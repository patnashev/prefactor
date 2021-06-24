
#include <cmath>
#include <iostream>

#include "gwnum.h"
#include "cpuid.h"
#include "stage1.h"
#include "exception.h"

using namespace arithmetic;

void report_factor(const Giant& f, InputNum& input);

void Stage1::init(InputNum& input, GWState& gwstate, int iterations, File* file, TaskState* state)
{
    Task::init(gwstate, iterations, file, state);
    _error_check = gwnear_fft_limit(gwstate.gwdata(), 1) == TRUE;
    _input = &input;
    _timer = getHighResTimer();
    _transforms = -(int)gwstate.handle.fft_count;
    _success = false;
}

void Stage1::done(const arithmetic::Giant& factor)
{
    _timer = (getHighResTimer() - _timer)/getHighResTimerFrequency();
    _transforms += (int)_gwstate->handle.fft_count;
    if (factor == 0 || factor == *_gwstate->N)
    {
        printf("All divisors of N < B1.\n");
        _success = true;
    }
    else if (factor != 1)
    {
        report_factor(factor, *_input);
        _success = true;
    }
    else
    {
        printf("No factors found, transforms: %d, time: %d s.\n", _transforms, (int)_timer);
    }
}

void Stage1::reinit_gwstate()
{
    _transforms += (int)_gwstate->handle.fft_count;
    _gwstate->done();
    _input->setup(*_gwstate);
    std::cout << "Using " << _gwstate->fft_description << std::endl;
    _error_check = gwnear_fft_limit(_gwstate->gwdata(), 1) == TRUE;
}

Giant Stage1::get_stage1_exp()
{
    int j, k;

    PrimeIterator it = primes().begin();
    for (; *it < 3; it++);

    int sqrtB1 = (int)sqrt(_B1);
    Giant tmp(GiantsArithmetic::default_arithmetic(), (int)(_B1/0.69/32) + 10);
    Giant tmp2(GiantsArithmetic::default_arithmetic(), tmp.capacity() < 8192 ? tmp.capacity() : 8192);
    tmp2 = 1;
    while (*it <= _B1)
    {
        // Building exponent with prime powers <= B1
        j = *it;
        if (*it <= sqrtB1)
        {
            k = _B1/(*it);
            while (j <= k)
                j *= *it;
        }
        tmp2 *= j;
        it++;
        if (*it > _B1 || tmp2.size() > 8190)
        {
            if (tmp == 0)
                tmp = tmp2;
            else
                tmp *= tmp2;
            tmp2 = 1;
        }
    }

    return tmp;
}


PM1Stage1::PM1Stage1(PrimeList& primes, int B1) : Stage1(primes, B1)
{
    _exp = get_stage1_exp();
}

void PM1Stage1::init(InputNum& input, GWState& gwstate, File* file)
{
    Stage1::init(input, gwstate, _exp.bitlen() - 1, file, read_state<State>(file));
    _state_update_period = MULS_PER_STATE_UPDATE;
    printf("%s, P-1 stage 1, B1 = %d%s", input.display_text().data(), _B1, input.gfn() != 0 ? ", GFN" : "");
    if (state() != nullptr)
        printf(", restarting at %.1f%%", 100.0*state()->iteration()/iterations());
    printf(".\n");
}

void PM1Stage1::execute()
{
    int i, len;

    GWNum X(gw());
    if (state() == nullptr)
    {
        i = 0;
        X = 3;
        gwset_carefully_count(gw().gwdata(), 30);
    }
    else
    {
        i = state()->iteration();
        X = state()->X();
    }
    gw().setmulbyconst(3);
    len = iterations();
    for (; i < len; i++, commit_execute<State>(i, X))
        gw().square(X, X, (_exp.bit(len - i - 1) ? GWMUL_MULBYCONST : 0) | GWMUL_STARTNEXTFFT_IF(!is_last(i)));
    if (i == len)
    {
        for (i = 2; i <= _B1; i <<= 1)
            gw().carefully().square(X, X, 0);
        for (i = 0; i < _input->gfn(); i++)
            gw().carefully().square(X, X, 0);
        i = len + 1;
        check();
        set_state(new State(i, X));
    }

    try
    {
        _V = (state()->X() + inv(state()->X(), gw().N()))%gw().N();
    }
    catch (const NoInverseException& e)
    {
        done(e.divisor);
        return;
    }
    done(gcd(state()->X() - 1, gw().N()));
}

PP1Stage1::PP1Stage1(PrimeList& primes, int B1, std::string& sP) : Stage1(primes, B1), _sP(sP)
{
    _Pa = stoi(sP);
    _Pb = 1;
    int i;
    for (i = sP[0] == '-' ? 1 : 0; isdigit(sP[i]); i++);
    if (sP[i] == '/')
        _Pb = stoi(sP.substr(i + 1));
}

void PP1Stage1::init(InputNum& input, GWState& gwstate, File* file)
{
    Stage1::init(input, gwstate, (int)std::expint(log(_B1)), file, read_state<State>(file));
    _state_update_period = MULS_PER_STATE_UPDATE/10;
    printf("%s, P+1 stage 1, B1 = %d, P = %s%s", input.display_text().data(), _B1, _sP.data(), input.gfn() != 0 ? ", GFN" : "");
    if (state() != nullptr)
        printf(", restarting at %.1f%%", 100.0*state()->iteration()/iterations());
    else
    {
        _P = _Pb;
        _P = _Pa*inv(std::move(_P), *gwstate.N)%(*gwstate.N);
    }
    printf(".\n");
}

void PP1Stage1::execute()
{
    int j, k;

    PrimeIterator it = primes().begin();

    LucasVArithmetic lucas(gw().carefully());
    LucasV V(lucas);
    if (state() == nullptr)
    {
        V.V() = _P;
        for (j = 2; j <= _B1; j <<= 1)
            lucas.dbl(V, V);
        for (j = 0; j < _input->gfn(); j++)
            lucas.dbl(V, V);
        it++;
        commit_execute<State>(1, V);
    }
    else
    {
        V.V() = state()->V();
        it += state()->iteration();
    }

    lucas.set_gw(gw());
    int sqrtB1 = (int)sqrt(_B1);
    while (*it <= _B1)
    {
        j = *it;
        lucas.mul(V, *it, it.pos(), V);
        if (*it <= sqrtB1)
        {
            k = _B1/(*it);
            while (j <= k)
            {
                lucas.mul(V, *it, it.pos(), V);
                j *= *it;
            }
        }
        it++;
        if (*it > _B1)
            _iterations = (int)it.pos();
        else if (_iterations == (int)it.pos())
            _iterations += 100;
        commit_execute<State>((int)it.pos(), V);
    }

    done(gcd(state()->V() - 2, gw().N()));
}

EdECMStage1::EdECMStage1(PrimeList& primes, int B1, int W) : Stage1(primes, B1), _W(W)
{
    if (_W < 2)
        _W = 2;
    if (_W > 16)
        _W = 16;
    Giant tmp = get_stage1_exp();
    get_NAF_W(_W, tmp, _NAF_W);
}

void EdECMStage1::init(InputNum& input, GWState& gwstate, File* file, arithmetic::Giant& X, arithmetic::Giant& Y, arithmetic::Giant& Z, arithmetic::Giant& T, arithmetic::Giant& EdD)
{
    Stage1::init(input, gwstate, (int)_NAF_W.size() - 1, file, read_state<State>(file));
    _state_update_period = MULS_PER_STATE_UPDATE/7/_W;
    printf("%s, EdECM stage 1, B1 = %d, W = %d", input.display_text().data(), _B1, _W);
    if (state() != nullptr)
        printf(", restarting at %.1f%%", 100.0*state()->iteration()/iterations());
    printf(".\n");

    _error_check = false;
    _X = X;
    _Y = Y;
    _Z = Z;
    _T = T;
    _EdD = EdD;
    ed.reset(new EdwardsArithmetic());
}

void EdECMStage1::setup()
{
    ed->set_gw(gw());
    if (!_ed_d)
    {
        _ed_d.reset(new GWNum(gw()));
        *_ed_d = _EdD;
    }

    if (_u.empty())
    {
        int i;
        std::vector<std::unique_ptr<arithmetic::EdPoint>> u;
        for (i = 0; i < (1 << (_W - 2)); i++)
            u.emplace_back(new EdPoint(*ed));
        u[0]->deserialize(_X, _Y, _Z, _T);
        ed->dbl(*u[0], *u[0], ed->ED_PROJECTIVE | GWMUL_STARTNEXTFFT);
        for (int i = 2; i <= _B1; i <<= 1)
            ed->dbl(*u[0], *u[0], ed->ED_PROJECTIVE | GWMUL_STARTNEXTFFT);
        ed->dbl(*u[0], *u[0], GWMUL_STARTNEXTFFT);

        if (_W > 2)
        {
            EdPoint u2(*ed);
            ed->dbl(*u[0], u2, GWMUL_STARTNEXTFFT);
            for (i = 1; i < (1 << (_W - 2)); i++)
                ed->add(*u[i - 1], u2, *u[i], GWMUL_STARTNEXTFFT);
        }
        if (!ed->on_curve(*u.back(), *_ed_d))
            throw TaskRestartException();

        try
        {
            ed->normalize(u.begin(), u.end(), 0);
        }
        catch (const NoInverseException& e)
        {
            done(e.divisor);
            return;
        }
        commit_setup();
        _u = std::move(u);
    }
}

void EdECMStage1::release()
{
    _ed_d.reset();
    _u.clear();
}

void EdECMStage1::execute()
{
    if (success())
        return;
    int i, j, len;

    ed->set_gw(gw());
    EdPoint p(*ed);
    if (state() == nullptr)
    {
        i = 0;
        p = *_u[_NAF_W.back()/2];
    }
    else
    {
        i = state()->iteration();
        p.deserialize(state()->X(), state()->Y(), state()->Z(), state()->T());
    }

    len = iterations();
    while (i < len)
    {
        int16_t cur = _NAF_W[len - i - 1];
        if (cur != 0)
        {
            for (j = 1; j < _W; j++)
                ed->dbl(p, p, GWMUL_STARTNEXTFFT | ed->ED_PROJECTIVE);
            ed->dbl(p, p, GWMUL_STARTNEXTFFT | (i < len - 1 ? 0 : ed->EDDBL_FOR_EXT_NORM_ADD));
            ed->add(p, *_u[abs(cur)/2], p, GWMUL_STARTNEXTFFT_IF(!is_last(i)) | (i < len - 1 ? ed->ED_PROJECTIVE : 0) | (cur < 0 ? ed->EDADD_NEGATIVE : 0));
        }
        else
            ed->dbl(p, p, GWMUL_STARTNEXTFFT_IF(!is_last(i)) | (i < len - 1 ? ed->ED_PROJECTIVE : 0));
        i++;
        //check();
        if (i%state_update_period() == 0 || i == len)
        {
            if (!ed->on_curve(p, *_ed_d))
                throw TaskRestartException();
            set_state(new State(i, p));
        }
    }
    
    done(gcd(state()->X(), gw().N()));
}
