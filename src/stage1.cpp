
#include <cmath>
#include <iostream>

#include "gwnum.h"
#include "cpuid.h"
#include "stage1.h"
#include "exception.h"

using namespace arithmetic;

void Stage1::init(InputNum* input, GWState* gwstate, File* file, TaskState* state, Logging* logging, int iterations)
{
    InputTask::init(input, gwstate, file, state, logging, iterations);
    _success = false;
}

void Stage1::done(const arithmetic::Giant& factor)
{
    InputTask::done();
    _logging->info("transforms: %d, time: %.1f s.\n", _transforms, _timer);
    if (factor == 0 || factor == *_gwstate->N)
    {
        _success = true;
        _logging->result(_success, "all divisors less than B1.\n");
        _logging->result_save("all divisors of " + _input->input_text()  + " are less than B1.\n");
    }
    else if (factor != 1)
    {
        _success = true;
        _logging->report_factor(*_input, factor);
    }
    else
    {
        //_logging->result(_success, "no factors found.\n");
    }
    _logging->set_prefix("");
}

void Stage1::reinit_gwstate()
{
    InputTask::reinit_gwstate();
    _error_check = gwnear_fft_limit(_gwstate->gwdata(), 1) == TRUE;
}

Giant Stage1::get_stage1_exp(int B1)
{
    return get_stage1_exp(2, B1, B1);
}

Giant Stage1::get_stage1_exp(uint64_t B1cur, uint64_t B1next, uint64_t B1max)
{
    std::vector<uint64_t> plist;
    PrimeIterator::get().sieve_range(B1cur + 1, B1next + 1, plist);
    auto it = plist.begin();

    uint64_t sqrtB1 = (uint64_t)(sqrt(B1max) + 0.5);
    Giant tmp(GiantsArithmetic::default_arithmetic(), (int)plist.size()*2 + 10);
    Giant tmp2(GiantsArithmetic::default_arithmetic(), 8192 < tmp.capacity() ? 8192 : tmp.capacity());
    tmp2 = 1;
    Giant g64;
    g64 = 1;
    g64 <<= 32;
    while (it != plist.end())
    {
        // Building exponent with prime powers <= B1
        uint64_t j = *it;
        if (*it <= sqrtB1)
        {
            uint64_t k = B1max/(*it);
            while (j <= k)
                j *= *it;
        }
        g64 = 1;
        if (j >= (1ULL << 32))
            g64 <<= 32;
        *(uint64_t*)g64.data() = j;
        tmp2 *= g64;
        it++;
        if (it == plist.end() || tmp2.size() > 8190)
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


PM1Stage1::PM1Stage1(int B1) : Stage1(B1)
{
    _squarings = (int)ceil(log2(B1));
    _exp = get_stage1_exp(B1);
    _exp_len = _exp.bitlen();
}

void PM1Stage1::init(InputNum* input, GWState* gwstate, File* file, Logging* logging)
{
    Stage1::init(input, gwstate, file, read_state<State>(file), logging, _exp_len - 1);
    _state_update_period = MULS_PER_STATE_UPDATE;
    _logging->set_prefix(input->display_text() + ", P-1 stage 1, ");
    _logging->info("B1 = %" PRId64 "%s", _B1, input->gfn() != 0 ? ", GFN" : "");
    if (state() != nullptr)
        _logging->info(", restarting at %.1f%%", 100.0*state()->iteration()/iterations());
    _logging->info(".\n");
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
        for (i = 0; i < _squarings + _input->gfn(); i++)
            gw().carefully().square(X, X, 0);
        i = len + 1;
        check();
        set_state<State>(i, X);
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
    if (_check_success)
        done(gcd(state()->X() - 1, gw().N()));
    else
        done(gw().popg() = 1);
}

PP1Stage1::PP1Stage1(int B1, std::string& sP) : Stage1(B1), _primes(B1), _sP(sP)
{
    _squarings = (int)ceil(log2(B1));
    _Pa = stoi(sP);
    _Pb = 1;
    int i;
    for (i = sP[0] == '-' ? 1 : 0; isdigit(sP[i]); i++);
    if (sP[i] == '/')
        _Pb = stoi(sP.substr(i + 1));
}

void PP1Stage1::init(InputNum* input, GWState* gwstate, File* file, Logging* logging)
{
    Stage1::init(input, gwstate, file, read_state<State>(file), logging, (int)std::expint(log(_B1)));
    _state_update_period = MULS_PER_STATE_UPDATE/10;
    _logging->set_prefix(input->display_text() + ", P+1 stage 1, ");
    _logging->info("B1 = %" PRId64 ", P = %s%s", _B1, _sP.data(), input->gfn() != 0 ? ", GFN" : "");
    if (state() != nullptr)
        _logging->info(", restarting at %.1f%%", 100.0*state()->iteration()/iterations());
    else
    {
        _P = _Pb;
        _P = _Pa*inv(std::move(_P), *gwstate->N)%(*gwstate->N);
    }
    _logging->info(".\n");
}

void PP1Stage1::execute()
{
    PrimeIterator it = _primes.begin();

    LucasVArithmetic lucas(gw().carefully());
    LucasV V(lucas);
    if (state() == nullptr)
    {
        V.V() = _P;
        for (int i = 0; i < _squarings + _input->gfn(); i++)
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
    uint64_t sqrtB1 = (uint64_t)sqrt(_B1);
    while (*it <= _B1)
    {
        lucas.mul(V, *it, (int)it.pos(), V);
        if (*it <= sqrtB1)
        {
            uint64_t j = *it;
            uint64_t k = _B1/(*it);
            while (j <= k)
            {
                lucas.mul(V, *it, (int)it.pos(), V);
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

    if (_check_success)
        done(gcd(state()->V() - 2, gw().N()));
    else
        done(gw().popg() = 1);
}

EdECMStage1::EdECMStage1(int B1, int W) : Stage1(B1), _W(W)
{
    _squarings = (int)ceil(log2(B1));
    if (_W < 2)
        _W = 2;
    if (_W > 16)
        _W = 16;
    Giant tmp = get_stage1_exp(B1);
    _exp_len = tmp.bitlen();
    get_NAF_W(_W, tmp, _NAF_W);
}

EdECMStage1::EdECMStage1(uint64_t B1cur, uint64_t B1next, uint64_t B1max, int max_size) : Stage1(B1next)
{
    if (B1cur < 2)
    {
        _squarings = (int)ceil(log2(B1max));
        B1cur = 2;
    }
    Giant tmp = get_stage1_exp(B1cur, B1next, B1max);
    _exp_len = tmp.bitlen();
    for (_W = 2; _W < 16 && 3 + 3*(1 << (_W - 1)) <= max_size && (15 << (_W - 2)) + _exp_len*(7 + 7/(_W + 1.0)) >(15 << (_W - 1)) + _exp_len*(7 + 7/(_W + 2.0)); _W++);
    get_NAF_W(_W, tmp, _NAF_W);
}

void EdECMStage1::init(InputNum* input, GWState* gwstate, File* file, Logging* logging, arithmetic::Giant* X, arithmetic::Giant* Y, arithmetic::Giant* Z, arithmetic::Giant* T, arithmetic::Giant* EdD)
{
    Stage1::init(input, gwstate, file, read_state<State>(file), logging, (int)_NAF_W.size() - 1);
    _state_update_period = MULS_PER_STATE_UPDATE/_W;
    _logging->set_prefix(input->display_text() + ", EdECM stage 1, ");
    _logging->info("B1 = %" PRId64 ", %d bits, W = %d", _B1, _exp_len, _W);
    if (state() != nullptr)
        _logging->info(", restarting at %.1f%%", 100.0*state()->iteration()/iterations());
    _logging->info(".\n");

    _error_check = false;
    _X = X;
    _Y = Y;
    _Z = Z;
    _T = T;
    _EdD = EdD;
}

void EdECMStage1::setup()
{
    if (!ed)
        ed.reset(new EdwardsArithmetic());
    ed->set_gw(gw());
    if (!_ed_d)
    {
        _ed_d.reset(new GWNum(gw()));
        *_ed_d = *_EdD;
    }

    if (_u.empty())
    {
        int i;
        std::vector<std::unique_ptr<arithmetic::EdPoint>> u;
        for (i = 0; i < (1 << (_W - 2)); i++)
            u.emplace_back(new EdPoint(*ed));
        u[0]->deserialize(*_X, *_Y, *_Z, *_T);
        for (i = 0; i < _squarings; i++)
            ed->dbl(*u[0], *u[0], (i < _squarings - 1 ? ed->ED_PROJECTIVE : 0) | GWMUL_STARTNEXTFFT);

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
    ed.release();
}

void EdECMStage1::execute()
{
    if (success())
        return;
    int i, j, len;

    ed->set_gw(gw());
    EdPoint p(*ed);
    if (iterations() < 0)
    {
        i = -1;
        set_state<State>(i, *_u[0]);
    }
    else if (state() == nullptr)
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
        if (i%state_update_period() == 0 || i == len || abort_flag())
        {
            if (!ed->on_curve(p, *_ed_d))
                throw TaskRestartException();
            set_state<State>(i, p);
        }
    }
    
    if (_check_success)
        done(gcd(state()->X(), gw().N()));
    else
        done(gw().popg() = 1);
}
