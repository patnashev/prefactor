
#include <deque>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <functional>
#include <string.h>

#include "gwnum.h"
#include "cpuid.h"
#include "stage2poly.h"
#include "exception.h"

using namespace arithmetic;

// http://cr.yp.to/arith/scaledmod-20040820.pdf

#ifdef _DEBUG
//#define DEBUG_STAGE2POLY_FAST
#define DEBUG_STAGE2POLY_INV
#define DEBUG_STAGE2POLY_MOD
#define DEBUG_STAGE2POLY_REM
#define DEBUG_STAGE2POLY_ROOT
#define DEBUG_STAGE2POLY_QUEUE
#define DEBUG_STAGE2POLY_CONVERT 0
#else
//#define DEBUG_STAGE2POLY_FAST
//#define DEBUG_STAGE2POLY_INV
#define DEBUG_STAGE2POLY_CONVERT 0
#endif

template<class Element>
void Stage2Poly<Element>::init(InputNum* input, GWState* gwstate, File* file, TaskState* state, Logging* logging)
{
#ifdef _DEBUG
    //_tinypoly_power = 0;
    //_poly_check = false;
#endif
    if (!gwstate->polymult)
    {
        gwstate->polymult = true;
        gwstate->done();
        input->setup(*gwstate);
    }
    if (PolyMult::max_polymult_output(*gwstate) < 2*(1 << poly_power()))
    {
        gwstate->done();
        gwstate->next_fft_count++;
        input->setup(*gwstate);
        std::string prefix = logging->prefix();
        logging->set_prefix("");
        logging->warning("Switching to %s\n", gwstate->fft_description.data());
        logging->set_prefix(prefix);
        logging->report_param("fft_desc", gwstate->fft_description);
        logging->report_param("fft_len", gwstate->fft_length);
    }

    _poly_gwstate.clear();
    GWState* cur = gwstate;
    for (int i = 1; i <= poly_power(); i++)
    {
        if (2*(1 << i) > PolyMult::max_polymult_output(*cur) || i == DEBUG_STAGE2POLY_CONVERT)
        {
            _poly_gwstate.emplace_back();
            _poly_gwstate.back().copy(*cur);
            _poly_gwstate.back().next_fft_count++;
            input->setup(_poly_gwstate.back());
            std::string prefix = logging->prefix();
            logging->set_prefix("");
            logging->warning("Using also %s\n", _poly_gwstate.back().fft_description.data());
            logging->set_prefix(prefix);
            cur = &_poly_gwstate.back();
        }
    }

    _first_D = (int)((_B1 + _D/2)/_D);
    _last_D = (int)((_B2 + _D/2)/_D);
    Stage2::init(input, gwstate, file, state, logging, _last_D - _first_D + 1);

    if (poly_check())
        _poly_check_index = Giant::rnd(32).data()[0]%_poly_mod_degree;

    //_state_update_period = poly_degree();
    logging->info("B2 = %" PRId64, _B2);
    if (this->state() != nullptr)
        logging->info(", restarting at %.1f%%", 100.0*this->state()->iteration()/iterations());
    logging->info(".\n");
    logging->info("polynomial mode, D = %d, degree %d/%d, %d steps in batches of %d.\n", _D, 1 << poly_power(), 1 << _tinypoly_power, _last_D - _first_D + 1, 1 << _smallpoly_power);
    if ((_last_D - _first_D + 1) < _poly_mod_degree)
        logging->info("recommended B2 = %" PRId64 ".\n", (_first_D + _poly_mod_degree - 1)*(uint64_t)_D + _D/2 - 1);
    else if ((_last_D - _first_D + 1 - _poly_mod_degree)%(1 << poly_power()) != 0)
        logging->info("recommended B2 = %" PRId64 ".\n", (_last_D + (1 << poly_power()) - (_last_D - _first_D + 1 - _poly_mod_degree)%(1 << poly_power()))*(uint64_t)_D + _D/2 - 1);
    if (_poly_threads > 1)
        logging->info("%d threads.\n", _poly_threads);
}

template<class Element>
void Stage2Poly<Element>::done(const arithmetic::Giant& factor)
{
    for (int i = 0; i < _workers.size(); i++)
        gwclone_merge_stats(_gwstate->gwdata(), _workers[i]->gw().gwdata());
    Stage2::done(factor);
}

template<class Element>
void Stage2Poly<Element>::smallpoly_init(int degree)
{
    int i, j, k;
    _smallpoly_count = ((degree + (1 << _smallpoly_power) - 1) >> _smallpoly_power);
    int last_tiny_count = (((degree - ((_smallpoly_count - 1) << _smallpoly_power)) + (1 << _tinypoly_power) - 1) >> _tinypoly_power);
    _smallpoly.resize(_smallpoly_count);
    _smallpoly_alloc.resize(_smallpoly_count);
    for (k = 0; k < _smallpoly_alloc.size(); k++)
    {
        _smallpoly_alloc[k].resize(_smallpoly_power + 1);
        for (j = 0; j <= _tinypoly_power; j++)
            _smallpoly_alloc[k][j].resize(k < _smallpoly_alloc.size() - 1 ? (size_t)1 << (_smallpoly_power - _tinypoly_power) : last_tiny_count, Poly(_poly_mult[j]));
        for (j = _tinypoly_power + 1; j <= _smallpoly_power; j++)
            _smallpoly_alloc[k][j].resize((_smallpoly_alloc[k][j - 1].size() + 1)/2, Poly(_poly_mult[j]));
        for (j = 0; j <= _smallpoly_power; j++)
            for (i = 0; i < _smallpoly_alloc[k][j].size(); i++)
                _poly_mult[j].free(_smallpoly_alloc[k][j][i]);
    }
}

template<class Element>
void Stage2Poly<Element>::smallpoly_init_level(gwarray data, int data_level)
{
    int i, j, k;
    j = data_level;
    for (k = 0; k < _smallpoly_alloc.size(); k++)
        for (i = 0; i < _smallpoly_alloc[k][j].size(); i++)
                _poly_mult[j].init(data + (k << _smallpoly_power) + (i << (j < _tinypoly_power ? _tinypoly_power : j)), (size_t)1 << (j < _tinypoly_power ? _tinypoly_power : j), false, true, _smallpoly_alloc[k][j][i]);
}

void build_tree(std::vector<std::vector<Poly>>& tree, std::vector<arithmetic::PolyMult>& poly_mult, int base_level, bool base_tree, std::function<bool(int)> preserve_predicate, std::vector<gwarray>& gwarrays, int level_size, std::vector<gwarray>& gwarrays_tmp)
{
    int i, j;
    for (j = base_level + 1; j < (base_tree ? base_level : 0) + tree.size(); j++)
    {
        PolyMult& pm = poly_mult[j - 1];
        PolyMult& pm_next = poly_mult[j];
        bool convert = &pm.gw() != &pm_next.gw();
        bool preserve = preserve_predicate(j);
        int options = !convert ? POLYMULT_NEXTFFT : 0;
        std::vector<Poly>& cur = tree[j - (base_tree ? base_level : 0) - 1];
        std::vector<Poly>& next = tree[j - (base_tree ? base_level : 0)];
        Poly tmp(pm);

        next.resize((cur.size() + 1)/2, Poly(pm_next));
        if ((convert || preserve) && next[0].empty())
        {
            if (gwarrays[j] == nullptr)
                gwarrays[j] = gwalloc_array(pm_next.gw().gwdata(), level_size);
            for (i = 0; i < next.size(); i++)
                pm_next.init(gwarrays[j] + (i << j), level_size >= ((i + 1) << j) ? (size_t)1 << j : level_size - (i << j), false, false, next[i]);
        }
        if (convert && preserve)
        {
            if (gwarrays_tmp[j - 1] == nullptr)
                gwarrays_tmp[j - 1] = gwalloc_array(pm.gw().gwdata(), 1 << j);
            pm.init(gwarrays_tmp[j - 1], (size_t)1 << j, false, false, tmp);
        }
        for (i = 0; i < next.size(); i++)
        {
            if (2*i + 1 == cur.size())
            {
                if (convert)
                    pm.convert(cur[2*i], pm_next, next[i]);
                else if (preserve)
                    pm.copy(cur[2*i], next[i]);
                else
                    next[i] = std::move(cur[2*i]);
                continue;
            }
            Poly& res = convert ? tmp : next[i];
            //if (preserve && workstage == 1)
            //    pm.preprocess_and_mul(cur[2*i], cur[2*i + 1], res, 1 << j, options | POLYMULT_PRE_FFT); else
            if (preserve)
                pm.mul(cur[2*i], cur[2*i + 1], res, options);
            else
                pm.mul(std::move(cur[2*i]), std::move(cur[2*i + 1]), res, options);
            if (convert)
                pm.convert(res, pm_next, next[i]);
            if (!preserve)
                pm.free(cur[2*i]);
            if (!preserve)
                pm.free(cur[2*i + 1]);
        }
    }
}

void build_tree_tinypoly(std::vector<std::vector<Poly>>& tree, int k, int tiny_power, std::vector<arithmetic::PolyMult>& poly_mult, int base_level, std::function<bool(int)> preserve_predicate, std::vector<gwarray>& gwarrays, int level_size, std::vector<gwarray>& gwarrays_tmp, GWNum* check_root, GWNum* check)
{
    int i, j;
    std::unique_ptr<GWNum> a1;
    for (j = base_level + 1; j <= tiny_power && j < tree.size(); j++)
    {
        PolyMult& pm = poly_mult[j - 1];
        PolyMult& pm_next = poly_mult[j];
        bool convert = &pm.gw() != &pm_next.gw();
        bool preserve = preserve_predicate(j);
        int options = !convert ? POLYMULT_NEXTFFT : 0;
        Poly& cur = tree[j - 1][k];
        Poly& next = tree[j][k];
        Poly tmp(pm);
        void* plan = NULL;

        if ((convert || preserve) && next.empty())
        {
            if (gwarrays[j] == nullptr)
                gwarrays[j] = gwalloc_array(pm_next.gw().gwdata(), level_size);
            pm_next.init(gwarrays[j] + (k << tiny_power), cur.size(), false, cur.monic(), next);
        }
        else if (!next.empty())
            pm_next.alloc(next, cur.size());
        else
            pm.init(cur.data(), cur.size(), false, cur.monic(), next);
        if (convert && preserve)
        {
            if (gwarrays_tmp[j - 1] == nullptr)
                gwarrays_tmp[j - 1] = gwalloc_array(pm.gw().gwdata(), 1 << tiny_power);
            pm.init(gwarrays_tmp[j - 1], cur.size(), false, cur.monic(), tmp);
        }
        else if (convert)
            pm.init(cur.data(), cur.size(), false, cur.monic(), tmp);

        Poly& res = convert ? tmp : next;
        int degree = (1 << (j - 1));
        for (int offset = 0; offset < cur.size(); offset += 2*degree)
        {
            if (offset + degree >= cur.size())
            {
                if (degree == 1 && check)
                    gwsubmul4(pm.gw().gwdata(), cur.data()[offset], **check_root, **check, **check, GWMUL_STARTNEXTFFT | GWMUL_FFT_S1);
                if (preserve)
                    for (i = offset; i < cur.size(); i++)
                        gwcopy(pm.gw().gwdata(), cur.data()[i], res.data()[i]);
                continue;
            }
            int sb = offset + 2*degree <= cur.size() ? degree : cur.size() - offset - degree;
            if (degree == 1)
            {
                if (preserve)
                {
                    gwadd3o(pm.gw().gwdata(), cur.data()[offset], cur.data()[offset + 1], res.data()[offset + 1], GWADD_FORCE_NORMALIZE);
                    if (!convert)
                        gwfft(pm.gw().gwdata(), res.data()[offset + 1], res.data()[offset + 1]);
                    gwmul3(pm.gw().gwdata(), cur.data()[offset], cur.data()[offset + 1], res.data()[offset], GWMUL_STARTNEXTFFT_IF(!convert) | GWMUL_FFT_S1 | GWMUL_FFT_S2);
                    if (!convert)
                        gwfft(pm.gw().gwdata(), res.data()[offset], res.data()[offset]);
                }
                else
                {
                    if (!a1)
                        a1.reset(new GWNum(pm.gw()));
                    gwadd3o(pm.gw().gwdata(), cur.data()[offset], cur.data()[offset + 1], **a1, GWADD_FORCE_NORMALIZE);
                    if (check)
                    {
                        gwsubmul4(pm.gw().gwdata(), cur.data()[offset], **check_root, **check, **check, GWMUL_STARTNEXTFFT | GWMUL_FFT_S1);
                        gwsubmul4(pm.gw().gwdata(), cur.data()[offset + 1], **check_root, **check, **check, GWMUL_STARTNEXTFFT | GWMUL_FFT_S1);
                    }
                    gwmul3(pm.gw().gwdata(), cur.data()[offset], cur.data()[offset + 1], res.data()[offset], GWMUL_STARTNEXTFFT_IF(!convert));
                    if (!convert)
                        gwfft(pm.gw().gwdata(), res.data()[offset], res.data()[offset]);
                    if (!convert)
                        gwfft(pm.gw().gwdata(), **a1, res.data()[offset + 1]);
                    else
                        gwcopy(pm.gw().gwdata(), **a1, res.data()[offset + 1]);
                }
            }
            else
            {
                polymult(pm.pmdata(), cur.data() + offset, degree, cur.data() + offset + degree, sb, res.data() + offset, degree + sb, options | POLYMULT_INVEC1_MONIC | POLYMULT_INVEC2_MONIC | (sb == degree ? POLYMULT_SAVE_PLAN : 0) | (sb == degree && offset > 0 ? POLYMULT_USE_PLAN : 0));
                if (sb == degree)
                    plan = pm.pmdata()->plan;
            }
        }
        if (plan)
            free(plan);
        if (convert)
            pm.convert(res, pm_next, next);
        if (!preserve)
            pm.free(cur);
    }
}

void rem_tree(std::vector<Poly>& rem, int k, std::vector<std::vector<Poly>>& tree, std::vector<arithmetic::PolyMult>& poly_mult, int cur_level, int base_level, int stop_level, std::vector<gwarray>& gwarrays_down, int level_size)
{
    int i, j;
    std::vector<Poly> rem_prev;
    rem_prev.reserve(rem.capacity());
    for (j = cur_level; gwarrays_down[j] == nullptr; j++);
    gwarray last_gwarray = gwarrays_down[j];
    for (j = cur_level - base_level; j >= stop_level; j--)
    {
        PolyMult& pm = poly_mult[base_level + j];
        PolyMult& pm_prev = poly_mult[base_level + j > 0 ? base_level + j - 1 : 0];
        bool convert = &pm.gw() != &pm_prev.gw();
        int options = !convert ? POLYMULT_NEXTFFT : 0;
        Poly tmp1(pm), tmp2(pm);
        gwarray tmp_gwarray;

        if (convert)
        {
            tmp_gwarray = last_gwarray;
            if (gwarrays_down[base_level + j - 1] == nullptr)
                gwarrays_down[base_level + j - 1] = gwalloc_array(pm_prev.gw().gwdata(), level_size);
            last_gwarray = gwarrays_down[base_level + j - 1];
        }
        for (i = 0; i < rem.size(); i++)
        {
            if (convert)
                pm.init(tmp_gwarray + (k << (cur_level + 1)), ((size_t)1 << (base_level + j)), false, true, tmp1);
            Poly& rem1 = rem_prev.emplace_back(pm_prev);
            pm_prev.init(last_gwarray + (k << (cur_level + 1)) + ((2*i) << (base_level + j)), ((size_t)1 << (base_level + j)), false, true, rem1);
            Poly& res1 = convert ? tmp1 : rem1;
            if (2*i + 1 == tree[j].size())
            {
                res1 = std::move(rem[i]);
                res1 >>= (1 << (base_level + j));
            }
            else
            {
                if (convert)
                    pm.init(tmp_gwarray + (k << (cur_level + 1)) + ((size_t)1 << (base_level + j)), ((size_t)1 << (base_level + j)), false, true, tmp2);
                Poly& rem2 = rem_prev.emplace_back(pm_prev);
                pm_prev.init(last_gwarray + (k << (cur_level + 1)) + ((2*i + 1) << (base_level + j)), ((size_t)1 << (base_level + j)), false, true, rem2);
                Poly& res2 = convert ? tmp2 : rem2;
                pm.mul_twohalf(rem[i], tree[j][2*i + 1], tree[j][2*i], res1, res2, 1 << (base_level + j), options);
                if (convert)
                    pm.convert(res2, pm_prev, rem2);
            }
            if (convert)
                pm.convert(res1, pm_prev, rem1);
        }
        std::swap(rem, rem_prev);
        rem_prev.clear();
    }
}

void rem_tree_tinypoly(std::vector<Poly>& poly_rem, int k, std::vector<std::vector<Poly>>& tree, int tiny_power, std::vector<arithmetic::PolyMult>& poly_mult, std::vector<gwarray>& gwarrays_down, int level_size)
{
    int i, j;
    GWNum a0(poly_mult[0].gw());
    Poly& rem = poly_rem[k];
    int size = rem.size();
    for (j = tiny_power - 1; gwarrays_down[j] == nullptr; j++);
    gwarray last_gwarray = gwarrays_down[j] + (k << tiny_power);
    for (; last_gwarray[0] != rem.data()[0]; last_gwarray++);
    poly_mult[tiny_power - 1].init(last_gwarray, (size_t)1 << tiny_power, false, rem.monic(), rem);
    for (j = tiny_power - 1; j >= 0; j--)
    {
        PolyMult& pm = poly_mult[j];
        PolyMult& pm_prev = poly_mult[j > 0 ? j - 1 : 0];
        bool convert = &pm.gw() != &pm_prev.gw();
        int options = !convert ? POLYMULT_NEXTFFT : 0;
        Poly rem_prev(pm_prev);
        void* plan = NULL;

        if (convert)
        {
            if (gwarrays_down[j - 1] == nullptr)
                gwarrays_down[j - 1] = gwalloc_array(pm_prev.gw().gwdata(), level_size);
            pm_prev.init(gwarrays_down[j - 1] + (k << tiny_power), (size_t)1 << tiny_power, false, true, rem_prev);
        }

        int degree = (1 << j);
        for (int offset = 0; offset < (1 << tiny_power); offset += 2*degree)
        {
            if (size == 1)
            {
                gwadd3o(pm.gw().gwdata(), tree[j][k].data()[offset + degree - 1], rem.data()[offset], rem.data()[offset + degree], GWADD_NORMALIZE_IF(convert));
                gwadd3o(pm.gw().gwdata(), tree[j][k].data()[offset + 2*degree - 1], rem.data()[offset], rem.data()[offset], GWADD_NORMALIZE_IF(convert));
            }
            else if (degree == 1)
            {
                gwmuladd4(pm.gw().gwdata(), tree[j][k].data()[offset + 1], rem.data()[offset + 1], rem.data()[offset], *a0, GWMUL_STARTNEXTFFT_IF(!convert) | GWMUL_FFT_S2 | GWMUL_FFT_S3);
                gwmuladd4(pm.gw().gwdata(), tree[j][k].data()[offset], rem.data()[offset + 1], rem.data()[offset], rem.data()[offset + 1], GWMUL_STARTNEXTFFT_IF(!convert));
                if (!convert)
                    gwfft(pm.gw().gwdata(), rem.data()[offset + 1], rem.data()[offset + 1]);
                if (!convert)
                    gwfft(pm.gw().gwdata(), *a0, rem.data()[offset]);
                else
                    gwcopy(pm.gw().gwdata(), *a0, rem.data()[offset]);
            }
            else
            {
                polymult_arg args[2];
                polymult_arg& arg_b = args[0];
                polymult_arg& arg_c = args[1];
                arg_b.invec2 = tree[j][k].data() + offset + degree;
                arg_b.invec2_size = degree;
                arg_b.outvec = rem.data() + offset;
                arg_b.outvec_size = size < degree ? size : degree;
                arg_b.fmavec = nullptr;
                arg_b.circular_size = (size > degree ? 2*degree : 0);
                arg_b.first_mulmid = 0;
                arg_b.options = (size > degree ? POLYMULT_CIRCULAR : 0) | POLYMULT_INVEC2_MONIC;
                arg_c.invec2 = tree[j][k].data() + offset;
                arg_c.invec2_size = degree;
                arg_c.outvec = rem.data() + offset + degree;
                arg_c.outvec_size = size < degree ? size : degree;
                arg_c.fmavec = nullptr;
                arg_c.circular_size = (size > degree ? 2*degree : 0);
                arg_c.first_mulmid = 0;
                arg_c.options = (size > degree ? POLYMULT_CIRCULAR : 0) | POLYMULT_INVEC2_MONIC;
                polymult_several(pm.pmdata(), rem.data() + offset, size < 2*degree ? size : 2*degree, args, 2, options | POLYMULT_MULHI | (size < 2*degree ? POLYMULT_INVEC1_MONIC : 0) | POLYMULT_SAVE_PLAN | (offset > 0 ? POLYMULT_USE_PLAN : 0));
                plan = pm.pmdata()->plan;
            }
            if (convert)
            {
                for (i = 0; i < size && i < degree; i++)
                    gwconvert(pm.gw().gwdata(), pm_prev.gw().gwdata(), rem.data()[offset + i], rem_prev.data()[offset + i]);
                for (i = 0; i < size && i < degree; i++)
                    gwconvert(pm.gw().gwdata(), pm_prev.gw().gwdata(), rem.data()[offset + degree + i], rem_prev.data()[offset + degree + i]);
            }
        }
        if (plan)
            free(plan);
        if (convert)
            rem = std::move(rem_prev);
    }
}

template<class Element>
void Stage2Poly<Element>::setup()
{
    int i, j, k;
    double timer = getHighResTimer();
    int transforms = -(int)_gwstate->gwdata()->fft_count;
    Element Xn(arithmetic());
    Element Xn1(arithmetic());
#ifdef DEBUG_STAGE2POLY_MOD
    Element X1(*_X1);
#endif

    std::vector<std::unique_ptr<Element>> precomp(_D/4);
    precomp[0].reset(new Element(arithmetic()));
    *precomp[0] = *_X1;
    Xn = *_X1;
    arithmetic().dbl(Xn, Xn); // V_2
    arithmetic().dbl(Xn, Xn1); // V_4

    int dist = 1;
    int precomp_size = 1;
    for (i = 1; i < _D/4; i++)
    {
        if (i == 2*dist && precomp_size > 20*_poly_threads)
            break;
        if (gcd(2*i + 1, 2*dist) != 1)
            continue;
        // V_{2i+1}
        precomp[i].reset(new Element(arithmetic()));
        precomp_size++;
        arithmetic().add(Xn, *precomp[i - dist], *precomp[i >= 2*dist ? i - 2*dist : 2*dist - i - 1], *precomp[i]);
        if ((i + 1)%dist == 0 && _D%(i + 1) == 0 && gcd((i + 1)/dist, 2*dist) == 1)
        {
            int pd = (i + 1)/dist;
            arithmetic().add(*precomp[dist*(pd - 1)/2 - 1], *precomp[dist*(pd + 1)/2], Xn1, Xn1);
            arithmetic().add(*precomp[dist*(pd - 1)/2], *precomp[dist*(pd + 1)/2], Xn, Xn);
            swap(Xn, Xn1);
            dist = i + 1;
            for (j = 0; j < i; j++)
                if (precomp[j] && gcd(2*j + 1, 2*dist) != 1)
                {
                    precomp[j].reset();
                    precomp_size--;
                }
        }
    }

    std::unique_lock<std::mutex> lock(_work_mutex);
    if (i < _D/4)
    {
        *_X1 = Xn;
        for (i = 0; i < dist; i++)
            if (precomp[i])
            {
                _workqueue.emplace_back();
                _workqueue.back().Xn = std::move(precomp[i]);
                _workqueue.back().Xn1 = std::move(precomp[i + dist]);
                _workqueue.back().count = _D/4/dist + 1;
                _workqueue.back().n = 2*i + 1;
                _workqueue.back().distance = 2*dist;
                while (_workqueue.back().count > 0 && ((_workqueue.back().count - 1)*_workqueue.back().distance + _workqueue.back().n >= _D/2 || gcd((_workqueue.back().count - 1)*_workqueue.back().distance + _workqueue.back().n, _D) != 1))
                    _workqueue.back().count--;
                if (_workqueue.back().count == 0)
                    _workqueue.pop_back();
            }
    }
    else
    {
        _X1 = nullptr;
        for (i = 0; i < _D/4; i++)
            if (precomp[i] && gcd(2*i + 1, _D) == 1)
            {
                _workqueue.emplace_back();
                _workqueue.back().Xn = std::move(precomp[i]);
                _workqueue.back().count = 1;
                _workqueue.back().n = 2*i + 1;
            }
    }
    
    _poly_mult.reserve(poly_power() + 1);
    _poly_mult.emplace_back(gw(), _poly_threads);
    int cur_gwstate = -1;
    for (int i = 1; i <= poly_power(); i++)
    {
        if (2*(1 << i) > _poly_mult[i - 1].max_output() || i == DEBUG_STAGE2POLY_CONVERT)
        {
            cur_gwstate++;
            _poly_gw.emplace_back(_poly_gwstate[cur_gwstate]);
            _poly_mult.emplace_back(_poly_gw.back(), _poly_threads);
        }
        else
            _poly_mult.emplace_back(_poly_mult[i - 1].gw(), _poly_threads);
    }

    GWASSERT((_poly_mod_degree & ((1 << _tinypoly_power) - 1)) == 0);
    _gwarrays.resize(poly_power() + 1);
    _gwarrays_tmp.resize(poly_power() + 1);
    _poly_mod.resize(poly_power() + 1);
    _element_map.resize(_poly_mod_degree);

    smallpoly_init(_poly_mod_degree);
    _gwarrays[0] = gwalloc_array(_poly_mult[0].gw().gwdata(), _poly_mod_degree);
    smallpoly_init_level(_gwarrays[0], 0);
    for (j = 1; j <= _smallpoly_power; j++)
    {
        PolyMult& pm = _poly_mult[j - 1];
        PolyMult& pm_next = _poly_mult[j];
        bool convert = &pm.gw() != &pm_next.gw();
        bool preserve = (j - 1)%_poly_optmem_small == 0;

        if (preserve || convert)
        {
            _gwarrays[j] = gwalloc_array(pm_next.gw().gwdata(), _poly_mod_degree);
            smallpoly_init_level(_gwarrays[j], j);
        }
    }
    
    _workstage = 1;
    lock.unlock();
    _workstage_signal.notify_all();
    _workers[0]->run();
    lock.lock();
    _workdone_signal.wait(lock, [&] {return _thread_exception || _smallpoly_count == 0; });
    lock.unlock();
    if (_thread_exception)
        std::rethrow_exception(_thread_exception);
    GWASSERT(_workqueue.empty());
    _logging->debug("all threads done.\n");

    _smallpoly_mod = std::move(_smallpoly);
    for (k = 0; k < _smallpoly_mod.size(); k++)
        _poly_mod[_smallpoly_power].emplace_back(std::move(_smallpoly_mod[k][_smallpoly_power][0]));

    build_tree(_poly_mod, _poly_mult, _smallpoly_power, false, [&](int j) { return (j - _smallpoly_power - 1)%_poly_optmem == 0 && _mem_model > -2; }, _gwarrays, _poly_mod_degree, _gwarrays_tmp);
    GWASSERT(_poly_mod[poly_power()][0].degree() == _poly_mod_degree);
#ifdef DEBUG_STAGE2POLY_FAST
    if (poly_check())
    {
        GWNum tmp(gw());
        gw().setmulbyconst(-1);
        gwunfft2(gw().gwdata(), *_smallpoly_mod[_poly_check_index >> _smallpoly_power][0][(_poly_check_index & ((1 << _smallpoly_power) - 1)) >> _tinypoly_power].at(_poly_check_index & ((1 << _tinypoly_power) - 1)), *tmp, GWMUL_MULBYCONST);
        if (_poly_mod[poly_power()][0].eval(tmp) != 0)
            _logging->error("incorrect modulus\n");
    }
#endif

#ifdef DEBUG_STAGE2POLY_MOD
    std::vector<std::unique_ptr<Element>> elements;
    for (i = 1, j = 0; i < _D/2; i++)
        if (gcd(_D, i) == 1)
            arithmetic().mul(X1, i, *elements.emplace_back(new Element(arithmetic())));
    GWASSERT(elements.size() == _poly_mod_degree);
    Poly poly(_poly_mult[0], _poly_mod_degree, false);
    _workers[0]->elements_to_gwnums(elements, elements.size(), poly.data());
    elements.clear();

    for (i = 0; i < _poly_mod_degree; i++)
    {
        gw().neg((GWNum&)poly.at(i), (GWNum&)poly.at(i));
        if (&_poly_mult[poly_power()].gw() != _gw)
        {
            GWNum tmp(_poly_mult[poly_power()].gw());
            gwconvert(gw().gwdata(), _poly_mult[poly_power()].gw().gwdata(), *poly.at(i), *tmp);
            GWASSERT(_poly_mod[poly_power()][0].eval(tmp) == 0);
        }
        else
            GWASSERT((_poly_mult[poly_power()].gw().popg() = _poly_mod[poly_power()][0].eval((GWNum&)poly.at(i)))%_poly_mult[poly_power()].gw().N() == 0);
    }
#endif
#ifdef DEBUG_STAGE2POLY_REM
    Element Xd(arithmetic());
    arithmetic().mul(X1, _D, Xd);
    arithmetic().mul(Xd, _first_D, Xn, Xn1);
    for (i = _first_D; i <= _last_D; i++)
    {
        elements.emplace_back(new Element(arithmetic()));
        if (i > 0)
            arithmetic().add(Xd, Xn1, Xn, *elements.back());
        else
            arithmetic().dbl(Xn1, *elements.back());
        swap(Xn, Xn1);
        swap(Xn1, *elements.back());
    }
    Poly poly_roots(_poly_mult[0], elements.size(), false);
    _workers[0]->elements_to_gwnums(elements, elements.size(), poly_roots.data());
    elements.clear();

    _gwarrays_tmp[0] = gwalloc_array(gw().gwdata(), 1 << poly_power());
    for (k = 0; k < _smallpoly_mod.size(); k++)
    {
        Poly poly_rem(_poly_mult[0]);
        _poly_mult[0].init(_gwarrays_tmp[0] + (k << _smallpoly_power), _smallpoly_mod[k][0].size() << _tinypoly_power, false, false, poly_rem);
        for (i = 0; i < poly_rem.size(); i++)
        {
            (GWNum&)poly_rem.at(i) = 1;
            for (j = 0; j < poly_roots.size(); j++)
                gw().submul((GWNum&)poly_roots.at(j), (GWNum&)_smallpoly_mod[k][0][i >> _tinypoly_power].at(i & ((1 << _tinypoly_power) - 1)), (GWNum&)poly_rem.at(i), (GWNum&)poly_rem.at(i), j < poly_roots.size() - 1 ? GWMUL_STARTNEXTFFT : 0);
        }
        _smallpoly_mod[k][_smallpoly_power][0] = std::move(poly_rem);
    }
#endif

    _logging->debug("calculating reciprocal...\n");
    PolyMult& pm = _poly_mult[poly_power()];
    Poly& modulus = _poly_mod[poly_power()][0];
    Poly reciprocal(pm);
    _gwarrays.push_back(gwalloc_array(pm.gw().gwdata(), 1 << poly_power()));
    pm.init(_gwarrays.back(), (1 << poly_power()) - 1, false, true, reciprocal);
    pm.reciprocal(modulus, reciprocal, POLYMULT_NEXTFFT);
    int degree = (1 << (poly_power() + 1)) - modulus.degree();
    if (reciprocal.degree() < degree)
        reciprocal <<= degree - reciprocal.degree();

    _logging->debug("preprocessing...\n");
    _modulus.reset(new Poly(pm));
    _reciprocal.reset(new Poly(pm));
    pm.preprocess(modulus, *_modulus, 1 << (poly_power() + 1), _mem_model >= 0 ? POLYMULT_PRE_FFT : 0);
    pm.preprocess(reciprocal, *_reciprocal, 1 << (poly_power() + 1), _mem_model >= 0 ? POLYMULT_PRE_FFT : 0);
#ifdef DEBUG_STAGE2POLY_MOD
    Poly polyT(pm);
    pm.mul_range(*_reciprocal, *_modulus, polyT, 1 << poly_power(), 1 << poly_power(), 0);
    for (i = (1 << poly_power()) - 1; i >= 1; i--)
        GWASSERT(polyT.at(i) == 0);
#endif
#ifdef DEBUG_STAGE2POLY_FAST
    if (poly_check())
    {
        Poly polyT(pm);
        pm.mul_range(*_reciprocal, *_modulus, polyT, 1 << poly_power(), 1 << poly_power(), 0);
        for (i = (1 << poly_power()) - 1; i >= 1; i--)
            if (polyT.at(i) != 0)
            {
                _logging->error("incorrect reciprocal\n");
                break;
            }
    }
#endif

    for (i = 0; i < _workers.size(); i++)
        gwclone_merge_stats(_gwstate->gwdata(), _workers[i]->gw().gwdata());
    transforms += (int)_gwstate->gwdata()->fft_count;
    _transforms -= transforms;
    timer = (getHighResTimer() - timer)/getHighResTimerFrequency();
    _logging->info("modulus degree %d, transforms: %d, time: %.3f s.\n", _modulus->degree(), transforms, timer);
    commit_setup();
}

template<class Element>
void Stage2Poly<Element>::execute()
{
    int i, j, k;
    double timer;
    double timer_merges = 0;
    int count_merges = 0;
    int degree;
    PolyMult& pm = _poly_mult[poly_power()];
    Poly Lo(pm);
    Poly Hi(pm);
    Poly Quot(pm);

#ifdef DEBUG_STAGE2POLY_QUEUE
    Element Xn(arithmetic());
    Element Xn1(arithmetic());
    for (i = 0; i < _workqueue.size(); i++)
    {
        arithmetic().mul(*_X1, _workqueue[i].n, Xn, Xn1);
        GWASSERT(*_workqueue[i].Xn == Xn);
        GWASSERT(!_workqueue[i].Xn1 || *_workqueue[i].Xn1 == Xn1);
        if (_workqueue[i].distance > (1 << _smallpoly_power))
        {
            arithmetic().mul(*_X1, _workqueue[i].n + (1 << _smallpoly_power), Xn, Xn1);
            GWASSERT(!_workqueue[i].Xdn || *_workqueue[i].Xdn == Xn);
            GWASSERT(!_workqueue[i].Xdn1 || *_workqueue[i].Xdn1 == Xn1);
        }
    }
#endif

    std::vector<gwarray> gwarrays_up(poly_power() + 1);
    gwarray last_gwarray = _gwarrays.back();
    _gwarrays.pop_back();
    for (j = poly_power(); j > 0 && _gwarrays[j] == nullptr; j--);
    gwarray mod_gwarray = _gwarrays[j];
    _gwarrays[j] = nullptr;
    for (j = poly_power(); j >= 0; j--)
    {
        if (j == 0 || &_poly_mult[j].gw() != &_poly_mult[j - 1].gw())
        {
            if (last_gwarray == nullptr)
                gwarrays_up[j] = gwalloc_array(_poly_mult[j].gw().gwdata(), 1 << poly_power());
            else
            {
                gwarrays_up[j] = last_gwarray;
                last_gwarray = nullptr;
            }
            if (j <= _smallpoly_power)
                break;
        }
    }
    std::vector<gwarray> gwarrays_down(poly_power() + 1);
    for (j = poly_power() + 1; j >= _smallpoly_power; j--)
    {
        if (j == poly_power() + 1 || &_poly_mult[j].gw() != &_poly_mult[j - 1].gw())
        {
            for (i = j - 1; i > 0 && gwarrays_up[i] == nullptr; i--);
            if (gwarrays_up[i] == nullptr)
                gwarrays_down[j - 1] = gwalloc_array(_poly_mult[j - 1].gw().gwdata(), 1 << poly_power());
            else
                gwarrays_down[j - 1] = gwarrays_up[i];
            j = i + 1;
        }
    }
    last_gwarray = gwarrays_down[poly_power()];

    std::vector<std::vector<arithmetic::Poly>> poly_prod(poly_power() + 1);
    int count = iterations() - state()->iteration();
    while (count > 0)
    {
        _logging->debug("%d steps left.\n", count);
        std::unique_lock<std::mutex> lock(_work_mutex);
        _smallpoly_alloc = std::move(_smallpoly);
        if (!_accumulator)
            degree = count >= _poly_mod_degree ? _poly_mod_degree : count;
        else
            degree = count >= (1 << poly_power()) ? (1 << poly_power()) : count;
        smallpoly_init(degree);
        for (j = 0; j <= _smallpoly_power && gwarrays_up[j] == nullptr; j++);
        smallpoly_init_level(gwarrays_up[j], j);

        _workstage = 2;
        lock.unlock();
        _workstage_signal.notify_all();
        _workers[0]->run();
        lock.lock();
        _workdone_signal.wait(lock, [&] {return _thread_exception || _smallpoly_count == 0; });
        lock.unlock();
        if (_thread_exception)
            std::rethrow_exception(_thread_exception);
        _logging->debug("all threads done.\n");

        poly_prod[_smallpoly_power].resize(_smallpoly.size(), Poly(_poly_mult[_smallpoly_power]));
        for (i = 0; i < _smallpoly.size(); i++)
            poly_prod[_smallpoly_power][i] = std::move(_smallpoly[i][_smallpoly_power][0]);

        for (j = _smallpoly_power + 1; j <= poly_power(); j++)
        {
            PolyMult* prev = &_poly_mult[j - 1];
            PolyMult* cur = &_poly_mult[j];
            bool convert = &cur->gw() != &prev->gw();
            int options = !convert ? POLYMULT_NEXTFFT : 0;
            Poly tmp(*prev, 0, true);

            poly_prod[j].resize((degree + (1 << j) - 1) >> j, Poly(_poly_mult[j]));
            if (convert)
            {
                GWASSERT(gwarrays_up[j] != nullptr);
                for (i = 0; i < poly_prod[j].size(); i++)
                {
                    GWASSERT(degree > (i << j));
                    _poly_mult[j].init(gwarrays_up[j] + (i << j), (size_t)1 << j, false, true, poly_prod[j][i]);
                }
            }

            for (i = 0; i < poly_prod[j].size(); i++)
            {
                GWASSERT(2*i < poly_prod[j - 1].size());
                if (2*i + 1 == poly_prod[j - 1].size())
                {
                    if (convert)
                        prev->convert(poly_prod[j - 1][2*i], *cur, poly_prod[j][i]);
                    else
                        poly_prod[j][i] = std::move(poly_prod[j - 1][2*i]);
                    continue;
                }
                Poly* res = convert ? &tmp : &poly_prod[j][i];
                if (!_accumulator && j == poly_power())
                {
                    if (poly_prod[j - 1][2*i].size() + poly_prod[j - 1][2*i + 1].size() == _poly_mod_degree)
                    {
                        if (convert)
                        {
                            prev->mul(std::move(poly_prod[j - 1][2*i]), std::move(poly_prod[j - 1][2*i + 1]), tmp, options);
                            prev->convert(tmp, *cur, poly_prod[j][0]);
                            for (i = 0; i < _poly_mod_degree; i++)
                                cur->gw().sub((GWNum&)poly_prod[j][0].at(i), (GWNum&)_poly_mod[poly_power()][0].at(i), (GWNum&)_poly_mod[poly_power()][0].at(i));
                            cur->init(mod_gwarray, _poly_mod_degree, false, false, poly_prod[j][0]);
                            cur->free(_poly_mod[poly_power()][0]);
                            break;
                        }
                        else
                        {
                            poly_prod[j][i] = std::move(_poly_mod[poly_power()][0]);
                            prev->fma_range(poly_prod[j - 1][2*i], poly_prod[j - 1][2*i + 1], poly_prod[j][i], poly_prod[j][i], 0, _poly_mod_degree, options | POLYMULT_FMSUB);
                        }
                    }
                    else
                    {
                        poly_prod[j][i] = std::move(_poly_mod[poly_power()][0]);
                        if (convert)
                            prev->mul(std::move(poly_prod[j - 1][2*i]), std::move(poly_prod[j - 1][2*i + 1]), *res, options);
                        else
                            prev->mul(poly_prod[j - 1][2*i], poly_prod[j - 1][2*i + 1], *res, options);
                    }
                }
                else
                    prev->mul(std::move(poly_prod[j - 1][2*i]), std::move(poly_prod[j - 1][2*i + 1]), *res, options);
                if (convert)
                    prev->convert(tmp, *cur, poly_prod[j][i]);
            }
        }
        count -= (int)poly_prod[poly_power()][0].size();

        if (!_accumulator && poly_prod[poly_power()][0].degree() < _poly_mod_degree)
            _accumulator.reset(new Poly(std::move(poly_prod[poly_power()][0])));
        else
        {
            _logging->debug("merging.\n");
            timer = getHighResTimer();
            if (!_accumulator)
                _accumulator.reset(new Poly(pm, 0, true));
            pm.init(mod_gwarray, _poly_mod_degree, false, false, Lo);
            pm.init(last_gwarray, (size_t)1 << poly_power(), false, false, Hi);

            pm.mul_split(*_accumulator, poly_prod[poly_power()][0], Lo, Hi, _poly_mod_degree, POLYMULT_NEXTFFT);
            if (Hi.degree() >= 0)
            {
                if (Hi.size() > 0)
                {
                    pm.init(last_gwarray, _poly_mod_degree, false, false, Quot);
                    pm.mul_range(Hi, *_reciprocal, Quot, (1 << (poly_power() + 1)) - _poly_mod_degree, _poly_mod_degree, POLYMULT_NEXTFFT);
                }
                else
                    pm.init(true, Quot);
                pm.init(mod_gwarray, _poly_mod_degree, false, false, *_accumulator);
                pm.fma_range(Quot, *_modulus, Lo, *_accumulator, 0, _poly_mod_degree, POLYMULT_FNMADD | POLYMULT_NEXTFFT);
            }
            else
                *_accumulator = std::move(Lo);

            count_merges++;
            timer_merges += (getHighResTimer() - timer)/getHighResTimerFrequency();
        }

        if (_file != nullptr)
            set_state<State>(iterations() - count);
    }
    poly_prod.clear();
    _logging->info("merges: %d, time: %.3f s.\n", count_merges, timer_merges);
    timer = getHighResTimer();

    if (poly_check())
        for (i = 1; i < _workers.size(); i++)
        {
            *_workers[0]->check() *= *_workers[i]->check();
            _workers[i]->check().reset(new GWNumWrapper(*_workers[0]->check()));
        }

#ifdef DEBUG_STAGE2POLY_REM
    if (poly_check())
        GWASSERT(*_workers[0]->check() == _smallpoly_mod[_poly_check_index >> _smallpoly_power][_smallpoly_power][0].at(_poly_check_index & ((1 << _smallpoly_power) - 1)));
#ifdef DEBUG_STAGE2POLY_ROOT
    GWNum root(_poly_mult[poly_power()].gw());
    for (k = 0; k < _smallpoly_mod.size(); k++)
        for (i = 0; i < (_smallpoly_mod[k][0].size() << _tinypoly_power); i++)
        {
            if (&_poly_mult[0].gw() != &_poly_mult[poly_power()].gw())
                gwconvert(_poly_mult[0].gw().gwdata(), _poly_mult[poly_power()].gw().gwdata(), _smallpoly_mod[k][0][i >> _tinypoly_power].data()[i & ((1 << _tinypoly_power) - 1)], *root);
            else
                gwunfft(_poly_mult[0].gw().gwdata(), _smallpoly_mod[k][0][i >> _tinypoly_power].data()[i & ((1 << _tinypoly_power) - 1)], *root);
            root.arithmetic().neg(root, root);
            GWNum val = _accumulator->eval(root);
            if (&_poly_mult[0].gw() != &_poly_mult[poly_power()].gw())
                gwconvert(_poly_mult[0].gw().gwdata(), _poly_mult[poly_power()].gw().gwdata(), _smallpoly_mod[k][_smallpoly_power][0].data()[i], *root);
            else
                gwunfft(_poly_mult[0].gw().gwdata(), _smallpoly_mod[k][_smallpoly_power][0].data()[i], *root);
            GWASSERT(root == val);
        }
#endif
#endif

    std::vector<Poly> poly_rem;
    poly_rem.reserve(_smallpoly_mod.size());
    poly_rem.emplace_back(pm);
    std::vector<Poly> poly_rem_base;
    poly_rem_base.reserve(_smallpoly_mod.size());

    pm.init(last_gwarray, (size_t)1 << poly_power(), false, false, poly_rem[0]);
    pm.mul_range(*_accumulator, *_reciprocal, poly_rem[0], 1 << poly_power(), 1 << poly_power(), gwarrays_up[poly_power()] == nullptr ? POLYMULT_NEXTFFT : 0);
    _modulus.reset();
    _reciprocal.reset();
    _accumulator.reset();
    
    if (gwarrays_up[poly_power()] != nullptr)
    {
        PolyMult& pm_prev = _poly_mult[poly_power() - 1];
        Poly tmp(pm_prev);
        pm_prev.init(gwarrays_down[poly_power() - 1], (size_t)1 << poly_power(), false, false, tmp);
        pm.convert(poly_rem[0], pm_prev, tmp);
        poly_rem[0] = std::move(tmp);
        gwfree_array(pm.gw().gwdata(), last_gwarray);
        gwfree_array(pm.gw().gwdata(), mod_gwarray);
    }
    else if (_gwarrays[poly_power() - 1] == nullptr)
        _gwarrays[poly_power() - 1] = mod_gwarray;
    else
        gwfree_array(pm.gw().gwdata(), mod_gwarray);
    gwarrays_up.resize(poly_power());

    if (_mem_model == -2)
    {
#ifdef DEBUG_STAGE2POLY_REM
        std::vector<Poly> rems;
        for (k = 0; k < _smallpoly_mod.size(); k++)
            rems.emplace_back(std::move(_smallpoly_mod[k][_smallpoly_power][0]));
#endif
        std::unique_lock<std::mutex> lock(_work_mutex);
        _smallpoly_count = _smallpoly_mod.size();
        _smallpoly.clear();
        _smallpoly.resize(_smallpoly_mod.size());
        _smallpoly_alloc = std::move(_smallpoly_mod);
        _gwarrays[1] = gwalloc_array(_poly_mult[1].gw().gwdata(), _poly_mod_degree);
        smallpoly_init_level(_gwarrays[1], 1);

        _workstage = 3;
        lock.unlock();
        _workstage_signal.notify_all();
        _workers[0]->run();
        lock.lock();
        _workdone_signal.wait(lock, [&] {return _thread_exception || _smallpoly_count == 0; });
        lock.unlock();
        if (_thread_exception)
            std::rethrow_exception(_thread_exception);
        _logging->debug("all threads done.\n");

        _smallpoly_mod = std::move(_smallpoly);
        for (k = 0; k < _smallpoly_mod.size(); k++)
            _poly_mod[_smallpoly_power][k] = std::move(_smallpoly_mod[k][_smallpoly_power][0]);
#ifdef DEBUG_STAGE2POLY_REM
        for (k = 0; k < _smallpoly_mod.size(); k++)
            _smallpoly_mod[k][_smallpoly_power][0] = std::move(rems[k]);
#endif
    }

    int cur_level = poly_power() - 1;
    while (cur_level > _smallpoly_power - 1)
    {
        int base_level = cur_level;
        while ((base_level - _smallpoly_power)%_poly_optmem != 0)
            base_level--;
        if (_mem_model == -2 && _poly_mod[base_level][0].empty())
        {
            _poly_mod.resize(base_level + 1);
            build_tree(_poly_mod, _poly_mult, _smallpoly_power, false, [&](int j) { return (j - _smallpoly_power - 1)%_poly_optmem == 0; }, _gwarrays, _poly_mod_degree, _gwarrays_tmp);
        }
        GWASSERT(!_poly_mod[base_level][0].empty());

        std::vector<std::vector<Poly>> tree(cur_level - base_level + 1);
        int base_size = (1 << tree.size());
        int level_size = (1 << (cur_level + 1)) <= _poly_mod_degree ? (1 << (cur_level + 1)) : _poly_mod_degree;
        for (k = 0; k < poly_rem.size(); k++)
        {
            tree[0].resize((k + 1)*base_size <= _poly_mod[base_level].size() ? base_size : _poly_mod[base_level].size() - k*base_size, Poly(_poly_mult[base_level]));
            for (i = 0; i < tree[0].size(); i++)
                tree[0][i] = std::move(_poly_mod[base_level][k*base_size + i]);

            build_tree(tree, _poly_mult, base_level, true, [](int) { return true; }, _gwarrays, level_size, _gwarrays_tmp);

            std::vector<Poly> rem;
            rem.reserve(tree[0].size());
            rem.emplace_back(std::move(poly_rem[k]));
            
            rem_tree(rem, k, tree, _poly_mult, cur_level, base_level, 0, gwarrays_down, 1 << poly_power());

            for (i = 0; i < rem.size(); i++)
                poly_rem_base.emplace_back(std::move(rem[i]));
        }
        for (j = base_level; j <= cur_level; j++)
        {
            if (_gwarrays[j] != nullptr)
                gwfree_array(_poly_mult[j].gw().gwdata(), _gwarrays[j]);
            if (_gwarrays_tmp[j] != nullptr)
                gwfree_array(_poly_mult[j].gw().gwdata(), _gwarrays_tmp[j]);
            if (gwarrays_up[j] != nullptr)
                gwfree_array(_poly_mult[j].gw().gwdata(), gwarrays_up[j]);
        }
        _gwarrays.resize(base_level);
        _gwarrays_tmp.resize(base_level);
        gwarrays_up.resize(base_level);
        _poly_mod.resize(base_level);

        std::swap(poly_rem, poly_rem_base);
        poly_rem_base.clear();
        cur_level = base_level - 1;
    }
    for (i = cur_level; gwarrays_down[i] == nullptr; i++);
    _gwarrays.push_back(gwarrays_down[i]);
    gwarrays_up.clear();
    gwarrays_down.clear();
    _poly_mod.clear();

    {
        std::unique_lock<std::mutex> lock(_work_mutex);
        _smallpoly_count = _smallpoly_mod.size();
        _smallpoly_rem = std::move(poly_rem);
        GWASSERT(_smallpoly_rem.size() == _smallpoly_mod.size());

        _workstage = 4;
        lock.unlock();
        _workstage_signal.notify_all();
        _workers[0]->run();
        lock.lock();
        _workdone_signal.wait(lock, [&] {return _thread_exception || _smallpoly_count == 0; });
        lock.unlock();
        if (_thread_exception)
            std::rethrow_exception(_thread_exception);
        _logging->debug("all threads done.\n");
    }

    if (_gwarrays_tmp[0] != nullptr)
        gwfree_array(gw().gwdata(), _gwarrays_tmp[0]);
    _gwarrays_tmp.clear();
    for (j = 0; j < _smallpoly_power; j++)
        if (_gwarrays[j] != nullptr)
            gwfree_array(_poly_mult[j].gw().gwdata(), _gwarrays[j]);
    gwfree_array(_poly_mult[_smallpoly_power - 1].gw().gwdata(), _gwarrays[_smallpoly_power]);
    _gwarrays.clear();

    timer = (getHighResTimer() - timer)/getHighResTimerFrequency();
    _logging->info("polymod time: %.3f s.\n", timer);

    for (i = 1; i < _workers.size(); i++)
        _workers[0]->G() *= _workers[i]->G();
    Giant tmp;
    tmp = _workers[0]->G();
    tmp %= gw().N();
    _res64 = tmp.to_res64();
    //GWASSERT(_res64 == "5E0857283607657F"); // PolyPower=7 "1023*2^3160+1" -edecm -curve seed 83988 -B1 10k -B2 2000k
    //GWASSERT(_res64 == "5CBB3D1C8B7819DA"); // PolyPower=9 "1023*2^3160+1" -edecm -curve seed 83988 -B1 10k -B2 2000k
    tmp = gcd(std::move(tmp), gw().N());

    done(tmp);
}

template<class Element>
void Stage2Poly<Element>::release()
{
    {
        std::lock_guard<std::mutex> lock(_work_mutex);
        _workstage = 10;
    }
    _workstage_signal.notify_all();
    _workqueue_signal.notify_all();
    while (!_threads.empty())
    {
        _threads.back().join();
        _threads.pop_back();
    }
    _workqueue.clear();
    _workers.clear();

    _smallpoly.clear();
    _smallpoly_alloc.clear();
    _smallpoly_mod.clear();
    _smallpoly_rem.clear();
    _poly_mod.clear();

    _modulus.reset();
    _reciprocal.reset();
    _accumulator.reset();

    _poly_mult.clear();
    _poly_gw.clear();

    if (_file != nullptr)
        _state.reset(::read_state<State>(_file));
}

template<class Element>
Stage2Poly<Element>::SmallPolyWorker::SmallPolyWorker(Stage2Poly<Element>& stage2) : _stage2(stage2), _gwstate(*stage2._gwstate), _gw(_gwstate)
{
    //_gwstate.copy(*stage2._gwstate);
    //stage2._input->setup(_gwstate);
    _G.reset(new GWNum(_gw));
    *_G = 1;
    _check.reset(new GWNum(_gw));
    *_check = 1;

    int power = _stage2._smallpoly_power;
    _poly_mult.reserve(power + 1);
    _poly_mult.emplace_back(_gw, 1);
    for (int i = 1; i <= power; i++)
    {
        if (2*(1 << i) > _poly_mult[i - 1].max_output() || i == DEBUG_STAGE2POLY_CONVERT)
        {
            _poly_gwstate.emplace_back(stage2._poly_gwstate[_poly_gwstate.size()]);
            //_poly_gwstate.emplace_back();
            //_poly_gwstate.back().copy(stage2._poly_gwstate[_poly_gwstate.size() - 1]);
            //stage2._input->setup(_poly_gwstate.back());
            _poly_gw.emplace_back(_poly_gwstate.back());
            _poly_mult.emplace_back(_poly_gw.back(), 1);
        }
        else
            _poly_mult.emplace_back(_poly_mult[i - 1].gw(), 1);
    }
}

template<class Element>
Stage2Poly<Element>::SmallPolyWorker::~SmallPolyWorker()
{
    _poly_mult.clear();
    _poly_gw.clear();
    _poly_gwstate.clear();
    _G.reset();
    _check.reset();
}

template<class Element>
void Stage2Poly<Element>::SmallPolyWorker::run()
{
    int i, j, k;
    int workstage = _stage2._workstage;
    int power = _stage2._smallpoly_power;
    int tiny_power = _stage2._tinypoly_power;
    int degree, index;
    std::vector<std::unique_ptr<Element>> elements;
    elements.reserve((size_t)1 << power);
    int elements_count = 0;
    std::vector<int> element_map;
    element_map.resize((size_t)1 << power);
    std::vector<std::vector<Poly>> poly_tree;
    std::vector<Poly> poly_rem;
    poly_rem.reserve((size_t)1 << (power - tiny_power));
    bool workitem_done = false;
    std::vector<gwarray> gwarrays(power + 1);
    std::vector<gwarray> gwarrays_tmp(power + 1);
    std::vector<gwarray> gwarrays_down(power + 1);
    std::unique_ptr<GWNum> check_root;

    try
    {
        while (workstage <= 2)
        {
            std::unique_ptr<Element> Xn;
            std::unique_ptr<Element> Xn1;
            std::unique_ptr<Element> TXn;
            std::unique_ptr<Element> TXn1;
            std::unique_ptr<Element> TXdn;
            std::unique_ptr<Element> TXdn1;
            int count;
            int n;
            int distance;

            {
                std::unique_lock<std::mutex> lock(_stage2._work_mutex);
                _stage2._workstage_signal.wait(lock, [&] {return (workstage = _stage2._workstage) >= 1; });

                if (!poly_tree.empty() && workitem_done)
                {
                    if (workstage == 1)
                        std::copy(element_map.begin(), element_map.begin() + degree, _stage2._element_map.begin() + (index << power));
                    _stage2._smallpoly[index] = std::move(poly_tree);
                    _stage2._smallpoly_count--;
                    workitem_done = false;
                    if (_stage2._smallpoly_count == 0 && !_main)
                        _stage2._workdone_signal.notify_one();
                }
                auto pred = [&] {return (workstage = _stage2._workstage) > 2 || !poly_tree.empty() || !_stage2._smallpoly_alloc.empty(); };
                if (_main && !pred())
                    break;
                _stage2._workstage_signal.wait(lock, pred);
                if (workstage > 2)
                    break;
                if (poly_tree.empty())
                {
                    poly_tree = std::move(_stage2._smallpoly_alloc.back());
                    _stage2._smallpoly_alloc.pop_back();
                    degree = poly_tree[0].size() << tiny_power;
                    index = _stage2._smallpoly_alloc.size();
                }

                _stage2._workqueue_signal.wait(lock, [&] {return (workstage = _stage2._workstage) > 2 || !_stage2._workqueue.empty(); });
                if (workstage > 2)
                    break;
                Xn = std::move(_stage2._workqueue.back().Xn);
                Xn1 = std::move(_stage2._workqueue.back().Xn1);
                TXdn = std::move(_stage2._workqueue.back().Xdn);
                TXdn1 = std::move(_stage2._workqueue.back().Xdn1);
                count = _stage2._workqueue.back().count;
                n = _stage2._workqueue.back().n;
                distance = _stage2._workqueue.back().distance;
                _stage2._workqueue.pop_back();
            }
            if (workstage == 2 && distance > (1 << power))
            {
                TXn.reset(new Element(arithmetic()));
                TXn1.reset(new Element(arithmetic()));
                *TXn = *Xn;
                *TXn1 = *Xn1;
                std::swap(TXn, Xn);
                std::swap(TXn1, Xn1);

                distance -= (1 << power);
                std::swap(TXn, TXdn);
                std::swap(TXn1, TXdn1);
                if (distance > (1 << power))
                {
                    if (n > 0)
                        arithmetic().add(*_stage2._Xd, *TXn, *Xn, *TXdn);
                    else
                        arithmetic().dbl(*TXn, *TXdn);
                }
                else
                    TXdn.reset();
                if (distance > (1 << power) + 1)
                    arithmetic().add(*_stage2._Xd, *TXn1, *Xn1, *TXdn1);
                else
                    TXdn1.reset();

                std::unique_lock<std::mutex> lock(_stage2._work_mutex);
                _stage2._workqueue.emplace_front();
                _stage2._workqueue.front().Xn = std::move(TXn);
                _stage2._workqueue.front().Xn1 = std::move(TXn1);
                _stage2._workqueue.front().Xdn = std::move(TXdn);
                _stage2._workqueue.front().Xdn1 = std::move(TXdn1);
                _stage2._workqueue.front().count = distance >= (1 << power) ? (1 << power) : distance;
                _stage2._workqueue.front().n = n + (1 << power);
                _stage2._workqueue.front().distance = distance;
                lock.unlock();
                _stage2._workqueue_signal.notify_one();

                distance = 1;
            }
            if (workstage == 2 && _stage2.poly_check() && !check_root)
                check_root.reset(new GWNumWrapper(_stage2._smallpoly_mod[_stage2._poly_check_index >> power][0][(_stage2._poly_check_index & ((1 << power) - 1)) >> tiny_power].at(_stage2._poly_check_index & ((1 << tiny_power) - 1))));
            
            while (elements_count < degree && count > 0)
            {
                if (elements_count >= elements.size())
                    elements.emplace_back(new Element(arithmetic()));
                std::unique_ptr<Element>& Xtmp = elements[elements_count];
                if (count > 2)
                {
                    if (n > 0)
                        arithmetic().add(*_stage2._X1, *Xn1, *Xn, *Xtmp);
                    else
                        arithmetic().dbl(*Xn1, *Xtmp);
                }
                std::swap(Xn, Xn1);
                std::swap(Xn1, Xtmp);
                if (workstage > 1 || gcd(_stage2.D(), n) == 1)
                {
                    element_map[elements_count] = n;
                    elements_count++;
                }
                count--;
                n += distance;
            }
            if (count > 0)
            {
                std::unique_lock<std::mutex> lock(_stage2._work_mutex);
                _stage2._workqueue.emplace_back();
                _stage2._workqueue.back().Xn = std::move(Xn);
                _stage2._workqueue.back().Xn1 = std::move(Xn1);
                _stage2._workqueue.back().count = count;
                _stage2._workqueue.back().n = n;
                _stage2._workqueue.back().distance = distance;
                lock.unlock();
                _stage2._workqueue_signal.notify_one();
                count = 0;
            }

            if (elements_count == degree || workstage == 2)
            {
                if (poly_tree[0][0].empty())
                {
                    if (gwarrays[0] == nullptr)
                        gwarrays[0] = gwalloc_array(_poly_mult[0].gw().gwdata(), 1 << power);
                    for (k = 0; k < poly_tree[0].size(); k++)
                        _poly_mult[0].init(gwarrays[0] + (k << tiny_power), (size_t)1 << tiny_power, false, true, poly_tree[0][k]);
                }
                if (elements_count < degree)
                {
                    for (j = 0; j <= tiny_power; j++)
                    {
                        poly_tree[j].erase(poly_tree[j].begin() + ((elements_count + ((size_t)1 << tiny_power) - 1) >> tiny_power), poly_tree[j].end());
                        if ((elements_count >> tiny_power) < poly_tree[j].size() && !poly_tree[j][elements_count >> tiny_power].empty())
                            _poly_mult[j].alloc(poly_tree[j][elements_count >> tiny_power], elements_count & ((1 << tiny_power) - 1));
                    }
                }
                std::vector<gwnum> roots(elements_count);
                for (i = 0; i < roots.size(); i++)
                    roots[i] = poly_tree[0][i >> tiny_power].data()[i & ((1 << tiny_power) - 1)];
                elements_to_gwnums(elements, elements_count, roots.data());
                elements_count = 0;
                
                if (tiny_power > 0)
                    for (k = 0; k < poly_tree[0].size(); k++)
                        build_tree_tinypoly(poly_tree, k, tiny_power, _poly_mult, 0, [&](int j) { return workstage == 1 && (j - 1)%_stage2._poly_optmem_small == 0; }, gwarrays, 1 << power, gwarrays_tmp, check_root ? check_root.get() : nullptr, check_root ? _check.get() : nullptr);
                if (tiny_power < power)
                    build_tree(poly_tree, _poly_mult, tiny_power, false, [&](int j) { return workstage == 1 && (j - 1)%_stage2._poly_optmem_small == 0; }, gwarrays, 1 << power, gwarrays_tmp);
                workitem_done = true;
            }
        }
        elements.clear();

        workitem_done = false;
        while (workstage == 3)
        {
            {
                std::unique_lock<std::mutex> lock(_stage2._work_mutex);

                if (workitem_done)
                {
                    _stage2._smallpoly[index] = std::move(poly_tree);
                    _stage2._smallpoly_count--;
                    workitem_done = false;
                    if (_stage2._smallpoly_count == 0 && !_main)
                        _stage2._workdone_signal.notify_one();
                }
                auto pred = [&] {return (workstage = _stage2._workstage) > 3 || !_stage2._smallpoly_alloc.empty(); };
                if (_main && !pred())
                    break;
                _stage2._workstage_signal.wait(lock, pred);
                if (workstage > 3)
                    break;
                poly_tree = std::move(_stage2._smallpoly_alloc.back());
                _stage2._smallpoly_alloc.pop_back();
                degree = poly_tree[0].size() << tiny_power;
                index = _stage2._smallpoly_alloc.size();
            }

            if (tiny_power > 0)
                for (k = 0; k < poly_tree[0].size(); k++)
                    build_tree_tinypoly(poly_tree, k, tiny_power, _poly_mult, 0, [&](int j) { return j == 1; }, gwarrays, 1 << power, gwarrays_tmp, nullptr, nullptr);
            if (tiny_power < power)
                build_tree(poly_tree, _poly_mult, tiny_power, false, [&](int j) { return j == 1; }, gwarrays, 1 << power, gwarrays_tmp);
            workitem_done = true;
        }

        check_root.reset();
        workitem_done = false;
        while (workstage == 4)
        {
            {
                std::unique_lock<std::mutex> lock(_stage2._work_mutex);

                if (workitem_done)
                {
                    _stage2._smallpoly_count--;
                    workitem_done = false;
                }
                if (_stage2._smallpoly_rem.empty())
                {
                    if (_stage2._smallpoly_count == 0 && !_main)
                        _stage2._workdone_signal.notify_one();
                    break;
                }
                if ((workstage = _stage2._workstage) > 4)
                    break;
                poly_tree = std::move(_stage2._smallpoly_mod.back());
                _stage2._smallpoly_mod.pop_back();
                poly_rem.clear();
                poly_rem.emplace_back(std::move(_stage2._smallpoly_rem.back()));
                _stage2._smallpoly_rem.pop_back();
                gwarrays_down[power - 1] = _stage2._gwarrays.back() + (_stage2._smallpoly_rem.size() << power);
                if (_stage2.poly_check() && _stage2._smallpoly_rem.size() == (_stage2._poly_check_index >> power))
                    check_root.reset(new GWNumWrapper(*_check));
            }

#ifdef DEBUG_STAGE2POLY_REM
            Poly debug_poly_rem = std::move(poly_tree[power][0]);
#endif

            poly_tree.resize(power);
            if (tiny_power > 0)
                for (k = 0; k < poly_tree[0].size(); k++)
                    build_tree_tinypoly(poly_tree, k, tiny_power, _poly_mult, 0, [](int) { return true; }, gwarrays, 1 << power, gwarrays_tmp, nullptr, nullptr);
            if (tiny_power < power - 1)
                build_tree(poly_tree, _poly_mult, tiny_power, false, [](int) { return true; }, gwarrays, 1 << power, gwarrays_tmp);
            if (tiny_power < power)
                rem_tree(poly_rem, 0, poly_tree, _poly_mult, power - 1, 0, tiny_power, gwarrays_down, 1 << power);
            if (tiny_power > 0)
                for (k = 0; k < poly_rem.size(); k++)
                    rem_tree_tinypoly(poly_rem, k, poly_tree, tiny_power, _poly_mult, gwarrays_down, 1 << power);

#ifdef DEBUG_STAGE2POLY_REM
            GWASSERT(debug_poly_rem.size() == (poly_rem.size() << tiny_power));
            for (i = 0; i < debug_poly_rem.size(); i++)
                GWASSERT(gw().cmp((GWNum&)debug_poly_rem.at(i), (GWNum&)poly_rem[i >> tiny_power].at(i & ((1 << tiny_power) - 1))) == 0);
#endif
            if (check_root)
            {
                if ((gw().popg() = *check_root)%gw().N() != (gw().popg() = poly_rem[(_stage2._poly_check_index & ((1 << power) - 1)) >> tiny_power].at(_stage2._poly_check_index & ((1 << tiny_power) - 1)))%gw().N())
                {
                    _stage2._logging->error("check failed.\n");
#ifdef DEBUG_STAGE2POLY_REM
                    GWASSERT(0);
#endif
                    throw TaskAbortException();
                }
                else
                    _stage2._logging->info("check passed.\n");
                check_root.reset();
            }

            for (k = 0; k < poly_rem.size(); k++)
                for (i = 0; i < poly_rem[k].size(); i++)
                    gw().mul((GWNum&)poly_rem[k].at(i), G(), G(), GWMUL_STARTNEXTFFT);

            workitem_done = true;
        }

        gwarrays_down[power - 1] = nullptr;
        for (j = 0; j < gwarrays.size(); j++)
            if (gwarrays[j] != nullptr)
                gwfree_array(_poly_mult[j].gw().gwdata(), gwarrays[j]);
        gwarrays.clear();
        for (j = 0; j < gwarrays_tmp.size(); j++)
            if (gwarrays_tmp[j] != nullptr)
                gwfree_array(_poly_mult[j].gw().gwdata(), gwarrays_tmp[j]);
        gwarrays_tmp.clear();
        for (j = 0; j < gwarrays_down.size(); j++)
            if (gwarrays_down[j] != nullptr)
                gwfree_array(_poly_mult[j].gw().gwdata(), gwarrays_down[j]);
        gwarrays_down.clear();
        GWASSERT(gw().gwdata()->array_list == nullptr);
    }
    catch (const std::exception&)
    {
        if (_main)
            throw;
        {
            std::unique_lock<std::mutex> lock(_stage2._work_mutex);
            _stage2._workstage = 11;
            _stage2._thread_exception = std::current_exception();
        }
        _stage2._workstage_signal.notify_all();
        _stage2._workqueue_signal.notify_all();
        _stage2._workdone_signal.notify_one();
    }
}

template<class Element>
void Stage2Poly<Element>::write_state(Writer* writer)
{
    int i, j, k;
    Giant tmp;
    if (poly_check())
    {
        GWNum check(gw());
        check = *_workers[0]->check();
        for (i = 1; i < _workers.size(); i++)
            check *= *_workers[i]->check();
        writer->write(_element_map[_poly_check_index]);
        writer->write(tmp = check);
    }
    else
    {
        writer->write(0);
        writer->write(0);
    }

    int count = (1 << _smallpoly_power);
    int files = ((int)_accumulator->size() + count - 1)/count;
    writer->write(_accumulator->monic() ? -files : files);
    for (i = 0, k = 0; k < files; k++)
    {
        File* subf = _file->add_child(std::to_string(k), File::unique_fingerprint(_file->fingerprint(), std::to_string(state()->iteration())));
        std::unique_ptr<Writer> subwriter(subf->get_writer(State::SUB_TYPE, 0));
        count = (1 << _smallpoly_power);
        if (i + count > _accumulator->size())
            count = _accumulator->size() - i;
        subwriter->write(count);
        for (j = 0; j < count; j++, i++)
            subwriter->write(tmp = _accumulator->at(i));
        subf->commit_writer(*subwriter);
        subf->free_buffer();
    }
}

template<class Element>
bool Stage2Poly<Element>::read_state(Reader* reader)
{
    int i, j, k;
    Giant tmp;
    if (!reader->read(j))
        return false;
    if (poly_check() && j != 0)
        for (_poly_check_index = 0; _poly_check_index < _element_map.size() && _element_map[_poly_check_index] != j; _poly_check_index++);
    if (!reader->read(tmp))
        return false;
    if (poly_check() && tmp != 0)
        *_workers[0]->check() = tmp;
    _logging->info("%d, %x\n", _poly_check_index, _element_map[_poly_check_index]);
    //execute();
    
    bool monic = false;
    int files = 0;
    if (!reader->read(files))
        return false;
    if (files <= 0)
    {
        monic = true;
        files *= -1;
    }
    PolyMult& pm = _poly_mult[poly_power()];
    Poly& res = _poly_mod[poly_power()][0];
    for (i = 0, k = 0; k < files; k++)
    {
        File* subf = _file->add_child(std::to_string(k), File::unique_fingerprint(_file->fingerprint(), std::to_string(state()->iteration())));
        std::unique_ptr<Reader> subreader(subf->get_reader());
        if (!subreader)
            break;
        int count;
        if (!subreader->read(count))
            break;
        for (j = 0; j < count; j++, i++)
        {
            if (!subreader->read(tmp))
                break;
            (GWNum&)res.at(i) = tmp;
        }
        if (j < count)
            break;
        subf->free_buffer();
    }
    if (k < files)
    {
        if (i == 0)
            return false;
        Poly one(pm, 0, true);
        pm.mul(*_modulus, one, res, POLYMULT_NEXTFFT);
        return false;
    }
    _accumulator.reset(new Poly(pm));
    pm.init(res.data(), i, false, monic, *_accumulator);
    pm.free(res);
    return true;
}


void PP1Stage2Poly::init(InputNum* input, GWState* gwstate, File* file, Logging* logging, Giant& P, bool minus1)
{
    logging->set_prefix(input->display_text() + (minus1 ? ", P-1 stage 2, " : ", P+1 stage 2, "));
    Stage2Poly::init(input, gwstate, file, ::read_state<State>(file), logging);
    _lucas.reset(new LucasVArithmetic());
    _P = P;
}

template void Stage2Poly<LucasV>::done(const arithmetic::Giant& factor);

void PP1Stage2Poly::SmallPolyWorker::elements_to_gwnums(std::vector<std::unique_ptr<arithmetic::LucasV>>& elements, int count, gwnum* data)
{
    for (int i = 0; i < count; i++)
        gw().unfft(elements[i]->V(), (GWNum&)gw().wrap(data[i]));
}

void PP1Stage2Poly::write_state()
{
    if (_file == nullptr || state() == nullptr || state()->iteration() == 0)
        return;
    _logging->debug("saving state to disk.\n");
    Giant tmp;
    std::unique_ptr<Writer> writer(_file->get_writer(state()->type(), state()->version()));
    auto write_LucasV = [&](std::unique_ptr<LucasV>& p)
    {
        if (p)
            writer->write(tmp = p->V());
        else
            writer->write(0);
    };

    writer->write(state()->iteration());
    writer->write((int)_workqueue.size());
    for (auto it = _workqueue.begin(); it != _workqueue.end(); it++)
    {
        writer->write(it->n);
        writer->write(it->count);
        writer->write(it->distance);
        write_LucasV(it->Xn);
        write_LucasV(it->Xn1);
        write_LucasV(it->Xdn);
        write_LucasV(it->Xdn1);
    }

    Stage2Poly::write_state(writer.get());
    _file->commit_writer(*writer);
    state()->set_written();
    _logging->debug("state saved.\n");
}

void PP1Stage2Poly::setup()
{
    if (!_lucas)
        _lucas.reset(new LucasVArithmetic());
    if (!_safe_gw)
        _safe_gw.reset(new ThreadSafeGWArithmetic(*_gwstate));
    _lucas->set_gw(*_safe_gw);
    if (!_W)
        _W.reset(new LucasV(arithmetic()));
    if (!_Wd)
        _Wd.reset(new LucasV(arithmetic()));

    _workstage = 0;
    _workers.reserve(_poly_threads);
    _workers.emplace_back(new SmallPolyWorker(*this));
    _workers[0]->set_main(true);
    for (int i = 1; i < _poly_threads; i++)
    {
        _workers.emplace_back(new SmallPolyWorker(*this));
        _threads.emplace_back(&SmallPolyWorker::run, _workers.back().get());
    }

    _gwstate->gwdata()->gwnum_max_free_count = ((1 << _smallpoly_power) + 4) + 16;

    if (!_modulus)
    {
        _W->V() = _P;
        gw().fft(_W->V(), _W->V());
        _X1 = _W.get();
        Stage2Poly::setup();
    }

    int v;
    int count = iterations();
    int distance = (count >> _smallpoly_power)/_poly_threads;
    if (distance == 0)
        distance = 1;

    _W->V() = _P;
    arithmetic().mul(*_W, _D, *_W);
    gw().fft(_W->V(), _W->V());
    _X1 = _W.get();
    *_Wd = *_W;
    for (int i = 0; i < _smallpoly_power; i++)
        arithmetic().dbl(*_Wd, *_Wd);
    gw().fft(_Wd->V(), _Wd->V());
    _Xd = _Wd.get();
    std::unique_ptr<LucasV> Wdist;

    if (state() != nullptr)
    {
        std::unique_ptr<Reader> reader(_file->get_reader());
        reader->read(v);
        reader->read(count);
        Giant tmp;
        for (; count > 0; count--)
        {
            _workqueue.emplace_back();
            reader->read(_workqueue.back().n);
            reader->read(_workqueue.back().count);
            reader->read(_workqueue.back().distance);
            reader->read(tmp);
            if (tmp != 0)
            {
                _workqueue.back().Xn.reset(new LucasV(arithmetic()));
                _workqueue.back().Xn->V() = tmp;
            }
            reader->read(tmp);
            if (tmp != 0)
            {
                _workqueue.back().Xn1.reset(new LucasV(arithmetic()));
                _workqueue.back().Xn1->V() = tmp;
            }
            reader->read(tmp);
            if (tmp != 0)
            {
                _workqueue.back().Xdn.reset(new LucasV(arithmetic()));
                _workqueue.back().Xdn->V() = tmp;
            }
            reader->read(tmp);
            if (tmp != 0)
            {
                _workqueue.back().Xdn1.reset(new LucasV(arithmetic()));
                _workqueue.back().Xdn1->V() = tmp;
            }
        }
        read_state(reader.get());
    }

    if (state() == nullptr)
    {
        reset_state<State>();
        v = _first_D;

        _workqueue.emplace_back();
        _workqueue.back().n = v;
        _workqueue.back().count = count >= (1 << _smallpoly_power) ? (1 << _smallpoly_power) : count;
        _workqueue.back().distance = count >= (distance << _smallpoly_power) ? (distance << _smallpoly_power) : count;

        LucasV Pn(arithmetic());
        LucasV Pn1(arithmetic());
        if (v > 1)
            arithmetic().mul(*_W, v, Pn, Pn1);
        else if (v > 0)
        {
            Pn = *_W;
            arithmetic().dbl(Pn, Pn1);
        }
        else
        {
            arithmetic().init(Pn);
            Pn1 = *_W;
        }
        _workqueue.back().Xn.reset(new LucasV(Pn));
        if (_workqueue.back().count > 1)
            _workqueue.back().Xn1.reset(new LucasV(Pn1));

        if (_workqueue.back().distance > (1 << _smallpoly_power))
        {
            arithmetic().mul(*_W, v + (1 << _smallpoly_power), Pn, Pn1);
            _workqueue.back().Xdn.reset(new LucasV(Pn));
            if (_workqueue.back().distance > (1 << _smallpoly_power) + 1)
                _workqueue.back().Xdn1.reset(new LucasV(Pn1));
        }

        int total = _workqueue.back().distance;
        while (total < count)
        {
            _workqueue.emplace_back();
            _workqueue.back().n = v + total;
            _workqueue.back().count = count - total >= (1 << _smallpoly_power) ? (1 << _smallpoly_power) : count - total;
            _workqueue.back().distance = count - total >= (distance << _smallpoly_power) ? (distance << _smallpoly_power) : count - total;
            total += _workqueue.back().distance;

            if (_workqueue.size() > 2)
            {
                if (!Wdist)
                {
                    Wdist.reset(new LucasV(arithmetic()));
                    arithmetic().mul(*_Wd, distance, *Wdist);
                }
                _workqueue.back().Xn.reset(new LucasV(arithmetic()));
                arithmetic().add(*Wdist, *_workqueue[_workqueue.size() - 2].Xn, *_workqueue[_workqueue.size() - 3].Xn, *_workqueue.back().Xn);
                if (_workqueue.back().count > 1)
                {
                    _workqueue.back().Xn1.reset(new LucasV(arithmetic()));
                    arithmetic().add(*Wdist, *_workqueue[_workqueue.size() - 2].Xn1, *_workqueue[_workqueue.size() - 3].Xn1, *_workqueue.back().Xn1);
                }

                if (_workqueue.back().distance > (1 << _smallpoly_power))
                {
                    _workqueue.back().Xdn.reset(new LucasV(arithmetic()));
                    arithmetic().add(*Wdist, *_workqueue[_workqueue.size() - 2].Xdn, *_workqueue[_workqueue.size() - 3].Xdn, *_workqueue.back().Xdn);
                    if (_workqueue.back().distance > (1 << _smallpoly_power) + 1)
                    {
                        _workqueue.back().Xdn1.reset(new LucasV(arithmetic()));
                        arithmetic().add(*Wdist, *_workqueue[_workqueue.size() - 2].Xdn1, *_workqueue[_workqueue.size() - 3].Xdn1, *_workqueue.back().Xdn1);
                    }
                }
            }
            else
            {
                arithmetic().mul(*_W, _workqueue.back().n, Pn, Pn1);
                _workqueue.back().Xn.reset(new LucasV(Pn));
                if (_workqueue.back().count > 1)
                    _workqueue.back().Xn1.reset(new LucasV(Pn1));

                if (_workqueue.back().distance > (1 << _smallpoly_power))
                {
                    arithmetic().mul(*_W, _workqueue.back().n + (1 << _smallpoly_power), Pn, Pn1);
                    _workqueue.back().Xdn.reset(new LucasV(Pn));
                    if (_workqueue.back().distance > (1 << _smallpoly_power) + 1)
                        _workqueue.back().Xdn1.reset(new LucasV(Pn1));
                }
            }
        }

        commit_setup();
    }
}

void PP1Stage2Poly::release()
{
    Stage2Poly::release();
    _Wd.reset();
    _W.reset();
    _lucas.reset();
    _safe_gw.reset();
}

template void Stage2Poly<LucasV>::execute();

void EdECMStage2Poly::init(InputNum* input, GWState* gwstate, File* file, Logging* logging, arithmetic::Giant& X, arithmetic::Giant& Y, arithmetic::Giant& Z, arithmetic::Giant& T, arithmetic::Giant& EdD)
{
    logging->set_prefix(input->display_text() + ", EdECM stage 2, ");
    Stage2Poly::init(input, gwstate, file, ::read_state<State>(file), logging);
    _X = X;
    _Y = Y;
    _Z = Z;
    _T = T;
    _EdD = EdD;
}

template void Stage2Poly<EdY>::done(const arithmetic::Giant& factor);

void EdECMStage2Poly::SmallPolyWorker::elements_to_gwnums(std::vector<std::unique_ptr<arithmetic::EdY>>& elements, int count, gwnum* data)
{
    int i;
    for (i = 0; i < count; i++)
    {
        if (!elements[i]->Z)
        {
            elements[i]->Z.reset(new GWNum(gw()));
            *elements[i]->Z = 1;
        }
        if (!elements[i]->ZpY)
            elements[i]->ZpY.reset(new GWNum(gw()));
    }
    std::swap(elements[0]->ZpY, elements[0]->Z);
    for (i = 1; i < count; i++)
        gw().mul(*elements[i - 1]->ZpY, *elements[i]->Z, *elements[i]->ZpY, GWMUL_STARTNEXTFFT_IF(i + 1 < count));
    try
    {
#ifdef DEBUG_STAGE2POLY_INV
        GWNum test(gw());
        test = *elements[count - 1]->ZpY;
#endif
        gw().inv(*elements[count - 1]->ZpY, *elements[count - 1]->ZpY);
#ifdef DEBUG_STAGE2POLY_INV
        gw().carefully().mul(*elements[count - 1]->ZpY, test, test);
        if ((gw().popg() = test)%gw().N() != 1)
        {
            static_cast<EdECMStage2Poly&>(_stage2)._logging->error("invg() failure.\n");
            throw TaskAbortException();
        }
#endif
    }
    catch (const ArithmeticException&)
    {
        std::swap(elements[0]->ZpY, elements[0]->Z);
        throw;
    }
    for (i = count - 1; i >= 0; i--)
    {
        if (i > 0)
        {
            gw().mul(*elements[i]->ZpY, *elements[i - 1]->ZpY, *elements[i - 1]->ZpY, GWMUL_STARTNEXTFFT);
            std::swap(elements[i]->ZpY, elements[i - 1]->ZpY);
            gw().mul(*elements[i]->Z, *elements[i - 1]->ZpY, *elements[i - 1]->ZpY, GWMUL_STARTNEXTFFT);
        }
        if (elements[i]->Y)
            gw().mul(*elements[i]->ZpY, *elements[i]->Y, (GWNum&)gw().wrap(data[i]), 0);
    }
}

void EdECMStage2Poly::write_state()
{
    if (_file == nullptr || state() == nullptr || state()->iteration() == 0)
        return;
    _logging->debug("saving state to disk.\n");
    Giant tmp;
    std::unique_ptr<Writer> writer(_file->get_writer(state()->type(), state()->version()));
    auto write_EdY = [&](std::unique_ptr<EdY>& p)
    {
        if (p)
        {
            writer->write((tmp = *p->Y));
            if (p->Z)
                writer->write((tmp = *p->Z));
            else
            {
                writer->write(1);
                writer->write(1);
            }
        }
        else
        {
            writer->write(0);
            writer->write(0);
        }
    };

    writer->write(state()->iteration());
    writer->write((int)_workqueue.size());
    for (auto it = _workqueue.begin(); it != _workqueue.end(); it++)
    {
        writer->write(it->n);
        writer->write(it->count);
        writer->write(it->distance);
        write_EdY(it->Xn);
        write_EdY(it->Xn1);
        write_EdY(it->Xdn);
        write_EdY(it->Xdn1);
    }

    Stage2Poly::write_state(writer.get());
    _file->commit_writer(*writer);
    state()->set_written();
    _logging->debug("state saved.\n");
}

/*EdY* to_EdY(MontgomeryArithmetic& arithmetic, EdPoint& a)
{
    std::unique_ptr<EdY> res(new EdY(arithmetic, a));
    if (!res->Z)
    {
        res->Z.reset(new GWNum(arithmetic.gw()));
        *res->Z = 1;
    }
    arithmetic.optimize(*res);
    return res.release();
}*/
#define to_EdY new EdY

void EdECMStage2Poly::setup()
{
    if (!_ed_d)
    {
        _ed_d.reset(new GWNum(gw()));
        *_ed_d = _EdD;
    }
    if (!_montgomery)
        _montgomery.reset(new MontgomeryArithmetic(*_ed_d));
    if (!_safe_gw)
        _safe_gw.reset(new ThreadSafeGWArithmetic(*_gwstate));
    _montgomery->set_gw(*_safe_gw);
    if (!_W)
        _W.reset(new EdY(arithmetic()));
    if (!_Wd)
        _Wd.reset(new EdY(arithmetic()));

    Giant tmp;
    EdwardsArithmetic ed(gw());
    EdPoint EdP(ed);
    EdPoint EdP1(ed);
    EdPoint EdW(ed);
    EdPoint EdWd(ed);
    EdPoint EdWdist(ed);
    int v, count, distance;

    _workstage = 0;
    _workers.reserve(_poly_threads);
    _workers.emplace_back(new SmallPolyWorker(*this));
    _workers[0]->set_main(true);
    for (int i = 1; i < _poly_threads; i++)
    {
        _workers.emplace_back(new SmallPolyWorker(*this));
        _threads.emplace_back(&SmallPolyWorker::run, _workers.back().get());
    }

    _gwstate->gwdata()->gwnum_max_free_count = 4*((1 << _smallpoly_power) + 4) + 16;

    try
    {
        if (!_modulus)
        {
            _W->deserialize(_Y, _Z);
            _X1 = _W.get();
            Stage2Poly::setup();
        }

        count = iterations();
        distance = (count >> _smallpoly_power)/_poly_threads;
        if (distance == 0)
            distance = 1;

        EdP.deserialize(_X, _Y, _Z, _T);
        tmp = _D;
        ed.mul(EdP, tmp, EdW);
        tmp = (1 << _smallpoly_power);
        ed.mul(EdW, tmp, EdWd);
        tmp = distance;
        ed.mul(EdWd, tmp, EdWdist);
        std::vector<EdPoint*> to_norm;
        to_norm.push_back(&EdW);
        to_norm.push_back(&EdWd);
        to_norm.push_back(&EdWdist);
        ed.normalize(to_norm.begin(), to_norm.end(), GWMUL_STARTNEXTFFT);
        *_W = EdW;
        *_Wd = EdWd;
        _X1 = _W.get();
        _Xd = _Wd.get();
    }
    catch (const NoInverseException& e)
    {
        done(e.divisor);
        return;
    }

    if (state() != nullptr)
    {
        auto read_files = [&]()
        {
            std::unique_ptr<Reader> reader(_file->get_reader());
            if (!reader)
                return false;
            if (!reader->read(v))
                return false;
            if (!reader->read(count))
                return false;
            Giant Y, Z;
            for (; count > 0; count--)
            {
                _workqueue.emplace_back();
                reader->read(_workqueue.back().n);
                reader->read(_workqueue.back().count);
                reader->read(_workqueue.back().distance);
                if (!reader->read(Y))
                    return false;
                if (!reader->read(Z))
                    return false;
                if (Y != 0)
                {
                    _workqueue.back().Xn.reset(new EdY(arithmetic()));
                    _workqueue.back().Xn->deserialize(Y, Z);
                }
                if (!reader->read(Y))
                    return false;
                if (!reader->read(Z))
                    return false;
                if (Y != 0)
                {
                    _workqueue.back().Xn1.reset(new EdY(arithmetic()));
                    _workqueue.back().Xn1->deserialize(Y, Z);
                }
                if (!reader->read(Y))
                    return false;
                if (!reader->read(Z))
                    return false;
                if (Y != 0)
                {
                    _workqueue.back().Xdn.reset(new EdY(arithmetic()));
                    _workqueue.back().Xdn->deserialize(Y, Z);
                }
                if (!reader->read(Y))
                    return false;
                if (!reader->read(Z))
                    return false;
                if (Y != 0)
                {
                    _workqueue.back().Xdn1.reset(new EdY(arithmetic()));
                    _workqueue.back().Xdn1->deserialize(Y, Z);
                }
            }
            return read_state(reader.get());
        };
        if (!read_files())
        {
            _logging->warning("restart failed.\n");
            _state.reset();
            _workqueue.clear();
            _accumulator.reset();
            if (poly_check())
                *_workers[0]->check() = 1;
            count = iterations();
        }
    }

    if (state() == nullptr)
    {
        reset_state<State>();
        v = _first_D;

        _workqueue.emplace_back();
        _workqueue.back().n = v;
        _workqueue.back().count = count >= (1 << _smallpoly_power) ? (1 << _smallpoly_power) : count;
        _workqueue.back().distance = count >= (distance << _smallpoly_power) ? (distance << _smallpoly_power) : count;

        *EdP.X = 0;
        *EdP.Y = 1;
        EdP.Z.reset();
        if (v > 1)
        {
            if (v == (1 << _smallpoly_power))
                EdP = EdWd;
            else if (v == (distance << _smallpoly_power))
                EdP = EdWdist;
            else if (v == 2)
                ed.dbl(EdW, EdP);
            else
            {
                tmp = v;
                ed.mul(EdW, tmp, EdP);
            }
            ed.add(EdP, EdW, EdP1, ed.ED_PROJECTIVE);
        }
        else if (v > 0)
        {
            EdP = EdW;
            ed.dbl(EdP, EdP1, ed.ED_PROJECTIVE);
        }
        else
            EdP1 = EdW;
        _workqueue.back().Xn.reset(to_EdY(arithmetic(), EdP));
        if (_workqueue.back().count > 1)
            _workqueue.back().Xn1.reset(to_EdY(arithmetic(), EdP1));

        if (_workqueue.back().distance > (1 << _smallpoly_power))
        {
            if (v == 0)
                EdP1 = EdWd;
            else if (v == ((distance - 1) << _smallpoly_power))
                EdP1 = EdWdist;
            else if (v == (1 << _smallpoly_power))
                ed.dbl(EdWd, EdP1);
            else
                ed.add(EdP, EdWd, EdP1);
            _workqueue.back().Xdn.reset(to_EdY(arithmetic(), EdP1));

            if (_workqueue.back().distance > (1 << _smallpoly_power) + 1)
            {
                ed.add(EdP1, EdW, EdP1, ed.ED_PROJECTIVE);
                _workqueue.back().Xdn1.reset(to_EdY(arithmetic(), EdP1));
            }
        }

        int total = _workqueue.back().distance;
        while (total < count)
        {
            _workqueue.emplace_back();
            _workqueue.back().n = v + total;
            _workqueue.back().count = count - total >= (1 << _smallpoly_power) ? (1 << _smallpoly_power) : count - total;
            _workqueue.back().distance = count - total >= (distance << _smallpoly_power) ? (distance << _smallpoly_power) : count - total;
            total += _workqueue.back().distance;

            if (_workqueue.back().n == (1 << _smallpoly_power))
                EdP = EdWd;
            else if (_workqueue.back().n == (distance << _smallpoly_power))
                EdP = EdWdist;
            else if (_workqueue.back().n == 2*(distance << _smallpoly_power))
                ed.dbl(EdWdist, EdP);
            else
                ed.add(EdP, EdWdist, EdP);
            _workqueue.back().Xn.reset(to_EdY(arithmetic(), EdP));

            if (_workqueue.back().count > 1)
            {
                ed.add(EdP, EdW, EdP1, ed.ED_PROJECTIVE);
                _workqueue.back().Xn1.reset(to_EdY(arithmetic(), EdP1));
            }

            if (_workqueue.back().distance > (1 << _smallpoly_power))
            {
                if (_workqueue.back().n == (1 << _smallpoly_power))
                    ed.dbl(EdWd, EdP1);
                else
                    ed.add(EdP, EdWd, EdP1);
                _workqueue.back().Xdn.reset(to_EdY(arithmetic(), EdP1));

                if (_workqueue.back().distance > (1 << _smallpoly_power) + 1)
                {
                    ed.add(EdP1, EdW, EdP1, ed.ED_PROJECTIVE);
                    _workqueue.back().Xdn1.reset(to_EdY(arithmetic(), EdP1));
                }
            }
        }

#ifdef DEBUG_STAGE2POLY_QUEUE
        GWASSERT(total == count);
        EdY P1(arithmetic());
        P1.deserialize(_Y, _Z);
        arithmetic().mul(P1, _D, P1);
        GWASSERT(P1 == *_X1);
        P1.deserialize(_Y, _Z);
        arithmetic().mul(P1, _D << _smallpoly_power, P1);
        GWASSERT(P1 == *_Xd);
#endif

        commit_setup();
    }
}

void EdECMStage2Poly::release()
{
    Stage2Poly::release();
    _W.reset();
    _Wd.reset();
    _montgomery.reset();
    _safe_gw.reset();
    _ed_d.reset();
}

void EdECMStage2Poly::execute()
{
    if (success()) // NoInverseException in setup()
        return;

    try
    {
        Stage2Poly::execute();
    }
    catch (const NoInverseException& e)
    {
        done(e.divisor);
    }
}
