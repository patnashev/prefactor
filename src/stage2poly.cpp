
#include <deque>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <string.h>

#include "gwnum.h"
#include "cpuid.h"
#include "stage2poly.h"
#include "exception.h"
#include "windows.h"

using namespace arithmetic;

// http://cr.yp.to/arith/scaledmod-20040820.pdf

#ifdef _DEBUG
#define DEBUG_STAGE2POLY_MOD
#define DEBUG_STAGE2POLY_REM
#define DEBUG_STAGE2POLY_ROOT
#define DEBUG_STAGE2POLY_QUEUE
#define DEBUG_STAGE2POLY_CONVERT 0
#else
#define DEBUG_STAGE2POLY_CONVERT 0
#endif

template<class Element>
void Stage2Poly<Element>::init(InputNum* input, GWState* gwstate, File* file, TaskState* state, Logging* logging)
{
#ifdef _DEBUG
    //_tinypoly_power = 0;
#endif
    if (!gwstate->polymult)
    {
        gwstate->polymult = true;
        gwstate->done();
        input->setup(*gwstate);
    }
    if (gwstate->max_polymult_output() < 2*(1 << poly_power()))
    {
        gwstate->done();
        gwstate->next_fft_count++;
        input->setup(*gwstate);
        std::string prefix = logging->prefix();
        logging->set_prefix("");
        logging->warning("Switching to %s\n", gwstate->fft_description.data());
        logging->set_prefix(prefix);
    }
    _poly_gwstate.clear();
    GWState* cur = gwstate;
    for (int i = 1; i <= poly_power(); i++)
    {
        if (2*(1 << i) > cur->max_polymult_output() || i == DEBUG_STAGE2POLY_CONVERT)
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
    //if (state() != nullptr)
    //    logging->info(", restarting at %.1f%%", 100.0*state()->iteration()/iterations());
    logging->info(".\n");
    logging->info("polynomial mode, D = %d, degree %d, %d steps in batches of %d.\n", _D, 1 << poly_power(), _last_D - _first_D + 1, 1 << _smallpoly_power);
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
void Stage2Poly<Element>::smallpoly_init(int count, gwarray data, int data_level)
{
    int i, j, k;
    _smallpoly_count = count;
    _smallpoly_alloc.resize(_smallpoly_count);
    for (k = 0; k < _smallpoly_alloc.size(); k++)
    {
        _smallpoly_alloc[k].resize(_smallpoly_power + 1);
        for (j = 0; j <= _tinypoly_power; j++)
        {
            _smallpoly_alloc[k][j].resize((size_t)1 << (_smallpoly_power - _tinypoly_power), Poly(_poly_mult[j]));
            for (i = 0; i < _smallpoly_alloc[k][j].size(); i++)
                _poly_mult[j].free(_smallpoly_alloc[k][j][i]);
        }
        for (j = _tinypoly_power + 1; j <= _smallpoly_power; j++)
            _smallpoly_alloc[k][j].clear();
    }
    smallpoly_init_level(data, data_level);
}

template<class Element>
void Stage2Poly<Element>::smallpoly_init_level(gwarray data, int data_level)
{
    int i, j, k;
    j = data_level;
    for (k = 0; k < _smallpoly_alloc.size(); k++)
    {
        if (j < _tinypoly_power)
        {
            for (i = 0; i < _smallpoly_alloc[k][j].size(); i++)
                _poly_mult[j].init(data + ((_smallpoly_count - 1 - k) << _smallpoly_power) + (i << _tinypoly_power), (size_t)1 << _tinypoly_power, false, true, _smallpoly_alloc[k][j][i]);
        }
        else
        {
            _smallpoly_alloc[k][j].resize((size_t)1 << (_smallpoly_power - j), Poly(_poly_mult[j]));
            for (i = 0; i < _smallpoly_alloc[k][j].size(); i++)
                _poly_mult[j].init(data + ((_smallpoly_count - 1 - k) << _smallpoly_power) + (i << j), (size_t)1 << j, false, true, _smallpoly_alloc[k][j][i]);
        }
    }
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

        if ((convert || preserve) && next.empty())
        {
            next.resize((cur.size() + 1)/2, Poly(pm_next));
            if (gwarrays[j] == nullptr)
                gwarrays[j] = gwalloc_array(pm_next.gw().gwdata(), level_size);
            for (i = 0; i < next.size(); i++)
                pm_next.init(gwarrays[j] + (i << j), level_size >= ((i + 1) << j) ? (size_t)1 << j : level_size - (i << j), false, false, next[i]);
        }
        else
            next.resize((cur.size() + 1)/2, Poly(pm_next));
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
        }
        if (!preserve)
            cur.clear();
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
                polymult(pm.pmdata(), cur.data() + offset, degree, cur.data() + offset + degree, sb, res.data() + offset, degree + sb, options | POLYMULT_INVEC1_MONIC | POLYMULT_INVEC2_MONIC);
        }
        if (convert)
            pm.convert(res, pm_next, next);
        if (!preserve)
            pm.free(cur);
    }
}

void rem_tree(std::vector<Poly>& rem, int k, std::vector<std::vector<Poly>>& tree, int tree_size, std::vector<arithmetic::PolyMult>& poly_mult, int cur_level, int base_level, int stop_level, std::vector<gwarray>& gwarrays_down, int level_size)
{
    int i, j;
    std::vector<Poly> rem_prev;
    rem_prev.reserve(tree_size);
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
                //res1 = std::move(rem[i]);
                //res1 >>= (1 << (base_level + j));
                rem[i] >>= (1 << (base_level + j));
                res1 = rem[i];
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
    gwarray last_gwarray = gwarrays_down[j];
    poly_mult[tiny_power - 1].init(last_gwarray + (k << tiny_power), (size_t)1 << tiny_power, false, rem.monic(), rem);
    for (j = tiny_power - 1; j >= 0; j--)
    {
        PolyMult& pm = poly_mult[j];
        PolyMult& pm_prev = poly_mult[j > 0 ? j - 1 : 0];
        bool convert = &pm.gw() != &pm_prev.gw();
        int options = !convert ? POLYMULT_NEXTFFT : 0;
        Poly rem_prev(pm_prev);

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
                polymult_several(pm.pmdata(), rem.data() + offset, size < 2*degree ? size : 2*degree, args, 2, options | POLYMULT_MULHI | (size < 2*degree ? POLYMULT_INVEC1_MONIC : 0));
            }
            if (convert)
            {
                for (i = 0; i < size && i < degree; i++)
                    gwconvert(pm.gw().gwdata(), pm_prev.gw().gwdata(), rem.data()[offset + i], rem_prev.data()[offset + i]);
                for (i = 0; i < size && i < degree; i++)
                    gwconvert(pm.gw().gwdata(), pm_prev.gw().gwdata(), rem.data()[offset + degree + i], rem_prev.data()[offset + degree + i]);
            }
        }
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
        _X1 = &Xn;
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

    GWASSERT((_poly_mod_degree & ((1 << _smallpoly_power) - 1)) == 0);
    _gwarrays.resize(poly_power() + 1);
    _gwarrays_tmp.resize(poly_power() + 1);
    _poly_mod.resize(poly_power() + 1);

    _gwarrays[0] = gwalloc_array(_poly_mult[0].gw().gwdata(), _poly_mod_degree);
    smallpoly_init(_poly_mod_degree >> _smallpoly_power, _gwarrays[0], 0);
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
    _workdone_signal.wait(lock, [&] {return _thread_exception || _smallpoly.size() == _smallpoly_count; });
    lock.unlock();
    if (_thread_exception)
        std::rethrow_exception(_thread_exception);
    GWASSERT(_workqueue.empty());

    _smallpoly_mod = std::move(_smallpoly);
    for (k = 0; k < _smallpoly_mod.size(); k++)
        _poly_mod[_smallpoly_power].emplace_back(std::move(_smallpoly_mod[k][_smallpoly_power][0]));

    build_tree(_poly_mod, _poly_mult, _smallpoly_power, false, [&](int j) { return (j - _smallpoly_power - 1)%_poly_optmem == 0; }, _gwarrays, _poly_mod_degree, _gwarrays_tmp);
    GWASSERT(_poly_mod[poly_power()][0].degree() == _poly_mod_degree);

#ifdef DEBUG_STAGE2POLY_MOD
    std::vector<std::unique_ptr<Element>> elements;
    for (i = 1, j = 0; i < _D/2; i++)
        if (gcd(_D, i) == 1)
            arithmetic().mul(X1, i, *elements.emplace_back(new Element(arithmetic())));
    GWASSERT(elements.size() == _poly_mod_degree);
    Poly poly(_poly_mult[0], _poly_mod_degree, false);
    _workers[0]->elements_to_gwnums(elements, poly.data());
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
            GWASSERT(_poly_mod[poly_power()][0].eval((GWNum&)poly.at(i)) == 0);
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
    _workers[0]->elements_to_gwnums(elements, poly_roots.data());
    elements.clear();

    _gwarrays_tmp[0] = gwalloc_array(gw().gwdata(), _poly_mod_degree);
    for (k = 0; k < _smallpoly_count; k++)
    {
        Poly poly_rem(_poly_mult[0]);
        _poly_mult[0].init(_gwarrays_tmp[0] + (k << _smallpoly_power), (size_t)1 << _smallpoly_power, false, false, poly_rem);
        for (i = 0; i < (1 << _smallpoly_power); i++)
        {
            (GWNum&)poly_rem.at(i) = 1;
            for (j = 0; j < poly_roots.size(); j++)
                gw().submul((GWNum&)poly_roots.at(j), (GWNum&)_smallpoly_mod[k][0][i >> _tinypoly_power].at(i & ((1 << _tinypoly_power) - 1)), (GWNum&)poly_rem.at(i), (GWNum&)poly_rem.at(i), j < poly_roots.size() - 1 ? GWMUL_STARTNEXTFFT : 0);
        }
        _smallpoly_mod[k][_smallpoly_power][0] = std::move(poly_rem);
    }
#endif

    _logging->info("calculating reciprocal...\n");
    PolyMult& pm = _poly_mult[poly_power()];
    Poly& modulus = _poly_mod[poly_power()][0];
    Poly reciprocal(pm);
    _gwarrays.push_back(gwalloc_array(pm.gw().gwdata(), 1 << poly_power()));
    pm.init(_gwarrays.back(), (1 << poly_power()) - 1, false, true, reciprocal);
    pm.reciprocal(modulus, reciprocal, POLYMULT_NEXTFFT);
    int degree = (1 << (poly_power() + 1)) - modulus.degree();
    if (reciprocal.degree() < degree)
        reciprocal <<= degree - reciprocal.degree();

    _logging->info("preprocessing...\n");
    _modulus.reset(new Poly(pm));
    _reciprocal.reset(new Poly(pm));
    pm.preprocess(modulus, *_modulus, 1 << (poly_power() + 1));
    pm.preprocess(reciprocal, *_reciprocal, 1 << (poly_power() + 1));
#ifdef DEBUG_STAGE2POLY_MOD
    Poly polyT(pm);
    pm.mul_range(*_reciprocal, *_modulus, polyT, 1 << poly_power(), 1 << poly_power(), 0);
    for (i = (1 << poly_power()) - 1; i >= 1; i--)
        GWASSERT(polyT.at(i) == 0);
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
        std::unique_lock<std::mutex> lock(_work_mutex);
        _smallpoly_alloc = std::move(_smallpoly);
        if (!_accumulator)
            degree = count >= _poly_mod_degree ? _poly_mod_degree : count;
        else
            degree = count >= (1 << poly_power()) ? (1 << poly_power()) : count;
        for (j = 0; j <= _smallpoly_power && gwarrays_up[j] == nullptr; j++);
        smallpoly_init(((degree + (1 << _smallpoly_power) - 1) >> _smallpoly_power), gwarrays_up[j], j);

        _workstage = 2;
        lock.unlock();
        _workstage_signal.notify_all();
        _workers[0]->run();
        lock.lock();
        _workdone_signal.wait(lock, [&] {return _thread_exception || _smallpoly.size() == _smallpoly_count; });
        lock.unlock();
        if (_thread_exception)
            std::rethrow_exception(_thread_exception);

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
        for (i = 0; i < (1 << _smallpoly_power); i++)
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

    _smallpoly_count = _smallpoly_mod.size();
    std::vector<Poly> poly_rem;
    poly_rem.reserve(_smallpoly_count);
    poly_rem.emplace_back(pm);
    std::vector<Poly> poly_rem_base;
    poly_rem_base.reserve(_smallpoly_count);

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

    int cur_level = poly_power() - 1;
    while (cur_level > _smallpoly_power - 1)
    {
        int base_level = cur_level;
        while ((base_level - _smallpoly_power)%_poly_optmem != 0)
            base_level--;
        GWASSERT(!_poly_mod[base_level][0].empty());

        std::vector<std::vector<Poly>> tree(cur_level - base_level + 1);
        int base_size = (1 << tree.size());
        int level_size = (1 << (cur_level + 1)) <= _poly_mod_degree ? (1 << (cur_level + 1)) : _poly_mod_degree;
        for (k = 0; k < poly_rem.size(); k++)
        {
            tree[0].resize((k + 1)*base_size <= _poly_mod[base_level].size() ? base_size : _poly_mod[base_level].size() - k*base_size, Poly(_poly_mult[base_level]));
            for (i = 0; i < tree[0].size(); i++)
                tree[0][i] = std::move(_poly_mod[base_level][k*base_size + i]);

            build_tree(tree, _poly_mult, base_level, true, [](int j) { return true; }, _gwarrays, level_size, _gwarrays_tmp);

            std::vector<Poly> rem;
            rem.reserve(tree[0].size());
            rem.emplace_back(std::move(poly_rem[k]));
            
            rem_tree(rem, k, tree, tree[0].size(), _poly_mult, cur_level, base_level, 0, gwarrays_down, 1 << poly_power());

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
        _smallpoly_rem = std::move(poly_rem);
        GWASSERT(_smallpoly_rem.size() == _smallpoly_mod.size());

        _workstage = 3;
        lock.unlock();
        _workstage_signal.notify_all();
        _workers[0]->run();
        lock.lock();
        _workdone_signal.wait(lock, [&] {return _thread_exception || _smallpoly_count == 0; });
        lock.unlock();
        if (_thread_exception)
            std::rethrow_exception(_thread_exception);
    }

    if (_gwarrays_tmp[0] != nullptr)
        gwfree_array(gw().gwdata(), _gwarrays_tmp[0]);
    _gwarrays_tmp.clear();
    for (j = 0; j < _smallpoly_power; j++)
        if (_gwarrays[j] != nullptr)
            gwfree_array(_poly_mult[j].gw().gwdata(), _gwarrays[j]);
    gwfree_array(_poly_mult[_smallpoly_power - 1].gw().gwdata(), _gwarrays[_smallpoly_power]);
    _gwarrays.clear();
    GWASSERT(gw().gwdata()->array_list == nullptr);

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
    while (!_threads.empty())
    {
        _threads.back().join();
        _threads.pop_back();
    }
    _workers.clear();
    _workqueue.clear();
    
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
    _gwstate.gwdata()->gwnum_max_free_count = 4*(1 << power) + 10;
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
    std::vector<std::unique_ptr<Element>> elements;
    elements.reserve((size_t)1 << power);
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
                    _stage2._smallpoly.push_back(std::move(poly_tree));
                    workitem_done = false;
                    if (_stage2._smallpoly.size() == _stage2._smallpoly_count && !_main)
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
                }

                _stage2._workqueue_signal.wait(lock, [&] {return !_stage2._workqueue.empty(); });
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
                _stage2._workqueue.emplace_back();
                _stage2._workqueue.back().Xn = std::move(TXn);
                _stage2._workqueue.back().Xn1 = std::move(TXn1);
                _stage2._workqueue.back().Xdn = std::move(TXdn);
                _stage2._workqueue.back().Xdn1 = std::move(TXdn1);
                _stage2._workqueue.back().count = distance >= (1 << power) ? (1 << power) : distance;
                _stage2._workqueue.back().n = n + (1 << power);
                _stage2._workqueue.back().distance = distance;
                lock.unlock();
                _stage2._workqueue_signal.notify_one();

                distance = 1;
            }
            if (workstage == 2 && _stage2.poly_check() && !check_root)
                check_root.reset(new GWNumWrapper(_stage2._smallpoly_mod[_stage2._poly_check_index >> power][0][(_stage2._poly_check_index & ((1 << power) - 1)) >> tiny_power].at(_stage2._poly_check_index & ((1 << tiny_power) - 1))));
            
            while (count > 0)
            {
                while (elements.size() < ((size_t)1 << power) && count > 0)
                {
                    std::unique_ptr<Element> Xtmp;
                    if (count > 2)
                    {
                        Xtmp.reset(new Element(arithmetic()));
                        if (n > 0)
                            arithmetic().add(*_stage2._X1, *Xn1, *Xn, *Xtmp);
                        else
                            arithmetic().dbl(*Xn1, *Xtmp);
                    }
                    std::swap(Xn, Xn1);
                    std::swap(Xn1, Xtmp);
                    if (workstage > 1 || gcd(_stage2.D(), n) == 1)
                        elements.push_back(std::move(Xtmp));
                    count--;
                    n += distance;
                }
                
                if (elements.size() == ((size_t)1 << power) || workstage == 2)
                {
                    poly_tree[0].resize((size_t)1 << (power - tiny_power), Poly(_poly_mult[0]));
                    if (poly_tree[0][0].empty())
                    {
                        if (gwarrays[0] == nullptr)
                            gwarrays[0] = gwalloc_array(_poly_mult[0].gw().gwdata(), 1 << power);
                        for (k = 0; k < poly_tree[0].size(); k++)
                            _poly_mult[0].init(gwarrays[0] + (k << tiny_power), (size_t)1 << tiny_power, false, true, poly_tree[0][k]);
                    }
                    if (elements.size() < ((size_t)1 << power))
                    {
                        for (j = 0; j <= tiny_power; j++)
                        {
                            poly_tree[j].erase(poly_tree[j].begin() + ((elements.size() + ((size_t)1 << tiny_power) - 1) >> tiny_power), poly_tree[j].end());
                            if ((elements.size() >> tiny_power) < poly_tree[j].size() && !poly_tree[j][elements.size() >> tiny_power].empty())
                                _poly_mult[j].alloc(poly_tree[j][elements.size() >> tiny_power], elements.size() & ((1 << tiny_power) - 1));
                        }
                    }
                    std::vector<gwnum> roots(elements.size());
                    for (i = 0; i < roots.size(); i++)
                        roots[i] = poly_tree[0][i >> tiny_power].data()[i & ((1 << tiny_power) - 1)];
                    elements_to_gwnums(elements, roots.data());
                    elements.clear();

                    if (tiny_power > 0)
                        for (k = 0; k < poly_tree[0].size(); k++)
                            build_tree_tinypoly(poly_tree, k, tiny_power, _poly_mult, 0, [&](int j) { return workstage == 1 && (j - 1)%_stage2._poly_optmem_small == 0; }, gwarrays, 1 << power, gwarrays_tmp, check_root ? check_root.get() : nullptr, check_root ? _check.get() : nullptr);
                    if (tiny_power < power)
                        build_tree(poly_tree, _poly_mult, tiny_power, false, [&](int j) { return workstage == 1 && (j - 1)%_stage2._poly_optmem_small == 0; }, gwarrays, 1 << power, gwarrays_tmp);
                    workitem_done = true;
                    
                    if (count > 0)
                    {
                        std::unique_lock<std::mutex> lock(_stage2._work_mutex);
                        _stage2._smallpoly.push_back(std::move(poly_tree));
                        workitem_done = false;
                        if (!_stage2._smallpoly_alloc.empty())
                        {
                            poly_tree = std::move(_stage2._smallpoly_alloc.back());
                            _stage2._smallpoly_alloc.pop_back();
                        }
                        else
                        {
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
                    }
                }
            }
            GWASSERT(!Xn);
            GWASSERT(!Xn1);
        }

        check_root.reset();
        workitem_done = false;
        while (workstage == 3)
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
                if ((workstage = _stage2._workstage) > 3)
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
            if (tiny_power < power)
                poly_tree[tiny_power].resize((size_t)1 << (power - tiny_power), Poly(_poly_mult[tiny_power]));
            if (tiny_power > 0)
                for (k = 0; k < (1 << (power - tiny_power)); k++)
                    build_tree_tinypoly(poly_tree, k, tiny_power, _poly_mult, 0, [](int j) { return true; }, gwarrays, 1 << power, gwarrays_tmp, nullptr, nullptr);
            if (tiny_power < power - 1)
                build_tree(poly_tree, _poly_mult, tiny_power, false, [](int j) { return true; }, gwarrays, 1 << power, gwarrays_tmp);

            rem_tree(poly_rem, 0, poly_tree, 1 << (power - tiny_power), _poly_mult, power - 1, 0, tiny_power, gwarrays_down, 1 << power);

            if (tiny_power > 0)
            {
                GWASSERT(poly_rem.size() == poly_tree[tiny_power - 1].size());
                for (k = 0; k < poly_rem.size(); k++)
                    rem_tree_tinypoly(poly_rem, k, poly_tree, tiny_power, _poly_mult, gwarrays_down, 1 << power);
            }

#ifdef DEBUG_STAGE2POLY_REM
            GWASSERT(debug_poly_rem.size() == (poly_rem.size() << tiny_power));
            for (i = 0; i < debug_poly_rem.size(); i++)
                GWASSERT(gw().cmp((GWNum&)debug_poly_rem.at(i), (GWNum&)poly_rem[i >> tiny_power].at(i & ((1 << tiny_power) - 1))) == 0);
#endif
            if (check_root)
            {
                if (gw().cmp(*check_root, (GWNum&)poly_rem[(_stage2._poly_check_index & ((1 << power) - 1)) >> tiny_power].at(_stage2._poly_check_index & ((1 << tiny_power) - 1))) != 0)
                {
                    _stage2._logging->error("check failed.\n");
#ifdef DEBUG_STAGE2POLY_REM
                    GWASSERT(0);
#endif
                    throw TaskAbortException();
                }
                _stage2._logging->error("check passed.\n");
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
        _stage2._thread_exception = std::current_exception();
        _stage2._workdone_signal.notify_one();
    }
}


void PP1Stage2Poly::init(InputNum* input, GWState* gwstate, File* file, Logging* logging, Giant& P, bool minus1)
{
    logging->set_prefix(input->display_text() + (minus1 ? ", P-1 stage 2, " : ", P+1 stage 2, "));
    Stage2Poly::init(input, gwstate, file, /*read_state<State>(file)*/nullptr, logging);
    _lucas.reset(new LucasVArithmetic());
    _P = P;
}

template void Stage2Poly<LucasV>::done(const arithmetic::Giant& factor);

void PP1Stage2Poly::SmallPolyWorker::elements_to_gwnums(std::vector<std::unique_ptr<arithmetic::LucasV>>& elements, gwnum* data)
{
    for (int i = 0; i < elements.size(); i++)
        gw().unfft(elements[i]->V(), (GWNum&)GWNumWrapper(gw(), data[i]));
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

    _workers.reserve(_poly_threads);
    _workers.emplace_back(new SmallPolyWorker(*this));
    _workers[0]->set_main(true);
    for (int i = 1; i < _poly_threads; i++)
    {
        _workers.emplace_back(new SmallPolyWorker(*this));
        _threads.emplace_back(&SmallPolyWorker::run, _workers.back().get());
    }

    if (!_modulus)
    {
        _W->V() = _P;
        gw().fft(_W->V(), _W->V());
        _X1 = _W.get();
        Stage2Poly::setup();
    }

    if (state() == nullptr)
        reset_state<State>();
    int v = state()->iteration() + _first_D;
    int count = iterations() - state()->iteration();
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
    Stage2Poly::init(input, gwstate, file, /*read_state<State>(file)*/nullptr, logging);
    _X = X;
    _Y = Y;
    _Z = Z;
    _T = T;
    _EdD = EdD;
}

template void Stage2Poly<EdY>::done(const arithmetic::Giant& factor);

void EdECMStage2Poly::SmallPolyWorker::elements_to_gwnums(std::vector<std::unique_ptr<arithmetic::EdY>>& elements, gwnum* data)
{
    using Iter = typename std::remove_reference<decltype(elements)>::type::iterator;
    Iter begin = elements.begin();
    Iter end = elements.end();
    Iter it;
    Iter first = end;
    Iter last = end;
    for (it = begin; it != end; it++)
    {
        if (!(*it))
            continue;
        if (first == end)
            first = it;
        last = it;
        if (!(*it)->Z)
        {
            (*it)->Z.reset(new GWNum(gw()));
            *(*it)->Z = 1;
        }
        if (!(*it)->ZpY)
            (*it)->ZpY.reset(new GWNum(gw()));
        (*it)->ZmY.reset();
    }
    if (first == end)
        return;
    std::swap((*first)->ZpY, (*first)->Z);
    Iter prev = first;
    for ((it = first)++; it != end; it++)
        if (*it)
        {
            gw().mul(*(*prev)->ZpY, *(*it)->Z, *(*it)->ZpY, it != last ? GWMUL_STARTNEXTFFT : 0);
            prev = it;
        }
    try
    {
        //GWNum tmp(gw());
        //tmp = *(*last)->ZpY;
        gw().inv(*(*last)->ZpY, *(*last)->ZpY);
        //gw().mul(tmp, *(*last)->ZpY, tmp);
        //GWASSERT(tmp == 1);
    }
    catch (const ArithmeticException&)
    {
        std::swap((*first)->ZpY, (*first)->Z);
        throw;
    }
    size_t i = 0;
    for ((it = last)++, prev = last; it != first;)
    {
        it = prev;
        if (it != first)
        {
            for ((prev = it)--; !(*prev); prev--);
            gw().mul(*(*it)->ZpY, *(*prev)->ZpY, *(*prev)->ZpY, GWMUL_STARTNEXTFFT);
            std::swap((*it)->ZpY, (*prev)->ZpY);
            gw().mul(*(*it)->Z, *(*prev)->ZpY, *(*prev)->ZpY, GWMUL_STARTNEXTFFT);
        }
        if ((*it)->Y)
            gw().mul(*(*it)->ZpY, *(*it)->Y, (GWNum&)GWNumWrapper(gw(), data[i++]), 0);
        (*it)->Z.reset();
        (*it)->ZpY.reset();
    }
}

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

    _workers.reserve(_poly_threads);
    _workers.emplace_back(new SmallPolyWorker(*this));
    _workers[0]->set_main(true);
    for (int i = 1; i < _poly_threads; i++)
    {
        _workers.emplace_back(new SmallPolyWorker(*this));
        _threads.emplace_back(&SmallPolyWorker::run, _workers.back().get());
    }

    try
    {
        if (!_modulus)
        {
            _W->deserialize(_Y, _Z);
            _X1 = _W.get();
            Stage2Poly::setup();
        }

        if (state() == nullptr)
            reset_state<State>();
        v = state()->iteration() + _first_D;
        count = iterations() - state()->iteration();
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
    _workqueue.back().Xn.reset(new EdY(arithmetic(), EdP));
    arithmetic().optimize(*_workqueue.back().Xn);
    if (_workqueue.back().count > 1)
    {
        _workqueue.back().Xn1.reset(new EdY(arithmetic(), EdP1));
        arithmetic().optimize(*_workqueue.back().Xn1);
    }

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
        _workqueue.back().Xdn.reset(new EdY(arithmetic(), EdP1));
        arithmetic().optimize(*_workqueue.back().Xdn);

        if (_workqueue.back().distance > (1 << _smallpoly_power) + 1)
        {
            ed.add(EdP1, EdW, EdP1, ed.ED_PROJECTIVE);
            _workqueue.back().Xdn1.reset(new EdY(arithmetic(), EdP1));
            arithmetic().optimize(*_workqueue.back().Xdn1);
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
        _workqueue.back().Xn.reset(new EdY(arithmetic(), EdP));
        arithmetic().optimize(*_workqueue.back().Xn);

        if (_workqueue.back().count > 1)
        {
            ed.add(EdP, EdW, EdP1, ed.ED_PROJECTIVE);
            _workqueue.back().Xn1.reset(new EdY(arithmetic(), EdP1));
            arithmetic().optimize(*_workqueue.back().Xn1);
        }

        if (_workqueue.back().distance > (1 << _smallpoly_power))
        {
            if (_workqueue.back().n == (1 << _smallpoly_power))
                ed.dbl(EdWd, EdP1);
            else
                ed.add(EdP, EdWd, EdP1);
            _workqueue.back().Xdn.reset(new EdY(arithmetic(), EdP1));
            arithmetic().optimize(*_workqueue.back().Xdn);
            
            if (_workqueue.back().distance > (1 << _smallpoly_power) + 1)
            {
                ed.add(EdP1, EdW, EdP1, ed.ED_PROJECTIVE);
                _workqueue.back().Xdn1.reset(new EdY(arithmetic(), EdP1));
                arithmetic().optimize(*_workqueue.back().Xdn1);
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
