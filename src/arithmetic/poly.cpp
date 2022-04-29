#define GDEBUG
#include <algorithm>
#include <stdlib.h>
#include <string.h>
#include "gwnum.h"
#include "poly.h"
#include "exception.h"

namespace arithmetic
{
    GWNum Poly::eval(GWNum& x)
    {
        GWNum res(pm().gw());
        if (size() == 0 && monic())
            res = 1;
        else if (size() == 0)
            res = 0;
        else
            res = at(0);
        std::unique_ptr<GWNum> t;
        std::unique_ptr<GWNum> p;
        int i, d;
        d = degree();
        for (i = 1; i <= d; i++)
        {
            if (i == 2)
            {
                p.reset(new GWNum(x));
                pm().gw().mul(x, x, *p, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT);
            }
            if (i > 2)
                pm().gw().mul(x, *p, *p, GWMUL_STARTNEXTFFT);
            if (i < size() && _poly[i] != nullptr)
            {
                GWNumWrapper a = at(i);
                if (!t)
                    t.reset(new GWNum(pm().gw()));
                if (i == 1)
                    pm().gw().mul(x, a, *t, GWMUL_FFT_S1 | GWMUL_PRESERVE_S2);
                else
                    pm().gw().mul(*p, a, *t, GWMUL_PRESERVE_S2);
                pm().gw().add(res, *t, res, i < d ? GWADD_DELAY_NORMALIZE : GWADD_FORCE_NORMALIZE);
            }
        }
        if (monic() && d > 1)
            pm().gw().add(res, *p, res, GWADD_FORCE_NORMALIZE);
        if (monic() && d == 1)
            pm().gw().add(res, x, res, GWADD_FORCE_NORMALIZE);
        return res;
    }

    Poly Poly::reciprocal(int precision, int options)
    {
        Poly res(pm(), precision, monic());
        pm().reciprocal(*this, res, options);
        return res;
    }

    PolyMult::PolyMult(GWArithmetic& gw, int max_threads) : _gw(gw)
    {
        _max_output = gw.state().max_polymult_output();
        polymult_init(pmdata(), gw.gwdata());
#ifdef polymult_set_max_num_threads
        polymult_set_max_num_threads(pmdata(), max_threads);
#endif
    }

    PolyMult::~PolyMult()
    {
        polymult_done(pmdata());
    }

    void PolyMult::set_threads(int threads)
    {
#ifdef polymult_set_max_num_threads
        polymult_set_num_threads(pmdata(), threads);
#endif
    }

    void PolyMult::alloc(Poly& a)
    {
        for (auto it = a._poly.begin(); it != a._poly.end(); it++)
            if (*it == nullptr)
                *it = gwalloc(gw().gwdata());
        if (a._cache != nullptr)
            gwfree_array(gw().gwdata(), a._cache);
        a._cache = nullptr;
    }

    void PolyMult::alloc(Poly& a, int size)
    {
        for (size_t i = size; i < a.size(); i++)
            gwfree(gw().gwdata(), a._poly[i]);
        a._poly.resize(size);
        for (auto it = a._poly.begin(); it != a._poly.end(); it++)
            if (*it == nullptr)
                *it = gwalloc(gw().gwdata());
        if (a._cache != nullptr)
            gwfree_array(gw().gwdata(), a._cache);
        a._cache = nullptr;
    }

    void PolyMult::free(Poly& a)
    {
        for (auto it = a._poly.begin(); it != a._poly.end(); it++)
            gwfree(gw().gwdata(), *it);
        a._poly.clear();
        if (a._cache != nullptr)
            gwfree_array(gw().gwdata(), a._cache);
        a._cache = nullptr;
        a._monic = false;
    }

    void PolyMult::copy(const Poly& a, Poly& res)
    {
        alloc(res, a.size());
        for (size_t i = 0; i < a.size(); i++)
            gwcopy(gw().gwdata(), a._poly[i], res._poly[i]);
        res._monic = a._monic;
    }

    void PolyMult::move(Poly&& a, Poly& res)
    {
        if (!res.empty())
            free(res);
        res._poly = std::move(a._poly);
        res._monic = a._monic;
        res._cache = a._cache;
        res._cache_size = a._cache_size;
        a._cache = nullptr;
    }

    void PolyMult::fft(const Poly& a, Poly& res)
    {
        alloc(res, a.size());
        for (size_t i = 0; i < a.size(); i++)
            gwfft(gw().gwdata(), a._poly[i], res._poly[i]);
        res._monic = a._monic;
    }

    void PolyMult::init(bool monic, Poly& res)
    {
        free(res);
        res._monic = monic;
    }

    void PolyMult::init(const GWNum& a, bool monic, Poly& res)
    {
        alloc(res, 1);
        gwcopy(gw().gwdata(), *a, res._poly[0]);
        res._monic = monic;
    }

    void PolyMult::init(GWNum&& a, bool monic, Poly& res)
    {
        if (!res.empty())
            free(res);
        res._poly.push_back(a._gwnum);
        a._gwnum = nullptr;
        res._monic = monic;
    }

    void PolyMult::mul(Poly& a, Poly& b, Poly& res, int options)
    {
        if (a.degree() < b.degree() || a.size() < b.size())
        {
            mul(b, a, res, options);
            return;
        }
        if (b.degree() < 0)
        {
            free(res);
            return;
        }
        if (b.size() == 0) // monic
        {
            if (&a != &res)
                copy(a, res);
            return;
        }

        int sa = a.size();
        int sb = b.size();
        alloc(res, sa + sb - (a.monic() || b.monic() ? 0 : 1));

        if (sb == 1)
        {
            if (a.monic() && b.monic())
                gwadd3o(gw().gwdata(), a._poly[sa - 1], b._poly[0], res._poly[sa], GWADD_FORCE_NORMALIZE);
            else if (a.monic())
                gwcopy(gw().gwdata(), b._poly[0], res._poly[sa]);
            else if (b.monic())
                gwcopy(gw().gwdata(), a._poly[sa - 1], res._poly[sa]);
            for (int i = sa - 1; i > 0; i--)
                if (b.monic())
                    gwmuladd4(gw().gwdata(), a._poly[i], b._poly[0], a._poly[i - 1], res._poly[i], (options & POLYMULT_STARTNEXTFFT ? GWMUL_STARTNEXTFFT : 0) | GWMUL_FFT_S1 | GWMUL_FFT_S2);
                else
                    gwmul3(gw().gwdata(), a._poly[i], b._poly[0], res._poly[i], (options & POLYMULT_STARTNEXTFFT ? GWMUL_STARTNEXTFFT : 0) | GWMUL_FFT_S1 | GWMUL_FFT_S2);
            gwmul3(gw().gwdata(), a._poly[0], b._poly[0], res._poly[0], (options & POLYMULT_STARTNEXTFFT ? GWMUL_STARTNEXTFFT : 0) | GWMUL_FFT_S1 | GWMUL_FFT_S2);
            res._monic = a.monic() && b.monic();
            return;
        }

        polymult(pmdata(), a.data(), sa, b.data(), sb, res.data(), res.size(), options | (a.monic() ? POLYMULT_INVEC1_MONIC : 0) | (b.monic() ? POLYMULT_INVEC2_MONIC : 0));
        res._monic = a.monic() && b.monic();
    }

    void PolyMult::mul(Poly&& a, Poly&& b, Poly& res, int options)
    {
        GWASSERT(&a.pm().gw() == &res.pm().gw());
        GWASSERT(&b.pm().gw() == &res.pm().gw());
        if (a.degree() < b.degree() || a.size() < b.size())
        {
            mul(std::move(b), std::move(a), res, options);
            return;
        }
        if (b.degree() < 0)
        {
            free(res);
            return;
        }
        if (&b == &res)
        {
            Poly tmp(*this, 0, false);
            move(std::move(b), tmp);
            mul(std::move(a), std::move(tmp), res, options);
            return;
        }

        int sa = a.size();
        int sb = b.size();
        if (&a != &res)
            move(std::move(a), res);
        else
            alloc(res, sa); // free cache

        if (sb == 0) // monic
            return;
        if (sb == 1)
        {
            mul(res, b, res, options);
            free(b);
            return;
        }

        res._poly.resize(sa + sb - (res.monic() || b.monic() ? 0 : 1));
        for (int i = 0; sa + i < res.size(); i++)
            res._poly[sa + i] = b._poly[i];
        polymult(pmdata(), res._poly.data(), sa, b._poly.data(), sb, res._poly.data(), res.size(), options | (res.monic() ? POLYMULT_INVEC1_MONIC : 0) | (b.monic() ? POLYMULT_INVEC2_MONIC : 0));
        if (!res.monic() && !b.monic())
            gwfree(gw().gwdata(), b._poly[sb - 1]);
        b._poly.clear();
        res._monic = res.monic() && b.monic();
    }

    void PolyMult::mul_lohi(Poly&& a, Poly&& b, Poly& res_lo, Poly& res_hi, int half, int options)
    {
        int sa = a.size();
        int sb = b.size();
        mul(std::move(a), std::move(b), res_lo, options);
        free(res_hi);

        int sr = sa + sb - (a.monic() || b.monic() ? 0 : 1);
        if (sr >= half)
        {
            res_hi._poly.insert(res_hi._poly.begin(), res_lo._poly.begin() + half, res_lo._poly.end());
            res_lo._poly.erase(res_lo._poly.begin() + half, res_lo._poly.end());
            res_hi._monic = res_lo._monic;
            res_lo._monic = false;
        }
    }

    void PolyMult::preprocess(Poly& res, int size)
    {
        if (res.size() > 1 && res._cache == nullptr)
        {
            res._cache = polymult_preprocess(pmdata(), res._poly.data(), res._poly.size(), size, size, POLYMULT_CIRCULAR | POLYMULT_FFT | (res.monic() ? POLYMULT_INVEC1_MONIC : 0));
            res._cache_size = res._poly.size();
            for (auto it = res._poly.begin(); it != res._poly.end(); it++)
                gwfree(gw().gwdata(), *it);
            res._poly.clear();
        }
    }

    void PolyMult::preprocess_and_mul(Poly& a, Poly& b, Poly& res, int size, int options)
    {
        GWASSERT(&a.pm().gw() == &res.pm().gw());
        GWASSERT(&b.pm().gw() == &res.pm().gw());
        GWASSERT(&a != &res);
        GWASSERT(&b != &res);
        if (a.degree() < b.degree() || a.size() < b.size())
        {
            preprocess_and_mul(b, a, res, size, options);
            return;
        }
        if (b.degree() < 0)
        {
            free(res);
            return;
        }

        int sa = a.size();
        int sb = b.size();
        if (sa > 1)
        {
            move(std::move(a), res);
            a._cache = polymult_preprocess(pmdata(), res._poly.data(), sa, size, size, POLYMULT_CIRCULAR | POLYMULT_FFT | (res.monic() ? POLYMULT_INVEC1_MONIC : 0));
            a._cache_size = sa;
        }
        else
            copy(a, res);

        if (sb == 0) // monic
            return;
        if (sb == 1)
        {
            mul(res, b, res, options);
            return;
        }
        b._cache = polymult_preprocess(pmdata(), b._poly.data(), sb, size, size, POLYMULT_CIRCULAR | POLYMULT_FFT| (b.monic() ? POLYMULT_INVEC1_MONIC : 0));
        b._cache_size = sb;

        res._poly.resize(sa + sb - (res.monic() || b.monic() ? 0 : 1));
        for (int i = 0; sa + i < res.size(); i++)
            res._poly[sa + i] = b._poly[i];
        int padding = size - res._poly.size();
        res._poly.insert(res._poly.end(), padding, nullptr);
        polymult(pmdata(), a._cache, sa, b._cache, sb, res._poly.data(), size, options | POLYMULT_CIRCULAR | (res.monic() ? POLYMULT_INVEC1_MONIC : 0) | (b.monic() ? POLYMULT_INVEC2_MONIC : 0));
        res._poly.erase(res._poly.end() - padding, res._poly.end());
        if (!res.monic() && !b.monic())
            gwfree(gw().gwdata(), b._poly[sb - 1]);
        b._poly.clear();
        res._monic = res.monic() && b.monic();
        if (res._monic && padding == 0)
        {
            gwunfft(gw().gwdata(), res._poly[0], res._poly[0]);
            gwaddsmall(gw().gwdata(), res._poly[0], -1);
        }
        if (res._cache != nullptr)
            gwfree_array(gw().gwdata(), res._cache);
        res._cache = nullptr;
    }

    void PolyMult::mul_hi(Poly& a, Poly& b, Poly& res, int options)
    {
        if (a.degree() < b.degree() || a.size() < b.size())
        {
            mul_hi(b, a, res, options);
            return;
        }
        if (b.degree() < 0 || a.degree() + b.degree() < res.size())
        {
            free(res);
            return;
        }
        if (b.size() == 0) // monic
        {
            for (int i = 0; i < res.size() && i + res.size() < a.size(); i++)
                gwcopy(gw().gwdata(), a._poly[i + res.size()], res._poly[i]);
            alloc(res, a.size() - res.size());
            res._monic = a.monic();
            return;
        }

        int sa = a.size();
        int sb = b.size();
        int half = res.size();
        int full = sa + sb - (a.monic() || b.monic() ? 0 : 1);
        if (full < 2*half)
            alloc(res, full - half);

        if (sb == 1)
        {
            if (sa >= half && sa < 2*half)
            {
                if (a.monic() && b.monic())
                    gwadd3o(gw().gwdata(), a._poly[sa - 1], b._poly[0], res._poly[sa - half], GWADD_FORCE_NORMALIZE);
                else if (a.monic())
                    gwcopy(gw().gwdata(), b._poly[0], res._poly[sa - half]);
                else if (b.monic())
                    gwcopy(gw().gwdata(), a._poly[sa - 1], res._poly[sa - half]);
            }
            for (int i = sa - 1; i >= half; i--)
                if (b.monic())
                    gwmuladd4(gw().gwdata(), a._poly[i], b._poly[0], a._poly[i - 1], res._poly[i - half], (options & POLYMULT_STARTNEXTFFT ? GWMUL_STARTNEXTFFT : 0) | GWMUL_FFT_S1 | GWMUL_FFT_S2);
                else
                    gwmul3(gw().gwdata(), a._poly[i], b._poly[0], res._poly[i - half], (options & POLYMULT_STARTNEXTFFT ? GWMUL_STARTNEXTFFT : 0) | GWMUL_FFT_S1 | GWMUL_FFT_S2);
            res._monic = a.monic() && b.monic() && sa + sb < 2*half;
            return;
        }
        
        res._poly.insert(res._poly.begin(), half, nullptr);
        int padding = 2*half - res._poly.size();
        res._poly.insert(res._poly.end(), padding, nullptr);
        // (full > 2*half ? POLYMULT_CIRCULAR : 0)
        polymult(pmdata(), a.data(), sa, b.data(), sb, res.data(), res.size(), options | POLYMULT_CIRCULAR | (a.monic() ? POLYMULT_INVEC1_MONIC : 0) | (b.monic() ? POLYMULT_INVEC2_MONIC : 0));
        res._poly.erase(res._poly.end() - padding, res._poly.end());
        res._poly.erase(res._poly.begin(), res._poly.begin() + half);
        res._monic = a.monic() && b.monic() && sa + sb < 2*half;
    }

    void PolyMult::mul_hi(Poly&& a, Poly& b, Poly& res, int half, int options)
    {
        GWASSERT(&a.pm().gw() == &res.pm().gw());
        //GWASSERT(b.preprocessed() || b.size() <= 1);
        if (b.degree() < 0 || a.degree() + b.degree() < half)
        {
            free(res);
            return;
        }
        if (b.size() == 0) // monic
        {
            if (&a != &res)
                move(std::move(a), res);
            res >>= half;
            return;
        }

        int sa = a.size();
        int sb = b.size();
        int full = sa + sb - (a.monic() || b.monic() ? 0 : 1);
        int size = full < 2*half ? full - half : half;
        Poly to_free(a.pm());
        std::vector<gwnum> poly;
        poly.reserve(size);
        for (int i = 0; i < a.size(); i++)
            if (a._poly[i] != nullptr)
            {
                if (poly.size() < size)
                    poly.push_back(a._poly[i]);
                else
                    to_free._poly.push_back(a._poly[i]);
            }
        while (poly.size() < size)
            poly.push_back(gwalloc(a.pm().gw().gwdata()));

        if (sb == 1)
        {
            if (sa >= half && sa < 2*half)
            {
                if (a.monic() && b.monic())
                    gwadd3o(gw().gwdata(), a._poly[sa - 1], b._poly[0], poly[sa - half], GWADD_FORCE_NORMALIZE);
                else if (a.monic())
                    gwcopy(gw().gwdata(), b._poly[0], poly[sa - half]);
                else if (b.monic())
                    gwcopy(gw().gwdata(), a._poly[sa - 1], poly[sa - half]);
            }
            for (int i = sa - 1; i >= half; i--)
                if (b.monic())
                    gwmuladd4(gw().gwdata(), a._poly[i], b._poly[0], a._poly[i - 1], poly[i - half], (options & POLYMULT_STARTNEXTFFT ? GWMUL_STARTNEXTFFT : 0) | GWMUL_FFT_S1 | GWMUL_FFT_S2);
                else
                    gwmul3(gw().gwdata(), a._poly[i], b._poly[0], poly[i - half], (options & POLYMULT_STARTNEXTFFT ? GWMUL_STARTNEXTFFT : 0) | GWMUL_FFT_S1 | GWMUL_FFT_S2);
            a._poly.clear();
            res._poly = std::move(poly);
            res._monic = a.monic() && b.monic() && sa + sb < 2*half;
            return;
        }

        poly.insert(poly.begin(), half, nullptr);
        int padding = 2*half - poly.size();
        if (padding > 0)
            poly.insert(poly.end(), padding, nullptr);
        polymult(pmdata(), a.data(), sa, b.data(), sb, poly.data(), 2*half, options | POLYMULT_CIRCULAR | (a.monic() ? POLYMULT_INVEC1_MONIC : 0) | (b.monic() ? POLYMULT_INVEC2_MONIC : 0));
        if (padding > 0)
            poly.erase(poly.end() - padding, poly.end());
        poly.erase(poly.begin(), poly.begin() + half);
        a._poly.clear();
        res._poly = std::move(poly);
        res._monic = a.monic() && b.monic() && sa + sb < 2*half;
    }

    void PolyMult::mul_hi(Poly&& a, Poly& b, Poly& c, Poly& res1, Poly& res2, int half, int options)
    {
        GWASSERT(&a.pm().gw() == &res1.pm().gw());
        GWASSERT(&a.pm().gw() == &res2.pm().gw());
        //GWASSERT(b.preprocessed() || b.size() <= 1);
        //GWASSERT(c.preprocessed() || c.size() <= 1);
        if (b.degree() < 0 || a.degree() + b.degree() < half)
        {
            free(res1);
            mul_hi(std::move(a), c, res2, half, options);
            return;
        }
        if (c.degree() < 0 || a.degree() + c.degree() < half)
        {
            free(res2);
            mul_hi(std::move(a), b, res1, half, options);
            return;
        }
        if (b.size() == 0 && c.size() == 0) // monic
        {
            move(std::move(a), res1);
            free(res2);
            res2._poly.insert(res2._poly.begin(), res1._poly.begin() + half, res1._poly.end());
            res1._poly.resize(half);
            res2._monic = res1._monic;
            copy(res2, res1);
            return;
        }
        if (b.size() == 0) // monic
        {
            alloc(res1, a.size() - half);
            for (int i = 0; i < res1.size(); i++)
                gwcopy(gw().gwdata(), a._poly[i + half], res1._poly[i]);
            res1._monic = a.monic() && a.size() < 2*half;
            mul_hi(std::move(a), c, res2, half, options);
            return;
        }
        if (c.size() == 0) // monic
        {
            alloc(res2, a.size() - half);
            for (int i = 0; i < res2.size(); i++)
                gwcopy(gw().gwdata(), a._poly[i + half], res2._poly[i]);
            res2._monic = a.monic() && a.size() < 2*half;
            mul_hi(std::move(a), b, res1, half, options);
            return;
        }

        int sa = a.size();
        int sb = b.size();
        int sc = c.size();
        int full1 = sa + sb - (a.monic() || b.monic() ? 0 : 1);
        int full2 = sa + sc - (a.monic() || c.monic() ? 0 : 1);
        int size1 = full1 < 2*half ? full1 - half : half;
        int size2 = full2 < 2*half ? full2 - half : half;

        if (sb == 1)
        {
            alloc(res1, half);
            mul_hi(a, b, res1, options);
            mul_hi(std::move(a), c, res2, half, options);
            return;
        }
        if (sc == 1)
        {
            alloc(res2, half);
            mul_hi(a, c, res2, options);
            mul_hi(std::move(a), b, res1, half, options);
            return;
        }

        free(res1);
        free(res2);
        res1._poly = a._poly;
        if (sa > size1)
        {
            res2._poly.insert(res2._poly.begin(), res1._poly.begin() + size1, res1._poly.end());
            res1._poly.resize(size1);
        }
        if (sa > size1 + size2)
            res2._poly.resize(size2);
        alloc(res1, size1);
        alloc(res2, size2);

#if !defined(POLYMULT_VECTOR)
        a._cache = polymult_preprocess(pmdata(), a._poly.data(), a._poly.size(), 2*half, 2*half, POLYMULT_CIRCULAR | POLYMULT_FFT | (a.monic() ? POLYMULT_INVEC1_MONIC : 0));

        res1._poly.insert(res1._poly.begin(), half, nullptr);
        int padding1 = 2*half - res1._poly.size();
        if (padding1 > 0)
            res1._poly.insert(res1._poly.end(), padding1, nullptr);
        polymult(pmdata(), a._cache, sa, b.data(), sb, res1.data(), 2*half, options | POLYMULT_CIRCULAR | (a.monic() ? POLYMULT_INVEC1_MONIC : 0) | (b.monic() ? POLYMULT_INVEC2_MONIC : 0));
        if (padding1 > 0)
            res1._poly.erase(res1._poly.end() - padding1, res1._poly.end());
        res1._poly.erase(res1._poly.begin(), res1._poly.begin() + half);

        res2._poly.insert(res2._poly.begin(), half, nullptr);
        int padding2 = 2*half - res2._poly.size();
        if (padding2 > 0)
            res2._poly.insert(res2._poly.end(), padding2, nullptr);
        polymult(pmdata(), a._cache, sa, c.data(), sc, res2.data(), 2*half, options | POLYMULT_CIRCULAR | (a.monic() ? POLYMULT_INVEC1_MONIC : 0) | (c.monic() ? POLYMULT_INVEC2_MONIC : 0));
        if (padding2 > 0)
            res2._poly.erase(res2._poly.end() - padding2, res2._poly.end());
        res2._poly.erase(res2._poly.begin(), res2._poly.begin() + half);

        gwfree_array(gw().gwdata(), a._cache);
        a._cache = nullptr;
#else
        res1._poly.insert(res1._poly.begin(), half, nullptr);
        int padding1 = 2*half - res1._poly.size();
        if (padding1 > 0)
            res1._poly.insert(res1._poly.end(), padding1, nullptr);
        res2._poly.insert(res2._poly.begin(), half, nullptr);
        int padding2 = 2*half - res2._poly.size();
        if (padding2 > 0)
            res2._poly.insert(res2._poly.end(), padding2, nullptr);

        pmarg args[2];
        args[0].invec2 = b.data();
        args[0].invec2_size = sb;
        args[0].outvec = res1.data();
        args[0].outvec_size = 2*half;
        args[0].options = b.monic() ? POLYMULT_INVEC2_MONIC : 0;
        args[1].invec2 = c.data();
        args[1].invec2_size = sc;
        args[1].outvec = res2.data();
        args[1].outvec_size = 2*half;
        args[1].options = c.monic() ? POLYMULT_INVEC2_MONIC : 0;
        polymult_vector(pmdata(), a.data(), sa, args, 2, options | POLYMULT_CIRCULAR | (a.monic() ? POLYMULT_INVEC1_MONIC : 0));

        if (padding1 > 0)
            res1._poly.erase(res1._poly.end() - padding1, res1._poly.end());
        res1._poly.erase(res1._poly.begin(), res1._poly.begin() + half);
        if (padding2 > 0)
            res2._poly.erase(res2._poly.end() - padding2, res2._poly.end());
        res2._poly.erase(res2._poly.begin(), res2._poly.begin() + half);
#endif

        for (int i = size1 + size2; i < a.size(); i++)
            gwfree(gw().gwdata(), a._poly[i]);
        a._poly.clear();
        res1._monic = a.monic() && b.monic() && sa + sb < 2*half;
        res2._monic = a.monic() && c.monic() && sa + sc < 2*half;
    }

    void PolyMult::reciprocal(Poly& a, Poly& res, int options)
    {
        GWASSERT(&a.pm().gw() == &res.pm().gw());
        int i, j;
        int d;
        for (d = 2; d < res.size() + (a.monic() ? 1 : 0); d <<= 1);

        gwnum* g = res._poly.data() + res.size() - 1;
        std::vector<gwnum> tmp(d);
        std::vector<gwnum> tmp2(d);
        std::vector<gwnum> f(d);
        for (j = d/2; j < d; j++)
            tmp[j] = gwalloc(gw().gwdata());

        int sa = a.size();
        if (a.monic())
        {
            dbltogw(gw().gwdata(), 0, *g);
            gwsub3o(gw().gwdata(), *g, a._poly[sa - 1], *g, GWADD_GUARANTEED_OK);
            for (i = 2; i < d; g -= i, i <<= 1)
            {
                for (j = 0; j < 2*i; j++)
                    f[j] = (sa - 2*i + j) >= 0 ? a._poly[sa - 2*i + j] : nullptr;
                polymult(pmdata(), g, i - 1, f.data(), 2*i, tmp.data() + d/2 - i, 2*i, POLYMULT_CIRCULAR | POLYMULT_INVEC1_MONIC);
                for (j = 0; j < i; j++)
                    tmp2[d/2 + j] = (g + j - i) >= res._poly.data() ? g[j - i] : nullptr;
                polymult(pmdata(), g, i - 1, tmp.data() + d/2 - 1, i + 1, tmp2.data() + d/2 - i, 2*i, POLYMULT_CIRCULAR | POLYMULT_INVEC1_MONIC | POLYMULT_INVEC2_NEGATE);
            }
        }
        else
        {
            Giant tmp = gw().popg();
            gwtobinary(gw().gwdata(), a._poly[sa - 1], tmp.data(), tmp.capacity());
            tmp.inv(gw().N());
        }

        for (j = d/2; j < d; j++)
            gwfree(gw().gwdata(), tmp[j]);
        res._monic = a.monic();
    }

    void PolyMult::shiftleft(Poly& a, int b, Poly& res)
    {
        GWASSERT(&a.pm().gw() == &res.pm().gw());
        if (&a != &res)
            copy(a, res);
        res._poly.insert(res._poly.begin(), b, nullptr);
    }

    void PolyMult::shiftright(Poly& a, int b, Poly& res)
    {
        GWASSERT(&a.pm().gw() == &res.pm().gw());
        if (&a == &res)
        {
            for (int i = 0; i < b && i < res.size(); i++)
                if (res._poly[i] != nullptr)
                    gwfree(gw().gwdata(), res._poly[i]);
            res._poly.erase(res._poly.begin(), res._poly.begin() + b);
        }
        else
        {
            alloc(res, (int)a.size() - b);
            for (int i = 0; i < res.size(); i++)
                gwcopy(gw().gwdata(), a._poly[b + i], res._poly[i]);
        }
    }

    void PolyMult::convert(const Poly& a, Poly& res)
    {
        res._monic = a.monic();
        res.pm().alloc(res, a.size());
        for (int i = 0; i < a.size(); i++)
            if (gwconvert(a.pm().gw().gwdata(), res.pm().gw().gwdata(), a._poly[i], res._poly[i]) != 0)
                throw InvalidFFTDataException();
    }

    void PolyMult::insert(GWNum&& a, Poly& res, size_t pos)
    {
        res._poly.insert(res._poly.begin() + pos, *a);
        a._gwnum = nullptr;
    }

    GWNum PolyMult::remove(Poly& a, size_t pos)
    {
        gwnum res = a._poly.at(pos);
        a._poly.erase(a._poly.begin() + pos);
        return GWNum(gw(), res);
    }
}
