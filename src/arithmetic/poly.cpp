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
        GWNum res(gw());
        if (front().is_small())
            res = front().small();
        else
            res = front().value();
        std::unique_ptr<GWNum> t;
        std::unique_ptr<GWNum> p;
        size_t i, l;
        l = degree();
        for (i = 1; i <= l; i++)
        {
            if (i == 2)
            {
                p.reset(new GWNum(x));
                gw().mul(x, x, *p, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT);
            }
            if (i > 2)
                gw().mul(x, *p, *p, GWMUL_STARTNEXTFFT);
            if (!at(i).is_zero())
            {
                if (!t)
                    t.reset(new GWNum(gw()));
                if (at(i).is_small())
                    *t = at(i).small();
                if (i == 1)
                    gw().mul(x, at(i).is_small() ? *t : at(i).value(), *t, GWMUL_FFT_S1 | GWMUL_PRESERVE_S2);
                else
                    gw().mul(*p, at(i).is_small() ? *t : at(i).value(), *t, GWMUL_PRESERVE_S2);
                gw().add(res, *t, res, i < l ? GWADD_DELAY_NORMALIZE : GWADD_FORCE_NORMALIZE);
            }
        }
        return res;
    }

    Poly Poly::reciprocal(int size)
    {
        PolyMulKaratsuba mul(gw());
        return mul.reciprocal(*this, size);
    }

    Poly Poly::mul(Poly& b)
    {
        PolyMul poly_mul(gw());
        return poly_mul.mul(*this, b);
    }

    Poly Poly::mul_half(Poly& b, int half)
    {
        PolyMul poly_mul(gw());
        return poly_mul.mul_half(*this, b, half);
    }

    void Poly::do_fft()
    {
        for (int i = degree(); i >= 0; i--)
            at(i).do_fft();
    }
    
    void Poly::clear_fft()
    {
        for (auto it = begin(); it != end(); it++)
            it->clear_fft();
        _fft.reset();
    }

    SubPoly::SubPoly(const Poly& poly, int offset, int count, bool transpose) : Poly(poly.gw(), 0)
    {
        for (int i = offset; (count < 0 && i < (int)poly.size()) || (count > 0 && i < offset + count); i++)
            if (i >= 0 && i < (int)poly.size())
                emplace(transpose ? begin() : end(), poly[i]);
            else
                emplace(transpose ? begin() : end());
    }

    void SubPolyFFT::init(bool force_fft)
    {
        clear();
        for (int i = _offset; (_count < 0 && i < (int)_poly.size()) || (_count > 0 && i < _offset + _count); i++)
            if (i >= 0 && i < (int)_poly.size())
            {
                if (force_fft)
                    _poly[i].do_fft();
                emplace(_transpose ? begin() : end(), _poly[i]);
            }
            else
                emplace(_transpose ? begin() : end());
    }

    Poly PolyMul::reciprocal(Poly& a, int size)
    {
        int i, j;
        int d = a.degree();
        Poly g(gw(), 1, true);
        Poly tmp(gw(), size/2, false);
        Poly tmp2(gw(), size/2, false);

        if (a[d].is_small() && abs(a[d].small()) == 1)
            g[0].set_small(a[d].small());
        else
        {
            Giant tmp;
            if (a[d].is_small())
                tmp = a[d].small();
            else
                tmp = a[d].value();
            g[0].own_set(gw(), tmp.inv(gw().N()));
        }

        for (i = 1; 2*i <= size; i <<= 1)
        {
            SubPolyFFT f(a, d - 2*i, 2*i);
            mul_half(g, f, tmp, i);
            tmp.emplace(tmp.begin());
            mul_half(g, tmp, tmp2, i);
            tmp.erase(tmp.begin());
            g.fft().reset();
            for (j = i - 1; j >= 0; j--)
            {
                g.emplace(g.begin());
                if (tmp2[j].is_small())
                    g.front().set_small(-tmp2[j].small());
                else
                {
                    g.front().own(gw());
                    gw().sub(g.front().value(), tmp2[j].value(), g.front().value(), GWADD_FORCE_NORMALIZE);
                }
            }
        }
        return g;
    }

    Poly PolyMul::mul(Poly& a, Poly& b)
    {
        Poly res(gw(), a.degree() + b.degree() + 1);
        mul(a, b, res);
        return res;
    }

    Poly PolyMul::mul_half(Poly& a, Poly& b, int half)
    {
        Poly res(gw(), half);
        mul_half(a, b, res, half);
        return res;
    }

    void PolyMul::mul(Poly& a, Poly& b, Poly& res)
    {
        if (a.preserve_fft())
            a.do_fft();
        if (b.preserve_fft())
            b.do_fft();
        res.clear_fft();
        res.set_zero();
        if (res.size() < a.degree() + b.degree() + 1)
            res.resize(a.degree() + b.degree() + 1);
        mul(a, b, res.data(), 0, (int)res.size());
    }

    void PolyMul::mul_half(Poly& a, Poly& b, Poly& res, int half)
    {
        GWASSERT(a.degree() + b.degree() < 3*half);
        if (a.preserve_fft())
            a.do_fft();
        if (b.preserve_fft())
            b.do_fft();
        res.clear_fft();
        res.set_zero();
        if (res.size() < half)
            res.resize(half);
        mul(a, b, res.data(), half, half);
    }

    void PolyMul::mul(const Poly& a, const Poly& b, PolyCoeff* res, int offset, int count)
    {
        int i, j, k;

        std::vector<int> smalls;
        std::vector<std::vector<GWNum*>> adds;
        std::vector<std::vector<GWNum*>> subs;
        std::vector<std::vector<std::pair<GWNum*, int>>> smallmuls;
        std::vector<std::vector<std::pair<GWNum*, GWNum*>>> muls;
        smalls.resize(count);
        adds.resize(count);
        subs.resize(count);
        smallmuls.resize(count);
        muls.resize(count);
        for (i = a.degree(); i >= 0; i--)
            for (j = b.degree(); j >= 0; j--)
                if (i + j >= offset && i + j < offset + count)
                {
                    k = i + j - offset;
                    if (a[i].is_zero() || b[j].is_zero())
                        continue;
                    else if (a[i].is_small() && b[j].is_small())
                        smalls[k] += a[i].small()*b[j].small();
                    else if (a[i].is_small() && a[i].small() == 1)
                        adds[k].push_back(&b[j].value());
                    else if (a[i].is_small() && a[i].small() == -1)
                        subs[k].push_back(&b[j].value());
                    else if (b[j].is_small() && b[j].small() == 1)
                        adds[k].push_back(&a[i].value());
                    else if (b[j].is_small() && b[j].small() == -1)
                        subs[k].push_back(&a[i].value());
                    else if (a[i].is_small())
                        smallmuls[k].emplace_back(b[j].has_fft() ? &b[j].fft() : &b[j].value(), a[i].small());
                    else if (b[j].is_small())
                        smallmuls[k].emplace_back(a[i].has_fft() ? &a[i].fft() : &a[i].value(), b[j].small());
                    else
                        muls[k].emplace_back(a[i].has_fft() ? &a[i].fft() : &a[i].value(), b[j].has_fft() ? &b[j].fft() : &b[j].value());
                }
        for (i = 0; i < count; i++)
        {
            if (abs(smalls[i]) < PolyCoeff::POLY_MAX_SMALL)
                res[i].set_small(smalls[i]);
            else
                res[i].own_set(gw(), smalls[i]);

            bool last = subs[i].empty() && smallmuls[i].empty() && muls[i].empty();
            for (auto it = adds[i].begin(); it != adds[i].end(); it++)
                if (res[i].is_zero())
                    res[i].own_set(gw(), **it);
                else
                {
                    res[i].own(gw());
                    gw().add(res[i].value(), **it, res[i].value(), last && it + 1 == adds[i].end() ? GWADD_FORCE_NORMALIZE : GWADD_DELAY_NORMALIZE);
                }
            last = smallmuls[i].empty() && muls[i].empty();
            for (auto it = subs[i].begin(); it != subs[i].end(); it++)
            {
                res[i].own(gw());
                gw().sub(res[i].value(), **it, res[i].value(), last && it + 1 == subs[i].end() ? GWADD_FORCE_NORMALIZE : GWADD_DELAY_NORMALIZE);
            }
            if (!res[i].is_small() && !gwnum_is_not_ffted(gw().gwdata(), *res[i].value()))
                gw().unfft(res[i].value(), res[i].value());
        }
        for (i = 0; i < count; i++)
        {
            bool last = muls[i].empty();
            for (auto it = smallmuls[i].begin(); it != smallmuls[i].end(); it++)
            {
                _tmp = it->second;
                gw().mul(*it->first, _tmp, _tmp, GWMUL_FFT_S1);
                if (res[i].is_zero())
                    res[i].own_swap(_tmp);
                else
                {
                    res[i].own(gw());
                    gw().add(res[i].value(), _tmp, res[i].value(), last && it + 1 == smallmuls[i].end() ? GWADD_FORCE_NORMALIZE : GWADD_DELAY_NORMALIZE);
                }
            }

            for (auto it = muls[i].begin(); it != muls[i].end(); it++)
            {
                auto itp = it + 1;
                if (!mul_safe(gw().gwdata(), (unnorms(**it->first) + 1)*(unnorms(**it->second) + 1) - 1, 0))
                    printf("1");
                //if (itp != muls[i].end() && !mul_safe(gw().gwdata(), (unnorms(**it->first) + 1)*(unnorms(**it->second) + 1) + (unnorms(**itp->first) + 1)*(unnorms(**itp->second) + 1) - 1, 0))
                //    printf("1");
                if (itp != muls[i].end() && mul_safe(gw().gwdata(), (unnorms(**it->first) + 1)*(unnorms(**it->second) + 1) + (unnorms(**itp->first) + 1)*(unnorms(**itp->second) + 1) - 1, 0))
                {
                    gw().mulmuladd(*it->first, *it->second, *itp->first, *itp->second, _tmp, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_FFT_S3 | GWMUL_FFT_S4);
                    it = itp;
                }
                else
                    gw().mul(*it->first, *it->second, _tmp, GWMUL_FFT_S1 | GWMUL_FFT_S2);
                if (res[i].is_zero())
                    res[i].own_swap(_tmp);
                else
                {
                    res[i].own(gw());
                    gw().add(res[i].value(), _tmp, res[i].value(), it + 1 == muls[i].end() ? GWADD_FORCE_NORMALIZE : GWADD_DELAY_NORMALIZE);
                }
            }
        }
    }
    
    void PolyMulKaratsuba::mul(Poly& a, Poly& b, Poly& res)
    {
        a.do_fft();
        b.do_fft();
        res.clear_fft();
        res.set_zero();
        if (res.size() < a.degree() + b.degree() + 1)
            res.resize(a.degree() + b.degree() + 1);
        karatsuba(a, b, res.data(), res.size());
    }

    void PolyMulKaratsuba::mul_half(Poly& a, Poly& b, Poly& res, int half)
    {
        GWASSERT(a.degree() + b.degree() < 3*half);
        a.do_fft();
        b.do_fft();
        res.clear_fft();
        res.set_zero();
        if (res.size() < half)
            res.resize(half);
        karatsuba_half(a, b, res.data(), res.size(), half);
    }

    void PolyMulKaratsuba::poly_sum(Poly& a, Poly& b, Poly& res)
    {
        int i;
        res.set_zero();
        int da = a.degree();
        int db = b.degree();
        for (i = 0; i <= da && i <= db; i++)
        {
            if (a[i].is_zero() && b[i].is_zero())
                continue;
            else if (a[i].is_zero())
                res[i] = b[i];
            else if (b[i].is_zero())
                res[i] = a[i];
            else if (a[i].is_small() && b[i].is_small())
            {
                int small = a[i].small() + b[i].small();
                if (abs(small) < PolyCoeff::POLY_MAX_SMALL)
                    res[i].set_small(small);
                else
                {
                    res[i].own_set(gw(), small);
                    res[i].do_fft();
                }
            }
            else if (a[i].is_small())
            {
                gw().add(b[i].value(), a[i].small(), _tmp);
                res[i].own_set(gw(), _tmp);
                res[i].do_fft();
            }
            else if (b[i].is_small())
            {
                gw().add(a[i].value(), b[i].small(), _tmp);
                res[i].own_set(gw(), _tmp);
                res[i].do_fft();
            }
            else
            {
                if (!gwnum_is_not_ffted(gw().gwdata(), *a[i].value()))
                    printf("1");
                if (!gwnum_is_not_ffted(gw().gwdata(), *b[i].value()))
                    printf("1");
                if (mul_safe(gw().gwdata(), (unnorms(*a[i].fft()) + 1 + unnorms(*b[i].fft()) + 1)*(unnorms(*a[i].fft()) + 1 + unnorms(*b[i].fft()) + 1) - 1, 0))
                {
                    gw().add(a[i].value(), b[i].value(), _tmp, GWADD_DELAY_NORMALIZE);
                    gw().add(a[i].fft(), b[i].fft(), _tmp_fft, GWADD_DELAY_NORMALIZE);
                }
                else
                {
                    gw().add(a[i].value(), b[i].value(), _tmp, GWADD_FORCE_NORMALIZE);
                    gw().fft(_tmp, _tmp_fft);
                }
                res[i].own_swap(_tmp, _tmp_fft);
            }
        }
        for (; i <= da; i++)
            res[i] = a[i];
        for (; i <= db; i++)
            res[i] = b[i];
    }

    void PolyMulKaratsuba::poly_add(const PolyCoeff& a, PolyCoeff& res)
    {
        if (a.is_zero())
            return;
        else if (res.is_zero() && a.is_small())
            res.set_small(a.small());
        else if (res.is_zero())
            res.own_swap(a.value());
        else if (a.is_small() && res.is_small())
        {
            int small = a.small() + res.small();
            if (abs(small) < PolyCoeff::POLY_MAX_SMALL)
                res.set_small(small);
            else
                res.own_set(gw(), small);
        }
        else if (a.is_small())
        {
            res.own(gw());
            gw().add(res.value(), a.small(), res.value());
        }
        else if (res.is_small())
        {
            int small = res.small();
            res.own_swap(a.value());
            gw().add(res.value(), small, res.value());
        }
        else
        {
            res.own(gw());
            gw().add(res.value(), a.value(), res.value(), GWADD_FORCE_NORMALIZE);
        }
    }

    void PolyMulKaratsuba::poly_sub(const PolyCoeff& a, PolyCoeff& res)
    {
        if (a.is_zero())
            return;
        else if (res.is_zero() && a.is_small())
            res.set_small(-a.small());
        else if (res.is_zero())
        {
            res.own(gw());
            gw().sub(res.value(), a.value(), res.value(), GWADD_GUARANTEED_OK);
        }
        else if (a.is_small() && res.is_small())
        {
            int small = res.small() - a.small();
            if (abs(small) < PolyCoeff::POLY_MAX_SMALL)
                res.set_small(small);
            else
                res.own_set(gw(), small);
        }
        else if (a.is_small())
        {
            res.own(gw());
            gw().sub(res.value(), a.small(), res.value());
        }
        else
        {
            res.own(gw());
            gw().sub(res.value(), a.value(), res.value(), GWADD_FORCE_NORMALIZE);
        }
    }

    Poly& PolyMulKaratsuba::get_tmp_poly(int level, int id, int size)
    {
        if (_tmp_polys.size() <= level)
            _tmp_polys.resize(level + 1);
        std::vector<std::unique_ptr<Poly>>& tmp_level = _tmp_polys[level];
        if (tmp_level.size() <= id)
            tmp_level.resize(id + 1);
        std::unique_ptr<Poly>& tmp_ptr = tmp_level[id];
        if (!tmp_ptr)
            tmp_ptr.reset(new Poly(gw(), 0));
        Poly* tmp_poly = tmp_ptr.get();
        for (auto it = tmp_poly->begin(); it != tmp_poly->end(); it++)
            it->set_zero();
        if (tmp_poly->size() < size)
            tmp_poly->resize(size);
        return *tmp_poly;
    }

    void PolyMulKaratsuba::karatsuba(const Poly& a, const Poly& b, PolyCoeff* res, size_t size, int level)
    {
        int i;
        int da = a.degree();
        int db = b.degree();
        if (da == -1 || db == -1)
            return;
        if (da < 4 || db < 4)
        {
            GWASSERT(da + db < size);
            PolyMul::mul(a, b, res, 0, da + db + 1);
            return;
        }
        int d = (da + 1)/2;

        SubPoly a0(a, 0, d);
        SubPoly a1(a, d);
        SubPoly b0(b, 0, d);
        SubPoly b1(b, d);

        Poly& poly0 = get_tmp_poly(level, 0, 2*d - 1);
        Poly& poly1 = get_tmp_poly(level, 1, da + db - 2*d + 1);
        karatsuba(a0, b0, poly0.data(), poly0.size(), level + 1);
        karatsuba(a1, b1, poly1.data(), poly1.size(), level + 1);

        Poly& as = get_tmp_poly(level, 2, std::max(da - d, d - 1) + 1);
        Poly& bs = get_tmp_poly(level, 3, std::max(db - d, d - 1) + 1);
        poly_sum(a0, a1, as);
        poly_sum(b0, b1, bs);
        karatsuba(as, bs, res + d, size - d, level + 1);

        GWASSERT(poly0.degree() + d < size);
        GWASSERT(poly1.degree() + 2*d < size);
        for (i = poly0.degree(); i >= 0; i--)
            poly_sub(poly0[i], res[i + d]);
        for (i = poly1.degree(); i >= 0; i--)
            poly_sub(poly1[i], res[i + d]);
        for (i = poly0.degree(); i >= 0; i--)
            poly_add(poly0[i], res[i]);
        for (i = poly1.degree(); i >= 0; i--)
            poly_add(poly1[i], res[i + 2*d]);
    }

    void PolyMulKaratsuba::karatsuba_half(const Poly& a, const Poly& b, PolyCoeff* res, size_t size, int half, int level)
    {
        int i;
        int da = a.degree();
        int db = b.degree();
        if (da == -1 || db == -1)
            return;
        if (da < 4 || db < 4)
        {
            GWASSERT(half <= size);
            PolyMul::mul(a, b, res, half, half);
            return;
        }

        Poly& poly0 = get_tmp_poly(level, 0, half);
        Poly& poly1 = get_tmp_poly(level, 1, half - 1);
        karatsuba_halfother(SubPoly(a, 0, half), SubPoly(b, 1, half), res, size, half, level + 1);
        karatsuba_halfother(SubPoly(a, half, half, true), SubPoly(b, 0, half, true), poly0.data(), poly0.size(), half, level + 1);
        karatsuba_halfother(SubPoly(a, 0, half - 1, true), SubPoly(b, half + 1, half - 1, true), poly1.data(), poly1.size(), half - 1, level + 1);

        GWASSERT(poly0.degree() < half);
        GWASSERT(poly1.degree() < half);
        GWASSERT(half <= size);
        for (i = poly0.degree(); i >= 0; i--)
            poly_add(poly0[i], res[half - 1 - i]);
        for (i = poly1.degree(); i >= 0; i--)
            poly_add(poly1[i], res[half - 1 - i]);
    }

    void PolyMulKaratsuba::karatsuba_halfother(const Poly& a, const Poly& b, PolyCoeff* res, size_t size, int half, int level)
    {
        int i;
        int da = a.degree();
        int db = b.degree();
        if (da == -1 || db == -1)
            return;
        if (da < 4 || db < 4)
        {
            GWASSERT(half <= size);
            PolyMul::mul(a, b, res, half - 1, half);
            return;
        }
        int d = half/2;

        Poly& poly01 = get_tmp_poly(level, 0, d);
        Poly& poly10 = get_tmp_poly(level, 1, d);
        karatsuba(SubPoly(a, d, half - d), SubPoly(b, d, half - d), res + (1 - (half & 1)), size - (1 - (half & 1)), level + 1);
        karatsuba_halfother(SubPoly(a, 0, d), SubPoly(b, half - d, d), poly01.data(), poly01.size(), d, level + 1);
        karatsuba_halfother(SubPoly(b, 0, d), SubPoly(a, half - d, d), poly10.data(), poly10.size(), d, level + 1);

        GWASSERT(poly01.degree() < (int)size);
        GWASSERT(poly10.degree() < (int)size);
        for (i = poly01.degree(); i >= 0; i--)
            poly_add(poly01[i], res[i]);
        for (i = poly10.degree(); i >= 0; i--)
            poly_add(poly10[i], res[i]);
    }
    
    PolyMulFFT::PolyMulFFT(GWArithmetic& gw, int size) : PolyMul(gw), _size(size), _gwstate(), _gwpoly(_gwstate)
    {
        _coeff_size = (int)(2*gw.state().bit_length + 16 + 31)/32;
        _gwstate.safety_margin = 0.25;
        _gwstate.setup(_coeff_size*size*32);
        _tmp_g.arithmetic().alloc(_tmp_g, _gwstate.giants.capacity());
    }

    PolyMulFFT::~PolyMulFFT()
    {
    }

    void PolyMulFFT::poly_fft(Poly& a)
    {
        if (a.fft() && &a.fft()->arithmetic() == &gwpoly())
            return;
        a.fft().reset(new GWNum(gwpoly()));

        int i;
        memset(_tmp_g.data(), 0, 4*_coeff_size*(a.degree() + 1));
        for (i = a.degree(); i >= 0; i--)
        {
            uint32_t* coeff = _tmp_g.data() + _coeff_size*i;
            if (a[i].is_small())
                coeff[0] = a[i].small();
            else if (gwtobinary(gw().gwdata(), *a[i].value(), coeff, _coeff_size) < 0)
                throw InvalidFFTDataException();
        }
        binarytogw(gwdata(), _tmp_g.data(), _coeff_size*(a.degree() + 1), **a.fft());
    }

    extern "C" void emulate_mod(gwhandle *gwdata, gwnum	s);

    void PolyMulFFT::reduce_coeff(uint32_t* data, int count, PolyCoeff& res)
    {
        int reduced_size = (int)(gw().state().bit_length + 31)/32;
        giantstruct g0;
        setmaxsize(&g0, reduced_size);
        giantstruct g1;
        setmaxsize(&g1, reduced_size);

        if (gw().gwdata()->GENERAL_MOD)
        {
            setmaxsize(&g0, _coeff_size);
            g0.n = data;
            g0.sign = count;
            while (g0.sign > 0 && g0.n[g0.sign - 1] == 0)
                g0.sign--;
            gianttogw(gw().gwdata(), &g0, *_tmp);
            emulate_mod(gw().gwdata(), *_tmp);
            res.own_swap(_tmp);
        }
        else if (gw().gwdata()->k == 1.0 && gw().gwdata()->b == 2 && gw().gwdata()->c == 1 && (gw().gwdata()->n & 31) == 0)
        {
            g0.n = data;
            g0.sign = reduced_size;
            if (g0.sign > count)
                g0.sign = count;
            while (g0.sign > 0 && g0.n[g0.sign - 1] == 0)
                g0.sign--;
            g1.n = data + reduced_size;
            g1.sign = reduced_size;
            if (g1.sign > count - reduced_size)
                g1.sign = count - reduced_size;
            if (g1.sign < 0)
                g1.sign = 0;
            while (g1.sign > 0 && g1.n[g1.sign - 1] == 0)
                g1.sign--;
            subg(&g1, &g0);
            setmaxsize(&g0, 2*reduced_size);
            if (2*reduced_size < count)
                uladdg(data[2*reduced_size], &g0);
            if (g0.sign < 0)
            {
                g0.sign = reduced_size;
                for (int j = 0; j < g0.sign; j++)
                    g0.n[j] = ~g0.n[j];
                while (g0.sign > 0 && g0.n[g0.sign - 1] == 0)
                    g0.sign--;
                uladdg(2, &g0);
            }
            if (g0.sign == 0)
                res.set_zero();
            else if (g0.sign == 1 && g0.n[0] < PolyCoeff::POLY_MAX_SMALL)
                res.set_small(g0.n[0]);
            else
            {
                gianttogw(gw().gwdata(), &g0, *_tmp);
                res.own_swap(_tmp);
            }
        }
        else
        {
            setmaxsize(&g0, _coeff_size);
            g0.n = data;
            g0.sign = count;
            while (g0.sign > 0 && g0.n[g0.sign - 1] == 0)
                g0.sign--;
            specialmodg(gw().gwdata(), &g0);
            gianttogw(gw().gwdata(), &g0, *_tmp);
            res.own_swap(_tmp);
        }
    }

    void PolyMulFFT::mul(Poly& a, Poly& b, Poly& res)
    {
        GWASSERT(a.degree() + b.degree() < size());
        res.clear_fft();
        res.set_zero();
        if (res.size() < a.degree() + b.degree() + 1)
            res.resize(a.degree() + b.degree() + 1);

        poly_fft(a);
        poly_fft(b);
        if (a.preserve_fft() && b.preserve_fft())
        {
            GWNum tmp(gwpoly());
            gwpoly().mul(*a.fft(), *b.fft(), tmp, 0);
            _tmp_g = tmp;
        }
        else if (b.preserve_fft())
        {
            gwpoly().mul(*a.fft(), *b.fft(), *a.fft(), 0);
            _tmp_g = *a.fft();
            a.fft().reset();
        }
        else
        {
            gwpoly().mul(*a.fft(), *b.fft(), *b.fft(), 0);
            _tmp_g = *b.fft();
            b.fft().reset();
            if (!a.preserve_fft())
                a.fft().reset();
        }

        int len = _tmp_g.size();
        for (int i = 0; len > 0; i++, len -= _coeff_size)
            reduce_coeff(_tmp_g.data() + _coeff_size*i, len > _coeff_size ? _coeff_size : len, res[i]);
    }

    void PolyMulFFT::mul_half(Poly& a, Poly& b, Poly& res, int half)
    {
        GWASSERT(a.degree() + b.degree() < 3*half);
        GWASSERT(2*half <= size());
        res.clear_fft();
        res.set_zero();
        if (res.size() < half)
            res.resize(half);

        poly_fft(a);
        poly_fft(b);
        if (a.preserve_fft() && b.preserve_fft())
        {
            GWNum tmp(gwpoly());
            gwpoly().mul(*a.fft(), *b.fft(), tmp, 0);
            _tmp_g = tmp;
        }
        else if (b.preserve_fft())
        {
            gwpoly().mul(*a.fft(), *b.fft(), *a.fft(), 0);
            _tmp_g = *a.fft();
            a.fft().reset();
        }
        else
        {
            gwpoly().mul(*a.fft(), *b.fft(), *b.fft(), 0);
            _tmp_g = *b.fft();
            b.fft().reset();
            if (!a.preserve_fft())
                a.fft().reset();
        }

        int len = _tmp_g.size() - _coeff_size*half;
        for (int i = 0; len > 0 && i < half; i++, len -= _coeff_size)
            reduce_coeff(_tmp_g.data() + _coeff_size*(i + half), len > _coeff_size ? _coeff_size : len, res[i]);
    }

    PolyMul& PolyMulOptimal::get_optimal(const Poly& a, const Poly& b, int half)
    {
        int da = a.degree();
        int db = b.degree();
        int sda, sdb;
        for (sda = da; sda >= 0 && a[sda].is_small(); sda--);
        for (sdb = db; sdb >= 0 && b[sdb].is_small(); sdb--);
        if (std::min(sda, sdb) < 4)
        {
            if (!_mul)
                _mul.reset(new PolyMul(gw()));
            return *_mul;
        }
        if (std::max(sda, sdb) < 32 || (!gw().gwdata()->GENERAL_MOD && std::max(sda, sdb) < 1024))
        {
            if (!_karatsuba)
                _karatsuba.reset(new PolyMulKaratsuba(gw()));
            return *_karatsuba;
        }
        int i = da + db + 1;
        if (half > 0)
            i = 2*half;
        if (_ffts.size() < i + 1)
            _ffts.resize(i + 1);
        if (!_ffts[i])
            _ffts[i].reset(new PolyMulFFT(gw(), i));
        return *_ffts[i];
    }

    void PolyMulOptimal::mul(Poly& a, Poly& b, Poly& res)
    {
        PolyMul& optimal = get_optimal(a, b, 0);
        optimal.mul(a, b, res);
    }

    void PolyMulOptimal::mul_half(Poly& a, Poly& b, Poly& res, int half)
    {
        PolyMul& optimal = get_optimal(a, b, half);
        optimal.mul_half(a, b, res, half);
    }

    /*PolyMulPrime::PolyMulPrime(GWArithmetic& gw, int size) : PolyMul(gw), _size(size)
    {
        int i, j;
        Giant tmp;
        int len, depth, a;
        std::vector<GWNum> roots;

        GWASSERT((size & (size - 1)) == 0);
        for (depth = 0; (1 << depth) < size; depth++);
        a = 3; // Jacobi check!!!
        roots.emplace_back(gw);
        roots.front() = a;
        gw.setmulbyconst(a);
        gwset_carefully_count(gw.gwdata(), 30);
        len = gw.N().bitlen() - 1;
        for (j = 1; !gw.N().bit(j); j++);
        for (i = 1; i <= len - j; i++)
            gw.square(roots.front(), roots.front(), (gw.N().bit(len - i) ? GWMUL_MULBYCONST : 0) | GWMUL_STARTNEXTFFT);
        for (; i <= len - depth; i++)
            gw.square(roots.front(), roots.front(), i < len - depth ? GWMUL_STARTNEXTFFT : 0);
        for (; i <= len; i++)
        {
            roots.emplace(roots.begin(), gw);
            gw.carefully().square(roots[1], roots[0], 0);
        }
        if (roots.size() != depth + 1)
            throw ArithmeticException("Not enough roots of unity.");
        if (roots[0] != 1)
            throw ArithmeticException("The number is not prime.");
        if (roots[1] != -1)
            throw ArithmeticException("Invalid root. Jacobi check missing?");
        _roots.emplace_back(std::move(roots[0]));
        for (i = 2; i <= depth; i++)
        {
            _roots.emplace_back(std::move(roots[i]));
            for (j = 1; i > 2 && j < (1 << (i - 2)); j++)
            {
                _roots.emplace_back(gw);
                gw.carefully().mul(_roots[j], _roots[_roots.size() - 1 - j], _roots.back(), 0);
            }
        }
        GWASSERT(_roots.size() == size/2);
        for (i = 1; i < _roots.size(); i++)
        {
            if (i%2 == 0)
                GWASSERT(_roots[i]*_roots[i] - _roots[i >> 1] == 0);
            if (i%2 == 1)
                GWASSERT(_roots[i]*_roots[i] + _roots[i >> 1] == 0);
            for (j = 0; j < _roots.size(); j++)
                GWASSERT(i == j || (_roots[i] != _roots[j] && _roots[i] + _roots[j] != 0));
        }
        for (auto it = _roots.begin(); it != _roots.end(); it++)
        {
            _inv_roots.emplace_back(gw);
            tmp = *it;
            if (it == _roots.begin())
                tmp = size;
            _inv_roots.back() = inv(tmp, gw.N());
        }
    }

    Poly PolyMulPrime::mul(Poly&a, Poly&b)
    {
        Poly fft_a(gw(), size());
        Poly fft_b(gw(), size());
        transform(a, fft_a);
        transform(b, fft_b);
        for (int i = 0; i < size(); i++)
            fft_a.a(i) *= fft_b.a(i);
        inv_transform(fft_a);
        return fft_a;
    }

    Poly PolyMulPrime::mul_half(Poly&a, Poly&b, int half)
    {
        Poly fft_a(gw(), size());
        Poly fft_b(gw(), size());
        transform(a, fft_a);
        transform(b, fft_b);
        for (int i = 0; i < size(); i++)
            fft_a.a(i) *= fft_b.a(i);
        inv_transform(fft_a);
        GWASSERT(0);
        return fft_a;
    }

    void PolyMulPrime::transform(Poly& src, Poly& dst)
    {
        transform(src, dst, 0, _size, 0);
    }

    void PolyMulPrime::transform(Poly& src, Poly& dst, int offset, int count, int root)
    {
        if (count == 1)
            return;
        int i, m;
        m = count/2;
        for (i = 0; i < m; i++)
        {
            if (!src[offset + i] && !src[offset + m + i])
            {
                dst[offset + i].reset();
                dst[offset + m + i].reset();
                continue;
            }
            if (!src[offset + m + i])
            {
                dst.a(offset + i) = src.a(offset + i);
                dst.a(offset + m + i) = src.a(offset + i);
                continue;
            }
            if (root != 0)
                gw().mul(src.a(offset + m + i), _roots[root], dst.a(offset + m + i), GWMUL_FFT_S2);
            else if (&dst.a(offset + m + i) != &src.a(offset + m + i))
                dst.a(offset + m + i) = src.a(offset + m + i);
            if (!src[offset + i])
            {
                dst.a(offset + i) = dst.a(offset + m + i);
                gw().neg(dst.a(offset + m + i), dst.a(offset + m + i));
                continue;
            }
            gw().addsub(src.a(offset + i), dst.a(offset + m + i), dst.a(offset + i), dst.a(offset + m + i), GWADD_FORCE_NORMALIZE);
        }
        transform(dst, dst, offset, count/2, root*2);
        transform(dst, dst, offset + m, count/2, root*2 + 1);
    }

    void PolyMulPrime::inv_transform(Poly& dst)
    {
        inv_transform(dst, 0, _size, 0);
        for (auto it = dst.begin(); it != dst.end(); it++)
            if (*it)
                gw().mul(_inv_roots[0], **it, **it, GWMUL_FFT_S1);
    }

    void PolyMulPrime::inv_transform(Poly& dst, int offset, int count, int root)
    {
        if (count == 1)
            return;
        int i, m;
        m = count/2;
        inv_transform(dst, offset, count/2, root*2);
        inv_transform(dst, offset + m, count/2, root*2 + 1);
        for (i = 0; i < m; i++)
        {
            gw().addsub(dst.a(offset + i), dst.a(offset + m + i), dst.a(offset + i), dst.a(offset + m + i), GWADD_FORCE_NORMALIZE);
            if (root != 0)
                gw().mul(dst.a(offset + m + i), _inv_roots[root], dst.a(offset + m + i), GWMUL_FFT_S2);
        }
    }*/
}
