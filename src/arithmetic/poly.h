#pragma once

#include <memory>
#include <vector>

#include "arithmetic.h"

namespace arithmetic
{
    class PolyCoeff
    {
    public:
        static const int POLY_NOT_SMALL = 0x01;
        static const int POLY_OWN = 0x02;
        static const int POLY_FFT = 0x04;
        static const int POLY_MAX_SMALL = 32;

    public:
        PolyCoeff() : _flag(0), _value(nullptr), _fft(nullptr) { }
        PolyCoeff(GWArithmetic& gw) : _flag(POLY_OWN), _value(new GWNum(gw)), _fft(nullptr) { }
        ~PolyCoeff() { reset(); }
        PolyCoeff(PolyCoeff&& a) noexcept : _flag(a._flag), _value(a._value), _fft(a._fft) { a._flag = 0; a._value = nullptr; a._fft = nullptr; }
        PolyCoeff& operator = (PolyCoeff&& a) noexcept { reset(); _flag = a._flag; _value = a._value; _fft = a._fft; a._flag = 0; a._value = nullptr; a._fft = nullptr; return *this; }

        bool operator == (const PolyCoeff& a) { if (is_small() && a.is_small()) return small() == a.small(); if (is_small() && !a.is_small()) return small() == a.value();  if (!is_small() && a.is_small()) return value() == a.small(); return value() == a.value(); }
        bool operator != (const PolyCoeff& a) { return !(*this == a); }

        void set_zero() { _flag &= (POLY_OWN | POLY_FFT); }
        void set_small(int value) { GWASSERT(abs(value) < POLY_MAX_SMALL); _flag = (_flag & (POLY_OWN | POLY_FFT)) | (value << 8); }
        void reset() { reset_fft(); if (_flag & POLY_OWN) delete _value; _value = nullptr; _flag &= ~(POLY_OWN | POLY_NOT_SMALL); }
        void reset_fft() { if ((_flag & POLY_OWN) && _fft != nullptr) delete _fft; _fft = nullptr; _flag &= ~POLY_FFT; }
        void clear_fft() { if (!(_flag & POLY_OWN)) _fft = nullptr; _flag &= ~POLY_FFT; }

        void set(const PolyCoeff& a)
        {
            if (!has_own())
                _flag = (a._flag & ~POLY_OWN);
            else
                _flag = (a._flag & ~POLY_FFT) | (_flag & POLY_FFT) | POLY_OWN;
            if (!is_small())
            {
                if (!has_own())
                {
                    _value = a._value;
                    if (has_fft())
                        _fft = a._fft;
                    else
                        _fft = nullptr;
                }
                else
                {
                    value() = a.value();
                    if (has_fft() && a.has_fft())
                        fft() = a.fft();
                    if (has_fft() && !a.has_fft())
                        value().arithmetic().fft(value(), fft());
                }
            }
        }

        void set(GWNum* a)
        {
            if (!has_own())
                _value = a;
            else
            {
                value() = *a;
                if (has_fft())
                    value().arithmetic().fft(value(), fft());
            }
            _flag |= POLY_NOT_SMALL;
        }

        void do_fft()
        {
            if (!has_own() || has_fft())
                return;
            if (_fft == nullptr)
                _fft = new GWNum(value().arithmetic());
            _flag |= POLY_FFT;
            if (!is_small())
                value().arithmetic().fft(value(), fft());
        }

        void own(GWArithmetic& gw)
        {
            GWNum* newval = has_own() ? _value : new GWNum(gw);
            if (is_small())
                *newval = small();
            else if (!has_own())
                *newval = *_value;
            _value = newval;
            if (has_fft())
            {
                newval = has_own() ? _fft : new GWNum(gw);
                if (is_small())
                    gw.fft(value(), *newval);
                else if (!has_own())
                    *newval = *_fft;
                _fft = newval;
            }
            _flag = POLY_OWN | POLY_NOT_SMALL | (_flag & POLY_FFT);
        }

        template<class T>
        void own_set(GWArithmetic& gw, const T& value)
        {
            if (!has_own())
                _value = new GWNum(gw);
            if (!has_own() && has_fft())
                _fft = new GWNum(gw);
            *_value = value;
            if (has_fft())
                gw.fft(*_value, *_fft);
            _flag = POLY_OWN | POLY_NOT_SMALL | (_flag & POLY_FFT);
        }

        void own_swap(GWNum& value)
        {
            if (!has_own())
                _value = new GWNum(value.arithmetic());
            if (!has_own() && has_fft())
                _fft = new GWNum(value.arithmetic());
            swap(value, *_value);
            if (has_fft())
                value.arithmetic().fft(*_value, *_fft);
            _flag = POLY_OWN | POLY_NOT_SMALL | (_flag & POLY_FFT);
        }

        void own_swap(GWNum& value, GWNum& fft)
        {
            if (!has_own())
                _value = new GWNum(value.arithmetic());
            if (!has_own() || _fft == nullptr)
                _fft = new GWNum(fft.arithmetic());
            swap(value, *_value);
            swap(fft, *_fft);
            _flag = POLY_OWN | POLY_NOT_SMALL | POLY_FFT;
        }

        bool is_zero() const { return (_flag & ~(POLY_OWN | POLY_FFT)) == 0; }
        bool is_small() const { return (_flag & POLY_NOT_SMALL) == 0; }
        bool has_own() const { return (_flag & POLY_OWN) != 0; }
        bool has_fft() const { return (_flag & POLY_FFT) != 0; }
        int32_t small() const { return _flag >> 8; }
        GWNum& value() const { return *_value; }
        GWNum& fft() const { return *_fft; }

    private:
        int _flag;
        GWNum* _value;
        GWNum* _fft;
    };

    class Poly : public std::vector<PolyCoeff>
    {
    public:
        Poly(GWArithmetic& gw, int size, bool preserve_fft = false) : _gw(gw), _preserve_fft(preserve_fft)
        {
            for (int i = 0; i < size; i++)
                emplace_back();
        }
        Poly(Poly&& a) noexcept : vector(std::move(a)), _gw(a.gw()), _fft(std::move(a._fft)), _preserve_fft(a._preserve_fft){ }
        virtual ~Poly() { }
        Poly& operator = (Poly&& a) noexcept { std::vector<PolyCoeff>::operator=(std::move(a)); _fft = std::move(a._fft); _preserve_fft = a._preserve_fft; return *this; }

        void set_zero() { for (auto it = begin(); it != end(); it++) it->set_zero(); }
        virtual void do_fft();
        void clear_fft();
        GWNum eval(GWNum& x);
        Poly mul(Poly& b);
        Poly mul_half(Poly& b, int half);
        Poly reciprocal(int size);

        GWArithmetic& gw() const { return _gw; }
        int degree() const { for (size_t i = size(); i > 0; i--) if (!at(i - 1).is_zero()) return (int)(i - 1); return -1; }
        bool empty() const { return degree() == -1; }
        std::unique_ptr<GWNum>& fft() { return _fft; }
        bool preserve_fft() const { return _preserve_fft; }
        void set_preserve_fft(bool value) { _preserve_fft = value; }

    private:
        GWArithmetic& _gw;
        std::unique_ptr<GWNum> _fft;
        bool _preserve_fft;
    };

    class SubPoly : public Poly
    {
    public:
        SubPoly(const Poly& poly, int offset, int count = -1, bool transpose = false);
    };

    class SubPolyFFT : public Poly
    {
    public:
        SubPolyFFT(Poly& poly, int offset, int count = -1, bool transpose = false) : Poly(poly.gw(), 0), _poly(poly), _offset(offset), _count(count), _transpose(transpose) { init(false); }

        void init(bool force_fft);
        virtual void do_fft() override { init(true); }

    private:
        Poly& _poly;
        int _offset;
        int _count;
        bool _transpose;
    };

    class PolyMul
    {
    public:
        PolyMul(GWArithmetic& gw) : _gw(gw), _tmp(gw) { }
        virtual ~PolyMul() { }

        virtual Poly mul(Poly& a, Poly& b);
        virtual Poly mul_half(Poly& a, Poly& b, int half);
        virtual void mul(Poly& a, Poly& b, Poly& res);
        virtual void mul_half(Poly& a, Poly& b, Poly& res, int half);

        void mul(const Poly& a, const Poly& b, PolyCoeff* res, int offset, int count);

        Poly reciprocal(Poly& a, int size);

        GWArithmetic& gw() { return _gw; }

    protected:
        GWArithmetic& _gw;
        GWNum _tmp;
    };

    class PolyMulKaratsuba : public PolyMul
    {
    public:
        PolyMulKaratsuba(GWArithmetic& gw) : PolyMul(gw), _tmp_fft(gw) { }

        using PolyMul::mul;
        using PolyMul::mul_half;
        virtual void mul(Poly& a, Poly& b, Poly& res) override;
        virtual void mul_half(Poly& a, Poly& b, Poly& res, int half) override;

        void karatsuba(const Poly& a, const Poly& b, PolyCoeff* res, size_t size, int level = 0);
        void karatsuba_half(const Poly& a, const Poly& b, PolyCoeff* res, size_t size, int half, int level = 0);
        void karatsuba_halfother(const Poly& a, const Poly& b, PolyCoeff* res, size_t size, int half, int level = 0);

    private:
        void poly_sum(Poly& a, Poly& b, Poly& res);
        void poly_add(const PolyCoeff& a, PolyCoeff& res);
        void poly_sub(const PolyCoeff& a, PolyCoeff& res);
        Poly& get_tmp_poly(int level, int id, int size);

    private:
        GWNum _tmp_fft;
        std::vector<std::vector<std::unique_ptr<Poly>>> _tmp_polys;
    };

    class PolyMulFFT : public PolyMul
    {
    public:
        PolyMulFFT(GWArithmetic& gw, int size);
        ~PolyMulFFT();

        using PolyMul::mul;
        using PolyMul::mul_half;
        virtual void mul(Poly& a, Poly& b, Poly& res) override;
        virtual void mul_half(Poly& a, Poly& b, Poly& res, int half) override;

        int size() { return _size; }
        gwhandle* gwdata() { return _gwstate.gwdata(); }
        GWArithmetic& gwpoly() { return _gwpoly; }

    private:
        void poly_fft(Poly& a);
        void reduce_coeff(uint32_t* data, int count, PolyCoeff& res);

    private:
        int _size;
        int _coeff_size;
        GWState _gwstate;
        GWArithmetic _gwpoly;
        Giant _tmp_g;
    };

    class PolyMulOptimal : public PolyMul
    {
    public:
        PolyMulOptimal(GWArithmetic& gw) : PolyMul(gw) { }

        using PolyMul::mul;
        using PolyMul::mul_half;
        virtual void mul(Poly& a, Poly& b, Poly& res) override;
        virtual void mul_half(Poly& a, Poly& b, Poly& res, int half) override;

        PolyMul& get_optimal(const Poly& a, const Poly& b, int half);

    private:
        std::unique_ptr<PolyMul> _mul;
        std::unique_ptr<PolyMulKaratsuba> _karatsuba;
        std::vector<std::unique_ptr<PolyMulFFT>> _ffts;
    };

    /*class PolyMulPrime : public PolyMul
    {
    public:
        PolyMulPrime(GWArithmetic& gw, int size);

        virtual Poly mul(Poly&a, Poly&b) override;
        virtual Poly mul_half(Poly&a, Poly&b, int half) override;

        void transform(Poly& src, Poly& dst);
        void inv_transform(Poly& dst);

        int size() const { return _size; }

    private:
        void transform(Poly& src, Poly& dst, int offset, int count, int depth);
        void inv_transform(Poly& dst, int offset, int count, int depth);

    private:
        int _size;
        std::vector<GWNum> _roots;
        std::vector<GWNum> _inv_roots;
    };*/
}
