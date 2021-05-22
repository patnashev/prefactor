
#include "gwnum.h"
#include "montgomery.h"

namespace arithmetic
{
    void MontgomeryArithmetic::copy(const EdY& a, EdY& res)
    {
        if (a.Y)
        {
            if (!res.Y)
                res.Y.reset(new GWNum(gw()));
            *res.Y = *a.Y;
        }
        else
            res.Y.reset();
        if (a.Z)
        {
            if (!res.Z)
                res.Z.reset(new GWNum(gw()));
            *res.Z = *a.Z;
        }
        else
            res.Z.reset();
        if (a.ZpY)
        {
            if (!res.ZpY)
                res.ZpY.reset(new GWNum(gw()));
            *res.ZpY = *a.ZpY;
        }
        else
            res.ZpY.reset();
        if (a.ZmY)
        {
            if (!res.ZmY)
                res.ZmY.reset(new GWNum(gw()));
            *res.ZmY = *a.ZmY;
        }
        else
            res.ZmY.reset();
    }

    void MontgomeryArithmetic::move(EdY&& a, EdY& res)
    {
        res.Y = std::move(a.Y);
        res.Z = std::move(a.Z);
        res.ZpY = std::move(a.ZpY);
        res.ZmY = std::move(a.ZmY);
    }

    void MontgomeryArithmetic::init(EdY& a)
    {
        a.Y.reset();
        a.Z.reset();
        a.ZpY.reset();
        a.ZmY.reset();
    }

    void MontgomeryArithmetic::init(const EdPoint& a, EdY& res)
    {
        res.Y.reset(new GWNum(gw()));
        *res.Y = *a.Y;
        if (a.Z)
        {
            res.Z.reset(new GWNum(gw()));
            *res.Z = *a.Z;
        }
    }

    void MontgomeryArithmetic::add(EdY& a, EdY& b, EdY& a_minus_b, EdY& res)
    {
        bool safe1 = square_safe(gw().gwdata(), 1);
        bool safe11 = mul_safe(gw().gwdata(), 1, 1);

        if (!a_minus_b.Y)
        {
            dbl(a, res);
            return;
        }
        if (!res.Y)
            res.Y.reset(new GWNum(gw()));
        if (!res.Z)
            res.Z.reset(new GWNum(gw()));
        if (!res.ZpY)
            res.ZpY.reset(new GWNum(gw()));
        if (!res.ZmY)
            res.ZmY.reset(new GWNum(gw()));
        if (!a_minus_b.ZpY)
        {
            a_minus_b.ZpY.reset(new GWNum(gw()));
            a_minus_b.ZmY.reset(new GWNum(gw()));
            if (a_minus_b.Z)
                *a_minus_b.ZpY = *a_minus_b.Z;
            else
                *a_minus_b.ZpY = 1;
            *a_minus_b.ZmY = *a_minus_b.Y;
            gw().addsub(*a_minus_b.ZpY, *a_minus_b.ZmY, *a_minus_b.ZpY, *a_minus_b.ZmY, GWADD_DELAYNORM_IF(safe11));
        }

        if (b.Z)
            gw().mul(*a.Y, *b.Z, *res.Y, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT_IF(safe1));
        else if (&a != &res)
            gw().copy(*a.Y, *res.Y);
        if (a.Z)
            gw().mul(*a.Z, *b.Y, *res.Z, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT_IF(safe1));
        else
            gw().copy(*b.Y, *res.Z);
        gw().addsub(*res.Y, *res.Z, *res.Z, *res.Y, GWADD_DELAYNORM_IF(safe1));
        gw().square(*res.Z, *res.Z, GWMUL_STARTNEXTFFT);
        gw().square(*res.Y, *res.Y, GWMUL_STARTNEXTFFT);
        gw().mulmuladd(*a_minus_b.ZmY, *res.Z, *a_minus_b.ZpY, *res.Y, *res.ZpY, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_FFT_S3 | GWMUL_FFT_S4 | GWMUL_STARTNEXTFFT_IF(safe11));
        gw().mulmulsub(*a_minus_b.ZmY, *res.Z, *a_minus_b.ZpY, *res.Y, *res.ZmY, GWMUL_STARTNEXTFFT_IF(safe11 && (b.Z || safe1)));
        if (gwnum_is_partially_ffted(gw().gwdata(), **res.ZpY))
            gw().fft(*res.ZpY, *res.ZpY);
        if (gwnum_is_partially_ffted(gw().gwdata(), **res.ZmY))
            gw().fft(*res.ZmY, *res.ZmY);
        gw().copy(*res.ZpY, *res.Z);
        gw().copy(*res.ZmY, *res.Y);
        gw().addsub(*res.ZpY, *res.ZmY, *res.ZpY, *res.ZmY, GWADD_DELAYNORM_IF(safe11));
    }

    void MontgomeryArithmetic::dbl(EdY& a, EdY& res)
    {
        bool safe1 = square_safe(gw().gwdata(), 1);
        bool safe11 = mul_safe(gw().gwdata(), 1, 1);

        if (!res.Y)
            res.Y.reset(new GWNum(gw()));
        if (!res.Z)
            res.Z.reset(new GWNum(gw()));
        if (!res.ZpY)
            res.ZpY.reset(new GWNum(gw()));
        if (!res.ZmY)
            res.ZmY.reset(new GWNum(gw()));

        // d = 1 - 1/Ad4
        // t1 = zz*(yy - d*yy)
        // t2 = (zz - yy)*(zz - d*yy)
        // y_2 = t1 - t2
        // z_2 = t1 + t2

        gw().square(*a.Y, *res.Y, GWMUL_STARTNEXTFFT_IF(safe11)); // yy
        gw().square(*a.Z, *res.Z, GWMUL_STARTNEXTFFT_IF(safe11)); // zz
        gw().sub(*res.Z, *res.Y, *res.ZmY, GWADD_DELAYNORM_IF(safe11)); // zz - yy
        gw().mul(ed_d(), *res.Y, *res.ZpY, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT); // d*yy
        gw().submul(*res.Z, *res.ZpY, *res.ZmY, *res.ZmY, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT_IF(safe1)); // t2
        gw().submul(*res.Y, *res.ZpY, *res.Z, *res.Z, GWMUL_STARTNEXTFFT_IF(safe1)); // t1
        gw().copy(*res.ZmY, *res.Y);
        gw().copy(*res.Z, *res.ZpY);
        gw().addsub(*res.Z, *res.Y, *res.Z, *res.Y, GWADD_DELAYNORM_IF(safe1));
    }

    void MontgomeryArithmetic::optimize(EdY& a)
    {
        //normalize(a);
    }

    void MontgomeryArithmetic::normalize(EdY& a)
    {
        std::vector<EdY*> tmp;
        tmp.push_back(&a);
        normalize(tmp.begin(), tmp.end());
    }

    template <typename Iter>
    void MontgomeryArithmetic::normalize(Iter begin, Iter end)
    {
        Iter it;
        Iter first = end;
        Iter last = end;
        for (it = begin; it != end; it++)
        {
            if (*it == nullptr)
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
        swap(*(*first)->ZpY, *(*first)->Z);
        Iter prev = first;
        for ((it = first)++; it != end; it++)
            if (*it != nullptr)
            {
                gw().mul(*(*prev)->ZpY, *(*it)->Z, *(*it)->ZpY, it != last ? GWMUL_STARTNEXTFFT : 0);
                prev = it;
            }
        try
        {
            (*last)->ZpY->inv();
        }
        catch (const ArithmeticException&)
        {
            swap(*(*first)->ZpY, *(*first)->Z);
            throw;
        }
        for ((it = last)++; it != first;)
        {
            it--;
            if (it != first)
            {
                for ((prev = it)--; *prev == nullptr; prev--);
                gw().mul(*(*it)->ZpY, *(*prev)->ZpY, *(*prev)->ZpY, GWMUL_STARTNEXTFFT);
                swap(*(*it)->ZpY, *(*prev)->ZpY);
                gw().mul(*(*it)->Z, *(*prev)->ZpY, *(*prev)->ZpY, GWMUL_STARTNEXTFFT);
            }
            if ((*it)->Y)
                gw().mul(*(*it)->ZpY, *(*it)->Y, *(*it)->Y, GWMUL_STARTNEXTFFT);
            (*it)->Z.reset();
            (*it)->ZpY.reset();
        }
    }
}
