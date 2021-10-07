
#include <vector>
#include <stdlib.h>
#include "gwnum.h"
#include "arithmetic.h"

namespace arithmetic
{
    void GWState::setup(int k, int b, int n, int c)
    {
        gwset_num_threads(gwdata(), thread_count);
        gwset_larger_fftlen_count(gwdata(), next_fft_count);
        gwset_safety_margin(gwdata(), safety_margin);
        gwset_maxmulbyconst(gwdata(), maxmulbyconst);
        if (will_error_check)
            gwset_will_error_check(gwdata());
        if (large_pages)
            gwset_use_large_pages(gwdata());
        if (gwsetup(gwdata(), k, b, n, c))
            throw std::exception();
        bit_length = (int)gwdata()->bit_length;
        giants._capacity = (bit_length >> 5) + 10;
        if (N != nullptr)
            delete N;
        N = new Giant(giants);
        *N = k*power(std::move(*N = b), n) + c;
        if (gwdata()->GENERAL_MOD)
            bit_length = N->bitlen();
        fingerprint = *N%3417905339UL;
        char buf[200];
        gwfft_description(gwdata(), buf);
        fft_description = buf;
        fft_length = gwfftlen(gwdata());
    }

    void GWState::setup(const Giant& g)
    {
        gwset_num_threads(gwdata(), thread_count);
        gwset_larger_fftlen_count(gwdata(), next_fft_count);
        gwset_safety_margin(gwdata(), safety_margin);
        gwset_maxmulbyconst(gwdata(), maxmulbyconst);
        if (will_error_check)
            gwset_will_error_check(gwdata());
        if (large_pages)
            gwset_use_large_pages(gwdata());
        if (gwsetup_general_mod(gwdata(), g.data(), g.size()))
            throw std::exception();
        bit_length = g.bitlen();
        giants._capacity = (bit_length >> 5) + 10;
        if (N != nullptr)
            delete N;
        N = new Giant(giants);
        *N = g;
        fingerprint = *N%3417905339UL;
        char buf[200];
        gwfft_description(gwdata(), buf);
        fft_description = buf;
        fft_length = gwfftlen(gwdata());
    }

    void GWState::setup(int bitlen)
    {
        gwset_num_threads(gwdata(), thread_count);
        gwset_larger_fftlen_count(gwdata(), next_fft_count);
        gwset_safety_margin(gwdata(), safety_margin);
        gwset_maxmulbyconst(gwdata(), maxmulbyconst);
        if (will_error_check)
            gwset_will_error_check(gwdata());
        if (large_pages)
            gwset_use_large_pages(gwdata());
        if (gwsetup_without_mod(gwdata(), bitlen))
            throw std::exception();
        bit_length = (int)gwdata()->bit_length;;
        giants._capacity = (bit_length >> 5) + 10;
        if (N != nullptr)
            delete N;
        fingerprint = 0;
        char buf[200];
        gwfft_description(gwdata(), buf);
        fft_description = buf;
        fft_length = gwfftlen(gwdata());
    }

    void GWState::done()
    {
        if (N != nullptr)
            delete N;
        N = nullptr;
        fft_description.clear();
        fft_length = 0;
        gwdone(&handle);
        gwinit(&handle);
    }

    GWArithmetic::GWArithmetic(GWState& state) : _state(state)
    {
        _careful = new CarefulGWArithmetic(state);
    }

    GWArithmetic::~GWArithmetic()
    {
        if (_careful != nullptr && (GWArithmetic*)_careful != this)
        {
            delete _careful;
            _careful = nullptr;
        }
    }

    void GWArithmetic::alloc(GWNum& a)
    {
        a._gwnum = gwalloc(gwdata());
    }

    void GWArithmetic::free(GWNum& a)
    {
        gwfree(gwdata(), a._gwnum);
        a._gwnum = nullptr;
    }

    void GWArithmetic::copy(const GWNum& a, GWNum& res)
    {
        if (*res == nullptr)
            alloc(res);
        gwcopy(gwdata(), *a, *res);
    }

    void GWArithmetic::move(GWNum&& a, GWNum& res)
    {
        if (res._gwnum != nullptr)
            free(res);
        res._gwnum = a._gwnum;
        a._gwnum = nullptr;
    }

    void GWArithmetic::init(int32_t a, GWNum& res)
    {
        dbltogw(gwdata(), a, *res);
    }

    void GWArithmetic::init(const std::string& a, GWNum& res)
    {
        (popg() = a).to_GWNum(res);
    }

    int GWArithmetic::cmp(const GWNum& a, const GWNum& b)
    {
        return _state.giants.cmp(popg() = a, popg() = b);
    }

    int GWArithmetic::cmp(const GWNum& a, int32_t b)
    {
        if (b >= 0)
            return _state.giants.cmp(popg() = a, b);
        else
            return _state.giants.cmp(popg() = a, N() + b);
    }

    int GWArithmetic::cmp(const GWNum& a, const Giant& b)
    {
        return _state.giants.cmp(popg() = a, b);
    }

    void GWArithmetic::add(GWNum& a, GWNum& b, GWNum& res)
    {
        add(a, b, res, GWADD_DELAY_NORMALIZE);
    }

    void GWArithmetic::add(GWNum& a, int32_t b, GWNum& res)
    {
        if (*res != *a)
            copy(a, res);
        gwsmalladd(gwdata(), b, *res);
    }

    void GWArithmetic::add(GWNum& a, GWNum& b, GWNum& res, int options)
    {
        gwadd3o(gwdata(), *a, *b, *res, options);
    }

    void GWArithmetic::sub(GWNum& a, GWNum& b, GWNum& res)
    {
        sub(a, b, res, GWADD_DELAY_NORMALIZE);
    }

    void GWArithmetic::sub(GWNum& a, int32_t b, GWNum& res)
    {
        if (*res != *a)
            copy(a, res);
        gwsmalladd(gwdata(), -b, *res);
    }

    void GWArithmetic::sub(GWNum& a, GWNum& b, GWNum& res, int options)
    {
        gwsub3o(gwdata(), *a, *b, *res, options);
    }

    void GWArithmetic::neg(GWNum& a, GWNum& res)
    {
        if (*res != *a)
            copy(a, res);
        gwsmallmul(gwdata(), -1, *res);
    }

    void GWArithmetic::mul(GWNum& a, GWNum& b, GWNum& res)
    {
        mul(a, b, res, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
    }

    void GWArithmetic::mul(GWNum& a, int32_t b, GWNum& res)
    {
        if (*res != *a)
            copy(a, res);
        gwsmallmul(gwdata(), b, *res);
    }

    void GWArithmetic::mul(GWNum& a, GWNum& b, GWNum& res, int options)
    {
        gwmul3(gwdata(), *a, *b, *res, options);
    }

    void GWArithmetic::addsub(GWNum& a, GWNum& b, GWNum& res1, GWNum& res2, int options)
    {
        gwaddsub4o(gwdata(), *a, *b, *res1, *res2, options);
    }

    void GWArithmetic::addmul(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options)
    {
        gwaddmul4(gwdata(), *a, *b, *c, *res, options);
    }

    void GWArithmetic::submul(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options)
    {
        gwsubmul4(gwdata(), *a, *b, *c, *res, options);
    }

    void GWArithmetic::muladd(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options)
    {
        gwmuladd4(gwdata(), *a, *b, *c, *res, options);
    }

    void GWArithmetic::mulsub(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options)
    {
        gwmulsub4(gwdata(), *a, *b, *c, *res, options);
    }

    void GWArithmetic::mulmuladd(GWNum& a, GWNum& b, GWNum& c, GWNum& d, GWNum& res, int options)
    {
        gwmulmuladd5(gwdata(), *a, *b, *c, *d, *res, options);
    }

    void GWArithmetic::mulmulsub(GWNum& a, GWNum& b, GWNum& c, GWNum& d, GWNum& res, int options)
    {
        gwmulmulsub5(gwdata(), *a, *b, *c, *d, *res, options);
    }

    void GWArithmetic::fft(GWNum& a, GWNum& res)
    {
        gwfft(gwdata(), *a, *res);
    }

    void GWArithmetic::unfft(GWNum& a, GWNum& res)
    {
        gwunfft(gwdata(), *a, *res);
    }

    void GWArithmetic::div(GWNum& a, GWNum& b, GWNum& res)
    {
        GWNum inv_b = b;
        inv_b.inv();
        mul(a, inv_b, res);
    }

    void GWArithmetic::div(GWNum& a, int32_t b, GWNum& res)
    {
        GWNum inv_b(a.arithmetic());
        (popg() = b).inv(N()).to_GWNum(inv_b);
        mul(a, inv_b, res);
    }

    void GWArithmetic::gcd(GWNum& a, GWNum& b, GWNum& res)
    {
        (popg() = a).gcd(popg() = b).to_GWNum(res);
    }

    void GWArithmetic::inv(GWNum& a, GWNum& n, GWNum& res)
    {
        (popg() = a).inv(popg() = n).to_GWNum(res);
    }

    void GWArithmetic::inv(GWNum& a, GWNum& res)
    {
        (popg() = a).inv(N()).to_GWNum(res);
    }

    void GWArithmetic::mod(GWNum& a, GWNum& b, GWNum& res)
    {
        ((popg() = a)%(popg() = b)).to_GWNum(res);
    }

    std::string GWNum::to_string() const
    {
        return (arithmetic().popg() = *this).to_string();
    }

    void CarefulGWArithmetic::add(GWNum& a, int32_t b, GWNum& res)
    {
        GWNum tmp(*this);
        tmp = b;
        unfft(a, a);
        gwadd3o(gwdata(), *a, *tmp, *res, GWADD_FORCE_NORMALIZE);
    }

    void CarefulGWArithmetic::sub(GWNum& a, int32_t b, GWNum& res)
    {
        GWNum tmp(*this);
        tmp = b;
        unfft(a, a);
        gwsub3o(gwdata(), *a, *tmp, *res, GWADD_FORCE_NORMALIZE);
    }

    void CarefulGWArithmetic::add(GWNum& a, GWNum& b, GWNum& res, int options)
    {
        unfft(a, a);
        unfft(b, b);
        gwadd3o(gwdata(), *a, *b, *res, GWADD_FORCE_NORMALIZE);
    }

    void CarefulGWArithmetic::sub(GWNum& a, GWNum& b, GWNum& res, int options)
    {
        unfft(a, a);
        unfft(b, b);
        gwsub3o(gwdata(), *a, *b, *res, GWADD_FORCE_NORMALIZE);
    }

    void CarefulGWArithmetic::mul(GWNum& a, GWNum& b, GWNum& res, int options)
    {
        unfft(a, a);
        unfft(b, b);
        gwmul3_carefully(gwdata(), *a, *b, *res, options);
    }

    void CarefulGWArithmetic::addsub(GWNum& a, GWNum& b, GWNum& res1, GWNum& res2, int options)
    {
        unfft(a, a);
        unfft(b, b);
        gwaddsub4o(gwdata(), *a, *b, *res1, *res2, GWADD_FORCE_NORMALIZE);
    }

    void CarefulGWArithmetic::addmul(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options)
    {
        if (*c == *res)
        {
            GWNum tmp = c;
            add(a, b, res);
            mul(res, tmp, res, options);
        }
        else
        {
            add(a, b, res);
            mul(res, c, res, options);
        }
    }

    void CarefulGWArithmetic::submul(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options)
    {
        if (*c == *res)
        {
            GWNum tmp = c;
            sub(a, b, res);
            mul(res, tmp, res, options);
        }
        else
        {
            sub(a, b, res);
            mul(res, c, res, options);
        }
    }

    void CarefulGWArithmetic::muladd(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options)
    {
        if (*c == *res || (options & GWMUL_MULBYCONST))
        {
            GWNum tmp = c;
            if (options & GWMUL_MULBYCONST)
                gwsmallmul(gwdata(), gwdata()->mulbyconst, *tmp);
            mul(a, b, res, options);
            add(res, tmp, res);
        }
        else
        {
            mul(a, b, res, options);
            add(res, c, res);
        }
    }

    void CarefulGWArithmetic::mulsub(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options)
    {
        if (*c == *res || (options & GWMUL_MULBYCONST))
        {
            GWNum tmp = c;
            if (options & GWMUL_MULBYCONST)
                gwsmallmul(gwdata(), gwdata()->mulbyconst, *tmp);
            mul(a, b, res, options);
            sub(res, tmp, res);
        }
        else
        {
            mul(a, b, res, options);
            sub(res, c, res);
        }
    }

    void CarefulGWArithmetic::mulmuladd(GWNum& a, GWNum& b, GWNum& c, GWNum& d, GWNum& res, int options)
    {
        GWNum tmp(*this);
        mul(c, d, tmp, options & (~GWMUL_ADDINCONST));
        mul(a, b, res, options);
        add(res, tmp, res);
    }

    void CarefulGWArithmetic::mulmulsub(GWNum& a, GWNum& b, GWNum& c, GWNum& d, GWNum& res, int options)
    {
        GWNum tmp(*this);
        mul(c, d, tmp, options & (~GWMUL_ADDINCONST));
        mul(a, b, res, options);
        sub(res, tmp, res);
    }

    void ReliableGWArithmetic::reset()
    {
        _op = 0;
        _suspect_ops.clear();
        _restart_flag = false;
        _failure_flag = false;
    }

    void ReliableGWArithmetic::restart(int op)
    {
        _op = op;
        _restart_flag = false;
    }

    void ReliableGWArithmetic::mul(GWNum& a, GWNum& b, GWNum& res, int options)
    {
        bool suspect = _suspect_ops.count(_op) != 0;
        if (!suspect)
        {
            // Try a normal operation.
            gwerror_checking(gwdata(), true);
            gwmul3(gwdata(), *a, *b, *res, options);
            gwerror_checking(gwdata(), false);
            if (gw_get_maxerr(gwdata()) > _max_roundoff)
            {
                // If roundoff exceeds maximum, mark the operation as suspect.
                gw_clear_maxerr(gwdata());
                _suspect_ops.insert(_op);
                // If either source is the same as destination, we can't redo the operation.
                // Need to restart from a checkpoint.
                if (*res == *a || *res == *b)
                    _restart_flag = true;
                else
                    suspect = true;
            }
        }
        if (suspect)
        {
            // If operation is marked as suspect, make a copy of the source arguments.
            GWNum s1 = a;
            GWNum s2 = b;
            // Retry it with original arguments in case it's a hardware problem.
            gwerror_checking(gwdata(), true);
            gwmul3(gwdata(), *a, *b, *res, options);
            gwerror_checking(gwdata(), false);
            if (gw_get_maxerr(gwdata()) > _max_roundoff)
            {
                // Not a hardware problem, retry carefully.
                gw_clear_maxerr(gwdata());
                gwerror_checking(gwdata(), true);
                gwmul3_carefully(gwdata(), *s1, (*a == *b) ? *s1 : *s2, *res, options & (~GWMUL_PRESERVE_S1) & (~GWMUL_PRESERVE_S2));
                gwerror_checking(gwdata(), false);
                if (gw_get_maxerr(gwdata()) > _max_roundoff)
                {
                    // All is lost, need larger FFT.
                    gw_clear_maxerr(gwdata());
                    _failure_flag = true;
                }
            }
            else // It's a hardware problem after all.
                _suspect_ops.erase(_op);
        }
        // Increase operation counter.
        _op++;
    }

    void ReliableGWArithmetic::addmul(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options)
    {
        bool suspect = _suspect_ops.count(_op) != 0;
        if (!suspect)
        {
            gwerror_checking(gwdata(), true);
            gwaddmul4(gwdata(), *a, *b, *c, *res, options);
            gwerror_checking(gwdata(), false);
            if (gw_get_maxerr(gwdata()) > _max_roundoff)
            {
                gw_clear_maxerr(gwdata());
                _suspect_ops.insert(_op);
                if (*res == *a || *res == *b || *res == *c)
                    _restart_flag = true;
                else
                    suspect = true;
            }
        }
        if (suspect)
        {
            GWNum s1 = a;
            GWNum s2 = b;
            GWNum s3 = c;
            gwerror_checking(gwdata(), true);
            gwaddmul4(gwdata(), *a, *b, *c, *res, options);
            gwerror_checking(gwdata(), false);
            if (gw_get_maxerr(gwdata()) > _max_roundoff)
            {
                gw_clear_maxerr(gwdata());
                gwerror_checking(gwdata(), true);
                gwunfft(gwdata(), *s1, *s1);
                gwunfft(gwdata(), *s2, *s2);
                gwadd3o(gwdata(), *s1, *s2, *s1, GWADD_FORCE_NORMALIZE);
                gwmul3_carefully(gwdata(), *s1, *s3, *res, options & (~GWMUL_PRESERVE_S1) & (~GWMUL_PRESERVE_S2));
                gwerror_checking(gwdata(), false);
                if (gw_get_maxerr(gwdata()) > _max_roundoff)
                {
                    gw_clear_maxerr(gwdata());
                    _failure_flag = true;
                }
            }
            else
                _suspect_ops.erase(_op);
        }
        _op++;
    }

    void ReliableGWArithmetic::submul(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options)
    {
        bool suspect = _suspect_ops.count(_op) != 0;
        if (!suspect)
        {
            gwerror_checking(gwdata(), true);
            gwsubmul4(gwdata(), *a, *b, *c, *res, options);
            gwerror_checking(gwdata(), false);
            if (gw_get_maxerr(gwdata()) > _max_roundoff)
            {
                gw_clear_maxerr(gwdata());
                _suspect_ops.insert(_op);
                if (*res == *a || *res == *b || *res == *c)
                    _restart_flag = true;
                else
                    suspect = true;
            }
        }
        if (suspect)
        {
            GWNum s1 = a;
            GWNum s2 = b;
            GWNum s3 = c;
            gwerror_checking(gwdata(), true);
            gwsubmul4(gwdata(), *a, *b, *c, *res, options);
            gwerror_checking(gwdata(), false);
            if (gw_get_maxerr(gwdata()) > _max_roundoff)
            {
                gw_clear_maxerr(gwdata());
                gwerror_checking(gwdata(), true);
                gwunfft(gwdata(), *s1, *s1);
                gwunfft(gwdata(), *s2, *s2);
                gwsub3o(gwdata(), *s1, *s2, *s1, GWADD_FORCE_NORMALIZE);
                gwmul3_carefully(gwdata(), *s1, *s3, *res, options & (~GWMUL_PRESERVE_S1) & (~GWMUL_PRESERVE_S2));
                gwerror_checking(gwdata(), false);
                if (gw_get_maxerr(gwdata()) > _max_roundoff)
                {
                    gw_clear_maxerr(gwdata());
                    _failure_flag = true;
                }
            }
            else
                _suspect_ops.erase(_op);
        }
        _op++;
    }

    void ReliableGWArithmetic::muladd(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options)
    {
        bool suspect = _suspect_ops.count(_op) != 0;
        if (!suspect)
        {
            gwerror_checking(gwdata(), true);
            gwmuladd4(gwdata(), *a, *b, *c, *res, options);
            gwerror_checking(gwdata(), false);
            if (gw_get_maxerr(gwdata()) > _max_roundoff)
            {
                gw_clear_maxerr(gwdata());
                _suspect_ops.insert(_op);
                if (*res == *a || *res == *b || *res == *c)
                    _restart_flag = true;
                else
                    suspect = true;
            }
        }
        if (suspect)
        {
            GWNum s1 = a;
            GWNum s2 = b;
            GWNum s3 = c;
            gwerror_checking(gwdata(), true);
            gwmuladd4(gwdata(), *a, *b, *c, *res, options);
            gwerror_checking(gwdata(), false);
            if (gw_get_maxerr(gwdata()) > _max_roundoff)
            {
                gw_clear_maxerr(gwdata());
                gwerror_checking(gwdata(), true);
                gwmul3_carefully(gwdata(), *s1, (*a == *b) ? *s1 : *s2, *res, options & (~GWMUL_PRESERVE_S1) & (~GWMUL_PRESERVE_S2) & (~GWMUL_STARTNEXTFFT));
                gwunfft(gwdata(), *s3, *s3);
                if (options & GWMUL_MULBYCONST)
                    gwsmallmul(gwdata(), gwdata()->mulbyconst, *s3);
                gwadd3o(gwdata(), *res, *s3, *res, GWADD_FORCE_NORMALIZE);
                gwerror_checking(gwdata(), false);
                if (gw_get_maxerr(gwdata()) > _max_roundoff)
                {
                    gw_clear_maxerr(gwdata());
                    _failure_flag = true;
                }
            }
            else
                _suspect_ops.erase(_op);
        }
        _op++;
    }

    void ReliableGWArithmetic::mulsub(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options)
    {
        bool suspect = _suspect_ops.count(_op) != 0;
        if (!suspect)
        {
            gwerror_checking(gwdata(), true);
            gwmulsub4(gwdata(), *a, *b, *c, *res, options);
            gwerror_checking(gwdata(), false);
            if (gw_get_maxerr(gwdata()) > _max_roundoff)
            {
                gw_clear_maxerr(gwdata());
                _suspect_ops.insert(_op);
                if (*res == *a || *res == *b || *res == *c)
                    _restart_flag = true;
                else
                    suspect = true;
            }
        }
        if (suspect)
        {
            GWNum s1 = a;
            GWNum s2 = b;
            GWNum s3 = c;
            gwerror_checking(gwdata(), true);
            gwmulsub4(gwdata(), *a, *b, *c, *res, options);
            gwerror_checking(gwdata(), false);
            if (gw_get_maxerr(gwdata()) > _max_roundoff)
            {
                gw_clear_maxerr(gwdata());
                gwerror_checking(gwdata(), true);
                gwmul3_carefully(gwdata(), *s1, (*a == *b) ? *s1 : *s2, *res, options & (~GWMUL_PRESERVE_S1) & (~GWMUL_PRESERVE_S2) & (~GWMUL_STARTNEXTFFT));
                gwunfft(gwdata(), *s3, *s3);
                if (options & GWMUL_MULBYCONST)
                    gwsmallmul(gwdata(), gwdata()->mulbyconst, *s3);
                gwsub3o(gwdata(), *res, *s3, *res, GWADD_FORCE_NORMALIZE);
                gwerror_checking(gwdata(), false);
                if (gw_get_maxerr(gwdata()) > _max_roundoff)
                {
                    gw_clear_maxerr(gwdata());
                    _failure_flag = true;
                }
            }
            else
                _suspect_ops.erase(_op);
        }
        _op++;
    }

    void ReliableGWArithmetic::mulmuladd(GWNum& a, GWNum& b, GWNum& c, GWNum& d, GWNum& res, int options)
    {
        bool suspect = _suspect_ops.count(_op) != 0;
        if (!suspect)
        {
            gwerror_checking(gwdata(), true);
            gwmulmuladd5(gwdata(), *a, *b, *c, *d, *res, options);
            gwerror_checking(gwdata(), false);
            if (gw_get_maxerr(gwdata()) > _max_roundoff)
            {
                gw_clear_maxerr(gwdata());
                _suspect_ops.insert(_op);
                if (*res == *a || *res == *b || *res == *c || *res == *d)
                    _restart_flag = true;
                else
                    suspect = true;
            }
        }
        if (suspect)
        {
            GWNum s1 = a;
            GWNum s2 = b;
            GWNum s3 = c;
            GWNum s4 = d;
            gwerror_checking(gwdata(), true);
            gwmulmuladd5(gwdata(), *a, *b, *c, *d, *res, options);
            gwerror_checking(gwdata(), false);
            if (gw_get_maxerr(gwdata()) > _max_roundoff)
            {
                gw_clear_maxerr(gwdata());
                gwerror_checking(gwdata(), true);
                gwmul3_carefully(gwdata(), *s1, (*a == *b) ? *s1 : *s2, *res, options & (~GWMUL_PRESERVE_S1) & (~GWMUL_PRESERVE_S2) & (~GWMUL_STARTNEXTFFT));
                gwmul3_carefully(gwdata(), *s3, (*c == *d) ? *s3 : *s4, *s3, options & (~GWMUL_PRESERVE_S1) & (~GWMUL_PRESERVE_S2) & (~GWMUL_STARTNEXTFFT) & (~GWMUL_ADDINCONST));
                gwadd3o(gwdata(), *res, *s3, *res, GWADD_FORCE_NORMALIZE);
                gwerror_checking(gwdata(), false);
                if (gw_get_maxerr(gwdata()) > _max_roundoff)
                {
                    gw_clear_maxerr(gwdata());
                    _failure_flag = true;
                }
            }
            else
                _suspect_ops.erase(_op);
        }
        _op++;
    }

    void ReliableGWArithmetic::mulmulsub(GWNum& a, GWNum& b, GWNum& c, GWNum& d, GWNum& res, int options)
    {
        bool suspect = _suspect_ops.count(_op) != 0;
        if (!suspect)
        {
            gwerror_checking(gwdata(), true);
            gwmulmulsub5(gwdata(), *a, *b, *c, *d, *res, options);
            gwerror_checking(gwdata(), false);
            if (gw_get_maxerr(gwdata()) > _max_roundoff)
            {
                gw_clear_maxerr(gwdata());
                _suspect_ops.insert(_op);
                if (*res == *a || *res == *b || *res == *c || *res == *d)
                    _restart_flag = true;
                else
                    suspect = true;
            }
        }
        if (suspect)
        {
            GWNum s1 = a;
            GWNum s2 = b;
            GWNum s3 = c;
            GWNum s4 = d;
            gwerror_checking(gwdata(), true);
            gwmulmulsub5(gwdata(), *a, *b, *c, *d, *res, options);
            gwerror_checking(gwdata(), false);
            if (gw_get_maxerr(gwdata()) > _max_roundoff)
            {
                gw_clear_maxerr(gwdata());
                gwerror_checking(gwdata(), true);
                gwmul3_carefully(gwdata(), *s1, (*a == *b) ? *s1 : *s2, *res, options & (~GWMUL_PRESERVE_S1) & (~GWMUL_PRESERVE_S2) & (~GWMUL_STARTNEXTFFT));
                gwmul3_carefully(gwdata(), *s3, (*c == *d) ? *s3 : *s4, *s3, options & (~GWMUL_PRESERVE_S1) & (~GWMUL_PRESERVE_S2) & (~GWMUL_STARTNEXTFFT) & (~GWMUL_ADDINCONST));
                gwsub3o(gwdata(), *res, *s3, *res, GWADD_FORCE_NORMALIZE);
                gwerror_checking(gwdata(), false);
                if (gw_get_maxerr(gwdata()) > _max_roundoff)
                {
                    gw_clear_maxerr(gwdata());
                    _failure_flag = true;
                }
            }
            else
                _suspect_ops.erase(_op);
        }
        _op++;
    }
}
