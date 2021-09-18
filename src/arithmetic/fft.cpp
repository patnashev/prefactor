#include <stdlib.h>
#include "gwnum.h"
#include "fft.h"
#include "exception.h"

namespace arithmetic
{
    FFT::FFT(GWArithmetic& gw, int size) : _gw(&gw), _size(size)
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
                GWASSERT(i == j || _roots[i] != _roots[j]);
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

    void FFT::transform(Poly& src, Poly& dst)
    {
        transform(src, dst, 0, _size, 0);
    }

    void FFT::transform(Poly& src, Poly& dst, int offset, int count, int root)
    {
        if (count == 1)
            return;
        int i, m;
        m = count/2;
        for (i = 0; i < m; i++)
        {
            _gw->mul(src[offset + m + i], _roots[root], dst[offset + m + i], GWMUL_FFT_S2);
            _gw->addsub(src[offset + i], dst[offset + m + i], dst[offset + i], dst[offset + m + i], GWADD_FORCE_NORMALIZE);
        }
        transform(src, dst, offset, count/2, root*2);
        transform(src, dst, offset + m, count/2, root*2 + 1);
    }

    void FFT::inv_transform(Poly& src, Poly& dst)
    {
        inv_transform(src, dst, 0, _size, 0);
        for (auto it = dst.begin(); it != dst.end(); it++)
            _gw->mul(_inv_roots[0], *it, *it, GWMUL_FFT_S1);
    }

    void FFT::inv_transform(Poly& src, Poly& dst, int offset, int count, int root)
    {
        if (count == 1)
            return;
        int i, m;
        m = count/2;
        inv_transform(src, dst, offset, count/2, root*2);
        inv_transform(src, dst, offset + m, count/2, root*2 + 1);
        for (i = 0; i < m; i++)
        {
            _gw->addsub(src[offset + i], src[offset + m + i], dst[offset + i], dst[offset + m + i], GWADD_FORCE_NORMALIZE);
            if (root != 0)
               _gw->mul(dst[offset + m + i], _inv_roots[root], dst[offset + m + i], GWMUL_FFT_S2);
        }
    }
}
