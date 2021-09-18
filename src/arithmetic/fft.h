#pragma once

#include <memory>
#include <vector>

#include "arithmetic.h"

namespace arithmetic
{
    class Poly;

    class FFT
    {
    public:
        FFT(GWArithmetic& gw, int size);

        void transform(Poly& src, Poly& dst);
        void inv_transform(Poly& src, Poly& dst);

        GWArithmetic& gw() { return *_gw; }
        int size() const { return _size; }

    private:
        void transform(Poly& src, Poly& dst, int offset, int count, int depth);
        void inv_transform(Poly& src, Poly& dst, int offset, int count, int depth);

    private:
        GWArithmetic* _gw;
        int _size;
        std::vector<GWNum> _roots;
        std::vector<GWNum> _inv_roots;
    };

    class Poly : public std::vector<GWNum>
    {
    public:
        Poly(FFT& fft) : _fft(fft)
        {
            for (int i = 0; i < fft.size(); i++)
            {
                emplace_back(fft.gw());
                back() = 0;
            }
        }

    private:
        FFT& _fft;
    };
}
