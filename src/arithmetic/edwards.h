#pragma once

#include "group.h"

namespace arithmetic
{
    class EdPoint;

    class EdwardsArithmetic : public GroupArithmetic<EdPoint>
    {
        friend class EdPoint;

    public:
        static const int ED_PROJECTIVE = 0x100000;
        static const int EDADD_NEGATIVE = 0x200000;
        static const int EDDBL_FOR_EXT_NORM_ADD = 0x400000;

    public:
        EdwardsArithmetic(GWArithmetic& gw) : _gw(&gw) { }

        virtual void copy(const EdPoint& a, EdPoint& res) override;
        virtual void move(EdPoint&& a, EdPoint& res) override;
        virtual void init(EdPoint& res) override;
        virtual void init(const GWNum& X, const GWNum& Y, EdPoint& res);
        virtual void init(const GWNum& X, const GWNum& Y, const GWNum& Z, const GWNum& T, EdPoint& res);
        virtual void add(EdPoint& a, EdPoint& b, EdPoint& res) override;
        virtual void sub(EdPoint& a, EdPoint& b, EdPoint& res) override;
        virtual void add(EdPoint& a, EdPoint& b, EdPoint& res, int options);
        virtual void neg(EdPoint& a, EdPoint& res) override;
        virtual void dbl(EdPoint& a, EdPoint& res) override;
        virtual void dbl(EdPoint& a, EdPoint& res, int options);
        virtual void mul(EdPoint& a, Giant& b, EdPoint& res);
        virtual void mul(EdPoint& a, int W, std::vector<int16_t>& naf_w, EdPoint& res) override;

        virtual void normalize(EdPoint& a, int options);
        template <typename Iter>
        void normalize(Iter begin, Iter end, int options);
        EdPoint gen_curve(int seed, GWNum* ed_d);
        EdPoint from_small(int32_t xa, int32_t xb, int32_t ya, int32_t yb, GWNum* ed_d);
        GWNum jinvariant(GWNum& ed_d);
        bool on_curve(EdPoint& a, GWNum& ed_d);
        void d_ratio(EdPoint& a, GWNum& ed_d_a, GWNum& ed_d_b);

        GWArithmetic& gw() { return *_gw; }
        void set_gw(GWArithmetic& gw) { _gw = &gw; }

    private:
        GWArithmetic* _gw;
    };

    class EdPoint : public GroupElement<EdwardsArithmetic, EdPoint>
    {
        friend class EdwardsArithmetic;

    public:
        EdPoint(EdwardsArithmetic& arithmetic) : GroupElement<EdwardsArithmetic, EdPoint>(arithmetic)
        {
            arithmetic.init(*this);
        }
        EdPoint(EdwardsArithmetic& arithmetic, const GWNum& X, const GWNum& Y) : GroupElement<EdwardsArithmetic, EdPoint>(arithmetic)
        {
            arithmetic.init(X, Y, *this);
        }
        EdPoint(EdwardsArithmetic& arithmetic, const GWNum& X, const GWNum& Y, const GWNum& Z, const GWNum& T) : GroupElement<EdwardsArithmetic, EdPoint>(arithmetic)
        {
            arithmetic.init(X, Y, Z, T, *this);
        }
        ~EdPoint()
        {
        }
        EdPoint(const EdPoint& a) : GroupElement<EdwardsArithmetic, EdPoint>(a.arithmetic())
        {
            arithmetic().copy(a, *this);
        }
        EdPoint(EdPoint&& a) : GroupElement<EdwardsArithmetic, EdPoint>(a.arithmetic())
        {
            arithmetic().move(std::move(a), *this);
        }

        EdPoint& operator = (const EdPoint& a)
        {
            arithmetic().copy(a, *this);
            return *this;
        }
        EdPoint& operator = (EdPoint&& a)
        {
            arithmetic().move(std::move(a), *this);
            return *this;
        }

        friend bool operator == (EdPoint& a, EdPoint& b)
        {
            if (!a.Z && !b.Z)
                return (*a.X == *b.X) && (*a.Y == *b.Y);
            GWNum tmp = *a.X;
            if (b.Z)
                tmp *= *b.Z;
            GWNum tmp2 = *b.X;
            if (a.Z)
                tmp2 *= *a.Z;
            if (tmp != tmp2)
                return false;
            tmp = *a.X;
            if (b.Z)
                tmp *= *b.Z;
            tmp2 = *b.X;
            if (a.Z)
                tmp2 *= *a.Z;
            return (tmp == tmp2);
        }
        friend bool operator != (EdPoint& a, EdPoint& b)
        {
            return !(a == b);
        }
        EdPoint& normalize()
        {
            arithmetic().normalize(*this, 0);
            return *this;
        }
        EdPoint& extend()
        {
            if (!T)
                T.reset(new GWNum(arithmetic().gw()));
            *T = *X * (*Y);
            if (Z)
            {
                *X *= *Z;
                *Y *= *Z;
                *Z *= *Z;
            }
            return *this;
        }

    public:
        std::unique_ptr<GWNum> X;
        std::unique_ptr<GWNum> Y;
        std::unique_ptr<GWNum> Z;
        std::unique_ptr<GWNum> T;
    };

}
