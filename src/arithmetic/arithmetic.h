#pragma once

#include <unordered_set>

#include "field.h"
#include "giant.h"

namespace arithmetic
{
    class GWNum;

    class GWState
    {
    public:
        GWState() : giants(&handle)
        {
            gwinit(&handle);
        }
        ~GWState()
        {
            done();
        }

        void setup(int k, int b, int n, int c);
        void setup(const Giant& g);
        void done();

        gwhandle* gwdata() { return &handle; }
        Giant popg() { return Giant(giants); }

        int thread_count = 1;
        int next_fft_count = 0;
        double safety_margin = 0;
        int maxmulbyconst = 3;
        bool will_error_check = false;
        bool large_pages = false;

        gwhandle handle;
        GWGiantsArithmetic giants;
        Giant *N = nullptr;
        uint32_t fingerprint;
        std::string fft_description;
        int fft_length;
    };

    class CarefulGWArithmetic;
    class ReliableGWArithmetic;

    class GWArithmetic : public FieldArithmetic<GWNum>
    {
        friend class GWNum;

    protected:
        GWArithmetic(GWState& state, CarefulGWArithmetic *careful) : _state(state), _careful(careful) { }
    public:
        GWArithmetic(GWState& state);
        virtual ~GWArithmetic();

        virtual void alloc(GWNum& a) override;
        virtual void free(GWNum& a) override;
        virtual void copy(const GWNum& a, GWNum& res) override;
        virtual void move(GWNum&& a, GWNum& res) override;
        virtual void init(int32_t a, GWNum& res) override;
        virtual void init(const std::string& a, GWNum& res) override;
        virtual int cmp(const GWNum& a, const GWNum& b) override;
        virtual int cmp(const GWNum& a, int32_t b) override;
        virtual int cmp(const GWNum& a, const Giant& b);
        virtual void add(GWNum& a, GWNum& b, GWNum& res) override;
        virtual void add(GWNum& a, int32_t b, GWNum& res) override;
        virtual void add(GWNum& a, GWNum& b, GWNum& res, int options);
        virtual void sub(GWNum& a, GWNum& b, GWNum& res) override;
        virtual void sub(GWNum& a, int32_t b, GWNum& res) override;
        virtual void sub(GWNum& a, GWNum& b, GWNum& res, int options);
        virtual void neg(GWNum& a, GWNum& res) override;
        virtual void mul(GWNum& a, GWNum& b, GWNum& res) override;
        virtual void mul(GWNum& a, int32_t b, GWNum& res) override;
        virtual void mul(GWNum& a, GWNum& b, GWNum& res, int options);
        virtual void div(GWNum& a, GWNum& b, GWNum& res) override;
        virtual void div(GWNum& a, int32_t b, GWNum& res) override;
        virtual void gcd(GWNum& a, GWNum& b, GWNum& res) override;
        virtual void inv(GWNum& a, GWNum& n, GWNum& res) override;
        virtual void inv(GWNum& a, GWNum& res);
        virtual void mod(GWNum& a, GWNum& b, GWNum& res) override;

        virtual void addsub(GWNum& a, GWNum& b, GWNum& res1, GWNum& res2, int options);
        virtual void addmul(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options);
        virtual void submul(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options);
        virtual void muladd(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options);
        virtual void mulsub(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options);
        virtual void mulmuladd(GWNum& a, GWNum& b, GWNum& c, GWNum& d, GWNum& res, int options);
        virtual void mulmulsub(GWNum& a, GWNum& b, GWNum& c, GWNum& d, GWNum& res, int options);
        virtual void fft(GWNum& a, GWNum& res);
        virtual void unfft(GWNum& a, GWNum& res);

        void square(GWNum& a, GWNum& res, int options) { mul(a, a, res, options); }
        void setmulbyconst(int32_t a) { gwsetmulbyconst(gwdata(), a); }
        void setaddin(int32_t a) { gwsetaddin(gwdata(), a); }

        GWState& state() { return _state; }
        gwhandle* gwdata() { return &_state.handle; }
        Giant popg() { return Giant(_state.giants); }
        Giant& N() { GWASSERT(_state.N != nullptr); return *_state.N; }
        CarefulGWArithmetic& carefully() { return *_careful; }

    protected:
        GWState& _state;
        CarefulGWArithmetic* _careful = nullptr;
    };

    class CarefulGWArithmetic : public GWArithmetic
    {
    public:
        CarefulGWArithmetic(GWState& state) : GWArithmetic(state, this) { }
        virtual ~CarefulGWArithmetic() { }

        using GWArithmetic::add;
        using GWArithmetic::sub;
        using GWArithmetic::mul;
        virtual void add(GWNum& a, int32_t b, GWNum& res) override;
        virtual void sub(GWNum& a, int32_t b, GWNum& res) override;
        virtual void add(GWNum& a, GWNum& b, GWNum& res, int options) override;
        virtual void sub(GWNum& a, GWNum& b, GWNum& res, int options) override;
        virtual void mul(GWNum& a, GWNum& b, GWNum& res, int options) override;
        virtual void addsub(GWNum& a, GWNum& b, GWNum& res1, GWNum& res2, int options) override;
        virtual void addmul(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options) override;
        virtual void submul(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options) override;
        virtual void muladd(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options) override;
        virtual void mulsub(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options) override;
        virtual void mulmuladd(GWNum& a, GWNum& b, GWNum& c, GWNum& d, GWNum& res, int options) override;
        virtual void mulmulsub(GWNum& a, GWNum& b, GWNum& c, GWNum& d, GWNum& res, int options) override;
    };

    class ReliableGWArithmetic : public GWArithmetic
    {
    public:
        ReliableGWArithmetic(GWState& state) : GWArithmetic(state) { }
        virtual ~ReliableGWArithmetic() { }

        virtual void mul(GWNum& a, GWNum& b, GWNum& res, int options) override;
        virtual void addmul(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options) override;
        virtual void submul(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options) override;
        virtual void muladd(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options) override;
        virtual void mulsub(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options) override;
        virtual void mulmuladd(GWNum& a, GWNum& b, GWNum& c, GWNum& d, GWNum& res, int options) override;
        virtual void mulmulsub(GWNum& a, GWNum& b, GWNum& c, GWNum& d, GWNum& res, int options) override;

        void reset();
        void restart(int op = 0);

        int op() { return _op; }
        bool restart_flag() { return _restart_flag; }
        bool failure_flag() { return _failure_flag; }

    private:
        int _op = 0;
        std::unordered_set<int> _suspect_ops;
        bool _restart_flag = false;
        bool _failure_flag = false;
        double _max_roundoff = 0.4;
    };

    class GWNum : public FieldElement<GWArithmetic, GWNum>
    {
        friend class GWArithmetic;

    public:
        GWNum(GWArithmetic& arithmetic) : FieldElement<GWArithmetic, GWNum>(arithmetic)
        {
            arithmetic.alloc(*this);
        }
        ~GWNum()
        {
            if (_gwnum != nullptr)
                arithmetic().free(*this);
        }
        GWNum(const GWNum& a) : FieldElement<GWArithmetic, GWNum>(a.arithmetic())
        {
            arithmetic().copy(a, *this);
        }
        GWNum(GWNum&& a) : FieldElement<GWArithmetic, GWNum>(a.arithmetic())
        {
            arithmetic().move(std::move(a), *this);
        }

        GWNum& operator = (const GWNum& a)
        {
            arithmetic().copy(a, *this);
            return *this;
        }
        GWNum& operator = (GWNum&& a)
        {
            arithmetic().move(std::move(a), *this);
            return *this;
        }
        using FieldElement<GWArithmetic, GWNum>::operator=;

        virtual std::string to_string() const override;
        virtual GWNum& operator = (const Giant& a)
        {
            a.to_GWNum(*this);
            return *this;
        }

        friend gwnum operator *(const GWNum& a) { return a._gwnum; }

        GWNum& operator /= (GWNum&& a)
        {
            *this *= a.inv();
            return *this;
        }
        using FieldElement<GWArithmetic, GWNum>::operator/=;
        friend GWNum operator / (GWNum& a, GWNum& b)
        {
            GWNum t(b);
            t.inv();
            return a * std::move(t);
        }
        friend GWNum operator / (GWNum& a, GWNum&& b)
        {
            return a * std::move(b.inv());
        }
        friend GWNum operator / (GWNum&& a, GWNum& b)
        {
            GWNum res(std::move(a));
            res /= ((void*)&a == (void*)&b) ? res : b;
            return res;
        }
        friend GWNum operator / (GWNum&& a, GWNum&& b)
        {
            return std::move(a) * b.inv();
        }
        friend GWNum operator / (int32_t a, GWNum& b)
        {
            GWNum t(b);
            t.inv();
            return a * std::move(t);
        }
        friend GWNum operator / (int32_t a, GWNum&& b)
        {
            return a * std::move(b.inv());
        }
        friend Giant gcd(GWNum& a, Giant& b)
        {
            Giant res(b.arithmetic());
            res = a;
            res.gcd(b);
            return res;
        }
        GWNum& inv()
        {
            arithmetic().inv(*this, *this);
            return *this;
        }
        using FieldElement<GWArithmetic, GWNum>::inv;
        friend GWNum inv(GWNum& a)
        {
            GWNum t(a);
            t.inv();
            return t;
        }
        friend GWNum inv(GWNum&& a)
        {
            return a.inv();
        }

    private:
        gwnum _gwnum = nullptr;
    };
}
