#pragma once

#include <string>
#include <exception>
#include <unordered_set>

namespace arithmetic
{
    template<class Element>
    class FieldArithmetic
    {
    public:
        virtual void alloc(Element& a) = 0;
        virtual void free(Element& a) = 0;
        virtual void copy(const Element& a, Element& res) = 0;
        virtual void move(Element&& a, Element& res) = 0;
        virtual void init(int32_t a, Element& res) = 0;
        virtual void init(const std::string& a, Element& res) = 0;
        virtual int cmp(const Element& a, const Element& b) = 0;
        virtual int cmp(const Element& a, int32_t b) = 0;
        virtual void add(Element& a, Element& b, Element& res) = 0;
        virtual void add(Element& a, int32_t b, Element& res) = 0;
        virtual void sub(Element& a, Element& b, Element& res) = 0;
        virtual void sub(Element& a, int32_t b, Element& res) = 0;
        virtual void neg(Element& a, Element& res) = 0;
        virtual void mul(Element& a, Element& b, Element& res) = 0;
        virtual void mul(Element& a, int32_t b, Element& res) = 0;
        virtual void div(Element& a, Element& b, Element& res) = 0;
        virtual void div(Element& a, int32_t b, Element& res) = 0;
        virtual void gcd(Element& a, Element& b, Element& res) = 0;
        virtual void inv(Element& a, Element& n, Element& res) = 0;
        virtual void mod(Element& a, Element& b, Element& res) = 0;
    };

    template<class Arithmetic, class Element>
    class FieldElement
    {
    public:
        FieldElement(Arithmetic& arithmetic) : _arithmetic(arithmetic) { }

        Arithmetic& arithmetic() const { return _arithmetic; }

        virtual std::string to_string() const = 0;

        Element& operator = (int32_t a)
        {
            arithmetic().init(a, (Element&)*this);
            return (Element&)*this;
        }

        Element& operator = (const std::string& a)
        {
            arithmetic().init(a, (Element&)*this);
            return (Element&)*this;
        }

        friend bool operator == (const Element& a, const Element& b)
        {
            return a.arithmetic().cmp(a, b) == 0;
        }

        friend bool operator != (const Element& a, const Element& b)
        {
            return a.arithmetic().cmp(a, b) != 0;
        }

        friend bool operator > (const Element& a, const Element& b)
        {
            return a.arithmetic().cmp(a, b) > 0;
        }

        friend bool operator < (const Element& a, const Element& b)
        {
            return a.arithmetic().cmp(a, b) < 0;
        }

        friend bool operator >= (const Element& a, const Element& b)
        {
            return a.arithmetic().cmp(a, b) >= 0;
        }

        friend bool operator <= (const Element& a, const Element& b)
        {
            return a.arithmetic().cmp(a, b) <= 0;
        }

        friend bool operator == (const Element& a, int b)
        {
            return a.arithmetic().cmp(a, b) == 0;
        }

        friend bool operator != (const Element& a, int b)
        {
            return a.arithmetic().cmp(a, b) != 0;
        }

        friend bool operator > (const Element& a, int b)
        {
            return a.arithmetic().cmp(a, b) > 0;
        }

        friend bool operator < (const Element& a, int b)
        {
            return a.arithmetic().cmp(a, b) < 0;
        }

        friend bool operator >= (const Element& a, int b)
        {
            return a.arithmetic().cmp(a, b) >= 0;
        }

        friend bool operator <= (const Element& a, int b)
        {
            return a.arithmetic().cmp(a, b) <= 0;
        }

        friend bool operator == (int a, const Element& b)
        {
            return b.arithmetic().cmp(b, a) == 0;
        }

        friend bool operator != (int a, const Element& b)
        {
            return b.arithmetic().cmp(b, a) != 0;
        }

        friend bool operator > (int a, const Element& b)
        {
            return b.arithmetic().cmp(b, a) < 0;
        }

        friend bool operator < (int a, const Element& b)
        {
            return b.arithmetic().cmp(b, a) > 0;
        }

        friend bool operator >= (int a, const Element& b)
        {
            return b.arithmetic().cmp(b, a) <= 0;
        }

        friend bool operator <= (int a, const Element& b)
        {
            return b.arithmetic().cmp(b, a) >= 0;
        }

        Element& add(Element& a)
        {
            arithmetic().add((Element&)*this, a, (Element&)*this);
            return (Element&)*this;
        }

        Element& operator += (Element& a)
        {
            arithmetic().add((Element&)*this, a, (Element&)*this);
            return (Element&)*this;
        }

        friend Element operator + (Element& a, Element& b)
        {
            Element res(a.arithmetic());
            res.arithmetic().add(a, b, res);
            return res;
        }

        friend Element operator + (Element&& a, Element& b)
        {
            Element res(std::move(a));
            res += ((void*)&a == (void*)&b) ? res : b;
            return res;
        }

        friend Element operator + (Element& a, Element&& b)
        {
            Element res(std::move(b));
            res += ((void*)&a == (void*)&b) ? res : a;
            return res;
        }

        friend Element operator + (Element&& a, Element&& b)
        {
            return std::move(a) + b;
        }

        Element& operator += (int32_t a)
        {
            arithmetic().add((Element&)*this, a, (Element&)*this);
            return (Element&)*this;
        }

        friend Element operator + (Element& a, int32_t b)
        {
            Element res(a.arithmetic());
            res.arithmetic().add(a, b, res);
            return res;
        }

        friend Element operator + (int32_t a, Element& b)
        {
            Element res(b.arithmetic());
            res.arithmetic().add(b, a, res);
            return res;
        }

        friend Element operator + (Element&& a, int32_t b)
        {
            Element res(std::move(a));
            res += b;
            return res;
        }

        friend Element operator + (int32_t a, Element&& b)
        {
            Element res(std::move(b));
            res += a;
            return res;
        }

        Element& operator -= (Element& a)
        {
            arithmetic().sub((Element&)*this, a, (Element&)*this);
            return (Element&)*this;
        }

        friend Element operator - (Element& a, Element& b)
        {
            Element res(a.arithmetic());
            res.arithmetic().sub(a, b, res);
            return res;
        }

        friend Element operator - (Element&& a, Element& b)
        {
            Element res(std::move(a));
            res -= ((void*)&a == (void*)&b) ? res : b;
            return res;
        }

        friend Element operator - (Element& a, Element&& b)
        {
            Element res(std::move(b));
            res -= ((void*)&a == (void*)&b) ? res : a;
            res.arithmetic().neg(res, res);
            return res;
        }

        friend Element operator - (Element&& a, Element&& b)
        {
            return std::move(a) - b;
        }

        Element& operator -= (int32_t a)
        {
            arithmetic().sub((Element&)*this, a, (Element&)*this);
            return (Element&)*this;
        }

        friend Element operator - (Element& a, int32_t b)
        {
            Element res(a.arithmetic());
            res.arithmetic().sub(a, b, res);
            return res;
        }

        friend Element operator - (int32_t a, Element& b)
        {
            Element res(b.arithmetic());
            res.arithmetic().init(a, res);
            res.arithmetic().sub(res, b, res);
            return res;
        }

        friend Element operator - (Element&& a, int32_t b)
        {
            Element res(std::move(a));
            res -= b;
            return res;
        }

        friend Element operator - (int32_t a, Element&& b)
        {
            Element res(std::move(b));
            res -= a;
            res.arithmetic().neg(res, res);
            return res;
        }

        friend Element operator - (Element& a)
        {
            Element res(a.arithmetic());
            res.arithmetic().neg(a, res);
            return res;
        }

        friend Element operator - (Element&& a)
        {
            Element res(std::move(a));
            res.arithmetic().neg(res, res);
            return res;
        }

        Element& square()
        {
            arithmetic().mul((Element&)*this, (Element&)*this, (Element&)*this);
            return (Element&)*this;
        }

        friend Element square(Element& a)
        {
            Element res(a.arithmetic());
            res.arithmetic().mul(a, a, res);
            return res;
        }

        friend Element square(Element&& a)
        {
            Element res(std::move(a));
            res.square();
            return res;
        }

        Element& operator *= (Element& a)
        {
            arithmetic().mul((Element&)*this, a, (Element&)*this);
            return (Element&)*this;
        }

        Element& operator *= (Element&& a)
        {
            *this *= a;
            return (Element&)*this;
        }

        friend Element operator * (Element& a, Element& b)
        {
            Element res(a.arithmetic());
            res.arithmetic().mul(a, b, res);
            return res;
        }

        friend Element operator * (Element&& a, Element& b)
        {
            Element res(std::move(a));
            res *= ((void*)&a == (void*)&b) ? res : b;
            return res;
        }

        friend Element operator * (Element& a, Element&& b)
        {
            Element res(std::move(b));
            res *= ((void*)&a == (void*)&b) ? res : a;
            return res;
        }

        friend Element operator * (Element&& a, Element&& b)
        {
            return std::move(a) * b;
        }

        Element& operator *= (int32_t a)
        {
            arithmetic().mul((Element&)*this, a, (Element&)*this);
            return (Element&)*this;
        }

        friend Element operator * (Element& a, int32_t b)
        {
            Element res(a.arithmetic());
            res.arithmetic().mul(a, b, res);
            return res;
        }

        friend Element operator * (int32_t a, Element& b)
        {
            Element res(b.arithmetic());
            res.arithmetic().mul(b, a, res);
            return res;
        }

        friend Element operator * (Element&& a, int32_t b)
        {
            Element res(std::move(a));
            res *= b;
            return res;
        }

        friend Element operator * (int32_t a, Element&& b)
        {
            Element res(std::move(b));
            res *= a;
            return res;
        }

        Element& operator /= (Element& a)
        {
            arithmetic().div((Element&)*this, a, (Element&)*this);
            return (Element&)*this;
        }

        Element& operator /= (int32_t a)
        {
            arithmetic().div((Element&)*this, a, (Element&)*this);
            return (Element&)*this;
        }

        friend Element operator / (Element& a, int32_t b)
        {
            Element res(a.arithmetic());
            res.arithmetic().div(a, b, res);
            return res;
        }

        friend Element operator / (Element&& a, int32_t b)
        {
            Element res(std::move(a));
            res /= b;
            return res;
        }

        Element& operator %= (Element& a)
        {
            arithmetic().mod((Element&)*this, a, (Element&)*this);
            return (Element&)*this;
        }

        friend Element operator % (Element& a, Element& b)
        {
            Element res(a.arithmetic());
            res.arithmetic().mod(a, b, res);
            return res;
        }

        friend Element operator % (Element&& a, Element& b)
        {
            Element res(std::move(a));
            res %= ((void*)&a == (void*)&b) ? res : b;
            return res;
        }

        friend Element operator % (Element& a, Element&& b)
        {
            return a % b;
        }

        friend Element operator % (Element&& a, Element&& b)
        {
            return std::move(a) % b;
        }

        Element& gcd(Element& a)
        {
            arithmetic().gcd((Element&)*this, a, (Element&)*this);
            return (Element&)*this;
        }

        friend Element gcd(Element& a, Element& b)
        {
            Element res(a.arithmetic());
            res.arithmetic().gcd(a, b, res);
            return res;
        }

        friend Element gcd(Element&& a, Element& b)
        {
            Element res(std::move(a));
            res.gcd(((void*)&a == (void*)&b) ? res : b);
            return res;
        }

        friend Element gcd(Element& a, Element&& b)
        {
            Element res(std::move(b));
            res.gcd(((void*)&a == (void*)&b) ? res : a);
            return res;
        }

        friend Element gcd(Element&& a, Element&& b)
        {
            return gcd(std::move(a), b);
        }

        Element& inv(Element& n)
        {
            arithmetic().inv((Element&)*this, n, (Element&)*this);
            return (Element&)*this;
        }

        friend Element inv(Element& a, Element& n)
        {
            Element res(a.arithmetic());
            res.arithmetic().inv(a, n, res);
            return res;
        }

        friend Element inv(Element&& a, Element& b)
        {
            Element res(std::move(a));
            res.inv(((void*)&a == (void*)&b) ? res : b);
            return res;
        }

        friend Element inv(Element& a, Element&& b)
        {
            return inv(a, b);
        }

        friend Element inv(Element&& a, Element&& b)
        {
            return inv(std::move(a), b);
        }

        friend void swap(Element& a, Element& b)
        {
            Element t(std::move(a));
            a = std::move(b);
            b = std::move(t);
        }

    private:
        Arithmetic& _arithmetic;
    };

    class Giant;
    class GWNum;

    class GiantsArithmetic : public FieldArithmetic<Giant>
    {
        friend class Giant;

    public:
        GiantsArithmetic() { }

        virtual void alloc(Giant& a) override;
        virtual void alloc(Giant& a, int size);
        virtual void free(Giant& a) override;
        virtual void copy(const Giant& a, Giant& res) override;
        virtual void move(Giant&& a, Giant& res) override;
        virtual void init(int32_t a, Giant& res) override;
        virtual void init(uint32_t a, Giant& res);
        virtual void init(const std::string& a, Giant& res) override;
        virtual int cmp(const Giant& a, const Giant& b) override;
        virtual int cmp(const Giant& a, int32_t b) override;
        virtual void add(Giant& a, Giant& b, Giant& res) override;
        virtual void add(Giant& a, int32_t b, Giant& res) override;
        virtual void sub(Giant& a, Giant& b, Giant& res) override;
        virtual void sub(Giant& a, int32_t b, Giant& res) override;
        virtual void neg(Giant& a, Giant& res) override;
        virtual void mul(Giant& a, Giant& b, Giant& res) override;
        virtual void mul(Giant& a, int32_t b, Giant& res) override;
        virtual void mul(Giant& a, uint32_t b, Giant& res);
        virtual void div(Giant& a, Giant& b, Giant& res) override;
        virtual void div(Giant& a, int32_t b, Giant& res) override;
        virtual void gcd(Giant& a, Giant& b, Giant& res) override;
        virtual void inv(Giant& a, Giant& n, Giant& res) override;
        virtual void mod(Giant& a, Giant& b, Giant& res) override;
        virtual void shiftleft(Giant& a, int b, Giant& res);
        virtual void shiftright(Giant& a, int b, Giant& res);
        virtual int bitlen(const Giant& a);
        virtual bool bit(const Giant& a, int b);
        virtual void power(Giant& a, int32_t b, Giant& res);

        static GiantsArithmetic& default_arithmetic();
    };

    class GWGiantsArithmetic : public GiantsArithmetic
    {
        friend class GWState;

    public:
        GWGiantsArithmetic(gwhandle *gwdata) : _gwdata(gwdata) { }

        virtual void alloc(Giant& a) override;
        virtual void alloc(Giant& a, int size) override;
        virtual void free(Giant& a) override;

        int size() const { return _size; }

    private:
        gwhandle *_gwdata;
        int _size;
    };

    class Giant : public FieldElement<GiantsArithmetic, Giant>
    {
        friend class GiantsArithmetic;
        friend class GWGiantsArithmetic;

    public:
        Giant() : FieldElement<GiantsArithmetic, Giant>(GiantsArithmetic::default_arithmetic())
        {
            arithmetic().alloc(*this);
        }
        Giant(GiantsArithmetic& arithmetic) : FieldElement<GiantsArithmetic, Giant>(arithmetic)
        {
            arithmetic.alloc(*this);
        }
        Giant(GiantsArithmetic& arithmetic, int size) : FieldElement<GiantsArithmetic, Giant>(arithmetic)
        {
            arithmetic.alloc(*this, size);
        }
        ~Giant()
        {
            if (_giant != nullptr)
                arithmetic().free(*this);
        }
        Giant(const Giant& a) : FieldElement<GiantsArithmetic, Giant>(a.arithmetic())
        {
            arithmetic().copy(a, *this);
        }
        Giant(Giant&& a) : FieldElement<GiantsArithmetic, Giant>(a.arithmetic())
        {
            arithmetic().move(std::move(a), *this);
        }

        Giant& operator = (const Giant& a)
        {
            arithmetic().copy(a, *this);
            return *this;
        }
        Giant& operator = (Giant&& a)
        {
            arithmetic().move(std::move(a), *this);
            return *this;
        }
        using FieldElement<GiantsArithmetic, Giant>::operator=;

        virtual std::string to_string() const override;
        virtual void to_GWNum(GWNum& a) const;
        virtual Giant& operator = (const GWNum& a);

        giant to_giant() { return _giant; }
        int size() const { return _size; }

        Giant& operator = (uint32_t a)
        {
            arithmetic().init(a, *this);
            return *this;
        }
        Giant& operator *= (uint32_t a)
        {
            arithmetic().mul(*this, a, *this);
            return *this;
        }
        using FieldElement<GiantsArithmetic, Giant>::operator*=;
        Giant& operator <<= (int a)
        {
            arithmetic().shiftleft(*this, a, *this);
            return *this;
        }
        Giant& operator >>= (int a)
        {
            arithmetic().shiftright(*this, a, *this);
            return *this;
        }
        Giant& operator /= (Giant&& a)
        {
            *this /= a;
            return *this;
        }
        using FieldElement<GiantsArithmetic, Giant>::operator/=;
        friend Giant operator / (Giant& a, Giant& b)
        {
            Giant res(a.arithmetic());
            res.arithmetic().div(a, b, res);
            return res;
        }
        friend Giant operator / (Giant&& a, Giant& b)
        {
            Giant res(std::move(a));
            res /= ((void*)&a == (void*)&b) ? res : b;
            return res;
        }
        friend Giant operator / (Giant& a, Giant&& b)
        {
            return a / b;
        }
        friend Giant operator / (Giant&& a, Giant&& b)
        {
            return std::move(a) / b;
        }
        int bitlen() const { return arithmetic().bitlen(*this); }
        bool bit(int b) const { return arithmetic().bit(*this, b); }
        Giant& power(int32_t a)
        {
            arithmetic().power(*this, a, *this);
            return *this;
        }
        friend Giant power(Giant& a, int32_t b)
        {
            Giant res(a.arithmetic());
            res.arithmetic().power(a, b, res);
            return res;
        }
        friend Giant power(Giant&& a, int32_t b)
        {
            Giant res(std::move(a));
            res.power(b);
            return res;
        }

    private:
        giant _giant = nullptr;
        int _size = 0;
    };

    class GWState
    {
    public:
        GWState() : giants(&handle)
        {
            gwinit(&handle);
        }
        ~GWState()
        {
            if (N != nullptr)
                delete N;
            gwdone(&handle);
        }

        void setup(int k, int b, int n, int c);
        void setup(Giant& g);

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
        ~GWArithmetic();

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

        using GWArithmetic::add;
        using GWArithmetic::sub;
        using GWArithmetic::mul;
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

        virtual void mul(GWNum& a, GWNum& b, GWNum& res, int options) override;
        virtual void addmul(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options) override;
        virtual void submul(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options) override;
        virtual void muladd(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options) override;
        virtual void mulsub(GWNum& a, GWNum& b, GWNum& c, GWNum& res, int options) override;
        virtual void mulmuladd(GWNum& a, GWNum& b, GWNum& c, GWNum& d, GWNum& res, int options) override;
        virtual void mulmulsub(GWNum& a, GWNum& b, GWNum& c, GWNum& d, GWNum& res, int options) override;

        void reset();
        void restart(int op = 0);

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
            return (a.arithmetic().popg() = a).gcd(b);
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

        virtual std::string to_string() const override;

    private:
        gwnum _gwnum = nullptr;
    };

    class ArithmeticException : public std::exception
    {
    public:
        ArithmeticException() { }
        ArithmeticException(const char *message) : std::exception(message) { }
    };

    class InvalidFFTDataException : public ArithmeticException
    {
    public:
        InvalidFFTDataException() : ArithmeticException("Invalid FFT data.") { }
    };

    class NoInverseException : public ArithmeticException
    {
    public:
        NoInverseException(Giant& divisor) : ArithmeticException("The inverse does not exist."), divisor(GiantsArithmetic::default_arithmetic())
        {
            this->divisor = divisor;
            if (this->divisor < 0)
                this->divisor.arithmetic().neg(this->divisor, this->divisor);
        }
    public:
        Giant divisor;
    };
}
