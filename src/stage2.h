#pragma once

#include "arithmetic.h"
#include "group.h"
#include "lucas.h"
#include "montgomery.h"
#include "edwards.h"
#include "primelist.h"
#include "inputnum.h"
#include "task.h"
#include "file.h"
#include "logging.h"
#include "poly.h"

class Stage2 : public Task
{
public:
    template<int TL>
    struct prime
    {
        int value;
        int match;
        int adjacency[TL];

        prime(int val) { value = val; match = 0; }
    };

    struct Pairing
    {
        int total;
        int pairs;
        int first_D;
        int last_D;
        std::vector<int> distances;
    };

    class State : public TaskState
    {
    public:
        State() : TaskState(10) { _G = 1; }
        void set(int iteration, arithmetic::GWNum& G) { TaskState::set(iteration); _G = G; }
        arithmetic::Giant& G() { return _G; }
        bool read(Reader& reader) override { return TaskState::read(reader) && reader.read(_G); }
        void write(Writer& writer) override { TaskState::write(writer); writer.write(_G); }

    private:
        arithmetic::Giant _G;
    };

public:
    Stage2(Logging& logging, PrimeList& primes, uint64_t B1, uint64_t B2, int D, int A, int L, int poly_power) : Task(true), _primes(primes), _B1(B1), _B2(B2), _D(D), _A(A), _L(L), _poly_power(poly_power)
    {
        if (poly_power == 0)
            _pairing = get_pairing(logging, primes, (int)B1, (int)B2, D, A, L, true);
        else
        {
            _pairing.first_D = (int)(B1/D + 1);
            _pairing.last_D = (int)((B2 + D - 1)/D);
        }
    }

    static Pairing get_pairing(Logging& logging, PrimeList& primes, int B1, int B2, int D, int A, int L, bool with_distances);

    PrimeList& primes() { return _primes; }
    bool success() { return _success; }
    State* state() { return static_cast<State*>(Task::state()); }

protected:
    template<class Element>
    int precompute(arithmetic::DifferentialGroupArithmetic<Element>& arithmetic, Element& X1, Element& XD, Element& XDA, std::vector<std::unique_ptr<Element>>& precomp);
    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, TaskState* state, Logging* logging);
    void reinit_gwstate() override;
    void done(const arithmetic::Giant& factor);

    bool is_poly() { return _poly_power > 0; }
    int poly_power() { return _poly_power; }
    void poly_setup(std::vector<arithmetic::GWNum*>& roots, int LN);
    void poly_release();
    void poly_execute(std::vector<arithmetic::GWNum*>& roots, arithmetic::GWNum& G);

private:
    PrimeList& _primes;

protected:
    uint64_t _B1;
    uint64_t _B2;
    int _D;
    int _A;
    int _L;
    Pairing _pairing;

    InputNum* _input = nullptr;
    double _timer = 0;
    int _transforms = 0;
    bool _success = false;

    int _poly_power;
    std::unique_ptr<arithmetic::PolyMulOptimal> _poly_mul;
    std::vector<std::vector<arithmetic::Poly>> _poly_mod;
    std::vector<std::vector<arithmetic::Poly>> _poly_prod;
};

class PP1Stage2 : public Stage2
{
public:
    PP1Stage2(Logging& logging, PrimeList& primes, uint64_t B1, uint64_t B2, int D, int A, int L, int poly_power) : Stage2(logging, primes, B1, B2, D, A, L, poly_power) { }

    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging, arithmetic::Giant& P, bool minus1);

protected:
    void setup() override;
    void execute() override;
    void release() override;

private:
    arithmetic::Giant _P;

    std::unique_ptr<arithmetic::LucasVArithmetic> lucas;
    std::vector<std::unique_ptr<arithmetic::LucasV>> _precomp;
    std::unique_ptr<arithmetic::LucasV> _W;
    std::unique_ptr<arithmetic::LucasV> _Wa;
    std::unique_ptr<arithmetic::LucasV> _Vn;
    std::unique_ptr<arithmetic::LucasV> _Vn1;
    std::unique_ptr<arithmetic::LucasV> _Va;
    std::unique_ptr<arithmetic::LucasV> _Va1;
};

class EdECMStage2 : public Stage2
{
public:
    EdECMStage2(Logging& logging, PrimeList& primes, uint64_t B1, uint64_t B2, int D, int L, int LN, int poly_power) : Stage2(logging, primes, B1, B2, D, 1, L, poly_power), _LN(LN) { }

    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging, arithmetic::Giant& X, arithmetic::Giant& Y, arithmetic::Giant& Z, arithmetic::Giant& T, arithmetic::Giant& EdD);

protected:
    void setup() override;
    void execute() override;
    void release() override;

private:
    int _LN;
    arithmetic::Giant _X;
    arithmetic::Giant _Y;
    arithmetic::Giant _Z;
    arithmetic::Giant _T;
    arithmetic::Giant _EdD;

    std::unique_ptr<arithmetic::GWNum> _ed_d;
    std::unique_ptr<arithmetic::MontgomeryArithmetic> montgomery;
    std::vector<std::unique_ptr<arithmetic::EdY>> _precomp;
    std::unique_ptr<arithmetic::EdY> _W;
    std::unique_ptr<arithmetic::EdY> _Pn;
    std::unique_ptr<arithmetic::EdY> _Pn1;
};
