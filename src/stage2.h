#pragma once

#include "arithmetic.h"
#include "group.h"
#include "lucas.h"
#include "montgomery.h"
#include "edwards.h"
#include "primelist.h"
#include "inputnum.h"
#include "task.h"

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
        State() : TaskState(0) { _G = 1; }
        State(int iteration, arithmetic::GWNum& G) : TaskState(iteration) { _G = G; }
        arithmetic::Giant& G() { return _G; }

    private:
        arithmetic::Giant _G;
    };

public:
    Stage2(PrimeList& primes, int B1, int B2, int D, int A, int L) : Task(true), _primes(primes), _B1(B1), _B2(B2), _D(D), _A(A), _L(L)
    {
        _pairing = get_pairing(primes, B1, B2, D, A, L, true);
    }

    static Pairing get_pairing(PrimeList& primes, int B1, int B2, int D, int A, int L, bool with_distances);
    template<class Element>
    int precompute(arithmetic::DifferentialGroupArithmetic<Element>& arithmetic, Element& X1, Element& XD, Element& XDA, std::vector<std::unique_ptr<Element>>& precomp);
    void init(InputNum& input, arithmetic::GWState& gwstate);
    void reinit_gwstate() override;
    void done(const arithmetic::Giant& factor);

    PrimeList& primes() { return _primes; }
    bool success() { return _success; }
    State* state() { return static_cast<State*>(Task::state()); }

private:
    PrimeList& _primes;

protected:
    int _B1;
    int _B2;
    int _D;
    int _A;
    int _L;
    Pairing _pairing;

    InputNum* _input = nullptr;
    double _timer = 0;
    int _transforms = 0;
    bool _success = false;
};

class PP1Stage2 : public Stage2
{
public:
    PP1Stage2(PrimeList& primes, int B1, int B2, int D, int A, int L) : Stage2(primes, B1, B2, D, A, L) { }

    void init(InputNum& input, arithmetic::GWState& gwstate, arithmetic::Giant& P, bool minus1);
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
    EdECMStage2(PrimeList& primes, int B1, int B2, int D, int L, int LN) : Stage2(primes, B1, B2, D, 1, L), _LN(LN) { }

    void init(InputNum& input, arithmetic::GWState& gwstate, arithmetic::Giant& X, arithmetic::Giant& Y, arithmetic::Giant& Z, arithmetic::Giant& T, arithmetic::Giant& EdD);
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
