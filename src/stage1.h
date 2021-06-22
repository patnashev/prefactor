#pragma once

#include "arithmetic.h"
#include "group.h"
#include "lucas.h"
#include "montgomery.h"
#include "edwards.h"
#include "primelist.h"
#include "inputnum.h"
#include "task.h"

class Stage1 : public Task
{
public:

public:
    Stage1(PrimeList& primes, int B1) : Task(true), _primes(primes), _B1(B1)
    {
    }

    void setup() override { }
    void release() override { }
    void reinit_gwstate() override;
    void init(InputNum& input, arithmetic::GWState& gwstate, int iterations);
    void done(const arithmetic::Giant& factor);
    arithmetic::Giant get_stage1_exp();

    PrimeList& primes() { return _primes; }
    bool success() { return _success; }

private:
    PrimeList& _primes;

protected:
    int _B1;

    InputNum* _input = nullptr;
    double _timer = 0;
    int _transforms = 0;
    bool _success = false;
};

class PM1Stage1 : public Stage1
{
public:
    class State : public TaskState
    {
    public:
        State(int iteration, arithmetic::GWNum& X) : TaskState(iteration) { _X = X; }
        arithmetic::Giant& X() { return _X; }

    private:
        arithmetic::Giant _X;
    };

public:
    PM1Stage1(PrimeList& primes, int B1);

    void init(InputNum& input, arithmetic::GWState& gwstate);
    void execute() override;

    State* state() { return static_cast<State*>(Task::state()); }
    arithmetic::Giant& V() { return _V; }

private:
    arithmetic::Giant _exp;
    arithmetic::Giant _V;
};

class PP1Stage1 : public Stage1
{
public:
    class State : public TaskState
    {
    public:
        State(arithmetic::Giant& P) : TaskState(0) { _V = P; }
        State(int iteration, arithmetic::LucasV& V) : TaskState(iteration) { _V = V.V(); }
        arithmetic::Giant& V() { return _V; }

    private:
        arithmetic::Giant _V;
    };
    
public:
    PP1Stage1(PrimeList& primes, int B1, std::string& sP);

    void init(InputNum& input, arithmetic::GWState& gwstate);
    void execute() override;

    State* state() { return static_cast<State*>(Task::state()); }

private:
    std::string _sP;
    int _Pa;
    int _Pb;
    arithmetic::Giant _P;
};

class EdECMStage1 : public Stage1
{
public:
    class State : public TaskState
    {
    public:
        State(int iteration, arithmetic::EdPoint& p) : TaskState(iteration) { p.serialize(_X, _Y, _Z, _T); }
        arithmetic::Giant& X() { return _X; }
        arithmetic::Giant& Y() { return _Y; }
        arithmetic::Giant& Z() { return _Z; }
        arithmetic::Giant& T() { return _T; }

    private:
        arithmetic::Giant _X;
        arithmetic::Giant _Y;
        arithmetic::Giant _Z;
        arithmetic::Giant _T;
    };

public:
    EdECMStage1(PrimeList& primes, int B1, int W);

    void init(InputNum& input, arithmetic::GWState& gwstate, arithmetic::Giant& X, arithmetic::Giant& Y, arithmetic::Giant& Z, arithmetic::Giant& T, arithmetic::Giant& EdD);
    void setup() override;
    void release() override;
    void execute() override;

    State* state() { return static_cast<State*>(Task::state()); }

private:
    int _W;
    std::vector<int16_t> _NAF_W;
    arithmetic::Giant _X;
    arithmetic::Giant _Y;
    arithmetic::Giant _Z;
    arithmetic::Giant _T;
    arithmetic::Giant _EdD;
    std::unique_ptr<arithmetic::EdwardsArithmetic> ed;
    std::vector<std::unique_ptr<arithmetic::EdPoint>> _u;
};
