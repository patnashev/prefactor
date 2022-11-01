#pragma once

#include "arithmetic.h"
#include "group.h"
#include "lucas.h"
#include "montgomery.h"
#include "edwards.h"
#include "integer.h"
#include "inputnum.h"
#include "task.h"
#include "file.h"

class Stage1 : public InputTask
{
public:

public:
    Stage1(uint64_t B1) : _B1(B1)
    {
    }

    static arithmetic::Giant get_stage1_exp(int B1);
    static arithmetic::Giant get_stage1_exp(uint64_t B1cur, uint64_t B1next, uint64_t B1max);

    void set_no_check_success() { _check_success = false; }
    bool success() { return _success; }
    int exp_len() { return _exp_len; }

protected:
    void setup() override { }
    void release() override { }
    void reinit_gwstate() override;
    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, TaskState* state, Logging* logging, int iterations);
    void done(const arithmetic::Giant& factor);

protected:
    uint64_t _B1;
    int _squarings = 0;
    int _exp_len = 0;

    double _timer = 0;
    int _transforms = 0;
    bool _success = false;
    bool _check_success = true;
};

class PM1Stage1 : public Stage1
{
public:
    class State : public TaskState
    {
    public:
        State() : TaskState(1) { }
        void set(int iteration, arithmetic::GWNum& X) { TaskState::set(iteration); _X = X; }
        arithmetic::Giant& X() { return _X; }
        bool read(Reader& reader) override { return TaskState::read(reader) && reader.read(_X); }
        void write(Writer& writer) override { TaskState::write(writer); writer.write(_X); }

    private:
        arithmetic::Giant _X;
    };

public:
    PM1Stage1(int B1);

    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging);

    State* state() { return static_cast<State*>(Task::state()); }
    arithmetic::Giant& V() { return _V; }

protected:
    void execute() override;

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
        State() : TaskState(2) { }
        State(arithmetic::Giant& P) : TaskState(2) { _V = P; }
        void set(int iteration, arithmetic::LucasV& V) { TaskState::set(iteration); _V = V.V(); }
        arithmetic::Giant& V() { return _V; }
        bool read(Reader& reader) override { return TaskState::read(reader) && reader.read(_V); }
        void write(Writer& writer) override { TaskState::write(writer); writer.write(_V); }

    private:
        arithmetic::Giant _V;
    };
    
public:
    PP1Stage1(int B1, std::string& sP);

    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging);

    State* state() { return static_cast<State*>(Task::state()); }

protected:
    void execute() override;

protected:
    arithmetic::PrimeList _primes;

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
        State() : TaskState(3) { }
        void set(int iteration, arithmetic::EdPoint& p) { TaskState::set(iteration); p.serialize(_X, _Y, _Z, _T); }
        arithmetic::Giant& X() { return _X; }
        arithmetic::Giant& Y() { return _Y; }
        arithmetic::Giant& Z() { return _Z; }
        arithmetic::Giant& T() { return _T; }
        bool read(Reader& reader) override { return TaskState::read(reader) && reader.read(_X) && reader.read(_Y) && reader.read(_Z) && reader.read(_T); }
        void write(Writer& writer) override { TaskState::write(writer); writer.write(_X); writer.write(_Y); writer.write(_Z); writer.write(_T); }

    private:
        arithmetic::Giant _X;
        arithmetic::Giant _Y;
        arithmetic::Giant _Z;
        arithmetic::Giant _T;
    };

public:
    EdECMStage1(int B1, int W);
    EdECMStage1(uint64_t B1cur, uint64_t B1next, uint64_t B1max, int max_size);

    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging, arithmetic::Giant* X, arithmetic::Giant* Y, arithmetic::Giant* Z, arithmetic::Giant* T, arithmetic::Giant* EdD);

    State* state() { return static_cast<State*>(Task::state()); }
    int W() { return _W; }

protected:
    void setup() override;
    void release() override;
    void execute() override;

private:
    int _W;
    std::vector<int16_t> _NAF_W;
    arithmetic::Giant* _X;
    arithmetic::Giant* _Y;
    arithmetic::Giant* _Z;
    arithmetic::Giant* _T;
    arithmetic::Giant* _EdD;
    std::unique_ptr<arithmetic::EdwardsArithmetic> ed;
    std::unique_ptr<arithmetic::GWNum> _ed_d;
    std::vector<std::unique_ptr<arithmetic::EdPoint>> _u;
};
