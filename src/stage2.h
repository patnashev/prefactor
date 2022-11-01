#pragma once

#include <deque>
#include "arithmetic.h"
#include "group.h"
#include "lucas.h"
#include "montgomery.h"
#include "edwards.h"
#include "integer.h"
#include "inputnum.h"
#include "task.h"
#include "file.h"
#include "logging.h"

class Stage2 : public InputTask
{
protected:
    Stage2(uint64_t B1, uint64_t B2, int D) : _B1(B1), _B2(B2), _D(D) { }

public:
    bool success() { return _success; }
    std::string& res64() { return _res64; }

protected:
    virtual void init(InputNum* input, arithmetic::GWState* gwstate, File* file, TaskState* state, Logging* logging, int iterations);
    virtual void done(const arithmetic::Giant& factor);

protected:
    uint64_t _B1;
    uint64_t _B2;
    int _D;

    double _timer = 0;
    int _transforms = 0;
    bool _success = false;
    std::string _res64;
};

class Stage2Pairing : public Stage2
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

protected:
    Stage2Pairing(uint64_t B1, uint64_t B2, int D, int A, int L, Logging& logging) : Stage2(B1, B2, D)
    {
        _A = A;
        _L = L;
        arithmetic::PrimeList primes((int)_B2 + 100);
        _pairing = get_pairing(logging, primes, (int)_B1, (int)_B2, D, A, L, true);
    }

public:
    static Pairing get_pairing(Logging& logging, arithmetic::PrimeList& primes, int B1, int B2, int D, int A, int L, bool with_distances);

    State* state() { return static_cast<State*>(Task::state()); }

protected:
    template<class Element>
    int precompute(arithmetic::DifferentialGroupArithmetic<Element>& arithmetic, Element& X1, Element& XD, Element& XDA, std::vector<std::unique_ptr<Element>>& precomp);

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
    std::string _res64;
};

class IPP1Stage2
{
public:
    virtual void init(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging, arithmetic::Giant& P, bool minus1) = 0;
};

class PP1Stage2 : public Stage2Pairing, public IPP1Stage2
{
public:
    PP1Stage2(uint64_t B1, uint64_t B2, int D, int A, int L, Logging& logging) : Stage2Pairing(B1, B2, D, A, L, logging) { }

    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging, arithmetic::Giant& P, bool minus1) override;

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

class IEdECMStage2
{
public:
    virtual void init(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging, arithmetic::Giant& X, arithmetic::Giant& Y, arithmetic::Giant& Z, arithmetic::Giant& T, arithmetic::Giant& EdD) = 0;
};

class EdECMStage2 : public Stage2Pairing, public IEdECMStage2
{
public:
    EdECMStage2(uint64_t B1, uint64_t B2, int D, int L, int LN, Logging& logging) : Stage2Pairing(B1, B2, D, 1, L, logging)
    {
        _LN = LN;
    }

    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging, arithmetic::Giant& X, arithmetic::Giant& Y, arithmetic::Giant& Z, arithmetic::Giant& T, arithmetic::Giant& EdD) override;

protected:
    void setup() override;
    void execute() override;
    void release() override;

protected:
    int _LN;

private:
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
