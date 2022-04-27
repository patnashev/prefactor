#pragma once

#include <deque>
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

#ifdef _DEBUG
#define DEBUG_POLY_STAGE2
#endif

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

    struct Stage2Thread
    {
        gwthread id;
        std::unique_ptr<Stage2> stage2;
        std::unique_ptr<arithmetic::GWState> gwstate;
        File* file;
    };

protected:
    Stage2(uint64_t B1, uint64_t B2) : Task(true), _B1(B1), _B2(B2) { }

    void stage2_pairing(int D, int A, int L, Logging& logging)
    {
        _D = D;
        _A = A;
        _L = L;
        PrimeList primes((int)_B2 + 100);
        _pairing = get_pairing(logging, primes, (int)_B1, (int)_B2, D, A, L, true);
    }

    template<class T>
    void stage2_poly(int D, int L, int poly_degree, int poly_power, int poly_threads)
    {
        _D = D;
        _A = 1;
        _L = L;
        _poly_degree = poly_degree;
        _poly_power = poly_power;
        _poly_mod_degree = phi(D)/2;
        _poly_threads = poly_threads;
        _poly_thread_helpers.resize(poly_threads - 1);
        for (auto it = _poly_thread_helpers.begin(); it != _poly_thread_helpers.end(); it++)
            it->stage2.reset(new T(static_cast<T*>(this)));
        for (int i = 1; i < _poly_threads; i <<= 1)
        {
            _poly_thread_mergers.push_back(_poly_thread_helpers[i - 1].stage2.get());
            for (int j = 2*i; j + i < _poly_threads; j += 2*i)
                _poly_thread_helpers[j - 1].stage2->_poly_thread_mergers.push_back(_poly_thread_helpers[j + i - 1].stage2.get());
        }
    }

public:
    static Pairing get_pairing(Logging& logging, PrimeList& primes, int B1, int B2, int D, int A, int L, bool with_distances);

    bool success() { return _success; }
    State* state() { return static_cast<State*>(Task::state()); }

protected:
    template<class Element>
    int precompute(arithmetic::DifferentialGroupArithmetic<Element>& arithmetic, Element& X1, Element& XD, Element& XDA, std::vector<std::unique_ptr<Element>>& precomp);
    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, TaskState* state, Logging* logging);
    void reinit_gwstate() override;
    void done(const arithmetic::Giant& factor);

    bool is_poly() { return _poly_power > 0; }
    bool is_poly_threaded() { return _poly_threads > 1; }
    bool is_poly_helper() { return _poly_thread_main != nullptr; }
    int poly_degree() { return _poly_accumulator ? _poly_degree : _poly_mod_degree - 1; }
    int poly_power() { return _poly_power; }
    void poly_init();
    void poly_setup(std::vector<arithmetic::GWNum*>& roots);
    void poly_release();
    void poly_execute(std::vector<arithmetic::GWNum>& roots);
    void poly_merge(arithmetic::Poly& G);
    void poly_final(arithmetic::GWNum& G);
    void poly_threads_init();
    void poly_helper_done();
    void poly_threads_merge();
    friend void poly_thread(void* data);

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

    int _poly_degree = 0;
    int _poly_power = 0;
    int _poly_mod_degree = 0;
    std::deque<arithmetic::GWState> _poly_gwstate;
    std::deque<arithmetic::GWArithmetic> _poly_gw;
    std::vector<arithmetic::PolyMult> _poly_mult;
    std::vector<std::vector<arithmetic::Poly>> _poly_mod;
    std::unique_ptr<arithmetic::Poly> _poly_reciprocal;
    std::vector<std::vector<arithmetic::Poly>> _poly_prod;
    std::unique_ptr<arithmetic::Poly> _poly_accumulator;
    double _poly_timer;
#ifdef DEBUG_POLY_STAGE2
    std::vector<arithmetic::GWNum> _poly_rem;
#endif

    int _poly_threads = 1;
    Stage2* _poly_thread_main = nullptr;
    std::vector<Stage2Thread> _poly_thread_helpers;
    std::vector<Stage2*> _poly_thread_mergers;
    std::exception_ptr _poly_thread_exception;
    gwevent _poly_done;
    gwevent _poly_merged;
};

class PP1Stage2 : public Stage2
{
public:
    PP1Stage2(uint64_t B1, uint64_t B2) : Stage2(B1, B2) { }
    PP1Stage2(PP1Stage2* main) : Stage2(main->_B1, main->_B2)
    {
        _poly_thread_main = main;
        stage2_poly(main->_D, main->_L, main->_poly_degree, main->_poly_power, 1);
    }

    void stage2_pairing(int D, int A, int L, Logging& logging)
    {
        Stage2::stage2_pairing(D, A, L, logging);
    }

    void stage2_poly(int D, int L, int poly_degree, int poly_power, int poly_threads)
    {
        Stage2::stage2_poly<PP1Stage2>(D, L, poly_degree, poly_power, poly_threads);
    }

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
    EdECMStage2(uint64_t B1, uint64_t B2) : Stage2(B1, B2) { }
    EdECMStage2(EdECMStage2* main) : Stage2(main->_B1, main->_B2)
    {
        _poly_thread_main = main;
        stage2_poly(main->_D, main->_L, main->_LN, main->_poly_degree, main->_poly_power, 1);
    }

    void stage2_pairing(int D, int L, int LN, Logging& logging)
    {
        _LN = LN;
        Stage2::stage2_pairing(D, 1, L, logging);
    }

    void stage2_poly(int D, int L, int LN, int poly_degree, int poly_power, int poly_threads)
    {
        _LN = LN;
        Stage2::stage2_poly<EdECMStage2>(D, L, poly_degree, poly_power, poly_threads);
    }

    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging, arithmetic::Giant& X, arithmetic::Giant& Y, arithmetic::Giant& Z, arithmetic::Giant& T, arithmetic::Giant& EdD);

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
