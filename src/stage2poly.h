#pragma once

#include <condition_variable>
#include <thread>

#include "stage2.h"
#include "poly.h"

template<class Element>
class Stage2Poly : public Stage2
{
public:

    class State : public TaskState
    {
    public:
        State() : TaskState(11) { _G = 1; }
        void set(int iteration, arithmetic::GWNum& G) { TaskState::set(iteration); _G = G; }
        arithmetic::Giant& G() { return _G; }
        bool read(Reader& reader) override { return TaskState::read(reader) && reader.read(_G); }
        void write(Writer& writer) override { TaskState::write(writer); writer.write(_G); }

    private:
        arithmetic::Giant _G;
    };

protected:

    class SmallPolyWorker
    {
    public:
        SmallPolyWorker(Stage2Poly& stage2);
        virtual ~SmallPolyWorker();

        void run();
        virtual void elements_to_gwnums(std::vector<std::unique_ptr<Element>>& elements, gwnum* data) = 0;

        virtual typename Element::Arithmetic& arithmetic() = 0;
        arithmetic::GWArithmetic& gw() { return _gw; }
        arithmetic::PolyMult* poly_mult() { return _poly_mult.data(); }
        arithmetic::GWNum& G() { return *_G; }
        std::unique_ptr<arithmetic::GWNum>& check() { return _check; }
        void set_main(bool value) { _main = value; }

    protected:
        bool _main = false;
        Stage2Poly<Element>& _stage2;
        arithmetic::GWState _gwstate;
        arithmetic::ThreadSafeGWArithmetic _gw;
        std::unique_ptr<arithmetic::GWNum> _G;
        std::unique_ptr<arithmetic::GWNum> _check;

        std::deque<arithmetic::GWState> _poly_gwstate;
        std::deque<arithmetic::GWArithmetic> _poly_gw;
        std::vector<arithmetic::PolyMult> _poly_mult;
    };

protected:
    Stage2Poly(uint64_t B1, uint64_t B2, int D, int poly_power, int smallpoly_power, int threads, bool check) : Stage2(B1, B2, D)
    {
        _poly_power = poly_power;
        _poly_mod_degree = phi(D)/2;
        _poly_threads = threads;
        _poly_check = check;
        _smallpoly_power = smallpoly_power;
        _tinypoly_power = _smallpoly_power;
    }

public:
    State* state() { return static_cast<State*>(Task::state()); }
    int D() { return _D; }
    virtual typename Element::Arithmetic& arithmetic() = 0;

protected:

    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, TaskState* state, Logging* logging);
    void setup() override;
    void execute() override;
    void release() override;
    void done(const arithmetic::Giant& factor) override;

    int poly_power() { return _poly_power; }
    bool poly_check() { return _poly_check; }

private:
    void smallpoly_init(int count, gwarray data, int data_level);
    void smallpoly_init_level(gwarray data, int data_level);

protected:
    int _first_D;
    int _last_D;

    int _poly_power;
    int _smallpoly_power;
    int _tinypoly_power;
    int _poly_mod_degree;
    int _poly_optmem = 10;
    int _poly_optmem_small = 10;
    int _poly_threads = 1;
    bool _poly_check = false;
    int _poly_check_index;

    std::deque<arithmetic::GWState> _poly_gwstate;
    std::deque<arithmetic::GWArithmetic> _poly_gw;
    std::vector<arithmetic::PolyMult> _poly_mult;

    std::vector<gwarray> _gwarrays;
    std::vector<gwarray> _gwarrays_tmp;
    int _smallpoly_count;
    std::vector<std::vector<std::vector<arithmetic::Poly>>> _smallpoly_alloc;
    std::vector<std::vector<std::vector<arithmetic::Poly>>> _smallpoly;
    std::vector<std::vector<arithmetic::Poly>> _poly_mod;
    std::vector<std::vector<std::vector<arithmetic::Poly>>> _smallpoly_mod;
    std::vector<arithmetic::Poly> _smallpoly_rem;

    std::unique_ptr<arithmetic::Poly> _modulus;
    std::unique_ptr<arithmetic::Poly> _reciprocal;
    std::unique_ptr<arithmetic::Poly> _accumulator;

    std::vector<std::thread> _threads;
    std::exception_ptr _thread_exception;
    std::vector<std::unique_ptr<SmallPolyWorker>> _workers;
    int _workstage = 0;
    std::mutex _work_mutex;
    std::condition_variable _workstage_signal;
    std::condition_variable _workqueue_signal;
    std::condition_variable _workdone_signal;

    struct WorkItem
    {
        std::unique_ptr<Element> Xn;
        std::unique_ptr<Element> Xn1;
        std::unique_ptr<Element> Xdn;
        std::unique_ptr<Element> Xdn1;
        int count;
        int n;
        int distance;
    };
    std::vector<WorkItem> _workqueue;
    Element* _X1;
    Element* _Xd;
};

class PP1Stage2Poly : public Stage2Poly<arithmetic::LucasV>, public IPP1Stage2
{
public:
    class SmallPolyWorker : public Stage2Poly<arithmetic::LucasV>::SmallPolyWorker
    {
    public:
        SmallPolyWorker(PP1Stage2Poly& stage2) : Stage2Poly<arithmetic::LucasV>::SmallPolyWorker(stage2), _lucas(_gw) { }

        virtual void elements_to_gwnums(std::vector<std::unique_ptr<arithmetic::LucasV>>& elements, gwnum* data) override;
        virtual arithmetic::LucasVArithmetic& arithmetic() override { return _lucas; }

    protected:
        arithmetic::LucasVArithmetic _lucas;
    };

public:
    PP1Stage2Poly(uint64_t B1, uint64_t B2, int D, int poly_power, int smallpoly_power, int threads, bool check = false) : Stage2Poly(B1, B2, D, poly_power, smallpoly_power, threads, check)
    {
    }

    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging, arithmetic::Giant& P, bool minus1);

    virtual arithmetic::LucasVArithmetic& arithmetic() override { return *_lucas; }

protected:
    void setup() override;
    void release() override;

private:
    arithmetic::Giant _P;

    std::unique_ptr<arithmetic::ThreadSafeGWArithmetic> _safe_gw;
    std::unique_ptr<arithmetic::LucasVArithmetic> _lucas;
    std::unique_ptr<arithmetic::LucasV> _W;
    std::unique_ptr<arithmetic::LucasV> _Wd;
};

class EdECMStage2Poly : public Stage2Poly<arithmetic::EdY>, public IEdECMStage2
{
public:
    class SmallPolyWorker : public Stage2Poly<arithmetic::EdY>::SmallPolyWorker
    {
    public:
        SmallPolyWorker(EdECMStage2Poly& stage2) : Stage2Poly<arithmetic::EdY>::SmallPolyWorker(stage2), _montgomery(_gw, *stage2._ed_d) { }

        virtual void elements_to_gwnums(std::vector<std::unique_ptr<arithmetic::EdY>>& elements, gwnum* data) override;
        virtual arithmetic::MontgomeryArithmetic& arithmetic() override { return _montgomery; }

    protected:
        arithmetic::MontgomeryArithmetic _montgomery;
    };

public:
    EdECMStage2Poly(uint64_t B1, uint64_t B2, int D, int poly_power, int smallpoly_power, int threads, bool check = false) : Stage2Poly(B1, B2, D, poly_power, smallpoly_power, threads, check)
    {
    }

    void init(InputNum* input, arithmetic::GWState* gwstate, File* file, Logging* logging, arithmetic::Giant& X, arithmetic::Giant& Y, arithmetic::Giant& Z, arithmetic::Giant& T, arithmetic::Giant& EdD);

    virtual arithmetic::MontgomeryArithmetic& arithmetic() override { return *_montgomery; }

protected:
    void setup() override;
    void execute() override;
    void release() override;

private:
    arithmetic::Giant _X;
    arithmetic::Giant _Y;
    arithmetic::Giant _Z;
    arithmetic::Giant _T;
    arithmetic::Giant _EdD;

    std::unique_ptr<arithmetic::GWNum> _ed_d;
    std::unique_ptr<arithmetic::ThreadSafeGWArithmetic> _safe_gw;
    std::unique_ptr<arithmetic::MontgomeryArithmetic> _montgomery;
    std::unique_ptr<arithmetic::EdY> _W;
    std::unique_ptr<arithmetic::EdY> _Wd;
};
