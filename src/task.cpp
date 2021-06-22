
#include "gwnum.h"
#include "task.h"
#include "exception.h"

using namespace arithmetic;

void Task::init(arithmetic::GWState& gwstate, int iterations)
{
    _gwstate = &gwstate;
    _iterations = iterations;
    _gw = nullptr;
}

void Task::run()
{
    ReliableGWArithmetic* reliable;
    std::unique_ptr<arithmetic::GWArithmetic> gw_setup(_error_check ? new ReliableGWArithmetic(*_gwstate) : new GWArithmetic(*_gwstate));
    std::unique_ptr<arithmetic::GWArithmetic> gw_execute(_error_check ? new ReliableGWArithmetic(*_gwstate) : new GWArithmetic(*_gwstate));

    _restart_op = 0;
    int restart_count = 0;
    while (true)
    {
        _gw = gw_setup.get();
        reliable = dynamic_cast<ReliableGWArithmetic*>(_gw);
        if (reliable)
            reliable->reset();
        int setup_restart_count = 0;
        bool setup_failed = false;
        while (true)
        {
            try
            {
                setup();
                break;
            }
            catch (const TaskRestartException&)
            {
                if (!reliable)
                {
                    release();
                    _error_check = true;
                    gw_setup.reset(new ReliableGWArithmetic(*_gwstate));
                    gw_execute.reset(new ReliableGWArithmetic(*_gwstate));
                    _gw = gw_setup.get();
                    reliable = dynamic_cast<ReliableGWArithmetic*>(_gw);
                }
            }
            catch (const ArithmeticException& e)
            {
                printf("ArithmeticException: %s\n", e.what());
            }
            //GWASSERT(0);
            if (reliable && reliable->restart_flag() && !reliable->failure_flag())
            {
                reliable->restart();
                continue;
            }
            if ((!reliable || !reliable->failure_flag()) && setup_restart_count < 5)
            {
                setup_restart_count++;
                continue;
            }
            setup_failed = true;
            break;
        }

        _gw = gw_execute.get();
        reliable = dynamic_cast<ReliableGWArithmetic*>(_gw);

        if (!setup_failed)
        {
            try
            {
                execute();
                break;
            }
            catch (const TaskRestartException&)
            {
                if (!reliable)
                {
                    release();
                    _error_check = true;
                    gw_setup.reset(new ReliableGWArithmetic(*_gwstate));
                    gw_execute.reset(new ReliableGWArithmetic(*_gwstate));
                }
            }
            catch (const ArithmeticException& e)
            {
                printf("ArithmeticException: %s\n", e.what());
            }
            //GWASSERT(0);
            if (reliable && reliable->restart_flag() && !reliable->failure_flag())
            {
                reliable->restart(_restart_op);
                continue;
            }
            if ((!reliable || !reliable->failure_flag()) && restart_count < 5)
            {
                restart_count++;
                continue;
            }
        }

        release();
        if (_gwstate->next_fft_count >= 5)
        {
            printf("Maximum FFT increment reached.\n");
            throw TaskAbortException();
        }
        _gwstate->next_fft_count++;
        reinit_gwstate();
        restart_count = 0;
        _restart_op = 0;
        if (reliable)
            reliable->reset();
    }
    release();
    _gw = nullptr;
}

void Task::check()
{
    if (_error_check)
    {
        ReliableGWArithmetic* reliable = dynamic_cast<ReliableGWArithmetic*>(_gw);
        if (reliable->restart_flag() || reliable->failure_flag())
            throw TaskRestartException();
    }
}

void Task::set_state(TaskState* state)
{
    _state.reset(state);
    if (_error_check && _gw != nullptr)
    {
        ReliableGWArithmetic* reliable = dynamic_cast<ReliableGWArithmetic*>(_gw);
        _restart_op = reliable->op();
    }
}

void Task::commit_setup()
{
    check();
    if (_error_check)
    {
        ReliableGWArithmetic* reliable = dynamic_cast<ReliableGWArithmetic*>(_gw);
        _restart_op = reliable->op();
    }
}
