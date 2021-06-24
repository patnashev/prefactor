#pragma once

#include <chrono>
#include <stdexcept>
#include <memory>
#include "arithmetic.h"
#include "inputnum.h"
#include "file.h"

class TaskRestartException : public std::exception
{
public:
    TaskRestartException() { }
};

class TaskAbortException : public std::exception
{
public:
    TaskAbortException() { }
};

class TaskState
{
public:
    TaskState(int iteration) : _iteration(iteration) { }
    virtual ~TaskState() { }
    
    virtual bool read(Reader& reader);
    virtual void write(Writer& writer);

    int iteration() { return _iteration; }
    char version() { return 0; }

protected:
    int _iteration;
};

template<class State>
State* read_state(File* file)
{
    State* state = nullptr;
    if (file != nullptr)
    {
        state = new State();
        if (!file->read(*state))
        {
            delete state;
            state = nullptr;
        }
    }
    return state;
}

class Task
{
public:
    const int MULS_PER_STATE_UPDATE = 10000;
    const int DISK_WRITE_TIME = 1;

public:
    Task(bool error_check) : _error_check(error_check){ }
    virtual ~Task() { }

    virtual void run();

    arithmetic::GWArithmetic& gw() { return *_gw; }
    TaskState* state() { return _state.get(); }
    int iterations() { return _iterations; }
    int state_update_period() { return _state_update_period; }
    bool is_last(int iteration) { return iteration + 1 == _iterations || (iteration + 1)%_state_update_period == 0; }

protected:
    virtual void init(arithmetic::GWState& gwstate, int iterations, File* file, TaskState* state);
    virtual void set_state(TaskState* state);
    virtual void setup() = 0;
    virtual void execute() = 0;
    virtual void reinit_gwstate() = 0;
    virtual void release() = 0;

    void check();
    void commit_setup();
    template<class TState, class... Args>
    void commit_execute(int iteration, Args&&... args)
    {
        if (iteration%state_update_period() == 0 || iteration == iterations())
        {
            check();
            set_state(new TState(iteration, std::forward<Args>(args)...));
        }
    }

protected:
    bool _error_check;
    arithmetic::GWState* _gwstate = nullptr;
    arithmetic::GWArithmetic* _gw = nullptr;
    std::unique_ptr<TaskState> _state;
    int _iterations = 0;
    int _state_update_period = MULS_PER_STATE_UPDATE;
    File* _file;
    std::chrono::system_clock::time_point _last_write;
private:
    int _restart_op = 0;
};