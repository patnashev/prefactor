#pragma once

#include <memory>
#include <vector>
#include "inputnum.h"

class Progress
{
public:
    Progress() { }

    void add_stage(int cost) { _total_cost += cost; _costs.push_back(cost); }
    void next_stage() { _cur_stage++; update(0, 0); }
    void update(double progress, int fft_count);
    void time_init(double elapsed);

    double progress_stage() { return _cur_progress; }
    double progress_total() { int cost = 0; for (int i = 0; i < _cur_stage; cost += _costs[i], i++); return (cost + _costs[_cur_stage]*_cur_progress)/_total_cost; }
    int cost_total() { return _total_cost; }
    double time_total() { return _time_total; }
    double time_fft() { return _time_fft; }
    int fft_count() { return _fft_count; }

private:
    std::vector<int> _costs;
    int _total_cost = 0;
    int _cur_stage = 0;
    double _cur_progress = 0;
    double _time_total = 0;
    double _time_fft = 0;
    double _timer = 0;
    int _fft_count = 0;
};

class Logging
{
public:
    static const int LEVEL_DEBUG = 0;
    static const int LEVEL_INFO = 1;
    static const int LEVEL_WARNING = 2;
    static const int LEVEL_ERROR = 3;

public:
    Logging(int level = LEVEL_INFO) : _level(level) { }

    void debug(const char* message...);
    void info(const char* message...);
    void warning(const char* message...);
    void error(const char* message...);
    void result(const char* message...);

    void report(const std::string& message);
    void report_progress();
    void report_factor(InputNum& input, const arithmetic::Giant& f);
    void report_result(const std::string& message);

    int level() { return _level; }
    Progress& progress() { return _progress; }
    const std::string& prefix() { return _prefix; }
    void set_prefix(const std::string& prefix) { _prefix = prefix; }

private:
    int _level;
    Progress _progress;
    std::string _prefix;
    bool _print_prefix = true;
};
