#pragma once

#include <memory>
#include <vector>
#include "arithmetic.h"
#include "edwards.h"
#include "inputnum.h"
#include "logging.h"
#include "file.h"

#define FERMAT_APPID 2

int fermat_main(int argc, char *argv[]);

class Fermat
{
public:
    Fermat(int exponent, arithmetic::GWState& gwstate, Logging& logging) : _exponent(exponent), _input(1, 2, exponent, 1), _gwstate(gwstate), _logging(logging), _gw(gwstate), _ed(_gw)
    {
        int j;
        for (j = 1, _exponent_c = 0; j < exponent; j <<= 1, _exponent_c++);
    }

    bool read_points(File& file);
    bool read_state(File& file, uint64_t B1);
    std::string verify(bool verify_curve);
    void modulus(int curve, File& file_result);
    void stage1(uint64_t B1, File& file_state, File& file_result);

    std::vector<std::unique_ptr<arithmetic::EdPoint>>& points() { return _points; }
    uint64_t B0() { return _B0; }

private:
    void write_file(File& file, uint64_t B1, std::vector<std::unique_ptr<arithmetic::EdPoint>>& points);
    bool read_file(File& file, int& seed, uint64_t& B0, std::vector<std::unique_ptr<arithmetic::EdPoint>>& points);

private:
    int _exponent;
    char _exponent_c;
    InputNum _input;
    arithmetic::GWState& _gwstate;
    Logging& _logging;
    arithmetic::GWArithmetic _gw;
    arithmetic::EdwardsArithmetic _ed;

    int _seed;
    uint64_t _B0;
    std::vector<std::unique_ptr<arithmetic::EdPoint>> _points;
    std::vector<std::unique_ptr<arithmetic::EdPoint>> _state;
};
