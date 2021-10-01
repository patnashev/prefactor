#pragma once

#include <vector>

class ProbSmooth
{
public:
    ProbSmooth();

    double F(double a);
    double rho(double u) { return F(1/u); }
    double smooth(double log_number, double log_largest_divisor) { return F(log_largest_divisor/log_number); }
    double factoring(double log_B1, double log_B2, double log_sieving_depth, double log_known_divisors);
    double factoring_fixed(double log_B1, double log_B2, double log_factor, double log_known_divisors);

private:
    std::vector<double> _F;
};
