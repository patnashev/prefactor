#pragma once

#include <vector>
#include "prob.h"

class Params
{
protected:
    Params() { }
    Params(int B1_, int B2_) : B1(B1_), B2(B2_) { }
    virtual ~Params() { }

    double get_profit(int B1, int B2, int max_size, ProbSmooth& prob, double log_sieving_depth, double log_known_divisors, double primality_cost);
    void find_pp1_optimal_bounds(int max_size, ProbSmooth& prob, double log_sieving_depth, double log_known_divisors, double primality_cost);
    void find_pp1_optimal_B2(int max_size, ProbSmooth& prob, double log_sieving_depth, double log_known_divisors, double primality_cost);
    void find_pp1_stage2_optimal_params(int max_size);

public:
    virtual int stage1_cost();
    virtual int stage1_size();
    virtual int stage2_cost();
    virtual int stage2_size();

public:
    int B1;
    int B2;
    int D;
    int A;
    int L;
    double pairing;
};

class PM1Params : public Params
{
public:
    PM1Params(int B1, int B2, int max_size);
    PM1Params(int B1, int max_size, ProbSmooth& prob, double log_sieving_depth, double log_known_divisors, double primality_cost);
    PM1Params(int max_size, ProbSmooth& prob, double log_sieving_depth, double log_known_divisors, double primality_cost);

};

class PP1Params : public Params
{
public:
    PP1Params(int B1, int B2, int max_size);
    PP1Params(int B1, int max_size, ProbSmooth& prob, double log_sieving_depth, double log_known_divisors, double primality_cost);
    PP1Params(int max_size, ProbSmooth& prob, double log_sieving_depth, double log_known_divisors, double primality_cost);

    virtual int stage1_cost() override;
};

class EdECMParams : public Params
{
public:
    EdECMParams(int B1, int B2, int max_size);

    virtual int stage1_cost() override;
    virtual int stage1_size() override;
    virtual int stage2_cost() override;
    virtual int stage2_size() override;

public:
    int W;
    int LN;
};
