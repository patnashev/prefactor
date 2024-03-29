#pragma once

#include <string>
#include <vector>
#include <cstdint>
#include "prob.h"

class Params
{
protected:
    Params() { }
    Params(uint64_t B1_, uint64_t B2_) : B1(B1_), B2(B2_) { }
    virtual ~Params() { }

    double get_profit(uint64_t B1, uint64_t B2, int max_size, ProbSmooth& prob, double log_sieving_depth, double log_known_divisors, double primality_cost);
    void find_pp1_optimal_bounds(int max_size, ProbSmooth& prob, double log_sieving_depth, double log_known_divisors, double primality_cost);
    void find_pp1_optimal_B2(int max_size, ProbSmooth& prob, double log_sieving_depth, double log_known_divisors, double primality_cost);
    void find_pp1_stage2_optimal_params(int max_size);
    void find_poly_stage2_optimal_params(int max_size, int threads, int mem_model);

public:
    virtual double stage1_cost();
    virtual int stage1_size();
    virtual double stage2_cost();
    virtual int stage2_size();
    virtual std::string unique_id();

public:
    uint64_t B1;
    uint64_t B2;
    int D;
    int A;
    int L;
    int PolyPower;

protected:
    double pairing;
    int PolyMemModel;
};

class PM1Params : public Params
{
public:
    PM1Params(uint64_t B1, uint64_t B2, int max_size, bool poly, int threads, int mem_model);
    PM1Params(uint64_t B1, int max_size, ProbSmooth& prob, double log_sieving_depth, double log_known_divisors, double primality_cost);
    PM1Params(int max_size, ProbSmooth& prob, double log_sieving_depth, double log_known_divisors, double primality_cost);

};

class PP1Params : public Params
{
public:
    PP1Params(uint64_t B1, uint64_t B2, int max_size, bool poly, int threads, int mem_model);
    PP1Params(uint64_t B1, int max_size, ProbSmooth& prob, double log_sieving_depth, double log_known_divisors, double primality_cost);
    PP1Params(int max_size, ProbSmooth& prob, double log_sieving_depth, double log_known_divisors, double primality_cost);

    virtual double stage1_cost() override;
};

class EdECMParams : public Params
{
public:
    EdECMParams(uint64_t B1, uint64_t B2, int max_size, bool poly, int threads, int mem_model);

    virtual double stage1_cost() override;
    virtual int stage1_size() override;
    virtual double stage2_cost() override;
    virtual int stage2_size() override;

public:
    int W;
    int LN;
};
