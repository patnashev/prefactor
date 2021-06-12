#pragma once

#include "arithmetic.h"
#include "group.h"
#include "edwards.h"
#include "primelist.h"
#include "inputnum.h"

class Stage2
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
        std::vector<int> distances;
    };

public:
    Stage2(arithmetic::GWArithmetic& gw, PrimeList& primes, int B1, int B2, int D, int A, int L) : _gw(gw), _primes(primes), _B1(B1), _B2(B2), _D(D), _A(A), _L(L)
    {
        _pairing = get_pairing(primes, B1, B2, D, A, L, true);
    }

    static Pairing get_pairing(PrimeList& primes, int B1, int B2, int D, int A, int L, bool with_distances);
    template<class Element>
    int precompute(arithmetic::DifferentialGroupArithmetic<Element>& arithmetic, Element& X1, Element& XD, Element& XDA, std::vector<std::unique_ptr<Element>>& precomp);

    arithmetic::GWArithmetic& gw() { return _gw; }
    PrimeList& primes() { return _primes; }

protected:
    int _B1;
    int _B2;
    int _D;
    int _A;
    int _L;
    Pairing _pairing;

private:
    arithmetic::GWArithmetic& _gw;
    PrimeList& _primes;
};

class PP1Stage2 : public Stage2
{
public:
    PP1Stage2(arithmetic::GWArithmetic& gw, PrimeList& primes, int B1, int B2, int D, int A, int L) : Stage2(gw, primes, B1, B2, D, A, L) { }

    void run(InputNum& input, arithmetic::Giant& P, bool minus1);
};

class EdECMStage2 : public Stage2
{
public:
    EdECMStage2(arithmetic::GWArithmetic& gw, PrimeList& primes, int B1, int B2, int D, int L, int LN) : Stage2(gw, primes, B1, B2, D, 1, L), _LN(LN) { }

    void run(InputNum& input, arithmetic::GWNum& ed_d, arithmetic::EdPoint& P);

private:
    int _LN;
};
