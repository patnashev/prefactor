#include "prob.h"
#include "math.h"

#define PRECISION 50
#define F_SIZE 1000

// Eric Bach and Rene Peralta
// Asymptotic semismoothness probabilities
// https://doi.org/10.1090/S0025-5718-96-00775-2
// Section 4.
// rho(100) == 1e-229
ProbSmooth::ProbSmooth()
{
    int i, j, k, fi;
    double p, t;

    _F.resize(F_SIZE + 1);
    _F[F_SIZE] = 1.0;
    fi = F_SIZE - 1;
    std::vector<double> c_prev(PRECISION);
    std::vector<double> c_cur(PRECISION);
    for (k = 2; k <= 125 && fi > 0; k++)
    {
        if (k == 2)
        {
            for (c_cur[0] = 1 - log(2), i = 1; i < c_cur.size(); i++)
                c_cur[i] = 1.0/i/pow(2.0, i);
        }
        else
        {
            for (i = 1; i < c_cur.size(); i++)
                for (c_cur[i] = 0, p = k*i, j = i - 1; j >= 0; j--, p *= k)
                    c_cur[i] += c_prev[j]/p;
            for (c_cur[0] = 0, j = 1; j < c_cur.size(); j++)
                c_cur[0] += c_cur[j]/(j + 1);
            c_cur[0] /= k - 1;
        }
        for (; fi >= (double)F_SIZE/k; fi--)
            for (_F[fi] = 0, p = 1, t = k - (double)F_SIZE/fi, i = 0; i < c_cur.size(); i++, p *= t)
                _F[fi] += c_cur[i]*p;
        c_cur.swap(c_prev);
    }
    for (; fi > 0; fi--)
        _F[fi] = 0;
}

double ProbSmooth::F(double a)
{
    if (a >= 1.0)
        return 1.0;
    int i = (int)(a*F_SIZE);
    double t = a*F_SIZE - i;
    //return _F[i]*(1 - t) + _F[i + 1]*t;
    return exp(log(_F[i])*(1 - t) + log(_F[i + 1])*t);
}

double ProbSmooth::factoring(double log_B1, double log_B2, double log_sieving_depth, double log_known_divisors)
{
    int i, j;
    double delta = 0.1;
    double sum = 0;
    for (i = 0; i < 1000; i++)
    {
        double l = log_sieving_depth + i*delta; // log(factor)
        double d = log_sieving_depth*(1/l - 1/(l + delta)); // probability of one factor in the range
        double ds = d; // probability of up to 5 factors in the range
        for (j = 0; j < 5; j++)
            ds = ds*d + d;
        double a = log_B1/(l - log_known_divisors); // B1 = (factor/knownDivisors)^a
        double b = log_B2/(l - log_known_divisors); // B2 = (factor/knownDivisors)^b
        double stage1 = F(a); // probability of B1-smooth factor
        double stage2 = 0; // integrating stage 2 probabilities
        for (j = (int)(a*F_SIZE); j < (int)(b*F_SIZE); j++)
            stage2 += F(a*F_SIZE/(double)(F_SIZE - j))/j;
        sum += (stage1 + stage2)*ds; // probability of successful factorization
    }

    return sum;
}

double ProbSmooth::factoring_fixed(double log_B1, double log_B2, double log_factor, double log_known_divisors)
{
    int j;
    double l = log_factor; // log(factor)
    double a = log_B1/(log_factor - log_known_divisors); // B1 = (factor/knownDivisors)^a
    double b = log_B2/(log_factor - log_known_divisors); // B2 = (factor/knownDivisors)^b
    double stage1 = F(a); // probability of B1-smooth factor
    double stage2 = 0; // integrating stage 2 probabilities
    for (j = (int)(a*F_SIZE); j < (int)(b*F_SIZE); j++)
        stage2 += F(a*F_SIZE/(double)(F_SIZE - j))/j;
    return stage1 + stage2; // probability of successful factorization
}
