
#include "lucas.h"
#include "edwards.h"

// Stage 1 exponent with prime powers <= B1
Giant get_stage1_exp(PrimeList primes, int B0, int B1)
{
    int j, k;

    PrimeIterator it = primes.begin();
    for (; *it < B0; it++);

    int sqrtB1 = (int)sqrt(B1);
    Giant tmp(GiantsArithmetic::default_arithmetic(), (int)(B1/0.69/32) + 10);
    Giant tmp2(GiantsArithmetic::default_arithmetic(), tmp.capacity() < 8192 ? tmp.capacity() : 8192);
    tmp2 = 1;
    while (*it <= B1)
    {
        // Building exponent with prime powers <= B1
        j = *it;
        if (*it <= sqrtB1)
        {
            k = B1/(*it);
            while (j <= k)
                j *= *it;
        }
        tmp2 *= j;
        it++;
        if (*it > B1 || tmp2.size() > 8190)
        {
            if (tmp == 0)
                tmp = tmp2;
            else
                tmp *= tmp2;
            tmp2 = 1;
        }
    }

    return tmp;
}

// P-1 factoring stage 1.
int do_minus1stage1(PrimeList primes, InputNum& input, GWArithmetic& gw, int B1, Giant& res)
{
    int i, j;
    int len;
    int retval = TRUE;

    GWNum X(gw);

    printf("%s, P-1 stage 1, B1 = %d%s\n", input.display_text().data(), B1, input.gfn() != 0 ? ", GFN" : "");
    double timer = getHighResTimer();
    int transforms = -(int)gw.gwdata()->fft_count;

    Giant tmp = get_stage1_exp(primes, 3, B1);
    gw.setmulbyconst(3);
    X = 3;
    len = tmp.bitlen() - 1;
    gwset_carefully_count(gw.gwdata(), 30);
    for (i = 1; i <= len; i++)
        gw.square(X, X, (tmp.bit(len - i) ? GWMUL_MULBYCONST : 0) | (i < len ? GWMUL_STARTNEXTFFT : 0));
    for (j = 2; j <= B1; j <<= 1)
        gw.carefully().square(X, X, 0);
    for (j = 0; j < input.gfn(); j++)
        gw.carefully().square(X, X, 0);

    tmp = X;
    res = tmp + inv(tmp, gw.N());
    tmp = gcd(std::move(tmp) - 1, gw.N());

    timer = (getHighResTimer() - timer)/getHighResTimerFrequency();
    transforms += (int)gw.gwdata()->fft_count;
    if (tmp == 0 || tmp == gw.N())
    {
        printf("All divisors of N < B1.\n");
        retval = FALSE;
    }
    else if (tmp != 1)
    {
        report_factor(tmp, input);
        retval = FALSE;
    }
    else
    {
        printf("No factors found, transforms: %d, time: %d s.\n", transforms, (int)timer);
    }

    return retval;
}

// P+1 factoring stage 1.
int do_plus1stage1(PrimeList primes, InputNum& input, GWArithmetic& gw, int B0, int B1, Giant& P, std::string& sP, Giant& res)
{
    int j, k;
    int retval = TRUE;
    Giant tmp;

    printf("%s, P+1 stage 1, B1 = %d, P = %s%s.\n", input.display_text().data(), B1, sP.data(), input.gfn() != 0 ? ", GFN" : "");
    double timer = getHighResTimer();
    int transforms = -(int)gw.gwdata()->fft_count;

    LucasVArithmetic lucas(gw.carefully());
    LucasV V(lucas);
    // V_tmp with careful start
    V.V() = P;
    if (B0 <= 2)
    {
        for (j = 2; j <= B1; j <<= 1)
            lucas.dbl(V, V);
        for (j = 0; j < input.gfn(); j++)
            lucas.dbl(V, V);
    }
    if (B0 < 3)
        B0 = 3;

    PrimeIterator it = primes.begin();
    for (; *it < B0; it++);

    lucas.set_gw(gw);
    int sqrtB1 = (int)sqrt(B1);
    while (*it <= B1)
    {
        j = *it;
        lucas.mul(V, *it, it.pos(), V);
        if (*it <= sqrtB1)
        {
            k = B1/(*it);
            while (j <= k)
            {
                lucas.mul(V, *it, it.pos(), V);
                j *= *it;
            }
        }
        it++;
    }

    tmp = V.V();
    res = tmp;
    tmp = gcd(std::move(tmp) - 2, gw.N());

    timer = (getHighResTimer() - timer)/getHighResTimerFrequency();
    transforms += (int)gw.gwdata()->fft_count;
    if (tmp == 0 || tmp == gw.N())
    {
        printf("All divisors of N < B1.\n");
        retval = FALSE;
    }
    else if (tmp != 1)
    {
        report_factor(tmp, input);
        retval = FALSE;
    }
    else
    {
        printf("No factors found, transforms: %d, time: %d s.\n", transforms, (int)timer);
    }

    return retval;
}

// EdECM factoring stage 1 using signed window.
int do_edecm_stage1(PrimeList primes, InputNum& input, GWArithmetic& gw, int B1, int W, EdPoint& P)
{
    int j;
    int retval = TRUE;
    Giant tmp;
    std::vector<int16_t> naf_w;

    if (W > 16)
        W = 16;
    tmp = get_stage1_exp(primes, 3, B1);
    get_NAF_W(W, tmp, naf_w);

    printf("%s, EdECM stage 1, B1 = %d, W = %d.\n", input.display_text().data(), B1, W);
    double timer = getHighResTimer();
    int transforms = -(int)gw.gwdata()->fft_count;

    EdwardsArithmetic& ed = P.arithmetic();
    ed.set_gw(gw.carefully());
    for (j = 2; j <= B1; j <<= 1)
        ed.dbl(P, P, ed.ED_PROJECTIVE);
    ed.dbl(P, P, 0);
    ed.set_gw(gw);
    ed.mul(P, W, naf_w, P);

    tmp = *P.X;
    tmp = gcd(std::move(tmp), gw.N());

    timer = (getHighResTimer() - timer)/getHighResTimerFrequency();
    transforms += (int)gw.gwdata()->fft_count;
    if (tmp == 0 || tmp == gw.N())
    {
        printf("All divisors of N < B1.\n");
        retval = FALSE;
    }
    else if (tmp != 1)
    {
        report_factor(tmp, input);
        retval = FALSE;
    }
    else
    {
        printf("No factors found, transforms: %d, time: %d s.\n", transforms, (int)timer);
    }

    return retval;
}
