#define PREFACTOR_VERSION "0.6.0"

#include <ctype.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "gwnum/gwnum.h"
#include "gwnum/cpuid.h"

// The number being tested
giant N = NULL;
// N mod 3417905339
uint32_t fingerprint;
// Bitlen of base
double log2_b = 1;
// Human-readable representation
char Nstr[100];
// Gwnum handle
gwhandle* gwdata = NULL;
// Returns giant for operations modulo N
#define getg() popg(&gwdata->gdata, ((int)gwdata->bit_length >> 5) + 10)
#define freeg() pushg(&gwdata->gdata, 1)
// Normalization required after each subtract
int must_norm = 1;
// Max number of gwnums
int maxSize = 0;
// All primes for the test
int *primes = NULL;
// Edwards curve d parameter
gwnum EdD = NULL;

#include "md5.c"
#include "utils.c"
#include "file.c"
#include "lucas.c"
#include "edwards.c"
#include "stage1.c"
#include "stage2.c"
#include "prob.c"

// Initializes Gwnum in N=k*b^n+c mode. Computes N.
int initKBNC(int threadCount, double k, unsigned int b, unsigned int n, int c)
{
    if (b != 2)
        log2_b = log(b)/log(2.0);

    N = allocgiant(((int)(log2(b)*n + log2(k)) >> 5) + 10);
    ultog(b, N);
    power(N, n);
    dblmulg(k, N);
    sladdg(c, N);
    fingerprint = gmodul(N, 3417905339);

    gwdata = (gwhandle*)malloc(sizeof(gwhandle));
    gwinit(gwdata);
    //gwset_safety_margin(gwdata, 0.5);
    //gwset_larger_fftlen_count(gwdata, 1);
    //gwset_maxmulbyconst(gwdata, 19*33);
    gwset_num_threads(gwdata, threadCount);
    if (gwsetup(gwdata, k, b, n, c))
    {
        printf("gwsetup failed.\n");
        return FALSE;
    }
    must_norm = gwdata->EXTRA_BITS < 2.0;

    if ((int)k != 1)
        snprintf(Nstr, 100, "%d*%d^%d%c%d", (int)k, b, n, c < 0 ? '-' : '+', abs(c));
    else
        snprintf(Nstr, 100, "%d^%d%c%d", b, n, c < 0 ? '-' : '+', abs(c));

    printf("Starting factoring of %s\n", Nstr);
    char buf[200];
    gwfft_description(gwdata, buf);
    printf("Using %s\n", buf);

    return TRUE;
}

// Initializes Gwnum in generic mode. Computes N.
int initKBNCstrings(int threadCount, unsigned int k, unsigned int b, unsigned int n, int c, char *sk, char *sb)
{
    giant gk;
    if (k == 0)
    {
        gk = allocgiant(((int)(strlen(sk)*3.33) >> 5) + 10);
        ctog(sk, gk);
    }
    else
    {
        gk = allocgiant(1);
        ultog(k, gk);
    }
    giant gb;
    if (b == 0)
    {
        gb = allocgiant(((int)(strlen(sb)*3.33) >> 5) + 10);
        ctog(sb, gb);
        log2_b = bitlen(gb);
    }
    else
    {
        gb = allocgiant(1);
        ultog(b, gb);
        log2_b = log(b)/log(2.0);
    }

    N = allocgiant(((bitlen(gb)*n + bitlen(gk)) >> 5) + 10);
    gtog(gb, N);
    power(N, n);
    mulg(gk, N);
    sladdg(c, N);
    fingerprint = gmodul(N, 3417905339);

    gwdata = (gwhandle*)malloc(sizeof(gwhandle));
    gwinit(gwdata);
    //gwset_safety_margin(gwdata, 0.5);
    //gwset_larger_fftlen_count(gwdata, 1);
    gwset_num_threads(gwdata, threadCount);
    if (gwsetup_general_mod_giant(gwdata, N))
    {
        printf("gwsetup failed.\n");
        return FALSE;
    }
    must_norm = gwdata->EXTRA_BITS < 2.0;

    if (k != 1)
        sprintf(Nstr, "%d*", k);
    if (k == 0)
        snprintf(Nstr, 50, "%s*", sk);
    int len = strlen(Nstr);
    if (b != 0)
        snprintf(Nstr + len, 50 - len, "%d^", b);
    else
        snprintf(Nstr + len, 50 - len, "%s^", sb);
    len += strlen(Nstr + len);
    if (len + snprintf(Nstr + len, 50 - len, "%d%c%d", n, c < 0 ? '-' : '+', abs(c) >= 50))
        strcat(Nstr + len, "...");

    printf("Starting factoring of %s\n", Nstr);
    char buf[200];
    gwfft_description(gwdata, buf);
    printf("Using %s\n", buf);

    return TRUE;
}

// Initializes Gwnum in generic mode with generic N.
int initGeneric(int threadCount)
{
    fingerprint = gmodul(N, 3417905339);

    gwdata = (gwhandle*)malloc(sizeof(gwhandle));
    gwinit(gwdata);
    //gwset_safety_margin(gwdata, 0.5);
    //gwset_larger_fftlen_count(gwdata, 1);
    gwset_num_threads(gwdata, threadCount);
    if (gwsetup_general_mod_giant(gwdata, N))
    {
        printf("gwsetup failed.\n");
        return FALSE;
    }
    must_norm = gwdata->EXTRA_BITS < 2.0;

    snprintf(Nstr, 100, "%d-bit number", bitlen(N));
    printf("Starting factoring of %s\n", Nstr);
    char buf[200];
    gwfft_description(gwdata, buf);
    printf("Using %s\n", buf);

    return TRUE;
}

// Releases memory
void cleanup()
{
    if (N)
        free(N);
    N = NULL;
    if (gwdata)
    {
        gwdone(gwdata);
        free(gwdata);
    }
    gwdata = NULL;
}

// Tries to parse the string as "k*b^n+c".
int parseKBNC(const char *s, int *type, int *num, char *sk, char *sb, int max_input)
{
    int j;
    char buf[10];
    int ret_type = 1;

    j = 0;
    while (s[j] && isdigit(s[j]))
        j++;
    if (!j || j >= max_input)
        return FALSE;
    if (s[j] == '*')
    {
        strncpy(sk, s, j);
        sk[j] = 0;
        if (j < 10)
            num[0] = atoi(sk);
        else
        {
            num[0] = 0;
            ret_type = 2;
        }
        s += j + 1;
        j = 0;
        while (s[j] && isdigit(s[j]))
            j++;
        if (!j || j >= max_input)
            return FALSE;
    }
    else
        num[0] = 1;
    if (s[j] == '^')
    {
        strncpy(sb, s, j);
        sb[j] = 0;
        if (j < 10)
            num[1] = atoi(sb);
        else
        {
            num[1] = 0;
            ret_type = 2;
        }
    }
    else
        return FALSE;
    s += j + 1;
    j = 0;
    while (s[j] && isdigit(s[j]))
        j++;
    if (!j || j > 9)
        return FALSE;
    strncpy(buf, s, j);
    buf[j] = 0;
    num[2] = atoi(buf);
    if (s[j] == '+')
        num[3] = 1;
    else if (s[j] == '-')
        num[3] = -1;
    else
        return FALSE;
    s += j + 1;
    j = 0;
    while (s[j] && isdigit(s[j]))
        j++;
    if (!j || j > 9)
        return FALSE;
    strncpy(buf, s, j);
    buf[j] = 0;
    num[3] *= atoi(buf);

    *type = ret_type;
    return TRUE;
}

int parseFraction(char *s, giant res)
{
    int j;
    char buf[10];
    int a, b = 1;

    j = 0;
    while (s[j] && isdigit(s[j]))
        j++;
    if (!j || j >= 10)
        return FALSE;
    strncpy(buf, s, j);
    buf[j] = 0;
    a = atoi(buf);
    if (s[j] == '/')
    {
        s += j + 1;
        j = 0;
        while (s[j] && isdigit(s[j]))
            j++;
        if (!j || j >= 10)
            return FALSE;
        strncpy(buf, s, j);
        buf[j] = 0;
        b = atoi(buf);
    }
    itog(b, res);
    invg(N, res);
    if (res->sign < 0)
        return FALSE;
    imulg(a, res);
    modg(N, res);

    return TRUE;
}

int main(int argc, char	*argv[])
{
    int i, j;
    char *s;

    int ThreadCount = 1;
    int B1 = 0;
    int B2 = 0;
    int minus1 = 0;
    int plus1 = 0;
    int edecm = 0;
    int gfn = 0;
    double sievingDepth = 0;
    int maxMem = 2048;
    int D, A, L;
    char *sP = NULL;
    int curveType = 0;
    int curveSeed = 0;
    char *curveX = NULL;
    char *curveY = NULL;
    int K = 0;
    int InputType = 0;
    int InputKBNC[4] = {0};
    #define MAX_INPUT 10000
    char sk[MAX_INPUT];
    char sb[MAX_INPUT];

    for (i = 1; i < argc; i++)
        if (argv[i][0] == '-' && argv[i][1])
        {
            switch (argv[i][1])
            {
            case 't':
                if (argv[i][2])
                    ThreadCount = atoi(argv[i] + 2);
                else if (i < argc - 1)
                {
                    i++;
                    ThreadCount = atoi(argv[i]);
                }
                if (ThreadCount == 0 || ThreadCount > 64)
                    ThreadCount = 1;
                break;

            case 'q':
                s = argv[i] + 2;
                if (*s == '"')
                    s++;
                if (!parseKBNC(s, &InputType, InputKBNC, sk, sb, MAX_INPUT))
                {
                    printf("Invalid number format.\n");
                    return 1;
                }
                break;

            default:
                if (i < argc - 1 && strcmp(argv[i], "-B1") == 0)
                {
                    i++;
                    B1 = atoi(argv[i]);
                    if (B1 < 100)
                        B1 = 100;
                }
                else if (i < argc - 1 && strcmp(argv[i], "-B2") == 0)
                {
                    i++;
                    B2 = atoi(argv[i]);
                }
                else if (i < argc - 1 && strcmp(argv[i], "-S") == 0)
                {
                    i++;
                    if (sscanf(argv[i], "%lf", &sievingDepth))
                    {
                        j = strlen(argv[i]);
                        if (argv[i][j - 1] == 'P')
                            sievingDepth *= 1e15;
                        if (argv[i][j - 1] == 'T')
                            sievingDepth *= 1e12;
                        if (argv[i][j - 1] == 'G')
                            sievingDepth *= 1e9;
                        if (argv[i][j - 1] == 'M')
                            sievingDepth *= 1e6;
                    }
                    if (sievingDepth > 100)
                        sievingDepth = log(sievingDepth)/log(2);
                }
                else if (i < argc - 1 && strcmp(argv[i], "-M") == 0)
                {
                    i++;
                    maxMem = atoi(argv[i]);
                }
                else if (i < argc - 1 && strcmp(argv[i], "-P") == 0)
                {
                    i++;
                    sP = argv[i];
                    plus1 = 1;
                }
                else if (i < argc - 1 && strcmp(argv[i], "-curve") == 0)
                {
                    i++;
                    if (strcmp(argv[i], "curve2x8") == 0)
                        curveType = 0;
                    else if (strcmp(argv[i], "curve12") == 0)
                        curveType = 1;
                    else if (strcmp(argv[i], "seed") == 0)
                    {
                        i++;
                        curveType = 2;
                        curveSeed = atoi(argv[i]);
                    }
                    else if (strcmp(argv[i], "random") == 0)
                    {
                        curveType = 2;
                        double timer = getHighResTimer();
                        curveSeed = *(int *)&timer;
                    }
                    else if (strcmp(argv[i], "xy") == 0)
                    {
                        i++;
                        curveType = 3;
                        curveX = argv[i];
                        i++;
                        curveY = argv[i];
                    }

                    edecm = 1;
                }
                else if (strcmp(argv[i], "-minus") == 0 || strcmp(argv[i], "-minus1") == 0)
                    minus1 = 1;
                else if (strcmp(argv[i], "-plus") == 0 || strcmp(argv[i], "-plus1") == 0)
                    plus1 = 1;
                else if (strcmp(argv[i], "-ecm") == 0 || strcmp(argv[i], "-eecm") == 0 || strcmp(argv[i], "-edecm") == 0)
                    edecm = 1;
                else if (strcmp(argv[i], "-v") == 0)
                {
                    printf("Prefactor version "PREFACTOR_VERSION", Gwnum library version "GWNUM_VERSION"\n");
                    return 0;
                }
                break;
            }
        }
        else
        {
            s = argv[i];
            if (*s == '"')
                s++;
            if (!parseKBNC(s, &InputType, InputKBNC, sk, sb, MAX_INPUT))
            {
                if (!readFromFile(argv[i], 0, &j, NULL, NULL))
                {
                    printf("File %s is missing or corrupted.\n", argv[i]);
                    return 1;
                }
                InputType = 3;
            }
        }
    if (InputType == 0)
    {
        printf("Usage: prefactor {-B1 10000 -B2 100000 | -S sievingDepth [-B1 10000] [-B2 100000]} [-minus1] [-plus1] [-edecm] options {\"K*B^N+C\" | file}\n");
        printf("Options: [-M maxMemory] [-t Threads] [-P 2/7] [-curve {curve2x8 | curve12 | random | seed 123 | xy 17/19 17/33}]\n");
        return 0;
    }

    if ((InputType == 1 || InputType == 2) && InputKBNC[0] == 1 && InputKBNC[3] == 1 && (InputKBNC[2] & (InputKBNC[2] - 1)) == 0)
        gfn = InputKBNC[2];
    if (!minus1 && !plus1 && !edecm)
        minus1 = 1;

    // KBNC test
    if (InputType == 1)
        if (!initKBNC(ThreadCount, InputKBNC[0], InputKBNC[1], InputKBNC[2], InputKBNC[3]))
            return 1;
    if (InputType == 2)
        if (!initKBNCstrings(ThreadCount, InputKBNC[0], InputKBNC[1], InputKBNC[2], InputKBNC[3], sk, sb))
            return 1;
    if (InputType == 3)
        if (!initGeneric(ThreadCount))
            return 1;

    double primalityCost = 0;
    if (InputType == 1 || InputType == 2)
    {
        if (InputKBNC[1] == 2)
            primalityCost = InputKBNC[2];
        else
            primalityCost = InputKBNC[2]*log2_b*1.2;
    }
    if (InputType == 3)
        primalityCost = bitlen(N)*1.2;

    double knownDivisors = 0;
    if (minus1 && gfn > 0)
        knownDivisors = log(gfn)/log(2);
    if (!minus1 && plus1 && !edecm && (sP == NULL || strcmp(sP, "2/7") == 0))
        knownDivisors = log(3)/log(2);
    if (!minus1 && plus1 && !edecm && (sP != NULL && strcmp(sP, "6/5") == 0))
        knownDivisors = 1;
    if (!minus1 && !plus1 && edecm)
        knownDivisors = log(8)/log(2);

    int maxSize = (int)(maxMem/(gwnum_size(gwdata)/1048576.0));

    if (sievingDepth != 0 && B1 == 0)
    {
        get_optimal_bounds(&B1, &B2, maxSize, sievingDepth, knownDivisors, primalityCost, minus1 ? 0 : plus1 ? 1 : 5);
        printf("Optimal B1 = %d, B2 = %d.\n", B1, B2);
    }
    if (sievingDepth != 0 && B1 != 0 && B2 == 0)
    {
        get_optimal_B2(B1, &B2, maxSize, sievingDepth, knownDivisors, primalityCost, minus1 ? 0 : plus1 ? 1 : 5);
    }
    if (B1 == 0)
    {
        printf("B1 parameter is missing and no sievingDepth set to calculate it.\n");
        return 1;
    }
    if (B2 < B1)
        B2 = B1;

    if (edecm)
        get_edecm_stage1_params(B1, maxSize, &K);
    double pairing;
    get_stage2_params(B1, B2, maxSize, &D, &A, &L, &pairing);
    int cost = get_stage1_cost(B1, minus1 ? 0 : plus1 ? 1 : K) + get_stage2_cost(B1, B2, D, A, L, pairing);
    int size = get_stage2_size(D, A, L);
    if (edecm && size < get_edecm_stage1_size(K))
        size = get_edecm_stage1_size(K);
    printf("Running at 1/%.0f cost of a primality test, using %.0f MB.\n", primalityCost/cost, gwnum_size(gwdata)/1048576.0*size);
    if (sievingDepth != 0)
    {
        double value = prob_smooth(B1, B2, sievingDepth, knownDivisors);
        printf("Probability of a factor 1/%.0f, overall speedup %.2f%%.\n", 1/value, 100*(value - cost/primalityCost));
    }
    costInit(cost);

    sieve(B2 + 100);

    giant P = getg();
    if (minus1)
    {
        if (do_minus1stage1(B1, gfn, P) && B2 > B1)
            do_pm1stage2(B1, B2, P, 1, D, A, L);
    }
    if (plus1)
    {
        if (sP == NULL)
            //sP = "6/5";
            sP = "2/7";
        if (!parseFraction(sP, P))
        {
            printf("Invalid P+1 parameter.\n");
            return 1;
        }
        if (do_plus1stage1(0, B1, P, sP, P) && B2 > B1)
            do_pm1stage2(B1, B2, P, 0, D, A, L);
    }
    if (edecm)
    {
        ed_point EdP = ed_alloc();
        if (curveType == 0)
            ed_from_small(17, 19, 17, 33, EdP);
        else if (curveType == 1)
            ed_from_small(5, 23, -1, 7, EdP);
        else if (curveType == 2)
        {
            if (!gen_ed_curve(curveSeed, EdP))
            {
                printf("Invalid curve.\n");
                return 1;
            }
        }
        else if (curveType == 3)
        {
            int xa = atoi(curveX);
            int xb = 1;
            for (i = curveX[0] == '-' ? 1 : 0; isdigit(curveX[i]); i++);
            if (curveX[i] == '/')
                xb = atoi(curveX + i + 1);
            int ya = atoi(curveY);
            int yb = 1;
            for (i = curveY[0] == '-' ? 1 : 0; isdigit(curveY[i]); i++);
            if (curveY[i] == '/')
                yb = atoi(curveY + i + 1);
            ed_from_small(xa, xb, ya, yb, EdP);
        }
        get_j_invariant(P);
        if (abs(P->sign) > 1)
            printf("Curve j-invariant RES64: %08X%08X\n", P->n[1], P->n[0]);
        else if (abs(P->sign) > 0)
            printf("Curve j-invariant RES64: %08X%08X\n", 0, P->n[0]);
        else
            printf("Curve j-invariant RES64: %08X%08X\n", 0, 0);
        if (do_edecm_stage1(B1, K, EdP) && B2 > B1)
            do_edecm_stage2(B1, B2, EdP, 210, 5, 2);
        ed_free(EdP);
    }
    freeg();

    cleanup();

    return 0;
}
