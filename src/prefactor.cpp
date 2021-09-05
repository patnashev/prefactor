#define PREFACTOR_VERSION "0.8.0"

#include <ctype.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <math.h>
#include <array>
#include <iostream>
#include <fstream>
#include <unordered_set>

#include "gwnum.h"
#include "cpuid.h"
#include "arithmetic.h"
#include "exception.h"
#include "inputnum.h"
#include "primelist.h"
#include "stage1.h"
#include "stage2.h"
#include "file.h"
#include "logging.h"
#include "prob.h"
#ifdef FACTORING
#include "factoring.h"
#endif
#ifdef NETPF
int net_main(int argc, char *argv[]);
#endif

using namespace arithmetic;

#include "precompute.h"

#include "prob.c"

void sigterm_handler(int signo)
{
    Task::abort();
    printf("Terminating...\n");
    signal(signo, sigterm_handler);
}

int main(int argc, char *argv[])
{
    signal(SIGTERM, sigterm_handler);
    signal(SIGINT, sigterm_handler);

    int i, j;
    GWState gwstate;
    GWArithmetic gw(gwstate);
    int B1 = 0;
    int B2 = 0;
    int minus1 = 0;
    int plus1 = 0;
    int edecm = 0;
    double sievingDepth = 0;
    int maxMem = 2048;
    int D, A, L;
    std::string sP;
    int curveType = 0;
    int curveSeed = 0;
    std::string curveX;
    std::string curveY;
    int K = 0;
    InputNum input;
    std::string toFile;
    int log_level = Logging::LEVEL_INFO;
    bool success = false;

    for (i = 1; i < argc; i++)
        if (argv[i][0] == '-' && argv[i][1])
        {
            switch (argv[i][1])
            {
            case 't':
                if (argv[i][2] && isdigit(argv[i][2]))
                    gwstate.thread_count = atoi(argv[i] + 2);
                else if (!argv[i][2] && i < argc - 1)
                {
                    i++;
                    gwstate.thread_count = atoi(argv[i]);
                }
                else
                    break;
                if (gwstate.thread_count == 0 || gwstate.thread_count > 64)
                    gwstate.thread_count = 1;
                continue;

            case 'q':
                if (argv[i][2] != '\"' && !isdigit(argv[i][2]))
                    break;
                if (!input.parse(argv[i] + 2))
                {
                    printf("Invalid number format.\n");
                    return 1;
                }
                continue;
            }

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
                    j = (int)strlen(argv[i]);
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
                    sievingDepth = log2(sievingDepth);
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
#ifdef FACTORING
            else if (strcmp(argv[i], "-factoring") == 0)
                return factoring_main(argc, argv);
#endif
#ifdef NETPF
            else if (strcmp(argv[i], "-net") == 0)
                return net_main(argc, argv);
#endif
            else if (i < argc - 1 && strcmp(argv[i], "-file") == 0)
            {
                i++;
                toFile = argv[i];
            }
            else if (strcmp(argv[i], "-v") == 0)
            {
                printf("Prefactor version " PREFACTOR_VERSION ", Gwnum library version " GWNUM_VERSION "\n");
                return 0;
            }
        }
        else
        {
            if (!input.parse(argv[i]))
            {
                File file(argv[i], 0);
                if (!input.read(file))
                {
                    printf("File %s is missing or corrupted.\n", argv[i]);
                    return 1;
                }
            }
        }
    if (input.empty())
    {
        printf("Usage: prefactor {-B1 10000 -B2 100000 | -S sievingDepth [-B1 10000] [-B2 100000]} [-minus1] [-plus1] [-edecm] options {\"K*B^N+C\" | file}\n");
        printf("Options: [-M maxMemory] [-t Threads] [-P 2/7] [-curve {curve2x8 | curve12 | random | seed 123 | xy 17/19 17/33}]\n");
        return 0;
    }
    if (!toFile.empty())
    {
        File file(toFile, 0);
        input.write(file);
    }

    if (!minus1 && !plus1 && !edecm)
        minus1 = 1;

    Logging logging(log_level);
    logging.progress().time_init(0);
    input.setup(gwstate);
    logging.info("Using %s.\n", gwstate.fft_description.data());

    double primalityCost = 0;
    if (input.b() == 2)
        primalityCost = input.n();
    else
        primalityCost = gwstate.N->bitlen()*1.2;

    ProbSmooth prob;
    double knownDivisors = 1;
    if (minus1 && input.gfn() > 0)
        knownDivisors = input.gfn() + 1;
    if (!minus1 && plus1 && !edecm && (sP.empty() || sP == "2/7"))
        knownDivisors = log2(6);
    if (!minus1 && plus1 && !edecm && (sP == "6/5"))
        knownDivisors = log2(4);
    if (!minus1 && !plus1 && edecm)
        knownDivisors = log2(31); // by Yves

    int maxSize = (int)(maxMem/(gwnum_size(gwstate.gwdata())/1048576.0));

    if (sievingDepth != 0 && B1 == 0)
    {
        get_optimal_bounds(&B1, &B2, maxSize, sievingDepth, knownDivisors, primalityCost, minus1 ? 0 : plus1 ? 1 : 5);
        logging.info("Optimal B1 = %d, B2 = %d.\n", B1, B2);
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

    PrimeList primes(B2 + 100);

    if (edecm)
        get_edecm_stage1_params(B1, maxSize, &K);
    double pairing;
    get_stage2_params(B1, B2, maxSize, &D, &A, &L, &pairing);
    int size = 0;
    if (minus1)
    {
        logging.progress().add_stage(get_stage1_cost(B1, 0));
        logging.progress().add_stage(get_stage2_cost(B1, B2, D, A, L, pairing));
        size = get_stage2_size(D, A, L);
    }
    if (plus1)
    {
        logging.progress().add_stage(get_stage1_cost(B1, 1));
        logging.progress().add_stage(get_stage2_cost(B1, B2, D, A, L, pairing));
        size = get_stage2_size(D, A, L);
    }
    if (edecm)
    {
        logging.progress().add_stage(get_stage1_cost(B1, K));
        logging.progress().add_stage(get_stage2_cost(B1, B2, D, A, L, pairing)); // TODO:
        int size_ed1 = get_edecm_stage1_size(K);
        if (size < size_ed1)
            size = size_ed1;
    }
    int cost = logging.progress().cost_total();
    if (2*cost > primalityCost)
        logging.info("Running at %.0f%% cost of a primality test, using %.0f MB.\n", cost/primalityCost*100, gwnum_size(gwstate.gwdata())/1048576.0*size);
    else
        logging.info("Running at 1/%.0f cost of a primality test, using %.0f MB.\n", primalityCost/cost, gwnum_size(gwstate.gwdata())/1048576.0*size);
    if (sievingDepth != 0)
    {
        double value = prob.factoring(log2(B1), log2(B2), sievingDepth, knownDivisors);
        logging.info("Probability of a factor 1/%.0f, overall speedup %.2f%%.\n", 1/value, 100*(value - cost/primalityCost));
    }

    try
    {
        if (minus1 && !success)
        {
            File file1(std::to_string(gwstate.fingerprint) + ".m1", gwstate.fingerprint);
            File file12(std::to_string(gwstate.fingerprint) + ".m12", gwstate.fingerprint);
            File file2(std::to_string(gwstate.fingerprint) + ".m2", gwstate.fingerprint);
            PP1Stage1::State* interstate = read_state<PP1Stage1::State>(&file12);
            if (interstate == nullptr)
            {
                PM1Stage1 stage1(primes, B1);
                stage1.init(&input, &gwstate, &file1, &logging);
                stage1.run();
                success = stage1.success();
                if (!success && B2 > B1)
                {
                    interstate = new PP1Stage1::State();
                    interstate->V() = std::move(stage1.V());
                    file12.write(*interstate);
                }
            }
            logging.progress().next_stage();
            if (interstate != nullptr)
            {
                PP1Stage2 stage2(logging, primes, B1, B2, D, A, L);
                stage2.init(&input, &gwstate, &file2, &logging, interstate->V(), true);
                stage2.run();
                success = stage2.success();
            }
            logging.progress().next_stage();
            file1.clear();
            file12.clear();
            file2.clear();
        }
        if (plus1 && !success)
        {
            if (sP.empty())
                //sP = "6/5";
                sP = "2/7";
            File file1(std::to_string(gwstate.fingerprint) + ".p1", gwstate.fingerprint);
            File file12(std::to_string(gwstate.fingerprint) + ".p12", gwstate.fingerprint);
            File file2(std::to_string(gwstate.fingerprint) + ".p2", gwstate.fingerprint);
            PP1Stage1::State* interstate = read_state<PP1Stage1::State>(&file12);
            if (interstate == nullptr)
            {
                PP1Stage1 stage1(primes, B1, sP);
                stage1.init(&input, &gwstate, &file1, &logging);
                stage1.run();
                success = stage1.success();
                if (!success && B2 > B1)
                {
                    interstate = new PP1Stage1::State(std::move(*stage1.state()));
                    file12.write(*interstate);
                }
            }
            logging.progress().next_stage();
            if (interstate != nullptr)
            {
                PP1Stage2 stage2(logging, primes, B1, B2, D, A, L);
                stage2.init(&input, &gwstate, &file2, &logging, interstate->V(), false);
                stage2.run();
                success = stage2.success();
            }
            logging.progress().next_stage();
            file1.clear();
            file12.clear();
            file2.clear();
        }
        if (edecm && !success)
        {
            Giant X, Y, Z, T, EdD;
            {
                EdwardsArithmetic ed(gw.carefully());
                EdPoint P(ed);
                GWNum ed_d(gw.carefully());

                if (curveType == 0)
                    P = ed.from_small(17, 19, 17, 33, &ed_d);
                else if (curveType == 1)
                    P = ed.from_small(5, 23, -1, 7, &ed_d);
                else if (curveType == 2)
                {
                    try
                    {
                        P = ed.gen_curve(curveSeed, &ed_d);
                    }
                    catch (const NoInverseException& e)
                    {
                        logging.set_prefix(input.display_text() + ", EdECM, ");
                        logging.report_factor(input, e.divisor);
                        return 0;
                    }
                    catch (const ArithmeticException&)
                    {
                        printf("Invalid curve.\n");
                        return 1;
                    }
                }
                else if (curveType == 3)
                {
                    int xa = stoi(curveX);
                    int xb = 1;
                    for (i = curveX[0] == '-' ? 1 : 0; isdigit(curveX[i]); i++);
                    if (curveX[i] == '/')
                        xb = stoi(curveX.substr(i + 1));
                    int ya = stoi(curveY);
                    int yb = 1;
                    for (i = curveY[0] == '-' ? 1 : 0; isdigit(curveY[i]); i++);
                    if (curveY[i] == '/')
                        yb = stoi(curveY.substr(i + 1));
                    P = ed.from_small(xa, xb, ya, yb, &ed_d);
                }
                Giant tmp;
                tmp = ed.jinvariant(ed_d);
                if (tmp.size() > 1)
                    logging.info("Curve j-invariant RES64: %08X%08X\n", tmp.data()[1], tmp.data()[0]);
                else if (tmp.size() > 0)
                    logging.info("Curve j-invariant RES64: %08X%08X\n", 0, tmp.data()[0]);
                else
                    logging.info("Curve j-invariant RES64: %08X%08X\n", 0, 0);
                P.serialize(X, Y, Z, T);
                EdD = ed_d;
            };

            File file1(std::to_string(gwstate.fingerprint) + ".ed1", gwstate.fingerprint);
            File file12(std::to_string(gwstate.fingerprint) + ".ed12", gwstate.fingerprint);
            File file2(std::to_string(gwstate.fingerprint) + ".ed2", gwstate.fingerprint);
            EdECMStage1::State* interstate = read_state<EdECMStage1::State>(&file12);
            if (interstate == nullptr)
            {
                EdECMStage1 stage1(primes, B1, K);
                stage1.init(&input, &gwstate, &file1, &logging, X, Y, Z, T, EdD);
                stage1.run();
                success = stage1.success();
                if (!success && B2 > B1)
                {
                    interstate = new EdECMStage1::State(std::move(*stage1.state()));
                    file12.write(*interstate);
                }
            }
            logging.progress().next_stage();
            if (interstate != nullptr)
            {
                EdECMStage2 stage2(logging, primes, B1, B2, 210, 5, 20);
                stage2.init(&input, &gwstate, &file2, &logging, interstate->X(), interstate->Y(), interstate->Z(), interstate->T(), EdD);
                stage2.run();
                success = stage2.success();
            }
            logging.progress().next_stage();
            file1.clear();
            file12.clear();
            file2.clear();
        }
        if (!success)
        {
            logging.info("%s, no factors found.\n", input.input_text().data(), logging.progress().time_total());
            logging.result("%s, no factors found, time: %.1f s.\n", input.input_text().data(), logging.progress().time_total());
        }
    }
    catch (const TaskAbortException&)
    {
    }

    gwstate.done();

    return 0;
}
