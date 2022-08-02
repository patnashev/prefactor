#define PREFACTOR_VERSION "0.9.0"

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
#include "params.h"
#include "poly.h"
#include "stage2poly.h"
#ifdef FACTORING
#include "factoring.h"
#endif
#ifdef NETPF
int net_main(int argc, char *argv[]);
#endif

using namespace arithmetic;

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
    uint64_t B1 = 0;
    uint64_t B2 = 0;
    int minus1 = 0;
    int plus1 = 0;
    int edecm = 0;
    bool poly = false;
    int polyThreads = 1;
    bool polyCheck = false;
    double sievingDepth = 0;
    uint64_t maxMem = 2048*1048576ULL;
    std::unique_ptr<PM1Params> params_pm1;
    std::unique_ptr<PP1Params> params_pp1;
    std::unique_ptr<EdECMParams> params_edecm;
    std::string sP;
    int curveType = 0;
    int curveSeed = 0;
    std::string curveX;
    std::string curveY;
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

            case 'f':
                if (argv[i][2] && isdigit(argv[i][2]))
                    gwstate.known_factors = argv[i] + 2;
                else if (!argv[i][2] && i < argc - 1)
                {
                    i++;
                    gwstate.known_factors = argv[i];
                }
                else
                    break;
                continue;
            }

            if (i < argc - 1 && strcmp(argv[i], "-B1") == 0)
            {
                i++;
                B1 = InputNum::parse_numeral(argv[i]);
                //if (B1 < 100)
                //    B1 = 100;
            }
            else if (i < argc - 1 && strcmp(argv[i], "-B2") == 0)
            {
                i++;
                B2 = InputNum::parse_numeral(argv[i]);
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
                maxMem = InputNum::parse_numeral(argv[i]);
            }
            else if (strncmp(argv[i], "-fft", 4) == 0 && ((!argv[i][4] && i < argc - 1) || argv[i][4] == '+'))
            {
                if (argv[i][4] == '+')
                    gwstate.next_fft_count = atoi(argv[i] + 5);
                else if (argv[i + 1][0] == '+')
                {
                    i++;
                    gwstate.next_fft_count = atoi(argv[i] + 1);
                }
            }
            else if (strcmp(argv[i], "-generic") == 0)
                gwstate.force_general_mod = true;
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
            else if (strcmp(argv[i], "-poly") == 0)
            {
                poly = true;
                while (true)
                    if (i < argc - 2 && strcmp(argv[i + 1], "threads") == 0)
                    {
                        i += 2;
                        polyThreads = atoi(argv[i]);
                    }
                    else if (i < argc - 1 && argv[i + 1][0] == 't')
                    {
                        i++;
                        polyThreads = atoi(argv[i] + 1);
                    }
                    else if (i < argc - 1 && strcmp(argv[i + 1], "check") == 0)
                    {
                        i++;
                        polyCheck = true;
                    }
                    else
                        break;
            }
            else if (strcmp(argv[i], "-time") == 0)
            {
                while (true)
                    if (i < argc - 2 && strcmp(argv[i + 1], "write") == 0)
                    {
                        i += 2;
                        Task::DISK_WRITE_TIME = atoi(argv[i]);
                    }
                    else if (i < argc - 2 && strcmp(argv[i + 1], "progress") == 0)
                    {
                        i += 2;
                        Task::PROGRESS_TIME = atoi(argv[i]);
                    }
                    else
                        break;
            }
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
            else if (i < argc - 1 && strcmp(argv[i], "-log") == 0)
            {
                i++;
                if (strcmp(argv[i], "debug") == 0)
                    log_level = Logging::LEVEL_DEBUG;
                if (strcmp(argv[i], "info") == 0)
                    log_level = Logging::LEVEL_INFO;
                if (strcmp(argv[i], "warning") == 0)
                    log_level = Logging::LEVEL_WARNING;
                if (strcmp(argv[i], "error") == 0)
                    log_level = Logging::LEVEL_ERROR;
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
        printf("Options: [-M maxMemory] [-t Threads] [-P 2/7] [-curve {curve2x8 | curve12 | random | seed 123 | xy 17/19 17/33}] [-poly [tThreads]] [-fft+1]\n");
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
    double knownDivisors_pm1 = 1;
    if (minus1 && input.gfn() > 0)
        knownDivisors_pm1 = input.gfn() + 1;
    double knownDivisors_pp1 = 1;
    if (plus1 && (sP.empty() || sP == "2/7"))
        knownDivisors_pp1 = log2(6);
    if (plus1 && (sP == "6/5"))
        knownDivisors_pp1 = log2(4);
    double knownDivisors_edecm = 1;
    if (edecm && (curveType == 0 || curveType == 2))
        knownDivisors_edecm = log2(31); // by Yves
    if (edecm && (curveType == 1))
        knownDivisors_edecm = log2(12);

    int maxSize = (int)(maxMem/(gwnum_size(gwstate.gwdata())));

    if (sievingDepth != 0 && B1 == 0)
    {
        if (minus1)
        {
            params_pm1.reset(new PM1Params(maxSize, prob, sievingDepth, knownDivisors_pm1, primalityCost));
            logging.info("Optimal P-1 B1 = %" PRId64 ", B2 = %" PRId64 ".\n", params_pm1->B1, params_pm1->B2);
        }
        if (plus1)
        {
            params_pp1.reset(new PP1Params(maxSize, prob, sievingDepth, knownDivisors_pp1, primalityCost));
            logging.info("Optimal P+1 B1 = %" PRId64 ", B2 = %" PRId64 ".\n", params_pp1->B1, params_pp1->B2);
        }
        if (edecm)
        {
            logging.error("EdECM needs explicit B1 and B2\n");
            return 1;
        }
    }
    else if (sievingDepth != 0 && B1 != 0 && B2 == 0)
    {
        if (minus1)
        {
            params_pm1.reset(new PM1Params(B1, maxSize, prob, sievingDepth, knownDivisors_pm1, primalityCost));
            logging.info("Optimal P-1 B2 = %" PRId64 ".\n", params_pm1->B2);
        }
        if (plus1)
        {
            params_pp1.reset(new PP1Params(B1, maxSize, prob, sievingDepth, knownDivisors_pp1, primalityCost));
            logging.info("Optimal P+1 B2 = %" PRId64 ".\n", params_pp1->B2);
        }
        if (edecm)
        {
            logging.error("EdECM needs explicit B1 and B2\n");
            return 1;
        }
    }
    else
    {
        if (B1 == 0)
        {
            printf("B1 parameter is missing and no sievingDepth set to calculate it.\n");
            return 1;
        }
        if (B2 < B1)
            B2 = B1;
        if (minus1)
            params_pm1.reset(new PM1Params(B1, B2, maxSize, poly, polyThreads));
        if (plus1)
            params_pp1.reset(new PP1Params(B1, B2, maxSize, poly, polyThreads));
        if (edecm)
            params_edecm.reset(new EdECMParams(B1, B2, maxSize, poly, polyThreads));
    }

    int size = 0;
    if (params_pm1)
    {
        logging.progress().add_stage(params_pm1->stage1_cost());
        logging.progress().add_stage(params_pm1->stage2_cost());
        size = params_pm1->stage2_size();
    }
    if (params_pp1)
    {
        logging.progress().add_stage(params_pp1->stage1_cost());
        logging.progress().add_stage(params_pp1->stage2_cost());
        int size_pp1 = params_pp1->stage2_size();
        if (size < size_pp1)
            size = size_pp1;
    }
    if (params_edecm)
    {
        logging.progress().add_stage(params_edecm->stage1_cost());
        logging.progress().add_stage(params_edecm->stage2_cost());
        int size_ed1 = params_edecm->stage1_size();
        if (size < size_ed1)
            size = size_ed1;
        int size_ed2 = params_edecm->stage2_size();
        if (size < size_ed2)
            size = size_ed2;
    }
    double cost = logging.progress().cost_total();
    if (2*cost > primalityCost)
        logging.info("Running at %.0f%% cost of a primality test, using %.0f MB.\n", cost/primalityCost*100, gwnum_size(gwstate.gwdata())/1048576.0*size);
    else
        logging.info("Running at 1/%.0f cost of a primality test, using %.0f MB.\n", primalityCost/cost, gwnum_size(gwstate.gwdata())/1048576.0*size);
    if (sievingDepth != 0)
    {
        double value = 0;
        if (params_pm1)
            value += prob.factoring(log2(params_pm1->B1), log2(params_pm1->B2), sievingDepth, knownDivisors_pm1);
        if (params_pp1)
            value += prob.factoring(log2(params_pp1->B1), log2(params_pp1->B2), sievingDepth, knownDivisors_pp1);
        if (params_edecm)
            value += prob.factoring(log2(params_edecm->B1), log2(params_edecm->B2), sievingDepth, knownDivisors_edecm);
        logging.info("Probability of a factor 1/%.0f, overall speedup %.2f%%.\n", 1/value, 100*(value - cost/primalityCost));
    }

    try
    {
        if (params_pm1 && !success)
        {
            uint32_t fingerprint = File::unique_fingerprint(gwstate.fingerprint, std::to_string(params_pm1->B1));
            File file1(std::to_string(fingerprint) + ".m1", fingerprint);
            File file12(std::to_string(fingerprint) + ".m12", fingerprint);
            File file2(std::to_string(fingerprint) + ".m2", File::unique_fingerprint(fingerprint, std::to_string(params_pm1->B2)));
            PP1Stage1::State* interstate = read_state<PP1Stage1::State>(&file12);
            if (interstate == nullptr)
            {
                PM1Stage1 stage1((int)params_pm1->B1);
                stage1.init(&input, &gwstate, &file1, &logging);
                stage1.run();
                success = stage1.success();
                if (!success && params_pm1->B2 > params_pm1->B1)
                {
                    interstate = new PP1Stage1::State();
                    interstate->V() = std::move(stage1.V());
                    file12.write(*interstate);
                }
            }
            logging.progress().next_stage();
            if (interstate != nullptr)
            {
                std::unique_ptr<Stage2> stage2;
                if (params_pm1->PolyPower == 0)
                    stage2.reset(new PP1Stage2(params_pm1->B1, params_pm1->B2, params_pm1->D, params_pm1->A, params_pm1->L, logging));
                else
                    stage2.reset(new PP1Stage2Poly(params_pm1->B1, params_pm1->B2, params_pm1->D, params_pm1->PolyPower, params_pm1->PolySmallPower, params_pm1->PolyThreads, polyCheck));
                dynamic_cast<IPP1Stage2*>(stage2.get())->init(&input, &gwstate, &file2, &logging, interstate->V(), true);
                stage2->run();
                success = stage2->success();
            }
            logging.progress().next_stage();
            file1.clear();
            file12.clear();
            file2.clear();
        }
        if (params_pp1 && !success)
        {
            if (sP.empty())
                //sP = "6/5";
                sP = "2/7";
            uint32_t fingerprint = File::unique_fingerprint(gwstate.fingerprint, sP + "." + std::to_string(params_pp1->B1));
            File file1(std::to_string(fingerprint) + ".p1", fingerprint);
            File file12(std::to_string(fingerprint) + ".p12", fingerprint);
            File file2(std::to_string(fingerprint) + ".p2", File::unique_fingerprint(fingerprint, std::to_string(params_pp1->B2)));
            PP1Stage1::State* interstate = read_state<PP1Stage1::State>(&file12);
            if (interstate == nullptr)
            {
                PP1Stage1 stage1((int)params_pp1->B1, sP);
                stage1.init(&input, &gwstate, &file1, &logging);
                stage1.run();
                success = stage1.success();
                if (!success && params_pp1->B2 > params_pp1->B1)
                {
                    interstate = new PP1Stage1::State(std::move(*stage1.state()));
                    file12.write(*interstate);
                }
            }
            logging.progress().next_stage();
            if (interstate != nullptr)
            {
                std::unique_ptr<Stage2> stage2;
                if (params_pp1->PolyPower == 0)
                    stage2.reset(new PP1Stage2(params_pp1->B1, params_pp1->B2, params_pp1->D, params_pp1->A, params_pp1->L, logging));
                else
                    stage2.reset(new PP1Stage2Poly(params_pp1->B1, params_pp1->B2, params_pp1->D, params_pp1->PolyPower, params_pp1->PolySmallPower, params_pp1->PolyThreads, polyCheck));
                dynamic_cast<IPP1Stage2*>(stage2.get())->init(&input, &gwstate, &file2, &logging, interstate->V(), false);
                stage2->run();
                success = stage2->success();
            }
            logging.progress().next_stage();
            file1.clear();
            file12.clear();
            file2.clear();
        }
        if (params_edecm && !success)
        {
            std::string jinvariant;
            Giant X, Y, Z, T, EdD;
            {
                EdwardsArithmetic ed(gw.carefully());
                EdPoint P(ed);
                GWNum ed_d(gw.carefully());

                try
                {
                    if (curveType == 0)
                        P = ed.from_small(17, 19, 17, 33, &ed_d);
                    else if (curveType == 1)
                        P = ed.from_small(5, 23, -1, 7, &ed_d);
                    else if (curveType == 2)
                    {
                        P = ed.gen_curve(curveSeed, &ed_d);
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
                    jinvariant = tmp.to_res64();
                    logging.info("Curve j-invariant RES64: %s\n", jinvariant.data());
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

                P.serialize(X, Y, Z, T);
                EdD = ed_d;
            };

            uint32_t fingerprint = File::unique_fingerprint(gwstate.fingerprint, jinvariant + "." + std::to_string(params_edecm->B1));
            File file1(std::to_string(fingerprint) + ".ed1", File::unique_fingerprint(fingerprint, std::to_string(params_edecm->W)));
            File file12(std::to_string(fingerprint) + ".ed12", fingerprint);
            File file2(std::to_string(fingerprint) + ".ed2", File::unique_fingerprint(fingerprint, std::to_string(params_edecm->B2)));
            EdECMStage1::State* interstate = read_state<EdECMStage1::State>(&file12);
            if (interstate == nullptr)
            {
                EdECMStage1 stage1((int)params_edecm->B1, params_edecm->W);
                stage1.init(&input, &gwstate, &file1, &logging, &X, &Y, &Z, &T, &EdD);
                stage1.run();
                success = stage1.success();
                if (!success && params_edecm->B2 > params_edecm->B1)
                {
                    interstate = new EdECMStage1::State(std::move(*stage1.state()));
                    file12.write(*interstate);
                }
            }
            logging.progress().next_stage();
            //for (i = 0; i < (3 << params_edecm->PolyPower); i++)
            if (interstate != nullptr)
            {
                //file2.clear();
                //params_edecm->B2 = params_edecm->B1 + i*params_edecm->D;
                std::unique_ptr<Stage2> stage2;
                if (params_edecm->PolyPower == 0)
                    stage2.reset(new EdECMStage2(params_edecm->B1, params_edecm->B2, params_edecm->D, params_edecm->L, params_edecm->LN, logging));
                else
                    stage2.reset(new EdECMStage2Poly(params_edecm->B1, params_edecm->B2, params_edecm->D, params_edecm->PolyPower, params_edecm->PolySmallPower, params_edecm->PolyThreads, polyCheck));
                dynamic_cast<IEdECMStage2*>(stage2.get())->init(&input, &gwstate, &file2, &logging, interstate->X(), interstate->Y(), interstate->Z(), interstate->T(), EdD);
                stage2->run();
                success = stage2->success();
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
