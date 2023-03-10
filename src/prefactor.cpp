
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
#include "config.h"
#include "inputnum.h"
#include "stage1.h"
#include "stage2.h"
#include "file.h"
#include "logging.h"
#include "prob.h"
#include "params.h"
#include "poly.h"
#include "stage2poly.h"
#include "version.h"

#ifdef FACTORING
int factoring_main(int argc, char *argv[]);
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
    setvbuf(stdout, NULL, _IONBF, 0);
#if defined(_MSC_VER) && !defined(_DEBUG)
    _set_error_mode(_OUT_TO_STDERR);
    _set_abort_behavior(0, _CALL_REPORTFAULT);
#endif

    File::FILE_APPID = 1;

    int i;
    GWState gwstate;
    GWArithmetic gw(gwstate);
    uint64_t B1 = 0;
    uint64_t B2 = 0;
    int minus1 = 0;
    int plus1 = 0;
    int edecm = 0;
    bool poly = false;
    int polyThreads = 1;
    int polyMemModel = 0;
    bool polyCheck = true;
    bool polyWrite = false;
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
    std::string log_file;
    bool success = false;

    auto known_factor_handler = [&](const char* factor) {
        if (!isdigit(factor[0]))
            return false;
        if (gwstate.known_factors.empty())
            gwstate.known_factors = factor;
        else
        {
            Giant f;
            f = factor;
            gwstate.known_factors *= f;
        }
        return true;
    };

    Config cnfg;
    cnfg.value_number("-t", 0, gwstate.thread_count, 1, 256)
        .value_number("-t", ' ', gwstate.thread_count, 1, 256)
        .value_number("-spin", ' ', gwstate.spin_threads, 1, 256)
        .value_enum("-cpu", ' ', gwstate.instructions, Enum<std::string>().add("SSE2", "SSE2").add("AVX", "AVX").add("FMA3", "FMA3").add("AVX512F", "AVX512F"))
        .value_number("-fft", '+', gwstate.next_fft_count, 0, 5)
        .group("-fft")
            .value_number("+", 0, gwstate.next_fft_count, 0, 5)
            .value_number("safety", ' ', gwstate.safety_margin, -10.0, 10.0)
            .check("generic", gwstate.force_general_mod, true)
            .end()
        .value_number("-M", ' ', maxMem)
        .value_number("-L3", ' ', PolyMult::L3_CACHE_MB, 0, INT_MAX)
        .value_code("-f", 0, known_factor_handler)
        .value_code("-f", ' ', known_factor_handler)
        .value_number("-B1", ' ', B1)
        .value_number("-B2", ' ', B2)
        .value_number("-S", ' ', sievingDepth, 0.0, 1e100)
            .on_code([&] { if (sievingDepth > 100) sievingDepth = log2(sievingDepth); })
        .check("-minus", minus1, 1)
        .check("-minus1", minus1, 1)
        .check("-plus", plus1, 1)
        .check("-plus1", plus1, 1)
        .value_string("-P", ' ', sP)
            .on_check(plus1, 1)
        .check("-ecm", edecm, 1)
        .check("-eecm", edecm, 1)
        .check("-edecm", edecm, 1)
        .group("-curve")
            .exclusive()
                .ex_case().check("curve2x8", curveType, 0).end()
                .ex_case().check("curve12", curveType, 1).end()
                .ex_case()
                    .value_number("seed", ' ', curveSeed, 1, INT_MAX)
                        .on_check(curveType, 2)
                    .end()
                .ex_case()
                    .check("random", curveType, 2)
                        .on_code([&] {
                                double timer = getHighResTimer();
                                curveSeed = abs(*(int *)&timer);
                            })
                    .end()
                .ex_case()
                    .list("xy", ' ', ' ')
                        .value_string(curveX)
                        .value_string(curveY)
                        .end()
                        .on_check(curveType, 3)
                    .end()
                .end()
                .on_check(edecm, 1)
            .end()
        .group("-poly")
            .value_number("t", 0, polyThreads, 1, 256)
            .value_number("t", ' ', polyThreads, 1, 256)
            .value_number("threads", ' ', polyThreads, 1, 256)
            .value_enum("mem", ' ', polyMemModel, Enum<int>().add("lowest", -2).add("low", -1).add("normal", 0).add("high", 1).add("highest", 2))
            .check("nocheck", polyCheck, false)
            .check("write", polyWrite, true)
            .end()
        .group("-time")
            .value_number("write", ' ', Task::DISK_WRITE_TIME, 1, INT_MAX)
            .value_number("progress", ' ', Task::PROGRESS_TIME, 1, INT_MAX)
            .end()
        .group("-log")
            .exclusive()
                .ex_case().check("debug", log_level, Logging::LEVEL_DEBUG).end()
                .ex_case().check("info", log_level, Logging::LEVEL_INFO).end()
                .ex_case().check("warning", log_level, Logging::LEVEL_WARNING).end()
                .ex_case().check("error", log_level, Logging::LEVEL_ERROR).end()
                .end()
            .value_string("file", ' ', log_file)
            .end()
        .check("-d", log_level, Logging::LEVEL_INFO)
#ifdef FACTORING
        .check_code("-factoring", [&] { exit(factoring_main(argc, argv)); })
#endif
#ifdef NETPF
        .check_code("-net", [&] { exit(net_main(argc, argv)); })
#endif
        .check_code("-v", [&] {
                print_banner();
                exit(0);
            })
        .value_code("-ini", ' ', [&](const char* param) {
                File ini_file(param, 0);
                ini_file.read_buffer();
                if (ini_file.buffer().empty())
                    printf("ini file not found: %s.\n", param);
                else
                    cnfg.parse_ini(ini_file);
                return true;
            })
        .value_string("-file", ' ', toFile)
        .value_code("-q", 0, [&](const char* param) {
                if (param[0] != '\"' && !isdigit(param[0]))
                    return false;
                if (!input.parse(param))
                    return false;
                return true;
            })
        .default_code([&](const char* param) {
                if (!input.parse(param))
                {
                    File file(param, 0);
                    if (!input.read(file))
                        printf("Unknown option %s.\n", param);
                }
            })
        .parse_args(argc, argv);

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
    if (edecm)
        gwstate.maxmulbyconst = 8;

    Logging logging(log_level);
    if (!log_file.empty())
        logging.file_log(log_file);
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
            params_pm1.reset(new PM1Params(B1, B2, maxSize, poly, polyThreads, polyMemModel));
        if (plus1)
            params_pp1.reset(new PP1Params(B1, B2, maxSize, poly, polyThreads, polyMemModel));
        if (edecm)
            params_edecm.reset(new EdECMParams(B1, B2, maxSize, poly, polyThreads, polyMemModel));
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
        File file_progress(std::to_string(gwstate.fingerprint), gwstate.fingerprint);
        file_progress.hash = false;
        logging.file_progress(&file_progress);

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
                    stage2.reset(new PP1Stage2Poly(params_pm1->B1, params_pm1->B2, params_pm1->D, params_pm1->PolyPower, polyThreads, polyMemModel, polyCheck));
                dynamic_cast<IPP1Stage2*>(stage2.get())->init(&input, &gwstate, &file2, &logging, interstate->V(), true);
                stage2->run();
                success = stage2->success();
            }
            logging.progress().next_stage();
            file1.clear();
            file12.clear();
            file2.clear(true);
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
                    stage2.reset(new PP1Stage2Poly(params_pp1->B1, params_pp1->B2, params_pp1->D, params_pp1->PolyPower, polyThreads, polyMemModel, polyCheck));
                dynamic_cast<IPP1Stage2*>(stage2.get())->init(&input, &gwstate, params_pp1->PolyPower == 0 || polyWrite ? &file2 : nullptr, &logging, interstate->V(), false);
                stage2->run();
                success = stage2->success();
            }
            logging.progress().next_stage();
            file1.clear();
            file12.clear();
            file2.clear(true);
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
                    stage2.reset(new EdECMStage2Poly(params_edecm->B1, params_edecm->B2, params_edecm->D, params_edecm->PolyPower, polyThreads, polyMemModel, polyCheck));
                dynamic_cast<IEdECMStage2*>(stage2.get())->init(&input, &gwstate, params_edecm->PolyPower == 0 || polyWrite ? &file2 : nullptr, &logging, interstate->X(), interstate->Y(), interstate->Z(), interstate->T(), EdD);
                stage2->run();
                success = stage2->success();
            }
            logging.progress().next_stage();
            file1.clear();
            file12.clear();
            file2.clear(true);
        }
        if (!success)
        {
            logging.result(false, "%s no factors found.\n", input.display_text().data());
            logging.result_save(input.input_text() + " no factors found, time: " + std::to_string((int)logging.progress().time_total()) + " s.\n");
        }
        file_progress.clear();
    }
    catch (const TaskAbortException&)
    {
    }

    gwstate.done();

    return 0;
}
