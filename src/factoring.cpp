#define PREFACTOR_FACTORING_VERSION "1.0.0"

#include <iostream>
#include <cmath>
#include <string.h>
#include <stdlib.h>
#include "gwnum.h"
#include "cpuid.h"
#include "factoring.h"
#include "exception.h"
#include "task.h"
#include "primelist.h"
#include "params.h"
#include "stage1.h"
#include "stage2.h"

using namespace arithmetic;

void write_dhash(const std::string& filename, const std::string& dhash)
{
    FILE* fp = fopen(filename.data(), "w");
    if (fp)
    {
        fwrite(dhash.data(), 1, dhash.length(), fp);
        fclose(fp);
    }
}

void Factoring::write_file(File& file, char type, uint64_t B1, std::vector<Curve>& points)
{
    std::unique_ptr<Writer> writer(file.get_writer());
    writer->write(File::MAGIC_NUM);
    writer->write(FACTORING_APPID + (0 << 8) + (type << 16) + (_file_version << 24));
    writer->write(_gwstate.fingerprint);
    writer->write(_seed);
    writer->write((int)points.size());
    writer->write(B1);
    for (auto it = points.begin(); it != points.end(); it++)
    {
        writer->write(it->X);
        writer->write(it->Y);
        writer->write(it->Z);
        writer->write(it->T);
    }
    file.commit_writer(*writer);
}

bool Factoring::read_file(File& file, char type, int& seed, uint64_t& B1, std::vector<Curve>& points)
{
    file.appid = FACTORING_APPID;
    std::unique_ptr<Reader> reader(file.get_reader());
    if (!reader)
    {
        file.appid = 2;
        reader.reset(file.get_reader());
        if (!reader)
            return false;
        _file_version = 0;
    }
    else
    {
        if (!reader || reader->type() != type)
            return false;
        _file_version = reader->version();
        uint32_t fingerprint;
        if (!reader->read(fingerprint) || fingerprint != _gwstate.fingerprint)
            return false;
    }
    if (!reader->read(seed))
        return false;
    int count;
    if (!reader->read(count))
        return false;
    if (!reader->read(B1))
        return false;
    points.resize(count);
    for (int i = 0; i < count; i++)
    {
        if (!reader->read(points[i].X) || !reader->read(points[i].Y) || !reader->read(points[i].Z) || !reader->read(points[i].T))
            return false;
    }
    return true;
}

void Factoring::write_points(File& file)
{
    write_file(file, 0, _B1, _points);
}

bool Factoring::read_points(File& file)
{
    int seed;
    uint64_t B1;
    std::vector<Curve> points;
    if (!read_file(file, 0, seed, B1, points))
        return false;
    _seed = seed;
    _B1 = B1;
    _points = std::move(points);
    _logging.info("%d curve%s starting with #%d, B1 = %" PRId64 ".\n", _points.size(), _points.size() > 1 ? "s" : "", _seed, _B1);
    return true;
}

bool Factoring::read_state(File& file, uint64_t B1)
{
    int seed;
    uint64_t b1;
    std::vector<Curve> points;
    points.reserve(_points.size());
    if (!read_file(file, 1, seed, b1, points))
        return false;
    if (_seed != seed)
        return false;
    if (B1 != b1)
        return false;
    _state = std::move(points);
    _logging.info("resuming from curve %d.\n", (int)_state.size());
    return true;
}

std::string Factoring::generate(int seed, int count)
{
    _seed = seed;
    _B1 = 1;
    _points.clear();
    _state.clear();
    _points.reserve(count);

    GWArithmetic& gw = _gw.carefully();
    Writer hasher;

    int i, len;
    Giant tmp;
    GWNum s(gw), t(gw);
    s = 12;
    t = 40;
    tmp = seed;
    len = tmp.bitlen() - 1;
    for (i = 1; i <= len; i++)
    {
        GWNum lambda = (3*square(s) - 8)/(2*t);
        GWNum snext = lambda*lambda - s - s;
        t = lambda*(std::move(s) - snext) - t;
        s = std::move(snext);
        if (tmp.bit(len - i))
        {
            lambda = (t - 40)/(s - 12);
            snext = lambda*lambda - 12 - s;
            t = lambda*(std::move(s) - snext) - t;
            s = std::move(snext);
        }
    }

    Giant gx, gy, gx8, gisdx8, gd;
    Giant n1, nx8, nisdx8;
    for (i = 0; i < count; i++)
    {
        // Asserting T^2 = S^3 - 8S - 32
        GWASSERT((gw.popg() = square(t))%gw.N() == (gw.popg() = (square(s) - 8)*s - 32)%gw.N());

        // SqrtD = (8A^2 - 1) * (8A^2 + 8A + 1) / (8A^2 + 4A + 1)^2
        // B = A * 2(4A + 1) / (8A^2 - 1)
        GWNum alpha = 1/((t + 25)/(s - 9) + 1);
        GWNum alpha8alpha(gw);
        gw.setmulbyconst(8);
        gw.square(alpha, alpha8alpha, GWMUL_MULBYCONST);
        alpha += alpha;
        GWNum sqrt_d = alpha8alpha - 1;
        alpha8alpha += alpha;
        GWNum beta = alpha8alpha/sqrt_d;
        alpha8alpha += 1;
        alpha8alpha += alpha;
        alpha += alpha;
        sqrt_d *= alpha8alpha + std::move(alpha);
        sqrt_d /= square(std::move(alpha8alpha));

        GWNum x8 = 2*beta - 1;
        GWNum isdx8 = inv(x8*sqrt_d);
        GWNum x = x8*(4*beta - 3)/(6*beta - 5);
        GWNum y = x8*(t*(t + 50) - 104 - square(s)*(2*s - 27))/((t - 2 + 3*s)*(t + 16 + s));

        (gx = x) %= gw.N();
        (gy = y) %= gw.N();
        (gx8 = x8) %= gw.N();
        (gisdx8 = isdx8) %= gw.N();
        (gd = square(sqrt_d)) %= gw.N();
        (n1 = gw.N()) -= 1;
        (nx8 = gw.N()) -= gx8;
        (nisdx8 = gw.N()) -= gisdx8;

        // Checking torsion points
        ArithmeticException e("Torsion point.");
        if (gx == 0 && gy == 1)
            throw e;
        if (gx == 0 && gy == n1)
            throw e;
        if (gx == 1 && gy == 0)
            throw e;
        if (gx == n1 && gy == 0)
            throw e;
        if (gx == gx8 && gy == gx8)
            throw e;
        if (gx == nx8 && gy == gx8)
            throw e;
        if (gx == gx8 && gy == nx8)
            throw e;
        if (gx == nx8 && gy == nx8)
            throw e;
        if (gx == gisdx8 && gy == gisdx8)
            throw e;
        if (gx == nisdx8 && gy == gisdx8)
            throw e;
        if (gx == gisdx8 && gy == nisdx8)
            throw e;
        if (gx == nisdx8 && gy == nisdx8)
            throw e;

        // Asserting X^2 + Y^2 = 1 + d * X^2 * Y^2
        GWASSERT((gw.popg() = square(x) + square(y))%gw.N() == (gw.popg() = 1 + square(sqrt_d)*square(x)*square(y))%gw.N());

        hasher.write((char*)gd.data(), gd.size()*4);

        _points.emplace_back();
        _points.back().X = gx;
        _points.back().Y = gy;
        _points.back().Z = 1;
        _points.back().T = 0;
        _points.back().D = gd;

        if (i + seed == 1)
        {
            GWNum lambda = (3*square(s) - 8)/(2*t);
            GWNum snext = lambda*lambda - s - s;
            t = lambda*(std::move(s) - snext) - t;
            s = std::move(snext);
        }
        else
        {
            GWNum lambda = (t - 40)/(s - 12);
            GWNum snext = lambda*lambda - 12 - s;
            t = lambda*(std::move(s) - snext) - t;
            s = std::move(snext);
        }

        if ((i + 1)%1000 == 0)
            printf("%d\n", i + 1);
    }

    return hasher.hash_str();
}

void Factoring::compute_d(bool stateless)
{
    int i;
    std::vector<std::unique_ptr<arithmetic::EdPoint>> ed_points(_points.size());
    for (i = stateless ? (int)_state.size() : 0; i < _points.size(); i++)
        if (_points[i].D.empty())
        {
            ed_points[i].reset(new EdPoint(_ed));
            ed_points[i]->deserialize(_points[i].X, _points[i].Y, _points[i].Z, _points[i].T);
            _ed.d_ratio(*ed_points[i], *ed_points[i]->X, *ed_points[i]->Y);
            ed_points[i]->Z.reset(ed_points[i]->Y.release());
        }
    _ed.normalize(ed_points.begin(), ed_points.end(), _ed.ED_PROJECTIVE);
    for (i = 0; i < _points.size(); i++)
        if (ed_points[i])
            (_points[i].D = *ed_points[i]->X) %= _gw.N();
}

void Factoring::normalize(std::vector<Curve>& points)
{
    int i;
    std::vector<std::unique_ptr<arithmetic::EdPoint>> ed_points(points.size());
    for (i = 0; i < points.size(); i++)
    {
        ed_points[i].reset(new EdPoint(_ed));
        ed_points[i]->deserialize(points[i].X, points[i].Y, points[i].Z, points[i].T);
    }
    _ed.normalize(ed_points.begin(), ed_points.end(), _ed.ED_PROJECTIVE);
    for (i = 0; i < points.size(); i++)
        ed_points[i]->serialize(points[i].X, points[i].Y, points[i].Z, points[i].T);
}

std::string Factoring::verify(bool verify_curve)
{
    int i;
    std::vector<std::unique_ptr<arithmetic::EdPoint>> ed_points(_points.size());
    for (i = 0; i < _points.size(); i++)
    {
        ed_points[i].reset(new EdPoint(_ed));
        ed_points[i]->deserialize(_points[i].X, _points[i].Y, _points[i].Z, _points[i].T);
        _ed.d_ratio(*ed_points[i], *ed_points[i]->X, *ed_points[i]->Y);
        ed_points[i]->Z.reset(ed_points[i]->Y.release());
    }
    try
    {
        _ed.normalize(ed_points.begin(), ed_points.end(), _ed.ED_PROJECTIVE);
    }
    catch (const NoInverseException& e)
    {
        if (!verify_curve)
            throw;
        _logging.report_factor(_input, e.divisor);
        for (i = 0; i < _points.size(); i++)
        {
            try
            {
                ed_points[i]->normalize();
            }
            catch (const NoInverseException&)
            {
                _logging.warning("curve #%d found the factor.\n", _seed + i);
                _ed.gen_curve(_seed + i, ed_points[i]->X.get());
                ed_points[i]->Z.reset();
            }
        }
    }

    Giant gd;
    Writer hasher;
    for (i = 0; i < _points.size(); i++)
    {
        (gd = *ed_points[i]->X) %= _gw.N();
        hasher.write((char*)gd.data(), gd.size()*4);
    }

    std::string dhash = hasher.hash_str();
    _logging.info("dhash: %s\n", dhash.data());

    return dhash;
}

void Factoring::modulus(int curve, File& file_result)
{
    int i, j;
    if (curve != 0 && curve - _seed >= 0 && curve - _seed < 1024)
    {
        i = curve - _seed;
        j = i + 1;
        _logging.info("applying modulus to curve #%d (offset %d).\n", curve, i);
    }
    else if (curve != 0)
    {
        _logging.error("curve #%d not found.\n", curve);
        return;
    }
    else
    {
        i = 0;
        j = (int)_points.size();
        _logging.info("applying modulus to all curves.\n");
    }
    for (; i < j; i++)
    {
        _points[i].X %= _gw.N();
        _points[i].Y %= _gw.N();
    }
    write_file(file_result, 0, _B1, _points);
}

bool Factoring::split(int offset, int count, Factoring& result)
{
    if (offset < 0 || offset >= _points.size())
    {
        _logging.error("offset %d outside boundaries.\n", offset);
        return false;
    }
    if (count <= 0 || offset + count > _points.size())
    {
        _logging.error("can't find %d curves at offset %d.\n", count, offset);
        return false;
    }
    _logging.info("splitting %d curve%s starting with #%d.\n", count, count > 1 ? "s" : "", _seed + offset);

    result._seed = _seed + offset;
    result._B1 = _B1;
    result._state.clear();
    if (&result != this)
    {
        result._points.resize(count);
        for (int i = 0; i < count; i++)
            result._points[i] = _points[offset + i];
    }
    else
    {
        result._points.erase(result._points.begin(), result._points.begin() + offset);
        result._points.resize(count);
    }
    return true;
}

bool Factoring::merge(Factoring& other)
{
    if (other._B1 != _B1)
    {
        _logging.error("B1 mismatch.\n");
        return false;
    }
    if (other._seed != _seed + _points.size())
    {
        _logging.error("first curve of merging file must follow the last curve of source file.\n");
        return false;
    }
    _logging.info("merging %d curve%s starting with #%d.\n", other._points.size(), other._points.size() > 1 ? "s" : "", other._seed);

    _points.reserve(_points.size() + other._points.size());
    for (auto it = other._points.begin(); it != other._points.end(); it++)
        _points.push_back(std::move(*it));
    return true;
}

void Factoring::copy(Factoring& result)
{
    if (&result == this)
        return;

    result._seed = _seed;
    result._B1 = _B1;
    result._state.clear();
    result._points.resize(_points.size());
    for (size_t i = 0; i < _points.size(); i++)
        result._points[i] = _points[i];
}

void Factoring::stage1(uint64_t B1next, uint64_t B1max, uint64_t maxMem, File& file_state, File& file_result)
{
    int i;

    int maxSize = (int)(maxMem/(gwnum_size(_gwstate.gwdata())));
    EdECMStage1 stage1(_B1, B1next, B1max, maxSize);
    stage1.set_no_check_success();

    std::string prefix = _logging.prefix();
    Logging* logging = &_logging;
    SubLogging sub_logging(_logging, _logging.level() + 1);
    for (i = 0; i < _points.size(); i++)
        _logging.progress().add_stage(stage1.exp_len());
    for (i = 0; i < _state.size(); i++)
        _logging.progress().next_stage();
    if (_points.size() > 1)
    {
        _logging.info("%d bits, W = %d\n", stage1.exp_len(), stage1.W());
        _logging.progress().update(0, 0);
        sub_logging.progress().add_stage(stage1.exp_len());
        logging = &sub_logging;
    }

    compute_d(true);
    _state.reserve(_points.size());
    time_t last_write = time(NULL) - 270;
    while (_state.size() < _points.size())
    {
        i = (int)_state.size();
        if (_points.size() > 1)
            _logging.set_prefix(prefix + "#" + std::to_string(_seed + i) + ", ");
        File* file_stage1 = file_state.add_child(std::to_string(B1next) + "." + std::to_string(stage1.W()) + "." + std::to_string(i));
        //double timer = getHighResTimer();
        stage1.init(&_input, &_gwstate, file_stage1, logging, &_points[i].X, &_points[i].Y, &_points[i].Z, &_points[i].T, &_points[i].D);
        if (_points.size() == 1)
            _logging.set_prefix(prefix);
        if (stage1.state() != nullptr)
            last_write = 0;
        stage1.run();
        //timer = (getHighResTimer() - timer)/getHighResTimerFrequency();
        //_logging.info("%.1f%% done, %.3f ms per kilobit.\n", i/10.24, 1000000*timer/len);
        file_stage1->clear();
        _logging.set_prefix(prefix);

        _state.emplace_back();
        _state[i].X = std::move(stage1.state()->X());
        _state[i].Y = std::move(stage1.state()->Y());
        _state[i].Z = std::move(stage1.state()->Z());
        _state[i].T = std::move(stage1.state()->T());

        if (_points.size() > 1)
            _logging.progress().update(1, stage1.exp_len()/1000);
        if ((time(NULL) - last_write > 300 || Task::abort_flag()) && _state.size() < _points.size())
        {
            _logging.report_progress();
            write_file(file_state, 1, B1next, _state);
            last_write = time(NULL);
        }
        _logging.progress().next_stage();
        if (Task::abort_flag())
            throw TaskAbortException();
    }
    if (_points.size() > 1)
        _logging.report_progress();

    normalize(_state);
    write_file(file_result, 0, B1next, _state);
    _B1 = B1next;
    _points = std::move(_state);
}

void Factoring::stage2(uint64_t B2, uint64_t maxMem, bool poly, int threads, File& file_state)
{
    int i;

    i = 0;
    file_state.appid = FACTORING_APPID;
    std::unique_ptr<Reader> reader(file_state.get_reader());
    if (reader && reader->type() == 2)
        reader->read(i);
    reader.reset();

    int maxSize = (int)(maxMem/(gwnum_size(_gwstate.gwdata())));
    EdECMParams params_edecm(_B1, B2, maxSize, poly, threads);
    if (params_edecm.PolyPower > 0 && _gwstate.max_polymult_output() < 2*(1 << params_edecm.PolyPower))
    {
        _logging.error("FFT size too small for polynomial stage 2. Set -fft+1.\n");
        return;
    }
    EdECMStage2 stage2(params_edecm.B1, params_edecm.B2);
    if (params_edecm.PolyPower == 0)
        stage2.stage2_pairing(params_edecm.D, params_edecm.L, params_edecm.LN, _logging);
    else
        stage2.stage2_poly(params_edecm.D, params_edecm.L, params_edecm.LN, params_edecm.PolyDegree, params_edecm.PolyPower, params_edecm.PolyThreads);
    _logging.info("stage 2, B2 = %" PRId64 ", D = %d, degree %d.\n", B2, params_edecm.D, params_edecm.PolyDegree);

    SubLogging logging(_logging, _logging.level() + 1);
    logging.progress().add_stage((int)((params_edecm.B2 - params_edecm.B1)/params_edecm.D));
    _logging.progress().update(i/(double)_points.size(), i);
    std::string prefix = _logging.prefix();

    compute_d(true);
    time_t last_write = time(NULL) - 270;
    while (i < _points.size())
    {
        _logging.set_prefix(prefix + "#" + std::to_string(_seed + i) + ", ");
        //double timer = getHighResTimer();
        stage2.init(&_input, &_gwstate, nullptr, &logging, _points[i].X, _points[i].Y, _points[i].Z, _points[i].T, _points[i].D);
        stage2.run();
        //if (stage2.success())
        //    _logging.warning("curve #%d found the factor.\n", _seed + i);
        //timer = (getHighResTimer() - timer)/getHighResTimerFrequency();
        //_logging.info("%.1f%% done, %.3f ms per kilobit.\n", i/10.24, 1000000*timer/len);
        _logging.set_prefix(prefix);

        i++;
        if ((time(NULL) - last_write > 300 || Task::abort_flag()) && i < _points.size())
        {
            _logging.progress().update(i/(double)_points.size(), i);
            _logging.report_progress();
            std::unique_ptr<Writer> writer(file_state.get_writer());
            writer->write(File::MAGIC_NUM);
            writer->write(FACTORING_APPID + (0 << 8) + (2 << 16) + (0 << 24));
            writer->write(i);
            file_state.commit_writer(*writer);
            last_write = time(NULL);
        }
        if (Task::abort_flag())
            throw TaskAbortException();
    }
}

int factoring_main(int argc, char *argv[])
{
    int i;
    GWState gwstate;
    GWArithmetic gw(gwstate);
    EdwardsArithmetic ed(gw);
    uint64_t B1next = 0;
    uint64_t B1max = 0;
    uint64_t B2 = 0;
    uint64_t maxMem = 2048*1048576ULL;
    bool poly = false;
    int polyThreads = 1;
    int seed = 0;
    int count = 0;
    std::string filename;
    int generate = 0;
    bool range = false;
    int rangeOffset = 0;
    int rangeCount = 0;
    int verify = 0;
    int verifyCurve = 0;
    int modulus = 0;
    int modCurve = 0;
    int split = 0;
    int splitOffset = 0;
    int splitCount = 0;
    std::string splitName;
    int merge = 0;
    std::string mergeName;
    int log_level = Logging::LEVEL_INFO;
    InputNum input;
    Giant tmp;

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
                std::string sB1 = argv[i];
                size_t sep = sB1.find("/");
                if (sep != std::string::npos)
                {
                    B1next = InputNum::parse_numeral(sB1.substr(0, sep));
                    B1max = InputNum::parse_numeral(sB1.substr(sep + 1));
                }
                else
                {
                    B1next = InputNum::parse_numeral(sB1);
                    B1max = B1next;
                }
                //if (B1 < 100)
                //    B1 = 100;
            }
            else if (i < argc - 1 && strcmp(argv[i], "-B2") == 0)
            {
                i++;
                B2 = InputNum::parse_numeral(argv[i]);
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
            else if (strcmp(argv[i], "-poly") == 0)
            {
                poly = true;
                if (i < argc - 2 && strcmp(argv[i + 1], "threads") == 0)
                {
                    i += 2;
                    polyThreads = atoi(argv[i]);
                }
                if (i < argc - 1 && argv[i + 1][0] == 't')
                {
                    i++;
                    polyThreads = atoi(argv[i] + 1);
                }
            }
            else if (i < argc - 2 && strcmp(argv[i], "-generate") == 0)
            {
                generate = 1;
                i++;
                seed = atoi(argv[i]);
                i++;
                count = atoi(argv[i]);
            }
            else if (strcmp(argv[i], "-verify") == 0)
            {
                verify = 1;
                if (i < argc - 1 && strcmp(argv[i + 1], "curve") == 0)
                {
                    verifyCurve = 1;
                    i++;
                }
            }
            else if (strcmp(argv[i], "-mod") == 0)
            {
                modulus = 1;
                if (i < argc - 2 && strcmp(argv[i + 1], "curve") == 0)
                {
                    modCurve = atoi(argv[i + 2]);
                    i += 2;
                }
            }
            else if (i < argc - 2 && strcmp(argv[i], "-range") == 0)
            {
                range = true;
                i++;
                rangeOffset = atoi(argv[i]);
                i++;
                rangeCount = atoi(argv[i]);
            }
            else if (i < argc - 3 && strcmp(argv[i], "-split") == 0)
            {
                split = 1;
                i++;
                splitOffset = atoi(argv[i]);
                i++;
                splitCount = atoi(argv[i]);
                i++;
                splitName = argv[i];
            }
            else if (i < argc - 1 && strcmp(argv[i], "-merge") == 0)
            {
                merge = 1;
                i++;
                mergeName = argv[i];
            }
            else if (strcmp(argv[i], "-p54") == 0)
            {
                tmp = "568630647535356955169033410940867804839360742060818433";
                gwstate.known_factors *= tmp;
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
                printf("Prefactor-Factoring version " PREFACTOR_FACTORING_VERSION ", Gwnum library version " GWNUM_VERSION "\n");
                return 0;
            }
        }
        else
        {
            if (i < argc - 1 && strcmp(argv[i], "fermat") == 0)
            {
                i++;
                int fermat_n = atoi(argv[i]);
                input.init(1, 2, 1 << fermat_n, 1);
                if (fermat_n == 12)
                    gwstate.known_factors = "45477879701734570611058964078361695337745924097";
                if (fermat_n == 13)
                    gwstate.known_factors = "8314626596650587038214450998145116566054205961594349402971227571059785400321";
                if (fermat_n == 15)
                    gwstate.known_factors = "476875482933546652582154243389358929222507544696763973633";
                if (fermat_n == 16)
                    gwstate.known_factors = "156052367171184321737113706005266433";
                if (fermat_n == 17)
                    gwstate.known_factors = "31065037602817";
                if (fermat_n == 18)
                    gwstate.known_factors = "1107895052308076834874643709953";
                if (fermat_n == 19)
                    gwstate.known_factors = "1714509847183606156843894401498451927424901089206317613057";
                if (fermat_n == 21)
                    gwstate.known_factors = "4485296422913";
                if (fermat_n == 22)
                    gwstate.known_factors = "64658705994591851009055774868504577";
                if (fermat_n == 23)
                    gwstate.known_factors = "167772161";
            }
            else if (!input.empty())
                filename = argv[i];
            else if (!input.parse(argv[i]))
            {
                File file(argv[i], 0);
                if (!input.read(file))
                {
                    filename = argv[i];
                }
            }
        }
    if (input.empty() || filename.empty())
    {
        printf("Usage: prefactor -factoring {\"NUMBER\" | FILE | fermat N} [-generate CURVE COUNT] [-B1 10000] [-mod [curve CURVE]] [-verify [curve]] [-split OFFSET COUNT FILE] [-merge FILE] [-f FACTOR] FILE\n");
        return 0;
    }

    Logging logging(log_level);
    input.setup(gwstate);
    logging.info("Using %s.\n", gwstate.fft_description.data());

    try
    {
        Factoring factoring(input, gwstate, logging);
        logging.set_prefix(filename + ", ");
        File file_points(filename, gwstate.fingerprint);
        file_points.hash = false;
        std::string filename_state = filename + ".tmp";
        File file_state(filename_state, gwstate.fingerprint);
        file_state.hash = false;

        if (!factoring.read_points(file_points))
        {
            if (generate)
            {
                if (seed <= 0)
                {
                    logging.error("invalid curve #.\n");
                    return 1;
                }
                if (count <= 0)
                {
                    logging.error("invalid # of curves.\n");
                    return 1;
                }
                std::string dhash = factoring.generate(seed, count);
                factoring.write_points(file_points);
                write_dhash(filename + ".dhash", dhash);
                logging.info("generated %d curve%s starting with #%d.\n", count, count > 1 ? "s" : "", seed);
                logging.info("dhash: %s\n", dhash.data());
            }
            else
            {
                logging.error("file is missing or corrupted.\n");
                return 1;
            }
        }
        else if (generate)
            logging.warning("file exists, -generate ignored.\n");

        if (range)
            factoring.split(rangeOffset, rangeCount, factoring);

        if (B1next > factoring.B1())
        {
            if (range)
            {
                Factoring factoring_range(input, gwstate, logging);
                factoring.copy(factoring_range);
                write_dhash(filename + ".dhash", factoring_range.verify(false));
            }
            factoring.read_state(file_state, B1next);
            factoring.stage1(B1next, B1max, maxMem, file_state, file_points);
            remove(filename_state.data());
        }
        else
        {
            if (B1next != 0)
                logging.warning("file is at a higher B1.\n");
        }

        if (B2 > factoring.B1())
        {
            logging.progress().add_stage((int)factoring.points().size());
            factoring.stage2(B2, maxMem, poly, polyThreads, file_state);
            remove(filename_state.data());
            logging.progress().next_stage();
        }

        if (modulus)
        {
            factoring.modulus(modCurve, file_points);
        }

        if (split)
        {
            Factoring factoring_split(input, gwstate, logging);
            if (factoring.split(splitOffset, splitCount, factoring_split))
            {
                File file_split(splitName, gwstate.fingerprint);
                file_split.hash = false;
                factoring_split.write_points(file_split);
                write_dhash(splitName + ".dhash", factoring_split.verify(false));
            }
        }

        if (merge)
        {
            File file_merge(mergeName, gwstate.fingerprint);
            file_merge.hash = false;
            Factoring factoring_merge(input, gwstate, logging);
            if (factoring_merge.read_points(file_merge) && factoring.merge(factoring_merge))
            {
                factoring.write_points(file_points);
                write_dhash(filename + ".dhash", factoring.verify(false));
                verify = 0;
            }
        }

        if (verify)
        {
            std::string dhash = factoring.verify(verifyCurve);
            std::string dhashfile = filename + ".dhash";
            FILE* fp = fopen(dhashfile.data(), "r");
            if (fp)
            {
                char hash[33];
                if (fread(hash, 1, 32, fp) == 32)
                {
                    hash[32] = 0;
                    if (dhash != hash)
                        logging.error("dhash mismatch!\n");
                    else
                        logging.info("dhash ok.\n");
                }
                fclose(fp);
            }
        }
    }
    catch (const NoInverseException& e)
    {
        logging.report_factor(input, e.divisor);
    }
    catch (const TaskAbortException&)
    {
    }

    gwstate.done();

    return 0;
}
