#define PREFACTOR_FACTORING_VERSION "1.0.0"

#include <iostream>
#include <string.h>
#include <stdlib.h>
#include "gwnum.h"
#include "cpuid.h"
#include "factoring.h"
#include "exception.h"
#include "task.h"
#include "primelist.h"

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

void Factoring::write_file(File& file, char type, uint64_t B1, std::vector<std::unique_ptr<EdPoint>>& points)
{
    std::unique_ptr<Writer> writer(file.get_writer());
    writer->write(File::MAGIC_NUM);
    writer->write(FACTORING_APPID + (0 << 8) + (type << 16) + (0 << 24));
    writer->write(_gwstate.fingerprint);
    writer->write(_seed);
    writer->write((int)points.size());
    writer->write(B1);
    Giant X, Y, Z, T;
    for (auto it = points.begin(); it != points.end(); it++)
    {
        (*it)->serialize(X, Y, Z, T);
        writer->write(X);
        writer->write(Y);
        writer->write(Z);
        writer->write(T);
    }
    file.commit_writer(*writer);
}

bool Factoring::read_file(File& file, char type, int& seed, uint64_t& B0, std::vector<std::unique_ptr<EdPoint>>& points)
{
    file.appid = FACTORING_APPID;
    std::unique_ptr<Reader> reader(file.get_reader());
    if (!reader)
    {
        file.appid = 2;
        reader.reset(file.get_reader());
        if (!reader)
            return false;
    }
    else
    {
        if (!reader || reader->type() != type)
            return false;
        uint32_t fingerprint;
        if (!reader->read(fingerprint) || fingerprint != _gwstate.fingerprint)
            return false;
    }
    if (!reader->read(seed))
        return false;
    int count;
    if (!reader->read(count))
        return false;
    if (!reader->read(B0))
        return false;
    points.resize(count);
    for (int i = 0; i < count; i++)
    {
        points[i].reset(new EdPoint(_ed));
        Giant X, Y, Z, T;
        if (!reader->read(X) || !reader->read(Y) || !reader->read(Z) || !reader->read(T))
            return false;
        points[i]->deserialize(X, Y, Z, T);
    }
    return true;
}

void Factoring::write_points(File& file)
{
    write_file(file, 0, _B0, _points);
}

bool Factoring::read_points(File& file)
{
    int seed;
    uint64_t B0;
    std::vector<std::unique_ptr<EdPoint>> points;
    if (!read_file(file, 0, seed, B0, points))
        return false;
    _seed = seed;
    _B0 = B0;
    _points = std::move(points);
    _logging.info("%d curve%s starting with #%d, B1 = %" PRId64 ".\n", _points.size(), _points.size() > 1 ? "s" : "", _seed, _B0);
    return true;
}

bool Factoring::read_state(File& file, uint64_t B1)
{
    int seed;
    uint64_t b1;
    std::vector<std::unique_ptr<EdPoint>> points;
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
    _B0 = 1;
    _points.clear();
    _state.clear();

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

        Giant gx = (gw.popg() = x)%gw.N();
        Giant gy = (gw.popg() = y)%gw.N();
        Giant gx8 = (gw.popg() = x8)%gw.N();
        Giant gisdx8 = (gw.popg() = isdx8)%gw.N();
        Giant n1 = gw.N() - 1;
        Giant nx8 = gw.N() - gx8;
        Giant nisdx8 = gw.N() - gisdx8;

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

        Giant gd = (gw.popg() = square(sqrt_d))%gw.N();
        hasher.write((char*)gd.data(), gd.size()*4);

        _points.emplace_back(new EdPoint(_ed));
        _points.back()->X.reset(new GWNum(std::move(x)));
        _points.back()->Y.reset(new GWNum(std::move(y)));

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

std::string Factoring::verify(bool verify_curve)
{
    int i;
    for (i = 0; i < _points.size(); i++)
    {
        _ed.d_ratio(*_points[i], *_points[i]->X, *_points[i]->Y);
        _points[i]->Z.reset(_points[i]->Y.release());
    }
    try
    {
        _ed.normalize(_points.begin(), _points.end(), _ed.ED_PROJECTIVE);
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
                _points[i]->normalize();
            }
            catch (const NoInverseException&)
            {
                _logging.warning("curve #%d found the factor.\n", _seed + i);
                _ed.gen_curve(_seed + i, _points[i]->X.get());
                _points[i]->Z.reset();
            }
        }
    }

    Writer hasher;
    for (i = 0; i < _points.size(); i++)
    {
        Giant gd = (_gw.popg() = *_points[i]->X)%_gw.N();
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
        *_points[i]->X = (_gw.popg() = *_points[i]->X)%_gw.N();
        *_points[i]->Y = (_gw.popg() = *_points[i]->Y)%_gw.N();
    }
    write_file(file_result, 0, _B0, _points);
}

bool Factoring::split(int offset, int count, Factoring& result)
{
    if (offset < 0 || offset >= _points.size())
    {
        _logging.error("offset %d outside boundaries.\n", offset);
        return false;
    }
    _logging.info("splitting %d curve%s starting with #%d.\n", count, count > 1 ? "s" : "", _seed + offset);

    result._seed = _seed + offset;
    result._B0 = _B0;
    result._points.clear();
    result._state.clear();
    for (int i = 0; i < count; i++)
        result._points.emplace_back(new EdPoint(*_points[offset + i]));
    return true;
}

bool Factoring::merge(Factoring& other)
{
    if (other._B0 != _B0)
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

    for (auto it = other._points.begin(); it != other._points.end(); it++)
        _points.emplace_back(it->release());
    return true;
}

Giant get_exp(std::vector<uint64_t>& primes)
{
    uint64_t j, k;

    auto it = primes.begin();

    uint64_t sqrtB1 = (uint64_t)1e9;
    Giant tmp(GiantsArithmetic::default_arithmetic(), (int)primes.size()*2 + 10);
    Giant tmp2(GiantsArithmetic::default_arithmetic(), 8192 < tmp.capacity() ? 8192 : tmp.capacity());
    tmp2 = 1;
    Giant g64;
    g64 = 1;
    g64 <<= 32;
    while (it != primes.end())
    {
        // Building exponent with prime powers <= B1
        j = *it;
        if (*it <= sqrtB1)
        {
            k = ((uint64_t)1e18)/(*it);
            while (j <= k)
                j *= *it;
        }
        *(uint64_t*)g64.data() = j;
        tmp2 *= g64;
        it++;
        if (it == primes.end() || tmp2.size() > 8190)
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

void Factoring::stage1(uint64_t B1, File& file_state, File& file_result)
{
    int i, j;

    if (_B0 == 1)
    {
        _ed.set_gw(_gw.carefully());
        for (i = 0; i < _points.size(); i++)
            for (j = 0; j < 60; j++)
                _ed.dbl(*_points[i], *_points[i], _ed.ED_PROJECTIVE);
        _ed.normalize(_points.begin(), _points.end(), _ed.ED_PROJECTIVE);
        _B0 = 2;
        _ed.set_gw(_gw);
        if (_B0 == B1)
        {
            write_file(file_result, 0, B1, _points);
            return;
        }
    }

    PrimeList primes(65536);
    std::vector<uint64_t> plist;
    primes.begin().sieve_range(_B0 + 1, B1 + 1, plist);
    Giant tmp = get_exp(plist);
    plist.clear();

    int W;
    int len = tmp.bitlen();
    for (W = 2; W < 15 && (14 << (W - 2)) + len*(7 + 7/(W + 1.0)) > (14 << (W - 1)) + len*(7 + 7/(W + 2.0)); W++);
    std::vector<int16_t> naf_w;
    get_NAF_W(W, tmp, naf_w);

    _logging.info("%d bits, W = %d\n", len, W);
    _logging.progress().update(_state.size()/(double)_points.size(), (int)_state.size()*len/1000);

    time_t last_write = time(NULL);
    while (_state.size() < _points.size())
    {
        i = (int)_state.size();
        _state.emplace_back(new EdPoint(_ed));
        //double timer = getHighResTimer();
        _ed.mul(*_points[i], W, naf_w, *_state[i]);
        //timer = (getHighResTimer() - timer)/getHighResTimerFrequency();
        //_logging.info("%.1f%% done, %.3f ms per kilobit.\n", i/10.24, 1000000*timer/len);

        if ((time(NULL) - last_write > 300 || Task::abort_flag()) && _state.size() < _points.size())
        {
            _logging.progress().update(_state.size()/(double)_points.size(), (int)_state.size()*len/1000);
            _logging.report_progress();
            write_file(file_state, 1, B1, _state);
            last_write = time(NULL);
        }
        if (Task::abort_flag())
            throw TaskAbortException();
    }

    _ed.normalize(_state.begin(), _state.end(), _ed.ED_PROJECTIVE);
    write_file(file_result, 0, B1, _state);
    _B0 = B1;
    _points = std::move(_state);
}

int factoring_main(int argc, char *argv[])
{
    int i, j;
    GWState gwstate;
    GWArithmetic gw(gwstate);
    EdwardsArithmetic ed(gw);
    uint64_t B1 = 0;
    int seed = 0;
    int count = 0;
    std::string filename;
    int generate = 0;
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
    Giant factors;
    factors = "1";
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
                    factors = argv[i] + 2;
                else if (!argv[i][2] && i < argc - 1)
                {
                    i++;
                    factors = argv[i];
                }
                else
                    break;
                continue;
            }

            if (i < argc - 1 && strcmp(argv[i], "-B1") == 0)
            {
                i++;
                B1 = atoll(argv[i]);
                j = (int)strlen(argv[i]);
                if (argv[i][j - 1] == 'P')
                    B1 *= 1e15;
                if (argv[i][j - 1] == 'T')
                    B1 *= 1e12;
                if (argv[i][j - 1] == 'G')
                    B1 *= 1e9;
                if (argv[i][j - 1] == 'M')
                    B1 *= 1e6;
                if (argv[i][j - 1] == 'K')
                    B1 *= 1e3;
                //if (B1 < 100)
                //    B1 = 100;
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
                factors *= tmp;
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
                    factors = "45477879701734570611058964078361695337745924097";
            }
            else if (!input.parse(argv[i]))
            {
                File file(argv[i], 0);
                if (!input.read(file))
                {
                    filename = argv[i];
                }
            }
        }
    if (filename.empty())
    {
        printf("Usage: prefactor -factoring {\"NUMBER\" | FILE | fermat N} [-generate CURVE COUNT] [-B1 10000] [-mod [curve CURVE]] [-verify [curve]] [-split OFFSET COUNT FILE] [-merge FILE] [-f FACTOR] FILE\n");
        return 0;
    }

    Logging logging(log_level);
    input.setup(gwstate);
    logging.info("Using %s.\n", gwstate.fft_description.data());
    *gwstate.N /= factors;

    try
    {
        Factoring factoring(input, gwstate, logging);
        logging.set_prefix(filename + ", ");
        File file_points(filename, gwstate.fingerprint);
        file_points.hash = false;
        std::string filename_state = filename + ".tmp";
        File file_state(filename_state, gwstate.fingerprint);
        file_state.hash = false;

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
        }
        else if (!factoring.read_points(file_points))
        {
            logging.error("file is missing or corrupted.\n");
            return 1;
        }

        if (B1 > factoring.B0())
        {
            logging.progress().add_stage((int)factoring.points().size());
            factoring.read_state(file_state, B1);
            factoring.stage1(B1, file_state, file_points);
            remove(filename_state.data());
            logging.progress().next_stage();
        }
        else
        {
            if (B1 != 0)
                logging.warning("file is at a higher B1.\n");
        }

        if (modulus)
        {
            factoring.modulus(modCurve, file_points);
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
