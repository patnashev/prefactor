#define PREFACTOR_FERMAT_VERSION "0.8.0"

#include <stdlib.h>
#include "gwnum.h"
#include "cpuid.h"
#include "fermat.h"
#include "exception.h"
#include "task.h"
#include "primelist.h"

using namespace arithmetic;

void writeToFileFermat(const std::string& filename, char exponent, int id, uint64_t b0, std::vector<std::unique_ptr<EdPoint>>& points)
{
    Writer writer;
    writer.write(File::MAGIC_NUM);
    writer.write(FERMAT_APPID + (0 << 8) + ((int)exponent << 16));
    writer.write(id);
    writer.write((int)points.size());
    writer.write(b0);
    Giant X, Y, Z, T;
    for (auto it = points.begin(); it != points.end(); it++)
    {
        (*it)->serialize(X, Y, Z, T);
        writer.write(X);
        writer.write(Y);
        writer.write(Z);
        writer.write(T);
    }
    File file(filename, 0);
    file.hash = false;
    file.commit_writer(writer);
}

void gen_F12_ed_curve(const std::string& prefix, char exponent, int seed, int count, EdwardsArithmetic& ed)
{
    GWArithmetic& gw = ed.gw().carefully();
    std::vector<std::unique_ptr<EdPoint>> points;
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

    for (i = 0; i < 1024*count; i++)
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

        if (i < 1024)
            points.emplace_back(new EdPoint(ed));
        points[i & 1023]->X.reset(new GWNum(std::move(x)));
        points[i & 1023]->Y.reset(new GWNum(std::move(y)));

        if (((i + 1) & 1023) == 0)
        {
            //ed.normalize(points.begin(), points.end(), 0);

            std::string filename = prefix + std::to_string((seed + i - 1)/1024);
            writeToFileFermat(filename, exponent, seed + i - 1023, 1, points);

            std::string hash = hasher.hash_str();
            filename += ".dhash";
            FILE *fp = fopen(filename.data(), "a");
            if (fp)
            {
                fwrite(hash.data(), 1, hash.length(), fp);
                fclose(fp);
            }
            hasher.buffer().clear();
        }

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
}

void ed_mul(EdPoint& p, Giant& exp)
{
    int i, len;
    EdPoint p0 = p;
    len = exp.bitlen() - 1;
    for (i = 1; i <= len; i++)
    {
        p.arithmetic().dbl(p, p);
        if (exp.bit(len - i))
            p.arithmetic().add(p, p0, p);
    }
}

void Fermat::write_file(File& file, uint64_t B1, std::vector<std::unique_ptr<EdPoint>>& points)
{
    std::unique_ptr<Writer> writer(file.get_writer());
    writer->write(File::MAGIC_NUM);
    writer->write(FERMAT_APPID + (0 << 8) + ((int)_exponent_c << 16) + (0 << 24));
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

bool Fermat::read_file(File& file, int& seed, uint64_t& B0, std::vector<std::unique_ptr<EdPoint>>& points)
{
    file.appid = FERMAT_APPID;
    std::unique_ptr<Reader> reader(file.get_reader());
    if (!reader)
    {
        file.appid = FILE_APPID;
        reader.reset(file.get_reader());
        if (!reader)
            return false;
    }
    if (reader->type() != _exponent_c)
        return false;
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
        if (!reader->read(X) || !reader->read(Y) || !reader->read(Z) || (file.appid == FERMAT_APPID && !reader->read(T)))
            return false;
        if (file.appid != FERMAT_APPID)
        {
            if (Z == 0)
                Z = 1;
            T = 0;
        }
        points[i]->deserialize(X, Y, Z, T);
    }
    return true;
}

bool Fermat::read_points(File& file)
{
    int seed;
    uint64_t B0;
    std::vector<std::unique_ptr<EdPoint>> points;
    if (!read_file(file, seed, B0, points))
        return false;
    _seed = seed;
    _B0 = B0;
    _points = std::move(points);
    _logging.info("%s, %d curves, B1 = %" PRId64 ".\n", _input.display_text().data(), _points.size(), _B0);
    return true;
}

bool Fermat::read_state(File& file, uint64_t B1)
{
    int seed;
    uint64_t b1;
    std::vector<std::unique_ptr<EdPoint>> points;
    if (!read_file(file, seed, b1, points))
        return false;
    if (_seed != seed)
        return false;
    if (B1 != b1)
        return false;
    _state = std::move(points);
    _logging.info("Resuming from curve %d.\n", (int)_state.size());
    return true;
}

std::string Fermat::verify(bool verify_curve)
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
                _logging.warning("Curve #%d found the factor.\n", _seed + i);
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

void Fermat::modulus(int curve, File& file_result)
{
    int i, j;
    if (curve != 0 && curve - _seed >= 0 && curve - _seed < 1024)
    {
        i = curve - _seed;
        j = i + 1;
        _logging.info("Applying modulus to curve #%d (offset %d).\n", curve, i);
    }
    else if (curve != 0)
    {
        _logging.error("Curve #%d not found.\n", curve);
        return;
    }
    else
    {
        i = 0;
        j = (int)_points.size();
        _logging.info("Applying modulus to all curves.\n");
    }
    for (; i < j; i++)
    {
        *_points[i]->X = (_gw.popg() = *_points[i]->X)%_gw.N();
        *_points[i]->Y = (_gw.popg() = *_points[i]->Y)%_gw.N();
    }
    write_file(file_result, _B0, _points);
}

Giant get_exp(std::vector<int>& primes)
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

void Fermat::stage1(uint64_t B1, File& file_state, File& file_result)
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
            write_file(file_result, B1, _points);
            return;
        }
    }

    PrimeList primes(65536);
    std::vector<int> plist;
    primes.sieve_range((int)_B0 + 1, (int)B1 + 1, plist);
    Giant tmp = get_exp(plist);
    plist.clear();

    int W;
    int len = tmp.bitlen();
    for (W = 2; W < 15 && (14 << (W - 2)) + len*(7 + 7/(W + 1.0)) > (14 << (W - 1)) + len*(7 + 7/(W + 2.0)); W++);
    std::vector<int16_t> naf_w;
    get_NAF_W(W, tmp, naf_w);

    _logging.info("%d bits, W = %d\n", len, W);

    time_t last_write = time(NULL);
    double timer = getHighResTimer();
    int timer_i = (int)_state.size();
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
            timer = (getHighResTimer() - timer)/getHighResTimerFrequency();
            write_file(file_state, B1, _state);
            last_write = time(NULL);
            _logging.info("%.1f%% done, %.3f ms per kilobit.\n", i/10.24, 1000000*timer/len/(i + 1 - timer_i));
            timer = getHighResTimer();
            timer_i = i + 1;
        }
        if (Task::abort_flag())
            throw TaskAbortException();
    }

    _ed.normalize(_state.begin(), _state.end(), _ed.ED_PROJECTIVE);
    write_file(file_result, B1, _state);
    _B0 = B1;
    _points = std::move(_state);
}

int fermat_main(int argc, char *argv[])
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
    int exponent = 4096;
    char exponent_c = 12;
    int log_level = Logging::LEVEL_INFO;
    Giant factors;
    factors = "45477879701734570611058964078361695337745924097";
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

            case 'n':
                if (argv[i][2] && isdigit(argv[i][2]))
                    exponent = atoi(argv[i] + 2);
                else if (!argv[i][2] && i < argc - 1)
                {
                    i++;
                    exponent = atoi(argv[i]);
                }
                else
                    break;
                for (j = 1, exponent_c = 0; j < exponent; j <<= 1, exponent_c++);
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
                //if (B1 < 100)
                //    B1 = 100;
            }
            else if (i < argc - 2 && strcmp(argv[i], "-generate") == 0)
            {
                generate = 1;
                i++;
                seed = atoi(argv[i])*1024 + 1;
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
            else if (strcmp(argv[i], "-p54") == 0)
            {
                tmp = "568630647535356955169033410940867804839360742060818433";
                factors *= tmp;
            }
            else if (strcmp(argv[i], "-v") == 0)
            {
                printf("Prefactor-Fermat version " PREFACTOR_FERMAT_VERSION ", Gwnum library version " GWNUM_VERSION "\n");
                return 0;
            }
        }
        else
        {
            filename = argv[i];
        }
    if (filename.empty())
    {
        printf("Usage: prefactor -fermat [-n 4096] [-B1 10000] [-generate OFFSET COUNT PREFIX] [-mod [curve 123456]] [-verify [curve]] [-f knownFactors] [-p54] file\n");
        return 0;
    }

    InputNum input(1, 2, exponent, 1);
    Logging logging(log_level);
    input.setup(gwstate);
    logging.info("Using %s.\n", gwstate.fft_description.data());
    *gwstate.N /= factors;

    try
    {
        if (generate)
        {
            gen_F12_ed_curve(filename, exponent_c, seed, count, ed);
            return 0;
        }

        if (0)
        {
            factors = "568630647535356955169033410940867804839360742060818433";
            EdPoint p0 = ed.gen_curve(7636607, nullptr);
            EdPoint p(ed);
            PrimeList primes(3000000);
            std::vector<int> powers(primes.size());
            uint64_t sqrtB1 = (uint64_t)1e9;
            for (i = 0; i < primes.size(); i++)
            {
                uint64_t pp, k;
                // Building exponent with prime powers <= B1
                pp = primes[i];
                powers[i] = 1;
                if (primes[i] <= sqrtB1)
                {
                    k = ((uint64_t)1e18)/primes[i];
                    while (pp <= k)
                    {
                        pp *= primes[i];
                        powers[i]++;
                    }
                }
            }

            std::vector<int> rprimes;
            std::vector<int> rpowers;
            int r = (int)primes.size();
            EdPoint pl = p0;
            while (r > 1)
            {
                int l = 0;
                while (l < r - 1)
                {
                    p = pl;
                    int t = (l + r)/2;
                    for (i = l; i < t; i++)
                    {
                        tmp = primes[i];
                        for (j = 0; j < powers[i]; j++)
                            ed_mul(p, tmp);
                    }
                    if ((gw.popg() = *p.X)%factors == 0)
                    {
                        r = t;
                    }
                    else
                    {
                        pl = p;
                        l = t;
                    }
                }
                rprimes.push_back(primes[l]);
                rpowers.push_back(powers[l]);
                pl = p0;
                for (i = 0; i < rprimes.size(); i++)
                {
                    tmp = rprimes[i];
                    for (j = 0; j < rpowers[i]; j++)
                        ed_mul(pl, tmp);
                }
            }
            for (r = 0; r < rprimes.size(); r++)
            {
                while (1)
                {
                    rpowers[r]--;
                    pl = p0;
                    for (i = 0; i < rprimes.size(); i++)
                    {
                        tmp = rprimes[i];
                        for (j = 0; j < rpowers[i]; j++)
                            ed_mul(pl, tmp);
                    }
                    if ((gw.popg() = *pl.X)%factors != 0)
                    {
                        rpowers[r]++;
                        break;
                    }
                }
            }
            for (i = 0; i < rprimes.size(); i++)
            {
                if (i > 0)
                    printf(" * ");
                else
                    printf("#E = ");
                if (rpowers[i] > 1)
                    printf("%d^%d", rprimes[i], rpowers[i]);
                else
                    printf("%d", rprimes[i]);
            }
            printf("\n");
        }

        Fermat fermat(exponent, filename, gwstate, logging);
        File file_points(filename, 0);
        file_points.hash = false;
        std::string filename_state = filename + ".tmp";
        File file_state(filename_state, 0);
        file_state.hash = false;

        if (!fermat.read_points(file_points))
        {
            logging.warning("File %s is missing or corrupted.\n", filename.data());
            return 1;
        }

        if (B1 > fermat.B0())
        {
            fermat.read_state(file_state, B1);
            fermat.stage1(B1, file_state, file_points);
            remove(filename_state.data());
        }
        else
        {
            if (B1 != 0)
                logging.warning("File %s is at a higher B1.\n", filename.data());
        }

        if (modulus)
        {
            fermat.modulus(modCurve, file_points);
        }

        if (verify)
        {
            std::string dhash = fermat.verify(verifyCurve);
            std::string dhashfile = filename + ".dhash";
            FILE* fp = fopen(dhashfile.data(), "r");
            if (fp)
            {
                char hash[33];
                if (fread(hash, 1, 32, fp) == 32)
                {
                    hash[32] = 0;
                    if (dhash != hash)
                        logging.error("File %s dhash mismatch!\n", filename.data());
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
