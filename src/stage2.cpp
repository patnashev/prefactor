
#include <deque>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <string.h>

#include "gwnum.h"
#include "cpuid.h"
#include "stage2.h"
#include "exception.h"

using namespace arithmetic;

#ifdef _DEBUG
#define DEBUG_STAGE2_PRECOMP
#endif

void Stage2::init(InputNum* input, GWState* gwstate, File* file, TaskState* state, Logging* logging, int iterations)
{
    Task::init(gwstate, file, state, logging, iterations);
    _error_check = gwnear_fft_limit(gwstate->gwdata(), 1) == TRUE;
    _input = input;
    _timer = getHighResTimer();
    _transforms = -(int)gwstate->handle.fft_count;
    _success = false;
}

void Stage2::reinit_gwstate()
{
    double fft_count = _gwstate->handle.fft_count;
    _gwstate->done();
    _input->setup(*_gwstate);
    _gwstate->handle.fft_count = fft_count;
    std::string prefix = _logging->prefix();
    _logging->set_prefix("");
    _logging->error("Restarting using %s\n", _gwstate->fft_description.data());
    _logging->set_prefix(prefix);
    _logging->report_param("fft_desc", _gwstate->fft_description);
    _logging->report_param("fft_len", _gwstate->fft_length);
}

void Stage2::done(const arithmetic::Giant& factor)
{
    _timer = (getHighResTimer() - _timer)/getHighResTimerFrequency();
    _transforms += (int)_gwstate->handle.fft_count;
    _logging->progress().update(1, (int)_gwstate->handle.fft_count/2);
    _logging->info("transforms: %d, time: %.1f s.\n", _transforms, _timer);
    if (factor == 0 || factor == *_gwstate->N)
    {
        _success = true;
        _logging->result(_success, "all divisors less than B1.\n");
        _logging->result_save("all divisors of " + _input->input_text()  + " are less than B1.\n");
    }
    else if (factor != 1)
    {
        _success = true;
        _logging->report_factor(*_input, factor);
    }
    else
    {
        _logging->result(_success, "RES64: %s.\n", _res64.data());
    }
    _logging->set_prefix("");
}

template<int TL>
Stage2Pairing::Pairing get_pairing_L(Logging& logging, PrimeList& primes, int B1, int B2, int D, int A, int L, bool with_distances)
{
    int i, j, k;
    int d, p, q;
    Stage2Pairing::Pairing ret;

    double timer = getHighResTimer();

    int second_base = D/A;
    //int second_base = D/A*(A-1)/2;

    std::vector<int> dist;
    /*for (i = 1; i < D; i++)
        if (gcd(D, i) == 1)
        {
            dist[dist_count] = i;
            dist_count++;
        }
    k = dist_count/2;
    for (j = 1; j < L - 1; j++)
        for (i = 0; i < k; i++)
        {
            dist[dist_count] = dist[i] + ((1 << j) - 1)*D;
            dist_count++;
        }*/
        //for (i = 3; i < D*L + (A > 1 ? D/A : 0); i++)
    for (i = 1; i < D/2*L; i++)
        if (gcd(D/A, i%(D/2)) == 1)
        {
            //d = i;
            //d = i + i/D*D;
            //d = i%D + ((1 << (i/D*1)) - 1)*D;
            d = i%(D/2) + ((1 << (i/(D/2)*1)) - 1)*D;
            //d = i%(D/2) + (i/(D/2) == 1 ? -1 : i/(D/2) == 2 ? 2 : i/(D/2) == 3 ? 6 : i/(D/2) == 4 ? 14 : i/(D/2) == 5 ? 30 : 0)*D;
            //d = i%(D/2) + (i/(D/2) == 1 ? 1 : i/(D/2) == 2 ? 7 : i/(D/2) == 3 ? 31 : i/(D/2) == 4 ? 127 : i/(D/2) == 5 ? 511 : 0)*D;
            //d = i%(D/2) + (i/(D/2) == 1 ? 1 : i/(D/2) == 2 ? 4 : i/(D/2) == 3 ? 16 : i/(D/2) == 4 ? 64 : i/(D/2) == 5 ? 256 : 0)*D;
            dist.push_back(d);
        }
    for (k = (int)dist.size(), i = 0; i < k; i++)
        dist.push_back(-dist[i]);
    int max_dist = d;

    std::vector<int> relocs;
    for (j = 3; B1*j <= B2; j += 2)
        if (gcd(D/A, j) == 1)
            relocs.push_back(j);

    int *map = new int[(B2 + max_dist)/2 + 1];
    memset(map, 0, sizeof(int)*((B2 + max_dist)/2 + 1));

    PrimeIterator it = primes.begin();
    for (; *it <= B1; it++);

    std::vector<Stage2Pairing::prime<TL>> plist;
    std::vector<Stage2Pairing::prime<1>> srclist;
    plist.emplace_back(0);
    srclist.emplace_back(0);
    for (ret.total = 0; *it <= B2; it++, ret.total++)
    {
        if (!relocs.empty() && *it*relocs[0] <= B2)
        {
            srclist.emplace_back((int)plist.size());
            srclist.back().adjacency[0] = (int)plist.size();
            for (i = 0; i < relocs.size() && *it*relocs[i] <= B2 + max_dist; i++)
                if (*it*relocs[i] > B2/relocs[0] - max_dist)
                {
                    if (*it*relocs[i] <= B2)
                        srclist.back().value = (int)plist.size();
                    map[*it*relocs[i]/2] = (int)plist.size();
                    plist.emplace_back(*it*relocs[i]);
                    plist.back().match = 1 - (int)srclist.size();
                }
            if (srclist.back().adjacency[0] == plist.size())
                throw ArithmeticException("Invalid B1/B2.");
            if (*it > B2/relocs[0] - max_dist)
            {
                map[*it/2] = (int)plist.size();
                plist.emplace_back(*it);
                plist.back().match = 1 - (int)srclist.size();
            }
        }
        else
        {
            map[*it/2] = (int)plist.size();
            plist.emplace_back(*it);
        }
    }

    logging.info("Pairing %d primes, D=%d, L=%d", ret.total, D, L);
    if (A > 1)
        logging.info(", A=%d", A);
    logging.info("... ");

    std::vector<std::vector<int>> dist_rem(D);
    for (i = 1; i < D; i++)
        if (gcd(D/A, i) == 1)
        {
            for (j = 0; j < dist.size(); j++)
            {
                d = i + dist[j];
                if ((d%D != 0 && (d%D + D)%D != second_base))
                    continue;
                dist_rem[i].push_back(dist[j]*2);
            }
        }
    if (!relocs.empty())
        d = B2/relocs[0]/D*D;
    else
        d = B1/D*D;
    for (auto itp = plist.begin(); itp != plist.end(); itp++)
    {
        j = 0;
        p = itp->value;
        auto& p_dists = dist_rem[p%D];
        for (auto itd = p_dists.begin(); itd != p_dists.end(); itd++)
        {
            q = p + *itd;
            if (q <= B1 || q > B2 + max_dist || map[q/2] == 0 || p + q < 2*d || p + q > 2*B2)
                continue;
            if (itp->match < 0 && itp->match == plist[map[q/2]].match)
                continue;
            itp->adjacency[j] = map[q/2];
            j++;
        }
        if (j > TL)
            throw ArithmeticException("Incorrect L.");
        if (j < TL)
            itp->adjacency[j] = 0;
    }
    delete map;

    for (p = 0; p < plist.size(); p++)
    {
        if (plist[p].match > 0 || (plist[p].match < 0 && srclist[-plist[p].match].match > 0))
            continue;
        int* adj = plist[p].adjacency;
        for (i = 0; i < TL && adj[i] != 0; i++)
        {
            q = adj[i];
            if (plist[p].value > plist[q].value)
                continue;
            if (plist[q].match == 0 || (plist[q].match < 0 && srclist[-plist[q].match].match == 0))
            {
                if (plist[p].match < 0)
                {
                    srclist[-plist[p].match].value = p;
                    srclist[-plist[p].match].match = q;
                }
                else
                    plist[p].match = q;
                if (plist[q].match < 0)
                {
                    srclist[-plist[q].match].value = q;
                    srclist[-plist[q].match].match = p;
                }
                else
                    plist[q].match = p;
                break;
            }
        }
    }

    std::deque<int> queue;
    int *plink = new int[plist.size()];
    int *srclink = new int[srclist.size()];
    k = 0;
    int kk = 0;
    int pairs = 0;
    bool flag = true;
    while (flag)
    {
        flag = false;
        kk = k%plist.size();
        memset(plink, 0, sizeof(int)*plist.size());
        memset(srclink, 0, sizeof(int)*srclist.size());
        for (k = kk; k < plist.size() + kk && !flag; k++)
        {
            p = k < (int)plist.size() ? k : k - (int)plist.size();
            if (p == 0 || plist[p].match > 0 || (plist[p].match < 0 && srclist[-plist[p].match].match > 0))
                continue;
            queue.clear();
            queue.push_back(p);
            plink[p] = -1;
            if (plist[p].match < 0)
                srclink[-plist[p].match] = p;

            for (; !queue.empty() && !flag; queue.pop_front())
            {
                int* adj = plist[queue.front()].adjacency;
                for (i = 0; i < TL && adj[i] != 0; i++)
                {
                    p = queue.front();
                    q = adj[i];
                    if ((plist[q].match >= 0 && plink[q] != 0) || (plist[q].match < 0 && srclink[-plist[q].match] != 0))
                        continue;
                    if (plist[q].match == 0 || (plist[q].match < 0 && srclist[-plist[q].match].match == 0))
                    {
                        while (1)
                        {
                            if (plist[p].match < 0)
                            {
                                srclist[-plist[p].match].value = p;
                                srclist[-plist[p].match].match = q;
                            }
                            else
                                plist[p].match = q;
                            if (plist[q].match < 0)
                            {
                                srclist[-plist[q].match].value = q;
                                srclist[-plist[q].match].match = p;
                            }
                            else
                                plist[q].match = p;
                            if (plist[p].match < 0)
                                p = srclink[-plist[p].match];
                            if (plink[p] == -1)
                                break;
                            q = plink[p];
                            if (plist[q].match < 0)
                                q = srclink[-plist[q].match];
                            p = plink[q];
                        }
                        pairs++;
                        flag = true;
                        break;
                    }

                    plink[q] = p;
                    if (plist[q].match < 0)
                    {
                        srclink[-plist[q].match] = q;
                        q = srclist[-plist[q].match].value;
                        p = srclist[-plist[q].match].match;
                    }
                    else
                        p = plist[q].match;
                    plink[p] = q;
                    if (plist[p].match < 0)
                    {
                        srclink[-plist[p].match] = p;
                        for (j = srclist[-plist[p].match].adjacency[0]; j < plist.size() && plist[j].match == plist[p].match; j++)
                            queue.push_back(j);
                    }
                    else
                        queue.push_back(p);
                }
            }

            if (Task::abort_flag())
                throw TaskAbortException();
        }
    }
    delete plink;
    delete srclink;

    ret.pairs = 0;
    if (!with_distances)
    {
        ret.first_D = B1/D;
        if (!relocs.empty())
            ret.first_D = B2/relocs[0]/D;
        ret.last_D = B2/D;
        for (auto itp = plist.begin(); itp != plist.end(); itp++)
            if (itp->match > 0 && plist[itp->match].value > itp->value)
                ret.pairs++;
        for (auto its = srclist.begin(); its != srclist.end(); its++)
            if (its->match > 0 && plist[its->match].value > plist[its->value].value)
                ret.pairs++;
    }
    else
    {
        std::vector<std::pair<int, int>> D_distances;
        for (auto itp = plist.begin() + 1; itp != plist.end(); itp++)
            if (itp->match == 0 && itp->value%D < D/2)
                D_distances.emplace_back(itp->value - itp->value%D, itp->value%D);
            else if (itp->match == 0)
                D_distances.emplace_back(itp->value + D - itp->value%D, D - itp->value%D);
            else if (itp->match > 0 && plist[itp->match].value > itp->value)
            {
                D_distances.emplace_back((plist[itp->match].value + itp->value)/2, (plist[itp->match].value - itp->value)/2);
                ret.pairs++;
            }
        for (auto its = srclist.begin() + 1; its != srclist.end(); its++)
            if (its->match == 0 && plist[its->value].value%D < D/2)
                D_distances.emplace_back(plist[its->value].value - plist[its->value].value%D, plist[its->value].value%D);
            else if (its->match == 0)
                D_distances.emplace_back(plist[its->value].value + D - plist[its->value].value%D, D - plist[its->value].value%D);
            else if (its->match > 0 && plist[its->match].value > plist[its->value].value)
            {
                D_distances.emplace_back((plist[its->match].value + plist[its->value].value)/2, (plist[its->match].value - plist[its->value].value)/2);
                ret.pairs++;
            }
        
        std::sort(D_distances.begin(), D_distances.end(), [](const std::pair<int, int>& a, const std::pair<int, int>& b) { return a.first < b.first; });
        ret.first_D = D_distances.front().first/D;
        int cur = ret.first_D*D;
        auto it = D_distances.begin();
        while (it != D_distances.end())
        {
            for (; it != D_distances.end() && cur == it->first; it++)
                ret.distances.push_back(it->second);
            ret.distances.push_back(0);
            if (A > 1)
            {
                for (; it != D_distances.end() && cur + second_base == it->first; it++)
                    ret.distances.push_back(it->second);
                ret.distances.push_back(0);
            }
            cur += D;
        }
        ret.last_D = cur/D - 1;
    }

    timer = (getHighResTimer() - timer)/getHighResTimerFrequency();
    logging.info("%d pairs (%.1f%%), time: %.1f s.\n", ret.pairs, 200.0*ret.pairs/ret.total, timer);

    return ret;
}

Stage2Pairing::Pairing Stage2Pairing::get_pairing(Logging& logging, PrimeList& primes, int B1, int B2, int D, int A, int L, bool with_distances)
{
    if (A > 1)
    {
        if (L == 1)
            return get_pairing_L<2>(logging, primes, B1, B2, D, A, L, with_distances);
        if (L == 2)
            return get_pairing_L<4>(logging, primes, B1, B2, D, A, L, with_distances);
        if (L == 3)
            return get_pairing_L<6>(logging, primes, B1, B2, D, A, L, with_distances);
        if (L == 4)
            return get_pairing_L<8>(logging, primes, B1, B2, D, A, L, with_distances);
        if (L == 5)
            return get_pairing_L<10>(logging, primes, B1, B2, D, A, L, with_distances);
        if (L == 6)
            return get_pairing_L<12>(logging, primes, B1, B2, D, A, L, with_distances);
        if (L == 7)
            return get_pairing_L<14>(logging, primes, B1, B2, D, A, L, with_distances);
        if (L == 8)
            return get_pairing_L<16>(logging, primes, B1, B2, D, A, L, with_distances);
        if (L == 9)
            return get_pairing_L<18>(logging, primes, B1, B2, D, A, L, with_distances);
        if (L == 10)
            return get_pairing_L<20>(logging, primes, B1, B2, D, A, L, with_distances);
        if (L == 11)
            return get_pairing_L<22>(logging, primes, B1, B2, D, A, L, with_distances);
        if (L == 12)
            return get_pairing_L<24>(logging, primes, B1, B2, D, A, L, with_distances);
        if (L == 13)
            return get_pairing_L<26>(logging, primes, B1, B2, D, A, L, with_distances);
        if (L == 14)
            return get_pairing_L<28>(logging, primes, B1, B2, D, A, L, with_distances);
        if (L == 15)
            return get_pairing_L<30>(logging, primes, B1, B2, D, A, L, with_distances);
    }
    else
    {
        if (L == 1)
            return get_pairing_L<1>(logging, primes, B1, B2, D, A, L, with_distances);
        if (L == 2)
            return get_pairing_L<2>(logging, primes, B1, B2, D, A, L, with_distances);
        if (L == 3)
            return get_pairing_L<3>(logging, primes, B1, B2, D, A, L, with_distances);
        if (L == 4)
            return get_pairing_L<4>(logging, primes, B1, B2, D, A, L, with_distances);
        if (L == 5)
            return get_pairing_L<5>(logging, primes, B1, B2, D, A, L, with_distances);
        if (L == 6)
            return get_pairing_L<6>(logging, primes, B1, B2, D, A, L, with_distances);
        if (L == 7)
            return get_pairing_L<7>(logging, primes, B1, B2, D, A, L, with_distances);
        if (L == 8)
            return get_pairing_L<8>(logging, primes, B1, B2, D, A, L, with_distances);
        if (L == 9)
            return get_pairing_L<9>(logging, primes, B1, B2, D, A, L, with_distances);
        if (L == 10)
            return get_pairing_L<10>(logging, primes, B1, B2, D, A, L, with_distances);
        if (L == 11)
            return get_pairing_L<11>(logging, primes, B1, B2, D, A, L, with_distances);
        if (L == 12)
            return get_pairing_L<12>(logging, primes, B1, B2, D, A, L, with_distances);
        if (L == 13)
            return get_pairing_L<13>(logging, primes, B1, B2, D, A, L, with_distances);
        if (L == 14)
            return get_pairing_L<14>(logging, primes, B1, B2, D, A, L, with_distances);
        if (L == 15)
            return get_pairing_L<15>(logging, primes, B1, B2, D, A, L, with_distances);
    }
    throw ArithmeticException("L not supported.");
}

template<class Element>
int Stage2Pairing::precompute(DifferentialGroupArithmetic<Element>& arithmetic, Element& X1, Element& XD, Element& XDA, std::vector<std::unique_ptr<Element>>& precomp)
{
    int i, j;
    int v;
    Element Xn(X1.arithmetic());
    Element Xn1(X1.arithmetic());

    precomp.resize((1 << (_L - 1))*_D);
    precomp[0].reset(new Element(X1.arithmetic()));
    *precomp[0] = X1;
    Xn = X1;
    arithmetic.dbl(Xn, Xn); // V_2
    arithmetic.dbl(Xn, Xn1); // V_4

    if (_A > 1)
        _D /= _A;
    int dist = 1;
    int precomp_size = 1;
    for (i = 1; i < (_L == 1 ? _D/4*_A : _D/2*_A); i++)
    {
        if (gcd(2*i + 1, 2*dist) != 1)
            continue;
        // V_{2i+1}
        precomp[i].reset(new Element(X1.arithmetic()));
        precomp_size++;
        arithmetic.add(Xn, *precomp[i - dist], *precomp[i >= 2*dist ? i - 2*dist : 2*dist - i - 1], *precomp[i]);
        if ((i + 1)%dist == 0 && _D%(i + 1) == 0 && gcd((i + 1)/dist, 2*dist) == 1)
        {
            int pd = (i + 1)/dist;
            arithmetic.add(*precomp[dist*(pd - 1)/2 - 1], *precomp[dist*(pd + 1)/2], Xn1, Xn1);
            arithmetic.add(*precomp[dist*(pd - 1)/2], *precomp[dist*(pd + 1)/2], Xn, Xn);
            swap(Xn, Xn1);
            dist = i + 1;
            for (j = 0; j < i; j++)
                if (precomp[j] && gcd(2*j + 1, 2*dist) != 1)
                {
                    precomp[j].reset();
                    precomp_size--;
                }
        }
    }
    // V_{D/A}
    if (_A > 1)
    {
        XDA = Xn;
        if (_D/2/dist != 1)
        {
            int pd = _D/2/dist;
            for (j = 0; (pd & 1) == 0; pd >>= 1, j++);
            if (pd > 1)
                arithmetic.add(*precomp[dist*(pd - 1)/2 - 1], *precomp[dist*(pd + 1)/2], Xn1, XDA);
            for (; j > 0; j--)
                arithmetic.dbl(XDA, XDA);
        }
        _D *= _A;
#ifdef DEBUG_STAGE2_PRECOMP
        XD = X1;
        XD *= _D/_A;
        GWASSERT(XD == XDA);
#endif
    }
    // V_D
    swap(XD, Xn);
    if (_L == 1 && _D%4 != 0)
        arithmetic.init(XD); // No way to compute XD
    else
    {
        if (_D/2/dist != 1)
        {
            int pd = _D/2/dist;
            for (j = 0; (pd & 1) == 0; pd >>= 1, j++);
            if (pd > 1)
                arithmetic.add(*precomp[dist*(pd - 1)/2 - 1], *precomp[dist*(pd + 1)/2], Xn1, XD);
            for (; j > 0; j--)
                arithmetic.dbl(XD, XD);
        }
#ifdef DEBUG_STAGE2_PRECOMP
        Xn = X1;
        Xn *= _D;
        GWASSERT(Xn == XD);
#endif
    }
    // Cleanup
    for (j = 0; j < (_L == 1 ? _D/4 : _D/2); j++)
        if (precomp[j] && gcd(2*j + 1, _D) != 1)
        {
            precomp[j].reset();
            precomp_size--;
        }

    // Irregular L
    Xn = XD;
    if (_L > 2)
        for (j = _D/4 + 0; j < _D/2; j++)
            if (precomp[j])
                arithmetic.optimize(*precomp[j]);
    for (i = 1; i < _L; i++)
    {
        for (j = 0; j < _D/4 + 0; j++)
            if (precomp[j])
            {
                v = ((1 << i) - 1)*_D/2 + j;
                precomp[v].reset(new Element(X1.arithmetic()));
                precomp_size++;
                arithmetic.add(Xn, *precomp[v - (1 << (i - 1))*_D/2], *precomp[_D/2 - j - 1], *precomp[v]);
            }
        arithmetic.dbl(Xn, Xn);
    }
    for (j = _D/4 + 0; j < _D/2; j++)
        if (precomp[j])
        {
            precomp[j].reset();
            precomp_size--;
        }
#ifdef DEBUG_STAGE2_PRECOMP
    for (i = 1, j = 0; i < _D/2*_L; i++)
        if (gcd(_D/_A, i%(_D/2)) == 1)
        {
            int d = i%(_D/2) + ((1 << (i/(_D/2)*1)) - 1)*_D;
            Xn = X1;
            Xn *= d;
            GWASSERT(Xn == *precomp[d/2]);
            j++;
        }
    GWASSERT(j == precomp_size);
#endif

    return precomp_size;
}

void PP1Stage2::init(InputNum* input, GWState* gwstate, File* file, Logging* logging, Giant& P, bool minus1)
{
    Stage2::init(input, gwstate, file, read_state<State>(file), logging, _pairing.last_D - _pairing.first_D + 1);
    _state_update_period = MULS_PER_STATE_UPDATE/10;
    _logging->set_prefix(input->display_text() + (minus1 ? ", P-1 stage 2, " : ", P+1 stage 2, "));
    _logging->info("B2 = %" PRId64, _B2);
    if (state() != nullptr)
        _logging->info(", restarting at %.1f%%", 100.0*state()->iteration()/iterations());
    _logging->info(".\n");
    lucas.reset(new LucasVArithmetic());
    _P = P;
}

void PP1Stage2::setup()
{
    lucas->set_gw(gw());

    if (!_Vn)
        _Vn.reset(new LucasV(*lucas));
    if (!_Vn1)
        _Vn1.reset(new LucasV(*lucas));
    if (!_W)
        _W.reset(new LucasV(*lucas));
    if (!_Va && _A > 1)
        _Va.reset(new LucasV(*lucas));
    if (!_Va1 && _A > 1)
        _Va1.reset(new LucasV(*lucas));
    if (!_Wa && _A > 1)
        _Wa.reset(new LucasV(*lucas));

    if (_precomp.empty())
    {
        int transforms = -(int)gw().gwdata()->fft_count;
        _Vn->V() = _P;
        std::vector<std::unique_ptr<LucasV>> precomp;
        int precomp_size = precompute<LucasV>(*lucas, *_Vn, *_W, _Wa ? *_Wa : *_Vn1, precomp);
        if (_L == 1 && _D%4 != 0)
            lucas->mul(*_Vn, _D, *_W, *_Vn1);
        commit_setup();
        _precomp = std::move(precomp);

        transforms += (int)gw().gwdata()->fft_count;
        _transforms -= transforms;
        _logging->info("%d precomputed values (%d transforms), %d steps.\n", precomp_size, transforms, _pairing.last_D - _pairing.first_D + 1);
    }

    if (state() == nullptr)
        reset_state<State>();
    int v = state()->iteration() + _pairing.first_D;

    if (_A > 1)
    {
        lucas->init(*_Va);
        *_Va1 = *_Wa;
        if (v > 0)
            lucas->mul(*_Wa, v*_A, *_Va, *_Va1);
        swap(*_Vn1, *_Va);
        *_Vn = *_Vn1;
        *_Va = *_Va1;
        for (int i = 0; i < _A; i++)
        {
            lucas->add(*_Wa, *_Va1, *_Vn1, *_Vn1);
            swap(*_Vn1, *_Va1);
        }
    }
    else
    {
        lucas->init(*_Vn);
        *_Vn1 = *_W;
        if (v > 0)
            lucas->mul(*_W, v, *_Vn, *_Vn1);
    }
    commit_setup();
}

void PP1Stage2::release()
{
    _precomp.clear();
    _Vn.reset();
    _Vn1.reset();
    _W.reset();
    _Va.reset();
    _Va1.reset();
    _Wa.reset();
}

void PP1Stage2::execute()
{
    int v;
    Giant tmp;

    lucas->set_gw(gw());
    GWNum G(gw());
    G = state()->G();

    v = _pairing.first_D;
    auto it = _pairing.distances.begin();
    while (v - _pairing.first_D < state()->iteration() && it != _pairing.distances.end())
    {
        for (; *it != 0; it++);
        it++;
        v++;
        if (_A > 1)
        {
            for (; *it != 0; it++);
            it++;
        }
    }
    v = state()->iteration() + _pairing.first_D;

    while (it != _pairing.distances.end())
    {
        for (; *it != 0; it++)
            gw().submul(_Vn->V(), _precomp[*it/2]->V(), G, G, GWMUL_STARTNEXTFFT_IF(!is_last(v - _pairing.first_D)));
        it++;
        if (v > 0)
            lucas->add(*_W, *_Vn1, *_Vn, *_Vn);
        else
            lucas->dbl(*_Vn1, *_Vn);
        swap(*_Vn, *_Vn1);
        if (_A > 1)
        {
            for (; *it != 0; it++)
                gw().submul(_Va->V(), _precomp[*it/2]->V(), G, G, GWMUL_STARTNEXTFFT_IF(!is_last(v - _pairing.first_D)));
            it++;
            lucas->add(*_W, *_Va1, *_Va, *_Va);
            swap(*_Va, *_Va1);
        }
        v++;
        commit_execute<State>(v - _pairing.first_D, G);
    }

    tmp = G;
    tmp %= gw().N();
    _res64 = tmp.to_res64();
    tmp = gcd(std::move(tmp), gw().N());

    done(tmp);
}

void EdECMStage2::init(InputNum* input, GWState* gwstate, File* file, Logging* logging, arithmetic::Giant& X, arithmetic::Giant& Y, arithmetic::Giant& Z, arithmetic::Giant& T, arithmetic::Giant& EdD)
{
    Stage2::init(input, gwstate, file, read_state<State>(file), logging, _pairing.last_D - _pairing.first_D + 1);
    _state_update_period = MULS_PER_STATE_UPDATE*10/_D;
    if (_LN > 0 && _state_update_period%_LN != 0)
        _state_update_period += _LN - _state_update_period%_LN;
    _logging->set_prefix(input->display_text() + ", EdECM stage 2, ");
    _logging->info("B2 = %" PRId64, _B2);
    if (state() != nullptr)
        _logging->info(", restarting at %.1f%%", 100.0*state()->iteration()/iterations());
    _logging->info(".\n");
    _X = X;
    _Y = Y;
    _Z = Z;
    _T = T;
    _EdD = EdD;
}

void EdECMStage2::setup()
{
    if (!_ed_d)
    {
        _ed_d.reset(new GWNum(gw()));
        *_ed_d = _EdD;
    }
    if (!montgomery)
        montgomery.reset(new MontgomeryArithmetic(*_ed_d));
    montgomery->set_gw(gw());

    if (!_Pn)
        _Pn.reset(new EdY(*montgomery));
    if (!_Pn1)
        _Pn1.reset(new EdY(*montgomery));
    if (!_W)
        _W.reset(new EdY(*montgomery));

    EdwardsArithmetic ed(gw());
    EdPoint EdP(ed);
    EdP.deserialize(_X, _Y, _Z, _T);
    EdPoint EdW(ed);
    Giant tmp;
    tmp = _D;
    ed.mul(EdP, tmp, EdW);

    if (_precomp.empty())
    {
        int transforms = -(int)gw().gwdata()->fft_count;
        *_Pn = EdP;
        std::vector<std::unique_ptr<EdY>> precomp;
        int precomp_size = precompute<EdY>(*montgomery, *_Pn, *_W, *_Pn1, precomp);
        if (_L == 1 && _D%4 != 0)
            *_W = EdW;

        precomp.emplace_back(_W.release());
        try
        {
            montgomery->normalize(precomp.begin(), precomp.end());
        }
        catch (const NoInverseException& e)
        {
            done(e.divisor);
            return;
        }
        _W.reset(precomp.back().release());
        //GWASSERT((gw().popg() = *EdW.Y)%gw().N() == (gw().popg() = (*_W->Y)*(*EdW.Z))%gw().N());
        commit_setup();
        _precomp = std::move(precomp);

        transforms += (int)gw().gwdata()->fft_count;
        _transforms -= transforms;
        _logging->info("%d precomputed values (%d transforms), %d steps", precomp_size, transforms, _pairing.last_D - _pairing.first_D + 1);
        if (_LN > 0)
            _logging->info(" in batches of %d", _LN);
        _logging->info(".\n");
    }

    if (state() == nullptr)
        reset_state<State>();
    int v = state()->iteration() + _pairing.first_D;

    montgomery->init(*_Pn);
    *_Pn1 = *_W;
    if (v > 0)
    {
        tmp = v;
        ed.mul(EdW, tmp, EdP);
        *_Pn = EdP;
        montgomery->optimize(*_Pn);
        if (v > 1)
            ed.add(EdP, EdW, EdP, ed.ED_PROJECTIVE);
        else
            ed.dbl(EdP, EdP, ed.ED_PROJECTIVE);
        *_Pn1 = EdP;
        montgomery->optimize(*_Pn1);
    }
    commit_setup();
}

void EdECMStage2::release()
{
    _precomp.clear();
    _Pn.reset();
    _Pn1.reset();
    _W.reset();
    _ed_d.reset();
    montgomery.reset();
}

void EdECMStage2::execute()
{
    if (success())
        return;
    int i;
    Giant tmp;

    montgomery->set_gw(gw());

    std::vector<std::unique_ptr<EdY>> norm_bases;
    norm_bases.reserve(_LN);
    int cur_base = 0;

    GWNum G(gw());
    G = state()->G();
    GWNum TG(gw());

    int v = _pairing.first_D;
    auto it = _pairing.distances.begin();
    while (v - _pairing.first_D < state()->iteration() && it != _pairing.distances.end())
    {
        for (; *it != 0; it++);
        it++;
        v++;
    }
    v = state()->iteration() + _pairing.first_D;

    try
    {
        while (v <= _pairing.last_D || cur_base < (int)norm_bases.size())
        {
            EdY* iD = _Pn.get();
            if (_LN > 0)
            {
                if (cur_base == (int)norm_bases.size())
                {
                    for (i = 0; i < _LN && v <= _pairing.last_D; i++, v++)
                    {
                        while (i >= norm_bases.size())
                            norm_bases.emplace_back(new EdY(*montgomery));
                        if (v + 2 <= _pairing.last_D)
                        {
                            if (v > 0)
                                montgomery->add(*_W, *_Pn1, *_Pn, *norm_bases[i]);
                            else
                                montgomery->dbl(*_Pn1, *norm_bases[i]);
                        }
                        swap(*_Pn, *_Pn1);
                        swap(*_Pn1, *norm_bases[i]);
                        norm_bases[i]->ZmY.reset();
                    }
                    norm_bases.resize(i);
                    montgomery->normalize(norm_bases.begin(), norm_bases.end());
                    cur_base = 0;
                }
                iD = norm_bases[cur_base].get();
            }
            for (; *it != 0; it++)
                if (iD->Z)
                {
                    gw().mul(*iD->Z, *_precomp[*it/2]->Y, TG, GWMUL_STARTNEXTFFT);
                    gw().submul(*iD->Y, TG, G, G, GWMUL_STARTNEXTFFT_IF(!is_last(v - _pairing.first_D - ((int)norm_bases.size() - cur_base))));
                }
                else
                    gw().submul(*iD->Y, *_precomp[*it/2]->Y, G, G, GWMUL_STARTNEXTFFT_IF(!is_last(v - _pairing.first_D - ((int)norm_bases.size() - cur_base))));
            it++;
            if (_LN > 0)
                cur_base++;
            else
            {
                if (v > 0)
                    montgomery->add(*_W, *_Pn1, *_Pn, *_Pn);
                else
                    montgomery->dbl(*_Pn1, *_Pn);
                swap(*_Pn, *_Pn1);
                v++;
            }
            commit_execute<State>(v - _pairing.first_D - ((int)norm_bases.size() - cur_base), G);
            GWASSERT(state()->iteration() != v - _pairing.first_D - ((int)norm_bases.size() - cur_base) || ((int)norm_bases.size() == cur_base)); // deterministic restart
        }
        norm_bases.clear();

        tmp = G;
        tmp %= gw().N();
        _res64 = tmp.to_res64();
        tmp = gcd(std::move(tmp), gw().N());

        done(tmp);
    }
    catch (const NoInverseException& e)
    {
        done(e.divisor);
    }
}
