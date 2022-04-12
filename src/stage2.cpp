
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

template<int TL>
Stage2::Pairing get_pairing_L(Logging& logging, PrimeList& primes, int B1, int B2, int D, int A, int L, bool with_distances)
{
    int i, j, k;
    int d, p, q;
    Stage2::Pairing ret;

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

    std::vector<Stage2::prime<TL>> plist;
    std::vector<Stage2::prime<1>> srclist;
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

Stage2::Pairing Stage2::get_pairing(Logging& logging, PrimeList& primes, int B1, int B2, int D, int A, int L, bool with_distances)
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
int Stage2::precompute(DifferentialGroupArithmetic<Element>& arithmetic, Element& X1, Element& XD, Element& XDA, std::vector<std::unique_ptr<Element>>& precomp)
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
        /*XD = X1;
        XD *= _D/_A;
        GWASSERT(XD == XDA);*/
    }
    // V_D
    swap(XD, Xn);
    if (_L == 1 && _D%4 != 0)
    {
        arithmetic.init(XD);
        for (j = 0; j < _D/4; j++)
            if (precomp[j] && gcd(2*j + 1, _D) != 1)
            {
                precomp[j].reset();
                precomp_size--;
            }
    }
    else if (_D/2/dist != 1)
    {
        int pd = _D/2/dist;
        for (j = 0; (pd & 1) == 0; pd >>= 1, j++);
        if (pd > 1)
            arithmetic.add(*precomp[dist*(pd - 1)/2 - 1], *precomp[dist*(pd + 1)/2], Xn1, XD);
        for (; j > 0; j--)
            arithmetic.dbl(XD, XD);
    }
    /*Xn = X1;
    Xn *= _D;
    GWASSERT(Xn == XD);*/

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
    /*for (i = 1; i < _D/2*_L; i++)
        if (gcd(_D/_A, i%(_D/2)) == 1)
        {
            int d = i%(_D/2) + ((1 << (i/(_D/2)*1)) - 1)*_D;
            Xn = X1;
            Xn *= d;
            GWASSERT(Xn == *precomp[d/2]);
        }*/

    return precomp_size;
}

void Stage2::poly_init()
{
    int i;

    _poly_mult.reserve(poly_power() + 1);
    _poly_mult.emplace_back(gw());
    for (i = 1; i <= poly_power(); i++)
    {
        if (2*(1 << i) > _poly_mult[i - 1].max_output())
        {
            _poly_gwstate.emplace_back();
            if (_poly_gwstate.size() == 1)
                _poly_gwstate.back().copy(gw().state());
            else
                _poly_gwstate.back().copy(_poly_gwstate[_poly_gwstate.size() - 2]);
            _poly_gwstate.back().next_fft_count++;
            _input->setup(_poly_gwstate.back());
            _poly_gw.emplace_back(_poly_gwstate.back());
            _poly_mult.emplace_back(_poly_gw.back());
        }
        else
            _poly_mult.emplace_back(_poly_mult[i - 1].gw());
        //if (i > 5)
        //polymult_set_num_threads(_poly_mult[i].pmdata(), 2);
    }
}

// http://cr.yp.to/arith/scaledmod-20040820.pdf
void Stage2::poly_setup(std::vector<GWNum*>& roots)
{
    int i, j;
    double timer = getHighResTimer();

    poly_init();
    if (!_poly_mod)
    {
        _poly_mod_container.resize(poly_power() + 1);
        _poly_mod = _poly_mod_container.data();
    }

    GWASSERT((int)roots.size() <= (1 << poly_power()));
    for (i = 0; i < (1 << poly_power()); i++)
    {
        _poly_mod[0].emplace_back(_poly_mult[0], 0, true);
        if (i < roots.size())
            _poly_mult[0].init(std::move(*roots[i]), true, _poly_mod[0][i]);
    }
    for (j = 1; j <= poly_power(); j++)
        if (&_poly_mult[j].gw() != &_poly_mult[j - 1].gw())
        {
            Poly tmp(_poly_mult[j - 1], 0, true);
            for (i = 0; i < (1 << (poly_power() - j)); i++)
            {
                _poly_mult[j - 1].preprocess_and_mul(_poly_mod[j - 1][2*i], _poly_mod[j - 1][2*i + 1], tmp, 1 << j, 0);
                _poly_mod[j].emplace_back(_poly_mult[j], 0, true);
                _poly_mult[j - 1].convert(tmp, _poly_mod[j][i]);
            }
        }
        else
        {
            for (i = 0; i < (1 << (poly_power() - j)); i++)
            {
                _poly_mod[j].emplace_back(_poly_mult[j], 0, true);
                _poly_mult[j - 1].preprocess_and_mul(_poly_mod[j - 1][2*i], _poly_mod[j - 1][2*i + 1], _poly_mod[j][i], 1 << j, POLYMULT_STARTNEXTFFT);
            }
        }
    for (i = 0; i < roots.size(); i++)
        gwfft(gw().gwdata(), *_poly_mod[0][i].data(), *_poly_mod[0][i].data());
#ifdef _DEBUG
/*    std::vector<GWNum> polyTest;
    for (i = 0; i < (1 << poly_power()) && _poly_mod[0][i].degree() > 0; i++)
    {
        polyTest.emplace_back(gw());
        GWNum a(gw());
        a = 0;
        gw().sub(a, *roots[i], a, GWADD_GUARANTEED_OK);
        polyTest.back() = _poly_mod[poly_power()][0].eval(a);
        GWASSERT(polyTest.back() == 0);
    }*/
#endif
    Poly polyR = _poly_mod[poly_power()][0].reciprocal(poly_degree(), POLYMULT_STARTNEXTFFT);
    int degree = (1 << (poly_power() + 1)) - _poly_mod[poly_power()][0].degree();
    if (polyR.degree() < degree)
        polyR <<= degree - polyR.degree();
#ifdef _DEBUG
/*    Poly polyT(_poly_mult[0], 1 << poly_power(), false);
    _poly_mult[0].mul_half(polyR, _poly_mod[poly_power()][0], polyT, 0);
    for (i = 0; i < _LN; i++)
        GWASSERT(polyT.at((1 << poly_power()) - i - 1) == 0);*/
#endif
    _poly_mult[poly_power()].preprocess(polyR, 1 << (poly_power() + 1));
    _poly_mod[poly_power()][0] = std::move(polyR);

    commit_setup();
    timer = (getHighResTimer() - timer)/getHighResTimerFrequency();
    _logging->info("polynomial mode, D=%d, degree %d, setup time: %.3f s.\n", _D, poly_degree(), timer);
}

void Stage2::poly_release()
{
    _poly_mod = nullptr;
    _poly_mod_container.clear();
    _poly_prod.clear();
    _poly_mult.clear();
    _poly_gw.clear();
    _poly_gwstate.clear();
}

void Stage2::poly_execute(std::vector<GWNum>& roots, GWNum& G)
{
    int i, j;
    for (i = 0; i < _poly_mult.size(); i++)
        _poly_mult[i].gw().gwdata()->gwnum_max_free_count = 2*(1 << poly_power());

    GWASSERT((int)roots.size() <= (1 << poly_power()));
    if (_poly_prod.size() < poly_power() + 1)
        _poly_prod.resize(poly_power() + 1);
    for (i = 0; i < (1 << poly_power()); i++)
    {
        while (_poly_prod[0].size() <= i)
            _poly_prod[0].emplace_back(_poly_mult[0], 0, true);
        if (i < roots.size())
#ifdef _DEBUG
            _poly_mult[0].init(roots[i], true, _poly_prod[0][i]);
#else
            _poly_mult[0].init(std::move(roots[i]), true, _poly_prod[0][i]);
#endif
        else
            _poly_mult[0].init(true, _poly_prod[0][i]);
    }
    for (j = 1; j <= poly_power(); j++)
        if (&_poly_mult[j].gw() != &_poly_mult[j - 1].gw())
        {
            Poly tmp(_poly_mult[j - 1], 0, true);
            for (i = 0; i < (1 << (poly_power() - j)); i++)
            {
                _poly_mult[j - 1].mul(std::move(_poly_prod[j - 1][2*i]), std::move(_poly_prod[j - 1][2*i + 1]), tmp, 0);
                while (_poly_prod[j].size() <= i)
                    _poly_prod[j].emplace_back(_poly_mult[j], 0, true);
                _poly_mult[j - 1].convert(tmp, _poly_prod[j][i]);
            }
        }
        else
        {
            for (i = 0; i < (1 << (poly_power() - j)); i++)
            {
                while (_poly_prod[j].size() <= i)
                    _poly_prod[j].emplace_back(_poly_mult[j], 0, true);
                _poly_mult[j - 1].mul(std::move(_poly_prod[j - 1][2*i]), std::move(_poly_prod[j - 1][2*i + 1]), _poly_prod[j][i], POLYMULT_STARTNEXTFFT);
            }
        }
#ifdef _DEBUG
/*    for (i = 0; i < roots.size(); i++)
    {
        GWNum a(gw());
        a = 0;
        gw().sub(a, *roots[i], a, GWADD_GUARANTEED_OK);
        GWNum b(_poly_mult[pm].gw());
        b = a;
        GWASSERT(_poly_prod[poly_power()][0].eval(b) == 0);
    }*/
#endif

    _poly_mult[poly_power()].mul_half_preprocessed(std::move(_poly_prod[poly_power()][0]), _poly_mod[poly_power()][0], _poly_prod[poly_power()][0], 1 << poly_power(), &_poly_mod[poly_power()][0].pm().gw() == &_poly_mod[poly_power() - 1][0].pm().gw() ? POLYMULT_STARTNEXTFFT : 0);

    for (j = poly_power() - 1; j >= 0; j--)
    {
        PolyMult* prev = &_poly_mult[j + 1];
        PolyMult* cur = &_poly_mult[j];
        PolyMult* next = j > 0 ? &_poly_mult[j - 1] : cur;
        bool convert = &prev->gw() != &cur->gw();
        int options = !convert ? POLYMULT_STARTNEXTFFT : 0;

        for (i = 0; i < (1 << (poly_power() - j)); i += 2)
        {
            if (_poly_mod[j][i + 1].degree() == 0 && _poly_mod[j][i].degree() == 0)
            {
                prev->free(_poly_prod[j + 1][i/2]);
                continue;
            }
            if (convert)
            {
                Poly tmp(*cur, 0, false);
                prev->convert(_poly_prod[j + 1][i/2], tmp);
                prev->free(_poly_prod[j + 1][i/2]);

                /*cur->alloc(_poly_prod[j][i], 1 << j);
                cur->mul_half(tmp, _poly_mod[j][i + 1], _poly_prod[j][i], options);
                cur->alloc(_poly_prod[j][i + 1], 1 << j);
                cur->mul_half(tmp, _poly_mod[j][i], _poly_prod[j][i + 1], options);*/
                if (_poly_mod[j][i + 1].degree() == 0)
                    cur->mul_half_preprocessed(std::move(tmp), _poly_mod[j][i + 1], _poly_prod[j][i], 1 << j, options);
                else
                    cur->mul_half_preprocessed(std::move(tmp), _poly_mod[j][i + 1], _poly_mod[j][i], _poly_prod[j][i], _poly_prod[j][i + 1], 1 << j, options);
            }
            else
            {
                /*cur->alloc(_poly_prod[j][i], 1 << j);
                cur->mul_half(_poly_prod[j + 1][i/2], _poly_mod[j][i + 1], _poly_prod[j][i], options);
                cur->alloc(_poly_prod[j][i + 1], 1 << j);
                cur->mul_half(_poly_prod[j + 1][i/2], _poly_mod[j][i], _poly_prod[j][i + 1], options);
                prev->free(_poly_prod[j + 1][i/2]);*/
                if (_poly_mod[j][i + 1].degree() == 0)
                    cur->mul_half_preprocessed(std::move(_poly_prod[j + 1][i/2]), _poly_mod[j][i + 1], _poly_prod[j][i], 1 << j, options);
                else
                    cur->mul_half_preprocessed(std::move(_poly_prod[j + 1][i/2]), _poly_mod[j][i + 1], _poly_mod[j][i], _poly_prod[j][i], _poly_prod[j][i + 1], 1 << j, options);
            }
        }
    }

#ifdef _DEBUG
    GWNum polyTest(gw());
    for (i = 0; i < (1 << poly_power()) && _poly_mod[0][i].degree() > 0; i++)
    {
        polyTest = 1;
        GWNumWrapper a = _poly_mod[0][i].at(0);
        for (int j = 0; j < roots.size(); j++)
            polyTest *= roots[j] - a;
        GWASSERT(polyTest == _poly_prod[0][i].at(0));
    }
#endif

    for (i = (1 << poly_power()) - 1; i >= 0; i--)
        if (_poly_mod[0][i].degree() > 0)
        {
            GWNumWrapper a = _poly_prod[0][i].at(0);
            gw().mul(a, G, G, GWMUL_STARTNEXTFFT_IF(i > 0));
        }
#if !_DEBUG
    for (i = 0; i < roots.size(); i++)
        if (_poly_prod[0][i].size() > 0)
            roots[i] = std::move(_poly_prod[0][i].pop_back());
        else
            roots[i].arithmetic().alloc(roots[i]);
#endif
/*
    for (i = 0; i < roots.size(); i++)
        if (_poly_prod[poly_power()][0].size() > 0)
            *roots[i] = std::move(_poly_prod[poly_power()][0].pop_back());
        else
            roots[i]->arithmetic().alloc(*roots[i]);*/
}

void poly_thread(void* data)
{
    Stage2* stage2 = (Stage2*)data;
    try
    {
        stage2->run();
    }
    catch (const std::exception&)
    {
    }
    gwmutex_lock(&stage2->_poly_thread_main->_poly_mutex);
    int num_active = --stage2->_poly_thread_main->_poly_threads_active;
    gwmutex_unlock(&stage2->_poly_thread_main->_poly_mutex);
    if (num_active == 0)
        gwevent_signal(&stage2->_poly_thread_main->_poly_done);
}

void Stage2::poly_threads_init()
{
    gwmutex_init(&_poly_mutex);
    gwevent_init(&_poly_done);
    _poly_G.reset(new GWNum(gw()));
    *_poly_G = 1;
    _poly_threads_active = _poly_thread_helpers.size();
    for (auto it = _poly_thread_helpers.begin(); it != _poly_thread_helpers.end(); it++)
        gwthread_create_waitable(&it->id, &poly_thread, (void *)it->stage2.get());
}

void Stage2::poly_helper_done(arithmetic::GWNum& G)
{
    gwmutex_lock(&_poly_thread_main->_poly_mutex);
    gw().mul(G, *_poly_thread_main->_poly_G, *_poly_thread_main->_poly_G, GWMUL_STARTNEXTFFT);
    gwmutex_unlock(&_poly_thread_main->_poly_mutex);
}

void Stage2::poly_threads_wait(arithmetic::GWNum& G)
{
    gwevent_wait(&_poly_done, 0);
    gw().mul(*_poly_G, G, G, 0);
    _poly_G.reset();
    gwmutex_destroy(&_poly_mutex);
    gwevent_destroy(&_poly_done);
}

void Stage2::init(InputNum* input, GWState* gwstate, File* file, TaskState* state, Logging* logging)
{
    while (is_poly() && gwstate->max_polymult_output() < 2*(1 << poly_power()))
    {
        gwstate->done();
        gwstate->next_fft_count++;
        input->setup(*gwstate);
        logging->warning("Switching to %s\n", gwstate->fft_description.data());
    }
    if (is_poly_threaded())
    {
        int Ds_per_thread = (_pairing.last_D - _pairing.first_D + _poly_threads)/_poly_threads;
        int Ds_per_thread_last = _pairing.last_D - _pairing.first_D + 1 - Ds_per_thread*(_poly_threads - 1);
        _pairing.last_D = _pairing.first_D + Ds_per_thread - 1;
        for (int i = 0; i < _poly_thread_helpers.size(); i++)
        {
            _poly_thread_helpers[i].stage2->_pairing.first_D = _pairing.first_D + Ds_per_thread*(i + 1);
            _poly_thread_helpers[i].stage2->_pairing.last_D = _poly_thread_helpers[i].stage2->_pairing.first_D + (i < _poly_threads - 2 ? Ds_per_thread : Ds_per_thread_last) - 1;
            _poly_thread_helpers[i].gwstate.reset(new GWState());
            _poly_thread_helpers[i].gwstate->copy(*gwstate);
            input->setup(*_poly_thread_helpers[i].gwstate);
            _poly_thread_helpers[i].file = file->add_child(std::to_string(i + 2), gwstate->fingerprint);
        }
    }

    Task::init(gwstate, file, state, logging, _pairing.last_D - _pairing.first_D + 1);
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
}

void Stage2::done(const arithmetic::Giant& factor)
{
    _timer = (getHighResTimer() - _timer)/getHighResTimerFrequency();
    _transforms += (int)_gwstate->handle.fft_count;
    _logging->progress().update(1, (int)_gwstate->handle.fft_count/2);
    _logging->info("transforms: %d, time: %.1f s.\n", _transforms, _timer);
    if (factor == 0 || factor == *_gwstate->N)
    {
        _logging->warning("all divisors less than B1.\n");
        _logging->result("all divisors less than B1, time: %.1f s.\n", _logging->progress().time_total());
        _success = true;
    }
    else if (factor != 1)
    {
        _logging->report_factor(*_input, factor);
        _success = true;
    }
    else
    {
        //_logging->info("No factors found.\n");
    }
    _logging->set_prefix("");

    if (is_poly_threaded())
    {
        for (auto it = _poly_thread_helpers.begin(); it != _poly_thread_helpers.end(); it++)
            it->file->clear();
        _poly_thread_helpers.clear();
    }
}

void PP1Stage2::init(InputNum* input, GWState* gwstate, File* file, Logging* logging, Giant& P, bool minus1)
{
    Stage2::init(input, gwstate, file, read_state<State>(file), logging);
    _state_update_period = MULS_PER_STATE_UPDATE/10;
    if (!is_poly_helper())
    {
        _logging->set_prefix(input->display_text() + (minus1 ? ", P-1 stage 2, " : ", P+1 stage 2, "));
        _logging->info("B2 = %" PRId64, _B2);
        if (is_poly())
            _logging->info(", suggested B2 = %" PRId64 "", (_pairing.first_D + (iterations() + poly_degree() - 1 - (iterations() + poly_degree() - 1)%poly_degree())*_poly_threads - 1)*_D + _D/2 - 1);
        if (state() != nullptr)
            _logging->info(", restarting at %.1f%%", 100.0*state()->iteration()/iterations());
        _logging->info(".\n");
        if (_poly_threads > 1)
            _logging->info("%d threads.\n", _poly_threads);
    }
    lucas.reset(new LucasVArithmetic());
    _P = P;

    for (auto it = _poly_thread_helpers.begin(); it != _poly_thread_helpers.end(); it++)
        static_cast<PP1Stage2*>(it->stage2.get())->init(input, it->gwstate.get(), it->file, logging, P, minus1);
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

    if (is_poly_helper())
    {
        //*_W = *static_cast<PP1Stage2*>(_poly_thread_main)->_W;
        _Vn->V() = _P;
        lucas->mul(*_Vn, _D, *_W, *_Vn1);
    }
    else if (_precomp.empty())
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

    if (is_poly_helper() && _poly_mod == nullptr)
    {
        poly_init();
        _poly_mod = static_cast<PP1Stage2*>(_poly_thread_main)->_poly_mod;
    }
    else if (is_poly() && _poly_mod == nullptr)
    {
        std::vector<GWNum*> poly_r;
        for (auto it = _precomp.begin(); it != _precomp.end(); it++)
            if (*it)
                poly_r.push_back(&(*it)->V());
        poly_setup(poly_r);
    }
}

void PP1Stage2::release()
{
    if (is_poly())
        poly_release();
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

    if (is_poly_threaded())
        poly_threads_init();

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

    if (!is_poly())
    {
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
    }
    else
    {
        int i;
        std::vector<GWNum> bases;
        bases.reserve(poly_degree());
        for (i = 0; i < poly_degree(); i++)
            bases.emplace_back(gw());
        LucasV Vtmp(*lucas);

        while (v <= _pairing.last_D)
        {
            for (i = 0; i < poly_degree() && v <= _pairing.last_D; i++, v++)
            {
                if (v + 2 <= _pairing.last_D)
                {
                    if (v > 0)
                        lucas->add(*_W, *_Vn1, *_Vn, Vtmp);
                    else
                        lucas->dbl(*_Vn1, Vtmp);
                }
                swap(*_Vn, *_Vn1);
                swap(*_Vn1, Vtmp);
                swap(Vtmp.V(), bases[i]);
            }
            bases.erase(bases.begin() + i, bases.end());
            poly_execute(bases, G);
            commit_execute<State>(v - _pairing.first_D, G);
        }
    }

    if (is_poly_helper())
        poly_helper_done(G);
    else
    {
        if (is_poly_threaded())
            poly_threads_wait(G);

        tmp = G;
        tmp = gcd(std::move(tmp), gw().N());

        done(tmp);
    }
}

void EdECMStage2::init(InputNum* input, GWState* gwstate, File* file, Logging* logging, arithmetic::Giant& X, arithmetic::Giant& Y, arithmetic::Giant& Z, arithmetic::Giant& T, arithmetic::Giant& EdD)
{
    Stage2::init(input, gwstate, file, read_state<State>(file), logging);
    if (is_poly())
        _state_update_period = poly_degree();
    else
    {
        _state_update_period = MULS_PER_STATE_UPDATE*10/_D;
        if (_LN > 0 && _state_update_period%_LN != 0)
            _state_update_period += _LN - _state_update_period%_LN;
    }
    if (!is_poly_helper())
    {
        _logging->set_prefix(input->display_text() + ", EdECM stage 2, ");
        _logging->info("B2 = %" PRId64, _B2);
        if (is_poly())
            _logging->info(", suggested B2 = %" PRId64 "", (_pairing.first_D + (iterations() + poly_degree() - 1 - (iterations() + poly_degree() - 1)%poly_degree())*_poly_threads - 1)*_D + _D/2 - 1);
        if (state() != nullptr)
            _logging->info(", restarting at %.1f%%", 100.0*state()->iteration()/iterations());
        _logging->info(".\n");
        if (_poly_threads > 1)
            _logging->info("%d threads.\n", _poly_threads);
    }
    _X = X;
    _Y = Y;
    _Z = Z;
    _T = T;
    _EdD = EdD;

    for (auto it = _poly_thread_helpers.begin(); it != _poly_thread_helpers.end(); it++)
        static_cast<EdECMStage2*>(it->stage2.get())->init(input, it->gwstate.get(), it->file, logging, X, Y, Z, T, EdD);
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

    if (is_poly_helper())
    {
        *_W = EdW;
    }
    else if (_precomp.empty())
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

    if (is_poly_helper() && _poly_mod == nullptr)
    {
        poly_init();
        _poly_mod = static_cast<EdECMStage2*>(_poly_thread_main)->_poly_mod;
    }
    else if (is_poly() && _poly_mod == nullptr)
    {
        std::vector<GWNum*> poly_r;
        for (auto it = _precomp.begin(); it != _precomp.end(); it++)
            if (*it)
                poly_r.push_back((*it)->Y.get());
        poly_setup(poly_r);
    }
}

void EdECMStage2::release()
{
    if (is_poly())
        poly_release();
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

    if (is_poly_threaded())
        poly_threads_init();

    std::vector<GWNum> poly_bases;
    poly_bases.reserve(poly_degree());
    int cur_degree = 0;

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
                    for (i = 0; i < _LN && (poly_degree() == 0 || cur_degree + i < poly_degree()) && v <= _pairing.last_D; i++, v++)
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
            if (it != _pairing.distances.end())
            {
                for (; *it != 0; it++)
                    if (iD->Z)
                    {
                        gw().mul(*iD->Z, *_precomp[*it/2]->Y, TG, GWMUL_STARTNEXTFFT);
                        gw().submul(*iD->Y, TG, G, G, GWMUL_STARTNEXTFFT_IF(!is_last(v - _pairing.first_D - ((int)norm_bases.size() - cur_base))));
                    }
                    else
                        gw().submul(*iD->Y, *_precomp[*it/2]->Y, G, G, GWMUL_STARTNEXTFFT_IF(!is_last(v - _pairing.first_D - ((int)norm_bases.size() - cur_base))));
                it++;
            }
            if (is_poly())
            {
                while (cur_degree >= poly_bases.size())
                    poly_bases.emplace_back(gw());
                if (iD->Z)
                {
                    poly_bases[cur_degree] = *iD->Z;
                    poly_bases[cur_degree].inv();
                    gw().mul(*iD->Y, poly_bases[cur_degree], poly_bases[cur_degree], GWMUL_STARTNEXTFFT);
                }
                else
                    swap(poly_bases[cur_degree], *iD->Y);
                cur_degree++;
            }
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
            if (is_poly() && (cur_degree == poly_degree() || (v > _pairing.last_D && cur_base == (int)norm_bases.size())))
            {
                poly_bases.erase(poly_bases.begin() + cur_degree, poly_bases.end());
                poly_execute(poly_bases, G);
                cur_degree = 0;
            }
            commit_execute<State>(v - _pairing.first_D - ((int)norm_bases.size() - cur_base), G);
            GWASSERT(state()->iteration() != v - _pairing.first_D - ((int)norm_bases.size() - cur_base) || ((int)norm_bases.size() == cur_base)); // deterministic restart
        }

        if (is_poly_helper())
            poly_helper_done(G);
        else
        {
            if (is_poly_threaded())
                poly_threads_wait(G);

            tmp = G;
            tmp = gcd(std::move(tmp), gw().N());

            done(tmp);
        }
    }
    catch (const NoInverseException& e)
    {
        if (is_poly_helper())
        {
            G = e.divisor;
            poly_helper_done(G);
        }
        else
            done(e.divisor);
    }
}
