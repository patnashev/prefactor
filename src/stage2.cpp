
#include <deque>
#include <algorithm>

#include "gwnum.h"
#include "cpuid.h"
#include "stage2.h"
#include "exception.h"
#include "lucas.h"
#include "montgomery.h"

using namespace arithmetic;

void report_factor(const Giant& f, InputNum& input);

template<int TL>
Stage2::Pairing get_pairing_L(PrimeList& primes, int B1, int B2, int D, int A, int L, bool with_distances)
{
    int i, j, k;
    int d, p, q;
    Stage2::Pairing ret;

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
            //d = i%(D/2) + (i/(D/2) == 1 ? -1 : i/(D/2) == 2 ? 2 : i/(D/2) == 3 ? 6 : i/(D/2) == 4 ? 14 : i/(D/2) == 5 ? 30 : 0)*D;
            d = i%(D/2) + ((1 << (i/(D/2)*1)) - 1)*D;
            //d = i;
            //d = i + i/D*D;
            //d = i%D + ((1 << (i/D*1)) - 1)*D;
            //d = i%D + (i/D == 1 ? 1 : i/D == 2 ? 2 : i/D == 3 ? 8 : i/D == 4 ? 21 : i/D == 5 ? 48 : 0)*D;
            //d = i%D + (i/D == 1 ? 1 : i/D == 2 ? 2 : i/D == 3 ? 3 : i/D == 4 ? 5 : i/D == 5 ? 7 : 0)*D;
            dist.push_back(d);
        }
    for (k = dist.size(), i = 0; i < k; i++)
        dist.push_back(-dist[i]);

    std::vector<std::vector<int>> dist_rem(D);
    for (i = 1; i < D; i++)
        if (gcd(D/A, i) == 1)
        {
            for (j = 0; j < dist.size(); j++)
            {
                d = i + dist[j];
                if ((d%D != 0 && (d%D + D)%D != D/A))
                    continue;
                dist_rem[i].push_back(dist[j]*2);
            }
        }

    std::vector<int> relocs;
    for (j = 3; B1*j <= B2; j += 2)
        if (gcd(D/A, j) == 1)
            relocs.push_back(j);

    int *map = new int[B2/2];
    memset(map, 0, sizeof(int)*B2/2);

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
            srclist.back().adjacency[0] = plist.size();
            for (i = 0; i < relocs.size() && *it*relocs[i] <= B2; i++)
                if (*it*relocs[i]*relocs[0] > B2)
                {
                    map[*it*relocs[i]/2] = plist.size();
                    plist.emplace_back(*it*relocs[i]);
                    plist.back().match = 1 - srclist.size();
                }
            if (srclist.back().adjacency[0] == plist.size())
                throw ArithmeticException("Invalid B1/B2.");
        }
        else
        {
            map[*it/2] = plist.size();
            plist.emplace_back(*it);
        }
    }
    for (auto itp = plist.begin(); itp != plist.end(); itp++)
    {
        j = 0;
        p = itp->value;
        auto& p_dists = dist_rem[p%D];
        for (auto itd = p_dists.begin(); itd != p_dists.end(); itd++)
        {
            q = p + *itd;
            if (q <= B1 || q > B2 || map[q/2] == 0)
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
            p = k < plist.size() ? k : k - plist.size();
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
        
        std::sort(D_distances.begin(), D_distances.end(), [](std::pair<int, int>& a, std::pair<int, int>& b) { return a.first < b.first; });
        ret.first_D = D_distances[0].first/D;
        int cur = ret.first_D*D;
        auto it = D_distances.begin();
        while (it != D_distances.end())
        {
            for (; it != D_distances.end() && cur == it->first; it++)
                ret.distances.push_back(it->second);
            ret.distances.push_back(0);
            if (A > 1)
            {
                for (; it != D_distances.end() && cur + D/A == it->first; it++)
                    ret.distances.push_back(it->second);
                ret.distances.push_back(0);
            }
            cur += D;
        }
    }

    return ret;
}

Stage2::Pairing Stage2::get_pairing(PrimeList& primes, int B1, int B2, int D, int A, int L, bool with_distances)
{
    if (A > 1)
    {
        if (L == 1)
            return get_pairing_L<2>(primes, B1, B2, D, A, L, with_distances);
        if (L == 2)
            return get_pairing_L<4>(primes, B1, B2, D, A, L, with_distances);
        if (L == 3)
            return get_pairing_L<6>(primes, B1, B2, D, A, L, with_distances);
        if (L == 4)
            return get_pairing_L<8>(primes, B1, B2, D, A, L, with_distances);
        if (L == 5)
            return get_pairing_L<10>(primes, B1, B2, D, A, L, with_distances);
        if (L == 6)
            return get_pairing_L<12>(primes, B1, B2, D, A, L, with_distances);
        if (L == 7)
            return get_pairing_L<14>(primes, B1, B2, D, A, L, with_distances);
        if (L == 8)
            return get_pairing_L<16>(primes, B1, B2, D, A, L, with_distances);
        if (L == 9)
            return get_pairing_L<18>(primes, B1, B2, D, A, L, with_distances);
        if (L == 10)
            return get_pairing_L<20>(primes, B1, B2, D, A, L, with_distances);
        if (L == 11)
            return get_pairing_L<22>(primes, B1, B2, D, A, L, with_distances);
        if (L == 12)
            return get_pairing_L<24>(primes, B1, B2, D, A, L, with_distances);
        if (L == 13)
            return get_pairing_L<26>(primes, B1, B2, D, A, L, with_distances);
        if (L == 14)
            return get_pairing_L<28>(primes, B1, B2, D, A, L, with_distances);
        if (L == 15)
            return get_pairing_L<30>(primes, B1, B2, D, A, L, with_distances);
    }
    else
    {
        if (L == 1)
            return get_pairing_L<1>(primes, B1, B2, D, A, L, with_distances);
        if (L == 2)
            return get_pairing_L<2>(primes, B1, B2, D, A, L, with_distances);
        if (L == 3)
            return get_pairing_L<3>(primes, B1, B2, D, A, L, with_distances);
        if (L == 4)
            return get_pairing_L<4>(primes, B1, B2, D, A, L, with_distances);
        if (L == 5)
            return get_pairing_L<5>(primes, B1, B2, D, A, L, with_distances);
        if (L == 6)
            return get_pairing_L<6>(primes, B1, B2, D, A, L, with_distances);
        if (L == 7)
            return get_pairing_L<7>(primes, B1, B2, D, A, L, with_distances);
        if (L == 8)
            return get_pairing_L<8>(primes, B1, B2, D, A, L, with_distances);
        if (L == 9)
            return get_pairing_L<9>(primes, B1, B2, D, A, L, with_distances);
        if (L == 10)
            return get_pairing_L<10>(primes, B1, B2, D, A, L, with_distances);
        if (L == 11)
            return get_pairing_L<11>(primes, B1, B2, D, A, L, with_distances);
        if (L == 12)
            return get_pairing_L<12>(primes, B1, B2, D, A, L, with_distances);
        if (L == 13)
            return get_pairing_L<13>(primes, B1, B2, D, A, L, with_distances);
        if (L == 14)
            return get_pairing_L<14>(primes, B1, B2, D, A, L, with_distances);
        if (L == 15)
            return get_pairing_L<15>(primes, B1, B2, D, A, L, with_distances);
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
    for (i = 1; i < _D/2*_A; i++)
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
    if (_D/2/dist != 1)
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

// P+-1 factoring stage 2.
void PP1Stage2::run(InputNum& input, Giant& P, bool minus1)
{
    int i;
    int v;
    Giant tmp;

    if (_B2 <= _B1)
        return;
    printf("%s, P%c1 stage 2, B2 = %d.\n", input.display_text().data(), minus1 ? '-' : '+', _B2);

    double timer = getHighResTimer();
    int transforms = -(int)gw().gwdata()->fft_count;

    LucasVArithmetic lucas(gw());
    LucasV Vn(lucas);
    LucasV Vn1(lucas);
    LucasV W(lucas);
    GWNum G(gw());
    LucasV Va(lucas);
    LucasV Va1(lucas);
    std::vector<std::unique_ptr<LucasV>> precomp;

    Vn.V() = P;
    int precomp_size = precompute<LucasV>(lucas, Vn, W, Va, precomp);

    v = _pairing.first_D;

    if (_A > 1)
    {
        LucasV Wa(lucas);
        swap(Wa, Va);
        Va1 = Wa;
        if (v > 0)
            lucas.mul(Wa, v*_A, Va, Va1);
        swap(Vn1, Va);
        Vn = Vn1;
        Va = Va1;
        for (i = 0; i < _A; i++)
        {
            lucas.add(Wa, Va1, Vn1, Vn1);
            swap(Vn1, Va1);
        }
    }
    else
    {
        // V_{v*D} V_{(v+1)*D}
        lucas.init(Vn);
        Vn1 = W;
        if (v > 0)
            lucas.mul(W, v, Vn, Vn1);
    }

    transforms += (int)gw().gwdata()->fft_count;
    printf("%d precomputed values (%d transforms), %d%% pairing, D = %d, L = %d", precomp_size, transforms, (200*_pairing.pairs + _pairing.total/2)/_pairing.total, _D, _L);
    if (_A > 1)
        printf(", A = %d.\n", _A);
    else
        printf(".\n");

    transforms = -(int)gw().gwdata()->fft_count;
    G = 1;
    auto it = _pairing.distances.begin();
    while (it != _pairing.distances.end())
    {
        for (; *it != 0; it++)
            gw().submul(Vn.V(), precomp[*it/2]->V(), G, G, GWMUL_STARTNEXTFFT);
        it++;
        if (v > 0)
            lucas.add(W, Vn1, Vn, Vn);
        else
            lucas.dbl(Vn1, Vn);
        swap(Vn, Vn1);
        v++;
        if (_A > 1)
        {
            for (; *it != 0; it++)
                gw().submul(Va.V(), precomp[*it/2]->V(), G, G, GWMUL_STARTNEXTFFT);
            it++;
            lucas.add(W, Va1, Va, Va);
            swap(Va, Va1);
        }
    }

    tmp = G;
    tmp = gcd(std::move(tmp), gw().N());

    timer = (getHighResTimer() - timer)/getHighResTimerFrequency();
    transforms += (int)gw().gwdata()->fft_count;
    if (tmp == 0 || tmp == gw().N())
    {
        printf("All divisors of N < B1.\n");
    }
    else if (tmp != 1)
    {
        report_factor(tmp, input);
    }
    else
    {
        printf("No factors found, transforms: %d, time: %d s.\n", transforms, (int)timer);
    }
}

// EdECM factoring stage 2.
void EdECMStage2::run(InputNum& input, GWNum& ed_d, EdPoint& P)
{
    int i;
    int v;
    Giant tmp;

    if (_B2 <= _B1)
        return;
    printf("%s, EdECM stage 2, B2 = %d.\n", input.display_text().data(), _B2);
    
    double timer = getHighResTimer();
    int transforms = -(int)gw().gwdata()->fft_count;

    MontgomeryArithmetic montgomery(gw(), ed_d);
    EdY Pn(montgomery);
    EdY Pn1(montgomery);
    EdY W(montgomery);
    GWNum G(gw());
    GWNum TG(gw());
    std::vector<std::unique_ptr<EdY>> precomp;

    Pn = P;
    int precomp_size = precompute<EdY>(montgomery, Pn, W, Pn1, precomp);
    precomp.emplace_back(&W);

    try
    {
        montgomery.normalize(precomp.begin(), precomp.end());

        v = _pairing.first_D;

        // V_{v*D} V_{(v+1)*D}
        montgomery.init(Pn);
        Pn1 = W;
        if (v > 0)
        {
            EdwardsArithmetic edwards = P.arithmetic();
            EdPoint EdW(edwards);
            tmp = _D;
            edwards.mul(P, tmp, EdW);
            //GWASSERT(*EdW.Y == *W.Y*(*EdW.Z));
            tmp = v;
            EdPoint EdPn(edwards);
            edwards.mul(EdW, tmp, EdPn);
            Pn = EdPn;
            montgomery.optimize(Pn);
            if (v > 1)
                edwards.add(EdPn, EdW, EdPn, edwards.ED_PROJECTIVE);
            else
                edwards.dbl(EdPn, EdPn, edwards.ED_PROJECTIVE);
            Pn1 = EdPn;
            montgomery.optimize(Pn1);
        }
        std::deque<std::unique_ptr<EdY>> norm_bases;

        transforms += (int)gw().gwdata()->fft_count;
        printf("%d precomputed values (%d transforms), %d%% pairing, D = %d, L = %d", precomp_size, transforms, (200*_pairing.pairs + _pairing.total/2)/_pairing.total, _D, _L);
        if (_LN > 0)
            printf(", LN = %d.\n", _LN);
        else
            printf(".\n");

        transforms = -(int)gw().gwdata()->fft_count;
        G = 1;
        auto it = _pairing.distances.begin();
        while (it != _pairing.distances.end())
        {
            EdY* iD = &Pn;
            if (_LN > 0)
            {
                if (norm_bases.empty())
                {
                    for (i = 0; i < _LN && (v - 1)*_D < _B2; i++, v++)
                    {
                        norm_bases.emplace_back(new EdY(montgomery));
                        if ((v + 1)*_D < _B2)
                            if (v > 0)
                                montgomery.add(W, Pn1, Pn, *norm_bases.back());
                            else
                                montgomery.dbl(Pn1, *norm_bases.back());
                        swap(Pn, Pn1);
                        swap(Pn1, *norm_bases.back());
                    }
                    montgomery.normalize(norm_bases.begin(), norm_bases.end());
                }
                iD = norm_bases.front().get();
            }
            for (; *it != 0; it++)
                if (iD->Z)
                {
                    gw().mul(*iD->Z, *precomp[*it/2]->Y, TG, GWMUL_STARTNEXTFFT);
                    gw().submul(*iD->Y, TG, G, G, GWMUL_STARTNEXTFFT);
                }
                else
                    gw().submul(*iD->Y, *precomp[*it/2]->Y, G, G, GWMUL_STARTNEXTFFT);
            it++;
            if (_LN > 0)
                norm_bases.pop_front();
            else
            {
                if (v > 0)
                    montgomery.add(W, Pn1, Pn, Pn);
                else
                    montgomery.dbl(Pn1, Pn);
                swap(Pn, Pn1);
                v++;
            }
        }

        tmp = G;
        tmp = gcd(std::move(tmp), gw().N());
    }
    catch (const NoInverseException& e)
    {
        tmp = e.divisor;
    }
    precomp.back().release();

    timer = (getHighResTimer() - timer)/getHighResTimerFrequency();
    transforms += (int)gw().gwdata()->fft_count;
    if (tmp == 0 || tmp == gw().N())
    {
        printf("All divisors of N < B1.\n");
    }
    else if (tmp != 1)
    {
        report_factor(tmp, input);
    }
    else
    {
        printf("No factors found, transforms: %d, time: %d s.\n", transforms, (int)timer);
    }
}
