
#include <cmath>
#include "group.h"
using namespace arithmetic;

int get_stage2_size(int D, int A, int L)
{
    int i;
    int size = 0;

    size += 4;
    if (A > 1)
        size += 2;

    int precomp = 0;
    for (i = 1; i < D/2; i++)
        if (gcd(D/A, i) == 1)
            precomp++;
    size += precomp*L;

    return size;
}

int get_stage2_cost_best(int B1, int B2, int D, int A, int L)
{
    int i, j;
    int stage2cost = 0;

    double num_primes = std::expint(log(B2)) - std::expint(log(B1));
    stage2cost += (int)(num_primes/2);

    PrimeList primes(B2/B1 + 100);
    for (PrimeIterator it = primes.begin(); *it*B1 < B2; it++)
        if (D/A%(*it) != 0)
        {
            B1 = B2/(*it);
            break;
        }

    if (A > 1)
    {
        D /= A;
    }

    stage2cost++;
    stage2cost++;
    int dist = 1;
    int precomp = 1;
    std::vector<char> precomp_V(D*A/2, 0);
    precomp_V[0] = 1;
    //stage2cost += D/2*L/2;
    for (i = 1; i < D*A/4; i++)
    {
        if (gcd(2*i + 1, 2*dist) != 1)
            continue;
        stage2cost++;
        precomp++;
        precomp_V[i] = 1;
        if ((i + 1)%dist == 0 && D%(i + 1) == 0 && gcd((i + 1)/dist, 2*dist) == 1)
        {
            stage2cost++;
            stage2cost++;
            dist = i + 1;
            for (j = 0; j < i; j++)
                if (precomp_V[j] != 0 && gcd(2*j + 1, 2*dist) != 1)
                {
                    precomp_V[j] = 0;
                    precomp--;
                }
        }
    }

    if (A > 1)
    {
        if (D/2/dist != 1)
        {
            int pd = D/2/dist;
            for (j = 0; (pd & 1) == 0; pd >>= 1, j++);
            if (pd > 1)
                stage2cost++;
            stage2cost += j;
        }
        D *= A;
    }
    if (D/2/dist != 1)
    {
        int pd = D/2/dist;
        for (j = 0; (pd & 1) == 0; pd >>= 1, j++);
        if (pd > 1)
            stage2cost++;
        stage2cost += j;
    }
    stage2cost += (precomp + 1)*L;

    int v = B1/D;

    if (A > 1)
    {
        if (v > 0)
            stage2cost += 2*(int)log2(v*A);
        stage2cost += A;
    }
    else
    {
        if (v > 0)
            stage2cost += 2*(int)log2(v);
    }

    stage2cost += (B2 - B1)/D;
    if (A > 1)
        stage2cost += (B2 - B1)/D;

    return stage2cost;
}

void get_stage2_cost_and_pairing(PrimeList& primes, int B1, int B2, int D, int A, int L, int *cost, double *pairing)
{
    int i, j;
    int stage2cost = 0;

    Logging logging(Logging::LEVEL_WARNING);
    Stage2::Pairing PP = Stage2::get_pairing(logging, primes, B1, B2, D, A, L, false);
    stage2cost += PP.total - PP.pairs;
    *pairing = 2*PP.pairs/(double)PP.total;

    if (A > 1)
    {
        D /= A;
    }

    stage2cost++;
    stage2cost++;
    int dist = 1;
    int precomp = 1;
    std::vector<char> precomp_V(D*A/2, 0);
    precomp_V[0] = 1;
    //stage2cost += D/2*L/2;
    for (i = 1; i < D*A/4; i++)
    {
        if (gcd(2*i + 1, 2*dist) != 1)
            continue;
        stage2cost++;
        precomp++;
        precomp_V[i] = 1;
        if ((i + 1)%dist == 0 && D%(i + 1) == 0 && gcd((i + 1)/dist, 2*dist) == 1)
        {
            stage2cost++;
            stage2cost++;
            dist = i + 1;
            for (j = 0; j < i; j++)
                if (precomp_V[j] != 0 && gcd(2*j + 1, 2*dist) != 1)
                {
                    precomp_V[j] = 0;
                    precomp--;
                }
        }
    }

    if (A > 1)
    {
        if (D/2/dist != 1)
        {
            int pd = D/2/dist;
            for (j = 0; (pd & 1) == 0; pd >>= 1, j++);
            if (pd > 1)
                stage2cost++;
            stage2cost += j;
        }
        D *= A;
    }
    if (D/2/dist != 1)
    {
        int pd = D/2/dist;
        for (j = 0; (pd & 1) == 0; pd >>= 1, j++);
        if (pd > 1)
            stage2cost++;
        stage2cost += j;
    }
    stage2cost += (precomp + 1)*L;

    int v = PP.first_D;

    if (A > 1)
    {
        if (v > 0)
            stage2cost += 2*(int)log2(v*A);
        stage2cost += A;
    }
    else
    {
        if (v > 0)
            stage2cost += 2*(int)log2(v);
    }

    stage2cost += PP.last_D - PP.first_D;
    if (A > 1)
        stage2cost += PP.last_D - PP.first_D;

    *cost = stage2cost;
}

#include <omp.h>

extern const short precomputed_stage2_params[];
extern const int precomputed_stage2_params_len;

void do_stage2_params(int B1, int B2)
{
    int i, j, k;
    std::vector<std::array<int,6>> best;
    
    PrimeList primes(B2 + 100);
    int totalp;
    for (totalp = 0; primes[totalp] <= B1; totalp++);
    int num_primes;
    for (num_primes = 0; primes[totalp] <= B2; totalp++, num_primes++);
    totalp++;

#ifndef _DEBUG
    omp_set_num_threads(1);
#endif

    std::vector<std::array<int,6>> savedCosts;
    std::ifstream fIn("costs." + std::to_string(B1/1000) + ".txt");
    fIn >> i >> j;
    while (!fIn.fail())
    {
        std::array<int, 6> tmp;
        fIn >> tmp[0] >> tmp[1] >> tmp[2] >> tmp[3] >> tmp[4] >> tmp[5];
        if (i == B1 && j == B2)
            savedCosts.push_back(tmp);
        fIn >> i >> j;
    }
    fIn.close();

    std::ofstream fOut("costs." + std::to_string(B1/1000) + ".txt", std::fstream::app);

    std::vector<std::array<int, 3>> params;
    for (i = 0; i < precomputed_stage2_params_len - 7; i += 7)
    {
        for (j = 0; j < params.size() && (precomputed_stage2_params[i + 3] != params[j][0] || precomputed_stage2_params[i + 4] != params[j][1] || precomputed_stage2_params[i + 5] != params[j][2]); j++);
        if (j >= params.size())
            params.push_back({ precomputed_stage2_params[i + 3], precomputed_stage2_params[i + 4], precomputed_stage2_params[i + 5]});
    }
    {
        for (int D = 30; D < 2150; D += 6)
            for (int iA = 0; iA < 10; iA++)
            {
                int A = 1;
                if (iA > 0)
                    A = primes[iA];
                if (D%A != 0 || A == 2)
                    continue;

                for (int L = 1; L <= 15; L++)
                {
                    for (j = 0; j < params.size() && (D != params[j][0] || A != params[j][1] || L != params[j][2]); j++);
                    if (j >= params.size())
                        params.push_back({D, A, L});
                }
            }
    }

    int it;
    #pragma omp parallel for schedule(dynamic) default(none) private(i) shared(best,totalp,savedCosts)
    for (it = 0; it < params.size(); it++)
    {
        int D = params[it][0];
        int A = params[it][1];
        int L = params[it][2];

        int size = get_stage2_size(D, A, L);
        if (size > 1000)
            continue;
        double pairing = 1.0;
        int cost = get_stage2_cost_best(B1, B2, D, A, L);
        int bad = 0;
        #pragma omp critical
        while (1)
        {
            for (i = 0; i < best.size() && ((best[i][0] > size && best[i][1] < cost) || (best[i][0] < size && best[i][1] > cost)); i++);
            if (i < best.size() && best[i][0] <= size && best[i][1] <= cost)
                bad = 1;
            break;
        }
        if (bad)
            continue;

        bool saved = false;
        for (auto it = savedCosts.begin(); it != savedCosts.end() && !saved; it++)
            if ((*it)[2] == D && (*it)[3] == A && (*it)[4] == L)
            {
                cost = (*it)[1];
                pairing = 2.0*(*it)[5]/num_primes;
                saved = true;
            }
        if (!saved)
            get_stage2_cost_and_pairing(primes, B1, B2, D, A, L, &cost, &pairing);

        #pragma omp critical
        while (1)
        {
            if (!saved)
            {
                fOut << B1 << " " << B2 << " " << size << " " << cost << " " << D << " " << A << " " << L << " " << (int)(num_primes*pairing/2) << std::endl;
                fOut.flush();
            }
            for (i = 0; i < best.size() && ((best[i][0] > size && best[i][1] < cost) || (best[i][0] < size && best[i][1] > cost)); i++);
            if (i < best.size() && best[i][0] <= size && best[i][1] <= cost)
                break;
            //printf("size=%d, D=%d, A=%d, L=%d, pairs=%d\n", size, D, A, L, (int)(num_primes*pairing/2));
            if (i < best.size())
                best[i] = { size, cost, D, A, L, (int)(pairing*1000) };
            else
                best.push_back({ size, cost, D, A, L, (int)(pairing*1000) });
            i++;
            for (; i < best.size(); i++)
                if (best[i][0] >= size && best[i][1] >= cost)
                {
                    best.erase(best.begin() + i);
                    i--;
                }
            break;
        }
    }
    fOut.close();

    /*std::sort(best.begin(), best.end(), [](std::array<int, 6>& a, std::array<int, 6>& b) { return a[0] < b[0]; });
    printf("size    cost    D  A  L %%pair\n");
    printf("---- ------- ---- -- -- -----\n");
    for (i = 0; i < best.size(); i++)
        printf("%4d %7d %4d %2d %2d %5d\n", best[i][0], best[i][1], best[i][2], best[i][3], best[i][4], best[i][5]);*/

    int prevk = -1;
    int bounds[] = { 50, 100, 150, 200, 250, 350, 450, 550, 750, 1000 };
    for (j = 0; j < 10; j++)
    {
        k = -1;
        int min = B2;
        for (i = 0; i < best.size(); i++)
            if (min > best[i][1] && best[i][0] < bounds[j])
            {
                min = best[i][1];
                k = i;
            }
        if (k != prevk)
        {
            printf("  %d,%d,%d,%d,%d,%d,%d, \n", B1/1000, B2/B1, best[k][0], best[k][2], best[k][3], best[k][4], best[k][5]);
            std::ofstream fRes("result.txt", std::fstream::app);
            fRes << "  " << B1/1000 << "," << B2/B1 << "," << best[k][0] << "," << best[k][2] << "," << best[k][3] << "," << best[k][4] << "," << best[k][5] << ", " << std::endl;
            fRes.flush();
        }
        prevk = k;
    }
}

void precompute_stage2_params()
{
    int b2mul[] = { 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 23, 26, 30, 35, 40, 50, 60, 75, 100 };

    int b1, b2;
    for (b1 = 1000; b1 < 10000; b1 += 1000)
        for (b2 = 0; b2 < 20; b2++)
            do_stage2_params(b1, b2mul[b2]*b1);
    for (b1 = 10000; b1 < 100000; b1 += 10000)
        for (b2 = 0; b2 < 20; b2++)
            do_stage2_params(b1, b2mul[b2]*b1);
    for (b1 = 100000; b1 <= 1000000/2; b1 += 100000)
        for (b2 = 0; b2 < 20; b2++)
            do_stage2_params(b1, b2mul[b2]*b1);
}

void precompute_DAC_S_d()
{
    PrimeList primes(65536);
    for (PrimeIterator it = primes.begin(); *it < 10000000; it++)
    {
        int Slen = 60;
        int Sd = get_DAC_S_d(*it, 1, *it/2 + 1, &Slen);
        printf("%d,\n", Sd);
    }
}
