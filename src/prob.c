
#include <cmath>

#define F_len 10000
double F_val[F_len + 1];
int F_inited = 0;
double rho_val[] = { 1, 1, 3.0685282e-1, 4.8608388e-2, 4.9109256e-3, 3.5472470e-4, 1.9649696e-5, 8.7456700e-7, 3.2320693e-8, 1.0162483e-9, 2.7701718e-11, 0.664480e-12, 0.141971e-13, 0.272918e-15, 0.476063e-17, 0.758990e-19 };

double F(double x)
{
    if (x >= 1.0)
        return 1.0;
    int i = (int)(x*F_len);
    double t = x - i/(double)F_len;
    return F_val[i]*(1 - t) + F_val[i + 1]*t;
}

void F_init()
{
    int i, j;
    F_val[F_len] = 1.0;
    for (i = F_len - 1; i >= F_len/4; i--)
        F_val[i] = F_val[i + 1] - F(i/(double)(F_len - i))/i;
    double a = 0;
    double b = log(F_val[F_len/4]);
    for (j = 5; j <= 15; j++)
    {
        a = b;
        b = log(rho_val[j]);
        for (; i >= (double)F_len/j; i--)
            F_val[i] = exp(a + (b - a)*((double)F_len/i - (j - 1)));
    }
    for (; i > 0; i--)
        F_val[i] = exp(a + (b - a)*((double)F_len/i - (j - 2)));
    F_val[0] = 0.0;
    F_inited = 1;
}

double prob_smooth(int B1, int B2, double sievingDepth, double knownDivisors)
{
    if (!F_inited)
        F_init();

    int i, j;
    double logB1 = log(B1)/log(2.0);
    double logB2 = log(B2)/log(2.0);
    double delta = 0.1;
    double sum = 0;
    for (i = 0; i < 1000; i++)
    {
        double l = sievingDepth + i*delta; // log(factor)
        double d = sievingDepth*(1/l - 1/(l + delta)); // probability of one factor in the range
        double ds = d; // probability of up to 5 factors in the range
        for (j = 0; j < 5; j++)
            ds = ds*d + d;
        double a = logB1/(l - 1 - knownDivisors); // B1 = (factor/2/knownDivisors)^a
        double b = logB2/(l - 1 - knownDivisors); // B2 = (factor/2/knownDivisors)^b
        double stage1 = F(a); // probability of B1-smooth factor
        double stage2 = 0; // integrating stage 2 probabilities
        for (j = (int)(a*F_len); j < (int)(b*F_len); j++)
            stage2 += F(a*F_len/(double)(F_len - j))/j;
        sum += (stage1 + stage2)*ds; // probability of successful factorization
    }

    return sum;
}

int get_stage1_cost(int B1, int stage1param)
{
    int stage1cost = (int)(B1/0.69);
    if (stage1param == 1)
        stage1cost = (int)(stage1cost*1.5);
    if (stage1param > 1)
        stage1cost = (15 << (stage1param - 2)) + (int)(stage1cost*(7 + 7/(stage1param + 1.0)));
    return stage1cost;
}

void get_edecm_stage1_params(int B1, int maxSize, int *K)
{
    int k;
    for (k = 2; k < 16 && 3 + 3*(1 << (k - 1)) <= maxSize && (15 << (k - 2)) + B1/0.69*(7 + 7/(k + 1.0)) > (15 << (k - 1)) + B1/0.69*(7 + 7/(k + 2.0)); k++);
    *K = k;
}

int get_edecm_stage1_size(int K)
{
    return 3 + 3*(1 << (K - 2));
}

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

int get_stage2_cost(int B1, int B2, int D, int A, int L, double pairing)
{
    int i, j;
    int stage2cost = 0;

    double num_primes = std::expint(log(B2)) - std::expint(log(B1));
    stage2cost += (int)(num_primes - pairing*num_primes/2);

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
            stage2cost += 2*(int)(log(v*A)/log(2.0));
        stage2cost += A;
    }
    else
    {
        if (v > 0)
            stage2cost += 2*(int)(log(v)/log(2.0));
    }

    stage2cost += (B2 - B1)/D;
    if (A > 1)
        stage2cost += (B2 - B1)/D;

    return stage2cost;
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
            stage2cost += 2*(int)(log(v*A)/log(2.0));
        stage2cost += A;
    }
    else
    {
        if (v > 0)
            stage2cost += 2*(int)(log(v)/log(2.0));
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
    B1 = PP.first_D*D;
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

    int v = B1/D;

    if (A > 1)
    {
        if (v > 0)
            stage2cost += 2*(int)(log(v*A)/log(2.0));
        stage2cost += A;
    }
    else
    {
        if (v > 0)
            stage2cost += 2*(int)(log(v)/log(2.0));
    }

    stage2cost += (B2 - B1)/D;
    if (A > 1)
        stage2cost += (B2 - B1)/D;

    *cost = stage2cost;
}

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

    std::vector<std::array<int,6>> savedCosts;
    std::ifstream fIn("costs.txt");
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

    std::ofstream fOut("costs.txt", std::fstream::app);

    int D;
    for (D = 30; D < 2150; D += 6)
    {
        for (int L = 1; L < 15 && L < 4300/D && D*(L + 1) < 16000; L++)
        {
            for (int iA = 0; iA < 10; iA++)
            {
                int A = 1;
                if (iA > 0)
                    A = primes[iA];
                if (D%A != 0 || A == 2)
                    continue;

                int size = get_stage2_size(D, A, L);
                if (size > 1000)
                    continue;
                double pairing = 1.0;
                int cost = get_stage2_cost_best(B1, B2, D, A, L);
                int bad = 0;
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
            printf("  %d,%d,%d,%d,%d,%d,%d, \n", B1/1000, B2/B1, best[k][0], best[k][2], best[k][3], best[k][4], best[k][5]);
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
    /*for (b1 = 10000; b1 < 100000; b1 += 10000)
        for (b2 = 0; b2 < 20; b2++)
            do_stage2_params(b1, b2mul[b2]*b1);*/
    /*for (b1 = 1000000; b1 <= 1000000; b1 += 100000)
        for (b2 = 4; b2 < 20; b2 += 5)
            do_stage2_params(b1, b2mul[b2]*b1);*/
}

void get_stage2_params(int B1, int B2, int maxSize, int *D, int *A, int *L, double *pairing)
{
    int i;
    int b2mul[] = { 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 23, 26, 30, 35, 40, 50, 60, 75, 100 };
    B2 /= B1;
    for (i = 0; i < 19 && B2 >= b2mul[i + 1]; i++);
    B2 = b2mul[i];
    B1 /= 1000;
    if (B1 < 1)
        B1 = 1;
    if (B1 > 10)
        B1 -= B1%10;
    if (B1 > 100)
        B1 -= B1%100;
    if (B2 > 100)
        B2 = 100;

    int l = 0;
    int r = precomputed_stage2_params_len/7 - 1;
    while (l < r - 1)
    {
        int m = (l + r)/2;
        if (precomputed_stage2_params[m*7 + 0] > B1)
            r = m;
        else if (precomputed_stage2_params[m*7 + 0] < B1)
            l = m;
        else if (precomputed_stage2_params[m*7 + 1] > B2)
            r = m;
        else if (precomputed_stage2_params[m*7 + 1] < B2)
            l = m;
        else if (precomputed_stage2_params[m*7 + 2] > maxSize)
            r = m;
        else
            l = m;
    }
    *D = precomputed_stage2_params[l*7 + 3];
    *A = precomputed_stage2_params[l*7 + 4];
    *L = precomputed_stage2_params[l*7 + 5];
    if (pairing != NULL)
        *pairing = precomputed_stage2_params[l*7 + 6]/1000.0;
}

double get_profit(int B1, int B2, int maxSize, double sievingDepth, double knownDivisors, double primalityCost, int stage1param)
{
    int D, A, L;
    double pairing;
    get_stage2_params(B1, B2, maxSize, &D, &A, &L, &pairing);
    int stage1cost = get_stage1_cost(B1, stage1param);
    int stage2cost = get_stage2_cost(B1, B2, D, A, L, pairing);
    double value = prob_smooth(B1, B2, sievingDepth, knownDivisors);
    return value - (stage1cost + stage2cost)/primalityCost;
}

void get_optimal_bounds(int *B1, int *B2, int maxSize, double sievingDepth, double knownDivisors, double primalityCost, int stage1param)
{
    int i, j;
    int b1, b2;
    int b1min = 0;
    int b1max = 1024000;
    int b2min = 0;
    int b2max = 10240000;
    double p0, p1;

    for (i = 0; i < 10; i++)
    {
        b1 = (b1min + b1max)/2;
        b2min = 0;
        b2max = 10240000;
        for (j = 0; j < 10; j++)
        {
            b2 = (b2min + b2max)/2;
            p0 = get_profit(b1, b2, maxSize, sievingDepth, knownDivisors, primalityCost, stage1param);
            p1 = get_profit(b1, b2 + 10000, maxSize, sievingDepth, knownDivisors, primalityCost, stage1param);
            if (p1 > p0)
                b2min = b2;
            else
                b2max = b2;
        }
        b2 = b2max;
        p0 = get_profit(b1, b2, maxSize, sievingDepth, knownDivisors, primalityCost, stage1param);
        p1 = get_profit(b1 + 1000, b2, maxSize, sievingDepth, knownDivisors, primalityCost, stage1param);
        if (p1 > p0)
            b1min = b1;
        else
            b1max = b1;
    }
    b1 = b1max;
    *B1 = b1;
    *B2 = b2;
}

void get_optimal_B2(int B1, int *B2, int maxSize, double sievingDepth, double knownDivisors, double primalityCost, int stage1param)
{
    int j;
    int b2;
    int b2min = 0;
    int b2max = 10240000;
    double p0, p1;

    b2min = 0;
    b2max = 10240000;
    for (j = 0; j < 10; j++)
    {
        b2 = (b2min + b2max)/2;
        p0 = get_profit(B1, b2, maxSize, sievingDepth, knownDivisors, primalityCost, stage1param);
        p1 = get_profit(B1, b2 + 10000, maxSize, sievingDepth, knownDivisors, primalityCost, stage1param);
        if (p1 > p0)
            b2min = b2;
        else
            b2max = b2;
    }
    b2 = b2max;
    *B2 = b2;
}

int get_edecm_stage2_size(int D, int L, int LR)
{
    int i, j;
    int size = 0;

    size++;
    size++;
    size += 4;
    size += 4;
    size += 2;
    size += 2;
    char *precomp_V = new char[D/2*L];
    memset(precomp_V, 0, D/2*L);

    size += 1;
    int dist = 1;
    for (i = 1; i < D/2*L; i++)
    {
        if (gcd(2*i + 1, 2*dist) != 1)
            continue;
        precomp_V[i] = 1;
        size += 1;
        if ((i + 1)%dist == 0 && D%(i + 1) == 0 && gcd((i + 1)/dist, 2*dist) == 1)
        {
            dist = i + 1;
            for (j = 0; j < i; j++)
                if (precomp_V[j] != 0 && gcd(2*j + 1, 2*dist) != 1)
                {
                    precomp_V[j] = 0;
                    size -= 1;
                }
        }
    }

    if (LR > 0)
        size += L + LR;
    else
        size += 2*L;

    delete precomp_V;
    return size;
}

int get_edecm_stage2_cost(int B1, int B2, int D, int L, int LR, double pairing)
{
    int i, j;
    int stage2cost = 0;
    double num_primes = std::expint(log(B2)) - std::expint(log(B1));

    PrimeList primes(B2/B1 + 100);
    for (PrimeIterator it = primes.begin(); *it*B1 < B2; it++)
        if (D%(*it) != 0)
        {
            B1 = B2/(*it);
            break;
        }

    stage2cost++;
    stage2cost++;
    int dist = 1;
    //stage2cost += D/2*L/2;
    for (i = 1; i < D/2*L; i++)
    {
        if (gcd(2*i + 1, 2*dist) != 1)
            continue;
        stage2cost++;
        if ((i + 1)%dist == 0 && D%(i + 1) == 0 && gcd((i + 1)/dist, 2*dist) == 1)
        {
            stage2cost++;
            stage2cost++;
            dist = i + 1;
        }
    }

    if (D/2/dist != 1)
    {
        int pd = D/2/dist;
        for (j = 0; (pd & 1) == 0; pd >>= 1, j++);
        if (pd > 1)
            stage2cost++;
        stage2cost += j;
    }

    int v = B1/D;

    if (v > 0)
        stage2cost += 2*(int)(log(v)/log(2.0));

    stage2cost += L;
    stage2cost += (int)(num_primes - num_primes*pairing/2);
    stage2cost += (B2 - B1)/D;

    return stage2cost;
}

