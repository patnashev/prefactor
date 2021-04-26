
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
        stage1cost *= 1.5;
    if (stage1param > 1)
        stage1cost = (15 << (stage1param - 2)) + stage1cost*(7 + 7/(stage1param + 1.0));
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
    int i, j;
    int size = 0;

    size++;
    size++;
    size++;

    if (A > 1)
    {
        D /= A;
        L = L*A + 1;
    }
    size++;
    size++;
    char *precomp_V = malloc(D/2*L);
    memset(precomp_V, 0, D/2*L);
    size++;
    if (A > 1)
    {
        size++;
        size++;
    }

    size++;
    int dist = 1;
    for (i = 1; i < D/2*L; i++)
    {
        if (gcd(2*i + 1, 2*dist) != 1)
            continue;
        precomp_V[i] = 1;
        size++;
        if ((i + 1)%dist == 0 && D%(i + 1) == 0 && gcd((i + 1)/dist, 2*dist) == 1)
        {
            dist = i + 1;
            for (j = 0; j < i; j++)
                if (precomp_V[j] != 0 && gcd(2*j + 1, 2*dist) != 1)
                {
                    precomp_V[j] = 0;
                    size--;
                }
        }
    }
    if (A > 1)
    {
        D *= A;
        L /= A;
    }

    if (A > 1)
        size += L + 1;
    size += L + 1;

    free(precomp_V);
    return size;
}

int get_stage2_cost(int B1, int B2, int D, int A, int L, double pairing)
{
    int i, j;
    int stage2cost = 0;
    double num_primes = Ei(log(B2)) - Ei(log(B1));

    if (A > 1)
    {
        D /= A;
        L = L*A + 1;
    }

    if (primes == NULL)
        sieve(B2/B1 + 100);
    for (i = 0; primes[i]*B1 < B2; i++)
        if (D%primes[i] != 0)
        {
            B1 = B2/primes[i];
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
        L /= A;
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

    if (A > 1)
    {
        if (v > 0)
            stage2cost += 2*(int)(log(v*A)/log(2.0));
        stage2cost += A;
        stage2cost += L - 1;
    }
    else
    {
        if (v > 0)
            stage2cost += 2*(int)(log(v)/log(2.0));
    }

    stage2cost += L;
    stage2cost += (int)(num_primes - num_primes*pairing/2);
    stage2cost += (B2 - B1)/D;
    if (A > 1)
        stage2cost += (B2 - B1)/D;

    return stage2cost;
}

int get_stage2_cost_best(int B1, int B2, int D, int A, int L, int num_primes)
{
    int i, j;
    int stage2cost = 0;

    if (A > 1)
    {
        D /= A;
        L = L*A + 1;
    }

    for (i = 0; primes[i]*B1 < B2; i++)
        if (D%primes[i] != 0)
        {
            B1 = B2/primes[i];
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
        L /= A;
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

    if (A > 1)
    {
        if (v > 0)
            stage2cost += 2*(int)(log(v*A)/log(2.0));
        stage2cost += A;
        stage2cost += L - 1;
    }
    else
    {
        if (v > 0)
            stage2cost += 2*(int)(log(v)/log(2.0));
    }

    stage2cost += L;

    stage2cost += num_primes/2;
    stage2cost += (B2 - B1)/D;
    if (A > 1)
        stage2cost += (B2 - B1)/D;

    return stage2cost;
}

void get_stage2_cost_and_pairing(int* primes, int B1, int B2, int D, int A, int L, int *cost, double *pairing)
{
    int i, j;
    int stage2cost = 0;

    if (A > 1)
    {
        D /= A;
        L = L*A + 1;
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
        L /= A;
    }
    if (D/2/dist != 1)
    {
        int pd = D/2/dist;
        for (j = 0; (pd & 1) == 0; pd >>= 1, j++);
        if (pd > 1)
            stage2cost++;
        stage2cost += j;
    }

    int it = 0;
    while (primes[it] <= B1)
        it++;
    if (A > 1)
        D /= A;
    for (i = 0; primes[i]*B1 < B2; i++)
        if (D%primes[i] != 0)
        {
            for (j = it; primes[j]*primes[i] <= B2; j++)
                while (primes[j]*primes[i] <= B2)
                    primes[j] *= primes[i];
            for (; primes[j] <= B2; j++);
            qsort(primes + it, j - it, sizeof(int), cmp);
            break;
        }
    if (A > 1)
        D *= A;

    int paired = 0;
    short* distances;
    if (A > 1)
        distances = get_distances_A(primes, B1, B2, D, A, L);
    else
        distances = get_distances_simple(primes, B1, B2, D, L);
    for (j = 0; distances[j] != -1; j++)
        if (distances[j] == 0)
            paired++;
    free(distances);

    int v = primes[it]/D;

    if (A > 1)
    {
        if (v > 0)
            stage2cost += 2*(int)(log(v*A)/log(2.0));
        stage2cost += A;
        stage2cost += L - 1;
    }
    else
    {
        if (v > 0)
            stage2cost += 2*(int)(log(v)/log(2.0));
    }

    stage2cost += L;

    stage2cost += j - paired;
    stage2cost += (B2 - primes[it])/D;
    if (A > 1)
        stage2cost += (B2 - primes[it])/D;

    *cost = stage2cost;
    *pairing = 2*paired/(double)j;
}

void do_stage2_params(int B1, int B2)
{
    int i, j, k;
    int best[100][6];
    best[0][0] = 0;
    int totalp;
    for (totalp = 0; primes[totalp] <= B1; totalp++);
    int num_primes;
    for (num_primes = 0; primes[totalp] <= B2; totalp++, num_primes++);
    totalp++;

    int Ds[] = {66, 294, 588};

    int D;
    for (D = 30; D < 2150; D += 6)
    {
        int *plist = malloc(sizeof(int)*totalp);
        for (int L = 1; L < 20 && L < 4300/D && D*(L + 1) < 16000; L++)
        {
            for (int iA = 0; iA < 10; iA++)
            {
                int A = 1;
                if (iA > 0)
                    A = primes[iA];
                if (D%A != 0 || A == 2)
                    continue;

                memcpy(plist, primes, sizeof(int)*totalp);
                int size = get_stage2_size(D, A, L);
                if (size > 1000)
                    continue;
                double pairing = 1.0;
                int cost = get_stage2_cost_best(B1, B2, D, A, L, num_primes);
                int bad = 0;
                while (1)
                {
                    for (i = 0; i < 100 && best[i][0] != 0 && ((best[i][0] > size && best[i][1] < cost) || (best[i][0] < size && best[i][1] > cost)); i++);
                    if (i == 100)
                        printf("overflow\n");
                    if (best[i][0] != 0 && best[i][0] <= size && best[i][1] <= cost)
                        bad = 1;
                    break;
                }
                if (bad)
                    continue;

                get_stage2_cost_and_pairing(plist, B1, B2, D, A, L, &cost, &pairing);

                while (1)
                {
                    for (i = 0; i < 100 && best[i][0] != 0 && ((best[i][0] > size && best[i][1] < cost) || (best[i][0] < size && best[i][1] > cost)); i++);
                    if (i == 100)
                        printf("overflow\n");
                    if (best[i][0] != 0 && best[i][0] <= size && best[i][1] <= cost)
                        break;
                    if (best[i][0] == 0 && i < 100 - 1)
                        best[i + 1][0] = 0;
                    best[i][0] = size;
                    best[i][1] = cost;
                    best[i][2] = D;
                    best[i][3] = A;
                    best[i][4] = L;
                    best[i][5] = (int)(pairing*1000);
                    i++;
                    for (; i < 100 && best[i][0] != 0; i++)
                        if (best[i][0] >= size && best[i][1] >= cost)
                        {
                            for (j = i + 1; j < 100 && best[j][0] != 0; j++)
                                for (k = 0; k < 6; k++)
                                    best[j - 1][k] = best[j][k];
                            best[j - 1][0] = 0;
                            i--;
                        }
                    break;
                }
            }
        }
        free(plist);
    }
    int prevk = -1;
    int bounds[] = { 50, 100, 150, 200, 250, 350, 450, 550, 750, 1000 };
    for (j = 0; j < 10; j++)
    {
        k = -1;
        int min = B2;
        for (i = 0; i < 100 && best[i][0] != 0; i++)
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
    if (primes != NULL)
        free(primes);
    sieve(100*1000000 + 1000);

    int b2mul[] = { 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 23, 26, 30, 35, 40, 50, 60, 75, 100 };

    int b1, b2;
    /*for (b1 = 1000; b1 < 10000; b1 += 1000)
        for (b2 = 0; b2 < 20; b2++)
            do_stage2_params(b1, b2mul[b2]*b1);*/
    do_stage2_params(30000, 100*30000);
    for (b1 = 40000; b1 <= 100000; b1 += 10000)
        for (b2 = 0; b2 < 20; b2++)
            do_stage2_params(b1, b2mul[b2]*b1);
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
    char *precomp_V = malloc(D/2*L);
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

    free(precomp_V);
    return size;
}

int get_edecm_stage2_cost(int B1, int B2, int D, int L, int LR, double pairing)
{
    int i, j;
    int stage2cost = 0;
    double num_primes = Ei(log(B2)) - Ei(log(B1));

    if (primes == NULL)
        sieve(B2/B1 + 100);
    for (i = 0; primes[i]*B1 < B2; i++)
        if (D%primes[i] != 0)
        {
            B1 = B2/primes[i];
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

