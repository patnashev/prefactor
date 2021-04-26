
int cmp(const void *a, const void *b) { return *(int*)a - *(int*)b; }
int cmpu64(const void *a, const void *b) { return *(uint64_t*)a > *(uint64_t*)b ? 1 : *(uint64_t*)a == *(uint64_t*)b ? 0 : -1; }

int gcd(int a, int b)
{
    while (a)
    {
        int t = a;
        a = b;
        b = t;
        a = a%b;
    }
    return b;
}

uint32_t gmodul(giant num, uint32_t den)
{
    uint64_t res = 0;
    for (int i = abs(num->sign) - 1; i >= 0; i--)
    {
        res <<= 32;
        res += num->n[i];
        res %= den;
    }
    return (uint32_t)res;
}

void get_nafw(giant exp, int W, short **res_nafw, int *res_len)
{
    int i, len;
    len = bitlen(exp);
    uint32_t *giant_n = exp->n;
    giant cur = allocgiant(2);
    short *nafw = malloc(sizeof(short)*(len + 1));
    int nafwMask = (1 << W) - 1;
    int nafwBit = (1 << (W - 1));
    int bit = 0;
    for (i = 0; !isZero(exp); i++)
    {
        if (bitval(exp, bit))
        {
            cur->n[0] = exp->n[0];
            cur->sign = 1;
            if (exp->sign > 1 && exp->n[1] != 0 && bit + W >= 32)
            {
                cur->n[1] = exp->n[1];
                cur->sign = 2;
            }
            gshiftright(bit, cur);
            nafw[i] = cur->n[0] & nafwMask;
            if (nafw[i] & nafwBit)
                nafw[i] -= nafwMask + 1;
            itog(-nafw[i], cur);
            gshiftleft(bit, cur);
            addg(cur, exp);
            bit += W;
        }
        else
        {
            nafw[i] = 0;
            bit++;
        }
        if (bit >= 32 && exp->sign > 0)
        {
            bit -= 32;
            exp->n++;
            exp->sign--;
        }
    }
    exp->n = giant_n;
    *res_nafw = nafw;
    *res_len = i;
}

int gwtogiantVerbose(gwhandle* gwdata, gwnum gg, giant v)
{
    int ret = gwtogiant(gwdata, gg, v);
    if (ret < 0)
        printf("Invalid FFT data.\n");
    return ret;
}

// Sieves primes < B
void sieve(int B)
{
    int i, j, k;
    char *bitmap = malloc(B/2);
    memset(bitmap, 0, B/2);
    k = 0;
    for (i = 1; i < B/2; i++)
        if (!bitmap[i])
        {
            k++;
            if (i > 16384)
                continue;
            for (j = (i*2 + 1)*(i*2 + 1)/2; j < B/2; j += i*2 + 1)
                bitmap[j] = 1;
        }
    if (primes != NULL)
        free(primes);
    primeCount = k + 1;
    primes = malloc(sizeof(int)*(k + 1));
    primes[0] = 2;
    k = 1;
    for (i = 1; i < B/2; i++)
        if (!bitmap[i])
        {
            primes[k] = i*2 + 1;
            k++;
        }
    free(bitmap);
}

double Ei(double x)
{
    double res = 0.57721566490153286 + log(x);
    double t = x;
    for (int i = 1; t > i; i++, t = t*x/i)
        res += t/i;
    return res;
}

void sieveRange(int start, int end, int **list, int *count)
{
    int i, j;
    if (!(start & 1))
        start++;
    int bitmapLen = (end - start)/2;
    char *bitmap = malloc(bitmapLen);
    memset(bitmap, 0, bitmapLen);
    for (i = 1; i < primeCount && primes[i]*primes[i] < end; i++)
    {
        j = (primes[i] - start%primes[i])%primes[i];
        for (j = ((j & 1) != 0 ? (j + primes[i]) : j)/2; j < bitmapLen; j += primes[i])
            bitmap[j] = 1;
    }
    *count = (Ei(log(end)) - Ei(log(start)))*1.1 + 100;
    *list = malloc(sizeof(int)*(*count));
    i = 0;
    if (start < sqrt(end))
    {
        for (j = 0; primes[j] < start; j++);
        for (; primes[j]*primes[j] < end; j++)
        {
            (*list)[i] = primes[j];
            i++;
        }
    }
    for (j = 0; j < bitmapLen; j++)
        if (!bitmap[j])
        {
            if (i >= *count)
            {
                printf("error\n");
                break;
            }
            (*list)[i] = start + j*2;
            i++;
        }
    *count = i;
    free(bitmap);
}

int costTotal = 0;
int costCurrent = 0;
int costLast = 0;
double costTimer = 0;

void costInit(int cost)
{
    costTotal = cost*2;
    costTimer = getHighResTimer();
}

void costAdd(int cost)
{
    costCurrent += cost*2;
    if (costCurrent - costLast >= 10000)
    {
        costTimer = (getHighResTimer() - costTimer)/getHighResTimerFrequency();
        printf("Progress : %d / %d [%.1f%%]. Speed : %.3f ms/op.\n", costCurrent, costTotal, costCurrent*100.0/costTotal, costTimer*1000.0/(costCurrent - costLast));
        costLast = costCurrent;
        costTimer = getHighResTimer();
    }
}
