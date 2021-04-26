
// Stage 1 exponent with prime powers <= B1
giant get_stage1_exp(int B0, int B1)
{
    int j, k;
    int p, it;

    for (it = 0; primes[it] < B0; it++);
    p = primes[it++];

    int sqrtB1 = (int)sqrt(B1);
    giant tmp = popg(&gwdata->gdata, ((int)(gwdata->bit_length > B1/0.69 ? gwdata->bit_length : B1/0.69) >> 5) + 10);
    tmp->sign = 0;
    giant tmp2 = popg(&gwdata->gdata, ((int)(262144 < B1/0.69 ? 262144 : B1/0.69) >> 5) + 10);
    ultog(1, tmp2);
    while (p <= B1)
    {
        // Building exponent with prime powers <= B1
        j = p;
        if (p <= sqrtB1)
        {
            k = B1/p;
            while (j <= k)
                j *= p;
        }
        ulmulg(j, tmp2);
        p = primes[it++];
        if (p > B1 || tmp2->sign > 8192)
        {
            if (tmp->sign == 0)
                gtog(tmp2, tmp);
            else
                mulg(tmp2, tmp);
            ultog(1, tmp2);
        }
    }
    freeg();

    return tmp;
}

// P-1 factoring stage 1.
int do_minus1stage1(int B1, int gfn, giant res)
{
    int i, j;
    int len;
    int retval = TRUE;

    gwnum X = gwalloc(gwdata);

    printf("%s, P-1 stage 1, B1 = %d%s\n", Nstr, B1, gfn != 0 ? ", GFN" : "");
    double timer = getHighResTimer();
    int transforms = -(int)gwdata->fft_count;

    giant tmp = get_stage1_exp(3, B1);
    gwsetmulbyconst(gwdata, 3);
    dbltogw(gwdata, 3, X);
    len = bitlen(tmp) - 1;
    gwset_carefully_count(gwdata, 30);
    for (i = 1; i <= len; i++, costAdd(1))
        gwsquare2(gwdata, X, X, (bitval(tmp, len - i) ? GWMUL_MULBYCONST : 0) | (i < len ? GWMUL_STARTNEXTFFT : 0));
    for (j = 2; j <= B1; j <<= 1, costAdd(1))
        gwsquare_carefully(gwdata, X);
    for (j = 2; j <= gfn; j <<= 1, costAdd(1))
        gwsquare_carefully(gwdata, X);

    gwtogiant(gwdata, X, tmp);
    gtog(tmp, res);
    invg(N, res);
    if (res->sign > 0)
        addg(tmp, res);
    ulsubg(1, tmp);
    gcdg(N, tmp);

    timer = (getHighResTimer() - timer)/getHighResTimerFrequency();
    transforms += (int)gwdata->fft_count;
    if (isZero(tmp) || gcompg(tmp, N) == 0)
    {
        printf("All divisors of N < B1.\n");
        retval = FALSE;
    }
    else if (!isone(tmp))
    {
        report_factor(tmp);
        retval = FALSE;
    }
    else
    {
        printf("No factors found, transforms: %d, time: %d s.\n", transforms, (int)timer);
    }

    gwfree(gwdata, X);
    freeg();
    return retval;
}

// P+1 factoring stage 1.
int do_plus1stage1(int B0, int B1, int gfn, giant P, char *sP, giant res)
{
    int j, k;
    int p, it;
    int retval = TRUE;
    giant tmp = NULL;

    gwnum V = gwalloc(gwdata);

    printf("%s, P+1 stage 1, B1 = %d, P = %s%s.\n", Nstr, B1, sP, gfn != 0 ? ", GFN" : "");
    double timer = getHighResTimer();
    int transforms = -(int)gwdata->fft_count;


    // V_tmp with careful start
    gianttogw(gwdata, P, V);
    if (B0 <= 2)
    {
        int power2 = 0;
        for (j = 2; j <= B1; j <<= 1)
            power2++;
        for (j = 2; j <= gfn; j <<= 1)
            power2++;
        gwset_carefully_count(gwdata, power2);
        lucas_V_shiftleft(power2, V);
    }
    if (B0 < 3)
        B0 = 3;

    for (it = 0; primes[it] < B0; it++);
    p = primes[it];

    int sqrtB1 = (int)sqrt(B1);
    while (p <= B1)
    {
        j = p;
        lucas_V_mul_prime(it, &V);
        if (p <= sqrtB1)
        {
            k = B1/p;
            while (j <= k)
            {
                lucas_V_mul_prime(it, &V);
                j *= p;
            }
        }
        it++;
        p = primes[it];
    }

    tmp = getg();
    gwtogiant(gwdata, V, tmp);
    gtog(tmp, res);
    ulsubg(2, tmp);
    gcdg(N, tmp);

    timer = (getHighResTimer() - timer)/getHighResTimerFrequency();
    transforms += (int)gwdata->fft_count;
    if (isZero(tmp) || gcompg(tmp, N) == 0)
    {
        printf("All divisors of N < B1.\n");
        retval = FALSE;
    }
    else if (!isone(tmp))
    {
        report_factor(tmp);
        retval = FALSE;
    }
    else
    {
        printf("No factors found, transforms: %d, time: %d s.\n", transforms, (int)timer);
    }

    gwfree(gwdata, V);
    freeg();
    return retval;
}

// EdECM factoring stage 1, mul&add small point.
int do_edecm_stage1_simple(int B1)
{
    int i;
    int len;
    int retval = TRUE;

    ///int xa = 5;
    //int xb = 23;
    //int ya = -1;
    //int yb = 7;
    int xa = 17;
    int xb = 19;
    int ya = 17;
    int yb = 33;

    ed_point P = ed_alloc();
    ed_from_small(xa, xb, ya, yb, P);

    printf("%s, EdECM stage 1, B1 = %d.\n", Nstr, B1);
    double timer = getHighResTimer();

    gwset_carefully_count(gwdata, 50);
    giant tmp = get_stage1_exp(2, B1);
    len = bitlen(tmp) - 1;
    for (i = 1; i <= len; i++)
        if (bitval(tmp, len - i))
            ed_mul2_addsmall(xa, xb, ya, yb, P, i < len ? GWMUL_STARTNEXTFFT : 0);
        else
            ed_mul2(P, i < len ? GWMUL_STARTNEXTFFT : 0);

    gwtogiant(gwdata, P->X, tmp);
    gcdg(N, tmp);

    timer = (getHighResTimer() - timer)/getHighResTimerFrequency();
    if (isZero(tmp) || gcompg(tmp, N) == 0)
    {
        printf("All divisors of N < B1.\n");
        retval = FALSE;
    }
    else if (!isone(tmp))
    {
        report_factor(tmp);
        retval = FALSE;
    }
    else
    {
        printf("No factors found, time: %d s.\n", (int)timer);
    }

    ed_free(P);
    freeg();
    return retval;
}

// EdECM factoring stage 1 using signed window.
int do_edecm_stage1(int B1, int W, ed_point P)
{
    int j;
    int len;
    int retval = TRUE;
    giant tmp = NULL;
    short *nafw = NULL;

    if (W > 16)
        W = 16;
    tmp = get_stage1_exp(3, B1);
    get_nafw(tmp, W, &nafw, &len);
    freeg();
    tmp = NULL;

    printf("%s, EdECM stage 1, B1 = %d, W = %d.\n", Nstr, B1, W);
    double timer = getHighResTimer();
    int transforms = -(int)gwdata->fft_count;

    for (j = 2; j <= B1; j <<= 1)
        ed_mul2_carefully(P);
    ed_mul2(P, ED_MUL_FULL);

    tmp = ed_mul_nafw(W, nafw, len, P);
    free(nafw);

    if (tmp == NULL)
    {
        tmp = getg();
        gwtogiant(gwdata, P->X, tmp);
        gcdg(N, tmp);
    }
    timer = (getHighResTimer() - timer)/getHighResTimerFrequency();
    transforms += (int)gwdata->fft_count;
    if (isZero(tmp) || gcompg(tmp, N) == 0)
    {
        printf("All divisors of N < B1.\n");
        retval = FALSE;
    }
    else if (!isone(tmp))
    {
        report_factor(tmp);
        retval = FALSE;
    }
    else
    {
        printf("No factors found, transforms: %d, time: %.3f s.\n", transforms, timer);
    }

    freeg();
    return retval;
}
