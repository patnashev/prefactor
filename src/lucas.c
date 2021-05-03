
#define GW_FMA

void lucas_V_shiftleft(int x, gwnum V)
{
    gwsetaddin(gwdata, -2);
    for (; x > 0; x--, costAdd(1))
        gwsquare2(gwdata, V, V, GWMUL_ADDINCONST | GWMUL_STARTNEXTFFT);
    gwsetaddin(gwdata, 0);
}

void lucas_V_inc(gwnum Vinc, gwnum Vprev, gwnum Vcur, gwnum Vnext)
{
#ifdef GW_FMA
    gwmulsub4(gwdata, Vinc, Vcur, Vprev, Vnext, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_FFT_S3 | GWMUL_STARTNEXTFFT);
#else
    gwmul3(gwdata, Vinc, Vcur, Vnext, GWMUL_FFT_S1 | options);
    gwsub3o(gwdata, Vnext, Vprev, Vnext, GWADD_FORCE_NORMALIZE);
#endif
    costAdd(1);
}

#define lucas_V_mul_2(V) lucas_V_shiftleft(1, V)
#define lucas_V_add(Va, Vb, Vbma, Vbpa) lucas_V_inc(Va, Vbma, Vb, Vbpa)

void lucas_V_mul_giant(giant g, gwnum V, gwnum Vres, gwnum Vres1)
{
    // Vres = V_1
    gwcopy(gwdata, V, Vres);
    // Vres1 = V_2
    gwcopy(gwdata, V, Vres1);
    gwsetaddin(gwdata, -2);
    gwsquare2(gwdata, Vres1, Vres1, GWMUL_ADDINCONST);
    gwsetaddin(gwdata, 0);
    costAdd(1);
#ifdef GW_FMA
    gwnum Vtmp = gwalloc(gwdata);
    gwfft_for_fma(gwdata, V, Vtmp);
#endif

    // Each iteration computes V_n = V_{2n+bit}, V_{n+1} = V_{2n+bit+1}
    int len = bitlen(g) - 1;
    gwsetaddin(gwdata, -2);
    for (int i = 1; i <= len; i++, costAdd(2))
    {
        if (bitval(g, len - i))
        {
#ifdef GW_FMA
            gwmulsub4(gwdata, Vres1, Vres, Vtmp, Vres, i < len ? GWMUL_STARTNEXTFFT : 0);
#else
            gwmul3(gwdata, Vres1, Vres, Vres, 0);
            gwsub3o(gwdata, Vres, V, Vres, i == len ? GWADD_FORCE_NORMALIZE : !bitval(g, len - i - 1) ? GWADD_SQUARE_INPUT : GWADD_MUL_INPUT);
#endif
            gwsquare2(gwdata, Vres1, Vres1, GWMUL_ADDINCONST + (i < len ? GWMUL_STARTNEXTFFT : 0));
        }
        else
        {
#ifdef GW_FMA
            gwmulsub4(gwdata, Vres, Vres1, Vtmp, Vres1, i < len ? GWMUL_STARTNEXTFFT : 0);
#else
            gwmul3(gwdata, Vres, Vres1, Vres1, 0);
            gwsub3o(gwdata, Vres1, V, Vres1, i == len || (must_norm && bitval(g, len - i - 1)) ? GWADD_SQUARE_INPUT : GWADD_MUL_INPUT);
#endif
            gwsquare2(gwdata, Vres, Vres, GWMUL_ADDINCONST + (i < len ? GWMUL_STARTNEXTFFT : 0));
        }
    }
    gwsetaddin(gwdata, 0);

#ifdef GW_FMA
    gwfree(gwdata, Vtmp);
#endif
}

void lucas_V_mul_int(int x, gwnum V, gwnum Vres, gwnum Vres1)
{
    if (x == 0)
    {
        dbltogw(gwdata, 2, Vres);
        if (Vres1 != NULL)
            gwcopy(gwdata, V, Vres1);
        return;
    }
    int len;
    for (len = 1; x >= (1 << len); len++);
    int i = x;
    if (len > 2)
        i >>= len - 2;
    len -= 2;
    gwnum Vn = Vres;
    gwnum Vn1 = Vres1;
    if (len > 0 && Vn1 == NULL)
        Vn1 = gwalloc(gwdata);
    gwsetaddin(gwdata, -2);
    if (i == 1)
    {
        gwcopy(gwdata, V, Vn); // V_1
        if (Vn1 != NULL)
        {
            gwsquare2(gwdata, V, Vn1, GWMUL_ADDINCONST | GWMUL_STARTNEXTFFT); // V_2
            costAdd(1);
        }
    }
    if (i == 2)
    {
        gwsquare2(gwdata, V, Vn, GWMUL_ADDINCONST | GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT); // V_2
        costAdd(1);
        if (Vn1 != NULL)
        {
            gwmulsub4(gwdata, V, Vn, V, Vn1, GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT); // V_3
            costAdd(1);
        }
    }
    if (i == 3)
    {
        if (Vn1 != NULL)
        {
            gwsquare2(gwdata, V, Vn1, GWMUL_ADDINCONST | GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT); // V_2
            gwmulsub4(gwdata, V, Vn1, V, Vn, GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT); // V_3
            gwsquare2(gwdata, Vn1, Vn1, GWMUL_ADDINCONST | GWMUL_STARTNEXTFFT); // V_4
            costAdd(3);
        }
        else
        {
            gwsquare2(gwdata, V, Vn, GWMUL_ADDINCONST | GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT); // V_2
            gwmulsub4(gwdata, V, Vn, V, Vn, GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT); // V_3
            costAdd(2);
        }
    }
    gwsetaddin(gwdata, 0);
    if (len <= 0)
        return;
    gwnum Vtmp = V;
    if (len > 1)
    {
        Vtmp = gwalloc(gwdata);
        gwfft_for_fma(gwdata, V, Vtmp);
    }

    // Each iteration computes V_n = V_{2n+bit}, V_{n+1} = V_{2n+bit+1}
    gwsetaddin(gwdata, -2);
    for (i = (1 << (len - 1)); i; i >>= 1, costAdd(2))
    {
        if (x & i)
        {
            gwmulsub4(gwdata, Vn1, Vn, Vtmp, Vn, GWMUL_STARTNEXTFFT);
            if (i > 1 || Vres1 != NULL)
                gwsquare2(gwdata, Vn1, Vn1, GWMUL_ADDINCONST | GWMUL_STARTNEXTFFT);
        }
        else
        {
            if (i > 1 || Vres1 != NULL)
                gwmulsub4(gwdata, Vn, Vn1, Vtmp, Vn1, GWMUL_STARTNEXTFFT);
            gwsquare2(gwdata, Vn, Vn, GWMUL_ADDINCONST | GWMUL_STARTNEXTFFT);
        }
    }
    gwsetaddin(gwdata, 0);

    if (Vres1 == NULL)
        gwfree(gwdata, Vn1);
    if (Vtmp != V)
        gwfree(gwdata, Vtmp);
}

int get_DAC_S_d(int e, int start, int end, int *maxlen)
{
    int minlen = *maxlen;
    int mind = -1;
    if (start < 1)
        start = 1;
    if (end > e)
        end = e;
    for (int d = start; d < end; d++)
    {
        if (gcd(d, e) != 1)
            continue;
        int te = e;
        int td = d;
        int len = 0;
        while (td != 0 && len < minlen)
        {
            if (te/2 < td)
                td = te - td;
            else if ((te & 1) == 0 && td < te/4)
            {
                len += 2;
                te = te/2;
            }
            else
            {
                len++;
                te = te - td;
            }
        }
        if (minlen > len)
        {
            if (td != 0 || te != 1)
                printf("error");
            minlen = len;
            mind = d;
        }
    }
    *maxlen = minlen;
    return mind;
}

int get_DAC_M_d(int e, int start, int end, int *maxlen)
{
    int minlen = *maxlen;
    int mind = -1;
    if (start < e/2 + 1)
        start = e/2 + 1;
    if (end > e)
        end = e;
    for (int d = start; d < end; d++)
    {
        if (gcd(d, e) != 1)
            continue;
        int te = e - d;
        int td = d - te;
        int len = 2;
        while (td != te && len < minlen)
        {
            if (td < te)
            {
                int t = td;
                td = te;
                te = t;
            }
            if (100*td <= 296*te || ((td & 1) == 1 && (te & 1) == 0))
            {
                len++;
                td = td - te;
            }
            else if ((td & 1) == (te & 1))
            {
                len += 2;
                td = (td - te) >> 1;
            }
            else
            {
                len += 2;
                td = td >> 1;
            }
        }
        if (minlen > len)
        {
            minlen = len;
            mind = d;
        }
    }
    *maxlen = minlen;
    return mind;
}

void precompute_DAC_S_d()
{
    for (int i = 1000; i < primeCount; i++)
    {
        int Slen = 60;
        int Sd = get_DAC_S_d(primes[i], 1, primes[i]/2 + 1, &Slen);
        printf("%d,\n", Sd);

        /*int Mlen = 60;
        int Md = get_DAC_M_d(primes[i], primes[i]/2 + 1, primes[i], &Mlen);
        if (Slen < Mlen)
            printf("S %d\n", primes[i]);
        if (Mlen < Slen)
            printf("M %d\n", primes[i]);*/
    }
}

void lucas_V_mul_prime(int index, gwnum *V)
{
    if (index < 5)
    {
        gwnum Vres = gwalloc(gwdata);
        lucas_V_mul_int(primes[index], *V, Vres, NULL);
        gwswap(*V, Vres);
        gwfree(gwdata, Vres);
        return;
    }
    int i;
    int chain[100];
    int len = 60;
    int e = primes[index];
    int d = 1;
    if (index < precomputed_DAC_S_d_len)
        d = precomputed_DAC_S_d[index];
    else
        d = get_DAC_S_d(e, e*2/(1 + sqrt(5)) - 100, e*2/(1 + sqrt(5)) + 100, &len);
    len = 0;
    while (d != 0)
    {
        if (e/2 < d)
        {
            d = e - d;
            chain[len] = 0;
            len++;
        }
        else if ((e & 1) == 0 && d < e/4)
        {
            chain[len] = 2;
            len++;
            e = e/2;
        }
        else
        {
            chain[len] = 1;
            len++;
            e = e - d;
        }
    }
    if (len > 100)
        printf("error");
    d = 1;
    e = 2;
    int ed = 1;
    gwnum Ve = gwalloc(gwdata);
    gwnum Vd = *V;
    gwnum Ved = gwalloc(gwdata);
    gwsetaddin(gwdata, -2);
    gwmul3(gwdata, Vd, Vd, Ve, GWMUL_FFT_S1 | GWMUL_ADDINCONST | GWMUL_STARTNEXTFFT);
    gwsetaddin(gwdata, 0);
    costAdd(1);
    gwcopy(gwdata, Vd, Ved);
    for (i = len - 3; i >= 0; i--)
    {
        if (ed != e - d)
            printf("error");
        if (chain[i] == 2)
        {
            ed = ed + e;
            e = 2*e;
            lucas_V_inc(Ved, Vd, Ve, Ved, GWMUL_STARTNEXTFFT);
            lucas_V_mul_2(Ve);
        }
        else if (chain[i] == 1)
        {
            ed = e;
            e = e + d;
            lucas_V_inc(Vd, Ved, Ve, Ved, GWMUL_STARTNEXTFFT);
            gwswap(Ve, Ved);
        }
        else
        {
            ed = d;
            d = e - d;
            gwswap(Vd, Ved);
        }
    }
    if (e != primes[index])
        printf("error");
    *V = Ve;
    gwfree(gwdata, Vd);
    gwfree(gwdata, Ved);
}
