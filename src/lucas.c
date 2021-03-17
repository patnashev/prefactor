

void lucas_V_mul_giant(giant g, gwnum V, gwnum Vres, gwnum Vres1)
{
    // Vres = V_1
    gwcopy(gwdata, V, Vres);
    // Vres1 = V_2
    gwcopy(gwdata, V, Vres1);
    gwsetaddin(gwdata, -2);
    gwsquare(gwdata, Vres1);
    gwsetaddin(gwdata, 0);
    costAdd(1);

    // Each iteration computes V_n = V_{2n+bit}, V_{n+1} = V_{2n+bit+1}
    int len = bitlen(g) - 1;
    for (int i = 1; i <= len; i++, costAdd(2))
    {
        if (bitval(g, len - i))
        {
            gwmul3(gwdata, Vres1, Vres, Vres, 0);
            if (i == len || (must_norm && !bitval(g, len - i - 1)))
                force_normalize(Vres);
            gwsub(gwdata, V, Vres);
            gwsetaddin(gwdata, -2);
            gwsquare2(gwdata, Vres1, Vres1, i < len ? GWMUL_STARTNEXTFFT : 0);
        }
        else
        {
            gwmul3(gwdata, Vres, Vres1, Vres1, 0);
            if (i == len || (must_norm && bitval(g, len - i - 1)))
                force_normalize(Vres1);
            gwsub(gwdata, V, Vres1);
            gwsetaddin(gwdata, -2);
            gwsquare2(gwdata, Vres, Vres, i < len ? GWMUL_STARTNEXTFFT : 0);
        }
        gwsetaddin(gwdata, 0);
    }
}

void lucas_V_mul_int(int x, gwnum V, gwnum Vres, gwnum Vres1)
{
    giant g = allocgiant(1);
    itog(x, g);
    lucas_V_mul_giant(g, V, Vres, Vres1);
    free(g);
}

void lucas_V_shiftleft(int x, gwnum V)
{
    gwsetaddin(gwdata, -2);
    for (; x > 0; x--, costAdd(1))
        gwsquare(gwdata, V);
    gwsetaddin(gwdata, 0);
}

void lucas_V_inc(gwnum Vinc, gwnum Vprev, gwnum Vcur, gwnum Vnext, int options)
{
    gwmul3(gwdata, Vinc, Vcur, Vnext, GWMUL_FFT_S1 |options);
    force_normalize(Vnext);
    gwsub(gwdata, Vprev, Vnext);
    costAdd(1);
}

#define lucas_V_mul_2(V) lucas_V_shiftleft(1, V)
#define lucas_V_add(Va, Vb, Vbma, Vbpa, options) lucas_V_inc(Va, Vbma, Vb, Vbpa, options)

