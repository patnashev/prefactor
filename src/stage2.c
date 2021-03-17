
// Get distances to multiples of D, look ahead for L sections.
// Simple pairing with the first match.
short *get_distances_simple(int *primes, int B1, int B2, int D, int L)
{
    int i, j;
    char* map = malloc(B2);
    memset(map, 0, B2);
    for (j = 0; primes[j] <= B1; j++);
    int* plist = primes + j;
    for (j = 0; plist[j] <= B2; j++)
        map[plist[j]] = 1;

    short* distances = malloc(sizeof(short)*(j + 1));
    int d, p;
    for (j = 0; plist[j] <= B2; j++)
    {
        p = plist[j];
        distances[j] = 0;
        if (map[p] != 1)
            continue;
        d = p - p%D + D;
        short distance = d - p;
        for (i = 1; i <= L; i++)
        {
            d = p - p%D + i*D;
            if (2*d - p <= B2 && map[2*d - p] == 1)
            {
                distance = d - p;
                map[2*d - p] = 0;
                break;
            }
        }
        distances[j] = distance;
    }
    distances[j] = -1;

    free(map);
    return distances;
}

// Get distances to multiples of D and i*D + D/A, look ahead for L sections.
// Trying to find the best pairing.
short* get_distances_A(int *primes, int B1, int B2, int D, int A, int L)
{
    int i, j, k;
    short* bases = malloc(sizeof(short)*(L + 1)*2);
    for (i = 0; i <= L; i++)
    {
        bases[2*i] = (i + 1)*D;
        bases[2*i + 1] = (i + 1)*D + D/A;
    }
    short *map = malloc(sizeof(short)*B2);
    memset(map, 0, sizeof(short)*B2);
    for (j = 0; primes[j] <= B1; j++);
    int* plist = primes + j;
    for (j = 0; plist[j] <= B2; j++)
        map[plist[j]] = 1;
    int ptotal = j;

    int d, p;
    for (j = 0; plist[j] <= B2; j++)
    {
        p = plist[j];
        if (map[p] != 1)
            continue;
        for (i = 0; i < 2*L; i++)
        {
            d = p - p%D + bases[i];
            if (2*d - p <= B2 && map[2*d - p] == 1 && d - p < 16000)
            {
                map[p] = 2*d - 2*p;
                map[2*d - p] = 2*p - 2*d;
                break;
            }
        }
    }

    #define qlen 1000
    int queue[qlen];
    int qstart = 0, qend = 0;
    short *link = malloc(sizeof(short)*B2);
    k = 0;
    int kk = 0;
    int pairs = 0;
    int flag = 1;
    while (flag)
    {
        flag = 0;
        kk += k;
        memset(link, 0, sizeof(short)*B2);
        for (k = 0; plist[k] <= B2 && !flag; k++)
        {
            p = plist[(k + kk)%ptotal];
            if (map[p] != 1)
                continue;
            qstart = 0; qend = 1; queue[0] = p;
            link[p] = 1;

            while (qstart != qend)
            {
                qstart = qstart%qlen + 1;
                p = queue[qstart - 1];
                for (i = 0; i < 2*2*L + 2; i++)
                {
                    d = p - p%D + (i < 2*L ? bases[i] : bases[i - 2*L] - (L + 1)*D);
                    if (abs(d - p) >= D*L + D/A || abs(d - p) > 16000)
                        continue;
                    j = 2*d - p;
                    if (j <= B1 || j > B2 || map[j] == 0 || link[j] != 0)
                        continue;
                    if (map[j] == 1)
                    {
                        while (1)
                        {
                            map[j] = p - j;
                            map[p] = j - p;
                            if (link[p] == 1)
                                break;
                            j = p + link[p];
                            p = j + link[j];
                            if (p < 0)
                                printf("error");
                        }
                        qstart = qend;
                        pairs++;
                        //if (pairs%100 == 0)
                        //    printf("%d\n", pairs);
                        flag = 1;
                        break;
                    }

                    link[j] = p - j;
                    j += map[j];
                    link[j] = map[j];
                    if (j < B1 || j > B2)
                        printf("error");
                    qend = qend%qlen + 1;
                    queue[qend - 1] = j;
                    if (qstart == qend)
                        printf("error");
                }
            }
        }
    }
    free(link);

    short* distances = malloc(sizeof(short)*(ptotal + 1));
    for (j = 0; plist[j] <= B2; j++)
    {
        p = plist[j];
        if (map[p] < 0)
            distances[j] = 0;
        else if (map[p] == 1)
            distances[j] = D - p%D;
        else
            distances[j] = map[p]/2;
    }
    distances[j] = -1;
    free(map);
    free(bases);
    return distances;
}

// P+-1 factoring stage 2.
void do_pm1stage2(int B1, int B2, giant P, int minus1, int D, int A, int L)
{
    int i, j;
    int p, it;

    if (B2 <= B1)
        return;
    for (it = 0; primes[it] <= B1; it++);
    p = primes[it++];

    printf("%s, P%c1 stage 2, B2 = %d.\n", Nstr, minus1 ? '-' : '+', B2);
    double timer = getHighResTimer();
    int transforms = -(int)gwdata->fft_count;

    // p = (v+1)*D - u
    if (A > 1)
    {
        D /= A;
        L = L*A + 1;
    }

    for (i = 0; primes[i]*B1 < B2; i++)
        if (D%primes[i] != 0)
        {
            it--;
            for (j = it; primes[j]*primes[i] <= B2; j++)
                while (primes[j]*primes[i] <= B2)
                    primes[j] *= primes[i];
            for (; primes[j] <= B2; j++);
            qsort(primes + it, j - it, sizeof(int), cmp);
            p = primes[it++];
            break;
        }

    gwnum Vn = gwalloc(gwdata);
    gwnum Vn1 = gwalloc(gwdata);
    gwnum W = gwalloc(gwdata);
    gwnum G = gwalloc(gwdata);
    gwnum Vtmp = gwalloc(gwdata);
    gwnum *precomp_V = malloc(sizeof(gwnum)*D/2*L);
    memset(precomp_V, 0, sizeof(gwnum)*D/2*L);
    gwnum* VL = malloc(sizeof(gwnum)*(L/A + 1));
    gwnum Va = NULL;
    gwnum Va1 = NULL;
    gwnum* VA = NULL;
    if (A > 1)
    {
        Va = gwalloc(gwdata);
        Va1 = gwalloc(gwdata);
        VA = malloc(sizeof(gwnum)*(L/A + 1));
    }

    // Precomputation of V_u(V_n) where gcd(u,D)=1
    gianttogw(gwdata, P, Vn);
    precomp_V[0] = gwalloc(gwdata);
    gwcopy(gwdata, Vn, precomp_V[0]);
    lucas_V_mul_2(Vn); // V_2
    gwcopy(gwdata, Vn, Vn1);
    lucas_V_mul_2(Vn1); // V_4
    gwfft(gwdata, Vn, Vtmp);
    int dist = 1;
    int precomp = 1;
    for (i = 1; i < D/2*L; i++)
    {
        if (gcd(2*i + 1, 2*dist) != 1)
            continue;
        // V_{2i+1}
        precomp_V[i] = gwalloc(gwdata);
        precomp++;
        lucas_V_inc(Vtmp, precomp_V[i >= 2*dist ? i - 2*dist : 2*dist - i - 1], precomp_V[i - dist], precomp_V[i], GWMUL_PRESERVE_S2);
        if ((i + 1)%dist == 0 && D%(i + 1) == 0 && gcd((i + 1)/dist, 2*dist) == 1)
        {
            int pd = (i + 1)/dist;
            gwswap(Vn, Vtmp);
            lucas_V_add(precomp_V[dist*(pd - 1)/2 - 1], precomp_V[dist*(pd + 1)/2], Vn1, Vn, GWMUL_PRESERVE_S1 | GWMUL_PRESERVE_S2);
            lucas_V_add(precomp_V[dist*(pd - 1)/2], precomp_V[dist*(pd + 1)/2], Vtmp, Vn1, GWMUL_PRESERVE_S1 | GWMUL_PRESERVE_S2);
            gwfft(gwdata, Vn, Vtmp);
            dist = i + 1;
            for (j = 0; j < i; j++)
                if (precomp_V[j] != NULL && gcd(2*j + 1, 2*dist) != 1)
                {
                    gwfree(gwdata, precomp_V[j]);
                    precomp_V[j] = NULL;
                    precomp--;
                }
        }
    }
    // V_{D/A}
    if (A > 1)
    {
        gwcopy(gwdata, Vn, Va);
        if (D/2/dist != 1)
        {
            int pd = D/2/dist;
            for (j = 0; (pd & 1) == 0; pd >>= 1, j++);
            if (pd > 1)
                lucas_V_add(precomp_V[dist*(pd - 1)/2 - 1], precomp_V[dist*(pd + 1)/2], Vn1, Va, 0);
            lucas_V_shiftleft(j, Va);
        }
        D *= A;
        L /= A;
    }
    // W = V_D
    gwswap(W, Vn);
    if (D/2/dist != 1)
    {
        int pd = D/2/dist;
        for (j = 0; (pd & 1) == 0; pd >>= 1, j++);
        if (pd > 1)
            lucas_V_add(precomp_V[dist*(pd - 1)/2 - 1], precomp_V[dist*(pd + 1)/2], Vn1, W, 0);
        lucas_V_shiftleft(j, W);
    }

    int v = p/D;

    if (A > 1)
    {
        gwswap(G, Va);
        dbltogw(gwdata, 2, Va);
        gwcopy(gwdata, G, Va1);
        if (v > 0)
            lucas_V_mul_int(v*A, G, Va, Va1);
        gwswap(Vn1, Va);
        gwcopy(gwdata, Vn1, Vn);
        gwcopy(gwdata, Va1, Va);
        for (i = 0; i < A; i++)
        {
            lucas_V_inc(G, Vn1, Va1, Vtmp, GWMUL_PRESERVE_S2);
            gwswap(Vn1, Va1);
            gwswap(Va1, Vtmp);
        }

        VA[0] = gwalloc(gwdata);
        gwfft(gwdata, Va, VA[0]);
        VA[1] = gwalloc(gwdata);
        gwfft(gwdata, Va1, VA[1]);
        for (i = 2; i <= L; i++)
        {
            lucas_V_inc(W, Va, VA[i - 1], Vtmp, 0);
            gwswap(Va, Va1);
            gwswap(Va1, Vtmp);
            VA[i] = gwalloc(gwdata);
            gwfft(gwdata, Va1, VA[i]);
        }
    }
    else
    {
        // V_{v*D} V_{(v+1)*D}
        dbltogw(gwdata, 2, Vn);
        gwcopy(gwdata, W, Vn1);
        if (v > 0)
            lucas_V_mul_int(v, W, Vn, Vn1);
    }

    VL[0] = gwalloc(gwdata);
    gwfft(gwdata, Vn1, VL[0]);
    for (i = 1; i <= L; i++)
    {
        lucas_V_inc(W, Vn, VL[i - 1], Vtmp, 0);
        gwswap(Vn, Vn1);
        gwswap(Vn1, Vtmp);
        VL[i] = gwalloc(gwdata);
        gwfft(gwdata, Vn1, VL[i]);
    }

    int paired = 0;
    short* distances;
    if (A > 1)
        distances = get_distances_A(primes, B1, B2, D, A, L);
    else
        distances = get_distances_simple(primes, B1, B2, D, L);
    for (j = 0; distances[j] != -1; j++)
        if (distances[j] == 0)
            paired++;
    transforms += (int)gwdata->fft_count;
    printf("%d precomputed values (%d transforms), %d%% pairing, D = %d, L = %d", precomp, transforms, (200*paired + j/2)/j, D, L);
    if (A > 1)
        printf(", A = %d.\n", A);
    else
        printf(".\n");

    transforms = -(int)gwdata->fft_count;
    dbltogw(gwdata, 1, G);
    int last;
    for (last = j - 1; distances[last] == 0; last--);
    for (j = 0; j <= last; j++)
    {
        while (v < p/D)
        {
            // n += w
            lucas_V_inc(W, Vn, VL[L], Vtmp, 0);
            gwswap(Vn, Vn1);
            gwswap(Vn1, Vtmp);
            for (i = 1; i <= L; i++)
                gwswap(VL[i - 1], VL[i]);
            gwfft(gwdata, Vn1, VL[L]);
            // nA += w
            if (A > 1)
            {
                lucas_V_inc(W, Va, VA[L], Vtmp, 0);
                gwswap(Va, Va1);
                gwswap(Va1, Vtmp);
                for (i = 1; i <= L; i++)
                    gwswap(VA[i - 1], VA[i]);
                gwfft(gwdata, Va1, VA[L]);
            }
            v++;
        }

        if (distances[j] != 0)
        {
            int d = distances[j] + p%D;
            gwsubmul4(gwdata, d%D == 0 ? VL[d/D - 1] : VA[d/D], precomp_V[distances[j]/2], G, G, j < last ? GWMUL_STARTNEXTFFT : 0);
            costAdd(1);
        }

        p = primes[it++];
    }

    giant tmp = getg();
    gwtogiant(gwdata, G, tmp);
    gcdg(N, tmp);

    timer = (getHighResTimer() - timer)/getHighResTimerFrequency();
    transforms += (int)gwdata->fft_count;
    if (isZero(tmp) || gcompg(tmp, N) == 0)
    {
        printf("All divisors of N < B2.\n");
    }
    else if (!isone(tmp))
    {
        report_factor(tmp);
    }
    else
    {
        printf("No factors found, transforms: %d, time: %d s.\n", transforms, (int)timer);
    }

    free(distances);
    for (i = 0; i < D/2; i++)
        if (precomp_V[i] != NULL)
            gwfree(gwdata, precomp_V[i]);
    free(precomp_V);
    for (i = 0; i <= L; i++)
        gwfree(gwdata, VL[i]);
    free(VL);
    if (A > 1)
    {
        for (i = 0; i <= L; i++)
            gwfree(gwdata, VA[i]);
        free(VA);
        gwfree(gwdata, Va);
        gwfree(gwdata, Va1);
    }
    gwfree(gwdata, G);
    gwfree(gwdata, Vtmp);
    gwfree(gwdata, Vn);
    gwfree(gwdata, Vn1);
    gwfree(gwdata, W);
    freeg();
}

// EdECM factoring stage 2.
void do_edecm_stage2(int B1, int B2, ed_point P, int D, int A, int L)
{
    int i, j;
    int p, it;

    if (B2 <= B1)
        return;
    for (it = 0; primes[it] <= B1; it++);
    p = primes[it++];

    printf("%s, EdECM stage 2, B2 = %d.\n", Nstr, B2);
    double timer = getHighResTimer();
    int transforms = -(int)gwdata->fft_count;

    // p = (v+1)*D - u
    if (A > 1)
    {
        D /= A;
        L = L*A + 1;
    }

    for (i = 0; primes[i]*B1 < B2; i++)
        if (D%primes[i] != 0)
        {
            it--;
            for (j = it; primes[j]*primes[i] <= B2; j++)
                while (primes[j]*primes[i] <= B2)
                    primes[j] *= primes[i];
            for (; primes[j] <= B2; j++);
            qsort(primes + it, j - it, sizeof(int), cmp);
            p = primes[it++];
            break;
        }

    gwnum G = gwalloc(gwdata);
    gwnum TG = gwalloc(gwdata);
    ed_point Pn = ed_alloc();
    ed_point Pn1 = ed_alloc();
    ed_point W = ed_alloc();
    ed_point Ptmp = ed_alloc();
    int precomp_size = D/2*L + 2;
    ed_point *precomp_P = malloc(sizeof(ed_point)*precomp_size);
    memset(precomp_P, 0, sizeof(ed_point)*precomp_size);
    ed_point* PL = malloc(sizeof(ed_point)*(L/A + 1));
    ed_point Pa = NULL;
    ed_point Pa1 = NULL;
    ed_point Wa = NULL;
    ed_point* PA = NULL;
    if (A > 1)
    {
        Pa = ed_alloc();
        Pa1 = ed_alloc();
        Wa = ed_alloc();
        PA = malloc(sizeof(ed_point)*(L/A + 1));
    }

    // Precomputation of V_u(V_n) where gcd(u,D)=1
    ed_copy(P, Pn);
    precomp_P[0] = ed_alloc();
    ed_copy(Pn, precomp_P[0]);
    ed_y_mul2(Pn); // V_2
    ed_copy(Pn, Pn1);
    ed_y_mul2(Pn1); // V_4
    ed_y_optimize(Pn, Ptmp);
    int dist = 1;
    int precomp = 1;
    for (i = 1; i < D/2*L; i++)
    {
        if (gcd(2*i + 1, 2*dist) != 1)
            continue;
        // V_{2i+1}
        precomp_P[i] = ed_alloc();
        precomp++;
        ed_y_inc(Ptmp, precomp_P[i >= 2*dist ? i - 2*dist : 2*dist - i - 1], precomp_P[i - dist], precomp_P[i]);
        if ((i + 1)%dist == 0 && D%(i + 1) == 0 && gcd((i + 1)/dist, 2*dist) == 1)
        {
            int pd = (i + 1)/dist;
            ed_swap(Ptmp, Pn);
            ed_y_add(precomp_P[dist*(pd - 1)/2 - 1], precomp_P[dist*(pd + 1)/2], Pn1, Pn);
            ed_y_add(precomp_P[dist*(pd - 1)/2], precomp_P[dist*(pd + 1)/2], Ptmp, Pn1);
            ed_y_optimize(Pn, Ptmp);
            dist = i + 1;
            for (j = 0; j < i; j++)
                if (precomp_P[j] != NULL && gcd(2*j + 1, 2*dist) != 1)
                {
                    ed_free(precomp_P[j]);
                    precomp_P[j] = NULL;
                    precomp--;
                }
        }
    }
    // V_{D/A}
    /*if (A > 1)
    {
        gwcopy(gwdata, Vn, Va);
        if (D/2/dist != 1)
        {
            int pd = D/2/dist;
            for (j = 0; (pd & 1) == 0; pd >>= 1, j++);
            if (pd > 1)
                lucas_V_add(precomp_V[dist*(pd - 1)/2 - 1], precomp_V[dist*(pd + 1)/2], Vn1, Va, 0);
            lucas_V_shiftleft(j, Va);
        }
        D *= A;
        L /= A;
    }*/
    // W = V_D
    /*ed_swap(W, Pn);
    if (D/2/dist != 1)
    {
        int pd = D/2/dist;
        for (j = 0; (pd & 1) == 0; pd >>= 1, j++);
        if (pd > 1)
            ed_y_add(precomp_P[dist*(pd - 1)/2 - 1], precomp_P[dist*(pd + 1)/2], Pn1, W);
        ed_y_shiftleft(j, W, GWMUL_STARTNEXTFFT1);
    }*/
    if (A > 1)
    {
        ed_mul_int_add(D, P, Wa, NULL);
        precomp_P[precomp_size - 2] = Wa;
        ed_mul_int_add(A, Wa, W, NULL);
        D *= A;
        L /= A;
    }
    else
        ed_mul_int_add(D, P, W, NULL);
    /*{
        ed_normalize(W);
        ed_point T = ed_alloc();
        ed_mul_int_add(D, P, T, NULL);
        ed_normalize(T);
        giant t1 = getg();
        gwunfft(gwdata, W->Y, W->Y);
        gwtogiant(gwdata, W->Y, t1);
        giant t2 = getg();
        gwunfft(gwdata, T->Y, T->Y);
        gwtogiant(gwdata, T->Y, t2);
        printf("1");
    }*/
    precomp_P[precomp_size - 1] = W;
    ed_normalize_pool(precomp_P, precomp_size);
    precomp_P[precomp_size - 1] = NULL;
    precomp_P[precomp_size - 2] = NULL;

    int v = p/D;

    if (A > 1)
    {
        ed_zero(Pa);
        ed_copy(Wa, Pa1);
        if (v > 0)
            ed_mul_int_add(v*A, Wa, Pa, Pa1);
        ed_swap(Pn1, Pa);
        ed_copy(Pn1, Pn);
        ed_copy(Pa1, Pa);
        ed_add(W, Pn1, GWMUL_STARTNEXTFFT1);
        ed_add(W, Pa1, GWMUL_STARTNEXTFFT1);

        PA[0] = ed_alloc();
        ed_y_optimize(Pa, PA[0]);
        PA[1] = ed_alloc();
        ed_y_optimize(Pa1, PA[1]);
        for (i = 2; i <= L; i++)
        {
            ed_y_inc(W, Pa, PA[i - 1], Ptmp);
            ed_swap(Pa, Pa1);
            ed_swap(Pa1, Ptmp);
            PA[i] = ed_alloc();
            ed_y_optimize(Pa1, PA[i]);
        }
    }
    else
    {
        // V_{v*D} V_{(v+1)*D}
        ed_zero(Pn);
        ed_copy(W, Pn1);
        if (v > 0)
            ed_mul_int_add(v, W, Pn, Pn1);
    }

    PL[0] = ed_alloc();
    ed_copy(Pn1, PL[0]);
    for (i = 1; i <= L; i++)
    {
        ed_y_inc(W, Pn, PL[i - 1], Ptmp);
        ed_swap(Pn, Pn1);
        ed_swap(Pn1, Ptmp);
        PL[i] = ed_alloc();
        ed_y_optimize(Pn1, PL[i]);
    }

    int paired = 0;
    short* distances;
    if (A > 1)
        distances = get_distances_A(primes, B1, B2, D, A, L);
    else
        distances = get_distances_simple(primes, B1, B2, D, L);
    for (j = 0; distances[j] != -1; j++)
        if (distances[j] == 0)
            paired++;
    transforms += (int)gwdata->fft_count;
    printf("%d precomputed values (%d transforms), %d%% pairing, D = %d, L = %d", precomp, transforms, (200*paired + j/2)/j, D, L);
    if (A > 1)
        printf(", A = %d.\n", A);
    else
        printf(".\n");

    transforms = -(int)gwdata->fft_count;
    dbltogw(gwdata, 1, G);
    int last;
    for (last = j - 1; distances[last] == 0; last--);
    for (j = 0; j <= last; j++)
    {
        while (v < p/D)
        {
            // n += w
            ed_y_inc(W, Pn, Pn1, Ptmp);
            ed_swap(Pn, Pn1);
            ed_swap(Pn1, Ptmp);
            for (i = 1; i <= L; i++)
                ed_swap(PL[i - 1], PL[i]);
            ed_y_optimize(Pn1, PL[L]);
            // nA += w
            if (A > 1)
            {
                ed_y_inc(W, Pa, Pa1, Ptmp);
                ed_swap(Pa, Pa1);
                ed_swap(Pa1, Ptmp);
                for (i = 1; i <= L; i++)
                    ed_swap(PA[i - 1], PA[i]);
                ed_y_optimize(Pa1, PA[L]);
            }
            v++;
        }

        if (distances[j] != 0)
        {
            int d = distances[j] + p%D;
            /*{
                //ed_normalize(W);
                ed_point T = ed_alloc();
                ed_mul_int_add(distances[j], P, T, NULL);
                ed_normalize(T);
                giant t1 = getg();
                gwunfft(gwdata, precomp_P[distances[j]/2]->Y, precomp_P[distances[j]/2]->Y);
                gwtogiant(gwdata, precomp_P[distances[j]/2]->Y, t1);
                giant t2 = getg();
                gwunfft(gwdata, T->Y, T->Y);
                gwtogiant(gwdata, T->Y, t2);
                printf("1");
            }*/
            ed_point PD = d%D == 0 ? PL[d/D - 1] : PA[d/D];
            if (PD->Z != NULL)
            {
                gwmul3(gwdata, precomp_P[distances[j]/2]->Y, PD->Z, TG, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT);
                gwsubmul4(gwdata, PD->Y, TG, G, G, j < last ? GWMUL_STARTNEXTFFT : 0);
            }
            else
                gwsubmul4(gwdata, PD->Y, precomp_P[distances[j]/2]->Y, G, G, j < last ? GWMUL_STARTNEXTFFT : 0);
            costAdd(1);
        }

        p = primes[it++];
    }

    giant tmp = getg();
    gwtogiant(gwdata, G, tmp);
    gcdg(N, tmp);

    timer = (getHighResTimer() - timer)/getHighResTimerFrequency();
    transforms += (int)gwdata->fft_count;
    if (isZero(tmp) || gcompg(tmp, N) == 0)
    {
        printf("All divisors of N < B2.\n");
    }
    else if (!isone(tmp))
    {
        report_factor(tmp);
    }
    else
    {
        printf("No factors found, transforms: %d, time: %d s.\n", transforms, (int)timer);
    }

    free(distances);
    for (i = 0; i < precomp_size; i++)
        if (precomp_P[i] != NULL)
            ed_free(precomp_P[i]);
    free(precomp_P);
    for (i = 0; i <= L; i++)
        ed_free(PL[i]);
    free(PL);
    if (A > 1)
    {
        for (i = 0; i <= L; i++)
            ed_free(PA[i]);
        free(PA);
        ed_free(Wa);
        ed_free(Pa);
        ed_free(Pa1);
    }
    ed_free(Ptmp);
    ed_free(Pn);
    ed_free(Pn1);
    ed_free(W);
    gwfree(gwdata, TG);
    gwfree(gwdata, G);
    freeg();
}