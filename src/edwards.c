
// Edwards curve point in extended coordinates with several temporary values
typedef struct ed_point_struct
{
    gwnum X;
    gwnum Y;
    gwnum Z;
    gwnum T;
    gwnum YpX;
    gwnum YmX;
    gwnum ZpY;
    gwnum ZmY;
} *ed_point;

// Allocs Edwards curve point.
ed_point ed_alloc();
// Frees Edwards curve point.
void ed_free(ed_point p);
// Copies Edwards curve point.
void ed_copy(ed_point s, ed_point d);
// Swaps two Edwards curve points.
#define ed_swap(s, d) {ed_point t; t = s; s = d; d = t;}
// Zeroes Edwards curve point.
void ed_zero(ed_point p);
// Sets Edwards curve point to (xa/xb, ya/yb).
void ed_from_small(int xa, int xb, int ya, int yb, ed_point p);
// Doubles Edwards curve point.
void ed_mul2(ed_point p, int options);
// Doubles Edwards curve point using _carefully operations.
void ed_mul2_carefully(ed_point p);
// Adds two Edwards curve points.
void ed_add(ed_point padd, ed_point p, int options);
// Doubles Edwards curve point and adds another one to it.
void ed_mul2_add(ed_point padd, ed_point p, int options);
// Doubles Edwards curve point and adds a numerically small (xa/xb, ya/yb) to it.
void ed_mul2_addsmall(int xa, int xb, int ya, int yb, ed_point p, int options);
// Option for signed adds. ED_SIGNED_ADD may be necessary if followed by ED_SIGNED_SUB.
#define ED_SIGNED_ADD 0x200000
#define ED_SIGNED_SUB 0x400000
// Option to get the full set of coordinates. Necessary if the result is used for addition.
#define ED_MUL_FULL 0x800000
// Pooled normalization.
giant ed_normalize_pool(ed_point *pool, int count);
// Normalizes Edwards curve point.
#define ed_normalize(p) ed_normalize_pool(&(p), 1)
// Multiplies Edwards curve point by a giant scalar.
void ed_mul_giant_add(giant g, ed_point p, ed_point pres, ed_point pres1);
// Multiplies Edwards curve point by a scalar.
void ed_mul_int_add(int x, ed_point p, ed_point pres, ed_point pres1);
// Multiplies Edwards curve point by a NAF_w.
giant ed_mul_nafw(int W, short *nafw, int len, ed_point P);
// Doubles y coordinate of Edwards curve point.
void ed_y_mul2(ed_point p);
// Multiplies y coordinate of Edwards curve point by power of 2.
void ed_y_shiftleft(int x, ed_point p);
// Increments y coordinate of Edwards curve point in sequence.
void ed_y_inc(ed_point pinc, ed_point pprev, ed_point pcur, ed_point pnext);
// Adds y coordinates of two Edwards curve points with known difference.
#define ed_y_add(pa, pb, pbma, pbpa) ed_y_inc(pa, pbma, pb, pbpa)
// Optimizes memory use of Edwards curve point.
void ed_y_optimize(ed_point s, ed_point d);
// Returns j-invariant = 16*(1 + 14*d + d^2)^3/(d*(1 - d)^4)
void get_j_invariant(giant res);
// Generates Edwards curve and returns a non-torsion point.
int gen_ed_curve(int seed, ed_point P);

ed_point ed_alloc()
{
    ed_point p = malloc(sizeof(struct ed_point_struct));
    p->X = NULL;
    p->Y = NULL;
    p->Z = NULL;
    p->T = NULL;
    p->YpX = NULL;
    p->YmX = NULL;
    p->ZpY = NULL;
    p->ZmY = NULL;
    return p;
}

void ed_free(ed_point p)
{
    ed_zero(p);
    free(p);
}

void ed_zero(ed_point p)
{
    if (p->X != NULL)
        gwfree(gwdata, p->X);
    if (p->Y != NULL)
        gwfree(gwdata, p->Y);
    if (p->Z != NULL)
        gwfree(gwdata, p->Z);
    if (p->T != NULL)
        gwfree(gwdata, p->T);
    if (p->YpX != NULL)
        gwfree(gwdata, p->YpX);
    if (p->YmX != NULL)
        gwfree(gwdata, p->YmX);
    if (p->ZpY != NULL)
        gwfree(gwdata, p->ZpY);
    if (p->ZmY != NULL)
        gwfree(gwdata, p->ZmY);
    p->X = NULL;
    p->Y = NULL;
    p->Z = NULL;
    p->T = NULL;
    p->YpX = NULL;
    p->YmX = NULL;
    p->ZpY = NULL;
    p->ZmY = NULL;
}

void ed_copy(ed_point s, ed_point d)
{
    if (s->X != NULL)
    {
        if (d->X == NULL)
            d->X = gwalloc(gwdata);
        gwcopy(gwdata, s->X, d->X);
    }
    else if (d->X != NULL)
    {
        gwfree(gwdata, d->X);
        d->X = NULL;
    }
    if (s->Y != NULL)
    {
        if (d->Y == NULL)
            d->Y = gwalloc(gwdata);
        gwcopy(gwdata, s->Y, d->Y);
    }
    else if (d->Y != NULL)
    {
        gwfree(gwdata, d->Y);
        d->Y = NULL;
    }
    if (s->Z != NULL)
    {
        if (d->Z == NULL)
            d->Z = gwalloc(gwdata);
        gwcopy(gwdata, s->Z, d->Z);
    }
    else if (d->Z != NULL)
    {
        gwfree(gwdata, d->Z);
        d->Z = NULL;
    }
    if (s->T != NULL)
    {
        if (d->T == NULL)
            d->T = gwalloc(gwdata);
        gwcopy(gwdata, s->T, d->T);
    }
    else if (d->T != NULL)
    {
        gwfree(gwdata, d->T);
        d->T = NULL;
    }
    if (s->YpX != NULL)
    {
        if (d->YpX == NULL)
            d->YpX = gwalloc(gwdata);
        gwcopy(gwdata, s->YpX, d->YpX);
    }
    else if (d->YpX != NULL)
    {
        gwfree(gwdata, d->YpX);
        d->YpX = NULL;
    }
    if (s->YmX != NULL)
    {
        if (d->YmX == NULL)
            d->YmX = gwalloc(gwdata);
        gwcopy(gwdata, s->YmX, d->YmX);
    }
    else if (d->YmX != NULL)
    {
        gwfree(gwdata, d->YmX);
        d->YmX = NULL;
    }
    if (s->ZpY != NULL)
    {
        if (d->ZpY == NULL)
            d->ZpY = gwalloc(gwdata);
        gwcopy(gwdata, s->ZpY, d->ZpY);
    }
    else if (d->ZpY != NULL)
    {
        gwfree(gwdata, d->ZpY);
        d->ZpY = NULL;
    }
    if (s->ZmY != NULL)
    {
        if (d->ZmY == NULL)
            d->ZmY = gwalloc(gwdata);
        gwcopy(gwdata, s->ZmY, d->ZmY);
    }
    else if (d->ZmY != NULL)
    {
        gwfree(gwdata, d->ZmY);
        d->ZmY = NULL;
    }
}

void ed_from_small(int xa, int xb, int ya, int yb, ed_point p)
{
    giant tmp = getg();
    giant tmp2 = getg();
    giant tmp3 = getg();

    if (p->X == NULL)
        p->X = gwalloc(gwdata);
    if (p->Y == NULL)
        p->Y = gwalloc(gwdata);
    if (p->Z == NULL)
        p->Z = gwalloc(gwdata);
    if (p->T == NULL)
        p->T = gwalloc(gwdata);
    if (p->YpX == NULL)
        p->YpX = gwalloc(gwdata);
    if (p->YmX == NULL)
        p->YmX = gwalloc(gwdata);
    if (p->ZpY != NULL)
        gwfree(gwdata, p->ZpY);
    p->ZpY = NULL;
    if (p->ZmY != NULL)
        gwfree(gwdata, p->ZmY);
    p->ZmY = NULL;
    if (EdD == NULL)
        EdD = gwalloc(gwdata);

    if (xb < 0)
    {
        xa = -xa;
        xb = -xb;
    }
    if (yb < 0)
    {
        ya = -ya;
        yb = -yb;
    }
    ultog(abs(xa), tmp);
    ulmulg(abs(xa), tmp);
    ultog(xb, tmp2);
    ulmulg(xb, tmp2);
    subg(tmp2, tmp);
    ultog(yb, tmp2);
    ulmulg(yb, tmp2);
    mulg(tmp2, tmp);
    ultog(abs(ya), tmp2);
    ulmulg(abs(ya), tmp2);
    ulmulg(xb, tmp2);
    ulmulg(xb, tmp2);
    addg(tmp2, tmp);
    if (tmp->sign < 0)
        addg(N, tmp);
    gianttogw(gwdata, tmp, EdD);
    ultog(abs(xa), tmp2);
    ulmulg(abs(xa), tmp2);
    ulmulg(abs(ya), tmp2);
    ulmulg(abs(ya), tmp2);
    invg(N, tmp2);
    if (tmp2->sign < 0)
    {
        gwfree(gwdata, EdD);
        EdD = NULL;
    }
    gianttogw(gwdata, tmp2, p->T);
    gwmul_carefully(gwdata, p->T, EdD);

    itog(xa, tmp);
    ulmulg(yb, tmp);
    gtog(tmp, tmp3);
    if (tmp->sign < 0)
        addg(N, tmp);
    gianttogw(gwdata, tmp, p->X);
    itog(ya, tmp2);
    ulmulg(xb, tmp2);
    addg(tmp2, tmp);
    if (tmp2->sign < 0)
        addg(N, tmp2);
    gianttogw(gwdata, tmp2, p->Y);
    if (tmp->sign < 0)
        addg(N, tmp);
    gianttogw(gwdata, tmp, p->YpX);
    subg(tmp3, tmp2);
    if (tmp2->sign < 0)
        addg(N, tmp2);
    gianttogw(gwdata, tmp2, p->YmX);
    ultog(xb, tmp);
    ulmulg(yb, tmp);
    gianttogw(gwdata, tmp, p->Z);
    ultog(abs(xa), tmp);
    ulmulg(abs(ya), tmp);
    gianttogw(gwdata, tmp, p->T);

    freeg();
    freeg();
}

void ed_mul2(ed_point p, int options)
{
    if (p->T == NULL)
        p->T = gwalloc(gwdata);
    if (p->YpX != NULL)
        gwfree(gwdata, p->YpX);
    p->YpX = NULL;
    if (p->YmX != NULL)
        gwfree(gwdata, p->YmX);
    p->YmX = NULL;

    gwsetmulbyconst(gwdata, 2);
    gwmul3(gwdata, p->X, p->Y, p->T, GWMUL_FFT_S1 | GWMUL_MULBYCONST | GWMUL_STARTNEXTFFT); // E = 2X1Y1
    if (p->Z != NULL)
    {
        gwsquare2(gwdata, p->Z, p->Z, GWMUL_MULBYCONST | GWMUL_STARTNEXTFFT); // C = 2Z1^2
        costAdd(1);
    }
    else
    {
        p->Z = gwalloc(gwdata);
        dbltogw(gwdata, 2, p->Z);
    }
    gwsquare(gwdata, p->X); // A = X1^2
    gwsquare(gwdata, p->Y); // B = Y1^2
    force_normalize(p->Y);
    gwaddsub(gwdata, p->Y, p->X); // G = B + A, H = B - A
    //gwsub(gwdata, p->Y, p->Z); // F = C - G
    gwnum tmp = NULL;
    if (options & ED_MUL_FULL)
    {
        tmp = gwalloc(gwdata);
        gwmul3(gwdata, p->X, p->T, tmp, GWMUL_FFT_S1 | (options & (~ED_MUL_FULL))); // T3 = E * H
    }
    //gwmul3(gwdata, p->Z, p->T, p->T, options & (~ED_MUL_FULL)); // X3 = E * F
    //gwmul3(gwdata, p->Y, p->Z, p->Z, options & (~ED_MUL_FULL)); // Z3 = F * G
    gwsubmul4(gwdata, p->Z, p->Y, p->T, p->T, options & (~ED_MUL_FULL)); // F = C - G, X3 = F * E
    gwsubmul4(gwdata, p->Z, p->Y, p->Y, p->Z, options & (~ED_MUL_FULL)); // F = C - G, Z3 = F * G
    gwmul3(gwdata, p->X, p->Y, p->Y, options & (~ED_MUL_FULL)); // Y3 = G * H
    gwswap(p->X, p->T);
    if (options & ED_MUL_FULL)
    {
        gwswap(p->T, tmp);
        costAdd(1);
        gwfree(gwdata, tmp);
    }
    else
    {
        gwfree(gwdata, p->T);
        p->T = NULL;
    }
    costAdd(6);
}

void ed_mul2_carefully(ed_point p)
{
    if (p->T == NULL)
        p->T = gwalloc(gwdata);
    if (p->YpX != NULL)
        gwfree(gwdata, p->YpX);
    p->YpX = NULL;
    if (p->YmX != NULL)
        gwfree(gwdata, p->YmX);
    p->YmX = NULL;

    gwsetnormroutine(gwdata, 0, 0, 1);
    gwsetmulbyconst(gwdata, 2);
    gwcopy(gwdata, p->X, p->T);
    gwmul_carefully(gwdata, p->Y, p->T); // E = 2X1Y1
    if (p->Z != NULL)
        gwsquare_carefully(gwdata, p->Z); // C = 2Z1^2
    else
    {
        p->Z = gwalloc(gwdata);
        dbltogw(gwdata, 2, p->Z);
    }
    gwsetnormroutine(gwdata, 0, 0, 0);
    gwsquare_carefully(gwdata, p->X); // A = X1^2
    gwsquare_carefully(gwdata, p->Y); // B = Y1^2
    force_normalize(p->Y);
    gwaddsub(gwdata, p->Y, p->X); // G = B + A, H = B - A
    force_normalize(p->Z);
    gwsub(gwdata, p->Y, p->Z); // F = C - G
    gwmul_carefully(gwdata, p->Z, p->T); // X3 = E * F
    gwmul_carefully(gwdata, p->Y, p->Z); // Z3 = F * G
    gwmul_carefully(gwdata, p->X, p->Y); // Y3 = G * H
    gwswap(p->X, p->T);
    gwfree(gwdata, p->T);
    p->T = NULL;
    costAdd(7);
}

void ed_add(ed_point padd, ed_point p, int options)
{
    if (p->T == NULL)
        p->T = gwalloc(gwdata);
    if (p->YpX != NULL)
        gwfree(gwdata, p->YpX);
    p->YpX = NULL;
    if (p->YmX != NULL)
        gwfree(gwdata, p->YmX);
    p->YmX = NULL;

    if (p->Z != NULL)
    {
        gwmul3(gwdata, padd->T, p->Z, p->Z, !must_norm ? GWMUL_STARTNEXTFFT : 0); // C = Z1 * T2
        costAdd(1);
    }
    else
    {
        p->Z = gwalloc(gwdata);
        if (must_norm)
            gwunfft(gwdata, padd->T, p->Z);
        else
            gwcopy(gwdata, padd->T, p->Z);
    }
    if (padd->Z != NULL)
    {
        gwmul3(gwdata, padd->Z, p->T, p->T, !must_norm ? GWMUL_STARTNEXTFFT : 0); // D = T1 * Z2
        costAdd(1);
    }
    if (must_norm)
        force_normalize(p->T);
    gwaddsub(gwdata, p->T, p->Z); // E = D + C, H = D - C
    gwnum tmp = gwalloc(gwdata);
    gwswap(tmp, p->X);
    gwsub3(gwdata, tmp, p->Y, p->X); // X1 - Y1
    if (padd->YpX == NULL)
    {
        padd->YpX = gwalloc(gwdata);
        gwcopy(gwdata, padd->X, padd->YpX);
        force_normalize(padd->YpX);
        gwadd(gwdata, padd->Y, padd->YpX);
    }
    gwmul3(gwdata, padd->YpX, p->X, p->X, 0); // (X1 - Y1) * (X2 + Y2)
    gwmul3(gwdata, padd->X, tmp, tmp, 0); // A = X1 * X2
    gwmul3(gwdata, padd->Y, p->Y, p->Y, 0); // B = Y1 * Y2
    if (must_norm)
        force_normalize(p->Y);
    gwaddsub(gwdata, p->Y, tmp); // G = B + A, B - A
    force_normalize(p->X);
    gwadd(gwdata, tmp, p->X); // F = (X1 - Y1) * (X2 + Y2) + B - A
    gwmul3(gwdata, p->Z, p->T, tmp, GWMUL_FFT_S1 | (!must_norm ? options : 0)); // T3 = E * H
    gwmul3(gwdata, p->X, p->T, p->T, options); // X3 = E * F
    gwmul3(gwdata, p->Y, p->Z, p->Z, options); // Y3 = G * H
    gwmul3(gwdata, p->X, p->Y, p->Y, options); // Z3 = F * G
    gwswap(p->Y, p->Z);
    gwswap(p->X, p->T);
    gwswap(p->T, tmp);
    costAdd(7);
    gwfree(gwdata, tmp);
}

void ed_mul2_add(ed_point padd, ed_point p, int options)
{
    if (p->T == NULL)
        p->T = gwalloc(gwdata);
    if (p->YpX != NULL)
        gwfree(gwdata, p->YpX);
    p->YpX = NULL;
    if (p->YmX != NULL)
        gwfree(gwdata, p->YmX);
    p->YmX = NULL;

    gwsetmulbyconst(gwdata, 2);
    gwmul3(gwdata, p->X, p->Y, p->T, GWMUL_FFT_S1 | GWMUL_MULBYCONST | GWMUL_STARTNEXTFFT); // E = 2X1Y1
    if (p->Z != NULL)
    {
        gwsquare2(gwdata, p->Z, p->Z, GWMUL_MULBYCONST | GWMUL_STARTNEXTFFT); // C = 2Z1^2
        costAdd(1);
    }
    else
    {
        p->Z = gwalloc(gwdata);
        dbltogw(gwdata, 2, p->Z);
    }
    gwsquare(gwdata, p->X); // A = X1^2
    gwsquare(gwdata, p->Y); // B = Y1^2
    force_normalize(p->Y);
    gwaddsub(gwdata, p->Y, p->X); // G = B + A, H = B - A
    gwsub(gwdata, p->Y, p->Z); // F = C - G
    gwnum tmp = gwalloc(gwdata);
    gwmul3(gwdata, p->Z, p->T, tmp, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT); // X3 = E * F
    gwmul3(gwdata, p->Y, p->Z, p->Z, GWMUL_STARTNEXTFFT); // Z3 = F * G
    if (options & ED_SIGNED_SUB)
        gwsetmulbyconst(gwdata, -1);
    gwmul3(gwdata, padd->T, p->Z, p->Z, ((options & ED_SIGNED_SUB) ? GWMUL_MULBYCONST : 0) | ((options & ED_MUL_FULL) ? GWMUL_STARTNEXTFFT1 : GWMUL_STARTNEXTFFT)); // C = Z1 * T2
    gwmul3(gwdata, p->X, p->T, p->T, (options & ED_MUL_FULL) && padd->Z == NULL ? GWMUL_STARTNEXTFFT1 : GWMUL_STARTNEXTFFT); // T3 = E * H
    if (padd->Z != NULL)
    {
        gwmul3(gwdata, padd->Z, p->T, p->T, (options & ED_MUL_FULL) ? GWMUL_STARTNEXTFFT1 : GWMUL_STARTNEXTFFT); // D = T1 * Z2
        costAdd(1);
    }
    if ((options & ED_MUL_FULL) && must_norm)
        force_normalize(p->T);
    gwaddsub(gwdata, p->T, p->Z); // E = D + C, H = D - C
    gwmul3(gwdata, p->X, p->Y, p->Y, GWMUL_STARTNEXTFFT); // Y3 = G * H
    gwfft(gwdata, tmp, tmp);
    gwfft(gwdata, p->Y, p->Y);
    gwsub3(gwdata, tmp, p->Y, p->X); // X1 - Y1
    if (!(options & (ED_SIGNED_ADD | ED_SIGNED_SUB)) && padd->YpX == NULL)
    {
        padd->YpX = gwalloc(gwdata);
        gwcopy(gwdata, padd->X, padd->YpX);
        force_normalize(padd->YpX);
        gwadd(gwdata, padd->Y, padd->YpX);
    }
    if ((options & (ED_SIGNED_ADD | ED_SIGNED_SUB)) && padd->YpX == NULL)
    {
        padd->YpX = gwalloc(gwdata);
        if (padd->YmX == NULL)
            padd->YmX = gwalloc(gwdata);
        gwcopy(gwdata, padd->Y, padd->YpX);
        gwcopy(gwdata, padd->X, padd->YmX);
        force_normalize(padd->YpX);
        gwaddsub(gwdata, padd->YpX, padd->YmX);
    }
    if (options & ED_SIGNED_SUB)
        gwmul3(gwdata, padd->YmX, p->X, p->X, 0); // (X1 - Y1) * (-X2 + Y2)
    else
        gwmul3(gwdata, padd->YpX, p->X, p->X, 0); // (X1 - Y1) * (X2 + Y2)
    gwmul3(gwdata, padd->X, tmp, tmp, (options & ED_SIGNED_SUB) ? GWMUL_MULBYCONST : 0); // A = X1 * X2
    gwmul3(gwdata, padd->Y, p->Y, p->Y, 0); // B = Y1 * Y2
    if (must_norm)
        force_normalize(p->Y);
    gwaddsub(gwdata, p->Y, tmp); // G = B + A, B - A
    force_normalize(p->X);
    gwadd(gwdata, tmp, p->X); // F = (X1 - Y1) * (X2 + Y2) + B - A
    if (options & ED_MUL_FULL)
        gwmul3(gwdata, p->Z, p->T, tmp, GWMUL_FFT_S1 | (options & (~ED_MUL_FULL))); // T3 = E * H
    gwmul3(gwdata, p->X, p->T, p->T, options & (~ED_MUL_FULL)); // X3 = E * F
    gwmul3(gwdata, p->Y, p->Z, p->Z, options & (~ED_MUL_FULL)); // Y3 = G * H
    gwmul3(gwdata, p->X, p->Y, p->Y, options & (~ED_MUL_FULL)); // Z3 = F * G
    gwswap(p->Y, p->Z);
    gwswap(p->X, p->T);
    if (options & ED_MUL_FULL)
    {
        gwswap(p->T, tmp);
        costAdd(1);
    }
    else
    {
        gwfree(gwdata, p->T);
        p->T = NULL;
    }
    costAdd(14);
    gwfree(gwdata, tmp);
}

void ed_mul2_addsmall(int xa, int xb, int ya, int yb, ed_point p, int options)
{
    GWASSERT(xa*ya <= gwdata->maxmulbyconst);
    GWASSERT(xb*yb <= gwdata->maxmulbyconst);

    gwsetnormroutine(gwdata, 0, 0, 1);
    gwsetmulbyconst(gwdata, 2);
    gwmul3(gwdata, p->X, p->Y, p->T, GWMUL_GLOBALMULBYCONST | GWMUL_STARTNEXTFFT); // E = 2X1Y1
    gwsetmulbyconst(gwdata, -2);
    gwsquare2(gwdata, p->Z, p->Z, GWMUL_GLOBALMULBYCONST | GWMUL_STARTNEXTFFT); // -C = -2Z1^2
    gwsetnormroutine(gwdata, 0, 0, 0);
    gwsquare(gwdata, p->X); // A = X1^2
    gwsquare(gwdata, p->Y); // B = Y1^2
    force_normalize(p->X);
    gwaddsub(gwdata, p->X, p->Y); // G = A + B, H = A - B
    gwadd(gwdata, p->X, p->Z); // F = G - C
    gwnum tmp = gwalloc(gwdata);
    gwmul3(gwdata, p->Z, p->T, tmp, 0); // X3 = E * F
    gwsetnormroutine(gwdata, 0, 0, 1);
    gwsetmulbyconst(gwdata, xa*ya);
    gwmul3(gwdata, p->X, p->Z, p->Z, GWMUL_GLOBALMULBYCONST | GWMUL_STARTNEXTFFT); // Z3 = F * G, C = Z1 * T2
    gwsetmulbyconst(gwdata, xb*yb);
    gwmul3(gwdata, p->Y, p->T, p->T, GWMUL_GLOBALMULBYCONST | GWMUL_STARTNEXTFFT); // T3 = E * H, D = T1 * Z2
    gwsetnormroutine(gwdata, 0, 0, 0);
    gwaddsub(gwdata, p->T, p->Z); // E = D + C, H = D - C
    gwmul3(gwdata, p->X, p->Y, p->Y, 0); // Y3 = G * H
    gwsub3quick(gwdata, tmp, p->Y, p->X); // X1 - Y1
    gwsmallmul(gwdata, xa*yb + ya*xb, p->X); // (X1 - Y1) * (X2 + Y2)
    gwsmallmul(gwdata, xa*yb, tmp); // A = X1 * X2
    gwsmallmul(gwdata, ya*xb, p->Y); // B = Y1 * Y2
    if (must_norm)
        force_normalize(p->Y);
    gwaddsub(gwdata, p->Y, tmp); // G = B + A, B - A
    force_normalize(p->X);
    gwadd(gwdata, tmp, p->X); // F = (X1 - Y1) * (X2 + Y2) + B - A
    gwmul3(gwdata, p->X, p->T, p->T, options); // X3 = E * F
    gwmul3(gwdata, p->Y, p->Z, p->Z, options); // Y3 = G * H
    gwmul3(gwdata, p->X, p->Y, p->Y, options); // Z3 = F * G
    gwswap(p->Y, p->Z);
    gwswap(p->X, p->T);
    costAdd(11);
    gwfree(gwdata, tmp);
}

giant ed_normalize_pool(ed_point* pool, int count)
{
    int i;
    int first = -1;
    int last;
    for (i = 0; i < count; i++)
    {
        if (pool[i] == NULL)
            continue;
        if (first < 0)
            first = i;
        last = i;
        if (pool[i]->T == NULL)
            pool[i]->T = gwalloc(gwdata);
        if (pool[i]->YpX != NULL)
            gwfree(gwdata, pool[i]->YpX);
        pool[i]->YpX = NULL;
        if (pool[i]->ZpY != NULL)
            gwfree(gwdata, pool[i]->ZpY);
        pool[i]->ZpY = NULL;
        if (pool[i]->ZmY != NULL)
            gwfree(gwdata, pool[i]->ZmY);
        pool[i]->ZmY = NULL;
    }
    if (first < 0)
        return NULL;
    gwswap(pool[first]->T, pool[first]->Z);
    int prev = first;
    for (i = first + 1; i <= last; i++)
        if (pool[i] != NULL)
        {
            gwmul3(gwdata, pool[prev]->T, pool[i]->Z, pool[i]->T, i != last ? GWMUL_STARTNEXTFFT : 0);
            prev = i;
            costAdd(1);
        }
    giant tmp = getg();
    gwtogiant(gwdata, pool[last]->T, tmp);
    invg(N, tmp);
    if (tmp->sign <= 0)
    {
        tmp->sign = -tmp->sign;
        return tmp;
    }
    gianttogw(gwdata, tmp, pool[last]->T);
    for (i = last; i >= first; i = prev)
    {
        if (i > first)
        {
            for (prev = i - 1; pool[prev] == NULL; prev--);
            gwmul3(gwdata, pool[i]->T, pool[prev]->T, pool[prev]->T, GWMUL_STARTNEXTFFT);
            gwswap(pool[i]->T, pool[prev]->T);
            gwmul3(gwdata, pool[i]->Z, pool[prev]->T, pool[prev]->T, GWMUL_STARTNEXTFFT);
            costAdd(2);
        }
        else
            prev = -1;
        if (pool[i]->X != NULL)
            gwmul3(gwdata, pool[i]->T, pool[i]->X, pool[i]->X, 0);
        if (pool[i]->Y != NULL)
            gwmul3(gwdata, pool[i]->T, pool[i]->Y, pool[i]->Y, 0);
        gwswap(pool[i]->YpX, pool[i]->Z);
        pool[i]->Z = NULL;
        if (pool[i]->X != NULL && pool[i]->Y != NULL)
        {
            if (pool[i]->YmX == NULL)
                pool[i]->YmX = gwalloc(gwdata);
            gwcopy(gwdata, pool[i]->Y, pool[i]->YpX);
            gwcopy(gwdata, pool[i]->X, pool[i]->YmX);
            force_normalize(pool[i]->YpX);
            gwaddsub(gwdata, pool[i]->YpX, pool[i]->YmX);
            gwmul3(gwdata, pool[i]->X, pool[i]->Y, pool[i]->T, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT);
            costAdd(3);
        }
        else
        {
            gwfree(gwdata, pool[i]->T);
            pool[i]->T = NULL;
            gwfree(gwdata, pool[i]->YpX);
            pool[i]->YpX = NULL;
            if (pool[i]->YmX != NULL)
                gwfree(gwdata, pool[i]->YmX);
            pool[i]->YmX = NULL;
            costAdd(1);
        }
    }
    freeg();
    return NULL;
}

void ed_y_mul2(ed_point p)
{
    if (p->X != NULL)
        gwfree(gwdata, p->X);
    p->X = NULL;
    if (p->T != NULL)
        gwfree(gwdata, p->T);
    p->T = NULL;
    if (p->YpX != NULL)
        gwfree(gwdata, p->YpX);
    p->YpX = NULL;
    if (p->YmX != NULL)
        gwfree(gwdata, p->YmX);
    p->YmX = NULL;
    if (p->ZpY == NULL)
        p->ZpY = gwalloc(gwdata);
    if (p->ZmY == NULL)
        p->ZmY = gwalloc(gwdata);

    gwsquare2(gwdata, p->Y, p->Y, GWMUL_STARTNEXTFFT1);
    gwsquare2(gwdata, p->Z, p->Z, GWMUL_STARTNEXTFFT1);
    gwcopy(gwdata, p->Z, p->ZmY);
    if (must_norm)
        force_normalize(p->ZmY);
    gwsub(gwdata, p->Y, p->ZmY);
    gwmul3(gwdata, EdD, p->Y, p->ZpY, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT);
    gwsubmul4(gwdata, p->Z, p->ZpY, p->ZmY, p->ZmY, GWMUL_STARTNEXTFFT1);
    gwsubmul4(gwdata, p->Y, p->ZpY, p->Z, p->Z, GWMUL_STARTNEXTFFT1);
    gwcopy(gwdata, p->ZmY, p->Y);
    gwcopy(gwdata, p->Z, p->ZpY);
    if (must_norm)
        force_normalize(p->Z);
    gwaddsub(gwdata, p->Z, p->Y);
}

void ed_y_inc(ed_point pinc, ed_point pprev, ed_point pcur, ed_point pnext)
{
    if (pprev->Y == NULL)
    {
        ed_copy(pcur, pnext);
        ed_y_mul2(pnext);
        return;
    }
    if (pnext->X != NULL)
        gwfree(gwdata, pnext->X);
    pnext->X = NULL;
    if (pnext->T != NULL)
        gwfree(gwdata, pnext->T);
    pnext->T = NULL;
    if (pnext->YpX != NULL)
        gwfree(gwdata, pnext->YpX);
    pnext->YpX = NULL;
    if (pnext->YmX != NULL)
        gwfree(gwdata, pnext->YmX);
    pnext->YmX = NULL;
    if (pnext->Y == NULL)
        pnext->Y = gwalloc(gwdata);
    if (pnext->Z == NULL)
        pnext->Z = gwalloc(gwdata);
    if (pnext->ZpY == NULL)
        pnext->ZpY = gwalloc(gwdata);
    if (pnext->ZmY == NULL)
        pnext->ZmY = gwalloc(gwdata);
    if (pprev->ZpY == NULL)
    {
        pprev->ZpY = gwalloc(gwdata);
        pprev->ZmY = gwalloc(gwdata);
        if (pprev->Z != NULL)
            gwcopy(gwdata, pprev->Z, pprev->ZpY);
        else
            dbltogw(gwdata, 1, pprev->ZpY);
        gwcopy(gwdata, pprev->Y, pprev->ZmY);
        if (must_norm)
            force_normalize(pprev->ZpY);
        gwaddsub(gwdata, pprev->ZpY, pprev->ZmY);
    }

    if (pinc->Z != NULL)
        gwmul3(gwdata, pcur->Y, pinc->Z, pnext->Z, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT1);
    else
        gwcopy(gwdata, pcur->Y, pnext->Z);
    if (pcur->Z != NULL)
        gwmul3(gwdata, pcur->Z, pinc->Y, pnext->Y, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT1);
    else
        gwcopy(gwdata, pinc->Y, pnext->Y);
    if (must_norm || pinc->Z == NULL || pcur->Z == NULL)
        force_normalize(pnext->Z);
    gwaddsub(gwdata, pnext->Z, pnext->Y);
    gwsquare2(gwdata, pnext->Z, pnext->Z, GWMUL_STARTNEXTFFT);
    gwsquare2(gwdata, pnext->Y, pnext->Y, GWMUL_STARTNEXTFFT);
    gwmul3(gwdata, pprev->ZmY, pnext->Z, pnext->Z, GWMUL_STARTNEXTFFT1);
    gwmul3(gwdata, pprev->ZpY, pnext->Y, pnext->Y, GWMUL_STARTNEXTFFT1);
    gwcopy(gwdata, pnext->Z, pnext->ZpY);
    gwcopy(gwdata, pnext->Y, pnext->ZmY);
    if (must_norm)
        force_normalize(pnext->Z);
    gwaddsub(gwdata, pnext->Z, pnext->Y);
}

void ed_mul_giant_add(giant g, ed_point p, ed_point pres, ed_point pres1)
{
    int K;
    //for (K = 1; 5 + 4*(1 << K) <= maxSize && (15 << (K - 1)) + bitlen(g)*(7 + 8/(K + 1.0)) > (15 << K) + bitlen(g)*(7 + 8/(K + 2.0)); K++);
    for (K = 1; 5 + 4*(1 << K) <= maxSize && (9 << (K - 1)) + bitlen(g)*(7 + 9/(K + 1.0)) > (9 << K) + bitlen(g)*(7 + 9/(K + 2.0)); K++);

    int i, j;
    int len;
    giant tmp = NULL;
    
    ed_point* u = malloc(sizeof(ed_point) << (K - 1));
    for (i = 0; i < (1 << (K - 1)); i++)
        u[i] = ed_alloc();

    // Dictionary
    ed_copy(p, u[0]);
    if (K > 1)
    {
        ed_copy(p, pres);
        ed_mul2(pres, ED_MUL_FULL);
        //tmp = ed_normalize(P);
        for (i = 1; i < (1 << (K - 1)) && tmp == NULL; i++)
        {
            ed_copy(u[i - 1], u[i]);
            ed_add(pres, u[i], GWMUL_STARTNEXTFFT);
        }
    }
    //if (tmp == NULL)
    //    tmp = ed_normalize_pool(u, 1 << (K - 1));

    if (tmp == NULL)
    {
        len = bitlen(g) - 1;

        // Sliding window
        int init = 0;
        for (i = len; i >= 0; i--)
        {
            if (bitval(g, i) == 0)
                ed_mul2(pres, i > 0 ? GWMUL_STARTNEXTFFT : ED_MUL_FULL);
            else
            {
                j = i - K + 1;
                if (j < 0)
                    j = 0;
                for (; bitval(g, j) == 0; j++);
                int ui = 0;
                while (i > j)
                {
                    if (init)
                        ed_mul2(pres, GWMUL_STARTNEXTFFT);
                    ui <<= 1;
                    ui += bitval(g, i) ? 1 : 0;
                    i--;
                }
                if (init)
                    ed_mul2_add(u[ui], pres, i > 0 ? GWMUL_STARTNEXTFFT : ED_MUL_FULL);
                else
                {
                    ed_copy(u[ui], pres);
                    init = 1;
                }
            }
        }

        if (pres1 != NULL)
        {
            ed_copy(pres, pres1);
            ed_add(p, pres1, 0);
        }
    }

    for (i = 0; i < (1 << (K - 1)); i++)
        ed_free(u[i]);
    free(u);
}

void ed_mul_int_add(int x, ed_point p, ed_point pres, ed_point pres1)
{
    giant g = allocgiant(1);
    itog(x, g);
    ed_mul_giant_add(g, p, pres, pres1);
    free(g);
}

giant ed_mul_nafw(int W, short *nafw, int len, ed_point P)
{
    int i, j;
    giant factor = NULL;

    ed_point* u = malloc(sizeof(ed_point) << (W - 2));
    for (i = 0; i < (1 << (W - 2)); i++)
        u[i] = ed_alloc();

    // Dictionary
    ed_copy(P, u[0]);
    if (W > 2)
    {
        ed_mul2(P, ED_MUL_FULL);
        //factor = ed_normalize(P);
        for (i = 1; i < (1 << (W - 2)) && factor == NULL; i++)
        {
            ed_copy(u[i - 1], u[i]);
            ed_add(P, u[i], GWMUL_STARTNEXTFFT);
        }
    }
    if (factor == NULL)
        factor = ed_normalize_pool(u, 1 << (W - 2));

    if (factor == NULL)
    {
        // Signed window
        ed_copy(u[nafw[len - 1]/2], P);
        for (i = len - 2; i >= 0; i--)
        {
            if (nafw[i] != 0)
                for (j = 1; j < W; j++)
                    ed_mul2(P, GWMUL_STARTNEXTFFT);
            if (nafw[i] > 0)
                ed_mul2_add(u[nafw[i]/2], P, (i > 0 ? GWMUL_STARTNEXTFFT : ED_MUL_FULL) | ED_SIGNED_ADD);
            else if (nafw[i] < 0)
                ed_mul2_add(u[-nafw[i]/2], P, (i > 0 ? GWMUL_STARTNEXTFFT : ED_MUL_FULL) | ED_SIGNED_SUB);
            else
                ed_mul2(P, i > 0 ? GWMUL_STARTNEXTFFT : ED_MUL_FULL);
        }
    }

    for (i = 0; i < (1 << (W - 2)); i++)
        ed_free(u[i]);
    free(u);

    return factor;
}

void ed_y_shiftleft(int x, ed_point p)
{
    for (; x > 0; x--, costAdd(1))
        ed_y_mul2(p);
}

void ed_y_optimize(ed_point s, ed_point d)
{
    if (s != d)
    {
        ed_zero(d);
        d->Y = gwalloc(gwdata);
        gwcopy(gwdata, s->Y, d->Y);
        d->Z = gwalloc(gwdata);
        gwcopy(gwdata, s->Z, d->Z);
    }
    //ed_normalize(d);
}

void get_j_invariant(giant res)
{
    gwnum M = gwalloc(gwdata);
    gwnum R = gwalloc(gwdata);

    gwcopy(gwdata, EdD, M);
    gwsmalladd(gwdata, 14, M);
    gwmul_carefully(gwdata, EdD, M);
    dbltogw(gwdata, 1, R);
    force_normalize(M);
    gwadd(gwdata, R, M);
    gwcopy(gwdata, M, R);
    gwsquare_carefully(gwdata, M);
    gwmul_carefully(gwdata, R, M);
    gwsmallmul(gwdata, 16, M);
    dbltogw(gwdata, 1, R);
    force_normalize(R);
    gwsub(gwdata, EdD, R);
    gwsquare_carefully(gwdata, R);
    gwsquare_carefully(gwdata, R);
    gwmul_carefully(gwdata, EdD, R);
    gwtogiant(gwdata, R, res);
    invg(N, res);
    if (res->sign > 0)
    {
        gianttogw(gwdata, res, R);
        gwmul_carefully(gwdata, R, M);
        gwtogiant(gwdata, M, res);
    }

    gwfree(gwdata, R);
    gwfree(gwdata, M);
}

int gen_ed_curve(int seed, ed_point P)
{
    gwnum S = gwalloc(gwdata);
    gwnum T = gwalloc(gwdata);
    gwnum X = gwalloc(gwdata);
    gwnum Y = gwalloc(gwdata);
    gwnum M = gwalloc(gwdata);
    gwnum R = gwalloc(gwdata);
    gwnum A = gwalloc(gwdata);
    gwnum B = gwalloc(gwdata);
    giant tmp = getg();
    giant r = getg();
    giant x8 = getg();
    giant sqrdx8 = getg();
    giant x = getg();
    giant y = getg();

    int i, len;
    int retval = TRUE;
    ultog(seed, tmp);

    // T^2 = S^3 - 8S - 32
    // (S, T) = id * (12, 40)
    dbltogw(gwdata, 12, S);
    dbltogw(gwdata, 40, T);
    len = bitlen(tmp) - 1;
    for (i = 1; i <= len; i++)
    {
        gwcopy(gwdata, S, M);
        gwsquare_carefully(gwdata, M);
        gwsmallmul(gwdata, 3, M);
        gwsmalladd(gwdata, -8, M);
        gwcopy(gwdata, T, R);
        gwsmallmul(gwdata, 2, R);
        gwtogiant(gwdata, R, r);
        invg(N, r);
        if (r->sign < 0)
            goto error;
        gianttogw(gwdata, r, R);
        gwmul_carefully(gwdata, R, M);
        gwcopy(gwdata, M, X);
        gwsquare_carefully(gwdata, X);
        gwsub(gwdata, S, X);
        force_normalize(X);
        gwsub(gwdata, S, X);
        gwcopy(gwdata, S, Y);
        gwsub(gwdata, X, Y);
        gwmul_carefully(gwdata, M, Y);
        force_normalize(Y);
        gwsub(gwdata, T, Y);
        gwswap(S, X);
        gwswap(T, Y);

        if (bitval(tmp, len - i))
        {
            dbltogw(gwdata, 40, M);
            gwsub(gwdata, T, M);
            dbltogw(gwdata, 12, R);
            gwsub(gwdata, S, R);
            gwtogiant(gwdata, R, r);
            invg(N, r);
            if (r->sign < 0)
                goto error;
            gianttogw(gwdata, r, R);
            gwmul_carefully(gwdata, R, M);
            gwcopy(gwdata, M, X);
            gwsquare_carefully(gwdata, X);
            gwsmalladd(gwdata, -12, X);
            force_normalize(X);
            gwsub(gwdata, S, X);
            gwcopy(gwdata, S, Y);
            gwsub(gwdata, X, Y);
            gwmul_carefully(gwdata, M, Y);
            force_normalize(Y);
            gwsub(gwdata, T, Y);
            gwswap(S, X);
            gwswap(T, Y);
        }
    }

    // Asserting T^2 = S^3 - 8S - 32
    /*giant s = popg(&gwdata->gdata, ((int)gwdata->bit_length >> 4) + 10);
    gwtogiant(gwdata, S, s);
    gtog(s, r);
    squareg(s);
    modg(N, s);
    sladdg(-8, s);
    mulg(r, s);
    sladdg(-32, s);
    modg(N, s);
    giant t = popg(&gwdata->gdata, ((int)gwdata->bit_length >> 4) + 10);
    gwtogiant(gwdata, T, t);
    squareg(t);
    modg(N, t);
    subg(t, s);
    GWASSERT(s->sign == 0);
    freeg();
    freeg();*/

    // A = 1/((T + 25)/(S - 9) + 1)
    gwcopy(gwdata, S, R);
    gwsmalladd(gwdata, -9, R);
    gwtogiant(gwdata, R, r);
    invg(N, r);
    if (r->sign < 0)
        goto error;
    gianttogw(gwdata, r, R);
    gwcopy(gwdata, T, A);
    gwsmalladd(gwdata, 25, A);
    gwmul_carefully(gwdata, R, A);
    gwsmalladd(gwdata, 1, A);
    gwtogiant(gwdata, A, r);
    invg(N, r);
    if (r->sign < 0)
        goto error;
    gianttogw(gwdata, r, A);

    // SqrD = X/Y
    // X = (8A^2 - 1) * (8A^2 + 8A + 1)
    // Y = (8A^2 + 4A + 1)^2
    // B = A * 2(4A + 1) / (8A^2 - 1)
    gwcopy(gwdata, A, X);
    gwsmallmul(gwdata, 8, X);
    gwcopy(gwdata, X, B);
    gwsmalladd(gwdata, 1, X);
    gwsmalladd(gwdata, 2, B);
    gwmul_carefully(gwdata, A, B);
    gwcopy(gwdata, A, Y);
    gwsmallmul(gwdata, 4, Y);
    gwsmalladd(gwdata, 1, Y);
    gwcopy(gwdata, A, R);
    gwsquare_carefully(gwdata, R);
    gwsmallmul(gwdata, 8, R);
    force_normalize(X);
    gwadd(gwdata, R, X);
    force_normalize(Y);
    gwadd(gwdata, R, Y);
    gwsquare_carefully(gwdata, Y);
    gwsmalladd(gwdata, -1, R);
    gwmul_carefully(gwdata, R, X);
    gwtogiant(gwdata, R, r);
    invg(N, r);
    if (r->sign < 0)
        goto error;
    gianttogw(gwdata, r, R);
    gwmul_carefully(gwdata, R, B);

    // D = SqrD^2
    gwtogiant(gwdata, Y, r);
    invg(N, r);
    if (r->sign < 0)
        goto error;
    gianttogw(gwdata, r, Y);
    gwmul_carefully(gwdata, Y, X);
    gwcopy(gwdata, X, Y);
    gwsquare_carefully(gwdata, X);
    EdD = gwalloc(gwdata);
    gwcopy(gwdata, X, EdD);

    // x8 = 2B - 1
    // 1/(SqrD * x8)
    gwcopy(gwdata, B, A);
    gwsmallmul(gwdata, 2, A);
    //force_normalize(A);
    gwsmalladd(gwdata, -1, A);
    gwtogiant(gwdata, A, x8);
    gwmul_carefully(gwdata, A, Y);
    gwtogiant(gwdata, Y, r);
    invg(N, r);
    if (r->sign < 0)
        goto error;
    gtog(r, sqrdx8);

    // X = x8 * (4B - 3) / (6B - 5)
    gwcopy(gwdata, B, R);
    gwsmallmul(gwdata, 6, R);
    gwsmalladd(gwdata, -5, R);
    gwtogiant(gwdata, R, r);
    invg(N, r);
    if (r->sign < 0)
        goto error;
    gianttogw(gwdata, r, R);
    gwcopy(gwdata, B, X);
    gwsmallmul(gwdata, 4, X);
    gwsmalladd(gwdata, -3, X);
    gwmul_carefully(gwdata, R, X);
    gwmul_carefully(gwdata, A, X);

    // Y = x8 * (T^2 + 50T - 2S^3 + 27S^2 - 104) / ((T + 3S - 2) * (T + S + 16))
    gwcopy(gwdata, S, M);
    gwsmallmul(gwdata, 3, M);
    gwsmalladd(gwdata, -2, M);
    force_normalize(M);
    gwadd(gwdata, T, M);
    gwcopy(gwdata, S, R);
    gwsmalladd(gwdata, 16, R);
    force_normalize(R);
    gwadd(gwdata, T, R);
    gwmul_carefully(gwdata, M, R);
    gwtogiant(gwdata, R, r);
    invg(N, r);
    if (r->sign < 0)
        goto error;
    gianttogw(gwdata, r, R);
    gwcopy(gwdata, S, M);
    gwsmallmul(gwdata, -2, M);
    gwsmalladd(gwdata, 27, M);
    gwmul_carefully(gwdata, S, M);
    gwmul_carefully(gwdata, S, M);
    gwsmalladd(gwdata, -104, M);
    gwcopy(gwdata, T, Y);
    gwsmalladd(gwdata, 50, Y);
    gwmul_carefully(gwdata, T, Y);
    gwadd(gwdata, M, Y);
    gwmul_carefully(gwdata, R, Y);
    gwmul_carefully(gwdata, A, Y);

    gwtogiant(gwdata, X, x);
    gwtogiant(gwdata, Y, y);

    // Checking torsion points
    gtog(N, r);
    sladdg(-1, r);
    if (isZero(x) && isone(y))
        retval = FALSE;
    if (isZero(x) && !gcompg(y, r))
        retval = FALSE;
    if (isone(x) && isZero(y))
        retval = FALSE;
    if (!gcompg(x, r) && isZero(y))
        retval = FALSE;
    gtog(N, r);
    subg(x8, r);
    if (!gcompg(x, x8) && !gcompg(y, x8))
        retval = FALSE;
    if (!gcompg(x, r) && !gcompg(y, x8))
        retval = FALSE;
    if (!gcompg(x, x8) && !gcompg(y, r))
        retval = FALSE;
    if (!gcompg(x, r) && !gcompg(y, r))
        retval = FALSE;
    gtog(N, r);
    subg(sqrdx8, r);
    if (!gcompg(x, sqrdx8) && !gcompg(y, sqrdx8))
        retval = FALSE;
    if (!gcompg(x, r) && !gcompg(y, sqrdx8))
        retval = FALSE;
    if (!gcompg(x, sqrdx8) && !gcompg(y, r))
        retval = FALSE;
    if (!gcompg(x, r) && !gcompg(y, r))
        retval = FALSE;

    // Asserting X^2 + Y^2 = 1 + d * X^2 * Y^2
    /*gwcopy(gwdata, X, A);
    gwsquare_carefully(gwdata, A);
    gwcopy(gwdata, Y, B);
    gwsquare_carefully(gwdata, B);
    gwcopy(gwdata, A, M);
    gwadd(gwdata, B, M);
    gwmul_carefully(gwdata, B, A);
    gwmul_carefully(gwdata, EdD, A);
    gwsmalladd(gwdata, 1, A);
    gwsub(gwdata, A, M);
    gwtogiant(gwdata, M, tmp);
    GWASSERT(tmp->sign == 0);*/

    if (retval)
    {
        P->X = X;
        X = NULL;
        P->Y = Y;
        Y = NULL;
        P->T = gwalloc(gwdata);
        gwcopy(gwdata, P->X, P->T);
        gwmul_carefully(gwdata, P->Y, P->T);
    }

    goto cleanup;
error:
    report_factor(r);
    retval = FALSE;
cleanup:
    gwfree(gwdata, S);
    gwfree(gwdata, T);
    gwfree(gwdata, R);
    gwfree(gwdata, M);
    gwfree(gwdata, A);
    gwfree(gwdata, B);
    if (X != NULL)
        gwfree(gwdata, X);
    if (Y != NULL)
        gwfree(gwdata, Y);
    freeg();
    freeg();
    freeg();
    freeg();
    freeg();
    freeg();

    return retval;
}