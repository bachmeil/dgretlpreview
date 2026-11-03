/* Port of lib/src/random.c: gretl's pseudo-random number generators. The engine is xoshiro256+ (seeded through
   splitmix64), both included here as the C includes them (rng/splitmix64.c and rng/xoshiro256plus.c, public
   domain, by Sebastiano Vigna and David Blackman); on top of it are the Ziggurat normal generator, uniform,
   gamma, chi-square, Student, F, binomial, Poisson, Weibull, exponential, logistic, GED, Laplace, beta,
   beta-binomial, truncated normal, discrete (alias method), Dirichlet, Halton sequences, inverse Wishart
   matrices and random hex strings. The header random.h has only prototypes.

   Differences from the C: the MPI initialisation (and its stand-in) and the 32-bit-x86 Ziggurat variant are
   dropped, and the Windows seeding branch is dropped. */

module gretl.random;


import core.stdc.errno : errno;
import core.stdc.float_ : DBL_MAX, DBL_MIN;
import core.stdc.limits : INT_MAX;
import core.stdc.math : cos, exp, expm1, fabs, floor, log, log1p, pow, sin, sqrt;
import core.stdc.stdio : fprintf, printf, sprintf, stderr;
import core.stdc.stdlib : calloc, free, getenv, malloc, realloc, strtoull;
import core.stdc.string : memcpy, memset;
import core.stdc.time : time;
version (Posix) import core.sys.posix.fcntl : O_RDONLY, open; else import gretl.winport : O_RDONLY, open;
version (Posix) import core.sys.posix.sys.types : ssize_t; else import gretl.winport : ssize_t;
version (Posix) import core.sys.posix.unistd : read; else import gretl.winport : read;

import std.algorithm : max, min;

import gretl.gretl_errorsh : GretlError;
import gretl.gretl_matrix;
import gretl.gretl_matrixh;
import gretl.libgretlh;
import gretl.missingh : NADBL, na;

import gretl.gretl_utils : ijton;
import gretl.pvalues : gammafun, normal_cdf, normal_cdf_inverse;
import gretl.strutils : gretl_strdup;

/* splitmix64 (rng/splitmix64.c): a fixed-increment version of Java 8's SplittableRandom generator. See
   http://dx.doi.org/10.1145/2714064.2660195 */

private __gshared ulong smx; /* The state can be seeded with any value. */

ulong splitmix64_next()
{
    ulong z = (smx += 0x9e3779b97f4a7c15UL);
        
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9UL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebUL;
    return z ^ (z >> 31);
}

/* xoshiro256+ 1.0 (rng/xoshiro256plus.c), written in 2018 by David Blackman and Sebastiano Vigna (public domain). */

/* This is xoshiro256+ 1.0, our best and fastest generator for floating-point
   numbers. We suggest to use its upper bits for floating-point
   generation, as it is slightly faster than xoshiro256++/xoshiro256**. It
   passes all tests we are aware of except for the lowest three bits,
   which might fail linearity tests (and just those), so if low linear
   complexity is not considered an issue (as it is usually the case) it
   can be used to generate 64-bit outputs, too.

   We suggest to use a sign test to extract a random Boolean value, and
   right shifts to extract subsets of bits.

   The state must be seeded so that it is not everywhere zero. If you have
   a 64-bit seed, we suggest to seed a splitmix64 generator and use its
   output to fill s.
*/

ulong xor_rotl(const ulong x, int k)
{
    return (x << k) | (x >> (64 - k));
}

private __gshared ulong[4] xor_state;

void set_xor_state(ulong u)
{
    xor_state[0] = smx = u;
    xor_state[1] = splitmix64_next();
    xor_state[2] = splitmix64_next();
    xor_state[3] = splitmix64_next();
}

/* A replacement for the engine, local to the thread (module-level variables are thread-local in D unless __gshared): while a source is set, every draw
   of this thread, by every generator of this file, takes its 64 random bits from it instead of from the xoshiro state, which is then left alone. This is
   not in the C. It lets another generator, for example one stream of a family of parallel generators, drive all the distributions of this file, one thread
   for each source, with no change to the global generator. */

alias RandSource = extern (D) ulong function(void *);

private RandSource source_fn;           // thread-local on purpose
private void * source_ctx;

void gretl_rand_set_source(RandSource fn, void * ctx)
{
    source_fn = fn;
    source_ctx = ctx;
}

void gretl_rand_get_source(RandSource * fn, void ** ctx)
{
    *fn = source_fn;
    *ctx = source_ctx;
}

void gretl_rand_clear_source()
{
    source_fn = null;
    source_ctx = null;
}

/* The tables of the Ziggurat generator are made on first use: this makes them (once, with no draw) so that threads that use the generator
   at the same time do not race to do it. */

void gretl_rand_prepare()
{
    if (initt) {
        create_ziggurat_tables();
    }
}

/* get the next pseudo-random ulong: from the source of this thread if it has one, else from the xoshiro engine */

ulong xor_i64()
{
    if (source_fn !is null) {
        return source_fn(source_ctx);
    }
    return xoshiro_next();
}

private ulong xoshiro_next()
{
    const ulong ret = xor_state[0] + xor_state[3];
    const ulong t = xor_state[1] << 17;

    xor_state[2] ^= xor_state[0];
    xor_state[3] ^= xor_state[1];
    xor_state[1] ^= xor_state[2];
    xor_state[0] ^= xor_state[3];

    xor_state[2] ^= t;
    xor_state[3] = xor_rotl(xor_state[3], 45);

    return ret;
}

/* This is the jump function for the generator. It is equivalent to
   2^128 calls to xor_i64(); it can be used to generate 2^128
   non-overlapping subsequences for parallel computations.
*/

void xor_jump()
{
    static immutable ulong[4] JUMP = [
        0x180ec6d33cfd0aba, 0xd5a61266f0c9392c,
        0xa9582618e03fc9aa, 0x39abdc4529b1661c
    ];
    ulong s0 = 0;
    ulong s1 = 0;
    ulong s2 = 0;
    ulong s3 = 0;
    int i, b;

    for (i = 0; i < JUMP.length; i++) {
        for (b = 0; b < 64; b++) {
            if (JUMP[i] & (1UL << b)) {
                s0 ^= xor_state[0];
                s1 ^= xor_state[1];
                s2 ^= xor_state[2];
                s3 ^= xor_state[3];
            }
            xoshiro_next();
        }
    }

    xor_state[0] = s0;
    xor_state[1] = s1;
    xor_state[2] = s2;
    xor_state[3] = s3;
}



/**
 * SECTION:random
 * @short_description: generate pseudo-random values
 * @title: PRNG
 * @include: libgretl.h
 *
 * Libgretl uses xoshiro256+ as its underlying engine for
 * uniform random values, but offers added value in the
 * form of generators for several distributions commonly
 * used in econometrics.
 *
 * Note that before using the libgretl PRNG you must call
 * either libgretl_init() or the specific initialization
 * function shown below, gretl_rand_init().
 */

private __gshared ulong xor_seed;

double double_from_uint64(ulong u)
{
    /* Set the exponent to 0x3FF (for 1.0) and the sign bit to 0;
       the remaining 52 bits are from the random ulong.
    */
    const ulong MASK = 0x3FF0000000000000UL;
    ulong bits = MASK | (u >> 12);
    double ret;

    /* create a double in the range [1.0, 2.0) */
    memcpy(&ret, &bits, double.sizeof);

    return ret - 1.0;
}

uint randi32()
{
    return cast(uint) (xor_i64() >> 32);
}

ulong get_auto_seed()
{
    ulong u;

    int fd = open("/dev/urandom", O_RDONLY);

    if (fd == -1) {
        u = time(null);
    } else {
        ssize_t sz = read(fd, &u, ulong.sizeof);

        if (sz != ulong.sizeof) {
            u = time(null);
        }
    }

    return u;
}

/**
 * gretl_rand_init:
 *
 * Initialize gretl's PRNG, using the system time as seed.
 */

void gretl_rand_init()
{
    char * fseed = getenv("GRETL_FORCE_SEED");
    int done = 0;

    if (fseed != null) {
        char * test = null;
        ulong u;

        errno = 0;
        u = strtoull(fseed, &test, 10);
        if (*test == '\0' && errno == 0 && u > 0) {
            xor_seed = u;
            done = 1;
        }
    }

    if (!done) {
        xor_seed = get_auto_seed();
    }

    set_xor_state(xor_seed);
}



/**
 * gretl_rand_get_seed:
 *
 * Returns: the value of the seed for gretl's PRNG.
 */

ulong gretl_rand_get_seed()
{
    return xor_seed;
}

/**
 * gretl_rand_set_seed:
 * @seed: the chosen seed value.
 *
 * Set a specific (and hence reproducible) seed for gretl's PRNG.
 * But if the value 0 is given for @seed, set the seed
 * automatically.
 */

void gretl_rand_set_seed(ulong seed)
{
    if (seed > 0) {
        xor_seed = cast(ulong) seed;
        set_xor_state(xor_seed);
    } else {
        xor_seed = get_auto_seed();
        set_xor_state(xor_seed);
    }
}

/**
 * gretl_rand_01:
 *
 * Returns: the next random double, equally distributed over
 * the range [0, 1).
 */

double gretl_rand_01()
{
    return double_from_uint64(xor_i64());
}


/* 53 bits for mantissa + 1 bit sign */

ulong randi54()
{
    const ulong u = xor_i64();
    const ulong mask = (1UL << 54) - 1;

    return u & mask;
}


/* generates a uniform random double on (0,1) with 53-bit resolution */

double randu53()
{
    const ulong u = xor_i64() >> 11;

    return cast(double) u * 0x1.0p-53;
}

/* Ziggurat normal generator: this Ziggurat code here is shamelessly
   filched from GNU Octave (with minor modifications), since it turned
   out to be faster than what we had (until 2016-10-08), with no loss
   in quality of the random normal variates (in fact, perhaps a gain).

   See randmtzig.c in the Octave code-base. It appears that David
   Bateman was the main author; anyway, thanks to whoever wrote the
   implementation!
*/

enum ZIGGURAT_TABLE_SIZE = 256;
enum ZIGGURAT_NOR_R = 3.6541528853610088;
enum ZIGGURAT_NOR_INV_R = 0.27366123732975828;
enum NOR_SECTION_AREA = 0.00492867323399;

enum NMANTISSA = 9007199254740992.0;  /* 53 bit mantissa */

private __gshared ulong[ZIGGURAT_TABLE_SIZE] ki = 0;
private __gshared double[ZIGGURAT_TABLE_SIZE] wi = 0;
private __gshared double[ZIGGURAT_TABLE_SIZE] fi = 0;

private __gshared int initt = 1;

void create_ziggurat_tables()
{
    int i;
    double x, x1;

    x1 = ZIGGURAT_NOR_R;
    wi[255] = x1 / NMANTISSA;
    fi[255] = exp(-0.5 * x1 * x1);

    /* Index zero is special for tail strip, where Marsaglia and Tsang
       defines this as:
       k_0 = 2^31 * r * f(r) / v, w_0 = 0.5^31 * v / f(r), f_0 = 1,
       where v is the area of each strip of the ziggurat.
    */
    ki[0] = cast(ulong) (x1 * fi[255] / NOR_SECTION_AREA * NMANTISSA);
    wi[0] = NOR_SECTION_AREA / fi[255] / NMANTISSA;
    fi[0] = 1.;

    for (i = 254; i > 0; i--) {
    /* New x is given by x = f^{-1}(v/x_{i+1} + f(x_{i+1})), thus
       need inverse operator of y = exp(-0.5*x*x) -> x = sqrt(-2*ln(y))
    */
    x = sqrt (-2. * log(NOR_SECTION_AREA / x1 + fi[i+1]));
    ki[i+1] = cast(ulong) (x / x1 * NMANTISSA);
    wi[i] = x / NMANTISSA;
    fi[i] = exp(-0.5 * x * x);
    x1 = x;
    }

    ki[1] = 0;

    initt = 0;
}

/**
 * gretl_one_snormal:
 *
 * Returns: a single drawing from the standard normal distribution.
 */

double gretl_one_snormal()
{
    if (initt) {
    create_ziggurat_tables();
    }

    while (1) {
    const ulong r = randi54();
    const long rabs = r >> 1;
    const int idx = cast(int) (rabs & 0xFF);
    const(double) x = ((r & 1) ? -rabs : rabs) * wi[idx];

    if (rabs < cast(long) (ki[idx])) {
        return x; /* 99.3% of the time we return here 1st try */
    } else if (idx == 0) {
        /* As stated in Marsaglia and Tsang:
           For the normal tail, the method of Marsaglia provides:
           generate x = -ln(U_1)/r, y = -ln(U_2), until y+y > x*x,
           then return r+x. Except that r+x is always in the positive
           tail. Anything random might be used to determine the
           sign, but as we already have r we might as well use it.
        */
        double xx, yy;

        do {
        xx = - ZIGGURAT_NOR_R * log(randu53());
        yy = - log(randu53());
            } while (yy+yy <= xx*xx);
        return (rabs & 0x100) ? -ZIGGURAT_NOR_R-xx : ZIGGURAT_NOR_R+xx;
        } else if ((fi[idx-1] - fi[idx]) * randu53() + fi[idx] < exp(-0.5*x*x)) {
        return x;
    }
    }
}

/**
 * gretl_rand_normal:
 * @a: target array
 * @t1: start of the fill range
 * @t2: end of the fill range
 *
 * Fill the selected range of array @a with pseudo-random drawings
 * from the standard normal distribution, using the Mersenne Twister
 * for uniform input and the Ziggurat method for converting to the
 * normal distribution.
 */

void gretl_rand_normal(double * a, int t1, int t2)
{
    int t;

    for (t=t1; t<=t2; t++) {
    a[t] = gretl_one_snormal();
    }
}

/**
 * gretl_rand_normal_full:
 * @a: target array
 * @t1: start of the fill range
 * @t2: end of the fill range
 * @mean: mean of the distribution
 * @sd: standard deviation
 *
 * Fill the selected range of array @a with pseudo-random drawings
 * from the normal distribution with the given mean and standard
 * deviation, using the Mersenne Twister for uniform input and
 * the Ziggurat method for converting to the normal distribution.
 *
 * Returns: 0 on success, 1 on invalid input.
 */

int gretl_rand_normal_full(double * a, int t1, int t2, double mean, double sd)
{
    int t;

    if (na(mean) && na(sd)) {
    mean = 0.0;
    sd = 1.0;
    } else if (na(mean) || na(sd) || sd <= 0.0) {
    return GretlError.invarg;
    }

    gretl_rand_normal(a, t1, t2);

    if (mean != 0.0 || sd != 1.0) {
    for (t=t1; t<=t2; t++) {
        a[t] = mean + a[t] * sd;
    }
    }

    return 0;
}

uint rand_int_range(uint begin, uint end)
{
    uint dist = end - begin;
    uint rval = 0;

    if (dist > 0) {
    /* maxval is set to the predecessor of the greatest
       multiple of dist less than or equal to 2^32
    */
    uint maxval;

    if (dist <= 0x80000000u) { /* 2^31 */
        /* maxval = 2^32 - 1 - (2^32 % dist) */
        uint rem = (0x80000000u % dist) * 2;

        if (rem >= dist) rem -= dist;
        maxval = 0xffffffffu - rem;
    } else {
        maxval = dist - 1;
    }

        do {
            rval = randi32();
        } while (rval > maxval);

    rval %= dist;
    }

    return begin + rval;
}

/**
 * gretl_rand_uniform_minmax:
 * @a: target array.
 * @t1: start of the fill range.
 * @t2: end of the fill range.
 * @min: lower bound of range.
 * @max: upper bound of range.
 *
 * Fill the selected subset of array @a with pseudo-random drawings
 * from the uniform distribution on @min to @max, using the Mersenne
 * Twister.
 *
 * Returns: 0 on success, 1 on invalid input.
 */

int gretl_rand_uniform_minmax(double * a, int t1, int t2, double min, double max)
{
    int t;

    if (na(min) && na(max)) {
    min = 0.0;
    max = 1.0;
    } else if (na(min) || na(max) || max <= min) {
    return GretlError.invarg;
    }

    for (t=t1; t<=t2; t++) {
        a[t] = double_from_uint64(xor_i64()) * (max - min) + min;
    }

    return 0;
}

/**
 * gretl_rand_int_minmax:
 * @a: target array.
 * @n: length of array.
 * @min: lower closed bound of range.
 * @max: upper closed bound of range.
 *
 * Fill array @a of length @n with pseudo-random drawings
 * from the uniform distribution on [@min, @max].
 *
 * Returns: 0 on success, 1 on invalid input.
 */

int gretl_rand_int_minmax(int * a, int n, int min, int max)
{
    int i, err = 0;

    if (max < min) {
    err = GretlError.invarg;
    } else if (min == max) {
    for (i=0; i<n; i++) {
        a[i] = min;
    }
    } else {
    int offset = 0;

    if (min < 0) {
        offset = -min;
        max += offset;
        min += offset;
    }

    for (i=0; i<n; i++) {
        a[i] = rand_int_range(min, max + 1) - offset;
    }
    }

    return err;
}

int already_selected(double * a, int n, double val, int offset)
{
    int i;

    for (i=0; i<n; i++) {
    if (offset > 0) {
        if (a[i] == val - offset) {
        return 1;
        }
    } else if (a[i] == val) {
        return 1;
    }
    }

    return 0;
}

/**
 * gretl_rand_uniform_int_minmax:
 * @a: target array.
 * @t1: start of the fill range.
 * @t2: end of the fill range.
 * @min: lower closed bound of range.
 * @max: upper closed bound of range.
 * @opt: include OPT_O for sampling without replacement.
 *
 * Fill the selected subset of array @a with pseudo-random drawings
 * from the uniform distribution on [@min, @max], using the
 * Mersenne Twister.
 *
 * Returns: 0 on success, 1 on invalid input.
 */

int gretl_rand_uniform_int_minmax(double * a, int t1, int t2, int min, int max, Gretlopt opt)
{
    int t, err = 0;

    if (max < min) {
    err = GretlError.invarg;
    } else if (max == min) {
    for (t=t1; t<=t2; t++) {
        a[t] = min;
    }
    } else {
    int i = 0, offset = 0;
    double x;

    if (min < 0) {
        offset = -min;
        max += offset;
        min += offset;
    }

    for (t=t1; t<=t2; t++) {
        x = rand_int_range(min, max + 1);
        if (opt & Gretlopt.OPT_O) {
        while (already_selected(a, i, x, offset)) {
            x = rand_int_range(min, max + 1);
        }
        }
        a[t] = x - offset;
        i++;
    }
    }

    return err;
}

/**
 * gretl_rand_uniform:
 * @a: target array
 * @t1: start of the fill range
 * @t2: end of the fill range
 *
 * Fill the selected range of array @a with pseudo-random drawings
 * from the uniform distribution on [0-1), using the Mersenne
 * Twister.
 */

void gretl_rand_uniform(double * a, int t1, int t2)
{
    int t;

    for (t=t1; t<=t2; t++) {
        a[t] = double_from_uint64(xor_i64());
    }
}

double gretl_rand_uniform_one()
{
    return double_from_uint64(xor_i64());
}

double gretl_rand_gamma_one(double shape, double scale)
{
    double k = shape;
    double d, c, x, v, u, dv;

    if (shape <= 0 || scale <= 0) {
    return NADBL;
    }

    if (shape < 1) {
    k = shape + 1.0;
    }

    d = k - 1.0/3;
    c = 1.0 / sqrt(9*d);

    while (1) {
    x = gretl_one_snormal();
    v = pow(1 + c*x, 3);
    if (v > 0.0) {
        dv = d * v;
        u = gretl_rand_01();
        /* apply squeeze */
        if (u < 1 - 0.0331 * pow(x, 4) ||
        log(u) < 0.5*x*x + d*(1-v+log(v))) {
        break;
        }
    }
    }
    if (shape < 1) {
    u = gretl_rand_01();
    dv *= pow(u, 1/shape);
    }

    return dv * scale;
}

/* Marsaglia-Tsang, "A Simple Method for Generating Gamma Variables",
   ACM Transactions on Mathematical Software, Vol. 26, No. 3,
   September 2000, Pages 363­372.

   (1) Setup: d=a-1/3, c=1/sqrt(9*d).
   (2) Generate v=(1+c*x)^3 with x normal; repeat if v <= 0.
   (3) Generate uniform U.
   (4) If U < 1-0.0331*x^4 return d*v.
   (5) If log(U) < 0.5*x^2+d*(1-v+log(v)) return d*v.
   (6) Go to step 2.
*/

/**
 * gretl_rand_gamma:
 * @a: target array.
 * @t1: start of the fill range.
 * @t2: end of the fill range.
 * @shape: shape parameter.
 * @scale: scale parameter.
 *
 * Fill the selected range of array @a with pseudo-random drawings
 * from the specified gamma distribution.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_rand_gamma(double * a, int t1, int t2, double shape, double scale)
{
    double k = shape;
    double d, c, x, v, u, dv;
    int t;

    if (shape <= 0 || scale <= 0) {
    return GretlError.data;
    }

    if (shape < 1) {
    k = shape + 1.0;
    }

    d = k - 1.0/3;
    c = 1.0 / sqrt(9*d);

    for (t=t1; t<=t2; t++) {
    while (1) {
        x = gretl_one_snormal();
        v = pow(1 + c*x, 3);
        if (v > 0.0) {
        dv = d * v;
        u = gretl_rand_01();
        /* apply squeeze */
        if (u < 1 - 0.0331 * pow(x, 4) ||
            log(u) < 0.5*x*x + d*(1-v+log(v))) {
            break;
        }
        }
    }
    if (shape < 1) {
        u = gretl_rand_01();
        dv *= pow(u, 1/shape);
    }
    a[t] = dv * scale;
    }

    return 0;
}

/**
 * gretl_rand_chisq:
 * @a: target array.
 * @t1: start of the fill range.
 * @t2: end of the fill range.
 * @v: degrees of freedom.
 *
 * Fill the selected range of array @a with pseudo-random drawings
 * from the Chi-squared distribution with @v degrees of freedom,
 * using the gamma r.v. generator.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_rand_chisq(double * a, int t1, int t2, int v)
{
    return gretl_rand_gamma(a, t1, t2, 0.5*v, 2);
}

/**
 * gretl_rand_student:
 * @a: target array.
 * @t1: start of the fill range.
 * @t2: end of the fill range.
 * @v: degrees of freedom.
 *
 * Fill the selected range of array @a with pseudo-random drawings
 * from the Student t distribution with @v degrees of freedom,
 * using the Mersenne Twister for uniform input and the ziggurat
 * method for converting to the normal distribution.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_rand_student(double * a, int t1, int t2, double v)
{
    double * X2 = null;
    int T = t2 - t1 + 1;
    int t;

    if (v <= 0) {
    return GretlError.invarg;
    }

    X2 = cast(typeof(X2)) malloc(T * (*X2).sizeof);
    if (X2 == null) {
    return GretlError.alloc;
    }

    gretl_rand_normal(a, t1, t2);
    gretl_rand_gamma(X2, 0, T-1, 0.5*v, 2);

    for (t=0; t<T; t++) {
    a[t + t1] /= sqrt(X2[t] / v);
    }

    free(X2);

    return 0;
}

/**
 * gretl_rand_F:
 * @a: target array.
 * @t1: start of the fill range.
 * @t2: end of the fill range.
 * @v1: numerator degrees of freedom.
 * @v2: denominator degrees of freedom.
 *
 * Fill the selected range of array @a with pseudo-random drawings
 * from the F distribution with @v1 and @v2 degrees of freedom,
 * using the Mersenne Twister for uniform input and the Ziggurat
 * method for converting to the normal distribution.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_rand_F(double * a, int t1, int t2, int v1, int v2)
{
    double * b = null;
    int T = t2 - t1 + 1;
    int s, t;

    if (v1 < 1 || v2 < 1) {
    return GretlError.invarg;
    }

    b = cast(typeof(b)) malloc(T * (*b).sizeof);
    if (b == null) {
    return GretlError.alloc;
    }

    gretl_rand_chisq(a, t1, t2, v1);
    gretl_rand_chisq(b, 0, T-1, v2);

    for (t=0; t<T; t++) {
    s = t + t1;
    a[s] = (a[s]/v1) / (b[t]/v2);
    }

    free(b);

    return 0;
}

/**
 * gretl_rand_binomial:
 * @a: target array.
 * @t1: start of the fill range.
 * @t2: end of the fill range.
 * @n: number of trials.
 * @p: success probability per trial.
 *
 * Fill the selected range of array @a with pseudo-random drawings
 * from the binomial distribution with parameters @n and @p.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_rand_binomial(double * a, int t1, int t2, int n, double p)
{
    int t;

    if (n < 0 || p < 0 || p > 1) {
    return GretlError.invarg;
    }

    if (n == 0 || p == 0.0) {
    for (t=t1; t<=t2; t++) {
        a[t] = 0.0;
    }
    } else if (p == 1.0) {
    for (t=t1; t<=t2; t++) {
        a[t] = n;
    }
    } else {
    double * b = cast(double *) malloc(n * double.sizeof);
    int i;

    if (b == null) {
        return GretlError.alloc;
    }

    for (t=t1; t<=t2; t++) {
        a[t] = 0.0;
        gretl_rand_uniform(b, 0, n - 1);
        for (i=0; i<n; i++) {
        if (b[i] <= p) {
            a[t] += 1;
        }
        }
    }

    free(b);
    }

    return 0;
}

double gretl_rand_binomial_one(int n, double p, double * b)
{
    double ret;

    if (n < 0 || p < 0 || p > 1) {
    return NADBL;
    }

    if (n == 0 || p == 0.0) {
    ret = 0.0;
    } else if (p == 1.0) {
    ret = n;
    } else {
    int i;

    ret = 0.0;
    gretl_rand_uniform(b, 0, n - 1);
    for (i=0; i<n; i++) {
        if (b[i] <= p) {
        ret += 1;
        }
    }
    }

    return ret;
}

/* Poisson rv with mean m */

double genpois(const(double) m)
{
    double x;

    if (m > 200) {
    x = (m + 0.5) + sqrt(m) * gretl_one_snormal();
    x = floor(x);
    } else {
    int y = 0;

    x = exp(m) * gretl_rand_01();
    while (x > 1) {
        y++;
        x *= gretl_rand_01();
    }
    x = cast(double) y;
    }

    return x;
}

/**
 * gretl_rand_poisson:
 * @a: target array.
 * @t1: start of the fill range.
 * @t2: end of the fill range.
 * @m: mean (see below).
 * @vec: should be 1 if @m is an array, else 0.
 *
 * Fill the selected range of array @a with pseudo-random drawings
 * from the Poisson distribution with a mean determined by
 * @m, which can either be a pointer to a scalar, or an array
 * of length greater than or equal to @t2 + 1.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_rand_poisson(double * a, int t1, int t2, const(double) * m, int vec)
{
    double mt;
    int t;

    for (t=t1; t<=t2; t++) {
    mt = (vec)? m[t] : *m;
    a[t] = (mt <= 0)? NADBL : genpois(mt);
    }

    return 0;
}

/* f(x; k, \lambda) = (k/\lambda) (x/\lambda)^{k-1} e^{-(x/\lambda)^k} */

/**
 * gretl_rand_weibull:
 * @a: target array.
 * @t1: start of the fill range.
 * @t2: end of the fill range.
 * @shape: shape parameter > 0.
 * @scale: scale parameter > 0.
 *
 * Fill the selected range of array @a with pseudo-random drawings
 * from the Weibull distribution with shape @k and scale @lambda.
 *
 * Returns: 0 on success, non-zero if either parameter is out of
 * bounds.
 */

int gretl_rand_weibull(double * a, int t1, int t2, double shape, double scale)
{
    int err = 0;

    if (shape <= 0 || scale <= 0) {
    err = GretlError.data;
    } else {
    double u, kinv = 1.0 / shape;
    int t;

    for (t=t1; t<=t2; t++) {
        u = gretl_rand_01();
        while (u == 0.0) {
        u = gretl_rand_01();
        }
        a[t] = scale * pow(-log(u), kinv);
    }
    }

    return err;
}

/**
 * gretl_rand_exponential:
 * @a: target array.
 * @t1: start of the fill range.
 * @t2: end of the fill range.
 * @mu: scale parameter > 0.
 *
 * Fill the selected range of array @a with pseudo-random drawings
 * from the exponential distribution with scale @mu.
 *
 * Returns: 0 on success, non-zero if @mu is out of
 * bounds.
 */

int gretl_rand_exponential(double * a, int t1, int t2, double mu)
{
    int err = 0;

    if (mu <= 0) {
    err = GretlError.data;
    } else {
    double u;
    int t;

    for (t=t1; t<=t2; t++) {
        u = gretl_rand_01();
        while (u == 0.0) {
        u = gretl_rand_01();
        }
        a[t] = -mu * log(u);
    }
    }

    return err;
}

/**
 * gretl_rand_logistic:
 * @a: target array.
 * @t1: start of the fill range.
 * @t2: end of the fill range.
 * @loc: location parameter.
 * @scale: scale parameter > 0.
 *
 * Fill the selected range of array @a with pseudo-random drawings
 * from the logistic distribution with location @loc and scale @scale.
 *
 * Returns: 0 on success, non-zero if @scale is out of
 * bounds.
 */

int gretl_rand_logistic(double * a, int t1, int t2, double loc, double scale)
{
    int err = 0;

    if (scale <= 0) {
    err = GretlError.data;
    } else {
    double u;
    int t;

    for (t=t1; t<=t2; t++) {
        u = gretl_rand_01();
        while (u == 0.0) {
        u = gretl_rand_01();
        }
        a[t] = loc + scale * log(u / (1 - u));
    }
    }

    return err;
}

/**
 * gretl_rand_GED:
 * @a: target array.
 * @t1: start of the fill range.
 * @t2: end of the fill range.
 * @nu: shape parameter > 0.
 *
 * Fill the selected range of array @a with pseudo-random drawings
 * from the GED distribution with shape @nu. We exploit the fact that
 * if x ~ GED(n), then |x/k|^n is a Gamma rv.
 *
 * Returns: 0 on success, non-zero if @nu is out of bounds.
 */

int gretl_rand_GED(double * a, int t1, int t2, double nu)
{
    int err, t;
    double p, scale;

    if (nu < 0) {
    return GretlError.invarg;
    }

    p = 1.0/nu;
    scale = pow(0.5, p) * sqrt(gammafun(p) / gammafun(3.0*p));
    err = gretl_rand_gamma(a, t1, t2, p, 2);

    if (!err) {
    for (t=t1; t<=t2; t++) {
        a[t] = scale * pow(a[t], p);
        if (gretl_rand_01() < 0.5) {
        a[t] = -a[t];
        }
    }
    }

    return err;
}

/**
 * gretl_rand_laplace:
 * @a: target array.
 * @t1: start of the fill range.
 * @t2: end of the fill range.
 * @mu: mean.
 * @b: shape parameter > 0.
 *
 * Fill the selected range of array @a with pseudo-random drawings
 * from the Laplace distribution with mean @mu and scale @b.
 *
 * Returns: 0 on success, non-zero if @b is out of bounds.
 */

int gretl_rand_laplace(double * a, int t1, int t2, double mu, double b)
{
    int t, sgn;
    double U;

    if (b < 0) {
    return GretlError.invarg;
    }

    /* uniform on [0,1) */
    gretl_rand_uniform(a, t1, t2);

    for (t=t1; t<=t2; t++) {
    /* convert to (-1/2,1/2] */
    U = 0.5 - a[t];
    sgn = U < 0 ? -1 : 1;
    a[t] = mu - b*sgn * log(1 - 2*fabs(U));
    }

    return 0;
}

/**
 * gretl_rand_beta:
 * @x: target array.
 * @t1: start of the fill range.
 * @t2: end of the fill range.
 * @s1: shape parameter > 0.
 * @s2: shape parameter > 0.
 *
 * Fill the selected range of array @a with pseudo-random drawings
 * from the beta distribution with shape parameters @s1 and @s2.
 * The code here is adapted from http://www.netlib.org/random/random.f90
 * which implements the method of R.C.H. Cheng, "Generating beta
 * variates with nonintegral shape parameters", Communications of the
 * ACM, 21(4), April 1978.
 *
 * Returns: 0 on success, non-zero if @s1 or @s2 are out of bounds.
 */

int gretl_rand_beta(double * x, int t1, int t2, double s1, double s2)
{
    double aln4 = 1.3862944;
    double a, b, s, u, v, y, z;
    double d, f, h, u0, c;
    int t, swap;
    double val;

    if (s1 <= 0 || s2 <= 0) {
    return GretlError.data;
    }

    /* initialization */
    a = s1;
    b = s2;
    swap = b > a;
    if (swap) {
    f = b;
    b = a;
    a = f;
    }
    d = a/b;
    f = a+b;
    if (b > 1) {
    h = sqrt((2*a*b - f)/(f - 2.0));
    u0 = 1.0;
    } else {
    h = b;
    u0 = 1.0/(1.0 + pow(a/(DBL_MAX*b), b));
    }
    c = a+h;

    /* generation */
    for (t=t1; t<=t2; t++) {
    while (1) {
        u = gretl_rand_uniform_one();
        v = gretl_rand_uniform_one();
        s = u * u * v;
        if (u < DBL_MIN || s <= 0) continue;
        if (u < u0) {
        v = log(u/(1.0 - u))/h;
        y = d*exp(v);
        z = c*v + f*log((1.0 + d)/(1.0 + y)) - aln4;
        if (s - 1.0 > z) {
            if (s - s*z > 1.0) continue;
            if (log(s) > z) continue;
        }
        val = y/(1.0 + y);
        } else {
        if (4.0*s > pow(1.0 + 1.0/d, f)) continue;
        val = 1.0;
        }
        break;
    }
    x[t] = swap ? (1.0 - val) : val;
    }

    return 0;
}

/**
 * gretl_rand_beta_binomial:
 * @x: target array.
 * @t1: start of the fill range.
 * @t2: end of the fill range.
 * @n: number of trials.
 * @s1: beta shape parameter > 0.
 * @s2: beta shape parameter > 0.
 *
 * Fill the selected range of array @x with pseudo-random drawings
 * from the binomial distribution with @n trials and success
 * probability distributed according to the beta distribution with
 * shape parameters @s1 and @s2.
 *
 * Returns: 0 on success, non-zero if @n, @s1 or @s2 are out of bounds.
 */

int gretl_rand_beta_binomial(double * x, int t1, int t2, int n, double s1, double s2)
{
    int t, err;

    err = gretl_rand_beta(x, t1, t2, s1, s2);

    if (!err) {
    double * b = cast(double *) malloc(n * double.sizeof);

    if (b == null) {
        return GretlError.alloc;
    } else {
        for (t=t1; t<=t2; t++) {
        x[t] = gretl_rand_binomial_one(n, x[t], b);
        }
        free(b);
    }
    }

    return err;
}

/* Samples from the standard normal distribution truncated over the
   region [l,u], where 0 < l < u, using acceptance-rejection from
   Rayleigh distribution, similar to Marsaglia (1964).

   Use u = NADBL for setting a right-unbounded region (ie u = \infty)
*/

double ntail(double l, double u)
{
    double c, f, x, u01;
    int acc, iter = 0;

    c = 0.5 * l*l;
    f = na(u) ? -1.0 : expm1(c - 0.5 * u*u);

    /* sample from Rayleigh */
    x = c - log1p(gretl_rand_01() * f);
    u01 = gretl_rand_01();
    acc = (x*u01*u01) < c;

    while (!acc && iter < 1000) {
        /* find the threshold */
    // x = c - logp1(gretl_rand_01() * f);
    x = c - log1p(gretl_rand_01() * f);
    u01 = gretl_rand_01();
    acc = (x*u01*u01) < c;
        if (acc) {
        /* found! */
        break;
    }
        iter++;
    }

    if (iter == 1000) {
    x = NADBL;
    } else {
    /* this Rayleigh transform can be delayed till the end */
    x = sqrt(2*x);
    }

    return x;
}

/**
 * gretl_rand_tnormal:
 * @x: target array.
 * @t1: start of the fill range.
 * @t2: end of the fill range.
 * @lo: left truncation (NADBL for -infty).
 * @hi: right truncation (NADBL for +infty).
 *
 * Fill the selected range of array @x with pseudo-random drawings
 * from the truncated standard normal distribution between @lo and
 * @hi.
 *
 * Returns: 0 on success.
 */

int gretl_rand_tnormal(double * x, int t1, int t2, double l, double u)
{
    /*
      Efficient generator of a double from the standard normal
      distribution, truncated over the region [l,u]. Adapted from
      Zdravko Botev's Matlab code at
      https://web.maths.unsw.edu.au/~zdravkobotev/

      Remark: If you wish to simulate a random variable 'Z' from the
      non-standard Gaussian N(m,s^2) conditional on l < Z < u, then
      first simulate X = tmnormal((l-m)/s, (u-m)/s) and set Z = m+s*X

      Reference: Botev, Z. I. (2016). "The normal law under linear
      restrictions: simulation and estimation via minimax
      tilting". Journal of the Royal Statistical Society: Series B
      (Statistical Methodology).  doi:10.1111/rssb.12162
    */
    double a = 0.66;
    double tol = 2.0;
    double ret;
    int t;

    int ok_l = !na(l);
    int ok_u = !na(u);

    if (l >= u) {
    return GretlError.invarg;
    }

    /* thresholds for switching between methods; they can be tuned
       for maximum speed */

    if (!ok_l && !ok_u) {
    /* no limits! */
    for (t=t1; t<=t2; t++) {
        x[t] = gretl_one_snormal();
    }
    } else if (ok_l && (l > a)) {
    /* right tail */
    for (t=t1; t<=t2; t++) {
        x[t] = ntail(l, u);
    }
    } else if (ok_u && (u < -a)) {
    /* left tail */
    for (t=t1; t<=t2; t++) {
        x[t] = -ntail(-u, -l);
    }
    } else {
    if (ok_l && ok_u && fabs(l-u) > tol) {
        /* simulate from truncated normal via acceptance-rejection */
        for (t=t1; t<=t2; t++) {
        int out_ = 1;
        double z;

        while (out_) {
            z = gretl_one_snormal();
            out_ = (z < l) || (z > u);
        }
        x[t] = z;
        }
    } else {
        /* simulate from truncated normal via inversion */
        double P0 = ok_l ? normal_cdf(l) : 0.0;
        double P1 = ok_u ? normal_cdf(u) : 1.0;

        for (t=t1; t<=t2; t++) {
        x[t] = normal_cdf_inverse(P1 - (P1-P0) * gretl_rand_01());
        }
    }
    }

    return 0;
}

/* Start helper functions for gretl_rand_discrete(), using the alias
   method. */

double alias_pop(int * x, int * pn)
{
    *pn -= 1;
    return x[*pn];
}

void alias_push(int * x, int * pn, double xi)
{
    int n = *pn + 1;

    x[n-1] = cast(int) xi;
    *pn = n;
}

/* Construct lookup arrays for the alias method. For details see Luc
   Devroye's Non-Uniform Random Variate Generation (Springer-Verlag, New
   York, 1986) [ https://luc.devroye.org/chapter_three.pdf ].

   The originator of this method was A.J. Walker; see “New fast method
   for generating discrete random numbers wlth arbitrary frequency
   distributions,” Electronics Letters, vol. 10, pp. 127-128, 1974,
   and “An efficient method for generating discrete random
   variables with general distributions,” ACM Transactions on
   Mathematical Software, vol. 3, pp. 253-256, 1977.
*/

int alias_setup(const(GretlMatrixRaw) * p, int ** pA, double ** pR)
{
    int * A = null;
    double * R = null;
    int * smaller = null;
    int * greater = null;
    int sm, gr;
    int ns = 0;
    int ng = 0;
    int k, i;

    k = gretl_vector_get_length(p);
    *pA = A = cast(typeof(A)) malloc(k * (*A).sizeof);
    *pR = R = cast(typeof(R)) malloc(k * (*R).sizeof);
    smaller = cast(typeof(smaller)) malloc(2 * k * (*smaller).sizeof);
    greater = smaller + k;

    for (i=0; i<k; i++) {
    A[i] = 0;
    }

    /* scale probabilities by k */
    for (i=0; i<k; i++) {
    R[i] = k * p.val[i];
    if (R[i] < 1.0) {
        alias_push(smaller, &ns, i);
    } else {
        alias_push(greater, &ng, i);
    }
    }

    /* balance probabilities into pairs */
    while (ns > 0 && ng > 0) {
    sm = cast(int) alias_pop(smaller, &ns);
    gr = cast(int) alias_pop(greater, &ng);
    A[sm] = gr;
    R[gr] -= 1.0 - R[sm];
    if (R[gr] < 1.0) {
        alias_push(smaller, &ns, gr);
    } else {
        alias_push(greater, &ng, gr);
    }
    }

    free(smaller);

    return 0;
}

int alias_draw(const(int) * A, const(double) * R, int k)
{
    /* pick a random bucket */
    int i = cast(int) floor(gretl_rand_01() * k);
    double u = gretl_rand_01();

    /* decide between the primary item and the alias */
    return u < R[i] ? i : A[i];
}

/* End helper functions for gretl_rand_discrete(). */

/**
 * gretl_rand_discrete:
 * @x: target vector.
 * @t1: start of the fill range.
 * @t2: end of the fill range.
 * @p: probabilities vector.
 *
 * Fill the selected range of array @x with pseudo-random drawings
 * from the discrete distribution defined by @p, with support from 1
 * to @n, using the alias method.
 *
 * Returns: 0 on success, non-zero if @p is not a proper probability vector.
 */

int gretl_rand_discrete(double * x, int t1, int t2, const(GretlVector) * p)
{
    double psum = 0.0;
    int i, k;
    int err = 0;

    k = gretl_vector_get_length(p);

    /* sanity check on @p */
    for (i=0; i<k; i++) {
        if (p.val[i] < 0.0 || p.val[i] > 1.0) {
            fprintf(stderr, "gretl_rand_discrete: p is outside[0,1]\n".ptr);
            err = GretlError.invarg;
        break;
        }
    psum += p.val[i];
    }
    if (!err && fabs(psum - 1.0) > 1.0e-7) {
    err = GretlError.invarg;
    }

    if (!err) {
    int * A = null;
    double * R = null;
    int t;

    err = alias_setup(p, &A, &R);
    if (!err) {
        for (t=t1; t<=t2; t++) {
        x[t] = cast(double) (1 + alias_draw(A, R, k));
        }
    }
    free(A);
    free(R);
    }

    return err;
}

/**
 * gretl_rand_dirichlet:
 * @a: parameter vector, length k.
 * @n: number of rows (replications) in return matrix.
 * @err: location to receive error code.
 *
 * Returns: on success, an @n x @k matrix matrix containing
 * @n drawings from the Dirichlet distribution of order k.
 */

GretlMatrixRaw * gretl_rand_dirichlet(const(GretlVector) * a, int n, int * err)
{
    int k = gretl_vector_get_length(a);
    GretlMatrixRaw * D = null;

    if (k < 2) {
    *err = GretlError.nonconf;
    return null;
    }

    D = gretl_matrix_alloc(n, k);

    if (D == null) {
    *err = GretlError.alloc;
    } else {
    double dij, rsum;
    int t1 = 0, t2 = n-1;
    int i, j;

    for (j=0; j<k && !*err; j++) {
        *err = gretl_rand_gamma(D.val, t1, t2, a.val[j], 1.0);
        t1 += n;
        t2 += n;
    }
    if (*err) {
        gretl_matrix_free(D);
        D = null;
    } else {
        for (i=0; i<n; i++) {
        rsum = 0.0;
        for (j=0; j<k; j++) {
            rsum += gretl_matrix_get(D, i, j);
        }
        for (j=0; j<k; j++) {
            dij = gretl_matrix_get(D, i, j);
            gretl_matrix_set(D, i, j, dij / rsum);
        }
        }
    }
    }

    return D;
}

/**
 * gretl_rand_int_max:
 * @max: the maximum value (open)
 *
 * Returns: a pseudo-random unsigned int in the interval
 * [0, max-1].
 */

uint gretl_rand_int_max(uint max)
{
    return rand_int_range(0, max);
}

/**
 * gretl_rand_int:
 *
 * Returns: a pseudo-random unsigned int on the interval
 * [0, 2^32-1].
 */

uint gretl_rand_int()
{
    return randi32();
}

double halton(int i, int base)
{
    double f = 1.0 / base;
    double h = 0.0;

    while (i > 0) {
    h += f * (i % base);
    i = cast(int) floor(cast(double) (i / base));
    f /= base;
    }

    return h;
}

/* This is simple, and no doubt could be more efficient, but it's fast
   up to the number of primes that might reasonably be expected in the
   context of Halton sequences. The incoming @p is prime, and
   therefore odd. Starting at the next odd number q and skipping
   evens, we check for odd factors (since odd numbers cannot have even
   factors) from 3 to a maximum of q/3. If none are found then q is
   prime and we return it. There's a guard against running out of
   positive 32-bit integers; in that case we return 0.
*/

int next_prime(int p)
{
    int f, fmax, pr, q;

    for (q=p+2; q<INT_MAX-1; q+=2) {
    fmax = q/3;
    pr = 1;
    for (f=3; f<=fmax; f+=2) {
        if (q % f == 0) {
        pr = 0;
        break;
        }
    }
    if (pr) {
        return q;
    }
    }

    return 0;
}

/**
 * halton_matrix:
 * @m: number of rows (sequences).
 * @r: number of columns (elements in each sequence).
 * @offset: the number of initial values to discard.
 * @err: location to receive error code.
 *
 * Returns: an @m x @r matrix containing @m Halton
 * sequences. The sequences are contructed using the first
 * @m primes, and the first @offset elements of each sequence
 * are discarded.
 */

GretlMatrixRaw * halton_matrix(int m, int r, int offset, int * err)
{
    static immutable int[50] bases = [
    /* the first 50 primes */
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
    53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107,
    109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167,
    173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229
    ];
    GretlMatrixRaw * H;
    double hij;
    int i, j, k, n, p = 0;

    if (offset < 0 || m <= 0 || r <= 0) {
    *err = GretlError.invarg;
    return null;
    }

    H = gretl_matrix_alloc(m, r);
    if (H == null) {
    *err = GretlError.alloc;
    return null;
    }

    /* we'll discard the first @offset elements */
    n = r + offset;

    for (i=0; i<m; i++) {
    if (i < 50) {
        p = bases[i];
    } else {
        p = next_prime(p);
        if (p == 0) {
        *err = GretlError.invarg;
        break;
        }
    }
    j = 0;
    for (k=1; k<n; k++) {
        hij = halton(k, p);
        if (k >= offset) {
        gretl_matrix_set(H, i, j++, hij);
        }
    }
    }

    if (*err) {
    gretl_matrix_free(H);
    H = null;
    }

    return H;
}

int wishart_workspace(GretlMatrixRaw ** pW, GretlMatrixRaw ** pB, double ** pZ, int p)
{
    int err = 0;

    *pW = gretl_matrix_alloc(p, p);
    *pB = gretl_matrix_alloc(p, p);

    if (*pW == null || *pB == null) {
        err = GretlError.alloc;
    } else {
    int n = p * (p + 1) / 2;

    *pZ = cast(typeof(*pZ)) malloc(n * (**pZ).sizeof);
    if (*pZ == null) {
        err = GretlError.alloc;
    } else {
        gretl_rand_normal(*pZ, 0, n - 1);
    }
    }

    if (err) {
    gretl_matrix_free(*pW);
    gretl_matrix_free(*pB);
    }

    return err;
}

int wishart_inverse_finalize(const(GretlMatrixRaw) * C, GretlMatrixRaw * B, GretlMatrixRaw * W)
{
    gretl_matrix_qform(C, GretlMatrixMod.none,
               B, W, GretlMatrixMod.none);
    return gretl_invpd(W);
}

void odell_feiveson_compute(GretlMatrixRaw * B, double * Z, int v)
{
    double Xi, Zri, Zrj;
    double wii, wij;
    int p = B.rows;
    int i, j, r;

    for (i=0; i<p; i++) {
    gretl_rand_chisq(&Xi, 0, 0, v - i);
    wii = Xi;
    for (r=0; r<i; r++) {
        Zri = Z[ijton(r, i, p)];
        wii += Zri * Zri;
    }
    gretl_matrix_set(B, i, i, wii);
    for (j=i+1; j<p; j++) {
        wij = Z[ijton(i, j, p)] * sqrt(Xi);
        for (r=0; r<i; r++) {
        Zri = Z[ijton(r, i, p)];
        Zrj = Z[ijton(r, j, p)];
        wij += Zri * Zrj;
        }
        gretl_matrix_set(B, i, j, wij);
        gretl_matrix_set(B, j, i, wij);
    }
    }
}

/**
 * inverse_wishart_matrix:
 * @S: p x p positive definite scale matrix.
 * @v: degrees of freedom: v >= p.
 * @err: location to receive error code.
 *
 * Computes a draw from the Inverse Wishart distribution
 * with @v degrees of freedom, using the method of Odell and
 * Feiveson, "A numerical procedure to generate a sample
 * covariance matrix", JASA 61, pp. 199­203, 1966.
 *
 * Returns: a p x p Inverse Wishart matrix, or NULL on error.
 */

GretlMatrixRaw * inverse_wishart_matrix(const(GretlMatrixRaw) * S, int v, int * err)
{
    GretlMatrixRaw * C, B = null, W = null;
    double * Z = null;

    if (S == null || S.cols != S.rows || v < S.rows) {
    *err = GretlError.invarg;
    return null;
    }

    C = gretl_matrix_copy(S);
    if (C == null) {
        *err = GretlError.alloc;
        return null;
    }

    /* invert and decompose C */
    *err = cholesky_factor_of_inverse(C);

    if (!*err) {
    *err = wishart_workspace(&W, &B, &Z, S.rows);
    }

    if (*err) {
    gretl_matrix_free(C);
    return null;
    }

    odell_feiveson_compute(B, Z, v);
    *err = wishart_inverse_finalize(C, B, W);

    if (*err) {
    gretl_matrix_free(W);
    W = null;
    }

    gretl_matrix_free(B);
    gretl_matrix_free(C);
    free(Z);

    return W;
}

void vech_into_row(GretlMatrixRaw * targ, int row, const(GretlMatrixRaw) * m)
{
    double mij;
    int n = m.rows;
    int i, j, jj = 0;

    for (i=0; i<n; i++) {
    for (j=i; j<n; j++) {
        mij = gretl_matrix_get(m, i, j);
        gretl_matrix_set(targ, row, jj++, mij);
    }
    }
}

/**
 * inverse_wishart_sequence:
 * @S: p x p positive definite scale matrix.
 * @v: degrees of freedom: v >= p.
 * @n: number of replications.
 * @err: location to receive error code.
 *
 * Computes @n draws from the Inverse Wishart distribution
 * with @v degrees of freedom, using the method of Odell and
 * Feiveson, "A numerical procedure to generate a sample
 * covariance matrix", JASA 61, pp. 199­203, 1966.
 *
 * Returns: an n x k matrix, where k = p*(p+1)/2, each row
 * being the vech of an Inverse Wishart matrix, or NULL on error.
 */

GretlMatrixRaw * inverse_wishart_sequence(const(GretlMatrixRaw) * S, int v, int n, int * err)
{
    GretlMatrixRaw * C, B = null, W = null;
    GretlMatrixRaw * Seq = null;
    double * Z = null;
    int i, k = 0;

    if (S == null || S.cols != S.rows || v < S.rows) {
    *err = GretlError.invarg;
    return null;
    }

    if (n < 1) {
    *err = GretlError.invarg;
    return null;
    }

    C = gretl_matrix_copy(S);
    if (C == null) {
        *err = GretlError.alloc;
        return null;
    }

    /* invert and decompose C */
    *err = cholesky_factor_of_inverse(C);

    if (!*err) {
    *err = wishart_workspace(&W, &B, &Z, S.rows);
    }

    if (!*err) {
    k = C.rows * (C.rows + 1) / 2;
    Seq = gretl_matrix_alloc(n, k);
    if (Seq == null) {
        *err = GretlError.alloc;
    }
    }

    for (i=0; i<n && !*err; i++) {
    odell_feiveson_compute(B, Z, v);
    *err = wishart_inverse_finalize(C, B, W);
        if (!*err) {
        /* write vech of W into row i of Seq */
        vech_into_row(Seq, i, W);
        if (i < n - 1) {
        /* refresh Normal draws */
        gretl_rand_normal(Z, 0, k - 1);
        }
    }
    }

    gretl_matrix_free(C);
    gretl_matrix_free(W);
    gretl_matrix_free(B);
    free(Z);

    if (*err && Seq != null) {
    gretl_matrix_free(Seq);
    Seq = null;
    }

    return Seq;
}

char * gretl_rand_hex_string(int len, int * err)
{
    const(char) * src = "0123456789abcdef";
    int * ivals = null;
    char * ret = null;
    int i;

    if (len < 0) {
    *err = GretlError.invarg;
    return null;
    } else if (len == 0) {
    return gretl_strdup("");
    }

    ivals = cast(typeof(ivals)) malloc(len * (*ivals).sizeof);
    ret = cast(typeof(ret)) malloc(len + 1);

    if (ivals == null || ret == null) {
    *err = GretlError.alloc;
    return ret;
    }

    gretl_rand_int_minmax(ivals, len, 0, 15);

    for (i=0; i<len; i++) {
    ret[i] = src[ivals[i]];
    }

    ret[i] = '\0';
    free(ivals);

    return ret;
}

