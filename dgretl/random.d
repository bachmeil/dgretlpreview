/* Random numbers: the generators of the library (gretl.random: the xoshiro256+ engine with Ziggurat normals, gamma, binomial, Poisson, beta, ... generators) and
   a family of parallel streams that can drive all of them.

   import dgretl;
   setSeed(42);                                     // the library's own generator: sequential, reproducible from the seed
   auto x = rnorm(1000, mean: 1, sd: 2);            // double[]
   auto z = randomMatrix!rgamma(100, 3, 2.0);       // a 100 x 3 GretlMatrix of gamma(2) draws

   auto streams = makeStreams(8, seed: 2024);       // eight non-overlapping streams of L'Ecuyer's MRG31k3p generator
   auto draws = parallelDraws(streams, (size_t i) => rgamma(10_000, 2.0));       // each task draws from its own stream, in parallel
   auto more = streams[3].rnorm(100);               // one stream, one distribution

   Sequential use is the library's engine (xoshiro256+): fast, one global state, so not for use by several threads at once. For parallel simulation the library's engine is
   replaced, for the thread and for the time of a draw, by one stream of the MRG31k3p generator (a port of Pierre L'Ecuyer's generator as in SSJ: dgretl.prng), whose
   streams start 2^134 steps apart and so do not overlap; every distribution of the library then draws from the stream, and the global generator is not touched.
   The result of a stream does not depend on the thread that runs it or on the other streams, so a parallel simulation is reproducible from the seed whatever the number of
   processors. A stream must be used by one thread at a time.

   Invalid parameters (a negative standard deviation, a probability outside [0, 1], ...) throw GretlException. Draws are `double[]` (`int[]` for the integer-valued uniform and the
   discrete distribution), or a GretlMatrix for the multivariate ones and for randomMatrix. */

module dgretl.random;

public import dgretl.matrix;
public import dgretl.exception;

import std.conv: to;
import std.exception: enforce;
import std.math: isFinite, isNaN;
import std.parallelism: parallel, taskPool;
import std.range: iota;

import dgretl.internal;
import dgretl.prng: Stream, StreamGen;
import gretl.gretl_matrix: gretl_matrix_free;
import gretl.gretl_matrixh: GretlMatrixRaw;
import gretl.libgretlh: Gretlopt;
import gretl.missingh: NADBL;
import gretl.random: RandSource, gretl_rand_beta, gretl_rand_beta_binomial, gretl_rand_binomial, gretl_rand_chisq, gretl_rand_clear_source, gretl_rand_dirichlet,
    gretl_rand_discrete, gretl_rand_exponential, gretl_rand_F, gretl_rand_gamma, gretl_rand_GED, gretl_rand_get_seed, gretl_rand_get_source, gretl_rand_int_minmax,
    gretl_rand_laplace, gretl_rand_logistic, gretl_rand_normal_full, gretl_rand_poisson, gretl_rand_prepare, gretl_rand_set_seed, gretl_rand_set_source,
    gretl_rand_student, gretl_rand_tnormal, gretl_rand_uniform_int_minmax, gretl_rand_uniform_minmax, gretl_rand_weibull, inverse_wishart_matrix,
    inverse_wishart_sequence;

/* ---- the library's generator ---- */

/* setSeed(seed): seed the library's generator (the one that the functions below use outside a stream) so that what follows is reproducible. 0 chooses a seed from the system. */
void setSeed(ulong seed) @trusted {
    gretl_rand_set_seed(seed);
}

/* getSeed(): the seed of the library's generator. */
ulong getSeed() @trusted {
    return gretl_rand_get_seed();
}

/* ---- the distributions ---- */

private double[] filled(string who, size_t n, scope int delegate(double *, int) @system f) @trusted {
    enforce!GretlException(n < int.max, who ~ ": too many draws");
    auto a = new double[n];

    if (n == 0) return a;
    a[] = 0;
    const err = f(a.ptr, cast(int) n - 1);

    enforce!GretlException(err == 0, who ~ ": invalid parameters");
    return a;
}

private void need(bool ok, string who, string what) {
    enforce!GretlException(ok, who ~ ": " ~ what);
}

/* runif(n, min, max): n draws from the uniform distribution on [min, max) (max above min). */
double[] runif(size_t n, double min = 0, double max = 1) {
    need(isFinite(min) && isFinite(max) && max > min, "runif", "need finite min < max");
    return filled("runif", n, (double * a, int last) => gretl_rand_uniform_minmax(a, 0, last, min, max));
}

/* runifInt(n, min, max, withoutReplacement): n integers uniformly from min .. max (both included). With withoutReplacement the draws are distinct
   (n must not exceed max - min + 1), like a sample of that many from the range. */
int[] runifInt(size_t n, int min, int max, bool withoutReplacement = false) @trusted {
    need(max >= min, "runifInt", "need min <= max");
    enforce!GretlException(n < int.max, "runifInt: too many draws");
    auto r = new int[n];

    if (n == 0) return r;
    if (withoutReplacement) {
        need(n <= cast(size_t) (cast(long) max - min + 1), "runifInt", "cannot draw " ~ to!string(n) ~ " distinct integers from " ~ to!string(min) ~ " .. " ~ to!string(max));
        auto a = new double[n];

        a[] = 0;
        const err = gretl_rand_uniform_int_minmax(a.ptr, 0, cast(int) n - 1, min, max, Gretlopt.OPT_O);

        need(err == 0, "runifInt", "invalid parameters");
        foreach (i; 0 .. n) r[i] = cast(int) a[i];
    } else {
        const err = gretl_rand_int_minmax(r.ptr, cast(int) n, min, max);

        need(err == 0, "runifInt", "invalid parameters");
    }
    return r;
}

/* rnorm(n, mean, sd): the normal distribution (Ziggurat method). */
double[] rnorm(size_t n, double mean = 0, double sd = 1) {
    need(isFinite(mean) && isFinite(sd) && sd > 0, "rnorm", "need a finite mean and a positive sd");
    return filled("rnorm", n, (double * a, int last) => gretl_rand_normal_full(a, 0, last, mean, sd));
}

/* rtnorm(n, lower, upper, mean, sd): the normal distribution with the given mean and standard deviation truncated to [lower, upper] (use -double.infinity and
   double.infinity for no bound on that side). */
double[] rtnorm(size_t n, double lower = -double.infinity, double upper = double.infinity, double mean = 0, double sd = 1) {
    need(isFinite(mean) && isFinite(sd) && sd > 0, "rtnorm", "need a finite mean and a positive sd");
    need(!isNaN(lower) && !isNaN(upper) && lower < upper, "rtnorm", "need lower < upper");
    // the library draws the standard normal between its bounds (NADBL for none)
    const lo = isFinite(lower) ? (lower - mean) / sd : NADBL, hi = isFinite(upper) ? (upper - mean) / sd : NADBL;
    auto r = filled("rtnorm", n, (double * a, int last) => gretl_rand_tnormal(a, 0, last, lo, hi));

    foreach (ref v; r) v = mean + sd * v;
    return r;
}

/* rgamma(n, shape, scale): the gamma distribution with mean shape * scale. */
double[] rgamma(size_t n, double shape, double scale = 1) {
    need(shape > 0 && scale > 0 && isFinite(shape) && isFinite(scale), "rgamma", "need a positive shape and scale");
    return filled("rgamma", n, (double * a, int last) => gretl_rand_gamma(a, 0, last, shape, scale));
}

/* rchisq(n, df): the chi-square distribution with df (integer) degrees of freedom. */
double[] rchisq(size_t n, int df) {
    need(df > 0, "rchisq", "need positive degrees of freedom");
    return filled("rchisq", n, (double * a, int last) => gretl_rand_chisq(a, 0, last, df));
}

/* rt(n, df): Student's t distribution with df degrees of freedom. */
double[] rt(size_t n, double df) {
    need(df > 0 && isFinite(df), "rt", "need positive degrees of freedom");
    return filled("rt", n, (double * a, int last) => gretl_rand_student(a, 0, last, df));
}

/* rF(n, df1, df2): the F distribution (integer degrees of freedom). */
double[] rF(size_t n, int df1, int df2) {
    need(df1 > 0 && df2 > 0, "rF", "need positive degrees of freedom");
    return filled("rF", n, (double * a, int last) => gretl_rand_F(a, 0, last, df1, df2));
}

/* rbinom(n, size, p): the binomial distribution with `size` trials and success probability p (the draws are doubles holding integers). */
double[] rbinom(size_t n, int size, double p) {
    need(size >= 0 && p >= 0 && p <= 1, "rbinom", "need size >= 0 and 0 <= p <= 1");
    return filled("rbinom", n, (double * a, int last) => gretl_rand_binomial(a, 0, last, size, p));
}

/* rbetabinom(n, size, a, b): the beta-binomial distribution: `size` trials, success probability drawn from the beta(a, b) distribution. */
double[] rbetabinom(size_t n, int size, double a, double b) {
    need(size >= 0 && a > 0 && b > 0, "rbetabinom", "need size >= 0 and positive shape parameters");
    return filled("rbetabinom", n, (double * x, int last) => gretl_rand_beta_binomial(x, 0, last, size, a, b));
}

/* rpois(n, lambda): the Poisson distribution with mean lambda. */
double[] rpois(size_t n, double lambda) {
    need(lambda >= 0 && isFinite(lambda), "rpois", "need a non-negative mean");
    double m = lambda;

    return filled("rpois", n, (double * a, int last) => gretl_rand_poisson(a, 0, last, &m, 0));
}

/* rpois(lambdas): one Poisson draw for each mean in the array. */
double[] rpois(const(double)[] lambdas) @trusted {
    foreach (l; lambdas) need(l >= 0 && isFinite(l), "rpois", "need non-negative means");
    const(double)[] m = lambdas;

    return filled("rpois", lambdas.length, (double * a, int last) => gretl_rand_poisson(a, 0, last, m.ptr, 1));
}

/* rweibull(n, shape, scale): the Weibull distribution. */
double[] rweibull(size_t n, double shape, double scale = 1) {
    need(shape > 0 && scale > 0 && isFinite(shape) && isFinite(scale), "rweibull", "need a positive shape and scale");
    return filled("rweibull", n, (double * a, int last) => gretl_rand_weibull(a, 0, last, shape, scale));
}

/* rexp(n, mean): the exponential distribution with the given mean (the library's `scale`). */
double[] rexp(size_t n, double mean = 1) {
    need(mean > 0 && isFinite(mean), "rexp", "need a positive mean");
    return filled("rexp", n, (double * a, int last) => gretl_rand_exponential(a, 0, last, mean));
}

/* rlogis(n, location, scale): the logistic distribution. */
double[] rlogis(size_t n, double location = 0, double scale = 1) {
    need(isFinite(location) && scale > 0 && isFinite(scale), "rlogis", "need a finite location and a positive scale");
    return filled("rlogis", n, (double * a, int last) => gretl_rand_logistic(a, 0, last, location, scale));
}

/* rged(n, nu): the generalized error distribution with shape nu (the library's parameterization; nu = 2 is the normal distribution). */
double[] rged(size_t n, double nu) {
    need(nu > 0 && isFinite(nu), "rged", "need a positive shape");
    return filled("rged", n, (double * a, int last) => gretl_rand_GED(a, 0, last, nu));
}

/* rlaplace(n, mu, b): the Laplace distribution with mean mu and scale b. */
double[] rlaplace(size_t n, double mu = 0, double b = 1) {
    need(isFinite(mu) && b > 0 && isFinite(b), "rlaplace", "need a finite mean and a positive scale");
    return filled("rlaplace", n, (double * a, int last) => gretl_rand_laplace(a, 0, last, mu, b));
}

/* rbeta(n, a, b): the beta distribution with shape parameters a and b. */
double[] rbeta(size_t n, double a, double b) {
    need(a > 0 && b > 0 && isFinite(a) && isFinite(b), "rbeta", "need positive shape parameters");
    return filled("rbeta", n, (double * x, int last) => gretl_rand_beta(x, 0, last, a, b));
}

/* rdiscrete(n, probs): n draws from the distribution on 0 .. probs.length - 1 with these probabilities (non-negative, adding to 1), by the alias method. (The library's
   values run from 1; these are one less.) */
size_t[] rdiscrete(size_t n, const(double)[] probs) @trusted {
    need(probs.length >= 1, "rdiscrete", "need at least one probability");
    double total = 0;

    foreach (p; probs) {
        need(p >= 0 && isFinite(p), "rdiscrete", "probabilities must be non-negative");
        total += p;
    }
    need(total > 1 - 1e-9 && total < 1 + 1e-9, "rdiscrete", "the probabilities must add to 1");
    auto pv = fresh(probs.length, 1);

    foreach (i, p; probs) pv.ptr.val[i] = p;
    auto d = filled("rdiscrete", n, (double * a, int last) => gretl_rand_discrete(a, 0, last, pv.ptr));
    auto r = new size_t[n];

    foreach (i; 0 .. n) r[i] = cast(size_t) d[i] - 1;
    return r;
}

/* rdirichlet(alpha, n): n draws from the Dirichlet distribution with parameter vector alpha (positive), as the rows of an n x k matrix. */
GretlMatrix rdirichlet(const(double)[] alpha, size_t n) @trusted {
    need(alpha.length >= 2 && n >= 1, "rdirichlet", "need at least two parameters and one draw");
    foreach (a; alpha) need(a > 0 && isFinite(a), "rdirichlet", "parameters must be positive");
    auto av = fresh(alpha.length, 1);

    foreach (i, a; alpha) av.ptr.val[i] = a;
    int err = 0;
    auto m = gretl_rand_dirichlet(av.ptr, cast(int) n, &err);

    need(err == 0 && m !is null, "rdirichlet", "invalid parameters");
    scope (exit) gretl_matrix_free(m);
    return copyRaw(m);
}

/* rinvwishart(S, df): one draw from the inverse Wishart distribution with p x p scale matrix S (positive definite) and df >= p degrees of freedom (Odell and Feiveson's
   method), as a p x p matrix. */
GretlMatrix rinvwishart(M)(M scale, int df) @trusted if (isMatrixLike!M) {
    auto s = cbind(scale);

    need(s.nrow == s.ncol && s.nrow >= 1 && df >= cast(int) s.nrow, "rinvwishart", "need a square scale matrix and df >= its size");
    requireFinite(s, "rinvwishart: scale");
    int err = 0;
    auto m = inverse_wishart_matrix(s.ptr, df, &err);

    need(err == 0 && m !is null, "rinvwishart", "the scale matrix must be positive definite");
    scope (exit) gretl_matrix_free(m);
    return copyRaw(m);
}

/* rinvwishart(S, df, n): n draws, as the rows of an n x p(p + 1)/2 matrix, each row the vech (the lower triangle, column by column) of an inverse Wishart matrix. */
GretlMatrix rinvwishart(M)(M scale, int df, size_t n) @trusted if (isMatrixLike!M) {
    auto s = cbind(scale);

    need(s.nrow == s.ncol && s.nrow >= 1 && df >= cast(int) s.nrow && n >= 1, "rinvwishart", "need a square scale matrix, df >= its size and n >= 1");
    requireFinite(s, "rinvwishart: scale");
    int err = 0;
    auto m = inverse_wishart_sequence(s.ptr, df, cast(int) n, &err);

    need(err == 0 && m !is null, "rinvwishart", "the scale matrix must be positive definite");
    scope (exit) gretl_matrix_free(m);
    return copyRaw(m);
}

private GretlMatrix copyRaw(const(GretlMatrixRaw) * m) @trusted {
    auto r = fresh(m.rows, m.cols);

    foreach (i; 0 .. cast(size_t) m.rows * m.cols) r.ptr.val[i] = m.val[i];
    return r;
}

/* randomMatrix!rnorm(nrow, ncol, mean, sd): an nrow x ncol GretlMatrix filled (column by column) with draws of any of the functions above that return a double[]
   (the arguments after nrow and ncol are theirs after n). */
GretlMatrix randomMatrix(alias f, Args...)(size_t nrow, size_t ncol, Args args) @trusted {
    auto d = f(nrow * ncol, args);
    auto m = fresh(nrow, ncol);

    foreach (i; 0 .. nrow * ncol) m.ptr.val[i] = d[i];
    return m;
}

/* ---- parallel streams: L'Ecuyer's MRG31k3p ---- */

/* The streams are those of prng.d (dgretl.prng), a port of MRG31k3p.java of the SSJ project (https://github.com/umontreal-simul/ssj): a combined multiple recursive
   generator with two components of order 3 and moduli 2^31 - 1 and 2^31 - 21069, period about 2^185. A family of streams is made by jumping 2^134 steps from the last
   one, so the streams never overlap in practice. The licence of the original is Apache 2.0. */

/* One stream of random numbers: a prng.d `Stream` (six integers of state, `gen()` for U(0, 1)) that can also give 64 random bits and drive the library's generators. Every
   distribution function of this module is also a method of it (stream.rnorm(100), stream.rgamma(50, 2.0), ...) and then draws from this stream. */
final class RandomStream : Stream {
    this(int[6] state) {
        super(state);
    }

    /* the state (the six integers), to start another stream from it */
    int[6] state() const {
        return [x11, x12, x13, x21, x22, x23];
    }

    /* the next U(0, 1) draw (never 0 or 1; 31 bits of resolution) */
    double uniform() {
        return gen();
    }

    /* 64 random bits, from three steps of the generator (each step has 31 bits): what the library's generators ask of their engine */
    ulong next64() {
        const a = cast(ulong) (gen() * 2147483648.0), b = cast(ulong) (gen() * 2147483648.0), c = cast(ulong) (gen() * 2147483648.0);

        return a | (b << 31) | ((c & 3) << 62);
    }

    private static ulong thunk(void * p) {
        return (cast(RandomStream) p).next64();
    }

    /* run f with this stream as the source of every random draw of the library, in this thread (the draws of the global generator are not used meanwhile) */
    void use(scope void delegate() f) {
        auto guard = StreamGuard(this);

        f();
    }

    /* the distribution functions of this module, drawing from this stream */
    auto opDispatch(string name, Args...)(Args args) if (isDistribution!name) {
        auto guard = StreamGuard(this);

        return mixin(name ~ "(args)");
    }
}

private enum distributions = ["runif", "runifInt", "rnorm", "rtnorm", "rgamma", "rchisq", "rt", "rF", "rbinom", "rbetabinom", "rpois", "rweibull", "rexp", "rlogis", "rged",
                              "rlaplace", "rbeta", "rdiscrete", "rdirichlet", "rinvwishart", "randomMatrix"];

private template isDistribution(string name) {
    private bool has() {
        foreach (d; distributions) if (d == name) return true;
        return false;
    }
    enum isDistribution = has();
}

private struct StreamGuard {
    RandSource prevFn;
    void * prevCtx;

    @disable this(this);

    this(RandomStream s) @trusted {
        gretl_rand_get_source(&prevFn, &prevCtx);
        gretl_rand_prepare();
        gretl_rand_set_source(&RandomStream.thunk, cast(void *) s);
    }

    ~this() @trusted {
        gretl_rand_set_source(prevFn, prevCtx);
    }
}

/* a valid starting state for a nonzero seed: six integers, the first three in 1 .. M1 - 1 and the last three in 1 .. M2 - 1, from splitmix64 */
private int[6] startState(ulong seed) {
    int[6] state;
    ulong z = seed;

    ulong next() {
        z += 0x9e3779b97f4a7c15UL;
        ulong r = z;

        r = (r ^ (r >> 30)) * 0xbf58476d1ce4e5b9UL;
        r = (r ^ (r >> 27)) * 0x94d049bb133111ebUL;
        return r ^ (r >> 31);
    }
    foreach (i; 0 .. 6) {
        const m = i < 3 ? StreamGen.M1 : StreamGen.M2;

        state[i] = cast(int) (1 + next() % (m - 1));
    }
    return state;
}

/* makeStreams(n, seed): n non-overlapping streams, made by prng.d's StreamGen. With seed 0 (the default) the first starts from the generator's original state (all six
   integers 12345); another seed starts the family from a state derived from it. */
RandomStream[] makeStreams(size_t n, ulong seed = 0) {
    auto g = (seed == 0) ? StreamGen.init : StreamGen(startState(seed));
    RandomStream[] r;

    foreach (i; 0 .. n) {
        auto s = g.createStream();

        r ~= new RandomStream([s.x11, s.x12, s.x13, s.x21, s.x22, s.x23]);
    }
    gretl_rand_prepare();
    return r;
}

/* parallelDraws(streams, f): calls f(i) for i = 0, 1, ... in parallel on the task pool of std.parallelism, task i with stream i as the source of every random draw
   (the functions of this module, called inside f, draw from it), and returns the results in order. The results do not depend on the number of processors or on how the tasks
   are scheduled. */
T[] parallelDraws(T)(RandomStream[] streams, scope T delegate(size_t) f) {
    auto result = new T[streams.length];

    gretl_rand_prepare();
    foreach (i; parallel(iota(streams.length))) {
        auto guard = StreamGuard(streams[i]);

        result[i] = f(i);
    }
    return result;
}

/// ditto: for a function that needs the stream itself as well
T[] parallelDraws(T)(RandomStream[] streams, scope T delegate(RandomStream, size_t) f) {
    return parallelDraws!T(streams, (size_t i) => f(streams[i], i));
}
