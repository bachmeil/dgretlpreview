/* Tests for dgretl/random.d: the library's generators as D functions, the thread-local source hook of gretl/random.d, the MRG31k3p streams and the parallel draws. */

import core.memory: GC;
import core.stdc.stdlib;
import std.conv: to;
import std.math: PI, abs, exp, isNaN, sqrt;
import std.mathspecial: normalDistribution, gamma;
import std.parallelism: totalCPUs;
import std.stdio: writefln, writeln;
import dgretl;
import dgretl.prng: Stream, StreamGen, genDouble, genInteger;
import gretl.random: double_from_uint64, gretl_rand_01, gretl_rand_clear_source, gretl_rand_set_source, randi32, RandSource;

extern (C) {
    struct mallinfo_t { int arena, ordblks, smblks, hblks, hblkhd, usmblks, fsmblks, uordblks, fordblks, keepcost; }
    mallinfo_t mallinfo();
}

int nTests, nFailed;

void check(bool ok, string what) {
    nTests++;
    if (!ok) nFailed++;
    writefln("%s %s", ok ? "ok:  " : "FAIL:", what);
}

bool close(double a, double b, double tol = 1e-9) { return abs(a - b) <= tol * (1.0 + abs(b)); }

bool throws(void delegate() f) {
    try f(); catch (GretlException) return true;
    return false;
}

double meanOf(const(double)[] x) { double s = 0; foreach (v; x) s += v; return s / x.length; }
double varOf(const(double)[] x) { const m = meanOf(x); double s = 0; foreach (v; x) s += (v - m) ^^ 2; return s / (x.length - 1); }

/* the sample mean and variance are within 5 standard errors of the theory (the variance's standard error from the fourth moment, estimated) */
bool moments(const(double)[] x, double mean, double variance, string what) {
    const n = x.length;
    const m = meanOf(x), v = varOf(x);
    double m4 = 0;
    foreach (e; x) m4 += (e - m) ^^ 4;
    m4 /= n;
    const seMean = sqrt(variance / n), seVar = sqrt((m4 - v * v) / n);
    const ok = abs(m - mean) < 5 * seMean && abs(v - variance) < 5 * seVar;
    if (!ok) writefln("  %s: mean %g (theory %g, se %g), variance %g (theory %g, se %g)", what, m, mean, seMean, v, variance, seVar);
    return ok;
}

/* a source for the hook test: returns the values of an array in turn */
__gshared ulong[] hookValues;
__gshared size_t hookPos;
extern (D) ulong hookNext(void * ctx) { return hookValues[hookPos++ % hookValues.length]; }

void main() {
    enum N = 200_000;

    // ---- the source hook of gretl.random ----
    {
        setSeed(5);
        const g1 = gretl_rand_01(), g2 = gretl_rand_01();
        setSeed(5);
        const h1 = gretl_rand_01();
        hookValues = [0x8000000000000000UL, 0xFFFFFFFFFFFFFFFFUL, 0x0000000000000000UL];
        hookPos = 0;
        gretl_rand_set_source(&hookNext, null);
        const a = gretl_rand_01(), b = gretl_rand_01(), c = gretl_rand_01(), d32 = randi32();
        gretl_rand_clear_source();
        const h2 = gretl_rand_01();
        check(g1 == h1 && g2 == h2, "the hook leaves the global generator's state alone: the sequence continues as if it had not been used");
        check(a == double_from_uint64(0x8000000000000000UL) && close(a, 0.5, 1e-15) && b == double_from_uint64(0xFFFFFFFFFFFFFFFFUL) && b < 1 && c == 0.0, "while a source is set, gretl_rand_01 takes its bits from it");
        check(d32 == 0x80000000u, "randi32 takes the top 32 bits of the source");
    }

    // ---- the global generator ----
    setSeed(42);
    const r1 = rnorm(10);
    setSeed(42);
    const r2 = rnorm(10);
    check(r1 == r2 && getSeed() == 42, "setSeed: the same seed gives the same draws, getSeed returns it");
    setSeed(43);
    check(rnorm(10) != r1, "another seed, other draws");
    check(rnorm(0).length == 0 && runif(0).length == 0 && rgamma(0, 2.0).length == 0, "n = 0 gives an empty array");

    setSeed(2024);
    {
        auto u = runif(N, 2, 5);
        check(moments(u, 3.5, 9.0 / 12, "runif") && minimum(u) >= 2 && maximum(u) < 5, "runif(2, 5): moments, within [2, 5)");
        auto i = runifInt(N, -3, 4);
        double[] id;
        foreach (v; i) id ~= v;
        check(minimum(id) == -3 && maximum(id) == 4 && moments(id, 0.5, (8.0 ^^ 2 - 1) / 12, "runifInt"), "runifInt(-3, 4): range and moments");
        auto nr = runifInt(10, 1, 10, true);
        bool[11] seen;
        bool distinct = true;
        foreach (v; nr) { distinct = distinct && !seen[v]; seen[v] = true; }
        check(distinct && nr.length == 10, "runifInt without replacement: ten distinct values from 1 .. 10");
        check(throws({ runifInt(11, 1, 10, true); }) && throws({ runifInt(3, 5, 4); }), "runifInt: too many distinct values, or min above max, throw");
        check(moments(rnorm(N, 1, 2), 1, 4, "rnorm"), "rnorm(1, 2)");
        auto tn = rtnorm(N, 0, double.infinity);
        const tm = sqrt(2 / PI);
        check(minimum(tn) >= 0 && moments(tn, tm, 1 - 2 / PI, "rtnorm one-sided"), "rtnorm(0, inf): non-negative, half-normal moments");
        auto tn2 = rtnorm(N, -1, 2, 1, 2);       // N(1, 2^2) between -1 and 2: standardized bounds -1 and 0.5
        const al = -1.0, be = 0.5;
        const phi = (double z) => exp(-z * z / 2) / sqrt(2 * PI);
        const zz = normalDistribution(be) - normalDistribution(al);
        const mt = 1 + 2 * (phi(al) - phi(be)) / zz;
        const vt = 4 * (1 + (al * phi(al) - be * phi(be)) / zz - ((phi(al) - phi(be)) / zz) ^^ 2);
        check(minimum(tn2) >= -1 && maximum(tn2) <= 2 && moments(tn2, mt, vt, "rtnorm two-sided"), "rtnorm(-1, 2, mean 1, sd 2): within the bounds, truncated-normal moments");
        check(moments(rgamma(N, 3, 2), 6, 12, "rgamma 3") && moments(rgamma(N, 0.5, 2), 1, 2, "rgamma 0.5"), "rgamma(3, 2) and rgamma(0.5, 2)");
        check(moments(rchisq(N, 4), 4, 8, "rchisq"), "rchisq(4)");
        check(moments(rt(N, 10), 0, 10.0 / 8, "rt"), "rt(10)");
        check(moments(rF(N, 5, 12), 1.2, 2.0 * 144 * 15 / (5 * 100 * 8), "rF"), "rF(5, 12)");
        check(moments(rbinom(N, 20, 0.3), 6, 4.2, "rbinom"), "rbinom(20, 0.3)");
        check(moments(rbetabinom(N, 10, 2, 3), 4, 10.0 * 6 * 15 / (25 * 6), "rbetabinom"), "rbetabinom(10, 2, 3)");
        check(moments(rpois(N, 4.5), 4.5, 4.5, "rpois") && moments(rpois(N, 40), 40, 40, "rpois large"), "rpois(4.5) and rpois(40)");
        double[] lams;
        foreach (k; 0 .. 100_000) lams ~= (k % 2 == 0) ? 1.0 : 10.0;
        auto pv = rpois(lams);
        double[] odd, even;
        foreach (k, v; pv) (k % 2 == 0 ? even : odd) ~= v;
        check(pv.length == 100_000 && moments(even, 1, 1, "rpois vec 1") && moments(odd, 10, 10, "rpois vec 10"), "rpois with an array of means: each element has its own mean");
        const wm = 3 * gamma(1.5), wv = 9 * (gamma(2.0) - gamma(1.5) ^^ 2);
        check(moments(rweibull(N, 2, 3), wm, wv, "rweibull"), "rweibull(2, 3)");
        check(moments(rexp(N, 2), 2, 4, "rexp"), "rexp(mean 2)");
        check(moments(rlogis(N, 1, 2), 1, PI * PI * 4 / 3, "rlogis"), "rlogis(1, 2)");
        check(moments(rlaplace(N, 1, 2), 1, 8, "rlaplace"), "rlaplace(1, 2)");
        check(moments(rbeta(N, 2, 5), 2.0 / 7, 10.0 / (49 * 8), "rbeta"), "rbeta(2, 5)");
        auto ged2 = rged(N, 2), ged1 = rged(N, 1);
        check(abs(meanOf(ged2)) < 0.02 && abs(meanOf(ged1)) < 0.02 && varOf(ged2) > 0.5 && varOf(ged2) < 2 && kurtosis(ged1) > kurtosis(ged2) + 0.5 && abs(kurtosis(ged2)) < 0.15, "rged: symmetric, nu = 2 has normal kurtosis (0), nu = 1 heavier tails");
    }
    {
        auto d = rdiscrete(N, [0.2, 0.5, 0.3]);
        size_t[3] cnt;
        foreach (v; d) cnt[v]++;
        check(d.length == N && abs(cnt[0] / cast(double) N - 0.2) < 0.01 && abs(cnt[1] / cast(double) N - 0.5) < 0.01 && abs(cnt[2] / cast(double) N - 0.3) < 0.01, "rdiscrete: values 0 .. 2 with the given frequencies");
        auto dm = rdirichlet([1, 2, 3], 20000);
        bool sums = true;
        foreach (i; 0 .. dm.nrow) sums = sums && abs(dm[i, 0] + dm[i, 1] + dm[i, 2] - 1) < 1e-9;
        auto cm = mean(dm);
        check(dm.nrow == 20000 && dm.ncol == 3 && sums && abs(cm[0, 0] - 1.0 / 6) < 0.01 && abs(cm[0, 1] - 2.0 / 6) < 0.01 && abs(cm[0, 2] - 3.0 / 6) < 0.01, "rdirichlet([1, 2, 3]): rows add to 1, column means alpha / sum");
        auto S = fresh2(2, 2);
        S[0, 0] = 2; S[0, 1] = 0.5; S[1, 0] = 0.5; S[1, 1] = 1;
        auto one = rinvwishart(S, 10);
        check(one.nrow == 2 && one.ncol == 2 && abs(one[0, 1] - one[1, 0]) < 1e-12 && one[0, 0] > 0 && one[0, 0] * one[1, 1] > one[0, 1] ^^ 2, "rinvwishart: a symmetric positive definite matrix");
        auto seq = rinvwishart(S, 10, 20000);
        auto sm = mean(seq);
        check(seq.nrow == 20000 && seq.ncol == 3 && abs(sm[0, 0] - 2.0 / 7) < 0.02 && abs(sm[0, 1] - 0.5 / 7) < 0.01 && abs(sm[0, 2] - 1.0 / 7) < 0.01, "rinvwishart(S, 10, n): vech rows with mean S / (df - p - 1)");
    }

    // ---- invalid parameters ----
    check(throws({ runif(3, 1, 1); }) && throws({ rnorm(3, 0, -1); }) && throws({ rnorm(3, double.nan, 1); }) && throws({ rtnorm(3, 2, 1); }) && throws({ rgamma(3, -1); }) && throws({ rgamma(3, 1, 0); }), "invalid uniform, normal, truncated normal and gamma parameters throw");
    check(throws({ rchisq(3, 0); }) && throws({ rt(3, -1); }) && throws({ rF(3, 0, 5); }) && throws({ rbinom(3, 10, 1.5); }) && throws({ rbinom(3, -1, 0.5); }) && throws({ rbetabinom(3, 10, -1, 2); }), "invalid chi-square, t, F, binomial and beta-binomial parameters throw");
    check(throws({ rpois(3, -1); }) && throws({ rpois([1.0, -2.0]); }) && throws({ rweibull(3, 0); }) && throws({ rexp(3, 0); }) && throws({ rlogis(3, 0, -1); }) && throws({ rged(3, 0); }) && throws({ rlaplace(3, 0, 0); }) && throws({ rbeta(3, 0, 1); }), "invalid Poisson, Weibull, exponential, logistic, GED, Laplace and beta parameters throw");
    check(throws({ rdiscrete(3, [0.5, 0.6]); }) && throws({ rdiscrete(3, [-0.5, 1.5]); }) && throws({ rdiscrete(3, null); }) && throws({ rdirichlet([1.0], 3); }) && throws({ rdirichlet([1.0, -1.0], 3); }) && throws({ rinvwishart(fresh2(2, 3), 5); }) && throws({ rinvwishart(fresh2(2, 2), 1); }), "invalid discrete, Dirichlet and inverse Wishart parameters throw");

    // ---- matrices ----
    {
        auto m = randomMatrix!rnorm(300, 4, 1.0, 0.5);
        check(m.nrow == 300 && m.ncol == 4 && moments(m.ptr.val[0 .. 1200], 1.0, 0.25, "randomMatrix"), "randomMatrix!rnorm(300, 4, 1, 0.5): shape and moments");
        auto g = randomMatrix!rgamma(50, 3, 2.0, 1.0);
        check(g.nrow == 50 && g.ncol == 3 && minimum(g)[0, 0] > 0, "randomMatrix!rgamma");
    }

    // ---- the MRG31k3p streams: against an independent implementation of the recurrences ----
    enum long M1 = 2147483647, M2 = 2147462579;
    {
        auto ss = makeStreams(1);
        auto s = ss[0];
        check(s.state == [12345, 12345, 12345, 12345, 12345, 12345], "the first default stream starts from the original generator's state (all 12345, as prng.d)");
        long[3] h1 = [12345, 12345, 12345], h2 = [12345, 12345, 12345];
        bool same = true;
        foreach (step; 0 .. 5000) {
            const y1 = (4194304 * h1[1] + 129 * h1[2]) % M1;
            const y2 = (32768 * h2[0] + 32769 * h2[2]) % M2;
            h1 = [y1, h1[0], h1[1]];
            h2 = [y2, h2[0], h2[1]];
            const z = (y1 <= y2) ? y1 - y2 + M1 : y1 - y2;
            same = same && s.uniform() == z * 4.656612873077392578125e-10;
        }
        check(same, "5000 draws of a stream equal the two recurrences x1(n) = 2^22 x1(n-2) + 129 x1(n-3) mod 2^31 - 1, x2(n) = 2^15 x2(n-1) + 32769 x2(n-3) mod 2^31 - 21069, combined");
    }
    long[3][3] matMul(long[3][3] a, long[3][3] b, long m) {
        long[3][3] c;
        foreach (i; 0 .. 3) foreach (j; 0 .. 3) { long t = 0; foreach (k; 0 .. 3) t = (t + a[i][k] * b[k][j]) % m; c[i][j] = t; }
        return c;
    }
    long[3][3] A1 = [[0, 4194304, 129], [1, 0, 0], [0, 1, 0]], A2 = [[32768, 0, 32769], [1, 0, 0], [0, 1, 0]];
    foreach (_; 0 .. 134) { A1 = matMul(A1, A1, M1); A2 = matMul(A2, A2, M2); }
    {
        // the jump of 2^134 steps made by squaring 134 times is what separates the streams
        auto ss = makeStreams(4);
        long[3] v1 = [12345, 12345, 12345], v2 = [12345, 12345, 12345];
        bool ok = true;
        foreach (k; 0 .. 4) {
            ok = ok && ss[k].state == [cast(int) v1[0], cast(int) v1[1], cast(int) v1[2], cast(int) v2[0], cast(int) v2[1], cast(int) v2[2]];
            long[3] n1, n2;
            foreach (i; 0 .. 3) { n1[i] = (A1[i][0] * v1[0] + A1[i][1] * v1[1] + A1[i][2] * v1[2]) % M1; n2[i] = (A2[i][0] * v2[0] + A2[i][1] * v2[1] + A2[i][2] * v2[2]) % M2; }
            v1 = n1; v2 = n2;
        }
        check(ok, "stream k starts 2^134 steps after stream k - 1 (the jump matrices are A^(2^134), checked by squaring the transition matrices 134 times)");
        check(ss[0].state != ss[1].state && ss[1].state != ss[2].state && ss[2].state != ss[3].state, "the streams have different states");
    }
    {
        // a seed other than 0 gives a valid, different, reproducible family
        auto a = makeStreams(3, 7), b = makeStreams(3, 7), c = makeStreams(3, 8), d = makeStreams(3);
        bool valid = true;
        foreach (st; a) {
            auto x = st.state;
            valid = valid && x[0] > 0 && x[1] > 0 && x[2] > 0 && x[0] < M1 && x[1] < M1 && x[2] < M1 && x[3] > 0 && x[4] > 0 && x[5] > 0 && x[3] < M2 && x[4] < M2 && x[5] < M2;
        }
        check(valid && a[0].state == b[0].state && a[2].state == b[2].state && a[0].state != c[0].state && a[0].state != d[0].state, "seeds: valid states, reproducible, different seeds give different families");
        check(a.length == 3 && makeStreams(0).length == 0, "makeStreams(n) makes n streams");
    }
    {
        // quality of the uniforms and of the 64-bit words
        auto ss = makeStreams(2, 99);
        enum n = 100_000;
        double[] u0, u1;
        foreach (_; 0 .. n) { u0 ~= ss[0].uniform(); u1 ~= ss[1].uniform(); }
        size_t[10] bins;
        foreach (v; u0) bins[cast(size_t) (v * 10)]++;
        double chi = 0;
        foreach (b; bins) chi += (b - n / 10.0) ^^ 2 / (n / 10.0);
        check(minimum(u0) > 0 && maximum(u0) < 1 && chi < 21.7 && moments(u0, 0.5, 1.0 / 12, "stream uniform"), "stream uniforms: in (0, 1), uniform over ten bins (chi-square below its 1% value), moments");
        check(abs(correlation(u0, u1)) < 4 / sqrt(cast(double) n) && abs(correlation(u0[0 .. n - 1], u0[1 .. n])) < 4 / sqrt(cast(double) n), "no correlation between two streams, nor at lag 1 within one");
        auto c = new RandomStream(makeStreams(1, 5)[0].state);
        double[64] ones = 0;
        foreach (_; 0 .. 20_000) {
            const w = c.next64();
            foreach (b; 0 .. 64) if (w & (1UL << b)) ones[b]++;
        }
        bool bits = true;
        foreach (b; 0 .. 64) bits = bits && abs(ones[b] / 20_000 - 0.5) < 5 * 0.5 / sqrt(20_000.0);
        check(bits, "next64: each of the 64 bits is 1 about half the time (the low and the top two bits too)");
        auto e1 = new RandomStream([1, 2, 3, 4, 5, 6]), e2 = new RandomStream([1, 2, 3, 4, 5, 6]);
        const w3 = e1.next64();
        const a1 = cast(ulong) (e2.uniform() * 2147483648.0), a2 = cast(ulong) (e2.uniform() * 2147483648.0), a3 = cast(ulong) (e2.uniform() * 2147483648.0);
        check(w3 == (a1 | (a2 << 31) | ((a3 & 3) << 62)), "next64 is three steps combined: 31 + 31 + 2 bits");
    }

    // ---- prng.d itself (the pure D generator that the streams are built on) ----
    {
        StreamGen sg;
        auto p0 = sg.createStream(), p1 = sg.createStream(), p2 = sg.createStream();
        auto ms = makeStreams(3);
        check(p0.x11 == ms[0].state[0] && [p1.x11, p1.x12, p1.x13, p1.x21, p1.x22, p1.x23] == ms[1].state && [p2.x11, p2.x12, p2.x13, p2.x21, p2.x22, p2.x23] == ms[2].state,
              "prng.d: StreamGen makes the streams that makeStreams wraps");
        auto q = new Stream([12345, 12345, 12345, 12345, 12345, 12345]);
        check(genDouble(q) == makeStreams(1)[0].uniform() && makeStreams(1)[0].uniform() > 0, "prng.d: genDouble is the stream's U(0, 1) draw");
        size_t[6] faces;
        foreach (_; 0 .. 60_000) faces[genInteger(q, 1, 6) - 1]++;
        bool even = true;
        foreach (f; faces) even = even && abs(f - 10_000.0) < 5 * sqrt(60_000.0 * (1.0 / 6) * (5.0 / 6));
        check(even, "prng.d: genInteger(1, 6) gives the six values equally often");
        auto custom = StreamGen([11, 22, 33, 44, 55, 66]);
        auto c0 = custom.createStream();
        auto c1 = StreamGen([11, 22, 33, 44, 55, 66]).createStream();
        check([c0.x11, c0.x12, c0.x13, c0.x21, c0.x22, c0.x23] == [11, 22, 33, 44, 55, 66] && c0.gen() == c1.gen(), "prng.d: StreamGen(start) begins from the given state, reproducibly");
        Stream asBase = makeStreams(1)[0];
        check(asBase.gen() == makeStreams(1)[0].uniform(), "a RandomStream is a prng.d Stream");
    }

    // ---- streams drive the library's distributions ----
    {
        // a stream's draws equal the module function run with the stream as source, and are reproducible
        auto s = makeStreams(3, 11);
        auto t = makeStreams(3, 11);
        const x = s[1].rnorm(20);
        const y = t[1].rnorm(20);
        double[] z;
        auto c = new RandomStream(makeStreams(3, 11)[1].state);
        c.use({ z = rnorm(20); });
        check(x == y && x == z, "stream.rnorm equals rnorm run under use(), and the same seed gives the same stream");
        check(s[0].rnorm(20) != s[2].rnorm(20), "different streams give different draws");
        // the global generator is not touched
        setSeed(9);
        const g1 = rnorm(3);
        setSeed(9);
        auto u = makeStreams(2, 3);
        u[0].rgamma(50, 2.0);
        u[1].rbinom(50, 10, 0.3);
        const g2 = rnorm(3);
        check(g1 == g2, "draws from streams do not use the global generator: its sequence is the same with them in between");
        // outside use() the global generator is back
        setSeed(9);
        const g3 = rnorm(3);
        check(g3 == g1, "after a stream's draw the module functions use the global generator again");
        // nested use restores the outer source
        auto o1 = makeStreams(2, 21), o2 = makeStreams(2, 21);
        double[] outer;
        o1[0].use({
            auto a1 = rnorm(2);
            auto inner = o1[1].rnorm(2);
            auto a2 = rnorm(2);
            outer = a1 ~ a2;
        });
        auto expect = o2[0].rnorm(4);
        check(outer == expect, "a nested stream restores the outer one afterwards: the outer stream's draws are unbroken");
        // statistics of the distributions driven by a stream
        auto st = makeStreams(1, 77)[0];
        check(moments(st.rnorm(N, 1, 2), 1, 4, "stream rnorm") && moments(st.rgamma(N, 3, 2), 6, 12, "stream rgamma") && moments(st.rbinom(N, 20, 0.3), 6, 4.2, "stream rbinom") && moments(st.rpois(N, 4.5), 4.5, 4.5, "stream rpois"), "stream-driven normal, gamma, binomial and Poisson draws have the right moments");
        check(moments(st.rbeta(N, 2, 5), 2.0 / 7, 10.0 / (49 * 8), "stream rbeta") && moments(st.rt(N, 10), 0, 10.0 / 8, "stream rt") && moments(st.rweibull(N, 2, 3), 3 * gamma(1.5), 9 * (gamma(2.0) - gamma(1.5) ^^ 2), "stream rweibull"), "stream-driven beta, t and Weibull draws");
        auto sm = st.rdirichlet([1.0, 2.0, 3.0], 5000);
        check(sm.nrow == 5000 && abs(mean(sm)[0, 2] - 0.5) < 0.02, "a multivariate draw (Dirichlet) from a stream");
        check(throws({ st.rnorm(3, 0, -1); }), "invalid parameters throw under a stream too");
        check(st.rnorm(2).length == 2, "and the source is released after an error: the stream still works");
    }

    // ---- parallel draws ----
    {
        static double[] work(size_t i) {
            double[] r;
            r ~= rnorm(500, cast(double) i, 1.5);
            r ~= rgamma(300, 2.0 + i, 1.0);
            r ~= rbinom(200, 15, 0.05 + 0.03 * i);
            r ~= rpois(200, 3.0 + i);
            r ~= rtnorm(100, -1, 2);
            r ~= rbeta(100, 2, 3);
            r ~= rt(100, 5);
            auto d = rdirichlet([1.0, 2.0, 3.0], 20);
            foreach (v; d.ptr.val[0 .. 60]) r ~= v;
            foreach (v; rdiscrete(50, [0.1, 0.4, 0.5])) r ~= v;
            return r;
        }
        enum K = 24;
        auto seqStreams = makeStreams(K, 2025);
        double[][] sequential;
        foreach (i; 0 .. K) {
            double[] r;
            seqStreams[i].use({ r = work(i); });
            sequential ~= r;
        }
        auto parStreams = makeStreams(K, 2025);
        auto parallel = parallelDraws(parStreams, (size_t i) => work(i));
        check(parallel == sequential, "parallelDraws: 24 streams over 9 distributions give exactly the draws of the same streams run one after the other");
        auto again = parallelDraws(makeStreams(K, 2025), (size_t i) => work(i));
        check(again == parallel, "and a second parallel run gives the same results (it does not depend on the scheduling; " ~ totalCPUs.to!string ~ " processors)");
        auto withStream = parallelDraws(makeStreams(K, 2025), (RandomStream s, size_t i) => s.rnorm(500, cast(double) i, 1.5));
        check(withStream[5] == sequential[5][0 .. 500] && withStream.length == K, "the form of parallelDraws that passes the stream: the same stream, the same draws");
        bool diff = true;
        foreach (i; 1 .. K) diff = diff && parallel[i][0 .. 50] != parallel[i - 1][0 .. 50];
        check(diff, "different streams, different draws");
        // means of the parallel pieces behave
        double[] means;
        foreach (i; 0 .. K) means ~= meanOf(parallel[i][0 .. 500]);
        bool near = true;
        foreach (i; 0 .. K) near = near && abs(means[i] - i) < 5 * 1.5 / sqrt(500.0);
        check(near, "the normal draws of each task have the mean of that task's parameter");
        auto none = parallelDraws(makeStreams(0), (size_t i) => rnorm(3));
        check(none.length == 0, "no streams: no results");
        // the global generator is not disturbed by a parallel run
        setSeed(31);
        const q1 = rnorm(4);
        setSeed(31);
        parallelDraws(makeStreams(8, 1), (size_t i) => rgamma(100, 2.0));
        check(rnorm(4) == q1, "a parallel run leaves the global generator alone");
    }

    // ---- no growth of the heap ----
    {
        void round() {
            auto s = makeStreams(4, 3);
            auto p = parallelDraws(s, (size_t i) => rgamma(200, 2.0));
            auto d = s[0].rdirichlet([1.0, 1.0, 1.0], 20);
            auto q = rbinom(100, 10, 0.5);
            auto w = rinvwishart(eye(2), 5, 10);
        }
        foreach (_; 0 .. 3) round();
        GC.collect();
        const before = mallinfo().uordblks;
        foreach (_; 0 .. 40) round();
        GC.collect();
        const after = mallinfo().uordblks;
        writefln("  heap in use: %d bytes before, %d after 40 rounds", before, after);
        check(after - before < 40000, "streams, parallel draws and the multivariate generators do not leak");
    }

    writefln("\n%d tests, %d failed", nTests, nFailed);
    if (nFailed) throw new Exception("failed");
}

/* a matrix of zeros: the library's allocator leaves the values uninitialized, and the tests fill in only the nonzero elements */
GretlMatrix fresh2(size_t r, size_t c) {
    import dgretl.internal: fresh;
    auto m = fresh(r, c);
    foreach (i; 0 .. r * c) m.ptr.val[i] = 0;
    return m;
}
