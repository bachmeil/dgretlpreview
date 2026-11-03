/* Tests for dgretl/distributions.d: densities, distribution functions and quantiles against closed forms, std.mathspecial, direct sums, numerical integration and the
   defining properties of a quantile; shapes (scalar, array, matrix), tails, edges, parameter errors, leaks. */

import core.memory: GC;
import std.algorithm: map;
import std.array: array;
import std.conv: to;
import std.math: PI, abs, atan, exp, floor, isNaN, log, sqrt;
import std.mathspecial: betaIncomplete, gammaIncomplete, logGamma, normalDistribution, normalDistributionInverse;
import std.stdio: writefln;
import dgretl;
import dgretl.internal: makeDataset;

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
bool rel(double a, double b, double tol = 1e-9) { return abs(a - b) <= tol * abs(b); }

bool throws(void delegate() f) {
    try f(); catch (GretlException) return true;
    return false;
}

/* Simpson's rule */
double integrate(double delegate(double) f, double a, double b, int n = 20000) {
    const h = (b - a) / n;
    double s = f(a) + f(b);

    foreach (i; 1 .. n) s += f(a + i * h) * ((i % 2) ? 4 : 2);
    return s * h / 3;
}

/* the density integrates to the distribution function, from `from` (below the support or far in the tail) */
bool densityMatchesCdf(double delegate(double) pdf, double delegate(double) cdf, double from, double[] at, double tol = 1e-6) {
    foreach (x; at) {
        const got = cdf(from) + integrate(pdf, from, x);

        if (abs(got - cdf(x)) > tol) {
            writefln("  at %g: integral %g, cdf %g", x, got, cdf(x));
            return false;
        }
    }
    return true;
}

/* quantile undoes the distribution function, in both tails */
bool roundTrip(double delegate(double) cdf, double delegate(double) upper, double delegate(double, bool) quant, double[] ps, double tol = 1e-7) {
    foreach (p; ps) {
        const x = quant(p, true), y = quant(p, false);

        if (!close(cdf(x), p, tol) || !close(upper(y), p, tol)) {
            writefln("  p = %g: cdf(q) = %g, upper(q_upper) = %g", p, cdf(x), upper(y));
            return false;
        }
    }
    return true;
}

immutable double[] probs = [1e-8, 1e-4, 0.01, 0.1, 0.3, 0.5, 0.77, 0.95, 0.999];

void main() {
    // ---- normal ----
    {
        check(rel(pnorm(1.96), normalDistribution(1.96)), "pnorm against std.mathspecial");
        check(close(pnorm(0), 0.5, 1e-15) && close(pnorm(1.0, 1.0, 2.0), 0.5, 1e-15), "pnorm at the mean is 1/2; location and scale");
        check(rel(pnorm(10.0, 0, 1, false), 7.619853024160527e-24, 1e-9), "upper tail of the normal in the far tail (7.6e-24)");
        check(rel(pnorm(-10.0), 7.619853024160527e-24, 1e-9), "lower tail of the normal in the far tail");
        check(close(dnorm(0), 1 / sqrt(2 * PI), 1e-14) && close(dnorm(3.0, 1, 2), exp(-0.5) / (2 * sqrt(2 * PI)), 1e-14), "dnorm");
        check(close(qnorm(0.975), 1.959963984540054, 1e-9) && close(qnorm(0.975, 10.0, 3.0), 10 + 3 * 1.959963984540054, 1e-9), "qnorm");
        check(close(qnorm(0.025, lower: false), 1.959963984540054, 1e-9), "qnorm upper tail");
        check(roundTrip(x => pnorm(x, 2, 3), x => pnorm(x, 2, 3, false), (p, lw) => qnorm(p, 2, 3, lw), probs.dup), "normal quantile round trip, both tails");
        check(pnorm(double.infinity) == 1 && pnorm(-double.infinity) == 0 && pnorm(double.infinity, 0, 1, false) == 0 && dnorm(double.infinity) == 0, "infinite arguments");
        check(qnorm(0.0) == -double.infinity && qnorm(1.0) == double.infinity && qnorm(0.0, lower: false) == double.infinity, "quantiles at 0 and 1");
        check(densityMatchesCdf(x => dnorm(x), x => pnorm(x), -12, [-2, 0, 1.5, 4]), "dnorm integrates to pnorm");
    }

    // ---- t ----
    {
        check(close(pt(2.0, 1), 0.5 + atan(2.0) / PI, 1e-12), "pt(df = 1) is the Cauchy distribution");
        check(close(pt(2.0, 5), 0.949030260585, 1e-9), "pt(2, 5) against a known value");
        check(close(qt(0.975, 5), 2.570581835636, 1e-9) && close(qt(0.975, 5.0), -qt(0.025, 5.0), 1e-12), "qt(.975, 5); symmetry");
        check(close(pt(-3.0, 7) + pt(3.0, 7), 1, 1e-12) && close(pt(3.0, 7, lower: false), pt(-3.0, 7), 1e-12), "pt symmetry and upper tail");
        check(close(dt(0, 1), 1 / PI, 1e-12) && close(dt(1.5, 3), 2 / (PI * sqrt(3.0)) / (1.75 * 1.75), 1e-12), "dt(0, 1) is 1/pi; dt(1.5, 3)");
        check(roundTrip(x => pt(x, 6), x => pt(x, 6, 0, false), (p, lw) => qt(p, 6, 0, lw), probs.dup), "t quantile round trip, both tails");
        check(densityMatchesCdf(x => dt(x, 4), x => pt(x, 4), -400, [-3, 0, 2]), "dt integrates to pt");
        // non-central: ncp = 0 is the central distribution; the density integrates to the distribution function; the quantile undoes it
        check(close(pt(1.2, 8, ncp: 0), pt(1.2, 8), 1e-15), "ncp = 0 is central");
        check(densityMatchesCdf(x => dt(x, 5, 1.0), x => pt(x, 5, 1.0), -60, [-1, 1, 3], 1e-4), "non-central t: density integrates to the distribution function");
        check(roundTrip(x => pt(x, 5, 1.0), x => pt(x, 5, 1.0, false), (p, lw) => qt(p, 5, 1.0, lw), [0.05, 0.5, 0.9], 1e-6), "non-central t quantile round trip");
    }

    // ---- chi-square and gamma ----
    {
        check(close(pchisq(3.0, 2), 1 - exp(-1.5), 1e-12) && close(dchisq(3.0, 2), 0.5 * exp(-1.5), 1e-12), "chi-square(2) is exponential");
        check(close(qchisq(0.95, 3), 7.814727903251, 1e-9) && close(qchisq(0.05, 3, lower: false), 7.814727903251, 1e-9), "qchisq(.95, 3)");
        check(rel(pchisq(2.7, 2.5), gammaIncomplete(1.25, 1.35), 1e-9), "chi-square with a non-integer df against the incomplete gamma");
        check(close(qchisq(0.9, 4.5), qgamma(0.9, 2.25, 2), 1e-9), "non-integer chi-square is gamma(df/2, 2)");
        check(roundTrip(x => pchisq(x, 4), x => pchisq(x, 4, 0, false), (p, lw) => qchisq(p, 4, 0, lw), probs.dup), "chi-square round trip");
        check(pchisq(-1.0, 3) == 0 && pchisq(-1.0, 3, 0, false) == 1 && dchisq(-1.0, 3) == 0, "below the support");
        check(rel(pchisq(1e-3, 3, 0, false), 1 - pchisq(1e-3, 3), 1e-12), "chi-square tails add up");
        check(densityMatchesCdf(x => dchisq(x, 4, 2.0), x => pchisq(x, 4, 2.0), 0, [1, 4, 9], 1e-5), "non-central chi-square density integrates");
        check(roundTrip(x => pchisq(x, 4, 2.0), x => pchisq(x, 4, 2.0, false), (p, lw) => qchisq(p, 4, 2.0, lw), [0.1, 0.5, 0.9], 1e-6), "non-central chi-square round trip");
        check(rel(pgamma(4.0, 2, 3), gammaIncomplete(2, 4.0 / 3), 1e-12), "pgamma against the incomplete gamma");
        check(close(dgamma(4.0, 2, 3), 4.0 / 9 * exp(-4.0 / 3), 1e-12), "dgamma(shape 2, scale 3)");
        check(roundTrip(x => pgamma(x, 2.5, 1.5), x => pgamma(x, 2.5, 1.5, false), (p, lw) => qgamma(p, 2.5, 1.5, lw), probs.dup), "gamma round trip");
    }

    // ---- F and beta ----
    {
        const x = 2.0;
        check(rel(pF(x, 4, 9), betaIncomplete(2.0, 4.5, 4 * x / (4 * x + 9)), 1e-9), "pF against the incomplete beta");
        check(close(qF(0.95, 4, 9), 3.633088568, 1e-7), "qF(.95, 4, 9)");
        check(close(qF(0.95, 1, 7), qt(0.975, 7) ^^ 2, 1e-9), "F(1, n) is t squared");
        check(roundTrip(a => pF(a, 5, 12), a => pF(a, 5, 12, 0, false), (p, lw) => qF(p, 5, 12, 0, lw), probs.dup), "F round trip");
        check(densityMatchesCdf(a => dF(a, 5, 12), a => pF(a, 5, 12), 0, [0.5, 1, 3]), "dF integrates");
        check(throws({ pF(1.0, 2.5, 4); }) && throws({ pF(1.0, 2, 0); }), "F rejects fractional or zero degrees of freedom");
        check(densityMatchesCdf(a => dF(a, 4, 9, 2.0), a => pF(a, 4, 9, 2.0), 0, [1, 3], 1e-5), "non-central F density integrates");
        check(rel(pbeta(0.3, 2, 5), betaIncomplete(2, 5, 0.3), 1e-12) && close(dbeta(0.3, 2, 5), 30 * 0.3 * 0.7 ^^ 4, 1e-12), "beta");
        check(close(pbeta(0.3, 2, 5, lower: false), 1 - pbeta(0.3, 2, 5), 1e-12), "beta upper tail");
        check(roundTrip(a => pbeta(a, 2.5, 0.7), a => pbeta(a, 2.5, 0.7, false), (p, lw) => qbeta(p, 2.5, 0.7, lw), probs.dup), "beta round trip");
        check(pbeta(-0.5, 2, 3) == 0 && pbeta(1.5, 2, 3) == 1 && dbeta(1.5, 2, 3) == 0, "outside [0, 1]");
    }

    // ---- exponential, Weibull, logistic, uniform, Laplace, GED ----
    {
        check(close(pexp(1.0, 2), 1 - exp(-2.0), 1e-14) && close(dexp(1.0, 2), 2 * exp(-2.0), 1e-14) && close(qexp(0.4, 2), -log(0.6) / 2, 1e-14), "exponential with a rate");
        check(rel(pexp(50.0, 1, false), exp(-50.0), 1e-12), "exponential far tail");
        check(close(pweibull(2.0, 2, 3), 1 - exp(-(2.0 / 3) ^^ 2), 1e-14) && close(dweibull(2.0, 2, 3), 2.0 / 3 * (2.0 / 3) * exp(-(2.0 / 3) ^^ 2), 1e-14), "Weibull(shape 2, scale 3)");
        check(roundTrip(a => pweibull(a, 1.7, 2), a => pweibull(a, 1.7, 2, false), (p, lw) => qweibull(p, 1.7, 2, lw), probs.dup, 1e-9), "Weibull round trip");
        check(close(plogis(1.0), 1 / (1 + exp(-1.0)), 1e-14) && close(dlogis(0), 0.25, 1e-14) && close(qlogis(0.8), log(4.0), 1e-12), "logistic");
        check(roundTrip(a => plogis(a, 1, 2), a => plogis(a, 1, 2, false), (p, lw) => qlogis(p, 1, 2, lw), probs.dup), "logistic round trip");
        check(close(punif(0.25), 0.25, 1e-15) && close(dunif(2.0, 1, 5), 0.25, 1e-15) && close(qunif(0.5, 1, 5), 3, 1e-15) && punif(7.0, 1, 5) == 1 && dunif(7.0, 1, 5) == 0, "uniform");
        check(close(qunif(0.2, 1, 5, lower: false), 4.2, 1e-14), "uniform quantile, upper tail");
        check(close(plaplace(2.0, 1, 2), 1 - 0.5 * exp(-0.5), 1e-12) && close(dlaplace(2.0, 1, 2), 0.25 * exp(-0.5), 1e-12), "Laplace");
        check(roundTrip(a => plaplace(a, 1, 2), a => plaplace(a, 1, 2, false), (p, lw) => qlaplace(p, 1, 2, lw), probs.dup), "Laplace round trip");
        check(densityMatchesCdf(a => dged(a, 1.5), a => pged(a, 1.5), -30, [-1, 0.5, 2], 1e-5), "GED density integrates");
        check(roundTrip(a => pged(a, 1.5), a => pged(a, 1.5, false), (p, lw) => qged(p, 1.5, lw), [0.01, 0.3, 0.5, 0.9], 1e-6), "GED round trip");
    }

    // ---- truncated normal ----
    {
        const a = -1.0, b = 2.0;
        const mass = normalDistribution(b) - normalDistribution(a);
        check(close(ptnorm(1.0, 0, 1, a, b), (normalDistribution(1.0) - normalDistribution(a)) / mass, 1e-12), "truncated normal distribution function");
        check(close(dtnorm(1.0, 0, 1, a, b), exp(-0.5) / sqrt(2 * PI) / mass, 1e-12), "truncated normal density");
        check(ptnorm(-2.0, 0, 1, a, b) == 0 && ptnorm(3.0, 0, 1, a, b) == 1 && dtnorm(3.0, 0, 1, a, b) == 0, "outside the truncation limits");
        check(roundTrip(x => ptnorm(x, 0.5, 2, a, b), x => ptnorm(x, 0.5, 2, a, b, false), (p, lw) => qtnorm(p, 0.5, 2, a, b, lw), probs.dup), "truncated normal round trip");
        check(qtnorm(0.0, 0, 1, a, b) == a && qtnorm(1.0, 0, 1, a, b) == b, "quantiles at the limits");
        // one-sided truncation and an interval far in the right tail
        check(close(ptnorm(0.5, 0, 1, 0), 2 * normalDistribution(0.5) - 1, 1e-12), "half normal");
        check(densityMatchesCdf(x => dtnorm(x, 0, 1, 5, 6), x => ptnorm(x, 0, 1, 5, 6), 5, [5.2, 5.5, 5.9], 1e-8), "truncated normal on [5, 6] integrates");
        check(roundTrip(x => ptnorm(x, 0, 1, 5, 6), x => ptnorm(x, 0, 1, 5, 6, false), (p, lw) => qtnorm(p, 0, 1, 5, 6, lw), [0.01, 0.5, 0.99], 1e-8), "truncated normal on [5, 6] round trip");
        check(close(ptnorm(0.3, 0, 1), pnorm(0.3), 1e-12), "no truncation is the normal");
    }

    // ---- binomial and Poisson ----
    {
        const n = 10; const pr = 0.3;
        double mass(int k) { return exp(logGamma(n + 1.0) - logGamma(k + 1.0) - logGamma(n - k + 1.0) + k * log(pr) + (n - k) * log(1 - pr)); }
        double cum(int k) { double s = 0; foreach (j; 0 .. k + 1) s += mass(j); return s; }
        bool all = true, allP = true, allQ = true, allU = true;
        foreach (k; 0 .. n + 1) {
            all = all && close(dbinom(k, n, pr), mass(k), 1e-12);
            allP = allP && close(pbinom(k, n, pr), cum(k), 1e-10);
            allU = allU && close(pbinom(k, n, pr, false), 1 - cum(k), 1e-10);
        }
        check(all, "binomial mass against the closed form");
        check(allP && allU, "binomial distribution function and upper tail against the sums");
        foreach (p; [0.01, 0.2, 0.5, 0.7, 0.95, 0.999]) {
            const q = cast(int) qbinom(p, n, pr);
            // the smallest k with cdf(k) >= p
            allQ = allQ && cum(q) >= p - 1e-12 && (q == 0 || cum(q - 1) < p);
            const qu = cast(int) qbinom(p, n, pr, false);
            allQ = allQ && pbinom(qu, n, pr, false) <= p + 1e-12 && (qu == 0 || pbinom(qu - 1, n, pr, false) > p);
        }
        check(allQ, "binomial quantile is the smallest k with cdf >= p (and the upper-tail form)");
        check(dbinom(2.5, n, pr) == 0 && close(pbinom(2.5, n, pr), pbinom(2.0, n, pr), 1e-15) && pbinom(-1.0, n, pr) == 0 && pbinom(10.0, n, pr) == 1 && pbinom(50.0, n, pr) == 1, "non-integer and out-of-range arguments");
        check(qbinom(0.0, n, pr) == 0 && qbinom(1.0, n, pr) == n, "binomial quantiles at 0 and 1");

        const lam = 4.0;
        double pm(int k) { return exp(-lam + k * log(lam) - logGamma(k + 1.0)); }
        double pc(int k) { double s = 0; foreach (j; 0 .. k + 1) s += pm(j); return s; }
        bool pAll = true, pQ = true;
        foreach (k; 0 .. 20) {
            pAll = pAll && close(dpois(k, lam), pm(k), 1e-12) && close(ppois(k, lam), pc(k), 1e-10) && close(ppois(k, lam, false), 1 - pc(k), 1e-10);
        }
        check(pAll, "Poisson mass, distribution function and upper tail against the sums");
        foreach (p; [0.01, 0.2, 0.5, 0.7, 0.95, 0.999]) {
            const q = cast(int) qpois(p, lam);
            pQ = pQ && pc(q) >= p - 1e-12 && (q == 0 || pc(q - 1) < p);
        }
        check(pQ, "Poisson quantile is the smallest k with cdf >= p");
        check(qpois(1.0, lam) == double.infinity && qpois(0.0, lam) == 0, "Poisson quantiles at 0 and 1");
        check(ppois(double.infinity, lam) == 1 && dpois(-1.0, lam) == 0, "Poisson edges");
    }

    // ---- shapes, missing values, errors ----
    {
        double[] x = [-1.0, 0, 1.5, double.nan];
        auto r = pnorm(x);
        check(is(typeof(r) == double[]) && r.length == 4 && close(r[0], pnorm(-1.0), 1e-15) && close(r[2], pnorm(1.5), 1e-15) && isNaN(r[3]), "arrays give arrays; NaN stays NaN");
        auto qi = pnorm([1, 2, 3]);
        check(qi.length == 3 && close(qi[1], pnorm(2.0), 1e-15), "an int array is accepted");
        check(close(pnorm(1), pnorm(1.0), 1e-15), "an int argument is accepted");

        auto m = cbind([0.1, 0.4], [0.2, 0.5], [0.3, 0.6]);
        auto qm = qnorm(m);
        bool same = qm.nrow == m.nrow && qm.ncol == m.ncol;
        foreach (j; 0 .. m.ncol) foreach (i; 0 .. m.nrow) same = same && close(qm[i, j], qnorm(m[i, j]), 1e-15);
        check(same, "a matrix gives a matrix of the same shape, element by element");
        auto sub = pt(m.cols[0 .. 2], 5);
        check(sub.nrow == 2 && sub.ncol == 2 && close(sub[1, 1], pt(0.5, 5), 1e-15), "a submatrix is accepted");

        check(isNaN(qnorm(-0.1)) && isNaN(qnorm(1.1)) && isNaN(qchisq(double.nan, 3)) && isNaN(qbinom(2.0, 10, 0.3)), "probabilities outside [0, 1] and missing ones give NaN");
        check(close(dnorm(x: 0.0, mean: 0.0, sd: 1.0), 1 / sqrt(2 * PI), 1e-14) && close(pt(2.0, df: 8, ncp: 0.0, lower: false), pt(2.0, 8, 0, false), 1e-15), "named arguments");

        check(throws({ dnorm(0, 0, 0); }) && throws({ pnorm(0, 0, -1); }) && throws({ pnorm(0, double.nan, 1); }), "normal parameter errors");
        check(throws({ pt(1.0, 0); }) && throws({ pt(1.0, -2); }) && throws({ pt(1.0, double.infinity); }), "t parameter errors");
        check(throws({ pchisq(1.0, 0); }) && throws({ pchisq(1.0, 3, -1); }) && throws({ pgamma(1.0, 0); }) && throws({ pgamma(1.0, 1, -2); }), "chi-square and gamma parameter errors");
        check(throws({ pbeta(0.5, 0, 1); }) && throws({ pexp(1.0, 0); }) && throws({ pweibull(1.0, -1); }) && throws({ plogis(0.0, 0, 0); }), "beta, exponential, Weibull, logistic parameter errors");
        check(throws({ punif(0.5, 1, 1); }) && throws({ plaplace(0.0, 0, 0); }) && throws({ pged(0.0, 0); }), "uniform, Laplace, GED parameter errors");
        check(throws({ ptnorm(0.0, 0, 1, 2, 1); }) && throws({ ptnorm(0.0, 0, 0); }), "truncated normal parameter errors");
        check(throws({ pbinom(1.0, -1, 0.5); }) && throws({ pbinom(1.0, 5, 1.5); }) && throws({ ppois(1.0, 0); }) && throws({ ppois(1.0, -1); }), "binomial and Poisson parameter errors");
    }

    // ---- the series of a Dataset ----
    {
        double[] z = [-1.5, 0, 0.7, 2.2, double.nan];
        double[] u = [0.1, 0.5, 0.9, 0.99, 0.2];
        auto ds = makeDataset([z, u], ["const", "z", "u"]);
        auto a = pnorm(ds, "z", 1, 2), b = dnorm(ds, "z"), c = qt(ds, "u", 8, lower: false), e = pF(ds, "z", 3, 7, 1.5);
        bool ok = a.length == 5 && b.length == 5 && c.length == 5 && e.length == 5;
        foreach (i; 0 .. 5) {
            ok = ok && (isNaN(z[i]) ? isNaN(a[i]) && isNaN(b[i]) && isNaN(e[i]) : close(a[i], pnorm(z[i], 1, 2), 1e-15) && close(b[i], dnorm(z[i]), 1e-15) && close(e[i], pF(z[i], 3, 7, 1.5), 1e-15));
            ok = ok && close(c[i], qt(u[i], 8, 0, false), 1e-15);
        }
        check(ok, "a series of a Dataset: d, p and q forms, with parameters, ncp, lower and a missing value");
        check(close(pbinom(ds, "u", 5, 0.3)[1], pbinom(0.5, 5, 0.3), 1e-15), "Dataset form of the binomial (integer size)");
        check(throws({ pnorm(ds, "nosuchseries"); }) && throws({ pnorm(ds, "z", 0, -1); }), "unknown series and invalid parameters throw");
    }

    // ---- parameters that vary by element ----
    {
        double[] x = [-1.0, 0.3, 1.2, 2.5];
        double[] mu = [0.0, 0.5, 1.0, 1.5], sg = [1.0, 2.0, 0.5, 3.0];
        auto a = pnorm(x, mu, sg), b = dnorm(x, mu, sg), c = qnorm([0.1, 0.5, 0.9, 0.99], mu, sg, false);
        bool ok = a.length == 4 && b.length == 4 && c.length == 4;
        foreach (i; 0 .. 4) {
            ok = ok && close(a[i], pnorm(x[i], mu[i], sg[i]), 1e-15) && close(b[i], dnorm(x[i], mu[i], sg[i]), 1e-15)
                && close(c[i], qnorm([0.1, 0.5, 0.9, 0.99][i], mu[i], sg[i], false), 1e-15);
        }
        check(ok, "array parameters: element by element (mean and sd, d, p and q, lower: false)");

        auto s1 = pnorm(x, 1.0, sg), s2 = pnorm(x, mu, 2);
        check(close(s1[1], pnorm(0.3, 1, 2), 1e-15) && close(s2[3], pnorm(2.5, 1.5, 2), 1e-15), "a number and an array parameter: the number is repeated");
        auto sx = pnorm(0.4, mu, sg);
        check(sx.length == 4 && close(sx[2], pnorm(0.4, 1.0, 0.5), 1e-15), "a number with array parameters gives an array");

        double[] dfs = [3, 5, 8, 30];
        auto tq = qt(0.975, dfs), tp = pt(x, dfs, ncp: 0.0);
        check(close(tq[0], qt(0.975, 3), 1e-15) && close(tq[3], qt(0.975, 30), 1e-15) && close(tp[2], pt(1.2, 8), 1e-15), "the first parameter (df) as an array; ncp keeps its default");
        double[] ncps = [0, 0.5, 1, 2];
        auto nc = pt(x, 8, ncps);
        check(close(nc[0], pt(-1.0, 8), 1e-15) && close(nc[3], pt(2.5, 8, 2), 1e-15), "the second parameter (ncp) as an array, df a number");
        auto ln = pnorm(x, 0, 1, lower: false);
        check(close(ln[1], pnorm(0.3, 0, 1, false), 1e-15), "numbers still work (and lower by name)");
        auto bs = pbinom([2.0, 3, 4], [5, 10, 20], 0.3);
        check(close(bs[1], pbinom(3.0, 10, 0.3), 1e-15) && close(bs[2], pbinom(4.0, 20, 0.3), 1e-15), "binomial size as an array");
        auto tn = ptnorm(x, 0, 1, [-2.0, -1, 0, 1], 3);
        check(close(tn[1], ptnorm(0.3, 0, 1, -1, 3), 1e-15), "truncated normal with an array of lower limits");

        // matrices
        auto m = cbind([0.1, 0.4], [0.2, 0.5]);
        auto mm = cbind([0.0, 1.0], [2.0, 3.0]);
        auto r = qnorm(m, mm, 2);
        check(r.nrow == 2 && r.ncol == 2 && close(r[1, 0], qnorm(0.4, 1, 2), 1e-15) && close(r[0, 1], qnorm(0.2, 2, 2), 1e-15), "a matrix parameter for a matrix argument (same shape)");
        auto r2 = qnorm(m, [0.0, 1.0, 2.0, 3.0], 2);
        check(close(r2[1, 0], qnorm(0.4, 1, 2), 1e-15) && close(r2[1, 1], qnorm(0.5, 3, 2), 1e-15), "an array parameter for a matrix argument (column by column)");

        // missing values and errors
        double[] mu2 = [0.0, double.nan, 1.0, 1.5];
        auto nn = pnorm(x, mu2, 1);
        check(isNaN(nn[1]) && close(nn[2], pnorm(1.2, 1, 1), 1e-15), "a missing value in an array parameter gives NaN for that element");
        check(throws({ pnorm(x, [0.0, 1], 1); }) && throws({ pnorm([1.0, 2], mu, 1); }), "arrays of different lengths throw");
        check(throws({ pnorm(x, 0, [1.0, 2, -1, 1]); }) && throws({ pt(x, [3.0, 0, 5, 5]); }) && throws({ pbinom(1.0, [5, 2.5], 0.5); }), "an invalid element throws");
        check(throws({ pnorm(x, double.nan, 1); }), "a missing single parameter still throws");
        check(pnorm(new double[0], new double[0], 1).length == 0, "empty arrays");

        // Dataset: parameters as series names
        auto ds = makeDataset([x, mu, sg], ["const", "y", "m", "s"]);
        auto dn = pnorm(ds, "y", "m", "s");
        check(close(dn[1], pnorm(0.3, 0.5, 2), 1e-15) && close(dn[3], pnorm(2.5, 1.5, 3), 1e-15), "Dataset: parameters named by series");
        auto dm = pnorm(ds, "y", "m", 1);
        auto dd = pt(ds, "y", [3.0, 5, 8, 30]);
        auto dq = qnorm(ds, "y", 0, "s", false);
        check(close(dm[2], pnorm(1.2, 1, 1), 1e-15) && close(dd[3], pt(2.5, 30), 1e-15) && close(dq[1], qnorm(0.3, 0, 2, false), 1e-15), "Dataset: a series name, an array and a number mixed");
        check(throws({ pnorm(ds, "y", "nosuch", 1); }), "Dataset: an unknown series for a parameter throws");
    }

    // ---- array parameters by name, skipping the ones before ----
    {
        double[] x = [-1.0, 0.3, 1.2, 2.5];
        double[] sg = [1.0, 2.0, 0.5, 3.0], mu = [0.0, 0.5, 1.0, 1.5];
        auto a = pnorm(x, sd: sg);
        auto b = dnorm(x, sd: sg, mean: 1.0);
        auto c = qnorm([0.1, 0.5, 0.9, 0.99], mean: mu, lower: false);
        bool ok = true;
        foreach (i; 0 .. 4) {
            ok = ok && close(a[i], pnorm(x[i], 0, sg[i]), 1e-15) && close(b[i], dnorm(x[i], 1, sg[i]), 1e-15) && close(c[i], qnorm([0.1, 0.5, 0.9, 0.99][i], mu[i], 1, false), 1e-15);
        }
        check(ok, "pnorm(x, sd: array) skips the mean; named in any order; mean: array with lower: false");
        auto n2 = ptnorm(x, hi: [3.0, 3, 4, 4]);
        auto n3 = ptnorm(x, 0, 2, lo: [-3.0, -2, -1, 0]);
        check(close(n2[1], ptnorm(0.3, 0, 1, -double.infinity, 3), 1e-15) && close(n3[3], ptnorm(2.5, 0, 2, 0, double.infinity), 1e-15), "truncated normal: a later limit as an array by name");
        auto g = pgamma(x.map!(v => v * v + 0.1).array, 2, scale: sg);
        check(close(g[1], pgamma(0.19, 2, 2), 1e-15), "gamma: scale as an array by name, shape a number");
        auto ncs = pt(x, 8, ncp: [0.0, 0.5, 1, 2]);
        auto ncs2 = pt(x, df: 8, ncp: [0.0, 0.5, 1, 2], lower: false);
        check(close(ncs[3], pt(2.5, 8, 2), 1e-15) && close(ncs2[2], pt(1.2, 8, 1, false), 1e-15), "ncp as an array by name");
        auto lo = punif(x, max: [1.0, 2, 3, 4]);
        check(close(lo[1], punif(0.3, 0, 2), 1e-15), "uniform: max as an array, min the default");
        auto bn = pbinom([1.0, 2, 3], 10, prob: [0.1, 0.2, 0.3]);
        check(close(bn[2], pbinom(3.0, 10, 0.3), 1e-15), "binomial: prob as an array by name");
        auto bsz = pbinom([1.0, 2, 3], size: [5, 6, 7], prob: 0.4);
        check(close(bsz[1], pbinom(2.0, 6, 0.4), 1e-15), "binomial: size as an array by name");
        // scalar calls are untouched
        check(close(pnorm(1.0, sd: 2), pnorm(1.0, 0, 2), 1e-15) && close(pt(1.0, 5, lower: false), pt(1.0, 5, 0, false), 1e-15), "named numbers still work");
        // the Dataset forms
        auto ds = makeDataset([x, sg], ["const", "y", "s"]);
        auto dd = pnorm(ds, "y", sd: "s");
        auto de = pnorm(ds, "y", sd: [1.0, 2, 0.5, 3], mean: 1.0);
        check(close(dd[3], pnorm(2.5, 0, 3), 1e-15) && close(de[2], pnorm(1.2, 1, 0.5), 1e-15), "Dataset: a later parameter by name, as a series name or an array");
        // matrices still work as parameters
        auto m = cbind([0.1, 0.4], [0.2, 0.5]);
        auto mp = cbind([1.0, 2.0], [3.0, 4.0]);
        auto mr = qnorm(m, 0, mp);
        check(close(mr[1, 1], qnorm(0.5, 0, 4), 1e-15) && close(mr[0, 0], qnorm(0.1, 0, 1), 1e-15), "a matrix parameter, in order");
        check(throws({ pnorm(x, sd: [1.0, 2]); }) && throws({ pnorm(x, sd: [1.0, 2, -1, 1]); }), "named array parameters: wrong length and invalid element throw");
    }

    // ---- matrix and submatrix parameters, by position and by name ----
    {
        auto m = cbind([0.1, 0.4], [0.2, 0.5]);
        auto mu = cbind([0.0, 1.0], [2.0, 3.0]);
        auto sd = cbind([1.0, 2.0], [3.0, 4.0]);
        auto byName = qnorm(m, sd: sd);
        auto both = qnorm(m, mean: mu, sd: sd);
        auto swapped = qnorm(m, sd: sd, mean: mu);
        bool ok = byName.nrow == 2 && byName.ncol == 2;
        foreach (j; 0 .. 2) foreach (i; 0 .. 2) {
            ok = ok && close(byName[i, j], qnorm(m[i, j], 0, sd[i, j]), 1e-15) && close(both[i, j], qnorm(m[i, j], mu[i, j], sd[i, j]), 1e-15) && close(swapped[i, j], both[i, j], 1e-15);
        }
        check(ok, "matrix parameters by name: sd alone, mean and sd, in either order");
        auto lowerM = pnorm(m, sd: sd, lower: false);
        check(close(lowerM[1, 1], pnorm(0.5, 0, 4, false), 1e-15), "a matrix parameter by name with lower: false");
        auto wide = cbind([1.0, 2.0], [3.0, 4.0], [5.0, 6.0]);
        auto sub = pnorm(m, 0, wide.cols[1 .. 3]);
        check(close(sub[0, 0], pnorm(0.1, 0, 3), 1e-15) && close(sub[1, 1], pnorm(0.5, 0, 6), 1e-15), "a Submatrix parameter, by position");
        auto subN = pnorm(m, sd: wide.cols[0 .. 2]);
        check(close(subN[1, 0], pnorm(0.4, 0, 2), 1e-15) && close(subN[0, 1], pnorm(0.2, 0, 3), 1e-15), "a Submatrix parameter, by name");
        auto mixed = pnorm(m, mean: mu, sd: [1.0, 2, 3, 4]);
        auto mixed2 = ptnorm(m, hi: wide.cols[0 .. 2], lo: -1);
        check(close(mixed[1, 1], pnorm(0.5, 3, 4), 1e-15) && close(mixed2[1, 0], ptnorm(0.4, 0, 1, -1, 2), 1e-15), "matrices, arrays and numbers mixed, by name");
        auto vecx = pnorm([0.1, 0.4, 0.2, 0.5], sd: sd);
        check(vecx.length == 4 && close(vecx[2], pnorm(0.2, 0, 3), 1e-15), "an array argument with a matrix parameter (column by column)");
        check(throws({ pnorm(m, sd: wide); }) && throws({ pnorm(m, sd: cbind([1.0, -2.0], [3.0, 4.0])); }), "a matrix parameter of the wrong size, or with an invalid element, throws");
    }

    // ---- Dataset forms with matrix and submatrix parameters ----
    {
        auto y = [0.1, 0.4, 0.2, 0.5];
        auto ds = makeDataset([y], ["const", "y"]);
        auto sdM = cbind([1.0, 2.0], [3.0, 4.0]);
        auto wide = cbind([1.0, 2.0], [3.0, 4.0], [5.0, 6.0]);
        auto a = pnorm(ds, "y", sd: sdM);
        auto b = pnorm(ds, "y", 0, wide.cols[0 .. 2]);
        auto c = qnorm(ds, "y", mean: sdM, sd: [1.0, 2, 3, 4], lower: false);
        auto d = ptnorm(ds, "y", hi: wide.cols[1 .. 3], lo: -1);
        check(a.length == 4 && close(a[2], pnorm(0.2, 0, 3), 1e-15) && close(b[3], pnorm(0.5, 0, 4), 1e-15), "Dataset: matrix and submatrix parameters");
        check(close(c[1], qnorm(0.4, 2, 2, false), 1e-15) && close(d[3], ptnorm(0.5, 0, 1, -1, 6), 1e-15), "Dataset: matrix, array and number mixed, by name");
        check(throws({ pnorm(ds, "y", sd: wide); }), "Dataset: a matrix parameter of the wrong size throws");
    }

    // ---- no leaks ----
    {
        auto m = cbind([0.1, 0.2], [0.3, 0.4]);

        foreach (i; 0 .. 200) { auto q = qt(m, 7); auto d = dgamma(m, 2, 1.5); auto r = pchisq([1.0, 2, 3], 4, 1.5); }
        GC.collect();
        const before = mallinfo().uordblks;
        foreach (i; 0 .. 2000) {
            auto q = qt(m, 7);
            auto d = dgamma(m, 2, 1.5);
            auto r = pchisq([1.0, 2, 3], 4, 1.5);
            auto b = pbinom(3.0, 10, 0.3);
            auto v = pnorm(m, m, [1.0, 2, 3, 4]);
        }
        GC.collect();
        const grown = mallinfo().uordblks - before;
        check(grown < 200_000, "no leaks over 2000 rounds (" ~ to!string(grown) ~ " bytes)");
    }

    writefln("%d checks, %d failed", nTests, nFailed);
    if (nFailed) {
        import core.stdc.stdlib: exit;
        exit(1);
    }
}
