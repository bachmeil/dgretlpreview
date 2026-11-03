/* Probability distributions: density (or mass), distribution function and quantile function of the distributions of the library (the generators of dgretl.random
   draw from the same family), with the names and the argument order of R.

   import dgretl;
   double p = pnorm(1.96);                      // 0.975
   double c = qchisq(0.95, 3);                  // 95th percentile of chi-square(3)
   double u = pt(2.5, 10, lower: false);        // upper tail P(T > 2.5): accurate in the far tail
   auto d = dnorm(x, mean: 1, sd: 2);           // x: double[] (or a GretlMatrix): the same shape back
   double pw = pt(2.0, 8, ncp: 1.5);            // non-central t

   For each distribution X there are three functions: dX (density; the probability mass of a discrete distribution), pX (P(X <= q), or P(X > q) with lower: false) and
   qX (the smallest x with P(X <= x) >= p for a discrete distribution, the x with P(X <= x) = p for a continuous one; with lower: false, p is the upper-tail probability).
   The first argument is one number, a numeric array (the result is a double[]), or a GretlMatrix or Submatrix (the result is a matrix of the same shape), or a Dataset and the
   name of one of its series (the result is a double[], one value per observation: pnorm(d, "z"), qt(d, "p", 8)); the parameters
   are single numbers (or arrays, matrices, series names: see `distribution` below). Missing values (NaN) in the first argument give NaN. A probability outside [0, 1] gives NaN. Invalid parameters (a negative standard deviation,

   The parameters can vary by element: any of them can be an array (or a matrix) with one value per element of the first argument, and a number is repeated:
       pnorm(x, mu, sigma)            // x, mu, sigma: double[] of the same length
       qnorm(0.975, mu, sd)           // a number with array parameters gives an array (one value per element of the parameters)
       pt(x, dfArray);  pt(x, 8, ncpArray);  pnorm(d, "y", "mu", "sigma")     // a Dataset: a parameter can name a series
   An array parameter can be given by position or by name, in any combination with number parameters, and the parameters that are not given take their defaults
   (pnorm(x, sd: sigma); pt(x, 8, ncp: ncps)); the same for a matrix or a submatrix (pnorm(M, mean: Mu, sd: S)). With a matrix the elements of an array parameter are taken
   column by column. Arrays of different lengths throw; an invalid element throws; a missing value (NaN) in an array
   parameter gives NaN for that element (a missing single-number parameter throws).
   a non-positive number of degrees of freedom, ...) throw GretlException.

       distribution             parameters (defaults)                       functions
       normal                   mean = 0, sd = 1                            dnorm pnorm qnorm
       Student's t              df, ncp = 0                                 dt pt qt
       chi-square               df, ncp = 0                                 dchisq pchisq qchisq
       F                        df1, df2 (integers), ncp = 0                dF pF qF
       gamma                    shape, scale = 1                            dgamma pgamma qgamma
       beta                     shape1, shape2                              dbeta pbeta qbeta
       exponential              rate = 1                                    dexp pexp qexp
       Weibull                  shape, scale = 1                            dweibull pweibull qweibull
       logistic                 location = 0, scale = 1                     dlogis plogis qlogis
       uniform                  min = 0, max = 1                            dunif punif qunif
       Laplace                  mean = 0, scale = 1                         dlaplace plaplace qlaplace
       generalized error (GED)  shape                                       dged pged qged
       truncated normal         mean = 0, sd = 1, lo = -inf, hi = inf       dtnorm ptnorm qtnorm
       binomial                 size (integer), prob                        dbinom pbinom qbinom
       Poisson                  lambda                                      dpois ppois qpois

   The gamma, Weibull and Laplace scales are scales (the mean of the exponential with rate r is 1/r, that of the gamma is shape * scale), as in the library. The chi-square
   distribution accepts any positive df; the F distribution only integer degrees of freedom (the library's restriction). The numerical work is the library's
   (the cephes routines behind gretl.pvalues), except for the exponential, Weibull, logistic, uniform and beta quantiles and tails, and the truncated normal, which are
   closed forms here. */

module dgretl.distributions;

public import dgretl.matrix;
public import dgretl.exception;

import std.exception: enforce;
import std.math: exp, expm1, floor, fabs, isFinite, isInfinity, isNaN, log, log1p, pow, sqrt, PI;
import std.mathspecial: betaIncomplete, betaIncompleteInverse;
import std.conv: to;
import std.meta: allSatisfy, anySatisfy;
import std.traits: isArray, isNumeric, Unqual;

import dgretl.dataset: Dataset;
import dgretl.internal: fresh;
import gretl.pvalues: gretl_get_cdf, gretl_get_cdf_inverse, gretl_get_critval, gretl_get_pdf, gretl_get_pvalue, normal_cdf, normal_cdf_comp, normal_cdf_inverse;
import gretl.pvaluesh: DistCode;

private enum inf = double.infinity;
private enum nan = double.nan;

/* ---- one element at a time, and the shape of the argument ---- */

private enum isData(T) = isNumeric!T || (isArray!T && isNumeric!(typeof(T.init[0]))) || isMatrixLike!T;

/* a parameter is a number, or one value per element: an array or a matrix with as many elements as the first argument (column by column) */
private enum isParam(T) = isData!T;
private enum isArr(T) = isArray!T && isNumeric!(typeof(T.init[0]));
private struct Unset {}     // the type of an array parameter that was not given: the overload is not viable
/* the same for a matrix and a submatrix parameter, and the default that converts to the matrix type: the default of a parameter that follows one with a default */
private struct UnsetM { GretlMatrix m; alias m this; }
private struct UnsetS { Submatrix s; alias s this; }

private size_t countOf(T)(T a) {
    static if (isNumeric!T) {
        return 1;
    } else static if (isArray!T) {
        return a.length;
    } else {
        return a.nrow * a.ncol;
    }
}

private double elementOf(T)(T a, size_t i) @trusted {
    static if (isNumeric!T) {
        return cast(double) a;
    } else static if (isArray!T) {
        return cast(double) a[i];
    } else {
        return a[i % a.nrow, i / a.nrow];
    }
}

/* f(value, parameters) for every element of x, the parameters being numbers or one value per element. The result has the shape of x (a number if x and all the
   parameters are numbers; an array if x is a number and a parameter is not). A missing value (NaN) in a parameter that is not a single number gives NaN for
   that element; one in a single number is an error, which f reports. */
private auto broadcast(X, P...)(X x, scope double delegate(double, const(double)[]) @system f, P ps) @trusted {
    static assert(P.length <= 8);
    enum scalarResult = isNumeric!X && allSatisfy!(isNumeric, P);

    size_t n = 1;
    bool joined = false;

    void join(size_t m) {
        if (!joined) {
            n = m;
            joined = true;
        } else {
            enforce!GretlException(m == n, "the arguments must have the same number of elements (" ~ to!string(n) ~ " and " ~ to!string(m) ~ ")");
        }
    }

    static if (!isNumeric!X) {
        join(countOf(x));
    }
    static foreach (k, T; P) {
        static if (!isNumeric!T) {
            join(countOf(ps[k]));
        }
    }

    double one(size_t i) {
        double[8] buf;
        bool missing = false;

        static foreach (k, T; P) {
            buf[k] = elementOf(ps[k], i);
            static if (!isNumeric!T) {
                if (isNaN(buf[k])) {
                    missing = true;
                }
            }
        }
        return missing ? nan : f(elementOf(x, i), buf[0 .. P.length]);
    }

    static if (scalarResult) {
        return one(0);
    } else static if (isMatrixLike!X) {
        auto r = fresh(x.nrow, x.ncol);

        foreach (i; 0 .. n) {
            r[i % x.nrow, i / x.nrow] = one(i);
        }
        return r;
    } else {
        auto r = new double[n];

        foreach (i; 0 .. n) {
            r[i] = one(i);
        }
        return r;
    }
}

/* the Dataset forms: a parameter may name a series */
private auto seriesOrValue(A)(Dataset ds, A a) {
    static if (is(A == string)) {
        return ds[a];
    } else {
        return a;
    }
}

/* ---- the library's generic functions ---- */

private struct Law {
    DistCode code;
    double a = 0, b = 0, c = 0;

    double pdf(double x) @trusted {
        double[3] p = [a, b, c];
        return gretl_get_pdf(code, p.ptr, x);
    }

    double cdf(double x) @trusted {
        double[3] p = [a, b, c];
        return gretl_get_cdf(code, p.ptr, x);
    }

    /* P(X > x); the library has it (accurately in the tail) for the common distributions */
    double upper(double x) @trusted {
        switch (code) {
            case DistCode.normal, DistCode.student, DistCode.chisq, DistCode.snedecor, DistCode.gamma, DistCode.binomial, DistCode.poisson, DistCode.ged, DistCode.laplace: {
                double[3] p = [a, b, c];
                return gretl_get_pvalue(code, p.ptr, x);
            }
            default:
                return 1 - cdf(x);
        }
    }

    double inverse(double p) @trusted {
        double[3] q = [a, b, c];
        return gretl_get_cdf_inverse(code, q.ptr, p);
    }

    /* the x with P(X > x) = p */
    double critical(double p) @trusted {
        switch (code) {
            case DistCode.normal, DistCode.student, DistCode.chisq, DistCode.snedecor: {
                double[3] q = [a, b, c];
                return gretl_get_critval(code, q.ptr, p);
            }
            default:
                return inverse(1 - p);
        }
    }
}

/* ---- support, ties, edges ---- */

private struct Support {
    double lo = -inf, hi = inf;
    bool discrete;
}

private double densityOn(double x, Support s, scope double delegate(double) @system f) @trusted {
    if (isNaN(x)) {
        return nan;
    }
    if (x < s.lo || x > s.hi || isInfinity(x)) {
        return 0;
    }
    if (s.discrete && x != floor(x)) {
        return 0;
    }
    return f(x);
}

private double tailOn(double x, Support s, bool lower, scope double delegate(double) @system lowerF, scope double delegate(double) @system upperF) @trusted {
    if (isNaN(x)) {
        return nan;
    }
    if (s.discrete) {
        x = floor(x);
    }
    const below = s.discrete ? (x < s.lo) : (x <= s.lo);
    const above = x >= s.hi || (s.discrete && x > int.max);

    if (below || x == -inf) {
        return lower ? 0 : 1;
    }
    if (above) {
        return lower ? 1 : 0;
    }
    return lower ? lowerF(x) : upperF(x);
}

private double quantileOn(double p, Support s, bool lower, scope double delegate(double, bool) @system f) @trusted {
    if (isNaN(p) || p < 0 || p > 1) {
        return nan;
    }
    if (p == 0) {
        return lower ? s.lo : s.hi;
    }
    if (p == 1) {
        return lower ? s.hi : s.lo;
    }
    return f(p, lower);
}

/* ---- parameter checks ---- */

private void positive(string who, string name, double v) {
    enforce!GretlException(isFinite(v) && v > 0, who ~ ": " ~ name ~ " must be positive and finite");
}

private void finite(string who, string name, double v) {
    enforce!GretlException(isFinite(v), who ~ ": " ~ name ~ " must be finite");
}

private void wholeNumber(string who, string name, double v, double least) {
    enforce!GretlException(isFinite(v) && v == floor(v) && v >= least && v <= int.max, who ~ ": " ~ name ~ " must be an integer of at least " ~ to!string(least));
}

/* ---- the distributions ---- */

/* normal, with the standard normal of the library at the centre */

private void impl_normCheck(double mean, double sd) {
    finite("normal", "mean", mean);
    positive("normal", "sd", sd);
}

private double impl_normPdf(double x, double mean, double sd) {
    if (isNaN(x)) {
        return nan;
    }
    const z = (x - mean) / sd;

    return exp(-0.5 * z * z) / (sd * sqrt(2 * PI));
}

private double impl_normCdf(double x, double mean, double sd, bool lower) {
    if (isNaN(x)) {
        return nan;
    }
    const z = (x - mean) / sd;

    if (z == inf || z == -inf) {
        return ((z > 0) == lower) ? 1 : 0;
    }
    return lower ? normal_cdf(z) : normal_cdf_comp(z);
}

private double impl_normQuant(double p, double mean, double sd, bool lower) {
    return quantileOn(p, Support(), lower, (double q, bool lw) => mean + sd * (lw ? normal_cdf_inverse(q) : -normal_cdf_inverse(q)));
}

/* Student's t (and the non-central t) */

private Law tLaw(double df, double ncp) {
    return (ncp == 0) ? Law(DistCode.student, df) : Law(DistCode.nc_t, df, ncp);
}

private void impl_tCheck(double df, double ncp) {
    positive("t", "df", df);
    finite("t", "ncp", ncp);
}

private double impl_tPdf(double x, double df, double ncp) {
    auto g = tLaw(df, ncp);

    return densityOn(x, Support(), (double v) => g.pdf(v));
}

private double impl_tCdf(double x, double df, double ncp, bool lower) {
    auto g = tLaw(df, ncp);

    return tailOn(x, Support(), lower, (double v) => g.cdf(v), (double v) => g.upper(v));
}

private double impl_tQuant(double p, double df, double ncp, bool lower) {
    auto g = tLaw(df, ncp);

    return quantileOn(p, Support(), lower, (double q, bool lw) => lw ? g.inverse(q) : g.critical(q));
}

/* chi-square: the library's chi-square for an integer df, the gamma distribution with shape df/2 and scale 2 otherwise */

private Law chisqLaw(double df, double ncp) {
    if (ncp != 0) {
        return Law(DistCode.nc_chisq, df, ncp);
    }
    if (df == floor(df) && df <= int.max) {
        return Law(DistCode.chisq, df);
    }
    return Law(DistCode.gamma, df / 2, 2);
}

private void impl_chisqCheck(double df, double ncp) {
    positive("chisq", "df", df);
    enforce!GretlException(isFinite(ncp) && ncp >= 0, "chisq: ncp must be non-negative and finite");
}

private double impl_chisqPdf(double x, double df, double ncp) {
    auto g = chisqLaw(df, ncp);

    return densityOn(x, Support(0, inf), (double v) => g.pdf(v));
}

private double impl_chisqCdf(double x, double df, double ncp, bool lower) {
    auto g = chisqLaw(df, ncp);

    return tailOn(x, Support(0, inf), lower, (double v) => g.cdf(v), (double v) => g.upper(v));
}

private double impl_chisqQuant(double p, double df, double ncp, bool lower) {
    auto g = chisqLaw(df, ncp);

    return quantileOn(p, Support(0, inf), lower, (double q, bool lw) => lw ? g.inverse(q) : g.critical(q));
}

/* F */

private Law fLaw(double d1, double d2, double ncp) {
    return (ncp == 0) ? Law(DistCode.snedecor, d1, d2) : Law(DistCode.nc_f, d1, d2, ncp);
}

private void impl_FCheck(double df1, double df2, double ncp) {
    wholeNumber("F", "df1", df1, 1);
    wholeNumber("F", "df2", df2, 1);
    enforce!GretlException(isFinite(ncp) && ncp >= 0, "F: ncp must be non-negative and finite");
}

private double impl_FPdf(double x, double df1, double df2, double ncp) {
    auto g = fLaw(df1, df2, ncp);

    return densityOn(x, Support(0, inf), (double v) => g.pdf(v));
}

private double impl_FCdf(double x, double df1, double df2, double ncp, bool lower) {
    auto g = fLaw(df1, df2, ncp);

    return tailOn(x, Support(0, inf), lower, (double v) => g.cdf(v), (double v) => g.upper(v));
}

private double impl_FQuant(double p, double df1, double df2, double ncp, bool lower) {
    auto g = fLaw(df1, df2, ncp);

    return quantileOn(p, Support(0, inf), lower, (double q, bool lw) => lw ? g.inverse(q) : g.critical(q));
}

/* gamma */

private void impl_gammaCheck(double shape, double scale) {
    positive("gamma", "shape", shape);
    positive("gamma", "scale", scale);
}

private double impl_gammaPdf(double x, double shape, double scale) {
    auto g = Law(DistCode.gamma, shape, scale);

    return densityOn(x, Support(0, inf), (double v) => g.pdf(v));
}

private double impl_gammaCdf(double x, double shape, double scale, bool lower) {
    auto g = Law(DistCode.gamma, shape, scale);

    return tailOn(x, Support(0, inf), lower, (double v) => g.cdf(v), (double v) => g.upper(v));
}

private double impl_gammaQuant(double p, double shape, double scale, bool lower) {
    auto g = Law(DistCode.gamma, shape, scale);

    return quantileOn(p, Support(0, inf), lower, (double q, bool lw) => lw ? g.inverse(q) : g.inverse(1 - q));
}

/* beta */

private void impl_betaCheck(double shape1, double shape2) {
    positive("beta", "shape1", shape1);
    positive("beta", "shape2", shape2);
}

private double impl_betaPdf(double x, double a, double b) {
    auto g = Law(DistCode.beta, a, b);

    return densityOn(x, Support(0, 1), (double v) => g.pdf(v));
}

private double impl_betaCdf(double x, double a, double b, bool lower) {
    return tailOn(x, Support(0, 1), lower, (double v) => betaIncomplete(a, b, v), (double v) => betaIncomplete(b, a, 1 - v));
}

private double impl_betaQuant(double p, double a, double b, bool lower) {
    return quantileOn(p, Support(0, 1), lower, (double q, bool lw) => lw ? betaIncompleteInverse(a, b, q) : 1 - betaIncompleteInverse(b, a, q));
}

/* exponential, with the rate */

private void impl_expCheck(double rate) {
    positive("exp", "rate", rate);
}

private double impl_expPdf(double x, double rate) {
    return densityOn(x, Support(0, inf), (double v) => rate * exp(-rate * v));
}

private double impl_expCdf(double x, double rate, bool lower) {
    return tailOn(x, Support(0, inf), lower, (double v) => -expm1(-rate * v), (double v) => exp(-rate * v));
}

private double impl_expQuant(double p, double rate, bool lower) {
    return quantileOn(p, Support(0, inf), lower, (double q, bool lw) => (lw ? -log1p(-q) : -log(q)) / rate);
}

/* Weibull */

private void impl_weibullCheck(double shape, double scale) {
    positive("weibull", "shape", shape);
    positive("weibull", "scale", scale);
}

private double impl_weibullPdf(double x, double shape, double scale) {
    return densityOn(x, Support(0, inf), (double v) => (shape / scale) * pow(v / scale, shape - 1) * exp(-pow(v / scale, shape)));
}

private double impl_weibullCdf(double x, double shape, double scale, bool lower) {
    return tailOn(x, Support(0, inf), lower, (double v) => -expm1(-pow(v / scale, shape)), (double v) => exp(-pow(v / scale, shape)));
}

private double impl_weibullQuant(double p, double shape, double scale, bool lower) {
    return quantileOn(p, Support(0, inf), lower, (double q, bool lw) => scale * pow(lw ? -log1p(-q) : -log(q), 1 / shape));
}

/* logistic */

private void impl_logisCheck(double location, double scale) {
    finite("logis", "location", location);
    positive("logis", "scale", scale);
}

private double impl_logisPdf(double x, double location, double scale) {
    return densityOn(x, Support(), (double v) {
        const e = exp(-fabs((v - location) / scale));

        return e / (scale * (1 + e) * (1 + e));
    });
}

private double impl_logisCdf(double x, double location, double scale, bool lower) {
    return tailOn(x, Support(), lower, (double v) => 1 / (1 + exp(-(v - location) / scale)), (double v) => 1 / (1 + exp((v - location) / scale)));
}

private double impl_logisQuant(double p, double location, double scale, bool lower) {
    return quantileOn(p, Support(), lower, (double q, bool lw) => location + scale * (lw ? log(q) - log1p(-q) : log1p(-q) - log(q)));
}

/* uniform */

private void impl_unifCheck(double min, double max) {
    finite("unif", "min", min);
    finite("unif", "max", max);
    enforce!GretlException(min < max, "unif: min must be less than max");
}

private double impl_unifPdf(double x, double min, double max) {
    return densityOn(x, Support(min, max), (double v) => 1 / (max - min));
}

private double impl_unifCdf(double x, double min, double max, bool lower) {
    return tailOn(x, Support(min, max), lower, (double v) => (v - min) / (max - min), (double v) => (max - v) / (max - min));
}

private double impl_unifQuant(double p, double min, double max, bool lower) {
    return quantileOn(p, Support(min, max), lower, (double q, bool lw) => lw ? min + q * (max - min) : max - q * (max - min));
}

/* Laplace */

private void impl_laplaceCheck(double mean, double scale) {
    finite("laplace", "mean", mean);
    positive("laplace", "scale", scale);
}

private double impl_laplacePdf(double x, double mean, double scale) {
    auto g = Law(DistCode.laplace, mean, scale);

    return densityOn(x, Support(), (double v) => g.pdf(v));
}

private double impl_laplaceCdf(double x, double mean, double scale, bool lower) {
    auto g = Law(DistCode.laplace, mean, scale);

    return tailOn(x, Support(), lower, (double v) => g.cdf(v), (double v) => g.upper(v));
}

private double impl_laplaceQuant(double p, double mean, double scale, bool lower) {
    auto g = Law(DistCode.laplace, mean, scale);

    return quantileOn(p, Support(), lower, (double q, bool lw) => lw ? g.inverse(q) : g.inverse(1 - q));
}

/* generalized error distribution */

private void impl_gedCheck(double shape) {
    positive("ged", "shape", shape);
}

private double impl_gedPdf(double x, double shape) {
    auto g = Law(DistCode.ged, shape);

    return densityOn(x, Support(), (double v) => g.pdf(v));
}

private double impl_gedCdf(double x, double shape, bool lower) {
    auto g = Law(DistCode.ged, shape);

    return tailOn(x, Support(), lower, (double v) => g.cdf(v), (double v) => g.upper(v));
}

private double impl_gedQuant(double p, double shape, bool lower) {
    auto g = Law(DistCode.ged, shape);

    return quantileOn(p, Support(), lower, (double q, bool lw) => lw ? g.inverse(q) : g.inverse(1 - q));
}

/* truncated normal: the normal(mean, sd) conditioned on [lo, hi]. Standardized limits a, b; the tail that is the smaller one is
   computed from the upper tails when the interval lies on the right, so that no precision is lost. */

private void impl_tnormCheck(double mean, double sd, double lo, double hi) {
    finite("tnorm", "mean", mean);
    positive("tnorm", "sd", sd);
    enforce!GretlException(!isNaN(lo) && !isNaN(hi) && lo < hi, "tnorm: lo must be less than hi");
}

private double Phi(double z) {
    return (z == inf) ? 1 : (z == -inf) ? 0 : normal_cdf(z);
}

private double Qhi(double z) {
    return (z == inf) ? 0 : (z == -inf) ? 1 : normal_cdf_comp(z);
}

private double impl_tnormPdf(double x, double mean, double sd, double lo, double hi) {
    const a = (lo - mean) / sd, b = (hi - mean) / sd;
    const mass = (a > 0) ? Qhi(a) - Qhi(b) : Phi(b) - Phi(a);

    return densityOn(x, Support(lo, hi), (double v) {
        const z = (v - mean) / sd;

        return exp(-0.5 * z * z) / (sd * sqrt(2 * PI) * mass);
    });
}

private double impl_tnormCdf(double x, double mean, double sd, double lo, double hi, bool lower) {
    const a = (lo - mean) / sd, b = (hi - mean) / sd;
    const right = a > 0;
    const mass = right ? Qhi(a) - Qhi(b) : Phi(b) - Phi(a);

    return tailOn(x, Support(lo, hi), lower,
        (double v) { const z = (v - mean) / sd; return right ? (Qhi(a) - Qhi(z)) / mass : (Phi(z) - Phi(a)) / mass; },
        (double v) { const z = (v - mean) / sd; return right ? (Qhi(z) - Qhi(b)) / mass : (Phi(b) - Phi(z)) / mass; });
}

private double impl_tnormQuant(double p, double mean, double sd, double lo, double hi, bool lower) {
    const a = (lo - mean) / sd, b = (hi - mean) / sd;
    const right = a > 0;
    const mass = right ? Qhi(a) - Qhi(b) : Phi(b) - Phi(a);

    return quantileOn(p, Support(lo, hi), lower, (double q, bool lw) {
        /* the probability, as a lower tail of the normal (or an upper tail when the interval is on the right) */
        double z;

        if (right) {
            const u = lw ? Qhi(a) - q * mass : Qhi(b) + q * mass;      // upper-tail probability of the standard normal
            z = -normal_cdf_inverse(u);
        } else {
            const u = lw ? Phi(a) + q * mass : Phi(b) - q * mass;
            z = normal_cdf_inverse(u);
        }
        const x = mean + sd * z;

        return (x < lo) ? lo : (x > hi) ? hi : x;
    });
}

/* binomial */

private void impl_binomCheck(double size, double prob) {
    wholeNumber("binom", "size", size, 0);
    enforce!GretlException(prob >= 0 && prob <= 1, "binom: prob must be in [0, 1]");
}

private double impl_binomPdf(double x, double size, double prob) {
    auto g = Law(DistCode.binomial, prob, size);

    return densityOn(x, Support(0, size, true), (double v) => g.pdf(v));
}

private double impl_binomCdf(double x, double size, double prob, bool lower) {
    auto g = Law(DistCode.binomial, prob, size);

    return tailOn(x, Support(0, size, true), lower, (double v) => g.cdf(v), (double v) => g.upper(v));
}

private double impl_binomQuant(double p, double size, double prob, bool lower) {
    auto g = Law(DistCode.binomial, prob, size);

    return quantileOn(p, Support(0, size, true), lower, (double q, bool lw) => g.inverse(lw ? q : 1 - q));
}

/* Poisson */

private void impl_poisCheck(double lambda) {
    positive("pois", "lambda", lambda);
}

private double impl_poisPdf(double x, double lambda) {
    auto g = Law(DistCode.poisson, lambda);

    return densityOn(x, Support(0, inf, true), (double v) => g.pdf(v));
}

private double impl_poisCdf(double x, double lambda, bool lower) {
    auto g = Law(DistCode.poisson, lambda);

    return tailOn(x, Support(0, inf, true), lower, (double v) => g.cdf(v), (double v) => g.upper(v));
}

private double impl_poisQuant(double p, double lambda, bool lower) {
    auto g = Law(DistCode.poisson, lambda);

    return quantileOn(p, Support(0, inf, true), lower, (double q, bool lw) => g.inverse(lw ? q : 1 - q));
}

/* ---- the public functions: dX, pX, qX for each distribution, generated from the pieces above ---- */

/* The functions of one distribution. `names` and `defaults` are its parameters (a default of "" means that the parameter is required). For every one of dX, pX and qX:
     - the parameters are plain numbers, with their defaults: pnorm(x), pnorm(x, sd: 2), pt(x, 8, lower: false);
     - some of the parameters are arrays, matrices or submatrices, one value per element of x (and the others numbers), by position or by name, in any combination: pnorm(x, mu), pnorm(x, sd: sigma),
       pnorm(x, 0, sigma), pt(x, 8, ncp: ncps). A parameter that is an array is recognised by its type, so every overload below has a set of parameters that are arrays
       (`varying`) and the others are numbers; a parameter that is not given takes its default;
     - the same with a Dataset and the name of a series, where an array parameter may also be the name of a series: pnorm(d, "y", "mu", "sigma");
     (An array, a matrix and a submatrix are different kinds of parameter, so for each parameter an overload exists for each kind, and the overloads are made for every
     combination of kinds. The type of an argument selects the kind, and there is no ambiguity.) */
private string distribution(string key, string[] names, string[] defaults)() {
    const k = names.length;
    size_t required = 0;

    foreach (d; defaults) {
        required += d.length ? 0 : 1;
    }
    string r;

    foreach (kind; ["d:Pdf", "p:Cdf", "q:Quant"]) {
        const fn = kind[0 .. 1], impl = kind[2 .. $];
        const tail = (fn == "d") ? "" : ", bool lower = true";
        const lowerArg = (fn == "d") ? "" : ", lower";
        const arg0 = (fn == "d") ? "x" : (fn == "p") ? "q" : "p";
        const lam = "(double v, const(double)[] a) { impl_" ~ key ~ "Check(" ~ unpack(k) ~ "); return impl_" ~ key ~ impl ~ "(v, " ~ unpack(k) ~ lowerArg ~ "); }";

        // numbers
        string decl, args;

        foreach (n, nm; names) {
            decl ~= ", double " ~ nm ~ (defaults[n].length ? " = " ~ defaults[n] : "");
            args ~= ", " ~ nm;
        }
        r ~= "auto " ~ fn ~ key ~ "(X)(X " ~ arg0 ~ decl ~ tail ~ ") @trusted if (isData!X) { return broadcast(" ~ arg0 ~ ", " ~ lam ~ args ~ "); }\n";
        r ~= "double[] " ~ fn ~ key ~ "(Dataset ds, string series" ~ decl ~ tail ~ ") @trusted { return " ~ fn ~ key ~ "(ds[series]" ~ args ~ lowerArg ~ "); }\n";

        // some parameters are arrays, matrices or submatrices (the kind of each is its digit 1, 2 or 3 below; 0 is a number): every combination with at least one
        foreach (code; 1 .. 4 ^^ k) {
            string tmpl = "X", d2, con = "isData!X", dsTmpl, dsCon = "true", forward;
            bool anyDefault = false;
            size_t rest = code;

            foreach (n, nm; names) {
                const kindOf = rest % 4;
                const A = "A" ~ to!string(n);

                rest /= 4;
                if (kindOf == 0) {
                    d2 ~= ", double " ~ nm ~ (defaults[n].length ? " = " ~ defaults[n] : "");
                    anyDefault = anyDefault || defaults[n].length > 0;
                    forward ~= ", " ~ nm;
                } else {
                    const unset = (kindOf == 1) ? "Unset" : (kindOf == 2) ? "UnsetM" : "UnsetS";

                    tmpl ~= ", " ~ A ~ " = " ~ unset;
                    d2 ~= ", " ~ A ~ " " ~ nm ~ (anyDefault ? " = " ~ ((kindOf == 1) ? "null" : unset ~ "()") : "");
                    con ~= " && " ~ ((kindOf == 1) ? "isArr!" ~ A : (kindOf == 2) ? "is(Unqual!" ~ A ~ " == GretlMatrix)" : "is(Unqual!" ~ A ~ " == Submatrix)");
                    dsTmpl ~= (dsTmpl.length ? ", " : "") ~ A ~ " = " ~ unset;
                    dsCon ~= " && " ~ ((kindOf == 1) ? "(isArr!" ~ A ~ " || is(" ~ A ~ " == string))" : (kindOf == 2) ? "is(Unqual!" ~ A ~ " == GretlMatrix)" : "is(Unqual!" ~ A ~ " == Submatrix)");
                    forward ~= ", seriesOrValue(ds, " ~ nm ~ ")";
                }
            }
            r ~= "auto " ~ fn ~ key ~ "(" ~ tmpl ~ ")(X " ~ arg0 ~ d2 ~ tail ~ ") @trusted if (" ~ con ~ ") { return broadcast(" ~ arg0 ~ ", " ~ lam ~ args ~ "); }\n";
            r ~= "double[] " ~ fn ~ key ~ "(" ~ dsTmpl ~ ")(Dataset ds, string series" ~ d2 ~ tail ~ ") @trusted if (" ~ dsCon ~ ") { return " ~ fn ~ key ~ "(ds[series]" ~ forward ~ lowerArg ~ "); }\n";
        }
    }
    return r;
}

private string unpack(size_t n) {
    string r;

    foreach (k; 0 .. n) {
        r ~= (k ? ", " : "") ~ "a[" ~ to!string(k) ~ "]";
    }
    return r;
}

/* normal */
mixin(distribution!("norm", ["mean", "sd"], ["0", "1"]));
/* Student's t */
mixin(distribution!("t", ["df", "ncp"], ["", "0"]));
/* chi-square */
mixin(distribution!("chisq", ["df", "ncp"], ["", "0"]));
/* F */
mixin(distribution!("F", ["df1", "df2", "ncp"], ["", "", "0"]));
/* gamma */
mixin(distribution!("gamma", ["shape", "scale"], ["", "1"]));
/* beta */
mixin(distribution!("beta", ["shape1", "shape2"], ["", ""]));
/* exponential */
mixin(distribution!("exp", ["rate"], ["1"]));
/* Weibull */
mixin(distribution!("weibull", ["shape", "scale"], ["", "1"]));
/* logistic */
mixin(distribution!("logis", ["location", "scale"], ["0", "1"]));
/* uniform */
mixin(distribution!("unif", ["min", "max"], ["0", "1"]));
/* Laplace */
mixin(distribution!("laplace", ["mean", "scale"], ["0", "1"]));
/* generalized error */
mixin(distribution!("ged", ["shape"], [""]));
/* truncated normal */
mixin(distribution!("tnorm", ["mean", "sd", "lo", "hi"], ["0", "1", "-double.infinity", "double.infinity"]));
/* binomial */
mixin(distribution!("binom", ["size", "prob"], ["", ""]));
/* Poisson */
mixin(distribution!("pois", ["lambda"], [""]));
