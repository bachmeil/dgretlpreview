/* Descriptive statistics: the functionality of describe.d (the library's summary statistics, quantiles, moments, correlations, frequency distributions,
   normality tests, two-sample tests, distances) in a D-friendly form. Every function takes plain data as an array (double[], int[], ...), a GretlMatrix or Submatrix,
   or one row or column of a matrix (Col(m, j), Row(m, i)), and the results are plain structs and numbers; nothing is freed by hand and every function throws
   GretlException if the library reports an error.

   import dgretl;
   double m = mean(x);                      // x: double[], a column Col(M, 0) ...
   auto means = mean(M);                    // M: a matrix: the statistic of each column, as a 1 x k GretlMatrix
   auto q = quantile(x, 0.9);
   auto s = summary(M, names: ["a", "b"]);  // everything at once: s.mean, s.sd, s.min, s.max, s.median, ...
   auto c = correlation(M);                 // k x k correlation matrix

   Missing values are NaN, and (as in the library) the statistics skip them: mean([1.0, double.nan, 3.0]) is 2. A statistic of no valid observation is NaN.
   The names are those of the library, except that `total` is the sum and `validCount` the number of valid values (so that they cannot clash with std.algorithm's `sum` and
   `count`), `rowDistances` is the library's `distance`, and `kurtosis` is the EXCESS kurtosis
   (zero for a normal distribution), as in the library. */

module dgretl.descriptive;

public import dgretl.matrix;
public import dgretl.exception;

import std.conv: to;
import std.exception: enforce;
import std.format: format;
import std.math: isNaN, sqrt;
import std.range: ElementType;
import std.array: split;
import std.string: fromStringz, indexOf, strip;

import dgretl.dataset: Dataset;
import dgretl.internal;
import gretl.libset: pop_program_state, push_program_state;
import gretl.describe: libDistanceFn = distance, free_mahal_dist, free_xtab, get_mahal_distances, mahal_dist_get_distances, mahal_dist_get_n, multivariate_normality_test,
    gretl_normtest_matrix, means_test, vars_test, satterthwaite_df, single_crosstab, free_freq, free_summary, get_freq, get_summary, get_summary_weighted, gretl_corr, gretl_covar, gretl_gini, gretl_kurtosis, gretl_long_run_variance, gretl_max, gretl_mean, gretl_median, gretl_min,
    gretl_quantile, gretl_skewness, gretl_sst, gretl_stddev, gretl_sum, gretl_variance, real_gretl_variance;
import gretl.gretl_errors: gretl_errmsg_get;
import gretl.gretl_list: gretl_list_new;
import gretl.libgretlh: Gretlopt;
import gretl.nonparam: rank_sum_test, sign_test, signed_rank_test;
import gretl.gretl_matrix: gretl_matrix_free;
import gretl.gretl_matrixh: GretlMatrixRaw;
import gretl.gretl_prn: Prn, gretl_print_destroy, gretl_print_new;
import gretl.gretl_prnh: PrnType;
import gretl.gretl_utils: get_last_pvalue, get_last_test_statistic;
import gretl.pvalues: binomial_cdf, chisq_cdf_comp, normal_cdf, normal_pvalue_1, snedecor_pvalue_2, student_pvalue_2;
import gretl.libset: pop_program_state, push_program_state;

/* ---- what the functions take ---- */

/* a series: an array of numbers, or one row or column of a matrix */
enum isSeries(T) = isNumericArray!T || is(T == MatrixLine);

private const(double)[] asDoubles(S)(S x) @trusted if (isSeries!S) {
    static if (is(S == MatrixLine)) {
        return x.dup;
    } else static if (is(ElementType!S == double)) {
        return x;
    } else {
        double[] a;

        foreach (v; x) a ~= v;
        return a;
    }
}

/* the statistic of each column of a matrix, as a 1 x k matrix */
private GretlMatrix perColumn(M)(M m, scope double delegate(const(double)[]) @system f) @trusted if (isMatrixLike!M) {
    auto v = cbind(m);                  // a contiguous copy
    auto r = fresh(1, v.ncol);

    foreach (j; 0 .. v.ncol) r.ptr.val[j] = f(v.ptr.val[j * v.nrow .. (j + 1) * v.nrow]);
    return r;
}

private string lastError() @trusted {
    return fromStringz(gretl_errmsg_get()).idup;
}

/* the pairs of functions: for a series the number, for a matrix the 1 x k matrix of the columns' numbers */
private string pair(string name, string core, string extraParams = "", string extraArgs = "") {
    return "double " ~ name ~ "(S)(S x" ~ extraParams ~ ") @trusted if (isSeries!S) { return " ~ core ~ "(asDoubles(x)" ~ extraArgs ~ "); }\n"
         ~ "GretlMatrix " ~ name ~ "(M)(M m" ~ extraParams ~ ") @trusted if (isMatrixLike!M) { return perColumn(m, (const(double)[] a) => " ~ core ~ "(a" ~ extraArgs ~ ")); }\n";
}

/* ---- location, spread and shape ---- */

private double meanA(const(double)[] a) @trusted { return a.length ? gretl_mean(0, cast(int) a.length - 1, a.ptr) : double.nan; }
private double medianA(const(double)[] a) @trusted { return a.length ? gretl_median(0, cast(int) a.length - 1, a.ptr) : double.nan; }
private double totalA(const(double)[] a) @trusted { return a.length ? gretl_sum(0, cast(int) a.length - 1, a.ptr) : double.nan; }
private double minimumA(const(double)[] a) @trusted { return a.length ? gretl_min(0, cast(int) a.length - 1, a.ptr) : double.nan; }
private double maximumA(const(double)[] a) @trusted { return a.length ? gretl_max(0, cast(int) a.length - 1, a.ptr) : double.nan; }
private double sstA(const(double)[] a) @trusted { return a.length ? gretl_sst(0, cast(int) a.length - 1, a.ptr) : double.nan; }
private double skewnessA(const(double)[] a) @trusted { return a.length ? gretl_skewness(0, cast(int) a.length - 1, a.ptr) : double.nan; }
private double kurtosisA(const(double)[] a) @trusted { return a.length ? gretl_kurtosis(0, cast(int) a.length - 1, a.ptr) : double.nan; }
private double giniA(const(double)[] a) @trusted { return a.length ? gretl_gini(0, cast(int) a.length - 1, a.ptr) : double.nan; }
private double countA(const(double)[] a) @trusted {
    size_t n = 0;

    foreach (v; a) if (!isNaN(v)) n++;
    return n;
}

private double varianceA(const(double)[] a, bool population) @trusted {
    return a.length ? real_gretl_variance(0, cast(int) a.length - 1, a.ptr, population ? 1 : 0) : double.nan;
}

private double stddevA(const(double)[] a, bool population) @trusted {
    const v = varianceA(a, population);

    return isNaN(v) ? v : sqrt(v);
}

/* The rule for the quantile of a sample (Hyndman and Fan's types 6, 7 and 8; the library's setting `quantile_type`, default Q6). */
enum QuantileType {
    q6,     // position (n + 1) p: what the library uses by default
    q7,     // position (n - 1) p + 1: the default of R and numpy
    q8      // position (n + 1/3) p + 1/3: approximately median-unbiased
}

private double quantileA(const(double)[] a, double p, QuantileType type = QuantileType.q6) @trusted {
    enforce!GretlException(p >= 0 && p <= 1, "quantile: the probability must be between 0 and 1");
    if (p == 0) return minimumA(a);
    if (p == 1) return maximumA(a);
    if (a.length == 0) return double.nan;
    enforce!GretlException(push_program_state() == 0, "could not save the library settings");
    scope (exit) pop_program_state();
    setLibrary("quantile_type", ["Q6", "Q7", "Q8"][type]);
    int err = 0;

    return gretl_quantile(0, cast(int) a.length - 1, a.ptr, p, Gretlopt.OPT_Q, &err);
}

private double lrvA(const(double)[] a, int lag, double mean) @trusted {
    return a.length ? gretl_long_run_variance(0, cast(int) a.length - 1, a.ptr, lag, mean) : double.nan;
}

/* mean(x): the arithmetic mean (of each column, for a matrix). */
mixin(pair("mean", "meanA"));

/* median(x): the median (the 0.5 quantile). */
mixin(pair("median", "medianA"));

/* total(x): the sum of the values. */
mixin(pair("total", "totalA"));

/* minimum(x) and maximum(x): the smallest and the largest value. */
mixin(pair("minimum", "minimumA"));
/// ditto
mixin(pair("maximum", "maximumA"));

/* sst(x): the sum of squared deviations from the mean. */
mixin(pair("sst", "sstA"));

/* validCount(x): the number of valid (not missing) values, as a double. */
mixin(pair("validCount", "countA"));

/* variance(x): the sample variance (divisor n - 1), or with population: true the one with divisor n. */
mixin(pair("variance", "varianceA", ", bool population = false", ", population"));

/* stddev(x): the square root of variance(x, population). */
mixin(pair("stddev", "stddevA", ", bool population = false", ", population"));

/* skewness(x): the sample skewness, the third standardized moment (divisors n; NaN for a constant series). */
mixin(pair("skewness", "skewnessA"));

/* kurtosis(x): the EXCESS kurtosis, the fourth standardized moment minus 3 (divisors n; NaN for a constant series). */
mixin(pair("kurtosis", "kurtosisA"));

/* gini(x): the Gini coefficient of inequality. */
mixin(pair("gini", "giniA"));

/* quantile(x, p, type): the p quantile, 0 <= p <= 1 (p = 0 and p = 1 give the minimum and the maximum), by interpolation between the order statistics at the position
   the QuantileType gives. The result is NaN when the sample is too small for p (for the default type, when (n + 1) p is below 1 or above n). */
mixin(pair("quantile", "quantileA", ", double p, QuantileType type = QuantileType.q6", ", p, type"));

/* longRunVariance(x, lag, mean): the long-run variance with Bartlett weights of `lag` lags (-1, the default, takes the library's choice n^(1/3)); `mean`
   is the mean to deviate from (NaN, the default, uses the sample mean). Unlike the other functions it does not skip missing values: they give NaN. */
mixin(pair("longRunVariance", "lrvA", ", int lag = -1, double mean = double.nan", ", lag, mean"));

/* quantiles(x, ps): the quantiles at the probabilities ps, as an array. */
double[] quantiles(S)(S x, const(double)[] ps) if (isSeries!S) {
    double[] r;
    auto a = asDoubles(x);

    foreach (p; ps) r ~= quantileA(a, p);
    return r;
}

/* the five-number summary and more: min, quartiles, max */
struct FiveNumbers {
    double min, q1, median, q3, max;
}

FiveNumbers fiveNumbers(S)(S x) if (isSeries!S) {
    auto a = asDoubles(x);

    return FiveNumbers(minimumA(a), quantileA(a, 0.25), medianA(a), quantileA(a, 0.75), maximumA(a));
}

/* ---- two series: covariance and correlation ---- */

private void checkPair(size_t a, size_t b, string who) {
    enforce!GretlException(a == b, who ~ ": the series have " ~ to!string(a) ~ " and " ~ to!string(b) ~ " observations");
}

/* covariance(x, y): the sample covariance (divisor n - 1) of two series of the same length, using the observations where neither is missing. */
double covariance(S, T)(S x, T y) @trusted if (isSeries!S && isSeries!T) {
    auto a = asDoubles(x), b = asDoubles(y);

    checkPair(a.length, b.length, "covariance");
    return a.length ? gretl_covar(0, cast(int) a.length - 1, a.ptr, b.ptr, null) : double.nan;
}

/* correlation(x, y): the correlation coefficient of two series, using the observations where neither is missing. */
double correlation(S, T)(S x, T y) @trusted if (isSeries!S && isSeries!T) {
    auto a = asDoubles(x), b = asDoubles(y);

    checkPair(a.length, b.length, "correlation");
    return a.length ? gretl_corr(0, cast(int) a.length - 1, a.ptr, b.ptr, null) : double.nan;
}

private GretlMatrix pairwise(M)(M m, bool corr) @trusted if (isMatrixLike!M) {
    auto v = cbind(m);
    const k = v.ncol, n = v.nrow;
    auto r = fresh(k, k);

    foreach (i; 0 .. k) {
        foreach (j; i .. k) {
            auto x = v.ptr.val + i * n, y = v.ptr.val + j * n;
            double c = double.nan;

            if (n > 0) c = corr ? gretl_corr(0, cast(int) n - 1, x, y, null) : gretl_covar(0, cast(int) n - 1, x, y, null);
            r.ptr.val[j * k + i] = c;
            r.ptr.val[i * k + j] = c;
        }
    }
    return r;
}

/* covariance(m): the k x k sample covariance matrix of the columns of a matrix, each entry from the rows where both columns are valid. */
GretlMatrix covariance(M)(M m) @trusted if (isMatrixLike!M) {
    return pairwise(m, false);
}

/* correlation(m): the k x k correlation matrix of the columns of a matrix, each entry from the rows where both columns are valid. */
GretlMatrix correlation(M)(M m) @trusted if (isMatrixLike!M) {
    return pairwise(m, true);
}

/* ---- the summary table ---- */

/* The summary statistics of k series: each field is a 1 x k matrix, in the order of the series. A series with no valid observation has NaN in every field.
   `missing` counts the missing values among the observations that count (with weights, those with a non-zero weight). */
struct SummaryStats {
    string[] names;
    GretlMatrix valid;          // number of valid observations (weighted: with a non-zero weight)
    GretlMatrix missing;        // number of missing observations
    GretlMatrix mean;
    GretlMatrix median;
    GretlMatrix sd;             // standard deviation (divisor n - 1; weighted: the library's weighted formula)
    GretlMatrix min;
    GretlMatrix max;
    GretlMatrix cv;             // coefficient of variation |sd / mean| (NaN if the mean is 0)
    GretlMatrix skewness;
    GretlMatrix excessKurtosis;
    GretlMatrix p5;             // 5th and 95th percentiles
    GretlMatrix p95;
    GretlMatrix iqr;            // interquartile range: the 0.75 quantile minus the 0.25 quantile

    string toString() {
        string s = format("%-10s %6s %10s %10s %10s %10s %10s %10s %10s\n", "", "n", "mean", "median", "min", "max", "s.d.", "skew", "ex.kurt");

        foreach (j, nm; names) {
            s ~= format("%-10s %6.0f %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g\n", nm, valid[0, j], mean[0, j], median[0, j], min[0, j], max[0, j], sd[0, j], skewness[0, j],
                        excessKurtosis[0, j]);
        }
        return s;
    }
}

/* summary(x): the summary statistics (the library's `summary` command) of a series, or of each column of a matrix.
       names     the names of the series (default x1, x2, ...; `x` for a single series)
       weights   one weight per observation: the statistics are weighted; observations with weight 0 or NaN are left out (the same weights for every column) */
SummaryStats summary(X)(X x, string[] names = null, const(double)[] weights = null) @trusted if (isMatrixLike!X || isSeries!X) {
    static if (isMatrixLike!X) {
        auto m = cbind(x);
    } else {
        auto m = cbind(asDoubles(x).dup);
    }
    const k = m.ncol, n = m.nrow;

    enforce!GretlException(k >= 1 && n >= 1, "summary: no data");
    enforce!GretlException(names is null || names.length == k, "summary: " ~ to!string(k) ~ " series need " ~ to!string(k) ~ " names");
    enforce!GretlException(weights is null || weights.length == n, "summary: the weights need one value per observation");
    SummaryStats r;

    r.names = names !is null ? names.dup : (k == 1 ? ["x"] : numbered("x", k));
    foreach (f; [&r.valid, &r.missing, &r.mean, &r.median, &r.sd, &r.min, &r.max, &r.cv, &r.skewness, &r.excessKurtosis, &r.p5, &r.p95, &r.iqr]) {
        *f = fresh(1, k);
        foreach (j; 0 .. k) f.ptr.val[j] = double.nan;
    }
    auto cols = columnsOf(m);
    const(double)[][] all = cols;

    if (weights !is null) all ~= weights;
    auto ds = makeDataset(all, ["const"] ~ r.names ~ (weights !is null ? ["weights"] : null));
    auto list = gretl_list_new(cast(int) k);

    foreach (j; 0 .. k) list[j + 1] = cast(int) (j + 1);
    scope (exit) {
        import core.stdc.stdlib: free;
        free(list);
    }
    int err = 0;
    auto sm = weights !is null ? get_summary_weighted(list, ds.ptr, cast(int) (k + 1), Gretlopt.OPT_NONE, null, &err)
                               : get_summary(list, ds.ptr, Gretlopt.OPT_NONE, null, &err);

    enforce!GretlException(err == 0 && sm !is null, "summary: " ~ lastError());
    scope (exit) free_summary(sm);
    foreach (i; 0 .. cast(size_t) sm.list[0]) {
        const j = cast(size_t) (sm.list[i + 1] - 1);

        size_t used = 0;       // the observations that count: a valid value and, if weighted, a weight that is neither 0 nor missing

        foreach (t; 0 .. n) {
            if (!isNaN(cols[j][t]) && (weights is null || (!isNaN(weights[t]) && weights[t] != 0))) used++;
        }
        r.valid.ptr.val[j] = used;
        r.missing.ptr.val[j] = sm.misscount[i];
        r.mean.ptr.val[j] = sm.mean[i];
        r.median.ptr.val[j] = sm.median[i];
        r.sd.ptr.val[j] = sm.sd[i];
        r.min.ptr.val[j] = sm.low[i];
        r.max.ptr.val[j] = sm.high[i];
        r.cv.ptr.val[j] = sm.cv[i];
        r.skewness.ptr.val[j] = sm.skew[i];
        r.excessKurtosis.ptr.val[j] = sm.xkurt[i];
        r.p5.ptr.val[j] = sm.perc05[i];
        r.p95.ptr.val[j] = sm.perc95[i];
        r.iqr.ptr.val[j] = sm.iqr[i];
    }
    return r;
}

/// ditto, for the named series of a Dataset
SummaryStats summary(Dataset d, string[] names, const(double)[] weights = null) @trusted {
    enforce!GretlException(names.length > 0, "summary: no series named");
    double[][] cols;

    foreach (nm; names) cols ~= d[nm].dup;
    return summary(matrixFrom(cols), names, weights);
}

/* ---- frequency distributions ---- */

/* A frequency distribution: one row per bin (or, if `discrete`, per distinct value), each field nbins x 1. */
struct FrequencyTable {
    bool discrete;              // one row for each distinct value of a series that takes few of them (midpoint is the value)
    GretlMatrix midpoint;       // the middle of each bin, or the value
    GretlMatrix lower;          // the left end of each bin (the right end is the next one's `lower`; see `upper`); equal to the value if discrete
    GretlMatrix upper;          // the right end of each bin; equal to the value if discrete
    GretlMatrix count;          // number of observations in each bin
    GretlMatrix fraction;       // count / nobs
    GretlMatrix cumulative;     // the cumulative fractions
    int nobs;                   // observations used (not missing)
    double mean, sd;            // of the series
    double normalTest = double.nan;     // Doornik-Hansen chi-square(2) statistic for normality, if asked for (and the series has at least 8 values)
    double normalPvalue = double.nan;
}

/* frequencies(x, bins, min, width, normalTest): the frequency distribution of a series (the library's `freq`). The bins are of equal width: by default their number is chosen
   from the sample size (5 to 29, always odd); `bins` fixes the number (at least 2); `min` and `width` (given together) fix the left end of the first bin and the width,
   and the number follows. A series with few distinct values is tabulated by value instead (`discrete`). The first and last bins also hold anything below or above the
   outer ends. A series that is constant, or has fewer than 8 valid values, throws. `normalTest` adds the Doornik-Hansen test of normality. */
FrequencyTable frequencies(S)(S x, int bins = 0, double min = double.nan, double width = double.nan, bool normalTest = false) @trusted if (isSeries!S) {
    enforce!GretlException(bins == 0 || bins >= 2, "frequencies: the number of bins must be at least 2");
    enforce!GretlException(isNaN(min) == isNaN(width), "frequencies: give the minimum and the width together");
    enforce!GretlException(isNaN(width) || width > 0, "frequencies: the width must be positive");
    enforce!GretlException(bins == 0 || isNaN(width), "frequencies: give either the number of bins or the minimum and the width");
    auto a = asDoubles(x);

    enforce!GretlException(a.length >= 3, "frequencies: too few observations");
    return makeFrequencies(a, bins, min, width, normalTest, false);
}

/* valueCounts(x): how often each distinct value of a series occurs, as a FrequencyTable with discrete = true. */
FrequencyTable valueCounts(S)(S x) @trusted if (isSeries!S) {
    auto a = asDoubles(x);

    enforce!GretlException(a.length >= 3, "valueCounts: too few observations");
    return makeFrequencies(a, 0, double.nan, double.nan, false, true);
}

private FrequencyTable makeFrequencies(const(double)[] a, int bins, double min, double width, bool normalTest, bool discrete) @trusted {
    auto ds = makeDataset([a], ["const", "x"]);
    auto opt = Gretlopt.OPT_S;

    if (normalTest) opt |= Gretlopt.OPT_Z;
    if (discrete) opt |= Gretlopt.OPT_D;
    int err = 0;
    auto fd = get_freq(1, ds.ptr, min, width, bins, 1, opt, &err);

    enforce!GretlException(err == 0 && fd !is null, "frequencies: " ~ lastError());
    scope (exit) free_freq(fd);
    FrequencyTable t;
    const k = fd.numbins;

    t.discrete = fd.discrete != 0;
    t.nobs = fd.n;
    t.mean = fd.xbar;
    t.sd = fd.sdx;
    t.midpoint = fresh(k, 1);
    t.lower = fresh(k, 1);
    t.upper = fresh(k, 1);
    t.count = fresh(k, 1);
    t.fraction = fresh(k, 1);
    t.cumulative = fresh(k, 1);
    double cum = 0;

    foreach (i; 0 .. k) {
        t.midpoint.ptr.val[i] = fd.midpt[i];
        t.lower.ptr.val[i] = t.discrete ? fd.midpt[i] : fd.endpt[i];
        t.upper.ptr.val[i] = t.discrete ? fd.midpt[i] : fd.endpt[i + 1];
        t.count.ptr.val[i] = fd.f[i];
        t.fraction.ptr.val[i] = cast(double) fd.f[i] / fd.n;
        cum += t.fraction.ptr.val[i];
        t.cumulative.ptr.val[i] = cum;
    }
    if (!isNaN(fd.test)) {
        t.normalTest = fd.test;
        t.normalPvalue = chisq_cdf_comp(2, fd.test);
    }
    return t;
}

/* ---- normality tests ---- */

/* The result of a test: the statistic and its p-value. */
struct StatTest {
    double statistic = double.nan;
    double pvalue = double.nan;
}

/* The four tests of normality of the library's `normtest --all`; each is NaN if it cannot be computed for the sample. */
struct NormalityTests {
    StatTest doornikHansen;     // chi-square(2), from the skewness and kurtosis
    StatTest shapiroWilk;       // W (the p-value is approximate)
    StatTest jarqueBera;        // chi-square(2)
    StatTest lilliefors;        // Kolmogorov-Smirnov with estimated mean and variance (the p-value is approximate)
}

/* normalityTests(x): the tests of the null that a series is normally distributed (missing values are left out; at least 8 valid values, not constant). */
NormalityTests normalityTests(S)(S x) @trusted if (isSeries!S) {
    double[] a;

    foreach (v; asDoubles(x)) if (!isNaN(v)) a ~= v;
    enforce!GretlException(a.length >= 8, "normalityTests: too few valid observations (at least 8 are needed)");
    int err = 0;
    auto m = gretl_normtest_matrix(a.ptr, 0, cast(int) a.length - 1, Gretlopt.OPT_A, &err);

    enforce!GretlException(err == 0 && m !is null && m.rows == 4, "normalityTests: " ~ lastError());
    scope (exit) gretl_matrix_free(m);
    NormalityTests r;

    r.doornikHansen = StatTest(m.val[0], m.val[4]);
    r.shapiroWilk = StatTest(m.val[1], m.val[5]);
    r.jarqueBera = StatTest(m.val[2], m.val[6]);
    r.lilliefors = StatTest(m.val[3], m.val[7]);
    return r;
}

/* multivariateNormality(e, sigma): the Doornik-Hansen test of multivariate normality of the rows of e (n x p, for example the residuals of a system), chi-square with 2 p
   degrees of freedom. sigma is the covariance matrix of the columns (p x p); by default e'e / n. */
StatTest multivariateNormality(E)(E e, GretlMatrix sigma = GretlMatrix.init) @trusted if (isMatrixLike!E) {
    auto em = cbind(e);

    requireFinite(em, "multivariateNormality: e");
    enforce!GretlException(em.ncol >= 2 && em.nrow > em.ncol + 8, "multivariateNormality: need at least two columns and more than ncol + 8 rows");
    auto sg = sigma.nrow == 0 ? (t(em) * em) / cast(double) em.nrow : sigma;

    enforce!GretlException(sg.nrow == em.ncol && sg.ncol == em.ncol, "multivariateNormality: sigma must be " ~ to!string(em.ncol) ~ " x " ~ to!string(em.ncol));
    auto pr = Printer(0);
    const err = multivariate_normality_test(em.ptr, sg.ptr, Gretlopt.OPT_NONE, pr.prn);

    enforce!GretlException(err == 0, "multivariateNormality: " ~ lastError());
    return StatTest(get_last_test_statistic(), get_last_pvalue());
}

/* ---- two-sample tests ---- */

/* The result of a two-sample test. */
struct TwoSampleTest {
    string test;
    double statistic;
    double pvalue;              // two-sided
    double df;                  // degrees of freedom (the numerator's, for an F test)
    double df2 = double.nan;    // the denominator's, for an F test
    double difference = double.nan;     // for a test of means: the difference of the sample means, and its standard error
    double se = double.nan;
    size_t n1, n2;              // valid observations in each sample (the paired differences: n1)
}

/* How the means of two samples are compared. */
enum MeansTest {
    equalVariances,         // the pooled-variance t test
    unequalVariances,       // Welch's t test with the Satterthwaite degrees of freedom (the library rounds them down to an integer)
    paired                  // the t test of the mean of the differences of paired observations (same length; a pair with a missing value is left out)
}

private double[] validOf(const(double)[] v) {
    double[] a;

    foreach (e; v) if (!isNaN(e)) a ~= e;
    return a;
}

/* meansTest(x, y, kind): the test of the null that two populations have the same mean. The samples may have different lengths (not for `paired`). */
TwoSampleTest meansTest(S, T)(S x, T y, MeansTest kind = MeansTest.equalVariances) @trusted if (isSeries!S && isSeries!T) {
    auto xa = asDoubles(x), ya = asDoubles(y);

    if (kind == MeansTest.paired) {
        enforce!GretlException(xa.length == ya.length, "meansTest: paired samples need the same length");
    }
    TwoSampleTest r;
    double[] xv, yv;

    r.test = ["Equality of means, equal variances", "Equality of means, unequal variances", "Paired differences"][kind];
    if (kind == MeansTest.paired) {
        foreach (i; 0 .. xa.length) {
            if (!isNaN(xa[i]) && !isNaN(ya[i])) {
                xv ~= xa[i];
                yv ~= ya[i];
            }
        }
    } else {
        xv = validOf(xa);
        yv = validOf(ya);
    }
    enforce!GretlException(xv.length >= 2 && yv.length >= 2, "meansTest: each sample needs at least two valid observations");
    const n = xv.length > yv.length ? xv.length : yv.length;
    auto xp = xv.dup, yp = yv.dup;

    xp.length = n; yp.length = n;
    foreach (i; xv.length .. n) xp[i] = double.nan;
    foreach (i; yv.length .. n) yp[i] = double.nan;
    auto ds = makeDataset([xp, yp], ["const", "x", "y"]);
    auto list = makeList([1, 2]);
    auto opt = Gretlopt.OPT_Q;

    if (kind == MeansTest.unequalVariances) opt |= Gretlopt.OPT_O;
    if (kind == MeansTest.paired) opt |= Gretlopt.OPT_P;
    enforce!GretlException(push_program_state() == 0, "could not save the library settings");
    scope (exit) {
        pop_program_state();
        import core.stdc.stdlib: free;
        free(list);
    }
    auto pr = Printer(0);
    const err = means_test(list, ds.ptr, opt, pr.prn);

    enforce!GretlException(err == 0, "meansTest: " ~ lastError());
    r.statistic = get_last_test_statistic();
    r.pvalue = get_last_pvalue();
    r.n1 = xv.length;
    r.n2 = kind == MeansTest.paired ? 0 : yv.length;
    if (kind == MeansTest.paired) {
        double[] d;

        foreach (i; 0 .. xv.length) d ~= xv[i] - yv[i];
        r.difference = meanA(d);
        r.df = cast(double) xv.length - 1;
    } else {
        r.difference = meanA(xv) - meanA(yv);
        const v1 = varianceA(xv, false), v2 = varianceA(yv, false);

        r.df = kind == MeansTest.unequalVariances ? satterthwaite_df(v1, cast(int) xv.length, v2, cast(int) yv.length) : cast(double) (xv.length + yv.length - 2);
    }
    r.se = r.difference / r.statistic;
    return r;
}

/* variancesTest(x, y): the F test of the null that two populations have the same variance, F = var(x) / var(y) with n1 - 1 and n2 - 1 degrees of freedom. `pvalue` is the
   two-sided p-value, 2 min(P(F' > F), P(F' < F)); `pvalueUpper` is the one-sided P(F' > F) that the library itself records (see bugs.md #124). */
struct VarianceTest {
    TwoSampleTest test;
    double pvalueUpper;
}

VarianceTest variancesTest(S, T)(S x, T y) @trusted if (isSeries!S && isSeries!T) {
    auto xv = validOf(asDoubles(x)), yv = validOf(asDoubles(y));

    enforce!GretlException(xv.length >= 2 && yv.length >= 2, "variancesTest: each sample needs at least two valid observations");
    const n = xv.length > yv.length ? xv.length : yv.length;
    auto xp = xv.dup, yp = yv.dup;

    xp.length = n; yp.length = n;
    foreach (i; xv.length .. n) xp[i] = double.nan;
    foreach (i; yv.length .. n) yp[i] = double.nan;
    auto ds = makeDataset([xp, yp], ["const", "x", "y"]);
    auto list = makeList([1, 2]);

    enforce!GretlException(push_program_state() == 0, "could not save the library settings");
    scope (exit) {
        pop_program_state();
        import core.stdc.stdlib: free;
        free(list);
    }
    auto pr = Printer(0);
    const err = vars_test(list, ds.ptr, Gretlopt.OPT_Q, pr.prn);

    enforce!GretlException(err == 0, "variancesTest: " ~ lastError());
    VarianceTest r;

    r.test.test = "Equality of variances (F)";
    r.test.statistic = get_last_test_statistic();
    r.pvalueUpper = get_last_pvalue();
    r.test.df = cast(double) xv.length - 1;
    r.test.df2 = cast(double) yv.length - 1;
    r.test.pvalue = snedecor_pvalue_2(r.test.df, r.test.df2, r.test.statistic);
    r.test.n1 = xv.length;
    r.test.n2 = yv.length;
    return r;
}

/* ---- rank and sign tests ---- */

/* The sign test of the null that two paired series have the same median: W is the number of pairs with x > y among those that differ. */
struct SignTest {
    size_t n;                   // pairs with valid values that differ
    size_t w;                   // of them, the number with x > y
    double pvalueUpper;         // P(W' >= w) with W' binomial(n, 1/2): the one the library records
    double pvalueLower;         // P(W' <= w)
    double pvalue;              // two-sided: min(1, 2 min(upper, lower))
    string report;
}

/* signTest(x, y): the sign test for series of the same length (pairs where either value is missing, or the values are equal, are left out). */
SignTest signTest(S, T)(S x, T y) @trusted if (isSeries!S && isSeries!T) {
    auto xa = asDoubles(x), ya = asDoubles(y);

    enforce!GretlException(xa.length == ya.length, "signTest: the series have " ~ to!string(xa.length) ~ " and " ~ to!string(ya.length) ~ " observations");
    auto pr = Printer(0);
    SignTest r;

    const res = runDiffTest(xa, ya, NonparametricKind.sign, "signTest", pr);

    r.report = pr.text;
    r.w = cast(size_t) res[0];
    r.pvalueUpper = res[1];
    foreach (i; 0 .. xa.length) {
        if (!isNaN(xa[i]) && !isNaN(ya[i]) && xa[i] != ya[i]) r.n++;
    }
    r.pvalueLower = binomial_cdf(0.5, cast(int) r.n, cast(int) r.w);
    const m = r.pvalueUpper < r.pvalueLower ? r.pvalueUpper : r.pvalueLower;

    r.pvalue = 2 * m > 1 ? 1 : 2 * m;
    return r;
}

/* The Wilcoxon signed-rank test of the null that the median of the differences x - y is zero (Cureton's normal approximation, with ties and zero differences). */
struct SignedRankTest {
    size_t pairs;               // pairs with valid values
    size_t n;                   // of them, those with a non-zero difference
    size_t zeroDifferences;
    size_t nonZeroTies;         // groups of equal absolute differences
    double wPlus = 0, wMinus = 0;       // the sums of the ranks of the positive and the negative differences
    double expected, variance;  // of W+ under the null (NaN when n is 8 or less)
    double statistic = double.nan;      // z with the continuity correction; NaN when n <= 8 (the normal approximation is not used: see `report`, which gives 5% critical values for n 6 to 8)
    double pvalue = double.nan;         // two-sided
    double pvalueOneSided = double.nan; // the tail beyond z: the one the library records
    string report;
}

/* signedRankTest(x, y): the Wilcoxon signed-rank test for paired series of the same length. Zero differences take the lowest ranks. If every difference is zero the
   library does not complain: n is 0, the statistic NaN and the report says that H0 is not rejected. */
SignedRankTest signedRankTest(S, T)(S x, T y) @trusted if (isSeries!S && isSeries!T) {
    auto xa = asDoubles(x), ya = asDoubles(y);

    enforce!GretlException(xa.length == ya.length, "signedRankTest: the series have " ~ to!string(xa.length) ~ " and " ~ to!string(ya.length) ~ " observations");
    auto pr = Printer(0);
    SignedRankTest r;

    const res = runDiffTest(xa, ya, NonparametricKind.signedRank, "signedRankTest", pr);

    r.report = pr.text;
    // the sums of ranks, by the library's rules: zero differences first, ties get the average rank
    double[] ad;
    bool[] pos;

    foreach (i; 0 .. xa.length) {
        if (isNaN(xa[i]) || isNaN(ya[i])) continue;
        r.pairs++;
        const d = xa[i] - ya[i];

        if (d != 0) {
            ad ~= d < 0 ? -d : d;
            pos ~= d > 0;
        }
    }
    r.n = ad.length;
    r.zeroDifferences = r.pairs - r.n;
    auto idx = new size_t[ad.length];

    foreach (i, ref v; idx) v = i;
    import std.algorithm.sorting: stdSort = sort;
    stdSort!((a, b) => ad[a] < ad[b])(idx);
    double tsum = 0;

    for (size_t i = 0; i < idx.length;) {
        size_t j = i;

        while (j + 1 < idx.length && ad[idx[j + 1]] == ad[idx[i]]) j++;
        const m = j - i;                                           // extra members of the tie group
        const rank = (r.zeroDifferences + i + 1 + r.zeroDifferences + j + 1) / 2.0;

        foreach (k; i .. j + 1) {
            if (pos[idx[k]]) r.wPlus += rank; else r.wMinus += rank;
        }
        if (m > 0) {
            tsum += (cast(double) m) ^^ 3 - m;
            r.nonZeroTies++;
        }
        i = j + 1;
    }
    if (r.n > 8) {
        const N = cast(double) r.pairs, Z = cast(double) r.zeroDifferences;

        r.expected = (N * (N + 1) - Z * (Z + 1)) / 4.0;
        r.variance = (N * (N + 1) * (2 * N + 1) - Z * (Z + 1) * (2 * Z + 1) - tsum / 2) / 24.0;
        r.statistic = res[0];
        r.pvalueOneSided = res[1];
        r.pvalue = isNaN(r.pvalueOneSided) ? double.nan : (2 * r.pvalueOneSided > 1 ? 1 : 2 * r.pvalueOneSided);
    } else {
        r.expected = double.nan;
        r.variance = double.nan;
    }
    return r;
}

/* Critical values of the rank sum of the smaller sample (the library's table, for smaller samples of 4 to 9 and larger ones of up to 12 observations): reject at the level
   if the rank sum of the smaller sample is at most the lower or at least the upper value. NaN where the table has no value. */
struct RankSumCritical {
    bool available;
    double lower1 = double.nan, lower5 = double.nan, lower10 = double.nan;      // lower tail, 1%, 5%, 10%
    double upper10 = double.nan, upper5 = double.nan, upper1 = double.nan;      // upper tail, 10%, 5%, 1%
}

/* The Wilcoxon rank-sum (Mann-Whitney) test of the null that two samples come from populations with the same median. */
struct RankSumTest {
    size_t nx, ny;              // valid observations in each sample
    double rankSumX = 0, rankSumY = 0;  // sums of the (mid)ranks of each sample in the pooled ranking: they add to (nx + ny)(nx + ny + 1) / 2
    double rankSumSmaller = double.nan; // the rank sum of the smaller sample (x if the sizes are equal)
    double statistic = double.nan;      // z for the SMALLER sample's rank sum (x if the sizes are equal); NaN unless both samples have 10 or more observations
    double pvalue = double.nan;         // two-sided
    double pvalueOneSided = double.nan; // the tail beyond z: the one the library records
    RankSumCritical critical;           // for small samples (smaller 4 to 9, larger at most 12): the critical values of rankSumSmaller
    string report;
}

/* rankSumTest(x, y): the Wilcoxon rank-sum test for two independent samples (any lengths, missing values left out). The normal approximation is used when both samples have
   at least 10 observations. For smaller samples there is no p-value (`statistic` and the p-values are NaN): when the smaller has 4 to 9 observations and the larger at most 12,
   `critical` gives the critical values of the table for the rank sum of the smaller sample; below that the sample is too small for a test. */
RankSumTest rankSumTest(S, T)(S x, T y) @trusted if (isSeries!S && isSeries!T) {
    auto xv = validOf(asDoubles(x)), yv = validOf(asDoubles(y));

    enforce!GretlException(xv.length >= 1 && yv.length >= 1, "rankSumTest: each sample needs at least one valid observation");
    const n = xv.length > yv.length ? xv.length : yv.length;
    auto xp = xv.dup, yp = yv.dup;

    xp.length = n; yp.length = n;
    foreach (i; xv.length .. n) xp[i] = double.nan;
    foreach (i; yv.length .. n) yp[i] = double.nan;
    auto pr = Printer(0);

    const res = runDiffTest(xp, yp, NonparametricKind.rankSum, "rankSumTest", pr);
    RankSumTest r;

    r.report = pr.text;
    r.nx = xv.length;
    r.ny = yv.length;
    // midranks of the pooled sample
    auto all = xv ~ yv;
    auto idx = new size_t[all.length];

    foreach (i, ref v; idx) v = i;
    import std.algorithm.sorting: stdSort = sort;
    stdSort!((a, b) => all[a] < all[b])(idx);
    for (size_t i = 0; i < idx.length;) {
        size_t j = i;

        while (j + 1 < idx.length && all[idx[j + 1]] == all[idx[i]]) j++;
        const rank = (i + 1 + j + 1) / 2.0;

        foreach (k; i .. j + 1) {
            if (idx[k] < xv.length) r.rankSumX += rank; else r.rankSumY += rank;
        }
        i = j + 1;
    }
    const smallN = xv.length <= yv.length ? xv.length : yv.length, largeN = xv.length <= yv.length ? yv.length : xv.length;

    r.rankSumSmaller = xv.length <= yv.length ? r.rankSumX : r.rankSumY;
    if (smallN >= 4 && smallN <= 9 && largeN >= smallN && largeN <= 12) {
        // the library prints the table row ("lower tail: 1% a, 5% b, 10% c"); the upper values follow by symmetry (the printed upper tail of the first row is wrong, bugs.md #128)
        double[3] lower = double.nan;
        const marker = "lower tail:";
        const at = r.report.indexOf(marker);

        if (at >= 0) {
            auto line = r.report[at + marker.length .. $];
            const eol = line.indexOf('\n');

            if (eol >= 0) line = line[0 .. eol];
            foreach (part; line.split(",")) {
                auto bits = part.strip().split("%");

                if (bits.length == 2) {
                    const level = to!int(bits[0].strip());
                    const value = to!int(bits[1].strip());

                    if (level == 1) lower[0] = value;
                    else if (level == 5) lower[1] = value;
                    else if (level == 10) lower[2] = value;
                }
            }
        }
        const total = cast(double) smallN * (smallN + largeN + 1);

        r.critical.available = !isNaN(lower[1]);
        r.critical.lower1 = lower[0];
        r.critical.lower5 = lower[1];
        r.critical.lower10 = lower[2];
        r.critical.upper1 = isNaN(lower[0]) ? double.nan : total - lower[0];
        r.critical.upper5 = isNaN(lower[1]) ? double.nan : total - lower[1];
        r.critical.upper10 = isNaN(lower[2]) ? double.nan : total - lower[2];
    }
    if (xv.length >= 10 && yv.length >= 10) {
        r.statistic = res[0];
        r.pvalueOneSided = res[1];
        r.pvalue = isNaN(r.pvalueOneSided) ? double.nan : (2 * r.pvalueOneSided > 1 ? 1 : 2 * r.pvalueOneSided);
    }
    return r;
}

/* The three tests of gretl.nonparam, which diff_test (the `difftest` command) calls with arrays that it has compacted to the valid values but that the tests index over the
   whole sample range: with a missing value that reads beyond the arrays and misaligns the pairs (bugs.md #125). Here they are called directly with full-length arrays in which a
   missing value is NaN, as they expect. */
private enum NonparametricKind { sign, signedRank, rankSum }

private double[2] runDiffTest(const(double)[] x, const(double)[] y, NonparametricKind kind, string who, ref Printer pr) @trusted {
    enforce!GretlException(x.length >= 2 && x.length == y.length, who ~ ": too few observations");
    auto ds = makeDataset([x, y], ["const", "x", "y"]);
    double[2] result = double.nan;

    enforce!GretlException(push_program_state() == 0, "could not save the library settings");
    scope (exit) pop_program_state();
    int err;

    final switch (kind) {
        case NonparametricKind.sign: err = sign_test(x.ptr, y.ptr, 1, 2, ds.ptr, result.ptr, Gretlopt.OPT_NONE, pr.prn); break;
        case NonparametricKind.signedRank: err = signed_rank_test(x.ptr, y.ptr, 1, 2, ds.ptr, result.ptr, Gretlopt.OPT_NONE, pr.prn); break;
        case NonparametricKind.rankSum: err = rank_sum_test(x.ptr, y.ptr, 1, 2, ds.ptr, result.ptr, Gretlopt.OPT_NONE, pr.prn); break;
    }
    enforce!GretlException(err == 0, who ~ ": " ~ (err == 0 ? "" : (lastError().length ? lastError() : "no valid observations")));
    return result;
}

/// ditto, for series of a Dataset
SignTest signTest(Dataset d, string first, string second) @trusted {
    return signTest(d[first], d[second]);
}

/// ditto
SignedRankTest signedRankTest(Dataset d, string first, string second) @trusted {
    return signedRankTest(d[first], d[second]);
}

/// ditto
RankSumTest rankSumTest(Dataset d, string first, string second) @trusted {
    return rankSumTest(d[first], d[second]);
}

/* ---- distances ---- */

/* The distance between observations (rows). */
enum DistanceMetric {
    euclidean,
    manhattan,
    hamming,        // the number (the proportion in the library) of coordinates that differ
    chebyshev,      // the largest absolute difference
    cosine,         // 1 minus the cosine of the angle between the two vectors
    mahalanobis     // the distance with the covariance matrix of the rows of X
}

/* rowDistances(X, Y, metric): the distances between the rows of X (m x k) and the rows of Y (p x k): an m x p matrix. */
GretlMatrix rowDistances(MX, MY)(MX x, MY y, DistanceMetric metric = DistanceMetric.euclidean) @trusted if (isMatrixLike!MX && isMatrixLike!MY) {
    auto xm = cbind(x), ym = cbind(y);

    enforce!GretlException(xm.ncol == ym.ncol, "rowDistances: X and Y must have the same number of columns");
    requireFinite(xm, "rowDistances: X");
    requireFinite(ym, "rowDistances: Y");
    return libDistance(xm, ym.ptr, metric);
}

/* rowDistanceMatrix(X, metric): the m x m symmetric matrix of the distances between the rows of X (zeros on the diagonal). */
GretlMatrix rowDistanceMatrix(MX)(MX x, DistanceMetric metric = DistanceMetric.euclidean) @trusted if (isMatrixLike!MX) {
    auto xm = cbind(x);

    requireFinite(xm, "rowDistanceMatrix: X");
    enforce!GretlException(xm.nrow >= 2, "rowDistanceMatrix: need at least two rows");
    auto v = libDistance(xm, null, metric);
    const m = xm.nrow;
    auto r = zerosMatrix(m, m);
    size_t k = 0;

    foreach (j; 0 .. m) {
        foreach (i; j + 1 .. m) {
            r.ptr.val[j * m + i] = v.ptr.val[k];
            r.ptr.val[i * m + j] = v.ptr.val[k];
            k++;
        }
    }
    return r;
}

private GretlMatrix zerosMatrix(size_t r, size_t c) @trusted {
    auto m = fresh(r, c);

    foreach (i; 0 .. r * c) m.ptr.val[i] = 0;
    return m;
}

private GretlMatrix libDistance(GretlMatrix xm, GretlMatrixRaw * yraw, DistanceMetric metric) @trusted {
    int err = 0;
    auto mraw = libDistanceFn(xm.ptr, yraw, ["euclidean", "manhattan", "hamming", "chebyshev", "cosine", "mahalanobis"][metric].ptr, &err);

    enforce!GretlException(err == 0 && mraw !is null, "rowDistances: " ~ lastError());
    auto r = fresh(mraw.rows, mraw.cols);

    foreach (i; 0 .. cast(size_t) mraw.rows * mraw.cols) r.ptr.val[i] = mraw.val[i];
    gretl_matrix_free(mraw);
    return r;
}

/* mahalanobis(X): the Mahalanobis distance of each row of X (n x k, k >= 1) from the centroid of the rows, as an n x 1 matrix: sqrt((x - mean)' S^-1 (x - mean)) with S the
   sample covariance matrix. A row with a missing value gets NaN (and is not used for the centroid or S). */
GretlMatrix mahalanobis(MX)(MX x) @trusted if (isMatrixLike!MX) {
    auto xm = cbind(x);
    const k = xm.ncol, n = xm.nrow;

    enforce!GretlException(k >= 1 && n > k, "mahalanobis: need more rows than columns");
    auto ds = makeDataset(columnsOf(xm), ["const"] ~ numbered("x", k));
    int[] spec;

    foreach (j; 0 .. k) spec ~= cast(int) (1 + j);
    auto list = makeList(spec);

    scope (exit) {
        import core.stdc.stdlib: free;
        free(list);
    }
    int err = 0;
    auto md = get_mahal_distances(list, ds.ptr, Gretlopt.OPT_NONE, null, &err);

    enforce!GretlException(err == 0 && md !is null, "mahalanobis: " ~ lastError());
    scope (exit) free_mahal_dist(md);
    auto r = fresh(n, 1);
    const d = mahal_dist_get_distances(md);

    foreach (i; 0 .. n) r.ptr.val[i] = d[i];
    return r;
}

/* ---- cross-tabulation ---- */

/* A cross-tabulation of two discrete series: counts[i, j] is the number of observations with the i-th value of the first series and the j-th of the second. */
struct CrossTab {
    GretlMatrix counts;         // rows x cols
    GretlMatrix rowValues;      // rows x 1: the distinct values of the first series, increasing
    GretlMatrix colValues;      // cols x 1
    GretlMatrix rowTotals;      // rows x 1
    GretlMatrix colTotals;      // cols x 1
    int nobs;                   // observations counted (both values present)
    double pearson;             // Pearson's chi-square statistic of independence: the sum of (O - E)^2 / E
    int df;                     // (rows - 1)(cols - 1)
    double pvalue;              // chi-square(df) p-value (NaN if df is 0)
}

/* crosstab(x, y): the table of counts of the pairs of values of two discrete series (integer valued, or taking few distinct values) of the same length, with the Pearson
   chi-square test of independence. Observations with a missing value in either series are left out. */
CrossTab crosstab(S, T)(S x, T y) @trusted if (isSeries!S && isSeries!T) {
    auto xa = asDoubles(x), ya = asDoubles(y);

    enforce!GretlException(xa.length == ya.length, "crosstab: the series have " ~ to!string(xa.length) ~ " and " ~ to!string(ya.length) ~ " observations");
    enforce!GretlException(xa.length >= 2, "crosstab: too few observations");
    auto ds = makeDataset([xa, ya], ["const", "x", "y"]);
    auto list = makeList([1, 2]);

    scope (exit) {
        import core.stdc.stdlib: free;
        free(list);
    }
    int err = 0;
    auto xt = single_crosstab(list, ds.ptr, Gretlopt.OPT_NONE, null, &err);

    enforce!GretlException(err == 0 && xt !is null, "crosstab: " ~ lastError());
    scope (exit) free_xtab(xt);
    CrossTab r;
    const nr = xt.rows, nc = xt.cols;

    r.counts = fresh(nr, nc);
    r.rowValues = fresh(nr, 1);
    r.colValues = fresh(nc, 1);
    r.rowTotals = fresh(nr, 1);
    r.colTotals = fresh(nc, 1);
    r.nobs = xt.n;
    foreach (i; 0 .. nr) {
        r.rowValues.ptr.val[i] = xt.rval[i];
        r.rowTotals.ptr.val[i] = xt.rtotal[i];
        foreach (j; 0 .. nc) r.counts.ptr.val[j * nr + i] = xt.f[i][j];
    }
    foreach (j; 0 .. nc) {
        r.colValues.ptr.val[j] = xt.cval[j];
        r.colTotals.ptr.val[j] = xt.ctotal[j];
    }
    r.df = (nr - 1) * (nc - 1);
    r.pearson = 0;
    foreach (i; 0 .. nr) {
        foreach (j; 0 .. nc) {
            const e = r.rowTotals[i, 0] * r.colTotals[j, 0] / r.nobs, o = r.counts[i, j];

            if (e > 0) r.pearson += (o - e) * (o - e) / e;
        }
    }
    r.pvalue = r.df > 0 ? chisq_cdf_comp(r.df, r.pearson) : double.nan;
    return r;
}

/* ---- on the series of a Dataset ---- */

/* The functions for one series take a series of a Dataset as it is: mean(d["gdp"]), quantile(d["gdp"], 0.9), frequencies(d["gdp"]). Those that take several
   series have forms with the names: */

private GretlMatrix columnsOfDataset(Dataset d, string[] names, string who) @trusted {
    enforce!GretlException(names.length > 0, who ~ ": no series named");
    double[][] cols;

    foreach (nm; names) cols ~= d[nm].dup;
    return matrixFrom(cols);
}

/// ditto
GretlMatrix covariance(Dataset d, string[] names) @trusted {
    return covariance(columnsOfDataset(d, names, "covariance"));
}

/// ditto
GretlMatrix correlation(Dataset d, string[] names) @trusted {
    return correlation(columnsOfDataset(d, names, "correlation"));
}

/// ditto
GretlMatrix mahalanobis(Dataset d, string[] names) @trusted {
    return mahalanobis(columnsOfDataset(d, names, "mahalanobis"));
}

/// ditto: the cross-tabulation of two series
CrossTab crosstab(Dataset d, string rows, string cols) @trusted {
    return crosstab(d[rows], d[cols]);
}

/// ditto
TwoSampleTest meansTest(Dataset d, string first, string second, MeansTest kind = MeansTest.equalVariances) @trusted {
    return meansTest(d[first], d[second], kind);
}

/// ditto
VarianceTest variancesTest(Dataset d, string first, string second) @trusted {
    return variancesTest(d[first], d[second]);
}
