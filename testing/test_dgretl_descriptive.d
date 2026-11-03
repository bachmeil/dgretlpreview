/* Tests for dgretl/descriptive.d. */

import core.stdc.stdlib;
import std.conv: to;
import std.math: abs, isNaN, sqrt;
import std.stdio: writefln, writeln;
import dgretl;
import dgretl.descriptive;

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

ulong state = 777;
double urand() {
    state = state * 6364136223846793005UL + 1442695040888963407UL;
    return ((state >> 11) / cast(double) (1UL << 53)) - 0.5;
}

void main() {
    enum n = 50;
    double[] x = new double[n], y = new double[n], z = new double[n];

    foreach (i; 0 .. n) {
        x[i] = 10 * urand() + 3;
        y[i] = 0.5 * x[i] + 4 * urand();
        z[i] = urand() * urand();
    }
    double mu = 0;
    foreach (v; x) mu += v;
    mu /= n;
    double ss = 0, s3 = 0, s4 = 0;
    foreach (v; x) { ss += (v - mu) ^^ 2; }
    const popvar = ss / n;
    foreach (v; x) { s3 += ((v - mu) / sqrt(popvar)) ^^ 3; s4 += ((v - mu) / sqrt(popvar)) ^^ 4; }

    check(close(mean(x), mu) && close(total(x), mu * n), "mean and total of an array");
    check(close(variance(x), ss / (n - 1)) && close(variance(x, true), popvar) && close(stddev(x), sqrt(ss / (n - 1))) && close(stddev(x, true), sqrt(popvar)), "variance and stddev, sample and population");
    check(close(sst(x), ss), "sst");
    check(close(skewness(x), s3 / n) && close(kurtosis(x), s4 / n - 3), "skewness and the EXCESS kurtosis");
    double lo = x[0], hi = x[0];
    foreach (v; x) { if (v < lo) lo = v; if (v > hi) hi = v; }
    check(minimum(x) == lo && maximum(x) == hi && validCount(x) == n, "minimum, maximum, count");

    // quantiles against the order statistics
    auto sorted = x.dup;
    import std.algorithm.sorting: stdSort = sort;
    stdSort(sorted);
    check(close(median(x), (sorted[24] + sorted[25]) / 2), "median of an even number of values");
    check(close(quantile(x, 0), sorted[0]) && close(quantile(x, 1), sorted[n - 1]), "quantile at 0 and 1 are the minimum and the maximum");
    double qfor(double p, double offset, double slope) {          // position (n + offset-ish): h = slope * p + offset, 1-based
        const hh = slope * p + offset - 1;
        const lo = cast(size_t) hh;
        return sorted[lo] + (hh - lo) * (sorted[lo + 1] - sorted[lo]);
    }
    check(close(quantile(x, 0.9), qfor(0.9, 0, n + 1.0)) && close(quantile(x, 0.9, QuantileType.q6), qfor(0.9, 0, n + 1.0)), "quantile: type 6 by default, (n + 1) p");
    check(close(quantile(x, 0.9, QuantileType.q7), qfor(0.9, 1, n - 1.0)) && close(quantile(x, 0.3, QuantileType.q8), qfor(0.3, 1.0 / 3, n + 1.0 / 3)), "quantile: types 7 and 8");
    check(isNaN(quantile(x, 0.001)) && isNaN(quantile(x[0 .. 3], 0.99)) && isNaN(quantile(new double[0], 0.5)), "quantile: NaN when the sample is too small for p");
    check(close(quantile(x, 0.9), quantile(x, 0.9, QuantileType.q6)), "quantile: the setting is restored afterwards");
    auto qs = quantiles(x, [0.25, 0.5, 0.75]);
    auto f5 = fiveNumbers(x);
    check(qs.length == 3 && close(qs[1], median(x)) && f5.min == lo && f5.max == hi && close(f5.q1, qs[0]) && close(f5.q3, qs[2]) && f5.q1 < f5.median && f5.median < f5.q3, "quantiles and fiveNumbers");
    check(throws({ quantile(x, 1.5); }) && throws({ quantile(x, -0.1); }), "quantile: p outside [0, 1] throws");

    // Gini of a constant-sum two-point distribution: (0, 0, 1, 1) has Gini 0.5
    check(close(gini([0.0, 0.0, 1.0, 1.0]), 0.5, 1e-9) && close(gini([1.0, 1.0, 1.0, 1.0]), 0, 1e-12), "gini");

    // missing values are skipped
    auto xm = x.dup;
    xm[3] = double.nan;
    xm[10] = double.nan;
    double m2 = 0;
    foreach (i, v; xm) if (!isNaN(v)) m2 += v;
    m2 /= n - 2;
    check(close(mean(xm), m2) && validCount(xm) == n - 2 && !isNaN(variance(xm)) && !isNaN(median(xm)), "missing values are skipped");
    check(isNaN(mean(new double[0])) && isNaN(mean([double.nan, double.nan])) && isNaN(median([double.nan])) && isNaN(minimum([double.nan])) && variance([double.nan, 2.0]) == 0, "no valid observation gives NaN");
    check(variance([1.0]) == 0 && stddev([1.0]) == 0, "one valid observation: the library gives variance 0");

    // int arrays, and lines of a matrix
    int[] ii = [3, 1, 4, 1, 5, 9, 2, 6];
    check(close(mean(ii), 31.0 / 8) && close(median(ii), 3.5) && maximum(ii) == 9, "an int array");
    auto M = cbind(x, y, z);
    check(close(mean(Col(M, 1)), mean(y)) && close(median(Row(M, 0)), median([x[0], y[0], z[0]])), "a column and a row of a matrix");
    check(close(variance(Col(M, 2)), variance(z)), "variance of a column of a matrix");

    // matrices: column by column
    auto mm = mean(M);
    check(mm.nrow == 1 && mm.ncol == 3 && close(mm[0, 0], mean(x)) && close(mm[0, 1], mean(y)) && close(mm[0, 2], mean(z)), "mean of a matrix: one value per column");
    auto sd = stddev(M, true), sk = skewness(M), ku = kurtosis(M), md = median(M), mn = minimum(M), mx = maximum(M), qq = quantile(M, 0.75), tt = total(M), vv = variance(M);
    check(close(sd[0, 1], stddev(y, true)) && close(sk[0, 2], skewness(z)) && close(ku[0, 0], kurtosis(x)) && close(md[0, 1], median(y)), "stddev, skewness, kurtosis and median of a matrix");
    check(mn[0, 0] == minimum(x) && mx[0, 2] == maximum(z) && close(qq[0, 1], quantile(y, 0.75)) && close(tt[0, 2], total(z)) && close(vv[0, 0], variance(x)), "minimum, maximum, quantile, total, variance of a matrix");
    auto sub = M.sub[10 .. 30, 0 .. 2];
    check(mean(sub).ncol == 2 && close(mean(sub)[0, 1], mean(y[10 .. 30])), "a submatrix");
    double[] px = new double[n], pz = new double[n];
    foreach (i; 0 .. n) { px[i] = abs(x[i]); pz[i] = abs(z[i]); }
    auto gm = gini(cbind(px, pz)), cm = validCount(M);
    check(gm.ncol == 2 && close(gm[0, 1], gini(pz)) && cm[0, 0] == n, "gini and count of a matrix");

    // long run variance with lag 0 is the population variance
    double sx2 = 0;
    foreach (v; x) sx2 += v * v;
    check(close(longRunVariance(x, 0), popvar, 1e-9) && longRunVariance(x, 4) > 0, "long-run variance: lag 0 is the population variance");
    check(!isNaN(longRunVariance(x)) && close(longRunVariance(x, 0, 0.0), sx2 / n, 1e-9), "long-run variance: automatic lag and a given mean");

    // covariance and correlation
    double my = mean(y), cxy = 0, sy = 0;
    foreach (i; 0 .. n) { cxy += (x[i] - mu) * (y[i] - my); sy += (y[i] - my) ^^ 2; }
    check(close(covariance(x, y), cxy / (n - 1)) && close(correlation(x, y), cxy / sqrt(ss * sy)), "covariance and correlation of two series");
    check(close(correlation(x, x), 1.0) && close(covariance(x, x), variance(x)), "a series with itself");
    auto cv = covariance(M), cr = correlation(M);
    check(cv.nrow == 3 && cv.ncol == 3 && close(cv[0, 1], covariance(x, y)) && close(cv[1, 0], cv[0, 1]) && close(cv[2, 2], variance(z)), "covariance matrix");
    check(close(cr[0, 0], 1) && close(cr[0, 1], correlation(x, y)) && close(cr[2, 1], cr[1, 2]), "correlation matrix");
    auto ym = y.dup;
    ym[5] = double.nan;
    auto crm = correlation(cbind(x, ym));
    double[] xs2, ys2;
    foreach (i; 0 .. n) if (!isNaN(ym[i])) { xs2 ~= x[i]; ys2 ~= ym[i]; }
    check(close(crm[0, 1], correlation(xs2, ys2)), "correlation matrix: pairwise deletion of missing values");
    check(throws({ covariance(x, y[0 .. 10]); }) && throws({ correlation(x, y[0 .. 10]); }), "series of different length throw");

    stage2(x, y, z);
    stage3(x, y);
    stage4(x, y, z);
    stage5(x, y, z);
    stage6();

    writefln("\n%d tests, %d failed", nTests, nFailed);
    if (nFailed) throw new Exception("failed");
}


void stage2(double[] x, double[] y, double[] z) {
    enum n = 50;
    auto M = cbind(x, y, z);
    auto sm = summary(M, ["a", "b", "c"]);

    check(sm.names == ["a", "b", "c"] && sm.mean.ncol == 3 && sm.valid[0, 0] == n && sm.missing[0, 2] == 0, "summary: shape, names and counts");
    bool ok = true;
    foreach (j, col; [x, y, z]) {
        ok = ok && close(sm.mean[0, j], mean(col)) && close(sm.median[0, j], median(col)) && close(sm.sd[0, j], stddev(col)) && close(sm.min[0, j], minimum(col))
                && close(sm.max[0, j], maximum(col)) && close(sm.skewness[0, j], skewness(col)) && close(sm.excessKurtosis[0, j], kurtosis(col))
                && close(sm.p5[0, j], quantile(col, 0.05)) && close(sm.p95[0, j], quantile(col, 0.95)) && close(sm.iqr[0, j], quantile(col, 0.75) - quantile(col, 0.25))
                && close(sm.cv[0, j], abs(stddev(col) / mean(col)));
    }
    check(ok, "summary: every statistic equals the function of the same name");
    auto xm = x.dup;
    xm[2] = double.nan; xm[9] = double.nan;
    auto sx = summary(xm);
    check(sx.names == ["x"] && sx.valid[0, 0] == n - 2 && sx.missing[0, 0] == 2 && close(sx.mean[0, 0], mean(xm)), "summary of a single series with missing values");
    double[] allnan = new double[n];
    allnan[] = double.nan;
    auto sn = summary(cbind(x, allnan), ["a", "b"]);
    check(sn.valid[0, 0] == n && isNaN(sn.mean[0, 1]) && isNaN(sn.valid[0, 1]), "summary: a series without valid observations has NaN everywhere");
    // weights: integer weights are replication
    double[] w = new double[n];
    foreach (i; 0 .. n) w[i] = i % 3 == 0 ? 2.0 : 1.0;
    double sw = 0, swx = 0;
    foreach (i; 0 .. n) { sw += w[i]; swx += w[i] * x[i]; }
    auto sw1 = summary(x, null, w);
    check(close(sw1.mean[0, 0], swx / sw, 1e-9), "summary: the weighted mean");
    w[4] = 0;
    auto sw2 = summary(x, null, w);
    check(sw2.valid[0, 0] == n - 1, "summary: observations with weight 0 are left out");
    check(throws({ summary(M, ["a", "b"]); }) && throws({ summary(x, null, w[0 .. 5]); }) && throws({ summary(new double[0]); }), "summary: bad arguments throw");
    check(sm.toString().length > 0, "summary: a printable table");
    {
        auto d = Dataset(3, n);
        import core.stdc.string: strcpy;
        strcpy(d.ptr.varname[1], "p"); strcpy(d.ptr.varname[2], "q");
        foreach (i; 0 .. n) { d.ptr.Z[1][i] = x[i]; d.ptr.Z[2][i] = y[i]; }
        auto sd2 = summary(d, ["p", "q"]);
        check(close(sd2.mean[0, 1], mean(y)) && sd2.names == ["p", "q"], "summary on a Dataset");
    }
}

void stage3(double[] x, double[] y) {
    enum n = 50;
    auto f = frequencies(x);

    check(!f.discrete && f.nobs == n && f.count.nrow == f.midpoint.nrow && f.count.nrow % 2 == 1, "frequencies: an odd number of equal-width bins chosen from n");
    double tot = 0;
    foreach (i; 0 .. f.count.nrow) tot += f.count[i, 0];
    check(tot == n && close(f.cumulative[f.count.nrow - 1, 0], 1.0) && close(f.fraction[0, 0], f.count[0, 0] / n), "frequencies: the counts add to n, the cumulative fraction to 1");
    bool contiguous = true;
    foreach (i; 0 .. f.count.nrow - 1) contiguous = contiguous && close(f.upper[i, 0], f.lower[i + 1, 0]) && close(f.midpoint[i, 0], (f.lower[i, 0] + f.upper[i, 0]) / 2);
    check(contiguous && f.lower[0, 0] < minimum(x) && f.upper[f.count.nrow - 1, 0] > maximum(x), "frequencies: contiguous bins covering the data, midpoints in the middle");
    // each count against a direct count
    bool counts = true;
    foreach (i; 0 .. f.count.nrow) {
        double c = 0;
        foreach (v; x) {
            const first = i == 0, last = i == f.count.nrow - 1;
            if ((first || v >= f.lower[i, 0]) && (last || v < f.upper[i, 0])) c++;
        }
        counts = counts && c == f.count[i, 0];
    }
    check(counts, "frequencies: every count agrees with a direct count");
    check(close(f.mean, mean(x)) && close(f.sd, stddev(x)), "frequencies: mean and standard deviation");
    auto f9 = frequencies(x, 9);
    check(f9.count.nrow == 9, "frequencies: a given number of bins");
    auto fw = frequencies(x, 0, 0.0, 2.0);
    check(close(fw.lower[0, 0], 0.0) && close(fw.upper[0, 0] - fw.lower[0, 0], 2.0) && fw.count.nrow == cast(size_t) ((maximum(x) - 0.0) / 2.0 + 1) , "frequencies: a given minimum and width");
    auto fn = frequencies(x, 0, double.nan, double.nan, true);
    check(!isNaN(fn.normalTest) && fn.normalPvalue > 0 && fn.normalPvalue <= 1 && isNaN(f.normalTest), "frequencies: the normality test only when asked for");
    // a roughly normal sample is not rejected, a strongly skewed one is
    double[] skewed = new double[200];
    foreach (i; 0 .. 200) skewed[i] = (urand() + 0.5) ^^ 6;
    check(frequencies(skewed, 0, double.nan, double.nan, true).normalPvalue < 0.01, "frequencies: a skewed sample is rejected as normal");

    // discrete data
    double[] dd = new double[100];
    foreach (i; 0 .. 100) dd[i] = cast(int) ((urand() + 0.5) * 5);
    auto fd = valueCounts(dd);
    double tot2 = 0;
    foreach (i; 0 .. fd.count.nrow) tot2 += fd.count[i, 0];
    check(fd.discrete && fd.count.nrow <= 5 && tot2 == 100 && fd.midpoint[0, 0] == minimum(dd), "valueCounts: one row per value, counts add to n");
    bool vals = true;
    foreach (i; 0 .. fd.count.nrow) {
        double c = 0;
        foreach (v; dd) if (v == fd.midpoint[i, 0]) c++;
        vals = vals && c == fd.count[i, 0] && fd.lower[i, 0] == fd.upper[i, 0];
    }
    check(vals, "valueCounts: counts agree with direct counts");
    check(frequencies(dd).discrete, "frequencies: a series with few distinct values is tabulated by value");

    // missing values, errors
    auto xm = x.dup;
    xm[0] = double.nan;
    check(frequencies(xm).nobs == n - 1, "frequencies: missing values are left out");
    check(throws({ frequencies(new double[n]); }) && throws({ frequencies(x[0 .. 5]); }) && throws({ frequencies(x, 1); }) && throws({ frequencies(x, 5, 0.0, 1.0); })
          && throws({ frequencies(x, 0, 1.0); }) && throws({ frequencies(x, 0, 0.0, -1.0); }), "frequencies: constant, short, or inconsistent arguments throw");
}

void stage4(double[] x, double[] y, double[] z) {
    import gretl.pvalues: chisq_cdf_comp, snedecor_cdf_comp;
    enum n = 50;

    // normality: a sum of uniforms is close to normal, a sixth power of a uniform is not
    enum m = 400;
    double[] g = new double[m], sk = new double[m];
    foreach (i; 0 .. m) {
        g[i] = urand() + urand() + urand() + urand() + urand() + urand();
        sk[i] = (urand() + 0.5) ^^ 6;
    }
    auto ng = normalityTests(g), ns = normalityTests(sk);
    check(ng.doornikHansen.pvalue > 0.01 && ng.shapiroWilk.pvalue > 0.01 && ng.jarqueBera.pvalue > 0.01 && ng.lilliefors.pvalue > 0.01, "normality tests: a near-normal sample is not rejected");
    check(ns.doornikHansen.pvalue < 1e-6 && ns.shapiroWilk.pvalue < 1e-6 && ns.jarqueBera.pvalue < 1e-6 && ns.lilliefors.pvalue < 0.01, "normality tests: a skewed sample is rejected by all four");
    const jb = m / 6.0 * (skewness(sk) ^^ 2 + kurtosis(sk) ^^ 2 / 4);
    check(close(ns.jarqueBera.statistic, jb, 1e-8) && close(ns.jarqueBera.pvalue, chisq_cdf_comp(2, jb), 1e-8), "normality tests: Jarque-Bera from the skewness and kurtosis");
    auto gm = g.dup;
    gm[3] = double.nan;
    check(normalityTests(gm).shapiroWilk.statistic > 0.9 && throws({ normalityTests(g[0 .. 5]); }) && throws({ normalityTests(new double[20]); }), "normality tests: missing values are left out; too few or constant samples throw");

    // multivariate normality
    auto E = cbind(g[0 .. 200], g[200 .. 400]);
    auto Es = cbind(sk[0 .. 200], sk[200 .. 400]);
    auto mn = multivariateNormality(E), ms = multivariateNormality(Es);
    check(mn.pvalue > 0.01 && ms.pvalue < 1e-4 && mn.statistic > 0, "multivariate normality: a normal sample is not rejected, a skewed one is");
    check(throws({ multivariateNormality(cbind(g[0 .. 20])); }) && throws({ multivariateNormality(E, eye(3)); }), "multivariate normality: bad arguments throw");

    // means tests
    auto x1 = x[0 .. 30], y1 = y[0 .. 40];
    double mx = mean(x1), my = mean(y1), vx = variance(x1), vy = variance(y1);
    auto eq = meansTest(x1, y1);
    const sp = (29 * vx + 39 * vy) / 68, tq = (mx - my) / sqrt(sp / 30 + sp / 40);
    check(eq.n1 == 30 && eq.n2 == 40 && eq.df == 68 && close(eq.statistic, tq, 1e-9) && close(eq.difference, mx - my) && close(eq.se, sqrt(sp / 30 + sp / 40), 1e-9), "meansTest: the pooled t statistic");
    import gretl.pvalues: student_pvalue_2;
    check(close(eq.pvalue, student_pvalue_2(68, tq), 1e-9), "meansTest: the two-sided p-value");
    auto un = meansTest(x1, y1, MeansTest.unequalVariances);
    const tu = (mx - my) / sqrt(vx / 30 + vy / 40);
    const dfu = ((vx / 30 + vy / 40) ^^ 2) / ((vx / 30) ^^ 2 / 29 + (vy / 40) ^^ 2 / 39);
    check(close(un.statistic, tu, 1e-9) && un.df == cast(double) cast(int) dfu, "meansTest: Welch's statistic and the Satterthwaite degrees of freedom (rounded down)");
    double[] diffs;
    foreach (i; 0 .. n) diffs ~= x[i] - y[i];
    auto pd = meansTest(x, y, MeansTest.paired);
    check(pd.n1 == n && pd.df == n - 1 && close(pd.statistic, mean(diffs) / (stddev(diffs) / sqrt(cast(double) n)), 1e-9) && close(pd.difference, mean(diffs)), "meansTest: the paired t statistic");
    auto xn = x1.dup;
    xn[4] = double.nan;
    check(meansTest(xn, y1).n1 == 29, "meansTest: missing values are left out");
    check(throws({ meansTest(x, y[0 .. 10], MeansTest.paired); }) && throws({ meansTest([1.0], y1); }), "meansTest: bad samples throw");

    // variance test
    auto vt = variancesTest(x1, y1);
    const F = vx / vy;
    check(vt.test.n1 == 30 && vt.test.n2 == 40 && vt.test.df == 29 && vt.test.df2 == 39 && close(vt.test.statistic, F, 1e-9), "variancesTest: the F statistic and its degrees of freedom");
    check(close(vt.pvalueUpper, snedecor_cdf_comp(29, 39, F), 1e-9) && close(vt.test.pvalue, 2 * (vt.pvalueUpper < 0.5 ? vt.pvalueUpper : 1 - vt.pvalueUpper), 1e-6), "variancesTest: the two-sided p-value, and the library's one-sided one");

    // distances
    auto X = cbind(x[0 .. 6], y[0 .. 6], z[0 .. 6]);
    auto Y = cbind(x[10 .. 14], y[10 .. 14], z[10 .. 14]);
    auto de = rowDistances(X, Y), dm = rowDistances(X, Y, DistanceMetric.manhattan), dc = rowDistances(X, Y, DistanceMetric.chebyshev);
    bool ok = de.nrow == 6 && de.ncol == 4;
    foreach (i; 0 .. 6) foreach (j; 0 .. 4) {
        double e2 = 0, e1 = 0, ei = 0;
        foreach (c; 0 .. 3) {
            const dd = abs(X[i, c] - Y[j, c]);
            e2 += dd * dd; e1 += dd; if (dd > ei) ei = dd;
        }
        ok = ok && close(de[i, j], sqrt(e2), 1e-9) && close(dm[i, j], e1, 1e-9) && close(dc[i, j], ei, 1e-9);
    }
    check(ok, "rowDistances: Euclidean, Manhattan and Chebyshev distances between rows");
    auto dmat = rowDistanceMatrix(X);
    ok = dmat.nrow == 6 && dmat.ncol == 6;
    foreach (i; 0 .. 6) foreach (j; 0 .. 6) {
        double e2 = 0;
        foreach (c; 0 .. 3) e2 += (X[i, c] - X[j, c]) ^^ 2;
        ok = ok && close(dmat[i, j], sqrt(e2), 1e-9);
    }
    check(ok, "rowDistanceMatrix: symmetric, zero diagonal, the pairwise Euclidean distances");
    auto dcos = rowDistances(X, X, DistanceMetric.cosine);
    check(close(dcos[2, 2], 0, 1e-9), "rowDistances: cosine distance of a row to itself is 0");
    check(throws({ rowDistances(X, cbind(x[0 .. 5], y[0 .. 5])); }) && throws({ rowDistanceMatrix(cbind(x[0 .. 1])); }), "rowDistances: mismatched columns or a single row throw");

    // Mahalanobis
    auto Xm = cbind(x, y);
    auto md = mahalanobis(Xm);
    auto S = covariance(Xm);
    const mux = mean(x), muy = mean(y);
    auto Si = inv(S);
    bool mok = md.nrow == n;
    foreach (i; 0 .. n) {
        const d0 = x[i] - mux, d1 = y[i] - muy;
        const q = d0 * (Si[0, 0] * d0 + Si[0, 1] * d1) + d1 * (Si[1, 0] * d0 + Si[1, 1] * d1);
        mok = mok && close(md[i, 0], sqrt(q), 1e-8);
    }
    check(mok, "mahalanobis: the distance from the centroid with the sample covariance matrix");
    auto ym2 = y.dup;
    ym2[7] = double.nan;
    auto md2 = mahalanobis(cbind(x, ym2));
    check(isNaN(md2[7, 0]) && !isNaN(md2[8, 0]), "mahalanobis: a row with a missing value gets NaN");

    // crosstab
    double[] a = new double[100], b = new double[100];
    foreach (i; 0 .. 100) { a[i] = cast(int) ((urand() + 0.5) * 3); b[i] = cast(int) ((urand() + 0.5) * 2) + (a[i] == 2 ? 1 : 0); }
    auto ct = crosstab(a, b);
    double tot = 0;
    foreach (i; 0 .. ct.counts.nrow) foreach (j; 0 .. ct.counts.ncol) tot += ct.counts[i, j];
    bool cok = tot == 100 && ct.nobs == 100 && ct.df == (ct.counts.nrow - 1) * (ct.counts.ncol - 1);
    foreach (i; 0 .. ct.counts.nrow) foreach (j; 0 .. ct.counts.ncol) {
        double c = 0;
        foreach (k; 0 .. 100) if (a[k] == ct.rowValues[i, 0] && b[k] == ct.colValues[j, 0]) c++;
        cok = cok && c == ct.counts[i, j];
    }
    check(cok, "crosstab: every cell agrees with a direct count");
    double rt = 0, chi = 0;
    foreach (i; 0 .. ct.counts.nrow) {
        double r = 0;
        foreach (j; 0 .. ct.counts.ncol) r += ct.counts[i, j];
        rt += abs(r - ct.rowTotals[i, 0]);
        foreach (j; 0 .. ct.counts.ncol) {
            double c = 0;
            foreach (ii; 0 .. ct.counts.nrow) c += ct.counts[ii, j];
            const e = r * c / 100;
            chi += (ct.counts[i, j] - e) ^^ 2 / e;
        }
    }
    check(rt == 0 && close(ct.pearson, chi, 1e-9) && close(ct.pvalue, chisq_cdf_comp(ct.df, chi), 1e-9) && ct.pvalue < 0.05, "crosstab: the marginal totals and Pearson's chi-square");
    auto an = a.dup;
    an[5] = double.nan;
    check(crosstab(an, b).nobs == 99 && throws({ crosstab(a, b[0 .. 10]); }), "crosstab: missing values are left out; different lengths throw");
}

void stage5(double[] x, double[] y, double[] z) {
    enum n = 50;
    auto d = Dataset(4, n);
    import core.stdc.string: strcpy;
    strcpy(d.ptr.varname[1], "a"); strcpy(d.ptr.varname[2], "b"); strcpy(d.ptr.varname[3], "c");
    foreach (i; 0 .. n) { d.ptr.Z[1][i] = x[i]; d.ptr.Z[2][i] = y[i]; d.ptr.Z[3][i] = z[i]; }

    check(close(mean(d["a"]), mean(x)) && close(quantile(d["b"], 0.75), quantile(y, 0.75)) && close(skewness(d["c"]), skewness(z)), "a series of a Dataset goes straight into the functions");
    check(frequencies(d["a"]).nobs == n && normalityTests(d["a"]).jarqueBera.pvalue >= 0, "frequencies and normalityTests on a series of a Dataset");
    auto cv = covariance(d, ["a", "c"]), cr = correlation(d, ["a", "b", "c"]);
    check(cv.nrow == 2 && close(cv[0, 1], covariance(x, z)) && cr.nrow == 3 && close(cr[0, 1], correlation(x, y)), "covariance and correlation by names");
    check(mahalanobis(d, ["a", "b"]).nrow == n && close(mahalanobis(d, ["a", "b"])[3, 0], mahalanobis(cbind(x, y))[3, 0]), "mahalanobis by names");
    check(close(meansTest(d, "a", "b").statistic, meansTest(x, y).statistic) && close(variancesTest(d, "a", "b").test.statistic, variancesTest(x, y).test.statistic), "two-sample tests by names");
    double[] da = new double[n], db = new double[n];
    foreach (i; 0 .. n) { da[i] = i % 3; db[i] = i % 2; }
    auto d2 = Dataset(3, n);
    strcpy(d2.ptr.varname[1], "p"); strcpy(d2.ptr.varname[2], "q");
    foreach (i; 0 .. n) { d2.ptr.Z[1][i] = da[i]; d2.ptr.Z[2][i] = db[i]; }
    check(crosstab(d2, "p", "q").counts.nrow == 3 && throws({ covariance(d, ["a", "nothere"]); }) && throws({ correlation(d, null); }), "crosstab by names; an unknown series or no names throw");
}

void stage6() {
    import gretl.pvalues: binomial_cdf, normal_cdf;
    enum n = 40;
    double[] x = new double[n], y = new double[n];

    foreach (i; 0 .. n) {
        y[i] = cast(int) (urand() * 20);                        // integers: ties and zero differences
        x[i] = y[i] + cast(int) (urand() * 8 + 1.5);            // x tends to be larger
    }
    x[3] = y[3]; x[11] = y[11]; x[20] = y[20];                  // zero differences

    // sign test
    auto st = signTest(x, y);
    size_t nn = 0, w = 0;
    foreach (i; 0 .. n) if (x[i] != y[i]) { nn++; if (x[i] > y[i]) w++; }
    check(st.n == nn && st.w == w && st.report.length > 0, "signTest: n and W count the pairs that differ");
    check(close(st.pvalueLower, binomial_cdf(0.5, cast(int) nn, cast(int) w), 1e-12) && close(st.pvalueUpper, 1 - binomial_cdf(0.5, cast(int) nn, cast(int) w - 1), 1e-9), "signTest: the binomial tail probabilities");
    check(st.pvalue < 0.001 && close(st.pvalue, 2 * (st.pvalueUpper < st.pvalueLower ? st.pvalueUpper : st.pvalueLower)), "signTest: x is larger, the two-sided p-value is small");
    auto sym = signTest(y, x);
    check(sym.w == nn - w && close(sym.pvalue, st.pvalue, 1e-9), "signTest: swapping the series swaps W and keeps the two-sided p-value");
    auto xm = x.dup, ym = y.dup;
    xm[0] = double.nan;
    check(signTest(xm, ym).n == (x[0] != y[0] ? nn - 1 : nn), "signTest: a missing value drops its pair");
    check(throws({ signTest(x, y[0 .. 10]); }) && throws({ signTest(y, y); }), "signTest: different lengths, or no differences at all, throw");

    // signed-rank test: ranks of |d| with zeros first, the sums of ranks, Cureton's variance
    auto sr = signedRankTest(x, y);
    double[] ad; bool[] pos;
    foreach (i; 0 .. n) if (x[i] != y[i]) { ad ~= abs(x[i] - y[i]); pos ~= x[i] > y[i]; }
    const Z = n - ad.length;
    auto order = new size_t[ad.length];
    foreach (i, ref v; order) v = i;
    import std.algorithm.sorting: stdSort = sort;
    stdSort!((a, b) => ad[a] < ad[b])(order);
    double wp = 0, wm = 0, tsum = 0;
    size_t ties = 0;
    for (size_t i = 0; i < order.length;) {
        size_t j = i;
        while (j + 1 < order.length && ad[order[j + 1]] == ad[order[i]]) j++;
        const rk = (Z + i + 1 + Z + j + 1) / 2.0;
        foreach (k; i .. j + 1) { if (pos[order[k]]) wp += rk; else wm += rk; }
        if (j > i) { tsum += (cast(double) (j - i)) ^^ 3 - (j - i); ties++; }
        i = j + 1;
    }
    check(sr.pairs == n && sr.n == ad.length && sr.zeroDifferences == Z && sr.nonZeroTies == ties, "signedRankTest: counts of pairs, differences, zero differences and ties");
    check(close(sr.wPlus, wp) && close(sr.wMinus, wm) && close(wp + wm, (cast(double) n * (n + 1) - Z * (Z + 1)) / 2.0, 1e-9), "signedRankTest: W+ and W- add up to the sum of the ranks");
    const xe = (cast(double) n * (n + 1) - Z * (Z + 1)) / 4.0;
    const s2 = (cast(double) n * (n + 1) * (2 * n + 1) - Z * (Z + 1) * (2 * Z + 1) - tsum / 2) / 24.0;
    double num = wp - xe;
    num = num > 0.25 ? num - 0.5 : num + 0.5;                      // the library's continuity correction
    const zz = num / sqrt(s2);
    check(close(sr.expected, xe) && close(sr.variance, s2) && close(sr.statistic, zz, 1e-9), "signedRankTest: the expected value, Cureton's variance and the z statistic");
    check(close(sr.pvalueOneSided, zz > 0 ? 1 - normal_cdf(zz) : normal_cdf(zz), 1e-9) && close(sr.pvalue, 2 * sr.pvalueOneSided, 1e-9) && sr.pvalue < 0.001, "signedRankTest: the tail probability and the two-sided p-value");
    auto xn2 = x.dup, yn2 = y.dup;
    xn2[5] = double.nan;
    yn2[17] = double.nan;
    auto srm = signedRankTest(xn2, yn2);
    check(srm.pairs == n - 2 && srm.n <= sr.n && !isNaN(srm.statistic), "signedRankTest: pairs with a missing value in either series are dropped, the others stay paired");
    auto small = signedRankTest(x[0 .. 8], y[0 .. 8]);
    check(isNaN(small.statistic) && isNaN(small.pvalue) && small.report.length > 0 && small.wPlus + small.wMinus > 0, "signedRankTest: a sample of 8 or fewer differences has no normal approximation");
    check(throws({ signedRankTest(x, y[0 .. 3]); }) && throws({ signedRankTest([double.nan, double.nan, 1.0], [1.0, double.nan, double.nan]); }), "signedRankTest: different lengths or no valid pairs throw");
    auto zero = signedRankTest(y, y);
    check(zero.n == 0 && zero.zeroDifferences == n && isNaN(zero.statistic) && zero.report.length > 0, "signedRankTest: identical series give n = 0 and no statistic, not an error");

    // rank-sum test
    double[] a = new double[25], b = new double[18];
    foreach (i; 0 .. 25) a[i] = cast(int) (urand() * 15 + 8);
    foreach (i; 0 .. 18) b[i] = cast(int) (urand() * 15 + 2);
    auto rs = rankSumTest(a, b);
    auto all = a ~ b;
    auto ord = new size_t[all.length];
    foreach (i, ref v; ord) v = i;
    stdSort!((p, q) => all[p] < all[q])(ord);
    double ra = 0, rb = 0;
    for (size_t i = 0; i < ord.length;) {
        size_t j = i;
        while (j + 1 < ord.length && all[ord[j + 1]] == all[ord[i]]) j++;
        const rk = (i + 1 + j + 1) / 2.0;
        foreach (k; i .. j + 1) { if (ord[k] < 25) ra += rk; else rb += rk; }
        i = j + 1;
    }
    check(rs.nx == 25 && rs.ny == 18 && close(rs.rankSumX, ra) && close(rs.rankSumY, rb) && close(ra + rb, 43.0 * 44 / 2), "rankSumTest: the rank sums of the pooled ranking");
    // z is for the smaller sample (b)
    const mb = 18.0 * 44 / 2, sb = sqrt(25.0 * 18 * 44 / 12);
    const zb = (rb - mb) / sb;
    check(close(rs.statistic, zb, 1e-9) && close(rs.pvalue, 2 * (zb > 0 ? 1 - normal_cdf(zb) : normal_cdf(zb)), 1e-9) && rs.pvalue < 0.01, "rankSumTest: z for the smaller sample and the two-sided p-value");
    auto rsw = rankSumTest(b, a);
    check(close(rsw.statistic, zb, 1e-9) && close(rsw.pvalue, rs.pvalue, 1e-9), "rankSumTest: swapping the samples gives the same z (the smaller sample's)");
    auto rss = rankSumTest(a[0 .. 6], b[0 .. 5]);
    check(isNaN(rss.statistic) && isNaN(rss.pvalue) && rss.nx == 6 && rss.ny == 5 && close(rss.rankSumX + rss.rankSumY, 66.0), "rankSumTest: small samples have no normal approximation, but the rank sums are there");
    check(rss.critical.available && close(rss.rankSumSmaller, rss.rankSumY) && rss.critical.lower1 == 17 && rss.critical.lower5 == 20 && rss.critical.lower10 == 22
          && rss.critical.upper10 == 38 && rss.critical.upper5 == 40 && rss.critical.upper1 == 43, "rankSumTest: critical values of the rank sum of the smaller sample (5 and 6 observations)");
    auto r44 = rankSumTest(a[0 .. 4], b[0 .. 4]);
    check(r44.critical.available && isNaN(r44.critical.lower1) && r44.critical.lower5 == 11 && r44.critical.lower10 == 13 && r44.critical.upper10 == 23 && r44.critical.upper5 == 25 && isNaN(r44.critical.upper1)
          && close(r44.rankSumSmaller, r44.rankSumX), "rankSumTest: 4 and 4 observations: only 5% and 10%, with the right upper values (the library prints them wrongly)");
    check(!rankSumTest(a[0 .. 3], b[0 .. 8]).critical.available && !rankSumTest(a[0 .. 6], b[0 .. 14]).critical.available && !rankSumTest(a, b).critical.available, "rankSumTest: no critical values outside the table");
    auto am = a.dup;
    am[2] = double.nan;
    check(rankSumTest(am, b).nx == 24 && throws({ rankSumTest(new double[0], b); }) && throws({ rankSumTest([double.nan], b); }), "rankSumTest: missing values are left out; an empty sample throws");

    // Dataset forms
    auto d = Dataset(3, n);
    import core.stdc.string: strcpy;
    strcpy(d.ptr.varname[1], "u"); strcpy(d.ptr.varname[2], "v");
    foreach (i; 0 .. n) { d.ptr.Z[1][i] = x[i]; d.ptr.Z[2][i] = y[i]; }
    check(signTest(d, "u", "v").w == w && close(signedRankTest(d, "u", "v").statistic, sr.statistic) && close(rankSumTest(d, "u", "v").statistic, rankSumTest(x, y).statistic), "the Dataset forms of the sign, signed-rank and rank-sum tests");
}
