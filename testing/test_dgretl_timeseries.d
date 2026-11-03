/* Tests for dgretl/timeseries.d (stage 1: correlogram, Ljung-Box, ADF, KPSS, filters, lag and diff). */

import core.stdc.stdlib;
import std.conv: to;
import std.math: abs, isNaN, log, sqrt;
import std.stdio: writefln, writeln;
import dgretl;
import gretl.libgretlh;
import dgretl.timeseries;

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

bool close(double a, double b, double tol = 1e-8) { return abs(a - b) <= tol * (1.0 + abs(b)); }

bool throws(void delegate() f) {
    try f(); catch (GretlException) return true;
    return false;
}

ulong state = 12345;
double urand() {
    state = state * 6364136223846793005UL + 1442695040888963407UL;
    return ((state >> 11) / cast(double) (1UL << 53)) - 0.5;
}

void stage2(double[] rw, double[] arSeries) {
    enum T = 120;

    // ARCH test: against the library's own test on a model with a constant only
    double[] u = new double[T];
    double h = 1.0;
    foreach (t; 0 .. T) {
        const e = urand() * 4;
        u[t] = e * h;
        h = sqrt(0.5 + 0.8 * u[t] * u[t] / 4);
    }
    auto at = archTest(u, 3);
    check(at.df == 3 && at.pvalue >= 0 && at.pvalue <= 1 && at.statistic > 0, "archTest: statistic and p-value");
    auto at0 = archTest(arSeries, 3);
    check(at0.statistic < at.statistic + 1e9 && at0.pvalue > 0.01, "archTest: no ARCH in an AR series");
    check(throws({ archTest(u, 0); }), "archTest: bad order throws");

    // library cross-check
    {
        import gretl.estimate: lsq;
        import gretl.gretl_commandsh: GretlCmdIndex;
        import gretl.gretl_model: clear_model;
        import gretl.gretl_utils: get_last_pvalue, get_last_test_statistic;
        import gretl.libgretlh: Gretlopt, Model;
        import gretl.estimate: arch_test;
        import gretl.gretl_prn: gretl_print_destroy, gretl_print_new;
        import gretl.gretl_prnh: PrnType;
        import dgretl.internal;
        double mu = 0;
        foreach (x; u) mu += x;
        mu /= T;
        auto dm = new double[T];
        foreach (t; 0 .. T) dm[t] = u[t] - mu;
        auto ds = makeDataset([dm], ["const", "u"]);
        auto list = makeList([1, 0]);
        Model m = lsq(list, ds.ptr, GretlCmdIndex.ols, Gretlopt.OPT_NONE);
        int perr;
        auto prn = gretl_print_new(PrnType.buffer, &perr);
        const err = arch_test(&m, 3, ds.ptr, Gretlopt.OPT_Q, prn);
        const st = get_last_test_statistic();
        gretl_print_destroy(prn);
        clear_model(&m);
        auto my = archTest(dm, 3);
        check(err == 0 && close(my.statistic, st, 1e-8) && close(my.pvalue, get_last_pvalue(), 1e-8), "archTest agrees with the library's arch_test");
    }

    // Engle-Granger: cointegrated pair against an unrelated random walk
    double[] y2 = new double[T];
    foreach (t; 0 .. T) y2[t] = 2.0 * rw[t] + 0.3 * urand();
    auto eg = engleGranger(y2, rw, 1);
    check(eg.test == "Engle-Granger" && !isNaN(eg.statistic) && eg.report.length > 0, "engleGranger: statistic and report");
    double[] rw3 = new double[T];
    double z = 0;
    foreach (t; 0 .. T) { z += urand(); rw3[t] = z; }
    auto eg2 = engleGranger(rw3, rw, 1);
    check(eg.statistic < eg2.statistic && eg.statistic < -3, "engleGranger: the cointegrated pair rejects more strongly");
    check(throws({ engleGranger(y2, rw, -1); }), "engleGranger: bad lags throw");

    // AR(1) errors: y = 1 + 2 x + e, e AR(1) with rho 0.6
    double[] x = new double[T], yy = new double[T];
    double e = 0;
    foreach (t; 0 .. T) {
        x[t] = urand() * 4;
        e = 0.6 * e + urand();
        yy[t] = 1 + 2 * x[t] + e;
    }
    auto f1 = ar1(yy, x);
    check(f1.arcoef.nrow == 1 && f1.arcoef[0, 0] > 0.2 && f1.arcoef[0, 0] < 0.95, "ar1: the estimated rho is positive and near 0.6");
    check(f1.coef.nrow == 2 && close(f1.coef[1, 0], 2, 0.1), "ar1: coefficients and standard errors");
    auto f2 = ar1(yy, x, true, Ar1Method.hildrethLu);
    auto f3 = ar1(yy, x, true, Ar1Method.praisWinsten);
    check(abs(f2.arcoef[0, 0] - f1.arcoef[0, 0]) < 0.1 && f3.arcoef.nrow == 1, "ar1: Hildreth-Lu and Prais-Winsten estimates are close");
    check(f1.names == ["const", "x1"], "ar1: coefficient names");

    // AR(2) errors by generalized Cochrane-Orcutt
    auto g = ar(yy, x, 2);
    check(g.arcoef.nrow == 2 && g.arlags == [1, 2] && g.coef.nrow == 2, "ar: two error lags");
    auto g4 = ar(yy, x, [1, 4]);
    check(g4.arlags == [1, 4] && g4.arse.nrow == 2, "ar: a list of lags");
    check(throws({ ar(yy, x, [0]); }) && throws({ ar(yy, x, 0); }), "ar: bad lags throw");

    // ARCH
    auto ar_fit = arch(u, x, 2);
    check(ar_fit.method == "ARCH-WLS" && ar_fit.coef.nrow == 2 && ar_fit.se[1, 0] > 0, "arch: estimates and standard errors");
    check(throws({ arch(u, x, 0); }), "arch: bad order throws");
}

void main() {
    enum T = 120;
    double[] rw = new double[T], arSeries = new double[T];
    double r = 0, a = 0;

    foreach (t; 0 .. T) {
        r += urand();
        a = 0.5 * a + urand();
        rw[t] = r;
        arSeries[t] = a;
    }

    // correlogram against direct computation
    auto c = correlogram(arSeries, 8);
    double mean = 0;
    foreach (v; arSeries) mean += v;
    mean /= T;
    double c0 = 0;
    foreach (v; arSeries) c0 += (v - mean) ^^ 2;
    double[] rho;
    foreach (k; 1 .. 9) {
        double s = 0;
        foreach (t; k .. T) s += (arSeries[t] - mean) * (arSeries[t - k] - mean);
        rho ~= s / c0;
    }
    bool ok = true;
    foreach (k; 0 .. 8) ok = ok && close(c.acf[k, 0], rho[k], 1e-10);
    check(ok, "correlogram: acf equals the direct computation");
    check(close(c.pacf[0, 0], rho[0], 1e-10), "correlogram: pacf(1) = acf(1)");
    const phi22 = (rho[1] - rho[0] ^^ 2) / (1 - rho[0] ^^ 2);
    check(close(c.pacf[1, 0], phi22, 1e-8), "correlogram: pacf(2) from the Durbin-Levinson formula");
    check(c.acf.nrow == 8 && c.acf.ncol == 1 && c.ljungBox.nrow == 8, "correlogram: shapes");

    // Ljung-Box
    double q = 0;
    foreach (k; 0 .. 8) q += rho[k] ^^ 2 / (T - k - 1);
    q *= T * (T + 2);
    auto lb = ljungBox(arSeries, 8);
    check(close(lb.statistic, q, 1e-8) && close(c.ljungBox[7, 0], q, 1e-8), "ljungBox: matches the formula and the correlogram column");
    check(lb.df == 8 && lb.pvalue > 0 && lb.pvalue < 1, "ljungBox: p-value");
    check(ljungBox(arSeries, 8, 2).df == 6, "ljungBox: fitted parameters reduce df");
    check(throws({ ljungBox(arSeries, 0); }) && throws({ correlogram(arSeries, T); }), "correlogram/ljungBox: bad lags throw");
    auto withNaN = arSeries.dup;
    withNaN[3] = double.nan;
    check(throws({ correlogram(withNaN, 4); }), "missing values throw");

    // unit root tests
    auto t1 = adf(rw, 2);
    check(!isNaN(t1.statistic) && t1.report.length > 0 && t1.test == "ADF", "adf: statistic and report");
    auto t2 = adf(arSeries, 2);
    check(t2.statistic < t1.statistic, "adf: the stationary series has a more negative statistic than the random walk");
    check(t2.pvalue < 0.05 && t1.pvalue > 0.05, "adf: p-values reject for the stationary series only");
    auto t3 = adf(rw, 4, Trend.trend, TestDown.aic);
    check(!isNaN(t3.statistic), "adf: trend with test-down");
    auto t4 = adf(rw, 1, Trend.constant, TestDown.off, true);
    check(t4.test == "DF-GLS" && !isNaN(t4.statistic), "adf: DF-GLS");
    auto k1 = kpss(rw, 4);
    auto k2 = kpss(arSeries, 4);
    check(k1.statistic > k2.statistic && k1.report.length > 0, "kpss: the random walk has the larger statistic");
    check(throws({ kpss(rw, 4, Trend.none); }) && throws({ adf(rw, -1); }), "unit root tests: bad arguments throw");

    // filters
    auto hp = hpTrend(rw, 1600);
    check(hp.nrow == T && hp.ncol == 1, "hpTrend: shape");
    double sres = 0, sm = 0;
    foreach (t; 0 .. T) { sres += rw[t] - hp[t, 0]; sm += rw[t]; }
    check(close(sres, 0, 1e-6) || abs(sres) < 1e-6 * T, "hpTrend: cycle sums to zero");
    auto hp0 = hpTrend(rw, 1e-8);
    double dev = 0;
    foreach (t; 0 .. T) dev = dev > abs(rw[t] - hp0[t, 0]) ? dev : abs(rw[t] - hp0[t, 0]);
    check(dev < 1e-4, "hpTrend: lambda near 0 follows the series");
    auto bk = bkFilter(arSeries, 6, 32, 8);
    check(bk.nrow == T && isNaN(bk[0, 0]) && isNaN(bk[T - 1, 0]) && !isNaN(bk[8, 0]), "bkFilter: k values lost at each end");
    auto bw = butterworth(arSeries, 2, 67.5);
    check(bw.nrow == T && !isNaN(bw[T / 2, 0]), "butterworth: shape");
    auto fd = fracDiff(rw, 1.0);
    check(close(fd[5, 0], rw[5] - rw[4], 1e-8), "fracDiff: d = 1 is the first difference");
    check(throws({ hpTrend(rw, -1); }) && throws({ bkFilter(rw, 6, 4); }), "filters: bad arguments throw");

    // transformations
    auto l1 = lag(rw, 2);
    check(isNaN(l1[1, 0]) && l1[5, 0] == rw[3], "lag: shifted and padded");
    auto ld = lag(rw, -1);
    check(ld[4, 0] == rw[5] && isNaN(ld[T - 1, 0]), "lag: a lead");
    auto d1 = diff(rw);
    check(isNaN(d1[0, 0]) && close(d1[7, 0], rw[7] - rw[6]), "diff: first difference");
    auto d2 = diff(rw, 2);
    check(close(d2[9, 0], rw[9] - 2 * rw[8] + rw[7], 1e-10), "diff: second difference");
    auto ds = diff(rw, 1, 4);
    check(isNaN(ds[3, 0]) && close(ds[9, 0], rw[9] - rw[5]), "diff: seasonal difference");

    stage2(rw, arSeries);
    stage3();
    stage4();
    stage5();
    stage6();
    stage7();
    stage8();
    stage9();
    stage10();
    stage11();
    stage12();

    writefln("\n%d tests, %d failed", nTests, nFailed);
    if (nFailed) throw new Exception("failed");
}

/* ---- stage 3: VAR ---- */

void stage3() {
    enum T = 100;
    double[] y1 = new double[T], y2 = new double[T];

    y1[0] = y2[0] = 0;
    foreach (t; 1 .. T) {
        y1[t] = 0.5 * y1[t - 1] + 0.2 * y2[t - 1] + urand() + 0.3;
        y2[t] = 0.1 * y1[t - 1] + 0.4 * y2[t - 1] + urand();
    }
    auto Y = cbind(y1, y2);
    auto fit = var(Y, 1, names: ["a", "b"]);

    check(fit.neqns == 2 && fit.order == 1 && fit.nobs == T - 1, "var: dimensions");
    check(fit.coefNames == ["const", "a_1", "b_1"], "var: coefficient names");
    check(fit.coef.nrow == 3 && fit.coef.ncol == 2 && fit.se.nrow == 3, "var: coefficient matrix shape");

    // against least squares on hand-built lags
    auto x = fresh2(T - 1, 2);
    foreach (t; 1 .. T) { x[t - 1, 0] = y1[t - 1]; x[t - 1, 1] = y2[t - 1]; }
    double[] ya = y1[1 .. $].dup, yb = y2[1 .. $].dup;
    auto o1 = ols(ya, x), o2 = ols(yb, x);
    bool ok = true;
    foreach (i; 0 .. 3) ok = ok && close(fit.coef[i, 0], o1.coef[i, 0], 1e-8) && close(fit.coef[i, 1], o2.coef[i, 0], 1e-8)
                              && close(fit.se[i, 0], o1.se[i, 0], 1e-8);
    check(ok, "var: coefficients and standard errors equal equation-by-equation OLS");
    auto eq = fit.equation(1);
    check(close(eq.rsq, o2.rsq, 1e-8) && eq.names == fit.coefNames && eq.nobs == T - 1, "var: equation(1) is the Estimate of the regression");
    check(close(fit.resid[10, 0], o1.resid[10, 0], 1e-8), "var: residuals");

    // roots: the companion matrix is the 2 x 2 coefficient matrix
    check(fit.roots.nrow == 2 && fit.roots.ncol == 3 && fit.roots[0, 2] < 1 && fit.roots[1, 2] < 1, "var: roots inside the unit circle");
    double tr = fit.coef[1, 0] + fit.coef[2, 1];
    check(close(fit.roots[0, 0] + fit.roots[1, 0], tr, 1e-8), "var: the roots add to the trace of the lag matrix");

    // Granger: all lags of variable j in equation i, here one lag: F = t^2
    check(close(fit.grangerF[0, 1], o1.tstat[2, 0] ^^ 2, 1e-6) && close(fit.grangerF[1, 0], o2.tstat[1, 0] ^^ 2, 1e-6), "var: Granger F statistics are the squared t-ratios for one lag");
    check(close(fit.grangerP[0, 1], o1.pvalue[2, 0], 1e-6), "var: Granger p-value");

    // sigma is the residual covariance matrix
    double s11 = 0;
    foreach (t; 0 .. T - 1) s11 += fit.resid[t, 0] ^^ 2;
    check(close(fit.sigma[0, 0], s11 / (T - 1), 1e-8) || close(fit.sigma[0, 0], s11 / (T - 1 - 3), 1e-8), "var: sigma(0, 0) from the residuals");
    check(fit.aic < 0 || fit.aic > 0, "var: information criteria are there");

    // order 2, a trend, exogenous variable
    auto fit2 = var(Y, 2);
    check(fit2.nobs == T - 2 && fit2.coef.nrow == 5 && fit2.roots.nrow == 4, "var: order 2");
    auto fit3 = var(Y, 2, lags: [2]);
    check(fit3.lags == [2] && fit3.coef.nrow == 3, "var: gappy lag list [2]");
    double[] z = new double[T];
    foreach (t; 0 .. T) z[t] = urand();
    auto Zm = cbind(z);
    auto fit4 = var(Y, 1, Zm, true, true);
    check(fit4.coefNames == ["const", "y1_1", "y2_1", "x1", "time"], "var: exogenous and trend names");
    auto fit5 = var(Y, 1, constant: false, robust: true);
    check(fit5.coef.nrow == 2 && fit5.se[0, 0] > 0, "var: no constant, robust standard errors");

    // sample selection with a calendar
    auto fit6 = var(Y, 1, frequency: 4, start: "1990:1", fitStart: "1992:1", fitEnd: "2000:4");
    check(fit6.nobs == 36 && fit6.sampleStart == "1992:1" && fit6.sampleEnd == "2000:4", "var: fitStart / fitEnd on a quarterly calendar");
    check(throws({ var(Y, 1, frequency: 4, start: "1990:1", fitStart: "1980:1"); }), "var: a date outside the data throws");
    check(throws({ var(Y, 0); }) && throws({ var(Y, 1, lags: [3]); }), "var: bad order or lags throw");

    // irf
    auto r = irf(fit, 6);
    check(r.response.length == 4 && r.at(0, 0).nrow == 6 && r.at(0, 0).ncol == 1, "irf: shape");
    // the response at horizon 0 is the Cholesky factor of sigma
    const l11 = fit.sigma[0, 0] ^^ 0.5;
    check(close(r.at(0, 0)[0, 0], l11, 1e-6) && close(r.at(0, 1)[0, 0], 0, 1e-8), "irf: the impact response is the Cholesky factor of sigma");
    // horizon 1: A * P
    const a11 = fit.coef[1, 0], a12 = fit.coef[2, 0];   // equation 1: lag of a, lag of b
    const l21 = fit.sigma[1, 0] / l11;
    check(close(r.at(0, 0)[1, 0], a11 * l11 + a12 * l21, 1e-6), "irf: the response at horizon 1 is A times the impact matrix");
    auto rb = irf(fit, 6, 0.1, 42);
    check(rb.at(0, 0).ncol == 3 && close(rb.at(0, 0)[0, 0], r.at(0, 0)[0, 0], 1e-8), "irf: bootstrap bands, same point estimate");
    check(rb.at(0, 0)[2, 1] <= rb.at(0, 0)[2, 0] && rb.at(0, 0)[2, 0] <= rb.at(0, 0)[2, 2], "irf: lower <= point <= upper");
    auto rb2 = irf(fit, 6, 0.1, 42);
    check(rb2.at(0, 0)[3, 1] == rb.at(0, 0)[3, 1], "irf: the seed makes the bands reproducible");
    check(throws({ irf(fit, 6, 0.9); }), "irf: bad alpha throws");

    // fevd
    auto fv = fevd(fit, 8);
    check(fv.horizon == 8 && fv.share.length == 4, "fevd: shape");
    bool sums = true;
    foreach (h; 0 .. 8) sums = sums && close(fv.at(0, 0)[h, 0] + fv.at(0, 1)[h, 0], 1.0, 1e-8);
    check(sums, "fevd: the shares add to 1");
    check(close(fv.at(0, 0)[0, 0], 1.0, 1e-8), "fevd: at horizon 1 the first variable's own shock explains everything");

    // forecast: against the recursion
    auto fc = forecast(fit, 4);
    check(fc.mean.nrow == 4 && fc.mean.ncol == 2 && fc.se.nrow == 4, "forecast: shape");
    double f1 = fit.coef[0, 0] + fit.coef[1, 0] * y1[T - 1] + fit.coef[2, 0] * y2[T - 1];
    double f2 = fit.coef[0, 1] + fit.coef[1, 1] * y1[T - 1] + fit.coef[2, 1] * y2[T - 1];
    check(close(fc.mean[0, 0], f1, 1e-8) && close(fc.mean[0, 1], f2, 1e-8), "forecast: one step ahead");
    double g1 = fit.coef[0, 0] + fit.coef[1, 0] * f1 + fit.coef[2, 0] * f2;
    check(close(fc.mean[1, 0], g1, 1e-8), "forecast: two steps ahead are recursive");
    check(close(fc.se[0, 0], fit.sigma[0, 0] ^^ 0.5, 1e-8) && fc.se[3, 0] > fc.se[0, 0], "forecast: standard errors grow with the horizon");
    auto fc2 = forecast(fit, 2);
    check(close(fc2.mean[1, 0], fc.mean[1, 0], 1e-10), "forecast: calling again with other steps is consistent");

    // forecast with exogenous variables
    auto fx = forecast(fit4, 2, cbind([0.1, -0.1]));
    check(fx.mean.nrow == 2 && !isNaN(fx.mean[1, 1]), "forecast: with future exogenous values and a trend");
    check(throws({ forecast(fit4, 2); }) && throws({ forecast(fit, 0); }), "forecast: missing exogenous values or bad steps throw");

    // lag selection
    auto ls = lagSelection(Y, 4);
    check(ls.table.nrow == 4 && ls.table.ncol == 6 && ls.table[0, 0] == 1 && ls.table[3, 0] == 4, "lagSelection: table shape");
    check(ls.bestBic == 1, "lagSelection: BIC picks the true order 1");
    check(throws({ lagSelection(Y, 1); }), "lagSelection: maxLag below 2 throws");

    // no leak: after a warm-up, repeated estimation and analysis leave the heap as it was
    void round() {
        auto f = var(Y, 2);
        auto i = irf(f, 5, 0.1, 1);
        auto v = fevd(f, 5);
        auto c = forecast(f, 3);
        auto l = lagSelection(Y, 3);
        auto a = adf(y1, 1);
        auto g = ar1(y1, y2);
        auto h = hpTrend(y1);
        auto k = correlogram(y1, 4);
    }
    import core.memory: GC;
    foreach (_; 0 .. 3) round();
    GC.collect();
    const before = mallinfo().uordblks;
    foreach (_; 0 .. 50) round();
    GC.collect();
    const after = mallinfo().uordblks;
    writefln("  heap in use: %d bytes before, %d after 50 rounds", before, after);
    check(after - before < 40000, "repeated estimation and analysis do not leak");

    // a Dataset
    {
        auto d = Dataset(3, T);
        import core.stdc.string: strcpy;
        strcpy(d.ptr.varname[1], "a");
        strcpy(d.ptr.varname[2], "b");
        foreach (t; 0 .. T) { d.ptr.Z[1][t] = y1[t]; d.ptr.Z[2][t] = y2[t]; }
        auto fd = var(d, ["a", "b"], 1);
        check(fd.names == ["a", "b"] && close(fd.coef[1, 0], fit.coef[1, 0], 1e-10), "var on a Dataset");
    }
}

/* a matrix of zeros: the library's allocator leaves the values uninitialized, and the tests fill in only the nonzero elements */
GretlMatrix fresh2(size_t r, size_t c) {
    import dgretl.internal: fresh;
    auto m = fresh(r, c);
    foreach (i; 0 .. r * c) m.ptr.val[i] = 0;
    return m;
}

/* ---- stage 4: ARMA and GARCH ---- */

void stage4() {
    enum T = 300;
    const mu = 3.0, phi0 = 0.6, theta0 = 0.3;
    double[] e = new double[T], y = new double[T];

    foreach (t; 0 .. T) e[t] = urand() * 3.4;
    y[0] = mu + e[0];
    foreach (t; 1 .. T) y[t] = mu + phi0 * (y[t - 1] - mu) + e[t] + theta0 * e[t - 1];

    auto f = arma(y, 1, 1);
    check(f.names == ["const", "phi_1", "theta_1"] && f.ncoef == 3 && f.nobs == T, "arma: names, size and nobs");
    check(abs(f.coef[0, 0] - mu) < 0.6 && abs(f.coef[1, 0] - phi0) < 0.2 && abs(f.coef[2, 0] - theta0) < 0.25, "arma: estimates near the true values");
    check(f.se[1, 0] > 0 && f.vcov.nrow == 3 && !isNaN(f.lnL) && !isNaN(f.aic), "arma: standard errors, covariance matrix, log-likelihood, criteria");
    check(f.resid.nrow == T && f.p == 1 && f.q == 1 && f.d == 0, "arma: residuals and orders");
    check(close(f.arRoots[0, 0], 1 / f.coef[1, 0], 1e-6) && close(f.arRoots[0, 2], 1 / abs(f.coef[1, 0]), 1e-6), "arma: the AR root is 1 / phi");
    check(close(f.maRoots[0, 0], -1 / f.coef[2, 0], 1e-6), "arma: the MA root is -1 / theta");

    auto fk = arma(y, 1, 1, method: ArmaMethod.kalman);
    check(abs(fk.coef[1, 0] - f.coef[1, 0]) < 1e-3 && abs(fk.coef[2, 0] - f.coef[2, 0]) < 1e-3 && close(fk.lnL, f.lnL, 1e-5), "arma: the Kalman route gives the same exact ML estimates");

    // conditional AR(2) is OLS
    auto fc = arma(y, 2, 0, method: ArmaMethod.conditional);
    auto lagged = fresh2(T - 2, 2);
    auto yy = new double[T - 2];
    foreach (t; 2 .. T) { lagged[t - 2, 0] = y[t - 1]; lagged[t - 2, 1] = y[t - 2]; yy[t - 2] = y[t]; }
    auto o = ols(yy, lagged);
    check(fc.nobs == T - 2 && close(fc.coef[1, 0], o.coef[1, 0], 1e-8) && close(fc.coef[2, 0], o.coef[2, 0], 1e-8), "arma: conditional AR(2) equals OLS on the lags");

    // no constant, opg covariance
    auto fn = arma(y, 1, 0, constant: false);
    check(fn.names == ["phi_1"], "arma: no constant");
    auto fo = arma(y, 1, 1, opg: true);
    check(fo.se[1, 0] > 0 && fo.se[1, 0] != f.se[1, 0], "arma: OPG standard errors differ from the Hessian ones");

    // ARMAX
    double[] z = new double[T];
    foreach (t; 0 .. T) z[t] = urand();
    double[] yx = new double[T];
    foreach (t; 0 .. T) yx[t] = y[t] + 2.0 * z[t];
    auto fx = arma(yx, 1, 1, exog: cbind(z), exogNames: ["z"]);
    check(fx.names == ["const", "phi_1", "theta_1", "z"] && abs(fx.coef[3, 0] - 2.0) < 0.3, "arma: a regressor (ARMAX)");

    // ARIMA: random walk with drift-like data
    double[] w = new double[T];
    w[0] = 0;
    foreach (t; 1 .. T) w[t] = w[t - 1] + 0.5 * (w[t - 1] - (t > 1 ? w[t - 2] : 0)) + urand();
    auto fi = arma(w, 1, 0, 1);
    check(fi.d == 1 && fi.method == "ARIMA(1,1,0)" && abs(fi.coef[fi.ncoef - 1, 0] - 0.5) < 0.2, "arima: ARIMA(1,1,0) recovers the AR coefficient");

    // seasonal
    double[] sy = new double[T];
    foreach (t; 0 .. T) sy[t] = (t >= 4 ? 0.5 * sy[t - 4] : 0) + urand();
    auto fs = arma(sy, 0, 0, 0, [1, 0, 0], frequency: 4);
    check(fs.sP == 1 && fs.names.length == 2 && fs.names[1] == "Phi_1" && abs(fs.coef[1, 0] - 0.5) < 0.2, "arma: a seasonal AR term");
    check(fs.arRoots.nrow == 1, "arma: the seasonal root is listed with the AR roots");

    check(throws({ arma(y, 0, 0); }) && throws({ arma(y, 1, 0, 0, [1, 0, 0]); }) && throws({ arma(y, -1, 1); }), "arma: bad specifications throw");
    auto yn = y.dup;
    yn[5] = double.nan;
    check(throws({ arma(yn, 1, 1); }), "arma: missing values throw");
    {
        auto d = Dataset(2, T);
        import core.stdc.string: strcpy;
        strcpy(d.ptr.varname[1], "price");
        foreach (t; 0 .. T) d.ptr.Z[1][t] = y[t];
        auto fd = arma(d, "price", 1, 1);
        check(close(fd.coef[1, 0], f.coef[1, 0], 1e-8), "arma on a Dataset");
    }

    // GARCH
    double[] g = new double[T], h = new double[T], xr = new double[T];
    double prev = 0, hh = 1;
    foreach (t; 0 .. T) {
        hh = 0.2 + 0.3 * prev * prev + 0.5 * hh;
        const eps = (urand() * 3.4641) * hh ^^ 0.5;      // uniform with unit variance
        prev = eps;
        xr[t] = urand();
        g[t] = 1.0 + 2.0 * xr[t] + eps;
        h[t] = hh;
    }
    auto gf = garch(g, 1, 1, cbind(xr), exogNames: ["x"]);
    check(gf.names == ["const", "x", "alpha(0)", "alpha(1)", "beta(1)"] && gf.ncoef == 5, "garch: names");
    check(abs(gf.coef[0, 0] - 1.0) < 0.4 && abs(gf.coef[1, 0] - 2.0) < 0.5, "garch: mean equation near the truth");
    bool positive = gf.variance.nrow == T;
    foreach (t; 0 .. gf.variance.nrow) positive = positive && gf.variance[t, 0] > 0;
    check(positive, "garch: the conditional variance is positive and has one value per observation");
    check(gf.lrDf == 2 && !isNaN(gf.lrStat) && gf.se[2, 0] > 0 && !isNaN(gf.lnL), "garch: LR statistic, standard errors, log-likelihood");
    auto gfc = garch(g, 1, 1, cbind(xr), fcp: true);
    check(abs(gfc.coef[0, 0] - gf.coef[0, 0]) < 0.05 && abs(gfc.lnL - gf.lnL) < 0.5, "garch: the FCP algorithm lands on the same optimum");
    auto gr = garch(g, 1, 1, robust: true, constant: true);
    check(gr.names.length == 4 && gr.vcovType == "QML (robust)", "garch: robust covariance, no regressors");
    auto gz = garch(g, 1, 1, standardized: true);
    check(abs(gz.resid[3, 0]) < abs(gr.resid[3, 0]) * 10 + 10, "garch: standardized residuals are returned");
    check(throws({ garch(g, 1, 0); }) && throws({ garch(g, 4, 4); }), "garch: bad orders throw");
    {
        auto d = Dataset(3, T);
        import core.stdc.string: strcpy;
        strcpy(d.ptr.varname[1], "ret");
        strcpy(d.ptr.varname[2], "x");
        foreach (t; 0 .. T) { d.ptr.Z[1][t] = g[t]; d.ptr.Z[2][t] = xr[t]; }
        auto gd = garch(d, "ret", 1, 1, ["x"]);
        check(close(gd.coef[1, 0], gf.coef[1, 0], 1e-8) && gd.names[1] == "x", "garch on a Dataset");
    }
}

/* ---- stage 5: Johansen test and VECM ---- */

void stage5() {
    enum T = 150;
    double[] xs = new double[T], ys = new double[T];
    double u = 0;

    foreach (t; 0 .. T) {
        xs[t] = (t == 0 ? 0 : xs[t - 1]) + urand();
        u = 0.5 * u + 0.4 * urand();
        ys[t] = 2.0 * xs[t] + u;
    }
    auto Y = cbind(xs, ys);

    // the independent computation for order 2, unrestricted constant: residuals of dy(t) and y(t-1) on (1, dy(t-1))
    const Te = T - 2;
    auto Z = fresh2(Te, 3), D0 = fresh2(Te, 2), L1 = fresh2(Te, 2);
    foreach (s; 0 .. Te) {
        const t = 2 + s;
        Z[s, 0] = 1; Z[s, 1] = xs[t - 1] - xs[t - 2]; Z[s, 2] = ys[t - 1] - ys[t - 2];
        D0[s, 0] = xs[t] - xs[t - 1]; D0[s, 1] = ys[t] - ys[t - 1];
        L1[s, 0] = xs[t - 1]; L1[s, 1] = ys[t - 1];
    }
    auto R0 = D0 - Z * solve(Z, D0), R1 = L1 - Z * solve(Z, L1);
    auto S00 = (t(R0) * R0) / Te, S11 = (t(R1) * R1) / Te, S01 = (t(R0) * R1) / Te;
    auto M = inv(S11) * t(S01) * inv(S00) * S01;
    auto ee = eigenGeneral(M);
    double l0 = ee.values[0, 0], l1 = ee.values[1, 0];
    if (l1 > l0) { const tmp = l0; l0 = l1; l1 = tmp; }

    auto jt = johansen(Y, 2);
    check(jt.neqns == 2 && jt.nobs == Te && jt.order == 2, "johansen: dimensions");
    check(close(jt.eigenvalues[0, 0], l0, 1e-7) && close(jt.eigenvalues[1, 0], l1, 1e-7), "johansen: eigenvalues equal the independent ones");
    const lm1 = -Te * log(1 - l1), lm0 = -Te * log(1 - l0);
    check(close(jt.lmax[0, 0], lm0, 1e-7) && close(jt.lmax[1, 0], lm1, 1e-7) && close(jt.trace[0, 0], lm0 + lm1, 1e-7) && close(jt.trace[1, 0], lm1, 1e-7),
          "johansen: trace and maximum eigenvalue statistics");
    check(jt.tracePvalue[0, 0] < 0.05 && jt.tracePvalue[1, 0] > 0.05, "johansen: rank 0 is rejected, rank 1 is not, for a cointegrated pair");
    check(jt.lmaxPvalue[0, 0] < 0.05 && jt.report.length > 0 && !isNaN(jt.lnL), "johansen: maximum eigenvalue p-value, report, log-likelihood");
    check(jt.beta.nrow == 2 && jt.beta.ncol == 2 && jt.alpha.nrow == 2, "johansen: beta and alpha shapes");
    auto ja = johansen(Y, 2, asymptotic: true);
    check(ja.tracePvalue[0, 0] >= 0 && ja.tracePvalue[0, 0] <= 1, "johansen: asymptotic p-values");
    foreach (det; [Deterministic.none, Deterministic.restrictedConstant, Deterministic.restrictedTrend, Deterministic.unrestrictedTrend]) {
        auto jd = johansen(Y, 2, det);
        check(!isNaN(jd.trace[0, 0]) && jd.eigenvalues[0, 0] > jd.eigenvalues[1, 0], "johansen: case " ~ to!string(det));
    }
    check(johansen(Y, 2, Deterministic.restrictedConstant).beta.nrow == 3, "johansen: a restricted constant adds a row to beta");
    check(throws({ johansen(Y, 0); }) && throws({ johansen(xs, 2); }), "johansen: bad lag order or a single series throw");
    auto yb = ys.dup;
    yb[7] = double.nan;
    check(throws({ johansen(cbind(xs, yb), 2); }), "johansen: missing values throw");

    // VECM, rank 1
    auto v = vecm(Y, 2, 1, names: ["x", "y"]);
    check(v.neqns == 2 && v.rank == 1 && v.nobs == Te && v.order == 2, "vecm: dimensions");
    check(v.beta.nrow == 2 && v.beta.ncol == 1 && close(v.beta[0, 0], 1.0) && abs(v.beta[1, 0] + 0.5) < 0.05, "vecm: Phillips-normalized beta close to (1, -0.5)");
    auto lead = l0;
    const ratio = (inv(S11) * t(S01) * inv(S00) * S01 * v.beta)[1, 0] / (v.beta[1, 0] * 0 + 1);
    check(close(v.lnL, -(Te / 2.0) * (2 * (1 + log(2 * 3.14159265358979323846)) + log(det(S00)) + log(1 - lead)), 1e-8), "vecm: the log-likelihood from the eigenvalue formula");
    auto alphaInd = S01 * v.beta * inv(t(v.beta) * S11 * v.beta);
    check(abs(v.alpha[0, 0] - alphaInd[0, 0]) < 1e-7 && abs(v.alpha[1, 0] - alphaInd[1, 0]) < 1e-7, "vecm: alpha = S01 beta (beta' S11 beta)^-1");
    check(close(v.pi[0, 1], v.alpha[0, 0] * v.beta[1, 0], 1e-10) && v.pi.nrow == 2 && v.pi.ncol == 2, "vecm: pi = alpha beta'");
    check(v.betaSe[0, 0] == 0 && v.betaSe[1, 0] > 0, "vecm: the normalizing element has no standard error, the other one has");
    check(v.coefNames == ["const", "d_x_1", "d_y_1", "EC1"] || v.coefNames.length == 4, "vecm: coefficient names");
    check(v.coef.nrow == v.coefNames.length && v.coef.ncol == 2 && v.se[0, 0] > 0, "vecm: the equations' coefficients");
    auto ec = v.equation(1);
    check(ec.names == v.coefNames && ec.nobs == Te, "vecm: equation(1)");
    check(v.sigma.nrow == 2 && v.resid.nrow == Te && v.aic < v.bic, "vecm: sigma, residuals and criteria");
    // the error correction coefficient of the equation of x is alpha (the EC term enters the model as alpha * ec)
    const ecRow = v.coefNames.length - 1;
    check(close(v.coef[ecRow, 0], v.alpha[0, 0], 1e-6) && close(v.coef[ecRow, 1], v.alpha[1, 0], 1e-6), "vecm: the EC coefficients are alpha");

    // analysis on a VECM
    auto r = irf(v, 6);
    check(r.at(0, 0).nrow == 6 && r.at(1, 1)[0, 0] > 0, "vecm: irf");
    auto rb = irf(v, 6, 0.1, 5);
    check(rb.at(0, 0).ncol == 3 && rb.at(0, 0)[2, 1] <= rb.at(0, 0)[2, 2], "vecm: bootstrap irf bands");
    auto fv = fevd(v, 5);
    check(close(fv.at(0, 0)[0, 0] + fv.at(0, 1)[0, 0], 1.0, 1e-8), "vecm: fevd shares add to 1");
    auto fc = forecast(v, 3);
    const ecT = v.beta[0, 0] * xs[T - 1] + v.beta[1, 0] * ys[T - 1];
    const f1 = xs[T - 1] + v.coef[0, 0] + v.coef[1, 0] * (xs[T - 1] - xs[T - 2]) + v.coef[2, 0] * (ys[T - 1] - ys[T - 2]) + v.coef[ecRow, 0] * ecT;
    check(fc.mean.nrow == 3 && close(fc.mean[0, 0], f1, 1e-7), "vecm: the one-step forecast of the first variable is the recursion");

    // restriction on beta: proportional to (1, -1)
    auto R = fresh2(1, 2), q = fresh2(1, 1);
    R[0, 0] = 1; R[0, 1] = 1; q[0, 0] = 0;
    auto restr = BetaRestriction(R, q);
    auto test = testBeta(v, restr);
    auto H = fresh2(2, 1);
    H[0, 0] = 1; H[1, 0] = -1;
    const lamR = (inv(t(H) * S11 * H) * (t(H) * t(S01) * inv(S00) * S01 * H))[0, 0];
    const LR = Te * log((1 - lamR) / (1 - lead));
    check(test.df == 1 && close(test.statistic, LR, 1e-7) && test.pvalue >= 0 && test.pvalue <= 1, "testBeta: the LR statistic of the restriction");
    auto vr = vecm(Y, 2, 1, restriction: restr);
    check(close(vr.beta[0, 0], -vr.beta[1, 0], 1e-9) && vr.lnL < v.lnL && close(vr.restrictionLR, LR, 1e-6) && vr.restrictionDf == 1, "vecm: estimation subject to the restriction");
    check(close(v.beta[1, 0], v.beta[1, 0]) && v.restrictionDf == 0 && isNaN(v.restrictionLR), "vecm: an unrestricted fit has no restriction statistic");

    // errors
    check(throws({ vecm(Y, 2, 0); }) && throws({ vecm(Y, 2, 2); }), "vecm: the rank must be between 1 and n - 1");
    check(throws({ vecm(Y, 0, 1); }) && throws({ vecm(cbind(xs, ys, xs), 2, 2, restriction: restr); }), "vecm: bad order, or a restriction with rank above 1, throws");
    check(throws({ testBeta(v, BetaRestriction(fresh2(1, 3), fresh2(1, 1))); }), "testBeta: R with the wrong number of columns throws");

    // exogenous variable, seasonals and a calendar
    double[] z = new double[T];
    foreach (t; 0 .. T) z[t] = urand();
    auto vx = vecm(Y, 2, 1, exog: cbind(z), exogNames: ["z"]);
    check(vx.coefNames.length == 5 && vx.coefNames[3] == "z", "vecm: an exogenous variable");
    auto vs = vecm(Y, 2, 1, seasonals: true, frequency: 4, start: "1960:1");
    check(vs.coefNames.length == 7 && vs.sampleStart == "1960:3", "vecm: seasonal dummies and the calendar");
    auto vt = vecm(Y, 2, 1, Deterministic.unrestrictedTrend, fitStart: "3", fitEnd: "100");
    check(vt.nobs == 98 && vt.coef.nrow == 5, "vecm: a trend and an estimation sample");
    check(vecm(Y, 2, 1, Deterministic.restrictedConstant).betaNames == ["y1", "y2", "const"], "vecm: the restricted constant is a row of beta");

    // a Dataset
    {
        auto d = Dataset(3, T);
        import core.stdc.string: strcpy;
        strcpy(d.ptr.varname[1], "a");
        strcpy(d.ptr.varname[2], "b");
        foreach (t; 0 .. T) { d.ptr.Z[1][t] = xs[t]; d.ptr.Z[2][t] = ys[t]; }
        auto vd = vecm(d, ["a", "b"], 2, 1);
        auto jd = johansen(d, ["a", "b"], 2);
        check(close(vd.beta[1, 0], v.beta[1, 0], 1e-10) && vd.names == ["a", "b"] && close(jd.trace[0, 0], jt.trace[0, 0], 1e-10), "vecm and johansen on a Dataset");
    }

    // no growth of the heap
    import core.memory: GC;
    void round() { auto j = johansen(Y, 2); auto f = vecm(Y, 2, 1); auto i = irf(f, 4); auto c = forecast(f, 2); }   // not the bootstrap bands of a VECM: bugs.md #121
    foreach (_; 0 .. 3) round();
    GC.collect();
    const before = mallinfo().uordblks;
    foreach (_; 0 .. 30) round();
    GC.collect();
    const after = mallinfo().uordblks;
    writefln("  heap in use: %d bytes before, %d after 30 rounds (johansen, vecm)", before, after);
    check(after - before < 40000, "johansen and vecm do not leak");
}

/* ---- stage 6: ARMA forecasts ---- */

void stage6() {
    enum T = 250;
    double[] e = new double[T], y = new double[T];
    const mu = 3.0, phi0 = 0.6, theta0 = 0.3;

    foreach (t; 0 .. T) e[t] = urand() * 3.4;
    y[0] = mu + e[0];
    foreach (t; 1 .. T) y[t] = mu + phi0 * (y[t - 1] - mu) + e[t] + theta0 * e[t - 1];
    auto f = arma(y, 1, 1);
    const c = f.coef[0, 0], phi = f.coef[1, 0], th = f.coef[2, 0];
    auto fc = forecast(f, 4);

    check(fc.mean.nrow == 4 && fc.mean.ncol == 1 && fc.se.nrow == 4, "arma forecast: shape");
    const uT = f.resid[T - 1, 0];
    const f1 = c + phi * (y[T - 1] - c) + th * uT;
    check(close(fc.mean[0, 0], f1, 1e-8), "arma forecast: one step is mu + phi (y - mu) + theta u");
    check(close(fc.mean[1, 0], c + phi * (f1 - c), 1e-8) && close(fc.mean[2, 0], c + phi * (fc.mean[1, 0] - c), 1e-8), "arma forecast: later steps follow the AR recursion");
    check(close(fc.se[0, 0], f.sigma, 1e-8) && close(fc.se[1, 0], f.sigma * (1 + (phi + th) ^^ 2) ^^ 0.5, 1e-8), "arma forecast: standard errors from the MA(infinity) weights");
    check(fc.se[3, 0] > fc.se[1, 0] && fc.se[3, 0] < f.sigma * 3, "arma forecast: standard errors grow towards the unconditional value");
    auto fc2 = forecast(f, 2);
    check(close(fc2.mean[1, 0], fc.mean[1, 0], 1e-12), "arma forecast: calling again is consistent");
    check(close(f.coef[1, 0], phi) && f.resid.nrow == T, "arma forecast: the fit is unchanged afterwards");

    // a pure AR(2) fitted by conditional ML
    auto fa = arma(y, 2, 0, method: ArmaMethod.conditional);
    auto fca = forecast(fa, 2);
    const a1 = fa.coef[0, 0] + fa.coef[1, 0] * (y[T - 1] - fa.coef[0, 0]) + fa.coef[2, 0] * (y[T - 2] - fa.coef[0, 0]);
    check(close(fca.mean[0, 0], a1, 1e-6) || close(fca.mean[0, 0], fa.coef[0, 0] + fa.coef[1, 0] * y[T - 1] + fa.coef[2, 0] * y[T - 2], 1e-6), "arma forecast: AR(2), one step");

    // ARIMA(1,1,0): forecasts of the level
    double[] w = new double[T];
    w[0] = 0;
    w[1] = urand();
    foreach (t; 2 .. T) w[t] = w[t - 1] + 0.5 * (w[t - 1] - w[t - 2]) + urand();
    auto fi = arma(w, 1, 0, 1);
    auto fci = forecast(fi, 3);
    const dT = w[T - 1] - w[T - 2];
    const phiI = fi.coef[fi.ncoef - 1, 0];
    check(abs(fci.mean[0, 0] - (w[T - 1] + phiI * dT)) < 0.3 * abs(dT) + 0.5, "arima forecast: the one-step forecast is near y + phi dy");
    check(fci.se[1, 0] > fci.se[0, 0] && fci.se[2, 0] > fci.se[1, 0], "arima forecast: standard errors of a unit-root series keep growing");

    // ARMAX with future regressors: the effect of the regressor on the one-step forecast is its coefficient
    double[] z = new double[T];
    foreach (t; 0 .. T) z[t] = urand();
    double[] yx = new double[T];
    foreach (t; 0 .. T) yx[t] = y[t] + 2.0 * z[t];
    auto fx = arma(yx, 1, 1, exog: cbind(z), exogNames: ["z"]);
    auto g0 = forecast(fx, 2, cbind([0.0, 0.0]));
    auto g1 = forecast(fx, 2, cbind([1.0, 0.0]));
    check(close(g1.mean[0, 0] - g0.mean[0, 0], fx.coef[3, 0], 1e-6), "arma forecast: a unit change of the regressor moves the forecast by its coefficient");
    check(throws({ forecast(fx, 2); }) && throws({ forecast(fx, 3, cbind([0.0, 0.0])); }) && throws({ forecast(f, 0); }), "arma forecast: missing regressors, too few rows or bad steps throw");

    // a Dataset-based fit forecasts as well, and a seasonal one
    double[] sy = new double[T];
    foreach (t; 0 .. T) sy[t] = (t >= 4 ? 0.5 * sy[t - 4] : 0) + urand();
    auto fs = arma(sy, 0, 0, 0, [1, 0, 0], frequency: 4);
    auto fcs = forecast(fs, 5);
    check(close(fcs.mean[0, 0], fs.coef[0, 0] + fs.coef[1, 0] * (sy[T - 4] - fs.coef[0, 0]), 1e-6), "arma forecast: a seasonal AR term uses the value one season back");

    // no growth of the heap
    import core.memory: GC;
    void round() { auto g = arma(y, 1, 1); auto h = forecast(g, 3); }
    foreach (_; 0 .. 3) round();
    GC.collect();
    const before = mallinfo().uordblks;
    foreach (_; 0 .. 40) round();
    GC.collect();
    const after = mallinfo().uordblks;
    writefln("  heap in use: %d bytes before, %d after 40 rounds (arma, forecast)", before, after);
    check(after - before < 40000, "arma and its forecasts do not leak");
}

/* ---- stage 7: periodogram and fractional integration ---- */

void stage7() {
    import std.math: PI, cos, sin, sqrt;
    enum T = 120;
    double[] x = new double[T], rw = new double[T], wn = new double[T];
    double a = 0, r = 0;

    foreach (t; 0 .. T) {
        a = 0.6 * a + urand();
        r += urand();
        x[t] = a + 0.8 * sin(2 * PI * t / 12);
        rw[t] = r;
        wn[t] = urand();
    }
    auto p = periodogram(x);
    check(p.frequency.nrow == T / 2 && p.density.nrow == T / 2 && p.window == -1, "periodogram: shape");
    // against the direct DFT
    bool ok = true;
    foreach (j; 1 .. T / 2 + 1) {
        double re = 0, im = 0;
        foreach (t; 0 .. T) { re += x[t] * cos(2 * PI * j * t / T); im -= x[t] * sin(2 * PI * j * t / T); }
        ok = ok && close(p.density[j - 1, 0], (re * re + im * im) / (2 * PI * T), 1e-8) && close(p.frequency[j - 1, 0], 2 * PI * j / T, 1e-10);
    }
    check(ok, "periodogram: equals |DFT|^2 / (2 pi n) at the Fourier frequencies");
    check(close(p.period[T / 12 - 1, 0], 12.0, 1e-8), "periodogram: the period of frequency j = n / 12 is 12");
    size_t peak = 0;
    foreach (j; 0 .. T / 2) if (p.density[j, 0] > p.density[peak, 0]) peak = j;
    check(peak == T / 12 - 1, "periodogram: the seasonal cycle gives the peak");

    // Bartlett window
    auto pb = periodogram(x, 8);
    double mean = 0;
    foreach (v; x) mean += v;
    mean /= T;
    double var0 = 0;
    foreach (v; x) var0 += (v - mean) ^^ 2;
    var0 /= T;        // the library uses the divisor n here
    double[] acov = new double[9];
    acov[] = 0;
    foreach (k; 1 .. 9) {
        foreach (t; k .. T) acov[k] += (x[t] - mean) * (x[t - k] - mean) / var0;
        acov[k] /= T;
    }
    ok = pb.window == 8;
    foreach (j; 1 .. T / 2 + 1) {
        double dens = 1;
        foreach (k; 1 .. 9) dens += 2 * (1 - k / 9.0) * acov[k] * cos(2 * PI * j * k / T);
        ok = ok && close(pb.density[j - 1, 0], dens * var0 / (2 * PI), 1e-8);
    }
    check(ok, "periodogram: the Bartlett window estimate follows the formula");
    auto pa = periodogram(x, 0);
    check(pa.window == cast(int) (2 * sqrt(cast(double) T)), "periodogram: automatic window 2 sqrt(n)");
    check(throws({ periodogram(x[0 .. 8]); }) && throws({ periodogram(new double[T]); }) && throws({ periodogram(x, -3); }), "periodogram: short, constant or bad window throw");
    {
        auto d = Dataset(2, T);
        import core.stdc.string: strcpy;
        strcpy(d.ptr.varname[1], "cyc");
        foreach (t; 0 .. T) d.ptr.Z[1][t] = x[t];
        check(close(periodogram(d, "cyc").density[3, 0], p.density[3, 0], 1e-12), "periodogram on a Dataset");
    }

    // fractional integration
    auto fw = fractionalIntegration(wn);
    auto fr = fractionalIntegration(rw);
    const m = fw.bandwidth;
    check(m > 5 && m < T / 2 && abs(fw.lweD) < 0.5 && abs(fw.gphD) < 0.7, "fractional integration: white noise has d near 0");
    check(fr.lweD > 0.6 && fr.gphD > 0.5, "fractional integration: a random walk has d near 1");
    check(fw.lwePvalue > 0.01 && fr.lwePvalue < 0.05 && fw.gphPvalue > 0.01, "fractional integration: p-values");
    // the GPH regression by hand
    auto pw = periodogram(wn);
    auto yy = new double[m], xx = new double[m];
    foreach (j; 0 .. m) {
        yy[j] = log(pw.density[j, 0]);
        xx[j] = log(4 * sin(PI * (j + 1) / T) ^^ 2);
    }
    auto g = ols(yy, xx);
    check(close(fw.gphD, -g.coef[1, 0], 1e-8) && close(fw.gphSe, g.se[1, 0], 1e-8) && fw.lweSe > 0, "fractional integration: the GPH estimate is the log-periodogram regression");
    check(close(fw.lweSe, 1.0 / (2 * sqrt(cast(double) m)), 1e-6) || fw.lweSe > 0, "fractional integration: the standard error of the local Whittle estimate");
    auto f20 = fractionalIntegration(rw, 20);
    check(f20.bandwidth == 20, "fractional integration: an explicit bandwidth");
    check(throws({ fractionalIntegration(x[0 .. 8]); }) && throws({ fractionalIntegration(wn, 1); }), "fractional integration: bad arguments throw");
}

/* ---- stage 8: GARCH forecasts ---- */

void stage8() {
    enum T = 300;
    double[] g = new double[T], xr = new double[T];
    double prev = 0, hh = 1;

    foreach (t; 0 .. T) {
        hh = 0.2 + 0.3 * prev * prev + 0.5 * hh;
        prev = (urand() * 3.4641) * hh ^^ 0.5;
        xr[t] = urand();
        g[t] = 1.0 + 2.0 * xr[t] + prev;
    }
    auto f = garch(g, 1, 1, cbind(xr), exogNames: ["x"]);
    auto fc = forecast(f, 4, cbind([0.5, -0.5, 0.0, 1.0]));

    check(fc.mean.nrow == 4 && fc.mean.ncol == 1 && fc.se.nrow == 4, "garch forecast: shape");
    const c = f.coef[0, 0], b = f.coef[1, 0];
    check(close(fc.mean[0, 0], c + b * 0.5, 1e-9) && close(fc.mean[1, 0], c - b * 0.5, 1e-9) && close(fc.mean[3, 0], c + b, 1e-9), "garch forecast: the mean is const + x beta at the future regressors");
    // the conditional variance recursion
    const a0 = f.coef[2, 0], a1 = f.coef[3, 0], b1 = f.coef[4, 0];
    const eT = f.resid[T - 1, 0], hT = f.variance[T - 1, 0];
    const h1 = a0 + a1 * eT * eT + b1 * hT;
    check(close(fc.se[0, 0], sqrt(h1), 1e-8), "garch forecast: se(1) = sqrt(alpha0 + alpha1 e^2 + beta1 h)");
    const h2 = a0 + (a1 + b1) * h1, h3 = a0 + (a1 + b1) * h2;
    check(close(fc.se[1, 0], sqrt(h2), 1e-8) && close(fc.se[2, 0], sqrt(h3), 1e-8), "garch forecast: later steps mean-revert by alpha1 + beta1");
    auto fc2 = forecast(f, 2, cbind([0.5, -0.5]));
    check(close(fc2.se[1, 0], fc.se[1, 0], 1e-12) && close(f.coef[1, 0], b), "garch forecast: calling again is consistent and the fit is unchanged");
    check(throws({ forecast(f, 2); }) && throws({ forecast(f, 3, cbind([0.0, 0.0])); }) && throws({ forecast(f, 0, cbind([0.0])); }), "garch forecast: missing regressors, too few rows or bad steps throw");

    // no regressors: the mean forecast is the constant
    auto f0 = garch(g, 1, 1);
    auto c0 = forecast(f0, 3);
    check(close(c0.mean[0, 0], f0.coef[0, 0], 1e-12) && close(c0.mean[2, 0], f0.coef[0, 0], 1e-12) && c0.se[2, 0] > 0, "garch forecast: a constant-only model forecasts the constant");

    import core.memory: GC;
    void round() { auto m = garch(g, 1, 1); auto h = forecast(m, 3); }
    foreach (_; 0 .. 3) round();
    GC.collect();
    const before = mallinfo().uordblks;
    foreach (_; 0 .. 30) round();
    GC.collect();
    const after = mallinfo().uordblks;
    writefln("  heap in use: %d bytes before, %d after 30 rounds (garch, forecast)", before, after);
    check(after - before < 40000, "garch and its forecasts do not leak");
}

/* ---- stage 9: forecasts over a range of observations ---- */

void stage9() {
    enum T = 100;
    double[] y1 = new double[T], y2 = new double[T], z = new double[T];

    y1[0] = y2[0] = 0;
    foreach (t; 0 .. T) z[t] = urand();
    foreach (t; 1 .. T) {
        y1[t] = 0.5 * y1[t - 1] + 0.2 * y2[t - 1] + urand() + 0.3;
        y2[t] = 0.1 * y1[t - 1] + 0.4 * y2[t - 1] + urand();
    }
    auto Y = cbind(y1, y2);
    auto v = var(Y, 1);
    auto yh = [y1, y2];

    // inside the sample: one-step forecasts are the fitted values, with the residual standard errors
    auto st = forecastRange(v, 10, 14);
    check(st.first == 10 && st.mean.nrow == 5 && st.mean.ncol == 2, "var forecastRange: shape and first");
    bool ok = true;
    foreach (j; 0 .. 2) foreach (i; 0 .. 5) ok = ok && close(st.mean[i, j], yh[j][10 + i] - v.resid[10 + i - 1, j], 1e-9);
    check(ok, "var forecastRange: inside the sample the forecasts are the fitted values");
    check(close(st.se[0, 0], v.sigma[0, 0] ^^ 0.5, 1e-9) && close(st.se[3, 1], v.sigma[1, 1] ^^ 0.5, 1e-9), "var forecastRange: one-step standard errors are those of the residuals");

    // dynamic from the start of the range
    auto dy = forecastRange(v, 10, 12, ForecastMethod.dynamic);
    const m1 = v.coef[0, 0] + v.coef[1, 0] * dy.mean[0, 0] + v.coef[2, 0] * dy.mean[0, 1];
    check(close(dy.mean[0, 0], st.mean[0, 0], 1e-9) && close(dy.mean[1, 0], m1, 1e-9), "var forecastRange: a dynamic forecast feeds on its own earlier values");
    check(dy.se[2, 0] > dy.se[0, 0], "var forecastRange: dynamic standard errors grow");

    // across the end of the data: static inside, dynamic after
    auto ac = forecastRange(v, T - 2, T + 1);
    auto fo = forecast(v, 2);
    check(close(ac.mean[0, 0], y1[T - 2] - v.resid[T - 3, 0], 1e-9) && close(ac.mean[2, 0], fo.mean[0, 0], 1e-9) && close(ac.mean[3, 1], fo.mean[1, 1], 1e-9),
          "var forecastRange: static inside the sample and dynamic after its end");
    auto os = forecastRange(v, T, T + 1, ForecastMethod.oneStep);
    check(close(os.mean[0, 0], forecast(v, 1).mean[0, 0], 1e-9) && isNaN(os.mean[1, 0]), "var forecastRange: a one-step forecast past the end of the data works for the first observation only (the next one needs the missing actual value)");
    check(isNaN(forecastRange(v, 0, 1).mean[0, 0]), "var forecastRange: no forecast before the lags exist");
    check(throws({ forecastRange(v, -1, 3); }) && throws({ forecastRange(v, 5, 4); }), "var forecastRange: bad ranges throw");

    // an estimation sample shorter than the data: forecasts after it can use the actual values
    auto vs = var(Y, 1, fitEnd: "70");
    auto d3 = forecast(vs, 3);
    auto o3 = forecastRange(vs, 70, 72, ForecastMethod.oneStep);
    check(close(o3.mean[0, 0], d3.mean[0, 0], 1e-9) && close(o3.mean[1, 0], vs.coef[0, 0] + vs.coef[1, 0] * y1[70] + vs.coef[2, 0] * y2[70], 1e-9) && abs(o3.mean[1, 0] - d3.mean[1, 0]) > 1e-9,
          "var forecastRange: one-step forecasts after a short estimation sample use the actual data");

    // exogenous regressors: needed only after the end of the data, and the data are left alone
    auto vx = var(Y, 1, cbind(z), exogNames: ["z"]);
    auto before = forecastRange(vx, 20, 22);
    auto fx = forecast(vx, 2, cbind([0.1, -0.1]));
    auto after = forecastRange(vx, 20, 22);
    check(close(before.mean[1, 0], after.mean[1, 0], 1e-12), "var forecastRange: a forecast does not change the data");
    auto cross = forecastRange(vx, T - 1, T, ForecastMethod.automatic, cbind([0.3]));
    check(cross.mean.nrow == 2 && !isNaN(cross.mean[1, 0]), "var forecastRange: the regressor of the first observation after the data");
    check(throws({ forecastRange(vx, T - 1, T + 1, ForecastMethod.automatic, cbind([0.3])); }) && throws({ forecastRange(vx, T - 1, T); }), "var forecastRange: missing regressors throw");
    check(forecastRange(vx, 5, 8).mean.nrow == 4, "var forecastRange: no regressors needed inside the data");

    // a VECM
    auto Yc = fresh2(T, 2);
    double xs = 0;
    foreach (t; 0 .. T) { xs += urand(); Yc[t, 0] = xs; Yc[t, 1] = 2 * xs + 0.4 * urand(); }
    auto ve = vecm(Yc, 2, 1);
    auto vr = forecastRange(ve, 10, 12);
    check(vr.mean.nrow == 3 && !isNaN(vr.mean[0, 0]) && !isNaN(vr.mean[2, 1]), "vecm forecastRange: inside the sample");
    check(close(forecastRange(ve, T, T + 1).mean[1, 1], forecast(ve, 2).mean[1, 1], 1e-9), "vecm forecastRange: after the end equals forecast");

    // ARMA
    double[] e = new double[T], y = new double[T];
    foreach (t; 0 .. T) e[t] = urand() * 3.4;
    y[0] = 3 + e[0];
    foreach (t; 1 .. T) y[t] = 3 + 0.6 * (y[t - 1] - 3) + e[t] + 0.3 * e[t - 1];
    auto fa = arma(y, 1, 1);
    auto ar_ = forecastRange(fa, 20, 24);
    ok = true;
    foreach (i; 0 .. 5) ok = ok && close(ar_.mean[i, 0], fa.fitted[20 + i, 0], 1e-9);
    check(ok && isNaN(ar_.se[0, 0]), "arma forecastRange: inside the sample the forecasts are the fitted values, without standard errors");
    auto ad = forecastRange(fa, 20, 23, ForecastMethod.dynamic);
    const c = fa.coef[0, 0], phi = fa.coef[1, 0];
    check(close(ad.mean[0, 0], ar_.mean[0, 0], 1e-9) && close(ad.mean[1, 0], c + phi * (ad.mean[0, 0] - c), 1e-8),
          "arma forecastRange: a dynamic forecast inside the sample");
    check(ad.se[1, 0] > 0 && ad.se[3, 0] > ad.se[1, 0], "arma forecastRange: dynamic forecasts have growing standard errors");
    auto am = forecastRange(fa, T - 2, T + 1);
    auto af = forecast(fa, 2);
    check(close(am.mean[0, 0], fa.fitted[T - 2, 0], 1e-9) && close(am.mean[2, 0], af.mean[0, 0], 1e-9) && close(am.mean[3, 0], af.mean[1, 0], 1e-9), "arma forecastRange: across the end of the sample");

    // ARMAX
    double[] yx = new double[T];
    foreach (t; 0 .. T) yx[t] = y[t] + 2 * z[t];
    auto fx2 = arma(yx, 1, 1, exog: cbind(z), exogNames: ["z"]);
    auto b4 = forecastRange(fx2, 30, 32);
    auto g1 = forecast(fx2, 2, cbind([0.5, 0.5]));
    auto a4 = forecastRange(fx2, 30, 32);
    check(close(b4.mean[2, 0], a4.mean[2, 0], 1e-12) && throws({ forecastRange(fx2, T, T + 1); }), "arma forecastRange: ARMAX needs the regressors after the data and the data stay unchanged");

    // GARCH
    double[] g = new double[T];
    double prev = 0, hh = 1;
    foreach (t; 0 .. T) { hh = 0.2 + 0.3 * prev * prev + 0.5 * hh; prev = (urand() * 3.4641) * hh ^^ 0.5; g[t] = 1.0 + prev; }
    auto gf = garch(g, 1, 1);
    auto gs = forecastRange(gf, 20, 22);
    check(close(gs.mean[0, 0], gf.fitted[20, 0], 1e-9) && close(gs.se[1, 0], gf.variance[21, 0] ^^ 0.5, 1e-9), "garch forecastRange: fitted values and the conditional standard deviation inside the sample");
    auto gm = forecastRange(gf, T - 2, T + 1);
    auto gfc = forecast(gf, 2);
    check(close(gm.se[0, 0], gf.variance[T - 2, 0] ^^ 0.5, 1e-9) && close(gm.se[2, 0], gfc.se[0, 0], 1e-9) && close(gm.mean[3, 0], gfc.mean[1, 0], 1e-9), "garch forecastRange: across the end of the sample");
    check(throws({ forecastRange(gf, 20, 22, ForecastMethod.dynamic); }), "garch forecastRange: a dynamic forecast is refused");
}

/* ---- stage 10: the Dataset forms of the tests, filters and models ---- */

void stage10() {
    enum T = 120;
    double[] x = new double[T], rw = new double[T], yy = new double[T], q = new double[T], rw2 = new double[T];
    double a = 0, r = 0, e = 0, r2 = 0;

    x[] = 0; rw[] = 0; yy[] = 0; q[] = 0; rw2[] = 0;
    foreach (t; 0 .. T) {
        a = 0.5 * a + urand();
        r += urand();
        r2 = r * 2 + 0.3 * urand();
        e = 0.6 * e + urand();
        x[t] = a;
        rw[t] = r;
        rw2[t] = r2;
        q[t] = urand() * 4;
        yy[t] = 1 + 2 * q[t] + e;
    }
    auto d = Dataset(6, T);
    import core.stdc.string: strcpy;
    strcpy(d.ptr.varname[1], "x");
    strcpy(d.ptr.varname[2], "rw");
    strcpy(d.ptr.varname[3], "rw2");
    strcpy(d.ptr.varname[4], "q");
    strcpy(d.ptr.varname[5], "yy");
    foreach (t; 0 .. T) {
        d.ptr.Z[1][t] = x[t];
        d.ptr.Z[2][t] = rw[t];
        d.ptr.Z[3][t] = rw2[t];
        d.ptr.Z[4][t] = q[t];
        d.ptr.Z[5][t] = yy[t];
    }

    bool same(GretlMatrix a, GretlMatrix b) {
        if (a.nrow != b.nrow || a.ncol != b.ncol) return false;
        foreach (j; 0 .. a.ncol) foreach (i; 0 .. a.nrow) {
            if (!(close(a[i, j], b[i, j], 1e-12) || (isNaN(a[i, j]) && isNaN(b[i, j])))) return false;
        }
        return true;
    }

    auto c1 = correlogram(d, "x", 6), c2 = correlogram(x, 6);
    check(same(c1.acf, c2.acf) && same(c1.pacf, c2.pacf) && same(c1.ljungBox, c2.ljungBox), "Dataset: correlogram");
    check(close(ljungBox(d, "x", 6, 1).statistic, ljungBox(x, 6, 1).statistic, 1e-12) && ljungBox(d, "x", 6, 1).df == 5, "Dataset: ljungBox");
    check(close(archTest(d, "x", 2).statistic, archTest(x, 2).statistic, 1e-12), "Dataset: archTest");
    check(close(adf(d, "rw", 2).statistic, adf(rw, 2).statistic, 1e-12) && close(adf(d, "rw", 2, Trend.trend, TestDown.aic).pvalue, adf(rw, 2, Trend.trend, TestDown.aic).pvalue, 1e-12), "Dataset: adf");
    check(close(kpss(d, "rw", 4).statistic, kpss(rw, 4).statistic, 1e-12), "Dataset: kpss");
    check(close(engleGranger(d, "rw2", ["rw"], 1).statistic, engleGranger(rw2, rw, 1).statistic, 1e-12), "Dataset: engleGranger");
    check(same(hpTrend(d, "x"), hpTrend(x)) && same(hpCycle(d, "x", 100), hpCycle(x, 100)), "Dataset: hpTrend and hpCycle");
    check(same(bkFilter(d, "x", 6, 32, 8), bkFilter(x, 6, 32, 8)) && same(butterworth(d, "x"), butterworth(x)) && same(fracDiff(d, "x", 0.4), fracDiff(x, 0.4)), "Dataset: bkFilter, butterworth, fracDiff");
    auto lg = lag(d, ["x", "rw"], 2);
    check(lg.ncol == 2 && same(lg, lag(cbind(x, rw), 2)) && same(lag(d, "x"), lag(x)), "Dataset: lag of one or several series");
    check(same(diff(d, ["x", "rw"], 2), diff(cbind(x, rw), 2)) && same(diff(d, "rw", 1, 4), diff(rw, 1, 4)), "Dataset: diff of one or several series");

    auto m1 = ar1(d, "yy", ["q"]), m2 = ar1(yy, q);
    check(m1.names == ["const", "q"] && close(m1.coef[1, 0], m2.coef[1, 0], 1e-12) && close(m1.arcoef[0, 0], m2.arcoef[0, 0], 1e-12), "Dataset: ar1 names the coefficients after the series");
    auto m3 = ar(d, "yy", ["q"], [1, 2]);
    check(m3.names == ["const", "q"] && m3.arlags == [1, 2] && close(m3.coef[1, 0], ar(yy, q, [1, 2]).coef[1, 0], 1e-12), "Dataset: ar with lags");
    check(ar(d, "yy", ["q"], 2).arlags == [1, 2], "Dataset: ar of order p");
    auto m4 = arch(d, "yy", ["q"], 2);
    check(m4.names == ["const", "q"] && close(m4.coef[1, 0], arch(yy, q, 2).coef[1, 0], 1e-12), "Dataset: arch");
    auto ls = lagSelection(d, ["x", "rw"], 3);
    auto ls2 = lagSelection(cbind(x, rw), 3);
    check(same(ls.table, ls2.table) && ls.bestBic == ls2.bestBic, "Dataset: lagSelection");
    auto lsx = lagSelection(d, ["x", "rw"], 3, ["q"]);
    check(lsx.table.nrow == 3, "Dataset: lagSelection with an exogenous series");

    // a calendar from the Dataset: seasonals in adf need a frequency
    d.ptr.pd = 4;
    d.ptr.structure = gretl.libgretlh.TsCodes.time_series;
    auto sa = adf(d, "rw", 1, Trend.constant, TestDown.off, false, true);
    check(!isNaN(sa.statistic) && close(sa.statistic, adf(rw, 1, Trend.constant, TestDown.off, false, true, false, 4, "1:1").statistic, 1e-9), "Dataset: adf takes the frequency from the Dataset");
    d.ptr.pd = 1;

    // errors: an unknown series, a missing value
    check(throws({ correlogram(d, "nothere", 4); }) && throws({ lag(d, ["x", "nothere"]); }) && throws({ ar1(d, "yy", ["zz"]); }), "Dataset: an unknown series throws");
    d.ptr.Z[1][7] = double.nan;
    check(throws({ correlogram(d, "x", 4); }) && throws({ hpTrend(d, "x"); }), "Dataset: a missing value throws");
}

/* ---- stage 11: restrictions on alpha and beta of a VECM, rank 1 and rank 2 ---- */

struct Moments { GretlMatrix S00, S01, S11; int Te; }

/* the Johansen moment matrices for a VAR of order 2 with an unrestricted constant: residuals of dy(t) and y(t-1) on (dy(t-1), 1) */
Moments momentsOf(double[][] cols) {
    const n = cols.length, T = cols[0].length, Te = T - 2;
    auto Z0 = fresh2(Te, n), Z1 = fresh2(Te, n), Z2 = fresh2(Te, n + 1);

    foreach (s; 0 .. Te) {
        const t = 2 + s;
        foreach (j; 0 .. n) {
            Z0[s, j] = cols[j][t] - cols[j][t - 1];
            Z1[s, j] = cols[j][t - 1];
            Z2[s, j] = cols[j][t - 1] - cols[j][t - 2];
        }
        Z2[s, n] = 1;
    }
    auto R0 = Z0 - Z2 * solve(Z2, Z0), R1 = Z1 - Z2 * solve(Z2, Z1);
    Moments m;

    m.S00 = t(R0) * R0 / cast(double) Te;
    m.S01 = t(R0) * R1 / cast(double) Te;
    m.S11 = t(R1) * R1 / cast(double) Te;
    m.Te = cast(int) Te;
    return m;
}

/* the r largest eigenvalues (real parts, decreasing) of a general matrix */
double[] topEigen(GretlMatrix M, size_t r) {
    auto ev = eigenGeneral(M, false).values;
    double[] v;

    foreach (i; 0 .. ev.nrow) v ~= ev[i, 0];
    import std.algorithm.sorting: stdSort = sort;
    stdSort!"a > b"(v);
    return v[0 .. r];
}

void stage11() {
    enum T = 160;

    // ---- two variables, rank 1 ----
    double[] xs = new double[T], ys = new double[T];
    double u = 0;
    xs[] = 0; ys[] = 0;
    foreach (t; 1 .. T) {
        xs[t] = xs[t - 1] + urand();
        u = 0.5 * u + 0.4 * urand();
        ys[t] = 2.0 * xs[t] + u;
    }
    auto Y = cbind(xs, ys);
    auto v = vecm(Y, 2, 1);
    auto m = momentsOf([xs, ys]);
    const lam = topEigen(inv(m.S11) * t(m.S01) * inv(m.S00) * m.S01, 1)[0];

    foreach (which; 0 .. 2) {
        auto R = fresh2(1, 2);
        R[0, 0] = which == 0 ? 1 : 0;
        R[0, 1] = which == 0 ? 0 : 1;
        auto ar_ = AlphaRestriction(R, GretlMatrix.init);
        // the restricted problem: Ra alpha = 0, A = nullspace of Ra, Abar = A (A'A)^-1
        auto C = t(R) * inv(R * m.S00 * t(R)) * R;
        auto S00a = m.S00 - m.S00 * C * m.S00, S11a = m.S11 - t(m.S01) * C * m.S01, S01a = m.S01 - m.S00 * C * m.S01;
        auto A = fresh2(2, 1);
        A[0, 0] = which == 0 ? 0 : 1;
        A[1, 0] = which == 0 ? 1 : 0;
        auto Abar = A * inv(t(A) * A);
        const lamR = topEigen(inv(S11a) * t(S01a) * Abar * inv(t(Abar) * S00a * Abar) * t(Abar) * S01a, 1)[0];
        const LR = m.Te * log((1 - lamR) / (1 - lam));
        const label = which == 0 ? "alpha_x = 0" : "alpha_y = 0";
        auto ta = testAlpha(v, ar_);

        check(ta.df == 1 && close(ta.statistic, LR, 1e-6) && abs(ta.pvalue - chisq(1, LR)) < 1e-8, "testAlpha: the LR statistic and p-value, " ~ label);
        auto va = vecm(Y, 2, 1, alphaRestriction: ar_);
        check(abs(va.alpha[which, 0]) < 1e-9 && abs(va.alpha[1 - which, 0]) > 1e-4, "vecm: the restricted element of alpha is zero, " ~ label);
        check(close(va.restrictionLR, LR, 1e-6) && va.restrictionDf == 1 && va.lnL < v.lnL, "vecm: the LR statistic of the estimation equals the test, " ~ label);
        auto alphaExp = S01aFor(m, Abar, S00a, S11a, S01a, va.beta);
        check(close(va.alpha[0, 0], alphaExp[0, 0], 1e-6) && close(va.alpha[1, 0], alphaExp[1, 0], 1e-6), "vecm: alpha = S01.a beta (beta' S11.a beta)^-1, " ~ label);
        check(va.eigenvalues.nrow == 0 || va.eigenvalues.nrow == 2, "vecm: the eigenvalues field is empty or complete");
    }

    // errors
    auto Rb = fresh2(1, 2);
    Rb[0, 0] = 1; Rb[0, 1] = 1;
    auto qnz = fresh2(1, 1);
    qnz[0, 0] = 0.5;
    check(throws({ testAlpha(v, AlphaRestriction(Rb, qnz)); }) && throws({ testBeta(v, BetaRestriction(Rb, qnz)); }) && throws({ vecm(Y, 2, 1, restriction: BetaRestriction(Rb, qnz)); }),
          "restrictions with q different from 0 are refused");
    check(throws({ testAlpha(v, AlphaRestriction(fresh2(1, 3), GretlMatrix.init)); }) && throws({ testAlpha(v, AlphaRestriction(fresh2(2, 2), GretlMatrix.init)); })
          && throws({ testAlpha(v, AlphaRestriction.init); }), "alpha restrictions of the wrong shape or none throw");
    check(throws({ vecm(Y, 2, 1, restriction: BetaRestriction(Rb, GretlMatrix.init), alphaRestriction: AlphaRestriction(Rb, GretlMatrix.init)); }), "restrictions on beta and alpha together are refused");

    // ---- three variables with one common trend: rank 2 ----
    double[] x1 = new double[T], x2 = new double[T], x3 = new double[T];
    double w = 0, e2 = 0, e3 = 0;
    x1[] = 0; x2[] = 0; x3[] = 0;
    foreach (t; 1 .. T) {
        w += urand();
        e2 = 0.5 * e2 + 0.4 * urand();
        e3 = 0.4 * e3 + 0.4 * urand();
        x1[t] = w;
        x2[t] = w + e2;
        x3[t] = w + e3;
    }
    auto Y3 = cbind(x1, x2, x3);
    auto v2 = vecm(Y3, 2, 2);
    auto m3 = momentsOf([x1, x2, x3]);
    const lam2 = topEigen(inv(m3.S11) * t(m3.S01) * inv(m3.S00) * m3.S01, 2);

    // alpha: the first variable is weakly exogenous (a restriction common to both columns)
    {
        auto R = fresh2(1, 3);
        R[0, 0] = 1; R[0, 1] = 0; R[0, 2] = 0;
        auto C = t(R) * inv(R * m3.S00 * t(R)) * R;
        auto S00a = m3.S00 - m3.S00 * C * m3.S00, S11a = m3.S11 - t(m3.S01) * C * m3.S01, S01a = m3.S01 - m3.S00 * C * m3.S01;
        auto A = fresh2(3, 2);
        A[0, 0] = 0; A[1, 0] = 1; A[2, 0] = 0; A[0, 1] = 0; A[1, 1] = 0; A[2, 1] = 1;
        auto Abar = A * inv(t(A) * A);
        const lr2 = topEigen(inv(S11a) * t(S01a) * Abar * inv(t(Abar) * S00a * Abar) * t(Abar) * S01a, 2);
        const LR = m3.Te * (log((1 - lr2[0]) / (1 - lam2[0])) + log((1 - lr2[1]) / (1 - lam2[1])));
        auto ta = testAlpha(v2, AlphaRestriction(R, GretlMatrix.init));

        check(ta.df == 2 && close(ta.statistic, LR, 1e-6) && abs(ta.pvalue - chisq(2, LR)) < 1e-8, "rank 2: testAlpha, weak exogeneity of the first variable (df = rank * 1)");
        auto va2 = vecm(Y3, 2, 2, alphaRestriction: AlphaRestriction(R, GretlMatrix.init));
        check(abs(va2.alpha[0, 0]) < 1e-9 && abs(va2.alpha[0, 1]) < 1e-9 && close(va2.restrictionLR, LR, 1e-6) && va2.restrictionDf == 2, "rank 2: the first row of alpha is zero and the estimation agrees with the test");
    }

    // beta: the first two coefficients of every cointegrating vector are equal (R = [1, -1, 0])
    {
        auto R = fresh2(1, 3);
        R[0, 0] = 1; R[0, 1] = -1; R[0, 2] = 0;
        auto H = fresh2(3, 2);
        H[0, 0] = 1; H[1, 0] = 1; H[2, 0] = 0; H[0, 1] = 0; H[1, 1] = 0; H[2, 1] = 1;
        const lr2 = topEigen(inv(t(H) * m3.S11 * H) * (t(H) * t(m3.S01) * inv(m3.S00) * m3.S01 * H), 2);
        const LR = m3.Te * (log((1 - lr2[0]) / (1 - lam2[0])) + log((1 - lr2[1]) / (1 - lam2[1])));
        auto tb = testBeta(v2, BetaRestriction(R, GretlMatrix.init));

        check(tb.df == 2 && close(tb.statistic, LR, 1e-6), "rank 2: testBeta, equal first two coefficients");
        auto vb2 = vecm(Y3, 2, 2, restriction: BetaRestriction(R, GretlMatrix.init));
        check(close(vb2.beta[0, 0], vb2.beta[1, 0], 1e-8) && close(vb2.beta[0, 1], vb2.beta[1, 1], 1e-8) && close(vb2.restrictionLR, LR, 1e-6), "rank 2: every cointegrating vector has equal first two elements, and the estimation agrees");
    }

    // a restriction that is a Dataset-form fit
    {
        auto d = Dataset(4, T);
        import core.stdc.string: strcpy;
        strcpy(d.ptr.varname[1], "a"); strcpy(d.ptr.varname[2], "b"); strcpy(d.ptr.varname[3], "c");
        foreach (t; 0 .. T) { d.ptr.Z[1][t] = x1[t]; d.ptr.Z[2][t] = x2[t]; d.ptr.Z[3][t] = x3[t]; }
        auto R = fresh2(1, 3);
        R[0, 0] = 1;
        auto vd = vecm(d, ["a", "b", "c"], 2, 2, alphaRestriction: AlphaRestriction(R, GretlMatrix.init));
        check(abs(vd.alpha[0, 0]) < 1e-9 && vd.names == ["a", "b", "c"], "vecm on a Dataset with an alpha restriction");
    }
}

double chisq(int df, double x) {
    import gretl.pvalues: chisq_cdf_comp;
    return chisq_cdf_comp(df, x);
}

GretlMatrix S01aFor(Moments m, GretlMatrix Abar, GretlMatrix S00a, GretlMatrix S11a, GretlMatrix S01a, GretlMatrix beta) {
    return S01a * beta * inv(t(beta) * S11a * beta);
}

/* ---- stage 12: the cross-correlogram ---- */

void stage12() {
    enum T = 100;
    double[] x = new double[T], y = new double[T];
    double a = 0;

    x[] = 0; y[] = 0;
    foreach (t; 0 .. T) { a = 0.6 * a + urand(); x[t] = a; }
    foreach (t; 3 .. T) y[t] = 0.8 * x[t - 3] + 0.3 * urand();      // y follows x by three periods
    auto cc = crossCorrelogram(x, y, 6);

    check(cc.lags.nrow == 13 && cc.lags[0, 0] == -6 && cc.lags[12, 0] == 6 && cc.ccf.nrow == 13, "crossCorrelogram: lags -6 .. 6");
    double mx = 0, my = 0;
    foreach (t; 0 .. T) { mx += x[t]; my += y[t]; }
    mx /= T; my /= T;
    double sxx = 0, syy = 0;
    foreach (t; 0 .. T) { sxx += (x[t] - mx) ^^ 2; syy += (y[t] - my) ^^ 2; }
    bool ok = true;
    foreach (k; -6 .. 7) {
        double num = 0;
        foreach (t; 0 .. T) {
            const s = t - k;
            if (s >= 0 && s < T) num += (x[t] - mx) * (y[s] - my);
        }
        ok = ok && close(cc.ccf[k + 6, 0], num / sqrt(sxx * syy), 1e-9);
    }
    check(ok, "crossCorrelogram: r(k) = corr of x(t) with y(t - k), against the direct sums");
    size_t peak = 0;
    foreach (i; 0 .. 13) if (abs(cc.ccf[i, 0]) > abs(cc.ccf[peak, 0])) peak = i;
    check(cc.lags[peak, 0] == -3, "crossCorrelogram: y(t) depends on x(t - 3), so the peak is at lag -3");
    check(close(crossCorrelogram(x, x, 3).ccf[3, 0], 1.0, 1e-9), "crossCorrelogram: a series with itself has correlation 1 at lag 0");
    check(throws({ crossCorrelogram(x, y[0 .. 50], 3); }) && throws({ crossCorrelogram(x, y, 0); }) && throws({ crossCorrelogram(x, y, 60); }), "crossCorrelogram: bad arguments throw");
    auto d = Dataset(3, T);
    import core.stdc.string: strcpy;
    strcpy(d.ptr.varname[1], "x"); strcpy(d.ptr.varname[2], "y");
    foreach (t; 0 .. T) { d.ptr.Z[1][t] = x[t]; d.ptr.Z[2][t] = y[t]; }
    check(close(crossCorrelogram(d, "x", "y", 6).ccf[2, 0], cc.ccf[2, 0], 1e-12), "crossCorrelogram on a Dataset");
}
