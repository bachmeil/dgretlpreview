/* Tests for dgretl.estimate: ols on matrices and arrays, checked against the formulas written out with the matrix operators
   (solve, inv, t, the element-wise forms), for the classical, HC0-HC3, HAC and cluster-robust covariance matrices.

   Build and run with testing/run_tests.sh. */

import core.stdc.stdio : printf;

import core.stdc.string : strcpy;
import std.exception : collectException;
import std.math : abs, exp, fabs, log, sqrt;
import std.mathspecial : erf;

import dgretl;

struct MallInfo {
    int arena, ordblks, smblks, hblks, hblkhd, usmblks, fsmblks, uordblks, fordblks, keepcost;
}
extern (C) MallInfo mallinfo();

int n_tests, n_failed;

void check(bool ok, string what)
{
    n_tests++;
    if (!ok) {
        n_failed++;
        printf("FAIL: %.*s\n", cast(int) what.length, what.ptr);
    } else {
        printf("ok:   %.*s\n", cast(int) what.length, what.ptr);
    }
}

bool close(double a, double b, double tol = 1e-8)
{
    return fabs(a - b) <= tol * (1 + fabs(b));
}

ulong lcg = 12345;

double urand()
{
    lcg = lcg * 6364136223846793005UL + 1442695040888963407UL;
    return ((lcg >> 11) / cast(double) (1UL << 53)) - 0.5;
}

GretlMatrix mat(size_t r, size_t c)
{
    return GretlMatrix(gretl_matrix_alloc(cast(int) r, cast(int) c));
}

bool vclose(GretlMatrix a, GretlMatrix b, double tol = 1e-8)
{
    if (a.nrow != b.nrow || a.ncol != b.ncol) return false;
    foreach (j; 0 .. a.ncol) foreach (i; 0 .. a.nrow) if (!close(a[i, j], b[i, j], tol)) return false;
    return true;
}

int main()
{
    enum n = 80;
    auto X = mat(n, 2), y = mat(n, 1);
    foreach (i; 0 .. n) {
        X[i, 0] = urand();
        X[i, 1] = urand() + 0.5 * X[i, 0];
        y[i, 0] = 1.0 + 2.0 * X[i, 0] - 1.5 * X[i, 1] + 0.5 * urand() * (1 + 3 * fabs(X[i, 0]));
    }
    auto Z = cbind(1, X);                                  // with the constant, for the by-hand formulas
    auto beta = solve(Z, y);
    auto u = y - Z * beta;
    auto XtXi = inv(t(Z) * Z);
    const k = 3;

    /* classical */
    auto f = ols(y, X);
    check(f.method == "OLS" && f.nobs == n && f.ncoef == 3 && f.df == n - 3 && f.names == ["const", "x1", "x2"], "ols: method, sizes and names");
    check(vclose(f.coef, beta), "ols: coefficients agree with solve(Z, y)");
    double ess = 0;
    foreach (i; 0 .. n) ess += u[i, 0] * u[i, 0];
    const s2 = ess / (n - k);
    auto V = XtXi * s2;
    check(close(f.ess, ess) && close(f.sigma, sqrt(s2)) && vclose(f.vcov, V, 1e-7), "ols: ESS, sigma and the classical covariance matrix");
    bool seOk = true;
    foreach (i; 0 .. k) seOk = seOk && close(f.se[i, 0], sqrt(V[i, i])) && close(f.tstat[i, 0], f.coef[i, 0] / f.se[i, 0]);
    check(seOk && f.pvalue[1, 0] > 0 && f.pvalue[1, 0] < 1e-6, "ols: standard errors, t ratios and p-values");
    check(f.resid.nrow == n && f.fitted.nrow == n && close(f.resid[5, 0], u[5, 0]) && close(f.fitted[5, 0] + f.resid[5, 0], y[5, 0]), "ols: residuals and fitted values");
    double ybar = 0;
    foreach (i; 0 .. n) ybar += y[i, 0] / n;
    double tss = 0;
    foreach (i; 0 .. n) tss += (y[i, 0] - ybar) ^^ 2;
    check(close(f.tss, tss) && close(f.rsq, 1 - ess / tss) && close(f.adjrsq, 1 - (ess / (n - k)) / (tss / (n - 1))) && f.fstat > 0, "ols: TSS, R-squared, adjusted R-squared, F");
    check(f.aic == f.aic && f.bic > f.aic - 100 && f.lnL == f.lnL, "ols: information criteria and log-likelihood are numbers");

    /* the other ways in: arrays, views, no constant */
    auto x1 = new double[n], yy = new double[n];
    foreach (i; 0 .. n) { x1[i] = X[i, 0]; yy[i] = y[i, 0]; }
    auto fa = ols(yy, x1);
    auto b1 = solve(cbind(1, x1), y);
    check(fa.ncoef == 2 && close(fa.coef[0, 0], b1[0, 0]) && close(fa.coef[1, 0], b1[1, 0]) && fa.names == ["const", "x1"], "ols with arrays");
    check(vclose(ols(y.cols[0 .. 1], X.sub[0 .. n, 0 .. 2]).coef, beta), "ols with views");
    auto fnc = ols(y, X, constant: false);
    check(fnc.ncoef == 2 && fnc.names == ["x1", "x2"] && vclose(fnc.coef, solve(X, y)), "ols without a constant (named argument)");
    auto fz = ols(y, Z, constant: false);
    check(vclose(fz.coef, beta) && fz.ncoef == 3, "ols with an explicit column of ones");

    /* HC0-HC3 */
    auto h = mat(n, 1);
    foreach (i; 0 .. n) {
        double hi = 0;
        foreach (a; 0 .. k) foreach (b; 0 .. k) hi += Z[i, a] * XtXi[a, b] * Z[i, b];
        h[i, 0] = hi;
    }
    GretlMatrix sandwich(double delegate(size_t) wt)
    {
        auto meat = mat(k, k);
        foreach (a; 0 .. k) foreach (b; 0 .. k) {
            double acc = 0;
            foreach (i; 0 .. n) acc += Z[i, a] * Z[i, b] * u[i, 0] * u[i, 0] * wt(i);
            meat[a, b] = acc;
        }
        return XtXi * meat * XtXi;
    }
    auto hc0 = sandwich((i) => 1.0);
    auto hc1 = sandwich((i) => cast(double) n / (n - k));
    auto hc2 = sandwich((i) => 1.0 / (1.0 - h[i, 0]));
    auto hc3 = sandwich((i) => 1.0 / ((1.0 - h[i, 0]) ^^ 2));
    check(vclose(ols(y, X, vcov: Vcov.hc0).vcov, hc0, 1e-6), "ols: HC0 covariance matrix");
    check(vclose(ols(y, X, vcov: Vcov.hc1).vcov, hc1, 1e-6) && ols(y, X, vcov: Vcov.hc1).vcovType == "HC1", "ols: HC1 covariance matrix");
    check(vclose(ols(y, X, vcov: Vcov.hc2).vcov, hc2, 1e-6), "ols: HC2 covariance matrix");
    check(vclose(ols(y, X, vcov: Vcov.hc3).vcov, hc3, 1e-6), "ols: HC3 covariance matrix");
    auto fr = ols(y, X, vcov: Vcov.hc1);
    check(vclose(fr.coef, beta) && fr.se[1, 0] != f.se[1, 0], "robust standard errors differ from the classical ones, the coefficients do not");
    check(ols(y, X).se[1, 0] == f.se[1, 0], "the library settings are restored after a robust call (a classical call is unchanged)");

    /* clusters */
    auto grp = new double[n];
    foreach (i; 0 .. n) grp[i] = i % 8;
    auto fc = ols(y, X, cluster: grp);
    check(fc.vcovType == "cluster-robust" && fc.df == 7 && vclose(fc.coef, beta), "ols: cluster-robust call runs, df is the number of clusters minus one");
    auto meatC = mat(k, k);
    foreach (a; 0 .. k) foreach (b; 0 .. k) {
        double acc = 0;
        foreach (g; 0 .. 8) {
            double sa = 0, sb = 0;
            foreach (i; 0 .. n) if (i % 8 == g) { sa += Z[i, a] * u[i, 0]; sb += Z[i, b] * u[i, 0]; }
            acc += sa * sb;
        }
        meatC[a, b] = acc;
    }
    auto vcl = XtXi * meatC * XtXi;
    const ratio = fc.vcov[1, 1] / vcl[1, 1];
    printf("  cluster vcov / by-hand (no small-sample factor) = %g\n", ratio);
    check(close(ratio, 8.0 / 7 * (n - 1.0) / (n - k), 1e-6), "ols: cluster-robust covariance is the by-hand one times G/(G-1) * (n-1)/(n-k)");
    check(collectException!GretlException(ols(y, X, cluster: grp[0 .. 10])) !is null, "a cluster series of the wrong length throws");

    /* HAC */
    auto fh = ols(y, X, vcov: Vcov.hac, lag: 3);
    auto meatH = mat(k, k);
    foreach (a; 0 .. k) foreach (b; 0 .. k) {
        double acc = 0;
        foreach (l; 0 .. 4) {
            double gam = 0;
            foreach (i; l .. n) gam += Z[i, a] * u[i, 0] * Z[i - l, b] * u[i - l, 0];
            double gam2 = 0;
            foreach (i; l .. n) gam2 += Z[i, b] * u[i, 0] * Z[i - l, a] * u[i - l, 0];
            const w = 1.0 - l / 4.0;
            acc += (l == 0) ? gam : w * (gam + gam2);
        }
        meatH[a, b] = acc;
    }
    auto vhac = XtXi * meatH * XtXi;
    const hr = fh.vcov[1, 1] / vhac[1, 1];
    printf("  HAC vcov / by-hand Bartlett (lag 3) = %g, %g\n", hr, fh.vcov[0, 0] / vhac[0, 0]);
    check(fh.vcovType == "HAC" && vclose(fh.coef, beta) && close(hr, 1.0, 1e-6) && vclose(fh.vcov, vhac, 1e-6), "ols: HAC covariance matrix with a fixed lag is the Bartlett estimator, no degrees-of-freedom factor");
    check(ols(y, X, vcov: Vcov.hac).se[1, 0] > 0 && ols(y, X, vcov: Vcov.hac, kernel: HACKernel.parzen, lag: 4).se[1, 0] > 0 && ols(y, X, vcov: Vcov.hac, prewhiten: true).se[1, 0] > 0, "ols: HAC with the default lag, Parzen, and pre-whitening");

    /* errors and edge cases */
    check(collectException!GretlException(ols(y.rows[0 .. 10], X)) !is null, "ols with different numbers of rows throws");
    check(collectException!GretlException(ols(cbind(y, y), X)) !is null, "ols with a two-column y throws");
    auto bad = y.dup;
    bad[3, 0] = double.nan;
    check(collectException!GretlException(ols(bad, X)) !is null, "ols with a missing value throws");
    check(collectException!GretlException(ols(y.rows[0 .. 3], X.rows[0 .. 3])) !is null, "ols with too few observations throws");
    auto Xc = cbind(X, X.cols[0 .. 1]);                    // a duplicated column
    auto fdup = ols(y, Xc);
    check(fdup.ncoef == 3 && fdup.names.length == 3 && vclose(fdup.coef, beta), "a perfectly collinear regressor is dropped");
    auto text = f.toString();
    check(text.length > 50 && text[0 .. 3] == "OLS", "toString gives a table");

    /* WLS: the weighted problem by hand, with weights w and w^(1/2) scaling of the rows of [1 X] and y */
    auto wts = new double[n];
    foreach (i; 0 .. n) wts[i] = 0.5 + 2.0 * (urand() + 0.5);
    auto fw = wls(y, X, wts);
    auto Zs = mat(n, k), ys = mat(n, 1);
    foreach (i; 0 .. n) {
        foreach (j; 0 .. k) Zs[i, j] = Z[i, j] * sqrt(wts[i]);
        ys[i, 0] = y[i, 0] * sqrt(wts[i]);
    }
    auto bw = solve(Zs, ys);
    auto uw = ys - Zs * bw;
    double essw = 0;
    foreach (i; 0 .. n) essw += uw[i, 0] ^^ 2;
    auto Vw = inv(t(Zs) * Zs) * (essw / (n - k));
    check(fw.method == "WLS" && vclose(fw.coef, bw) && vclose(fw.vcov, Vw, 1e-6) && fw.nobs == n, "wls: coefficients and covariance matrix from the weighted problem");
    check(!vclose(fw.coef, beta, 1e-12) && close(fw.ess, essw, 1e-7), "wls: differs from ols, and ESS is the weighted one");
    auto ones = new double[n];
    ones[] = 1.0;
    check(vclose(wls(y, X, ones).coef, beta) && vclose(wls(y, X, ones).vcov, V, 1e-7), "wls with equal weights is ols");
    auto dummyw = ones.dup;
    dummyw[0] = 0;
    check(wls(y, X, dummyw).nobs == n - 1 && vclose(wls(y, X, dummyw).coef, ols(y.rows[1 .. n], X.rows[1 .. n]).coef), "a 0/1 weight series leaves out the observations with weight 0");
    auto zerow = wts.dup;
    zerow[0] = 0;
    check(collectException!GretlException(wls(y, X, zerow)) !is null, "a zero weight among other weights throws");
    auto negw = wts.dup;
    negw[2] = -1;
    check(collectException!GretlException(wls(y, X, negw)) !is null && collectException!GretlException(wls(y, X, wts[0 .. 5])) !is null, "wls: negative weights and weights of the wrong length throw");
    check(vclose(wls(y, X, wts, vcov: Vcov.hc1).coef, bw) && wls(y, X, wts, vcov: Vcov.hc1).vcovType == "HC1", "wls with robust standard errors");

    /* 2SLS: x1 is endogenous (correlated with the error), z1 and z2 are instruments, w is an exogenous regressor */
    auto zi = mat(n, 2), xe = mat(n, 1), we = mat(n, 1), yi = mat(n, 1);
    foreach (i; 0 .. n) {
        const e1 = urand(), e2 = urand();
        zi[i, 0] = urand();
        zi[i, 1] = urand();
        we[i, 0] = urand();
        xe[i, 0] = zi[i, 0] + 0.7 * zi[i, 1] + 0.5 * we[i, 0] + 0.8 * e1 + 0.3 * e2;
        yi[i, 0] = 1.0 + 2.0 * xe[i, 0] - 1.0 * we[i, 0] + e1;
    }
    auto fi = tsls(yi, xe, zi, exog: we);
    auto Xr = cbind(1, we, xe);                                      // regressors in the order const, exog, endog
    auto Zr = cbind(1, we, zi);                                      // instruments
    auto Pz = Zr * inv(t(Zr) * Zr) * t(Zr);
    auto bi = inv(t(Xr) * Pz * Xr) * (t(Xr) * Pz * yi);
    auto ui = yi - Xr * bi;
    double essi = 0;
    foreach (i; 0 .. n) essi += ui[i, 0] ^^ 2;
    check(fi.method == "2SLS" && fi.names == ["const", "w1", "x1"] && vclose(fi.coef, bi, 1e-6), "tsls: coefficients agree with (X'PX)^-1 X'P y");
    auto Vi = inv(t(Xr) * Pz * Xr);
    const ratioI = fi.vcov[2, 2] / Vi[2, 2];
    printf("  2SLS sigma^2 = %g (by hand ess / n = %g, ess / (n - k) = %g)\n", ratioI, essi / n, essi / (n - 3));
    check(close(ratioI, essi / (n - 3), 1e-6), "tsls: covariance matrix is ess / (n - k) times (X'PX)^-1");
    check(fi.se[2, 0] > 0 && fi.pvalue[2, 0] < 1e-6 && fi.resid.nrow == n, "tsls: standard errors, p-values, residuals");
    auto fo = ols(yi, Xr.cols[1 .. 3]);                              // OLS is biased here: the endogenous coefficient differs from 2SLS
    check(abs(fo.coef[2, 0] - 2.0) > abs(fi.coef[2, 0] - 2.0), "tsls is closer to the truth than ols when x is endogenous");
    check(tsls(yi, xe, zi, exog: we, vcov: Vcov.hc1).vcovType == "HC1" && tsls(yi, xe, zi).ncoef == 2 && tsls(yi, xe, zi, constant: false).ncoef == 1, "tsls: robust option, no exog, no constant");
    check(collectException!GretlException(tsls(yi, cbind(xe, we), zi.cols[0 .. 1])) !is null, "tsls: too few instruments throws");
    check(collectException!GretlException(tsls(yi.rows[0 .. 9], xe, zi)) !is null, "tsls: sizes that do not fit throw");

    /* logit and probit */
    auto yb = mat(n, 1);
    foreach (i; 0 .. n) {
        const uu = urand() + 0.5;
        yb[i, 0] = (0.3 + 1.2 * X[i, 0] - 0.8 * X[i, 1] + log((uu + 1e-9) / (1 - uu + 1e-9)) > 0) ? 1.0 : 0.0;
    }
    auto fl = logit(yb, X);
    auto pr = mat(n, 1);
    auto Wd = mat(n, 1);
    foreach (i; 0 .. n) {
        double eta = 0;
        foreach (j; 0 .. k) eta += Z[i, j] * fl.coef[j, 0];
        pr[i, 0] = 1.0 / (1.0 + exp(-eta));
    }
    auto grad = t(Z) * (yb - pr);
    double lnl = 0;
    foreach (i; 0 .. n) lnl += yb[i, 0] * log(pr[i, 0]) + (1 - yb[i, 0]) * log(1 - pr[i, 0]);
    auto Zw = mat(n, k);
    foreach (i; 0 .. n) foreach (j; 0 .. k) Zw[i, j] = Z[i, j] * sqrt(pr[i, 0] * (1 - pr[i, 0]));
    auto Vl = inv(t(Zw) * Zw);
    check(fl.method == "Logit" && fl.nobs == n && fl.ncoef == 3 && fl.names == ["const", "x1", "x2"], "logit: method, sizes and names");
    check(norm(grad, Norm.max) < 1e-5 && close(fl.lnL, lnl, 1e-8), "logit: the score is zero at the estimates and lnL agrees");
    check(vclose(fl.vcov, Vl, 1e-4) && fl.df == 0 && fl.pvalue[1, 0] > 0 && fl.pvalue[1, 0] < 1, "logit: covariance matrix is (X'WX)^-1; z-ratio p-values");
    check(fl.aic == fl.aic && fl.fitted.nrow == n && close(fl.fitted[3, 0], pr[3, 0], 1e-6), "logit: fitted values are the probabilities");
    auto fp = probit(yb, X);
    double score0 = 0, score1 = 0;
    foreach (i; 0 .. n) {
        double eta = 0;
        foreach (j; 0 .. k) eta += Z[i, j] * fp.coef[j, 0];
        const phi = exp(-0.5 * eta * eta) / sqrt(2 * 3.141592653589793), cdf = 0.5 * (1 + erf(eta / sqrt(2.0)));
        score0 += phi * (yb[i, 0] - cdf) / (cdf * (1 - cdf));
        score1 += phi * (yb[i, 0] - cdf) / (cdf * (1 - cdf)) * Z[i, 1];
    }
    check(fp.method == "Probit" && abs(score0) < 1e-4 && abs(score1) < 1e-4 && fp.lnL < 0 && fp.se[1, 0] > 0, "probit: the score is zero at the estimates");
    check(abs(fp.coef[1, 0] / fl.coef[1, 0] - 0.55) < 0.15, "probit coefficients are about 0.55 times the logit ones");
    check(logit(yb, X, robust: true).vcovType == "QML (sandwich)" && vclose(logit(yb, X, robust: true).coef, fl.coef, 1e-6), "logit with robust standard errors");
    check(collectException!GretlException(logit(y, X)) !is null && collectException!GretlException(probit(yb.rows[0 .. 10], X)) !is null, "logit: y that is not 0/1 and sizes that do not fit throw");
    check(logit(yb, X, constant: false).ncoef == 2, "logit without a constant");

    /* Dataset-based wrappers */
    auto ds = Dataset(10, n);                       // the constant plus nine series
    {
        auto raw = ds.ptr;
        foreach (v, nm; ["const", "gdp", "unemp", "rate", "grp", "wt", "yb", "xe", "ye", "zin"]) strcpy(raw.varname[v], (nm ~ "\0").ptr);
        foreach (i; 0 .. n) {
            raw.Z[1][i] = y[i, 0];
            raw.Z[2][i] = X[i, 0];
            raw.Z[3][i] = X[i, 1];
            raw.Z[4][i] = grp[i];
            raw.Z[5][i] = wts[i];
            raw.Z[6][i] = yb[i, 0];
            raw.Z[7][i] = xe[i, 0];
            raw.Z[8][i] = yi[i, 0];
            raw.Z[9][i] = zi[i, 0];
        }
    }
    auto fd = ols(ds, "gdp", ["unemp", "rate"]);
    check(fd.names == ["const", "unemp", "rate"] && vclose(fd.coef, beta) && vclose(fd.vcov, V, 1e-7) && fd.nobs == n, "ols on a Dataset: names are the series names, results match the matrix form");
    check(vclose(ols(ds, "gdp", ["unemp", "rate"], vcov: Vcov.hc1).vcov, hc1, 1e-6) && ols(ds, "gdp", ["rate"], constant: false).names == ["rate"], "ols on a Dataset with options (vcov, constant)");
    check(vclose(ols(ds, "gdp", ["unemp", "rate"], cluster: "grp").vcov, fc.vcov, 1e-8), "ols on a Dataset with a cluster series given by name");
    check(vclose(wls(ds, "gdp", ["unemp", "rate"], "wt").coef, bw) && wls(ds, "gdp", ["unemp", "rate"], "wt").names == ["const", "unemp", "rate"], "wls on a Dataset");
    auto fdi = tsls(ds, "ye", ["xe"], ["zin"]);
    auto fdm = tsls(yi, xe, zi.cols[0 .. 1]);
    check(fdi.names == ["const", "xe"] && vclose(fdi.coef, fdm.coef, 1e-8), "tsls on a Dataset");
    auto fdl = logit(ds, "yb", ["unemp", "rate"]), fdp = probit(ds, "yb", ["unemp"]);
    check(fdl.names == ["const", "unemp", "rate"] && vclose(fdl.coef, fl.coef, 1e-8) && fdp.method == "Probit" && fdp.names == ["const", "unemp"], "logit and probit on a Dataset");
    check(collectException!GretlException(ols(ds, "gdp", ["nosuch"])) !is null && collectException!GretlException(ols(ds, "nosuch", ["rate"])) !is null, "a series that is not there throws");

    /* a Dataset with gaps: observations with a missing value in a series used are left out */
    {
        auto raw = ds.ptr;
        raw.Z[2][4] = double.nan;
        raw.Z[1][10] = double.nan;
    }
    auto fg = ols(ds, "gdp", ["unemp", "rate"]);
    auto keepRows = mat(n - 2, 3), keepY = mat(n - 2, 1);
    size_t rr = 0;
    foreach (i; 0 .. n) {
        if (i == 4 || i == 10) continue;
        keepRows[rr, 0] = 1; keepRows[rr, 1] = X[i, 0]; keepRows[rr, 2] = X[i, 1]; keepY[rr, 0] = y[i, 0];
        rr++;
    }
    check(fg.nobs == n - 2 && fg.resid.nrow == n - 2 && vclose(fg.coef, solve(keepRows, keepY), 1e-8), "a Dataset with missing values: those observations are left out");
    check(ols(ds, "rate", ["grp"]).nobs == n, "series that have no gaps keep every observation");

    /* no leak: after a warm-up, many estimations leave the heap as it was */
    /* (a long warm-up: the allocator of the Windows C library takes one-time growth of the heap in the first rounds) */
    foreach (i; 0 .. 250) { auto w1 = ols(y, X); auto w2 = logit(yb, X); auto w3 = tsls(yi, xe, zi, exog: we); auto w4 = wls(y, X, wts, vcov: Vcov.hac); }
    int before = mallinfo().uordblks;
    foreach (i; 0 .. 200) { auto w1 = ols(y, X, vcov: Vcov.hc1); auto w2 = logit(yb, X); auto w3 = tsls(yi, xe, zi, exog: we); auto w4 = wls(y, X, wts, vcov: Vcov.hac); }
    int after = mallinfo().uordblks;
    printf("  heap in use: %d bytes before, %d after 800 estimations\n", before, after);
    check(after - before < 64 * 1024, "estimating 800 times does not leak");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
