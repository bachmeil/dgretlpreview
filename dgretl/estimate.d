/* Estimation: regressions on matrices (or arrays), as a D-friendly layer over the estimators of the port. An estimator takes the data as
   GretlMatrix, Submatrix or arrays, runs the library's code (gretl.estimate) on a temporary dataset made from them, and hands back an
   Estimate: a plain D struct holding GretlMatrix fields and numbers, so nothing has to be cleared or freed by hand.

   import dgretl;
   auto fit = ols(y, X);                                       // a constant is added; classical standard errors
   auto fit = ols(y, X, constant: false, vcov: Vcov.hc1);      // robust (heteroskedasticity-consistent) standard errors
   writeln(fit);                                               // a table of the coefficients
   fit.coef, fit.se, fit.vcov, fit.resid, fit.rsq, ...

   Every estimator throws GretlException if the library reports an error. The library settings the options change (HC version, HAC kernel
   and lag, pre-whitening) are restored afterwards. */

module dgretl.estimate;

public import dgretl.matrix;
public import dgretl.exception;
public import gretl.gretl_modelh: HACKernel;

import core.stdc.string: strcpy;
import std.conv: to;
import std.exception: enforce;
import std.format: format;
import std.math: fabs, isFinite, isNaN;
import std.string: fromStringz, toStringz;

import cephes.ndtr: ndtr;
import cephes.stdtr: stdtr;
import dgretl.dataset: Dataset;
import dgretl.internal;
import gretl.discrete: binary_logit, binary_probit;
import gretl.estimate: ivreg, lsq;
import gretl.gretl_commandsh: GretlCmdIndex;
import gretl.gretl_errors: gretl_errmsg_get;
import gretl.gretl_list: gretl_list_new;
import gretl.gretl_listh: LISTSEP;
import gretl.gretl_model: clear_model;
import gretl.libgretlh: DatasetRaw, Gretlopt, InfoCriterion, Model, TsCodes;
import gretl.libset: execute_set, pop_program_state, push_program_state;
import gretl.matrix_extra: gretl_vcv_matrix_from_model;
import gretl.options: set_optval_string;

/* The covariance matrix of the coefficients. */
enum Vcov {
    classical,      /* sigma^2 (X'X)^-1 */
    hc0,            /* White's heteroskedasticity-consistent estimator, no correction */
    hc1,            /* HC0 times n / (n - k): the usual "robust" standard errors */
    hc2,            /* squared residuals divided by 1 - h (the leverage) */
    hc3,            /* squared residuals divided by (1 - h)^2: the jackknife-like version */
    hac             /* heteroskedasticity and autocorrelation consistent (Newey-West and relatives), the data in time order */
}

/* The result of an estimation. Vectors of coefficient statistics are k x 1 GretlMatrices in the order of `names`; resid and fitted are n x 1. */
struct Estimate {
    string method;              // "OLS", "WLS", "2SLS", "Logit", "Probit"
    string vcovType;            // what the standard errors are: "classical", "HC1", "HAC", ...
    string[] names;             // coefficient names: "const", then the names given or "x1", "x2", ...
    GretlMatrix coef;           // k x 1 coefficients
    GretlMatrix se;             // k x 1 standard errors
    GretlMatrix tstat;          // k x 1 coef / se (a z-ratio for the maximum likelihood estimators)
    GretlMatrix pvalue;         // k x 1 two-sided p-values (Student t with `df` degrees of freedom, or the normal)
    GretlMatrix vcov;           // k x k covariance matrix of the coefficients
    GretlMatrix resid;          // n x 1 residuals (y - fitted for the least squares estimators)
    GretlMatrix fitted;         // n x 1 fitted values
    size_t nobs;                // number of observations used
    size_t ncoef;               // number of coefficients (a perfectly collinear regressor is dropped by the library)
    int df;                     // residual degrees of freedom, nobs - ncoef for least squares
    double sigma = double.nan;  // standard error of the regression
    double ess = double.nan;    // error (residual) sum of squares
    double tss = double.nan;    // total sum of squares
    double rsq = double.nan;    // R squared
    double adjrsq = double.nan; // adjusted R squared
    double fstat = double.nan;  // overall F statistic
    double lnL = double.nan;    // log-likelihood
    double aic = double.nan, bic = double.nan, hqc = double.nan;    // information criteria (gretl's definitions)
    double dw = double.nan;     // Durbin-Watson statistic
    double rho = double.nan;    // first-order autocorrelation of the residuals
    GretlMatrix arcoef;         // autoregressive errors (ar1, ar): the AR coefficients of the errors (p x 1), otherwise empty
    GretlMatrix arse;           // their standard errors
    int[] arlags;               // the lags they belong to

    /* a table of the coefficients, with the usual summary lines */
    string toString() {
        string s = format("%s, n = %d, df = %d, standard errors: %s\n", method, nobs, df, vcovType);

        s ~= format("%-12s %12s %12s %10s %9s\n", "", "coef", "std.err", "t-ratio", "p-value");
        foreach (i, nm; names) {
            s ~= format("%-12s %12.6g %12.6g %10.4f %9.4f\n", nm, coef[i, 0], se[i, 0], tstat[i, 0], pvalue[i, 0]);
        }
        if (!isNaN(rsq)) s ~= format("R-squared %.6f, adjusted %.6f, sigma %.6g, F %.4g\n", rsq, adjrsq, sigma, fstat);
        if (!isNaN(lnL)) s ~= format("log-likelihood %.6f, AIC %.4f, BIC %.4f\n", lnL, aic, bic);
        return s;
    }
}

/* the data as an n x k GretlMatrix (a copy), whether it comes as a matrix, a view or an array (one column) */
private GretlMatrix asMatrix(T)(T x) {
    return cbind(x);
}

/* copy what the user needs out of a library model; names[v] is the name of series v of the temporary dataset */
package(dgretl) Estimate fromModel(ref Model m, string method, string vcovType, string[] names, size_t n, bool needVcov = true) @trusted {
    Estimate e;
    const k = cast(size_t) m.ncoeff;

    e.method = method;
    e.vcovType = vcovType;
    e.nobs = m.nobs;
    e.ncoef = k;
    e.df = m.dfd;
    foreach (i; 0 .. k) {
        const id = (m.list !is null) ? m.list[i + 2] : cast(int) i;

        e.names ~= (id == 0) ? "const" : ((id > 0 && id < cast(int) names.length) ? names[id] : "x" ~ to!string(id));
    }
    e.coef = columnOf(m.coeff, 0, k);
    e.se = columnOf(m.sderr, 0, k);
    e.tstat = fresh(k, 1);
    e.pvalue = fresh(k, 1);
    foreach (i; 0 .. k) {
        const t = e.coef.ptr.val[i] / e.se.ptr.val[i];

        e.tstat.ptr.val[i] = t;
        e.pvalue.ptr.val[i] = (e.se.ptr.val[i] > 0 && isFinite(t)) ? 2.0 * stdtr(e.df, -fabs(t)) : double.nan;
    }
    if (needVcov || m.vcv !is null) {
        int err = 0;
        auto vraw = gretl_vcv_matrix_from_model(&m, null, &err);

        enforce!GretlException(err == 0 && vraw !is null, "could not get the covariance matrix of the estimates");
        e.vcov = GretlMatrix(vraw);
    }
    e.resid = columnOf(m.uhat, m.t1, n);
    e.fitted = columnOf(m.yhat, m.t1, n);
    e.sigma = m.sigma;
    e.ess = m.ess;
    e.tss = m.tss;
    e.rsq = m.rsq;
    e.adjrsq = m.adjrsq;
    e.fstat = m.fstt;
    e.lnL = m.lnL;
    e.aic = m.criterion[InfoCriterion.aic];
    e.bic = m.criterion[InfoCriterion.bic];
    e.hqc = m.criterion[InfoCriterion.hqc];
    e.dw = m.dw;
    e.rho = m.rho;
    if (m.arinfo !is null && m.arinfo.arlist !is null) {
        const p = cast(size_t) m.arinfo.arlist[0];

        e.arcoef = columnOf(m.arinfo.rho, 0, p);
        e.arse = columnOf(m.arinfo.sderr, 0, p);
        foreach (i; 0 .. p) e.arlags ~= m.arinfo.arlist[i + 1];
    }
    return e;
}

/* the options that choose the covariance matrix, shared by the estimators; `apply` sets the library for one estimation, inside a
   push_program_state / pop_program_state pair of the caller, and returns the option flags */
private struct Robust {
    Vcov vcov;
    int lag;
    HACKernel kernel;
    bool prewhiten;
    const(double)[] cluster;

    Gretlopt apply(DatasetRaw * raw, string clusterName, out string label) @trusted {
        Gretlopt opt = Gretlopt.OPT_NONE;

        label = "classical";
        final switch (vcov) {
            case Vcov.classical:
                break;
            case Vcov.hc0, Vcov.hc1, Vcov.hc2, Vcov.hc3:
                const v = cast(int) vcov - cast(int) Vcov.hc0;

                opt = Gretlopt.OPT_R;
                label = "HC" ~ to!string(v);
                setLibrary("hc_version", to!string(v));
                setLibrary("force_hc", "1");
                break;
            case Vcov.hac:
                opt = Gretlopt.OPT_R;
                label = "HAC";
                raw.structure = TsCodes.time_series;
                raw.pd = 1;
                setLibrary("force_hc", "0");
                setLibrary("hac_kernel", ["bartlett", "parzen", "qs"][kernel]);
                setLibrary("hac_prewhiten", prewhiten ? "1" : "0");
                if (lag > 0) setLibrary("hac_lag", to!string(lag));
                break;
        }
        if (cluster !is null) {
            opt |= Gretlopt.OPT_R | Gretlopt.OPT_C;
            label = "cluster-robust";
            set_optval_string(GretlCmdIndex.ols, Gretlopt.OPT_C, clusterName.toStringz);
            set_optval_string(GretlCmdIndex.ivreg, Gretlopt.OPT_C, clusterName.toStringz);
        }
        return opt;
    }
}

private void checkCluster(const(double)[] cluster, size_t n, string who) {
    enforce!GretlException(cluster is null || cluster.length == n, who ~ ": the cluster series needs one value per observation");
}

/* ols(y, x): least squares of y (n x 1) on the columns of x (n x k) and, unless constant is false, a constant.
   y and x are GretlMatrix or Submatrix, or arrays (an array is one column).
       vcov       the covariance matrix of the coefficients: Vcov.classical (default), hc0 .. hc3, or hac
       lag        for Vcov.hac, the number of lags; 0 (default) lets the library choose from the sample size
       kernel     for Vcov.hac: HACKernel.bartlett (default), parzen or qs
       prewhiten  for Vcov.hac: pre-whiten the data
       cluster    a series (one value per observation) that identifies clusters: the standard errors are made robust to correlation
                  within them (cluster-robust), whatever `vcov` is
   A perfectly collinear regressor is dropped by the library (the result then has fewer coefficients than columns, named in `names`).
   Missing or non-finite values throw GretlException. */
Estimate ols(Y, X)(Y y, X x, bool constant = true, Vcov vcov = Vcov.classical, int lag = 0, HACKernel kernel = HACKernel.bartlett,
                   bool prewhiten = false, const(double)[] cluster = null) @trusted if (isOperand!Y && isOperand!X) {
    auto ym = asMatrix(y), xm = asMatrix(x);
    const n = ym.nrow, k = xm.ncol;

    enforce!GretlException(ym.ncol == 1, "ols: y must have one column, it has " ~ to!string(ym.ncol));
    enforce!GretlException(xm.nrow == n, "ols: y has " ~ to!string(n) ~ " observations but x has " ~ to!string(xm.nrow) ~ " rows");
    checkCluster(cluster, n, "ols");
    enforce!GretlException(n > k + (constant ? 1 : 0), "ols: too few observations (" ~ to!string(n) ~ ") for " ~ to!string(k + (constant ? 1 : 0)) ~ " coefficients");
    requireFinite(ym, "ols: y");
    requireFinite(xm, "ols: x");

    auto names = ["const", "y"] ~ numbered("x", k);
    auto cols = columnsOf(ym) ~ columnsOf(xm);

    if (cluster !is null) {
        cols ~= cluster;
        names ~= "cluster";
    }
    auto ds = makeDataset(cols, names);
    int[] spec = [1];

    if (constant) spec ~= 0;
    foreach (j; 0 .. k) spec ~= cast(int) (2 + j);
    return runLsq(ds, spec, GretlCmdIndex.ols, "OLS", Robust(vcov, lag, kernel, prewhiten, cluster), names, n);
}

private Estimate runLsq(Dataset ds, int[] spec, GretlCmdIndex ci, string method, Robust rb, string[] names, size_t n) @trusted {
    auto raw = ds.ptr;
    auto list = makeList(spec);
    string vtype;

    enforce!GretlException(push_program_state() == 0, "could not save the library settings");
    scope (exit) pop_program_state();
    auto opt = rb.apply(raw, names[$ - 1], vtype);
    Model m;

    scope (exit) {
        clear_model(&m);
        import core.stdc.stdlib: free;
        free(list);
    }
    m = (ci == GretlCmdIndex.ivreg) ? ivreg(list, raw, opt) : lsq(list, raw, ci, opt);
    enforce!GretlException(m.errcode == 0, method ~ ": the estimation failed, " ~ errorText(m));
    return fromModel(m, method, vtype, names, n);
}

/* wls(y, x, w): weighted least squares of y on x (and a constant), where the weight series w (one positive value per observation) is the
   inverse of the error variance of the observation: the sum of w * residual^2 is minimized, so that an observation with weight 4 counts as four
   observations of weight 1. Weights must be positive, except that a series of only 0s and 1s is accepted as well and leaves out the observations with
   weight 0; any other zero weight, and a negative weight, throws GretlException. The other arguments are those of ols. The residuals and fitted values are the unweighted ones, y - x * coef. */
Estimate wls(Y, X)(Y y, X x, const(double)[] w, bool constant = true, Vcov vcov = Vcov.classical, int lag = 0, HACKernel kernel = HACKernel.bartlett,
                    bool prewhiten = false, const(double)[] cluster = null) @trusted if (isOperand!Y && isOperand!X) {
    auto ym = asMatrix(y), xm = asMatrix(x);
    const n = ym.nrow, k = xm.ncol;

    enforce!GretlException(ym.ncol == 1, "wls: y must have one column, it has " ~ to!string(ym.ncol));
    enforce!GretlException(xm.nrow == n, "wls: y has " ~ to!string(n) ~ " observations but x has " ~ to!string(xm.nrow) ~ " rows");
    enforce!GretlException(w.length == n, "wls: the weights need one value per observation");
    checkCluster(cluster, n, "wls");
    enforce!GretlException(n > k + (constant ? 1 : 0), "wls: too few observations");
    requireFinite(ym, "wls: y");
    requireFinite(xm, "wls: x");
    foreach (i, wi; w) {
        enforce!GretlException(isFinite(wi) && wi >= 0, "wls: the weight in position " ~ to!string(i) ~ " is negative or not finite");
    }
    auto names = ["const", "y"] ~ numbered("x", k) ~ ["weight"];
    auto cols = columnsOf(ym) ~ columnsOf(xm) ~ [w];

    if (cluster !is null) {
        cols ~= cluster;
        names ~= "cluster";
    }
    auto ds = makeDataset(cols, names);
    int[] spec = [cast(int) (2 + k), 1];            // the weight series comes first in a WLS list

    if (constant) spec ~= 0;
    foreach (j; 0 .. k) spec ~= cast(int) (2 + j);
    return runLsq(ds, spec, GretlCmdIndex.wls, "WLS", Robust(vcov, lag, kernel, prewhiten, cluster), names, n);
}

/* tsls(y, endog, instruments): two-stage least squares of y on the endogenous regressors `endog` (n x k1) and, if given, the included exogenous
   regressors `exog` (n x k2) and a constant, using `instruments` (n x m, the EXCLUDED instruments) and also the exogenous regressors and the constant
   as instruments. The model is identified if m >= k1; the library then reports the Sargan overidentification test in the model it builds (not copied here).
   Coefficients come in the order: constant, exog, endog, with the names "const", "w1", ..., "x1", ... Other arguments as for ols (the
   cluster-robust form is not available for 2SLS: it throws). */
Estimate tsls(Y, E, Z)(Y y, E endog, Z instruments, GretlMatrix exog = GretlMatrix.init, bool constant = true, Vcov vcov = Vcov.classical,
                       int lag = 0, HACKernel kernel = HACKernel.bartlett, bool prewhiten = false)
                       @trusted if (isOperand!Y && isOperand!E && isOperand!Z) {
    auto ym = asMatrix(y), em = asMatrix(endog), zm = asMatrix(instruments);
    const n = ym.nrow, k1 = em.ncol, mz = zm.ncol, k2 = exog.nrow == 0 ? 0 : exog.ncol;

    enforce!GretlException(ym.ncol == 1, "tsls: y must have one column, it has " ~ to!string(ym.ncol));
    enforce!GretlException(em.nrow == n && zm.nrow == n && (k2 == 0 || exog.nrow == n), "tsls: y, endog, instruments and exog need the same number of rows");
    enforce!GretlException(mz >= k1, "tsls: " ~ to!string(mz) ~ " instruments are not enough for " ~ to!string(k1) ~ " endogenous regressors");
    enforce!GretlException(n > k1 + k2 + mz + (constant ? 1 : 0), "tsls: too few observations");
    requireFinite(ym, "tsls: y");
    requireFinite(em, "tsls: endog");
    requireFinite(zm, "tsls: instruments");
    if (k2 > 0) requireFinite(exog, "tsls: exog");

    auto xw = (k2 > 0) ? columnsOf(exog) : (const(double)[][]).init;
    auto cols = columnsOf(ym) ~ xw ~ columnsOf(em) ~ columnsOf(zm);
    auto names = ["const", "y"] ~ numbered("w", k2) ~ numbered("x", k1) ~ numbered("z", mz);
    auto ds = makeDataset(cols, names);
    int[] spec = [1];

    if (constant) spec ~= 0;
    foreach (j; 0 .. k2 + k1) spec ~= cast(int) (2 + j);
    spec ~= LISTSEP;
    if (constant) spec ~= 0;
    foreach (j; 0 .. k2) spec ~= cast(int) (2 + j);
    foreach (j; 0 .. mz) spec ~= cast(int) (2 + k2 + k1 + j);
    return runLsq(ds, spec, GretlCmdIndex.ivreg, "2SLS", Robust(vcov, lag, kernel, prewhiten, null), names, n);
}

/* logit(y, x) and probit(y, x): maximum likelihood estimation of a binary choice model: P(y = 1) = F(x * coef) with F the logistic (logit) or
   normal (probit) distribution function; y holds 0s and 1s (n values), x (n x k) the regressors, and a constant is added unless constant is false.
   `robust` gives the quasi-ML (sandwich) covariance matrix. The result has the log-likelihood, the information criteria, McFadden's R squared (in rsq),
   the z-ratios and normal p-values (df is 0), `fitted` is the estimated probability of y = 1 and `resid` the generalized residuals of the library. */
Estimate logit(Y, X)(Y y, X x, bool constant = true, bool robust = false) @trusted if (isOperand!Y && isOperand!X) {
    return binaryChoice(y, x, constant, robust, true);
}

/// ditto
Estimate probit(Y, X)(Y y, X x, bool constant = true, bool robust = false) @trusted if (isOperand!Y && isOperand!X) {
    return binaryChoice(y, x, constant, robust, false);
}

private Estimate binaryChoice(Y, X)(Y y, X x, bool constant, bool robust, bool isLogit) @trusted {
    const who = isLogit ? "logit" : "probit";
    auto ym = asMatrix(y), xm = asMatrix(x);
    const n = ym.nrow, k = xm.ncol;

    enforce!GretlException(ym.ncol == 1, who ~ ": y must have one column, it has " ~ to!string(ym.ncol));
    enforce!GretlException(xm.nrow == n, who ~ ": y has " ~ to!string(n) ~ " observations but x has " ~ to!string(xm.nrow) ~ " rows");
    enforce!GretlException(n > k + (constant ? 1 : 0), who ~ ": too few observations");
    requireFinite(ym, who ~ ": y");
    requireFinite(xm, who ~ ": x");
    foreach (i; 0 .. n) {
        enforce!GretlException(ym[i, 0] == 0.0 || ym[i, 0] == 1.0, who ~ ": y must be 0 or 1, but observation " ~ to!string(i) ~ " is " ~ to!string(ym[i, 0]));
    }
    auto names = ["const", "y"] ~ numbered("x", k);
    auto ds = makeDataset(columnsOf(ym) ~ columnsOf(xm), names);
    auto raw = ds.ptr;
    int[] spec = [1];

    if (constant) spec ~= 0;
    foreach (j; 0 .. k) spec ~= cast(int) (2 + j);
    auto list = makeList(spec);

    enforce!GretlException(push_program_state() == 0, "could not save the library settings");
    scope (exit) pop_program_state();
    Model m = isLogit ? binary_logit(list, raw, robust ? Gretlopt.OPT_R : Gretlopt.OPT_NONE, null)
                      : binary_probit(list, raw, robust ? Gretlopt.OPT_R : Gretlopt.OPT_NONE, null);

    scope (exit) {
        clear_model(&m);
        import core.stdc.stdlib: free;
        free(list);
    }
    enforce!GretlException(m.errcode == 0, who ~ ": the estimation failed, " ~ errorText(m));
    auto e = fromModel(m, isLogit ? "Logit" : "Probit", robust ? "QML (sandwich)" : "ML (Hessian)", names, n);

    // maximum likelihood: z-ratios with the normal distribution
    foreach (i; 0 .. e.ncoef) {
        const z = e.tstat.ptr.val[i];

        e.pvalue.ptr.val[i] = (e.se.ptr.val[i] > 0 && isFinite(z)) ? 2.0 * ndtr(-fabs(z)) : double.nan;
    }
    e.df = 0;
    return e;
}

/* ---- estimators on a Dataset: series are chosen by name ---- */

/* replace the generated names (x1, x2, ... and w1, w2, ...) by the names of the series */
private void renameCoefficients(ref Estimate e, string[] xnames, string[] wnames = null) {
    foreach (ref nm; e.names) {
        foreach (j, s; xnames) {
            if (nm == "x" ~ to!string(j + 1)) nm = s;
        }
        foreach (j, s; wnames) {
            if (nm == "w" ~ to!string(j + 1)) nm = s;
        }
    }
}

/* ols(d, "y", ["x1", "x2"]): the regression of the series y on the series x1 and x2 of the Dataset d and a constant, with the options of the matrix form
   (constant, vcov, lag, kernel, prewhiten) and `cluster`, here the NAME of a series. Observations with a missing value in any series used are left out
   (the result says how many were used, in nobs); the coefficients carry the names of the series. */
Estimate ols(Dataset d, string y, string[] x, bool constant = true, Vcov vcov = Vcov.classical, int lag = 0, HACKernel kernel = HACKernel.bartlett,
             bool prewhiten = false, string cluster = null) @trusted {
    auto s = selectSeries(d, [[y], x] ~ (cluster is null ? null : [[cluster]]));
    auto e = ols(matrixFrom(s.cols[0 .. 1]), matrixFrom(s.cols[1 .. 1 + x.length]), constant, vcov, lag, kernel, prewhiten,
                 cluster is null ? null : s.cols[$ - 1]);

    renameCoefficients(e, x);
    return e;
}

/* wls(d, "y", ["x1", "x2"], "w"): weighted least squares with the series w as the weights (see the matrix form). */
Estimate wls(Dataset d, string y, string[] x, string weights, bool constant = true, Vcov vcov = Vcov.classical, int lag = 0,
             HACKernel kernel = HACKernel.bartlett, bool prewhiten = false, string cluster = null) @trusted {
    auto s = selectSeries(d, [[y], x, [weights]] ~ (cluster is null ? null : [[cluster]]));
    auto e = wls(matrixFrom(s.cols[0 .. 1]), matrixFrom(s.cols[1 .. 1 + x.length]), s.cols[1 + x.length], constant, vcov, lag, kernel, prewhiten,
                 cluster is null ? null : s.cols[$ - 1]);

    renameCoefficients(e, x);
    return e;
}

/* tsls(d, "y", ["endog"], ["z1", "z2"], exog: ["w"]): two-stage least squares with the series named: the endogenous regressors, the excluded
   instruments and the included exogenous regressors (see the matrix form). */
Estimate tsls(Dataset d, string y, string[] endog, string[] instruments, string[] exog = null, bool constant = true, Vcov vcov = Vcov.classical,
              int lag = 0, HACKernel kernel = HACKernel.bartlett, bool prewhiten = false) @trusted {
    auto s = selectSeries(d, [[y], endog, instruments, exog]);
    const ne = endog.length, nz = instruments.length;
    auto exogM = exog.length ? matrixFrom(s.cols[1 + ne + nz .. $]) : GretlMatrix.init;
    auto e = tsls(matrixFrom(s.cols[0 .. 1]), matrixFrom(s.cols[1 .. 1 + ne]), matrixFrom(s.cols[1 + ne .. 1 + ne + nz]), exogM, constant, vcov, lag, kernel, prewhiten);

    renameCoefficients(e, endog, exog);
    return e;
}

/* logit(d, "y", ["x1", "x2"]) and probit(...): binary choice models on the series of a Dataset (see the matrix forms). */
Estimate logit(Dataset d, string y, string[] x, bool constant = true, bool robust = false) @trusted {
    auto s = selectSeries(d, [[y], x]);
    auto e = logit(matrixFrom(s.cols[0 .. 1]), matrixFrom(s.cols[1 .. $]), constant, robust);

    renameCoefficients(e, x);
    return e;
}

/// ditto
Estimate probit(Dataset d, string y, string[] x, bool constant = true, bool robust = false) @trusted {
    auto s = selectSeries(d, [[y], x]);
    auto e = probit(matrixFrom(s.cols[0 .. 1]), matrixFrom(s.cols[1 .. $]), constant, robust);

    renameCoefficients(e, x);
    return e;
}
