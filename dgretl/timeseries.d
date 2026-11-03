/* Time series: correlograms, unit-root and stationarity tests, filters and transformations, as a D-friendly layer over the port (gretl.describe,
   gretl.adf_kpss, gretl.genfuncs). Series come as arrays, GretlMatrix or Submatrix (a series is one column), results are plain structs of
   GretlMatrix and numbers, and every function throws GretlException when the library reports an error. Missing (NaN) values are refused.

   import dgretl;
   auto c = correlogram(x, 12);                  // c.acf, c.pacf: 12 x 1
   auto t = adf(x, lags: 2, trend: Trend.constant);
   writeln(t.statistic, " ", t.pvalue);
   auto cycle = hpCycle(x, 1600);                // filters return the filtered series as an n x 1 GretlMatrix

   Time information, where a function needs it (seasonal dummies in adf, the calendar of a VAR), is given by the optional arguments `frequency`
   (observations per year: 1, 4, 12, ...) and `start` (the first observation, as "1990:1" or "1990Q1"); without them the series is annual from 1 (or from 1:1 when only a frequency is given). */

module dgretl.timeseries;

public import dgretl.matrix;
public import dgretl.exception;

import std.conv: to;
import std.typecons: RefCountedAutoInitialize, SafeRefCounted;
import std.exception: enforce;
import std.format: format;
import std.math: isFinite, isNaN;
import std.string: fromStringz, toStringz;

import dgretl.dataset: Dataset;
import dgretl.estimate: Estimate, fromModel, ols;
import dgretl.internal;
import gretl.adf_kpss: adf_test, engle_granger_test, kpss_test;
import gretl.estimate: ar1_model, ar_model, arch_model;
import gretl.gretl_model: clear_model, gretl_model_get_double;
import gretl.libgretlh: Model;
import gretl.gretl_listh: LISTSEP;
import gretl.describe: acf_matrix, ljung_box;
import gretl.genfuncs: bkbp_filter, butterworth_filter, fracdiff_series, hp_filter;
import gretl.gretl_commandsh: GretlCmdIndex;
import gretl.gretl_errors: gretl_errmsg_get;
import gretl.gretl_list: gretl_list_new;
import gretl.gretl_matrix: gretl_matrix_free;
import gretl.gretl_matrixh: GretlMatrixRaw;
import gretl.gretl_prn: Prn, gretl_print_destroy, gretl_print_get_buffer, gretl_print_new;
import gretl.gretl_prnh: PrnType;
import gretl.gretl_utils: get_last_pvalue, get_last_test_statistic, simple_set_obs;
import gretl.libgretlh: Gretlopt;
import gretl.libset: pop_program_state, push_program_state;
import gretl.options: set_optval_string;
import gretl.pvalues: chisq_cdf_comp, snedecor_cdf_comp;
import gretl.dataio: dateton;
import gretl.dataset: dataset_add_observations;
import gretl.libgretlh: DatasetRaw, TsCodes;
import gretl.random: gretl_rand_set_seed;
import gretl.describe: xcf_vec, auto_spectrum_order, fractint, get_last_result_data, periodogram_matrix;
import gretl.libgretlh: GretlType;
import gretl.pvalues: student_pvalue_2;
import gretl.plugin.arma: arma_model;
import gretl.plugin.johansen: vecm_test_restriction;
import gretl.gretl_restricth: GretlRestriction;
import gretl.gretl_utils: get_last_pvals_matrix, get_last_test_matrix;
import gretl.var: VarTask, johansen_wrapper, johansen_test;
import gretl.johansenh: JohansenCode;
import gretl.plugin.garch: garch_model;
import cephes.mconfh: CephesCmplx = Cmplx;
import gretl.gretl_model: gretl_model_get_data, gretl_model_get_int;
import gretl.var: default_VAR_horizon, gretl_VAR, gretl_VAR_free, gretl_VAR_get_FEVD_matrix, gretl_VAR_get_forecast_matrix,
    gretl_VAR_get_impulse_response, gretl_VAR_get_model, gretl_VAR_get_roots, gretl_VAR_param_names, gui_VAR_lagsel;
import gretl.varh: GretlVar, var_n_lags;

/* ---- helpers ---- */

/* a series (one column) as a plain array, checked for finite values */
private double[] seriesOf(X)(X x, string who) @trusted if (isOperand!X) {
    auto m = cbind(x);

    enforce!GretlException(m.ncol == 1, who ~ ": a series has one column, this has " ~ to!string(m.ncol));
    requireFinite(m, who);
    return m.ptr.val[0 .. m.nrow].dup;
}

private GretlMatrix columnFrom(const(double)[] v) @trusted {
    auto m = fresh(v.length, 1);

    foreach (i, x; v) m.ptr.val[i] = x;
    return m;
}

/* a temporary dataset holding the series `cols` (named "x1", "x2", ...) with the calendar given */
private Dataset timeDataset(const(double)[][] cols, int frequency, string start, string who, string[] names = null) @trusted {
    enforce!GretlException(frequency >= 1, who ~ ": the frequency must be at least 1");
    auto ds = makeDataset(cols, ["const"] ~ (names is null ? numbered("x", cols.length) : names));
    if (start == "1" && frequency > 1) start = "1:1";     // the default start of a sub-annual series
    const err = simple_set_obs(ds.ptr, frequency, start.toStringz, Gretlopt.OPT_T);

    enforce!GretlException(err == 0, who ~ ": cannot use the starting observation \"" ~ start ~ "\" with frequency " ~ to!string(frequency) ~ ": " ~ lastError());
    return ds;
}

private string lastError() @trusted {
    return fromStringz(gretl_errmsg_get()).idup;
}

/* ---- correlograms ---- */

/* The autocorrelation and partial autocorrelation functions at lags 1 .. maxLag, each maxLag x 1. */
struct Correlogram {
    GretlMatrix acf;
    GretlMatrix pacf;
    GretlMatrix ljungBox;   // maxLag x 1: the Ljung-Box statistic Q(k) for k = 1 .. maxLag
}

/* correlogram(x, maxLag): the sample ACF and PACF of the series x at lags 1 .. maxLag (library: acf_matrix), and the Ljung-Box statistics Q(1) .. Q(maxLag).
   maxLag must be between 1 and nobs - 1. */
Correlogram correlogram(X)(X x, int maxLag) @trusted if (isOperand!X) {
    auto v = seriesOf(x, "correlogram");
    const n = v.length;

    enforce!GretlException(maxLag >= 1 && maxLag < cast(int) n, "correlogram: maxLag must be between 1 and " ~ to!string(n - 1));
    int err = 0;
    auto raw = acf_matrix(v.ptr, maxLag, null, cast(int) n, &err);

    enforce!GretlException(err == 0 && raw !is null, "correlogram: " ~ lastError());
    scope (exit) gretl_matrix_free(raw);
    Correlogram c;

    c.acf = fresh(maxLag, 1);
    c.pacf = fresh(maxLag, 1);
    c.ljungBox = fresh(maxLag, 1);
    foreach (k; 0 .. maxLag) {
        c.acf.ptr.val[k] = raw.val[k];
        c.pacf.ptr.val[k] = raw.val[maxLag + k];
        int e2;

        c.ljungBox.ptr.val[k] = ljung_box(k + 1, 0, cast(int) n - 1, v.ptr, &e2);
        enforce!GretlException(e2 == 0, "correlogram: the Ljung-Box statistic failed");
    }
    return c;
}

/* The result of a test with a statistic, a p-value and degrees of freedom (the chi-square tests of this module). */
struct TestResult {
    double statistic;
    double pvalue;
    int df;
}

alias PortmanteauTest = TestResult;

/* ljungBox(x, lags): the Ljung-Box test of no autocorrelation up to `lags`; the p-value is from chi-square with `lags - fitted` degrees of freedom
   (set `fitted` to the number of ARMA parameters estimated when x holds residuals of a fitted model). */
PortmanteauTest ljungBox(X)(X x, int lags, int fitted = 0) @trusted if (isOperand!X) {
    auto v = seriesOf(x, "ljungBox");
    const n = v.length;

    enforce!GretlException(lags >= 1 && lags < cast(int) n, "ljungBox: lags must be between 1 and " ~ to!string(n - 1));
    enforce!GretlException(fitted >= 0 && fitted < lags, "ljungBox: fitted must be between 0 and lags - 1");
    int err = 0;
    const q = ljung_box(lags, 0, cast(int) n - 1, v.ptr, &err);

    enforce!GretlException(err == 0 && !isNaN(q), "ljungBox: " ~ lastError());
    PortmanteauTest t;

    t.statistic = q;
    t.df = lags - fitted;
    t.pvalue = chisq_cdf_comp(t.df, q);
    return t;
}

/* The cross-correlation function of two series at the lags -maxLag .. maxLag. */
struct CrossCorrelogram {
    GretlMatrix lags;           // (2 maxLag + 1) x 1: -maxLag .. maxLag
    GretlMatrix ccf;            // (2 maxLag + 1) x 1: the correlation of x(t) with y(t - k) at lag k
}

/* crossCorrelogram(x, y, maxLag): the sample cross-correlations r(k) = sum (x(t) - xbar) (y(t - k) - ybar) / sqrt(sum (x - xbar)^2 sum (y - ybar)^2) for k = -maxLag .. maxLag
   (the library's `xcorrgm`). The series must have the same length, no missing values and not be constant; 2 maxLag must not exceed n. */
CrossCorrelogram crossCorrelogram(X, Y)(X x, Y y, int maxLag) @trusted if (isOperand!X && isOperand!Y) {
    auto a = seriesOf(x, "crossCorrelogram"), b = seriesOf(y, "crossCorrelogram");

    enforce!GretlException(a.length == b.length, "crossCorrelogram: the series have " ~ to!string(a.length) ~ " and " ~ to!string(b.length) ~ " observations");
    enforce!GretlException(maxLag >= 1 && 2 * maxLag <= cast(int) a.length, "crossCorrelogram: maxLag must be between 1 and n / 2");
    int err = 0;
    auto m = xcf_vec(a.ptr, b.ptr, maxLag, null, cast(int) a.length, &err);

    enforce!GretlException(err == 0 && m !is null, "crossCorrelogram: " ~ lastError());
    scope (exit) gretl_matrix_free(m);
    CrossCorrelogram r;
    const h = 2 * maxLag + 1;

    r.lags = fresh(h, 1);
    r.ccf = fresh(h, 1);
    foreach (i; 0 .. h) {
        r.lags.ptr.val[i] = i - maxLag;
        r.ccf.ptr.val[i] = m.val[i];
    }
    return r;
}

/// ditto, for two series of a Dataset
CrossCorrelogram crossCorrelogram(Dataset d, string xname, string yname, int maxLag) @trusted {
    return crossCorrelogram(named(d, xname, "crossCorrelogram"), named(d, yname, "crossCorrelogram"), maxLag);
}

/* ---- unit-root and stationarity tests ---- */

/* The deterministic terms of the test regression. */
enum Trend {
    none,           // no constant
    constant,       // a constant
    trend,          // constant and linear trend
    quadratic       // constant, linear and quadratic trend (adf only)
}

/* How the lag order of the augmented regression is chosen by `testDown`. */
enum TestDown {
    off,            // use `lags` as given
    tstat,          // start from `lags` and drop the last lag while its t-ratio is not significant (10% level)
    aic,            // the lag order, from 0 to `lags`, that minimizes AIC
    bic             // ... BIC
}

/* The result of a unit-root or stationarity test. `report` is the text the library prints for it. */
struct UnitRootTest {
    string test;
    double statistic;
    double pvalue;      // NaN when the library has no p-value for the case (the KPSS test gives critical values in the report only)
    string report;
}

private Gretlopt trendOpt(Trend t) {
    final switch (t) {
        case Trend.none: return Gretlopt.OPT_N;
        case Trend.constant: return Gretlopt.OPT_C;
        case Trend.trend: return Gretlopt.OPT_T;
        case Trend.quadratic: return Gretlopt.OPT_R;
    }
}

/* adf(x, lags): the augmented Dickey-Fuller test of a unit root in the series x with `lags` lagged differences.
       trend      the deterministic terms (Trend.constant by default)
       testDown   choose the lag order from 0 .. lags by a criterion (the report shows the choice)
       gls        Elliott-Rothenberg-Stock GLS detrending (DF-GLS); trend must be constant or trend
       seasonals  add seasonal dummies (needs frequency 4 or 12 etc.)
       difference apply the test to the first difference of x
       frequency, start   the calendar of the series (see the module comment) */
UnitRootTest adf(X)(X x, int lags = 0, Trend trend = Trend.constant, TestDown testDown = TestDown.off, bool gls = false, bool seasonals = false,
                    bool difference = false, int frequency = 1, string start = "1") @trusted if (isOperand!X) {
    auto v = seriesOf(x, "adf");

    enforce!GretlException(lags >= 0, "adf: lags must not be negative");
    enforce!GretlException(v.length > cast(size_t) (lags + 10), "adf: too few observations for " ~ to!string(lags) ~ " lags");
    auto ds = timeDataset([v], frequency, start, "adf");
    auto opt = trendOpt(trend);

    if (testDown != TestDown.off) {
        opt |= Gretlopt.OPT_E;
        set_optval_string(GretlCmdIndex.adf, Gretlopt.OPT_E, ["", "tstat", "AIC", "BIC"][testDown].toStringz);
    }
    if (gls) opt |= Gretlopt.OPT_G;
    if (seasonals) opt |= Gretlopt.OPT_D;
    if (difference) opt |= Gretlopt.OPT_F;
    auto list = makeList([1, 1]);
    auto pr = Printer(0);

    enforce!GretlException(push_program_state() == 0, "could not save the library settings");
    scope (exit) {
        pop_program_state();
        import core.stdc.stdlib: free;
        free(list);
    }
    const err = adf_test(lags, list, ds.ptr, opt, pr.prn);

    enforce!GretlException(err == 0, "adf: the test failed, " ~ lastError());
    return UnitRootTest(gls ? "DF-GLS" : "ADF", get_last_test_statistic(), get_last_pvalue(), pr.text);
}

/* kpss(x, lags): the KPSS test of the null that x is stationary (around a constant, or around a trend with Trend.trend), with the Bartlett
   window of `lags` lags. Use at least 1; a common choice is 4 * (n / 100)^0.25. Only Trend.constant and Trend.trend are meaningful. */
UnitRootTest kpss(X)(X x, int lags, Trend trend = Trend.constant, bool difference = false) @trusted if (isOperand!X) {
    auto v = seriesOf(x, "kpss");

    enforce!GretlException(trend == Trend.constant || trend == Trend.trend, "kpss: trend must be Trend.constant or Trend.trend");
    enforce!GretlException(lags >= 0 && v.length > cast(size_t) (lags + 2), "kpss: lags must be between 0 and nobs - 3");
    auto ds = timeDataset([v], 1, "1", "kpss");
    auto opt = Gretlopt.OPT_NONE;

    if (trend == Trend.trend) opt |= Gretlopt.OPT_T;
    if (difference) opt |= Gretlopt.OPT_F;
    auto list = makeList([1, 1]);
    auto pr = Printer(0);

    enforce!GretlException(push_program_state() == 0, "could not save the library settings");
    scope (exit) {
        pop_program_state();
        import core.stdc.stdlib: free;
        free(list);
    }
    const err = kpss_test(lags, list, ds.ptr, opt, pr.prn);

    enforce!GretlException(err == 0, "kpss: the test failed, " ~ lastError());
    return UnitRootTest("KPSS", get_last_test_statistic(), get_last_pvalue(), pr.text);
}

/* ---- filters ---- */

private GretlMatrix runFilter(string who, const(double)[] v, int frequency, scope int delegate(const(double) * x, double * y, void * ds) @system call) @trusted {
    auto ds = timeDataset([v], frequency, "1", who);
    auto y = new double[v.length];

    y[] = double.nan;
    const err = call(v.ptr, y.ptr, ds.ptr);

    enforce!GretlException(err == 0, who ~ ": the filter failed, " ~ lastError());
    return columnFrom(y);
}

/* hpTrend(x, lambda): the trend of the Hodrick-Prescott filter, an n x 1 matrix; the cycle is x - trend. lambda is the smoothing parameter (1600 for
   quarterly data, 100 for annual, 14400 for monthly are the usual choices). */
GretlMatrix hpTrend(X)(X x, double lambda = 1600) @trusted if (isOperand!X) {
    auto v = seriesOf(x, "hpTrend");

    enforce!GretlException(lambda > 0 && v.length >= 4, "hpTrend: lambda must be positive and there must be at least 4 observations");
    return runFilter("hpTrend", v, 1, (const(double) * p, double * y, void * ds) {
        import gretl.libgretlh: DatasetRaw;
        return hp_filter(p, y, cast(DatasetRaw *) ds, lambda, Gretlopt.OPT_T);   // OPT_T: the trend (the default is the cycle)
    });
}

/* hpCycle(x, lambda): the cyclical component x - trend of the Hodrick-Prescott filter (see hpTrend). */
GretlMatrix hpCycle(X)(X x, double lambda = 1600) @trusted if (isOperand!X) {
    auto v = seriesOf(x, "hpCycle");
    auto t = hpTrend(v, lambda);
    auto r = fresh(v.length, 1);

    foreach (i; 0 .. v.length) r.ptr.val[i] = v[i] - t.ptr.val[i];
    return r;
}

/* bkFilter(x, lower, upper, k): the Baxter-King band-pass filter keeping the cycles with a period between `lower` and `upper` observations
   (for quarterly data 6 and 32 are the usual choices) with a moving average of k leads and lags (12 for quarterly data). The first and last k values are
   NaN. */
GretlMatrix bkFilter(X)(X x, int lower = 6, int upper = 32, int k = 12) @trusted if (isOperand!X) {
    auto v = seriesOf(x, "bkFilter");

    enforce!GretlException(lower >= 2 && upper > lower, "bkFilter: need 2 <= lower < upper");
    enforce!GretlException(k >= 1 && v.length > cast(size_t) (2 * k + 1), "bkFilter: k must be at least 1 and the series longer than 2k + 1");
    return runFilter("bkFilter", v, 1, (const(double) * p, double * y, void * ds) {
        import gretl.libgretlh: DatasetRaw;
        return bkbp_filter(p, y, cast(DatasetRaw *) ds, lower, upper, k);
    });
}

/* butterworth(x, order, cutoff): the low-pass Butterworth filter of the given order; cutoff is the cutoff frequency in degrees (0 < cutoff < 180). */
GretlMatrix butterworth(X)(X x, int order = 2, double cutoff = 67.5) @trusted if (isOperand!X) {
    auto v = seriesOf(x, "butterworth");

    enforce!GretlException(order >= 1 && cutoff > 0 && cutoff < 180, "butterworth: order must be positive and cutoff between 0 and 180 degrees");
    return runFilter("butterworth", v, 1, (const(double) * p, double * y, void * ds) {
        import gretl.libgretlh: DatasetRaw;
        return butterworth_filter(p, y, cast(DatasetRaw *) ds, order, cutoff);
    });
}

/* fracDiff(x, d): the fractional difference (1 - L)^d of x (d may be negative to integrate); the series keeps its length. */
GretlMatrix fracDiff(X)(X x, double d) @trusted if (isOperand!X) {
    auto v = seriesOf(x, "fracDiff");

    return runFilter("fracDiff", v, 1, (const(double) * p, double * y, void * ds) {
        import gretl.libgretlh: DatasetRaw;
        return fracdiff_series(p, y, d, 1, -1, cast(DatasetRaw *) ds);
    });
}

/* ---- transformations ---- */

/* lag(x, k): x shifted k observations later (a lag for k > 0, a lead for k < 0), the vacated positions are NaN; same length as x. */
GretlMatrix lag(X)(X x, int k = 1) @trusted if (isOperand!X) {
    auto m = cbind(x);
    auto r = fresh(m.nrow, m.ncol);

    foreach (j; 0 .. m.ncol) {
        foreach (i; 0 .. m.nrow) {
            const s = cast(long) i - k;

            r.ptr.val[j * m.nrow + i] = (s >= 0 && s < cast(long) m.nrow) ? m[cast(size_t) s, j] : double.nan;
        }
    }
    return r;
}

/* diff(x, order, seasonal): the difference x(t) - x(t - s), repeated `order` times, where s = `seasonal` (1 for the ordinary first difference);
   the first order * s values are NaN; same length as x. */
GretlMatrix diff(X)(X x, int order = 1, int seasonal = 1) @trusted if (isOperand!X) {
    enforce!GretlException(order >= 1 && seasonal >= 1, "diff: order and seasonal must be positive");
    auto m = cbind(x);
    auto r = fresh(m.nrow, m.ncol);

    foreach (j; 0 .. m.ncol) {
        auto col = new double[m.nrow];

        foreach (i; 0 .. m.nrow) col[i] = m[i, j];
        foreach (_; 0 .. order) {
            auto next = new double[m.nrow];

            foreach (i; 0 .. m.nrow) next[i] = (i >= cast(size_t) seasonal) ? col[i] - col[i - seasonal] : double.nan;
            col = next;
        }
        foreach (i; 0 .. m.nrow) r.ptr.val[j * m.nrow + i] = col[i];
    }
    return r;
}

/* ---- tests for conditional heteroskedasticity and cointegration ---- */

/* archTest(u, order): the LM test of Engle for ARCH effects up to `order` in the series u (usually residuals): regress u^2 on a constant and its own
   `order` lags; the statistic is (n - order) times the R squared, chi-square with `order` degrees of freedom under the null of none. */
TestResult archTest(U)(U u, int order) @trusted if (isOperand!U) {
    auto v = seriesOf(u, "archTest");
    const n = v.length;

    enforce!GretlException(order >= 1 && order < cast(int) n - 2, "archTest: order must be between 1 and nobs - 3");
    const m = n - order;
    auto y = new double[m];
    auto x = fresh(m, order);

    foreach (t; 0 .. m) {
        y[t] = v[t + order] ^^ 2;
        foreach (i; 0 .. order) x.ptr.val[i * m + t] = v[t + order - 1 - i] ^^ 2;
    }
    auto fit = ols(y, x);
    TestResult r;

    r.df = order;
    r.statistic = m * fit.rsq;
    r.pvalue = chisq_cdf_comp(order, r.statistic);
    return r;
}

/* engleGranger(y, x, lags): the Engle-Granger test of cointegration between y (n x 1) and the columns of x: the cointegrating regression of y on x
   (with the deterministic terms of `trend`: Trend.none, constant (default), trend or quadratic), then an ADF test of its residuals with `lags` lags.
   The p-value is the one of the library (MacKinnon); the report has the regression and the unit-root tests of the individual series. */
UnitRootTest engleGranger(Y, X)(Y y, X x, int lags = 0, Trend trend = Trend.constant, TestDown testDown = TestDown.off) @trusted
        if (isOperand!Y && isOperand!X) {
    auto ym = cbind(y), xm = cbind(x);

    enforce!GretlException(ym.ncol == 1 && xm.nrow == ym.nrow, "engleGranger: y must be a column and x must have as many rows");
    requireFinite(ym, "engleGranger: y");
    requireFinite(xm, "engleGranger: x");
    enforce!GretlException(lags >= 0 && ym.nrow > cast(size_t) (lags + 10 + xm.ncol), "engleGranger: too few observations");
    auto ds = timeDataset(columnsOf(ym) ~ columnsOf(xm), 1, "1", "engleGranger");
    int[] spec = [1];

    foreach (j; 0 .. xm.ncol) spec ~= cast(int) (2 + j);
    auto list = makeList(spec);
    auto opt = Gretlopt.OPT_S;

    final switch (trend) {
        case Trend.none: opt |= Gretlopt.OPT_N; break;
        case Trend.constant: break;
        case Trend.trend: opt |= Gretlopt.OPT_T; break;
        case Trend.quadratic: opt |= Gretlopt.OPT_R; break;
    }
    if (testDown != TestDown.off) {
        opt |= Gretlopt.OPT_E;
        set_optval_string(GretlCmdIndex.coint, Gretlopt.OPT_E, ["", "tstat", "AIC", "BIC"][testDown].toStringz);
    }
    auto pr = Printer(0);

    enforce!GretlException(push_program_state() == 0, "could not save the library settings");
    scope (exit) {
        pop_program_state();
        import core.stdc.stdlib: free;
        free(list);
    }
    const err = engle_granger_test(lags, list, ds.ptr, opt, pr.prn);

    enforce!GretlException(err == 0, "engleGranger: the test failed, " ~ lastError());
    return UnitRootTest("Engle-Granger", get_last_test_statistic(), get_last_pvalue(), pr.text);
}

/* ---- single-equation models with autoregressive errors or ARCH ---- */

/* The estimator for the AR(1) error model. */
enum Ar1Method {
    cochraneOrcutt,     // iterated Cochrane-Orcutt
    praisWinsten,       // iterated Prais-Winsten (keeps the first observation)
    hildrethLu          // Hildreth-Lu grid search
}

private struct ModelInput {
    Dataset ds;
    int[] spec;
    string[] names;
    size_t n;
}

/* y on x (and a constant) in a temporary time-series dataset; spec is the list of the regression: y, [0,] x... */
private ModelInput modelInput(Y, X)(Y y, X x, bool constant, string who) @trusted {
    auto ym = cbind(y), xm = cbind(x);
    const n = ym.nrow, k = xm.ncol;

    enforce!GretlException(ym.ncol == 1, who ~ ": y must have one column");
    enforce!GretlException(xm.nrow == n, who ~ ": y has " ~ to!string(n) ~ " observations but x has " ~ to!string(xm.nrow) ~ " rows");
    enforce!GretlException(n > k + 5, who ~ ": too few observations");
    requireFinite(ym, who ~ ": y");
    requireFinite(xm, who ~ ": x");
    ModelInput mi;

    mi.names = ["const", "y"] ~ numbered("x", k);
    mi.ds = timeDataset(columnsOf(ym) ~ columnsOf(xm), 1, "1", who);
    mi.spec = [1];
    if (constant) mi.spec ~= 0;
    foreach (j; 0 .. k) mi.spec ~= cast(int) (2 + j);
    mi.n = n;
    return mi;
}

/* ar1(y, x): regression of y on x (and a constant) with first-order autoregressive errors, e(t) = rho e(t-1) + u(t). The rho is in `arcoef` (the
   library gives no standard error for it, `arse` is NaN); the other fields are those of the regression on the transformed data. method: Cochrane-Orcutt (default), Prais-Winsten or
   Hildreth-Lu. */
Estimate ar1(Y, X)(Y y, X x, bool constant = true, Ar1Method method = Ar1Method.cochraneOrcutt) @trusted if (isOperand!Y && isOperand!X) {
    auto mi = modelInput(y, x, constant, "ar1");
    auto list = makeList(mi.spec);
    auto opt = Gretlopt.OPT_NONE;

    final switch (method) {
        case Ar1Method.cochraneOrcutt: break;
        case Ar1Method.praisWinsten: opt = Gretlopt.OPT_P; break;
        case Ar1Method.hildrethLu: opt = Gretlopt.OPT_H; break;
    }
    auto pr = Printer(0);

    enforce!GretlException(push_program_state() == 0, "could not save the library settings");
    scope (exit) pop_program_state();
    Model m = ar1_model(list, mi.ds.ptr, opt, pr.prn);

    scope (exit) {
        clear_model(&m);
        import core.stdc.stdlib: free;
        free(list);
    }
    enforce!GretlException(m.errcode == 0, "ar1: the estimation failed, " ~ errorText(m));
    auto e = fromModel(m, ["Cochrane-Orcutt", "Prais-Winsten", "Hildreth-Lu"][method], "classical", mi.names, mi.n);

    // the library keeps the rho of an AR(1) model in the model's data, without a standard error (`rho`, the residual autocorrelation, is something else)
    e.arcoef = fresh(1, 1);
    e.arcoef.ptr.val[0] = gretl_model_get_double(&m, "rho_gls");
    e.arse = fresh(1, 1);
    e.arse.ptr.val[0] = double.nan;
    e.arlags = [1];
    return e;
}

/* ar(y, x, lags): regression of y on x (and a constant) with autoregressive errors at the given lags (an array such as [1, 2] or [1, 4], or a single
   int p for lags 1 .. p), estimated by the generalized Cochrane-Orcutt iteration. The error coefficients are in `arcoef` and `arse`, their lags in `arlags`. */
Estimate ar(Y, X)(Y y, X x, int[] lags, bool constant = true) @trusted if (isOperand!Y && isOperand!X) {
    enforce!GretlException(lags.length >= 1, "ar: give at least one lag");
    foreach (l; lags) enforce!GretlException(l >= 1, "ar: lags must be positive");
    auto mi = modelInput(y, x, constant, "ar");
    int[] full = lags ~ [LISTSEP] ~ mi.spec;
    auto list = makeList(full);
    auto pr = Printer(0);

    enforce!GretlException(push_program_state() == 0, "could not save the library settings");
    scope (exit) pop_program_state();
    Model m = ar_model(list, mi.ds.ptr, Gretlopt.OPT_NONE, pr.prn);

    scope (exit) {
        clear_model(&m);
        import core.stdc.stdlib: free;
        free(list);
    }
    enforce!GretlException(m.errcode == 0, "ar: the estimation failed, " ~ errorText(m));
    return fromModel(m, "AR", "classical", mi.names, mi.n);
}

/// ditto
Estimate ar(Y, X)(Y y, X x, int p, bool constant = true) @trusted if (isOperand!Y && isOperand!X) {
    enforce!GretlException(p >= 1, "ar: the order must be at least 1");
    int[] lags;

    foreach (i; 1 .. p + 1) lags ~= i;
    return ar(y, x, lags, constant);
}

/* arch(y, x, order): the two-step weighted least squares estimate of ARCH(order): OLS of y on x (and a constant), a regression of the squared
   residuals on their `order` lags, then WLS with the fitted variances as weights. The result is the final WLS regression. */
Estimate arch(Y, X)(Y y, X x, int order, bool constant = true) @trusted if (isOperand!Y && isOperand!X) {
    enforce!GretlException(order >= 1, "arch: the order must be at least 1");
    auto mi = modelInput(y, x, constant, "arch");

    enforce!GretlException(mi.n > cast(size_t) (order + mi.spec.length + 3), "arch: too few observations for this order");
    auto list = makeList(mi.spec);

    enforce!GretlException(push_program_state() == 0, "could not save the library settings");
    scope (exit) pop_program_state();
    Model m = arch_model(list, order, mi.ds.ptr, Gretlopt.OPT_NONE);

    scope (exit) {
        clear_model(&m);
        import core.stdc.stdlib: free;
        free(list);
    }
    enforce!GretlException(m.errcode == 0, "arch: the estimation failed, " ~ errorText(m));
    return fromModel(m, "ARCH-WLS", "classical", mi.names, m.nobs);
}

/* ---- vector autoregressions ---- */

/* The estimate of a VAR: coefficient matrices as GretlMatrix, and the library's structure kept (reference counted) so that the analysis functions
   irf, fevd and forecast can use it. Copies share it; it is freed when the last copy goes away. */
struct VarFit {
    private struct Payload {
        GretlVar * raw;
        Dataset ds;                 // the dataset the VAR was estimated on (the forecasts need it)
        int nOrig;                  // its number of observations before any forecast extended it

        this(GretlVar * v, Dataset d) {
            raw = v;
            ds = d;
            nOrig = d.ptr.n;
        }

        ~this() {
            if (raw !is null) {
                gretl_VAR_free(raw);
                raw = null;
            }
        }

        @disable this(this);
    }

    private SafeRefCounted!(Payload, RefCountedAutoInitialize.no) _p;

    string[] names;             // the endogenous variables
    string[] coefNames;         // the rows of `coef`: "const", then each variable's lags (y1_1, y1_2, y2_1, ...), the exogenous variables, "S1" ..., "time"
    int[] lags;                 // the lags used (1 .. order unless a gappy list was given)
    size_t order;               // the maximum lag
    size_t nobs;                // observations used by each equation
    size_t neqns;
    int df;                     // residual degrees of freedom of each equation
    GretlMatrix coef;           // ncoef x neqns, one column per equation
    GretlMatrix se;             // standard errors, ncoef x neqns (robust if requested)
    GretlMatrix sigma;          // neqns x neqns covariance matrix of the residuals (divisor: the equations' df average, as the library does it)
    GretlMatrix resid;          // nobs x neqns
    GretlMatrix roots;          // neqns * order x 3: real part, imaginary part, modulus of the eigenvalues of the companion matrix; stable if all moduli < 1
    GretlMatrix grangerF;       // neqns x neqns: [i, j] is the F statistic that the lags of variable j are all zero in equation i (NaN on the diagonal when the test is not available)
    GretlMatrix grangerP;       // their p-values: F(number of lags, df)
    double lnL = double.nan;
    double aic = double.nan, bic = double.nan, hqc = double.nan;    // the VAR criteria of the library (per observation)
    double ldet = double.nan;   // log determinant of sigma
    string sampleStart, sampleEnd;  // labels of the first and last observation used

    package(dgretl) GretlVar * raw() @trusted {
        return _p.refCountedStore.isInitialized ? _p.refCountedPayload.raw : null;
    }

    package(dgretl) Dataset dataset() @trusted {
        return _p.refCountedPayload.ds;
    }

    package(dgretl) int nOrig() @trusted {
        return _p.refCountedPayload.nOrig;
    }

    /* equation i (0-based) as the Estimate of a regression: coefficients, standard errors, t-ratios, R squared and so on */
    Estimate equation(size_t i) @trusted {
        enforce!GretlException(i < neqns, "VarFit.equation: there are only " ~ to!string(neqns) ~ " equations");
        auto m = *cast(Model *) gretl_VAR_get_model(raw, cast(int) i);
        auto e = fromModel(m, "VAR equation for " ~ names[i], "classical", ["const"] ~ names, nobs);

        e.names = coefNames.dup;
        return e;
    }

    string toString() {
        string s = format("VAR(%d) of %s, %d observations (%s to %s)\n", order, names, nobs, sampleStart, sampleEnd);

        s ~= format("log-likelihood %.4f, AIC %.5f, BIC %.5f, HQC %.5f\n", lnL, aic, bic, hqc);
        foreach (j; 0 .. neqns) {
            s ~= format("\nequation for %s:\n%-14s %12s %12s\n", names[j], "", "coef", "std.err");
            foreach (i, nm; coefNames) s ~= format("%-14s %12.6g %12.6g\n", nm, coef[i, j], se[i, j]);
        }
        return s;
    }
}

/* where a name is not given, "y1", "y2", ... and "x1", ... */
private string[] namesOr(string[] given, size_t count, string stem, string who) {
    if (given is null) return numbered(stem, count);
    enforce!GretlException(given.length == count, who ~ ": " ~ to!string(given.length) ~ " names for " ~ to!string(count) ~ " series");
    return given.dup;
}

/* the observation number (0-based) of a date such as "1990:3" in the dataset's calendar */
private int observation(Dataset ds, string date, string who) @trusted {
    const t = dateton(date.toStringz, ds.ptr);

    enforce!GretlException(t >= 0, who ~ ": the observation \"" ~ date ~ "\" is not in the sample");
    return t;
}

private string observationLabel(Dataset ds, int t) @trusted {
    import gretl.dataio: ntolabel;
    char[64] buf = 0;

    ntolabel(buf.ptr, t, ds.ptr);
    return fromStringz(buf.ptr).idup;
}

/* var(y, order): the vector autoregression of the columns of y (n x k) with lags 1 .. order.
       exog       exogenous regressors (n x m GretlMatrix), same rows as y
       constant   include a constant (default); trend: add a linear trend; seasonals: add seasonal dummies (frequency > 1)
       lags       a gappy list of lags such as [1, 4] (then `order` must be the largest)
       robust     heteroskedasticity-consistent standard errors
       names      names of the columns of y (default y1, y2, ...); exogNames those of exog (x1, ...)
       frequency, start   the calendar of the rows: observations per year (1, 4, 12, ...) and the first observation ("1990:1" or "1990Q1")
       fitStart, fitEnd   the estimation sample as observations of that calendar (default: all rows; the first `order` rows are used as presample) */
VarFit var(Y)(Y y, int order, GretlMatrix exog = GretlMatrix.init, bool constant = true, bool trend = false, bool seasonals = false,
              int[] lags = null, bool robust = false, string[] names = null, string[] exogNames = null, int frequency = 1,
              string start = "1", string fitStart = null, string fitEnd = null) @trusted if (isOperand!Y) {
    auto ym = cbind(y);
    const nx = exog.nrow == 0 ? 0 : exog.ncol;

    enforce!GretlException(ym.ncol >= 1, "var: y has no columns");
    requireFinite(ym, "var: y");
    if (nx > 0) {
        enforce!GretlException(exog.nrow == ym.nrow, "var: y has " ~ to!string(ym.nrow) ~ " rows but exog has " ~ to!string(exog.nrow));
        requireFinite(exog, "var: exog");
    }
    return estimateVar(columnsOf(ym), nx > 0 ? columnsOf(exog) : null, order, lags, constant, trend, seasonals, robust, namesOr(names, ym.ncol, "y", "var"),
                       namesOr(exogNames, nx, "x", "var"), frequency, start, fitStart, fitEnd);
}

/* var(d, ["y1", "y2"], order): a VAR on the named series of a Dataset, with `exog` naming the exogenous series; the calendar of the Dataset is used if it has one
   (frequency 0 and the default start take it from there), and the series must have no missing values in the sample (leading or trailing NaNs can be
   left out by choosing fitStart and fitEnd). */
VarFit var(Dataset d, string[] endog, int order, string[] exog = null, bool constant = true, bool trend = false, bool seasonals = false,
           int[] lags = null, bool robust = false, int frequency = 0, string start = null, string fitStart = null, string fitEnd = null) @trusted {
    const(double)[][] ycols, xcols;

    foreach (nm; endog) ycols ~= d[nm];
    foreach (nm; exog) xcols ~= d[nm];
    calendarOf(d, frequency, start);
    return estimateVar(ycols, xcols, order, lags, constant, trend, seasonals, robust, endog.dup, exog.dup, frequency, start, fitStart, fitEnd);
}

/* the calendar of a Dataset where the caller did not give one: frequency 0 and a null start mean "take it from the data" (annual from 1 if it has none) */
private void calendarOf(Dataset d, ref int frequency, ref string start) @trusted {
    auto raw = d.ptr;
    const dated = raw.structure == TsCodes.time_series;

    if (frequency == 0) frequency = dated ? raw.pd : 1;
    if (start is null) start = dated ? fromStringz(raw.stobs.ptr).idup : "1";
}

private VarFit estimateVar(const(double)[][] ycols, const(double)[][] xcols, int order, int[] lags, bool constant, bool trend, bool seasonals, bool robust,
                           string[] ynames, string[] xnames, int frequency, string start, string fitStart, string fitEnd) @trusted {
    const k = ycols.length, nx = xcols.length;

    enforce!GretlException(order >= 1, "var: the lag order must be at least 1");
    enforce!GretlException(seasonals ? frequency > 1 : true, "var: seasonal dummies need a frequency above 1");
    if (lags !is null) {
        enforce!GretlException(lags.length > 0, "var: empty lag list");
        foreach (l; lags) enforce!GretlException(l >= 1 && l <= order, "var: the lags must be between 1 and the order");
    }
    foreach (c; xcols) enforce!GretlException(c.length == ycols[0].length, "var: all series need the same number of observations");
    auto ds = timeDataset(ycols ~ xcols, frequency, start, "var", ynames ~ xnames);
    auto raw = ds.ptr;

    if (fitStart !is null) raw.t1 = observation(ds, fitStart, "var");
    if (fitEnd !is null) raw.t2 = observation(ds, fitEnd, "var");
    enforce!GretlException(raw.t1 < raw.t2, "var: the estimation sample is empty");

    int[] spec;

    foreach (j; 0 .. k) spec ~= cast(int) (1 + j);
    if (nx > 0) {
        spec ~= LISTSEP;
        foreach (j; 0 .. nx) spec ~= cast(int) (1 + k + j);
    }
    auto list = makeList(spec);
    int[] lagSpec = lags is null ? null : lags.dup;
    int * laglist = null;

    if (lagSpec !is null) laglist = makeList(lagSpec);
    auto opt = Gretlopt.OPT_Q;

    if (!constant) opt |= Gretlopt.OPT_N;
    if (trend) opt |= Gretlopt.OPT_T;
    if (seasonals) opt |= Gretlopt.OPT_D;
    if (robust) opt |= Gretlopt.OPT_R;
    enforce!GretlException(push_program_state() == 0, "could not save the library settings");
    scope (exit) {
        pop_program_state();
        import core.stdc.stdlib: free;
        free(list);
        if (laglist !is null) free(laglist);
    }
    int err = 0;
    auto v = gretl_VAR(order, laglist, list, raw, opt, null, &err);

    enforce!GretlException(err == 0 && v !is null, "var: the estimation failed, " ~ lastError());
    return collectVar(v, ds, ynames);
}

private VarFit collectVar(GretlVar * v, Dataset ds, string[] ynames) @trusted {
    VarFit r;

    r._p = SafeRefCounted!(VarFit.Payload, RefCountedAutoInitialize.no)(v, ds);
    r.names = ynames.dup;
    r.neqns = v.neqns;
    r.order = v.order;
    r.nobs = v.T;
    r.df = v.df;
    r.lnL = v.ll;
    r.aic = v.AIC;
    r.bic = v.BIC;
    r.hqc = v.HQC;
    r.ldet = v.ldet;
    if (v.lags !is null) {
        foreach (i; 0 .. v.lags[0]) r.lags ~= v.lags[i + 1];
    } else {
        foreach (i; 1 .. v.order + 1) r.lags ~= i;
    }
    const n = v.neqns, g = v.ncoeff;

    r.coef = copyOf(v.B);
    r.sigma = copyOf(v.S);
    r.resid = copyOf(v.E);
    r.se = fresh(g, n);
    foreach (j; 0 .. n) {
        const m = gretl_VAR_get_model(v, cast(int) j);

        foreach (i; 0 .. g) r.se.ptr.val[j * g + i] = m.sderr[i];
    }
    // coefficient names
    {
        import core.stdc.stdlib: calloc, free;
        auto params = cast(char **) calloc(g, (char *).sizeof);

        foreach (i; 0 .. g) params[i] = cast(char *) calloc(64, 1);
        gretl_VAR_param_names(v, params, ds.ptr);
        foreach (i; 0 .. g) {
            r.coefNames ~= fromStringz(params[i]).idup;
            free(params[i]);
        }
        free(params);
    }
    // roots of the companion matrix
    int err = 0;
    auto L = gretl_VAR_get_roots(v, &err);

    if (err == 0 && L !is null) {
        const nr = L.rows;

        r.roots = fresh(nr, 3);
        foreach (i; 0 .. nr) {
            const re = L.val[i], im = L.val[nr + i];

            r.roots.ptr.val[i] = re;
            r.roots.ptr.val[nr + i] = im;
            r.roots.ptr.val[2 * nr + i] = (re * re + im * im) ^^ 0.5;
        }
    }
    // Granger-type F tests: the library's Fvals hold, for each equation, n tests (all lags of variable j) and then one more
    const p = var_n_lags(v);

    r.grangerF = fresh(n, n);
    r.grangerP = fresh(n, n);
    if (v.Fvals !is null) {
        const dfd = gretl_VAR_get_model(v, 0).dfd;

        foreach (i; 0 .. n) {
            foreach (j; 0 .. n) {
                const f = v.Fvals[i * (n + 1) + j];

                r.grangerF.ptr.val[j * n + i] = f;
                r.grangerP.ptr.val[j * n + i] = isNaN(f) ? double.nan : snedecor_cdf_comp(p, dfd, f);
            }
        }
    }
    r.sampleStart = observationLabel(ds, v.t1);
    r.sampleEnd = observationLabel(ds, v.t2);
    return r;
}

private GretlMatrix copyOf(const(GretlMatrixRaw) * m) @trusted {
    auto r = fresh(m.rows, m.cols);

    foreach (i; 0 .. cast(size_t) m.rows * m.cols) r.ptr.val[i] = m.val[i];
    return r;
}

/* ---- lag selection ---- */

/* The result of lagSelection: one row per lag order, columns lag, loglik, pLR, aic, bic, hqc (see `columns`); the `best...` fields are the lag orders
   chosen by each criterion. */
struct LagSelection {
    GretlMatrix table;
    string[] columns = ["lags", "loglik", "p(LR)", "AIC", "BIC", "HQC"];
    int bestAic, bestBic, bestHqc;
}

/* lagSelection(y, maxLag): estimate VARs of order 1 .. maxLag on the same sample (the one of the longest VAR) and tabulate the log-likelihood, the p-value of the
   likelihood-ratio test of each lag against the one before, and the information criteria (Akaike, Schwarz-Bayesian, Hannan-Quinn: smaller is better). The options of `var` that
   describe the model (exog, constant, trend, seasonals, names, calendar and sample) apply. */
LagSelection lagSelection(Y)(Y y, int maxLag, GretlMatrix exog = GretlMatrix.init, bool constant = true, bool trend = false, bool seasonals = false,
                             int frequency = 1, string start = "1", string fitStart = null, string fitEnd = null) @trusted if (isOperand!Y) {
    auto ym = cbind(y);
    const nx = exog.nrow == 0 ? 0 : exog.ncol;

    enforce!GretlException(maxLag >= 2, "lagSelection: maxLag must be at least 2");
    requireFinite(ym, "lagSelection: y");
    if (nx > 0) {
        enforce!GretlException(exog.nrow == ym.nrow, "lagSelection: y and exog need the same number of rows");
        requireFinite(exog, "lagSelection: exog");
    }
    auto ds = timeDataset(columnsOf(ym) ~ (nx > 0 ? columnsOf(exog) : null), frequency, start, "lagSelection");
    auto raw = ds.ptr;

    if (fitStart !is null) raw.t1 = observation(ds, fitStart, "lagSelection");
    if (fitEnd !is null) raw.t2 = observation(ds, fitEnd, "lagSelection");
    int[] spec;

    foreach (j; 0 .. ym.ncol) spec ~= cast(int) (1 + j);
    if (nx > 0) {
        spec ~= LISTSEP;
        foreach (j; 0 .. nx) spec ~= cast(int) (1 + ym.ncol + j);
    }
    auto list = makeList(spec);
    auto opt = Gretlopt.OPT_S | Gretlopt.OPT_L;   // L: lag selection only

    if (!constant) opt |= Gretlopt.OPT_N;
    if (trend) opt |= Gretlopt.OPT_T;
    if (seasonals) opt |= Gretlopt.OPT_D;
    enforce!GretlException(push_program_state() == 0, "could not save the library settings");
    scope (exit) {
        pop_program_state();
        import core.stdc.stdlib: free;
        free(list);
    }
    GretlMatrixRaw * pm = null;
    const err = gui_VAR_lagsel(maxLag, list, raw, opt, null, &pm);

    enforce!GretlException(err == 0 && pm !is null, "lagSelection: it failed, " ~ lastError());
    scope (exit) gretl_matrix_free(pm);
    LagSelection ls;

    ls.table = copyOf(pm);
    const rows = pm.rows;

    foreach (c; 0 .. 3) {
        size_t best = 0;

        foreach (i; 1 .. rows) {
            if (ls.table[i, 3 + c] < ls.table[best, 3 + c]) best = i;
        }
        const lag = cast(int) ls.table[best, 0];

        if (c == 0) ls.bestAic = lag; else if (c == 1) ls.bestBic = lag; else ls.bestHqc = lag;
    }
    return ls;
}

/* ---- impulse responses, variance decompositions, forecasts ---- */

/* what irf, fevd and forecast need from a fit: the library's structure, its dataset and the number of equations (VarFit and VecmFit have them) */
enum hasVarHandle(F) = is(typeof(F.init.raw()) == GretlVar *) && is(typeof(F.init.neqns) : size_t);

/* The impulse responses of a VAR: response[targ * neqns + shock] is a periods x 1 matrix (the response of variable `targ` to a one-standard-deviation
   orthogonalized shock to variable `shock`, the Cholesky factor of sigma with the variables in the order of the VAR), or periods x 3 (point, lower, upper)
   when bootstrap bands were asked for. */
struct Irf {
    size_t neqns;
    int periods;
    double alpha;                   // 0: no bands; otherwise the bands cover 1 - alpha
    GretlMatrix[] response;

    GretlMatrix at(size_t targ, size_t shock) {
        enforce!GretlException(targ < neqns && shock < neqns, "Irf.at: variable out of range");
        return response[targ * neqns + shock];
    }
}

/* irf(fit, periods): the orthogonalized impulse responses (of a VarFit or a VecmFit). alpha in [0.01, 0.6] adds bootstrap bands (percentile intervals covering 1 - alpha, from the library's
   residual bootstrap); `seed` makes the bootstrap reproducible. periods = 0 takes the library's default horizon for the data frequency. */
Irf irf(F)(ref F fit, int periods = 0, double alpha = 0.0, long seed = -1) @trusted if (hasVarHandle!F) {
    enforce!GretlException(alpha == 0.0 || (alpha >= 0.01 && alpha <= 0.6), "irf: alpha must be 0 or between 0.01 and 0.6");
    enforce!GretlException(periods >= 0, "irf: periods must not be negative");
    auto raw = fit.raw;
    auto ds = fit.dataset;
    const h = periods == 0 ? default_VAR_horizon(ds.ptr) : periods;
    Irf r;

    r.neqns = fit.neqns;
    r.periods = h;
    r.alpha = alpha;
    if (seed >= 0) gretl_rand_set_seed(cast(ulong) seed);
    foreach (targ; 0 .. fit.neqns) {
        foreach (shock; 0 .. fit.neqns) {
            int err = 0;
            auto m = gretl_VAR_get_impulse_response(raw, cast(int) targ, cast(int) shock, h, alpha, ds.ptr, &err);

            enforce!GretlException(err == 0 && m !is null, "irf: failed, " ~ lastError());
            r.response ~= copyOf(m);
            gretl_matrix_free(m);
        }
    }
    return r;
}

/* The forecast error variance decomposition: share[targ * neqns + shock] is a horizon x 1 matrix with the fraction of the variance of the
   horizon-step forecast error of variable `targ` that comes from the orthogonalized shock to `shock` (the fractions over the shocks add to 1). */
struct Fevd {
    size_t neqns;
    int horizon;
    GretlMatrix[] share;

    GretlMatrix at(size_t targ, size_t shock) {
        enforce!GretlException(targ < neqns && shock < neqns, "Fevd.at: variable out of range");
        return share[targ * neqns + shock];
    }
}

Fevd fevd(F)(ref F fit, int horizon = 0) @trusted if (hasVarHandle!F) {
    enforce!GretlException(horizon >= 0, "fevd: horizon must not be negative");
    auto ds = fit.dataset;
    int err = 0;
    auto m = gretl_VAR_get_FEVD_matrix(fit.raw, -1, -1, horizon, ds.ptr, &err);

    enforce!GretlException(err == 0 && m !is null, "fevd: failed, " ~ lastError());
    scope (exit) gretl_matrix_free(m);
    Fevd r;

    r.neqns = fit.neqns;
    r.horizon = m.rows;
    foreach (c; 0 .. fit.neqns * fit.neqns) {
        auto col = fresh(m.rows, 1);

        foreach (i; 0 .. m.rows) col.ptr.val[i] = m.val[c * m.rows + i];
        r.share ~= col;
    }
    return r;
}

/* How a forecast over a range of observations is made: `automatic` (the library's rule) is one-step-ahead (static) forecasts from the actual data inside the estimation
   sample and dynamic forecasts after it; `oneStep` and `dynamic` force one or the other over the whole range (a dynamic forecast uses its own earlier
   forecasts from the start of the range on). */
enum ForecastMethod {
    automatic,
    oneStep,
    dynamic
}

private Gretlopt methodOpt(ForecastMethod m) {
    final switch (m) {
        case ForecastMethod.automatic: return Gretlopt.OPT_NONE;
        case ForecastMethod.oneStep: return Gretlopt.OPT_S;
        case ForecastMethod.dynamic: return Gretlopt.OPT_D;
    }
}

/* Puts the future values of the regressors into the dataset (extending it with missing values as far as `last`) for the time of a forecast; the cells that held actual
   data are restored afterwards, so that one forecast does not alter the data the next one sees. Row s of `values` goes to observation exogStart + s. */
private struct ExogGuard {
    DatasetRaw * dset;
    int[] ids;
    int t0;
    double[] saved;         // [column * rows + s], for the cells that existed
    size_t rows;
    int nBefore;

    @disable this(this);

    this(DatasetRaw * d, int[] regressorIds, GretlMatrix values, int exogStart, size_t rowCount, int last) @trusted {
        dset = d;
        ids = regressorIds;
        t0 = exogStart;
        rows = rowCount;
        nBefore = d.n;
        if (last >= d.n) {
            const err = dataset_add_observations(d, last + 1 - d.n, Gretlopt.OPT_A);

            enforce!GretlException(err == 0, "forecast: cannot extend the data, " ~ lastError());
        }
        saved = new double[ids.length * rows];
        foreach (k, id; ids) {
            foreach (s; 0 .. rows) {
                saved[k * rows + s] = (t0 + cast(int) s < nBefore) ? d.Z[id][t0 + s] : double.nan;
                d.Z[id][t0 + s] = values[s, k];
            }
        }
    }

    ~this() @trusted {
        foreach (k, id; ids) {
            foreach (s; 0 .. rows) {
                if (t0 + cast(int) s < nBefore) dset.Z[id][t0 + s] = saved[k * rows + s];
            }
        }
    }
}

/* the regressors a forecast needs: a check of `futureExog`, which has one row for each observation from exogStart on up to `last` */
private void checkExog(size_t needed, size_t have, size_t columns, GretlMatrix futureExog, int steps) @trusted {
    if (columns == 0) return;
    enforce!GretlException(needed == 0 || (futureExog.nrow >= needed && futureExog.ncol == columns),
                           "forecast: this model has " ~ to!string(columns) ~ " exogenous variables, so futureExog must have " ~ to!string(needed) ~ " rows and that many columns");
    if (needed > 0) requireFinite(futureExog, "forecast: futureExog");
}

/* The forecasts of a VAR or VECM: mean and standard error of each variable at each observation of the range (rows x neqns). */
struct VarForecast {
    int first;              // the observation (0-based) of row 0
    GretlMatrix mean;
    GretlMatrix se;
}

/* forecast(fit, steps): dynamic forecasts for the `steps` observations after the end of the estimation sample, with standard errors of the forecast errors
   (the ones of the library: they ignore the uncertainty of the coefficients). A VAR with exogenous regressors needs their future values in
   `futureExog` (steps x m). Seasonal dummies and the trend continue by themselves. */
VarForecast forecast(F)(ref F fit, int steps, GretlMatrix futureExog = GretlMatrix.init) @trusted if (hasVarHandle!F) {
    enforce!GretlException(steps >= 1, "forecast: steps must be at least 1");
    const t1 = fit.raw.t2 + 1;

    return varRange(fit.raw, fit.dataset, cast(int) fit.neqns, t1, t1 + steps - 1, ForecastMethod.automatic, futureExog, t1, steps);
}

/* forecastRange(fit, first, last): forecasts of a VAR or VECM for the observations first .. last (0-based rows of the data the model was estimated on, counting from
   the beginning of the series), inside the sample, after it, or across its end. By default (ForecastMethod.automatic) forecasts inside the estimation sample are
   one-step-ahead from the actual data and later ones are dynamic; `method` forces either over the whole range. The standard errors of one-step forecasts are those of the
   residuals. Observations beyond the end of the data (rows of the original matrix) need `futureExog`, if the model has exogenous variables: row 0 is the first
   observation after the data. The first `order` observations have no forecast (NaN). */
VarForecast forecastRange(F)(ref F fit, int first, int last, ForecastMethod method = ForecastMethod.automatic, GretlMatrix futureExog = GretlMatrix.init)
        @trusted if (hasVarHandle!F) {
    enforce!GretlException(first >= 0 && last >= first, "forecastRange: need 0 <= first <= last");
    const nOrig = fit.nOrig;

    return varRange(fit.raw, fit.dataset, cast(int) fit.neqns, first, last, method, futureExog, nOrig, last >= nOrig ? last - nOrig + 1 : 0);
}

private VarForecast varRange(GretlVar * raw, Dataset ds, int neqns, int t1, int t2, ForecastMethod method, GretlMatrix futureExog, int exogStart,
                             size_t exogRows) @trusted {
    auto dset = ds.ptr;
    int[] ids;

    if (raw.xlist !is null) {
        foreach (j; 1 .. raw.xlist[0] + 1) ids ~= raw.xlist[j];
    }
    checkExog(exogRows, 0, ids.length, futureExog, 0);
    auto guard = ExogGuard(dset, ids, futureExog, exogStart, ids.length > 0 ? exogRows : 0, t2);
    int err = 0;
    auto F = gretl_VAR_get_forecast_matrix(raw, t1, t2, dset, methodOpt(method), &err);

    enforce!GretlException(err == 0 && F !is null, "forecast: failed, " ~ lastError());
    const n = neqns, steps = t2 - t1 + 1;
    VarForecast f;

    f.first = t1;
    f.mean = fresh(steps, n);
    f.se = fresh(steps, n);
    foreach (j; 0 .. n) {
        foreach (s; 0 .. steps) {
            f.mean.ptr.val[j * steps + s] = F.val[j * steps + s];
            f.se.ptr.val[j * steps + s] = F.val[(n + j) * steps + s];
        }
    }
    return f;
}

/* ---- ARMA / ARIMA ---- */

/* How an ARMA model is estimated. */
enum ArmaMethod {
    exact,          // exact Gaussian maximum likelihood (the library's default: AS 154 / AS 197 routines)
    kalman,         // exact maximum likelihood by the Kalman filter (the library's `--kalman`)
    conditional     // conditional maximum likelihood (pre-sample innovations set to zero)
}

/* The estimate of an ARMA or ARIMA model: an Estimate (coefficients const, phi_1.., theta_1.., Phi_1.., Theta_1.., then the regressors; the residuals
   `resid` and `fitted` are those of the differenced model's one-step predictions) plus the orders and the roots of the polynomials. */
struct ArmaFit {
    private struct Payload {
        Model m;                    // the library's model, kept for the forecasts
        Dataset ds;                 // and the data it was estimated on
        int nOrig;                  // its number of observations before any forecast extended it

        this(ref Model src, Dataset d) {
            m = src;
            ds = d;
            nOrig = d.ptr.n;
        }

        ~this() {
            clear_model(&m);
        }

        @disable this(this);
    }

    private SafeRefCounted!(Payload, RefCountedAutoInitialize.no) _p;

    package(dgretl) Model * model() @trusted {
        return _p.refCountedStore.isInitialized ? &_p.refCountedPayload.m : null;
    }

    package(dgretl) Dataset dataset() @trusted {
        return _p.refCountedPayload.ds;
    }

    package(dgretl) int nOrig() @trusted {
        return _p.refCountedPayload.nOrig;
    }

    Estimate fit;
    alias fit this;
    int p, d, q;                // non-seasonal orders
    int sP, sD, sQ;             // seasonal orders
    GretlMatrix arRoots;        // roots of the AR polynomial (non-seasonal then seasonal): real part, imaginary part, modulus; outside the unit circle for stationarity
    GretlMatrix maRoots;        // roots of the MA polynomial in the same layout; outside the unit circle for invertibility
}

private string arimaLabel(int p, int d, int q, int sp, int sd, int sq, int pd) {
    string s = format("ARIMA(%d,%d,%d)", p, d, q);

    if (sp + sd + sq > 0) s ~= format("(%d,%d,%d)[%d]", sp, sd, sq, pd);
    return s;
}

/* arma(y, p, q): the model y(t) = const + phi_1 y(t-1) + .. + phi_p y(t-p) + e(t) + theta_1 e(t-1) + .. + theta_q e(t-q), or with d > 0 (or a seasonal difference) the ARIMA
   model for the differenced series.
       d          order of differencing (ARIMA)
       seasonal   [P, D, Q]: seasonal AR order, seasonal differences, seasonal MA order (needs frequency > 1)
       exog       regressors (ARMAX), n x m; exogNames their names (default x1, ..)
       constant   a constant (default); for an integrated model there is none to estimate and the regressors are differenced too
       method     ArmaMethod.exact (default), kalman or conditional
       opg        covariance matrix from the outer product of the gradient instead of the Hessian (the library gives no covariance matrix, and so
                  this throws, if the OPG matrix cannot be inverted or for a pure AR model on the exact route: use the default then)
       frequency, start, fitStart, fitEnd   as in `var`
   Missing values throw. The result is an ArmaFit (everything of an Estimate, through alias this, plus the orders and roots). */
ArmaFit arma(Y)(Y y, int p, int q, int d = 0, int[] seasonal = null, GretlMatrix exog = GretlMatrix.init, bool constant = true,
                ArmaMethod method = ArmaMethod.exact, bool opg = false, string[] exogNames = null, int frequency = 1, string start = "1",
                string fitStart = null, string fitEnd = null) @trusted if (isOperand!Y) {
    auto ym = cbind(y);
    const nx = exog.nrow == 0 ? 0 : exog.ncol;

    enforce!GretlException(ym.ncol == 1, "arma: y must have one column");
    requireFinite(ym, "arma: y");
    if (nx > 0) {
        enforce!GretlException(exog.nrow == ym.nrow, "arma: y and exog need the same number of rows");
        requireFinite(exog, "arma: exog");
    }
    return estimateArma(columnsOf(ym)[0], nx > 0 ? columnsOf(exog) : null, p, d, q, seasonal, constant, method, opg,
                        namesOr(exogNames, nx, "x", "arma"), frequency, start, fitStart, fitEnd);
}

/// ditto, for the series `y` (and `exog`) of a Dataset, using its calendar
ArmaFit arma(Dataset data, string y, int p, int q, int d = 0, int[] seasonal = null, string[] exog = null, bool constant = true,
             ArmaMethod method = ArmaMethod.exact, bool opg = false, int frequency = 0, string start = null, string fitStart = null,
             string fitEnd = null) @trusted {
    const(double)[][] xcols;

    foreach (nm; exog) xcols ~= data[nm];
    calendarOf(data, frequency, start);
    foreach (v; data[y]) enforce!GretlException(isFinite(v), "arma: " ~ y ~ " has missing values");
    return estimateArma(data[y], xcols, p, d, q, seasonal, constant, method, opg, exog.dup, frequency, start, fitStart, fitEnd);
}

private ArmaFit estimateArma(const(double)[] y, const(double)[][] xcols, int p, int d, int q, int[] seasonal, bool constant, ArmaMethod method, bool opg,
                             string[] xnames, int frequency, string start, string fitStart, string fitEnd) @trusted {
    int sP = 0, sD = 0, sQ = 0;

    if (seasonal !is null) {
        enforce!GretlException(seasonal.length == 3, "arma: seasonal is [P, D, Q]");
        sP = seasonal[0]; sD = seasonal[1]; sQ = seasonal[2];
    }
    enforce!GretlException(p >= 0 && q >= 0 && d >= 0 && sP >= 0 && sD >= 0 && sQ >= 0, "arma: orders must not be negative");
    enforce!GretlException(p + q + sP + sQ > 0, "arma: the model has no AR or MA terms");
    const seasonal_ = sP + sD + sQ > 0;

    enforce!GretlException(!seasonal_ || frequency > 1, "arma: a seasonal specification needs a frequency above 1");
    enforce!GretlException(y.length > cast(size_t) (p + q + d + frequency * (sP + sQ + sD) + 10), "arma: too few observations");
    auto ds = timeDataset([y] ~ xcols, frequency, start, "arma", ["y"] ~ xnames);
    auto raw = ds.ptr;

    if (fitStart !is null) raw.t1 = observation(ds, fitStart, "arma");
    if (fitEnd !is null) raw.t2 = observation(ds, fitEnd, "arma");
    const arima = d > 0 || sD > 0;
    int[] spec;

    spec ~= arima ? [p, d, q] : [p, q];
    spec ~= LISTSEP;
    if (seasonal_) spec ~= (arima ? [sP, sD, sQ] : [sP, sQ]) ~ [LISTSEP];
    spec ~= 1;
    if (xcols.length > 0 && constant) spec ~= 0;
    foreach (j; 0 .. xcols.length) spec ~= cast(int) (2 + j);
    auto list = makeList(spec);
    auto opt = Gretlopt.OPT_NONE;

    final switch (method) {
        case ArmaMethod.exact: break;
        case ArmaMethod.kalman: opt |= Gretlopt.OPT_K; break;
        case ArmaMethod.conditional: opt |= Gretlopt.OPT_C; break;
    }
    if (!constant) opt |= Gretlopt.OPT_N;
    if (opg) opt |= Gretlopt.OPT_G;
    enforce!GretlException(push_program_state() == 0, "could not save the library settings");
    scope (exit) pop_program_state();
    Model m = arma_model(list, null, raw, opt, null);
    bool keep = false;

    scope (exit) {
        if (!keep) clear_model(&m);
        import core.stdc.stdlib: free;
        free(list);
    }
    enforce!GretlException(m.errcode == 0, "arma: the estimation failed, " ~ errorText(m));
    ArmaFit r;
    const label = arimaLabel(p, d, q, sP, sD, sQ, frequency);

    r.fit = fromModel(m, label, "ML", ["const", "y"] ~ xnames, m.nobs);
    r.fit.names.length = 0;
    foreach (i; 0 .. cast(size_t) m.ncoeff) r.fit.names ~= fromStringz(m.params[i]).idup;
    r.p = p; r.d = d; r.q = q; r.sP = sP; r.sD = sD; r.sQ = sQ;
    r._p = SafeRefCounted!(ArmaFit.Payload, RefCountedAutoInitialize.no)(m, ds);
    keep = true;                // the payload owns the model now
    auto rts = cast(CephesCmplx *) gretl_model_get_data(&m, "roots");
    const nAr = p + sP, nMa = q + sQ;

    r.arRoots = fresh(nAr, 3);
    r.maRoots = fresh(nMa, 3);
    if (rts !is null) {
        foreach (kind; 0 .. 2) {
            const n = kind == 0 ? nAr : nMa, first = kind == 0 ? 0 : nAr;
            auto dst = kind == 0 ? r.arRoots : r.maRoots;

            foreach (i; 0 .. n) {
                const re = rts[first + i].r, im = rts[first + i].i;

                dst.ptr.val[i] = re;
                dst.ptr.val[n + i] = im;
                dst.ptr.val[2 * n + i] = (re * re + im * im) ^^ 0.5;
            }
        }
    }
    return r;
}

/* ---- GARCH ---- */

/* The estimate of a GARCH model: the mean equation and the variance parameters in one Estimate (names: const, the regressors, alpha(0), alpha(1) .. alpha(q),
   beta(1) .. beta(p)); `variance` is the conditional variance series (in the units of the data). */
struct GarchFit {
    private SafeRefCounted!(ArmaFit.Payload, RefCountedAutoInitialize.no) _p;     // the model and data, kept for the forecasts

    package(dgretl) Model * model() @trusted {
        return _p.refCountedStore.isInitialized ? &_p.refCountedPayload.m : null;
    }

    package(dgretl) Dataset dataset() @trusted {
        return _p.refCountedPayload.ds;
    }

    package(dgretl) int nOrig() @trusted {
        return _p.refCountedPayload.nOrig;
    }

    Estimate fit;
    alias fit this;
    int p, q;                   // GARCH (beta) order p, ARCH (alpha) order q
    GretlMatrix variance;       // nobs x 1: the conditional variance h(t)
    double lrStat = double.nan; // likelihood-ratio statistic of the GARCH terms against the plain regression (its p-value is not the usual chi-square: the null is on the boundary)
    int lrDf;
}

/* garch(y, p, q): y(t) = const + x(t) b + e(t) with e(t) | past ~ N(0, h(t)), h(t) = alpha0 + alpha_1 e(t-1)^2 + .. + alpha_q e(t-q)^2 + beta_1 h(t-1) + .. + beta_p h(t-p).
   p is the number of GARCH lags (of h) and q the number of ARCH lags (of e^2), as in the library: `garch(y, 1, 1)` is the usual GARCH(1,1). A GARCH(p, q) with p > 0 needs q > 0, and p + q is at most 7.
       exog       regressors of the mean equation (n x m GretlMatrix); constant: a constant (default)
       robust     quasi-ML (robust) covariance matrix
       fcp        use the Fiorentini-Calzolari-Panattoni algorithm (the default is the library's BFGS-based estimator)
       standardized   residuals divided by the conditional standard deviation
       armaInit   start the variance parameters from an ARMA fit of the squared residuals */
GarchFit garch(Y)(Y y, int p, int q, GretlMatrix exog = GretlMatrix.init, bool constant = true, bool robust = false, bool fcp = false,
                  bool standardized = false, bool armaInit = false, string[] exogNames = null) @trusted if (isOperand!Y) {
    auto ym = cbind(y);
    const nx = exog.nrow == 0 ? 0 : exog.ncol;

    enforce!GretlException(ym.ncol == 1, "garch: y must have one column");
    requireFinite(ym, "garch: y");
    if (nx > 0) {
        enforce!GretlException(exog.nrow == ym.nrow, "garch: y and exog need the same number of rows");
        requireFinite(exog, "garch: exog");
    }
    return estimateGarch(columnsOf(ym)[0], nx > 0 ? columnsOf(exog) : null, p, q, constant, robust, fcp, standardized, armaInit,
                         namesOr(exogNames, nx, "x", "garch"));
}

/// ditto, for the series `y` (and `exog`) of a Dataset
GarchFit garch(Dataset data, string y, int p, int q, string[] exog = null, bool constant = true, bool robust = false, bool fcp = false,
               bool standardized = false, bool armaInit = false) @trusted {
    const(double)[][] xcols;

    foreach (nm; exog) xcols ~= data[nm];
    foreach (v; data[y]) enforce!GretlException(isFinite(v), "garch: " ~ y ~ " has missing values");
    return estimateGarch(data[y], xcols, p, q, constant, robust, fcp, standardized, armaInit, exog.dup);
}

private GarchFit estimateGarch(const(double)[] y, const(double)[][] xcols, int p, int q, bool constant, bool robust, bool fcp, bool standardized, bool armaInit,
                               string[] xnames) @trusted {
    enforce!GretlException(p >= 0 && q >= 1, "garch: q must be at least 1 (and p at least 0)");
    enforce!GretlException(p + q <= 7, "garch: p + q must not exceed 7");
    enforce!GretlException(y.length > cast(size_t) (p + q + xcols.length + 20), "garch: too few observations");
    auto ds = timeDataset([y] ~ xcols, 1, "1", "garch", ["y"] ~ xnames);
    int[] spec = [p, q, LISTSEP, 1];

    foreach (j; 0 .. xcols.length) spec ~= cast(int) (2 + j);
    auto list = makeList(spec);
    auto opt = Gretlopt.OPT_NONE;

    if (!constant) opt |= Gretlopt.OPT_N;
    if (robust) opt |= Gretlopt.OPT_R;
    if (fcp) opt |= Gretlopt.OPT_F;
    if (standardized) opt |= Gretlopt.OPT_Z;
    if (armaInit) opt |= Gretlopt.OPT_A;
    enforce!GretlException(push_program_state() == 0, "could not save the library settings");
    scope (exit) pop_program_state();
    Model m = garch_model(list, ds.ptr, null, opt);
    bool keep = false;

    scope (exit) {
        if (!keep) clear_model(&m);
        import core.stdc.stdlib: free;
        free(list);
    }
    enforce!GretlException(m.errcode == 0, "garch: the estimation failed, " ~ errorText(m));
    GarchFit r;

    r._p = SafeRefCounted!(ArmaFit.Payload, RefCountedAutoInitialize.no)(m, ds);
    keep = true;                // the payload owns the model now

    r.fit = fromModel(m, "GARCH", robust ? "QML (robust)" : "ML", ["const", "y"] ~ xnames, m.nobs);
    r.fit.names.length = 0;
    foreach (i; 0 .. cast(size_t) m.ncoeff) r.fit.names ~= fromStringz(m.params[i]).idup;
    r.p = p;
    r.q = q;
    auto h = cast(const(double) *) gretl_model_get_data(&m, "garch_h");

    if (h !is null) r.variance = columnOf(h, m.t1, m.nobs);
    r.lrStat = gretl_model_get_double(&m, "garch_LR");
    r.lrDf = gretl_model_get_int(&m, "garch_LR_df");
    return r;
}

/* ---- cointegration: the Johansen test and the VECM ---- */

/* The deterministic terms of a cointegration model (Johansen's five cases). */
enum Deterministic {
    none,                   // no constant
    restrictedConstant,     // a constant restricted to the cointegrating space
    unrestrictedConstant,   // an unrestricted constant (default): drift in the levels
    restrictedTrend,        // an unrestricted constant and a trend restricted to the cointegrating space
    unrestrictedTrend       // an unrestricted constant and trend
}

private Gretlopt deterministicOpt(Deterministic d) {
    final switch (d) {
        case Deterministic.none: return Gretlopt.OPT_N;
        case Deterministic.restrictedConstant: return Gretlopt.OPT_R;
        case Deterministic.unrestrictedConstant: return Gretlopt.OPT_NONE;
        case Deterministic.restrictedTrend: return Gretlopt.OPT_A;
        case Deterministic.unrestrictedTrend: return Gretlopt.OPT_T;
    }
}

/* The Johansen cointegration test: row r of each n x 1 matrix is for the null hypothesis that the cointegration rank is at most r. */
struct JohansenTest {
    size_t neqns;
    size_t nobs;
    int order;                  // lag order of the VAR in levels
    GretlMatrix eigenvalues;    // n x 1, decreasing
    GretlMatrix trace;          // trace statistics
    GretlMatrix lmax;           // maximum eigenvalue statistics
    GretlMatrix tracePvalue;    // p-values of the trace tests (corrected for the sample size, unless `asymptotic`; the library's gamma approximation)
    GretlMatrix lmaxPvalue;     // asymptotic p-values of the maximum eigenvalue tests
    GretlMatrix beta;           // n(+) x n: eigenvectors, normalized so that beta' S11 beta = I (the restricted terms add rows)
    GretlMatrix alpha;          // n x n: the matching adjustment coefficients
    double lnL = double.nan;
    string report;              // the library's printed table
}

/* johansen(y, order): the Johansen trace and maximum eigenvalue tests on the columns of y (n x k), for a VAR in levels of lag order `order` (>= 1: order - 1
   lagged differences in the VECM form).
       deterministic   which of the five cases (Deterministic.unrestrictedConstant by default)
       exog            unrestricted exogenous regressors (n x m); seasonals: centered seasonal dummies (frequency > 1)
       asymptotic      asymptotic p-values for the trace test
       names, frequency, start, fitStart, fitEnd   as in `var` */
JohansenTest johansen(Y)(Y y, int order, Deterministic deterministic = Deterministic.unrestrictedConstant, GretlMatrix exog = GretlMatrix.init,
                         bool seasonals = false, bool asymptotic = false, string[] names = null, string[] exogNames = null, int frequency = 1,
                         string start = "1", string fitStart = null, string fitEnd = null) @trusted if (isOperand!Y) {
    auto ym = cbind(y);
    const nx = exog.nrow == 0 ? 0 : exog.ncol;

    enforce!GretlException(ym.ncol >= 2, "johansen: need at least two series");
    requireFinite(ym, "johansen: y");
    if (nx > 0) {
        enforce!GretlException(exog.nrow == ym.nrow, "johansen: y and exog need the same number of rows");
        requireFinite(exog, "johansen: exog");
    }
    return runJohansen(columnsOf(ym), nx > 0 ? columnsOf(exog) : null, order, deterministic, seasonals, asymptotic, namesOr(names, ym.ncol, "y", "johansen"),
                       namesOr(exogNames, nx, "x", "johansen"), frequency, start, fitStart, fitEnd);
}

/// ditto, for the named series of a Dataset (using its calendar)
JohansenTest johansen(Dataset d, string[] endog, int order, Deterministic deterministic = Deterministic.unrestrictedConstant, string[] exog = null,
                      bool seasonals = false, bool asymptotic = false, int frequency = 0, string start = null, string fitStart = null,
                      string fitEnd = null) @trusted {
    const(double)[][] ycols, xcols;

    foreach (nm; endog) ycols ~= d[nm];
    foreach (nm; exog) xcols ~= d[nm];
    calendarOf(d, frequency, start);
    return runJohansen(ycols, xcols, order, deterministic, seasonals, asymptotic, endog.dup, exog.dup, frequency, start, fitStart, fitEnd);
}

/* a dataset with the series y then x, the list for the library, and the common checks of the cointegration functions */
private struct CointInput {
    Dataset ds;
    int * list;
    string[] names;
}

private CointInput cointInput(const(double)[][] ycols, const(double)[][] xcols, int order, bool seasonals, string[] ynames, string[] xnames, int frequency,
                              string start, string fitStart, string fitEnd, string who) @trusted {
    enforce!GretlException(order >= 1, who ~ ": the lag order must be at least 1");
    enforce!GretlException(!seasonals || frequency > 1, who ~ ": seasonal dummies need a frequency above 1");
    foreach (c; xcols) enforce!GretlException(c.length == ycols[0].length, who ~ ": all series need the same number of observations");
    CointInput ci;

    ci.names = ynames ~ xnames;
    ci.ds = timeDataset(ycols ~ xcols, frequency, start, who, ci.names);
    auto raw = ci.ds.ptr;

    if (fitStart !is null) raw.t1 = observation(ci.ds, fitStart, who);
    if (fitEnd !is null) raw.t2 = observation(ci.ds, fitEnd, who);
    enforce!GretlException(raw.t1 < raw.t2, who ~ ": the estimation sample is empty");
    int[] spec;

    foreach (j; 0 .. ycols.length) spec ~= cast(int) (1 + j);
    if (xcols.length > 0) {
        spec ~= LISTSEP;
        foreach (j; 0 .. xcols.length) spec ~= cast(int) (1 + ycols.length + j);
    }
    ci.list = makeList(spec);
    return ci;
}

private JohansenTest runJohansen(const(double)[][] ycols, const(double)[][] xcols, int order, Deterministic det, bool seasonals, bool asymptotic,
                                 string[] ynames, string[] xnames, int frequency, string start, string fitStart, string fitEnd) @trusted {
    auto ci = cointInput(ycols, xcols, order, seasonals, ynames, xnames, frequency, start, fitStart, fitEnd, "johansen");
    auto opt = deterministicOpt(det);

    if (seasonals) opt |= Gretlopt.OPT_D;
    if (asymptotic) opt |= Gretlopt.OPT_Y;
    auto pr = Printer(0);

    enforce!GretlException(push_program_state() == 0, "could not save the library settings");
    scope (exit) {
        pop_program_state();
        import core.stdc.stdlib: free;
        free(ci.list);
    }
    auto v = johansen_test(order, ci.list, ci.ds.ptr, opt, pr.prn);

    scope (exit) if (v !is null) gretl_VAR_free(v);
    enforce!GretlException(v !is null && v.err == 0, "johansen: the test failed, " ~ lastError());
    JohansenTest r;
    const n = v.neqns;

    r.neqns = n;
    r.nobs = v.T;
    r.order = order;
    r.lnL = v.ll;
    r.report = pr.text;
    r.eigenvalues = copyOf(v.jinfo.evals);
    r.beta = copyOf(v.jinfo.Beta);
    r.alpha = copyOf(v.jinfo.Alpha);
    int err = 0;
    auto tests = get_last_test_matrix(&err);
    auto pvals = get_last_pvals_matrix(&err);

    enforce!GretlException(err == 0 && tests !is null && pvals !is null, "johansen: the test statistics were not recorded");
    scope (exit) {
        gretl_matrix_free(tests);
        gretl_matrix_free(pvals);
    }
    r.trace = fresh(n, 1);
    r.lmax = fresh(n, 1);
    r.tracePvalue = fresh(n, 1);
    r.lmaxPvalue = fresh(n, 1);
    foreach (i; 0 .. n) {
        r.trace.ptr.val[i] = tests.val[i];
        r.lmax.ptr.val[i] = tests.val[n + i];
        r.tracePvalue.ptr.val[i] = pvals.val[i];
        r.lmaxPvalue.ptr.val[i] = pvals.val[n + i];
    }
    return r;
}

/* The estimate of a vector error correction model
       dy(t) = alpha beta' y(t-1) + Gamma_1 dy(t-1) + .. + Gamma_(p-1) dy(t-p+1) + deterministic terms + e(t).
   Like VarFit it keeps the library's structure, so irf, fevd and forecast work on it. */
struct VecmFit {
    private SafeRefCounted!(VarFit.Payload, RefCountedAutoInitialize.no) _p;

    string[] names;             // the endogenous variables
    string[] coefNames;         // the rows of `coef`: const, d_y1_1 .., the exogenous variables, EC1 ..
    string[] betaNames;         // the rows of `beta`: the variables, then "const" or "trend" when restricted to the cointegrating space
    size_t neqns;
    size_t rank;
    size_t order;               // the lag order of the VAR in levels
    size_t nobs;
    int df;
    GretlMatrix beta;           // cointegrating vectors, (neqns or more) x rank; normalized (Phillips) unless restricted
    GretlMatrix alpha;          // adjustment coefficients, neqns x rank
    GretlMatrix betaSe;         // standard errors of beta (0 for the normalizing elements), same shape; NaN if the library gives none
    GretlMatrix alphaSe;
    GretlMatrix pi;             // alpha * beta', the long-run matrix
    GretlMatrix coef;           // the equations' coefficients, ncoef x neqns, given beta
    GretlMatrix se;             // their standard errors
    GretlMatrix sigma;          // covariance matrix of the residuals
    GretlMatrix resid;          // nobs x neqns
    GretlMatrix eigenvalues;    // the eigenvalues of the Johansen problem, all n (empty for a model with a restricted beta)
    double lnL = double.nan;
    double aic = double.nan, bic = double.nan, hqc = double.nan;
    double restrictionLR = double.nan;  // for a restricted beta or alpha: the LR statistic against the unrestricted model, and its degrees of freedom
    int restrictionDf;
    string sampleStart, sampleEnd;

    package(dgretl) GretlVar * raw() @trusted {
        return _p.refCountedStore.isInitialized ? _p.refCountedPayload.raw : null;
    }

    package(dgretl) Dataset dataset() @trusted {
        return _p.refCountedPayload.ds;
    }

    package(dgretl) int nOrig() @trusted {
        return _p.refCountedPayload.nOrig;
    }

    /* equation i (0-based) of the error correction form as an Estimate (the dependent variable is the difference); the library keeps no covariance
       matrix for these equations, `vcov` is empty */
    Estimate equation(size_t i) @trusted {
        enforce!GretlException(i < neqns, "VecmFit.equation: there are only " ~ to!string(neqns) ~ " equations");
        auto m = *cast(Model *) gretl_VAR_get_model(raw, cast(int) i);
        auto e = fromModel(m, "VECM equation for d_" ~ names[i], "classical", ["const"] ~ names, nobs, false);

        e.names = coefNames.dup;
        return e;
    }

    string toString() {
        string s = format("VECM of %s, lag order %d, rank %d, %d observations (%s to %s)\n", names, order, rank, nobs, sampleStart, sampleEnd);

        s ~= format("log-likelihood %.4f, AIC %.5f, BIC %.5f, HQC %.5f\n\nbeta:\n", lnL, aic, bic, hqc);
        foreach (i, nm; betaNames) {
            s ~= format("%-12s", nm);
            foreach (j; 0 .. rank) s ~= format(" %12.6g (%.4g)", beta[i, j], betaSe[i, j]);
            s ~= "\n";
        }
        s ~= "\nalpha:\n";
        foreach (i, nm; names) {
            s ~= format("%-12s", nm);
            foreach (j; 0 .. rank) s ~= format(" %12.6g (%.4g)", alpha[i, j], alphaSe[i, j]);
            s ~= "\n";
        }
        return s;
    }
}

/* A homogeneous linear restriction on the cointegrating vectors: R beta = 0, with R (g x rows of beta; q, g x 1, may be left out and must be zero if given). It applies
   to every cointegrating vector (every column of beta), whatever the rank: R = [1 1] says that each is proportional to (1, -1). Restrictions with q different from 0 are
   refused (the library returns NaN for the test of one and fails to estimate). */
struct BetaRestriction {
    GretlMatrix R;
    GretlMatrix q;
}

/* A homogeneous linear restriction on the adjustment coefficients: R alpha = 0, with R (g x neqns; q, g x 1, may be left out and must be zero if given). It applies to
   every column of alpha: R = [1 0] says that the first variable does not adjust (weak exogeneity), R = [1 0; 0 1] that neither of the first two does. */
struct AlphaRestriction {
    GretlMatrix R;
    GretlMatrix q;
}

/* the library's restriction structure for at most one of a beta and an alpha restriction, with its matrices allocated (release with freeRestriction) */
private GretlRestriction makeRestriction(BetaRestriction b, AlphaRestriction a, size_t betaRows, size_t neqns, int rank, string who) @trusted {
    const hasB = b.R.nrow > 0, hasA = a.R.nrow > 0;
    GretlRestriction rs;

    enforce!GretlException(!(hasB && hasA), who ~ ": restrictions on beta and alpha together are not supported (give one of them)");
    if (!hasB && !hasA) return rs;
    auto R = hasB ? b.R : a.R;
    auto q = hasB ? b.q : a.q;
    const cols = hasB ? betaRows : neqns;
    const side = hasB ? "beta" : "alpha";

    enforce!GretlException(R.ncol == cols, who ~ ": R for " ~ side ~ " must have " ~ to!string(cols) ~ " columns, it has " ~ to!string(R.ncol));
    enforce!GretlException(R.nrow < cols, who ~ ": R for " ~ side ~ " has as many rows as columns: that fixes " ~ side ~ " completely");
    enforce!GretlException(q.nrow == 0 || (q.nrow == R.nrow && q.ncol == 1), who ~ ": q must have one row for each row of R");
    requireFinite(R, who ~ ": R");
    bool homogeneous = true;

    foreach (i; 0 .. q.nrow) {
        if (q[i, 0] != 0) homogeneous = false;
    }
    enforce!GretlException(homogeneous, who ~ ": only homogeneous restrictions (q = 0) are supported: the library's tests and estimation of the non-homogeneous case do not work");
    auto qm = q.nrow == 0 ? zerosMatrix(R.nrow, 1) : q;
    auto rraw = copyRaw(R), qraw = copyRaw(qm);

    if (hasB) {
        rs.R = rraw;
        rs.q = qraw;
    } else {
        rs.Ra = rraw;
        rs.qa = qraw;
    }
    rs.g = cast(int) R.nrow;
    rs.gmax = cast(int) cols;
    return rs;
}

private void freeRestriction(ref GretlRestriction rs) @trusted {
    gretl_matrix_free(rs.R);
    gretl_matrix_free(rs.q);
    gretl_matrix_free(rs.Ra);
    gretl_matrix_free(rs.qa);
}

private GretlMatrix zerosMatrix(size_t r, size_t c) @trusted {
    auto m = fresh(r, c);

    foreach (i; 0 .. r * c) m.ptr.val[i] = 0;
    return m;
}

/* vecm(y, order, rank): estimate the VECM of the columns of y (n x k) with lag order `order` in levels and cointegration rank `rank` (1 .. k - 1), by Johansen's
   maximum likelihood. Options as in `johansen`, plus
       restriction       a BetaRestriction on the cointegrating vectors: beta is estimated subject to it, and `restrictionLR` is the test of it
       alphaRestriction  an AlphaRestriction on the adjustment coefficients (not together with a restriction on beta) */
VecmFit vecm(Y)(Y y, int order, int rank, Deterministic deterministic = Deterministic.unrestrictedConstant, GretlMatrix exog = GretlMatrix.init,
                bool seasonals = false, BetaRestriction restriction = BetaRestriction.init, AlphaRestriction alphaRestriction = AlphaRestriction.init, string[] names = null,
                string[] exogNames = null,
                int frequency = 1, string start = "1", string fitStart = null, string fitEnd = null) @trusted if (isOperand!Y) {
    auto ym = cbind(y);
    const nx = exog.nrow == 0 ? 0 : exog.ncol;

    enforce!GretlException(ym.ncol >= 2, "vecm: need at least two series");
    requireFinite(ym, "vecm: y");
    if (nx > 0) {
        enforce!GretlException(exog.nrow == ym.nrow, "vecm: y and exog need the same number of rows");
        requireFinite(exog, "vecm: exog");
    }
    return estimateVecm(columnsOf(ym), nx > 0 ? columnsOf(exog) : null, order, rank, deterministic, seasonals, restriction, alphaRestriction, namesOr(names, ym.ncol, "y", "vecm"),
                        namesOr(exogNames, nx, "x", "vecm"), frequency, start, fitStart, fitEnd);
}

/// ditto, for the named series of a Dataset (using its calendar)
VecmFit vecm(Dataset d, string[] endog, int order, int rank, Deterministic deterministic = Deterministic.unrestrictedConstant, string[] exog = null,
             bool seasonals = false, BetaRestriction restriction = BetaRestriction.init, AlphaRestriction alphaRestriction = AlphaRestriction.init, int frequency = 0,
             string start = null, string fitStart = null,
             string fitEnd = null) @trusted {
    const(double)[][] ycols, xcols;

    foreach (nm; endog) ycols ~= d[nm];
    foreach (nm; exog) xcols ~= d[nm];
    calendarOf(d, frequency, start);
    return estimateVecm(ycols, xcols, order, rank, deterministic, seasonals, restriction, alphaRestriction, endog.dup, exog.dup, frequency, start, fitStart, fitEnd);
}

private VecmFit estimateVecm(const(double)[][] ycols, const(double)[][] xcols, int order, int rank, Deterministic det, bool seasonals, BetaRestriction restriction,
                             AlphaRestriction alphaRestriction, string[] ynames, string[] xnames, int frequency, string start, string fitStart, string fitEnd) @trusted {
    const n = ycols.length;

    enforce!GretlException(rank >= 1 && rank < cast(int) n, "vecm: the rank must be between 1 and " ~ to!string(n - 1));
    auto ci = cointInput(ycols, xcols, order, seasonals, ynames, xnames, frequency, start, fitStart, fitEnd, "vecm");
    auto opt = deterministicOpt(det) | Gretlopt.OPT_Q;

    if (seasonals) opt |= Gretlopt.OPT_D;
    const betaRows = n + (det == Deterministic.restrictedConstant || det == Deterministic.restrictedTrend ? 1 : 0);
    auto rs = makeRestriction(restriction, alphaRestriction, betaRows, n, rank, "vecm");
    const restricted = restriction.R.nrow > 0 || alphaRestriction.R.nrow > 0;

    enforce!GretlException(push_program_state() == 0, "could not save the library settings");
    scope (exit) {
        pop_program_state();
        import core.stdc.stdlib: free;
        free(ci.list);
        freeRestriction(rs);
    }
    int err = 0;
    GretlVar * v;

    if (restricted) {
        v = johansen_wrapper(VarTask.vecm_estimate, order, rank, null, ci.list, &rs, ci.ds.ptr, opt, null, &err);
    } else {
        v = johansen_wrapper(VarTask.vecm_estimate, order, rank, null, ci.list, null, ci.ds.ptr, opt, null, &err);
    }
    if (v is null || err != 0 || v.err != 0) {
        if (v !is null) gretl_VAR_free(v);
        throw new GretlException("vecm: the estimation failed, " ~ lastError());
    }
    VecmFit r;

    r._p = SafeRefCounted!(VarFit.Payload, RefCountedAutoInitialize.no)(v, ci.ds);
    r.names = ynames.dup;
    r.neqns = v.neqns;
    r.rank = rank;
    r.order = order;
    r.nobs = v.T;
    r.df = v.df;
    r.lnL = v.ll;
    r.aic = v.AIC;
    r.bic = v.BIC;
    r.hqc = v.HQC;
    auto jv = v.jinfo;

    r.beta = copyOf(jv.Beta);
    r.alpha = copyOf(jv.Alpha);
    r.eigenvalues = jv.evals !is null ? copyOf(jv.evals) : GretlMatrix.init;   // the library keeps none for a restricted model
    r.betaSe = jv.Bse !is null ? copyOf(jv.Bse) : nanLike(r.beta);
    r.alphaSe = jv.Ase !is null ? copyOf(jv.Ase) : nanLike(r.alpha);
    r.pi = r.alpha * t(r.beta);
    r.sigma = copyOf(v.S);
    r.resid = copyOf(v.E);
    r.betaNames = ynames.dup;
    if (r.beta.nrow > n) {
        if (det == Deterministic.restrictedConstant) r.betaNames ~= "const";
        else if (det == Deterministic.restrictedTrend) r.betaNames ~= "trend";
    }
    while (r.betaNames.length < r.beta.nrow) r.betaNames ~= "x" ~ to!string(r.betaNames.length + 1);
    if (restricted) {
        r.restrictionLR = 2 * (jv.ll0 - v.ll);
        r.restrictionDf = jv.lrdf;
    }
    const g = v.models[0].ncoeff;

    r.coef = fresh(g, n);
    r.se = fresh(g, n);
    foreach (j; 0 .. n) {
        const m = gretl_VAR_get_model(v, cast(int) j);

        foreach (i; 0 .. g) {
            r.coef.ptr.val[j * g + i] = m.coeff[i];
            r.se.ptr.val[j * g + i] = m.sderr[i];
        }
    }
    {
        import core.stdc.stdlib: calloc, free;
        auto params = cast(char **) calloc(v.ncoeff, (char *).sizeof);

        foreach (i; 0 .. cast(size_t) v.ncoeff) params[i] = cast(char *) calloc(64, 1);
        gretl_VAR_param_names(v, params, ci.ds.ptr);
        foreach (i; 0 .. cast(size_t) v.ncoeff) {
            r.coefNames ~= fromStringz(params[i]).idup;
            free(params[i]);
        }
        free(params);
    }
    r.sampleStart = observationLabel(ci.ds, v.t1);
    r.sampleEnd = observationLabel(ci.ds, v.t2);
    return r;
}

private GretlMatrix nanLike(GretlMatrix m) @trusted {
    auto r = fresh(m.nrow, m.ncol);

    foreach (i; 0 .. m.nrow * m.ncol) r.ptr.val[i] = double.nan;
    return r;
}

private GretlMatrixRaw * copyRaw(GretlMatrix m) @trusted {
    auto r = gretl_matrix_alloc(cast(int) m.nrow, cast(int) m.ncol);

    foreach (j; 0 .. m.ncol) {
        foreach (i; 0 .. m.nrow) r.val[j * m.nrow + i] = m[i, j];
    }
    return r;
}

/* testBeta(fit, restriction): the likelihood-ratio test of a linear restriction R beta = q on the cointegrating vectors (chi-square with rank * rows(R) degrees of
   freedom for a restriction that applies to every column of beta). The fit is left as it was. */
TestResult testBeta(ref VecmFit fit, BetaRestriction restriction) @trusted {
    enforce!GretlException(restriction.R.nrow > 0, "testBeta: no restriction given");
    auto rs = makeRestriction(restriction, AlphaRestriction.init, fit.beta.nrow, fit.neqns, cast(int) fit.rank, "testBeta");

    return runRestrictionTest(fit, rs, "testBeta");
}

/* testAlpha(fit, restriction): the likelihood-ratio test of a linear restriction R alpha = q on the adjustment coefficients (chi-square with rank * rows(R) degrees of
   freedom for a restriction that applies to every column of alpha); for example R = [1 0] tests that the first variable is weakly exogenous. The fit is left as it was. */
TestResult testAlpha(ref VecmFit fit, AlphaRestriction restriction) @trusted {
    enforce!GretlException(restriction.R.nrow > 0, "testAlpha: no restriction given");
    auto rs = makeRestriction(BetaRestriction.init, restriction, fit.beta.nrow, fit.neqns, cast(int) fit.rank, "testAlpha");

    return runRestrictionTest(fit, rs, "testAlpha");
}

private TestResult runRestrictionTest(ref VecmFit fit, GretlRestriction rs, string who) @trusted {
    scope (exit) freeRestriction(rs);
    auto pr = Printer(0);
    auto ds = fit.dataset;
    const err = vecm_test_restriction(fit.raw, &rs, ds.ptr, Gretlopt.OPT_NONE, pr.prn);

    enforce!GretlException(err == 0, who ~ ": the test failed, " ~ lastError());
    TestResult t;

    t.statistic = rs.test;
    t.pvalue = rs.pval;
    t.df = rs.g * cast(int) fit.rank;
    return t;
}

/* ---- ARMA and GARCH forecasts ---- */

/* The forecasts of an ARMA or GARCH model: mean and standard error of the forecast error at each observation (rows x 1 each). */
struct SeriesForecast {
    int first;                  // the observation (0-based) of row 0
    GretlMatrix mean;
    GretlMatrix se;
}

alias ArmaForecast = SeriesForecast;
alias GarchForecast = SeriesForecast;

/* forecast(fit, steps): dynamic forecasts of the series for the `steps` observations after the end of the estimation sample, with the standard errors of the
   forecast errors (from the MA(infinity) weights of the model; they ignore the uncertainty of the coefficients). An ARIMA model gives forecasts of the level
   of the series. A model with regressors (ARMAX) needs their future values in `futureExog` (steps x m, in the order of `exog`). */
SeriesForecast forecast(ref ArmaFit fit, int steps, GretlMatrix futureExog = GretlMatrix.init) @trusted {
    enforce!GretlException(steps >= 1, "forecast: steps must be at least 1");
    const t1 = fit.model.t2 + 1;

    return modelForecast(fit.model, fit.dataset, t1, t1 + steps - 1, ForecastMethod.automatic, futureExog, t1, steps);
}

/* forecast(garchFit, steps): forecasts of the mean equation of a GARCH model for the `steps` observations after the sample (the regressors' future values in
   `futureExog`, steps x m, if the model has any), with `se` the forecast of the conditional standard deviation sqrt(h) at each step (the library's: it is the
   standard deviation of the error of the mean forecast given the GARCH process, ignoring parameter uncertainty). */
SeriesForecast forecast(ref GarchFit fit, int steps, GretlMatrix futureExog = GretlMatrix.init) @trusted {
    enforce!GretlException(steps >= 1, "forecast: steps must be at least 1");
    const t1 = fit.model.t2 + 1;

    return modelForecast(fit.model, fit.dataset, t1, t1 + steps - 1, ForecastMethod.automatic, futureExog, t1, steps);
}

/* forecastRange(fit, first, last): forecasts of an ARMA or GARCH model for the observations first .. last (0-based rows of the data the model was estimated
   on), inside the estimation sample, after it, or across its end. With ForecastMethod.automatic the forecasts inside the sample are the one-step-ahead (fitted)
   values, the later ones dynamic; `oneStep` and `dynamic` force one or the other (a dynamic forecast started inside the sample uses the model's own forecasts
   from `first` on). One-step forecasts of an ARMA model have no standard error (NaN); for a GARCH model it is the conditional standard deviation. Observations
   after the end of the data need `futureExog` if there are regressors: row 0 is the first observation after the data. */
SeriesForecast forecastRange(ref ArmaFit fit, int first, int last, ForecastMethod method = ForecastMethod.automatic, GretlMatrix futureExog = GretlMatrix.init)
        @trusted {
    enforce!GretlException(first >= 0 && last >= first, "forecastRange: need 0 <= first <= last");
    const nOrig = fit.nOrig;

    return modelForecast(fit.model, fit.dataset, first, last, method, futureExog, nOrig, last >= nOrig ? last - nOrig + 1 : 0);
}

/// ditto
SeriesForecast forecastRange(ref GarchFit fit, int first, int last, ForecastMethod method = ForecastMethod.automatic, GretlMatrix futureExog = GretlMatrix.init)
        @trusted {
    enforce!GretlException(first >= 0 && last >= first, "forecastRange: need 0 <= first <= last");
    enforce!GretlException(method != ForecastMethod.dynamic, "forecastRange: a dynamic forecast is not applicable to a GARCH model (its mean equation has no lagged dependent variable)");
    const nOrig = fit.nOrig;
    auto f = modelForecast(fit.model, fit.dataset, first, last, method, futureExog, nOrig, last >= nOrig ? last - nOrig + 1 : 0);
    const t1m = fit.model.t1, t2m = fit.model.t2;

    // inside the sample the standard error is that of the conditional distribution, sqrt(h(t)): the library leaves it out, or (when the range runs past the
    // sample) gives h itself instead of its square root
    foreach (s; 0 .. cast(size_t) (last - first + 1)) {
        const t = first + cast(int) s;

        if (t >= t1m && t <= t2m) f.se.ptr.val[s] = fit.variance[t - t1m, 0] ^^ 0.5;
    }
    return f;
}

private SeriesForecast modelForecast(Model * m, Dataset ds, int t1, int t2, ForecastMethod method, GretlMatrix futureExog, int exogStart, size_t exogRows)
        @trusted {
    import gretl.forecast: free_fit_resid, get_forecast;
    auto dset = ds.ptr;
    int[] ids;                  // the regressors of the model other than the constant

    if (m.xlist !is null) {
        foreach (i; 1 .. m.xlist[0] + 1) {
            if (m.xlist[i] != 0) ids ~= m.xlist[i];
        }
    }
    checkExog(exogRows, 0, ids.length, futureExog, 0);
    auto guard = ExogGuard(dset, ids, futureExog, exogStart, ids.length > 0 ? exogRows : 0, t2);
    int err = 0;
    auto fr = get_forecast(m, t1, t2, 0, dset, methodOpt(method), &err);

    enforce!GretlException(err == 0 && fr !is null, "forecast: failed, " ~ lastError());
    scope (exit) free_fit_resid(fr);
    const steps = t2 - t1 + 1;
    SeriesForecast f;

    f.first = t1;
    f.mean = fresh(steps, 1);
    f.se = fresh(steps, 1);
    foreach (s; 0 .. steps) {
        f.mean.ptr.val[s] = fr.fitted[t1 + s];
        f.se.ptr.val[s] = fr.sderr !is null ? fr.sderr[t1 + s] : double.nan;
    }
    return f;
}

/* ---- the periodogram and fractional integration ---- */

/* The sample spectrum of a series at the Fourier frequencies 2 pi j / n, j = 1 .. n / 2 (each n/2 x 1). */
struct Periodogram {
    GretlMatrix frequency;      // radians, between 0 and pi
    GretlMatrix period;         // 2 pi / frequency, in observations
    GretlMatrix density;        // the spectral density estimate (the library's scaling: the periodogram is |DFT|^2 / (2 pi n))
    int window;                 // the Bartlett window width used (-1: none, the plain sample periodogram)
}

/* periodogram(x, window): the sample periodogram of x (at least 12 observations, not constant). window = -1 (default) gives the plain periodogram; window = 0 smooths it
   with the Bartlett lag window of the library's automatic width (2 sqrt(n)); window > 0 smooths with that width (at most n / 2). */
Periodogram periodogram(X)(X x, int window = -1) @trusted if (isOperand!X) {
    auto v = seriesOf(x, "periodogram");
    const n = v.length;

    enforce!GretlException(n >= 12, "periodogram: too few observations (at least 12 are needed)");
    enforce!GretlException(window >= -1, "periodogram: window must be -1 (none), 0 (automatic) or positive");
    int err = 0;
    auto m = periodogram_matrix(v.ptr, 0, cast(int) n - 1, window, &err);

    enforce!GretlException(err == 0 && m !is null, "periodogram: " ~ lastError());
    scope (exit) gretl_matrix_free(m);
    Periodogram r;
    const h = m.rows;

    r.frequency = fresh(h, 1);
    r.period = fresh(h, 1);
    r.density = fresh(h, 1);
    foreach (i; 0 .. h) {
        r.frequency.ptr.val[i] = m.val[i];
        r.period.ptr.val[i] = 2 * 3.14159265358979323846 / m.val[i];
        r.density.ptr.val[i] = m.val[h + i];
    }
    r.window = window == 0 ? auto_spectrum_order(cast(int) n, Gretlopt.OPT_O) : (window < 0 ? -1 : (window > cast(int) n / 2 ? cast(int) n / 2 : window));
    return r;
}

/// ditto, for a series of a Dataset
Periodogram periodogram(Dataset d, string name, int window = -1) @trusted {
    return periodogram(d[name].dup, window);
}

/* Estimates of the degree of fractional integration d of a series (d = 0: stationary, d = 1: a unit root), from the frequencies up to `bandwidth`. */
struct FractionalIntegration {
    int bandwidth;              // the number of Fourier frequencies used
    double lweD, lweSe;         // the local Whittle estimator (Robinson), its standard error 1 / (2 sqrt(m))
    double lweZ, lwePvalue;     // z = d / se for the null of d = 0, and the two-sided normal p-value
    double gphD, gphSe;         // the Geweke and Porter-Hudak log-periodogram regression estimate and its standard error
    double gphT, gphPvalue;     // t = d / se with bandwidth - 2 degrees of freedom, and the two-sided p-value
}

/* fractionalIntegration(x, bandwidth): both estimators. bandwidth = 0 (default) takes the library's choice, about n^0.6; otherwise at most n / 2. Needs 12 or more observations. */
FractionalIntegration fractionalIntegration(X)(X x, int bandwidth = 0) @trusted if (isOperand!X) {
    auto v = seriesOf(x, "fractionalIntegration");
    const n = v.length;

    enforce!GretlException(n >= 12, "fractionalIntegration: too few observations (at least 12 are needed)");
    enforce!GretlException(bandwidth >= 0 && bandwidth != 1, "fractionalIntegration: the bandwidth must be 0 (automatic) or at least 2");
    auto ds = timeDataset([v], 1, "1", "fractionalIntegration");
    auto pr = Printer(0);
    auto list = makeList([1]);

    enforce!GretlException(push_program_state() == 0, "could not save the library settings");
    scope (exit) {
        pop_program_state();
        import core.stdc.stdlib: free;
        free(list);
    }
    const err = fractint(1, bandwidth, ds.ptr, Gretlopt.OPT_A | Gretlopt.OPT_Q, pr.prn);

    enforce!GretlException(err == 0, "fractionalIntegration: it failed, " ~ lastError());
    GretlType type;
    int e2 = 0;
    auto res = cast(GretlMatrixRaw *) get_last_result_data(&type, &e2);

    enforce!GretlException(e2 == 0 && res !is null && res.rows == 2, "fractionalIntegration: no result recorded");
    scope (exit) gretl_matrix_free(res);
    FractionalIntegration r;

    r.bandwidth = bandwidth <= 0 ? auto_spectrum_order(cast(int) n, Gretlopt.OPT_NONE) : (bandwidth > cast(int) n / 2 ? cast(int) n / 2 : bandwidth);
    r.lweD = res.val[0];
    r.lweSe = res.val[2];
    r.gphD = res.val[1];
    r.gphSe = res.val[3];
    r.lweZ = r.lweD / r.lweSe;
    r.lwePvalue = get_last_pvalue();
    r.gphT = r.gphD / r.gphSe;
    r.gphPvalue = student_pvalue_2(r.bandwidth - 2, r.gphT);
    return r;
}

/* ---- the tests, filters and models on the series of a Dataset ---- */

/* The functions below take a Dataset and the NAMES of the series instead of an array or matrix, and do what the forms for matrices do: correlogram(d, "gdp", 12),
   adf(d, "gdp", 2), hpTrend(d, "gdp"), ar1(d, "y", ["x1", "x2"]), lagSelection(d, ["y1", "y2"], 8), .... Coefficients of the models carry the names of the series.
   The series must have no missing values (a NaN throws GretlException: choose the rows first); `adf` takes the calendar from the Dataset if it has one. */

private double[] named(Dataset d, string name, string who) @trusted {
    return seriesOf(d[name].dup, who);
}

private GretlMatrix columnsNamed(Dataset d, string[] names, string who) @trusted {
    enforce!GretlException(names.length > 0, who ~ ": no series named");
    double[][] cols;

    foreach (nm; names) cols ~= named(d, nm, who);
    return matrixFrom(cols);
}

/// ditto
Correlogram correlogram(Dataset d, string name, int maxLag) @trusted {
    return correlogram(named(d, name, "correlogram"), maxLag);
}

/// ditto
TestResult ljungBox(Dataset d, string name, int lags, int fitted = 0) @trusted {
    return ljungBox(named(d, name, "ljungBox"), lags, fitted);
}

/// ditto
TestResult archTest(Dataset d, string name, int order) @trusted {
    return archTest(named(d, name, "archTest"), order);
}

/// ditto
UnitRootTest adf(Dataset d, string name, int lags = 0, Trend trend = Trend.constant, TestDown testDown = TestDown.off, bool gls = false, bool seasonals = false,
                 bool difference = false, int frequency = 0, string start = null) @trusted {
    calendarOf(d, frequency, start);
    return adf(named(d, name, "adf"), lags, trend, testDown, gls, seasonals, difference, frequency, start);
}

/// ditto
UnitRootTest kpss(Dataset d, string name, int lags, Trend trend = Trend.constant, bool difference = false) @trusted {
    return kpss(named(d, name, "kpss"), lags, trend, difference);
}

/// ditto
UnitRootTest engleGranger(Dataset d, string y, string[] x, int lags = 0, Trend trend = Trend.constant, TestDown testDown = TestDown.off) @trusted {
    return engleGranger(named(d, y, "engleGranger"), columnsNamed(d, x, "engleGranger"), lags, trend, testDown);
}

/// ditto
GretlMatrix hpTrend(Dataset d, string name, double lambda = 1600) @trusted {
    return hpTrend(named(d, name, "hpTrend"), lambda);
}

/// ditto
GretlMatrix hpCycle(Dataset d, string name, double lambda = 1600) @trusted {
    return hpCycle(named(d, name, "hpCycle"), lambda);
}

/// ditto
GretlMatrix bkFilter(Dataset d, string name, int lower = 6, int upper = 32, int k = 12) @trusted {
    return bkFilter(named(d, name, "bkFilter"), lower, upper, k);
}

/// ditto
GretlMatrix butterworth(Dataset d, string name, int order = 2, double cutoff = 67.5) @trusted {
    return butterworth(named(d, name, "butterworth"), order, cutoff);
}

/// ditto
GretlMatrix fracDiff(Dataset d, string name, double dd) @trusted {
    return fracDiff(named(d, name, "fracDiff"), dd);
}

/// ditto: one column for each series named
GretlMatrix lag(Dataset d, string[] names, int k = 1) @trusted {
    return lag(columnsNamed(d, names, "lag"), k);
}

/// ditto
GretlMatrix lag(Dataset d, string name, int k = 1) @trusted {
    return lag(d, [name], k);
}

/// ditto
GretlMatrix diff(Dataset d, string[] names, int order = 1, int seasonal = 1) @trusted {
    return diff(columnsNamed(d, names, "diff"), order, seasonal);
}

/// ditto
GretlMatrix diff(Dataset d, string name, int order = 1, int seasonal = 1) @trusted {
    return diff(d, [name], order, seasonal);
}

/* the generated names x1, x2, ... of the regressors of a model replaced by the names of the series */
private void nameRegressors(ref Estimate e, string[] xnames) {
    foreach (ref nm; e.names) {
        foreach (j, s; xnames) {
            if (nm == "x" ~ to!string(j + 1)) nm = s;
        }
    }
}

/// ditto
Estimate ar1(Dataset d, string y, string[] x, bool constant = true, Ar1Method method = Ar1Method.cochraneOrcutt) @trusted {
    auto e = ar1(named(d, y, "ar1"), columnsNamed(d, x, "ar1"), constant, method);

    nameRegressors(e, x);
    return e;
}

/// ditto
Estimate ar(Dataset d, string y, string[] x, int[] lags, bool constant = true) @trusted {
    auto e = ar(named(d, y, "ar"), columnsNamed(d, x, "ar"), lags, constant);

    nameRegressors(e, x);
    return e;
}

/// ditto
Estimate ar(Dataset d, string y, string[] x, int p, bool constant = true) @trusted {
    auto e = ar(named(d, y, "ar"), columnsNamed(d, x, "ar"), p, constant);

    nameRegressors(e, x);
    return e;
}

/// ditto
Estimate arch(Dataset d, string y, string[] x, int order, bool constant = true) @trusted {
    auto e = arch(named(d, y, "arch"), columnsNamed(d, x, "arch"), order, constant);

    nameRegressors(e, x);
    return e;
}

/// ditto: `exog` names the exogenous series; the calendar of the Dataset is used
LagSelection lagSelection(Dataset d, string[] endog, int maxLag, string[] exog = null, bool constant = true, bool trend = false, bool seasonals = false,
                          int frequency = 0, string start = null, string fitStart = null, string fitEnd = null) @trusted {
    calendarOf(d, frequency, start);
    return lagSelection(columnsNamed(d, endog, "lagSelection"), maxLag, exog.length > 0 ? columnsNamed(d, exog, "lagSelection") : GretlMatrix.init, constant, trend,
                        seasonals, frequency, start, fitStart, fitEnd);
}
