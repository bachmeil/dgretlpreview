/* Port of adf_kpss.c */

module gretl.adf_kpss;


/*
 * Implementations of the (Augmented) Dickey-Fuller test
 * and the Kwiatkowski, Phillips, Schmidt and Shin test for
 * the presence of a unit root in a time series, along with
 * the Engle-Granger test for cointegration of two or more
 * time series.
 *
 * The Johansen cointegration test is also provided in
 * libgretl; see johansen_test() and johansen_test_simple().
 */

import core.stdc.math : floor, log, pow, sqrt;
import core.stdc.stdio : fprintf, sprintf, stderr;
import core.stdc.stdlib : free, getenv, malloc;
import core.stdc.string : strcmp, strcpy;
import std.algorithm.comparison : min;

import gretl.compareh;
import gretl.dataset : dataset_add_allocated_series, dataset_add_series, dataset_drop_last_variables, dataset_drop_variable, multi_unit_panel_sample, series_set_label;
import gretl.dataseth;
import gretl.estimate : lsq;
import gretl.glib;
import gretl.gmem;
import gretl.gretl_commandsh;
import gretl.gretl_errorsh;
import gretl.gretl_intlh;
import gretl.gretl_matrix : gretl_matrix_alloc, gretl_matrix_free, gretl_matrix_ols,
    gretl_matrix_set_rownames;
import gretl.gretl_matrixh;
import gretl.gretl_model : clear_model, gretl_model_init, gretl_model_set_double, gretl_model_set_int;
import gretl.gretl_prn : Prn, pprintf, pputc, pputs;
import gretl.gstrfuncs : g_strdup_printf;
import gretl.libgretlh;
import gretl.missingh;
import gretl.pvaluesh;
import gretl.pvalues : chisq_cdf_comp, normal_cdf_inverse, normal_pvalue_1, normal_pvalue_2,
    snedecor_cdf_comp, student_pvalue_1;
import gretl.gretl_list : gretl_list_copy, gretl_list_delete_at_pos, gretl_list_new;
import gretl.strutils : gretl_strdup, strings_array_free, strings_array_new;
import gretl.options : get_optval_string, incompatible_options, options_incompatible_with;
import gretl.gretl_errors : gretl_errmsg_set, gretl_errmsg_sprintf;
import gretl.missing : list_adjust_sample, series_adjust_sample;
import gretl.transforms : diffgenr, laggenr;
import gretl.gretl_utils : gretl_isconst, record_matrix_test_result, record_test_result;
import gretl.dataio : ntolabel;
import gretl.plugins : get_plugin_function;
import gretl.genmain : series_index;
import gretl.modelprint : printmodel;
import gretl.compare : wald_omit_F;
import gretl.genfuncs : seasonals_list;

/* ---- declarations of not-yet-ported libgretl items used below ----
   Replace each with an import when its defining file is ported. */

extern (C) {
    /* printf-style variadics */

    /* models */

    /* matrices */

    /* lists and strings */

    /* dataset */

    /* transforms and genfuncs */

    /* distributions */

    /* options, errors, test records, plugins */
}

/* ---- end of declarations ---- */

/* codes for deterministic regressors */
enum DetCode {
    no_const = 1,
    const_,
    trend,
    quad_trend,
    max
}

/* flags for "special stuff" going on */
enum AdfFlags {
    eg_test   = 1 << 0, /* doing Engle-Granger test */
    eg_resids = 1 << 1, /* final stage of the above */
    panel     = 1 << 2, /* working on panel data */
    gls       = 1 << 3, /* GLS: Elliott, Rothenberg, Stock */
    pq        = 1 << 4, /* Perron-Qu, 2007 */
    test_down = 1 << 5  /* Testing down for lag order */
}

/* automatic lag selection methods */
enum LagSelection {
    aic = 1,
    bic,
    tstat
}

/* demean/detrend method */
enum DemeanMethod {
    gls = 1,
    ols
}

/* p-value type */
enum PvalueType {
    finite = 1,
    asy,
    approx
}

struct AdfInfo {
    int v;           /* ID number of series to test (in/out) */
    int order;       /* lag order for ADF (in/out) */
    int kmax;        /* max. order (for testing down) */
    int altv;        /* ID of modified series (detrended) */
    int niv;         /* number of (co-)integrated vars (Engle-Granger) */
    AdfFlags flags = cast(AdfFlags) 0; /* bitflags: see above */
    DetCode det = cast(DetCode) 0;     /* code for deterministics */
    int pvtype;      /* code for p-value type (see PvalueType; 0 = unset) */
    int nseas;       /* number of seasonal dummies */
    int T;           /* number of obs used in test */
    int df;          /* degrees of freedom, test regression */
    double b0 = 0;   /* coefficient on lagged level */
    double tau = 0;  /* test statistic */
    double pval = 0; /* p-value of test stat */
    int * list;      /* regression list */
    int * slist;     /* list of seasonal dummies, if applicable */
    const(char) * vname; /* name of series tested */
    GretlMatrixRaw * g; /* GLS coefficients (if applicable) */
    int verbosity;
}

struct KpssInfo {
    int T;
    int order;
    int pmin;
    int pmax;
    double test;
    double pval;
}

/* replace @y with demeaned or detrended y via GLS */

private int GLS_demean_detrend(DatasetRaw * dset, int v,
                               int offset, int T,
                               DetCode det,
                               AdfInfo * ainfo,
                               Prn * prn)
{
    double * y = dset.Z[v];
    GretlMatrixRaw * yd = null;
    GretlMatrixRaw * Xd = null;
    GretlMatrixRaw * b = null;
    double c;
    int i, t, xcols;
    int err = 0;

    xcols = (det == DetCode.trend) ? 2 : 1;

    if (T - xcols <= 0) {
        return GretlError.df;
    }

    yd = gretl_column_vector_alloc(T);
    Xd = gretl_matrix_alloc(T, xcols);
    b = gretl_column_vector_alloc(xcols);

    if (yd == null || Xd == null || b == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    c = (det == DetCode.const_) ? (1.0 - 7.0/T) : (1.0 - 13.5/T);

    /* invalidate pre-offset observations */
    for (t=0; t<offset; t++) {
        y[t] = NADBL;
    }
    y += offset;

    gretl_vector_set(yd, 0, y[0]);
    // gretl_vector_set(yd, 0, (1 - c) * y[0]);
    for (t=1; t<T; t++) {
        gretl_vector_set(yd, t, y[t] - c * y[t-1]);
    }

    gretl_matrix_set(Xd, 0, 0, 1);
    if (xcols == 2) {
        gretl_matrix_set(Xd, 0, 1, 1);
    }

    for (t=1; t<T; t++) {
        gretl_matrix_set(Xd, t, 0, 1 - c);
        if (xcols == 2) {
            gretl_matrix_set(Xd, t, 1, t+1 - t*c);
        }
    }

    err = gretl_matrix_ols(yd, Xd, b, null, null, null);

    if (!err && ainfo.verbosity > 1) {
        char[OBSLEN] date1;
        char[OBSLEN] date2;

        ntolabel(date1.ptr, offset, dset);
        ntolabel(date2.ptr, offset + T - 1, dset);
        pprintf(prn, _("\nGLS demean/detrend: using %s-%s (T = %d)\n").ptr,
                date1.ptr, date2.ptr, yd.rows);
        for (i=0; i<xcols; i++) {
            pprintf(prn, "  %s = %#.8g\n", ((i == 0) ? "mean" : "trend").ptr, b.val[i]);
        }
    }

    if (!err) {
        for (t=0; t<T; t++) {
            y[t] -= b.val[0];
            if (xcols == 2) {
                y[t] -= b.val[1] * (t+1);
            }
        }
    }

 bailout:

    gretl_matrix_free(yd);
    gretl_matrix_free(Xd);

    if (err) {
        gretl_matrix_free(b);
    } else {
        ainfo.g = b;
    }

    return err;
}

/* replace @y with demeaned or detrended y via OLS */

private int OLS_demean_detrend(double * y, int offset,
                               int T, DetCode det,
                               Prn * prn)
{
    GretlMatrixRaw * yd = null;
    GretlMatrixRaw * Xd = null;
    GretlMatrixRaw * b = null;
    int t, xcols = 1;
    int err = 0;

    if (det == DetCode.quad_trend) {
        xcols = 3;
    } else if (det == DetCode.trend) {
        xcols = 2;
    }

    if (T - xcols <= 0) {
        return GretlError.df;
    }

    yd = gretl_column_vector_alloc(T);
    Xd = gretl_matrix_alloc(T, xcols);
    b = gretl_column_vector_alloc(xcols);

    if (yd == null || Xd == null || b == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    /* invalidate pre-offset observations */
    for (t=0; t<offset; t++) {
        y[t] = NADBL;
    }
    y += offset;

    for (t=0; t<T; t++) {
        gretl_vector_set(yd, t, y[t]);
        gretl_matrix_set(Xd, t, 0, 1.0);
        if (xcols > 1) {
            gretl_matrix_set(Xd, t, 1, t+1);
        }
        if (xcols > 2) {
            gretl_matrix_set(Xd, t, 2, (t+1) * (t+1));
        }
    }

    err = gretl_matrix_ols(yd, Xd, b, null, null, null);

    if (!err && prn != null) {
        pprintf(prn, _("\nOLS demean/detrend: (T = %d)\n").ptr, yd.rows);
        for (t=0; t<xcols; t++) {
            pprintf(prn, "  %s = %#.8g\n", ((t == 0) ? "mean" : "trend").ptr, b.val[t]);
        }
    }

    if (!err) {
        for (t=0; t<T; t++) {
            y[t] -= b.val[0];
            if (xcols > 1) {
                y[t] -= b.val[1] * (t+1);
            }
            if (xcols > 2) {
                y[t] -= b.val[2] * (t+1) * (t+1);
            }
        }
    }

 bailout:

    gretl_matrix_free(yd);
    gretl_matrix_free(Xd);
    gretl_matrix_free(b);

    return err;
}

private int real_adf_form_list(AdfInfo * ainfo,
                               DatasetRaw * dset)
{
    int v, save_t1 = dset.t1;
    int k, j, err = 0;

    /* using the original var, or transformed? */
    v = ainfo.altv > 0 ? ainfo.altv : ainfo.v;

    /* temporararily reset sample */
    dset.t1 = 0;

    /* generate the first difference of series @v:
       this will be the LHS variable in the test
    */
    ainfo.list[1] = diffgenr(v, GretlCmdIndex.diff, dset);
    if (ainfo.list[1] < 0) {
        dset.t1 = save_t1;
        return GretlError.data;
    }

    /* generate lag 1 of series @v: the basic RHS series */
    ainfo.list[2] = laggenr(v, 1, dset);
    if (ainfo.list[2] < 0) {
        dset.t1 = save_t1;
        return GretlError.data;
    }

    /* generate lagged differences for augmented test */
    j = 3;
    for (k=1; k<=ainfo.order && !err; k++) {
        int vk = laggenr(ainfo.list[1], k, dset);

        if (vk < 0) {
            fprintf(stderr, "Error generating lags\n");
            err = GretlError.data;
        } else {
            ainfo.list[j++] = vk;
        }
    }

    if (!err && ainfo.nseas > 0) {
        /* note: centered */
        ainfo.slist = seasonals_list(dset, dset.pd, 1, &err);
        if (err) {
            ainfo.nseas = 0;
        } else {
            ainfo.nseas = ainfo.slist[0];
        }
    }

    /* restore incoming sample */
    dset.t1 = save_t1;

    return err;
}

/* The "offset" will determine where the data start for the
   initial detrending regression. I suppose we don't want
   to include excessive "pre-sample" data if the user has
   explicitly moved the sample start off zero, but we need
   the data to start early enough to provide @order + 1
   lags, to be in sync with the ADF regressions.
*/

private int adf_y_offset(AdfInfo * ainfo, int v, DatasetRaw * dset)
{
    int min_offset = dset.t1 - (ainfo.order + 1);
    int t, offset = 0;

    /* copy original data into series v */
    for (t=0; t<=dset.t2; t++) {
        dset.Z[v][t] = dset.Z[ainfo.v][t];
        if (na(dset.Z[v][t])) {
            offset = t+1;
        }
    }

    if (offset < min_offset) {
        offset = min_offset;
    }

    return offset;
}

/* Generate the various differences and lags required for
   the ADF test, detrending first if required.
*/

private int adf_prepare_vars(AdfInfo * ainfo, DatasetRaw * dset,
                             Gretlopt opt, int demean_mode,
                             Prn * prn)
{
    int err = 0;

    if (ainfo.v == 0) {
        return GretlError.data;
    }

    if (ainfo.list == null) {
        /* the max number of terms (quadratic trend case) */
        int len = 5 + ainfo.nseas + ainfo.order;

        ainfo.list = gretl_list_new(len);
        if (ainfo.list == null) {
            return GretlError.alloc;
        }
    }

    if (demean_mode == DemeanMethod.ols) {
        /* OLS adjustment is wanted (first pass) */
        DetCode det = DetCode.const_;
        int v = dset.v;

        if (opt & Gretlopt.OPT_R) {
            det = DetCode.quad_trend;
        } else if (opt & Gretlopt.OPT_T) {
            det = DetCode.trend;
        }

        err = dataset_add_series(dset, 1);

        if (!err) {
            Prn * vprn = ainfo.verbosity > 1 ? prn : null;
            int offset = adf_y_offset(ainfo, v, dset);
            int T = dset.t2 - offset + 1;

            err = OLS_demean_detrend(dset.Z[v], offset, T, det, vprn);
        }

        if (!err) {
            /* replace with OLS-detrended version */
            strcpy(dset.varname[v], "ydetols");
            ainfo.altv = v;
        }
    } else if (demean_mode == DemeanMethod.gls) {
        /* GLS adjustment is wanted */
        DetCode det = (opt & Gretlopt.OPT_T) ? DetCode.trend : DetCode.const_;
        int v;

        if (ainfo.altv > 0) {
            /* if we already did detrending, re-use the
               storage we added earlier */
            v = ainfo.altv;
        } else {
            v = dset.v;
            err = dataset_add_series(dset, 1);
        }

        if (!err) {
            int offset = adf_y_offset(ainfo, v, dset);
            int T = dset.t2 - offset + 1;

            err = GLS_demean_detrend(dset, v, offset, T,
                                     det, ainfo, prn);
        }

        if (!err) {
            /* replace with GLS-detrended version */
            strcpy(dset.varname[v], "ydetrend");
            ainfo.altv = v;
        }
    }

    if (!err) {
        err = real_adf_form_list(ainfo, dset);
    }

    return err;
}

/* Peter Sephton's response-surface based critical values for
   ADF-GLS with testing down of lag order as per Ng and Perron,
   or Perron and Qu. See PS's paper in Computational Economics,
   https://doi.org/10.1007/s10614-020-10082-6
*/

private void sephton_adf_critvals(AdfInfo * ainfo, double * c)
{
    /* Each row contains 6 coeffs: const, 4 powers of 1/T,
       and maxlag / T. The rows are in blocks of four:
       (1) Ng-Perron constant; (2) Perron-Qu constant;
       (3) Ng-Perron trend;    (3) Perron-Qu trend.
       Within each block are coefficients for 4 critical
       values: 0.01, 0.025, 0.05 and 0.1.
    */
    static immutable double[6][16] coef = [
        /* Ng-Perron constant */
        [-2.5512, 1.084, -720.43, 18279, -175060, 0.3174],
        [-2.2186, -3.7919, -333.25, 6665.6, -52297, 0.31557],
        [-1.9341, -4.7232, -124.1, 1441.1, -6754.7, 0.0096449],
        [-1.6136, -8.9613, 120.77, -3103.9, 24355, -0.15212],
        /* Perron-Qu constant */
        [-2.5449, -3.2252,  -626, 14049, -107060, 0.33716],
        [-2.2154, -9.9681, -250.85, 3955.5, -12315, 0.5985],
        [-1.932, -7.5717, -259.69, 5407.7, -31513, 0.2145],
        [-1.6108, -8.9303, -83.401, 2341.3, -13132, -0.10669],
        /* Ng-Perron trend */
        [-3.3772, -0.2059, -1809.6, 46527, -415670,  1.27],
        [-3.0761, 2.9114, -1466.4, 35993, -312430, 0.7397],
        [-2.828, 7.4946, -1579.2, 39383, -341880, 0.47643],
        [-2.545, 5.5607, -1133.2, 26203, -213800, 0.24246],
        /* Perron-Qu trend */
        [-3.3719, -14.027, -1194.7, 32450, -301000, 1.6079],
        [-3.0688, -4.2253, -1135, 29133, -257020, 0.73835],
        [-2.8196, 2.4706, -1346.6, 35398, -314840, 0.37564],
        [-2.5401, -2.7442, -950.57, 25683, -233160, 0.52035]
    ];
    int PQ = (ainfo.flags & AdfFlags.pq) ? 1 : 0;
    int trend = (ainfo.det == DetCode.trend);
    /* NB: precedence as in the C source: (4 * PQ + trend) ? 8 : 0 */
    int i, j, rmin = 4 * PQ + trend ? 8 : 0;
    const(double) * b;
    double[6] X;

    X[0] = 1.0;
    X[1] = 1.0 / ainfo.T;
    X[2] = X[1] * X[1];
    X[3] = X[2] * X[1];
    X[4] = X[3] * X[1];
    X[5] = ainfo.kmax * X[1];

    for (i=0; i<4; i++) {
        b = coef[rmin+i].ptr;
        c[i] = 0;
        for (j=0; j<6; j++) {
            c[i] += X[j] * b[j];
        }
    }
}

/* display an F-test for the joint significance of the lagged
   \delta y terms in ADF test */

private void show_lags_test(Model * pmod, int order, Prn * prn)
{
    int * llist = gretl_list_new(order);
    double F;
    int i;

    if (llist != null) {
        for (i=0; i<order; i++) {
            /* lagged differences */
            llist[i+1] = pmod.list[pmod.ifc + 3 + i];
        }

        F = wald_omit_F(llist, pmod);

        if (!na(F)) {
            pprintf(prn, "  %s: F(%d, %d) = %.3f [%.4f]\n",
                    _("lagged differences").ptr, order, pmod.dfd, F,
                    snedecor_cdf_comp(order, pmod.dfd, F));
        }

        free(llist);
    }
}

private string test_down_string(int i, Gretlopt opt)
{
    if (opt & Gretlopt.OPT_U) {
        /* perron-qu */
        if (i == LagSelection.bic) {
            return _("modified BIC, Perron-Qu");
        } else {
            return _("modified AIC, Perron-Qu");
        }
    } else {
        int gls = (opt & Gretlopt.OPT_G);

        if (i == LagSelection.bic) {
            return gls ? _("modified BIC") : _("BIC");
        } else if (i == LagSelection.tstat) {
            return _("t-statistic");
        } else {
            /* the default */
            return gls ? _("modified AIC") : _("AIC");
        }
    }
}

private void DF_header(const(char) * s, int p, int pmax,
                       int test_down, Gretlopt opt,
                       Prn * prn)
{
    pputc(prn, '\n');

    if (p <= 0 && pmax == 0) {
        if (opt & Gretlopt.OPT_G) {
            pprintf(prn, _("Dickey-Fuller (GLS) test for %s\n").ptr, s);
        } else {
            pprintf(prn, _("Dickey-Fuller test for %s\n").ptr, s);
        }
    } else {
        if (opt & Gretlopt.OPT_G) {
            pprintf(prn, _("Augmented Dickey-Fuller (GLS) test for %s\n").ptr, s);
        } else {
            pprintf(prn, _("Augmented Dickey-Fuller test for %s\n").ptr, s);
        }
        if (test_down) {
            string critstr = test_down_string(test_down, opt);

            pprintf(prn, _("testing down from %d lags, criterion %s").ptr,
                    pmax, critstr.ptr);
        } else {
            if (p == 1) {
                pprintf(prn, _("including one lag of (1-L)%s").ptr, s);
            } else {
                pprintf(prn, _("including %d lags of (1-L)%s").ptr, p, s);
            }
        }
        pputc(prn, '\n');
    }
}

private string DF_model_string(int i)
{
    static immutable string[4] models = [
        "(1-L)y = (a-1)*y(-1) + e",
        "(1-L)y = b0 + (a-1)*y(-1) + e",
        "(1-L)y = b0 + b1*t + (a-1)*y(-1) + e",
        "(1-L)y = b0 + b1*t + b2*t^2 + (a-1)*y(-1) + e"
    ];

    if (i >= 0 && i < 4) {
        return models[i];
    } else {
        return "";
    }
}

private string ADF_model_string(int i)
{
    static immutable string[4] models = [
        "(1-L)y = (a-1)*y(-1) + ... + e",
        "(1-L)y = b0 + (a-1)*y(-1) + ... + e",
        "(1-L)y = b0 + b1*t + (a-1)*y(-1) + ... + e",
        "(1-L)y = b0 + b1*t + b2*t^2 + (a-1)*y(-1) + ... + e"
    ];

    if (i >= 0 && i < 4) {
        return models[i];
    } else {
        return "";
    }
}

private string DF_test_string(int i)
{
    static immutable string[4] tests = [
        N_("test without constant"),
        N_("test with constant"),
        N_("with constant and trend"),
        N_("with constant, linear and quadratic trend")
    ];

    if (i >= 0 && i < 4) {
        return tests[i];
    } else {
        return "";
    }
}

private void print_adf_results(AdfInfo * ainfo, Model * dfmod,
                               int * blurb_done, Gretlopt opt,
                               int test_down, Prn * prn)
{
    static immutable string[4] urcstrs = [
        "nc", "c", "ct", "ctt"
    ];
    gchar * pvstr = null;
    char[16] taustr;
    int i;

    if (prn == null) return;

    /* convert deterministics code to 0-base */
    i = ainfo.det - 1;

    if (na(ainfo.pval)) {
        pvstr = g_strdup_printf("%s %s", _("p-value").ptr, _("unknown").ptr);
    } else if (ainfo.pvtype == PvalueType.asy) {
        pvstr = g_strdup_printf("%s %.4g", _("asymptotic p-value").ptr, ainfo.pval);
    } else if (ainfo.pvtype == PvalueType.approx) {
        pvstr = g_strdup_printf("%s %.3f", _("approximate p-value").ptr, ainfo.pval);
    } else {
        pvstr = g_strdup_printf("%s %.4g", _("p-value").ptr, ainfo.pval);
    }

    if (*blurb_done == 0) {
        DF_header(ainfo.vname, ainfo.order, ainfo.kmax,
                  test_down, opt, prn);
        pprintf(prn, _("sample size %d\n").ptr, ainfo.T);
        if (ainfo.flags & AdfFlags.panel) {
            pputc(prn, '\n');
        } else {
            pputs(prn, _("unit-root null hypothesis: a = 1").ptr);
            pputs(prn, "\n\n");
        }
        *blurb_done = 1;
    }

    if (ainfo.flags & AdfFlags.eg_resids) {
        /* last step of Engle-Granger test */
        pprintf(prn, "  %s ", _(DF_test_string(0)).ptr);
        pputc(prn, '\n');
        if (test_down) {
            pputs(prn, "  ");
            if (ainfo.order == 1) {
                pprintf(prn, _("including one lag of (1-L)%s").ptr, ainfo.vname);
            } else {
                pprintf(prn, _("including %d lags of (1-L)%s").ptr,
                        ainfo.order, ainfo.vname);
            }
            pputc(prn, '\n');
        }
        pprintf(prn, "  %s: %s\n", _("model").ptr,
                ((ainfo.order > 0) ? ADF_model_string(0) :
                 DF_model_string(0)).ptr);
    } else {
        pprintf(prn, "  %s ", _(DF_test_string(i)).ptr);
        if (ainfo.nseas > 0) {
            pputs(prn, _("plus seasonal dummies").ptr);
        }
        pputc(prn, '\n');
        if (test_down) {
            pputs(prn, "  ");
            if (ainfo.order == 1) {
                pprintf(prn, _("including one lag of (1-L)%s").ptr, ainfo.vname);
            } else {
                pprintf(prn, _("including %d lags of (1-L)%s").ptr,
                        ainfo.order, ainfo.vname);
            }
            pputc(prn, '\n');
        }
        pprintf(prn, "  %s: %s\n", _("model").ptr,
                ((ainfo.order > 0) ? ADF_model_string(i) :
                 DF_model_string(i)).ptr);
    }

    if (ainfo.flags & AdfFlags.gls) {
        strcpy(taustr.ptr, "tau");
    } else {
        sprintf(taustr.ptr, "tau_%s(%d)", urcstrs[i].ptr, ainfo.niv);
    }

    pprintf(prn, "  %s: %g\n"
            ~ "  %s: %s = %g\n",
            _("estimated value of (a - 1)").ptr, ainfo.b0,
            _("test statistic").ptr, taustr.ptr, ainfo.tau);

    if ((ainfo.flags & AdfFlags.gls) && na(ainfo.pval)) {
        double[4] c;

        sephton_adf_critvals(ainfo, c.ptr);
        pprintf(prn, "\n  %*s    ", cast(int) TRANSLATED_WIDTH(_("Critical values")), " ".ptr);
        pprintf(prn, "%g%%     %g%%     %g%%     %g%%\n", 10.0, 5.0, 2.5, 1.0);
        pprintf(prn, "  %s: %.2f   %.2f   %.2f   %.2f\n",
                _("Critical values").ptr, c[3], c[2], c[1], c[0]);
    } else {
        pprintf(prn, "  %s\n", pvstr);
    }

    g_free(pvstr);

    if (!na(dfmod.rho)) {
        pprintf(prn, "  %s: %.3f\n", _("1st-order autocorrelation coeff. for e").ptr,
                dfmod.rho);
    }

    if (ainfo.order > 1) {
        show_lags_test(dfmod, ainfo.order, prn);
    }
}

private Gretlopt ur_lsq_opt(const(DatasetRaw) * dset, Gretlopt opt)
{
    if (dset.pd == 52 && (opt & Gretlopt.OPT_D)) {
        /* weekly seasonals present */
        return Gretlopt.OPT_A;
    } else {
        return Gretlopt.OPT_A | Gretlopt.OPT_Z;
    }
}

/* test the lag order down using the t-statistic criterion */

private int t_adjust_order(AdfInfo * ainfo, DatasetRaw * dset,
                           Gretlopt opt, int * err, Prn * prn)
{
    Gretlopt kmod_opt;
    Model kmod;
    int kmax = ainfo.kmax;
    double tstat, pval;
    int k, pos;

    kmod_opt = ur_lsq_opt(dset, opt);

    for (k=kmax; k>0; k--) {
        kmod = lsq(ainfo.list, dset, GretlCmdIndex.ols, kmod_opt);
        if (!kmod.errcode && kmod.dfd == 0) {
            kmod.errcode = GretlError.df;
        }
        if (kmod.errcode) {
            fprintf(stderr, "t_adjust_order: k = %d, err = %d\n", k,
                    kmod.errcode);
            *err = kmod.errcode;
            clear_model(&kmod);
            k = -1;
            break;
        }
        pos = k + kmod.ifc;
        tstat = kmod.coeff[pos] / kmod.sderr[pos];
        clear_model(&kmod);
        pval = normal_pvalue_2(tstat);

        if (pval > 0.10) {
            gretl_list_delete_at_pos(ainfo.list, k + 2);
        } else {
            break;
        }
    }

    return k;
}

/* compute modified information criterion */

private double get_MIC(Model * pmod, int k, double sum_ylag2,
                       LagSelection kmethod, const(DatasetRaw) * dset)
{
    double g, CT, ttk, s2k = 0;
    int t, T = pmod.nobs;

    g = pmod.coeff[pmod.ifc];

    for (t=pmod.t1; t<=pmod.t2; t++) {
        s2k += pmod.uhat[t] * pmod.uhat[t];
    }

    s2k /= pmod.nobs;
    ttk = g * g * sum_ylag2 / s2k;
    CT = kmethod == LagSelection.bic ? log(T) : 2.0;

    return log(s2k) + CT * (ttk + k)/T;
}

/* component calculation for modified IC methods */

private double get_sum_y2(AdfInfo * ainfo, Model * pmod,
                          const(DatasetRaw) * dset)
{
    const(double) * ylag = dset.Z[ainfo.list[2]];
    double sumy2 = 0;
    int t;

    for (t=pmod.t1; t<=pmod.t2; t++) {
        sumy2 += ylag[t] * ylag[t];
    }

    return sumy2;
}

/* Using modified information criterion, as per Ng and Perron,
   "Lag Length Selection and the Construction of Unit Root Tests
   with Good Size and Power", Econometrica 69/6, Nov 2001, pp.
   1519-1554, for the GLS case -- otherwise plain IC (as of
   2015-03-31).
*/

private int ic_adjust_order(AdfInfo * ainfo, LagSelection kmethod,
                            DatasetRaw * dset, Gretlopt opt,
                            int test_num, int * err,
                            Prn * prn)
{
    Model kmod;
    Gretlopt kmod_opt;
    double IC, ICmin = 0;
    double sum_ylag2 = 0;
    int kmax = ainfo.kmax;
    int k, kopt = kmax;
    int save_t1 = dset.t1;
    int save_t2 = dset.t2;
    int use_MIC = 0;
    int * tmplist;

    kmod_opt = ur_lsq_opt(dset, opt);

    tmplist = gretl_list_copy(ainfo.list);
    if (tmplist == null) {
        *err = GretlError.alloc;
        return -1;
    }

    if (ainfo.flags & AdfFlags.gls) {
        /* modified criterion wanted */
        use_MIC = 1;
    }

    for (k=kmax; k>=0; k--) {
        kmod = lsq(tmplist, dset, GretlCmdIndex.ols, kmod_opt);
        if (!kmod.errcode && kmod.dfd == 0) {
            kmod.errcode = GretlError.df;
        }
        if (kmod.errcode) {
            fprintf(stderr, "ic_adjust_order: k = %d, err = %d\n", k,
                    kmod.errcode);
            *err = kmod.errcode;
            clear_model(&kmod);
            kopt = -1;
            break;
        }
        if (use_MIC) {
            if (k == kmax) {
                /* this need only be done once */
                sum_ylag2 = get_sum_y2(ainfo, &kmod, dset);
            }
            IC = get_MIC(&kmod, k, sum_ylag2, kmethod, dset);
        } else if (kmethod == LagSelection.bic) {
            IC = kmod.criterion[InfoCriterion.bic];
        } else {
            IC = kmod.criterion[InfoCriterion.aic];
        }
        if (k == kmax) {
            /* ensure a uniform sample */
            dset.t1 = kmod.t1;
            dset.t2 = kmod.t2;
            ICmin = IC;
        } else if (IC < ICmin) {
            ICmin = IC;
            kopt = k;
        }
        if (ainfo.verbosity > 1) {
            printmodel(&kmod, dset, Gretlopt.OPT_S, prn);
        }
        if (ainfo.verbosity) {
            const(char) * tag;

            if (use_MIC) {
                tag = (kmethod == LagSelection.bic) ? "MBIC".ptr : "MAIC".ptr;
            } else {
                tag = (kmethod == LagSelection.bic) ? "BIC".ptr : "AIC".ptr;
            }

            if (0 && k == kmax && test_num == 1) {
                /* FIXME */
                pputc(prn, '\n');
            }
            pprintf(prn, "  k = %2d: %s = %#g\n", k, tag, IC);
        }
        clear_model(&kmod);
        gretl_list_delete_at_pos(tmplist, k + 2);
    }

    if (ainfo.verbosity && test_num > 1) {
        pputc(prn, '\n');
    }

    free(tmplist);

    if (kopt >= 0) {
        /* now trim the "real" list to @kopt lags */
        for (k=kmax; k>kopt; k--) {
            gretl_list_delete_at_pos(ainfo.list, k + 2);
        }
    }

    dset.t1 = save_t1;
    dset.t2 = save_t2;

    return kopt;
}

/* targ must be big enough to accept all of src! */

private void copy_list_values(int * targ, const(int) * src)
{
    int i;

    for (i=0; i<=src[0]; i++) {
        targ[i] = src[i];
    }
}

/*
 * get_urc_pvalue:
 * @tau: test statistic.
 * @n: sample size (or 0 for asymptotic result).
 * @niv: number of potentially cointegrated variables
 * (1 for simple unit-root test).
 * @itv: code: 1, 2, 3, 4 for nc, c, ct, ctt models
 * respectively.
 *
 * Retrieves the p-value for @tau from the Dickey–Fuller
 * unit-root test or the Engle–Granger cointegration
 * test, as per James MacKinnon (1996).
 *
 * Returns: p-value, or %NADBL on failure.
 */

double get_urc_pvalue(double tau, int n, int niv, int itv)
{
    alias PvFunc = extern (C) double function(double, int, int, int);
    PvFunc pvfunc;
    double pval = NADBL;

    pvfunc = cast(PvFunc) get_plugin_function("mackinnon_pvalue");
    if (pvfunc == null) {
        return pval;
    }

    pval = pvfunc(tau, n, niv, itv);

    return pval;
}

private double get_mackinnon_pvalue(AdfInfo * ainfo)
{
    int asy = (ainfo.kmax > 0 || ainfo.order > 0);
    int niv, T = asy ? 0 : ainfo.T;
    double pval;

    if (ainfo.flags & AdfFlags.eg_resids) {
        niv = ainfo.niv;
    } else {
        niv = 1;
    }

    pval = get_urc_pvalue(ainfo.tau, T, niv, ainfo.det);
    if (!na(pval)) {
        ainfo.pvtype = asy ? PvalueType.asy : PvalueType.finite;
    }

    return pval;
}

private double get_dfgls_pvalue(AdfInfo * ainfo)
{
    alias DfglsFunc = extern (C) double function(double, int, int, int, int, int *);
    DfglsFunc pvfunc;
    double pval = NADBL;
    int T = ainfo.T;
    int PQ, trend, err = 0;

    pvfunc = cast(DfglsFunc) get_plugin_function("dfgls_pvalue");
    if (pvfunc == null) {
        return pval;
    }

    if (ainfo.kmax == 0 && ainfo.order > 0) {
        /* fixed lag order > 0: asymptotic */
        T = 0;
    }
    trend = (ainfo.det == DetCode.trend);
    PQ = (ainfo.flags & AdfFlags.pq) ? 1 : 0;

    pval = pvfunc(ainfo.tau, T, trend, ainfo.kmax, PQ, &err);
    if (!na(pval)) {
        ainfo.pvtype = (T == 0) ? PvalueType.asy : PvalueType.approx;
    }

    return pval;
}

private bool test_opt_not_set(Gretlopt o)
{
    return (!(o & Gretlopt.OPT_N) && !(o & Gretlopt.OPT_C) &&
            !(o & Gretlopt.OPT_T) && !(o & Gretlopt.OPT_R));
}

private int test_wanted(DetCode det, Gretlopt opt)
{
    int ret = 0;

    switch (det) {
    case DetCode.no_const:
        ret = (opt & Gretlopt.OPT_N);
        break;
    case DetCode.const_:
        ret = (opt & Gretlopt.OPT_C);
        break;
    case DetCode.trend:
        ret = (opt & Gretlopt.OPT_T);
        break;
    case DetCode.quad_trend:
        ret = (opt & Gretlopt.OPT_R);
        break;
    default:
        break;
    }

    return ret;
}

private DetCode engle_granger_itv(Gretlopt opt)
{
    DetCode itv = DetCode.const_;

    if (opt & Gretlopt.OPT_N) {
        itv = DetCode.no_const;
    } else if (opt & Gretlopt.OPT_T) {
        itv = DetCode.trend;
    } else if (opt & Gretlopt.OPT_R) {
        itv = DetCode.quad_trend;
    }

    return itv;
}

private int gettrend(DatasetRaw * dset, int square)
{
    int idx, t, v = dset.v;
    double x;

    idx = series_index(dset, (square) ? "timesq" : "time");

    if (idx < v) {
        return idx;
    }

    if (dataset_add_series(dset, 1)) {
        return 0; /* error: valid value cannot == 0 */
    }

    for (t=0; t<dset.n; t++) {
        x = cast(double) t + 1;
        dset.Z[v][t] = (square) ? x * x : x;
    }

    if (square) {
        strcpy(dset.varname[v], "timesq");
        series_set_label(dset, v, _("squared time trend variable").ptr);
    } else {
        strcpy(dset.varname[v], "time");
        series_set_label(dset, v, _("time trend variable").ptr);
    }

    return idx;
}

private void print_df_model(AdfInfo * ainfo, Model * pmod,
                            int dfnum, DatasetRaw * dset,
                            Prn * prn)
{
    pmod.aux = (ainfo.order > 0) ? ModelAuxCode.adf : ModelAuxCode.df;

    if (!na(ainfo.pval)) {
        gretl_model_set_int(pmod, "dfnum", dfnum);
        gretl_model_set_double(pmod, "dfpval", ainfo.pval);
    }

    if (ainfo.flags & AdfFlags.eg_resids) {
        gretl_model_set_int(pmod, "eg-resids", 1);
    }

    if (ainfo.g != null) {
        gretl_model_set_int(pmod, "dfgls", 1);
    }

    printmodel(pmod, dset, Gretlopt.OPT_NONE, prn);

    if (ainfo.g != null) {
        pputs(prn, "  ");
        if (ainfo.g.rows == 2) {
            pprintf(prn, _("GLS detrending: b0 = %g, b1 = %g\n").ptr,
                    ainfo.g.val[0], ainfo.g.val[1]);
        } else {
            pprintf(prn, _("GLS estimate of b0: %g\n").ptr,
                    ainfo.g.val[0]);
        }
        pputc(prn, '\n');
    }
}

private int set_deterministic_terms(AdfInfo * ainfo,
                                    DatasetRaw * dset)
{
    int i, j;

    /* Note that list[1] and list[2], plus the @order
       lagged differences, are in common for all
       specifications. So we start adding deterministic
       terms at position 3 + ainfo->order.
    */

    ainfo.list[0] = 1 + ainfo.order + ainfo.det + ainfo.nseas;
    i = 3 + ainfo.order;

    if (ainfo.det >= DetCode.trend) {
        ainfo.list[i] = gettrend(dset, 0);
        if (ainfo.list[i] == 0) {
            return GretlError.alloc;
        }
        i++;
    }

    if (ainfo.det == DetCode.quad_trend) {
        ainfo.list[i] = gettrend(dset, 1);
        if (ainfo.list[i] == 0) {
            return GretlError.alloc;
        }
        i++;
    }

    if (ainfo.nseas > 0) {
        for (j=0; j<ainfo.nseas; j++) {
            ainfo.list[i++] = ainfo.slist[j+1];
        }
    }

    if (ainfo.det != DetCode.no_const) {
        /* stick constant on end of list */
        ainfo.list[i] = 0;
    }

    return 0;
}

/* When we're done with testing down in the DF-GLS context
   we may need to regenerate the detrended data: this
   applies if (a) we used Perron-Qu OLS detrending in the
   test-down phase or (b) the sample start is currently
   set after the start of the dataset.
*/

private int reset_detrended_data(AdfInfo * ainfo,
                                 DatasetRaw * dset)
{
    if (ainfo.flags & AdfFlags.pq) {
        /* testing down using Perron-Qu */
        return 1;
    } else if ((ainfo.flags & AdfFlags.gls) && dset.t1 > 0) {
        /* GLS + testing down + t1 > 0 */
        return 1;
    } else {
        return 0;
    }
}

private int handle_test_down_option(AdfInfo * ainfo,
                                    Gretlopt opt,
                                    int * err)
{
    int kmethod = 0;
    const(char) * s;

    if (ainfo.flags & (AdfFlags.eg_test | AdfFlags.eg_resids)) {
        s = get_optval_string(GretlCmdIndex.coint, Gretlopt.OPT_E);
    } else {
        s = get_optval_string(GretlCmdIndex.adf, Gretlopt.OPT_E);
    }

    if (s == null || *s == '\0') {
        /* the default */
        kmethod = LagSelection.aic;
    } else if (!strcmp(s, "MAIC") || !strcmp(s, "AIC")) {
        kmethod = LagSelection.aic;
    } else if (!strcmp(s, "MBIC") || !strcmp(s, "BIC")) {
        kmethod = LagSelection.bic;
    } else if (!strcmp(s, "tstat")) {
        kmethod = LagSelection.tstat;
    } else {
        gretl_errmsg_set(_("Invalid option").ptr);
        *err = GretlError.data;
    }

    if (!*err) {
        /* take the given order to be the max */
        ainfo.kmax = ainfo.order;
        if (ainfo.flags & AdfFlags.pq) {
            /* --perron-qu */
            if (kmethod != LagSelection.aic && kmethod != LagSelection.bic) {
                *err = GretlError.badopt;
            }
        }
    }

    if (kmethod > 0) {
        ainfo.flags |= AdfFlags.test_down;
    }

    return kmethod;
}

private int check_adf_options(AdfInfo * ainfo, Gretlopt opt)
{
    int err = 0;

    if (opt & Gretlopt.OPT_G) {
        /* we can only have one of the basic deterministics options */
        err = incompatible_options(opt, Gretlopt.OPT_N | Gretlopt.OPT_C | Gretlopt.OPT_T | Gretlopt.OPT_R);
        /* options incompatible with --gls: no-const, seasonals,
           and quadratic trend
        */
        if (!err && (opt & (Gretlopt.OPT_N | Gretlopt.OPT_D | Gretlopt.OPT_R))) {
            err = GretlError.badopt;
        }
        if (!err) {
            /* --quiet and --verbose are icompatible */
            err = incompatible_options(opt, Gretlopt.OPT_Q | Gretlopt.OPT_V);
        }
        if (!err) {
            ainfo.flags |= AdfFlags.gls;
            if (opt & Gretlopt.OPT_U) {
                /* Perron-Qu */
                ainfo.flags |= AdfFlags.pq;
            }
        }
    } else if (opt & Gretlopt.OPT_U) {
        /* option dependent on --gls: Perron-Qu modified AIC/BIC */
        err = GretlError.badopt;
    }

    return err;
}

private int real_adf_test(AdfInfo * ainfo, DatasetRaw * dset,
                          Gretlopt opt, Prn * prn)
{
    Model dfmod;
    Gretlopt eg_opt = Gretlopt.OPT_NONE;
    Gretlopt df_mod_opt;
    int * biglist = null;
    int orig_nvars = dset.v;
    int blurb_done = 0;
    int demean_mode = 0;
    int test_down = 0;
    int test_num = 0;
    int i, err;

    df_mod_opt = ur_lsq_opt(dset, opt);

    /* (most of) this may have been done already
       but it won't hurt to check here */
    err = check_adf_options(ainfo, opt);
    if (err) {
        return err;
    }

    /* safety-first initializations */
    ainfo.nseas = ainfo.kmax = ainfo.altv = 0;
    ainfo.vname = dset.varname[ainfo.v];
    ainfo.list = ainfo.slist = null;
    ainfo.det = cast(DetCode) 0;

    if (gretl_isconst(dset.t1, dset.t2, dset.Z[ainfo.v])) {
        gretl_errmsg_sprintf(_("%s is a constant").ptr, ainfo.vname);
        return GretlError.data;
    }

    if (opt & Gretlopt.OPT_E) {
        /* --test-down[=...] */
        test_down = handle_test_down_option(ainfo, opt, &err);
        if (err) {
            return err;
        }
    }

    if (ainfo.order < 0) {
        test_down = LagSelection.aic;
        ainfo.order = ainfo.kmax = -ainfo.order;
    }

    if (ainfo.flags & AdfFlags.eg_resids) {
        /* Final step of Engle-Granger test: the (A)DF test
           regression will contain no deterministic terms,
           but the selection of the p-value is based on the
           deterministic terms in the cointegrating
           regression, represented by @eg_opt.
        */
        int verbose = (opt & Gretlopt.OPT_V);
        int silent = (opt & Gretlopt.OPT_I);

        eg_opt = opt;
        opt = Gretlopt.OPT_N;
        if (silent) {
            opt |= Gretlopt.OPT_I;
        } else if (verbose) {
            opt |= Gretlopt.OPT_V;
        }
    }

    if ((opt & Gretlopt.OPT_D) && dset.pd > 1) {
        /* arrange to add seasonal dummies */
        ainfo.nseas = dset.pd - 1;
    }

    if (test_opt_not_set(opt)) {
        /* default model(s) */
        if (ainfo.flags & AdfFlags.gls) {
            opt |= Gretlopt.OPT_C;
        } else {
            opt |= (Gretlopt.OPT_C | Gretlopt.OPT_T);
        }
    }

    if (ainfo.flags & AdfFlags.pq) {
        demean_mode = DemeanMethod.ols;
    } else if (ainfo.flags & AdfFlags.gls) {
        demean_mode = DemeanMethod.gls;
    }

    err = adf_prepare_vars(ainfo, dset, opt, demean_mode, prn);
    if (err) {
        return err;
    }

    if (test_down) {
        ainfo.list[0] = ainfo.order + 5;
        biglist = gretl_list_copy(ainfo.list);
        if (biglist == null) {
            err = GretlError.alloc;
            goto bailout;
        }
    }

    gretl_model_init(&dfmod, dset);

    /* Now loop across the wanted deterministics cases:
       in many instances we'll actually be doing only one
       case.
    */

    for (i=DetCode.no_const; i<DetCode.max; i++) {
        int b0pos = (i > DetCode.no_const);

        ainfo.det = cast(DetCode) i;
        if (!test_wanted(ainfo.det, opt)) {
            continue;
        }

        if (test_down) {
            /* re-establish max order before testing down */
            ainfo.order = ainfo.kmax;
            copy_list_values(ainfo.list, biglist);
        }

        if (ainfo.flags & AdfFlags.gls) {
            /* skip deterministics */
            ainfo.list[0] = ainfo.order + 2;
            b0pos = 0;
        } else {
            err = set_deterministic_terms(ainfo, dset);
            if (err) {
                goto bailout;
            }
        }

        test_num++;

        if (test_down) {
            /* determine the optimal lag order */
            if (test_down == LagSelection.tstat) {
                ainfo.order = t_adjust_order(ainfo, dset, opt,
                                             &err, prn);
            } else {
                ainfo.order = ic_adjust_order(ainfo, cast(LagSelection) test_down,
                                              dset, opt, test_num,
                                              &err, prn);
            }
            if (err) {
                clear_model(&dfmod);
                goto bailout;
            } else if (reset_detrended_data(ainfo, dset)) {
                /* swap out the detrended data */
                err = adf_prepare_vars(ainfo, dset, opt, DemeanMethod.gls, prn);
            }
        }

        /* run the actual test regression */
        dfmod = lsq(ainfo.list, dset, GretlCmdIndex.ols, df_mod_opt);

        if (!dfmod.errcode && dfmod.dfd == 0) {
            /* we can't tolerate an exact fit here */
            dfmod.errcode = GretlError.df;
        }
        if (dfmod.errcode) {
            fprintf(stderr, "adf_test: dfmod.errcode = %d\n",
                    dfmod.errcode);
            err = dfmod.errcode;
            clear_model(&dfmod);
            goto bailout;
        }

        /* transcribe info from test regression */
        ainfo.b0 = dfmod.coeff[b0pos];
        ainfo.tau = ainfo.b0 / dfmod.sderr[b0pos];
        ainfo.T = dfmod.nobs;
        ainfo.df = dfmod.dfd;

        if (ainfo.flags & AdfFlags.eg_resids) {
            ainfo.det = engle_granger_itv(eg_opt);
        }

        /* obtain p-value if wanted */
        if (getenv("DFGLS_NO_PVALUE")) {
            /* skip it, to speed up monte carlo stuff */
            ainfo.pval = NADBL;
        } else if (ainfo.flags & AdfFlags.gls) {
            /* use Cottrell/Komashko or Sephton */
            ainfo.pval = get_dfgls_pvalue(ainfo);
        } else {
            /* no GLS adjustment applied: use MacKinnon p-values,
               either finite sample or asymptotic in the case
               of an augmented test
            */
            ainfo.pval = get_mackinnon_pvalue(ainfo);
        }

        if (!(opt & (Gretlopt.OPT_Q | Gretlopt.OPT_I)) && !(ainfo.flags & AdfFlags.panel)) {
            print_adf_results(ainfo, &dfmod, &blurb_done,
                              opt, test_down, prn);
        }

        if ((opt & Gretlopt.OPT_V) && !(ainfo.flags & AdfFlags.panel)) {
            /* verbose */
            print_df_model(ainfo, &dfmod, b0pos, dset, prn);
        } else if (!(opt & Gretlopt.OPT_Q) && !(ainfo.flags & (AdfFlags.eg_resids | AdfFlags.panel))) {
            pputc(prn, '\n');
        }

        clear_model(&dfmod);
    }

    if (!err) {
        if (!(ainfo.flags & (AdfFlags.eg_test | AdfFlags.panel)) ||
            (ainfo.flags & AdfFlags.eg_resids)) {
            record_test_result(ainfo.tau, ainfo.pval);
        }
    }

 bailout:

    free(ainfo.list);
    ainfo.list = null;
    free(ainfo.slist);
    ainfo.slist = null;
    gretl_matrix_free(ainfo.g);
    ainfo.g = null;
    free(biglist);

    dataset_drop_last_variables(dset, dset.v - orig_nvars);

    return err;
}

/* print critical value for ADF-IPS or KPSS */

private void print_critical_values(double * a, double * cv,
                                   GretlCmdIndex ci, Prn * prn)
{
    string label = N_("Critical values");
    int figs = (ci == GretlCmdIndex.adf) ? 2 : 3;

    pprintf(prn, "%*s    ", cast(int) TRANSLATED_WIDTH(_(label)), " ".ptr);
    pprintf(prn, "%g%%      %g%%      %g%%\n",
            100*a[0], 100*a[1], 100*a[2]);
    pprintf(prn, "%s: %.*f   %.*f   %.*f\n",
            _(label).ptr, figs, cv[0], figs, cv[1], figs, cv[2]);
}

private int panel_adjust_ADF_opt(Gretlopt * opt)
{
    int err = 0;

    /* Has the user selected an option governing the
       deterministic terms to be included? If so, don't
       mess with it.
    */
    if (*opt & (Gretlopt.OPT_N | Gretlopt.OPT_C | Gretlopt.OPT_R | Gretlopt.OPT_T)) {
        ; /* no-op */
    } else {
        /* panel default: test with constant */
        *opt |= Gretlopt.OPT_C;
    }

    return err;
}

private int DF_index(Gretlopt opt)
{
    if (opt & Gretlopt.OPT_N) {
        return 0;
    } else if (opt & Gretlopt.OPT_C) {
        return 1;
    } else if (opt & Gretlopt.OPT_T) {
        return 2;
    } else {
        return 3;
    }
}

/* See Im, Pesaran and Shin, "Testing for unit roots in
   heterogeneous panels", Journal of Econometrics 115 (2003),
   53-74.
*/

private int do_IPS_test(double tbar, int n, const(int) * Ti,
                        int order, const(int) * Oi,
                        Gretlopt opt, Prn * prn)
{
    alias CritvalsFunc = extern (C) int function(int, int, int, double *);
    alias MomentsFunc = extern (C) int function(int, double *, double *);
    alias RhoMomentsFunc = extern (C) int function(int, int, int, double *, double *);
    CritvalsFunc get_IPS_critvals;
    MomentsFunc IPS_tbar_moments;
    RhoMomentsFunc IPS_tbar_rho_moments;
    int Tmin = Ti[1], Tmax = Ti[1];
    int i, T, err = 0;

    for (i=2; i<=n; i++) {
        if (Ti[i] > Tmax) {
            Tmax = Ti[i];
        }
        if (Ti[i] < Tmin) {
            Tmin = Ti[i];
        }
    }

    if (Oi != null || order > 0) {
        /* non-zero lag order: use IPS's W_{tbar} statistic */
        double E, V, Wtbar, Etbar = 0, Vtbar = 0;
        int order_i;

        IPS_tbar_rho_moments = cast(RhoMomentsFunc) get_plugin_function("IPS_tbar_rho_moments");

        if (IPS_tbar_rho_moments != null) {
            for (i=0; i<n && !err; i++) {
                T = Ti[i+1];
                order_i = (Oi != null) ? Oi[i+1] : order;
                err = IPS_tbar_rho_moments(order_i, T, (opt & Gretlopt.OPT_T), &E, &V);
                Etbar += E;
                Vtbar += V;
            }
            if (!err) {
                Etbar /= n;
                Vtbar /= n;
                Wtbar = sqrt(cast(double) n) * (tbar - Etbar) / sqrt(Vtbar);
                pprintf(prn, "N = %d, Tmin = %d, Tmax = %d\n", n, Tmin, Tmax);
                pprintf(prn, "Im-Pesaran-Shin W_tbar = %g [%.4f]\n", Wtbar,
                        normal_pvalue_1(-Wtbar));
            } else {
                pputs(prn, _("Im-Pesaran-Shin test: cannot be computed\n").ptr);
            }
        }
    } else if (Tmax > Tmin) {
        /* sample sizes differ: use IPS's Z_{tbar} */
        double E, V, Ztbar, Etbar = 0, Vtbar = 0;

        IPS_tbar_moments = cast(MomentsFunc) get_plugin_function("IPS_tbar_moments");

        if (IPS_tbar_moments != null) {
            for (i=0; i<n && !err; i++) {
                T = Ti[i+1];
                err = IPS_tbar_moments(T, &E, &V);
                Etbar += E;
                Vtbar += V;
            }
            if (!err) {
                Etbar /= n;
                Vtbar /= n;
                Ztbar = sqrt(cast(double) n) * (tbar - Etbar) / sqrt(Vtbar);
                pprintf(prn, "N = %d, Tmin = %d, Tmax = %d\n", n, Tmin, Tmax);
                pprintf(prn, "Im-Pesaran-Shin Z_tbar = %g [%.4f]\n", Ztbar,
                        normal_pvalue_1(-Ztbar));
            } else {
                pputs(prn, _("Im-Pesaran-Shin test: cannot be computed\n").ptr);
            }
        }
    } else {
        /* simple case: use tbar with exact critical values */
        pprintf(prn, "N,T = (%d,%d)\n", n, Tmax);
        pprintf(prn, "Im-Pesaran-Shin t-bar = %g\n", tbar);

        get_IPS_critvals = cast(CritvalsFunc) get_plugin_function("get_IPS_critvals");

        if (get_IPS_critvals != null) {
            double[3] a = [0.1, 0.05, 0.01];
            double[3] cv;

            err = get_IPS_critvals(n, Tmax, (opt & Gretlopt.OPT_T), cv.ptr);
            if (!err) {
                print_critical_values(a.ptr, cv.ptr, GretlCmdIndex.adf, prn);
            }
        } else {
            pputs(prn, _("Im-Pesaran-Shin test: cannot be computed\n").ptr);
        }
    }

    return err;
}

private void record_choi_test(double * ct, double * cp)
{
    GretlMatrixRaw * tests = gretl_column_vector_alloc(3);
    GretlMatrixRaw * pvals = gretl_column_vector_alloc(3);
    char ** S = strings_array_new(3);
    int i;

    if (tests != null && pvals != null) {
        for (i=0; i<3; i++) {
            tests.val[i] = ct[i];
            pvals.val[i] = cp[i];
        }
        if (S != null) {
            S[0] = gretl_strdup(_("Inverse chi-square").ptr);
            S[1] = gretl_strdup(_("Inverse normal test").ptr);
            S[2] = gretl_strdup(_("Logit test").ptr);
            gretl_matrix_set_rownames(tests, S);
        }
        record_matrix_test_result(tests, pvals);
    } else {
        gretl_matrix_free(tests);
        gretl_matrix_free(pvals);
        strings_array_free(S, 3);
    }
}

/* See In Choi, "Unit root tests for panel data", Journal of
   International Money and Finance 20 (2001), 249-272.
*/

private void do_choi_test(double ppv, double zpv, double lpv,
                          int n, Prn * prn)
{
    double[3] choi_test;
    double[3] choi_pval;
    int tdf = 5 * n + 4;
    double k = (3.0*tdf)/(M_PI*M_PI*n*(5*n+2));

    /* inverse chi^2 */
    choi_test[0] = -2 * ppv;
    choi_pval[0] = chisq_cdf_comp(2*n, choi_test[0]);

    /* inverse normal */
    choi_test[1] = zpv / sqrt(cast(double) n);
    choi_pval[1] = normal_pvalue_1(-choi_test[1]);

    /* logit */
    choi_test[2] = sqrt(k) * lpv;
    choi_pval[2] = student_pvalue_1(tdf, -choi_test[2]);

    if (prn != null) {
        pprintf(prn, "%s\n", _("Choi meta-tests:").ptr);
        pprintf(prn, "   %s(%d) = %g [%.4f]\n", _("Inverse chi-square").ptr,
                2*n, choi_test[0], choi_pval[0]);
        pprintf(prn, "   %s = %g [%.4f]\n", _("Inverse normal test").ptr,
                choi_test[1], choi_pval[1]);
        pprintf(prn, "   %s: t(%d) = %g [%.4f]\n", _("Logit test").ptr,
                tdf, choi_test[2], choi_pval[2]);
    }

    record_choi_test(choi_test.ptr, choi_pval.ptr);
}

private void panel_unit_DF_print(AdfInfo * ainfo, int i, Prn * prn)
{
    if (ainfo.flags & AdfFlags.test_down) {
        pputc(prn, '\n');
    }
    pprintf(prn, "%s %d, T = %d, %s = %d\n", _("Unit").ptr, i,
            ainfo.T, _("lag order").ptr, ainfo.order);
    pprintf(prn, "   %s: %g\n"
            ~ "   %s = %g",
            _("estimated value of (a - 1)").ptr, ainfo.b0,
            _("test statistic").ptr, ainfo.tau);
    if (na(ainfo.pval)) {
        pputs(prn, "\n\n");
    } else {
        pprintf(prn, " [%.4f]\n\n", ainfo.pval);
    }
}

/* When the test order is given as -1 we default to L_{12}: see
   G. W. Schwert, "Tests for Unit Roots: A Monte Carlo Investigation",
   Journal of Business and Economic Statistics, 7(2), 1989,
   pp. 5-17. However, in the case of panel data with a short time-series
   dimension we somewhat arbitrarily cap the lag order at 5 in an effort
   to enable the Im-Pesaran-Shinn test.
*/

private int auto_lag_order(const(DatasetRaw) * dset)
{
    int T = sample_size(dset);
    int ret;

    ret = cast(int) floor(12.0 * pow(T/100.0, 0.25));
    if (dataset_is_panel(dset) && T < 24) {
        ret = min(ret, 5);
    }

    return ret;
}

private void set_ainfo_verbosity(AdfInfo * ainfo,
                                 Gretlopt opt)
{
    ainfo.verbosity = 0;

    if (opt & Gretlopt.OPT_V) {
        const(char) * s = get_optval_string(GretlCmdIndex.adf, Gretlopt.OPT_V);

        if (s != null && !strcmp(s, "max")) {
            ainfo.verbosity = 2;
        } else {
            ainfo.verbosity = 1;
        }
    }
}

private int panel_DF_test(int v, int order, DatasetRaw * dset,
                          Gretlopt opt, Prn * prn)
{
    AdfInfo ainfo;
    int u0 = dset.t1 / dset.pd;
    int uN = dset.t2 / dset.pd;
    int quiet = (opt & Gretlopt.OPT_Q);
    int verbose = (opt & Gretlopt.OPT_V);
    double ppv = 0.0, zpv = 0.0, lpv = 0.0;
    double pval, tbar = 0.0;
    int * Ti = null, Oi = null;
    int autolag = 0;
    int i, n, err;

    err = panel_adjust_ADF_opt(&opt);
    if (err) {
        return err;
    }

    if (opt & Gretlopt.OPT_G) {
        /* GLS option: can't do IPS t-bar test */
        tbar = NADBL;
    } else {
        Ti = gretl_list_new(uN - u0 + 1);
        if (Ti == null) {
            return GretlError.alloc;
        }
        if ((opt & Gretlopt.OPT_E) || order < 0) {
            /* testing down: lag order may vary by unit */
            Oi = gretl_list_new(uN - u0 + 1);
            if (Oi == null) {
                free(Ti);
                return GretlError.alloc;
            }
        }
    }

    if (!quiet) {
        int j = DF_index(opt);

        DF_header(dset.varname[v], (opt & Gretlopt.OPT_E) ? 0 : order,
                  0, 0, opt, prn);
        pprintf(prn, "   %s ", _(DF_test_string(j)).ptr);
        pputc(prn, '\n');
        pprintf(prn, "   %s: %s\n\n", _("model").ptr,
                ((order > 0) ? ADF_model_string(j) : DF_model_string(j)).ptr);
    }

    /* number of units in sample range */
    n = uN - u0 + 1;

    /* initialize @ainfo */
    ainfo.v = v;
    ainfo.niv = 1;
    ainfo.flags = AdfFlags.panel;
    set_ainfo_verbosity(&ainfo, opt);

    if (order == -1) {
        /* automatic order selection */
        autolag = 1;
    }

    /* run a Dickey-Fuller test for each unit and record the
       results */

    for (i=u0; i<=uN && !err; i++) {
        ainfo.order = order;
        dset.t1 = i * dset.pd;
        dset.t2 = dset.t1 + dset.pd - 1;
        err = series_adjust_sample(dset.Z[v], &dset.t1, &dset.t2);
        if (!err) {
            if (autolag) {
                ainfo.order = auto_lag_order(dset);
            }
            err = real_adf_test(&ainfo, dset, opt, prn);
            if (!err && verbose) {
                panel_unit_DF_print(&ainfo, i+1, prn);
            }
        }
        if (!err) {
            if (Ti != null) {
                Ti[i-u0+1] = ainfo.T;
            }
            if (Oi != null) {
                Oi[i-u0+1] = ainfo.order;
            }
            if (!na(tbar)) {
                tbar += ainfo.tau;
            }
            pval = ainfo.pval;
            if (na(pval)) {
                ppv = zpv = lpv = NADBL;
            } else if (!na(ppv)) {
                ppv += log(pval);
                zpv += normal_cdf_inverse(pval);
                lpv += log(pval / (1-pval));
            }
        }
    }

    /* process the results as per Im-Pesaran-Shin and/or Choi */

    if (!err) {
        if (!quiet) {
            pprintf(prn, "%s\n\n", _("H0: all groups have unit root").ptr);
        }
        if (!quiet && !na(tbar)) {
            tbar /= n;
            do_IPS_test(tbar, n, Ti, order, Oi, opt, prn);
        }
        if (!na(ppv)) {
            if (quiet) {
                do_choi_test(ppv, zpv, lpv, n, null);
            } else {
                pputc(prn, '\n');
                do_choi_test(ppv, zpv, lpv, n, prn);
            }
        }
        if (!quiet) {
            pputc(prn, '\n');
        }
    }

    free(Ti);
    free(Oi);

    return err;
}

/*
 * levin_lin_test:
 * @vnum: ID number of variable to test.
 * @plist: list of ADF lag orders.
 * @dset: data information struct.
 * @opt: option flags.
 * @prn: gretl printing struct.
 *
 * Carries out and prints the results of the Levin-Lin-Chu test
 * for a unit root in panel data.
 *
 * The list @plist should contain either a single lag order
 * to be applied to all units, or a set of unit-specific
 * orders; in the latter case the length of the list must
 * equal the number of panel units in the current sample
 * range. (This is a gretl list: the first element holds
 * a count of the number of elements following.)
 *
 * By default a test with constant is performed, but the
 * (mutually exclusive) options OPT_N and OPT_T in @opt switch to
 * the case of no constant or constant plus trend respectively.
 * The OPT_Q flag may be used to suppress printed output.
 *
 * Returns: 0 on successful completion, non-zero on error.
 */

int levin_lin_test(int vnum, const(int) * plist,
                   DatasetRaw * dset, Gretlopt opt,
                   Prn * prn)
{
    alias LevinLinFunc = extern (C) int function(int, const(int) *, DatasetRaw *,
                                                 Gretlopt, Prn *);
    LevinLinFunc real_levin_lin;
    int panelmode;
    int err = 0;

    panelmode = multi_unit_panel_sample(dset);

    if (!panelmode || incompatible_options(opt, Gretlopt.OPT_N | Gretlopt.OPT_T)) {
        return GretlError.badopt;
    }

    real_levin_lin = cast(LevinLinFunc) get_plugin_function("real_levin_lin");

    if (real_levin_lin == null) {
        err = GretlError.fopen;
    } else {
        int save_t1 = dset.t1;
        int save_t2 = dset.t2;

        err = real_levin_lin(vnum, plist, dset, opt, prn);

        dset.t1 = save_t1;
        dset.t2 = save_t2;
    }

    return err;
}

/* Get the ID number of the series to be tested, allowing for the
   possibility that the --difference option (OPT_F) has been specified.
*/

private int get_real_v(int v0, DatasetRaw * dset,
                       Gretlopt opt, int * err)
{
    int v = v0;

    if (opt & Gretlopt.OPT_F) {
        /* --difference */
        int save_t1 = dset.t1;

        dset.t1 = 0;
        v = diffgenr(v0, GretlCmdIndex.diff, dset);
        dset.t1 = save_t1;
        if (v < 0) {
            *err = GretlError.data; /* ? */
        }
    }

    return v;
}

/*
 * adf_test:
 * @order: lag order for the (augmented) test.
 * @list: list of variables to test.
 * @dset: dataset struct.
 * @opt: option flags.
 * @prn: gretl printing struct.
 *
 * Carries out and prints the results of the Augmented Dickey-Fuller
 * test for a unit root.
 *
 * By default two tests are performed, one for a model
 * including a constant and one including a linear trend. The
 * deterministic components of the model can be controlled via
 * flags in @opt as follows: OPT_N, omit the constant; OPT_C,
 * run just one test using the constant; OPT_T, one test including
 * linear trend; OPT_R, one test including a quadratic trend;
 * OPT_D, include seasonal dummy variables.
 *
 * Additional flags that may be given in @opt include:
 * OPT_V for verbose operation; OPT_F to apply first-differencing
 * before testing; OPT_G for GLS preprocessing as in Elliott, Rothenberg
 * and Stock (incompatible with OPT_N, OPT_R, OPT_D); OPT_E to
 * "test down" from a given maximum lag order (see the entry for
 * "adf" in the Gretl Command Reference for details).
 *
 * Returns: 0 on successful completion, non-zero on error.
 */

int adf_test(int order, const(int) * list, DatasetRaw * dset,
             Gretlopt opt, Prn * prn)
{
    int save_t1 = dset.t1;
    int save_t2 = dset.t2;
    int orig_v = dset.v;
    int panelmode;
    int err = 0;

    if (opt & Gretlopt.OPT_G) {
        /* GLS is incompatible with no const (OPT_N), quadratic
           trend (OPT_R) and seasonals (OPT_D) */
        err = options_incompatible_with(opt, Gretlopt.OPT_G,
                                        Gretlopt.OPT_N | Gretlopt.OPT_R | Gretlopt.OPT_D);
        if (!err) {
            /* and you have to choose between OPT_C nd OPT_T */
            err = incompatible_options(opt, Gretlopt.OPT_C | Gretlopt.OPT_T);
        }
    }

    panelmode = multi_unit_panel_sample(dset);

    if (panelmode) {
        int v = get_real_v(list[1], dset, opt, &err);

        if (!err) {
            err = panel_DF_test(v, order, dset, opt, prn);
        }
    } else {
        /* the regular time series case */
        int i;
        int[2] vlist = [1, 0];
        AdfInfo ainfo;

        ainfo.niv = 1;
        set_ainfo_verbosity(&ainfo, opt);

        for (i=1; i<=list[0] && !err; i++) {
            ainfo.v = vlist[1] = get_real_v(list[i], dset, opt, &err);
            if (err) {
                break;
            }
            ainfo.order = order;
            err = list_adjust_sample(vlist.ptr, &dset.t1, &dset.t2, dset, null);
            if (!err && order == -1) {
                ainfo.order = auto_lag_order(dset);
            }
            if (!err) {
                err = real_adf_test(&ainfo, dset, opt, prn);
            }
            dset.t1 = save_t1;
            dset.t2 = save_t2;
        }
    }

    dset.t1 = save_t1;
    dset.t2 = save_t2;

    if (dset.v > orig_v) {
        dataset_drop_last_variables(dset, dset.v - orig_v);
    }

    return err;
}

/* See Peter S. Sephton, "Response surface estimates of the KPSS
   stationarity test", Economics Letters 47 (1995) 255-261.

   The estimates below of \beta_\infty and \beta_1 (based on
   a bigger replication using Sephton's methodology) allow the
   construction of better critical values for finite samples
   than the values given in the original KPSS article.
*/

private void kpss_parms(double a, int trend, double * b)
{
    static immutable double[3] b0_level = [ 0.74404, 0.46158, 0.34742 ];
    static immutable double[3] b1_level = [ -0.99120, 0.01642, 0.19814 ];
    static immutable double[3] b0_trend = [ 0.21787, 0.14797, 0.11925 ];
    static immutable double[3] b1_trend = [ -0.25128, 0.03270, 0.10244 ];
    int i = (a == .01) ? 0 : (a == .05) ? 1 : 2;

    if (trend) {
        b[0] = b0_trend[i];
        b[1] = b1_trend[i];
    } else {
        b[0] = b0_level[i];
        b[1] = b1_level[i];
    }
}

private double kpss_critval(double alpha, int T, int trend)
{
    double[2] b;

    kpss_parms(alpha, trend, b.ptr);

    return b[0] + b[1]/T;
}

GretlMatrixRaw * kpss_critvals(int T, int trend, int * err)
{
    GretlMatrixRaw * m = null;

    if (T < 5) {
        *err = GretlError.toofew;
    } else {
        m = gretl_matrix_alloc(1, 3);
        if (m == null) {
            *err = GretlError.alloc;
        } else {
            m.val[0] = kpss_critval(0.10, T, trend);
            m.val[1] = kpss_critval(0.05, T, trend);
            m.val[2] = kpss_critval(0.01, T, trend);
        }
    }

    return m;
}

enum PV_GT10 = 1.1;
enum PV_LT01 = -1.0;

private double kpss_interp(double s, int T, int trend)
{
    double c10, c05, c01;
    double pv;

    c10 = kpss_critval(.10, T, trend);
    if (s < c10) {
        return PV_GT10;
    }

    c01 = kpss_critval(.01, T, trend);
    if (s > c01) {
        return PV_LT01;
    }

    /* OK, p-value must lie between .01 and .10 */

    c05 = kpss_critval(.05, T, trend);
    if (s > c05) {
        pv = .01 + .04 * (c01 - s) / (c01 - c05);
    } else {
        pv = .05 + .05 * (c05 - s) / (c05 - c10);
    }

    return pv;
}

private int
real_kpss_test(int order, int varno, DatasetRaw * dset,
               Gretlopt opt, KpssInfo * kinfo,
               Prn * prn)
{
    Model KPSSmod;
    int * list = null;
    int hastrend = 0, hasseas = 0;
    double et, s2 = 0.0;
    double cumsum = 0.0, cumsum2 = 0.0;
    double teststat, pval = NADBL;
    double * autocov;
    int t1, t2, T;
    int i, t, ndum, nreg;
    int err = 0;

    /* sanity check */
    if (varno <= 0 || varno >= dset.v) {
        return GretlError.data;
    }

    if (gretl_isconst(dset.t1, dset.t2, dset.Z[varno])) {
        gretl_errmsg_sprintf(_("%s is a constant").ptr, dset.varname[varno]);
        return GretlError.data;
    }

    if (opt & Gretlopt.OPT_F) {
        /* difference the variable before testing */
        varno = diffgenr(varno, GretlCmdIndex.diff, dset);
        if (varno < 0) {
            return GretlError.data;
        }
    }

    hastrend = (opt & Gretlopt.OPT_T) ? 1 : 0;
    hasseas = (opt & Gretlopt.OPT_D) ? 1 : 0;
    ndum = hasseas ? (dset.pd - 1) : 0;
    nreg = 1 + hastrend + ndum;

    list = gretl_list_new(nreg + 1);
    if (list == null) {
        return GretlError.alloc;
    }

    list[1] = varno;
    list[2] = 0;

    if (hastrend) {
        list[3] = gettrend(dset, 0);
        if (list[3] == 0) {
            return GretlError.alloc;
        }
    }

    if (hasseas) {
        int * slist = null;

        slist = seasonals_list(dset, dset.pd, 1, &err);
        if (err) {
            free(list);
            return err;
        } else {
            for (i=0; i<ndum; i++) {
                list[3 + hastrend + i] = slist[i+1];
            }
            free(slist);
        }
    }

    /* OPT_M: reject missing values within sample range */
    KPSSmod = lsq(list, dset, GretlCmdIndex.ols, Gretlopt.OPT_A | Gretlopt.OPT_M);
    if (KPSSmod.errcode) {
        clear_model(&KPSSmod);
        return KPSSmod.errcode;
    }

    t1 = KPSSmod.t1;
    t2 = KPSSmod.t2;
    T = KPSSmod.nobs;

    if (order < 0) {
        order = cast(int) (4.0 * pow(T / 100.0, 0.25));
    }

    /* record order for panel case */
    if (kinfo != null) {
        kinfo.order = order;
        if (order < kinfo.pmin) {
            kinfo.pmin = order;
        }
        if (order > kinfo.pmax) {
            kinfo.pmax = order;
        }
    }

    if (kinfo == null && (opt & Gretlopt.OPT_V)) {
        KPSSmod.aux = ModelAuxCode.kpss;
        printmodel(&KPSSmod, dset, Gretlopt.OPT_NONE, prn);
    }

    autocov = cast(double *) malloc(order * double.sizeof);
    if (autocov == null) {
        return GretlError.alloc;
    }

    for (i=0; i<order; i++) {
        autocov[i] = 0.0;
    }

    for (t=t1; t<=t2; t++) {
        et = KPSSmod.uhat[t];
        if (na(et)) {
            continue;
        }
        cumsum += et;
        cumsum2 += cumsum * cumsum;
        s2 += et * et;
        for (i=0; i<order; i++) {
            int s = i + 1;

            if (t - s >= t1) {
                autocov[i] += et * KPSSmod.uhat[t - s];
            }
        }
    }

    for (i=0; i<order; i++) {
        double wt = 1.0 - (cast(double) (i + 1)) / (order + 1);

        s2 += 2.0 * wt * autocov[i];
    }

    s2 /= T;

    if (s2 <= 0.0) {
        teststat = pval = NADBL;
    } else {
        teststat = cumsum2 / (s2 * T * T);
        pval = kpss_interp(teststat, T, hastrend);
    }

    if (kinfo != null) {
        /* storing info for panel test */
        kinfo.T = T;
        kinfo.test = teststat;
        kinfo.pval = pval;
    } else {
        /* testing individual time series */
        if (opt & Gretlopt.OPT_V) {
            pprintf(prn, "  %s: %g\n", _("Robust estimate of variance").ptr, s2);
            pprintf(prn, "  %s: %g\n", _("Sum of squares of cumulated residuals").ptr,
                    cumsum2);
        }
    }

    if (!(opt & Gretlopt.OPT_Q)) {
        double[3] a = [ 0.1, 0.05, 0.01 ];
        double[3] cv;

        cv[0] = kpss_critval(a[0], T, hastrend);
        cv[1] = kpss_critval(a[1], T, hastrend);
        cv[2] = kpss_critval(a[2], T, hastrend);

        pprintf(prn, _("\nKPSS test for %s").ptr, dset.varname[varno]);
        if (hastrend) {
            if (hasseas) {
                pputs(prn, _(" (including trend and seasonals)\n\n").ptr);
            } else {
                pputs(prn, _(" (including trend)\n\n").ptr);
            }
        } else {
            if (hasseas) {
                pputs(prn, _(" (including seasonals)\n\n").ptr);
            } else {
                pputs(prn, "\n\n");
            }
        }

        pprintf(prn, "T = %d\n", T);
        pprintf(prn, _("Lag truncation parameter = %d\n").ptr, order);
        pprintf(prn, "%s = %g\n\n", _("Test statistic").ptr, teststat);
        print_critical_values(a.ptr, cv.ptr, GretlCmdIndex.kpss, prn);
        if (pval == PV_GT10) {
            pprintf(prn, "%s > .10\n", _("P-value").ptr);
        } else if (pval == PV_LT01) {
            pprintf(prn, "%s < .01\n", _("P-value").ptr);
        } else if (!na(pval)) {
            pprintf(prn, "%s %.3f\n", _("Interpolated p-value").ptr, pval);
        }
        pputc(prn, '\n');
    }

    if (kinfo == null) {
        /* not doing the panel case */
        if (pval == PV_GT10 || pval == PV_LT01) {
            /* invalidate for record_test_result */
            pval = NADBL;
        }
        record_test_result(teststat, pval);
    }

    clear_model(&KPSSmod);
    free(list);
    free(autocov);

    return err;
}

private void panel_kpss_header(const(char) * vname, KpssInfo * kinfo,
                               Gretlopt opt, Prn * prn)
{
    pprintf(prn, _("\nKPSS test for %s %s\n").ptr, vname,
            ((opt & Gretlopt.OPT_T) ? _("(including trend)") : _("(without trend)")).ptr);
    if (kinfo.pmin > 0 && kinfo.pmax == kinfo.pmin) {
        pprintf(prn, _("Lag truncation parameter = %d\n").ptr, kinfo.pmin);
    } else if (kinfo.pmin > 0 && kinfo.pmax > kinfo.pmin) {
        pprintf(prn, _("Lag truncation parameter: %d - %d\n").ptr,
                kinfo.pmin, kinfo.pmax);
    } else {
        pputs(prn, _("Lag truncation automatic\n").ptr);
    }
    pputc(prn, '\n');
}

private int panel_kpss_test(int order, int v, DatasetRaw * dset,
                            Gretlopt opt, Prn * prn)
{
    KpssInfo kinfo;
    int u0 = dset.t1 / dset.pd;
    int uN = dset.t2 / dset.pd;
    int n = uN - u0 + 1;
    int verbose = (opt & Gretlopt.OPT_V);
    double ppv = 0.0, zpv = 0.0, lpv = 0.0;
    int gt_10 = 0, lt_01 = 0;
    double pval;
    int i, err = 0;

    if (order < 0) {
        kinfo.pmin = 1000000;
        kinfo.pmax = -1;
    } else {
        kinfo.pmin = kinfo.pmax = order;
    }

    if (verbose) {
        panel_kpss_header(dset.varname[v], &kinfo, opt, prn);
    }

    /* run a KPSS test for each unit and record the
       results */

    for (i=u0; i<=uN && !err; i++) {
        dset.t1 = i * dset.pd;
        dset.t2 = dset.t1 + dset.pd - 1;
        err = series_adjust_sample(dset.Z[v], &dset.t1, &dset.t2);
        if (!err) {
            err = real_kpss_test(order, v, dset, opt | Gretlopt.OPT_Q, &kinfo, prn);
            if (!err && verbose) {
                if (order < 0) {
                    pprintf(prn, "Unit %d, T = %d, order %d\n", i + 1,
                            kinfo.T, kinfo.order);
                } else {
                    pprintf(prn, "Unit %d, T = %d\n", i + 1, kinfo.T);
                }
                if (na(kinfo.pval)) {
                    pputs(prn, "\n\n");
                } else {
                    pprintf(prn, "test = %g, ", kinfo.test);
                    if (kinfo.pval == PV_GT10) {
                        pprintf(prn, "%s > .10\n", _("p-value").ptr);
                    } else if (kinfo.pval == PV_LT01) {
                        pprintf(prn, "%s < .01\n", _("p-value").ptr);
                    } else {
                        pprintf(prn, "%s %.3f\n", _("interpolated p-value").ptr,
                                kinfo.pval);
                    }
                    pputc(prn, '\n');
                }
            }
        }

        if (!err) {
            pval = kinfo.pval;

            if (pval == PV_GT10) {
                gt_10++;
                if (lt_01 == 0) {
                    /* record lower bound */
                    pval = .10;
                } else {
                    pval = NADBL;
                }
            } else if (pval == PV_LT01) {
                lt_01++;
                if (gt_10 == 0) {
                    /* record upper bound */
                    pval = .01;
                } else {
                    pval = NADBL;
                }
            }

            if (na(pval)) {
                ppv = zpv = lpv = NADBL;
            } else if (!na(ppv)) {
                ppv += log(pval);
                zpv += normal_cdf_inverse(pval);
                lpv += log(pval / (1-pval));
            }
        }
    }

    if (!err && !verbose) {
        panel_kpss_header(dset.varname[v], &kinfo, opt, prn);
    }

    if (!err && !na(ppv)) {
        /* process the results as per Choi, as best we can */
        pprintf(prn, "%s\n\n", _("H0: all groups are stationary").ptr);
        do_choi_test(ppv, zpv, lpv, n, prn);
        if (gt_10 > 0) {
            pprintf(prn, "   %s\n", _("Note: these are LOWER BOUNDS on the true p-values").ptr);
            pprintf(prn, "   (%s: %d)\n", _("Individual p-values > .10, and recorded as .10").ptr,
                    gt_10);
        } else if (lt_01 > 0) {
            pprintf(prn, "   %s\n", _("Note: these are UPPER BOUNDS on the true p-values").ptr);
            pprintf(prn, "   (%s: %d)\n", _("Individual p-values < .01, and recorded as .01").ptr,
                    lt_01);
        }
        pputc(prn, '\n');
    } else {
        pprintf(prn, _("Choi test: cannot be calculated\n").ptr);
    }

    return err;
}

/*
 * kpss_test:
 * @order: window size for Bartlett smoothing.
 * @list: list of variables to test.
 * @dset: dataset struct.
 * @opt: option flags.
 * @prn: gretl printing struct.
 *
 * Carries out and prints the results of the KPSS test for
 * stationarity. Flags that may be given in @opt include:
 * OPT_T to include a linear trend; OPT_F to apply
 * first-differencing before testing; OPT_V for verbose
 * operation.
 *
 * Returns: 0 on successful completion, non-zero on error.
 */

int kpss_test(int order, const(int) * list, DatasetRaw * dset,
              Gretlopt opt, Prn * prn)
{
    int save_t1 = dset.t1;
    int save_t2 = dset.t2;
    int orig_nvars = dset.v;
    int err = 0;

    if (multi_unit_panel_sample(dset)) {
        err = panel_kpss_test(order, list[1], dset, opt, prn);
    } else {
        /* regular time series case */
        int i, v;
        int[2] vlist = [1, 0];

        for (i=1; i<=list[0] && !err; i++) {
            v = list[i];
            vlist[1] = v;
            err = list_adjust_sample(vlist.ptr, &dset.t1, &dset.t2,
                                     dset, null);
            if (!err) {
                err = real_kpss_test(order, v, dset, opt, null, prn);
            }
            dset.t1 = save_t1;
            dset.t2 = save_t2;
        }
    }

    dset.t1 = save_t1;
    dset.t2 = save_t2;

    /* added 2012-03-22 for consistency with adf test */
    dataset_drop_last_variables(dset, dset.v - orig_nvars);

    return err;
}

private int * make_coint_list(const(int) * list, DetCode detcode,
                              Gretlopt opt, int * nv,
                              DatasetRaw * dset, int * err)
{
    int * clist = null;
    int nseas = 0;
    int ifc = 0;
    int i, j = 1;

    /* does the incoming list contain a constant? */
    for (i=1; i<=list[0]; i++) {
        if (list[i] == 0) {
            ifc = 1;
            break;
        }
    }

    /* check for sufficient arguments */
    *nv = list[0] - ifc;
    if (*nv < 2) {
        *err = GretlError.args;
        return null;
    }

    /* space for seasonals, if applicable */
    if (!*err && (opt & Gretlopt.OPT_D) && dset.pd > 1) {
        nseas = dset.pd - 1;
    }

    /* allocate list for cointegrating regression */
    clist = gretl_list_new(*nv + detcode - 1 + nseas);
    if (clist == null) {
        *err = GretlError.alloc;
        return null;
    }

    /* transcribe original vars */
    for (i=1; i<=list[0]; i++) {
        if (list[i] != 0) {
            clist[j++] = list[i];
        }
    }

    /* add trend, if wanted */
    if (detcode >= DetCode.trend) {
        clist[j] = gettrend(dset, 0);
        if (clist[j++] == 0) {
            *err = GretlError.alloc;
        }
    }

    /* add trend-squared, if wanted */
    if (!*err && detcode == DetCode.quad_trend) {
        clist[j] = gettrend(dset, 1);
        if (clist[j++] == 0) {
            *err = GretlError.alloc;
        }
    }

    /* add const, if wanted */
    if (!*err && detcode != DetCode.no_const) {
        clist[j++] = 0;
    }

    /* add centered seasonals, if wanted */
    if (nseas > 0) {
        int * slist = seasonals_list(dset, dset.pd, 1, err);

        if (!*err) {
            for (i=0; i<nseas; i++) {
                clist[j++] = slist[i+1];
            }
            free(slist);
        }
    }

    return clist;
}

private int
coint_check_opts(Gretlopt opt, DetCode * detcode, Gretlopt * adf_opt)
{
    if (opt & Gretlopt.OPT_N) {
        if ((opt & Gretlopt.OPT_T) || (opt & Gretlopt.OPT_R)) {
            return GretlError.badopt;
        } else {
            *detcode = DetCode.no_const;
            *adf_opt = Gretlopt.OPT_N;
        }
    } else if (opt & Gretlopt.OPT_T) {
        if (opt & Gretlopt.OPT_R) {
            return GretlError.badopt;
        } else {
            *detcode = DetCode.trend;
            *adf_opt = Gretlopt.OPT_T;
        }
    } else if (opt & Gretlopt.OPT_R) {
        *detcode = DetCode.quad_trend;
        *adf_opt = Gretlopt.OPT_R;
    }

    if (opt & Gretlopt.OPT_D) {
        *adf_opt |= Gretlopt.OPT_D;
    }

    if (opt & Gretlopt.OPT_E) {
        *adf_opt |= Gretlopt.OPT_E;
    }

    return 0;
}

/* Engle-Granger: try to ensure a uniform sample for the individual
   (A)DF tests and the test on the cointegrating regression
*/

private int coint_set_sample(const(int) * list, int nv, int order,
                             DatasetRaw * dset)
{
    int anymiss;
    int i, v, t;

    for (t=dset.t1; t<dset.t2; t++) {
        anymiss = 0;
        for (i=1; i<=nv; i++) {
            v = list[i];
            if (na(dset.Z[v][t])) {
                anymiss = 1;
                break;
            }
        }
        if (!anymiss) {
            break;
        }
    }

    dset.t1 = t + order + 1;

    for (t=dset.t2; t>dset.t1; t--) {
        anymiss = 0;
        for (i=1; i<=nv; i++) {
            v = list[i];
            if (na(dset.Z[v][t])) {
                anymiss = 1;
                break;
            }
        }
        if (!anymiss) {
            break;
        }
    }

    dset.t2 = t;

    return 0;
}

/*
 * engle_granger_test:
 * @order: lag order for the test.
 * @list: specifies the variables to use.
 * @dset: dataset struct.
 * @opt: option flags.
 * @prn: gretl printing struct.
 *
 * Carries out the Engle-Granger test for cointegration.
 * Flags that may be given in @opt include: OPT_N, do
 * not an include a constant in the cointegrating regression;
 * OPT_T include constant and linear trend; OPT_R, include
 * quadratic trend; OPT_S, skip DF tests for individual variables;
 * OPT_E, test down from maximum lag order (see the entry for
 * "adf" in the Gretl Command Reference for details); OPT_V,
 * verbose operation.
 *
 * Returns: 0 on successful completion, non-zero code
 * on error.
 */

int engle_granger_test(int order, const(int) * list, DatasetRaw * dset,
                       Gretlopt opt, Prn * prn)
{
    int orig_t1 = dset.t1;
    int orig_t2 = dset.t2;
    AdfInfo ainfo;
    Gretlopt adf_opt = Gretlopt.OPT_C;
    Model cmod;
    DetCode detcode = DetCode.const_;
    int i, nv, k = 0;
    int step = 1;
    int skip = 0;
    int silent = 0;
    int * clist = null;
    int err = 0;

    if (multi_unit_panel_sample(dset)) {
        gretl_errmsg_set(_("Sorry, this command is not yet available "
                           ~ "for panel data").ptr);
        return GretlError.data;
    }

    err = coint_check_opts(opt, &detcode, &adf_opt);
    if (err) {
        return err;
    }

    clist = make_coint_list(list, detcode, opt, &nv, dset, &err);
    if (err) {
        return err;
    }

    /* backward compatibility: let a negative lag order
       indicate that we should test down */
    if (order < 0) {
        order = -order;
        adf_opt |= Gretlopt.OPT_E;
    }

    /* verbosity? */
    if (opt & Gretlopt.OPT_V) {
        adf_opt |= Gretlopt.OPT_V;
    }

    /* or silence? */
    if (opt & Gretlopt.OPT_I) {
        adf_opt |= Gretlopt.OPT_I;
        silent = skip = 1;
    } else if (opt & Gretlopt.OPT_S) {
        skip = 1;
    }

    gretl_model_init(&cmod, dset);

    if (!skip) {
        /* start by testing all candidate vars for unit root */
        int uniform_order = (opt & Gretlopt.OPT_E) ? 0 : order;

        coint_set_sample(clist, nv, uniform_order, dset);
        for (i=1; i<=nv; i++) {
            ainfo.v = clist[i];
            ainfo.order = order;
            ainfo.niv = 1;
            ainfo.flags = AdfFlags.eg_test;

            if (step == 1) {
                pputc(prn, '\n');
            }
            pprintf(prn, _("Step %d: testing for a unit root in %s\n").ptr,
                    step++, dset.varname[ainfo.v]);
            real_adf_test(&ainfo, dset, adf_opt, prn);
        }
    }

    if (!silent) {
        if (step == 1) {
            pputc(prn, '\n');
        }
        pprintf(prn, _("Step %d: cointegrating regression\n").ptr, step++);
    }

    dset.t1 = orig_t1;
    dset.t2 = orig_t2;

    cmod = lsq(clist, dset, GretlCmdIndex.ols, Gretlopt.OPT_NONE);
    err = cmod.errcode;
    if (err) {
        goto bailout;
    }

    if (!silent) {
        cmod.aux = ModelAuxCode.coint;
        printmodel(&cmod, dset, Gretlopt.OPT_NONE, prn);
    }

    /* add residuals from cointegrating regression to data set */
    err = dataset_add_allocated_series(dset, cmod.uhat);
    if (err) {
        goto bailout;
    }

    k = dset.v - 1;
    strcpy(dset.varname[k], "uhat");
    cmod.uhat = null;

    if (!silent) {
        pprintf(prn, _("Step %d: testing for a unit root in %s\n").ptr,
                step, dset.varname[k]);
    }

    ainfo.v = k;
    ainfo.order = order;
    ainfo.niv = nv;
    ainfo.flags = AdfFlags.eg_test | AdfFlags.eg_resids;

    /* Run (A)DF test on the residuals */
    real_adf_test(&ainfo, dset, adf_opt, prn);

    if (!silent) {
        pputs(prn, _("\nThere is evidence for a cointegrating relationship if:\n"
                     ~ "(a) The unit-root hypothesis is not rejected for the individual"
                     ~ " variables, and\n(b) the unit-root hypothesis is rejected for the "
                     ~ "residuals (uhat) from the \n    cointegrating regression.\n").ptr);
        pputc(prn, '\n');
    }

 bailout:

    clear_model(&cmod);
    free(clist);
    if (k > 0) {
        dataset_drop_variable(k, dset);
    }

    dset.t1 = orig_t1;
    dset.t2 = orig_t2;

    return err;
}
