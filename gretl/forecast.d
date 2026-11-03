/* Port of lib/src/forecast.c: forecasts based on estimated models (static, dynamic and k-step, with standard
   errors) for OLS-type, AR, ARMA, GARCH, NLS, binary/ordered/count models, VARs and equation systems, the
   fitted-and-residual values of a model and the user-level "fcast" command (the header forecast.h is forecasth.d).

   Differences from the C:
   - the AR_DEBUG, NLS_DEBUG, MIDAS_DEBUG, GFC_DEBUG and KSTEP_DEBUG code (compiled out in the C), the DPRINTF
     macro and the `#if 0` alternative in nls_fcast() are not ported;
   - the file-internal struct Forecast_ is Forecast, the anonymous out-of-sample enum is OsCode;
   - the function-like macros are functions.
   Dependencies not yet ported are declared extern (C) below, in the style of the other ports. */

module gretl.forecast;


import core.stdc.errno : errno;
import core.stdc.math : exp, fabs, floor, isnan, log, pow, sqrt, ceil;
import core.stdc.stdio : fprintf, sprintf, stderr;
import core.stdc.stdlib : calloc, free, malloc, qsort, realloc;
import core.stdc.string : memcpy, memset, strcmp, strcpy, strlen;

import std.algorithm : max, min;

import gretl.forecasth;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_matrix;
import gretl.gretl_matrixh;
import gretl.libgretlh;
import gretl.missingh : NADBL, na;
import gretl.systemh;
import gretl.varh;
import core.stdc.stdio : sscanf;
import core.stdc.string : strncat;
import gretl.dataio : check_varname, dateton, get_precision, ntolabel;
import gretl.dataset : dataset_add_NA_series, dataset_add_allocated_series, dataset_add_series, dataset_drop_last_variables, is_periodic_dummy, is_trend_variable, rename_series, series_get_lag, series_get_parent_id, series_get_transform, series_is_log, series_set_label;
import gretl.dataseth : dataset_is_panel, dataset_is_time_series;
import gretl.describe : doornik_chisq, series_get_skew_kurt;
import gretl.discrete : mn_logit_prediction, mn_logit_probabilities, ordered_model_prediction, ordered_probabilities;
import gretl.estimate : lsq;
import gretl.genmain : current_series_index, generate, generate_scalar, gretl_reserved_word, series_index;
import gretl.genmainh : ModelDataIndex;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_errors : gretl_errmsg_set, gretl_errmsg_sprintf, gretl_error_clear;
import gretl.gretl_list : gretl_consecutive_list_new, gretl_list_new, gretl_null_list;
import gretl.gretl_model : arima_delta_coeffs, arma_model_AR_MA_coeffs, arma_model_get_x_coeffs, arma_model_max_AR_lag, arma_model_max_MA_lag, clear_model, gretl_is_arima_model, gretl_model_ahat_vec, gretl_model_get_data, gretl_model_get_depvar, gretl_model_get_double, gretl_model_get_int, gretl_model_get_x_list, regarma_model_AR_coeffs;
import gretl.gretl_modelh : AR_MODEL, ASYMPTOTIC_MODEL, ArmaFlags, COUNT_MODEL, LIMDEP, SIMPLE_AR_MODEL;
import gretl.gretl_prn : Prn, pprintf, pputc;
import gretl.gretl_utils : copyvec, gretl_int_from_string, gretl_isdummy;
import gretl.libset : gretl_echo_on, gretl_messages_on;
import gretl.matrix_extra : gretl_coeff_vector_from_model, gretl_matrix_print_to_prn, gretl_vcv_matrix_from_model;
import gretl.objstackh : GretlObjType;
import gretl.options : get_optval_int, inapplicable_option_error, incompatible_options;
import gretl.pvalues : chisq_cdf_comp, normal_cdf;
import gretl.strutils : gretl_strdup, strings_array_new;
import gretl.transforms : is_standard_diff;
import gretl.printout : maybe_trim_varname, print_fcast_stats_matrix, text_print_forecast;
import gretl.var : gretl_VAR_get_forecast_matrix;
import gretl.subsample : fetch_full_dataset, get_dataset_submask_size, same_dataset;
import gretl.nls : nl_model_run_aux_genrs;
import gretl.system : system_get_forecast_matrix;
import gretl.uservar : destroy_private_uvars, gretl_is_scalar, gretl_scalar_get_value, user_var_add_or_replace;
import gretl.genfuncs : forecast_stats;

/* unported dependencies: per-file declarations from the 2026c headers */
extern (C) {
    void * get_last_model(GretlObjType * type);
    int midas_forecast_setup(const(Model) * pmod, DatasetRaw * dset, FcastMethod method, char ** pformula);
}
/**
 * SECTION:forecast
 * @short_description: Forecasting based on models
 * @title: Forecasting
 * @include: gretl/libgretl.h, gretl/forecast.h
 *
 * Write-up needed here.
 *
 */

/* estimators where a simple X*b does _not_ give the
   predicted value of the dependent variable */

bool FCAST_SPECIAL(int c)
{
    return (c == GretlCmdIndex.logit ||
            c == GretlCmdIndex.logistic ||
            c == GretlCmdIndex.midasreg ||
            c == GretlCmdIndex.negbin ||
            c == GretlCmdIndex.nls ||
            c == GretlCmdIndex.poisson ||
            c == GretlCmdIndex.probit ||
            c == GretlCmdIndex.tobit);
}

struct Forecast {
    FcastMethod method; /* static, dynamic or auto */
    double * yhat; /* array of forecast values */
    double * sderr; /* array of forecast standard errors */
    double * eps; /* array of estimated forecast errors */
    int * dvlags; /* info on lags of dependent variable */
    int t1; /* start of forecast range */
    int t2; /* end of forecast range */
    int model_t2; /* end of period over which model was estimated */
}

__gshared GretlMatrixRaw * fcast_matrix;
__gshared GretlMatrixRaw * fcerr_matrix;

void set_fcast_matrices(GretlMatrixRaw * fc, GretlMatrixRaw * se)
{
    gretl_matrix_free(fcast_matrix);
    gretl_matrix_free(fcerr_matrix);

    fcast_matrix = fc;
    fcerr_matrix = se;
}

/* Create an empty, dummy AR info structure for use with models
   that don't have an explicit AR error process, but that do
   have a lagged dependent variable that in effect produces
   an AR error, for forecasting purposes. Or add a single-
   element AR info structure for the AR1 case.
*/

int dummy_ar_info_init(Model * pmod)
{
    int err = 0;

    pmod.arinfo = cast(typeof(pmod.arinfo)) malloc((*pmod.arinfo).sizeof);

    if (pmod.arinfo == null) {
    err = GretlError.alloc;
    } else {
    pmod.arinfo.rho = null;
    pmod.arinfo.sderr = null;
    }

    if (!err && pmod.ci == GretlCmdIndex.ar1) {
    pmod.arinfo.arlist = gretl_list_new(1);
    pmod.arinfo.rho = cast(typeof(pmod.arinfo.rho)) malloc(double.sizeof);
    pmod.arinfo.arlist[1] = 1;
    pmod.arinfo.rho[0] = gretl_model_get_double(pmod, "rho_gls");
    } else if (!err) {
    pmod.arinfo.arlist = gretl_null_list();
    }

    if (err) {
    free(pmod.arinfo);
    pmod.arinfo = null;
    }

    return err;
}

int fit_resid_add_sderr(FitResid * fr)
{
    int t, err = 0;

    fr.sderr = cast(typeof(fr.sderr)) malloc(fr.nobs * (*fr.sderr).sizeof);

    if (fr.sderr == null) {
    err = GretlError.alloc;
    } else {
    for (t=0; t<fr.nobs; t++) {
        fr.sderr[t] = NADBL;
    }
    }

    return err;
}

/**
 * free_fit_resid:
 * @fr: the pointer to be freed.
 *
 * Frees all resources associated with @fr, then frees the pointer
 * itself.
 */

void free_fit_resid(FitResid * fr)
{
    if (fr != null) {
    free(fr.actual);
    free(fr.fitted);
    free(fr.resid);
    free(fr.sderr);
    free(fr);
    }
}

/**
 * fit_resid_new_with_length:
 * @n: the number of observations to allow for, or 0 if this
 * information will be added later.
 *
 * Allocates a #FITRESID struct for holding fitted values and
 * residuals from a model (or out-of-sample forecasts).  If
 * @n is greater than 0 the arrays required for that number
 * of observations will be allocated.
 *
 * Returns: pointer to allocated structure, or %NULL on failure.
 */

FitResid * fit_resid_new_with_length(int n, int add_errs)
{
    FitResid * f = cast(FitResid *) malloc((FitResid).sizeof);

    if (f == null) {
    return null;
    }

    f.method = cast(FcastMethod) 0;
    f.model_ID = 0;
    f.asymp = 0;
    f.std = 0;
    f.model_t1 = 0;
    f.t0 = 0;
    f.t1 = 0;
    f.t2 = 0;
    f.df = 0;
    f.k = 0;
    f.nobs = 0;
    f.pmax = PMAX_NOT_AVAILABLE;

    f.sigma = NADBL;
    f.alpha = 0.05;
    f.a0 = NADBL;
    f.a0meth = 0;

    f.actual = null;
    f.fitted = null;
    f.resid = null;
    f.sderr = null;

    f.depvar[0] = '\0';
    f.ylabel[0] = '\0';

    f.actual = cast(typeof(f.actual)) malloc(n * (*f.actual).sizeof);
    f.fitted = cast(typeof(f.fitted)) malloc(n * (*f.fitted).sizeof);
    f.resid = cast(typeof(f.resid)) malloc(n * (*f.resid).sizeof);
    if (add_errs) {
    f.sderr = cast(typeof(f.sderr)) malloc(n * (*f.sderr).sizeof);
    }

    if (f.actual == null || f.fitted == null || f.resid == null ||
    (add_errs && f.sderr == null)) {
    free_fit_resid(f);
    f = null;
    } else {
    int t;

    for (t=0; t<n; t++) {
        f.actual[t] = f.fitted[t] = f.resid[t] = NADBL;
        if (f.sderr != null) {
        f.sderr[t] = NADBL;
        }
    }
    f.nobs = n;
    }

    return f;
}

FitResid * fit_resid_new_for_model(const(Model) * pmod, const(DatasetRaw) * dset, int t1, int t2, int pre_n, Gretlopt opt, int * err)
{
    FitResid * fr;

    if (t1 < 0 || t2 < 0 || t2 < t1) {
    *err = GretlError.obs;
    return null;
    }

    if (dset != null) {
        fr = fit_resid_new_with_length(dset.n, 0);
    } else {
        fr = fit_resid_new_with_length(t2 - t1 + 1, 0);
    }

    if (fr == null) {
    *err = GretlError.alloc;
    return null;
    }

    fr.t1 = t1;
    fr.t2 = t2;

    if (pre_n > 0) {
    fr.t0 = fr.t1 - pre_n;
    } else {
    fr.t0 = t1;
    }

    fr.model_ID = pmod.ID;
    fr.model_t1 = pmod.t1;

    fr.asymp = ASYMPTOTIC_MODEL(cast(GretlCmdIndex) pmod.ci);
    fr.opt = opt;

    return fr;
}

FitResid * fit_resid_new_for_system(int asy, const(DatasetRaw) * dset, int t1, int t2, int pre_n, int * err)
{
    FitResid * fr;

    if (t1 < 0 || t2 < 0 || t2 < t1) {
    *err = GretlError.obs;
    return null;
    }

    fr = fit_resid_new_with_length(dset.n, 1);

    if (fr == null) {
    *err = GretlError.alloc;
    return null;
    }

    fr.t1 = t1;
    fr.t2 = t2;

    if (pre_n > 0) {
    fr.t0 = fr.t1 - pre_n;
    } else {
    fr.t0 = t1;
    }

    fr.asymp = asy;

    return fr;
}

void fit_resid_set_dec_places(FitResid * fr)
{
    int n = fr.t2 - fr.t0 + 1;

    if (gretl_isdummy(fr.t0, fr.t2, fr.actual) > 0) {
    fr.pmax = get_precision(fr.fitted + fr.t0, n, 6);
    } else {
    fr.pmax = get_precision(fr.actual + fr.t0, n, 6);
    }

    if (fr.pmax == 0) {
    fr.pmax = 2;
    }
}

int special_depvar(const(Model) * pmod)
{
    int ret = 0;

    if (pmod.ci == GretlCmdIndex.intreg) {
    ret = 1;
    } else if (pmod.ci == GretlCmdIndex.panel) {
    if (pmod.opt & Gretlopt.OPT_B) {
        /* panel "between" model */
        ret = 1;
    }
    }

    return ret;
}

/**
 * get_fit_resid:
 * @pmod: the model for which actual and fitted values
 * are wanted.
 * @dset: dataset struct.
 * @err: location to receive error code.
 *
 * Allocates a #FITRESID structure and fills it out with
 * the actual and predicted values of the dependent variable
 * in @pmod.
 *
 * Returns: pointer to allocated structure, or %NULL on failure.
 */

FitResid * get_fit_resid(const(Model) * pmod, const(DatasetRaw) * dset, int * err)
{
    double * uhat = cast(double *) pmod.uhat;
    double * yhat = cast(double *) pmod.yhat;
    FitResid * fr;
    int t, dv = -1;

    if (!special_depvar(pmod)) {
    dv = gretl_model_get_depvar(pmod);
    if (dv < 0 || dv >= dset.v) {
        *err = GretlError.data;
        return null;
    }
    }

    if (uhat == null || yhat == null) {
    GretlMatrixRaw * uh = cast(GretlMatrixRaw *) gretl_model_get_data(pmod, "uhat");
    GretlMatrixRaw * yh = cast(GretlMatrixRaw *) gretl_model_get_data(pmod, "yhat");

    if (uh == null || yh == null) {
        *err = GretlError.data;
        return null;
    } else {
        uhat = uh.val;
        yhat = yh.val;
    }
    }

    fr = fit_resid_new_for_model(pmod, dset, pmod.t1, pmod.t2,
                 0, Gretlopt.OPT_NONE, err);
    if (*err) {
    return null;
    }

    if (LIMDEP(cast(GretlCmdIndex) pmod.ci)) {
    fr.sigma = NADBL;
    } else if (pmod.ci == GretlCmdIndex.garch && (pmod.opt & Gretlopt.OPT_Z)) {
    fr.sigma = 1.0;
    fr.std = 1;
    } else {
    fr.sigma = gretl_model_get_double(pmod, "sigma_orig");
    if (na(fr.sigma)) {
        fr.sigma = pmod.sigma;
    }
    }

    for (t=0; t<fr.nobs; t++) {
    if (dv < 0) {
        if (na(yhat[t]) || na(uhat[t])) {
        fr.actual[t] = NADBL;
        } else {
        fr.actual[t] = yhat[t] + uhat[t];
        }
    } else {
        fr.actual[t] = dset.Z[dv][t];
    }
    fr.fitted[t] = yhat[t];
    fr.resid[t] = uhat[t];
    }

    if (dv < 0) {
    fr.pmax = PMAX_NOT_AVAILABLE;
    } else {
    fit_resid_set_dec_places(fr);
    }

    if (dv < 0) {
    const(char) * dvname = pmod.depvar;

    if (dvname != null) {
        strcpy(fr.depvar.ptr, dvname);
    } else {
        strcpy(fr.depvar.ptr, "implicit y");
    }
    } else {
    strcpy(fr.depvar.ptr, dset.varname[dv]);
    }

    return fr;
}

int xy_equal(double x, double y)
{
    return x == y || (na(x) && na(y));
}

int detect_lag(int vy, int vx, int t1, int t2, int pmax, const(DatasetRaw) * dset)
{
    const(double) * y = dset.Z[vy];
    const(double) * x = dset.Z[vx];
    int p, t, ok;
    int ret = 0;

    for (p=1; p<=pmax; p++) {
        ok = 1;
        for (t=t1+p; t<=t2; t++) {
            if (!xy_equal(x[t], y[t-p])) {
                ok = 0;
                break;
            }
        }
        if (ok) {
            ret = p;
            break;
        }
    }

    return ret;
}

/* Makes a list to keep track of any "independent variables" that are
   actually lags of the dependent variable.  This list has as many
   elements as the model has independent variables, and in each place we
   write a zero if the regressor is not a lag of the dependent variable
   or a positive integer recording the lag order.  However, if the list
   of regressors contains no lags of the dependent var, the function
   returns NULL.
*/

int * process_lagged_depvar(Model * pmod, const(DatasetRaw) * dset)
{
    int * dvlags = null;
    int nl = 0;

    if (pmod.xlist != null) {
        dvlags = cast(typeof(dvlags)) malloc(pmod.xlist[0] * (*dvlags).sizeof);
    }

    if (dvlags != null) {
        int yv = gretl_model_get_depvar(pmod);
    int i, vi, lag;

    for (i=1; i<=pmod.xlist[0]; i++) {
            lag = 0;
            vi = pmod.xlist[i];
            if (series_get_transform(dset, vi) == GretlCmdIndex.lags) {
                if (series_get_parent_id(dset, vi) == yv) {
                    lag = series_get_lag(dset, vi);
                }
            } else if (vi > 0) {
                int pmax = dset.pd < 4 ? 4 : dset.pd;

                lag = detect_lag(yv, vi, dset.t1, dset.t2,
                                 pmax, dset);
            }
            dvlags[i-1] = lag;
            if (lag > 0) {
                nl++;
            }
        }
    }

    if (nl == 0) {
        /* didn't find any */
        free(dvlags);
        dvlags = null;
    }

    return dvlags;
}

/* Get a value for a lag of the dependent variable.  If method is
   dynamic we prefer lagged prediction to lagged actual.  If method is
   static, we never want the lagged prediction, only the actual.  If
   method is "auto", which value we prefer depends on whether we're in
   or out of sample (actual within, lagged prediction without).
*/

double fcast_get_ldv(Forecast * fc, int i, int t, int lag, const(DatasetRaw) * dset)
{
    int p = t - lag;
    double ldv;

    /* initialize to actual lagged value, if available */
    if (p < 0) {
    ldv = NADBL;
    } else {
    ldv = dset.Z[i][p];
    }

    if (fc.method != FcastMethod.static_ && p >= 0) {
    double yhat = fc.yhat[p];
    int fc_ok = !na(yhat);

    if (fc.method == FcastMethod.dynamic) {
        if (fc_ok) {
        ldv = yhat;
        }
    } else if (fc.method == FcastMethod.auto_) {
        if (t > fc.model_t2 + lag || na(ldv)) {
        if (fc_ok) {
            ldv = yhat;
        }
        }
    }
    }

    return ldv;
}

/* Get forecasts, plus standard errors for same, for models without
   autoregressive errors and without "special requirements"
   (e.g. nonlinearity).

   The forecast standard errors include both uncertainty over the
   error process and parameter uncertainty (Davidson and MacKinnon
   method).
*/

int static_fcast_with_errs(Forecast * fc, Model * pmod, const(DatasetRaw) * dset, Gretlopt opt)
{
    GretlMatrixRaw * V = null;
    GretlVector * Xt = null;
    GretlVector * b = null;
    double s2 = pmod.sigma * pmod.sigma;
    double vyh, xval;
    int k = pmod.ncoeff;
    int same_data = 0;
    int i, vi, t;
    int err = 0;

    V = gretl_vcv_matrix_from_model(pmod, null, &err);
    if (err) {
    goto bailout;
    }

    b = gretl_coeff_vector_from_model(pmod, null, &err);
    if (err) {
    goto bailout;
    }

    Xt = gretl_vector_alloc(k);
    if (Xt == null) {
    err = GretlError.alloc;
    goto bailout;
    }

    same_data = same_dataset(pmod, dset);

    for (t=fc.t1; t<=fc.t2 && !err; t++) {
    int missing = 0;

    /* skip if we can't compute forecast */
    if (same_data && t >= pmod.t1 && t <= pmod.t2) {
        missing = na(pmod.yhat[t]);
    }

    /* populate Xt vector for observation */
    for (i=0; i<k && !missing; i++) {
        vi = pmod.list[i + 2];
        xval = dset.Z[vi][t];
        if (na(xval)) {
        fc.sderr[t] = fc.yhat[t] = NADBL;
        missing = 1;
        } else {
        gretl_vector_set(Xt, i, xval);
        }
    }

    if (missing) {
        fc.sderr[t] = fc.yhat[t] = NADBL;
        continue;
    }

    /* forecast value */
    fc.yhat[t] = gretl_vector_dot_product(Xt, b, null);

    /* forecast variance */
    vyh = gretl_scalar_qform(Xt, V, &err);
    if (na(vyh)) {
        err = 1;
    } else {
        vyh += s2;
        if (vyh >= 0.0) {
        fc.sderr[t] = sqrt(vyh);
        } else {
        fc.sderr[t] = NADBL;
        err = 1;
        }
    }
    }

 bailout:

    gretl_matrix_free(V);
    gretl_vector_free(Xt);
    gretl_vector_free(b);

    return err;
}

bool individual_effects_model(const(Model) * m)
{
    return (m != null && m.ci == GretlCmdIndex.panel &&
            (m.opt & (Gretlopt.OPT_F | Gretlopt.OPT_U)));
}

void special_set_fcast_matrix(GretlMatrixRaw * fc)
{
    size_t sz;

    /* copy forecast from second column to first */
    sz = fc.rows * double.sizeof;
    memcpy(fc.val, fc.val + fc.rows, sz);

    /* and release some storage */
    gretl_matrix_realloc(fc, fc.rows, 1);

    set_fcast_matrices(fc, null);
}

void special_print_fc_stats(GretlMatrixRaw * fc, Prn * prn)
{
    const(double) * y = fc.val;
    const(double) * f = y + fc.rows;
    GretlMatrixRaw * fcs;
    int n_used;
    int err = 0;

    fcs = forecast_stats(y, f, 0, fc.rows - 1, &n_used,
             Gretlopt.OPT_O, &err);
    if (fcs != null) {
    print_fcast_stats_matrix(fcs, n_used, Gretlopt.OPT_NONE, prn);
    gretl_matrix_free(fcs);
    }
}

void transcribe_to_matrix(GretlMatrixRaw * fc, double yt, double fct, char ** S, const(DatasetRaw) * dset, int t, int j)
{
    if (S != null) {
    char[OBSLEN] obsstr;

    ntolabel(obsstr.ptr, t, dset);
    S[j] = gretl_strdup(obsstr.ptr);
    }
    if (fc.cols == 2) {
    gretl_matrix_set(fc, j, 0, yt);
    gretl_matrix_set(fc, j, 1, fct);
    } else {
    gretl_matrix_set(fc, j, 0, fct);
    }
}

/* Specific to "dpanel": handles both within-sample and out of sample
   forecasts.
*/

int real_dpanel_fcast(double * yhat, GretlMatrixRaw * fc, char ** rlabels, Model * pmod, const(char) * mask, const(DatasetRaw) * dset, Gretlopt opt)
{
    const(double) * y, xi, b;
    double * yi = null;
    int * xlist, ylags = null;
    int yno, free_ylags = 0;
    int fc_static = (opt & Gretlopt.OPT_S);
    int dpdstyle = (pmod.opt & Gretlopt.OPT_X);
    int t1min, ntdum, ifc, tdt = 0;
    int i, vi, s, t, k, u;
    int T, ti, row;

    ylags = cast(typeof(ylags)) gretl_model_get_data(pmod, "ylags");
    if (ylags == null) {
    ylags = gretl_consecutive_list_new(1, pmod.list[1]);
    free_ylags = 1;
    }

    /* per unit time series */
    T = dset.pd;
    yi = cast(typeof(yi)) malloc(T * (*yi).sizeof);

    if (ylags == null || yi == null) {
    return GretlError.alloc;
    }

    if (mask == null && !(opt & Gretlopt.OPT_D)) {
    /* within-sample: default to static */
    fc_static = 1;
    }

    /* first usable time-series observation, zero based */
    t1min = gretl_model_get_int(pmod, "t1min");
    if (t1min > 0) {
    t1min--;
    } else {
    /* fallback? problematic if we have time dummies */
    t1min = 1 + ylags[ylags[0]];
    }

    /* is a constant included? */
    ifc = gretl_model_get_int(pmod, "ifc");

    /* do we have time dummies? */
    ntdum = gretl_model_get_int(pmod, "ntdum");
    if (ntdum > 0) {
    /* period to which first time dummy refers */
    tdt = t1min + ifc;
    }

    yno = gretl_model_get_depvar(pmod);
    xlist = gretl_model_get_x_list(pmod);
    y = dset.Z[yno];
    b = pmod.coeff;
    t = u = row = 0;

    for (s=0; s<=dset.t2; s++) {
    double fct = NADBL;
    int missing = 0;
    int fc_skip = 0;
    int j = 0;

    if (t == 0) {
        /* fill dependent var vector for unit */
        for (ti=0; ti<T; ti++) {
        yi[ti] = y[s+ti];
        }
    }

    if (s < dset.t1) {
        /* we haven't got to the requested sample yet */
        goto next_t;
    }

    if (mask != null && mask[s] == 1) {
        /* we're doing out-of-sample forecast and this is
           an in-sample observation
        */
        fc_skip = 1;
        goto transcribe;
    }

    if (t < t1min) {
        /* required lags not available */
        goto transcribe;
    }

    /* get lagged level for integration */
    if (na(yi[t-1])) {
        missing = 1;
    } else {
        fct = yi[t-1];
    }

    /* handle autoregressive terms */
    for (i=1; i<=ylags[0] && !missing; i++) {
        k = ylags[i];
        if (t - k - 1 < 0 || na(yi[t-k]) || na(yi[t-k-1])) {
        missing = 1;
        } else {
        fct += b[j++] * (yi[t-k] - yi[t-k-1]);
        }
    }

    /* then regular regressors, differenced */
    for (i=1; i<=xlist[0] && !missing; i++) {
        vi = xlist[i];
        if (vi == 0) {
        fct += b[j++];
        } else {
        xi = dset.Z[vi];
        if (t - 1 < 0 || na(xi[s]) || na(xi[s-1])) {
            missing = 1;
        } else {
            fct += b[j++] * (xi[s] - xi[s-1]);
        }
        }
    }

    /* then a time effect, if applicable */
    if (!missing && ntdum > 0 && t >= tdt) {
        fct += b[j + t - tdt];
        if (!dpdstyle) {
        /* time dummies differenced */
        if (t - tdt > 0) {
            fct -= b[j + t - tdt - 1];
        }
        }
    }

    if (missing) {
        fct = NADBL;
    } else if (!fc_static) {
        /* dynamic: replace observed with forecast */
        yi[t] = fct;
    }

    transcribe:

    /* transcribe result? */
    if (yhat != null) {
        yhat[s] = fct;
    } else if (fc_skip) {
        ; /* skip */
    } else if (row >= fc.rows) {
        fprintf(stderr, "dpanel_fcast out of bounds! s=%d, row=%d\n".ptr, s, row);
    } else {
        transcribe_to_matrix(fc, y[s], fct, rlabels, dset, s, row++);
    }

    next_t:

    if (t == dset.pd - 1) {
        t = 0;
        u++;
    } else {
        t++;
    }
    }

    free(xlist);
    free(yi);
    if (free_ylags) {
    free(ylags);
    }

    return 0;
}

/* For pooled OLS, fixed effects and random effects: handles both
   within-sample and out of sample.
*/

int real_panel_fcast(double * yhat, GretlMatrixRaw * fc, char ** rlabels, Model * pmod, const(char) * mask, const(DatasetRaw) * dset, Gretlopt opt)
{
    const(double) * y, b;
    GretlVector * a_vec = null;
    int i, vi, s, t;
    int T, row, err = 0;

    if (individual_effects_model(pmod)) {
    a_vec = gretl_model_ahat_vec(pmod, &err);
    if (err) {
        fprintf(stderr, "real_panel_fcast: ahat vector is missing\n".ptr);
        return err;
    }
    }

    y = dset.Z[pmod.list[1]];
    b = pmod.coeff;
    s = 0; /* panel unit index */
    T = dset.pd;
    row = 0;

    for (t=0; t<=dset.t2; t++) {
    double xit, fct = NADBL;
    int missing = 0;
    int fc_skip = 0;
    int j = 0;

    if (t > 0 && t % T == 0) {
        /* advance unit index */
        s++;
    }

    if (t < dset.t1) {
        continue;
    }

    if (mask != null && mask[t] == 1) {
        /* we're doing out-of-sample forecast and this is
           an in-sample observation
        */
        fc_skip = 1;
        goto transcribe;
    }

    if (a_vec != null) {
        /* individual effect */
        if (na(a_vec.val[s])) {
        missing = 1;
        } else {
        fct = a_vec.val[s];
        }
    } else {
        fct = 0.0;
    }

    /* regular regressors */
    for (i=2, j=0; i<=pmod.list[0] && !missing; i++, j++) {
        vi = pmod.list[i];
        if (vi == 0) {
        if (!(pmod.opt & Gretlopt.OPT_F)) {
            fct += b[j];
        }
        } else {
        xit = dset.Z[vi][t];
        if (na(xit)) {
            missing = 1;
        } else {
            fct += b[j] * xit;
        }
        }
    }

    /* FIXME handling of time dummies? */

    if (missing) {
        fct = NADBL;
    }

    transcribe:

    /* transcribe result? */
    if (yhat != null) {
        yhat[t] = fct;
    } else if (fc_skip) {
        ; /* skip */
    } else if (row >= fc.rows) {
        fprintf(stderr, "panel_fcast out of bounds! t=%d, row=%d\n".ptr, t, row);
    } else {
        transcribe_to_matrix(fc, y[t], fct, rlabels, dset, t, row++);
    }
    }

    gretl_matrix_free(a_vec);

    return 0;
}

int panel_fcast(Forecast * fc, Model * pmod, DatasetRaw * dset, Gretlopt opt)
{
    int save_t1 = dset.t1;
    int save_t2 = dset.t2;
    int err;

    dset.t1 = fc.t1;
    dset.t2 = fc.t2;
    if (pmod.ci == GretlCmdIndex.dpanel) {
    err = real_dpanel_fcast(fc.yhat, null, null, pmod,
                null, dset, opt);
    } else {
    err = real_panel_fcast(fc.yhat, null, null, pmod,
                   null, dset, opt);
    }
    dset.t1 = save_t1;
    dset.t2 = save_t2;

    return err;
}

/* Out-of-sample (in the time dimension) forecast routine
   for panel data, supporting pooled OLS, fixed effects,
   random effects and dpanel. This is the top-level
   driver.
*/

int panel_os_special(Model * pmod, DatasetRaw * dset, const(char) * vname, Gretlopt opt, Prn * prn)
{
    DatasetRaw * fset = fetch_full_dataset();
    int named_save = vname[0] != '\0';
    GretlMatrixRaw * fc = null;
    double * yhat = null;
    char * mask = dset.submask;
    char ** Sr = null;
    int vi, sn;
    int err = 0;

    sn = get_dataset_submask_size(dset);
    if (sn != fset.n) {
    fprintf(stderr, "fullset->n = %d but submask size = %d\n".ptr,
        fset.n, sn);
    return GretlError.data;
    }

    if (named_save) {
    yhat = cast(typeof(yhat)) malloc(fset.n * (*yhat).sizeof);
    if (yhat == null) {
        return GretlError.alloc;
    }
    } else {
    int nfc = fset.n - dset.n;
    int cols = opt & Gretlopt.OPT_Q ? 1 : 2;

    fc = gretl_matrix_alloc(nfc, cols);
    if (fc == null) {
        return GretlError.alloc;
    }
    Sr = strings_array_new(nfc);
    }

    if (pmod.ci == GretlCmdIndex.dpanel) {
    err = real_dpanel_fcast(yhat, fc, Sr, pmod, mask, fset, opt);
    } else {
    err = real_panel_fcast(yhat, fc, Sr, pmod, mask, fset, opt);
    }

    if (!err && yhat != null) {
    err = dataset_add_NA_series(dset, 1);
    if (!err) {
            vi = dset.v - 1;
        err = rename_series(dset, vi, vname, Gretlopt.OPT_NONE);
        series_set_label(dset, vi, _("predicted values").ptr);
    }
    if (!err) {
        /* keep number of series in sync */
        err = dataset_add_allocated_series(fset, yhat);
    }
    yhat = null;
    } else if (!err) {
    gretl_matrix_set_rownames(fc, Sr);
    if (opt & Gretlopt.OPT_Q) {
        /* @fc should be suitable to purpose */
        set_fcast_matrices(fc, null);
    } else {
        if (!(opt & Gretlopt.OPT_T)) {
        /* not --stats-only */
        char ** Sc = strings_array_new(2);
        int yno = gretl_model_get_depvar(pmod);

        Sc[0] = gretl_strdup(fset.varname[yno]);
        Sc[1] = gretl_strdup(_("prediction").ptr);
        gretl_matrix_set_colnames(fc, Sc);
        gretl_matrix_print_to_prn(fc, null, prn);
        }
        if (!(opt & Gretlopt.OPT_N)) {
        /* not --no-stats */
        special_print_fc_stats(fc, prn);
        }
        /* @fc needs some adjustment */
        special_set_fcast_matrix(fc);
    }
    fc = null;
    }

    gretl_matrix_free(fc);
    free(yhat);

    return err;
}

/* Generate forecasts from nonlinear least squares model, using the
   string specification of the regression function that was saved as
   data on the model (see nls.c).  If the NLS formula is dynamic and
   the user has not requested a static forecast, we do an
   autoregressive genr out of sample.  If we're doing a static
   forecast we add a simple-minded forecast error, namely the standard
   error of the NLS regression.
*/

int nls_fcast(Forecast * fc, const(Model) * pmod, DatasetRaw * dset)
{
    int oldt1 = dset.t1;
    int oldt2 = dset.t2;
    int origv = dset.v;
    int yno = 0;
    const(char) * nlfunc;
    double * y = null;
    char[MAXLINE] formula;
    int t, err = 0;

    nlfunc = cast(typeof(nlfunc)) gretl_model_get_data(pmod, "nl_regfunc");
    if (nlfunc == null) {
    err = GretlError.data;
    }

    if (fc.method == FcastMethod.auto_) {
    yno = pmod.list[1];
    y = copyvec(dset.Z[yno], dset.n);
    if (y == null) {
        err = GretlError.alloc;
    }
    }

    /* The computation of $nl_y below may depend on computing
       some auxiliary quantities used within the nls block, so
       try that first below: nl_model_run_aux_genrs().
    */

    if (!err) {
    /* upper limit of static forecast */
    int t2 = (fc.method == FcastMethod.static_)? fc.t2 : pmod.t2;

    if (t2 >= fc.t1) {
        /* non-null static range */
        int fcv;

        dset.t1 = fc.t1;
        dset.t2 = t2;
        err = nl_model_run_aux_genrs(pmod, dset);
        fcv = dset.v;
        if (!err) {
        sprintf(formula.ptr, "$nl_y = %s".ptr, nlfunc);
        err = generate(formula.ptr, dset, GretlType.series,
                   Gretlopt.OPT_P, null);
        }
        if (!err) {
        for (t=dset.t1; t<=dset.t2; t++) {
            fc.yhat[t] = dset.Z[fcv][t];
        }
        }
    }

    if (!err && fc.method == FcastMethod.auto_ && fc.t2 > pmod.t2) {
        /* dynamic forecast out of sample: in this context
           pmod->depvar is the expression to generate
           the series in question under its own name
        */
        dset.t1 = pmod.t2 + 1;
        dset.t2 = fc.t2;
        err = nl_model_run_aux_genrs(pmod, dset);
        if (!err) {
        strcpy(formula.ptr, pmod.depvar);
        err = generate(formula.ptr, dset, GretlType.series,
                   Gretlopt.OPT_P, null);
        }
        if (!err) {
        for (t=dset.t1; t<=dset.t2; t++) {
            fc.yhat[t] = dset.Z[yno][t];
        }
        }
    }
    }

    /* restore dataset state */
    if (dset.v > origv) {
    err = dataset_drop_last_variables(dset, dset.v - origv);
    }
    dset.t1 = oldt1;
    dset.t2 = oldt2;

    if (y != null) {
    /* restore original dependent variable */
    memcpy(dset.Z[yno], y, dset.n * (*y).sizeof);
    free(y);
    }

    if (!err && fc.method == FcastMethod.static_ && fc.sderr != null) {
    /* kinda simple-minded */
    for (t=fc.t1; t<=fc.t2; t++) {
        fc.sderr[t] = pmod.sigma;
    }
    }

    return err;
}

int midas_fcast(Forecast * fc, const(Model) * pmod, DatasetRaw * dset)
{
    int oldt1 = dset.t1;
    int oldt2 = dset.t2;
    int origv = dset.v;
    char[MAXLINE] formula;
    char * mdsfunc = null;
    double * y = null;
    int yno = pmod.list[1];
    int t, err;

    err = midas_forecast_setup(pmod, dset, fc.method, &mdsfunc);

    if (!err && fc.method == FcastMethod.auto_) {
    y = copyvec(dset.Z[yno], dset.n);
    if (y == null) {
        err = GretlError.alloc;
    }
    }

    if (!err) {
    /* upper limit of static forecast */
    int t2 = (fc.method == FcastMethod.static_)? fc.t2 : pmod.t2;

    if (t2 >= fc.t1) {
        /* non-null static range */
        int fcv = dset.v;

        dset.t1 = fc.t1;
        dset.t2 = t2;
        sprintf(formula.ptr, "$midas_y=%s".ptr, mdsfunc);
        err = generate(formula.ptr, dset, GretlType.series,
               Gretlopt.OPT_P, null);
        if (!err) {
        for (t=dset.t1; t<=dset.t2; t++) {
            fc.yhat[t] = dset.Z[fcv][t];
        }
        }
    }

    if (!err && fc.method == FcastMethod.auto_ && fc.t2 > pmod.t2) {
        /* dynamic forecast out of sample */
        dset.t1 = pmod.t2 + 1;
        dset.t2 = fc.t2;
        sprintf(formula.ptr, "%s=%s".ptr, dset.varname[yno], mdsfunc);
        err = generate(formula.ptr, dset, GretlType.series,
               Gretlopt.OPT_P, null);
        if (!err) {
        for (t=dset.t1; t<=dset.t2; t++) {
            fc.yhat[t] = dset.Z[yno][t];
        }
        }
    }
    }

    /* restore dataset state */
    if (dset.v > origv) {
    err = dataset_drop_last_variables(dset, dset.v - origv);
    }
    dset.t1 = oldt1;
    dset.t2 = oldt2;

    if (y != null) {
    /* restore original dependent variable */
    memcpy(dset.Z[yno], y, dset.n * (*y).sizeof);
    free(y);
    }

    if (!err && fc.method == FcastMethod.static_ && fc.sderr != null) {
    for (t=fc.t1; t<=fc.t2; t++) {
        fc.sderr[t] = pmod.sigma;
    }
    }

    /* clean up any uservars the MIDAS mechanism created */
    destroy_private_uvars();

    free(mdsfunc);

    return err;
}

int depvar_lag(const(Forecast) * f, int i)
{
    return (f.dvlags != null)? f.dvlags[i] : 0;
}

/*
   GARCH error variance process:

   h(t) = a(0) + sum(i=1 to q) a(i) * u(t-i)^2
                   + sum(j=1 to p) b(j) * h(t-j)

   This is then complexified if the model includes lags of the
   dependent variable among the regressors, but in the following
   function we just calculate the successive h's.
*/

double * garch_h_hat(const(Model) * pmod, int t1, int t2, int xvars)
{
    const(double) * alpha;
    const(double) * beta;
    double hlag;
    double * h, mh;
    int nf, q, p, lmax;
    int i, s, t, ti;
    int err = 0;

    mh = cast(typeof(mh)) gretl_model_get_data(pmod, "garch_h");
    if (mh == null) {
    return null;
    }

    nf = t2 - t1 + 1;
    if (nf <= 0) {
    return null;
    }

    h = cast(typeof(h)) malloc(nf * (*h).sizeof);
    if (h == null) {
    return null;
    }

    q = pmod.list[1];
    p = pmod.list[2];
    alpha = pmod.coeff + xvars;
    beta = alpha + q + 1;

    lmax = max(p, q);

    for (t=t1; t<=t2 && !err; t++) {
    s = t - t1;

    h[s] = alpha[0];
    for (i=1; i<=lmax; i++) {
        ti = t - i;
        if (ti < 0) {
        break;
        }
        hlag = 0.0;
        if (ti <= pmod.t2) {
        if (!na(mh[ti])) {
            hlag = mh[ti];
        }
        } else {
        hlag = h[s-i];
        }
        if (i <= q) {
        if (ti <= pmod.t2 && !na(pmod.uhat[ti])) {
            h[s] += alpha[i] * pmod.uhat[ti] * pmod.uhat[ti];
        } else {
            /* lagged uhat^2 not available: substitute its
               expectation, h at the same date */
            h[s] += alpha[i] * hlag;
        }
        }
        if (i <= p) {
        h[s] += beta[i-1] * hlag;
        }
    }
    if (h[s] < 0.0) {
        err = 1;
    }
    }

    if (err) {
    free(h);
    h = null;
    }

    return h;
}

/* Infinite MA representation, in case GARCH model has lags of the
   dependent var among the regressors (the autoregressive coeffs are
   recorded in the array phi).
*/

double * garch_psi(const(double) * phi, int p, int nf)
{
    double * psi;
    int i, s;

    if (phi == null) {
    return null;
    }

    psi = cast(typeof(psi)) malloc(nf * (*psi).sizeof);
    if (psi == null) {
    return null;
    }

    psi[0] = 1.0;

    /* Are we OK with regard to signs below?  Do we need to keep a
       record of the squares of the psi elements too? (Or instead?)
    */

    for (s=1; s<nf; s++) {
    psi[s] = 0.0;

    /* add psi[s] components derived from dep var lags */
    for (i=0; i<p && i<s; i++) {
        psi[s] += phi[i] * psi[s-i-1];
    }
    }

    return psi;
}

/* Construct an array, "phi", of autoregressive coefficients, if a
   GARCH model has any lags of the dependent variable among the
   regressors; the array has dimension equal to the highest-order lag
   of the dependent variable.
*/

double * garch_ldv_phi(const(Model) * pmod, const(int) * dvlags, int * ppmax)
{
    double * phi = null;
    int i, pmax = 0;

    if (dvlags != null && pmod.xlist != null) {
    for (i=0; i<pmod.xlist[0]; i++) {
        if (dvlags[i] > pmax) {
        pmax = dvlags[i];
        }
    }

    phi = cast(typeof(phi)) malloc(pmax * (*phi).sizeof);
    if (phi != null) {
        for (i=0; i<pmax; i++) {
        phi[i] = 0.0;
        }
        for (i=0; i<pmod.xlist[0]; i++) {
        if (dvlags[i] > 0) {
            phi[dvlags[i] - 1] = pmod.coeff[i];
        }
        }
    }
    }

    if (phi != null) {
    *ppmax = pmax;
    }

    return phi;
}

/* Compute forecast standard errors for a GARCH model that includes
   lags of the dependent variable among the regressors.  This could do
   with some scrutiny (or the calculation, above, of the "psi" array
   could do with checking).
*/

double garch_ldv_sderr(const(double) * h, const(double) * psi, int s)
{
    double ss, vs = h[s];
    int i;

    for (i=1; i<=s; i++) {
    vs += h[s-i] * psi[i] * psi[i];
    }

    if (vs >= 0.0) {
    ss = sqrt(vs);
    } else {
    ss = NADBL;
    }

    return ss;
}

int garch_fcast(Forecast * fc, Model * pmod, const(DatasetRaw) * dset)
{
    double xval;
    int xvars, yno;
    double * h = null;
    double * mh = null;
    double * phi = null;
    double * psi = null;
    int i, v, t;

    yno = gretl_model_get_depvar(pmod);

    if (pmod.xlist != null) {
    xvars = pmod.xlist[0];
    } else {
    xvars = 0;
    }

    if (fc.sderr != null) {
    h = garch_h_hat(pmod, pmod.t2 + 1, fc.t2, xvars);
    mh = cast(typeof(mh)) gretl_model_get_data(pmod, "garch_h");
    }

    if (fc.dvlags != null) {
    int ns = fc.t2 - pmod.t2;
    int pmax = 0;

    if (ns > 0) {
        phi = garch_ldv_phi(pmod, fc.dvlags, &pmax);
        psi = garch_psi(phi, pmax, ns);
        free(phi);
    }
    }

    for (t=fc.t1; t<=fc.t2; t++) {
    int lag, miss = 0;
    double yh = 0.0;

    if (fc.method != FcastMethod.dynamic &&
        t >= pmod.t1 && t <= pmod.t2) {
        fc.yhat[t] = pmod.yhat[t];
        if (fc.sderr != null && mh != null) {
        fc.sderr[t] = mh[t];
        }
        continue;
    }

    for (i=0; i<xvars; i++) {
        v = pmod.xlist[i+1];
        if ((lag = depvar_lag(fc, i)) != 0) {
        xval = fcast_get_ldv(fc, yno, t, lag, dset);
        } else {
        xval = dset.Z[v][t];
        }
        if (na(xval)) {
        miss = 1;
        } else {
        yh += pmod.coeff[i] * xval;
        }
    }

    if (miss) {
        fc.yhat[t] = NADBL;
    } else {
        fc.yhat[t] = yh;
    }

    if (h != null) {
        if (t > pmod.t2) {
        if (psi != null) {
            /* build in effect of lagged dependent var */
            fc.sderr[t] = garch_ldv_sderr(h, psi, t - pmod.t2 - 1);
        } else {
            /* no lagged dependent variable */
            fc.sderr[t] = sqrt(h[t - pmod.t2 - 1]);
        }
        } else {
        fc.sderr[t] = NADBL;
        }
    }
    }

    if (h != null) {
    free(h);
    }
    if (psi != null) {
    free(psi);
    }

    return 0;
}

/* Compute ARMA forecast error variance (ignoring parameter
   uncertainty, as is common), via recursion. Cf. Box and Jenkins,
   1976, p. 508, "Program 4", V(l) algorithm (with the sign of theta
   changed).
*/

double arma_variance(const(double) * phi, int p, const(double) * theta, int q, double * psi, int npsi, int l)
{
    /* the sum of squared psi's */
    static double sspsi;

    if (l == 1) {
    int i, j;

    psi[0] = 1.0;
    for (j=1; j<npsi; j++) {
        psi[j] = 0.0;
        for (i=1; i<=j; i++) {
        if (i <= p) {
            psi[j] += phi[i] * psi[j-i];
        }
        }
        if (theta != null && j <= q) {
        psi[j] += theta[j];
        }
    }
    sspsi = 1.0;
    } else {
    sspsi += psi[l-1] * psi[l-1];
    }

    return sspsi;
}

double arima_difference_obs(const(double) * x, int t, int * c, int k)
{
    double dx = x[t];

    if (!na(dx)) {
    int i, s;

    for (i=0; i<k; i++) {
        if (c[i] != 0) {
        s = t - i - 1;
        if (s < 0 || na(x[s])) {
            dx = NADBL;
            break;
        } else {
            dx -= c[i] * x[s];
        }
        }
    }
    }

    return dx;
}

/* When forecasting based on an armax model estimated using X12A,
   or via the Kalman filter, we need to form the series X\beta
   so that we can subtract X\beta_{t-i} from y_{t-i} in
   computing the AR portion of the forecast. "beta" below
   is the array of ARMAX coefficients, not including the
   constant.

   If the model is ARIMAX and was NOT estimated with the
   --y-diff-only option (OPT_Y) then the regressors X have
   to be differenced in this context.
*/

double * create_Xb_series(Forecast * fc, const(Model) * pmod, const(double) * beta, const(double) ** Z)
{
    double * Xb;
    int * delta = null;
    double x;
    int miss, diff = 0;
    int d, D, s = 0, k = 0;
    int i, j, t, vi;

    Xb = cast(typeof(Xb)) malloc((fc.t2 + 1) * (*Xb).sizeof);
    if (Xb == null) {
    return null;
    }

    d = gretl_model_get_int(pmod, "arima_d");
    D = gretl_model_get_int(pmod, "arima_D");

    if ((d > 0 || D > 0) && !(pmod.opt & Gretlopt.OPT_Y)) {
    s = gretl_model_get_int(pmod, "arma_pd");
    diff = 1;
    }

    if (diff) {
    delta = arima_delta_coeffs(d, D, s);
    if (delta == null) {
        free(Xb);
        return null;
    }
    k = d + s * D;
    }

    for (t=0; t<=fc.t2; t++) {
    Xb[t] = 0.0;
    miss = 0;
    j = 0;
    for (i=1; i<=pmod.xlist[0] && !miss; i++) {
        vi = pmod.xlist[i];
        if (vi == 0) {
        Xb[t] += pmod.coeff[0];
        } else {
        if (diff) {
            x = arima_difference_obs(Z[vi], t, delta, k);
        } else {
            x = Z[vi][t];
        }
        if (na(x)) {
            Xb[t] = NADBL;
            miss = 1;
        } else {
            Xb[t] += beta[j++] * x;
        }
        }
    }
    }

    if (delta != null) {
    free(delta);
    }

    return Xb;
}

int want_x_beta_prep(const(Model) * pmod)
{
    int ret = 0;

    if (pmod.xlist != null) {
    int aflags = gretl_model_get_int(pmod, "arma_flags");

    if (aflags & (ArmaFlags.exact | ArmaFlags.x12a)) {
        ret = 1;
    }
    }

    return ret;
}

/* generate forecasts for AR(I)MA (or ARMAX) models, including
   forecast standard errors if we're doing out-of-sample
   forecasting
*/

int arma_fcast(Forecast * fc, Model * pmod, const(DatasetRaw) * dset)
{
    double * psi = null;
    GretlVector * phi = null;
    GretlVector * theta = null;
    double * phi0 = null;
    double * Xb = null;
    const(double) * beta;
    const(double) * y;
    double xval, yval, vl;
    double mu = NADBL;
    int xvars, yno;
    int p, q, px = 0, npsi = 0;
    int fcstart = fc.t1;
    int ar_smax, ma_smax;
    int regarma = 0;
    int i, s, t;
    int err = 0;

    if (fc.method != FcastMethod.dynamic) {
    /* use pre-calculated fitted values over model estimation range,
       and don't bother calculating forecast error variance */
    for (t=fc.t1; t<=pmod.t2; t++) {
        fc.yhat[t] = pmod.yhat[t];
        if (fc.sderr != null) {
        fc.sderr[t] = NADBL;
        }
    }
    if (fc.t2 <= pmod.t2) {
        /* no "real" forecasts were called for, we're done */
        return 0;
    }
    fcstart = pmod.t2 + 1;
    }

    p = arma_model_max_AR_lag(pmod);
    q = arma_model_max_MA_lag(pmod);

    yno = gretl_model_get_depvar(pmod);

    xvars = (pmod.xlist != null)? pmod.xlist[0] : 0;

    err = arma_model_AR_MA_coeffs(pmod, &phi, &theta, Gretlopt.OPT_I);
    if (err) {
    goto bailout;
    }

    beta = arma_model_get_x_coeffs(pmod);

    if (want_x_beta_prep(pmod)) {
    if (gretl_is_arima_model(pmod)) {
        regarma = 1;
        err = regarma_model_AR_coeffs(pmod, &phi0, &px);
    }
    if (pmod.xlist[0] == 1 && pmod.xlist[1] == 0) {
        /* just a const, no ARMAX */
        if (phi0 != null) {
        mu = 1.0;
        for (i=1; i<=px; i++) {
            mu -= phi0[i];
        }
        mu *= pmod.coeff[0];
        } else {
        mu = pmod.coeff[0];
        }
    } else {
        /* we have ARMAX terms */
        Xb = create_Xb_series(fc, pmod, beta,
                  cast(const(double) **) dset.Z);
        if (Xb == null) {
        err = GretlError.alloc;
        }
    }
    if (err) {
        goto bailout;
    }
    }

    /* setup for forecast error variance */
    if (fc.sderr != null) {
    npsi = fc.t2 - fcstart + 1;
    psi = cast(typeof(psi)) malloc(npsi * (*psi).sizeof);
    }

    /* cut-off points for using actual rather than forecast
       values of y in generating further forecasts (FIXME??) */
    if (fc.method == FcastMethod.static_) {
    ar_smax = fc.t2;
    ma_smax = pmod.t2;
    } else if (fc.method == FcastMethod.dynamic) {
    ar_smax = max(fc.t1 - 1, p - 1);
    ma_smax = max(fc.t1 - 1, q);
    } else {
    ar_smax = pmod.t2;
    ma_smax = pmod.t2;
    }

    /* dependent variable */
    y = dset.Z[yno];

    /* do real forecast */
    for (t=fcstart; t<=fc.t2 && !err; t++) {
    double yh = 0.0;
    int miss = 0;

    /* contribution of const and/or independent variables */

    if (!na(mu)) {
        yh = mu;
    } else if (Xb != null) {
        /* X\beta series is pre-computed */
        if (na(Xb[t])) {
        miss = 1;
        } else {
        yh = Xb[t];
        }
    } else {
        int j = 0;

        for (i=1; i<=xvars; i++) {
        if (pmod.xlist[i] == 0) {
            yh += pmod.coeff[0];
        } else {
            xval = dset.Z[pmod.xlist[i]][t];
            if (na(xval)) {
            miss = 1;
            } else {
            yh += beta[j++] * xval;
            }
        }
        }
    }

    if (phi0 != null && Xb != null) {
        /* the regarma case */
        for (i=1; i<=px; i++) {
        s = t - i;
        if (s < 0) {
            miss = 1;
        } else if (na(Xb[s])) {
            miss = 1;
        } else {
            yh -= phi0[i] * Xb[s];
        }
        }
    }

    /* AR contribution (incorporating any differencing) */

    for (i=1; i<=p && !miss; i++) {
        if (phi.val[i] == 0.0) {
        continue;
        }
        s = t - i;
        if (s < 0) {
        yval = NADBL;
        } else if (s <= ar_smax) {
        yval = y[s];
        } else {
        yval = fc.yhat[s];
        }
        if (na(yval)) {
        miss = 1;
        } else {
        if (!regarma && Xb != null) {
            if (na(Xb[s])) {
            miss = 1;
            } else {
            yh += phi.val[i] * (yval - Xb[s]);
            }
        } else if (!regarma && !na(mu)) {
            yh += phi.val[i] * (yval - mu);
        } else {
            yh += phi.val[i] * yval;
        }
        }
    }

    /* MA contribution */

    for (i=1; i<=q && !miss; i++) {
        if (theta.val[i] == 0.0) {
        continue;
        }
        s = t - i;
        if (s >= pmod.t1 && s <= ma_smax) {
        if (na(pmod.uhat[s])) {
            char[OBSLEN] obsstr;

            ntolabel(obsstr.ptr, s, dset);
            fprintf(stderr, "Uh oh, pmod->uhat[%d] is missing! (%s)\n".ptr,
                s, obsstr.ptr);
        } else {
            yh += theta.val[i] * pmod.uhat[s];
        }
        } else if (fc.eps != null) {
        yh += theta.val[i] * fc.eps[s];
        }
    }

    if (miss) {
        fc.yhat[t] = NADBL;
    } else {
        fc.yhat[t] = yh;
        if (fc.eps != null) {
        /* form estimated error: we do this only in
           the case of a "static" forecast */
        fc.eps[t] = y[t] - yh;
        }
    }

    /* forecast error variance */
    if (psi != null) {
        const(double) * arvec = (phi != null)? phi.val : null;
        const(double) * mavec = (theta != null)? theta.val : null;

        vl = arma_variance(arvec, p, mavec, q,
                   psi, npsi, t - fcstart + 1);
        fc.sderr[t] = pmod.sigma * sqrt(vl);
    }
    }

 bailout:

    free(psi);
    free(Xb);
    free(phi0);

    gretl_vector_free(phi);
    gretl_vector_free(theta);

    return err;
}

/* construct the "phi" array of AR coefficients, based on the ARINFO
   that was added to the model at estimation time.  The latter's rho
   member may be a compacted array, with zero elements omitted, but
   here we need a full-length array with zeros inserted as required.
*/

double * make_phi_from_arinfo(const(ArInfo) * arinfo, int pmax)
{
    double * phi = cast(double *) malloc((pmax + 1) * (double).sizeof);

    if (phi != null) {
    int i, lag;

    for (i=0; i<=pmax; i++) {
        phi[i] = 0.0;
    }

    for (i=1; i<=arinfo.arlist[0]; i++) {
        lag = arinfo.arlist[i];
        phi[lag] = arinfo.rho[i-1];
    }
    }

    return phi;
}

/* Determine the greatest lag order for a model, either via explicit
   AR error process or via inclusion of lagged dependent var as
   regressor.
*/

int max_ar_lag(Forecast * fc, const(Model) * pmod, int p)
{
    int i, pmax = p; /* explicit AR order */

    if (fc.dvlags != null) {
    for (i=0; i<pmod.ncoeff; i++) {
        if (fc.dvlags[i] > pmax) {
        pmax = fc.dvlags[i];
        }
    }
    }

    return pmax;
}

/* Set things up for computing forecast error variance for ar models
   (AR, AR1).  This is complicated by the fact that there
   may be a lagged dependent variable in the picture.  If there is,
   the effective AR coefficients have to be incremented, for the
   purposes of calculating forecast variance.  But I'm not sure this
   is quite right yet.
*/

int set_up_ar_fcast_variance(const(Model) * pmod, int pmax, int npsi, double ** pphi, double ** ppsi, double ** perrphi)
{
    double * errphi = null;
    double * psi = null;
    double * phi = null;
    int err = 0;

    errphi = make_phi_from_arinfo(pmod.arinfo, pmax);
    psi = cast(typeof(psi)) malloc(npsi * (*psi).sizeof);
    phi = cast(typeof(phi)) malloc((pmax + 1) * (*phi).sizeof);

    if (errphi == null || psi == null || phi == null) {
    free(errphi);
    free(psi);
    free(phi);
    err = GretlError.alloc;
    } else {
    *perrphi = errphi;
    *pphi = phi;
    *ppsi = psi;
    }

    return err;
}

/*
   The code below generates forecasts that incorporate the
   predictable portion of an AR error term:

       u_t = r1 u_{t-1} + r2 u_{t-1} + ... + e_t

   where e_t is white noise.  The forecasts are based on the
   representation of a model with such an error term as

       (1 - r(L)) y_t = (1 - r(L)) X_t b + e_t

   or

       y_t = r(L) y_t + (1 - r(L)) X_t b + e_t

   where r(L) is a polynomial in the lag operator.

   We also attempt to calculate forecast error variance for
   out-of-sample forecasts.  These calculations, like those for
   ARMA, do not take into account parameter uncertainty.

   This code is used for AR and AR1 models; it is also used for
   dynamic forecasting with models that do not have an explicit AR
   error process but that have one or more lagged values of the
   dependent variable as regressors.
*/

int ar_fcast(Forecast * fc, Model * pmod, const(DatasetRaw) * dset)
{
    const(int) * arlist;
    double * phi = null;
    double * psi = null;
    double * errphi = null;
    double * rho;
    double xval, yh, vl;
    double ylag, xlag;
    int miss, yno;
    int i, k, v, t, tk;
    int p, dvlag, pmax = 0;
    int pwe, npsi = 0;
    int err = 0;

    yno = pmod.list[1];

    if (pmod.ci == GretlCmdIndex.ar1 && pmod.arinfo == null) {
    dummy_ar_info_init(pmod);
    }

    arlist = pmod.arinfo.arlist;
    p = arlist[arlist[0]]; /* AR order of error term */
    rho = pmod.arinfo.rho;

    if (fc.t2 > pmod.t2 && fc.sderr != null) {
    /* we'll compute variance only if we're forecasting out of
       sample */
    pmax = max_ar_lag(fc, pmod, p);
    npsi = fc.t2 - fc.t1 + 1;
    /* npsi = fc->t2 - pmod->t2; */
    set_up_ar_fcast_variance(pmod, pmax, npsi, &phi, &psi, &errphi);
    }

    pwe = (pmod.opt & Gretlopt.OPT_P);

    for (t=fc.t1; t<=fc.t2; t++) {
    miss = 0;
    yh = 0.0;

    if (t < p) {
        fc.yhat[t] = NADBL;
        if (fc.sderr != null) {
        fc.sderr[t] = NADBL;
        }
        continue;
    }

        if (pwe && t == pmod.t1) {
            /* PWE first obs is special */
            fc.yhat[t] = pmod.yhat[t];
            continue;
        }

    if (phi != null) {
        /* initialize the phi's based on the AR error process
           alone */
        for (i=1; i<=pmax; i++) {
        phi[i] = errphi[i];
        }
    }

    /* r(L) y_t */
    for (k=1; k<=arlist[0]; k++) {
        tk = t - arlist[k];
        ylag = fcast_get_ldv(fc, yno, tk, 0, dset);
        if (na(ylag)) {
        miss = 1;
        } else {
        yh += rho[k-1] * ylag;
        }
    }

    /* (1 - r(L)) X_t b */
    for (i=0; i<pmod.ncoeff && !miss; i++) {
        v = pmod.list[i+2];
        if ((dvlag = depvar_lag(fc, i)) != 0) {
        xval = fcast_get_ldv(fc, yno, t, dvlag, dset);
        } else {
        xval = dset.Z[v][t];
        }
        if (na(xval)) {
        miss = 1;
        } else {
        if (dvlag > 0 && phi != null) {
            /* augment phi for computation of variance */
            phi[dvlag] += pmod.coeff[i];
        }
        for (k=1; k<=arlist[0]; k++) {
            tk = t - arlist[k];
            if (dvlag > 0) {
            xlag = fcast_get_ldv(fc, yno, tk, dvlag, dset);
            } else {
            xlag = dset.Z[v][tk];
            }
            if (!na(xlag)) {
            xval -= rho[k-1] * xlag;
            }
        }
        yh += xval * pmod.coeff[i];
        }
    }

    if (miss) {
        fc.yhat[t] = NADBL;
    } else {
        fc.yhat[t] = yh;
    }

    /* forecast error variance */
    if (phi != null && pmod.ci != GretlCmdIndex.garch) {
        if (t > pmod.t2) {
        vl = arma_variance(phi, pmax, null, 0, psi, npsi, t - pmod.t2);
        fc.sderr[t] = pmod.sigma * sqrt(vl);
        } else {
        fc.sderr[t] = NADBL;
        }
    }
    }

    if (psi != null) {
    free(phi);
    free(psi);
    free(errphi);
    }

    return err;
}

/* Calculates the transformation required to get from xb (= X*b) to
   the actual prediction for the dependent variable, for models of
   type LOGISTIC, LOGIT, PROBIT, TOBIT, POISSON and NEGBIN.
 */

double fcast_transform(double xb, const(Model) * pmod, int t, const(double) * offset, double lmax)
{
    int ymin, ci = pmod.ci;
    double yf = xb;

    if (ci == GretlCmdIndex.tobit) {
    if (xb < 0.0) {
        yf = 0.0;
    }
    } else if (ci == GretlCmdIndex.logit) {
    if (gretl_model_get_int(pmod, "ordered")) {
        ymin = gretl_model_get_int(pmod, "ymin");
        yf = ordered_model_prediction(pmod, xb, ymin);
    } else {
        yf = exp(xb) / (1.0 + exp(xb));
    }
    } else if (ci == GretlCmdIndex.probit) {
    if (gretl_model_get_int(pmod, "ordered")) {
        ymin = gretl_model_get_int(pmod, "ymin");
        yf = ordered_model_prediction(pmod, xb, ymin);
    } else {
        yf = normal_cdf(xb);
    }
    } else if (ci == GretlCmdIndex.logistic) {
    if (na(lmax)) {
        yf = 1.0 / (1.0 + exp(-xb));
    } else {
        yf = lmax / (1.0 + exp(-xb));
    }
    } else if (ci == GretlCmdIndex.poisson || ci == GretlCmdIndex.negbin) {
    if (offset != null) {
        if (na(offset[t])) {
        yf = NADBL;
        } else {
        yf = exp(xb + log(offset[t]));
        }
    } else {
        yf = exp(xb);
    }
    }

    return yf;
}

int integrated_fcast(Forecast * fc, const(Model) * pmod, int yno, const(DatasetRaw) * dset)
{
    double s2 = NADBL;
    int dyn_start; /* start of dynamic forecast */
    double xval, yht = 0.0;
    int i, vi, t;

    if (fc.method == FcastMethod.static_) {
    dyn_start = dset.n; /* i.e., never */
    } else if (fc.method == FcastMethod.dynamic) {
    dyn_start = fc.t1;
    } else {
    /* FC_AUTO: dynamic out of sample */
    dyn_start = pmod.t2 + 1;
    }

    if (fc.sderr != null) {
    s2 = pmod.sigma * pmod.sigma;
    fc.sderr[dyn_start - 1] = 0.0;
    }

    for (t=fc.t1; t<=fc.t2; t++) {
    int miss = 0;

    if (t <= dyn_start) {
        yht = dset.Z[yno][t-1];
        if (na(yht)) {
        miss = 1;
        }
    }

    for (i=0; i<pmod.ncoeff && !miss; i++) {
        vi = pmod.list[i+2];
        xval = dset.Z[vi][t];
        if (na(xval)) {
        miss = 1;
        } else {
        yht += xval * pmod.coeff[i];
        }
    }

    if (miss) {
        if (t >= dyn_start) {
        /* all subsequent values are NA */
        int s;

        for (s=t; s<=fc.t2; s++) {
            fc.yhat[s] = NADBL;
        }
        break;
        } else {
        fc.yhat[t] = NADBL;
        }
    } else {
        fc.yhat[t] = yht;
        if (fc.sderr != null && t >= dyn_start) {
        fc.sderr[t] = fc.sderr[t-1] + s2;
        }
    }
    }

    if (fc.sderr != null) {
    fc.sderr[dyn_start - 1] = NADBL;
    for (t=dyn_start; t<=fc.t2; t++) {
        if (!na(fc.sderr[t])) {
        fc.sderr[t] = sqrt(fc.sderr[t]);
        }
    }
    }

    return 0;
}

int mlogit_allprobs(Forecast * fc, const(Model) * pmod, const(DatasetRaw) * dset)
{
    GretlMatrixRaw * P;
    int err = 0;

    P = mn_logit_probabilities(pmod, fc.t1, fc.t2, dset, &err);
    if (P != null) {
    set_fcast_matrices(P, null);
    }

    return err;
}

int mlogit_fcast(Forecast * fc, const(Model) * pmod, const(DatasetRaw) * dset)
{
    const(GretlMatrixRaw) * yvals;
    GretlMatrixRaw * Xt;
    int k, i, vi, t;

    yvals = cast(typeof(yvals)) gretl_model_get_data(pmod, "yvals");
    if (yvals == null) {
    return GretlError.data;
    }

    k = pmod.list[0] - 1;
    Xt = gretl_matrix_alloc(1, k);
    if (Xt == null) {
    return GretlError.alloc;
    }

    for (t=fc.t1; t<=fc.t2; t++) {
    double xti, yht;
    int miss = 0;

    for (i=0; i<k && !miss; i++) {
        vi = pmod.list[i+2];
        xti = dset.Z[vi][t];
        if (na(xti)) {
        miss = 1;
        } else {
        Xt.val[i] = xti;
        }
    }

    if (miss) {
        fc.yhat[t] = NADBL;
    } else {
        yht = mn_logit_prediction(Xt, pmod.coeff, yvals);
        fc.yhat[t] = yht;
    }
    }

    gretl_matrix_free(Xt);

    return 0;
}

/* Compute forecasts for linear models without autoregressive errors.
   We don't compute forecast standard errors, unless we're integrating
   the forecast from a static model estimated via OLS.
*/

int linear_fcast(Forecast * fc, const(Model) * pmod, int yno, const(DatasetRaw) * dset, const(GretlMatrixRaw) * X, Gretlopt opt)
{
    const(double) * offvar = null;
    double xval, lmax = NADBL;
    int all_probs = 0;
    int k = pmod.ncoeff;
    int i, vi, t;
    int err = 0;

    if (COUNT_MODEL(cast(GretlCmdIndex) pmod.ci) && dset != null) {
    /* special for "offset" variable */
    int offnum = gretl_model_get_int(pmod, "offset_var");

    if (offnum > 0) {
        offvar = dset.Z[offnum];
    }
    } else if (pmod.ci == GretlCmdIndex.logistic) {
    lmax = gretl_model_get_double(pmod, "lmax");
    } else if (pmod.ci == GretlCmdIndex.logit || pmod.ci == GretlCmdIndex.probit) {
    if (gretl_model_get_int(pmod, "ordered")) {
        /* we need to know how many coeffs are not
           just estimated cut-points */
        k = gretl_model_get_int(pmod, "nx");
        if (k <= 0) {
        return GretlError.missdata;
        } else if (opt & Gretlopt.OPT_L) {
        all_probs = 1;
        }
    }
    }

    for (t=fc.t1; t<=fc.t2; t++) {
    double yht = 0.0;
    int miss = 0;

    for (i=0; i<k && !miss; i++) {
        int lag = depvar_lag(fc, i);

        if (lag != 0) {
        xval = fcast_get_ldv(fc, yno, t, lag, dset);
        } else if (dset != null) {
                vi = pmod.list[i+2];
        xval = dset.Z[vi][t];
        } else {
                xval = X.val[i * X.rows + t];
            }
        if (na(xval)) {
        miss = 1;
        } else {
        yht += xval * pmod.coeff[i];
        }
    }

    if (miss) {
        fc.yhat[t] = NADBL;
    } else if (FCAST_SPECIAL(pmod.ci) && !all_probs) {
        /* special handling for LOGIT and others */
        fc.yhat[t] = fcast_transform(yht, pmod, t, offvar, lmax);
    } else {
        fc.yhat[t] = yht;
    }
    }

    if (all_probs && dset != null) {
    GretlMatrixRaw * P;

    P = ordered_probabilities(pmod, fc.yhat, fc.t1, fc.t2,
                  dset, &err);
    set_fcast_matrices(P, null);
    }

    return err;
}

bool dynamic_nls(const(Model) * m)
{
    return (m.ci == GretlCmdIndex.nls || m.ci == GretlCmdIndex.midasreg)
        && gretl_model_get_int(m, "dynamic");
}

int get_forecast_method(Forecast * fc, Model * pmod, const(DatasetRaw) * dset, Gretlopt opt)
{
    int dyn_ok = 0;
    int dyn_errs_ok = 0;
    int err;

    err = incompatible_options(opt, Gretlopt.OPT_D | Gretlopt.OPT_S);
    if (err) {
    return err;
    }

    fc.dvlags = null;
    fc.method = FcastMethod.static_;

    /* do setup for possible lags of the dependent variable,
       unless OPT_S for "static" has been given
    */
    if (dataset_is_time_series(dset) && !(opt & Gretlopt.OPT_S) && pmod.ci != GretlCmdIndex.arma) {
    fc.dvlags = process_lagged_depvar(pmod, dset);
    }

    if (!(opt & Gretlopt.OPT_S)) {
    /* user didn't give the "static" option */
    if (pmod.ci == GretlCmdIndex.arma || fc.dvlags != null ||
        dynamic_nls(pmod) || (opt & Gretlopt.OPT_I)) {
        /* dynamic forecast is possible */
        dyn_ok = 1;
    }
    if (SIMPLE_AR_MODEL(cast(GretlCmdIndex) pmod.ci) || pmod.ci == GretlCmdIndex.garch) {
        dyn_errs_ok = 1;
    }
    }

    if (!dyn_ok && (opt & Gretlopt.OPT_D)) {
    /* "dynamic" option given, but can't be honored */
    return inapplicable_option_error(GretlCmdIndex.fcast, Gretlopt.OPT_D);
    }

    /* NLS: we can only do dynamic out of sample (fc->t1 > pmod->t2) */
    if (pmod.ci == GretlCmdIndex.nls && dyn_ok && (opt & Gretlopt.OPT_D)) {
    if (fc.t1 <= pmod.t2) {
        return inapplicable_option_error(GretlCmdIndex.fcast, Gretlopt.OPT_D);
    }
    }

    if (opt & Gretlopt.OPT_D) {
    /* user requested dynamic forecast and it seems OK */
    fc.method = FcastMethod.dynamic;
    } else if ((dyn_ok || dyn_errs_ok) && fc.t2 > pmod.t2) {
    /* do dynamic f'cast out of sample */
    fc.method = FcastMethod.auto_;
    }

    if (fc.method == FcastMethod.dynamic || fc.method == FcastMethod.auto_) {
    if (fc.t1 > fc.model_t2 + 1) {
        /* we'll probably need intervening forecasts */
        fc.t1 = fc.model_t2 + 1;
    }
    }

    return err;
}

void forecast_init(Forecast * fc)
{
    fc.method = FcastMethod.auto_;

    fc.t1 = 0;
    fc.t2 = 0;
    fc.model_t2 = 0;

    fc.yhat = null;
    fc.sderr = null;
    fc.eps = null;
    fc.dvlags = null;
}

void forecast_free(Forecast * fc)
{
    if (fc.dvlags != null) {
    free(fc.dvlags);
    }

    if (fc.eps != null) {
    free(fc.eps);
    }
}

int fc_add_eps(Forecast * fc, int n)
{
    int t, err = 0;

    fc.eps = cast(typeof(fc.eps)) malloc(n * (*fc.eps).sizeof);

    if (fc.eps == null) {
    err = GretlError.alloc;
    } else {
    for (t=0; t<n; t++) {
        fc.eps[t] = 0.0;
    }
    }

    return err;
}

/* Find the earliest starting point at which we have a previous valid
   observation on the level of series @yno with which to initialize an
   integrated forecast.
*/

int revise_fr_start(FitResid * fr, int yno, const(DatasetRaw) * dset)
{
    const(double) * y = dset.Z[yno];
    int t, t0;

    if (fr.t0 == 0) {
    /* certainly can't start before obs 1 */
    fr.t0 = 1;
    }

    /* previous obs */
    t0 = fr.t0 - 1;

    /* find first non-missing value */
    for (t=t0; t<dset.n; t++) {
    if (!na(y[t])) {
        break;
    }
    }

    t0 = t;

    if (t0 >= dset.n - 2) {
    return GretlError.missdata;
    }

    fr.t0 = t0 + 1;

    if (fr.t1 < fr.t0) {
    fr.t1 = fr.t0;
    }
    if (fr.t2 < fr.t1) {
    fr.t2 = fr.t1;
    }

    fr.fitted[t0] = y[t0];
    fr.resid[t0] = 0.0;

    return 0;
}

/* check whether we really can do an integrated forecast */

int check_integrated_forecast_option(Model * pmod, DatasetRaw * dset, int * pyno)
{
    int err = 0;

    if (pmod.ci != GretlCmdIndex.ols) {
    err = inapplicable_option_error(GretlCmdIndex.fcast, Gretlopt.OPT_I);
    } else {
    int yno = gretl_model_get_depvar(pmod);
    int d, parent;

    d = is_standard_diff(yno, dset, &parent);
    if (!d) {
        err = inapplicable_option_error(GretlCmdIndex.fcast, Gretlopt.OPT_I);
    } else if (pyno != null) {
        /* in caller, make yno = ID of level variable */
        *pyno = parent;
    }
    }

    return err;
}

int check_all_probs_option(Model * pmod)
{
    int err = GretlError.badopt;

    if (pmod.ci == GretlCmdIndex.logit) {
    if (gretl_model_get_int(pmod, "multinom") ||
        gretl_model_get_int(pmod, "ordered")) {
        err = 0;
    }
    } else if (pmod.ci == GretlCmdIndex.probit) {
    if (gretl_model_get_int(pmod, "ordered")) {
        err = 0;
    }
    }

    return err;
}

/* Here we respond to the --exponentiate option for the "fcast" command.
   If the optional int parameter is not given, or given as 0, we apply
   the normal correction if the Doornik-Hansen normality test on the
   residual series fails to reject at the 10 percent level, or otherwise
   Duan's smearing estimator. If the int parameter is given we respond
   as below:

   1 -> apply the normal correction unconditionally
   2 -> apply the smearing correction unconditionally
   3 -> don't apply any correction (just exponentiate)
*/

int add_expon_factor(FitResid * fr, Model * pmod)
{
    double test, skew, xkurt;
    double normpv = NADBL;
    int mode = 0;
    int n, err = 0;

    mode = get_optval_int(GretlCmdIndex.fcast, Gretlopt.OPT_X, &err);
    if (!err && (mode < 0 || mode > 3)) {
    err = GretlError.invarg;
    }

    if (err) {
    return err;
    }

    if (mode == 0) {
    /* automatic: normal or Duan */
    err = series_get_skew_kurt(pmod.t1, pmod.t2, pmod.uhat,
                   &skew, &xkurt, &n);
    if (err) {
        return err;
    }
    test = doornik_chisq(skew, xkurt, n);
    if (na(test)) {
        return GretlError.data;
    } else {
        normpv = chisq_cdf_comp(2, test);
        if (na(normpv)) {
        return GretlError.data;
        } else if (normpv > 0.10) {
        fr.a0meth = 1;
        } else {
        fr.a0meth = 2;
        }
    }
    } else if (mode == 1) {
    /* normal correction forced */
    fr.a0meth = 1;
    } else if (mode == 2) {
    /* Duan forced */
    fr.a0meth = 2;
    } else {
    /* naive exponentiation forced */
    fr.a0meth = 0;
    fr.a0 = 1.0;
    }

    if (fr.a0meth == 1) {
    /* normal correction */
    fr.a0 = exp(pmod.sigma * pmod.sigma / 2.0);
    } else if (fr.a0meth == 2) {
    /* Duan's smearing estimator */
    double sum = 0.0;
    int t;

    for (t=pmod.t1, n=0; t<=pmod.t2; t++) {
        if (!na(pmod.uhat[t])) {
        sum += exp(pmod.uhat[t]);
        n++;
        }
    }
    fr.a0 = sum / n;
    }

    return 0;
}

double log2lev(FitResid * fr, int t)
{
    if (!na(fr.a0)) {
    return fr.a0 * exp(fr.fitted[t]);
    } else {
    return exp(fr.fitted[t]);
    }
}

/* Write to @label (which should be of length SHORTLEN or greater) a
   string to identify the series subject to forecast in @fr, allowing
   for the possibility that the forecast was originally for a series in
   log form and has been exponentiated.
*/

void fr_set_ylabel(FitResid * fr, const(DatasetRaw) * dset)
{
    if (fr.opt & Gretlopt.OPT_X) {
    char[VNAMELEN] parent = 0;
    int v = current_series_index(dset, fr.depvar.ptr);

    if (v > 0 && v < dset.v) {
        series_is_log(dset, v, parent.ptr);
    }
    if (parent[0] != '\0') {
        maybe_trim_varname(fr.ylabel.ptr, parent.ptr);
        return;
    }
    }

    maybe_trim_varname(fr.ylabel.ptr, fr.depvar.ptr);
}

/* Driver for various functions that compute forecasts for different
   sorts of models.
*/

int real_get_fcast(FitResid * fr, Model * pmod, DatasetRaw * dset, Gretlopt opt)
{
    Forecast fc;
    const(double) * y = null;
    int yno = gretl_model_get_depvar(pmod);
    int integrate = (opt & Gretlopt.OPT_I);
    int expon = (opt & Gretlopt.OPT_X);
    int dummy_AR = 0;
    int DM_errs = 0;
    int dyn_errs = 0;
    int asy_errs = 0;
    int int_errs = 0;
    int same_data = 0;
    int nf = 0;
    int t, err;

    if (opt & Gretlopt.OPT_L) {
    err = check_all_probs_option(pmod);
    if (err) {
        return err;
    }
    }

    forecast_init(&fc);

    if (integrate) {
    err = check_integrated_forecast_option(pmod, dset, &yno);
    if (!err) {
        err = revise_fr_start(fr, yno, dset);
    }
    if (err) {
        return err;
    }
    }

    fc.t1 = fr.t1;
    fc.t2 = fr.t2;
    fc.model_t2 = pmod.t2;

    err = get_forecast_method(&fc, pmod, dset, opt);
    if (err) {
    return err;
    }

    if (pmod.ci == GretlCmdIndex.nls && fc.method == FcastMethod.static_) {
    asy_errs = 1;
    }

    if (pmod.ci == GretlCmdIndex.panel || pmod.ci == GretlCmdIndex.dpanel) {
    ; /* don't do the things below */
    } else if (!FCAST_SPECIAL(pmod.ci) && !integrate) {
    if (!AR_MODEL(cast(GretlCmdIndex) pmod.ci) && fc.dvlags == null) {
        /* we'll do Davidson-MacKinnon error variance */
        DM_errs = 1;
    } else if (fc.method == FcastMethod.dynamic) {
        /* we'll do dynamic forecast errors throughout */
        dyn_errs = 1;
    } else if (fc.method == FcastMethod.auto_ && fc.t2 > pmod.t2) {
        /* do dynamic forecast errors out of sample */
        dyn_errs = 1;
    }
    }

    if (integrate && pmod.ci == GretlCmdIndex.ols &&
    fc.dvlags == null && fc.method != FcastMethod.static_) {
    int_errs = 1;
    }

    if (dyn_errs && !AR_MODEL(cast(GretlCmdIndex) pmod.ci) && fc.dvlags != null) {
    /* create dummy AR info structure for model with lagged
       dependent variable */
    int derr = dummy_ar_info_init(pmod);

    if (derr) {
        dyn_errs = 0;
    } else {
        dummy_AR = 1;
    }
    }

    if (DM_errs || dyn_errs || asy_errs || int_errs) {
    err = fit_resid_add_sderr(fr);
    }

    if (err) {
    return err;
    }

    fc.yhat = fr.fitted;
    fc.sderr = fr.sderr;

    if (pmod.ci == GretlCmdIndex.arma && fc.method == FcastMethod.static_) {
    fc_add_eps(&fc, dset.n);
    }

    /* compute the actual forecast */
    if (pmod.ci == GretlCmdIndex.dpanel || individual_effects_model(pmod)) {
    err = panel_fcast(&fc, pmod, dset, opt);
    } else if (DM_errs) {
    err = static_fcast_with_errs(&fc, pmod, dset, opt);
    } else if (pmod.ci == GretlCmdIndex.nls) {
    err = nls_fcast(&fc, pmod, dset);
    } else if (pmod.ci == GretlCmdIndex.midasreg) {
    err = midas_fcast(&fc, pmod, dset);
    } else if (SIMPLE_AR_MODEL(cast(GretlCmdIndex) pmod.ci) || dummy_AR) {
    err = ar_fcast(&fc, pmod, dset);
    } else if (pmod.ci == GretlCmdIndex.arma) {
    err = arma_fcast(&fc, pmod, dset);
    } else if (pmod.ci == GretlCmdIndex.garch) {
    err = garch_fcast(&fc, pmod, dset);
    } else if (integrate) {
    err = integrated_fcast(&fc, pmod, yno, dset);
    } else if (pmod.ci == GretlCmdIndex.logit && gretl_model_get_int(pmod, "multinom")) {
    if (opt & Gretlopt.OPT_L) {
        /* --all-probs */
        err = mlogit_allprobs(&fc, pmod, dset);
    } else {
        err = mlogit_fcast(&fc, pmod, dset);
    }
    } else {
    err = linear_fcast(&fc, pmod, yno, dset, null, opt);
    }

    /* free any auxiliary info */
    forecast_free(&fc);

    if (opt & Gretlopt.OPT_L) {
        /* --all-probs: we're done */
    return err;
    }

    if (dummy_AR) {
    free(pmod.arinfo.arlist);
    free(pmod.arinfo);
    pmod.arinfo = null;
    }

    if (expon) {
    add_expon_factor(fr, pmod);
    }

    same_data = same_dataset(pmod, dset);
    y = dset.Z[yno];

    for (t=0; t<fr.nobs; t++) {
    if (t >= fr.t1 && t <= fr.t2) {
        /* in the forecast range */
        if (!na(fr.fitted[t])) {
        nf++;
        if (expon) {
            fr.fitted[t] = log2lev(fr, t);
        }
        }
    } else if (same_data && t >= fr.t0 && t <= fr.t2 &&
           t >= pmod.t1 && t <= pmod.t2) {
        /* within sample observations */
        if (integrate) {
        fr.fitted[t] = fr.fitted[t-1] + pmod.yhat[t];
        fr.resid[t] = fr.resid[t-1] + pmod.uhat[t];
        } else if (expon) {
        if (!na(fr.a0)) {
            fr.fitted[t] = fr.a0 * exp(pmod.yhat[t]);
        } else {
            fr.fitted[t] = exp(pmod.yhat[t]);
        }
        fr.resid[t] = exp(y[t]) - fr.fitted[t];
        } else {
        fr.fitted[t] = pmod.yhat[t];
        fr.resid[t] = pmod.uhat[t];
        }
    }
    fr.actual[t] = expon ? exp(y[t]) : y[t];
    }

    if (nf == 0) {
    err = GretlError.missdata;
    } else {
    fit_resid_set_dec_places(fr);
    strcpy(fr.depvar.ptr, dset.varname[yno]);
    fr_set_ylabel(fr, dset);
    fr.df = pmod.dfd;
    }

    return err;
}

int matrix_forecast_supported(Model * pmod)
{
    if (pmod.errcode) {
        return 0;
    } else if (pmod.ncoeff < 1) {
        return 0;
    } else if (gretl_model_get_int(pmod, "ldepvar")) {
        return 0;
    } else if (pmod.ci != GretlCmdIndex.ols &&
               pmod.ci != GretlCmdIndex.logit &&
               pmod.ci != GretlCmdIndex.probit) {
        /* FIXME relax this! */
        return 0;
    } else if (pmod.xlist == null) {
        return 0;
    }

    if (pmod.ncoeff != pmod.xlist[0]) {
        return 0;
    } else {
        fprintf(stderr, "HERE ncoeff = nx = %d, OK?\n".ptr, pmod.ncoeff);
    }

    return 1;
}

/* April 2025: simple forecast interface: X-data given in matrix
   form and prediction returned in matrix form.
*/

GretlMatrixRaw * matrix_forecast(Model * pmod, const(GretlMatrixRaw) * X, int * err)
{
    GretlMatrixRaw * F = null;
    FitResid * fr = null;
    Forecast fc;
    int n = 0;

    if (!matrix_forecast_supported(pmod)) {
        *err = GretlError.data;
    } else if (X == null || X.rows < 1) {
        *err = GretlError.data;
    } else {
        n = X.rows;
        fr = fit_resid_new_for_model(pmod, null, 0, n, 0, Gretlopt.OPT_NONE, err);
    }

    if (!*err) {
        int yno = pmod.list[1]; /* FIXME */

        forecast_init(&fc);
        fc.method = FcastMethod.static_;
        fc.yhat = fr.fitted;
        fc.sderr = fr.sderr;
        fc.t1 = 0;
        fc.t2 = n - 1;
        *err = linear_fcast(&fc, pmod, yno, null, X, Gretlopt.OPT_NONE);
    }

    if (!*err) {
        int cols = fc.sderr == null ? 1 : 2;
        int i;

        F = gretl_zero_matrix_new(n, cols);
        if (F == null) {
            *err = GretlError.alloc;
        } else {
            for (i=0; i<n; i++) {
                gretl_matrix_set(F, i, 0, fc.yhat[i]);
                if (cols > 1) {
                    gretl_matrix_set(F, i, 1, fc.sderr[i]);
                }
            }
        }
    }

    free_fit_resid(fr);

    return F;
}

int fcast_get_limit(const(char) * s, DatasetRaw * dset)
{
    double x;
    int t, err = 0;

    if (gretl_is_scalar(s)) {
    x = gretl_scalar_get_value(s, null);
    } else {
    x = generate_scalar(s, dset, &err);
    }

    if (x < 1 || x > dset.n) {
    gretl_error_clear();
    gretl_errmsg_set(_("Observation number out of bounds").ptr);
    t = -1;
    } else {
    t = cast(int) (x - 1);
    }

    return t;
}

enum OsCode {
    ok = 0,
    err,
    panel
}

int out_of_sample_check(Model * pmod, DatasetRaw * dset)
{
    int panel, ret;

    if (pmod.ci == GretlCmdIndex.panel || pmod.ci == GretlCmdIndex.dpanel ||
    (pmod.ci == GretlCmdIndex.ols && dataset_is_panel(dset))) {
    panel = 1;
    ret = OsCode.err;
    if (gretl_model_get_int(pmod, "ntdum") > 0) {
        gretl_errmsg_set(_("Specification includes time dummies: cannot forecast out of sample").ptr);
        return ret;
    }
    } else {
    panel = 0;
    ret = OsCode.ok;
    }

    if (pmod.smpl.t1 == 0 && pmod.smpl.t2 == dset.n - 1) {
    /* no out-of-sample obs reachable via t1, t2 */
    ret = OsCode.err;
    if (!panel) {
        if (pmod.t2 < dset.n - 1) {
        /* 2019-11-15: may be OK for dynamic model */
        ret = OsCode.ok;
        }
    } else {
        /* the panel case */
        DatasetRaw * fullset = fetch_full_dataset();

        if (dataset_is_panel(fullset)) {
        int Tfull = fullset.pd;
        int Tcurr = dset.pd;
        int Nfull = fullset.n / Tfull;
        int Ncurr = dset.n / Tcurr;

        if (Ncurr == Nfull && Tcurr < Tfull &&
            dset.v == fullset.v) {
            /* sub-sampled in the time dimension only, and
               no series added or deleted (we trust)
            */
            ret = OsCode.panel;
        }
        }
    }
    if (ret == OsCode.err) {
        gretl_errmsg_set(_("No out-of-sample observations are available").ptr);
    }
    } else if (panel) {
    gretl_errmsg_set(_("The --out-of-sample option is only supported in the time dimension").ptr);
    }

    return ret;
}

int parse_forecast_string(const(char) * s, Gretlopt opt, Model * pmod, int t2est, DatasetRaw * dset, int * pt1, int * pt2, int * pk, char * vname, int * os_case)
{
    char[32][4] f = 0;
    char * t1str = null, t2str = null;
    char * kstr = null, vstr = null;
    int nmax, nmin = (opt & Gretlopt.OPT_R)? 1 : 0;
    int t1 = 0, t2 = 0;
    int nf = 0;
    int err = 0;

    /* "static" and "recursive" can't be combined */
    err = incompatible_options(opt, Gretlopt.OPT_S | Gretlopt.OPT_R);
    if (err) {
    return err;
    }

    vname[0] = '\0';

    /* How many fields should we be looking for in the user input?
       If OPT_R ("recursive") is given, the max is 4:

       t1, t2, k (steps-ahead), vname

       Otherwise the max is 3 (no k-value):

       t1, t2, vname

       However, if OPT_O ("out-of-sample") is given, we should not
       expect t1 and t2, so the max becomes 2 (if OPT_R), else 1.

       Also note: t1 and t2 need not be given (even in the absence
       of OPT_O) since these default to the current sample range.
    */

    if (opt & Gretlopt.OPT_O) {
    nmax = (opt & Gretlopt.OPT_R)? 2 : 1;
    } else {
    nmax = (opt & Gretlopt.OPT_R)? 4 : 3;
    }

    if (s[0] == '\0') {
    nf = 0;
    } else {
    nf = sscanf(s, "%31s %31s %31s %31s".ptr, f[0].ptr, f[1].ptr, f[2].ptr, f[3].ptr);
    if (nf > nmax) {
        /* try for parenthesized t1, t2 terms? */
        nf = sscanf(s, "(%31[^)]) (%31[^)]) %31s %31s".ptr, f[0].ptr, f[1].ptr, f[2].ptr, f[3].ptr);
    }
    }

    if (nf < nmin || nf > nmax) {
    gretl_errmsg_sprintf("fcast: expected %d to %d fields in input, got %d".ptr,
                 nmin, nmax, nf);
    return GretlError.parse;
    }

    if (opt & Gretlopt.OPT_R) {
    if (nf == 4) {
        /* t1, t2, k, vname */
        t1str = f[0].ptr;
        t2str = f[1].ptr;
        kstr = f[2].ptr;
        vstr = f[3].ptr;
    } else if (nf == 3) {
        /* t1, t2, k */
        t1str = f[0].ptr;
        t2str = f[1].ptr;
        kstr = f[2].ptr;
    } else if (nf == 2) {
        /* k, vname */
        kstr = f[0].ptr;
        vstr = f[1].ptr;
    } else if (nf == 1) {
        /* just k */
        kstr = f[0].ptr;
    }
    } else {
    if (nf == 3) {
        /* t1, t2, vname */
        t1str = f[0].ptr;
        t2str = f[1].ptr;
        vstr = f[2].ptr;
    } else if (nf == 2) {
        /* t1, t2 */
        t1str = f[0].ptr;
        t2str = f[1].ptr;
    } else if (nf == 1) {
        /* just vname */
        vstr = f[0].ptr;
    }
    }

    if (t1str != null || t2str != null) {
    if (opt & Gretlopt.OPT_O) {
        /* t1, t2 should not be given with --out-of-sample */
        gretl_errmsg_set("fcast: unexpected t1 and/or t2 field in input");
        return GretlError.data;
    }
    }

    if (kstr != null && pk != null) {
    *pk = gretl_int_from_string(kstr, &err);
    if (err) {
        return err;
    }
    }

    if (vstr != null) {
    strncat(vname, vstr, VNAMELEN - 1);
    }

    if (t1str != null && t2str != null) {
    t1 = dateton(t1str, dset);
    if (t1 < 0) {
        t1 = fcast_get_limit(t1str, dset);
    }
    t2 = dateton(t2str, dset);
    if (t2 < 0) {
        t2 = fcast_get_limit(t2str, dset);
    }
    if (t1 < 0 || t2 < 0 || t2 < t1) {
        err = GretlError.data;
    }
    } else if (opt & Gretlopt.OPT_O) {
    /* out of sample, if possible */
    if (pmod != null) {
        *os_case = out_of_sample_check(pmod, dset);
    }
    if (*os_case == OsCode.err) {
        err = GretlError.data;
    } else if (*os_case == OsCode.ok) {
        if (dset.n - t2est - 1 > 0) {
        t1 = t2est + 1;
        t2 = dset.n - 1;
        } else {
        err = GretlError.obs;
        }
    }
    } else {
    /* default: current sample range */
    t1 = dset.t1;
    t2 = dset.t2;
    }

    if (!err) {
    /* in case we hit any "temporary" errors above */
    gretl_error_clear();
    *pt1 = t1;
    *pt2 = t2;
    }

    return err;
}

/**
 * get_forecast:
 * @pmod: the model from which forecasts are wanted.
 * @t1: start of forecast range.
 * @t2: end of forecast range.
 * @pre_n: number of pre-forecast observations to include.
 * @dset: dataset struct.
 * @opt: if %OPT_D, force a dynamic forecast; if %OPT_S, force
 * a static forecast.  By default, the forecast is static within
 * the data range over which the model was estimated, and dynamic
 * out of sample (in cases where a dynamic forecast is meaningful).
 * If @opt includes %OPT_I, integrate the forecast (only relevant
 * if the dependent variable in the model in question is recognized
 * as the first difference of another variable).
 * @err: location to receive error code.
 *
 * Allocates a #FITRESID structure and fills it out with forecasts
 * based on @pmod, over the specified range of observations.
 * For some sorts of models forecast standard errors are also
 * computed (these appear in the %sderr member of the structure
 * to which a pointer is returned; otherwise the %sderr member is
 * %NULL).
 *
 * The calculation of forecast errors, where applicable, is based
 * on Davidson and MacKinnon, Econometric Theory and Methods,
 * chapter 3 (p. 104), which shows how the variance of forecast errors
 * can be computed given the covariance matrix of the parameter
 * estimates, provided the error term may be assumed to be serially
 * uncorrelated.
 *
 * Returns: pointer to allocated structure, or %NULL on failure,
 * in which case an error code is assigned via @err.
 */

FitResid * get_forecast(Model * pmod, int t1, int t2, int pre_n, DatasetRaw * dset, Gretlopt opt, int * err)
{
    FitResid * fr;

    if (pmod.errcode) {
    *err = GretlError.data;
    return null;
    }

    fr = fit_resid_new_for_model(pmod, dset, t1, t2, pre_n, opt, err);

    if (!*err) {
    *err = real_get_fcast(fr, pmod, dset, opt);
    if (*err) {
        free_fit_resid(fr);
        fr = null;
    }
    }

    return fr;
}

void forecast_matrix_cleanup()
{
    set_fcast_matrices(null, null);
}

GretlMatrixRaw * get_forecast_matrix(int idx, int * err)
{
    GretlMatrixRaw * ret = null;
    GretlMatrixRaw * M = null;

    if (idx == GretlCmdIndex.fcast) {
    M = fcast_matrix;
    } else if (idx == ModelDataIndex.fcse) {
    M = fcerr_matrix;
    }

    if (M == null) {
    *err = GretlError.badstat;
    } else {
    ret = gretl_matrix_copy(M);
    if (ret == null) {
        *err = GretlError.alloc;
    }
    }

    return ret;
}

void fr_adjust_sample(FitResid * fr, int * ft1, int * ft2, int * et1, int * et2)
{
    int t;

    t = fr.t0;
    while (t<=fr.t2 && na(fr.fitted[t])) {
    *ft1 = ++t;
    }

    t = fr.t2;
    while (t>=*ft1 && na(fr.fitted[t])) {
    *ft2 = --t;
    }

    if (fr.sderr != null) {
    t = fr.t0;
    while (t<=fr.t2 && na(fr.sderr[t])) {
        *et1 = ++t;
    }

    t = fr.t2;
    while (t>=*et1 && na(fr.sderr[t])) {
        *et2 = --t;
    }
    }
}

int set_forecast_matrices_from_fr(FitResid * fr)
{
    GretlMatrixRaw * f = null;
    GretlMatrixRaw * e = null;
    int ft1 = fr.t0;
    int ft2 = fr.t2;
    int et1 = fr.t0;
    int et2 = fr.t2;
    int T, t, s;
    int err = 0;

    fr_adjust_sample(fr, &ft1, &ft2, &et1, &et2);

    T = ft2 - ft1 + 1;
    if (T <= 0) {
    return 0;
    }

    f = gretl_matrix_alloc(T, 1);
    if (f == null) {
    return GretlError.alloc;
    }

    gretl_matrix_set_t1(f, ft1);
    gretl_matrix_set_t2(f, ft2);

    if (fr.sderr != null) {
    T = et2 - et1 + 1;
    if (T > 0) {
        e = gretl_matrix_alloc(T, 1);
        if (e == null) {
        err = GretlError.alloc;
        } else {
        gretl_matrix_set_t1(e, et1);
        gretl_matrix_set_t2(e, et2);
        }
    }
    }

    s = 0;
    for (t=ft1; t<=ft2; t++) {
    f.val[s++] = fr.fitted[t];
    }

    if (e != null) {
    s = 0;
    for (t=et1; t<=et2; t++) {
        e.val[s++] = fr.sderr[t];
    }
    }

    set_fcast_matrices(f, e);

    return err;
}

int matrix_first_ok_row(const(GretlMatrixRaw) * a, int jmin, int cols)
{
    int jmax = jmin + cols;
    double x;
    int i, j, miss;

    for (i=0; i<a.rows; i++) {
    miss = 0;
    for (j=jmin; j<jmax && !miss; j++) {
        x = gretl_matrix_get(a, i, j);
        if (na(x)) {
        miss = 1;
        }
    }
    if (!miss) {
        break;
    }
    }

    return i;
}

int matrix_last_ok_row(const(GretlMatrixRaw) * a, int jmin, int cols)
{
    int jmax = jmin + cols;
    double x;
    int i, j, miss;

    for (i=a.rows-1; i>=0; i--) {
    miss = 0;
    for (j=jmin; j<jmax && !miss; j++) {
        x = gretl_matrix_get(a, i, j);
        if (na(x)) {
        miss = 1;
        }
    }
    if (!miss) {
        break;
    }
    }

    return i;
}

void F_matrix_adjust_sample(const(GretlMatrixRaw) * F, int * f0, int * fn, int * e0, int * en)
{
    int k = F.cols / 2;

    *f0 = matrix_first_ok_row(F, 0, k);
    *fn = matrix_last_ok_row(F, 0, k);

    *e0 = matrix_first_ok_row(F, k, k);
    *en = matrix_last_ok_row(F, k, k);
}

int set_forecast_matrices_from_F(const(GretlMatrixRaw) * F, int imin, int imax)
{
    GretlMatrixRaw * F1 = null;
    GretlMatrixRaw * f = null;
    GretlMatrixRaw * e = null;
    int n = F.cols;
    int k = n / 2;
    int fT = F.rows;
    int eT = F.rows;
    int f0 = 0, fn = fT;
    int e0 = 0, en = eT;
    double x;
    int i, j, Ft1, mt1;
    int err = 0;

    Ft1 = gretl_matrix_get_t1(F);

    if (imin == imax) {
    /* extract one pair of columns */
    F1 = gretl_matrix_alloc(fT, 2);
    if (F1 == null) {
        return GretlError.alloc;
    }
    for (i=0; i<fT; i++) {
        x = gretl_matrix_get(F, i, imin);
        gretl_matrix_set(F1, i, 0, x);
        x = gretl_matrix_get(F, i, imin + k);
        gretl_matrix_set(F1, i, 1, x);
    }
    F = F1;
    k = 1;
    }

    F_matrix_adjust_sample(F, &f0, &fn, &e0, &en);

    fT = fn - f0 + 1;
    if (fT <= 0) {
    return GretlError.missdata;
    }

    f = gretl_matrix_alloc(fT, k);
    if (f == null) {
    return GretlError.alloc;
    }

    mt1 = Ft1 + f0;
    gretl_matrix_set_t1(f, mt1);
    gretl_matrix_set_t2(f, mt1 + fT - 1);

    eT = en - e0 + 1;

    if (eT > 0) {
    e = gretl_matrix_alloc(eT, k);
    if (e == null) {
        err = GretlError.alloc;
    } else {
        mt1 = Ft1 + e0;
        gretl_matrix_set_t1(e, mt1);
        gretl_matrix_set_t2(e, mt1 + eT - 1);
    }
    }

    for (j=0; j<k; j++) {
    for (i=0; i<fT; i++) {
        x = gretl_matrix_get(F, i + f0, j);
        gretl_matrix_set(f, i, j, x);
    }
    }

    if (e != null) {
    for (j=0; j<k; j++) {
        for (i=0; i<eT; i++) {
        x = gretl_matrix_get(F, i + e0, k + j);
        gretl_matrix_set(e, i, j, x);
        }
    }
    }

    set_fcast_matrices(f, e);

    if (F1 != null) {
    gretl_matrix_free(F1);
    }

    return err;
}

/* compatibility with behaviour of old fcast command */

int add_fcast_to_dataset(FitResid * fr, const(char) * vname, DatasetRaw * dset, Prn * prn)
{
    int oldv = dset.v;
    int v, err = 0;

    v = series_index(dset, vname);

    if (v == dset.v) {
    /* new variable */
    if (check_varname(vname)) {
        err = GretlError.data;
    } else {
        err = dataset_add_series(dset, 1);
        if (!err) {
        strcpy(dset.varname[v], vname);
        }
    }
    }

    if (!err) {
    int t;

    for (t=0; t<dset.n; t++) {
        dset.Z[v][t] = fr.fitted[t];
    }

    series_set_label(dset, v, _("predicted values").ptr);

    if (gretl_messages_on()) {
        if (v < oldv) {
        pprintf(prn, _("Replaced series %s (ID %d)").ptr,
            vname, v);
        } else {
        pprintf(prn, _("Generated series %s (ID %d)").ptr,
            vname, v);
        }
        pputc(prn, '\n');
    }
    }

    return err;
}

void print_probs_matrix(Prn * prn)
{
    if (fcast_matrix != null) {
        if (gretl_echo_on()) {
            pputc(prn, '\n');
        }
        pprintf(prn, "%s\n\n".ptr, _("Estimated outcome probabilities").ptr);
        gretl_matrix_print_to_prn(fcast_matrix, null, prn);
    }
}

int add_fcast_matrix(const(char) * name, DatasetRaw * dset, Prn * prn)
{
    int err = 0;

    if (gretl_reserved_word(name)) {
        err = GretlError.types;
    } else if (current_series_index(dset, name) >= 0) {
    err = GretlError.types;
    } else {
    err = user_var_add_or_replace(name, GretlType.matrix,
                      fcast_matrix);
    fcast_matrix = null;
    }

    return err;
}

int model_do_forecast(const(char) * str, Model * pmod, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    char[VNAMELEN] vname;
    FitResid * fr = null;
    int t1, t2, k = -1;
    int all_probs = (opt & Gretlopt.OPT_L);
    int os_case = 0;
    int quiet = 0;
    int err;

    /* FIXME relax some of these limitations? */
    if (pmod.ci == GretlCmdIndex.heckit || pmod.ci == GretlCmdIndex.duration) {
    return GretlError.notimp;
    } else if (pmod.ci == GretlCmdIndex.panel && !(pmod.opt & Gretlopt.OPT_P)) {
    if (opt & Gretlopt.OPT_R) {
        return GretlError.notimp;
    }
    } else if (pmod.ci == GretlCmdIndex.probit && (pmod.opt & Gretlopt.OPT_E)) {
    return GretlError.notimp;
    }

    /* can't use --plot with --all-probs */
    if (all_probs && (opt & Gretlopt.OPT_U)) {
        return GretlError.badopt;
    }

    /* OPT_I for integrate: reject for non-OLS, or if the dependent
       variable is not recognized as a first difference */
    if (opt & Gretlopt.OPT_I) {
    err = check_integrated_forecast_option(pmod, dset, null);
    if (err) {
        return err;
    }
    }

    err = parse_forecast_string(str, opt, pmod, pmod.t2, dset,
                &t1, &t2, &k, vname.ptr, &os_case);
    if (err) {
    return err;
    }

    if (vname[0] != '\0') {
    /* saving to named series or matrix: be quiet */
    opt |= (Gretlopt.OPT_A | Gretlopt.OPT_Q);
    }

    quiet = (opt & Gretlopt.OPT_Q);

    if (os_case == OsCode.panel) {
    return panel_os_special(pmod, dset, vname.ptr, opt, prn);
    } else if (opt & Gretlopt.OPT_R) {
    fr = recursive_OLS_k_step_fcast(pmod, dset, t1, t2,
                    k, 0, &err);
    } else {
    fr = get_forecast(pmod, t1, t2, 0, dset, opt, &err);
    }

    if (err) {
        goto bailout;
    }

    if (!err && (opt & Gretlopt.OPT_A)) {
        /* named output */
    if (opt & Gretlopt.OPT_L) {
            /* add forecast matrix directly */
        err = add_fcast_matrix(vname.ptr, dset, prn);
    } else {
            /* add forecast series directly */
        err = add_fcast_to_dataset(fr, vname.ptr, dset, prn);
    }
    }

    if (!quiet && (opt & Gretlopt.OPT_L)) {
        print_probs_matrix(prn);
    } else if (!quiet) {
        Gretlopt printopt = opt;

        if (opt & Gretlopt.OPT_U) {
            /* do graph (from command line) */
            printopt |= Gretlopt.OPT_P;
        }
        err = text_print_forecast(fr, dset, printopt, prn);
    }

    if (!err && !(opt & Gretlopt.OPT_L)) {
    set_forecast_matrices_from_fr(fr);
    }

 bailout:

    free_fit_resid(fr);

    return err;
}

int get_sys_fcast_var(const(int) * ylist, const(char) * vname, DatasetRaw * dset)
{
    int i, vi;

    for (i=0; i<ylist[0]; i++) {
    vi = ylist[i+1];
    if (!strcmp(vname, dset.varname[vi])) {
        return i;
    }
    }

    return -1;
}

/* grab the relevant range from system forecast matrix and write it
   into a FITRESID struct for printing or graphing
*/

int fill_system_forecast(FitResid * fr, int i, int yno, GretlVar * var, EquationSystem * sys, const(GretlMatrixRaw) * F, DatasetRaw * dset)
{
    int m = F.cols / 2;
    int s, t, nf;
    int err = 0;

    strcpy(fr.depvar.ptr, dset.varname[yno]);

    /* "pre-forecast" observations */
    for (t=fr.t0; t<fr.t1; t++) {
    fr.actual[t] = dset.Z[yno][t];
    if (sys != null) {
        if (i < sys.neqns && sys.lists[i][1] == yno) {
        /* Note: right now we can only handle endogenous
           variables that appear on the LHS of stochastic
           equations, because only such variables have an
           associated column in sys->yhat.
        */
        if (t >= sys.t1 && t <= sys.t2) {
            s = t - sys.t1;
            fr.fitted[t] = gretl_matrix_get(sys.yhat, s, i);
        }
        }
    } else if (var != null) {
        if (t >= var.t1 && t <= var.t2) {
        s = t - var.t1;
        fr.fitted[t] = fr.actual[t] - gretl_matrix_get(var.E, s, i);
        }
    }
    }

    /* actual forecasts */
    nf = 0;
    for (t=fr.t1, s=0; t<=fr.t2; t++, s++) {
    fr.actual[t] = dset.Z[yno][t];
    fr.fitted[t] = gretl_matrix_get(F, s, i);
    if (!na(fr.fitted[t])) {
        nf++;
    }
    if (fr.sderr != null) {
        fr.sderr[t] = gretl_matrix_get(F, s, i + m);
    }
    }

    if (nf == 0) {
    err = GretlError.missdata;
    } else {
    fit_resid_set_dec_places(fr);
    }

    return err;
}

int system_do_forecast(const(char) * str, void * ptr, int type, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    char[VNAMELEN] vname;
    GretlVar * var = null;
    EquationSystem * sys = null;
    const(int) * ylist = null;
    const(GretlMatrixRaw) * F = null;
    int t1 = 0, t2 = 0;
    int t2est, ci;
    int imax, imin = 0;
    int os_case = 0;
    int df = 0;
    int err = 0;

    if (type == GretlObjType.var) {
    var = cast(GretlVar *) ptr;
    imax = var.neqns - 1;
    t2est = var.t2;
    ci = var.ci;
    df = var.df;
    ylist = var.ylist;
    } else {
    sys = cast(EquationSystem *) ptr;
    imax = sys.neqns + sys.nidents - 1;
    t2est = sys.t2;
    ci = GretlCmdIndex.system;
    df = sys.df;
    ylist = sys.ylist;
    }

    err = parse_forecast_string(str, opt, null, t2est, dset,
                &t1, &t2, null, vname.ptr, &os_case);

    if (!err) {
    if (var != null) {
        F = gretl_VAR_get_forecast_matrix(var, t1, t2, dset,
                          opt, &err);
    } else {
        F = system_get_forecast_matrix(sys, t1, t2, dset,
                       opt, &err);
    }
    }

    if (!err && vname[0] != '\0') {
    imin = get_sys_fcast_var(ylist, vname.ptr, dset);
    if (imin < 0) {
        err = GretlError.data;
    } else {
        imax = imin;
    }
    }

    if (!err) {
    /* arrange to save forecast and errors */
    err = set_forecast_matrices_from_F(F, imin, imax);
    }

    if (err) {
    return err;
    }

    if (!(opt & Gretlopt.OPT_Q)) {
    /* assemble and print per-equation forecasts */
    Gretlopt printopt = (opt & Gretlopt.OPT_N)? Gretlopt.OPT_N : Gretlopt.OPT_NONE;
    FitResid * fr = null;
    int i, asy;

    asy = (ci == GretlCmdIndex.vecm);

    fr = fit_resid_new_for_system(asy, dset, t1, t2, 0, &err);
    if (err) {
        return err;
    }

        if (asy) {
        /* asymptotic normal */
        fr.df = var.T;
    } else {
        fr.df = df;
    }

    if ((opt & Gretlopt.OPT_U) && vname[0] != '\0') {
        /* produce a plot */
        printopt |= Gretlopt.OPT_P;
    }

    for (i=imin; i<=imax && !err; i++) {
        err = fill_system_forecast(fr, i, ylist[i+1], var, sys,
                       F, dset);
        if (!err) {
        err = text_print_forecast(fr, dset, printopt, prn);
        }
        printopt |= Gretlopt.OPT_Q;
    }

    free_fit_resid(fr);
    }

    return err;
}

/**
 * do_forecast:
 * @str: command string, which may include a starting
 * observation and ending observation, and/or the name of a
 * variable for saving the forecast values.
 * @dset: dataset struct.
 * @opt: if %OPT_D, force a dynamic forecast; if %OPT_S, force
 * a static forecast.  By default, the forecast is static within
 * the data range over which the model was estimated, and dynamic
 * out of sample (in cases where this distinction is meaningful).
 * %OPT_R: do recursive forecast.
 * %OPT_Q: suppress printing of the forecast;
 * %OPT_P: ensure that the values are printed.
 * %OPT_U: produce gnuplot plot.
 * @prn: gretl printing struct.
 *
 * In the case of "simple" models with an autoregressive error term
 * (%AR, %AR1) the forecast values incorporate the predictable
 * component of the error.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int do_forecast(const(char) * str, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    void * ptr;
    GretlObjType type;
    int err;

    ptr = get_last_model(&type);
    if (ptr == null) {
    return GretlError.badstat;
    }

    if ((opt & Gretlopt.OPT_U) && dataset_is_panel(dset)) {
    gretl_errmsg_set(_("Forecast plot not implemented for panel data").ptr);
    err = GretlError.notimp;
    } else if ((opt & (Gretlopt.OPT_R | Gretlopt.OPT_I)) && type != GretlObjType.eqn) {
    /* "recursive", "integrate": single equations only */
    err = GretlError.badopt;
    } else if (type == GretlObjType.eqn) {
    err = model_do_forecast(str, cast(Model *) ptr, dset, opt, prn);
    } else if (type == GretlObjType.sys || type == GretlObjType.var) {
    err = system_do_forecast(str, ptr, type, dset, opt, prn);
    } else {
    err = GretlError.data;
    }

    return err;
}

/* try to determine in advance how far we can go with a forecast,
   either dynamic or static (@ftype) */

int fcast_get_t2max(const(int) * dvlags, const(Model) * pmod, const(DatasetRaw) * dset, int ftype)
{
    const(double) * ay = null;
    int i, vi, t, l0;

    if (pmod.ci == GretlCmdIndex.arma && ftype == FcastMethod.static_) {
    int yno = gretl_model_get_depvar(pmod);

    ay = dset.Z[yno];
    }

    l0 = pmod.xlist == null ? 0 : pmod.xlist[0];

    for (t=pmod.t2; t<dset.n; t++) {
    int p, vj, all_ok = 1;

    if (ay != null && na(ay[t-1])) {
        /* FIXME? */
        break;
    }

    for (i=1; i<=l0; i++) {
        vi = pmod.xlist[i];
        if (vi == 0) {
        continue;
        } else if (dvlags != null && dvlags[i-1] != 0) {
        continue;
        } else if (is_trend_variable(dset.Z[vi], dset.n)) {
        continue;
        } else if (is_periodic_dummy(dset.Z[vi], dset)) {
        continue;
        } else if ((p = series_get_lag(dset, vi)) > 0) {
        vj = series_get_parent_id(dset, vi);
        if (vj > 0 && !na(dset.Z[vj][t-p])) {
            continue;
        }
        }
        if (na(dset.Z[vi][t])) {
        all_ok = 0;
        break;
        }
    }

    if (!all_ok) {
        t--;
        break;
    } else if (t == dset.n - 1) {
        break;
    }
    }

    return t;
}

/**
 * get_system_forecast:
 * @p: pointer to the VAR or equation system from which
 * forecasts are wanted.
 * @ci: command index for system (%VAR, %VECM or %SYSTEM)
 * @i: 0-based index for the variable to forecast, within
 * the equation system.
 * @t1: start of forecast range.
 * @t2: end of forecast range.
 * @pre_n: number of pre-forecast observations to include.
 * @dset: dataset struct.
 * @opt: if %OPT_D, force a dynamic forecast; if %OPT_S, force
 * a static forecast.  By default, the forecast is static within
 * the data range over which the model was estimated, and dynamic
 * out of sample.
 * @err: location to receive error code.
 *
 * Allocates a #FITRESID structure and fills it out with forecasts
 * based on the system at location @p, over the specified range of
 * observations.
 *
 * Returns: pointer to allocated structure, or %NULL on failure.
 */

FitResid * get_system_forecast(void * p, int ci, int i, int t1, int t2, int pre_n, DatasetRaw * dset, Gretlopt opt, int * err)
{
    FitResid * fr;
    GretlVar * var = null;
    EquationSystem * sys = null;
    const(GretlMatrixRaw) * F = null;
    int nf = t2 - t1 + 1;
    int asy, df = 0, yno = 0;

    if (nf <= 0) {
    *err = GretlError.data;
    return null;
    }

    if (ci == GretlCmdIndex.var || ci == GretlCmdIndex.vecm) {
    var = cast(GretlVar *) p;
    yno = var.ylist[i+1];
    df = var.df;
    F = gretl_VAR_get_forecast_matrix(var, t1, t2, dset, opt, err);
    } else if (ci == GretlCmdIndex.system) {
    sys = cast(EquationSystem *) p;
    yno = sys.ylist[i+1];
    df = sys.df;
    F = system_get_forecast_matrix(sys, t1, t2, dset, opt, err);
    } else {
    *err = GretlError.data;
    }

    if (*err) {
    fprintf(stderr, "get_system_forecast: matrix F is NULL\n".ptr);
    return null;
    }

    asy = (ci == GretlCmdIndex.vecm);

    fr = fit_resid_new_for_system(asy, dset, t1, t2, pre_n, err);
    if (*err) {
    return null;
    }

    if (asy) {
    /* asymptotic normal */
    fr.df = var.T;
    } else {
    fr.df = df;
    }

    *err = fill_system_forecast(fr, i, yno, var, sys, F, dset);

    if (*err) {
    free_fit_resid(fr);
    fr = null;
    }

    return fr;
}

/* Try to determine whether adding observations to the dataset (without
   actually adding more "real" data) can serve to extend the range of
   out-of-sample prediction. The answer will be No if the specification
   includes exogenous regressors other than deterministic terms that can
   be extended automatically (trend or periodic dummies).  However, if
   an exogenous regressor is a lag series of order p we can extend it
   automatically for p periods.

   The @dvlags argument will be non-NULL only if the specification
   includes at least one lag of the dependent variable.
*/

int addobs_can_help(Model * pmod, const(int) * dvlags, const(DatasetRaw) * dset)
{
    int i, vi, ret = 1;

    if (pmod.xlist != null) {
    for (i=1; i<=pmod.xlist[0]; i++) {
        vi = pmod.xlist[i];
        if (vi == 0) {
        /* constant, OK */
        continue;
        } else if (dvlags != null && dvlags[i-1] != 0) {
        /* a lag of the dependent variable, OK */
        continue;
        } else if (is_trend_variable(dset.Z[vi], dset.n)) {
        /* trend, OK */
        continue;
        } else if (is_periodic_dummy(dset.Z[vi], dset)) {
        /* a predictable dummy, OK */
        continue;
        } else if (series_get_lag(dset, vi) &&
               series_get_parent_id(dset, vi)) {
        /* a lagged exogenous term, OK up to a point */
        continue;
        } else {
        /* none of the above, not OK */
        ret = 0;
        break;
        }
    }
    }

    return ret;
}

/**
 * forecast_options_for_model:
 * @pmod: the model from which forecasts are wanted.
 * @dset: dataset struct.
 * @flags: location to receive flags from among #FcastFlags.
 * @dt2max: location to receive the last observation that can
 * be supported for a dynamic forecast.
 * @st2max: location to receive the last observation that can
 * be supported for a static forecast.
 *
 * Examines @pmod and determines which forecasting options are
 * applicable.
 */

void forecast_options_for_model(Model * pmod, const(DatasetRaw) * dset, FcastFlags * flags, int * dt2max, int * st2max)
{
    int * dvlags = null;
    int dvcheck = 0;
    int dv;

    *flags = cast(FcastFlags) 0;
    dv = gretl_model_get_depvar(pmod);

    if (pmod.ci == GretlCmdIndex.ols) {
    if (is_standard_diff(dv, dset, null)) {
        *flags |= FcastFlags.integrate_ok;
    } else if (series_is_log(dset, dv, null)) {
        *flags |= FcastFlags.expon_ok;
    }
    } else if (pmod.ci == GretlCmdIndex.nls) {
    /* we'll try winging it! */
    if (gretl_model_get_int(pmod, "dynamic") && pmod.t2 < dset.n - 1) {
        *flags |= FcastFlags.auto_ok;
    }
    return;
    }

    *dt2max = *st2max = pmod.t2;

    if (pmod.ci == GretlCmdIndex.arma || pmod.ci == GretlCmdIndex.garch) {
    *flags |= FcastFlags.dynamic_ok;
    } else if (AR_MODEL(cast(GretlCmdIndex) pmod.ci)) {
    *flags |= FcastFlags.dynamic_ok;
    } else if (dataset_is_time_series(dset)) {
        /* assume OK until shown otherwise */
    *flags |= FcastFlags.dynamic_ok;
        dvcheck = 1;
    }

    if (*flags & FcastFlags.dynamic_ok) {
    dvlags = process_lagged_depvar(pmod, dset);
        if (dvcheck && dvlags == null) {
            /* plain time-series, no dvlags found */
            *flags ^= FcastFlags.dynamic_ok;
        }
    }

    if (addobs_can_help(pmod, dvlags, dset)) {
    *flags |= FcastFlags.addobs_ok;
    }

    *dt2max = fcast_get_t2max(dvlags, pmod, dset, FcastMethod.dynamic);
    *st2max = fcast_get_t2max(dvlags, pmod, dset, FcastMethod.static_);

    if (dvlags != null) {
    free(dvlags);
    }
}

int ols_y_lag(int v, int parent_id, const(DatasetRaw) * dset)
{
    int lag = 0;

    if (series_get_transform(dset, v) == GretlCmdIndex.lags) {
    int pv = series_get_parent_id(dset, v);

    if (pv == parent_id) {
        lag = series_get_lag(dset, v);
    }
    } else {
        int pmax = dset.pd < 4 ? 4 : dset.pd;

        lag = detect_lag(parent_id, v, dset.t1, dset.t2,
                         pmax, dset);
    }

    return lag;
}

/* initialize for recursive k-step-ahead forecast */

int k_step_init(Model * pmod, const(DatasetRaw) * dset, int k, double ** py, int ** pllist)
{
    int * llist = null;
    int vy = pmod.list[1];
    int i, p;
    int nl = 0;
    int err = 0;

    llist = gretl_list_new(pmod.list[0] - 1);
    if (llist == null) {
        return GretlError.alloc;
    }

    for (i=2; i<=pmod.list[0]; i++) {
        if (pmod.list[i] > 0) {
            p = ols_y_lag(pmod.list[i], vy, dset);
            if (p > 0) {
                /* vi is lagged dependent variable */
                llist[i-1] = p;
                nl++;
            }
        }
    }

    if (nl > 0) {
        double * y = cast(double *) malloc(dset.n * (double).sizeof);

        if (y == null) {
            err = GretlError.alloc;
        } else {
            for (i=0; i<dset.n; i++) {
                y[i] = dset.Z[vy][i];
            }
            *py = y;
            *pllist = llist;
            llist = null;
        }
    }

    free(llist);

    return err;
}

int recursive_fcast_adjust_obs(Model * pmod, int * t1, int t2, int k)
{
    /* the earliest possible forecast start */
    int t1min = pmod.t1 + pmod.ncoeff + k - 1;
    int err = 0;

    if (*t1 < t1min) {
    *t1 = t1min;
    }

    /* minimal forecast range = 1, otherwise error */
    if (t2 < *t1) {
    err = GretlError.obs;
    }

    return err;
}

/* recursive k-step ahead forecasts, for models estimated via OLS */

FitResid * recursive_OLS_k_step_fcast(Model * pmod, DatasetRaw * dset, int t1, int t2, int k, int pre_n, int * err)
{
    FitResid * fr;
    int orig_t1 = dset.t1;
    int orig_t2 = dset.t2;
    double * y = null;
    int * llist = null;
    double xit, yf = NADBL;
    Model mod;
    int j, p, vi, nf;
    int i, s, t;

    if (pmod.ci != GretlCmdIndex.ols) {
    *err = GretlError.olsonly;
    return null;
    }

    if (k < 1) {
    gretl_errmsg_set("recursive forecast: steps-ahead must be >= 1");
    *err = GretlError.data;
    return null;
    }

    /* check feasibility of forecast range */
    *err = recursive_fcast_adjust_obs(pmod, &t1, t2, k);
    if (*err) {
    return null;
    }

    if (k > 1) {
        /* check the regressors for y-lags */
        *err = k_step_init(pmod, dset, k, &y, &llist);
        if (*err) {
            return null;
        }
    }

    fr = fit_resid_new_for_model(pmod, dset, t1, t2, pre_n, Gretlopt.OPT_NONE, err);
    if (*err) {
    free(y);
    free(llist);
    return null;
    }

    fr.method = FcastMethod.kstep;
    fr.k = k;

    /* sample range for initial estimation */
    dset.t1 = pmod.t1;
    dset.t2 = t1 - k; /* start of fcast range minus k */

    /* number of forecasts */
    nf = t2 - t1 + 1;

    for (t=0; t<dset.n; t++) {
    fr.actual[t] = dset.Z[pmod.list[1]][t];
    }

    /* estimate over successive sample ranges */
    for (s=0; s<nf; s++) {
    mod = lsq(pmod.list, dset, GretlCmdIndex.ols, Gretlopt.OPT_A | Gretlopt.OPT_Z);
    if (mod.errcode) {
        *err = mod.errcode;
        clear_model(&mod);
        break;
    }
    /* the first obs following the estimation sample */
    t = dset.t2 + 1;
        /* iterate over steps */
    for (j=0; j<k; j++) {
        yf = 0.0;
        for (i=0; i<mod.ncoeff; i++) {
        vi = mod.list[i+2];
        p = (llist != null)? llist[i+1] : 0;
        if (p > 0 && p <= j) {
            xit = y[t-p];
        } else {
            xit = dset.Z[vi][t];
        }
        if (na(xit)) {
            yf = NADBL;
            break;
        } else {
            yf += mod.coeff[i] * xit;
        }
        }
        if (y != null && j < k - 1) {
        y[t++] = yf;
        }
    }
    fr.fitted[t] = yf;
    if (!na(fr.actual[t]) && !na(fr.fitted[t])) {
        fr.resid[t] = fr.actual[t] - fr.fitted[t];
    }
    clear_model(&mod);
    dset.t2 += 1;
    }

    dset.t1 = orig_t1;
    dset.t2 = orig_t2;

    if (*err) {
    free_fit_resid(fr);
    fr = null;
    } else {
    fit_resid_set_dec_places(fr);
    strcpy(fr.depvar.ptr, dset.varname[pmod.list[1]]);
    }

    free(y);
    free(llist);

    return fr;
}

void fcast_get_continuous_range(const(FitResid) * fr, int * pt1, int * pt2)
{
    int t, t1 = fr.t1, t2 = fr.t2;

    for (t=t1; t<=t2; t++) {
    if (na(fr.actual[t]) || na(fr.fitted[t])) {
        t1++;
    } else {
        break;
    }
    }

    for (t=t2; t>=t1; t--) {
    if (na(fr.actual[t]) || na(fr.fitted[t])) {
        t2--;
    } else {
        break;
    }
    }

    *pt1 = t1;
    *pt2 = t2;
}
