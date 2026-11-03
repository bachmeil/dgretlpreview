/* Port of lib/src/describe.c: descriptive statistics plus some tests (the header describe.h is
   describeh.d).

   Differences from the C:
   - the anonymous enums of the file are the private CorrCov, DistanceMetric and PergmUsage;
   - function-level statics are private __gshared module variables;
   - qsort comparators are extern (C) functions;
   - C's N_("...") label arrays are arrays of D strings;
   - the Shapiro-Wilk code keeps the C's float arithmetic;
   - the `#if 0` alternatives of the C are dropped.
   Dependencies not yet ported are declared extern (C) below, in the style of the other ports. */

module gretl.describe;


import core.stdc.limits : INT_MAX;
import core.stdc.math : asin, ceil, cos, exp, floor, fmod, fabs, log, log10, pow, sin, sqrt;
import core.stdc.stdio : FILE, fprintf, fputs, sprintf, stderr;
import core.stdc.stdlib : atoi, calloc, free, malloc, qsort;
import core.stdc.string : memcpy, strcmp, strcpy, strlen, strncmp;
import core.stdc.ctype : isdigit;

import gretl.describeh;
import gretl.dataio : ntolabel;
import gretl.dataset;
import gretl.dataseth : dataset_is_panel, dataset_is_time_series, sample_size;
import gretl.genmainh : ModelDataIndex;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_errors : gretl_error_clear, gretl_errmsg_get, gretl_errmsg_set, gretl_errmsg_sprintf;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_intl : gretl_pop_c_numeric_locale, gretl_push_c_numeric_locale;
import gretl.gretl_intlh : UTF_WIDTH;
import gretl.gretl_list;
import gretl.gretl_matrix;
import gretl.gretl_matrixh;
import gretl.gretl_model : clear_model, exact_fit_check, gretl_model_get_normality_test, gretl_model_init;
import gretl.gretl_modelh : LIMDEP;
import gretl.gretl_prn : Prn, pprintf, pputc, pputs;
import gretl.gretl_string_table : SeriesTable, series_table_get_string, series_table_get_strings;
import gretl.gretl_utils;
import gretl.gretl_utilsh : floateq;
import gretl.glib : gchar;
import gretl.gmem : g_free;
import gretl.graphingh : GptFlags, PlotType;
import gretl.gstrfuncs : g_strdup, g_strdup_printf, g_utf8_strlen;
import gretl.estimate : lsq;
import gretl.libgretlh;
import gretl.libset : libset_get_int;
import gretl.libseth : SetKey;
import gretl.missing : list_adjust_sample, series_adjust_sample;
import gretl.missingh : NADBL, na;
import gretl.options : get_optval_double, get_optval_int, get_optval_string, incompatible_options;
import gretl.plotspech : GpKeySpec;
import gretl.printouth : SHORTLEN;
import gretl.pvalues;
import gretl.pvaluesh : DistCode;
import gretl.strutils;
import gretl.genmain : current_series_index, generate_matrix, generate_series;
import gretl.matrix_extra : gretl_covariance_matrix_from_varlist, gretl_matrix_print_to_prn, gretl_vector_from_series;
import gretl.gretl_cmatrix : gretl_matrix_fft;
import gretl.gretl_panel : panel_variance_info;
import gretl.modelprint : printmodel;
import gretl.printout : bufspace, get_gretl_digits, gretl_fix_exponent, max_namelen_in_list, max_obs_marker_length, maybe_trim_varname, print_freq, print_freq_test, print_obs_marker, print_xtab, text_print_vmatrix;
import gretl.plotspec : print_keypos_string;
import gretl.graphing : correlogram_plot, finalize_plot_input_file, gnuplot_graph_wanted, open_plot_input_file, periodogram_plot, plot_corrmat, plot_freq, plotname;
import gretl.gretl_bundle : gretl_bundle_copy, gretl_bundle_destroy;
import gretl.genfuncs : genr_fit_resid;
import gretl.nonparam : lockes_test, rank_sum_test, sign_test, signed_rank_test, vge_gamma_test;
import gretl.usermat : get_matrix_by_name;


/**
 * SECTION:describe
 * @short_description: descriptive statistics plus some tests
 * @title: Descriptive statistics
 * @include: libgretl.h
 *
 * Computation and printing of numerous descriptive statistics
 * along with some hypothesis tests, for example regarding the
 * normality of a data series.
 */

/* return the number of "good" (non-missing) observations in series x;
   also (optionally) write to x0 the first non-missing value */

private int good_obs(const(double) * x, int n, double * x0)
{
    int t, good = n;

    if (x0 != null) {
        *x0 = NADBL;
    }

    for (t=0; t<n; t++) {
        if (na(x[t])) {
            good--;
        } else if (x0 != null && na(*x0)) {
            *x0 = x[t];
        }
    }

    return good;
}

/**
 * gretl_minmax:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 * @min: location to receive minimum value.
 * @max: location to receive maximum value.
 *
 * Puts the minimum and maximum values of the series @x,
 * from obs @t1 to obs @t2, into the variables @min and @max.
 *
 * Returns: the number of valid observations in the given
 * data range.
 */

extern (C) int gretl_minmax(int t1, int t2, const(double) * x,
                            double * min, double * max)
{
    int t, n = 0;

    *max = *min = NADBL;

    for (t=t1; t<=t2; t++) {
        if (!na(x[t])) {
            if (n == 0) {
                *max = *min = x[t];
            } else {
                if (x[t] > *max) *max = x[t];
                if (x[t] < *min) *min = x[t];
            }
            n++;
        }
    }

    return n;
}

private int summary_minmax(int t1, int t2, const(double) * x,
                           Summary * summ, int i)
{
    int n = gretl_minmax(t1, t2, x, &summ.low[i], &summ.high[i]);

    if (summ.high[i] > 1.0e5 && summ.high[i] < 1.0e10) {
        /* 2024-05-29: add check for big integers */
        summ.ival[i] = gretl_isint(t1, t2, x);
    }

    return n;
}

/**
 * gretl_min:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 *
 * Returns: the minimum value of @x over the given range,
 * or #NADBL if no valid vaues are found.
 */

double gretl_min(int t1, int t2, const(double) * x)
{
    double min, max;

    gretl_minmax(t1, t2, x, &min, &max);
    return min;
}

/**
 * gretl_max:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 *
 * Returns: the maximum value of @x over the given range,
 * or #NADBL if no valid vaues are found.
 */

double gretl_max(int t1, int t2, const(double) * x)
{
    double min, max;

    gretl_minmax(t1, t2, x, &min, &max);
    return max;
}

/**
 * gretl_sum:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 *
 * Returns: the sum of the series @x from obs
 * @t1 to obs @t2, skipping any missing values, or #NADBL
 * in case there are no valid observations.
 */

double gretl_sum(int t1, int t2, const(double) * x)
{
    double sum = 0.0;
    int t, n = 0;

    for (t=t1; t<=t2; t++) {
        if (!(na(x[t]))) {
            sum += x[t];
            n++;
        }
    }

    return (n == 0)? NADBL : sum;
}

/**
 * gretl_mean:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 *
 * Returns: the arithmetic mean of the series @x from obs
 * @t1 to obs @t2, skipping any missing values, or #NADBL
 * in case there are no valid observations.
 */

double gretl_mean(int t1, int t2, const(double) * x)
{
    double xbar, sum = 0.0;
    int t, n = 0;

    for (t=t1; t<=t2; t++) {
        if (!(na(x[t]))) {
            sum += x[t];
            n++;
        }
    }

    if (n == 0) {
        return NADBL;
    }

    xbar = sum / n;
    sum = 0.0;

    for (t=t1; t<=t2; t++) {
        if (!(na(x[t]))) {
            sum += (x[t] - xbar);
        }
    }

    return xbar + sum / n;
}

/**
 * eval_ytest:
 * @y: reference numerical value.
 * @op: operator.
 * @test: numerical test value.
 *
 * Returns: 1 if the expression @y @yop @test (for example
 * "y = 2" or "y <= 45") evaluates as true, else 0.
 */

int eval_ytest(double y, GretlOp op, double test)
{
    int ret = 0;

    switch (op) {
    case GretlOp.eq:
        if (na(test)) {
            ret = na(y);
        } else {
            ret = (y == test);
        }
        break;
    case GretlOp.gt:
        ret = (y > test);
        break;
    case GretlOp.lt:
        ret = (y < test);
        break;
    case GretlOp.neq:
        if (na(test)) {
            ret = !na(y);
        } else {
            ret = (y != test);
        }
        break;
    case GretlOp.gte:
        ret = (y >= test);
        break;
    case GretlOp.lte:
        ret = (y <= test);
        break;
    default:
        break;
    }

    return ret;
}

/**
 * gretl_restricted_mean:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 * @y: criterion series.
 * @yop: criterion operator.
 * @yval: criterion value.
 *
 * Returns: the arithmetic mean of the series @x in the
 * range @t1 to @t2 (inclusive), but including only
 * observations where the criterion variable @y bears the
 * relationship @yop to the value @yval -- or #NADBL in case
 * there are no observations that satisfy the restriction.
 */

double gretl_restricted_mean(int t1, int t2, const(double) * x,
                             const(double) * y, GretlOp yop,
                             double yval)
{
    double xbar, sum = 0.0;
    int t, n = 0;

    for (t=t1; t<=t2; t++) {
        if (!na(x[t]) && eval_ytest(y[t], yop, yval)) {
            sum += x[t];
            n++;
        }
    }

    if (n == 0) {
        return NADBL;
    }

    xbar = sum / n;
    sum = 0.0;

    for (t=t1; t<=t2; t++) {
        if (!na(x[t]) && eval_ytest(y[t], yop, yval)) {
            sum += (x[t] - xbar);
        }
    }

    return xbar + sum / n;
}

/*
   For an array a with n elements, find the element which would be a[h]
   if the array were sorted from smallest to largest (without the need
   to do a full sort).  The algorithm is due to C. A. R. Hoare; the C
   implementation is shamelessly copied, with minor adaptations, from
   http://www.astro.northwestern.edu/~ato/um/programming/sorting.html
   (Atakan G\"urkan: link now dead).
*/

private double find_hoare(double * a, int n, int h)
{
    double w, x;
    int l = 0;
    int r = n - 1;
    int i, j;

    while (l < r) {
        x = a[h];
        i = l;
        j = r;
        while (j >= i) {
            while (a[i] < x) i++;
            while (x < a[j]) j--;
            if (i <= j) {
                w = a[i];
                a[i++] = a[j];
                a[j--] = w;
            }
        }
        if (j < h) l = i;
        if (h < i) r = j;
    }

    return a[h];
}

private double find_hoare_inexact(double * a, double p,
                                  double xmin, double xmax,
                                  double frac, int n,
                                  int hf, int hc)
{
    double tmp, high, low;
    int i;

    if (p < 0.5) {
        high = find_hoare(a, n, hc);
        tmp = xmin;
        for (i=0; i<hc; i++) {
            if (a[i] > tmp) {
                tmp = a[i];
            }
        }
        low = tmp;
    } else {
        low = find_hoare(a, n, hf);
        tmp = xmax;
        for (i=hc; i<n; i++) {
            if (a[i] < tmp) {
                tmp = a[i];
            }
        }
        high = tmp;
    }

    return low + frac * (high - low);
}

/* See https://en.wikipedia.org/wiki/Quantile, also Hyndman and Fan,
   "Sample Quantiles in Statistical Packages" (The American
   Statistician, 1996).  Return the "index" of the @p-quantile for a
   sample of size @n (which may not be an integer).
*/

private double quantile_index(int n, double p)
{
    int Qtype = libset_get_int(SetKey.quantile_type);
    double h;

    if (Qtype == 0) {
        /* \hat{Q}_6 */
        h = (n + 1) * p;
    } else if (Qtype == 1) {
        /* \hat{Q}_7 */
        h = (n - 1) * p + 1;
    } else {
        /* \hat{Q}_8 */
        h = (n + 1.0/3) * p + 1.0/3;
    }

    return h - 1; /* 0-based */
}

/**
 * gretl_quantile:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 * @p: probability.
 * @opt: may include OPT_Q to hush warning when sample is
 * too small.
 * @err: location to receive error code.
 *
 * Returns: the @p quantile of the series @x from obs
 * @t1 to obs @t2, skipping any missing values, or #NADBL
 * on failure.
 */

double gretl_quantile(int t1, int t2, const(double) * x, double p,
                      Gretlopt opt, int * err)
{
    double * a = null;
    double xmin, xmax;
    double h, ret;
    int hf, hc;
    int t, n = 0;

    if (*err) {
        /* don't compound a prior error */
        return NADBL;
    }

    if (p <= 0 || p >= 1) {
        /* sanity check */
        *err = GretlError.data;
    } else {
        n = gretl_minmax(t1, t2, x, &xmin, &xmax);
        if (n == 0) {
            *err = GretlError.data;
        }
    }

    if (!*err) {
        h = quantile_index(n, p);
        hf = cast(int) floor(h);
        hc = cast(int) ceil(h);

        if (hc == 0 || hc == n) {
            /* too few usable observations for the specified
               quantile; don't treat as fatal error */
            return NADBL;
        } else {
            a = cast(double *) malloc(n * (*a).sizeof);
            if (a == null) {
                *err = GretlError.alloc;
            }
        }
    }

    if (*err) {
        return NADBL;
    }

    n = 0;
    for (t=t1; t<=t2; t++) {
        if (!na(x[t])) {
            a[n++] = x[t];
        }
    }

    if (hf == hc) {
        /* "exact" */
        ret = find_hoare(a, n, hf);
    } else {
        ret = find_hoare_inexact(a, p, xmin, xmax, h - hf,
                                 n, hf, hc);
    }

    free(a);

    return ret;
}

/**
 * gretl_array_quantiles:
 * @a: data array (this gets re-ordered).
 * @n: length of array.
 * @p: array of probabilities (over-written by quantiles).
 * @k: number of probabilities.
 *
 * Computes @k quantiles (given by the elements of @p) for the
 * first n elements of the array @a, which is re-ordered in
 * the process. On successful exit, @p contains the quantiles.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_array_quantiles(double * a, int n, double * p, int k)
{
    double h, xmin = 0, xmax = NADBL;
    int hf, hc, i;
    int err = 0;

    if (n <= 0 || k <= 0) {
        return GretlError.data;
    }

    for (i=0; i<k; i++) {
        if (p[i] <= 0.0 || p[i] >= 1.0) {
            p[i] = NADBL;
            err = GretlError.invarg;
            break;
        }

        h = quantile_index(n, p[i]);
        hf = cast(int) floor(h);
        hc = cast(int) ceil(h);

        if (hc == 0 || hc == n) {
            p[i] = NADBL;
        } else if (hf == hc) {
            p[i] = find_hoare(a, n, hf);
        } else {
            if (na(xmax)) {
                gretl_minmax(0, n-1, a, &xmin, &xmax);
            }
            p[i] = find_hoare_inexact(a, p[i], xmin, xmax, h - hf,
                                      n, hf, hc);
        }
    }

    return err;
}

/**
 * gretl_array_quantile:
 * @a: array on which to operate.
 * @n: number of elements in @a.
 * @p: probability.
 *
 * Returns: the @p quantile of the first @n elements in @a,
 * which is re-ordered in the process, or #NADBL on
 * failure.
 */

double gretl_array_quantile(double * a, int n, double p)
{
    gretl_array_quantiles(a, n, &p, 1);
    return p;
}

/**
 * gretl_median:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 *
 * Returns: the median value of the series @x from obs
 * @t1 to obs @t2, skipping any missing values, or #NADBL
 * on failure.
 */

double gretl_median(int t1, int t2, const(double) * x)
{
    int err = 0;

    return gretl_quantile(t1, t2, x, 0.5, Gretlopt.OPT_NONE, &err);
}

/**
 * gretl_sst:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 *
 * Returns: the sum of squared deviations from the mean for
 * the series @x from obs @t1 to obs @t2, skipping any missing
 * values, or #NADBL on failure.
 */

double gretl_sst(int t1, int t2, const(double) * x)
{
    int t, n = t2 - t1 + 1;
    double sumsq, xx, xbar;

    if (n == 0) {
        return NADBL;
    }

    xbar = gretl_mean(t1, t2, x);
    if (na(xbar)) {
        return NADBL;
    }

    sumsq = 0.0;

    for (t=t1; t<=t2; t++) {
        if (!na(x[t])) {
            xx = x[t] - xbar;
            sumsq += xx * xx;
        }
    }

    return sumsq;
}

/* allows for computing MLE rather than using df correction,
   via the @asy argument */

double real_gretl_variance(int t1, int t2, const(double) * x,
                           int asy)
{
    int t, n = t2 - t1 + 1;
    double v, dx, xbar;

    if (n == 0) {
        /* null sample */
        return NADBL;
    }

    xbar = gretl_mean(t1, t2, x);
    if (na(xbar)) {
        return NADBL;
    }

    v = 0.0;
    n = 0;

    for (t=t1; t<=t2; t++) {
        if (!na(x[t])) {
            dx = x[t] - xbar;
            v += dx * dx;
            n++;
        }
    }

    if (asy && n > 0) {
        v /= n;
    } else if (n > 1) {
        v /= (n - 1);
    } else {
        v = 0.0;
    }

    return (v >= 0)? v : NADBL;
}

/**
 * gretl_variance:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 *
 * Returns: the variance of the series @x from obs
 * @t1 to obs @t2, skipping any missing values, or #NADBL
 * on failure.
 */

double gretl_variance(int t1, int t2, const(double) * x)
{
    return real_gretl_variance(t1, t2, x, 0);
}

/**
 * gretl_restricted_variance:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 * @y: criterion series.
 * @yop: criterion operator.
 * @yval: criterion value.
 *
 * Returns: the variance of the series @x from obs
 * @t1 to obs @t2, skipping any missing values and
 * observations where the series @y does not bear the
 * relationship @yop to the value @yval, or #NADBL on
 * failure.
 */

double gretl_restricted_variance(int t1, int t2, const(double) * x,
                                 const(double) * y, GretlOp yop,
                                 double yval)
{
    double sumsq, xx, xbar;
    int t, n = 0;

    xbar = gretl_restricted_mean(t1, t2, x, y, yop, yval);
    if (na(xbar)) {
        return NADBL;
    }

    sumsq = 0.0;

    for (t=t1; t<=t2; t++) {
        if (!na(x[t]) && eval_ytest(y[t], yop, yval)) {
            xx = x[t] - xbar;
            sumsq += xx * xx;
            n++;
        }
    }

    sumsq = (n > 1)? sumsq / (n - 1) : 0.0;

    return (sumsq >= 0)? sumsq : NADBL;
}

/**
 * gretl_stddev:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 *
 * Returns: the standard deviation of the series @x from obs
 * @t1 to obs @t2, skipping any missing values, or #NADBL
 * on failure.
 */

double gretl_stddev(int t1, int t2, const(double) * x)
{
    double v = gretl_variance(t1, t2, x);

    return (na(v))? v : sqrt(v);
}

/**
 * gretl_restricted_stddev:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 * @y: criterion series.
 * @yop: criterion operator.
 * @yval: criterion value.
 *
 * Returns: the standard deviation of the series @x from obs
 * @t1 to obs @t2, skipping any missing values and observations
 * where the series @y does not bear the relationship @yop to
 * the value @yval, or #NADBL on failure.
 */

double gretl_restricted_stddev(int t1, int t2, const(double) * x,
                               const(double) * y, GretlOp yop,
                               double yval)
{
    double xx = gretl_restricted_variance(t1, t2, x, y, yop, yval);

    return (na(xx))? xx : sqrt(xx);
}

/**
* gretl_long_run_variance:
* @t1: starting observation.
* @t2: ending observation.
* @x: data series.
* @m: bandwidth (<= 0 for automatic).
* @mu: mean (or NADBL to use sample mean).
*
* Returns: the long-run variance of the series @x from obs
* @t1 to obs @t2, using Bartlett kernel weights, or #NADBL
* on failure (which includes encountering missing values).
*/

double gretl_long_run_variance(int t1, int t2, const(double) * x,
                               int m, double mu)
{
    double zt, wi, xbar, s2 = 0.0;
    int use_mu = !na(mu);
    int i, t, n, order;

    if (series_adjust_sample(x, &t1, &t2) != 0) {
        return NADBL;
    }
    if ((n = t2 - t1 + 1) < 2) {
        return NADBL;
    }

    if (m < 0) {
        order = cast(int) exp(log(cast(double) n) / 3.0);
    } else {
        order = m;
    }

    if (use_mu) {
        xbar = mu;
    }

    if (use_mu && mu == 0.0) {
        for (t=t1; t<=t2; t++) {
            s2 += x[t] * x[t];
        }
    } else {
        if (!use_mu) {
            xbar = 0.0;
            for (t=t1; t<=t2; t++) {
                xbar += x[t];
            }
            xbar /= n;
        }
        for (t=t1; t<=t2; t++) {
            zt = x[t] - xbar;
            s2 += zt * zt;
        }
    }

    for (i=1; i<=order; i++) {
        wi = 1.0 - i /(cast(double) order + 1);
        for (t=t1+i; t<=t2; t++) {
            zt = x[t] - xbar;
            s2 += 2 * wi * zt * (x[t-i] - xbar);
        }
    }

    return s2 / n;
}

/**
 * gretl_covar:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 * @y: data series.
 * @missing: location to receive information on the number
 * of missing observations that were skipped, or %NULL.
 *
 * Returns: the covariance of the series @x and @y from obs
 * @t1 to obs @t2, skipping any missing values, or #NADBL
 * on failure.
 */

double gretl_covar(int t1, int t2, const(double) * x, const(double) * y,
                   int * missing)
{
    double sx, sy, sxy, xbar, ybar;
    int t, nn, n = t2 - t1 + 1;

    if (n == 0) {
        /* void sample */
        return NADBL;
    }

    nn = 0;
    sx = sy = 0.0;

    for (t=t1; t<=t2; t++) {
        if (!na(x[t]) && !na(y[t])) {
            sx += x[t];
            sy += y[t];
            nn++;
        }
    }

    if (nn < 2) {
        return NADBL;
    }

    xbar = sx / nn;
    ybar = sy / nn;
    sxy = 0.0;

    for (t=t1; t<=t2; t++) {
        if (!na(x[t]) && !na(y[t])) {
            sx = x[t] - xbar;
            sy = y[t] - ybar;
            sxy = sxy + (sx * sy);
        }
    }

    if (missing != null) {
        /* record number of missing obs */
        *missing = n - nn;
    }

    return sxy / (nn - 1);
}

/**
 * gretl_corr:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 * @y: data series.
 * @missing: location to receive information on the number
 * of missing observations that were skipped, or %NULL.
 *
 * Returns: the correlation coefficient for the series @x and @y
 * from obs @t1 to obs @t2, skipping any missing values, or #NADBL
 * on failure.
 */

double gretl_corr(int t1, int t2, const(double) * x, const(double) * y,
                  int * missing)
{
    int t, nn, n = t2 - t1 + 1;
    double sx, sy, sxx, syy, sxy, den, xbar, ybar;
    double cval = 0.0;

    if (n == 0) {
        /* void sample */
        return NADBL;
    }

    if (gretl_isconst(t1, t2, x) || gretl_isconst(t1, t2, y)) {
        /* correlation between any x and a constant is
           undefined */
        return NADBL;
    }

    nn = 0;
    sx = sy = 0.0;

    for (t=t1; t<=t2; t++) {
        if (!na(x[t]) && !na(y[t])) {
            sx += x[t];
            sy += y[t];
            nn++;
        }
    }

    if (nn < 2) {
        /* zero or 1 observations: no go! */
        return NADBL;
    }

    xbar = sx / nn;
    ybar = sy / nn;
    sxx = syy = sxy = 0.0;

    for (t=t1; t<=t2; t++) {
        if (!na(x[t]) && !na(y[t])) {
            sx = x[t] - xbar;
            sy = y[t] - ybar;
            sxx += sx * sx;
            syy += sy * sy;
            sxy += sx * sy;
        }
    }

    if (sxy != 0.0) {
        den = sxx * syy;
        if (den > 0.0) {
            cval = sxy / sqrt(den);
        } else {
            cval = NADBL;
        }
    }

    if (missing != null) {
        /* record number of missing observations */
        *missing = n - nn;
    }

    return cval;
}

/**
 * gretl_corr_rsq:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 * @y: data series.
 *
 * Returns: the square of the correlation coefficient for the series
 * @x and @y from obs @t1 to obs @t2, skipping any missing values,
 * or #NADBL on failure.  Used as alternative value for R^2 in a
 * regression without an intercept.
 */

double gretl_corr_rsq(int t1, int t2, const(double) * x, const(double) * y)
{
    double r = gretl_corr(t1, t2, x, y, null);

    return (na(r))? r : r * r;
}

/* we're supposing a variance smaller than this is just noise */
private enum TINYVAR = 1.0e-36;

double gretl_skewness(int t1, int t2, const(double) * x)
{
    double s, sd, xx, xbar, ret;
    int t, n = 0;

    xbar = gretl_mean(t1, t2, x);
    if (na(xbar)) {
        return NADBL;
    }

    s = 0.0;

    for (t=t1; t<=t2; t++) {
        if (!na(x[t])) {
            xx = x[t] - xbar;
            s += xx * xx;
            n++;
        }
    }

    s /= n;

    if (s <= TINYVAR) {
        ret = NADBL;
    } else {
        sd = sqrt(s);
        s = 0.0;

        for (t=t1; t<=t2; t++) {
            if (!na(x[t])) {
                xx = (x[t] - xbar) / sd;
                s += xx * xx * xx;
            }
        }

        ret = s/n;
    }

    return ret;
}

/* note: the following computes EXCESS kurtosis */

double gretl_kurtosis(int t1, int t2, const(double) * x)
{
    double s, sd, xx, xbar, ret;
    int t, n = 0;

    xbar = gretl_mean(t1, t2, x);
    if (na(xbar)) {
        return NADBL;
    }

    s = 0.0;

    for (t=t1; t<=t2; t++) {
        if (!na(x[t])) {
            xx = x[t] - xbar;
            s += xx * xx;
            n++;
        }
    }

    s /= n;

    if (s <= TINYVAR) {
        ret = NADBL;
    } else {
        sd = sqrt(s);
        s = 0.0;

        for (t=t1; t<=t2; t++) {
            if (!na(x[t])) {
                xx = (x[t] - xbar) / sd;
                s += xx * xx * xx * xx;
            }
        }

        ret = s/n - 3.0;
    }

    return ret;
}

/**
 * gretl_moments:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 * @wts: weights (may be NULL).
 * @xbar: pointer to receive mean.
 * @sd: pointer to receive standard deviation.
 * @skew: pointer to receive skewness (may be NULL).
 * @kurt: pointer to receive excess kurtosis (may be NULL).
 * @k: degrees of freedom loss (generally 1).
 *
 * Calculates sample moments for series @x from obs @t1 to obs
 * @t2.
 *
 * Returns: 0 on success, 1 on error.
 */

int gretl_moments(int t1, int t2, const(double) * x,
                  const(double) * wts,
                  double * xbar, double * std,
                  double * skew, double * kurt, int k)
{
    int t, n;
    double dev, var;
    double s, s2, s3, s4;
    int allstats = 1;
    int weighted = (wts != null);
    double wt, wn = 0;

    if (skew == null && kurt == null) {
        allstats = 0;
    }

    while (na(x[t1]) && t1 <= t2) {
        /* skip missing values at start of sample */
        t1++;
    }

    if (gretl_isconst(t1, t2, x)) {
        /* no variation in x */
        *xbar = x[t1];
        *std = 0.0;
        if (allstats) {
            *skew = *kurt = NADBL;
        }
        return 1;
    }

    /* get number and sum of valid observations */
    s = 0.0;
    n = 0;
    for (t=t1; t<=t2; t++) {
        if (!na(x[t])) {
            if (weighted) {
                wt = wts[t];
                if (!na(wt) && wt != 0.0) {
                    s += wt * x[t];
                    wn += wt;
                    n++;
                }
            } else {
                s += x[t];
                n++;
            }
        }
    }

    if (n == 0) {
        /* no valid data */
        *xbar = *std = NADBL;
        if (allstats) {
            *skew = *kurt = 0.0;
        }
        return 1;
    }

    if (!weighted) {
        wn = n;
    }

    /* calculate mean and initialize other stats */
    *xbar = s / wn;
    var = 0.0;
    if (allstats) {
        *skew = *kurt = 0.0;
    }

    /* compute quantities needed for higher moments */
    s2 = s3 = s4 = 0.0;
    for (t=t1; t<=t2; t++) {
        if (!na(x[t])) {
            dev = x[t] - *xbar;
            if (weighted) {
                wt = wts[t];
                if (!na(wt) && wt != 0.0) {
                    s2 += wt * dev * dev;
                    if (allstats) {
                        s3 += wt * pow(dev, 3);
                        s4 += wt * pow(dev, 4);
                    }
                }
            } else {
                s2 += dev * dev;
                if (allstats) {
                    s3 += pow(dev, 3);
                    s4 += pow(dev, 4);
                }
            }
        }
    }


    /* variance */
    if (weighted) {
        /* as per Stata's -summary- with aweights */
        var = n * s2 / (wn * (n-1));
    } else {
        var = s2 / (wn - k);
    }

    if (var < 0.0) {
        /* in principle impossible, but this is a computer */
        *std = NADBL;
        if (allstats) {
            *skew = *kurt = NADBL;
        }
        return 1;
    }

    if (var > TINYVAR) {
        *std = sqrt(var);
    } else {
        /* screen out minuscule variance */
        *std = 0.0;
    }

    if (allstats) {
        if (var > TINYVAR) {
            *skew = (s3 / wn) / pow(s2 / wn, 1.5);
            *kurt = ((s4 / wn) / pow(s2 / wn, 2)) - 3.0; /* excess kurtosis */
        } else {
            /* if variance is effectively zero, these should be undef'd */
            *skew = *kurt = NADBL;
        }
    }

    return 0;
}

/**
 * gretl_sorted_series:
 * @v: ID number of input series.
 * @dset: dataset struct.
 * @opt: may include OPT_M to flag an error in case
 * missing values are found.
 * @n: on input, the minimum acceptable number of
 * non-missing observations; on output, the number
 * of valid observations.
 * @err: location to receive error code.
 *
 * Returns: an array containing the valid values of the
 * input series over the sample range given in @dset,
 * sorted from smallest to largest, or NULL on error.
 * An error is flagged if the number of valid observations
 * is less than that given in @n on input, or if OPT_M
 * is given and the input contains missing values.
 */

double * gretl_sorted_series(int v, const(DatasetRaw) * dset,
                             Gretlopt opt, int * n,
                             int * err)
{
    double * y = null;
    int t, k = 0;

    if (*n < 1) {
        *n = 1;
    }

    for (t=dset.t1; t<=dset.t2; t++) {
        if (!na(dset.Z[v][t])) {
            k++;
        } else if (opt & Gretlopt.OPT_M) {
            *err = GretlError.missdata;
            return null;
        }
    }

    if (k < *n) {
        gretl_errmsg_set(_("Insufficient data").ptr);
        *err = GretlError.data;
        return null;
    }

    y = cast(double *) malloc(k * (*y).sizeof);
    if (y == null) {
        *err = GretlError.alloc;
        return null;
    }

    *n = k;

    k = 0;
    for (t=dset.t1; t<=dset.t2; t++) {
        if (!na(dset.Z[v][t])) {
            y[k++] = dset.Z[v][t];
        }
    }

    qsort(y, k, (*y).sizeof, &gretl_compare_doubles);

    return y;
}

/**
 * free_freq:
 * @freq: pointer to gretl frequency distribution struct
 *
 * Frees all resources associated with @freq, and the
 * pointer itself.
 */

void free_freq(FreqDist * freq)
{
    if (freq == null) {
        return;
    }

    free(freq.midpt);
    free(freq.endpt);
    free(freq.f);

    if (freq.S != null) {
        strings_array_free(freq.S, freq.numbins);
    }

    free(freq);
}

private FreqDist * freq_new(const(DatasetRaw) * dset, int v)
{
    FreqDist * freq;

    freq = cast(FreqDist *) malloc((*freq).sizeof);
    if (freq == null) return null;

    strcpy(freq.varname.ptr, dset.varname[v]);
    strcpy(freq.gname.ptr, plotname(dset, v, 1));

    freq.midpt = null;
    freq.endpt = null;
    freq.S = null;
    freq.f = null;

    freq.dist = cast(DistCode) 0;
    freq.discrete = 0;
    freq.strvals = 0;

    freq.xbar = NADBL;
    freq.sdx = NADBL;
    freq.test = NADBL;

    return freq;
}

/*
 * dh_root_b1_to_z1:
 * @rb1: square root b1, skewness.
 * @n: number of observations.
 *
 * Performs the transformation from skewness, root b1, to the
 * normal score, z1, as set out in Doornik and Hansen, "An Omnibus
 * Test for Normality", 1994.  The transformation is originally
 * due to D'Agostino (Biometrika, 1970).
 *
 * Returns: the z1 value.
 */

private double dh_root_b1_to_z1(double rb1, double n)
{
    double b, w2, d, y, z1;

    b = 3.0 * (n*n + 27*n - 70) * (n+1) * (n+3) /
        ((n-2) * (n+5) * (n+7) * (n+9));

    w2 = -1.0 + sqrt(2 * (b-1));
    d = 1.0 / sqrt(log(sqrt(w2)));
    y = rb1 * sqrt(((w2-1.0)/2.0) * ((n+1.0)*(n+3.0))/(6.0*(n-2)));
    z1 = d * log(y + sqrt(y*y + 1));

    return z1;
}

/*
 * dh_b2_to_z2:
 * @b1: skewness.
 * @b2: kurtosis.
 * @n: number of observations.
 *
 * Performs the transformation from kurtosis, b2, to the
 * normal score, z2, as set out in Doornik and Hansen, "An Omnibus
 * Test for Normality", 1994.
 *
 * Returns: the z2 value, or #NADBL on failure.
 */

private double dh_b2_to_z2(double b1, double b2, double n)
{
    double d, a, c, k, alpha, chi, z2;
    double n2 = n * n;

    d = (n-3) * (n+1) * (n2 + 15*n - 4.0);
    a = ((n-2) * (n+5) * (n+7) * (n2 + 27*n - 70.0)) / (6.0 * d);
    c = ((n-7) * (n+5) * (n+7) * (n2 + 2*n - 5.0)) / (6.0 * d);
    k = ((n+5) * (n+7) * (n2*n + 37*n2 + 11*n - 313.0)) / (12.0 * d);

    alpha = a + b1 * c;
    z2 =  (1.0 / (9.0 * alpha)) - 1.0;
    chi = (b2 - 1.0 - b1) * 2.0 * k;

    if (chi > 0.0) {
       z2 += pow(chi/(2*alpha), 1.0 / 3.0);
    }

    z2 *= sqrt(9.0*alpha);

    return z2;
}

/**
 * doornik_chisq:
 * @skew: skewness.
 * @xkurt: excess kurtosis.
 * @n: number of observations.
 *
 * Calculates the Chi-square test for normality as set out by
 * Doornik and Hansen, "An Omnibus Test for Normality", 1994.
 * This is a modified version of the test proposed by Bowman and
 * Shenton (Biometrika, 1975).
 *
 * Returns: the Chi-square value, which has 2 degrees of freedom.
 */

double doornik_chisq(double skew, double xkurt, int n)
{
    double z1, z2, x2 = NADBL;

    z1 = dh_root_b1_to_z1(skew, cast(double) n);
    z2 = dh_b2_to_z2(skew * skew, xkurt + 3.0, cast(double) n);

    if (!na(z2)) {
        x2 = z1*z1 + z2*z2;
    }

    return x2;
}

int series_get_skew_kurt(int t1, int t2, const(double) * x,
                         double * skew, double * xkurt,
                         int * pn)
{
    double dev;
    double[4] s = 0.0;
    int t, n = 0;
    int err = 0;

    for (t=t1; t<=t2; t++) {
        if (!na(x[t])) {
            s[0] += x[t];
            n++;
        }
    }

    *pn = n;
    if (n < 3) {
        return GretlError.data;
    }

    s[0] /= n;

    for (t=t1; t<=t2; t++) {
        if (!na(x[t])) {
            dev = x[t] - s[0];
            s[1] += dev * dev;
            s[2] += pow(dev, 3);
            s[3] += pow(dev, 4);
        }
    }

    s[1] /= n;

    if (s[1] > 0.0) {
        *skew = (s[2] / n) / pow(s[1], 1.5);
        *xkurt = ((s[3] / n) / pow(s[1], 2)) - 3.0;
    } else {
        /* if variance is effectively zero, these should be undef'd */
        *skew = *xkurt = NADBL;
        err = 1;
    }

    return err;
}

private int get_moments(const(GretlMatrixRaw) * M, int row, double * skew, double * kurt)
{
    int j, n = gretl_matrix_cols(M);
    double dev;
    double[4] s = 0.0;
    int err = 0;

    for (j=0; j<n; j++) {
        s[0] += gretl_matrix_get(M, row, j);
    }

    s[0] /= n;

    for (j=0; j<n; j++) {
        dev = gretl_matrix_get(M, row, j) - s[0];
        s[1] += dev * dev;
        s[2] += pow(dev, 3);
        s[3] += pow(dev, 4);
    }

    s[1] /= n;

    if (s[1] > 0.0) {
        *skew = (s[2] / n) / pow(s[1], 1.5);
        *kurt = ((s[3] / n) / pow(s[1], 2));
    } else {
        /* if variance is effectively zero, these should be undef'd */
        *skew = *kurt = NADBL;
        err = 1;
    }

    return err;
}

int multivariate_normality_test(const(GretlMatrixRaw) * E,
                                const(GretlMatrixRaw) * Sigma,
                                Gretlopt opt, Prn * prn)
{
    GretlMatrixRaw * S = null;
    GretlMatrixRaw * V = null;
    GretlMatrixRaw * C = null;
    GretlMatrixRaw * X = null;
    GretlMatrixRaw * R = null;
    GretlMatrixRaw * evals = null;
    GretlMatrixRaw * tmp = null;
    /* convenience pointers: do not free! */
    GretlMatrixRaw * H;
    GretlVector * Z1;
    GretlVector * Z2;
    double x, skew, kurt;
    double X2 = NADBL;
    int n, p;
    int i, j;
    int err = 0;

    if (E == null || Sigma == null) {
        err = GretlError.data;
        goto bailout;
    }

    p = gretl_matrix_cols(E);
    n = gretl_matrix_rows(E);

    clear_gretl_matrix_err();

    S = gretl_matrix_copy(Sigma);
    V = gretl_vector_alloc(p);
    C = gretl_matrix_alloc(p, p);
    X = gretl_matrix_copy_transpose(E);
    R = gretl_matrix_alloc(p, n);
    tmp = gretl_matrix_alloc(p, p);

    err = get_gretl_matrix_err();
    if (err) {
        goto bailout;
    }

    for (i=0; i<p; i++) {
        x = gretl_matrix_get(S, i, i);
        gretl_vector_set(V, i, 1.0 / sqrt(x));
    }

    err = gretl_matrix_diagonal_sandwich(V, S, C);
    if (err) {
        goto bailout;
    }

    if (!(opt & Gretlopt.OPT_Q)) {
        pputc(prn, '\n');
        gretl_matrix_print_to_prn(C, _("Residual correlation matrix, C").ptr, prn);
    }

    evals = gretl_symmetric_matrix_eigenvals(C, 1, &err);
    if (err) {
        goto bailout;
    }

    if (!(opt & Gretlopt.OPT_Q)) {
        pprintf(prn, "%s\n\n", _("Eigenvalues of C").ptr);
        for (i=0; i<p; i++) {
            pprintf(prn, " %10g\n", evals.val[i]);
        }
        pputc(prn, '\n');
    }

    /* C should now contain eigenvectors of the original C:
       relabel as 'H' for perspicuity */
    H = C;
    gretl_matrix_copy_values(tmp, H);

    /* make "tmp" into $H \Lambda^{-1/2}$ */
    for (i=0; i<p; i++) {
        for (j=0; j<p; j++) {
            x = gretl_matrix_get(tmp, i, j);
            x *= 1.0 / sqrt(evals.val[j]);
            gretl_matrix_set(tmp, i, j, x);
        }
    }

    /* make S into $H \Lambda^{-1/2} H'$ */
    gretl_matrix_multiply_mod(tmp, GretlMatrixMod.none,
                              H, GretlMatrixMod.transpose,
                              S, GretlMatrixMod.none);

    gretl_matrix_demean_by_row(X);

    /* compute VX', in X (don't need to subtract means, because these
       are OLS residuals)
    */
    for (i=0; i<p; i++) {
        for (j=0; j<n; j++) {
            x = gretl_matrix_get(X, i, j);
            x *= gretl_vector_get(V, i);
            gretl_matrix_set(X, i, j, x);
        }
    }

    /* finally, compute $R' = H  \Lambda^{-1/2} H' V X'$
       Doornik and Hansen, 1994, section 3 */
    gretl_matrix_multiply(S, X, R);

    /* Z_1 and Z_2 are p-vectors: use existing storage */
    Z1 = V;
    Z2 = gretl_matrix_reuse(tmp, p, 1);

    for (i=0; i<p && !err; i++) {
        get_moments(R, i, &skew, &kurt);
        if (na(skew) || na(kurt)) {
            err = 1;
        } else {
            double z1i = dh_root_b1_to_z1(skew, n);
            double z2i = dh_b2_to_z2(skew * skew, kurt, n);

            if (na(z2i)) {
                X2 = NADBL;
                err = GretlError.nan;
            } else {
                gretl_vector_set(Z1, i, z1i);
                gretl_vector_set(Z2, i, z2i);
            }
        }
    }

    if (!err) {
        X2 = gretl_vector_dot_product(Z1, Z1, &err);
        X2 += gretl_vector_dot_product(Z2, Z2, &err);
    }

    if (na(X2)) {
        pputs(prn, "Calculation of test statistic failed\n");
    } else {
        /* print and record result */
        double pv = chisq_cdf_comp(2 * p, X2);

        if (!(opt & Gretlopt.OPT_I)) {
            pputs(prn, _("Doornik-Hansen test").ptr);
            pprintf(prn, "\n %s(%d) = %g [%.4f]\n\n", _("Chi-square").ptr, 2 * p,
                    X2, pv);
        }
        record_test_result(X2, pv);
    }

 bailout:

    gretl_matrix_free(S);
    gretl_matrix_free(V);
    gretl_matrix_free(C);
    gretl_matrix_free(X);
    gretl_matrix_free(R);
    gretl_matrix_free(evals);
    gretl_matrix_free(tmp);

    return err;
}

private int freq_add_arrays(FreqDist * freq, int n)
{
    int err = 0;

    if (!freq.discrete) {
        freq.endpt = cast(double *) malloc((n + 1) * (*freq.endpt).sizeof);
        if (freq.endpt == null) {
            err = GretlError.alloc;
        }
    }

    if (!err) {
        if (freq.strvals) {
            freq.S = strings_array_new(n);
            if (freq.S == null) {
                err = GretlError.alloc;
            }
        } else {
            freq.midpt = cast(double *) malloc(n * (*freq.midpt).sizeof);
            if (freq.midpt == null) {
                err = GretlError.alloc;
            }
        }
    }

    if (!err) {
        freq.f = cast(int *) malloc(n * (*freq.f).sizeof);
        if (freq.f == null) {
            err = GretlError.alloc;
        }
    }

    if (!err) {
        freq.numbins = n;
    }

    return err;
}

private int freq_n_bins(int k0, int n, double s, double r)
{
    int k;

    if (k0 > 0) {
        /* got a pre-specified value */
        k = k0;
    } else if (n < 16) {
        k = 5;
    } else if (n < 50) {
        k = 7;
    } else if (n > 850) {
        k = 29;
    } else {
        k = cast(int) sqrt(cast(double) n);
        if (k > 99) {
            k = 99;
        }
    }
    if (k % 2 == 0) {
        k++;
    }

    return k;
}

int freq_setup(int v, const(DatasetRaw) * dset, int * pn,
               double * pxmax, double * pxmin, int * nbins,
               double * binwidth)
{
    const(double) * x = dset.Z[v];
    GretlMatrixRaw * xvals;
    int all_n = sample_size(dset);
    int nv, missvals = 0;
    int n = 0;
    int err = 0;

    xvals = gretl_matrix_values_full(x + dset.t1, all_n, Gretlopt.OPT_S,
                                     &missvals, &err);
    if (!err) {
        nv = xvals.rows;
        n = all_n - missvals;
        if (nv < 2) {
            gretl_errmsg_sprintf(_("%s is a constant").ptr, dset.varname[v]);
            err = GretlError.data;
        } else if (n < 8) {
            gretl_errmsg_sprintf(_("Insufficient data to build frequency "
                                   ~ "distribution for variable %s").ptr,
                                 dset.varname[v]);
            err = GretlError.toofew;
        }
    }

    if (!err) {
        double xmin = xvals.val[0];
        double xmax = xvals.val[nv-1];
        double xrange = xmax - xmin;
        double s = gretl_stddev(dset.t1, dset.t2, x);
        int k = freq_n_bins(*nbins, n, s, xrange);

        if (nv < k) {
            /* too few distinct values: don't do binning */
            *nbins = 0;
        } else {
            *nbins = k;
            *binwidth = xrange / (k - 1);
        }
        *pn = n;
        *pxmax = xmax;
        *pxmin = xmin;
    }

    gretl_matrix_free(xvals);

    return err;
}

/* calculate test stat for distribution, if the sample
   is big enough */

private void freq_dist_stat(FreqDist * freq, const(double) * x, Gretlopt opt, int k)
{
    double skew, kurt;

    freq.test = NADBL;
    freq.dist = DistCode.none;

    gretl_moments(freq.t1, freq.t2, x, null,
                  &freq.xbar, &freq.sdx,
                  &skew, &kurt, k);

    if (freq.n > 7) {
        if (opt & Gretlopt.OPT_O) {
            int err = 0;

            if (freq.n > 500) {
                freq.test = vge_gamma_test(x, freq.t1, freq.t2, &err);
            } else {
                freq.test = lockes_test(x, freq.t1, freq.t2, &err);
            }
            freq.dist = DistCode.gamma;
        } else if (opt & Gretlopt.OPT_Z) {
            freq.test = doornik_chisq(skew, kurt, freq.n);
            freq.dist = DistCode.normal;
        }
    }
}

private struct StrvalSorter {
    char * s;
    int n;
}

extern (C) int compare_strvals_freq(const(void) * a, const(void) * b)
{
    const(StrvalSorter) * sa = cast(const(StrvalSorter) *) a;
    const(StrvalSorter) * sb = cast(const(StrvalSorter) *) b;

    return strcmp(sa.s, sb.s);
}

FreqDist * get_string_freq(int v, const(DatasetRaw) * dset,
                           int * err)
{
    FreqDist * freq;
    StrvalSorter * ss;
    const(double) * x = dset.Z[v];
    SeriesTable * st;
    char ** S;
    int i, t, n, ns;

    freq = freq_new(dset, v);
    if (freq == null) {
        *err = GretlError.alloc;
        return null;
    }

    st = series_get_string_table(dset, v);
    S = series_table_get_strings(st, &ns);

    ss = cast(StrvalSorter *) malloc(ns * (*ss).sizeof);
    if (ss == null) {
        *err = GretlError.alloc;
        free(freq);
        return null;
    }

    for (i=0; i<ns; i++) {
        ss[i].s = S[i];
        ss[i].n = 0;
    }

    n = 0;
    for (t=dset.t1; t<=dset.t2; t++) {
        if (!na(x[t])) {
            i = cast(int) (x[t] - 1);
            if (i >= 0 && i < ns) {
                ss[i].n += 1;
                n++;
            }
        }
    }

    qsort(ss, ns, (*ss).sizeof, &compare_strvals_freq);

    freq.t1 = dset.t1;
    freq.t2 = dset.t2;
    freq.n = n;
    freq.discrete = 1;
    freq.strvals = 1;

    if (freq_add_arrays(freq, ns)) {
        *err = GretlError.alloc;
    } else {
        for (i=0; i<ns; i++) {
            freq.S[i] = gretl_strdup(ss[i].s);
            freq.f[i] = ss[i].n;
        }
    }

    free(ss);

    if (*err && freq != null) {
        free_freq(freq);
        freq = null;
    }

    return freq;
}

FreqDist * get_discrete_freq(int v, const(DatasetRaw) * dset,
                             Gretlopt opt, int * err)
{
    FreqDist * freq;
    const(double) * x = dset.Z[v];
    int * ifreq = null;
    double * ivals = null;
    double * sorted = null;
    double last;
    int i, t, nv;

    freq = freq_new(dset, v);
    if (freq == null) {
        *err = GretlError.alloc;
        return null;
    }

    freq.t1 = dset.t1;
    freq.t2 = dset.t2;

    freq.n = 0;
    for (t=freq.t1; t<=freq.t2; t++) {
        if (!na(x[t])) {
            freq.n += 1;
        }
    }

    if (freq.n < 3) {
        gretl_errmsg_sprintf(_("Insufficient data to build frequency "
                               ~ "distribution for variable %s").ptr,
                             dset.varname[v]);
        *err = GretlError.toofew;
        goto bailout;
    }

    freq.discrete = 1;
    freq.test = NADBL;
    freq.dist = DistCode.none;

    sorted = cast(double *) malloc(freq.n * (*sorted).sizeof);
    if (sorted == null) {
        *err = GretlError.alloc;
        goto bailout;
    }

    i = 0;
    for (t=freq.t1; t<=freq.t2; t++) {
        if (!na(x[t])) {
            sorted[i++] = x[t];
        }
    }

    qsort(sorted, freq.n, (*sorted).sizeof, &gretl_compare_doubles);
    nv = count_distinct_values(sorted, freq.n);

    if (nv >= 10 && !(opt & Gretlopt.OPT_X)) {
        freq_dist_stat(freq, x, opt, 1);
    } else if (opt & (Gretlopt.OPT_Z | Gretlopt.OPT_O)) {
        freq.xbar = gretl_mean(freq.t1, freq.t2, x);
        freq.sdx = gretl_stddev(freq.t1, freq.t2, x);
    }

    ifreq = cast(int *) malloc(nv * (*ifreq).sizeof);
    ivals = cast(double *) malloc(nv * (*ivals).sizeof);
    if (ifreq == null || ivals == null) {
        *err = GretlError.alloc;
        goto bailout;
    }

    ivals[0] = last = sorted[0];
    ifreq[0] = i = 1;

    for (t=1; t<freq.n; t++) {
        if (sorted[t] != last) {
            last = sorted[t];
            ifreq[i] = 1;
            ivals[i++] = last;
        } else {
            ifreq[i-1] += 1;
        }
    }

    if (freq_add_arrays(freq, nv)) {
        *err = GretlError.alloc;
    } else {
        int allints = 1;

        for (i=0; i<nv; i++) {
            if (allints && ivals[i] != floor(ivals[i])) {
                allints = 0;
            }
            freq.midpt[i] = ivals[i];
            freq.f[i] = ifreq[i];
        }

        if (allints) {
            freq.discrete = 2;
        }
    }

 bailout:

    free(sorted);
    free(ivals);
    free(ifreq);

    if (*err && freq != null) {
        free_freq(freq);
        freq = null;
    }

    return freq;
}

/**
 * get_freq:
 * @varno: ID number of variable to process.
 * @dset: dataset struct.
 * @fmin: lower limit of left-most bin (or #NADBL for automatic).
 * @fwid: bin width (or #NADBL for automatic).
 * @nbins: number of bins to use (or 0 for automatic).
 * @params: degrees of freedom loss (generally = 1 unless we're dealing
 * with the residual from a regression).
 * @opt: if includes %OPT_Z, set up for comparison with normal dist;
 * if includes %OPT_O, compare with gamma distribution;
 * if includes %OPT_S we're not printing results; if includes %OPT_D,
 * treat the variable as discrete; %OPT_X indicates that this function
 * is called as part of a cross-tabulation; if includes %OPT_G,
 * add the arrays that are needed for a plot even when we're not
 * printing (%OPT_S).
 * @err: location to receive error code.
 *
 * Calculates the frequency distribution for the specified variable.
 *
 * Returns: pointer to struct containing the distribution.
 */

FreqDist * get_freq(int varno, const(DatasetRaw) * dset,
                    double fmin, double fwid,
                    int nbins, int params,
                    Gretlopt opt, int * err)
{
    FreqDist * freq;
    const(double) * x;
    double xx, xmin, xmax;
    double binwidth = fwid;
    int t, k, n;

    if (is_string_valued(dset, varno)) {
        return get_string_freq(varno, dset, err);
    } else if (series_is_discrete(dset, varno) || (opt & Gretlopt.OPT_D)) {
        return get_discrete_freq(varno, dset, opt, err);
    } else if (gretl_isdiscrete(dset.t1, dset.t2, dset.Z[varno]) > 1) {
        return get_discrete_freq(varno, dset, opt, err);
    }

    freq = freq_new(dset, varno);
    if (freq == null) {
        *err = GretlError.alloc;
        return null;
    }

    *err = freq_setup(varno, dset, &n, &xmax, &xmin, &nbins, &binwidth);

    if (!*err && nbins == 0) {
        /* no binning wanted: switch to discrete variant */
        free_freq(freq);
        return get_discrete_freq(varno, dset, opt, err);
    }

    if (*err) {
        goto bailout;
    }

    if (!na(fmin) && !na(fwid)) {
        /* endogenous implied number of bins */
        nbins = cast(int) ceil((xmax - fmin) / fwid);
        if (nbins <= 0 || nbins > 5000) {
            *err = GretlError.invarg;
            goto bailout;
        } else {
            binwidth = fwid;
        }
    }

    freq.t1 = dset.t1;
    freq.t2 = dset.t2;
    freq.n = n;

    x = dset.Z[varno];
    freq_dist_stat(freq, x, opt, params);

    if (freq_add_arrays(freq, nbins)) {
        *err = GretlError.alloc;
        goto bailout;
    }

    if (na(fmin) || na(fwid)) {
        freq.endpt[0] = xmin - .5 * binwidth;
        if (xmin > 0.0 && freq.endpt[0] < 0.0) {
            double rshift;

            freq.endpt[0] = 0.0;
            rshift = 1.0 - xmin / binwidth;
            freq.endpt[freq.numbins] = xmax + rshift * binwidth;
        } else {
            freq.endpt[freq.numbins] = xmax + .5 * binwidth;
        }
    } else {
        freq.endpt[0] = fmin;
        freq.endpt[freq.numbins] = fmin + nbins * fwid;
    }

    for (k=0; k<freq.numbins; k++) {
        freq.f[k] = 0;
        if (k > 0) {
            freq.endpt[k] = freq.endpt[k-1] + binwidth;
        }
        freq.midpt[k] = freq.endpt[k] + .5 * binwidth;
    }

    for (t=dset.t1; t<=dset.t2; t++) {
        xx = x[t];
        if (na(xx)) {
            continue;
        }
        if (xx < freq.endpt[1]) {
            freq.f[0] += 1;
        } else if (xx >= freq.endpt[freq.numbins]) {
            freq.f[freq.numbins - 1] += 1;
            continue;
        } else {
            for (k=1; k<freq.numbins; k++) {
                if (freq.endpt[k] <= xx && xx < freq.endpt[k+1]) {
                    freq.f[k] += 1;
                    break;
                }
            }
        }
    }

 bailout:

    if (*err && freq != null) {
        free_freq(freq);
        freq = null;
    }

    return freq;
}

private void record_freq_test(const(FreqDist) * freq)
{
    double pval = NADBL;

    if (freq.dist == DistCode.normal) {
        pval = chisq_cdf_comp(2, freq.test);
    } else if (freq.dist == DistCode.gamma) {
        pval = normal_pvalue_2(freq.test);
    }

    if (!na(pval)) {
        record_test_result(freq.test, pval);
    }
}

private int check_freq_opts(Gretlopt opt, int * n_bins,
                            double * fmin, double * fwid)
{
    double x;
    int err = 0;

    if (opt & Gretlopt.OPT_N) {
        /* number of bins given */
        if (opt & (Gretlopt.OPT_M | Gretlopt.OPT_W)) {
            /* can't give min, width spec */
            return GretlError.badopt;
        } else {
            int n = get_optval_int(GretlCmdIndex.freq, Gretlopt.OPT_N, &err);

            if (!err) {
                if (n < 2 || n > 10000) {
                    err = GretlError.invarg;
                } else {
                    *n_bins = n;
                }
            }
            if (err) {
                return err;
            }
        }
    }

    if (opt & Gretlopt.OPT_M) {
        /* minimum specified */
        if (!(opt & Gretlopt.OPT_W)) {
            /* but no width given */
            return GretlError.args;
        } else {
            x = get_optval_double(GretlCmdIndex.freq, Gretlopt.OPT_M, &err);
            if (err) {
                return err;
            } else if (na(x)) {
                return GretlError.args;
            } else {
                *fmin = x;
            }
        }
    }

    if (opt & Gretlopt.OPT_W) {
        /* width specified */
        if (!(opt & Gretlopt.OPT_M)) {
            /* but no min given */
            return GretlError.args;
        } else {
            x = get_optval_double(GretlCmdIndex.freq, Gretlopt.OPT_W, &err);
            if (err) {
                return err;
            } else if (na(x)) {
                return GretlError.args;
            } else if (x <= 0) {
                return GretlError.invarg;
            } else {
                *fwid = x;
            }
        }
    }

    return 0;
}

private void record_freq_matrix(FreqDist * fd)
{
    GretlMatrixRaw * m = null;
    char ** Sr = null;
    int i, n = fd.numbins;

    if (fd.S != null) {
        m = gretl_matrix_alloc(n, 1);
        Sr = strings_array_dup(fd.S, n);
    } else {
        m = gretl_matrix_alloc(n, 2);
    }

    if (m != null) {
        char ** Sc = null;

        if (fd.S == null) {
            Sc = strings_array_new(2);
            if (Sc != null) {
                Sc[0] = gretl_strdup("midpoint");
                Sc[1] = gretl_strdup("count");
            }
        }
        for (i=0; i<n; i++) {
            if (fd.S != null) {
                gretl_matrix_set(m, i, 0, fd.f[i]);
            } else {
                gretl_matrix_set(m, i, 0, fd.midpt[i]);
                gretl_matrix_set(m, i, 1, fd.f[i]);
            }
        }
        if (Sc != null) {
            gretl_matrix_set_colnames(m, Sc);
        }
        if (Sr != null) {
            gretl_matrix_set_rownames(m, Sr);
        }
        set_last_result_data(m, GretlType.matrix);
    }
}

/* Wrapper function: get the distribution, print it if wanted, graph
   it if wanted, then free stuff.
*/

int freqdist(int varno, const(DatasetRaw) * dset,
             Gretlopt opt, Prn * prn)
{
    FreqDist * freq = null;
    DistCode dist = DistCode.none;
    PlotType ptype = cast(PlotType) 0;
    double fmin = NADBL;
    double fwid = NADBL;
    int n_bins = 0;
    int do_graph = 0;
    int err = 0;

    if (opt & Gretlopt.OPT_O) {
        dist = DistCode.gamma;
        ptype = PlotType.freq_gamma;
    } else if (opt & Gretlopt.OPT_Z) {
        dist = DistCode.normal;
        ptype = PlotType.freq_normal;
    } else {
        ptype = PlotType.freq_simple;
    }

    do_graph = gnuplot_graph_wanted(ptype, opt, &err);

    if (!err) {
        err = check_freq_opts(opt, &n_bins, &fmin, &fwid);
    }

    if (!err) {
        Gretlopt fopt = opt;

        if (do_graph) {
            fopt |= Gretlopt.OPT_G;
        }
        freq = get_freq(varno, dset, fmin, fwid, n_bins, 1, fopt, &err);
    }

    if (!err) {
        if (!(opt & Gretlopt.OPT_S)) {
            print_freq(freq, varno, dset, prn);
        } else if (dist) {
            record_freq_test(freq);
        }
        record_freq_matrix(freq);

        if (do_graph && freq.numbins < 2) {
            do_graph = 0;
        }
        if (do_graph) {
            int gerr = plot_freq(freq, dist);

            if (gerr) {
                const(char) * msg = gretl_errmsg_get();

                if (msg != null) {
                    pprintf(prn, "%s\n", msg);
                } else {
                    pputs(prn, _("gnuplot command failed\n").ptr);
                }
                do_graph = 0;
            }
        }
        free_freq(freq);
    }

    return err;
}

/* Wrapper function: get the frequency distribution and write it into
   a gretl_matrix: the first column holds the mid-points of the bins
   and the second holds the frequencies.
*/

GretlMatrixRaw * freqdist_matrix(const(double) * x, int t1, int t2, int * err)
{
    DatasetRaw * dset = null;
    FreqDist * freq = null;
    GretlMatrixRaw * m = null;
    int i, t, T = t2 - t1 + 1;

    dset = create_auxiliary_dataset(1, T, Gretlopt.OPT_NONE);
    if (dset == null) {
        *err = GretlError.alloc;
    }

    if (!*err) {
        i = 0;
        for (t=t1; t<=t2; t++) {
            dset.Z[0][i++] = x[t];
        }
        freq = get_freq(0, dset, NADBL, NADBL, 0, 1, Gretlopt.OPT_NONE, err);
    }

    if (!*err) {
        m = gretl_matrix_alloc(freq.numbins, 2);
        if (m == null) {
            *err = GretlError.alloc;
        }
    }

    if (!*err) {
        for (i=0; i<freq.numbins; i++) {
            gretl_matrix_set(m, i, 0, freq.midpt[i]);
            gretl_matrix_set(m, i, 1, freq.f[i]);
        }
    }

    destroy_dataset(dset);
    free_freq(freq);

    return m;
}

private double alt_mget(const(GretlMatrixRaw) * m,
                        int i, int j, int tr)
{
    return tr ? m.val[i*m.rows+j] : m.val[j*m.rows+i];
}

private int corresp_status(const(GretlMatrixRaw) * H, int tr)
{
    int dim1 = tr ? H.rows : H.cols;
    int dim2 = tr ? H.cols : H.rows;
    int sum1, max1 = 0;
    int sum2, max2 = 0;
    int i, j;
    int ret = 0;

    /* find the max sum in dimension 1 */
    for (j=0; j<dim1 && max1 < 2; j++) {
        sum1 = 0;
        for (i=0; i<dim2 && sum1 < 2; i++) {
            sum1 += alt_mget(H, i, j, tr) > 0;
        }
        if (sum1 > max1) {
            max1 = sum1;
        }
    }
    if (max1 == 1) {
        ret = tr ? -1 : 1;
        /* find the max sum in dimension 2 */
        for (i=0; i<dim2 && max2 < 2; i++) {
            sum2 = 0;
            for (j=0; j<dim1 && sum2 < 2; j++) {
                sum2 += alt_mget(H, i, j, tr) > 0;
            }
            if (sum2 > max2) {
                max2 = sum2;
            }
        }
        if (max2 == 1) {
            ret = 2;
        }
    }

    return ret;
}

int correspondence(const(double) * x, const(double) * y,
                   int n, int * err)
{
    GretlMatrixRaw * H = null;
    int status = 0;

    /* status: 2 means 1-to-1 x:y relationship
               1 means 1-to-n x:y relationship
              -1 means n-to-1 x:y relationship
               0 means no relationship
    */

    H = gretl_matrix_xtab(x, y, n, err);

    if (!*err) {
        status = corresp_status(H, 0);
        if (status == 0) {
            status = corresp_status(H, 1);
        }
    }

    gretl_matrix_free(H);

    return status;
}

GretlMatrixRaw * xtab_to_matrix(const(Xtab) * tab)
{
    GretlMatrixRaw * m;
    double x;
    int i, j;

    if (tab == null) {
        return null;
    }

    m = gretl_matrix_alloc(tab.rows, tab.cols);
    if (m == null) {
        return null;
    }

    for (j=0; j<tab.cols; j++) {
        for (i=0; i<tab.rows; i++) {
            x = cast(double) tab.f[i][j];
            gretl_matrix_set(m, i, j, x);
        }
    }

    return m;
}

__gshared private void * last_result;
__gshared private GretlType last_result_type;

void * get_last_result_data(GretlType * type, int * err)
{
    void * ret = null;

    if (last_result == null) {
        *type = GretlType.none;
        *err = GretlError.badstat;
    } else {
        *type = last_result_type;
        if (*type == GretlType.matrix) {
            ret = gretl_matrix_copy(cast(GretlMatrixRaw *) last_result);
        } else {
            ret = gretl_bundle_copy(cast(GretlBundle *) last_result, err);
        }
    }

    return ret;
}

void set_last_result_data(void * data, GretlType type)
{
    if (last_result != null) {
        if (last_result_type == GretlType.matrix) {
            gretl_matrix_free(cast(GretlMatrixRaw *) last_result);
        } else if (last_result_type == GretlType.bundle) {
            gretl_bundle_destroy(cast(GretlBundle *) last_result);
        }
    }

    last_result = data;
    last_result_type = type;
}

void last_result_cleanup()
{
    if (last_result != null) {
        if (last_result_type == GretlType.matrix) {
            gretl_matrix_free(cast(GretlMatrixRaw *) last_result);
        } else if (last_result_type == GretlType.bundle) {
            gretl_bundle_destroy(cast(GretlBundle *) last_result);
        }
    }

    last_result = null;
    last_result_type = cast(GretlType) 0;
}

/**
 * free_xtab:
 * @tab: pointer to gretl crosstab struct.
 *
 * Frees all resources associated with @tab, and the
 * pointer itself.
 */

void free_xtab(Xtab * tab)
{
    int i;

    if (tab == null) {
        return;
    }

    free(tab.rname);
    free(tab.cname);

    free(tab.rtotal);
    free(tab.ctotal);
    free(tab.rval);
    free(tab.cval);

    if (tab.f != null) {
        for (i=0; i<tab.rows; i++) {
            free(tab.f[i]);
        }
        free(tab.f);
    }

    if (tab.Sr != null) {
        strings_array_free(tab.Sr, tab.rows);
    }
    if (tab.Sc != null) {
        strings_array_free(tab.Sc, tab.cols);
    }

    free(tab);
}

private Xtab * xtab_new(int n, int t1, int t2)
{
    Xtab * tab = cast(Xtab *) malloc(Xtab.sizeof);

    if (tab == null) return null;

    tab.rtotal = null;
    tab.ctotal = null;
    tab.rval = null;
    tab.cval = null;
    tab.f = null;

    tab.n = n;
    tab.t1 = t1;
    tab.t2 = t2;
    tab.missing = 0;

    tab.rname = null;
    tab.cname = null;
    tab.Sr = null;
    tab.Sc = null;
    tab.rstrs = 0;
    tab.cstrs = 0;

    return tab;
}

/* allocate the required arrays (apart from S, which must be
   handled separately) and initialize counters to zero
*/

private int xtab_allocate_arrays(Xtab * tab)
{
    int i, j, err = 0;

    if (!tab.rstrs && tab.rval == null) {
        tab.rval = cast(double *) malloc(tab.rows * (*tab.rval).sizeof);
        if (tab.rval == null) {
            err = GretlError.alloc;
        }
    }

    if (!err && !tab.cstrs && tab.cval == null) {
        tab.cval = cast(double *) malloc(tab.cols * (*tab.cval).sizeof);
        if (tab.cval == null) {
            err = GretlError.alloc;
        }
    }

    if (!err) {
        tab.rtotal = cast(int *) malloc(tab.rows * (*tab.rtotal).sizeof);
        tab.ctotal = cast(int *) malloc(tab.cols * (*tab.ctotal).sizeof);
        if (tab.rtotal == null || tab.ctotal == null) {
            err = GretlError.alloc;
        } else {
            for (i=0; i<tab.rows; i++) {
                tab.rtotal[i] = 0;
            }
            for (j=0; j<tab.cols; j++) {
                tab.ctotal[j] = 0;
            }
        }
    }

    if (!err) {
        tab.f = cast(int **) malloc(tab.rows * (*tab.f).sizeof);
        if (tab.f == null) {
            err = GretlError.alloc;
        } else {
            for (i=0; i<tab.rows; i++) {
                tab.f[i] = null;
            }
        }
    }

    for (i=0; i<tab.rows && !err; i++) {
        tab.f[i] = cast(int *) malloc(tab.cols * (*tab.f[i]).sizeof);
        if (tab.f[i] == null) {
            err = GretlError.alloc;
        } else {
            for (j=0; j<tab.cols; j++) {
                tab.f[i][j] = 0;
            }
        }
    }

    return err;
}

/* also used in gretl_matrix.c */

extern (C) int compare_xtab_rows(const(void) * a, const(void) * b)
{
    const(double) ** da = cast(const(double) **) a;
    const(double) ** db = cast(const(double) **) b;
    double ret = da[0][0] - db[0][0];

    if (ret == 0) {
        ret = da[0][1] - db[0][1];
    }

    return ret < 0 ? -1 : ret > 0 ? 1: 0;
}

private int xtab_get_data(Xtab * tab, int v, int cols,
                          const(DatasetRaw) * dset,
                          SeriesTable ** pst)
{
    const(double) * x = dset.Z[v] + dset.t1;
    int n = sample_size(dset);
    GretlMatrixRaw * u;
    double ** xtarg;
    char *** Starg;
    int * itarg;
    int * ttarg;
    int err = 0;

    if (cols) {
        xtarg = &tab.cval;
        Starg = &tab.Sc;
        itarg = &tab.cols;
        ttarg = &tab.cstrs;
    } else {
        xtarg = &tab.rval;
        Starg = &tab.Sr;
        itarg = &tab.rows;
        ttarg = &tab.rstrs;
    }

    u = gretl_matrix_values(x, n, Gretlopt.OPT_S, &err);

    if (!err && is_string_valued(dset, v)) {
        SeriesTable * st;
        int ns, nv = u.rows;
        char ** S = null;
        char ** S0;

        *pst = st = series_get_string_table(dset, v);
        S0 = series_table_get_strings(st, &ns);

        if (ns > nv) {
            int i, k;

            S = strings_array_new(nv);
            if (S != null) {
                for (i=0; i<nv; i++) {
                    k = cast(int) (u.val[i] - 1);
                    if (k < 0 || k >= ns) {
                        /* the series' numeric codes are inconsistent
                           with its attached string table
                        */
                        err = GretlError.data;
                        break;
                    }
                    S[i] = gretl_strdup(S0[k]);
                }
            }
        } else {
            S = strings_array_dup(S0, ns);
        }
        if (!err && S == null) {
            err = GretlError.alloc;
        }
        if (!err) {
            *ttarg = 1;
            *itarg = nv;
            *Starg = S;
        }
    } else if (!err) {
        *itarg = u.rows;
        *xtarg = u.val;
        u.val = null;
    }

    gretl_matrix_free(u);

    return err;
}

private int xtab_row_match(Xtab * tab, int i, const(char) * s, double x)
{
    if (tab.rstrs) {
        return s != null && strcmp(s, tab.Sr[i]) == 0;
    } else {
        return x == tab.rval[i];
    }
}

private int xtab_col_match(Xtab * tab, int j, const(char) * s, double x)
{
    if (tab.cstrs) {
        return s != null && strcmp(s, tab.Sc[j]) == 0;
    } else {
        return x == tab.cval[j];
    }
}

private char * xtab_dim_string(int v, const(DatasetRaw) * dset,
                               Gretlopt opt)
{
    if (opt & Gretlopt.OPT_T) {
        return gretl_strdup(dset.varname[v]);
    } else {
        const(char) * s = series_get_display_name(dset, v);

        if (s != null && *s != '\0') {
            return gretl_strdup(s);
        } else {
            return gretl_strdup(dset.varname[v]);
        }
    }
}

private bool complete_obs(const(double) * x, const(double) * y, int t)
{
    return !na(x[t]) && !na(y[t]);
}

/* crosstab struct creation functions */

private Xtab * get_new_xtab(int rv, int cv, const(DatasetRaw) * dset,
                            Gretlopt opt, int * err)
{
    SeriesTable * sti = null;
    SeriesTable * stj = null;
    const(char) * s1 = null;
    const(char) * s2 = null;
    double x1 = 0, x2 = 0;
    Xtab * tab = null;
    int imatch, jmatch;
    int i, j, t, n = 0;

    for (t=dset.t1; t<=dset.t2; t++) {
        if (complete_obs(dset.Z[rv], dset.Z[cv], t)) {
            n++;
        }
    }

    if (n == 0) {
        *err = GretlError.missdata;
    } else {
        tab = xtab_new(n, dset.t1, dset.t2);
        if (tab == null) {
            *err = GretlError.alloc;
        }
    }

    if (!*err) {
        /* assemble row data */
        *err = xtab_get_data(tab, rv, 0, dset, &sti);
    }
    if (!*err) {
        /* assemble column data */
        *err = xtab_get_data(tab, cv, 1, dset, &stj);
    }
    if (!*err) {
        *err = xtab_allocate_arrays(tab);
    }

    if (*err) goto bailout;

    tab.missing = (dset.t2 - dset.t1 + 1) - n;
    tab.rname = xtab_dim_string(rv, dset, opt);
    tab.cname = xtab_dim_string(cv, dset, opt);

    /* The following could be made more efficient by substituting
       sorted arrays for dset->Z[rv] and dset->Z[cv] but it's not
       clear that the fixed cost would be recouped.
    */

    for (t=dset.t1; t<=dset.t2; t++) {
        if (!complete_obs(dset.Z[rv], dset.Z[cv], t)) {
            continue;
        }
        x1 = dset.Z[rv][t];
        if (tab.rstrs) {
            s1 = series_table_get_string(sti, x1);
        }
        x2 = dset.Z[cv][t];
        if (tab.cstrs) {
            s2 = series_table_get_string(stj, x2);
        }
        for (i=0; i<tab.rows; i++) {
            imatch = xtab_row_match(tab, i, s1, x1);
            if (imatch) {
                jmatch = 0;
                for (j=0; j<tab.cols && !jmatch; j++) {
                    jmatch = xtab_col_match(tab, j, s2, x2);
                    if (jmatch) {
                        tab.f[i][j] += 1;
                        tab.rtotal[i] += 1;
                        tab.ctotal[j] += 1;
                    }
                }
                break;
            }
        }
    }

 bailout:

    if (*err) {
        free_xtab(tab);
        tab = null;
    }

    return tab;
}

private void record_xtab(const(Xtab) * tab, const(DatasetRaw) * dset,
                         Gretlopt opt)
{
    GretlMatrixRaw * X = null;
    char ** Sc = null;
    char ** Sr = null;
    double xij, cj, ri;
    int totals, rows, cols;
    int i, j;

    totals = (opt & Gretlopt.OPT_N)? 0 : 1;
    rows = tab.rows + totals;
    cols = tab.cols + totals;

    X = gretl_zero_matrix_new(rows, cols);
    if (X == null) {
        return;
    }

    /* column labels */
    Sc = strings_array_new(cols);
    if (Sc != null) {
        for (j=0; j<tab.cols; j++) {
            if (tab.cstrs) {
                Sc[j] = gretl_strdup(tab.Sc[j]);
            } else {
                cj = tab.cval[j];
                Sc[j] = gretl_strdup_printf("%4g", cj);
            }
        }
        if (totals) {
            Sc[cols-1] = gretl_strdup("TOTAL");
        }
    }

    /* row labels */
    Sr = strings_array_new(rows);
    if (Sr != null) {
        for (i=0; i<tab.rows; i++) {
            if (tab.rstrs) {
                Sr[i] = gretl_strdup(tab.Sr[i]);
            } else {
                ri = tab.rval[i];
                Sr[i] = gretl_strdup_printf("%4g", ri);
            }
        }
        if (totals) {
            Sr[rows-1] = gretl_strdup("TOTAL");
        }
    }

    /* body of table */

    for (i=0; i<tab.rows; i++) {
        if (tab.rtotal[i] > 0) {
            /* row counts */
            for (j=0; j<tab.cols; j++) {
                if (tab.ctotal[j] > 0) {
                    if (opt & (Gretlopt.OPT_C | Gretlopt.OPT_R)) {
                        if (opt & Gretlopt.OPT_C) {
                            xij = 100.0 * tab.f[i][j] / tab.ctotal[j];
                        } else {
                            xij = 100.0 * tab.f[i][j] / tab.rtotal[i];
                        }
                    } else {
                        xij = tab.f[i][j];
                    }
                    gretl_matrix_set(X, i, j, xij);
                }
            }
            if (totals) {
                /* row totals */
                if (opt & Gretlopt.OPT_C) {
                    xij = 100.0 * tab.rtotal[i] / tab.n;
                } else {
                    xij = tab.rtotal[i];
                }
                gretl_matrix_set(X, i, cols-1, xij);
            }
        }
    }

    if (totals) {
        /* column totals */
        for (j=0; j<tab.cols; j++) {
            if (opt & Gretlopt.OPT_R) {
                xij = 100.0 * tab.ctotal[j] / tab.n;
            } else {
                xij = tab.ctotal[j];
            }
            gretl_matrix_set(X, rows-1, j, xij);
        }
        gretl_matrix_set(X, rows-1, cols-1, tab.n);
    }

    gretl_matrix_set_colnames(X, Sc);
    gretl_matrix_set_rownames(X, Sr);
    set_last_result_data(X, GretlType.matrix);
}

/* For use in the context of "xtab" with --quiet option:
   compute and record the Pearson chi-square value and its
   p-value.
*/

private int xtab_do_pearson(const(Xtab) * tab)
{
    double x, y, ymin = 1.0e-7;
    double pearson = 0.0;
    int i, j, err = 0;

    for (i=0; i<tab.rows && !err; i++) {
        if (tab.rtotal[i] > 0) {
            for (j=0; j<tab.cols && !err; j++) {
                y = (cast(double) tab.rtotal[i] * tab.ctotal[j]) / tab.n;
                if (y < ymin) {
                    err = GretlError.data;
                } else {
                    x = cast(double) tab.f[i][j] - y;
                    pearson += x * x / y;
                }
            }
        }
    }

    if (!err) {
        int df = (tab.rows - 1) * (tab.cols - 1);
        double pval = chisq_cdf_comp(df, pearson);

        if (!na(pval)) {
            record_test_result(pearson, pval);
        } else {
            err = GretlError.data;
        }
    }

    if (err) {
        record_test_result(NADBL, NADBL);
    }

    return err;
}

int crosstab_from_matrix(Gretlopt opt, Prn * prn)
{
    const(char) * mname;
    GretlMatrixRaw * m;
    Xtab * tab = null;
    int i, j, nvals, n = 0;
    int free_m = 0;
    double x;
    int err = 0;

    mname = get_optval_string(GretlCmdIndex.xtab, Gretlopt.OPT_X);
    if (mname == null) {
        return GretlError.data;
    }

    m = get_matrix_by_name(mname);
    if (m == null) {
        m = generate_matrix(mname, null, &err);
        if (m == null) {
            return GretlError.unkvar;
        } else {
            free_m = 1;
        }
    }

    if (m.rows < 2 || m.cols < 2) {
        err = GretlError.data;
    }

    nvals = m.rows * m.cols;

    for (i=0; i<nvals && !err; i++) {
        x = m.val[i];
        if (x < 0 || x != floor(x) || x > INT_MAX) {
            err = GretlError.data;
        }
        n += cast(int) x;
    }

    if (err) {
        gretl_errmsg_sprintf(_("Matrix %s does not represent a "
                               ~ "contingency table").ptr, mname);
    } else {
        tab = xtab_new(n, 0, 0);
        if (tab == null) {
            err = GretlError.alloc;
        }
    }

    if (err) {
        goto bailout;
    }

    tab.rows = m.rows;
    tab.cols = m.cols;
    tab.Sc = cast(char **) gretl_matrix_get_colnames(m);
    tab.Sr = cast(char **) gretl_matrix_get_rownames(m);

    if (xtab_allocate_arrays(tab)) {
        err = GretlError.alloc;
        goto bailout;
    }

    for (i=0; i<m.rows; i++) {
        tab.rval[i] = i + 1;
        tab.rtotal[i] = 0;
        for (j=0; j<m.cols; j++) {
            tab.f[i][j] = cast(int) gretl_matrix_get(m, i, j);
            tab.rtotal[i] += tab.f[i][j];
        }
    }

    for (j=0; j<m.cols; j++) {
        tab.cval[j] = j + 1;
        tab.ctotal[j] = 0;
        for (i=0; i<m.rows; i++) {
            tab.ctotal[j] += tab.f[i][j];
        }
    }

    if (opt & Gretlopt.OPT_Q) {
        xtab_do_pearson(tab);
    } else {
        print_xtab(tab, null, opt | Gretlopt.OPT_S, prn);
    }

 bailout:

    if (tab != null) {
        tab.Sc = null;
        tab.Sr = null;
        free_xtab(tab);
    }
    if (free_m) {
        gretl_matrix_free(m);
    }

    return err;
}

int crosstab(const(int) * list, const(DatasetRaw) * dset,
             Gretlopt opt, Prn * prn)
{
    Xtab * tab;
    int * rowvar = null;
    int * colvar = null;
    int i, j, vi, vj, k;
    int pos, onelist;
    int nrv, ncv;
    int err = 0;

    pos = gretl_list_separator_position(list);
    onelist = (pos == 0);

    if (pos == 0) {
        /* single list case */
        nrv = list[0];
        ncv = nrv - 1;
    } else {
        /* double list case */
        nrv = pos - 1;
        ncv = list[0] - pos;
    }

    if (nrv == 0 || ncv == 0) {
        return GretlError.parse;
    }

    rowvar = gretl_list_new(nrv);
    if (rowvar == null) {
        return GretlError.alloc;
    }

    j = 1;
    for (i=1; i<=nrv; i++) {
        k = list[i];
        if (accept_as_discrete(dset, k, 0)) {
            rowvar[j++] = k;
        } else {
            pprintf(prn, _("dropping %s: not a discrete variable\n").ptr,
                    dset.varname[k]);
            rowvar[0] -= 1;
        }
    }

    if (rowvar[0] == 0 || (onelist && rowvar[0] == 1)) {
        gretl_errmsg_set("xtab: variables must be discrete");
        free(rowvar);
        return GretlError.types;
    }

    if (onelist && rowvar[0] == 2) {
        /* the bivariate case */
        tab = get_new_xtab(rowvar[1], rowvar[2], dset, opt, &err);
        if (!err) {
            /* make $result matrix available */
            record_xtab(tab, dset, opt);
            if (opt & Gretlopt.OPT_Q) {
                /* quiet: run and record the Pearson test */
                xtab_do_pearson(tab);
            } else {
                /* print, and record Pearson test */
                print_xtab(tab, dset, opt | Gretlopt.OPT_S, prn);
            }
            free_xtab(tab);
        }
        goto finish;
    }

    if (!onelist) {
        /* construct the second list */
        colvar = gretl_list_new(ncv);
        if (colvar == null) {
            err = GretlError.alloc;
        } else {
            j = 1;
            for (i=1; i<=ncv; i++) {
                k = list[pos+i];
                if (accept_as_discrete(dset, k, 0)) {
                    colvar[j++] = k;
                } else {
                    colvar[0] -= 1;
                }
            }
            if (colvar[0] == 0) {
                err = GretlError.types;
            }
        }
    }

    for (i=1; i<=rowvar[0] && !err; i++) {
        vi = rowvar[i];
        if (onelist) {
            /* single list case */
            for (j=1; j<i && !err; j++) {
                vj = rowvar[j];
                tab = get_new_xtab(vj, vi, dset, opt, &err);
                if (!err) {
                    print_xtab(tab, dset, opt, prn);
                    free_xtab(tab);
                }
            }
        } else {
            /* double list case */
            for (j=1; j<=colvar[0] && !err; j++) {
                vj = colvar[j];
                tab = get_new_xtab(vi, vj, dset, opt, &err);
                if (!err) {
                    print_xtab(tab, dset, opt, prn);
                    free_xtab(tab);
                }
            }
        }
    }

 finish:

    free(rowvar);
    free(colvar);

    return err;
}

Xtab * single_crosstab(const(int) * list, const(DatasetRaw) * dset,
                       Gretlopt opt, Prn * prn, int * err)
{
    Xtab * tab = null;
    int rv, cv;

    if (list[0] != 2) {
        *err = GretlError.data;
        return null;
    }

    rv = list[1];
    cv = list[2];

    if (accept_as_discrete(dset, rv, 0) &&
        accept_as_discrete(dset, cv, 0)) {
        tab = get_new_xtab(rv, cv, dset, opt, err);
    } else {
        *err = GretlError.types;
    }

    if (!*err) {
        print_xtab(tab, dset, opt, prn);
    }

    return tab;
}

private int just_record_freq_test(const(FreqDist) * freq)
{
    double pval = NADBL;

    if (freq.dist == DistCode.normal) {
        pval = chisq_cdf_comp(2, freq.test);
    } else if (freq.dist == DistCode.gamma) {
        pval = normal_pvalue_2(freq.test);
    }

    if (na(pval)) {
        return GretlError.nan;
    } else {
        record_test_result(freq.test, pval);
        return 0;
    }
}

int model_error_dist(const(Model) * pmod, DatasetRaw * dset,
                     Gretlopt opt, Prn * prn)
{
    FreqDist * freq = null;
    Gretlopt fopt = Gretlopt.OPT_Z; /* show normality test */
    int save_t1 = dset.t1;
    int save_t2 = dset.t2;
    int err = 0;

    if (pmod == null || pmod.uhat == null) {
        return GretlError.data;
    }

    err = gretl_model_get_normality_test(pmod, prn);

    if (!err) {
        return 0;
    } else if (LIMDEP(cast(GretlCmdIndex) pmod.ci)) {
        return err;
    } else {
        err = 0;
    }

    if (exact_fit_check(pmod, prn)) {
        return 0;
    }

    if (genr_fit_resid(pmod, dset, ModelDataIndex.uhat)) {
        return GretlError.alloc;
    }

    if (!err) {
        dset.t1 = pmod.t1;
        dset.t2 = pmod.t2;
        freq = get_freq(dset.v - 1, dset, NADBL, NADBL, 0,
                        pmod.ncoeff, fopt, &err);
    }

    if (!err) {
        if (opt & Gretlopt.OPT_I) {
            err = just_record_freq_test(freq);
        } else if (opt & Gretlopt.OPT_Q) {
            print_freq_test(freq, prn);
        } else {
            print_freq(freq, 0, null, prn);
        }
        free_freq(freq);
    }

    dataset_drop_last_variables(dset, 1);

    dset.t1 = save_t1;
    dset.t2 = save_t2;

    return err;
}

/* Compute PACF via Durbin-Levinson algorithm and write
   the values into @A following the ACF values already
   present.
*/

private int get_pacf(GretlMatrixRaw * A)
{
    int m = gretl_matrix_rows(A);
    GretlMatrixRaw * phi;
    double * acf = A.val;
    double * pacf = acf + m;
    double x, num, den;
    int i, j;

    phi = gretl_matrix_alloc(m, m);
    if (phi == null) {
        return GretlError.alloc;
    }

    gretl_matrix_set(A, 0, 1, acf[0]);
    gretl_matrix_set(phi, 0, 0, acf[0]);

    for (i=1; i<m; i++) {
        num = acf[i];
        for (j=0; j<i; j++) {
            num -= gretl_matrix_get(phi, i-1, j) * acf[i-j-1];
        }
        den = 1.0;
        for (j=0; j<i; j++) {
            den -= gretl_matrix_get(phi, i-1, j) * acf[j];
        }
        pacf[i] = num / den;
        gretl_matrix_set(phi, i, i, pacf[i]);
        for (j=0; j<i; j++) {
            x = gretl_matrix_get(phi, i-1, j);
            x -= pacf[i] * gretl_matrix_get(phi, i-1, i-j-1);
            gretl_matrix_set(phi, i, j, x);
        }
    }

    gretl_matrix_free(phi);

    return 0;
}

/* also used in GUI, library.c */

int auto_acf_order(int T)
{
    int p = cast(int) (10 * log10(cast(double) T));

    if (p > T / 5) {
        /* restrict to 20 percent of data (Tadeusz) */
        p = T / 5;
    }

    return p;
}

/**
 * gretl_acf:
 * @k: lag order.
 * @t1: starting observation.
 * @t2: ending observation.
 * @y: data series.
 * @ybar: mean of @y over range @t1 to @t2.
 *
 * Returns: the autocorrelation at lag @k for the series @y over
 * the range @t1 to @t2, or #NADBL on failure.
 */

private double gretl_acf(int k, int t1, int t2, const(double) * y,
                         double ybar)
{
    double z, num = 0, den = 0;
    int t, n = 0;

    for (t=t1; t<=t2; t++) {
        if (na(y[t])) {
            return NADBL;
        }
        z = y[t] - ybar;
        den += z * z;
        if (t - k >= t1) {
            if (na(y[t-k])) {
                return NADBL;
            }
            num += z * (y[t-k] - ybar);
            n++;
        }
    }

    if (n == 0) {
        return NADBL;
    }

    return num / den;
}

/**
 * ljung_box:
 * @m: maximum lag.
 * @t1: starting observation.
 * @t2: ending observation.
 * @y: data series.
 * @err: location to receive error code.
 *
 * Returns: the Ljung-Box statistic for lag order @m for
 * the series @y over the sample @t1 to @t2, or #NADBL
 * on failure.
 */

double ljung_box(int m, int t1, int t2, const(double) * y, int * err)
{
    double acf, ybar = 0.0, LB = 0.0;
    int k, n = t2 - t1 + 1;

    *err = 0;

    if (n == 0 || gretl_isconst(t1, t2, y)) {
        *err = GretlError.data;
    } else if (m <= 0) {
        gretl_errmsg_sprintf(_("Invalid lag order %d").ptr, m);
        *err = GretlError.data;
    } else {
        ybar = gretl_mean(t1, t2, y);
        if (na(ybar)) {
            *err = GretlError.data;
        }
    }

    if (*err) {
        return NADBL;
    }

    /* calculate acf up to lag m, cumulating LB */
    for (k=1; k<=m; k++) {
        acf = gretl_acf(k, t1, t2, y, ybar);
        if (na(acf)) {
            *err = GretlError.missdata;
            break;
        }
        LB += acf * acf / (n - k);
    }

    if (*err) {
        LB = NADBL;
    } else {
        LB *= n * (n + 2.0);
    }

    return LB;
}

/**
 * gretl_xcf:
 * @k: lag order (or lead order if < 0).
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: first data series.
 * @y: second data series.
 * @xbar: mean of first series.
 * @ybar: mean of second series.
 *
 * Returns: the cross-correlation at lag (or lead) @k for the
 * series @x and @y over the range @t1 to @t2, or #NADBL on failure.
 */

private double gretl_xcf(int k, int t1, int t2, const(double) * x, const(double) * y,
                         double xbar, double ybar)
{
    double num = 0, den1 = 0, den2 = 0;
    double zx, zy;
    int t;

    for (t=t1; t<=t2; t++) {
        if (na(x[t]) || na(y[t])) {
            return NADBL;
        }
        zx = x[t] - xbar;
        zy = y[t] - ybar;
        den1 += zx * zx;
        den2 += zy * zy;
        if ((k >= 0 && t - k >= t1) || (k < 0 && t - k <= t2)) {
            num += zx * (y[t-k] - ybar);
        }
    }

    return num / sqrt(den1 * den2);
}

/* Full Bartlett calculation for three z-values, called from
   corrgm_plot() if we're printing ACF output. BPM will have three
   columns.
*/

private void do_acf_bartlett(GretlMatrixRaw * BPM, int i,
                             int T, double ssr,
                             const(double) * z)
{
    double c = sqrt((1.0/T) * (1 + 2*ssr));
    int j;

    for (j=0; j<BPM.cols; j++) {
        gretl_matrix_set(BPM, i, j, z[j] * c);
    }
}

/* Basic Bartlett if we're not printing ACF output. BPM will have a
   single column.
*/

private void corrgm_basic_bartlett(const(double) * acf,
                                   GretlMatrixRaw * BPM,
                                   int T)
{
    double Tm1 = 1.0 / T;
    double z = 1.96;
    double ssr = 0;
    int k;

    for (k=0; k<BPM.rows; k++) {
        if (na(acf[k])) {
            continue;
        }
        if (k > 0 && !na(acf[k-1])) {
            ssr += acf[k-1] * acf[k-1];
        }
        BPM.val[k] = z * sqrt(Tm1 * (1 + 2*ssr));
    }
}

private void corrgm_do_ljung_box(const(double) * acf,
                                 int m, int T,
                                 int nparam)
{
    double lbox = 0.0;
    int k, dfQ = 0;

    for (k=0; k<m; k++) {
        if (!na(acf[k])) {
            lbox += (T * (T + 2.0)) * acf[k] * acf[k] / (T - (k + 1));
            dfQ++;
        }
    }

    dfQ -= nparam;
    if (lbox > 0 && dfQ > 0) {
        double pval = chisq_cdf_comp(dfQ, lbox);

        if (!na(pval)) {
            record_test_result(lbox, pval);
        }
    }
}

/* The following is called if the --silent option is not given to
   "corrgm". Note that it embeds Bartlett calculations for three
   confidence levels. We make do with a simpler Bartlett variant if
   we're not printing the ACF table.
*/

private void corrgm_print(const(char) * vname,
                          const(double) * acf,
                          const(double) * pacf,
                          GretlMatrixRaw * BPM,
                          int m, int T,
                          int nparam,
                          Gretlopt opt,
                          Prn * prn)
{
    immutable double[3] z = [1.65, 1.96, 2.58];
    double[3] pm = 0;
    double pval = NADBL;
    double lbox = 0.0;
    double ssr = 0.0;
    int i, k, dfQ = 1;

    if (opt & Gretlopt.OPT_R) {
        pprintf(prn, "\n%s\n", _("Residual autocorrelation function").ptr);
    } else {
        pputc(prn, '\n');
        pprintf(prn, _("Autocorrelation function for %s").ptr, vname);
        pputc(prn, '\n');
    }

    pputs(prn, _("***, **, * indicate significance at the 1%, 5%, 10% levels\n").ptr);

    if (opt & Gretlopt.OPT_B) {
        pputs(prn, _("using Bartlett standard errors for ACF").ptr);
    } else {
        pprintf(prn, _("using standard error 1/T^%.1f").ptr, 0.5);
    }

    pputs(prn, "\n\n");
    if (pacf != null) {
        pputs(prn, _("  LAG      ACF          PACF         Q-stat. [p-value]").ptr);
    } else {
        pputs(prn, _("  LAG      ACF          Q-stat. [p-value]").ptr);
    }
    pputs(prn, "\n\n");

    if (BPM == null || pacf != null) {
        /* if we're not doing Bartlett confidence bands, or we're
           showing PACF as well as ACF
        */
        for (i=0; i<3; i++) {
            pm[i] = z[i] / sqrt(cast(double) T);
            if (pm[i] > 0.5) {
                pm[i] = 0.5;
            }
        }
    }

    for (k=0; k<m; k++) {
        double pm0, pm1, pm2;

        if (na(acf[k])) {
            pprintf(prn, "%5d\n", k + 1);
            continue;
        }

        /* ACF */
        pprintf(prn, "%5d%9.4f ", k + 1, acf[k]);
        if (BPM != null) {
            if (k > 0 && !na(acf[k-1])) {
                ssr += acf[k-1] * acf[k-1];
            }
            do_acf_bartlett(BPM, k, T, ssr, z.ptr);
            pm0 = gretl_matrix_get(BPM, k, 0);
            pm1 = gretl_matrix_get(BPM, k, 1);
            pm2 = gretl_matrix_get(BPM, k, 2);
        } else {
            pm0 = pm[0];
            pm1 = pm[1];
            pm2 = pm[2];
        }
        if (fabs(acf[k]) > pm2) {
            pputs(prn, " ***");
        } else if (fabs(acf[k]) > pm1) {
            pputs(prn, " ** ");
        } else if (fabs(acf[k]) > pm0) {
            pputs(prn, " *  ");
        } else {
            pputs(prn, "    ");
        }

        /* PACF, if present */
        if (pacf != null) {
            if (na(pacf[k])) {
                bufspace(13, prn);
            } else {
                /* PACF */
                pprintf(prn, "%9.4f", pacf[k]);
                if (fabs(pacf[k]) > pm[2]) {
                    pputs(prn, " ***");
                } else if (fabs(pacf[k]) > pm[1]) {
                    pputs(prn, " ** ");
                } else if (fabs(pacf[k]) > pm[0]) {
                    pputs(prn, " *  ");
                } else {
                    pputs(prn, "    ");
                }
            }
        }

        /* Ljung-Box Q */
        lbox += (T * (T + 2.0)) * acf[k] * acf[k] / (T - (k + 1));

        if (k >= nparam) {
            /* i.e., if the real df is > 0 */
            pprintf(prn, "%12.4f", lbox);
            pval = chisq_cdf_comp(dfQ++, lbox);
            pprintf(prn, "  [%5.3f]", pval);
        }

        pputc(prn, '\n');
    }
    pputc(prn, '\n');

    if (lbox > 0 && !na(pval)) {
        record_test_result(lbox, pval);
    }
}

/**
 * corrgram:
 * @varno: ID number of variable to process.
 * @order: integer order for autocorrelation function.
 * @nparam: number of estimated parameters (e.g. for the
 * case of ARMA), used to correct the degrees of freedom
 * for Q test.
 * @dset: dataset struct.
 * @opt: if includes OPT_R, variable in question is a model
 * residual generated "on the fly"; OPT_U can be used to
 * specify a plot option; OPT_S says not to print anything;
 * OPT_A means ACF only (no PACF); OPT_B means Bartlett.
 * @prn: gretl printing struct.
 *
 * Computes the autocorrelation function and plots the correlogram for
 * the variable specified by @varno.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int corrgram(int varno, int order, int nparam, DatasetRaw * dset,
             Gretlopt opt, Prn * prn)
{
    GretlMatrixRaw * BPM = null;
    GretlMatrixRaw * A = null;
    double ybar;
    double * acf = null;
    double * pacf = null;
    const(char) * vname;
    int BPMcols = 3;
    int k, m, T;
    int do_plot = 0;
    int do_print = 1;
    int do_pacf = 1;
    int t1 = dset.t1;
    int t2 = dset.t2;
    int quiet = 0;
    int err = 0;

    if (!dataset_is_time_series(dset)) {
        return GretlError.tsonly;
    } else if (order < 0) {
        gretl_errmsg_sprintf(_("Invalid lag order %d").ptr, order);
        return GretlError.data;
    }

    if (opt & (Gretlopt.OPT_Q | Gretlopt.OPT_S)) {
        quiet = 1;
    }

    gretl_error_clear();

    err = series_adjust_sample(dset.Z[varno], &t1, &t2);
    if (err) {
        return err;
    }

    if ((T = t2 - t1 + 1) < 4) {
        return GretlError.toofew;
    }

    if (gretl_isconst(t1, t2, dset.Z[varno])) {
        gretl_errmsg_sprintf(_("%s is a constant").ptr, dset.varname[varno]);
        return GretlError.data;
    }

    ybar = gretl_mean(t1, t2, dset.Z[varno]);
    if (na(ybar)) {
        return GretlError.data;
    }

    if (quiet) {
        do_print = 0;
        BPMcols = 1;
    }
    if (opt & Gretlopt.OPT_A) {
        do_pacf = 0;
    }

    /* lag order for acf */
    m = order;
    if (m == 0) {
        m = auto_acf_order(T);
    } else if (m > T - dset.pd) {
        int mmax = T - 1;

        if (m > mmax) {
            m = mmax;
        }
    }

    /* allocate space for ACF, and perhaps PACF */
    A = gretl_matrix_alloc(m, do_pacf ? 2 : 1);
    if (A == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    if (opt & Gretlopt.OPT_B) {
        /* Bartlett */
        BPM = gretl_matrix_alloc(m, BPMcols);
        if (BPM == null) {
            err = GretlError.alloc;
            goto bailout;
        }
    }

    /* graphing? */
    do_plot = gnuplot_graph_wanted(PlotType.correlogram, opt, &err);
    if (err) {
        goto bailout;
    }

    /* convenience pointer into @A */
    acf = A.val;

    /* calculate acf up to order @m */
    for (k=0; k<m; k++) {
        acf[k] = gretl_acf(k+1, t1, t2, dset.Z[varno], ybar);
    }

    /* try adding pacf into @A if wanted */
    if (do_pacf) {
        int pacf_err = get_pacf(A);

        pacf = pacf_err ? null : acf + m;
    }

    vname = plotname(dset, varno, 1);

    if (do_print) {
        corrgm_print(vname, acf, pacf, BPM, m, T, nparam, opt, prn);
    } else {
        if (opt & Gretlopt.OPT_B) {
            corrgm_basic_bartlett(acf, BPM, T);
        }
        corrgm_do_ljung_box(acf, m, T, nparam);
    }

    if (do_plot) {
        double pm = 1.96 / sqrt(cast(double) T);

        if (pm > 0.5) pm = 0.5;
        err = correlogram_plot(vname, acf, pacf, BPM, m, pm, opt);
    }

 bailout:

    gretl_matrix_free(A);
    gretl_matrix_free(BPM);

    return err;
}

/**
 * acf_matrix:
 * @x: series to analyse.
 * @order: maximum lag for autocorrelation function.
 * @dset: information on the data set, or %NULL.
 * @n: length of series (required if @dset is %NULL).
 * @err: location to receive error code.
 *
 * Computes the autocorrelation function for series @x with
 * maximum lag @order.
 *
 * Returns: two-column matrix containing the values of the
 * ACF and PACF at the successive lags, or %NULL on error.
 */

GretlMatrixRaw * acf_matrix(const(double) * x, int order,
                         const(DatasetRaw) * dset, int n,
                         int * err)
{
    GretlMatrixRaw * A = null;
    double xbar;
    int m, k, t, T;
    int t1, t2;

    if (dset != null) {
        t1 = dset.t1;
        t2 = dset.t2;

        while (na(x[t1])) t1++;
        while (na(x[t2])) t2--;

        T = t2 - t1 + 1;
    } else {
        t1 = 0;
        t2 = n - 1;
        T = n;
    }

    if (T < 4) {
        *err = GretlError.toofew;
        return null;
    }

    for (t=t1; t<=t2; t++) {
        if (na(x[t])) {
            *err = GretlError.missdata;
            return null;
        }
    }

    if (gretl_isconst(t1, t2, x)) {
        gretl_errmsg_set(_("Argument is a constant").ptr);
        *err = GretlError.data;
        return null;
    }

    xbar = gretl_mean(t1, t2, x);
    if (na(xbar)) {
        *err = GretlError.data;
        return null;
    }

    m = order;

    if (dset == null) {
        if (m < 1 || m > T) {
            *err = GretlError.data;
            return null;
        }
    } else {
        if (m == 0) {
            m = auto_acf_order(T);
        } else if (m > T - dset.pd) {
            int mmax = T - 1;

            if (m > mmax) {
                m = mmax;
            }
        }
    }

    A = gretl_matrix_alloc(m, 2);
    if (A == null) {
        *err = GretlError.alloc;
        return null;
    }

    /* calculate ACF up to order m */
    for (k=0; k<m && !*err; k++) {
        A.val[k] = gretl_acf(k+1, t1, t2, x, xbar);
        if (na(A.val[k])) {
            *err = GretlError.data;
        }
    }

    /* add PACF */
    if (!*err) {
        *err = get_pacf(A);
    }

    if (*err) {
        gretl_matrix_free(A);
        A = null;
    }

    return A;
}

private int xcorrgm_graph(const(char) * xname, const(char) * yname,
                          double * xcf, int m, double * pm,
                          int allpos)
{
    char[16] crit_string = 0;
    gchar * title;
    FILE * fp;
    int k, err = 0;

    fp = open_plot_input_file(PlotType.xcorrelogram, cast(GptFlags) 0, &err);
    if (err) {
        return err;
    }

    sprintf(crit_string.ptr, "%.2f/T^%.1f", 1.96, 0.5);

    gretl_push_c_numeric_locale();

    fputs("set xzeroaxis\n", fp);
    fputs("set yzeroaxis\n", fp);
    print_keypos_string(GpKeySpec.right | GpKeySpec.top, fp);
    fprintf(fp, "set xlabel '%s'\n", _("lag").ptr);
    if (allpos) {
        fputs("set yrange [-0.1:1.1]\n", fp);
    } else {
        fputs("set yrange [-1.1:1.1]\n", fp);
    }
    title = g_strdup_printf(_("Correlations of %s and lagged %s").ptr,
                            xname, yname);
    fprintf(fp, "set title '%s'\n", title);
    g_free(title);
    fprintf(fp, "set xrange [%d:%d]\n", -(m + 1), m + 1);
    if (allpos) {
        fprintf(fp, "plot \\\n"
                ~ "'-' using 1:2 notitle w impulses lw 5, \\\n"
                ~ "%g title '%s' lt 2\n", pm[1], crit_string.ptr);
    } else {
        fprintf(fp, "plot \\\n"
                ~ "'-' using 1:2 notitle w impulses lw 5, \\\n"
                ~ "%g title '+- %s' lt 2, \\\n"
                ~ "%g notitle lt 2\n", pm[1], crit_string.ptr, -pm[1]);
    }

    for (k=-m; k<=m; k++) {
        fprintf(fp, "%d %g\n", k, xcf[k+m]);
    }
    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    return finalize_plot_input_file(fp);
}

/* We assume here that all data issues have already been assessed (lag
   length, missing values etc.) and we just get on with the job.
*/

private GretlMatrixRaw * real_xcf_vec(const(double) * x, const(double) * y,
                                   int p, int T, int * err)
{
    GretlMatrixRaw * xcf;
    double xbar, ybar;
    int i;

    xbar = gretl_mean(0, T-1, x);
    if (na(xbar)) {
        *err = GretlError.data;
        return null;
    }

    ybar = gretl_mean(0, T-1, y);
    if (na(ybar)) {
        *err = GretlError.data;
        return null;
    }

    xcf = gretl_column_vector_alloc(p * 2 + 1);
    if (xcf == null) {
        *err = GretlError.alloc;
        return null;
    }

    for (i=-p; i<=p; i++) {
        xcf.val[i+p] = gretl_xcf(i, 0, T - 1, x, y, xbar, ybar);
    }

    return xcf;
}

/* for arrays @x and @y, check that there are no missing values
   and that neither series has a constant value
*/

private int xcf_data_check(const(double) * x, const(double) * y, int T,
                           int * badvar)
{
    int xconst = 1, yconst = 1;
    int t;

    if (T < 5) {
        return GretlError.toofew;
    }

    for (t=0; t<T; t++) {
        if (na(x[t]) || na(y[t])) {
            return GretlError.missdata;
        }
        if (t > 0 && x[t] != x[0]) {
            xconst = 0;
        }
        if (t > 0 && y[t] != y[0]) {
            yconst = 0;
        }
    }

    if (xconst) {
        *badvar = 1;
        return GretlError.data;
    } else if (yconst) {
        *badvar = 2;
        return GretlError.data;
    }

    return 0;
}

/**
 * xcf_vec:
 * @x: first series.
 * @y: second series.
 * @p: maximum lag for cross-correlation function.
 * @dset: information on the data set, or NULL.
 * @n: length of series (required only if @dset is NULL).
 * @err: location to receive error code.
 *
 * Computes the cross-correlation function for series @x with
 * series @y up to maximum lag @order.
 *
 * Returns: column vector containing the values of the
 * cross-correlation function, or NULL on error.
 */

GretlMatrixRaw * xcf_vec(const(double) * x, const(double) * y,
                      int p, const(DatasetRaw) * dset,
                      int n, int * err)
{
    GretlMatrixRaw * xcf = null;
    int t1, t2;
    int T, badvar = 0;

    if (p <= 0) {
        *err = GretlError.data;
        return null;
    }

    if (dset != null) {
        int yt1, yt2;

        t1 = yt1 = dset.t1;
        t2 = yt2 = dset.t2;

        while (na(x[t1])) t1++;
        while (na(y[yt1])) yt1++;

        while (na(x[t2])) t2--;
        while (na(y[yt2])) yt2--;

        t1 = (yt1 > t1)? yt1 : t1;
        t2 = (yt2 < t2)? yt2 : t2;

        T = t2 - t1 + 1;
    } else {
        t1 = 0;
        t2 = n - 1;
        T = n;
    }

    if (dset != null) {
        if (2 * p > T - dset.pd) { /* ?? */
            *err = GretlError.data;
        }
    } else if (2 * p > T) {
        *err = GretlError.data;
    }

    if (!*err) {
        *err = xcf_data_check(x + t1, y + t1, T, &badvar);
        if (badvar) {
            gretl_errmsg_sprintf(_("Argument %d is a constant").ptr,
                                 badvar);
        }
    }

    if (!*err) {
        xcf = real_xcf_vec(x + t1, y + t1, p, T, err);
    }

    return xcf;
}

/**
 * xcorrgram:
 * @list: should contain ID numbers of two variables.
 * @order: integer order for autocorrelation function.
 * @dset: dataset struct.
 * @opt: may include OPT_U for plot options, OPT_Q to
 * skip printing results.
 * @prn: gretl printing struct.
 *
 * Computes the cross-correlation function and plots the
 * cross-correlogram for the specified variables.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int xcorrgram(const(int) * list, int order, DatasetRaw * dset,
              Gretlopt opt, Prn * prn)
{
    GretlMatrixRaw * xcf = null;
    double[3] pm;
    const(char) * xname;
    const(char) * yname;
    const(double) * x;
    const(double) * y;
    int t1 = dset.t1, t2 = dset.t2;
    int do_print = 1;
    int do_plot = 0;
    int k, p, badvar = 0;
    int T, err = 0;

    if (!dataset_is_time_series(dset)) {
        return GretlError.tsonly;
    } else if (order < 0) {
        gretl_errmsg_sprintf(_("Invalid lag order %d").ptr, order);
        return GretlError.data;
    }

    gretl_error_clear();

    if (list[0] != 2) {
        return GretlError.data;
    }

    x = dset.Z[list[1]];
    y = dset.Z[list[2]];

    err = list_adjust_sample(list, &t1, &t2, dset, null);

    if (!err) {
        T = t2 - t1 + 1;
        err = xcf_data_check(x + t1, y + t1, T, &badvar);
    }

    if (err) {
        if (badvar) {
            gretl_errmsg_sprintf(_("%s is a constant").ptr,
                                 dset.varname[list[badvar]]);
        }
        return err;
    }

    if (opt & (Gretlopt.OPT_Q | Gretlopt.OPT_S)) {
        /* suppress printing */
        do_print = 0;
    }

    /* graphing? */
    do_plot = gnuplot_graph_wanted(PlotType.xcorrelogram, opt, &err);
    if (err) {
        return err;
    }

    xname = dset.varname[list[1]];
    yname = dset.varname[list[2]];

    p = order;
    if (p == 0) {
        p = auto_acf_order(T) / 2;
    } else if (2 * p > T - dset.pd) {
        p = (T - 1) / 2; /* ?? */
    }

    xcf = real_xcf_vec(x + t1, y + t1, p, T, &err);
    if (err) {
        return err;
    }

    /* for confidence bands */
    pm[0] = 1.65 / sqrt(cast(double) T);
    pm[1] = 1.96 / sqrt(cast(double) T);
    pm[2] = 2.58 / sqrt(cast(double) T);

    if (do_print) {
        pputc(prn, '\n');
        pprintf(prn, _("Cross-correlation function for %s and %s").ptr,
                xname, yname);
        pputs(prn, "\n\n");
        pputs(prn, _("  LAG      XCF").ptr);
        pputs(prn, "\n\n");

        for (k=-p; k<=p; k++) {
            double xk = xcf.val[k + p];

            pprintf(prn, "%5d%9.4f", k, xk);
            if (fabs(xk) > pm[2]) {
                pputs(prn, " ***");
            } else if (fabs(xk) > pm[1]) {
                pputs(prn, " **");
            } else if (fabs(xk) > pm[0]) {
                pputs(prn, " *");
            }
            pputc(prn, '\n');
        }
        pputc(prn, '\n');
    }

    if (do_plot) {
        int allpos = 1;

        for (k=-p; k<=p; k++) {
            if (xcf.val[k+p] < 0) {
                allpos = 0;
                break;
            }
        }
        err = xcorrgm_graph(xname, yname, xcf.val, p, pm.ptr, allpos);
    }

    if (!err) {
        /* save a matrix accessible via $result */
        GretlMatrixRaw * res;
        char ** S = null;
        int r = 2*p + 1;

        res = gretl_matrix_alloc(r, 2);
        for (k=0; k<r; k++) {
            gretl_matrix_set(res, k, 0, k-p);
            gretl_matrix_set(res, k, 1, xcf.val[k]);
        }
        S = strings_array_new(2);
        S[0] = gretl_strdup("LAG");
        S[1] = gretl_strdup("XCF");
        gretl_matrix_set_colnames(res, S);
        set_last_result_data(res, GretlType.matrix);
    }

    gretl_matrix_free(xcf);

    return err;
}

private struct FractintTest {
    double d;    /* estimated degree of integration */
    double se;   /* standard error of the above */
}

private int fract_int_GPH(int T, int m, const(double) * dens,
                          FractintTest * ft)
{
    GretlMatrixRaw * y = null;
    GretlMatrixRaw * X = null;
    GretlMatrixRaw * b = null;
    GretlMatrixRaw * V = null;
    double x, w;
    int t, err = 0;

    ft.d = ft.se = NADBL;

    y = gretl_column_vector_alloc(m);
    X = gretl_unit_matrix_new(m, 2);
    b = gretl_column_vector_alloc(2);
    V = gretl_matrix_alloc(2, 2);

    if (y == null || X == null || b == null || V == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    /* Test from Geweke and Porter-Hudak, as set out in
       Greene, Econometric Analysis 4e, p. 787 */

    for (t=0; t<m; t++) {
        y.val[t] = log(dens[t]);
        w = M_2PI * (t + 1) / cast(double) T;
        x = sin(w / 2);
        gretl_matrix_set(X, t, 1, log(4 * x * x));
    }

    err = gretl_matrix_ols(y, X, b, V, null, &x);

    if (!err) {
        ft.d = -b.val[1];
        ft.se = sqrt(gretl_matrix_get(V, 1, 1));
    }

 bailout:

    gretl_matrix_free(y);
    gretl_matrix_free(X);
    gretl_matrix_free(b);
    gretl_matrix_free(V);

    return err;
}

private GretlMatrixRaw * gretl_matrix_pergm(const(GretlMatrixRaw) * x, int m,
                                         int * err)
{
    GretlMatrixRaw * p = null;
    GretlMatrixRaw * f = null;

    f = gretl_matrix_fft(x, err);
    if (*err) {
        return null;
    }

    p = gretl_column_vector_alloc(m);

    if (p == null) {
        *err = GretlError.alloc;
    } else {
        int T = gretl_vector_get_length(x);
        double re, im, scale = M_2PI * T;
        Cmplx z;
        int i;

        for (i=0; i<m; i++) {
            z = gretl_cmatrix_get(f, i+1, 0);
            re = z.r;
            im = z.i;
            p.val[i] = (re*re + im*im) / scale;
        }
    }

    gretl_matrix_free(f);

    return p;
}

private struct LweHelper {
    GretlMatrixRaw * lambda;
    GretlMatrixRaw * lpow;
    GretlMatrixRaw * I1;
    GretlMatrixRaw * I2;
    double lcm = 0;
}

private void LWE_free(LweHelper * L)
{
    gretl_matrix_free(L.lambda);
    gretl_matrix_free(L.lpow);
    gretl_matrix_free(L.I1);
    gretl_matrix_free(L.I2);
}

private double LWE_obj_func(LweHelper * L, double d)
{
    double dd = 2.0 * d;
    int i;

    gretl_matrix_copy_values(L.lpow, L.lambda);
    gretl_matrix_raise(L.lpow, dd);

    for (i=0; i<L.I1.rows; i++) {
        L.I2.val[i] = L.I1.val[i] * L.lpow.val[i];
    }

    return -(log(gretl_vector_mean(L.I2)) - dd * L.lcm);
}

private GretlMatrixRaw * LWE_lambda(const(GretlMatrixRaw) * I1, int n)
{
    GretlMatrixRaw * lambda;
    int i, m = gretl_vector_get_length(I1);

    lambda = gretl_column_vector_alloc(m);

    if (lambda != null) {
        for (i=0; i<m; i++) {
            gretl_vector_set(lambda, i, (M_2PI / n) * (i + 1));
        }
    }

    return lambda;
}

private int LWE_init(LweHelper * L, const(GretlMatrixRaw) * X, int m)
{
    int i, err = 0;

    L.I2 = L.lpow = null;

    L.I1 = gretl_matrix_pergm(X, m, &err);
    if (err) {
        return err;
    }

    L.lambda = LWE_lambda(L.I1, X.rows);
    if (L.lambda == null) {
        gretl_matrix_free(L.I1);
        return GretlError.alloc;
    }

    L.lpow = gretl_matrix_copy(L.lambda);
    L.I2 = gretl_matrix_copy(L.I1);

    if (L.lpow == null || L.I2 == null) {
        err = GretlError.alloc;
        LWE_free(L);
    } else {
        L.lcm = 0.0;
        for (i=0; i<m; i++) {
            L.lcm += log(L.lambda.val[i]);
        }
        L.lcm /= m;
    }

    return err;
}

private enum LWE_MAXITER = 100;

private double LWE_calc(const(GretlMatrixRaw) * X, int m, int * err)
{
    LweHelper L;
    double d = 0, dd = 1.0;
    double eps = 1.0e-05;
    double f, incr, incl, deriv, h;
    int iter = 0;

    *err = LWE_init(&L, X, m);
    if (*err) {
        return NADBL;
    }

    while (fabs(dd) > 1.0e-06 && iter++ < LWE_MAXITER) {
        f = LWE_obj_func(&L, d);
        incr = LWE_obj_func(&L, d + eps) / eps;
        incl = LWE_obj_func(&L, d - eps) / eps;

        deriv = (incr - incl) / 2.0;
        h = (0.5 * (incr + incl) - f / eps) / eps;
        dd = (h < 0)? (-deriv / h) : deriv;

        if (fabs(dd) > 1) {
            dd = (dd > 0)? 1 : -1;
        }

        d += 0.5 * dd;
    }

    if (*err) {
        d = NADBL;
    } else if (iter == LWE_MAXITER) {
        fprintf(stderr, "LWE: max iterations reached\n");
        d = NADBL;
    }

    LWE_free(&L);

    return d;
}

int auto_spectrum_order(int T, Gretlopt opt)
{
    int m;

    if (opt & Gretlopt.OPT_O) {
        /* Bartlett */
        m = cast(int) (2.0 * sqrt(cast(double) T));
    } else {
        /* fractional integration test */
        double m1 = floor(cast(double) T / 2.0);
        double m2 = floor(pow(cast(double) T, 0.6));

        m = cast(int) ((m1 < m2)? m1 : m2);
        m--;
    }

    return m;
}

private int fract_int_LWE(const(double) * x, int m, int t1, int t2,
                          FractintTest * ft)
{
    GretlMatrixRaw * X;
    int err = 0;

    X = gretl_vector_from_series(x, t1, t2);

    if (X == null) {
        err = GretlError.alloc;
    } else {
        ft.d = LWE_calc(X, m, &err);
        if (!err) {
            ft.se = 1.0 / (2 * sqrt(cast(double) m));
        }
        gretl_matrix_free(X);
    }

    return err;
}

/* called by pergm_or_fractint(), when responding to the
   "pergm" command only, without its --silent option
*/

private void pergm_print(const(char) * vname, const(double) * d,
                         int T, int L, Gretlopt opt, Prn * prn)
{
    char[32] xstr = 0;
    double dt, yt;
    int t;

    if (vname == null) {
        pprintf(prn, "\n%s\n", _("Residual periodogram").ptr);
    } else {
        pprintf(prn, _("\nPeriodogram for %s\n").ptr, vname);
    }

    pprintf(prn, _("Number of observations = %d\n").ptr, T);
    if (opt & Gretlopt.OPT_O) {
        pprintf(prn, _("Using Bartlett lag window, length %d\n\n").ptr, L);
    } else {
        pputc(prn, '\n');
    }

    if (opt & Gretlopt.OPT_L) {
        pputs(prn, _(" omega  scaled frequency  periods  log spectral density\n\n").ptr);
    } else {
        pputs(prn, _(" omega  scaled frequency  periods  spectral density\n\n").ptr);
    }

    for (t=1; t<=T/2; t++) {
        dt = (opt & Gretlopt.OPT_L)? log(d[t-1]) : d[t-1];
        if (opt & Gretlopt.OPT_R) {
            yt = 2 * M_PI * cast(double) t / T;
            pprintf(prn, " %.5f%8d%16.2f", yt, t, cast(double) T / t);
        } else if (opt & Gretlopt.OPT_D) {
            yt = 360 * t / cast(double) T;
            pprintf(prn, " %7.3f%8d%16.2f", yt, t, cast(double) T / t);
        } else {
            yt = M_2PI * t / cast(double) T;
            pprintf(prn, " %.5f%8d%16.2f", yt, t, cast(double) T / t);
        }
        dt = (opt & Gretlopt.OPT_L)? log(d[t-1]) : d[t-1];
        sprintf(xstr.ptr, "%#.5g", dt);
        gretl_fix_exponent(xstr.ptr);
        pprintf(prn, "%16s\n", xstr.ptr);
    }

    pputc(prn, '\n');
}

private int finalize_fractint(const(double) * x,
                              const(double) * dens,
                              int t1, int t2, int width,
                              const(char) * vname,
                              Gretlopt opt,
                              Prn * prn)
{
    FractintTest ft;
    GretlMatrixRaw * result = null;
    int do_GPH = (opt & (Gretlopt.OPT_G | Gretlopt.OPT_A)) != 0;
    int do_LWE = !(opt & Gretlopt.OPT_G);
    int do_print = !(opt & Gretlopt.OPT_Q);
    int T = t2 - t1 + 1;
    int m, err = 0;

    /* order for test */
    if (width <= 0) {
        m = auto_spectrum_order(T, Gretlopt.OPT_NONE);
    } else {
        m = (width > T / 2)? T / 2 : width;
    }

    if (do_print && vname != null) {
        pprintf(prn, "\n%s, T = %d\n\n", vname, T);
    }

    if (do_GPH && do_LWE) {
        result = gretl_matrix_alloc(2, 2);
    } else {
        result = gretl_matrix_alloc(1, 2);
    }

    if (do_LWE) {
        /* do Local Whittle if wanted */
        err = fract_int_LWE(x, m, t1, t2, &ft);
        if (!err) {
            double z = ft.d / ft.se;
            double pv = normal_pvalue_2(z);

            record_test_result(z, pv);
            gretl_matrix_set(result, 0, 0, ft.d);
            gretl_matrix_set(result, 0, 1, ft.se);

            if (do_print) {
                pprintf(prn, "%s (m = %d)\n"
                        ~ "  %s = %g (%g)\n"
                        ~ "  %s: z = %g, %s %.4f\n\n",
                        _("Local Whittle Estimator").ptr, m,
                        _("Estimated degree of integration").ptr, ft.d, ft.se,
                        _("test statistic").ptr, z,
                        _("with p-value").ptr, pv);
            }
        }
    }

    if (!err && do_GPH) {
        /* do GPH if wanted (options --all or --gph) */
        int row = do_LWE ? 1 : 0;

        err = fract_int_GPH(T, m, dens, &ft);
        gretl_matrix_set(result, row, 0, ft.d);
        gretl_matrix_set(result, row, 1, ft.se);

        if (!err && (!do_LWE || do_print)) {
            double tval = ft.d / ft.se;
            int df = m - 2;
            double pv = student_pvalue_2(df, tval);

            if (!do_LWE) {
                record_test_result(tval, pv);
            }
            if (do_print) {
                pprintf(prn, "%s (m = %d)\n"
                        ~ "  %s = %g (%g)\n"
                        ~ "  %s: t(%d) = %g, %s %.4f\n\n",
                        _("GPH test").ptr, m,
                        _("Estimated degree of integration").ptr, ft.d, ft.se,
                        _("test statistic").ptr, df, tval,
                        _("with p-value").ptr, pv);
            }
        }
    }

    if (err) {
        gretl_matrix_free(result);
    } else {
        char ** colnames = strings_array_new(2);

        colnames[0] = gretl_strdup("d");
        colnames[1] = gretl_strdup("se");
        gretl_matrix_set_colnames(result, colnames);
        set_last_result_data(result, GretlType.matrix);
    }

    return err;
}

private double * pergm_bartlett_density(const(double) * x, int t1, int t2,
                                        int L, int * err)
{
    double * sdy;
    double * acov;
    double * dens;
    double xx, yy, vx, sx, w;
    int k, t, T = t2 - t1 + 1;

    sdy = cast(double *) malloc(T * (*sdy).sizeof);
    acov = cast(double *) malloc((L + 1) * (*acov).sizeof);
    dens = cast(double *) malloc(T/2 * (*dens).sizeof);

    if (sdy == null || acov == null || dens == null) {
        *err = GretlError.alloc;
        free(sdy);
        free(acov);
        free(dens);
        return null;
    }

    xx = gretl_mean(t1, t2, x);
    vx = real_gretl_variance(t1, t2, x, 1);
    sx = sqrt(vx);

    for (t=t1; t<=t2; t++) {
        sdy[t-t1] = (x[t] - xx) / sx;
    }

    /* autocovariances */
    for (k=1; k<=L; k++) {
        acov[k] = 0.0;
        for (t=k; t<T; t++) {
            acov[k] += sdy[t] * sdy[t-k];
        }
        acov[k] /= T;
    }

    vx /= M_2PI;

    for (t=0; t<T/2; t++) {
        yy = M_2PI * (t+1) / cast(double) T;
        xx = 1.0;
        for (k=1; k<=L; k++) {
            w = 1.0 - cast(double) k/(L + 1);
            xx += 2.0 * w * acov[k] * cos(yy * k);
        }
        dens[t] = xx * vx;
    }

    free(sdy);
    free(acov);

    return dens;
}

private enum PergmUsage {
    pergm_cmd,
    pergm_func,
    fractint_cmd
}

/* The following is a little complicated because we're supporting
   three use cases:

   1) Implementing the "pergm" command, with or without the
   option of using a Bartlett window (OPT_O); in this case the
   final matrix-location argument will be NULL.

   2) Implementing the user-space pergm function: in this case
   the final matrix-location argument will be non-NULL, and the
   @vname argument will be NULL.

   3) Implementing the fractint command, computing the Local
   Whittle Estimator and/or the GPH test -- this is flagged by
   opt & OPT_F. If we're doing Whittle only, we don't need to compute
   the spectral density via the autocovariances approach; that is
   handled via FFT in the LWE functions above.
*/

private int pergm_or_fractint(PergmUsage usage, const(double) * x, int t1, int t2,
                              int width, const(char) * vname, Gretlopt opt,
                              Prn * prn, GretlMatrixRaw ** pmat)
{
    double * dens = null;
    int bartlett = (opt & Gretlopt.OPT_O) != 0;
    int whittle_only = 0;
    int t, T, L = 0;
    int do_plot = 0;
    int err = 0;

    gretl_error_clear();

    /* common to all uses: check for data problems */

    err = series_adjust_sample(x, &t1, &t2);
    if (!err && (T = t2 - t1 + 1) < 12) {
        err = GretlError.toofew;
    }
    if (!err && gretl_isconst(t1, t2, x)) {
        if (vname != null) {
            gretl_errmsg_sprintf(_("%s is a constant").ptr, vname);
        }
        err = GretlError.data;
    }
    if (err) {
        return err;
    }

    if (usage == PergmUsage.pergm_cmd) {
        do_plot = gnuplot_graph_wanted(PlotType.periodogram, opt, &err);
        if (err) {
            return err;
        }
    }

    if (usage == PergmUsage.fractint_cmd) {
        /* doing fractint */
        if (!(opt & (Gretlopt.OPT_A | Gretlopt.OPT_G))) {
            /* not --all, not --gph */
            whittle_only = 1;
        }
    }

    if (!whittle_only) {
        /* Chatfield (1996); William Greene, 4e, p. 772 */
        if (bartlett) {
            if (width <= 0) {
                L = auto_spectrum_order(T, opt);
            } else {
                L = (width > T / 2)? T / 2 : width;
            }
        } else {
            /* use full sample */
            L = T - 1;
        }

        if (bartlett) {
            dens = pergm_bartlett_density(x, t1, t2, L, &err);
        } else {
            GretlMatrixRaw * pg;
            GretlVector vt;

            gretl_matrix_init(&vt);
            vt.rows = t2 - t1 + 1;
            vt.cols = 1;
            vt.val = cast(double *) x + t1;

            pg = gretl_matrix_pergm(&vt, T/2, &err);
            if (!err) {
                dens = gretl_matrix_steal_data(pg);
            }
            gretl_matrix_free(pg);
        }
        if (err) {
            return err;
        }
    }

    if (usage == PergmUsage.pergm_func) {
        /* make matrix for pergm() function */
        int T2 = T / 2;
        GretlMatrixRaw * pm;

        *pmat = pm = gretl_matrix_alloc(T2, 2);

        if (pm == null) {
            err = GretlError.alloc;
        } else {
            for (t=0; t<T2; t++) {
                gretl_matrix_set(pm, t, 0, M_2PI * (t+1) / cast(double) T);
                gretl_matrix_set(pm, t, 1, dens[t]);
            }
        }
    } else if (usage == PergmUsage.fractint_cmd) {
        /* supporting the "fractint" command */
        err = finalize_fractint(x, dens, t1, t2, width,
                                vname, opt, prn);
    } else {
        /* supporting the "pergm" command */
        if (!(opt & Gretlopt.OPT_S)) {
            /* not --silent */
            pergm_print(vname, dens, T, L, opt, prn);
        }
        if (do_plot) {
            err = periodogram_plot(vname, T, L, dens, opt);
        }
    }

    free(dens);

    return err;
}

/**
 * periodogram:
 * @varno: ID number of variable to process.
 * @width: width of window.
 * @dset: dataset struct.
 * @opt: if includes OPT_O, use Bartlett lag window for periodogram;
 * OPT_L, use log scale; OPT_S, don't print anything (just do plot).
 * @prn: gretl printing struct.
 *
 * Computes and displays the periodogram for the series specified
 * by @varno.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int periodogram(int varno, int width, const(DatasetRaw) * dset,
                Gretlopt opt, Prn * prn)
{
    if (!dataset_is_time_series(dset)) {
        return GretlError.tsonly;
    } else {
        return pergm_or_fractint(PergmUsage.pergm_cmd, dset.Z[varno],
                                 dset.t1, dset.t2,
                                 width, dset.varname[varno],
                                 opt, prn, null);
    }
}

/**
 * residual_periodogram:
 * @x: series to process.
 * @width: width of window.
 * @dset: dataset struct.
 * @opt: if includes OPT_O, use Bartlett lag window for periodogram;
 * OPT_L, use log scale.
 * @prn: gretl printing struct.
 *
 * Computes and displays the periodogram for @x, which is presumed to
 * be a model residual series.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int residual_periodogram(const(double) * x, int width, const(DatasetRaw) * dset,
                         Gretlopt opt, Prn * prn)
{
    return pergm_or_fractint(PergmUsage.pergm_cmd, x,
                             dset.t1, dset.t2,
                             width, null,
                             opt, prn, null);
}

/**
 * periodogram_matrix:
 * @x: the series to process.
 * @t1: starting observation in @x.
 * @t2: ending observation in @x.
 * @width: width of Bartlett window, or -1 for plain sample
 * periodogram.
 * @err: location to receive error code.
 *
 * Implements the userspace gretl pergm function, which can
 * be used on either a series from the dataset or a gretl
 * vector.
 *
 * Returns: allocated matrix on success, NULL on failure.
 */

GretlMatrixRaw * periodogram_matrix(const(double) * x, int t1, int t2,
                                 int width, int * err)
{
    Gretlopt opt = (width < 0)? Gretlopt.OPT_NONE : Gretlopt.OPT_O;
    GretlMatrixRaw * m = null;

    *err = pergm_or_fractint(PergmUsage.pergm_func, x, t1, t2, width,
                             null, opt, null, &m);

    return m;
}

/**
 * fractint:
 * @varno: ID number of variable to process.
 * @order: lag order / window size.
 * @dset: dataset struct.
 * @opt: option flags.
 * @prn: gretl printing struct.
 *
 * Computes and prints a test for fractional integration of the
 * series specified by @varno. By default the test uses the
 * Local Whittle Estimator but if @opt includes OPT_G then
 * the Geweke and Porter-Hudak test is done instead, or if
 * OPT_A then both tests are shown. If OPT_Q is given the
 * test results are not printed, just recorded (with
 * preference given to the LWE in case of OPT_A).
 *
 * Returns: 0 on successful completion, error code on error.
 */

int fractint(int varno, int order, const(DatasetRaw) * dset,
             Gretlopt opt, Prn * prn)
{
    int err = incompatible_options(opt, (Gretlopt.OPT_G | Gretlopt.OPT_A));

    if (!err && !dataset_is_time_series(dset)) {
        err = GretlError.tsonly;
    }
    if (!err) {
        err = pergm_or_fractint(PergmUsage.fractint_cmd, dset.Z[varno],
                                dset.t1, dset.t2,
                                order,
                                dset.varname[varno],
                                opt, prn, null);
    }

    return err;
}

private void printfw(int w, double x, int d, int ival, Prn * prn)
{
    w = w - 1; /* allow for leading space */
    pputc(prn, ' ');
    if (na(x)) {
        pprintf(prn, "%*s", UTF_WIDTH(_("NA"), w), _("NA").ptr);
    } else if (ival) {
        pprintf(prn, "%*.0f", w, x);
    } else if (x > 999 && x < 100000) {
        int p = cast(int) (1 + floor(log10(x)));

        d -= p;
        if (d < 0) d = 0;
        pprintf(prn, "%*.*f", w, d, x);
    } else {
        pprintf(prn, "%#*.*g", w, d, x);
    }
}

private void output_line(const(char) * str, Prn * prn, int dblspc)
{
    pputs(prn, str);
    if (dblspc) {
        pputs(prn, "\n\n");
    } else {
        pputc(prn, '\n');
    }
}

private void prhdr(const(char) * str, const(DatasetRaw) * dset,
                   int missing, Prn * prn)
{
    char[OBSLEN] date1 = 0;
    char[OBSLEN] date2 = 0;
    gchar * tmp;

    ntolabel(date1.ptr, dset.t1, dset);
    ntolabel(date2.ptr, dset.t2, dset);

    pputc(prn, '\n');

    tmp = g_strdup_printf(_("%s, using the observations %s - %s").ptr,
                          str, date1.ptr, date2.ptr);
    output_line(tmp, prn, 0);
    g_free(tmp);

    if (missing) {
        output_line(_("(missing values were skipped)").ptr, prn, 1);
    }
}

private void summary_print_val(int w, double x, int digits,
                               int places, Prn * prn)
{
    pputc(prn, ' ');
    if (na(x)) {
        pprintf(prn, "%*s", UTF_WIDTH(_("NA"), w), _("NA").ptr);
    } else if (digits < 0) {
        pprintf(prn, "%*.0f", w, x);
    } else if (digits > 0) {
        pprintf(prn, "%#*.*g", w, digits, x);
    } else if (places > 0) {
        pprintf(prn, "%*.*f", w, places, x);
    } else {
        /* the default */
        pprintf(prn, "%#*.5g", w, x);
    }
}

private enum NSUMM = 12;

void print_summary_single(const(Summary) * s,
                          int digits, int places,
                          const(DatasetRaw) * dset,
                          Prn * prn)
{
    string[NSUMM] labels = [
        N_("Mean"),
        N_("Median"),
        N_("Minimum"),
        N_("Maximum"),
        N_("Standard deviation"),
        N_("C.V."),
        N_("Skewness"),
        N_("Ex. kurtosis"),
        /* xgettext:no-c-format */
        N_("5% percentile"),
        /* xgettext:no-c-format */
        N_("95% percentile"),
        N_("Interquartile range"),
        N_("Missing obs.")
    ];
    string wstr = N_("Within s.d.");
    string bstr = N_("Between s.d.");
    double[NSUMM] vals;
    int[NSUMM] simple_skip = [0,1,0,0,0,1,1,1,1,1,1,0];
    int skip0595 = 0;
    int offset = 2;
    int slen = 0, i = 0;

    if (s.opt & Gretlopt.OPT_B) {
        offset = 4;
    } else {
        const(char) * vname = dset.varname[s.list[1]];
        char[OBSLEN] obs1 = 0;
        char[OBSLEN] obs2 = 0;
        gchar * tmp = null;

        ntolabel(obs1.ptr, dset.t1, dset);
        ntolabel(obs2.ptr, dset.t2, dset);

        prhdr(_("Summary statistics").ptr, dset, 0, prn);

        if (isdigit(*vname)) {
            const(char) * mname = dataset_get_matrix_name(dset);

            if (mname != null) {
                tmp = g_strdup_printf(_("for column %d of %s (%d valid observations)").ptr,
                                      atoi(vname), mname, s.n);
            } else {
                tmp = g_strdup_printf(_("for column %d (%d valid observations)").ptr,
                                      atoi(vname), s.n);
            }
        } else {
            tmp = g_strdup_printf(_("for the variable '%s' (%d valid observations)").ptr,
                                  dset.varname[s.list[1]], s.n);
        }
        output_line(tmp, prn, 1);
        g_free(tmp);
    }

    vals[0]  = s.mean[0];
    vals[1]  = s.median[0];
    vals[2]  = s.low[0];
    vals[3]  = s.high[0];
    vals[4]  = s.sd[0];
    vals[5]  = s.cv[0];
    vals[6]  = s.skew[0];
    vals[7]  = s.xkurt[0];
    vals[8]  = s.perc05[0];
    vals[9]  = s.perc95[0];
    vals[10] = s.iqr[0];
    vals[11] = s.misscount[0];

    if (na(vals[8]) && na(vals[9])) {
        skip0595 = 1;
    }

    for (i=0; i<NSUMM; i++) {
        if ((s.opt & Gretlopt.OPT_S) && simple_skip[i]) {
            continue;
        } else if (skip0595 && (i == 8 || i == 9)) {
            continue;
        }
        if (strlen(_(labels[i]).ptr) > slen) {
            slen = cast(int) g_utf8_strlen(_(labels[i]).ptr, -1);
        }
    }
    slen++;

    for (i=0; i<NSUMM; i++) {
        if ((s.opt & Gretlopt.OPT_S) && simple_skip[i]) {
            continue;
        } else if (skip0595 && (i == 8 || i == 9)) {
            continue;
        }
        bufspace(offset, prn);
        if (i == 8 || i == 9) {
            /* the strings contain '%' */
            gchar * pcstr = g_strdup(_(labels[i]).ptr);
            int n = cast(int) (slen - g_utf8_strlen(pcstr, -1));

            pputs(prn, pcstr);
            if (n > 0) {
                bufspace(n, prn);
            }
            g_free(pcstr);
        } else {
            pprintf(prn, "%-*s", UTF_WIDTH(_(labels[i]), slen), _(labels[i]).ptr);
        }
        if (i == NSUMM - 1 || s.ival[0]) {
            summary_print_val(14, vals[i], -1, places, prn);
        } else {
            summary_print_val(14, vals[i], digits, places, prn);
        }
        pputc(prn, '\n');
    }

    if (!na(s.sw) && !na(s.sb)) {
        pputc(prn, '\n');
        bufspace(offset, prn);
        pprintf(prn, "%-*s", UTF_WIDTH(_(wstr), slen), _(wstr).ptr);
        summary_print_val(14, s.sw, digits, places, prn);
        pputc(prn, '\n');
        bufspace(offset, prn);
        pprintf(prn, "%-*s", UTF_WIDTH(_(bstr), slen), _(bstr).ptr);
        summary_print_val(14, s.sb, digits, places, prn);
    }

    pputc(prn, '\n');
}

private void summary_print_varname(const(char) * src,
                                   int len, Prn * prn)
{
    char[SHORTLEN] vname = 0;

    maybe_trim_varname(vname.ptr, src);
    pprintf(prn, "%-*s", len, vname.ptr);
}

/**
 * print_summary:
 * @summ: pointer to gretl summary statistics struct.
 * @dset: information on the data set.
 * @prn: gretl printing struct.
 *
 * Prints the summary statistics for a given variable.
 */

void print_summary(const(Summary) * summ,
                   const(DatasetRaw) * dset,
                   Prn * prn)
{
    int dmax, d = get_gretl_digits();
    int nv = summ.list[0];
    int len = 0;
    int i, vi;

    if (summ.list == null || summ.list[0] == 0) {
        return;
    }

    if (summ.weight_var > 0) {
        pputc(prn, '\n');
        pprintf(prn, _("Weighting variable: %s\n").ptr,
                dset.varname[summ.weight_var]);
    }

    if (nv == 1 && !(summ.opt & (Gretlopt.OPT_B | Gretlopt.OPT_M))) {
        print_summary_single(summ, 6, 0, dset, prn);
        return;
    }

    /* number of significant figures to use */
    dmax = 5; // was: (summ->opt & OPT_S)? 4 : 5;
    d = d > dmax ? dmax : d;

    if (nv > 1 || (summ.opt & Gretlopt.OPT_M)) {
        int maxlen = max_namelen_in_list(summ.list, dset);

        len = maxlen + 1;
    }

    pputc(prn, '\n');

    if (summ.opt & Gretlopt.OPT_S) {
        /* the "simple" option */
        string[5] h = [
            N_("Mean"),
            N_("Median"),
            /* TRANSLATORS: 'S.D.' means Standard Deviation */
            N_("S.D."),
            N_("Min"),
            N_("Max"),
        ];
        int w = 12;
        int ival;

        if (len > 0) {
            pprintf(prn, "%*s", len, " ".ptr);
        }
        pprintf(prn, "%*s%*s%*s%*s%*s\n",
                UTF_WIDTH(_(h[0]), w), _(h[0]).ptr,
                UTF_WIDTH(_(h[1]), w), _(h[1]).ptr,
                UTF_WIDTH(_(h[2]), w), _(h[2]).ptr,
                UTF_WIDTH(_(h[3]), w), _(h[3]).ptr,
                UTF_WIDTH(_(h[4]), w), _(h[4]).ptr);

        for (i=0; i<summ.list[0]; i++) {
            vi = summ.list[i+1];
            ival = summ.ival[i];
            if (nv > 1 || (summ.opt & Gretlopt.OPT_M)) {
                summary_print_varname(dset.varname[vi], len, prn);
            }
            printfw(w, summ.mean[i], d, ival, prn);
            printfw(w, summ.median[i], d, ival, prn);
            printfw(w, summ.sd[i], d, ival, prn);
            printfw(w, summ.low[i], d, ival, prn);
            printfw(w, summ.high[i], d, ival, prn);
            pputc(prn, '\n');
        }
        pputc(prn, '\n');
    } else {
        /* print all available stats */
        string[4] ha = [
            N_("Mean"),
            N_("Median"),
            N_("Minimum"),
            N_("Maximum"),
        ];
        string[4] hb = [
            N_("Std. Dev."),
            N_("C.V."),
            N_("Skewness"),
            N_("Ex. kurtosis")
        ];
        string[4] hc = [
            /* xgettext:no-c-format */
            N_("5% perc."),
            /* xgettext:no-c-format */
            N_("95% perc."),
            N_("IQ range"),
            N_("Missing obs.")
        ];
        /* cases where 0.05 and 0.95 quantiles are OK */
        int npct = 0;
        int w = 15;
        int ival;

        if (len > 0) {
            pprintf(prn, "%*s", len, " ".ptr);
        }
        pprintf(prn, "%*s%*s%*s%*s\n",
                UTF_WIDTH(_(ha[0]), w), _(ha[0]).ptr,
                UTF_WIDTH(_(ha[1]), w), _(ha[1]).ptr,
                UTF_WIDTH(_(ha[2]), w), _(ha[2]).ptr,
                UTF_WIDTH(_(ha[3]), w), _(ha[3]).ptr);

        for (i=0; i<summ.list[0]; i++) {
            vi = summ.list[i+1];
            ival = summ.ival[i];
            if (nv > 1) {
                summary_print_varname(dset.varname[vi], len, prn);
            }
            printfw(w, summ.mean[i], d, ival, prn);
            printfw(w, summ.median[i], d, ival, prn);
            printfw(w, summ.low[i], d, ival, prn);
            printfw(w, summ.high[i], d, ival, prn);
            pputc(prn, '\n');
            /* while we're at it, register cases where we can
               show the 0.05 and 0.95 quantiles
            */
            if (!na(summ.perc05[i]) && !na(summ.perc95[i])) {
                npct++;
            }
        }
        pputc(prn, '\n');

        if (len > 0) {
            pprintf(prn, "%*s", len, " ".ptr);
        }
        pprintf(prn, "%*s%*s%*s%*s\n",
                UTF_WIDTH(_(hb[0]), w), _(hb[0]).ptr,
                UTF_WIDTH(_(hb[1]), w), _(hb[1]).ptr,
                UTF_WIDTH(_(hb[2]), w), _(hb[2]).ptr,
                UTF_WIDTH(_(hb[3]), w), _(hb[3]).ptr);

        for (i=0; i<summ.list[0]; i++) {
            double cv;

            vi = summ.list[i+1];
            ival = summ.ival[i];
            if (nv > 1) {
                summary_print_varname(dset.varname[vi], len, prn);
            }
            if (floateq(summ.mean[i], 0.0)) {
                cv = NADBL;
            } else if (floateq(summ.sd[i], 0.0)) {
                cv = 0.0;
            } else {
                cv = fabs(summ.sd[i] / summ.mean[i]);
            }
            printfw(w, summ.sd[i], d, ival, prn);
            printfw(w, cv, d, 0, prn);
            printfw(w, summ.skew[i], d, 0, prn);
            printfw(w, summ.xkurt[i], d, 0, prn);
            pputc(prn, '\n');
        }
        pputc(prn, '\n');

        if (npct > 0) {
            /* note: use pputs for strings containing literal '%' */
            gchar * hc0 = g_strdup(_(hc[0]).ptr);
            gchar * hc1 = g_strdup(_(hc[1]).ptr);
            int n;

            if (len > 0) {
                pprintf(prn, "%*s", len, " ".ptr);
            }
            n = cast(int) (w - g_utf8_strlen(hc0, -1));
            if (n > 0) bufspace(n, prn);
            pputs(prn, hc0);
            n = cast(int) (w - g_utf8_strlen(hc1, -1));
            if (n > 0) bufspace(n, prn);
            pputs(prn, hc1);
            g_free(hc0); g_free(hc1);

            pprintf(prn, "%*s%*s\n",
                    UTF_WIDTH(_(ha[2]), w), _(hc[2]).ptr,
                    UTF_WIDTH(_(ha[3]), w), _(hc[3]).ptr);
        } else {
            /* not showing any 0.05, 0.95 quantiles */
            if (len > 0) {
                pprintf(prn, "%*s", len, " ".ptr);
            }
            pprintf(prn, "%*s%*s\n",
                    UTF_WIDTH(_(ha[2]), w), _(hc[2]).ptr,
                    UTF_WIDTH(_(ha[3]), w), _(hc[3]).ptr);
        }

        for (i=0; i<summ.list[0]; i++) {
            vi = summ.list[i+1];
            ival = summ.ival[i];
            if (nv > 1) {
                summary_print_varname(dset.varname[vi], len, prn);
            }
            if (!na(summ.perc05[i]) && !na(summ.perc95[i])) {
                printfw(w, summ.perc05[i], d, ival, prn);
                printfw(w, summ.perc95[i], d, ival, prn);
            } else if (npct > 0) {
                pprintf(prn, "%*s", w, "NA".ptr);
                pprintf(prn, "%*s", w, "NA".ptr);
            }
            printfw(w, summ.iqr[i], d, ival, prn);
            pprintf(prn, "%*d", w, cast(int) summ.misscount[i]);
            pputc(prn, '\n');
        }
        pputc(prn, '\n');
    }
}

/**
 * free_summary:
 * @summ: pointer to gretl summary statistics struct
 *
 * Frees all resources associated with @summ, and
 * the pointer itself.
 */

void free_summary(Summary * summ)
{
    free(summ.list);
    free(summ.misscount);
    free(summ.stats);
    free(summ.ival);

    free(summ);
}

/***
 * summary_new:
 * @list: list of variables we want summary statistics for
 * @wv: weighting variable (0=no weights)
 * @opt: options
 *
 * Allocates a new "Summary" struct and initializes a few
 * things inside it
 */

private Summary * summary_new(const(int) * list, int wv,
                              Gretlopt opt, int * err)
{
    Summary * s;
    int nv;

    if (list == null) {
        fprintf(stderr, "summary_new: list is NULL\n");
        *err = GretlError.data;
        return null;
    }

    nv = list[0];

    s = cast(Summary *) malloc((*s).sizeof);
    if (s == null) {
        *err = GretlError.alloc;
        return null;
    }

    s.list = gretl_list_copy(list);
    if (s.list == null) {
        *err = GretlError.alloc;
        free(s);
        return null;
    }

    s.weight_var = wv;
    s.opt = opt;
    s.n = 0;
    s.misscount = cast(int *) malloc(nv * (*s.misscount).sizeof);
    s.ival = cast(int *) calloc(nv, (*s.ival).sizeof);

    s.stats = cast(double *) malloc(11 * nv * (*s.stats).sizeof);
    if (s.stats == null) {
        *err = GretlError.alloc;
        free_summary(s);
        return null;
    }

    s.mean   = s.stats;
    s.median = s.mean + nv;
    s.sd     = s.median + nv;
    s.skew   = s.sd + nv;
    s.xkurt  = s.skew + nv;
    s.low    = s.xkurt + nv;
    s.high   = s.low + nv;
    s.cv     = s.high + nv;
    s.perc05 = s.cv + nv;
    s.perc95 = s.perc05 + nv;
    s.iqr    = s.perc95 + nv;

    s.sb = s.sw = NADBL;

    return s;
}

int summary_has_missing_values(const(Summary) * summ)
{
    if (summ.misscount != null) {
        int i, nv = summ.list[0];

        for (i=0; i<nv; i++) {
            if (summ.misscount[i] > 0) {
                return 1;
            }
        }
    }

    return 0;
}

private extern (C) int compare_wgt_ord_rows(const(void) * a, const(void) * b)
{
    const(double) ** da = cast(const(double) **) a;
    const(double) ** db = cast(const(double) **) b;

    return (da[0][0] > db[0][0]) - (da[0][0] < db[0][0]);
}

private int weighted_order_stats(const(double) * y, const(double) * w,
                                 double * ostats, int n, int k)
{
    /* on input "ostats" contains the quantiles to compute;
       on output it holds the results
    */
    double p, q, wsum = 0.0;
    double ** X;
    int t, i, err = 0;

    X = doubles_array_new(n, 2);
    if (X == null) {
        return GretlError.alloc;
    }

    i = 0;
    for (t=0; t<n; t++) {
        if (na(y[t]) || na(w[t]) || w[t] == 0.0) {
            continue;
        } else {
            X[i][0] = y[t];
            X[i][1] = w[t];
            wsum += w[t];
            i++;
        }
    }

    qsort(X, i, (*X).sizeof, &compare_wgt_ord_rows);

    for (i=0; i<k; i++) {
        p = ostats[i] * wsum;
        q = X[0][1];
        t = 0;
        while (q <= p) {
            q += X[t][1];
            if (q < p) t++;
        }

        if (t == 0 || t >= n - 1) {
            ostats[i] = NADBL;
            continue;
        }

        if (X[t-1][0] == X[t][0]) {
            ostats[i] = X[t][0];
        } else {
            double a = (q - p) / X[t][1];

            ostats[i] = a * X[t+1][0] + (1-a) * X[t][0];
        }
    }

    doubles_array_free(X, n);

    return err;
}

/**
 * get_summary_weighted:
 * @list: list of variables to process.
 * @dset: dataset struct.
 * @wgt: series to use as weights.
 * @opt: may include OPT_S for "simple" version.
 * @prn: gretl printing struct.
 * @err: location to receive error code.
 *
 * Calculates descriptive summary statistics for the specified
 * variables, weighting the observations by @rv. The series @rv must
 * be of full length (dset->n).
 *
 * Returns: #Summary object containing the summary statistics, or NULL
 * on failure.
 */

Summary * get_summary_weighted(const(int) * list, const(DatasetRaw) * dset,
                               int wtvar, Gretlopt opt,
                               Prn * prn, int * err)
{
    const(double) * wts;
    double * x;
    double[5] ostats;
    int t1 = dset.t1;
    int t2 = dset.t2;
    Summary * s;
    int i, t;

    s = summary_new(list, wtvar, opt, err);
    if (s == null) {
        return null;
    }

    x = cast(double *) malloc(dset.n * (*x).sizeof);
    if (x == null) {
        *err = GretlError.alloc;
        free_summary(s);
        return null;
    }

    wts = dset.Z[wtvar];

    for (i=0; i<list[0]; i++)  {
        double * pskew = null;
        double * pkurt = null;
        int vi = s.list[i+1];
        int ni = 0;
        int ntot = 0;
        double wt;

        /* prepare the series for weighting: substitute NAs for values
           at which the weights are invalid or zero
        */
        for (t=t1; t<=t2; t++) {
            wt = wts[t];
            if (!na(wt) && wt != 0.0) {
                ntot++;
                x[t] = dset.Z[vi][t];
                if (!na(x[t])) {
                    ni++;
                }
            } else {
                x[t] = NADBL;
            }
        }

        s.misscount[i] = ntot - ni;

        if (ni > s.n) {
            s.n = ni;
        }

        if (ni == 0) {
            pprintf(prn, _("Dropping %s: sample range contains no valid "
                           ~ "observations\n").ptr, dset.varname[vi]);
            gretl_list_delete_at_pos(s.list, i + 1);
            if (s.list[0] == 0) {
                return s;
            } else {
                i--;
                continue;
            }
        }

        if (opt & Gretlopt.OPT_S) {
            s.skew[i] = NADBL;
            s.xkurt[i] = NADBL;
            s.cv[i] = NADBL;
            s.median[i] = NADBL;
        } else {
            pskew = &s.skew[i];
            pkurt = &s.xkurt[i];
        }

        summary_minmax(t1, t2, x, s, i);
        gretl_moments(t1, t2, x, wts, &s.mean[i], &s.sd[i], pskew, pkurt, 1);

        if (!(opt & Gretlopt.OPT_S)) {
            int werr;       /* C: "err", shadowing the argument */

            if (floateq(s.mean[i], 0.0)) {
                s.cv[i] = NADBL;
            } else if (floateq(s.sd[i], 0.0)) {
                s.cv[i] = 0.0;
            } else {
                s.cv[i] = fabs(s.sd[i] / s.mean[i]);
            }

            ostats[0] = 0.05;
            ostats[1] = 0.25;
            ostats[2] = 0.5;
            ostats[3] = 0.75;
            ostats[4] = 0.95;

            werr = weighted_order_stats(x, wts, ostats.ptr, ni, 5);

            if (!werr) {
                s.median[i] = ostats[2];
                s.perc05[i] = ostats[0];
                s.perc95[i] = ostats[4];
                if (!na(ostats[1]) && !na(ostats[3])) {
                    s.iqr[i] = ostats[3] - ostats[1];
                }
            }
        }
    }

    free(x);

    return s;
}

/* Get the additional statistics wanted if the --simple
   flag was not given to the summary command.
*/

private int get_extra_stats(Summary * s, int i,
                            int t1, int t2,
                            const(double) * x)
{
    int err = 0;

    if (floateq(s.mean[i], 0.0)) {
        s.cv[i] = NADBL;
    } else if (floateq(s.sd[i], 0.0)) {
        s.cv[i] = 0.0;
    } else {
        s.cv[i] = fabs(s.sd[i] / s.mean[i]);
    }

    s.perc05[i] = gretl_quantile(t1, t2, x, 0.05, Gretlopt.OPT_Q, &err);
    s.perc95[i] = gretl_quantile(t1, t2, x, 0.95, Gretlopt.OPT_Q, &err);
    s.iqr[i] = gretl_quantile(t1, t2, x, 0.75, Gretlopt.OPT_NONE, &err);
    if (!na(s.iqr[i])) {
        s.iqr[i] -= gretl_quantile(t1, t2, x, 0.25, Gretlopt.OPT_NONE, &err);
    }

    return err;
}

/**
 * get_summary_restricted:
 * @list: list of variables to process.
 * @dset: dataset struct.
 * @rv: series to use as restriction dummy.
 * @opt: may include OPT_S for "simple" version.
 * @prn: gretl printing struct.
 * @err: location to receive error code.
 *
 * Calculates descriptive summary statistics for the specified variables,
 * with the observations restricted to those for which @rv has a non-zero
 * (and non-missing) value. The series @rv must be of full length (dset->n).
 *
 * Returns: #Summary object containing the summary statistics, or NULL
 * on failure.
 */

Summary * get_summary_restricted(const(int) * list, const(DatasetRaw) * dset,
                                 const(double) * rv, Gretlopt opt,
                                 Prn * prn, int * err)
{
    int t1 = dset.t1;
    int t2 = dset.t2;
    Summary * s;
    double * x;
    int i, t;

    s = summary_new(list, 0, opt, err);
    if (s == null) {
        return null;
    }

    x = cast(double *) malloc(dset.n * (*x).sizeof);
    if (x == null) {
        *err = GretlError.alloc;
        free_summary(s);
        return null;
    }

    for (i=0; i<s.list[0]; i++)  {
        double * pskew = null;
        double * pkurt = null;
        int vi = s.list[i+1];
        int strvals;
        int ni = 0;
        int ntot = 0;

        strvals = is_string_valued(dset, vi);

        if (!strvals) {
            /* create the restricted series: substitute NAs
               for values at which the restriction dummy is
               invalid or zero
            */
            for (t=t1; t<=t2; t++) {
                if (!na(rv[t]) && rv[t] != 0.0) {
                    ntot++;
                    x[t] = dset.Z[vi][t];
                    if (!na(x[t])) {
                        ni++;
                    }
                } else {
                    x[t] = NADBL;
                }
            }
            s.misscount[i] = ntot - ni;
            if (ni > s.n) {
                s.n = ni;
            }
        }

        if (ni == 0) {
            if (strvals) {
                pprintf(prn, _("Dropping %s: string-valued series\n").ptr,
                        dset.varname[vi]);
            } else {
                pprintf(prn, _("Dropping %s: sample range contains no valid "
                               ~ "observations\n").ptr, dset.varname[vi]);
            }
            gretl_list_delete_at_pos(s.list, i + 1);
            if (s.list[0] == 0) {
                return s;
            } else {
                i--;
                continue;
            }
        }

        if (opt & Gretlopt.OPT_S) {
            s.skew[i] = NADBL;
            s.xkurt[i] = NADBL;
            s.cv[i] = NADBL;
        } else {
            pskew = &s.skew[i];
            pkurt = &s.xkurt[i];
        }

        summary_minmax(t1, t2, x, s, i);
        gretl_moments(t1, t2, x, null, &s.mean[i], &s.sd[i], pskew, pkurt, 1);
        s.median[i] = gretl_median(t1, t2, x);

        if (!(opt & Gretlopt.OPT_S)) {
            *err = get_extra_stats(s, i, t1, t2, x);
        }

        if (dataset_is_panel(dset) && list[0] == 1) {
            panel_variance_info(x, dset, s.mean[0], &s.sw, &s.sb);
        }
    }

    free(x);

    return s;
}

/**
 * get_summary:
 * @list: list of variables to process.
 * @dset: dataset struct.
 * @opt: may include OPT_S for "simple" version.
 * @prn: gretl printing struct.
 * @err: location to receive error code.
 *
 * Calculates descriptive summary statistics for the specified variables.
 *
 * Returns: #Summary object containing the summary statistics, or NULL
 * on failure.
 */

Summary * get_summary(const(int) * list, const(DatasetRaw) * dset,
                      Gretlopt opt, Prn * prn, int * err)
{
    int t1 = dset.t1;
    int t2 = dset.t2;
    Summary * s;
    int i, nmax;

    s = summary_new(list, 0, opt, err);
    if (s == null) {
        return null;
    }

    nmax = sample_size(dset);

    for (i=0; i<s.list[0]; i++)  {
        double * pskew = null;
        double * pkurt = null;
        const(double) * x;
        double x0;
        int vi = s.list[i+1];
        int strvals;
        int ni = 0;

        strvals = is_string_valued(dset, vi);
        if (!strvals) {
            x = dset.Z[vi];
            ni = good_obs(x + t1, nmax, &x0);
            s.misscount[i] = nmax - ni;
            if (ni > s.n) {
                s.n = ni;
            }
        }
        if (ni == 0) {
            if (strvals) {
                pprintf(prn, _("Dropping %s: string-valued series\n").ptr,
                        dset.varname[vi]);
            } else {
                pprintf(prn, _("Dropping %s: sample range contains no valid "
                               ~ "observations\n").ptr, dset.varname[vi]);
            }
            gretl_list_delete_at_pos(s.list, i + 1);
            if (s.list[0] == 0) {
                return s;
            } else {
                i--;
                continue;
            }
        }

        if (opt & Gretlopt.OPT_S) {
            s.skew[i] = NADBL;
            s.xkurt[i] = NADBL;
            s.cv[i] = NADBL;
            s.median[i] = NADBL;
        } else {
            pskew = &s.skew[i];
            pkurt = &s.xkurt[i];
        }

        summary_minmax(t1, t2, x, s, i);

        if (s.weight_var == 0) {
            gretl_moments(t1, t2, x, null, &s.mean[i], &s.sd[i],
                          pskew, pkurt, 1);
        } else {
            gretl_moments(t1, t2, x, dset.Z[s.weight_var], &s.mean[i], &s.sd[i],
                          pskew, pkurt, 0);
        }

        /* included in both simple and full variants */
        s.median[i] = gretl_median(t1, t2, x);

        if (!(opt & Gretlopt.OPT_S)) {
            *err = get_extra_stats(s, i, t1, t2, x);
        }

        if (dataset_is_panel(dset) && list[0] == 1) {
            panel_variance_info(x, dset, s.mean[0], &s.sw, &s.sb);
        }
    }

    return s;
}

/* basic summary stats column heads */
private immutable string[6] hb = [
    N_("Byval"), N_("Mean"), N_("Median"), N_("S.D."), N_("Min"), N_("Max")
];

/* full summary stats column heads */
private immutable string[13] hf = [
    N_("Byval"), N_("Mean"), N_("Median"), N_("Minimum"), N_("Maximum"),
    N_("Std. Dev."), N_("C.V."), N_("Skewness"), N_("Ex. kurtosis"),
    /* xgettext:no-c-format */
    N_("5% perc."),
    /* xgettext:no-c-format */
    N_("95% perc."),
    N_("IQ range"),
    N_("Missing obs.")
];

private int add_summary_colnames(GretlMatrixRaw * m,
                                 Gretlopt opt)
{
    char ** S = strings_array_new(m.cols);

    if (S == null) {
        return GretlError.alloc;
    } else {
        immutable(string)[] h = (opt & Gretlopt.OPT_S)? hb[] : hf[];
        int k = (opt & Gretlopt.OPT_B) ? 0 : 1;
        int i;

        for (i=0; i<m.cols; i++) {
            S[i] = gretl_strdup(_(h[k++]).ptr);
        }
        gretl_matrix_set_colnames(m, S);
    }

    return 0;
}

/* the factorized statistics case */

private int add_summary_by_rownames(GretlMatrixRaw * m,
                                    const(GretlMatrixRaw) * xvals,
                                    const(DatasetRaw) * dset,
                                    SeriesTable * st)
{
    char ** S = strings_array_new(m.rows);
    gchar * tmp;
    double xvi;
    int i;

    if (S == null) {
        return GretlError.alloc;
    } else {
        for (i=0; i<m.rows; i++) {
            xvi = xvals.val[i];
            if (st != null) {
                S[i] = gretl_strdup(series_table_get_string(st, xvi));
            } else {
                tmp = gretl_strdup_printf("%g", xvi);
                S[i] = gretl_strdup(tmp);
                g_free(tmp);
            }
        }
        gretl_matrix_set_rownames(m, S);
    }

    return 0;
}

/* the (relatively) simple case */

private int add_summary_rownames(GretlMatrixRaw * m,
                                 const(int) * list,
                                 const(DatasetRaw) * dset)
{
    char ** S = strings_array_new(m.rows);
    int i, vi;

    if (S == null) {
        return GretlError.alloc;
    } else {
        for (i=0; i<m.rows; i++) {
            vi = list[i+1];
            S[i] = gretl_strdup(dset.varname[vi]);
        }
        gretl_matrix_set_rownames(m, S);
    }

    return 0;
}

private GretlMatrixRaw * make_summary_recorder(const(int) * list,
                                            Gretlopt opt,
                                            const(DatasetRaw) * dset,
                                            const(GretlMatrixRaw) * xvals,
                                            SeriesTable * st,
                                            int * err)
{
    GretlMatrixRaw * m = null;
    int rows, cols;

    cols = (opt & Gretlopt.OPT_S)? 5 : 12;

    if (opt & Gretlopt.OPT_B) {
        /* add an extra column for factor value */
        cols++;
        /* number of distinct factor values */
        rows = gretl_vector_get_length(xvals);
    } else {
        rows = list[0];
    }

    m = gretl_zero_matrix_new(rows, cols);

    if (m == null) {
        *err = GretlError.alloc;
    } else {
        *err = add_summary_colnames(m, opt);
        if (!*err) {
            if (xvals != null) {
                *err = add_summary_by_rownames(m, xvals, dset, st);
            } else {
                *err = add_summary_rownames(m, list, dset);
            }
        }
    }

    if (*err && m != null) {
        gretl_matrix_free(m);
        m = null;
    }

    return m;
}

/* If @ival is NULL, write the content of @summ into matrix @m,
   otherwise
*/

private void transcribe_summary(Summary * summ, GretlMatrixRaw * m,
                                int * ival, double byval)
{
    int factorized = ival != null;
    int nv = summ.list[0];
    int offset = 0;
    int i, r;

    if (factorized) {
        /* writing to a single row */
        gretl_matrix_set(m, *ival, 0, byval);
        offset = 1;
    }

    if (summ.opt & Gretlopt.OPT_S) {
        /* the "simple" case */
        for (i=0; i<nv; i++) {
            r = factorized ? *ival : i;
            gretl_matrix_set(m, r, 0 + offset, summ.mean[i]);
            gretl_matrix_set(m, r, 1 + offset, summ.median[i]);
            gretl_matrix_set(m, r, 2 + offset, summ.sd[i]);
            gretl_matrix_set(m, r, 3 + offset, summ.low[i]);
            gretl_matrix_set(m, r, 4 + offset, summ.high[i]);
        }
    } else {
        /* record all available stats */
        double cv;

        for (i=0; i<nv; i++) {
            r = factorized ? *ival : i;
            gretl_matrix_set(m, r, 0 + offset, summ.mean[i]);
            gretl_matrix_set(m, r, 1 + offset, summ.median[i]);
            gretl_matrix_set(m, r, 2 + offset, summ.low[i]);
            gretl_matrix_set(m, r, 3 + offset, summ.high[i]);
            gretl_matrix_set(m, r, 4 + offset, summ.sd[i]);
            if (floateq(summ.mean[i], 0.0)) {
                cv = NADBL;
            } else if (floateq(summ.sd[i], 0.0)) {
                cv = 0.0;
            } else {
                cv = fabs(summ.sd[i] / summ.mean[i]);
            }
            gretl_matrix_set(m, r,  5 + offset, cv);
            gretl_matrix_set(m, r,  6 + offset, summ.skew[i]);
            gretl_matrix_set(m, r,  7 + offset, summ.xkurt[i]);
            gretl_matrix_set(m, r,  8 + offset, summ.perc05[i]);
            gretl_matrix_set(m, r,  9 + offset, summ.perc95[i]);
            gretl_matrix_set(m, r, 10 + offset, summ.iqr[i]);
            gretl_matrix_set(m, r, 11 + offset, summ.misscount[i]);
        }
    }
}

private void record_summary(Summary * summ, const(DatasetRaw) * dset)
{
    GretlMatrixRaw * m;
    int err = 0;

    m = make_summary_recorder(summ.list, summ.opt, dset,
                              null, null, &err);
    if (!err) {
        transcribe_summary(summ, m, null, 0);
        set_last_result_data(m, GretlType.matrix);
    }
}

/**
 * list_summary:
 * @list: list of series to process.
 * @dset: dataset struct.
 * @opt: may include %OPT_S for "simple" version.
 * @prn: gretl printing struct.
 *
 * Prints descriptive statistics for the listed variables.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int list_summary(const(int) * list, int wgtvar,
                 const(DatasetRaw) * dset,
                 Gretlopt opt, Prn * prn)
{
    Summary * summ = null;
    int err = 0;

    if (list != null) {
        if (list[0] == 0) {
            return 0;
        }
        if (wgtvar == 0) {
            /* no weights */
            summ = get_summary(list, dset, opt, prn, &err);
        } else {
            summ = get_summary_weighted(list, dset, wgtvar,
                                        opt, prn, &err);
        }
    } else {
        int * tmplist = full_var_list(dset, null);

        if (tmplist != null) {
            if (wgtvar == 0) {
                /* no weights */
                summ = get_summary(tmplist, dset, opt, prn, &err);
            } else {
                summ = get_summary_weighted(tmplist, dset, wgtvar,
                                            opt, prn, &err);
            }
            free(tmplist);
        }
    }

    if (summ != null) {
        if (!(opt & Gretlopt.OPT_Q)) {
            print_summary(summ, dset, prn);
        }
        record_summary(summ, dset);
        free_summary(summ);
    }

    return err;
}

int summary_statistics_by(const(int) * list, DatasetRaw * dset,
                          Gretlopt opt, Prn * prn)
{
    const(char) * byname = get_optval_string(GretlCmdIndex.summary, Gretlopt.OPT_B);
    int quiet = (opt & Gretlopt.OPT_Q) != 0;
    int byvar;
    SeriesTable * st = null;
    GretlMatrixRaw * xvals = null;
    GretlMatrixRaw * m = null;
    int * mylist = null;
    const(double) * x;
    int i, nvals = 0;
    int single, err = 0;

    if (dset == null || byname == null) {
        return GretlError.data;
    }

    /* FIXME accept "unit" and "time"/"period" in place of actual
       variables for panel data? */

    byvar = current_series_index(dset, byname);
    if (byvar < 0) {
        return GretlError.unkvar;
    } else if (!accept_as_discrete(dset, byvar, 0)) {
        gretl_errmsg_sprintf(_("The variable '%s' is not discrete").ptr, byname);
        return GretlError.data;
    }

    x = cast(const(double) *) dset.Z[byvar];

    if (list == null) {
        /* compose full series list, minus the "by" variable */
        int pos;

        mylist = full_var_list(dset, null);
        if (mylist == null) {
            return GretlError.alloc;
        }
        pos = in_gretl_list(mylist, byvar);
        if (pos > 0) {
            gretl_list_delete_at_pos(mylist, pos);
        }
        if (mylist[0] == 0) {
            free(mylist);
            return GretlError.data;
        }
        list = mylist;
    }

    single = (list[0] == 1);

    xvals = gretl_matrix_values(x + dset.t1, dset.t2 - dset.t1 + 1,
                                Gretlopt.OPT_S, &err);

    if (!err) {
        nvals = gretl_vector_get_length(xvals);
        if (nvals == 0) {
            err = GretlError.data;
        } else {
            st = series_get_string_table(dset, byvar);
        }
    }

    if (!err && single) {
        m = make_summary_recorder(list, opt, dset, xvals, st, &err);
        if (!quiet) {
            pputc(prn, '\n');
            pprintf(prn, _("Summary statistics for %s, by value of %s").ptr,
                    dset.varname[list[1]], byname);
            pputc(prn, '\n');
        }
    }

    for (i=0; i<nvals && !err; i++) {
        Summary * summ = null;
        char[64] genline = 0;
        double xi = gretl_vector_get(xvals, i);
        double * rv = null;

        gretl_push_c_numeric_locale();
        sprintf(genline.ptr, "%s == %g", byname, xi);
        gretl_pop_c_numeric_locale();
        rv = generate_series(genline.ptr, dset, prn, &err);

        if (!err) {
            summ = get_summary_restricted(list, dset, rv,
                                          opt, prn, &err);
        }
        if (!err && !quiet) {
            if (i == 0) {
                pputc(prn, '\n');
            }
            if (st != null) {
                const(char) * si = series_table_get_string(st, xi);

                pprintf(prn, "%s (n = %d):", si, summ.n);
            } else {
                pprintf(prn, "%s = %g (n = %d):\n", byname, xi, summ.n);
            }
            print_summary(summ, dset, prn);
        }
        if (!err && m != null) {
            transcribe_summary(summ, m, &i, xi);
        }
        if (summ != null) {
            free_summary(summ);
        }
        free(rv);
    }

    if (!err && m != null) {
        set_last_result_data(m, GretlType.matrix);
    }

    gretl_matrix_free(xvals);
    if (mylist != null) {
        free(mylist);
    }

    return err;
}

/**
 * vmatrix_new:
 *
 * Returns: an allocated and initialized #VMatrix, or
 * %NULL on failure.
 */

VMatrix * vmatrix_new()
{
    VMatrix * v = cast(VMatrix *) malloc(VMatrix.sizeof);

    if (v != null) {
        v.vec = null;
        v.xbar = null;
        v.ssx = null;
        v.list = null;
        v.names = null;
        v.ci = 0;
        v.dim = 0;
        v.t1 = 0;
        v.t2 = 0;
        v.nmin = 0;
        v.nmax = 0;
        v.ncrit = 0;
        v.missing = 0;
    }

    return v;
}

/**
 * free_vmatrix:
 * @vmat: pointer to gretl correlation matrix struct
 *
 * Frees all resources associated with @vmat, and
 * the pointer itself.
 */

void free_vmatrix(VMatrix * vmat)
{
    if (vmat != null) {
        if (vmat.names != null) {
            strings_array_free(vmat.names, vmat.dim);
        }
        if (vmat.vec != null) {
            free(vmat.vec);
        }
        if (vmat.xbar != null) {
            free(vmat.xbar);
        }
        if (vmat.ssx != null) {
            free(vmat.ssx);
        }
        if (vmat.list != null) {
            free(vmat.list);
        }
        free(vmat);
    }
}

private enum CorrCov {
    corrmat = 0,
    covmat
}

/* Compute correlation or covariance matrix, using the maximum
   available sample for each coefficient.
*/

private int max_corrcov_matrix(VMatrix * v, const(DatasetRaw) * dset,
                               int flag)
{
    const(double **) Z = dset.Z;
    int i, j, vi, vj, idx;
    int vn = v.t2 - v.t1 + 1;
    int m = v.dim;
    int nij, nmiss;

    v.nmin = vn;
    v.nmax = 0;

    for (i=0; i<m; i++) {
        vi = v.list[i+1];
        for (j=i; j<m; j++)  {
            vj = v.list[j+1];
            idx = ijton(i, j, m);
            nij = vn;
            if (i == j && flag == CorrCov.corrmat) {
                v.vec[idx] = 1.0;
            } else {
                nmiss = 0;
                if (flag == CorrCov.covmat) {
                    v.vec[idx] = gretl_covar(v.t1, v.t2, Z[vi], Z[vj],
                                             &nmiss);
                } else {
                    v.vec[idx] = gretl_corr(v.t1, v.t2, Z[vi], Z[vj],
                                            &nmiss);
                }
                if (nmiss > 0) {
                    v.missing += 1;
                    nij -= nmiss;
                }
                if (nij > v.nmax) {
                    v.nmax = nij;
                }
                if (nij < v.nmin) {
                    v.nmin = nij;
                }
            }
        }
    }

    /* We'll record an "ncrit" value if there's something resembling
       a common number of observations across the coefficients.
       Specifically, we require that the difference between the max
       and min number of observations is less than 10 percent of the
       maximum; in that case we set "ncrit" to the minimum.
    */
    if (v.missing > 0) {
        v.ncrit = 0;
        if (v.nmax > 0) {
            double d = (v.nmax - v.nmin) / cast(double) v.nmax;

            if (d < 0.10) {
                v.ncrit = v.nmin;
            }
        }
    } else {
        v.ncrit = v.nmax;
    }

    return 0;
}

/* Compute correlation or covariance matrix, ensuring we use the same
   sample for all coefficients. We may be doing this in the context of
   the "corr" command or in the context of "pca". In the latter case
   we want to save the "xbar" and "ssx" values for later use when
   calculating the component series.
*/

private int uniform_corrcov_matrix(VMatrix * v, const(DatasetRaw) * dset,
                                   int ci, int flag)
{
    const(double **) Z = dset.Z;
    double * xbar = null;
    double * ssx = null;
    const(double) * x;
    int m = v.dim;
    int mm = (m * (m + 1)) / 2;
    double d1, d2;
    int i, j, jmin, nij, t;
    int miss, n = 0;

    xbar = cast(double *) malloc(m * (*xbar).sizeof);
    if (xbar == null) {
        return GretlError.alloc;
    }

    if (ci == GretlCmdIndex.pca || flag == CorrCov.corrmat) {
        ssx = cast(double *) malloc(m * (*ssx).sizeof);
        if (ssx == null) {
            free(xbar);
            return GretlError.alloc;
        }
    }

    for (i=0; i<m; i++) {
        xbar[i] = 0.0;
    }

    if (ssx != null) {
        for (i=0; i<m; i++) {
            ssx[i] = 0.0;
        }
    }

    v.missing = 0;

    /* first pass: get sample size and sums */

    for (t=v.t1; t<=v.t2; t++) {
        miss = 0;
        for (i=0; i<m; i++) {
            if (na(Z[v.list[i+1]][t])) {
                miss = 1;
                v.missing += 1;
                break;
            }
        }
        if (!miss) {
            n++;
            for (i=0; i<m; i++) {
                xbar[i] += Z[v.list[i+1]][t];
            }
        }
    }

    if (n < 2) {
        free(xbar);
        free(ssx);
        return GretlError.missdata;
    }

    for (i=0; i<m; i++) {
        xbar[i] /= n;
    }

    for (i=0; i<mm; i++) {
        v.vec[i] = 0.0;
    }

    /* second pass: get deviations from means and cumulate */

    for (t=v.t1; t<=v.t2; t++) {
        miss = 0;
        for (i=0; i<m; i++) {
            x = Z[v.list[i+1]];
            if (na(x[t])) {
                miss = 1;
                break;
            }
        }
        if (!miss) {
            for (i=0; i<m; i++) {
                x = Z[v.list[i+1]];
                d1 = x[t] - xbar[i];
                if (ssx != null) {
                    ssx[i] += d1 * d1;
                }
                jmin = (flag == CorrCov.covmat)? i : i + 1;
                for (j=jmin; j<m; j++) {
                    x = Z[v.list[j+1]];
                    nij = ijton(i, j, m);
                    d2 = x[t] - xbar[j];
                    v.vec[nij] += d1 * d2;
                }
            }
        }
    }

    /* finalize: compute correlations or covariances */

    if (flag == CorrCov.corrmat) {
        /* correlations */
        for (i=0; i<m; i++) {
            for (j=i; j<m; j++)  {
                nij = ijton(i, j, m);
                if (i == j) {
                    v.vec[nij] = 1.0;
                } else if (ssx[i] == 0.0 || ssx[j] == 0.0) {
                    v.vec[nij] = NADBL;
                } else {
                    v.vec[nij] /= sqrt(ssx[i] * ssx[j]);
                }
            }
        }
    } else {
        /* covariances */
        for (i=0; i<mm; i++) {
            v.vec[i] /= (n - 1);
        }
    }

    v.nmax = v.nmin = v.ncrit = n;

    if (ci == GretlCmdIndex.pca) {
        v.xbar = xbar;
        v.ssx = ssx;
    } else {
        free(xbar);
        free(ssx);
    }

    return 0;
}

private int corr_skip_obs(int k, int n_ok, int uniform)
{
    if (uniform) {
        return n_ok < k;
    } else {
        return n_ok < 2;
    }
}

private void corrlist_adjust_sample(const(int) * list,
                                    int * t1, int * t2,
                                    Gretlopt opt,
                                    const(DatasetRaw) * dset)
{
    int uniform = (opt & Gretlopt.OPT_N)? 1 : 0;
    int t1min = *t1;
    int t2max = *t2;
    int k = list[0];
    int n_ok;
    int i, t;

    /* Note: In this context what counts as a "usable" observation
       depends on whether the --uniform option (OPT_N) is given. If
       so, an observation is usable only if it's "complete" (no NAs);
       otherwise it's potentially usable so long as there at least two
       non-missing values among the members of @list.
    */

    /* advance start of sample range to skip any unusable obs */
    for (t=t1min; t<t2max; t++) {
        n_ok = k;
        for (i=1; i<=k; i++) {
            if (na(dset.Z[list[i]][t])) {
                n_ok--;
            }
        }
        if (corr_skip_obs(k, n_ok, uniform)) {
            t1min++;
        } else {
            break;
        }
    }

    /* retard end of sample range to skip any unusable obs */
    for (t=t2max; t>t1min; t--) {
        n_ok = k;
        for (i=1; i<=k; i++) {
            if (na(dset.Z[list[i]][t])) {
                n_ok--;
            }
        }
        if (corr_skip_obs(k, n_ok, uniform)) {
            t2max--;
        } else {
            break;
        }
    }

    *t1 = t1min;
    *t2 = t2max;
}

/**
 * corrlist:
 * @ci: command index.
 * @list: list of variables to process, by ID number.
 * @dset: dataset struct.
 * @opt: option flags.
 * @err: location to receive error code.
 *
 * Computes pairwise correlation coefficients for the variables
 * specified in @list, skipping any constants.  If the option
 * flags contain OPT_N, a uniform sample is ensured: only those
 * observations for which all the listed variables have valid
 * values are used.  If OPT_C is included, we actually calculate
 * covariances rather than correlations.
 *
 * Returns: gretl correlation matrix struct, or NULL on failure.
 */

VMatrix * corrlist(int ci, int * list, const(DatasetRaw) * dset,
                   Gretlopt opt, int * err)
{
    int flag = (opt & Gretlopt.OPT_C)? CorrCov.covmat : CorrCov.corrmat;
    VMatrix * v = null;
    int i, m, mm;

    if (list == null) {
        fprintf(stderr, "corrlist: list is NULL\n");
        *err = GretlError.data;
        return null;
    }

    v = vmatrix_new();
    if (v == null) {
        *err = GretlError.alloc;
        return null;
    }

    /* drop any constants from list */
    for (i=list[0]; i>0; i--) {
        if (gretl_isconst(dset.t1, dset.t2, dset.Z[list[i]])) {
            gretl_list_delete_at_pos(list, i);
        }
    }

    if (list[0] < 2) {
        gretl_errmsg_set(_("corr: needs at least two non-constant arguments").ptr);
        *err = GretlError.data;
        goto bailout;
    }

    v.t1 = dset.t1;
    v.t2 = dset.t2;

    /* adjust sample range if need be */
    corrlist_adjust_sample(list, &v.t1, &v.t2, opt, dset);

    if (v.t2 - v.t1 + 1 < 3) {
        *err = GretlError.toofew;
        goto bailout;
    }

    v.dim = m = list[0];
    mm = (m * (m + 1)) / 2;

    v.names = strings_array_new(m);
    v.vec = cast(double *) malloc(mm * (*v.vec).sizeof);
    v.list = gretl_list_copy(list);

    if (v.names == null || v.vec == null || v.list == null) {
        *err = GretlError.alloc;
        goto bailout;
    }

    if (ci == GretlCmdIndex.pca) {
        opt |= Gretlopt.OPT_N;
    }

    if (opt & Gretlopt.OPT_N) {
        /* impose uniform sample size */
        *err = uniform_corrcov_matrix(v, dset, ci, flag);
    } else {
        /* sample sizes may differ */
        *err = max_corrcov_matrix(v, dset, flag);
    }

    if (!*err) {
        for (i=0; i<m; i++) {
            v.names[i] = gretl_strdup(dset.varname[list[i+1]]);
            if (v.names[i] == null) {
                *err = GretlError.alloc;
                goto bailout;
            }
        }
    }

    v.ci = (flag == CorrCov.corrmat)? GretlCmdIndex.corr : 0; /* FIXME? */

 bailout:

    if (*err) {
        free_vmatrix(v);
        v = null;
    }

    return v;
}

private void printcorr(const(VMatrix) * v, Prn * prn)
{
    double r = v.vec[1];

    pprintf(prn, "\ncorr(%s, %s)", v.names[0], v.names[1]);

    if (na(r)) {
        pprintf(prn, ": %s\n\n", _("undefined").ptr);
    } else {
        pprintf(prn, " = %.8f\n", r);
        if (fabs(r) < 1.0) {
            int nm2 = v.ncrit - 2;
            double tval = r * sqrt(nm2 / (1 - r*r));

            pputs(prn, _("Under the null hypothesis of no correlation:\n ").ptr);
            pprintf(prn, _("t(%d) = %g, with two-tailed p-value %.4f\n").ptr,
                    nm2, tval, student_pvalue_2(nm2, tval));
            pputc(prn, '\n');
        } else {
            pprintf(prn, _("Two-tailed critical values for n = %d").ptr, v.ncrit);
            pprintf(prn, ": 5%% %.4f, 1%% %.4f\n", rhocrit(v.ncrit, 0.05),
                    rhocrit(v.ncrit, 0.01));
            pputc(prn, '\n');
        }
    }
}

/**
 * print_corrmat:
 * @corr: gretl correlation matrix.
 * @dset: dataset information.
 * @prn: gretl printing struct.
 *
 * Prints a gretl correlation matrix to @prn.
 */

void print_corrmat(VMatrix * corr, const(DatasetRaw) * dset, Prn * prn)
{
    if (corr.dim == 2) {
        printcorr(corr, prn);
    } else {
        char[OBSLEN] date1 = 0;
        char[OBSLEN] date2 = 0;
        gchar * tmp = null;

        ntolabel(date1.ptr, corr.t1, dset);
        ntolabel(date2.ptr, corr.t2, dset);

        pputc(prn, '\n');

        if (corr.nmin == corr.t2 - corr.t1 + 1) {
            tmp = g_strdup_printf(_("%s, using the %d observations %s - %s").ptr,
                                  _("Correlation Coefficients").ptr, corr.ncrit,
                                  date1.ptr, date2.ptr);
        } else if (corr.nmin == corr.nmax) {
            tmp = g_strdup_printf(_("%s, using %d observations from %s - %s").ptr,
                                  _("Correlation Coefficients").ptr, corr.ncrit,
                                  date1.ptr, date2.ptr);
        } else {
            tmp = g_strdup_printf(_("%s, using samples of size %d to %d").ptr,
                                  _("Correlation Coefficients").ptr,
                                  corr.nmin, corr.nmax);
        }
        output_line(tmp, prn, 0);
        g_free(tmp);

        if (corr.missing > 0) {
            output_line(_("(missing values were skipped)").ptr, prn, 1);
        }

        if (corr.ncrit > 0) {
            pprintf(prn, _("Two-tailed critical values for n = %d").ptr, corr.ncrit);
            pprintf(prn, ": 5%% %.4f, 1%% %.4f\n", rhocrit(corr.ncrit, 0.05),
                    rhocrit(corr.ncrit, 0.01));
            pputc(prn, '\n');
        }

        text_print_vmatrix(corr, prn);
    }
}

private void record_corr_matrix(VMatrix * c)
{
    GretlMatrixRaw * m = null;
    int i, j, k, n = c.dim;
    double cij;

    m = gretl_matrix_alloc(n, n);

    if (m != null) {
        k = 0;
        for (i=0; i<n; i++) {
            for (j=i; j<n; j++) {
                cij = c.vec[k];
                gretl_matrix_set(m, i, j, cij);
                gretl_matrix_set(m, j, i, cij);
                k++;
            }
        }

        if (c.names) {
            char ** rlab = strings_array_new(n);
            char ** clab = strings_array_new(n);

            if (rlab != null && clab != null) {
                for (i=0; i<n; i++) {
                    rlab[i] = gretl_strdup(c.names[i]);
                    clab[i] = gretl_strdup(c.names[i]);
                }
                gretl_matrix_set_rownames(m, rlab);
                gretl_matrix_set_colnames(m, clab);
            }
        }

        set_last_result_data(m, GretlType.matrix);
    }
}

/**
 * gretl_corrmx:
 * @list: gives the ID numbers of the variables to process.
 * @dset: dataset struct.
 * @opt: option flags: OPT_N means use uniform sample size,
 * OPT_U controls plotting options.
 * @prn: gretl printing struct.
 *
 * Computes and prints the correlation matrix for the specified list
 * of variables.
 *
 * Returns: 0 on successful completion, 1 on error.
 */

int gretl_corrmx(int * list, const(DatasetRaw) * dset,
                 Gretlopt opt, Prn * prn)
{
    VMatrix * corr = null;
    int err = 0;

    if (list != null) {
        if (list[0] == 0) {
            return 0;
        } else {
            corr = corrlist(GretlCmdIndex.corr, list, dset, opt, &err);
        }
    } else {
        int * tmplist = full_var_list(dset, null);

        if (tmplist != null) {
            corr = corrlist(GretlCmdIndex.corr, tmplist, dset, opt, &err);
            free(tmplist);
        }
    }

    if (corr != null) {
        if (!(opt & Gretlopt.OPT_Q)) {
            print_corrmat(corr, dset, prn);
        }
        if (corr.dim > 2 && gnuplot_graph_wanted(PlotType.heatmap, opt, null)) {
            err = plot_corrmat(corr, dset, opt);
        }

        record_corr_matrix(corr);
        free_vmatrix(corr);
    }

    return err;
}

/*
  Satterthwaite, F.E. (1946). "An Approximate Distribution of Estimates
  of Variance Components". Biometrics Bulletin, 2, 6, pp. 110–114.

  v = (s^2_1/n_1 + s^2_2/n_2)^2 /
        [(s^2_1/n1)^2 / (n_1 - 1) + (s^2_2/n_2)^2 / (n_2 - 1)]
*/

int satterthwaite_df(double v1, int n1,
                     double v2, int n2)
{
    double rtnum = v1 / n1 + v2 / n2;
    double d1 = (v1/n1) * (v1/n1) / (n1 - 1);
    double d2 = (v2/n2) * (v2/n2) / (n2 - 1);
    double v = (rtnum * rtnum) / (d1 + d2);

    return cast(int) floor(v);
}

private int test_data_from_dummy(const(DatasetRaw) * dset,
                                 int v1, int v2,
                                 double ** px,
                                 double ** py,
                                 int * pn1, int * pn2)
{
    double * x = null;
    double * y = null;
    int n1 = 0;
    int n2 = 0;
    int t, ix, iy;

    if (!gretl_isdummy(dset.t1, dset.t2, dset.Z[v2])) {
        return GretlError.invarg;
    }

    for (t=dset.t1; t<=dset.t2; t++) {
        if (na(dset.Z[v1][t])) {
            continue;
        } else if (dset.Z[v2][t] == 0.0) {
            n1++;
        } else if (dset.Z[v2][t] == 1.0) {
            n2++;
        }
    }
    if (n1 < 2 || n2 < 2) {
        return GretlError.toofew;
    }

    x = cast(double *) malloc(n1 * (*x).sizeof);
    y = cast(double *) malloc(n2 * (*y).sizeof);
    if (x == null || y == null) {
        free(x); free(y);
        return GretlError.alloc;
    }

    ix = iy = 0;
    for (t=dset.t1; t<=dset.t2; t++) {
        if (na(dset.Z[v1][t])) {
            continue;
        } else if (dset.Z[v2][t] == 0.0) {
            x[ix++] = dset.Z[v1][t];
        } else if (dset.Z[v2][t] == 1.0) {
            y[iy++] = dset.Z[v1][t];
        }
    }

    *px = x;
    *py = y;
    *pn1 = n1;
    *pn2 = n2;

    return 0;
}

private double * transcribe_ok_obs(const(DatasetRaw) * dset, int v,
                                   int * pn, int * err)
{
    const(double) * src = dset.Z[v];
    double * ret = null;
    int t, n = 0;

    for (t=dset.t1; t<=dset.t2; t++) {
        if (!na(src[t])) {
            n++;
        }
    }

    if (n < 2) {
        *err = GretlError.toofew;
    } else {
        ret = cast(double *) malloc(n * (*ret).sizeof);
        if (ret == null) {
            *err = GretlError.alloc;
        }
    }
    if (!*err) {
        *pn = n;
        n = 0;
        for (t=dset.t1; t<=dset.t2; t++) {
            if (!na(src[t])) {
                ret[n++] = src[t];
            }
        }
    }

    return ret;
}

private double * transcribe_ok_pairs(const(DatasetRaw) * dset,
                                     int v1, int v2,
                                     int * pn, int * err)
{
    const(double) * x = dset.Z[v1];
    const(double) * y = dset.Z[v2];
    double * ret = null;
    int t, n = 0;

    for (t=dset.t1; t<=dset.t2; t++) {
        if (!na(x[t]) && !na(y[t])) {
            n++;
        }
    }

    if (n < 2) {
        *err = GretlError.toofew;
    } else {
        ret = cast(double *) malloc(n * (*ret).sizeof);
        if (ret == null) {
            *err = GretlError.alloc;
        }
    }
    if (!*err) {
        *pn = n;
        n = 0;
        for (t=dset.t1; t<=dset.t2; t++) {
            if (!na(x[t]) && !na(y[t])) {
                ret[n++] = x[t] - y[t];
            }
        }
    }

    return ret;
}

/* HAC difference of means test using OLS on stacked data */

private int ols_means_test(const(DatasetRaw) * dset,
                           const(double) * x, int n1,
                           const(double) * y, int n2,
                           double * m1, double * m2,
                           double * se, int * df,
                           Prn * prn)
{
    Model mod;
    DatasetRaw * aset;
    int[4] list = [3, 1, 0, 2];
    int nobs = n1 + n2;
    int i, err = 0;

    aset = create_auxiliary_dataset(list[0], nobs, Gretlopt.OPT_NONE);
    if (aset == null) {
        return GretlError.alloc;
    }
    aset.pd = dset.pd;
    aset.structure = dset.structure;

    /* stacked dependent variable, @x | @y */
    memcpy(aset.Z[1], x, n1 * (*x).sizeof);
    memcpy(aset.Z[1] + n1, y, n2 * (*y).sizeof);
    strcpy(aset.varname[1], "x_y");

    /* dummy to select the @y observations */
    for (i=0; i<nobs; i++) {
        aset.Z[2][i] = i >= n1;
    }

    gretl_model_init(&mod, aset);
    mod = lsq(list.ptr, aset, GretlCmdIndex.ols, Gretlopt.OPT_A | Gretlopt.OPT_R);
    err = mod.errcode;

    if (!err) {
        if (prn != null) {
            printmodel(&mod, aset, Gretlopt.OPT_S, prn);
        }
        *m1 = mod.coeff[0];
        *m2 = mod.coeff[0] + mod.coeff[1];
        *se = mod.sderr[1];
        *df = mod.dfd;
    }

    clear_model(&mod);
    destroy_dataset(aset);

    return err;
}

private int test_data_from_list(const(DatasetRaw) * dset,
                                int v1, int v2,
                                double ** px,
                                double ** py,
                                int * pn1, int * pn2)
{
    double * x;
    double * y;
    int n1, n2;
    int err = 0;

    x = transcribe_ok_obs(dset, v1, &n1, &err);
    if (!err) {
        y = transcribe_ok_obs(dset, v2, &n2, &err);
    }
    if (!err) {
        *px = x;
        *py = y;
        *pn1 = n1;
        *pn2 = n2;
    }

    return err;
}

private int paired_data_from_list(const(DatasetRaw) * dset,
                                  int v1, int v2,
                                  double ** px, int * pn)
{
    double * x;
    int n, err = 0;

    x = transcribe_ok_pairs(dset, v1, v2, &n, &err);
    if (!err) {
        *px = x;
        *pn = n;
    }

    return err;
}

private void paired_diff_print(int n,
                               double mdiff,
                               double se,
                               int df,
                               double t,
                               double pval,
                               Prn * prn)
{
    pputc(prn, '\n');
    pputs(prn, _("Paired difference test").ptr);
    pputs(prn, "\n\n");
    pputs(prn, "   ");
    pprintf(prn, _("Number of observations = %d\n").ptr, n);
    pputs(prn, "   ");
    pprintf(prn, _("Sample mean paired difference = %g").ptr, mdiff);
    pputc(prn, '\n');
    pputs(prn, _("   Null hypothesis: The population mean difference is zero.").ptr);
    pputc(prn, '\n');
    pprintf(prn, _("   Estimated standard error = %g\n").ptr, se);
    pprintf(prn, _("   Test statistic: t(%d) = %g\n").ptr, df, t);
    pprintf(prn, _("   p-value (two-tailed) = %g\n\n").ptr, pval);
}

private int splitby_varnum(GretlCmdIndex ci, const(DatasetRaw) * dset)
{
    const(char) * vname = get_optval_string(ci, Gretlopt.OPT_D);

    return current_series_index(dset, vname);
}

/**
 * diff_test:
 * @list: list containing 2 variables.
 * @dset: dataset struct.
 * @opt: %OPT_G, sign test; %OPT_R, rank sum; %OPT_I,
 * signed rank; %OPT_V, verbose (for rank tests).
 * %OPT_D (--split-by) requires OPT_R.
 * @prn: gretl printing struct.
 *
 * Performs, and prints the results of, a non-parametric
 * test for a difference between two variables or groups.
 * The specific test performed depends on @opt, the sign
 * test being the default.
 *
 * Returns: 0 on successful completion, non-zero on error.
 */

int diff_test(const(int) * list, const(DatasetRaw) * dset,
              Gretlopt opt, Prn * prn)
{
    double * x = null;
    double * y = null;
    double[2] result = 0;
    int usedum = (opt & Gretlopt.OPT_D)? 1 : 0;
    int v1, v2;
    int n1 = 0;
    int n2 = 0;
    int err = 0;

    if (usedum && !(opt & Gretlopt.OPT_R)) {
        /* OPT_D is supported only in conjunction with OPT_R */
        return GretlError.badopt;
    }

    err = incompatible_options(opt, Gretlopt.OPT_G | Gretlopt.OPT_R | Gretlopt.OPT_I);
    if (err) {
        return err;
    }

    v1 = list[1];
    if (usedum) {
        v2 = splitby_varnum(GretlCmdIndex.difftest, dset);
        if (v2 < 0) {
            err = GretlError.unkvar;
        }
    } else {
        v2 = list[2];
    }

    if (err) {
        return err;
    }

    if (usedum) {
        /* get @x and @y as portions of series v1 */
        err = test_data_from_dummy(dset, v1, v2, &x, &y, &n1, &n2);
    } else {
        /* get @x and @y from the series v1 and v2 */
        err = test_data_from_list(dset, v1, v2, &x, &y, &n1, &n2);
    }

    if (opt & Gretlopt.OPT_R) {
        err = rank_sum_test(x, y, v1, v2, dset, result.ptr, opt, prn);
    } else if (opt & Gretlopt.OPT_I) {
        err = signed_rank_test(x, y, v1, v2, dset, result.ptr, opt, prn);
    } else {
        err = sign_test(x, y, v1, v2, dset, result.ptr, opt, prn);
    }

    record_test_result(result[0], result[1]);

    return err;
}

/**
 * means_test:
 * @list: gives the ID numbers of the series to compare.
 * @dset: dataset struct.
 * @opt: if OPT_O, assume that the population variances differ;
 * if OPT_D (--split-by) interpret the second series as a factor (dummy);
 * if OPT_R, carry out a robust, regression-based test; if OPT_P
 * (--paired) carry out a paired-difference test; if OPT_Q, don't print
 * the resuults.
 * @prn: gretl printing struct.
 *
 * Carries out test of the null hypothesis that the means of two
 * series are equal.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int means_test(const(int) * list, const(DatasetRaw) * dset,
               Gretlopt opt, Prn * prn)
{
    double m1, m2, s1, s2, se, mdiff, t, pval;
    double * x = null;
    double * y = null;
    int vardiff = (opt & Gretlopt.OPT_O)? 1 : 0;
    int usedum = (opt & Gretlopt.OPT_D)? 1 : 0;
    int robust = (opt & Gretlopt.OPT_R)? 1 : 0;
    int paired = (opt & Gretlopt.OPT_P)? 1 : 0;
    int v1, v2;
    int n1 = 0;
    int n2 = 0;
    int df, err;

    if (paired && usedum) {
        /* --split-by and --paired are incompatible */
        return GretlError.badopt;
    } else if (!usedum && list[0] == 1) {
        /* We need two series */
        return GretlError.args;
    }

    err = incompatible_options(opt, Gretlopt.OPT_O | Gretlopt.OPT_P | Gretlopt.OPT_R);
    if (!err && (opt & Gretlopt.OPT_R) && !dataset_is_time_series(dset)) {
        gretl_errmsg_set("meantest: the --robust option requires time-series data");
        err = GretlError.tsonly;
    }

    if (!err) {
        v1 = list[1];
        if (usedum) {
            v2 = splitby_varnum(GretlCmdIndex.meantest, dset);
            if (v2 < 0) {
                err = GretlError.unkvar;
            }
        } else {
            v2 = list[2];
        }
    }

    if (err) {
        return err;
    }

    if (paired) {
        /* fill @x with paired differences between series v1 and v2 */
        err = paired_data_from_list(dset, v1, v2, &x, &n1);
    } else if (usedum) {
        /* get @x and @y as portions of series v1 */
        err = test_data_from_dummy(dset, v1, v2, &x, &y, &n1, &n2);
    } else {
        /* get @x and @y from the series v1 and v2 */
        err = test_data_from_list(dset, v1, v2, &x, &y, &n1, &n2);
    }

    if (err) {
        free(x); free(y);
        return err;
    }

    if (paired) {
        /* @x actually holds the x,y differences */
        gretl_moments(0, n1-1, x, null, &m1, &se, null, null, 1);
        se /= sqrt(cast(double) n1);
        mdiff = m1;
        df = n1 - 1;
    } else if (robust) {
        /* use regression method */
        err = ols_means_test(dset, x, n1, y, n2, &m1, &m2, &se, &df, prn);
        if (err) {
            goto bailout;
        }
        mdiff = m1 - m2;
    } else {
        double var1, var2;

        df = n1 + n2 - 2;
        gretl_moments(0, n1-1, x, null, &m1, &s1, null, null, 1);
        gretl_moments(0, n2-1, y, null, &m2, &s2, null, null, 1);
        mdiff = m1 - m2;
        var1 = s1 * s1;
        var2 = s2 * s2;

        if (vardiff) {
            /* assuming unequal variances */
            se = sqrt((var1 / n1) + (var2 / n2));
            df = satterthwaite_df(var1, n1, var2, n2);
        } else {
            /* form pooled estimate of variance */
            double s2p;

            s2p = ((n1-1) * var1 + (n2-1) * var2) / df;
            se = sqrt(s2p / n1 + s2p / n2);
            df = n1 + n2 - 2;
        }
    }

    t = mdiff / se;
    pval = student_pvalue_2(df, t);
    record_test_result(t, pval);

    if (opt & Gretlopt.OPT_Q) {
        /* quiet */
        goto bailout;
    }

    if (paired) {
        paired_diff_print(n1, mdiff, se, df, t, pval, prn);
    } else {
        if (robust) {
            pputs(prn, _("Robust equality of means test").ptr);
            pputs(prn, "\n\n");
        } else {
            pprintf(prn, _("\nEquality of means test "
                           ~ "(assuming %s variances)\n\n").ptr,
                    vardiff ? _("unequal").ptr : _("equal").ptr);
        }
        if (usedum) {
            pprintf(prn, "   %s (%s==0): ", dset.varname[v1], dset.varname[v2]);
            pprintf(prn, _("Number of observations = %d\n").ptr, n1);
            pprintf(prn, "   %s (%s==1): ", dset.varname[v1], dset.varname[v2]);
            pprintf(prn, _("Number of observations = %d\n").ptr, n2);
        } else {
            pprintf(prn, "   %s: ", dset.varname[v1]);
            pprintf(prn, _("Number of observations = %d\n").ptr, n1);
            pprintf(prn, "   %s: ", dset.varname[v2]);
            pprintf(prn, _("Number of observations = %d\n").ptr, n2);
        }
        pprintf(prn, _("   Difference between sample means = %g - %g = %g\n").ptr,
                m1, m2, mdiff);
        pputs(prn, _("   Null hypothesis: The two population means are the same.\n").ptr);
        pprintf(prn, _("   Estimated standard error = %g\n").ptr, se);
        pprintf(prn, _("   Test statistic: t(%d) = %g\n").ptr, df, t);
        pprintf(prn, _("   p-value (two-tailed) = %g\n\n").ptr, pval);
    }

 bailout:

    free(x);
    free(y);

    return 0;
}

private double trimmed_mean(const(double) * x, double p,
                            int n, int * err)
{
    double ret;
    double * sx;
    int skip;

    sx = cast(double *) malloc(n * (*sx).sizeof);
    memcpy(sx, x, n * (*sx).sizeof);
    qsort(sx, n, (*sx).sizeof, &gretl_compare_doubles);
    skip = cast(int) floor(p * n);
    ret = gretl_mean(skip, n-skip-1, sx);
    free(sx);

    return ret;
}

private double levene_center(const(double) * x, int n,
                             int f, double p,
                             int * err)
{
    if (f == 1) {
        return gretl_mean(0, n-1, x);
    } else if (f == 2) {
        return gretl_median(0, n-1, x);
    } else {
        return trimmed_mean(x, p, n, err);
    }
}

/* Levene's test for difference of variance, with the options added by
   Brown and Forsythe. Refs:

   Levene, H. (1960) Robust Tests for Equality of Variances. In: Olkin,
   I. (ed), Contributions to Probability and Statistics, Stanford
   University Press. pp. 278-292.

   Brown, Morton B and Forsythe, Alan B. (1974). "Robust tests for the
   equality of variances". Journal of the American Statistical
   Association 69 (346), pp. 364–367.
*/

private double levene_test(const(double) * x, int n1,
                           const(double) * y, int n2,
                           int f, double p,
                           int * err)
{
    double zxbar, zybar, zbar, lc;
    double num, den, W;
    double d1, d2;
    double * zx;
    double * zy;
    int i, dfd;

    zx = cast(double *) malloc(n1 * (*zx).sizeof);
    zy = cast(double *) malloc(n2 * (*zy).sizeof);
    if (zx == null || zy == null) {
        *err = GretlError.alloc;
        return NADBL;
    }

    /* absolute deviations, x */
    lc = levene_center(x, n1, f, p, err);
    zxbar = 0;
    for (i=0; i<n1; i++) {
        zx[i] = fabs(x[i] - lc);
        zxbar += zx[i];
    }
    zxbar /= n1;

    /* absolute deviations, y */
    lc = levene_center(y, n2, f, p, err);
    zybar = 0;
    for (i=0; i<n2; i++) {
        zy[i] = fabs(y[i] - lc);
        zybar += zy[i];
    }
    zybar /= n2;

    /* grand mean */
    zbar = (n1 * zxbar + n2 * zybar) / (n1+n2);
    dfd = n1 + n2 - 2;

    /* numerator of W */
    d1 = zxbar - zbar;
    d2 = zybar - zbar;
    num = n1 * d1 * d1 + n2 * d2 * d2;

    /* denominator of W */
    den = 0;
    for (i=0; i<n1; i++) {
        d1 = zx[i] - zxbar;
        den += d1 * d1;
    }
    for (i=0; i<n2; i++) {
        d2 = zy[i] - zybar;
        den += d2 * d2;
    }

    W = dfd * num / den;

    free(zx);
    free(zy);

    return W;
}

private int handle_levene_options(double * p, int * err)
{
    const(char) * s = get_optval_string(GretlCmdIndex.vartest, Gretlopt.OPT_R);
    double default_trim = 0.10;
    int f = 1;

    if (s != null) {
        if (!strcmp(s, "mean")) {
            /* OK */
        } else if (!strcmp(s, "median")) {
            f = 2;
        } else if (!strcmp(s, "trimmed")) {
            f = 3;
            *p = default_trim;
        } else if (!strncmp(s, "trimmed,", 8)) {
            int k = positive_int_from_string(s + 8);

            if (k > 0 && k < 50) {
                *p = k / 100.0;
                f = 3;
            } else {
                *err = GretlError.invarg;
            }
        }
    }

    return f;
}

private void levene_print(double W, double pval,
                          int dfn, int dfd, int f,
                          double p, Prn * prn)
{
    pputs(prn, _("\nRobust equality of variances test\n").ptr);
    pputs(prn, _("Centering: ").ptr);
    if (f == 1) {
        pputs(prn, _("mean").ptr);
    } else if (f == 2) {
        pputs(prn, _("median").ptr);
    } else {
        pprintf(prn, _("%d percent trimmed mean").ptr, cast(int) (p*100));
    }
    pputc(prn, '\n');
    pprintf(prn, "F(%d,%d) = %g ", dfn, dfd, W);
    pprintf(prn, "%s %g\n", _("with p-value").ptr, pval);
}

private void vartest_print(double F, double pval,
                           int dfn, int dfd,
                           const(DatasetRaw) * dset,
                           int v1, int v2,
                           int n1, int n2,
                           int usedum, Prn * prn)
{
    pputs(prn, _("\nEquality of variances test\n\n").ptr);
    if (usedum) {
        pprintf(prn, "   %s (%s==0): ", dset.varname[v1], dset.varname[v2]);
        pprintf(prn, _("Number of observations = %d\n").ptr, n1);
        pprintf(prn, "   %s (%s==1): ", dset.varname[v1], dset.varname[v2]);
        pprintf(prn, _("Number of observations = %d\n").ptr, n2);
    } else {
        pprintf(prn, "   %s: ", dset.varname[v1]);
        pprintf(prn, _("Number of observations = %d\n").ptr, n1);
        pprintf(prn, "   %s: ", dset.varname[v2]);
        pprintf(prn, _("Number of observations = %d\n").ptr, n2);
    }
    pprintf(prn, _("   Ratio of sample variances = %g\n").ptr, F);
    pprintf(prn, "   %s: %s\n", _("Null hypothesis").ptr,
            _("The two population variances are equal").ptr);
    pprintf(prn, "   %s: F(%d,%d) = %g\n", _("Test statistic").ptr, dfn, dfd, F);
    pprintf(prn, _("   p-value (two-tailed) = %g\n\n").ptr, pval);
}

/**
 * vars_test:
 * @list: gives the ID numbers of the variables to compare.
 * @dset: pointer to dataset.
 * @opt: if includes OPT_R, carry out a robust (Levene-type) test;
 * if OPT_D, interpret the second series as a factor (dummy); if
 * OPT_Q (quiet), don't print the results.
 * @prn: gretl printing struct.
 *
 * Carries out test of the null hypothesis that the variances of two
 * variables are equal.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int vars_test(const(int) * list, const(DatasetRaw) * dset,
              Gretlopt opt, Prn * prn)
{
    double F, pval;
    double p = NADBL;
    double * x = null;
    double * y = null;
    int robust = (opt & Gretlopt.OPT_R)? 1 : 0;
    int usedum = (opt & Gretlopt.OPT_D)? 1 : 0;
    int dfn, dfd;
    int v1, v2;
    int n1 = 0;
    int n2 = 0;
    int f = 0;
    int err = 0;

    v1 = list[1];
    if (usedum) {
        v2 = splitby_varnum(GretlCmdIndex.vartest, dset);
        if (v2 < 0) {
            err = GretlError.unkvar;
        }
    } else {
        v2 = list[2];
    }

    if (!err && robust) {
        f = handle_levene_options(&p, &err);
        if (err) {
            return err;
        }
    }

    if (usedum) {
        /* get @x and @y as portions of dset->Z[v1] */
        err = test_data_from_dummy(dset, v1, v2, &x, &y, &n1, &n2);
    } else {
        /* get @x and @y as dset->Z[v1], dset->Z[v2] */
        err = test_data_from_list(dset, v1, v2, &x, &y, &n1, &n2);
    }

    if (err) {
        free(x); free(y);
        return err;
    }

    if (robust) {
        F = levene_test(x, n1, y, n2, f, p, &err);
        dfn = 1; /* = 2 - 1 */
        dfd = n1 + n2 - 2;
        pval = snedecor_cdf_comp(dfn, dfd, F);
    } else {
        double mx, my, sx, sy;

        gretl_moments(0, n1-1, x, null, &mx, &sx, null, null, 1);
        gretl_moments(0, n2-1, y, null, &my, &sy, null, null, 1);
        F = (sx * sx) / (sy * sy);
        dfn = n1 - 1;
        dfd = n2 - 1;
        pval = snedecor_pvalue_2(dfn, dfd, F);
    }

    pval = snedecor_cdf_comp(dfn, dfd, F);
    record_test_result(F, pval);

    if (!(opt & Gretlopt.OPT_Q)) {
        /* not quiet */
        if (robust) {
            levene_print(F, pval, dfn, dfd, f, p, prn);
        } else {
            vartest_print(F, pval, dfn, dfd, dset, v1, v2,
                          n1, n2, usedum, prn);
        }
    }

    free(x);
    free(y);

    return err;
}

struct MahalDist {
    int * list;
    int n;
    GretlMatrixRaw * d;
}

const(double) * mahal_dist_get_distances(const(MahalDist) * md)
{
    return md.d.val;
}

int mahal_dist_get_n(const(MahalDist) * md)
{
    return md.n;
}

const(int) * mahal_dist_get_varlist(const(MahalDist) * md)
{
    return md.list;
}

void free_mahal_dist(MahalDist * md)
{
    if (md != null) {
        free(md.list);
        gretl_matrix_free(md.d);
        free(md);
    }
}

private MahalDist * mahal_dist_new(const(int) * list, int n)
{
    MahalDist * md = cast(MahalDist *) malloc(MahalDist.sizeof);

    if (md != null) {
        md.d = gretl_column_vector_alloc(n);
        if (md.d == null) {
            free(md);
            md = null;
        } else {
            md.list = gretl_list_copy(list);
            if (md.list == null) {
                gretl_matrix_free(md.d);
                free(md);
                md = null;
            } else {
                md.n = n;
            }
        }
    }

    if (md != null) {
        int t;

        for (t=0; t<n; t++) {
            md.d.val[t] = NADBL;
        }
    }

    return md;
}

private int mdist_saver(DatasetRaw * dset)
{
    int t, v, err;

    err = dataset_add_series(dset, 1);

    if (err) {
        v = 0;
    } else {
        v = dset.v - 1;

        for (t=0; t<dset.n; t++) {
            dset.Z[v][t] = NADBL;
        }

        strcpy(dset.varname[v], "mdist");
        make_varname_unique(dset.varname[v], v, dset);
        series_set_label(dset, v, _("Mahalanobis distances").ptr);
    }

    return v;
}

private int real_mahalanobis_distance(const(int) * list, DatasetRaw * dset,
                                      Gretlopt opt, MahalDist * md,
                                      Prn * prn)
{
    GretlMatrixRaw * S = null;
    GretlVector * means = null;
    GretlVector * xdiff;
    int orig_t1 = dset.t1;
    int orig_t2 = dset.t2;
    int n, err = 0;

    list_adjust_sample(list, &dset.t1, &dset.t2, dset, null);

    n = sample_size(dset);
    if (n < 2) {
        err = GretlError.data;
        goto bailout;
    }

    xdiff = gretl_column_vector_alloc(list[0]);
    if (xdiff == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    S = gretl_covariance_matrix_from_varlist(list,
                                             dset,
                                             &means,
                                             &err);

    if (!err) {
        if (opt & Gretlopt.OPT_V) {
            gretl_matrix_print_to_prn(S, _("Covariance matrix").ptr, prn);
        }
        err = gretl_invert_symmetric_matrix(S);
        if (err) {
            fprintf(stderr, "error inverting covariance matrix\n");
        } else if (opt & Gretlopt.OPT_V) {
            gretl_matrix_print_to_prn(S, _("Inverse of covariance matrix").ptr, prn);
        }
    }

    if (!err) {
        int k = gretl_vector_get_length(means);
        int miss, obslen = 0, savevar = 0;
        double m, x, xbar;
        int i, t, vi;

        if (opt & Gretlopt.OPT_S) {
            /* save the results to a data series */
            savevar = mdist_saver(dset);
        }

        if (prn != null) {
            pprintf(prn, "%s\n", _("Mahalanobis distances from the centroid").ptr);
            pprintf(prn, "%s\n", _("using the variables:").ptr);
            for (i=1; i<=list[0]; i++) {
                pprintf(prn, " %s\n", dset.varname[list[i]]);
            }
            pputc(prn, '\n');

            obslen = max_obs_marker_length(dset);
            if (obslen < 8) {
                obslen = 8;
            }
        }

        for (t=dset.t1; t<=dset.t2; t++) {
            miss = 0;

            /* write vector of deviations from centroid for
               observation t */
            for (i=0; i<k; i++) {
                vi = list[i+1];
                xbar = gretl_vector_get(means, i);
                x = dset.Z[vi][t];
                miss += na(x);
                if (!miss) {
                    gretl_vector_set(xdiff, i,  x - xbar);
                }
            }

            if (miss) {
                m = NADBL;
            } else {
                m = gretl_scalar_qform(xdiff, S, &err);
                if (!err) {
                    m = sqrt(m);
                }
            }

            if (prn != null) {
                print_obs_marker(t, dset, obslen, prn);
                if (err || miss) {
                    pprintf(prn, "NA\n");
                } else {
                    pprintf(prn, "%9.6f\n", m);
                }
            }

            if (savevar > 0) {
                dset.Z[savevar][t] = m;
            } else if (md != null) {
                md.d.val[t] = m;
            }
        }

        if (savevar > 0 && prn != null) {
            pputc(prn, '\n');
            pprintf(prn, _("Distances saved as '%s'").ptr,
                    dset.varname[savevar]);
            pputc(prn, '\n');
        }

        if (prn != null) {
            pputc(prn, '\n');
        }
    }

    gretl_matrix_free(xdiff);
    gretl_matrix_free(means);
    gretl_matrix_free(S);

 bailout:

    dset.t1 = orig_t1;
    dset.t2 = orig_t2;

    return err;
}

int mahalanobis_distance(const(int) * list, DatasetRaw * dset,
                         Gretlopt opt, Prn * prn)
{
    return real_mahalanobis_distance(list, dset, opt, null,
                                     (opt & Gretlopt.OPT_Q) ? null : prn);
}

/* for use in the GUI -- see gui/library.c */

MahalDist * get_mahal_distances(const(int) * list, DatasetRaw * dset,
                                Gretlopt opt, Prn * prn, int * err)
{
    MahalDist * md = mahal_dist_new(list, dset.n);

    if (md == null) {
        *err = GretlError.alloc;
    } else {
        *err = real_mahalanobis_distance(list, dset, opt,
                                         md, prn);
        if (*err) {
            free_mahal_dist(md);
            md = null;
        }
    }

    return md;
}

private enum DistanceMetric {
    euclidean,
    manhattan,
    hamming,
    chebyshev,
    cosine,
    mahalanobis
}

private int distance_metric(const(char) * s)
{
    int n = cast(int) strlen(s);

    if (n == 0 || (n == 1 && (*s == 'c' || *s == 'm'))) {
        return -1;
    }

    if (!strncmp(s, "euclidean", n)) {
        return DistanceMetric.euclidean;
    } else if (!strncmp(s, "manhattan", n)) {
        return DistanceMetric.manhattan;
    } else if (!strncmp(s, "hamming", n)) {
        return DistanceMetric.hamming;
    } else if (!strncmp(s, "chebyshev", n) ||
               !strncmp(s, "chebychev", n)) {
        return DistanceMetric.chebyshev;
    } else if (!strncmp(s, "cosine", n)) {
        return DistanceMetric.cosine;
    } else if (!strncmp(s, "mahalanobis", n)) {
        return DistanceMetric.mahalanobis;
    } else {
        return -1;
    }
}

private GretlMatrixRaw * mahal_cholesky_factor(const(GretlMatrixRaw) * X,
                                            int * err)
{
    GretlMatrixRaw * C;

    C = gretl_covariance_matrix(X, 0, 1, err);
    if (!*err) {
        *err = gretl_matrix_cholesky_decomp(C);
        if (!*err) {
            gretl_square_matrix_transpose(C);
            *err = gretl_invert_triangular_matrix(C, 'U');
        }
    }

    if (*err) {
        gretl_matrix_free(C);
        C = null;
    }

    return C;
}

/* Produce scaled versions of @X and @Y using the inverse Cholesky
   factor of @X */

private int mahalanobis_prep(const(GretlMatrixRaw) * X,
                             const(GretlMatrixRaw) * Y,
                             GretlMatrixRaw ** XC,
                             GretlMatrixRaw ** YC)
{
    GretlMatrixRaw * C;
    int err = 0;

    C = mahal_cholesky_factor(X, &err);
    if (err) {
        return err;
    }

    *XC = gretl_matrix_alloc(X.rows, X.cols);
    if (*XC == null) {
        err = GretlError.alloc;
    } else {
        gretl_matrix_multiply(X, C, *XC);
    }

    if (!err && Y != null) {
        *YC = gretl_matrix_alloc(Y.rows, Y.cols);
        if (*YC == null) {
            gretl_matrix_free(*XC);
            err = GretlError.alloc;
        } else {
            gretl_matrix_multiply(Y, C, *YC);
        }
    }

    gretl_matrix_free(C);

    return err;
}

/* distance(): produces a set of pairwise distances, either between
   the rows of @X, or if @Y is non-NULL, between the rows of @X and
   the rows of @Y.

   If @X is m x n and @Y is NULL, the result is a vector @v of length
   m * (m - 1) / 2, comprising the non-redundant pairwise distances
   with the zeros on the diagonal suppressed. This can be turned into
   a square matrix via unvech(v, 0).

   If @X is m x n and @Y is p x n, the result is a matrix with m rows
   and p columns.

   If @X and @Y have different column dimensions an error is flagged.
*/

GretlMatrixRaw * distance(const(GretlMatrixRaw) * X,
                       const(GretlMatrixRaw) * Y,
                       const(char) * type,
                       int * err)
{
    GretlMatrixRaw * ret;
    const(GretlMatrixRaw) * P = X;
    const(GretlMatrixRaw) * Q = Y;
    GretlMatrixRaw * XC = null;
    GretlMatrixRaw * YC = null;
    double d, dij, x, y;
    double den1, den2;
    int metric, pos;
    int m, n, p, jmin;
    int i, j, k;
    int colvec = (Y == null);

    if (gretl_is_null_matrix(X)) {
        fprintf(stderr, "distance(): X is null matrix\n");
        *err = GretlError.invarg;
    } else if (gretl_is_complex(X) || gretl_is_complex(Y)) {
        *err = GretlError.cmplx;
    } else if (Y != null && Y.cols != X.cols) {
        fprintf(stderr, "distance(): Y is not conformable\n");
        *err = GretlError.nonconf;
    } else if ((metric = distance_metric(type)) < 0) {
        fprintf(stderr, "distance(): invalid metric value '%s'\n", type);
        *err = GretlError.invarg;
    }

    if (*err) {
        return null;
    }

    m = X.rows;
    n = X.cols;

    if (colvec) {
        p = m;
        /* column vector for results */
        ret = gretl_matrix_alloc(m * (m-1) / 2, 1);
    } else {
        p = Y.rows;
        /* matrix for results */
        ret = gretl_matrix_alloc(m, p);
    }

    if (ret == null) {
        *err = GretlError.alloc;
        return null;
    }

    if (metric == DistanceMetric.mahalanobis) {
        *err = mahalanobis_prep(X, Y, &XC, &YC);
        if (*err) {
            free(ret);
            return null;
        }
        P = XC;
        Q = YC;
        metric = DistanceMetric.euclidean;
    }

    pos = 0;
    for (i=0; i<m; i++) {
        jmin = colvec ? i + 1 : 0;
        for (j=jmin; j<p; j++) {
            den1 = den2 = dij = 0;
            for (k=0; k<n; k++) {
                x = gretl_matrix_get(P, i, k);
                if (colvec) {
                    y = gretl_matrix_get(P, j, k);
                } else {
                    y = gretl_matrix_get(Q, j, k);
                }
                if (metric == DistanceMetric.manhattan) {
                    dij += fabs(x - y);
                } else if (metric == DistanceMetric.hamming) {
                    dij += (x - y != 0);
                } else if (metric == DistanceMetric.chebyshev) {
                    d = fabs(x - y);
                    if (d > dij) {
                        dij = d;
                    }
                } else if (metric == DistanceMetric.cosine) {
                    dij  += x * y;
                    den1 += x * x;
                    den2 += y * y;
                } else {
                    /* Euclidean or Mahalanobis */
                    d = x - y;
                    dij += d * d;
                }
            }
            if (metric == DistanceMetric.euclidean) {
                dij = sqrt(dij);
            } else if (metric == DistanceMetric.hamming) {
                dij /= n;
            } else if (metric == DistanceMetric.cosine) {
                dij = 1.0 - dij / sqrt(den1 * den2);
            }
            if (colvec) {
                ret.val[pos++] = dij;
            } else {
                gretl_matrix_set(ret, i, j, dij);
            }
        }
    }

    gretl_matrix_free(XC);
    gretl_matrix_free(YC);

    return ret;
}

/*
   G = [(2 * \sum_i^n (i * x_i)) / (n * \sum_i^n x_i )] - [(n + 1)/n]

   where x is sorted: x_i <= x_{i+1}
*/

private double gini_coeff(const(double) * x, int t1, int t2, double ** plz,
                          int * pn, int * err)
{
    int m = t2 - t1 + 1;
    double * sx = null;
    double csx = 0.0, sumx = 0.0, sisx = 0.0;
    double G;
    int t, n = 0;

    gretl_error_clear();

    sx = cast(double *) malloc(m * (*sx).sizeof);

    if (sx == null) {
        *err = GretlError.alloc;
        return NADBL;
    }

    n = 0;
    for (t=t1; t<=t2; t++) {
        if (na(x[t])) {
            continue;
        } else if (x[t] < 0.0) {
            gretl_errmsg_set(_("illegal negative value").ptr);
            *err = GretlError.data;
            break;
        } else {
            sx[n++] = x[t];
            sumx += x[t];
        }
    }

    if (!*err && (n == 0 || sumx == 0.0)) {
        *err = GretlError.data;
    }

    if (*err) {
        free(sx);
        return NADBL;
    }

    qsort(sx, n, (*sx).sizeof, &gretl_compare_doubles);

    for (t=0; t<n; t++) {
        csx += sx[t];
        sisx += (t + 1) * sx[t];
        sx[t] = csx / sumx; /* sx now = Lorenz curve */
    }

    G = 2.0 * sisx / (n * sumx) - (cast(double) n + 1) / n;

    if (plz != null) {
        *plz = sx;
        *pn = n;
    } else {
        free(sx);
    }

    return G;
}

/**
 * gretl_gini:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series.
 *
 * Returns: the Gini coefficient for the series @x from obs
 * @t1 to obs @t2, skipping any missing values, or #NADBL
 * on failure.
 */

double gretl_gini(int t1, int t2, const(double) * x)
{
    double g;
    int err = 0;

    g = gini_coeff(x, t1, t2, null, null, &err);

    if (err) {
        g = NADBL;
    }

    return g;
}

private int lorenz_graph(const(char) * vname, double * lz, int n)
{
    FILE * fp;
    double idx;
    int downsample = 0;
    int t, err = 0;

    fp = open_plot_input_file(PlotType.regular, cast(GptFlags) 0, &err);
    if (err) {
        return err;
    }

    print_keypos_string(GpKeySpec.left | GpKeySpec.top, fp);

    fprintf(fp, "set title '%s'\n", vname);
    fprintf(fp, "plot \\\n"
            ~ "'-' using 1:2 title '%s' w lines, \\\n"
            ~ "'-' using 1:2 notitle w lines\n",
            _("Lorenz curve").ptr);

    gretl_push_c_numeric_locale();

    if (n > 4000) {
        downsample = cast(int) (n / 1000.0);
    }

    for (t=0; t<n; t++) {
        if (downsample && t % downsample) {
            continue;
        }
        idx = t + 1;
        fprintf(fp, "%g %g\n", idx / n, lz[t]);
    }
    fputs("e\n", fp);

    for (t=0; t<n; t++) {
        if (downsample && t % downsample) {
            continue;
        }
        idx = (cast(double) t + 1) / n;
        fprintf(fp, "%g %g\n", idx, idx);
    }
    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    return finalize_plot_input_file(fp);
}

/**
 * gini:
 * @varno: ID number of variable to examine.
 * @dset: dataset struct.
 * @opt: unused at present.
 * @prn: gretl printing struct.
 *
 * Graphs the Lorenz curve for variable @vnum and prints the
 * Gini coefficient.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int gini(int varno, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    double * lz = null;
    double gini_val;
    int n, fulln;
    int err = 0;

    gini_val = gini_coeff(dset.Z[varno], dset.t1, dset.t2,
                          &lz, &n, &err);

    if (err) {
        return err;
    }

    fulln = dset.t2 - dset.t1 - 1;
    pprintf(prn, "%s\n", dset.varname[varno], n);
    pprintf(prn, _("Number of observations = %d\n").ptr, n);

    if (n < fulln) {
        pputs(prn, _("Warning: there were missing values\n").ptr);
    }

    pputc(prn, '\n');

    pprintf(prn, "%s = %g\n", _("Sample Gini coefficient").ptr, gini_val);
    pprintf(prn, "%s = %g\n", _("Estimate of population value").ptr,
            gini_val * cast(double) n / (n - 1));

    err = lorenz_graph(dset.varname[varno], lz, n);

    free(lz);

    return err;
}

/* Following: apparatus for Shapiro-Wilk normality test.
   Thanks to Marcin Blazejowski for the contribution.
*/

private int min(int a, int b)
{
    return (a > b) ? b : a;
}

private int sign(int a)
{
    return (a == 0)? 0 : ((a < 0)? -1 : 1);
}

private extern (C) int compare_floats(const(void) * a, const(void) * b)
{
    const(float) * fa = cast(const(float) *) a;
    const(float) * fb = cast(const(float) *) b;

    return (*fa > *fb) - (*fa < *fb);
}

/* Algorithm AS 181.2 Appl. Statist. (1982) Vol. 31, No. 2
   Calculates the algebraic polynomial of order n-1 with
   array of coefficients cc.  Zero-order coefficient is
   cc(1) = cc[0]
*/

private double poly(const(float) * cc, int n, float x)
{
    double p, ret = cc[0]; /* preserve precision! */
    int j;

    if (n > 1) {
        p = x * cc[n-1];
        for (j=n-2; j>0; j--) {
            p = (p + cc[j]) * x;
        }
        ret += p;
    }

    return ret;
}

/* Calculate coefficients for the Shapiro-Wilk test */

private void sw_coeff(int n, int n2, float * a)
{
    immutable float[6] c1 = [ 0.0f, 0.221157f, -0.147981f, -2.07119f, 4.434685f, -2.706056f ];
    immutable float[6] c2 = [ 0.0f, 0.042981f, -0.293762f, -1.752461f, 5.682633f, -3.582633f ];
    float summ2, ssumm2;
    float a0, a1, an = cast(float) n;
    float fac, an25, rsn;
    int i, i1, nn2 = n / 2;

    if (n == 3) {
        a[0] = cast(float) sqrt(0.5);
    } else {
        an25 = an + 0.25f;
        summ2 = 0.0;
        for (i=0; i<n2; i++) {
            a[i] = cast(float) normal_cdf_inverse((i + 1 - 0.375f) / an25);
            summ2 += a[i] * a[i];
        }
        summ2 *= 2.0;
        ssumm2 = cast(float) sqrt(cast(double) summ2);
        rsn = cast(float) (1.0 / sqrt(cast(double) an));
        a0 = cast(float) (poly(c1.ptr, 6, rsn) - a[0] / ssumm2);

        /* Normalize array a */
        if (n > 5) {
            i1 = 2;
            a1 = cast(float) (-a[1] / ssumm2 + poly(c2.ptr, 6, rsn));
            fac = cast(float) sqrt((summ2 - 2.0 * (a[0] * a[0]) - 2.0 * (a[1] * a[1])) /
                                   (1.0 - 2.0 * (a0 * a0) - 2.0 * (a1 * a1)));
            a[1] = a1;
        } else {
            i1 = 1;
            fac = cast(float) sqrt((summ2 - 2.0 * (a[0] * a[0])) / (1.0 - 2.0 * (a0 * a0)));
        }
        a[0] = a0;
        for (i=i1; i<nn2; i++) {
            a[i] /= -fac;
        }
    }
}

/* ALGORITHM AS R94 APPL. STATIST. (1995) vol. 44, no. 4, 547-551.
   Calculates the Shapiro-Wilk W test and its significance level,
   Rendered into C by Marcin Blazejowski, May 2008.
*/

private int sw_w(float * x, int n, int n1, float * a,
                 double * W, double * pval)
{
    immutable float z90 = 1.2816f;
    immutable float z95 = 1.6449f;
    immutable float z99 = 2.3263f;
    immutable float zm = 1.7509f;
    immutable float zss = 0.56268f;
    immutable float bf1 = 0.8378f;
    immutable double xx90 = 0.556;
    immutable double xx95 = 0.622;
    immutable float little = 1e-19f; /* "small" is reserved on win32 */
    immutable float pi6 = 1.909859f;
    immutable float stqr = 1.047198f;

    /* polynomial coefficients */
    immutable float[2] g = [ -2.273f, 0.459f ];
    immutable float[4] c3 = [ 0.544f, -0.39978f, 0.025054f, -6.714e-4f ];
    immutable float[4] c4 = [ 1.3822f, -0.77857f, 0.062767f, -0.0020322f ];
    immutable float[4] c5 = [ -1.5861f, -0.31082f, -0.083751f, 0.0038915f ];
    immutable float[3] c6 = [ -0.4803f, -0.082676f, 0.0030302f ];
    immutable float[2] c7 = [ 0.164f, 0.533f ];
    immutable float[2] c8 = [ 0.1736f, 0.315f ];
    immutable float[2] c9 = [ 0.256f, -0.00635f ];

    float r1, zbar, ssassx, gamma, range;
    float bf, ld, m, s, sa, xi, sx, xx, y, w1;
    float asa, ssa, z90f, sax, zfm, z95f, zsd, z99f, ssx, xsx;
    float an = cast(float) n;
    int i, j, ncens = n - n1;

    /* Check for zero range */
    range = x[n1 - 1] - x[0];
    if (range < little) {
        fprintf(stderr, "sw_w: range is too small\n");
        return 1;
    }

    /* Check for correct sort order on range-scaled X */
    /* In fact we don't need do this, since we use qsort() from libc
       and we can suppose, that qsort() can sort... */
    xx = x[0] / range;
    sx = xx;
    sa = -a[0];
    j = n - 1;
    for (i = 1; i < n1; --j) {
        xi = x[i] / range;
        sx += xi;
        ++i;
        if (i != j) {
            sa += sign(i - j) * a[min(i,j) - 1];
        }
        xx = xi;
    }

    /* Calculate W statistic as squared correlation between data and
       coefficients */
    sa /= n1;
    sx /= n1;
    ssa = ssx = sax = 0.0;
    j = n - 1;
    for (i=0; i<n1; i++, j--) {
        if (i != j) {
            asa = sign(i - j) * a[min(i,j)] - sa;
        } else {
            asa = -sa;
        }
        xsx = x[i] / range - sx;
        ssa += asa * asa;
        ssx += xsx * xsx;
        sax += asa * xsx;
    }

    /* W1 equals (1-W) calculated to avoid excessive rounding error
       for W very near 1 (a potential problem in very large samples)
    */
    ssassx = cast(float) sqrt(cast(double) (ssa * ssx));
    w1 = (ssassx - sax) * (ssassx + sax) / (ssa * ssx);
    *W = 1 - w1;

    /* Calculate significance level for W */
    if (n == 3) {
        /* exact P-value */
        *pval = pi6 * (asin(sqrt(*W)) - stqr);
        return 0;
    }

    y = cast(float) log(cast(double) w1);
    xx = cast(float) log(cast(double) an);

    if (n <= 11) {
        gamma = cast(float) poly(g.ptr, 2, an);
        if (y >= gamma) {
            /* FIXME: rather use an even smaller value, or NA? */
            *pval = little;
            return 0;
        }
        y = cast(float) -log(cast(double) (gamma - y));
        m = cast(float) poly(c3.ptr, 4, an);
        s = cast(float) exp(poly(c4.ptr, 4, an));
    } else {
        /* n >= 12 */
        m = cast(float) poly(c5.ptr, 4, xx);
        s = cast(float) exp(poly(c6.ptr, 3, xx));
    }

    if (ncens > 0) {
        /* Censoring by proportion ncens/n.
           Calculate mean and sd of normal equivalent deviate of W
        */
        float delta = cast(float) ncens / an;

        ld = cast(float) -log(cast(double) delta);
        bf = 1.0f + xx * bf1;
        r1 = cast(float) pow(xx90, cast(double) xx);
        z90f = cast(float) (z90 + bf * pow(poly(c7.ptr, 2, r1), cast(double) ld));
        r1 = cast(float) pow(xx95, cast(double) xx);
        z95f = cast(float) (z95 + bf * pow(poly(c8.ptr, 2, r1), cast(double) ld));
        z99f = cast(float) (z99 + bf * pow(poly(c9.ptr, 2, xx), cast(double) ld));

        /* Regress Z90F,...,Z99F on normal deviates Z90,...,Z99 to get
           pseudo-mean and pseudo-sd of z as the slope and intercept
        */
        zfm = cast(float) ((z90f + z95f + z99f) / 3.0);
        zsd = cast(float) ((z90 * (z90f - zfm) + z95 * (z95f - zfm) + z99 * (z99f - zfm)) / zss);
        zbar = zfm - zsd * zm;
        m += zbar * s;
        s *= zsd;
    }

    *pval = normal_cdf_comp((cast(double) y - cast(double) m) / cast(double) s);

    return 0;
}

private int sw_sample_check(int n, int n1)
{
    int ncens = n - n1;
    float delta;
    float an = n;

    if (n < 3 || n1 < 3) {
        fprintf(stderr, "There is no way to run SW test for less then 3 obs\n");
        return GretlError.data;
    }

    if (ncens < 0 || (ncens > 0 && n < 20)) {
        fprintf(stderr, "sw_w: not enough uncensored obserations\n");
        return GretlError.data;
    }

    delta = cast(float) ncens / an;
    if (delta > 0.8f) {
        fprintf(stderr, "sw_w: too many censored obserations\n");
        return GretlError.data;
    }

    return 0;
}

/**
 * shapiro_wilk:
 * @x: data array.
 * @t1: starting observation.
 * @t2: ending observation.
 * @W: location to receive test statistic.
 * @pval: location to receive p-value.
 *
 * Computes the Shapiro-Wilk W statistic as a test for
 * normality of the data @x, and also the p-value for
 * the test. These are written into the pointer
 * arguments @W and @pval.
 *
 * Returns: 0 on success, non-zero on failure.
*/

int shapiro_wilk(const(double) * x, int t1, int t2, double * W, double * pval)
{
    int n1 = 0;                 /* number of uncensored obs */
    float * a = null;           /* array of coefficients */
    float * xf = null;          /* copy of x in float format */
    int i, t, n2, n = 0;
    int err = 0;

    *W = *pval = NADBL;

    for (t=t1; t<=t2; t++) {
        if (!na(x[t])) n++;
    }

    /* for now we assume all obs are uncensored */
    n1 = n;

    err = sw_sample_check(n, n1);
    if (err) {
        return err;
    }

    /* How many coeffs should be computed? */
    n2 = cast(int) fmod(cast(double) n, 2.0);
    n2 = (n2 == 0)?  n / 2 : (n - 1) / 2;

    xf = cast(float *) malloc(n * (*xf).sizeof);
    a = cast(float *) malloc(n2 * (*a).sizeof);

    if (xf == null || a == null) {
        err = GretlError.alloc;
    } else {
        i = 0;
        for (t=t1; t<=t2; t++) {
            if (!na(x[t])) {
                xf[i++] = cast(float) x[t];
            }
        }
        qsort(xf, n, (*xf).sizeof, &compare_floats);
        /* Main job: compute W stat and its p-value */
        sw_coeff(n, n2, a);
        err = sw_w(xf, n, n1, a, W, pval);
    }

    free(a);
    free(xf);

    return err;
}

private double round2(double x)
{
    x *= 100.0;
    x = (x - floor(x) < 0.5)? floor(x) : ceil(x);
    return x / 100;
}

/* Polynomial approximation to the p-value for the Lilliefors test,
   said to be good to 2 decimal places.  See Abdi and Molin,
   "Lilliefors/Van Soest's test of normality", in Salkind (ed.)
   Encyclopedia of Measurement and Statistics, Sage, 2007; also
   http://www.utd.edu/~herve/Abdi-Lillie2007-pretty.pdf
 */

private double lilliefors_pval(double L, int N)
{
    double b0 = 0.37872256037043;
    double b1 = 1.30748185078790;
    double b2 = 0.08861783849346;
    double b1pN = b1 + N;
    double Lm2 = 1.0 / (L * L);
    double A, pv;

    A = (-b1pN + sqrt(b1pN * b1pN - 4 * b2 * (b0 - Lm2))) / (2 * b2);

    pv = -0.37782822932809 + 1.67819837908004*A
        - 3.02959249450445*A*A + 2.80015798142101*pow(A, 3.)
        - 1.39874347510845*pow(A, 4.) + 0.40466213484419*pow(A, 5.)
        - 0.06353440854207*pow(A, 6.) + 0.00287462087623*pow(A, 7.)
        + 0.00069650013110*pow(A, 8.) - 0.00011872227037*pow(A, 9.)
        + 0.00000575586834*pow(A, 10.);

    if (pv < 0) {
        pv = 0;
    } else if (pv > 1) {
        pv = 1;
    } else {
        /* round to 2 digits */
        pv = round2(pv);
    }

    return pv;
}

private int lilliefors_test(const(double) * x, int t1, int t2,
                            double * L, double * pval)
{
    double * zx;   /* copy of x for sorting, z-scoring */
    int i, t, n = 0;
    int err = 0;

    *L = *pval = NADBL;

    for (t=t1; t<=t2; t++) {
        if (!na(x[t])) n++;
    }

    if (n < 5) {
        /* we need more than 4 data points */
        return GretlError.data;
    }

    zx = cast(double *) malloc(n * (*zx).sizeof);

    if (zx == null) {
        err = GretlError.alloc;
    } else {
        double dx, mx = 0.0, sx = 0.0;
        double Phi, Dp = 0.0, Dm = 0.0;
        double Dmax = 0.0;

        /* compute sample mean and standard deviation plus
           sorted z-scores */

        i = 0;
        for (t=t1; t<=t2; t++) {
            if (!na(x[t])) {
                zx[i++] = x[t];
                mx += x[t];
            }
        }
        mx /= n;

        for (t=t1; t<=t2; t++) {
            if (!na(x[t])) {
                dx = x[t] - mx;
                sx += dx * dx;
            }
        }
        sx = sqrt(sx / (n - 1));

        qsort(zx, n, (*zx).sizeof, &gretl_compare_doubles);

        for (i=0; i<n; i++) {
            zx[i] = (zx[i] - mx) / sx;
        }

        /* The following is based on the formula given in the
           "nortest" package for GNU R, function lillie.test (written
           by Juergen Gross).  It differs from the algorithm given by
           Abdi and Molin (whose p-value approximation is used above).
           The A & M account of the statistic itself seems to be
           wrong: it consistently over-rejects (i.e. around 8 percent
           of the time for a nominal 5 percent significance level).
        */

        for (i=0; i<n; i++) {
            Phi = normal_cdf(zx[i]);
            Dp = cast(double) (i+1) / n - Phi;
            Dm = Phi - cast(double) i / n;
            if (Dp > Dmax) Dmax = Dp;
            if (Dm > Dmax) Dmax = Dm;
        }

        *L = Dmax;
        *pval = lilliefors_pval(Dmax, n);
    }

    free(zx);

    return err;
}

private int skew_kurt_test(const(double) * x, int t1, int t2,
                           double * test, double * pval,
                           Gretlopt opt)
{
    double skew, xkurt;
    int n, err = 0;

    *test = *pval = NADBL;

    err = series_get_skew_kurt(t1, t2, x, &skew, &xkurt, &n);

    if (!err) {
        if (opt & Gretlopt.OPT_J) {
            /* Jarque-Bera */
            *test = (n / 6.0) * (skew * skew + xkurt * xkurt/4.0);
        } else {
            /* Doornik-Hansen */
            *test = doornik_chisq(skew, xkurt, n);
        }
    }

    if (!err && na(*test)) {
        err = GretlError.nan;
    }

    if (!na(*test)) {
        *pval = chisq_cdf_comp(2, *test);
    }

    return err;
}

private void print_normality_stat(double test, double pval,
                                  Gretlopt opt, Prn * prn)
{
    string[4] tstrs = [
        N_("Shapiro-Wilk W"),
        N_("Jarque-Bera test"),
        N_("Lilliefors test"),
        N_("Doornik-Hansen test")
    ];
    int i = (opt & Gretlopt.OPT_W)? 0 : (opt & Gretlopt.OPT_J)? 1 :
        (opt & Gretlopt.OPT_L)? 2 : 3;

    if (na(test) || na(pval)) {
        return;
    }

    if (opt & Gretlopt.OPT_L) {
        pprintf(prn, " %s = %g, %s ~= %g\n\n",
                _(tstrs[i]).ptr, test, _("with p-value").ptr, pval);
    } else {
        pprintf(prn, " %s = %g, %s %g\n\n",
                _(tstrs[i]).ptr, test, _("with p-value").ptr, pval);
    }
}

/**
 * gretl_normality_test:
 * @varno: ID number of the variable to process.
 * @dset: dataset struct.
 * @opt: OPT_A: all tests; OPT_D: Doornik-Hansen; OPT_W: Shapiro-Wilk;
 * OPT_J: Jarque-Bera; OPT_L: Lilliefors; default is Doornik-Hansen.
 * Also use OPT_Q for quiet.
 * @prn: gretl printing struct.
 *
 * Performs, and prints the results of, the specified test(s) randomness
 * for the variable specified by @v.
 *
 * Returns: 0 on successful completion, non-zero on error.
 */

int gretl_normality_test(int varno, const(DatasetRaw) * dset,
                         Gretlopt opt, Prn * prn)
{
    Gretlopt alltests = Gretlopt.OPT_D | Gretlopt.OPT_W | Gretlopt.OPT_J | Gretlopt.OPT_L;
    double test = NADBL;
    double pval = NADBL;
    double trec = NADBL;
    double pvrec = NADBL;
    int err;

    if (varno < 0 || varno >= dset.v) {
        return GretlError.data;
    }

    err = incompatible_options(opt, Gretlopt.OPT_A | alltests);
    if (err) {
        return err;
    }

    if (opt & Gretlopt.OPT_A) {
        /* show all tests */
        opt |= alltests;
    }

    if (!(opt & alltests)) {
        /* no method selected: use Doornik-Hansen */
        opt |= Gretlopt.OPT_D;
    }

    if (!(opt & Gretlopt.OPT_Q)) {
        pprintf(prn, _("Test for normality of %s:").ptr, dset.varname[varno]);
        if (opt & Gretlopt.OPT_A) {
            pputs(prn, "\n\n");
        } else {
            pputc(prn, '\n');
        }
    }

    if (opt & Gretlopt.OPT_D) {
        err = skew_kurt_test(dset.Z[varno], dset.t1, dset.t2,
                             &test, &pval, Gretlopt.OPT_D);
        if (!err && !(opt & Gretlopt.OPT_Q)) {
            print_normality_stat(test, pval, Gretlopt.OPT_D, prn);
        }
        if (!err) {
            /* record Doornik-Hansen result by default */
            trec = test;
            pvrec = pval;
        }
    }

    if (opt & Gretlopt.OPT_W) {
        err = shapiro_wilk(dset.Z[varno], dset.t1, dset.t2,
                           &test, &pval);
        if (!err && !(opt & Gretlopt.OPT_Q)) {
            print_normality_stat(test, pval, Gretlopt.OPT_W, prn);
        }
    }

    if (opt & Gretlopt.OPT_L) {
        err = lilliefors_test(dset.Z[varno], dset.t1, dset.t2,
                              &test, &pval);
        if (!err && !(opt & Gretlopt.OPT_Q)) {
            print_normality_stat(test, pval, Gretlopt.OPT_L, prn);
        }
    }

    if (opt & Gretlopt.OPT_J) {
        err = skew_kurt_test(dset.Z[varno], dset.t1, dset.t2,
                             &test, &pval, Gretlopt.OPT_J);
        if (!err && !(opt & Gretlopt.OPT_Q)) {
            print_normality_stat(test, pval, Gretlopt.OPT_J, prn);
        }
    }

    if (na(trec) && !na(test)) {
        trec = test;
    }

    if (na(pvrec) && !na(pval)) {
        pvrec = pval;
    }

    if (!na(trec) && !na(pvrec)) {
        record_test_result(trec, pvrec);
    }

    return err;
}

GretlMatrixRaw * gretl_normtest_matrix(const(double) * y,
                                    int t1, int t2,
                                    Gretlopt opt,
                                    int * err)
{
    GretlMatrixRaw * ret = null;
    double test = NADBL;
    double pval = NADBL;
    int do_all = 0;

    if (opt & Gretlopt.OPT_J) {
        /* Jarque-Bera */
        *err = skew_kurt_test(y, t1, t2, &test, &pval, opt);
    } else if (opt & Gretlopt.OPT_W) {
        /* Shapiro-Wilk */
        *err = shapiro_wilk(y, t1, t2, &test, &pval);
    } else if (opt & Gretlopt.OPT_L) {
        /* Lilliefors */
        *err = lilliefors_test(y, t1, t2, &test, &pval);
    } else if (opt & Gretlopt.OPT_A) {
        /* all tests */
        do_all = 1;
    } else {
        /* default: Doornik-Hansen */
        *err = skew_kurt_test(y, t1, t2, &test, &pval, opt);
    }

    if (do_all) {
        ret = gretl_matrix_alloc(4, 2);
        if (ret == null) {
            *err = GretlError.alloc;
        } else {
            char ** Sc;
            char ** Sr;

            Sc = strings_array_new(2);
            Sr = strings_array_new(4);
            if (Sc != null) {
                Sc[0] = gretl_strdup("test");
                Sc[1] = gretl_strdup("p-value");
                gretl_matrix_set_colnames(ret, Sc);
            }
            if (Sr != null) {
                Sr[0] = gretl_strdup("Doornik-Hansen");
                Sr[1] = gretl_strdup("Shapiro-Wilk");
                Sr[2] = gretl_strdup("Jarque-Bera");
                Sr[3] = gretl_strdup("Lilliefors");
                gretl_matrix_set_rownames(ret, Sr);
            }

            /* Doornik-Hansen */
            skew_kurt_test(y, t1, t2, &test, &pval, Gretlopt.OPT_NONE);
            gretl_matrix_set(ret, 0, 0, test);
            gretl_matrix_set(ret, 0, 1, pval);
            /* Shapiro-Wilk */
            shapiro_wilk(y, t1, t2, &test, &pval);
            gretl_matrix_set(ret, 1, 0, test);
            gretl_matrix_set(ret, 1, 1, pval);
            /* Jarque-Bera */
            skew_kurt_test(y, t1, t2, &test, &pval, Gretlopt.OPT_J);
            gretl_matrix_set(ret, 2, 0, test);
            gretl_matrix_set(ret, 2, 1, pval);
            /* Lilliefors */
            lilliefors_test(y, t1, t2, &test, &pval);
            gretl_matrix_set(ret, 3, 0, test);
            gretl_matrix_set(ret, 3, 1, pval);
        }
        return ret;
    }

    if (!*err) {
        ret = gretl_matrix_alloc(1, 2);
        if (ret == null) {
            *err = GretlError.alloc;
        } else {
            ret.val[0] = test;
            ret.val[1] = pval;
        }
    }

    return ret;
}
