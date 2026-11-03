/* Smoke tests for describe.d: the basic statistics (mean, variance, quantiles, moments, covariance and
   correlation), the Gini coefficient, frequency distributions, summary statistics (computed and printed),
   correlation matrices, ACF/Ljung-Box/cross-correlation, the periodogram, the normality tests, the
   difference-of-means and variance tests, distances and the other small helpers. Reference values were
   worked out by hand or in a few lines of Python.

   Not exercised: functions that print through or draw with unported code (frequency tables, cross-tabulations,
   graphs) and the periodogram (it needs gretl_matrix_fft, not ported yet).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs, sqrt;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;
import core.stdc.stdlib : free;
import core.stdc.string : strcmp, strstr;

import gretl.dataioh : DataClearCode;
import gretl.dataset : clear_datainfo, create_new_dataset, destroy_dataset;
import gretl.describe;
import gretl.describeh : FreqDist, Summary;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_matrixh : GretlMatrixRaw;
import gretl.gretl_matrix : gretl_matrix_alloc, gretl_matrix_free;
import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_get_buffer, gretl_print_new;
import gretl.gretl_prnh : PrnType;
import gretl.libgretlh;
import gretl.missingh : NADBL, na;

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

bool close(double a, double b, double tol = 1e-9)
{
    return fabs(a - b) <= tol * (1.0 + fabs(b));
}

void setname(DatasetRaw * d, int i, string name)
{
    d.varname[i][0 .. name.length] = name[];
    d.varname[i][name.length] = 0;
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    immutable double[8] x = [2, 4, 4, 4, 5, 5, 7, 9];
    immutable double[8] y = [1, 3, 2, 5, 4, 6, 8, 7];
    double[9] xm = [2, 4, double.nan, 4, 5, 5, 7, 9, 4];
    int err;

    /* basic statistics */
    check(gretl_sum(0, 7, x.ptr) == 40 && gretl_mean(0, 7, x.ptr) == 5, "sum and mean");
    check(gretl_min(0, 7, x.ptr) == 2 && gretl_max(0, 7, x.ptr) == 9, "min and max");
    double lo, hi;
    check(gretl_minmax(0, 8, xm.ptr, &lo, &hi) == 8 && lo == 2 && hi == 9, "minmax skips missing values");
    check(close(gretl_variance(0, 7, x.ptr), 32.0 / 7) && close(gretl_stddev(0, 7, x.ptr), 2.138089935299395),
          "variance and standard deviation");
    check(close(real_gretl_variance(0, 7, x.ptr, 1), 4.0), "variance with the asymptotic divisor");
    check(gretl_sst(0, 7, x.ptr) == 32, "sum of squared deviations");
    check(gretl_median(0, 7, x.ptr) == 4.5, "median");
    err = 0;
    check(gretl_quantile(0, 7, x.ptr, 0.75, Gretlopt.OPT_NONE, &err) == 6.5 && err == 0, "75th percentile");
    err = 0;
    check(na(gretl_quantile(0, 7, x.ptr, 0.99, Gretlopt.OPT_NONE, &err)) && err == 0, "too few observations for a quantile: NA, no error");
    err = 0;
    gretl_quantile(0, 7, x.ptr, 1.5, Gretlopt.OPT_NONE, &err);
    check(err == GretlError.data, "quantile of p > 1 is an error");
    double[8] a = [9, 2, 7, 4, 5, 4, 5, 4];
    check(gretl_array_quantile(a.ptr, 8, 0.5) == 4.5, "array quantile (reorders the array)");
    check(close(gretl_skewness(0, 7, x.ptr), 0.65625) && close(gretl_kurtosis(0, 7, x.ptr), -0.21875), "skewness and excess kurtosis");
    double xbar, sd, skew, kurt;
    check(gretl_moments(0, 7, x.ptr, null, &xbar, &sd, &skew, &kurt, 1) == 0 && xbar == 5 && close(sd, 2.138089935299395)
          && close(skew, 0.65625) && close(kurt, -0.21875), "gretl_moments");
    check(close(gretl_long_run_variance(0, 7, x.ptr, 0, NADBL), 4.0), "long-run variance with no lags is the ML variance");

    /* covariance and correlation */
    int miss;
    check(close(gretl_covar(0, 7, x.ptr, y.ptr, &miss), 4.428571428571429) && miss == 0, "covariance");
    check(close(gretl_corr(0, 7, x.ptr, y.ptr, &miss), 0.8455943246644706), "correlation");
    check(close(gretl_corr_rsq(0, 7, x.ptr, y.ptr), 0.8455943246644706 ^^ 2), "squared correlation");
    check(close(gretl_corr(0, 7, x.ptr, x.ptr, null), 1.0), "corr(x,x) = 1");
    immutable double[4] cst = [3, 3, 3, 3];
    check(na(gretl_corr(0, 3, cst.ptr, cst.ptr, null)), "correlation with a constant is NA");

    /* restricted statistics and comparisons */
    check(eval_ytest(2, GretlOp.gte, 2) == 1 && eval_ytest(1, GretlOp.gt, 2) == 0 && eval_ytest(NADBL, GretlOp.eq, NADBL) == 1,
          "eval_ytest");
    check(close(gretl_restricted_mean(0, 7, x.ptr, y.ptr, GretlOp.gt, 4), 6.25), "restricted mean (x where y > 4)");

    /* Gini and small helpers */
    immutable double[4] g4 = [1, 2, 3, 4];
    check(close(gretl_gini(0, 3, g4.ptr), 0.25), "Gini coefficient");
    check(satterthwaite_df(4.0, 10, 9.0, 12) == 19, "Satterthwaite degrees of freedom");

    /* ACF family */
    double[8] ys = y;
    err = 0;
    GretlMatrixRaw * A = acf_matrix(ys.ptr, 2, null, 8, &err);
    check(err == 0 && A !is null && close(A.val[0], 0.49404761904761907) && close(A.val[1], 0.2857142857142857),
          "acf_matrix: ACF at lags 1 and 2");
    check(A !is null && close(A.val[2], 0.49404761904761907), "acf_matrix: PACF(1) = ACF(1)");
    gretl_matrix_free(A);
    err = 0;
    check(close(ljung_box(2, 0, 7, y.ptr, &err), 3.8779559442824754), "Ljung-Box statistic");
    err = 0;
    GretlMatrixRaw * X = xcf_vec(x.ptr, y.ptr, 1, null, 8, &err);
    check(err == 0 && X !is null && X.rows == 3 && close(X.val[1], 0.8455943246644705), "cross-correlation at lag 0");
    gretl_matrix_free(X);
    check(auto_acf_order(100) == 20, "automatic ACF order");

    /* distances */
    GretlMatrixRaw * Z = gretl_matrix_alloc(2, 2);
    Z.val[0] = 0; Z.val[1] = 3; Z.val[2] = 0; Z.val[3] = 4;       /* rows: (0,0) and (3,4) */
    err = 0;
    GretlMatrixRaw * D = distance(Z, null, "euclidean", &err);
    check(err == 0 && D !is null && D.rows == 1 && close(D.val[0], 5.0), "Euclidean distance");
    gretl_matrix_free(D);
    D = distance(Z, null, "manhattan", &err);
    check(D !is null && close(D.val[0], 7.0), "Manhattan distance");
    gretl_matrix_free(D);
    D = distance(Z, null, "chebyshev", &err);
    check(D !is null && close(D.val[0], 4.0), "Chebyshev distance");
    gretl_matrix_free(D);
    err = 0;
    D = distance(Z, null, "bogus", &err);
    check(D is null && err == GretlError.invarg, "invalid distance metric");
    gretl_matrix_free(Z);

    /* normality tests */
    double[40] nz;
    double[40] sk;
    foreach (i; 0 .. 40) {
        import gretl.gretl_matrix : gretl_matrix_alloc;
        import gretl.pvalues : normal_cdf_inverse;

        nz[i] = normal_cdf_inverse((i + 0.5) / 40.0);
        sk[i] = (i + 1.0) ^^ 4;
    }
    double W, pv;
    check(shapiro_wilk(nz.ptr, 0, 39, &W, &pv) == 0 && W > 0.97 && W <= 1.0 && pv > 0.5, "Shapiro-Wilk: normal quantiles look normal");
    check(shapiro_wilk(sk.ptr, 0, 39, &W, &pv) == 0 && W < 0.9 && pv < 0.01, "Shapiro-Wilk: strongly skewed data rejected");
    err = 0;
    GretlMatrixRaw * N = gretl_normtest_matrix(sk.ptr, 0, 39, Gretlopt.OPT_A, &err);
    check(err == 0 && N !is null && N.rows == 4 && N.cols == 2, "gretl_normtest_matrix: all tests");
    if (N !is null && N.rows == 4) {
        check(N.val[4 + 0] < 0.01 && N.val[4 + 1] < 0.01 && N.val[4 + 2] < 0.01 && N.val[4 + 3] < 0.05,
              "all four normality tests reject the skewed data");
    }
    gretl_matrix_free(N);
    check(shapiro_wilk(x.ptr, 0, 1, &W, &pv) != 0, "Shapiro-Wilk needs at least 3 observations");

    /* a dataset for the rest */
    DatasetRaw * d = create_new_dataset(4, 8, 0);

    d.pd = 1;
    d.structure = TsCodes.cross_section;
    d.t1 = 0;
    d.t2 = 7;
    setname(d, 1, "x");
    setname(d, 2, "y");
    setname(d, 3, "g");
    foreach (t; 0 .. 8) {
        d.Z[1][t] = x[t];
        d.Z[2][t] = y[t];
        d.Z[3][t] = t < 4 ? 0 : 1;
    }

    /* sorted series, frequency distribution */
    int nn = 1;

    err = 0;
    double * srt = gretl_sorted_series(1, d, Gretlopt.OPT_NONE, &nn, &err);
    check(err == 0 && srt !is null && nn == 8 && srt[0] == 2 && srt[7] == 9 && srt[3] == 4, "gretl_sorted_series");
    free(srt);
    err = 0;
    FreqDist * fd = get_discrete_freq(1, d, Gretlopt.OPT_NONE, &err);
    check(err == 0 && fd !is null && fd.numbins == 5 && fd.f[0] == 1 && fd.f[1] == 3 && fd.f[2] == 2, "get_discrete_freq: counts 1,3,2,... for 2,4,5,7,9");
    free_freq(fd);
    err = 0;
    GretlMatrixRaw * fm = freqdist_matrix(y.ptr, 0, 7, &err);
    check(err == 0 && fm !is null && fm.cols == 2, "freqdist_matrix has midpoints and counts");
    gretl_matrix_free(fm);
    immutable double[6] u = [1, 2, 3, 4, 5, 6];
    immutable double[6] v = [10, 10, 20, 20, 30, 30];
    err = 0;
    check(correspondence(u.ptr, v.ptr, 6, &err) == -1 && err == 0, "correspondence: n-to-1 relationship");

    /* summary statistics */
    int[2] lst = [1, 1];
    Prn * prn = gretl_print_new(PrnType.buffer, null);

    err = 0;
    Summary * sm = get_summary(lst.ptr, d, Gretlopt.OPT_NONE, prn, &err);
    check(err == 0 && sm !is null && sm.n == 8 && sm.mean[0] == 5 && sm.low[0] == 2 && sm.high[0] == 9 && sm.median[0] == 4.5,
          "get_summary: n, mean, min, max, median");
    check(sm !is null && close(sm.sd[0], 2.138089935299395) && close(sm.skew[0], 0.65625) && summary_has_missing_values(sm) == 0,
          "get_summary: sd, skewness, no missing values");
    print_summary(sm, d, prn);
    const(char) * out_ = gretl_print_get_buffer(prn);
    check(strstr(out_, "Summary statistics") !is null && strstr(out_, "Mean") !is null && strstr(out_, "Median") !is null,
          "print_summary writes the table");
    free_summary(sm);

    /* restricted (by dummy) summary */
    int[2] lst2 = [1, 1];
    err = 0;
    sm = get_summary_restricted(lst2.ptr, d, d.Z[3], Gretlopt.OPT_NONE, prn, &err);
    check(err == 0 && sm !is null && sm.n == 4 && sm.mean[0] == 6.5, "get_summary_restricted: mean over the g==1 observations");
    free_summary(sm);

    /* correlation matrix */
    int[3] clist = [2, 1, 2];
    err = 0;
    VMatrix * vm = corrlist(GretlCmdIndex.corr, clist.ptr, d, Gretlopt.OPT_NONE, &err);
    check(err == 0 && vm !is null && vm.dim == 2 && close(vm.vec[1], 0.8455943246644706), "corrlist: correlation coefficient");
    gretl_print_destroy(prn);
    prn = gretl_print_new(PrnType.buffer, null);
    if (vm !is null) {
        print_corrmat(vm, d, prn);
        check(strstr(gretl_print_get_buffer(prn), "corr(x, y)") !is null, "print_corrmat for two variables");
    }
    free_vmatrix(vm);

    /* tests: difference of means and variances */
    gretl_print_destroy(prn);
    prn = gretl_print_new(PrnType.buffer, null);
    int[3] pair = [2, 1, 2];

    err = means_test(pair.ptr, d, Gretlopt.OPT_NONE, prn);
    out_ = gretl_print_get_buffer(prn);
    check(err == 0 && strstr(out_, "Equality of means test") !is null && strstr(out_, "Difference between sample means = 5 - 4.5 = 0.5") !is null,
          "means_test");
    gretl_print_destroy(prn);
    prn = gretl_print_new(PrnType.buffer, null);
    err = vars_test(pair.ptr, d, Gretlopt.OPT_NONE, prn);
    out_ = gretl_print_get_buffer(prn);
    check(err == 0 && strstr(out_, "Equality of variances test") !is null && strstr(out_, "Ratio of sample variances = 0.761905") !is null,
          "vars_test");
    gretl_print_destroy(prn);
    prn = gretl_print_new(PrnType.buffer, null);
    err = gretl_normality_test(1, d, Gretlopt.OPT_NONE, prn);
    check(err == 0 && strstr(gretl_print_get_buffer(prn), "Test for normality of x") !is null, "gretl_normality_test (Doornik-Hansen)");

    /* last-result store */
    GretlType rtype;

    err = 0;
    void * lr = get_last_result_data(&rtype, &err);
    check(lr !is null || err != 0, "get_last_result_data does not crash");
    last_result_cleanup();

    gretl_print_destroy(prn);
    destroy_dataset(d);
    printf("\n%d checks, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
