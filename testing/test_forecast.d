/* Smoke tests for forecast.d on an OLS model of a small dataset: fitted values and residuals (get_fit_resid),
   static out-of-sample forecasts with standard errors (get_forecast), matrix_forecast, the recursive k-step
   forecast and the argument handling of get_forecast.

   Not exercised: dynamic forecasts with lagged dependent variables, ARMA/GARCH/NLS/limited-dependent and VAR/system
   forecasts, and the "fcast" command (do_forecast), which needs the genr machinery.

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs, sqrt;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import gretl.dataset : create_new_dataset, destroy_dataset;
import gretl.estimate : lsq;
import gretl.forecast;
import gretl.forecasth : FitResid;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_matrix : gretl_matrix_alloc, gretl_matrix_free;
import gretl.gretl_matrixh : GretlMatrixRaw;
import gretl.gretl_model : clear_model;
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

bool close(double a, double b, double tol = 1e-8)
{
    return fabs(a - b) <= tol * (1.0 + fabs(b));
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    enum n = 14;
    immutable double[n] noise = [0.3, -0.2, 0.5, -0.6, 0.1, 0.4, -0.3, -0.2, -0.5, 0.2, 0.0, 0.0, 0.0, 0.0];
    DatasetRaw * d = create_new_dataset(3, n, 0);

    d.pd = 1;
    d.structure = TsCodes.time_series;
    foreach (i, nm; ["const", "y", "x"]) {
        d.varname[i][0 .. nm.length] = nm[];
        d.varname[i][nm.length] = 0;
    }
    foreach (t; 0 .. n) {
        d.Z[0][t] = 1;
        d.Z[2][t] = t + 1.0;
        d.Z[1][t] = 1.0 + 2.0 * (t + 1.0) + noise[t];
    }
    d.t1 = 0;
    d.t2 = 9;
    int[4] list = [3, 1, 0, 2];
    Model m = lsq(list.ptr, d, GretlCmdIndex.ols, Gretlopt.OPT_NONE);

    check(m.errcode == 0 && m.ncoeff == 2, "OLS model estimated over the first ten observations");

    int err = 0;
    FitResid * fr = get_fit_resid(&m, d, &err);

    check(err == 0 && fr !is null && fr.nobs >= 10, "get_fit_resid");
    check(close(fr.fitted[3], m.coeff[0] + m.coeff[1] * 4.0), "fitted value = b0 + b1 x");
    check(close(fr.resid[3], fr.actual[3] - fr.fitted[3]), "residual = actual - fitted");
    check(close(fr.sigma, m.sigma), "standard error of regression carried over");
    free_fit_resid(fr);

    /* static forecast over the four observations after the sample */
    fr = get_forecast(&m, 10, 13, 0, d, Gretlopt.OPT_NONE, &err);
    check(err == 0 && fr !is null && fr.t1 == 10 && fr.t2 == 13, "get_forecast over 10..13");
    double[4] expect;

    foreach (i; 0 .. 4) expect[i] = m.coeff[0] + m.coeff[1] * (11.0 + i);
    bool allok = true;

    foreach (i; 0 .. 4) allok = allok && close(fr.fitted[10 + i], expect[i]);
    check(allok, "forecasts are the fitted line at the new x values");
    check(fr.sderr !is null && close(fr.sderr[10], m.sigma * sqrt(1 + 1.0 / 10 + (11 - 5.5) * (11 - 5.5) / 82.5), 1e-6) && fr.sderr[13] > fr.sderr[10], "forecast standard errors follow the textbook formula");
    free_fit_resid(fr);

    /* matrix forecast */
    GretlMatrixRaw * X = gretl_matrix_alloc(3, 2);

    X.val[0] = 1; X.val[1] = 1; X.val[2] = 1;
    X.val[3] = 2; X.val[4] = 4; X.val[5] = 6;
    GretlMatrixRaw * F = matrix_forecast(&m, X, &err);

    check(err == 0 && F !is null && F.rows == 3 && close(F.val[1], m.coeff[0] + m.coeff[1] * 4), "matrix_forecast");
    gretl_matrix_free(F);
    gretl_matrix_free(X);

    /* recursive one-step-ahead forecasts */
    fr = recursive_OLS_k_step_fcast(&m, d, 6, 13, 1, 0, &err);
    check(err == 0 && fr !is null && fr.k == 1, "recursive_OLS_k_step_fcast runs");
    if (fr !is null) free_fit_resid(fr);

    /* a range that starts before the estimation sample ends still works; a reversed range does not */
    fr = get_forecast(&m, 13, 10, 0, d, Gretlopt.OPT_NONE, &err);
    check(fr is null && err != 0, "a reversed forecast range is an error");
    err = 0;

    clear_model(&m);
    destroy_dataset(d);
    printf("\n%d checks, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
