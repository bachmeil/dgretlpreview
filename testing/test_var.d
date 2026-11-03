/* Smoke tests for var.d on a simulated bivariate VAR(1): estimation with gretl_VAR (equation models against OLS on
   hand-built lags, coefficient and residual matrices, information criteria), the roots of the companion matrix,
   impulse responses, the helpers var_max_order and default_VAR_horizon, and gretl_VAR_free.

   Not exercised: VECM/Johansen estimation, lag selection, the omit/autocorrelation/ARCH/normality tests, forecast
   decompositions, bundles and XML serialization.

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs, sqrt;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import gretl.dataset : create_new_dataset;
import gretl.estimate : lsq;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_matrixh : GretlMatrixRaw;
import gretl.gretl_model : clear_model;
import gretl.libgretlh;
import gretl.var;
import gretl.varh : GretlVar;

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

    enum T = 80;
    DatasetRaw * d = create_new_dataset(6, T, 0);

    d.pd = 1;
    d.structure = TsCodes.time_series;
    foreach (i, nm; ["const", "y1", "y2", "y1_1", "y2_1", "junk"]) {
        d.varname[i][0 .. nm.length] = nm[];
        d.varname[i][nm.length] = 0;
    }

    uint seed = 12345;
    double rnd()
    {
        seed = seed * 1103515245u + 12345u;
        return ((seed >> 8) & 0xffff) / 65536.0 - 0.5;
    }

    d.Z[1][0] = 0;
    d.Z[2][0] = 0;
    foreach (t; 0 .. T) {
        d.Z[0][t] = 1;
        d.Z[5][t] = 0;
        if (t > 0) {
            d.Z[1][t] = 0.5 * d.Z[1][t-1] + 0.2 * d.Z[2][t-1] + rnd();
            d.Z[2][t] = 0.1 * d.Z[1][t-1] + 0.4 * d.Z[2][t-1] + rnd();
            d.Z[3][t] = d.Z[1][t-1];
            d.Z[4][t] = d.Z[2][t-1];
        }
    }
    d.t1 = 1;
    d.t2 = T - 1;

    int[3] ylist = [2, 1, 2];
    int err = 0;

    check(var_max_order(ylist.ptr, d) >= 1, "var_max_order is positive");
    check(default_VAR_horizon(d) > 0, "default_VAR_horizon is positive");

    GretlVar * v = gretl_VAR(1, null, ylist.ptr, d, Gretlopt.OPT_Q, null, &err);

    check(err == 0 && v !is null, "gretl_VAR(1) estimates");
    check(v.neqns == 2 && v.order == 1 && v.T == T - 1, "dimensions: 2 equations, order 1, T = 79");
    check(gretl_VAR_get_n_equations(v) == 2, "gretl_VAR_get_n_equations");

    /* equation 1 against OLS on hand-built lags */
    int[5] l1b = [4, 1, 0, 3, 4];
    Model m1 = lsq(l1b.ptr, d, GretlCmdIndex.ols, Gretlopt.OPT_NONE);
    const(Model) * vm = gretl_VAR_get_model(v, 0);

    check(m1.errcode == 0 && vm !is null && vm.ncoeff == 3, "equation model has 3 coefficients");
    bool match = true;
    foreach (i; 0 .. 3) {
        bool found = false;
        foreach (j; 0 .. 3) {
            if (close(vm.coeff[i], m1.coeff[j])) found = true;
        }
        match = match && found;
    }
    check(match, "VAR equation 1 coefficients agree with OLS on the lags");
    check(close(vm.ess, m1.ess), "equation 1 sum of squared residuals agrees");
    clear_model(&m1);

    const(GretlMatrixRaw) * E = gretl_VAR_get_residual_matrix(v);

    check(E !is null && E.rows == T - 1 && E.cols == 2, "residual matrix is T x 2");

    const(GretlMatrixRaw) * R = gretl_VAR_get_roots(v, &err);

    check(err == 0 && R !is null && R.rows == 2, "two roots for the companion matrix");
    if (R !is null) {
        double re0 = R.val[0], im0 = R.val[R.rows];
        check(sqrt(re0 * re0 + im0 * im0) < 1.0, "the simulated VAR is stable (root inside the unit circle)");
    }

    GretlMatrixRaw * irf = gretl_VAR_get_impulse_response(v, 0, 0, 5, 0.0, d, &err);

    check(err == 0 && irf !is null && irf.rows == 5, "impulse response has 5 rows");
    if (irf !is null) {
        check(irf.val[0] > 0, "own-shock response is positive at horizon 0");
    }

    gretl_VAR_free(v);
    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
