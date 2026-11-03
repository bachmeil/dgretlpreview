/* Smoke tests for vartest.d on a simulated bivariate VAR(1): the lag-omission F statistics and portmanteau statistic
   computed during estimation, the multivariate autocorrelation, ARCH and normality tests, and the Wald omission
   test.

   Not exercised: VECM-specific tests, lag-selection tables, the likelihood-ratio omit test with restricted VECMs and
   the exact test statistics (only their presence and ranges are checked).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs, sqrt;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import gretl.dataset : create_new_dataset;
import gretl.estimate : lsq;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_matrixh : GretlMatrixRaw;
import gretl.gretl_model : clear_model;
import gretl.libgretlh;
import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_get_buffer, gretl_print_new;
import gretl.gretl_prnh : PrnType;
import core.stdc.string : strlen;
import gretl.var;
import gretl.vartest;
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
    Prn * prn = gretl_print_new(PrnType.buffer, &err);

    GretlVar * v = gretl_VAR(1, null, ylist.ptr, d, Gretlopt.OPT_Q, null, &err);

    check(err == 0 && v !is null, "gretl_VAR(1) estimates");
    check(v.Fvals !is null, "lag-omission F statistics were computed");
    if (v.Fvals !is null) {
        check(v.Fvals[0] >= 0 && v.Fvals[0] == v.Fvals[0], "first omission F statistic is a number");
    }
    check(v.LB > 0 && v.LBs > 0, "portmanteau statistic and its order are set");

    err = gretl_VAR_autocorrelation_test(v, 2, d, Gretlopt.OPT_NONE, prn);
    const(char) * buf = gretl_print_get_buffer(prn);
    check(err == 0 && buf !is null && strlen(buf) > 0, "multivariate autocorrelation test (order 2) prints");
    gretl_print_destroy(prn);

    prn = gretl_print_new(PrnType.buffer, &err);
    err = gretl_VAR_arch_test(v, 1, d, Gretlopt.OPT_NONE, prn);
    buf = gretl_print_get_buffer(prn);
    check(err == 0 && buf !is null && strlen(buf) > 0, "multivariate ARCH test (order 1) prints");
    gretl_print_destroy(prn);

    prn = gretl_print_new(PrnType.buffer, &err);
    err = gretl_VAR_normality_test(v, Gretlopt.OPT_NONE, prn);
    buf = gretl_print_get_buffer(prn);
    check(err == 0 && buf !is null && strlen(buf) > 0, "multivariate normality test prints");
    gretl_print_destroy(prn);

    int[2] omit = [1, 2];
    prn = gretl_print_new(PrnType.buffer, &err);
    err = gretl_VAR_wald_omit_test(v, omit.ptr, d, Gretlopt.OPT_NONE, prn);
    buf = gretl_print_get_buffer(prn);
    check(err == 0 && buf !is null && strlen(buf) > 0, "Wald test for omitting variable 2's lags prints");
    gretl_print_destroy(prn);

    gretl_VAR_free(v);
    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
