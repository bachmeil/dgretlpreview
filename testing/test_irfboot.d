/* Smoke tests for irfboot.d on a simulated bivariate VAR(1): bootstrap confidence bands for impulse responses obtained
   through gretl_VAR_get_impulse_response with a non-zero significance level.

   Not exercised: VECM bootstraps (they need the johansen plugin), studentized intervals and the exact values (the
   bootstrap is random, so only the ordering of the bands and their coverage of the point estimate are checked).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs, sqrt;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import gretl.dataset : create_new_dataset;
import gretl.estimate : lsq;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_matrixh : GretlMatrixRaw;
import gretl.gretl_model : clear_model;
import gretl.random : gretl_rand_set_seed;
import gretl.libgretlh;
import gretl.irfboot;
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
    gretl_rand_set_seed(20260101);

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

    GretlVar * v = gretl_VAR(1, null, ylist.ptr, d, Gretlopt.OPT_Q, null, &err);

    check(err == 0 && v !is null, "gretl_VAR(1) estimates");

    GretlMatrixRaw * pt = gretl_VAR_get_impulse_response(v, 0, 0, 5, 0.0, d, &err);

    check(err == 0 && pt !is null && pt.rows == 5, "point impulse response has 5 rows");

    GretlMatrixRaw * bt = gretl_VAR_get_impulse_response(v, 0, 0, 5, 0.1, d, &err);

    check(err == 0 && bt !is null && bt.rows == 5 && bt.cols == 3, "bootstrapped response has 5 rows and 3 columns (point, lower, upper)");
    if (bt !is null && bt.cols == 3) {
        bool ordered = true, point_ok = true;
        foreach (i; 0 .. 5) {
            double p = bt.val[i], lo = bt.val[i + 5], hi = bt.val[i + 10];
            ordered = ordered && lo <= hi;
            point_ok = point_ok && close(p, pt.val[i]);
        }
        check(ordered, "lower band <= upper band at every horizon");
        check(point_ok, "first column is the point estimate");
        check(bt.val[5] <= bt.val[0] + 1e-9 && bt.val[0] <= bt.val[10] + 1e-9, "the band covers the point estimate at horizon 0");
    }

    gretl_VAR_free(v);
    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
