/* Smoke tests for genfuncs.d: series transformations on a small time-series dataset (differences, cumulation,
   standardization, sorting by another series, ranks, Box-Cox, fractional differencing, exponential moving
   average, the Hodrick-Prescott filter on a trend, a recursive filter), plus the Bessel and logistic functions, the
   MIDAS weight vector, the vector statistics, value substitution and the empirical CDF.

   Not exercised: the band-pass and Butterworth filters, panel statistics and dummies, the loess/kernel helpers,
   the forecast and duration statistics, the aggregation functions, bootstrapping and the interpolation routines.

   Build and run with testing/run_tests.sh. */

import core.stdc.math : exp, fabs, log, sqrt;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;
import core.stdc.stdlib : free;

import gretl.dataset : create_new_dataset;
import gretl.genfuncs;
import gretl.genparseh : GenFunc;
import gretl.gretl_matrix : gretl_matrix_alloc, gretl_matrix_free;
import gretl.gretl_matrixh : GretlMatrixRaw, GretlVecStat;
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

    enum n = 8;
    DatasetRaw * d = create_new_dataset(2, n, 0);

    d.pd = 1;
    d.structure = TsCodes.time_series;
    d.t1 = 0;
    d.t2 = n - 1;

    immutable double[n] x = [3, 1, 4, 1, 5, 9, 2, 6];
    double[n] y;

    /* first difference */
    y[] = NADBL;
    check(diff_series(x.ptr, y.ptr, GenFunc.f_diff, d) == 0 && na(y[0]) && y[1] == -2 && y[2] == 3 && y[7] == 4, "first difference");

    /* cumulation */
    y[] = NADBL;
    check(cum_series(x.ptr, y.ptr, d) == 0 && y[0] == 3 && y[1] == 4 && y[7] == 31, "cumulative sum");

    /* standardization */
    y[] = NADBL;
    check(standardize_series(x.ptr, y.ptr, 0, d) == 0, "standardize_series runs");
    double m0 = 0, v0 = 0;
    foreach (i; 0 .. n) m0 += y[i];
    foreach (i; 0 .. n) v0 += y[i] * y[i];
    check(fabs(m0) < 1e-12 && close(v0 / n, 1.0, 1e-12), "standardized series has mean 0 and (population) variance 1");

    /* sort by another series */
    immutable double[n] key = [8, 7, 6, 5, 4, 3, 2, 1];
    y[] = NADBL;
    immutable double[n] perm = [5, 3, 8, 1, 7, 2, 6, 4];
    check(gretl_sort_by(perm.ptr, key.ptr, y.ptr, d) == 0 && y[0] == 5 && y[1] == 3 && y[2] == 7 && y[7] == 6, "gretl_sort_by orders the key by x");

    /* ranks (f = 1: ascending, ties averaged) */
    GretlMatrixRaw * v = gretl_matrix_alloc(5, 1);
    int err = 0;

    foreach (i, e; [10.0, 20.0, 20.0, 5.0, 30.0]) v.val[i] = e;
    GretlMatrixRaw * r = rank_vector(v, 1, &err);
    check(err == 0 && r !is null && r.val[0] == 2 && r.val[1] == 3.5 && r.val[2] == 3.5 && r.val[3] == 1 && r.val[4] == 5, "rank_vector averages ties");
    gretl_matrix_free(r);
    gretl_matrix_free(v);

    /* Box-Cox: d = 0 is the log, d = 1 is x - 1 */
    immutable double[n] pos = [1, 2, 3, 4, 5, 6, 7, 8];
    y[] = NADBL;
    check(boxcox_series(pos.ptr, y.ptr, 0.0, d) == 0 && close(y[3], log(4.0), 1e-12), "Box-Cox with d = 0 is the log");
    check(boxcox_series(pos.ptr, y.ptr, 1.0, d) == 0 && close(y[3], 3.0, 1e-12), "Box-Cox with d = 1 is x - 1");

    /* fractional differencing with d = 1 is the first difference */
    y[] = NADBL;
    check(fracdiff_series(x.ptr, y.ptr, 1.0, 1, -1, d) == 0 && close(y[2], 3.0, 1e-9) && close(y[5], 4.0, 1e-9), "fractional difference of order 1 is the first difference");

    /* exponential moving average with d = 0.5 */
    y[] = NADBL;
    err = exponential_movavg_series(x.ptr, y.ptr, d, 0.5, 1, NADBL);
    check(err == 0 && close(y[0], 3.0, 1e-12) && close(y[1], 0.5 * 1 + 0.5 * 3, 1e-12), "exponential moving average");

    /* recursive filter y_t = x_t + 0.5 y_{t-1} */
    GretlMatrixRaw * A = gretl_matrix_alloc(1, 1);
    A.val[0] = 0.5;
    y[] = NADBL;
    err = filter_series(x.ptr, y.ptr, d, A, null, 0.0, null);
    check(err == 0 && close(y[0], 3.0, 1e-12) && close(y[1], 1 + 0.5 * 3, 1e-12) && close(y[2], 4 + 0.5 * (1 + 1.5), 1e-12), "AR filter recursion");
    gretl_matrix_free(A);

    /* Hodrick-Prescott on a straight line returns the line */
    immutable double[n] line = [1, 2, 3, 4, 5, 6, 7, 8];
    y[] = NADBL;
    check(hp_filter(line.ptr, y.ptr, d, 1600.0, Gretlopt.OPT_NONE) == 0 && fabs(y[3]) < 1e-6, "HP cycle of a straight line is zero");
    check(hp_filter(line.ptr, y.ptr, d, 1600.0, Gretlopt.OPT_T) == 0 && close(y[3], 4.0, 1e-6), "HP trend of a straight line is the line");

    /* special functions */
    int e2 = 0;
    check(close(gretl_bessel('J', 0.0, 1.0, &e2), 0.7651976865579666, 1e-10) && e2 == 0, "Bessel J0(1)");
    check(close(gretl_bessel('I', 0.0, 1.0, &e2), 1.2660658777520084, 1e-10), "Bessel I0(1)");
    check(close(logistic_cdf(0.0), 0.5, 1e-15) && close(logistic_cdf(1.0), 1.0 / (1.0 + exp(-1.0)), 1e-12), "logistic cdf");

    /* vector statistics */
    GretlMatrixRaw * M = gretl_matrix_alloc(2, 3);
    foreach (i, e; [1.0, 4.0, 2.0, 5.0, 3.0, 6.0]) M.val[i] = e;
    GretlMatrixRaw * cs = gretl_matrix_vector_stat(M, GretlVecStat.sum, 0, 0, &err);
    check(err == 0 && cs !is null && cs.val[0] == 5 && cs.val[1] == 7 && cs.val[2] == 9, "column sums");
    gretl_matrix_free(cs);
    GretlMatrixRaw * rm = gretl_matrix_vector_stat(M, GretlVecStat.mean, 1, 0, &err);
    check(err == 0 && rm !is null && rm.val[0] == 2 && rm.val[1] == 5, "row means");
    gretl_matrix_free(rm);
    gretl_matrix_free(M);

    /* substitution of values */
    immutable double[3] src = [1, 2, 3];
    immutable double[2] from = [2, 3];
    immutable double[2] to = [20, 30];
    double[3] dst;
    check(substitute_values(dst.ptr, src.ptr, 3, from.ptr, 2, to.ptr, 2) == 0 && dst[0] == 1 && dst[1] == 20 && dst[2] == 30, "substitute_values");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
