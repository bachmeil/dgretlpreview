/* Smoke tests for nonparam.d: Spearman's rho and Kendall's tau (matrix form), the sign, signed-rank and
   rank-sum tests in quiet mode, sort_pairs_by_x, loess_fit on exactly linear data, and the kernel bandwidth
   rule. Reference values were worked out by hand.

   Not exercised: the printing forms (spearman_rho, kendall_tau, runs_test), Locke's and the VGE tests (random
   shuffling), and the robust loess variant.

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs, pow;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import gretl.dataset : create_new_dataset;
import gretl.gretl_matrix : gretl_matrix_alloc, gretl_matrix_free;
import gretl.gretl_matrixh : GretlMatrixRaw;
import gretl.libgretlh;
import gretl.nonparam;

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

    int err = 0;

    /* Spearman */
    immutable double[8] x8 = [1, 2, 3, 4, 5, 6, 7, 8];
    immutable double[8] sq = [1, 4, 9, 16, 25, 36, 49, 64];
    immutable double[8] rev = [8, 7, 6, 5, 4, 3, 2, 1];
    GretlMatrixRaw * m = spearman_rho_func(x8.ptr, sq.ptr, 8, &err);
    check(err == 0 && m !is null && close(m.val[0], 1.0), "Spearman rho of a monotone relation is 1");
    gretl_matrix_free(m);
    m = spearman_rho_func(x8.ptr, rev.ptr, 8, &err);
    check(err == 0 && m !is null && close(m.val[0], -1.0), "Spearman rho of a reversed relation is -1");
    gretl_matrix_free(m);

    /* Kendall: 15 pairs, 13 concordant, 2 discordant */
    immutable double[6] kx = [1, 2, 3, 4, 5, 6];
    immutable double[6] ky = [1, 3, 2, 4, 6, 5];
    m = kendall_tau_func(kx.ptr, ky.ptr, 6, &err);
    check(err == 0 && m !is null && close(m.val[0], 11.0 / 15.0), "Kendall tau = (13 - 2) / 15");
    gretl_matrix_free(m);

    /* sorting pairs */
    GretlMatrixRaw * X = gretl_matrix_alloc(3, 1);
    GretlMatrixRaw * Y = gretl_matrix_alloc(3, 1);
    X.val[0] = 3; X.val[1] = 1; X.val[2] = 2;
    Y.val[0] = 30; Y.val[1] = 10; Y.val[2] = 20;
    int * order = null;
    err = sort_pairs_by_x(X, Y, &order, null);
    check(err == 0 && X.val[0] == 1 && X.val[1] == 2 && X.val[2] == 3 && Y.val[0] == 10 && Y.val[2] == 30, "sort_pairs_by_x sorts both vectors");
    check(order !is null && order[0] == 1 && order[1] == 2 && order[2] == 0, "sort_pairs_by_x returns the permutation");
    import core.stdc.stdlib : free;
    free(order);
    gretl_matrix_free(X);
    gretl_matrix_free(Y);

    /* two-sample tests, quiet; the series live in a small dataset */
    enum n = 10;
    DatasetRaw * d = create_new_dataset(3, n, 0);
    d.t1 = 0;
    d.t2 = n - 1;
    immutable double[n] a = [1.1, 2.3, 3.2, 4.8, 5.1, 6.7, 7.2, 8.9, 9.4, 10.5];
    immutable double[n] b = [0.5, 1.9, 2.9, 4.1, 4.9, 6.0, 6.9, 8.1, 9.0, 9.9];
    double[2] res = 0;

    err = sign_test(a.ptr, b.ptr, 1, 2, d, res.ptr, Gretlopt.OPT_Q, null);
    check(err == 0 && res[0] == 10 && close(res[1], 1.0 / 1024.0, 1e-9), "sign_test: all 10 differences positive, p = 2^-10");
    err = signed_rank_test(a.ptr, b.ptr, 1, 2, d, res.ptr, Gretlopt.OPT_Q, null);
    check(err == 0 && close(res[0], 2.75209, 1e-5) && close(res[1], 0.00296077, 1e-4), "signed_rank_test: z and p-value");
    err = rank_sum_test(a.ptr, b.ptr, 1, 2, d, res.ptr, Gretlopt.OPT_Q, null);
    check(err == 0 && close(res[0], 0.377964, 1e-5) && close(res[1], 0.352728, 1e-4), "rank_sum_test: z and p-value");

    /* loess on exactly linear data reproduces it */
    enum N = 12;
    GretlMatrixRaw * lx = gretl_matrix_alloc(N, 1);
    GretlMatrixRaw * ly = gretl_matrix_alloc(N, 1);
    foreach (i; 0 .. N) {
        lx.val[i] = i + 1;
        ly.val[i] = 2.0 * (i + 1) + 1.0;
    }
    GretlMatrixRaw * yh = loess_fit(lx, ly, 1, 0.8, Gretlopt.OPT_NONE, &err);
    bool okl = err == 0 && yh !is null;
    if (okl) {
        foreach (i; 0 .. N) okl = okl && close(yh.val[i], ly.val[i], 1e-8);
    }
    check(okl, "loess with d = 1 reproduces linear data");
    gretl_matrix_free(yh);
    gretl_matrix_free(lx);
    gretl_matrix_free(ly);

    /* kernel bandwidth: 0.9 * min(sd, IQR / 1.349) * n^(-1/5) */
    double h = kernel_bandwidth(x8.ptr, 8);
    check(h > 0 && h < 2.0, "kernel_bandwidth is positive and sensible");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
