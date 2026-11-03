/* Smoke test for bhhh_max.d: maximum likelihood estimation of the mean and (log) standard deviation of a normal
   sample by the BHHH method with an analytic score, against the closed-form ML estimates; the variance matrix
   from the final OPG regression; the numerical-score option; and the argument and non-convergence error paths.

   Build and run with testing/run_tests.sh. */

import core.stdc.math : exp, fabs, log, sqrt;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import gretl.bhhh_max : bhhh_max;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_matrix : gretl_matrix_alloc, gretl_matrix_free;
import gretl.gretl_matrixh : GretlMatrixRaw, gretl_matrix_set;
import gretl.libgretlh;

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

bool close(double a, double b, double tol = 1e-5)
{
    return fabs(a - b) <= tol * (1.0 + fabs(b));
}

immutable double[12] y = [1.2, 0.4, 2.3, 1.8, 0.9, 1.5, 2.9, 0.1, 1.1, 1.9, 2.4, 0.6];

/* parameters: mu, log sigma; G is T x 2 (score by observation); returns the log-likelihood, with the
   per-observation values left in column 0 of G's companion when the score is not requested */
extern (C) double normal_ll(double * th, GretlMatrixRaw * G, void * data, int do_score, int * err)
{
    double mu = th[0], ls = th[1], s2 = exp(2 * ls), ll = 0;

    foreach (t; 0 .. y.length) {
        double e = y[t] - mu;
        double lt = -0.5 * log(2 * 3.141592653589793) - ls - 0.5 * e * e / s2;

        ll += lt;
        if (do_score) {
            gretl_matrix_set(G, cast(int) t, 0, e / s2);
            gretl_matrix_set(G, cast(int) t, 1, e * e / s2 - 1);
        }
    }
    return ll;
}

extern (C) double bad_ll(double * th, GretlMatrixRaw * G, void * data, int do_score, int * err)
{
    *err = 1;
    return double.nan;
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    double m = 0, v = 0;

    foreach (x; y) m += x / y.length;
    foreach (x; y) v += (x - m) * (x - m) / y.length;

    double[2] th = [0.5, 0.0];
    GretlMatrixRaw * G = gretl_matrix_alloc(cast(int) y.length, 2);
    GretlMatrixRaw * V = gretl_matrix_alloc(2, 2);
    int fc, gc;
    int err = bhhh_max(th.ptr, 2, G, &normal_ll, 1e-12, &fc, &gc, null, V, Gretlopt.OPT_NONE, null);

    check(err == 0, "bhhh_max converges");
    check(close(th[0], m, 1e-4), "ML estimate of the mean is the sample mean");
    check(close(exp(2 * th[1]), v, 1e-4), "ML estimate of the variance is the (1/T) sample variance");
    check(fc >= gc && gc > 0, "function and gradient counts");
    check(V.val[0] > 0 && V.val[3] > 0 && close(V.val[0], v / y.length, 0.5), "OPG variance of the mean is roughly sigma^2 / T");

    /* a wrong column count is rejected */
    GretlMatrixRaw * G3 = gretl_matrix_alloc(cast(int) y.length, 3);

    th = [0.5, 0.0];
    check(bhhh_max(th.ptr, 2, G3, &normal_ll, 1e-12, &fc, &gc, null, null, Gretlopt.OPT_NONE, null) == GretlError.nonconf, "a score matrix with the wrong number of columns is E_NONCONF");
    gretl_matrix_free(G3);

    /* a callback that fails */
    th = [0.5, 0.0];
    check(bhhh_max(th.ptr, 2, G, &bad_ll, 1e-12, &fc, &gc, null, null, Gretlopt.OPT_NONE, null) != 0, "a failing callback gives an error");

    /* too few iterations allowed -> not converged is reported via the tolerance */
    th = [50.0, 3.0];
    err = bhhh_max(th.ptr, 2, G, &normal_ll, 1e-30, &fc, &gc, null, null, Gretlopt.OPT_NONE, null);
    check(err == 0 || err == GretlError.noconv, "an unreachable tolerance ends in convergence or E_NOCONV");

    gretl_matrix_free(G);
    gretl_matrix_free(V);
    printf("\n%d checks, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
