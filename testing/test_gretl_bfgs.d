/* Smoke tests for gretl_bfgs.d on small smooth problems with known solutions: BFGS (numerical and analytic
   gradient, maximizing and in the NA-criterion error case), Newton-Raphson with an analytic Hessian, the
   numerical Hessian and gradient, Nelder-Mead, the amoeba wrapper and the bracketed root finder.

   Not exercised here: L-BFGS (see test_lbfgsb.d), simulated annealing, the user-level
   functions (user_BFGS, user_NR, deriv_free_optimize, user_fdjac), which need the function-call machinery.

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs, sqrt;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import gretl.gretl_bfgs;
import gretl.gretl_bfgsh : BFGS_CRIT_FUNC;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_matrix : gretl_matrix_alloc, gretl_matrix_free;
import gretl.gretl_matrixh : GretlMatrixRaw;
import gretl.libgretlh;
import gretl.nlsh : CriterionType;

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

/* concave criterion with maximum at (1, 2) */

extern (C) double crit(const(double) * b, void * data)
{
    double x = b[0] - 1, y = b[1] - 2;

    return -(x * x + 3 * y * y + 0.5 * x * y);
}

extern (C) int grad(double * b, double * g, int n, BFGS_CRIT_FUNC f, void * data)
{
    double x = b[0] - 1, y = b[1] - 2;

    g[0] = -(2 * x + 0.5 * y);
    g[1] = -(6 * y + 0.5 * x);
    return 0;
}

extern (C) int hess(double * b, GretlMatrixRaw * H, void * data)
{
    /* positive definite "information" form, as the Newton code expects */
    H.val[0] = 2;
    H.val[1] = 0.5;
    H.val[2] = 0.5;
    H.val[3] = 6;
    return 0;
}

extern (C) double cmin(const(double) * b, void * data)
{
    return -crit(b, data);
}

extern (C) double zf(double x, void * data)
{
    return x * x * x - 2;
}

extern (C) double nacrit(const(double) * b, void * data)
{
    return double.nan;
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    double[2] b = [0, 0];
    int fc, gc;

    int err = BFGS_max(b.ptr, 2, 200, 1e-12, &fc, &gc, &crit, CriterionType.other, null, null, null, Gretlopt.OPT_NONE, null);

    check(err == 0 && close(b[0], 1, 1e-4) && close(b[1], 2, 1e-4), "BFGS with a numerical gradient finds the maximum");
    b = [0, 0];
    err = BFGS_max(b.ptr, 2, 200, 1e-12, &fc, &gc, &crit, CriterionType.other, &grad, null, null, Gretlopt.OPT_NONE, null);
    check(err == 0 && close(b[0], 1, 1e-6) && close(b[1], 2, 1e-6) && fc > 0 && gc > 0, "BFGS with an analytic gradient");

    b = [0, 0];
    int iters;

    err = newton_raphson_max(b.ptr, 2, 100, 1e-12, 1e-9, &iters, CriterionType.other, &crit, &grad, &hess, null, Gretlopt.OPT_NONE, null);
    check(err == 0 && close(b[0], 1, 1e-6) && close(b[1], 2, 1e-6), "Newton-Raphson with an analytic Hessian");
    check(iters >= 1 && iters <= 5, "a quadratic needs very few Newton iterations");

    /* numerical derivatives at the maximum and elsewhere */
    double[2] at = [1.5, 1];
    double[2] g;

    check(numeric_gradient(at.ptr, g.ptr, 2, &crit, null) == 0, "numeric_gradient runs");
    double[2] ga;

    grad(at.ptr, ga.ptr, 2, &crit, null);
    check(close(g[0], ga[0], 1e-4) && close(g[1], ga[1], 1e-4), "numeric gradient matches the analytic one");
    GretlMatrixRaw * H = gretl_matrix_alloc(2, 2);

    check(numerical_hessian(at.ptr, H, &crit, null, 1, 0.01) == 0, "numerical_hessian runs");
    check(close(H.val[0], 2, 1e-3) && close(H.val[3], 6, 1e-3) && close(H.val[1], 0.5, 1e-3), "numerical Hessian (negated) of a quadratic");
    gretl_matrix_free(H);
    GretlMatrixRaw * Hi = numerical_hessian_inverse(at.ptr, 2, &crit, null, 0.01, &err);

    check(err == 0 && close(Hi.val[0] * 2 + Hi.val[1] * 0.5, 1, 1e-3), "numerical_hessian_inverse");
    gretl_matrix_free(Hi);

    /* derivative-free methods */
    double[2] st = [0, 0];
    double[2] xmin;
    double[2] step = [0.5, 0.5];
    double ynew;
    int ncalls, nresets;

    err = nelder_mead(&crit, 2, st.ptr, xmin.ptr, &ynew, 1e-14, step.ptr, 2000, &ncalls, &nresets, null, Gretlopt.OPT_NONE, null);
    check(err == 0 && close(xmin[0], 1, 1e-3) && close(xmin[1], 2, 1e-3), "Nelder-Mead maximizes the criterion");
    b = [0, 0];
    err = gretl_amoeba(b.ptr, 2, 2000, &crit, null, Gretlopt.OPT_NONE, null);
    check(err == 0 && close(b[0], 1, 1e-3) && close(b[1], 2, 1e-3), "gretl_amoeba");

    /* root finder: x^3 = 2 */
    double[2] br = [0, 2];
    double root;

    err = gretl_fzero(br.ptr, 1e-12, &zf, null, &root, Gretlopt.OPT_NONE, null);
    check(err == 0 && close(root, 1.2599210498948732, 1e-8), "gretl_fzero finds the cube root of 2");

    /* error handling */
    b = [0, 0];
    err = BFGS_max(b.ptr, 2, 50, 1e-12, &fc, &gc, &nacrit, CriterionType.other, null, null, null, Gretlopt.OPT_NONE, null);
    check(err != 0, "a criterion that is NA at the start is an error");

    int maxit;
    double tol;

    BFGS_defaults(&maxit, &tol, 0);
    check(maxit > 0, "BFGS_defaults gives a positive iteration limit");

    printf("\n%d checks, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
