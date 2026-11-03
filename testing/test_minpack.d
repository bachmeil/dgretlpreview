/* Tests for gretl/minpack/*.d (enorm, qrfac, qrsolv, lmpar, fdjac2, lmder, lmdif): MINPACK, the Levenberg-Marquardt least-squares routines.

   Each routine is checked against an independent computation with the matrix operators of dgretl or against a problem with a known solution:
   enorm on scaled vectors (over- and underflow), qrfac by R'R = (AP)'(AP) and the column norms and the ordering of the pivots, qrsolv and lmpar
   against the normal equations, fdjac2 against the analytic Jacobian (the three "quality" settings), lmdif and lmder on exponential
   regression with exact data, a linear regression (the OLS solution), and Rosenbrock's function written as two residuals.

   Build and run with testing/run_tests.sh. */

import core.stdc.math : exp, fabs, sqrt, pow;
import core.stdc.stdio : printf, setvbuf, stdout, _IOLBF;

import dgretl;

import gretl.minpack.enorm : enorm_;
import gretl.minpack.fdjac2 : fdjac2_;
import gretl.minpack.lmder : lmder_;
import gretl.minpack.lmdif : lmdif_;
import gretl.minpack.lmpar : lmpar_;
import gretl.minpack.qrfac : qrfac_;
import gretl.minpack.qrsolv : qrsolv_;

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

bool close(double a, double b, double tol)
{
    return fabs(a - b) <= tol * (1 + fabs(b));
}

/* ---- test problems ---- */

/* exponential regression: y_i = a exp(b t_i), residuals f_i = a exp(b t_i) - y_i */
struct ExpData {
    double[] t, y;
}

extern (C) int expFcn(int m, int n, double * x, double * fvec, int * iflag, void * p)
{
    auto d = cast(ExpData *) p;

    foreach (i; 0 .. m) fvec[i] = x[0] * exp(x[1] * d.t[i]) - d.y[i];
    return 0;
}

extern (C) int expFcnJ(int m, int n, double * x, double * fvec, double * fjac, int ldfjac, int * iflag, void * p)
{
    auto d = cast(ExpData *) p;

    if (*iflag == 1) {
        foreach (i; 0 .. m) fvec[i] = x[0] * exp(x[1] * d.t[i]) - d.y[i];
    } else {
        foreach (i; 0 .. m) {
            fjac[i] = exp(x[1] * d.t[i]);                                  // column 0
            fjac[i + ldfjac] = x[0] * d.t[i] * exp(x[1] * d.t[i]);         // column 1
        }
    }
    return 0;
}

/* linear regression y = b0 + b1 x */
struct LinData {
    double[] x, y;
}

extern (C) int linFcn(int m, int n, double * b, double * fvec, int * iflag, void * p)
{
    auto d = cast(LinData *) p;

    foreach (i; 0 .. m) fvec[i] = b[0] + b[1] * d.x[i] - d.y[i];
    return 0;
}

extern (C) int linFcnJ(int m, int n, double * b, double * fvec, double * fjac, int ldfjac, int * iflag, void * p)
{
    auto d = cast(LinData *) p;

    if (*iflag == 1) {
        foreach (i; 0 .. m) fvec[i] = b[0] + b[1] * d.x[i] - d.y[i];
    } else {
        foreach (i; 0 .. m) {
            fjac[i] = 1.0;
            fjac[i + ldfjac] = d.x[i];
        }
    }
    return 0;
}

/* Rosenbrock as residuals */
extern (C) int rosFcn(int m, int n, double * x, double * fvec, int * iflag, void * p)
{
    fvec[0] = 10 * (x[1] - x[0] * x[0]);
    fvec[1] = 1 - x[0];
    return 0;
}

extern (C) int rosFcnJ(int m, int n, double * x, double * fvec, double * fjac, int ldfjac, int * iflag, void * p)
{
    if (*iflag == 1) {
        fvec[0] = 10 * (x[1] - x[0] * x[0]);
        fvec[1] = 1 - x[0];
    } else {
        fjac[0] = -20 * x[0];
        fjac[1] = -1;
        fjac[0 + ldfjac] = 10;
        fjac[1 + ldfjac] = 0;
    }
    return 0;
}

/* a smooth vector function for the Jacobian test: f(x) = (x0^2 + x1, sin(x0) * x1, exp(x0 + x1)) */
extern (C) int vecFcn(int m, int n, double * x, double * fvec, int * iflag, void * p)
{
    import core.stdc.math : sin;

    fvec[0] = x[0] * x[0] + x[1];
    fvec[1] = sin(x[0]) * x[1];
    fvec[2] = exp(x[0] + x[1]);
    return 0;
}

GretlMatrix colMajor(double[] a, int m, int n)
{
    auto M = GretlMatrix(gretl_matrix_alloc(m, n));

    foreach (j; 0 .. n) foreach (i; 0 .. m) M[i, j] = a[i + j * m];
    return M;
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    /* ---- enorm ---- */
    {
        double[2] a = [3.0, 4.0];
        double[3] big = [1e200, 1e200, 1e200];
        double[3] tiny = [1e-200, 1e-200, 1e-200];
        double[4] z = [0, 0, 0, 0];
        double[5] mixed = [1e-170, 2.0, -3.0, 1e150, 1e-5];

        check(enorm_(2, a.ptr) == 5.0, "enorm of (3, 4) is 5");
        check(close(enorm_(3, big.ptr), sqrt(3.0) * 1e200, 1e-14), "enorm of three 1e200 does not overflow");
        check(close(enorm_(3, tiny.ptr), sqrt(3.0) * 1e-200, 1e-14), "enorm of three 1e-200 does not underflow");
        check(enorm_(4, z.ptr) == 0.0 && enorm_(0, a.ptr) == 0.0, "enorm of zeros and of an empty vector is 0");
        check(close(enorm_(5, mixed.ptr), 1e150, 1e-14), "enorm of a vector with a huge and tiny components is the huge component");
    }

    /* ---- qrfac ---- */
    enum M = 6, N = 3;
    double[M * N] A = [1.0, 2, 3, 4, 5, 6,   2, 1, 0, 1, 2, 3,   10, 1, 4, 2, 7, 5];     // column-major, the third column has the largest norm
    double[N] rdiag, acnorm;
    int[N] ipvt;
    double[N] wa;
    double[M * N] QR = A;

    {
        qrfac_(M, N, QR.ptr, M, ipvt.ptr, rdiag.ptr, acnorm.ptr, wa.ptr);
        // the norms of the columns of A
        bool nrm = true;

        foreach (j; 0 .. N) {
            double s = 0;

            foreach (i; 0 .. M) s += A[i + j * M] ^^ 2;
            nrm = nrm && close(acnorm[j], sqrt(s), 1e-14);
        }
        check(nrm, "qrfac: acnorm holds the column norms of A");
        check(ipvt[0] == 2 && fabs(rdiag[0]) >= fabs(rdiag[1]) && fabs(rdiag[1]) >= fabs(rdiag[2]), "qrfac: pivoting puts the largest column first and |diag R| is non-increasing");
        // R: upper triangle of QR with rdiag on the diagonal; check R'R = (AP)'(AP)
        auto R = GretlMatrix(gretl_matrix_alloc(N, N));

        foreach (i; 0 .. N) foreach (j; 0 .. N) R[i, j] = (i < j) ? QR[i + j * M] : (i == j ? rdiag[i] : 0.0);
        auto AP = GretlMatrix(gretl_matrix_alloc(M, N));

        foreach (j; 0 .. N) foreach (i; 0 .. M) AP[i, j] = A[i + ipvt[j] * M];
        auto lhs = t(R) * R, rhs = t(AP) * AP;
        bool ok = true;

        foreach (i; 0 .. N) foreach (j; 0 .. N) ok = ok && close(lhs[i, j], rhs[i, j], 1e-12);
        check(ok, "qrfac: R'R = (AP)'(AP) for the permuted matrix");
    }

    /* ---- qrsolv and lmpar, from the factorization above ---- */
    {
        double[M] b = [1.0, 0.5, -1, 2, 3, 0.25];
        double[N] diag = [1.0, 2.0, 0.5];
        // Q'b = R^-T (AP)' b
        auto R = GretlMatrix(gretl_matrix_alloc(N, N));

        foreach (i; 0 .. N) foreach (j; 0 .. N) R[i, j] = (i < j) ? QR[i + j * M] : (i == j ? rdiag[i] : 0.0);
        auto AP = GretlMatrix(gretl_matrix_alloc(M, N));

        foreach (j; 0 .. N) foreach (i; 0 .. M) AP[i, j] = A[i + ipvt[j] * M];
        auto bv = GretlMatrix(gretl_matrix_alloc(M, 1));

        foreach (i; 0 .. M) bv[i, 0] = b[i];
        auto qtbv = inv(t(R)) * (t(AP) * bv);
        double[N] qtb;

        foreach (i; 0 .. N) qtb[i] = qtbv[i, 0];

        // qrsolv solves min ||A x - b||^2 + ||D x||^2 in the permuted variables: (AP'AP + P'D^2P) z = (AP)'b, x = P z
        double[N * N] r;                                                       // upper triangle of R, column-major, ldr = N
        double[N] x, sdiag;

        foreach (j; 0 .. N) foreach (i; 0 .. N) r[i + j * N] = (i <= j) ? ((i == j) ? rdiag[i] : QR[i + j * M]) : 0.0;
        double[N] wa2;

        qrsolv_(N, r.ptr, N, ipvt.ptr, diag.ptr, qtb.ptr, x.ptr, sdiag.ptr, wa2.ptr);
        auto Am = colMajor(A[], M, N);                                            // the original matrix; x comes back in the original variable order (x = P z)
        auto Dm = GretlMatrix(gretl_matrix_alloc(N, N));

        foreach (i; 0 .. N) foreach (j; 0 .. N) Dm[i, j] = (i == j) ? diag[i] ^^ 2 : 0.0;
        auto z = inv(t(Am) * Am + Dm) * (t(Am) * bv);
        bool ok = true;

        foreach (j; 0 .. N) ok = ok && close(x[j], z[j, 0], 1e-10);
        check(ok, "qrsolv: the least-squares solution with the diagonal D equals the normal-equations solution");

        // lmpar: find the Levenberg-Marquardt parameter for a trust-region radius delta
        double[N * N] r2;
        double[N] x2, sdiag2, wa1, wa3;
        double par = 0;
        const delta = 0.3;

        foreach (j; 0 .. N) foreach (i; 0 .. N) r2[i + j * N] = (i <= j) ? ((i == j) ? rdiag[i] : QR[i + j * M]) : 0.0;
        lmpar_(N, r2.ptr, N, ipvt.ptr, diag.ptr, qtb.ptr, delta, &par, x2.ptr, sdiag2.ptr, wa1.ptr, wa3.ptr);
        // x2 is the solution in the original variable order: ||D x|| is within 10% of delta unless par = 0
        double dxn = 0;

        foreach (j; 0 .. N) dxn += (diag[j] * x2[j]) ^^ 2;
        dxn = sqrt(dxn);
        check(par > 0 && fabs(dxn - delta) <= 0.1 * delta + 1e-12, "lmpar: the parameter makes ||D x|| within 10% of the radius (the full step is longer than the radius here)");
        // x2 solves (A'A + par D^2) x = A'b
        auto Dp = GretlMatrix(gretl_matrix_alloc(N, N));

        foreach (i; 0 .. N) foreach (j; 0 .. N) Dp[i, j] = (i == j) ? par * diag[i] ^^ 2 : 0.0;
        auto zz = inv(t(Am) * Am + Dp) * (t(Am) * bv);
        bool ok2 = true;

        foreach (j; 0 .. N) ok2 = ok2 && close(x2[j], zz[j, 0], 1e-8);
        check(ok2, "lmpar: x solves the damped normal equations with the parameter it found");
        // a large radius: par = 0 and x is the Gauss-Newton step
        double par0 = 0;
        double[N] x3, sd3, w1, w3;
        double[N * N] r3;

        foreach (j; 0 .. N) foreach (i; 0 .. N) r3[i + j * N] = (i <= j) ? ((i == j) ? rdiag[i] : QR[i + j * M]) : 0.0;
        lmpar_(N, r3.ptr, N, ipvt.ptr, diag.ptr, qtb.ptr, 100.0, &par0, x3.ptr, sd3.ptr, w1.ptr, w3.ptr);
        auto gn = inv(t(Am) * Am) * (t(Am) * bv);
        bool ok3 = par0 == 0;

        foreach (j; 0 .. N) ok3 = ok3 && close(x3[j], gn[j, 0], 1e-9);
        check(ok3, "lmpar: for a large radius the parameter is 0 and x is the Gauss-Newton step");
    }

    /* ---- fdjac2 ---- */
    {
        double[2] x = [0.7, -0.4];
        double[3] fv;
        int iflag = 1;

        vecFcn(3, 2, x.ptr, fv.ptr, &iflag, null);
        double[6] J;
        double[3] wa3;
        // analytic Jacobian, column-major 3 x 2
        import core.stdc.math : cos, sin;
        double[6] Ja = [2 * x[0], cos(x[0]) * x[1], exp(x[0] + x[1]),   1.0, sin(x[0]), exp(x[0] + x[1])];
        double[3] errs;

        foreach (q; 0 .. 3) {
            iflag = 1;
            fdjac2_(&vecFcn, 3, 2, q, x.ptr, fv.ptr, J.ptr, 3, &iflag, 0.0, wa3.ptr, null);
            double e = 0;

            foreach (k; 0 .. 6) e = e > fabs(J[k] - Ja[k]) ? e : fabs(J[k] - Ja[k]);
            errs[q] = e;
        }
        check(errs[0] < 1e-6 && errs[1] < 1e-7 && errs[2] < 1e-8, "fdjac2: forward, central and five-point differences approximate the analytic Jacobian");
        check(errs[2] <= errs[0], "fdjac2: the highest quality is at least as accurate as the plain forward difference");
    }

    /* ---- the Levenberg-Marquardt routines ---- */
    {
        // exponential regression with exact data a = 2.5, b = 0.4
        ExpData ed;

        foreach (i; 0 .. 10) {
            ed.t ~= 0.2 * i;
            ed.y ~= 2.5 * exp(0.4 * 0.2 * i);
        }
        enum m = 10, n = 2;
        double[n] x = [1.0, 0.1], diag;
        double[m] fvec;
        double[m * n] fjac;
        int[n] ip;
        double[n] qtf, w1, w2, w3;
        double[m] w4;
        int info = 0, nfev = 0, njev = 0;

        lmdif_(&expFcn, m, n, x.ptr, fvec.ptr, 1e-12, 1e-12, 0.0, 400, 0.0, diag.ptr, 1, 100.0, 0, &info, &nfev, fjac.ptr, m, ip.ptr, qtf.ptr, w1.ptr, w2.ptr, w3.ptr, w4.ptr, &ed);
        check(info >= 1 && info <= 4 && close(x[0], 2.5, 1e-7) && close(x[1], 0.4, 1e-7), "lmdif: exponential regression recovers (a, b) = (2.5, 0.4)");
        double ss = 0;

        foreach (i; 0 .. m) ss += fvec[i] ^^ 2;
        check(ss < 1e-14 && nfev > 0 && nfev < 400, "lmdif: the residuals vanish and the number of evaluations is reasonable");

        double[n] xj = [1.0, 0.1], diag2;
        double[m] fvec2;
        double[m * n] fjac2;
        int[n] ip2;
        double[n] qtf2, a1, a2, a3;
        double[m] a4;

        lmder_(&expFcnJ, m, n, xj.ptr, fvec2.ptr, fjac2.ptr, m, 1e-12, 1e-12, 0.0, 400, diag2.ptr, 1, 100.0, 0, &info, &nfev, &njev, ip2.ptr, qtf2.ptr, a1.ptr, a2.ptr, a3.ptr, a4.ptr, &ed);
        check(info >= 1 && info <= 4 && close(xj[0], 2.5, 1e-8) && close(xj[1], 0.4, 1e-8) && njev > 0, "lmder: the same problem with the analytic Jacobian");
        check(nfev <= 400 && close(xj[0], x[0], 1e-6) && close(xj[1], x[1], 1e-6), "lmder and lmdif reach the same solution");

        // linear regression: the least-squares solution is the OLS estimate
        LinData ld;
        double[] noise = [0.3, -0.1, 0.25, -0.4, 0.05, 0.2, -0.35, 0.15, -0.2, 0.1, 0.0, -0.05];

        foreach (i; 0 .. 12) {
            ld.x ~= i * 0.5;
            ld.y ~= 1.0 + 2.0 * ld.x[i] + noise[i];
        }
        auto X = GretlMatrix(gretl_matrix_alloc(12, 2));
        auto Y = GretlMatrix(gretl_matrix_alloc(12, 1));

        foreach (i; 0 .. 12) { X[i, 0] = 1; X[i, 1] = ld.x[i]; Y[i, 0] = ld.y[i]; }
        auto ols = inv(t(X) * X) * (t(X) * Y);
        double[2] b = [0.0, 0.0], d3;
        double[12] fv3;
        double[24] fj3;
        int[2] ip3;
        double[2] q3, c1, c2, c3;
        double[12] c4;
        int info3, nf3, nj3;

        lmdif_(&linFcn, 12, 2, b.ptr, fv3.ptr, 1e-14, 1e-14, 0.0, 400, 0.0, d3.ptr, 1, 100.0, 0, &info3, &nf3, fj3.ptr, 12, ip3.ptr, q3.ptr, c1.ptr, c2.ptr, c3.ptr, c4.ptr, &ld);
        check(close(b[0], ols[0, 0], 1e-7) && close(b[1], ols[1, 0], 1e-7), "lmdif on a linear regression gives the OLS estimates");
        double[2] bj = [0.0, 0.0];

        lmder_(&linFcnJ, 12, 2, bj.ptr, fv3.ptr, fj3.ptr, 12, 1e-14, 1e-14, 0.0, 400, d3.ptr, 1, 100.0, 0, &info3, &nf3, &nj3, ip3.ptr, q3.ptr, c1.ptr, c2.ptr, c3.ptr, c4.ptr, &ld);
        check(close(bj[0], ols[0, 0], 1e-9) && close(bj[1], ols[1, 0], 1e-9), "lmder on a linear regression gives the OLS estimates");

        // Rosenbrock from the classic start (-1.2, 1)
        double[2] xr = [-1.2, 1.0], dr;
        double[2] fvr, wr1, wr2, wr3, qr_;
        double[2] wr4;
        double[4] fjr;
        int[2] ipr;
        int ir, nfr, njr;

        lmdif_(&rosFcn, 2, 2, xr.ptr, fvr.ptr, 1e-12, 1e-12, 0.0, 600, 0.0, dr.ptr, 1, 100.0, 0, &ir, &nfr, fjr.ptr, 2, ipr.ptr, qr_.ptr, wr1.ptr, wr2.ptr, wr3.ptr, wr4.ptr, null);
        check(close(xr[0], 1.0, 1e-7) && close(xr[1], 1.0, 1e-7), "lmdif finds the minimum (1, 1) of Rosenbrock's function");
        double[2] xs = [-1.2, 1.0];

        lmder_(&rosFcnJ, 2, 2, xs.ptr, fvr.ptr, fjr.ptr, 2, 1e-12, 1e-12, 0.0, 600, dr.ptr, 1, 100.0, 0, &ir, &nfr, &njr, ipr.ptr, qr_.ptr, wr1.ptr, wr2.ptr, wr3.ptr, wr4.ptr, null);
        check(close(xs[0], 1.0, 1e-8) && close(xs[1], 1.0, 1e-8), "lmder finds the minimum (1, 1) of Rosenbrock's function");

        // limits: too few evaluations are reported through info = 5 (maxfev reached)
        double[2] xl = [-1.2, 1.0];

        lmdif_(&rosFcn, 2, 2, xl.ptr, fvr.ptr, 1e-14, 1e-14, 0.0, 5, 0.0, dr.ptr, 1, 100.0, 0, &ir, &nfr, fjr.ptr, 2, ipr.ptr, qr_.ptr, wr1.ptr, wr2.ptr, wr3.ptr, wr4.ptr, null);
        check(ir == 5 && nfr >= 5, "lmdif: stopping at the evaluation limit gives info = 5");
        // invalid input: n > m gives info = 0
        double[2] xi = [1.0, 1.0];
        double[1] fv1;
        int[2] ipi;

        lmdif_(&rosFcn, 1, 2, xi.ptr, fv1.ptr, 1e-8, 1e-8, 0.0, 100, 0.0, dr.ptr, 1, 100.0, 0, &ir, &nfr, fjr.ptr, 1, ipi.ptr, qr_.ptr, wr1.ptr, wr2.ptr, wr3.ptr, wr4.ptr, null);
        check(ir == 0, "lmdif: n > m is an input error (info = 0)");
    }

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
