/* Tests for plugin/kalman_arma.d: exact ML estimation of ARMA models with the gretl-native Kalman filter.

   kalman_arma() is the entry point that arma.c calls; here it is called directly with a hand-made ArmaInfo (the model-level tests
   that go through arma_model are in test_arma.d). The functions it uses from arma.c and arma_init.c are the ported ones.

   The estimates are checked against an independent computation of the exact Gaussian log-likelihood of a stationary ARMA(p, q) model
   from the autocovariance matrix of the sample (psi-weights, then the T x T Toeplitz matrix), concentrated with respect to the
   innovation variance as the Kalman code does: ll = -0.5 (T (1 + log 2 pi + log(q/T)) + log |S|), q = e' S^-1 e, S the matrix of
   autocovariances for unit innovation variance. The test checks that the reported log-likelihood is that of the estimated
   coefficients, that those coefficients are a maximum (the numerical score of the independent likelihood vanishes), that sigma^2 is
   q/T, and that the Hessian-based covariance matrix is the inverse of the numerical Hessian of the independent likelihood.

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs, log;
import core.stdc.stdlib : free;
import core.stdc.stdio : printf, setvbuf, stdout, _IOLBF;

import std.math : PI;

import dgretl;

import gretl.dataset : create_new_dataset, destroy_dataset;
import gretl.gretl_matrixh : GretlMatrixRaw;
import gretl.gretl_list : gretl_list_new;
import gretl.gretl_listh : LISTSEP;
import gretl.gretl_model : clear_model, gretl_model_init;
import gretl.gretl_modelh : ArmaFlags;
import gretl.libgretlh : DatasetRaw, Gretlopt, Model;
import gretl.matrix_extra : gretl_vcv_matrix_from_model, gretl_vector_from_series;
import gretl.plugin.arma_priv;
import gretl.plugin.kalman_arma : kalman_arma;

import gretl.gretl_prn : Prn;

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

ulong lcg = 12345;

double urand()
{
    lcg = lcg * 6364136223846793005UL + 1442695040888963407UL;
    // sum of 4 uniforms, rescaled to unit variance: close enough to normal for a likelihood test
    double s = 0;

    foreach (i; 0 .. 4) {
        lcg = lcg * 6364136223846793005UL + 1442695040888963407UL;
        s += ((lcg >> 11) / cast(double) (1UL << 53)) - 0.5;
    }
    return s * 1.7320508;
}

/* full lag polynomials: phi*(L) = (1 - phi_1 L - ...)(1 - Phi_1 L^s - ...) returned as coefficients c_k with x_t = sum c_k x_{t-k} + ...,
   and theta*(L) = (1 + theta_1 L + ...)(1 + Theta_1 L^s + ...) as coefficients m_k (k >= 1) */
struct Polys { double[] ar, ma; }

Polys polys(double[] phi, double[] Phi, double[] theta, double[] Theta, int pd)
{
    double[] arA = [1.0], maA = [1.0];

    foreach (i, x; phi) { arA.length = i + 2; arA[i + 1] = -x; }
    foreach (i, x; theta) { maA.length = i + 2; maA[i + 1] = x; }
    double[] mul(double[] A, double[] B)
    {
        auto R = new double[A.length + B.length - 1];

        R[] = 0;
        foreach (i, x; A) foreach (j, y; B) R[i + j] += x * y;
        return R;
    }
    double[] arB = [1.0], maB = [1.0];

    foreach (i, x; Phi) { arB.length = (i + 1) * pd + 1; arB[] = 0; arB[0] = 1; arB[(i + 1) * pd] = -x; }
    foreach (i, x; Theta) { maB.length = (i + 1) * pd + 1; maB[] = 0; maB[0] = 1; maB[(i + 1) * pd] = x; }
    auto ar = mul(arA, arB), ma = mul(maA, maB);

    foreach (ref x; ar[1 .. $]) x = -x;      // now x_t = sum_k ar[k] x_{t-k} + ...
    return Polys(ar, ma);
}

/* simulate y_t = mu + phi(L)^-1 theta(L) e_t */
double[] simulate(int T, double mu, Polys pl)
{
    const burn = 300, n = T + burn;
    double[] e = new double[n], x = new double[n];

    x[] = 0;

    foreach (t; 0 .. n) e[t] = urand();
    foreach (t; 0 .. n) {
        double v = 0;

        foreach (k, m; pl.ma) if (t >= k) v += m * e[t - k];
        foreach (k; 1 .. pl.ar.length) if (t >= k) v += pl.ar[k] * x[t - k];
        x[t] = v;
    }
    auto y = new double[T];

    foreach (t; 0 .. T) y[t] = mu + x[burn + t];
    return y;
}

/* the independent concentrated log-likelihood; b = (mu, phi..., Phi..., theta..., Theta...) */
double indepLL(double[] y, int p, int P, int q, int Q, int pd, double[] b, double * q_out = null)
{
    const T = cast(int) y.length;
    auto pl = polys(b[1 .. 1 + p], b[1 + p .. 1 + p + P], b[1 + p + P .. 1 + p + P + q], b[1 + p + P + q .. 1 + p + P + q + Q], pd);
    enum J = 900;
    double[] psi = new double[J];

    foreach (j; 0 .. J) {
        double v = (j < pl.ma.length) ? pl.ma[j] : 0.0;

        foreach (k; 1 .. pl.ar.length) if (j >= k) v += pl.ar[k] * psi[j - k];
        psi[j] = v;
    }
    auto gam = new double[T];

    foreach (k; 0 .. T) {
        double g = 0;

        foreach (j; 0 .. J - k) g += psi[j] * psi[j + k];
        gam[k] = g;
    }
    auto S = GretlMatrix(gretl_matrix_alloc(T, T));

    foreach (i; 0 .. T) foreach (j; 0 .. T) S[i, j] = gam[i > j ? i - j : j - i];
    auto e = GretlMatrix(gretl_matrix_alloc(T, 1));

    foreach (t; 0 .. T) e[t, 0] = y[t] - b[0];
    const qf = (t(e) * inv(S) * e)[0, 0];

    if (q_out !is null) *q_out = qf;
    return -0.5 * (T * (1 + log(2 * PI) + log(qf / T)) + log(det(S)));
}

void runModel(string name, int p, int P, int q, int Q, int pd, double mu, double[] truth, double[] start)
{
    enum T = 150;
    auto pl = polys(truth[0 .. p], truth[p .. p + P], truth[p + P .. p + P + q], truth[p + P + q .. p + P + q + Q], pd);
    auto y = simulate(T, mu, pl);
    auto dset = create_new_dataset(2, T, 0);

    foreach (t; 0 .. T) dset.Z[1][t] = y[t];

    ArmaInfo ainfo;

    ainfo.yno = 1;
    ainfo.flags = ArmaFlags.exact;
    ainfo.ifc = 1;
    ainfo.p = p;
    ainfo.np = p;
    ainfo.P = P;
    ainfo.q = q;
    ainfo.nq = q;
    ainfo.Q = Q;
    ainfo.nc = 1 + p + P + q + Q;
    ainfo.t1 = 0;
    ainfo.t2 = T - 1;
    ainfo.pd = pd;
    ainfo.T = T;
    ainfo.fullT = T;
    ainfo.y = dset.Z[1];
    ainfo.yscale = 1.0;
    ainfo.yshift = 0.0;
    if (P > 0 || Q > 0) set_arma_has_seasonal(&ainfo);

    Model pmod;

    gretl_model_init(&pmod, dset);
    ainfo.alist = gretl_list_new(4);        // p, q, ;, y (what the arma command would have recorded)
    ainfo.alist[1] = ainfo.p;
    ainfo.alist[2] = ainfo.q;
    ainfo.alist[3] = LISTSEP;
    ainfo.alist[4] = 1;
    const err = kalman_arma(start.ptr, dset, &ainfo, &pmod, Gretlopt.OPT_NONE);

    check(err == 0 && pmod.errcode == 0, name ~ ": kalman_arma returns no error");
    if (err == 0) {
        double[] b = new double[ainfo.nc];

        foreach (i; 0 .. ainfo.nc) b[i] = pmod.coeff[i];
        double qf;
        double LL(double[] bb, double * qo = null) { return indepLL(y, p, P, q, Q, pd, bb, qo); }
        const ll = LL(b, &qf);

        check(close(pmod.lnL, ll, 1e-7), name ~ ": reported log-likelihood is the exact ARMA log-likelihood of the estimates");
        check(close(pmod.sigma * pmod.sigma, qf / T, 1e-6), name ~ ": sigma^2 = q/T");

        // the estimates are a maximum of the independent likelihood
        bool zero = true;

        foreach (i; 0 .. b.length) {
            const h = 1e-4;
            auto bp = b.dup, bm = b.dup;

            bp[i] += h;
            bm[i] -= h;
            zero = zero && fabs((LL(bp) - LL(bm)) / (2 * h)) < 1e-3;
        }
        check(zero, name ~ ": the numerical score of the independent likelihood is zero at the estimates");

        bool max = true;

        foreach (i; 0 .. b.length) {
            auto bp = b.dup;

            bp[i] += 0.02;
            max = max && LL(bp) < ll + 1e-9;
            bp[i] -= 0.04;
            max = max && LL(bp) < ll + 1e-9;
        }
        check(max, name ~ ": the likelihood is lower for perturbed coefficients");

        // covariance matrix: inverse of minus the numerical Hessian of the independent likelihood
        int e2 = 0;
        auto V = GretlMatrix(gretl_vcv_matrix_from_model(&pmod, null, &e2));
        const k = b.length;
        auto H = GretlMatrix(gretl_matrix_alloc(cast(int) k, cast(int) k));
        const h = 1e-3;

        foreach (i; 0 .. k) foreach (j; 0 .. k) {
            double f(double di, double dj) { auto bb = b.dup; bb[i] += di; bb[j] += dj; return LL(bb); }

            H[i, j] = -(f(h, h) - f(h, -h) - f(-h, h) + f(-h, -h)) / (4 * h * h);
        }
        auto Vi = inv(H);
        bool vok = e2 == 0 && V.nrow == cast(int) k;

        foreach (i; 0 .. k) foreach (j; 0 .. k) vok = vok && fabs(V[i, j] - Vi[i, j]) < 0.02 * (fabs(Vi[i, j]) + (Vi[i, i] * Vi[j, j]) ^^ 0.5 * 0.1);
        check(vok, name ~ ": Hessian covariance matrix agrees with the numerical Hessian of the independent likelihood");

        // the true parameters are within a few standard errors of the estimates
        bool near = fabs(b[0] - mu) < 5 * (V[0, 0] ^^ 0.5);

        foreach (i; 1 .. k) near = near && fabs(b[i] - truth[i - 1]) < 5 * (V[i, i] ^^ 0.5);
        check(near, name ~ ": estimates are within five standard errors of the true values");
    }
    clear_model(&pmod);
    free(ainfo.alist);
    destroy_dataset(dset);
}

/* ARIMA(1,1,0) by the levels formulation of kalman_arma.c (which exists for ARIMA with within-sample missing values): y_t = y_{t-1} + w_t,
   w_t = phi w_{t-1} + e_t, y_0 known. The observed levels y_obs = S (y_0 + C w) with C the cumulative-sum matrix and S the selection of the
   observed periods, so for unit innovation variance their covariance is Omega = S C Sigma C' S' (Sigma: AR(1) autocovariances), and the
   concentrated log-likelihood is -0.5 (n (1 + log 2 pi + log(q/n)) + log|Omega|), q = (y_obs - y_0)' Omega^-1 (y_obs - y_0). */
double arimaLL(double[] y, size_t[] obs, double phi, double * q_out = null)
{
    const T = cast(int) y.length - 1;                     // y[0] is the pre-sample value
    auto Sig = GretlMatrix(gretl_matrix_alloc(T, T));

    foreach (i; 0 .. T) foreach (j; 0 .. T) {
        const k = i > j ? i - j : j - i;

        Sig[i, j] = phi ^^ k / (1 - phi * phi);
    }
    const n = cast(int) obs.length;
    auto SC = GretlMatrix(gretl_matrix_alloc(n, T));
    auto d = GretlMatrix(gretl_matrix_alloc(n, 1));

    foreach (a, t; obs) {                                 // t in 1 .. T: y_t = y_0 + w_1 + ... + w_t
        foreach (j; 0 .. T) SC[cast(int) a, j] = (j < t) ? 1.0 : 0.0;
        d[cast(int) a, 0] = y[t] - y[0];
    }
    auto Om = SC * Sig * t(SC);
    const qf = (t(d) * inv(Om) * d)[0, 0];

    if (q_out !is null) *q_out = qf;
    return -0.5 * (n * (1 + log(2 * PI) + log(qf / n)) + log(det(Om)));
}

void runArima(string name, size_t[] missing)
{
    enum T = 120;
    const phi0 = 0.5;
    double[] w = new double[T + 1], y = new double[T + 1];

    w[] = 0;

    y[0] = 10.0;
    foreach (t; 1 .. T + 1) {
        w[t] = phi0 * w[t - 1] + urand();
        y[t] = y[t - 1] + w[t];
    }
    size_t[] obs;

    foreach (t; 1 .. T + 1) {
        bool miss = false;

        foreach (m; missing) if (m == t) miss = true;
        if (miss) y[t] = double.nan; else obs ~= t;
    }
    auto dset = create_new_dataset(2, T + 1, 0);

    foreach (t; 0 .. T + 1) dset.Z[1][t] = y[t];

    ArmaInfo ainfo;

    ainfo.yno = 1;
    ainfo.flags = ArmaFlags.exact;
    ainfo.ifc = 0;
    ainfo.p = 1;
    ainfo.np = 1;
    ainfo.d = 1;
    ainfo.nc = 1;
    ainfo.t1 = 1;
    ainfo.t2 = T;
    ainfo.pd = 1;
    ainfo.T = cast(int) obs.length;
    ainfo.fullT = T;
    ainfo.y = dset.Z[1];
    ainfo.yscale = 1.0;
    ainfo.yshift = 0.0;
    set_arima_levels(&ainfo);
    set_arma_is_arima(&ainfo);
    set_arma_na_ok(&ainfo);
    if (missing.length > 0) set_arma_missvals(&ainfo);

    Model pmod;
    double[1] start = [0.2];

    gretl_model_init(&pmod, dset);
    ainfo.alist = gretl_list_new(4);        // p, q, ;, y (what the arma command would have recorded)
    ainfo.alist[1] = ainfo.p;
    ainfo.alist[2] = ainfo.q;
    ainfo.alist[3] = LISTSEP;
    ainfo.alist[4] = 1;
    const err = kalman_arma(start.ptr, dset, &ainfo, &pmod, Gretlopt.OPT_NONE);

    check(err == 0 && pmod.errcode == 0, name ~ ": kalman_arma returns no error");
    if (err == 0) {
        double qf;
        const phi = pmod.coeff[0];
        const ll = arimaLL(y, obs, phi, &qf);
        const h = 1e-4;

        check(close(pmod.lnL, ll, 1e-6), name ~ ": reported log-likelihood is the exact log-likelihood of the observed levels");
        check(fabs((arimaLL(y, obs, phi + h) - arimaLL(y, obs, phi - h)) / (2 * h)) < 1e-3, name ~ ": the estimate is a maximum of the independent likelihood");
        check(close(pmod.sigma * pmod.sigma, qf / obs.length, 1e-5), name ~ ": sigma^2 = q/n");
        check(fabs(phi - phi0) < 0.3, name ~ ": the estimate of phi is near the true value");
    }
    clear_model(&pmod);
    free(ainfo.alist);
    destroy_dataset(dset);
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    //       name           p P q Q pd  mu     true coefficients (phi, Phi, theta, Theta)  start values
    runModel("ARMA(1,1)",    1, 0, 1, 0, 1,  2.0, [0.6, 0.3],            [1.5, 0.3, 0.1]);
    runModel("MA(1)",        0, 0, 1, 0, 1, -1.0, [0.5],                 [-0.5, 0.1]);
    runModel("AR(2)",        2, 0, 0, 0, 1,  0.5, [0.5, 0.25],           [0.3, 0.2, 0.1]);
    runModel("AR(0)xSAR(1)_4", 0, 1, 0, 0, 4, 1.0, [0.5], [0.5, 0.2]);
    runModel("MA(0)xSMA(1)_4", 0, 0, 0, 1, 4, 1.0, [0.5], [0.5, 0.2]);
    runModel("AR(1)xSAR(1)_4", 1, 1, 0, 0, 4, 1.0, [0.5, 0.4], [0.5, 0.2, 0.2]);
    runModel("MA(1)xSMA(1)_4", 0, 0, 1, 1, 4, 1.0, [0.4, 0.4], [0.5, 0.2, 0.2]);
    runModel("ARMA(1,1)x(1,1)_4", 1, 1, 1, 1, 4, 1.0, [0.5, 0.4, 0.3, 0.3], [0.5, 0.2, 0.2, 0.1, 0.1]);

    runArima("ARIMA(1,1,0) levels, complete data", []);
    runArima("ARIMA(1,1,0) levels, missing observations", [30, 31, 77]);

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
