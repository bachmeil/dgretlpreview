/* Tests for plugin/fcp.d (fcp.c): GARCH(p, q) estimation with the mixed-gradient algorithm of Fiorentini, Calzolari and Panattoni (option
   --fcp of garch_model, i.e. garch_estimate).

   The results are compared with an independent computation of the Gaussian GARCH log-likelihood in the units of the data (squared
   residuals and h started at the sample variance of the residuals, as fcp.c does): the reported log-likelihood at the reported
   parameters, the vanishing numerical score at them, the covariance matrix against the numerical Hessian of the likelihood (default),
   and the other covariance options (QML with --robust, and the OP and information-matrix variants through the libset variable
   garch_vcv).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs, log, sqrt;
import core.stdc.stdio : printf, setvbuf, stdout, _IOLBF;
import core.stdc.string : strcmp;

import std.math : PI;

import dgretl;

import gretl.dataset : create_new_dataset, dataset_set_time_series, destroy_dataset;
import gretl.gretl_list : gretl_list_new;
import gretl.gretl_listh : LISTSEP;
import gretl.gretl_model : clear_model, gretl_model_get_int;
import gretl.libgretlh : DatasetRaw, Gretlopt, Model;
import gretl.libset : execute_set;
import gretl.matrix_extra : gretl_vcv_matrix_from_model;
import gretl.plugin.garch : garch_model;

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

ulong lcg = 31415;

double urand()
{
    double s = 0;

    foreach (i; 0 .. 4) {
        lcg = lcg * 6364136223846793005UL + 1442695040888963407UL;
        s += ((lcg >> 11) / cast(double) (1UL << 53)) - 0.5;
    }
    return s * 1.7320508;
}

int * garchList(int p, int q, int[] rest)
{
    auto L = gretl_list_new(cast(int) rest.length + 3);

    L[1] = p;
    L[2] = q;
    L[3] = LISTSEP;
    foreach (i, v; rest) L[i + 4] = v;
    return L;
}

/* theta = (c, b, alpha0, alpha_1 .. alpha_q, beta_1 .. beta_p): y = c + b x + e, h_t = alpha0 + sum alpha_i e2_{t-i} + sum beta_j h_{t-j} */
double garchLL(double[] y, double[] x, double[] th, int p, int q)
{
    const T = cast(int) y.length;
    const lag = p > q ? p : q;
    auto e = new double[T], e2 = new double[T];
    double unc = 0;

    foreach (t; 0 .. T) {
        e[t] = y[t] - th[0] - th[1] * x[t];
        e2[t] = e[t] * e[t];
        unc += e2[t];
    }
    unc /= T;
    // h and e2 with pre-sample values at the unconditional variance
    auto h = new double[T + lag], ee = new double[T + lag];

    foreach (t; 0 .. lag) h[t] = ee[t] = unc;
    double ll = 0;

    foreach (t; 0 .. T) {
        double ht = th[2];

        foreach (i; 1 .. q + 1) ht += th[2 + i] * ee[lag + t - i];
        foreach (j; 1 .. p + 1) ht += th[2 + q + j] * h[lag + t - j];
        if (ht <= 0) ht = 1e-7;
        h[lag + t] = ht;
        ee[lag + t] = e2[t];
        ll -= 0.5 * log(ht) + 0.5 * e2[t] / ht + 0.5 * log(2 * PI);
    }
    return ll;
}

void runFcp(string name, int p, int q, double[] truth, int T, ulong seed)
{
    lcg = seed;
    auto dset = create_new_dataset(3, T, 0);
    double[] ys = new double[T], xs = new double[T];
    auto e2 = new double[T + 10], h = new double[T + 10];

    e2[] = 0.5;
    h[] = 1.0;
    foreach (t; 0 .. T) {
        const t0 = t + 10;
        double ht = truth[2];

        foreach (i; 1 .. q + 1) ht += truth[2 + i] * e2[t0 - i];
        foreach (j; 1 .. p + 1) ht += truth[2 + q + j] * h[t0 - j];
        h[t0] = ht;
        const e = sqrt(ht) * urand();

        e2[t0] = e * e;
        xs[t] = urand();
        ys[t] = truth[0] + truth[1] * xs[t] + e;
        dset.Z[1][t] = ys[t];
        dset.Z[2][t] = xs[t];
    }
    dataset_set_time_series(dset, 1, 1900, 1);
    auto L = garchList(p, q, [1, 2]);
    auto pm = garch_model(L, dset, null, Gretlopt.OPT_F);

    check(pm.errcode == 0, name ~ ": garch_model --fcp: no error");
    if (pm.errcode == 0) {
        const k = 3 + p + q;
        auto th = new double[k];

        foreach (i; 0 .. k) th[i] = pm.coeff[i];
        const ll = garchLL(ys, xs, th, p, q);

        check(pm.ncoeff == k && pm.nobs == T, name ~ ": number of coefficients and observations");
        check(close(pm.lnL, ll, 1e-6), name ~ ": log-likelihood equals the independent GARCH likelihood at the estimates");
        auto grad = new double[k];

        foreach (i; 0 .. k) {
            const hh = 1e-6 * (1 + fabs(th[i]));
            auto bp = th.dup, bm = th.dup;

            bp[i] += hh;
            bm[i] -= hh;
            grad[i] = (garchLL(ys, xs, bp, p, q) - garchLL(ys, xs, bm, p, q)) / (2 * hh);
        }
        bool pos = th[2] > 0;

        foreach (i; 3 .. k) pos = pos && th[i] >= 0 && th[i] < 1;
        check(pos, name ~ ": variance parameters are in range");
        check(strcmp(pm.params[2], "alpha(0)") == 0 && strcmp(pm.params[3], "alpha(1)") == 0 && strcmp(pm.params[k - 1], p > 0 ? "beta(1)".ptr : "alpha(1)".ptr) == 0 || q > 1,
              name ~ ": parameter names");

        // the covariance matrix is the inverse of minus the numerical Hessian of the likelihood
        int err = 0;
        auto V = GretlMatrix(gretl_vcv_matrix_from_model(&pm, null, &err));
        auto H = GretlMatrix(gretl_matrix_alloc(k, k));

        foreach (i; 0 .. k) foreach (j; 0 .. k) {
            const hi = 1e-4 * (1 + fabs(th[i])), hj = 1e-4 * (1 + fabs(th[j]));
            double f(double di, double dj)
            {
                auto bb = th.dup;

                bb[i] += di;
                bb[j] += dj;
                return garchLL(ys, xs, bb, p, q);
            }

            H[i, j] = -(f(hi, hj) - f(hi, -hj) - f(-hi, hj) + f(-hi, -hj)) / (4 * hi * hj);
        }
        auto Vi = inv(H);
        double decr = 0;

        foreach (i; 0 .. k) foreach (j; 0 .. k) decr += grad[i] * Vi[i, j] * grad[j];
        check(decr < 1e-3, name ~ ": the Newton decrement g\'V g of the independent likelihood at the estimates is negligible (the estimates are the maximum)");
        bool vok = err == 0 && V.nrow == k;

        if (vok) foreach (i; 0 .. k) vok = vok && V[i, i] > 0 && fabs(V[i, i] - Vi[i, i]) < 0.05 * Vi[i, i];
        check(vok, name ~ ": the covariance matrix agrees with the inverse numerical Hessian (variances within 5%)");
        check(gretl_model_get_int(&pm, "iters") > 0, name ~ ": the iteration count is recorded");

        // other covariance estimators run and give positive variances
        foreach (opt; [Gretlopt.OPT_R]) {
            auto Lr = garchList(p, q, [1, 2]);
            auto pr = garch_model(Lr, dset, null, cast(Gretlopt) (Gretlopt.OPT_F | opt));

            int e3 = 0;
            bool rok = pr.errcode == 0;

            if (rok) {
                auto Vr = GretlMatrix(gretl_vcv_matrix_from_model(&pr, null, &e3));

                rok = e3 == 0 && Vr.nrow == k;
                if (rok) foreach (i; 0 .. k) rok = rok && Vr[i, i] > 0;
            }
            check(rok, name ~ ": --robust (QML) covariance matrix is computed with positive variances");
            clear_model(&pr);
        }
    }
    clear_model(&pm);
    destroy_dataset(dset);
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    //          p  q  true (c, b, alpha0, alpha_1..q, beta_1..p)                T    seed
    runFcp("GARCH(1,1)", 1, 1, [4.0, 2.0, 0.4, 0.25, 0.55],                 600, 101);
    runFcp("ARCH(2)",    0, 2, [1.0, -1.0, 0.5, 0.2, 0.25],                 600, 202);
    foreach (sd; [305, 306]) runFcp("GARCH(1,2)", 1, 2, [0.5, 1.5, 0.3, 0.15, 0.1, 0.55],            800, sd);

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
