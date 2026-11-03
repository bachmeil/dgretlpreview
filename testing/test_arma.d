/* Tests for plugin/arma.d (arma.c with arma_common.c): the ARMA / ARIMA model driver arma_model(), its list and option checking, the
   model statistics, parameter names and roots it writes, and its registration in the plugin table.

   Only the exact-ML route through the Kalman filter (option --kalman, OPT_K) with user-supplied initial values can be run end to end for
   now: the other routes need arma_init.c (initialization, conditional least squares), arma_cml.c (conditional ML by BHHH) or as_driver.c
   (AS 154 / AS 197), which are not ported yet. arma_select() runs arma_model() for each order and so waits for the same files.

   The results are compared with an independent exact Gaussian likelihood of a stationary ARMA(p, q) model, computed from the
   autocovariance matrix of the sample with the matrix operators of dgretl (see test_kalman_arma.d); ARIMA(p,1,q) is the ARMA model for the
   first differences. 

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs, log, sqrt;
import core.stdc.stdio : printf, setvbuf, stdout, _IOLBF;
import core.stdc.string : strcmp, strcpy;

import std.math : PI;

import dgretl;

import gretl.dataset : create_new_dataset, dataset_set_time_series, destroy_dataset;
import gretl.gretl_list : gretl_list_new;
import gretl.gretl_matrix : gretl_matrix_free;
import gretl.gretl_matrixh : GretlMatrixRaw, gretl_matrix_get;
import gretl.gretl_listh : LISTSEP;
import gretl.gretl_model : clear_model, gretl_model_get_int;
import gretl.gretl_modelh : ArmaFlags;
import gretl.libgretlh : DatasetRaw, Gretlopt, GretlType, Model;
import gretl.libset : execute_set;
import gretl.matrix_extra : gretl_vcv_matrix_from_model;
import gretl.plugin.arma : arma_model, gui_arma_select;
import gretl.plugin.arma_priv : ArmaInfo;
import gretl.plugins : get_plugin_function;
import gretl.uservar : user_var_add_or_replace, user_var_delete_by_name;



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

ulong lcg = 4242;

double urand()
{
    double s = 0;

    foreach (i; 0 .. 4) {
        lcg = lcg * 6364136223846793005UL + 1442695040888963407UL;
        s += ((lcg >> 11) / cast(double) (1UL << 53)) - 0.5;
    }
    return s * 1.7320508;
}

/* x_t = phi x_{t-1} + e_t + theta e_{t-1} */
double[] simulateArma(int T, double phi, double theta)
{
    const burn = 200, n = T + burn;
    double[] e = new double[n], x = new double[n];

    x[] = 0;
    foreach (t; 0 .. n) {
        e[t] = urand();
        x[t] = e[t] + ((t > 0) ? theta * e[t - 1] + phi * x[t - 1] : 0);
    }
    return x[burn .. $].dup;
}

/* independent concentrated exact log-likelihood of ARMA(p, q) for zero-mean data x (general orders) */
double indepLLg(double[] x, double[] phi, double[] theta, double * qout = null)
{
    const T = cast(int) x.length;
    enum J = 700;
    double[] psi = new double[J];

    foreach (j; 0 .. J) {
        double v = (j == 0) ? 1.0 : (j <= theta.length ? theta[j - 1] : 0.0);

        foreach (i, ph; phi) if (j > i) v += ph * psi[j - i - 1];
        psi[j] = v;
    }
    auto S = GretlMatrix(gretl_matrix_alloc(T, T));
    auto gam = new double[T];

    foreach (k; 0 .. T) {
        double g = 0;

        foreach (j; 0 .. J - k) g += psi[j] * psi[j + k];
        gam[k] = g;
    }
    foreach (i; 0 .. T) foreach (j; 0 .. T) S[i, j] = gam[i > j ? i - j : j - i];
    auto e = GretlMatrix(gretl_matrix_alloc(T, 1));

    foreach (t; 0 .. T) e[t, 0] = x[t];
    const qf = (t(e) * inv(S) * e)[0, 0];

    if (qout !is null) *qout = qf;
    return -0.5 * (T * (1 + log(2 * PI) + log(qf / T)) + log(det(S)));
}

/* independent concentrated exact log-likelihood of ARMA(1,1) for zero-mean data x */
double indepLL(double[] x, double phi, double theta, double * qout = null)
{
    const T = cast(int) x.length;
    enum J = 600;
    double[] psi = new double[J];

    psi[0] = 1;
    foreach (j; 1 .. J) psi[j] = (j == 1 ? theta : 0) + phi * psi[j - 1];
    auto S = GretlMatrix(gretl_matrix_alloc(T, T));
    auto gam = new double[T];

    foreach (k; 0 .. T) {
        double g = 0;

        foreach (j; 0 .. J - k) g += psi[j] * psi[j + k];
        gam[k] = g;
    }
    foreach (i; 0 .. T) foreach (j; 0 .. T) S[i, j] = gam[i > j ? i - j : j - i];
    auto e = GretlMatrix(gretl_matrix_alloc(T, 1));

    foreach (t; 0 .. T) e[t, 0] = x[t];
    const qf = (t(e) * inv(S) * e)[0, 0];

    if (qout !is null) *qout = qf;
    return -0.5 * (T * (1 + log(2 * PI) + log(qf / T)) + log(det(S)));
}

int * armaList(int[] orders, int yno)
{
    auto L = gretl_list_new(cast(int) orders.length + 2);

    foreach (i, o; orders) L[i + 1] = o;
    L[orders.length + 1] = LISTSEP;
    L[orders.length + 2] = yno;
    return L;
}

void setInitvals(double[] v)
{
    auto m = gretl_matrix_alloc(cast(int) v.length, 1);

    foreach (i, x; v) m.val[i] = x;
    user_var_add_or_replace("iv_test", GretlType.matrix, m);
    execute_set("initvals", "iv_test", null, Gretlopt.OPT_NONE, null);
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    enum T = 140;
    const mu = 3.0, phi0 = 0.6, theta0 = 0.3;
    auto x = simulateArma(T, phi0, theta0);
    double[] y = new double[T];

    foreach (t; 0 .. T) y[t] = mu + x[t];
    auto dset = create_new_dataset(2, T, 0);

    foreach (t; 0 .. T) dset.Z[1][t] = y[t];
    dataset_set_time_series(dset, 1, 1900, 1);

    /* ---- ARMA(1,1) with a constant: exact ML by the Kalman filter, user initial values ---- */
    {
        auto L = armaList([1, 1], 1);

        setInitvals([2.5, 0.4, 0.1]);
        auto pm = arma_model(L, null, dset, cast(Gretlopt) (Gretlopt.OPT_K | Gretlopt.OPT_U), null);

        check(pm.errcode == 0, "arma_model --kalman with initvals: no error");
        if (pm.errcode == 0) {
            check(pm.ncoeff == 3 && pm.nobs == T && pm.t1 == 0 && pm.t2 == T - 1, "ncoeff, nobs and the sample range");
            check(pm.ci == 0 || pm.ci > 0, "model command index set");
            // the log-likelihood and the estimates agree with the independent exact likelihood
            double[] xs = new double[T];

            foreach (t; 0 .. T) xs[t] = y[t] - pm.coeff[0];
            double qf;
            const ll = indepLL(xs, pm.coeff[1], pm.coeff[2], &qf);

            check(close(pm.lnL, ll, 1e-6), "log-likelihood equals the exact ARMA(1,1) likelihood at the estimates");
            check(close(pm.sigma * pm.sigma, qf / T, 1e-6), "sigma^2 = q/T");
            bool zero = true;

            foreach (i; 0 .. 3) {
                const h = 1e-4;
                double f(double d)
                {
                    double[3] b = [pm.coeff[0], pm.coeff[1], pm.coeff[2]];
                    double[] xx = new double[T];

                    b[i] += d;
                    foreach (t; 0 .. T) xx[t] = y[t] - b[0];
                    return indepLL(xx, b[1], b[2]);
                }

                zero = zero && fabs((f(h) - f(-h)) / (2 * h)) < 1e-3;
            }
            check(zero, "the estimates are a maximum of the independent likelihood");
            check(pm.nparams == 3 && strcmp(pm.params[0], "const") == 0 && strcmp(pm.params[1], "phi_1") == 0 && strcmp(pm.params[2], "theta_1") == 0,
                  "parameter names const, phi_1, theta_1");
            check(gretl_model_get_int(&pm, "arma_flags") == ArmaFlags.exact, "arma_flags records exact ML");
            check(pm.criterion[0] == pm.criterion[0] && pm.criterion[1] == pm.criterion[1], "information criteria are computed");
            int err = 0;
            auto V = GretlMatrix(gretl_vcv_matrix_from_model(&pm, null, &err));

            check(err == 0 && V.nrow == 3 && V[0, 0] > 0 && V[1, 1] > 0 && V[2, 2] > 0, "covariance matrix of the coefficients");
            check(fabs(pm.coeff[1] - phi0) < 5 * sqrt(V[1, 1]) && fabs(pm.coeff[2] - theta0) < 5 * sqrt(V[2, 2]), "estimates are near the true values");
            check(pm.uhat !is null && pm.uhat[pm.t1] == pm.uhat[pm.t1], "residuals stored");
        }
        clear_model(&pm);
    }

    /* ---- exact ML by the Kalman filter with the program's own initialization (no initvals) ---- */
    foreach (spec; [[1, 1], [1, 2], [2, 1]]) {
        const p = spec[0], q = spec[1];
        auto L = armaList([p, q], 1);
        auto pm = arma_model(L, null, dset, Gretlopt.OPT_K, null);
        const name = "ARMA(" ~ cast(char) ('0' + p) ~ "," ~ cast(char) ('0' + q) ~ ") --kalman, default initialization";

        check(pm.errcode == 0, name ~ ": no error");
        if (pm.errcode == 0) {
            auto b = new double[pm.ncoeff];

            foreach (i; 0 .. pm.ncoeff) b[i] = pm.coeff[i];
            double LL(double[] bb)
            {
                auto xx = new double[T];

                foreach (t; 0 .. T) xx[t] = y[t] - bb[0];
                return indepLLg(xx, bb[1 .. 1 + p], bb[1 + p .. 1 + p + q]);
            }
            check(pm.ncoeff == 1 + p + q && close(pm.lnL, LL(b), 1e-6), name ~ ": log-likelihood equals the exact likelihood at the estimates");
            bool zero = true;

            foreach (i; 0 .. b.length) {
                const h = 1e-4;
                auto bp = b.dup, bm = b.dup;

                bp[i] += h;
                bm[i] -= h;
                zero = zero && fabs((LL(bp) - LL(bm)) / (2 * h)) < 2e-3;
            }
            check(zero, name ~ ": the estimates are a maximum of the independent likelihood");
        }
        clear_model(&pm);
    }

    /* ---- conditional ML for a pure AR(2) with a constant: OLS ---- */
    {
        auto L = armaList([2, 0], 1);
        auto pm = arma_model(L, null, dset, Gretlopt.OPT_C, null);

        check(pm.errcode == 0, "AR(2) --conditional: no error");
        if (pm.errcode == 0) {
            const n = T - 2;
            auto X = GretlMatrix(gretl_matrix_alloc(n, 3));
            auto yy = GretlMatrix(gretl_matrix_alloc(n, 1));

            foreach (t; 2 .. T) {
                X[t - 2, 0] = 1;
                X[t - 2, 1] = y[t - 1];
                X[t - 2, 2] = y[t - 2];
                yy[t - 2, 0] = y[t];
            }
            auto beta = inv(t(X) * X) * (t(X) * yy);
            auto u = yy - X * beta;
            const ess = (t(u) * u)[0, 0];
            
            check(pm.nobs == n && pm.t1 == 2 && pm.ncoeff == 3, "AR(2) --conditional: sample and number of coefficients");
            check(close(pm.coeff[1], beta[1, 0], 1e-8) && close(pm.coeff[2], beta[2, 0], 1e-8), "AR coefficients are the OLS ones");
            check(close(pm.coeff[0], beta[0, 0], 1e-8), "the constant is the OLS intercept (conditional least squares)");
            check(close(pm.ess, ess, 1e-8), "residual sum of squares");
            check(gretl_model_get_int(&pm, "arma_flags") == ArmaFlags.ls, "arma_flags records conditional (least squares)");
        }
        clear_model(&pm);
    }

    /* ---- exact ML by AS 197 (the default route, arma_model without --kalman) ---- */
    foreach (spec; [[1, 1], [2, 1], [1, 2]]) {
        const p = spec[0], q = spec[1];
        auto L = armaList([p, q], 1);
        auto pm = arma_model(L, null, dset, Gretlopt.OPT_NONE, null);
        auto pk = arma_model(L, null, dset, Gretlopt.OPT_K, null);
        const name = "ARMA(" ~ cast(char) ('0' + p) ~ "," ~ cast(char) ('0' + q) ~ ") exact ML by AS 197";

        check(pm.errcode == 0 && pk.errcode == 0, name ~ ": no error");
        if (pm.errcode == 0 && pk.errcode == 0) {
            auto b = new double[pm.ncoeff];

            foreach (i; 0 .. pm.ncoeff) b[i] = pm.coeff[i];
            double LL(double[] bb)
            {
                auto xx = new double[T];

                foreach (t; 0 .. T) xx[t] = y[t] - bb[0];
                return indepLLg(xx, bb[1 .. 1 + p], bb[1 + p .. 1 + p + q]);
            }
            check(pm.ncoeff == 1 + p + q && close(pm.lnL, LL(b), 1e-6), name ~ ": log-likelihood equals the exact likelihood at the estimates");
            bool zero = true;

            foreach (i; 0 .. b.length) {
                const h = 1e-4;
                auto bp = b.dup, bm = b.dup;

                bp[i] += h;
                bm[i] -= h;
                zero = zero && fabs((LL(bp) - LL(bm)) / (2 * h)) < 2e-3;
            }
            check(zero, name ~ ": the estimates are a maximum of the independent likelihood");
            bool same = close(pm.lnL, pk.lnL, 1e-6);

            foreach (i; 0 .. b.length) same = same && fabs(pm.coeff[i] - pk.coeff[i]) < 2e-3 * (1 + fabs(pk.coeff[i]));
            check(same, name ~ ": same estimates and log-likelihood as the Kalman filter");
            int e1 = 0, e2 = 0;
            auto V1 = GretlMatrix(gretl_vcv_matrix_from_model(&pm, null, &e1));
            auto V2 = GretlMatrix(gretl_vcv_matrix_from_model(&pk, null, &e2));
            bool vok = e1 == 0 && e2 == 0 && V1.nrow == V2.nrow;

            if (vok) foreach (i; 0 .. V1.nrow) vok = vok && fabs(V1[i, i] - V2[i, i]) < 0.05 * V2[i, i];
            check(vok, name ~ ": the Hessian variances agree with those of the Kalman route");
            check(gretl_model_get_int(&pm, "arma_flags") == ArmaFlags.exact, name ~ ": arma_flags records exact ML");
        }
        clear_model(&pm);
        clear_model(&pk);
    }

    /* AS 197 for a seasonal model and an ARMAX model, with initial values */
    {
        const phi1 = 0.4, Phi1 = 0.5;
        double[] xs = new double[T + 200], ys = new double[T];

        xs[] = 0;
        foreach (t; 0 .. T + 200) {
            double v = urand();

            if (t >= 1) v += phi1 * xs[t - 1];
            if (t >= 4) v += Phi1 * xs[t - 4];
            if (t >= 5) v -= phi1 * Phi1 * xs[t - 5];
            xs[t] = v;
        }
        foreach (t; 0 .. T) ys[t] = 1.0 + xs[200 + t];
        auto ds = create_new_dataset(2, T, 0);

        foreach (t; 0 .. T) ds.Z[1][t] = ys[t];
        dataset_set_time_series(ds, 4, 1900, 1);
        auto L = gretl_list_new(7);

        L[1] = 1; L[2] = 0; L[3] = LISTSEP; L[4] = 1; L[5] = 0; L[6] = LISTSEP; L[7] = 1;
        setInitvals([1.0, 0.3, 0.3]);
        auto pm = arma_model(L, null, ds, Gretlopt.OPT_NONE, null);
        setInitvals([1.0, 0.3, 0.3]);
        auto pk = arma_model(L, null, ds, cast(Gretlopt) (Gretlopt.OPT_K | Gretlopt.OPT_U), null);

        check(pm.errcode == 0 && pk.errcode == 0, "seasonal ARMA(1,0)x(1,0)_4 by AS 197: no error");
        if (pm.errcode == 0 && pk.errcode == 0) {
            bool same = close(pm.lnL, pk.lnL, 1e-6);

            foreach (i; 0 .. 3) same = same && fabs(pm.coeff[i] - pk.coeff[i]) < 2e-3 * (1 + fabs(pk.coeff[i]));
            check(same, "seasonal ARMA: AS 197 and the Kalman filter agree");
        }
        clear_model(&pm);
        clear_model(&pk);
        destroy_dataset(ds);

        auto xr = new double[T], yr = new double[T];
        double u = 0;

        foreach (t; 0 .. T) {
            xr[t] = urand();
            u = 0.5 * u + urand();
            yr[t] = 2.0 + 1.5 * xr[t] + u;
        }
        auto dx = create_new_dataset(3, T, 0);

        foreach (t; 0 .. T) { dx.Z[1][t] = yr[t]; dx.Z[2][t] = xr[t]; }
        dataset_set_time_series(dx, 1, 1900, 1);
        auto L2 = gretl_list_new(6);

        L2[1] = 1; L2[2] = 0; L2[3] = LISTSEP; L2[4] = 1; L2[5] = 0; L2[6] = 2;
        setInitvals([1.5, 0.3, 1.0]);
        auto pa = arma_model(L2, null, dx, Gretlopt.OPT_NONE, null);
        setInitvals([1.5, 0.3, 1.0]);
        auto pb = arma_model(L2, null, dx, cast(Gretlopt) (Gretlopt.OPT_K | Gretlopt.OPT_U), null);

        check(pa.errcode == 0 && pb.errcode == 0, "ARMAX(1,0) by AS 197: no error");
        if (pa.errcode == 0 && pb.errcode == 0) {
            bool same = close(pa.lnL, pb.lnL, 1e-6);

            foreach (i; 0 .. 3) same = same && fabs(pa.coeff[i] - pb.coeff[i]) < 2e-3 * (1 + fabs(pb.coeff[i]));
            check(same, "ARMAX: AS 197 and the Kalman filter agree");
        }
        clear_model(&pa);
        clear_model(&pb);
        destroy_dataset(dx);
    }

    /* ---- exact ML by AS 154 (option --as154 = OPT_A, and the default for samples with missing observations) ---- */
    foreach (spec; [[1, 1], [2, 1], [1, 2]]) {
        const p = spec[0], q = spec[1];
        auto L = armaList([p, q], 1);
        auto pa = arma_model(L, null, dset, Gretlopt.OPT_A, null);
        auto pk = arma_model(L, null, dset, Gretlopt.OPT_K, null);
        const name = "ARMA(" ~ cast(char) ('0' + p) ~ "," ~ cast(char) ('0' + q) ~ ") exact ML by AS 154";

        check(pa.errcode == 0 && pk.errcode == 0, name ~ ": no error");
        if (pa.errcode == 0 && pk.errcode == 0) {
            auto b = new double[pa.ncoeff];

            foreach (i; 0 .. pa.ncoeff) b[i] = pa.coeff[i];
            double LL(double[] bb)
            {
                auto xx = new double[T];

                foreach (t; 0 .. T) xx[t] = y[t] - bb[0];
                return indepLLg(xx, bb[1 .. 1 + p], bb[1 + p .. 1 + p + q]);
            }
            check(close(pa.lnL, LL(b), 1e-6), name ~ ": log-likelihood equals the exact likelihood at the estimates");
            bool zero = true;

            foreach (i; 0 .. b.length) {
                const h = 1e-4;
                auto bp = b.dup, bm = b.dup;

                bp[i] += h;
                bm[i] -= h;
                zero = zero && fabs((LL(bp) - LL(bm)) / (2 * h)) < 2e-3;
            }
            check(zero, name ~ ": the estimates are a maximum of the independent likelihood");
            bool same = close(pa.lnL, pk.lnL, 1e-6);

            foreach (i; 0 .. b.length) same = same && fabs(pa.coeff[i] - pk.coeff[i]) < 2e-3 * (1 + fabs(pk.coeff[i]));
            check(same, name ~ ": same estimates and log-likelihood as the Kalman filter");
        }
        clear_model(&pa);
        clear_model(&pk);
    }

    /* OPG matrix through as154_llt_callback, and the sample with missing observations (AS 154 is then the default) */
    {
        auto L = armaList([1, 1], 1);
        auto pa = arma_model(L, null, dset, cast(Gretlopt) (Gretlopt.OPT_A | Gretlopt.OPT_G), null);
        auto pk = arma_model(L, null, dset, cast(Gretlopt) (Gretlopt.OPT_K | Gretlopt.OPT_G), null);

        check(pa.errcode == 0 && pk.errcode == 0, "ARMA(1,1) by AS 154 with --opg: no error");
        if (pa.errcode == 0 && pk.errcode == 0) {
            int e1 = 0, e2 = 0;
            auto V1 = GretlMatrix(gretl_vcv_matrix_from_model(&pa, null, &e1));
            auto V2 = GretlMatrix(gretl_vcv_matrix_from_model(&pk, null, &e2));
            bool vok = e1 == 0 && e2 == 0 && V1.nrow == 3 && V2.nrow == 3;

            if (vok) foreach (i; 0 .. 3) vok = vok && V1[i, i] > 0 && fabs(V1[i, i] - V2[i, i]) < 0.1 * V2[i, i];
            check(vok, "the AS 154 OPG variances agree with the Kalman ones");
        }
        clear_model(&pa);
        clear_model(&pk);

        auto dm = create_new_dataset(2, T, 0);

        foreach (t; 0 .. T) dm.Z[1][t] = y[t];
        dm.Z[1][30] = double.nan;
        dm.Z[1][31] = double.nan;
        dm.Z[1][77] = double.nan;
        dataset_set_time_series(dm, 1, 1900, 1);
        auto pm = arma_model(L, null, dm, Gretlopt.OPT_NONE, null);             // missing values: AS 197 not usable, so AS 154
        auto pq = arma_model(L, null, dm, Gretlopt.OPT_K, null);

        check(pm.errcode == 0 && pq.errcode == 0, "ARMA(1,1) with missing observations: no error (AS 154 and Kalman)");
        if (pm.errcode == 0 && pq.errcode == 0) {
            bool same = close(pm.lnL, pq.lnL, 1e-6) && pm.nobs == T - 3;

            foreach (i; 0 .. 3) same = same && fabs(pm.coeff[i] - pq.coeff[i]) < 5e-3 * (1 + fabs(pq.coeff[i]));
            check(same, "missing observations: AS 154 and the Kalman filter agree (log-likelihood, estimates, nobs)");
        }
        clear_model(&pm);
        clear_model(&pq);
        destroy_dataset(dm);
    }

    /* ---- OPG covariance matrix (option --opg) through the AS 197 and Kalman callbacks ---- */
    {
        auto L = armaList([1, 1], 1);
        auto pm = arma_model(L, null, dset, Gretlopt.OPT_G, null);
        auto pk = arma_model(L, null, dset, cast(Gretlopt) (Gretlopt.OPT_K | Gretlopt.OPT_G), null);

        check(pm.errcode == 0 && pk.errcode == 0, "ARMA(1,1) with --opg: no error (AS 197 and Kalman)");
        if (pm.errcode == 0 && pk.errcode == 0) {
            int e1 = 0, e2 = 0;
            auto V1 = GretlMatrix(gretl_vcv_matrix_from_model(&pm, null, &e1));
            auto V2 = GretlMatrix(gretl_vcv_matrix_from_model(&pk, null, &e2));
            bool vok = e1 == 0 && e2 == 0 && V1.nrow == 3 && V2.nrow == 3;

            if (vok) foreach (i; 0 .. 3) vok = vok && V1[i, i] > 0 && fabs(V1[i, i] - V2[i, i]) < 0.1 * V2[i, i];
            check(vok, "the OPG variances of the two routes agree");
        }
        clear_model(&pm);
        clear_model(&pk);
    }

    /* ---- conditional ML by BHHH (arma_cml.d): ARMA(1,1) and MA(1) with a constant ---- */
    foreach (spec; [[1, 1], [0, 1]]) {
        const p = spec[0], q = spec[1];
        auto L = armaList([p, q], 1);
        auto pm = arma_model(L, null, dset, Gretlopt.OPT_C, null);
        const name = "ARMA(" ~ cast(char) ('0' + p) ~ "," ~ cast(char) ('0' + q) ~ ") --conditional";

        check(pm.errcode == 0, name ~ ": no error");
        if (pm.errcode == 0) {
            // e_t = y_t - c - phi y_{t-1} - theta e_{t-1}, starting at observation p, with e_{p-1} = 0
            double cll(double[] b, double * ssrOut = null)
            {
                double ssr = 0, eprev = 0;
                const n = T - p;

                foreach (t; p .. T) {
                    double e = y[t] - b[0] - (p > 0 ? b[1] * y[t - 1] : 0) - b[p + 1] * eprev;

                    ssr += e * e;
                    eprev = e;
                }
                if (ssrOut !is null) *ssrOut = ssr;
                return -0.5 * n * (log(2 * PI) + 1 + log(ssr / n));
            }
            auto b = new double[pm.ncoeff];

            foreach (i; 0 .. pm.ncoeff) b[i] = pm.coeff[i];
            double ssr;
            const ll = cll(b, &ssr);

            check(pm.ncoeff == 2 + p && pm.nobs == T - p && pm.t1 == p, name ~ ": coefficients and sample");
            check(close(pm.lnL, ll, 1e-8), name ~ ": log-likelihood equals the conditional likelihood at the estimates");
            bool zero = true;

            foreach (i; 0 .. b.length) {
                const h = 1e-5;
                auto bp = b.dup, bm = b.dup;

                bp[i] += h;
                bm[i] -= h;
                zero = zero && fabs((cll(bp) - cll(bm)) / (2 * h)) < 1e-2;
            }
            check(zero, name ~ ": the estimates maximize the conditional likelihood");
            check(close(pm.sigma * pm.sigma, ssr / (T - p), 1e-6), name ~ ": sigma^2 = SSR/n");
            check(gretl_model_get_int(&pm, "arma_flags") == ArmaFlags.ls || gretl_model_get_int(&pm, "arma_flags") == 0, "arma_flags is not exact");
            int err = 0;
            auto V = GretlMatrix(gretl_vcv_matrix_from_model(&pm, null, &err));

            check(err == 0 && V.nrow == pm.ncoeff && V[0, 0] > 0 && V[p + 1, p + 1] > 0, name ~ ": covariance matrix of the coefficients");
        }
        clear_model(&pm);
    }

    /* ---- conditional ML with a seasonal MA term: (1,0)x(0,1)_4 ---- */
    {
        auto ds = create_new_dataset(2, T, 0);

        foreach (t; 0 .. T) ds.Z[1][t] = y[t];
        dataset_set_time_series(ds, 4, 1900, 1);
        auto L = gretl_list_new(7);

        L[1] = 1; L[2] = 0; L[3] = LISTSEP; L[4] = 0; L[5] = 1; L[6] = LISTSEP; L[7] = 1;      // arma 1 0 ; 0 1 ; y
        auto pm = arma_model(L, null, ds, Gretlopt.OPT_C, null);

        check(pm.errcode == 0, "ARMA(1,0)x(0,1)_4 --conditional: no error");
        if (pm.errcode == 0) {
            double cll(double[] b)
            {
                double ssr = 0;
                auto e = new double[T];

                e[] = 0;
                foreach (t; 1 .. T) {
                    e[t] = y[t] - b[0] - b[1] * y[t - 1] - (t - 4 >= 1 ? b[2] * e[t - 4] : 0);
                    ssr += e[t] * e[t];
                }
                const n = T - 1;

                return -0.5 * n * (log(2 * PI) + 1 + log(ssr / n));
            }
            auto b = new double[pm.ncoeff];

            foreach (i; 0 .. pm.ncoeff) b[i] = pm.coeff[i];
            check(pm.ncoeff == 3 && close(pm.lnL, cll(b), 1e-8), "seasonal MA: log-likelihood equals the conditional likelihood at the estimates");
            bool zero = true;

            foreach (i; 0 .. b.length) {
                const h = 1e-5;
                auto bp = b.dup, bm = b.dup;

                bp[i] += h;
                bm[i] -= h;
                zero = zero && fabs((cll(bp) - cll(bm)) / (2 * h)) < 5e-2;
            }
            foreach (i; 0 .. b.length) {
                auto bp = b.dup;

                bp[i] += 0.01;
                zero = zero && cll(bp) < cll(b) + 1e-9;
                bp[i] -= 0.02;
                zero = zero && cll(bp) < cll(b) + 1e-9;
            }
            check(zero, "seasonal MA: the estimates maximize the conditional likelihood (to the BHHH tolerance)");
        }
        clear_model(&pm);
        destroy_dataset(ds);
    }

    /* ---- selection of the AR order by information criteria (conditional) ---- */
    {
        auto L = armaList([3, 0], 1);                 // p = 0..3, q = 0
        GretlMatrixRaw * sel = null;
        auto err = gui_arma_select(L, dset, cast(Gretlopt) (Gretlopt.OPT_C | Gretlopt.OPT_Q), null, &sel);

        check(err == 0 && sel !is null, "gui_arma_select returns the criteria matrix");
        if (sel !is null) {
            check(sel.rows == 4 && sel.cols == 6, "one row per order, columns p q AIC BIC HQC loglik");
            bool ok = true;

            dset.t1 = 3;                                // the common sample of the selection: the first max-order observations are lost
            foreach (i; 0 .. 4) {
                auto Li = armaList([i, 0], 1);
                auto pm = arma_model(Li, null, dset, cast(Gretlopt) (Gretlopt.OPT_C | Gretlopt.OPT_Q), null);

                ok = ok && pm.errcode == 0 && close(gretl_matrix_get(sel, i, 5), pm.lnL, 1e-8) && close(gretl_matrix_get(sel, i, 2), pm.criterion[0], 1e-8);
                clear_model(&pm);
            }
            dset.t1 = 0;
            check(ok, "each row agrees with arma_model of that order on the common sample (log-likelihood, AIC)");
            gretl_matrix_free(sel);
        }
    }

    /* ---- seasonal ARMA and ARMAX by exact ML (Kalman), default initialization ---- */
    {
        // seasonal AR: (1 - phi L)(1 - Phi L^4) x_t = e_t, quarterly data
        const phi1 = 0.4, Phi1 = 0.5;
        double[] xs = new double[T + 200], ys = new double[T];

        xs[] = 0;
        foreach (t; 0 .. T + 200) {
            double v = urand();

            if (t >= 1) v += phi1 * xs[t - 1];
            if (t >= 4) v += Phi1 * xs[t - 4] - phi1 * Phi1 * xs[t - 5 < 0 ? 0 : t - 5] * (t >= 5 ? 1 : 0);
            xs[t] = v;
        }
        foreach (t; 0 .. T) ys[t] = 1.0 + xs[200 + t];
        auto ds = create_new_dataset(2, T, 0);

        foreach (t; 0 .. T) ds.Z[1][t] = ys[t];
        dataset_set_time_series(ds, 4, 1900, 1);
        auto L = gretl_list_new(7);

        L[1] = 1; L[2] = 0; L[3] = LISTSEP; L[4] = 1; L[5] = 0; L[6] = LISTSEP; L[7] = 1;      // arma 1 0 ; 1 0 ; y
        strcpy(ds.varname[1], "y");
        // the default initialization of this model is by NLS: a genr expression for the regression function, parsed by genlex / gensyntax, evaluated
        // by the ported evaluator and minimized by MINPACK
        auto pm = arma_model(L, null, ds, Gretlopt.OPT_K, null);

        check(pm.errcode == 0, "ARMA(1,0)x(1,0)_4 --kalman with the default (NLS) initialization: no error");
        if (pm.errcode == 0) {
            auto b = new double[pm.ncoeff];

            foreach (i; 0 .. pm.ncoeff) b[i] = pm.coeff[i];
            // full AR polynomial (1 - phi L)(1 - Phi L^4) = 1 - phi L - Phi L^4 + phi Phi L^5
            double[] phiFull = [b[1], 0, 0, b[2], -b[1] * b[2]];
            double LL(double[] bb)
            {
                auto xx = new double[T];

                foreach (t; 0 .. T) xx[t] = ys[t] - bb[0];
                return indepLLg(xx, [bb[1], 0, 0, bb[2], -bb[1] * bb[2]], []);
            }
            check(pm.ncoeff == 3 && close(pm.lnL, LL(b), 1e-6), "seasonal ARMA: log-likelihood equals the exact likelihood at the estimates");
            bool zero = true;

            foreach (i; 0 .. b.length) {
                const h = 1e-4;
                auto bp = b.dup, bm = b.dup;

                bp[i] += h;
                bm[i] -= h;
                zero = zero && fabs((LL(bp) - LL(bm)) / (2 * h)) < 2e-3;
            }
            check(zero, "seasonal ARMA: the estimates are a maximum of the independent likelihood");
            check(strcmp(pm.params[2], "Phi_1") == 0, "seasonal parameter name Phi_1");
        }
        clear_model(&pm);
        destroy_dataset(ds);

        // ARMAX: y_t = c + b x_t + u_t, u an AR(1)
        auto xr = new double[T], yr = new double[T];
        double u = 0;

        foreach (t; 0 .. T) {
            xr[t] = urand();
            u = 0.5 * u + urand();
            yr[t] = 2.0 + 1.5 * xr[t] + u;
        }
        auto dx = create_new_dataset(3, T, 0);

        foreach (t; 0 .. T) { dx.Z[1][t] = yr[t]; dx.Z[2][t] = xr[t]; }
        dataset_set_time_series(dx, 1, 1900, 1);
        auto L2 = gretl_list_new(6);

        L2[1] = 1; L2[2] = 0; L2[3] = LISTSEP; L2[4] = 1; L2[5] = 0; L2[6] = 2;       // arma 1 0 ; y 0 x  (the 0 is the constant)
        strcpy(dx.varname[1], "y");
        strcpy(dx.varname[2], "x");
        auto pa = arma_model(L2, null, dx, Gretlopt.OPT_K, null);

        check(pa.errcode == 0, "ARMAX(1,0) --kalman with the default (NLS) initialization: no error");
        if (pa.errcode == 0) {
            auto b = new double[pa.ncoeff];

            foreach (i; 0 .. pa.ncoeff) b[i] = pa.coeff[i];       // const, phi, beta
            double LL(double[] bb)
            {
                auto xx = new double[T];

                foreach (t; 0 .. T) xx[t] = yr[t] - bb[0] - bb[2] * xr[t];
                return indepLLg(xx, [bb[1]], []);
            }
            check(pa.ncoeff == 3 && close(pa.lnL, LL(b), 1e-6), "ARMAX: log-likelihood equals the exact likelihood at the estimates");
            bool zero = true;

            foreach (i; 0 .. b.length) {
                const h = 1e-4;
                auto bp = b.dup, bm = b.dup;

                bp[i] += h;
                bm[i] -= h;
                zero = zero && fabs((LL(bp) - LL(bm)) / (2 * h)) < 2e-3;
            }
            check(zero, "ARMAX: the estimates are a maximum of the independent likelihood");
            check(fabs(b[2] - 1.5) < 0.3, "ARMAX: the regressor coefficient is near its true value");
        }
        clear_model(&pa);
        destroy_dataset(dx);
    }

    /* ---- ARIMA(1,1,0) without a constant: exact ML by the Kalman filter on the first differences ---- */
    {
        double[] z = new double[T];

        z[0] = 50;
        foreach (t; 1 .. T) z[t] = z[t - 1] + x[t];          // levels whose differences are the ARMA(1,1) series x
        auto d2 = create_new_dataset(2, T, 0);

        foreach (t; 0 .. T) d2.Z[1][t] = z[t];
        dataset_set_time_series(d2, 1, 1900, 1);
        auto L = armaList([1, 1, 0], 1);

        setInitvals([0.3]);
        auto pm = arma_model(L, null, d2, cast(Gretlopt) (Gretlopt.OPT_K | Gretlopt.OPT_N | Gretlopt.OPT_U), null);   // --nc: no constant; OPT_U: use initvals

        check(pm.errcode == 0, "ARIMA(1,1,0): no error");
        if (pm.errcode == 0) {
            double[] w = new double[T - 1];

            foreach (t; 1 .. T) w[t - 1] = z[t] - z[t - 1];
            check(pm.ncoeff == 1 && pm.nobs == T - 1, "one coefficient on T - 1 differenced observations");
            double qf;
            const ll = indepLL(w, pm.coeff[0], 0.0, &qf);

            check(close(pm.lnL, ll, 1e-6), "log-likelihood equals the exact AR(1) likelihood of the differences");
            const h = 1e-4;
            check(fabs((indepLL(w, pm.coeff[0] + h, 0) - indepLL(w, pm.coeff[0] - h, 0)) / (2 * h)) < 1e-3, "the estimate maximizes the independent likelihood");
            check(gretl_model_get_int(&pm, "arima_d") == 1, "arima_d recorded");
        }
        clear_model(&pm);
        destroy_dataset(d2);
    }

    /* ---- list and option checking ---- */
    {
        auto L1 = armaList([200, 1], 1);
        auto pm = arma_model(L1, null, dset, Gretlopt.OPT_K, null);

        check(pm.errcode != 0, "an AR order above the maximum is an error");
        clear_model(&pm);
        auto L2 = armaList([-1, 1], 1);
        pm = arma_model(L2, null, dset, Gretlopt.OPT_K, null);
        check(pm.errcode != 0, "a negative order is an error");
        clear_model(&pm);
        auto L3 = armaList([1, 3, 1], 1);              // d = 3 > 2
        pm = arma_model(L3, null, dset, Gretlopt.OPT_K, null);
        check(pm.errcode != 0, "an integration order above 2 is an error");
        clear_model(&pm);
        auto L4 = armaList([1, 1], 1);
        pm = arma_model(L4, null, dset, cast(Gretlopt) (Gretlopt.OPT_C | Gretlopt.OPT_X | Gretlopt.OPT_K), null);
        check(pm.errcode != 0, "incompatible options are an error");
        clear_model(&pm);
    }

    /* ---- the plugin table ---- */
    check(get_plugin_function("arma_model") !is null && get_plugin_function("arma_select") !is null && get_plugin_function("gui_arma_select") !is null,
          "arma_model, arma_select and gui_arma_select are found through the plugin table");

    destroy_dataset(dset);
    user_var_delete_by_name("iv_test", null);
    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
