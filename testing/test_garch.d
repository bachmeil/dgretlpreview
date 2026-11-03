/* Tests for plugin/garch.d (garch.c): the GARCH driver garch_model(): list checking, the initial OLS, the scaling of y and its undoing,
   the model statistics (parameter names, sigma, residuals, the LR test against the OLS fit, standardized residuals), the ARMA-based
   initialization of the variance parameters, verbose output, and the registration in the plugin table.

   The estimators themselves are tested in test_fcp.d (fcp.c, option --fcp) and test_mod_garch.d (mod_garch.c, the default).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs, log, sqrt;
import core.stdc.stdio : printf, setvbuf, stdout, _IOLBF;
import core.stdc.string : strcmp, strcpy, strstr;

import std.math : PI;

import dgretl;

import gretl.dataset : create_new_dataset, dataset_set_time_series, destroy_dataset;
import gretl.gretl_list : gretl_list_new;
import gretl.gretl_listh : LISTSEP;
import gretl.gretl_model : clear_model, gretl_model_get_data, gretl_model_get_double, gretl_model_get_int;
import gretl.gretl_prn : gretl_print_destroy, gretl_print_get_buffer, gretl_print_new;
import gretl.gretl_prnh : PrnType;
import gretl.libgretlh : DatasetRaw, Gretlopt, Model;
import gretl.plugin.garch : garch_model;
import gretl.plugins : get_plugin_function;

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

ulong lcg = 777;

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

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    enum T = 500;
    auto dset = create_new_dataset(3, T, 0);
    double[] ys = new double[T], xs = new double[T];
    double e = 0, h = 1;

    foreach (t; 0 .. T) {
        h = 0.5 + 0.2 * e * e + 0.6 * h;
        e = sqrt(h) * urand();
        xs[t] = urand();
        ys[t] = 4.0 + 2.0 * xs[t] + e;
        dset.Z[1][t] = ys[t];
        dset.Z[2][t] = xs[t];
    }
    dataset_set_time_series(dset, 1, 1900, 1);
    strcpy(dset.varname[1], "y");
    strcpy(dset.varname[2], "x");

    /* ---- list checking ---- */
    {
        auto L1 = garchList(1, 1, [1]);
        L1[3] = 2;                                          // the separator is missing
        auto m = garch_model(L1, dset, null, Gretlopt.OPT_NONE);

        check(m.errcode != 0, "a list without the separator is an error");
        clear_model(&m);
        auto L2 = garchList(-1, 1, [1, 2]);
        m = garch_model(L2, dset, null, Gretlopt.OPT_NONE);
        check(m.errcode != 0, "a negative order is an error");
        clear_model(&m);
        auto L3 = garchList(1, 0, [1, 2]);
        m = garch_model(L3, dset, null, Gretlopt.OPT_NONE);
        check(m.errcode != 0, "p > 0 with q = 0 (unidentified) is an error");
        clear_model(&m);
        auto L4 = garchList(4, 4, [1, 2]);
        m = garch_model(L4, dset, null, Gretlopt.OPT_NONE);
        check(m.errcode != 0, "p + q above 7 is an error");
        clear_model(&m);
    }

    /* ---- a GARCH(1,1) with y on a constant and x ---- */
    {
        auto L = garchList(1, 1, [1, 2]);                    // garch 1 1 ; y x  (the constant is added)
        auto m = garch_model(L, dset, null, Gretlopt.OPT_NONE);

        check(m.errcode == 0, "garch_model: no error");
        if (m.errcode == 0) {
            bool restored = true;

            foreach (t; 0 .. T) restored = restored && close(dset.Z[1][t], ys[t], 1e-12);
            check(restored, "the dependent variable of the dataset is restored after estimation");
            check(m.ncoeff == 5 && m.nobs == T && m.t1 == 0 && m.t2 == T - 1, "five coefficients (const, x, alpha0, alpha1, beta1) on the whole sample");
            check(strcmp(m.params[0], "const") == 0 && strcmp(m.params[1], "x") == 0 && strcmp(m.params[2], "alpha(0)") == 0
                  && strcmp(m.params[3], "alpha(1)") == 0 && strcmp(m.params[4], "beta(1)") == 0, "parameter names");
            check(gretl_model_get_int(&m, "garch_p") == 1 && gretl_model_get_int(&m, "garch_q") == 1, "garch_p and garch_q recorded");
            // sigma is the unconditional standard deviation, residuals are in the units of the data
            const den = 1 - m.coeff[3] - m.coeff[4];

            check(close(m.sigma * m.sigma, m.coeff[2] / den, 1e-9), "sigma^2 = alpha0 / (1 - alpha1 - beta1)");
            double ess = 0, ssr = 0;

            foreach (t; 0 .. T) {
                ess += m.uhat[t] * m.uhat[t];
                const r = ys[t] - m.coeff[0] - m.coeff[1] * xs[t];

                ssr += (r - m.uhat[t]) * (r - m.uhat[t]);
            }
            check(close(m.ess, ess, 1e-9) && ssr < 1e-12 * T, "ess is the sum of squared residuals; the residuals are y - Xb in the original units");
            // the LR test of the GARCH terms against the OLS fit
            auto X = GretlMatrix(gretl_matrix_alloc(T, 2));
            auto yv = GretlMatrix(gretl_matrix_alloc(T, 1));

            foreach (t; 0 .. T) {
                X[t, 0] = 1;
                X[t, 1] = xs[t];
                yv[t, 0] = ys[t];
            }
            auto beta = inv(t(X) * X) * (t(X) * yv);
            auto u = yv - X * beta;
            const llols = -0.5 * T * (log(2 * PI) + 1 + log((t(u) * u)[0, 0] / T));
            const LR = 2 * (m.lnL - llols);

            check(close(gretl_model_get_double(&m, "garch_LR"), LR, 1e-8) && gretl_model_get_int(&m, "garch_LR_df") == 2, "the LR statistic against OLS and its 2 degrees of freedom");
            check(gretl_model_get_data(&m, "garch_h") !is null, "the conditional variance series is attached");

            // standardized residuals (--stdresid = OPT_Z)
            auto Lz = garchList(1, 1, [1, 2]);
            auto mz = garch_model(Lz, dset, null, Gretlopt.OPT_Z);

            check(mz.errcode == 0, "garch_model --stdresid: no error");
            if (mz.errcode == 0) {
                const hz = cast(double *) gretl_model_get_data(&mz, "garch_h");
                bool ok = hz !is null;

                if (ok) foreach (t; 0 .. T) ok = ok && close(mz.uhat[t], m.uhat[t] / sqrt(hz[t]), 1e-6);
                check(ok, "the residuals are the unstandardized ones divided by the conditional standard deviation");
            }
            clear_model(&mz);
        }
        clear_model(&m);
    }

    /* ---- initialization of the variance parameters from an ARMA fit of the squared residuals (--arma-init = OPT_A) ---- */
    {
        auto L = garchList(1, 1, [1, 2]);
        auto m1 = garch_model(L, dset, null, Gretlopt.OPT_NONE);
        auto m2 = garch_model(L, dset, null, Gretlopt.OPT_A);

        check(m1.errcode == 0 && m2.errcode == 0, "garch_model --arma-init: no error");
        if (m1.errcode == 0 && m2.errcode == 0) {
            bool same = close(m1.lnL, m2.lnL, 1e-6);

            foreach (i; 0 .. 5) same = same && fabs(m1.coeff[i] - m2.coeff[i]) < 5e-3 * (1 + fabs(m1.coeff[i]));
            check(same, "the ARMA-based initialization reaches the same maximum");
        }
        clear_model(&m1);
        clear_model(&m2);
    }

    /* ---- verbose output, and the constant given explicitly or suppressed ---- */
    {
        int perr = 0;
        auto prn = gretl_print_new(PrnType.buffer, &perr);
        auto L = garchList(1, 1, [1, 2]);
        auto m = garch_model(L, dset, prn, Gretlopt.OPT_V);

        check(m.errcode == 0 && strstr(gretl_print_get_buffer(prn), "Iteration") !is null, "--verbose prints the iterations");
        clear_model(&m);
        gretl_print_destroy(prn);

        auto L2 = garchList(1, 1, [1, 0, 2]);                // constant given explicitly, between y and x
        auto m2 = garch_model(L2, dset, null, Gretlopt.OPT_NONE);

        check(m2.errcode == 0 && m2.ncoeff == 5 && strcmp(m2.params[0], "const") == 0, "an explicit constant is the first coefficient");
        clear_model(&m2);
        auto d0 = create_new_dataset(3, T, 0);                 // the same data with the intercept removed
        double ymean = 4.0;

        foreach (t; 0 .. T) { d0.Z[1][t] = ys[t] - ymean; d0.Z[2][t] = xs[t]; }
        dataset_set_time_series(d0, 1, 1900, 1);
        strcpy(d0.varname[1], "y");
        strcpy(d0.varname[2], "x");
        auto L3 = garchList(1, 1, [1, 2]);
        auto m3 = garch_model(L3, d0, null, Gretlopt.OPT_N);

        check(m3.errcode == 0 && m3.ncoeff == 4 && strcmp(m3.params[0], "x") == 0, "--nc: no constant is added");
        clear_model(&m3);
        destroy_dataset(d0);
    }

    check(get_plugin_function("garch_model") !is null, "garch_model is found through the plugin table");

    destroy_dataset(dset);
    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
