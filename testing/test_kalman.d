/* Tests for kalman.d: the Kalman filter and smoother through the "Kalman bundle" interface of the library (kalman_bundle_new,
   maybe_set_kalman_element, kalman_bundle_filter, kalman_bundle_smooth, kalman_bundle_simulate, maybe_retrieve_kalman_element,
   kalman_bundle_copy, print_kalman_bundle_info, kalman_serialize) and the older direct interface used for ARMA estimation
   (kalman_new, kfilter_standard, kalman_get_loglik).

   Everything is checked against an independent computation written here with the matrix operators of dgretl: the filter
   v_t = y_t - Z a_t, F_t = Z P_t Z' + VY, K_t = T P_t Z' inv(F_t), a_{t+1} = T a_t + K_t v_t, P_{t+1} = T P_t T' + VS - K_t F_t K_t'
   with the Gaussian log-likelihood; the fixed-interval (RTS) state smoother; the disturbance smoother (Durbin-Koopman);
   simulation; and the exact diffuse initialization, for which the log-likelihood of the local level model is known in closed
   form (the first observation only fixes the state).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs, log;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;
import core.stdc.string : strstr;

import std.math : PI;

import dgretl;

import gretl.gretl_bundle : gretl_bundle_destroy, gretl_bundle_get_private_data, kalman_bundle_new;
import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_get_buffer, gretl_print_new;
import gretl.gretl_prnh : PrnType;
import gretl.kalman : kalman_bundle_copy, kalman_bundle_filter, kalman_bundle_n_members, kalman_bundle_simdata, kalman_bundle_simulate,
                      kalman_bundle_smooth, kalman_serialize, maybe_retrieve_kalman_element, maybe_set_kalman_element, print_kalman_bundle_info;
import gretl.libgretlh : GretlBundle, GretlType;

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

bool close(double a, double b, double tol = 1e-9)
{
    return fabs(a - b) <= tol * (1 + fabs(b));
}

bool mclose(GretlMatrix a, GretlMatrix b, double tol = 1e-9)
{
    if (a.nrow != b.nrow || a.ncol != b.ncol) return false;
    foreach (j; 0 .. a.ncol) foreach (i; 0 .. a.nrow) if (!close(a[i, j], b[i, j], tol)) return false;
    return true;
}

GretlMatrix mat(int r, int c, double[] v)
{
    auto m = GretlMatrix(gretl_matrix_alloc(r, c));

    foreach (i, x; v) m.ptr.val[i] = x;          // column by column
    return m;
}

GretlMatrix rowsOf(double[][] a)
{
    auto m = GretlMatrix(gretl_matrix_alloc(cast(int) a.length, cast(int) a[0].length));

    foreach (i, row; a) foreach (j, x; row) m[i, j] = x;
    return m;
}

ulong lcg = 99;

double urand()
{
    lcg = lcg * 6364136223846793005UL + 1442695040888963407UL;
    return ((lcg >> 11) / cast(double) (1UL << 53)) - 0.5;
}

/* a Kalman bundle and what is needed to set and read its members */
struct KB {
    GretlBundle * b;
    void * k;

    bool setMatrix(string key, GretlMatrix m)
    {
        int err = 0;

        maybe_set_kalman_element(k, key.ptr, m.ptr, GretlType.matrix, 1, &err);
        return err == 0;
    }

    bool setScalar(string key, double x)
    {
        int err = 0;

        maybe_set_kalman_element(k, key.ptr, &x, GretlType.double_, 1, &err);
        return err == 0;
    }

    GretlMatrix get(string key)
    {
        GretlType ty;
        int priv, err = 0;
        auto p = cast(GretlMatrixRaw *) maybe_retrieve_kalman_element(k, key.ptr, &ty, &priv, &err);

        return (p is null) ? GretlMatrix.init : GretlMatrix(gretl_matrix_copy(p));
    }

    double scalar(string key)
    {
        GretlType ty;
        int priv, err = 0;
        auto p = cast(double *) maybe_retrieve_kalman_element(k, key.ptr, &ty, &priv, &err);

        return (p is null) ? double.nan : *p;
    }

    int filter()
    {
        int err = 0;

        kalman_bundle_filter(b, null, &err);
        return err;
    }
}

/* ysT: y (N x n), ZT: r x n (Z'), T: r x r, VS: r x r; the observation variance and initial values are set separately */
KB makeBundle(GretlMatrix y, GretlMatrix ZT, GretlMatrix Tm, GretlMatrix VS, int * err)
{
    GretlMatrixRaw *[4] M = [y.ptr, ZT.ptr, Tm.ptr, VS.ptr];
    int[4] copy = [1, 1, 1, 1];
    KB kb;

    kb.b = kalman_bundle_new(M.ptr, copy.ptr, 4, 0, err);
    if (kb.b !is null) kb.k = gretl_bundle_get_private_data(kb.b);
    return kb;
}

/* the independent filter */
struct Filt {
    GretlMatrix V, F, A, P, K;           // row t holds v_t', vech(F_t), a_t', vech(P_t), vec(K_t)
    double ll;
    GretlMatrix aNext, PNext;
}

GretlMatrix row(GretlMatrix a, size_t t, size_t c)
{
    auto r = GretlMatrix(gretl_matrix_alloc(cast(int) c, 1));

    foreach (j; 0 .. c) r[j, 0] = a[t, j];
    return r;
}

/* vech, column by column of the lower triangle */
double[] vech(GretlMatrix m)
{
    double[] r;

    foreach (j; 0 .. m.ncol) foreach (i; j .. m.nrow) r ~= m[i, j];
    return r;
}

Filt indepFilter(GretlMatrix y, GretlMatrix ZT, GretlMatrix Tm, GretlMatrix VS, GretlMatrix VY, GretlMatrix a0, GretlMatrix P0, size_t first = 0)
{
    const N = y.nrow, n = y.ncol, r = Tm.nrow;
    Filt f;
    f.ll = 0;

    f.V = GretlMatrix(gretl_matrix_alloc(cast(int) N, cast(int) n));
    f.F = GretlMatrix(gretl_matrix_alloc(cast(int) N, cast(int) (n * (n + 1) / 2)));
    f.A = GretlMatrix(gretl_matrix_alloc(cast(int) N, cast(int) r));
    f.P = GretlMatrix(gretl_matrix_alloc(cast(int) N, cast(int) (r * (r + 1) / 2)));
    f.K = GretlMatrix(gretl_matrix_alloc(cast(int) N, cast(int) (r * n)));
    auto a = a0.dup, P = P0.dup;

    foreach (tt; 0 .. N) {
        auto v = row(y, tt, n) - t(ZT) * a;
        auto F = t(ZT) * P * ZT + VY;
        auto Kt = Tm * P * ZT * inv(F);

        foreach (j; 0 .. n) f.V[tt, j] = v[j, 0];
        foreach (j, x; vech(F)) f.F[tt, j] = x;
        foreach (j; 0 .. r) f.A[tt, j] = a[j, 0];
        foreach (j, x; vech(P)) f.P[tt, j] = x;
        size_t q = 0;
        foreach (j; 0 .. n) foreach (i; 0 .. r) f.K[tt, q++] = Kt[i, j];
        if (tt >= first) {
            f.ll += -0.5 * (n * log(2 * PI) + log(det(F)) + (t(v) * inv(F) * v)[0, 0]);
        }
        a = Tm * a + Kt * v;
        P = Tm * P * t(Tm) + VS - Kt * F * t(Kt);
    }
    f.aNext = a;
    f.PNext = P;
    return f;
}

/* model 1: local level, y_t = mu_t + e_t, mu_{t+1} = mu_t + eta_t */
struct Data1 {
    GretlMatrix y;
    int N;
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);
    int err = 0;

    /* ---- local level ---- */
    enum N = 40;
    double[] ys;
    double mu = 0;

    foreach (t; 0 .. N) { mu += 0.6 * urand(); ys ~= mu + urand(); }
    auto y = mat(N, 1, ys);
    auto kb = makeBundle(y, mat(1, 1, [1.0]), mat(1, 1, [1.0]), mat(1, 1, [0.4]), &err);

    check(kb.b !is null && err == 0, "kalman_bundle_new: a local level model");
    check(kb.setMatrix("obsvar", mat(1, 1, [1.0])) && kb.setMatrix("inistate", mat(1, 1, [0.0])) && kb.setMatrix("inivar", mat(1, 1, [10.0])), "setting obsvar, inistate and inivar");
    check(kb.filter() == 0, "kalman_bundle_filter runs");
    auto f1 = indepFilter(y, mat(1, 1, [1.0]), mat(1, 1, [1.0]), mat(1, 1, [0.4]), mat(1, 1, [1.0]), mat(1, 1, [0.0]), mat(1, 1, [10.0]));
    check(mclose(kb.get("prederr"), f1.V) && mclose(kb.get("pevar"), f1.F) && mclose(kb.get("state"), f1.A) && mclose(kb.get("stvar"), f1.P), "prediction errors, their variances, predicted states and their variances");
    check(mclose(kb.get("gain"), f1.K) && close(kb.scalar("lnl"), f1.ll, 1e-10), "gain matrices and the log-likelihood");
    double s2 = 0;
    foreach (t; 0 .. N) s2 += f1.V[t, 0] * f1.V[t, 0] / f1.F[t, 0];
    check(close(kb.scalar("s2"), s2 / N, 1e-9), "s2 is the mean of the standardized squared prediction errors");
    check(kb.scalar("r") == 1 && kb.scalar("n") == 1 && kb.scalar("N") == N && kb.scalar("diffuse") == 0, "scalar members: r, n, N, diffuse");

    /* ---- fixed-interval smoothing of the states ---- */
    {
        // RTS smoother: r_{t-1} = Z'F^-1 v + L' r_t, N_{t-1} = Z'F^-1 Z + L'N_t L, ahat_t = a_t + P_t r_{t-1}, V_t = P_t - P_t N_{t-1} P_t, L_t = T - K_t Z
        auto rr = mat(1, 1, [0.0]), NN = mat(1, 1, [0.0]);
        auto ahat = GretlMatrix(gretl_matrix_alloc(N, 1)), Vt = GretlMatrix(gretl_matrix_alloc(N, 1));
        foreach_reverse (t; 0 .. N) {
            const F = f1.F[t, 0], v = f1.V[t, 0], P = f1.P[t, 0], a = f1.A[t, 0], K = f1.K[t, 0];
            const L = 1.0 - K;

            const rprev = v / F + L * rr[0, 0];
            const Nprev = 1.0 / F + L * NN[0, 0] * L;

            ahat[t, 0] = a + P * rprev;
            Vt[t, 0] = P - P * Nprev * P;
            rr[0, 0] = rprev;
            NN[0, 0] = Nprev;
        }
        foreach (variant; ["standard", "sequential", "dejong"]) {
            auto kb2 = makeBundle(y, mat(1, 1, [1.0]), mat(1, 1, [1.0]), mat(1, 1, [0.4]), &err);

            kb2.setMatrix("obsvar", mat(1, 1, [1.0]));
            kb2.setMatrix("inistate", mat(1, 1, [0.0]));
            kb2.setMatrix("inivar", mat(1, 1, [10.0]));
            if (variant != "standard") kb2.setScalar(variant, 1);
            kb2.filter();
            kalman_bundle_smooth(kb2.b, 0, null);
            check(mclose(kb2.get("state"), ahat, 1e-8) && mclose(kb2.get("stvar"), Vt, 1e-8), "state smoothing (" ~ variant ~ "): smoothed states and their variances agree with the RTS recursions");
        }

        /* disturbance smoothing: eps_hat = VY (v/F - K r), eta_hat = VS r (r = r_t of the Durbin-Koopman recursions) */
        auto epsh = GretlMatrix(gretl_matrix_alloc(N, 1)), etah = GretlMatrix(gretl_matrix_alloc(N, 1));
        double rn = 0;

        foreach_reverse (t; 0 .. N) {
            const F = f1.F[t, 0], v = f1.V[t, 0], K = f1.K[t, 0];

            epsh[t, 0] = 1.0 * (v / F - K * rn);
            etah[t, 0] = 0.4 * rn;
            rn = v / F + (1.0 - K) * rn;
        }
        foreach (variant; ["standard", "sequential", "dejong"]) {
            auto kb3 = makeBundle(y, mat(1, 1, [1.0]), mat(1, 1, [1.0]), mat(1, 1, [0.4]), &err);

            kb3.setMatrix("obsvar", mat(1, 1, [1.0]));
            kb3.setMatrix("inistate", mat(1, 1, [0.0]));
            kb3.setMatrix("inivar", mat(1, 1, [10.0]));
            if (variant != "standard") kb3.setScalar(variant, 1);
            kb3.filter();
            kalman_bundle_smooth(kb3.b, 1, null);
            auto U = kb3.get("smdist");

            check(U.nrow == N && U.ncol == 2, "disturbance smoothing (" ~ variant ~ ") returns an N x (r + n) matrix");
            // the first row of the de Jong smoother is wrong in the C code (bugs.md #116: "if (t > 1)" should be "t > 0"), so start at row 1 there
            const t0 = (variant == "dejong") ? 1 : 0;
            bool dist = true;
            foreach (t; t0 .. N) dist = dist && close(U[t, 0], etah[t, 0], 1e-8) && close(U[t, 1], epsh[t, 0], 1e-8);
            check(dist, "smoothed state and observation disturbances (" ~ variant ~ ") agree with the Durbin-Koopman recursions");
        }
    }

    /* ---- the sequential and de Jong variants give the same (non-diffuse) result ---- */
    foreach (variant; ["sequential", "dejong"]) {
        auto kb4 = makeBundle(y, mat(1, 1, [1.0]), mat(1, 1, [1.0]), mat(1, 1, [0.4]), &err);

        kb4.setMatrix("obsvar", mat(1, 1, [1.0]));
        kb4.setMatrix("inistate", mat(1, 1, [0.0]));
        kb4.setMatrix("inivar", mat(1, 1, [10.0]));
        const okSet = kb4.setScalar(variant, 1);
        const rc = kb4.filter();

        check(okSet && rc == 0 && close(kb4.scalar("lnl"), f1.ll, 1e-9), "the " ~ variant ~ " filter gives the same log-likelihood");
        check(mclose(kb4.get("prederr"), f1.V, 1e-8) && mclose(kb4.get("pevar"), f1.F, 1e-8), "the " ~ variant ~ " filter gives the same prediction errors and variances");
    }

    /* ---- exact diffuse initialization of the local level model ---- */
    {
        // known result: the first observation fixes the state, a_2 = y_1, P_2 = VY + VS; the log-likelihood is that of observations 2 .. N
        const VYv = 1.0, VSv = 0.4;
        double a = ys[0], P = VYv + VSv, ll = 0;

        foreach (t; 1 .. N) {
            const F = P + VYv, v = ys[t] - a, K = P / F;

            ll += -0.5 * (log(2 * PI) + log(F) + v * v / F);
            a += K * v;
            P = P * (1 - K) + VSv;
        }
        foreach (variant; ["sequential", "dejong"]) {
            auto kb5 = makeBundle(y, mat(1, 1, [1.0]), mat(1, 1, [1.0]), mat(1, 1, [VSv]), &err);

            kb5.setMatrix("obsvar", mat(1, 1, [VYv]));
            kb5.setScalar("diffuse", 2);
            kb5.setScalar(variant, 1);
            const rc = kb5.filter();

            // the de Jong code starts from P*_1 = VS, which adds 0.5 log(F*_1) = 0.5 log(VY + VS) (bugs.md #117), as ported
            const expected = (variant == "dejong") ? ll + 0.5 * log(VYv + VSv) : ll;

            check(rc == 0 && close(kb5.scalar("lnl"), expected, 1e-9), "exact diffuse log-likelihood of the local level model (" ~ variant ~ ")");
            check(kb5.scalar("diffuse") == 2 && kb5.scalar("exact") == 1, "the initialization is flagged exact diffuse (" ~ variant ~ ")");
        }
    }

    /* ---- a two-state model with a vector of observables: local linear trend with two series sharing the state ---- */
    {
        enum N2 = 30;
        auto Z2T = rowsOf([[1.0, 0.5], [0.0, 1.0]]);                // r x n = 2 x 2: Z' of the two observation equations
        auto T2 = rowsOf([[1.0, 1.0], [0.0, 1.0]]);
        auto VS2 = rowsOf([[0.3, 0.05], [0.05, 0.2]]);
        auto VY2 = rowsOf([[1.0, 0.2], [0.2, 0.8]]);
        auto a02 = mat(2, 1, [0.5, -0.2]);
        auto P02 = rowsOf([[2.0, 0.3], [0.3, 1.5]]);
        auto y2 = GretlMatrix(gretl_matrix_alloc(N2, 2));
        double l = 0, s = 0;

        foreach (t; 0 .. N2) { l += 0.3 * urand(); s += 0.1 * urand() + 0.05; y2[t, 0] = l + 0.5 * s + urand(); y2[t, 1] = s + 0.8 * urand(); }
        auto kb6 = makeBundle(y2, Z2T, T2, VS2, &err);

        check(kb6.b !is null && err == 0, "kalman_bundle_new: two states, two observables");
        kb6.setMatrix("obsvar", VY2);
        kb6.setMatrix("inistate", a02);
        kb6.setMatrix("inivar", P02);
        check(kb6.filter() == 0, "the filter runs for n = 2, r = 2");
        auto f2 = indepFilter(y2, Z2T, T2, VS2, VY2, a02, P02);
        check(mclose(kb6.get("prederr"), f2.V, 1e-8) && mclose(kb6.get("pevar"), f2.F, 1e-8), "prediction errors and their (vech) variances, n = 2");
        check(mclose(kb6.get("state"), f2.A, 1e-8) && mclose(kb6.get("stvar"), f2.P, 1e-8), "predicted states and their (vech) variances, r = 2");
        check(mclose(kb6.get("gain"), f2.K, 1e-8) && close(kb6.scalar("lnl"), f2.ll, 1e-9), "gain matrices (vec of the r x n matrix) and the log-likelihood");

        /* the sequential filter (needs a diagonal VY: use independent observables) */
        auto VY3 = rowsOf([[1.0, 0.0], [0.0, 0.8]]);
        auto kb7 = makeBundle(y2, Z2T, T2, VS2, &err);

        kb7.setMatrix("obsvar", VY3);
        kb7.setMatrix("inistate", a02);
        kb7.setMatrix("inivar", P02);
        kb7.setScalar("sequential", 1);
        auto f3 = indepFilter(y2, Z2T, T2, VS2, VY3, a02, P02);
        check(kb7.filter() == 0 && close(kb7.scalar("lnl"), f3.ll, 1e-9), "the sequential filter, n = 2 with independent observables: log-likelihood");

        /* smoothing, two states: the standard, sequential and de Jong smoothers agree (states and disturbances) */
        {
            GretlMatrix[3] st, ds;

            foreach (i, variant; ["standard", "sequential", "dejong"]) {
                foreach (dist; 0 .. 2) {
                    auto kbs = makeBundle(y2, Z2T, T2, VS2, &err);

                    kbs.setMatrix("obsvar", VY3);
                    kbs.setMatrix("inistate", a02);
                    kbs.setMatrix("inivar", P02);
                    if (i > 0) kbs.setScalar(variant, 1);
                    kbs.filter();
                    kalman_bundle_smooth(kbs.b, dist, null);
                    if (dist) {
                        ds[i] = kbs.get("smdist");
                    } else {
                        st[i] = kbs.get("state");
                    }
                }
            }
            check(mclose(st[0], st[1], 1e-7) && mclose(st[0], st[2], 1e-7), "two states: smoothed states agree across the standard, sequential and de Jong smoothers");
            check(mclose(ds[0], ds[1], 1e-7), "two states: smoothed disturbances of the standard and sequential smoothers agree");
            check(ds[2].nrow == N2 && ds[2].ncol == 4, "two states: de Jong disturbance smoothing has r + n columns");
        }

        /* the cross-correlated disturbances form (H, G) and the Durbin-Koopman form (R, Q) of the same model */
        {
            // DK: state variance R Q R' with R = I and Q = VS gives the standard model
            GretlMatrixRaw *[5] M5 = [y2.ptr, Z2T.ptr, T2.ptr, VS2.ptr, mat(2, 2, [1.0, 0.0, 0.0, 1.0]).ptr];
            int[5] copy5 = [1, 1, 1, 1, 1];
            auto Ridentity = mat(2, 2, [1.0, 0.0, 0.0, 1.0]);

            M5[4] = Ridentity.ptr;
            auto bdk = kalman_bundle_new(M5.ptr, copy5.ptr, 5, 1, &err);
            KB kdk;

            check(bdk !is null && err == 0, "kalman_bundle_new with Q and R (the Durbin-Koopman form)");
            kdk.b = bdk;
            kdk.k = gretl_bundle_get_private_data(bdk);
            kdk.setMatrix("obsvar", VY2);
            kdk.setMatrix("inistate", a02);
            kdk.setMatrix("inivar", P02);
            check(kdk.filter() == 0 && close(kdk.scalar("lnl"), f2.ll, 1e-9), "the Q, R form with R = I gives the same log-likelihood");
            check(kdk.scalar("dkvar") == 1, "the dkvar member is 1");
        }

        /* simulation with given disturbances: y_t = Z a_t + eps_t, a_{t+1} = T a_t + eta_t  (a_1 = initial state) */
        {
            auto kb8 = makeBundle(y2, Z2T, T2, VS2, &err);

            kb8.setMatrix("obsvar", VY2);
            kb8.setMatrix("inistate", a02);
            kb8.setMatrix("inivar", P02);
            auto Ud = GretlMatrix(gretl_matrix_alloc(N2, 4));            // eta (2) then eps (2) per period
            foreach (t; 0 .. N2) foreach (j; 0 .. 4) Ud[t, j] = urand();
            int e3 = 0;
            auto sim = kalman_bundle_simulate(kb8.b, Ud.ptr, 0, null, &e3);

            check(sim !is null && e3 == 0 && sim.rows == N2 && sim.cols == 2, "kalman_bundle_simulate returns an N x n matrix of simulated observables");
            if (sim !is null && sim.rows == N2) {
                auto S = GretlMatrix(sim);
                auto a = a02.dup;
                bool ok = true;

                foreach (tt; 0 .. N2) {
                    // U holds the state disturbance (r columns) and then the observation disturbance (n columns), unscaled
                    auto eta = mat(2, 1, [Ud[tt, 0], Ud[tt, 1]]);
                    auto eps = mat(2, 1, [Ud[tt, 2], Ud[tt, 3]]);
                    auto yt = t(Z2T) * a + eps;

                    ok = ok && close(S[tt, 0], yt[0, 0], 1e-9) && close(S[tt, 1], yt[1, 0], 1e-9);
                    a = T2 * a + eta;
                }
            check(ok, "simulated observables follow the state-space recursions with the given disturbances");
            }
        }
    }

    /* ---- copies, member count, printing ---- */
    {
        int e = 0;
        auto cp = kalman_bundle_copy(kb.b, &e);

        check(cp !is null && e == 0, "kalman_bundle_copy");
        if (cp !is null) {
            KB c;

            c.b = cp;
            c.k = gretl_bundle_get_private_data(cp);
            check(c.filter() == 0 && close(c.scalar("lnl"), f1.ll, 1e-10) && mclose(c.get("state"), f1.A), "the copy filters to the same result");
            gretl_bundle_destroy(cp);
        }
        check(kalman_bundle_n_members(kb.b) > 10, "kalman_bundle_n_members counts the matrices and scalars");
        Prn * prn = gretl_print_new(PrnType.buffer, &err);

        print_kalman_bundle_info(kb.k, prn);
        check(strstr(gretl_print_get_buffer(prn), "obsy") !is null, "print_kalman_bundle_info lists the members");
        gretl_print_destroy(prn);
        prn = gretl_print_new(PrnType.buffer, &err);
        kalman_serialize(kb.k, prn);
        check(strstr(gretl_print_get_buffer(prn), "gretl-kalman") !is null, "kalman_serialize writes XML");
        gretl_print_destroy(prn);
    }

    /* ---- errors ---- */
    {
        GretlMatrixRaw *[4] M = [y.ptr, null, mat(1, 1, [1.0]).ptr, mat(1, 1, [0.4]).ptr];
        int[4] copy = [1, 1, 1, 1];
        int e = 0;
        auto bad = kalman_bundle_new(M.ptr, copy.ptr, 4, 0, &e);

        check(bad is null && e != 0, "a missing obsymat is an error");
        auto kbBad = makeBundle(y, mat(1, 1, [1.0]), mat(1, 1, [1.0]), mat(1, 1, [0.4]), &e);
                check(!kbBad.setMatrix("statemat", mat(2, 2, [1.0, 0.0, 0.0, 1.0])), "replacing the state matrix with one of the wrong size is refused");
    }

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
