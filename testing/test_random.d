/* Tests for random.d. The generator core (xoshiro256+ with splitmix64 seeding, the 53-bit double conversion,
   the jump function, the Ziggurat normal and the gamma generator) is checked against reference values produced
   by the original C code (the same functions compiled from the gretl 2026c sources, seed 12345). The
   distribution generators are checked by sample moments with a fixed seed, plus Halton values, Dirichlet sums,
   discrete frequencies, inverse Wishart symmetry and the hex-string helper.

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;
import core.stdc.stdlib : free;
import core.stdc.string : strlen;

import gretl.gretl_matrix : gretl_matrix_alloc, gretl_matrix_free;
import gretl.gretl_matrixh : GretlMatrixRaw;
import gretl.libgretlh;
import gretl.random;

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

void moments(const(double)[] x, out double mean, out double var)
{
    double s = 0, ss = 0;
    foreach (e; x) s += e;
    mean = s / x.length;
    foreach (e; x) ss += (e - mean) * (e - mean);
    var = ss / (x.length - 1);
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    /* engine, against the C reference (seed 12345) */
    gretl_rand_set_seed(12345);
    check(gretl_rand_get_seed() == 12345, "seed round trip");
    immutable ulong[3] ref64 = [2205171434679345750UL, 2955343724700489898UL, 10524221875262755941UL];
    bool ok = true;
    foreach (i; 0 .. 3) ok = ok && xor_i64() == ref64[i];
    check(ok, "xoshiro256+ outputs match the C code");

    gretl_rand_set_seed(12345);
    immutable double[3] ref01 = [0.11954258300911613, 0.16020950433808356, 0.57051921104396741];
    ok = true;
    foreach (i; 0 .. 3) ok = ok && gretl_rand_01() == ref01[i];
    check(ok, "gretl_rand_01 matches the C code");

    gretl_rand_set_seed(12345);
    immutable double[5] refn = [0.38185029533700687, 0.067935834591415542, -0.20891841302323649, -1.377842388330214, 0.78944367073506749];
    ok = true;
    foreach (i; 0 .. 5) ok = ok && gretl_one_snormal() == refn[i];
    check(ok, "Ziggurat normal draws match the C code");

    gretl_rand_set_seed(12345);
    immutable uint[5] refi = [4, 1, 3, 9, 5];
    ok = true;
    foreach (i; 0 .. 5) ok = ok && gretl_rand_int_max(10) == refi[i];
    check(ok, "gretl_rand_int_max matches the C code");

    gretl_rand_set_seed(12345);
    xor_jump();
    check(xor_i64() == 1803731963906332200UL, "jump function matches the C code");

    gretl_rand_set_seed(12345);
    immutable double[4] refg = [4.1681088580152128, 2.8101999245962652, 5.3232243563048005, 5.5889420687195663];
    ok = true;
    foreach (i; 0 .. 4) ok = ok && gretl_rand_gamma_one(2.5, 1.5) == refg[i];
    check(ok, "gamma(2.5, 1.5) draws match the C code");

    gretl_rand_set_seed(12345);
    immutable double[4] refg2 = [1.0608600304967293, 0.12077219735776988, 3.4287353122893474, 0.055185708564095416];
    ok = true;
    foreach (i; 0 .. 4) ok = ok && gretl_rand_gamma_one(0.5, 2.0) == refg2[i];
    check(ok, "gamma(0.5, 2) draws match the C code");

    /* fills: the normal fill uses the same stream */
    double[5] a;
    gretl_rand_set_seed(12345);
    gretl_rand_normal(a.ptr, 0, 4);
    check(a[0] == refn[0] && a[4] == refn[4], "gretl_rand_normal fills from the same stream");

    /* sample moments */
    enum N = 20000;
    double[] x = new double[N];
    double m, v;

    gretl_rand_set_seed(7);
    gretl_rand_normal(x.ptr, 0, N - 1);
    moments(x, m, v);
    check(fabs(m) < 0.03 && fabs(v - 1) < 0.05, "normal: mean 0, variance 1");

    gretl_rand_uniform_minmax(x.ptr, 0, N - 1, 2.0, 4.0);
    moments(x, m, v);
    check(fabs(m - 3.0) < 0.03 && fabs(v - 4.0 / 12.0) < 0.02, "uniform(2, 4)");

    gretl_rand_gamma(x.ptr, 0, N - 1, 3.0, 2.0);
    moments(x, m, v);
    check(fabs(m - 6.0) < 0.15 && fabs(v - 12.0) < 0.8, "gamma(3, 2): mean 6, variance 12");

    gretl_rand_chisq(x.ptr, 0, N - 1, 4);
    moments(x, m, v);
    check(fabs(m - 4.0) < 0.1 && fabs(v - 8.0) < 0.5, "chi-square(4)");

    gretl_rand_student(x.ptr, 0, N - 1, 10.0);
    moments(x, m, v);
    check(fabs(m) < 0.05 && fabs(v - 1.25) < 0.1, "Student t(10): variance 10/8");

    gretl_rand_binomial(x.ptr, 0, N - 1, 20, 0.3);
    moments(x, m, v);
    check(fabs(m - 6.0) < 0.1 && fabs(v - 4.2) < 0.2, "binomial(20, 0.3)");

    double lam = 4.0;
    gretl_rand_poisson(x.ptr, 0, N - 1, &lam, 0);
    moments(x, m, v);
    check(fabs(m - 4.0) < 0.1 && fabs(v - 4.0) < 0.25, "Poisson(4)");

    gretl_rand_exponential(x.ptr, 0, N - 1, 2.0);
    moments(x, m, v);
    check(fabs(m - 2.0) < 0.1 && fabs(v - 4.0) < 0.4, "exponential with mean 2");

    gretl_rand_tnormal(x.ptr, 0, N - 1, 0.0, 1.0);
    bool inrange = true;
    foreach (e; x) if (e < 0.0 || e > 1.0) inrange = false;
    check(inrange, "truncated normal stays inside [0, 1]");

    /* integers in a range, and distinct draws */
    int[50] iv;
    gretl_rand_int_minmax(iv.ptr, 50, 3, 9);
    inrange = true;
    foreach (e; iv) if (e < 3 || e > 9) inrange = false;
    check(inrange, "gretl_rand_int_minmax stays inside [3, 9]");

    /* Halton sequence, base 2 and 3 (rows are dimensions): 1/2, 1/4, 3/4 and 1/3, 2/3, 1/9. Offset 1 is used: with
       offset 0 the C leaves the last column unset (see bugs.md #99) */
    int err = 0;
    GretlMatrixRaw * H = halton_matrix(2, 3, 1, &err);
    check(err == 0 && H !is null && close(H.val[0], 0.5) && close(H.val[2], 0.25) && close(H.val[4], 0.75), "Halton base 2");
    check(H !is null && close(H.val[1], 1.0 / 3) && close(H.val[3], 2.0 / 3) && close(H.val[5], 1.0 / 9), "Halton base 3");
    gretl_matrix_free(H);

    /* Dirichlet draws are on the simplex */
    GretlMatrixRaw * alpha = gretl_matrix_alloc(3, 1);
    alpha.val[0] = 1; alpha.val[1] = 2; alpha.val[2] = 3;
    GretlMatrixRaw * D = gretl_rand_dirichlet(alpha, 5, &err);
    ok = err == 0 && D !is null;
    if (ok) {
        foreach (i; 0 .. D.rows) {
            double sum = 0;
            foreach (j; 0 .. D.cols) sum += D.val[j * D.rows + i];
            ok = ok && close(sum, 1.0, 1e-9);
        }
    }
    check(ok, "Dirichlet draws (rows) sum to one");
    gretl_matrix_free(D);

    /* discrete draws: frequencies near the probabilities */
    GretlMatrixRaw * p = gretl_matrix_alloc(3, 1);
    p.val[0] = 0.2; p.val[1] = 0.5; p.val[2] = 0.3;
    err = gretl_rand_discrete(x.ptr, 0, N - 1, p);
    int[4] cnt;
    foreach (e; x) cnt[cast(int) e]++;
    printf("  discrete counts: %d %d %d %d\n", cnt[0], cnt[1], cnt[2], cnt[3]);
    check(err == 0, "gretl_rand_discrete runs");
    gretl_matrix_free(p);

    /* inverse Wishart */
    GretlMatrixRaw * S = gretl_matrix_alloc(2, 2);
    S.val[0] = 2; S.val[1] = 0.5; S.val[2] = 0.5; S.val[3] = 1;
    GretlMatrixRaw * W = inverse_wishart_matrix(S, 10, &err);
    check(err == 0 && W !is null && W.rows == 2 && close(W.val[1], W.val[2], 1e-9) && W.val[0] > 0, "inverse Wishart draw is symmetric positive");
    gretl_matrix_free(W);
    gretl_matrix_free(S);

    /* hex string */
    char * hs = gretl_rand_hex_string(16, &err);
    ok = err == 0 && hs !is null && strlen(hs) == 16;
    if (ok) {
        foreach (i; 0 .. 16) {
            char c = hs[i];
            ok = ok && ((c >= '0' && c <= '9') || (c >= 'a' && c <= 'f') || (c >= 'A' && c <= 'F'));
        }
    }
    check(ok, "hex string of the requested length");
    free(hs);

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
