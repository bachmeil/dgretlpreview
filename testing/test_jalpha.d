/* Tests for plugin/jalpha.d: restrictions on alpha, the adjustment coefficients of a VECM (Johansen 1995, section 8.2.1, "the same
   restriction on all alpha"), through the entry points that use it: vecm_test_restriction (the LR test) and johansen_estimate
   (estimation subject to the restriction). The restrictions here are R alpha = 0 with R = (1, 0) (x does not adjust) and R = (0, 1) (y does not adjust), on simulated
   data (x a random walk, y = 2 x + AR(1) noise).

   Everything is compared with an independent computation with the matrix operators of dgretl: the moment matrices conditional on
   the restriction (S00.a, S11.a, S01.a), the eigenvalues of S11.a^-1 S10.a Abar (Abar' S00.a Abar)^-1 Abar' S01.a, the LR
   statistic T sum ln((1 - lambda~)/(1 - lambda)) with r (p - s) degrees of freedom, alpha = S01.a beta (beta' S11.a beta)^-1.

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs, log;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;
import core.stdc.stdlib : calloc;
import core.stdc.string : strstr;

import std.math : PI;

import dgretl;

import gretl.dataset : create_new_dataset;
import gretl.gretl_list : gretl_list_new;
import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_get_buffer, gretl_print_new;
import gretl.gretl_prnh : PrnType;
import gretl.gretl_restricth : GretlRestriction;
import gretl.libgretlh : DatasetRaw, Gretlopt, TsCodes;
import gretl.plugin.jalpha : alt_get_eigenvalues, vecm_alpha_test;
import gretl.plugin.johansen : vecm_test_restriction;
import gretl.pvalues : chisq_cdf_comp;
import gretl.var : VarTask, gretl_VAR_free, gretl_VECM, johansen_wrapper;
import gretl.varh : GretlVar;

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

bool mclose(GretlMatrix a, GretlMatrix b, double tol = 1e-8)
{
    if (a.nrow != b.nrow || a.ncol != b.ncol) return false;
    foreach (j; 0 .. a.ncol) foreach (i; 0 .. a.nrow) if (!close(a[i, j], b[i, j], tol)) return false;
    return true;
}

GretlMatrix copyOf(const(GretlMatrixRaw) * p)
{
    return GretlMatrix(gretl_matrix_copy(p));
}

GretlMatrix resid(GretlMatrix Y, GretlMatrix X)
{
    return Y - X * solve(X, Y);
}

enum T = 120;
double[T] xs, ys;

GretlRestriction * alphaRestriction(double r0, double r1)
{
    GretlRestriction * rs = cast(GretlRestriction *) calloc(1, GretlRestriction.sizeof);

    rs.Ra = gretl_matrix_alloc(1, 2);
    rs.qa = gretl_matrix_alloc(1, 1);
    rs.Ra.val[0] = r0;
    rs.Ra.val[1] = r1;
    rs.qa.val[0] = 0;
    rs.g = 1;
    rs.gmax = 2;
    return rs;
}

/* the independent computation for case 3 (unrestricted constant), VAR of order 2; the restriction is Ra alpha = 0 with the
   given 1 x 2 Ra, and Abar for it (the nullspace A is (0, 1)' for Ra = (1, 0) and (1, 0)' for Ra = (0, 1)) */
struct Indep {
    GretlMatrix S00, S01, S11, S00a, S01a, S11a;
    double eig, eigR;
    int Teff;
}

Indep independent(double r0, double r1)
{
    Indep r;
    const t0 = 2, Te = T - t0;
    auto Z0 = GretlMatrix(gretl_matrix_alloc(Te, 2));
    auto Z1 = GretlMatrix(gretl_matrix_alloc(Te, 2));
    auto Z2 = GretlMatrix(gretl_matrix_alloc(Te, 3));

    foreach (s; 0 .. Te) {
        const t = t0 + s;
        Z0[s, 0] = xs[t] - xs[t - 1];
        Z0[s, 1] = ys[t] - ys[t - 1];
        Z1[s, 0] = xs[t - 1];
        Z1[s, 1] = ys[t - 1];
        Z2[s, 0] = xs[t - 1] - xs[t - 2];
        Z2[s, 1] = ys[t - 1] - ys[t - 2];
        Z2[s, 2] = 1;
    }
    auto R0 = resid(Z0, Z2), R1 = resid(Z1, Z2);
    r.S00 = t(R0) * R0 / cast(double) Te;
    r.S01 = t(R0) * R1 / cast(double) Te;
    r.S11 = t(R1) * R1 / cast(double) Te;
    r.Teff = Te;
    auto Ap = GretlMatrix(gretl_matrix_alloc(1, 2));
    Ap[0, 0] = r0;
    Ap[0, 1] = r1;
    auto C = t(Ap) * inv(Ap * r.S00 * t(Ap)) * Ap;
    r.S00a = r.S00 - r.S00 * C * r.S00;
    r.S11a = r.S11 - t(r.S01) * C * r.S01;
    r.S01a = r.S01 - r.S00 * C * r.S01;
    auto A = GretlMatrix(gretl_matrix_alloc(2, 1));
    A[0, 0] = r1 != 0 ? 1 : 0;
    A[1, 0] = r0 != 0 ? 1 : 0;                        // nullspace of (r0, r1) when one of them is zero
    auto Abar = A * inv(t(A) * A);
    auto Mr = inv(r.S11a) * t(r.S01a) * Abar * inv(t(Abar) * r.S00a * Abar) * t(Abar) * r.S01a;
    auto ev = eigenGeneral(Mr, false).values;
    r.eigR = ev[0, 0] > ev[1, 0] ? ev[0, 0] : ev[1, 0];
    auto M = inv(r.S11) * t(r.S01) * inv(r.S00) * r.S01;
    auto ev0 = eigenGeneral(M, false).values;
    r.eig = ev0[0, 0] > ev0[1, 0] ? ev0[0, 0] : ev0[1, 0];
    return r;
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    uint seed = 4242;
    double rnd()
    {
        seed = seed * 1103515245u + 12345u;
        return ((seed >> 8) & 0xffff) / 65536.0 - 0.5;
    }
    double u = 0;
    DatasetRaw * d = create_new_dataset(3, T, 0);

    d.pd = 1;
    d.structure = TsCodes.time_series;
    foreach (i, nm; ["const", "x", "y"]) {
        d.varname[i][0 .. nm.length] = nm[];
        d.varname[i][nm.length] = 0;
    }
    foreach (t; 0 .. T) {
        xs[t] = (t == 0 ? 0 : xs[t - 1]) + rnd();
        u = 0.5 * u + 0.4 * rnd();
        ys[t] = 2.0 * xs[t] + u;
        d.Z[0][t] = 1;
        d.Z[1][t] = xs[t];
        d.Z[2][t] = ys[t];
    }
    int * list = gretl_list_new(2);

    list[1] = 1;
    list[2] = 2;
    int err = 0;

    foreach (which; 0 .. 2) {
        const r0 = which == 0 ? 1.0 : 0.0, r1 = which == 0 ? 0.0 : 1.0;
        const label = which == 0 ? "alpha_x = 0" : "alpha_y = 0";
        auto ind = independent(r0, r1);
        const LR = ind.Teff * log((1 - ind.eigR) / (1 - ind.eig));

        GretlVar * v = gretl_VECM(2, 1, list, d, Gretlopt.OPT_Q, null, &err);

        check(v !is null && err == 0, "gretl_VECM (rank 1) for the restriction test: " ~ label);
        auto jr = v.jinfo;
        auto beta0 = copyOf(jr.Beta), alpha0 = copyOf(jr.Alpha);

        /* the LR test */
        GretlRestriction * rs = alphaRestriction(r0, r1);
        Prn * prn = gretl_print_new(PrnType.buffer, &err);

        err = vecm_test_restriction(v, rs, d, Gretlopt.OPT_NONE, prn);
        check(err == 0, "vecm_test_restriction runs for a restriction on alpha: " ~ label);
        check(close(rs.test, LR, 1e-6) && close(rs.pval, chisq_cdf_comp(1, LR), 1e-8), "LR statistic and p-value (1 df): " ~ label);
        char[64] lrbuf;
        import core.stdc.stdio : snprintf;
        snprintf(lrbuf.ptr, lrbuf.length, "%g", rs.test);
        check(strstr(gretl_print_get_buffer(prn), "Restricted loglikelihood") !is null && strstr(gretl_print_get_buffer(prn), lrbuf.ptr) !is null, "the output shows the restricted log-likelihood and the statistic: " ~ label);
        check(mclose(copyOf(jr.Beta), beta0) && mclose(copyOf(jr.Alpha), alpha0), "beta and alpha are restored after the test: " ~ label);
        gretl_print_destroy(prn);

        gretl_VAR_free(v);

        /* estimation subject to the restriction, from scratch as var.d does it */
        prn = gretl_print_new(PrnType.buffer, &err);
        v = johansen_wrapper(VarTask.vecm_estimate, 2, 1, null, list, rs, d, Gretlopt.OPT_NONE, prn, &err);
        check(v !is null && err == 0 && v.err == 0, "a VECM estimated subject to a restriction on alpha (johansen_wrapper): " ~ label);
        jr = v.jinfo;
        auto Ar = copyOf(jr.Alpha), Br = copyOf(jr.Beta);
        const fixedIdx = which == 0 ? 0 : 1;
        check(fabs(Ar[fixedIdx, 0]) < 1e-9 && fabs(Ar[1 - fixedIdx, 0]) > 1e-4, "the restricted element of alpha is zero, the other is not: " ~ label);
        auto alphaExp = ind.S01a * Br * inv(t(Br) * ind.S11a * Br);
        check(mclose(Ar, alphaExp, 1e-7), "alpha = S01.a beta (beta' S11.a beta)^-1: " ~ label);
        const llr = -(ind.Teff / 2.0) * (2 * (1 + log(2 * PI)) + log(det(ind.S00)) + log(1 - ind.eigR));
        check(close(v.ll, llr, 1e-9), "restricted log-likelihood: " ~ label);
        check(jr.lrdf == 1 && jr.Ra !is null && jr.qa !is null, "the df and the restriction matrices are recorded: " ~ label);
        gretl_print_destroy(prn);
        gretl_VAR_free(v);
    }

    /* the eigenvalue helper alone: with a trivial restriction (none) it gives the unrestricted eigenvalues */
    auto ind = independent(1.0, 0.0);
    auto AS00 = inv(ind.S00);
    auto Mtmp = GretlMatrix(gretl_matrix_alloc(2, 2));
    auto Mv = GretlMatrix(gretl_matrix_alloc(2, 2));
    GretlMatrixRaw * evals = null;
    auto tmp = GretlMatrix(gretl_matrix_alloc(2, 2));
    err = alt_get_eigenvalues(AS00.ptr, ind.S01.ptr, ind.S11.ptr, Mv.ptr, &evals, tmp.ptr, 1);
    check(err == 0 && evals !is null && close(evals.val[0], ind.eig, 1e-8), "alt_get_eigenvalues with S00^-1 gives the unrestricted eigenvalue");
    gretl_matrix_free(evals);

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
