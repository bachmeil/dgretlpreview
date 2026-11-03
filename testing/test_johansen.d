/* Tests for plugin/johansen.d: the Johansen cointegration test and VECM estimation through the library entry points johansen_test,
   gretl_VECM (-> johansen_coint_test, johansen_estimate), vecm_test_restriction and the bootstrap impulse responses
   (johansen_boot_round).

   Everything is checked against an independent computation written here with the matrix operators of dgretl: the residual
   matrices R0 and R1 of the two auxiliary regressions, S_ij = R_i' R_j / T, the eigenvalues of
   inv(S11) S10 inv(S00) S01, the trace statistics, the log-likelihood, alpha = S01 beta inv(beta' S11 beta), the standard errors of
   Phillips-normalized beta, and the LR test of a restriction on beta through the restricted eigenproblem, for the five
   deterministic cases (no constant, restricted constant, unrestricted constant, restricted trend, unrestricted trend).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs, log, sqrt;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;
import core.stdc.stdlib : calloc;
import core.stdc.string : strstr;

import std.conv : to;
import std.format : format;
import std.math : PI;

import dgretl;

import gretl.dataset : create_new_dataset;
import gretl.gretl_list : gretl_list_new;
import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_get_buffer, gretl_print_new;
import gretl.gretl_prnh : PrnType;
import gretl.gretl_restricth : GretlRestriction;
import gretl.libgretlh : DatasetRaw, Gretlopt, TsCodes;
import gretl.libset : execute_set;
import gretl.plugin.johansen : johansen_boot_round, trace_pvalue, vecm_test_restriction;
import gretl.pvalues : chisq_cdf_comp;
import gretl.var : VarTask, gretl_VAR_free, gretl_VAR_get_impulse_response, gretl_VECM, johansen_test, johansen_wrapper;
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

/* the residuals of the columns of Y regressed on the columns of X */
GretlMatrix resid(GretlMatrix Y, GretlMatrix X)
{
    return Y - X * solve(X, Y);
}

enum T = 120;
double[T] xs, ys;

enum Case { noConst, restConst, unrestConst, restTrend, unrestTrend }

/* the independent computation, for a VAR of order p = 2 in levels (one lagged difference), n = 2 */
struct Indep {
    GretlMatrix R0, R1, S00, S01, S11;
    double[] eig;            // descending
    int Teff;
}

Indep independent(Case c)
{
    Indep r;
    const t0 = 2, Te = T - t0;
    auto Z0 = GretlMatrix(gretl_matrix_alloc(Te, 2));
    auto Z1base = GretlMatrix(gretl_matrix_alloc(Te, 2));
    auto D1 = GretlMatrix(gretl_matrix_alloc(Te, 2));
    foreach (s; 0 .. Te) {
        const t = t0 + s;
        Z0[s, 0] = xs[t] - xs[t - 1];
        Z0[s, 1] = ys[t] - ys[t - 1];
        Z1base[s, 0] = xs[t - 1];
        Z1base[s, 1] = ys[t - 1];
        D1[s, 0] = xs[t - 1] - xs[t - 2];
        D1[s, 1] = ys[t - 1] - ys[t - 2];
    }
    auto one = new double[Te], trend = new double[Te];
    foreach (s; 0 .. Te) { one[s] = 1; trend[s] = t0 + s; }
    GretlMatrix Z1, Z2;
    final switch (c) {
        case Case.noConst:
            Z1 = Z1base; Z2 = D1; break;
        case Case.restConst:
            Z1 = cbind(Z1base, one); Z2 = D1; break;
        case Case.unrestConst:
            Z1 = Z1base; Z2 = cbind(D1, one); break;
        case Case.restTrend:
            Z1 = cbind(Z1base, trend); Z2 = cbind(D1, one); break;
        case Case.unrestTrend:
            Z1 = Z1base; Z2 = cbind(D1, one, trend); break;
    }
    r.R0 = resid(Z0, Z2);
    r.R1 = resid(Z1, Z2);
    r.S00 = t(r.R0) * r.R0 / cast(double) Te;
    r.S01 = t(r.R0) * r.R1 / cast(double) Te;
    r.S11 = t(r.R1) * r.R1 / cast(double) Te;
    auto M = inv(r.S11) * t(r.S01) * inv(r.S00) * r.S01;
    auto ev = eigenGeneral(M, false).values;
    foreach (i; 0 .. ev.nrow) r.eig ~= ev[i, 0];
    r.eig.sort!"a > b";
    r.Teff = Te;
    return r;
}

import std.algorithm.sorting : sort;

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    /* a cointegrated pair: x a random walk, y = 2 x + AR(1) noise */
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

    /* ---- the cointegration test in the five cases ---- */
    static immutable Gretlopt[5] opts = [Gretlopt.OPT_N, Gretlopt.OPT_R, Gretlopt.OPT_NONE, Gretlopt.OPT_A, Gretlopt.OPT_T];
    static immutable string[5] caseNames = ["no constant", "restricted constant", "unrestricted constant", "restricted trend", "unrestricted trend"];

    foreach (ci; 0 .. 5) {
        auto ind = independent(cast(Case) ci);
        Prn * prn = gretl_print_new(PrnType.buffer, &err);
        GretlVar * jv = johansen_test(2, list, d, opts[ci], prn);
        string text = gretl_print_get_buffer(prn).to!string;

        check(jv !is null && jv.err == 0 && jv.T == ind.Teff, "johansen_test runs: " ~ caseNames[ci]);
        if (jv is null || jv.err != 0) {
            gretl_print_destroy(prn);
            continue;
        }
        auto jr = jv.jinfo;
        check(mclose(copyOf(jr.R0), ind.R0, 1e-9) && mclose(copyOf(jr.R1), ind.R1, 1e-9), "residual matrices R0 and R1: " ~ caseNames[ci]);
        check(mclose(copyOf(jr.S00), ind.S00, 1e-9) && mclose(copyOf(jr.S01), ind.S01, 1e-9) && mclose(copyOf(jr.S11), ind.S11, 1e-9), "S00, S01, S11: " ~ caseNames[ci]);
        bool evOk = jr.evals.rows * jr.evals.cols == 2;
        foreach (i; 0 .. 2) evOk = evOk && close(jr.evals.val[i], ind.eig[i], 1e-8);
        check(evOk, "eigenvalues: " ~ caseNames[ci]);

        /* trace statistics and the printed table */
        const lmax1 = -ind.Teff * log(1 - ind.eig[1]), lmax0 = -ind.Teff * log(1 - ind.eig[0]);
        const trace0 = lmax0 + lmax1;
        char[64] tbuf, ebuf;
        snprintf_g(tbuf, trace0);
        snprintf_g(ebuf, ind.eig[0]);
        check(strstr(text.ptr, tbuf.ptr) !is null && strstr(text.ptr, ebuf.ptr) !is null, "the table shows the trace statistic and eigenvalue: " ~ caseNames[ci]);

        /* log-likelihood: -(T/2) (n (1 + ln 2 pi) + ln|S00| + sum ln(1 - lambda)) */
        const ll = -(ind.Teff / 2.0) * (2 * (1 + log(2 * PI)) + log(det(ind.S00)) + log(1 - ind.eig[0]) + log(1 - ind.eig[1]));
        check(close(jv.ll, ll, 1e-9), "log-likelihood: " ~ caseNames[ci]);

        /* beta solves the generalized eigenproblem and is normalized by beta' S11 beta = I */
        auto B = copyOf(jr.Beta);
        auto A = copyOf(jr.Alpha);
        auto S11 = ind.S11, S01 = ind.S01, S00 = ind.S00;
        auto lhs = S11 * B * diag2(jr.evals);
        auto rhs = t(S01) * inv(S00) * S01 * B;
        check(mclose(lhs, rhs, 1e-7), "beta columns are eigenvectors of S11^-1 S10 S00^-1 S01: " ~ caseNames[ci]);
        auto bsb = t(B) * S11 * B;
        check(close(bsb[0, 0], 1, 1e-7) && close(bsb[1, 1], 1, 1e-7) && fabs(bsb[0, 1]) < 1e-7, "beta' S11 beta is the identity: " ~ caseNames[ci]);
        gretl_VAR_free(jv);
        gretl_print_destroy(prn);
    }

    /* ---- VECM estimation, rank 1, unrestricted constant ---- */
    auto ind = independent(Case.unrestConst);
    Prn * prn = gretl_print_new(PrnType.buffer, &err);
    GretlVar * v = gretl_VECM(2, 1, list, d, Gretlopt.OPT_NONE, prn, &err);

    check(v !is null && err == 0 && v.err == 0 && v.jinfo.rank == 1, "gretl_VECM (rank 1) runs");
    auto Bt = copyOf(v.jinfo.Beta), At = copyOf(v.jinfo.Alpha);
    check(Bt.nrow == 2 && Bt.ncol == 1 && close(Bt[0, 0], 1.0) && Bt[1, 0] < 0, "Phillips normalization: beta = (1, b)");

    /* the cointegrating vector is the leading eigenvector, rescaled: compare directions with the independent eigenvector */
    auto M = inv(ind.S11) * t(ind.S01) * inv(ind.S00) * ind.S01;
    auto ee = eigenGeneral(M);
    size_t lead = ee.values[0, 0] > ee.values[1, 0] ? 0 : 1;
    const ratio = ee.vectorsRe[1, lead] / ee.vectorsRe[0, lead];
    check(close(Bt[1, 0], ratio, 1e-7), "beta is the leading eigenvector divided by its first element");
    check(fabs(Bt[1, 0] + 0.5) < 0.05, "the estimated cointegrating relation is close to y = 2 x");

    /* alpha = S01 beta (beta' S11 beta)^-1 */
    auto alpha = ind.S01 * Bt * inv(t(Bt) * ind.S11 * Bt);
    check(mclose(At, alpha, 1e-7), "alpha = S01 beta (beta' S11 beta)^-1");

    /* the equations: OLS of (dy - alpha beta' y_{t-1}) on [d_{t-1}, 1] -- compare through the models */
    const Te = ind.Teff;
    auto ec = GretlMatrix(gretl_matrix_alloc(Te, 1));
    auto Dl = GretlMatrix(gretl_matrix_alloc(Te, 4));
    auto Y0 = GretlMatrix(gretl_matrix_alloc(Te, 2));
    foreach (s; 0 .. Te) {
        const t = 2 + s;
        ec[s, 0] = Bt[0, 0] * xs[t - 1] + Bt[1, 0] * ys[t - 1];
        Dl[s, 0] = 1;
        Dl[s, 1] = xs[t - 1] - xs[t - 2];
        Dl[s, 2] = ys[t - 1] - ys[t - 2];
        Dl[s, 3] = ec[s, 0];
        Y0[s, 0] = xs[t] - xs[t - 1];
        Y0[s, 1] = ys[t] - ys[t - 1];
    }
    auto olsFit = solve(Dl, Y0);
    bool eqOk = true;
    foreach (e; 0 .. 2) {
        foreach (k; 0 .. 4) eqOk = eqOk && close(v.models[e].coeff[k], olsFit[k, e], 1e-7);
    }
    check(eqOk, "the equations' coefficients (const, lagged differences, EC term) are the OLS ones given beta");

    /* standard errors of beta: sqrt of the diagonal of (alpha' Omega^-1 alpha (x) H'S11 H)^-1 / df, H'S11H = S11[1,1] */
    auto Om = copyOf(v.S);
    auto aOa = t(At) * inv(Om) * At;
    const varb = 1.0 / (aOa[0, 0] * ind.S11[1, 1]) / v.df;
    check(close(v.jinfo.Bse.val[1], sqrt(varb), 1e-7) && v.jinfo.Bse.val[0] == 0, "standard error of the free element of beta");

    /* the rank-1 log-likelihood */
    const ll1 = -(ind.Teff / 2.0) * (2 * (1 + log(2 * PI)) + log(det(ind.S00)) + log(1 - ind.eig[0]));
    check(close(v.ll, ll1, 1e-9), "VECM log-likelihood for rank 1");
    check(v.jinfo.Gamma !is null && v.A !is null && v.S !is null && v.AIC < 0 && v.BIC > v.AIC, "Gamma, the VAR form, the covariance matrix and information criteria are filled in");
    /* the VAR form: A_1 = I + Pi + Gamma_1, A_2 = -Gamma_1 (companion block) */
    auto Pi = At * t(Bt);
    check(v.A.rows == 4 && v.A.cols == 4, "companion matrix of the VAR form is 4 x 4");
    gretl_print_destroy(prn);

    /* ---- test of a restriction on beta: beta proportional to (1, -1) ---- */
    GretlRestriction * rs = cast(GretlRestriction *) calloc(1, GretlRestriction.sizeof);
    rs.R = gretl_matrix_alloc(1, 2);
    rs.q = gretl_matrix_alloc(1, 1);
    rs.R.val[0] = 1;
    rs.R.val[1] = 1;
    rs.q.val[0] = 0;
    rs.g = 1;
    rs.gmax = 2;
    prn = gretl_print_new(PrnType.buffer, &err);
    auto B0 = copyOf(v.jinfo.Beta);
    err = vecm_test_restriction(v, rs, d, Gretlopt.OPT_NONE, prn);
    check(err == 0, "vecm_test_restriction runs for a simple restriction on beta");
    /* restricted problem: H = (1, -1)', one eigenvalue of (H'S11H)^-1 H'S10 S00^-1 S01 H */
    auto H = GretlMatrix(gretl_matrix_alloc(2, 1));
    H[0, 0] = 1;
    H[1, 0] = -1;
    auto Mr = inv(t(H) * ind.S11 * H) * (t(H) * t(ind.S01) * inv(ind.S00) * ind.S01 * H);
    const lamR = Mr[0, 0];
    const LR = ind.Teff * log((1 - lamR) / (1 - ind.eig[0]));
    check(close(rs.test, LR, 1e-7) && close(rs.pval, chisq_cdf_comp(1, LR), 1e-9) && LR > 0, "the LR statistic and p-value of the restriction (1 df)");
    check(mclose(copyOf(v.jinfo.Beta), B0), "vecm_test_restriction leaves the model's beta as it was");
    gretl_print_destroy(prn);

    /* ---- estimation subject to the restriction (from scratch, as var.d does it): beta is (1, -1) up to scale ---- */
    prn = gretl_print_new(PrnType.buffer, &err);
    GretlVar * vr = johansen_wrapper(VarTask.vecm_estimate, 2, 1, null, list, rs, d, Gretlopt.OPT_NONE, prn, &err);
    check(vr !is null && err == 0 && vr.err == 0, "a VECM estimated subject to a simple restriction on beta (johansen_wrapper)");
    auto Br = copyOf(vr.jinfo.Beta);
    check(Br.nrow == 2 && close(Br[0, 0], -Br[1, 0], 1e-9) && close(Br[0, 0], 1.0, 1e-9), "the restricted beta is (1, -1)");
    const llr1 = -(ind.Teff / 2.0) * (2 * (1 + log(2 * PI)) + log(det(ind.S00)) + log(1 - lamR));
    check(close(vr.ll, llr1, 1e-9) && vr.ll < ll1, "the restricted log-likelihood is lower than the unrestricted one and equals the eigenvalue formula");
    check(close(vr.jinfo.ll0, ll1, 1e-9) && vr.jinfo.lrdf == 1, "the unrestricted log-likelihood and the df of the restriction are recorded");
    check(vr.models !is null && vr.models[0].ncoeff == 4 && vr.A !is null && vr.jinfo.Bse !is null, "the equations, the VAR form and the standard errors are filled in");
    gretl_VAR_free(vr);
    gretl_print_destroy(prn);

    /* ---- bootstrap impulse responses of a VECM (johansen_boot_round) ---- */
    GretlVar * vb = gretl_VECM(2, 1, list, d, Gretlopt.OPT_Q, null, &err);

    execute_set("boot_iters", "499", null, Gretlopt.OPT_NONE, null);
    err = 0;
    auto point = gretl_VAR_get_impulse_response(vb, 0, 0, 6, 0.0, d, &err);
    check(point !is null && err == 0 && point.rows == 6, "point impulse responses of the VECM");
    auto bands = gretl_VAR_get_impulse_response(vb, 0, 0, 6, 0.1, d, &err);
    check(bands !is null && err == 0 && bands.rows == 6 && bands.cols == 3, "bootstrap impulse responses of the VECM (uses johansen_boot_round)");
    if (bands !is null && bands.rows == 6) {
        bool within = true;
        foreach (h; 0 .. 6) {
            const pt = gretl_matrix_get(bands, h, 0), lo = gretl_matrix_get(bands, h, 1), hi = gretl_matrix_get(bands, h, 2);

            within = within && lo <= hi && lo <= pt + 1e-9 && pt <= hi + 1e-9;
        }
        check(within, "the bootstrap band contains the point response");
        check(close(gretl_matrix_get(bands, 0, 0), gretl_matrix_get(point, 0, 0), 1e-9), "the first column of the bands is the point response");
    }
    check(johansen_boot_round(vb, d) == 0, "johansen_boot_round can be called directly on an estimated VECM");
    gretl_VAR_free(vb);
    gretl_VAR_free(v);

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}

/* diag of the eigenvalue vector as a matrix */
GretlMatrix diag2(const(GretlMatrixRaw) * ev)
{
    const n = ev.rows * ev.cols;
    auto m = GretlMatrix(gretl_matrix_alloc(n, n));

    foreach (j; 0 .. n) foreach (i; 0 .. n) m[i, j] = (i == j) ? ev.val[i] : 0.0;
    return m;
}

/* the text of "%#.5g" for the table checks */
void snprintf_g(ref char[64] buf, double x)
{
    import core.stdc.stdio : snprintf;
    snprintf(buf.ptr, buf.length, "%#.5g", x);
}
