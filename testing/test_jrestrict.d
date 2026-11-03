/* Tests for plugin/jrestrict.d: general restrictions on beta (and alpha) in a VECM (Boswijk and Doornik 2004): the identification
   analysis through the rank of the Jacobian, the switching algorithm and the LBFGS alternative, removal of column scalings,
   restrictions on alpha combined with restrictions on beta, the variances from the information matrix, the printed result, and
   the bootstrap use.

   The entry point is general_vecm_analysis(jvar, rset, dset, prn); with the option OPT_F in the restriction set it replaces beta,
   alpha, Omega and the log-likelihood in the VECM (the old log-likelihood goes to jinfo.ll0, the degrees of freedom of the
   restriction to jinfo.lrdf), which is what the tests look at. Each result is compared with an independent computation with the
   matrix operators of dgretl:
     - homogeneous restrictions common to all columns of beta: the eigenvalue problem for H' S11 H (the "simple" method of
       johansen.d), whose log-likelihood the general algorithm must reach;
     - a restriction that just identifies beta (normalizations, df = 0): the unrestricted estimates and log-likelihood;
     - a fully specified beta: alpha = S01 beta (beta' S11 beta)^-1, Omega, the log-likelihood and the LR statistic;
     - alpha_x = 0 together with a normalization of beta: the log-likelihood of the alpha-restricted eigenproblem
       (Johansen 1995, section 8.2.1);
     - the standard error of the free element of beta from the information matrix.

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs, floor, log, sqrt;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;
import core.stdc.stdlib : calloc;
import core.stdc.string : strstr;

import std.conv : to;
import std.math : PI;

import dgretl;

import gretl.dataset : create_new_dataset;
import gretl.gretl_list : gretl_list_new;
import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_get_buffer, gretl_print_new;
import gretl.gretl_prnh : PrnType;
import gretl.gretl_restricth : GretlRestriction;
import gretl.libgretlh : DatasetRaw, Gretlopt, TsCodes;
import gretl.libset : execute_set;
import gretl.plugin.jrestrict : general_vecm_analysis;
import gretl.pvalues : chisq_cdf_comp;
import gretl.var : VarTask, gretl_VAR_free, gretl_VAR_get_impulse_response, gretl_VECM, johansen_wrapper;
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

/* the moment matrices for a VAR of order 2 in levels with an unrestricted constant (case 3), for the series in `cols` */
struct Mom {
    GretlMatrix S00, S01, S11;
    double[] eig;                 // descending
    int Teff, n;
}

Mom moments(const(double)[][] cols)
{
    Mom m;
    const n = cols.length, T = cols[0].length, t0 = 2, Te = T - t0;
    auto Z0 = GretlMatrix(gretl_matrix_alloc(cast(int) Te, cast(int) n));
    auto Z1 = GretlMatrix(gretl_matrix_alloc(cast(int) Te, cast(int) n));
    auto Z2 = GretlMatrix(gretl_matrix_alloc(cast(int) Te, cast(int) n + 1));

    foreach (s; 0 .. Te) {
        const t = t0 + s;

        foreach (i; 0 .. n) {
            Z0[s, i] = cols[i][t] - cols[i][t - 1];
            Z1[s, i] = cols[i][t - 1];
            Z2[s, i] = cols[i][t - 1] - cols[i][t - 2];
        }
        Z2[s, n] = 1;
    }
    auto R0 = resid(Z0, Z2), R1 = resid(Z1, Z2);
    m.S00 = t(R0) * R0 / cast(double) Te;
    m.S01 = t(R0) * R1 / cast(double) Te;
    m.S11 = t(R1) * R1 / cast(double) Te;
    m.Teff = cast(int) Te;
    m.n = cast(int) n;
    auto ev = eigenGeneral(inv(m.S11) * t(m.S01) * inv(m.S00) * m.S01, false).values;
    foreach (i; 0 .. ev.nrow) m.eig ~= ev[i, 0];
    m.eig.sort!"a > b";
    return m;
}

import std.algorithm.sorting : sort;

double llFromEig(Mom m, const(double)[] eig, int r)
{
    double s = 0;

    foreach (i; 0 .. r) s += log(1 - eig[i]);
    return -(m.Teff / 2.0) * (m.n * (1 + log(2 * PI)) + log(det(m.S00)) + s);
}

/* the log-likelihood of the model with beta = H phi (a common homogeneous restriction): eigenvalues of (H'S11H)^-1 H'S10 S00^-1 S01 H */
double llHomogeneous(Mom m, GretlMatrix H, int r)
{
    auto M = inv(t(H) * m.S11 * H) * (t(H) * t(m.S01) * inv(m.S00) * m.S01 * H);
    auto ev = eigenGeneral(M, false).values;
    double[] e;

    foreach (i; 0 .. ev.nrow) e ~= ev[i, 0];
    e.sort!"a > b";
    return llFromEig(m, e, r);
}

GretlRestriction * makeSet(GretlMatrix R, GretlMatrix q, GretlMatrix Ra = GretlMatrix.init, GretlMatrix qa = GretlMatrix.init, Gretlopt opt = Gretlopt.OPT_NONE)
{
    GretlRestriction * rs = cast(GretlRestriction *) calloc(1, GretlRestriction.sizeof);

    if (R.nrow > 0) {
        rs.R = gretl_matrix_copy(R.ptr);
        rs.q = gretl_matrix_copy(q.ptr);
        rs.g = cast(int) R.nrow;
        rs.gmax = cast(int) R.ncol;
    }
    if (Ra.nrow > 0) {
        rs.Ra = gretl_matrix_copy(Ra.ptr);
        rs.qa = gretl_matrix_copy(qa.ptr);
    }
    rs.opt = opt;
    return rs;
}

GretlMatrix rowvec(double[] a)
{
    auto m = GretlMatrix(gretl_matrix_alloc(1, cast(int) a.length));

    foreach (j, x; a) m[0, j] = x;
    return m;
}

GretlMatrix colvec(double[] a)
{
    auto m = GretlMatrix(gretl_matrix_alloc(cast(int) a.length, 1));

    foreach (i, x; a) m[i, 0] = x;
    return m;
}

GretlMatrix rows(double[][] a)
{
    auto m = GretlMatrix(gretl_matrix_alloc(cast(int) a.length, cast(int) a[0].length));

    foreach (i, row; a) foreach (j, x; row) m[i, j] = x;
    return m;
}

DatasetRaw * makeData(const(double)[][] cols, string[] names)
{
    DatasetRaw * d = create_new_dataset(cast(int) cols.length + 1, cast(int) cols[0].length, 0);

    d.pd = 1;
    d.structure = TsCodes.time_series;
    foreach (i, nm; ["const"] ~ names) {
        d.varname[i][0 .. nm.length] = nm[];
        d.varname[i][nm.length] = 0;
    }
    foreach (t; 0 .. cols[0].length) {
        d.Z[0][t] = 1;
        foreach (i; 0 .. cols.length) d.Z[i + 1][t] = cols[i][t];
    }
    return d;
}

int * listOf(int n)
{
    int * list = gretl_list_new(n);

    foreach (i; 1 .. n + 1) list[i] = i;
    return list;
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

    /* ---- the two-variable system: x a random walk, y = 2 x + AR(1) noise (r = 1) ---- */
    enum T2 = 120;
    auto xs = new double[T2], ys = new double[T2];
    double u = 0;

    foreach (t; 0 .. T2) {
        xs[t] = (t == 0 ? 0 : xs[t - 1]) + rnd();
        u = 0.5 * u + 0.4 * rnd();
        ys[t] = 2.0 * xs[t] + u;
    }
    auto d2 = makeData([xs, ys], ["x", "y"]);
    int * list2 = listOf(2);
    int err = 0;
    auto m2 = moments([xs, ys]);
    const ll_u1 = llFromEig(m2, m2.eig, 1);

    GretlVar * vnew() { int e; return gretl_VECM(2, 1, list2, d2, Gretlopt.OPT_Q, null, &e); }

    /* beta_x = 1: a normalization, so the model is just identified (df = 0) and nothing changes */
    {
        auto v = vnew();
        auto Bphillips = copyOf(v.jinfo.Beta), Aun = copyOf(v.jinfo.Alpha);
        auto rs = makeSet(rowvec([1, 0]), colvec([1]), GretlMatrix.init, GretlMatrix.init, Gretlopt.OPT_F);
        Prn * prn = gretl_print_new(PrnType.buffer, &err);

        err = general_vecm_analysis(v, rs, d2, prn);
        check(err == 0, "general_vecm_analysis: beta_x = 1 runs");
        check(strstr(gretl_print_get_buffer(prn), "Model is fully identified") !is null && v.jinfo.lrdf == 0, "the model is reported fully identified, df = 0");
        check(close(v.ll, ll_u1, 1e-6) && close(v.jinfo.ll0, ll_u1, 1e-9), "log-likelihood unchanged by a just-identifying normalization");
        auto Bn = copyOf(v.jinfo.Beta), An = copyOf(v.jinfo.Alpha);
        check(close(Bn[0, 0], 1, 1e-9) && close(Bn[1, 0], Bphillips[1, 0], 1e-5) && mclose(An, Aun, 1e-4), "beta = (1, b) and alpha agree with the Phillips-normalized VECM");
        /* standard error of b from the information matrix: df = T - floor(K / p), K = (order p + ifc) p + blen + alen */
        auto Om = copyOf(v.S);
        auto aOa = t(An) * inv(Om) * An;
        const df = m2.Teff - cast(int) floor(((1 * 2 + 1) * 2 + 1 + 2) / 2.0);
        const se = sqrt(1.0 / (aOa[0, 0] * m2.S11[1, 1]) / df);
        check(v.jinfo.Bse !is null && v.jinfo.Bse.val[0] == 0 && close(v.jinfo.Bse.val[1], se, 1e-5), "standard error of the free element of beta from the information matrix");
        check(v.jinfo.Ase !is null && v.jinfo.Ase.val[0] > 0 && v.jinfo.Ase.val[1] > 0, "standard errors of alpha");
        gretl_print_destroy(prn);
        gretl_VAR_free(v);
    }

    /* beta fully specified: (1, -2)' */
    {
        auto v = vnew();
        auto rs = makeSet(rows([[1, 0], [0, 1]]), colvec([1, -2]), GretlMatrix.init, GretlMatrix.init, Gretlopt.OPT_F);
        Prn * prn = gretl_print_new(PrnType.buffer, &err);

        err = general_vecm_analysis(v, rs, d2, prn);
        check(err == 0, "general_vecm_analysis: beta fully specified runs");
        auto B = colvec([1, -2]);
        auto alpha = m2.S01 * B / (t(B) * m2.S11 * B)[0, 0];
        auto Om = m2.S00 - alpha * t(alpha) * (t(B) * m2.S11 * B)[0, 0];
        const llr = -(m2.Teff / 2.0) * (2 * (1 + log(2 * PI)) + log(det(Om)));
        check(close(v.ll, llr, 1e-9) && v.jinfo.lrdf == 1 && close(v.jinfo.ll0, ll_u1, 1e-9), "beta = (1, -2): log-likelihood, df = 1 and the unrestricted log-likelihood recorded");
        check(mclose(copyOf(v.jinfo.Alpha), alpha, 1e-9) && mclose(copyOf(v.jinfo.Beta), B, 1e-12) && mclose(copyOf(v.S), Om, 1e-9), "alpha = S01 beta / (beta' S11 beta), beta and Omega");
        gretl_print_destroy(prn);
        /* the test only (no OPT_F): the printed LR statistic and the recorded results */
        rs = makeSet(rows([[1, 0], [0, 1]]), colvec([1, -2]));
        prn = gretl_print_new(PrnType.buffer, &err);
        v = vnew();
        err = general_vecm_analysis(v, rs, d2, prn);
        const LR = 2 * (ll_u1 - llr);
        check(err == 0 && close(rs.test, LR, 1e-7) && close(rs.pval, chisq_cdf_comp(1, LR), 1e-8) && close(v.ll, ll_u1, 1e-9), "the LR test of beta = (1, -2) (1 df) and the VECM left alone");
        check(strstr(gretl_print_get_buffer(prn), "Restricted loglikelihood") !is null, "the output shows the restricted log-likelihood");
        gretl_print_destroy(prn);
        gretl_VAR_free(v);
    }

    /* homogeneous restriction beta proportional to (1, -1): the general algorithm reaches the eigenvalue-method likelihood */
    {
        auto H = colvec([1, -1]);
        const llr = llHomogeneous(m2, H, 1);

        foreach (variant; 0 .. 3) {
            auto v = vnew();
            Gretlopt o = Gretlopt.OPT_F;
            string label = "switching algorithm";

            if (variant == 1) { o |= Gretlopt.OPT_J; label = "switching algorithm with simulated annealing initialization"; }
            if (variant == 2) { o |= Gretlopt.OPT_N; label = "switching algorithm without scale removal"; }
            auto rs = makeSet(rowvec([1, 1]), colvec([0]), GretlMatrix.init, GretlMatrix.init, o);
            Prn * prn = gretl_print_new(PrnType.buffer, &err);

            err = general_vecm_analysis(v, rs, d2, prn);
            check(err == 0, "general_vecm_analysis: beta ~ (1, -1), " ~ label);
            check(close(v.ll, llr, 1e-5), "log-likelihood equals the eigenvalue method: " ~ label);
            auto B = copyOf(v.jinfo.Beta);
            check(close(B[0, 0] + B[1, 0], 0, 1e-6) && v.jinfo.lrdf == 1, "beta is proportional to (1, -1), df = 1: " ~ label);
            gretl_print_destroy(prn);
            gretl_VAR_free(v);
        }
    }

    /* LBFGS needs an identified model: with beta ~ (1, -1) and alpha unrestricted the scale of beta is not identified, which is an error */
    {
        auto v = vnew();
        auto rs = makeSet(rowvec([1, 1]), colvec([0]), GretlMatrix.init, GretlMatrix.init, Gretlopt.OPT_F | Gretlopt.OPT_L);
        Prn * prn = gretl_print_new(PrnType.buffer, &err);

        err = general_vecm_analysis(v, rs, d2, prn);
        check(err != 0 && strstr(gretl_print_get_buffer(prn), "Model is not fully identified") !is null, "LBFGS on a model that is not fully identified is refused");
        gretl_print_destroy(prn);
        gretl_VAR_free(v);
    }

    /* alpha_x = 0 together with the normalization beta_x = 1: identified, so the LBFGS alternative (option OPT_L, L-BFGS-B with
       analytical derivatives) can be used as well as the switching algorithm */
    foreach (useLbfgs; 0 .. 2) {
        const label = useLbfgs ? "LBFGS" : "switching algorithm";
        auto v = vnew();
        auto rs = makeSet(rowvec([1, 0]), colvec([1]), rowvec([1, 0]), colvec([0]), useLbfgs ? Gretlopt.OPT_F | Gretlopt.OPT_L : Gretlopt.OPT_F);
        Prn * prn = gretl_print_new(PrnType.buffer, &err);

        err = general_vecm_analysis(v, rs, d2, prn);
        check(err == 0, "general_vecm_analysis: alpha_x = 0 and beta_x = 1 runs (" ~ label ~ ")");
        // the alpha-restricted eigenproblem
        auto Ap = rowvec([1, 0]);
        auto C = t(Ap) * inv(Ap * m2.S00 * t(Ap)) * Ap;
        auto S00a = m2.S00 - m2.S00 * C * m2.S00, S11a = m2.S11 - t(m2.S01) * C * m2.S01, S01a = m2.S01 - m2.S00 * C * m2.S01;
        auto A = colvec([0, 1]);
        auto Abar = A * inv(t(A) * A);
        auto ev = eigenGeneral(inv(S11a) * t(S01a) * Abar * inv(t(Abar) * S00a * Abar) * t(Abar) * S01a, false).values;
        const lam = ev[0, 0] > ev[1, 0] ? ev[0, 0] : ev[1, 0];
        const llr = -(m2.Teff / 2.0) * (2 * (1 + log(2 * PI)) + log(det(m2.S00)) + log(1 - lam));
        check(close(v.ll, llr, useLbfgs ? 1e-7 : 1e-5), "log-likelihood equals that of the alpha-restricted eigenproblem (" ~ label ~ ")");
        auto An = copyOf(v.jinfo.Alpha), Bn = copyOf(v.jinfo.Beta);
        check(fabs(An[0, 0]) < 1e-8 && fabs(An[1, 0]) > 1e-4 && close(Bn[0, 0], 1, 1e-9), "alpha_x is zero, alpha_y is not, beta_x = 1 (" ~ label ~ ")");
        check(v.jinfo.lrdf == 1, "one degree of freedom for the alpha restriction (" ~ label ~ ")");
        gretl_print_destroy(prn);
        gretl_VAR_free(v);
    }

    /* a non-homogeneous restriction on alpha is refused */
    {
        auto v = vnew();
        auto rs = makeSet(rowvec([1, 0]), colvec([1]), rowvec([1, 0]), colvec([0.5]), Gretlopt.OPT_F);
        Prn * prn = gretl_print_new(PrnType.buffer, &err);

        err = general_vecm_analysis(v, rs, d2, prn);
        check(err != 0, "a non-homogeneous restriction on alpha is an error");
        gretl_print_destroy(prn);
        gretl_VAR_free(v);
    }

    /* ---- three variables, r = 2: y = 2 x + u1, z = -x + u2 ---- */
    enum T3 = 160;
    auto x3 = new double[T3], y3 = new double[T3], z3 = new double[T3];
    double u1 = 0, u2 = 0;

    foreach (t; 0 .. T3) {
        x3[t] = (t == 0 ? 0 : x3[t - 1]) + rnd();
        u1 = 0.5 * u1 + 0.3 * rnd();
        u2 = 0.4 * u2 + 0.3 * rnd();
        y3[t] = 2.0 * x3[t] + u1;
        z3[t] = -1.0 * x3[t] + u2;
    }
    auto d3 = makeData([x3, y3, z3], ["x", "y", "z"]);
    int * list3 = listOf(3);
    auto m3 = moments([x3, y3, z3]);
    const ll_u2 = llFromEig(m3, m3.eig, 2);

    GretlVar * vnew3() { int e; return gretl_VECM(2, 2, list3, d3, Gretlopt.OPT_Q, null, &e); }

    {
        auto v = vnew3();
        check(v !is null && v.jinfo.rank == 2 && close(v.ll, ll_u2, 1e-9), "gretl_VECM with r = 2 on three variables, log-likelihood from the eigenvalues");
        gretl_VAR_free(v);
    }

    /* a restriction common to both columns of beta: beta_x + beta_y = 0 (the same row for each column) */
    {
        // H: 3 x 2 with beta' (1, 1, 0) = 0
        auto H = rows([[1, 0], [-1, 0], [0, 1]]);
        const llr = llHomogeneous(m3, H, 2);
        auto v = vnew3();
        auto rs = makeSet(rowvec([1, 1, 0]), colvec([0]), GretlMatrix.init, GretlMatrix.init, Gretlopt.OPT_F);
        Prn * prn = gretl_print_new(PrnType.buffer, &err);

        err = general_vecm_analysis(v, rs, d3, prn);
        check(err == 0, "general_vecm_analysis: a common restriction with r = 2 runs");
        check(close(v.ll, llr, 1e-5) && v.jinfo.lrdf == 2, "log-likelihood equals the eigenvalue method; df = r (p - s) = 2");
        auto B = copyOf(v.jinfo.Beta);
        check(B.nrow == 3 && B.ncol == 2 && close(B[0, 0] + B[1, 0], 0, 1e-6) && close(B[0, 1] + B[1, 1], 0, 1e-6), "both columns of beta satisfy the restriction");
        gretl_print_destroy(prn);
        gretl_VAR_free(v);
    }

    /* beta = [I | free]: four normalizations that just identify beta */
    {
        auto v = vnew3();
        auto R = rows([[1, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 1, 0]]);
        // vec(beta) is column by column: beta_11 = 1, beta_21 = 0, beta_12 = 0, beta_22 = 1
        auto rs = makeSet(R, colvec([1, 0, 0, 1]), GretlMatrix.init, GretlMatrix.init, Gretlopt.OPT_F);
        Prn * prn = gretl_print_new(PrnType.buffer, &err);

        err = general_vecm_analysis(v, rs, d3, prn);
        check(err == 0, "general_vecm_analysis: beta = [I | free] with r = 2 runs");
        check(close(v.ll, ll_u2, 1e-5) && v.jinfo.lrdf == 0, "just-identifying normalizations leave the log-likelihood and df = 0");
        auto B = copyOf(v.jinfo.Beta);
        check(close(B[0, 0], 1, 1e-9) && fabs(B[1, 0]) < 1e-9 && fabs(B[0, 1]) < 1e-9 && close(B[1, 1], 1, 1e-9), "the normalized elements are as specified");
        // the free elements: the third row; y = 2x + u1 and z = -x + u2 give beta_31 = -(1/2)... check that beta' y_t is stationary-ish: compare with the Phillips-normalized VECM
        auto vp = vnew3();
        auto Bp = copyOf(vp.jinfo.Beta);
        check(mclose(B, Bp, 1e-3), "the free elements agree with the Phillips-normalized VECM");
        gretl_VAR_free(vp);
        gretl_print_destroy(prn);
        gretl_VAR_free(v);
    }

    /* identification: a single normalization with r = 2 leaves the model under-identified */
    {
        auto v = vnew3();
        auto rs = makeSet(rowvec([1, 0, 0, 0, 0, 0]), colvec([1]), GretlMatrix.init, GretlMatrix.init, Gretlopt.OPT_NONE);
        Prn * prn = gretl_print_new(PrnType.buffer, &err);

        err = general_vecm_analysis(v, rs, d3, prn);
        check(strstr(gretl_print_get_buffer(prn), "Model is not fully identified") !is null, "an under-identified model is reported as such");
        gretl_print_destroy(prn);
        gretl_VAR_free(v);
    }

    /* ---- estimation from scratch subject to a general restriction (the way var.d does it), and the bootstrap use ---- */
    {
        auto rs = makeSet(rowvec([1, 0]), colvec([1]), GretlMatrix.init, GretlMatrix.init, Gretlopt.OPT_F);   // OPT_F: full estimation
        Prn * prn = gretl_print_new(PrnType.buffer, &err);
        GretlVar * v = johansen_wrapper(VarTask.vecm_estimate, 2, 1, null, list2, rs, d2, Gretlopt.OPT_NONE, prn, &err);

        check(v !is null && err == 0 && v.err == 0, "a VECM estimated subject to beta_x = 1 from scratch (johansen_wrapper)");
        auto Bp = copyOf(vnew().jinfo.Beta);
        check(close(v.jinfo.Beta.val[0], 1, 1e-9) && close(v.jinfo.Beta.val[1], Bp[1, 0], 1e-4) && close(v.ll, ll_u1, 1e-5), "its beta and log-likelihood (df = 0)");
        check(v.jinfo.R !is null && v.jinfo.q !is null, "the restriction is recorded in the VECM");
        check(v.models !is null && v.models[0].ncoeff == 4 && v.A !is null && v.A.rows == 4, "the equations and the VAR form have the right size");
        gretl_print_destroy(prn);
        execute_set("boot_iters", "499", null, Gretlopt.OPT_NONE, null);
        err = 0;
        auto point = gretl_VAR_get_impulse_response(v, 0, 0, 5, 0.0, d2, &err);
        auto bands = gretl_VAR_get_impulse_response(v, 0, 0, 5, 0.1, d2, &err);
        check(point !is null && bands !is null && err == 0 && bands.rows == 5 && bands.cols == 3, "bootstrap impulse responses of a restricted VECM (the OPT_B path of the general algorithm)");
        if (bands !is null && bands.rows == 5) {
            bool within = true;
            foreach (h; 0 .. 5) {
                const pt = gretl_matrix_get(bands, h, 0), lo = gretl_matrix_get(bands, h, 1), hi = gretl_matrix_get(bands, h, 2);

                within = within && lo <= hi && lo <= pt + 1e-9 && pt <= hi + 1e-9;
            }
            check(within, "the bands contain the point response");
        }
        gretl_VAR_free(v);
    }

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
