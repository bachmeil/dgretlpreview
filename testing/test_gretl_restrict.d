/* Smoke tests for gretl_restrict.d on a 30-observation OLS model: a hand-built restriction set (R and q matrices) for
   an OLS model, the accessors, and the finalizing Wald/F test (the printed F for a single restriction against the
   by-hand t-ratio squared, and a joint test).

   Not exercised: parsing of restriction lines (needs the genr parser), system and VECM restrictions, nonlinear restrictions through a function, bootstrap options and
   gretl_sum_test.

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs;
import core.stdc.stdio : _IOLBF, printf, setvbuf, snprintf, stdout;
import core.stdc.string : strlen, strstr;

import gretl.dataset : create_new_dataset;
import gretl.estimate : lsq;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_matrixh : GretlMatrixRaw;
import gretl.gretl_model : clear_model, gretl_model_get_vcv_element;
import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_get_buffer, gretl_print_new;
import gretl.gretl_prnh : PrnType;
import gretl.gretl_restrict;
import gretl.gretl_restricth : GretlRestriction;
import gretl.libgretlh;
import gretl.objstackh : GretlObjType;

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

bool has(const(char) * buf, const(char)[] needle)
{
    import std.string : toStringz;
    return buf !is null && strstr(buf, toStringz(needle)) !is null;
}

GretlRestriction * make_set(Model * pmod, const(double)[] rvals, const(double)[] qvals, int g = 1)
{
    import core.stdc.stdlib : calloc;
    import gretl.gretl_matrix : gretl_matrix_alloc;

    GretlRestriction * rs = cast(GretlRestriction *) calloc(1, GretlRestriction.sizeof);
    GretlMatrixRaw * R = gretl_matrix_alloc(g, 3);
    GretlMatrixRaw * q = gretl_matrix_alloc(g, 1);

    foreach (i, v; rvals) {
        /* row-major input into the column-major matrix */
        R.val[(i % 3) * g + i / 3] = v;
    }
    foreach (i, v; qvals) {
        q.val[i] = v;
    }
    rs.g = g;
    rs.gmax = 3;
    rs.R = R;
    rs.q = q;
    rs.obj = pmod;
    rs.otype = GretlObjType.eqn;
    return rs;
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    enum n = 30;
    DatasetRaw * d = create_new_dataset(4, n, 0);

    foreach (i, nm; ["const", "y", "x1", "x2"]) {
        d.varname[i][0 .. nm.length] = nm[];
        d.varname[i][nm.length] = 0;
    }
    uint seed = 777;
    double rnd()
    {
        seed = seed * 1103515245u + 12345u;
        return ((seed >> 8) & 0xffff) / 65536.0 - 0.5;
    }
    foreach (t; 0 .. n) {
        d.Z[0][t] = 1;
        d.Z[2][t] = 3.0 * rnd();
        d.Z[3][t] = 3.0 * rnd();
        d.Z[1][t] = 1.0 + 2.0 * d.Z[2][t] + 0.5 * d.Z[3][t] + rnd();
    }
    d.t1 = 0;
    d.t2 = n - 1;

    int[5] list = [4, 1, 0, 2, 3];
    Model m = lsq(list.ptr, d, GretlCmdIndex.ols, Gretlopt.OPT_NONE);

    check(m.errcode == 0 && m.ncoeff == 3, "OLS y on const x1 x2");

    /* a set built by hand (parsing a restriction line needs the genr parser): b[2] = 2 */
    GretlRestriction * rs = make_set(&m, [0.0, 1.0, 0.0], [2.0]);

    check(rs !is null, "restriction set built");

    const(GretlMatrixRaw) * R = rset_get_R_matrix(rs);
    const(GretlMatrixRaw) * q = rset_get_q_matrix(rs);

    check(R !is null && R.rows == 1 && R.cols == 3 && q !is null && q.val[0] == 2, "R and q accessors");
    check(gretl_restriction_get_type(rs) == GretlObjType.eqn, "the set refers to an equation");

    int err = 0;
    Prn * prn = gretl_print_new(PrnType.buffer, &err);
    err = gretl_restriction_finalize(rs, d, Gretlopt.OPT_NONE, prn);
    check(err == 0, "finalize runs the test");

    double se = sqrt_(gretl_model_get_vcv_element(&m, 1, 1, 3));
    double t = (m.coeff[1] - 2.0) / se;
    char[64] fs;
    snprintf(fs.ptr, fs.length, "F(1, 27) = %g", t * t);
    check(has(gretl_print_get_buffer(prn), fs[0 .. strlen(fs.ptr)]), "the printed F equals the squared t-ratio");
    gretl_print_destroy(prn);

    /* two restrictions: b[2] = 2, b[3] = 0.5 */
    rs = make_set(&m, [0.0, 1.0, 0.0, 0.0, 0.0, 1.0], [2.0, 0.5], 2);
    prn = gretl_print_new(PrnType.buffer, &err);
    err = gretl_restriction_finalize(rs, d, Gretlopt.OPT_NONE, prn);
    check(err == 0 && has(gretl_print_get_buffer(prn), "F(2, 27)"), "joint test has two numerator degrees of freedom");
    gretl_print_destroy(prn);

    /* a null set is an error */
    check(gretl_restriction_finalize(null, d, Gretlopt.OPT_NONE, null) != 0, "finalizing a null set is an error");

    clear_model(&m);
    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}

double sqrt_(double x)
{
    import core.stdc.math : sqrt;
    return sqrt(x);
}
