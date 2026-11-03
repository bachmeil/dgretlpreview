/* Smoke tests for compare.d on a 30-observation OLS model: the F test for omitting a variable (the printed F against the
   by-hand restricted-model F), adding a variable, RESET, the non-linearity (squares) test, the LM autocorrelation
   test and the Chow test (against the by-hand split-sample F).

   Not exercised: stepwise elimination (plugin), QLR and CUSUM (plugins/gnuplot), panel tests, leverage, VIF and
   BKW, IV autocorrelation, and the "compare models" functions.

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs;
import core.stdc.stdio : _IOLBF, printf, setvbuf, snprintf, stdout;
import core.stdc.string : strlen, strstr;

import gretl.compare;
import gretl.dataset : create_new_dataset;
import gretl.estimate : lsq;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_model : clear_model;
import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_get_buffer, gretl_print_new;
import gretl.gretl_prnh : PrnType;
import gretl.libgretlh;
import gretl.compareh : ModelAuxCode;

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

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    enum n = 30;
    DatasetRaw * d = create_new_dataset(5, n, 0);

    foreach (i, nm; ["const", "y", "x1", "x2", "x3"]) {
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
        d.Z[4][t] = 3.0 * rnd();
        d.Z[1][t] = 1.0 + 2.0 * d.Z[2][t] + 0.5 * d.Z[3][t] + rnd();
    }
    d.t1 = 0;
    d.t2 = n - 1;

    int[5] list = [4, 1, 0, 2, 3];
    Model m = lsq(list.ptr, d, GretlCmdIndex.ols, Gretlopt.OPT_NONE);

    check(m.errcode == 0 && m.ncoeff == 3, "OLS y on const x1 x2");

    int err = 0;
    Prn * prn = gretl_print_new(PrnType.buffer, &err);

    /* omit x2: F against the restricted regression */
    int[4] rlist = [3, 1, 0, 2];
    Model r = lsq(rlist.ptr, d, GretlCmdIndex.ols, Gretlopt.OPT_NONE);
    double Fhand = ((r.ess - m.ess) / 1.0) / (m.ess / m.dfd);
    int[2] omitv = [1, 3];

    err = omit_test(&m, omitv.ptr, d, Gretlopt.OPT_NONE, prn);
    check(err == 0, "omit_test runs");
    char[64] fs;
    snprintf(fs.ptr, fs.length, "F(1, 27) = %g", Fhand);
    check(has(gretl_print_get_buffer(prn), fs[0 .. strlen(fs.ptr)]), "omit_test prints the by-hand F statistic");
    clear_model(&r);
    gretl_print_destroy(prn);

    /* add x3 */
    prn = gretl_print_new(PrnType.buffer, &err);
    int[2] addv = [1, 4];
    err = add_test(&m, addv.ptr, d, Gretlopt.OPT_NONE, prn);
    const(char) * buf = gretl_print_get_buffer(prn);
    check(err == 0 && buf !is null && buf[0] != 0, "add_test runs and prints");
    gretl_print_destroy(prn);

    /* RESET */
    prn = gretl_print_new(PrnType.buffer, &err);
    err = reset_test(&m, d, Gretlopt.OPT_NONE, prn);
    buf = gretl_print_get_buffer(prn);
    check(err == 0 && buf !is null && buf[0] != 0, "reset_test runs and prints");
    gretl_print_destroy(prn);

    /* non-linearity (squares) */
    prn = gretl_print_new(PrnType.buffer, &err);
    err = nonlinearity_test(&m, d, ModelAuxCode.sq, Gretlopt.OPT_NONE, prn);
    buf = gretl_print_get_buffer(prn);
    check(err == 0 && buf !is null && buf[0] != 0, "nonlinearity_test (squares) runs and prints");
    gretl_print_destroy(prn);

    /* LM autocorrelation test */
    prn = gretl_print_new(PrnType.buffer, &err);
    err = autocorr_test(&m, 1, d, Gretlopt.OPT_NONE, prn);
    buf = gretl_print_get_buffer(prn);
    check(err == 0 && buf !is null && buf[0] != 0, "autocorr_test (order 1) runs and prints");
    gretl_print_destroy(prn);

    /* Chow test at observation 16 against the by-hand split F */
    d.t2 = 15;
    Model s1 = lsq(list.ptr, d, GretlCmdIndex.ols, Gretlopt.OPT_NONE);
    d.t1 = 16;
    d.t2 = n - 1;
    Model s2 = lsq(list.ptr, d, GretlCmdIndex.ols, Gretlopt.OPT_NONE);
    d.t1 = 0;
    double Fchow = ((m.ess - s1.ess - s2.ess) / 3.0) / ((s1.ess + s2.ess) / (n - 6.0));

    prn = gretl_print_new(PrnType.buffer, &err);
    err = chow_test(16, &m, d, Gretlopt.OPT_NONE, prn);
    check(err == 0, "chow_test runs");
    snprintf(fs.ptr, fs.length, "F(3, 24) = %g", Fchow);
    check(has(gretl_print_get_buffer(prn), fs[0 .. strlen(fs.ptr)]), "chow_test prints the by-hand split-sample F");
    gretl_print_destroy(prn);
    clear_model(&s1);
    clear_model(&s2);
    clear_model(&m);

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
