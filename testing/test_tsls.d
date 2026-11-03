/* Smoke tests for tsls.d: two-stage least squares with one endogenous regressor and two instruments (coefficients
   against the by-hand two-step procedure with lsq, residuals from the structural equation, the over-identification
   test), and the list helpers ivreg_process_lists and tsls_make_endolist.

   Not exercised: ivreg_list_omit/add, LIML, GMM options, robust or panel variants, the Hausman test and the multiple-endogenous case.

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;
import core.stdc.stdlib : free;

import gretl.dataset : create_new_dataset;
import gretl.estimate : lsq;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_listh : LISTSEP;
import gretl.gretl_model : clear_model;
import gretl.libgretlh;
import gretl.tsls : ivreg_process_lists, tsls, tsls_make_endolist;

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

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    enum n = 40;
    DatasetRaw * d = create_new_dataset(8, n, 0);

    foreach (i, nm; ["const", "y", "x", "z1", "z2", "xhat", "u", "v"]) {
        d.varname[i][0 .. nm.length] = nm[];
        d.varname[i][nm.length] = 0;
    }
    uint seed = 4242;
    double rnd()
    {
        seed = seed * 1103515245u + 12345u;
        return ((seed >> 8) & 0xffff) / 65536.0 - 0.5;
    }
    foreach (t; 0 .. n) {
        double z1 = rnd() * 4, z2 = rnd() * 4, v = rnd(), u = rnd() + 0.5 * v;

        d.Z[0][t] = 1;
        d.Z[3][t] = z1;
        d.Z[4][t] = z2;
        d.Z[2][t] = 1.0 + z1 + 0.5 * z2 + v;
        d.Z[1][t] = 1.0 + 2.0 * d.Z[2][t] + u;
    }
    d.t1 = 0;
    d.t2 = n - 1;

    /* y const x ; const z1 z2 */
    int[8] list = [7, 1, 0, 2, LISTSEP, 0, 3, 4];
    Model m = tsls(list.ptr, d, Gretlopt.OPT_NONE);

    check(m.errcode == 0 && m.ncoeff == 2 && m.ci == GretlCmdIndex.ivreg, "tsls estimates a model with two coefficients");

    /* by hand: first stage x on const, z1, z2; second stage y on const, xhat */
    int[5] l1 = [4, 2, 0, 3, 4];
    Model s1 = lsq(l1.ptr, d, GretlCmdIndex.ols, Gretlopt.OPT_NONE);
    foreach (t; 0 .. n) {
        d.Z[5][t] = s1.coeff[0] + s1.coeff[1] * d.Z[3][t] + s1.coeff[2] * d.Z[4][t];
    }
    int[4] l2 = [3, 1, 0, 5];
    Model s2 = lsq(l2.ptr, d, GretlCmdIndex.ols, Gretlopt.OPT_NONE);

    check(s1.errcode == 0 && s2.errcode == 0, "by-hand first and second stages estimate");
    check(close(m.coeff[0], s2.coeff[0], 1e-6) && close(m.coeff[1], s2.coeff[1], 1e-6), "2SLS coefficients equal the two-step OLS coefficients");
    check(fabs(m.coeff[1] - 2.0) < 0.5, "slope estimate is near its true value 2");

    /* residuals use the actual (not fitted) x */
    double e0 = d.Z[1][0] - (m.coeff[0] + m.coeff[1] * d.Z[2][0]);
    check(close(m.uhat[0], e0, 1e-6), "residuals come from the structural equation");
    check(m.dfd == n - 2, "degrees of freedom n - 2");

    /* the over-identification test (one extra instrument) is attached to the model */
    check(m.ntests >= 1, "a test (Sargan over-identification) is attached");

    clear_model(&s1);
    clear_model(&s2);
    clear_model(&m);

    /* list helpers */
    int * reglist, instlist;
    int err = ivreg_process_lists(list.ptr, &reglist, &instlist);

    check(err == 0 && reglist !is null && reglist[0] == 3, "ivreg_process_lists splits the regressors");
    check(instlist !is null && instlist[0] == 3, "ivreg_process_lists splits the instruments");

    int * endo = tsls_make_endolist(reglist, instlist, &err);

    check(err == 0 && endo !is null && endo[0] == 1 && endo[1] == 2, "tsls_make_endolist finds x as the only endogenous regressor");
    free(endo);
    free(reglist);
    free(instlist);

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
