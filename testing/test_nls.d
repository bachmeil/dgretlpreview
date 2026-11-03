/* Smoke tests for nls.d: creating and destroying a nonlinear-model specification, adding a parameter list, the
   coefficient updating, setting the sample range, the Gauss-Newton regression (against plain OLS) and the default
   tolerance.

   Not exercised: NLS/MLE/GMM estimation proper and parsing of "nls"/"mle"/"deriv" blocks, which evaluate genr
   expressions and so need the parser (not ported), and the bootstrap helpers.

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import gretl.dataset : create_new_dataset;
import gretl.estimate : lsq;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_model : clear_model;
import gretl.libgretlh;
import gretl.nls;
import gretl.uservar : gretl_scalar_get_value;
import gretl.nlspech : Nlspec;

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

    enum n = 20;
    DatasetRaw * d = create_new_dataset(3, n, 0);

    foreach (i, nm; ["const", "y", "x"]) {
        d.varname[i][0 .. nm.length] = nm[];
        d.varname[i][nm.length] = 0;
    }
    uint seed = 99;
    double rnd()
    {
        seed = seed * 1103515245u + 12345u;
        return ((seed >> 8) & 0xffff) / 65536.0 - 0.5;
    }
    foreach (t; 0 .. n) {
        d.Z[0][t] = 1;
        d.Z[2][t] = 3.0 * rnd();
        d.Z[1][t] = 1.0 + 2.0 * d.Z[2][t] + rnd();
    }
    d.t1 = 0;
    d.t2 = n - 1;

    double tol = get_default_nls_toler();

    check(tol > 0 && tol < 1, "default NLS tolerance is a small positive number");

    Nlspec * spec = nlspec_new(GretlCmdIndex.nls, d);

    check(spec !is null && spec.ci == GretlCmdIndex.nls && spec.nparam == 0 && spec.t1 == 0 && spec.t2 == n - 1, "nlspec_new takes the sample range from the dataset");

    nlspec_set_t1_t2(spec, 2, 17);
    check(spec.t1 == 2 && spec.t2 == 17, "nlspec_set_t1_t2");

    double[2] vals = [1.5, 2.5];
    const(char) *[2] names = ["alpha".ptr, "beta".ptr];
    int err = nlspec_add_param_list(spec, 2, vals.ptr, cast(char **) names.ptr);

    check(err == 0 && spec.nparam == 2 && spec.ncoeff == 2, "nlspec_add_param_list adds two scalar parameters");
    check(spec.params !is null && spec.params[0].name[0] == 'a' && spec.params[1].name[0] == 'b', "parameter names are kept");
    check(spec.coeff !is null && spec.coeff[0] == 1.5 && spec.coeff[1] == 2.5, "initial values are stored");

    double[2] nb = [3.0, 4.0];
    err = update_coeff_values(nb.ptr, spec);
    int e1 = 0, e2 = 0;
    check(err == 0 && gretl_scalar_get_value("alpha", &e1) == 3.0 && gretl_scalar_get_value("beta", &e2) == 4.0, "update_coeff_values sets the parameter scalars");

    nlspec_destroy(spec);

    /* Gauss-Newton regression is OLS on the Jacobian columns */
    int[4] glist = [3, 1, 0, 2];
    Model g = GNR(glist.ptr, d, Gretlopt.OPT_NONE, null);
    Model o = lsq(glist.ptr, d, GretlCmdIndex.ols, Gretlopt.OPT_A);

    check(g.errcode == 0 && o.errcode == 0, "GNR and OLS estimate");
    check(close(g.coeff[0], o.coeff[0], 1e-9) && close(g.coeff[1], o.coeff[1], 1e-9), "GNR coefficients equal OLS");
    clear_model(&g);
    clear_model(&o);

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
