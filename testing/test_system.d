/* Smoke tests for system.d: building a two-equation system with equation_system_new/append, the method-name
   helpers, list accessors, sample adjustment, naming and destroy, and estimation through the sysest plugin module: OLS
   equation by equation against lsq, and SUR with identical regressors (which equals OLS).

   Not exercised: restrictions, forecasts, identities, FIML/LIML and the XML/bundle serialisation.

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;
import core.stdc.string : strcmp;

import gretl.dataset : create_new_dataset;
import gretl.estimate : lsq;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_model : clear_model;
import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_new;
import gretl.gretl_prnh : PrnType;
import gretl.libgretlh;
import gretl.system;
import gretl.systemh : EquationSystem, GretlSystemMethods;

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

    check(system_method_from_string("sur") == GretlSystemMethods.sur && system_method_from_string("3sls") == GretlSystemMethods.sls3
          && system_method_from_string("ols") == GretlSystemMethods.ols, "system_method_from_string");
    check(system_method_from_string("nonsense") == GretlSystemMethods.max, "unknown method name gives the maximum code");
    check(strcmp(system_method_short_string(GretlSystemMethods.tsls), "tsls") == 0, "system_method_short_string");
    check(strcmp(system_method_full_string(GretlSystemMethods.sur), "Seemingly Unrelated Regressions") == 0, "system_method_full_string");

    enum n = 40;
    DatasetRaw * d = create_new_dataset(5, n, 0);

    foreach (i, nm; ["const", "y1", "y2", "x", "z"]) {
        d.varname[i][0 .. nm.length] = nm[];
        d.varname[i][nm.length] = 0;
    }
    uint seed = 31337;
    double rnd()
    {
        seed = seed * 1103515245u + 12345u;
        return ((seed >> 8) & 0xffff) / 65536.0 - 0.5;
    }
    foreach (t; 0 .. n) {
        d.Z[0][t] = 1;
        d.Z[3][t] = 3.0 * rnd();
        d.Z[4][t] = 3.0 * rnd();
        d.Z[1][t] = 1.0 + 2.0 * d.Z[3][t] + rnd();
        d.Z[2][t] = -1.0 + 0.5 * d.Z[3][t] + d.Z[4][t] + rnd();
    }
    d.t1 = 0;
    d.t2 = n - 1;

    int err = 0;
    EquationSystem * sys = equation_system_new(GretlSystemMethods.ols, "mysys", &err);

    check(err == 0 && sys !is null, "equation_system_new");

    int[4] l1 = [3, 1, 0, 3];
    int[5] l2 = [4, 2, 0, 3, 4];

    check(equation_system_append(sys, l1.ptr) == 0 && equation_system_append(sys, l2.ptr) == 0, "two equations appended");
    check(sys.neqns == 2 && system_get_method(sys) == GretlSystemMethods.ols, "two equations, method OLS");
    check(system_get_depvar(sys, 0) == 1 && system_get_depvar(sys, 1) == 2, "dependent variables of the equations");

    check(system_get_list_length(sys, 0) == 3 && system_get_list_length(sys, 1) == 4, "regression list lengths");
    int * l = system_get_list(sys, 1);
    check(l !is null && l[0] == 4 && l[1] == 2, "list of equation 2 starts with its dependent variable");
    check(system_n_indep_vars(sys) >= 4, "independent variables counted across equations");

    err = system_adjust_t1t2(sys, d);
    check(err != 0, "system_adjust_t1t2 before the lists have been checked is an error");

    equation_system_set_name(sys, "renamed");
    check(strcmp(sys.name, "renamed") == 0, "equation_system_set_name");

    char * full = system_get_full_string(sys);
    check(full !is null && full[0] != 0, "system_get_full_string describes the system");
    import core.stdc.stdlib : free;
    free(full);

    /* estimation by OLS, equation by equation */
    Prn * prn = gretl_print_new(PrnType.buffer, &err);

    err = equation_system_estimate(sys, d, Gretlopt.OPT_NONE, prn);
    check(err == 0, "equation_system_estimate by OLS");

    Model m1 = lsq(l1.ptr, d, GretlCmdIndex.ols, Gretlopt.OPT_NONE);
    Model m2 = lsq(l2.ptr, d, GretlCmdIndex.ols, Gretlopt.OPT_NONE);
    /* the per-equation models are freed after estimation unless the system is saved; the coefficients stay in b */
    check(sys.b !is null && sys.b.rows == 5, "the stacked coefficient vector has 5 elements");
    if (sys.b !is null && sys.b.rows == 5) {
        const(double) * b = sys.b.val;

        check(close(b[0], m1.coeff[0], 1e-8) && close(b[1], m1.coeff[1], 1e-8), "equation 1 coefficients equal OLS");
        check(close(b[2], m2.coeff[0], 1e-8) && close(b[3], m2.coeff[1], 1e-8) && close(b[4], m2.coeff[2], 1e-8),
              "equation 2 coefficients equal OLS");
    }
    clear_model(&m1);
    clear_model(&m2);
    gretl_print_destroy(prn);
    equation_system_destroy(sys);

    /* SUR with identical regressors in both equations equals OLS */
    EquationSystem * sur = equation_system_new(GretlSystemMethods.sur, "sursys", &err);
    int[4] la = [3, 1, 0, 3];
    int[4] lb = [3, 2, 0, 3];

    equation_system_append(sur, la.ptr);
    equation_system_append(sur, lb.ptr);
    prn = gretl_print_new(PrnType.buffer, &err);
    err = equation_system_estimate(sur, d, Gretlopt.OPT_NONE, prn);
    check(err == 0, "equation_system_estimate by SUR");

    Model oa = lsq(la.ptr, d, GretlCmdIndex.ols, Gretlopt.OPT_NONE);
    Model ob = lsq(lb.ptr, d, GretlCmdIndex.ols, Gretlopt.OPT_NONE);
    if (sur.b !is null && sur.b.rows == 4) {
        const(double) * b = sur.b.val;

        check(close(b[1], oa.coeff[1], 1e-6) && close(b[3], ob.coeff[1], 1e-6), "SUR with identical regressors gives the OLS slopes");
    } else {
        check(false, "SUR coefficient vector has 4 elements");
    }
    clear_model(&oa);
    clear_model(&ob);
    gretl_print_destroy(prn);
    equation_system_destroy(sur);
    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
