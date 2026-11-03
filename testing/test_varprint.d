/* Smoke tests for varprint.d on a simulated bivariate VAR(1): the full text printout of the VAR (via gretl_VAR and via
   gretl_VAR_print), the contemporaneous covariance matrix, the impulse responses and the forecast-error variance
   decompositions.

   Not exercised: VECM printing, the TeX/RTF/CSV variants and the numerical values (only the presence of expected
   headings is checked).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs, sqrt;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import gretl.dataset : create_new_dataset;
import gretl.estimate : lsq;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_matrixh : GretlMatrixRaw;
import gretl.gretl_model : clear_model;
import gretl.libgretlh;
import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_get_buffer, gretl_print_new;
import gretl.gretl_prnh : PrnType;
import core.stdc.string : strlen, strstr;
import gretl.var;
import gretl.varprint;
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

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    enum T = 80;
    DatasetRaw * d = create_new_dataset(6, T, 0);

    d.pd = 1;
    d.structure = TsCodes.time_series;
    foreach (i, nm; ["const", "y1", "y2", "y1_1", "y2_1", "junk"]) {
        d.varname[i][0 .. nm.length] = nm[];
        d.varname[i][nm.length] = 0;
    }

    uint seed = 12345;
    double rnd()
    {
        seed = seed * 1103515245u + 12345u;
        return ((seed >> 8) & 0xffff) / 65536.0 - 0.5;
    }

    d.Z[1][0] = 0;
    d.Z[2][0] = 0;
    foreach (t; 0 .. T) {
        d.Z[0][t] = 1;
        d.Z[5][t] = 0;
        if (t > 0) {
            d.Z[1][t] = 0.5 * d.Z[1][t-1] + 0.2 * d.Z[2][t-1] + rnd();
            d.Z[2][t] = 0.1 * d.Z[1][t-1] + 0.4 * d.Z[2][t-1] + rnd();
            d.Z[3][t] = d.Z[1][t-1];
            d.Z[4][t] = d.Z[2][t-1];
        }
    }
    d.t1 = 1;
    d.t2 = T - 1;

    int[3] ylist = [2, 1, 2];
    int err = 0;
    Prn * prn = gretl_print_new(PrnType.buffer, &err);

    GretlVar * v = gretl_VAR(1, null, ylist.ptr, d, Gretlopt.OPT_NONE, prn, &err);
    const(char) * buf = gretl_print_get_buffer(prn);

    check(err == 0 && v !is null, "gretl_VAR(1) estimates and prints");
    check(has(buf, "VAR system"), "heading 'VAR system'");
    check(has(buf, "Equation 1") && has(buf, "Equation 2"), "both equations are printed");
    check(has(buf, "y1_1") || has(buf, "y1(-1)") || has(buf, "y1_"), "lagged variable names appear");
    check(has(buf, "Akaike") || has(buf, "AIC"), "information criteria are printed");
    gretl_print_destroy(prn);

    prn = gretl_print_new(PrnType.buffer, &err);
    err = gretl_VAR_print_sigma(v, prn);
    buf = gretl_print_get_buffer(prn);
    check(err == 0 && buf !is null && strlen(buf) > 0, "gretl_VAR_print_sigma prints the covariance matrix");
    gretl_print_destroy(prn);

    prn = gretl_print_new(PrnType.buffer, &err);
    err = gretl_VAR_print_all_impulse_responses(v, d, 4, prn);
    buf = gretl_print_get_buffer(prn);
    check(err == 0 && has(buf, "esponse"), "impulse responses are printed");
    gretl_print_destroy(prn);

    prn = gretl_print_new(PrnType.buffer, &err);
    err = gretl_VAR_print_all_fcast_decomps(v, d, 4, prn);
    buf = gretl_print_get_buffer(prn);
    check(err == 0 && has(buf, "ecomposition"), "forecast-error variance decompositions are printed");
    gretl_print_destroy(prn);

    gretl_VAR_free(v);
    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}

bool has(const(char) * buf, string needle)
{
    import std.string : toStringz;
    return buf !is null && strstr(buf, toStringz(needle)) !is null;
}
