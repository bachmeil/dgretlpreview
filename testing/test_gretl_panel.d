/* Smoke tests for gretl_panel.d: panel structure from unit/time variables, panel_isconst, the orientation switch,
   and fixed- and random-effects, between and pooled estimation on a small balanced panel (slopes compared with
   the within and between estimators computed here by hand).

   Not exercised: unbalanced panels, clustered and Driscoll-Kraay covariance matrices, WLS/ML groupwise models, the
   Hausman and Breusch-Pagan output, time dummies and the panel-dataset utilities beyond those above.

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import gretl.dataset : create_new_dataset, destroy_dataset;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_model : clear_model;
import gretl.gretl_panel;
import gretl.libgretlh;

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

enum N = 5;
enum T = 4;

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    /* series: 0 const, 1 y, 2 x, 3 unit, 4 time; stacked time series */
    DatasetRaw * d = create_new_dataset(5, N * T, 0);

    foreach (i, nm; ["const", "y", "x", "unit", "time"]) {
        d.varname[i][0 .. nm.length] = nm[];
        d.varname[i][nm.length] = 0;
    }
    immutable double[N * T] noise = [0.3, -0.2, 0.5, -0.6, 0.1, 0.4, -0.3, -0.2, -0.5, 0.2, 0.6, -0.3, 0.0, 0.2, -0.1, 0.4, -0.4, 0.3, 0.1, -0.1];
    foreach (i; 0 .. N) {
        foreach (t; 0 .. T) {
            int s = i * T + t;
            double x = (i * 1.7 + t * 2.3 + (i * t) % 3) % 7;

            d.Z[0][s] = 1;
            d.Z[2][s] = x;
            d.Z[1][s] = 3.0 * i + 2.0 * x + noise[s];
            d.Z[3][s] = i + 1;
            d.Z[4][s] = 2000 + t;
        }
    }
    d.pd = 1;
    d.structure = TsCodes.cross_section;

    check(set_panel_structure_from_vars(3, 4, d) == 0, "set_panel_structure_from_vars");
    check(d.structure == TsCodes.stacked_time_series && d.pd == T, "the dataset is a stacked time-series panel with T periods");
    check(panel_isconst(0, N * T - 1, T, d.Z[3], 0) == 1 && panel_isconst(0, N * T - 1, T, d.Z[3], 1) == 0 && panel_isconst(0, N * T - 1, T, d.Z[4], 1) == 1, "panel_isconst: the unit variable is constant within units but varies across them");

    /* hand-computed within and between slopes */
    double[N] ybar = 0, xbar = 0;
    double sxy = 0, sxx = 0;

    foreach (i; 0 .. N) {
        foreach (t; 0 .. T) { ybar[i] += d.Z[1][i * T + t] / T; xbar[i] += d.Z[2][i * T + t] / T; }
    }
    foreach (i; 0 .. N) {
        foreach (t; 0 .. T) {
            double dx = d.Z[2][i * T + t] - xbar[i];
            double dy = d.Z[1][i * T + t] - ybar[i];

            sxy += dx * dy;
            sxx += dx * dx;
        }
    }
    double bwithin = sxy / sxx;
    double mx = 0, my = 0, bxy = 0, bxx = 0;

    foreach (i; 0 .. N) { mx += xbar[i] / N; my += ybar[i] / N; }
    foreach (i; 0 .. N) { bxy += (xbar[i] - mx) * (ybar[i] - my); bxx += (xbar[i] - mx) * (xbar[i] - mx); }
    double bbetween = bxy / bxx;

    int[4] list = [3, 1, 0, 2];
    Model m = real_panel_model(list.ptr, d, Gretlopt.OPT_F, null);

    check(m.errcode == 0, "fixed-effects model runs");
    check(m.ncoeff >= 1 && close(m.coeff[m.ncoeff - 1], bwithin, 1e-6), "fixed-effects slope equals the within estimator");
    clear_model(&m);

    m = real_panel_model(list.ptr, d, Gretlopt.OPT_B, null);
    check(m.errcode == 0 && close(m.coeff[1], bbetween, 1e-6), "between slope equals the group-means regression slope");
    clear_model(&m);

    m = real_panel_model(list.ptr, d, Gretlopt.OPT_U, null);
    check(m.errcode == 0 && m.coeff[1] >= (bwithin < bbetween ? bwithin : bbetween) - 1e-6 && m.coeff[1] <= (bwithin > bbetween ? bwithin : bbetween) + 1e-6, "random-effects slope lies between the within and between slopes");
    clear_model(&m);

    m = real_panel_model(list.ptr, d, Gretlopt.OPT_NONE, null);
    check(m.errcode == 0 && m.ncoeff == 2, "pooled OLS through the panel interface");
    clear_model(&m);

    /* orientation switch: stacked time series -> stacked cross sections and back */
    double first = d.Z[1][1];

    check(switch_panel_orientation(d) == 0 && d.structure == TsCodes.stacked_cross_section, "switch_panel_orientation");
    check(d.Z[1][N] == first, "the second observation of unit 1 moves to the second block");
    check(switch_panel_orientation(d) == 0 && d.Z[1][1] == first, "switching back restores the order");
    check(panel_padding_rows(d) == 0, "no padding rows in a balanced panel");

    destroy_dataset(d);
    printf("\n%d checks, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
