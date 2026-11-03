/* Smoke tests for boxplots.d: the quartiles, whiskers and outlier detection (quartiles_etc, with the tiny-sample
   special cases), the bootstrap confidence interval for the median, the printed summaries (six_numbers,
   five_numbers_with_interval) and the flag helpers.

   Not exercised: writing a complete gnuplot boxplot (needs the gnuplot support of graphing.c), factorized and
   panel boxplots, the numerical summary read back from a boxplot file and the tic parsing.

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;
import core.stdc.string : strstr;

import gretl.boxplots;
import gretl.gretl_matrix : gretl_matrix_free;
import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_get_buffer, gretl_print_new;
import gretl.gretl_prnh : PrnType;
import gretl.random : gretl_rand_set_seed;
import gretl.libgretlh;
import gretl.missingh : NADBL, na;

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
    gretl_rand_set_seed(20260101);

    /* a sample with one outlier */
    double[12] x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 100];
    Boxplot bp;

    boxplot_init(&bp);
    quartiles_etc(x.ptr, 12, &bp, 1.5);
    check(close(bp.mean, 166.0 / 12.0, 1e-12), "mean");
    check(bp.min == 1 && bp.max == 100, "data minimum and maximum");
    check(!na(bp.lq) && !na(bp.uq) && bp.lq < bp.median && bp.median < bp.uq, "quartiles are ordered around the median");
    check(close(bp.median, 6.5, 1e-9), "median is 6.5");
    check(bp.wmin == 1 && bp.wmax == 11, "the whiskers stop at the last non-outlying values");
    check(bp.outliers !is null && bp.outliers.rows * bp.outliers.cols == 1 && bp.outliers.val[0] == 100.0, "100 is the only outlier");
    gretl_matrix_free(bp.outliers);

    /* without an outlier limit the whiskers reach the extremes */
    Boxplot b2;

    boxplot_init(&b2);
    quartiles_etc(x.ptr, 12, &b2, 0.0);
    check(b2.wmin == 1 && b2.wmax == 100 && b2.outliers is null, "limit 0: whiskers at the extremes, no outliers");

    /* tiny samples */
    Boxplot t1, t2, t3;

    boxplot_init(&t1);
    boxplot_init(&t2);
    boxplot_init(&t3);
    quartiles_etc(x.ptr, 1, &t1, 1.5);
    check(t1.min == 1 && t1.max == 1 && na(t1.median) && na(t1.lq), "one observation");
    quartiles_etc(x.ptr, 2, &t2, 1.5);
    check(t2.mean == 1.5 && na(t2.median), "two observations: only the mean");
    quartiles_etc(x.ptr, 3, &t3, 1.5);
    check(t3.median == 2.0 && na(t3.lq) && close(t3.mean, 2.0, 1e-12), "three observations: median and mean, no quartiles");
    gretl_matrix_free(t1.outliers);
    gretl_matrix_free(t2.outliers);
    gretl_matrix_free(t3.outliers);

    /* bootstrap interval for the median */
    double[11] y = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11];
    double lo, hi;

    check(median_interval(y.ptr, 11, &lo, &hi) == 0 && lo <= 6.0 && 6.0 <= hi && lo >= 1.0 && hi <= 11.0, "the bootstrap interval brackets the median");

    /* missing-value markers and printing */
    check(bpna(NADBL) && bpna(double.infinity) && !bpna(2.5), "bpna");
    int err = 0;
    Prn * prn = gretl_print_new(PrnType.buffer, &err);
    bp.outliers = null;
    six_numbers(&bp, 1, prn);
    const(char) * buf = gretl_print_get_buffer(prn);
    check(buf !is null && strstr(buf, "100") !is null && strstr(buf, "6.5") !is null, "six_numbers prints the extremes and the median");
    gretl_print_destroy(prn);

    /* flags */
    PlotGroup g;

    g.flags = 0;
    set_show_mean(&g);
    set_do_intervals(&g);
    check(show_mean(&g) && do_intervals(&g) && !use_xtics(&g) && !factorized(&g), "flag setters and getters");
    unset_show_mean(&g);
    check(!show_mean(&g) && do_intervals(&g), "unset_show_mean");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
