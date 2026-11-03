/* Smoke tests for adf_kpss.d: each function is called on small
   simulated series, and we check only that the call returns and
   produces output (not that the numbers are right).

   Build and run with testing/run_tests.sh. */

import core.stdc.stdio : _IOLBF, fflush, printf, setvbuf, stdout;
import core.stdc.stdlib : calloc;
import core.stdc.string : strcpy, strlen;

import gretl.adf_kpss;
import gretl.dataset : create_new_dataset;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_list : gretl_list_new;
import gretl.gretl_matrixh;
import gretl.gretl_matrix : gretl_matrix_free;
import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_get_buffer, gretl_print_new;
import gretl.gretl_prnh : PrnType;
import gretl.gretl_utils : get_last_pvalue, get_last_test_statistic;
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

/* deterministic pseudo-random numbers in (-0.5, 0.5) */

ulong lcg_state = 12345;

double urand()
{
    lcg_state = lcg_state * 6364136223846793005UL + 1442695040888963407UL;
    return ((lcg_state >> 11) / cast(double) (1UL << 53)) - 0.5;
}

/* series: 0 const, 1 random walk, 2 stationary AR(1), 3 cointegrated with 1, 4 trend */
enum RW = 1, AR1 = 2, COINT = 3;

DatasetRaw * make_dataset(int T)
{
    DatasetRaw * dset = create_new_dataset(5, T, 0);
    string[] names = ["const", "rw", "ar1", "rw2", "trend"];

    foreach (i, nm; names) {
        strcpy(dset.varname[i], (nm ~ "\0").ptr);
    }
    dset.pd = 4;
    dset.structure = TsCodes.time_series;
    dset.sd0 = 1990.1;
    strcpy(dset.stobs.ptr, "1990:1");
    dset.t1 = 0;
    dset.t2 = T - 1;

    double rw = 0, ar = 0;

    foreach (t; 0 .. T) {
        rw += urand();
        ar = 0.5 * ar + urand();
        dset.Z[1][t] = rw;
        dset.Z[2][t] = ar;
        dset.Z[3][t] = 2.0 * rw + 0.3 * urand();
        dset.Z[4][t] = t + 1;
    }

    return dset;
}

int * list1(int v)
{
    int * l = gretl_list_new(1);

    l[1] = v;
    return l;
}

int * list2(int v1, int v2)
{
    int * l = gretl_list_new(2);

    l[1] = v1;
    l[2] = v2;
    return l;
}

void run_test(string name, int delegate(Prn *) fn)
{
    int err;
    Prn * prn = gretl_print_new(PrnType.buffer, &err);
    int ret = fn(prn);
    const(char) * buf = gretl_print_get_buffer(prn);

    printf("---- %.*s (return value %d) ----\n%s\n", cast(int) name.length, name.ptr, ret, buf);
    check(ret == 0 && strlen(buf) > 0, name ~ ": returns 0 and prints something");
    gretl_print_destroy(prn);
}

void main()
{
    setvbuf(stdout, null, _IOLBF, 4096);   /* keep output ordered if a stub aborts */
    DatasetRaw * dset = make_dataset(100);
    int err;

    /* p-value lookup (no dataset needed) */
    foreach (itv; [1, 2, 3, 4]) {
        double p = get_urc_pvalue(-2.5, 100, 1, itv);

        printf("get_urc_pvalue(-2.5, 100, 1, %d) = %g\n", itv, p);
        /* the p-value comes from the (not yet available) mackinnon plugin,
           so NA is the expected result here: we only check that it returns */
        check(!na(p) && p >= 0 && p <= 1, "get_urc_pvalue returns a probability");
    }

    /* KPSS critical values */
    GretlMatrixRaw * kc = kpss_critvals(100, 0, &err);

    check(err == 0 && kc !is null, "kpss_critvals (no trend)");
    gretl_matrix_free(kc);
    kc = kpss_critvals(100, 1, &err);
    check(err == 0 && kc !is null, "kpss_critvals (trend)");
    gretl_matrix_free(kc);

    /* ADF, with the various deterministic-term options */
    run_test("adf_test: default (constant)", delegate(Prn * prn) {
        return adf_test(2, list1(RW), dset, Gretlopt.OPT_C, prn);
    });
    run_test("adf_test: no constant", delegate(Prn * prn) {
        return adf_test(2, list1(RW), dset, Gretlopt.OPT_N, prn);
    });
    run_test("adf_test: constant and trend", delegate(Prn * prn) {
        return adf_test(1, list1(RW), dset, Gretlopt.OPT_T, prn);
    });
    run_test("adf_test: stationary series", delegate(Prn * prn) {
        return adf_test(2, list1(AR1), dset, Gretlopt.OPT_C, prn);
    });
    run_test("adf_test: verbose, test-down", delegate(Prn * prn) {
        return adf_test(4, list1(RW), dset, Gretlopt.OPT_C | Gretlopt.OPT_V | Gretlopt.OPT_E, prn);
    });
    printf("last ADF statistic %g, p-value %g\n", get_last_test_statistic(), get_last_pvalue());

    /* KPSS */
    run_test("kpss_test: level", delegate(Prn * prn) {
        return kpss_test(4, list1(RW), dset, Gretlopt.OPT_NONE, prn);
    });
    run_test("kpss_test: trend", delegate(Prn * prn) {
        return kpss_test(4, list1(AR1), dset, Gretlopt.OPT_T, prn);
    });
    run_test("kpss_test: verbose", delegate(Prn * prn) {
        return kpss_test(2, list1(RW), dset, Gretlopt.OPT_V, prn);
    });

    /* Engle-Granger */
    run_test("engle_granger_test: cointegrated pair", delegate(Prn * prn) {
        return engle_granger_test(2, list2(COINT, RW), dset, Gretlopt.OPT_NONE, prn);
    });

    printf("\n%d checks, %d failed\n", n_tests, n_failed);
    fflush(stdout);
    if (n_failed) {
        throw new Exception("some tests failed");
    }
}
