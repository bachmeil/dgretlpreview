/* Smoke tests for subsample.d: restricting the sample by a mask, by a dummy variable, by dropping missing
   values and at random, restoring the full sample, the full-dataset bookkeeping, set_sample on a cross-section
   and the sub-sample mask helpers.

   Not exercised: boolean-expression restrictions (need the genr parser), panel samples, permanent sub-samples
   and string-valued series recoding, sample setting on time-series and daily data.

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;
import core.stdc.stdlib : free, malloc;

import gretl.dataset : create_new_dataset;
import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_new;
import gretl.gretl_prnh : PrnType;
import gretl.libgretlh;
import gretl.missingh : NADBL;
import gretl.subsample;
import gretl.subsampleh : SUBMASK_SENTINEL;

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

DatasetRaw * make_data()
{
    enum n = 10;
    DatasetRaw * d = create_new_dataset(4, n, 0);

    foreach (i, nm; ["const", "x", "dum", "m"]) {
        d.varname[i][0 .. nm.length] = nm[];
        d.varname[i][nm.length] = 0;
    }
    foreach (t; 0 .. n) {
        d.Z[0][t] = 1;
        d.Z[1][t] = 10.0 * (t + 1);
        d.Z[2][t] = (t >= 5)? 1 : 0;
        d.Z[3][t] = (t == 2 || t == 7)? NADBL : t;
    }
    d.t1 = 0;
    d.t2 = n - 1;
    return d;
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    int err = 0, nd = 0;
    Prn * prn = gretl_print_new(PrnType.buffer, &err);
    DatasetRaw * d = make_data();

    check(!dataset_is_subsampled(d) && complex_subsampled() == 0, "a fresh dataset is not sub-sampled");

    /* mask */
    char * mask = cast(char *) malloc(11);   /* gretl masks end with a sentinel byte */
    foreach (t; 0 .. 10) {
        mask[t] = (t >= 2 && t <= 6) ? 1 : 0;
    }
    mask[10] = SUBMASK_SENTINEL;
    err = restrict_sample_from_mask(mask, d, Gretlopt.OPT_NONE);
    check(err == 0 && d.n == 5, "mask keeps observations 2..6");
    check(dataset_is_subsampled(d) && get_full_length_n() == 10, "sub-sampled; full length is 10");
    check(close(d.Z[1][0], 30.0) && close(d.Z[1][4], 70.0), "series values are those of the kept observations");
    DatasetRaw * full = fetch_full_dataset();
    check(full !is null && full.n == 10, "fetch_full_dataset gives the full dataset");

    char * m2 = copy_dataset_submask(d, &err);
    check(err == 0 && m2 !is null && submask_cmp(m2, d.submask) == 0, "copy_dataset_submask matches the dataset's mask");
    free(m2);

    err = restore_full_sample(d, null);
    check(err == 0 && d.n == 10 && !dataset_is_subsampled(d), "restore_full_sample");
    check(close(d.Z[1][9], 100.0), "full data are back");
    free(mask);

    /* dummy */
    err = restrict_sample("dum", null, d, null, Gretlopt.OPT_O, prn, &nd);
    check(err == 0 && nd == 5 && dataset_is_subsampled(d) && d.t1 == 5 && d.t2 == 9, "dummy restriction keeps the contiguous range 5..9 (dum = 1)");
    check(close(d.Z[1][d.t1], 60.0), "first kept observation is the sixth");
    restore_full_sample(d, null);

    /* drop missing */
    int[3] list = [2, 1, 3];
    err = restrict_sample(null, list.ptr, d, null, Gretlopt.OPT_M, prn, &nd);
    check(err == 0 && d.n == 8, "no-missing restriction on x and m drops the two incomplete observations");
    restore_full_sample(d, null);

    /* random */
    err = restrict_sample("4", null, d, null, Gretlopt.OPT_N, prn, &nd);
    check(err == 0 && d.n == 4, "random restriction keeps four observations");
    restore_full_sample(d, null);

    /* incompatible options */
    err = restrict_sample("dum", null, d, null, Gretlopt.OPT_O | Gretlopt.OPT_M, prn, &nd);
    check(err != 0, "incompatible options are rejected");

    /* set_sample on a cross-section */
    err = set_sample("3", "8", d, Gretlopt.OPT_NONE);
    check(err == 0 && d.t1 == 2 && d.t2 == 7, "set_sample 3 8 gives t1 = 2, t2 = 7");
    err = set_sample("1", "99", d, Gretlopt.OPT_NONE);
    check(err != 0, "set_sample beyond the data is an error");

    gretl_print_destroy(prn);
    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
