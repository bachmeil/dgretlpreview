/* Smoke tests for tsjoin.d: the calendar helpers (q2m, m2q, mofq), use_tsjoin, and fill_ts_joiner for
   monthly <- annual and quarterly <- annual and monthly <- quarterly joins, and tj_continue.

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import gretl.dataset : create_new_dataset;
import gretl.libgretlh;
import gretl.tsjoin;

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

DatasetRaw * make_ts(int n, int pd, string start)
{
    import gretl.dataio : get_date_x;

    DatasetRaw * d = create_new_dataset(1, n, 0);

    d.pd = pd;
    d.structure = TsCodes.time_series;
    d.stobs[0 .. start.length] = start[];
    d.stobs[start.length] = 0;
    d.sd0 = get_date_x(pd, d.stobs.ptr);
    d.t1 = 0;
    d.t2 = n - 1;
    return d;
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    check(q2m(1) == 1 && q2m(2) == 4 && q2m(4) == 10, "q2m");
    check(m2q(1) == 1 && m2q(3) == 1 && m2q(4) == 2 && m2q(12) == 4, "m2q");
    check(mofq(1) == 1 && mofq(2) == 2 && mofq(3) == 3 && mofq(4) == 1 && mofq(12) == 3, "mofq");

    DatasetRaw * m24 = make_ts(24, 12, "2020:01");
    DatasetRaw * q8 = make_ts(8, 4, "2020:1");
    DatasetRaw * a2 = make_ts(2, 1, "2020");

    check(use_tsjoin(m24, a2) == 1 && use_tsjoin(m24, q8) == 1 && use_tsjoin(q8, a2) == 1, "use_tsjoin: lower frequency on the right");
    check(use_tsjoin(a2, m24) == 0 && use_tsjoin(q8, m24) == 0 && use_tsjoin(m24, m24) == 0, "use_tsjoin: otherwise not");

    TsJoiner tj;

    check(fill_ts_joiner(m24, a2, &tj) == 0, "fill_ts_joiner monthly <- annual");
    check(tj.pdr == 12 && tj.t1 == 0 && tj.rt1 == 0 && tj.t2 == 23 && tj.rminor == 1, "ratio 12, whole sample, minor index 1");

    m24.t1 = 3;   /* start in 2020:04 */
    check(fill_ts_joiner(m24, a2, &tj) == 0 && tj.t1 == 3 && tj.rt1 == 0 && tj.rminor == 4, "start moved to 2020:04 gives minor index 4");
    m24.t1 = 0;

    check(fill_ts_joiner(q8, a2, &tj) == 0, "fill_ts_joiner quarterly <- annual");
    check(tj.pdr == 4 && tj.t1 == 0 && tj.t2 == 7 && tj.rminor == 1, "ratio 4, whole sample");

    check(fill_ts_joiner(m24, q8, &tj) == 0, "fill_ts_joiner monthly <- quarterly");
    check(tj.pdr == 3 && tj.t1 == 0 && tj.rt1 == 0 && tj.t2 == 23 && tj.rminor == 1, "ratio 3, whole sample");

    int s = 0;
    tj.pdr = 3;
    check(tj_continue(&tj, 1, &s) == 2 && s == 0, "tj_continue within a cycle");
    check(tj_continue(&tj, 3, &s) == 1 && s == 1, "tj_continue wraps and advances the source index");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
