/* Smoke tests for stata_export.d: a dataset is written with get_plugin_function("stata_export") and read
   back with dta_get_data. Covers the storage types chosen by make_types_array() (byte for a dummy, int,
   long, double), missing values, a string-valued series (fixed-length string type), variable labels, the
   --lcnames option, a list of series, and the added time variable for quarterly time series.

   Build and run with testing/run_tests.sh. */

import core.stdc.stdio : _IOLBF, printf, remove, setvbuf, stdout;
import core.stdc.string : strcmp;

import gretl.dataioh : DataClearCode;
import gretl.dataset : clear_datainfo, create_new_dataset, destroy_dataset, is_string_valued, series_get_label,
    series_get_string_for_obs, series_set_label, series_set_string_val;
import gretl.gretl_string_table : series_table_new;
import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_new;
import gretl.gretl_prnh : PrnType;
import gretl.libgretlh;
import gretl.missingh : na;
import gretl.plugins : get_plugin_function;

alias Exporter = extern (C) int function(const(char) *, const(int) *, Gretlopt, const(DatasetRaw) *);
alias Importer = extern (C) int function(const(char) *, DatasetRaw *, Gretlopt, Prn *);

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

void setname(DatasetRaw * d, int i, string name)
{
    d.varname[i][0 .. name.length] = name[];
    d.varname[i][name.length] = 0;
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    Exporter exp = cast(Exporter) get_plugin_function("stata_export");
    Importer imp = cast(Importer) get_plugin_function("dta_get_data");

    check(exp !is null && imp !is null, "get_plugin_function returns the D exporter and importer");

    enum N = 4;
    DatasetRaw * d = create_new_dataset(7, N, 0);

    setname(d, 1, "Dummy");
    setname(d, 2, "small");
    setname(d, 3, "big");
    setname(d, 4, "frac");
    setname(d, 5, "word");
    setname(d, 6, "huge");
    double[N] dummy = [0, 1, 1, 0];
    double[N] small = [1, -250, 9000, double.nan];
    double[N] big = [100000, -2000000, 3, 4];
    double[N] frac = [1.5, double.nan, 3.25, -4];
    double[N] huge = [1e10, 2, 3, 4];
    foreach (t; 0 .. N) {
        d.Z[1][t] = dummy[t];
        d.Z[2][t] = small[t];
        d.Z[3][t] = big[t];
        d.Z[4][t] = frac[t];
        d.Z[6][t] = huge[t];
    }
    string[N] words = ["alpha", "beta", "alpha", "gamma"];
    int terr;

    d.varinfo[5].st = series_table_new(null, 0, &terr);
    foreach (t; 0 .. N) {
        series_set_string_val(d, 5, t, (words[t] ~ "\0").ptr);
    }
    series_set_label(d, 3, "a big series");

    int err = exp("build/export_test.dta", null, cast(Gretlopt) 0, d);

    check(err == 0, "stata_export returns 0");

    Prn * prn = gretl_print_new(PrnType.buffer, null);
    DatasetRaw r;

    err = imp("build/export_test.dta", &r, Gretlopt.OPT_Q, prn);
    check(err == 0, "the file reads back");
    if (!err) {
        check(r.v == 7 && r.n == N, "dimensions");
        check(strcmp(r.varname[1], "Dummy") == 0 && strcmp(r.varname[5], "word") == 0, "names");
        check(r.Z[1][1] == 1 && r.Z[1][3] == 0, "dummy (byte) values");
        check(r.Z[2][1] == -250 && r.Z[2][2] == 9000 && na(r.Z[2][3]), "small integer (int) values and a missing");
        check(r.Z[3][1] == -2000000 && r.Z[3][0] == 100000, "larger integer (long) values");
        check(r.Z[4][0] == 1.5 && na(r.Z[4][1]) && r.Z[4][3] == -4, "non-integer (double) values and a missing");
        check(r.Z[6][0] == 1e10, "very large integers as double");
        check(is_string_valued(&r, 5) != 0, "string-valued series comes back string-valued");
        const(char) * s = series_get_string_for_obs(&r, 5, 3);
        check(s !is null && strcmp(s, "gamma") == 0, "string value");
        const(char) * lab = series_get_label(&r, 3);
        check(lab !is null && strcmp(lab, "a big series") == 0, "variable label");
        clear_datainfo(&r, DataClearCode.full);
    }

    /* lower-case names and a list of series */
    int[3] list = [2, 1, 4];
    DatasetRaw r2;

    err = exp("build/export_test.dta", list.ptr, Gretlopt.OPT_L, d);
    check(err == 0, "export of a list with --lcnames");
    err = imp("build/export_test.dta", &r2, Gretlopt.OPT_Q, prn);
    check(err == 0 && r2.v == 3 && strcmp(r2.varname[1], "dummy") == 0 && strcmp(r2.varname[2], "frac") == 0,
          "only the listed series, in lower case");
    if (!err) clear_datainfo(&r2, DataClearCode.full);

    /* quarterly time series: a time variable is added */
    d.pd = 4;
    d.structure = TsCodes.time_series;
    d.sd0 = 2000.1;
    d.stobs[0 .. 7] = "2000:1\0";
    d.t1 = 0;
    d.t2 = N - 1;
    DatasetRaw r3;

    err = exp("build/export_test.dta", null, cast(Gretlopt) 0, d);
    check(err == 0, "export of a quarterly time series");
    err = imp("build/export_test.dta", &r3, Gretlopt.OPT_Q, prn);
    check(err == 0, "time-series file reads back");
    if (!err) {
        check(r3.v == 8 && strcmp(r3.varname[1], "quarter") == 0, "time variable 'quarter' added");
        check(r3.Z[1][0] == 160 && r3.Z[1][3] == 163, "quarters since 1960q1");
        check(r3.structure == TsCodes.time_series && r3.pd == 4 && strcmp(r3.stobs.ptr, "2000:1") == 0,
              "time-series structure recovered from the %tq format");
        clear_datainfo(&r3, DataClearCode.full);
    }

    err = exp("no_such_dir/x.dta", null, cast(Gretlopt) 0, d);
    check(err != 0, "unwritable path is an error");

    remove("build/export_test.dta");
    destroy_dataset(d);
    gretl_print_destroy(prn);
    printf("\n%d checks, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
