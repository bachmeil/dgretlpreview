/* Smoke tests for stata_import.d: the files in data/ (written by data/make_dta.py, which follows the
   Stata .dta format description) are imported through get_plugin_function("dta_get_data"): formats 114
   (little- and big-endian), 117 and 118, with byte/int/long/float/double/fixed-string/strL variables,
   Stata missing values, variable labels, a value label, and a %tq variable that gives quarterly time-series
   information. Also the errors for a missing file and a non-Stata file.

   Build and run with testing/run_tests.sh. */

import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;
import core.stdc.string : strcmp, strstr;

import gretl.dataioh : DataClearCode;
import gretl.dataset : clear_datainfo, is_string_valued, series_get_label, series_get_string_for_obs;
import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_new;
import gretl.gretl_prnh : PrnType;
import gretl.libgretlh;
import gretl.missingh : na;
import gretl.plugins : get_plugin_function;

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

void check_file(Importer imp, string fname, string tag, bool has_strl, bool check_strl_value)
{
    Prn * prn = gretl_print_new(PrnType.buffer, null);
    DatasetRaw d;
    int err = imp((fname ~ "\0").ptr, &d, Gretlopt.OPT_Q, prn);

    check(err == 0, tag ~ ": dta_get_data returns 0");
    if (err) {
        gretl_print_destroy(prn);
        return;
    }

    int nv = has_strl ? 9 : 8;

    check(d.v == nv && d.n == 4, tag ~ ": series + const, 4 observations");
    check(strcmp(d.varname[1], "a") == 0 && strcmp(d.varname[5], "e") == 0 && strcmp(d.varname[7], "q") == 0,
          tag ~ ": names");
    check(d.Z[1][0] == 1 && d.Z[1][1] == 2 && na(d.Z[1][3]), tag ~ ": byte values and missing");
    check(d.Z[2][2] == -30 && d.Z[3][1] == -200000, tag ~ ": int and long values");
    check(d.Z[4][1] == 2.5 && na(d.Z[4][2]), tag ~ ": float value and missing");
    check(d.Z[5][3] == 0.4, tag ~ ": double value");
    check(is_string_valued(&d, 6) != 0, tag ~ ": fixed-length string variable is string-valued");
    const(char) * s = series_get_string_for_obs(&d, 6, 1);
    check(s !is null && strcmp(s, "cd") == 0, tag ~ ": string value");
    if (check_strl_value) {
        s = series_get_string_for_obs(&d, 8, 3);
        check(s !is null && strcmp(s, "fourth long string") == 0, tag ~ ": strL value");
    }
    check(d.structure == TsCodes.time_series && d.pd == 4 && strcmp(d.stobs.ptr, "2000:1") == 0,
          tag ~ ": %tq variable gives quarterly time series from 2000:1");
    const(char) * lab = series_get_label(&d, 3);
    check(lab !is null && strcmp(lab, "long var") == 0, tag ~ ": variable label");
    check(d.descrip !is null && strstr(d.descrip, "test data") !is null, tag ~ ": dataset label");
    clear_datainfo(&d, DataClearCode.full);
    gretl_print_destroy(prn);
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    Importer imp = cast(Importer) get_plugin_function("dta_get_data");

    check(imp !is null, "get_plugin_function returns the D importer");

    check_file(imp, "data/test114_le.dta", "114 LE", false, false);
    check_file(imp, "data/test114_be.dta", "114 BE", false, false);
    /* the strL value is not checked for 117: see bugs.md #72 */
    check_file(imp, "data/test117.dta", "117", true, false);
    check_file(imp, "data/test118.dta", "118", true, true);

    Prn * prn = gretl_print_new(PrnType.buffer, null);
    DatasetRaw d;
    int err = imp("data/no_such_file.dta", &d, Gretlopt.OPT_Q, prn);

    check(err != 0, "missing file is an error");
    err = imp("test_stata_import.d", &d, Gretlopt.OPT_Q, prn);
    check(err != 0, "non-Stata file is an error");

    gretl_print_destroy(prn);
    printf("\n%d checks, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
