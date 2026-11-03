/* Smoke tests for excel_import.d, workbook.d, ms_ole.d and biffh.d: data/test_import.xls (written by
   LibreOffice: sheets "Data" (numbers with one blank), "Strs" (a string column and a number column)
   and "Dates" (quarterly dates)) is imported through get_plugin_function("xls_get_data").

   Build and run with testing/run_tests.sh. */

import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;
import core.stdc.string : strcmp;

import gretl.dataioh : DataClearCode;
import gretl.dataset : clear_datainfo, is_string_valued, series_get_string_for_obs;
import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_new;
import gretl.gretl_prnh : PrnType;
import gretl.libgretlh;
import gretl.missingh : na;
import gretl.plugins : get_plugin_function;

alias Importer = extern (C) int function(const(char) *, int *, char *, DatasetRaw *, Gretlopt, Prn *, void *);

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

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    Importer imp = cast(Importer) get_plugin_function("xls_get_data");

    check(imp !is null, "get_plugin_function returns the D importer");

    Prn * prn = gretl_print_new(PrnType.buffer, null);
    DatasetRaw d;
    int err = imp("data/test_import.xls", null, null, &d, cast(Gretlopt) 0, prn, null);

    check(err == 0, "xls_get_data returns 0 (first sheet)");
    if (!err) {
        check(d.v == 4 && d.n == 8, "3 series + const, 8 observations");
        check(strcmp(d.varname[1], "x") == 0 && strcmp(d.varname[3], "z") == 0, "names");
        check(d.Z[1][0] == 1 && d.Z[2][1] == 3 && d.Z[3][1] == 4, "numeric values");
        check(na(d.Z[3][2]), "blank cell is NA");
        clear_datainfo(&d, DataClearCode.full);
    }

    DatasetRaw d2;
    char[8] sname = 0;

    sname[0 .. 4] = "Strs";
    err = imp("data/test_import.xls", null, sname.ptr, &d2, cast(Gretlopt) 0, prn, null);
    check(err == 0, "sheet selected by name");
    if (!err) {
        check(d2.n == 4, "four observations");
        if (d2.S !is null) {
            check(strcmp(d2.S[2], "c") == 0 && d2.Z[1][2] == 3, "string column used as obs labels");
        } else {
            check(d2.v == 3 && is_string_valued(&d2, 1) != 0, "string column imported as string-valued");
        }
        clear_datainfo(&d2, DataClearCode.full);
    }

    DatasetRaw d3;
    sname[0 .. 5] = "Date\0";
    sname[0 .. 5] = "Dates";
    err = imp("data/test_import.xls", null, sname.ptr, &d3, cast(Gretlopt) 0, prn, null);
    check(err == 0, "dates sheet");
    if (!err) {
        check(d3.n == 8 && d3.pd == 4 && d3.structure == TsCodes.time_series, "quarterly time series detected");
        check(strcmp(d3.stobs.ptr, "2000:1") == 0, "starting observation");
        check(d3.Z[1][7] == 73.5, "last value");
        clear_datainfo(&d3, DataClearCode.full);
    }

    DatasetRaw d4;
    sname[0 .. 5] = "Nope\0";
    err = imp("data/test_import.xls", null, sname.ptr, &d4, cast(Gretlopt) 0, prn, null);
    check(err != 0, "unknown sheet name is an error");

    err = imp("build/no_such_file.xls", null, null, &d4, cast(Gretlopt) 0, prn, null);
    check(err != 0, "missing file is an error");

    err = imp("test_excel_import.d", null, null, &d4, cast(Gretlopt) 0, prn, null);
    check(err != 0, "non-OLE file is an error");

    gretl_print_destroy(prn);
    printf("\n%d checks, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
