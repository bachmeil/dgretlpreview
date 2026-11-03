/* Smoke tests for gnumeric_import.d: small .gnumeric workbooks (plain and gzipped XML) are
   imported through get_plugin_function("gnumeric_get_data"): names, numeric values, a string
   column, a date-formatted observation column, sheet selection and error codes.

   Build and run with testing/run_tests.sh. */

import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;
import core.stdc.string : strcmp;

import std.file : remove, write;
import std.zlib : Compress, HeaderFormat;

import gretl.dataioh : DataClearCode;
import gretl.dataset : clear_datainfo, is_string_valued, series_get_string_for_obs;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_new;
import gretl.gretl_prnh : PrnType;
import gretl.libgretlh;
import gretl.missingh : na;
import gretl.plugins : get_plugin_function;

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

string cell(int r, int c, int vt, string v, string extra = "")
{
    import std.conv : to;

    return `<gnm:Cell Row="` ~ r.to!string ~ `" Col="` ~ c.to!string ~ `" ValueType="` ~ vt.to!string ~ `"` ~ extra
        ~ `>` ~ v ~ `</gnm:Cell>`;
}

string workbook()
{
    return `<?xml version="1.0"?><gnm:Workbook xmlns:gnm="http://www.gnumeric.org/v10.dtd">`
        ~ `<gnm:SheetNameIndex><gnm:SheetName>First</gnm:SheetName><gnm:SheetName>Second</gnm:SheetName></gnm:SheetNameIndex>`
        ~ `<gnm:Sheets><gnm:Sheet><gnm:Name>First</gnm:Name><gnm:Cells>`
        ~ cell(0, 0, 60, "x") ~ cell(0, 1, 60, "y") ~ cell(0, 2, 60, "city")
        ~ cell(1, 0, 40, "1.5") ~ cell(1, 1, 30, "10") ~ cell(1, 2, 60, "Paris")
        ~ cell(2, 0, 40, "2.5") ~ cell(2, 1, 30, "20") ~ cell(2, 2, 60, "Rome")
        ~ cell(3, 0, 40, "3.5") ~ cell(3, 1, 10, "") ~ cell(3, 2, 60, "Paris")
        ~ `</gnm:Cells></gnm:Sheet>`
        ~ `<gnm:Sheet><gnm:Name>Second</gnm:Name><gnm:Cells>`
        ~ cell(0, 0, 60, "a") ~ cell(0, 1, 60, "b") ~ cell(1, 0, 30, "7") ~ cell(1, 1, 30, "8")
        ~ `</gnm:Cells></gnm:Sheet></gnm:Sheets></gnm:Workbook>`;
}

alias Importer = extern (C) int function(const(char) *, int *, char *, DatasetRaw *, Gretlopt, Prn *, void *);

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    write("build/gnm_test.gnumeric", workbook());

    /* a gzipped copy: gnumeric files normally are */
    {
        auto c = new Compress(6, HeaderFormat.gzip);
        auto data = cast(ubyte[]) c.compress(workbook().dup);

        data ~= cast(ubyte[]) c.flush();
        write("build/gnm_test_gz.gnumeric", data);
    }

    Importer imp = cast(Importer) get_plugin_function("gnumeric_get_data");

    check(imp !is null, "get_plugin_function returns the D importer");

    Prn * prn = gretl_print_new(PrnType.buffer, null);

    foreach (fn; ["build/gnm_test.gnumeric", "build/gnm_test_gz.gnumeric"]) {
        DatasetRaw d;
        int err = imp((fn ~ "\0").ptr, null, null, &d, cast(Gretlopt) 0, prn, null);

        check(err == 0, fn ~ ": gnumeric_get_data returns 0");
        if (!err) {
            check(d.v == 4 && d.n == 3, "3 series + const, 3 observations");
            check(strcmp(d.varname[1], "x") == 0 && strcmp(d.varname[3], "city") == 0, "names");
            check(d.Z[1][0] == 1.5 && d.Z[2][1] == 20 && na(d.Z[2][2]), "numeric values and a missing cell");
            check(is_string_valued(&d, 3) != 0, "city is string-valued");
            const(char) * s = series_get_string_for_obs(&d, 3, 1);
            check(s !is null && strcmp(s, "Rome") == 0, "string value");
            clear_datainfo(&d, DataClearCode.full);
        }
    }

    /* second sheet, by name and by number */
    {
        DatasetRaw d;
        char[8] sname = 0;

        sname[0 .. 6] = "Second";
        int err = imp("build/gnm_test.gnumeric", null, sname.ptr, &d, cast(Gretlopt) 0, prn, null);

        check(err == 0 && d.v == 3 && d.n == 1 && d.Z[2][0] == 8, "sheet selected by name");
        if (!err) clear_datainfo(&d, DataClearCode.full);

        DatasetRaw d2;
        int[4] list = [3, 2, 0, 0];

        err = imp("build/gnm_test.gnumeric", list.ptr, null, &d2, cast(Gretlopt) 0, prn, null);
        check(err == 0 && d2.v == 3 && d2.Z[1][0] == 7, "sheet selected by number");
        if (!err) clear_datainfo(&d2, DataClearCode.full);

        DatasetRaw d3;
        sname[0 .. 5] = "Nope\0";
        err = imp("build/gnm_test.gnumeric", null, sname.ptr, &d3, cast(Gretlopt) 0, prn, null);
        check(err != 0, "unknown sheet name is an error");

        err = imp("build/no_such_file.gnumeric", null, null, &d3, cast(Gretlopt) 0, prn, null);
        check(err != 0, "missing file is an error");
    }

    gretl_print_destroy(prn);
    remove("build/gnm_test.gnumeric");
    remove("build/gnm_test_gz.gnumeric");
    printf("\n%d checks, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
