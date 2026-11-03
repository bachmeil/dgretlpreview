/* Smoke tests for ods_import.d: a small .ods (a zip with mimetype and content.xml) is imported through
   get_plugin_function("ods_get_data"): names, numeric values, a blank cell, a string column, a repeated
   row, table selection by name and number, the empty-table pruning and error codes.

   Build and run with testing/run_tests.sh. */

import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;
import core.stdc.string : strcmp;

import std.file : chdir, exists, getcwd, mkdirRecurse, remove, rmdirRecurse, write;
import std.zip : ArchiveMember, CompressionMethod, ZipArchive;

import gretl.dataioh : DataClearCode;
import gretl.dataset : clear_datainfo, is_string_valued, series_get_string_for_obs;
import gretl.gretl_paths : gretl_set_path_by_name;
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

string scell(string s)
{
    return `<table:table-cell office:value-type="string"><text:p>` ~ s ~ `</text:p></table:table-cell>`;
}

string ncell(string v)
{
    return `<table:table-cell office:value-type="float" office:value="` ~ v ~ `"><text:p>` ~ v ~ `</text:p></table:table-cell>`;
}

string row(string cells, int repeat = 1)
{
    import std.conv : to;

    return `<table:table-row` ~ (repeat > 1 ? ` table:number-rows-repeated="` ~ repeat.to!string ~ `"` : "") ~ `>` ~ cells ~ `</table:table-row>`;
}

void make_ods(string fname)
{
    auto zip = new ZipArchive();
    auto mt = new ArchiveMember();

    mt.name = "mimetype";
    mt.expandedData = cast(ubyte[]) "application/vnd.oasis.opendocument.spreadsheet".dup;
    mt.compressionMethod = CompressionMethod.none;
    zip.addMember(mt);

    auto c = new ArchiveMember();

    c.name = "content.xml";
    c.expandedData = cast(ubyte[]) (`<?xml version="1.0"?><office:document-content xmlns:office="o" xmlns:table="t" xmlns:text="x">`
        ~ `<office:body><office:spreadsheet>`
        ~ `<table:table table:name="Data">`
        ~ row(scell("x") ~ scell("y") ~ scell("city"))
        ~ row(ncell("1.5") ~ ncell("10") ~ scell("Paris"))
        ~ row(ncell("2.5") ~ `<table:table-cell/>` ~ scell("Rome"))
        ~ row(ncell("3.5") ~ ncell("30") ~ scell("Paris"), 2)
        ~ `</table:table>`
        ~ `<table:table table:name="Blank"><table:table-row><table:table-cell/></table:table-row></table:table>`
        ~ `<table:table table:name="Second">`
        ~ row(scell("a") ~ scell("b")) ~ row(ncell("7") ~ ncell("8"))
        ~ `</table:table>`
        ~ `</office:spreadsheet></office:body></office:document-content>`).dup;
    c.compressionMethod = CompressionMethod.deflate;
    zip.addMember(c);
    write(fname, zip.build());
}

alias Importer = extern (C) int function(const(char) *, int *, char *, DatasetRaw *, Gretlopt, Prn *, void *);

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    string home = getcwd();

    if (exists("build/ods_dot")) {
        rmdirRecurse("build/ods_dot");
    }
    mkdirRecurse("build/ods_dot");
    gretl_set_path_by_name("dotdir", (home ~ "/build/ods_dot/\0").ptr);
    make_ods("build/ods_test.ods");

    Importer imp = cast(Importer) get_plugin_function("ods_get_data");

    check(imp !is null, "get_plugin_function returns the D importer");

    Prn * prn = gretl_print_new(PrnType.buffer, null);
    DatasetRaw d;
    int err = imp("build/ods_test.ods", null, null, &d, cast(Gretlopt) 0, prn, null);

    check(err == 0, "ods_get_data returns 0");
    if (!err) {
        check(d.v == 4 && d.n == 4, "3 series + const, 4 observations (one row repeated)");
        check(strcmp(d.varname[1], "x") == 0 && strcmp(d.varname[3], "city") == 0, "names");
        check(d.Z[1][0] == 1.5 && d.Z[2][0] == 10 && na(d.Z[2][1]), "numeric values and a blank cell");
        check(d.Z[1][3] == 3.5 && d.Z[2][3] == 30, "repeated row copied");
        check(is_string_valued(&d, 3) != 0, "city is string-valued");
        const(char) * s = series_get_string_for_obs(&d, 3, 1);
        check(s !is null && strcmp(s, "Rome") == 0, "string value");
        clear_datainfo(&d, DataClearCode.full);
    }

    chdir(home);
    DatasetRaw d2;
    char[8] sname = 0;

    sname[0 .. 6] = "Second";
    err = imp("build/ods_test.ods", null, sname.ptr, &d2, cast(Gretlopt) 0, prn, null);
    check(err == 0 && d2.v == 3 && d2.Z[2][0] == 8, "table selected by name (the blank table was pruned)");
    if (!err) clear_datainfo(&d2, DataClearCode.full);

    chdir(home);
    DatasetRaw d3;
    sname[0 .. 5] = "Nope\0";
    err = imp("build/ods_test.ods", null, sname.ptr, &d3, cast(Gretlopt) 0, prn, null);
    check(err != 0, "unknown table name is an error");

    chdir(home);
    err = imp("build/no_such_file.ods", null, null, &d3, cast(Gretlopt) 0, prn, null);
    check(err != 0, "missing file is an error");

    gretl_print_destroy(prn);
    chdir(home);
    rmdirRecurse("build/ods_dot");
    remove("build/ods_test.ods");
    printf("\n%d checks, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
