/* Smoke tests for xlsx_import.d: a small workbook is written with std.zip (workbook, relationships,
   one sheet with numeric and shared-string cells, a shared-strings table) and imported through
   get_plugin_function("xlsx_get_data"), as gretl's open command does. gretl_unzip() is a stand-in
   in stubs.d until gretl_zip.c is ported.

   Build and run with testing/run_tests.sh. */

import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;
import core.stdc.string : strcmp;

import std.file : exists, mkdirRecurse, read, rmdirRecurse, write;
import std.zip : ArchiveMember, CompressionMethod, ZipArchive;

import gretl.dataioh : DataClearCode;
import gretl.dataset : clear_datainfo, is_string_valued, series_get_string_for_obs;
import gretl.gretl_errorsh : GretlError;
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

void add(ZipArchive zip, string name, string content)
{
    auto m = new ArchiveMember();

    m.name = name;
    m.expandedData = cast(ubyte[]) content.dup;
    m.compressionMethod = CompressionMethod.deflate;
    zip.addMember(m);
}

void make_xlsx(string fname)
{
    auto zip = new ZipArchive();

    add(zip, "xl/workbook.xml",
        `<?xml version="1.0"?><workbook xmlns:r="x"><sheets>`
        ~ `<sheet name="Data" sheetId="1" r:id="rId1"/><sheet name="Empty" sheetId="2" r:id="rId2"/>`
        ~ `</sheets></workbook>`);
    add(zip, "xl/_rels/workbook.xml.rels",
        `<?xml version="1.0"?><Relationships>`
        ~ `<Relationship Id="rId1" Target="worksheets/sheet1.xml"/>`
        ~ `<Relationship Id="rId2" Target="worksheets/sheet2.xml"/></Relationships>`);
    add(zip, "xl/sharedStrings.xml",
        `<?xml version="1.0"?><sst count="5" uniqueCount="5"><si><t>x</t></si><si><t>y</t></si>`
        ~ `<si><t>city</t></si><si><t>Paris</t></si><si><t>Rome</t></si></sst>`);
    add(zip, "xl/worksheets/sheet1.xml",
        `<?xml version="1.0"?><worksheet><sheetData>`
        ~ `<row r="1"><c r="A1" t="s"><v>0</v></c><c r="B1" t="s"><v>1</v></c><c r="C1" t="s"><v>2</v></c></row>`
        ~ `<row r="2"><c r="A2"><v>1.5</v></c><c r="B2"><v>10</v></c><c r="C2" t="s"><v>3</v></c></row>`
        ~ `<row r="3"><c r="A3"><v>2.5</v></c><c r="B3"><v>20</v></c><c r="C3" t="s"><v>4</v></c></row>`
        ~ `<row r="4"><c r="A4"><v>3.5</v></c><c r="C4" t="s"><v>3</v></c></row>`
        ~ `</sheetData></worksheet>`);
    add(zip, "xl/worksheets/sheet2.xml", `<?xml version="1.0"?><worksheet><sheetData/></worksheet>`);
    write(fname, zip.build());
}

alias Importer = extern (C) int function(const(char) *, int *, char *, DatasetRaw *, Gretlopt, Prn *, void *);

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    import std.file : chdir, getcwd;

    string home = getcwd();
    string dot = home ~ "/build/xlsx_dot/\0";

    if (exists("build/xlsx_dot")) {
        rmdirRecurse("build/xlsx_dot");
    }
    mkdirRecurse("build/xlsx_dot");
    gretl_set_path_by_name("dotdir", dot.ptr);
    make_xlsx("build/xlsx_test.xlsx");

    Importer imp = cast(Importer) get_plugin_function("xlsx_get_data");

    check(imp !is null, "get_plugin_function returns the D importer");

    DatasetRaw d;
    Prn * prn = gretl_print_new(PrnType.buffer, null);
    int err = imp("build/xlsx_test.xlsx", null, null, &d, cast(Gretlopt) 0, prn, null);

    check(err == 0, "xlsx_get_data returns 0");
    if (!err) {
        check(d.v == 4 && d.n == 3, "3 series + const, 3 observations");
        check(strcmp(d.varname[1], "x") == 0 && strcmp(d.varname[2], "y") == 0 && strcmp(d.varname[3], "city") == 0,
              "names x, y, city");
        check(d.Z[1][0] == 1.5 && d.Z[1][2] == 3.5 && d.Z[2][1] == 20, "numeric values");
        check(na(d.Z[2][2]), "missing cell is NA");
        check(is_string_valued(&d, 3) != 0, "city is string-valued");
        const(char) * s = series_get_string_for_obs(&d, 3, 1);
        check(s !is null && strcmp(s, "Rome") == 0, "string value read through the shared strings");
        check(d.descrip !is null, "import note added");
        clear_datainfo(&d, DataClearCode.full);
    }

    /* choose the sheet by name: the empty sheet has no data */
    DatasetRaw d2;
    char[8] sname = 0;
    sname[0 .. 5] = "Empty";
    chdir(home);
    err = imp("build/xlsx_test.xlsx", null, sname.ptr, &d2, cast(Gretlopt) 0, prn, null);
    check(err == GretlError.data, "asking for the empty sheet fails with E_DATA");

    sname[0 .. 5] = "Nope\0";
    chdir(home);
    err = imp("build/xlsx_test.xlsx", null, sname.ptr, &d2, cast(Gretlopt) 0, prn, null);
    check(err == GretlError.data, "unknown sheet name -> E_DATA");

    chdir(home);
    err = imp("build/no_such_file.xlsx", null, null, &d2, cast(Gretlopt) 0, prn, null);
    check(err == GretlError.fopen, "missing file -> E_FOPEN");

    gretl_print_destroy(prn);
    chdir(home); /* the importer leaves the process in the dotdir, as the C does */
    rmdirRecurse("build/xlsx_dot");
    import std.file : remove;
    remove("build/xlsx_test.xlsx");
    printf("\n%d checks, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
