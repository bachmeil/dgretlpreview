/* Tests for dgretl.d, the D front end: no library start-up call appears anywhere in this file (dgretl starts the library
   before main), datasets and matrices are reference counted and free themselves, options are optional named
   arguments, and failures are GretlException. A leak check compares the heap in use before and after many imports.

   Build and run with testing/run_tests.sh. */

import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import std.algorithm : canFind;
import std.file : remove, tempDir, write;
import std.process : thisProcessID;
import std.conv : to;

import gretl.dataset : create_new_dataset;
import dgretl;
import gretl.gretl_utils : libgretl_init;

struct MallInfo {
    int arena, ordblks, smblks, hblks, hblkhd, usmblks, fsmblks, uordblks, fordblks, keepcost;
}
extern (C) MallInfo mallinfo();

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

string tmp(string name)
{
    return tempDir() ~ "/test_dgretl_" ~ to!string(thisProcessID) ~ "_" ~ name;
}

/* the exception message of a failing call, or "" if it did not throw */
string failure(void delegate() dg)
{
    try {
        dg();
    } catch (GretlException e) {
        return e.msg;
    }
    return "";
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    string f1 = tmp("a.csv"), f2 = tmp("b.csv"), f3 = tmp("c.csv");
    write(f1, "a,b\n1,2\n3,4\n5,6\n");
    write(f2, "year,gdp\n2001,1.5\n2002,2.5\n2003,3.5\n");
    write(f3, "1,2\n3,4\n");
    scope (exit) { remove(f1); remove(f2); remove(f3); }

    /* a plain import */
    auto d = importCsv(f1);
    check(d.nobs == 3 && d.nseries == 2 && d.names == ["a", "b"], "importCsv: observations, series and names");
    check(d["a"] == [1.0, 3.0, 5.0] && d["b"] == [2.0, 4.0, 6.0] && d[1] == [2.0, 4.0, 6.0], "series by name and by number");
    check(d.ptr !is null && d.ptr.n == 3, "the raw pointer is available");

    /* a first column of years is the observation index, unless allColumns is given */
    auto y = importCsv(f2);
    auto y2 = importCsv(f2, "allColumns");
    check(y.names == ["gdp"] && y2.names == ["year", "gdp"] && y2["year"] == [2001.0, 2002.0, 2003.0], "the label column and the allColumns option");

    /* no header row */
    auto nh = importCsv(f3, "noHeader");
    check(nh.nobs == 2 && nh.names == ["v1", "v2"] && nh[0] == [1.0, 3.0], "the noHeader option names the series v1, v2");

    /* copies share the data, which outlives the original */
    Dataset keep;
    check(keep.nobs == 0 && keep.names.length == 0 && keep.ptr is null, "a default Dataset is empty");
    {
        auto tmpd = importCsv(f1);
        keep = tmpd;
        auto another = keep;     /* three handles on one dataset */
        check(another.nobs == 3, "a copy reads the same data");
    }
    check(keep.nobs == 3 && keep["a"][2] == 5.0, "the data survive when the other handles are gone");

    /* failures */
    string m1 = failure({ importCsv(tmp("no_such_file.csv")); });
    check(m1.canFind("importCsv") && m1.canFind("no_such_file.csv") && m1.canFind("Couldn't open"), "a missing file throws GretlException with the name and the importer's message");
    check(failure({ auto s = d["nonesuch"]; }).canFind("nonesuch") && failure({ auto s = d[9]; }).canFind("out of range"), "unknown series name and number");

    /* adopting a raw dataset made with the low-level functions */
    {
        auto raw = create_new_dataset(2, 4, 0);
        raw.Z[1][0] = 7.0;
        auto ad = Dataset(raw);
        check(ad.nobs == 4 && ad.nseries == 1 && ad[0][0] == 7.0, "Dataset(raw) takes over a raw dataset");
    }

    /* the constructor that makes a new dataset: 3 series counting the constant, 5 observations */
    {
        auto nd = Dataset(3, 5);
        nd.ptr.Z[1][4] = 1.5;
        check(nd.nobs == 5 && nd.nseries == 2 && nd[0].length == 5 && nd[0][4] == 1.5, "Dataset(nvar, nobs) makes a new dataset");
    }

    /* several options at once, as strings (and "verbose" shows the importer's messages, here to the test's output) */
    check(importCsv(f3, "noHeader", "allColumns").nobs == 2 && importCsv(f3, "noHeader", "verbose").nseries == 2, "several string options, including verbose");

    /* the same options given at compile time: importCsv!("noHeader", "allColumns")(file) */
    auto ct = importCsv!("noHeader", "allColumns")(f3);
    check(ct.nobs == 2 && ct.names == ["v1", "v2"] && ct[1] == [2.0, 4.0], "compile-time options give the same result");
    check(importCsv!"verbose"(f1).nobs == 3 && importCsv!("allColumns")(f2).names == ["year", "gdp"], "a single compile-time option; allColumns keeps the year");
    static assert(__traits(compiles, importCsv!("noHeader", "allColumns")(f3)));
    static assert(!__traits(compiles, importCsv!("noHeaders")(f3)));      /* a misspelt option does not compile */
    static assert(!__traits(compiles, importCsv!(42)(f3)));               /* options must be strings */

    /* spreadsheets: errors only (the Kansas example reads a real worksheet) */
    check(failure({ importSpreadsheet(tmp("nothing.csv")); }).canFind("must end in"), "importSpreadsheet refuses other file types");
    check(failure({ importSpreadsheet(tmp("nothing.xlsx")); }).canFind("importSpreadsheet"), "importSpreadsheet reports a missing workbook");

    /* matrices */
    auto m = importCsvAsMatrix(f1);
    check(m.nrow == 3 && m.ncol == 2 && m[0, 0] == 1.0 && m[2, 1] == 6.0 && m[1, 0] == 3.0, "importCsvAsMatrix: size and elements (row, column)");
    m[1, 1] = 40.0;
    auto m2 = m;
    check(m2[1, 1] == 40.0 && m.ptr == m2.ptr, "assignment through an index; a copy shares the matrix");
    check(failure({ auto q = m[3, 0]; }).canFind("out of range") && failure({ importCsvAsMatrix(tmp("no_such_file.csv")); }).canFind("importCsvAsMatrix"), "matrix index and file errors");
    auto my = importCsvAsMatrix(f2);
    check(my.ncol == 2 && my[0, 0] == 2001.0, "the matrix keeps the label column");

    /* the library start-up is idempotent: calling it again changes nothing */
    libgretl_init();
    check(importCsv(f1).nobs == 3, "libgretl_init() a second time is harmless");

    /* heap check: after a warm-up, many imports must not grow the heap in use */
    foreach (i; 0 .. 5) { auto w = importCsv(f1); auto wm = importCsvAsMatrix(f1); }
    int before = mallinfo().uordblks;
    foreach (i; 0 .. 300) {
        auto t = importCsv(f1);
        auto u = importCsvAsMatrix(f2);
    }
    int after = mallinfo().uordblks;
    printf("  heap in use: %d bytes before, %d after 300 imports\n", before, after);
    check(after - before < 64 * 1024, "importing and releasing 300 times does not leak");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
