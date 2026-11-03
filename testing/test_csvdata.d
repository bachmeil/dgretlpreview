/* Smoke tests for csvdata.d: CSV files written to temporary files are
   imported with import_csv() and the resulting dataset is examined
   (dimensions, names, values, date detection, delimiters, missing
   values, header-less files), plus probe_csv() and peek_at_csv().

   Build and run with testing/run_tests.sh. */

import core.stdc.stdio : _IOLBF, FILE, fclose, fopen, fprintf, printf, remove, setvbuf, stdout;
import core.stdc.string : strcmp;

import gretl.csvdata;
import gretl.dataset : clear_datainfo, is_string_valued, series_get_string_for_obs;
import gretl.dataioh : DataClearCode;
import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_get_buffer, gretl_print_new;
import gretl.gretl_prnh : PrnType;
import gretl.libgretlh;
import gretl.missingh : na;

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

const(char) * fname = "build/csvdata_test.csv";

void write_file(string content)
{
    FILE * fp = fopen(fname, "w");

    fprintf(fp, "%.*s", cast(int) content.length, content.ptr);
    fclose(fp);
}

void test_dated()
{
    DatasetRaw d;
    Prn * prn = gretl_print_new(PrnType.buffer, null);

    write_file("obs,x,y\n1990:1,1.5,10\n1990:2,2.5,20\n1990:3,3.5,30\n1990:4,4.5,40\n"
               ~ "1991:1,5.5,50\n1991:2,6.5,60\n");
    int err = import_csv(fname, &d, cast(Gretlopt) 0, prn);

    check(err == 0, "import_csv (quarterly dates) returns 0");
    if (!err) {
        check(d.v == 3 && d.n == 6, "2 series, 6 obs");
        check(strcmp(d.varname[1], "x") == 0 && strcmp(d.varname[2], "y") == 0, "names x, y");
        check(d.pd == 4 && strcmp(d.stobs.ptr, "1990:1") == 0, "quarterly frequency, start 1990:1");
        check(d.Z[1][1] == 2.5 && d.Z[2][5] == 60, "values");
        check(d.structure == TsCodes.time_series, "time-series structure");
        clear_datainfo(&d, DataClearCode.full);
    }
    gretl_print_destroy(prn);
}

void test_undated_missing()
{
    DatasetRaw d;
    Prn * prn = gretl_print_new(PrnType.buffer, null);

    write_file("a,b\n1,2\n3,NA\n5,6\n,8\n");
    int err = import_csv(fname, &d, cast(Gretlopt) 0, prn);

    check(err == 0, "import_csv (undated, NA) returns 0");
    if (!err) {
        check(d.v == 3 && d.n == 4, "2 series, 4 obs");
        check(d.Z[1][2] == 5 && na(d.Z[2][1]) && na(d.Z[1][3]), "NA and blank cells are missing");
        check(d.pd == 1 && d.structure == TsCodes.cross_section, "undated cross-section");
        clear_datainfo(&d, DataClearCode.full);
    }
    gretl_print_destroy(prn);
}

void test_semicolon_comma()
{
    DatasetRaw d;
    Prn * prn = gretl_print_new(PrnType.buffer, null);

    write_file("p;q\n1,5;2,25\n3,5;4,75\n");
    int err = import_csv(fname, &d, cast(Gretlopt) 0, prn);

    check(err == 0, "import_csv (semicolons, decimal commas) returns 0");
    if (!err) {
        check(d.n == 2 && d.Z[1][0] == 1.5 && d.Z[2][1] == 4.75, "decimal-comma values read");
        clear_datainfo(&d, DataClearCode.full);
    }
    gretl_print_destroy(prn);
}

void test_noheader()
{
    DatasetRaw d;
    Prn * prn = gretl_print_new(PrnType.buffer, null);

    write_file("1,2,3\n4,5,6\n7,8,9\n");
    int err = import_csv(fname, &d, Gretlopt.OPT_H, prn);

    check(err == 0, "import_csv (no header) returns 0");
    if (!err) {
        check(d.v == 4 && d.n == 3 && strcmp(d.varname[1], "v1") == 0, "autonamed v1..v3, 3 obs");
        check(d.Z[3][2] == 9, "values");
        clear_datainfo(&d, DataClearCode.full);
    }
    gretl_print_destroy(prn);
}

void test_thousands()
{
    DatasetRaw d;
    Prn * prn = gretl_print_new(PrnType.buffer, null);

    write_file("v;w\n1,234.5;1\n2,345.5;2\n10,000.25;3\n");
    int err = import_csv(fname, &d, cast(Gretlopt) 0, prn);

    check(err == 0, "import_csv (thousands separators) returns 0");
    if (!err) {
        check(d.Z[1][0] == 1234.5 && d.Z[1][2] == 10000.25, "thousands separators removed");
        clear_datainfo(&d, DataClearCode.full);
    }
    gretl_print_destroy(prn);
}

void test_strings()
{
    DatasetRaw d;
    Prn * prn = gretl_print_new(PrnType.buffer, null);

    write_file("city,pop\nParis,2.1\nRome,2.8\nParis,2.2\nOslo,0.7\n");
    int err = import_csv(fname, &d, cast(Gretlopt) 0, prn);

    check(err == 0, "import_csv (string column) returns 0");
    if (!err) {
        check(d.v == 3 && d.n == 4, "2 series, 4 obs");
        check(is_string_valued(&d, 1) != 0, "city is string-valued");
        const(char) * s0 = series_get_string_for_obs(&d, 1, 2);
        check(s0 !is null && strcmp(s0, "Paris") == 0, "obs 3 is Paris");
        check(d.Z[1][0] == d.Z[1][2] && d.Z[1][0] != d.Z[1][1], "equal strings share a code");
        clear_datainfo(&d, DataClearCode.full);
    }
    gretl_print_destroy(prn);
}

void test_probe_peek()
{
    char ** names = null;
    int nv = 0;
    Gretlopt opt = cast(Gretlopt) 0;

    write_file("alpha,beta,gamma\n1,2,3\n4,5,6\n");
    int err = probe_csv(fname, &names, &nv, &opt);

    check(err == 0 && nv == 4 && names !is null && strcmp(names[3], "gamma") == 0, "probe_csv finds variable names");

    Prn * prn = gretl_print_new(PrnType.buffer, null);

    err = peek_at_csv(fname, 2, prn);
    const(char) * buf = gretl_print_get_buffer(prn);

    check(err == 0 && buf !is null && strcmp(buf, "alpha,beta,gamma\n1,2,3\n") == 0, "peek_at_csv returns the first lines");
    gretl_print_destroy(prn);
}

void test_helpers()
{
    check(import_na_string("NA") == 1 && import_na_string("n/a") == 1 && import_na_string("3.5") == 0,
          "import_na_string defaults");
    check(import_obs_label("obs") == 1 && import_obs_label("Date") == 1 && import_obs_label("price") == 0,
          "import_obs_label");
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);
    import_na_init();
    test_helpers();
    test_dated();
    test_undated_missing();
    test_semicolon_comma();
    test_noheader();
    test_strings();
    test_thousands();
    test_probe_peek();
    remove(fname);
    printf("\n%d checks, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
