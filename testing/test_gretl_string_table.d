/* Smoke tests for gretl_string_table.d: series tables (lookup, adding strings, mapping, copying),
   gretl string tables (building, validation, attaching to a dataset), the built-in strings,
   file content retrieval and "@name" substitution.

   Build and run with testing/run_tests.sh. */

import core.stdc.stdio : _IOLBF, FILE, fclose, fopen, fputs, printf, remove, setvbuf, stdout;
import core.stdc.stdlib : free;
import core.stdc.string : strcmp;
version (Posix) import core.sys.posix.stdlib : setenv; else import gretl.winport : setenv;

import gretl.dataset : clear_datainfo, create_new_dataset, destroy_dataset, is_string_valued,
    series_get_string_for_obs;
import gretl.dataioh : DataClearCode;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_list : gretl_list_new;
import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_new;
import gretl.gretl_prnh : PrnType;
import gretl.gretl_string_table;
import gretl.libgretlh;
import gretl.missingh : NADBL, na;
import gretl.strutils : gretl_strdup, strings_array_new;

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

void test_series_table()
{
    int err = 0;
    char ** strs = strings_array_new(3);

    strs[0] = gretl_strdup("red");
    strs[1] = gretl_strdup("green");
    strs[2] = gretl_strdup("blue");
    SeriesTable * st = series_table_new(strs, 3, &err);

    check(err == 0 && st !is null && series_table_get_n_strings(st) == 3, "series_table_new");
    check(series_table_get_value(st, "green") == 2.0 && na(series_table_get_value(st, "pink")), "series_table_get_value");
    check(strcmp(series_table_get_string(st, 3.0), "blue") == 0 && series_table_get_string(st, 4.0) is null &&
          series_table_get_string(st, NADBL) is null, "series_table_get_string");

    check(series_table_add_string(st, "\"quoted\"") == 4 && series_table_get_value(st, "quoted") == 4.0,
          "series_table_add_string strips quotes");

    SeriesTable * st2 = series_table_copy(st);
    int[2] dummy;
    const(char) *[2] more = ["violet", "red"];

    check(st2 !is null && series_table_add_strings(st2, more.ptr, 2) == 0 && series_table_get_n_strings(st2) == 6,
          "series_table_copy / add_strings");

    int * map = series_table_map(st, st2);

    check(map !is null && map[0] == 4 && map[1] == 6 && map[4] == 4, "series_table_map (a duplicate string takes the later code)");
    free(map);
    series_table_destroy(st);
    series_table_destroy(st2);
}

void test_string_table()
{
    int[3] l = [2, 1, 2];
    Prn * prn = gretl_print_new(PrnType.buffer, null);
    GretlStringTable * gst = gretl_string_table_new(null);

    check(gst !is null, "gretl_string_table_new (empty)");
    /* build the table for series 1 from scratch, as the CSV reader does */
    check(gretl_string_table_index(gst, "apple", 1, 1, prn) == 1, "index: first string gets 1");
    check(gretl_string_table_index(gst, "pear", 1, 1, prn) == 2 && gretl_string_table_index(gst, "apple", 1, 1, prn) == 1,
          "index: new strings and repeats");
    check(in_string_table(gst, 1) == 1 && in_string_table(gst, 2) == 0, "in_string_table");
    check(gretl_string_table_validate(gst, cast(Gretlopt) 0) == 0, "validate: ordinary strings are fine");

    DatasetRaw * dset = create_new_dataset(2, 3, 0);

    dset.Z[1][0] = 1; dset.Z[1][1] = 2; dset.Z[1][2] = 1;
    check(gretl_string_table_finalize(gst, dset) == 1, "finalize attaches one table");
    check(is_string_valued(dset, 1) != 0, "series is string-valued");
    const(char) * s = series_get_string_for_obs(dset, 1, 1);

    check(s !is null && strcmp(s, "pear") == 0, "string for obs 2 is pear");
    gretl_string_table_destroy(gst);

    /* numeric-looking strings must be rejected by validate */
    gst = gretl_string_table_new(null);
    gretl_string_table_index(gst, "12", 1, 1, prn);
    gretl_string_table_index(gst, "34.5", 1, 1, prn);
    check(gretl_string_table_validate(gst, cast(Gretlopt) 0) == GretlError.data, "validate rejects quasi-numeric strings");
    gretl_string_table_destroy(gst);

    gretl_print_destroy(prn);
    destroy_dataset(dset);
}

void test_builtin_and_env()
{
    gretl_insert_builtin_string("workdir", "/home/user/work/");
    gretl_insert_builtin_string("pkgdir", "/opt/pkg");
    check(strcmp(get_built_in_string_by_name("workdir"), "/home/user/work") == 0, "built-in string drops the trailing slash");
    check(strcmp(get_built_in_string_by_name("pkgdir"), "/opt/pkg") == 0, "built-in string");
    check(strcmp(get_built_in_string_by_name("lang"), "") == 0 && get_built_in_string_by_name("nonesuch") is null,
          "unset built-in is empty, unknown is NULL");

    setenv("GRETL_ST_TEST", "value1", 1);
    int defined = 0, err = 0;
    char * v = gretl_getenv("GRETL_ST_TEST", &defined, &err);

    check(err == 0 && defined == 1 && strcmp(v, "value1") == 0, "gretl_getenv (defined)");
    free(v);
    v = gretl_getenv("GRETL_ST_NOPE_XYZ", &defined, &err);
    check(defined == 0 && strcmp(v, "") == 0, "gretl_getenv (undefined)");
    free(v);

    /* "@name" substitution */
    char[128] line = 0;
    string src = "print @pkgdir/file.txt\0";

    line[0 .. src.length] = src[];
    int subst = 0;

    check(substitute_named_strings(line.ptr, &subst) == 0 && subst == 1 && strcmp(line.ptr, "print /opt/pkg/file.txt") == 0,
          "substitute_named_strings");
}

void test_content()
{
    const(char) * fname = "build/st_content.txt";
    FILE * fp = fopen(fname, "w");

    fputs("hello\nworld\n", fp);
    fclose(fp);

    int err = 0;
    char * c = retrieve_file_content(fname, null, &err);

    check(err == 0 && c !is null && strcmp(c, "hello\nworld\n") == 0, "retrieve_file_content");
    free(c);
    err = 0;
    c = retrieve_file_content("build/no_such_file_xyz.txt", null, &err);
    check(err != 0 && c is null, "retrieve_file_content: missing file is an error");
    err = 0;
    c = retrieve_file_content("http://example.com/x.txt", null, &err);
    check(err != 0 && c is null, "retrieve_file_content: web resources not supported");
    remove(fname);

    err = 0;
    DatasetRaw * dset = create_new_dataset(2, 3, 0);
    char * d = retrieve_date_string(2, dset, &err);

    check(err == 0 && d !is null, "retrieve_date_string (no markers)");
    free(d);
    check(retrieve_date_string(9, dset, &err) is null && err == GretlError.data, "retrieve_date_string out of range");
    destroy_dataset(dset);
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);
    test_series_table();
    test_string_table();
    test_builtin_and_env();
    test_content();
    printf("\n%d checks, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
