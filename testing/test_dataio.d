/* Smoke tests for dataio.d: date <-> observation-number conversion,
   number-format checks, CSV/octave round trips through temporary
   files, observation markers, transpose and nulldata creation.

   Build and run with testing/run_tests.sh. */

import core.stdc.stdio : _IOLBF, setvbuf, stdout, FILE, fclose, fgets, fopen, fprintf, printf, remove;
import core.stdc.stdlib : free;
import core.stdc.string : strcmp, strcpy, strlen, strstr;

import gretl.dataio;
import gretl.dataset : create_new_dataset, destroy_dataset, clear_datainfo, datainfo_new;
import gretl.dataioh : DataClearCode, GretlFileType;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_get_buffer, gretl_print_new;
import gretl.gretl_prnh : PrnType;
import gretl.libgretlh;
import gretl.missingh : NADBL, na;

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

DatasetRaw * make_dataset(int T)
{
    DatasetRaw * dset = create_new_dataset(3, T, 0);

    strcpy(dset.varname[0], "const");
    strcpy(dset.varname[1], "x");
    strcpy(dset.varname[2], "y");
    dset.pd = 4;
    dset.structure = TsCodes.time_series;
    dset.sd0 = 1990.1;
    strcpy(dset.stobs.ptr, "1990:1");
    dset.t1 = 0;
    dset.t2 = T - 1;
    strcpy(dset.endobs.ptr, "1991:4");

    foreach (t; 0 .. T) {
        dset.Z[1][t] = t + 0.5;
        dset.Z[2][t] = 10.0 * (t + 1);
    }

    return dset;
}

void write_file(const(char) * fname, string content)
{
    FILE * fp = fopen(fname, "w");

    fprintf(fp, "%.*s", cast(int) content.length, content.ptr);
    fclose(fp);
}

void test_dates()
{
    DatasetRaw * dset = make_dataset(8);
    char[OBSLEN] lbl = 0;

    check(dateton("1990:1", dset) == 0, "dateton 1990:1 == 0");
    check(dateton("1991:2", dset) == 5, "dateton 1991:2 == 5");
    ntolabel(lbl.ptr, 6, dset);
    check(strcmp(lbl.ptr, "1991:3") == 0, "ntolabel(6) == 1991:3");
    check(get_date_x(4, "1991:3") == 1991.3, "get_date_x 1991:3");

    int[2] mm;
    date_maj_min(6, dset, &mm[0], &mm[1]);
    check(mm[0] == 1991 && mm[1] == 3, "date_maj_min(6) == 1991, 3");

    destroy_dataset(dset);
}

void test_checks()
{
    double[4] x = [1.5, 2.25, 3.125, 4.0];

    check(get_precision(x.ptr, 4, 10) == 3, "get_precision == 3");
    check(check_varname("abc_1") == 0, "check_varname accepts abc_1");
    check(check_varname("1abc") != 0, "check_varname rejects 1abc");
    check(check_atof("3.25") == 0, "check_atof 3.25");
    check(check_atof("3.2x") != 0, "check_atof rejects 3.2x");
    check(check_atof("") == 0, "check_atof accepts blank");
    check(check_atoi("42") == 0, "check_atoi 42");
    check(check_atoi("4.2") != 0, "check_atoi rejects 4.2");
    check(check_atoi("99999999999") != 0, "check_atoi out of range");
}

void test_nulldata()
{
    DatasetRaw dset;
    Prn * prn = gretl_print_new(PrnType.buffer, null);

    int err = open_nulldata(&dset, 0, 25, cast(Gretlopt) 0, prn);
    check(err == 0 && dset.n == 25, "open_nulldata(25)");
    clear_datainfo(&dset, DataClearCode.full);
    gretl_print_destroy(prn);
}

void test_csv_write()
{
    DatasetRaw * dset = make_dataset(8);
    Prn * prn = gretl_print_new(PrnType.buffer, null);
    int[3] list = [2, 1, 2];
    char[256] line = 0;
    const(char) * fname = "build/dataio_test.csv";

    int err = write_data(fname, list.ptr, dset, Gretlopt.OPT_C, prn);
    check(err == 0, "write_data csv returns 0");

    FILE * fp = fopen(fname, "r");
    bool got_header = false, got_row = false;

    if (fp !is null) {
        while (fgets(line.ptr, cast(int) line.sizeof, fp)) {
            if (strstr(line.ptr, "obs") && strstr(line.ptr, "x") && strstr(line.ptr, "y")) {
                got_header = true;
            }
            if (strstr(line.ptr, "1991Q3") != null) {
                got_row = true;
            }
        }
        fclose(fp);
    }
    check(got_header, "csv header has obs, x, y");
    check(got_row, "csv uses ISO quarterly label 1991Q3");
    check(strcmp(detect_filetype(cast(char *) "build/dataio_test.csv".ptr, cast(Gretlopt) 0) == GretlFileType.csv ? "csv" : "other",
                 "csv") == 0, "detect_filetype -> csv");

    remove(fname);
    gretl_print_destroy(prn);
    destroy_dataset(dset);
}

void test_octave()
{
    const(char) * fname = "build/dataio_test.m";
    DatasetRaw dset;
    Prn * prn = gretl_print_new(PrnType.buffer, null);

    write_file(fname,
        "# Created by test\n"
        ~ "# name: aa\n# type: matrix\n# rows: 3\n# columns: 2\n"
        ~ " 1 2\n 3 4\n 5 6\n\n"
        ~ "# name: bb\n# type: matrix\n# rows: 3\n# columns: 1\n"
        ~ " 7\n 8\n 9\n");

    int err = import_other(fname, GretlFileType.octave, &dset, cast(Gretlopt) 0, prn);
    check(err == 0, "import_other octave returns 0");
    check(dset.v == 4 && dset.n == 3, "octave: 4 series (incl. const), 3 obs");
    if (dset.v == 4 && dset.n == 3) {
        check(strcmp(dset.varname[1], "aa1") == 0 && strcmp(dset.varname[3], "bb") == 0,
              "octave: variable names aa1, aa2, bb");
        check(dset.Z[2][1] == 4.0 && dset.Z[3][2] == 9.0, "octave: values");
    }

    remove(fname);
    clear_datainfo(&dset, DataClearCode.full);
    gretl_print_destroy(prn);
}

void test_markers_transpose()
{
    DatasetRaw * dset = make_dataset(4);
    const(char) * fname = "build/dataio_test_markers.txt";

    dset.structure = TsCodes.cross_section;
    dset.pd = 1;
    dset.sd0 = 1;
    write_file(fname, "alpha\nbeta\ngamma\ndelta\n");

    int err = add_obs_markers_from_file(dset, fname);
    check(err == 0 && dset.markers != 0, "add_obs_markers_from_file");
    check(err == 0 && strcmp(dset.S[2], "gamma") == 0, "marker 2 == gamma");
    remove(fname);

    err = transpose_data(dset);
    check(err == 0 && dset.v == 5 && dset.n == 2, "transpose_data dimensions");
    if (err == 0 && dset.v == 5) {
        check(strcmp(dset.varname[1], "alpha") == 0 && strcmp(dset.varname[4], "delta") == 0,
              "transposed names come from markers");
        check(dset.Z[2][0] == 1.5 && dset.Z[2][1] == 20.0, "transposed values");
    }
    destroy_dataset(dset);
}

void test_import_info()
{
    DatasetRaw * dset = make_dataset(4);

    dataset_add_import_info(dset, "/some/dir/foo.csv", GretlFileType.csv);
    check(dset.descrip !is null && strstr(dset.descrip, "foo.csv") != null &&
          strstr(dset.descrip, "CSV") != null, "dataset_add_import_info note");
    destroy_dataset(dset);
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);
    test_dates();
    test_checks();
    test_nulldata();
    test_csv_write();
    test_octave();
    test_markers_transpose();
    test_import_info();

    printf("\n%d checks, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
