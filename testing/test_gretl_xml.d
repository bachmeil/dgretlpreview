/* Smoke tests for gretl_xml.d: the writers for matrices, lists and
   strings, the gdt data-file round trip (write, read, read series
   names, read the description) and a few of the XML readers.

   Build and run with testing/run_tests.sh. */

import core.stdc.stdio : _IOLBF, FILE, fclose, fopen, fprintf, printf, remove, setvbuf, stdout;
import core.stdc.stdlib : free;
import core.stdc.string : strcmp, strcpy, strlen, strstr;

import gretl.dataset : clear_datainfo, create_new_dataset, destroy_dataset, series_get_label,
    series_set_label;
import gretl.dataioh : DataClearCode;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_matrix : gretl_matrix_alloc, gretl_matrix_free;
import gretl.gretl_matrixh;
import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_get_buffer, gretl_print_new;
import gretl.gretl_prnh : PrnType;
import gretl.gretl_xml;
import gretl.libgretlh;
import gretl.missingh : NADBL, na;
import gretl.xmlh;

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
    strcpy(dset.endobs.ptr, "1991:4");
    dset.t1 = 0;
    dset.t2 = T - 1;
    series_set_label(dset, 1, "the x series");

    foreach (t; 0 .. T) {
        dset.Z[1][t] = t + 0.25;
        dset.Z[2][t] = (t == 3) ? NADBL : 10.0 * (t + 1);
    }

    return dset;
}

void test_serialize()
{
    Prn * prn = gretl_print_new(PrnType.buffer, null);
    GretlMatrixRaw * m = gretl_matrix_alloc(2, 2);
    int[4] list = [3, 1, 2, 5];

    m.val[0] = 1; m.val[1] = 2; m.val[2] = 3; m.val[3] = 4.5;

    check(gretl_matrix_serialize(m, "mm", prn) == 0, "gretl_matrix_serialize returns 0");
    gretl_list_serialize(list.ptr, "ll", prn);
    gretl_finite_scalar_serialize(2.5, "sc", prn);
    gretl_string_serialize("hello", "ss", prn);
    gretl_xml_put_tagged_string("tag", "a < b", prn);
    gretl_xml_put_int("n", 7, prn);

    const(char) * buf = gretl_print_get_buffer(prn);

    check(strstr(buf, "<gretl-matrix name=\"mm\" rows=\"2\" cols=\"2\">") != null, "matrix header");
    check(strstr(buf, "<list name=\"ll\">") != null && strstr(buf, "3 1 2 5 ") != null, "list serialized");
    check(strstr(buf, "value=\"2.5\"") != null, "scalar serialized");
    check(strstr(buf, "<tag>a &lt; b</tag>") != null, "xml special characters escaped");
    check(strstr(buf, "n=\"7\" ") != null, "put_int attribute");

    /* read the matrix back through libxml */
    xmlDocPtr doc;
    xmlNodePtr node;
    char * xml = null;

    const(char) * start = strstr(buf, "<gretl-matrix");
    const(char) * end = strstr(start, "</gretl-matrix>") + strlen("</gretl-matrix>");
    char[512] mbuf = 0;

    mbuf[0 .. end - start] = start[0 .. end - start];

    int err = gretl_xml_read_buffer(mbuf.ptr, "gretl-matrix", &doc, &node);
    check(err == 0, "gretl_xml_read_buffer parses matrix");
    if (!err) {
        int e2 = 0;
        GretlMatrixRaw * m2 = gretl_xml_get_matrix(node, doc, &e2);

        check(e2 == 0 && m2 !is null && m2.rows == 2 && m2.cols == 2 &&
              m2.val[3] == 4.5, "gretl_xml_get_matrix round trip");
        gretl_matrix_free(m2);
        xmlFreeDoc(doc);
    }

    gretl_matrix_free(m);
    gretl_print_destroy(prn);
}

void test_gdt_roundtrip()
{
    const(char) * fname = "build/xml_test.gdt";
    DatasetRaw * dset = make_dataset(8);
    int err = gretl_write_gdt(fname, null, dset, cast(Gretlopt) 0, 0);

    check(err == 0, "gretl_write_gdt returns 0");

    DatasetRaw in_;
    Prn * prn = gretl_print_new(PrnType.buffer, null);

    err = gretl_read_gdt(fname, &in_, cast(Gretlopt) 0, prn);
    check(err == 0, "gretl_read_gdt returns 0");
    if (!err) {
        check(in_.n == 8 && in_.v == 3, "read: 8 obs, 3 series");
        check(strcmp(in_.varname[1], "x") == 0 && strcmp(in_.varname[2], "y") == 0, "read: series names");
        check(in_.pd == 4 && strcmp(in_.stobs.ptr, "1990:1") == 0, "read: frequency and start obs");
        check(in_.Z[1][5] == 5.25 && in_.Z[2][1] == 20.0, "read: values");
        check(na(in_.Z[2][3]), "read: NA preserved");
        const(char) * lab = series_get_label(&in_, 1);
        check(lab !is null && strcmp(lab, "the x series") == 0, "read: series label");
    }

    char ** names = null;
    int nv = 0;

    err = gretl_read_native_varnames(fname, &names, &nv);
    check(err == 0 && nv == 3 && names !is null && strcmp(names[2], "y") == 0,
          "gretl_read_native_varnames");

    clear_datainfo(&in_, DataClearCode.full);
    gretl_print_destroy(prn);
    destroy_dataset(dset);
    remove(fname);
}

void test_description()
{
    const(char) * fname = "build/xml_test_desc.gdt";
    DatasetRaw * dset = make_dataset(4);

    dset.descrip = cast(char *) strcpy(cast(char *) malloc_(32), "A <small> test");
    gretl_write_gdt(fname, null, dset, cast(Gretlopt) 0, 0);

    int err = 0;
    char * d = gretl_get_gdt_description(fname, &err);

    check(err == 0 && d !is null && strcmp(d, "A <small> test") == 0, "gretl_get_gdt_description");
    free(d);
    destroy_dataset(dset);
    remove(fname);
}

import core.stdc.stdlib : malloc;
alias malloc_ = malloc;

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);
    test_serialize();
    test_gdt_roundtrip();
    test_description();
    printf("\n%d checks, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
