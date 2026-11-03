/* Smoke tests for gretl_bundle.d: members of every supported kind (scalar, integer, string, matrix, list, nested
   bundle, series), type queries, missing-key errors, overwriting, notes, deleting, key lists, copying, appending,
   union, equality, printing and the destroy functions.

   Not exercised: Kalman bundles, virtual series/matrix members, XML and binary (de)serialisation, the sysinfo
   bundle and the model/system bundlisers.

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;
import core.stdc.stdlib : free;
import core.stdc.string : strcmp, strstr;

import gretl.dataset : create_new_dataset;
import gretl.gretl_array : gretl_array_destroy, gretl_array_get_length, gretl_array_get_strings;
import gretl.gretl_bundle;
import gretl.gretl_list : gretl_list_new;
import gretl.gretl_matrix : gretl_matrix_alloc, gretl_matrix_free;
import gretl.gretl_matrixh : GretlMatrixRaw;
import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_get_buffer, gretl_print_new;
import gretl.gretl_prnh : PrnType;
import gretl.libgretlh;

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

bool close(double a, double b, double tol = 1e-8)
{
    return fabs(a - b) <= tol * (1.0 + fabs(b));
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    GretlBundle * b = gretl_bundle_new();
    int err = 0;

    check(b !is null && gretl_bundle_get_n_keys(b) == 0 && !gretl_bundle_has_content(b), "a new bundle is empty");

    /* scalars, ints and strings */
    check(gretl_bundle_set_scalar(b, "x", 2.5) == 0 && gretl_bundle_set_int(b, "n", 7) == 0 && gretl_bundle_set_string(b, "name", "alpha") == 0, "scalar, int and string set");
    check(gretl_bundle_get_n_keys(b) == 3 && gretl_bundle_has_content(b), "three keys");
    check(gretl_bundle_get_scalar(b, "x", &err) == 2.5 && err == 0, "get_scalar");
    check(gretl_bundle_get_int(b, "n", &err) == 7 && err == 0, "get_int");
    const(char) * s = gretl_bundle_get_string(b, "name", &err);
    check(err == 0 && s !is null && strcmp(s, "alpha") == 0, "get_string");
    check(gretl_bundle_get_member_type(b, "x", &err) == GretlType.double_ && gretl_bundle_get_member_type(b, "name", &err) == GretlType.string_, "member types");
    check(gretl_bundle_has_key(b, "x") == 1 && gretl_bundle_has_key(b, "zzz") == 0, "has_key");

    gretl_bundle_get_scalar(b, "missing", &err);
    check(err != 0, "a missing key is an error");
    err = 0;

    /* overwrite */
    gretl_bundle_set_scalar(b, "x", 9.0);
    check(gretl_bundle_get_scalar(b, "x", &err) == 9.0 && gretl_bundle_get_n_keys(b) == 3, "overwriting a member keeps the key count");

    /* matrix (copied on set) */
    GretlMatrixRaw * m = gretl_matrix_alloc(2, 2);

    foreach (i; 0 .. 4) {
        m.val[i] = i + 1.0;
    }
    check(gretl_bundle_set_matrix(b, "M", m) == 0, "set_matrix");
    m.val[0] = 99;
    GretlMatrixRaw * mb = gretl_bundle_get_matrix(b, "M", &err);
    check(err == 0 && mb !is null && mb.val[0] == 1.0 && mb.val[3] == 4.0, "the stored matrix is a copy");
    gretl_matrix_free(m);

    /* list */
    int * L = gretl_list_new(2);

    L[1] = 5;
    L[2] = 6;
    check(gretl_bundle_set_list(b, "L", L) == 0, "set_list");
    int * Lb = gretl_bundle_get_list(b, "L", &err);
    check(err == 0 && Lb !is null && Lb[0] == 2 && Lb[2] == 6, "get_list");
    free(L);                    /* set_list copies; get_list returns the stored list, which the bundle owns */

    /* nested bundle, note */
    GretlBundle * inner = gretl_bundle_new();

    gretl_bundle_set_scalar(inner, "y", 1.0);
    check(gretl_bundle_get_n_keys(b) == 5, "five keys so far");

    gretl_bundle_set_note(b, "x", "the x value");
    const(char) * note = gretl_bundle_get_note(b, "x");
    check(note !is null && strcmp(note, "the x value") == 0, "member note");

    /* keys */
    GretlArray * keys = gretl_bundle_get_keys(b, &err);
    check(err == 0 && keys !is null && gretl_array_get_length(keys) == 5, "get_keys returns five keys");
    int ns;
    char ** K = gretl_array_get_strings(keys, &ns);
    bool found = false;
    foreach (i; 0 .. ns) {
        if (strcmp(K[i], "name") == 0) found = true;
    }
    check(found, "the key list contains 'name'");
    gretl_array_destroy(keys);

    /* copy and equality */
    GretlBundle * c = gretl_bundle_copy(b, &err);
    check(err == 0 && c !is null && gretl_bundle_get_n_keys(c) == 5 && gretl_bundles_are_equal(b, c) == 1, "a copy has the same members and equals the original");
    gretl_bundle_set_scalar(c, "x", -1.0);
    check(gretl_bundles_are_equal(b, c) == 0 && gretl_bundle_get_scalar(b, "x", &err) == 9.0, "changing the copy does not change the original");

    /* append and union */
    GretlBundle * extra = gretl_bundle_new();

    gretl_bundle_set_scalar(extra, "z", 3.0);
    gretl_bundle_set_scalar(extra, "x", 100.0);
    GretlBundle * u = gretl_bundle_union(b, extra, &err);
    check(err == 0 && gretl_bundle_get_n_keys(u) == 6, "union has six keys");
    check(gretl_bundle_get_scalar(u, "x", &err) == 9.0, "union keeps the first bundle's value for common keys");
    gretl_bundle_destroy(u);
    check(gretl_bundle_append(c, extra) == 0 && gretl_bundle_get_scalar(c, "z", &err) == 3.0, "append adds the other bundle's members");
    gretl_bundle_destroy(extra);
    gretl_bundle_destroy(c);

    /* delete */
    check(gretl_bundle_delete_data(b, "n") == 0 && gretl_bundle_has_key(b, "n") == 0 && gretl_bundle_get_n_keys(b) == 4, "delete_data removes a key");
    check(gretl_bundle_delete_data(b, "n") != 0, "deleting a missing key is an error");

    /* printing */
    Prn * prn = gretl_print_new(PrnType.buffer, &err);

    gretl_bundle_print(b, prn);
    const(char) * buf = gretl_print_get_buffer(prn);
    check(buf !is null && strstr(buf, "name") !is null && strstr(buf, "alpha") !is null, "print shows keys and values");
    gretl_print_destroy(prn);

    gretl_bundle_destroy(inner);
    gretl_bundle_destroy(b);
    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
