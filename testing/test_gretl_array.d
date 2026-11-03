/* Smoke tests for gretl_array.d: arrays of strings (creation, setting, appending, copying, sorting, joining,
   union, intersection, concatenation, flattening, position and membership, dropping and deleting) and arrays of
   matrices (the common-block allocation, appending, element access), plus equality and type queries.

   Not exercised: arrays of bundles, lists and arrays, sorting with a user function, serialisation, printing and
   the stringify and subspec functions.

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;
import core.stdc.stdlib : free;
import core.stdc.string : strcmp;

import gretl.gretl_array;
import gretl.gretl_matrix : gretl_matrix_alloc;
import gretl.gretl_matrixh : GretlMatrixRaw;
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

bool iseq(GretlArray * a, int i, string s)
{
    import std.string : toStringz;
    int ns;
    char ** S = gretl_array_get_strings(a, &ns);

    return S !is null && i < ns && S[i] !is null && strcmp(S[i], toStringz(s)) == 0;
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    int err = 0;
    GretlArray * a = gretl_array_new(GretlType.strings, 3, &err);

    check(err == 0 && a !is null && gretl_array_get_length(a) == 3 && gretl_array_get_type(a) == GretlType.strings, "a new strings array of length 3");
    check(gretl_array_get_content_type(a) == GretlType.string_, "content type is string");

    gretl_array_set_string(a, 0, cast(char *) "pear".ptr, 1);
    gretl_array_set_string(a, 1, cast(char *) "apple".ptr, 1);
    gretl_array_set_string(a, 2, cast(char *) "fig".ptr, 1);
    check(iseq(a, 0, "pear") && iseq(a, 1, "apple") && iseq(a, 2, "fig"), "strings set by copy");

    gretl_array_append_string(a, cast(char *) "kiwi".ptr, 1);
    check(gretl_array_get_length(a) == 4 && iseq(a, 3, "kiwi"), "append_string");

    check(gretl_strings_array_includes(a, "fig") == 1 && gretl_strings_array_includes(a, "plum") == 0, "strings_array_includes");

    GretlMatrixRaw * pos = gretl_strings_array_pos(a, "apple", &err);
    check(err == 0 && pos !is null && pos.val[0] == 2.0, "strings_array_pos is 1-based");
    free(pos);

    /* sorting */
    GretlArray * s = gretl_strings_sort(a, 0, &err);
    check(err == 0 && iseq(s, 0, "apple") && iseq(s, 1, "fig") && iseq(s, 2, "kiwi") && iseq(s, 3, "pear"), "ascending sort");
    GretlArray * sd = gretl_strings_sort(a, 1, &err);
    check(iseq(sd, 0, "pear") && iseq(sd, 3, "apple"), "descending sort");
    gretl_array_destroy(sd);

    /* copy and equality */
    GretlArray * c = gretl_array_copy(a, &err);
    check(err == 0 && gretl_arrays_are_equal(a, c) == 1, "a copy equals the original");
    gretl_array_set_string(c, 0, cast(char *) "plum".ptr, 1);
    check(gretl_arrays_are_equal(a, c) == 0 && iseq(a, 0, "pear"), "changing the copy does not change the original");
    gretl_array_destroy(c);

    /* set operations */
    /* explicit terminators: LDC 1.43 on Windows packs the literals of this initializer without their NULs */
    char *[3] sb = [cast(char *) "fig\0".ptr, cast(char *) "date\0".ptr, cast(char *) "apple\0".ptr];
    GretlArray * b = gretl_array_from_strings(sb.ptr, 3, 1, &err);

    check(err == 0 && gretl_array_get_length(b) == 3 && iseq(b, 1, "date"), "array from strings");

    GretlArray * j = gretl_arrays_join(a, b, &err);
    check(err == 0 && gretl_array_get_length(j) == 7, "join keeps all elements");
    gretl_array_destroy(j);

    GretlArray * u = gretl_arrays_union(a, b, &err);
    check(err == 0 && gretl_array_get_length(u) == 5, "union has five distinct strings");
    gretl_array_destroy(u);

    GretlArray * it = gretl_arrays_intersection(a, b, &err);
    check(err == 0 && gretl_array_get_length(it) == 2, "intersection has the two common strings");
    gretl_array_destroy(it);

    GretlArray * cc = gretl_arrays_concat(a, b, &err);   /* element-wise concatenation needs equal lengths */
    check(err != 0 || cc is null, "concat of arrays of different lengths is an error");
    err = 0;
    if (cc !is null) gretl_array_destroy(cc);

    char * flat = gretl_strings_array_flatten(b, ",", &err);
    check(err == 0 && flat !is null && strcmp(flat, "fig,date,apple") == 0, "flatten with a separator");
    free(flat);

    /* dropping and deleting */
    check(gretl_array_drop_string(a, "fig") == 0 && gretl_array_get_length(a) == 3 && !gretl_strings_array_includes(a, "fig"), "drop_string");
    check(gretl_array_delete_element(a, 0) == 0 && gretl_array_get_length(a) == 2 && iseq(a, 0, "apple"), "delete_element shifts the rest down");
    gretl_array_destroy(a);
    gretl_array_destroy(b);
    gretl_array_destroy(s);

    /* matrices */
    GretlArray * M = gretl_matrix_array_sized(3, 2, 2, &err);

    check(err == 0 && M !is null && gretl_array_get_length(M) == 3 && gretl_array_get_type(M) == GretlType.matrices, "matrix array of three 2 x 2 matrices");
    GretlMatrixRaw * m1 = cast(GretlMatrixRaw *) gretl_array_get_data(M, 1);
    check(m1 !is null && m1.rows == 2 && m1.cols == 2, "elements have the requested size");
    m1.val[3] = 7.5;
    GretlType ty;
    GretlMatrixRaw * m1b = cast(GretlMatrixRaw *) gretl_array_get_element(M, 1, &ty, &err);
    check(err == 0 && ty == GretlType.matrix && m1b.val[3] == 7.5, "get_element returns the same matrix");
    GretlMatrixRaw * extra = gretl_matrix_alloc(1, 1);
    extra.val[0] = 3.0;
    check(gretl_array_append_matrix(M, extra, 0) == 0 && gretl_array_get_length(M) == 4, "append_matrix");
    gretl_array_destroy(M);

    GretlArray * sing = gretl_singleton_array(cast(void *) "solo".ptr, GretlType.strings, 1, &err);
    check(err == 0 && sing !is null && gretl_array_get_length(sing) == 1 && iseq(sing, 0, "solo"), "singleton array");
    gretl_array_destroy(sing);

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
