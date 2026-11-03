/* Tests for gretl_typemap.d: singular/plural and plain/reference type relations, type names and their parsing,
   the "gen type" keywords, the type ordering and predicates, the mismatch test and the name-conflict message.

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;
import core.stdc.stdlib : free;
import core.stdc.string : strcmp;

import gretl.gretl_typemap;
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

    check(gretl_type_get_plural(GretlType.matrix) == GretlType.matrices && gretl_type_get_plural(GretlType.string_) == GretlType.strings, "plural types");
    check(gretl_type_get_plural(GretlType.double_) == GretlType.none && gretl_type_get_plural(GretlType.none) == GretlType.none, "scalars and null have no plural");
    check(gretl_type_get_singular(GretlType.bundles) == GretlType.bundle && gretl_type_get_singular(GretlType.arrays) == GretlType.array, "singular types");
    check(gretl_type_get_ref_type(GretlType.list) == GretlType.list_ref && gretl_type_get_ref_type(GretlType.lists) == GretlType.lists_ref, "reference types of plain and plural types");
    check(gretl_type_get_ref_type(GretlType.matrix_ref) == GretlType.matrix_ref, "a reference type passes through");
    check(gretl_type_get_plain_type(GretlType.series_ref) == GretlType.series && gretl_type_get_plain_type(GretlType.strings_ref) == GretlType.strings, "plain types of references");
    check(gretl_type_get_plain_type(GretlType.matrix) == GretlType.matrix, "a plain type passes through");

    check(gretl_is_array_type(GretlType.matrices) == 1 && gretl_is_array_type(GretlType.matrix) == 0, "gretl_is_array_type");
    check(gretl_is_array_ref_type(GretlType.bundles_ref) == 1 && gretl_is_array_ref_type(GretlType.bundles) == 0, "gretl_is_array_ref_type");
    check(gretl_is_arrayable_type(GretlType.matrix) == 1 && gretl_is_arrayable_type(GretlType.series) == 0, "arrayable types");

    check(strcmp(gretl_type_get_name(GretlType.double_), "scalar") == 0 && strcmp(gretl_type_get_name(GretlType.matrix_ref), "matrix *") == 0
          && strcmp(gretl_type_get_name(GretlType.none), "null") == 0, "type names");
    check(strcmp(gretl_type_get_name(cast(GretlType) 9999), "invalid") == 0, "an unknown type is invalid");

    check(gretl_type_from_string("matrix") == GretlType.matrix && gretl_type_from_string("matrix *") == GretlType.matrix_ref
          && gretl_type_from_string("matrixref") == GretlType.matrix_ref, "parsing matrix and its reference forms");
    check(gretl_type_from_string("bundles") == GretlType.bundles && gretl_type_from_string("strings *") == GretlType.strings_ref, "parsing plural forms");
    check(gretl_type_from_string("boolean") == GretlType.bool_ && gretl_type_from_string("unsigned64") == GretlType.uint64 && gretl_type_from_string("numeric") == GretlType.numeric, "parsing aliases");
    check(gretl_type_from_string("bogus") == GretlType.none && gretl_type_from_string("matrixx") == GretlType.none, "unknown names give none");

    check(gretl_get_gen_type("series") == GretlType.series && gretl_get_gen_type("matrices") == GretlType.matrices && gretl_get_gen_type("nosuch") == GretlType.none, "gen type keywords");
    gretl_typemap_cleanup();
    check(gretl_get_gen_type("list") == GretlType.list, "the keyword table is rebuilt after clean-up");

    check(gretl_type_get_order(GretlType.double_) == 1 && gretl_type_get_order(GretlType.series) == 2 && gretl_type_get_order(GretlType.list) == 7
          && gretl_type_get_order(GretlType.void_) == 0, "typeof() ordering");
    check(gretl_is_scalar_type(GretlType.int_) == 1 && gretl_is_scalar_type(GretlType.matrix) == 0, "gretl_is_scalar_type");
    check(gretl_is_series_type(GretlType.useries) == 1 && gretl_is_series_type(GretlType.list) == 0, "gretl_is_series_type");
    check(gretl_type_mismatch(GretlType.matrix, GretlType.matrix_ref) == 0 && gretl_type_mismatch(GretlType.matrix, GretlType.list) == 1, "type mismatch ignores reference-ness");

    char * msg = name_conflict_message("x", GretlType.matrix);
    check(msg !is null && strcmp(msg, "A matrix named x already exists") == 0, "name-conflict message");
    free(msg);
    check(name_conflict_message("x", GretlType.void_) is null, "no message for other types");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
