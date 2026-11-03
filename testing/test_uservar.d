/* Smoke tests for uservar.d: adding scalars, matrices, strings and lists as named user variables, lookups by
   name, type and data pointer, scalar get/set (including the fast scalar accessors), replacing, deleting,
   counting by type and the private ("$name") variables.

   Not exercised: function-level scoping (the second hash table), shell variables, the GUI callbacks and the
   XML serialisation of the user variables.

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;
import core.stdc.stdlib : malloc;
import core.stdc.string : strcmp, strdup;

import gretl.gretl_list : gretl_list_new;
import gretl.gretl_matrix : gretl_matrix_alloc;
import gretl.gretl_matrixh : GretlMatrixRaw;
import gretl.libgretlh;
import gretl.uservar;
import gretl.uservarh : UserVar;

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

    check(gretl_is_scalar("a") == 0 && gretl_is_user_var("a") == 0, "nothing is defined at the start");

    /* scalars */
    check(gretl_scalar_add("a", 2.5) == 0, "gretl_scalar_add");
    int err = 0;
    check(gretl_is_scalar("a") == 1 && gretl_is_user_var("a") == 1 && gretl_scalar_get_value("a", &err) == 2.5 && err == 0, "the scalar can be looked up and read");
    check(gretl_scalar_set_value("a", 7.0) == 0 && gretl_scalar_get_value("a", &err) == 7.0, "gretl_scalar_set_value");
    gretl_scalar_get_value("nosuch", &err);
    check(err != 0, "reading an undefined scalar is an error");
    err = 0;

    UserVar * u = get_user_var_by_name("a");
    check(u !is null && user_var_get_type(u) == GretlType.double_ && strcmp(user_var_get_name(u), "a") == 0, "get_user_var_by_name gives a scalar named a");
    check(user_var_get_scalar_value(u) == 7.0 && user_var_set_scalar_value(u, 8.0) == 0 && gretl_scalar_get_value("a", &err) == 8.0, "scalar accessors on the variable");
    check(get_user_var_by_data(u.ptr) is u, "lookup by data pointer");

    /* matrix */
    GretlMatrixRaw * m = gretl_matrix_alloc(2, 2);
    m.val[0] = 1;
    m.val[1] = 2;
    m.val[2] = 3;
    m.val[3] = 4;
    check(user_var_add("M", GretlType.matrix, m) == 0 && user_var_get_type_by_name("M") == GretlType.matrix, "a matrix variable");
    check(get_user_var_of_type_by_name("M", GretlType.matrix) !is null && get_user_var_of_type_by_name("M", GretlType.double_) is null, "lookup by name and type");
    check(get_user_var_by_data(m) !is null, "the matrix is found by its pointer");

    /* string and list */
    check(user_var_add("s", GretlType.string_, strdup("hello")) == 0, "a string variable");
    const(char) * sv = get_string_by_name("s");
    check(sv !is null && strcmp(sv, "hello") == 0 && gretl_is_string("s") == 1, "the string can be read back");
    int * L = gretl_list_new(2);
    L[1] = 3;
    L[2] = 4;
    check(user_var_add("L", GretlType.list, L) == 0, "a list variable");
    int * Lb = get_list_by_name("L");
    check(Lb !is null && Lb[0] == 2 && Lb[2] == 4, "the list can be read back");

    check(user_var_count_for_type(GretlType.double_) == 1 && user_var_count_for_type(GretlType.matrix) == 1, "counts by type");

    /* duplicates, replace, private */
    check(gretl_scalar_add("a", 1.0) != 0 || gretl_scalar_get_value("a", &err) == 1.0, "adding an existing name is refused or replaces");
    double * pd = cast(double *) malloc(double.sizeof);
    *pd = 42.0;
    check(user_var_add_or_replace("a", GretlType.double_, pd) == 0 && gretl_scalar_get_value("a", &err) == 42.0, "add_or_replace changes the value");

    check(gretl_scalar_add("$tmp", 3.0) == 0 && gretl_is_scalar("$tmp") == 1, "a private variable can be added");

    /* delete */
    check(user_var_delete_by_name("a", null) == 0 && gretl_is_scalar("a") == 0, "delete_by_name removes the scalar");
    check(user_var_delete_by_name("a", null) != 0, "deleting a missing variable is an error");

    destroy_user_vars();
    check(gretl_is_user_var("M") == 0 && gretl_is_user_var("s") == 0, "destroy_user_vars empties the table");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
