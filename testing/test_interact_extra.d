/* Tests for interact_extra.d: check_stringvar_name (new, existing and wrongly typed names) and the plot-produced
   flag.

   Build and run with testing/run_tests.sh. */

import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import gretl.gretl_errorsh : GretlError;
import gretl.interact_extra;
import gretl.libgretlh;
import gretl.uservar : get_string_by_name, user_var_add;

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

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    check(check_stringvar_name("nosuchvar", 0, null) == GretlError.data, "missing variable is an error when new ones are not allowed");
    check(check_stringvar_name("outstr", 1, null) == 0 && get_string_by_name("outstr") !is null, "a new string variable is created on request");
    check(check_stringvar_name("outstr", 0, null) == 0, "an existing string variable is accepted");
    check(check_stringvar_name("1bad", 1, null) != 0, "an invalid identifier is refused");

    double x = 1.0;
    user_var_add("xs", GretlType.double_, &x);
    check(check_stringvar_name("xs", 1, null) == GretlError.types, "a scalar is not a string variable");

    check(plot_was_produced() == 0, "flag starts clear");
    set_plot_produced();
    check(plot_was_produced() == 1, "flag set");
    clear_plot_produced();
    check(plot_was_produced() == 0, "flag cleared");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
