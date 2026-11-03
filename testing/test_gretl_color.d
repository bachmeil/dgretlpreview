/* Tests for gretl_color.d: channel decomposition and colormix_array() (without plotting). Build and run with
   testing/run_tests.sh. */

import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;
import core.stdc.string : strcmp;

import gretl.gretl_array : gretl_array_destroy, gretl_array_get_data;
import gretl.gretl_color;
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

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    ubyte a, r, g, b;
    decompose_argb(0x80ff7f01, &a, &r, &g, &b);
    check(a == 0x80 && r == 0xff && g == 0x7f && b == 0x01, "decompose_argb splits the channels");

    int err = 0;
    immutable double[3] f = [0.0, 0.5, 1.0];
    GretlArray * A = colormix_array(0xff0000, 0x0000ff, f.ptr, 3, 0, &err);
    check(err == 0 && A !is null, "colormix_array succeeds");
    auto s0 = cast(const(char) *) gretl_array_get_data(A, 0);
    auto s1 = cast(const(char) *) gretl_array_get_data(A, 1);
    auto s2 = cast(const(char) *) gretl_array_get_data(A, 2);
    check(strcmp(s0, "0xff0000") == 0, "f = 0 gives the first color");
    check(strcmp(s1, "0x80007f") == 0 || strcmp(s1, "0x800080") == 0, "f = 0.5 gives the midpoint");
    check(strcmp(s2, "0x0000ff") == 0, "f = 1 gives the second color");
    gretl_array_destroy(A);

    A = colormix_array(0x80ff0000, 0x800000ff, f.ptr, 1, 0, &err);
    check(strcmp(cast(const(char) *) gretl_array_get_data(A, 0), "0x80ff0000") == 0, "non-zero alpha gives 8 hex digits");
    gretl_array_destroy(A);

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
