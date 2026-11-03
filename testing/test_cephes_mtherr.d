/* Tests for cephes/mtherr.d: the messages written to stderr (captured through a temporary file), the recorded
   error code and its clearing by get_cephes_errno(), out-of-range codes, and set_cephes_hush(). The C
   comparison of the original mtherr.c (the same messages and codes for all codes and for both entry points) was
   run once when the module was ported (see log.md).

   Build and run with testing/run_tests.sh. */

import core.stdc.stdio : FILE, _IOLBF, fflush, fileno, fread, printf, rewind, setvbuf, stderr, stdout, tmpfile;
import core.stdc.string : memset, strcmp;
version (Posix) import core.sys.posix.unistd : dup, dup2, close; else import gretl.winport : dup, dup2, close;

import cephes.libprobh : CephesError;
import cephes.mtherr;

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

/* run dg with stderr redirected to a temporary file and return what was written */
string capture(void delegate() dg)
{
    FILE * tmp = tmpfile();
    fflush(stderr);
    int saved = dup(2);
    dup2(fileno(tmp), 2);
    dg();
    fflush(stderr);
    dup2(saved, 2);
    close(saved);
    char[512] buf;
    memset(buf.ptr, 0, buf.length);
    rewind(tmp);
    size_t n = fread(buf.ptr, 1, buf.length - 1, tmp);
    return buf[0 .. n].idup;
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    get_cephes_errno();   // clear any earlier state

    /* messages for each code */
    check(capture({ mtherr("fn", CephesError.domain); }) == "fn domain error\n" && get_cephes_errno() == CephesError.domain, "domain error message and code");
    check(capture({ mtherr("fn", CephesError.sing); }) == "fn singularity error\n" && get_cephes_errno() == CephesError.sing, "singularity");
    check(capture({ mtherr("fn", CephesError.overflow); }) == "fn overflow error\n", "overflow");
    check(capture({ mtherr("fn", CephesError.underflow); }) == "fn underflow error\n", "underflow");
    check(capture({ mtherr("fn", CephesError.tloss); }) == "fn total loss of precision error\n", "total loss of precision");
    check(capture({ mtherr("fn", CephesError.ploss); }) == "fn partial loss of precision error\n", "partial loss of precision");
    check(capture({ mtherr("fn", CephesError.unknown); }) == "fn unknown error\n", "unknown");
    get_cephes_errno();

    /* with an argument */
    check(capture({ mtherr_with_arg("lgam", CephesError.sing, -5.0); }) == "lgam singularity error (arg = -5)\n" && get_cephes_errno() == CephesError.sing, "message with the offending argument");
    check(capture({ mtherr_with_arg("ndtri", CephesError.domain, 2.5); }) == "ndtri domain error (arg = 2.5)\n", "argument printed with %g");

    /* out-of-range codes are reported as unknown */
    check(capture({ mtherr("fn", 0); }) == "fn unknown error\n" && get_cephes_errno() == CephesError.unknown, "code 0 is unknown");
    check(capture({ mtherr("fn", 99); }) == "fn unknown error\n" && get_cephes_errno() == CephesError.unknown, "code 99 is unknown");
    check(capture({ mtherr("fn", -3); }) == "fn unknown error\n", "negative code is unknown");

    /* the code is cleared by reading it, and holds the last error */
    capture({ mtherr("a", CephesError.domain); mtherr("b", CephesError.overflow); });
    check(get_cephes_errno() == CephesError.overflow && get_cephes_errno() == 0, "the last code is returned once, then 0");

    /* hush: no output and no recording */
    set_cephes_hush(1);
    string quiet = capture({ mtherr("fn", CephesError.domain); mtherr_with_arg("fn", CephesError.sing, 1.0); });
    int quietcode = get_cephes_errno();
    set_cephes_hush(0);
    check(quiet == "" && quietcode == 0, "set_cephes_hush(1) silences messages and the recorded code");
    check(capture({ mtherr("fn", CephesError.domain); }) == "fn domain error\n" && get_cephes_errno() == CephesError.domain, "set_cephes_hush(0) restores them");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
