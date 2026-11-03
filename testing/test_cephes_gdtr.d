/* Tests for cephes/gdtr.d: the gamma distribution function gdtr and its upper tail gdtrc, against closed forms
   (shape 1 and 2), the relation to the chi-square distribution (rate 1/2, shape df/2), the complement, and the
   domain error. A bit-for-bit comparison with the original C code over 1,152 combinations was run once when the
   module was ported (see log.md).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : exp, fabs;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import cephes.chdtr : chdtr, chdtrc;
import cephes.gdtr;

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

bool close(double a, double b, double tol = 1e-12)
{
    return fabs(a - b) <= tol * (fabs(b) + 1e-300);
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    /* shape 1 is the exponential distribution, shape 2 the Erlang-2 */
    bool ok = true;
    foreach (a; [0.3, 1.0, 4.0]) {
        foreach (x; [0.05, 0.7, 2.0, 6.0]) {
            ok = ok && close(gdtr(a, 1, x), 1 - exp(-a * x)) && close(gdtrc(a, 1, x), exp(-a * x));
            ok = ok && close(gdtrc(a, 2, x), exp(-a * x) * (1 + a * x));
        }
    }
    check(ok, "shape 1 and shape 2 match their closed forms");

    /* chi-square with df degrees of freedom is the gamma distribution with rate 1/2 and shape df/2 */
    ok = true;
    foreach (df; [1.0, 3.0, 8.0, 25.0]) {
        foreach (x; [0.5, 3.0, 10.0, 40.0]) {
            ok = ok && close(gdtr(0.5, df / 2, x), chdtr(df, x), 1e-13) && close(gdtrc(0.5, df / 2, x), chdtrc(df, x), 1e-13);
        }
    }
    check(ok, "gdtr(1/2, df/2, x) = chdtr(df, x)");

    check(close(gdtr(2.5, 3.5, 1.7) + gdtrc(2.5, 3.5, 1.7), 1.0, 1e-14) && gdtr(2, 3, 0.0) == 0.0 && gdtrc(2, 3, 0.0) == 1.0, "complement and the limits at 0");

    check(gdtr(2, 3, -1.0) == 0.0 && gdtrc(2, 3, -1.0) == 0.0, "negative x gives 0 (with an error message)");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
