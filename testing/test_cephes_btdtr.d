/* Tests for cephes/btdtr.d: the beta distribution function: equal to incbet, closed forms and the domain case.

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs, pow;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import cephes.btdtr;
import cephes.incbet : incbet;

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

    bool same = true;
    foreach (a; [0.5, 1.0, 3.2, 20.0]) foreach (b; [0.7, 2.0, 9.5]) foreach (x; [0.0, 0.05, 0.4, 0.9, 1.0]) same = same && btdtr(a, b, x) == incbet(a, b, x);
    check(same, "btdtr is incbet");
    check(fabs(btdtr(1, 1, 0.3) - 0.3) < 1e-15 && fabs(btdtr(2, 1, 0.3) - 0.09) < 1e-15 && fabs(btdtr(1, 2, 0.3) - (1 - 0.49)) < 1e-15, "uniform, x^2 and 1 - (1-x)^2 cases");
    check(btdtr(-1, 2, 0.5) == 0.0, "a <= 0 gives 0 (with an error message)");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
