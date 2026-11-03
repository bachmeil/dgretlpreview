/* Tests for cephes/constants.d: the constants have the C values, are consistent with each other and with the
   machine epsilon, and are used by the cephes functions (igami, which reads MACHEP, MAXNUM and MAXLOG, is called
   from here).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs, signbit, sqrt, log;
import core.stdc.float_ : DBL_EPSILON, DBL_MAX, DBL_MIN;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import cephes.constants;
import cephes.igami : igami;

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

bool close(double a, double b, double tol = 1e-15)
{
    return fabs(a - b) <= tol * (1.0 + fabs(b));
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    check(MACHEP == DBL_EPSILON, "MACHEP is the double epsilon");
    check(MAXNUM == DBL_MAX, "MAXNUM is the largest double");
    /* UFLOWTHRESH from the C is DBL_MIN plus one unit in the last place (bugs.md #101); nothing uses it */
    check(close(UFLOWTHRESH, DBL_MIN, 1e-15) && UFLOWTHRESH > DBL_MIN, "UFLOWTHRESH is the C value, one ulp above the smallest normal double");
    check(close(MAXLOG, 708.396418532264106224) && MINLOG == -MAXLOG, "MAXLOG and MINLOG are +-708.396...");
    check(close(MAXLOG, log(DBL_MIN) * -1.0, 1e-12), "MAXLOG is -log of the smallest normal double");

    check(close(PI, 3.141592653589793) && PIO2 * 2 == PI && PIO4 * 4 == PI, "PI, PIO2 and PIO4 agree");
    check(close(THPIO4, 3 * PI / 4) && close(TWOOPI, 2 / PI), "THPIO4 and TWOOPI");
    check(close(SQRT2 * SQRT2, 2.0, 1e-15) && close(SQRTH * 2, SQRT2) && close(SQ2OPI, sqrt(2 / PI)), "SQRT2, SQRTH and SQ2OPI");
    check(close(LOGE2, log(2.0)) && close(LOGSQ2 * 2, LOGE2) && close(LOG2E * LOGE2, 1.0), "LOGE2, LOGSQ2 and LOG2E");
    check(NEGZERO == 0.0 && signbit(NEGZERO) != 0, "NEGZERO is minus zero");

    /* igami reads MACHEP, MAXNUM and MAXLOG: igami(1, y) solves 1 - exp(-x) = 1 - y, so igami(1, exp(-2)) = 2 */
    check(close(igami(1.0, 0.1353352832366127), 2.0, 1e-10), "igami works with the constants");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
