/* Tests for plugin/johansen.d: trace_pvalue, the p-value of Johansen's trace test by Doornik's gamma approximation (the only part of
   johansen.c ported so far). 80 p-values
   are compared with an independent computation (a Python program that parses the same tables from gretl's johansen.c and evaluates the
   formula with its own incomplete gamma function: agreement to a relative 1e-9); the published 5 percent critical values of the
   trace test give p-values near 0.05; and the plugin look-up and johansen_trace_pval (pvalues.d) are checked.
   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import std.conv : to;
import std.math : isNaN;

import gretl.plugin.johansen;
import gretl.plugins : get_plugin_function;
import gretl.pvalues : johansen_trace_pval;

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

bool close(double a, double b, double tol)
{
    return fabs(a - b) <= tol * (fabs(b) > 1e-300 ? fabs(b) : 1.0);
}

struct Ref {
    double trace;
    int n, det, T;
    double pval;
}

/* from the independent Python computation: trace, n, det, T, p-value */
immutable Ref[] refs = [
    Ref(8, 1, 0, 0, 0.0049024176410796426),
    Ref(25, 1, 0, 0, 3.5397635078355089e-07),
    Ref(8, 1, 0, 60, 0.0052837939994625658),
    Ref(25, 1, 0, 60, 4.8398231618807399e-07),
    Ref(8, 2, 0, 0, 0.23822729898142866),
    Ref(25, 2, 0, 0, 0.00018137204663110697),
    Ref(8, 2, 0, 60, 0.25088337048116799),
    Ref(25, 2, 0, 60, 0.00020411950182701769),
    Ref(8, 4, 0, 0, 0.99998716874327109),
    Ref(25, 4, 0, 0, 0.65103074285442974),
    Ref(8, 4, 0, 60, 0.99999124146797302),
    Ref(25, 4, 0, 60, 0.69859120804529029),
    Ref(8, 7, 0, 0, 1),
    Ref(25, 7, 0, 0, 0.99999999999999956),
    Ref(8, 7, 0, 60, 1),
    Ref(25, 7, 0, 60, 1),
    Ref(8, 1, 1, 0, 0.083815846026574076),
    Ref(25, 1, 1, 0, 1.6029215144635856e-05),
    Ref(8, 1, 1, 60, 0.084387928685340716),
    Ref(25, 1, 1, 60, 1.6735750488757384e-05),
    Ref(8, 2, 1, 0, 0.82160511311528373),
    Ref(25, 2, 1, 0, 0.0088941327616109105),
    Ref(8, 2, 1, 60, 0.8309781723904387),
    Ref(25, 2, 1, 60, 0.010522905064608623),
    Ref(8, 4, 1, 0, 0.9999999998863075),
    Ref(25, 4, 1, 0, 0.98349931280501579),
    Ref(8, 4, 1, 60, 0.9999999999329946),
    Ref(25, 4, 1, 60, 0.98849165997715494),
    Ref(8, 7, 1, 0, 1),
    Ref(25, 7, 1, 0, 1),
    Ref(8, 7, 1, 60, 1),
    Ref(25, 7, 1, 60, 1),
    Ref(8, 1, 2, 0, 0.0046777349810472636),
    Ref(25, 1, 2, 0, 5.7330314375838793e-07),
    Ref(8, 1, 2, 60, 0.0056353423871736957),
    Ref(25, 1, 2, 60, 9.8535310727783693e-07),
    Ref(8, 2, 2, 0, 0.47265649112752439),
    Ref(25, 2, 2, 0, 0.0010503969361228944),
    Ref(8, 2, 2, 60, 0.48819052542936248),
    Ref(25, 2, 2, 60, 0.0012900738573974372),
    Ref(8, 4, 2, 0, 0.99999997142397379),
    Ref(25, 4, 2, 0, 0.91473082397354966),
    Ref(8, 4, 2, 60, 0.99999998345198726),
    Ref(25, 4, 2, 60, 0.93582255659860714),
    Ref(8, 7, 2, 0, 1),
    Ref(25, 7, 2, 0, 1),
    Ref(8, 7, 2, 60, 1),
    Ref(25, 7, 2, 60, 1),
    Ref(8, 1, 3, 0, 0.25929023784799793),
    Ref(25, 1, 3, 0, 0.00015917828928989415),
    Ref(8, 1, 3, 60, 0.2601284372700583),
    Ref(25, 1, 3, 60, 0.00014912042426179373),
    Ref(8, 2, 3, 0, 0.97653213639169256),
    Ref(25, 2, 3, 0, 0.062142243675612871),
    Ref(8, 2, 3, 60, 0.97895923547857266),
    Ref(25, 2, 3, 60, 0.06989626176434037),
    Ref(8, 4, 3, 0, 0.99999999999999745),
    Ref(25, 4, 3, 0, 0.99957565342719434),
    Ref(8, 4, 3, 60, 0.99999999999999889),
    Ref(25, 4, 3, 60, 0.99976236321428802),
    Ref(8, 7, 3, 0, 1),
    Ref(25, 7, 3, 0, 1),
    Ref(8, 7, 3, 60, 1),
    Ref(25, 7, 3, 60, 1),
    Ref(8, 1, 4, 0, 0.004677734981047268),
    Ref(25, 1, 4, 0, 5.7330314375838634e-07),
    Ref(8, 1, 4, 60, 0.0055563661083403767),
    Ref(25, 1, 4, 60, 9.5352868041683563e-07),
    Ref(8, 2, 4, 0, 0.68728385746467757),
    Ref(25, 2, 4, 0, 0.0043013910747010673),
    Ref(8, 2, 4, 60, 0.68875098704829185),
    Ref(25, 2, 4, 60, 0.004424684177777442),
    Ref(8, 4, 4, 0, 0.99999999995639055),
    Ref(25, 4, 4, 0, 0.98808924657551411),
    Ref(8, 4, 4, 60, 0.99999999996540068),
    Ref(25, 4, 4, 60, 0.9900629408679712),
    Ref(8, 7, 4, 0, 1),
    Ref(25, 7, 4, 0, 1),
    Ref(8, 7, 4, 60, 1),
    Ref(25, 7, 4, 60, 1),
];

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    /* the independent computation */
    bool refOk = true;
    foreach (r; refs) {
        const p = trace_pvalue(r.trace, r.n, r.det, r.T);

        if (!close(p, r.pval, 1e-9)) {
            printf("  mismatch: trace %g n %d det %d T %d: %.17g, reference %.17g\n", r.trace, r.n, r.det, r.T, p, r.pval);
            refOk = false;
        }
    }
    check(refOk, "trace_pvalue agrees with the independent computation on 80 values");

    /* published asymptotic 5 percent critical values of the trace test (Osterwald-Lenum 1992; Doornik's approximation is close) */
    check(fabs(trace_pvalue(4.1, 1, 0, 0) - 0.05) < 2e-3 && fabs(trace_pvalue(12.3, 2, 0, 0) - 0.05) < 2e-3, "no constant: 4.1 (n = 1) and 12.3 (n = 2) at 5 percent");
    check(fabs(trace_pvalue(15.41, 2, 2, 0) - 0.05) < 2e-3 && fabs(trace_pvalue(29.68, 3, 2, 0) - 0.05) < 4e-3, "unrestricted constant: 15.41 (n = 2) and 29.68 (n = 3) at 5 percent");

    /* the shape of the function */
    check(trace_pvalue(10, 3, 2, 0) > trace_pvalue(20, 3, 2, 0) && trace_pvalue(20, 3, 2, 0) > trace_pvalue(40, 3, 2, 0), "p decreases as the statistic grows");
    check(trace_pvalue(25, 3, 2, 50) < trace_pvalue(25, 3, 2, 0) + 0.2 && trace_pvalue(25, 3, 2, 50) != trace_pvalue(25, 3, 2, 0), "a finite T changes the result");
    check(trace_pvalue(25, 3, 2, 10000) == trace_pvalue(25, 3, 2, 0) && trace_pvalue(25, 3, 2, 20000) == trace_pvalue(25, 3, 2, 0), "T of 10000 or more is the asymptotic result");
    check(trace_pvalue(0.5, 4, 3, 0) > 0.99 && trace_pvalue(500, 4, 3, 0) < 1e-10, "extreme statistics give p near 1 and near 0");
    check(isNaN(trace_pvalue(10, 0, 2, 0)) && isNaN(trace_pvalue(10, -1, 2, 0)) && isNaN(trace_pvalue(10, 3, -1, 0)) && isNaN(trace_pvalue(10, 3, 5, 0)), "n < 1 or det outside 0 .. 4 gives NA");

    /* the library finds it */
    check(get_plugin_function("trace_pvalue") !is null, "get_plugin_function finds trace_pvalue");
    check(johansen_trace_pval(3, 2, 0, 29.68) == trace_pvalue(29.68, 3, 2, 0) && !isNaN(johansen_trace_pval(2, 1, 100, 20.0)), "johansen_trace_pval (pvalues.d) uses it");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
