/* Tests for plugin/urcdist.d: MacKinnon p-values for unit-root tests and the DF-GLS p-values. The tables compiled into the
   program (gretl.plugin.urcdata) are compared bit for bit with the binary files they were generated from (gretl/data); the
   p-values are compared with the original C code (urcdist.c's fpval, gls, cholx and eval_all_crit compiled verbatim and fed
   from urcdata.bin) on a grid of 2720 values, 44 of which are listed here (agreement to a relative 1e-12: the two use a
   different normal distribution function); and the published asymptotic critical values give p-values near 1, 5, 10 percent.
   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import std.file : read;
import std.math : isNaN;

import gretl.adf_kpss : get_urc_pvalue;
import gretl.plugin.urcdata : embedded_matrices, urcdata;
import gretl.plugin.urcdist : dfgls_pvalue, mackinnon_pvalue;
import gretl.plugins : get_plugin_function;

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
    return fabs(a - b) <= tol * (fabs(b) > 1e-300 ? fabs(b) : 1.0);
}

struct Ref {
    int niv, itv, T;
    double tau, pval;
}

/* from the C reference program (niv, itv, T, tau, p-value) */
immutable Ref[] refs = [
    Ref(1, 1, 0, -3, 0.0026319074981709365),
    Ref(1, 2, 0, -3, 0.034908423328518366),
    Ref(1, 3, 0, -3, 0.13210496961829848),
    Ref(1, 4, 0, -3, 0.29615320493385489),
    Ref(2, 1, 0, -3, 0.027125736444118635),
    Ref(2, 2, 0, -3, 0.11025743917844427),
    Ref(2, 3, 0, -3, 0.26430803059933988),
    Ref(2, 4, 0, -3, 0.45943248645575357),
    Ref(3, 1, 0, -3, 0.10044208984043021),
    Ref(3, 2, 0, -3, 0.24168258333954593),
    Ref(3, 3, 0, -3, 0.43010327918878788),
    Ref(3, 4, 0, -3, 0.62086585901522295),
    Ref(4, 1, 0, -3, 0.23347257078724265),
    Ref(4, 2, 0, -3, 0.41246682731574696),
    Ref(4, 3, 0, -3, 0.59920554638858237),
    Ref(4, 4, 0, -3, 0.75706842643030481),
    Ref(5, 1, 0, -3, 0.40622520609149176),
    Ref(5, 2, 0, -3, 0.58706401160409682),
    Ref(5, 3, 0, -3, 0.74280993592520073),
    Ref(5, 4, 0, -3, 0.85689013122559088),
    Ref(6, 1, 0, -3, 0.58372633015331388),
    Ref(6, 2, 0, -3, 0.73621288949458141),
    Ref(6, 3, 0, -3, 0.84887398892051036),
    Ref(6, 4, 0, -3, 0.92197826671940764),
    Ref(7, 1, 0, -3, 0.73467269643127131),
    Ref(7, 2, 0, -3, 0.84538015762905905),
    Ref(7, 3, 0, -3, 0.91786370749082358),
    Ref(7, 4, 0, -3, 0.96014916380791715),
    Ref(8, 1, 0, -3, 0.8452308111972775),
    Ref(8, 2, 0, -3, 0.91637145605802672),
    Ref(8, 3, 0, -3, 0.95829641210917949),
    Ref(8, 4, 0, -3, 0.98085067316240493),
    Ref(1, 2, 0, -6, 1.1980077375303474e-07),
    Ref(1, 2, 0, -4.5, 0.00019178664955760337),
    Ref(1, 2, 0, 0, 0.95755963299045144),
    Ref(1, 2, 0, 5, 0.99999999999638867),
    Ref(1, 2, 25, -6, 4.3182821425567276e-05),
    Ref(1, 2, 25, -4.5, 0.0016109005925058406),
    Ref(1, 2, 25, 0, 0.95000972635418568),
    Ref(1, 2, 25, 5, 0.99999998531354961),
    Ref(1, 2, 250, -6, 3.0568955541124626e-07),
    Ref(1, 2, 250, -4.5, 0.00025483835515272002),
    Ref(1, 2, 250, 0, 0.95685656836410593),
    Ref(1, 2, 250, 5, 0.99999999999035105)
];

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    /* the embedded tables are the files */
    auto raw = cast(const(ubyte)[]) read("../gretl/data/urcdata.bin");
    bool same = raw.length == urcdata.length * double.sizeof;
    if (same) {
        auto fromFile = cast(const(double)[]) raw;
        foreach (i; 0 .. urcdata.length) {
            if (*cast(const(ulong) *) &fromFile[i] != *cast(const(ulong) *) &urcdata[i]) {
                same = false;
                break;
            }
        }
    }
    check(same && urcdata.length == 32045, "urcdata is bit for bit the content of urcdata.bin");
    bool matricesOk = embedded_matrices.length == 7;
    enum header = "gretl_binary_matrix".length;
    foreach (ref em; embedded_matrices) {
        auto bytes = cast(const(ubyte)[]) read("../gretl/data/" ~ em.name);
        const rows = *cast(const(int) *) &bytes[header], cols = *cast(const(int) *) &bytes[header + 4];
        auto vals = cast(const(double)[]) bytes[header + 8 .. $];
        bool ok = rows == em.rows && cols == em.cols && vals.length == em.val.length;
        if (ok) {
            foreach (i; 0 .. vals.length) {
                if (*cast(const(ulong) *) &vals[i] != *cast(const(ulong) *) &em.val[i]) ok = false;
            }
        }
        matricesOk = matricesOk && ok;
    }
    check(matricesOk, "the seven embedded DF-GLS matrices are bit for bit the content of their .bin files");

    /* agreement with the original C code */
    bool refOk = true;
    foreach (r; refs) {
        const p = mackinnon_pvalue(r.tau, r.T, r.niv, r.itv);

        if (!close(p, r.pval)) {
            printf("  mismatch: niv %d itv %d T %d tau %g: %.17g, C gives %.17g\n", r.niv, r.itv, r.T, r.tau, p, r.pval);
            refOk = false;
        }
    }
    check(refOk, "mackinnon_pvalue agrees with the original C code on the reference values");

    /* published asymptotic critical values (constant, one integrated variable) */
    check(fabs(mackinnon_pvalue(-3.4336, 0, 1, 2) - 0.01) < 5e-4 && fabs(mackinnon_pvalue(-2.8621, 0, 1, 2) - 0.05) < 5e-4
          && fabs(mackinnon_pvalue(-2.5671, 0, 1, 2) - 0.10) < 5e-4, "asymptotic 1, 5 and 10 percent critical values of the constant case");
    check(fabs(mackinnon_pvalue(-2.8913, 100, 1, 2) - 0.05) < 5e-4, "finite-sample (T = 100) 5 percent critical value");
    check(fabs(mackinnon_pvalue(-1.9416, 0, 1, 1) - 0.05) < 5e-4 && fabs(mackinnon_pvalue(-3.4104, 0, 1, 3) - 0.05) < 1e-3, "no-constant and trend cases at 5 percent");
    check(mackinnon_pvalue(-10, 0, 1, 2) < 1e-6 && mackinnon_pvalue(10, 0, 1, 2) > 0.999, "extreme statistics give p near 0 and 1");
    check(mackinnon_pvalue(-3, 100, 1, 2) < mackinnon_pvalue(-2, 100, 1, 2), "the p-value increases with the statistic");
    check(isNaN(mackinnon_pvalue(-2, 0, 0, 2)) && isNaN(mackinnon_pvalue(-2, 0, 9, 2)) && isNaN(mackinnon_pvalue(-2, 0, 1, 0)) && isNaN(mackinnon_pvalue(-2, 0, 1, 5)), "invalid niv or itv give NA");

    /* DF-GLS */
    int err = 0;
    check(fabs(dfgls_pvalue(-1.94, 0, 0, 0, 0, &err) - 0.05) < 2e-3 && err == 0, "DF-GLS asymptotic 5 percent critical value, constant (-1.94)");
    check(fabs(dfgls_pvalue(-2.85, 0, 1, 0, 0, &err) - 0.05) < 6e-3 && err == 0, "DF-GLS asymptotic 5 percent critical value, trend (about -2.89)");
    check(dfgls_pvalue(-1.94, 100, 0, 0, 0, &err) > dfgls_pvalue(-1.94, 0, 0, 0, 0, &err) && dfgls_pvalue(-1.94, 40, 0, 0, 0, &err) > dfgls_pvalue(-1.94, 100, 0, 0, 0, &err),
          "the finite-sample p-value is larger and falls towards the asymptotic one as T grows");
    foreach (trend; 0 .. 2) {
        foreach (pq; 0 .. 2) {
            const p = dfgls_pvalue(-2.2, 100, trend, 4, pq, &err);

            check(err == 0 && p > 0 && p < 1, "DF-GLS with test-down (Sephton surfaces)");
        }
    }
    check(dfgls_pvalue(-50, 100, 0, 0, 0, &err) == 1.0e-5 && dfgls_pvalue(5, 100, 0, 0, 0, &err) <= 1.0, "DF-GLS extremes: the conservative 1e-5 and a value not above 1");

    /* the library finds the functions, and adf_kpss uses them */
    check(get_plugin_function("mackinnon_pvalue") !is null && get_plugin_function("dfgls_pvalue") !is null, "get_plugin_function finds both functions");
    const pu = get_urc_pvalue(-2.8621, 0, 1, 2);
    check(!isNaN(pu) && fabs(pu - 0.05) < 5e-4, "get_urc_pvalue (used by adf_test) returns the real p-value");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
