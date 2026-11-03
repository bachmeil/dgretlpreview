/* Tests for minpack/lbfgsb.d: the L-BFGS-B routines (setulb_, the f2c translation of version 3.0 of Nocedal's code), and the
   LBFGS option of gretl_bfgs.d that uses them (BFGS_max with OPT_L).

   setulb_ is driven as the original does: START, then FG (compute f and g) and NEW_X until a CONVERGENCE or error message. The
   problems are Rosenbrock's function in two variables, free and with the bound x <= 0.5, the 25-variable problem of the
   original driver1.f (bounds on every variable), a separable quadratic with bounds, and the 25-variable problem again with a tight
   tolerance and a memory of 8. The reference results are those of the original C code (minpack/lbfgsb.c, compiled with a
   small C driver and linked to the same BLAS): the final f, the final x, the number of function evaluations and iterations, and
   the first six characters of the final message agree exactly, bit for bit.
   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;
import core.stdc.string : memset, strcpy, strncmp;

import gretl.gretl_bfgs : BFGS_max;
import gretl.gretl_bfgsh : BFGS_CRIT_FUNC;
import gretl.gretl_matrix;
import gretl.libgretlh : Gretlopt;
import gretl.minpack.lbfgsb : setulb_;
import gretl.nlsh : CriterionType;

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

alias FgFn = void function(int n, const(double) * x, double * f, double * g);

void rosen2(int n, const(double) * x, double * f, double * g)
{
    const a = 1 - x[0], b = x[1] - x[0] * x[0];

    *f = a * a + 100 * b * b;
    g[0] = -2 * a - 400 * x[0] * b;
    g[1] = 200 * b;
}

void rosen25(int n, const(double) * x, double * f, double * g)
{
    double t1, t2;

    *f = .25 * (x[0] - 1) * (x[0] - 1);
    foreach (i; 1 .. n) *f += (x[i] - x[i-1] * x[i-1]) * (x[i] - x[i-1] * x[i-1]);
    *f *= 4;
    t1 = x[1] - x[0] * x[0];
    g[0] = 2 * (x[0] - 1) - 1.6e1 * x[0] * t1;
    foreach (i; 1 .. n - 1) {
        t2 = t1;
        t1 = x[i+1] - x[i] * x[i];
        g[i] = 8 * t2 - 1.6e1 * x[i] * t1;
    }
    g[n-1] = 8 * t1;
}

void quad5(int n, const(double) * x, double * f, double * g)
{
    *f = 0;
    foreach (i; 0 .. n) {
        *f += (x[i] - (i + 1)) * (x[i] - (i + 1));
        g[i] = 2 * (x[i] - (i + 1));
    }
}

struct Result {
    char[6] task;
    int nfg, iter;
    double f;
    double[] x;
}

Result run(FgFn fn, int n, int m, double[] x, double[] l, double[] u, int[] nbd, double factr, double pgtol)
{
    char[60] task = ' ', csave = 0;
    int[4] lsave;
    int[44] isave;
    double[29] dsave;
    double f = 0;
    auto g = new double[n];
    auto wa = new double[(2 * m + 5) * n + 11 * m * m + 8 * m];
    auto iwa = new int[3 * n];
    int nfg = 0;
    Result r;

    strcpy(task.ptr, "START");
    while (true) {
        setulb_(&n, &m, x.ptr, l.ptr, u.ptr, nbd.ptr, &f, g.ptr, &factr, &pgtol, wa.ptr, iwa.ptr, task.ptr, csave.ptr, lsave.ptr, isave.ptr, dsave.ptr);
        if (!strncmp(task.ptr, "FG", 2)) {
            fn(n, x.ptr, &f, g.ptr);
            nfg++;
        } else if (!strncmp(task.ptr, "NEW_X", 5)) {
            continue;
        } else {
            break;
        }
    }
    r.task[] = task[0 .. 6];
    r.nfg = nfg;
    r.iter = isave[29];
    r.f = f;
    r.x = x.dup;
    return r;
}

struct Ref {
    string name, task;
    int nfg, iter;
    double f;
    double[] x;
}

/* from the C reference program (the final message's first six characters, evaluations, iterations, f, x) */
immutable Ref[] refs = [
    Ref("rosen2", "CONVER", 47, 38, 0x1.fe9b12b85f7f8p-43, [0x1.fffff67ee8801p-1, 0x1.ffffee46a0819p-1]),
    Ref("rosen2_bounded", "CONVER", 10, 6, 0x1p-2, [0x1p-1, 0x1p-2]),
    Ref("rosen25", "CONVER", 28, 23, 0x1.29d3d6251e542p-30, [0x1.0000270591022p+0, 0x1.000033149fd7cp+0, 0x1.00004bbd2b3d6p+0, 0x1.0000552ea900fp+0, 0x1.000061a1bafc9p+0, 0x1.000075ddb7f6dp+0, 0x1.00008bfc64fa4p+0, 0x1.0000ca7855743p+0, 0x1.00014a59e0b7cp+0, 0x1.0002747bc7d09p+0, 0x1.00050c071aabbp+0, 0x1.000a613341de2p+0, 0x1.00152dbbd5de9p+0, 0x1.002a8d2739ce5p+0, 0x1.00553ca946553p+0, 0x1.00aa85620362fp+0, 0x1.015586fdd9c4dp+0, 0x1.02acd6e43fb71p+0, 0x1.0560e1f9e1a1bp+0, 0x1.0adea89939691p+0, 0x1.1633827655bf5p+0, 0x1.2e53f6f48ffa2p+0, 0x1.650a8671498d6p+0, 0x1.f1f62e4ce2a79p+0, 0x1.e44f0ed5981b2p+1]),
    Ref("quad5", "CONVER", 3, 2, 0x1.18p+3, [0x1p+0, 0x1p+1, 0x1.4p+1, 0x1.4p+1, 0x1.4p+1]),
    Ref("rosen25_m8", "CONVER", 39, 33, 0x1.0c67510fc5c94p-47, [0x1.00000136504e9p+0, 0x1.00000297cd9b5p+0, 0x1.000005578038ap+0, 0x1.00000ab01601p+0, 0x1.0000156f547b7p+0, 0x1.00002ac444aeep+0, 0x1.00005580d8ff4p+0, 0x1.0000aae9eb663p+0, 0x1.000155d476739p+0, 0x1.0002ab89f22d8p+0, 0x1.00055721fde13p+0, 0x1.000aae4e91e54p+0, 0x1.00155d0478c52p+0, 0x1.002abbc807b6p+0, 0x1.00557ebcf14bbp+0, 0x1.00ab19f5910a9p+0, 0x1.0156a64a57e0ep+0, 0x1.02af1735ede61p+0, 0x1.05656288127b4p+0, 0x1.0ae7e304a57eep+0, 0x1.1646b5d9b671p+0, 0x1.2e7da68870d6bp+0, 0x1.656cb2722c888p+0, 0x1.f30857ff5c1ddp+0, 0x1.e664c3aef92dap+1]),
];

bool sameBits(double a, double b)
{
    return *cast(ulong *) &a == *cast(ulong *) &b;
}

version (Windows) version = Tolerant;
version (LDC) version = Tolerant;

bool matches(Result r, immutable Ref ref_)
{
    version (Tolerant) {
        /* the reference values come from the C code built by gcc (x87/SSE, glibc): the C library of Windows and LDC's code
           generation (fused multiply-adds) round some results differently, and the last-digit differences can change the
           number of iterations: compare the solution to a tight tolerance */
        bool near(double a, double b) { return fabs(a - b) <= 1e-9 + 1e-6 * fabs(b); }
        bool okw = r.task[] == ref_.task && near(r.f, ref_.f) && r.x.length == ref_.x.length;

        if (okw) {
            foreach (i; 0 .. r.x.length) okw = okw && near(r.x[i], ref_.x[i]);
        }
        return okw;
    } else {
        bool ok = r.task[] == ref_.task && r.nfg == ref_.nfg && r.iter == ref_.iter && sameBits(r.f, ref_.f) && r.x.length == ref_.x.length;

        if (ok) {
            foreach (i; 0 .. r.x.length) ok = ok && sameBits(r.x[i], ref_.x[i]);
        }
        return ok;
    }
}

/* the log-likelihood of a normal regression y = a + b x with unit variance is minus half the sum of squares: BFGS_max maximizes it */
struct RegData {
    double[] xs, ys;
}

extern (C) double regLoglik(const(double) * b, void * data)
{
    auto d = cast(RegData *) data;
    double ss = 0;

    foreach (i; 0 .. d.xs.length) {
        const e = d.ys[i] - b[0] - b[1] * d.xs[i];

        ss += e * e;
    }
    return -0.5 * ss;
}

extern (C) int regScore(double * b, double * g, int n, BFGS_CRIT_FUNC f, void * data)
{
    auto d = cast(RegData *) data;

    g[0] = g[1] = 0;
    foreach (i; 0 .. d.xs.length) {
        const e = d.ys[i] - b[0] - b[1] * d.xs[i];

        g[0] += e;
        g[1] += e * d.xs[i];
    }
    return 0;
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    /* Rosenbrock, free */
    auto r1 = run(&rosen2, 2, 5, [-1.2, 1.0], [0.0, 0.0], [0.0, 0.0], [0, 0], 1e7, 1e-5);
    check(matches(r1, refs[0]), "Rosenbrock (2 variables, no bounds): result identical to the original C code");
    check(fabs(r1.x[0] - 1) < 1e-5 && fabs(r1.x[1] - 1) < 1e-5, "Rosenbrock: the minimum is at (1, 1)");

    /* Rosenbrock with x <= 0.5 */
    auto r2 = run(&rosen2, 2, 5, [0.0, 0.0], [0.0, 0.0], [0.5, 0.0], [3, 0], 1e7, 1e-5);
    check(matches(r2, refs[1]), "Rosenbrock with an upper bound: result identical to the original C code");
    check(r2.x[0] == 0.5 && r2.x[1] == 0.25 && r2.f == 0.25, "Rosenbrock with x <= 0.5: the constrained minimum is (0.5, 0.25), f = 0.25");

    /* the 25-variable problem of driver1.f */
    double[] l25, u25, x25;
    int[] nbd25;

    foreach (i; 0 .. 25) {
        nbd25 ~= 2;
        l25 ~= (i % 2 == 0) ? 1.0 : -100.0;
        u25 ~= 100.0;
        x25 ~= 3.0;
    }
    auto r3 = run(&rosen25, 25, 5, x25.dup, l25, u25, nbd25, 1e7, 1e-5);
    check(matches(r3, refs[2]), "the 25-variable bounded problem of driver1.f: result identical to the original C code");

    /* separable quadratic with bounds */
    auto r4 = run(&quad5, 5, 3, [0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0], [2.5, 2.5, 2.5, 2.5, 2.5], [2, 2, 2, 2, 2], 1e7, 1e-5);
    check(matches(r4, refs[3]), "separable quadratic with bounds: result identical to the original C code");
    check(r4.x == [1.0, 2.0, 2.5, 2.5, 2.5] && r4.f == 8.75, "the solution is (1, 2, 2.5, 2.5, 2.5) with f = 8.75");

    /* tight tolerance, memory 8 */
    auto r5 = run(&rosen25, 25, 8, x25.dup, l25, u25, nbd25, 10, 1e-10);
    check(matches(r5, refs[4]), "the 25-variable problem with factr = 10, pgtol = 1e-10 and m = 8: identical to the original C code");

    /* an invalid problem is reported through the task string */
    {
        char[60] task = ' ', csave = 0;
        int n = 2, m = 5, nbd = 0;
        int[4] lsave;
        int[44] isave;
        double[29] dsave;
        double[2] x = [0.0, 0.0], l = 0, u = 0, g = 0;
        double f = 0, factr = -1, pgtol = 0;
        auto wa = new double[(2 * m + 5) * n + 11 * m * m + 8 * m];
        auto iwa = new int[3 * n];
        int[2] nbd2 = 0;

        strcpy(task.ptr, "START");
        setulb_(&n, &m, x.ptr, l.ptr, u.ptr, nbd2.ptr, &f, g.ptr, &factr, &pgtol, wa.ptr, iwa.ptr, task.ptr, csave.ptr, lsave.ptr, isave.ptr, dsave.ptr);
        check(!strncmp(task.ptr, "ERROR: FACTR .LT. 0", 19), "a negative factr is reported as an error");
    }

    /* the LBFGS option of BFGS_max: a normal regression with known OLS solution */
    {
        RegData d;
        double[2] b = [0.0, 0.0];
        int fncount, grcount;

        foreach (i; 0 .. 30) {
            d.xs ~= i * 0.1;
            d.ys ~= 1.5 + 2.0 * d.xs[i] + ((i * 7) % 5 - 2) * 0.05;
        }
        const err = BFGS_max(b.ptr, 2, 100, 1e-12, &fncount, &grcount, &regLoglik, CriterionType.loglik, &regScore, &d, null, Gretlopt.OPT_L, null);
        double sx = 0, sy = 0, sxx = 0, sxy = 0;
        const nn = d.xs.length;

        foreach (i; 0 .. nn) { sx += d.xs[i]; sy += d.ys[i]; sxx += d.xs[i] * d.xs[i]; sxy += d.xs[i] * d.ys[i]; }
        const slope = (nn * sxy - sx * sy) / (nn * sxx - sx * sx), icpt = (sy - slope * sx) / nn;
        check(err == 0 && fabs(b[0] - icpt) < 1e-6 && fabs(b[1] - slope) < 1e-6, "BFGS_max with OPT_L (L-BFGS-B) maximizes a regression likelihood to the OLS solution");
        check(fncount > 0 && grcount > 0, "function and gradient evaluations are counted");
    }

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
