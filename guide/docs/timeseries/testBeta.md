# testBeta

## Signature

```d
TestResult testBeta(ref VecmFit fit, BetaRestriction restriction);

struct BetaRestriction { GretlMatrix R, q; }
struct TestResult {
    double statistic;
    double pvalue;
    int df;
}
```

## Parameters

- `fit`: a `VecmFit`. It is left as it was.
- `restriction`: a `BetaRestriction`: a homogeneous linear restriction `R beta = 0`, with `R` a matrix with one column for each row of `beta` (and `q` zero or left out), applied to every cointegrating vector. For example `R = [1 1]` says that each cointegrating vector is proportional to `(1, -1)`.

## Return value

A `TestResult` for the likelihood-ratio test of the restriction on the cointegrating vectors: `statistic` is the LR statistic, `df` is `rank` times the number of rows of `R` and `pvalue` is the upper tail of the chi-square distribution. A small p-value is evidence against the restriction. A `GretlException` is thrown if no restriction is given or the restriction is not valid. (`vecm` with a `restriction` estimates the restricted model.)

## Example

```d
import std.stdio;
import dgretl;

double[] arSeries(size_t n, double phi) {     // an AR(1) series from the library's random numbers
    auto e = rnorm(n);
    auto y = new double[n];

    y[0] = e[0];
    foreach (t; 1 .. n) {
        y[t] = phi * y[t - 1] + e[t];
    }
    return y;
}

void main() {
    setSeed(41);
    auto n = 250;
    auto walk = arSeries(n, 1.0);
    auto e = arSeries(n, 0.3);
    auto y1 = walk.dup;
    auto y2 = new double[n];
    foreach (t; 0 .. n) {
        y2[t] = walk[t] + e[t];
    }
    auto fit = vecm(cbind(y1, y2), 2, 1);

    auto equal = BetaRestriction(cbind([1.0], [1.0]));       // beta proportional to (1, -1): true here
    auto t = testBeta(fit, equal);
    writefln("(1, -1): LR = %.3f, df = %d, p = %.4f", t.statistic, t.df, t.pvalue);

    auto wrong = BetaRestriction(cbind([1.0], [-2.0]));      // beta proportional to (2, 1): false
    auto u = testBeta(fit, wrong);
    writefln("(2, 1):  LR = %.3f, df = %d, p = %.6f", u.statistic, u.df, u.pvalue);
}
```

Output:

```
(1, -1): LR = 0.001, df = 1, p = 0.9697
(2, 1):  LR = 93.470, df = 1, p = 0.000000
```
