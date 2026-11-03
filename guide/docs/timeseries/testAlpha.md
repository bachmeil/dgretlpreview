# testAlpha

## Signature

```d
TestResult testAlpha(ref VecmFit fit, AlphaRestriction restriction);

struct AlphaRestriction { GretlMatrix R, q; }
struct TestResult {
    double statistic;
    double pvalue;
    int df;
}
```

## Parameters

- `fit`: a `VecmFit`. It is left as it was.
- `restriction`: an `AlphaRestriction`: a homogeneous linear restriction `R alpha = 0`, with `R` a matrix with `neqns` columns (and `q` zero or left out), applied to every column of `alpha`. For example `R = [1 0]` tests that the first variable is weakly exogenous (does not adjust to the long-run relationship), and `R = [1 0; 0 1]` that neither of the first two does.

## Return value

A `TestResult` for the likelihood-ratio test of the restriction on the adjustment coefficients: `statistic`, `df` (`rank` times the number of rows of `R`) and the chi-square `pvalue`. A small p-value is evidence against the restriction. A `GretlException` is thrown if no restriction is given or it is not valid.

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

    auto first = AlphaRestriction(cbind([1.0], [0.0]));      // y1 is weakly exogenous: alpha[0] = 0
    auto t = testAlpha(fit, first);
    writefln("y1 does not adjust: LR = %.3f, df = %d, p = %.4f", t.statistic, t.df, t.pvalue);

    auto second = AlphaRestriction(cbind([0.0], [1.0]));     // y2 is weakly exogenous: alpha[1] = 0
    auto u = testAlpha(fit, second);
    writefln("y2 does not adjust: LR = %.3f, p = %.6f", u.statistic, u.pvalue);
}
```

Output:

```
y1 does not adjust: LR = 0.336, df = 1, p = 0.5619
y2 does not adjust: LR = 42.257, p = 0.000000
```
