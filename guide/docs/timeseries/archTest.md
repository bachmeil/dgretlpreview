# archTest

## Signature

```d
TestResult archTest(U)(U u, int order);
TestResult archTest(Dataset d, string name, int order);

struct TestResult {
    double statistic;
    double pvalue;
    int df;
}
```

## Parameters

- `u`: the series to test, usually the residuals of a regression (an array or a one-column matrix). No missing values.
- `order`: the number of lags of the squared series in the test regression, from 1 to `n - 3`.

## Return value

A `TestResult` for Engle's LM test of the null hypothesis of no ARCH effects: the squared series is regressed on a constant and its own `order` lags, `statistic` is `(n - order)` times the R squared of that regression, `df` is `order`, and `pvalue` is the upper tail of the chi-square distribution. A small p-value is evidence of autoregressive conditional heteroskedasticity. In the `Dataset` form, `name` is the name of the series.

## Example

```d
import std.stdio;
import std.math: sqrt;
import dgretl;

void main() {
    setSeed(22);
    auto z = rnorm(500);
    auto u = new double[500];                        // an ARCH(1) process: the variance depends on the last shock
    u[0] = z[0];
    foreach (t; 1 .. 500) {
        u[t] = sqrt(0.2 + 0.7 * u[t - 1] ^^ 2) * z[t];
    }

    auto t = archTest(u, 4);
    writefln("ARCH data:   LM = %.3f, df = %d, p = %.6f", t.statistic, t.df, t.pvalue);
    auto w = archTest(z, 4);
    writefln("white noise: LM = %.3f, p = %.4f", w.statistic, w.pvalue);
}
```

Output:

```
ARCH data:   LM = 219.999, df = 4, p = 0.000000
white noise: LM = 0.351, p = 0.9863
```
