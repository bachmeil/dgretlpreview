# longRunVariance

## Signature

```d
double longRunVariance(S)(S x, int lag = -1, double mean = double.nan);
GretlMatrix longRunVariance(M)(M m, int lag = -1, double mean = double.nan);
```

## Parameters

- `x`: a series: an array of numbers (`double[]`, `int[]`, ...) or one row or column of a matrix (`Row(m, i)`, `Col(m, j)`). Unlike the other statistics, a missing value is not skipped: it makes the result NaN.
- `m`: a `GretlMatrix` or `Submatrix`. The statistic is computed for each column.
- `lag`: the number of lags in the Bartlett (triangular) weights. The default -1 takes the library's choice, `n^(1/3)`.
- `mean`: the mean to take deviations from. The default NaN uses the sample mean.

## Return value

For a series, a `double`: the long-run variance, `gamma_0 + 2 * sum_{j=1..lag} (1 - j / (lag + 1)) gamma_j`, with `gamma_j` the autocovariances about `mean` (this is the variance of the sample mean times `n`, allowing for serial correlation). For an uncorrelated series it is close to the ordinary variance. The result is NaN if the series has no values or contains a missing value.

For a matrix, a 1 x k `GretlMatrix` with the value for each column.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    setSeed(3);
    auto e = rnorm(500);                          // uncorrelated: the long-run variance is near 1
    writefln("%.3f", longRunVariance(e));
    writefln("%.3f", longRunVariance(e, 0));        // no lags: just the variance ...
    writefln("%.3f", variance(e, population: true));  // ... with divisor n

    auto ar = new double[500];                      // a positively correlated series
    ar[0] = e[0];
    foreach (t; 1 .. 500) {
        ar[t] = 0.7 * ar[t - 1] + e[t];
    }
    writefln("variance %.3f, long-run variance %.3f", variance(ar), longRunVariance(ar, 12));   // the long-run value is larger
}
```

Output:

```
0.961
1.107
1.107
variance 1.933, long-run variance 7.472
```
