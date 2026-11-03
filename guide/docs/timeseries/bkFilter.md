# bkFilter

## Signature

```d
GretlMatrix bkFilter(X)(X x, int lower = 6, int upper = 32, int k = 12);
GretlMatrix bkFilter(Dataset d, string name, int lower = 6, int upper = 32, int k = 12);
```

## Parameters

- `x`: the series: an array of numbers, or a `GretlMatrix` or `Submatrix` with one column. Missing values (NaN) are not allowed and throw a `GretlException`. It must be longer than `2k + 1`.
- `lower`: the shortest period (in observations) of the cycles to keep, at least 2. The default is 6, the usual choice for quarterly data.
- `upper`: the longest period of the cycles to keep, greater than `lower`. The default is 32.
- `k`: the number of leads and lags in the moving average, at least 1. The default is 12.

## Return value

An n x 1 `GretlMatrix` with the Baxter-King band-pass filtered series: the components of `x` with periods between `lower` and `upper` observations. The first `k` and the last `k` values are NaN, because the moving average needs `k` observations on each side. In the `Dataset` form, `name` is the name of a series of the `Dataset` instead of the data themselves (the series must have no missing values).

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
    setSeed(12);
    auto x = arSeries(100, 0.9);
    auto bk = bkFilter(x);                           // periods of 6 to 32 observations, k = 12

    writefln("first value %s, value 12: %.4f, last value %s", bk[0, 0], bk[12, 0], bk[99, 0]);

    size_t missing = 0;
    foreach (t; 0 .. 100) {
        if (bk[t, 0] != bk[t, 0]) missing++;
    }
    writeln(missing);                                // 2 * k
}
```

Output:

```
first value nan, value 12: 0.8659, last value nan
24
```
