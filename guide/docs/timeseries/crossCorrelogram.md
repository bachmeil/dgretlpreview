# crossCorrelogram

## Signature

```d
CrossCorrelogram crossCorrelogram(X, Y)(X x, Y y, int maxLag);
CrossCorrelogram crossCorrelogram(Dataset d, string xname, string yname, int maxLag);

struct CrossCorrelogram {
    GretlMatrix lags;
    GretlMatrix ccf;
}
```

## Parameters

- `x`, `y`: two series of the same length (arrays, or one-column matrices). They must not be constant or contain missing values.
- `maxLag`: the largest lag, from 1 to `n / 2`.

## Return value

A `CrossCorrelogram` struct with two `(2 * maxLag + 1)` x 1 matrices: `lags`, the lags from `-maxLag` to `maxLag`, and `ccf`, the sample cross-correlation at each lag, `r(k) = sum (x(t) - xbar)(y(t - k) - ybar) / sqrt(sum (x - xbar)^2 sum (y - ybar)^2)`: the correlation of `x(t)` with `y(t - k)`. If `x` follows `y` with a delay of `j` observations, the peak is at `k = j`; if `y` follows `x`, it is at `k = -j`. In the `Dataset` form `xname` and `yname` name the series.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    setSeed(16);
    auto x = rnorm(300);
    auto y = new double[300];                        // y(t) = x(t - 2) plus a little noise
    auto e = rnorm(300);
    foreach (t; 0 .. 300) {
        y[t] = (t >= 2 ? x[t - 2] : 0.0) + 0.3 * e[t];
    }

    auto c = crossCorrelogram(x, y, 4);
    foreach (i; 0 .. 9) {
        writefln("lag %2.0f: %7.4f", c.lags[i, 0], c.ccf[i, 0]);   // the peak shows the delay
    }
}
```

Output:

```
lag -4: -0.0248
lag -3: -0.0520
lag -2:  0.9596
lag -1: -0.0529
lag  0: -0.0257
lag  1:  0.0016
lag  2:  0.0411
lag  3:  0.0150
lag  4:  0.0007
```
