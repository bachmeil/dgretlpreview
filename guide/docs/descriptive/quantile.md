# quantile

## Signature

```d
double quantile(S)(S x, double p, QuantileType type = QuantileType.q6);
GretlMatrix quantile(M)(M m, double p, QuantileType type = QuantileType.q6);

enum QuantileType { q6, q7, q8 }
```

## Parameters

- `x`: a series: an array of numbers (`double[]`, `int[]`, ...) or one row or column of a matrix (`Row(m, i)`, `Col(m, j)`). Missing values (NaN) are skipped.
- `m`: a `GretlMatrix` or `Submatrix`. The statistic is computed for each column.
- `p`: the probability, from 0 to 1. A value outside that range throws a `GretlException`.
- `type`: the rule for interpolating between the order statistics (Hyndman and Fan's types): `QuantileType.q6` (the default, the library's own) uses position `(n + 1) p`; `q7` (the default of R and NumPy) uses `(n - 1) p + 1`; `q8` (approximately median-unbiased) uses `(n + 1/3) p + 1/3`.

## Return value

For a series, a `double`: the `p` quantile of the valid values, interpolated linearly between the two neighbouring order statistics at the position the `type` gives. `p = 0` gives the minimum and `p = 1` the maximum. The result is NaN when the sample is too small for `p` (for `q6`, when `(n + 1) p` is below 1 or above `n`), or when there are no valid values. Missing values are skipped.

For a matrix, a 1 x k `GretlMatrix` with the quantile of each column. `quantiles` computes several probabilities at once and `fiveNumbers` the usual five.

## Example

```d
import std.stdio;
import dgretl;

void show(M)(M m) {
    foreach (i; 0 .. m.nrow) {
        foreach (j; 0 .. m.ncol) {
            writef("%9.4f", m[i, j]);
        }
        writeln();
    }
}

void main() {
    auto x = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
    writeln(quantile(x, 0.5));                      // the median
    writeln(quantile(x, 0.8));
    writeln(quantile(x, 0.8, QuantileType.q7));     // as R's quantile(x, 0.8)
    writeln(quantile(x, 0.8, QuantileType.q8));
    writeln(quantile(x, 0.9));                      // 8 values are too few for this p with q6: NaN
    writeln(quantile(x, 0.9, QuantileType.q7));     // but q7 can interpolate it

    auto m = cbind(x, [8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0]);
    show(quantile(m, 0.25));
}
```

Output:

```
4.5
7.4
6.2
7
nan
7.6
   4.0000   2.2500
```
