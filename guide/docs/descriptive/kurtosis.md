# kurtosis

## Signature

```d
double kurtosis(S)(S x);
GretlMatrix kurtosis(M)(M m);
```

## Parameters

- `x`: a series: an array of numbers (`double[]`, `int[]`, ...) or one row or column of a matrix (`Row(m, i)`, `Col(m, j)`). Missing values (NaN) are skipped.
- `m`: a `GretlMatrix` or `Submatrix`. The statistic is computed for each column.

## Return value

For a series, a `double`: the EXCESS kurtosis: the fourth standardized moment (divisors `n`) minus 3, so that it is 0 for a normal distribution. It is NaN for a constant series. Missing values (NaN) are skipped. The statistic of a series with no valid observation (or an empty series) is NaN.

For a matrix, a 1 x k `GretlMatrix` with the statistic of each of the k columns.

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
    writeln(kurtosis(x));

    auto gap = [1.0, double.nan, 3.0, 5.0];
    writeln(kurtosis(gap));                  // the missing value is skipped

    auto m = cbind(x, [8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0]);
    show(kurtosis(m));                       // one value per column
    writeln(kurtosis(Col(m, 1)));            // one column of a matrix
}
```

Output:

```
-0.21875
-1.5
  -0.2188  -1.2381
-1.2381
```
