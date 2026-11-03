# sst

## Signature

```d
double sst(S)(S x);
GretlMatrix sst(M)(M m);
```

## Parameters

- `x`: a series: an array of numbers (`double[]`, `int[]`, ...) or one row or column of a matrix (`Row(m, i)`, `Col(m, j)`). Missing values (NaN) are skipped.
- `m`: a `GretlMatrix` or `Submatrix`. The statistic is computed for each column.

## Return value

For a series, a `double`: the sum of squared deviations from the mean, `sum((x - mean)^2)`. Missing values (NaN) are skipped. The statistic of a series with no valid observation (or an empty series) is NaN.

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
    writeln(sst(x));                            // 32

    auto gap = [1.0, double.nan, 3.0, 5.0];
    writeln(sst(gap));                  // the missing value is skipped

    auto m = cbind(x, [8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0]);
    show(sst(m));                       // one value per column
    writeln(sst(Col(m, 1)));            // one column of a matrix
}
```

Output:

```
32
8
  32.0000  42.0000
42
```
