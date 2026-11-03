# stddev

## Signature

```d
double stddev(S)(S x, bool population = false);
GretlMatrix stddev(M)(M m, bool population = false);
```

## Parameters

- `x`: a series: an array of numbers (`double[]`, `int[]`, ...) or one row or column of a matrix (`Row(m, i)`, `Col(m, j)`). Missing values (NaN) are skipped.
- `m`: a `GretlMatrix` or `Submatrix`. The statistic is computed for each column.
- `population`: if `true`, divide by `n` instead of `n - 1`. The default is `false`.

## Return value

For a series, a `double`: the standard deviation, the square root of `variance(x, population)`. Missing values (NaN) are skipped. The statistic of a series with no valid observation (or an empty series) is NaN.

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
    writeln(stddev(x));                        // sqrt(32 / 7)

    auto gap = [1.0, double.nan, 3.0, 5.0];
    writeln(stddev(gap));                  // the missing value is skipped

    auto m = cbind(x, [8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0]);
    show(stddev(m));                       // one value per column
    writeln(stddev(Col(m, 1)));            // one column of a matrix
}
```

Output:

```
2.13809
2
   2.1381   2.4495
2.44949
```
