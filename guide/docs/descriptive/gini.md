# gini

## Signature

```d
double gini(S)(S x);
GretlMatrix gini(M)(M m);
```

## Parameters

- `x`: a series: an array of numbers (`double[]`, `int[]`, ...) or one row or column of a matrix (`Row(m, i)`, `Col(m, j)`). Missing values (NaN) are skipped.
- `m`: a `GretlMatrix` or `Submatrix`. The statistic is computed for each column.

## Return value

For a series, a `double`: the Gini coefficient of inequality, from 0 (all values equal) towards 1 (one value has everything). It is meant for non-negative data such as incomes. Missing values (NaN) are skipped. The statistic of a series with no valid observation (or an empty series) is NaN.

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
    writeln(gini(x));

    auto gap = [1.0, double.nan, 3.0, 5.0];
    writeln(gini(gap));                  // the missing value is skipped

    auto m = cbind(x, [8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0]);
    show(gini(m));                       // one value per column
    writeln(gini(Col(m, 1)));            // one column of a matrix
}
```

Output:

```
0.2125
0.296296
   0.2125   0.2917
0.291667
```
