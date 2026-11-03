# correlation

## Signature

```d
double correlation(S, T)(S x, T y);
GretlMatrix correlation(M)(M m);
GretlMatrix correlation(Dataset d, string[] names);
```

## Parameters

- `x`, `y`: two series of the same length (arrays or a row or column of a matrix). The observations where either value is missing are left out. Different lengths throw a `GretlException`.
- `m`: a `GretlMatrix` or `Submatrix` whose columns are the series.
- `d`, `names`: a `Dataset` and the names of the series to use, in order.

## Return value

For two series, a `double`: the correlation coefficient over the observations where both are valid; NaN if there is none or one of the series is constant.

For a matrix or a `Dataset` with names, a k x k symmetric `GretlMatrix` of the correlations of the columns (ones on the diagonal), each entry computed from the rows where both of its columns are valid (pairwise deletion). `covariance` gives the covariances.

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
    auto x = [1.0, 2.0, 3.0, 4.0, 5.0];
    auto y = [2.0, 4.0, 5.0, 4.0, 5.0];
    writefln("%.4f", correlation(x, y));

    auto z = [5.0, 4.0, 3.0, 2.0, 1.0];
    show(correlation(cbind(x, y, z)));            // x and z are perfectly negatively correlated
}
```

Output:

```
0.7746
   1.0000   0.7746  -1.0000
   0.7746   1.0000  -0.7746
  -1.0000  -0.7746   1.0000
```
