# covariance

## Signature

```d
double covariance(S, T)(S x, T y);
GretlMatrix covariance(M)(M m);
GretlMatrix covariance(Dataset d, string[] names);
```

## Parameters

- `x`, `y`: two series of the same length (arrays or a row or column of a matrix). The observations where either value is missing are left out. Different lengths throw a `GretlException`.
- `m`: a `GretlMatrix` or `Submatrix` whose columns are the series.
- `d`, `names`: a `Dataset` and the names of the series to use, in order.

## Return value

For two series, a `double`: the sample covariance (divisor `n - 1`) over the observations where both are valid; NaN if there is none.

For a matrix or a `Dataset` with names, a k x k symmetric `GretlMatrix` of the covariances of the columns (variances on the diagonal), each entry computed from the rows where both of its columns are valid (pairwise deletion). `correlation` gives the correlations.

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
    writefln("%.4f", covariance(x, y));

    auto m = cbind(x, y);
    show(covariance(m));                          // variances on the diagonal

    auto gap = [2.0, 4.0, double.nan, 4.0, 5.0];
    writefln("%.4f", covariance(x, gap));          // the third observation is left out
}
```

Output:

```
1.5000
   2.5000   1.5000
   1.5000   1.5000
2.0000
```
