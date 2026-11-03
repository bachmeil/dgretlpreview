# mahalanobis

## Signature

```d
GretlMatrix mahalanobis(MX)(MX x);
GretlMatrix mahalanobis(Dataset d, string[] names);
```

## Parameters

- `x`: an n x k `GretlMatrix` or `Submatrix` with `k >= 1` and more rows than columns; the rows are the observations. In the `Dataset` form, `d` and the `names` of the series to use.

## Return value

An n x 1 `GretlMatrix` with the Mahalanobis distance of each row of `x` from the centroid (the mean) of the rows, `sqrt((x - mean)' S^-1 (x - mean))` with `S` the sample covariance matrix. A row with a missing value gets NaN and is not used for the centroid or for `S`. Too few rows throws a `GretlException`. Large distances flag outliers.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    setSeed(8);
    auto x = randomMatrix!rnorm(100, 2);
    x[0, 0] = 6;                                         // make the first observation an outlier
    auto d = mahalanobis(x);
    writefln("%d distances; the first is %.3f, the second %.3f", d.nrow, d[0, 0], d[1, 0]);
}
```

Output:

```
100 distances; the first is 5.049, the second 1.846
```
