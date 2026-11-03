# rowDistances

## Signature

```d
GretlMatrix rowDistances(MX, MY)(MX x, MY y, DistanceMetric metric = DistanceMetric.euclidean);

enum DistanceMetric { euclidean, manhattan, hamming, chebyshev, cosine, mahalanobis }
```

## Parameters

- `x`: an m x k `GretlMatrix` or `Submatrix`; its rows are the first set of points.
- `y`: a p x k matrix with the same number of columns; its rows are the second set of points.
- `metric`: the distance: `euclidean` (the default), `manhattan` (sum of absolute differences), `hamming` (the proportion of coordinates that differ), `chebyshev` (the largest absolute difference), `cosine` (1 minus the cosine of the angle between the two vectors) or `mahalanobis` (with the covariance matrix of the rows of `x`). Missing values are not allowed.

## Return value

An m x p `GretlMatrix` whose element `[i, j]` is the distance between row `i` of `x` and row `j` of `y`. A different number of columns, or a missing value, throws a `GretlException`. (This is the library's `distance` function; it is named `rowDistances` here.)

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
    auto x = cbind([0.0, 3.0], [0.0, 4.0]);              // two points: (0, 0) and (3, 4)
    auto y = cbind([0.0], [0.0]);                        // the origin

    show(rowDistances(x, y));                            // 0 and 5
    show(rowDistances(x, y, DistanceMetric.manhattan));  // 0 and 7
    show(rowDistances(x, y, DistanceMetric.chebyshev));  // 0 and 4
}
```

Output:

```
   0.0000
   5.0000
   0.0000
   7.0000
   0.0000
   4.0000
```
