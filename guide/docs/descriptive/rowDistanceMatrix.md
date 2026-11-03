# rowDistanceMatrix

## Signature

```d
GretlMatrix rowDistanceMatrix(MX)(MX x, DistanceMetric metric = DistanceMetric.euclidean);
```

## Parameters

- `x`: an m x k `GretlMatrix` or `Submatrix` with at least two rows; its rows are the points. Missing values are not allowed.
- `metric`: the distance, as for `rowDistances`. The default is `DistanceMetric.euclidean`.

## Return value

A symmetric m x m `GretlMatrix` whose element `[i, j]` is the distance between rows `i` and `j` of `x`, with zeros on the diagonal. Fewer than two rows or a missing value throws a `GretlException`.

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
    auto x = cbind([0.0, 3.0, 0.0], [0.0, 4.0, 1.0]);   // the points (0, 0), (3, 4) and (0, 1)
    show(rowDistanceMatrix(x));
    show(rowDistanceMatrix(x, DistanceMetric.manhattan));
}
```

Output:

```
   0.0000   5.0000   1.0000
   5.0000   0.0000   4.2426
   1.0000   4.2426   0.0000
   0.0000   7.0000   1.0000
   7.0000   0.0000   6.0000
   1.0000   6.0000   0.0000
```
