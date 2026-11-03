# rank

## Signature

```d
size_t rank(M)(M m, double tol = 0.0);
```

## Parameters

- `m`: a `GretlMatrix` or a `Submatrix` (a view of part of a matrix). It is not changed.
- `tol`: the relative tolerance. Singular values not larger than `tol * max(nrow, ncol) * (largest singular value)` count as zero. The default 0 means machine epsilon (2.2e-16), as in NumPy and Matlab.

## Return value

A `size_t`, the numerical rank of `m`: the number of singular values above the tolerance. It is at most `min(nrow, ncol)`, and 0 for an empty matrix. The singular values are computed by `svd`, whatever the shape of the matrix. It can also be called as a method: `m.rank`.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    auto m = cbind([1.0, 2.0, 3.0], [2.0, 4.0, 6.0], [1.0, 0.0, 1.0]);   // the second column is twice the first
    writeln(rank(m));                   // 2
    writeln(rank(eye(4)));              // 4
}
```

Output:

```
2
4
```
