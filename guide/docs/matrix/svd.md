# svd

## Signature

```d
SVD svd(M)(M m, bool vectors = true, bool full = false);

struct SVD {
    GretlMatrix u;
    GretlMatrix s;
    GretlMatrix v;
}
```

## Parameters

- `m`: a `GretlMatrix` or a `Submatrix` (a view of part of a matrix). It is not changed.
- `vectors`: if `false`, compute only the singular values, which is cheaper. The default is `true`.
- `full`: if `true`, return the full `u` (nrow x nrow) and `v` (ncol x ncol) when `m` is not square, instead of the first `min(nrow, ncol)` columns. The default is `false`.

## Return value

An `SVD` struct such that `m = u * diag(s) * t(v)`, with k = `min(nrow, ncol)`. `s` is a k x 1 matrix of the singular values, largest first. `u` is nrow x k (nrow x nrow if `full`), with orthonormal columns: the left singular vectors. `v` is ncol x k (ncol x ncol if `full`), with orthonormal columns: the right singular vectors (this is V, not its transpose). `u` and `v` are empty if `vectors` is `false`, and all three are empty for an empty matrix. The signs of the singular vectors are LAPACK's.

For a matrix with more rows than columns the library's own routine would square the condition number and lose accuracy in the small singular values; this function decomposes the transpose instead and swaps `u` and `v`, which gives the same decomposition computed accurately. It can also be called as a method: `m.svd`.

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
    auto m = cbind([3.0, 0.0, 0.0], [0.0, 2.0, 0.0]);  // 3 x 2
    auto d = svd(m);
    show(d.s);                          // 3 and 2

    auto S = eye(2);
    diag(S) = Col(d.s, 0);              // diag(s) as a matrix
    show(d.u * S * t(d.v));             // equals m
}
```

Output:

```
   3.0000
   2.0000
   3.0000   0.0000
   0.0000   2.0000
   0.0000   0.0000
```
