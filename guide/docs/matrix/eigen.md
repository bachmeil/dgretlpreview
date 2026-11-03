# eigen

## Signature

```d
Eigen eigen(M)(M m, bool vectors = true);

struct Eigen {
    GretlMatrix values;
    GretlMatrix vectors;
}
```

## Parameters

- `m`: a real symmetric matrix: a `GretlMatrix` or a `Submatrix` (a view of part of a matrix). It is not changed.
- `vectors`: if `false`, compute only the eigenvalues, which is cheaper. The default is `true`.

## Return value

An `Eigen` struct with two `GretlMatrix` fields. `values` is an n x 1 matrix of the eigenvalues in ascending order. `vectors` is the n x n matrix whose columns are the matching orthonormal eigenvectors, so that `m * v = v * diag(values)` with `v = vectors` and `t(v) * v` the identity (it is empty if `vectors` was `false`). The sign of each eigenvector is whatever LAPACK returns, so compare them only up to sign.

The matrix must be square and symmetric (to 1e-12 relative to its largest element), otherwise a `GretlException` is thrown. For a non-symmetric matrix use `eigenGeneral`. It can also be called as a method: `m.eigen`.

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
    auto m = cbind([2.0, 1.0], [1.0, 2.0]);
    auto e = eigen(m);
    show(e.values);                     // 1 and 3, ascending
    show(e.vectors);

    auto v = e.vectors;
    show(t(v) * v);                     // orthonormal: the identity
    show(m * v);                        // equals v times diag(values)
}
```

Output:

```
   1.0000
   3.0000
  -0.7071   0.7071
   0.7071   0.7071
   1.0000   0.0000
   0.0000   1.0000
  -0.7071   2.1213
   0.7071   2.1213
```
