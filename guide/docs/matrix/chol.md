# chol

## Signature

```d
GretlMatrix chol(M)(M m, Triangle tri = Triangle.upper);
```

## Parameters

- `m`: a symmetric positive definite matrix: a `GretlMatrix` or a `Submatrix` (a view of part of a matrix). It is not changed.
- `tri`: which factor to return. `Triangle.upper` (the default) returns the upper triangle `U` with `t(U) * U = m`, as R and Matlab do; `Triangle.lower` returns the lower triangle `L` with `L * t(L) = m`, as gretl's own `cholesky()` does.

## Return value

A new `GretlMatrix`, the Cholesky factor of `m`: a triangular matrix of the same size. The two choices are transposes of each other. Only the lower triangle of `m` is read, and symmetry is not checked. A matrix that is not square or not positive definite throws a `GretlException`. It can also be called as a method: `m.chol`.

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
    auto m = cbind([4.0, 2.0], [2.0, 5.0]);
    auto U = chol(m);
    show(U);
    show(t(U) * U);                     // back to m

    auto L = chol(m, Triangle.lower);
    show(L);
}
```

Output:

```
   2.0000   1.0000
   0.0000   2.0000
   4.0000   2.0000
   2.0000   5.0000
   2.0000   0.0000
   1.0000   2.0000
```
