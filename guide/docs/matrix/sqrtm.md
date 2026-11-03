# sqrtm

## Signature

```d
GretlMatrix sqrtm(M)(M m);
```

## Parameters

- `m`: a square matrix: a `GretlMatrix` or a `Submatrix` (a view of part of a matrix). It is not changed.

## Return value

A new `GretlMatrix`, the principal matrix square root: the matrix `s` with `s * s = m` whose eigenvalues all have a non-negative real part. It is called `sqrtm` so that it is not confused with the square root of each element (that is `m.el ^^ 0.5`).

- A symmetric matrix (to 1e-12 relative) must be positive semidefinite, and the root is the symmetric positive semidefinite one, found from `eigen`. A singular symmetric matrix is fine.
- Any other square matrix uses the Denman-Beavers iteration, which needs `m` nonsingular and converges when `m` has no eigenvalue on the closed negative real axis. Otherwise there is no real principal root and a `GretlException` is thrown.

A non-square matrix throws a `GretlException`; a 0 x 0 matrix gives a 0 x 0 matrix. It can also be called as a method: `m.sqrtm`.

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
    auto m = cbind([4.0, 0.0], [0.0, 9.0]);
    show(sqrtm(m));                     // 2 and 3 on the diagonal

    auto p = cbind([2.0, 1.0], [1.0, 2.0]);
    auto r = sqrtm(p);
    show(r * r);                        // back to p
}
```

Output:

```
   2.0000   0.0000
   0.0000   3.0000
   2.0000   1.0000
   1.0000   2.0000
```
