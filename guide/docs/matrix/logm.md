# logm

## Signature

```d
GretlMatrix logm(M)(M m);
```

## Parameters

- `m`: a square matrix: a `GretlMatrix` or a `Submatrix` (a view of part of a matrix). It is not changed.

## Return value

A new `GretlMatrix`, the principal matrix logarithm: the matrix `l` with `expm(l) = m` whose eigenvalues have imaginary parts between -pi and pi. It is called `logm` so that it is not confused with the logarithm of each element.

- A symmetric matrix (to 1e-12 relative) must be positive definite, and the logarithm is found from `eigen`.
- Any other square matrix uses inverse scaling and squaring, and needs `m` nonsingular with no eigenvalue on the closed negative real axis, which is what a real principal logarithm requires. Otherwise a `GretlException` is thrown.

A non-square matrix throws a `GretlException`; a 0 x 0 matrix gives a 0 x 0 matrix. It can also be called as a method: `m.logm`.

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
    auto l = logm(m);
    show(l);
    show(expm(l));                      // back to m
}
```

Output:

```
   0.5493   0.5493
   0.5493   0.5493
   2.0000   1.0000
   1.0000   2.0000
```
