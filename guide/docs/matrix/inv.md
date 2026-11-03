# inv

## Signature

```d
GretlMatrix inv(M)(M m);
```

## Parameters

- `m`: a `GretlMatrix` or a `Submatrix` (a view of part of a matrix). It is not changed.

## Return value

A new `GretlMatrix`, the inverse of the square matrix `m`, so that `m * inv(m)` is the identity to rounding error. A symmetric positive definite matrix is inverted by a Cholesky factorization and any other by an LU factorization. The matrix must be square; otherwise a `GretlException` is thrown. A singular matrix throws a `GretlException`. To solve a system, `solve` is faster and more accurate than forming the inverse; for a matrix that may be singular or rectangular see `pinv`. It can also be called as a method: `m.inv`.

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
    auto m = cbind([4.0, 2.0], [7.0, 6.0]);
    auto mi = inv(m);
    show(mi);
    show(m * mi);                       // the identity
}
```

Output:

```
   0.6000  -0.7000
  -0.2000   0.4000
   1.0000   0.0000
   0.0000   1.0000
```
