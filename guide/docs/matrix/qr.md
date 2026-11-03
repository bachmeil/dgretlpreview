# qr

## Signature

```d
QR qr(M)(M m);

struct QR {
    GretlMatrix q;
    GretlMatrix r;
}
```

## Parameters

- `m`: a matrix with at least as many rows as columns: a `GretlMatrix` or a `Submatrix` (a view of part of a matrix). It is not changed.

## Return value

A `QR` struct, the thin QR decomposition `m = q * r` by Householder reflections (without pivoting). `q` has the shape of `m` and orthonormal columns; `r` is the `ncol` x `ncol` upper triangular matrix, whose diagonal elements can have either sign. A `GretlException` is thrown if `m` has more columns than rows. It can also be called as a method: `m.qr`.

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
    auto m = cbind([1.0, 1.0, 1.0], [1.0, 2.0, 3.0]);  // 3 x 2
    auto f = qr(m);
    show(f.r);                          // upper triangular
    show(t(f.q) * f.q);                 // orthonormal columns: the identity
    show(f.q * f.r);                    // equals m
}
```

Output:

```
  -1.7321  -3.4641
   0.0000  -1.4142
   1.0000  -0.0000
  -0.0000   1.0000
   1.0000   1.0000
   1.0000   2.0000
   1.0000   3.0000
```
