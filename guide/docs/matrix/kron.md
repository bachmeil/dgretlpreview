# kron

## Signature

```d
GretlMatrix kron(MA, MB)(MA a, MB b);
```

## Parameters

- `a`: a `GretlMatrix` or `Submatrix`, the left factor.
- `b`: a `GretlMatrix` or `Submatrix`, the right factor.

## Return value

A new `GretlMatrix` with `a.nrow * b.nrow` rows and `a.ncol * b.ncol` columns, the Kronecker product: each element `a[i, j]` is replaced by the block `a[i, j] * b`. It can also be called as a method: `m.kron(b)`.

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
    auto a = cbind([1.0, 2.0], [3.0, 4.0]);
    auto b = eye(2);
    show(kron(a, b));                   // 4 x 4: each element of a times the identity
}
```

Output:

```
   1.0000   0.0000   3.0000   0.0000
   0.0000   1.0000   0.0000   3.0000
   2.0000   0.0000   4.0000   0.0000
   0.0000   2.0000   0.0000   4.0000
```
