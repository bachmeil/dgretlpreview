# expm

## Signature

```d
GretlMatrix expm(M)(M m);
```

## Parameters

- `m`: a square matrix: a `GretlMatrix` or a `Submatrix` (a view of part of a matrix). It is not changed.

## Return value

A new `GretlMatrix`, the matrix exponential `exp(m) = I + m + m^2/2! + ...`, computed by scaling and squaring with a Pade approximation. It is called `expm` so that it is not confused with `std.math.exp`. For a diagonal matrix it is the exponential of each diagonal element, and for any `m`, `det(expm(m))` equals `exp(trace(m))`. A non-square matrix throws a `GretlException`. It can also be called as a method: `m.expm`.

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
    auto m = cbind([1.0, 0.0], [0.0, 2.0]);            // diagonal
    show(expm(m));                      // exp(1) and exp(2) on the diagonal

    auto a = cbind([0.0, 1.0], [-1.0, 0.0]);            // exp of a rotation generator is a rotation
    show(expm(a));
}
```

Output:

```
   2.7183   0.0000
   0.0000   7.3891
   0.5403  -0.8415
   0.8415   0.5403
```
