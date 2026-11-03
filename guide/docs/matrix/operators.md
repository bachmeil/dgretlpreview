# operators

## Signature

```d
GretlMatrix opBinary(op)(...)   // m + mm, m - mm, m * mm, m / mm, m + a, a + m, m * a, a / m, ...
```

## Parameters

- `m`, `mm`: a `GretlMatrix` or a `Submatrix`.
- `a`: a number.
- `m + mm`, `m - mm`: element by element; the sizes must be equal.
- `m + a`, `a + m`, `m - a`, `a - m`, `m * a`, `a * m`, `m / a`: with a number, applied to every element. `a / m` divides the number by each element.
- `m * mm`: the matrix product (the number of columns of `m` must equal the number of rows of `mm`).
- `m / mm`: matrix right division as in Matlab, the `x` with `x * mm = m`: `m * inv(mm)` for a square `mm`, the least-squares solution (by QR) when `mm` has more columns than rows.
- `m.el * mm`, `m * mm.el`, `m.el / mm`, `m.el ^^ a`, `a ^^ m.el`, `m.el ^^ mm`: element by element, through `el`. A plain `m ^^ n` is refused (use `power`).
- `m ~ x`, `x ~ m`: `cbind(m, x)`, where `x` is a matrix, number or array.
- `m ~= x`: add the columns of `x` to the right of `m`, in place if no other handle refers to the data and otherwise by rebinding `m` to a new matrix.

## Return value

Every operator returns a new `GretlMatrix` (except `~=`, which changes `m`). Sizes that do not fit throw a `GretlException`. The operators are available for a `GretlMatrix` and for a `Submatrix`, in either argument.

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
    auto b = cbind([0.0, 1.0], [1.0, 0.0]);

    show(a + b);
    show(a * 2 + 1);
    show(a * b);                              // matrix product
    show(a.el * b);                           // element by element
    show(a / a);                              // right division: the identity
    show(1 / a.el);                           // the reciprocal of each element

    a ~= [5.0, 6.0];                          // add a column
    show(a);
    show(a ~ 7);                              // a column of 7s on the right
}
```

Output:

```
   1.0000   4.0000
   3.0000   4.0000
   3.0000   7.0000
   5.0000   9.0000
   3.0000   1.0000
   4.0000   2.0000
   0.0000   3.0000
   2.0000   0.0000
   1.0000  -0.0000
  -0.0000   1.0000
   1.0000   0.3333
   0.5000   0.2500
   1.0000   3.0000   5.0000
   2.0000   4.0000   6.0000
   1.0000   3.0000   5.0000   7.0000
   2.0000   4.0000   6.0000   7.0000
```
