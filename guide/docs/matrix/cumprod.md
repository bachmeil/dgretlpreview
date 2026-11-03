# cumprod

## Signature

```d
GretlMatrix cumprod(M)(M m, Along dir = Along.column);
```

## Parameters

- `m`: a `GretlMatrix` or a `Submatrix` (a view of part of a matrix). It is not changed.
- `dir`: the direction. `Along.column` (the default) multiplies down each column; `Along.row` across each row.

## Return value

A new `GretlMatrix` of the same size as `m` holding the cumulative products: with the default direction, element `[i, j]` is the product of `m[0, j] .. m[i, j]`. A NaN makes the rest of its column (or row) NaN. An empty matrix gives an empty matrix. `cumsum` accumulates sums. It can also be called as a method: `m.cumprod`.

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
    auto m = cbind([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]);
    show(cumprod(m));
    show(cumprod(m, Along.row));
}
```

Output:

```
   1.0000   4.0000
   2.0000  20.0000
   6.0000 120.0000
   1.0000   4.0000
   2.0000  10.0000
   3.0000  18.0000
```
