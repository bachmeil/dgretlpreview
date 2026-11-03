# sort

## Signature

```d
GretlMatrix sort(M)(M m, size_t first, size_t[] rest...);
```

## Parameters

- `m`: a `GretlMatrix` or a `Submatrix` (a view of part of a matrix). It is not changed.
- `first`: the (0-based) column that decides the order.
- `rest`: further columns, used to break ties in the order given.

## Return value

A new `GretlMatrix` with the rows of `m` in ascending order of column `first`; rows that tie on it are ordered by the next column given, and so on. The sort is stable (rows that tie on every column given keep their order) and a NaN sorts after every number. The matrix `m` is not changed. A column number outside the matrix throws a `GretlException`. `order` returns the permutation instead of the sorted matrix. It can also be called as a method: `m.sort(1, 0)`.

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
    auto m = cbind([3.0, 1.0, 2.0, 1.0], [10.0, 20.0, 30.0, 5.0]);
    show(m);
    show(sort(m, 0));                   // by the first column; the two rows with 1 keep their order
    show(m.sort(0, 1));                 // by the first column, then the second
}
```

Output:

```
   3.0000  10.0000
   1.0000  20.0000
   2.0000  30.0000
   1.0000   5.0000
   1.0000  20.0000
   1.0000   5.0000
   2.0000  30.0000
   3.0000  10.0000
   1.0000   5.0000
   1.0000  20.0000
   2.0000  30.0000
   3.0000  10.0000
```
