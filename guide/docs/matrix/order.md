# order

## Signature

```d
size_t[] order(M)(M m, size_t first, size_t[] rest...);
```

## Parameters

- `m`: a `GretlMatrix` or a `Submatrix` (a view of part of a matrix). It is not changed.
- `first`: the (0-based) column that decides the order.
- `rest`: further columns, used to break ties in the order given.

## Return value

A `size_t[]` with one entry for each row of `m`: the 0-based row numbers in the order that sorts the rows by column `first` (then by the columns in `rest`), like R's `order`. So `m.order(0)[0]` is the row with the smallest first column. Ties keep their original order and a NaN sorts last. `sort` returns the sorted matrix itself. A column number outside the matrix throws a `GretlException`. It can also be called as a method: `m.order(0)`.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    auto m = cbind([3.0, 1.0, 2.0, 1.0], [10.0, 20.0, 30.0, 5.0]);
    writeln(order(m, 0));               // [1, 3, 2, 0]: row 1 first (the first of the two 1s), then row 3
    writeln(m.order(0, 1));             // [3, 1, 2, 0]: the tie is broken by the second column
}
```

Output:

```
[1, 3, 2, 0]
[3, 1, 2, 0]
```
