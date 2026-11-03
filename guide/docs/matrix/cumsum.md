# cumsum

## Signature

```d
GretlMatrix cumsum(M)(M m, Along dir = Along.column);
```

## Parameters

- `m`: a `GretlMatrix` or a `Submatrix` (a view of part of a matrix). It is not changed.
- `dir`: the direction. `Along.column` (the default) accumulates down each column, so element `[i, j]` is the sum of `m[0, j] .. m[i, j]`; `Along.row` accumulates across each row.

## Return value

A new `GretlMatrix` of the same size as `m` holding the cumulative sums. For a column vector the default gives its running total. A NaN makes the rest of its column (or row) NaN. An empty matrix gives an empty matrix. `cumprod` accumulates products. It can also be called as a method: `m.cumsum`.

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
    show(cumsum(m));                    // down each column
    show(cumsum(m, Along.row));         // across each row
}
```

Output:

```
   1.0000   4.0000
   3.0000   9.0000
   6.0000  15.0000
   1.0000   5.0000
   2.0000   7.0000
   3.0000   9.0000
```
