# cbind

## Signature

```d
GretlMatrix cbind(Args...)(Args args);
```

## Parameters

- `args`: one or more arguments, each a `GretlMatrix` or `Submatrix` (all with the same number of rows), a number (repeated to fill a column of that height), or an array of numbers (one column, whose length must equal the number of rows). A 0 x 0 matrix is left out, so `cbind(acc, x)` can start from `GretlMatrix.init`.

## Return value

A new `GretlMatrix` holding the arguments side by side, as in R: the columns of the first argument, then those of the second, and so on. Its number of rows is that of the matrix and array arguments; with only numbers it has one row. A `GretlException` is thrown if the row counts do not match.

The operator `m ~ x` is `cbind(m, x)`, and `m ~= x` adds the columns of `x` to `m` in place (see `operators`). `rbind` stacks rows in the same way.

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
    auto a = cbind([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]);   // two arrays become two columns
    show(a);

    auto b = cbind(1, a);                                // a column of ones in front
    show(b);

    auto c = cbind(a, a.cols[0 .. 1]);                   // a matrix and a submatrix
    writefln("%d x %d", c.nrow, c.ncol);
}
```

Output:

```
   1.0000   4.0000
   2.0000   5.0000
   3.0000   6.0000
   1.0000   1.0000   4.0000
   1.0000   2.0000   5.0000
   1.0000   3.0000   6.0000
3 x 3
```
