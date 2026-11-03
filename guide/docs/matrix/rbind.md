# rbind

## Signature

```d
GretlMatrix rbind(Args...)(Args args);
```

## Parameters

- `args`: one or more arguments, each a `GretlMatrix` or `Submatrix` (all with the same number of columns), a number (repeated to fill a row of that width), or an array of numbers (one row, whose length must equal the number of columns). A 0 x 0 matrix is left out, so `rbind(acc, row)` can start from `GretlMatrix.init`.

## Return value

A new `GretlMatrix` holding the arguments one below the other, as in R: the rows of the first argument, then those of the second, and so on. Its number of columns is that of the matrix and array arguments; with only numbers it has one column. A `GretlException` is thrown if the column counts do not match. `cbind` puts the arguments side by side.

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
    auto a = rbind([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]);   // two arrays become two rows
    show(a);

    auto b = rbind(a, [7.0, 8.0, 9.0]);                  // add a row at the bottom
    show(b);

    GretlMatrix acc;                                     // start from an empty matrix
    foreach (i; 0 .. 3) {
        acc = rbind(acc, [i, i * i]);
    }
    show(acc);
}
```

Output:

```
   1.0000   2.0000   3.0000
   4.0000   5.0000   6.0000
   1.0000   2.0000   3.0000
   4.0000   5.0000   6.0000
   7.0000   8.0000   9.0000
   0.0000   0.0000
   1.0000   1.0000
   2.0000   4.0000
```
