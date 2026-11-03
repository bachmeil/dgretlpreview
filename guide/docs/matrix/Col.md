# Col

## Signature

```d
MatrixLine Col(M)(M m, size_t j);
```

## Parameters

- `m`: a `GretlMatrix` or a `Submatrix`.
- `j`: the (0-based) column number.

## Return value

A `MatrixLine`: a reference to column `j` of `m` (not a copy), which behaves like a `double[]`. It has `length`, can be indexed and iterated (`foreach (i, val; Col(m, 2))`, with `ref` to change the matrix), and `dup` returns a `double[]` copy. Assigning to it writes into the matrix, from a `double[]`, another line of the same length, or a single number for every element. A column number outside the matrix throws a `GretlException`. `Row` is the same for a row and `diag` for the diagonal.

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
    auto m = cbind([1.0, 2.0], [3.0, 4.0], [5.0, 6.0]);
    writeln(Col(m, 2).dup);                   // [5, 6]

    Col(m, 0) = Col(m, 2);                    // copy a column onto another
    show(m);

    Col(m, 1) = 0;                            // set a whole column to a number
    show(m);
}
```

Output:

```
[5, 6]
   5.0000   3.0000   5.0000
   6.0000   4.0000   6.0000
   5.0000   0.0000   5.0000
   6.0000   0.0000   6.0000
```
