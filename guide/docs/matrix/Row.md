# Row

## Signature

```d
MatrixLine Row(M)(M m, size_t i);
```

## Parameters

- `m`: a `GretlMatrix` or a `Submatrix`.
- `i`: the (0-based) row number.

## Return value

A `MatrixLine`: a reference to row `i` of `m` (not a copy), which behaves like a `double[]`. It has `length`, can be indexed and iterated (`foreach (val; Row(m, 1))`, with `ref` to change the matrix), and `dup` returns a `double[]` copy. Assigning to it writes into the matrix, from a `double[]`, another line of the same length, or a single number for every element. A row number outside the matrix throws a `GretlException`. `Col` is the same for a column and `diag` for the diagonal.

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
    writeln(Row(m, 1).dup);                   // [2, 4, 6]

    foreach (ref v; Row(m, 0)) {
        v *= 10;                              // changes m
    }
    Row(m, 1) = [0.0, 0.0, 0.0];
    show(m);
}
```

Output:

```
[2, 4, 6]
  10.0000  30.0000  50.0000
   0.0000   0.0000   0.0000
```
