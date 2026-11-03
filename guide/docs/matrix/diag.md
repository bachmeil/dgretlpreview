# diag

## Signature

```d
MatrixLine diag(M)(M m);
```

## Parameters

- `m`: a `GretlMatrix` or a `Submatrix` (a view of part of a matrix).

## Return value

A `MatrixLine`, a reference to the diagonal of `m` (not a copy) with `min(nrow, ncol)` elements. It is a random-access range of `double`, so it can be read, passed to `foreach` (with `ref` to change the matrix), indexed, and its `dup` is a `double[]` copy. Assigning to it writes into the matrix: `diag(m) = [1.1, 3.4]`, `diag(m) = Row(m, 1)` (a line of the same length) and `diag(m) = 0` (a number, for every element) all change `m` itself. A `Row` or `Col` is the same type.

The diagonal of a rectangular matrix starts at `[0, 0]` and has as many elements as the shorter side. It can also be called as a method: `m.diag`.

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
    auto m = cbind([1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]);
    writeln(diag(m).dup);               // [1, 5, 9]

    diag(m) = 0;                        // zero the diagonal of m itself
    show(m);

    diag(m) = [10.0, 20.0, 30.0];
    foreach (ref v; diag(m)) {
        v += 1;
    }
    show(m);
}
```

Output:

```
[1, 5, 9]
   0.0000   4.0000   7.0000
   2.0000   0.0000   8.0000
   3.0000   6.0000   0.0000
  11.0000   4.0000   7.0000
   2.0000  21.0000   8.0000
   3.0000   6.0000  31.0000
```
