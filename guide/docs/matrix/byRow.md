# byRow

## Signature

```d
ByRow byRow(M)(M m);
```

## Parameters

- `m`: a `GretlMatrix` or a `Submatrix`.

## Return value

A random-access range over the rows of `m`, like `byKey` for an associative array. Each element is a `Submatrix`, a 1 x ncol reference to one row, so `foreach (row; m.byRow) row[0, 2] = 0;` changes `m`. The range has `length`, can be indexed, and keeps the data of `m` alive. `byColumn` does the same for columns.

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
    size_t i = 0;
    foreach (row; m.byRow) {
        writefln("row %d has %d elements, first two sum to %.1f", i++, row.ncol, row[0, 0] + row[0, 1]);
    }

    foreach (row; m.byRow) {
        row[0, 0] = 0;                        // a reference into m
    }
    show(m);
}
```

Output:

```
row 0 has 2 elements, first two sum to 5.0
row 1 has 2 elements, first two sum to 7.0
row 2 has 2 elements, first two sum to 9.0
   0.0000   4.0000
   0.0000   5.0000
   0.0000   6.0000
```
