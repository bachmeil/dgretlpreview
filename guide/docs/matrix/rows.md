# rows

## Signature

```d
RowsIndexer rows(M)(M m);
```

## Parameters

- `m`: a `GretlMatrix` or a `Submatrix`.

## Return value

An object to index with one range, giving a `Submatrix` that refers to rows `a .. b - 1` of `m`, all columns: `m.rows[a .. b]`. It is a reference, not a copy: writing to it changes `m`. A range outside the matrix throws a `GretlException`. See `cols` for columns and `sub` for blocks.

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
    auto top = m.rows[0 .. 2];
    show(top);

    top[1, 1] = 0;                            // changes m
    show(m);
}
```

Output:

```
   1.0000   4.0000
   2.0000   5.0000
   1.0000   4.0000
   2.0000   0.0000
   3.0000   6.0000
```
