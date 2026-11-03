# byColumn

## Signature

```d
ByColumn byColumn(M)(M m);
```

## Parameters

- `m`: a `GretlMatrix` or a `Submatrix`.

## Return value

A random-access range over the columns of `m`. Each element is a `Submatrix`, an nrow x 1 reference to one column, so changing it changes `m`. The range has `length`, can be indexed, and keeps the data of `m` alive. `byRow` does the same for rows.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    auto m = cbind([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]);
    foreach (j; 0 .. m.byColumn.length) {
        auto col = m.byColumn[j];             // indexed like an array
        double s = 0;
        foreach (i; 0 .. col.nrow) {
            s += col[i, 0];
        }
        writefln("column %d sums to %.1f", j, s);
    }
    writeln(m.byColumn.length);
}
```

Output:

```
column 0 sums to 6.0
column 1 sums to 15.0
2
```
