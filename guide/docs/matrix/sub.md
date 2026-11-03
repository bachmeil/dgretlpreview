# sub

## Signature

```d
SubIndexer sub(M)(M m);
```

## Parameters

- `m`: a `GretlMatrix` or a `Submatrix`.

## Return value

An object to index with two ranges, giving a `Submatrix`, a reference to the block of `m` in rows `a .. b - 1` and columns `c .. d - 1`: `m.sub[a .. b, c .. d]`. Unlike `m[a .. b, c .. d]`, which copies, the `Submatrix` reads and writes the elements of `m` itself, and keeps `m`'s data alive. A range outside the matrix throws a `GretlException`. A `Submatrix` has `sub`, `rows` and `cols` of its own and can be passed to any function that takes a matrix. `.dup` makes a `GretlMatrix` copy.

See also `rows` and `cols`, for whole rows or columns.

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
    auto block = m.sub[1 .. 3, 0 .. 2];       // rows 1 and 2, columns 0 and 1
    show(block);

    block[0, 0] = 100;                        // writes into m
    show(m);

    auto copy = m[1 .. 3, 0 .. 2];            // plain slicing makes a copy
    copy[0, 0] = -1;
    writeln(m[1, 0]);                         // still 100
}
```

Output:

```
   2.0000   5.0000
   3.0000   6.0000
   1.0000   4.0000   7.0000
 100.0000   5.0000   8.0000
   3.0000   6.0000   9.0000
100
```
