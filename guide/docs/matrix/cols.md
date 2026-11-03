# cols

## Signature

```d
ColsIndexer cols(M)(M m);
```

## Parameters

- `m`: a `GretlMatrix` or a `Submatrix`.

## Return value

An object to index with one range, giving a `Submatrix` that refers to columns `a .. b - 1` of `m`, all rows: `m.cols[a .. b]`. It is a reference, not a copy: writing to it changes `m`. A range outside the matrix throws a `GretlException`. A `Submatrix` of columns can be passed anywhere a matrix is accepted, for example as the parameter of a distribution function. See `rows` for rows and `sub` for blocks.

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
    auto last2 = m.cols[1 .. 3];
    show(last2);
    show(last2 * last2.t);            // a Submatrix works in arithmetic
}
```

Output:

```
   3.0000   5.0000
   4.0000   6.0000
  34.0000  42.0000
  42.0000  52.0000
```
