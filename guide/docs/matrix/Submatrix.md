# Submatrix

## Signature

```d
struct Submatrix {
    size_t nrow();
    size_t ncol();
    size_t stride();
    bool isContiguous();
    double opIndex(size_t i, size_t j);
    double opIndexAssign(double x, size_t i, size_t j);
    GretlMatrix opIndex(MatrixRange rows, MatrixRange cols);  // a copy
    GretlMatrix dup();
}
```

## Parameters

- `nrow`, `ncol`: the size of the block.
- `m[i, j]`: the element in row `i` and column `j` of the block, counting from 0. Reading and assigning go through to the matrix the block was cut from.
- `m[a .. b, c .. d]`: a copy (a `GretlMatrix`) of part of the block.
- `dup`: a `GretlMatrix` copy of the whole block.
- `stride`: the distance in memory between the columns of the block (a block is not contiguous in the column-by-column storage), and `isContiguous()`: `true` if its elements are one run of memory.

## Return value

A `Submatrix` is a reference to a block of a matrix, made by `m.sub[a .. b, c .. d]`, `m.rows[a .. b]` or `m.cols[a .. b]` (and by the iterators `byRow` and `byColumn`). It reads and writes the elements of the matrix it was cut from, and keeps that matrix's data alive. It has `sub`, `rows` and `cols` of its own. Every function that takes a `GretlMatrix` also takes a `Submatrix`, and it works in arithmetic and as the parameter of the distribution functions. Slicing a `Submatrix` with `[a .. b, c .. d]` makes a copy, as for a `GretlMatrix`.

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
    auto s = m.sub[0 .. 2, 1 .. 3];           // rows 0, 1 and columns 1, 2
    writefln("%d x %d, stride %d, contiguous %s", s.nrow, s.ncol, s.stride, s.isContiguous);

    s[1, 1] = 0;                              // writes into m
    writeln(m[1, 2]);

    show(s * s);                              // works in arithmetic
    auto c = s.dup;                           // a GretlMatrix copy
    c[0, 0] = -1;
    writeln(m[0, 1]);                         // m is unchanged by changing the copy
}
```

Output:

```
2 x 2, stride 3, contiguous false
0
  51.0000  28.0000
  20.0000  35.0000
4
```
