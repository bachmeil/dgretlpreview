# GretlMatrix

## Signature

```d
struct GretlMatrix {
    size_t nrow();
    size_t ncol();
    double opIndex(size_t i, size_t j);                       // m[i, j]
    double opIndexAssign(double x, size_t i, size_t j);       // m[i, j] = x
    GretlMatrix opIndex(MatrixRange rows, MatrixRange cols);  // m[a .. b, c .. d]: a copy
    GretlMatrix dup();
    GretlMatrixRaw* ptr();
    this(GretlMatrixRaw* raw);
}
```

## Parameters

- `nrow`, `ncol`: the number of rows and columns. An empty (`GretlMatrix.init`) matrix has both equal to 0.
- `m[i, j]`: the element in row `i` and column `j`, counting from 0. An index outside the matrix throws a `GretlException`. Assigning to `m[i, j]` changes the element.
- `m[a .. b, c .. d]`: a copy (a new `GretlMatrix`) of rows `a` to `b - 1` and columns `c` to `d - 1`. Use `m.sub`, `m.rows` or `m.cols` for a reference instead of a copy; `$` is available inside the brackets (`m[0 .. $, 1 .. 2]`).
- `dup`: a copy of the whole matrix.
- `ptr`: the raw pointer (`GretlMatrixRaw*`), to call the low-level functions of the port with; null for an empty matrix.
- the constructor `this(GretlMatrixRaw* raw)` takes over a raw matrix made with the low-level functions; it is freed with `gretl_matrix_free` when the last copy goes away.

## Return value

A `GretlMatrix` is a matrix of doubles, stored column by column and indexed from 0. It is reference counted: copying a `GretlMatrix` (assignment, passing to a function) shares the same data, which is freed when the last copy goes away, so there is nothing to free by hand. Because copies share data, changing an element through one copy is seen through the others; use `dup` for an independent copy.

A matrix is made by `cbind`, `rbind`, `eye`, by any function that returns one, by the importers and estimators, or from a raw pointer. Functions that take a matrix accept either a `GretlMatrix` or a `Submatrix`. See `operators` for arithmetic, `sub`, `rows` and `cols` for views, and `byRow`, `byColumn`, `Row`, `Col`, `diag` and `el` for iteration.

If a program has a type of its own called `GretlMatrix`, use a selective import or the qualified name `dgretl.matrix.GretlMatrix`.

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
    writefln("%d x %d", m.nrow, m.ncol);
    writeln(m[2, 1]);                         // row 2, column 1: 6

    m[0, 0] = 10;
    auto other = m;                           // shares the data
    other[0, 0] = 20;
    writeln(m[0, 0]);                         // 20: the same matrix

    auto copy = m.dup;                        // an independent copy
    copy[0, 0] = 30;
    writeln(m[0, 0]);

    auto block = m[1 .. 3, 0 .. 2];           // a copy of rows 1 and 2
    show(block);
    GretlMatrix none;
    writefln("%d x %d", none.nrow, none.ncol);
}
```

Output:

```
3 x 2
6
20
20
   2.0000   5.0000
   3.0000   6.0000
0 x 0
```
