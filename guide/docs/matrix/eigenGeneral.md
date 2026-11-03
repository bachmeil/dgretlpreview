# eigenGeneral

## Signature

```d
EigenGeneral eigenGeneral(M)(M m, bool vectors = true);

struct EigenGeneral {
    GretlMatrix values;
    GretlMatrix vectorsRe;
    GretlMatrix vectorsIm;
    bool isReal();
}
```

## Parameters

- `m`: a real square matrix: a `GretlMatrix` or a `Submatrix` (a view of part of a matrix). It is not changed.
- `vectors`: if `false`, compute only the eigenvalues. The default is `true`.

## Return value

An `EigenGeneral` struct. The eigenvalues of a real non-symmetric matrix can be complex, and a `GretlMatrix` holds only real numbers, so each complex quantity comes as two real matrices. `values` is n x 2: the real parts in column 0 and the imaginary parts in column 1 (complex ones come in conjugate pairs, the one with the positive imaginary part first). `vectorsRe` and `vectorsIm` are n x n: the real and imaginary parts of the right eigenvectors, column `j` belonging to eigenvalue `j`, so that `m v = lambda v`. Each eigenvector has unit length. With `vectors = false` the two vector matrices are empty. `isReal()` is `true` if every eigenvalue is real.

The eigenvalues are not sorted (LAPACK's order). For a symmetric matrix `eigen` gives real, sorted values and orthogonal vectors. A non-square matrix throws a `GretlException`. It can also be called as a method: `m.eigenGeneral`.

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
    auto rot = cbind([0.0, 1.0], [-1.0, 0.0]);        // a rotation by 90 degrees: eigenvalues are +i and -i
    auto e = eigenGeneral(rot);
    show(e.values);                     // real parts, imaginary parts
    writeln(e.isReal());

    auto tri = cbind([2.0, 0.0], [1.0, 3.0]);           // triangular: eigenvalues 2 and 3
    auto f = eigenGeneral(tri, false);
    show(f.values);
    writeln(f.isReal());
}
```

Output:

```
   0.0000   1.0000
   0.0000  -1.0000
false
   2.0000   0.0000
   3.0000   0.0000
true
```
