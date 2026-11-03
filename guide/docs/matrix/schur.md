# schur

## Signature

```d
Schur schur(M)(M m, bool vectors = true);

struct Schur {
    GretlMatrix t;
    GretlMatrix z;
    GretlMatrix values;
}
```

## Parameters

- `m`: a real square matrix: a `GretlMatrix` or a `Submatrix` (a view of part of a matrix). It is not changed.
- `vectors`: if `false`, skip the Schur vectors `z`. The default is `true`.

## Return value

A `Schur` struct, the real Schur decomposition `m = z * t * t(z)`. `t` is n x n and quasi upper triangular: upper triangular except for 2 x 2 blocks on the diagonal, one for each pair of complex eigenvalues. `z` is n x n and orthogonal (the Schur vectors; empty if `vectors` is `false`). `values` is n x 2: the eigenvalues read off `t`, real parts in column 0 and imaginary parts in column 1, in the order of the diagonal of `t`, not sorted. A non-square matrix throws a `GretlException`, as does a failure of the iteration to converge. It can also be called as a method: `m.schur`.

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
    auto m = cbind([1.0, 3.0], [2.0, 4.0]);
    auto s = schur(m);
    show(s.t);                          // upper triangular: the eigenvalues are on its diagonal
    show(s.values);
    show(s.z * s.t * t(s.z));           // equals m
}
```

Output:

```
  -0.3723  -1.0000
   0.0000   5.3723
  -0.3723   0.0000
   5.3723   0.0000
   1.0000   2.0000
   3.0000   4.0000
```
