# lu

## Signature

```d
LU lu(M)(M m);

struct LU {
    GretlMatrix l;
    GretlMatrix u;
    GretlMatrix p;
    size_t[] perm;
    bool singular;
}
```

## Parameters

- `m`: a `GretlMatrix` or a `Submatrix` (a view of part of a matrix). It is not changed.

## Return value

An `LU` struct, the LU decomposition with partial pivoting, `p * m = l * u`, for a matrix of any shape (k = `min(nrow, ncol)`):

- `l`: nrow x k, unit lower triangular (ones on the diagonal, zeros above).
- `u`: k x ncol, upper triangular (upper trapezoidal if `m` is not square).
- `p`: nrow x nrow, the permutation matrix of the row interchanges.
- `perm`: the same permutation as 0-based indices: row `i` of `p * m` is row `perm[i]` of `m`.
- `singular`: `true` if an exact zero turned up on the diagonal of `u`. The factors are still valid, but `u` cannot be used to solve a system.

The input is copied, not changed. It can also be called as a method: `m.lu`.

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
    auto m = cbind([1.0, 3.0], [2.0, 4.0]);            // rows (1, 2) and (3, 4)
    auto f = lu(m);
    show(f.l);
    show(f.u);
    writeln(f.perm);                    // the rows were swapped to get the larger pivot
    show(f.p * m);                      // equals l * u
    show(f.l * f.u);
    writeln(f.singular);
}
```

Output:

```
   1.0000   0.0000
   0.3333   1.0000
   3.0000   4.0000
   0.0000   0.6667
[1, 0]
   3.0000   4.0000
   1.0000   2.0000
   3.0000   4.0000
   1.0000   2.0000
false
```
