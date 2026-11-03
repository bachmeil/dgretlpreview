# pinv

## Signature

```d
GretlMatrix pinv(M)(M m, double tol = 0.0);
```

## Parameters

- `m`: a `GretlMatrix` or a `Submatrix` (a view of part of a matrix). It is not changed.
- `tol`: the relative tolerance. Singular values not larger than `tol * max(nrow, ncol) * (largest singular value)` count as zero. The default 0 means machine epsilon (2.2e-16), as in NumPy and Matlab.

## Return value

A new `ncol` x `nrow` `GretlMatrix`, the Moore-Penrose pseudoinverse of `m`, computed from the singular value decomposition. For a square nonsingular matrix it is `inv(m)`, for a tall matrix of full column rank it is `inv(t(m) * m) * t(m)`, and for the zero matrix it is the zero matrix. Unlike `inv` it exists for singular and rectangular matrices. An empty matrix gives an empty `ncol` x `nrow` matrix. It can also be called as a method: `m.pinv`.

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
    auto m = cbind([1.0, 2.0, 3.0], [1.0, 1.0, 1.0]);  // 3 x 2, full column rank
    auto p = pinv(m);
    show(p);                                            // 2 x 3
    show(p * m);                                        // the 2 x 2 identity

    auto s = cbind([1.0, 2.0], [2.0, 4.0]);             // singular
    show(pinv(s));
}
```

Output:

```
  -0.5000  -0.0000   0.5000
   1.3333   0.3333  -0.6667
   1.0000  -0.0000
   0.0000   1.0000
   0.0400   0.0800
   0.0800   0.1600
```
