# norm

## Signature

```d
double norm(M)(M m, Norm kind = Norm.frobenius);
```

## Parameters

- `m`: a `GretlMatrix` or a `Submatrix` (a view of part of a matrix). It is not changed.
- `kind`: which norm: `Norm.frobenius` (the default: the square root of the sum of squares of all elements), `Norm.one` (the largest column sum of absolute values), `Norm.infinity` (the largest row sum of absolute values), `Norm.spectral` (the 2-norm: the largest singular value, which needs an SVD) or `Norm.max` (the largest absolute value of an element).

## Return value

A `double`, the chosen norm of `m`. For a vector (an n x 1 or 1 x n matrix) `Norm.frobenius` and `Norm.spectral` are the usual Euclidean length. The Frobenius norm is scaled by the largest element so that it does not overflow or underflow. An empty matrix has norm 0 and a matrix containing a NaN has norm NaN. It can also be called as a method: `m.norm(kind)`.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    auto m = cbind([1.0, -3.0], [2.0, 4.0]);           // rows (1, 2) and (-3, 4)
    writefln("frobenius %.4f", norm(m));
    writefln("one       %.4f", norm(m, Norm.one));
    writefln("infinity  %.4f", norm(m, Norm.infinity));
    writefln("spectral  %.4f", norm(m, Norm.spectral));
    writefln("max       %.4f", norm(m, Norm.max));
}
```

Output:

```
frobenius 5.4772
one       6.0000
infinity  7.0000
spectral  5.1167
max       4.0000
```
