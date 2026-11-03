# det

## Signature

```d
double det(M)(M m);
```

## Parameters

- `m`: a `GretlMatrix` or a `Submatrix` (a view of part of a matrix). It is not changed.

## Return value

A `double`, the determinant of the square matrix `m` (a closed form for 2 x 2 and 3 x 3 matrices, an LU factorization otherwise). The determinant of a 0 x 0 matrix is 1. The matrix must be square; otherwise a `GretlException` is thrown. For large matrices the determinant can overflow or underflow; `logdet` avoids that. It can also be called as a method: `m.det`.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    auto m = cbind([4.0, 2.0], [7.0, 6.0]);
    writeln(det(m));                    // 4 * 6 - 7 * 2 = 10
    writeln(det(eye(3)));
}
```

Output:

```
10
1
```
