# logdet

## Signature

```d
double logdet(M)(M m, bool absolute = false);
```

## Parameters

- `m`: a `GretlMatrix` or a `Submatrix` (a view of part of a matrix). It is not changed.
- `absolute`: if `true`, return the log of the absolute value of the determinant, which is finite for any nonsingular matrix. The default is `false`.

## Return value

A `double`, the natural logarithm of the determinant of the square matrix `m`, computed from one LU factorization as the sum of the logs of the pivots, so that it stays finite when the determinant itself would overflow or underflow. A determinant of zero gives `-double.infinity`, and a negative determinant gives NaN unless `absolute` is `true`. The log-determinant of a 0 x 0 matrix is 0. The matrix must be square; otherwise a `GretlException` is thrown. It can also be called as a method: `m.logdet`.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    auto m = cbind([4.0, 2.0], [7.0, 6.0]);
    writeln(logdet(m));                 // log(10)

    auto big = eye(500) * 10.0;         // det = 10^500 overflows a double
    writeln(logdet(big));               // but its log is fine: 500 * log(10)

    auto neg = cbind([0.0, 1.0], [1.0, 0.0]);           // det = -1
    writeln(logdet(neg));
    writeln(logdet(neg, true));
}
```

Output:

```
2.30259
1151.29
nan
0
```
