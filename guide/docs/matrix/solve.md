# solve

## Signature

```d
GretlMatrix solve(MX, MY)(MX X, MY y);
```

## Parameters

- `X`: the n x p matrix of regressors (a `GretlMatrix` or `Submatrix`). It needs n >= p and full column rank.
- `y`: the n x m matrix of right-hand sides (a `GretlMatrix` or `Submatrix`); a single n x 1 column is the usual case.

## Return value

A new p x m `GretlMatrix` `b` that minimizes the sum of squares of `y - X * b`, computed by a QR decomposition: Matlab's `X \ y` for data in columns. When `X` is square and nonsingular, `b` solves `X * b = y` exactly; for a vector `y`, `b` is the vector of ordinary-least-squares coefficients. A `GretlException` is thrown if `X` and `y` have different numbers of rows, if `X` has more columns than rows, or if `X` is rank deficient. (`m / mm` is the right-division operator, which solves `x * mm = m`.)

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
    auto X = cbind(1, [1.0, 2.0, 3.0, 4.0]);            // a constant and a regressor
    auto y = cbind([2.1, 3.9, 6.2, 7.8]);

    auto b = solve(X, y);                               // least-squares coefficients
    show(b);

    auto A = cbind([2.0, 1.0], [1.0, 3.0]);             // a square system A b = c
    auto c = cbind([3.0, 5.0]);
    show(solve(A, c));
}
```

Output:

```
   0.1500
   1.9400
   0.8000
   1.4000
```
