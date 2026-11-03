# rinvwishart

## Signature

```d
GretlMatrix rinvwishart(M)(M scale, int df);
GretlMatrix rinvwishart(M)(M scale, int df, size_t n);
```

## Parameters

- `scale`: the p x p scale matrix, symmetric and positive definite, as a `GretlMatrix`, `Submatrix` or anything else accepted as a matrix. All elements must be finite.
- `df`: the degrees of freedom, an integer of at least p.
- `n` (second form): the number of draws, at least 1.

## Return value

A `GretlMatrix`. In the first form it is one draw from the inverse Wishart distribution, a symmetric positive definite p x p matrix. In the second form it has `n` rows and p(p + 1)/2 columns, and each row is the `vech` of one draw: the lower triangle taken column by column. For `df > p + 1` the mean of a draw is `scale / (df - p - 1)`.

The draws come from the library's random number generator, so they change from call to call; call `setSeed(seed)` first to make them reproducible. A `GretlException` is thrown if the matrix is not square, if `df` is below p, or if the matrix is not positive definite.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    auto S = cbind([2.0, 0.5], [0.5, 1.0]);     // a 2 x 2 scale matrix

    setSeed(5);
    auto one = rinvwishart(S, 10);               // one 2 x 2 draw
    writefln("%d x %d, symmetric: %s", one.nrow, one.ncol, one[0, 1] == one[1, 0]);

    auto many = rinvwishart(S, 10, 20_000);      // rows are (v11, v21, v22)
    writefln("%d x %d", many.nrow, many.ncol);
    double s = 0;
    foreach (i; 0 .. many.nrow) {
        s += many[i, 0];
    }
    writefln("mean of v11 %.3f", s / many.nrow);  // near 2 / (10 - 2 - 1) = 0.286
}
```

Output:

```
2 x 2, symmetric: true
20000 x 3
mean of v11 0.286
```
