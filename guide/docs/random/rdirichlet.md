# rdirichlet

## Signature

```d
GretlMatrix rdirichlet(const(double)[] alpha, size_t n);
```

## Parameters

- `alpha`: the parameter vector of the Dirichlet distribution, at least two values, each finite and positive.
- `n`: the number of draws, at least 1.

## Return value

A `GretlMatrix` with `n` rows and `alpha.length` columns. Each row is one draw from the Dirichlet distribution: its elements are positive and add to 1. The expected value of column `j` is `alpha[j]` divided by the sum of `alpha`.

The draws come from the library's random number generator, so they change from call to call; call `setSeed(seed)` first to make them reproducible. A `GretlException` is thrown for invalid parameters (the message names the function and says what is needed).

The function can be called as a stream method (`streams[3].rdirichlet(alpha, 100)`).

## Example

```d
import std.stdio;
import dgretl;

void main() {
    setSeed(5);
    auto m = rdirichlet([1, 2, 3], 50_000);      // each row adds to 1

    writefln("%d x %d", m.nrow, m.ncol);
    writefln("row 0 sums to %.3f", m[0, 0] + m[0, 1] + m[0, 2]);
    foreach (j; 0 .. 3) {
        double s = 0;
        foreach (i; 0 .. m.nrow) {
            s += m[i, j];
        }
        writefln("column %d mean %.3f", j, s / m.nrow);      // near 1/6, 2/6, 3/6
    }
}
```

Output:

```
50000 x 3
row 0 sums to 1.000
column 0 mean 0.166
column 1 mean 0.334
column 2 mean 0.499
```
