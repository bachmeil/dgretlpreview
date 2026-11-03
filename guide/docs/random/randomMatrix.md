# randomMatrix

## Signature

```d
GretlMatrix randomMatrix(alias f, Args...)(size_t nrow, size_t ncol, Args args);
```

## Parameters

- `f` (a template argument): one of the distribution functions that returns a `double[]`, such as `rnorm`, `runif`, `rgamma`, `rbinom`, written `randomMatrix!rnorm`.
- `nrow`: the number of rows of the result.
- `ncol`: the number of columns of the result.
- `args`: the arguments of `f` that come after its `n`, in the same order and with the same defaults.

## Return value

A `GretlMatrix` with `nrow` rows and `ncol` columns, filled column by column with `nrow * ncol` draws of `f`, so the draws are the same as those of `f(nrow * ncol, args)` from the same generator state. Invalid parameters throw a `GretlException` as for `f` itself.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    setSeed(1);
    auto m = randomMatrix!rnorm(1000, 3, 10.0, 2.0);   // 1000 x 3, N(10, 2^2)
    writefln("%d x %d", m.nrow, m.ncol);

    double s = 0;
    foreach (i; 0 .. m.nrow) {
        s += m[i, 2];
    }
    writefln("mean of the last column %.3f", s / m.nrow);

    setSeed(1);
    auto u = randomMatrix!runif(2, 2);                  // 2 x 2, uniform on [0, 1)
    writeln(u[0, 0] >= 0 && u[1, 1] < 1);
}
```

Output:

```
1000 x 3
mean of the last column 9.934
true
```
