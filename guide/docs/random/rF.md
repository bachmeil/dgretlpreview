# rF

## Signature

```d
double[] rF(size_t n, int df1, int df2);
```

## Parameters

- `n`: the number of draws wanted. `n = 0` returns an empty array.
- `df1`: the numerator degrees of freedom, an integer greater than 0.
- `df2`: the denominator degrees of freedom, an integer greater than 0.

## Return value

A `double[]` of length `n`. Each element is a draw from the F distribution with `df1` and `df2` degrees of freedom. All draws are positive; the mean is `df2 / (df2 - 2)` when `df2 > 2`.

The draws come from the library's random number generator, so they change from call to call; call `setSeed(seed)` first to make them reproducible. A `GretlException` is thrown for invalid parameters (the message names the function and says what is needed).

To get an `nrow` by `ncol` `GretlMatrix` of draws instead of an array, use `randomMatrix!rF(nrow, ncol, ...)`. To draw from a particular parallel stream, call `rF` as a method of the stream (`streams[3].rF(100, ...)`, with the arguments after `n`).

## Example

```d
import std.stdio;
import std.algorithm: sum;
import dgretl;

void main() {
    setSeed(5);
    auto x = rF(100_000, 3, 20);

    writefln("mean %.3f", x.sum / x.length);    // near 20/18 = 1.111
}
```

Output:

```
mean 1.112
```
