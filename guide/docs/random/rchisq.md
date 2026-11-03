# rchisq

## Signature

```d
double[] rchisq(size_t n, int df);
```

## Parameters

- `n`: the number of draws wanted. `n = 0` returns an empty array.
- `df`: the degrees of freedom, an integer greater than 0.

## Return value

A `double[]` of length `n`. Each element is a draw from the chi-square distribution with `df` degrees of freedom (mean `df`, variance `2 * df`). All draws are positive.

The draws come from the library's random number generator, so they change from call to call; call `setSeed(seed)` first to make them reproducible. A `GretlException` is thrown for invalid parameters (the message names the function and says what is needed).

To get an `nrow` by `ncol` `GretlMatrix` of draws instead of an array, use `randomMatrix!rchisq(nrow, ncol, ...)`. To draw from a particular parallel stream, call `rchisq` as a method of the stream (`streams[3].rchisq(100, ...)`, with the arguments after `n`).

## Example

```d
import std.stdio;
import std.algorithm: sum;
import dgretl;

void main() {
    setSeed(5);
    auto x = rchisq(100_000, 4);

    writefln("mean %.3f", x.sum / x.length);    // near 4
}
```

Output:

```
mean 3.998
```
