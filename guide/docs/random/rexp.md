# rexp

## Signature

```d
double[] rexp(size_t n, double mean = 1);
```

## Parameters

- `n`: the number of draws wanted. `n = 0` returns an empty array.
- `mean`: the mean of the exponential distribution (the reciprocal of the rate). It must be finite and positive. The default is 1.

## Return value

A `double[]` of length `n`. Each element is a draw from the exponential distribution with the given mean. All draws are positive. Note that the parameter is the mean, not the rate as in R's `rexp`: for rate `r` pass `1 / r`.

The draws come from the library's random number generator, so they change from call to call; call `setSeed(seed)` first to make them reproducible. A `GretlException` is thrown for invalid parameters (the message names the function and says what is needed).

To get an `nrow` by `ncol` `GretlMatrix` of draws instead of an array, use `randomMatrix!rexp(nrow, ncol, ...)`. To draw from a particular parallel stream, call `rexp` as a method of the stream (`streams[3].rexp(100, ...)`, with the arguments after `n`).

## Example

```d
import std.stdio;
import std.algorithm: sum;
import dgretl;

void main() {
    setSeed(5);
    auto x = rexp(100_000, 0.5);

    writefln("mean %.3f", x.sum / x.length);    // near 0.5
}
```

Output:

```
mean 0.499
```
