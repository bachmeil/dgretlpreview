# rweibull

## Signature

```d
double[] rweibull(size_t n, double shape, double scale = 1);
```

## Parameters

- `n`: the number of draws wanted. `n = 0` returns an empty array.
- `shape`: the shape parameter. It must be finite and positive.
- `scale`: the scale parameter. It must be finite and positive. The default is 1.

## Return value

A `double[]` of length `n`. Each element is a draw from the Weibull distribution with the given shape and scale, whose distribution function is `1 - exp(-(x / scale)^shape)`. All draws are positive.

The draws come from the library's random number generator, so they change from call to call; call `setSeed(seed)` first to make them reproducible. A `GretlException` is thrown for invalid parameters (the message names the function and says what is needed).

To get an `nrow` by `ncol` `GretlMatrix` of draws instead of an array, use `randomMatrix!rweibull(nrow, ncol, ...)`. To draw from a particular parallel stream, call `rweibull` as a method of the stream (`streams[3].rweibull(100, ...)`, with the arguments after `n`).

## Example

```d
import std.stdio;
import std.algorithm: sum;
import dgretl;

void main() {
    setSeed(5);
    auto x = rweibull(100_000, 1, 2);   // shape 1 is the exponential distribution with mean 2

    writefln("mean %.3f", x.sum / x.length);    // near 2
}
```

Output:

```
mean 1.997
```
