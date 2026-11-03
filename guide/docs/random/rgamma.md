# rgamma

## Signature

```d
double[] rgamma(size_t n, double shape, double scale = 1);
```

## Parameters

- `n`: the number of draws wanted. `n = 0` returns an empty array.
- `shape`: the shape parameter. It must be finite and positive.
- `scale`: the scale parameter. It must be finite and positive. The default is 1.

## Return value

A `double[]` of length `n`. Each element is a draw from the gamma distribution with the given shape and scale, whose mean is `shape * scale` and variance `shape * scale^2`. All draws are positive.

The draws come from the library's random number generator, so they change from call to call; call `setSeed(seed)` first to make them reproducible. A `GretlException` is thrown for invalid parameters (the message names the function and says what is needed).

To get an `nrow` by `ncol` `GretlMatrix` of draws instead of an array, use `randomMatrix!rgamma(nrow, ncol, ...)`. To draw from a particular parallel stream, call `rgamma` as a method of the stream (`streams[3].rgamma(100, ...)`, with the arguments after `n`).

## Example

```d
import std.stdio;
import std.algorithm: sum;
import dgretl;

void main() {
    setSeed(5);
    auto x = rgamma(100_000, 3, 2);     // shape 3, scale 2: mean 6, variance 12

    const mean = x.sum / x.length;
    double ss = 0;
    foreach (v; x) {
        ss += (v - mean) ^^ 2;
    }
    writefln("mean %.3f, variance %.3f", mean, ss / (x.length - 1));
}
```

Output:

```
mean 5.998, variance 11.997
```
