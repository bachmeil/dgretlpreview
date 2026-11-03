# rged

## Signature

```d
double[] rged(size_t n, double nu);
```

## Parameters

- `n`: the number of draws wanted. `n = 0` returns an empty array.
- `nu`: the shape parameter, in the library's parameterization. It must be finite and positive. `nu = 2` is the normal distribution; smaller values give heavier tails.

## Return value

A `double[]` of length `n`. Each element is a draw from the generalized error distribution with shape `nu`, centred at 0.

The draws come from the library's random number generator, so they change from call to call; call `setSeed(seed)` first to make them reproducible. A `GretlException` is thrown for invalid parameters (the message names the function and says what is needed).

To get an `nrow` by `ncol` `GretlMatrix` of draws instead of an array, use `randomMatrix!rged(nrow, ncol, ...)`. To draw from a particular parallel stream, call `rged` as a method of the stream (`streams[3].rged(100, ...)`, with the arguments after `n`).

## Example

```d
import std.stdio;
import std.algorithm: sum;
import dgretl;

void main() {
    setSeed(5);
    auto x = rged(100_000, 2);          // the normal case

    double ss = 0;
    foreach (v; x) {
        ss += v * v;
    }
    writefln("mean %.3f, variance %.3f", x.sum / x.length, ss / x.length);
}
```

Output:

```
mean 0.002, variance 0.990
```
