# rlogis

## Signature

```d
double[] rlogis(size_t n, double location = 0, double scale = 1);
```

## Parameters

- `n`: the number of draws wanted. `n = 0` returns an empty array.
- `location`: the mean (and median) of the distribution. It must be finite. The default is 0.
- `scale`: the scale parameter. It must be finite and positive. The default is 1; the standard deviation is `scale * pi / sqrt(3)`.

## Return value

A `double[]` of length `n`. Each element is a draw from the logistic distribution with the given location and scale.

The draws come from the library's random number generator, so they change from call to call; call `setSeed(seed)` first to make them reproducible. A `GretlException` is thrown for invalid parameters (the message names the function and says what is needed).

To get an `nrow` by `ncol` `GretlMatrix` of draws instead of an array, use `randomMatrix!rlogis(nrow, ncol, ...)`. To draw from a particular parallel stream, call `rlogis` as a method of the stream (`streams[3].rlogis(100, ...)`, with the arguments after `n`).

## Example

```d
import std.stdio;
import std.algorithm: sum;
import dgretl;

void main() {
    setSeed(5);
    auto x = rlogis(100_000, 2, 0.5);

    writefln("mean %.3f", x.sum / x.length);    // near 2
}
```

Output:

```
mean 2.001
```
