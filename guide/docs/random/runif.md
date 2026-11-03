# runif

## Signature

```d
double[] runif(size_t n, double min = 0, double max = 1);
```

## Parameters

- `n`: the number of draws wanted. `n = 0` returns an empty array.
- `min`: the lower end of the interval. It must be finite. The default is 0.
- `max`: the upper end of the interval. It must be finite and greater than `min`. The default is 1.

## Return value

A `double[]` of length `n`. Each element is an independent draw from the continuous uniform distribution on [`min`, `max`): `min` can occur, `max` does not.

The draws come from the library's random number generator, so they change from call to call; call `setSeed(seed)` first to make them reproducible. A `GretlException` is thrown for invalid parameters (the message names the function and says what is needed).

To get an `nrow` by `ncol` `GretlMatrix` of draws instead of an array, use `randomMatrix!runif(nrow, ncol, ...)`. To draw from a particular parallel stream, call `runif` as a method of the stream (`streams[3].runif(100, ...)`, with the arguments after `n`).

## Example

```d
import std.stdio;
import std.algorithm: sum, minElement, maxElement;
import dgretl;

void main() {
    setSeed(1);
    auto u = runif(100_000, 2, 5);      // uniform on [2, 5)

    writefln("min %.3f, max %.3f, mean %.3f", u.minElement, u.maxElement, u.sum / u.length);   // the mean is near 3.5
}
```

Output:

```
min 2.000, max 5.000, mean 3.500
```
