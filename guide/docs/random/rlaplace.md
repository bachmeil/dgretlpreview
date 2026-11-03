# rlaplace

## Signature

```d
double[] rlaplace(size_t n, double mu = 0, double b = 1);
```

## Parameters

- `n`: the number of draws wanted. `n = 0` returns an empty array.
- `mu`: the mean (and median) of the distribution. It must be finite. The default is 0.
- `b`: the scale parameter. It must be finite and positive. The default is 1; the variance is `2 * b^2`.

## Return value

A `double[]` of length `n`. Each element is a draw from the Laplace (double exponential) distribution with mean `mu` and scale `b`.

The draws come from the library's random number generator, so they change from call to call; call `setSeed(seed)` first to make them reproducible. A `GretlException` is thrown for invalid parameters (the message names the function and says what is needed).

To get an `nrow` by `ncol` `GretlMatrix` of draws instead of an array, use `randomMatrix!rlaplace(nrow, ncol, ...)`. To draw from a particular parallel stream, call `rlaplace` as a method of the stream (`streams[3].rlaplace(100, ...)`, with the arguments after `n`).

## Example

```d
import std.stdio;
import std.algorithm: sum;
import dgretl;

void main() {
    setSeed(5);
    auto x = rlaplace(100_000, 1, 2);

    double ss = 0;
    foreach (v; x) {
        ss += (v - 1) ^^ 2;
    }
    writefln("mean %.3f, variance %.3f", x.sum / x.length, ss / x.length);   // near 1 and 8
}
```

Output:

```
mean 0.997, variance 8.012
```
