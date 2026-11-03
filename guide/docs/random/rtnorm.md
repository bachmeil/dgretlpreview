# rtnorm

## Signature

```d
double[] rtnorm(size_t n, double lower = -double.infinity, double upper = double.infinity, double mean = 0, double sd = 1);
```

## Parameters

- `n`: the number of draws wanted. `n = 0` returns an empty array.
- `lower`: the lower bound of the truncation. The default, `-double.infinity`, means no lower bound.
- `upper`: the upper bound of the truncation. It must be greater than `lower`. The default, `double.infinity`, means no upper bound.
- `mean`: the mean of the normal distribution before truncation. It must be finite. The default is 0.
- `sd`: the standard deviation of the normal distribution before truncation. It must be finite and positive. The default is 1.

## Return value

A `double[]` of length `n`. Each element is a draw from the normal distribution with the given `mean` and `sd`, conditioned to lie in [`lower`, `upper`]. Note that `mean` and `sd` describe the normal distribution before truncation, not the mean and standard deviation of the draws. All the draws are within the bounds.

The draws come from the library's random number generator, so they change from call to call; call `setSeed(seed)` first to make them reproducible. A `GretlException` is thrown for invalid parameters (the message names the function and says what is needed).

To get an `nrow` by `ncol` `GretlMatrix` of draws instead of an array, use `randomMatrix!rtnorm(nrow, ncol, ...)`. To draw from a particular parallel stream, call `rtnorm` as a method of the stream (`streams[3].rtnorm(100, ...)`, with the arguments after `n`).

## Example

```d
import std.stdio;
import std.algorithm: sum, minElement, maxElement;
import dgretl;

void main() {
    setSeed(3);
    auto x = rtnorm(100_000, 0, double.infinity);       // the standard normal truncated to be positive

    writefln("min %.3f, mean %.3f", x.minElement, x.sum / x.length);   // the mean is near sqrt(2/pi) = 0.798

    auto y = rtnorm(5, 1, 2, 1.5, 3);                   // N(1.5, 3^2) restricted to [1, 2]
    writeln(y.minElement >= 1 && y.maxElement <= 2);
}
```

Output:

```
min 0.000, mean 0.797
true
```
