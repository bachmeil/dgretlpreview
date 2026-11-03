# rt

## Signature

```d
double[] rt(size_t n, double df);
```

## Parameters

- `n`: the number of draws wanted. `n = 0` returns an empty array.
- `df`: the degrees of freedom. It must be finite and positive, and need not be an integer.

## Return value

A `double[]` of length `n`. Each element is a draw from Student's t distribution with `df` degrees of freedom. The distribution is symmetric about 0, with heavier tails than the normal distribution; its variance is `df / (df - 2)` when `df > 2`.

The draws come from the library's random number generator, so they change from call to call; call `setSeed(seed)` first to make them reproducible. A `GretlException` is thrown for invalid parameters (the message names the function and says what is needed).

To get an `nrow` by `ncol` `GretlMatrix` of draws instead of an array, use `randomMatrix!rt(nrow, ncol, ...)`. To draw from a particular parallel stream, call `rt` as a method of the stream (`streams[3].rt(100, ...)`, with the arguments after `n`).

## Example

```d
import std.stdio;
import std.algorithm: sum;
import dgretl;

void main() {
    setSeed(5);
    auto x = rt(100_000, 10);

    double ss = 0;
    foreach (v; x) {
        ss += v * v;
    }
    writefln("mean %.3f, variance %.3f", x.sum / x.length, ss / x.length);   // variance near 10/8 = 1.25
}
```

Output:

```
mean -0.005, variance 1.245
```
