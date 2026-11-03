# rbeta

## Signature

```d
double[] rbeta(size_t n, double a, double b);
```

## Parameters

- `n`: the number of draws wanted. `n = 0` returns an empty array.
- `a`: the first shape parameter. It must be finite and positive.
- `b`: the second shape parameter. It must be finite and positive.

## Return value

A `double[]` of length `n`. Each element is a draw from the beta distribution with shape parameters `a` and `b`, a number between 0 and 1 with mean `a / (a + b)`.

The draws come from the library's random number generator, so they change from call to call; call `setSeed(seed)` first to make them reproducible. A `GretlException` is thrown for invalid parameters (the message names the function and says what is needed).

To get an `nrow` by `ncol` `GretlMatrix` of draws instead of an array, use `randomMatrix!rbeta(nrow, ncol, ...)`. To draw from a particular parallel stream, call `rbeta` as a method of the stream (`streams[3].rbeta(100, ...)`, with the arguments after `n`).

## Example

```d
import std.stdio;
import std.algorithm: sum;
import dgretl;

void main() {
    setSeed(5);
    auto x = rbeta(100_000, 2, 6);

    writefln("mean %.3f", x.sum / x.length);    // near 2/8 = 0.25
}
```

Output:

```
mean 0.250
```
