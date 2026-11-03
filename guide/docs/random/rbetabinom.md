# rbetabinom

## Signature

```d
double[] rbetabinom(size_t n, int size, double a, double b);
```

## Parameters

- `n`: the number of draws wanted. `n = 0` returns an empty array.
- `size`: the number of trials. It must be 0 or greater.
- `a`: the first shape parameter of the beta distribution from which the success probability is drawn. It must be positive.
- `b`: the second shape parameter. It must be positive.

## Return value

A `double[]` of length `n`. Each element is the number of successes in `size` trials, whole numbers from 0 to `size` stored as doubles, where the success probability is itself drawn from the beta(`a`, `b`) distribution (once for each draw). The mean is `size * a / (a + b)`, and the variance is larger than that of a binomial with the same mean.

The draws come from the library's random number generator, so they change from call to call; call `setSeed(seed)` first to make them reproducible. A `GretlException` is thrown for invalid parameters (the message names the function and says what is needed).

To get an `nrow` by `ncol` `GretlMatrix` of draws instead of an array, use `randomMatrix!rbetabinom(nrow, ncol, ...)`. To draw from a particular parallel stream, call `rbetabinom` as a method of the stream (`streams[3].rbetabinom(100, ...)`, with the arguments after `n`).

## Example

```d
import std.stdio;
import std.algorithm: sum;
import dgretl;

void main() {
    setSeed(5);
    auto x = rbetabinom(100_000, 10, 2, 3);

    writefln("mean %.3f", x.sum / x.length);    // near 10 * 2/5 = 4
}
```

Output:

```
mean 3.999
```
