# rnorm

## Signature

```d
double[] rnorm(size_t n, double mean = 0, double sd = 1);
```

## Parameters

- `n`: the number of draws wanted. `n = 0` returns an empty array.
- `mean`: the mean of the distribution. It must be finite. The default is 0.
- `sd`: the standard deviation of the distribution. It must be finite and strictly positive. The default is 1.

## Return value

A `double[]` of length `n`. Each element is an independent draw from the normal distribution with the given mean and standard deviation, generated with the Ziggurat method. There are no missing values in the result.

The draws come from the library's random number generator, so they change from call to call; call `setSeed(seed)` first to make them reproducible. A `GretlException` is thrown for invalid parameters (the message names the function and says what is needed).

To get an `nrow` by `ncol` `GretlMatrix` of draws instead of an array, use `randomMatrix!rnorm(nrow, ncol, ...)`. To draw from a particular parallel stream, call `rnorm` as a method of the stream (`streams[3].rnorm(100, ...)`, with the arguments after `n`).

## Example

```d
import std.stdio;
import std.algorithm: sum;
import dgretl;

void main() {
    setSeed(42);
    auto x = rnorm(10_000, 5, 2);       // 10,000 draws from N(5, 2^2)

    double sumSq = 0;
    foreach (v; x) {
        sumSq += v * v;
    }
    const mean = x.sum / x.length;
    const sd = (sumSq / x.length - mean * mean) ^^ 0.5;
    writefln("draws: %d, mean %.3f, sd %.3f", x.length, mean, sd);

    auto z = rnorm(3);                  // standard normal: mean 0, sd 1
    writeln(z.length);

    auto m = randomMatrix!rnorm(100, 4, 1.0, 0.5);   // a 100 x 4 matrix of N(1, 0.5^2) draws
    writefln("%d x %d", m.nrow, m.ncol);

    try {
        rnorm(5, 0, -1);
    } catch (GretlException e) {
        writeln("error: ", e.msg);
    }
}
```

Output:

```
draws: 10000, mean 5.009, sd 2.004
3
100 x 4
error: rnorm: need a finite mean and a positive sd
```
