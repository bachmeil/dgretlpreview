# rbinom

## Signature

```d
double[] rbinom(size_t n, int size, double p);
```

## Parameters

- `n`: the number of draws wanted. `n = 0` returns an empty array.
- `size`: the number of trials. It must be 0 or greater.
- `p`: the probability of success on each trial, from 0 to 1.

## Return value

A `double[]` of length `n`. Each element is the number of successes in `size` independent trials, a whole number from 0 to `size`; the draws are stored as doubles, not integers (the mean is `size * p`).

The draws come from the library's random number generator, so they change from call to call; call `setSeed(seed)` first to make them reproducible. A `GretlException` is thrown for invalid parameters (the message names the function and says what is needed).

To get an `nrow` by `ncol` `GretlMatrix` of draws instead of an array, use `randomMatrix!rbinom(nrow, ncol, ...)`. To draw from a particular parallel stream, call `rbinom` as a method of the stream (`streams[3].rbinom(100, ...)`, with the arguments after `n`).

## Example

```d
import std.stdio;
import std.algorithm: sum, minElement, maxElement;
import dgretl;

void main() {
    setSeed(5);
    auto x = rbinom(100_000, 10, 0.3);

    writefln("mean %.3f", x.sum / x.length);    // near 3
    writeln(x.minElement >= 0 && x.maxElement <= 10);
}
```

Output:

```
mean 3.000
true
```
