# rpois

## Signature

```d
double[] rpois(size_t n, double lambda);
double[] rpois(const(double)[] lambdas);
```

## Parameters

- `n` (first form): the number of draws wanted. `n = 0` returns an empty array.
- `lambda` (first form): the mean of the Poisson distribution. It must be finite and not negative. A mean of exactly 0 gives NaN (a missing value) instead of the draw 0, which is how the library behaves (see `bugs.md`, item 133); use a positive mean.
- `lambdas` (second form): an array of means, one for each draw. Each must be finite and not negative; as in the first form, an element equal to 0 gives NaN.

## Return value

A `double[]`. In the first form it has length `n` and each element is a draw from the Poisson distribution with mean `lambda`. In the second form it has the length of `lambdas`, and element `i` is a draw from the Poisson distribution with mean `lambdas[i]`. The draws are whole numbers stored as doubles.

The draws come from the library's random number generator, so they change from call to call; call `setSeed(seed)` first to make them reproducible. A `GretlException` is thrown for invalid parameters (the message names the function and says what is needed).

The first form also works with `randomMatrix!rpois(nrow, ncol, lambda)` and as a stream method (`streams[3].rpois(100, 4.0)`).

## Example

```d
import std.stdio;
import std.algorithm: sum;
import dgretl;

void main() {
    setSeed(5);
    auto x = rpois(100_000, 4.5);
    writefln("mean %.3f", x.sum / x.length);    // near 4.5

    auto y = rpois([0.5, 1, 10, 100]);          // one draw for each mean
    writeln(y);
}
```

Output:

```
mean 4.502
[0, 0, 12, 108]
```
