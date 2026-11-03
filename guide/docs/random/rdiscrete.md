# rdiscrete

## Signature

```d
size_t[] rdiscrete(size_t n, const(double)[] probs);
```

## Parameters

- `n`: the number of draws wanted. `n = 0` returns an empty array.
- `probs`: the probabilities of the values 0, 1, ..., `probs.length - 1`. Each must be finite and not negative, and they must add to 1 (to within 1e-9). There must be at least one.

## Return value

A `size_t[]` of length `n`. Each element is one of the values `0 .. probs.length - 1`, drawn with the probabilities in `probs` by the alias method. The values start at 0 (the library's own function counts from 1; the front end subtracts one).

The draws come from the library's random number generator, so they change from call to call; call `setSeed(seed)` first to make them reproducible. A `GretlException` is thrown for invalid parameters (the message names the function and says what is needed).

The result is not a `double[]`, so `randomMatrix` does not apply; the stream method works (`streams[3].rdiscrete(100, probs)`).

## Example

```d
import std.stdio;
import dgretl;

void main() {
    setSeed(5);
    auto x = rdiscrete(100_000, [0.2, 0.5, 0.3]);

    size_t[3] count;
    foreach (v; x) {
        count[v]++;
    }
    foreach (k, c; count) {
        writefln("value %d: %.3f", k, c / cast(double) x.length);
    }
}
```

Output:

```
value 0: 0.199
value 1: 0.499
value 2: 0.302
```
