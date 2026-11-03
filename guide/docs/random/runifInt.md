# runifInt

## Signature

```d
int[] runifInt(size_t n, int min, int max, bool withoutReplacement = false);
```

## Parameters

- `n`: the number of draws wanted. `n = 0` returns an empty array.
- `min`: the smallest integer that can be drawn.
- `max`: the largest integer that can be drawn. It must not be below `min`.
- `withoutReplacement`: if `true`, the draws are distinct (a random sample from `min .. max`), so `n` must not exceed `max - min + 1`. The default is `false`, which draws with replacement.

## Return value

An `int[]` of length `n` (not a `double[]`, unlike the continuous draws). Each element is an integer from `min` to `max`, both included, with equal probability. With `withoutReplacement` no value occurs twice, in random order.

The draws come from the library's random number generator, so they change from call to call; call `setSeed(seed)` first to make them reproducible. A `GretlException` is thrown for invalid parameters (the message names the function and says what is needed). Asking for more distinct values than the range holds also throws.

`randomMatrix` does not apply to this function, because it returns `int[]`; the stream method (`streams[3].runifInt(...)`) does.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    import std.algorithm: sort;

    setSeed(7);
    auto dice = runifInt(10, 1, 6);                  // ten rolls of a die
    writeln(dice);

    auto lottery = runifInt(6, 1, 49, true);         // six distinct numbers from 1 .. 49
    sort(lottery);
    writeln(lottery);
}
```

Output:

```
[5, 2, 3, 2, 2, 2, 6, 5, 1, 6]
[11, 15, 27, 30, 35, 44]
```
