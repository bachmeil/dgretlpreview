# getSeed

## Signature

```d
ulong getSeed();
```

## Parameters

None.

## Return value

A `ulong`: the seed of the library's random number generator. After `setSeed(s)` with `s` nonzero it is `s`. If the generator was seeded from the system (a seed of 0, or no call to `setSeed` at all), it is the seed that was chosen, which can be passed to `setSeed` to repeat the run.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    setSeed(2024);
    writeln(getSeed());

    setSeed(0);                         // a seed chosen from the system
    const s = getSeed();
    auto a = rnorm(2);
    setSeed(s);
    writeln(rnorm(2) == a);             // the run can be repeated
}
```

Output:

```
2024
true
```
