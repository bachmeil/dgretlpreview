# setSeed

## Signature

```d
void setSeed(ulong seed);
```

## Parameters

- `seed`: the seed of the library's random number generator. A seed of 0 asks the library to choose one from the system, so the sequence that follows is not reproducible.

## Return value

Nothing (`void`). The call resets the library's own generator (the one used by `rnorm`, `runif` and the other draws when they are called outside a stream), so the same seed followed by the same sequence of calls gives the same draws. It does not affect the parallel streams made by `makeStreams`, which have their own seeds.

Use `getSeed` to read the seed back.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    setSeed(42);
    auto a = rnorm(3);

    setSeed(42);                        // start again from the same seed
    auto b = rnorm(3);

    writeln(a == b);                    // the same draws
    writeln(rnorm(3) == a);             // the generator has moved on: different draws
}
```

Output:

```
true
false
```
