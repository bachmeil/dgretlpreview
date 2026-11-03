# makeStreams

## Signature

```d
RandomStream[] makeStreams(size_t n, ulong seed = 0);
```

## Parameters

- `n`: the number of streams.
- `seed`: chooses where the family of streams starts. With the default 0 the first stream starts from the generator's original state; another value starts the family from a state derived from it.

## Return value

An array of `n` `RandomStream` objects, each a stream of L'Ecuyer's MRG31k3p generator. The streams start 2^134 steps apart, so they do not overlap in practice, and they do not depend on the library's own generator or on `setSeed`. The same `n` and `seed` give the same streams.

A stream draws through its methods (`streams[2].rnorm(100)`; see `RandomStream`), or is handed to a task by `parallelDraws`. A stream must be used by one thread at a time.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    auto streams = makeStreams(3, seed: 2024);
    auto a = streams[0].rnorm(3);
    auto b = streams[1].rnorm(3);
    writeln(a != b);                    // different streams, different draws

    auto again = makeStreams(3, seed: 2024);
    writeln(again[0].rnorm(3) == a);    // the same seed gives the same streams
}
```

Output:

```
true
true
```
