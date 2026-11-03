# parallelDraws

## Signature

```d
T[] parallelDraws(T)(RandomStream[] streams, scope T delegate(size_t) f);
T[] parallelDraws(T)(RandomStream[] streams, scope T delegate(RandomStream, size_t) f);
```

## Parameters

- `streams`: the streams, normally from `makeStreams`. One task is run for each stream.
- `f`: the task. In the first form it takes the index `i` of the stream; in the second it takes the stream and its index. Whatever the distribution functions of this module (`rnorm`, `rgamma`, ...) draw inside `f` comes from stream `i`. The result type `T` is that of `f`.

## Return value

A `T[]` of length `streams.length`: element `i` is the result of `f` for stream `i`. The tasks run in parallel on the task pool of `std.parallelism`, but the result does not depend on the number of processors or on how the tasks were scheduled, so a simulation is reproducible from the seed given to `makeStreams`. The library's own generator is not used or changed.

## Example

```d
import std.stdio;
import std.algorithm: sum;
import dgretl;

void main() {
    auto streams = makeStreams(4, seed: 2024);

    // each task draws 100,000 normals from its own stream and returns their mean
    auto means = parallelDraws(streams, (size_t i) => rnorm(100_000).sum / 100_000);
    foreach (i, m; means) {
        writefln("stream %d: %.3f", i, m);
    }

    // the same with the stream passed in
    auto first = parallelDraws(makeStreams(4, seed: 2024), (RandomStream s, size_t i) => s.uniform());
    writeln(first.length);

    // repeating the run gives the same results
    auto again = parallelDraws(makeStreams(4, seed: 2024), (size_t i) => rnorm(100_000).sum / 100_000);
    writeln(again == means);
}
```

Output:

```
stream 0: 0.003
stream 1: 0.003
stream 2: -0.002
stream 3: 0.000
4
true
```
