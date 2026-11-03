# RandomStream

## Signature

```d
final class RandomStream : Stream;
```

## Parameters

- `state()`: returns the six integers of the stream's state (`int[6]`), from which another stream can be started.
- `uniform()`: the next U(0, 1) draw, a `double` never equal to 0 or 1 (31 bits of resolution).
- `next64()`: 64 random bits as a `ulong`, from three steps of the generator.
- `use(f)`: runs the delegate `f` with this stream as the source of every random draw of the library, in the current thread. The library's global generator is not used meanwhile.
- any distribution function of this module (`runif`, `runifInt`, `rnorm`, `rtnorm`, `rgamma`, `rchisq`, `rt`, `rF`, `rbinom`, `rbetabinom`, `rpois`, `rweibull`, `rexp`, `rlogis`, `rged`, `rlaplace`, `rbeta`, `rdiscrete`, `rdirichlet`, `rinvwishart`, `randomMatrix`) can be called as a method, with the same arguments, and draws from the stream: `s.rnorm(100, 1, 2)`.

## Return value

A `RandomStream` is made by `makeStreams`, not constructed directly. The stream methods return what the functions of the same name return. Successive calls continue the same stream, so two calls give different draws, and the whole sequence depends only on the stream's starting state. A stream must be used by one thread at a time; to use several in parallel give each task its own, as `parallelDraws` does.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    auto s = makeStreams(2, seed: 11)[0];

    auto x = s.rnorm(5, 10, 2);         // five N(10, 2^2) draws from this stream
    auto g = s.rgamma(3, 2.0);          // then three gamma draws
    writeln(x.length, " ", g.length);

    writeln(s.uniform());               // the next U(0, 1) draw

    s.use({
        auto y = runif(3);              // inside use, the ordinary function draws from s
        writeln(y.length);
    });
}
```

Output:

```
5 3
0.689819
3
```
