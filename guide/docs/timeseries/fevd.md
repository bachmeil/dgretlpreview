# fevd

## Signature

```d
Fevd fevd(F)(ref F fit, int horizon = 0);

struct Fevd {
    size_t neqns;
    int horizon;
    GretlMatrix[] share;
    GretlMatrix at(size_t targ, size_t shock);
}
```

## Parameters

- `fit`: a `VarFit` or a `VecmFit`.
- `horizon`: the number of periods. The default 0 takes the library's default horizon.

## Return value

An `Fevd` struct with the forecast error variance decomposition. `at(targ, shock)` gives a `horizon` x 1 matrix: the fraction of the variance of the forecast error of variable `targ` (at each horizon, period 1 first) that comes from the orthogonalized shock to variable `shock`. For each `targ` and horizon the fractions over the shocks add to 1. (`share[targ * neqns + shock]` is the same matrix.)

## Example

```d
import std.stdio;
import dgretl;

void main() {
    setSeed(31);
    auto n = 200;
    auto e1 = rnorm(n);
    auto e2 = rnorm(n);
    auto y1 = new double[n];
    auto y2 = new double[n];
    y1[0] = e1[0];
    y2[0] = e2[0];
    foreach (t; 1 .. n) {
        y1[t] = 0.5 * y1[t - 1] + 0.2 * y2[t - 1] + e1[t];
        y2[t] = 0.3 * y1[t - 1] + 0.4 * y2[t - 1] + e2[t];
    }
    auto fit = var(cbind(y1, y2), 1);

    auto d = fevd(fit, 5);
    foreach (h; 0 .. 5) {
        writefln("horizon %d: y2's variance from shocks to y1 %.3f, to y2 %.3f", h + 1, d.at(1, 0)[h, 0], d.at(1, 1)[h, 0]);
    }
}
```

Output:

```
horizon 1: y2's variance from shocks to y1 0.002, to y2 0.998
horizon 2: y2's variance from shocks to y1 0.052, to y2 0.948
horizon 3: y2's variance from shocks to y1 0.091, to y2 0.909
horizon 4: y2's variance from shocks to y1 0.112, to y2 0.888
horizon 5: y2's variance from shocks to y1 0.122, to y2 0.878
```
