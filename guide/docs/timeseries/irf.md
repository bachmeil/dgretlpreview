# irf

## Signature

```d
Irf irf(F)(ref F fit, int periods = 0, double alpha = 0.0, long seed = -1);

struct Irf {
    size_t neqns;
    int periods;
    double alpha;
    GretlMatrix[] response;
    GretlMatrix at(size_t targ, size_t shock);
}
```

## Parameters

- `fit`: a `VarFit` or a `VecmFit`.
- `periods`: the number of periods of the response. The default 0 takes the library's default horizon for the frequency of the data.
- `alpha`: if not 0, add bootstrap bands covering `1 - alpha` (percentile intervals from the library's residual bootstrap). It must be 0 or between 0.01 and 0.6. The default is 0 (no bands).
- `seed`: a seed that makes the bootstrap reproducible. The default -1 leaves the library's generator as it is.

## Return value

An `Irf` struct with the orthogonalized impulse responses: the response of each variable to a one-standard-deviation shock to each variable, using the Cholesky factor of the residual covariance matrix with the variables in the order of the VAR. `at(targ, shock)` gives the response of variable `targ` to a shock to variable `shock` as a `periods` x 1 matrix (period 0 first), or a `periods` x 3 matrix with the point estimate, the lower and the upper band when `alpha` was given (`response[targ * neqns + shock]` is the same). `neqns`, `periods` and `alpha` record what was asked.

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

    auto r = irf(fit, 6);
    auto resp = r.at(1, 0);                          // the response of y2 to a shock to y1
    foreach (h; 0 .. 6) {
        writefln("period %d: %7.4f", h, resp[h, 0]);
    }

    auto banded = irf(fit, 6, 0.1, 99);              // with 90% bootstrap bands
    auto b = banded.at(1, 0);
    writefln("%d x %d: point %.3f, band [%.3f, %.3f] at period 1", b.nrow, b.ncol, b[1, 0], b[1, 1], b[1, 2]);
}
```

Output:

```
period 0: -0.0411
period 1:  0.2605
period 2:  0.2491
period 3:  0.1905
period 4:  0.1379
period 5:  0.0981
6 x 3: point 0.261, band [0.150, 0.372] at period 1
```
