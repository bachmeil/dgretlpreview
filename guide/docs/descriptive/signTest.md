# signTest

## Signature

```d
SignTest signTest(S, T)(S x, T y);
SignTest signTest(Dataset d, string first, string second);

struct SignTest {
    size_t n, w;
    double pvalueUpper, pvalueLower, pvalue;
    string report;
}
```

## Parameters

- `x`, `y`: two paired series of the same length (arrays or a row or column of a matrix), or in the `Dataset` form the names of two series. Different lengths throw a `GretlException`. Pairs where either value is missing, or where the two values are equal, are left out.

## Return value

A `SignTest` struct, the sign test of the null hypothesis that the paired series have the same median, so that `x > y` and `x < y` are equally likely:

- `n`: the number of pairs that count (valid and different).
- `w`: of those, the number with `x > y`.
- `pvalueUpper`: `P(W' >= w)` for `W'` binomial with `n` trials and probability 1/2 (the one the library records); `pvalueLower`: `P(W' <= w)`.
- `pvalue`: the two-sided p-value, `min(1, 2 min(pvalueUpper, pvalueLower))`.
- `report`: the text the library prints for the test.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    auto x = [8.2, 7.9, 9.1, 8.8, 7.5, 9.4, 8.0, 8.9, 9.6, 7.7];
    auto y = [7.9, 7.8, 8.5, 8.9, 7.1, 8.8, 8.1, 8.2, 9.0, 7.2];

    auto t = signTest(x, y);
    writefln("%d pairs counted, x > y in %d", t.n, t.w);
    writefln("one-sided p = %.4f, two-sided p = %.4f", t.pvalueUpper, t.pvalue);
}
```

Output:

```
10 pairs counted, x > y in 8
one-sided p = 0.0547, two-sided p = 0.1094
```
