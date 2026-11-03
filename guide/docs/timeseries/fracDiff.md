# fracDiff

## Signature

```d
GretlMatrix fracDiff(X)(X x, double d);
GretlMatrix fracDiff(Dataset d, string name, double dd);
```

## Parameters

- `x`: the series: an array of numbers, or a `GretlMatrix` or `Submatrix` with one column. Missing values (NaN) are not allowed and throw a `GretlException`.
- `d`: the order of fractional differencing. It may be negative, which integrates the series fractionally. `d = 1` is the ordinary first difference.

## Return value

An n x 1 `GretlMatrix` with the fractional difference `(1 - L)^d x`, computed from the binomial expansion of the operator, so the series keeps its length (the first element is `x[0]`). The Dataset form is `fracDiff(Dataset d, string name, double dd)`, where `dd` is the order. A random walk differenced with `d = 1` is white noise, and with a smaller `d` is only partly differenced.

## Example

```d
import std.stdio;
import dgretl;

void show(M)(M m) {
    foreach (i; 0 .. m.nrow) {
        foreach (j; 0 .. m.ncol) {
            writef("%10.4f", m[i, j]);
        }
        writeln();
    }
}

void main() {
    auto x = [1.0, 3.0, 6.0, 10.0, 15.0];
    show(fracDiff(x, 1));                            // the ordinary difference (the first value is kept)
    show(fracDiff(x, 0.5));
    show(fracDiff(x, 0));                            // unchanged
}
```

Output:

```
    1.0000
    2.0000
    3.0000
    4.0000
    5.0000
    1.0000
    2.5000
    4.3750
    6.5625
    9.0234
    1.0000
    3.0000
    6.0000
   10.0000
   15.0000
```
