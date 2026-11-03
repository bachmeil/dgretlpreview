# diff

## Signature

```d
GretlMatrix diff(X)(X x, int order = 1, int seasonal = 1);
GretlMatrix diff(Dataset d, string name, int order = 1, int seasonal = 1);
GretlMatrix diff(Dataset d, string[] names, int order = 1, int seasonal = 1);
```

## Parameters

- `x`: an array of numbers, or a `GretlMatrix` or `Submatrix`; each column is a series. (Missing values are allowed and propagate.)
- `order`: how many times to difference. The default is 1, the ordinary first difference. It must be at least 1.
- `seasonal`: the span `s` of the difference `x(t) - x(t - s)`: 1 (the default) for the ordinary difference, 4 for quarterly seasonal differences, 12 for monthly. It must be at least 1.
- `d`, `name` or `names`: in the `Dataset` forms, the series (or one column for each of several series) to difference.

## Return value

A `GretlMatrix` with the same number of rows and columns as `x`, holding the difference `x(t) - x(t - s)` taken `order` times, where `s` is `seasonal`. The first `order * seasonal` values of each column are NaN, because there is nothing to subtract.

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
    auto x = [1.0, 4.0, 9.0, 16.0, 25.0, 36.0];
    show(diff(x));                      // 3, 5, 7, ...
    show(diff(x, 2));                   // the second difference is constant
    show(diff(x, 1, 2));                // x(t) - x(t - 2)
}
```

Output:

```
       nan
    3.0000
    5.0000
    7.0000
    9.0000
   11.0000
       nan
       nan
    2.0000
    2.0000
    2.0000
    2.0000
       nan
       nan
    8.0000
   12.0000
   16.0000
   20.0000
```
