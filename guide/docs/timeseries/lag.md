# lag

## Signature

```d
GretlMatrix lag(X)(X x, int k = 1);
GretlMatrix lag(Dataset d, string name, int k = 1);
GretlMatrix lag(Dataset d, string[] names, int k = 1);
```

## Parameters

- `x`: an array of numbers, or a `GretlMatrix` or `Submatrix`; each column is a series. (Unlike the other functions of this section, missing values are allowed.)
- `k`: the number of observations to shift. A positive `k` gives a lag, a negative `k` a lead. The default is 1.
- `d`, `name` or `names`: in the `Dataset` forms, the series (or one column for each of several series) to shift.

## Return value

A `GretlMatrix` with the same number of rows and columns as `x`, where element `t` of each column is `x[t - k]` and the positions with no value to shift in (the first `k` for a lag, the last `-k` for a lead) are NaN.

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
    auto x = [1.0, 2.0, 3.0, 4.0, 5.0];
    show(lag(x));                       // the first value is missing
    show(lag(x, 2));
    show(lag(x, -1));                   // a lead
}
```

Output:

```
       nan
    1.0000
    2.0000
    3.0000
    4.0000
       nan
       nan
    1.0000
    2.0000
    3.0000
    2.0000
    3.0000
    4.0000
    5.0000
       nan
```
