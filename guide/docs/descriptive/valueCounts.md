# valueCounts

## Signature

```d
FrequencyTable valueCounts(S)(S x);
```

## Parameters

- `x`: a series: an array of numbers or one row or column of a matrix. Missing values are left out. At least 3 values are needed.

## Return value

A `FrequencyTable` (see `frequencies`) with `discrete = true`: one row for each distinct value of `x`, in increasing order. `midpoint` holds the value, `count` how often it occurs, `fraction` the proportion, and `cumulative` the running proportion.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    auto die = [1, 3, 2, 6, 3, 3, 5, 1, 6, 3, 2, 4];
    auto v = valueCounts(die);
    foreach (i; 0 .. v.count.nrow) {
        writefln("value %.0f: %.0f times (%.3f)", v.midpoint[i, 0], v.count[i, 0], v.fraction[i, 0]);
    }
}
```

Output:

```
value 1: 2 times (0.167)
value 2: 2 times (0.167)
value 3: 4 times (0.333)
value 4: 1 times (0.083)
value 5: 1 times (0.083)
value 6: 2 times (0.167)
```
