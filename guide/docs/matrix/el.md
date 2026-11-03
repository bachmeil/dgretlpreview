# el

## Signature

```d
Elements el(M)(M m);
```

## Parameters

- `m`: a `GretlMatrix` or a `Submatrix`.

## Return value

A reference to all the elements of `m`, used in two ways.

- Iteration: `foreach (val; m.el)` visits every element in memory order (column by column, so the order is not part of the contract), `foreach (ref val; m.el) val *= 2` changes `m`, and `foreach (e, val; m.el)` also gives `e.row` and `e.col`.
- Element-by-element arithmetic: `m.el * mm`, `m.el / mm`, `m.el * a`, `m.el ^^ 2`, `a ^^ m.el`, `m.el ^^ mm.el` and so on, each returning a new `GretlMatrix`. (`+` and `-` are already element by element, and `m * mm` alone is the matrix product; `m ^^ n` alone is refused, see `power`.)

## Example

```d
import std.stdio;
import dgretl;

void show(M)(M m) {
    foreach (i; 0 .. m.nrow) {
        foreach (j; 0 .. m.ncol) {
            writef("%9.4f", m[i, j]);
        }
        writeln();
    }
}

void main() {
    auto m = cbind([1.0, 2.0], [3.0, 4.0]);
    foreach (ref v; m.el) {
        v *= 2;                               // changes m
    }
    show(m);

    foreach (e, v; m.el) {
        if (v > 5) writefln("[%d, %d] = %.0f", e.row, e.col, v);
    }

    show(m.el * m);                           // squares each element of m (times m, element by element)
    show(m.el ^^ 0.5);                        // square root of each element
    show(2.0 ^^ m.el);
}
```

Output:

```
   2.0000   6.0000
   4.0000   8.0000
[0, 1] = 6
[1, 1] = 8
   4.0000  36.0000
  16.0000  64.0000
   1.4142   2.4495
   2.0000   2.8284
   4.0000  64.0000
  16.0000 256.0000
```
