# t

## Signature

```d
GretlMatrix t(M)(M m);
```

## Parameters

- `m`: a `GretlMatrix` or a `Submatrix` (a view of part of a matrix). It is not changed.

## Return value

A new `GretlMatrix`, the transpose of `m`: an `ncol` x `nrow` matrix whose element `[j, i]` is `m[i, j]`. The result is a copy, so changing it does not change `m`. It can also be called as a method: `m.t`.

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
    auto m = cbind([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]);   // 3 x 2
    show(m);
    show(t(m));                                          // 2 x 3
    show(m.t * m);                                       // m'm
}
```

Output:

```
   1.0000   4.0000
   2.0000   5.0000
   3.0000   6.0000
   1.0000   2.0000   3.0000
   4.0000   5.0000   6.0000
  14.0000  32.0000
  32.0000  77.0000
```
