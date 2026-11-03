# eye

## Signature

```d
GretlMatrix eye(size_t n);
GretlMatrix eye(size_t nr, size_t nc);
```

## Parameters

- `n`: the size of a square identity matrix.
- `nr`, `nc`: the number of rows and columns of the (possibly rectangular) result.

## Return value

A new `GretlMatrix`. `eye(n)` is the `n` x `n` identity matrix. `eye(nr, nc)` has `nr` rows and `nc` columns, with ones on the main diagonal and zeros everywhere else (so it is the identity padded with zeros when it is not square). A size of 0 gives an empty matrix.

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
    auto I = eye(3);
    show(I);

    auto E = eye(2, 4);
    show(E);
}
```

Output:

```
   1.0000   0.0000   0.0000
   0.0000   1.0000   0.0000
   0.0000   0.0000   1.0000
   1.0000   0.0000   0.0000   0.0000
   0.0000   1.0000   0.0000   0.0000
```
