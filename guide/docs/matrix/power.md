# power

## Signature

```d
GretlMatrix power(M)(M m, long n);
```

## Parameters

- `m`: a square matrix: a `GretlMatrix` or a `Submatrix` (a view of part of a matrix). It is not changed.
- `n`: the exponent, a whole number. A negative `n` gives a power of the inverse.

## Return value

A new `GretlMatrix`, the matrix power `m * m * ... * m` with `n` factors, computed by repeated squaring. `power(m, 2)` is `m * m`, `power(m, 0)` is the identity, and `power(m, -k)` is `power(inv(m), k)`. The matrix must be square; otherwise a `GretlException` is thrown. (The `^^` operator is refused for matrices because it would be ambiguous; use `m.el ^^ a` for element powers.) It can also be called as a method: `m.power(n)`.

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
    auto m = cbind([1.0, 1.0], [1.0, 0.0]);            // the Fibonacci matrix
    show(power(m, 10));                                 // contains Fibonacci numbers
    show(power(m, 0));
    show(power(m, -1));
}
```

Output:

```
  89.0000  55.0000
  55.0000  34.0000
   1.0000   0.0000
   0.0000   1.0000
   0.0000   1.0000
   1.0000  -1.0000
```
