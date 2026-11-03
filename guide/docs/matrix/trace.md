# trace

## Signature

```d
double trace(M)(M m);
```

## Parameters

- `m`: a `GretlMatrix` or a `Submatrix` (a view of part of a matrix). It is not changed.

## Return value

A `double`, the sum of the diagonal elements of the square matrix `m`. The matrix must be square; otherwise a `GretlException` is thrown. It can also be called as a method: `m.trace`.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    auto m = cbind([1.0, 2.0], [3.0, 4.0]);
    writeln(trace(m));                  // 1 + 4
}
```

Output:

```
5
```
