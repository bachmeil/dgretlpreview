# cond

## Signature

```d
double cond(M)(M m, Norm kind = Norm.spectral);
```

## Parameters

- `m`: a `GretlMatrix` or a `Submatrix` (a view of part of a matrix). It is not changed.
- `kind`: the norm: `Norm.spectral` (the default), `Norm.one`, `Norm.infinity` or `Norm.frobenius`. `Norm.max` does not define a condition number and throws a `GretlException`.

## Return value

A `double`, the condition number of `m`, `norm(m) * norm(inv(m))` in the chosen norm. The default is the ratio of the largest to the smallest singular value, which also works for a non-square matrix (as in Matlab). The other norms need a square matrix. A singular matrix has condition number `double.infinity` and an empty matrix 0. A value near 1 / 2.2e-16 means that `m` is singular to machine precision. A large condition number means that solving with `m` loses accuracy. It can also be called as a method: `m.cond`.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    auto good = eye(3);
    writeln(cond(good));                // 1

    auto bad = cbind([1.0, 1.0], [1.0, 1.0 + 1e-10]);   // nearly singular
    writefln("%.3e", cond(bad));

    writefln("%.3f", cond(bad, Norm.one));
}
```

Output:

```
1
4.000e+10
39999996694.385
```
