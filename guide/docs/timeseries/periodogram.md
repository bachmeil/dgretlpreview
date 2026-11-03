# periodogram

## Signature

```d
Periodogram periodogram(X)(X x, int window = -1);
Periodogram periodogram(Dataset d, string name, int window = -1);

struct Periodogram {
    GretlMatrix frequency;
    GretlMatrix period;
    GretlMatrix density;
    int window;
}
```

## Parameters

- `x`: the series: an array of numbers, or a `GretlMatrix` or `Submatrix` with one column. Missing values (NaN) are not allowed and throw a `GretlException`. At least 12 observations are needed, and the series must not be constant.
- `window`: -1 (the default) gives the plain sample periodogram; 0 smooths it with the Bartlett lag window of the library's automatic width (`2 sqrt(n)`); a positive value smooths with that width (at most `n / 2`).

## Return value

A `Periodogram` struct with three `n / 2` x 1 matrices at the Fourier frequencies `2 pi j / n`, `j = 1 .. n / 2`: `frequency` (in radians, between 0 and pi), `period` (`2 pi / frequency`, in observations) and `density` (the spectral density estimate: the periodogram is `|DFT|^2 / (2 pi n)` in the library's scaling). `window` is the window width used (-1 for none). A peak of the density marks a cycle with that period.

## Example

```d
import std.stdio;
import std.math: sin, PI;
import dgretl;

void main() {
    setSeed(20);
    auto noise = rnorm(120);
    auto x = new double[120];
    foreach (t; 0 .. 120) {
        x[t] = 2.0 * sin(2 * PI * t / 10.0) + noise[t];     // a cycle of period 10
    }
    auto p = periodogram(x);
    size_t best = 0;
    foreach (j; 0 .. p.density.nrow) {
        if (p.density[j, 0] > p.density[best, 0]) best = j;
    }
    writefln("%d frequencies; the peak is at period %.2f", p.density.nrow, p.period[best, 0]);

    auto s = periodogram(x, 0);                      // smoothed, with the automatic window
    writefln("window %d", s.window);
}
```

Output:

```
60 frequencies; the peak is at period 10.00
window 21
```
