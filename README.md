# dgretlpreview

**Disclosure:** Most of the work was done by Claude. It is a port with a more convenient interface to the low-level functionality. There will be changes in the future and more human involvement now that there is a foundation to build on.

This repo is a preview release of **dgretl**: [libgretl](https://gretl.sourceforge.net/) (the econometrics library behind gretl) ported to D, with an idiomatic D front end.

```d
import dgretl;
import std.stdio: writeln;

void main() {
    auto d = importCsv("example.csv");
    writeln(d.nobs, " observations: ", d.names);
}
```

One import, `import dgretl;`, gives you data sets, matrices, estimation (OLS and the rest of the regression family, VAR and VECM, ARMA, GARCH),
hypothesis tests, descriptive statistics, distributions and random numbers. There is nothing to start up or shut down. Failures are
`GretlException`s.

## Status

| Platform | State |
|---|---|
| Linux | Tested: the whole test suite passes (LDC 1.39 and 1.43, dmd) |
| Windows | Tested: the whole test suite passes (LDC 1.41 and 1.43, MSVC target, libraries from vcpkg) |
| macOS | **Not tested.** Expected to need small fixes; see the installation guide |

This is a preview: the API may still change. Not ported: the X-12-ARIMA and TRAMO wrappers and the panel unit-root tests.
The port is a translation of gretl 2023c, and the original C is not part of this repository.

## Use it

With dub (the libraries libxml2, zlib, LAPACK, BLAS and FFTW3 must be installed; LDC 1.39 or later, or dmd):

```sdl
dependency "dgretl" repository="git+https://github.com/bachmeil/dgretlpreview.git" version="a6133f8cb9cff8864db294ea87b14ab45155fb78"
```

dub wants a commit hash here, not a tag name; the one above is the commit tagged `v0.1.0`.

Or clone the repository and use `dependency "dgretl" path="/path/to/dgretlpreview"`.

## Where to read

- [`guide/installing-and-building.md`](guide/installing-and-building.md): what to install and how to build, for Linux, Windows and macOS
- [`guide/index.md`](guide/index.md): the guide, with a reference page per function in `guide/docs/`
- [`overview.md`](overview.md): the main functions at a glance
- [`testing/README.md`](testing/README.md): the test suite (`bash testing/run_tests_fast.sh`)

Notes made while porting, including the suspected bugs in the upstream C code, are in the development repository,
<https://github.com/bachmeil/gretlport> (`bugs.md`, `log.md`). Comments in the source that mention them refer to those files.

## License

GPL, version 3 or later, like gretl. See `LICENSE`. This code is derived from gretl (Allin Cottrell and Riccardo "Jack" Lucchetti and
contributors), the Cephes mathematical library (Stephen L. Moshier), MINPACK and L-BFGS-B (see `gretl/minpack/README`).
