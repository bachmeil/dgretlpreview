# Installing, building and running a first program

dgretl is a port of libgretl to D. There is no separate library to install: the port's source files are compiled
together with your program, which takes a few seconds. This page lists what you need, shows how to build and run a
program, and walks through a first one.

## What you need

- A D compiler: `ldmd2` (the LDC compiler with `dmd`'s command line) or `dmd` itself. LDC 1.39 or later is needed: older versions, such as the 1.36 that
  Ubuntu 24.04 packages, cannot compile the named arguments used by the front end's templates. LDC 1.43 is the version used for the latest testing.
- The development packages of the C libraries that the port links with: libxml2, zlib, LAPACK, BLAS and FFTW3 (OpenBLAS can replace LAPACK and BLAS, see below; Windows is covered further down). On
  Debian or Ubuntu they are called `libxml2-dev`, `zlib1g-dev`, `liblapack-dev`, `libblas-dev` and `libfftw3-dev`:

  ```bash
  sudo apt install libxml2-dev zlib1g-dev liblapack-dev libblas-dev libfftw3-dev
  ```

  This is the tested route on Linux. Windows is covered further down, and macOS below that (not tested).
- `gnuplot` on your `PATH`, only if you draw plots (the worked example [Loading a spreadsheet and plotting a
  series](loading-and-plotting.md) does).

## OpenBLAS instead of reference LAPACK and BLAS

OpenBLAS is a much faster implementation of BLAS, and it also contains LAPACK, so it can replace both `-llapack` and
`-lblas`. This matters for large matrix products, `solve`, `svd` and the like. Install `libopenblas-dev` (Debian,
Ubuntu) and replace `-L-llapack -L-lblas` with `-L-lopenblas` on the link line:

```bash
ldmd2 ... -L-lm -L-lxml2 -L-lz -L-lopenblas -L-lfftw3
```

A program built this way ran correctly with that change (checked on Linux with the Debian OpenBLAS packages).
With dub, change `"lapack" "blas"` to `"openblas"` in the `libs` line of `dub.sdl`. On Debian and Ubuntu you can
instead keep the build as is and switch the system's BLAS and LAPACK to OpenBLAS with
`sudo update-alternatives --config libblas.so.3-x86_64-linux-gnu` (and the same for `liblapack.so.3-x86_64-linux-gnu`);
programs then pick OpenBLAS up at run time without rebuilding.

## Windows

The port builds and runs natively on Windows 11 with LDC 1.41 (`x86_64-pc-windows-msvc`) and the Visual Studio 2022
Build Tools (for `link.exe`). `testing/build_win.ps1` is the working recipe; WSL2 (Ubuntu) also works and needs no changes.

**Dependencies** (vcpkg, triplet `x64-windows`, DLL builds):

```
vcpkg install libxml2 zlib openblas clapack fftw3 libiconv --triplet x64-windows
```

- The vcpkg OpenBLAS has BLAS but no LAPACK, so `clapack` (`lapack.lib`, `libf2c.lib`) supplies LAPACK.
- The `openblas.lib` import library also defines `printf`, `sprintf` and other CRT functions, which hijack the C runtime's
  and crash programs at start-up. `build_win.ps1` builds `openblas_clean.lib` from `dumpbin /exports openblas.dll`, keeping only
  the BLAS entry points; link that instead.
- Link the dynamic C runtime (`-mscrtlib=msvcrt` and `-L/NODEFAULTLIB:libcmt`), to match the vcpkg DLLs. Phobos already
  contains zlib, so `zlib.lib` is not linked. Link `libxml2.lib`, `fftw3.lib`, `iconv.lib`, `charset.lib`.
- Put `vcpkg\installed\x64-windows\bin` on `PATH` only when you run a program, **not while compiling**: with it on `PATH`, `ldc2` itself picks up vcpkg's `z.dll`/`iconv`/`libxml2` instead of its own copies and crashes (segmentation fault or "illegal instruction"). `build_win.ps1` adds it only around the test runs.

**Building and running the tests:** from PowerShell 7 (`pwsh`), in the repository root: `testing\build_win.ps1 [test_name ...]`
(set `VCPKG_ROOT` if vcpkg is not in `%USERPROFILE%\repos\vcpkg`; `-NoRun` only builds). It compiles every source file to
its own object file because compiling many modules in one run crashed LDC 1.41 (some of it was the `PATH` problem above; I have not retried the one-run build since), and the `dgretl` modules
together (without `-g`). Use absolute paths when linking by hand: LDC 1.41 crashes on some relative ones.
LDC 1.41 also crashes now and then ("illegal instruction") on a command that works the next time, so the script repeats a failed compile up to three times and gives each one 5 minutes.

**Source changes made for Windows:** `gretl/winport.d` supplies the POSIX functions that the port uses (`stat`, `opendir`,
`mkstemp`, `strptime`, `timegm`, `dlopen`, `mmap`, `waitpid`, ...) on top of the C runtime and Win32; `gretl/gspawn_win.d`
runs child processes (gnuplot) with `CreateProcess`. The existing Linux code is unchanged (`version (Posix)`).
Differences you may notice: C `long` is 32 bits (`c_long` is used where the C code passes a `long` to the C library), gzipped
XML is decompressed by the port because the vcpkg libxml2 does not do it, and `std.conv.to!double` can be 1 ulp off
(use `strtod`). `gnuplot` must be installed for plots (the port starts it with the path from the `gnuplot` line of the configuration file, or `gnuplot` from `PATH`; the full path may contain spaces). LDC 1.43 is steadier than 1.41 on Windows and is the better choice. `dub build` works on Windows with LDC 1.41: `dub.sdl` has the Windows library names and flags (the vcpkg `lib` directory must be on the `LIB` environment variable, and the vcpkg `bin` directory on `PATH` only when you run the program, not while dub compiles). It passes `-lowmem` to LDC, because LDC 1.41 on Windows otherwise fails at random (bogus template errors, crashes) when compiling this many modules in one run. `testing/build_win.ps1` is the other way to build.
## macOS (not tested)

Nobody has built dgretl on a Mac yet. The code is POSIX, so it is expected to need only small changes, and reports and pull
requests are welcome. What should be needed:

- A D compiler and the libraries from Homebrew: `brew install ldc libxml2 zlib fftw openblas gnuplot`. Use LDC 1.39 or later.
- BLAS and LAPACK: either OpenBLAS (see above; `-L-lopenblas`, and `brew --prefix openblas` for its location) or Apple's Accelerate
  framework (`-L-framework -LAccelerate` instead of `-llapack -lblas`).
- Tell the linker where Homebrew keeps the libraries: `/opt/homebrew/lib` on Apple Silicon, `/usr/local/lib` on Intel (for
  example `-L-L/opt/homebrew/lib` with `ldmd2`, or `export LIBRARY_PATH=/opt/homebrew/lib`). `dub.sdl` lists the libraries under
  `platform="posix"`, which includes macOS.
- Places likely to need a fix: calls in `core.sys.posix` whose layout differs on Darwin (`stat`, `dirent`), `mallinfo`, which some tests
  use to look for leaks and which macOS does not have (`testing/stubs.d` has a Windows replacement and would need one for macOS),
  `nl_langinfo`, and the way `gnuplot` is found.

## The source tree

The repository has these parts:

| Directory | What it holds |
|---|---|
| `gretl/` | the direct ports of the libgretl files, one D module for each, with the library's own names; `gretl/plugin/` and `gretl/minpack/` hold the plugins and the Minpack routines |
| `cephes/` | the special-function and probability library libgretl uses |
| `dgretl/` | the D-friendly front end: what `import dgretl;` brings in |
| `testing/` | the test programs, and two files of stand-ins for the parts of the library that are not ported yet (`stubs.d`, `stubs_unported.d`) |
| `guide/` | this guide, its reference pages (`guide/docs/`) and a build script |

A program uses `dgretl`. The modules in `gretl/` are there for the cases where the front end is not enough.

## Using dgretl from a dub project

The port has a `dub.sdl`, so a dub project can depend on it. With this repository cloned at `/path/to/dgretlpreview`, put this in
your project's `dub.sdl`:

```sdl
dependency "dgretl" path="/path/to/dgretlpreview"
```

(or `dependency "dgretl" repository="git+https://github.com/bachmeil/dgretlpreview.git" version="a6133f8cb9cff8864db294ea87b14ab45155fb78"`, where the version is a commit hash, here the one tagged `v0.1.0`) and write
`import dgretl;` in your code. dub also builds the `dgretl:stubs` sub-package, the aborting stand-ins for the unported
functions, and links libxml2, zlib, LAPACK, BLAS and FFTW3. The build works with `dub build --compiler=ldc2` and with
dmd; `-b release` makes a faster library.

## Building a program

The script `guide/build.sh` builds a program that is in the `guide/` directory. For a file `guide/first.d` (the
program below):

```bash
guide/build.sh first
```

makes the executable `guide/first`. The script does two things, both with the D compiler, and you can do them by hand to build a program
that lives anywhere. Set `PORT` to the directory of this repository:

```bash
PORT=/path/to/dgretlpreview
mkdir -p build

# 1. the stand-ins for the libgretl functions that are not ported yet, compiled on their own
ldmd2 -c -of=build/stubs_unported.o $PORT/testing/stubs_unported.d

# 2. your program together with the whole port
ldmd2 -I$PORT -of=first first.d \
      $PORT/gretl/*.d $PORT/gretl/plugin/*.d $PORT/gretl/minpack/*.d $PORT/cephes/*.d \
      $(ls $PORT/dgretl/*.d | grep -v -e /testgslprng.d -e /slowPrng.d) \
      $PORT/testing/stubs.d build/stubs_unported.o \
      -L-lm -L-lxml2 -L-lz -L-llapack -L-lblas -L-lfftw3
```

- `-I$PORT` lets the compiler find `import dgretl;`.
- The `ls ... | grep -v` leaves out two files of `dgretl/` that are stand-alone programs with their own `main`
  (`testgslprng.d` and `slowPrng.d`).
- The stand-ins: parts of libgretl are not ported yet (see "Not in version 1" in the [index](index.md)).
  `testing/stubs_unported.d` gives every function the ported code could call in them a stand-in that prints
  `*** call to unported libgretl function NAME()` and aborts, so a missing piece fails loudly instead of returning
  garbage. `testing/stubs.d` holds working stand-ins for the pieces the library calls in normal use, for example
  "no user functions are defined" or "single-threaded". Every program that uses the library links both. They are plain D, so no C compiler is needed. The stand-ins are compiled on their own (step 1) because the ported modules declare these functions with their real signatures and the stubs do not, which the compiler would otherwise reject.
- With `dmd` instead of `ldmd2` the command is the same.

## A first program

Save this as `first.d`:

```d
import std.stdio;
import dgretl;

void main() {
    setSeed(1);
    auto n = 100;
    auto x = rnorm(n);                         // 100 standard normal draws
    auto e = rnorm(n);
    auto y = new double[n];
    foreach (t; 0 .. n) {
        y[t] = 1.0 + 2.0 * x[t] + e[t];        // y = 1 + 2 x + noise
    }

    writefln("mean of y %.3f, standard deviation %.3f", mean(y), stddev(y));
    auto fit = ols(y, x);                      // least squares of y on x and a constant
    writeln(fit);
}
```

Build it as above and run `./first`. It prints:

```
mean of y 1.138, standard deviation 2.214
OLS, n = 100, df = 98, standard errors: classical
                     coef      std.err    t-ratio   p-value
const             1.00837    0.0909233    11.0903    0.0000
x1                2.07409     0.093547    22.1717    0.0000
R-squared 0.833781, adjusted 0.832085, sigma 0.907357, F 491.6
log-likelihood -131.161771, AIC 266.3235, BIC 271.5339
```

What the program does:

- `import dgretl;` is the only import a program needs for the library. There is nothing to start or shut down:
  `dgretl` starts libgretl before `main` runs (settings, the random number generator, the XML library, a
  per-process scratch directory, the gnuplot program and a default PNG font for plots) and shuts it down, removing
  the scratch directory, after `main` returns. A program that does not import `dgretl` has to call
  `libgretl_init()` and `libgretl_cleanup()` itself.
- `setSeed(1)` makes the random numbers, and so the output, the same on every run ([`setSeed`](docs/random/setSeed.md)).
- `rnorm(n)` returns a `double[]` ([`rnorm`](docs/random/rnorm.md)), and an array goes straight into `mean`,
  `stddev` ([`mean`](docs/descriptive/mean.md), [`stddev`](docs/descriptive/stddev.md)) and `ols`
  ([`ols`](docs/regression/ols.md)).
- `ols` returns an `Estimate`, a plain struct, and `writeln` prints it as a table of coefficients. The pieces are
  fields: `fit.coef`, `fit.se`, `fit.rsq`, and so on ([`Estimate`](docs/regression/Estimate.md)).
- Nothing is freed by hand.

## When something goes wrong

- A problem with an argument or a failure inside the library is a `GretlException` with a message that says what
  happened ([`GretlException`](docs/errors/GretlException.md)). An uncaught one ends the program with that message.
- `*** call to unported libgretl function NAME()` means the program reached a part of libgretl that has not been
  ported. It is not a bug in your program; see "Not in version 1" in the [index](index.md).
- A linking error that mentions `xml`, `lapack`, `blas` or `fftw3` means that the development package for that
  library is missing.

## Where to go next

- [Conventions](conventions.md): how data go in and come out, missing values, names and options, used throughout.
- [Loading a spreadsheet and plotting a series](loading-and-plotting.md): a complete worked example with real data.
- The [index](index.md), which lists every function with a link to its reference page.
