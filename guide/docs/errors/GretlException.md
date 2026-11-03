# GretlException

## Signature

```d
class GretlException : Exception {
    this(string msg, string file = __FILE__, size_t line = __LINE__);
}
```

## Parameters

- `msg`: the message, which has the details of the failure. For an error reported by the library it includes the library's own error text; the import functions add the importer's progress messages as well.
- `file`, `line`: where the exception was thrown, as for any `Exception` (the defaults are the place of the `throw`).

## Return value

`GretlException` is the one exception type that every failure in the dgretl modules is reported with, so a program can catch all of them with one `catch (GretlException e)`. It is an ordinary D `Exception`, so `e.msg`, `e.file`, `e.line` and `e.info` work as usual, and it is also caught by `catch (Exception e)`.

What throws it:

- invalid arguments: a negative standard deviation, a probability outside [0, 1], a wrong number of rows, arrays of different lengths, a lag order out of range, missing values where none are allowed;
- a failure inside the library: a singular matrix, a matrix that is not positive definite, a regression that cannot be estimated, a test that cannot be computed;
- an index outside a matrix, a submatrix or a dataset, and a series name that does not exist;
- a file that cannot be read.

The message usually names the function that failed (for example `rnorm: need a finite mean and a positive sd`). Nothing is left half done when it is thrown: the matrices, datasets and estimates already made are freed by the garbage collector or by their reference counts, and the library settings that a function changed (for a robust covariance matrix, a quantile rule and the like) are restored. A program can create and throw its own `GretlException` too, which is handy in code that wraps these functions.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    auto good = rnorm(3, 0, 1);

    try {
        rnorm(3, 0, -1);                              // a negative standard deviation
    } catch (GretlException e) {
        writeln("1: ", e.msg);
    }

    try {
        auto m = inv(cbind([1.0, 2.0], [2.0, 4.0]));  // a singular matrix
    } catch (GretlException e) {
        writeln("2: ", e.msg);
    }

    try {
        auto x = cbind([1.0, 2.0, 3.0])[5, 0];        // an index outside the matrix
    } catch (GretlException e) {
        writeln("3: ", e.msg);
    }

    try {
        pnorm([1.0, 2.0, 3.0], sd: [1.0, 2.0]);       // arrays of different lengths
    } catch (Exception e) {                            // GretlException is an Exception
        writeln("4: ", e.msg);
    }

    try {
        throw new GretlException("something of my own");
    } catch (GretlException e) {
        writeln("5: ", e.msg);
    }
}
```

Output:

```
1: rnorm: need a finite mean and a positive sd
2: inv: the matrix could not be inverted (singular?), error code 3
3: GretlMatrix: index [5, 0] is out of range
4: the arguments must have the same number of elements (3 and 2)
5: something of my own
```
