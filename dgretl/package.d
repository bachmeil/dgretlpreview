/* dgretl: the D-friendly front end of the library.

   `import dgretl;` is all a program needs. The library is started before main() runs and shut down after it
   returns (dgretl.startup: no libgretl_init / libgretl_cleanup calls), data live in reference-counted structs that
   free themselves (dgretl.dataset, dgretl.matrix), data files are read by the functions of dgretl.importers, and
   failures are GretlException (dgretl.exception):

       import dgretl;
       import std.stdio: writeln;

       void main() {
           auto d = importCsv("example.csv");
           writeln(d.nobs, " observations: ", d.names);
           writeln("gdp: ", d["gdp"][0 .. 2]);

           auto m = importCsvAsMatrix("example.csv");
           writeln(m.nrow, " x ", m.ncol, ", m[2, 1] = ", m[2, 1]);
       }

   The modules in the gretl package are the direct ports of the files of libgretl, one for one; the ones here are
   written for D programmers. The wrapper structs wrap the C-layout ones (DatasetRaw, GretlMatrixRaw): the `ptr`
   property of each gives the raw pointer for the functions of the ports. */

module dgretl;

public import dgretl.dataset;
public import dgretl.descriptive;
public import dgretl.distributions;
public import dgretl.estimate;
public import dgretl.exception;
public import dgretl.importers;
public import dgretl.random;
public import dgretl.matrix;
public import dgretl.startup;
public import dgretl.timeseries;
