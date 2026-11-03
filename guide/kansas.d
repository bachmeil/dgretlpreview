/* kansas.d: load the BEA state GDP growth spreadsheet and plot Kansas over 2005Q2-2026Q2. */

import dgretl;                                   // starts the library by itself and cleans up at exit

import std.stdio : writefln, writeln;
import core.stdc.string : strcmp;

import gretl.dataset : dataset_set_time_series, series_get_string_for_obs;
import gretl.graphing : gnuplot;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.libgretlh : Gretlopt;
import gretl.options : set_optval_string;

enum nquarters = 85;          // 2005Q2 ... 2026Q2

int main()
{
    // 1. read the worksheet: skip the 5 title rows so that the header row is the variable names. The series are
    //    GeoFIPS (number 0), GeoName (1) and then the 85 quarters (2 .. 86)
    auto raw = importSpreadsheet("gdp-by-state.xlsx", rowOffset: 5);

    // 2. find the Kansas row (the strings of a text series are reached through the low-level function)
    long row = -1;
    foreach (t; 0 .. raw.nobs) {
        const(char) * name = series_get_string_for_obs(raw.ptr, 2, cast(int) t);   // 2: GeoName, counting the constant
        if (name !is null && strcmp(name, "Kansas") == 0) {
            row = t;
            break;
        }
    }
    if (row < 0) {
        writeln("Kansas not found");
        return 1;
    }

    // 3. build a quarterly time-series dataset: one series, 85 observations starting at 2005Q2; the Dataset
    //    frees itself
    auto ks = Dataset(2, nquarters);
    dataset_set_time_series(ks.ptr, 4, 2005, 2);
    foreach (q; 0 .. nquarters) {
        ks.ptr.Z[1][q] = raw[2 + q][row];
    }
    ks.ptr.varname[1][0 .. 7] = "Kansas\0";

    writefln("Kansas: %d quarters, first %.1f, last %.1f", nquarters, ks.ptr.Z[1][0], ks.ptr.Z[1][nquarters - 1]);

    // 4. plot to a PNG
    set_optval_string(GretlCmdIndex.gnuplot, Gretlopt.OPT_U, "kansas.png");
    int[2] list = [1, 1];
    int err = gnuplot(list.ptr, null, ks.ptr, cast(Gretlopt) (Gretlopt.OPT_T | Gretlopt.OPT_O));
    if (err) {
        writefln("plot failed: err=%d", err);
        return 1;
    }
    writeln("wrote kansas.png");
    return 0;
}
