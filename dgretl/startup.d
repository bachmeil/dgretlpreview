/* Start-up and shut-down of the library. The module constructor runs libgretl_init() before main and sets up the
   things a program has no way to see; the module destructor runs libgretl_cleanup() after main and removes the scratch
   directory. Linking this module into a program is all it takes: nothing has to be called. */

module dgretl.startup;

import std.conv: to;
import std.file: exists, mkdirRecurse, rmdirRecurse, tempDir;
import std.process: thisProcessID;
import std.string: toStringz;

import gretl.gretl_paths: gretl_set_path_by_name, set_gretl_png_font;
import gretl.gretl_utils: libgretl_cleanup, libgretl_init;

private string dotDirectory;

shared static this() {
    libgretl_init();

    /* libgretl keeps scratch files (unpacked workbooks, plot scripts) in a "dot directory": one per process, under the
       temporary directory, removed again at exit */
    dotDirectory = tempDir() ~ "/dgretl-" ~ to!string(thisProcessID) ~ "/";
    mkdirRecurse(dotDirectory);
    gretl_set_path_by_name("dotdir", dotDirectory.toStringz);

    /* defaults a program has no way to see: the gnuplot program (found through PATH) and a PNG font */
    gretl_set_path_by_name("gnuplot", "gnuplot");
    set_gretl_png_font("sans 10");
}

shared static ~this() {
    libgretl_cleanup();

    try {
        if (exists(dotDirectory)) {
            rmdirRecurse(dotDirectory);
        }
    } catch (Exception) {
        /* leave the scratch directory behind rather than fail at exit */
    }
}
