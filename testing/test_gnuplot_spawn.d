/* Runs gnuplot through the port's spawn code: gretl_spawn() on a quoted command line (the path may contain spaces, as under
   C:\Program Files), then gnuplot_make_image() with the path read from a configuration file, producing a PNG.
   Needs gnuplot: the test is skipped (and passes) when the GNUPLOT_EXE environment variable is not set; build_win.ps1 sets it
   when it finds gnuplot.

   Build and run with testing/run_tests.sh or testing\build_win.ps1. */

import core.stdc.stdio : FILE, fclose, fopen, fputs, printf, remove, setvbuf, stdout, _IOLBF;
import core.stdc.stdlib : getenv;
import core.stdc.string : strlen;
import std.file : exists, getSize;
import std.format : format;
import std.string : toStringz;
version (Posix) import core.sys.posix.stdlib : setenv; else import gretl.winport : setenv;

import gretl.gretl_paths : cli_read_rc, gretl_gnuplot_path;
import gretl.gretl_utils : gnuplot_make_image, gretl_spawn;
import gretl.gretl_errors : gretl_errmsg_get;

int n_tests, n_failed;

void check(bool ok, string what)
{
    n_tests++;
    if (!ok) {
        n_failed++;
        printf("FAIL: %.*s\n", cast(int) what.length, what.ptr);
    } else {
        printf("ok:   %.*s\n", cast(int) what.length, what.ptr);
    }
}

void write(string name, string content)
{
    FILE * fp = fopen(name.toStringz, "w");

    fputs(content.toStringz, fp);
    fclose(fp);
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);
    const(char) * exe = getenv("GNUPLOT_EXE");

    if (exe is null) {
        printf("GNUPLOT_EXE is not set: skipping the gnuplot tests\n");
        return 0;
    }
    string gp = exe[0 .. strlen(exe)].idup;

    /* a quoted command line with a full path */
    string cmd = "\"" ~ gp ~ "\" --version\0";
    char[] cl = cmd.dup;
    int err = gretl_spawn(cl.ptr);
    check(err == 0, "gretl_spawn runs 'gnuplot --version' by its full (quoted) path");

    /* a bad command is reported as an error */
    cl = ("\"" ~ gp ~ "-no-such\" --version\0").dup;
    check(gretl_spawn(cl.ptr) != 0, "gretl_spawn reports a missing program");

    /* gnuplot_make_image with the path from the configuration file */
    string png = "build/gp_test.png";
    if (exists(png)) remove(png.toStringz);
    write("build/gp_test.gp", "set terminal png size 200,150\nset output 'build/gp_test.png'\nplot sin(x)\n");
    write("build/gp_test.rc", "gnuplot = " ~ gp ~ "\n");
    setenv("GRETL_CONFIG_FILE", "build/gp_test.rc", 1);
    setenv("HOME", "build", 1);
    cli_read_rc();
    check(gretl_gnuplot_path() !is null && strlen(gretl_gnuplot_path()) > 0, "the gnuplot path is read from the configuration file");
    err = gnuplot_make_image("build/gp_test.gp");
    check(err == 0, "gnuplot_make_image returns no error");
    check(exists(png) && getSize(png) > 100, "gnuplot wrote a PNG file");
    if (exists(png)) {
        ubyte[8] sig;
        import std.stdio : File;
        File(png, "rb").rawRead(sig[]);
        check(sig == [0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A], "the file has a PNG signature");
    }

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
