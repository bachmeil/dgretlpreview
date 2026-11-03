/* Smoke tests for plugins.d: a tiny C shared object named like a gretl plugin ("spss_import.so",
   with no libgretl dependencies) is built into a scratch directory, found through
   GRETL_PLUGIN_PATH, and called through get_plugin_function(); unknown functions, a missing
   plugin file and the GMP plugin are checked to fail cleanly.

   Build and run with testing/run_tests.sh. */

import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;
version (Posix) import core.sys.posix.stdlib : setenv; else import gretl.winport : setenv;

import std.file : exists, mkdirRecurse, rmdirRecurse, write;
import std.process : execute;

import gretl.plugins;

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

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    string dir = "build/plugtest";

    if (exists(dir)) {
        rmdirRecurse(dir);
    }
    mkdirRecurse(dir);
    write(dir ~ "/sav.c", "int sav_get_data(int a, int b) { return a * 10 + b; }\n");
    auto r = execute(["gcc", "-shared", "-fPIC", "-o", dir ~ "/spss_import.so", dir ~ "/sav.c"]);

    check(r.status == 0, "built the test plugin");

    /* gretl_plugin_path() takes GRETL_PLUGIN_PATH only on first use: set it before anything else */
    setenv("GRETL_PLUGIN_PATH", "build/plugtest", 1);

    alias Fn = extern (C) int function(int, int);
    Fn f = cast(Fn) get_plugin_function("sav_get_data");

    check(f !is null, "get_plugin_function finds a function in the plugin");
    check(f !is null && f(4, 2) == 42, "the plugin function runs");

    check(get_plugin_function("no_such_function") is null, "unknown function name -> NULL");
    check(get_plugin_function("wf1_get_data") is null, "function of a plugin file that does not exist -> NULL");
    check(get_plugin_function("mplsq") is null, "GMP plugin is refused");

    void * h = gretl_dlopen("build/plugtest/spss_import.so", 1);

    check(h !is null && gretl_dlsym(h, "sav_get_data") !is null && gretl_dlsym(h, "nope") is null, "gretl_dlopen/gretl_dlsym");
    close_plugin(h);

    plugins_cleanup();
    check(true, "plugins_cleanup");

    rmdirRecurse(dir);
    printf("\n%d checks, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
