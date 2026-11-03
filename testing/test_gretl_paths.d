/* Smoke tests for gretl_paths.d: path building and normalizing, file
   wrappers, directory creation/deletion, path searching, reading the
   configuration file and setting up the standard paths (in a scratch
   HOME directory).

   Build and run with testing/run_tests.sh. */

import core.stdc.stdio : FILE, _IOLBF, fclose, fopen, fprintf, fputs, printf, setvbuf, stdout;
import core.stdc.stdlib : free;
import std.file : exists, rmdirRecurse;
import std.string : fromStringz;
version (Posix) import core.sys.posix.stdlib : setenv; else import gretl.winport : setenv;
import core.stdc.string : strcmp, strlen, strstr;
version (Posix) import core.sys.posix.sys.stat : stat_t; else import gretl.winport : stat_t;

import gretl.gretl_errorsh : GretlError;
import gretl.gretl_paths;
import gretl.gfileutils : g_get_current_dir;
import gretl.gretl_pathsh : ConfigPaths, PkgType, SearchType;
import gretl.libgretlh;

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

void touch(const(char) * name, string content = "x\n")
{
    FILE * fp = fopen(name, "w");

    fputs((content ~ "\0").ptr, fp);
    fclose(fp);
}

void test_build_normalize()
{
    char[FILENAME_MAX] buf = 0;

    gretl_build_path(buf.ptr, "/usr/".ptr, "/share/".ptr, "gretl".ptr, null);
    check(strcmp(buf.ptr, "/usr/share/gretl") == 0, "gretl_build_path joins and collapses slashes");
    gretl_build_path(buf.ptr, "a".ptr, "".ptr, "b.txt".ptr, null);
    check(strcmp(buf.ptr, "a/b.txt") == 0, "gretl_build_path skips empty elements");

    strcpy_(buf, "/a/b/../c/./d");
    check(gretl_normalize_path(buf.ptr) == 0 && strcmp(buf.ptr, "/a/c/d") == 0, "gretl_normalize_path removes . and ..");

    strcpy_(buf, "/foo/bar");
    check(slash_terminate(buf.ptr) == 1 && strcmp(buf.ptr, "/foo/bar/") == 0 && slash_terminate(buf.ptr) == 0,
          "slash_terminate");

    strcpy_(buf, "data.gdt");
    check(gretl_path_prepend(buf.ptr, "/tmp/dir") == 0 && strcmp(buf.ptr, "/tmp/dir/data.gdt") == 0,
          "gretl_path_prepend");

    char[16] small = 0;

    check(gretl_path_compose(small.ptr, 16, "abcdefgh", "12345678") == GretlError.data, "gretl_path_compose rejects overlong");
    check(gretl_path_compose(small.ptr, 16, "abc", "def") == 0 && strcmp(small.ptr, "abcdef") == 0, "gretl_path_compose");
    check(utf8_encoded("plain") == 0 && utf8_encoded("café") == 1, "utf8_encoded");
}

void strcpy_(char[] dest, string s)
{
    dest[0 .. s.length] = s[];
    dest[s.length] = 0;
}

import core.stdc.stdio : remove;
alias remove_ = remove;

void test_files()
{
    /* absolute paths: gretl_deltree() chdir()s into the directory and back out, so it only works
       for a single relative component (see bugs.md #59) */
    char * cwd = g_get_current_dir();
    string base = cast(string) fromStringz(cwd).idup ~ "/build/paths_test";
    string[] names = [base, base ~ "/a/one.txt", base ~ "/a/two.txt", base ~ "/a/b", base ~ "/nope",
                      base ~ "/t.XXXXXX", base ~ "/xml.dat", base ~ "/a"];
    foreach (ref n; names) n ~= "\0";
    const(char) * dir = names[0].ptr;
    const(char) * f1 = names[1].ptr;
    const(char) * f2 = names[2].ptr;
    free(cwd);

    gretl_deltree(dir);
    check(gretl_mkdir(names[3].ptr) == 0 && gretl_isdir(names[3].ptr), "gretl_mkdir (with parents) and gretl_isdir");

    touch(f1);
    check(gretl_file_exists(f1) && !gretl_file_exists(names[4].ptr), "gretl_file_exists");
    check(gretl_test_fopen(f1, "r") == 0 && gretl_test_fopen(names[4].ptr, "r") != 0, "gretl_test_fopen");
    check(gretl_rename(f1, f2) == 0 && gretl_file_exists(f2) && !gretl_file_exists(f1), "gretl_rename");

    FILE * fp = gretl_fopen(f2, "r");

    check(fp !is null, "gretl_fopen");
    if (fp) fclose(fp);

    stat_t sb;

    check(gretl_stat(f2, &sb) == 0 && sb.st_size == 2, "gretl_stat");

    char[256] tmpl = 0;
    strcpy_(tmpl, base ~ "/t.XXXXXX");
    fp = gretl_mktemp(tmpl.ptr, "wb");
    check(fp !is null && strstr(tmpl.ptr, "XXXXXX") is null, "gretl_mktemp");
    if (fp) fclose(fp);

    touch(names[6].ptr, "<?xml version=\"1.0\"?>\n<a/>\n");
    check(gretl_is_xml_file(names[6].ptr) == 1 && gretl_is_xml_file(f2) == 0, "gretl_is_xml_file");
    check(gretl_write_access(cast(char *) names[7].ptr) == 0, "gretl_write_access");

    check(gretl_deltree(dir) == 0 && !gretl_file_exists(f2) && !gretl_isdir(dir), "gretl_deltree");
}

void test_search_and_config()
{
    /* a scratch HOME: gretl_set_paths() creates ~/.gretl, reads GRETL_HOME for gretldir */
    if (exists("build/paths_home")) rmdirRecurse("build/paths_home");
    gretl_mkdir("build/paths_home/gretl_share/data/mydata");
    gretl_mkdir("build/paths_home/work");
    touch("build/paths_home/gretl_share/COPYING");
    touch("build/paths_home/gretl_share/data/mydata/sample.gdt");
    touch("build/paths_local.csv");

    setenv("HOME", "build/paths_home", 1);
    setenv("GRETL_HOME", "build/paths_home/gretl_share", 1);

    /* a config file */
    touch("build/paths_home/rc", "gretldir = build/paths_home/gretl_share\nworkdir = build/paths_home/work\n"
          ~ "gnuplot = /opt/gp/bin/gnuplot\nshellok = true\n");
    setenv("GRETL_CONFIG_FILE", "build/paths_home/rc", 1);

    ConfigPaths cp;
    char[64] proxy = 0;
    int use_proxy = 0, updated = 0;
    char * theme = null;
    FILE * fp = fopen("build/paths_home/rc", "r");

    foreach (ref c; cp.gretldir) c = 0;
    foreach (ref c; cp.workdir) c = 0;
    foreach (ref c; cp.gnuplot) c = 0;
    get_gretl_config_from_file(fp, &cp, proxy.ptr, &use_proxy, &updated, &theme);
    fclose(fp);
    check(strcmp(cp.gretldir.ptr, "build/paths_home/gretl_share") == 0 &&
          strcmp(cp.gnuplot.ptr, "/opt/gp/bin/gnuplot") == 0, "get_gretl_config_from_file: gretldir, gnuplot");

    char[FILENAME_MAX] rcpath = 0;

    get_gretl_rc_path(rcpath.ptr);
    check(strcmp(rcpath.ptr, "build/paths_home/rc") == 0, "get_gretl_rc_path honours GRETL_CONFIG_FILE");

    check(cli_read_rc() == 0, "cli_read_rc (sets up the standard paths)");
    check(strstr(gretl_home(), "gretl_share/") !is null, "gretl_home");
    check(strstr(gretl_dotdir(), "paths_home/.gretl/") !is null && gretl_isdir(gretl_dotdir()), "gretl_dotdir created");
    check(gretl_isdir("build/paths_home/.gretl/functions") != 0, "dotdir/functions created");
    check(strstr(gretl_gnuplot_path(), "/opt/gp/bin/gnuplot") !is null, "gretl_gnuplot_path from the config file");

    char[FILENAME_MAX] name = 0;

    strcpy_(name, "paths_local.csv");
    check(gretl_addpath(name.ptr, 0) !is null && strstr(name.ptr, "build/paths_local.csv") !is null, "gretl_addpath finds a file in a subdirectory of workdir");
    strcpy_(name, "sample");
    char * found = gretl_addpath(name.ptr, 0);

    check(found !is null && strstr(name.ptr, "data/mydata/sample.gdt") !is null, "gretl_addpath finds a data file in a data subdirectory (adds .gdt)");
    strcpy_(name, "no_such_file.xyz");
    check(gretl_addpath(name.ptr, 0) is null && strcmp(name.ptr, "no_such_file.xyz") == 0, "gretl_addpath: not found leaves the name unchanged");

    check(get_full_filename("paths_local.csv", name.ptr, cast(Gretlopt) 0) == 0 && strstr(name.ptr, "build/paths_local.csv") !is null, "get_full_filename");
    check(get_full_filename("nothing_here.csv", name.ptr, cast(Gretlopt) 0) == GretlError.fopen, "get_full_filename: error code for a missing file");

    int n = 0;
    char ** dirs = get_plausible_search_dirs(SearchType.data_search, &n);

    check(dirs !is null && n >= 3 && strstr(dirs[0], ".gretl/data") !is null, "get_plausible_search_dirs");

    strcpy_(name, "x.gdt");
    gretl_maybe_prepend_dir(name.ptr);
    version (Windows) {
        /* a drive-letter path with backslashes */
        check(strstr(name.ptr, "x.gdt") !is null && (name[0] == '/' || name[0] == '\\' || name[1] == ':'), "gretl_maybe_prepend_dir (workdir)");
    } else {
        check(strstr(name.ptr, "/x.gdt") !is null && name[0] == '/', "gretl_maybe_prepend_dir (workdir)");
    }
    strcpy_(name, "~/x.gdt");
    gretl_maybe_prepend_dir(name.ptr);
    check(strstr(name.ptr, "paths_home/x.gdt") !is null, "gretl_maybe_prepend_dir expands ~/");

    gretl_set_script_dir("/a/b/s1.inp");
    gretl_set_script_dir("/c/d/s2.inp");
    gretl_set_script_dir("/a/b/s3.inp");
    gretl_script_dirs_cleanup();
    check(true, "gretl_set_script_dir / cleanup run");

    rmdirRecurse("build/paths_home");
    remove_("build/paths_local.csv");
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);
    test_build_normalize();
    test_files();
    test_search_and_config();
    printf("\n%d checks, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
