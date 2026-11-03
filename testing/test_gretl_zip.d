/* Smoke tests for gretl_zip.d: make a zip from a directory tree, unzip it elsewhere, zip/unzip
   a data.xml + data.bin pair, read the top directory of a session-like archive, reject unsafe
   entry names, and build a function-package zipfile with package_make_zipfile().

   Build and run with testing/run_tests.sh. */

import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;
import core.stdc.stdlib : free;
import core.stdc.string : strcmp, strstr;

import std.file : chdir, exists, getcwd, mkdirRecurse, read, rmdirRecurse, write;
import std.string : toStringz;
import std.zip : ArchiveMember, ZipArchive;

import gretl.gretl_paths : gretl_set_path_by_name;
import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_get_buffer, gretl_print_new;
import gretl.gretl_prnh : PrnType;
import gretl.gretl_zip;
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

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    string home = getcwd();
    string base = home ~ "/build/zip_test";

    if (exists(base)) {
        rmdirRecurse(base);
    }
    mkdirRecurse(base ~ "/src/sub");
    write(base ~ "/src/a.txt", "alpha\n");
    write(base ~ "/src/sub/b.txt", "beta\n");
    mkdirRecurse(base ~ "/src/empty");
    chdir(base);

    /* a tree: entries are named src/... */
    check(gretl_make_zipfile("tree.zip", "src") == 0 && exists("tree.zip"), "gretl_make_zipfile");
    mkdirRecurse("out1");
    check(gretl_unzip_into("tree.zip", "out1") == 0, "gretl_unzip_into");
    check(exists("out1/src/a.txt") && cast(string) read("out1/src/sub/b.txt") == "beta\n" && exists("out1/src/empty"),
          "tree restored (files, subdirectory, empty directory)");

    /* gretl_unzip: current directory */
    mkdirRecurse("out2");
    chdir("out2");
    check(gretl_unzip("../tree.zip") == 0 && exists("src/sub/b.txt"), "gretl_unzip (current directory)");
    chdir(base);

    /* data.xml + data.bin */
    mkdirRecurse("dat");
    write("dat/data.xml", "<gretldata/>\n");
    write("dat/data.bin", "BINARY");
    check(gretl_zip_datafile("data.zip", "dat", 6) == 0, "gretl_zip_datafile");
    mkdirRecurse("out3");
    check(gretl_unzip_into("data.zip", "out3") == 0 && cast(string) read("out3/data.bin") == "BINARY"
          && exists("out3/data.xml"), "data files round trip");
    write("dat2.tmp", "");
    mkdirRecurse("dat_missing");
    write("dat_missing/data.xml", "x");
    check(gretl_zip_datafile("bad.zip", "dat_missing", 6) != 0, "gretl_zip_datafile needs both files");

    /* session-style archive: top directory from session.xml */
    {
        auto zip = new ZipArchive();
        auto m = new ArchiveMember();

        m.name = "mysession/session.xml";
        m.expandedData = cast(ubyte[]) "<s/>".dup;
        zip.addMember(m);
        write("sess.zip", zip.build());
    }
    char * top = null;

    mkdirRecurse("out4");
    chdir("out4");
    check(gretl_unzip_session_file("../sess.zip", &top) == 0 && top !is null && strcmp(top, "mysession") == 0
          && exists("mysession/session.xml"), "gretl_unzip_session_file reports the top directory");
    free(top);
    chdir(base);

    /* zip-slip */
    {
        auto zip = new ZipArchive();
        auto m = new ArchiveMember();

        m.name = "../evil.txt";
        m.expandedData = cast(ubyte[]) "x".dup;
        zip.addMember(m);
        write("evil.zip", zip.build());
    }
    mkdirRecurse("out5");
    check(gretl_unzip_into("evil.zip", "out5") != 0 && !exists("evil.txt"), "unsafe entry name is rejected");
    check(gretl_unzip("no_such_file.zip") != 0, "missing zipfile is an error");

    /* package_make_zipfile */
    mkdirRecurse("dot");
    gretl_set_path_by_name("dotdir", toStringz(base ~ "/dot/"));
    mkdirRecurse("pkg");
    write("pkg/mypkg.gfn", "<gretl-function-package/>\n");
    Prn * prn = gretl_print_new(PrnType.buffer, null);
    char * zipname = null;

    int err = package_make_zipfile("pkg/mypkg.gfn", 0, null, 0, &zipname, null, cast(Gretlopt) 0, prn);

    check(err == 0 && zipname !is null && exists(base ~ "/dot/mypkg.zip"), "package_make_zipfile (zip left in the dotdir)");
    chdir(base);
    mkdirRecurse("out6");
    check(zipname !is null && gretl_unzip_into(zipname, "out6") == 0 && exists("out6/mypkg/mypkg.gfn"),
          "the package zip holds mypkg/mypkg.gfn");
    free(zipname);
    check(package_make_zipfile("pkg/mypkg.txt", 0, null, 0, null, "x.zip", cast(Gretlopt) 0, prn) != 0, "wrong extension is an error");
    check(strstr(gretl_print_get_buffer(prn), "OK") !is null, "progress messages were printed");
    gretl_print_destroy(prn);

    chdir(home);
    rmdirRecurse(base);
    printf("\n%d checks, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
