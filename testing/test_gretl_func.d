/* Tests for gretl_func0 .. gretl_func4.d (gretl_func.c): user-defined functions and function packages.

   What can run without an interpreter: the definition of functions (parsing of the signature, parameters with defaults, ranges, labels,
   descriptions, const / null / auto qualifiers, return types, name-collision handling, deleting), the introspection API, the call
   structures (fncall_new, push_function_arg, fncall_destroy), the printing and retrieval of the code, indentation and normalization of
   hansl code, and function packages (properties, writing a .gfn file, reading its header back, peeking at its dependencies).
   Executing the body of a function needs the command interpreter (interact.c, flow_control.c, monte_carlo.c), which is not ported yet.

   The functions of unported files that definition-time code calls are given stand-ins here, in the simplest form that works:
   get_command_index (classifies a line by its first word), gretl_cmd_init / gretl_cmd_free, statements_get_structure (no recursion). Remove them as those files are ported.

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs;
import core.stdc.stdio : FILE, fclose, fopen, fread, printf, setvbuf, stdout, _IOLBF;
import core.stdc.stdlib : free, malloc;
import core.stdc.string : memset, strcmp, strlen, strncmp, strstr;

import std.file : exists, remove, tempDir;
import std.string : toStringz;

import dgretl;

import gretl.cmd_privateh : Cmd, CmdFlags, ExecState, Stmt;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_func;
import gretl.gretl_errors : gretl_error_clear, gretl_errmsg_get;
import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_get_buffer, gretl_print_new;
import gretl.gretl_prnh : PrnType;
import gretl.libgretlh : DatasetRaw, GretlType, Gretlopt;

extern (C) int gretl_cmd_init(Cmd * cmd);   /* the stand-in is in stubs.d */

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

bool same(const(char) * a, string b)
{
    return a !is null && strcmp(a, b.toStringz) == 0;
}

/* define a function from its text: header line, body lines, "end function" */
int define(string[] lines)
{
    ExecState st;
    Cmd cmd;
    gretl_error_clear();
    int err = gretl_start_compiling_function(lines[0].toStringz, null, null);

    memset(&st, 0, ExecState.sizeof);
    gretl_cmd_init(&cmd);
    st.cmd = &cmd;
    foreach (l; lines[1 .. $]) {
        char[] buf = l.dup ~ '\0';

        st.line = buf.ptr;
        err = gretl_function_append_line(&st);
        if (err) break;
    }
    return err;
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);
    int err = 0;

    /* ---- indentation of hansl code ---- */
    {
        int ti = 0, ni = 0;

        adjust_indent("loop i=1..10", &ti, &ni);
        check(ti == 0 && ni == 1, "adjust_indent: a loop line indents the next line");
        ti = ni;
        adjust_indent("x += i", &ti, &ni);
        check(ti == 1 && ni == 1, "adjust_indent: an ordinary line keeps the indent");
        adjust_indent("endloop", &ti, &ni);
        check(ti == 0 && ni == 0, "adjust_indent: endloop outdents this line and the next");
        ti = ni = 1;
        adjust_indent("else", &ti, &ni);
        check(ti == 0 && ni == 1, "adjust_indent: else outdents this line only");
        ti = ni = 0;
        adjust_indent("if x > 0", &ti, &ni);
        check(ni == 1, "adjust_indent: if");
    }

    /* ---- normalize_hansl ---- */
    {
        int e = 0;
        auto prn = gretl_print_new(PrnType.buffer, &e);

        normalize_hansl("if x>0\nloop i=1..3\ny=i\nendloop\nendif\n", 4, prn);
        const text = gretl_print_get_buffer(prn);

        check(strstr(text, "    loop i=1..3".ptr) !is null && strstr(text, "        y=i".ptr) !is null, "normalize_hansl: nested blocks are indented by the tab width");
        gretl_print_destroy(prn);
    }

    /* ---- return types ---- */
    {
        int e1 = 0, e2 = 0, e3 = 0, e4 = 0, e5 = 0;

        check(return_type_from_string("scalar", &e1) == GretlType.double_ && return_type_from_string("matrix", &e2) == GretlType.matrix
              && return_type_from_string("void", &e3) == GretlType.void_ && return_type_from_string("bundle", &e4) == GretlType.bundle
              && e1 == 0 && e2 == 0 && e3 == 0 && e4 == 0, "return_type_from_string: scalar, matrix, void, bundle");
        return_type_from_string("nonsense", &e5);
        check(e5 != 0, "return_type_from_string: an unknown type is an error");
    }

    /* ---- defining functions ---- */
    err = define(["function scalar addn(scalar x, int n[1:10:3] \"how many\", bool flag[0], const matrix m[null])",
                  "  scalar y = x + n",
                  "  return y",
                  "end function"]);
    check(err == 0, "a function with a defaulted int with range, a bool and an optional const matrix is defined");
    auto fun = get_user_function_by_name("addn");

    check(fun !is null && is_user_function("addn") != 0, "the function is found by name");
    if (fun !is null) {
        check(fn_n_params(fun) == 4, "four parameters");
        check(same(fn_param_name(fun, 0), "x") && same(fn_param_name(fun, 1), "n") && same(fn_param_name(fun, 2), "flag") && same(fn_param_name(fun, 3), "m"), "parameter names");
        check(fn_param_type(fun, 0) == GretlType.double_ && fn_param_type(fun, 1) == GretlType.int_ && fn_param_type(fun, 2) == GretlType.bool_ && fn_param_type(fun, 3) == GretlType.matrix,
              "parameter types");
        check(fn_param_has_default(fun, 1) != 0 && fabs(fn_param_default(fun, 1) - 3.0) < 1e-12 && fabs(fn_param_minval(fun, 1) - 1.0) < 1e-12 && fabs(fn_param_maxval(fun, 1) - 10.0) < 1e-12,
              "n: default 3, minimum 1, maximum 10");
        check(same(fn_param_descrip(fun, 1), "how many"), "n: description");
        check(fn_param_has_default(fun, 2) != 0 && fn_param_default(fun, 2) == 0.0, "flag: default 0");
        check(fn_param_optional(fun, 3) != 0, "m: optional (null default)");
        check(user_func_get_return_type(fun) == GretlType.double_, "return type scalar");
        check(gretl_function_get_signature(fun) !is null && strstr(gretl_function_get_signature(fun), "addn".ptr) !is null, "the signature mentions the name");
        int nl = 0;
        auto code = gretl_function_retrieve_code(fun, &nl);

        check(nl == 2 && code !is null && strstr(code[0], "scalar y = x + n".ptr) !is null && strstr(code[1], "return y".ptr) !is null, "the code (two lines) is retrieved");
    }

    err = define(["function void noop()", "end function"]);
    check(err != 0, "a function with no statements is refused");

    err = define(["function scalar addn(scalar z)", "  return z", "end function"]);
    // redefinition replaces: the same name with a new signature is allowed in a session
    auto fun2 = get_user_function_by_name("addn");

    check(err == 0 && fun2 !is null && fn_n_params(fun2) == 1, "redefining a function replaces it");
    check(n_user_functions() == 1 && same(user_function_name_by_index(0), "addn"), "one function in memory; its name by index");

    check(define(["function nonsense f(scalar x)", "  return x", "end function"]) != 0, "an unknown return type is an error");
    check(define(["function scalar g(scalar x, scalar x)", "  return x", "end function"]) != 0, "duplicate parameter names are an error");
    check(define(["function scalar h(scalar x", "  return x", "end function"]) != 0, "a missing right parenthesis is an error");
    check(define(["function scalar q(foo x)", "  return x", "end function"]) != 0, "an unknown parameter type is an error");
    check(define(["function scalar r(scalar x)", "  quit", "end function"]) != 0, "quit inside a function is an error");

    err = define(["function matrix lst(list L, series y, string s[null], scalar *p, bundle b)",
                  "  return {1}", "end function"]);
    auto lf = get_user_function_by_name("lst");

    check(err == 0 && lf !is null && fn_n_params(lf) == 5, "a function with list, series, optional string, pointer and bundle parameters");
    if (lf !is null) {
        check(fn_param_type(lf, 0) == GretlType.list && fn_param_type(lf, 1) == GretlType.series && fn_param_type(lf, 2) == GretlType.string_ && fn_param_type(lf, 3) == GretlType.scalar_ref
              && fn_param_type(lf, 4) == GretlType.bundle, "types: list, series, string, scalar reference, bundle");
        check(fn_param_optional(lf, 2) != 0, "the string parameter is optional (null default)");
    }

    /* ---- docstring, deleting ---- */
    check(user_func_add_docstr(fun2, "adds nothing") == 0 && same(user_func_get_docstr(fun2), "adds nothing"), "docstring set and read");
    check(define(["function addn delete"]) == 0 && get_user_function_by_name("addn") is null, "'function addn delete' removes the function");

    /* ---- call structures ---- */
    {
        err = define(["function scalar plus(scalar a, scalar b)", "  return a + b", "end function"]);
        auto pf = get_user_function_by_name("plus");
        auto fc = fncall_new(pf, 0);
        double x = 3, y = 4;

        check(fc !is null && fncall_get_return_type(fc) == GretlType.double_, "fncall_new: a call with the function's return type");
        check(push_function_arg(fc, null, null, GretlType.double_, &x) == 0 && push_function_arg(fc, null, null, GretlType.double_, &y) == 0, "two arguments pushed");
        check(push_function_arg(fc, null, null, GretlType.double_, &y) != 0, "a third argument is refused");
        fncall_destroy(fc);
        auto fc2 = user_func_get_fncall(pf);

        check(fc2 !is null, "user_func_get_fncall: a preserved call");
    }

    /* ---- packages ---- */
    {
        err = define(["function scalar pkgfun(scalar a)", "  return a", "end function"]);
        const fname = tempDir ~ "/test_gretl_func_pkg.gfn";

        if (exists(fname)) remove(fname);
        char *[1] pub = [cast(char *) "pkgfun".ptr];
        int perr = 0;
        auto pkg = function_package_new(fname.toStringz, pub.ptr, 1, null, 0, &perr);

        check(pkg !is null && perr == 0, "function_package_new for one public function");
        if (pkg !is null) {
            const e1 = function_package_set_properties(pkg, "author".ptr, "A. Tester".ptr, "version".ptr, "1.2".ptr, "date".ptr, "2026-01-02".ptr,
                                                       "description".ptr, "A test package".ptr, "email".ptr, "a@example.org".ptr,
                                                       "data-requirement".ptr, cast(int) DataReq.nodata_ok, "min-version".ptr, 20260101, null);

            check(e1 == 0, "package properties set");
            char * auth, pver, descr;

            function_package_get_properties(pkg, "author".ptr, &auth, "version".ptr, &pver, "description".ptr, &descr, null);
            check(same(auth, "A. Tester") && same(pver, "1.2") && same(descr, "A test package"), "string properties read back");
            check(same(function_package_get_name(pkg), "test_gretl_func_pkg"), "the package name is taken from the file name");
            check(function_package_get_version(pkg) > 1.19 && function_package_get_version(pkg) < 1.21, "version as a number");
            const w = function_package_write_file(pkg);

            check(w == 0 && exists(fname), "the .gfn file is written");
            if (exists(fname)) {
                FILE * fp = fopen(fname.toStringz, "r");
                char[16384] buf;
                const n = fread(buf.ptr, 1, buf.length - 1, fp);

                fclose(fp);
                buf[n] = 0;
                check(strstr(buf.ptr, "<gretl-function-package".ptr) !is null && strstr(buf.ptr, ">A. Tester</author>".ptr) !is null
                      && strstr(buf.ptr, "<version>1.2</version>".ptr) !is null && strstr(buf.ptr, "<gretl-function name=\"pkgfun\"".ptr) !is null
                      && strstr(buf.ptr, "return a".ptr) !is null, "the file holds the package header, the properties and the function code");
                char * desc, ver, date, author;
                int pdf = 0;
                const h = get_function_file_header(fname.toStringz, &desc, &ver, &date, &author, &pdf);

                check(h == 0 && same(desc, "A test package") && same(ver, "1.2") && same(date, "2026-01-02") && same(author, "A. Tester"), "get_function_file_header reads the header back");
                check(fabs(gfn_file_get_version(fname.toStringz) - 1.2) < 1e-9, "gfn_file_get_version");
                int ndeps = -1;
                auto deps = package_peek_dependencies(fname.toStringz, &ndeps);

                check(deps is null, "a package without dependencies reports none (null)");
                // load the package back from the file: the function is defined again from the XML
                function_package_unload_full_by_filename(fname.toStringz);
                check(get_user_function_by_name("pkgfun") is null, "a full unload removes the package and its functions");
                int lerr = 0;
                auto lp = get_function_package_by_filename(fname.toStringz, &lerr);

                check(lp !is null && lerr == 0, "get_function_package_by_filename loads the package from the file");
                auto lf2 = get_user_function_by_name("pkgfun");

                check(lf2 !is null && fn_n_params(lf2) == 1 && same(fn_param_name(lf2, 0), "a") && fn_param_type(lf2, 0) == GretlType.double_, "the function and its parameter are read from the XML");
                const(char) * lver;

                check(function_package_is_loaded(fname.toStringz, &lver) != 0 && same(lver, "1.2"), "function_package_is_loaded reports it with its version");
                int ln2 = 0;
                auto lcode = gretl_function_retrieve_code(lf2, &ln2);

                check(ln2 == 1 && strstr(lcode[0], "return a".ptr) !is null, "the code is read back");
                function_package_unload_full_by_filename(fname.toStringz);
                check(function_package_is_loaded(fname.toStringz, &lver) == 0 && get_user_function_by_name("pkgfun") is null, "the package and its function are unloaded");
                remove(fname);
            }
        }
    }

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
