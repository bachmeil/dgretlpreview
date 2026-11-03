/* Tests for genlex.d (genlex.c): the lexer module of 'genr': the tokenizer `lex`, the symbol tables and their look-up functions, the
   "pointerized" math function table, the symbol-name function getsymb, function overrides.

   `lex` is driven through parser_init (geneval5.d) with an anonymous expression and a small dataset (series y and xx). The expected token
   sequences are the lexer rules read from the C: numbers (decimal, exponent, hex), strings, operators (single and double character, dot
   operators, ++ / --, the transpose quote), brackets, the ternary operator, names resolved to series, constants, accessors ($...), built-in
   functions, lags (y(-1)); an undefined name stops the token stream with the "unknown" symbol.

   gretl_get_object_by_name (objstack.c, not ported yet) is given a stand-in that knows no saved objects.

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs, sqrt, exp, log, sin, tgamma;
import core.stdc.stdio : printf, setvbuf, stdout, _IOLBF;
import core.stdc.string : strcmp, strcpy;

import std.string : toStringz;

import gretl.dataset : create_new_dataset, dataset_set_time_series, destroy_dataset;
import gretl.genlex;
import gretl.genmainh : BundleDataIndex, ModelDataIndex, RetrievalIndex;
import gretl.genparseh;
import gretl.geneval5 : parser_init;
import gretl.libgretlh : DatasetRaw;


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

struct Tok {
    int sym;
    int idnum;
    double xval;
    string idstr;
}

DatasetRaw * dset;

/* the tokens of an expression up to EOT (or until the lexer stops at an unknown symbol) */
Tok[] tokens(string src, int * perr = null)
{
    Parser p;
    int done = 0;
    Tok[] ret;

    parser_init(&p, src.toStringz, dset, null, GenFlags.anon, GenSym.unk, &done);
    for (int i = 0; i < 60 && p.sym != GenSym.eot && !p.err; i++) {
        lex(&p);
        string s = (p.idstr is null) ? "" : cast(string) p.idstr[0 .. strlen_(p.idstr)].idup;

        ret ~= Tok(p.sym, p.idnum, p.xval, s);
        if (p.sym == GenSym.unk) break;
    }
    if (perr !is null) *perr = p.err;
    return ret;
}

size_t strlen_(const(char) * s)
{
    size_t n = 0;

    while (s[n]) n++;
    return n;
}

bool seq(Tok[] t, int[] syms)
{
    if (t.length != syms.length) return false;
    foreach (i, s; syms) if (t[i].sym != s) return false;
    return true;
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    dset = create_new_dataset(3, 10, 0);
    strcpy(dset.varname[1], "y");
    strcpy(dset.varname[2], "xx");
    dataset_set_time_series(dset, 1, 1900, 1);

    /* ---- tokens ---- */
    {
        auto t = tokens("sqrt(y) >= 1e-3 && xx != 0");

        check(seq(t, [GenFunc.f_sqrt, GenSym.g_lpr, GenSym.series, GenSym.g_rpr, GenSym.b_gte, GenSym.cnum, GenSym.b_and, GenSym.series, GenSym.b_neq, GenSym.cnum, GenSym.eot]),
              "a function call, a comparison with an exponent number, && and !=");
        check(t.length == 11 && t[2].idnum == 1 && t[2].idstr == "y" && fabs(t[5].xval - 0.001) < 1e-15 && t[7].idnum == 2 && t[7].idstr == "xx", "series y and xx are resolved to their ID numbers; 1e-3 is 0.001");
        auto h = tokens("0x1F + 2.5e2 - .5");

        check(seq(h, [GenSym.cnum, GenSym.b_add, GenSym.cnum, GenSym.b_sub, GenSym.cnum, GenSym.eot]) && h[0].xval == 31 && h[2].xval == 250 && h[4].xval == 0.5, "hex 0x1F is 31, 2.5e2 is 250, .5 is 0.5");
        auto l = tokens("y(-1) * 2");

        check(seq(l, [GenSym.lag, GenSym.g_lpr, GenSym.b_sub, GenSym.cnum, GenSym.g_rpr, GenSym.b_mul, GenSym.cnum, GenSym.eot]) && l[0].idstr == "y" && l[3].xval == 1, "a lag y(-1) is a LAG symbol followed by the parenthesized offset");
        auto s = tokens("\"hello\" ~ \"w\"");

        check(seq(s, [GenSym.cstr, GenSym.b_hcat, GenSym.cstr, GenSym.eot]) && s[0].idstr == "hello" && s[2].idstr == "w", "string literals and the horizontal concatenation operator ~");
        auto c = tokens("$pi + $nobs - NA");

        check(seq(c, [GenSym.con, GenSym.b_add, GenSym.dvar, GenSym.b_sub, GenSym.con, GenSym.eot]) && c[0].idnum == ConstId.pi && c[2].idnum == RetrievalIndex.nobs && c[4].idnum == ConstId.na,
              "$pi and NA are constants, $nobs is a dataset accessor");
        auto q = tokens("y ? 1 : 2");

        check(seq(q, [GenSym.series, GenSym.query, GenSym.cnum, GenSym.p_col, GenSym.cnum, GenSym.eot]), "the ternary operator");
        auto u = tokens("undefined_name + 1");

        check(u.length == 1 && u[0].sym == 0, "an undefined name stops the token stream (symbol 0)");
        auto d = tokens("xx == y || xx <= 3");

        check(seq(d, [GenSym.series, GenSym.b_eq, GenSym.series, GenSym.b_or, GenSym.series, GenSym.b_lte, GenSym.cnum, GenSym.eot]), "== , || and <=");
    }

    /* ---- look-up functions ---- */
    check(const_lookup("$pi") == ConstId.pi && const_lookup("NA") == ConstId.na && const_lookup("TRUE") == ConstId.true_ && const_lookup("nonsense") == 0, "const_lookup: $pi, NA, TRUE, and 0 for an unknown name");
    check(same(constname(ConstId.pi), "$pi") && same(constname(ConstId.na), "NA"), "constname gives the names back");
    check(mvar_lookup("$ess") == ModelDataIndex.ess && mvar_lookup("$uhat") == ModelDataIndex.uhat && mvar_lookup("$nobs") == 0 && same(mvarname(ModelDataIndex.ess), "$ess"), "mvar_lookup / mvarname for model accessors; $nobs is not a model accessor");
    check(same(dvarname(RetrievalIndex.nobs), "$nobs"), "dvarname");
    check(bvar_lookup("$model") == BundleDataIndex.model && same(bvarname(BundleDataIndex.model), "$model") && bvar_lookup("$ess") == 0, "bvar_lookup / bvarname for bundle accessors");
    check(same(dumname(DumId.diag), "diag") && same(dumname(DumId.null_), "null"), "dumname");
    check(function_lookup("sqrt") == GenFunc.f_sqrt && function_lookup("log") == GenFunc.f_log && function_lookup("no_such_function") == 0 && genr_function_word("sqrt") == GenFunc.f_sqrt && genr_function_word("zzz") == 0,
          "function_lookup and genr_function_word");
    check(is_function_alias("fft2") != 0 && is_function_alias("sqrt") == 0, "is_function_alias: fft2 is an alias, sqrt is not");
    check(same(gretl_function_complete("sqr"), "sqrt"), "gretl_function_complete: 'sqr' completes to sqrt");
    check(gen_func_count() > 300 && same(gen_func_name(0), "abs") && model_var_count() > 40 && bundle_var_count() == 3 && data_var_count() > 20 && gretl_const_count() == 12 && same(gretl_const_name(0), "$pi"),
          "the table sizes and first names");
    check(is_gretl_accessor("$nobs") != 0 && is_gretl_accessor("nobs") == 0, "is_gretl_accessor");

    /* ---- the pointerized math functions ---- */
    {
        alias DF = extern (C) double function(double);
        struct Case { string name; double x; double value; }
        Case[] cases = [Case("sqrt", 2.5, sqrt(2.5)), Case("sin", 2.5, sin(2.5)), Case("exp", 2.5, exp(2.5)), Case("log", 2.5, log(2.5)), Case("abs", -2.5, 2.5),
                        Case("gammafun", 2.5, tgamma(2.5)), Case("round", 2.4, 2.0), Case("sgn", -3.0, -1.0), Case("cnorm", 0.0, 0.5), Case("logistic", 0.0, 0.5)];
        bool ok = true;

        foreach (c; cases) {
            auto f = cast(DF) get_genr_function_pointer(function_lookup(c.name.toStringz));

            ok = ok && f !is null && fabs(f(c.x) - c.value) < 1e-12 * (1 + fabs(c.value));
        }
        check(ok, "the function pointers of sqrt, sin, exp, log, abs, gammafun, round, sgn, cnorm and logistic compute the right values");
        check(get_genr_function_pointer(function_lookup("diff")) is null, "a function not in the pointer table has no pointer");
    }

    /* ---- symbol names ---- */
    check(same(getsymb(GenSym.b_add), "+") && same(getsymb(GenSym.b_eq), "==") && same(getsymb(GenSym.series), "SERIES") && same(getsymb(GenFunc.f_sqrt), "sqrt"), "getsymb: operators, types and functions");

    /* ---- function overrides ---- */
    check(install_function_override("sqrt", "testpkg", null) == 0, "install_function_override accepts an override");
    check(delete_function_override("sqrt", "testpkg") == 1 && delete_function_override("sqrt", "testpkg") == 0, "delete_function_override: deleted once, then gone");

    destroy_dataset(dset);
    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
