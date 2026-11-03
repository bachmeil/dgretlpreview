/* Smoke tests for genmain.d (and genparseh.d, usermath.d, uservarh.d): the series-index lookups, the reserved and
   special words, extract_varname and the parser-independent entry points. realgen() is now the real evaluator
   (geneval5.d), which needs the lexer and parser (genlex.c, genparse.c: not ported), so the generate_*() wrappers
   that used to be tested here against a fake realgen() are left until those exist.

   Build and run with testing/run_tests.sh. */

import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;
import core.stdc.string : strcmp;

import gretl.dataioh : DataClearCode;
import gretl.dataset : create_new_dataset, destroy_dataset;
import gretl.genmain;
import gretl.genmainh : Generator;
import gretl.genparseh;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_matrix : gretl_matrix_alloc, gretl_matrix_free;
import gretl.gretl_matrixh : GretlMatrixRaw;
import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_new;
import gretl.gretl_prnh : PrnType;
import gretl.libgretlh;
import gretl.missingh : NADBL, na;

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

/* ---- tests ---- */

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    DatasetRaw * d = create_new_dataset(4, 5, 0);

    d.pd = 1;
    d.structure = TsCodes.cross_section;
    foreach (i, nm; ["const", "alpha", "beta", "gamma"]) {
        d.varname[i][0 .. nm.length] = nm[];
        d.varname[i][nm.length] = 0;
    }

    /* series lookups */
    check(series_index(d, "alpha") == 1 && series_index(d, "gamma") == 3, "series_index finds a series by name");
    check(series_index(d, "nosuch") == 4, "series_index of an unknown name is the next ID");
    check(series_index(d, "const") == 0, "series_index of 'const' is 0");
    check(series_index(d, "7up") == 4 && series_index(d, "") == 4 && series_index(d, null) == 4, "numeric, empty and NULL names give the next ID");
    check(series_index(null, "alpha") == -1, "series_index with no dataset is -1");
    check(series_greatest_index(d, "beta") == 2, "series_greatest_index");
    check(current_series_index(d, "beta") == 2 && current_series_index(d, "nosuch") == -1 && current_series_index(d, "") == -1,
          "current_series_index");
    check(caller_series_index(d, "alpha") == -1, "caller_series_index outside a function is -1");
    check(gretl_is_series("gamma", d) == 1 && gretl_is_series("nosuch", d) == 0, "gretl_is_series");

    /* words */
    check(genr_special_word("dummy") && genr_special_word("weekday") && genr_special_word("dummy:3") && !genr_special_word("dummy:") && !genr_special_word("alpha"),
          "genr_special_word");
    check(gretl_reserved_word("scalar") != 0 && gretl_reserved_word("NA") != 0 && gretl_reserved_word("ols") != 0 && gretl_reserved_word("alpha") == 0,
          "gretl_reserved_word (types, constants and command words)");
    char[32] nm = 0;
    int len;

    check(extract_varname(nm.ptr, "abc_1 + 2", &len) == 0 && len == 5 && strcmp(nm.ptr, "abc_1") == 0, "extract_varname");
    check(extract_varname(nm.ptr, "a123456789012345678901234567890123", &len) == GretlError.unkvar, "extract_varname: name too long");

    /* generate() and the parser entry points that do not need the (not yet ported) lexer and parser:
       the generate_*() wrappers are exercised once genlex.c and genparse.c exist */
    Prn * prn = gretl_print_new(PrnType.buffer, null);
    int err = 0;

    check(generate(null, d, GretlType.double_, Gretlopt.OPT_NONE, prn) == GretlError.args, "generate(NULL) is E_ARGS");
    Generator * g = genr_compile("genr time", d, GretlType.series, Gretlopt.OPT_NONE, prn, &err);
    check(g is null && err == GretlError.eqn, "genr_compile refuses the 'genr time' specials");

    /* the parser flag helpers */
    Parser pp;

    genr_set_na_check(&pp);
    check((pp.flags & GenFlags.natest) != 0, "genr_set_na_check");
    genr_unset_na_check(&pp);
    check((pp.flags & GenFlags.natest) == 0, "genr_unset_na_check");
    pp.flags = GenFlags.discard;
    check(genr_no_assign(&pp) == 1 && genr_is_autoregressive(&pp) == 0, "genr_no_assign / genr_is_autoregressive");
    pp.targ = GenSym.series;
    check(genr_get_output_type(&pp) == GretlType.series, "genr_get_output_type");

    /* genparseh helpers */
    check(unary_op(GenSym.u_neg) && !unary_op(GenSym.b_add) && binary_op(GenSym.b_add) && bool_comp(GenSym.b_gte) && dot_op(GenSym.b_dotgt),
          "operator classification helpers");
    check(func1_symb(GenFunc.f_abs) && func2_symb(GenFunc.f_min) && func3_symb(GenFunc.f_wmean) && funcn_symb(GenFunc.f_pdf) && !func1_symb(GenFunc.f_min),
          "function-arity classification helpers");

    gretl_print_destroy(prn);
    destroy_dataset(d);
    printf("\n%d checks, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
