/* Port of options.c */

module gretl.options;


import core.stdc.errno : errno;
import core.stdc.stdio : fprintf, fputc, fputs, sprintf, stderr;
import core.stdc.stdlib : atoi, free, malloc, realloc;
import core.stdc.string : strcmp, strlen, strncmp, strspn;

import gretl.dataset : gretl_function_depth;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_modelh : MODEL_COMMAND;
import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_get_buffer, gretl_print_new,
    gretl_print_reset_buffer, pprintf, pputc, pputs;
import gretl.gretl_prnh : PrnType;
import gretl.libgretlh;
import gretl.missingh : NADBL;
import gretl.optionsh : OptStatus;
import gretl.strutils : gretl_strdup, gretl_strdup_printf, gretl_strndup, strings_array_free,
    strings_array_new, strings_array_sort;
import gretl.gretl_errors : gretl_errmsg_set, gretl_errmsg_sprintf;
import gretl.gretl_utils : gretl_double_from_string, gretl_int_from_string;
import gretl.gretl_intl : gretl_pop_c_numeric_locale, gretl_push_c_numeric_locale;
import gretl.genmain : generate_int, generate_scalar, generate_string;
import gretl.gretl_commands : gretl_command_number, gretl_command_word;

/* ---- declarations of not-yet-ported items used below ----
   Replace each with an import when its defining file is ported. */

extern (C) {
}

/* commands for which --vcv (OPT_O) is applicable */
private bool vcv_opt_ok(GretlCmdIndex c)
{
    return MODEL_COMMAND(c) || c == GretlCmdIndex.add || c == GretlCmdIndex.omit;
}

/* commands for which --window (OPT_W) is applicable */
private bool window_opt_ok(GretlCmdIndex c)
{
    return MODEL_COMMAND(c) || c == GretlCmdIndex.var || c == GretlCmdIndex.vecm;
}

/* commands which support the $result accessor */
private bool yields_result(GretlCmdIndex c)
{
    return c == GretlCmdIndex.corr || c == GretlCmdIndex.freq ||
        c == GretlCmdIndex.summary;
}

/* commands for which --quiet (OPT_Q) is applicable */
private bool quiet_opt_ok(GretlCmdIndex c)
{
    return (MODEL_COMMAND(c) ||
            yields_result(c) ||
            c == GretlCmdIndex.add ||
            c == GretlCmdIndex.adf ||
            c == GretlCmdIndex.anova ||
            c == GretlCmdIndex.append ||
            c == GretlCmdIndex.bds ||
            c == GretlCmdIndex.bkw ||
            c == GretlCmdIndex.coeffsum ||
            c == GretlCmdIndex.chow ||
            c == GretlCmdIndex.coint2 ||
            c == GretlCmdIndex.corrgm ||
            c == GretlCmdIndex.cusum ||
            c == GretlCmdIndex.data ||
            c == GretlCmdIndex.difftest ||
            c == GretlCmdIndex.estimate ||
            c == GretlCmdIndex.fcast ||
            c == GretlCmdIndex.foreign ||
            c == GretlCmdIndex.fractint ||
            c == GretlCmdIndex.kpss ||
            c == GretlCmdIndex.makepkg ||
            c == GretlCmdIndex.modtest ||
            c == GretlCmdIndex.leverage ||
            c == GretlCmdIndex.levinlin ||
            c == GretlCmdIndex.loop ||
            c == GretlCmdIndex.mahal ||
            c == GretlCmdIndex.meantest ||
            c == GretlCmdIndex.normtest ||
            c == GretlCmdIndex.ols ||
            c == GretlCmdIndex.omit ||
            c == GretlCmdIndex.open ||
            c == GretlCmdIndex.panspec ||
            c == GretlCmdIndex.pkg ||
            c == GretlCmdIndex.qlrtest ||
            c == GretlCmdIndex.rename ||
            c == GretlCmdIndex.reset ||
            c == GretlCmdIndex.restrict ||
            c == GretlCmdIndex.rmplot ||
            c == GretlCmdIndex.smpl ||
            c == GretlCmdIndex.system ||
            c == GretlCmdIndex.var ||
            c == GretlCmdIndex.vartest ||
            c == GretlCmdIndex.vecm ||
            c == GretlCmdIndex.vif ||
            c == GretlCmdIndex.xcorrgm ||
            c == GretlCmdIndex.xtab);
}

/* --output (OPT_U) as attached to GNUPLOT */
private bool plot_output_opt_ok(GretlCmdIndex c)
{
    return (c == GretlCmdIndex.gnuplot ||
            c == GretlCmdIndex.plot ||
            c == GretlCmdIndex.bxplot ||
            c == GretlCmdIndex.hfplot ||
            c == GretlCmdIndex.panplot ||
            c == GretlCmdIndex.qqplot ||
            c == GretlCmdIndex.kdplot ||
            c == GretlCmdIndex.rmplot ||
            c == GretlCmdIndex.scatters ||
            c == GretlCmdIndex.gridplot);
}

/* --plot (OPT_U) as attached to CORR */
private bool cmd_plot_opt_ok(GretlCmdIndex c)
{
    return (c == GretlCmdIndex.corr ||
            c == GretlCmdIndex.corrgm ||
            c == GretlCmdIndex.cusum ||
            c == GretlCmdIndex.fcast ||
            c == GretlCmdIndex.freq ||
            c == GretlCmdIndex.hurst ||
            c == GretlCmdIndex.leverage ||
            c == GretlCmdIndex.pergm ||
            c == GretlCmdIndex.qlrtest ||
            c == GretlCmdIndex.xcorrgm);
}

/* --outbuf (OPT_b) as attached to GNUPLOT */
private bool plot_outbuf_opt_ok(GretlCmdIndex c)
{
    return (plot_output_opt_ok(c) || cmd_plot_opt_ok(c));
}

struct GretlOption {
    GretlCmdIndex ci;    /* command index (gives context) */
    Gretlopt o;          /* index of integer type */
    string longopt;      /* "--"-style string representation of option */
    ubyte parminfo;      /* 0 = option can never take a parameter,
                            1 = option may take a parameter,
                            2 = option requires a parameter
                         */
}

private struct FlagMatch {
    Gretlopt o;
    ubyte c;
}

/* Below: This is used as a one-way mapping from the long form to the
   index (e.g. OPT_Q), so a given index can have more than one
   long-form counterpart depending on the context.  The last field
   indicates whether the given option does not accept (0), accepts
   (1), or requires (2) an associated parameter value.
*/

immutable GretlOption[] gretl_opts = [
    GretlOption(GretlCmdIndex.add, Gretlopt.OPT_A, "auto", 2),
    GretlOption(GretlCmdIndex.add, Gretlopt.OPT_B, "both", 0),
    GretlOption(GretlCmdIndex.add, Gretlopt.OPT_I, "silent", 0),
    GretlOption(GretlCmdIndex.add, Gretlopt.OPT_L, "lm", 0),
    GretlOption(GretlCmdIndex.adf, Gretlopt.OPT_N, "nc", 0),
    GretlOption(GretlCmdIndex.adf, Gretlopt.OPT_C, "c", 0),
    GretlOption(GretlCmdIndex.adf, Gretlopt.OPT_D, "seasonals", 0),
    GretlOption(GretlCmdIndex.adf, Gretlopt.OPT_R, "ctt", 0),
    GretlOption(GretlCmdIndex.adf, Gretlopt.OPT_T, "ct", 0),
    GretlOption(GretlCmdIndex.adf, Gretlopt.OPT_V, "verbose", 1),
    GretlOption(GretlCmdIndex.adf, Gretlopt.OPT_F, "difference", 0),
    GretlOption(GretlCmdIndex.adf, Gretlopt.OPT_E, "test-down", 1),
    GretlOption(GretlCmdIndex.adf, Gretlopt.OPT_G, "gls", 0),
    GretlOption(GretlCmdIndex.adf, Gretlopt.OPT_U, "perron-qu", 0),
    GretlOption(GretlCmdIndex.ar1, Gretlopt.OPT_B, "no-corc", 0),
    GretlOption(GretlCmdIndex.ar1, Gretlopt.OPT_H, "hilu", 0),
    GretlOption(GretlCmdIndex.ar1, Gretlopt.OPT_P, "pwe", 0),
    GretlOption(GretlCmdIndex.ar1, Gretlopt.OPT_L, "loose", 0),
    GretlOption(GretlCmdIndex.append, Gretlopt.OPT_A, "all-cols", 0),
    GretlOption(GretlCmdIndex.append, Gretlopt.OPT_T, "time-series", 0),
    GretlOption(GretlCmdIndex.append, Gretlopt.OPT_R, "rowoffset", 2),
    GretlOption(GretlCmdIndex.append, Gretlopt.OPT_C, "coloffset", 2),
    GretlOption(GretlCmdIndex.append, Gretlopt.OPT_F, "fixed-cols", 2),
    GretlOption(GretlCmdIndex.append, Gretlopt.OPT_M, "rowmask", 2),
    GretlOption(GretlCmdIndex.append, Gretlopt.OPT_L, "cols", 2),
    GretlOption(GretlCmdIndex.append, Gretlopt.OPT_S, "sheet", 2),
    GretlOption(GretlCmdIndex.append, Gretlopt.OPT_V, "verbose", 0),
    GretlOption(GretlCmdIndex.append, Gretlopt.OPT_U, "update-overlap", 0),
    GretlOption(GretlCmdIndex.append, Gretlopt.OPT_X, "fixed-sample", 0),
    GretlOption(GretlCmdIndex.append, Gretlopt.OPT_K, "frompkg", 2),
    GretlOption(GretlCmdIndex.arma, Gretlopt.OPT_A, "as154", 0),
    GretlOption(GretlCmdIndex.arma, Gretlopt.OPT_C, "conditional", 0),
    GretlOption(GretlCmdIndex.arma, Gretlopt.OPT_E, "save-ehat", 0),
    GretlOption(GretlCmdIndex.arma, Gretlopt.OPT_G, "opg", 0),
    GretlOption(GretlCmdIndex.arma, Gretlopt.OPT_H, "hessian", 0),
    GretlOption(GretlCmdIndex.arma, Gretlopt.OPT_K, "kalman", 0),
    GretlOption(GretlCmdIndex.arma, Gretlopt.OPT_L, "lbfgs", 0),
    GretlOption(GretlCmdIndex.arma, Gretlopt.OPT_N, "nc", 0),
    GretlOption(GretlCmdIndex.arma, Gretlopt.OPT_V, "verbose", 0),
    GretlOption(GretlCmdIndex.arma, Gretlopt.OPT_X, "x-13arima", 0),
    GretlOption(GretlCmdIndex.arma, Gretlopt.OPT_X, "x-12-arima", 0), /* legacy alias */
    GretlOption(GretlCmdIndex.arma, Gretlopt.OPT_Y, "y-diff-only", 0),
    GretlOption(GretlCmdIndex.arma, Gretlopt.OPT_R, "robust", 0),
    GretlOption(GretlCmdIndex.arma, Gretlopt.OPT_S, "stdx", 0),
    GretlOption(GretlCmdIndex.arma, Gretlopt.OPT_Z, "lagselect", 0),
    GretlOption(GretlCmdIndex.bds, Gretlopt.OPT_B, "boot", 2),
    GretlOption(GretlCmdIndex.bds, Gretlopt.OPT_C, "corr1", 2),
    GretlOption(GretlCmdIndex.bds, Gretlopt.OPT_S, "sdcrit", 2),
    GretlOption(GretlCmdIndex.bds, Gretlopt.OPT_X, "matrix", 2),
    GretlOption(GretlCmdIndex.biprobit, Gretlopt.OPT_G, "opg", 0),
    GretlOption(GretlCmdIndex.biprobit, Gretlopt.OPT_R, "robust", 0),
    GretlOption(GretlCmdIndex.biprobit, Gretlopt.OPT_V, "verbose", 0),
    GretlOption(GretlCmdIndex.biprobit, Gretlopt.OPT_C, "cluster", 2),
    GretlOption(GretlCmdIndex.biprobit, Gretlopt.OPT_X, "save-xbeta", 0),
    GretlOption(GretlCmdIndex.bxplot, Gretlopt.OPT_O, "notches", 0),
    GretlOption(GretlCmdIndex.bxplot, Gretlopt.OPT_L, "outliers", 2),
    GretlOption(GretlCmdIndex.bxplot, Gretlopt.OPT_P, "panel", 0),
    GretlOption(GretlCmdIndex.bxplot, Gretlopt.OPT_X, "matrix", 2),
    GretlOption(GretlCmdIndex.bxplot, Gretlopt.OPT_Z, "factorized", 0),
    GretlOption(GretlCmdIndex.bxplot, Gretlopt.OPT_B, "whiskerbars", 0),
    GretlOption(GretlCmdIndex.chow, Gretlopt.OPT_D, "dummy", 0),
    GretlOption(GretlCmdIndex.chow, Gretlopt.OPT_L, "limit-to", 2),
    GretlOption(GretlCmdIndex.clear, Gretlopt.OPT_A, "all", 0),
    GretlOption(GretlCmdIndex.clear, Gretlopt.OPT_D, "dataset", 0),
    GretlOption(GretlCmdIndex.clear, Gretlopt.OPT_F, "functions", 0),
    GretlOption(GretlCmdIndex.coint, Gretlopt.OPT_D, "seasonals", 0),
    GretlOption(GretlCmdIndex.coint, Gretlopt.OPT_E, "test-down", 1),
    GretlOption(GretlCmdIndex.coint, Gretlopt.OPT_N, "nc", 0),
    GretlOption(GretlCmdIndex.coint, Gretlopt.OPT_R, "ctt", 0),
    GretlOption(GretlCmdIndex.coint, Gretlopt.OPT_S, "skip-df", 0),
    GretlOption(GretlCmdIndex.coint, Gretlopt.OPT_T, "ct", 0),
    GretlOption(GretlCmdIndex.coint, Gretlopt.OPT_V, "verbose", 0),
    GretlOption(GretlCmdIndex.coint, Gretlopt.OPT_I, "silent", 0),
    GretlOption(GretlCmdIndex.coint2, Gretlopt.OPT_A, "crt", 0),
    GretlOption(GretlCmdIndex.coint2, Gretlopt.OPT_D, "seasonals", 0),
    GretlOption(GretlCmdIndex.coint2, Gretlopt.OPT_N, "nc", 0),
    GretlOption(GretlCmdIndex.coint2, Gretlopt.OPT_R, "rc", 0),
    GretlOption(GretlCmdIndex.coint2, Gretlopt.OPT_C, "uc", 0),
    GretlOption(GretlCmdIndex.coint2, Gretlopt.OPT_T, "ct", 0),
    GretlOption(GretlCmdIndex.coint2, Gretlopt.OPT_S, "silent", 0),
    GretlOption(GretlCmdIndex.coint2, Gretlopt.OPT_V, "verbose", 0),
    GretlOption(GretlCmdIndex.coint2, Gretlopt.OPT_Y, "asy", 0),
    GretlOption(GretlCmdIndex.corr, Gretlopt.OPT_K, "kendall", 0),
    GretlOption(GretlCmdIndex.corr, Gretlopt.OPT_S, "spearman", 0),
    GretlOption(GretlCmdIndex.corr, Gretlopt.OPT_N, "uniform", 0),
    GretlOption(GretlCmdIndex.corr, Gretlopt.OPT_V, "verbose", 0),
    GretlOption(GretlCmdIndex.corr, Gretlopt.OPT_U, "plot", 2),
    GretlOption(GretlCmdIndex.corr, Gretlopt.OPT_X, "matrix", 2),
    GretlOption(GretlCmdIndex.corr, Gretlopt.OPT_T, "triangle", 0),
    GretlOption(GretlCmdIndex.corrgm, Gretlopt.OPT_A, "acf-only", 0),
    GretlOption(GretlCmdIndex.corrgm, Gretlopt.OPT_B, "bartlett", 0),
    GretlOption(GretlCmdIndex.corrgm, Gretlopt.OPT_S, "silent", 0), /* legacy (= quiet) */
    GretlOption(GretlCmdIndex.cusum, Gretlopt.OPT_R, "squares", 0),
    GretlOption(GretlCmdIndex.data, Gretlopt.OPT_C, "compact", 2),
    GretlOption(GretlCmdIndex.data, Gretlopt.OPT_O, "odbc", 0),
    GretlOption(GretlCmdIndex.data, Gretlopt.OPT_N, "name", 2),
    GretlOption(GretlCmdIndex.data, Gretlopt.OPT_V, "verbose", 0),
    GretlOption(GretlCmdIndex.data, Gretlopt.OPT_F, "no-align", 0),
    GretlOption(GretlCmdIndex.datamod, Gretlopt.OPT_P, "preserve", 0),
    GretlOption(GretlCmdIndex.datamod, Gretlopt.OPT_T, "panel-time", 0),
    GretlOption(GretlCmdIndex.datamod, Gretlopt.OPT_W, "weekstart", 2),
    GretlOption(GretlCmdIndex.datamod, Gretlopt.OPT_R, "repday", 2),
    GretlOption(GretlCmdIndex.deleet, Gretlopt.OPT_D, "db", 0),
    GretlOption(GretlCmdIndex.deleet, Gretlopt.OPT_F, "force", 0),
    GretlOption(GretlCmdIndex.deleet, Gretlopt.OPT_L, "list", 0),
    GretlOption(GretlCmdIndex.deleet, Gretlopt.OPT_T, "type", 2),
    GretlOption(GretlCmdIndex.difftest, Gretlopt.OPT_D, "split-by", 2),
    GretlOption(GretlCmdIndex.difftest, Gretlopt.OPT_G, "sign", 0),
    GretlOption(GretlCmdIndex.difftest, Gretlopt.OPT_R, "rank-sum", 0),
    GretlOption(GretlCmdIndex.difftest, Gretlopt.OPT_I, "signed-rank", 0),
    GretlOption(GretlCmdIndex.difftest, Gretlopt.OPT_V, "verbose", 0),
    GretlOption(GretlCmdIndex.discrete, Gretlopt.OPT_R, "reverse", 0),
    GretlOption(GretlCmdIndex.dpanel, Gretlopt.OPT_A, "asymptotic", 0),
    GretlOption(GretlCmdIndex.dpanel, Gretlopt.OPT_D, "time-dummies", 1),
    GretlOption(GretlCmdIndex.dpanel, Gretlopt.OPT_K, "keep-extra", 0),
    GretlOption(GretlCmdIndex.dpanel, Gretlopt.OPT_L, "system", 0),
    GretlOption(GretlCmdIndex.dpanel, Gretlopt.OPT_T, "two-step", 0),
    GretlOption(GretlCmdIndex.dpanel, Gretlopt.OPT_V, "verbose", 0),
    GretlOption(GretlCmdIndex.dpanel, Gretlopt.OPT_X, "dpdstyle", 0),
    GretlOption(GretlCmdIndex.dpanel, Gretlopt.OPT_C, "collapse", 0),
    GretlOption(GretlCmdIndex.dummify, Gretlopt.OPT_F, "drop-first", 0),
    GretlOption(GretlCmdIndex.dummify, Gretlopt.OPT_L, "drop-last", 0),
    GretlOption(GretlCmdIndex.duration, Gretlopt.OPT_B, "weibull", 0),
    GretlOption(GretlCmdIndex.duration, Gretlopt.OPT_E, "exponential", 0),
    GretlOption(GretlCmdIndex.duration, Gretlopt.OPT_L, "loglogistic", 0),
    GretlOption(GretlCmdIndex.duration, Gretlopt.OPT_Z, "lognormal", 0),
    GretlOption(GretlCmdIndex.duration, Gretlopt.OPT_M, "medians", 0),
    GretlOption(GretlCmdIndex.duration, Gretlopt.OPT_G, "opg", 0),
    GretlOption(GretlCmdIndex.duration, Gretlopt.OPT_R, "robust", 0),
    GretlOption(GretlCmdIndex.duration, Gretlopt.OPT_C, "cluster", 2),
    GretlOption(GretlCmdIndex.duration, Gretlopt.OPT_V, "verbose", 0),
    GretlOption(GretlCmdIndex.eqnprint, Gretlopt.OPT_O, "complete", 0),
    GretlOption(GretlCmdIndex.eqnprint, Gretlopt.OPT_T, "t-ratios", 0),
    GretlOption(GretlCmdIndex.eqnprint, Gretlopt.OPT_U, "output", 2),
    GretlOption(GretlCmdIndex.tabprint, Gretlopt.OPT_O, "complete", 0),
    GretlOption(GretlCmdIndex.tabprint, Gretlopt.OPT_C, "csv", 0),
    GretlOption(GretlCmdIndex.tabprint, Gretlopt.OPT_R, "rtf", 0),
    GretlOption(GretlCmdIndex.tabprint, Gretlopt.OPT_T, "format", 2),
    GretlOption(GretlCmdIndex.tabprint, Gretlopt.OPT_U, "output", 2),
    GretlOption(GretlCmdIndex.equation, Gretlopt.OPT_M, "multi", 0),
    GretlOption(GretlCmdIndex.estimate, Gretlopt.OPT_I, "iterate", 0),
    GretlOption(GretlCmdIndex.estimate, Gretlopt.OPT_M, "geomean", 0),
    GretlOption(GretlCmdIndex.estimate, Gretlopt.OPT_N, "no-df-corr", 0),
    GretlOption(GretlCmdIndex.estimate, Gretlopt.OPT_V, "verbose", 0),
    GretlOption(GretlCmdIndex.estimate, Gretlopt.OPT_W, "window", 0),
    GretlOption(GretlCmdIndex.fcast, Gretlopt.OPT_L, "all-probs", 0),
    GretlOption(GretlCmdIndex.fcast, Gretlopt.OPT_D, "dynamic", 0),
    GretlOption(GretlCmdIndex.fcast, Gretlopt.OPT_N, "no-stats", 0),
    GretlOption(GretlCmdIndex.fcast, Gretlopt.OPT_T, "stats-only", 0),
    GretlOption(GretlCmdIndex.fcast, Gretlopt.OPT_S, "static", 0),
    GretlOption(GretlCmdIndex.fcast, Gretlopt.OPT_R, "recursive", 0),
    GretlOption(GretlCmdIndex.fcast, Gretlopt.OPT_R, "rolling", 0), /* legacy alias */
    GretlOption(GretlCmdIndex.fcast, Gretlopt.OPT_O, "out-of-sample", 0),
    GretlOption(GretlCmdIndex.fcast, Gretlopt.OPT_I, "integrate", 0),
    GretlOption(GretlCmdIndex.fcast, Gretlopt.OPT_X, "exponentiate", 1),
    GretlOption(GretlCmdIndex.foreign, Gretlopt.OPT_D, "send-data", 1),
    GretlOption(GretlCmdIndex.foreign, Gretlopt.OPT_V, "verbose", 0),
    GretlOption(GretlCmdIndex.foreign, Gretlopt.OPT_F, "frame", 0),
    GretlOption(GretlCmdIndex.foreign, Gretlopt.OPT_I, "io-funcs", 2),
    GretlOption(GretlCmdIndex.fractint, Gretlopt.OPT_G, "gph", 0),
    GretlOption(GretlCmdIndex.fractint, Gretlopt.OPT_A, "all", 0),
    GretlOption(GretlCmdIndex.freq, Gretlopt.OPT_Q, "quiet", 0), /* legacy */
    GretlOption(GretlCmdIndex.freq, Gretlopt.OPT_O, "gamma", 0),
    GretlOption(GretlCmdIndex.freq, Gretlopt.OPT_S, "silent", 0),
    GretlOption(GretlCmdIndex.freq, Gretlopt.OPT_Z, "normal", 0),
    GretlOption(GretlCmdIndex.freq, Gretlopt.OPT_N, "nbins", 2),
    GretlOption(GretlCmdIndex.freq, Gretlopt.OPT_M, "min", 2),
    GretlOption(GretlCmdIndex.freq, Gretlopt.OPT_W, "binwidth", 2),
    GretlOption(GretlCmdIndex.freq, Gretlopt.OPT_X, "matrix", 2),
    GretlOption(GretlCmdIndex.garch, Gretlopt.OPT_A, "arma-init", 0),
    GretlOption(GretlCmdIndex.garch, Gretlopt.OPT_F, "fcp", 0),
    GretlOption(GretlCmdIndex.garch, Gretlopt.OPT_N, "nc", 0),
    GretlOption(GretlCmdIndex.garch, Gretlopt.OPT_R, "robust", 0),
    GretlOption(GretlCmdIndex.garch, Gretlopt.OPT_V, "verbose", 0),
    GretlOption(GretlCmdIndex.garch, Gretlopt.OPT_Z, "stdresid", 0),
    GretlOption(GretlCmdIndex.gibbs, Gretlopt.OPT_K, "keep", 0),
    GretlOption(GretlCmdIndex.gibbs, Gretlopt.OPT_Q, "quiet", 0),
    GretlOption(GretlCmdIndex.gibbs, Gretlopt.OPT_T, "thinning", 2),
    GretlOption(GretlCmdIndex.gmm, Gretlopt.OPT_I, "iterate", 0),
    GretlOption(GretlCmdIndex.gmm, Gretlopt.OPT_L, "lbfgs", 1),
    GretlOption(GretlCmdIndex.gmm, Gretlopt.OPT_T, "two-step", 0),
    GretlOption(GretlCmdIndex.gmm, Gretlopt.OPT_V, "verbose", 0),
    GretlOption(GretlCmdIndex.gnuplot, Gretlopt.OPT_I, "input", 2),
    GretlOption(GretlCmdIndex.gnuplot, Gretlopt.OPT_O, "with-lines", 1),
    GretlOption(GretlCmdIndex.gnuplot, Gretlopt.OPT_F, "fit", 2),
    GretlOption(GretlCmdIndex.gnuplot, Gretlopt.OPT_M, "with-impulses", 1),
    GretlOption(GretlCmdIndex.gnuplot, Gretlopt.OPT_P, "with-lp", 1),
    GretlOption(GretlCmdIndex.gnuplot, Gretlopt.OPT_B, "with-boxes", 1),
    GretlOption(GretlCmdIndex.gnuplot, Gretlopt.OPT_Q, "with-steps", 1),
    GretlOption(GretlCmdIndex.gnuplot, Gretlopt.OPT_T, "time-series", 0),
    GretlOption(GretlCmdIndex.gnuplot, Gretlopt.OPT_Z, "factorized", 1),
    GretlOption(GretlCmdIndex.gnuplot, Gretlopt.OPT_Z, "dummy", 1), /* alias */
    GretlOption(GretlCmdIndex.gnuplot, Gretlopt.OPT_C, "control", 0),
    GretlOption(GretlCmdIndex.gnuplot, Gretlopt.OPT_U, "output", 2),
    GretlOption(GretlCmdIndex.gnuplot, Gretlopt.OPT_b, "outbuf", 2),
    GretlOption(GretlCmdIndex.gnuplot, Gretlopt.OPT_b, "buffer", 2), /* legacy alias */
    GretlOption(GretlCmdIndex.gnuplot, Gretlopt.OPT_i, "inbuf", 2),
    GretlOption(GretlCmdIndex.gnuplot, Gretlopt.OPT_Y, "single-yaxis", 0),
    GretlOption(GretlCmdIndex.gnuplot, Gretlopt.OPT_X, "matrix", 2),
    GretlOption(GretlCmdIndex.gnuplot, Gretlopt.OPT_N, "band", 2),
    GretlOption(GretlCmdIndex.gnuplot, Gretlopt.OPT_J, "band-style", 2),
    GretlOption(GretlCmdIndex.gnuplot, Gretlopt.OPT_a, "bands", 2),
    GretlOption(GretlCmdIndex.gnuplot, Gretlopt.OPT_W, "font", 2),
    GretlOption(GretlCmdIndex.gnuplot, Gretlopt.OPT_L, "ylogscale", 1),
    GretlOption(GretlCmdIndex.gnuplot, Gretlopt.OPT_D, "y2axis", 2),
    GretlOption(GretlCmdIndex.gnuplot, Gretlopt.OPT_E, "alpha", 2), /* experimental */
    GretlOption(GretlCmdIndex.graphpg, Gretlopt.OPT_M, "monochrome", 0),
    GretlOption(GretlCmdIndex.graphpg, Gretlopt.OPT_O, "output", 2),
    GretlOption(GretlCmdIndex.gridplot, Gretlopt.OPT_F, "fontsize", 2),
    GretlOption(GretlCmdIndex.gridplot, Gretlopt.OPT_W, "width", 2),
    GretlOption(GretlCmdIndex.gridplot, Gretlopt.OPT_H, "height", 2),
    GretlOption(GretlCmdIndex.gridplot, Gretlopt.OPT_R, "rows", 2),
    GretlOption(GretlCmdIndex.gridplot, Gretlopt.OPT_C, "cols", 2),
    GretlOption(GretlCmdIndex.gridplot, Gretlopt.OPT_L, "layout", 2),
    GretlOption(GretlCmdIndex.gridplot, Gretlopt.OPT_T, "title", 2),
    GretlOption(GretlCmdIndex.gridplot, Gretlopt.OPT_U, "output", 2),
    GretlOption(GretlCmdIndex.heckit, Gretlopt.OPT_M, "ml", 0),
    GretlOption(GretlCmdIndex.heckit, Gretlopt.OPT_G, "opg", 0),
    GretlOption(GretlCmdIndex.heckit, Gretlopt.OPT_R, "robust", 0),
    GretlOption(GretlCmdIndex.heckit, Gretlopt.OPT_C, "cluster", 2),
    GretlOption(GretlCmdIndex.heckit, Gretlopt.OPT_T, "two-step", 0),
    GretlOption(GretlCmdIndex.heckit, Gretlopt.OPT_V, "verbose", 0),
    GretlOption(GretlCmdIndex.help, Gretlopt.OPT_F, "func", 0),
    GretlOption(GretlCmdIndex.hfplot, Gretlopt.OPT_O, "with-lines", 0),
    GretlOption(GretlCmdIndex.hfplot, Gretlopt.OPT_T, "time-series", 0),
    GretlOption(GretlCmdIndex.hsk, Gretlopt.OPT_N, "no-squares", 0),
    GretlOption(GretlCmdIndex.include, Gretlopt.OPT_F, "force", 0),
    GretlOption(GretlCmdIndex.info, Gretlopt.OPT_F, "from-file", 2),
    GretlOption(GretlCmdIndex.info, Gretlopt.OPT_T, "to-file", 2),
    GretlOption(GretlCmdIndex.intreg, Gretlopt.OPT_G, "opg", 0),
    GretlOption(GretlCmdIndex.intreg, Gretlopt.OPT_R, "robust", 0),
    GretlOption(GretlCmdIndex.intreg, Gretlopt.OPT_C, "cluster", 2),
    GretlOption(GretlCmdIndex.intreg, Gretlopt.OPT_V, "verbose", 0),
    GretlOption(GretlCmdIndex.ivreg, Gretlopt.OPT_G, "gmm", 0),
    GretlOption(GretlCmdIndex.ivreg, Gretlopt.OPT_I, "iterate", 0),
    GretlOption(GretlCmdIndex.ivreg, Gretlopt.OPT_L, "liml", 0),
    GretlOption(GretlCmdIndex.ivreg, Gretlopt.OPT_N, "no-df-corr", 0),
    GretlOption(GretlCmdIndex.ivreg, Gretlopt.OPT_R, "robust", 0),
    GretlOption(GretlCmdIndex.ivreg, Gretlopt.OPT_S, "save", 0),
    GretlOption(GretlCmdIndex.ivreg, Gretlopt.OPT_T, "two-step", 0),
    GretlOption(GretlCmdIndex.ivreg, Gretlopt.OPT_H, "weights", 2),
    GretlOption(GretlCmdIndex.ivreg, Gretlopt.OPT_X, "no-tests", 1),
    GretlOption(GretlCmdIndex.ivreg, Gretlopt.OPT_C, "cluster", 2),
    GretlOption(GretlCmdIndex.ivreg, Gretlopt.OPT_M, "matrix-diff", 0),
    GretlOption(GretlCmdIndex.join, Gretlopt.OPT_I, "ikey", 2),
    GretlOption(GretlCmdIndex.join, Gretlopt.OPT_O, "okey", 2),
    GretlOption(GretlCmdIndex.join, Gretlopt.OPT_F, "filter", 2),
    GretlOption(GretlCmdIndex.join, Gretlopt.OPT_A, "aggr", 2),
    GretlOption(GretlCmdIndex.join, Gretlopt.OPT_D, "data", 2),
    GretlOption(GretlCmdIndex.join, Gretlopt.OPT_T, "tconv-fmt", 2),
    GretlOption(GretlCmdIndex.join, Gretlopt.OPT_K, "tkey", 2),
    GretlOption(GretlCmdIndex.join, Gretlopt.OPT_X, "tconvert", 2),
    GretlOption(GretlCmdIndex.join, Gretlopt.OPT_H, "no-header", 0),
    GretlOption(GretlCmdIndex.join, Gretlopt.OPT_V, "verbose", 0),
    GretlOption(GretlCmdIndex.join, Gretlopt.OPT_P, "pd", 2), /* undocumented: is it wanted? */
    GretlOption(GretlCmdIndex.join, Gretlopt.OPT_R, "frompkg", 2),
    GretlOption(GretlCmdIndex.kdplot, Gretlopt.OPT_O, "alt", 0),
    GretlOption(GretlCmdIndex.kdplot, Gretlopt.OPT_S, "scale", 2),
    GretlOption(GretlCmdIndex.kpss, Gretlopt.OPT_T, "trend", 0),
    GretlOption(GretlCmdIndex.kpss, Gretlopt.OPT_D, "seasonals", 0),
    GretlOption(GretlCmdIndex.kpss, Gretlopt.OPT_V, "verbose", 0),
    GretlOption(GretlCmdIndex.kpss, Gretlopt.OPT_F, "difference", 0),
    GretlOption(GretlCmdIndex.lags, Gretlopt.OPT_L, "bylag", 0),
    GretlOption(GretlCmdIndex.leverage, Gretlopt.OPT_S, "save", 0),
    GretlOption(GretlCmdIndex.leverage, Gretlopt.OPT_O, "overwrite", 0),
    GretlOption(GretlCmdIndex.levinlin, Gretlopt.OPT_N, "nc", 0),
    GretlOption(GretlCmdIndex.levinlin, Gretlopt.OPT_T, "ct", 0),
    GretlOption(GretlCmdIndex.levinlin, Gretlopt.OPT_V, "verbose", 0),
    GretlOption(GretlCmdIndex.makepkg, Gretlopt.OPT_D, "dtd", 2),
    GretlOption(GretlCmdIndex.makepkg, Gretlopt.OPT_I, "index", 0),
    GretlOption(GretlCmdIndex.makepkg, Gretlopt.OPT_T, "translations", 0),
    GretlOption(GretlCmdIndex.markers, Gretlopt.OPT_D, "delete", 0),
    GretlOption(GretlCmdIndex.markers, Gretlopt.OPT_F, "from-file", 2),
    GretlOption(GretlCmdIndex.markers, Gretlopt.OPT_T, "to-file", 2),
    GretlOption(GretlCmdIndex.markers, Gretlopt.OPT_A, "to-array", 2),
    GretlOption(GretlCmdIndex.markers, Gretlopt.OPT_R, "from-array", 2),
    GretlOption(GretlCmdIndex.markers, Gretlopt.OPT_S, "from-series", 2),
    GretlOption(GretlCmdIndex.modtest, Gretlopt.OPT_A, "autocorr", 0),
    GretlOption(GretlCmdIndex.modtest, Gretlopt.OPT_B, "breusch-pagan", 0),
    GretlOption(GretlCmdIndex.modtest, Gretlopt.OPT_C, "comfac", 0),
    GretlOption(GretlCmdIndex.modtest, Gretlopt.OPT_D, "xdepend", 0),
    GretlOption(GretlCmdIndex.modtest, Gretlopt.OPT_H, "arch", 0),
    GretlOption(GretlCmdIndex.modtest, Gretlopt.OPT_L, "logs", 0),
    GretlOption(GretlCmdIndex.modtest, Gretlopt.OPT_N, "normality", 0),
    GretlOption(GretlCmdIndex.modtest, Gretlopt.OPT_S, "squares", 0),
    GretlOption(GretlCmdIndex.modtest, Gretlopt.OPT_P, "panel", 0),
    GretlOption(GretlCmdIndex.modtest, Gretlopt.OPT_R, "robust", 0),
    GretlOption(GretlCmdIndex.modtest, Gretlopt.OPT_W, "white", 0),
    GretlOption(GretlCmdIndex.modtest, Gretlopt.OPT_X, "white-nocross", 0),
    GretlOption(GretlCmdIndex.modtest, Gretlopt.OPT_I, "silent", 0),
    GretlOption(GretlCmdIndex.modtest, Gretlopt.OPT_U, "univariate", 0),
    GretlOption(GretlCmdIndex.mpi, Gretlopt.OPT_D, "send-data", 1),
    GretlOption(GretlCmdIndex.mpi, Gretlopt.OPT_M, "send-metadata", 0),
    GretlOption(GretlCmdIndex.mpi, Gretlopt.OPT_R, "full-range", 0),
    GretlOption(GretlCmdIndex.mpi, Gretlopt.OPT_F, "send-functions", 0),
    GretlOption(GretlCmdIndex.mpi, Gretlopt.OPT_L, "local", 0),
    GretlOption(GretlCmdIndex.mpi, Gretlopt.OPT_N, "np", 2),
    GretlOption(GretlCmdIndex.mpi, Gretlopt.OPT_T, "omp-threads", 2),
    GretlOption(GretlCmdIndex.mpi, Gretlopt.OPT_Q, "quiet", 0),
    GretlOption(GretlCmdIndex.mpi, Gretlopt.OPT_V, "verbose", 1),
    GretlOption(GretlCmdIndex.mpi, Gretlopt.OPT_S, "single-rng", 0),
    GretlOption(GretlCmdIndex.mpi, Gretlopt.OPT_K, "key", 2),
    GretlOption(GretlCmdIndex.labels, Gretlopt.OPT_D, "delete", 0),
    GretlOption(GretlCmdIndex.labels, Gretlopt.OPT_F, "from-file", 2),
    GretlOption(GretlCmdIndex.labels, Gretlopt.OPT_T, "to-file", 2),
    GretlOption(GretlCmdIndex.labels, Gretlopt.OPT_A, "from-array", 2),
    GretlOption(GretlCmdIndex.labels, Gretlopt.OPT_R, "to-array", 2),
    GretlOption(GretlCmdIndex.lad, Gretlopt.OPT_N, "no-vcv", 0),
    GretlOption(GretlCmdIndex.logistic, Gretlopt.OPT_M, "ymax", 2),
    GretlOption(GretlCmdIndex.logistic, Gretlopt.OPT_R, "robust", 0),
    GretlOption(GretlCmdIndex.logistic, Gretlopt.OPT_C, "cluster", 2),
    GretlOption(GretlCmdIndex.logistic, Gretlopt.OPT_F, "fixed-effects", 0),
    GretlOption(GretlCmdIndex.logit, Gretlopt.OPT_M, "multinomial", 0),
    GretlOption(GretlCmdIndex.logit, Gretlopt.OPT_P, "p-values", 0),
    GretlOption(GretlCmdIndex.logit, Gretlopt.OPT_R, "robust", 0),
    GretlOption(GretlCmdIndex.logit, Gretlopt.OPT_C, "cluster", 2),
    GretlOption(GretlCmdIndex.logit, Gretlopt.OPT_V, "verbose", 0),
    GretlOption(GretlCmdIndex.logit, Gretlopt.OPT_S, "estrella", 0),
    GretlOption(GretlCmdIndex.logit, Gretlopt.OPT_N, "no-qr", 0),
    GretlOption(GretlCmdIndex.loop, Gretlopt.OPT_D, "decr", 0),
    GretlOption(GretlCmdIndex.loop, Gretlopt.OPT_P, "progressive", 0),
    GretlOption(GretlCmdIndex.loop, Gretlopt.OPT_V, "verbose", 0),
    GretlOption(GretlCmdIndex.mahal, Gretlopt.OPT_S, "save", 0),
    GretlOption(GretlCmdIndex.mahal, Gretlopt.OPT_V, "vcv", 0),
    GretlOption(GretlCmdIndex.meantest, Gretlopt.OPT_O, "unequal-vars", 0),
    GretlOption(GretlCmdIndex.meantest, Gretlopt.OPT_D, "split-by", 2),
    GretlOption(GretlCmdIndex.meantest, Gretlopt.OPT_P, "paired", 0),
    GretlOption(GretlCmdIndex.meantest, Gretlopt.OPT_R, "robust", 0),
    GretlOption(GretlCmdIndex.midasreg, Gretlopt.OPT_R, "robust", 0),
    GretlOption(GretlCmdIndex.midasreg, Gretlopt.OPT_V, "verbose", 0),
    GretlOption(GretlCmdIndex.midasreg, Gretlopt.OPT_L, "levenberg", 0),
    GretlOption(GretlCmdIndex.midasreg, Gretlopt.OPT_P, "print-spec", 0),
    GretlOption(GretlCmdIndex.midasreg, Gretlopt.OPT_B, "breaktest", 0),
    GretlOption(GretlCmdIndex.midasreg, Gretlopt.OPT_C, "clamp-beta", 0), /* legacy */
    GretlOption(GretlCmdIndex.mle, Gretlopt.OPT_A, "auxiliary", 0),
    GretlOption(GretlCmdIndex.mle, Gretlopt.OPT_G, "opg", 0),
    GretlOption(GretlCmdIndex.mle, Gretlopt.OPT_H, "hessian", 0),
    GretlOption(GretlCmdIndex.mle, Gretlopt.OPT_S, "no-gradient-check", 0),
    GretlOption(GretlCmdIndex.mle, Gretlopt.OPT_L, "lbfgs", 0),
    GretlOption(GretlCmdIndex.mle, Gretlopt.OPT_N, "numerical", 0),
    GretlOption(GretlCmdIndex.mle, Gretlopt.OPT_R, "robust", 1),
    GretlOption(GretlCmdIndex.mle, Gretlopt.OPT_C, "cluster", 2),
    GretlOption(GretlCmdIndex.mle, Gretlopt.OPT_V, "verbose", 0),
    GretlOption(GretlCmdIndex.modprint, Gretlopt.OPT_A, "addstats", 2),
    GretlOption(GretlCmdIndex.modprint, Gretlopt.OPT_O, "output", 2),
    GretlOption(GretlCmdIndex.modprint, Gretlopt.OPT_C, "complete", 0),
    GretlOption(GretlCmdIndex.modeltab, Gretlopt.OPT_O, "output", 2),
    GretlOption(GretlCmdIndex.modeltab, Gretlopt.OPT_C, "complete", 0),
    GretlOption(GretlCmdIndex.modeltab, Gretlopt.OPT_B, "options", 2),
    GretlOption(GretlCmdIndex.mpols, Gretlopt.OPT_S, "simple-print", 0),
    GretlOption(GretlCmdIndex.negbin, Gretlopt.OPT_G, "opg", 0),
    GretlOption(GretlCmdIndex.negbin, Gretlopt.OPT_M, "model1", 0),
    GretlOption(GretlCmdIndex.negbin, Gretlopt.OPT_R, "robust", 0),
    GretlOption(GretlCmdIndex.negbin, Gretlopt.OPT_C, "cluster", 2),
    GretlOption(GretlCmdIndex.negbin, Gretlopt.OPT_V, "verbose", 0),
    GretlOption(GretlCmdIndex.nls, Gretlopt.OPT_N, "numerical", 0),
    GretlOption(GretlCmdIndex.nls, Gretlopt.OPT_R, "robust", 0),
    GretlOption(GretlCmdIndex.nls, Gretlopt.OPT_V, "verbose", 0),
    GretlOption(GretlCmdIndex.nls, Gretlopt.OPT_S, "no-gradient-check", 0),
    GretlOption(GretlCmdIndex.normtest, Gretlopt.OPT_A, "all", 0),
    GretlOption(GretlCmdIndex.normtest, Gretlopt.OPT_D, "dhansen", 0),
    GretlOption(GretlCmdIndex.normtest, Gretlopt.OPT_W, "swilk", 0),
    GretlOption(GretlCmdIndex.normtest, Gretlopt.OPT_J, "jbera", 0),
    GretlOption(GretlCmdIndex.normtest, Gretlopt.OPT_L, "lillie", 0),
    GretlOption(GretlCmdIndex.nulldata, Gretlopt.OPT_N, "no-index", 0),
    GretlOption(GretlCmdIndex.nulldata, Gretlopt.OPT_P, "preserve", 0),
    GretlOption(GretlCmdIndex.ols, Gretlopt.OPT_F, "print-final", 0),
    GretlOption(GretlCmdIndex.ols, Gretlopt.OPT_J, "jackknife", 0),
    GretlOption(GretlCmdIndex.ols, Gretlopt.OPT_N, "no-df-corr", 0),
    GretlOption(GretlCmdIndex.ols, Gretlopt.OPT_O, "vcv", 0),
    GretlOption(GretlCmdIndex.ols, Gretlopt.OPT_R, "robust", 0),
    GretlOption(GretlCmdIndex.ols, Gretlopt.OPT_Q, "quiet", 0), /* note: for the sake of documentation */
    GretlOption(GretlCmdIndex.ols, Gretlopt.OPT_S, "simple-print", 0),
    GretlOption(GretlCmdIndex.ols, Gretlopt.OPT_V, "anova", 0),
    GretlOption(GretlCmdIndex.ols, Gretlopt.OPT_C, "cluster", 2),
    GretlOption(GretlCmdIndex.ols, Gretlopt.OPT_W, "window", 0),
    GretlOption(GretlCmdIndex.omit, Gretlopt.OPT_A, "auto", 1),
    GretlOption(GretlCmdIndex.omit, Gretlopt.OPT_B, "both", 0),
    GretlOption(GretlCmdIndex.omit, Gretlopt.OPT_X, "chi-square", 0),
    GretlOption(GretlCmdIndex.omit, Gretlopt.OPT_I, "silent", 0),
    GretlOption(GretlCmdIndex.omit, Gretlopt.OPT_W, "test-only", 0),
    GretlOption(GretlCmdIndex.omit, Gretlopt.OPT_T, "trend", 0), /* omit auto-trend: VAR only */
    GretlOption(GretlCmdIndex.omit, Gretlopt.OPT_E, "seasonals", 0), /* omit auto-seasonals: VAR only */
    GretlOption(GretlCmdIndex.open, Gretlopt.OPT_A, "all-cols", 0),
    GretlOption(GretlCmdIndex.open, Gretlopt.OPT_B, "progress-bar", 0),
    GretlOption(GretlCmdIndex.open, Gretlopt.OPT_D, "drop-empty", 0),
    GretlOption(GretlCmdIndex.open, Gretlopt.OPT_E, "select", 2),
    GretlOption(GretlCmdIndex.open, Gretlopt.OPT_F, "fixed-cols", 2),
    GretlOption(GretlCmdIndex.open, Gretlopt.OPT_O, "odbc", 0),
    GretlOption(GretlCmdIndex.open, Gretlopt.OPT_P, "preserve", 0),
    GretlOption(GretlCmdIndex.open, Gretlopt.OPT_R, "rowoffset", 2),
    GretlOption(GretlCmdIndex.open, Gretlopt.OPT_C, "coloffset", 2),
    GretlOption(GretlCmdIndex.open, Gretlopt.OPT_S, "sheet", 2),
    GretlOption(GretlCmdIndex.open, Gretlopt.OPT_W, "www", 0),
    GretlOption(GretlCmdIndex.open, Gretlopt.OPT_L, "cols", 2),
    GretlOption(GretlCmdIndex.open, Gretlopt.OPT_M, "rowmask", 2),
    GretlOption(GretlCmdIndex.open, Gretlopt.OPT_V, "verbose", 0),
    GretlOption(GretlCmdIndex.open, Gretlopt.OPT_K, "frompkg", 2),
    GretlOption(GretlCmdIndex.open, Gretlopt.OPT_H, "no-header", 0),
    GretlOption(GretlCmdIndex.open, Gretlopt.OPT_I, "ignore-quotes", 0),
    GretlOption(GretlCmdIndex.open, Gretlopt.OPT_U, "bundle", 2),
    GretlOption(GretlCmdIndex.outfile, Gretlopt.OPT_A, "append", 0),
    GretlOption(GretlCmdIndex.outfile, Gretlopt.OPT_Q, "quiet", 0),
    GretlOption(GretlCmdIndex.outfile, Gretlopt.OPT_D, "decpoint", 0),
    GretlOption(GretlCmdIndex.outfile, Gretlopt.OPT_B, "buffer", 1), /* note: 1 is just for backward compat */
    GretlOption(GretlCmdIndex.outfile, Gretlopt.OPT_T, "tempfile", 2),
    GretlOption(GretlCmdIndex.panel, Gretlopt.OPT_B, "between", 0),
    GretlOption(GretlCmdIndex.panel, Gretlopt.OPT_C, "cluster", 2),
    GretlOption(GretlCmdIndex.panel, Gretlopt.OPT_D, "time-dummies", 1),
    GretlOption(GretlCmdIndex.panel, Gretlopt.OPT_E, "nerlove", 0),
    GretlOption(GretlCmdIndex.panel, Gretlopt.OPT_F, "fixed-effects", 0),
    GretlOption(GretlCmdIndex.panel, Gretlopt.OPT_I, "iterate", 0),
    GretlOption(GretlCmdIndex.panel, Gretlopt.OPT_M, "matrix-diff", 0),
    GretlOption(GretlCmdIndex.panel, Gretlopt.OPT_N, "no-df-corr", 0),
    GretlOption(GretlCmdIndex.panel, Gretlopt.OPT_P, "pooled", 0),
    GretlOption(GretlCmdIndex.panel, Gretlopt.OPT_R, "robust", 0),
    GretlOption(GretlCmdIndex.panel, Gretlopt.OPT_U, "random-effects", 0),
    GretlOption(GretlCmdIndex.panel, Gretlopt.OPT_V, "verbose", 0),
    GretlOption(GretlCmdIndex.panel, Gretlopt.OPT_H, "unit-weights", 0),
    GretlOption(GretlCmdIndex.panel, Gretlopt.OPT_X, "unbalanced", 1),
    GretlOption(GretlCmdIndex.panplot, Gretlopt.OPT_M, "means", 0),
    GretlOption(GretlCmdIndex.panplot, Gretlopt.OPT_V, "overlay", 0),
    GretlOption(GretlCmdIndex.panplot, Gretlopt.OPT_S, "sequence", 0),
    GretlOption(GretlCmdIndex.panplot, Gretlopt.OPT_D, "grid", 0),
    GretlOption(GretlCmdIndex.panplot, Gretlopt.OPT_A, "stack", 0),
    GretlOption(GretlCmdIndex.panplot, Gretlopt.OPT_B, "boxplots", 0),
    GretlOption(GretlCmdIndex.panplot, Gretlopt.OPT_C, "boxplot", 0),
    GretlOption(GretlCmdIndex.panplot, Gretlopt.OPT_Y, "single-yaxis", 0),
    GretlOption(GretlCmdIndex.panplot, Gretlopt.OPT_W, "font", 2),
    GretlOption(GretlCmdIndex.panspec, Gretlopt.OPT_M, "matrix-diff", 0),
    GretlOption(GretlCmdIndex.panspec, Gretlopt.OPT_N, "nerlove", 0),
    GretlOption(GretlCmdIndex.poisson, Gretlopt.OPT_R, "robust", 0),
    GretlOption(GretlCmdIndex.poisson, Gretlopt.OPT_C, "cluster", 2),
    GretlOption(GretlCmdIndex.poisson, Gretlopt.OPT_V, "verbose", 0),
    GretlOption(GretlCmdIndex.pca, Gretlopt.OPT_C, "covariance", 0),
    GretlOption(GretlCmdIndex.pca, Gretlopt.OPT_A, "save-all", 0),
    GretlOption(GretlCmdIndex.pca, Gretlopt.OPT_O, "save", 1),
    GretlOption(GretlCmdIndex.pca, Gretlopt.OPT_Q, "quiet", 0),
    GretlOption(GretlCmdIndex.pergm, Gretlopt.OPT_O, "bartlett", 0),
    GretlOption(GretlCmdIndex.pergm, Gretlopt.OPT_L, "log", 0),
    GretlOption(GretlCmdIndex.pergm, Gretlopt.OPT_R, "radians", 0),
    GretlOption(GretlCmdIndex.pergm, Gretlopt.OPT_D, "degrees", 0),
    GretlOption(GretlCmdIndex.pergm, Gretlopt.OPT_S, "silent", 0),
    GretlOption(GretlCmdIndex.pkg, Gretlopt.OPT_L, "local", 0),
    GretlOption(GretlCmdIndex.pkg, Gretlopt.OPT_V, "verbose", 0),
    GretlOption(GretlCmdIndex.pkg, Gretlopt.OPT_S, "staging", 0),
    GretlOption(GretlCmdIndex.print, Gretlopt.OPT_O, "byobs", 0),
    GretlOption(GretlCmdIndex.print, Gretlopt.OPT_H, "head", 1),
    GretlOption(GretlCmdIndex.print, Gretlopt.OPT_T, "tail", 1),
    GretlOption(GretlCmdIndex.print, Gretlopt.OPT_L, "list", 0),
    GretlOption(GretlCmdIndex.print, Gretlopt.OPT_D, "no-dates", 0),
    GretlOption(GretlCmdIndex.print, Gretlopt.OPT_U, "numeric", 0),
    GretlOption(GretlCmdIndex.print, Gretlopt.OPT_M, "midas", 0),
    GretlOption(GretlCmdIndex.print, Gretlopt.OPT_C, "complex", 0),
    GretlOption(GretlCmdIndex.print, Gretlopt.OPT_Y, "tree", 0),
    GretlOption(GretlCmdIndex.print, Gretlopt.OPT_R, "range", 2),
    GretlOption(GretlCmdIndex.print, Gretlopt.OPT_X, "data-only", 0),
    GretlOption(GretlCmdIndex.probit, Gretlopt.OPT_P, "p-values", 0),
    GretlOption(GretlCmdIndex.probit, Gretlopt.OPT_R, "robust", 0),
    GretlOption(GretlCmdIndex.probit, Gretlopt.OPT_C, "cluster", 2),
    GretlOption(GretlCmdIndex.probit, Gretlopt.OPT_V, "verbose", 0),
    GretlOption(GretlCmdIndex.probit, Gretlopt.OPT_E, "random-effects", 0),
    GretlOption(GretlCmdIndex.probit, Gretlopt.OPT_G, "quadpoints", 2),
    GretlOption(GretlCmdIndex.probit, Gretlopt.OPT_B, "bootstrap", 1),
    GretlOption(GretlCmdIndex.probit, Gretlopt.OPT_S, "estrella", 0),
    GretlOption(GretlCmdIndex.qlrtest, Gretlopt.OPT_L, "limit-to", 2),
    GretlOption(GretlCmdIndex.qqplot, Gretlopt.OPT_R, "raw", 0),
    GretlOption(GretlCmdIndex.qqplot, Gretlopt.OPT_Z, "z-scores", 0),
    GretlOption(GretlCmdIndex.quantreg, Gretlopt.OPT_I, "intervals", 1),
    GretlOption(GretlCmdIndex.quantreg, Gretlopt.OPT_N, "no-df-corr", 0),
    GretlOption(GretlCmdIndex.quantreg, Gretlopt.OPT_R, "robust", 0),
    GretlOption(GretlCmdIndex.quit, Gretlopt.OPT_X, "exit", 0),
    GretlOption(GretlCmdIndex.rename, Gretlopt.OPT_C, "case", 2),
    GretlOption(GretlCmdIndex.reset, Gretlopt.OPT_C, "cubes-only", 0),
    GretlOption(GretlCmdIndex.reset, Gretlopt.OPT_R, "robust", 0),
    GretlOption(GretlCmdIndex.reset, Gretlopt.OPT_U, "squares-only", 0),
    GretlOption(GretlCmdIndex.reset, Gretlopt.OPT_I, "silent", 0),
    GretlOption(GretlCmdIndex.restrict, Gretlopt.OPT_B, "bootstrap", 1),
    GretlOption(GretlCmdIndex.restrict, Gretlopt.OPT_F, "full", 0),
    GretlOption(GretlCmdIndex.restrict, Gretlopt.OPT_J, "jitter", 0),
    GretlOption(GretlCmdIndex.restrict, Gretlopt.OPT_V, "verbose", 0),
    GretlOption(GretlCmdIndex.restrict, Gretlopt.OPT_L, "lbfgs", 0),
    GretlOption(GretlCmdIndex.restrict, Gretlopt.OPT_N, "no-scaling", 0),
    GretlOption(GretlCmdIndex.restrict, Gretlopt.OPT_S, "silent", 0),
    GretlOption(GretlCmdIndex.restrict, Gretlopt.OPT_W, "wald", 0),
    GretlOption(GretlCmdIndex.rmplot, Gretlopt.OPT_T, "trim", 0),
    GretlOption(GretlCmdIndex.runs, Gretlopt.OPT_D, "difference", 0),
    GretlOption(GretlCmdIndex.runs, Gretlopt.OPT_E, "equal", 0),
    GretlOption(GretlCmdIndex.scatters, Gretlopt.OPT_O, "with-lines", 0),
    GretlOption(GretlCmdIndex.scatters, Gretlopt.OPT_T, "time-series", 0),
    GretlOption(GretlCmdIndex.scatters, Gretlopt.OPT_U, "output", 2),
    GretlOption(GretlCmdIndex.scatters, Gretlopt.OPT_X, "matrix", 2),
    GretlOption(GretlCmdIndex.set, Gretlopt.OPT_F, "from-file", 2),
    GretlOption(GretlCmdIndex.set, Gretlopt.OPT_T, "to-file", 2),
    GretlOption(GretlCmdIndex.setinfo, Gretlopt.OPT_C, "continuous", 0),
    GretlOption(GretlCmdIndex.setinfo, Gretlopt.OPT_D, "discrete", 0),
    GretlOption(GretlCmdIndex.setinfo, Gretlopt.OPT_I, "description", 2),
    GretlOption(GretlCmdIndex.setinfo, Gretlopt.OPT_G, "graph-name", 2),
    GretlOption(GretlCmdIndex.setinfo, Gretlopt.OPT_M, "midas", 0),
    GretlOption(GretlCmdIndex.setinfo, Gretlopt.OPT_F, "coded", 0),
    GretlOption(GretlCmdIndex.setinfo, Gretlopt.OPT_N, "numeric", 0),
    GretlOption(GretlCmdIndex.setobs, Gretlopt.OPT_C, "stacked-cross-section", 0),
    GretlOption(GretlCmdIndex.setobs, Gretlopt.OPT_P, "panel-vars", 0),
    GretlOption(GretlCmdIndex.setobs, Gretlopt.OPT_R, "restructure", 0),
    GretlOption(GretlCmdIndex.setobs, Gretlopt.OPT_S, "stacked-time-series", 0),
    GretlOption(GretlCmdIndex.setobs, Gretlopt.OPT_T, "time-series", 0),
    GretlOption(GretlCmdIndex.setobs, Gretlopt.OPT_X, "cross-section", 0),
    GretlOption(GretlCmdIndex.setobs, Gretlopt.OPT_N, "special-time-series", 0),
    GretlOption(GretlCmdIndex.setobs, Gretlopt.OPT_G, "panel-groups", 0),
    GretlOption(GretlCmdIndex.setobs, Gretlopt.OPT_I, "panel-time", 0),
    GretlOption(GretlCmdIndex.smpl, Gretlopt.OPT_A, "no-all-missing", 0),
    GretlOption(GretlCmdIndex.smpl, Gretlopt.OPT_B, "preserve-panel", 0),
    GretlOption(GretlCmdIndex.smpl, Gretlopt.OPT_B, "balanced", 0), /* legacy alias */
    GretlOption(GretlCmdIndex.smpl, Gretlopt.OPT_C, "contiguous", 0),
    GretlOption(GretlCmdIndex.smpl, Gretlopt.OPT_D, "dates", 0),
    GretlOption(GretlCmdIndex.smpl, Gretlopt.OPT_F, "full", 0),
    GretlOption(GretlCmdIndex.smpl, Gretlopt.OPT_O, "dummy", 0),
    GretlOption(GretlCmdIndex.smpl, Gretlopt.OPT_M, "no-missing", 0),
    GretlOption(GretlCmdIndex.smpl, Gretlopt.OPT_N, "random", 0),
    GretlOption(GretlCmdIndex.smpl, Gretlopt.OPT_P, "replace", 0),
    GretlOption(GretlCmdIndex.smpl, Gretlopt.OPT_R, "restrict", 0),
    GretlOption(GretlCmdIndex.smpl, Gretlopt.OPT_T, "permanent", 0),
    GretlOption(GretlCmdIndex.smpl, Gretlopt.OPT_U, "unit", 0),
    GretlOption(GretlCmdIndex.smpl, Gretlopt.OPT_X, "time", 0),
    GretlOption(GretlCmdIndex.spearman, Gretlopt.OPT_V, "verbose", 0),
    GretlOption(GretlCmdIndex.square, Gretlopt.OPT_O, "cross", 0),
    GretlOption(GretlCmdIndex.stdize, Gretlopt.OPT_C, "center-only", 0),
    GretlOption(GretlCmdIndex.stdize, Gretlopt.OPT_N, "no-df-corr", 0),
    GretlOption(GretlCmdIndex.store, Gretlopt.OPT_A, "matrix", 2),
    GretlOption(GretlCmdIndex.store, Gretlopt.OPT_D, "database", 0),
    GretlOption(GretlCmdIndex.store, Gretlopt.OPT_E, "comment", 2),
    GretlOption(GretlCmdIndex.store, Gretlopt.OPT_F, "overwrite", 0),
    GretlOption(GretlCmdIndex.store, Gretlopt.OPT_G, "dat", 0),
    GretlOption(GretlCmdIndex.store, Gretlopt.OPT_I, "decimal-comma", 0),
    GretlOption(GretlCmdIndex.store, Gretlopt.OPT_J, "jmulti", 0),
    GretlOption(GretlCmdIndex.store, Gretlopt.OPT_L, "lcnames", 0),
    GretlOption(GretlCmdIndex.store, Gretlopt.OPT_M, "gnu-octave", 0),
    GretlOption(GretlCmdIndex.store, Gretlopt.OPT_N, "no-header", 0),
    GretlOption(GretlCmdIndex.store, Gretlopt.OPT_P, "preserve-strvals", 0),
    GretlOption(GretlCmdIndex.store, Gretlopt.OPT_R, "gnu-R", 0),
    GretlOption(GretlCmdIndex.store, Gretlopt.OPT_U, "functions", 0),
    GretlOption(GretlCmdIndex.store, Gretlopt.OPT_X, "omit-obs", 0),
    GretlOption(GretlCmdIndex.store, Gretlopt.OPT_Z, "gzipped", 1),
    GretlOption(GretlCmdIndex.summary, Gretlopt.OPT_B, "by", 2),
    GretlOption(GretlCmdIndex.summary, Gretlopt.OPT_M, "multi", 0), /* treat the 1-var case like list */
    GretlOption(GretlCmdIndex.summary, Gretlopt.OPT_S, "simple", 0),
    GretlOption(GretlCmdIndex.summary, Gretlopt.OPT_W, "weights", 2),
    GretlOption(GretlCmdIndex.summary, Gretlopt.OPT_X, "matrix", 2),
    GretlOption(GretlCmdIndex.system, Gretlopt.OPT_I, "iterate", 0),
    GretlOption(GretlCmdIndex.system, Gretlopt.OPT_V, "verbose", 0),
    GretlOption(GretlCmdIndex.system, Gretlopt.OPT_R, "robust", 0),
    GretlOption(GretlCmdIndex.system, Gretlopt.OPT_N, "no-df-corr", 0),
    GretlOption(GretlCmdIndex.textplot, Gretlopt.OPT_O, "one-scale", 0),
    GretlOption(GretlCmdIndex.textplot, Gretlopt.OPT_S, "time-series", 0),
    GretlOption(GretlCmdIndex.textplot, Gretlopt.OPT_T, "tall", 0),
    GretlOption(GretlCmdIndex.tobit, Gretlopt.OPT_L, "llimit", 2),
    GretlOption(GretlCmdIndex.tobit, Gretlopt.OPT_M, "rlimit", 2),
    GretlOption(GretlCmdIndex.tobit, Gretlopt.OPT_G, "opg", 0),
    GretlOption(GretlCmdIndex.tobit, Gretlopt.OPT_R, "robust", 0),
    GretlOption(GretlCmdIndex.tobit, Gretlopt.OPT_C, "cluster", 2),
    GretlOption(GretlCmdIndex.tobit, Gretlopt.OPT_V, "verbose", 0),
    GretlOption(GretlCmdIndex.var, Gretlopt.OPT_D, "seasonals", 0),
    GretlOption(GretlCmdIndex.var, Gretlopt.OPT_F, "variance-decomp", 0),
    GretlOption(GretlCmdIndex.var, Gretlopt.OPT_H, "robust-hac", 0),
    GretlOption(GretlCmdIndex.var, Gretlopt.OPT_I, "impulse-responses", 0),
    GretlOption(GretlCmdIndex.var, Gretlopt.OPT_L, "lagselect", 0),
    GretlOption(GretlCmdIndex.var, Gretlopt.OPT_N, "nc", 0),
    GretlOption(GretlCmdIndex.var, Gretlopt.OPT_R, "robust", 0),
    GretlOption(GretlCmdIndex.var, Gretlopt.OPT_T, "trend", 0),
    GretlOption(GretlCmdIndex.var, Gretlopt.OPT_S, "silent", 0),
    GretlOption(GretlCmdIndex.var, Gretlopt.OPT_M, "minlag", 2),
    GretlOption(GretlCmdIndex.varlist, Gretlopt.OPT_A, "accessors", 0),
    GretlOption(GretlCmdIndex.varlist, Gretlopt.OPT_S, "scalars", 0),
    GretlOption(GretlCmdIndex.varlist, Gretlopt.OPT_T, "type", 2),
    GretlOption(GretlCmdIndex.varlist, Gretlopt.OPT_D, "debug", 0),
    GretlOption(GretlCmdIndex.vartest, Gretlopt.OPT_D, "split-by", 2),
    GretlOption(GretlCmdIndex.vartest, Gretlopt.OPT_R, "robust", 1),
    GretlOption(GretlCmdIndex.vecm, Gretlopt.OPT_A, "crt", 0),
    GretlOption(GretlCmdIndex.vecm, Gretlopt.OPT_D, "seasonals", 0),
    GretlOption(GretlCmdIndex.vecm, Gretlopt.OPT_F, "variance-decomp", 0),
    GretlOption(GretlCmdIndex.vecm, Gretlopt.OPT_I, "impulse-responses", 1),
    GretlOption(GretlCmdIndex.vecm, Gretlopt.OPT_N, "nc", 0),
    GretlOption(GretlCmdIndex.vecm, Gretlopt.OPT_R, "rc", 0),
    GretlOption(GretlCmdIndex.vecm, Gretlopt.OPT_C, "uc", 0),
    GretlOption(GretlCmdIndex.vecm, Gretlopt.OPT_T, "ct", 0),
    GretlOption(GretlCmdIndex.vecm, Gretlopt.OPT_V, "verbose", 0),
    GretlOption(GretlCmdIndex.vecm, Gretlopt.OPT_S, "silent", 0),
    GretlOption(GretlCmdIndex.wls, Gretlopt.OPT_R, "robust", 0),
    GretlOption(GretlCmdIndex.wls, Gretlopt.OPT_C, "cluster", 2),
    GretlOption(GretlCmdIndex.wls, Gretlopt.OPT_Z, "allow-zeros", 0),
    GretlOption(GretlCmdIndex.xcorrgm, Gretlopt.OPT_S, "silent", 0), /* legacy (= quiet) */
    GretlOption(GretlCmdIndex.xtab, Gretlopt.OPT_C, "column", 0),
    GretlOption(GretlCmdIndex.xtab, Gretlopt.OPT_X, "matrix", 2),
    GretlOption(GretlCmdIndex.xtab, Gretlopt.OPT_R, "row", 0),
    GretlOption(GretlCmdIndex.xtab, Gretlopt.OPT_Z, "zeros", 0),
    GretlOption(GretlCmdIndex.xtab, Gretlopt.OPT_T, "tex", 1),
    GretlOption(GretlCmdIndex.xtab, Gretlopt.OPT_N, "no-totals", 0),
    GretlOption(GretlCmdIndex.xtab, Gretlopt.OPT_E, "equal", 0),
    GretlOption(GretlCmdIndex.xtab, Gretlopt.OPT_F, "no-fisher", 0),
    GretlOption(cast(GretlCmdIndex) 0, Gretlopt.OPT_NONE, null, 0)
];

private enum GretlCmdIndex no_ci = cast(GretlCmdIndex) 0;

private const(char) * get_longopt(GretlCmdIndex ci, Gretlopt opt)
{
    int i, got_ci = 0;

    for (i=0; gretl_opts[i].ci != no_ci; i++) {
        if (gretl_opts[i].ci == ci) {
            if (gretl_opts[i].o == opt) {
                return gretl_opts[i].longopt.ptr;
            }
            got_ci = 1;
        } else if (got_ci) {
            break;
        }
    }

    return "??";
}

int cluster_option_ok(GretlCmdIndex ci)
{
    int i, got_ci = 0;

    for (i=0; gretl_opts[i].ci != no_ci; i++) {
        if (gretl_opts[i].ci == ci) {
            if (gretl_opts[i].o == Gretlopt.OPT_C &&
                !strcmp(gretl_opts[i].longopt.ptr, "cluster")) {
                return 1;
            }
            got_ci = 1;
        } else if (got_ci) {
            break;
        }
    }

    return 0;
}

/* used in tokenize.c: detects the case where a command
   has an option of the form --matrix=foo, indicating that
   the columns of a named matrix replace the dataset
   series that the comand would otherwise require as
   arguments
*/

int matrix_data_option(GretlCmdIndex ci, Gretlopt opt)
{
    if (opt & Gretlopt.OPT_X) {
        int i, got_ci = 0;

        for (i=0; gretl_opts[i].ci != no_ci; i++) {
            if (gretl_opts[i].ci == ci) {
                if (gretl_opts[i].o == Gretlopt.OPT_X &&
                    !strcmp(gretl_opts[i].longopt.ptr, "matrix")) {
                    return 1;
                }
                got_ci = 1;
            } else if (got_ci) {
                break;
            }
        }
    }

    return 0;
}

/* this function is used in compiling the gretl reference
   manual */

char ** get_all_option_strings(int * pn)
{
    char ** optstrs;
    int i, n = 0;

    for (i=0; gretl_opts[i].ci != no_ci; i++) {
        n++;
    }

    optstrs = strings_array_new(n);

    if (optstrs != null) {
        for (i=0; i<n; i++) {
            optstrs[i] = gretl_strdup(gretl_opts[i].longopt.ptr);
            if (optstrs[i] == null) {
                strings_array_free(optstrs, n);
                optstrs = null;
                break;
            }
        }
    }

    if (optstrs != null) {
        strings_array_sort(&optstrs, &n, Gretlopt.OPT_U);
        *pn = n;
    }

    return optstrs;
}

/* used in checking command documentation */

const(char) ** get_opts_for_command(GretlCmdIndex ci, int * nopt)
{
    int i, j, got_ci, n = 0;
    const(char) ** ret = null;

    if (ci != GretlCmdIndex.ols) {
        /* widely applicable options attached to OLS */
        n += vcv_opt_ok(ci);
        n += quiet_opt_ok(ci);
        n += window_opt_ok(ci);
    }

    if (ci != GretlCmdIndex.gnuplot) {
        /* common plotting options attached to GNUPLOT */
        n += plot_output_opt_ok(ci);
        n += plot_outbuf_opt_ok(ci);
    }

    if (ci != GretlCmdIndex.corr) {
        /* auxiliary plotting options attached to CORR */
        n += cmd_plot_opt_ok(ci);
    }

    for (i=0; gretl_opts[i].ci != no_ci; i++) {
        n += (gretl_opts[i].ci == ci);
    }

    if (n == 0) {
        *nopt = 0;
        return null;
    }

    ret = cast(const(char) **) malloc(n * (*ret).sizeof);
    if (ret == null) return null;

    j = got_ci = 0;
    for (i=0; gretl_opts[i].ci != no_ci; i++) {
        if (gretl_opts[i].ci == ci) {
            ret[j++] = gretl_opts[i].longopt.ptr;
            got_ci = 1;
        } else if (got_ci) {
            break;
        }
    }

    if (ci != GretlCmdIndex.ols) {
        if (vcv_opt_ok(ci)) {
            ret[j++] = "vcv";
        }
        if (quiet_opt_ok(ci)) {
            ret[j++] = "quiet";
        }
        if (window_opt_ok(ci)) {
            ret[j++] = "window";
        }
    } else if (ci != GretlCmdIndex.gnuplot) {
        if (plot_output_opt_ok(ci)) {
            ret[j++] = "output";
        }
        if (plot_outbuf_opt_ok(ci)) {
            ret[j++] = "outbuf";
        }
    } else if (ci != GretlCmdIndex.corr) {
        if (cmd_plot_opt_ok(ci)) {
            ret[j++] = "plot";
        }
    }

    *nopt = n;

    return ret;
}

private immutable FlagMatch[] flag_matches = [
    FlagMatch(Gretlopt.OPT_A, 'a'),
    FlagMatch(Gretlopt.OPT_B, 'b'),
    FlagMatch(Gretlopt.OPT_C, 'c'),
    FlagMatch(Gretlopt.OPT_D, 'd'),
    FlagMatch(Gretlopt.OPT_E, 'e'),
    FlagMatch(Gretlopt.OPT_F, 'f'),
    FlagMatch(Gretlopt.OPT_G, 'g'),
    FlagMatch(Gretlopt.OPT_H, 'h'),
    FlagMatch(Gretlopt.OPT_I, 'i'),
    FlagMatch(Gretlopt.OPT_J, 'j'),
    FlagMatch(Gretlopt.OPT_K, 'k'),
    FlagMatch(Gretlopt.OPT_L, 'l'),
    FlagMatch(Gretlopt.OPT_M, 'm'),
    FlagMatch(Gretlopt.OPT_N, 'n'),
    FlagMatch(Gretlopt.OPT_O, 'o'),
    FlagMatch(Gretlopt.OPT_P, 'p'),
    FlagMatch(Gretlopt.OPT_Q, 'q'),
    FlagMatch(Gretlopt.OPT_R, 'r'),
    FlagMatch(Gretlopt.OPT_S, 's'),
    FlagMatch(Gretlopt.OPT_T, 't'),
    FlagMatch(Gretlopt.OPT_U, 'u'),
    FlagMatch(Gretlopt.OPT_V, 'v'),
    FlagMatch(Gretlopt.OPT_W, 'w'),
    FlagMatch(Gretlopt.OPT_X, 'x'),
    FlagMatch(Gretlopt.OPT_Y, 'y'),
    FlagMatch(Gretlopt.OPT_Z, 'z'),
    FlagMatch(Gretlopt.OPT_NONE, '\0')
];

/**
 * opt_from_flag:
 *
 * Returns: the gretl option value associated with a given
 * single-character flag.  Gives 0 if there is no associated
 * option.
 */

Gretlopt opt_from_flag(ubyte c)
{
    int i;

    for (i=0; flag_matches[i].c != '\0'; i++) {
        if (c == flag_matches[i].c) {
            return flag_matches[i].o;
        }
    }

    return Gretlopt.OPT_NONE;
}

private int opt_is_valid(Gretlopt opt, GretlCmdIndex ci, char c)
{
    int i, got_ci = 0;

    if (opt == Gretlopt.OPT_O && vcv_opt_ok(ci)) {
        return 1;
    } else if (opt == Gretlopt.OPT_Q && quiet_opt_ok(ci)) {
        return 1;
    } else if (opt == Gretlopt.OPT_W && window_opt_ok(ci)) {
        return 1;
    } else if (opt == Gretlopt.OPT_U && plot_output_opt_ok(ci)) {
        return 1;
    } else if (opt == Gretlopt.OPT_b && plot_outbuf_opt_ok(ci)) {
        return 1;
    } else if (opt == Gretlopt.OPT_U && cmd_plot_opt_ok(ci)) {
        return 1;
    }

    for (i=0; gretl_opts[i].ci != no_ci; i++) {
        if (ci == gretl_opts[i].ci) {
            if (opt == gretl_opts[i].o) {
                return 1;
            }
            got_ci = 1;
        } else if (got_ci) {
            break;
        }
    }

    if (c != 0) {
        gretl_errmsg_sprintf("Invalid option '-%c'", cast(int) c);
    }

    return 0;
}

private enum StoredOptFlags {
    setopt  = 1 << 0,
    persist = 1 << 1
}

/* The following apparatus is used for

   (a) setting and retrieving parameters associated with
       command options, as in --opt=val, and

   (b) storing options for a specified command via the
       "setopt" command (with or without parameters).
*/

private struct StoredOpt {
    GretlCmdIndex ci; /* index of the associated command */
    Gretlopt opt;     /* option flag */
    char * val;       /* option parameter value, or NULL */
    int flags;        /* may include OPT_SETOPT, OPT_PERSIST */
    int fd;           /* "function depth" at which registered */
}

private __gshared StoredOpt * optinfo;
private __gshared int n_stored_opts;

private void clear_one_option(StoredOpt * so)
{
    free(so.val);
    so.val = null;
    so.ci = no_ci;
    so.opt = Gretlopt.OPT_NONE;
    so.flags = 0;
}

/**
 * clear_stored_options_for_command:
 * @ci: target command index.
 *
 * Clears any (non-persistent) option information currently
 * associated with the given command, identified by its
 * index.
 */

void clear_stored_options_for_command(GretlCmdIndex ci)
{
    int i, fd = gretl_function_depth();

    for (i=0; i<n_stored_opts; i++) {
        if (optinfo[i].fd == fd && optinfo[i].ci == ci &&
            !(optinfo[i].flags & StoredOptFlags.persist)) {
            clear_one_option(&optinfo[i]);
        }
    }
}

/* called on exiting a user-defined function */

void destroy_option_params_at_level(int level)
{
    int i, n = n_stored_opts;

    for (i=0; i<n_stored_opts; i++) {
        if (optinfo[i].fd == level) {
            clear_one_option(&optinfo[i]);
            n--;
        }
    }

    if (n == 0) {
        free(optinfo);
        optinfo = null;
    }

    n_stored_opts = n;
}

/* unconditionally clean up the entire optinfo stack */

void stored_options_cleanup()
{
    int i;

    for (i=0; i<n_stored_opts; i++) {
        free(optinfo[i].val);
    }

    free(optinfo);
    optinfo = null;
    n_stored_opts = 0;
}

/* scrub just those options set via "setopt" */

void setopt_cleanup()
{
    int i, n = n_stored_opts;

    for (i=0; i<n_stored_opts; i++) {
        if (optinfo[i].flags & StoredOptFlags.setopt) {
            clear_one_option(&optinfo[i]);
            n--;
        }
    }

    if (n == 0) {
        free(optinfo);
        optinfo = null;
    }

    n_stored_opts = n;
}

private StoredOpt * matching_stored_opt(GretlCmdIndex ci, Gretlopt opt)
{
    int i, fd = gretl_function_depth();

    for (i=0; i<n_stored_opts; i++) {
        StoredOpt * so = &optinfo[i];

        if (so.ci == ci && so.opt == opt && so.fd == fd) {
            return so;
        }
    }

    return null;
}

private StoredOpt * empty_stored_opt_slot()
{
    int i;

    for (i=0; i<n_stored_opts; i++) {
        if (optinfo[i].ci == no_ci) {
            return &optinfo[i];
        }
    }

    return null;
}

/* for a given (@ci, @opt) pair, determine its status with regard
   to a parameter value: 0 = not allowed, 1 = optional, 2 = required
*/

private int option_parm_status(GretlCmdIndex ci, Gretlopt opt)
{
    int i, got_ci = 0;

    if (opt == Gretlopt.OPT_U) {
        if (plot_output_opt_ok(ci)) {
            ci = GretlCmdIndex.gnuplot;
        } else if (cmd_plot_opt_ok(ci)) {
            ci = GretlCmdIndex.corr;
        }
    } else if (opt == Gretlopt.OPT_b && plot_outbuf_opt_ok(ci)) {
        ci = GretlCmdIndex.gnuplot;
    }

    for (i=0; gretl_opts[i].ci != no_ci; i++) {
        if (gretl_opts[i].ci == ci) {
            if (gretl_opts[i].o == opt) {
                return gretl_opts[i].parminfo;
            }
            got_ci = 1;
        } else if (got_ci) {
            break;
        }
    }

    return 0;
}

/* Handle the case where we got input on the pattern

   --option-with-param=eval(foo)

   where "foo" is the name of a string variable and the
   option parameter should be set to the value of that
   variable.
*/

private int maybe_evaluate_optval(StoredOpt * so)
{
    char * s = so.val;
    int n = cast(int) strlen(s);
    int err = 0;

    if (!strncmp(s, "eval(", 5) && s[n-1] == ')') {
        char * sname = gretl_strndup(s + 5, n - 6);

        if (sname != null) {
            char * tmp = generate_string(sname, null, &err);

            if (tmp != null) {
                free(so.val);
                so.val = tmp;
            }
            free(sname);
        }
    }

    return err;
}

private int real_push_option(GretlCmdIndex ci, Gretlopt opt, char * val,
                             int checked, int flags)
{
    int fd = gretl_function_depth();
    int n, err = 0;
    int val_set = 0;
    StoredOpt * so;

    if (!checked && cast(int) ci > 0 &&
        option_parm_status(ci, opt) == OptStatus.no_parm) {
        return GretlError.data;
    }

    so = matching_stored_opt(ci, opt);

    if (so != null) {
        /* got a match for the (ci, opt) pair already */
        if (!(flags & StoredOptFlags.setopt)) {
            free(so.val);
            so.val = val;
            val_set = 1;
        }
        so.flags = flags;
        goto finish;
    }

    so = empty_stored_opt_slot();

    if (so != null) {
        /* re-use a vacant slot */
        so.ci = ci;
        so.opt = opt;
        so.val = val;
        val_set = 1;
        so.flags = flags;
        so.fd = fd;
        goto finish;
    }

    /* so we have to extend the array */
    n = n_stored_opts + 1;
    so = cast(StoredOpt *) realloc(optinfo, n * (*so).sizeof);

    if (so == null) {
        err = GretlError.alloc;
    } else {
        optinfo = so;
        so = &optinfo[n-1];
        so.ci = ci;
        so.opt = opt;
        so.val = val;
        val_set = 1;
        so.flags = flags;
        so.fd = fd;
        n_stored_opts = n;
    }

 finish:

    if (val_set && so.val != null) {
        maybe_evaluate_optval(so);
    }

    return err;
}

/**
 * push_option_param:
 * @ci: gretl command index.
 * @opt: gretl option flag.
 * @val: parameter value as string.
 *
 * Pushes onto an internal stack a record of the @val
 * to be associated with @opt for the current @ci. Note
 * that the command option apparatus takes ownership of
 * @val, so the value passed in should be copied if need
 * be.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int push_option_param(GretlCmdIndex ci, Gretlopt opt, char * val)
{
    return real_push_option(ci, opt, val, 0, 0);
}

/**
 * get_optval_string:
 * @ci: gretl command index.
 * @opt: gretl option value.
 *
 * Returns: the ancillary string value currently
 * associated with option @opt for command @ci, if any,
 * otherwise %NULL.
 */

const(char) * get_optval_string(GretlCmdIndex ci, Gretlopt opt)
{
    StoredOpt * so = matching_stored_opt(ci, opt);

    return (so != null)? so.val : null;
}

/**
 * get_optval_double:
 * @ci: gretl command index.
 * @opt: gretl option value.
 * @err: location to receive error code.
 *
 * Returns: the double-precision parameter value currently
 * associated with option @opt for command @ci, if any,
 * otherwise #NADBL. If @opt is an active option for
 * @ci but the parameter for this option cannot be
 * interpreted as a numerical value, E_INVARG is written
 * into @err.
 */

double get_optval_double(GretlCmdIndex ci, Gretlopt opt, int * err)
{
    StoredOpt * so = matching_stored_opt(ci, opt);
    double ret = NADBL;

    if (so != null && so.val != null) {
        ret = gretl_double_from_string(so.val, err);
        if (*err) {
            ret = generate_scalar(so.val, null, err);
        }
        if (*err) {
            gretl_errmsg_sprintf(_("%s: invalid option argument").ptr, so.val);
            *err = GretlError.invarg;
        }
    }

    return ret;
}

/**
 * get_optval_int:
 * @ci: gretl command index.
 * @opt: gretl option value.
 * @err: location to receive error code.
 *
 * Returns: the integer ancillary value currently
 * associated with option @opt for command @ci, if any,
 * otherwise 0. A non-zero value is written to @err if
 * such a value is required for the option in question
 * but is not present.
 */

int get_optval_int(GretlCmdIndex ci, Gretlopt opt, int * err)
{
    StoredOpt * so = matching_stored_opt(ci, opt);
    int status = option_parm_status(ci, opt);
    int ret = 0;

    if (so != null && so.val != null) {
        ret = gretl_int_from_string(so.val, err);
        if (*err) {
            ret = generate_int(so.val, null, err);
        }
        if (*err) {
            gretl_errmsg_sprintf(_("%s: invalid option argument").ptr, so.val);
            *err = GretlError.invarg;
        }
    } else if (status == 2 && err != null) {
        const(char) * longopt = get_longopt(ci, opt);

        gretl_errmsg_sprintf(_("The option '--%s' requires a parameter").ptr,
                             longopt);
        *err = GretlError.data;
    }

    return ret;
}

private void clear_options_for_command(GretlCmdIndex ci)
{
    int i, fd = gretl_function_depth();

    for (i=0; i<n_stored_opts; i++) {
        if (optinfo[i].fd == fd && optinfo[i].ci == ci) {
            clear_one_option(&optinfo[i]);
        }
    }
}

private void set_stored_options(GretlCmdIndex ci, Gretlopt opt, int flags)
{
    int i, got_ci = 0;

    if (vcv_opt_ok(ci) && (opt & Gretlopt.OPT_O)) {
        real_push_option(ci, Gretlopt.OPT_O, null, 1, flags);
        opt &= ~Gretlopt.OPT_O; /* handled */
    }

    if (quiet_opt_ok(ci) && (opt & Gretlopt.OPT_Q)) {
        real_push_option(ci, Gretlopt.OPT_Q, null, 1, flags);
        opt &= ~Gretlopt.OPT_Q; /* handled */
    }

    if (window_opt_ok(ci) && (opt & Gretlopt.OPT_W)) {
        real_push_option(ci, Gretlopt.OPT_W, null, 1, flags);
        opt &= ~Gretlopt.OPT_W; /* handled */
    }

    if (opt == 0) {
        return;
    }

    for (i=0; gretl_opts[i].o != 0; i++) {
        const(GretlOption) * gopt = &gretl_opts[i];

        if (ci == gopt.ci) {
            if (opt & gopt.o) {
                real_push_option(ci, gopt.o, null, 1, flags);
            }
            got_ci = cast(int) ci;
        } else if (got_ci > 0 && ci != got_ci) {
            break;
        }
    }
}

/* Apparatus for pre-selecting options for a specified command,
   using the "setopt" command.
*/

int set_options_for_command(const(char) * cmdword,
                            const(char) * param,
                            Gretlopt opt)
{
    GretlCmdIndex target_ci = cast(GretlCmdIndex) gretl_command_number(cmdword);
    int flags = StoredOptFlags.setopt;
    int clear = 0;
    int err = 0;

    if (target_ci == no_ci || target_ci == GretlCmdIndex.setopt) {
        gretl_errmsg_sprintf(_("field '%s' in command is invalid").ptr,
                             cmdword);
        return GretlError.data;
    }

    if (param != null && *param != '\0') {
        if (!strcmp(param, "persist")) {
            flags |= StoredOptFlags.persist;
        } else if (!strcmp(param, "clear")) {
            clear = 1;
        } else {
            gretl_errmsg_sprintf(_("field '%s' in command is invalid").ptr,
                                 param);
            return GretlError.data;
        }
    }

    if (clear) {
        clear_options_for_command(target_ci);
    } else if (opt == 0) {
        err = GretlError.args;
    } else {
        set_stored_options(target_ci, opt, flags);
    }

    return err;
}

void maybe_get_stored_options(GretlCmdIndex ci, Gretlopt * popt)
{
    int i, fd = gretl_function_depth();

    for (i=0; i<n_stored_opts; i++) {
        if (optinfo[i].fd == fd && optinfo[i].ci == ci &&
            (optinfo[i].flags & StoredOptFlags.setopt)) {
            *popt |= optinfo[i].opt;
        }
    }
}

/* When writing data as gdt this is called conditionally on OPT_Z
   being given, but when writing gdtb it is called unconditionally
   (since the format of the latter is a zipfile in all cases).
   The effects are then:

   gdt: no compression applied by default; if --gzipped given
   with no parameter, zlib level 1 used; if --gzipped given
   with parameter, param value respected within range 0 to 9.

   gdtb: zlib level 1 used by default, levels 0 to 9 can be
   specified via --gzipped + param.
*/

int get_compression_option(GretlCmdIndex ci)
{
    StoredOpt * so = matching_stored_opt(ci, Gretlopt.OPT_Z);
    int level = 0;

    if (so == null || so.val == null) {
        /* use zlib level 1 by default */
        level = 1;
    } else {
        level = atoi(so.val);
        if (level < 0) {
            level = 0;
        } else if (level > 9) {
            level = 9;
        }
    }

    return level;
}

/* below: called via GUI */

/**
 * set_optval_double:
 * @ci: gretl command index.
 * @opt: gretl option value.
 * @x: value to set.
 *
 * Sets a double-precision ancillary value to be associated
 * with option @opt for command @ci.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int set_optval_double(GretlCmdIndex ci, Gretlopt opt, double x)
{
    char[32] s;

    gretl_push_c_numeric_locale();
    sprintf(s.ptr, "%g", x);
    gretl_pop_c_numeric_locale();

    return real_push_option(ci, opt, gretl_strdup(s.ptr), 1, 0);
}

/**
 * set_optval_int:
 * @ci: gretl command index.
 * @opt: gretl option value.
 * @k: value to set.
 *
 * Sets a integer ancillary value to be associated
 * with option @opt for command @ci.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int set_optval_int(GretlCmdIndex ci, Gretlopt opt, int k)
{
    char * s = gretl_strdup_printf("%d", k);

    return real_push_option(ci, opt, s, 1, 0);
}

/**
 * set_optval_string:
 * @ci: gretl command index.
 * @opt: gretl option value.
 * @s: value to set.
 *
 * Sets a ancillary string value to be associated
 * with option @opt for command @ci.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int set_optval_string(GretlCmdIndex ci, Gretlopt opt, const(char) * s)
{
    if (s == null) {
        gretl_errmsg_set("set_optval_string: no string was supplied");
        return GretlError.data;
    } else {
        char * scpy = gretl_strdup(s);

        if (scpy == null) {
            return GretlError.alloc;
        } else {
            return real_push_option(ci, opt, scpy, 1, 0);
        }
    }
}

/* valid_long_opt: this is (semi-) public because we have need of
   it in the command tokenizer
*/

Gretlopt valid_long_opt(GretlCmdIndex ci, const(char) * s, OptStatus * status)
{
    Gretlopt opt = Gretlopt.OPT_NONE;
    int i, got_ci = 0;

    *status = OptStatus.no_parm;

    if (*s == '\0') {
        return Gretlopt.OPT_NONE;
    }

    /* common options without parameter */
    if (vcv_opt_ok(ci) && !strcmp(s, "vcv")) {
        return Gretlopt.OPT_O;
    }
    if (quiet_opt_ok(ci) && !strcmp(s, "quiet")) {
        return Gretlopt.OPT_Q;
    }
    if (window_opt_ok(ci) && !strcmp(s, "window")) {
        return Gretlopt.OPT_W;
    }

    /* common options with parameter: switch @ci to the
       command that "owns" the option in question
    */
    if (plot_output_opt_ok(ci) && !strcmp(s, "output")) {
        ci = GretlCmdIndex.gnuplot;
    }
    if (plot_outbuf_opt_ok(ci) && (!strcmp(s, "outbuf") ||
                                   !strcmp(s, "buffer"))) {
        ci = GretlCmdIndex.gnuplot;
    }
    if (cmd_plot_opt_ok(ci) && !strcmp(s, "plot")) {
        ci = GretlCmdIndex.corr;
    }

    /* start by looking for an exact match */
    for (i=0; gretl_opts[i].o != 0; i++) {
        if (ci == gretl_opts[i].ci) {
            if (!strcmp(s, gretl_opts[i].longopt.ptr)) {
                opt = gretl_opts[i].o;
                *status = cast(OptStatus) gretl_opts[i].parminfo;
                break;
            }
            got_ci = 1;
        } else if (got_ci) {
            break;
        }
    }

    /* if this failed, try for a unique abbreviation */
    if (opt == Gretlopt.OPT_NONE) {
        int optlen, slen = cast(int) strlen(s);
        int nmatch = 0;

        got_ci = 0;
        for (i=0; gretl_opts[i].o != 0; i++) {
            if (ci == gretl_opts[i].ci) {
                optlen = cast(int) strlen(gretl_opts[i].longopt.ptr);
                if (optlen > slen && !strncmp(s, gretl_opts[i].longopt.ptr, slen)) {
                    opt = gretl_opts[i].o;
                    *status = cast(OptStatus) gretl_opts[i].parminfo;
                    nmatch++;
                }
                got_ci = 1;
            } else if (got_ci) {
                break;
            }
        }
        if (nmatch > 1) {
            if (ci == GretlCmdIndex.setobs && !strcmp(s, "panel")) {
                /* backward compatibility: --panel-vars was there first */
                return Gretlopt.OPT_P;
            }
            *status = OptStatus.ambiguous;
            return Gretlopt.OPT_NONE;
        }
    }

    /* backward compatibility */
    if (opt == Gretlopt.OPT_NONE && !strcmp(s, "wald")) {
        opt = Gretlopt.OPT_W;
        *status = OptStatus.no_parm;
    }

    return opt;
}

/* see valid_long_opt() above */

Gretlopt valid_short_opt(GretlCmdIndex ci, char c)
{
    Gretlopt opt = Gretlopt.OPT_NONE;
    int i, ok;

    for (i=0; flag_matches[i].c != '\0'; i++) {
        if (c == flag_matches[i].c) {
            opt = flag_matches[i].o;
            break;
        }
    }

    if (opt) {
        ok = opt_is_valid(opt, ci, c);
        if (!ok) {
            opt = Gretlopt.OPT_NONE;
        }
    }

    return opt;
}

private void print_option_param(const(char) * s, Prn * prn)
{
    const(char) * qchars = "=%, ";
    const(char) * p = s;
    int wrap = 0;
    int escape = 0;

    while (*p) {
        if (strspn(p, qchars)) {
            wrap = 1;
        } else if (*p == '"') {
            escape = 1;
        }
        p++;
    }

    if (wrap) {
        if (escape) {
            pputs(prn, "=\"");
            while (*s) {
                if (*s == '"') {
                    pputs(prn, "\\\"");
                } else {
                    pputc(prn, *s);
                }
                s++;
            }
            pputs(prn, "\"\n");
        } else {
            pprintf(prn, "=\"%s\"", s);
        }
    } else {
        pprintf(prn, "=%s", s);
    }
}

private __gshared Prn * flagprn;

/**
 * print_flags:
 * @oflags: options.
 * @ci: command index, for context.
 *
 * Returns: a string representation of the options in @oflags,
 * or an empty string if no options are found.  The returned
 * value should not be modified in any way.
 */

const(char) * print_flags(Gretlopt oflags, GretlCmdIndex ci)
{
    const(char) * parm;
    Gretlopt opt, oprev;
    int i, got_ci;

    if (flagprn == null) {
        int err = 0;

        flagprn = gretl_print_new(PrnType.buffer, &err);
        if (err) {
            return "";
        }
    } else {
        gretl_print_reset_buffer(flagprn);
    }

    if (oflags == Gretlopt.OPT_NONE || ci == GretlCmdIndex.quit ||
        ci == GretlCmdIndex.genr) {
        /* no options, or only hidden ones */
        return "";
    }

    if ((oflags & Gretlopt.OPT_O) && vcv_opt_ok(ci)) {
        pputs(flagprn, " --vcv");
        oflags &= ~Gretlopt.OPT_O; /* handled */
    }
    if ((oflags & Gretlopt.OPT_Q) && quiet_opt_ok(ci)) {
        pputs(flagprn, " --quiet");
        oflags &= ~Gretlopt.OPT_Q; /* handled */
    }
    if ((oflags & Gretlopt.OPT_W) && window_opt_ok(ci)) {
        pputs(flagprn, " --window");
        oflags &= ~Gretlopt.OPT_W; /* handled */
    }

    if (plot_outbuf_opt_ok(ci)) {
        if (oflags & Gretlopt.OPT_b) {
            parm = get_optval_string(ci, Gretlopt.OPT_b);
            pprintf(flagprn, " --outbuf=%s\n", parm);
            oflags &= ~Gretlopt.OPT_b; /* handled */
        } else if (oflags & Gretlopt.OPT_U) {
            parm = get_optval_string(ci, Gretlopt.OPT_U);
            if (plot_output_opt_ok(ci)) {
                pprintf(flagprn, " --output=%s\n", parm);
            } else {
                pprintf(flagprn, " --plot=%s\n", parm);
            }
            oflags &= ~Gretlopt.OPT_U; /* handled */
        }
    }

    got_ci = 0;
    oprev = Gretlopt.OPT_NONE;
    for (i=0; gretl_opts[i].ci != no_ci; i++) {
        if (ci == gretl_opts[i].ci) {
            opt = gretl_opts[i].o;
            if ((oflags & opt) && opt != oprev) {
                pprintf(flagprn, " --%s", gretl_opts[i].longopt.ptr);
                if (gretl_opts[i].parminfo) {
                    parm = get_optval_string(ci, opt);
                    if (parm != null && *parm != '\0') {
                        print_option_param(parm, flagprn);
                    }
                }
                /* block printing an option alias */
                oprev = opt;
            }
            got_ci = 1;
        } else if (got_ci) {
            break;
        }
    }

    return gretl_print_get_buffer(flagprn);
}

void option_printing_cleanup()
{
    gretl_print_destroy(flagprn);
    flagprn = null;
}

/**
 * check_for_loop_only_options:
 * @ci: gretl command index.
 * @opt: option flag to be tested.
 * @prn: gretl printing struct.
 *
 * Returns: 1 if option @opt is applicable for command @ci only
 * in the context of a command loop (in which case a warning is
 * printed to @prn), otherwise 0.
 */

int check_for_loop_only_options(GretlCmdIndex ci, Gretlopt opt, Prn * prn)
{
    int ret = 0;

    if (ci == GretlCmdIndex.ols && (opt & Gretlopt.OPT_P)) {
        const(char) * flagstr = print_flags(Gretlopt.OPT_P, GretlCmdIndex.ols);

        pprintf(prn, _("Warning: option%s ignored outside of loop").ptr,
                flagstr);
        pputc(prn, '\n');
        ret = 1;
    }

    return ret;
}

/**
 * transcribe_options:
 * @targ: pointer to target option flags.
 * @src: source option flags.
 * @test: bitwise OR of flags that are of interest in context.
 *
 * If the intersection of the flags in @src and @test is non-
 * empty, set the corresponding flags in @targ.
 *
 * Returns: the (possibly modified) @targ.
 */

Gretlopt transcribe_option_flags(Gretlopt * targ, Gretlopt src,
                                 Gretlopt test)
{
    *targ |= (src & test);

    return *targ;
}

/**
 * delete_option_flags:
 * @targ: pointer to target option flags.
 * @test: bitwise OR of flags that to be deleted.
 *
 * If the intersection of the flags in @targ and @test is non-
 * empty, unset the corresponding flags in @targ.
 *
 * Returns: the (possibly modified) @targ.
 */

Gretlopt delete_option_flags(Gretlopt * targ, Gretlopt test)
{
    *targ &= ~(*targ & test);

    return *targ;
}

/**
 * incompatible_options:
 * @opt: option flags to be tested.
 * @test: bitwise OR of flags that are mutually incompatible in context.
 *
 * Returns: %E_BADOPT if @opt contains more than one of the flags
 * in @test, otherwise 0.
 */

int incompatible_options(Gretlopt opt, Gretlopt test)
{
    int optcount = 0;
    Gretlopt o;

    for (o=Gretlopt.OPT_A; o<=Gretlopt.OPT_i; o=cast(Gretlopt) (o<<1)) {
        if ((opt & o) && (test & o)) {
            optcount++;
            if (optcount > 1) {
                return GretlError.badopt;
            }
        }
    }

    return 0;
}

/**
 * options_incompatible_with:
 * @opt: option flags to be tested.
 * @base: "base" option.
 * @test: bitwise OR of flags that are incompatible with @base.
 *
 * Returns: %E_BADOPT if @opt contains both @base and one or more of
 * the flags in @test, otherwise 0.
 */

int options_incompatible_with(Gretlopt opt, Gretlopt base,
                              Gretlopt test)
{
    if (opt & base) {
        if (opt & test) {
            return GretlError.badopt;
        }
    }

    return 0;
}

/**
 * option_prereq_missing:
 * @opt: option flags to be tested.
 * @test: bitwise OR of flags that have a definite prequisite.
 * @prereq: bitwise OR of prerequisite flags.
 *
 * Returns: %E_BADOPT if @opt contains at least one element of
 * @test but no elements of @prereq, otherwise 0.
 */

int option_prereq_missing(Gretlopt opt, Gretlopt test,
                          Gretlopt prereq)
{
    if (opt & test) {
        if (!(opt & prereq)) {
            return GretlError.badopt;
        }
    }

    return 0;
}

/**
 * inapplicable_option_error:
 * @ci: command index.
 * @opt: bad option flag.
 *
 * Flags an error: to be used when @opt is not applicable
 * for command @ci, in context.
 *
 * Returns: %E_BADOPT.
 */

int inapplicable_option_error(GretlCmdIndex ci, Gretlopt opt)
{
    const(char) * s = print_flags(opt, ci);

    gretl_errmsg_sprintf(_("%s: inapplicable option").ptr, s);
    return GretlError.badopt;
}

void debug_print_option_flags(const(char) * msg, Gretlopt opt)
{
    if (msg != null && *msg != '\0') {
        fprintf(stderr, "%s: ", msg);
    }

    if (opt == 0) {
        fprintf(stderr, "opt=0\n");
    } else {
        char c = 'A';
        int i, started = 0;

        fprintf(stderr, "opt=%d (", cast(int) opt);

        for (i=Gretlopt.OPT_A; i<=Gretlopt.OPT_Y; i*=2) {
            if (opt & i) {
                if (started) {
                    fputc('|', stderr);
                }
                fprintf(stderr, "OPT_%c", cast(int) c);
                started = 1;
            }
            c++;
        }

        fputs(")\n", stderr);
    }
}
