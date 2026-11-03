/* Port of lib/src/gretl_commands.c: the table of command words, command-word lookup in both directions
   (index to word, word to index), the genr aliases, the help index and command completion. The header
   gretl_commands.h is gretl_commandsh.d.

   Differences from the C: the word-to-index hash table (a glib hash table) is the string-keyed table of
   ghashtable.d; the command indices are the GretlCmdIndex enum, so the functions take and return int as in the
   C and convert at the boundary where the enum members are used. */

module gretl.gretl_commands;


import core.stdc.string : strcmp, strlen, strncmp;

import gretl.ghashtable;
import gretl.glib : gpointer;
import gretl.gretl_commandsh : GretlCmdIndex;

private struct GretlCmd {
    int cnum;
    const(char) * cword;
}

private __gshared GretlCmd[] gretl_cmds = [
    GretlCmd(0, ""),
    GretlCmd(GretlCmdIndex.add, "add"),
    GretlCmd(GretlCmdIndex.adf, "adf"),
    GretlCmd(GretlCmdIndex.anova, "anova"),
    GretlCmd(GretlCmdIndex.append, "append"),
    GretlCmd(GretlCmdIndex.ar, "ar"),
    GretlCmd(GretlCmdIndex.ar1, "ar1"),
    GretlCmd(GretlCmdIndex.arch, "arch"),
    GretlCmd(GretlCmdIndex.arma, "arima"),
    GretlCmd(GretlCmdIndex.bds, "bds"),
    GretlCmd(GretlCmdIndex.biprobit, "biprobit"),
    GretlCmd(GretlCmdIndex.bkw, "bkw"),
    GretlCmd(GretlCmdIndex.break_, "break"),
    GretlCmd(GretlCmdIndex.bxplot, "boxplot"),
    GretlCmd(GretlCmdIndex.chow, "chow"),
    GretlCmd(GretlCmdIndex.clear, "clear"),
    GretlCmd(GretlCmdIndex.coeffsum, "coeffsum"),
    GretlCmd(GretlCmdIndex.coint, "coint"),
    GretlCmd(GretlCmdIndex.coint2, "johansen"),
    GretlCmd(GretlCmdIndex.continue_, "continue"),
    GretlCmd(GretlCmdIndex.corr, "corr"),
    GretlCmd(GretlCmdIndex.corrgm, "corrgm"),
    GretlCmd(GretlCmdIndex.cusum, "cusum"),
    GretlCmd(GretlCmdIndex.data, "data"),
    GretlCmd(GretlCmdIndex.datamod, "dataset"),
    GretlCmd(GretlCmdIndex.deleet, "delete"),
    GretlCmd(GretlCmdIndex.diff, "diff"),
    GretlCmd(GretlCmdIndex.difftest, "difftest"),
    GretlCmd(GretlCmdIndex.discrete, "discrete"),
    GretlCmd(GretlCmdIndex.dpanel, "dpanel"),
    GretlCmd(GretlCmdIndex.dummify, "dummify"),
    GretlCmd(GretlCmdIndex.duration, "duration"),
    GretlCmd(GretlCmdIndex.elif, "elif"),
    GretlCmd(GretlCmdIndex.else_, "else"),
    GretlCmd(GretlCmdIndex.end, "end"),
    GretlCmd(GretlCmdIndex.endif, "endif"),
    GretlCmd(GretlCmdIndex.endloop, "endloop"),
    GretlCmd(GretlCmdIndex.eqnprint, "eqnprint"),
    GretlCmd(GretlCmdIndex.equation, "equation"),
    GretlCmd(GretlCmdIndex.estimate, "estimate"),
    GretlCmd(GretlCmdIndex.eval, "eval"),
    GretlCmd(GretlCmdIndex.fcast, "fcast"),
    GretlCmd(GretlCmdIndex.flush, "flush"),
    GretlCmd(GretlCmdIndex.foreign, "foreign"),
    GretlCmd(GretlCmdIndex.fractint, "fractint"),
    GretlCmd(GretlCmdIndex.freq, "freq"),
    GretlCmd(GretlCmdIndex.func, "function"),
    GretlCmd(GretlCmdIndex.funcerr, "funcerr"),
    GretlCmd(GretlCmdIndex.garch, "garch"),
    GretlCmd(GretlCmdIndex.genr, "genr"),
    GretlCmd(GretlCmdIndex.gibbs, "gibbs"),
    GretlCmd(GretlCmdIndex.gmm, "gmm"),
    GretlCmd(GretlCmdIndex.gnuplot, "gnuplot"),
    GretlCmd(GretlCmdIndex.gpbuild, "gpbuild"),
    GretlCmd(GretlCmdIndex.graphpg, "graphpg"),
    GretlCmd(GretlCmdIndex.gridplot, "gridplot"),
    GretlCmd(GretlCmdIndex.heckit, "heckit"),
    GretlCmd(GretlCmdIndex.help, "help"),
    GretlCmd(GretlCmdIndex.hfplot, "hfplot"),
    GretlCmd(GretlCmdIndex.hsk, "hsk"),
    GretlCmd(GretlCmdIndex.hurst, "hurst"),
    GretlCmd(GretlCmdIndex.if_, "if"),
    GretlCmd(GretlCmdIndex.include, "include"),
    GretlCmd(GretlCmdIndex.info, "info"),
    GretlCmd(GretlCmdIndex.intreg, "intreg"),
    GretlCmd(GretlCmdIndex.join, "join"),
    GretlCmd(GretlCmdIndex.kdplot, "kdplot"),
    GretlCmd(GretlCmdIndex.kpss, "kpss"),
    GretlCmd(GretlCmdIndex.labels, "labels"),
    GretlCmd(GretlCmdIndex.lad, "lad"),
    GretlCmd(GretlCmdIndex.lags, "lags"),
    GretlCmd(GretlCmdIndex.ldiff, "ldiff"),
    GretlCmd(GretlCmdIndex.leverage, "leverage"),
    GretlCmd(GretlCmdIndex.levinlin, "levinlin"),
    GretlCmd(GretlCmdIndex.logistic, "logistic"),
    GretlCmd(GretlCmdIndex.logit, "logit"),
    GretlCmd(GretlCmdIndex.logs, "logs"),
    GretlCmd(GretlCmdIndex.loop, "loop"),
    GretlCmd(GretlCmdIndex.mahal, "mahal"),
    GretlCmd(GretlCmdIndex.makepkg, "makepkg"),
    GretlCmd(GretlCmdIndex.markers, "markers"),
    GretlCmd(GretlCmdIndex.meantest, "meantest"),
    GretlCmd(GretlCmdIndex.midasreg, "midasreg"),
    GretlCmd(GretlCmdIndex.mle, "mle"),
    GretlCmd(GretlCmdIndex.modeltab, "modeltab"),
    GretlCmd(GretlCmdIndex.modprint, "modprint"),
    GretlCmd(GretlCmdIndex.modtest, "modtest"),
    GretlCmd(GretlCmdIndex.mpi, "mpi"),
    GretlCmd(GretlCmdIndex.mpols, "mpols"),
    GretlCmd(GretlCmdIndex.negbin, "negbin"),
    GretlCmd(GretlCmdIndex.nls, "nls"),
    GretlCmd(GretlCmdIndex.normtest, "normtest"),
    GretlCmd(GretlCmdIndex.nulldata, "nulldata"),
    GretlCmd(GretlCmdIndex.ols, "ols"),
    GretlCmd(GretlCmdIndex.omit, "omit"),
    GretlCmd(GretlCmdIndex.open, "open"),
    GretlCmd(GretlCmdIndex.orthdev, "orthdev"),
    GretlCmd(GretlCmdIndex.outfile, "outfile"),
    GretlCmd(GretlCmdIndex.panel, "panel"),
    GretlCmd(GretlCmdIndex.panplot, "panplot"),
    GretlCmd(GretlCmdIndex.panspec, "panspec"),
    GretlCmd(GretlCmdIndex.pca, "pca"),
    GretlCmd(GretlCmdIndex.pergm, "pergm"),
    GretlCmd(GretlCmdIndex.plot, "plot"),
    GretlCmd(GretlCmdIndex.poisson, "poisson"),
    GretlCmd(GretlCmdIndex.print, "print"),
    GretlCmd(GretlCmdIndex.printf, "printf"),
    GretlCmd(GretlCmdIndex.probit, "probit"),
    GretlCmd(GretlCmdIndex.pval, "pvalue"),
    GretlCmd(GretlCmdIndex.quantreg, "quantreg"),
    GretlCmd(GretlCmdIndex.qlrtest, "qlrtest"),
    GretlCmd(GretlCmdIndex.qqplot, "qqplot"),
    GretlCmd(GretlCmdIndex.quit, "quit"),
    GretlCmd(GretlCmdIndex.rename, "rename"),
    GretlCmd(GretlCmdIndex.reset, "reset"),
    GretlCmd(GretlCmdIndex.restrict, "restrict"),
    GretlCmd(GretlCmdIndex.rmplot, "rmplot"),
    GretlCmd(GretlCmdIndex.run, "run"),
    GretlCmd(GretlCmdIndex.runs, "runs"),
    GretlCmd(GretlCmdIndex.scatters, "scatters"),
    GretlCmd(GretlCmdIndex.sdiff, "sdiff"),
    GretlCmd(GretlCmdIndex.set, "set"),
    GretlCmd(GretlCmdIndex.setinfo, "setinfo"),
    GretlCmd(GretlCmdIndex.setobs, "setobs"),
    GretlCmd(GretlCmdIndex.setopt, "setopt"),
    GretlCmd(GretlCmdIndex.setmiss, "setmiss"),
    GretlCmd(GretlCmdIndex.shell, "shell"),
    GretlCmd(GretlCmdIndex.smpl, "smpl"),
    GretlCmd(GretlCmdIndex.spearman, "spearman"),
    GretlCmd(GretlCmdIndex.square, "square"),
    GretlCmd(GretlCmdIndex.stdize, "stdize"),
    GretlCmd(GretlCmdIndex.store, "store"),
    GretlCmd(GretlCmdIndex.summary, "summary"),
    GretlCmd(GretlCmdIndex.system, "system"),
    GretlCmd(GretlCmdIndex.tabprint, "tabprint"),
    GretlCmd(GretlCmdIndex.textplot, "textplot"),
    GretlCmd(GretlCmdIndex.tobit, "tobit"),
    GretlCmd(GretlCmdIndex.ivreg, "tsls"),
    GretlCmd(GretlCmdIndex.var, "var"),
    GretlCmd(GretlCmdIndex.varlist, "varlist"),
    GretlCmd(GretlCmdIndex.vartest, "vartest"),
    GretlCmd(GretlCmdIndex.vecm, "vecm"),
    GretlCmd(GretlCmdIndex.vif, "vif"),
    GretlCmd(GretlCmdIndex.wls, "wls"),
    GretlCmd(GretlCmdIndex.xcorrgm, "xcorrgm"),
    GretlCmd(GretlCmdIndex.xtab, "xtab"),
    GretlCmd(GretlCmdIndex.funcret, "return"),
    GretlCmd(GretlCmdIndex.catch_, "catch"),
    GretlCmd(GretlCmdIndex.pkg, "pkg"),
    GretlCmd(GretlCmdIndex.nc, null)
];

private __gshared GretlCmd[] gretl_cmd_aliases = [
    GretlCmd(GretlCmdIndex.genr, "series"),
    GretlCmd(GretlCmdIndex.genr, "scalar"),
    GretlCmd(GretlCmdIndex.genr, "matrix"),
    GretlCmd(GretlCmdIndex.genr, "string"),
    GretlCmd(GretlCmdIndex.genr, "list"),
    GretlCmd(GretlCmdIndex.genr, "bundle"),
    GretlCmd(GretlCmdIndex.genr, "strings"),
    GretlCmd(GretlCmdIndex.genr, "matrices"),
    GretlCmd(GretlCmdIndex.genr, "bundles"),
    GretlCmd(GretlCmdIndex.genr, "lists"),
    GretlCmd(GretlCmdIndex.genr, "arrays"),
    GretlCmd(GretlCmdIndex.arma, "arma"),
    GretlCmd(GretlCmdIndex.coint2, "coint2"),
    GretlCmd(GretlCmdIndex.panspec, "hausman"),
    GretlCmd(GretlCmdIndex.nc, null)
];

int word_is_genr_alias(const(char) * s)
{
    int i;

    for (i=0; gretl_cmd_aliases[i].cnum == GretlCmdIndex.genr; i++) {
        if (!strcmp(s, gretl_cmd_aliases[i].cword)) {
            return 1;
        }
    }

    return 0;
}

const(char) * gretl_command_word(int i)
{
    if (i > 0 && i < GretlCmdIndex.nc) {
        return gretl_cmds[i].cword;
    } else {
        return "";
    }
}

private __gshared GHashTable * ht;

private void gretl_command_hash_init()
{
    int i;

    ht = g_hash_table_new(&g_str_hash, &g_str_equal);

    for (i=1; gretl_cmds[i].cword != null; i++) {
        g_hash_table_insert(ht, cast(gpointer) gretl_cmds[i].cword, cast(gpointer) cast(size_t) gretl_cmds[i].cnum);
    }

    for (i=0; gretl_cmd_aliases[i].cword != null; i++) {
        g_hash_table_insert(ht, cast(gpointer) gretl_cmd_aliases[i].cword, cast(gpointer) cast(size_t) gretl_cmd_aliases[i].cnum);
    }
}

int gretl_command_number(const(char) * s)
{
    gpointer p;
    int ret = 0;

    if (ht == null) {
        gretl_command_hash_init();
    }

    p = g_hash_table_lookup(ht, s);
    if (p != null) {
        ret = cast(int) cast(size_t) p;
    }

    return ret;
}

/* Note on gretl_help_index() as of 2023-10-26: the special treatment
   of "tsplots" here is kind of a nasty hack (to allow this alias for
   "scatters" to have its own help text). It will need rethinking if
   we ever add another such alias for any command.
*/

int gretl_help_index(const(char) * s)
{
    int ret = gretl_command_number(s);

    if (ret == 0 && !strcmp(s, "tsplots")) {
        /* an alias masquerading as a command */
        ret = GretlCmdIndex.nc;
    }

    return ret;
}

void gretl_command_hash_cleanup()
{
    if (ht != null) {
        g_hash_table_destroy(ht);
    }
}

const(char) * gretl_command_complete_next(const(char) * s, int idx)
{
    size_t n = strlen(s);
    int i;

    for (i=idx; i<GretlCmdIndex.nc; i++) {
        if (!strncmp(s, gretl_cmds[i].cword, n)) {
            return gretl_cmds[i].cword;
        }
    }

    return null;
}

const(char) * gretl_command_complete(const(char) * s)
{
    return gretl_command_complete_next(s, 0);
}
