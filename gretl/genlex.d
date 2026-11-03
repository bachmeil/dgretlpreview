/* Port of genlex.c: the lexer module for 'genr' and related commands: the symbol tables (constants, dummies, dataset / model / bundle accessors, built-in
   functions with their aliases and hidden ones), their look-up functions, function overrides by packages, the tokenizer `lex` and the symbol-name
   function `getsymb`.

   The C operator and function ids (F_ABS, R_NOBS, CONST_PI, DUM_DIAG, P_LISTDEF ...) are the members of the named enums GenSym, GenFunc, ConstId, DumId and
   GenFlags of genparseh.d (lower case, prefix dropped). The tables are __gshared arrays of StrTable / StrTableEx (the second carries a function pointer:
   the "pointerized" math functions are extern (C) doublecast(double) thunks here, since the library functions they stand for are D functions; the complex
   ones take a Cmplx). The function-like macros are functions. The debug output under LDEBUG (off), the code under #if 0 and the embedded-R branch
   (USE_RLIB: maybe_get_R_function is 0) are not in the port.

   Functions that other modules declared extern (C) while this file was missing (const_lookup, getsymb, mvar_lookup, constname, dvarname, mvarname,
   function_lookup, install_function_override, ...) are extern (C) and public here. */

module gretl.genlex;

import gretl.ghashtable : g_str_equal, g_str_hash;
import gretl.gstring : g_string_append, g_string_append_c, g_string_new;

import gretl.dataio : dateton;
import gretl.gretl_func1 : get_user_function_by_name, gretl_function_depth;
import gretl.genmain : current_series_index;
import gretl.geneval1 : newdbl;
import gretl.geneval5 : expr, parser_advance, parser_getc, parser_print_input, parser_ungetc;
import gretl.genfuncs : function_package_get_name, get_active_function_package, get_t_from_obs_string, gretl_round, gretl_sgn, logistic_cdf;
import gretl.genmainh : model_data_array, model_data_list, model_data_matrix, model_data_matrix_builder;
import gretl.ghashtable : g_hash_table_destroy, g_hash_table_insert, g_hash_table_lookup, g_hash_table_new, g_hash_table_new_full, g_hash_table_remove;
import gretl.glib : gchar, guint32;
import gretl.gmem : g_free;
import gretl.gretl_cmatrix : gretl_cquad;
import gretl.gretl_errors : gretl_errmsg_set, gretl_warnmsg_set;
import gretl.gretl_list : varname_match_any;
import gretl.gretl_normal : invmills;
import gretl.gretl_prnh : PrnType;
import gretl.gretl_prn : gretl_print_new, pprintf, pputc, pputs;
import gretl.gretl_utils : object_is_function_arg;
import gretl.gstrfuncs : g_strdup_printf;
import gretl.gstring : GString, g_string_free;
import gretl.pvalues : digamma, gammafun, lngamma, normal_cdf, normal_cdf_inverse, normal_pdf, trigamma;
import gretl.strutils : dot_atof, double_quote_position, gretl_charpos, gretl_strdup, gretl_strndup, is_greek_letter, shift_string_left;
import gretl.uservar : get_user_var_by_name;
import gretl.uservarh : UserVar;

import gretl.genmainh : BundleDataIndex, ModelDataIndex, RetrievalIndex;

import gretl.ghashtable : GHashTable;
import gretl.glib : gpointer;

import core.stdc.math : fabs, floor, isnan, log, sqrt, pow, exp, atan2, ceil, sin, cos, tan, asin, acos, atan, sinh, cosh, tanh, asinh, acosh, atanh, log10, log2;
import core.stdc.stdlib : free, malloc, calloc, realloc, atoi, atof, strtol, strtoul, strtod;
import core.stdc.stdio : fprintf, printf, sprintf, snprintf, sscanf, stderr;
import core.stdc.string : memcpy, memset, strcmp, strlen, strcpy, strncpy, strcat, strchr, strstr, strncmp, strspn, strcspn, strncat;
import core.stdc.ctype : isalpha, isdigit, isspace, isalnum, isupper, islower, isprint;

import gretl.genparseh;
import gretl.libgretlh;
import gretl.gretl_matrixh;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_prn : Prn;
import gretl.gretl_func0 : Ufunc, Fnpkg;
import gretl.gretl_string_table : get_built_in_string_by_name;

/* objstack.c, not ported yet */
extern (C) void * gretl_get_object_by_name(const(char) * name);

/* function pointers stored in the "pointerized" function table: extern (C) doublecast(double) thunks (the C stores the library functions themselves) */
extern (C) double ptr_acos(double x) { return acos(x); }
extern (C) double ptr_acosh(double x) { return acosh(x); }
extern (C) double ptr_asin(double x) { return asin(x); }
extern (C) double ptr_asinh(double x) { return asinh(x); }
extern (C) double ptr_atan(double x) { return atan(x); }
extern (C) double ptr_atanh(double x) { return atanh(x); }
extern (C) double ptr_cabs(Cmplx z) { return cabs(z); }
extern (C) double ptr_carg(Cmplx z) { return atan2(z.i, z.r); }
extern (C) double ptr_ceil(double x) { return ceil(x); }
extern (C) double ptr_cimag(Cmplx z) { return z.i; }
extern (C) double ptr_cos(double x) { return cos(x); }
extern (C) double ptr_cosh(double x) { return cosh(x); }
extern (C) double ptr_creal(Cmplx z) { return z.r; }
extern (C) double ptr_digamma(double x) { return digamma(x); }
extern (C) double ptr_exp(double x) { return exp(x); }
extern (C) double ptr_fabs(double x) { return fabs(x); }
extern (C) double ptr_floor(double x) { return floor(x); }
extern (C) double ptr_gammafun(double x) { return gammafun(x); }
extern (C) double ptr_gretl_cquad(Cmplx z) { return gretl_cquad(z); }
extern (C) double ptr_gretl_round(double x) { return gretl_round(x); }
extern (C) double ptr_gretl_sgn(double x) { return gretl_sgn(x); }
extern (C) double ptr_invmills(double x) { return invmills(x); }
extern (C) double ptr_lngamma(double x) { return lngamma(x); }
extern (C) double ptr_log(double x) { return log(x); }
extern (C) double ptr_log10(double x) { return log10(x); }
extern (C) double ptr_log2(double x) { return log2(x); }
extern (C) double ptr_logistic_cdf(double x) { return logistic_cdf(x); }
extern (C) double ptr_normal_cdf(double x) { return normal_cdf(x); }
extern (C) double ptr_normal_cdf_inverse(double x) { return normal_cdf_inverse(x); }
extern (C) double ptr_normal_pdf(double x) { return normal_pdf(x); }
extern (C) double ptr_sin(double x) { return sin(x); }
extern (C) double ptr_sinh(double x) { return sinh(x); }
extern (C) double ptr_sqrt(double x) { return sqrt(x); }
extern (C) double ptr_tan(double x) { return tan(x); }
extern (C) double ptr_tanh(double x) { return tanh(x); }
extern (C) double ptr_trigamma(double x) { return trigamma(x); }



/* lexer module for 'genr' and related commands */


enum NUMLEN = 32;
enum HEXLEN = 11;
enum MAXQUOTE = 64;

bool defining_list(const(Parser) * p) { return (p.flags & GenFlags.listdef) != 0; }

bool bare_data_type(int s) { return s > GenSym.punct_max && s < GenSym.dtype_max; }

bool closing_sym(int s) { return s == GenSym.g_rpr || s == GenSym.g_rbr || s == GenSym.g_rcb; }

__gshared const(char) * wordchars = ("abcdefghijklmnopqrstuvwxyz" ~
                                     "ABCDEFGHIJKLMNOPQRSTUVWXYZ" ~
                                     "0123456789_").ptr;

struct StrTable {
    int id;
    const(char) * str;
}

struct StrTableEx {
    int id;
    const(char) * str;
    void * ptr;
}

__gshared StrTable[] consts = [
    StrTable(ConstId.pi, "$pi".ptr),
    StrTable(ConstId.na, "NA".ptr),
    StrTable(ConstId.inf, "$inf".ptr),
    StrTable(ConstId.nan, "$nan".ptr),
    StrTable(ConstId.win32, "WIN32".ptr),
    StrTable(ConstId.eps, "$macheps".ptr),
    StrTable(ConstId.have_mpi, "$havempi".ptr),
    StrTable(ConstId.mpi_rank, "$mpirank".ptr),
    StrTable(ConstId.mpi_size, "$mpisize".ptr),
    StrTable(ConstId.n_proc, "$nproc".ptr),
    StrTable(ConstId.true_, "TRUE".ptr),
    StrTable(ConstId.false_, "FALSE".ptr),
    StrTable(0, null)];

__gshared StrTable[] dummies = [
    StrTable(DumId.null_, "null".ptr),
    StrTable(DumId.empty, "empty".ptr),
    StrTable(DumId.diag, "diag".ptr),
    StrTable(DumId.upper, "upper".ptr),
    StrTable(DumId.lower, "lower".ptr),
    StrTable(DumId.real_, "real".ptr),
    StrTable(DumId.imag, "imag".ptr),
    StrTable(DumId.end, "end".ptr),
    StrTable(DumId.dataset, "dataset".ptr),
    StrTable(0, null)];

/* Identify matrix-selection dummy constants:
   these can be valid only between '[' and ']'.
*/
bool MSEL_DUM(int d) { return d >= DumId.diag && d <= DumId.end; }

/* dvars: dataset- and test-related accessors */

__gshared StrTable[] dvars = [
    StrTable(RetrievalIndex.nobs, "$nobs".ptr),
    StrTable(RetrievalIndex.nvars, "$nvars".ptr),
    StrTable(RetrievalIndex.pd, "$pd".ptr),
    StrTable(RetrievalIndex.panel_pd, "$panelpd".ptr),
    StrTable(RetrievalIndex.t1, "$t1".ptr),
    StrTable(RetrievalIndex.t2, "$t2".ptr),
    StrTable(RetrievalIndex.tmax, "$tmax".ptr),
    StrTable(RetrievalIndex.datatype, "$datatype".ptr),
    StrTable(RetrievalIndex.test_stat, "$test".ptr),
    StrTable(RetrievalIndex.test_pval, "$pvalue".ptr),
    StrTable(RetrievalIndex.test_brk, "$qlrbreak".ptr),
    StrTable(RetrievalIndex.test_lnl, "$rlnl".ptr),
    StrTable(RetrievalIndex.stopwatch, "$stopwatch".ptr),
    StrTable(RetrievalIndex.time, "$time".ptr),
    StrTable(RetrievalIndex.punit, "$unit".ptr),
    StrTable(RetrievalIndex.obsmaj, "$obsmajor".ptr),
    StrTable(RetrievalIndex.obsmin, "$obsminor".ptr),
    StrTable(RetrievalIndex.obsmic, "$obsmicro".ptr),
    StrTable(RetrievalIndex.dates, "$obsdate".ptr),
    StrTable(RetrievalIndex.windows, "$windows".ptr),
    StrTable(RetrievalIndex.version_, "$version".ptr),
    StrTable(RetrievalIndex.errno, "$error".ptr),
    StrTable(RetrievalIndex.seed, "$seed".ptr),
    StrTable(RetrievalIndex.huge, "$huge".ptr),
    StrTable(RetrievalIndex.now, "$now".ptr),
    StrTable(RetrievalIndex.result, "$result".ptr),
    StrTable(RetrievalIndex.grpnames, "$grpnames".ptr),
    StrTable(RetrievalIndex.pngfont, "$pngfont".ptr),
    StrTable(RetrievalIndex.mapfile, "$mapfile".ptr),
    StrTable(RetrievalIndex.map, "$map".ptr),
    StrTable(RetrievalIndex.index, "obs".ptr),
    StrTable(RetrievalIndex.loglevel, "$loglevel".ptr),
    StrTable(RetrievalIndex.logstamp, "$logstamp".ptr),
    StrTable(0, null),];

/* mvars: model-related accessors */

__gshared StrTable[] mvars = [
    StrTable(ModelDataIndex.ess, "$ess".ptr),
    StrTable(ModelDataIndex.t, "$T".ptr),
    StrTable(ModelDataIndex.rsq, "$rsq".ptr),
    StrTable(ModelDataIndex.sigma, "$sigma".ptr),
    StrTable(ModelDataIndex.df, "$df".ptr),
    StrTable(ModelDataIndex.ncoeff, "$ncoeff".ptr),
    StrTable(ModelDataIndex.lnl, "$lnl".ptr),
    StrTable(ModelDataIndex.gmmcrit, "$gmmcrit".ptr),
    StrTable(ModelDataIndex.aic, "$aic".ptr),
    StrTable(ModelDataIndex.bic, "$bic".ptr),
    StrTable(ModelDataIndex.hqc, "$hqc".ptr),
    StrTable(ModelDataIndex.trsq, "$trsq".ptr),
    StrTable(ModelDataIndex.dw, "$dw".ptr),
    StrTable(ModelDataIndex.dwpval, "$dwpval".ptr),
    StrTable(ModelDataIndex.fstt, "$Fstat".ptr),
    StrTable(ModelDataIndex.chisq, "$chisq".ptr),
    StrTable(ModelDataIndex.diagtest, "$diagtest".ptr),
    StrTable(ModelDataIndex.diagpval, "$diagpval".ptr),
    StrTable(ModelDataIndex.pmanteau, "$pmanteau".ptr),
    StrTable(ModelDataIndex.uhat, "$uhat".ptr),
    StrTable(ModelDataIndex.yhat, "$yhat".ptr),
    StrTable(ModelDataIndex.llt, "$llt".ptr),
    StrTable(ModelDataIndex.ahat, "$ahat".ptr),
    StrTable(ModelDataIndex.sample, "$sample".ptr),
    StrTable(ModelDataIndex.h, "$h".ptr),
    StrTable(ModelDataIndex.coeff, "$coeff".ptr),
    StrTable(ModelDataIndex.se, "$stderr".ptr),
    StrTable(ModelDataIndex.vcv, "$vcv".ptr),
    StrTable(ModelDataIndex.rho, "$rho".ptr),
    StrTable(ModelDataIndex.compan, "$compan".ptr),
    StrTable(ModelDataIndex.xtxinv, "$xtxinv".ptr),
    StrTable(ModelDataIndex.vecg, "$vecGamma".ptr),
    StrTable(ModelDataIndex.vma, "$vma".ptr),
    StrTable(ModelDataIndex.fevd, "$fevd".ptr),
    StrTable(ModelDataIndex.evals, "$evals".ptr),
    StrTable(ModelDataIndex.jalpha, "$jalpha".ptr),
    StrTable(ModelDataIndex.jbeta, "$jbeta".ptr),
    StrTable(ModelDataIndex.jvbeta, "$jvbeta".ptr),
    StrTable(ModelDataIndex.js00, "$s00".ptr),
    StrTable(ModelDataIndex.js11, "$s11".ptr),
    StrTable(ModelDataIndex.js01, "$s01".ptr),
    StrTable(ModelDataIndex.ec, "$ec".ptr),
    StrTable(ModelDataIndex.hausman, "$hausman".ptr),
    StrTable(ModelDataIndex.sargan, "$sargan".ptr),
    StrTable(ModelDataIndex.sysgam, "$sysGamma".ptr),
    StrTable(ModelDataIndex.sysa, "$sysA".ptr),
    StrTable(ModelDataIndex.sysb, "$sysB".ptr),
    StrTable(ModelDataIndex.fcast, "$fcast".ptr),
    StrTable(ModelDataIndex.fcse, "$fcse".ptr),
    StrTable(ModelDataIndex.coeff_ci, "$coeff_ci".ptr),
    StrTable(ModelDataIndex.ehat, "$ehat".ptr),
    StrTable(ModelDataIndex.oddsratios, "$oddsratios".ptr),
    StrTable(ModelDataIndex.mnlprobs, "$mnlprobs".ptr), /* legacy */
    StrTable(ModelDataIndex.allprobs, "$allprobs".ptr),
    StrTable(ModelDataIndex.xlist, "$xlist".ptr),
    StrTable(ModelDataIndex.ylist, "$ylist".ptr),
    StrTable(ModelDataIndex.command, "$command".ptr),
    StrTable(ModelDataIndex.depvar, "$depvar".ptr),
    StrTable(ModelDataIndex.parnames, "$parnames".ptr),
    StrTable(0, null)];

/* bvars: bundle accessors */

__gshared StrTable[] bvars = [
    StrTable(BundleDataIndex.model, "$model".ptr),
    StrTable(BundleDataIndex.system, "$system".ptr),
    StrTable(BundleDataIndex.sysinfo, "$sysinfo".ptr),
    StrTable(0, null)];

/* Below, @ptrfuncs: table of functions for which we wish to
   attach function-pointers to the relevant Node. Nota bene:
   it's crucial that no function in @ptrfuncs is also listed
   in @funcs below!

   The order of function symbols in @ptrfuncs need not match
   the order in which they're listed in genparse.h, but it's
   crucial that every function with ID number <= FP_MAX has
   an entry in @ptrfuncs.

   "Crucial" -> certain crash on calling wrongly classified
   function!
*/

__gshared StrTableEx[] ptrfuncs = [
    StrTableEx(GenFunc.f_abs, "abs".ptr, cast(void *) &ptr_fabs),
    StrTableEx(GenFunc.f_ceil, "ceil".ptr, cast(void *) &ptr_ceil),
    StrTableEx(GenFunc.f_floor, "floor".ptr, cast(void *) &ptr_floor),
    StrTableEx(GenFunc.f_sin, "sin".ptr, cast(void *) &ptr_sin),
    StrTableEx(GenFunc.f_cos, "cos".ptr, cast(void *) &ptr_cos),
    StrTableEx(GenFunc.f_tan, "tan".ptr, cast(void *) &ptr_tan),
    StrTableEx(GenFunc.f_asin, "asin".ptr, cast(void *) &ptr_asin),
    StrTableEx(GenFunc.f_acos, "acos".ptr, cast(void *) &ptr_acos),
    StrTableEx(GenFunc.f_atan, "atan".ptr, cast(void *) &ptr_atan),
    StrTableEx(GenFunc.f_sinh, "sinh".ptr, cast(void *) &ptr_sinh),
    StrTableEx(GenFunc.f_cosh, "cosh".ptr, cast(void *) &ptr_cosh),
    StrTableEx(GenFunc.f_tanh, "tanh".ptr, cast(void *) &ptr_tanh),
    StrTableEx(GenFunc.f_asinh, "asinh".ptr, cast(void *) &ptr_asinh),
    StrTableEx(GenFunc.f_acosh, "acosh".ptr, cast(void *) &ptr_acosh),
    StrTableEx(GenFunc.f_atanh, "atanh".ptr, cast(void *) &ptr_atanh),
    StrTableEx(GenFunc.f_log, "log".ptr, cast(void *) &ptr_log),
    StrTableEx(GenFunc.f_log10, "log10".ptr, cast(void *) &ptr_log10),
    StrTableEx(GenFunc.f_log2, "log2".ptr, cast(void *) &ptr_log2),
    StrTableEx(GenFunc.f_exp, "exp".ptr, cast(void *) &ptr_exp),
    StrTableEx(GenFunc.f_sqrt, "sqrt".ptr, cast(void *) &ptr_sqrt),
    StrTableEx(GenFunc.f_gamma, "gammafun".ptr, cast(void *) &ptr_gammafun),
    StrTableEx(GenFunc.f_lngamma, "lngamma".ptr, cast(void *) &ptr_lngamma),
    StrTableEx(GenFunc.f_digamma, "digamma".ptr, cast(void *) &ptr_digamma),
    StrTableEx(GenFunc.f_trigamma, "trigamma".ptr, cast(void *) &ptr_trigamma),
    StrTableEx(GenFunc.f_invmills, "invmills".ptr, cast(void *) &ptr_invmills),
    StrTableEx(GenFunc.f_round, "round".ptr, cast(void *) &ptr_gretl_round),
    StrTableEx(GenFunc.f_sgn, "sgn".ptr, cast(void *) &ptr_gretl_sgn),
    StrTableEx(GenFunc.f_cnorm, "cnorm".ptr, cast(void *) &ptr_normal_cdf),
    StrTableEx(GenFunc.f_dnorm, "dnorm".ptr, cast(void *) &ptr_normal_pdf),
    StrTableEx(GenFunc.f_qnorm, "qnorm".ptr, cast(void *) &ptr_normal_cdf_inverse),
    StrTableEx(GenFunc.f_logistic, "logistic".ptr, cast(void *) &ptr_logistic_cdf),
    StrTableEx(GenFunc.f_real, "Re".ptr, cast(void *) &ptr_creal),
    StrTableEx(GenFunc.f_imag, "Im".ptr, cast(void *) &ptr_cimag),
    StrTableEx(GenFunc.f_carg, "carg".ptr, cast(void *) &ptr_carg),
    StrTableEx(GenFunc.f_cmod, "cmod".ptr, cast(void *) &ptr_cabs),
    StrTableEx(GenFunc.f_cquad, "cquad".ptr, cast(void *) &ptr_gretl_cquad),
    StrTableEx(0, null, null)];

__gshared StrTable[] funcs = [
    StrTable(GenFunc.f_atan2, "atan2".ptr),
    StrTable(GenFunc.f_diff, "diff".ptr),
    StrTable(GenFunc.f_ldiff, "ldiff".ptr),
    StrTable(GenFunc.f_sdiff, "sdiff".ptr),
    StrTable(GenFunc.f_llag, "lags".ptr),
    StrTable(GenFunc.f_hflag, "hflags".ptr),
    StrTable(GenFunc.f_dropcoll, "dropcoll".ptr),
    StrTable(GenFunc.f_toint, "int".ptr),
    StrTable(GenFunc.f_sort, "sort".ptr),
    StrTable(GenFunc.f_dsort, "dsort".ptr),
    StrTable(GenFunc.f_sortby, "sortby".ptr),
    StrTable(GenFunc.f_ranking, "ranking".ptr),
    StrTable(GenFunc.f_odev, "orthdev".ptr),
    StrTable(GenFunc.f_nobs, "nobs".ptr),
    StrTable(GenFunc.f_t1, "firstobs".ptr),
    StrTable(GenFunc.f_t2, "lastobs".ptr),
    StrTable(GenFunc.f_runiform, "uniform".ptr),
    StrTable(GenFunc.f_rnormal, "normal".ptr),
    StrTable(GenFunc.f_cum, "cum".ptr),
    StrTable(GenFunc.f_missing, "missing".ptr),
    StrTable(GenFunc.f_dataok, "ok".ptr),        /* opposite of missing */
    StrTable(GenFunc.f_misszero, "misszero".ptr),
    StrTable(GenFunc.f_lrvar, "lrvar".ptr),
    StrTable(GenFunc.f_lrcovar, "lrcovar".ptr),
    StrTable(GenFunc.f_fevd, "fevd".ptr),
    StrTable(GenFunc.f_quantile, "quantile".ptr),
    StrTable(GenFunc.f_median, "median".ptr),
    StrTable(GenFunc.f_gini, "gini".ptr),
    StrTable(GenFunc.f_zeromiss, "zeromiss".ptr),
    StrTable(GenFunc.f_sum, "sum".ptr),
    StrTable(GenFunc.f_sumall, "sumall".ptr),
    StrTable(GenFunc.f_mean, "mean".ptr),
    StrTable(GenFunc.f_min, "min".ptr),
    StrTable(GenFunc.f_max, "max".ptr),
    StrTable(GenFunc.f_sd, "sd".ptr),
    StrTable(GenFunc.f_vce, "var".ptr),
    StrTable(GenFunc.f_skewness, "skewness".ptr),
    StrTable(GenFunc.f_kurtosis, "kurtosis".ptr),
    StrTable(GenFunc.f_sst, "sst".ptr),
    StrTable(GenFunc.f_resample, "resample".ptr),
    StrTable(GenFunc.f_pnobs, "pnobs".ptr),     /* per-unit nobs in panels */
    StrTable(GenFunc.f_pmin, "pmin".ptr),      /* panel min */
    StrTable(GenFunc.f_pmax, "pmax".ptr),      /* panel max */
    StrTable(GenFunc.f_psum, "psum".ptr),      /* panel sum */
    StrTable(GenFunc.f_pmean, "pmean".ptr),     /* panel mean */
    StrTable(GenFunc.f_pxsum, "pxsum".ptr),     /* panel x-sectional sum */
    StrTable(GenFunc.f_pxnobs, "pxnobs".ptr),    /* panel x-sectional obs count */
    StrTable(GenFunc.f_pxmean, "pxmean".ptr),    /* panel x-sectional mean */
    StrTable(GenFunc.f_psd, "psd".ptr),       /* panel std dev */
    StrTable(GenFunc.f_pshrink, "pshrink".ptr),
    StrTable(GenFunc.f_pexpand, "pexpand".ptr),
    StrTable(GenFunc.f_hpfilt, "hpfilt".ptr),    /* Hodrick-Prescott filter */
    StrTable(GenFunc.f_bkfilt, "bkfilt".ptr),    /* Baxter-King filter */
    StrTable(GenFunc.f_bwfilt, "bwfilt".ptr),    /* Butterworth filter */
    StrTable(GenFunc.f_fracdiff, "fracdiff".ptr),  /* fractional difference */
    StrTable(GenFunc.f_boxcox, "boxcox".ptr),    /* Box-Cox transformation */
    StrTable(GenFunc.f_cov, "cov".ptr),
    StrTable(GenFunc.f_cor, "corr".ptr),
    StrTable(GenFunc.f_movavg, "movavg".ptr),
    StrTable(GenFunc.f_imat, "I".ptr),
    StrTable(GenFunc.f_zeros, "zeros".ptr),
    StrTable(GenFunc.f_ones, "ones".ptr),
    StrTable(GenFunc.f_seq, "seq".ptr),
    StrTable(GenFunc.f_replace, "replace".ptr),
    StrTable(GenFunc.f_munif, "muniform".ptr),
    StrTable(GenFunc.f_mnorm, "mnormal".ptr),
    StrTable(GenFunc.f_sumc, "sumc".ptr),
    StrTable(GenFunc.f_sumr, "sumr".ptr),
    StrTable(GenFunc.f_prodc, "prodc".ptr),
    StrTable(GenFunc.f_prodr, "prodr".ptr),
    StrTable(GenFunc.f_meanc, "meanc".ptr),
    StrTable(GenFunc.f_meanr, "meanr".ptr),
    StrTable(GenFunc.f_asort, "asort".ptr),
    StrTable(GenFunc.f_corresp, "corresp".ptr),
    StrTable(GenFunc.f_sdc, "sdc".ptr),
    StrTable(GenFunc.f_minc, "minc".ptr),
    StrTable(GenFunc.f_maxc, "maxc".ptr),
    StrTable(GenFunc.f_minr, "minr".ptr),
    StrTable(GenFunc.f_maxr, "maxr".ptr),
    StrTable(GenFunc.f_iminc, "iminc".ptr),
    StrTable(GenFunc.f_imaxc, "imaxc".ptr),
    StrTable(GenFunc.f_iminr, "iminr".ptr),
    StrTable(GenFunc.f_imaxr, "imaxr".ptr),
    StrTable(GenFunc.f_fft, "fft".ptr),
    StrTable(GenFunc.f_ffti, "ffti".ptr),
    StrTable(GenFunc.f_hdprod, "hdprod".ptr),
    StrTable(GenFunc.f_mcov, "mcov".ptr),
    StrTable(GenFunc.f_mcorr, "mcorr".ptr),
    StrTable(GenFunc.f_mxtab, "mxtab".ptr),
    StrTable(GenFunc.f_cdemean, "cdemean".ptr),
    StrTable(GenFunc.f_chol, "cholesky".ptr),
    StrTable(GenFunc.f_psdroot, "psdroot".ptr),
    StrTable(GenFunc.f_instrings, "instrings".ptr),
    StrTable(GenFunc.f_inv, "inv".ptr),
    StrTable(GenFunc.f_invpd, "invpd".ptr),
    StrTable(GenFunc.f_ginv, "ginv".ptr),
    StrTable(GenFunc.f_diag, "diag".ptr),
    StrTable(GenFunc.f_transp, "transp".ptr),
    StrTable(GenFunc.f_ctrans, "ctrans".ptr),
    StrTable(GenFunc.f_vec, "vec".ptr),
    StrTable(GenFunc.f_vech, "vech".ptr),
    StrTable(GenFunc.f_unvech, "unvech".ptr),
    StrTable(GenFunc.f_upper, "upper".ptr),
    StrTable(GenFunc.f_lower, "lower".ptr),
    StrTable(GenFunc.f_rows, "rows".ptr),
    StrTable(GenFunc.f_cols, "cols".ptr),
    StrTable(GenFunc.f_det, "det".ptr),
    StrTable(GenFunc.f_ldet, "ldet".ptr),
    StrTable(GenFunc.f_trace, "tr".ptr),
    StrTable(GenFunc.f_norm1, "onenorm".ptr),
    StrTable(GenFunc.f_infnorm, "infnorm".ptr),
    StrTable(GenFunc.f_rcond, "rcond".ptr),
    StrTable(GenFunc.f_rank, "rank".ptr),
    StrTable(GenFunc.f_qform, "qform".ptr),
    StrTable(GenFunc.f_mlag, "mlag".ptr),
    StrTable(GenFunc.f_qr, "qrdecomp".ptr),
    StrTable(GenFunc.f_eigsym, "eigensym".ptr),
    StrTable(GenFunc.f_eigen, "eigen".ptr),
    StrTable(GenFunc.f_cmult, "cmult".ptr),    /* legacy */
    StrTable(GenFunc.f_cdiv, "cdiv".ptr),     /* legacy */
    StrTable(GenFunc.f_schur, "schur".ptr),
    StrTable(GenFunc.f_eigsolve, "eigsolve".ptr),
    StrTable(GenFunc.f_nullspc, "nullspace".ptr),
    StrTable(GenFunc.f_princomp, "princomp".ptr),
    StrTable(GenFunc.f_mexp, "mexp".ptr),
    StrTable(GenFunc.f_mlog, "mlog".ptr),
    StrTable(GenFunc.f_fdjac, "fdjac".ptr),
    StrTable(GenFunc.f_bfgsmax, "BFGSmax".ptr),
    StrTable(GenFunc.f_bfgscmax, "BFGScmax".ptr),
    StrTable(GenFunc.f_nrmax, "NRmax".ptr),
    StrTable(GenFunc.f_numhess, "numhess".ptr),
    StrTable(GenFunc.f_obsnum, "obsnum".ptr),
    StrTable(GenFunc.f_isdiscr, "isdiscrete".ptr),
    StrTable(GenFunc.f_isdummy, "isdummy".ptr),
    StrTable(GenFunc.f_typeof, "typeof".ptr),
    StrTable(GenFunc.f_typename, "typename".ptr),
    StrTable(GenFunc.f_exists, "exists".ptr),
    StrTable(GenFunc.f_nelem, "nelem".ptr),
    StrTable(GenFunc.f_pdf, "pdf".ptr),
    StrTable(GenFunc.f_cdf, "cdf".ptr),
    StrTable(GenFunc.f_invcdf, "invcdf".ptr),
    StrTable(GenFunc.f_pval, "pvalue".ptr),
    StrTable(GenFunc.f_crit, "critical".ptr),
    StrTable(GenFunc.f_randgen, "randgen".ptr),
    StrTable(GenFunc.f_mrandgen, "mrandgen".ptr),
    StrTable(GenFunc.f_randgen1, "randgen1".ptr),
    StrTable(GenFunc.f_urcpval, "urcpval".ptr),
    StrTable(GenFunc.f_qlrpval, "qlrpval".ptr),
    StrTable(GenFunc.f_values, "values".ptr),
    StrTable(GenFunc.f_uniq, "uniq".ptr),
    StrTable(GenFunc.f_mshape, "mshape".ptr),
    StrTable(GenFunc.f_svd, "svd".ptr),
    StrTable(GenFunc.f_mols, "mols".ptr),
    StrTable(GenFunc.f_mpols, "mpols".ptr),
    StrTable(GenFunc.f_mrls, "mrls".ptr),
    StrTable(GenFunc.f_fols, "fols".ptr),
    StrTable(GenFunc.f_mread, "mread".ptr),
    StrTable(GenFunc.f_mwrite, "mwrite".ptr),
    StrTable(GenFunc.f_bread, "bread".ptr),
    StrTable(GenFunc.f_bwrite, "bwrite".ptr),
    StrTable(GenFunc.f_mcsel, "selifc".ptr),
    StrTable(GenFunc.f_mrsel, "selifr".ptr),
    StrTable(GenFunc.f_polroots, "polroots".ptr),
    StrTable(GenFunc.f_dumify, "dummify".ptr),
    StrTable(GenFunc.f_wmean, "wmean".ptr),
    StrTable(GenFunc.f_wvar, "wvar".ptr),
    StrTable(GenFunc.f_wsd, "wsd".ptr),
    StrTable(GenFunc.f_square, "square".ptr),
    StrTable(GenFunc.f_filter, "filter".ptr),
    StrTable(GenFunc.f_kfilter, "kfilter".ptr),
    StrTable(GenFunc.f_ksmooth, "ksmooth".ptr),
    StrTable(GenFunc.f_kdsmooth, "kdsmooth".ptr),
    StrTable(GenFunc.f_ksimul, "ksimul".ptr),
    StrTable(GenFunc.f_ksimdata, "ksimdata".ptr),
    StrTable(GenFunc.f_trimr, "trimr".ptr),
    StrTable(GenFunc.f_getenv, "getenv".ptr),
    StrTable(GenFunc.f_ngetenv, "ngetenv".ptr),
    StrTable(GenFunc.f_argname, "argname".ptr),
    StrTable(GenFunc.f_obslabel, "obslabel".ptr),
    StrTable(GenFunc.f_readfile, "readfile".ptr),
    StrTable(GenFunc.f_backtick, "grab".ptr),
    StrTable(GenFunc.f_strstr, "strstr".ptr),
    StrTable(GenFunc.f_instring, "instring".ptr),
    StrTable(GenFunc.f_strstrip, "strstrip".ptr),
    StrTable(GenFunc.f_strncmp, "strncmp".ptr),
    StrTable(GenFunc.f_strlen, "strlen".ptr),
    StrTable(GenFunc.f_printf, "printf".ptr),
    StrTable(GenFunc.f_sprintf, "sprintf".ptr),
    StrTable(GenFunc.f_sscanf, "sscanf".ptr),
    StrTable(GenFunc.f_varname, "varname".ptr),
    StrTable(GenFunc.f_varnames, "varnames".ptr),
    StrTable(GenFunc.f_varnum, "varnum".ptr),
    StrTable(GenFunc.f_tolower, "tolower".ptr),
    StrTable(GenFunc.f_toupper, "toupper".ptr),
    StrTable(GenFunc.f_cnameset, "cnameset".ptr),
    StrTable(GenFunc.f_rnameset, "rnameset".ptr),
    StrTable(GenFunc.f_ljungbox, "ljungbox".ptr),
    StrTable(GenFunc.f_msortby, "msortby".ptr),
    StrTable(GenFunc.f_lincomb, "lincomb".ptr),
    StrTable(GenFunc.f_imhof, "imhof".ptr),
    StrTable(GenFunc.f_toepsolv, "toepsolv".ptr),
    StrTable(GenFunc.f_tsolvepd, "tsolvepd".ptr),
    StrTable(GenFunc.f_kmeans, "kmeans".ptr),
    StrTable(GenFunc.f_rgbmix, "rgbmix".ptr),
    StrTable(GenFunc.f_dsum, "diagcat".ptr),
    StrTable(GenFunc.f_corrgm, "corrgm".ptr),
    StrTable(GenFunc.f_mcovg, "mcovg".ptr),
    StrTable(GenFunc.f_fcstats, "fcstats".ptr),
    StrTable(GenFunc.f_bessel, "bessel".ptr),
    StrTable(GenFunc.f_fraclag, "fraclag".ptr),
    StrTable(GenFunc.f_mrev, "mreverse".ptr),
    StrTable(GenFunc.f_deseas, "deseas".ptr),
    StrTable(GenFunc.f_tramolin, "linearize".ptr),
    StrTable(GenFunc.f_pergm, "pergm".ptr),
    StrTable(GenFunc.f_irr, "irr".ptr),
    StrTable(GenFunc.f_npv, "npv".ptr),
    StrTable(GenFunc.f_weekday, "weekday".ptr),
    StrTable(GenFunc.f_kdensity, "kdensity".ptr),
    StrTable(GenFunc.f_monthlen, "monthlen".ptr),
    StrTable(GenFunc.f_epochday, "epochday".ptr),
    StrTable(GenFunc.f_setnote, "setnote".ptr),
    StrTable(GenFunc.f_polyfit, "polyfit".ptr),
    StrTable(GenFunc.f_chowlin, "chowlin".ptr),
    StrTable(GenFunc.f_varsimul, "varsimul".ptr),
    StrTable(GenFunc.f_strsplit, "strsplit".ptr),
    StrTable(GenFunc.f_inlist, "inlist".ptr),
    StrTable(GenFunc.f_errmsg, "errmsg".ptr),
    StrTable(GenFunc.f_isconst, "isconst".ptr),
    StrTable(GenFunc.f_irf, "irf".ptr),
    StrTable(GenFunc.f_inbundle, "inbundle".ptr),
    StrTable(GenFunc.f_strsub, "strsub".ptr),
    StrTable(GenFunc.f_regsub, "regsub".ptr),
    StrTable(GenFunc.f_cnameget, "cnameget".ptr),
    StrTable(GenFunc.f_rnameget, "rnameget".ptr),
    StrTable(GenFunc.f_randint, "randint".ptr),
    StrTable(GenFunc.f_randstr, "randstr".ptr),
    StrTable(GenFunc.f_nadarwat, "nadarwat".ptr),   /* Nadaraya-Watson */
    StrTable(GenFunc.f_simann, "simann".ptr),     /* simulated annealing */
    StrTable(GenFunc.f_loess, "loess".ptr),
    StrTable(GenFunc.f_ghk, "ghk".ptr),
    StrTable(GenFunc.f_halton, "halton".ptr),
    StrTable(GenFunc.f_iwishart, "iwishart".ptr),
    StrTable(GenFunc.f_isnan, "isnan".ptr),
    StrTable(GenFunc.f_typestr, "typestr".ptr),
    StrTable(GenFunc.f_quadtab, "quadtable".ptr),
    StrTable(GenFunc.f_aggrby, "aggregate".ptr),
    StrTable(GenFunc.f_remove, "remove".ptr),
    StrTable(GenFunc.f_isodate, "isodate".ptr),
    StrTable(GenFunc.f_isoweek, "isoweek".ptr),
    StrTable(GenFunc.f_juldate, "juldate".ptr),
    StrTable(GenFunc.f_getline, "getline".ptr),
    StrTable(GenFunc.f_atof, "atof".ptr),
    StrTable(GenFunc.f_fixname, "fixname".ptr),
    StrTable(GenFunc.f_isoconv, "isoconv".ptr),
    StrTable(GenFunc.f_substr, "substr".ptr),
    StrTable(GenFunc.f_mpi_send, "mpisend".ptr),
    StrTable(GenFunc.f_mpi_recv, "mpirecv".ptr),
    StrTable(GenFunc.f_bcast, "mpibcast".ptr),
    StrTable(GenFunc.f_reduce, "mpireduce".ptr),
    StrTable(GenFunc.f_allreduce, "mpiallred".ptr),
    StrTable(GenFunc.f_scatter, "mpiscatter".ptr),
    StrTable(GenFunc.f_barrier, "mpibarrier".ptr),
    StrTable(GenFunc.f_easter, "easterday".ptr),
    StrTable(GenFunc.f_genseries, "genseries".ptr),
    StrTable(GenFunc.f_curl, "curl".ptr),
    StrTable(GenFunc.f_jsonget, "jsonget".ptr),
    StrTable(GenFunc.f_jsongetb, "jsongetb".ptr),
    StrTable(GenFunc.f_xmlget, "xmlget".ptr),
    StrTable(GenFunc.f_nlines, "nlines".ptr),
    StrTable(GenFunc.f_kpsscrit, "kpsscrit".ptr),
    StrTable(GenFunc.f_array, "array".ptr),
    StrTable(GenFunc.f_strvals, "strvals".ptr),
    StrTable(GenFunc.f_stringify, "stringify".ptr),
    StrTable(GenFunc.f_strvsort, "strvsort".ptr),
    StrTable(GenFunc.f_bootci, "bootci".ptr),
    StrTable(GenFunc.f_bootpval, "bootpval".ptr),
    StrTable(GenFunc.f_seasonals, "seasonals".ptr),
    StrTable(GenFunc.f_defarray, "defarray".ptr),
    StrTable(GenFunc.f_defbundle, "defbundle".ptr),
    StrTable(GenFunc.f_deflist, "deflist".ptr),
    StrTable(GenFunc.f_defargs, "_".ptr),
    StrTable(GenFunc.f_ksetup, "ksetup".ptr),
    StrTable(GenFunc.f_mweights, "mweights".ptr),
    StrTable(GenFunc.f_mgradient, "mgradient".ptr),
    StrTable(GenFunc.f_mlincomb, "mlincomb".ptr),
    StrTable(GenFunc.f_midasmult, "midasmult".ptr),
    StrTable(GenFunc.f_hfdiff, "hfdiff".ptr),
    StrTable(GenFunc.f_hfldiff, "hfldiff".ptr),
    StrTable(GenFunc.f_hflist, "hflist".ptr),
    StrTable(GenFunc.f_nmmax, "NMmax".ptr),
    StrTable(GenFunc.f_gssmax, "GSSmax".ptr),
    StrTable(GenFunc.f_cnumber, "cnumber".ptr),
    StrTable(GenFunc.f_naalen, "naalen".ptr),
    StrTable(GenFunc.f_kmeier, "kmeier".ptr),
    StrTable(GenFunc.f_normtest, "normtest".ptr),
    StrTable(GenFunc.f_ecdf, "ecdf".ptr),
    StrTable(GenFunc.f_npcorr, "npcorr".ptr),
    StrTable(GenFunc.f_dayspan, "dayspan".ptr),
    StrTable(GenFunc.f_smplspan, "smplspan".ptr),
    StrTable(GenFunc.f_sleep, "sleep".ptr),
    StrTable(GenFunc.f_getinfo, "getinfo".ptr),
    StrTable(GenFunc.f_cdumify, "cdummify".ptr),
    StrTable(GenFunc.f_svm, "svm".ptr),
    StrTable(GenFunc.f_getkeys, "getkeys".ptr),
    StrTable(GenFunc.f_feval, "feval".ptr),
    StrTable(GenFunc.f_fevalb, "fevalb".ptr),
    StrTable(GenFunc.f_binperms, "binperms".ptr),
    StrTable(GenFunc.f_brename, "brename".ptr),
    StrTable(GenFunc.f_ccode, "isocountry".ptr),
    StrTable(GenFunc.f_lsolve, "Lsolve".ptr),
    StrTable(GenFunc.f_hyp2f1, "hyp2f1".ptr),
    StrTable(GenFunc.f_strftime, "strftime".ptr),
    StrTable(GenFunc.f_strfday, "strfday".ptr),
    StrTable(GenFunc.f_strptime, "strptime".ptr),
    StrTable(GenFunc.f_strpday, "strpday".ptr),
    StrTable(GenFunc.f_bkw, "bkw".ptr),
    StrTable(GenFunc.f_fzero, "fzero".ptr),
    StrTable(GenFunc.f_conv2d, "conv2d".ptr),
    StrTable(GenFunc.f_msplitby, "msplitby".ptr),
    StrTable(GenFunc.f_flatten, "flatten".ptr),
    StrTable(GenFunc.f_errorif, "errorif".ptr),
    StrTable(GenFunc.f_iscmplx, "iscomplex".ptr),
    StrTable(GenFunc.f_complex, "complex".ptr),
    StrTable(GenFunc.f_conj, "conj".ptr),
    StrTable(GenFunc.f_cswitch, "cswitch".ptr),
    StrTable(GenFunc.f_randperm, "randperm".ptr),
    StrTable(GenFunc.f_stdize, "stdize".ptr),
    StrTable(GenFunc.f_stack, "stack".ptr),
    StrTable(GenFunc.f_geoplot, "geoplot".ptr),
    StrTable(GenFunc.f_bincoeff, "bincoeff".ptr),
    StrTable(GenFunc.f_tdisagg, "tdisagg".ptr),
    StrTable(GenFunc.f_assert, "assert".ptr),
    StrTable(GenFunc.f_vma, "vma".ptr),
    StrTable(GenFunc.f_bcheck, "bcheck".ptr),
    StrTable(GenFunc.f_contains, "contains".ptr),
    StrTable(GenFunc.f_lpsolve, "lpsolve".ptr),
    StrTable(GenFunc.f_distance, "distance".ptr),
    StrTable(GenFunc.f_interpol, "interpol".ptr),
    StrTable(GenFunc.f_mat2list, "mat2list".ptr),
    StrTable(GenFunc.f_dec2bin, "dec2bin".ptr),
    StrTable(GenFunc.f_bin2dec, "bin2dec".ptr),
    StrTable(GenFunc.f_commute, "commute".ptr),
    StrTable(GenFunc.f_access, "access".ptr),
    StrTable(GenFunc.f_sphcorr, "sphericorr".ptr),
    StrTable(GenFunc.f_threshold, "thresh".ptr),
    StrTable(GenFunc.f_lnmgamma, "lnmgamma".ptr),
    StrTable(GenFunc.f_which, "which".ptr),
    StrTable(GenFunc.f_gpcolors, "gpcolors".ptr),
    StrTable(GenFunc.f_permute, "permute".ptr),
    StrTable(0, null)];

__gshared StrTable[] func_alias = [
    StrTable(GenFunc.f_eigen, "eiggen2".ptr), /* deprecated */
    StrTable(GenFunc.f_fft, "fft2".ptr),    /* deprecated */
    StrTable(GenFunc.f_nmmax, "NMmin".ptr),
    StrTable(GenFunc.f_nrmax, "NRmin".ptr),
    StrTable(GenFunc.f_bfgsmax, "BFGSmin".ptr),
    StrTable(GenFunc.f_bfgscmax, "BFGScmin".ptr),
    StrTable(GenFunc.f_gssmax, "GSSmin".ptr),
    StrTable(GenFunc.f_gamma, "gammafunc".ptr),
    StrTable(GenFunc.f_gamma, "gamma".ptr),
    StrTable(GenFunc.f_log, "logs".ptr),
    StrTable(GenFunc.f_log, "ln".ptr),
    StrTable(GenFunc.f_square, "xpx".ptr),
    StrTable(GenFunc.f_backtick, "$".ptr),
    StrTable(GenFunc.f_cnameset, "colnames".ptr),
    StrTable(GenFunc.f_rnameset, "rownames".ptr),
    StrTable(GenFunc.f_cnameget, "colname".ptr),
    StrTable(GenFunc.f_rnameget, "rowname".ptr),
    StrTable(GenFunc.f_max, "xmax".ptr),
    StrTable(GenFunc.f_min, "xmin".ptr),
    StrTable(0, null)];

__gshared StrTable[] hidden_funcs = [
    StrTable(GenFunc.hf_jbterms, "_jbterms".ptr),
    StrTable(GenFunc.hf_listinfo, "_listinfo".ptr),
    StrTable(GenFunc.hf_regls, "_regls".ptr),
    StrTable(GenFunc.hf_felogitr, "_felogit_rec".ptr),
    StrTable(GenFunc.hf_fdepth, "_fdepth".ptr),
    StrTable(GenFunc.hf_glasso, "_glasso".ptr),
    StrTable(GenFunc.hf_vcnorm, "_vcnorm".ptr),
    StrTable(0, null)];

extern (C) int const_lookup(const(char) * s)
{
    int i;

    for (i=0; consts[i].id != 0; i++) {
        if (!strcmp(s, consts[i].str)) {
            return consts[i].id;
        }
    }

    return 0;
}

extern (C) const(char) * constname(int c)
{
    int i;

    for (i=0; consts[i].id != 0; i++) {
        if (c == consts[i].id) {
            return consts[i].str;
        }
    }

    return "unknown";
}

private GHashTable * gretl_function_hash_init()
{
    GHashTable *ht;
    int i;

    ht = g_hash_table_new(&g_str_hash, &g_str_equal);

    for (i=0; ptrfuncs[i].str != null; i++) {
        g_hash_table_insert(ht, cast(gpointer) ptrfuncs[i].str, &ptrfuncs[i]);
    }

    for (i=0; funcs[i].str != null; i++) {
        g_hash_table_insert(ht, cast(gpointer) funcs[i].str, &funcs[i]);
    }

    for (i=0; hidden_funcs[i].str != null; i++) {
        g_hash_table_insert(ht, cast(gpointer) hidden_funcs[i].str,
                            &hidden_funcs[i]);
    }

    return ht;
}

__gshared GHashTable * oht;

extern (C) int install_function_override(const(char) * funname, const(char) * pkgname, gpointer data)
{
    if (funname == null) {
        /* cleanup signal */
        if (oht != null) {
            g_hash_table_destroy(oht);
            oht = null;
        }
        return 0;
    }

    if (oht == null) {
        oht = g_hash_table_new_full(&g_str_hash, &g_str_equal,
                                    &g_free, null);
    }

    if (oht != null) {
        gchar *key = g_strdup_printf("%s::%s", pkgname, funname);

        g_hash_table_insert(oht, cast(gpointer) key, data);
    }

    return 0;
}

extern (C) int delete_function_override(const(char) * funname, const(char) * pkgname)
{
    int ret = 0;

    if (oht != null) {
        gchar *key = g_strdup_printf("%s::%s", pkgname, funname);

        if (g_hash_table_remove(oht, key)) {
            fprintf(stderr, "'%s': deleted override of built-in\n", key);
            ret = 1;
        }
        g_free(key);
    }

    return ret;
}

private Ufunc * get_function_override(const(char) * sf, gpointer p)
{
    const char *sp = function_package_get_name(cast(Fnpkg *) p);
    gchar *key = g_strdup_printf("%s::%s", sp, sf);
    Ufunc * uf = cast(Ufunc *) g_hash_table_lookup(oht, key);

    g_free(key);

    return uf;
}

/* Attention: this function may be called from function_loopkup()
   below, and in that context @p will be null -- so don't
   dereference @p without checking it for nullity first!
*/

private int real_function_lookup(const(char) * s, int aliases, Parser * p)
{
    static __gshared GHashTable * fht;
    gpointer fnp;

    if (s == null) {
        /* cleanup signal */
        if (fht != null) {
            g_hash_table_destroy(fht);
            fht = null;
        }
        return 0;
    }

    if (fht == null) {
        fht = gretl_function_hash_init();
    }

    fnp = g_hash_table_lookup(fht, s);
    if (fnp != null) {
        StrTable *st = cast(StrTable *) fnp;

        if (p != null && st.id > 0 && st.id < GenFunc.fp_max) {
            StrTableEx *sx = cast(StrTableEx *) fnp;

            p.data = sx.ptr;
        }
        if (p != null) {
            /* note: point d'appui for deprecation of built-in function */
            if (st.id == GenFunc.f_chowlin) {
                pprintf(p.prn, "*** Warning: %s() is obsolete, please use " ~
			"tdisagg() instead ***\n", st.str);
            }
        }
        return st.id;
    }

    if (aliases) {
        int i;

        for (i=0; func_alias[i].id != 0; i++) {
            if (!strcmp(s, func_alias[i].str)) {
                if (!strcmp(s, "fft2")) {
                    gretl_warnmsg_set(_("deprecated alias 'fft2': please call fft()").ptr);
                }
                if (p != null) {
                    p.flags |= GenFlags.aliased;
                }
                return func_alias[i].id;
            }
        }
    }

    return 0;
}

extern (C) void gretl_function_hash_cleanup()
{
    real_function_lookup(null, 0, null);
    install_function_override(null, null, null);
}

extern (C) int function_lookup(const(char) * s)
{
    return real_function_lookup(s, 0, null);
}

extern (C) int is_function_alias(const(char) * s)
{
    int i;

    for (i=0; func_alias[i].id != 0; i++) {
        if (!strcmp(s, func_alias[i].str)) {
            return 1;
        }
    }

    return 0;
}

extern (C) void * get_genr_function_pointer(int f)
{
    int i;

    for (i=0; ptrfuncs[i].str != null; i++) {
        if (ptrfuncs[i].id == f) {
            return ptrfuncs[i].ptr;
        }
    }

    return null;
}

private int function_lookup_with_alias(const(char) * s, Parser * p)
{
    if (oht != null) {
        /* we have a record of one or more package-private
           functions whose names collide with built-ins
        */
        gpointer pp = get_active_function_package(Gretlopt.OPT_O);

        if (pp != null) {
            Ufunc *uf = get_function_override(s, pp);

            if (uf != null) {
                p.idstr = gretl_strdup(s);
                p.data = uf;
                return GenSym.ufun;
            }
        }
    }

    return real_function_lookup(s, 1, p);
}

private const(char) * funname(int t)
{
    int i;

    for (i=0; ptrfuncs[i].id != 0; i++) {
        if (t == ptrfuncs[i].id) {
            return ptrfuncs[i].str;
        }
    }

    for (i=0; funcs[i].id != 0; i++) {
        if (t == funcs[i].id) {
            return funcs[i].str;
        }
    }

    for (i=0; hidden_funcs[i].id != 0; i++) {
        if (t == hidden_funcs[i].id) {
            return hidden_funcs[i].str;
        }
    }

    return "unknown";
}

private int show_alias(int i)
{
    if (strstr(func_alias[i].str, "min")) {
        return 1;
    } else {
        return 0;
    }
}

/* return the number of built-in functions */

extern (C) int gen_func_count()
{
    int i;
    int n = 0;

    for (i=0; ptrfuncs[i].id != 0; i++) {
        n++;
    }

    for (i=0; funcs[i].id != 0; i++) {
        n++;
    }

    for (i=0; func_alias[i].id != 0; i++) {
        if (show_alias(i)) {
            n++;
        }
    }

    return n;
}

/* return the name of function @i, including aliases */

extern (C) const(char) * gen_func_name(int i)
{
    int j;
    int seq = -1;

    for (j=0; ptrfuncs[j].id != 0; j++) {
        seq++;
        if (seq == i) {
            return ptrfuncs[i].str;
        }
    }

    for (j=0; funcs[j].id != 0; j++) {
        seq++;
        if (seq == i) {
            return funcs[j].str;
        }
    }

    for (j=0; func_alias[j].id != 0; j++) {
        if (show_alias(j)) {
            seq++;
        }
        if (seq == i) {
            return func_alias[j].str;
        }
    }

    return null;
}

extern (C) int model_var_count()
{
    int i;

    for (i=0; mvars[i].id != 0; i++) {}
    return i;
}

extern (C) const(char) * model_var_name(int i)
{
    return mvars[i].str;
}

extern (C) int bundle_var_count()
{
    int i;

    for (i=0; bvars[i].id != 0; i++) {}
    return i;
}

extern (C) const(char) * bundle_var_name(int i)
{
    return bvars[i].str;
}

extern (C) int data_var_count()
{
    int i;
    int n = 0;

    for (i=0; dvars[i].id != 0; i++) {
        if (dvars[i].str[0] == '$') {
            n++;
        }
    }

    return n;
}

extern (C) const(char) * data_var_name(int i)
{
    return dvars[i].str;
}

extern (C) const(char) * gretl_function_complete(const(char) * s)
{
    size_t n = strlen(s);
    int i;

    for (i=0; ptrfuncs[i].str != null; i++) {
        if (!strncmp(s, ptrfuncs[i].str, n)) {
            return ptrfuncs[i].str;
        }
    }

    for (i=0; funcs[i].str != null; i++) {
        if (!strncmp(s, funcs[i].str, n)) {
            return funcs[i].str;
        }
    }

    return null;
}

extern (C) int gretl_const_count()
{
    int i;

    for (i=0; consts[i].id != 0; i++) {}
    return i;
}

extern (C) const(char) * gretl_const_name(int i)
{
    return consts[i].str;
}

/* end external stuff */

private int dummy_lookup(const(char) * s, Parser * p)
{
    int i;
    int d = 0;

    for (i=0; dummies[i].id != 0; i++) {
        if (!strcmp(s, dummies[i].str)) {
            d = dummies[i].id;
            break;
        }
    }

    if (d == DumId.end) {
        ; /* let this pass: we'll check its validity later */
    } else if (MSEL_DUM(d) && parser_next_char(p) != ']') {
        /* most MSEL dummies are stand-alone; they can
           be valid only if followed by ']'
        */
        d = 0;
    } else if (d == DumId.empty && strcmp(p.rhs, "empty")) {
        /* if "empty" is not the only term on the right-hand
           side, convert it to "null"
        */
        p.sym = GenSym.empty;
    }

    return d;
}

extern (C) const(char) * dumname(int t)
{
    int i;

    for (i=0; dummies[i].id != 0; i++) {
        if (t == dummies[i].id) {
            return dummies[i].str;
        }
    }

    return "unknown";
}

private int dvar_lookup(const(char) * s)
{
    int i;

    for (i=0; dvars[i].id != 0; i++) {
        if (!strcmp(s, dvars[i].str)) {
            return dvars[i].id;
        }
    }

    return 0;
}

extern (C) const(char) * dvarname(int t)
{
    int i;

    for (i=0; dvars[i].id != 0; i++) {
        if (t == dvars[i].id) {
            return dvars[i].str;
        }
    }

    return "unknown";
}

extern (C) int mvar_lookup(const(char) * s)
{
    int i;

    for (i=0; mvars[i].id != 0; i++) {
        if (!strcmp(s, mvars[i].str)) {
            return mvars[i].id;
        }
    }

    /* aliases */

    if (!strcmp(s, "$nrsq")) {
        return ModelDataIndex.trsq;
    } else if (!strcmp(s, "$fcerr")) {
        return ModelDataIndex.fcse;
    }

    return 0;
}

extern (C) const(char) * mvarname(int t)
{
    int i;

    for (i=0; mvars[i].id != 0; i++) {
        if (t == mvars[i].id) {
            return mvars[i].str;
        }
    }

    return "unknown";
}

extern (C) int bvar_lookup(const(char) * s)
{
    int i;

    for (i=0; bvars[i].id != 0; i++) {
        if (!strcmp(s, bvars[i].str)) {
            return bvars[i].id;
        }
    }

    return 0;
}

extern (C) const(char) * bvarname(int t)
{
    int i;

    for (i=0; bvars[i].id != 0; i++) {
        if (t == bvars[i].id) {
            return bvars[i].str;
        }
    }

    return "unknown";
}

extern (C) int genr_function_word(const(char) * s)
{
    int ret = 0;

    ret = real_function_lookup(s, 0, null);
    if (!ret) {
        ret = dvar_lookup(s);
    }
    if (!ret) {
        ret = mvar_lookup(s);
    }
    if (!ret) {
        ret = bvar_lookup(s);
    }
    if (!ret) {
        ret = const_lookup(s);
    }

    return ret;
}

extern (C) int parser_ensure_error_buffer(Parser * p)
{
    if (p.prn == null && p.errprn == null) {
        p.errprn = gretl_print_new(PrnType.buffer, null);
        if (p.errprn != null) {
            p.prn = p.errprn;
            return 0;
        } else {
            return GretlError.alloc;
        }
    }

    return 0;
}

extern (C) void undefined_symbol_error(const(char) * s, Parser * p)
{
    parser_ensure_error_buffer(p);
    parser_print_input(p);

    if (p.ch == '.') {
        pprintf(p.prn, _("%s: no such object").ptr, s);
    } else {
        pprintf(p.prn, _("The symbol '%s' is undefined").ptr, s);
    }
    pputc(p.prn, '\n');
    p.err = GretlError.data;
}

private void function_noargs_error(const(char) * s, Parser * p)
{
    parser_ensure_error_buffer(p);
    parser_print_input(p);

    pprintf(p.prn, _("'%s': no argument was given").ptr, s);
    pputc(p.prn, '\n');
    p.err = GretlError.args;
}

extern (C) void context_error(int c, Parser * p, const(char) * func)
{
    parser_ensure_error_buffer(p);
    if (c != 0) {
        parser_print_input(p);
        pprintf(p.prn, _("The symbol '%c' is not valid in this context\n").ptr, c);
        if (c == '&') {
            pputs(p.prn, _("(for logical AND, use '&&')\n").ptr);
        } else if (c == '|') {
            pputs(p.prn, _("(for logical OR, use '||')\n").ptr);
        } else if (c == ',') {
            p.err = GretlError.parse;
        }
    } else if (p.sym == GenSym.eot) {
        parser_print_input(p);
        pputs(p.prn, _("Incomplete expression\n").ptr);
    } else {
        const char *s = getsymb_full(p.sym, p);

        if (s != null && *s != '\0' && strcmp(s, "unknown")) {
            pprintf(p.prn, _("The symbol '%s' is not valid in this context\n").ptr,
                    getsymb_full(p.sym, p));
        } else {
            parser_print_input(p);
        }
    }

    if (!p.err) {
        p.err = GretlError.parse;
    }
}

/* @parsing_query: we want to keep track of the case where we're
   lexing/parsing the branches of a ternary "query" expression. When
   such an expression is evaluated, it's OK if the branch _not_ taken
   contains an undefined symbol; indeed, this can occur by design, as
   in

     scalar y = exists(x) ? x : 0

   when "x" is in fact undefined.

   We therefore use the "UNDEF" node type to defuse the error that
   would otherwise arise on parsing. An error is triggered only if the
   branch that references the UNDEF node is selected (attempting to
   evaluate an UNDEF node automatically throws an error.)
*/

static __gshared int parsing_query;

extern (C) void set_parsing_query(int s)
{
    parsing_query = s;
}

/* 2023-09-15: apparatus for allowing \" as an embedded quote in string
   literals. Double-quote is taken as closing a quotation (rather than
   as an embedded literal) if (a) it immediately follows an opening
   quote, (b) it is not immediately preceded by a backslash, or (c) it
   is the last byte in a hansl statement.
*/

private int alt_closing_quote_pos(const(char) * s, int * esc)
{
    int i;
    int ret = -1;

    for (i=0; s[i]; i++) {
        if (s[i] == '"') {
            if (i == 0 || s[i-1] != '\\' || s[i+1] == '\0') {
                ret = i;
                break;
            } else {
                *esc = 1;
            }
        }
    }

    return ret;
}

private char * escape_quotes(char * s)
{
    int i;

    for (i=1; s[i]; i++) {
        if (s[i] == '"' && s[i-1] == '\\') {
            shift_string_left(s + i - 1, 1);
        }
    }

    return s;
}

private char * get_quoted_string(Parser * p, int prevsym)
{
    char * s = null;
    int esc = 0;
    int n;



    if (prevsym == GenFunc.f_sprintf || prevsym == GenFunc.f_printf) {
        /* look for a matching non-escaped double-quote,
           allowance made for "\\" as itself an escape
        */
        n = double_quote_position(p.point);
    } else {
        /* look for a matching non-escaped double-quote when
           backslash is special only when preceding a double
           quote
        */
        n = alt_closing_quote_pos(p.point, &esc);
    }

    if (n >= 0) {
        s = gretl_strndup(p.point, n);
        if (esc) {
            /* "\"" is accepted as an escape, but that's all: FIXME */
            escape_quotes(s);
        }
        parser_advance(p, n + 1);
    } else {
        parser_print_input(p);
        pprintf(p.prn, _("Unmatched '%c'\n").ptr, '"');
        p.err = GretlError.parse;
    }

    if (!p.err) {
        if (p.ch == '.' && *p.point == '$') {
            /* maybe quoted name of saved model followed by
               dollar variable? */
            p.sym = GenSym.mmemb;
        } else {
            p.sym = GenSym.cstr;
        }
    }

    return s;
}

private int might_be_date_string(const(char) * s, int n)
{
    char[12] test;
    int y;
    int m;
    int d;


    if (n > 10) {
        return 0;
    }

    test[0] = 0;
    strncat(test.ptr, s, n);

    if (strspn(s, "1234567890") == n) {
        /* plain integer (FIXME?) */
        return 1;
    } else if (sscanf(s, "%d:%d", &y, &m) == 2) {
        /* quarterly, monthly date */
        return 1;
    } else if (sscanf(s, "%d-%d-%d", &y, &m, &d) == 3) {
        /* daily date? */
        return 1;
    } else if (sscanf(s, "%d/%d/%d", &y, &m, &d) == 3) {
        /* daily date? */
        return 1;
    }

    return 0;
}

extern (C) Node * obs_node(Parser * p)
{
    Node *ret = null;
    char[OBSLEN + 2] word = 0;
    const(char) * s = p.point - 1;
    int close;
    int special = 0;
    int t = -1;

    close = gretl_charpos(']', s);


    if (close == 0) {
        pprintf(p.prn, _("Empty observation []\n").ptr);
        p.err = GretlError.parse;
    } else if (close < 0) {
        pprintf(p.prn, _("Unmatched '%c'\n").ptr, '[');
        p.err = GretlError.parse;
    } else if (*s == '"' && close < OBSLEN + 2 &&
               gretl_charpos('"', s+1) == close - 2) {
        /* quoted observation label? */
        strncat(word.ptr, s, close);
        special = 1;
    } else if (might_be_date_string(s, close)) {
        strncat(word.ptr, s, close);
        special = 1;
    }

    if (special && !p.err) {
        t = get_t_from_obs_string(word.ptr, p.dset);
        if (t >= 0) {
            /* convert to user-style 1-based index */
            t++;
        }
    }

    if (t > 0) {
        parser_advance(p, close - 1);
        lex(p);
        ret = newdbl(t);
    } else if (!p.err) {
        lex(p);
        ret = expr(p);
    }

    return ret;
}

extern (C) int is_gretl_accessor(const(char) * s)
{
    int i;
    int n;

    for (i=0; dvars[i].id != 0; i++) {
        n = cast(int) strlen(dvars[i].str);
        if (!strncmp(s, dvars[i].str, n)) {
            return !isalpha(s[n]);
        }
    }

    for (i=0; mvars[i].id != 0; i++) {
        n = cast(int) strlen(mvars[i].str);
        if (!strncmp(s, mvars[i].str, n)) {
            return !isalpha(s[n]);
        }
    }

    for (i=0; bvars[i].id != 0; i++) {
        n = cast(int) strlen(bvars[i].str);
        if (!strncmp(s, bvars[i].str, n)) {
            return !isalpha(s[n]);
        }
    }

    return 0;
}

private void look_up_dollar_word(const(char) * s, Parser * p)
{
    char * bstr;

    if ((p.idnum = dvar_lookup(s)) > 0) {
        p.sym = GenSym.dvar;
    } else if ((p.idnum = const_lookup(s)) > 0) {
        p.sym = GenSym.con;
    } else if ((p.idnum = mvar_lookup(s)) > 0) {
        p.sym = GenSym.mvar;
    } else if ((p.idnum = bvar_lookup(s)) > 0) {
        p.sym = GenSym.dbundle;
    } else if ((bstr = get_built_in_string_by_name(s+1)) != null) {
        p.sym = GenSym.cstr;
        p.idstr = gretl_strdup(bstr);
    } else {
        undefined_symbol_error(s, p);
    }

}

int maybe_get_R_function(const(char) * s) { return 0; }

static __gshared int doing_genseries;

extern (C) void set_doing_genseries(int s)
{
    doing_genseries = s;
}

/* Get the next non-space byte beyond what's already parsed:
   this will either be p.ch, or may be found at p.point
   or beyond.
*/

private int parser_next_char(Parser * p)
{
    if (p.ch != ' ') {
        return p.ch;
    } else {
        const(char) * s = p.point;

        while (*s) {
            if (!isspace(*s)) {
                return *s;
            }
            s++;
        }
        return 0;
    }
}

private int char_past_point(Parser * p)
{
    if (*p.point != '\0') {
        int i;

        for (i=1; p.point[i] != '\0'; i++) {
            if (!isspace(p.point[i])) {
                return p.point[i];
            }
        }
    }
    return 0;
}

/* helpers specific to the stack() function -- dispensable if
   we decide not to support the old-style syntax any more
*/

private int get_oldstyle_stack_args(const(char) * s, char ** arg, char ** opt1, char ** opt2)
{
    const(char) * p1;
    const(char) * p2;
    int len;

    /* Let the marker for old-style uage be the presence of
       one or more legacy option flags */

    p1 = strstr(s, "--length=");
    if (p1 != null) {
        len = cast(int) strcspn(p1 + 9, " \n");
        *opt1 = gretl_strndup(p1 + 9, len);
    }

    p2 = strstr(s, "--offset=");
    if (p2 != null) {
        len = cast(int) strcspn(p2 + 9, " \n");
        *opt2 = gretl_strndup(p2 + 9, len);
    }

    if (p1 != null || p2 != null) {
        len = cast(int) strcspn(s, ")");
        *arg = gretl_strndup(s, len);
        return 1;
    } else {
        return 0;
    }
}

private int stack_update_parser_input(Parser * p)
{
    char * opt1 = null;
    char * opt2 = null;
    char * arg = null;
    char * start = null;
    const(char) * s;
    char * tmp = null;
    GString * gs;
    int offset;

    offset = cast(int) (p.point - p.input);
    start = gretl_strndup(p.input, offset);
    gs = g_string_new(start);
    free(start);

    s = strstr(p.input, "stack(") + 6;
    get_oldstyle_stack_args(s, &arg, &opt1, &opt2);
    if (arg != null) {
        gs = g_string_append(gs, arg);
        free(arg);
    }
    if (opt1 != null) {
        gs = g_string_append_c(gs, ',');
        gs = g_string_append(gs, opt1);
        free(opt1);
    }
    if (opt2 != null) {
        gs = g_string_append_c(gs, ',');
        gs = g_string_append(gs, opt2);
        free(opt2);
    }
    gs = g_string_append_c(gs, ')');
    tmp = g_string_free(gs, 0);

    p.input = tmp;
    p.point = p.input + offset;
    p.flags |= GenFlags.altinp;

    return 0;
}

private void handle_lpnext(const(char) * s, Parser * p, int have_dset)
{
    Ufunc *u = get_user_function_by_name(s);
    int vnum = -1;

    if (have_dset) {
        vnum = current_series_index(p.dset, s);
    }

    if (u == null && vnum >= 0) {
        /* unambiguous: series */
        p.sym = GenSym.series;
    } else if (u != null && vnum < 0) {
        /* unambiguous: function */
        p.sym = GenSym.ufun;
    } else if (u != null) {
        /* ambiguous case! */
        if (gretl_function_depth() > 0) {
            /* function writers should avoid collisions
               when naming series
            */
            p.sym = GenSym.series;
        } else if (p.targ != GenSym.unk && p.targ != GenSym.list && p.targ != GenSym.series) {
            /* target not compatible with series lag? */
            p.sym = GenSym.ufun;
        } else if (defining_list(p) || p.targ == GenSym.series) {
            /* debatable */
            p.sym = GenSym.series;
        } else {
            /* debatable */
            p.sym = GenSym.ufun;
        }
    }

    if (p.sym != 0) {
        p.idstr = gretl_strdup(s);
        if (p.sym == GenSym.ufun) {
            p.data = u;
        } else {
            p.idnum = vnum;
            /* in case of any intervening space */
            while (p.ch == ' ') {
                parser_getc(p);
            }
        }
    }
}

private int is_function_word(const(char) * s)
{
    return function_lookup_with_alias(s, null) != 0 ||
        get_user_function_by_name(s) != null;
}

private void look_up_word(const(char) * s, Parser * p)
{
    int have_dset = (p.dset != null && p.dset.v > 0);
    int prevsym = p.sym;
    int lpnext;
    int err = 0;


    /* is the next (or next non-space) character left paren? */
    lpnext = parser_next_char(p) == '(';

    /* initialize */
    p.sym = 0;
    p.data = null;

    if (lpnext) {
        /* identifier is immediately followed by left paren:
           most likely a function call but could be the name
           of a series followed by a lag specifier
        */
        p.sym = function_lookup_with_alias(s, p);
        if (p.sym == 0) {
            handle_lpnext(s, p, have_dset);
        } else if (p.sym == GenFunc.f_stack) {
            /* special! */
            if (strstr(p.point, "--length") || strstr(p.point, "--offset")) {
                stack_update_parser_input(p);
            }
            p.flags |= GenFlags.stack;
        }
    }

    if (p.sym == 0) {
        p.idnum = const_lookup(s);
        if (p.idnum > 0) {
            p.sym = GenSym.con;
        } else {
            p.idnum = dummy_lookup(s, p);
            if (p.idnum > 0) {
                p.sym = GenSym.dum;
            } else if (have_dset &&
                       (p.idnum = current_series_index(p.dset, s)) >= 0) {
                p.sym = GenSym.series;
                p.idstr = gretl_strdup(s);
            } else if (have_dset && !strcmp(s, "time")) {
                p.sym = GenSym.dum;
                p.idnum = DumId.trend;
            } else if ((p.data = get_user_var_by_name(s)) != null) {
                UserVar * u = cast(UserVar *) p.data;

                if (u.type == GretlType.double_) {
                    p.sym = GenSym.num;
                } else if (u.type == GretlType.matrix) {
                    p.sym = GenSym.mat;
                } else if (u.type == GretlType.bundle) {
                    p.sym = GenSym.bundle;
                } else if (u.type == GretlType.array) {
                    p.sym = GenSym.array;
                } else if (u.type == GretlType.string_) {
                    p.sym = GenSym.str;
                } else if (u.type == GretlType.list) {
                    p.sym = GenSym.list;
                }
                p.idstr = gretl_strdup(s);
            } else if (defining_list(p) && varname_match_any(p.dset, s)) {
                p.sym = GenSym.wlist;
                p.idstr = gretl_strdup(s);
            } else if (have_dset && !strcmp(s, "t")) {
                /* if "t" has not been otherwise defined, treat it
                   as an alias for "obs"
                */
                p.sym = GenSym.dvar;
                p.idnum = RetrievalIndex.index;
            } else if (maybe_get_R_function(s)) {
                /* note: all "native" types take precedence over this */
                p.sym = GenSym.rfun;
                p.idstr = gretl_strdup(s + 2);
            } else if (parsing_query || prevsym == GenSym.b_and || prevsym == GenSym.b_or) {
                p.sym = GenSym.undef;
                p.idstr = gretl_strdup(s);
            } else if (p.flags & (GenFlags.and | GenFlags.or)) {
                p.sym = GenSym.undef;
                p.idstr = gretl_strdup(s);
            } else if (gretl_get_object_by_name(s)) {
                p.sym = GenSym.uobj;
                p.idstr = gretl_strdup(s);
            } else {
                err = GretlError.unkvar;
            }
        }
    }

    if (err) {
        if (is_function_word(s)) {
            /* @s is a function identifier with no
               following left paren */
            function_noargs_error(s, p);
        } else if (object_is_function_arg(s)) {
            /* @s is the name of a function parameter
               with a null value */
            p.sym = GenSym.nullarg;
            p.idstr = gretl_strdup(s);
        } else if (p.sym != DumId.empty) {
            undefined_symbol_error(s, p);
        }
    }

}

private void maybe_treat_as_postfix(Parser * p)
{
    if (p.sym == GenSym.num) {
        const(char) * ok = ")]}+-*/%,:";
        int c = char_past_point(p);

        /* Interpret as foo++ or foo-- ? Only if
           the following character is suitable.
        */
        if (c == 0 || strchr(ok, c)) {
            p.sym = p.ch == '+'? GenSym.num_p : GenSym.num_m;
            /* swallow the pluses or minuses */
            parser_advance(p, 1);
        }
    }
}

bool dollar_series(int t) { return t > RetrievalIndex.scalar_max && t < RetrievalIndex.series_max; }

bool could_be_matrix(int t)
{
    return model_data_matrix(t) || model_data_matrix_builder(t) || t == ModelDataIndex.uhat || t == ModelDataIndex.yhat || (t > RetrievalIndex.series_max && t < RetrievalIndex.max);
}

bool could_be_bundle(int t) { return t == RetrievalIndex.result; }

private void word_check_next_char(Parser * p)
{
    char[2] chk = [cast(char) p.ch, '\0'];

    p.upsym = 0;

    /* Here we're checking for validity and syntactical implications
       of one of "([.+-" immediately following a 'word' of some kind.
    */
    if (strspn(chk.ptr, "([.+-".ptr) == 0) {
        /* none of the above */
        return;
    }

    if (p.sym == GenSym.undef) {
        /* 2020-03-16: suspend disbelief, we might be in a branch
           of "cond ? x : y" that ends up not being taken */
        return;
    }

    if (p.ch == '(') {
        /* series or list (lag), or function */
        if (p.sym == GenSym.series) {
            if (p.idnum > 0 && p.idnum == p.lh.vnum) {
                p.flags |= GenFlags.autoreg;
            }
            p.upsym = p.sym;
            p.sym = GenSym.lag;
        } else if (p.sym == GenSym.list) {
            p.upsym = p.sym;
            p.sym = GenSym.lag;
        } else if (p.sym == GenSym.mvar && model_data_matrix(p.idnum)) {
            /* old-style "$coeff(x1)" etc. */
            p.sym = GenSym.dmstr;
        } else if (!func1_symb(p.sym) &&
                   !func2_symb(p.sym) &&
                   !func3_symb(p.sym) &&
                   !funcn_symb(p.sym) &&
                   p.sym != GenSym.ufun &&
                   p.sym != GenSym.rfun) {
            p.err = GretlError.parse;
        }
    } else if (p.ch == '[') {
        p.upsym = p.sym;
        if (p.sym == GenSym.mat || p.sym == GenSym.array ||
            p.sym == GenSym.list || p.sym == GenSym.str) {
            /* slice of sliceable object */
            p.sym = GenSym.osl;
        } else if (*p.point != '"' &&
                   (p.sym == GenSym.mvar || p.sym == GenSym.dvar) &&
                   could_be_matrix(p.idnum)) {
            /* slice of $-matrix? */
            p.sym = GenSym.osl;
        } else if (p.sym == GenSym.series) {
            /* observation from series */
            p.sym = GenSym.obs;
        } else if (p.sym == GenSym.dvar && dollar_series(p.idnum)) {
            /* observation from "dollar" series */
            p.sym = GenSym.obs;
        } else if (p.sym == GenSym.mvar && (model_data_list(p.idnum) ||
                                      model_data_array(p.idnum))) {
            /* element/range of accessor list or array */
            p.sym = GenSym.osl;
        } else if (p.sym == GenSym.bundle) {
            /* member from bundle */
            p.sym = GenSym.bmemb;
        } else if (p.sym == GenSym.dbundle) {
            /* member from $ bundle */
            p.sym = GenSym.dbmemb;
        } else if (p.sym == GenSym.dvar && could_be_bundle(p.idnum)) {
            p.sym = GenSym.dbmemb;
        } else {
            p.err = GretlError.parse;
        }
    } else if (p.ch == '.' && *p.point == '$') {
        if (p.sym == GenSym.uobj) {
            /* name of saved model followed by dollar variable? */
            p.sym = GenSym.mmemb;
        } else if (p.sym == GenSym.cstr) {
            /* maybe quoted name of saved object followed by
               dollar variable? */
            p.sym = GenSym.mmemb;
        } else {
            p.err = GretlError.parse;
        }
    } else if (p.ch == '.' && isalpha(*p.point)) {
        if (p.sym == GenSym.list) {
            p.sym = GenSym.listvar;
        } else if (p.sym == GenSym.bundle) {
            p.sym = GenSym.bmemb;
        } else if (p.sym == GenSym.dbundle) {
            p.sym = GenSym.dbmemb;
        } else if (p.sym == GenSym.dvar && could_be_bundle(p.idnum)) {
            p.sym = GenSym.dbmemb;
        } else {
            p.err = GretlError.parse;
        }
    } else if (p.ch == '+' && *p.point == '+') {
        maybe_treat_as_postfix(p);
    } else if (p.ch == '-' && *p.point == '-') {
        maybe_treat_as_postfix(p);
    }

    if (p.err) {
        context_error(p.ch, p, "word_check_next_char");
    }
}

private int is_word_char(Parser * p)
{
    if (strchr(wordchars, p.ch) != null) {
        return 1;
    } else if (defining_list(p) && !doing_genseries &&
               (p.ch == '*' || p.ch == '?')) {
        return 1;
    }

    return 0;
}

private void getword(Parser * p, int greek)
{
    char[32] word;
    int i = 0;

    if (greek) {
        /* we have a single (2-byte) UTF-8 greek letter in scope */
        for (i=0; i<2; i++) {
            word[i] = cast(char) p.ch;
            parser_getc(p);
        }
    } else {
        /* we know the first char is acceptable (and might be '$' or '_') */
        word[i++] = cast(char) p.ch;
        parser_getc(p);

        while (p.ch != 0 && is_word_char(p) && i < 31) {
            word[i++] = cast(char) p.ch;
            parser_getc(p);
        }
    }

    word[i] = '\0';


    while (p.ch != 0 && strchr(wordchars, p.ch) != null) {
        /* flush excess word characters */
        parser_getc(p);
    }

    if (p.flags & GenFlags.getstr) {
        /* uninterpreted string wanted */
        p.sym = GenSym.cstr;
        p.idstr = gretl_strdup(word.ptr);
        p.flags ^= GenFlags.getstr;
        return; /* FIXME bundle-member name */
    } else if ((*word.ptr == '$' && word[1]) || !strcmp(word.ptr, "obs".ptr)) {
        look_up_dollar_word(word.ptr, p);
    } else if (*word.ptr == '$' && word[1] == '\0' && (p.ch == '[' || p.ch == '.')) {
        p.sym = GenSym.dbundle;
        p.idnum = BundleDataIndex.model;
    } else {
        look_up_word(word.ptr, p);
    }

    if (!p.err) {
        word_check_next_char(p);
    }

}

private int colon_ok(Parser * p, char * s, int n)
{
    int i;

    if (p.flags & GenFlags.slicing) {
        /* calculating a matrix "slice": colon is a separator in
           this context, cannot be part of a time/panel observation
           string
        */
        return 0;
    }

    if (n != 1 && n != 3) {
        return 0;
    }

    for (i=0; i<=n; i++) {
        if (!isdigit(s[i])) {
            return 0;
        }
    }

    return 1;
}

/* Below: we're testing p.ch for validity, given what we've already
   packed into string @s up to element @i, and with some regard to
   the next character in line.
*/

private int ok_dbl_char(Parser * p, char * s, int i)
{
    int ret = 0;

    if (i < 0 || (p.ch >= '0' && p.ch <= '9')) {
        return 1;
    }

    switch (p.ch) {
    case '+':
    case '-':
        ret = s[i] == 'e' || s[i] == 'E';
        break;
    case '.':
        ret = !strchr(s, '.') && !strchr(s, ':') &&
            !strchr(s, 'e') && !strchr(s, 'E') &&
            *p.point != '.';
        break;
    case 'e':
    case 'E':
        ret = !strchr(s, 'e') && !strchr(s, 'E') &&
            !strchr(s, ':');
        break;
    case ':':
        /* allow for obs numbers in the form, e.g., "1995:10" */
        ret = colon_ok(p, s, i);
        break;
    default:
        break;
    }

    return ret;
}

private int ok_hex_char(Parser * p, char * s, int i)
{
    if (i < 1) {
        return 1;
    } else if (p.ch >= '0' && p.ch <= '9') {
        return 1;
    } else if (p.ch >= 'a' && p.ch <= 'f') {
        return 1;
    } else if (p.ch >= 'A' && p.ch <= 'F') {
        return 1;
    }

    return 0;
}

private void parse_number(Parser * p)
{
    char[NUMLEN] xstr = 0;
    int got_colon = 0;
    int got_hex = 0;
    int i = 0;

    if (p.ch == '0' && (p.point[0] == 'x' || p.point[0] == 'X')) {
        /* hexadecimal input */
        got_hex = 1;
        while (ok_hex_char(p, xstr.ptr, i - 1) && i < HEXLEN - 1) {
            xstr[i++] = cast(char) p.ch;
            parser_getc(p);
        }
    } else {
        /* decimal input */
        while (ok_dbl_char(p, xstr.ptr, i - 1) && i < NUMLEN - 1) {
            xstr[i++] = cast(char) p.ch;
            if (p.ch == ':') {
                got_colon = 1;
            }
            parser_getc(p);
        }
        while (p.ch >= '0' && p.ch <= '9') {
            /* flush excess numeric characters */
            parser_getc(p);
        }
    }


    if (got_colon) {
        if (p.dset == null || p.dset.n == 0) {
            p.err = GretlError.nodata;
        } else if (p.dset.pd == 1) {
            p.err = GretlError.pdwrong;
        } else if (dateton(xstr.ptr, p.dset) < 0) {
            p.err = GretlError.data;
        } else {
            p.idstr = gretl_strdup(xstr.ptr);
            p.sym = GenSym.cstr;
        }
    } else if (got_hex) {
        guint32 u = cast(guint32) strtoul(xstr.ptr, null, 16);

        p.xval = cast(double) u;
        p.sym = GenSym.cnum;
    } else {
        p.xval = dot_atof(xstr.ptr);
        p.sym = GenSym.cnum;
    }
}

private int wildcard_special(Parser * p)
{
    char cprev = *(p.point - 2);

    if (p.ch == '?') {
        char cnext = *p.point;

        if ((cprev == ' ' || cprev == ')') && cnext == ' ') {
            /* '?' is presumably ternary operator */
            return 0;
        }
    }

    if (cprev == ' ' &&
        (bare_data_type(p.sym) || closing_sym(p.sym) ||
         (p.sym == GenSym.lag))) {
        p.sym = GenSym.b_lcat;
    } else {
        getword(p, 0);
    }

    return 1;
}

/* 0xE2 0x88 0x92 = UTF-8 minus */

private void lex_try_utf8(Parser * p)
{
    if (cast(ubyte) *p.point == 0x88 &&
        cast(ubyte) *(p.point + 1) == 0x92) {
        p.sym = GenSym.b_sub;
        parser_getc(p);
        parser_getc(p);
        parser_getc(p);
    } else {
        pprintf(p.prn, _("Unexpected byte 0x%x\n").ptr,
                cast(ubyte) p.ch);
        p.err = GretlError.parse;
    }
}

bool word_start_special(int c) { return c == '$' || c == '_'; }

/* accept 'to', but only with spaces before and after */
bool lag_range_sym(const(Parser) * p)
{
    /* accept 'to', but only with spaces before and after */
    return (p.flags & GenFlags.lagprse) && p.ch == 't' && *p.point == 'o' && *(p.point - 2) == ' ' && *(p.point + 1) == ' ';
}

extern (C) void lex(Parser * p)
{
    static __gshared int[2] prevsyms;

    prevsyms[0] = prevsyms[1];
    prevsyms[1] = p.sym;

    if (p.ch == 0) {
        p.sym = GenSym.eot;
        return;
    }

    while (p.ch != 0) {
        if (cast(ubyte) p.ch == 0xE2) {
            lex_try_utf8(p);
            return;
        } else if (is_greek_letter(p.point - 1)) {
            getword(p, 1);
            return;
        }
        switch (p.ch) {
        case ' ':
        case '\t':
        case '\r':
        case '\n':
            parser_getc(p);
            break;
        case '+':
            p.sym = GenSym.b_add;
            parser_getc(p);
            return;
        case '-':
            p.sym = GenSym.b_sub;
            parser_getc(p);
            return;
        case '*':
            if (defining_list(p) && !doing_genseries) {
                /* allow for '*' as wildcard */
                wildcard_special(p);
                return;
            }
            parser_getc(p);
            if (p.ch == '*') {
                p.sym = GenSym.b_kron;
                parser_getc(p);
            } else {
                p.sym = GenSym.b_mul;
            }
            return;
        case '\'':
            p.sym = GenSym.b_trmul;
            parser_getc(p);
            return;
        case '/':
            p.sym = GenSym.b_div;
            parser_getc(p);
            return;
        case '\\':
            p.sym = GenSym.b_ldiv;
            parser_getc(p);
            return;
        case '%':
            p.sym = GenSym.b_mod;
            parser_getc(p);
            return;
        case '^':
            p.sym = GenSym.b_pow;
            parser_getc(p);
            return;
        case '&':
            parser_getc(p);
            if (p.ch == '&') {
                p.sym = GenSym.b_and;
                parser_getc(p);
            } else {
                p.sym = GenSym.u_addr;
            }
            return;
        case '|':
            parser_getc(p);
            if (p.ch == '|') {
                p.sym = GenSym.b_or;
                parser_getc(p);
            } else {
                p.sym = GenSym.b_vcat;
            }
            return;
        case '!':
            parser_getc(p);
            if (p.ch == '=') {
                p.sym = GenSym.b_neq;
                parser_getc(p);
            } else {
                p.sym = GenSym.u_not;
            }
            return;
        case '=':
            parser_getc(p);
            if (p.ch == '=') {
                parser_getc(p);
                p.sym = GenSym.b_eq;
            } else {
                gretl_errmsg_set(_("If you meant to test for " ~
			"equality, please use '=='").ptr);
                p.err = GretlError.parse;
            }
            return;
        case '>':
            parser_getc(p);
            if (p.ch == '=') {
                p.sym = GenSym.b_gte;
                parser_getc(p);
            } else {
                p.sym = GenSym.b_gt;
            }
            return;
        case '<':
            parser_getc(p);
            if (p.ch == '=') {
                p.sym = GenSym.b_lte;
                parser_getc(p);
            } else if (p.ch == '>') {
                p.sym = GenSym.b_neq;
                parser_getc(p);
            } else {
                p.sym = GenSym.b_lt;
            }
            return;
        case '(':
            p.sym = GenSym.g_lpr;
            parser_getc(p);
            return;
        case ')':
            p.sym = GenSym.g_rpr;
            parser_getc(p);
            return;
        case '[':
            p.sym = GenSym.g_lbr;
            parser_getc(p);
            return;
        case '{':
            p.sym = GenSym.g_lcb;
            parser_getc(p);
            return;
        case '}':
            p.sym = GenSym.g_rcb;
            parser_getc(p);
            return;
        case ']':
            p.sym = GenSym.g_rbr;
            parser_getc(p);
            return;
        case '~':
            p.sym = GenSym.b_hcat;
            parser_getc(p);
            return;
        case ',':
            p.sym = GenSym.p_com;
            parser_getc(p);
            return;
        case ';':
            if (defining_list(p)) {
                p.sym = GenSym.b_join;
            } else {
                /* used in matrix definition */
                p.sym = GenSym.p_semi;
            }
            parser_getc(p);
            return;
        case ':':
            p.sym = GenSym.p_col;
            parser_getc(p);
            return;
        case '?':
            if (defining_list(p) && !doing_genseries) {
                /* allow for '?' as wildcard */
                if (wildcard_special(p)) {
                    return;
                }
            }
            p.sym = GenSym.query;
            parser_getc(p);
            return;
        case '.':
            if (*p.point == '$') {
                p.sym = GenSym.p_dot;
                parser_getc(p);
                return;
            } else if (isalpha(*p.point)) {
                /* 2017-01-07 */
                p.sym = GenSym.bmemb;
                parser_getc(p);
                return;
            }
            parser_getc(p);
            if (p.ch == '*') {
                p.sym = GenSym.b_dotmult;
                parser_getc(p);
                return;
            } else if (p.ch == '/') {
                p.sym = GenSym.b_dotdiv;
                parser_getc(p);
                return;
            } else if (p.ch == '^') {
                p.sym = GenSym.b_dotpow;
                parser_getc(p);
                return;
            } else if (p.ch == '+') {
                p.sym = GenSym.b_dotadd;
                parser_getc(p);
                return;
            } else if (p.ch == '-') {
                p.sym = GenSym.b_dotsub;
                parser_getc(p);
                return;
            } else if (p.ch == '=') {
                p.sym = GenSym.b_doteq;
                parser_getc(p);
                return;
            } else if (p.ch == '>') {
                p.sym = GenSym.b_dotgt;
                parser_getc(p);
                if (p.ch == '=') {
                    p.sym = GenSym.b_dotgte;
                    parser_getc(p);
                }
                return;
            } else if (p.ch == '<') {
                p.sym = GenSym.b_dotlt;
                parser_getc(p);
                if (p.ch == '=') {
                    p.sym = GenSym.b_dotlte;
                    parser_getc(p);
                }
                return;
            } else if (p.ch == '!' && *p.point == '=') {
                p.sym = GenSym.b_dotneq;
                parser_getc(p);
                parser_getc(p);
                return;
            } else if (p.ch == '.') {
                p.sym = GenSym.b_ellip;
                parser_getc(p);
                return;
            } else {
                /* not a "dot operator", so back up */
                parser_ungetc(p);
            }
            goto default;   /* falls through */
        default:
            if (defining_list(p) && lag_range_sym(p)) {
                p.sym = GenSym.b_range;
                parser_getc(p);
                parser_getc(p);
                return;
            }
            if (defining_list(p) && !doing_genseries &&
                (bare_data_type(p.sym) || closing_sym(p.sym) ||
                 p.sym == GenSym.lag) && *(p.point - 2) == ' ') {
                /* may be forming a list, but only if there are
                   spaces between the terms
                */
                p.sym = GenSym.b_lcat;
                return;
            }
            if (isdigit(p.ch)) {
                parse_number(p);
                return;
            } else if (p.ch == '.' && isdigit(*p.point)) {
                parse_number(p);
                return;
            } else if (islower(p.ch) || isupper(p.ch) ||
                       word_start_special(p.ch)) {
                getword(p, 0);
                return;
            } else if (p.ch == '"') {
                p.idstr = get_quoted_string(p, prevsyms[0]);
                return;
            } else if (p.ch == '#' && prevsyms[0] == GenSym.cstr) {
                /* 2023-12-22: an uncaught inline comment? */
                p.ch = '\0';
                p.sym = GenSym.eot;
                return;
            } else {
                parser_print_input(p);
                if (isprint(p.ch)) {
                    pprintf(p.prn, _("Invalid character '%c'\n").ptr, p.ch);
                } else {
                    pprintf(p.prn, _("Unexpected byte 0x%x\n").ptr,
                            cast(ubyte) p.ch);
                }
                p.err = GretlError.parse;
                return;
            }
        } /* end ch switch */
    } /* end while ch != 0 */
}

extern (C) const(char) * getsymb_full(int t, const(Parser) * p)
{
    if (t == GenFunc.f_defmat) {
        /* pseudo function */
        return "DEFMAT";
    } else if ((t > GenFunc.f1_min && t < GenFunc.f1_max) ||
        (t > GenFunc.f1_max && t < GenFunc.f2_max) ||
        (t > GenFunc.f2_max && t < GenFunc.fn_max)) {
        return funname(t);
    }

    if (t == GenSym.eot) {
        return "EOT";
    }

    /* yes, well */
    if (t == GenSym.obs) {
        return "OBS";
    } else if (t == GenSym.osl) {
        return "OSL";
    } else if (t == GenSym.sub_addr) {
        return "SUB_ADDR";
    } else if (t == GenSym.dmstr) {
        return "DMSTR";
    } else if (t == GenSym.slraw) {
        return "SLRAW";
    } else if (t == GenSym.mspec) {
        return "MSPEC";
    } else if (t == GenSym.subsl) {
        return "SUBSL";
    } else if (t == GenSym.fargs) {
        return "FARGS";
    } else if (t == GenSym.list || t == GenSym.wlist) {
        return "LIST";
    } else if (t == GenSym.empty) {
        return "EMPTY";
    } else if (t == GenSym.listvar) {
        return "LISTVAR";
    } else if (t == GenSym.bmemb) {
        return "BMEMB";
    } else if (t == GenSym.series) {
        return "SERIES";
    } else if (t == GenSym.mat) {
        return "MAT";
    } else if (t == GenSym.undef) {
        return "UNDEF";
    } else if (t == GenSym.num) {
        return "NUM";
    } else if (t == GenSym.cnum) {
        return "CNUM";
    } else if (t == GenSym.ivec) {
        return "IVEC";
    } else if (t == GenSym.num_p) {
        return "NUM_P";
    } else if (t == GenSym.num_m) {
        return "NUM_M";
    } else if (t == GenSym.dbundle) {
        return "DBUNDLE";
    } else if (t == GenSym.dbmemb) {
        return "DBMEMB";
    } else if (t == GenSym.mmemb) {
        return "MMEMB";
    } else if (t == DumId.empty) {
        return "empty";
    }

    if (p != null) {
        if (t == GenSym.bundle) {
            return p.idstr;
        } else if (t == GenSym.array) {
            return p.idstr;
        } else if (t == GenSym.uobj) {
            return p.idstr;
        } else if (t == GenSym.con) {
            return constname(p.idnum);
        } else if (t == GenSym.dum) {
            return dumname(p.idnum);
        } else if (t == GenSym.dvar) {
            return dvarname(p.idnum);
        } else if (t == GenSym.mvar) {
            return mvarname(p.idnum);
        } else if (t == GenSym.ufun || t == GenSym.rfun) {
            return p.idstr;
        } else if (t == GenSym.str || t == GenSym.cstr) {
            return p.idstr;
        }
    } else {
        if (t == GenSym.bundle) {
            return "BUNDLE";
        } else if (t == GenSym.array) {
            return "ARRAY";
        } else if (t == GenSym.uobj) {
            return "UOBJ";
        } else if (t == GenSym.con) {
            return "CON";
        } else if (t == GenSym.dum) {
            return "dummy constant";
        } else if (t == GenSym.dvar) {
            return "DVAR";
        } else if (t == GenSym.mvar) {
            return "MVAR";
        } else if (t == GenSym.ufun) {
            return "UFUN";
        } else if (t == GenSym.rfun) {
            return "RFUN";
        } else if (t == GenSym.str) {
            return "STR";
        } else if (t == GenSym.cstr) {
            return "CSTR";
        }
    }

    switch (t) {
    case GenSym.b_asn:
        return "=";
    case GenSym.b_add:
    case GenSym.u_pos:
        return "+";
    case GenSym.b_sub:
    case GenSym.u_neg:
        return "-";
    case GenSym.b_mul:
        return "*";
    case GenSym.b_trmul:
        return "'";
    case GenSym.b_div:
        return "/";
    case GenSym.b_ldiv:
        return "\\";
    case GenSym.b_mod:
        return "%";
    case GenSym.b_pow:
        return "^";
    case GenSym.b_eq:
        return "==";
    case GenSym.b_neq:
        return "!=";
    case GenSym.b_gt:
        return ">";
    case GenSym.b_lt:
        return "<";
    case GenSym.b_gte:
        return ">=";
    case GenSym.b_lte:
        return "<=";
    case GenSym.b_and:
        return "&&";
    case GenSym.b_join:
        return "JOIN";
    case GenSym.b_range:
        return " to ";
    case GenSym.b_ellip:
        return "..";
    case GenSym.u_addr:
        return "&";
    case GenSym.b_or:
        return "||";
    case GenSym.u_not:
        return "!";
    case GenSym.g_lpr:
        return "(";
    case GenSym.g_rpr:
        return ")";
    case GenSym.g_lbr:
        return "[";
    case GenSym.g_rbr:
        return "]";
    case GenSym.g_lcb:
        return "{";
    case GenSym.g_rcb:
        return "}";
    case GenSym.b_dotmult:
        return ".*";
    case GenSym.b_dotdiv:
        return "./";
    case GenSym.b_dotpow:
        return ".^";
    case GenSym.b_dotadd:
        return ".+";
    case GenSym.b_dotsub:
        return ".-";
    case GenSym.b_doteq:
        return ".=";
    case GenSym.b_dotgt:
        return ".>";
    case GenSym.b_dotlt:
        return ".<";
    case GenSym.b_dotgte:
        return ".>=";
    case GenSym.b_dotlte:
        return ".<=";
    case GenSym.b_dotneq:
        return ".!=";
    case GenSym.b_dotasn:
        return "dot-assign";
    case GenSym.b_kron:
        return "**";
    case GenSym.b_hcat:
        return "~";
    case GenSym.b_vcat:
        return "|";
    case GenSym.b_lcat:
        return "LCAT";
    case GenSym.p_com:
        return ",";
    case GenSym.p_dot:
        return ".";
    case GenSym.p_semi:
        return ";";
    case GenSym.p_col:
        return ":";
    case GenSym.query:
        return "query";
    case GenSym.lag:
        return "lag";
    default:
        break;
    }

    return "unknown";
}

extern (C) const(char) * getsymb(int t)
{
    return getsymb_full(t, null);
}
