/* Port of libset.c */

module gretl.libset;


import core.stdc.ctype : isalpha;
import core.stdc.errno : ERANGE, errno;
import core.stdc.math : fabs, floor, log2, pow;
import core.stdc.stdio : FILE, fclose, fgets, fprintf, sprintf, sscanf, stderr;
import core.stdc.stdlib : atoi, calloc, free, getenv, labs, malloc, realloc, strtod, strtol, strtoull;
import core.stdc.string : memset, strchr, strcmp, strcpy, strlen, strncat, strncmp;
import core.stdc.time : time;
import std.algorithm.comparison : min;

import gretl.ghashtable;
import gretl.glib;
import gretl.gptrarray;
import gretl.gstrfuncs;
import gretl.gretl_commandsh;
import gretl.gretl_errorsh;
import gretl.gretl_matrix : get_blas_mnk_min, get_simd_k_max, get_simd_mn_min,
    gretl_matrix_copy, gretl_matrix_free, set_blas_mnk_min, set_simd_k_max,
    set_simd_mn_min;
import gretl.gretl_matrixh;
import gretl.gretl_modelh : HACKernel, MLVCVType;
import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_new_with_filename,
    pprintf, pputc, pputs;
import gretl.libgretlh;
import gretl.libseth : HACopts, OptimCode, PanelRobust, SetKey, ShowActivityFunc, VECMnorm;
import gretl.missingh;
import gretl.strutils : gretl_strdup, gretl_strndup, string_is_blank, tailstrip;
import gretl.options : get_optval_string, incompatible_options;
import gretl.gretl_errors : gretl_errmsg_set, gretl_errmsg_set_from_errno, gretl_errmsg_sprintf;
import gretl.gretl_utils : gretl_stopwatch;
import gretl.gretl_paths : gretl_addpath, gretl_fopen, gretl_maybe_prepend_dir, gretl_set_path_by_name, gretl_test_fopen, gretl_workdir;
import gretl.gretl_string_table : gretl_insert_builtin_string;
import gretl.gretl_intl : gretl_pop_c_numeric_locale, gretl_push_c_numeric_locale;
import gretl.genmain : current_series_index;
import gretl.matrix_extra : gretl_matrix_print_to_prn;
import gretl.printout : get_gretl_digits, set_gretl_digits;
import gretl.nls : get_default_nls_toler;
import gretl.graphing : get_plotstyle, set_plotstyle;
import gretl.uservar : get_scalar_value_by_name;
import gretl.usermat : get_matrix_by_name;
import gretl.random : gretl_rand_get_seed, gretl_rand_set_seed;

/* ---- declarations of not-yet-ported items used below ----
   Replace each with an import when its defining file is ported.
   NB: the headers for most of these (usermat.h, uservar.h, matrix_extra.h,
   gretl_func.h, gretl_string_table.h, gretl_mt.h, gretl_foreign.h, ...) were
   not provided, so the signatures are inferred from the call sites. */

extern (C) {
    int gretl_function_depth();
    int gretl_get_omp_threads();
    int gretl_set_omp_threads(int n);
}

/* ---- end of declarations ---- */

private enum AutoLag {
    stock_watson,
    wooldridge,
    neweywest
}

private enum SetCategory {
    behave = 1,
    numeric,
    rng,
    robust,
    ts,
    special
}

private enum SVType {
    all,
    int_,
    double_
}

/* for values that really want a non-negative integer */
private enum UNSET_INT = -9;

private bool is_unset(int i)
{
    return (i == UNSET_INT);
}

struct SetState {
    int flags;
    /* small integer values */
    gint8 optim;                /* code for preferred optimizer */
    gint8 vecm_norm;            /* VECM beta normalization */
    gint8 garch_vcv;            /* GARCH vcv variant */
    gint8 garch_alt_vcv;        /* GARCH vcv variant, robust estimation */
    gint8 arma_vcv;             /* ARMA vcv variant */
    gint8 wildboot_d;           /* distribution for wild bootstrap */
    gint8 fdjac_qual;           /* quality of "fdjac" function */
    gint8 use_qr;               /* off, on or pivot */
    gint8 max_verbose;          /* optimizer verbosity level */
    gint8 hc_version;           /* HCCME version */
    gint8 panel_robust;         /* panel robust vcv estimator */
    gint8 hac_kernel;           /* HAC kernel type */
    gint8 auto_hac_lag;         /* HAC automatic lag-length formula */
    gint8 user_hac_lag;         /* fixed user-set HAC lag length */
    gint8 lbfgs_mem;            /* memory span for L-BFGS-B */
    gint8 quantile_type;        /* Formula for computing quantiles */
    /* potentially larger integers */
    int horizon;                /* for VAR impulse responses */
    int bootrep;                /* bootstrap replications */
    int loop_maxiter;           /* max no. of iterations in non-for loops */
    int bfgs_maxiter;           /* max iterations, BFGS */
    int bfgs_verbskip;          /* BFGS: show one in n iterations  */
    int boot_iters;             /* max iterations, IRF bootstrap */
    int bhhh_maxiter;           /* max iterations, BHHH */
    int rq_maxiter;             /* max iterations for quantreg, simplex */
    int gmm_maxiter;            /* max iterations for iterated GMM */
    /* floating-point values */
    double conv_huge;           /* conventional value for $huge */
    double nls_toler;           /* NLS convergence criterion */
    double bfgs_toler;          /* convergence tolerance, BFGS */
    double bfgs_maxgrad;        /* max acceptable gradient norm, BFGS */
    double bhhh_toler;          /* convergence tolerance, BHHH */
    double qs_bandwidth;        /* bandwidth for QS HAC kernel */
    double nadarwat_trim;       /* multiple of h to use in nadarwat() for trimming */
    /* strings */
    char[8] csv_write_na;       /* representation of NA in CSV output */
    char[8] csv_read_na;        /* representation of NA in CSV input */
    /* matrices */
    GretlMatrixRaw * initvals;     /* parameter initializer */
    GretlMatrixRaw * initcurv;     /* initial curvature matrix for BFGS */
    GretlMatrixRaw * matmask;      /* mask for series -> matrix conversion */
}

private struct GlobalVars {
    gint8 gretl_debug;
    gint8 gretl_assert;
    gint8 datacols;
    gint8 plot_collect;
    gint8 R_functions;
    gint8 R_lib;
    gint8 loglevel;
    gint8 logstamp;
    gint8 csv_digits;
    gint8 hac_missvals;
    int gmp_bits;
}

private __gshared GlobalVars globals = GlobalVars(0, 0, 5, 0, 0, 1, 2, 0, UNSET_INT, HACopts.es, 256);

/* globals for internal use */
private __gshared int seed_is_set;
private __gshared int comments_on;
private __gshared char data_delim = ',';
private __gshared char export_decpoint = '.';
private __gshared int suppress_plots;

struct SetVar {
    SetKey key;       /* internal integer key */
    const(char) * name; /* userspace name */
    gint8 category;   /* for printing purposes */
    size_t offset;    /* byte offset into state or globals struct,
                         where applicable */
}

private enum BEHAVE = SetCategory.behave;
private enum NUMERIC = SetCategory.numeric;
private enum RNG = SetCategory.rng;
private enum ROBUST = SetCategory.robust;
private enum TS = SetCategory.ts;
private enum SPECIAL = SetCategory.special;

__gshared SetVar[85] setvars = [
    /* booleans (bitflags) */
    SetVar(SetKey.use_cwd,      "use_cwd",   BEHAVE),
    SetVar(SetKey.echo_on,      "echo",      BEHAVE),
    SetVar(SetKey.msgs_on,      "messages",  BEHAVE),
    SetVar(SetKey.force_decpoint, "force_decpoint", BEHAVE),
    SetVar(SetKey.use_svd,      "svd",       NUMERIC),
    SetVar(SetKey.prewhiten,    "hac_prewhiten", ROBUST),
    SetVar(SetKey.force_hc,     "force_hc",      ROBUST),
    SetVar(SetKey.use_lbfgs,    "lbfgs",        NUMERIC),
    SetVar(SetKey.shell_ok,     "shell_ok",     SPECIAL),
    SetVar(SetKey.warnings,     "warnings",     BEHAVE),
    SetVar(SetKey.skip_missing, "skip_missing", BEHAVE),
    SetVar(SetKey.bfgs_rstep,   "bfgs_richardson", NUMERIC),
    SetVar(SetKey.robust_z,     "robust_z", ROBUST),
    SetVar(SetKey.mwrite_g,     "mwrite_g", BEHAVE),
    SetVar(SetKey.mpi_use_smt,  "mpi_use_smt", BEHAVE),
    SetVar(SetKey.state_flag_max, null),
    /* small integers */
    SetVar(SetKey.gretl_optim,  "optimizer", NUMERIC, SetState.optim.offsetof),
    SetVar(SetKey.vecm_norm,    "vecm_norm", TS,      SetState.vecm_norm.offsetof),
    SetVar(SetKey.garch_vcv,    "garch_vcv", ROBUST,  SetState.garch_vcv.offsetof),
    SetVar(SetKey.garch_alt_vcv, "garch_alt_vcv", ROBUST, SetState.garch_alt_vcv.offsetof),
    SetVar(SetKey.arma_vcv,      "arma_vcv", ROBUST, SetState.arma_vcv.offsetof),
    SetVar(SetKey.wildboot_dist, "wildboot", BEHAVE, SetState.wildboot_d.offsetof),
    SetVar(SetKey.fdjac_qual,    "fdjac_quality", NUMERIC, SetState.fdjac_qual.offsetof),
    SetVar(SetKey.use_qr,        "force_qr", NUMERIC, SetState.use_qr.offsetof),
    SetVar(SetKey.max_verbose,   "max_verbose", BEHAVE, SetState.max_verbose.offsetof),
    SetVar(SetKey.hc_version,    "hc_version",  ROBUST, SetState.hc_version.offsetof),
    SetVar(SetKey.panel_robust,  "panel_robust", ROBUST, SetState.panel_robust.offsetof),
    SetVar(SetKey.hac_kernel,    "hac_kernel",  ROBUST, SetState.hac_kernel.offsetof),
    SetVar(SetKey.hac_lag,       "hac_lag",     ROBUST),
    SetVar(SetKey.user_hac_lag,  null),
    SetVar(SetKey.lbfgs_mem,     "lbfgs_mem",     NUMERIC, SetState.lbfgs_mem.offsetof),
    SetVar(SetKey.quantile_type, "quantile_type", BEHAVE, SetState.quantile_type.offsetof),
    SetVar(SetKey.state_small_int_max, null),
    /* larger integers */
    SetVar(SetKey.horizon,       "horizon", TS,     SetState.horizon.offsetof),
    SetVar(SetKey.bootrep,       "bootrep", BEHAVE, SetState.bootrep.offsetof),
    SetVar(SetKey.loop_maxiter,  "loop_maxiter",  BEHAVE,  SetState.loop_maxiter.offsetof),
    SetVar(SetKey.bfgs_maxiter,  "bfgs_maxiter",  NUMERIC, SetState.bfgs_maxiter.offsetof),
    SetVar(SetKey.bfgs_verbskip, "bfgs_verbskip", BEHAVE,  SetState.bfgs_verbskip.offsetof),
    SetVar(SetKey.boot_iters,    "boot_iters",    TS,      SetState.boot_iters.offsetof),
    SetVar(SetKey.bhhh_maxiter,  "bhhh_maxiter",  NUMERIC, SetState.bhhh_maxiter.offsetof),
    SetVar(SetKey.rq_maxiter,    "rq_maxiter",    NUMERIC, SetState.rq_maxiter.offsetof),
    SetVar(SetKey.gmm_maxiter,   "gmm_maxiter",   NUMERIC, SetState.gmm_maxiter.offsetof),
    SetVar(SetKey.state_int_max, null),
    /* doubles */
    SetVar(SetKey.conv_huge,     "huge",         BEHAVE,  SetState.conv_huge.offsetof),
    SetVar(SetKey.nls_toler,     "nls_toler",    NUMERIC, SetState.nls_toler.offsetof),
    SetVar(SetKey.bfgs_toler,    "bfgs_toler",   NUMERIC, SetState.bfgs_toler.offsetof),
    SetVar(SetKey.bfgs_maxgrad,  "bfgs_maxgrad", NUMERIC, SetState.bfgs_maxgrad.offsetof),
    SetVar(SetKey.bhhh_toler,    "bhhh_toler",   NUMERIC, SetState.bhhh_toler.offsetof),
    SetVar(SetKey.qs_bandwidth,  "qs_bandwidth", ROBUST,  SetState.qs_bandwidth.offsetof),
    SetVar(SetKey.nadarwat_trim, "nadarwat_trim", NUMERIC, SetState.nadarwat_trim.offsetof),
    SetVar(SetKey.state_float_max, null),
    /* strings */
    SetVar(SetKey.csv_write_na,  "csv_write_na", SPECIAL),
    SetVar(SetKey.csv_read_na,   "csv_read_na",  SPECIAL),
    /* matrices */
    SetVar(SetKey.initvals,      "initvals",    NUMERIC),
    SetVar(SetKey.initcurv,      "initcurv",    NUMERIC),
    SetVar(SetKey.matmask,       "matrix_mask", BEHAVE),
    SetVar(SetKey.state_vars_max, null),
    /* global ints */
    SetVar(SetKey.gretl_debug,   "debug",     BEHAVE, GlobalVars.gretl_debug.offsetof),
    SetVar(SetKey.gretl_assert,  "assert",    BEHAVE, GlobalVars.gretl_assert.offsetof),
    SetVar(SetKey.datacols,      "datacols",  BEHAVE, GlobalVars.datacols.offsetof),
    SetVar(SetKey.plot_collect,  "plot_collection", BEHAVE, GlobalVars.plot_collect.offsetof),
    SetVar(SetKey.r_functions,   "R_functions", BEHAVE, GlobalVars.R_functions.offsetof),
    SetVar(SetKey.r_lib,         "R_lib",       BEHAVE, GlobalVars.R_lib.offsetof),
    SetVar(SetKey.loglevel,      "loglevel",    BEHAVE, GlobalVars.loglevel.offsetof),
    SetVar(SetKey.logstamp,      "logstamp",    BEHAVE, GlobalVars.logstamp.offsetof),
    SetVar(SetKey.csv_digits,    "csv_digits",  BEHAVE, GlobalVars.csv_digits.offsetof),
    SetVar(SetKey.hac_missvals,  "hac_missvals", BEHAVE, GlobalVars.hac_missvals.offsetof),
    SetVar(SetKey.ns_small_int_max, null),
    SetVar(SetKey.gmp_bits,      "gmp_bits",    BEHAVE, GlobalVars.gmp_bits.offsetof),
    SetVar(SetKey.ns_max, null),
    /* delegated ints */
    SetVar(SetKey.blas_mnk_min,  "blas_mnk_min", BEHAVE),
    SetVar(SetKey.omp_n_threads, "omp_num_threads", SPECIAL),
    SetVar(SetKey.simd_k_max,    "simd_k_max",  BEHAVE),
    SetVar(SetKey.simd_mn_min,   "simd_mn_min", BEHAVE),
    /* specials */
    SetVar(SetKey.seed,          "seed",      RNG),
    SetVar(SetKey.csv_delim,     "csv_delim", SPECIAL),
    SetVar(SetKey.stopwatch,     "stopwatch", SPECIAL),
    SetVar(SetKey.verbose,       "verbose",   SPECIAL),
    SetVar(SetKey.sv_workdir,    "workdir",   SPECIAL),
    SetVar(SetKey.sv_logfile,    "logfile",   SPECIAL),
    SetVar(SetKey.graph_theme,   "graph_theme", SPECIAL),
    SetVar(SetKey.disp_digits,   "display_digits", SPECIAL),
    SetVar(SetKey.tex_plot_opts, "tex_plot_opts", SPECIAL)
];

private bool libset_boolvar(SetKey k)
{
    return (k < SetKey.state_flag_max || k == SetKey.r_functions ||
            k == SetKey.r_lib || k == SetKey.logstamp);
}

private bool libset_double(SetKey k)
{
    return (k > SetKey.state_int_max && k < SetKey.state_float_max);
}

private bool libset_int(SetKey k)
{
    return ((k > SetKey.state_flag_max && k < SetKey.state_int_max) ||
            (k > SetKey.state_vars_max && k < SetKey.ns_int_max));
}

private bool libset_small_int(SetKey k)
{
    return (k < SetKey.state_small_int_max ||
            (k > SetKey.state_vars_max && k < SetKey.ns_small_int_max));
}

private bool coded_intvar(SetKey k)
{
    return (k == SetKey.garch_vcv ||
            k == SetKey.garch_alt_vcv ||
            k == SetKey.arma_vcv ||
            k == SetKey.hac_lag ||
            k == SetKey.hac_kernel ||
            k == SetKey.hc_version ||
            k == SetKey.panel_robust ||
            k == SetKey.use_qr ||
            k == SetKey.vecm_norm ||
            k == SetKey.gretl_optim ||
            k == SetKey.max_verbose ||
            k == SetKey.wildboot_dist ||
            k == SetKey.quantile_type ||
            k == SetKey.gretl_assert ||
            k == SetKey.plot_collect ||
            k == SetKey.loglevel ||
            k == SetKey.hac_missvals);
}

/* the current set of state variables */
__gshared SetState * state;

/* In case we have a SetKey and want to print the
   associated userspace keyword. */

private const(char) * setkey_get_name(SetKey key)
{
    int i;

    for (i=0; i<setvars.length; i++) {
        if (setvars[i].key == key) {
            return setvars[i].name;
        }
    }
    return null;
}

/* Set up a hash table to map from userspace keywords
   to SetVar structs.
*/

private GHashTable * libset_hash_init()
{
    GHashTable * ht = g_hash_table_new(&g_str_hash, &g_str_equal);
    int i;

    for (i=0; i<setvars.length; i++) {
        if (setvars[i].name != null) {
            g_hash_table_insert(ht, cast(gpointer) setvars[i].name, &setvars[i]);
        }
    }

    return ht;
}

private __gshared GHashTable * svht;

private void libset_hash_cleanup()
{
    if (svht != null) {
        g_hash_table_destroy(svht);
        svht = null;
    }
}

private SetVar * get_setvar_by_name(const(char) * name)
{
    SetVar * ret = null;

    if (svht == null) {
        svht = libset_hash_init();
    }

    ret = cast(SetVar *) g_hash_table_lookup(svht, name);

    if (ret == null) {
        /* backward compatibility */
        if (!strcmp(name, "csv_na")) {
            ret = cast(SetVar *) g_hash_table_lookup(svht, "csv_write_na".ptr);
        }
    }

    return ret;
}

private void * setvar_get_target(SetVar * sv)
{
    void * p;

    if (sv.offset == 0 || sv.key > SetKey.gmp_bits) {
        /* FIXME criterion? (0 offset is legit for "debug") */
        if (sv.key != SetKey.gretl_debug) {
            return null;
        }
    }

    p = (sv.key < SetKey.state_vars_max) ? cast(void *) state : cast(void *) &globals;
    return (cast(ubyte *) p) + sv.offset;
}

private int INTS_OFFSET()
{
    return cast(int) (1 + log2(cast(double) SetKey.state_flag_max));
}

private void * setkey_get_target(SetKey key, SVType t)
{
    int i = INTS_OFFSET() + key - SetKey.gretl_optim;
    SetVar * sv = &setvars[i];

    if (sv.key != key) {
        fprintf(stderr, "*** internal error, looking for %s, found %s ***\n",
                setkey_get_name(key), sv.name);
        return null;
    } else if ((t == SVType.int_ && !libset_int(key)) ||
               (t == SVType.double_ && !libset_double(key))) {
        fprintf(stderr, "*** type mismatch in setkey_get_target for %s ***\n",
                sv.name);
        return null;
    } else {
        return setvar_get_target(sv);
    }
}

/* value strings for integer-coded variables */

private immutable(char) *[7] gvc_strs = ["unset", "hessian", "im", "op", "qml", "bw", null];
private immutable(char) *[3] gvr_strs = ["qml", "bw", null];
private immutable(char) *[3] avc_strs = ["hessian", "op", null];
private immutable(char) *[4] hkn_strs = ["bartlett", "parzen", "qs", null];
private immutable(char) *[6] hcv_strs = ["0", "1", "2", "3", "3a", null];
private immutable(char) *[5] vnm_strs = ["phillips", "diag", "first", "none", null];
private immutable(char) *[4] opt_strs = ["auto", "BFGS", "newton", null];
private immutable(char) *[4] mxv_strs = ["off", "on", "full", null];
private immutable(char) *[3] wbt_strs = ["rademacher", "mammen", null];
private immutable(char) *[4] qnt_strs = ["Q6", "Q7", "Q8", null];
private immutable(char) *[5] ast_strs = ["off", "warn", "stop", "fatal", null];
private immutable(char) *[4] plc_strs = ["off", "auto", "on", null];
private immutable(char) *[6] csv_strs = ["comma", "space", "tab", "semicolon", "pipe", null];
private immutable(char) *[4] ahl_strs = ["nw1", "nw2", "nw3", null];
private immutable(char) *[6] llv_strs = ["debug", "info", "warn", "error", "critical", null];
private immutable(char) *[4] qrp_strs = ["off", "on", "pivot", null];
private immutable(char) *[4] hmv_strs = ["off", "es", "am", null];
private immutable(char) *[4] pnr_strs = ["arellano", "pcse", "scc", null];

struct CodevarInfo {
    SetKey key;
    const(char *) * strvals;
}

/* look-up table for sets of value strings */

__gshared CodevarInfo[18] coded;

shared static this()
{
    coded = [
        CodevarInfo(SetKey.garch_vcv,     gvc_strs.ptr),
        CodevarInfo(SetKey.garch_alt_vcv, gvr_strs.ptr),
        CodevarInfo(SetKey.arma_vcv,      avc_strs.ptr),
        CodevarInfo(SetKey.hac_kernel,    hkn_strs.ptr),
        CodevarInfo(SetKey.hc_version,    hcv_strs.ptr),
        CodevarInfo(SetKey.panel_robust,  pnr_strs.ptr),
        CodevarInfo(SetKey.vecm_norm,     vnm_strs.ptr),
        CodevarInfo(SetKey.gretl_optim,   opt_strs.ptr),
        CodevarInfo(SetKey.max_verbose,   mxv_strs.ptr),
        CodevarInfo(SetKey.wildboot_dist, wbt_strs.ptr),
        CodevarInfo(SetKey.csv_delim,     csv_strs.ptr),
        CodevarInfo(SetKey.quantile_type, qnt_strs.ptr),
        CodevarInfo(SetKey.gretl_assert,  ast_strs.ptr),
        CodevarInfo(SetKey.plot_collect,  plc_strs.ptr),
        CodevarInfo(SetKey.hac_lag,       ahl_strs.ptr),
        CodevarInfo(SetKey.loglevel,      llv_strs.ptr),
        CodevarInfo(SetKey.use_qr,        qrp_strs.ptr),
        CodevarInfo(SetKey.hac_missvals,  hmv_strs.ptr)
    ];
}

private const(char *) * libset_option_strings(SetKey key)
{
    int i;

    for (i=0; i<coded.length; i++) {
        if (coded[i].key == key) {
            return coded[i].strvals;
        }
    }
    return null;
}

private void coded_var_show_opts(SetKey key, Prn * prn)
{
    const(char *) * S = libset_option_strings(key);

    if (S != null) {
        pputs(prn, "valid settings:");
        while (*S != null) {
            pprintf(prn, " %s", *S);
            S++;
        }
        pputc(prn, '\n');
    }
}

private const(char) * garch_alt_vcv_string()
{
    if (state.garch_alt_vcv == MLVCVType.qml) {
        return gvr_strs[0];
    } else if (state.garch_alt_vcv == MLVCVType.bw) {
        return gvr_strs[1];
    } else {
        return "unset";
    }
}

private const(char) * arma_vcv_string()
{
    if (state.arma_vcv == MLVCVType.hessian) {
        return avc_strs[0];
    } else if (state.arma_vcv == MLVCVType.op) {
        return avc_strs[1];
    } else {
        return "unset";
    }
}

private const(char) * libset_option_string(SetKey key)
{
    if (key == SetKey.hac_lag) {
        return hac_lag_string();          /* special */
    } else if (key == SetKey.garch_alt_vcv) {
        return garch_alt_vcv_string();    /* special */
    } else if (key == SetKey.arma_vcv) {
        return arma_vcv_string();         /* special */
    } else {
        const(char *) * strs = libset_option_strings(key);
        void * valp = setkey_get_target(key, SVType.int_);

        return strs[*cast(gint8 *) valp];
    }
}

private void print_state_matrix(SetKey key, Prn * prn, Gretlopt opt)
{
    GretlMatrixRaw * m;
    const(char) * name;

    if (key == SetKey.initvals) {
        name = "initvals";
        m = state.initvals;
    } else {
        name = "initcurv";
        m = state.initcurv;
    }

    if (opt & Gretlopt.OPT_D) {
        if (m == null) {
            pprintf(prn, " %s = auto\n", name);
        } else {
            gretl_matrix_print_to_prn(m, name, prn);
        }
    }
}

/* check_for_state() returns non-zero if the program options
   state is not readable */

private int check_for_state()
{
    if (state == null) {
        return libset_init();
    } else {
        return 0;
    }
}

private int flag_to_bool(SetState * sv, SetKey key)
{
    if (sv == null) {
        return 0;
    } else {
        return (sv.flags & key) ? 1 : 0;
    }
}

private void state_vars_copy(SetState * sv)
{
    /* copy everything */
    *sv = *state;
    /* but set matrix pointers to null */
    sv.initvals = null;
    sv.initcurv = null;
    sv.matmask = null;
}

private __gshared SetState default_state = SetState(
    SetKey.echo_on | SetKey.msgs_on | SetKey.warnings | SetKey.skip_missing, /* .flags */
    OptimCode.auto_,     /* .optim */
    VECMnorm.phillips,   /* .vecm_norm */
    MLVCVType.unset,     /* .garch_vcv */
    MLVCVType.unset,     /* .garch_alt_vcv */
    MLVCVType.hessian,   /* .arma_vcv */
    0,              /* .wildboot_dist */
    0,              /* .fdjac_qual */
    0,              /* .use_qr */
    0,              /* .max_verbose */
    1,              /* .hc_version */
    0,              /* .panel_robust */
    HACKernel.bartlett,       /* .hac_kernel */
    AutoLag.stock_watson,     /* .auto_hac_lag */
    UNSET_INT,             /* .user_hac_lag */
    8,             /* .lbfgs_mem */
    0,             /* .quantile_type */
    UNSET_INT,     /* .horizon */
    1000,          /* .bootrep */
    100000,        /* .loop_maxiter */
    UNSET_INT,     /* .bfgs_maxiter */
    1,       /* .bfgs_verbskip */
    1999,    /* .boot_iters */
    500,     /* .bhhh_maxiter */
    1000,    /* .rq_maxiter */
    250,     /* .gmm_maxiter */
    1.0e100, /* .conv_huge */
    NADBL,   /* .nls_toler */
    NADBL,   /* .bfgs_toler */
    5.0,     /* .bfgs_maxgrad */
    1.0e-6,  /* .bhhh_toler */
    2.0,     /* .qs_bandwidth */
    4.0,     /* .nadarwat_trim */
    "NA\0\0\0\0\0\0",      /* .csv_write_na */
    "default\0", /* .csv_read_na */
    null,  /* .initvals */
    null,  /* .initcurv */
    null   /* .matmask */
);

private void state_vars_init(SetState * sv)
{
    *sv = default_state;
}

int get_bkbp_k(const(DatasetRaw) * dset)
{
    if (dset.pd == 1) {
        return 3;
    } else if (dset.pd == 4) {
        return 12;
    } else if (dset.pd == 12) {
        return 36;
    } else {
        return 3;
    }
}

void get_bkbp_periods(const(DatasetRaw) * dset, int * l, int * u)
{
    *l = (dset.pd == 4) ? 6 :
        (dset.pd == 12) ? 18 : 2;

    *u = (dset.pd == 4) ? 32 :
        (dset.pd == 12) ? 96 : 8;
}

void set_gretl_echo(int e)
{
    if (check_for_state()) return;

    if (e) {
        state.flags |= SetKey.echo_on;
    } else {
        state.flags &= ~SetKey.echo_on;
    }
}

int gretl_echo_on()
{
    if (check_for_state()) {
        return 1;
    } else {
        return flag_to_bool(state, SetKey.echo_on);
    }
}

int gretl_comments_on()
{
    if (gretl_function_depth() > 0) {
        return 0;
    } else {
        return comments_on;
    }
}

void set_gretl_messages(int e)
{
    if (check_for_state()) return;

    if (e) {
        state.flags |= SetKey.msgs_on;
    } else {
        state.flags &= ~SetKey.msgs_on;
    }
}

int gretl_messages_on()
{
    if (check_for_state()) {
        return 1;
    } else {
        return flag_to_bool(state, SetKey.msgs_on);
    }
}

int gretl_warnings_on()
{
    if (check_for_state()) return 1;
    return flag_to_bool(state, SetKey.warnings);
}

int gretl_debugging_on()
{
    return globals.gretl_debug;
}

private enum DEFAULT_MP_BITS = 256;

private bool mp_bits_ok(int b)
{
    return (b >= 256 && b <= 8192);
}

/* Called from the mp_ols plugin, also gretl_utils.c and
   the GUI model specification dialog.
*/

int get_mp_bits()
{
    if (globals.gmp_bits > DEFAULT_MP_BITS) {
        return globals.gmp_bits;
    } else {
        char * s = getenv("GRETL_MP_BITS");
        int b;

        if (s != null) {
            b = atoi(s);
            if (mp_bits_ok(b)) {
                return b;
            }
        }
    }

    return DEFAULT_MP_BITS;
}

/* start accessors for libset matrices */

GretlMatrixRaw * get_initvals()
{
    GretlMatrixRaw * iv;

    /* note: we nullify initvals after first use */
    check_for_state();
    iv = state.initvals;
    state.initvals = null;
    return iv;
}

int n_initvals()
{
    check_for_state();
    if (state.initvals != null) {
        return gretl_vector_get_length(state.initvals);
    } else {
        return 0;
    }
}

GretlMatrixRaw * get_initcurv()
{
    GretlMatrixRaw * ic;

    /* note: like initvals, we nullify initcurv after first use */
    check_for_state();
    ic = state.initcurv;
    state.initcurv = null;
    return ic;
}

int n_initcurv()
{
    check_for_state();
    if (state.initcurv != null) {
        return gretl_vector_get_length(state.initcurv);
    } else {
        return 0;
    }
}

const(GretlMatrixRaw) * get_matrix_mask()
{
    check_for_state();
    return state.matmask;
}

int get_matrix_mask_nobs()
{
    int n = 0;

    check_for_state();

    if (state.matmask != null) {
        int i;

        for (i=0; i<state.matmask.rows; i++) {
            if (state.matmask.val[i] != 0.0) {
                n++;
            }
        }
    }

    return n;
}

/* end accessors for libset matrices */

int get_hac_lag(int T)
{
    int h = 0;

    check_for_state();

    /* Variants of Newey-West */

    if (state.user_hac_lag >= 0 && state.user_hac_lag < T - 2) {
        /* FIXME upper limit? */
        h = state.user_hac_lag;
    } else if (state.auto_hac_lag == AutoLag.wooldridge) {
        h = cast(int) floor(4 * pow(T/100.0, 2.0/9));
    } else {
        /* Stock-Watson default */
        h = cast(int) floor(0.75 * pow(cast(double) T, 1.0/3));
    }

    return h;
}

/* prewhitening implies nw3, but not vice versa */

int data_based_hac_bandwidth()
{
    if (is_unset(state.user_hac_lag)) {
        if (state.auto_hac_lag == AutoLag.neweywest ||
            (state.flags & SetKey.prewhiten)) {
            return 1;
        }
    }

    return 0;
}

private const(char) * hac_lag_string()
{
    check_for_state();

    if (state.user_hac_lag >= 0 && state.user_hac_lag < 127) {
        static __gshared char[6] lagstr;

        sprintf(lagstr.ptr, "%d", cast(int) state.user_hac_lag);
        return lagstr.ptr;
    } else {
        return ahl_strs[state.auto_hac_lag];
    }
}

/* set max lag for HAC estimation */

private int parse_hac_lag_variant(const(char) * s)
{
    int i, err = 0;

    for (i=0; ahl_strs[i] != null; i++) {
        if (!strcmp(s, ahl_strs[i])) {
            state.auto_hac_lag = cast(gint8) i;
            state.user_hac_lag = UNSET_INT;
            return 0;
        }
    }

    err = libset_get_scalar(SetKey.hac_lag, s, &i, null);
    if (!err) {
        state.user_hac_lag = cast(gint8) i;
    }

    return err;
}

private enum NumericStatus {
    numeric_ok,
    numeric_bad,
    non_numeric
}

/* Test @s for being a string representation of a numeric value
   in the C locale, using strtod() and/or strtol(). This can be
   used to retrieve a floating-point value (if @px is non-null)
   or an integer value (@pi non-null); exactly one of these
   pointers should be non-null.

   A return value of numeric_ok means that @s is indeed numeric
   and the converted value is within range for a double (if @px
   is non-null) or a 32-bit integer (if @pi is non-null).

   A return of numeric_bad means that @s is numeric but out of
   range for the target type.

   A return of non_numeric means that @s is not a numeric string;
   one may then proceed to test whether it's the name of a scalar
   variable.
*/

private NumericStatus
libset_numeric_test(const(char) * s, int * pi, double * px)
{
    NumericStatus ret = NumericStatus.numeric_ok;
    const(char) * test;

    if (!strcmp(s, "inf") || !strcmp(s, "nan")) {
        return NumericStatus.numeric_bad;
    } else if (isalpha(*s)) {
        return NumericStatus.non_numeric;
    }

    errno = 0;
    gretl_push_c_numeric_locale();

    if (px != null) {
        /* looking for a floating-point value */
        *px = strtod(s, &test);
        if (*test != '\0') {
            ret = NumericStatus.non_numeric;
        } else if (errno == ERANGE) {
            gretl_errmsg_set_from_errno(s, errno);
            ret = NumericStatus.numeric_bad;
        }
    } else {
        /* looking for an integer value */
        import core.stdc.config : c_long; c_long li = strtol(s, &test, 10);

        if (*test != '\0') {
            /* try for a floating-point value that's also a valid int? */
            const(char) * testx;
            double x;

            errno = 0;
            x = strtod(s, &testx);

            if (*testx != '\0') {
                ret = NumericStatus.non_numeric;
            } else if (errno == ERANGE) {
                ret = NumericStatus.numeric_bad;
            } else {
                /* numeric, but does it work as an int? */
                if (x == floor(x) && fabs(x) <= int.max) {
                    *pi = cast(int) x;
                } else {
                    ret = NumericStatus.numeric_bad;
                }
            }
        } else if (errno == ERANGE) {
            gretl_errmsg_set_from_errno(s, errno);
            ret = NumericStatus.numeric_bad;
        } else if (labs(li) > int.max) {
            /* OK as a long but too big for 32-bit int */
            ret = NumericStatus.numeric_bad;
        } else {
            *pi = cast(int) li;
        }
    }

    gretl_pop_c_numeric_locale();

    return ret;
}

private int negval_invalid(SetKey key)
{
    int ret = 1; /* presume invalid */

    if (key > 0) {
        if (key == SetKey.blas_mnk_min ||
            key == SetKey.simd_k_max || key == SetKey.simd_mn_min) {
            /* these can all be set to -1 */
            ret = 0;
        }
    }

    return ret;
}

private int libset_get_scalar(SetKey key, const(char) * arg,
                              int * pi, double * px)
{
    double x = NADBL;
    int err = 0;
    NumericStatus nstatus;

    if (arg == null || *arg == '\0') {
        return GretlError.args;
    }

    if (key == SetKey.omp_n_threads && !strcmp(arg, "default")) {
        /* (the OPENMP_BUILD variant, gretl_n_physical_cores(), is not part of this port) */
        *pi = 0;
        return 0;
    }

    nstatus = libset_numeric_test(arg, pi, px);

    if (nstatus == NumericStatus.numeric_bad) {
        return GretlError.invarg; /* handled */
    } else if (nstatus == NumericStatus.numeric_ok) {
        if (pi != null && negval_invalid(key) && *pi < 0) {
            err = GretlError.invarg;
        } else if (px != null && *px < 0.0) {
            err = GretlError.invarg;
        }
        return err; /* handled */
    }

    /* handle the non-numeric case */
    x = get_scalar_value_by_name(arg, &err);

    if (!err) {
        if (negval_invalid(key) && x < 0.0) {
            err = GretlError.invarg;
        } else if (px != null) {
            *px = x;
        } else if (pi != null) {
            if (na(x) || fabs(x) > cast(double) int.max) {
                err = GretlError.invarg;
            } else {
                *pi = cast(int) x;
            }
        }
    }

    return err;
}

/* Called only when setting the PRNG seed */

private int libset_get_seed(const(char) * arg, guint64 * pu, int * automatic)
{
    guint64 lu = 0;
    const(char) * test = null;
    double x = NADBL;
    int err = 0;

    if (!strcmp(arg, "auto")) {
        *automatic = 1;
        *pu = cast(guint64) time(null);
        return 0;
    }

    errno = 0;
    lu = strtoull(arg, &test, 10);

    if (*test == '\0') {
        if (errno) {
            return GretlError.invarg;
        } else {
            *pu = lu;
            return 0;
        }
    }

    /* fallback, when @arg is not numeric */
    x = get_scalar_value_by_name(arg, &err);
    if (err) {
        return err;
    }

    if (na(x) || x <= 0.0 || x > cast(double) ulong.max) {
        err = GretlError.invarg;
    } else {
        *pu = cast(guint64) floor(x);
    }

    return err;
}

private int n_strvals(const(char *) * s)
{
    int n = 0;

    while (*s != null) {
        n++; s++;
    }
    return n;
}

private int parse_libset_int_code(SetKey key, const(char) * val)
{
    int i, err = GretlError.data;

    if (key == SetKey.hac_lag) {
        err = parse_hac_lag_variant(val);
    } else if (coded_intvar(key)) {
        const(char *) * strs = libset_option_strings(key);
        void * valp = setkey_get_target(key, SVType.int_);
        int ival = -1;

        for (i=0; strs[i] != null; i++) {
            if (!g_ascii_strcasecmp(val, strs[i])) {
                ival = i;
                break;
            }
        }
        if (ival >= 0) {
            /* (the C source re-declares valp here, with the same value) */
            err = 0;
            if (key == SetKey.garch_alt_vcv) {
                ival = (ival == 1) ? MLVCVType.bw : MLVCVType.qml;
            } else if (key == SetKey.arma_vcv) {
                ival = (ival == 1) ? MLVCVType.op : MLVCVType.hessian;
            }
            *cast(gint8 *) valp = cast(gint8) ival;
        } else if (key == SetKey.max_verbose || key == SetKey.loglevel) {
            /* special: bare integers allowed? */
            int n = n_strvals(strs);

            for (i=0; i<n; i++) {
                if (val[0] == i + 48 && val[1] == '\0') {
                    *cast(gint8 *) valp = cast(gint8) i;
                    err = 0;
                    break;
                }
            }
        }
    }

    if (err) {
        gretl_errmsg_sprintf(_("%s: invalid value '%s'").ptr,
                             setkey_get_name(key), val);
    }

    return err;
}

/* start public functions called from gui/settings.c */

void set_xsect_hccme(const(char) * s)
{
    if (check_for_state()) return;

    if (!strncmp(s, "HC", 2)) {
        s += 2;
    }
    parse_libset_int_code(SetKey.hc_version, s);
}

void set_tseries_hccme(const(char) * s)
{
    if (check_for_state()) return;

    if (!strcmp(s, "HAC")) {
        libset_set_bool(SetKey.force_hc, 0);
    } else {
        if (!strncmp(s, "HC", 2)) {
            s += 2;
        }
        if (parse_libset_int_code(SetKey.hc_version, s) == 0) {
            /* non-HAC variant chosen */
            libset_set_bool(SetKey.force_hc, 1);
        }
    }
}

void set_panel_hccme(const(char) * s)
{
    if (check_for_state()) return;

    parse_libset_int_code(SetKey.panel_robust, s);
}

void set_garch_alt_vcv(const(char) * s)
{
    if (check_for_state()) return;

    parse_libset_int_code(SetKey.garch_alt_vcv, s);
}

/* end public functions called from gui/settings.c */

private int set_init_matrix(SetKey key, const(char) * name,
                            Prn * prn)
{
    GretlMatrixRaw ** targ;

    targ = (key == SetKey.initvals) ? &state.initvals : &state.initcurv;

    gretl_matrix_free(*targ);
    *targ = null;

    if (strcmp(name, "auto")) {
        GretlMatrixRaw * m = get_matrix_by_name(name);

        if (m == null) {
            pprintf(prn, _("'%s': no such matrix").ptr, name);
            pputc(prn, '\n');
            return GretlError.data;
        }
        *targ = gretl_matrix_copy(m);
        if (*targ == null) {
            return GretlError.alloc;
        }
    }

    return 0;
}

private int set_matrix_mask(const(char) * name, DatasetRaw * dset)
{
    gretl_matrix_free(state.matmask);
    state.matmask = null;

    if (strcmp(name, "null")) {
        int t, v = current_series_index(dset, name);

        if (v < 0) {
            return GretlError.unkvar;
        }
        state.matmask = gretl_column_vector_alloc(dset.n);
        if (state.matmask == null) {
            return GretlError.alloc;
        }
        for (t=0; t<dset.n; t++) {
            state.matmask.val[t] = dset.Z[v][t];
        }
    }

    return 0;
}

void destroy_matrix_mask()
{
    check_for_state();
    gretl_matrix_free(state.matmask);
    state.matmask = null;
}

private __gshared int function(char * s) workdir_callback;

void set_workdir_callback(int function(char * s) callback)
{
    workdir_callback = callback;
}

private int set_workdir(const(char) * s)
{
    int err = 0;

    if (gretl_function_depth() > 0) {
        gretl_errmsg_set("set workdir: cannot be done inside a function");
        return 1;
    } else if (*s == '\0') {
        return GretlError.data;
    } else {
        char[MAXLEN] workdir;

        workdir[0] = '\0';
        strncat(workdir.ptr, s, MAXLEN - 1);
        if (!err && workdir_callback != null) {
            err = workdir_callback(workdir.ptr);
        } else if (!err) {
            err = gretl_set_path_by_name("workdir", workdir.ptr);
        }
    }

    return err;
}

/* Apparatus for handling tex_plot_opts: we allow setting of this
   variable within a function, but only if it hasn't already been set
   at a "higher" level of execution. If setting within a function
   succeeds, we destroy the setting on exit so that it doesn't
   propagate "upwards".
*/

private enum TPO_UNSET = 1024;
private __gshared char * tex_plot_opts;
private __gshared int tpo_depth = TPO_UNSET;

private int set_tex_plot_opts(const(char) * arg)
{
    int fd = gretl_function_depth();

    if (fd > tpo_depth) {
        gretl_errmsg_sprintf("'%s': cannot be overwritten here",
                             "tex_plot_opts".ptr);
        return GretlError.invarg;
    } else {
        free(tex_plot_opts);
        if (!strcmp(arg, "null") || !strcmp(arg, "none")) {
            tex_plot_opts = null;
        } else {
            tex_plot_opts = gretl_strdup(arg);
            g_strstrip(tex_plot_opts);
        }
        tpo_depth = fd;
    }

    return 0;
}

/* if tex_plot_opts was set within a function, scrub it on exit */

private void maybe_destroy_tpo()
{
    if (gretl_function_depth() == tpo_depth) {
        free(tex_plot_opts);
        tex_plot_opts = null;
        tpo_depth = TPO_UNSET;
    }
}

/* accessor, called from graphing.c */

const(char) * get_tex_plot_opts()
{
    return tex_plot_opts;
}

/* end tex_plot_opts apparatus */

private int set_logfile(const(char) * s)
{
    int err = 0;

    if (gretl_function_depth() > 0) {
        gretl_errmsg_set("set logfile: cannot be done inside a function");
        return GretlError.data;
    } else if (*s == '\0' || !strcmp(s, "null")) {
        gretl_insert_builtin_string("logfile", "");
    } else if (!strcmp(s, "stdout") || !strcmp(s, "stderr")) {
        gretl_insert_builtin_string("logfile", s);
    } else if (strlen(s) >= FILENAME_MAX) {
        gretl_errmsg_set(_("set logfile: filename is too long").ptr);
        return GretlError.invarg;
    } else {
        char[FILENAME_MAX] outname;

        /* switch to workdir if needed */
        strcpy(outname.ptr, s);
        gretl_maybe_prepend_dir(outname.ptr);
        err = gretl_test_fopen(outname.ptr, "w");
        if (!err) {
            gretl_insert_builtin_string("logfile", outname.ptr);
        } else {
            gretl_errmsg_sprintf(_("Couldn't write to %s").ptr, outname.ptr);
            err = GretlError.fopen;
        }
    }

    return err;
}

__gshared const(char) * csv_delims = ", \t;|";

private char delim_from_arg(const(char) * s)
{
    int i;

    for (i=0; csv_strs[i] != null; i++) {
        if (!strcmp(s, csv_strs[i])) {
            return csv_delims[i];
        }
    }

    return 0;
}

private const(char) * arg_from_delim(char c)
{
    int i;

    for (i=0; csv_delims[i] != '\0'; i++) {
        if (c == csv_delims[i]) {
            return csv_strs[i];
        }
    }

    return "unset";
}

private void libset_print_bool(SetKey key, const(char) * s,
                               Prn * prn, Gretlopt opt)
{
    int v = libset_get_bool(key);

    if (gretl_function_depth() > 0 &&
        (key == SetKey.echo_on || key == SetKey.msgs_on || key == SetKey.verbose)) {
        return;
    }

    if (s == null) {
        s = setkey_get_name(key);
    }

    if (opt & Gretlopt.OPT_D) {
        pprintf(prn, " %s = %d\n", s, v);
    } else {
        pprintf(prn, "set %s %s\n", s, (v ? "on" : "off").ptr);
    }
}

private void libset_print_int(SetKey key, const(char) * s,
                              Prn * prn, Gretlopt opt)
{
    if (s == null) {
        s = setkey_get_name(key);
    }

    if (coded_intvar(key)) {
        if (opt & Gretlopt.OPT_D) {
            pprintf(prn, " %s = %s\n", s, libset_option_string(key));
        } else {
            pprintf(prn, "set %s %s\n", s, libset_option_string(key));
        }
    } else {
        int k = libset_get_int(key);

        if (opt & Gretlopt.OPT_D) {
            if (is_unset(k)) {
                pprintf(prn, " %s = auto\n", s);
            } else {
                pprintf(prn, " %s = %d\n", s, k);
            }
        } else if (!is_unset(k)) {
            pprintf(prn, "set %s %d\n", s, k);
        }
    }
}

private void libset_print_double(SetKey key, const(char) * s,
                                 Prn * prn, Gretlopt opt)
{
    double x = libset_get_double(key);

    if (s == null) {
        s = setkey_get_name(key);
    }

    if (opt & Gretlopt.OPT_D) {
        if (na(x)) {
            pprintf(prn, " %s = auto\n", s);
        } else {
            pprintf(prn, " %s = %.15g\n", s, x);
        }
    } else if (!na(x)) {
        pprintf(prn, "set %s %.15g\n", s, x);
    }
}

/* (the C code passes @s through gettext; _() is the identity here) */
private void libset_header(const(char) * s, Prn * prn, Gretlopt opt)
{
    if (opt & Gretlopt.OPT_D) {
        pputs(prn, "\n --- ");
        pputs(prn, s);
        pputs(prn, " ---\n");
    } else {
        pprintf(prn, "# %s\n", s);
    }
}

private void print_vars_for_category(int category, Prn * prn,
                                     Gretlopt opt)
{
    SetVar * v;
    int i;

    for (i=0; i<setvars.length; i++) {
        if (setvars[i].category != category) {
            continue;
        }
        v = &setvars[i];
        if (libset_boolvar(v.key)) {
            libset_print_bool(v.key, v.name, prn, opt);
        } else if (libset_int(setvars[i].key)) {
            libset_print_int(v.key, v.name, prn, opt);
        } else if (libset_double(v.key)) {
            libset_print_double(v.key, v.name, prn, opt);
        }
    }
}

/* print_settings: use OPT_D for "display", otherwise
   this gives script-type output */

private int print_settings(Prn * prn, Gretlopt opt)
{
    const(char) * workdir = gretl_workdir();
    guint64 u;

    gretl_push_c_numeric_locale();

    if (opt & Gretlopt.OPT_D) {
        pputs(prn, _("Variables that can be set using \"set\"").ptr);
        pputs(prn, " (");
        pputs(prn, _("\"help set\" for details").ptr);
        pputs(prn, "):\n");
    }

    libset_header(N_("Program interaction and behavior").ptr, prn, opt);

    if (opt & Gretlopt.OPT_D) {
        pprintf(prn, " workdir = \"%s\"\n", workdir);
    }
    /* (the C source has an `else if (0)` branch here, printing
       `set workdir "..."`: dead code, not ported) */

    if (opt & Gretlopt.OPT_D) {
        pprintf(prn, " csv_delim = %s\n", arg_from_delim(data_delim));
        pprintf(prn, " csv_write_na = %s\n", get_csv_na_write_string());
        pprintf(prn, " csv_read_na = %s\n", get_csv_na_read_string());
        pprintf(prn, " display_digits = %d\n", get_gretl_digits());
        pprintf(prn, " graph_theme = %s\n", get_plotstyle());
        pprintf(prn, " tex_plot_opts = \"%s\"\n",
                tex_plot_opts == null ? "null".ptr : tex_plot_opts);
    } else {
        const(char) * dl = arg_from_delim(data_delim);

        if (strcmp(dl, "unset")) {
            pprintf(prn, "set csv_delim %s\n", arg_from_delim(data_delim));
        }
        pprintf(prn, "set csv_write_na %s\n", get_csv_na_write_string());
        pprintf(prn, "set csv_read_na %s\n", get_csv_na_read_string());
        pprintf(prn, "set graph_theme %s\n", get_plotstyle());
        pprintf(prn, "set tex_plot_opts \"%s\"\n",
                tex_plot_opts == null ? "null".ptr : tex_plot_opts);
    }

    print_vars_for_category(SetCategory.behave, prn, opt);
    if (opt & Gretlopt.OPT_D) {
        /* display only */
        libset_print_bool(SetKey.shell_ok, null, prn, opt);
    }

    libset_header(N_("Numerical methods").ptr, prn, opt);
    print_vars_for_category(SetCategory.numeric, prn, opt);
    if ((opt & Gretlopt.OPT_D) && gretl_function_depth() == 0) {
        /* script version of this? */
        print_state_matrix(SetKey.initvals, prn, opt);
        print_state_matrix(SetKey.initcurv, prn, opt);
    }

    libset_header(N_("Random number generation").ptr, prn, opt);
    u = gretl_rand_get_seed();
    if (opt & Gretlopt.OPT_D) {
        pprintf(prn, " seed = %lu\n", u);
    } else {
        if (seed_is_set) {
            pprintf(prn, "set seed %lu\n", u);
        }
    }

    libset_header(N_("Robust estimation").ptr, prn, opt);
    print_vars_for_category(SetCategory.robust, prn, opt);

    libset_header(N_("Time series").ptr, prn, opt);
    print_vars_for_category(SetCategory.ts, prn, opt);

    pputc(prn, '\n');

    gretl_pop_c_numeric_locale();

    return 0;
}

private int libset_query_settings(SetVar * sv, Prn * prn)
{
    int err = 0;

    if (libset_boolvar(sv.key)) {
        pprintf(prn, "%s: boolean (on/off), currently %s\n",
                sv.name, (libset_get_bool(sv.key) ? "on" : "off").ptr);
    } else if (coded_intvar(sv.key)) {
        pprintf(prn, "%s: code, currently \"%s\"\n", sv.name,
                libset_option_string(sv.key));
        coded_var_show_opts(sv.key, prn);
    } else if (libset_int(sv.key)) {
        int k = libset_get_int(sv.key);

        if (is_unset(k)) {
            pprintf(prn, "%s: positive integer, currently unset\n", sv.name);
        } else if (libset_int_min(sv.key) == 0) {
            pprintf(prn, "%s: non-negative integer, currently %d\n", sv.name, k);
        } else {
            pprintf(prn, "%s: positive integer, currently %d\n", sv.name, k);
        }
    } else if (libset_double(sv.key)) {
        double x = libset_get_double(sv.key);

        if (na(x)) {
            pprintf(prn, "%s: positive floating-point value, "
                    ~ "currently automatic\n", sv.name);
        } else {
            pprintf(prn, "%s: positive floating-point value, "
                    ~ "currently %g\n", sv.name, x);
        }
    } else if (sv.key == SetKey.initvals || sv.key == SetKey.initcurv ||
               sv.key == SetKey.matmask) {
        GretlMatrixRaw * m =
            (sv.key == SetKey.initvals) ? state.initvals :
            (sv.key == SetKey.initcurv) ? state.initcurv : state.matmask;

        if (m != null) {
            pprintf(prn, "%s: matrix, currently\n", sv.name);
            gretl_matrix_print_to_prn(m, null, prn);
        } else {
            pprintf(prn, "%s: matrix, currently null\n", sv.name);
        }
    } else if (sv.key == SetKey.seed) {
        pprintf(prn, "%s: unsigned 64-bit int, "
                ~ "currently %lu (%s)\n",
                sv.name, gretl_rand_get_seed(),
                (seed_is_set ? "set by user" : "automatic").ptr);
    } else if (sv.key == SetKey.csv_delim) {
        pprintf(prn, "%s: named character, currently \"%s\"\n", sv.name,
                arg_from_delim(data_delim));
        coded_var_show_opts(sv.key, prn);
    } else if (sv.key == SetKey.sv_workdir) {
        pprintf(prn, "%s: string, currently \"%s\"\n", sv.name,
                gretl_workdir());
    } else if (sv.key == SetKey.csv_write_na) {
        pprintf(prn, "%s: string, currently \"%s\"\n", sv.name,
                state.csv_write_na.ptr);
    } else if (sv.key == SetKey.csv_read_na) {
        pprintf(prn, "%s: string, currently \"%s\"\n", sv.name,
                state.csv_read_na.ptr);
    } else if (sv.key == SetKey.disp_digits) {
        pprintf(prn, "%s: integer, currently %d\n", sv.name,
                get_gretl_digits());
    } else if (sv.key == SetKey.graph_theme) {
        pprintf(prn, "%s: keyword, currently \"%s\"\n", sv.name,
                get_plotstyle());
    } else if (sv.key == SetKey.tex_plot_opts) {
        pprintf(prn, "%s: string, currently \"%s\"\n", sv.name,
                tex_plot_opts == null ? "null".ptr : tex_plot_opts);
    } else if (sv.key == SetKey.verbose) {
        pprintf(prn, "%s: boolean (on/off), currently %s\n", sv.name,
                ((libset_get_bool(SetKey.echo_on) || libset_get_bool(SetKey.msgs_on)) ?
                "on" : "off").ptr);
        err = 0;
    } else {
        err = 1;
    }

    return err;
}

/* determine whether help text is available pertaining to
   putative libset variable @s.
*/

int libset_help_available(const(char) * s)
{
    SetVar * sv = get_setvar_by_name(s);
    int err = 1;

    if (sv != null) {
        err = libset_query_settings(sv, null);
    }

    return (err == 0);
}

private bool default_ok(SetKey k)
{
    return (k == SetKey.bfgs_toler || k == SetKey.bhhh_toler || k == SetKey.nls_toler);
}

private bool default_str(const(char) * s)
{
    return (!strcmp(s, "auto") || !strcmp(s, "default"));
}

private bool boolean_on(const(char) * s)
{
    return (!strcmp(s, "on") || !strcmp(s, "1") ||
            !strcmp(s, "true") || !strcmp(s, "TRUE"));
}

private bool boolean_off(const(char) * s)
{
    return (!strcmp(s, "off") || !strcmp(s, "0") ||
            !strcmp(s, "false") || !strcmp(s, "FALSE"));
}

private int write_or_read_settings(Gretlopt opt, Prn * prn)
{
    int err = incompatible_options(opt, (Gretlopt.OPT_T | Gretlopt.OPT_F));

    if (!err) {
        const(char) * fname = get_optval_string(GretlCmdIndex.set, opt);

        if (fname == null) {
            err = GretlError.data;
        } else if (opt == Gretlopt.OPT_T) {
            err = libset_write_script(fname);
        } else {
            err = real_libset_read_script(fname, prn);
        }
    }

    return err;
}

private int check_set_bool(SetKey key, const(char) * name,
                           const(char) * arg)
{
    if (boolean_on(arg)) {
        return libset_set_bool(key, 1);
    } else if (boolean_off(arg)) {
        return libset_set_bool(key, 0);
    } else {
        gretl_errmsg_sprintf(_("%s: invalid value '%s'").ptr, name, arg);
        return GretlError.parse;
    }
}

private int legacy_set_pcse(const(char) * arg)
{
    int err = 0;

    if (arg == null || *arg == '\0') {
        err = GretlError.parse;
    } else if (boolean_on(arg)) {
        libset_set_int(SetKey.panel_robust, PanelRobust.beck_katz);
    } else if (boolean_off(arg)) {
        libset_set_int(SetKey.panel_robust, PanelRobust.arellano);
    } else {
        gretl_errmsg_sprintf(_("%s: invalid value '%s'").ptr, "pcse".ptr, arg);
        err = GretlError.parse;
    }

    return err;
}

private int set_display_digits(const(char) * arg)
{
    if (gretl_function_depth() > 0) {
        gretl_errmsg_sprintf("'%s': cannot be set inside a function",
                             "display_digits".ptr);
        return GretlError.invarg;
    } else {
        return set_gretl_digits(atoi(arg));
    }
}

private int set_verbosity(const(char) * arg)
{
    int err = 0;

    if (!strcmp(arg, "on")) {
        set_gretl_messages(1);
        set_gretl_echo(1);
    } else if (!strcmp(arg, "off")) {
        set_gretl_messages(0);
        set_gretl_echo(0);
        comments_on = 0;
    } else if (!strcmp(arg, "comments")) {
        set_gretl_messages(0);
        set_gretl_echo(0);
        comments_on = 1;
    } else {
        err = GretlError.invarg;
    }

    return err;
}

int execute_set(const(char) * setobj, const(char) * setarg,
                DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    SetVar * sv = null;
    int k, argc, err;
    guint64 u;

    check_for_state();

    if (opt != Gretlopt.OPT_NONE) {
        return write_or_read_settings(opt, prn);
    }

    argc = (setobj != null) + (setarg != null);
    if (argc == 0) {
        return print_settings(prn, Gretlopt.OPT_D);
    }

    if (!strcmp(setobj, "pcse")) {
        return legacy_set_pcse(setarg);
    } else if (!strcmp(setobj, "dcmt")) {
        pputs(prn, "'dcmt': obsolete setting\n");
        return 0;
    } else if (!strcmp(setobj, "omp_mnk_min") ||
               !strcmp(setobj, "mp_mnk_min")) {
        pputs(prn, "'omp_mnk_min': obsolete setting, ignored\n");
        return 0;
    } else {
        sv = get_setvar_by_name(setobj);
    }

    if (sv == null) {
        gretl_errmsg_sprintf(_("set: unknown variable '%s'").ptr, setobj);
        return GretlError.data;
    }

    /* set error default */
    err = GretlError.parse;

    if (argc == 1) {
        if (sv.key == SetKey.stopwatch) {
            gretl_stopwatch();
            return 0;
        } else {
            return libset_query_settings(sv, prn);
        }
    } else if (argc == 2) {
        /* specials first */
        if (sv.key == SetKey.csv_write_na) {
            return set_csv_na_write_string(setarg);
        } else if (sv.key == SetKey.csv_read_na) {
            return set_csv_na_read_string(setarg);
        } else if (sv.key == SetKey.initvals || sv.key == SetKey.initcurv) {
            return set_init_matrix(sv.key, setarg, prn);
        } else if (sv.key == SetKey.matmask) {
            return set_matrix_mask(setarg, dset);
        } else if (sv.key == SetKey.sv_workdir) {
            return set_workdir(setarg);
        } else if (sv.key == SetKey.sv_logfile) {
            return set_logfile(setarg);
        } else if (sv.key == SetKey.graph_theme) {
            return set_plotstyle(setarg);
        } else if (sv.key == SetKey.disp_digits) {
            return set_display_digits(setarg);
        } else if (sv.key == SetKey.verbose) {
            return set_verbosity(setarg);
        } else if (sv.key == SetKey.tex_plot_opts) {
            return set_tex_plot_opts(setarg);
        }

        if (libset_boolvar(sv.key)) {
            if (sv.key == SetKey.shell_ok) {
                pprintf(prn, "'%s': this must be set via the gretl GUI\n", setobj);
                err = GretlError.data;
            } else {
                err = check_set_bool(sv.key, setobj, setarg);
            }
        } else if (libset_double(sv.key)) {
            if (default_ok(sv.key) && default_str(setarg)) {
                libset_set_double(sv.key, NADBL);
                err = 0;
            } else {
                double x;

                err = libset_get_scalar(sv.key, setarg, null, &x);
                if (!err) {
                    err = libset_set_double(sv.key, x);
                }
            }
        } else if (sv.key == SetKey.csv_delim) {
            char c = delim_from_arg(setarg);

            if (c > 0) {
                data_delim = c;
                err = 0;
            }
        } else if (sv.key == SetKey.seed) {
            int automatic = 0;

            err = libset_get_seed(setarg, &u, &automatic);
            if (!err) {
                gretl_rand_set_seed(u);
                if (gretl_messages_on()) {
                    pprintf(prn,
                            _("Pseudo-random number generator seeded with %u\n").ptr, u);
                }
                seed_is_set = !automatic;
            }
        } else if (sv.key == SetKey.horizon) {
            /* horizon for VAR impulse responses */
            if (!strcmp(setarg, "auto")) {
                state.horizon = UNSET_INT;
                err = 0;
            } else {
                err = libset_get_scalar(sv.key, setarg, &k, null);
                if (!err) {
                    state.horizon = k;
                } else {
                    state.horizon = UNSET_INT;
                }
            }
        } else if (coded_intvar(sv.key)) {
            err = parse_libset_int_code(sv.key, setarg);
        } else if (libset_int(sv.key)) {
            err = libset_get_scalar(sv.key, setarg, &k, null);
            if (!err) {
                err = libset_set_int(sv.key, k);
            }
        } else {
            gretl_errmsg_sprintf(_("set: unknown variable '%s'").ptr, setobj);
            err = GretlError.unkvar;
        }
    }

    return err;
}

double libset_get_double(SetKey key)
{
    void * valp;

    if (check_for_state()) {
        return NADBL;
    }

    valp = setkey_get_target(key, SVType.double_);
    if (valp != null) {
        double x = *cast(double *) valp;

        if (na(x) && (key == SetKey.nls_toler || key == SetKey.bfgs_toler)) {
            x = get_default_nls_toler();
        }
        return x;
    } else {
        fprintf(stderr, "libset_get_double: unrecognized "
                ~ "key %d\n", key);
        return 0;
    }
}

double libset_get_user_tolerance(SetKey key)
{
    if (key >= SetKey.nls_toler && key <= SetKey.bhhh_toler) {
        void * valp = setkey_get_target(key, SVType.all);

        return *cast(double *) valp;
    } else {
        return NADBL;
    }
}

int libset_set_double(SetKey key, double val)
{
    void * valp;
    int err = 0;

    if (check_for_state()) {
        return 1;
    }

    /* all the libset double vals must be positive */
    if (val <= 0.0 || na(val)) {
        return GretlError.data;
    }

    valp = setkey_get_target(key, SVType.double_);
    if (valp != null) {
        *cast(double *) valp = val;
    } else {
        fprintf(stderr, "libset_set_double: unrecognized key %d (%s)\n",
                key, setkey_get_name(key));
        err = GretlError.unkvar;
    }

    return err;
}

int libset_get_int(SetKey key)
{
    void * valp;

    if (check_for_state()) {
        return 0;
    }

    valp = setkey_get_target(key, SVType.int_);

    if (valp != null) {
        if (libset_small_int(key)) {
            return *cast(gint8 *) valp;
        } else {
            return *cast(int *) valp;
        }
    } else if (key == SetKey.blas_mnk_min) {
        return get_blas_mnk_min();
    } else if (key == SetKey.omp_n_threads) {
        return gretl_get_omp_threads();
    } else if (key == SetKey.simd_k_max) {
        return get_simd_k_max();
    } else if (key == SetKey.simd_mn_min) {
        return get_simd_mn_min();
    } else {
        fprintf(stderr, "libset_get_int: unrecognized "
                ~ "key %d\n", key);
        return 0;
    }
}

private struct IntLimits {
    SetKey key;
    int min;
    int max;
}

private int get_int_limits(SetKey key, int * min, int * max)
{
    static immutable IntLimits[15] ilims = [
        IntLimits(SetKey.hc_version, 0, 4),
        IntLimits(SetKey.panel_robust, 0, 2),
        IntLimits(SetKey.fdjac_qual, 0, 2),
        IntLimits(SetKey.use_qr, 0, 2),
        IntLimits(SetKey.lbfgs_mem,  3, 20),
        IntLimits(SetKey.gretl_debug, 0, 4),
        IntLimits(SetKey.datacols,    1, 15),
        IntLimits(SetKey.plot_collect, 0, 2),
        IntLimits(SetKey.csv_digits, 1, 25),
        IntLimits(SetKey.boot_iters, 499, 999999),
        IntLimits(SetKey.bfgs_verbskip, 0, 1000),
        IntLimits(SetKey.bootrep, 1, 99999),
        IntLimits(SetKey.horizon, 1, 1000),
        IntLimits(SetKey.loop_maxiter, 0, int.max - 1),
        IntLimits(SetKey.gmp_bits, 256, 8192)
    ];
    int i;

    for (i=0; i<ilims.length; i++) {
        if (ilims[i].key == key) {
            *min = ilims[i].min;
            *max = ilims[i].max;
            return 1;
        }
    }

    return 0;
}

private int libset_int_min(SetKey key)
{
    int m1, m0 = 1;

    get_int_limits(key, &m0, &m1);
    return m0;
}

/* Called from within libset.c and also from various places in
   libgretl. It's primarily designed for "real" integer variables
   (not int-coded categories), but for now we make an exception
   for HC_VERSION and PLOT_COLLECT, to support existing calls
   from lib/src/estimate.c and gui/settings.c.
*/

int libset_set_int(SetKey key, int val)
{
    int err = 0;

    if (check_for_state()) {
        return 1;
    }

    if (key == SetKey.blas_mnk_min) {
        set_blas_mnk_min(val);
    } else if (key == SetKey.simd_k_max) {
        set_simd_k_max(val);
    } else if (key == SetKey.simd_mn_min) {
        set_simd_mn_min(val);
    } else if (key == SetKey.omp_n_threads) {
        err = gretl_set_omp_threads(val);
    } else {
        int min = 1, max = 100000;
        void * valp;

        get_int_limits(key, &min, &max);
        if (val < min || val > max) {
            err = GretlError.data;
        } else {
            valp = setkey_get_target(key, SVType.int_);
            if (valp == null) {
                err = GretlError.data;
            } else if (libset_small_int(key)) {
                *cast(gint8 *) valp = cast(gint8) val;
            } else {
                *cast(int *) valp = val;
            }
        }
    }

    return err;
}

private void set_flag_from_env(SetKey flag, const(char) * s)
{
    char * e = getenv(s);

    if (e != null) {
        if (*e == '\0' || *e == '0') {
            state.flags &= ~flag;
        } else {
            state.flags |= flag;
        }
    }
}

private void set_int_from_env(SetKey key, const(char) * s)
{
    char * e = getenv(s);

    if (e != null) {
        if (*e == '\0' || *e == '0') {
            libset_set_int(key, 0);
        } else {
            libset_set_int(key, atoi(e));
        }
    }
}

private void maybe_check_env(SetKey key)
{
    if (key == SetKey.use_svd) {
        set_flag_from_env(SetKey.use_svd, "GRETL_USE_SVD");
    } else if (key == SetKey.use_lbfgs) {
        set_flag_from_env(SetKey.use_lbfgs, "GRETL_USE_LBFGS");
    } else if (key == SetKey.use_qr) {
        set_int_from_env(SetKey.use_qr, "GRETL_USE_QR");
    }
}

int libset_get_bool(SetKey key)
{
    /* global specials */
    if (key == SetKey.r_functions) {
        return globals.R_functions;
    } else if (key == SetKey.r_lib) {
        return globals.R_lib;
    } else if (key == SetKey.logstamp) {
        return globals.logstamp;
    }

    if (check_for_state()) {
        return 0;
    } else {
        maybe_check_env(key);
        return flag_to_bool(state, key);
    }
}

/* (In the C source the body of this function is under #ifdef ENABLE_NLS;
   NLS is not enabled in this port, so it is a no-op, as in a non-NLS build.) */
private void libset_set_decpoint(int on)
{
}

void set_data_export_decimal_comma(int s)
{
    if (s) {
        export_decpoint = ',';
    } else {
        export_decpoint = '.';
    }
}

char get_data_export_decpoint()
{
    char c = export_decpoint;

    /* revert to '.' on access */
    export_decpoint = '.';
    return c;
}

void set_data_export_delimiter(char c)
{
    data_delim = c;
}

char get_data_export_delimiter()
{
    return data_delim;
}

/* (USE_RLIB is not defined in this port: only the "not supported in this
   build" variant of this function is ported.) */
private int check_R_setting(gint8 * var, SetKey key, int val)
{
    int err = 0;

    if (val) {
        const(char) * s = (key == SetKey.r_functions) ? "R_functions" : "R_lib";

        gretl_errmsg_sprintf(_("%s: not supported in this build of gretl").ptr, s);
        err = GretlError.external;
    }

    return err;
}

int libset_set_bool(SetKey key, int val)
{
    if (check_for_state()) {
        return GretlError.alloc;
    }

    /* global specials */
    if (key == SetKey.r_functions) {
        return check_R_setting(&globals.R_functions, key, val);
    } else if (key == SetKey.r_lib) {
        return check_R_setting(&globals.R_lib, key, val);
    } else if (key == SetKey.logstamp) {
        globals.logstamp = cast(gint8) val;
        return 0;
    }

    if (val) {
        state.flags |= key;
    } else {
        state.flags &= ~key;
    }

    if (key == SetKey.force_decpoint) {
        libset_set_decpoint(val);
    }

    return 0;
}

/* Mechanism for pushing and popping program state for user-defined
   functions. push_program_state() is used when a function starts
   execution: the function gets a copy of the current program state,
   while that state is pushed onto the stack for restoration when the
   function exits.
*/

private __gshared int n_states;
private __gshared GPtrArray * state_stack;
private __gshared int state_idx = -1;

int push_program_state()
{
    SetState * newstate;
    int err = 0;

    if (n_states == 0) {
        state_stack = g_ptr_array_new();
    }

    state_idx++;

    if (state_idx < n_states) {
        newstate = cast(SetState *) state_stack.pdata[state_idx];
    } else {
        newstate = cast(SetState *) malloc((*newstate).sizeof);
        if (newstate == null) {
            err = GretlError.alloc;
        } else {
            g_ptr_array_add(state_stack, newstate);
            n_states++;
        }
    }

    if (newstate != null) {
        if (n_states == 1) {
            state_vars_init(newstate);
        } else {
            state_vars_copy(newstate);
        }
        state = newstate;
    }

    return err;
}

private void free_state(SetState * sv)
{
    if (sv != null) {
        gretl_matrix_free(sv.initvals);
        gretl_matrix_free(sv.initcurv);
        gretl_matrix_free(sv.matmask);
        free(sv);
    }
}

/* Called when a user-defined function exits: restores the program
   state that was in force when the function started executing.
*/

int pop_program_state()
{
    int err = 0;

    if (n_states < 2) {
        err = 1;
    } else {
        int fdp0, fdp = state.flags & SetKey.force_decpoint;

        /* restore prior stack level */
        state = cast(SetState *) state_stack.pdata[--state_idx];

        fdp0 = state.flags & SetKey.force_decpoint;
        if (fdp0 != fdp) {
            libset_set_decpoint(fdp0);
        }
        if (tex_plot_opts != null) {
            maybe_destroy_tpo();
        }
    }

    return err;
}

/* initialization of all user-settable settings */

int libset_init()
{
    static __gshared int done;
    int err = 0;

    if (!done) {
        err = push_program_state();
        done = 1;
    }

    if (getenv("CLI_NO_PLOTS")) {
        suppress_plots = 1;
    }

    return err;
}

/* state variables for looping */
private __gshared char * looping;
private __gshared int looplen;

void libset_cleanup()
{
    int i;

    for (i=0; i<n_states; i++) {
        free_state(cast(SetState *) state_stack.pdata[i]);
    }

    g_ptr_array_free(state_stack, 1);
    state_stack = null;
    n_states = 0;
    state_idx = -1;

    libset_hash_cleanup();

    free(tex_plot_opts);
    tex_plot_opts = null;
    tpo_depth = TPO_UNSET;

    free(looping);
    looping = null;
    looplen = 0;
}

/* switches for looping and batch mode: output: these depend on the
   state of the program calling libgretl, they are not user-settable
*/

private enum LMIN = 8;

void set_loop_on()
{
    int fd = gretl_function_depth();

    if (looping == null) {
        looplen = fd+1 < LMIN ? LMIN : fd+1;
        looping = cast(char *) calloc(1, looplen);
    } else if (looplen < fd+1) {
        int n = looplen * 2;

        n = n < fd+1 ? fd+1 : n;
        looping = cast(char *) realloc(looping, n);
        memset(looping + looplen, 0, n - looplen);
        looplen = n;
    }
    looping[fd] = 1;
}

void set_loop_off()
{
    if (looping != null) {
        looping[gretl_function_depth()] = 0;
    }
}

/* returns 1 if there's a loop going on anywhere in the "caller
   ancestry" of the current execution level, else 0.
*/

int gretl_looping()
{
    int i, fd = gretl_function_depth();
    int n = min(fd+1, looplen);

    for (i=0; i<n; i++) {
        if (looping[i]) {
            return 1;
        }
    }
    return 0;
}

/* returns 1 if there's a loop going on at the current execution
   stack level, else 0.
*/

int gretl_looping_currently()
{
    int fd = gretl_function_depth();

    return looplen >= fd+1 && looping[fd];
}

private __gshared int iter_depth;

void gretl_iteration_push()
{
    iter_depth++;
}

void gretl_iteration_pop()
{
    if (iter_depth > 0) {
        iter_depth--;
    }
}

int gretl_iteration_depth()
{
    return iter_depth;
}

/* portmanteau: a "loop" is active or (internal)
   iteration is in progress (e.g. in BFGS)
*/

int gretl_iterating()
{
    return gretl_looping() || iter_depth > 0;
}

private int batch_mode_switch(int set, int val)
{
    static __gshared int bmode;

    if (set) {
        bmode = val;
    }

    return bmode;
}

void gretl_set_batch_mode(int b)
{
    batch_mode_switch(1, b);
}

/* Returns 1 if we're running a script, otherwise 0.
   Note: a 0 return indicates that we're in an interactive
   mode, whether GUI or CLI.
*/

int gretl_in_batch_mode()
{
    return batch_mode_switch(0, 0);
}

private __gshared int gui_mode;

/* set by the GUI program at start-up */

void gretl_set_gui_mode()
{
    gui_mode = 1;
}

/* Returns 1 if we're running the GUI program. The current
   usage may be interactive (menu-driven or typing at the
   GUI "console") or script/batch. See also
   gretl_in_batch_mode().
*/

int gretl_in_gui_mode()
{
    return gui_mode;
}

/* "tool_mode" is set when gretlcli is being used as
   a build tool */

private __gshared int tool_mode;

void gretl_set_tool_mode()
{
    set_gretl_echo(0);
    set_gretl_messages(0);
    tool_mode = 1;
}

int gretl_in_tool_mode()
{
    return tool_mode;
}

/* "python_mode" is set when libgretl is being called from
   python, via gretl4py.
*/

private __gshared int python_mode;

void gretl_set_python_mode()
{
    python_mode = 1;
}

int gretl_in_python_mode()
{
    return python_mode;
}

/* for CLI use: suppress production of plots? */

void gretl_set_no_plots()
{
    suppress_plots = 1;
}

int gretl_no_plots()
{
    return gui_mode ? 0 : suppress_plots;
}

/* mechanism to support callback for representing ongoing
   activity in the GUI */

private __gshared ShowActivityFunc sfunc;

void set_show_activity_func(ShowActivityFunc func)
{
    sfunc = func;
}

int show_activity_func_installed()
{
    return sfunc != null;
}

void show_activity_callback()
{
    if (sfunc != null) {
        sfunc();
    }
}

/* Support for the GUI "Stop" button, to terminate execution
   of a long-running script.
*/

private __gshared int user_stop;

void set_user_stop(int s)
{
    user_stop = s;
}

int get_user_stop()
{
    return user_stop;
}

/* internal: set the value of a string-valued "setvar" */

private int set_string_setvar(char * targ, const(char) * s, int len)
{
    *targ = '\0';

    if (*s == '"') {
        /* @s is quoted */
        const(char) * p = strchr(s+1, '"');
        int n;

        if (p == null) {
            return GretlError.parse;
        } else {
            n = cast(int) (p - s - 1);
            if (n > len) {
                n = len;
            }
            strncat(targ, s+1, n);
        }
    } else {
        strncat(targ, s, len);
    }

    return 0;
}

/* for setting what we print for NAs on CSV output */

const(char) * get_csv_na_write_string()
{
    if (check_for_state()) {
        return "NA";
    } else {
        return state.csv_write_na.ptr;
    }
}

int set_csv_na_write_string(const(char) * s)
{
    if (check_for_state()) {
        return GretlError.data;
    } else {
        return set_string_setvar(state.csv_write_na.ptr, s, 7);
    }
}

/* and for setting what we read as NA on CSV input */

const(char) * get_csv_na_read_string()
{
    if (check_for_state()) {
        return "default";
    } else {
        return state.csv_read_na.ptr;
    }
}

int set_csv_na_read_string(const(char) * s)
{
    if (check_for_state()) {
        return GretlError.data;
    } else {
        return set_string_setvar(state.csv_read_na.ptr, s, 7);
    }
}

int libset_write_script(const(char) * fname)
{
    Prn * prn;
    int err = 0;

    /* FIXME maybe adjust path for fname? */

    prn = gretl_print_new_with_filename(fname, &err);

    if (!err) {
        gretl_push_c_numeric_locale();
        print_settings(prn, Gretlopt.OPT_NONE);
        gretl_print_destroy(prn);
        gretl_pop_c_numeric_locale();
    }

    return err;
}

private char * get_quoted_arg(const(char) * s, int * err)
{
    const(char) * p;
    int n = 0, matched = 0;
    char * ret = null;

    p = s = strchr(s, '"') + 1;

    while (*p) {
        if (*p == '"') {
            matched = 1;
            break;
        } else {
            n++;
        }
        p++;
    }

    if (!matched) {
        *err = GretlError.parse;
    } else {
        ret = gretl_strndup(s, n);
        if (ret == null) {
            *err = GretlError.alloc;
        }
    }

    return ret;
}

/* If @prn is non-null we're called via the "set" command. Otherwise
   we're called by the gretl GUI session apparatus.
*/

private int real_libset_read_script(const(char) * fname,
                                    Prn * prn)
{
    FILE * fp;
    int err = 0;

    fp = gretl_fopen(fname, "r");

    if (fp == null) {
        char[FILENAME_MAX] fullname;

        strcpy(fullname.ptr, fname);
        gretl_addpath(fullname.ptr, 0);
        fp = gretl_fopen(fname, "r");
        if (fp == null) {
            err = GretlError.fopen;
        }
    }

    if (!err) {
        char[32] setobj, setarg;
        char[1024] line;
        int nf;

        while (fgets(line.ptr, cast(int) line.sizeof, fp)) {
            if (line[0] == '#' || string_is_blank(line.ptr)) {
                continue;
            }
            tailstrip(line.ptr);
            nf = sscanf(line.ptr, "%*s %31s %31s", setobj.ptr, setarg.ptr);
            if (nf == 1) {
                err = execute_set(setobj.ptr, null, null, Gretlopt.OPT_NONE, prn);
            } else if (nf == 2) {
                if (setarg[0] == '"') {
                    char * q = get_quoted_arg(line.ptr, &err);

                    if (!err) {
                        err = execute_set(setobj.ptr, q, null, Gretlopt.OPT_NONE, prn);
                        free(q);
                    }
                } else {
                    err = execute_set(setobj.ptr, setarg.ptr, null, Gretlopt.OPT_NONE, prn);
                }
            }
            if (err && prn != null) {
                break;
            }
        }

        fclose(fp);
    }

    return err;
}

int libset_read_script(const(char) * fname)
{
    return real_libset_read_script(fname, null);
}
