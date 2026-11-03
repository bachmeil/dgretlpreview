/* Port of genparse.h: shared private header for all 'genr' related modules (the types and the flag
   enums; the prototypes of functions defined in genlex.c, geneval.c and genfuncs.c are not ported here: they
   are declared extern (C) where used, until those files are ported).

   Naming (decision for confirmation): the two anonymous enums of operators/types/punctuation and of functions
   share no common prefix, so they are the named enums GenSym and GenFunc whose members are the C names in
   lower case (U_NEG -> u_neg, NUM -> num, F_ABS -> f_abs). Node.t and Parser.sym hold values of both enums
   (the C keeps their numeric ranges disjoint), so they are plain ints/shorts.
   CONST_* and DUM_* become ConstId and DumId (prefix dropped; D keywords get a trailing underscore),
   enum node_flags becomes NodeFlags (the _NODE suffix dropped) and genflags becomes GenFlags (P_ dropped).
   The function-like macros are functions. */

module gretl.genparseh;


import gretl.glib : gint16, guint8;
import gretl.gretl_matrixh : GretlMatrixRaw;
import gretl.gretl_prn : Prn;
import gretl.gretl_string_table : SeriesTable;
import gretl.libgretlh;
import gretl.usermath : MatrixSubspec;
import gretl.uservarh : UserVar;

enum GENDEBUG = 0;

/* operators, types, punctuation */

enum GenSym {
    u_neg = 1,
    u_pos,
    u_not,
    u_addr,
    u_max,  /* SEPARATOR: end of unary operators */
    b_asn,
    b_add,
    b_sub,
    b_mul,
    b_div,
    b_mod,
    b_pow,
    b_eq,
    b_lt,
    b_gt,
    b_lte,
    b_gte,
    b_neq,
    b_and,
    b_or,
    b_trmul,
    b_range,
    b_dotmult,
    b_dotdiv,
    b_dotpow,
    b_dotadd,
    b_dotsub,
    b_doteq,
    b_dotlt,
    b_dotgt,
    b_dotlte,
    b_dotgte,
    b_dotneq,
    b_dotasn,
    b_kron,  /* Kronecker product */
    b_hcat,  /* horizontal concatenation */
    b_vcat,  /* vertical concatenation */
    b_lcat,  /* list concatentation */
    b_ldiv,  /* matrix left division */
    b_ellip,  /* list-generating ellipsis */
    b_join,  /* list-joining with separator */
    op_max,  /* SEPARATOR: end of binary operators */
    g_lpr,  /* grouping: left paren */
    g_rpr,  /* right paren */
    g_lbr,  /* left bracket */
    g_rbr,  /* right bracket */
    g_lcb,  /* left curly bracket */
    g_rcb,  /* right curly bracket */
    p_com,  /* punctuation: comma */
    p_dot,  /* period */
    p_semi,  /* semi-colon */
    p_col,  /* colon */
    punct_max,  /* SEPARATOR: end of grouping and punctuation marks */
    num,  /* scalar */
    series,  /* series */
    list,  /* list of series */
    mat,  /* matrix */
    bundle,  /* gretl bundle (hash table) */
    array,  /* generic array object */
    str,  /* string */
    cnum,  /* constant (literal) numeric value */
    cstr,  /* constant (literal) string */
    con,  /* named numeric constant */
    dum,  /* dummy variable */
    uobj,  /* user-defined object (e.g. model) */
    num_p,  /* user scalar++ */
    num_m,  /* user scalar-- */
    obs,  /* observation from a series */
    dmstr,  /* "dollar" matrix plus string subspec */
    slraw,  /* unevaluated "slice" specification */
    mspec,  /* evaluated matrix subspec */
    subsl,  /* row or column component of MSPEC */
    lag,  /* variable plus lag length */
    dvar,  /* $ "dataset" variable (mostly scalar or series) */
    mvar,  /* $ model var (scalar, series, matrix, or array) */
    listvar,  /* variable in list, dot syntax */
    dbundle,  /* $ bundle accessor */
    bmemb,  /* member of bundle */
    dbmemb,  /* member of $ bundle */
    mmemb,  /* member of named model */
    fargs,  /* set of n function arguments */
    wlist,  /* wildcard list spec */
    empty,  /* "null" or empty arg slot */
    undef,  /* undefined (allowed in "query" context only) */
    nullarg,  /* record of missing function argument */
    dtype_max,  /* SEPARATOR: end of "bare" types */
    ufun,  /* user-defined function */
    rfun,  /* GNU R function */
    ivec,  /* array of ints, not a varlist */
    osl,  /* "slice" of object (matrix, array, string) */
    useries,  /* named series (defined only for error reporting) */
    sub_addr,  /* "address" of (e.g.) array element */
    inc,  /* increment */
    dec,  /* decrement */
    query,  /* ternary "?" expression */
    eot,  /* end of transmission */
    unk}

/* functions: don't collide with the enumeration above */

enum GenFunc {
    f1_min = 1 << 8,
    f_abs,
    f_sgn,
    f_ceil,
    f_floor,
    f_sin,
    f_cos,
    f_tan,
    f_asin,
    f_acos,
    f_atan,
    f_sinh,
    f_cosh,
    f_tanh,
    f_asinh,
    f_acosh,
    f_atanh,
    f_log,
    f_log10,
    f_log2,
    f_exp,
    f_sqrt,
    f_gamma,
    f_lngamma,
    f_digamma,
    f_trigamma,
    f_invmills,
    f_round,
    f_cnorm,
    f_dnorm,
    f_qnorm,
    f_logistic,
    f_real,
    f_imag,
    f_carg,
    f_cmod,
    f_cquad,
    fp_max,  /* separator: end of pointerized functions */
    f_conj,
    f_toint,
    f_diff,  /* first difference */
    f_sdiff,  /* seasonal difference */
    f_sort,  /* ascending sort */
    f_dsort,  /* descending sort */
    f_ranking,
    f_odev,  /* orthogonal deviation */
    f_nobs,
    f_cum,
    f_missing,
    f_dataok,
    f_misszero,
    f_zeromiss,
    f_median,
    f_gini,
    f_sumall,
    f_skewness,
    f_kurtosis,
    f_sst,
    f_chol,
    f_inv,
    f_diag,
    f_transp,
    f_vec,
    f_rows,
    f_cols,
    f_det,
    f_ldet,
    f_trace,
    f_norm1,
    f_infnorm,
    f_rcond,
    f_obsnum,
    f_isdiscr,
    f_isdummy,
    f_typeof,
    f_typename,
    f_exists,
    f_nelem,
    f_values,
    f_uniq,
    f_mexp,
    f_fft,
    f_ffti,
    f_upper,
    f_lower,
    f_obslabel,
    f_backtick,
    f_strlen,
    f_varname,
    f_varnames,
    f_varnum,
    f_tolower,
    f_toupper,
    f_irr,
    f_errmsg,
    f_getenv,
    f_ngetenv,
    f_isnan,
    f_typestr,
    f_strstrip,
    f_remove,
    f_atof,
    f_mpi_recv,
    f_easter,
    f_curl,
    f_nlines,
    f_array,
    f_tramolin,
    f_cnumber,
    f_ecdf,
    f_sleep,
    f_getinfo,
    f_cdumify,
    f_getkeys,
    f_mcorr,
    f_iscmplx,
    f_ctrans,
    f_mlog,
    f_barrier,
    f_lpsolve,
    f_interpol,
    f_dec2bin,
    f_bin2dec,
    f_access,
    f_polroots,
    f_randstr,
    f_which,
    f_gpcolors,
    f_permute,
    hf_jbterms,
    hf_fdepth,
    f1_max,  /* SEPARATOR: end of single-arg functions */
    hf_listinfo,
    f_min,
    f_max,
    f_rank,
    f_ginv,
    f_sum,
    f_mean,
    f_vce,
    f_sd,
    f_argname,
    f_t1,
    f_t2,
    f_cov,
    f_mcov,
    f_dumify,
    f_sortby,
    f_runiform,
    f_rnormal,
    f_fracdiff,
    f_boxcox,
    f_zeros,
    f_ones,
    f_munif,
    f_qform,
    f_eigsym,
    f_quantile,
    f_hdprod,  /* horizontal direct product */
    f_mxtab,
    f_mrsel,
    f_mcsel,
    f_cnameset,
    f_rnameset,
    f_ljungbox,
    f_msortby,
    f_lincomb,
    f_imhof,
    f_fraclag,
    f_mrev,
    f_pergm,
    f_npv,
    f_dsum,
    f_polyfit,
    f_inlist,
    f_isconst,
    f_inbundle,
    f_cnameget,
    f_rnameget,
    f_pnobs,
    f_pmin,
    f_pmax,
    f_psum,
    f_pmean,
    f_pxsum,
    f_pxnobs,
    f_pxmean,
    f_psd,
    f_pshrink,
    f_randint,
    f_mread,
    f_bread,
    f_getline,
    f_isodate,
    f_juldate,
    f_readfile,
    f_printf,
    f_sprintf,
    f_mpi_send,
    f_bcast,
    f_allreduce,
    f_genseries,
    f_kpsscrit,
    f_stringify,
    f_square,
    f_seasonals,
    f_dropcoll,
    f_ksimdata,
    f_hfdiff,
    f_hfldiff,
    f_naalen,
    f_kmeier,
    f_normtest,
    f_cor,
    f_lrcovar,
    f_fixname,
    f_atan2,
    f_ccode,
    f_lsolve,
    f_strfday,
    f_strptime,
    f_strpday,
    f_conv2d,
    f_flatten,
    f_imat,
    f_complex,
    f_randperm,
    f_cswitch,
    f_psdroot,
    f_strvals,
    f_errorif,
    f_bincoeff,
    f_assert,
    f_contains,
    f_vech,
    f_unvech,
    f_deseas,
    f_mat2list,
    f_cmult,
    f_cdiv,
    f_pexpand,
    f_invpd,
    f_minc,
    f_maxc,
    f_minr,
    f_maxr,
    f_iminc,
    f_imaxc,
    f_iminr,
    f_imaxr,
    f_sumc,
    f_sumr,
    f_prodc,
    f_prodr,
    f_meanc,
    f_meanr,
    f_asort,
    f_corresp,
    f_strvsort,
    f_fevalb,
    f_binperms,
    f_mnorm,
    f_lnmgamma,
    f_nullspc,
    hf_vcnorm,
    hf_glasso,
    f2_max,  /* SEPARATOR: end of two-arg functions */
    f_wmean,
    f_wvar,
    f_wsd,
    f_llag,
    f_hflag,
    f_princomp,
    f_bfgsmax,
    f_mshape,
    f_svd,
    f_qr,
    f_trimr,
    f_corrgm,
    f_seq,
    f_replace,
    f_strncmp,
    f_bessel,
    f_weekday,
    f_monthlen,
    f_epochday,
    f_kdensity,
    f_setnote,
    f_bwfilt,
    f_varsimul,
    f_strsub,
    f_regsub,
    f_mlag,
    f_eigsolve,
    f_simann,
    f_halton,
    f_mwrite,
    f_bwrite,
    f_aggrby,
    f_iwishart,
    f_sscanf,
    f_substr,
    f_reduce,
    f_scatter,
    f_mweights,
    f_mgradient,
    f_mlincomb,
    f_hflist,
    f_nmmax,
    f_gssmax,
    f_npcorr,
    f_dayspan,
    f_smplspan,
    f_fdjac,
    f_numhess,
    f_strsplit,
    f_hpfilt,
    f_xmlget,
    f_jsonget,
    f_fevd,
    f_lrvar,
    f_brename,
    f_isoweek,
    f_bkw,
    f_fzero,
    f_eigen,
    f_schur,
    f_resample,
    f_stack,
    f_geoplot,
    f_vma,
    f_fcstats,
    f_bcheck,
    f_msplitby,
    f_distance,
    f_sphcorr,
    f_strstr,
    f_instring,
    f_strftime,
    f_ldiff,
    f_sdc,
    f_cdemean,
    f_stdize,
    f_instrings,
    f_threshold,
    f_jsongetb,
    f_tsolvepd,
    f_kmeans,
    hf_regls,
    f3_max,  /* SEPARATOR: end of three-arg functions */
    f_urcpval,
    f_randgen,
    f_mrandgen,
    f_randgen1,
    f_pdf,
    f_pval,
    f_cdf,
    f_invcdf,
    f_crit,
    f_bkfilt,
    f_mols,
    f_mpols,
    f_mrls,
    f_fols,
    f_filter,
    f_mcovg,
    f_kfilter,
    f_ksmooth,
    f_kdsmooth,
    f_ksimul,
    f_nrmax,
    f_loess,
    f_ghk,
    f_quadtab,
    f_isoconv,
    f_qlrpval,
    f_bootci,
    f_bootpval,
    f_movavg,
    f_defarray,
    f_defbundle,
    f_deflist,
    f_defargs,
    f_defmat,
    f_ksetup,
    f_bfgscmax,
    f_svm,
    f_irf,
    f_nadarwat,
    f_feval,
    f_chowlin,
    f_tdisagg,
    f_hyp2f1,
    f_midasmult,
    f_commute,
    f_toepsolv,
    f_rgbmix,
    hf_felogitr,
    fn_max,  /* SEPARATOR: end of n-arg functions */}

enum ConstId {
    pi = 1,
    na,
    inf,
    nan,
    win32,
    eps,
    have_mpi,
    mpi_rank,
    mpi_size,
    n_proc,
    true_,
    false_,
    sysinfo
}

enum DumId {
    null_ = 1,
    empty,
    diag,
    upper,
    lower,
    real_,
    imag,
    end,
    dataset,
    trend
}

enum GENSTRLEN = 128;
enum NO_VNUM = -1;

bool unary_op(int s)  { return s >= 1 && s < GenSym.u_max; }
bool binary_op(int s) { return s > GenSym.u_max && s < GenSym.op_max; }
bool bool_comp(int s) { return s >= GenSym.b_eq && s <= GenSym.b_or; }
bool dot_op(int s)    { return s >= GenSym.b_dotmult && s <= GenSym.b_dotneq; }

bool func1_symb(int s) { return s > GenFunc.f1_min && s < GenFunc.f1_max; }
bool func2_symb(int s) { return s > GenFunc.f1_max && s < GenFunc.f2_max; }
bool func3_symb(int s) { return s > GenFunc.f2_max && s < GenFunc.f3_max; }
bool funcn_symb(int s) { return s > GenFunc.f3_max && s < GenFunc.fn_max; }

bool bnsym(int s) { return s == GenSym.fargs; }

bool alias_reversed(const(Node) * n) { return (n.flags & NodeFlags.als) != 0; }

/* function with single string argument */
bool string_arg_func(int s)
{
    return (s == GenFunc.f_isdiscr || s == GenFunc.f_obsnum ||
            s == GenFunc.f_backtick || s == GenFunc.f_varnum ||
            s == GenFunc.f_remove || s == GenFunc.f_exists);
}

/* function with multiple args, string for first arg */
bool str0_func(int s)
{
    return (s == GenFunc.f_pval || s == GenFunc.f_cdf || s == GenFunc.f_invcdf ||
            s == GenFunc.f_crit || s == GenFunc.f_randgen || s == GenFunc.f_pdf ||
            s == GenFunc.f_bessel || s == GenFunc.f_mrandgen || s == GenFunc.f_randgen1);
}

/* functions taking a string arg in last position */
bool string_last_func(int s)
{
    return (s == GenFunc.f_aggrby || s == GenFunc.f_genseries ||
            s == GenFunc.f_printf || s == GenFunc.f_sprintf ||
            s == GenFunc.f_allreduce || s == GenFunc.f_normtest ||
            s == GenFunc.f_sscanf || s == GenFunc.f_npcorr ||
            s == GenFunc.f_inbundle || s == GenFunc.f_asort);
}

/* functions taking string arg in middle position */
bool string_mid_func(int s)
{
    return (s == GenFunc.f_reduce || s == GenFunc.f_scatter);
}

/* functions taking one or more "fncall" (string) arguments */
bool fncall_func(int s)
{
    return (s == GenFunc.f_bfgsmax || s == GenFunc.f_nrmax ||
            s == GenFunc.f_fdjac || s == GenFunc.f_simann ||
            s == GenFunc.f_bfgscmax || s == GenFunc.f_nmmax ||
            s == GenFunc.f_gssmax || s == GenFunc.f_numhess ||
            s == GenFunc.f_fzero);
}

/* functions with "reversing" aliases */
bool als_func(int s)
{
    return (s == GenFunc.f_bfgsmax || s == GenFunc.f_nrmax ||
            s == GenFunc.f_simann || s == GenFunc.f_bfgscmax ||
            s == GenFunc.f_nmmax || s == GenFunc.f_gssmax ||
            s == GenFunc.f_exists);
}

/* functions where the right-hand argument is actually a return
   location */
bool r_return(int s)
{
    return (s == GenFunc.f_qr || s == GenFunc.f_eigsym || s == GenFunc.f_eigen ||
            s == GenFunc.f_mols || s == GenFunc.f_mpols || s == GenFunc.f_svd);
}

/* functions where the middle argument is actually a return
   location */
bool m_return(int s)
{
    return (s == GenFunc.f_svd || s == GenFunc.f_eigen || s == GenFunc.f_qr);
}

bool undef_arg_ok(int s)
{
    return (s == GenFunc.f_typeof || s == GenFunc.f_typename ||
            s == GenFunc.f_iscmplx);
}

bool reusable(const(Parser) * p)
{
    return (p.flags & (GenFlags.compile | GenFlags.exec)) != 0;
}

struct BranchN {
    int n_nodes;
    Node ** n;
}

union Val {
    BranchN bn;
    int idnum;
    char * str;
    double xval;
    double * xvec;
    int * ivec;
    GretlMatrixRaw * m;
    MatrixSubspec * mspec;
    GretlBundle * b;
    GretlArray * a;
    void * ptr;
}

enum NodeFlags {
    aux = 1 << 0,   /* AUX_NODE: auxiliary: free on exit */
    tmp = 1 << 1,   /* TMP_NODE: temporary: free content on exit */
    prx = 1 << 2,   /* PRX_NODE: aux node is proxy (don't reuse!) */
    lht = 1 << 3,   /* LHT_NODE: node holds terminal of LHS */
    mut = 1 << 4,   /* MUT_NODE: node is inherently mutable in type */
    als = 1 << 5    /* ALS_NODE: function subject to "reversing" alias */
}

struct Node {
    gint16 t;        /* type identifier */
    guint8 flags;    /* AUX_NODE etc., see above */
    int vnum;        /* associated series ID number */
    char * vname;    /* associated variable name */
    UserVar * uv;    /* associated named variable */
    Val v;           /* value (of whatever type) */
    Node * L;        /* up to three child nodes */
    Node * M;
    Node * R;
    Node * aux;      /* auxiliary (result) node */
    Node * parent;   /* parent node (or NULL) */
    int refcount;    /* reference counter, used by aux nodes */
}

enum GenFlags {
    discard = 1 <<  0, /* compute and print, don't save */
    start   = 1 <<  1, /* first round of evaluation */
    autoreg = 1 <<  2, /* expression is autoregressive */
    decl    = 1 <<  3, /* statement is actually a declaration */
    priv    = 1 <<  4, /* generating a "private" or internal var */
    compile = 1 <<  5, /* compiling the parse tree */
    exec    = 1 <<  6, /* evaluating a compiled tree */
    natest  = 1 <<  7, /* testing for NAs in expression */
    ufret   = 1 <<  8, /* returning value generated by user function */
    quiet   = 1 <<  9, /* don't print any messages or labels */
    getstr  = 1 << 10, /* state: flag acceptance of plain strings */
    mmask   = 1 << 11, /* genr result is masked matrix */
    slicing = 1 << 12, /* state: calculating object slice (temporary) */
    lagprse = 1 << 13, /* state: parsing lag spec (temporary) */
    deltan  = 1 << 14, /* flag for change in series length */
    catch_  = 1 << 15, /* "catch" is in force */
    nodecl  = 1 << 16, /* type of result was not specified */
    listdef = 1 << 17, /* expression defines a list */
    anon    = 1 << 18, /* generating an anonymous object */
    void_   = 1 << 19, /* function call, no assignment */
    noexec  = 1 << 20, /* just compile, don't evaluate */
    msave   = 1 << 21, /* trying for reuse of an aux matrix */
    obsval  = 1 << 22, /* generating value of observation in series */
    aliased = 1 << 23, /* state: handling aliased object (temporary) */
    and     = 1 << 24, /* state: working on right-hand term of B_AND */
    or      = 1 << 25, /* state: working on right-hand term of B_OR */
    stack   = 1 << 26, /* executing stack() */
    altinp  = 1 << 27, /* the input string has been substituted */
    objqry  = 1 << 28, /* querying the existence of an object */
    prnlist = 1 << 29  /* defining a list for "print" */
}

struct LhInfo {
    int t;                  /* type of pre-existing LHS variable, if any */
    char[VNAMELEN] name;    /* name of LHS variable */
    char * label;           /* descriptive string for series */
    SeriesTable * stab;     /* holds string values for series */
    int vnum;               /* ID number of pre-existing LHS series */
    UserVar * uv;           /* address of pre-existing LHS variable */
    char * expr;            /* expression on left */
    GretlType gtype;        /* gretl type of LHS array, if any, or
                               of LHS bundle member */
    GretlMatrixRaw * mret;     /* matrix output (possibly under bundle or array) */
}

struct Parser {
    const(char) * input; /* complete input string */
    const(char) * point; /* remaining unprocessed input */
    const(char) * rhs;   /* for use in labelling */
    DatasetRaw * dset;      /* convenience pointer to dataset */
    Prn * prn;           /* for printing messages */
    Prn * errprn;        /* for storing error message in case @prn is NULL */
    GenFlags flags;      /* various attributes (see @GenFlags above) */
    int targ;            /* target type */
    int op;              /* assignment operator (possibly inflected) */
    LhInfo lh;           /* left-hand side info */
    Node * lhtree;       /* LHS syntax tree, if needed */
    Node * lhres;        /* result of eval() on @lhtree */
    Node * tree;         /* RHS syntax tree */
    Node * ret;          /* result of eval() on @tree */
    /* below: parser state variables */
    Node * aux;          /* convenience pointer to current auxiliary node */
    int callcount;
    int dset_n;
    int obs;
    int sym;
    int upsym;
    int ch;
    double xval;
    int idnum;
    char * idstr;
    void * data;
    int err;
}
