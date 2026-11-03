/* Port of gretl_bfgs.h: the MaxMethod enum and the callback types (the prototypes are not ported here; the C
   functions are declared extern (C) where used until gretl_bfgs.c is ported). The enum members keep their
   C names in lower case (BHHH_MAX -> bhhh_max, ROOT_FIND -> root_find) since they do not share a prefix. */

module gretl.gretl_bfgsh;


import gretl.gretl_matrixh : GretlMatrixRaw;

enum MaxMethod {
    bhhh_max,
    bfgs_max,
    lbfgs_max,
    simann_max,
    nm_max,
    gss_max,
    root_find
}

alias BFGS_CRIT_FUNC = extern (C) double function(const(double) *, void *);
alias BFGS_GRAD_FUNC = extern (C) int function(double *, double *, int, BFGS_CRIT_FUNC, void *);
alias BFGS_COMBO_FUNC = extern (C) double function(double *, double *, int, void *);
alias BFGS_LLT_FUNC = extern (C) const(double) * function(const(double) *, int, void *);
alias HESS_FUNC = extern (C) int function(double *, GretlMatrixRaw *, void *);
alias ZFUNC = extern (C) double function(double, void *);
