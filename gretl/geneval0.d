/* The common head of lib/src/geneval.c (see geneval_private.h): macros, small enums and shared settings used by all
   five parts geneval1.d .. geneval5.d, which port consecutive slices of geneval.c (geneval1.c .. geneval5.c).

   Notes on the port of geneval.c as a whole:
   - The C file is a very large file, split here into five D modules; there is no need for `static` at module
     boundaries in D, so functions that are static in the C are public in the D parts (the macros and the module
     state below are shared the same way). Names are unique across the parts because they were one C file.
   - The EDEBUG, LHDEBUG, AUX_NODES_DEBUG, UVDEBUG and END_DEBUG diagnostics (all compiled out in the C: they are
     0) and the `#if 0` blocks are not ported.
   - Build-time conditionals: HAVE_MPI, USE_CURL, USE_RLIB and WIN32 are taken as undefined (as in the earlier
     ports), so the non-MPI, non-curl, non-R and non-Windows branches are the ones ported.
   - The function-like macros are functions here (the macro that takes a parser or a node takes the pointer). */

module gretl.geneval0;


import gretl.genparseh;
import gretl.gretl_matrixh : gretl_matrix_is_cscalar, gretl_matrix_is_scalar;
import gretl.libgretlh;

import gretl.dataset : is_string_valued;

enum AGGR_PREFER_MATRIX = 1;
enum ONE_BY_ONE_CAST = 1;

/* the "code" argument of free_tree() */
enum FreeCode {
    tree = 1,
    ret,
    lhtree,
    lhres,
    aret
}

bool is_aux_node(const(Node) * n) { return n != null && (n.flags & NodeFlags.aux); }
bool is_tmp_node(const(Node) * n) { return n != null && (n.flags & NodeFlags.tmp); }
bool is_proxy_node(const(Node) * n) { return n != null && (n.flags & NodeFlags.prx); }

bool emptymat_ok(int f)
{
    return (f == GenFunc.f_ginv || f == GenFunc.f_diag || f == GenFunc.f_transp ||
            f == GenFunc.f_vec || f == GenFunc.f_vech || f == GenFunc.f_unvech ||
            f == GenFunc.f_chol || f == GenFunc.f_upper || f == GenFunc.f_lower ||
            f == GenFunc.f_sort || f == GenFunc.f_dsort || f == GenFunc.f_values ||
            f == GenFunc.f_mrev);
}

bool dataset_dum(const(Node) * n) { return n.t == GenSym.dum && n.v.idnum == DumId.dataset; }

bool postfix_node(const(Node) * n) { return n.t == GenSym.num_p || n.t == GenSym.num_m; }

bool uscalar_node(const(Node) * n) { return (n.t == GenSym.num && n.vname != null) || postfix_node(n); }

bool umatrix_node(const(Node) * n) { return n.t == GenSym.mat && n.vname != null; }
bool ubundle_node(const(Node) * n) { return n.t == GenSym.bundle && n.vname != null; }
bool uarray_node(const(Node) * n)  { return n.t == GenSym.array && n.vname != null; }
bool ustring_node(const(Node) * n) { return n.t == GenSym.str && n.vname != null; }
bool ulist_node(const(Node) * n)   { return n.t == GenSym.list && n.vname != null; }
bool useries_node(const(Node) * n) { return n.t == GenSym.series && n.vnum >= 0; }

bool uvar_node(const(Node) * n) { return n.vname != null; }
bool func_node(const(Node) * n)
{
    return n.t == GenSym.ufun || (n.t > GenFunc.f1_min && n.t < GenFunc.fn_max);
}

bool scalar_matrix_node(const(Node) * n) { return n.t == GenSym.mat && gretl_matrix_is_scalar(n.v.m); }
bool scalar_node(const(Node) * n) { return n.t == GenSym.num || scalar_matrix_node(n); }
bool ok_matrix_node(const(Node) * n) { return n.t == GenSym.mat || n.t == GenSym.num; }
bool complex_node(const(Node) * n) { return n.t == GenSym.mat && n.v.m.is_complex; }
bool cscalar_node(const(Node) * n) { return n.t == GenSym.mat && gretl_matrix_is_cscalar(n.v.m); }

bool strvals_node(const(Node) * n, const(Parser) * p)
{
    return n.vnum > 0 && is_string_valued(p.dset, n.vnum);
}

bool mutable_node(const(Node) * n) { return (n.flags & NodeFlags.mut) != 0; }

bool null_node(const(Node) * n) { return n == null || n.t == GenSym.empty; }
bool null_or_scalar(const(Node) * n) { return null_node(n) || scalar_node(n); }
bool null_or_string(const(Node) * n) { return null_node(n) || n.t == GenSym.str; }

bool ok_bundled_type(int t)
{
    return (t == GenSym.num || t == GenSym.str || t == GenSym.mat || t == GenSym.list ||
            t == GenSym.series || t == GenSym.bundle || t == GenSym.array);
}

bool compiled(const(Parser) * p) { return (p.flags & GenFlags.exec) != 0; }
bool starting(const(Parser) * p) { return (p.flags & GenFlags.start) != 0; }
bool autoreg(const(Parser) * p)  { return (p.flags & GenFlags.autoreg) != 0; }
bool exestart(const(Parser) * p)
{
    return (p.flags & GenFlags.exec) && (p.flags & GenFlags.start);
}

/* set by set_user_qsorting() (geneval1.d) and used in geneval5.d */
__gshared int user_qsorting;

/* the MIN and MAX macros of libgretl.h */
T MIN(T)(T a, T b) { return a < b ? a : b; }
T MAX(T)(T a, T b) { return a > b ? a : b; }

/* the types of gretl_func.c (gretl_func0.d): parts of geneval only ever pass them around by pointer */
public import gretl.gretl_func0 : Ufunc, Fncall, Fnpkg;
