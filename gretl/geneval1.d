/* Port of the first part (lines 177-5001) of lib/src/geneval.c: geneval1.c. See geneval0.d for the notes that
   apply to all five parts of the file. */

module gretl.geneval1;


import core.stdc.errno : errno;
import core.stdc.limits : INT_MAX;
import core.stdc.math : fmod, pow;
import core.stdc.stdio : fprintf, fputs, sprintf, stderr;
import core.stdc.stdlib : calloc, free, malloc;
import core.stdc.string : memcpy, memset, strcat, strchr, strcmp, strcpy, strlen;

import gretl.genparseh;
import gretl.geneval0;
import gretl.gretl_errorsh : GretlError;
import gretl.libgretlh;
import gretl.glib : guint32;
import gretl.gretl_matrixh;
import gretl.usermath : MatrixSubspec;
import gretl.uservarh : UserVar;
import gretl.gretl_bfgsh : MaxMethod;
import gretl.adf_kpss : get_urc_pvalue;
import gretl.dataio : dateton, get_date_x;
import gretl.dataset : dataset_add_series, is_string_valued, series_decode_string, series_get_string_for_obs, series_get_string_vals;
import gretl.dataseth : annual_data, dataset_is_panel, sample_size;
import gretl.describe : distance;
import gretl.geneval2 : node_get_list;
import gretl.geneval4 : node_type_error;
import gretl.genmain : caller_series_index, current_series_index;
import gretl.gretl_errors : check_gretl_warning, errmsg, gretl_errmsg_is_set, gretl_errmsg_set, gretl_errmsg_sprintf, gretl_error_clear, gretl_strerror, gretl_warnmsg_set, gretl_warnmsg_sprintf;
import gretl.gretl_list : gretl_list_new;
import gretl.gretl_matrix : gretl_cholesky_solve, gretl_cmatrix_new, gretl_covariance_matrix, gretl_invert_matrix, gretl_matrix_add, gretl_matrix_alloc, gretl_matrix_bool_sel, gretl_matrix_cholesky_decomp, gretl_matrix_col_concat, gretl_matrix_complex_divide, gretl_matrix_complex_multiply, gretl_matrix_cond_index, gretl_matrix_copy, gretl_matrix_copy_transpose, gretl_matrix_copy_values, gretl_matrix_cumcol, gretl_matrix_diag_qform, gretl_matrix_diffcol, gretl_matrix_direct_sum, gretl_matrix_divide, gretl_matrix_dot_op, gretl_matrix_exp, gretl_matrix_free, gretl_matrix_get_colnames, gretl_matrix_get_diagonal, gretl_matrix_get_rownames, gretl_matrix_get_t1, gretl_matrix_get_t2, gretl_matrix_hdproduct_new, gretl_matrix_infinity_norm, gretl_matrix_init, gretl_matrix_is_dated, gretl_matrix_is_symmetric, gretl_matrix_isfinite, gretl_matrix_kronecker_product_new, gretl_matrix_left_nullspace, gretl_matrix_minmax, gretl_matrix_moore_penrose, gretl_matrix_multiply, gretl_matrix_multiply_mod, gretl_matrix_one_norm, gretl_matrix_polroots, gretl_matrix_pow, gretl_matrix_psd_root, gretl_matrix_qform, gretl_matrix_rank, gretl_matrix_rcond, gretl_matrix_realloc, gretl_matrix_reverse_cols, gretl_matrix_reverse_rows, gretl_matrix_right_nullspace, gretl_matrix_row_concat, gretl_matrix_set_complex_full, gretl_matrix_set_equals_tolerance, gretl_matrix_set_t1, gretl_matrix_set_t2, gretl_matrix_sort_by_column, gretl_matrix_subtract, gretl_matrix_trace, gretl_matrix_transcribe_obs_info, gretl_matrix_transpose_in_place, gretl_matrix_unset_equals_tolerance, gretl_matrix_zero_lower, gretl_matrix_zero_upper, gretl_null_matrix_new, gretl_unit_matrix_new, gretl_zero_matrix_new;
import gretl.gretl_prn : Prn;
import gretl.gretl_string_table : SeriesTable, series_table_new;
import gretl.gretl_utils : bin2dec, dec2bin, gretl_int_from_double, gretl_permute, gretl_uint32_from_double;
import gretl.gstrfuncs : g_utf8_offset_to_pointer, g_utf8_strlen;
import gretl.libset : gretl_iteration_pop, gretl_iteration_push;
import gretl.missingh : NADBL, na;
import gretl.plugins : get_plugin_function;
import gretl.pvalues : dist_code_from_string, gretl_fill_pdf_array, gretl_fill_random_series, gretl_get_DW, gretl_get_cdf, gretl_get_cdf_inverse, gretl_get_critval, gretl_get_pdf, gretl_get_pvalue, gretl_get_random_matrix, gretl_get_random_scalar;
import gretl.pvaluesh : DistCode;
import gretl.strutils : gretl_strdup, gretl_string_split, strings_array_free, strings_array_new;
import gretl.uservarh : uvar_set_scalar_fast;
import gretl.matrix_extra : gretl_matrix_2d_convolution, gretl_matrix_which, gretl_matrix_write_to_file;
import gretl.gretl_cmatrix : gretl_cmatrix_AHB, gretl_cmatrix_add_sub, gretl_cmatrix_cholesky, gretl_cmatrix_determinant, gretl_cmatrix_divide, gretl_cmatrix_dot_op, gretl_cmatrix_exp, gretl_cmatrix_fft, gretl_cmatrix_ginv, gretl_cmatrix_hdprod, gretl_cmatrix_inverse, gretl_cmatrix_kronecker, gretl_cmatrix_multiply, gretl_cmatrix_rank, gretl_cmatrix_switch, gretl_cmatrix_trace, gretl_ctrans, gretl_ctrans_in_place, gretl_matrix_fft, gretl_matrix_ffti, gretl_matrix_log;
import gretl.discrete : bit_permutations;
import gretl.gretl_bfgs : deriv_free_optimize, user_BFGS, user_fdjac, user_numhess;
import gretl.gretl_normal : bvnorm_cdf;
import gretl.subsample : dataset_is_subsampled;
import gretl.gretl_array : get_array_by_name, gretl_array_copy, gretl_array_destroy, gretl_array_from_strings, gretl_array_get_data, gretl_array_get_length, gretl_array_get_type, gretl_matrix_split_by;
import gretl.gretl_bundle : get_bundle_by_name, gretl_bundle_destroy, gretl_bundle_write_to_file;
import gretl.uservar : get_list_by_name, get_string_by_name, get_user_var_by_name, get_user_var_of_type_by_name, gretl_scalar_set_value, user_list_append, user_list_replace, user_list_subtract, user_var_add, user_var_add_or_replace, user_var_get_type, user_var_replace_value;
import gretl.gretl_typemap : gretl_is_arrayable_type, gretl_is_scalar_type, gretl_type_get_name;
import gretl.genfuncs : gretl_matrix_vector_stat, imhof, rank_vector;
import gretl.usermat : get_matrix_by_name, umatrix_set_names_from_array, umatrix_set_names_from_list, umatrix_set_names_from_string, user_matrix_get_column_name, user_matrix_get_determinant, user_matrix_get_row_name, user_matrix_unvech, user_matrix_vec, user_matrix_vech;
import gretl.random : gretl_rand_dirichlet, gretl_rand_discrete;

/* unported dependencies: per-file declarations from the 2026c headers */
extern (C) {
    int const_lookup(const(char) * s);
    const(char) * getsymb(int t);
    int gretl_n_processors();
    Node * new_node(int t);
    Node * newdbl(double x);
    Node * newempty();
    void undefined_symbol_error(const(char) * s, Parser * p);
}

void set_user_qsorting(int s)
{
    if (s != 0) {
        user_qsorting = 1;
        gretl_iteration_push();
    } else {
        user_qsorting = 0;
        gretl_iteration_pop();
    }
}

/* ok_list_node: This is a first-pass assessment of whether a given
   node _may_ be interpretable as holding a LIST.  The follow-up is
   node_get_list(), which will determine whether the interpretation
   really works. The code below dates from 2023-05-11, when we
   reverted to an earlier version [a9ec40]. In the interim, version
   [b6ed06] had "if (null_node(n))" as the first test and "return
   scalar_node(n)" as the consequent of the last test, giving a
   somewhat different take on the case where @n is an EMPTY node (as
   opposed to NULL). It seems that the interim version broke some
   things.
*/

int ok_list_node(Node * n, Parser * p)
{
    if (n == null) {
        return 0;
    } else if (n.t == GenSym.list) {
        return 1;
    } else if (n.t == GenSym.series && n.vnum >= 0) {
        /* can interpret as singleton list */
        return 1;
    } else if (p.flags & GenFlags.listdef) {
        /* when defining a list we can be a bit more accommodating */
        return null_or_scalar(n);
    }

    return 0;
}

/* a more "lenient" version of the above, to accommodate list
   expressions such as (L - 0), indicating the list that results from
   dropping the constant from L
*/

int ok_list_node_plus(Node * n)
{
    if (n.t == GenSym.list) {
        return 1;
    } else if (n.t == GenSym.series && n.vnum >= 0) {
        return 1;
    } else if (n.t == GenSym.num) {
        return 1;
    } else {
        return 0;
    }
}

const(char) * typestr(int t)
{
    switch (t) {
    case GenSym.num:
        return "scalar";
    case GenSym.series:
        return "series";
    case GenSym.mat:
        return "matrix";
    case GenSym.str:
        return "string";
    case GenSym.u_addr:
        return "address";
    case GenSym.list:
        return "list";
    case GenSym.bundle:
        return "bundle";
    case GenSym.dbundle:
        return "$-bundle";
    case GenSym.array:
        return "array";
    case GenSym.useries:
        return "named series";
    case GenSym.empty:
        return "empty";
    case GenSym.undef:
        return "undefined";
    default:
        return "?";
    }
}

void free_mspec(MatrixSubspec * spec, Parser * p)
{
    if (spec != null) {
        free(spec.rslice);
        free(spec.cslice);
        free(spec);
    }
}

void clear_mspec(MatrixSubspec * spec, Parser * p)
{
    free(spec.rslice);
    free(spec.cslice);

    memset(spec, 0, (*spec).sizeof);
}

/* return 1 if @n equals @t or is found under @t, otherwise 0 */

int find_in_tree(Node * t, Node * n)
{
    if (t == null) {
        return 0;
    } else if (t == n) {
        return 1;
    } else if (find_in_tree(t.L, n)) {
        return 1;
    } else if (find_in_tree(t.R, n)) {
        return 1;
    } else if (find_in_tree(t.M, n)) {
        return 1;
    } else if (bnsym(t.t)) {
        int i;

        for (i=0; i<t.v.bn.n_nodes; i++) {
            if (find_in_tree(t.v.bn.n[i], n)) {
                return 1;
            }
        }
    }

    return 0;
}

/* used when we know that @t is a terminal node: skip
   the tests for attached tree */

void free_node(Node * t, Parser * p)
{
    if (t.refcount > 1) {
        t.refcount -= 1;
        return;
    }

    if (is_tmp_node(t)) {
        if (t.t == GenSym.series) {
            free(t.v.xvec);
        } else if (t.t == GenSym.list || t.t == GenSym.ivec) {
            free(t.v.ivec);
        } else if (t.t == GenSym.mat) {
            gretl_matrix_free(t.v.m);
        } else if (t.t == GenSym.mspec) {
            free_mspec(t.v.mspec, p);
        } else if (t.t == GenSym.bundle) {
            gretl_bundle_destroy(t.v.b);
        } else if (t.t == GenSym.array) {
            gretl_array_destroy(t.v.a);
        } else if (t.t == GenSym.str) {
            free(t.v.str);
        }
    }

    if (t.t == GenSym.uobj || t.t == GenSym.wlist) {
        free(t.v.str);
    }

    if (t.vname != null) {
        free(t.vname);
    }

    if (p != null && t == p.ret) {
        p.ret = null;
    }

    free(t);
}

/* A word on "aux" nodes. These come in two sorts, which
   might be described as "robust" and "fragile" respectively.

   A robust node (identified by the TMP_NODE flag) is one
   whose data pointer is independently allocated. With such
   a node it's OK simply to "pass on" the pointer in
   assignment, and if it's not passed on it should be freed
   on completion of "genr". (So nota bene: if it's assigned
   elsewhere, the pointer on the aux node itself must then
   be set to NULL to avoid double-freeing.)

   A fragile node is one whose data pointer is not
   independently allocated; it actually "belongs to someone
   else". In assignment, then, it must be deeply copied,
   and it must _not_ be freed on completion of genr.

   Obviously, it's necessary to be careful in handling
   fragile nodes, but the advantage of allowing them
   is that they cut down on wasteful deep-copying of objects
   that may be used in calculation, without being modified,
   on the fly.
*/

void free_tree(Node * t, Parser * p, int code)
{
    if (t == null) {
        return;
    }

    /* free recursively */
    if (bnsym(t.t)) {
        int i;

        for (i=0; i<t.v.bn.n_nodes; i++) {
            free_tree(t.v.bn.n[i], p, code);
        }
        free(t.v.bn.n);
    } if (!(t.flags & NodeFlags.lht)) {
        free_tree(t.L, p, code);
        free_tree(t.M, p, code);
        free_tree(t.R, p, code);
    }

    if (t.aux != null && t.aux != p.ret && t.aux != p.lhres) {
        free_node(t.aux, p);
    } else if (t.aux != null) {
        t.aux.refcount -= 1;
    }

    free_node(t, p);
}

void clear_uvnodes(Node * t)
{
    if (t == null) {
        return;
    }

    if (bnsym(t.t)) {
        int i;

        for (i=0; i<t.v.bn.n_nodes; i++) {
            clear_uvnodes(t.v.bn.n[i]);
        }
    } else {
        clear_uvnodes(t.L);
        clear_uvnodes(t.M);
        clear_uvnodes(t.R);
    }

    if (t.t == GenSym.series) {
        if (t.vnum >= 0 || t.vname != null) {
            t.v.xvec = null;
        }
    } else if (t.uv != null) {
        t.uv = null;
    }
}

void reset_p_aux(Parser * p, Node * n)
{
    p.aux = n;
}

double * na_array(int n)
{
    double * x = cast(double *) malloc(n * double.sizeof);
    int i;

    if (x != null) {
        for (i=0; i<n; i++) {
            x[i] = NADBL;
        }
    }

    return x;
}

/* new node to hold array of doubles */

Node * newseries(int n, int flags)
{
    Node * b = new_node(GenSym.series);

    if (b != null) {
        b.flags = cast(ubyte) flags;
        if (n > 0) {
            b.v.xvec = na_array(n);
            if (b.v.xvec == null) {
                free(b);
                b = null;
            }
        } else {
            b.v.xvec = null;
        }
    }

    return b;
}

/* new node to hold array of @n ints */

Node * newivec(int n)
{
    Node * b = new_node(GenSym.ivec);

    if (b != null) {
        b.flags = NodeFlags.tmp;
        if (n > 0) {
            b.v.ivec = cast(int *) malloc(n * int.sizeof);
            if (b.v.ivec == null) {
                free(b);
                b = null;
            }
        } else {
            b.v.ivec = null;
        }
    }

    return b;
}

/* new node to hold a gretl_matrix */

Node * newmat(int flags)
{
    Node * n = new_node(GenSym.mat);

    if (n != null) {
        n.flags = cast(ubyte) flags;
        n.v.m = null;
    }

    return n;
}

/* new node to hold a matrix specification */

Node * newmspec()
{
    Node * n = new_node(GenSym.mspec);

    if (n != null) {
        n.flags = NodeFlags.tmp;
        n.v.mspec = null;
    }

    return n;
}

/* new node to hold a list */

Node * newlist(int flags)
{
    Node * n = new_node(GenSym.list);

    if (n != null) {
        n.flags = cast(ubyte) flags;
        n.v.ivec = null;
    }

    return n;
}

Node * newstring(int flags)
{
    Node * n = new_node(GenSym.str);

    if (n != null) {
        n.flags = cast(ubyte) flags;
        n.v.str = null;
    }

    return n;
}

Node * newbundle(int flags)
{
    Node * n = new_node(GenSym.bundle);

    if (n != null) {
        n.flags = cast(ubyte) flags;
        n.v.b = null;
    }

    return n;
}

Node * newarray(int flags)
{
    Node * n = new_node(GenSym.array);

    if (n != null) {
        n.flags = cast(ubyte) flags;
        n.v.a = null;
    }

    return n;
}

void clear_tmp_node_data(Node * n, Parser * p)
{
    int nullify = 1;

    if (n.t == GenSym.list) {
        free(n.v.ivec);
    } else if (n.t == GenSym.mat) {
        /* (how) can we avoid doing this? */
        gretl_matrix_free(n.v.m);
    } else if (n.t == GenSym.mspec) {
        if (n.v.mspec != null) {
            clear_mspec(n.v.mspec, p);
        }
        nullify = 0;
    } else if (n.t == GenSym.bundle) {
        gretl_bundle_destroy(n.v.b);
    } else if (n.t == GenSym.array) {
        gretl_array_destroy(n.v.a);
    } else if (n.t == GenSym.str) {
        free(n.v.str);
    } else if (n.t == GenSym.series) {
        /* preserve any existing tmp series, unless the
           dataset series length has changed
        */
        if (p.flags & GenFlags.deltan) {
            free(n.v.xvec);
            n.v.xvec = null;
            if (p.dset_n > 0) {
                n.v.xvec = na_array(p.dset_n);
                if (n.v.xvec == null) {
                    p.err = GretlError.alloc;
                }
            }
        } else {
            /* scrub any pre-existing values in the current
               sample range */
            int t;

            for (t=p.dset.t1; t<=p.dset.t2; t++) {
                n.v.xvec[t] = NADBL;
            }
        }
        nullify = 0;
    } else {
        nullify = 0;
    }

    if (nullify) {
        n.v.ptr = null;
    }
}

int mutate_bundle_member_node(Node * n, int type,
                              int flags, Parser * p)
{
    int err = 0;

    if (is_tmp_node(n)) {
        /* some allocated storage should be freed */
        if (n.t == GenSym.series) {
            free(n.v.xvec);
        } else if (n.t == GenSym.mat) {
            gretl_matrix_free(n.v.m);
        }
        n.v.ptr = null;
    }

    if (type == GenSym.series) {
        /* switching to a series node : allocate xvec */
        n.v.xvec = cast(double *) malloc(p.dset.n * double.sizeof);
        if (n.v.xvec == null) {
            err = GretlError.alloc;
        }
    }

    if (!err) {
        n.t = cast(short) type;
        n.flags = cast(ubyte) flags;
    }

    return err;
}

/* We allow here for some equivocation in type between
   1 x 1 matrices and scalars in the course of executing
   a compiled parser with saved aux nodes.
*/

void maybe_switch_node_type(Node * n, int type,
                            int flags, Parser * p)
{
    if (mutable_node(n)) {
        /* bundle members only */
        p.err = mutate_bundle_member_node(n, type, flags, p);
    } else if (n.t == GenSym.mat && type == GenSym.num) {
        /* switch aux node @n from matrix to scalar */
        if (is_tmp_node(n)) {
            gretl_matrix_free(n.v.m);
        }
        n.t = GenSym.num;
        n.v.xval = NADBL;
        n.flags = 0;
        n.vnum = NO_VNUM;
        n.vname = null;
    } else if (n.t == GenSym.num && type == GenSym.mat) {
        /* switch @n from scalar to matrix */
        n.t = GenSym.mat;
        n.v.m = null;
        n.flags = cast(ubyte) flags;
    } else if (type == GenSym.empty) {
        /* LHS mechanism: OK */
    } else if (gen_type_is_arrayable(n.t) && type == GenSym.array) {
        n.t = GenSym.array;
        n.flags = cast(ubyte) flags;
    } else {
        /* any other discrepancy presumably means that
           things have gone badly wrong
        */
        fprintf(stderr, "aux node mismatch: expected %s but got %s, TMP_NODE=%d\n",
                getsymb(n.t), getsymb(type), (flags & NodeFlags.tmp)? 1 : 0);
        gretl_errmsg_set("internal genr error: aux node mismatch");
        p.err = GretlError.data;
    }
}

/* get an auxiliary node: if starting from scratch we allocate
   a new node, otherwise we look up an existing one */

Node * get_aux_node(Parser * p, int t, int n, int flags)
{
    Node * ret = p.aux;

    if (is_proxy_node(ret)) {
        /* this node will get freed later */
        ret = null;
    }

    if (ret != null) {
        /* got a pre-existing aux node */
        if (starting(p)) {
            if (ret.t != t) {
                maybe_switch_node_type(ret, t, flags, p);
            } else if (is_tmp_node(ret) && !(p.flags & GenFlags.msave)) {
                clear_tmp_node_data(ret, p);
            }
        }
    } else {
        /* we need to create a new aux node */
        if (t == GenSym.num) {
            ret = newdbl(NADBL);
        } else if (t == GenSym.series) {
            ret = newseries(n, flags);
        } else if (t == GenSym.ivec) {
            ret = newivec(n);
        } else if (t == GenSym.list) {
            ret = newlist(flags);
        } else if (t == GenSym.mat) {
            ret = newmat(flags);
        } else if (t == GenSym.mspec) {
            ret = newmspec();
        } else if (t == GenSym.str) {
            ret = newstring(flags);
        } else if (t == GenSym.bundle) {
            ret = newbundle(flags);
        } else if (t == GenSym.array) {
            ret = newarray(flags);
        } else if (t == GenSym.empty) {
            ret = newempty();
        } else {
            /* invalid aux node spec */
            p.err = GretlError.data;
        }

        if (!p.err && ret == null) {
            p.err = GretlError.alloc;
        }

        if (!p.err) {
            ret.flags |= NodeFlags.aux;
        }
    }

    return ret;
}

int e_types(Node * n)
{
    if (n != null) {
        GretlType t;
        const(char) * s;

        if (n.t == GenSym.array) {
            t = gretl_array_get_type(n.v.a);
            s = gretl_type_get_name(t);
        } else if (n.t == GenSym.u_addr) {
            s = "pointer-to-object";
        } else if (n.t == GenSym.empty) {
            s = "null";
        } else {
            t = gretl_type_from_gen_type(n.t);
            s = gretl_type_get_name(t);
            if (!strcmp(s, "invalid")) {
                s = getsymb(n.t);
            }
        }
        gretl_errmsg_sprintf("%s is %s, not acceptable in context",
                             n.vname == null ? "[anonymous]".ptr : n.vname, s);
    }

    return GretlError.types;
}

/* We come here by preference to the generic get_aux_node()
   (above) if we want an aux node holding an allocated matrix
   of known size (m x n). On the second or subsequent
   iterations of a loop, with any luck we may find that the
   aux node already holds a matrix of the required dimensions
   which can then be reused.
*/

Node * aux_sized_matrix_node(Parser * p, int m, int n,
                             int cmplx)
{
    Node * ret = p.aux;

    if (is_proxy_node(ret)) {
        /* this node will get freed later */
        ret = null;
    }

    if (ret != null) {
        /* got a pre-existing node */
        if (ret.t == GenSym.num) {
            /* switch @ret from scalar to matrix */
            ret.t = GenSym.mat;
            ret.v.m = null;
            ret.flags |= NodeFlags.tmp;
        } else if (ret.t != GenSym.mat) {
            p.err = e_types(ret);
        } else {
            /* check for reusable matrix */
            GretlMatrixRaw * a = ret.v.m;

            if (a != null) {
                if (a.is_complex + cmplx == 1) {
                    /* too difficult to reuse */
                    gretl_matrix_free(ret.v.m);
                    ret.v.m = null;
                } else if (a.rows != m || a.cols != n) {
                    p.err = gretl_matrix_realloc(ret.v.m, m, n);
                }
            }
        }
    } else {
        /* we need to create a new node */
        ret = newmat(NodeFlags.tmp | NodeFlags.aux);
        if (ret == null) {
            p.err = GretlError.alloc;
        }
    }

    if (!p.err && ret.v.m == null) {
        if (cmplx) {
            ret.v.m = gretl_cmatrix_new(m, n);
        } else {
            ret.v.m = gretl_matrix_alloc(m, n);
        }
        if (ret.v.m == null) {
            p.err = GretlError.alloc;
        }
    }

    return ret;
}

int no_data_error(Parser * p)
{
    p.err = GretlError.nodata;
    return GretlError.nodata;
}

Node * aux_series_node(Parser * p)
{
    if (p.dset == null || p.dset.n == 0) {
        no_data_error(p);
        return null;
    } else {
        return get_aux_node(p, GenSym.series, p.dset.n, NodeFlags.tmp);
    }
}

Node * aux_empty_series_node(Parser * p)
{
    if (p.dset == null || p.dset.n == 0) {
        no_data_error(p);
        return null;
    } else {
        return get_aux_node(p, GenSym.series, 0, NodeFlags.tmp);
    }
}

Node * aux_list_node(Parser * p)
{
    if (p.dset == null || p.dset.n == 0) {
        no_data_error(p);
        return null;
    } else {
        return get_aux_node(p, GenSym.list, 0, NodeFlags.tmp);
    }
}

Node * aux_scalar_node(Parser * p)        { return get_aux_node(p, GenSym.num, 0, 0); }
Node * aux_ivec_node(Parser * p, int n)   { return get_aux_node(p, GenSym.ivec, n, NodeFlags.tmp); }
Node * aux_matrix_node(Parser * p)        { return get_aux_node(p, GenSym.mat, 0, NodeFlags.tmp); }
Node * matrix_pointer_node(Parser * p)    { return get_aux_node(p, GenSym.mat, 0, 0); }
Node * aux_mspec_node(Parser * p)         { return get_aux_node(p, GenSym.mspec, 0, NodeFlags.tmp); }
Node * aux_string_node(Parser * p)        { return get_aux_node(p, GenSym.str, 0, NodeFlags.tmp); }
Node * string_pointer_node(Parser * p)    { return get_aux_node(p, GenSym.str, 0, 0); }
Node * aux_bundle_node(Parser * p)        { return get_aux_node(p, GenSym.bundle, 0, NodeFlags.tmp); }
Node * bundle_pointer_node(Parser * p)    { return get_aux_node(p, GenSym.bundle, 0, 0); }
Node * aux_array_node(Parser * p)         { return get_aux_node(p, GenSym.array, 0, NodeFlags.tmp); }
Node * array_pointer_node(Parser * p)     { return get_aux_node(p, GenSym.array, 0, 0); }
Node * aux_parent_node(Parser * p)        { return get_aux_node(p, GenSym.empty, 0, 0); }
Node * aux_any_node(Parser * p)           { return get_aux_node(p, 0, 0, 0); }
Node * aux_empty_node(Parser * p)         { return get_aux_node(p, GenSym.empty, 0, 0); }

/* Start of functions that probably should not be needed in
   their present full form, but testing is required before
   they're slimmed down. The general idea is that we should
   already have the user_var pointer we're in need of without
   having to look up it by name, again. However, we'll fall
   back to name look-up (and squawk about it on stderr) if
   need be.
*/

void * gen_get_lhs_var(Parser * p, GretlType type)
{
    void * data = null;

    if (p.lh.uv != null && p.lh.uv.type == type) {
        data = p.lh.uv.ptr;
    } else {
        if (p.lh.uv == null) {
            fprintf(stderr, "*** get: LHS %s '%s' is NULL!\n",
                    gretl_type_get_name(type), p.lh.name.ptr);
        } else {
            fprintf(stderr, "*** get: LHS uv '%s' of wrong type!\n",
                    p.lh.name.ptr);
        }
        if (type == GretlType.bundle) {
            data = get_bundle_by_name(p.lh.name.ptr);
        } else if (type == GretlType.array) {
            data = get_array_by_name(p.lh.name.ptr);
        } else if (type == GretlType.matrix) {
            data = get_matrix_by_name(p.lh.name.ptr);
        } else if (type == GretlType.string_) {
            data = get_string_by_name(p.lh.name.ptr);
        } else if (type == GretlType.list) {
            data = get_list_by_name(p.lh.name.ptr);
        }
    }

    return data;
}

struct TypeConv {
    GretlType t;
    int gen_t;
}

__gshared TypeConv[8] conversions = [
    TypeConv(GretlType.none,    GenSym.empty),
    TypeConv(GretlType.double_, GenSym.num),
    TypeConv(GretlType.series,  GenSym.series),
    TypeConv(GretlType.matrix,  GenSym.mat),
    TypeConv(GretlType.list,    GenSym.list),
    TypeConv(GretlType.string_, GenSym.str),
    TypeConv(GretlType.bundle,  GenSym.bundle),
    TypeConv(GretlType.array,   GenSym.array)
];

int gen_type_from_gretl_type(GretlType t)
{
    int i, n = cast(int) conversions.length;

    for (i=0; i<n; i++) {
        if (t == conversions[i].t) {
            return conversions[i].gen_t;
        }
    }

    if (gretl_is_scalar_type(t)) {
        return GenSym.num;
    }

    return GenSym.undef;
}

Node * aux_node_for_type(GretlType t, Parser * p)
{
    int gen_t = gen_type_from_gretl_type(t);

    return get_aux_node(p, gen_t, 0, 0);
}

GretlType gretl_type_from_gen_type(int gen_t)
{
    int i, n = cast(int) conversions.length;

    for (i=0; i<n; i++) {
        if (gen_t == conversions[i].gen_t) {
            return conversions[i].t;
        }
    }

    return GretlType.none;
}

int gen_type_is_arrayable(int gen_t)
{
    return gretl_is_arrayable_type(gretl_type_from_gen_type(gen_t));
}

int type_func_in_parentage(Node * n)
{
    while (n != null) {
        if (undef_arg_ok(n.t)) {
            return 1;
        }
        n = n.parent;
    }

    return 0;
}

Node * maybe_rescue_undef_node(Node * n, Parser * p)
{
    int v = current_series_index(p.dset, n.vname);
    UserVar * uv = null;

    if (v >= 0) {
        n.t = GenSym.series;
        n.vnum = v;
        n.v.xvec = p.dset.Z[v];
    } else if ((uv = get_user_var_by_name(n.vname)) != null) {
        GretlType type = user_var_get_type(uv);

        n.t = cast(short) gen_type_from_gretl_type(type);
        n.uv = uv;
        if (type == GretlType.double_) {
            n.v.xval = *cast(double *) uv.ptr;
        } else {
            n.v.ptr = uv.ptr;
        }
    } else {
        undefined_symbol_error(n.vname, p);
    }

    return n;
}

Node * revive_undef_node(Node * n, Parser * p)
{
    p.err = 0;
    gretl_error_clear();

    if (n == null) {
        n = aux_empty_node(p);
    } else {
        if (n.vname != null) {
            free(n.vname);
            n.vname = null;
        }
        n.t = GenSym.empty;
    }

    return n;
}

int gen_add_or_replace(Parser * p, GretlType type, void * data)
{
    int err;

    if (p.lh.uv != null) {
        err = user_var_replace_value(p.lh.uv, data, type);
    } else {
        err = user_var_add_or_replace(p.lh.name.ptr, type, data);
    }

    return err;
}

int gen_replace_lhs(Parser * p, GretlType type, void * data)
{
    if (p.lh.uv == null) {
        fputs("*** gen_replace_lhs: lhs user_var is NULL ***\n", stderr);
        fprintf(stderr, " (type is specified as %s)\n",
                gretl_type_get_name(type));
        return GretlError.data;
    } else {
        return user_var_replace_value(p.lh.uv, data, type);
    }
}

int gen_add_uvar(Parser * p, GretlType type, void * data)
{
    int err;

    err = user_var_add(p.lh.name.ptr, type, data);

    /* FIXME attach lh.uv pointer? */
    return err;
}

int gen_edit_list(Parser * p, int * list, int op)
{
    UserVar * u;
    int err;

    if (p.lh.uv != null && p.lh.uv.type == GretlType.list) {
        u = p.lh.uv;
    } else {
        if (p.lh.uv == null) {
            fprintf(stderr, "*** replace list: LHS uv is NULL!\n");
        } else {
            fprintf(stderr, "*** replace list: LHS uv of wrong type!\n");
        }
        u = get_user_var_of_type_by_name(p.lh.name.ptr, GretlType.list);
    }

    if (op == GenSym.b_asn) {
        err = user_list_replace(u, list);
    } else if (op == GenSym.b_add) {
        err = user_list_append(u, list);
    } else {
        /* must be B_SUB */
        err = user_list_subtract(u, list, p.dset);
    }

    return err;
}

int node_replace_scalar(Node * n, double x)
{
    int err = 0;

    if (n.uv != null && n.uv.type == GretlType.double_) {
        uvar_set_scalar_fast(n.uv, x);
    } else {
        if (n.uv == null) {
            fprintf(stderr, "*** node_replace scalar: node uv is NULL!\n");
        } else {
            fprintf(stderr, "*** node_replace scalar: node uv of wrong type!\n");
        }
        err = gretl_scalar_set_value(n.vname, x);
    }

    return err;
}

int gen_replace_scalar(Parser * p, double x)
{
    int err = 0;

    if (p.lh.uv != null && p.lh.uv.type == GretlType.double_) {
        uvar_set_scalar_fast(p.lh.uv, x);
    } else {
        if (p.lh.uv == null) {
            fprintf(stderr, "*** gen_replace scalar: LHS uv is NULL!\n");
        } else {
            fprintf(stderr, "***gen_ replace scalar: LHS uv of wrong type!\n");
        }
        err = gretl_scalar_set_value(p.lh.name.ptr, x);
    }

    return err;
}

/* end of functions that can probably be slimmed down */

void eval_warning(Parser * p, int op, int errnum)
{
    if (!check_gretl_warning()) {
        const(char) * w = (op == GenSym.b_pow)? "pow".ptr : getsymb(op);
        const(char) * s = (errnum)? gretl_strerror(errnum) : null;

        if (s != null) {
            gretl_warnmsg_sprintf("%s: %s", w, s);
        } else {
            gretl_warnmsg_set(w);
        }
    }
}

/* evaluation of binary operators (yielding x op y) for scalar operands
   (also increment/decrement operators)
*/

double xy_calc(double x, double y, int op, int targ, Parser * p)
{
    double z = NADBL;

    /* assignment */
    if (op == GenSym.b_asn || op == GenSym.b_dotasn) {
        return y;
    }

    /* testing for presence of NAs? */
    if ((p.flags & GenFlags.natest) && (na(x) || na(y))) {
        return NADBL;
    }

    /* 0 times anything (even NA) = 0 ? But let's not do this
       for matrices */
    if (targ != GenSym.mat && op == GenSym.b_mul && (x == 0 || y == 0)) {
        return 0;
    }

    /* logical OR: if x or y is valid and non-zero, ignore NA for
       the other term */
    if (op == GenSym.b_or && ((!na(x) && x != 0) || (!na(y) && y != 0))) {
        return 1.0;
    }

    /* logical AND: if either x or y is false, the logical product
       should be false, even if the other term is NA */
    if (op == GenSym.b_and && (x == 0 || y == 0)) {
        return 0;
    }

    if (na(x) || na(y)) {
        /* NaN always propagates to the result */
        return NADBL;
    }

    errno = 0;

    switch (op) {
    case GenSym.b_add:
    case GenSym.inc:
        return x + y;
    case GenSym.b_sub:
    case GenSym.dec:
        return x - y;
    case GenSym.b_mul:
        return x * y;
    case GenSym.b_div:
        return x / y;
    case GenSym.b_mod:
        return fmod(x, y);
    case GenSym.b_and:
        return x != 0 && y != 0;
    case GenSym.b_or:
        return x != 0 || y != 0;
    case GenSym.b_eq:
        return x == y;
    case GenSym.b_neq:
        return x != y;
    case GenSym.b_gt:
        return x > y;
    case GenSym.b_lt:
        return x < y;
    case GenSym.b_gte:
        return x >= y;
    case GenSym.b_lte:
        return x <= y;
    case GenSym.b_pow:
        z = pow(x, y);
        if (errno) {
            eval_warning(p, op, errno);
        }
        return z;
    default:
        return z;
    }
}

int rmatrix_xy_calc(GretlMatrixRaw * targ,
                    GretlMatrixRaw * src,
                    double x, int xleft,
                    int op, Parser * p)
{
    int i, n = targ.rows * targ.cols;

    if (xleft) {
        for (i=0; i<n; i++) {
            targ.val[i] = xy_calc(x, src.val[i], op, GenSym.mat, p);
        }
    } else {
        for (i=0; i<n; i++) {
            targ.val[i] = xy_calc(src.val[i], x, op, GenSym.mat, p);
        }
    }

    return p.err;
}

int operator_real_only(int op)
{
    gretl_errmsg_sprintf("'%s': %s", getsymb(op),
                         _("complex operands are not supported").ptr);
    return GretlError.cmplx;
}

int function_real_only(int f)
{
    gretl_errmsg_sprintf("%s: %s", getsymb(f),
                         _("complex arguments are not supported").ptr);
    return GretlError.cmplx;
}

Cmplx c_xy_calc(Cmplx x, Cmplx y, int op, Parser * p)
{
    if (op == GenSym.b_asn || op == GenSym.b_dotasn) {
        return y;
    }

    switch (op) {
    case GenSym.b_add:
        return x + y;
    case GenSym.b_sub:
        return x - y;
    case GenSym.b_mul:
        return x * y;
    case GenSym.b_div:
        return x / y;
    case GenSym.b_eq:
        return Cmplx(x == y, 0);
    case GenSym.b_neq:
        return Cmplx(x != y, 0);
    default:
        p.err = operator_real_only(op);
        return Cmplx(NADBL, 0);
    }
}

int cmatrix_xy_calc(GretlMatrixRaw * targ,
                    GretlMatrixRaw * src,
                    Cmplx x,
                    int xleft, int op,
                    Parser * p)
{
    int i, n = targ.rows * targ.cols;

    if (xleft) {
        for (i=0; i<n && !p.err; i++) {
            targ.z[i] = c_xy_calc(x, src.z[i], op, p);
        }
    } else {
        for (i=0; i<n && !p.err; i++) {
            targ.z[i] = c_xy_calc(src.z[i], x, op, p);
        }
    }

    return p.err;
}

double cmatrix_xy_comp(GretlMatrixRaw * m, double x,
                       int op, Parser * p)
{
    int i, n = m.rows * m.cols;
    Cmplx zcond;
    Cmplx z = Cmplx(x, 0);
    double ret = 1;

    for (i=0; i<n && ret==1; i++) {
        zcond = c_xy_calc(m.z[i], z, op, p);
        if (p.err) {
            ret = NADBL;
        } else if (zcond.r == 0 && zcond.i == 0) {
            ret = 0;
        }
    }

    return ret;
}

void n_args_error(int k, int nmin, int nmax, int f, Parser * p)
{
    if (nmax == nmin) {
        gretl_errmsg_sprintf(_("Number of arguments (%d) does not "
                               ~ "match the number of\nparameters for "
                               ~ "function %s (%d)").ptr, k, getsymb(f), nmin);
    } else {
        gretl_errmsg_sprintf(_("Number of arguments (%d) does not "
                               ~ "match the number of\nparameters for "
                               ~ "function %s (%d to %d)").ptr, k, getsymb(f),
                             nmin, nmax);
    }
    p.err = 1;
}

bool randgen(int f)
{
    return (f == GenFunc.f_randgen || f == GenFunc.f_mrandgen || f == GenFunc.f_randgen1);
}

int check_dist_count(int d, int f, int * np, int * argc)
{
    int err = 0;

    *np = *argc = 0;

    if (d == DistCode.nc_t || d == DistCode.nc_f || d == DistCode.nc_chisq) {
        /* non-central t, chisq and F: only CDF, PDF and INVCDF supported */
        if (f == GenFunc.f_pdf || f == GenFunc.f_cdf || f == GenFunc.f_invcdf) {
            *np = (d == DistCode.nc_f) ? 3 : 2;
        } else {
            err = GretlError.invarg;
        }
    } else if (d == DistCode.uniform || d == DistCode.udiscrt) {
        /* only randgen is supported */
        if (randgen(f)) {
            *np = 2; /* min, max */
        } else {
            err = GretlError.invarg;
        }
    } else if (d == DistCode.normal) {
        /* all functions supported */
        if (randgen(f)) {
            *np = 2; /* mu, sigma */
        } else {
            *np = 0; /* N(0,1) is assumed */
        }
    } else if (d == DistCode.student) {
        /* Student t: all functions supported */
        *np = 1; /* df */
    } else if (d == DistCode.chisq) {
        /* chi-square: all functions supported */
        *np = 1; /* df */
    } else if (d == DistCode.snedecor) {
        /* all functions supported */
        *np = 2; /* dfn, dfd */
    } else if (d == DistCode.gamma) {
        /* partial support */
        if (f == GenFunc.f_crit) {
            err = 1;
        } else {
            *np = 2; /* shape, scale */
        }
    } else if (d == DistCode.binomial) {
        *np = 2; /* prob, trials */
    } else if (d == DistCode.binorm) {
        /* bivariate normal: cdf only */
        if (f == GenFunc.f_cdf) {
            *np = 1; /* rho */
            *argc = 2; /* note: special */
        } else {
            err = GretlError.invarg;
        }
    } else if (d == DistCode.poisson) {
        *np = 1;
    } else if (d == DistCode.expon) {
        /* inverse cdf not supported */
        if (f == GenFunc.f_invcdf) {
            err = GretlError.invarg;
        } else {
            *np = 1; /* scale */
        }
    } else if (d == DistCode.weibull) {
        /* inverse cdf not supported */
        if (f == GenFunc.f_invcdf) {
            err = GretlError.invarg;
        } else {
            *np = 2; /* shape, scale */
        }
    } else if (d == DistCode.ged) {
        *np = 1; /* shape */
    } else if (d == DistCode.laplace) {
        *np = 2; /* mean, scale */
    } else if (d == DistCode.dw) {
        /* Durbin-Watson: only critical value */
        if (f == GenFunc.f_crit) {
            *np = 2; /* n, k */
        } else {
            err = GretlError.invarg;
        }
    } else if (d == DistCode.johansen) {
        /* Johansen trace test: only p-value */
        if (f == GenFunc.f_pval) {
            *np = 3;
        } else {
            err = GretlError.invarg;
        }
    } else if (d == DistCode.beta) {
        /* cdf, pdf, randgen only */
        if (f == GenFunc.f_cdf || f == GenFunc.f_pdf || randgen(f)) {
            *np = 2; /* shape1, shape2 */
        } else {
            err = GretlError.invarg;
        }
    } else if (d == DistCode.betabin) {
        /* randgen only */
        if (randgen(f)) {
            *np = 3; /* n, shape1, shape2 */
        } else {
            err = GretlError.invarg;
        }
    } else if (d == DistCode.logistic) {
        if (randgen(f)) {
            *np = 2; /* location, scale */
        } else if (f == GenFunc.f_cdf) {
            *np = 0; /* (0,1) assumed */
        } else {
            err = GretlError.invarg;
        }
    } else if (d == DistCode.dirichlet) {
        /* mrandgen only */
        if (randgen(f)) {
            *np = 1; /* alpha vector */
        } else {
            err = GretlError.invarg;
        }
    } else if (d == DistCode.discrete) {
        /* randgen and mrandgen only */
        if (randgen(f)) {
            *np = 1; /* prob vector */
        } else {
            err = GretlError.invarg;
        }
    } else if (d == DistCode.tnormal) {
        /* cdf, pdf, randgen only */
        if (f == GenFunc.f_cdf || f == GenFunc.f_pdf || randgen(f)) {
            *np = 2; /* truncation limits */
        } else {
            err = GretlError.invarg;
        }
    } else {
        err = GretlError.invarg;
    }

    if (!err && !randgen(f) && *argc == 0) {
        *argc = 1;
    }

    return err;
}

double scalar_pdist(int t, int d, const(double) * parm,
                    int np, double arg, Parser * p)
{
    double x = NADBL;
    int i;

    for (i=0; i<np; i++) {
        if (na(parm[i])) {
            return NADBL;
        }
    }

    if (t == GenFunc.f_pval) {
        x = gretl_get_pvalue(d, parm, arg);
    } else if (t == GenFunc.f_pdf) {
        x = gretl_get_pdf(d, parm, arg);
    } else if (t == GenFunc.f_cdf) {
        x = gretl_get_cdf(d, parm, arg);
    } else if (t == GenFunc.f_invcdf) {
        x = gretl_get_cdf_inverse(d, parm, arg);
    } else if (t == GenFunc.f_crit) {
        x = gretl_get_critval(d, parm, arg);
    } else {
        p.err = GretlError.parse;
    }

    return x;
}

/* @parm contains an array of scalar parameters;
   @argvec contains a series of argument values.
*/

int series_pdist(double * x, int f, int d,
                 double * parm, int np,
                 const(double) * argvec,
                 Parser * p)
{
    int t;

    if (f == GenFunc.f_pdf) {
        /* fast treatment, for pdf only at this point */
        int n = sample_size(p.dset);

        for (t=p.dset.t1; t<=p.dset.t2; t++) {
            x[t] = argvec[t];
        }
        gretl_fill_pdf_array(d, parm, x + p.dset.t1, n);
    } else {
        for (t=p.dset.t1; t<=p.dset.t2; t++) {
            x[t] = scalar_pdist(f, d, parm, np, argvec[t], p);
        }
    }

    return 0;
}

/* @parm contains an array of zero to two scalar parameters;
   @argmat contains an array of argument values.
*/

GretlMatrixRaw * matrix_pdist(int f, int d,
                           double * parm, int np,
                           GretlMatrixRaw * argmat,
                           Parser * p)
{
    GretlMatrixRaw * m;
    double x;
    int i, n;

    if (gretl_is_null_matrix(argmat)) {
        return gretl_null_matrix_new();
    }

    m = gretl_matrix_alloc(argmat.rows, argmat.cols);
    if (m == null) {
        p.err = GretlError.alloc;
        return null;
    }

    n = m.rows * m.cols;

    for (i=0; i<n && !p.err; i++) {
        x = scalar_pdist(f, d, parm, np, argmat.val[i], p);
        if (na(x)) {
            p.err = GretlError.missdata;
        } else {
            m.val[i] = x;
        }
    }

    if (p.err) {
        gretl_matrix_free(m);
        m = null;
    }

    return m;
}

/* Gets a matrix from a node of type MAT or NUM. In the
   latter case it's a static matrix, good for use in
   calculation, but it should NOT be passed on or
   freed.
*/

enum N_STATIC = 8;

/* the function-level statics of node_get_matrix() */
__gshared GretlMatrixRaw[N_STATIC] node_get_matrix_mm;
__gshared double[N_STATIC] node_get_matrix_mval = 0;
__gshared int node_get_matrix_started;

GretlMatrixRaw * node_get_matrix(Node * n, Parser * p,
                              int i, int argnum)
{
    if (p.err) {
        /* don't compound prior error */
        return null;
    } else if (n.t == GenSym.mat) {
        return n.v.m;
    } else if (n.t != GenSym.num) {
        if (argnum > 0) {
            gretl_errmsg_sprintf(_("arg %d is missing or of invalid type").ptr,
                                 argnum);
        }
        p.err = GretlError.invarg;
        return null;
    } else if (i < 0 || i >= N_STATIC) {
        fprintf(stderr, "node_get_matrix: invalid static matrix index\n");
        p.err = GretlError.data;
        return null;
    } else {
        /* "fake" a 1 x 1 matrix */
        GretlMatrixRaw * ret;

        if (!node_get_matrix_started) {
            int j;

            for (j=0; j<N_STATIC; j++) {
                gretl_matrix_init(&node_get_matrix_mm[j]);
                node_get_matrix_mm[j].rows = node_get_matrix_mm[j].cols = 1;
            }
            node_get_matrix_started = 1;
        }
        node_get_matrix_mval[i] = n.v.xval;
        ret = &node_get_matrix_mm[i];
        ret.val = &node_get_matrix_mval[i];

        return ret;
    }
}

GretlMatrixRaw * node_get_real_matrix(Node * n, Parser * p,
                                   int i, int argnum)
{
    GretlMatrixRaw * m = node_get_matrix(n, p, i, argnum);

    if (!p.err && m.is_complex) {
        p.err = GretlError.cmplx;
        return null;
    } else {
        return m;
    }
}

double node_get_scalar(Node * n, Parser * p)
{
    if (n.t == GenSym.num) {
        return n.v.xval;
    } else if (scalar_matrix_node(n)) {
        return n.v.m.val[0];
    } else if (n.t == GenSym.empty) {
        p.err = GretlError.args;
        return NADBL;
    } else {
        p.err = GretlError.invarg;
        return NADBL;
    }
}

int node_get_int(Node * n, Parser * p)
{
    double x = node_get_scalar(n, p);

    if (p.err) {
        return -1;
    } else {
        return gretl_int_from_double(x, &p.err);
    }
}

guint32 node_get_guint32(Node * n, Parser * p)
{
    double x = node_get_scalar(n, p);

    if (p.err) {
        return 0;
    } else {
        return gretl_uint32_from_double(x, &p.err);
    }
}

int node_get_bool(Node * n, Parser * p, int deflt)
{
    int ret = -1;

    if (!null_node(n)) {
        int k = node_get_int(n, p);

        if (!p.err) {
            ret = (k != 0);
        }
    } else if (deflt == 0 || deflt == 1) {
        ret = deflt;
    } else {
        p.err = GretlError.args;
    }

    return ret;
}

Node * dec2bin_node(Node * n, Parser * p)
{
    Node * ret = aux_matrix_node(p);

    if (!p.err) {
        GretlMatrixRaw * v = null;
        double x = NADBL;

        if (n.t == GenSym.num) {
            x = n.v.xval;
        } else {
            v = n.v.m;
        }
        if (!p.err) {
            ret.v.m = dec2bin(x, v, &p.err);
        }
    }

    return ret;
}

Node * bin2dec_node(Node * n, Parser * p)
{
    Node * ret = aux_matrix_node(p);

    if (!p.err) {
        ret.v.m = bin2dec(n.v.m, &p.err);
    }

    return ret;
}

Node * access_series_node(Node * n, Parser * p)
{
    Node * ret = aux_list_node(p);

    if (!p.err) {
        int v0 = caller_series_index(p.dset, n.v.str);

        if (v0 > 0) {
            int t, v1 = p.dset.v;

            p.err = dataset_add_series(p.dset, 1);
            if (!p.err) {
                for (t=0; t<p.dset.n; t++) {
                    p.dset.Z[v1][t] = p.dset.Z[v0][t];
                }
                strcpy(p.dset.varname[v1], n.v.str);
                ret.v.ivec = gretl_list_new(1);
                ret.v.ivec[1] = v1;
            }
        } else {
            ret.v.ivec = gretl_list_new(0);
        }
    }

    return ret;
}

Node * DW_node(Node * r, Parser * p)
{
    Node * ret = null;
    int i;
    int[2] parm = 0;

    for (i=0; i<2 && !p.err; i++) {
        parm[i] = node_get_int(r.v.bn.n[i+1], p);
    }

    if (!p.err && (parm[0] < 6 || parm[1] < 0)) {
        p.err = GretlError.invarg;
    }

    if (!p.err) {
        ret = aux_matrix_node(p);
        if (ret != null) {
            ret.v.m = gretl_get_DW(parm[0], parm[1], &p.err);
        }
    }

    return ret;
}

Node * eval_urcpval(Node * n, Node * r, Parser * p)
{
    Node * ret = null;

    if (starting(p)) {
        Node * e;
        int i, m = r.v.bn.n_nodes;
        int[3] iargs = 0;
        double tau = NADBL;

        if (m != 4) {
            n_args_error(m, 4, 4, n.t, p);
        }

        /* need double, int, int, int */
        for (i=0; i<4 && !p.err; i++) {
            e = r.v.bn.n[i];
            if (scalar_node(e)) {
                if (i == 0) {
                    tau = node_get_scalar(e, p);
                } else {
                    iargs[i-1] = node_get_int(e, p);
                }
            } else {
                p.err = e_types(e);
            }
        }

        if (!p.err) {
            int nobs = iargs[0];
            int niv = iargs[1];
            int itv = iargs[2];

            ret = aux_scalar_node(p);
            if (ret != null) {
                ret.v.xval = get_urc_pvalue(tau, nobs, niv, itv);
            }
        }
    } else {
        ret = aux_any_node(p);
    }

    return ret;
}

Node * binperm_node(Node * l, Node * r, Parser * p)
{
    Node * ret = aux_matrix_node(p);

    if (ret != null && starting(p)) {
        int n = node_get_int(l, p);
        int k = node_get_int(r, p);

        if (!p.err) {
            ret.v.m = bit_permutations(n, k, &p.err);
        }
    }

    return ret;
}

int get_matrix_size(GretlMatrixRaw * a, GretlMatrixRaw * b,
                    int * r, int * c)
{
    int err = 0;

    /* if both matrices are present, they must be the
       same size */

    if (a != null) {
        *r = a.rows;
        *c = a.cols;
        if (b != null && (b.rows != *r || b.cols != *c)) {
            err = GretlError.nonconf;
        }
    } else if (b != null) {
        *r = b.rows;
        *c = b.cols;
    } else {
        *r = *c = 0;
    }

    return err;
}

Node * bvnorm_node(Node * n, Parser * p)
{
    Node * ret = null;

    if (starting(p)) {
        double * avec = null;
        double * bvec = null;
        GretlMatrixRaw * amat = null;
        GretlMatrixRaw * bmat = null;
        double a, b;
        double[2] args;
        double rho = NADBL;
        Node * e;
        int i, mode = 0;

        for (i=0; i<3 && !p.err; i++) {
            e = n.v.bn.n[i+1];
            if (scalar_node(e)) {
                if (i == 0) {
                    rho = node_get_scalar(e, p);
                } else {
                    args[i-1] = node_get_scalar(e, p);
                }
            } else if (i == 1) {
                if (e.t == GenSym.series) {
                    avec = e.v.xvec;
                } else if (e.t == GenSym.mat) {
                    amat = e.v.m;
                }
            } else if (i == 2) {
                if (e.t == GenSym.series) {
                    bvec = e.v.xvec;
                } else if (e.t == GenSym.mat) {
                    bmat = e.v.m;
                }
            } else {
                node_type_error(GenFunc.f_cdf, i+1, GenSym.num, e, p);
            }
        }

        if (!p.err) {
            if ((avec != null && bmat != null) ||
                (bvec != null && amat != null)) {
                p.err = GretlError.invarg;
            } else if (avec != null || bvec != null) {
                mode = 1;
                ret = aux_series_node(p);
            } else if (amat != null || bmat != null) {
                mode = 2;
                ret = aux_matrix_node(p);
            } else {
                mode = 0;
                ret = aux_scalar_node(p);
            }
        }

        if (p.err) {
            return ret;
        }

        if (mode == 0) {
            /* a, b are both scalars */
            ret.v.xval = bvnorm_cdf(rho, args[0], args[1]);
        } else if (mode == 1) {
            /* a and/or b are series */
            int t;

            for (t=p.dset.t1; t<=p.dset.t2; t++) {
                a = (avec != null)? avec[t] : args[0];
                b = (bvec != null)? bvec[t] : args[1];
                if (na(a) || na(b)) {
                    ret.v.xvec[t] = NADBL;
                } else {
                    ret.v.xvec[t] = bvnorm_cdf(rho, a, b);
                }
            }
        } else if (mode == 2) {
            /* a and/or b are matrices */
            GretlMatrixRaw * m = null;
            int r, c;

            p.err = get_matrix_size(amat, bmat, &r, &c);

            if (!p.err && r > 0 && c > 0) {
                m = gretl_matrix_alloc(r, c);
            }
            if (m != null) {
                int ii, nn = r * c;

                for (ii=0; ii<nn && !p.err; ii++) {
                    a = (amat != null)? amat.val[ii] : args[0];
                    b = (bmat != null)? bmat.val[ii] : args[1];
                    m.val[ii] = bvnorm_cdf(rho, a, b);
                    if (na(m.val[ii])) {
                        /* matrix: change NAs to NaNs */
                        m.val[ii] = double.nan;
                    }
                }
            }
            if (ret.v.m != null) {
                gretl_matrix_free(ret.v.m);
            }
            ret.v.m = m;
        }
    } else {
        ret = aux_any_node(p);
    }

    return ret;
}

Node * dirichlet_node(Node * args, Parser * p)
{
    Node * ret = null;
    GretlMatrixRaw * alpha = null;
    int n = 0;

    if (args.v.bn.n[1].t == GenSym.mat) {
        alpha = args.v.bn.n[1].v.m;
    }
    if (args.v.bn.n[2].t == GenSym.num) {
        n = node_get_int(args.v.bn.n[2], p);
    }
    if (alpha == null || n <= 0) {
        p.err = GretlError.invarg;
    } else {
        ret = aux_matrix_node(p);
        if (!p.err) {
            ret.v.m = gretl_rand_dirichlet(alpha, n, &p.err);
        }
    }

    return ret;
}

Node * discrete_rand_node(Node * n, Node * args, Parser * p)
{
    Node * ret = null;
    GretlMatrixRaw * probs = null;
    int vsize = 0;

    if (args.v.bn.n[1].t == GenSym.mat) {
        probs = args.v.bn.n[1].v.m;
        if (gretl_vector_get_length(probs) == 0) {
            p.err = GretlError.invarg;
            return null;
        }
    }

    if (n.t == GenFunc.f_randgen1) {
        ret = aux_scalar_node(p);
    } else if (n.t == GenFunc.f_randgen) {
        ret = aux_series_node(p);
    } else {
        /* F_MRANDGEN */
        int r = node_get_int(args.v.bn.n[2], p);
        int c = node_get_int(args.v.bn.n[3], p);

        if (!p.err && (r <= 0 || c <= 0)) {
            p.err = GretlError.invarg;
        } else {
            ret = aux_sized_matrix_node(p, r, c, 0);
            vsize = r * c;
        }
    }

    if (!p.err && n.t == GenFunc.f_randgen1) {
        p.err = gretl_rand_discrete(&ret.v.xval, 0, 0, probs);
    } else if (!p.err && n.t == GenFunc.f_randgen) {
        p.err = gretl_rand_discrete(ret.v.xvec,
                                    p.dset.t1, p.dset.t2,
                                    probs);
    } else if (!p.err) {
        /* F_MRANDGEN */
        p.err = gretl_rand_discrete(ret.v.m.val, 0, vsize-1, probs);
    }

    return ret;
}

/* Return a node containing the evaluated result of a probability
   distribution function. Here, node @n is the function node and node
   @r holds the arguments, the number of which ir not known in
   advance.
*/

Node * eval_pdist(Node * n, Node * r, Parser * p)
{
    Node * ret = null;

    if (starting(p)) {
        Node * e;
        Node * s;
        int i, k, nr = r.v.bn.n_nodes;
        int rgen = (n.t == GenFunc.f_randgen);
        int mrgen = (n.t == GenFunc.f_mrandgen);
        int rgen1 = (n.t == GenFunc.f_randgen1);
        double[3] parm = 0;
        double argval = NADBL;
        double *[2] parmvec = [null, null];
        double * argvec = null;
        int[2] pvlen = 0;
        GretlMatrixRaw * argmat = null;
        int rows = 0, cols = 0;
        int np, argc, bb, d = 0;

        if (nr < 2) {
            p.err = GretlError.args;
            goto disterr;
        }

        /* the first argument must hold a string identifying
           the distribution to be employed
        */
        s = r.v.bn.n[0];
        if (s.t == GenSym.str) {
            char * dstr = s.v.str;

            d = dist_code_from_string(dstr);
            if (d == 0) {
                dstr = get_string_by_name(dstr);
                if (dstr != null) {
                    d = dist_code_from_string(dstr);
                }
            }
            if (d == 0) {
                p.err = GretlError.invarg;
            }
        } else {
            node_type_error(n.t, 0, GenSym.str, s, p);
        }

        if (mrgen) {
            if (d == DistCode.dirichlet) {
                if (nr != 3) {
                    n_args_error(nr, 3, 3, n.t, p);
                }
            } else if (d == DistCode.discrete) {
                if (nr != 4) {
                    n_args_error(nr, 4, 4, n.t, p);
                }
            } else if (nr < 4 || nr > 7) {
                n_args_error(nr, 4, 7, n.t, p);
            }
        } else if (nr > 5) {
            n_args_error(nr, 2, 5, n.t, p);
        } else {
            if (d == DistCode.dirichlet) {
                /* mrandgen only */
                p.err = GretlError.invarg;
            }
        }

        if (p.err) {
            goto disterr;
        } else {
            p.err = check_dist_count(d, n.t, &np, &argc);
        }
        if (!p.err && d != DistCode.dirichlet) {
            k = np + argc + 2 * mrgen;
            if (k != nr - 1) {
                p.err = GretlError.invarg;
            }
        }
        if (p.err) {
            goto disterr;
        }

        bb = (d == DistCode.betabin);

        if (d == DistCode.dw) {
            /* special: Durbin-Watson */
            return DW_node(r, p);
        } else if (d == DistCode.binorm) {
            /* special: bivariate normal */
            return bvnorm_node(r, p);
        } else if (d == DistCode.dirichlet) {
            /* special: dirichlet */
            return dirichlet_node(r, p);
        } else if (d == DistCode.discrete) {
            return discrete_rand_node(n, r, p);
        }

        for (i=1; i<=k && !p.err; i++) {
            e = r.v.bn.n[i];
            if (scalar_node(e)) {
                /* scalars always acceptable */
                if (mrgen) {
                    if (i == k) {
                        cols = node_get_int(e, p);
                    } else if (i == k-1) {
                        rows = node_get_int(e, p);
                    } else {
                        parm[i-1] = node_get_scalar(e, p);
                    }
                } else if (i == k && argc > 0) {
                    argval = node_get_scalar(e, p);
                } else {
                    parm[i-1] = node_get_scalar(e, p);
                }
            } else if (i == k && e.t == GenSym.series) {
                /* a series in the last place? */
                if (bb) {
                    node_type_error(n.t, i, GenSym.num, e, p);
                } else if (rgen) {
                    parmvec[i-1] = e.v.xvec;
                } else if (mrgen) {
                    node_type_error(n.t, i, GenSym.num, e, p);
                } else {
                    argvec = e.v.xvec;
                }
            } else if (i == k && e.t == GenSym.mat) {
                /* a matrix in the last place? */
                if (mrgen) {
                    parmvec[i-1] = e.v.m.val;
                    pvlen[i-1] = e.v.m.rows * e.v.m.cols;
                } else if (rgen) {
                    node_type_error(n.t, i, GenSym.num, e, p);
                } else {
                    argmat = e.v.m;
                }
            } else if (e.t == GenSym.series) {
                /* a series param for randgen? */
                if (rgen && !bb) {
                    parmvec[i-1] = e.v.xvec;
                } else {
                    node_type_error(n.t, i, GenSym.num, e, p);
                }
            } else if (e.t == GenSym.mat) {
                /* a matrix param for mrandgen? */
                if (mrgen && !bb) {
                    parmvec[i-1] = e.v.m.val;
                    pvlen[i-1] = e.v.m.rows * e.v.m.cols;
                } else {
                    node_type_error(n.t, i, GenSym.num, e, p);
                }
            } else {
                p.err = GretlError.invarg;
                fprintf(stderr, "eval_pdist: arg %d, bad type %d\n", i+1, e.t);
            }
        }

        if (mrgen) {
            int rlen = rows * cols;

            if ((parmvec[0] != null && pvlen[0] != rlen) ||
                (parmvec[1] != null && pvlen[1] != rlen)) {
                p.err = GretlError.nonconf;
            }
        }

        if (p.err) {
            goto disterr;
        }

        if (mrgen) {
            ret = aux_matrix_node(p);
        } else if (rgen || argvec != null) {
            ret = aux_series_node(p);
        } else if (argmat != null) {
            ret = aux_matrix_node(p);
        } else {
            ret = aux_scalar_node(p);
        }

        if (ret == null) {
            goto disterr;
        }

        if (rgen) {
            p.err = gretl_fill_random_series(ret.v.xvec, d, parm.ptr,
                                             parmvec[0], parmvec[1],
                                             p.dset);
        } else if (mrgen) {
            ret.v.m = gretl_get_random_matrix(d, parm.ptr,
                                              parmvec[0], parmvec[1],
                                              rows, cols,
                                              &p.err);
        } else if (rgen1) {
            ret.v.xval = gretl_get_random_scalar(d, parm.ptr, &p.err);
        } else if (argvec != null) {
            p.err = series_pdist(ret.v.xvec, n.t, d, parm.ptr, np,
                                 argvec, p);
        } else if (argmat != null) {
            ret.v.m = matrix_pdist(n.t, d, parm.ptr, np, argmat, p);
        } else {
            ret.v.xval = scalar_pdist(n.t, d, parm.ptr, np, argval, p);
        }
    } else {
        ret = aux_any_node(p);
    }

  disterr:

    return ret;
}

__gshared int mpi_rank = -1;
__gshared int mpi_size = 0;

void set_mpi_rank_and_size(int rank, int size)
{
    mpi_rank = rank;
    mpi_size = size;
}

double get_const_by_id(int id)
{
    if (id == ConstId.pi) {
        return M_PI;
    } else if (id == ConstId.eps) {
        /* see https://en.wikipedia.org/wiki/Machine_epsilon :
           we now use the (b) definition, as per Matlab, Gauss,
           R and others
        */
        return pow(2.0, -52);
    } else if (id == ConstId.inf) {
        return double.infinity;
    } else if (id == ConstId.nan) {
        return double.nan;
    } else if (id == ConstId.win32) {
        return 0;
    } else if (id == ConstId.have_mpi) {
        return 0;
    } else if (id == ConstId.mpi_rank) {
        return mpi_rank;
    } else if (id == ConstId.mpi_size) {
        return mpi_size;
    } else if (id == ConstId.n_proc) {
        return gretl_n_processors();
    } else if (id == ConstId.true_) {
        return 1;
    } else if (id == ConstId.false_) {
        return 0;
    } else {
        return NADBL;
    }
}

/* look up and return numerical values of symbolic constants */

Node * retrieve_const(Node * n, Parser * p)
{
    Node * ret = aux_scalar_node(p);

    if (ret != null && starting(p)) {
        ret.v.xval = get_const_by_id(n.v.idnum);
    }

    return ret;
}

double get_const_by_name(const(char) * name, int * err)
{
    int id = const_lookup(name);

    if (id > 0) {
        return get_const_by_id(id);
    } else {
        if (err != null) {
            *err = GretlError.data;
        }
        return NADBL;
    }
}

/* HAVE_MPI is not defined: the non-MPI versions */

Node * mpi_transfer_node(Node * l, Node * r, Node * r2,
                         int f, Parser * p)
{
    gretl_errmsg_set(_("MPI is not supported in this gretl build").ptr);
    p.err = 1;
    return null;
}

Node * mpi_barrier_node(Parser * p)
{
    gretl_errmsg_set(_("MPI is not supported in this gretl build").ptr);
    p.err = 1;
    return null;
}

Node * scalar_calc(Node * x, Node * y, int f, Parser * p)
{
    Node * ret = aux_scalar_node(p);

    if (ret != null && starting(p)) {
        ret.v.xval = xy_calc(x.v.xval, y.v.xval, f, GenSym.num, p);
    }

    return ret;
}

/* support addition of integer meaning subset of string */

Node * string_offset(Node * l, Node * r, Parser * p)
{
    Node * ret = aux_string_node(p);

    if (ret != null && starting(p)) {
        int n = cast(int) g_utf8_strlen(l.v.str, -1);
        int k = cast(int) r.v.xval;

        if (k < 0) {
            p.err = GretlError.data;
        } else if (k >= n) {
            ret.v.str = gretl_strdup("");
        } else {
            char * ps = g_utf8_offset_to_pointer(l.v.str, k);

            ret.v.str = gretl_strdup(ps);
        }
        if (!p.err && ret.v.str == null) {
            p.err = GretlError.alloc;
        }
    }

    return ret;
}

/* support multiplication by integer meaning repetition of a string */

Node * string_mult(Node * l, Node * r, Parser * p)
{
    Node * ret = aux_string_node(p);

    if (ret != null && starting(p)) {
        const(char) * s = null;
        int k = 0;

        if (l.t == GenSym.str && r.t == GenSym.num) {
            s = l.v.str;
            k = node_get_int(r, p);
        } else if (l.t == GenSym.num && r.t == GenSym.str) {
            s = r.v.str;
            k = node_get_int(l, p);
        } else {
            p.err = GretlError.types;
        }
        if (!p.err && k < 0) {
            p.err = GretlError.invarg;
        }
        if (!p.err && k == 0) {
            ret.v.str = gretl_strdup("");
        } else if (!p.err) {
            size_t n = strlen(s);
            long sz64 = k * n + 1;
            int i;

            if (sz64 > cast(long) (INT_MAX / 2)) {
                /* avoid potential wrap-around and overflow */
                p.err = GretlError.invarg;
            } else {
                ret.v.str = cast(char *) calloc(k * n + 1, 1);
                if (ret.v.str == null) {
                    p.err = GretlError.alloc;
                } else {
                    for (i=0; i<k; i++) {
                        strcat(ret.v.str, s);
                    }
                }
            }
        }
    }

    return ret;
}

Node * compare_strings(Node * l, Node * r, int f, Parser * p)
{
    Node * ret = aux_scalar_node(p);

    if (ret != null && starting(p)) {
        int s = strcmp(l.v.str, r.v.str);

        ret.v.xval = (f == GenSym.b_eq)? (s == 0) : (s != 0);
    }

    return ret;
}

/*
   We're looking at a comparison, with either a series on the left and
   a string on the right or vice versa.  This can work if the series
   in question is string-valued, as in

     Case 1: series foo = (x == "strval")

   It can also work if the string is an observation marker, as in

     Case 2: series foo = (obs >= "CA")
*/

Node * series_string_calc(Node * l, Node * r, int f, Parser * p)
{
    double xt = NADBL, yt = NADBL;
    double * x = null;
    double * y = null;
    double * alt;
    const(char) * strval;
    int vnum, t;
    Node * ret;

    if (r.t == GenSym.str) {
        strval = r.v.str;
        vnum = l.vnum;
        x = l.v.xvec;
        alt = &yt;
    } else {
        strval = l.v.str;
        vnum = r.vnum;
        y = r.v.xvec;
        alt = &xt;
    }

    ret = aux_series_node(p);
    if (p.err) {
        return ret;
    }

    if (vnum > 0 && is_string_valued(p.dset, vnum)) {
        /* we must be in Case 1 */
        *alt = series_decode_string(p.dset, vnum, strval);
        if (na(*alt)) {
            /* @strval is not a string value of the given series */
            double xval = (f == GenSym.b_eq)? 0 : (f == GenSym.b_neq)? 1 : NADBL;

            for (t=p.dset.t1; t<=p.dset.t2; t++) {
                ret.v.xvec[t] = xval;
            }
            return ret; /* NA case handled */
        }
    } else {
        /* try interpreting @strval as an observation string */
        if (annual_data(p.dset)) {
            *alt = get_date_x(p.dset.pd, strval);
        } else {
            t = dateton(strval, p.dset);
            if (t >= 0) {
                *alt = t + 1;
            }
        }
    }

    if (na(*alt)) {
        gretl_errmsg_sprintf(_("got invalid field '%s'").ptr, strval);
        p.err = GretlError.types;
        return null;
    }

    if (ret != null) {
        int t1 = autoreg(p) ? p.obs : p.dset.t1;
        int t2 = autoreg(p) ? p.obs : p.dset.t2;

        for (t=t1; t<=t2; t++) {
            if (x != null) {
                xt = x[t];
            } else if (y != null) {
                yt = y[t];
            }
            ret.v.xvec[t] = xy_calc(xt, yt, f, GenSym.series, p);
        }
    }

    return ret;
}

double * list_node_get_series(Node * n, Parser * p)
{
    if (n.v.ivec[0] == 1) {
        int v = n.v.ivec[1];

        if (v >= 0 && v < p.dset.v) {
            return p.dset.Z[v];
        }
    }

    p.err = GretlError.invarg;
    return null;
}

/* At least one of the nodes is a series; the other may be a
   scalar or 1 x 1 matrix */

Node * series_calc(Node * l, Node * r, int f, Parser * p)
{
    Node * ret = aux_series_node(p);
    const(double) * x = null;
    const(double) * y = null;
    double xt = 0, yt = 0;
    int tmax = p.dset.t2;

    if (ret == null) {
        return null;
    }

    if (p.dset.n > p.dset_n) {
        /* can arise when stack() is in the tree ->
           dataset gets extended on the fly
        */
        if (l.t == GenSym.series) {
            if (useries_node(l)) {
                l.v.xvec = p.dset.Z[l.vnum];
            } else {
                tmax = MIN(tmax, p.dset_n - 1);
            }
        }
        if (r.t == GenSym.series) {
            if (useries_node(r)) {
                r.v.xvec = p.dset.Z[r.vnum];
            } else {
                tmax = MIN(tmax, p.dset_n - 1);
            }
        }
    }

    if (l.t == GenSym.series) {
        x = l.v.xvec;
    } else if (l.t == GenSym.list) {
        x = list_node_get_series(l, p);
    } else if (l.t == GenSym.num) {
        xt = l.v.xval;
    } else if (l.t == GenSym.mat) {
        xt = l.v.m.val[0];
    }

    if (r.t == GenSym.series) {
        y = r.v.xvec;
    } else if (r.t == GenSym.list) {
        y = list_node_get_series(r, p);
    } else if (r.t == GenSym.num) {
        yt = r.v.xval;
    } else if (r.t == GenSym.mat) {
        yt = r.v.m.val[0];
    }

    if (!p.err) {
        int t1 = autoreg(p) ? p.obs : p.dset.t1;
        int t2 = autoreg(p) ? p.obs : tmax;
        int t;

        for (t=t1; t<=t2; t++) {
            if (x != null) {
                xt = x[t];
            }
            if (y != null) {
                yt = y[t];
            }
            ret.v.xvec[t] = xy_calc(xt, yt, f, GenSym.series, p);
        }
    }

    return ret;
}

/* Try to determine if it's OK for a function to return a
   string-valued series, given that such a series must be a member of
   a dataset. The node @f should hold the function in question. In
   case the node that's to hand in the caller is actually a top-level
   argument of the function, one can pass its "parent" pointer as @f.
   Returns an error code if this operation is not OK.
*/

int strvals_return_check(Node * f, Parser * p)
{
    if ((f != null && f != p.tree) || (p.targ != GenSym.series && p.targ != GenSym.unk)) {
        /* we must be at the top of the tree, with the target of
           assignment either known to be a series or undetermined
           as yet
        */
        gretl_errmsg_sprintf(_("%s: series usage is valid only in direct assignment").ptr,
                             f == null ? "?".ptr : getsymb(f.t));
        return GretlError.data;
    } else if (dataset_is_subsampled(p.dset) && p.lh.t == GenSym.series &&
               !is_string_valued(p.dset, p.lh.vnum)) {
        /* can't assign string values to just part of a numeric series */
        gretl_errmsg_set(_("Sorry, can't do this with a sub-sampled dataset").ptr);
        return GretlError.data;
    } else {
        /* OK, we'll try it */
        return 0;
    }
}

/* Create and store a series table to be attached to a string-valued
   series return.
*/

void prep_strvals_return(char ** S, int ns, Parser * p)
{
    SeriesTable * st = series_table_new(S, ns, &p.err);

    if (!p.err) {
        p.lh.stab = st;
    } else {
        strings_array_free(S, ns);
    }
}

/* Both nodes are string-valued series. We support a limited
   set of operations.
*/

Node * strvals_calc(Node * l, Node * r, Node * n, Parser * p)
{
    Node * ret = null;
    const(char) * sl;
    const(char) * sr;
    char ** Sx = null;
    int nr = 0, nx = 0;
    int vl, vr, f = n.t;
    int i, t, eq;

    if (f == GenSym.b_pow) {
        p.err = strvals_return_check(n, p);
    } else if (f != GenSym.b_eq && f != GenSym.b_neq) {
        p.err = GretlError.types;
    }
    if (p.err) {
        return null;
    }

    ret = aux_series_node(p);
    if (ret == null) {
        return null;
    }

    vl = l.vnum;
    vr = r.vnum;

    if (f == GenSym.b_pow) {
        /* "logical product" */
        char * slr;
        char ** Sl;
        char ** Sr;
        int nl, j, ll;

        Sl = series_get_string_vals(p.dset, vl, &nl, 1);
        Sr = series_get_string_vals(p.dset, vr, &nr, 1);
        nx = nl * nr;

        Sx = strings_array_new(nx);
        if (Sx == null) {
            p.err = GretlError.alloc;
            return ret;
        }

        for (i=0; i<nl; i++) {
            ll = cast(int) strlen(Sl[i]) + 2;
            for (j=0; j<nr; j++) {
                slr = cast(char *) calloc(ll + strlen(Sr[j]), 1);
                sprintf(slr, "%s.%s", Sl[i], Sr[j]);
                Sx[i*nr+j] = slr;
            }
        }
    }

    for (t=p.dset.t1; t<=p.dset.t2; t++) {
        sl = series_get_string_for_obs(p.dset, vl, t);
        sr = series_get_string_for_obs(p.dset, vr, t);
        if (sl == null || sr == null) {
            ret.v.xvec[t] = NADBL;
        } else if (f == GenSym.b_pow) {
            int il = cast(int) (p.dset.Z[vl][t] - 1);
            int ir = cast(int) (p.dset.Z[vr][t] - 1);

            ret.v.xvec[t] = il*nr + ir + 1;
        } else {
            eq = strcmp(sl, sr) == 0;
            ret.v.xvec[t] = (f == GenSym.b_eq)? eq : !eq;
        }
    }

    if (f == GenSym.b_pow && Sx != null) {
        prep_strvals_return(Sx, nx, p);
    }

    return ret;
}

int op_symbol(int op)
{
    switch (op) {
    case GenSym.b_dotmult: return '*';
    case GenSym.b_dotdiv:  return '/';
    case GenSym.b_dotpow:  return '^';
    case GenSym.b_dotadd:  return '+';
    case GenSym.b_dotsub:  return '-';
    case GenSym.b_doteq:   return '=';
    case GenSym.b_dotgt:   return '>';
    case GenSym.b_dotlt:   return '<';
    case GenSym.b_dotgte:  return ']';
    case GenSym.b_dotlte:  return '[';
    case GenSym.b_dotneq:  return '!';
    default: return 0;
    }
}

GretlMatrixRaw * nullmat_multiply(const(GretlMatrixRaw) * A,
                               const(GretlMatrixRaw) * B,
                               int op, int * err)
{
    GretlMatrixRaw * C = null;

    if (A.rows == 0 && A.cols == 0 &&
        B.rows == 0 && B.cols == 0) {
        C = gretl_null_matrix_new();
    } else {
        int Lc = op == GenSym.b_trmul ? A.rows : A.cols;
        int Cr = op == GenSym.b_trmul ? A.cols : A.rows;
        int Cc = B.cols;

        if (Lc != B.rows) {
            *err = GretlError.nonconf;
        } else {
            if (Cr > 0 && Cc > 0) {
                C = gretl_zero_matrix_new(Cr, Cc);
            } else {
                C = gretl_matrix_alloc(Cr, Cc);
            }
            if (C == null) {
                *err = GretlError.alloc;
            }
        }
    }

    return C;
}

GretlMatrixRaw * matrix_add_sub_scalar(const(GretlMatrixRaw) * A,
                                    const(GretlMatrixRaw) * B,
                                    int op)
{
    GretlMatrixRaw * C;
    double xval;
    const(double) * xvec;
    int r, c;

    if (gretl_matrix_is_scalar(A)) {
        r = B.rows;
        c = B.cols;
        xval = A.val[0];
        xvec = B.val;
    } else {
        r = A.rows;
        c = A.cols;
        xval = B.val[0];
        xvec = A.val;
    }

    C = gretl_matrix_alloc(r, c);

    if (C != null) {
        int i, n = r * c;

        if (op == GenSym.b_add) {
            for (i=0; i<n; i++) {
                C.val[i] = xvec[i] + xval;
            }
        } else {
            if (xvec == A.val) {
                for (i=0; i<n; i++) {
                    C.val[i] = xvec[i] - xval;
                }
            } else {
                for (i=0; i<n; i++) {
                    C.val[i] = xval - xvec[i];
                }
            }
        }
    }

    return C;
}

/* See if we can reuse an existing matrix on an
   auxiliary node. If so, return it; otherwise
   free it and return a newly allocated matrix.
*/

GretlMatrixRaw * calc_get_matrix(GretlMatrixRaw ** pM,
                              int r, int c)
{
    if (*pM == null) {
        /* allocate from scratch */
        return gretl_matrix_alloc(r, c);
    } else if ((*pM).rows == r && (*pM).cols == c) {
        /* reusable as-is */
        return *pM;
    } else if ((*pM).rows == c && (*pM).cols == r) {
        /* reusable if reoriented */
        (*pM).rows = r;
        (*pM).cols = c;
        return *pM;
    } else {
        /* new matrix needed */
        gretl_matrix_free(*pM);
        *pM = null;
        return gretl_matrix_alloc(r, c);
    }
}

bool op_no_complex(int o)
{
    return ((o >= GenSym.b_lt && o <= GenSym.b_gte) ||
            (o >= GenSym.b_dotlt && o <= GenSym.b_dotgte));
}

bool fn_no_complex(int f)
{
    return (f == GenFunc.f_qform || f == GenFunc.f_lsolve ||
            f == GenFunc.f_cmult || f == GenFunc.f_cdiv ||
            f == GenFunc.f_conv2d || f == GenFunc.f_sgn ||
            f == GenFunc.f_max || f == GenFunc.f_min);
}

/* return allocated result of binary operation performed on
   two matrices */

int real_matrix_calc(const(GretlMatrixRaw) * A,
                     const(GretlMatrixRaw) * B,
                     int op, GretlMatrixRaw ** pM)
{
    GretlMatrixMod mod;
    GretlMatrixRaw * C = null;
    int ra, ca;
    int rb, cb;
    int r, c;
    int err = 0;

    if (gretl_is_null_matrix(A) ||
        gretl_is_null_matrix(B)) {
        if (op != GenSym.b_hcat && op != GenSym.b_vcat && op != GenFunc.f_dsum &&
            op != GenSym.b_mul && op != GenSym.b_trmul) {
            return GretlError.nonconf;
        }
        if (op == GenSym.b_mul || op == GenSym.b_trmul) {
            C = nullmat_multiply(A, B, op, &err);
            goto finish;
        }
    }

    if (A.is_complex || B.is_complex) {
        /* gatekeeper for complex */
        if (op_no_complex(op)) {
            return operator_real_only(op);
        } else if (fn_no_complex(op)) {
            return function_real_only(op);
        }
    }

    switch (op) {
    case GenSym.b_add:
    case GenSym.b_sub:
        if (A.is_complex || B.is_complex) {
            int sgn = (op == GenSym.b_sub)? -1 : 1;

            C = gretl_cmatrix_add_sub(A, B, sgn, &err);
        } else if (gretl_matrix_is_scalar(A) ||
                   gretl_matrix_is_scalar(B)) {
            C = matrix_add_sub_scalar(A, B, op);
            if (C == null) {
                err = GretlError.alloc;
            }
        } else {
            C = calc_get_matrix(pM, A.rows, A.cols);
            if (C == null) {
                err = GretlError.alloc;
            } else if (op == GenSym.b_add) {
                err = gretl_matrix_add(A, B, C);
            } else {
                err = gretl_matrix_subtract(A, B, C);
            }
        }
        break;
    case GenSym.b_hcat:
    case GenSym.b_vcat:
        if (op == GenSym.b_hcat) {
            C = gretl_matrix_col_concat(A, B, &err);
        } else {
            C = gretl_matrix_row_concat(A, B, &err);
        }
        break;
    case GenFunc.f_dsum:
        C = gretl_matrix_direct_sum(A, B, &err);
        break;
    case GenSym.b_mul:
        if (A.is_complex || B.is_complex) {
            C = gretl_cmatrix_multiply(A, B, &err);
        } else {
            ra = gretl_matrix_rows(A);
            ca = gretl_matrix_cols(A);
            rb = gretl_matrix_rows(B);
            cb = gretl_matrix_cols(B);
            r = (ra == 1 && ca == 1)? rb : ra;
            c = (rb == 1 && cb == 1)? ca : cb;

            C = calc_get_matrix(pM, r, c);
            if (C == null) {
                err = GretlError.alloc;
            } else {
                err = gretl_matrix_multiply(A, B, C);
                if (!err) {
                    gretl_matrix_transcribe_obs_info(C, A);
                }
            }
        }
        break;
    case GenSym.b_trmul:
        if (A.is_complex || B.is_complex) {
            C = gretl_cmatrix_AHB(A, B, &err);
        } else {
            ra = gretl_matrix_cols(A);
            ca = gretl_matrix_rows(A);
            rb = gretl_matrix_rows(B);
            cb = gretl_matrix_cols(B);

            r = (ra == 1 && ca == 1)? rb : ra;
            c = (rb == 1 && cb == 1)? ca : cb;

            C = calc_get_matrix(pM, r, c);
            if (C == null) {
                err = GretlError.alloc;
            } else {
                err = gretl_matrix_multiply_mod(A, GretlMatrixMod.transpose,
                                                B, GretlMatrixMod.none,
                                                C, GretlMatrixMod.none);
            }
        }
        break;
    case GenFunc.f_qform:
        /* quadratic form, A * B * A', for symmetric B */
        ra = gretl_matrix_rows(A);
        ca = gretl_matrix_cols(A);
        rb = gretl_matrix_rows(B);
        cb = gretl_matrix_cols(B);

        /* input checks */
        if (rb == 1 || cb == 1) {
            /* the diagonal B case */
            if (ca != gretl_vector_get_length(B)) {
                err = GretlError.nonconf;
            }
        } else {
            /* the "full B" case */
            if (ca != rb || cb != rb) {
                err = GretlError.nonconf;
            } else {
                gretl_matrix_set_equals_tolerance(1.0e-7);
                if (!gretl_matrix_is_symmetric(B)) {
                    gretl_errmsg_set(_("Matrix is not symmetric").ptr);
                    err = GretlError.nonconf;
                }
                gretl_matrix_unset_equals_tolerance();
            }
        }
        if (!err) {
            mod = GretlMatrixMod.none;
            C = calc_get_matrix(pM, ra, ra);
            if (C == null) {
                err = GretlError.alloc;
            } else if (rb == 1 || cb == 1) {
                err = gretl_matrix_diag_qform(A, mod, B, C, mod);
            } else {
                err = gretl_matrix_qform(A, mod, B, C, mod);
            }
        }
        break;
    case GenSym.b_ldiv:
    case GenSym.b_div:
        /* Matrix left (A\B) or right (A/B) "division". Note that
           A/B = (B'\A')', which we handle by passing the transpose
           flag to gretl_{c}matrix_divide.
        */
        mod = (op == GenSym.b_ldiv)? GretlMatrixMod.none : GretlMatrixMod.transpose;
        if (A.is_complex || B.is_complex) {
            C = gretl_cmatrix_divide(A, B, mod, &err);
        } else {
            C = gretl_matrix_divide(A, B, mod, &err);
        }
        break;
    case GenFunc.f_lsolve:
        C = calc_get_matrix(pM, B.rows, B.cols);
        if (C == null) {
            err = GretlError.alloc;
        } else {
            gretl_matrix_copy_values(C, B);
            err = gretl_cholesky_solve(A, C);
        }
        break;
    case GenSym.b_dotmult:
    case GenSym.b_dotdiv:
    case GenSym.b_dotpow:
    case GenSym.b_dotadd:
    case GenSym.b_dotsub:
    case GenSym.b_doteq:
    case GenSym.b_dotgt:
    case GenSym.b_dotlt:
    case GenSym.b_dotgte:
    case GenSym.b_dotlte:
    case GenSym.b_dotneq:
        /* apply operator element-wise */
        if (A.is_complex || B.is_complex) {
            C = gretl_cmatrix_dot_op(A, B, cast(char) op_symbol(op), &err);
        } else {
            C = gretl_matrix_dot_op(A, B, cast(char) op_symbol(op), &err);
        }
        break;
    case GenSym.b_kron:
        if (A.is_complex || B.is_complex) {
            C = gretl_cmatrix_kronecker(A, B, &err);
        } else {
            C = gretl_matrix_kronecker_product_new(A, B, &err);
        }
        break;
    case GenFunc.f_hdprod:
        if (A.is_complex || B.is_complex) {
            C = gretl_cmatrix_hdprod(A, B, &err);
        } else {
            C = gretl_matrix_hdproduct_new(A, B, &err);
        }
        break;
    case GenFunc.f_cmult:
        C = gretl_matrix_complex_multiply(A, B, 0, &err);
        break;
    case GenFunc.f_cdiv:
        C = gretl_matrix_complex_divide(A, B, 0, &err);
        break;
    case GenFunc.f_mrsel:
        C = gretl_matrix_bool_sel(A, B, 1, &err);
        break;
    case GenFunc.f_mcsel:
        C = gretl_matrix_bool_sel(A, B, 0, &err);
        break;
    case GenFunc.f_conv2d:
        C = gretl_matrix_2d_convolution(A, B, &err);
        break;
    default:
        err = GretlError.types;
        break;
    }

    if (err) {
        if (C != null) {
            if (pM != null && *pM == C) {
                *pM = null;
            }
            gretl_matrix_free(C);
            C = null;
        }
    } else {
        /* preserve data-row info? */
        int At1 = gretl_matrix_get_t1(A);
        int At2 = gretl_matrix_get_t2(A);
        int Bt1 = gretl_matrix_get_t1(B);
        int Bt2 = gretl_matrix_get_t2(B);

        if (C.rows == A.rows && At1 >= 0 && At2 > At1) {
            gretl_matrix_set_t1(C, At1);
            gretl_matrix_set_t2(C, At2);
        } else if (C.rows == B.rows && Bt1 >= 0 && Bt2 > Bt1) {
            gretl_matrix_set_t1(C, Bt1);
            gretl_matrix_set_t2(C, Bt2);
        }
    }

 finish:

    if (*pM != null && *pM != C) {
        /* we neither freed nor reused *pM */
        gretl_matrix_free(*pM);
    }

    *pM = C;

    return err;
}

GretlMatrixRaw * tmp_matrix_from_series(Node * n, Parser * p)
{
    int T = sample_size(p.dset);
    const(double) * x = n.v.xvec;
    GretlMatrixRaw * m = null;

    m = gretl_column_vector_alloc(T);

    if (m == null) {
        p.err = GretlError.alloc;
    } else {
        memcpy(m.val, x + p.dset.t1, T * (*x).sizeof);
    }

    return m;
}

/* "Fake" a series using a column vector: the vector must be
   of the same length as the current dataset.
*/

const(double) * get_colvec_as_series(Node * n, int f, Parser * p)
{
    double * ret = null;

    if (n.t != GenSym.mat) {
        p.err = e_types(n);
    } else {
        const(GretlMatrixRaw) * m = n.v.m;

        if (m.rows == p.dset.n && m.cols == 1) {
            ret = cast(double *) m.val;
        } else {
            node_type_error(f, 1, GenSym.series, n, p);
        }
    }

    return ret;
}

/* Similar to get_colvec_as_series(), above, except that (a) we'll
   tolerate a column vector shorter than the length of the current
   dataset provided it equals the length of the current sample range,
   and (b) we return a _copy_ of the vector data, padded out to the
   current full series length with NAs if necessary.
*/

double * copy_colvec_as_series(Node * n, Parser * p)
{
    double * ret = null;

    if (n.t != GenSym.mat) {
        p.err = e_types(n);
    } else {
        const(GretlMatrixRaw) * m = n.v.m;
        size_t sz = p.dset.n * (*ret).sizeof;
        int ns = sample_size(p.dset);
        int i, t, r = m.rows;

        if (m.cols != 1) {
            p.err = GretlError.invarg;
        } else if (r != p.dset.n && r != ns) {
            p.err = GretlError.invarg;
        } else {
            ret = cast(double *) malloc(sz);
            if (ret == null) {
                p.err = GretlError.alloc;
            }
        }
        if (!p.err) {
            if (r == p.dset.n) {
                memcpy(ret, m.val, sz);
            } else {
                for (t=0, i=0; t<p.dset.n; t++) {
                    if (t < p.dset.t1 || t > p.dset.t2) {
                        ret[t] = NADBL;
                    } else {
                        ret[t] = m.val[i++];
                    }
                }
            }
        }
    }

    return ret;
}

/* One of the operands is a matrix (or scalar), the other
   a series: we "cast" the series to a matrix.
*/

Node * matrix_series_calc(Node * l, Node * r, int op, Parser * p)
{
    Node * ret = aux_matrix_node(p);

    if (ret != null && starting(p)) {
        GretlMatrixRaw * a;
        GretlMatrixRaw * b;
        GretlMatrixRaw * tmp;

        if (l.t == GenSym.series) {
            tmp = a = tmp_matrix_from_series(l, p);
            b = node_get_real_matrix(r, p, 0, 0);
        } else {
            a = node_get_real_matrix(l, p, 0, 0);
            tmp = b = tmp_matrix_from_series(r, p);
        }

        if (!p.err) {
            p.err = real_matrix_calc(a, b, op, &ret.v.m);
        }

        gretl_matrix_free(tmp);
    }

    return ret;
}

Node * array_str_calc(Node * l, Node * r, int op, Parser * p)
{
    Node * ret = aux_matrix_node(p);

    if (ret != null && starting(p)) {
        GretlArray * a = l.v.a;

        if (gretl_array_get_type(a) != GretlType.strings) {
            p.err = e_types(l);
        } else {
            int i, n = gretl_array_get_length(a);
            const(char) * si;

            ret.v.m = gretl_zero_matrix_new(1, n);
            if (ret.v.m == null) {
                p.err = GretlError.alloc;
            } else {
                for (i=0; i<n; i++) {
                    si = cast(const(char) *) gretl_array_get_data(a, i);
                    if (op == GenSym.b_doteq) {
                        if (si != null && !strcmp(si, r.v.str)) {
                            ret.v.m.val[i] = 1;
                        }
                    } else if (si == null || strcmp(si, r.v.str)) {
                        ret.v.m.val[i] = 1;
                    }
                }
            }
        }
    }

    return ret;
}

bool comparison_op(int o)
{
    return (o == GenSym.b_eq  || o == GenSym.b_neq ||
            o == GenSym.b_lt  || o == GenSym.b_gt ||
            o == GenSym.b_lte || o == GenSym.b_gte);
}

/* Here we know have a scalar and a 1 x 1 matrix to work with,
   in either order */

Node * matrix_scalar_calc2(Node * l, Node * r, int op,
                           Parser * p)
{
    Node * ret;

    if (scalar_node(l) && scalar_node(r) && comparison_op(op)) {
        ret = aux_scalar_node(p);
    } else if (l.t == GenSym.num && (op == GenSym.b_mod || op == GenSym.b_pow)) {
        /* the matrix on the right is functioning as
           a scalar argument, so produce a scalar
        */
        ret = aux_scalar_node(p);
    } else {
        /* one of the operands is a matrix, albeit 1 x 1,
           so it's safer to produce a matrix result
        */
        ret = aux_sized_matrix_node(p, 1, 1, 0);
    }

    if (!p.err) {
        double x, y;

        if (l.t == GenSym.num) {
            x = l.v.xval;
            y = r.v.m.val[0];
        } else {
            x = l.v.m.val[0];
            y = r.v.xval;
        }

        if (ret.t == GenSym.num) {
            ret.v.xval = xy_calc(x, y, op, GenSym.num, p);
        } else {
            ret.v.m.val[0] = xy_calc(x, y, op, GenSym.mat, p);
        }
    }

    return ret;
}

/* Mixed types: one of the operands is a matrix, the other a scalar,
   giving a matrix result unless we're looking at a comparison
   operator.
*/

Node * matrix_scalar_calc(Node * l, Node * r, int op, Parser * p)
{
    GretlMatrixRaw * m = null;
    int comp = comparison_op(op);
    double x;
    Node * ret = null;

    /* Check for the simple case of scalar and
       1 x 1 matrix, either way round
    */
    if ((l.t == GenSym.num && scalar_node(r)) ||
        (r.t == GenSym.num && scalar_node(l))) {
        return matrix_scalar_calc2(l, r, op, p);
    }

    /* get a scalar @x and matrix @m */
    x = (l.t == GenSym.num)? l.v.xval : r.v.xval;
    m = (l.t == GenSym.mat)? l.v.m : r.v.m;

    if (gretl_is_null_matrix(m)) {
        p.err = GretlError.data;
        return null;
    }

    /* mod, pow: the right-hand term must be scalar */
    if ((op == GenSym.b_mod || op == GenSym.b_pow) && !scalar_node(r)) {
        p.err = e_types(r);
        return null;
    }

    if (comp) {
        ret = aux_scalar_node(p);
    } else if (op == GenSym.b_pow) {
        ret = aux_matrix_node(p);
    } else {
        ret = aux_sized_matrix_node(p, m.rows, m.cols, m.is_complex);
    }

    if (ret == null) {
        return null;
    }

    if (op == GenSym.b_pow) {
        /* note: the (scalar, 1x1) and (1x1, scalar) cases are
           handled above
        */
        double s = node_get_scalar(r, p);

        if (!p.err) {
            ret.v.m = gretl_matrix_pow(m, s, &p.err);
        }
        return ret;
    } else {
        int i, n = m.rows * m.cols;

        if (comp && m.is_complex) {
            /* B_EQ and B_NOTEQ; needs special treatment */
            ret.v.xval = cmatrix_xy_comp(m, x, op, p);
        } else if (comp) {
            /* condition assumed true by until shown false */
            double cond = 1;

            for (i=0; i<n && cond==1; i++) {
                if (l.t == GenSym.num) {
                    cond = xy_calc(x, m.val[i], op, GenSym.mat, p);
                } else {
                    cond = xy_calc(m.val[i], x, op, GenSym.mat, p);
                }
            }
            ret.v.xval = cond;
        } else {
            int xleft = (l.t == GenSym.num);

            if (m.is_complex) {
                p.err = cmatrix_xy_calc(ret.v.m, m, Cmplx(x, 0), xleft, op, p);
            } else {
                p.err = rmatrix_xy_calc(ret.v.m, m, x, xleft, op, p);
            }
            if (gretl_matrix_is_dated(m)) {
                gretl_matrix_set_t1(ret.v.m, gretl_matrix_get_t1(m));
                gretl_matrix_set_t2(ret.v.m, gretl_matrix_get_t2(m));
            }
        }
    }

    return ret;
}

Node * matrix_transpose_node(Node * n, Parser * p)
{
    Node * ret = null;

    if (starting(p)) {
        if (is_tmp_node(n)) {
            /* transpose temp matrix in place */
            if (n.v.m.is_complex) {
                p.err = gretl_ctrans_in_place(n.v.m);
            } else {
                p.err = gretl_matrix_transpose_in_place(n.v.m);
            }
            ret = n;
        } else {
            /* create transpose as new matrix */
            ret = aux_matrix_node(p);
            if (!p.err) {
                if (n.v.m.is_complex) {
                    ret.v.m = gretl_ctrans(n.v.m, 1, &p.err);
                } else {
                    ret.v.m = gretl_matrix_copy_transpose(n.v.m);
                }
                if (ret.v.m == null) {
                    p.err = GretlError.alloc;
                }
            }
        }
    } else {
        ret = is_tmp_node(n) ? n : aux_matrix_node(p);
    }

    return ret;
}

/* We're looking at a string argument that is supposed to represent
   a function call: we'll do a rudimentary heuristic check here.
   FIXME this should be more rigorous?
*/

int is_function_call(const(char) * s)
{
    if (!strchr(s, '(') || !strchr(s, ')')) {
        gretl_errmsg_sprintf("%s: expected a function call", s);
        return 0;
    } else {
        return 1;
    }
}

Node * numeric_jacobian_or_hessian(Node * l, Node * m, Node * r,
                                   int f, Parser * p)
{
    Node * ret = null;

    if (starting(p)) {
        const(char) * s = m.v.str;
        double eps = 0.0;

        if (!is_function_call(s)) {
            p.err = GretlError.types;
            return null;
        }

        ret = aux_matrix_node(p);
        if (ret == null) {
            return null;
        }

        if (!null_node(r)) {
            eps = node_get_scalar(r, p);
        }

        if (!p.err) {
            if (f == GenFunc.f_fdjac) {
                ret.v.m = user_fdjac(l.v.m, s, eps, p.dset, &p.err);
            } else {
                ret.v.m = user_numhess(l.v.m, s, eps, p.dset, &p.err);
            }
        }
    } else {
        ret = aux_matrix_node(p);
    }

    return ret;
}

/* try to get a (named) uservar of type @t from a U_ADDR node */

UserVar * ptr_node_get_uvar(Node * n, int t, Parser * p)
{
    UserVar * uv = null;

    if (n.t == GenSym.u_addr) {
        if (n.L.t != t) {
            p.err = e_types(n.L);
        } else {
            uv = n.L.uv;
            if (uv == null && n.L.vname != null) {
                uv = n.L.uv = get_user_var_by_name(n.L.vname);
            }
            if (uv == null) {
                p.err = GretlError.data;
            }
        }
    } else {
        p.err = e_types(n);
    }

    return uv;
}

GretlMatrixRaw * ptr_node_get_matrix(Node * n, Parser * p)
{
    UserVar * uv = ptr_node_get_uvar(n, GenSym.mat, p);

    return uv != null ? cast(GretlMatrixRaw *) uv.ptr : null;
}

GretlBundle * ptr_node_get_bundle(Node * n, Parser * p)
{
    UserVar * uv = ptr_node_get_uvar(n, GenSym.bundle, p);

    return uv != null ? cast(GretlBundle *) uv.ptr : null;
}

/* To be called ONLY for a node @n whose verified type is U_ADDR */

Node * ptr_node_get_referent_node(Node * n, Parser * p)
{
    Node * ret = n.L;

    if (ret == null) {
        p.err = GretlError.data;
    } else if (ret.t == GenSym.empty) {
        return ret;
    } else if (ret.t == GenSym.series) {
        if (ret.v.xvec == null) {
            int id = current_series_index(p.dset, ret.vname);

            if (id > 0) {
                ret.v.xvec = p.dset.Z[id];
            }
        }
        if (ret.v.xvec == null) {
            p.err = GretlError.data;
            ret = null;
        }
    } else {
        if (ret.uv == null) {
            ret.uv = get_user_var_by_name(ret.vname);
        }
        if (ret.uv == null) {
            p.err = GretlError.data;
            ret = null;
        }
    }

    return ret;
}

/* note: allows @n to be either a regular matrix node or a
   matrix-pointer node
*/

GretlMatrixRaw * mat_node_get_real_matrix(Node * n, Parser * p)
{
    GretlMatrixRaw * m = null;

    if (n == null) {
        p.err = GretlError.data;
    } else if (n.t == GenSym.mat) {
        m = n.v.m;
    } else if (n.t == GenSym.u_addr) {
        m = ptr_node_get_matrix(n, p);
    } else {
        p.err = e_types(n);
    }

    if (!p.err) {
        if (m == null) {
            p.err = GretlError.data;
        } else if (m.is_complex) {
            p.err = GretlError.cmplx;
            m = null;
        }
    }

    return m;
}

const(char) * node_get_fncall(Node * n, Parser * p)
{
    const(char) * ret = null;

    if (n.t != GenSym.str) {
        p.err = e_types(n);
    } else {
        ret = n.v.str;
        if (!is_function_call(ret)) {
            p.err = GretlError.types;
        }
    }

    return ret;
}

Node * BFGS_constrained_max(Node * t, Node * n, Parser * p)
{
    Node * ret = null;
    Node * e = null;
    GretlMatrixRaw * b = null;
    GretlMatrixRaw * bounds = null;
    const(char) * sf = null;
    const(char) * sg = null;
    int i, k = n.v.bn.n_nodes;

    if (k < 3 || k > 4) {
        n_args_error(k, 3, 4, GenFunc.f_bfgscmax, p);
    }

    for (i=0; i<k && !p.err; i++) {
        e = n.v.bn.n[i];
        if (i == 0) {
            b = mat_node_get_real_matrix(e, p);
        } else if (i == 1) {
            bounds = mat_node_get_real_matrix(e, p);
        } else if (i == 2) {
            sf = node_get_fncall(e, p);
        } else if (i == 3 && !null_node(e)) {
            sg = node_get_fncall(e, p);
        }
    }

    if (!p.err) {
        ret = aux_scalar_node(p);
    }

    if (!p.err) {
        int minimize = alias_reversed(t) ? 1 : 0;

        ret.v.xval = user_BFGS(b, sf, sg, p.dset, bounds,
                               minimize, p.prn, &p.err);
    }

    return ret;
}

Node * BFGS_maximize(Node * l, Node * m, Node * r,
                     Parser * p, Node * t)
{
    Node * ret = null;

    if (starting(p)) {
        GretlMatrixRaw * b;
        const(char) * sf = m.v.str;
        const(char) * sg = null;

        b = mat_node_get_real_matrix(l, p);

        if (!p.err) {
            if (r.t == GenSym.str) {
                sg = r.v.str;
            } else if (r.t != GenSym.empty) {
                p.err = e_types(r);
            }
        }
        if (!p.err && !is_function_call(sf)) {
            p.err = GretlError.types;
        }
        if (!p.err && sg != null && !is_function_call(sg)) {
            p.err = GretlError.types;
        }
        if (!p.err && gretl_is_null_matrix(b)) {
            p.err = GretlError.data;
        }
        if (p.err) {
            return null;
        }

        ret = aux_scalar_node(p);
        if (ret != null) {
            int minimize = alias_reversed(t) ? 1 : 0;

            ret.v.xval = user_BFGS(b, sf, sg, p.dset, null,
                                   minimize, p.prn, &p.err);
        }
    } else {
        ret = aux_scalar_node(p);
    }

    return ret;
}

Node * deriv_free_node(Node * l, Node * m, Node * r,
                       Parser * p, Node * t)
{
    Node * ret = null;

    if (starting(p)) {
        GretlMatrixRaw * b = null;
        const(char) * fcall = m.v.str;
        double tol = NADBL;
        int maxit = 0;

        b = mat_node_get_real_matrix(l, p);
        if (!p.err) {
            if (gretl_is_null_matrix(b)) {
                p.err = GretlError.data;
            } else if (!is_function_call(fcall)) {
                p.err = GretlError.types;
            }
        }
        if (!p.err) {
            if (scalar_node(r)) {
                if (t.t == GenFunc.f_gssmax) {
                    tol = r.v.xval;
                } else {
                    maxit = node_get_int(r, p);
                }
            } else if (!null_node(r)) {
                p.err = e_types(r);
            }
        }
        if (!p.err) {
            ret = aux_scalar_node(p);
        }
        if (ret != null) {
            int minimize = alias_reversed(t) ? 1 : 0;
            MaxMethod method = MaxMethod.simann_max;

            if (t.t == GenFunc.f_nmmax) {
                method = MaxMethod.nm_max;
            } else if (t.t == GenFunc.f_gssmax) {
                method = MaxMethod.gss_max;
            }
            ret.v.xval = deriv_free_optimize(method, b, fcall, maxit, tol,
                                             minimize, p.dset, p.prn,
                                             &p.err);
        }
    } else {
        ret = aux_scalar_node(p);
    }

    return ret;
}

Node * fzero_node(Node * l, Node * m, Node * r, Parser * p)
{
    Node * ret = null;

    if (starting(p)) {
        const(char) * fcall = l.v.str;
        GretlMatrixRaw * b = null;
        double tol = NADBL;
        int free_b = 0;

        if (null_or_scalar(m)) {
            b = gretl_matrix_alloc(1, 2);
            b.val[0] = null_node(m) ? NADBL : node_get_scalar(m, p);
            b.val[1] = NADBL;
            free_b = 1;
        } else if (m.t == GenSym.mat) {
            b = m.v.m;
        } else {
            p.err = e_types(m);
        }
        if (!p.err) {
            if (gretl_is_null_matrix(b)) {
                p.err = GretlError.data;
            } else if (!is_function_call(fcall)) {
                p.err = GretlError.types;
            }
        }
        if (!p.err) {
            if (scalar_node(r)) {
                tol = node_get_scalar(r, p);
            } else if (!null_node(r)) {
                p.err = e_types(r);
            }
        }
        if (!p.err) {
            ret = aux_scalar_node(p);
        }
        if (ret != null) {
            ret.v.xval = deriv_free_optimize(MaxMethod.root_find, b, fcall, 0,
                                             tol, 0, p.dset, p.prn,
                                             &p.err);
        }
        if (free_b) {
            gretl_matrix_free(b);
        }
    } else {
        ret = aux_scalar_node(p);
    }

    return ret;
}

void lag_calc(double * y, const(double) * x,
              int k, int t1, int t2,
              int op, double mul,
              Parser * p)
{
    int s, t;

    for (t=t1; t<=t2; t++) {
        s = t - k;
        if (dataset_is_panel(p.dset)) {
            if (s / p.dset.pd != t / p.dset.pd) {
                /* s and t pertain to different units */
                s = -1;
            }
        }
        if (s >= 0 && s < p.dset.n) {
            if (op == GenSym.b_asn && mul == 1.0) {
                y[t] = x[s];
            } else if (op == GenSym.b_asn) {
                y[t] = mul * x[s];
            } else if (op == GenSym.b_add) {
                y[t] += mul * x[s];
            } else {
                p.err = GretlError.data;
            }
        }
    }
}

Node * matrix_file_write(Node * l, Node * m, Node * r, Parser * p)
{
    Node * ret = null;

    if (starting(p)) {
        GretlMatrixRaw * a = node_get_matrix(l, p, 0, 1);
        const(char) * fname = m.v.str;
        int export_ = node_get_bool(r, p, 0);

        if (!p.err) {
            ret = aux_scalar_node(p);
        }
        if (!p.err) {
            p.err = gretl_matrix_write_to_file(a, fname, export_);
            ret.v.xval = p.err;
        }
    } else {
        ret = aux_scalar_node(p);
    }

    return ret;
}

Node * bundle_file_write(Node * l, Node * m, Node * r, Parser * p)
{
    Node * ret = null;

    if (starting(p)) {
        const(char) * s = m.v.str;
        int control = 0;

        if (!null_node(r)) {
            control = node_get_int(r, p);
        }
        ret = aux_scalar_node(p);
        if (ret != null) {
            p.err = cast(int) (ret.v.xval = gretl_bundle_write_to_file(l.v.b, s, control));
        }
    } else {
        ret = aux_scalar_node(p);
    }

    return ret;
}

int check_cswitch_param(GretlMatrixRaw * m, int * k)
{
    int err = 0;

    if (*k < 0 || *k > 4) {
        /* out of bounds */
        err = GretlError.invarg;
    } else if (*k == 1) {
        /* real to complex, column-wise */
        if (m.cols % 2) {
            err = GretlError.nonconf;
        }
    } else if (*k == 3) {
        /* real to complex, row-wise */
        if (m.rows % 2) {
            err = GretlError.nonconf;
        }
    }

    return err;
}

/* matrix on left, scalar(s) on right: returns a matrix */

Node * matrix_scalar_func(Node * l, Node * r,
                          int f, Parser * p)
{
    Node * ret = null;

    if (starting(p)) {
        GretlMatrixRaw * m = l.v.m;
        int k = 0;

        if (f == GenFunc.f_cswitch) {
            k = 1; /* default */
            if (!null_node(r)) {
                k = node_get_int(r, p);
            }
            if (!p.err) {
                p.err = check_cswitch_param(m, &k);
            }
        } else {
            k = node_get_int(r, p);
        }
        if (!p.err && gretl_is_null_matrix(m)) {
            p.err = GretlError.invarg;
        }
        if (!p.err && f == GenFunc.f_msortby && m.is_complex) {
            p.err = GretlError.cmplx;
        }
        if (!p.err) {
            ret = aux_matrix_node(p);
        }
        if (p.err) {
            return null;
        }

        if (f == GenFunc.f_msortby) {
            ret.v.m = gretl_matrix_sort_by_column(m, k-1, &p.err);
        } else if (f == GenFunc.f_cswitch) {
            if (k > 2) {
                /* the old _setcmplx() */
                k = (k == 3)? 1 : 0;
                ret.v.m = gretl_matrix_copy(m);
                if (ret.v.m == null) {
                    p.err = GretlError.alloc;
                } else {
                    p.err = gretl_matrix_set_complex_full(ret.v.m, k);
                }
            } else {
                k = (k == 1)? 1 : 0;
                ret.v.m = gretl_cmatrix_switch(m, k, &p.err);
            }
        }
    } else {
        ret = aux_any_node(p);
    }

    return ret;
}

Node * matrix_vector_func(Node * l, Node * m, Node * r,
                          int f, Parser * p)
{
    Node * ret = null;

    /* at present only F_MSPLITBY comes here */

    if (starting(p)) {
        GretlMatrixRaw * a = node_get_matrix(l, p, 0, 1);
        GretlMatrixRaw * v = node_get_matrix(m, p, 1, 2);
        int colwise = 0;

        if (!p.err) {
            colwise = node_get_bool(r, p, 0);
        }
        if (!p.err) {
            ret = aux_array_node(p);
        }
        if (ret != null) {
            int chunks = m.t == GenSym.num; /* scalar */

            ret.v.a = gretl_matrix_split_by(a, v, colwise, chunks, &p.err);
        }
    } else {
        ret = aux_array_node(p);
    }

    return ret;
}

Node * matrix_string_func(Node * l, Node * m, Node * r,
                          int f, Parser * p)
{
    Node * ret = null;

    /* at present only F_DISTANCE comes here */

    if (starting(p)) {
        GretlMatrixRaw * X = l.v.m;
        GretlMatrixRaw * Y = null;
        const(char) * s = "euc"; /* euclidean default */

        if (!null_node(m)) {
            if (m.t == GenSym.str) {
                s = m.v.str;
            } else {
                node_type_error(f, 2, GenSym.str, m, p);
            }
        }
        if (!null_node(r)) {
            if (r.t == GenSym.mat) {
                Y = r.v.m;
            } else {
                node_type_error(f, 3, GenSym.mat, r, p);
            }
        }
        if (!p.err) {
            ret = aux_matrix_node(p);
        }
        if (ret != null) {
            ret.v.m = distance(X, Y, s, &p.err);
        }
    } else {
        ret = aux_matrix_node(p);
    }

    return ret;
}

/* both operands are known to be matrices or scalars */

Node * matrix_matrix_calc(Node * l, Node * r, int op, Parser * p)
{
    GretlMatrixRaw * ml = null;
    GretlMatrixRaw * mr = null;
    Node * ret = null;

    if ((op == GenSym.b_mul || op == GenSym.b_trmul || op == GenSym.b_add || op == GenSym.b_sub) &&
        l.t == GenSym.mat && r.t == GenSym.mat) {
        ml = l.v.m;
        mr = r.v.m;
        if (ml.is_complex || mr.is_complex) {
            ret = aux_matrix_node(p);
            if (!p.err) {
                p.err = real_matrix_calc(ml, mr, op, &ret.v.m);
            }
            return ret;
        }
    }

    if (op == GenSym.b_dotpow || op == GenSym.b_pow) {
        if (op == GenSym.b_pow) {
            if (scalar_node(l) && scalar_node(r)) {
                op = GenSym.b_dotpow;
            } else if (!scalar_node(r)) {
                p.err = e_types(r);
                return null;
            }
        }
        ret = aux_matrix_node(p);
    } else {
        /* experiment: try reusing aux matrix */
        p.flags |= GenFlags.msave;
        ret = get_aux_node(p, GenSym.mat, 0, NodeFlags.tmp);
        p.flags ^= GenFlags.msave;
    }

    if (ml == null) {
        ml = node_get_matrix(l, p, 0, 1);
        if (op != GenSym.b_pow) {
            mr = node_get_matrix(r, p, 1, 2);
        }
    }

    if (ret != null && starting(p)) {
        if (op == GenSym.b_dotpow) {
            if (ml.is_complex) {
                ret.v.m = gretl_cmatrix_dot_op(ml, mr, '^', &p.err);
            } else {
                ret.v.m = gretl_matrix_dot_op(ml, mr, '^', &p.err);
            }
        } else if (op == GenSym.b_pow) {
            int s = node_get_int(r, p);

            if (!p.err) {
                ret.v.m = gretl_matrix_pow(ml, s, &p.err);
            }
        } else {
            p.err = real_matrix_calc(ml, mr, op, &ret.v.m);
        }
    }

    return ret;
}

Node * matrix_and_or(Node * l, Node * r, int op, Parser * p)
{
    Node * ret = aux_matrix_node(p);

    if (ret != null && starting(p)) {
        const(GretlMatrixRaw) * a = l.v.m;
        const(GretlMatrixRaw) * b = r.v.m;
        int i, n = a.rows * a.cols;

        if (gretl_is_null_matrix(a) || gretl_is_null_matrix(b)) {
            p.err = GretlError.nonconf;
        } else if (a.rows != b.rows || a.cols != b.cols) {
            p.err = GretlError.nonconf;
        } else {
            ret.v.m = gretl_unit_matrix_new(a.rows, a.cols);
            if (ret.v.m == null) {
                p.err = GretlError.alloc;
                return null;
            }
            for (i=0; i<n; i++) {
                if (op == GenSym.b_and) {
                    if (a.val[i] == 0.0 || b.val[i] == 0.0) {
                        ret.v.m.val[i] = 0.0;
                    }
                } else if (op == GenSym.b_or) {
                    if (a.val[i] == 0.0 && b.val[i] == 0.0) {
                        ret.v.m.val[i] = 0.0;
                    }
                }
            }
        }
    }

    return ret;
}

/* both operands are matrices */

Node * matrix_bool(Node * l, Node * r, int op, Parser * p)
{
    Node * ret;

    if (op == GenSym.b_or || op == GenSym.b_and) {
        return matrix_and_or(l, r, op, p);
    }

    ret = aux_scalar_node(p);

    if (ret != null && starting(p)) {
        const(GretlMatrixRaw) * a = l.v.m;
        const(GretlMatrixRaw) * b = r.v.m;
        int i, n = a.rows * a.cols;

        if (gretl_is_null_matrix(a) || gretl_is_null_matrix(b)) {
            ret.v.xval = NADBL;
        } else if (a.rows != b.rows || a.cols != b.cols) {
            ret.v.xval = NADBL;
        } else {
            ret.v.xval = op == GenSym.b_neq ? 0 : 1;
            for (i=0; i<n; i++) {
                if (op == GenSym.b_eq && a.val[i] != b.val[i]) {
                    ret.v.xval = 0;
                    break;
                } else if (op == GenSym.b_lt && a.val[i] >= b.val[i]) {
                    ret.v.xval = 0;
                    break;
                } else if (op == GenSym.b_gt && a.val[i] <= b.val[i]) {
                    ret.v.xval = 0;
                    break;
                } else if (op == GenSym.b_lte && a.val[i] > b.val[i]) {
                    ret.v.xval = 0;
                    break;
                } else if (op == GenSym.b_gte && a.val[i] < b.val[i]) {
                    ret.v.xval = 0;
                    break;
                } else if (op == GenSym.b_neq && a.val[i] != b.val[i]) {
                    ret.v.xval = 1;
                    break;
                }
            }
        }
    }

    return ret;
}

void matrix_error(Parser * p)
{
    if (p.err == 0) {
        p.err = 1;
    }

    if (gretl_errmsg_is_set()) {
        errmsg(p.err, p.prn);
    }
}

/* functions taking a matrix argument and returning a
   scalar result */

Node * matrix_to_scalar_func(Node * l, Node * r, int f, Parser * p)
{
    Node * ret = aux_scalar_node(p);

    if (ret != null && starting(p)) {
        GretlMatrixRaw * m = node_get_matrix(l, p, 0, 0);

        if (m.is_complex && f != GenFunc.f_rows && f != GenFunc.f_cols && f != GenFunc.f_rank) {
            /* gatekeeper for complex */
            p.err = function_real_only(f);
            return ret;
        }

        switch (f) {
        case GenFunc.f_rows:
            ret.v.xval = m.rows;
            break;
        case GenFunc.f_cols:
            ret.v.xval = m.cols;
            break;
        case GenFunc.f_norm1:
            ret.v.xval = gretl_matrix_one_norm(m);
            break;
        case GenFunc.f_infnorm:
            ret.v.xval = gretl_matrix_infinity_norm(m);
            break;
        case GenFunc.f_rcond:
            ret.v.xval = gretl_matrix_rcond(m, &p.err);
            break;
        case GenFunc.f_cnumber:
            ret.v.xval = gretl_matrix_cond_index(m, &p.err);
            break;
        case GenFunc.f_rank:
            if (m.is_complex) {
                ret.v.xval = gretl_cmatrix_rank(m, &p.err);
            } else {
                double eps = NADBL;
                if (!null_node(r)) {
                    eps = node_get_scalar(r, p);
                }
                if (!p.err) {
                    ret.v.xval = gretl_matrix_rank(m, eps, &p.err);
                }
            }
            break;
        default:
            p.err = GretlError.parse;
            break;
        }

        if (p.err) {
            matrix_error(p);
        }
    }

    return ret;
}

/* Compute a value which will be a scalar for a real matrix
   but a complex scalar (2-vector) for a complex matrix:
   handles determinant and trace.
*/

Node * matrix_to_alt_node(Node * n, int f, Parser * p)
{
    GretlMatrixRaw * m = node_get_matrix(n, p, 0, 0);
    Node * ret = null;

    if (!p.err) {
        ret = m.is_complex ? aux_matrix_node(p) : aux_scalar_node(p);
    }

    if (!p.err) {
        if (m.is_complex) {
            if (f == GenFunc.f_trace) {
                ret.v.m = gretl_cmatrix_trace(m, &p.err);
            } else {
                ret.v.m = gretl_cmatrix_determinant(m, f==GenFunc.f_ldet, &p.err);
            }
        } else if (f == GenFunc.f_trace) {
            ret.v.xval = gretl_matrix_trace(m);
        } else {
            int tmpmat = (n.t == GenSym.mat && is_tmp_node(n));

            ret.v.xval =
                user_matrix_get_determinant(m, tmpmat, f==GenFunc.f_ldet, &p.err);
        }
    }

    return ret;
}

Node * matrix_add_names(Node * l, Node * r, int f, Parser * p)
{
    Node * ret = aux_scalar_node(p);

    if (ret != null && starting(p)) {
        GretlMatrixRaw * m = l.v.m;
        int byrow = (f == GenFunc.f_rnameset);

        if (m.is_complex) {
            /* we could set column names for a complex matrix
               but they wouldn't show up on printing
            */
            p.err = GretlError.cmplx;
            return ret;
        }

        if (r.t == GenSym.str) {
            ret.v.xval = umatrix_set_names_from_string(m, r.v.str,
                                                       byrow, p.prn);
        } else if (r.t == GenSym.array) {
            if (gretl_array_get_type(r.v.a) != GretlType.strings) {
                p.err = e_types(r);
            } else {
                ret.v.xval = umatrix_set_names_from_array(m, r.v.a,
                                                          byrow, p.prn);
            }
        } else {
            /* some sort of list-bearing node */
            int * list = node_get_list(r, p);

            if (p.err) {
                ret.v.xval = 1;
            } else {
                ret.v.xval = umatrix_set_names_from_list(m, list, p.dset,
                                                         byrow, p.prn);
            }
            free(list);
        }
    }

    return ret;
}

Node * matrix_get_col_or_row_name(int f, Node * l, Node * r,
                                  Parser * p)
{
    int get_all = null_node(r);
    Node * ret = get_all ? aux_array_node(p) : aux_string_node(p);

    if (ret != null && starting(p)) {
        if (get_all) {
            const(char) ** S;
            int n = 0;

            if (f == GenFunc.f_cnameget) {
                S = gretl_matrix_get_colnames(l.v.m);
                if (S != null) n = l.v.m.cols;
            } else {
                S = gretl_matrix_get_rownames(l.v.m);
                if (S != null) n = l.v.m.rows;
            }
            ret.v.a = gretl_array_from_strings(cast(char **) S, n,
                                               1, &p.err);
        } else {
            int i = node_get_int(r, p);

            if (f == GenFunc.f_cnameget) {
                ret.v.str = user_matrix_get_column_name(l.v.m, i, &p.err);
            } else {
                ret.v.str = user_matrix_get_row_name(l.v.m, i, &p.err);
            }
        }
    }

    return ret;
}

Node * matrix_imhof(Node * l, Node * r, Parser * p)
{
    Node * ret = aux_scalar_node(p);

    if (ret != null && starting(p)) {
        const(GretlMatrixRaw) * m = l.v.m;
        double arg = node_get_scalar(r, p);

        ret.v.xval = imhof(m, arg, &p.err);
    }

    return ret;
}

Node * bkw_node(Node * l, Node * m, Node * r, Parser * p)
{
    Node * ret = aux_matrix_node(p);

    if (ret != null && starting(p)) {
        const(GretlMatrixRaw) * V = l.v.m;
        GretlArray * pnames = null;
        Prn * vprn = null;
        int ns = 0;

        if (!null_node(m)) {
            if (m.t == GenSym.str) {
                /* For compatibility with Lee Adkins' bkw() we accept
                   space-separated parameter names.
                */
                char ** S = gretl_string_split(m.v.str, &ns, null);

                if (S == null) {
                    p.err = GretlError.data;
                } else {
                    pnames = gretl_array_from_strings(S, ns, 0, &p.err);
                }
            } else if (m.t == GenSym.array) {
                if (gretl_array_get_type(m.v.a) != GretlType.strings) {
                    p.err = e_types(m);
                } else {
                    pnames = gretl_array_copy(m.v.a, &p.err);
                }
            } else {
                p.err = e_types(m);
            }
        }

        if (node_get_bool(r, p, 0)) {
            /* optional verbose flag */
            vprn = p.prn;
        }

        if (!p.err) {
            alias BkwFunc = extern (C) GretlMatrixRaw * function(const(GretlMatrixRaw) *, GretlArray *,
                                                            Prn *, int *);
            BkwFunc bkwfunc;

            bkwfunc = cast(BkwFunc) get_plugin_function("bkw_matrix");
            if (bkwfunc == null) {
                p.err = GretlError.fopen;
            } else {
                ret.v.m = (*bkwfunc)(V, pnames, vprn, &p.err);
            }
        }

        gretl_array_destroy(pnames);
    }

    return ret;
}

/* Here we handle the case where the relevant libgretl C function
   overwrites its matrix argument. If @m is just an on-the-fly matrix
   it can be passed on as argument to the function in question, but if
   it's a named user-matrix we'll have to make a copy to pass.
*/

GretlMatrixRaw * do_matrix_ovwrite(GretlMatrixRaw * m,
                                int f, int parm,
                                double xparm,
                                int tmpmat,
                                int * err)
{
    GretlMatrixRaw * R = m;

    if (f == GenFunc.f_chol && !gretl_is_null_matrix(m) &&
        !gretl_matrix_is_symmetric(m)) {
        gretl_errmsg_set(_("Matrix is not symmetric").ptr);
        *err = GretlError.data;
        return null;
    }

    if (!tmpmat) {
        /* shouldn't overwrite @m */
        R = gretl_matrix_copy(m);
        if (R == null) {
            *err = GretlError.alloc;
        }
    }

    if (R != null) {
        if (f == GenFunc.f_chol) {
            *err = gretl_matrix_cholesky_decomp(R);
        } else if (f == GenFunc.f_psdroot) {
            *err = gretl_matrix_psd_root(R, parm);
        } else if (f == GenFunc.f_ginv) {
            *err = gretl_matrix_moore_penrose(R, xparm);
        } else if (f == GenFunc.f_inv) {
            *err = gretl_invert_matrix(R);
        } else if (f == GenFunc.f_upper) {
            *err = gretl_matrix_zero_lower(R);
        } else if (f == GenFunc.f_lower) {
            *err = gretl_matrix_zero_upper(R);
        } else {
            *err = GretlError.data;
        }
        if (*err && R != m) {
            gretl_matrix_free(R);
            R = null;
        }
    }

    return R;
}

void matrix_minmax_indices(int f, int * mm, int * rc, int * idx)
{
    /* mm: 1 for maximum functions, 0 for minimum */
    *mm = (f == GenFunc.f_maxr || f == GenFunc.f_maxc || f == GenFunc.f_imaxr || f == GenFunc.f_imaxc);
    /* rc: 1 for column-wise functions, 0 for row-wise */
    *rc = (f == GenFunc.f_minc || f == GenFunc.f_maxc || f == GenFunc.f_iminc || f == GenFunc.f_imaxc);
    /* idx: 1 for index functions, 0 for values */
    *idx = (f == GenFunc.f_iminr || f == GenFunc.f_iminc || f == GenFunc.f_imaxr || f == GenFunc.f_imaxc);
}

bool mmf_does_complex(int f)
{
    return (f==GenFunc.f_inv || f==GenFunc.f_upper || f==GenFunc.f_lower ||
            f==GenFunc.f_diag || f==GenFunc.f_transp || f==GenFunc.f_ctrans ||
            f==GenFunc.f_vec || f==GenFunc.f_vech || f==GenFunc.f_unvech ||
            f==GenFunc.f_mrev || f== GenFunc.f_fft || f==GenFunc.f_ffti ||
            f==GenFunc.f_cum || f==GenFunc.f_diff || f==GenFunc.f_sumc ||
            f==GenFunc.f_sumr || f==GenFunc.f_prodc || f==GenFunc.f_prodr ||
            f==GenFunc.f_meanc || f==GenFunc.f_meanr || f==GenFunc.f_ginv ||
            f==GenFunc.f_mlog || f==GenFunc.f_mexp || f==GenFunc.f_chol);
}

Node * matrix_to_matrix_func(Node * n, Node * r, int f, Parser * p)
{
    Node * ret = aux_matrix_node(p);

    if (ret != null && starting(p)) {
        GretlMatrixRaw * m = null;
        double xparm = NADBL;
        int tmpmat = 0;
        int parm = 0;
        int a = 0, b = 0, c = 0;

        /* Note: @parm is an integer parameter, required
           for some functions, optional for others. In the
           case of mcov() the default @parm value is 1,
           otherwise it's 0.
        */
        if (f == GenFunc.f_mcov) {
            parm = 1;
        }

        m = node_get_matrix(n, p, 0, 0);
        tmpmat = n.t == GenSym.mat && is_tmp_node(n);

        if (!p.err && m != null && m.is_complex) {
            /* gatekeeper for complex */
            if (!mmf_does_complex(f)) {
                p.err = function_real_only(f);
            }
        }

        if (p.err) {
            goto finalize;
        }

        if (f == GenFunc.f_mrev || f == GenFunc.f_mcov ||
            f == GenFunc.f_psdroot || f == GenFunc.f_vech) {
            /* if present, the @r node should hold an integer */
            if (!null_or_scalar(r)) {
                node_type_error(f, 2, GenSym.num, r, p);
            } else if (!null_node(r)) {
                parm = node_get_int(r, p);
            }
        } else if (f == GenFunc.f_unvech || f == GenFunc.f_ginv) {
            /* if present, the @r node should hold a scalar */
            if (!null_or_scalar(r)) {
                node_type_error(f, 2, GenSym.num, r, p);
            } else if (!null_node(r)) {
                xparm = node_get_scalar(r, p);
            }
        } else if (f == GenFunc.f_nullspc) {
            /* set the 'left' boolean */
            parm = node_get_bool(r, p, 0);
        } else if (f >= GenFunc.f_minc && f <= GenFunc.f_meanr) {
            /* set the skip_na boolean option */
            parm = node_get_bool(r, p, 0);
        } else if (f == GenFunc.f_ranking) {
            if (gretl_vector_get_length(m) == 0) {
                /* m must be a vector */
                p.err = GretlError.types;
            }
        }

        if (!p.err && gretl_is_null_matrix(m) && !emptymat_ok(f)) {
            p.err = GretlError.data;
        }

        if (p.err) {
            goto finalize;
        }

        gretl_error_clear();

        if (gretl_is_null_matrix(m)) {
            ret.v.m = gretl_null_matrix_new();
            goto finalize;
        }

        switch (f) {
        case GenFunc.f_sumc:
            ret.v.m = gretl_matrix_vector_stat(m, GretlVecStat.sum, 0, parm, &p.err);
            break;
        case GenFunc.f_sumr:
            ret.v.m = gretl_matrix_vector_stat(m, GretlVecStat.sum, 1, parm, &p.err);
            break;
        case GenFunc.f_prodc:
            ret.v.m = gretl_matrix_vector_stat(m, GretlVecStat.prod, 0, parm, &p.err);
            break;
        case GenFunc.f_prodr:
            ret.v.m = gretl_matrix_vector_stat(m, GretlVecStat.prod, 1, parm, &p.err);
            break;
        case GenFunc.f_meanc:
            ret.v.m = gretl_matrix_vector_stat(m, GretlVecStat.mean, 0, parm, &p.err);
            break;
        case GenFunc.f_meanr:
            ret.v.m = gretl_matrix_vector_stat(m, GretlVecStat.mean, 1, parm, &p.err);
            break;
        case GenFunc.f_mcov:
            ret.v.m = gretl_covariance_matrix(m, 0, parm, &p.err);
            break;
        case GenFunc.f_mcorr:
            ret.v.m = gretl_covariance_matrix(m, 1, 1, &p.err);
            break;
        case GenFunc.f_cum:
            ret.v.m = gretl_matrix_cumcol(m, &p.err);
            break;
        case GenFunc.f_diff:
            ret.v.m = gretl_matrix_diffcol(m, 0, &p.err);
            break;
        case GenFunc.f_dataok:
            ret.v.m = gretl_matrix_isfinite(m, &p.err);
            break;
        case GenFunc.f_inv:
            if (m.is_complex) {
                ret.v.m = gretl_cmatrix_inverse(m, &p.err);
            } else {
                ret.v.m = do_matrix_ovwrite(m, f, parm, xparm, tmpmat, &p.err);
            }
            break;
        case GenFunc.f_ginv:
            if (m.is_complex) {
                ret.v.m = gretl_cmatrix_ginv(m, &p.err);
            } else {
                ret.v.m = do_matrix_ovwrite(m, f, parm, xparm, tmpmat, &p.err);
            }
            break;
        case GenFunc.f_chol:
            if (m.is_complex) {
                ret.v.m = gretl_cmatrix_cholesky(m, &p.err);
            } else {
                ret.v.m = do_matrix_ovwrite(m, f, parm, xparm, tmpmat, &p.err);
            }
            break;
        case GenFunc.f_psdroot:
        case GenFunc.f_upper:
        case GenFunc.f_lower:
            ret.v.m = do_matrix_ovwrite(m, f, parm, xparm, tmpmat, &p.err);
            break;
        case GenFunc.f_diag:
            ret.v.m = gretl_matrix_get_diagonal(m, &p.err);
            break;
        case GenFunc.f_transp:
            if (m.is_complex) {
                ret.v.m = gretl_ctrans(m, 0, &p.err);
            } else {
                ret.v.m = gretl_matrix_copy_transpose(m);
            }
            break;
        case GenFunc.f_vec:
            ret.v.m = user_matrix_vec(m, &p.err);
            break;
        case GenFunc.f_vech:
            ret.v.m = user_matrix_vech(m, parm, &p.err);
            break;
        case GenFunc.f_unvech:
            ret.v.m = user_matrix_unvech(m, xparm, &p.err);
            break;
        case GenFunc.f_mrev:
            if (parm) {
                ret.v.m = gretl_matrix_reverse_cols(m, &p.err);
            } else {
                ret.v.m = gretl_matrix_reverse_rows(m, &p.err);
            }
            break;
        case GenFunc.f_nullspc:
            if (parm) {
                ret.v.m = gretl_matrix_left_nullspace(m, GretlMatrixMod.none, &p.err);
            } else {
                ret.v.m = gretl_matrix_right_nullspace(m, &p.err);
            }
            break;
        case GenFunc.f_mexp:
            if (m.is_complex) {
                ret.v.m = gretl_cmatrix_exp(m, &p.err);
            } else {
                ret.v.m = gretl_matrix_exp(m, &p.err);
            }
            break;
        case GenFunc.f_mlog:
            ret.v.m = gretl_matrix_log(m, &p.err);
            break;
        case GenFunc.f_fft:
            if (m.is_complex) {
                ret.v.m = gretl_cmatrix_fft(m, 0, &p.err);
            } else {
                ret.v.m = gretl_matrix_fft(m, &p.err);
            }
            break;
        case GenFunc.f_ffti:
            ret.v.m = gretl_matrix_ffti(m, &p.err);
            break;
        case GenFunc.f_polroots:
            ret.v.m = gretl_matrix_polroots(m, 1, 0, &p.err);
            break;
        case GenFunc.f_permute:
            ret.v.m = gretl_permute(m, &p.err);
            break;
        case GenFunc.f_ranking:
            ret.v.m = rank_vector(m, GenFunc.f_sort, &p.err);
            break;
        case GenFunc.f_minc:
        case GenFunc.f_maxc:
        case GenFunc.f_minr:
        case GenFunc.f_maxr:
        case GenFunc.f_iminc:
        case GenFunc.f_imaxc:
        case GenFunc.f_iminr:
        case GenFunc.f_imaxr:
            matrix_minmax_indices(f, &a, &b, &c);
            ret.v.m = gretl_matrix_minmax(m, a, b, c, parm, &p.err);
            break;
        case GenFunc.f_ctrans:
            ret.v.m = gretl_ctrans(m, 1, &p.err);
            break;
        case GenFunc.f_which:
            ret.v.m = gretl_matrix_which(m, &p.err);
            break;
        default:
            break;
        }

        if (ret.v.m == m && n.t == GenSym.mat) {
            /* input matrix was recycled: avoid double-freeing */
            n.v.m = null;
        }

    finalize:

        if (ret.v.m == null) {
            matrix_error(p);
        }
    }

    return ret;
}
