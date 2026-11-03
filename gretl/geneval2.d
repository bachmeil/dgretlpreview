/* Port of the second part (lines 5002-9656) of lib/src/geneval.c: geneval2.c. See geneval0.d for the notes that
   apply to all five parts of the file. */

module gretl.geneval2;


import core.stdc.errno : errno;
import core.stdc.limits : INT_MAX;
import core.stdc.math : fmod, pow, floor, ceil, log, lgamma, nearbyint, exp, atan2, isnan;
import core.stdc.stdio : FILE, fgets, fclose, fprintf, fputs, sprintf, sscanf, stderr;
import core.stdc.stdlib : abort, strtod;
import core.stdc.stdlib : calloc, free, malloc;
import core.stdc.string : memcpy, memset, strcat, strchr, strcmp, strcpy, strlen, strstr;

import gretl.genparseh;
import gretl.geneval0;
import gretl.gretl_errorsh : GretlError;
import gretl.libgretlh;
import core.thread : Thread;
import core.time : usecs;
import gretl.glib : guint32;
import gretl.gretl_matrixh;
import gretl.usermath : MatrixSubspec;
import gretl.uservarh : UserVar;
import gretl.gretl_bfgsh : MaxMethod;
import gretl.geneval1 : aux_scalar_node, aux_ivec_node, aux_matrix_node, aux_mspec_node, aux_string_node, aux_bundle_node, aux_array_node, aux_parent_node, aux_any_node, aux_empty_node;
import gretl.gretl_midash : MidasType;
import gretl.matrix_extrah : MMissingCode;
import gretl.geneval1 : mpi_rank;
import gretl.gstrfuncs : g_utf8_strdown, g_utf8_strup;
import gretl.libseth : SetKey;
import gretl.calendar : easterdate;
import gretl.csvdata : import_csv_as_matrix;
import gretl.dataio : set_dset_matrix_target;
import gretl.dataset : list_info_matrix, series_from_string_transform, series_get_string_table, series_get_string_vals, series_info_bundle, series_is_discrete;
import gretl.dataseth : dataset_is_seasonal, dataset_is_seasonal_panel, dataset_is_time_series, sample_size;
import gretl.describe : correspondence, gretl_corr, gretl_covar;
import gretl.geneval1 : aux_empty_series_node, aux_list_node, aux_node_for_type, aux_series_node, aux_sized_matrix_node, e_types, find_in_tree, gen_type_from_gretl_type, get_aux_node, get_colvec_as_series, gretl_type_from_gen_type, matrix_error, matrix_matrix_calc, node_get_bool, node_get_int, node_get_matrix, node_get_real_matrix, node_get_scalar, ptr_node_get_matrix, ptr_node_get_uvar, strvals_return_check;
import gretl.geneval3 : bundled_series_node;
import gretl.geneval4 : node_is_true, node_type_error;
import gretl.genmain : current_series_index, gretl_is_series;
import gretl.glib : gchar;
import gretl.gmem : g_free;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_errors : errmsg_get_with_default, gretl_errmsg_get, gretl_errmsg_set, gretl_errmsg_sprintf, gretl_error_clear;
import gretl.gretl_intl : gretl_is_ascii;
import gretl.gretl_list : full_var_list, gretl_list_append_list, gretl_list_append_term, gretl_list_cmp, gretl_list_copy, gretl_list_drop, gretl_list_from_vector, gretl_list_get_names, gretl_list_get_names_array, gretl_list_intersection, gretl_list_is_consecutive, gretl_list_new, gretl_list_plus, gretl_list_product, gretl_list_select, gretl_list_sublist, gretl_list_to_vector, gretl_list_union, gretl_null_list, in_gretl_list;
import gretl.gretl_listh : LISTSEP;
import gretl.gretl_matrix : correlated_normal_vec, gretl_invert_symmetric_matrix2, gretl_invpd, gretl_matrix_alloc, gretl_matrix_copy, gretl_matrix_fill, gretl_matrix_free, gretl_matrix_get_colnames, gretl_matrix_get_rownames, gretl_matrix_hdproduct_new, gretl_matrix_inscribe_I, gretl_matrix_lag, gretl_matrix_random_fill, gretl_matrix_set_complex, gretl_matrix_xtab, gretl_matrix_zero, gretl_null_matrix_new, gretl_zero_matrix_new;
import gretl.gretl_paths : gretl_fopen, gretl_make_dotpath, gretl_maybe_switch_dir, gretl_remove;
import gretl.gretl_prn : Prn, pprintf, pputc;
import gretl.gretl_string_table : SeriesTable, gretl_backtick, gretl_getenv, retrieve_date_string, retrieve_date_strings, retrieve_file_content, series_table_get_string;
import gretl.gretl_utils : gretl_mpi_initialized;
import gretl.gretl_xml : gretl_read_gdt;
import gretl.gstrfuncs : g_strdup, g_strdup_printf, g_utf8_strlen, g_utf8_validate;
import gretl.libset : libset_get_int;
import gretl.missingh : NADBL, na;
import gretl.plugins : get_plugin_function;
import gretl.pvalues : lnmgamma;
import gretl.pvaluesh : DistCode;
import gretl.strutils : count_fields, count_lines, gretl_literal_replace, gretl_normalize_varname, gretl_regexp_replace, gretl_strdup, gretl_string_split, gretl_strstrip_copy, gretl_substring, gretl_utf8_select, has_suffix, strings_array_free, strings_array_new;
import gretl.transforms : auto_dummify_list, dumgenr_with_oddval, hf_list_diffgenr, laggenr, laggenr_from_to, list_diffgenr, list_dropcoll, list_dumgenr, list_laggenr, list_loggenr, list_orthdev, list_resample, list_stdgenr, list_xpxgenr;
import gretl.usermath : MSEL_MAX, Msel, SelType, mspec_get_element, mspec_set_col_index, mspec_set_row_index;
import gretl.matrix_extra : gretl_matrix_data_subset, gretl_matrix_from_scalar, gretl_matrix_read_from_file;
import gretl.gretl_cmatrix : cmatrix_get_element, gretl_cmatrix_build, gretl_cmatrix_from_scalar, gretl_cmatrix_hdprod, matrix_is_complex;
import gretl.gretl_array : gretl_array_append_object, gretl_array_copy, gretl_array_copy_subspec, gretl_array_drop_string, gretl_array_from_strings, gretl_array_get_content_type, gretl_array_get_data, gretl_array_get_element, gretl_array_get_length, gretl_array_get_strings, gretl_array_get_type, gretl_array_new, gretl_array_qsort, gretl_array_set_string, gretl_arrays_are_equal, gretl_arrays_concat, gretl_arrays_intersection, gretl_arrays_join, gretl_arrays_union, gretl_matrix_array_flatten, gretl_strings_array_flatten, gretl_strings_array_includes, gretl_strings_array_pos, gretl_strings_sort;
import gretl.gretl_bundle : get_bundle_by_name, gretl_bundle_get_data, gretl_bundle_get_element, gretl_bundle_get_n_members, gretl_bundle_get_scalar, gretl_bundle_read_from_file, gretl_bundle_union, gretl_bundles_are_equal;
import gretl.uservar : user_var_get_type, user_var_get_type_by_name, user_var_set_scalar_value;
import gretl.gretl_typemap : gretl_type_get_name, gretl_type_get_order, gretl_type_get_ref_type;
import gretl.genfuncs : cross_sectional_stat, duration_func, fill_permutation_vector, forecast_stats, gen_time, get_observation_number, list_linear_combo, matrix_fc_stats, midas_linear_combo, seasonals_list, standardize_series, vector_to_midas_list, x_sectional_weighted_stat;
import gretl.nonparam : kendall_tau_func, spearman_rho_func;
import gretl.usermat : check_matrix_subspec, matrix_get_chunk, matrix_get_submatrix, matrix_subspec_new, mspec_get_string, mspec_make_list, user_matrix_eigensym;
import gretl.random : gretl_rand_hex_string, gretl_rand_normal_full, gretl_rand_uniform_minmax;

/* unported dependencies: per-file declarations from the 2026c headers */
extern (C) {
    void current_function_info(const(char) ** funcname, const(char) ** pkgname);
    char * gretl_func_get_arg_name(const(char) * argvar, int * err);
    int gretl_function_depth();
}

Node * glasso_node(Node * l, Node * r, Parser * p)
{
    alias GlassoFunc = extern (C) GretlMatrixRaw * function(const(GretlMatrixRaw) *, GretlBundle *,
                                                       Prn *, int *);
    GlassoFunc glasso;
    Node * ret = null;

    glasso = cast(GlassoFunc) get_plugin_function("gretl_glasso");
    if (glasso == null) {
        p.err = GretlError.fopen;
    } else {
        ret = aux_matrix_node(p);
    }
    if (ret != null) {
        ret.v.m = (*glasso)(l.v.m, r.v.b, p.prn, &p.err);
    }

    return ret;
}

Node * list_reverse_node(Node * n, Parser * p)
{
    Node * ret = aux_list_node(p);

    if (ret != null && starting(p)) {
        int i, nt = n.v.ivec[0];
        int * rev = gretl_list_new(nt);

        for (i=1; i<=nt; i++) {
            rev[i] = n.v.ivec[nt-i+1];
        }
        ret.v.ivec = rev;
    }

    return ret;
}

/* We come here if we got a ".csv" suffix for the argument
   to mread(). If we can open the file as-is, we check to
   make sure it's not a gretl-format matrix file with the
   wrong suffix.
*/

int check_matrix_file(const(char) * fname, int * csv)
{
    char[1024] line = 0;
    FILE * fp;
    int r, c, n;

    *csv = 1;

    fp = gretl_fopen(fname, "rb");
    if (fp == null) {
        /* just assume it's really CSV */
        return 0;
    }

    while (fgets(line.ptr, cast(int) line.length, fp)) {
        if (line[0] != '#') {
            /* heuristic: if the non-comment portion of the file
               starts with two tab-separated integers, it's
               actually a native gretl .mat file regardless of
               the filename suffix?
            */
            n = sscanf(line.ptr, "%d\t%d", &r, &c);
            if (n == 2 && count_fields(line.ptr, "\t") == 2) {
                *csv = 0;
            }
            break;
        }
    }

    fclose(fp);

    return 0;
}

Node * read_object_func(Node * n, Node * r, int f, Parser * p)
{
    Node * ret;

    if (f == GenFunc.f_mread) {
        ret = aux_matrix_node(p);
    } else {
        ret = aux_bundle_node(p);
    }

    if (ret != null && starting(p)) {
        const(char) * fname = n.v.str;
        const(char) * realpath = fname;
        gchar * tmp = null;
        int import_ = node_get_bool(r, p, 0);
        int csv = 0;
        int gdt = 0;

        gretl_error_clear();

        if (import_) {
            tmp = gretl_make_dotpath(fname);
            realpath = tmp;
        }

        if (has_suffix(realpath, ".csv")) {
            p.err = check_matrix_file(realpath, &csv);
            if (p.err) {
                return ret;
            }
        } else if (has_suffix(realpath, ".gdt") ||
                   has_suffix(realpath, ".gdtb")) {
            gdt = 1;
        }

        switch (f) {
        case GenFunc.f_mread:
            if (csv) {
                ret.v.m = import_csv_as_matrix(realpath, &p.err);
            } else if (gdt) {
                set_dset_matrix_target(&ret.v.m);
                p.err = gretl_read_gdt(realpath, null, Gretlopt.OPT_NONE, null);
                set_dset_matrix_target(null);
            } else {
                ret.v.m = gretl_matrix_read_from_file(realpath, 0, &p.err);
            }
            break;
        case GenFunc.f_bread:
            ret.v.b = gretl_bundle_read_from_file(realpath, 0, &p.err);
            break;
        default:
            break;
        }

        g_free(tmp);
    }

    return ret;
}

/* Build a node holding a complex matrix, given two scalars,
   two matrices, or matrix plus scalar.
*/

Node * complex_matrix_node(Node * l, Node * r, Parser * p)
{
    Node * ret = aux_matrix_node(p);

    if (ret != null) {
        GretlMatrixRaw * Re = l.t == GenSym.mat ? l.v.m : null;
        GretlMatrixRaw * Im = r.t == GenSym.mat ? r.v.m : null;
        double x = l.t == GenSym.num ? l.v.xval : 0;
        double y = r.t == GenSym.num ? r.v.xval : 0;

        if (l.t == GenSym.num && null_or_scalar(r)) {
            ret.v.m = gretl_cmatrix_from_scalar(Cmplx(x, y), &p.err);
        } else {
            ret.v.m = gretl_cmatrix_build(Re, Im, x, y, &p.err);
        }
    }

    return ret;
}

/* 2023-02-16: added optional scalar-pointer arg to invpd(),
   to retrieve the log determinant if wanted.
*/

Node * invpd_node(Node * l, Node * r, Parser * p)
{
    Node * ret = aux_matrix_node(p);

    if (ret != null && starting(p)) {
        GretlMatrixRaw * m = null;
        UserVar * uv = null;

        if (l.t == GenSym.mat) {
            if (gretl_is_null_matrix(l.v.m)) {
                p.err = GretlError.data;
            } else {
                m = l.v.m;
            }
        }
        if (!p.err && !null_node(r)) {
            uv = ptr_node_get_uvar(r, GenSym.num, p);
        }
        if (!p.err && l.t == GenSym.mat) {
            if (is_tmp_node(l)) {
                /* OK to overwrite @m */
                ret.v.m = m;
                l.v.m = null;
            } else {
                /* don't destroy @m */
                ret.v.m = gretl_matrix_copy(m);
                if (ret.v.m == null) {
                    p.err = GretlError.alloc;
                }
            }
        }

        if (!p.err) {
            if (l.t == GenSym.num) {
                /* handle scalar case first */
                double x = l.v.xval;

                if (x > 0) {
                    ret.v.m = gretl_matrix_from_scalar(1.0/x);
                    if (uv != null) {
                        user_var_set_scalar_value(uv, log(x));
                    }
                } else {
                    p.err = GretlError.singular;
                }
            } else if (uv != null) {
                double ldet;

                p.err = gretl_invert_symmetric_matrix2(ret.v.m, &ldet);
                user_var_set_scalar_value(uv, ldet);
            } else {
                p.err = gretl_invpd(ret.v.m);
            }
        }

        if (p.err && ret.v.m != null) {
            gretl_matrix_free(ret.v.m);
            ret.v.m = null;
        }
    }

    return ret;
}

Node * matrix_to_matrix2_func(Node * n, Node * r, int f, Parser * p)
{
    Node * ret = aux_matrix_node(p);

    if (ret != null && starting(p)) {
        GretlMatrixRaw * m1 = node_get_matrix(n, p, 0, 0);
        GretlMatrixRaw * m2 = null;

        if (!p.err && gretl_is_null_matrix(m1)) {
            p.err = GretlError.data;
        }
        if (!p.err && !null_node(r)) {
            m2 = ptr_node_get_matrix(r, p);
        }

        if (!p.err) {
            if (f == GenFunc.f_eigsym) {
                ret.v.m = user_matrix_eigensym(m1, m2, &p.err);
            } else if (f == GenFunc.f_hdprod) {
                if (m1.is_complex) {
                    ret.v.m = gretl_cmatrix_hdprod(m1, null, &p.err);
                } else {
                    ret.v.m = gretl_matrix_hdproduct_new(m1, null, &p.err);
                }
            }
        }

        if (ret.v.m == null) {
            matrix_error(p);
        }
    }

    return ret;
}

int ok_matrix_dim(int r, int c, int f)
{
    if (f == GenFunc.f_imat || f == GenFunc.f_zeros || f == GenFunc.f_ones ||
        f == GenFunc.f_munif || f == GenFunc.f_mnorm) {
        /* zero is OK for matrix creation functions, which then
           return an empty matrix
        */
        return (r >= 0 && c >= 0);
    } else {
        double xm = cast(double) r * cast(double) c;

        return (r > 0 && c > 0 && xm < INT_MAX);
    }
}

Node * matrix_fill_func(Node * l, Node * r, int f, Parser * p)
{
    int rows = node_get_int(l, p);
    int n = 0, cols = 0;
    Node * ret = null;

    if (!p.err) {
        if (f == GenFunc.f_randperm) {
            n = rows; /* switched interpretation of first arg */
            rows = 1; /* row vector, per Matlab */
            if (null_node(r)) {
                cols = n;
            } else {
                cols = node_get_int(r, p);
            }
        } else if (f == GenFunc.f_imat && null_node(r)) {
            /* default to square */
            cols = rows;
        } else if (null_node(r)) {
            /* default to a column vector */
            cols = 1;
        } else {
            cols = node_get_int(r, p);
        }
    }

    if (!p.err && !ok_matrix_dim(rows, cols, f)) {
        p.err = GretlError.invarg;
        matrix_error(p);
    }

    if (!p.err) {
        ret = aux_sized_matrix_node(p, rows, cols, 0);
    }

    if (p.err || rows * cols == 0) {
        return ret;
    }

    switch (f) {
    case GenFunc.f_imat:
        if (rows != cols) {
            gretl_matrix_zero(ret.v.m);
        }
        gretl_matrix_inscribe_I(ret.v.m, 0, 0, MIN(rows, cols));
        break;
    case GenFunc.f_zeros:
        gretl_matrix_fill(ret.v.m, 0.0);
        break;
    case GenFunc.f_ones:
        gretl_matrix_fill(ret.v.m, 1.0);
        break;
    case GenFunc.f_munif:
        gretl_matrix_random_fill(ret.v.m, DistCode.uniform);
        break;
    case GenFunc.f_mnorm:
        gretl_matrix_random_fill(ret.v.m, DistCode.normal);
        break;
    case GenFunc.f_randperm:
        p.err = fill_permutation_vector(ret.v.m, n);
        break;
    default:
        break;
    }

    return ret;
}

Node * mvnormal_node(Node * l, Node * r, Parser * p)
{
    Node * ret = null;
    GretlMatrixRaw * L = l.v.m;
    int prec;

    if (L == null || L.rows != L.cols) {
        p.err = GretlError.invarg;
    } else {
        prec = node_get_bool(r, p, -1);
    }
    if (!p.err) {
        ret = aux_sized_matrix_node(p, L.rows, 1, 0);
    }
    if (!p.err && L.rows > 0) {
        p.err = correlated_normal_vec(ret.v.m, L, prec);
    }

    return ret;
}

/* Putative row or column selection matrix: must be a vector;
   cannot contain zero; cannot have both positive and negative
   entries; and entries must be integer-valued.
*/

int set_sel_vector(MatrixSubspec * spec, int r,
                   GretlMatrixRaw * m)
{
    int i, n = gretl_vector_get_length(m);
    int err = 0;

    if (n == 0 && m != null) {
        /* (valid) empty selection */
        if (r) {
            spec.rsel.m = m;
            spec.rtype = SelType.matrix;
        } else {
            spec.lsel.m = m;
            spec.ltype = SelType.matrix;
        }
        return 0;
    }

    if (n > 0) {
        double x;
        int nneg = 0;

        for (i=0; i<n && !err; i++) {
            x = m.val[i];
            if (x == 0 || na(x) || x != floor(x)) {
                err = GretlError.data;
            }
            nneg += x < 0;
        }
        if (!err && nneg > 0 && nneg < n) {
            err = GretlError.data;
        }
    } else {
        err = GretlError.types;
    }

    if (err) {
        gretl_errmsg_set("Invalid selection vector");
    } else if (r) {
        spec.rsel.m = m;
        spec.rtype = SelType.matrix;
    } else {
        spec.lsel.m = m;
        spec.ltype = SelType.matrix;
    }

    return err;
}

/* Basic support for referring to a matrix row or column by
   name rather than by number.
*/

int matrix_index_from_string(const(char) * s,
                             const(GretlMatrixRaw) * m,
                             char pos, Parser * p)
{
    int i, n;
    int ret = 0;

    if (pos == 'L' || pos == 'R') {
        const(char) ** S = null;

        if (pos == 'L') {
            S = gretl_matrix_get_rownames(m);
            n = m.rows;
        } else {
            S = gretl_matrix_get_colnames(m);
            n = m.cols;
        }
        if (S != null) {
            for (i=0; i<n; i++) {
                if (!strcmp(s, S[i])) {
                    ret = i + 1;
                    break;
                }
            }
        }
    } else {
        /* the single-index case */
        const(char) ** SR = gretl_matrix_get_rownames(m);
        const(char) ** SC = gretl_matrix_get_colnames(m);

        if (SR != null) {
            for (i=0; i<m.rows; i++) {
                if (!strcmp(s, SR[i])) {
                    ret = i + 1;
                    break;
                }
            }
        }
        if (ret == 0 && SC != null) {
            for (i=0; i<m.cols; i++) {
                if (!strcmp(s, SC[i])) {
                    ret = i + 1;
                    break;
                }
            }
        }
    }

    if (ret == 0) {
        if (pos == 'L') {
            gretl_errmsg_sprintf(_("'%s': invalid row index").ptr, s);
        } else if (pos == 'R') {
            gretl_errmsg_sprintf(_("'%s': invalid column index").ptr, s);
        } else {
            gretl_errmsg_sprintf(_("'%s': invalid index").ptr, s);
        }
        p.err = GretlError.invarg;
    }

    return ret;
}

int get_scalar_index(Node * n, MatrixSubspec * spec,
                     const(GretlMatrixRaw) * m,
                     char pos, Parser * p)
{
    int idx = 0;

    if (n.t == GenSym.str) {
        if (m == null) {
            /* "can't happen" */
            p.err = GretlError.invarg;
        } else {
            idx = matrix_index_from_string(n.v.str, m, pos, p);
        }
    } else {
        idx = node_get_int(n, p);
        if (!p.err && idx == 0) {
            p.err = GretlError.invarg;
        }
    }

    return idx;
}

void build_mspec(Node * targ, Node * l, Node * r,
                 Node * obj, Parser * p)
{
    MatrixSubspec * spec = targ.v.mspec;
    GretlMatrixRaw * m = null;
    GretlBundle * b = null;
    int lscalar = 0;
    int rscalar = 0;
    int i = 0, j = 0;

    if (spec == null) {
        spec = matrix_subspec_new();
        if (spec == null) {
            p.err = GretlError.alloc;
            return;
        }
    }

    if (obj != null) {
        if (obj.t == GenSym.mat) {
            m = obj.v.m;
        } else if (obj.t == GenSym.bundle) {
            b = obj.v.b;
        }
    }

    /* special case: bundle membership */
    if (b != null && l.t == GenSym.str) {
        if (r == null) {
            spec.ltype = SelType.str;
            spec.rtype = SelType.null_;
            spec.lsel.str = l.v.str;
        } else {
            p.err = e_types(r);
        }
        goto finished;
    }

    lscalar = scalar_node(l) || l.t == GenSym.str;
    rscalar = (r != null && (scalar_node(r) || r.t == GenSym.str));

    if (lscalar) {
        char pos = r == null ? 'S' : 'L';

        i = get_scalar_index(l, spec, m, pos, p);
        if (!p.err && r == null && i > 0) {
            /* identify and flag the single index case */
            spec.ltype = SelType.single;
            spec.rtype = SelType.null_;
            mspec_set_row_index(spec, i);
            goto finished;
        }
    }

    if (!p.err && rscalar) {
        j = get_scalar_index(r, spec, m, 'R', p);
    }

    if (p.err) {
        goto finished;
    }

    if (l.t == GenSym.dum) {
        if (r != null) {
            p.err = GretlError.invarg;
        } else {
            spec.rtype = SelType.all;
            if (l.v.idnum == DumId.diag) {
                spec.ltype = SelType.diag;
            } else if (l.v.idnum == DumId.upper) {
                spec.ltype = SelType.upper;
            } else if (l.v.idnum == DumId.lower) {
                spec.ltype = SelType.lower;
            } else if (l.v.idnum == DumId.real_) {
                spec.ltype = SelType.real_;
            } else if (l.v.idnum == DumId.imag) {
                spec.ltype = SelType.imag;
            } else {
                p.err = GretlError.types;
            }
        }
        goto finished;
    } else if (i > 0 && j > 0) {
        spec.ltype = spec.rtype = SelType.element;
        mspec_set_row_index(spec, i);
        mspec_set_col_index(spec, j);
        goto finished;
    } else if (lscalar) {
        spec.ltype = i > 0 ? SelType.range : SelType.excl;
        mspec_set_row_index(spec, i);
    } else if (l.t == GenSym.ivec) {
        spec.ltype = SelType.range;
        spec.lsel.range[0] = l.v.ivec[0];
        spec.lsel.range[1] = l.v.ivec[1];
    } else if (l.t == GenSym.mat) {
        p.err = set_sel_vector(spec, 0, l.v.m);
    } else if (null_node(l)) {
        spec.ltype = SelType.all;
    } else {
        p.err = e_types(l);
        goto finished;
    }

    if (r == null) {
        spec.rtype = SelType.null_;
    } else if (rscalar) {
        spec.rtype = j > 0 ? SelType.range : SelType.excl;
        mspec_set_col_index(spec, j);
    } else if (r.t == GenSym.ivec) {
        spec.rtype = SelType.range;
        spec.rsel.range[0] = r.v.ivec[0];
        spec.rsel.range[1] = r.v.ivec[1];
    } else if (r.t == GenSym.mat) {
        p.err = set_sel_vector(spec, 1, r.v.m);
    } else if (null_node(r)) {
        spec.rtype = SelType.all;
    } else {
        p.err = e_types(r);
    }

 finished:

    if (p.err && spec != null) {
        free(spec);
        spec = null;
    }

    targ.v.mspec = spec;
}

/* Node holding evaluated result of sub-object specification.  The nodes
   @l and @r are the left-hand and right-hand parts of the slice
   specification @s, and they "typically" refer to the row and column
   dimensions of a matrix, respectively. But we have to allow for the
   possibility that @r is NULL (as in a slice of a vector), in which
   case @l could refer either to rows or columns. We also allow that
   the object to be sliced may be a bundle or array, in which case the
   slice must always be one-dimensional (@l only).

   To support use of column or row names as matrix indices we want to
   give the builder function build_mspec() access to the matrix in
   question, which requires finding it. It will surely be in the
   parentage of @s, but exactly where depends on whether the "genr"
   expression references the matrix directly, or indirectly as a member
   of a bundle (BMEMB) or element of an array of matrices (OSL).
*/

Node * object_slice_node(Node * s, Node * l, Node * r, Parser * p)
{
    Node * ret = null;

    if (l.t == GenSym.undef || (r != null && r.t == GenSym.undef)) {
        p.err = GretlError.data;
    } else {
        Node * obj = s.parent.L;

        if (obj != null) {
            if (obj.t == GenSym.mat || obj.t == GenSym.bundle || obj.t == GenSym.array) {
                /* OK */
            } else if (obj.t == GenSym.bmemb || obj.t == GenSym.osl || func_node(obj)) {
                if (obj.aux != null) {
                    obj = obj.aux;
                }
            }
        }
        ret = aux_mspec_node(p);
        if (ret != null && starting(p)) {
            build_mspec(ret, l, r, obj, p);
        }
    }

    return ret;
}

Node * submatrix_node(Node * l, Node * r, Parser * p)
{
    Node * ret = null;

    if (starting(p)) {
        MatrixSubspec * spec = r.v.mspec;
        GretlMatrixRaw * m = node_get_matrix(l, p, 0, 0);

        if (!p.err) {
            if (spec.ltype == SelType.str) {
                p.err = GretlError.types;
            } else {
                p.err = check_matrix_subspec(spec, m);
            }
        }
        if (!p.err) {
            ret = aux_matrix_node(p);
        }

        if (!p.err) {
            if (spec.ltype == SelType.contig) {
                ret.v.m = matrix_get_chunk(m, spec, &p.err);
            } else if (spec.ltype == SelType.element) {
                int i = mspec_get_element(spec);

                if (m.is_complex) {
                    ret.v.m = cmatrix_get_element(m, i, &p.err);
                } else {
                    /* 2020-12-29: don't collapse to scalar here */
                    ret.v.m = gretl_matrix_alloc(1,1);
                    ret.v.m.val[0] = m.val[i];
                }
            } else {
                ret.v.m = matrix_get_submatrix(m, spec, 1, &p.err);
            }
        }
    } else {
        ret = aux_any_node(p);
    }

    return ret;
}

/* Check a list that has been stored in a bundle or array to see
   if it can be interpreted as a list given the characteristics
   of the current dataset (or lack thereof). Return 1 if this
   is so, otherwise 0;
*/

int stored_list_ok(const(int) * list, const(DatasetRaw) * dset)
{
    int i;

    if (dset == null || dset.n == 0) {
        return 0;
    }

    for (i=1; i<=list[0]; i++) {
        if (list[i] >= dset.v ||
            (list[i] < 0 && list[i] != LISTSEP)) {
            return 0;
        }
    }

    return 1;
}

Node * array_element_node(GretlArray * a, int i,
                          Parser * p)
{
    Node * ret = null;
    GretlType type = cast(GretlType) 0;
    void * data;

    data = gretl_array_get_element(a, i-1, &type, &p.err);

    if (p.err == GretlError.invarg) {
        gretl_errmsg_sprintf(_("Index value %d is out of bounds").ptr, i);
    } else {
        ret = aux_node_for_type(type, p);
    }

    if (!p.err) {
        if (type == GretlType.string_) {
            ret.v.str = cast(char *) data;
        } else if (type == GretlType.matrix) {
            ret.v.m = cast(GretlMatrixRaw *) data;
        } else if (type == GretlType.bundle) {
            ret.v.b = cast(GretlBundle *) data;
        } else if (type == GretlType.array) {
            ret.v.a = cast(GretlArray *) data;
        } else if (type == GretlType.double_) {
            /* FIXME not reached? */
            ret.v.xval = *cast(double *) data;
        } else if (type == GretlType.list) {
            /* last revised 2022-05-06 */
            if (stored_list_ok(cast(const(int) *) data, p.dset)) {
                ret.v.ivec = cast(int *) data;
            } else {
                /* fallback: extract list as row vector */
                ret.t = GenSym.mat;
                ret.flags |= NodeFlags.tmp;
                ret.v.m = gretl_list_to_vector(cast(const(int) *) data,
                                               &p.err);
            }
        }
    }

    return ret;
}

Node * array_subspec_node(GretlArray * a, int * list,
                          Parser * p)
{
    Node * ret = aux_array_node(p);

    if (ret != null) {
        ret.v.a = gretl_array_copy_subspec(a, list, &p.err);
    }

    return ret;
}

Node * list_range_node(int * list, int r1, int r2, Parser * p)
{
    Node * ret = null;

    if (starting(p)) {
        ret = aux_list_node(p);
        if (ret != null) {
            ret.v.ivec = gretl_list_sublist(list, r1, r2);
            if (ret.v.ivec == null) {
                p.err = GretlError.alloc;
            }
        }
    } else {
        ret = aux_any_node(p);
    }

    return ret;
}

Node * string_range_node(const(char) * s, int r1, int r2, Parser * p)
{
    Node * ret = null;

    if (starting(p)) {
        ret = aux_string_node(p);
        if (ret != null) {
            ret.v.str = gretl_substring(s, r1, r2, &p.err);
        }
    } else {
        ret = aux_any_node(p);
    }

    return ret;
}

/* Return a node holding a shallow copy of series @vnum in p->dset */

Node * real_list_series_node(int vnum, Parser * p)
{
    Node * ret = aux_empty_series_node(p);

    if (!p.err) {
        /* overwrite the TMP_NODE flag on @ret, since we're
           borrowing a dset->Z member!
        */
        ret.flags = NodeFlags.aux;
        ret.vnum = vnum;
        ret.v.xvec = p.dset.Z[vnum];
    }

    return ret;
}

/* Handling a context where we have a list under node @l,
   and position within the list in @pos.
*/

Node * list_member_by_pos(Node * l, int pos, Parser * p)
{
    Node * ret = null;
    int * list = l.v.ivec;

    if (pos < 1 || pos > list[0]) {
        gretl_errmsg_sprintf(_("Index value %d is out of bounds").ptr, pos);
        p.err = GretlError.invarg;
    } else {
        int v = list[pos];

        if (v < 0 || v >= p.dset.v) {
            gretl_errmsg_sprintf(_("Variable number %d is out of bounds").ptr, v);
            p.err = GretlError.invarg;
        } else {
            ret = real_list_series_node(v, p);
        }
    }

    return ret;
}

/* Handling a context where we have a list under node @l,
   and the name of a (putative) list member under node @r.
*/

Node * list_member_by_name(Node * l, Node * r, Parser * p)
{
    Node * ret = null;
    int * list = l.v.ivec;
    const(char) * s = r.v.str;
    int i, v;

    if (p.dset == null) {
        p.err = GretlError.data;
        return null;
    }

    for (i=1; i<=list[0]; i++) {
        v = list[i];
        if (!strcmp(s, p.dset.varname[v])) {
            /* found the named member series */
            ret = real_list_series_node(v, p);
            break;
        }
    }

    if (ret == null) {
        p.err = GretlError.unkvar;
    }

    return ret;
}

int mspec_get_series_index(MatrixSubspec * s,
                           Parser * p)
{
    int t = -1;

    if (s.ltype == SelType.single) {
        t = s.lsel.range[0];
    } else if (s.ltype == SelType.range && s.rtype == SelType.null_) {
        if (s.lsel.range[0] == s.lsel.range[1]) {
            t = s.lsel.range[0];
        } else {
            /* allow for dates such as "2008:4" */
            gchar * tmp;

            tmp = g_strdup_printf("%d:%d", s.lsel.range[0],
                                  s.lsel.range[1]);
            t = get_observation_number(tmp, p.dset);
            g_free(tmp);
        }
    }

    if (t < 1 || t > p.dset.n) {
        p.err = GretlError.data;
    }

    return t;
}

int mspec_get_array_index(MatrixSubspec * spec,
                          int * err)
{
    int idx = 0;

    if (spec.ltype == SelType.single) {
        idx = spec.lsel.range[0];
    } else if (spec.ltype == SelType.range &&
               spec.rtype == SelType.null_ &&
               spec.lsel.range[0] == spec.lsel.range[1]) {
        idx = spec.lsel.range[0];
    } else {
        gretl_errmsg_set("Invalid left-hand side index value");
        *err = GretlError.types;
    }

    return idx;
}

/* stricter variant of test_for_single_range */

int get_single_element(MatrixSubspec * spec,
                       Parser * p)
{
    int ret = 0;

    if (spec.ltype == SelType.single) {
        ret = spec.lsel.range[0];
    } else {
        if (p != null) {
            p.err = GretlError.types;
        }
        ret = -1;
    }

    return ret;
}

void * sub_addr_get_data(Node * t, GretlType * ptype,
                         UserVar ** puv)
{
    Node * l = t.L;
    Node * r = t.R;
    GretlArray * a = cast(GretlArray *) l.v.ptr;
    GretlType type = cast(GretlType) 0;
    void * elem;
    int idx, err = 0;

    idx = get_single_element(r.v.mspec, null);
    elem = gretl_array_get_element(a, idx-1, &type, &err);
    *ptype = (elem == null)? GretlType.none :
        gretl_type_get_ref_type(type);
    *puv = l.uv;

    return elem;
}

Node * process_OSL_address(Node * t, Node * l, Node * r, Parser * p)
{
    int idx = get_single_element(r.v.mspec, p);
    Node * lb = l.L;
    Node * ret = null;

    if (lb.t != GenSym.osl || lb.uv == null || idx <= 0) {
        p.err = GretlError.types;
    } else {
        GretlType type = user_var_get_type(lb.uv);

        if (type == GretlType.array) {
            ret = aux_parent_node(p);
            if (ret != null) {
                ret.t = GenSym.sub_addr;
                ret.L = lb; /* extracted left-hand */
                ret.R = r;  /* evaluated right-hand */
                /* prevent double-freeing of children @l and @r */
                ret.flags |= NodeFlags.lht;
            }
        } else {
            p.err = GretlError.types;
        }
    }

    if (p.err) {
        gretl_errmsg_set(_("Wrong type of operand for unary '&'").ptr);
    }

    return ret;
}

int want_singleton_array(Node * n, Parser * p)
{
    if (p.aux != null && p.aux.t == GenSym.array) {
        GretlType t = gretl_array_get_type(n.v.a);

        /* We want to preserve the ARRAY type of the aux
           node associated with @n: this requires producing
           a singleton array unless we're looking at an
           array of arrays.
        */
        return t != GretlType.arrays;
    } else if (p.aux == null && p.targ == GenSym.array &&
               p.tree == n.parent) {
        return 1;
    }

    return 0;
}

int * array_subspec_list(Node * l, Node * r, Parser * p)
{
    MatrixSubspec * spec = r.v.mspec;
    int * list = null;

    if (spec.rtype != SelType.null_) {
        /* array selection must be one-dimensional */
        p.err = GretlError.invarg;
    } else {
        int len;

        if (l.t == GenSym.array) {
            len = gretl_array_get_length(l.v.a);
        } else if (l.t == GenSym.str) {
            len = cast(int) g_utf8_strlen(l.v.str, -1);
        } else if (l.t == GenSym.list) {
            len = l.v.ivec[0];
        } else {
            p.err = e_types(l);
            return null;
        }

        /* convert @spec to list of elements */
        list = mspec_make_list(spec.ltype, &spec.lsel,
                               len, &p.err);
    }

    return list;
}

Node * subobject_node(Node * l, Node * r, Parser * p)
{
    Node * ret = null;

    if (!starting(p)) {
        return aux_any_node(p);
    }

    if (r == null || r.t != GenSym.mspec) {
        p.err = e_types(r);
    } else if (l.t == GenSym.mat) {
        ret = submatrix_node(l, r, p);
    } else if (l.t == GenSym.array) {
        int * vlist = array_subspec_list(l, r, p);

        if (!p.err && vlist[0] == 1) {
            if (want_singleton_array(l, p)) {
                /* produce a 1-element array */
                ret = array_subspec_node(l.v.a, vlist, p);
            } else {
                /* extract an array element */
                ret = array_element_node(l.v.a, vlist[1], p);
            }
        } else if (!p.err) {
            ret = array_subspec_node(l.v.a, vlist, p);
        }
        free(vlist);
    } else if (l.t == GenSym.list || l.t == GenSym.str) {
        int * vlist = array_subspec_list(l, r, p);

        if (!p.err && vlist[0] == 1 && l.t == GenSym.list) {
            /* singleton list selection */
            ret = list_member_by_pos(l, vlist[1], p);
        } else if (!p.err && gretl_list_is_consecutive(vlist)) {
            /* selected elements are consecutive */
            int r1 = vlist[1];
            int r2 = vlist[vlist[0]];

            if (l.t == GenSym.list) {
                ret = list_range_node(l.v.ivec, r1, r2, p);
            } else {
                ret = string_range_node(l.v.str, r1, r2, p);
            }
        } else if (!p.err) {
            /* selection is not consecutive */
            if (l.t == GenSym.str) {
                ret = aux_string_node(p);
                ret.v.str = gretl_utf8_select(l.v.str, vlist);
            } else {
                ret = aux_list_node(p);
                ret.v.ivec = gretl_list_select(l.v.ivec, vlist, &p.err);
            }
        }
        free(vlist);
    } else if (l.t == GenSym.series) {
        int t = mspec_get_series_index(r.v.mspec, p);

        if (!p.err) {
            ret = aux_scalar_node(p);
            if (!p.err) {
                ret.v.xval = l.v.xvec[t-1];
            }
        }
    } else if (l.t == GenSym.bundle) {
        /* the "mspec" must hold a single key string */
        const(char) * key = mspec_get_string(r.v.mspec, 0);
        GretlType type = GretlType.none;
        void * val = null;
        int is_virtual = 0;
        int is_tmp = 0;

        if (key == null) {
            p.err = GretlError.types;
        } else {
            val = gretl_bundle_get_element(l.v.b, key, &type,
                                           &is_virtual, &p.err);
        }
        if (!p.err) {
            int t = gen_type_from_gretl_type(type);

            if (t == GenSym.num) {
                ret = aux_scalar_node(p);
                if (!p.err) {
                    ret.v.xval = gretl_bundle_get_scalar(l.v.b, key, null);
                }
            } else if (t == GenSym.series) {
                GretlMatrixRaw * m = cast(GretlMatrixRaw *) val;

                ret = bundled_series_node(l, m, is_virtual, &is_tmp, p);
            } else {
                ret = get_aux_node(p, t, 0, 0);
                if (!p.err) {
                    ret.v.ptr = val;
                }
            }
            if (is_tmp) {
                ret.flags |= NodeFlags.tmp;
            } else {
                ret.flags = cast(ubyte) (ret.flags & ~NodeFlags.tmp);
            }
            ret.flags |= NodeFlags.mut;
        }
    } else if (l.t == GenSym.num) {
        /* allow "indexing into" a scalar, but only for a
           single index with value 1
        */
        int i = get_single_element(r.v.mspec, p);

        if (i == 1) {
            ret = aux_scalar_node(p);
            if (!p.err) {
                ret.v.xval = l.v.xval;
            }
        } else {
            p.err = GretlError.types;
        }
    } else {
        p.err = e_types(l);
    }

    return ret;
}

Node * process_subslice(Node * l, Node * r, Parser * p)
{
    Node * ret = null;

    if (starting(p)) {
        if (!scalar_node(l)) {
            p.err = e_types(l);
        } else if (!null_or_scalar(r)) {
            p.err = e_types(r);
        } else {
            ret = aux_ivec_node(p, 2);
            if (ret != null) {
                ret.v.ivec[0] = node_get_int(l, p);
                if (null_node(r)) {
                    ret.v.ivec[1] = MSEL_MAX; /* placeholder */
                } else {
                    ret.v.ivec[1] = node_get_int(r, p);
                }
            }
        }
    } else {
        ret = aux_ivec_node(p, 2);
    }

    return ret;
}

/* Note: many standard and a few non-standard math functions
   are not included in the switch below, because pointers to
   the functions are saved, obviating the need for repeated
   lookup. See the @ptrfuncs mechanism in genlex.c.
*/

double real_apply_func(double x, int f, Parser * p)
{
    double y;

    errno = 0;

    if (na(x)) {
        switch (f) {
        case GenFunc.f_missing:
            return 1.0;
        case GenFunc.f_dataok:
        case GenFunc.f_misszero:
            return 0.0;
        default:
            return NADBL;
        }
    }

    switch (f) {
    case GenSym.u_neg:
        return -x;
    case GenSym.u_pos:
        return x;
    case GenSym.u_not:
        return x == 0;
    case GenFunc.f_toint:
        if (x > int.max || x < int.min) {
            return NADBL;
        } else {
            return cast(double) cast(int) x;
        }
    case GenFunc.f_missing:
        return 0.0;
    case GenFunc.f_dataok:
        return 1.0;
    case GenFunc.f_misszero:
        return x;
    case GenFunc.f_zeromiss:
        return (x == 0.0)? NADBL : x;
    case GenFunc.f_easter:
        y = easterdate(cast(int) x);
        return y;
        /* below: functions that should already be mapped;
           it should be possible to delete them
        */
    case GenFunc.f_log: /* in case it's aliased */
        return log(x);
    default:
        return 0.0;
    }
}

/* @n must be of type NUM, MAT or SERIES, pre-checked. If
   it's MAT and 1 x 1 we treat it as if it were a scalar.
*/

double node_get_double(Node * n, int i, Parser * p)
{
    if (n.t == GenSym.num) {
        return n.v.xval;
    } else if (n.t == GenSym.mat) {
        if (n.v.m.rows == 1 && n.v.m.cols == 1) {
            return n.v.m.val[0];
        } else {
            return n.v.m.val[i];
        }
    } else {
        return n.v.xvec[p.dset.t1 + i];
    }
}

/* @n must be of type NUM, MAT or SERIES, pre-checked */

void node_set_double(Node * n, int i, double x, Parser * p)
{
    if (n.t == GenSym.num) {
        n.v.xval = x;
    } else if (n.t == GenSym.mat) {
        n.v.m.val[i] = x;
    } else {
        n.v.xvec[p.dset.t1 + i] = x;
    }
}

int is_integral(double d)
{
    return d == floor(d) || d == ceil(d);
}

double bincoeff(double n, double k, int * err)
{
    double ret;

    if (n < k || k < 0) {
        *err = GretlError.invarg;
        return NADBL;
    }

    /* catch special cases first */
    if (n == k || k == 0) {
        ret = 1.0;
    } else if (n - k == 1 || k == 1) {
        ret = n;
    } else {
        ret = lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1);
        if (is_integral(n) && is_integral(k)) {
            ret = nearbyint(exp(ret));
        } else {
            ret = exp(ret);
        }
    }

    return ret;
}

/* hard or soft thresholding */

double threshold(double x, double a, int hard, int * err)
{
    double ret = 0.0;

    if (a < 0) {
        *err = GretlError.invarg;
        return NADBL;
    }

    if (x > a) {
        ret = hard ? x : x-a;
    } else if (x < -a) {
        ret = hard ? x : x+a;
    }

    return ret;
}

/* flexible_2arg_node() handles cases like atan2, where we have two
   possibly heterogeneous arguments (scalar, series, matrix) and the
   objective is to return a sensibly sized object. The @flag argument
   is used for minimal inflection of the function, as in thresh().
*/

Node * flexible_2arg_node(Node * l, Node * r, int f, int flag,
                          Parser * p)
{
    Node * ret = null;
    int rettype = 0;
    int nmin, nmax;
    int nl = 0, nr = 0;

    if (l.t == GenSym.num) {
        nl = 1;
    } else if (l.t == GenSym.mat) {
        nl = gretl_vector_get_length(l.v.m);
    } else {
        nl = sample_size(p.dset);
    }

    if (r.t == GenSym.num) {
        nr = 1;
    } else if (r.t == GenSym.mat) {
        nr = gretl_vector_get_length(r.v.m);
    } else {
        nr = sample_size(p.dset);
    }

    nmin = nr < nl ? nr : nl;
    nmax = nr > nl ? nr : nl;

    if (nmin == 0 || (nmin > 1 && nmax != nmin)) {
        p.err = GretlError.nonconf;
        return null;
    }

    /* ordering is MAT > SERIES > NUM */
    rettype = r.t > l.t ? r.t : l.t;

    if (rettype == GenSym.num) {
        ret = aux_scalar_node(p);
    } else if (rettype == GenSym.mat) {
        ret = aux_sized_matrix_node(p, nmax, 1, 0);
    } else {
        ret = aux_series_node(p);
    }

    if (ret != null && !p.err) {
        double x1, x2, y = NADBL;
        int i;

        for (i=0; i<nmax; i++) {
            x1 = node_get_double(l, i, p);
            x2 = node_get_double(r, i, p);
            if (f == GenFunc.f_atan2) {
                y = atan2(x1, x2);
            } else if (f == GenFunc.f_bincoeff) {
                y = bincoeff(x1, x2, &p.err);
                if (p.err) {
                    break;
                }
            } else if (f == GenFunc.f_lnmgamma) {
                y = lnmgamma(x1, cast(int) x2);
                if (p.err) {
                    break;
                }
            } else if (f == GenFunc.f_threshold) {
                y = threshold(x1, x2, flag, &p.err);
                if (p.err) {
                    break;
                }
            }
            node_set_double(ret, i, y, p);
        }
    }

    return ret;
}

/* the C type of the function pointers stored in the nodes of "pointerized" math functions */
alias DblFunc = extern (C) double function(double);

Node * apply_scalar_func(Node * n, Node * f, Parser * p)
{
    Node * ret = aux_scalar_node(p);

    if (ret != null) {
        DblFunc dfunc = cast(DblFunc) f.v.ptr;

        if (dfunc != null) {
            ret.v.xval = dfunc(n.v.xval);
        } else {
            ret.v.xval = real_apply_func(n.v.xval, f.t, p);
        }
    }

    return ret;
}

enum G_USEC_PER_SEC = 1000000;

Node * sleep_node(Node * n, int f, Parser * p)
{
    Node * ret = aux_scalar_node(p);

    if (ret != null) {
        double secs = node_get_scalar(n, p);

        if (secs < 0 || secs * G_USEC_PER_SEC >= cast(double) ulong.max) {
            p.err = GretlError.invarg;
        } else {
            /* g_usleep() of glib */
            Thread.sleep(usecs(cast(long) (G_USEC_PER_SEC * secs)));
            ret.v.xval = 0;
        }
    }

    return ret;
}

Node * scalar_isnan_node(Node * n, Parser * p)
{
    Node * ret = aux_scalar_node(p);

    if (ret != null) {
        double x = node_get_scalar(n, p);

        if (!p.err) {
            ret.v.xval = isnan(x) != 0;
        }
    }

    return ret;
}

Node * matrix_isnan_node(Node * n, Parser * p)
{
    Node * ret = aux_matrix_node(p);

    if (ret != null && starting(p)) {
        const(GretlMatrixRaw) * m = n.v.m;

        if (m.rows == 0 || m.cols == 0) {
            p.err = GretlError.data;
        } else {
            ret.v.m = gretl_matrix_alloc(m.rows, m.cols);
            if (ret.v.m == null) {
                p.err = GretlError.alloc;
            } else {
                int i, nn = m.rows * m.cols;

                for (i=0; i<nn; i++) {
                    ret.v.m.val[i] = isnan(m.val[i]) != 0;
                }
                gretl_matrix_set_complex(ret.v.m, m.is_complex);
            }
        }
    }

    return ret;
}

Node * apply_series_func(Node * n, Node * f, Parser * p)
{
    Node * ret = aux_series_node(p);
    int t;

    if (ret != null) {
        DblFunc dfunc = cast(DblFunc) f.v.ptr;
        const(double) * x;

        if (n.t == GenSym.series) {
            if (n.v.xvec == null) {
                fprintf(stderr, "apply_series_func: got NULL xvec\n");
                p.err = GretlError.data;
            } else {
                x = n.v.xvec;
            }
        } else {
            x = get_colvec_as_series(n, f.t, p);
        }

        if (!p.err) {
            if (autoreg(p)) {
                if (dfunc != null) {
                    ret.v.xvec[p.obs] = dfunc(x[p.obs]);
                } else {
                    ret.v.xvec[p.obs] = real_apply_func(x[p.obs], f.t, p);
                }
            } else if (dfunc != null) {
                for (t=p.dset.t1; t<=p.dset.t2; t++) {
                    ret.v.xvec[t] = dfunc(x[t]);
                }
            } else {
                for (t=p.dset.t1; t<=p.dset.t2; t++) {
                    ret.v.xvec[t] = real_apply_func(x[t], f.t, p);
                }
            }
        }
    }

    return ret;
}

/* argument is series or list; value returned is list */

Node * dummify_func(Node * l, Node * r, Parser * p)
{
    Node * ret = aux_list_node(p);

    if (ret != null && starting(p)) {
        int * list = null;
        double oddval = NADBL;

        if (!null_node(r)) {
            if (r.t != GenSym.num) {
                p.err = e_types(r);
                return ret;
            } else {
                oddval = r.v.xval;
            }
        }

        if (l.t == GenSym.list) {
            list = gretl_list_copy(l.v.ivec);
        } else if (useries_node(l)) {
            list = gretl_list_new(1);
            list[1] = l.vnum;
        } else {
            gretl_errmsg_set(_("The first argument must be a named series "
                               ~ "in the current dataset").ptr);
            p.err = GretlError.invarg;
        }

        if (p.err) {
            /* don't do anything more */
        } else if (list == null) {
            p.err = GretlError.alloc;
        } else if (null_node(r)) {
            /* got just one argument */
            p.err = list_dumgenr(&list, p.dset, Gretlopt.OPT_F);
            ret.v.ivec = list;
        } else if (list[0] > 1) {
            gretl_errmsg_set("dummify(x, y): first argument should be a single series");
            free(list);
            p.err = GretlError.data;
        } else {
            p.err = dumgenr_with_oddval(&list, p.dset, oddval);
            ret.v.ivec = list;
        }
    }

    return ret;
}

/* argument is list; value returned is list */

Node * cdummify_func(Node * n, Parser * p)
{
    Node * ret = aux_list_node(p);

    if (ret != null && starting(p)) {
        int * list = null;

        if (n.t == GenSym.list) {
            list = gretl_list_copy(n.v.ivec);
        } else if (useries_node(n)) {
            list = gretl_list_new(1);
            list[1] = n.vnum;
        } else {
            p.err = e_types(n);
        }

        if (p.err) {
            /* don't do anything more */
        } else if (list == null) {
            p.err = GretlError.alloc;
        } else {
            p.err = auto_dummify_list(&list, p.dset);
            ret.v.ivec = list;
        }
    }

    return ret;
}

Node * get_info_on_series(Node * n, Parser * p)
{
    Node * ret = aux_bundle_node(p);

    if (ret != null && starting(p)) {
        int v = 0;

        if (useries_node(n)) {
            v = n.vnum;
        } else {
            v = node_get_int(n, p);
        }

        if (!p.err) {
            ret.v.b = series_info_bundle(p.dset, v, &p.err);
        }
    }

    return ret;
}

Node * list_stdize(Node * l, Node * r, Parser * p)
{
    Node * ret = aux_list_node(p);

    if (ret != null && starting(p)) {
        int * list = null;
        int dfc = 1;

        if (!null_node(r)) {
            dfc = node_get_int(r, p);
        }
        if (!p.err) {
            list = node_get_list(l, p);
        }
        if (!p.err) {
            Gretlopt opt;

            opt = dfc < 0 ? Gretlopt.OPT_C : dfc == 0 ? Gretlopt.OPT_N : Gretlopt.OPT_NONE;
            if (list[0] > 0) {
                p.err = list_stdgenr(list, p.dset, opt);
            }
            ret.v.ivec = list;
        }
    }

    return ret;
}

Node * series_stdize(Node * l, Node * r, Parser * p)
{
    Node * ret = aux_series_node(p);

    if (ret != null && starting(p)) {
        int dfc = 1;

        if (!null_node(r)) {
            dfc = node_get_int(r, p);
        }
        if (!p.err) {
            p.err = standardize_series(l.v.xvec, ret.v.xvec,
                                       dfc, p.dset);
        }
    }

    return ret;
}

/* middle argument is series or list; value returned is list in
   either case */

Node * list_make_lags(Node * l, Node * m, Node * r, Parser * p)
{
    Node * ret = aux_list_node(p);

    if (ret != null && starting(p)) {
        Gretlopt opt = Gretlopt.OPT_NONE;
        GretlMatrixRaw * v = null;
        int * list = null;
        int k = 0;

        /* ordering of the results? */
        if (node_get_bool(r, p, 0) > 0 && !p.err) {
            opt = Gretlopt.OPT_L; /* by lags */
        }

        if (!p.err) {
            /* scalar lag order or vector */
            if (l.t == GenSym.num) {
                k = node_get_int(l, p);
            } else {
                v = l.v.m;
            }
        }

        if (!p.err) {
            list = node_get_list(m, p);
        }

        if (!p.err) {
            if (list[0] > 0) {
                p.err = list_laggenr(&list, 1, k, v,
                                     p.dset, 0, opt);
            }
            ret.v.ivec = list;
        }
    }

    return ret;
}

Node * matrix_make_lags(Node * l, Node * m, Node * r, Parser * p)
{
    Node * ret = aux_matrix_node(p);

    if (ret != null && starting(p)) {
        Gretlopt opt = Gretlopt.OPT_NONE;
        GretlMatrixRaw * kvec = null;
        GretlMatrixRaw * src = null;

        if (node_get_bool(r, p, 0) > 0 && !p.err) {
            opt = Gretlopt.OPT_L; /* by lags */
        }

        if (!p.err) {
            /* scalar max lag order or vector */
            if (l.t == GenSym.num) {
                int i, k = node_get_int(l, p);

                if (!p.err && k <= 0) {
                    /* scalar k must be positive */
                    p.err = GretlError.invarg;
                } else if (!p.err) {
                    kvec = gretl_vector_alloc(k);
                    if (kvec == null) {
                        p.err = GretlError.alloc;
                    } else {
                        for (i=0; i<k; i++) {
                            kvec.val[i] = i+1;
                        }
                    }
                }
            } else {
                kvec = l.v.m;
            }
        }

        if (!p.err) {
            src = m.v.m;
        }

        if (!p.err) {
            ret.v.m = gretl_matrix_lag(src, kvec, opt, 0.0);
        }

        if (kvec != l.v.m) {
            gretl_matrix_free(kvec);
        }
    }

    return ret;
}

/* args are minlag, maxlag, MIDAS-list-to-be-lagged */

Node * hf_list_make_lags(Node * l, Node * m, Node * r, Parser * p)
{
    Node * ret = aux_list_node(p);

    if (ret != null && starting(p)) {
        int * list = null;
        int lmin, lmax = 0;
        int cfac = 0;

        lmin = node_get_int(l, p);
        if (!p.err) {
            lmax = node_get_int(m, p);
        }
        if (!p.err) {
            list = node_get_list(r, p);
        }

        if (!p.err) {
            /* compaction factor for high-frequency data */
            cfac = list[0];
            if (cfac < 2) {
                fprintf(stderr, "hflags: not a MIDAS list\n");
                p.err = GretlError.invarg;
                free(list);
            }
        }

        if (!p.err) {
            if (list[0] > 0) {
                p.err = list_laggenr(&list, lmin, lmax, null,
                                     p.dset, cfac, Gretlopt.OPT_L);
            }
            ret.v.ivec = list;
        }
    }

    return ret;
}

/* In case we get either of the two optional arguments to
   ldiff(), implement them via post-processing.
*/

void post_process_ldiff(Node * ret, Node * m, Node * r,
                        Parser * p)
{
    double d = NADBL;
    int fill, mul = 0;

    if (ret.t != GenSym.series && ret.t != GenSym.list) {
        /* "can't happen" */
        p.err = GretlError.data;
        return;
    }

    fill = node_get_bool(m, p, 0);
    if (!p.err && !null_node(r)) {
        d = node_get_scalar(r, p);
        mul = 1;
    }

    if (!p.err) {
        double * x = null;
        int * list = null;
        int imax = 1;
        int i, t;

        if (ret.t == GenSym.series) {
            x = ret.v.xvec;
        } else {
            list = ret.v.ivec;
            imax = list[0];
        }
        for (i=0; i<imax; i++) {
            if (list != null) {
                x = p.dset.Z[list[i+1]];
            }
            for (t=p.dset.t1; t<=p.dset.t2; t++) {
                if (na(x[t])) {
                    if (t > 0 && fill) {
                        x[t] = 0.0;
                    }
                } else if (mul) {
                    x[t] *= d;
                }
            }
        }
    }
}

bool ok_list_func(int f)
{
    return (f == GenFunc.f_log || f == GenFunc.f_diff ||
            f == GenFunc.f_ldiff || f == GenFunc.f_sdiff ||
            f == GenFunc.f_square || f == GenFunc.f_odev ||
            f == GenFunc.f_resample || f == GenFunc.f_dropcoll ||
            f == GenFunc.f_hfdiff || f == GenFunc.f_hfldiff);
}

/* The following handles functions that are "basically" for series,
   but which can also be applied to lists -- except for F_DROPCOLL,
   F_HFDIFF and F_HDLDIFF, which require a list argument.
*/

Node * apply_list_func(Node * n, Node * r, int f, Parser * p)
{
    Node * ret = aux_list_node(p);

    if (!ok_list_func(f)) {
        p.err = GretlError.types;
        return ret;
    }

    if (ret != null && starting(p)) {
        int * list = node_get_list(n, p);
        Gretlopt opt = Gretlopt.OPT_NONE;
        double parm = NADBL;
        GretlCmdIndex t = cast(GretlCmdIndex) 0;

        if (f == GenFunc.f_square) {
            if (r != null && node_is_true(r, p)) {
                opt = Gretlopt.OPT_O;
            }
        } else if (f == GenFunc.f_dropcoll || f == GenFunc.f_hfdiff ||
                   f == GenFunc.f_hfldiff) {
            /* handle optional parameter */
            if (!null_node(r)) {
                parm = node_get_scalar(r, p);
                if (p.err) {
                    return ret;
                }
            }
        }

        /* note: @list is modified by the library functions
           called below */

        if (list != null) {
            /* note: an empty list argument produces an
               empty list return
            */
            if (list[0] > 0) {
                switch (f) {
                case GenFunc.f_log:
                    p.err = list_loggenr(list, p.dset);
                    break;
                case GenFunc.f_diff:
                case GenFunc.f_ldiff:
                case GenFunc.f_sdiff:
                    if (f == GenFunc.f_diff) t = GretlCmdIndex.diff;
                    else if (f == GenFunc.f_ldiff) t = GretlCmdIndex.ldiff;
                    else if (f == GenFunc.f_sdiff) t = GretlCmdIndex.sdiff;
                    p.err = list_diffgenr(list, t, p.dset);
                    break;
                case GenFunc.f_square:
                    p.err = list_xpxgenr(&list, p.dset, opt);
                    break;
                case GenFunc.f_odev:
                    p.err = list_orthdev(list, p.dset);
                    break;
                case GenFunc.f_resample:
                    p.err = list_resample(list, p.dset);
                    break;
                case GenFunc.f_dropcoll:
                    p.err = list_dropcoll(list, parm, p.dset);
                    break;
                case GenFunc.f_hfdiff:
                case GenFunc.f_hfldiff:
                    t = (f == GenFunc.f_hfdiff)? GretlCmdIndex.diff : GretlCmdIndex.ldiff;
                    p.err = hf_list_diffgenr(list, t, parm, p.dset);
                    break;
                default:
                    break;
                }
            }
            ret.v.ivec = list;
        }
    }

    return ret;
}

Node * hf_list_node(Node * l, Node * m, Node * r, Parser * p)
{
    GretlMatrixRaw * v = l.v.m;
    int f_ratio = node_get_int(m, p);
    char * pfx = r.v.str;
    Node * ret = null;

    if (!p.err) {
        int n = gretl_vector_get_length(v);

        if (n == 0) {
            p.err = GretlError.nonconf;
        } else if (f_ratio < 3) {
            p.err = GretlError.invarg;
        } else if (*pfx == '\0' || !gretl_is_ascii(pfx) ||
                   strlen(pfx) > 24) {
            p.err = GretlError.invarg;
        } else {
            int T = sample_size(p.dset);

            if (n != f_ratio * T) {
                p.err = GretlError.invarg;
            }
        }
    }

    if (!p.err) {
        ret = aux_list_node(p);
        if (ret != null) {
            ret.v.ivec = vector_to_midas_list(v, f_ratio, pfx,
                                              p.dset, &p.err);
        }
    }

    return ret;
}

Node * dataset_list_node(Parser * p)
{
    Node * ret = null;

    if (gretl_function_depth() > 0) {
        gretl_errmsg_set("'dataset' is not recognized as a list within functions");
        p.err = GretlError.data;
    } else {
        ret = aux_list_node(p);
    }

    if (ret != null && starting(p)) {
        int * list = full_var_list(p.dset, null);

        if (list == null) {
            list = gretl_null_list();
        }
        if (list == null) {
            p.err = GretlError.data;
        }
        ret.v.ivec = list;
    }

    return ret;
}

Node * trend_node(Parser * p)
{
    Node * ret = null;

    if (starting(p)) {
        ret = aux_empty_series_node(p);
        if (!p.err) {
            p.err = gen_time(p.dset, 1, &ret.vnum);
            if (!p.err) {
                ret.v.xvec = p.dset.Z[ret.vnum];
                /* not TMP_NODE because we're borrowing a Z column */
                ret.flags = cast(ubyte) (ret.flags & ~NodeFlags.tmp);
            }
        }
    }

    return ret;
}

int object_get_size(Node * n, int context, Parser * p)
{
    int ret = 0;

    if (n.t == GenSym.num) {
        ret = 1;
    } else if (n.t == GenSym.mat) {
        ret = n.v.m.rows * n.v.m.cols;
    } else if (n.t == GenSym.array) {
        ret = gretl_array_get_length(n.v.a);
    } else if (n.t == GenSym.bundle) {
        ret = gretl_bundle_get_n_members(n.v.b);
    } else if (n.t == GenSym.list) {
        ret = n.v.ivec[0];
    } else if (n.t == GenSym.str) {
        if (context == DumId.end) {
            ret = cast(int) g_utf8_strlen(n.v.str, -1);
        } else {
            ret = cast(int) strlen(n.v.str);
        }
    } else {
        p.err = e_types(n);
    }

    return ret;
}

/* Cash out 'end' where it is legit, namely indicating the last
   element of a matrix row or column, or an array, list or string. We
   need to find the object to which 'end' is supposed to refer (@obj),
   and we also need to determine whether 'end' occurs in the first or
   second of the two possible indexation slots under the SLRAW node in
   its parentage. The second slot (idx = 2) is valid only as the column
   index for a matrix.
*/

int object_end_index(Node * t, Parser * p)
{
    Node * pa = t.parent;
    Node * obj = null;
    int idx = 0;
    int ret = -1;

    while (pa != null) {
        /* examine the parentage of @t */
        if (pa.t == GenSym.slraw) {
            if (find_in_tree(pa.R, t)) {
                idx = 2;
            } else if (pa.R != null && find_in_tree(pa.L, t)) {
                idx = 1;
            }
        } else if (pa.t == GenSym.osl) {
            obj = (pa.L.aux != null)? pa.L.aux : pa.L;
            break;
        }
        pa = pa.parent;
    }

    if (obj == null || (idx != 0 && obj.t != GenSym.mat)) {
        /* either we didn't find a referent, or it's a
           case of invalid double indexation
        */
        p.err = GretlError.invarg;
    } else if (obj.t == GenSym.mat) {
        if (idx == 2) {
            /* index must refer to column */
            ret = obj.v.m.cols;
        } else if (idx == 1) {
            /* index must refer to row */
            ret = obj.v.m.rows;
        } else {
            /* singleton index: may refer to row or column
               if the MAT is a vector; otherwise the error
               will be caught later
            */
            ret = gretl_vector_get_length(obj.v.m);
            if (ret == 0) {
                ret = obj.v.m.rows;
            }
        }
    } else if (obj.t == GenSym.array || obj.t == GenSym.list || obj.t == GenSym.str) {
        ret = object_get_size(obj, DumId.end, p);
    } else {
        p.err = GretlError.invarg;
    }

    if (p.err) {
        gretl_errmsg_set("'end' is not defined in this context");
    }

    return ret;
}

Node * array_last_node(Node * t, Parser * p)
{
    Node * ret = null;

    if (starting(p)) {
        ret = aux_scalar_node(p);
        if (!p.err) {
            ret.v.xval = object_end_index(t, p);
            if (p.err) {
                ret.v.xval = NADBL;
            }
        }
    }

    return ret;
}

Node * seasonals_node(Node * l, Node * r, Parser * p)
{
    Node * ret = null;

    if (!dataset_is_seasonal(p.dset) &&
        !dataset_is_seasonal_panel(p.dset)) {
        p.err = GretlError.pdwrong;
    } else {
        int ref_ = 0, center = 0;

        if (!null_node(l)) {
            ref_ = node_get_int(l, p);
        }
        if (!null_node(r)) {
            center = node_is_true(r, p);
        }
        if (!p.err) {
            ret = aux_list_node(p);
        }
        if (ret != null) {
            ret.v.ivec = seasonals_list(p.dset, ref_, center, &p.err);
        }
    }

    return ret;
}

Node * get_lag_list(Node * l, Node * r, Parser * p)
{
    Node * ret = null;

    if (starting(p)) {
        int * list = null;
        int * srclist = null;
        int i, imin = 1, imax = 1;
        int lv = 0;

        if (!useries_node(l) && l.t != GenSym.list) {
            /* we need a named series or a list on the left */
            p.err = e_types(l);
        } else if (r.t != GenSym.ivec && r.t != GenSym.num) {
            /* we need one or more integers on the right */
            p.err = e_types(r);
        }

        if (p.err) {
            return null;
        }

        if (l.t == GenSym.list) {
            srclist = l.v.ivec;
            imax = srclist[0];
        } else {
            lv = l.vnum;
        }

        if (imax == 0) {
            /* empty list on input -> empty on output */
            list = gretl_list_new(0);
            if (list == null) {
                p.err = GretlError.alloc;
            }
            goto loopdone;
        }

        for (i=imin; i<=imax && !p.err; i++) {
            if (srclist != null) {
                lv = srclist[i];
            }
            if (r.t == GenSym.ivec) {
                int fromlag = r.v.ivec[0];
                int tolag = r.v.ivec[1];

                if (list == null) {
                    list = laggenr_from_to(lv, fromlag, tolag,
                                           p.dset, &p.err);
                } else {
                    int * tmp;

                    tmp = laggenr_from_to(lv, fromlag, tolag,
                                          p.dset, &p.err);
                    if (!p.err) {
                        p.err = gretl_list_append_list(&list, tmp);
                        free(tmp);
                    }
                }
            } else {
                int lag = cast(int) -r.v.xval;

                lv = laggenr(lv, lag, p.dset);
                if (lv > 0) {
                    list = gretl_list_append_term(&list, lv);
                    if (list == null) {
                        p.err = GretlError.alloc;
                    }
                }
            }
        }

    loopdone:

        if (list != null) {
            ret = aux_list_node(p);
            if (ret != null) {
                ret.v.ivec = list;
            } else {
                free(list);
            }
        }
    } else {
        ret = aux_any_node(p);
    }

    return ret;
}

int * list_from_strings_array(GretlArray * a, Parser * p)
{
    GretlType type = gretl_array_get_type(a);
    int * list = null;

    if (type != GretlType.strings) {
        p.err = GretlError.types;
    } else {
        int i, vi, n = 0;
        char ** S = gretl_array_get_strings(a, &n);

        for (i=0; i<n && !p.err; i++) {
            vi = current_series_index(p.dset, S[i]);
            if (vi < 0) {
                gretl_errmsg_sprintf("'%s' is not a known series", S[i]);
                p.err = GretlError.unkvar;
            }
        }

        if (!p.err) {
            list = gretl_list_new(n);
            if (list == null) {
                p.err = GretlError.alloc;
            } else {
                for (i=0; i<n; i++) {
                    list[i+1] = current_series_index(p.dset, S[i]);
                }
            }
        }
    }

    return list;
}

/* Get an *int LIST from node @n: note that the list is always newly
   allocated, and so should be freed by the caller if it's just for
   temporary use.
*/

int * node_get_list(Node * n, Parser * p)
{
    int * list = null;
    int v = 0;

    if (null_node(n)) {
        if (p.flags & GenFlags.prnlist) {
            p.err = e_types(n);
        } else {
            list = gretl_null_list();
        }
    } else if (n.t == GenSym.list) {
        list = gretl_list_copy(n.v.ivec);
    } else if (n.t == GenSym.series || n.t == GenSym.num) {
        v = (n.t == GenSym.series)? n.vnum : node_get_int(n, p);
        if (!p.err) {
            if (v < 0 || v >= p.dset.v) {
                p.err = GretlError.unkvar;
            } else {
                list = gretl_list_new(1);
                if (list == null) {
                    p.err = GretlError.alloc;
                } else {
                    list[1] = v;
                }
            }
        }
    } else if (dataset_dum(n)) {
        list = full_var_list(p.dset, null);
    } else if (n.t == GenSym.mat) {
        list = gretl_list_from_vector(n.v.m, p.dset, &p.err);
    } else {
        p.err = e_types(n);
    }

    if (!p.err && list == null) {
        p.err = GretlError.alloc;
    } else if (p.err == GretlError.unkvar && v != 0) {
        gretl_errmsg_sprintf(_("Variable number %d is out of bounds").ptr, v);
    }

    return list;
}

Node * eval_lcat(Node * l, Node * r, Parser * p)
{
    Node * ret = aux_list_node(p);

    if (ret != null && starting(p)) {
        int * list1;
        int * list2 = null;

        list1 = node_get_list(l, p); /* note, copied */
        if (list1 != null) {
            list2 = node_get_list(r, p); /* copied */
        }
        if (list2 != null) {
            p.err = gretl_list_append_list(&list1, list2);
        }
        ret.v.ivec = list1;
        free(list2);
    }

    return ret;
}

Node * list_list_op(Node * l, Node * r, int f, Parser * p)
{
    Node * ret = aux_list_node(p);

    if (ret != null && starting(p)) {
        int * llist;
        int * rlist = null;
        int * list = null;

        llist = node_get_list(l, p);
        if (llist != null) {
            rlist = node_get_list(r, p);
        }
        if (rlist != null) {
            if (f == GenSym.b_and) {
                list = gretl_list_intersection(llist, rlist, &p.err);
            } else if (f == GenSym.b_or) {
                list = gretl_list_union(llist, rlist, &p.err);
            } else if (f == GenSym.b_sub) {
                list = gretl_list_drop(llist, rlist, &p.err);
            } else if (f == GenSym.b_pow) {
                list = gretl_list_product(llist, rlist, p.dset, &p.err);
            } else if (f == GenSym.b_add) {
                list = gretl_list_plus(llist, rlist, &p.err);
            }
        }
        ret.v.ivec = list;
        free(llist);
        free(rlist);
    }

    return ret;
}

/* Binary operator applied to two bundles: at present only '+'
   (for union) is supported.
*/

Node * bundle_op(Node * l, Node * r, int f, Parser * p)
{
    Node * ret = null;

    if (starting(p)) {
        GretlBundle * bl = l.v.b;
        GretlBundle * br = r.v.b;

        if (f == GenSym.b_add) {
            ret = aux_bundle_node(p);
        } else if (f == GenSym.b_eq || f == GenSym.b_neq) {
            ret = aux_scalar_node(p);
        } else {
            p.err = GretlError.types;
        }
        if (!p.err && f == GenSym.b_add) {
            ret.v.b = gretl_bundle_union(bl, br, &p.err);
        } else if (!p.err) {
            int eq = gretl_bundles_are_equal(bl, br);

            ret.v.xval = (f == GenSym.b_eq)? eq : !eq;
        }
    } else {
        ret = aux_any_node(p);
    }

    return ret;
}

/* Binary operator applied to two arrays: '+' (append),
   '||' (union) or '&&' (intersection). But the latter
   two are only for strings arrays.
*/

Node * array_op(Node * l, Node * r, int f, Parser * p)
{
    Node * ret = null;

    if (starting(p)) {
        GretlArray * al = l.v.a;
        GretlArray * ar = r.v.a;

        if (!p.err) {
            if (f == GenSym.b_eq || f == GenSym.b_neq) {
                ret = aux_scalar_node(p);
            } else if (f == GenSym.b_add || f == GenSym.b_or || f == GenSym.b_and) {
                ret = aux_array_node(p);
            } else {
                p.err = GretlError.types;
            }
        }
        if (!p.err) {
            if (f == GenSym.b_add) {
                ret.v.a = gretl_arrays_join(al, ar, &p.err);
            } else if (f == GenSym.b_or) {
                ret.v.a = gretl_arrays_union(al, ar, &p.err);
            } else if (f == GenSym.b_and) {
                ret.v.a = gretl_arrays_intersection(al, ar, &p.err);
            } else {
                int eq = gretl_arrays_are_equal(al, ar);

                ret.v.xval = (f == GenSym.b_eq)? eq : !eq;
            }
        }
    } else {
        ret = aux_any_node(p);
    }

    return ret;
}

Node * augment_array_node(Node * l, Node * r, Parser * p)
{
    Node * ret = aux_array_node(p);

    if (ret != null && starting(p)) {
        GretlType lt = gretl_array_get_content_type(l.v.a);
        GretlType rt = gretl_type_from_gen_type(r.t);

        if (rt == lt) {
            ret.v.a = gretl_array_copy(l.v.a, &p.err);
            if (!p.err) {
                p.err = gretl_array_append_object(ret.v.a, r.v.ptr, 1);
            }
        } else {
            p.err = GretlError.types;
        }
    }

    return ret;
}

Node * subtract_from_array_node(Node * l, Node * r, Parser * p)
{
    Node * ret = aux_array_node(p);

    if (ret != null && starting(p)) {
        if (gretl_array_get_type(l.v.a) == GretlType.strings &&
            r.t == GenSym.str) {
            ret.v.a = gretl_array_copy(l.v.a, &p.err);
            if (!p.err) {
                p.err = gretl_array_drop_string(ret.v.a, r.v.str);
            }
        }
    } else {
        p.err = GretlError.types;
    }

    return ret;
}

/* in case we switched the LHS and RHS in a boolean comparison */

int reversed_comp(int f)
{
    if (f == GenSym.b_gt) {
        return GenSym.b_lt;
    } else if (f == GenSym.b_lt) {
        return GenSym.b_gt;
    } else if (f == GenSym.b_gte) {
        return GenSym.b_lte;
    } else if (f == GenSym.b_lte) {
        return GenSym.b_gte;
    } else {
        return f;
    }
}

/* Boolean test of all vars in list against a scalar or series, for
   each observation in the sample, hence generating a series.
   The list will always be on the left-hand node; the 'reversed'
   flag is set if the list was originally on the right.
*/

Node * list_bool_comp(Node * l, Node * r, int f, int reversed,
                      Parser * p)
{
    Node * ret = aux_series_node(p);

    if (ret != null && starting(p)) {
        int * list = node_get_list(l, p);
        double * x = ret.v.xvec;
        double xit, targ = NADBL;
        double * tvec = null;
        int i, t;

        if (r.t == GenSym.num) {
            targ = r.v.xval;
        } else {
            tvec = r.v.xvec;
        }

        if (reversed) {
            f = reversed_comp(f);
        }

        if (list != null) {
            for (t=p.dset.t1; t<=p.dset.t2; t++) {
                if (tvec != null) {
                    targ = tvec[t];
                }
                if (na(targ)) {
                    x[t] = NADBL;
                    continue;
                }
                x[t] = 1.0; /* assume 'true' */
                for (i=1; i<=list[0]; i++) {
                    xit = p.dset.Z[list[i]][t];
                    if (na(xit)) {
                        x[t] = NADBL;
                        break;
                    } else if (f == GenSym.b_eq && xit != targ) {
                        x[t] = 0.0;
                    } else if (f == GenSym.b_neq && xit == targ) {
                        x[t] = 0.0;
                    } else if (f == GenSym.b_lt && xit >= targ) {
                        x[t] = 0.0;
                    } else if (f == GenSym.b_gt && xit <= targ) {
                        x[t] = 0.0;
                    } else if (f == GenSym.b_lte && xit > targ) {
                        x[t] = 0.0;
                    } else if (f == GenSym.b_gte && xit < targ) {
                        x[t] = 0.0;
                    }
                }
            }
            free(list);
        }
    }

    return ret;
}

/* Test for whether or not two lists are identical.  Note that
   using gretl_list_cmp() the order of the members matters.
   Perhaps the order shouldn't matter?
*/

Node * list_list_comp(Node * l, Node * r, int f, Parser * p)
{
    Node * ret = aux_scalar_node(p);

    if (ret != null && starting(p)) {
        int * llist = node_get_list(l, p);
        int * rlist = node_get_list(r, p);

        if (llist != null && rlist != null) {
            int d = gretl_list_cmp(llist, rlist);

            if (f == GenSym.b_neq) {
                ret.v.xval = d;
            } else if (f == GenSym.b_eq) {
                ret.v.xval = !d;
            } else {
                p.err = GretlError.types;
            }
        }
        free(llist);
        free(rlist);
    }

    return ret;
}

/* argument is list; value returned is series */

Node * list_to_series_func(Node * n, int f, Parser * p)
{
    Node * ret = aux_series_node(p);

    if (ret != null && starting(p)) {
        int partial_ok = 0;
        int * list = null;

        if (f == GenFunc.f_min || f == GenFunc.f_max) {
            /* Accept partial observations? This is debatable */
            partial_ok = 1;
        }
        list = node_get_list(n, p);
        if (list != null) {
            p.err = cross_sectional_stat(ret.v.xvec, list,
                                         p.dset, f, partial_ok);
            free(list);
        }
    }

    return ret;
}

/* arguments are series on left, list on right: we add all members
   of list to series, or subtract all members */

Node * series_list_calc(Node * l, Node * r, int f, Parser * p)
{
    Node * ret = aux_series_node(p);

    if (ret != null && starting(p)) {
        int * list = node_get_list(r, p);

        if (list != null) {
            double xt, xi;
            int i, t;

            for (t=p.dset.t1; t<=p.dset.t2; t++) {
                xt = l.v.xvec[t];
                if (!na(xt)) {
                    for (i=1; i<=list[0]; i++) {
                        xi = p.dset.Z[list[i]][t];
                        if (na(xi)) {
                            xt = NADBL;
                            break;
                        } else if (f == GenSym.b_add) {
                            xt += xi;
                        } else {
                            xt -= xi;
                        }
                    }
                }
                ret.v.xvec[t] = xt;
            }
            free(list);
        }
    }

    return ret;
}

int node_get_midas_method(Node * n, Parser * p)
{
    int ret = -1;

    if (n.t == GenSym.str) {
        if (!strcmp(n.v.str, "umidas")) {
            ret = MidasType.u;
        } else if (!strcmp(n.v.str, "nealmon")) {
            ret = MidasType.nealmon;
        } else if (!strcmp(n.v.str, "beta0")) {
            ret = MidasType.beta0;
        } else if (!strcmp(n.v.str, "betan")) {
            ret = MidasType.betan;
        } else if (!strcmp(n.v.str, "almonp")) {
            ret = MidasType.almonp;
        } else if (!strcmp(n.v.str, "beta1")) {
            ret = MidasType.beta1;
        }
    } else {
        ret = node_get_int(n, p);
    }

    if (ret < MidasType.u || ret >= MidasType.max) {
        p.err = GretlError.invarg;
    }

    return ret;
}

Node * lincomb_func(Node * l, Node * m, Node * r, int f, Parser * p)
{
    Node * ret = aux_series_node(p);

    if (ret != null && starting(p)) {
        int * list = node_get_list(l, p);
        const(GretlMatrixRaw) * b = node_get_real_matrix(m, p, 0, 2);
        int k = 0;

        if (!p.err && (list == null || gretl_is_null_matrix(b))) {
            p.err = GretlError.data;
        }

        if (!p.err && f == GenFunc.f_mlincomb) {
            k = node_get_midas_method(r, p);
        }

        if (!p.err) {
            if (f == GenFunc.f_mlincomb) {
                p.err = midas_linear_combo(ret.v.xvec, list, b, k, p.dset);
            } else {
                p.err = list_linear_combo(ret.v.xvec, list, b, p.dset);
            }
        }

        free(list);
    }

    return ret;
}

Node * list_list_series_func(Node * l1, Node * l2, int f,
                             Node * o, Parser * p)
{
    Node * ret = aux_series_node(p);

    if (ret != null && starting(p)) {
        int partial_ok = node_get_bool(o, p, 0);
        int * llist = node_get_list(l1, p);
        int * rlist = node_get_list(l2, p);

        if (!p.err) {
            p.err = x_sectional_weighted_stat(ret.v.xvec, llist, rlist,
                                              p.dset, f, partial_ok);
        }
        free(llist);
        free(rlist);
    }

    return ret;
}

/* check for missing obs in a list of variables */

Node * list_ok_func(Node * n, int f, Parser * p)
{
    Node * ret = aux_series_node(p);

    if (ret != null && starting(p)) {
        int * list = n.v.ivec;
        int i, vi, t;
        double x;

        if (list[0] == 0) {
            return ret;
        }

        for (t=p.dset.t1; t<=p.dset.t2; t++) {
            x = (f == GenFunc.f_dataok)? 1 : 0;
            for (i=1; i<=list[0]; i++) {
                vi = list[i];
                if (na(p.dset.Z[vi][t])) {
                    x = (f == GenFunc.f_dataok)? 0 : 1;
                    break;
                }
            }
            ret.v.xvec[t] = x;
        }
    }

    return ret;
}

/* functions taking (up to) two scalars as arguments and
   returning a series result */

Node * series_fill_func(Node * l, Node * r, int f, Parser * p)
{
    Node * ret = aux_series_node(p);

    if (ret != null && starting(p)) {
        double x, y;

        x = null_node(l) ? NADBL : node_get_scalar(l, p);
        y = null_node(r) ? NADBL : node_get_scalar(r, p);

        switch (f) {
        case GenFunc.f_runiform:
            p.err = gretl_rand_uniform_minmax(ret.v.xvec,
                                              p.dset.t1,
                                              p.dset.t2,
                                              x, y);
            break;
        case GenFunc.f_rnormal:
            p.err = gretl_rand_normal_full(ret.v.xvec,
                                           p.dset.t1,
                                           p.dset.t2,
                                           x, y);
            break;
        default:
            break;
        }
    }

    return ret;
}

GretlMatrixRaw * fc_matrix_from_list(Node * n, int n1, Parser * p)
{
    GretlMatrixRaw * ret = null;

    if (sample_size(p.dset) != n1) {
        p.err = GretlError.nonconf;
    } else {
        ret = gretl_matrix_data_subset(n.v.ivec,
                                       p.dset,
                                       p.dset.t1,
                                       p.dset.t2,
                                       MMissingCode.ok,
                                       &p.err);
    }

    return ret;
}

Node * fcstats_node(Node * l, Node * m, Node * r, Parser * p)
{
    Node * ret = aux_matrix_node(p);

    if (starting(p)) {
        GretlMatrixRaw * Fmat = null;
        const(double) * x = null, y = null;
        Gretlopt Fopt = Gretlopt.OPT_D;
        int U2, free_Fmat = 0;
        int n = 0, n2 = 0;

        if (l.t == GenSym.series || m.t == GenSym.series || m.t == GenSym.list) {
            if (dataset_is_time_series(p.dset)) {
                Fopt |= Gretlopt.OPT_T;
            }
        }
        U2 = node_get_bool(r, p, (Fopt & Gretlopt.OPT_T) ? 1 : 0);
        if (U2) {
            Fopt |= Gretlopt.OPT_T;
        } else {
            Fopt = cast(Gretlopt)(Fopt & ~Gretlopt.OPT_T);
        }

        if (!p.err) {
            if (l.t == GenSym.series) {
                n = sample_size(p.dset);
                x = l.v.xvec + p.dset.t1;
            } else {
                n = gretl_vector_get_length(l.v.m);
                if (n == 0) {
                    p.err = GretlError.types;
                } else {
                    x = l.v.m.val;
                }
            }
        }

        if (!p.err) {
            if (m.t == GenSym.series) {
                n2 = sample_size(p.dset);
                if (n2 != n) {
                    p.err = GretlError.nonconf;
                } else {
                    y = m.v.xvec + p.dset.t1;
                }
            } else if (m.t == GenSym.mat) {
                n2 = gretl_vector_get_length(m.v.m);
                if (n2 != n) {
                    if (m.v.m.rows == n) {
                        Fmat = m.v.m;
                    } else {
                        p.err = GretlError.nonconf;
                    }
                } else {
                    y = m.v.m.val;
                }
            } else {
                Fmat = fc_matrix_from_list(m, n, p);
                if (!p.err) {
                    free_Fmat = 1;
                }
            }
        }

        if (!p.err) {
            if (Fmat != null) {
                ret.v.m = matrix_fc_stats(x, Fmat, Fopt, &p.err);
            } else {
                ret.v.m = forecast_stats(x, y, 0, n-1, null, Fopt, &p.err);
            }
            if (free_Fmat) {
                gretl_matrix_free(Fmat);
            }
        }
    }

    return ret;
}

/* Functions taking two series as arguments and returning a scalar
   or matrix result. We also accept as arguments two matrices if
   they are vectors of the same length. In the case of F_NAALEN
   and F_KMEIER we can accept input with no @r node argument
   (meaning no censoring).
*/

Node * series_2_func(Node * l, Node * r, int f, Parser * p)
{
    Node * ret = null;

    if (starting(p)) {
        const(double) * x = null, y = null;
        int n = 0, n2 = 0;

        if (!p.err) {
            if (l.t == GenSym.series) {
                n = sample_size(p.dset);
                x = l.v.xvec + p.dset.t1;
            } else if (l.t == GenSym.num) {
                n = 1;
                x = &l.v.xval;
            } else if (l.t == GenSym.mat) {
                n = gretl_vector_get_length(l.v.m);
                if (n == 0) {
                    p.err = GretlError.types;
                } else {
                    x = l.v.m.val;
                }
            } else {
                p.err = GretlError.invarg;
            }
        }

        if (!p.err) {
            if (null_node(r)) {
                ; /* OK for duration funcs only */
            } else if (r.t == GenSym.series) {
                /* series on right */
                n2 = sample_size(p.dset);
                if (n2 != n) {
                    p.err = GretlError.nonconf;
                } else {
                    y = r.v.xvec + p.dset.t1;
                }
            } else if (r.t == GenSym.num) {
                n2 = 1;
                if (n2 != n) {
                    p.err = GretlError.nonconf;
                } else {
                    y = &r.v.xval;
                }
            } else if (r.t == GenSym.mat) {
                /* matrix on right */
                n2 = gretl_vector_get_length(r.v.m);
                if (n2 != n) {
                    p.err = GretlError.nonconf;
                } else {
                    y = r.v.m.val;
                }
            } else {
                p.err = e_types(r);
            }
        }

        if (p.err) {
            return ret;
        } else if (f == GenFunc.f_naalen || f == GenFunc.f_kmeier) {
            ret = aux_matrix_node(p);
        } else {
            ret = aux_scalar_node(p);
        }
        if (ret == null) {
            return null;
        }

        switch (f) {
        case GenFunc.f_cor:
            ret.v.xval = gretl_corr(0, n-1, x, y, null);
            break;
        case GenFunc.f_cov:
            ret.v.xval = gretl_covar(0, n-1, x, y, null);
            break;
        case GenFunc.f_naalen:
            ret.v.m = duration_func(x, y, 0, n-1, Gretlopt.OPT_NONE, &p.err);
            break;
        case GenFunc.f_kmeier:
            ret.v.m = duration_func(x, y, 0, n-1, Gretlopt.OPT_K, &p.err);
            break;
        case GenFunc.f_corresp:
            ret.v.xval = correspondence(x, y, n, &p.err);
            break;
        default:
            break;
        }
    } else {
        ret = aux_any_node(p);
    }

    return ret;
}

Gretlopt get_npcorr_option(Node * n, Parser * p)
{
    Gretlopt opt = Gretlopt.OPT_NONE;

    if (null_node(n)) {
        ; /* OK */
    } else {
        /* screened already: must be string */
        const(char) * s = n.v.str;

        if (!strcmp(s, "kendall")) {
            opt = Gretlopt.OPT_K;
        } else if (!strcmp(s, "spearman")) {
            opt = Gretlopt.OPT_S;
        } else {
            p.err = GretlError.invarg;
        }
    }

    return opt;
}

Node * npcorr_node(Node * l, Node * m, Node * r, Parser * p)
{
    Node * ret = aux_matrix_node(p);

    if (ret != null && starting(p)) {
        const(double) * x = null, y = null;
        Gretlopt opt = Gretlopt.OPT_NONE;
        int n1 = 0, n2 = 0;

        if (l.t == GenSym.series) {
            x = l.v.xvec + p.dset.t1;
            n1 = sample_size(p.dset);
        } else {
            n1 = gretl_vector_get_length(l.v.m);
            if (n1 == 0) {
                p.err = GretlError.invarg;
            } else {
                x = l.v.m.val;
            }
        }

        if (!p.err && m.t == GenSym.series) {
            y = m.v.xvec + p.dset.t1;
            n2 = sample_size(p.dset);
        } else if (!p.err) {
            n2 = gretl_vector_get_length(m.v.m);
            if (n2 == 0) {
                p.err = GretlError.invarg;
            } else {
                y = m.v.m.val;
            }
        }

        if (!p.err && n1 != n2) {
            p.err = GretlError.nonconf;
        } else if (!p.err) {
            opt = get_npcorr_option(r, p);
        }

        if (!p.err) {
            if (opt & Gretlopt.OPT_S) {
                ret.v.m = spearman_rho_func(x, y, n1, &p.err);
            } else {
                ret.v.m = kendall_tau_func(x, y, n1, &p.err);
            }
        }
    }

    return ret;
}

/* takes two series or two vectors as arguments */

Node * mxtab_func(Node * l, Node * r, Parser * p)
{
    Node * ret = aux_matrix_node(p);

    if (ret != null && starting(p)) {
        const(double) * x = null;
        const(double) * y = null;
        int nx = 0;
        int ny = 0;

        if (l.t == GenSym.mat && r.t == GenSym.mat) {
            nx = gretl_vector_get_length(l.v.m);
            ny = gretl_vector_get_length(r.v.m);
            x = l.v.m.val;
            y = r.v.m.val;
        } else if (l.t == GenSym.series && r.t == GenSym.series) {
            nx = ny = sample_size(p.dset);
            x = l.v.xvec + p.dset.t1;
            y = r.v.xvec + p.dset.t1;
        } else {
            p.err = GretlError.types;
        }
        if (!p.err && (nx == 0 || ny != nx)) {
            p.err = GretlError.nonconf;
        }
        if (!p.err) {
            ret.v.m = gretl_matrix_xtab(x, y, nx, &p.err);
        }
    }

    return ret;
}

/* remedial code to support <bundle>.<key> as passed to exists() */

GretlType remedy_old_exists(const(char) * s)
{
    GretlType ret = GretlType.none;
    char ** S;
    int ns = 0;

    S = gretl_string_split(s, &ns, ".");

    if (ns == 2) {
        GretlBundle * b = get_bundle_by_name(S[0]);

        if (b != null) {
            gretl_bundle_get_data(b, S[1], &ret, null);
        }
    }

    strings_array_free(S, ns);

    return ret;
}

Node * object_status(Node * n, Node * func, Parser * p)
{
    Node * ret = aux_scalar_node(p);
    int f = func.t;

    if (ret != null && starting(p)) {
        const(char) * s = n.v.str;

        ret.v.xval = NADBL;

        if (f == GenFunc.f_exists) {
            GretlType type = user_var_get_type_by_name(s);

            if (type == 0 && gretl_is_series(s, p.dset)) {
                type = GretlType.series;
            } else if (type == 0 && strchr(s, '.')) {
                type = remedy_old_exists(s);
            }
            ret.v.xval = gretl_type_get_order(type);
        } else if (f == GenFunc.f_isdiscr) {
            int v = current_series_index(p.dset, s);

            if (v >= 0) {
                ret.v.xval = series_is_discrete(p.dset, v);
            }
        } else if (f == GenFunc.f_obsnum) {
            int t = get_observation_number(s, p.dset);

            if (t > 0) {
                ret.v.xval = t;
            }
        } else if (f == GenFunc.f_strlen) {
            /* ret->v.xval = strlen(s); */
            ret.v.xval = g_utf8_strlen(s, -1);
        } else if (f == GenFunc.f_nlines) {
            ret.v.xval = count_lines(s);
        } else if (f == GenFunc.f_remove) {
            gretl_maybe_switch_dir(s);
            ret.v.xval = gretl_remove(s);
        }
    }

    return ret;
}

Node * multi_str_node(Node * l, int f, Parser * p)
{
    Node * ret = null;

    if (l.t == GenSym.series) {
        if (!strvals_node(l, p)) {
            p.err = GretlError.types;
        } else {
            ret = aux_series_node(p);
        }
    } else if (l.t == GenSym.array) {
        if (gretl_array_get_type(l.v.a) != GretlType.strings) {
            p.err = e_types(l);
        } else {
            ret = aux_matrix_node(p);
        }
    } else {
        p.err = e_types(l);
    }

    if (!p.err && l.t == GenSym.series) {
        SeriesTable * st;
        const(char) * s;
        int t;

        st = series_get_string_table(p.dset, l.vnum);
        for (t=p.dset.t1; t<=p.dset.t2; t++) {
            s = series_table_get_string(st, l.v.xvec[t]);
            ret.v.xvec[t] = (s == null)? NADBL : g_utf8_strlen(s, -1);
        }
    } else if (!p.err) {
        GretlMatrixRaw * m = null;
        char ** S;
        int i, ns = 0;

        S = gretl_array_get_strings(l.v.a, &ns);
        m = (ns == 0)? gretl_null_matrix_new() :
            gretl_matrix_alloc(ns, 1);
        if (m == null) {
            p.err = GretlError.alloc;
        } else {
            for (i=0; i<ns; i++) {
                m.val[i] = g_utf8_strlen(S[i], -1);
            }
            ret.v.m = m;
        }
    }

    return ret;
}

Node * generic_typeof_node(Node * n, Node * func, Parser * p)
{
    Node * ret = null;

    if (func.t == GenFunc.f_typename) {
        ret = aux_string_node(p);
    } else {
        ret = aux_scalar_node(p);
    }

    if (ret != null) {
        GretlType t = gretl_type_from_gen_type(n.t);

        if (func.t == GenFunc.f_typename) {
            if (t == GretlType.array) {
                t = gretl_array_get_type(n.v.a);
            }
            ret.v.str = gretl_strdup(gretl_type_get_name(t));
        } else if (func.t == GenFunc.f_iscmplx) {
            if (t == GretlType.matrix) {
                ret.v.xval = matrix_is_complex(n.v.m);
            } else {
                ret.v.xval = NADBL;
            }
        } else if (alias_reversed(func)) {
            /* handle the "isnull" alias for exists() */
            ret.v.xval = (t == 0);
        } else {
            /* func->t == F_TYPEOF */
            ret.v.xval = gretl_type_get_order(t);
        }
    }

    return ret;
}

/* return scalar node holding the number of elements in
   the object associated with node @n
*/

Node * n_elements_node(Node * n, Parser * p)
{
    Node * ret = aux_scalar_node(p);

    if (ret != null && starting(p)) {
        ret.v.xval = object_get_size(n, 0, p);
    }

    return ret;
}

int look_up_vname(const(char) * s, const(DatasetRaw) * dset)
{
    int i;

    for (i=0; i<dset.v; i++) {
        if (!strcmp(s, dset.varname[i])) {
            return i;
        }
    }

    return -1;
}

/* return scalar node holding the position of the series
   associated with node @r in the list associated with node
   @l, or zero if the series is not present in the list
*/

Node * in_list_node(Node * l, Node * r, Parser * p)
{
    Node * ret = aux_scalar_node(p);

    if (p.err == 0 && (p.dset == null || p.dset.v == 0)) {
        p.err = GretlError.nodata;
    }

    if (ret != null && starting(p)) {
        int * list = node_get_list(l, p);

        if (list != null) {
            int k = -1;

            if (useries_node(r)) {
                k = r.vnum;
            } else if (r.t == GenSym.num) {
                if (r.v.xval >= 0 && r.v.xval < p.dset.v) {
                    k = cast(int) r.v.xval;
                }
            } else if (r.t == GenSym.str) {
                k = look_up_vname(r.v.str, p.dset);
                if (k < 0) {
                    ret.v.xval = 0;
                }
            } else {
                node_type_error(GenFunc.f_inlist, 2, GenSym.series, r, p);
            }
            if (k >= 0) {
                ret.v.xval = in_gretl_list(list, k);
            }
            free(list);
        }
    }

    return ret;
}

Node * list_info_node(Node * l, Node * r, Parser * p)
{
    Node * ret = aux_matrix_node(p);
    int k = 0;

    if (p.err == 0 && (p.dset == null || p.dset.v == 0)) {
        p.err = GretlError.nodata;
    }

    if (!p.err && !null_node(r)) {
        k = node_get_int(r, p);
    }

    if (!p.err) {
        const(int) * list = l.v.ivec;
        Gretlopt opt = Gretlopt.OPT_NONE;

        if (k & 1) {
            /* "condense" the result */
            opt |= Gretlopt.OPT_C;
        }
        if (k & 2) {
            /* rely on data rather than just labels */
            opt |= Gretlopt.OPT_B;
        }
        ret.v.m = list_info_matrix(list, p.dset, opt, &p.err);
    }

    return ret;
}

Node * argname_from_uvar(Node * n, Node * r, Parser * p)
{
    Node * ret = aux_string_node(p);

    if (ret != null && starting(p)) {
        const(char) * vname = null;

        if (!null_or_string(r)) {
            /* if @r is present it must hold a string */
            p.err = e_types(r);
            return ret;
        }

        if (n.t == GenSym.series) {
            vname = p.dset.varname[n.vnum];
        } else {
            vname = n.vname;
        }

        if (vname == null) {
            p.err = GretlError.data;
        } else {
            ret.v.str = gretl_func_get_arg_name(vname, &p.err);
        }

        if (ret.v.str == null || ret.v.str[0] == '\0') {
            if (!null_node(r)) {
                ret.v.str = gretl_strdup(r.v.str);
            }
        }
    }

    return ret;
}

Node * varnum_node(Node * n, Parser * p)
{
    Node * ret = aux_scalar_node(p);

    if (ret != null && starting(p)) {
        if (n.t == GenSym.str) {
            int v = current_series_index(p.dset, n.v.str);

            ret.v.xval = (v >= 0)? v : NADBL;
        } else {
            p.err = GretlError.data;
        }
    }

    return ret;
}

Node * int_to_string_func(Node * n, int f, Parser * p)
{
    Node * ret = null;

    if (starting(p)) {
        GretlMatrixRaw * v = null;
        int i = 0;

        if (scalar_node(n)) {
            i = node_get_int(n, p);
        } else if (n.t == GenSym.series && f == GenFunc.f_varname) {
            i = n.vnum;
        } else if (n.t == GenSym.mat && f == GenFunc.f_obslabel) {
            v = n.v.m;
        } else {
            node_type_error(f, 0, GenSym.num, n, p);
        }

        if (!p.err && v != null) {
            ret = aux_array_node(p);
        } else if (!p.err) {
            ret = aux_string_node(p);
        }

        if (!p.err) {
            if (f == GenFunc.f_obslabel && v != null) {
                ret.v.a = retrieve_date_strings(v, p.dset, &p.err);
            } else if (f == GenFunc.f_obslabel) {
                ret.v.str = retrieve_date_string(i, p.dset, &p.err);
            } else if (f == GenFunc.f_varname) {
                if (i >= 0 && i < p.dset.v) {
                    ret.v.str = gretl_strdup(p.dset.varname[i]);
                } else {
                    p.err = GretlError.invarg;
                }
            } else if (f == GenFunc.f_randstr) {
                ret.v.str = gretl_rand_hex_string(i, &p.err);
            } else {
                p.err = GretlError.data;
            }
        }

        if (!p.err && v == null && ret.v.str == null) {
            p.err = GretlError.alloc;
        }
    } else {
        ret = aux_any_node(p);
    }

    return ret;
}

Node * list_to_string_func(Node * n, int f, Parser * p)
{
    Node * ret;

    if (f == GenFunc.f_varnames) {
        ret = aux_array_node(p);
    } else {
        ret = aux_string_node(p);
    }

    if (ret != null && starting(p)) {
        int * list = node_get_list(n, p);

        if (p.err) {
            return ret;
        }

        if (f == GenFunc.f_varname) {
            ret.v.str = gretl_list_get_names(list, p.dset,
                                             &p.err);
        } else if (f == GenFunc.f_varnames) {
            char ** S = gretl_list_get_names_array(list, p.dset,
                                                   &p.err);
            int ns = list[0];

            if (!p.err) {
                ret.v.a = gretl_array_from_strings(S, ns, 0, &p.err);
            }
        } else {
            p.err = GretlError.data;
        }

        free(list);
    }

    return ret;
}

/* handles both getenv (string value of variable) and
   ngetenv (numerical value of variable)
*/

Node * do_getenv(Node * l, int f, Parser * p)
{
    Node * ret = (f == GenFunc.f_getenv)? aux_string_node(p) :
        aux_scalar_node(p);

    if (ret != null && starting(p)) {
        int defined = 0;
        char * estr;

        estr = gretl_getenv(l.v.str, &defined, &p.err);

        if (f == GenFunc.f_getenv) {
            ret.v.str = estr;
        } else {
            /* ngetenv */
            if (defined) {
                char * test = null;
                double x;

                errno = 0;
                x = strtod(estr, &test);
                if (*test == '\0' && errno == 0) {
                    ret.v.xval = x;
                }
            }
            free(estr);
        }
    }

    return ret;
}

/* HAVE_MPI is not defined in this port: nothing to write */

void write_mpi_errmsg(const(char) * funcname, const(char) * s)
{
    return;
}

Node * do_errorif(Node * l, Node * r, Parser * p)
{
    Node * ret = aux_scalar_node(p);
    int fd = gretl_function_depth();

    if (fd == 0 && !gretl_mpi_initialized()) {
        gretl_errmsg_sprintf("'%s': can only be used within a function".ptr,
                             "errorif".ptr);
        p.err = GretlError.data;
    } else {
        int cond = node_get_bool(l, p, 1);

        if (cond) {
            const(char) * funcname = null;

            if (fd > 0) {
                current_function_info(&funcname, null);
            }
            if (gretl_mpi_initialized()) {
                if (mpi_rank == 0) {
                    write_mpi_errmsg(funcname, r.v.str);
                    p.err = GretlError.funcerr;
                }
            } else {
                gretl_errmsg_sprintf(_("Error message from %s():\n %s").ptr,
                                     funcname, r.v.str);
                p.err = GretlError.funcerr;
            }
        }
    }

    if (ret != null) {
        ret.v.xval = 1;
    }

    return ret;
}

Node * do_assert(Node * l, Node * r, Parser * p)
{
    Node * ret = aux_scalar_node(p);
    int assert_val = libset_get_int(SetKey.gretl_assert);

    if (ret == null) {
        p.err = GretlError.alloc;
        return null;
    }

    if (l.v.xval != 0 && !na(l.v.xval)) {
        /* non-zero, non-missing: success */
        ret.v.xval = 1;
    } else if (assert_val == 0) {
        /* flag but ignore failure */
        ret.v.xval = 0;
    } else if (assert_val == 1) {
        /* warn about failure */
        pprintf(p.prn, _("Warning: assertion '%s' failed").ptr, r.v.str);
        pputc(p.prn, '\n');
        ret.v.xval = 0;
    } else if (assert_val == 2) {
        /* complain and throw an error on failure */
        p.err = 1;
        gretl_errmsg_sprintf(_("Assertion '%s' failed").ptr, r.v.str);
        ret.v.xval = l.v.xval;
    } else {
        /* abort execution altogether */
        fprintf(stderr, "Assertion '%s' failed\n".ptr, r.v.str);
        abort();
    }

    return ret;
}

Node * contains_node(Node * val, Node * set, Parser * p)
{
    Node * ret = null;

    if (starting(p)) {
        const(GretlMatrixRaw) * m = node_get_matrix(set, p, 0, 2);

        if (val.t == GenSym.num) {
            ret = aux_scalar_node(p);
        } else if (val.t == GenSym.series) {
            ret = aux_series_node(p);
        } else if (val.t == GenSym.mat) {
            ret = aux_matrix_node(p);
            if (!p.err) {
                int r = val.v.m.rows;
                int c = val.v.m.cols;

                ret.v.m = gretl_zero_matrix_new(r, c);
            }
        }
        if (!p.err) {
            int i, n = m.rows * m.cols;
            double x;

            if (val.t == GenSym.num) {
                x = val.v.xval;
                if (na(x)) {
                    ret.v.xval = NADBL;
                } else {
                    ret.v.xval = 0;
                    for (i=0; i<n; i++) {
                        if (m.val[i] == x) {
                            ret.v.xval = 1;
                            break;
                        }
                    }
                }
            } else if (val.t == GenSym.series) {
                int t;

                for (t=p.dset.t1; t<=p.dset.t2; t++) {
                    x = val.v.xvec[t];
                    if (na(x)) {
                        ret.v.xvec[t] = NADBL;
                    } else {
                        ret.v.xvec[t] = 0;
                        for (i=0; i<n; i++) {
                            if (m.val[i] == x) {
                                ret.v.xvec[t] = 1;
                                break;
                            }
                        }
                    }
                }
            } else if (val.t == GenSym.mat) {
                GretlMatrixRaw * v = val.v.m;
                int k, nv = v.rows * v.cols;

                for (k=0; k<nv; k++) {
                    x = v.val[k];
                    if (na(x)) {
                        ret.v.m.val[k] = NADBL;
                    } else {
                        ret.v.m.val[k] = 0;
                        for (i=0; i<n; i++) {
                            if (m.val[i] == x) {
                                ret.v.m.val[k] = 1;
                                break;
                            }
                        }
                    }
                }
            }
        }
    } else {
        ret = aux_any_node(p);
    }

    return ret;
}

Node * single_string_func(Node * n, Node * x, int f, Parser * p)
{
    Node * ret = aux_string_node(p);

    if (ret != null && starting(p)) {
        const(char) * s = n.v.str;

        if (f == GenFunc.f_argname) {
            char * deflt = null;

            if (!null_node(x)) {
                if (x.t == GenSym.str) {
                    deflt = x.v.str;
                } else {
                    p.err = e_types(x);
                }
            }
            s = (n.vname != null)? n.vname : s;
            ret.v.str = gretl_func_get_arg_name(s, &p.err);
            if (!p.err && ret.v.str[0] == '\0' && deflt != null) {
                ret.v.str = gretl_strdup(deflt);
            }
        } else if (f == GenFunc.f_backtick) {
            ret.v.str = gretl_backtick(s, &p.err);
        } else if (f == GenFunc.f_strstrip) {
            ret.v.str = gretl_strstrip_copy(s, &p.err);
        } else if (f == GenFunc.f_fixname) {
            int uscore = 0;

            if (!null_node(x)) {
                uscore = node_get_bool(x, p, 0);
            }
            ret.v.str = cast(char *) calloc(VNAMELEN, 1);
            gretl_normalize_varname(ret.v.str, s, uscore, 0);
        } else {
            p.err = GretlError.data;
        }
    }

    return ret;
}

Node * country_code_node(Node * n, Node * r, Parser * p)
{
    Node * ret = null;
    char * src = null;
    char * tmp = null;

    if (n.t == GenSym.str) {
        src = n.v.str;
        ret = aux_string_node(p);
    } else if (n.t == GenSym.array || n.t == GenSym.series || n.t == GenSym.mat) {
        ret = aux_array_node(p);
    } else if (n.t == GenSym.num) {
        int k = node_get_int(n, p);

        if (!p.err) {
            src = tmp = g_strdup_printf("%03d".ptr, k);
            ret = aux_string_node(p);
        }
    } else {
        p.err = GretlError.invarg;
    }

    if (!p.err) {
        int output = 0; /* default to automatic */

        if (!null_node(r)) {
            output = node_get_int(r, p);
        }
        if (!p.err && src != null) {
            alias CFunc = extern (C) char * function(const(char) *, int, Prn *, int *);
            CFunc cfunc = cast(CFunc) get_plugin_function("iso_country");

            if (cfunc == null) {
                p.err = GretlError.fopen;
            } else {
                ret.v.str = cfunc(src, output, p.prn, &p.err);
            }
        } else if (!p.err && n.t == GenSym.array) {
            alias CFunc = extern (C) GretlArray * function(GretlArray *, int, Prn *, int *);
            CFunc cfunc = cast(CFunc) get_plugin_function("iso_country_array");

            if (cfunc == null) {
                p.err = GretlError.fopen;
            } else {
                ret.v.a = cfunc(n.v.a, output, p.prn, &p.err);
            }
        } else if (!p.err) {
            alias CFunc = extern (C) GretlArray * function(const(double) *, int, int, Prn *, int *);
            CFunc cfunc = cast(CFunc) get_plugin_function("iso_country_series");

            if (cfunc == null) {
                p.err = GretlError.fopen;
            } else {
                const(double) * x;
                int nx;

                if (n.t == GenSym.series) {
                    nx = sample_size(p.dset);
                    x = n.v.xvec + p.dset.t1;
                } else {
                    /* matrix */
                    nx = gretl_vector_get_length(n.v.m);
                    if (nx == 0) {
                        p.err = GretlError.invarg;
                    } else {
                        x = n.v.m.val;
                    }
                }
                if (!p.err) {
                    ret.v.a = cfunc(x, nx, output, p.prn, &p.err);
                }
            }
        }
    }

    g_free(tmp);

    return ret;
}

Node * readfile_node(Node * l, Node * r, Parser * p)
{
    Node * ret = aux_string_node(p);

    if (ret != null && starting(p)) {
        const(char) * codeset = null;

        if (!null_node(r)) {
            if (r.t == GenSym.str) {
                codeset = r.v.str;
            } else {
                node_type_error(GenFunc.f_readfile, 2, GenSym.str, r, p);
            }
        }
        if (!p.err) {
            ret.v.str = retrieve_file_content(l.v.str, codeset, &p.err);
        }
    }

    return ret;
}

const(char) * strstr_insensitive(const(char) * sd1,
                                 const(char) * sd2,
                                 const(char) * sl)
{
    const(char) * ret = strstr(sd1, sd2);

    if (ret != null) {
        /* preserve case of @sl? */
        const(char) * chk = sl + (ret - sd1);

        if (g_utf8_validate(chk, -1, null)) {
            ret = chk;
        }
    }

    return ret;
}

Node * two_string_func(Node * l, Node * r, Node * x,
                       int f, Parser * p)
{
    Node * ret = null;

    if (starting(p)) {
        const(char) * sl = l.v.str;
        const(char) * sr = null;

        if (f == GenFunc.f_jsongetb) {
            /* more checks are done below */
            if (!null_node(x) && x.t != GenSym.array) {
                p.err = GretlError.types;
            }
        } else if (f == GenFunc.f_xmlget && r.t == GenSym.array) {
            if (gretl_array_get_type(r.v.a) != GretlType.strings) {
                p.err = GretlError.types;
            }
        } else if (r.t != GenSym.str) {
            p.err = e_types(r);
        }
        if (!p.err) {
            ret = (f == GenFunc.f_instring)? aux_scalar_node(p) :
                (f == GenFunc.f_jsongetb)? aux_bundle_node(p) :
                aux_string_node(p);
        }
        if (p.err) {
            return null;
        }

        if (r != null && r.t == GenSym.str) {
            sr = r.v.str;
        }

        if (f == GenFunc.f_strstr || f == GenFunc.f_instring) {
            int ignore_case = node_get_bool(x, p, 0);
            const(char) * sret = null;

            if (ignore_case) {
                char * sd1 = g_utf8_strdown(sl, -1);
                char * sd2 = g_utf8_strdown(sr, -1);

                sret = strstr_insensitive(sd1, sd2, sl);
                g_free(sd1); g_free(sd2);
            } else {
                sret = strstr(sl, sr);
            }
            if (f == GenFunc.f_instring) {
                ret.v.xval = sret != null;
            } else if (sret != null) {
                ret.v.str = gretl_strdup(sret);
            } else {
                ret.v.str = gretl_strdup("");
            }
        } else if (f == GenSym.b_hcat) {
            int n1 = cast(int) strlen(l.v.str);
            int n2 = cast(int) strlen(r.v.str);

            ret.v.str = cast(char *) malloc(n1 + n2 + 1);
            if (ret.v.str != null) {
                *ret.v.str = '\0';
                strcat(ret.v.str, l.v.str);
                strcat(ret.v.str, r.v.str);
            }
        } else if (f == GenFunc.f_jsonget) {
            alias JFunc = extern (C) char * function(const(char) *, const(char) *,
                                                     int *, int *);
            JFunc jfunc = null;
            UserVar * uv = null;

            if (!null_node(x)) {
                uv = ptr_node_get_uvar(x, GenSym.num, p);
            }
            if (!p.err) {
                jfunc = cast(JFunc) get_plugin_function("json_get_string");
                if (jfunc == null) {
                    p.err = GretlError.fopen;
                }
            }
            if (jfunc != null) {
                int nobj = 0;
                int * pnobj = uv == null ? null : &nobj;

                ret.v.str = jfunc(l.v.str, r.v.str, pnobj, &p.err);
                if (!p.err && uv != null) {
                    user_var_set_scalar_value(uv, cast(double) *pnobj);
                }
            }
        } else if (f == GenFunc.f_jsongetb) {
            alias JFunc = extern (C) GretlBundle * function(const(char) *, const(char) *,
                                                            GretlArray *, int *);
            const(char) * path = null_node(r) ? null: r.v.str;
            GretlArray * strvars = null_node(x) ? null : x.v.a;

            JFunc jfunc = cast(JFunc) get_plugin_function("json_get_bundle");
            if (jfunc == null) {
                p.err = GretlError.fopen;
            } else {
                ret.v.b = jfunc(l.v.str, path, strvars, &p.err);
            }
        } else if (f == GenFunc.f_xmlget) {
            alias XFunc = extern (C) char * function(const(char) *, void *, GretlType,
                                                     int *, int *);
            XFunc xfunc = null;
            UserVar * uv = null;

            if (!null_node(x)) {
                uv = ptr_node_get_uvar(x, GenSym.num, p);
            }
            if (!p.err) {
                xfunc = cast(XFunc) get_plugin_function("xml_get");
                if (xfunc == null) {
                    p.err = GretlError.fopen;
                }
            }
            if (xfunc != null) {
                int nobj = 0;
                int * pnobj = uv == null ? null : &nobj;

                if (r.t == GenSym.array) {
                    ret.v.str = xfunc(l.v.str, r.v.a, GretlType.array,
                                      pnobj, &p.err);
                } else {
                    ret.v.str = xfunc(l.v.str, r.v.str, GretlType.string_,
                                      pnobj, &p.err);
                }
                if (!p.err && uv != null) {
                    user_var_set_scalar_value(uv, cast(double) *pnobj);
                }
            }
        } else {
            p.err = GretlError.data;
        }

        if (!p.err && f != GenFunc.f_instring && ret.v.str == null) {
            p.err = GretlError.alloc;
        }
    } else {
        ret = aux_any_node(p);
    }

    return ret;
}

Node * strseq_node(Node * l, Node * r, int rstr, Parser * p)
{
    Node * ret = null;
    Node * mn = rstr ? l : r; /* the MAT, ARRAY or NUM node */
    Node * sn = rstr ? r : l; /* the string node */
    const(char) * s = sn.v.str;
    GretlArray * S = null;
    int dim = 0;

    if (mn.t == GenSym.array) {
        dim = gretl_array_get_length(mn.v.a);
    } else if (mn.t == GenSym.mat) {
        dim = gretl_vector_get_length(mn.v.m);
    } else if (mn.t == GenSym.num) {
        dim = 1;
    }
    if (dim == 0) {
        p.err = GretlError.invarg;
        return null;
    }

    if (mn.t == GenSym.num) {
        ret = aux_string_node(p);
    } else {
        ret = aux_array_node(p);
        if (ret != null) {
            S = gretl_array_new(GretlType.strings, dim, &p.err);
        }
    }

    if (!p.err) {
        char * tmp;
        char * ds;
        int i, d;

        for (i=0; i<dim && !p.err; i++) {
            ds = null;
            if (mn.t == GenSym.array) {
                d = 0;
                ds = g_strdup(cast(char *) gretl_array_get_data(mn.v.a, i));
            } else if (mn.t == GenSym.mat) {
                d = cast(int) floor(mn.v.m.val[i]);
                ds = g_strdup_printf("%d".ptr, d);
            } else {
                d = cast(int) floor(mn.v.xval);
                ds = g_strdup_printf("%d".ptr, d);
            }
            if (d < 0) {
                p.err = GretlError.invarg;
            } else {
                if (rstr) {
                    tmp = g_strdup_printf("%s%s".ptr, ds, s);
                } else {
                    tmp = g_strdup_printf("%s%s".ptr, s, ds);
                }
                if (S != null) {
                    gretl_array_set_string(S, i, tmp, 1);
                } else {
                    ret.v.str = gretl_strdup(tmp);
                }
                g_free(tmp);
            }
            g_free(ds);
        }
        if (ret.t == GenSym.array) {
            ret.v.a = S;
        }
    }

    return ret;
}

int is_strings_array_node(Node * n)
{
    return n.t == GenSym.array &&
        gretl_array_get_type(n.v.a) == GretlType.strings;
}

/* "horizontal" concatenation of two gretl objects: we try to be as
   flexible as reasonably possible.
*/

Node * horizontal_concat_node(Node * l, Node * r, Parser * p)
{
    Node * ret = null;

    if (l.t == GenSym.str && r.t == GenSym.str) {
        /* (a) stick one string onto another */
        ret = two_string_func(l, r, null, GenSym.b_hcat, p);
    } else if (ok_matrix_node(l) && ok_matrix_node(r)) {
        /* (b) stick two conformable matrices together */
        ret = matrix_matrix_calc(l, r, GenSym.b_hcat, p);
    } else if (l.t == GenSym.str && (ok_matrix_node(r) || is_strings_array_node(r))) {
        /* (c) "multi-media" concatenation, string on left */
        ret = strseq_node(l, r, 0, p);
    } else if ((ok_matrix_node(l) || is_strings_array_node(l)) && r.t == GenSym.str) {
        /* (a) the converse of (c) */
        ret = strseq_node(l, r, 1, p);
    } else if (is_strings_array_node(l) && is_strings_array_node(r)) {
        /* (e) concatenation of strings, element by element, if the
           two arrays are of the same length
        */
        ret = aux_array_node(p);
        if (ret != null) {
            ret.v.a = gretl_arrays_concat(l.v.a, r.v.a, &p.err);
        }
    } else {
        p.err = GretlError.types;
    }

    return ret;
}

int is_matcher(int f)
{
    return f == GenFunc.f_strsub || f == GenFunc.f_regsub;
}

/* helper function for variant_string_func(): multi-string case */

char ** get_transformed_strings(Node * l, Node * m, Node * r,
                                int f, Parser * p, int * ns)
{
    char ** S1 = null;
    char ** S0;
    char * match = null;
    char * repl = null;
    int ini = 0;
    int fin = 0;

    if (l.t == GenSym.array) {
        S0 = gretl_array_get_strings(l.v.a, ns);
    } else {
        S0 = series_get_string_vals(p.dset, l.vnum, ns, 0);
    }

    if (S0 == null) {
        p.err = GretlError.data;
    } else {
        S1 = strings_array_new(*ns);
        if (S1 == null) {
            p.err = GretlError.alloc;
        }
    }

    if (f == GenFunc.f_substr) {
        ini = node_get_int(m, p);
        fin = node_get_int(r, p);
    } else if (is_matcher(f)) {
        match = m.v.str;
        repl = r.v.str;
    }

    if (!p.err) {
        char * tmp = null;
        int i;

        /* FIXME compile and re-use regexps? */

        for (i=0; i<*ns && !p.err; i++) {
            if (f == GenFunc.f_substr) {
                S1[i] = gretl_substring(S0[i], ini, fin, &p.err);
            } else if (f == GenFunc.f_strsub) {
                S1[i] = gretl_literal_replace(S0[i], match, repl, &p.err);
            } else if (f == GenFunc.f_regsub) {
                tmp = gretl_regexp_replace(S0[i], match, repl, &p.err);
            } else if (f == GenFunc.f_tolower) {
                tmp = g_utf8_strdown(S0[i], -1);
            } else {
                tmp = g_utf8_strup(S0[i], -1);
            }
            if (tmp != null) {
                S1[i] = gretl_strdup(tmp);
                g_free(tmp);
            }
        }
    }

    return S1;
}

/* helper function for variant_string_func(): single string case */

char * get_transformed_string(const(char) * s,
                              Node * m, Node * r,
                              int f, Parser * p)
{
    char * match = is_matcher(f) ? m.v.str : null;
    char * repl =  is_matcher(f) ? r.v.str : null;
    char * ret = null;
    char * tmp = null;

    if (f == GenFunc.f_substr) {
        int ini = node_get_int(m, p);
        int fin = node_get_int(r, p);

        if (!p.err) {
            ret = gretl_substring(s, ini, fin, &p.err);
        }
    } else if (f == GenFunc.f_strsub) {
        ret = gretl_literal_replace(s, match, repl, &p.err);
    } else if (f == GenFunc.f_regsub) {
        tmp = gretl_regexp_replace(s, match, repl, &p.err);
    } else if (f == GenFunc.f_tolower) {
        tmp = g_utf8_strdown(s, -1);
    } else if (f == GenFunc.f_toupper) {
        tmp = g_utf8_strup(s, -1);
    }
    if (tmp != null) {
        ret = gretl_strdup(tmp);
        g_free(tmp);
    }

    return ret;
}

/* A function that can be applied to a single string, an array of strings,
   or a string-valued series. If @m and @r are non-NULL they have already
   been verified as of type STR.
*/

Node * variant_string_func(Node * l, Node * m, Node * r,
                           int f, Parser * p)
{
    Node * ret = null;

    if (l.t == GenSym.str) {
        ret = aux_string_node(p);
    } else if (l.t == GenSym.array) {
        if (gretl_array_get_type(l.v.a) != GretlType.strings) {
            p.err = GretlError.types;
        } else {
            ret = aux_array_node(p);
        }
    } else if (l.t == GenSym.series) {
        p.err = strvals_return_check(l.parent, p);
        if (!p.err) {
            ret = aux_series_node(p);
        }
    } else {
        p.err = GretlError.data;
    }

    if (ret == null) {
        return null;
    }

    if (l.t == GenSym.str) {
        /* a single string target */
        ret.v.str = get_transformed_string(l.v.str, m, r, f, p);
    } else {
        /* array or series target */
        char ** S = null;
        int ns = 0;

        S = get_transformed_strings(l, m, r, f, p, &ns);
        if (l.t == GenSym.array) {
            ret.v.a = gretl_array_from_strings(S, ns, 0, &p.err);
        } else {
            double * y = ret.v.xvec + p.dset.t1;
            double * x = l.v.xvec + p.dset.t1;
            int n = sample_size(p.dset);

            p.err = series_from_string_transform(y, x, n, S, ns,
                                                 &p.lh.stab);
        }
    }

    return ret;
}

char * escape_strsplit_sep(const(char) * s)
{
    char * ret = cast(char *) calloc(strlen(s) + 1, 1);
    int i = 0;

    while (*s) {
        if (*s == '\\') {
            if (*(s+1) == '\\') {
                ret[i++] = '\\';
                s++;
            } else if (*(s+1) == 'n') {
                ret[i++] = '\n';
                s++;
            } else if (*(s+1) == 'r') {
                ret[i++] = '\r';
                s++;
            } else if (*(s+1) == 't') {
                ret[i++] = '\t';
                s++;
            } else {
                ret[i++] = *s;
            }
        } else {
            ret[i++] = *s;
        }
        s++;
    }

    return ret;
}

Node * strsplit_node(int f, Node * l, Node * m, Node * r, Parser * p)
{
    Node * ret = null;

    if (starting(p)) {
        const(char) * s = l.v.str;
        const(char) * sep0 = null;
        char * sep = null;
        int k = 0;

        /* We'll accept the two trailing optional arguments,
           string separator and integer index, in either order.
        */
        if (m != null) {
            if (m.t == GenSym.str) {
                sep0 = m.v.str;
            } else if (m.t != GenSym.empty) {
                k = node_get_int(m, p);
                if (k < 1) {
                    p.err = GretlError.invarg;
                }
            }
        }
        if (!p.err && r != null) {
            if (r.t == GenSym.str) {
                if (sep0 == null) {
                    /* OK, didn't get @sep0 yet */
                    sep0 = r.v.str;
                } else {
                    p.err = GretlError.invarg;
                }
            } else if (r.t != GenSym.empty) {
                if (k == 0) {
                    /* OK, didn't get @k yet */
                    k = node_get_int(r, p);
                    if (k < 1) {
                        p.err = GretlError.invarg;
                    }
                } else {
                    p.err = GretlError.invarg;
                }
            }
        }

        if (!p.err) {
            ret = k > 0 ? aux_string_node(p) : aux_array_node(p);
        }

        if (!p.err) {
            if (sep0 == null) {
                /* default: split on white space */
                sep = gretl_strdup(" \t\r\n");
            } else if (strchr(sep0, '\\')) {
                sep = escape_strsplit_sep(sep0);
            } else {
                sep = gretl_strdup(sep0);
            }
        }

        if (!p.err) {
            char ** S = null;
            int ns = 0;

            S = gretl_string_split(s, &ns, sep);
            if (!p.err) {
                if (k > 0) {
                    ret.v.str = gretl_strdup(k > ns ? "" : S[k-1]);
                    strings_array_free(S, ns);
                } else {
                    ret.v.a = gretl_array_from_strings(S, ns, 0, &p.err);
                }
            }
        }

        free(sep);
    } else {
        ret = aux_any_node(p);
    }

    return ret;
}

Node * array_sort_node(Node * n, int f, Parser * p)
{
    Node * ret = null;

    if (gretl_array_get_type(n.v.a) != GretlType.strings) {
        p.err = e_types(n);
    } else {
        ret = aux_array_node(p);
        if (!p.err) {
            ret.v.a = gretl_strings_sort(n.v.a, f == GenFunc.f_dsort, &p.err);
        }
    }
    return ret;
}

Node * array_func_node(Node * l, Node * r, Node * b, int f, Parser * p)
{
    GretlArray * a = l.v.a;
    GretlType t = gretl_array_get_type(a);
    Node * ret = null;

    if (f == GenFunc.f_instrings) {
        if (t != GretlType.strings || r.t != GenSym.str) {
            p.err = GretlError.types;
        } else {
            int simple = node_get_bool(b, p, 0);
            const(char) * s = r.v.str;

            if (!p.err && simple) {
                ret = aux_scalar_node(p);
                if (!p.err) {
                    ret.v.xval = gretl_strings_array_includes(a, s);
                }
            } else if (!p.err) {
                ret = aux_matrix_node(p);
                if (!p.err) {
                    ret.v.m = gretl_strings_array_pos(a, s, &p.err);
                }
            }
        }
    } else if (f == GenFunc.f_asort) {
        if (r.t != GenSym.str) {
            p.err = GretlError.types;
        } else {
            ret = aux_scalar_node(p);
        }
        if (!p.err) {
            ret.v.xval = p.err = gretl_array_qsort(a, r.v.str,
                                                   p.dset, p.prn);
        }
    } else if (f == GenFunc.f_flatten) {
        if (t == GretlType.matrices) {
            int mode = 0;

            if (!null_node(r)) {
                mode = node_get_int(r, p);
                if (mode < 0 || mode > 2) {
                    p.err = GretlError.invarg;
                }
            }
            if (!p.err) {
                ret = aux_matrix_node(p);
            }
            if (!p.err) {
                ret.v.m = gretl_matrix_array_flatten(l.v.a, mode, &p.err);
            }
        } else if (t == GretlType.strings) {
            const(char) * sep = null;

            if (r.t == GenSym.str) {
                sep = r.v.str;
            } else {
                int space = node_get_bool(r, p, 0);

                if (!p.err) {
                    sep = space ? " " : "\n";
                }
            }
            if (!p.err) {
                ret = aux_string_node(p);
            }
            if (!p.err) {
                ret.v.str = gretl_strings_array_flatten(l.v.a, sep, &p.err);
            }
        }
    } else {
        p.err = GretlError.types;
    }

    return ret;
}

Node * errmsg_node(Node * l, Parser * p)
{
    Node * ret = aux_string_node(p);

    if (ret != null && starting(p)) {
        const(char) * src = null;

        if (null_node(l)) {
            src = gretl_errmsg_get();
        } else {
            int errval = node_get_int(l, p);

            if (errval < 0 || errval >= GretlError.max) {
                p.err = GretlError.data;
            } else {
                src = errmsg_get_with_default(errval);
            }
        }

        if (src != null) {
            ret.v.str = gretl_strdup(src);
            if (ret.v.str == null) {
                p.err = GretlError.alloc;
            }
        }
    }

    return ret;
}
