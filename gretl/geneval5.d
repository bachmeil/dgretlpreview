/* Port of the fifth part (lines 18375-23327) of lib/src/geneval.c: geneval5.c. See geneval0.d for the notes that
   apply to all five parts of the file. */

module gretl.geneval5;


import core.stdc.errno : errno;
import core.stdc.limits : INT_MAX;
import core.stdc.math : fabs, fmod, pow, floor, ceil, log, lgamma, nearbyint, exp, atan2, isnan;
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
import gretl.genmainh : RetrievalIndex, BundleDataIndex;
import core.stdc.string : strncat, strncmp;
import core.stdc.math : acos, asin, atan, cos, sin, tan, acosh, asinh, atanh, cosh, sinh, tanh, sqrt;
import core.stdc.time : tm, mktime;
version (Posix) import core.sys.posix.time : strptime; else import gretl.winport : strptime;
import core.stdc.stdio : sprintf, sscanf;
import gretl.geneval2 : node_get_list, node_get_midas_method;
import gretl.gretl_tdisaggh : TdisaggInfo;
import gretl.gretl_bundleh : BundleType;
import gretl.objstackh : GretlObjType;
import gretl.versionh : GRETL_VERSION;
import gretl.genmainh : ModelDataIndex, model_data_scalar, model_data_series, model_data_matrix, model_data_matrix_builder, model_data_list, model_data_string, model_data_array;
import core.stdc.stdlib : calloc;
import core.stdc.time : time, time_t;
version (Posix) import core.sys.posix.time : localtime_r; else import gretl.winport : localtime_r;
import core.stdc.ctype : isspace;
import core.stdc.string : strspn, strcspn;
import gretl.uservarh : UVFlags;
import gretl.dataio : check_identifier;
import gretl.dataset : assign_numeric_to_strvar, assign_strings_to_strvar, dataset_add_NA_series, is_string_valued, overwrite_err, series_from_strings, series_get_string_table, series_unset_orig_pd, set_dataset_is_changed;
import gretl.dataseth : dataset_is_panel, sample_size;
import gretl.geneval1 : BFGS_constrained_max, BFGS_maximize, access_series_node, array_str_calc, bin2dec_node, binperm_node, bkw_node, bundle_file_write, clear_uvnodes, cmatrix_xy_calc, compare_strings, dec2bin_node, deriv_free_node, e_types, eval_pdist, eval_urcpval, free_node, free_tree, fzero_node, gen_add_or_replace, gen_add_uvar, gen_edit_list, gen_get_lhs_var, gen_replace_lhs, gen_replace_scalar, gen_type_is_arrayable, gretl_type_from_gen_type, matrix_add_names, matrix_bool, matrix_file_write, matrix_get_col_or_row_name, matrix_imhof, matrix_matrix_calc, matrix_scalar_calc, matrix_scalar_func, matrix_series_calc, matrix_string_func, matrix_to_alt_node, matrix_to_matrix_func, matrix_to_scalar_func, matrix_transpose_node, matrix_vector_func, maybe_rescue_undef_node, mpi_barrier_node, mpi_transfer_node, no_data_error, node_get_bool, node_get_scalar, numeric_jacobian_or_hessian, ok_list_node, ok_list_node_plus, real_matrix_calc, retrieve_const, revive_undef_node, rmatrix_xy_calc, scalar_calc, series_calc, series_string_calc, string_mult, string_offset, strvals_calc, type_func_in_parentage, typestr, xy_calc;
import gretl.geneval2 : apply_list_func, apply_scalar_func, apply_series_func, argname_from_uvar, array_func_node, array_last_node, array_op, array_sort_node, augment_array_node, bundle_op, cdummify_func, complex_matrix_node, contains_node, country_code_node, dataset_list_node, do_assert, do_errorif, do_getenv, dummify_func, errmsg_node, eval_lcat, fcstats_node, flexible_2arg_node, generic_typeof_node, get_info_on_series, get_lag_list, glasso_node, hf_list_make_lags, hf_list_node, horizontal_concat_node, in_list_node, int_to_string_func, invpd_node, lincomb_func, list_bool_comp, list_info_node, list_list_comp, list_list_op, list_list_series_func, list_make_lags, list_member_by_name, list_ok_func, list_reverse_node, list_stdize, list_to_series_func, list_to_string_func, matrix_fill_func, matrix_isnan_node, matrix_make_lags, matrix_to_matrix2_func, multi_str_node, mvnormal_node, mxtab_func, n_elements_node, npcorr_node, object_slice_node, object_status, post_process_ldiff, process_OSL_address, process_subslice, read_object_func, readfile_node, scalar_isnan_node, seasonals_node, series_2_func, series_fill_func, series_list_calc, series_stdize, single_string_func, sleep_node, strsplit_node, subobject_node, subtract_from_array_node, trend_node, two_string_func, variant_string_func, varnum_node;
import gretl.geneval3 : apply_matrix_func, atof_node, bvar_get_bundle, curl_bundle_node, do_irr, do_panel_expand, do_panel_shrink, dollar_bundle_node, dot_assign_to_matrix, eval_3args_func, eval_Rfunc, eval_ufunc, function_depth_node, get_bundle_array, get_bundle_member, getline_node, isconst_or_dum_node, isodate_node, lpsolve_bundle_node, matrix_quantiles_node, pergm_node, plot_colors_node, pretend_matrix_scalar_func, series_lag, series_ljung_box, series_matrix_node, series_obs, series_polyfit, series_scalar_func, series_scalar_scalar_func, series_series_func, series_sort_by, set_array_value, set_bundle_value, set_list_value, set_matrix_chunk, set_series_obs_value, set_string_value, strftime_node, strptime_node, svm_driver_node, test_bundle_key, type_string_node, vector_sort, vector_values;
import gretl.geneval4 : between_minmax_node, bundle_pointer_arg0, cast_series_to_list, could_be_series_id, destroy_multi, dollar_str_node, dollar_var_node, ellipsis_list_node, empty_ok, eval_bessel_func, eval_epochday, eval_feval, eval_kalman_bundle_func, eval_nargs_func, eval_query, evaluate_multi_args, fevd_node, gen_array_node, gen_series_node, geoplot_node, get_series_stringvals, isoconv_node, kalman_data_node, kpss_crit_node, lhs_terminal_node, list_join_node, mat2list_node, matrix_def_node, model_var_node, node_is_false, node_is_true, node_reattach_data, node_type_error, object_def_node, print_func_node, replace_value, scalar_postfix_node, series_calc_nodes, sscanf_node, stringify_series, two_scalars_func, wildlist_node, within_minmax_node;
import gretl.genmain : current_series_index, gretl_reserved_word;
import gretl.gretl_errors : errmsg, gretl_errmsg_set, gretl_errmsg_sprintf, gretl_warnmsg_sprintf, set_gretl_warning;
import gretl.gretl_errorsh : GretlWarningCodes;
import gretl.gretl_list : gretl_list_from_matrix, gretl_list_print, gretl_list_to_compact_string, gretl_list_to_vector;
import gretl.gretl_matrix : gretl_matrix_copy, gretl_matrix_copy_data, gretl_matrix_free, gretl_matrix_get_t1, gretl_matrix_get_t2, gretl_matrix_set_t1, gretl_matrix_set_t2, gretl_null_matrix_new;
import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_get_buffer, pprintf, pputc, pputs;
import gretl.gretl_string_table : SeriesTable, series_table_destroy;
import gretl.gstrfuncs : g_strchug, g_utf8_offset_to_pointer, g_utf8_strlen, g_utf8_validate;
import gretl.libset : get_matrix_mask, get_matrix_mask_nobs, libset_get_bool;
import gretl.missing : series_adjust_sample;
import gretl.missingh : NADBL, na;
import gretl.strutils : gretl_namechar_spn, gretl_strdup, gretl_strndup, strings_array_free, tailstrip;
import gretl.uservarh : uvar_get_scalar_value, uvar_set_scalar_fast;
import gretl.matrix_extra : gretl_matrix_from_scalar, gretl_matrix_print_to_prn;
import gretl.gretl_cmatrix : gretl_cmatrix_print;
import gretl.subsample : dataset_is_subsampled;
import gretl.gretl_array : get_strings_array_from_series, gretl_array_append_object, gretl_array_copy, gretl_array_copy_into, gretl_array_drop_null, gretl_array_drop_string, gretl_array_get_content_type, gretl_array_get_length, gretl_array_get_strings, gretl_array_get_type, gretl_array_new, gretl_array_print, gretl_array_void_content, gretl_singleton_array;
import gretl.gretl_bundle : gretl_bundle_append, gretl_bundle_copy, gretl_bundle_get_member_type, gretl_bundle_new, gretl_bundle_print, gretl_bundle_void_content;
import gretl.uservar : get_list_by_name, get_user_var_by_name, gretl_scalar_add, gretl_scalar_add_mutable, gretl_scalar_convert_to_matrix, remember_list, user_var_add;
import gretl.gretl_typemap : gretl_type_get_singular;
import gretl.genfuncs : check_declarations;

/* unported dependencies: per-file declarations from the 2026c headers */
extern (C) {
    void context_error(int c, Parser * p, const(char) * func);
    Node * expr(Parser * s);
    int function_lookup(const(char) * s);
    int get_loop_renaming();
    Ufunc * get_user_function_by_name(const(char) * name);
    const(char) * getsymb(int t);
    int gretl_function_depth();
    void lex(Parser * s);
    int object_is_const(const(char) * name, int vnum);
    void undefined_symbol_error(const(char) * s, Parser * p);
    const(char) * bvarname(int t);
    const(char) * constname(int c);
    const(char) * dumname(int t);
    const(char) * dvarname(int t);
    const(char) * mvarname(int t);
}
Node * eval(Node * t, Parser * p)
{
    Node * l = null, m = null, r = null;
    Node * multi = null;
    Node * ret = null;
    int ldone = 0;
    int rdone = 0;

    if (t == null) {
        /* catch NULL node right away */
        return null;
    }


    /* handle terminals first */
    if (t.t == GenSym.mspec) {
        return t;
    } else if (t.t == GenSym.empty) {
        if (empty_ok(t, p)) {
            return t;
        } else {
            return null;
        }
    } else if (t.t >= GenSym.num && t.t <= GenSym.str) {
        if (exestart(p) && uvar_node(t)) {
            node_reattach_data(t, p);
        }
        return t;
    }

    if (t.t == GenSym.query) {
        /* needs special treatment, see eval_query() */
        l = eval(t.L, p);
        if (p.err) {
            goto bailout;
        } else {
            goto do_switch;
        }
    }

    /* handle multi-argument L or R subnodes */
    if (t.L != null && bnsym(t.L.t)) {
        t.L.parent = t;
        multi = evaluate_multi_args(t.L, p);
        ldone = 1;
    } else if (t.R != null && bnsym(t.R.t)) {
        t.R.parent = t;
        multi = evaluate_multi_args(t.R, p);
        rdone = 1;
    }

    if (!p.err && t.L != null && !ldone) {
        t.L.parent = t;
        if (type_func_in_parentage(t)) {
            p.flags |= GenFlags.objqry;
            l = eval(t.L, p);
            if (p.err) {
                l = revive_undef_node(l, p);
                t.L.aux = l;
            }
        } else {
        p.flags = cast(GenFlags)(p.flags & ~GenFlags.objqry);
            l = eval(t.L, p);
        }
        if (l == null && !p.err) {
            p.err = 1;
        }
    }

    if (!p.err && t.M != null) {
        t.M.parent = t;
        if (m_return(t.t)) {
            m = t.M;
        } else {
            m = eval(t.M, p);
            if (m == null && !p.err) {
                p.err = 1;
            }
        }
    }

    if (!p.err && t.R != null && !rdone) {
        t.R.parent = t;
        if (r_return(t.t)) {
            r = t.R;
        } else {
            if (t.t == GenSym.b_and || t.t == GenSym.b_or) {
                /* logical operators: avoid redundant evaluation */
                if (l != null && l.t == GenSym.num) {
                    if ((t.t == GenSym.b_and && node_is_false(l, p)) ||
                        (t.t == GenSym.b_or && node_is_true(l, p))) {
                        /* no need to evaluate @r */
                        ret = l;
                        if (is_aux_node(ret)) {
                            /* mark return node as proxy */
                            ret.flags |= NodeFlags.prx;
                        }
                        goto finish;
                    }
                }
            } else if (t.t == GenFunc.f_genseries && t.R.t == GenSym.ufun) {
                p.flags |= GenFlags.ufret;
            }
            if (r == null && !p.err) {
                r = eval(t.R, p);
                if (r == null && !p.err) {
                    p.err = 1;
                }
            }
        }
    }

    if (p.err) {
        goto bailout;
    }

 do_switch:

    /* establish convenience pointer */
    p.aux = t.aux;

    switch (t.t) {
    case GenSym.dbundle:
        if (t.vname != null) {
            /* pseudo-bundle: name of model */
            ret = t;
        } else {
            /* built-in bundle identifed by idnum */
            ret = dollar_bundle_node(t, p);
        }
        break;
    case GenSym.undef:
        if (type_func_in_parentage(t)) {
            ret = t;
        } else {
            ret = maybe_rescue_undef_node(t, p);
        }
        break;
    case GenSym.u_addr:
        if (!uvar_node(t.L) && t.L.t != GenSym.osl) {
            p.err = GretlError.data;
        } else if (exestart(p)) {
            node_reattach_data(t.L, p);
        }
        ret = t;
        break;
    case GenSym.wlist:
        ret = wildlist_node(t, p);
        break;
    case GenSym.dum:
        if (t.v.idnum == DumId.dataset) {
            ret = dataset_list_node(p);
        } else if (t.v.idnum == DumId.trend) {
            ret = trend_node(p);
        } else if (t.v.idnum == DumId.end) {
            ret = array_last_node(t, p);
        } else {
            /* otherwise treat as terminal */
            ret = t;
        }
        break;
    case GenSym.num_p:
    case GenSym.num_m:
        if (exestart(p)) {
            node_reattach_data(t, p);
        }
        ret = scalar_postfix_node(t, p);
        break;
    case GenSym.b_add:
    case GenSym.b_sub:
    case GenSym.b_mul:
    case GenSym.b_div:
    case GenSym.b_mod:
    case GenSym.b_pow:
    case GenSym.b_and:
    case GenSym.b_or:
    case GenSym.b_eq:
    case GenSym.b_neq:
    case GenSym.b_gt:
    case GenSym.b_lt:
    case GenSym.b_gte:
    case GenSym.b_lte:
        /* arithmetic and logical binary operators: be as
           flexible as possible with regard to argument types
        */
        if (t.t == GenSym.b_add && l.t == GenSym.str && r.t == GenSym.num) {
            ret = string_offset(l, r, p);
        } else if ((t.t == GenSym.b_eq || t.t == GenSym.b_neq) &&
                   l.t == GenSym.str && r.t == GenSym.str) {
            ret = compare_strings(l, r, t.t, p);
        } else if (l.t == GenSym.num && r.t == GenSym.num) {
            ret = scalar_calc(l, r, t.t, p);
        } else if (l.t == GenSym.bundle && r.t == GenSym.bundle) {
            ret = bundle_op(l, r, t.t, p);
        } else if (l.t == GenSym.array && r.t == GenSym.array) {
            ret = array_op(l, r, t.t, p);
        } else if (strvals_node(l, p) && strvals_node(r, p)) {
            ret = strvals_calc(l, r, t, p);
        } else if (series_calc_nodes(l, r)) {
            ret = series_calc(l, r, t.t, p);
        } else if (l.t == GenSym.mat && r.t == GenSym.mat) {
            if (bool_comp(t.t)) {
                ret = matrix_bool(l, r, t.t, p);
            } else {
                ret = matrix_matrix_calc(l, r, t.t, p);
            }
        } else if ((l.t == GenSym.mat && r.t == GenSym.num) ||
                   (l.t == GenSym.num && r.t == GenSym.mat)) {
            ret = matrix_scalar_calc(l, r, t.t, p);
        } else if ((l.t == GenSym.mat && r.t == GenSym.series) ||
                   (l.t == GenSym.series && r.t == GenSym.mat)) {
            ret = matrix_series_calc(l, r, t.t, p);
        } else if (t.t >= GenSym.b_eq && t.t <= GenSym.b_neq &&
                   ((l.t == GenSym.series && r.t == GenSym.str) ||
                    (l.t == GenSym.str && r.t == GenSym.series))) {
            ret = series_string_calc(l, r, t.t, p);
        } else if ((t.t == GenSym.b_and || t.t == GenSym.b_or || t.t == GenSym.b_add || t.t == GenSym.b_sub) &&
                   ok_list_node_plus(l) && ok_list_node_plus(r)) {
            ret = list_list_op(l, r, t.t, p);
        } else if (t.t == GenSym.b_pow && ok_list_node(l, p) && ok_list_node(r, p)) {
            ret = list_list_op(l, r, t.t, p);
        } else if (bool_comp(t.t)) {
            if (ok_list_node(l, p) && (r.t == GenSym.num || r.t == GenSym.series)) {
                ret = list_bool_comp(l, r, t.t, 0, p);
            } else if (ok_list_node(r, p) && (l.t == GenSym.num || l.t == GenSym.series)) {
                ret = list_bool_comp(r, l, t.t, 1, p);
            } else if (ok_list_node(l, p) && ok_list_node(r, p)) {
                ret = list_list_comp(r, l, t.t, p);
            } else {
                p.err = GretlError.types;
            }
        } else if ((t.t == GenSym.b_add || t.t == GenSym.b_sub) &&
                   l.t == GenSym.series && ok_list_node(r, p)) {
            ret = series_list_calc(l, r, t.t, p);
        } else if (t.t == GenSym.b_add && l.t == GenSym.array) {
            ret = augment_array_node(l, r, p);
        } else if (t.t == GenSym.b_sub && l.t == GenSym.array) {
            ret = subtract_from_array_node(l, r, p);
        } else if (t.t == GenSym.b_mul && (l.t == GenSym.str || r.t == GenSym.str)) {
            ret = string_mult(l, r, p);
        } else {
            p.err = GretlError.types;
        }
        break;
    case GenSym.b_trmul:
        /* matrix on left, otherwise be flexible */
        if (ok_matrix_node(l) && ok_matrix_node(r)) {
            ret = matrix_matrix_calc(l, r, t.t, p);
        } else if (l.t == GenSym.mat && r.t == GenSym.series) {
            ret = matrix_series_calc(l, r, t.t, p);
        } else if (l.t == GenSym.mat && r.t == GenSym.empty) {
            ret = matrix_transpose_node(l, p);
        } else if (l.t == GenSym.num && r.t == GenSym.empty) {
            ret = l;
        } else {
            p.err = GretlError.types;
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
        /* matrix-matrix or matrix-scalar binary operators:
           in addition we permit scalar-scalar to allow for
           the possibility that results that could be taken
           to be 1 x 1 matrix results have been registered
           internally as scalars.
        */
        if (ok_matrix_node(l) && ok_matrix_node(r)) {
            ret = matrix_matrix_calc(l, r, t.t, p);
        } else if ((ok_matrix_node(l) && r.t == GenSym.series) ||
                   (l.t == GenSym.series && ok_matrix_node(r))) {
            ret = matrix_series_calc(l, r, t.t, p);
        } else if ((t.t == GenSym.b_doteq || t.t == GenSym.b_dotneq) &&
                   l.t == GenSym.array && r.t == GenSym.str) {
            ret = array_str_calc(l, r, t.t, p);
        } else {
            node_type_error(t.t, (l.t == GenSym.mat)? 2 : 1,
                            GenSym.mat, (l.t == GenSym.mat)? r : l, p);
        }
        break;
    case GenSym.subsl:
    case GenSym.b_range:
        /* matrix sub-slice, x:y, or lag range, 'p to q' */
        ret = process_subslice(l, r, p);
        break;
    case GenSym.b_vcat:
    case GenFunc.f_qform:
    case GenFunc.f_hdprod:
    case GenFunc.f_cmult:
    case GenFunc.f_cdiv:
    case GenFunc.f_lsolve:
    case GenFunc.f_mrsel:
    case GenFunc.f_mcsel:
    case GenFunc.f_dsum:
    case GenSym.b_ldiv:
    case GenSym.b_kron:
    case GenFunc.f_conv2d:
        if (ok_matrix_node(l) && ok_matrix_node(r)) {
            /* matrix-only binary operators (but promote scalars) */
            ret = matrix_matrix_calc(l, r, t.t, p);
        } else if (ok_matrix_node(l) && null_node(r) && t.t == GenFunc.f_hdprod) {
            ret = matrix_to_matrix2_func(l, r, t.t, p);
        } else {
            node_type_error(t.t, (l.t == GenSym.mat)? 2 : 1,
                            GenSym.mat, (l.t == GenSym.mat)? r : l, p);
        }
        break;
    case GenSym.b_hcat:
        ret = horizontal_concat_node(l, r, p);
        break;
    case GenSym.b_ellip:
        /* list-making ellipsis */
        if ((useries_node(l) || could_be_series_id(l, p)) &&
            (useries_node(r) || could_be_series_id(r, p))) {
            ret = ellipsis_list_node(l, r, p);
        } else {
            p.err = GretlError.types;
        }
        break;
    case GenFunc.f_mat2list:
        /* parlay matrix into list of series */
        if (l.t != GenSym.mat) {
            node_type_error(t.t, 1, GenSym.mat, l, p);
        } else if (!null_node(r) && r.t != GenSym.str) {
            node_type_error(t.t, 2, GenSym.str, r, p);
        } else {
            ret = mat2list_node(l, r, p);
        }
        break;
    case GenSym.b_join:
        /* list join with separator */
        if (ok_list_node(l, p) && ok_list_node(r, p)) {
            ret = list_join_node(l, r, p);
        } else {
            p.err = GretlError.types;
        }
        break;
    case GenFunc.f_msortby:
    case GenFunc.f_cswitch:
        /* matrix on left, scalar on right */
        if (l.t == GenSym.mat && null_or_scalar(r)) {
            ret = matrix_scalar_func(l, r, t.t, p);
        } else if (l.t == GenSym.mat) {
            node_type_error(t.t, 2, GenSym.num, r, p);
        } else {
            node_type_error(t.t, 1, GenSym.mat, l, p);
        }
        break;
    case GenFunc.f_distance:
        /* matrix on left, plus optional terms */
        if (l.t == GenSym.mat) {
            ret = matrix_string_func(l, m, r, t.t, p);
        } else {
            node_type_error(t.t, 1, GenSym.mat, l, p);
        }
        break;
    case GenFunc.f_msplitby:
        /* matrix on left, vector, optional boolean */
        if (ok_matrix_node(l) && ok_matrix_node(m)) {
            ret = matrix_vector_func(l, m, r, t.t, p);
        } else {
            p.err = GretlError.types;
        }
        break;
    case GenFunc.f_llag:
        if (null_node(m)) {
            p.err = GretlError.args;
        } else if (ok_matrix_node(l) && m.t != GenSym.mat && ok_list_node(m, p)) {
            ret = list_make_lags(l, m, r, p);
        } else if (ok_matrix_node(l) && m.t == GenSym.mat) {
            ret = matrix_make_lags(l, m, r, p);
        } else {
            p.err = GretlError.types;
        }
        break;
    case GenFunc.f_stdize:
        if (!null_or_scalar(m)) {
            p.err = e_types(m);
        } else if (l.t == GenSym.series) {
            ret = series_stdize(l, m, p);
        } else if (l.t == GenSym.list) {
            ret = list_stdize(l, m, p);
        } else if (l.t == GenSym.mat) {
            ret = eval_3args_func(l, m, r, t.t, p);
        } else {
            p.err = GretlError.types;
        }
        break;
    case GenFunc.f_hflag:
        if (scalar_node(l) && scalar_node(m) && ok_list_node(r, p)) {
            ret = hf_list_make_lags(l, m, r, p);
        } else {
            p.err = GretlError.types;
        }
        break;
    case GenFunc.f_hflist:
        if (l.t == GenSym.mat && scalar_node(m) && r.t == GenSym.str) {
            ret = hf_list_node(l, m, r, p);
        } else {
            p.err = GretlError.types;
        }
        break;
    case GenFunc.f_hfdiff:
    case GenFunc.f_hfldiff:
        if (ok_list_node(l, p) && null_or_scalar(r)) {
            ret = apply_list_func(l, r, t.t, p);
        } else {
            p.err = GretlError.types;
        }
        break;
    case GenSym.u_neg:
    case GenSym.u_pos:
    case GenSym.u_not:
    case GenFunc.f_abs:
    case GenFunc.f_sgn:
    case GenFunc.f_toint:
    case GenFunc.f_ceil:
    case GenFunc.f_floor:
    case GenFunc.f_round:
    case GenFunc.f_sin:
    case GenFunc.f_cos:
    case GenFunc.f_tan:
    case GenFunc.f_asin:
    case GenFunc.f_acos:
    case GenFunc.f_atan:
    case GenFunc.f_sinh:
    case GenFunc.f_cosh:
    case GenFunc.f_tanh:
    case GenFunc.f_asinh:
    case GenFunc.f_acosh:
    case GenFunc.f_atanh:
    case GenFunc.f_log:
    case GenFunc.f_log10:
    case GenFunc.f_log2:
    case GenFunc.f_exp:
    case GenFunc.f_sqrt:
    case GenFunc.f_cnorm:
    case GenFunc.f_dnorm:
    case GenFunc.f_qnorm:
    case GenFunc.f_logistic:
    case GenFunc.f_gamma:
    case GenFunc.f_lngamma:
    case GenFunc.f_digamma:
    case GenFunc.f_trigamma:
    case GenFunc.f_invmills:
    case GenFunc.f_easter:
        /* functions taking one numeric argument, any type */
        if (l.t == GenSym.num) {
            ret = apply_scalar_func(l, t, p);
        } else if (l.t == GenSym.series) {
            if (cast_series_to_list(p, l, t.t)) {
                ret = apply_list_func(l, null, t.t, p);
            } else {
                ret = apply_series_func(l, t, p);
            }
        } else if (l.t == GenSym.mat) {
            ret = apply_matrix_func(l, t, p);
        } else if (ok_list_node(l, p) && t.t == GenFunc.f_log) {
            ret = apply_list_func(l, null, t.t, p);
        } else {
            p.err = e_types(l);
        }
        break;
    case GenFunc.f_real:
    case GenFunc.f_imag:
        if (l.t == GenSym.mat || l.t == GenSym.num) {
            ret = apply_matrix_func(l, t, p);
        } else {
            p.err = GretlError.types;
        }
        break;
    case GenFunc.f_carg:
    case GenFunc.f_conj:
    case GenFunc.f_cmod:
    case GenFunc.f_cquad:
        if (complex_node(l)) {
            ret = apply_matrix_func(l, t, p);
        } else {
            p.err = GretlError.types;
        }
        break;
    case GenFunc.f_atan2:
    case GenFunc.f_bincoeff:
    case GenFunc.f_lnmgamma:
        if ((l.t == GenSym.num || l.t == GenSym.mat || l.t == GenSym.series) &&
            (r.t == GenSym.num || r.t == GenSym.mat || r.t == GenSym.series)) {
            ret = flexible_2arg_node(l, r, t.t, 0, p);
        } else {
            p.err = GretlError.types;
        }
        break;
    case GenFunc.f_threshold:
        if ((l.t == GenSym.num || l.t == GenSym.mat || l.t == GenSym.series) &&
            (m.t == GenSym.num || m.t == GenSym.mat || m.t == GenSym.series)) {
        int hard = node_get_bool(r, p, 0);
        if (!p.err) {
        ret = flexible_2arg_node(l, m, t.t, hard, p);
        }
        } else {
            p.err = GretlError.types;
        }
        break;
    case GenFunc.f_dumify:
    case GenFunc.f_cdumify:
        /* series or list argument wanted */
        if (l.t == GenSym.series || l.t == GenSym.list) {
            if (t.t == GenFunc.f_dumify) {
                ret = dummify_func(l, r, p);
            } else {
                ret = cdummify_func(l, p);
            }
        } else {
            p.err = GretlError.invarg;
        }
        break;
    case GenFunc.f_getinfo:
        /* named series (or ID) argument wanted */
        if (useries_node(l) || l.t == GenSym.num) {
            ret = get_info_on_series(l, p);
        } else {
            p.err = e_types(l);
        }
        break;
    case GenFunc.f_seasonals:
        /* two optional args: int, bool */
        if (!null_or_scalar(l)) {
            node_type_error(t.t, 1, GenSym.num, l, p);
        } else if (!null_or_scalar(r)) {
            node_type_error(t.t, 2, GenSym.num, r, p);
        } else {
            ret = seasonals_node(l, r, p);
        }
        break;
    case GenFunc.f_misszero:
    case GenFunc.f_zeromiss:
        /* one series, scalar or matrix argument needed */
        if (l.t == GenSym.series) {
            ret = apply_series_func(l, t, p);
        } else if (l.t == GenSym.num) {
            ret = apply_scalar_func(l, t, p);
        } else if (l.t == GenSym.mat) {
            ret = apply_matrix_func(l, t, p);
        } else {
            node_type_error(t.t, 0, GenSym.series, l, p);
        }
        break;
    case GenFunc.f_missing:
    case GenFunc.f_dataok:
        /* series, scalar or list argument needed */
        if (l.t == GenSym.mat) {
            if (t.t == GenFunc.f_dataok) {
                ret = matrix_to_matrix_func(l, null, t.t, p);
            } else {
                ret = matrix_isnan_node(l, p);
            }
        } else if (l.t == GenSym.series) {
            ret = apply_series_func(l, t, p);
        } else if (l.t == GenSym.num) {
            ret = apply_scalar_func(l, t, p);
        } else if (l.t == GenSym.list) {
            ret = list_ok_func(l, t.t, p);
        } else {
            node_type_error(t.t, 0, GenSym.series, l, p);
        }
        break;
    case GenFunc.f_isnan:
        /* scalar or matrix */
        if (scalar_node(l)) {
            ret = scalar_isnan_node(l, p);
        } else if (l.t == GenSym.mat) {
            ret = matrix_isnan_node(l, p);
        } else {
            node_type_error(t.t, 0, GenSym.num, l, p);
        }
        break;
    case GenFunc.f_sleep:
        if (scalar_node(l)) {
            ret = sleep_node(l, t.t, p);
        } else {
            node_type_error(t.t, 0, GenSym.num, l, p);
        }
        break;
    case GenFunc.f_barrier:
        if (l.t == GenSym.empty) {
            ret = mpi_barrier_node(p);
        } else {
            node_type_error(t.t, 0, GenSym.empty, l, p);
        }
        break;
    case GenSym.lag:
        if (p.targ == GenSym.list) {
            ret = get_lag_list(l, r, p);
        } else if (l.t == GenSym.series && scalar_node(r)) {
            ret = series_lag(l, r, p);
        } else if (l.t != GenSym.series) {
            node_type_error(t.t, 1, GenSym.series, l, p);
        } else {
            node_type_error(t.t, 2, GenSym.num, r, p);
        }
        break;
    case GenFunc.f_ljungbox:
    case GenFunc.f_polyfit:
        /* series on left, scalar on right */
        if (l.t != GenSym.series) {
            node_type_error(t.t, 1, GenSym.series, l, p);
        } else if (!scalar_node(r)) {
            node_type_error(t.t, 2, GenSym.num, r, p);
        } else if (t.t == GenFunc.f_ljungbox) {
            ret = series_ljung_box(l, r, p);
        } else if (t.t == GenFunc.f_polyfit) {
            ret = series_polyfit(l, r, p);
        }
        break;
    case GenSym.obs:
        if (l.t != GenSym.series) {
            node_type_error(t.t, 1, GenSym.series, l, p);
        } else if (!scalar_node(r) && r.t != GenSym.str) {
            node_type_error(t.t, 2, GenSym.num, r, p);
        } else if (t.flags & NodeFlags.lht) {
            ret = lhs_terminal_node(t, l, r, p);
        } else {
            ret = series_obs(l, r, p);
        }
        break;
    case GenSym.osl:
        /* object plus subspec */
        if (t.flags & NodeFlags.lht) {
            ret = lhs_terminal_node(t, l, r, p);
        } else if (l.t == GenSym.u_addr) {
            ret = process_OSL_address(t, l, r, p);
        } else {
            ret = subobject_node(l, r, p);
        }
        break;
    case GenSym.slraw:
        /* unevaluated object slice spec */
        ret = object_slice_node(t, l, r, p);
        break;
    case GenSym.listvar:
        ret = list_member_by_name(l, r, p);
        break;
    case GenSym.bmemb:
    case GenFunc.f_inbundle:
        /* name of bundle plus string */
        if (l.t == GenSym.bundle && r.t == GenSym.str) {
            if (t.t == GenSym.bmemb) {
                if (t.flags & NodeFlags.lht) {
                    ret = lhs_terminal_node(t, l, r, p);
                } else {
                    ret = get_bundle_member(l, r, p);
                }
            } else {
                ret = test_bundle_key(l, r, p);
            }
        } else if (l.t == GenSym.bundle) {
            node_type_error(t.t, 2, GenSym.str, r, p);
        } else {
            node_type_error(t.t, 1, GenSym.bundle, l, p);
        }
        break;
    case GenFunc.f_getkeys:
    case GenFunc.hf_jbterms:
        if (l.t == GenSym.bundle) {
            ret = get_bundle_array(l, t.t, p);
        } else {
            node_type_error(t.t, 0, GenSym.bundle, l, p);
        }
        break;
    case GenSym.dbmemb:
        /* name of $-bundle plus string */
        if (l.t == GenSym.bundle && r.t == GenSym.str) {
            ret = get_bundle_member(l, r, p);
        } else if (l.t == GenSym.dbundle && r.t == GenSym.str) {
            ret = model_var_node(l, r, p);
        } else if (r.t != GenSym.str) {
            node_type_error(t.t, 2, GenSym.str, r, p);
        } else {
            node_type_error(t.t, 1, GenSym.bundle, l, p);
        }
        break;
    case GenFunc.hf_fdepth:
        ret = function_depth_node(l, p);
        break;
    case GenFunc.f_curl:
        ret = curl_bundle_node(l, p);
        break;
    case GenFunc.f_gpcolors:
    ret = plot_colors_node(l, p);
    break;
    case GenFunc.f_lpsolve:
        if (l.t == GenSym.bundle) {
            ret = lpsolve_bundle_node(l, p);
        } else {
            node_type_error(t.t, 0, GenSym.bundle, l, p);
        }
        break;
    case GenFunc.f_svm:
        ret = svm_driver_node(t, multi, p);
        break;
    case GenFunc.f_typestr:
        /* numerical type code to string */
        if (scalar_node(l)) {
            ret = type_string_node(l, p);
        } else {
            node_type_error(t.t, 0, GenSym.num, l, p);
        }
        break;
    case GenFunc.f_ldiff:
    case GenFunc.f_sdiff:
    case GenFunc.f_odev:
    case GenFunc.f_interpol:
        if (l.t == GenSym.series && cast_series_to_list(p, l, t.t)) {
            ret = apply_list_func(l, null, t.t, p);
        } else if (l.t == GenSym.series ||
                   (t.t != GenFunc.f_odev && t.t != GenFunc.f_interpol && l.t == GenSym.mat)) {
            ret = series_series_func(l, r, null, t.t, p);
        } else if (ok_list_node(l, p)) {
            ret = apply_list_func(l, null, t.t, p);
        } else {
            node_type_error(t.t, 0, GenSym.series, l, p);
        }
        break;
    case GenFunc.f_dropcoll:
        /* list argument is required on left, optional scalar
           on the right */
        if (l.t == GenSym.list) {
            ret = apply_list_func(l, r, t.t, p);
        } else {
            node_type_error(t.t, 0, GenSym.list, l, p);
        }
        break;
    case GenFunc.f_hpfilt:
    case GenFunc.f_fracdiff:
    case GenFunc.f_fraclag:
    case GenFunc.f_boxcox:
    case GenFunc.f_pnobs:
    case GenFunc.f_pmin:
    case GenFunc.f_pmax:
    case GenFunc.f_psum:
    case GenFunc.f_pmean:
    case GenFunc.f_pxsum:
    case GenFunc.f_pxnobs:
    case GenFunc.f_pxmean:
    case GenFunc.f_psd:
    case GenFunc.f_deseas:
    case GenFunc.f_tramolin:
        /* series argument needed */
        if (l.t == GenSym.series || l.t == GenSym.mat) {
            if (t.t == GenFunc.f_hpfilt) {
                ret = series_series_func(l, m, r, t.t, p);
            } else {
                ret = series_series_func(l, r, null, t.t, p);
            }
        } else {
            node_type_error(t.t, 0, GenSym.series, l, p);
        }
        break;
    case GenFunc.f_pshrink:
        if (l.t == GenSym.series) {
            int noskip = node_get_bool(r, p, 0);

            if (!p.err) {
                ret = do_panel_shrink(l, noskip, p);
            }
        } else {
            node_type_error(t.t, 0, GenSym.series, l, p);
        }
        break;
    case GenFunc.f_pexpand:
        if (l.t == GenSym.mat) {
            ret = do_panel_expand(l, r, p);
        } else {
            node_type_error(t.t, 0, GenSym.mat, l, p);
        }
        break;
    case GenFunc.f_cum:
    case GenFunc.f_diff:
    case GenFunc.f_resample:
    case GenFunc.f_ranking:
        /* series or matrix argument */
        if (l.t == GenSym.series) {
            if (cast_series_to_list(p, l, t.t)) {
                ret = apply_list_func(l, null, t.t, p);
            } else {
                ret = series_series_func(l, r, null, t.t, p);
            }
        } else if (l.t == GenSym.mat) {
            if (t.t == GenFunc.f_resample) {
                ret = eval_3args_func(l, m, r, t.t, p);
            } else {
                ret = matrix_to_matrix_func(l, r, t.t, p);
            }
        } else if (t.t == GenFunc.f_diff && ok_list_node(l, p)) {
            ret = apply_list_func(l, null, t.t, p);
        } else if (t.t == GenFunc.f_resample && ok_list_node(l, p) &&
                   null_node(r)) {
            ret = apply_list_func(l, null, t.t, p);
        } else {
            node_type_error(t.t, 0, GenSym.series, l, p);
        }
        break;
    case GenFunc.f_sort:
    case GenFunc.f_dsort:
        /* series or vector or string array argument needed */
        if (l.t == GenSym.series || l.t == GenSym.mat || l.t == GenSym.num) {
            ret = vector_sort(l, t.t, p);
        } else if (l.t == GenSym.array) {
            ret = array_sort_node(l, t.t, p);
        } else {
            node_type_error(t.t, 0, GenSym.series, l, p);
        }
        break;
    case GenFunc.f_flatten:
    case GenFunc.f_instrings:
    case GenFunc.f_asort:
        if (l.t == GenSym.array) {
            if (t.t == GenFunc.f_instrings) {
                ret = array_func_node(l, m, r, t.t, p);;
            } else {
                ret = array_func_node(l, r, null, t.t, p);
            }
        } else {
            node_type_error(t.t, 0, GenSym.array, l, p);
        }
        break;
    case GenFunc.f_values:
    case GenFunc.f_uniq:
    case GenFunc.f_pergm:
    case GenFunc.f_irr:
        /* series or vector argument needed */
        if (l.t == GenSym.series || l.t == GenSym.mat || l.t == GenSym.num) {
            if (t.t == GenFunc.f_pergm) {
                ret = pergm_node(l, r, p);
            } else if (t.t == GenFunc.f_values || t.t == GenFunc.f_uniq) {
                ret = vector_values(l, t.t, p);
            } else if (t.t == GenFunc.f_irr) {
                ret = do_irr(l, p);
            } else {
                ret = vector_sort(l, t.t, p);
            }
        } else {
            p.err = e_types(l);
        }
        break;
    case GenFunc.f_sum:
    case GenFunc.f_sumall:
    case GenFunc.f_mean:
    case GenFunc.f_sd:
    case GenFunc.f_vce:
    case GenFunc.f_sst:
    case GenFunc.f_skewness:
    case GenFunc.f_kurtosis:
    case GenFunc.f_median:
    case GenFunc.f_gini:
    case GenFunc.f_nobs:
        /* functions taking series arg (mostly), returning scalar */
        if (l.t == GenSym.series || l.t == GenSym.mat) {
            ret = series_scalar_func(l, t.t, r, p);
        } else if ((t.t == GenFunc.f_mean || t.t == GenFunc.f_sd ||
                    t.t == GenFunc.f_vce || t.t == GenFunc.f_sum ||
                    t.t == GenFunc.f_median || t.t == GenFunc.f_nobs)
                   && ok_list_node(l, p)) {
            /* list -> series also acceptable for these cases */
            ret = list_to_series_func(l, t.t, p);
        } else if (l.t == GenSym.num) {
            ret = pretend_matrix_scalar_func(l, t.t, p);
        } else {
            node_type_error(t.t, 0, GenSym.series, l, p);
        }
        break;
    case GenFunc.f_t1:
    case GenFunc.f_t2:
        /* firstobs() and lastobs() */
        if (l.t == GenSym.series) {
            ret = series_scalar_func(l, t.t, r, p);
        } else {
            node_type_error(t.t, 0, GenSym.series, l, p);
        }
        break;
    case GenFunc.f_min:
    case GenFunc.f_max:
        if (null_node(r)) {
            /* "within arg" min or max */
            ret = within_minmax_node(l, t.t, p);
        } else {
            /* "between args" min or max */
            ret = between_minmax_node(l, r, t.t, p);
        }
        break;
    case GenFunc.f_ecdf:
    case GenFunc.f_normtest:
        /* series or vector (plus optional string arg for
           normtest); returns matrix */
        if (l.t != GenSym.series && l.t != GenSym.mat) {
            p.err = e_types(l);
        } else if (null_or_string(r)) {
            ret = series_matrix_node(l, r, t.t, p);
        } else {
            p.err = e_types(r);
        }
        break;
    case GenFunc.f_lrvar:
    case GenFunc.f_isconst:
    case GenFunc.f_isdummy:
        /* takes series and scalar arg, returns scalar */
        if (l.t == GenSym.series || l.t == GenSym.mat) {
            if (t.t == GenFunc.f_isconst || t.t == GenFunc.f_isdummy ) {
                ret = isconst_or_dum_node(l, r, p, t.t);
            } else if (t.t == GenFunc.f_lrvar) {
                ret = series_scalar_scalar_func(l, m, r, t.t, p);
            }
        } else {
            node_type_error(t.t, 1, GenSym.series, l, p);
        }
        break;
    case GenFunc.f_npv:
        if (l.t != GenSym.series && l.t != GenSym.mat && l.t != GenSym.num) {
            node_type_error(t.t, 1, GenSym.series, l, p);
        } else if (!scalar_node(r)) {
            node_type_error(t.t, 2, GenSym.num, r, p);
        } else {
            ret = series_scalar_scalar_func(l, r, null, t.t, p);
        }
        break;
    case GenFunc.f_quantile:
        if (l.t == GenSym.series) {
            if (scalar_node(r)) {
                ret = series_scalar_scalar_func(l, r, null, t.t, p);
            } else {
                node_type_error(t.t, 2, GenSym.num, r, p);
            }
        } else if (l.t == GenSym.mat) {
            if (r.t == GenSym.mat || scalar_node(r)) {
                ret = matrix_quantiles_node(l, r, p);
            } else {
                node_type_error(t.t, 2, GenSym.mat, r, p);
            }
        } else {
            node_type_error(t.t, 1, (r.t == GenSym.mat)? GenSym.mat : GenSym.series,
                            l, p);
        }
        break;
    case GenFunc.f_runiform:
    case GenFunc.f_rnormal:
        /* functions taking zero or two scalars as args */
        if (scalar_node(l) && scalar_node(r)) {
            ret = series_fill_func(l, r, t.t, p);
        } else if (l.t == GenSym.empty && r.t == GenSym.empty) {
            ret = series_fill_func(l, r, t.t, p);
        } else {
            node_type_error(t.t, (l.t == GenSym.num)? 2 : 1,
                            GenSym.num, (l.t == GenSym.num)? r : l, p);
        }
        break;
    case GenFunc.f_cov:
    case GenFunc.f_cor:
    case GenFunc.f_naalen:
    case GenFunc.f_kmeier:
    case GenFunc.f_corresp:
        /* functions taking two series/vectors as args, mostly */
        if ((l.t == GenSym.series || l.t == GenSym.mat || l.t == GenSym.num) &&
            (r.t == GenSym.series || r.t == GenSym.mat || r.t == GenSym.num)) {
            ret = series_2_func(l, r, t.t, p);
        } else if ((l.t == GenSym.series || l.t == GenSym.mat) &&
                   null_node(r) &&
                   (t.t == GenFunc.f_naalen || t.t == GenFunc.f_kmeier)) {
            ret = series_2_func(l, null, t.t, p);
        } else {
            p.err = GretlError.invarg;
        }
        break;
    case GenFunc.f_fcstats:
        /* two series or vectors, plus optional boolean */
        if ((l.t == GenSym.series || l.t == GenSym.mat) &&
            (m.t == GenSym.series || m.t == GenSym.mat || m.t == GenSym.list)) {
            ret = fcstats_node(l, m, r, p);
        } else {
            p.err = GretlError.invarg;
        }
        break;
    case GenFunc.f_npcorr:
        /* two series or vectors, plus optional control string */
        if ((l.t == GenSym.series || l.t == GenSym.mat) &&
            (m.t == GenSym.series || m.t == GenSym.mat) &&
            null_or_string(r)) {
            ret = npcorr_node(l, m, r, p);
        } else {
            p.err = GretlError.invarg;
        }
        break;
    case GenFunc.f_mxtab:
        /* functions taking two series or matrices as args and returning
           a matrix */
        if ((l.t == GenSym.series && r.t == GenSym.series) || (l.t == GenSym.mat && r.t == GenSym.mat)) {
            ret = mxtab_func(l, r, p);
        } else {
            node_type_error(t.t, (l.t == GenSym.series)? 2 : 1,
                            GenSym.series, (l.t == GenSym.series)? r : l, p);
        }
        break;
    case GenFunc.f_sortby:
        /* takes two series as args, returns series */
        if (l.t == GenSym.series && r.t == GenSym.series) {
            ret = series_sort_by(l, r, p);
        } else {
            node_type_error(t.t, (l.t == GenSym.series)? 2 : 1,
                            GenSym.series, (l.t == GenSym.series)? r : l, p);
        }
        break;
    case GenFunc.f_imat:
    case GenFunc.f_zeros:
    case GenFunc.f_ones:
    case GenFunc.f_munif:
    case GenFunc.f_mnorm:
    case GenFunc.f_randperm:
        /* matrix-creation functions */
        if (scalar_node(l) && null_or_scalar(r)) {
            ret = matrix_fill_func(l, r, t.t, p);
        } else if (!scalar_node(l)) {
            node_type_error(t.t, 1, GenSym.num, l, p);
        } else {
            node_type_error(t.t, 2, GenSym.num, r, p);
        }
        break;
    case GenFunc.hf_vcnorm:
        if (l.t == GenSym.mat && scalar_node(r)) {
            ret = mvnormal_node(l, r, p);
        } else if (l.t != GenSym.mat) {
            node_type_error(t.t, 1, GenSym.mat, l, p);
        } else {
            node_type_error(t.t, 2, GenSym.num, r, p);
        }
        break;
    case GenFunc.f_sumc:
    case GenFunc.f_sumr:
    case GenFunc.f_prodc:
    case GenFunc.f_prodr:
    case GenFunc.f_meanc:
    case GenFunc.f_meanr:
    case GenFunc.f_mcov:
    case GenFunc.f_mcorr:
    case GenFunc.f_chol:
    case GenFunc.f_psdroot:
    case GenFunc.f_inv:
    case GenFunc.f_ginv:
    case GenFunc.f_diag:
    case GenFunc.f_transp:
    case GenFunc.f_mrev:
    case GenFunc.f_vec:
    case GenFunc.f_vech:
    case GenFunc.f_unvech:
    case GenFunc.f_upper:
    case GenFunc.f_lower:
    case GenFunc.f_nullspc:
    case GenFunc.f_mexp:
    case GenFunc.f_mlog:
    case GenFunc.f_minc:
    case GenFunc.f_maxc:
    case GenFunc.f_minr:
    case GenFunc.f_maxr:
    case GenFunc.f_iminc:
    case GenFunc.f_imaxc:
    case GenFunc.f_iminr:
    case GenFunc.f_imaxr:
    case GenFunc.f_fft:
    case GenFunc.f_ffti:
    case GenFunc.f_polroots:
    case GenFunc.f_permute:
    case GenFunc.f_ctrans:
    case GenFunc.f_which:
        /* matrix -> matrix functions */
        if (l.t == GenSym.mat || l.t == GenSym.num) {
            ret = matrix_to_matrix_func(l, r, t.t, p);
        } else if (t.t == GenFunc.f_mrev && l.t == GenSym.list) {
            ret = list_reverse_node(l, p);
        } else {
            node_type_error(t.t, 1, GenSym.mat, l, p);
        }
        break;
    case GenFunc.hf_glasso:
        if (l.t == GenSym.mat && r.t == GenSym.bundle) {
            ret = glasso_node(l, r, p);
        } else {
            p.err = GretlError.types;
        }
        break;
    case GenFunc.f_invpd:
        /* matrix -> matrix with optional scalar pointer */
        if (l.t == GenSym.mat || l.t == GenSym.num) {
            ret = invpd_node(l, r, p);
        } else {
            node_type_error(t.t, 1, GenSym.mat, l, p);
        }
        break;
    case GenFunc.f_rows:
    case GenFunc.f_cols:
    case GenFunc.f_norm1:
    case GenFunc.f_infnorm:
    case GenFunc.f_rcond:
    case GenFunc.f_cnumber:
    case GenFunc.f_rank:
        /* matrix -> scalar functions */
        if (l.t == GenSym.mat || l.t == GenSym.num) {
            ret = matrix_to_scalar_func(l, r, t.t, p);
        } else {
            node_type_error(t.t, 0, GenSym.mat, l, p);
        }
        break;
    case GenFunc.f_det:
    case GenFunc.f_ldet:
    case GenFunc.f_trace:
        if (l.t == GenSym.mat || l.t == GenSym.num) {
            ret = matrix_to_alt_node(l, t.t, p);
        }
        break;
    case GenFunc.f_mread:
    case GenFunc.f_bread:
        if (l.t != GenSym.str) {
            node_type_error(t.t, 1, GenSym.str, l, p);
        } else if (!null_or_scalar(r)) {
            node_type_error(t.t, 2, GenSym.num, r, p);
        } else {
            ret = read_object_func(l, r, t.t, p);
        }
        break;
    case GenFunc.f_eigsym:
        /* matrix -> matrix functions, with single indirect return */
        if (l.t != GenSym.mat && l.t != GenSym.num) {
            node_type_error(t.t, 1, GenSym.mat, l, p);
        } else if (r.t != GenSym.u_addr && r.t != GenSym.empty) {
            node_type_error(t.t, 2, GenSym.u_addr, r, p);
        } else {
            ret = matrix_to_matrix2_func(l, r, t.t, p);
        }
        break;
    case GenFunc.f_complex:
        if ((l.t == GenSym.mat || l.t == GenSym.num) &&
            (r.t == GenSym.mat || null_or_scalar(r))) {
            ret = complex_matrix_node(l, r, p);
        } else {
            p.err = GretlError.types;
        }
        break;
    case GenFunc.f_fdjac:
    case GenFunc.f_numhess:
        /* matrix, fncall, optional scalar */
        if (l.t == GenSym.mat && m.t == GenSym.str) {
            ret = numeric_jacobian_or_hessian(l, m, r, t.t, p);
        } else {
            p.err = GretlError.types;
        }
        break;
    case GenFunc.f_mwrite:
        /* matrix, with string as second arg */
        if ((l.t == GenSym.mat || l.t == GenSym.num) &&
            m.t == GenSym.str && null_or_scalar(r)) {
            ret = matrix_file_write(l, m, r, p);
        } else {
            p.err = GretlError.types;
        }
        break;
    case GenFunc.f_bwrite:
        /* bundle, with string as second arg */
        if (l.t == GenSym.bundle && m.t == GenSym.str && null_or_scalar(r)) {
            ret = bundle_file_write(l, m, r, p);
        } else {
            p.err = GretlError.types;
        }
        break;
    case GenFunc.f_bfgsmax:
        /* matrix-pointer, plus one or two string args */
        if ((l.t == GenSym.u_addr || l.t == GenSym.mat) && m.t == GenSym.str) {
            ret = BFGS_maximize(l, m, r, p, t);
        } else {
            p.err = GretlError.types;
        }
        break;
    case GenFunc.f_bfgscmax:
        ret = BFGS_constrained_max(t, multi, p);
        break;
    case GenFunc.f_simann:
    case GenFunc.f_nmmax:
    case GenFunc.f_gssmax:
        /* matrix(-pointer), plus string and scalar args */
        if ((l.t == GenSym.u_addr || l.t == GenSym.mat) && m.t == GenSym.str) {
            ret = deriv_free_node(l, m, r, p, t);
        } else {
            p.err = GretlError.types;
        }
        break;
    case GenFunc.f_fzero:
        if (l.t == GenSym.str) {
            ret = fzero_node(l, m, r, p);
        } else {
            p.err = e_types(l);
        }
        break;
    case GenFunc.f_imhof:
        /* matrix, scalar as second arg */
        if (l.t == GenSym.mat && scalar_node(r)) {
            ret = matrix_imhof(l, r, p);
        } else {
            p.err = GretlError.types;
        }
        break;
    case GenFunc.f_bkw:
        /* matrix, string(s) as second optional arg,
           quiet flag as optional third arg */
        if (l.t == GenSym.mat) {
            ret = bkw_node(l, m, r, p);
        } else {
            p.err = e_types(l);
        }
        break;
    case GenFunc.f_fevd:
        /* integer target, source plus optional bundle */
        if (scalar_node(l) && null_or_scalar(m)) {
            ret = fevd_node(l, m, r, p);
        } else {
            p.err = GretlError.types;
        }
        break;
    case GenFunc.f_cnameset:
    case GenFunc.f_rnameset:
        /* matrix, with (list, string or strings array) as second arg */
        if (l.t == GenSym.mat && (ok_list_node(r, p) || r.t == GenSym.str || r.t == GenSym.array)) {
            ret = matrix_add_names(l, r, t.t, p);
        } else {
            p.err = GretlError.types;
        }
        break;
    case GenFunc.f_cnameget:
    case GenFunc.f_rnameget:
        /* matrix, scalar as second arg */
        if (l.t == GenSym.mat && null_or_scalar(r)) {
            ret = matrix_get_col_or_row_name(t.t, l, r, p);
        } else {
            p.err = GretlError.types;
        }
        break;
    case GenFunc.f_randint:
    case GenFunc.f_kpsscrit:
        /* two scalars */
        if (scalar_node(l) && scalar_node(r)) {
            if (t.t == GenFunc.f_kpsscrit) {
                ret = kpss_crit_node(l, r, p);
            } else {
                ret = two_scalars_func(l, r, t.t, p);
            }
        } else {
            p.err = GretlError.types;
        }
        break;
    case GenFunc.f_mshape:
    case GenFunc.f_svd:
    case GenFunc.f_qr:
    case GenFunc.f_eigen:
    case GenFunc.f_schur:
    case GenFunc.f_trimr:
    case GenFunc.f_corrgm:
    case GenFunc.f_seq:
    case GenFunc.f_replace:
    case GenFunc.f_strncmp:
    case GenFunc.f_weekday:
    case GenFunc.f_dayspan:
    case GenFunc.f_smplspan:
    case GenFunc.f_monthlen:
    case GenFunc.f_epochday:
    case GenFunc.f_kdensity:
    case GenFunc.f_setnote:
    case GenFunc.f_bwfilt:
    case GenFunc.f_varsimul:
    case GenFunc.f_mlag:
    case GenFunc.f_eigsolve:
    case GenFunc.f_princomp:
    case GenFunc.f_halton:
    case GenFunc.f_aggrby:
    case GenFunc.f_iwishart:
    case GenFunc.f_mweights:
    case GenFunc.f_mgradient:
    case GenFunc.f_lrcovar:
    case GenFunc.f_brename:
    case GenFunc.f_isoweek:
    case GenFunc.f_stack:
    case GenFunc.f_vma:
    case GenFunc.f_bcheck:
    case GenFunc.f_sphcorr:
    case GenFunc.f_sdc:
    case GenFunc.f_cdemean:
    case GenFunc.f_tsolvepd:
    case GenFunc.f_kmeans:
    case GenFunc.hf_regls:
        /* built-in functions taking three args */
        if (t.t == GenFunc.f_replace) {
            ret = replace_value(l, m, r, p);
        } else if (t.t == GenFunc.f_epochday) {
            ret = eval_epochday(l, m, r, p);
        } else {
            ret = eval_3args_func(l, m, r, t.t, p);
        }
        break;
    case GenFunc.f_geoplot:
        ret = geoplot_node(l, m, r, p);
        break;
    case GenFunc.f_printf:
    case GenFunc.f_sprintf:
        if (l.t == GenSym.str && null_or_string(r)) {
            ret = print_func_node(l, r, t.t, p);
        } else {
            p.err = GretlError.invarg;
        }
        break;
    case GenFunc.f_sscanf:
        if ((l.t == GenSym.str || l.t == GenSym.array) &&
            (m.t == GenSym.str && r.t == GenSym.str)) {
            ret = sscanf_node(l, m, r, p);
        } else {
            p.err = GretlError.invarg;
        }
        break;
    case GenFunc.f_bessel:
        /* functions taking one char, one scalar and one
           scalar/series/matrix as args */
        if (l.t != GenSym.str) {
            node_type_error(t.t, 1, GenSym.str, l, p);
        } else if (!scalar_node(m)) {
            node_type_error(t.t, 2, GenSym.num, m, p);
        } else if (r.t != GenSym.num && r.t != GenSym.series && r.t != GenSym.mat) {
            node_type_error(t.t, 3, GenSym.num, r, p);
        } else {
            ret = eval_bessel_func(l, m, r, p);
        }
        break;
    case GenFunc.f_bkfilt:
    case GenFunc.f_mols:
    case GenFunc.f_mpols:
    case GenFunc.f_mrls:
    case GenFunc.f_fols:
    case GenFunc.f_filter:
    case GenFunc.f_mcovg:
    case GenFunc.f_nrmax:
    case GenFunc.f_loess:
    case GenFunc.f_ghk:
    case GenFunc.f_quadtab:
    case GenFunc.f_qlrpval:
    case GenFunc.f_bootci:
    case GenFunc.f_bootpval:
    case GenFunc.f_movavg:
    case GenFunc.f_irf:
    case GenFunc.f_nadarwat:
    case GenFunc.f_feval:
    case GenFunc.f_chowlin:
    case GenFunc.f_hyp2f1:
    case GenFunc.f_tdisagg:
    case GenFunc.f_midasmult:
    case GenFunc.f_commute:
    case GenFunc.f_toepsolv:
    case GenFunc.f_rgbmix:
    case GenFunc.hf_felogitr:
        /* built-in functions taking more than three args */
        if (multi == null) {
            fprintf(stderr, "INTERNAL ERROR: @multi is NULL\n".ptr);
            p.err = GretlError.data;
        } else if (t.t == GenFunc.f_feval) {
            multi.parent = t;
            ret = eval_feval(t.t, multi, null, p);
        } else {
            ret = eval_nargs_func(t, multi, p);
        }
        break;
    case GenFunc.f_fevalb:
        if (l.t == GenSym.str && r.t == GenSym.bundle) {
            l.parent = t;
            ret = eval_feval(t.t, l, r, p);
        } else {
            p.err = GretlError.types;
        }
        break;
    case GenFunc.f_defarray:
    case GenFunc.f_defbundle:
    case GenFunc.f_defargs:
    case GenFunc.f_deflist:
        ret = object_def_node(t, multi, p);
        break;
    case GenFunc.f_defmat:
        ret = matrix_def_node(t, multi, p);
        break;
    case GenFunc.f_ksetup:
    case GenFunc.f_kfilter:
    case GenFunc.f_ksmooth:
    case GenFunc.f_ksimul:
    case GenFunc.f_kdsmooth:
        if (t.t == GenFunc.f_ksetup || bundle_pointer_arg0(t)) {
            ret = eval_kalman_bundle_func(t, multi, p);
        } else {
            p.err = GretlError.types;
        }
        break;
    case GenFunc.f_ksimdata:
        if (r.t == GenSym.mat) {
            ret = kalman_data_node(l, r, p);
        } else {
            node_type_error(t.t, 2, GenSym.mat, r, p);
        }
        break;
    case GenFunc.f_isoconv:
        ret = isoconv_node(t, multi, p);
        break;
    case GenSym.mvar:
        /* variable "under" model */
        ret = model_var_node(t, null, p);
        break;
    case GenSym.dmstr:
        ret = dollar_str_node(t, null, p);
        break;
    case GenSym.dvar:
        /* dataset "dollar" variable */
        ret = dollar_var_node(t, p);
        break;
    case GenFunc.f_obsnum:
    case GenFunc.f_isdiscr:
    case GenFunc.f_nlines:
    case GenFunc.f_remove:
    case GenFunc.f_exists:
        if (l.t == GenSym.str) {
            ret = object_status(l, t, p);
        } else {
            node_type_error(t.t, 1, GenSym.str, l, p);
        }
        break;
    case GenFunc.f_strlen:
        if (l.t == GenSym.str) {
            ret = object_status(l, t, p);
        } else if (useries_node(l) || l.t == GenSym.array) {
            ret = multi_str_node(l, t.t, p);
        } else {
            node_type_error(t.t, 1, GenSym.str, l, p);
        }
        break;
    case GenFunc.f_typeof:
    case GenFunc.f_typename:
    case GenFunc.f_iscmplx:
        ret = generic_typeof_node(l, t, p);
        break;
    case GenFunc.f_nelem:
        ret = n_elements_node(l, p);
        break;
    case GenFunc.f_inlist:
    case GenFunc.hf_listinfo:
        if (l.t == GenSym.list && t.t == GenFunc.hf_listinfo) {
            ret = list_info_node(l, r, p);
        } else if (ok_list_node(l, p)) {
            ret = in_list_node(l, r, p);
        } else {
            node_type_error(t.t, 1, GenSym.list, l, p);
        }
        break;
    case GenFunc.f_pdf:
    case GenFunc.f_cdf:
    case GenFunc.f_invcdf:
    case GenFunc.f_crit:
    case GenFunc.f_pval:
    case GenFunc.f_randgen:
    case GenFunc.f_mrandgen:
    case GenFunc.f_randgen1:
    case GenFunc.f_urcpval:
        if (t.t == GenFunc.f_urcpval) {
            ret = eval_urcpval(t, multi, p);
        } else {
            ret = eval_pdist(t, multi, p);
        }
        break;
    case GenSym.con:
        /* built-in constant */
        ret = retrieve_const(t, p);
        break;
    case GenSym.ufun:
        ret = eval_ufunc(t, multi, p);
        break;
    case GenSym.rfun:
        ret = eval_Rfunc(t, multi, p);
        break;
    case GenSym.query:
        ret = eval_query(t, l, p);
        break;
    case GenSym.b_lcat:
        /* list concatenation */
        if (ok_list_node(l, p) && ok_list_node(r, p)) {
            ret = eval_lcat(l, r, p);
        } else {
            p.err = GretlError.types;
        }
        break;
    case GenFunc.f_square:
        if (ok_list_node(l, p) && null_or_scalar(r)) {
            ret = apply_list_func(l, r, t.t, p);
        } else {
            p.err = GretlError.types;
        }
        break;
    case GenFunc.f_wmean:
    case GenFunc.f_wvar:
    case GenFunc.f_wsd:
        /* two lists -> series, with optional boolean */
        if (ok_list_node(l, p) && ok_list_node(m, p)) {
            ret = list_list_series_func(l, m, t.t, r, p);
        } else {
            p.err = GretlError.types;
        }
        break;
    case GenFunc.f_lincomb:
        /* list + matrix -> series */
        if (ok_list_node(l, p) && ok_matrix_node(r)) {
            ret = lincomb_func(l, r, null, t.t, p);
        } else {
            p.err = GretlError.types;
        }
        break;
    case GenFunc.f_mlincomb:
        /* list + matrix + int -> series */
        if (ok_list_node(l, p) && m.t == GenSym.mat &&
            (scalar_node(r) || r.t == GenSym.str)) {
            ret = lincomb_func(l, m, r, t.t, p);
        } else {
            p.err = GretlError.types;
        }
        break;
    case GenFunc.f_argname:
    case GenFunc.f_backtick:
    case GenFunc.f_strstrip:
    case GenFunc.f_fixname:
        if (l.t == GenSym.str) {
            ret = single_string_func(l, r, t.t, p);
        } else if (t.t == GenFunc.f_argname && uvar_node(l)) {
            ret = argname_from_uvar(l, r, p);
        } else {
            node_type_error(t.t, 0, GenSym.str, l, p);
        }
        break;
    case GenFunc.f_ccode:
        ret = country_code_node(l, r, p);
        break;
    case GenFunc.f_readfile:
        if (l.t == GenSym.str) {
            ret = readfile_node(l, r, p);
        } else {
            node_type_error(t.t, 1, GenSym.str, l, p);
        }
        break;
    case GenFunc.f_getenv:
    case GenFunc.f_ngetenv:
        if (l.t == GenSym.str) {
            ret = do_getenv(l, t.t, p);
        } else {
            node_type_error(t.t, 0, GenSym.str, l, p);
        }
        break;
    case GenFunc.f_errorif:
        if (r.t != GenSym.str) {
            p.err = e_types(r);
        } else {
            ret = do_errorif(l, r, p);
        }
        break;
    case GenFunc.f_assert:
        if (l.t == GenSym.num && r.t == GenSym.str) {
            ret = do_assert(l, r, p);
        } else {
            p.err = GretlError.types;
        }
        break;
    case GenFunc.f_contains:
        if ((r.t == GenSym.mat || r.t == GenSym.num) &&
            (l.t == GenSym.num || l.t == GenSym.series || l.t == GenSym.mat)) {
            ret = contains_node(l, r, p);
        } else {
            p.err = GretlError.types;
        }
        break;
    case GenFunc.f_obslabel:
        if (l.t == GenSym.num || l.t == GenSym.mat) {
            ret = int_to_string_func(l, t.t, p);
        } else {
            node_type_error(t.t, 0, GenSym.num, l, p);
        }
        break;
    case GenFunc.f_varname:
        if (l.t == GenSym.num || l.t == GenSym.mat || l.t == GenSym.series) {
            ret = int_to_string_func(l, t.t, p);
        } else if (l.t == GenSym.list) {
            ret = list_to_string_func(l, t.t, p);
        } else {
            node_type_error(t.t, 0, GenSym.num, l, p);
        }
        break;
    case GenFunc.f_randstr:
        if (scalar_node(l)) {
            ret = int_to_string_func(l, t.t, p);
        } else {
            node_type_error(t.t, 0, GenSym.num, l, p);
        }
        break;
    case GenFunc.f_varnames:
        if (l.t == GenSym.list) {
            ret = list_to_string_func(l, t.t, p);
        } else {
            node_type_error(t.t, 0, GenSym.list, l, p);
        }
        break;
    case GenFunc.f_varnum:
        if (l.t == GenSym.str) {
            ret = varnum_node(l, p);
        } else {
            node_type_error(t.t, 0, GenSym.str, l, p);
        }
        break;
    case GenFunc.f_tolower:
    case GenFunc.f_toupper:
        if (l.t == GenSym.str || l.t == GenSym.array || useries_node(l)) {
            ret = variant_string_func(l, null, null, t.t, p);
        } else {
            p.err = GretlError.types;
        }
        break;
    case GenFunc.f_strsub:
    case GenFunc.f_regsub:
        if ((l.t == GenSym.str || l.t == GenSym.array || strvals_node(l, p)) &&
            m.t == GenSym.str && r.t == GenSym.str) {
            ret = variant_string_func(l, m, r, t.t, p);
        } else {
            p.err = GretlError.types;
        }
        break;
    case GenFunc.f_substr:
        if ((l.t == GenSym.str || l.t == GenSym.array || strvals_node(l, p)) &&
            scalar_node(m) && scalar_node(r)) {
            ret = variant_string_func(l, m, r, t.t, p);
        } else {
            p.err = GretlError.types;
        }
        break;
    case GenFunc.f_jsonget:
    case GenFunc.f_xmlget:
        if (l.t == GenSym.str && m.t == GenSym.str) {
            ret = two_string_func(l, m, r, t.t, p);
        } else if (l.t == GenSym.str && m.t == GenSym.array && t.t == GenFunc.f_xmlget) {
            ret = two_string_func(l, m, r, t.t, p);
        } else {
            node_type_error(t.t, (l.t == GenSym.str)? 2 : 1,
                            GenSym.str, (l.t == GenSym.str)? m : l, p);
        }
        break;
    case GenFunc.f_jsongetb:
        if (l.t == GenSym.str && null_or_string(m)) {
            ret = two_string_func(l, m, r, t.t, p);
        } else {
            p.err = GretlError.types;
        }
        break;
    case GenFunc.f_strstr:
    case GenFunc.f_instring:
        if (l.t == GenSym.str && m.t == GenSym.str) {
            ret = two_string_func(l, m, r, t.t, p);
        } else {
            node_type_error(t.t, (l.t == GenSym.str)? 2 : 1,
                            GenSym.str, (l.t == GenSym.str)? m : l, p);
        }
        break;
    case GenFunc.f_strsplit:
        if (l.t == GenSym.str) {
            ret = strsplit_node(t.t, l, m, r, p);
        } else {
            node_type_error(t.t, 1, GenSym.str, l, p);
        }
        break;
    case GenFunc.f_getline:
        if (l.t == GenSym.str && (null_or_string(r) || r.t == GenSym.u_addr)) {
            ret = getline_node(l, r, p);
        } else {
            node_type_error(t.t, (l.t == GenSym.str)? 2 : 1,
                            GenSym.str, (l.t == GenSym.str)? r : l, p);
        }
        break;
    case GenFunc.f_errmsg:
        if (null_or_scalar(l)) {
            ret = errmsg_node(l, p);
        } else {
            node_type_error(t.t, 0, GenSym.num, l, p);
        }
        break;
    case GenFunc.f_isodate:
    case GenFunc.f_juldate:
        ret = isodate_node(l, r, t, p);
        break;
    case GenFunc.f_strftime:
    case GenFunc.f_strfday:
        if (l.t == GenSym.num || l.t == GenSym.series || l.t == GenSym.mat) {
            if (t.t == GenFunc.f_strftime) {
                ret = strftime_node(l, m, r, t.t, 0, p);
            } else {
                ret = strftime_node(l, r, null, t.t, 0, p);
            }
        } else {
            node_type_error(t.t, 0, GenSym.num, l, p);
        }
        break;
    case GenFunc.f_strptime:
    case GenFunc.f_strpday:
        if (l.t == GenSym.str || l.t == GenSym.series || l.t == GenSym.array ||
            l.t == GenSym.mat || l.t == GenSym.num) {
            ret = strptime_node(l, r, t.t, p);
        } else {
            node_type_error(t.t, 0, GenSym.str, l, p);
        }
        break;
    case GenFunc.f_atof:
        if (l.t == GenSym.str || l.t == GenSym.series) {
            ret = atof_node(l, p);
        } else {
            node_type_error(t.t, 0, GenSym.str, l, p);
        }
        break;
    case GenFunc.f_mpi_recv:
        ret = mpi_transfer_node(l, null, null, t.t, p);
        break;
    case GenFunc.f_mpi_send:
    case GenFunc.f_bcast:
    case GenFunc.f_allreduce:
        if (t.t == GenFunc.f_allreduce && r.t != GenSym.str) {
            node_type_error(t.t, 2, GenSym.str, r, p);
        } else {
            ret = mpi_transfer_node(l, r, null, t.t, p);
        }
        break;
    case GenFunc.f_reduce:
    case GenFunc.f_scatter:
        if (m.t != GenSym.str) {
            node_type_error(t.t, 2, GenSym.str, m, p);
        } else if (!null_or_scalar(r)) {
            node_type_error(t.t, 3, GenSym.num, r, p);
        } else {
            ret = mpi_transfer_node(l, m, r, t.t, p);
        }
        break;
    case GenFunc.f_genseries:
        ret = gen_series_node(l, r, p);
        break;
    case GenFunc.f_array:
        ret = gen_array_node(t, l, p);
        break;
    case GenFunc.f_strvals:
        if (!useries_node(l)) {
            node_type_error(t.t, 0, GenSym.useries, l, p);
        } else {
            ret = get_series_stringvals(l, r, p);
        }
        break;
    case GenFunc.f_stringify:
    case GenFunc.f_strvsort:
        if (useries_node(l) && (r.t == GenSym.array || null_node(r))) {
            ret = stringify_series(l, r, t.t, p);
        } else {
            p.err = GretlError.types;
        }
        break;
    case GenFunc.f_dec2bin:
        if (l.t == GenSym.num || l.t == GenSym.mat) {
            ret = dec2bin_node(l, p);
        } else {
            node_type_error(t.t, 0, GenSym.num, l, p);
        }
        break;
    case GenFunc.f_bin2dec:
        if (l.t == GenSym.mat) {
            ret = bin2dec_node(l, p);
        } else {
            node_type_error(t.t, 0, GenSym.mat, l, p);
        }
        break;
    case GenFunc.f_binperms:
        if (scalar_node(l) && scalar_node(r)) {
           ret = binperm_node(l, r, p);
        } else {
            p.err = GretlError.types;
        }
        break;
    case GenFunc.f_access:
        if (l.t == GenSym.str) {
            ret = access_series_node(l, p);
        } else {
            node_type_error(t.t, 0, GenSym.str, l, p);
        }
        break;
    default:
        fprintf(stderr, "eval: weird node %s (t->t = %d)\n".ptr,
                getsymb(t.t), t.t);
        p.err = GretlError.types;
        break;
    }

 finish:

    if ((p.flags & GenFlags.prnlist) && ret != null &&
        scalar_node(ret) && is_aux_node(ret)) {
        /* don't try to interpret this scalar as a series ID */
        p.err = GretlError.data;
    }

    if (!p.err && ret != null && ret != t && is_aux_node(ret)) {
        if ((t.t == GenFunc.f_feval || t.t == GenFunc.f_fevalb) && ret.refcount > 0) {
            ; /* special case: don't attach here */
        } else {
            p.err = attach_aux_node(t, ret, p);
        }
    }

    if (!p.err && t.t == GenFunc.f_ldiff) {
        if (!null_node(m) || !null_node(r)) {
            post_process_ldiff(ret, m, r, p);
        }
    }

 bailout:

    if (multi != null) {
        destroy_multi(multi);
    }


    return ret;
}


/* non-debugging variant: easier to see what's going on */

int attach_aux_node(Node * t, Node * ret, Parser * p)
{
    if (t.aux == null) {
        t.aux = ret;
        ret.refcount += 1;
    } else if (t.aux != ret) {
        if (t.t == GenSym.query || t.t == GenSym.b_and || t.t == GenSym.b_or) {
            /* OK to switch aux node in these cases */
            free_node(t.aux, p);
            t.aux = ret;
            ret.refcount += 1;
        } else if (is_proxy_node(t.aux)) {
            /* an extension to the above */
            free_node(t.aux, p);
            t.aux = ret;
            ret.refcount += 1;
        } else {
            /* otherwise if we're trying to switch aux node,
               something must have gone wrong
            */
            fprintf(stderr, "! %s node already has aux node attached\n".ptr,
                    getsymb(t.t));
            fprintf(stderr, " orig aux %p (%s), incoming %p (%s)\n".ptr,
                    cast(void *) t.aux, getsymb(t.aux.t),
                    cast(void *) ret, getsymb(ret.t));
            return GretlError.data;
        }
    }

    return 0;
}


int more_input(const(char) * s)
{
    while (*s) {
        if (!isspace(cast(ubyte) *s)) {
            return 1;
        }
        s++;
    }

    return 0;
}

/* get the next input character for the lexer */

int parser_getc(Parser * p)
{

    p.ch = 0;

    if (more_input(p.point)) {
        p.ch = *p.point;
        p.point += 1;
    }


    return p.ch;
}

/* advance the read position by n characters */

void parser_advance(Parser * p, int n)
{
    p.point += n;
    p.ch = *p.point;
    p.point += 1;
}

/* throw back the last-read character */

void parser_ungetc(Parser * p)
{
    p.point -= 1;
    p.ch = *(p.point - 1);
}

/* Look ahead for the first occurrence of a given character in
   the remaining input stream; return its 0-based index or
   -1 if not found.
*/

int parser_char_index(Parser * p, int c)
{
    int i;

    for (i=0; p.point[i] != '\0'; i++) {
        if (p.point[i] == c) {
            return i;
        }
    }

    return -1;
}

/* For error reporting: print the input up to the current
   parse point, unless it's not valid UTF-8. Return 0
   if the input is printed OK, otherwise non-zero.
*/

int parser_print_input(Parser * p)
{
    int len = cast(int) (p.point - p.input);
    char * s = gretl_strndup(p.input, len);
    int err = 0;

    if (s != null) {
        if (g_utf8_validate(s, -1, null)) {
            pprintf(p.prn, "> %s\n".ptr, s);
        } else {
            err = 1;
        }
        free(s);
    } else {
        err = 1;
    }

    return err;
}

/* "pretty print" syntactic nodes and symbols */

void printsymb(int symb, Parser * p)
{
    pputs(p.prn, getsymb(symb));
}

void printnode(Node * t, Parser * p, int value)
{
    if (t == null) {
        pputs(p.prn, "NULL");
    } else if (!value && useries_node(t)) {
        pprintf(p.prn, "%s".ptr, p.dset.varname[t.vnum]);
    } else if (!value && uscalar_node(t)) {
        pprintf(p.prn, "%s".ptr, t.vname);
    } else if (t.t == GenSym.num) {
        if (na(t.v.xval)) {
            pputs(p.prn, "NA");
        } else {
            pprintf(p.prn, "%.8g".ptr, t.v.xval);
        }
    } else if (t.t == GenSym.series) {
        const(double) * x = t.v.xvec;
        int i, j = 1;

        if (p.lh.vnum > 0 && p.lh.vnum < p.dset.v) {
            pprintf(p.prn, "%s\n".ptr, p.dset.varname[p.lh.vnum]);
        }

        for (i=p.dset.t1; i<=p.dset.t2; i++, j++) {
            if (na(x[i])) {
                pputs(p.prn, "NA");
            } else {
                pprintf(p.prn, "%g".ptr, x[i]);
            }
            if (j % 8 == 0) {
                pputc(p.prn, '\n');
            } else if (i < p.dset.t2) {
                pputc(p.prn, ' ');
            }
        }
    } else if (t.t == GenSym.mat) {
        if (t.vname != null) {
            pputs(p.prn, t.vname);
        } else {
            gretl_matrix_print_to_prn(t.v.m, null, p.prn);
        }
    } else if (t.t == GenSym.bundle) {
        gretl_bundle_print(t.v.b, p.prn);
    } else if (t.t == GenSym.dbundle) {
        pputs(p.prn, bvarname(t.v.idnum));
    } else if (t.t == GenSym.array) {
        gretl_array_print(t.v.a, p.prn);
    } else if (t.t == GenSym.uobj) {
        pprintf(p.prn, "%s".ptr, t.v.str);
    } else if (t.t == GenSym.dvar) {
        pputs(p.prn, dvarname(t.v.idnum));
    } else if (t.t == GenSym.mvar) {
        pputs(p.prn, mvarname(t.v.idnum));
    } else if (t.t == GenSym.con) {
        pputs(p.prn, constname(t.v.idnum));
    } else if (t.t == GenSym.dum) {
        pputs(p.prn, dumname(t.v.idnum));
    } else if (binary_op(t.t)) {
        pputc(p.prn, '(');
        printnode(t.L, p, 0);
        printsymb(t.t, p);
        printnode(t.R, p, 0);
        pputc(p.prn, ')');
    } else if (t.t == GenSym.osl) {
        printnode(t.L, p, 0);
        pputc(p.prn, '[');
        printnode(t.R, p, 0);
        pputc(p.prn, ']');
    } else if (t.t == GenSym.slraw) {
        pputs(p.prn, "SLRAW");
    } else if (t.t == GenSym.subsl) {
        pputs(p.prn, "SUBSL");
    } else if (func1_symb(t.t)) {
        printsymb(t.t, p);
        pputc(p.prn, '(');
        printnode(t.L, p, 0);
        pputc(p.prn, ')');
    } else if (unary_op(t.t)) {
        printsymb(t.t, p);
        printnode(t.L, p, 0);
    } else if (func2_symb(t.t)) {
        printsymb(t.t, p);
        pputc(p.prn, '(');
        printnode(t.L, p, 0);
        if (t.R.t != GenSym.empty) {
            pputc(p.prn, ',');
        }
        printnode(t.R, p, 0);
        pputc(p.prn, ')');
    } else if (t.t == GenSym.str) {
        pprintf(p.prn, "%s".ptr, t.v.str);
    } else if (t.t == GenSym.ufun) {
        pprintf(p.prn, "%s".ptr, t.vname);
    } else if (t.t == GenSym.dmstr) {
        printnode(t.L, p, 0);
        pputc(p.prn, '(');
        printnode(t.R, p, 0);
        pputc(p.prn, ')');
    } else if (t.t == GenSym.listvar) {
        pprintf(p.prn, "%s.%s".ptr, t.L.v.str, t.R.v.str);
    } else if (t.t == GenSym.list) {
        pputs(p.prn, "LIST");
    } else if (t.t == GenSym.lag) {
        pputs(p.prn, "LAG");
    } else if (t.t != GenSym.empty) {
        pputs(p.prn, "weird tree - ");
        printsymb(t.t, p);
    }
}

/* which modified assignment operators of the type '+='
   will we accept, when generating various types of
   result? */
bool ok_matrix_op(int o)
{
    return (o == GenSym.b_asn  || o == GenSym.b_dotasn || 
           o == GenSym.b_add  || o == GenSym.b_sub || 
           o == GenSym.b_mul  || o == GenSym.b_div || 
           o == GenSym.b_hcat || o == GenSym.b_vcat);
}

bool ok_list_op(int o)
{
    return (o == GenSym.b_asn || o == GenSym.b_add || o == GenSym.b_sub);
}

bool ok_string_op(int o)
{
    return (o == GenSym.b_asn || o == GenSym.b_add || 
           o == GenSym.b_hcat || o == GenSym.inc);
}

bool ok_array_op(int o)
{
    return (o == GenSym.b_asn || o == GenSym.b_add || o == GenSym.b_sub);
}

bool ok_bundle_op(int o)
{
    return (o == GenSym.b_asn || o == GenSym.b_add);
}

struct ModAssign {
    int c;
    int op;
}

/* supported "inflections" of assignment */

immutable ModAssign[] m_assign = [

    ModAssign('+', GenSym.b_add),
    ModAssign('-', GenSym.b_sub),
    ModAssign('*', GenSym.b_mul),
    ModAssign('/', GenSym.b_div),
    ModAssign('%', GenSym.b_mod),
    ModAssign('^', GenSym.b_pow),
    ModAssign('~', GenSym.b_hcat),
    ModAssign('|', GenSym.b_vcat),
    ModAssign('.', GenSym.b_dotasn),
    ModAssign(0, 0)
];

/* read operator from formula: this is either
   simple assignment or something like '+=' */

int get_op(char * s)
{
    if (s[0] == '=') {
        s[1] = '\0';
        return GenSym.b_asn;
    } else if (!strcmp(s, "++")) {
        return GenSym.inc;
    } else if (!strcmp(s, "--")) {
        return GenSym.dec;
    } else if (s[1] == '=') {
        int i;

        for (i=0; m_assign[i].c; i++) {
            if (s[0] == m_assign[i].c) {
                return m_assign[i].op;
            }
        }
    }

    return 0;
}

const(char) * get_opstr(int op)
{
    static __gshared char[4] opstr = 0;

    if (op == GenSym.b_asn) {
        return "=";
    } else if (op == GenSym.inc) {
        return "++";
    } else if (op == GenSym.dec) {
        return "--";
    } else {
        int i;

        for (i=0; m_assign[i].c; i++) {
            if (op == m_assign[i].op) {
                opstr[0] = cast(char) m_assign[i].c;
                opstr[1] = '=';
                return opstr.ptr;
            }
        }
        return "??";
    }
}

/* implement the declaration of new variables */

void do_declaration(Parser * p)
{
    char ** S = null;
    int i, v, n;

    n = check_declarations(&S, p);

    if (n == 0) {
        return;
    }

    for (i=0; i<n && !p.err; i++) {
        if (S[i] != null) {
            if (p.targ == GenSym.series) {
                p.err = dataset_add_NA_series(p.dset, 1);
                if (!p.err) {
                    v = p.dset.v - 1;
                    strcpy(p.dset.varname[v], S[i]);
                }
            } else {
                GretlType type = GretlType.none;

                if (p.targ == GenSym.mat) {
                    type = GretlType.matrix;
                } else if (p.targ == GenSym.num) {
                    type = GretlType.double_;
                } else if (p.targ == GenSym.str) {
                    type = GretlType.string_;
                } else if (p.targ == GenSym.bundle) {
                    type = GretlType.bundle;
                } else if (p.targ == GenSym.list) {
                    type = GretlType.list;
                } else if (p.targ == GenSym.array) {
                    type = p.lh.gtype;
                } else {
                    p.err = GretlError.data;
                }
                if (!p.err) {
                    p.err = user_var_add(S[i], type, null);
                }
            }
        }
    }

    strings_array_free(S, n);
}

/* The expression supplied for evaluation does not contain an '=':
   can we interpret it as an implicit request to print the value
   of an existing variable?
*/

void parser_try_print(Parser * p, const(char) * s, int * done)
{
    if (p.lh.t != 0 && p.lh.expr == null) {
        p.flags |= GenFlags.discard;
        p.point = s;
    } else {
        p.err = GretlError.eqn;
    }
}

/* Here we try to parse out the LHS of the statement and also the
   operator. If we find a unitary LHS (simply an identifier) we write
   it into p->lh.name, but if we find a compound LHS (such as a
   sub-matrix specification) we save it as p->lh.expr. The content of
   @ps is advanced to the first position beyond the operator.
*/

int extract_lhs_and_op(const(char) ** ps, Parser * p, char * opstr)
{
    const(char) * s = *ps;
    int quoted = 0;
    int i, n = 0;
    int err = 0;


    if (p.targ != GenSym.unk && strchr(s, '=') == null) {
        /* we got a type specification but no assignment,
           so should be variable declaration(s) ?
        */
        p.flags |= GenFlags.decl;
        p.lh.expr = gretl_strdup(s);
        goto done;
    }

    /* Count bytes preceding first unquoted '='. Note that
       the "unquoted" condition is required only because
       a string-literal bundle key might contain an equals
       sign, as in b["foo=bar"] = ...
    */
    for (i=0; s[i] != '\0'; i++) {
        if (s[i] == '"') {
            quoted = !quoted;
        } else if (!quoted && s[i] == '=') {
            break;
        }
        n++;
    }

    if (n > 0) {
        char * lhs = null;
        int lhlen = n;

        if (s[n] == '=') {
            /* we actually reached an '=' */
            if (s[n+1] == '=') {
                /* 2022-04-05: should this be an error? */
                goto done;
            } else if (strspn(s + n - 1, "+-*/%^~|.") == 1) {
                /* preceded by a modifier: inflected assignment */
                lhlen--;
                opstr[0] = s[n-1];
                opstr[1] = '=';
            } else {
                /* straight assignment */
                opstr[0] = '=';
            }
            n++; /* plus 1 for '=' */
        }

        if (lhlen > 0) {
            lhs = gretl_strndup(s, lhlen);
            tailstrip(lhs);
            lhlen = cast(int) strlen(lhs);
        }

        if (opstr[0] == '\0' && lhlen > 2) {
            /* check for postfix operator */
            char * test = lhs + lhlen - 2;

            if (!strcmp(test, "++") || !strcmp(test, "--")) {
                strcpy(opstr, test);
                *test = '\0';
                lhlen -= 2;
            }
        }

        if (lhlen > 0) {
            if (lhlen == gretl_namechar_spn(lhs)) {
                /* a straight identifier? */
                if (lhlen >= VNAMELEN) {
                    pprintf(p.prn, _("'%s': name is too long (max %d characters)\n").ptr,
                            lhs, VNAMELEN - 1);
                    err = GretlError.parse;
                } else {
                    strcpy(p.lh.name.ptr, lhs);
                }
            } else if ((p.flags & GenFlags.priv) && (*lhs == '$' || *lhs == '_') &&
                       gretl_namechar_spn(lhs + 1) == lhlen - 1) {
                /* "private" genr of the form $foo=expr or _foo=expr */
                strcpy(p.lh.name.ptr, lhs);
            } else {
                /* treat as an expression to be evaluated */
                p.lh.expr = lhs;
                lhs = null; /* protect against freeing */
            }
        } else {
            /* nothing relevant found */
            err = GretlError.parse;
        }

        if (!err && opstr[0] != '\0') {
            p.op = get_op(opstr);
        }

        free(lhs);
        *ps = s + n;
    }

 done:


    if (!p.err && !(p.flags & GenFlags.decl) && p.lh.name[0] == '\0' && p.op == 0) {
        /* added 2021-05-29 */
        p.flags |= GenFlags.discard;
    }

    return err;
}

void maybe_do_type_errmsg(const(char) * name, int t)
{
    const(char) * tstr = typestr(t);

    if (tstr != null && strcmp(tstr, "?")) {
        if (name != null && *name != '\0') {
            gretl_errmsg_sprintf(_("The variable %s is of type %s, not acceptable in context").ptr,
                                 name, tstr);
        } else {
            gretl_errmsg_sprintf(_("A variable of type %s is not acceptable in context").ptr, tstr);
        }
    }
}

void assignment_type_errmsg(int targ, int rhs, int op)
{
    const(char) * rhstr = typestr(rhs);

    if (*rhstr == '?') {
        rhstr = getsymb(rhs);
    } else if (targ == GenSym.list && rhs == GenSym.series) {
        rhstr = "anonymous series";
    }
    gretl_errmsg_sprintf(_("Incompatible types in assignment: %s %s %s").ptr, typestr(targ), get_opstr(op),
                         rhstr);
}

int overwrite_type_check(Parser * p)
{
    int err = 0;

    /* FIXME check for series/function collision here */

    if (p.targ != p.lh.t) {
        /* don't overwrite one type with another */
        maybe_do_type_errmsg(p.lh.name.ptr, p.lh.t);
        err = GretlError.types;
    }

    return err;
}

int overwrite_const_check(const(char) * name, int vnum)
{
    if (object_is_const(name, vnum)) {
        return overwrite_err(name);
    } else {
        return 0;
    }
}

/* Check that we're not trying to modify a const object via a compound
   LHS expression; and while we're at it, check whether we should be
   generating a list (a list member of a bundle or an element of an
   array of lists).
*/

int compound_const_check(Node * lhs, Parser * p)
{
    Node * n = lhs;
    int i = 0, err = 0;

    if (n.t == GenSym.bmemb && n.R != null && n.R.t == GenSym.str) {
        GretlType t;

        t = gretl_bundle_get_member_type(n.L.v.b, n.R.v.str, null);
        if (t == GretlType.list) {
            p.flags |= GenFlags.listdef;
        }
    }

    while (n.t == GenSym.osl || n.t == GenSym.obs || n.t == GenSym.bmemb) {
        n = n.L;
        if (i == 0 && lhs.t == GenSym.osl && n.t == GenSym.array) {
            if (gretl_array_get_type(n.v.a) == GretlType.lists) {
                p.flags |= GenFlags.listdef;
            }
        }
        i++;
    }

    /* do we have a const object at the tip of the tree? */
    if (n.vname != null) {
        err = overwrite_const_check(n.vname, n.vnum);
    }

    return err;
}

int ok_array_decl(Parser * p, const(char) * s)
{
    p.lh.gtype = GretlType.none;

    if (!strncmp(s, "strings ", 8)) {
        p.lh.gtype = GretlType.strings;
    } else if (!strncmp(s, "matrices ", 9)) {
        p.lh.gtype = GretlType.matrices;
    } else if (!strncmp(s, "bundles ", 8)) {
        p.lh.gtype = GretlType.bundles;
    } else if (!strncmp(s, "lists ", 6)) {
        p.lh.gtype = GretlType.lists;
    } else if (!strncmp(s, "arrays ", 7)) {
    p.lh.gtype = GretlType.arrays;
    }

    return p.lh.gtype != 0;
}

/* Given an existing LHS variable, whose type is recorded in p->lh.t,
   check that the specified operator is supported for the type. Return
   error code if not.
*/

int check_operator_validity(Parser * p, const(char) * opstr)
{
    if (p.lh.t == GenSym.mat && !ok_matrix_op(p.op)) {
        /* matrices: we accept only a limited range of
           modified assignment operators */
        gretl_errmsg_sprintf(_("'%s' : not implemented for matrices").ptr, opstr);
        return GretlError.parse;
    } else if (p.lh.t == GenSym.list && !ok_list_op(p.op)) {
        /* lists: same story as matrices */
        gretl_errmsg_sprintf(_("'%s' : not implemented for lists").ptr, opstr);
        return GretlError.parse;
    } else if (p.lh.t == GenSym.str && !ok_string_op(p.op)) {
        /* strings: ditto */
        gretl_errmsg_sprintf(_("'%s' : not implemented for strings").ptr, opstr);
        return GretlError.parse;
    } else if (p.lh.t == GenSym.array && !ok_array_op(p.op)) {
        /* arrays: ditto */
        gretl_errmsg_sprintf(_("'%s' : not implemented for arrays").ptr, opstr);
        return GretlError.parse;
    } else if (p.lh.t == GenSym.bundle && !ok_bundle_op(p.op)) {
        /* bundles: ditto */
        gretl_errmsg_sprintf(_("'%s' : not implemented for this type").ptr, opstr);
        return GretlError.parse;
    } else if (p.lh.t != GenSym.mat && (p.op == GenSym.b_vcat || p.op == GenSym.b_dotasn)) {
        /* vertical concat: only OK for matrices */
        gretl_errmsg_sprintf(_("'%s' : only defined for matrices").ptr, opstr);
        return GretlError.parse;
    } else if (p.lh.t != GenSym.mat && p.lh.t != GenSym.str && p.op == GenSym.b_hcat) {
        /* horizontal concat: only OK for matrices, strings */
        gretl_errmsg_sprintf(_("'%s' : not implemented for this type").ptr, opstr);
        return GretlError.parse;
    }

    /* otherwise OK? */
    return 0;
}

/* Do we have an inline type specification preceding the statement
   proper? In most cases we shouldn't, since it will already have been
   handled by the tokenizer (and the type will now be recorded in
   p->targ). But we allow for finding a typespec here in case of "genrs"
   within a block command, where the statement bypasses the regular
   tokenizer. (FIXME?)
*/

void check_for_inline_typespec(const(char) ** ps, Parser * p)
{
    const(char) * s = *ps;

    if (!strncmp(s, "scalar ", 7)) {
        p.targ = GenSym.num;
        s += 7;
    } else if (!strncmp(s, "series ", 7)) {
        p.targ = GenSym.series;
        s += 7;
    } else if (!strncmp(s, "matrix ", 7)) {
        p.targ = GenSym.mat;
        s += 7;
    } else if (!strncmp(s, "list ", 5)) {
        p.targ = GenSym.list;
        s += 5;
    } else if (!strncmp(s, "string ", 7)) {
        p.targ = GenSym.str;
        s += 7;
    } else if (!strncmp(s, "bundle ", 7)) {
        p.targ = GenSym.bundle;
        s += 7;
    } else if (ok_array_decl(p, s)) {
        p.targ = GenSym.array;
        s += strcspn(s, " ") + 1;
    }

    /* advance pointer */
    *ps = s;
}

/* Check @p->lh.name for the name of an existing series or user_var of
   some kind. If found, record the relevant info in p->lh.t, and also
   either p->lh.vnum (series) or p->lh.uv (other types).
*/

int check_existing_lhs_type(Parser * p, int * newvar)
{
    UserVar * uvar;
    int v, err = 0;

    if ((err = gretl_reserved_word(p.lh.name.ptr)) != 0) {
        return err;
    }

    v = current_series_index(p.dset, p.lh.name.ptr);
    if (v >= 0) {
        p.lh.vnum = v;
        p.lh.t = GenSym.series;
        *newvar = 0;
        return 0;
    }

    uvar = get_user_var_by_name(p.lh.name.ptr);

    if (uvar != null) {
        GretlType vtype = uvar.type;

        p.lh.uv = uvar;
        *newvar = 0;

        if (vtype == GretlType.matrix) {
            p.lh.t = GenSym.mat;
        } else if (vtype == GretlType.double_) {
            if (uvar.flags & UVFlags.nodecl) {
                if (p.targ == GenSym.unk) {
                    p.flags |= GenFlags.nodecl;
                } else {
                    uvar.flags &= ~UVFlags.nodecl;
                }
            }
            p.lh.t = GenSym.num;
        } else if (vtype == GretlType.list) {
            p.lh.t = GenSym.list;
        } else if (vtype == GretlType.string_) {
            p.lh.t = GenSym.str;
        } else if (vtype == GretlType.bundle) {
            p.lh.t = GenSym.bundle;
        } else if (vtype == GretlType.array) {
            p.lh.gtype = gretl_array_get_type(cast(GretlArray *) uvar.ptr);
            p.lh.t = GenSym.array;
        }
    }

    return err;
}

/* pre-process a "genr" statement */

void gen_preprocess(Parser * p, int flags, int * done)
{
    const(char) * s = p.input;
    char[3] opstr = 0;
    int newvar = 1;

    while (isspace(*s)) s++;

    /* skip leading command word, if any */
    if (!strncmp(s, "genr ", 5)) {
        s += 5;
    } else if (!strncmp(s, "print ", 6)) {
        /* allow this within (e.g.) mle block */
        p.flags |= GenFlags.discard;
        s += 6;
    }

    while (isspace(*s)) s++;

    if (p.targ == GenSym.unk) {
        check_for_inline_typespec(&s, p);
    } else if (gretl_array_type(cast(GretlType) p.targ)) {
        /* record a plural type spec such as "matrices" under
           the "lh" member of @p.
        */
        p.lh.gtype = cast(GretlType) p.targ;
        p.targ = GenSym.array;
    }

    /* check for types that cannot be generated in the
       absence of a dataset */
    if ((p.targ == GenSym.series || p.targ == GenSym.list) &&
        (p.dset == null || p.dset_n == 0)) {
        no_data_error(p);
        return;
    }

    if (p.flags & GenFlags.discard) {
        /* doing a simple "eval" */
        p.point = s;
        return;
    }

    /* extract LHS expression and operator, and test for a declaration */
    p.err = extract_lhs_and_op(&s, p, opstr.ptr);
    if (p.err || (p.flags & (GenFlags.decl | GenFlags.discard))) {
        return;
    }

    /* record next read position */
    p.point = s;

    if (p.lh.expr != null) {
        /* create syntax tree for the LHS expression */
        const(char) * savepoint = p.point;

        p.point = p.lh.expr;
        p.ch = parser_getc(p);
        lex(p);
        p.lhtree = expr(p);
        p.point = savepoint;
        p.ch = 0;
        if (!p.err) {
            p.err = compound_const_check(p.lhtree, p);
        }
        if (p.err) {
            return;
        } else {
            goto get_rhs;
        }
    }

    /* find out if the LHS var already exists, and if
       so, what type it is */
    if (!p.err) {
        p.err = check_existing_lhs_type(p, &newvar);
    }


    if (p.err) {
        return;
    }

    if (newvar) {
        /* new variable: check name for legality */
        if (!(flags & GenFlags.priv)) {
            p.err = check_identifier(p.lh.name.ptr);
        }
    } else {
        /* pre-existing var: check for const-ness */
        p.err = overwrite_const_check(p.lh.name.ptr, p.lh.vnum);
    }

    if (p.err) {
        return;
    }

    if (p.lh.t != 0) {
        if (p.targ == GenSym.unk) {
            /* when a result type is not specified, set this
               from existing LHS variable, if present
            */
            p.targ = p.lh.t;
        } else if (overwrite_type_check(p)) {
            /* don't overwrite one type with another */
            p.err = GretlError.types;
            return;
        }
    }

 get_rhs:

    /* advance past white space */
    while (isspace(*s)) s++;
    p.point = p.rhs = s;

    if (p.lh.expr != null) {
        goto alt_set_targ;
    }

    /* expression ends here with no operator: a call to print? */
    if (*s == '\0' && p.op == 0) {
        parser_try_print(p, p.lh.name.ptr, done);
        return;
    }

    /* if the LHS variable does not already exist, then
       we can't do '+=' or anything of that sort, only
       simple assignment, B_ASN
    */
    if (newvar && p.op != GenSym.b_asn) {
        undefined_symbol_error(p.lh.name.ptr, p);
        return;
    }

    if (p.op) {
        p.err = check_operator_validity(p, opstr.ptr);
        if (p.err) {
            return;
        }
    }

 alt_set_targ:

    if (p.targ == GenSym.unk && *p.rhs == '{') {
        /* if the target type is still unknown and the RHS
           expression is wrapped in '{' and '}', make the target
           a matrix
        */
        p.targ = GenSym.mat;
    } else if (p.targ == GenSym.list) {
        /* flag presence of list target to parser */
        p.flags |= GenFlags.listdef;
    }
}

void preprocess_minimal(Parser * p)
{
    const(char) * s = p.input;

    while (isspace(*s)) s++;
    if (p.targ == GenSym.unk && *s == '{') {
        p.targ = GenSym.mat;
    }
}

/* tests for saving variable */

int matrix_may_be_masked(const(GretlMatrixRaw) * m, int n, Parser * p)
{
    int mt1 = gretl_matrix_get_t1(m);
    int mt2 = gretl_matrix_get_t2(m);
    int fullrows = mt2 - mt1 + 1;
    int nobs = get_matrix_mask_nobs();

    if (n == nobs && fullrows > n) {
        p.flags |= GenFlags.mmask;
        return 1;
    } else {
        return 0;
    }
}

/* check whether a matrix result can be assigned to a series
   on return */

int series_compatible_matrix(const(GretlMatrixRaw) * m, Parser * p)
{
    int n = gretl_vector_get_length(m);
    int mt2 = gretl_matrix_get_t2(m);
    int T = sample_size(p.dset);
    int ok = 0;

    if (mt2 > 0) {
        int mt1 = gretl_matrix_get_t1(m);

        if (n == mt2 - mt1 + 1) {
            /* sample is recorded on matrix */
            ok = 1;
        } else if (matrix_may_be_masked(m, n, p)) {
            ok = 1;
        }
    } else if (n == T) {
        /* length matches current sample */
        ok = 1;
    } else if (n == p.dset.n) {
        /* length matches full series length */
        ok = 1;
    } else if (n == 1) {
        /* scalar: can be expanded */
        ok = 1;
    }

    return ok;
}

/* Check for assignment of an array of strings to a series, which
   thereby becomes string-valued. The size of the array must equal the
   full length of the dataset, which must not be sub-sampled, and if
   there's a pre-existing series it must not be string-valued already.
*/

int series_compatible_array(GretlArray * a, Parser * p)
{
    if (gretl_array_get_type(a) != GretlType.strings) {
        return 0;
    } else if (gretl_array_get_length(a) != p.dset.n) {
        return 0;
    } else if (dataset_is_subsampled(p.dset)) {
        return 0;
    } else if (is_string_valued(p.dset, p.lh.vnum)) {
        return 0;
    } else {
        /* all tests passed */
        return 1;
    }
}

/* This function converts a series to a column vector, respecting the
   value of the set-variable "skip_missing".  In that respect it differs
   from tmp_matrix_from_series(), which always just grabs the entire
   sample range.
*/

GretlMatrixRaw * series_to_matrix(const(double) * x, Parser * p)
{
    int t, n = sample_size(p.dset);
    int t1 = p.dset.t1;
    int t2 = p.dset.t2;
    int skip_na, contiguous = 1;
    GretlMatrixRaw * v;

    skip_na = libset_get_bool(SetKey.skip_missing);

    if (skip_na) {
        int err = series_adjust_sample(x, &t1, &t2);

        if (!err) {
            /* no interior NAs, just use (possibly) revised sample */
            n = t2 - t1 + 1;
        } else {
            /* we have to count the non-missing values */
            n = 0;
            for (t=t1; t<=t2; t++) {
                if (!na(x[t])) n++;
            }
            /* the values we want are not contiguous */
            contiguous = 0;
        }
    }

    if (n == 0) {
        v = gretl_null_matrix_new();
    } else {
        v = gretl_column_vector_alloc(n);
    }

    if (v == null) {
        p.err = GretlError.alloc;
    } else if (n > 0) {
        if (contiguous) {
            memcpy(v.val, x + t1, n * (*x).sizeof);
        } else {
            int i = 0;

            for (t=t1; t<=t2; t++) {
                if (na(x[t])) {
                    if (!skip_na) {
                        v.val[i++] = x[t];
                    }
                } else {
                    v.val[i++] = x[t];
                }
            }
        }
        if (contiguous) {
            gretl_matrix_set_t1(v, t1);
            gretl_matrix_set_t2(v, t2);
        }
    }

    return v;
}

GretlMatrixRaw * retrieve_matrix_result(Parser * p)
{
    Node * r = p.ret;
    GretlMatrixRaw * m = null;


    if (r.t == GenSym.num) {
        m = gretl_matrix_from_scalar(r.v.xval);
        if (m == null) {
            p.err = GretlError.alloc;
        } else if (na(r.v.xval)) {
            set_gretl_warning(GretlWarningCodes.gennan);
        }
    } else if (r.t == GenSym.series) {
        m = series_to_matrix(r.v.xvec, p);
    } else if (r.t == GenSym.list) {
        m = gretl_list_to_vector(r.v.ivec, &p.err);
    } else if (r.t == GenSym.mat && is_tmp_node(r)) {
        /* result matrix is newly allocated, steal it */
        m = r.v.m;
        r.v.m = null; /* avoid double-freeing */
    } else if (r.t == GenSym.mat) {
        /* r->v.m is an existing user matrix (or bundled matrix):
           must make a copy to keep pointers distinct
        */
        m = gretl_matrix_copy(r.v.m);
        if (m == null) {
            p.err = GretlError.alloc;
        }
    } else if (null_node(r)) {
        m = gretl_null_matrix_new();
    } else {
        fprintf(stderr, "Looking for matrix, but r->t = %s\n".ptr, getsymb(r.t));
        p.err = GretlError.types;
    }

    return m;
}

/* Check to see if the existing LHS matrix is of the same dimensions as
   the RHS result */

int LHS_matrix_reusable(Parser * p, GretlMatrixRaw ** pm, GretlMatrixRaw * tmp)
{
    GretlMatrixRaw * m = cast(GretlMatrixRaw *) gen_get_lhs_var(p, GretlType.matrix);
    int ok = 0;

    if (m == null) {
        return 0;
    } else if (p.ret.t == GenSym.num) {
        ok = (m.rows == 1 && m.cols == 1);
    } else if (p.ret.t == GenSym.series) {
        ok = (m.rows == tmp.rows && m.cols == 1);
    } else if (p.ret.t == GenSym.mat) {
        GretlMatrixRaw * retm = p.ret.v.m;

        ok = (retm != null &&
              m.rows == retm.rows &&
              m.cols == retm.cols &&
              m.is_complex == retm.is_complex);
    }

    *pm = m;

    return ok;
}

/* Generating a matrix, and there's a pre-existing LHS matrix: we re-use
   the left-hand side matrix if possible.
*/

GretlMatrixRaw * assign_to_matrix(Parser * p)
{
    GretlMatrixRaw * m = null;
    GretlMatrixRaw * tmp = null;
    int free_tmp = 1;
    double x;

    if (p.ret.t == GenSym.series) {
        /* a legacy thing */
        tmp = series_to_matrix(p.ret.v.xvec, p);
        if (p.err) {
            return null;
        }
    }

    if (LHS_matrix_reusable(p, &m, tmp)) {
        /* The result is of the same dimensions as the LHS matrix:
           this means that we don't need to construct an RHS
           matrix if it doesn't already exist as such, nor do we
           need to copy it if it does already exist.
        */
        if (p.ret.t == GenSym.num) {
            /* using RHS scalar */
            m.val[0] = x = p.ret.v.xval;
        } else if (p.ret.t == GenSym.series) {
            /* using RHS series, converted to @tmp */
            p.err = gretl_matrix_copy_data(m, tmp);
        } else {
            /* using RHS matrix: just copy data across */
            p.err = gretl_matrix_copy_data(m, p.ret.v.m);
        }
    } else {
        /* Dimensions differ: replace the LHS matrix */
        if (tmp != null) {
            p.err = gen_replace_lhs(p, GretlType.matrix, tmp);
            free_tmp = 0; /* @tmp is the return value */
        } else {
            m = retrieve_matrix_result(p);
            if (!p.err) {
                p.err = gen_replace_lhs(p, GretlType.matrix, m);
            }
        }
    }

    if (tmp != null && free_tmp) {
        gretl_matrix_free(tmp);
    }

    return m;
}

/* Assigning to an existing (whole) LHS matrix, but using '+=' or some
   such modified/inflected assignment. Note that save_generated_var() is
   the only caller.
*/

GretlMatrixRaw * assign_to_matrix_mod(GretlMatrixRaw * m1, Parser * p)
{
    GretlMatrixRaw * m2 = null;
    int mcat;

    if (m1 == null) {
        p.err = GretlError.data;
    }

    /* In most cases we can take a shortcut when the RHS
       value is scalar, but we can't do that when the
       inflection is one of the matrix concatenation
       operators: here we record that fact.
    */
    mcat = (p.op == GenSym.b_hcat || p.op == GenSym.b_vcat);

    if (!p.err) {
        if (p.op == GenSym.b_dotasn) {
            p.err = dot_assign_to_matrix(m1, p);
            m2 = m1; /* no change in matrix pointer */
        } else if (!mcat && scalar_node(p.ret)) {
            double x = node_get_scalar(p.ret, p);

            if (m1.is_complex) {
                cmatrix_xy_calc(m1, m1, Cmplx(x, 0), 0, p.op, p);
            } else {
                rmatrix_xy_calc(m1, m1, x, 0, p.op, p);
            }
            m2 = m1; /* no change in matrix pointer */
        } else if (!mcat && cscalar_node(p.ret)) {
            if (m1.is_complex) {
                Cmplx z = p.ret.v.m.z[0];

                cmatrix_xy_calc(m1, m1, z, 0, p.op, p);
                m2 = m1; /* no change in matrix pointer */
            } else {
                p.err = GretlError.types;
            }
        } else {
            GretlMatrixRaw * tmp = retrieve_matrix_result(p);

            if (tmp != null) {
                p.err = real_matrix_calc(m1, tmp, p.op, &m2);
                gretl_matrix_free(tmp);
            }
        }
    }

    return m2;
}

int should_copy_into_array(Node * n, Parser * p)
{
    int donate = is_tmp_node(n);

    /* In general if @n is a tmp node we can 'donate' its content
       to an array rather than copying it in, but there's a tricky
       case: we're compiling a "genr", or running a compiled one,
       and @n holds a string literal. If we don't copy this, it
       will have disappeared by the next call to the genr. There
       may be a simpler and more elegant way of identifying this
       case, but for now the following works (2022-08-10).
    */
    if (donate && !is_aux_node(n) && n.t == GenSym.str &&
        (p.flags & (GenFlags.compile | GenFlags.exec))) {
        donate = 0; /* we _should_ copy! */
    }

    return !donate;
}

void do_array_append(Parser * p)
{
    GretlArray * A = null;
    GretlType atype;
    Node * rhs = p.ret;
    void * ptr = null;

    A = cast(GretlArray *) gen_get_lhs_var(p, GretlType.array);
    if (A == null) {
        p.err = GretlError.data;
        return;
    }

    atype = gretl_array_get_content_type(A);

    if (atype == GretlType.string_ && rhs.t == GenSym.str) {
        ptr = rhs.v.str;
    } else if (atype == GretlType.matrix && rhs.t == GenSym.mat) {
        ptr = rhs.v.m;
    } else if (atype == GretlType.bundle && rhs.t == GenSym.bundle) {
        ptr = rhs.v.b;
    } else if (atype == GretlType.list && rhs.t == GenSym.list) {
        ptr = rhs.v.ivec;
    } else if (atype == GretlType.array && rhs.t == GenSym.array) {
        ptr = rhs.v.a;
    } else if (rhs.t == GenSym.array) {
        /* special: not actually appending an _element_;
           stick rhs array onto end of lhs array
        */
        p.err = gretl_array_copy_into(A, rhs.v.a);
    } else {
        p.err = GretlError.types;
    }


    if (!p.err && ptr != null) {
        int copy = should_copy_into_array(rhs, p);

        p.err = gretl_array_append_object(A, ptr, copy);
        if (!copy && !p.err) {
            rhs.v.ptr = null;
        }
    }
}

void do_array_subtract(Parser * p)
{
    GretlArray * A;
    Node * rhs = p.ret;

    A = cast(GretlArray *) gen_get_lhs_var(p, GretlType.array);
    if (A == null) {
        p.err = GretlError.data;
    } else if (gretl_array_get_type(A) == GretlType.strings && rhs.t == GenSym.str) {
        p.err = gretl_array_drop_string(A, rhs.v.str);
    } else if (rhs.t == GenSym.empty) {
        p.err = gretl_array_drop_null(A);
    } else {
        p.err = GretlError.types;
    }
}

void do_bundle_append(Parser * p)
{
    GretlBundle * bl = null;
    GretlBundle * br = null;
    Node * rhs = p.ret;

    bl = cast(GretlBundle *) gen_get_lhs_var(p, GretlType.bundle);
    if (rhs.t == GenSym.bundle) {
        br = rhs.v.b;
    }
    if (bl == null || br == null) {
        p.err = GretlError.types;
    } else {
        p.err = gretl_bundle_append(bl, br);
    }
}

int create_or_edit_string(Parser * p)
{
    const(char) * src = null;
    const(char) * orig = null;
    char * newstr = null;
    UserVar * uvar;

    if (p.ret.t == GenSym.num) {
        /* OK only in case of "+=" */
        if (p.op != GenSym.b_add) {
            p.err = GretlError.types;
            return p.err;
        }
    } else if (p.ret.t == GenSym.list) {
        /* string representation of list: simple assign only */
        if (p.op != GenSym.b_asn) {
            p.err = GretlError.types;
        } else {
            newstr = gretl_list_to_compact_string(p.ret.v.ivec,
                                                  p.dset, 0, &p.err);
        }
        if (p.err) {
            return p.err;
        }
        src = g_strchug(newstr);
    } else if (null_node(p.ret) || p.ret.v.str == null) {
        src = "";
    } else {
        src = p.ret.v.str;
    }


    uvar = p.lh.uv;

    if (uvar != null) {
        orig = cast(const(char) *) uvar.ptr;
    } else if (p.op != GenSym.b_asn) {
        /* without an existing LHS string we can only assign */
        p.err = GretlError.data;
        return p.err;
    }

    if (p.ret.t == GenSym.num) {
        /* taking an offset into an existing string */
        int len = cast(int) g_utf8_strlen(orig, -1);
        int adj = cast(int) p.ret.v.xval;

        if (adj < 0) {
            p.err = GretlError.data;
        } else if (adj == 0) {
            ; /* no-op */
        } else {
            if (adj < len) {
                src = g_utf8_offset_to_pointer(orig, adj);
            } else {
                src = "";
            }
            newstr = gretl_strdup(src);
            if (newstr == null) {
                p.err = GretlError.alloc;
            } else {
                gen_replace_lhs(p, GretlType.string_, newstr);
            }
        }
    } else if (src == null) {
        ; /* no-op -- e.g. argname() didn't get anything */
    } else if (p.op == GenSym.b_asn) {
        /* simple assignment */
        if (newstr == null) {
            newstr = gretl_strdup(src);
            if (newstr == null) {
                p.err = GretlError.alloc;
            }
        }
        if (!p.err) {
            if (uvar == null) {
                gen_add_uvar(p, GretlType.string_, newstr);
            } else {
                gen_replace_lhs(p, GretlType.string_, newstr);
            }
        }
    } else if (p.op == GenSym.b_hcat || p.op == GenSym.b_add) {
        /* string concatenation */
        if (*src == '\0') {
            ; /* no-op */
        } else {
            newstr = cast(char *) malloc(strlen(orig) + strlen(src) + 1);
            if (newstr == null) {
                p.err = GretlError.alloc;
            } else {
                strcpy(newstr, orig);
                strcat(newstr, src);
                gen_replace_lhs(p, GretlType.string_, newstr);
            }
        }
    }

    return p.err;
}

int create_or_edit_list(Parser * p)
{
    int * list = null;

    if (p.ret.t == GenSym.mat && gretl_vector_get_length(p.ret.v.m) == 0) {
        /* special case, list from matrix */
        const(char) * prefix;

        prefix = p.ret.vname != null ? p.ret.vname : p.lh.name.ptr;
        list = gretl_list_from_matrix(p.ret.v.m, prefix,
                                      p.dset, &p.err);
    } else {
        list = node_get_list(p.ret, p); /* note: copied */
    }



    if (!p.err) {
        if (p.lh.t != GenSym.list) {
            /* no pre-existing LHS list: must be simple assignment */
            p.err = remember_list(list, p.lh.name.ptr, null);
        } else if (p.op == GenSym.b_asn || p.op == GenSym.b_add || p.op == GenSym.b_sub) {
            /* replace, append or subtract list members */
            p.err = gen_edit_list(p, list, p.op);
        } else {
            p.err = GretlError.types;
        }
    }


    free(list);

    return p.err;
}

bool ok_return_type(int t)
{
    return (t == GenSym.num || t == GenSym.series || t == GenSym.mat ||
           t == GenSym.list || t == GenSym.dum || t == GenSym.empty ||
           t == GenSym.str || t == GenSym.bundle || t == GenSym.array ||
           t == GenSym.u_addr || t == GenSym.dbundle);
}

/* Note: we're doing this only in relation to "primary" types (excluding
   bundle members, array elements, matrix sub-specs).
*/

int gen_check_return_type(Parser * p)
{
    Node * r = p.ret;
    int err = 0;

    if (r == null) {
        fprintf(stderr, "gen_check_return_type: p->ret = NULL!\n".ptr);
        return GretlError.data;
    }


    if (!ok_return_type(r.t)) {
        return GretlError.types;
    }

    if (r.t == GenSym.series && r.v.xvec == null) {
        fprintf(stderr, "got SERIES return with xvec = NULL!\n".ptr);
        return GretlError.data;
    }

    if (p.targ == GenSym.num) {
        if (r.t == GenSym.num || scalar_matrix_node(r)) {
            ; /* scalar or 1 x 1 matrix: OK */
        } else if (r.t == GenSym.mat && (p.flags & GenFlags.nodecl)) {
            ; /* morphing to matrix may be OK */
        } else {
            err = GretlError.types;
        }
    } else if (p.targ == GenSym.series) {
        /* result must be scalar, series, or conformable matrix
           or array of strings
        */
        if (r.t == GenSym.num || r.t == GenSym.series) {
            ; /* OK */
        } else if (r.t == GenSym.mat) {
            if (!series_compatible_matrix(r.v.m, p)) {
                err = GretlError.types;
            }
        } else if (r.t == GenSym.array) {
            if (p.op != GenSym.b_asn || !series_compatible_array(r.v.a, p)) {
                err = GretlError.types;
            }
        } else {
            err = GretlError.types;
        }
    } else if (p.targ == GenSym.mat) {
        ; /* no-op: handled later */
    } else if (p.targ == GenSym.list) {
        if (r.t != GenSym.empty && r.t != GenSym.mat && !ok_list_node(r, p)) {
            err = GretlError.types;
        }
    } else if (p.targ == GenSym.str) {
        if (r.t != GenSym.empty && r.t != GenSym.str && r.t != GenSym.num && r.t != GenSym.list) {
            err = GretlError.types;
        }
    } else if (p.targ == GenSym.bundle) {
        if (p.op == GenSym.b_asn) {
            /* plain assignment: bundle or null */
            if (r.t != GenSym.bundle && r.t != GenSym.dbundle && r.t != GenSym.empty) {
                err = GretlError.types;
            }
        } else {
            /* the only other assignment possibility is "+=",
               in which case we'll only accept a bundle
            */
            if (r.t != GenSym.bundle) {
                err = GretlError.types;
            }
        }
    } else if (p.targ == GenSym.array) {
        if (p.op == GenSym.b_asn) {
            /* plain assignment: array, null, or suitable series */
            if (gen_type_is_arrayable(r.t) || r.t == GenSym.empty) {
                ; /* OK */
            } else if (strvals_node(r, p)) {
                ; /* possibly OK */
            } else {
                err = GretlError.types;
            }
        } else {
            /* arrays: the only other assignment possibility is "+=",
               in which case we'll only accept an array or an
               object which matches the content type on the left
               (but the matching check is deferred)
            */
            if (!gen_type_is_arrayable(r.t) && r.t != GenSym.empty) {
                err = GretlError.types;
            }
        }
    }

    if (err == GretlError.types) {
        assignment_type_errmsg(p.targ, r.t, p.op);
    }


    return err;
}

/* Allocate storage if saving a new series to the dataset. We also check
   for collision with the name of a function and issue a warning if need
   be.
*/

int gen_allocate_storage(Parser * p)
{
    if (p.dset == null || p.dset.Z == null) {
        p.err = GretlError.data;
    } else {
        p.err = dataset_add_NA_series(p.dset, 1);
        if (!p.err) {
            p.lh.vnum = p.dset.v - 1;
        }
    }
    if (!p.err && gretl_function_depth() == 0 &&
        get_user_function_by_name(p.lh.name.ptr) != null) {
        gretl_warnmsg_sprintf(_("'%s' shadows a function of the same name").ptr,
                              p.lh.name.ptr);
    } else if (!p.err && function_lookup(p.lh.name.ptr)) {
        gretl_warnmsg_sprintf(_("'%s' shadows a function of the same name").ptr,
                              p.lh.name.ptr);
    }

    return p.err;
}

void align_matrix_to_series(double * y, const(GretlMatrixRaw) * m, Parser * p)
{
    const(GretlMatrixRaw) * mask = get_matrix_mask();
    int t, s = 0;

    if (mask == null || mask.rows != p.dset.n) {
        p.err = GretlError.data;
        return;
    }

    for (t=0; t<p.dset.n; t++) {
        if (mask.val[t] != 0.0) {
            if (t >= p.dset.t1 && t <= p.dset.t2) {
                y[t] = xy_calc(y[t], m.val[s], p.op, GenSym.series, p);
            }
            s++;
        }
    }
}

int assign_null_to_bundle(Parser * p)
{
    GretlBundle * b;
    int err = 0;

    if (p.lh.t == GenSym.bundle) {
        b = cast(GretlBundle *) gen_get_lhs_var(p, GretlType.bundle);
        gretl_bundle_void_content(b);
    } else {
        b = gretl_bundle_new();
        if (b == null) {
            err = GretlError.alloc;
        } else {
            err = gen_add_uvar(p, GretlType.bundle, b);
        }
    }

    return err;
}

int assign_null_to_array(Parser * p)
{
    GretlArray * a;
    int err = 0;

    if (p.lh.t == GenSym.array) {
        a = cast(GretlArray *) gen_get_lhs_var(p, GretlType.array);
        gretl_array_void_content(a);
    } else {
        a = gretl_array_new(p.lh.gtype, 0, &err);
        if (!err) {
            err = gen_add_uvar(p, p.lh.gtype, a);
        }
    }

    return err;
}

double get_lhs_scalar(Parser * p)
{
    if (p.lh.uv == null) {
        p.lh.uv = get_user_var_by_name(p.lh.name.ptr);
    }
    if (p.lh.uv != null) {
        return uvar_get_scalar_value(p.lh.uv);
    } else {
        p.err = GretlError.data;
        return NADBL;
    }
}

/* apply postfix '++' or '--' to LHS scalar, or '++' to
   LHS string (only) */

int do_incr_decr(Parser * p)
{
    if (p.lh.uv != null && p.lh.uv.type == GretlType.double_) {
        double x = get_lhs_scalar(p);

        if (!na(x)) {
            x += (p.op == GenSym.inc)? 1 : -1;
            uvar_set_scalar_fast(p.lh.uv, x);
        }
    } else if (p.lh.uv != null && p.lh.uv.type == GretlType.string_) {
        if (p.op == GenSym.dec) {
            p.err = GretlError.types;
        } else {
            char * s = cast(char *) p.lh.uv.ptr;

            if (*s != '\0') {
                char * smod = gretl_strdup(s + 1);

                gen_replace_lhs(p, GretlType.string_, smod);
            }
        }
    } else {
        p.err = GretlError.types;
    }

    return p.err;
}

bool has_aux_mat(Node * n)
{
    return (n.aux != null && n.aux.t == GenSym.mat);
}

int explore_node(Node * t, int lev, Node * prev, Parser * p)
{
    Node pms;
    int save_op;
    int err = 0;

    if (prev != null && (t.t == GenSym.mat || has_aux_mat(t))) {
        pms.t = GenSym.osl;
        /* pms.L: node holding target matrix */
        pms.L = t.t == GenSym.mat ? t : t.aux;
        /* pms.R: node holding mspec */
        pms.R = prev.R.aux;
        save_op = p.op;
        p.op = GenSym.b_asn;
        /* prev->aux holds the replacement matrix */
        err = set_matrix_chunk(&pms, prev.aux, p);
        p.op = save_op;
    }

    return err;
}

int traverse_left(Parser * p)
{
    Node * t = p.lhtree;
    Node * prev = null;
    int level = 0;
    int err = 0;

    while (t && !err) {
        err = explore_node(t, level, prev, p);
        if (t.aux != null && t.aux.t == GenSym.mat) {
            prev = t;
        } else {
            prev = null;
        }
        t = t.L;
        level++;
    }

    return err;
}

/* set_nested_matrix_value(): this and its helper above,
   traverse_left(), require a little comment.

   We come here when a hansl statement modifies a matrix that is "under"
   something else. That something could be a bundle or array, as in

   # Case 0
   b.m[diag] = x     # under bundle b
   a[3][1:2,1:2] = y # under array a

   In such cases the first invocation of set_matrix_chunk below is
   sufficient. However, we also come here when the matrix is "under"
   another matrix -- that is, we have a double index or subspec, as in
   these examples for a complex matrix, C:

   # Case 1
   C[real][1:2,1:2] = x
   C[3,3][real] = y

   and also in these more extended examples where the complex matrix is
   itself "under" something else:

   # Case 2
   b.C[real][1:2,1:2] = x
   a[3][i,j][imag] = y

   To handle such cases we have to crawl the parser's "lhtree"
   (left-hand side tree) to find the matrix that ultimately has to be
   modified, executing further calls to set_matrix_chunk(). In Case 1
   above, the matrix we're looking for will be at depth 1 in the lhtree,
   while in Case 2 it will be at depth 2.
*/

int set_nested_matrix_value(Node * lhs, Node * rhs, Parser * p)
{
    int err = set_matrix_chunk(lhs, rhs, p);

    if (!err) {
        err = traverse_left(p);
    }

    return err;
}

/* Support direct assignment from an array of strings to a series
   (which will become string-valued). The source array must hold
   a number of strings equal to the length of the dataset.
*/

int series_from_strings_array(DatasetRaw * dset, int v, GretlArray * a)
{
    char ** S;
    int ns;

    S = gretl_array_get_strings(a, &ns);
    return series_from_strings(dset, v, S, ns);
}

/* We come here is the genr target is SERIES but we got a matrix on the
   right-hand side. @y is the storage area for the series and @m is the
   source matrix.
*/

void series_from_matrix(double * y, const(GretlMatrixRaw) * m, Parser * p)
{
    int k = gretl_vector_get_length(m);
    int mt1 = gretl_matrix_get_t1(m);
    int mt2 = gretl_matrix_get_t2(m);
    int s, t;

    if (p.flags & GenFlags.mmask) {
        /* result needs special alignment */
        align_matrix_to_series(y, m, p);
    } else if (k == 1) {
        /* the result is effectively a scalar */
        for (t=p.dset.t1; t<=p.dset.t2; t++) {
            y[t] = xy_calc(y[t], m.val[0], p.op, GenSym.series, p);
        }
    } else if (k == p.dset.n) {
        /* treat the result as a full-length series */
        for (t=p.dset.t1; t<=p.dset.t2; t++) {
            y[t] = xy_calc(y[t], m.val[t], p.op, GenSym.series, p);
        }
    } else if (k == sample_size(p.dset)) {
        /* treat as a series of the current sample length */
        for (t=p.dset.t1, s=0; t<=p.dset.t2; t++, s++) {
            y[t] = xy_calc(y[t], m.val[s], p.op, GenSym.series, p);
        }
    } else if (mt1 >= 0 && mt2 >= mt1) {
        /* align using the matrix "t1" value */
        for (t=mt1; t<mt1 + k && t<=p.dset.t2; t++) {
            if (t >= p.dset.t1) {
                y[t] = xy_calc(y[t], m.val[t-mt1], p.op,
                   GenSym.series, p);
            }
        }
    }
}

enum {
    OVW_STRINGS = 1,
    OVW_NUMERIC
};

/* Determine if it seems OK to overwrite a string-valued series.  This
   action was banned altogether up to gretl 2023b, but we're now
   allowing it subject to two conditions:

   * It must be a case of simple (not inflected) assignment
   * The RHS must be a series

   If the RHS series is not string-valued we'll treat that as OK
   provisionally. The suitability of the RHS values will be checked
   later, when we call assign_numeric_to_strvar(): they must be
   compatible with the encoding on the LHS.

   If the case is not rejected as unsuitable, we return either
   OVW_STRINGS to indicate a string-valued RHS or OVW_NUMERIC to
   indicate an overwrite using plain numerical input.
*/

int strv_overwrite_ok(Parser * p)
{
    int lv = p.lh.vnum; /* LHS series ID */
    int rv = p.ret.vnum; /* RHS series ID */
    int ok = 1;
    int ret = 0;

    if (p.op != GenSym.b_asn) {
        /* must be a case of straight assignment */
        ok = 0;
    } else if (p.ret.t != GenSym.series) {
        /* RHS must be a series */
        ok = 0;
    } else if (lv == rv) {
        /* and not identical to LHS */
        ok = 0;
    } else if (is_string_valued(p.dset, rv)) {
        ret = OVW_STRINGS;
    } else if (p.lh.stab != null) {
        ret = OVW_STRINGS;
    } else {
        ret = OVW_NUMERIC;
    }

    if (!ok) {
        gretl_errmsg_set("Cannot overwrite string-valued series in this context");
        p.err = GretlError.types;
    }

    return ret;
}

int savegen_retval(int err)
{
    return err;
}


void handle_one_by_one_cast(Parser * p, Node * r, int * no_decl)
{
    if (p.targ == GenSym.unk) {
        if (scalar_matrix_node(r)) {
            /* "cast" 1 x 1 matrix to scalar */
            *no_decl = 1;
            p.targ = GenSym.num;
            p.flags |= GenFlags.nodecl;
        } else {
            p.targ = r.t;
        }
    } else if (p.targ == GenSym.num && r.t == GenSym.mat && (p.flags & GenFlags.nodecl)) {
        /* We're looking at a @targ that was previously
           set to NUM by the "auto-cast" mechanism: allow
           it to morph to matrix if need be.
        */
        if (scalar_matrix_node(r)) {
            ; /* not a problem */
        } else if (p.lh.t == 0) {
            /* no pre-existing scalar var */
            p.targ = GenSym.mat;
        } else if (p.lh.t == GenSym.num) {
            /* type-convert existing scalar */
            p.err = gretl_scalar_convert_to_matrix(p.lh.uv);
            if (!p.err) {
                p.targ = GenSym.mat;
            }
        }
    }
}


int handle_compound_target(Parser * p)
{
    Node * r = p.ret;
    int compound_t;

    p.lhtree.flags |= NodeFlags.lht;
    p.flags |= GenFlags.start;
    p.lhres = eval(p.lhtree, p);
    if (p.err) {
        return savegen_retval(p.err);
    }

    compound_t = p.lhres.t;
    if (compound_t == GenSym.bmemb) {
        p.err = set_bundle_value(p.lhres, r, p);
    } else if (compound_t == GenSym.obs) {
        p.err = set_series_obs_value(p.lhres, r, p);
    } else if (compound_t == GenSym.osl) {
        Node * lh1 = p.lhres.L;

        if (p.targ == GenSym.series && lh1.t != GenSym.list && lh1.t != GenSym.bundle) {
            /* a sub-object can be of series type only if the "parent"
               object is a list or bundle
            */
            p.err = GretlError.types;
        } else if (lh1.t == GenSym.array) {
            p.err = set_array_value(p.lhres, r, p);
        } else if (lh1.t == GenSym.list) {
            p.err = set_list_value(p.lhres, r, p);
        } else if (lh1.t == GenSym.str) {
            p.err = set_string_value(p.lhres, r, p);
        } else if (lh1.t == GenSym.bundle) {
            p.err = set_bundle_value(p.lhres, r, p);
        } else if (lh1.t == GenSym.series) {
            p.err = set_series_obs_value(p.lhres, r, p);
        } else if (lh1.t == GenSym.mat) {
            p.err = set_nested_matrix_value(p.lhres, r, p);
        } else {
            gretl_errmsg_set(_("Invalid left-hand side expression").ptr);
            p.err = GretlError.types;
        }
    } else {
        gretl_errmsg_set(_("Invalid left-hand side expression").ptr);
        p.err = GretlError.types;
    }

    return savegen_retval(p.err);
}

int handle_overwrite_strvals(Parser * p, int targ, Node * r)
{
    SeriesTable * st;
    int copy = 0;
    int err;

    if (r.vnum == NO_VNUM) {
        /* we're assigning from an on-the-fly series, whose
           temporary series_table is attached to @p->lh
        */
        st = p.lh.stab;
    } else {
        /* assigning from a string-valued member of the dataset,
           whose series_table must be preserved as is
        */
        st = series_get_string_table(p.dset, r.vnum);
        copy = 1;
    }

    err = assign_strings_to_strvar(p.dset, targ, r.v.xvec, st, copy);

    if (st == p.lh.stab) {
        /* @st has done its work, and must now be set to NULL */
        if (dataset_is_subsampled(p.dset)) {
            /* in the subsampled case @st is not "donated"
               to the LHS series, so we have to free it here
               to avoid leaking memory
            */
            series_table_destroy(p.lh.stab);
        }
        p.lh.stab = null;
    }

    return err;
}

int save_generated_var(Parser * p, Prn * prn)
{
    Node * r = p.ret;
    double ** Z = null;
    double x;
    int no_decl = 0;
    int strv_ovwrite = 0;
    int t, v = 0;


    if (p.lh.t == GenSym.series && is_string_valued(p.dset, p.lh.vnum) &&
        p.lhtree == null) {
        strv_ovwrite = strv_overwrite_ok(p);
        if (!strv_ovwrite) {
            return savegen_retval(p.err);
        }
    }

    /* some special cases */
    if (p.lhtree != null) {
        return handle_compound_target(p);
    } else if (p.op == GenSym.inc || p.op == GenSym.dec) {
        return do_incr_decr(p);
    }

    if (p.callcount < 2) {
        /* first exec: test for type mismatch errors */
        p.err = gen_check_return_type(p);
        if (p.err) {
            return savegen_retval(p.err);
        }
    }

    handle_one_by_one_cast(p, r, &no_decl);


    if (p.targ == GenSym.unk || p.targ == GenSym.empty) {
        p.err = GretlError.unkvar;
        return savegen_retval(p.err);
    }

    if (p.targ == GenSym.series && cast(ubyte) p.lh.name[0] > 126) {
        /* can't allow Greek letters for series names */
        gretl_errmsg_sprintf("Invalid series name '%s'".ptr, p.lh.name.ptr);
        p.err = GretlError.data;
        return savegen_retval(p.err);
    }

    /* allocate dataset storage, if needed */
    if (p.targ == GenSym.series && p.lh.vnum <= 0) {
        gen_allocate_storage(p);
        if (p.err) {
            return savegen_retval(p.err);
        }
    }

    if (p.dset != null && p.dset.Z != null) {
        /* convenience notation */
        Z = p.dset.Z;
        v = p.lh.vnum;
    }

    /* put the generated data into place */

    if (p.targ == GenSym.num) {
        if (p.lh.t == GenSym.num) {
            /* modifying an existing scalar */
            if (r.t == GenSym.num) {
                x = r.v.xval;
            } else if (scalar_matrix_node(r)) {
                x = r.v.m.val[0];
            } else {
                p.err = GretlError.types;
            }
            if (!p.err && p.op != GenSym.b_asn) {
                double x0 = get_lhs_scalar(p);

                x = xy_calc(x0, x, p.op, GenSym.num, p);
            }
            if (!p.err) {
                p.err = gen_replace_scalar(p, x);
            }
        } else {
            /* a new scalar */
            if (r.t == GenSym.num) {
                x = r.v.xval;
            } else if (scalar_matrix_node(r)) {
                x = r.v.m.val[0];
            } else {
                p.err = GretlError.types;
            }
            if (!p.err) {
                if (no_decl) {
                    p.err = gretl_scalar_add_mutable(p.lh.name.ptr, x);
                } else {
                    p.err = gretl_scalar_add(p.lh.name.ptr, x);
                }
            }
        }
    } else if (p.targ == GenSym.series) {
        /* writing a series */
        if (strv_ovwrite == OVW_STRINGS) {
            p.err = handle_overwrite_strvals(p, v, r);
        } else if (strv_ovwrite == OVW_NUMERIC) {
            p.err = assign_numeric_to_strvar(p.dset, v, r.v.xvec);
        } else if (r.t == GenSym.series) {
            const(double) * xs = r.v.xvec;

            if (p.op == GenSym.b_asn) {
                /* avoid multiple calls to xy_calc */
                if (Z[v] != xs) {
                    size_t sz = sample_size(p.dset) * (*xs).sizeof;

                    memcpy(Z[v] + p.dset.t1, xs + p.dset.t1, sz);
                }
            } else {
                for (t=p.dset.t1; t<=p.dset.t2; t++) {
                    Z[v][t] = xy_calc(Z[v][t], xs[t], p.op, GenSym.series, p);
                }
            }
        } else if (r.t == GenSym.num) {
            for (t=p.dset.t1; t<=p.dset.t2; t++) {
                Z[v][t] = xy_calc(Z[v][t], r.v.xval, p.op, GenSym.series, p);
            }
        } else if (r.t == GenSym.mat) {
            series_from_matrix(Z[v], r.v.m, p);
        } else if (r.t == GenSym.array) {
            p.err = series_from_strings_array(p.dset, v, r.v.a);
            if (p.err) {
                return savegen_retval(p.err);
            }
        }
        strcpy(p.dset.varname[v], p.lh.name.ptr);
        series_unset_orig_pd(p.dset, v); /* 2020-09-27 */
        if (!p.err) {
            /* (probably) changed or added a series */
            set_dataset_is_changed(p.dset, 1);
        }
    } else if (p.targ == GenSym.mat) {
        /* we're writing a matrix */
        GretlMatrixRaw * m = null;

        if (p.lh.uv == null) {
            /* there's no pre-existing left-hand side matrix */
            m = retrieve_matrix_result(p);
            if (!p.err) {
                p.err = gen_add_uvar(p, GretlType.matrix, m);
            }
        } else if (p.op == GenSym.b_asn) {
            /* uninflected assignment to an existing matrix */
            m = assign_to_matrix(p);
        } else {
            /* inflected assignment to entire existing matrix */
            GretlMatrixRaw * m1 = cast(GretlMatrixRaw *) gen_get_lhs_var(p, GretlType.matrix);

            m = assign_to_matrix_mod(m1, p);
            if (!p.err) {
                p.err = gen_replace_lhs(p, GretlType.matrix, m);
            }
        }
        /* note: for use by genr_get_output_matrix() */
        p.lh.mret = m;
    } else if (p.targ == GenSym.list) {
        create_or_edit_list(p);
    } else if (p.targ == GenSym.str) {
        create_or_edit_string(p);
    } else if (p.targ == GenSym.bundle) {
        if (null_node(r)) {
            /* as in "bundle b = null" */
            p.err = assign_null_to_bundle(p);
        } else if (p.op != GenSym.b_asn) {
            do_bundle_append(p);
        } else {
            /* full assignment of RHS bundle */
            GretlBundle * b;

            if (r.t == GenSym.dbundle) {
                b = bvar_get_bundle(r, p);
            } else if (is_tmp_node(r) || (p.flags & GenFlags.ufret)) {
                /* grabbing r->v.b is OK */
                b = r.v.b;
            } else {
                /* we need to make a copy */
                b = gretl_bundle_copy(r.v.b, &p.err);
            }

            if (!p.err) {
                p.err = gen_add_or_replace(p, GretlType.bundle, b);
                if (!p.err && r.t != GenSym.dbundle && b == r.v.b) {
                    /* avoid destroying the assigned bundle */
                    r.v.b = null;
                }
            }
        }
    } else if (p.targ == GenSym.array) {
        if (p.op == GenSym.b_add) {
            do_array_append(p);
        } else if (p.op == GenSym.b_sub) {
            do_array_subtract(p);
        } else if (null_node(r)) {
            /* as in, e.g., "strings A = empty" */
            p.err = assign_null_to_array(p);
        } else if (r.t == GenSym.array) {
            /* full assignment of RHS array */
            GretlType atype = gretl_array_get_type(r.v.a);
            GretlArray * a = null;

            if (p.lh.gtype > 0 && atype != p.lh.gtype) {
        p.err = GretlError.types;
            } else if (is_tmp_node(r) || (p.flags & GenFlags.ufret)) {
                /* grabbing r->v.a is OK */
                a = r.v.a;
            } else {
                /* we need to make a copy */
                a = gretl_array_copy(r.v.a, &p.err);
            }
            if (!p.err) {
                p.err = gen_add_or_replace(p, atype, a);
                if (!p.err && a == r.v.a) {
                    /* avoid destroying the assigned array */
                    r.v.a = null;
                }
            }
        } else if (r.t == GenSym.series) {
            /* verified in advance as string-valued */
            GretlType atype = GretlType.strings;

            if (p.lh.gtype > 0 && p.lh.gtype != atype) {
                p.err = GretlError.types;
            } else {
                GretlArray * a = get_strings_array_from_series(p.dset, r.vnum, &p.err);

                if (!p.err) {
                    p.err = gen_add_or_replace(p, atype, a);
                }
            }
        } else {
            /* Allow promotion of a single object to an array of
               size 1? Note 2021-08-12: not sure this is actually
               a good idea.
            */
            GretlType rtype = gretl_type_from_gen_type(r.t);
            GretlType atype = p.lh.gtype;

            if (rtype == gretl_type_get_singular(atype)) {
                GretlArray * a = gretl_singleton_array(r.v.ptr, atype, 1, &p.err);

                if (!p.err) {
                    p.err = gen_add_or_replace(p, atype, a);
                }
            } else {
                p.err = e_types(r);
            }
        }
    }

    return savegen_retval(p.err);
}

void maybe_update_lhs_uvar(Parser * p, GretlType * type)
{
    if (p.targ == GenSym.series) {
        /* targetting a series */
        int v = p.lh.vnum;

        if (get_loop_renaming() || v <= 0 || v >= p.dset.v) {
            p.lh.vnum = current_series_index(p.dset, p.lh.name.ptr);
        }
        if (p.lh.vnum < 0) {
            p.lh.vnum = 0;
        }
        return;
    }

    if (p.lh.uv == null) {
        p.lh.uv = get_user_var_by_name(p.lh.name.ptr);
    }

    if (p.lh.uv != null) {
        *type = p.lh.uv.type;
    }

    switch (*type) {
    case GretlType.double_:
        p.lh.t = GenSym.num;
        break;
    case GretlType.matrix:
        p.lh.t = GenSym.mat;
        if (p.targ == GenSym.num) {
            p.targ = GenSym.mat;
        }
        break;
    case GretlType.list:
        p.lh.t = GenSym.list;
        break;
    case GretlType.string_:
        p.lh.t = GenSym.str;
        break;
    case GretlType.bundle:
        p.lh.t = GenSym.bundle;
        break;
    case GretlType.array:
        p.lh.t = GenSym.array;
        break;
    default:
        p.lh.t = 0;
        break;
    }
}

void parser_reinit(Parser * p, DatasetRaw * dset, Prn * prn)
{
    /* flags that should be reinstated if they were
       set at compile time, or in previous execution
    */
    static immutable int[] saveflags = [GenFlags.natest, GenFlags.autoreg, GenFlags.discard, GenFlags.nodecl,
        GenFlags.listdef, GenFlags.priv, GenFlags.ufret, GenFlags.altinp, 0];
    int i, prevflags = p.flags;
    GretlType lhtype = GretlType.none;
    int dset_n;

    if (p.callcount > 1) {
        /* the flags should basically have stabilized by now */
        p.flags |= GenFlags.start;
        p.flags = cast(GenFlags)(p.flags & ~GenFlags.deltan);
    } else {
        p.flags = (GenFlags.start | GenFlags.exec);
        for (i=0; saveflags[i] > 0; i++) {
            if (prevflags & saveflags[i]) {
                p.flags |= saveflags[i];
            }
        }
    }

    p.dset = dset;
    p.prn = prn;

    p.obs = 0;
    p.sym = 0;
    p.ch = 0;
    p.xval = 0.0;
    p.idnum = 0;
    p.idstr = null;
    p.data = null;
    p.errprn = null;

    p.ret = null;
    p.lhres = null;
    p.err = 0;


    if (*p.lh.name.ptr != '\0') {
        maybe_update_lhs_uvar(p, &lhtype);
    }

    /* allow for change in length of dataset */
    dset_n = dset != null ? dset.n : 0;
    if (dset_n != p.dset_n) {
        p.dset_n = dset_n;
        p.flags |= GenFlags.deltan;
    }
}

void parser_init(Parser * p, const(char) * str, DatasetRaw * dset, Prn * prn, GenFlags flags, int targtype, int * done)
{
    p.flags = flags | GenFlags.start;
    if (p.flags & GenFlags.compile) {
        char * tmp = gretl_strdup(str);

        p.point = p.rhs = p.input = tmp;
        p.flags |= GenFlags.altinp;
    } else {
        p.point = p.rhs = p.input = str;
    }
    p.dset = dset;
    p.dset_n = dset != null ? dset.n : 0;
    p.prn = prn;
    p.errprn = null;
    p.targ = targtype;
    p.op = 0;

    p.lhtree = null;
    p.lhres = null;
    p.tree = null;
    p.ret = null;

    /* left-hand side info */
    p.lh.t = 0;
    p.lh.name[0] = '\0';
    p.lh.label = null;
    p.lh.stab = null;
    p.lh.vnum = 0;
    p.lh.uv = null;
    p.lh.expr = null;
    p.lh.gtype = GretlType.none;
    p.lh.mret = null;

    /* auxiliary apparatus */
    p.aux = null;

    p.callcount = 0;
    p.obs = 0;
    p.sym = 0;
    p.upsym = 0;
    p.ch = 0;
    p.xval = 0.0;
    p.idnum = 0;
    p.idstr = null;
    p.err = 0;

    if (p.input == null) {
        p.err = GretlError.data;
        return;
    }

    if (p.flags & GenFlags.void_) {
        p.flags |= GenFlags.discard;
    } else if (p.flags & GenFlags.anon) {
        preprocess_minimal(p);
    } else {
        gen_preprocess(p, flags, done);
    }
    if (p.targ == GenSym.list) {
        p.flags |= GenFlags.listdef;
    }

    if (!p.err) {
        p.ch = parser_getc(p);
    }
}

/* called from genmain.c (only!) */

void gen_save_or_print(Parser * p, Prn * prn)
{
    if (autoreg(p)) {
        /* no transcription required */
        return;
    }
    if (p.flags & GenFlags.discard) {
        /* doing "eval" */
        if (p.ret == null) {
            return;
        } else if (p.ret.t == GenSym.empty) {
            pputs(prn, "[null]\n");
        } else if (p.ret.t == GenSym.mat) {
            if (p.ret.v.m.is_complex) {
                gretl_cmatrix_print(p.ret.v.m, p.lh.name.ptr, p.prn);
            } else {
                gretl_matrix_print_to_prn(p.ret.v.m, p.lh.name.ptr, p.prn);
            }
        } else if (p.ret.t == GenSym.list) {
            if (p.lh.name[0] != '\0') {
                gretl_list_print(get_list_by_name(p.lh.name.ptr),
                                 p.dset, p.prn);
            } else {
                gretl_list_print(p.ret.v.ivec, p.dset, p.prn);
            }
        } else if (p.ret.t == GenSym.str) {
            if (p.lh.name[0] != '\0') {
                pprintf(p.prn, "%s\n".ptr, gen_get_lhs_var(p, GretlType.string_));
            } else {
                pprintf(p.prn, "%s\n".ptr, p.ret.v.str);
            }
        } else if (p.ret.t == GenSym.bundle) {
            gretl_bundle_print(p.ret.v.b, prn);
        } else if (p.ret.t == GenSym.array) {
            gretl_array_print(p.ret.v.a, prn);
        } else {
            /* scalar, series */
            printnode(p.ret, p, 1);
            pputc(p.prn, '\n');
        }
    } else if (p.flags & GenFlags.decl) {
        do_declaration(p);
    } else {
        save_generated_var(p, prn);
    }
}

void gen_cleanup(Parser * p)
{
    int save = (p.flags & (GenFlags.compile | GenFlags.exec));


    if (p.lh.label != null) {
        free(p.lh.label);
        p.lh.label = null;
    }

    if (p.lh.stab != null) {
        series_table_destroy(p.lh.stab);
        p.lh.stab = null;
    }

    if (p.flags & GenFlags.altinp) {
        free(cast(char *) p.input);
        p.input = null;
    }

    if (p.err && (p.flags & GenFlags.compile)) {
        save = 0;
    }

    if (!save) {
        if (p.lhtree != null) {
            if (p.lhtree != p.lhres) {
                /* we have to scrub the LHT_NODE flag on p->lhtree,
                   or else its children will not get freed and we'll
                   leak memory
                */
                p.lhtree.flags = cast(ubyte)(p.lhtree.flags & ~NodeFlags.lht);
                free_tree(p.lhtree, p, FreeCode.lhtree);
            }
            if (p.lhres != null) {
                free_tree(p.lhres, p, FreeCode.lhres);
            }
        }

        if (p.tree != p.ret) {
            free_tree(p.tree, p, FreeCode.tree);
        }

        free_tree(p.ret, p, FreeCode.ret);

        free(p.lh.expr);
    }

}

void real_reset_uvars(Parser * p)
{
    if (0 && p.err) {
        /* 2022-09-13: not sure about this */
        return;
    }


    clear_uvnodes(p.tree);

    if (p.lhtree != null) {
        clear_uvnodes(p.lhtree);
    }

    p.lh.uv = null;
    p.lh.vnum = 0;
}

void genr_reset_uvars(Parser * p)
{
    real_reset_uvars(p);
}

void maybe_set_return_flags(Parser * p)
{
    Node * t = p.tree;

    if (t != null && t.t == GenSym.ufun) {
        p.flags |= GenFlags.ufret;
    }
}

int decl_check(Parser * p, int flags)
{
    if (flags & GenFlags.compile) {
        p.err = GretlError.eqn;
        gretl_errmsg_sprintf("%s:\n> '%s'".ptr,
                             _("Bare declarations are not allowed here").ptr,
                             p.input);
    }

    return p.err;
}

void autoreg_error(Parser * p, int t)
{
    fprintf(stderr, "*** autoreg error at obs t = %d (t1 = %d):\n".ptr,
            t, p.dset.t1);

    if (p.ret != null && p.ret.t != GenSym.series) {
        fprintf(stderr, " ret type != SERIES (=%d), p->err = %d\n".ptr,
                p.ret.t, p.err);
    } else if (p.ret == null) {
        fprintf(stderr, " ret = NULL, p->err = %d\n".ptr, p.err);
    }

    fprintf(stderr, " input = '%s'\n".ptr, p.input);

    if (!p.err) {
        p.err = GretlError.data;
    }
}


int realgen(const(char) * s, Parser * p, DatasetRaw * dset, Prn * prn, int flags, int targtype)
{

    if (flags & GenFlags.exec) {
        parser_reinit(p, dset, prn);
        if (p.err) {
            fprintf(stderr, "error in parser_reinit\n".ptr);
            goto gen_finish;
        } else if (p.op == GenSym.inc || p.op == GenSym.dec) {
            /* more or less a no-op: the work is done by
               save_generated_var()
            */
            goto gen_finish;
        } else {
            goto starteval;
        }
    } else {
        int done = 0;

        parser_init(p, s, dset, prn, cast(GenFlags) flags, targtype, &done);
        if (p.err) {
            if (gretl_function_depth() == 0) {
                errmsg(p.err, prn);
            }
            goto gen_finish;
        } else if (done) {
            goto gen_finish;
        }
    }


    if (p.flags & GenFlags.decl) {
        /* check validity of declaration(s) */
        decl_check(p, flags);
        goto gen_finish;
    }

    if (p.op == GenSym.inc || p.op == GenSym.dec) {
        /* implemented via save_generated_var() */
        goto gen_finish;
    }

    /* fire up the lexer */
    lex(p);
    if (p.err) {
        goto gen_finish;
    }

    /* build the syntax tree */
    p.tree = expr(p);
    if (p.err) {
        goto gen_finish;
    }


    if (p.sym != GenSym.eot || p.ch != 0) {
        int c = p.ch;

        if (c == ' ') {
            c = 0;
        } else if (c != 0) {
            parser_ungetc(p);
            c = p.ch;
        }
        context_error(c, p, "realgen");
        goto gen_finish;
    }

    if (flags & GenFlags.noexec) {
        /* we're done at this point */
        goto gen_finish;
    }

    if (!p.err) {
        /* set P_UFRET here if relevant */
        maybe_set_return_flags(p);
    }

 starteval:


    if (autoreg(p)) {
        /* e.g. y = b*y(-1) : evaluate dynamically */
        double * y = p.dset.Z[p.lh.vnum];
        const(double) * x;
        int t, initted = 0;

        for (t=p.dset.t1; t<=p.dset.t2 && !p.err; t++) {
            /* initialize for this observation */
            p.obs = t;
            if (dataset_is_panel(p.dset) && t % p.dset.pd == 0) {
                initted = 0;
            }
            p.ret = eval(p.tree, p);
            if (p.ret != null && p.ret.t == GenSym.series) {
                x = p.ret.v.xvec;
                if (!initted && na(x[t])) {
                    ; /* don't overwrite initializer */
                } else {
                    if (p.op == GenSym.b_asn) {
                        y[t] = x[t];
                    } else {
                        y[t] = xy_calc(y[t], x[t], p.op, GenSym.series, p);
                    }
                    initted = 1;
                }
            } else {
                autoreg_error(p, t);
            }
            if (t == p.dset.t1) {
                p.flags = cast(GenFlags)(p.flags & ~GenFlags.start);
            }
        }
    } else {
        /* standard non-dynamic evaluation */
        p.ret = eval(p.tree, p);
    }

    if (p.flags & GenFlags.exec) {
        p.callcount += 1;
    }

 gen_finish:

    if (p.errprn != null) {
        /* Pick and forward any error message that may not be
           seen if realgen was invoked with a NULL value for
           the printer @prn.
        */
        const(char) * buf = gretl_print_get_buffer(p.errprn);

        if (buf != null && *buf != '\0') {
            gretl_errmsg_set(buf);
        }
        gretl_print_destroy(p.errprn);
        p.errprn = null;
    }


    return p.err;
}

