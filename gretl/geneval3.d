/* Port of the third part (lines 9657-14333) of lib/src/geneval.c: geneval3.c. See geneval0.d for the notes that
   apply to all five parts of the file. */

module gretl.geneval3;


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
import gretl.genmainh : RetrievalIndex, BundleDataIndex;
import core.stdc.string : strncat, strncmp;
import core.stdc.math : acos, asin, atan, cos, sin, tan, acosh, asinh, atanh, cosh, sinh, tanh, sqrt;
import core.stdc.time : tm, mktime;
version (Posix) import core.sys.posix.time : strptime; else import gretl.winport : strptime;
import core.stdc.stdio : sprintf, sscanf;
import gretl.geneval2 : node_get_list, node_get_midas_method;
import gretl.calendar : day_span, epoch_day_from_ymd, fill_monthlen_array, get_days_in_month, gretl_alt_strfdate, gretl_strfdate, gretl_strftime, gretl_strptime, iso_week_number, legacy_day_of_week, ymd_basic_from_epoch_day, ymd_extended_from_epoch_day;
import gretl.dataio : dateton;
import gretl.dataset : build_stacked_series, destroy_dataset, is_string_valued, overwrite_err, series_get_string_for_obs, series_set_string_val, set_dataset_is_changed, string_series_assign_value;
import gretl.dataseth : dataset_has_markers, dataset_is_panel, dataset_is_seasonal, dataset_is_time_series, sample_size;
import gretl.describe : acf_matrix, get_last_result_data, gretl_gini, gretl_kurtosis, gretl_long_run_variance, gretl_max, gretl_mean, gretl_median, gretl_min, gretl_normtest_matrix, gretl_quantile, gretl_skewness, gretl_sst, gretl_stddev, gretl_sum, gretl_variance, ljung_box, periodogram_matrix;
import gretl.geneval1 : aux_empty_series_node, aux_list_node, aux_node_for_type, aux_series_node, aux_sized_matrix_node, c_xy_calc, cmatrix_xy_calc, copy_colvec_as_series, e_types, function_real_only, gen_type_from_gretl_type, gretl_type_from_gen_type, lag_calc, mat_node_get_real_matrix, n_args_error, node_get_bool, node_get_guint32, node_get_int, node_get_matrix, node_get_real_matrix, node_get_scalar, prep_strvals_return, ptr_node_get_bundle, ptr_node_get_matrix, ptr_node_get_referent_node, ptr_node_get_uvar, real_matrix_calc, reset_p_aux, rmatrix_xy_calc, strvals_return_check, typestr, xy_calc;
import gretl.geneval2 : mspec_get_array_index, mspec_get_series_index, real_apply_func, stored_list_ok, sub_addr_get_data;
import gretl.geneval4 : node_type_error;
import gretl.geneval5 : eval, get_opstr, series_to_matrix;
import gretl.gretl_errors : gretl_errmsg_set, gretl_errmsg_sprintf, gretl_errmsg_sprintf_replace, gretl_warnmsg_set;
import gretl.gretl_list : gretl_consecutive_list_new, gretl_list_new, gretl_list_to_vector;
import gretl.gretl_matrix : gretl_dgeev, gretl_matrix_alloc, gretl_matrix_block_resample, gretl_matrix_center, gretl_matrix_column_sd, gretl_matrix_copy, gretl_matrix_fill, gretl_matrix_free, gretl_matrix_get_t1, gretl_matrix_get_t2, gretl_matrix_global_minmax, gretl_matrix_global_sum, gretl_matrix_init, gretl_matrix_lag, gretl_matrix_pca, gretl_matrix_quantiles, gretl_matrix_resample, gretl_matrix_seq, gretl_matrix_shape, gretl_matrix_standardize, gretl_matrix_trim_rows, gretl_matrix_values, gretl_matrix_varsimul, gretl_null_matrix_new, gretl_toeplitz_solve, gretl_vector_sort;
import gretl.gretl_prn : Prn, pprintf, pputc;
import gretl.gretl_string_table : SeriesTable;
import gretl.gretl_utils : accept_as_discrete, gretl_isconst, gretl_isdummy;
import gretl.gstrfuncs : g_strchomp, g_utf8_strlen;
import gretl.libset : get_matrix_mask, libset_get_bool;
import gretl.missing : series_adjust_sample;
import gretl.missingh : NADBL, na;
import gretl.plugins : get_plugin_function;
import gretl.qr_estimate : long_run_covariance;
import gretl.strutils : gretl_strdup, gretl_utf8_replace_char, strings_array_free, strings_array_new;
import gretl.usermath : SelType, is_sel_dummy, mspec_get_element;
import gretl.matrix_extra : gretl_matrix_data_subset, gretl_matrix_data_subset_special, gretl_matrix_from_scalar;
import gretl.gretl_cmatrix : CDFunc, apply_cmatrix_cfunc, apply_cmatrix_dfunc, apply_cmatrix_unary_op, gretl_cmatrix_QR_decomp, gretl_cmatrix_QR_pivot_decomp, gretl_cmatrix_fill, gretl_zgees, gretl_zgeev;
import gretl.gretl_panel : panel_isconst;
import gretl.printout : bufgets, bufgets_finalize, bufgets_peek_line_length, query_bufgets_init;
import gretl.var : vma_rep;
import gretl.graphing : get_user_colors;
import gretl.gretl_array : gretl_array_destroy, gretl_array_get_content_type, gretl_array_get_data, gretl_array_get_element, gretl_array_get_length, gretl_array_get_type, gretl_array_new, gretl_array_set_element, gretl_array_set_string;
import gretl.gretl_bundle : bundle_from_model, bundle_from_system, bundle_get_virtual_object, get_sysinfo_bundle, gretl_bundle_copy, gretl_bundle_delete_data, gretl_bundle_destroy, gretl_bundle_donate_data, gretl_bundle_donate_series, gretl_bundle_extract_args, gretl_bundle_get_element, gretl_bundle_get_keys, gretl_bundle_get_member_type, gretl_bundle_get_scalar, gretl_bundle_get_target, gretl_bundle_has_key, gretl_bundle_new, gretl_bundle_rekey_data, gretl_bundle_set_data, gretl_bundle_set_note, gretl_bundle_set_series;
import gretl.uservar : get_user_var_by_data, user_string_reset, user_string_resize, user_var_replace_value, user_var_set_scalar_value;
import gretl.gretl_typemap : gretl_is_array_type, gretl_is_scalar_type, gretl_type_get_name, gretl_type_get_order;
import gretl.genfuncs : R_from_omega, aggregate_by, block_resample_series, boxcox_matrix, boxcox_series, butterworth_filter, cum_series, diff_series, empirical_cdf, fill_day_of_week_array, fill_isoweek_array, fracdiff_series, gretl_irr, gretl_kmeans, gretl_npv, gretl_sort_by, hp_filter, interpolate_series, matrix_aggregate, maybe_get_values_map, midas_gradient, midas_weights, multi_acf, multi_xcf, omega_from_R, orthdev_series, oshp_filter, panel_expand, panel_shrink, panel_statistic, poly_trend, rank_series, resample_series, sample_span, seasonally_adjust_series, series_get_nobs, series_sum_all, sort_series, tramo_linearize_series;
import gretl.usermat : assign_scalar_to_submatrix, check_matrix_subspec, gretl_matrix_set_part, matrix_get_submatrix, matrix_replace_submatrix, user_gensymm_eigenvals, user_matrix_QR_decomp, user_matrix_SVD;
import gretl.random : halton_matrix, inverse_wishart_matrix, inverse_wishart_sequence;

/* unported dependencies: per-file declarations from the 2026c headers */
extern (C) {
    int fn_n_params(const(Ufunc) * fun);
    GretlType fncall_get_return_type(Fncall * fc);
    const(char) * getsymb(int t);
    int gretl_function_depth();
    int gretl_function_exec_full(Fncall * call, int rtype, DatasetRaw * dset, void * ret, char ** descrip, SeriesTable ** stab, Prn * prn);
    void gretl_pop_c_numeric_locale();
    void gretl_push_c_numeric_locale();
    Node * newb2(int t, Node * l, Node * r);
    Node * newempty();
    int object_is_const(const(char) * name, int vnum);
    int push_function_arg(Fncall * fc, const(char) * name, void * uvar, GretlType type, void * value);
    Fncall * user_func_get_fncall(Ufunc * fun);
    int user_func_get_return_type(const(Ufunc) * fun);
}

Node * aux_sized_array_node(GretlType type, int n, Parser * p)
{
    Node * ret = aux_array_node(p);

    if (ret != null) {
        ret.v.a = gretl_array_new(type, n, &p.err);
    }

    return ret;
}

/* strftime_node() implements both strftime() and strfday().  In both
   cases the @l node must hold a scalar, vector, or numeric series;
   for strftime() this is interpreted as a time_t, while for strfday
   it is interpreted as a GLib "julian day" or "epoch day".

   isodate() also piggy-backs here when a string-valued series is
   wanted on return.
*/

Node * strftime_node(Node * l, Node * r, Node * o, int f,
                     int julian, Parser * p)
{
    Node * ret = null;
    const(char) * fmt = null;
    int * vmap = null;
    char ** S = null;
    double offset = NADBL;
    int nv = 0;

    if (l.t == GenSym.series) {
        if (f != GenFunc.f_isodate && f != GenFunc.f_juldate) {
            /* this check has not been performed yet */
            p.err = strvals_return_check(null, p);
        } else {
            p.err = strvals_return_check(l.parent, p);
        }
        if (!p.err) {
            ret = aux_series_node(p);
        }
    } else if (l.t == GenSym.mat) {
        nv = gretl_vector_get_length(l.v.m);
        if (nv == 0) {
            p.err = GretlError.invarg;
        } else if (nv == 1) {
            ret = aux_string_node(p);
        } else {
            ret = aux_sized_array_node(GretlType.strings, nv, p);
        }
    } else {
        ret = aux_string_node(p);
    }

    if (ret == null) {
        return null;
    }

    /* if @r isn't null or empty it should hold a format string */
    if (r != null) {
        if (r.t == GenSym.str) {
            fmt = r.v.str;
        } else if (r.t != GenSym.empty) {
            p.err = e_types(r);
        }
    }

    /* if @o isn't null or empty it should hold an offset in seconds */
    if (!null_node(o)) {
        offset = node_get_scalar(o, p);
    }

    if (!p.err && l.t == GenSym.series) {
        const(double) * x = l.v.xvec + p.dset.t1;
        int n = sample_size(p.dset);

        vmap = maybe_get_values_map(x, n, &nv, &p.err);
        if (!p.err) {
            S = strings_array_new(nv);
            if (S == null) {
                p.err = GretlError.alloc;
            }
        }
    }

    if (!p.err) {
        char[64] buf = 0;
        char * dstr;
        int t, t1 = 0, t2 = 0;
        double tx;
        long tt;
        guint32 ed;
        int bytes;
        int i = 0;

        if (fmt == null && (f == GenFunc.f_strfday || f == GenFunc.f_isodate)) {
            /* default to ISO 8601 */
            fmt = "%Y-%m-%d";
        }
        if (l.t == GenSym.series) {
            t1 = p.dset.t1;
            t2 = p.dset.t2;
        } else if (l.t == GenSym.mat) {
            t2 = nv - 1;
        }

        for (t=t1; t<=t2 && !p.err; t++) {
            tx = (l.t == GenSym.series)? l.v.xvec[t] : (l.t == GenSym.mat)? l.v.m.val[t] :
                node_get_scalar(l, p);
            if (na(tx)) {
                if (l.t == GenSym.series) {
                    ret.v.xvec[t] = NADBL;
                } else if (l.t != GenSym.mat) {
                    ret.v.str = gretl_strdup("");
                }
            } else if (vmap != null && S[vmap[i]] != null) {
                /* a repeated date in series mode */
                ret.v.xvec[t] = vmap[i] + 1;
                i++;
            } else {
                bytes = 0;
                if (f == GenFunc.f_strfday) {
                    ed = cast(guint32) floor(tx);
                    bytes = gretl_strfdate(buf.ptr, buf.sizeof, fmt, ed);
                } else if (f == GenFunc.f_isodate) {
                    ed = cast(guint32) floor(tx);
                    bytes = gretl_alt_strfdate(buf.ptr, buf.sizeof, julian, ed);
                } else {
                    tt = cast(long) floor(tx);
                    bytes = gretl_strftime(buf.ptr, buf.sizeof, fmt, tt, offset);
                }
                if (bytes > 0) {
                    dstr = gretl_strdup(g_strchomp(buf.ptr));
                    if (l.t == GenSym.series) {
                        if (vmap != null) {
                            ret.v.xvec[t] = vmap[i] + 1;
                            S[vmap[i++]] = dstr;
                        } else {
                            ret.v.xvec[t] = i + 1;
                            S[i++] = dstr;
                        }
                    } else if (l.t == GenSym.mat && nv > 1) {
                        p.err = gretl_array_set_string(ret.v.a, t, dstr, 0);
                    } else {
                        ret.v.str = dstr;
                    }
                } else if (l.t == GenSym.series) {
                    ret.v.xvec[t] = NADBL;
                } else {
                    ret.v.str = gretl_strdup("");
                }
            }
        }
    }

    if (l.t == GenSym.series) {
        if (p.err && S != null) {
            strings_array_free(S, nv);
        } else if (!p.err) {
            prep_strvals_return(S, nv, p);
        }
        free(vmap);
    }

    return ret;
}

Node * isodate_node(Node * l, Node * r, Node * f, Parser * p)
{
    int julian = (f.t == GenFunc.f_juldate);
    int as_string;
    int n = 0;
    Node * ret = null;

    as_string = node_get_bool(r, p, 0);

    if (l.t == GenSym.mat) {
        n = gretl_vector_get_length(l.v.m);
        if (n == 0) {
            p.err = GretlError.nonconf;
        }
    } else if (!scalar_node(l) && l.t != GenSym.series) {
        node_type_error(f.t, 1, GenSym.num, l, p);
    }

    if (p.err) {
        return null;
    } else if ((l.t == GenSym.series || l.t == GenSym.mat) && as_string) {
        if (f.t == GenFunc.f_isodate) {
            p.err = strvals_return_check(f, p);
        }
        if (p.err) {
            return null;
        } else {
            return strftime_node(l, null, null, GenFunc.f_isodate, julian, p);
        }
    }

    if (l.t == GenSym.mat) {
        ret = aux_sized_matrix_node(p, n, 1, 0);
    } else if (l.t == GenSym.num) {
        ret = as_string ? aux_string_node(p) : aux_scalar_node(p);
    } else {
        ret = aux_series_node(p);
    }

    if (ret != null && scalar_node(l)) {
        double x = node_get_scalar(l, p);

        if (!as_string && na(x)) {
            ret.v.xval = NADBL;
        } else if (x < 1 || x > uint.max) {
            p.err = GretlError.invarg;
        } else if (as_string) {
            ret.v.str = ymd_extended_from_epoch_day(cast(guint32) x,
                                                    julian, &p.err);
        } else {
            ret.v.xval = ymd_basic_from_epoch_day(cast(guint32) x,
                                                  julian, &p.err);
        }
    } else if (ret != null) {
        double xt;
        double * targ;
        int t1 = l.t == GenSym.mat ? 0 : p.dset.t1;
        int t2 = l.t == GenSym.mat ? n-1 : p.dset.t2;
        int t;

        for (t=t1; t<=t2; t++) {
            targ = l.t == GenSym.mat ? &ret.v.m.val[t] : &ret.v.xvec[t];
            xt = l.t == GenSym.mat ? l.v.m.val[t] : l.v.xvec[t];
            if (na(xt)) {
                *targ = NADBL;
            } else if (xt >= 1 && xt <= uint.max) {
                *targ = ymd_basic_from_epoch_day(cast(guint32) xt,
                                                 julian, &p.err);
            } else {
                p.err = GretlError.invarg;
                break;
            }
        }
    }

    return ret;
}

int strptime_error(const(char) * s, double x, int t, int ntype)
{
    if (ntype == GenSym.str) {
        gretl_errmsg_sprintf("'%s': argument doesn't match format".ptr, s);
    } else if (ntype == GenSym.num) {
        gretl_errmsg_sprintf("Invalid date number %.12g".ptr, x);
    } else if (s != null) {
        gretl_errmsg_sprintf("'%s': argument doesn't match format at obs %d\n".ptr, s, t);
    } else if (!na(x)) {
        gretl_errmsg_sprintf("Invalid date number %.12g at obs %d\n".ptr,
                             x, t);
    }

    return GretlError.invarg;
}

int am_pm_error(const(char) * s)
{
    char[3] tmp = 0;

    strncat(tmp.ptr, s, 2);
    gretl_errmsg_sprintf(_("Unrecognized content '%s' in time string").ptr, tmp.ptr);
    return GretlError.invarg;
}

/* strptime_node() implements both strptime() and strpday(). In the
   primary use case the @l node contains a string, array of strings,
   or string-valued series. In the secondary case @l contains ISO 8601
   "basic" dates (scalar, series or vector).
*/

Node * strptime_node(Node * l, Node * r, int f, Parser * p)
{
    Node * ret = null;
    const(char) * fmt = null;
    const(char) * src = null;
    int str_series = 0;
    int v = -1;
    int n = 0;

    if (l.t == GenSym.series) {
        /* either string-valued or ISO basic */
        v = l.vnum;
        if (strvals_node(l, p)) {
            str_series = 1;
        }
        ret = aux_series_node(p);
    } else if (l.t == GenSym.array) {
        /* must be strings array */
        if (gretl_array_get_type(l.v.a) != GretlType.strings) {
            p.err = GretlError.invarg;
        } else {
            n = gretl_array_get_length(l.v.a);
            ret = aux_sized_matrix_node(p, n, 1, 0);
        }
    } else if (l.t == GenSym.mat) {
        /* must be ISO basic */
        n = gretl_vector_get_length(l.v.m);
        if (n == 0) {
            p.err = GretlError.invarg;
        } else {
            ret = aux_sized_matrix_node(p, n, 1, 0);
        }
    } else {
        /* single string or ISO basic */
        if (l.t == GenSym.str) {
            src = l.v.str;
        }
        ret = aux_scalar_node(p);
    }

    if (ret == null) {
        return null;
    }

    if (src == null && !str_series && l.t != GenSym.array) {
        /* without a string arg we won't accept a format string */
        if (r.t != GenSym.empty) {
            p.err = GretlError.invarg;
        }
    } else {
        /* if @r isn't empty it should hold a format string */
        if (r.t == GenSym.str) {
            fmt = r.v.str;
        } else if (r.t != GenSym.empty) {
            p.err = e_types(r);
        }
    }

    if (!p.err) {
        double * targ = null;
        int t, t1 = 0, t2 = 0;
        char[32] buf;
        const(char) * s;

        if (l.t == GenSym.series) {
            t1 = p.dset.t1;
            t2 = p.dset.t2;
        } else if (l.t == GenSym.array) {
            t2 = gretl_array_get_length(l.v.a) - 1;
        } else if (l.t == GenSym.mat) {
            t2 = n - 1;
        } else {
            targ = &ret.v.xval;
        }

        for (t=t1; t<=t2 && !p.err; t++) {
            tm tmv;
            int handled = 0;
            double x = NADBL;

            tmv.tm_mday = 1;
            tmv.tm_isdst = -1;

            if (l.t == GenSym.series) {
                if (str_series) {
                    src = series_get_string_for_obs(p.dset, v, t);
                } else {
                    x = l.v.xvec[t];
                }
                targ = &ret.v.xvec[t];
            } else if (l.t == GenSym.array) {
                if ((src = cast(const(char) *) gretl_array_get_data(l.v.a, t)) == null) {
                    p.err = GretlError.invarg;
                    break;
                }
                targ = &ret.v.m.val[t];
            } else if (l.t == GenSym.mat) {
                x = l.v.m.val[t];
                targ = &ret.v.m.val[t];
            } else if (l.t == GenSym.num) {
                x = l.v.xval;
            }

            if (src == null) {
                /* has to be ISO 8601 basic */
                if (na(x)) {
                    *targ = NADBL;
                    continue;
                } else if (x < 0 || x != floor(x) || x > 99991231) {
                    p.err = strptime_error(src, x, t, l.t);
                    break;
                } else {
                    sprintf(buf.ptr, "%08d", cast(int) x);
                    s = strptime(buf.ptr, "%Y%m%d", &tmv);
                }
            } else {
                if (fmt == null) {
                    /* default to ISO 8601 extended */
                    fmt = "%Y-%m-%d";
                }
                if (f == GenFunc.f_strpday) {
                    s = strptime(src, fmt, &tmv);
                } else {
                    s = gretl_strptime(src, fmt, targ);
                    handled = 1;
                }
            }

            if (s == null) {
                /* strptime() failed */
                p.err = strptime_error(src, x, t, l.t);
                break;
            } else if (!strncmp(s, "AM", 2) || !strncmp(s, "PM", 2)) {
                p.err = am_pm_error(s);
            } else if (f == GenFunc.f_strpday) {
                int y = tmv.tm_year + 1900;
                int m = tmv.tm_mon + 1;

                *targ = cast(double) epoch_day_from_ymd(y, m, tmv.tm_mday);
                if (*targ == 0) {
                    *targ = NADBL;
                }
            } else if (!handled) {
                /* F_STRPTIME */
                *targ = cast(double) mktime(&tmv);
            }
        }
    }

    return ret;
}

Node * atof_node(Node * l, Parser * p)
{
    Node * ret = null;
    const(char) * endptr = null;
    int warn = 0;

    errno = 0;

    gretl_push_c_numeric_locale();

    if (l.t == GenSym.str) {
        ret = aux_scalar_node(p);
        if (ret != null && starting(p)) {
            if (strchr(l.v.str, ',') != null) {
                warn = 1;
            }
            ret.v.xval = strtod(l.v.str, &endptr);
            if (errno || endptr == l.v.str) {
                errno = 0;
                ret.v.xval = NADBL;
            }
        }
    } else if (l.t == GenSym.series) {
        int v = l.vnum;

        if (!strvals_node(l, p)) {
            p.err = GretlError.types;
        } else {
            ret = aux_series_node(p);
        }
        if (ret != null && starting(p)) {
            const(char) * st;
            int t;

            for (t=p.dset.t1; t<=p.dset.t2; t++) {
                st = series_get_string_for_obs(p.dset, v, t);
                if (st == null) {
                    /* happens if obs @t is missing */
                    ret.v.xvec[t] = NADBL;
                } else {
                    if (!warn && strchr(st, ',') != null) {
                        warn = 1;
                    }
                    ret.v.xvec[t] = strtod(st, &endptr);
                    if (errno || endptr == st) {
                        errno = 0;
                        ret.v.xvec[t] = NADBL;
                    }
                }
            }
        }
    }

    if (warn) {
        gretl_warnmsg_set(_("comma(s) found in atof() input").ptr);
    }

    gretl_pop_c_numeric_locale();

    return ret;
}

void strip_newline(char * s)
{
    if (s != null && *s != '\0') {
        int i, len = cast(int) strlen(s);

        for (i=len-1; i>=0; i--) {
            if (s[i] == '\n' || s[i] == '\r') {
                s[i] = '\0';
            } else {
                break;
            }
        }
    }
}

Node * getline_node(Node * l, Node * r, Parser * p)
{
    Node * ret = aux_scalar_node(p);

    if (ret != null && starting(p)) {
        const(char) * buf = null;
        Node * rs = null;

        if (l.vname == null) {
            gretl_errmsg_set("getline: the source must be a named string variable");
            p.err = GretlError.invarg;
        } else {
            buf = l.v.str;
            if (null_node(r)) {
                /* clean-up only */
                bufgets_finalize(buf);
                ret.v.xval = 0;
                return ret;
            }
            if (r.t == GenSym.str && r.vname != null) {
                rs = r;
            } else if (r.t == GenSym.u_addr && r.L.t == GenSym.str) {
                rs = ptr_node_get_referent_node(r, p);
            } else {
                gretl_errmsg_set("getline: the target must be a named string variable");
                p.err = GretlError.invarg;
            }
        }

        if (!p.err) {
            p.err = query_bufgets_init(buf);
        }
        if (!p.err) {
            size_t len = bufgets_peek_line_length(buf);

            if (len == 0) {
                bufgets_finalize(buf);
                rs.v.str = user_string_reset(rs.vname, null, &p.err);
                ret.v.xval = 0;
            } else {
                rs.v.str = user_string_resize(rs.vname, len, &p.err);
                if (!p.err) {
                    bufgets(rs.v.str, len, buf);
                    strip_newline(rs.v.str);
                    ret.v.xval = 1;
                }
            }
        }
    }

    return ret;
}

int series_get_start(int t1, int t2, const(double) * x)
{
    int t;

    for (t=t1; t<=t2; t++) {
        if (!na(x[t])) {
            break;
        }
    }

    return t + 1;
}

int series_get_end(int t1, int t2, const(double) * x)
{
    int t;

    for (t=t2; t>=t1; t--) {
        if (!na(x[t])) {
            break;
        }
    }

    return t + 1;
}

void cast_to_series(Node * n, int f, GretlMatrixRaw ** tmp,
                    int * t1, int * t2, Parser * p)
{
    GretlMatrixRaw * m = n.v.m;
    int len = gretl_vector_get_length(m);

    if (gretl_is_null_matrix(m)) {
        p.err = GretlError.data;
    } else if (m.is_complex) {
        node_type_error(f, 1, GenSym.series, n, p);
    } else if (len > 0 && len == p.dset.n) {
        *tmp = m;
        n.v.xvec = m.val;
    } else if (len > 0 && t1 != null && t2 != null) {
        *tmp = m;
        n.v.xvec = m.val;
        *t1 = 0;
        *t2 = len - 1;
    } else {
        node_type_error(f, 1, GenSym.series, n, p);
    }
}

/* Functions taking a series (or matrix) as argument and returning a
   scalar; allowance is made for an additional boolean arg in some
   cases.
*/

Node * series_scalar_func(Node * n, int f,
                          Node * r, Parser * p)
{
    Node * ret = aux_scalar_node(p);

    if (ret != null && starting(p)) {
        GretlMatrixRaw * tmp = null;
        int t1 = p.dset.t1;
        int t2 = p.dset.t2;
        const(double) * x;

        if (n.t == GenSym.mat && gretl_is_complex(n.v.m)) {
            p.err = GretlError.cmplx;
        } else if (n.t == GenSym.mat) {
            /* we got a matrix rather than a series */
            if (f == GenFunc.f_sum || f == GenFunc.f_max || f == GenFunc.f_min) {
                /* we'll sum, max, or min all elements of a matrix */
                if (f == GenFunc.f_sum) {
                    ret.v.xval = gretl_matrix_global_sum(n.v.m, &p.err);
                } else {
                    ret.v.xval = gretl_matrix_global_minmax(n.v.m,
                                                            f == GenFunc.f_max,
                                                            &p.err);
                }
                return ret; /* handled */
            } else {
                cast_to_series(n, f, &tmp, &t1, &t2, p);
            }
        }

        if (f == GenFunc.f_t1 || f == GenFunc.f_t2) {
            /* matrix argument not supported */
            int insample = node_get_bool(r, p, 0);

            if (!p.err && !insample) {
                t1 = 0;
                t2 = p.dset.n - 1;
            }
        }

        if (p.err) {
            return ret;
        }

        /* pointer to series or matrix data */
        x = n.v.xvec;

        switch (f) {
        case GenFunc.f_sum:
            ret.v.xval = gretl_sum(t1, t2, x);
            break;
        case GenFunc.f_sumall:
            ret.v.xval = series_sum_all(t1, t2, x);
            break;
        case GenFunc.f_mean:
            ret.v.xval = gretl_mean(t1, t2, x);
            break;
        case GenFunc.f_sd:
            ret.v.xval = gretl_stddev(t1, t2, x);
            break;
        case GenFunc.f_vce:
            ret.v.xval = gretl_variance(t1, t2, x);
            break;
        case GenFunc.f_sst:
            ret.v.xval = gretl_sst(t1, t2, x);
            break;
        case GenFunc.f_skewness:
            ret.v.xval = gretl_skewness(t1, t2, x);
            break;
        case GenFunc.f_kurtosis:
            ret.v.xval = gretl_kurtosis(t1, t2, x);
            break;
        case GenFunc.f_min:
            ret.v.xval = gretl_min(t1, t2, x);
            break;
        case GenFunc.f_max:
            ret.v.xval = gretl_max(t1, t2, x);
            break;
        case GenFunc.f_median:
            ret.v.xval = gretl_median(t1, t2, x);
            break;
        case GenFunc.f_gini:
            ret.v.xval = gretl_gini(t1, t2, x);
            break;
        case GenFunc.f_nobs:
            ret.v.xval = series_get_nobs(t1, t2, x);
            break;
        case GenFunc.f_isconst:
            ret.v.xval = gretl_isconst(t1, t2, x);
            break;
        case GenFunc.f_isdummy:
            ret.v.xval = gretl_isdummy(t1, t2, x);
            break;
        case GenFunc.f_t1:
            ret.v.xval = series_get_start(t1, t2, x);
            break;
        case GenFunc.f_t2:
            ret.v.xval = series_get_end(t1, t2, x);
            break;
        default:
            break;
        }

        if (n.t == GenSym.mat) {
            n.v.m = tmp;
        }
    }

    return ret;
}

/* Functions normally taking a series or vector as argument and
   returning a scalar, but are evaluated on a scalar, so output is
   trivial.
*/

Node * pretend_matrix_scalar_func(Node * n, int f, Parser * p)
{
    Node * ret = aux_scalar_node(p);

    if (ret != null && starting(p)) {
        if (f == GenFunc.f_sum || f == GenFunc.f_sumall || f == GenFunc.f_mean || f == GenFunc.f_max ||
            f == GenFunc.f_min || f == GenFunc.f_median) {
            ret.v.xval = n.v.xval;
        } else if (f == GenFunc.f_nobs) {
            ret.v.xval = 1;
        } else if (f == GenFunc.f_sd || f == GenFunc.f_vce || f == GenFunc.f_sst || f == GenFunc.f_gini) {
            ret.v.xval = 0;
        } else if (f == GenFunc.f_skewness || f == GenFunc.f_kurtosis) {
            /* this is probably less intuitive than all the above:
               in "normal" cases we're returning NADBL when the variance
               is 0, so we're just being consistent here
            */
            ret.v.xval = NADBL;
        } else {
            /* any other cases not legit */
            node_type_error(f, 0, GenSym.series, n, p);
        }
    }

    return ret;
}

Gretlopt get_normtest_option(Node * n, Parser * p)
{
    Gretlopt opt = Gretlopt.OPT_NONE;

    if (null_node(n)) {
        ; /* OK */
    } else if (n.t == GenSym.str) {
        const(char) * s = n.v.str;

        if (!strcmp(s, "swilk")) {
            opt = Gretlopt.OPT_W;
        } else if (!strcmp(s, "lillie")) {
            opt = Gretlopt.OPT_L;
        } else if (!strcmp(s, "jbera")) {
            opt = Gretlopt.OPT_J;
        } else if (!strcmp(s, "dhansen")) {
            opt = Gretlopt.OPT_D;
        } else if (!strcmp(s, "all")) {
            opt = Gretlopt.OPT_A;
        } else {
            p.err = GretlError.invarg;
        }
    } else {
        p.err = e_types(n);
    }

    return opt;
}

Node * series_matrix_node(Node * l, Node * r, int f, Parser * p)
{
    Node * ret = null;

    if (starting(p)) {
        const(double) * x = null;
        int n = 0, t1 = 0, t2 = 0;
        Gretlopt opt = Gretlopt.OPT_NONE;

        if (f == GenFunc.f_normtest) {
            opt = get_normtest_option(r, p);
            if (p.err) {
                return null;
            }
        }

        if (l.t == GenSym.mat) {
            n = gretl_vector_get_length(l.v.m);
            if (n == 0) {
                p.err = GretlError.types;
            } else {
                x = l.v.m.val;
                t2 = n - 1;
            }
        } else if (f == GenFunc.f_normtest) {
            x = l.v.xvec;
            t1 = p.dset.t1;
            t2 = p.dset.t2;
        } else {
            x = l.v.xvec + p.dset.t1;
            n = sample_size(p.dset);
        }

        if (!p.err) {
            ret = aux_matrix_node(p);
        }

        if (!p.err) {
            if (f == GenFunc.f_normtest) {
                ret.v.m = gretl_normtest_matrix(x, t1, t2, opt, &p.err);
            } else {
                /* F_ECDF */
                ret.v.m = empirical_cdf(x, n, &p.err);
            }
        }
    }

    return ret;
}

/* There must be a matrix in @l; @r may hold a vector or
   a scalar value */

Node * matrix_quantiles_node(Node * l, Node * r, Parser * p)
{
    Node * ret = null;

    if (starting(p)) {
        GretlMatrixRaw * pmat = node_get_real_matrix(r, p, 0, 0);

        if (!p.err) {
            ret = aux_matrix_node(p);
        }
        if (ret != null) {
            ret.v.m = gretl_matrix_quantiles(l.v.m, pmat, &p.err);
        }
    }

    return ret;
}

/* functions taking a series and a scalar as arguments and returning
   a scalar
*/

Node * series_scalar_scalar_func(Node * l, Node * r,
                                 Node * r2, int f,
                                 Parser * p)
{
    Node * ret = null;

    if (f == GenFunc.f_lrvar && l.t == GenSym.series && dataset_is_panel(p.dset)) {
        gretl_errmsg_set("lrvar(): panel data not supported");
        p.err = GretlError.data;
        return ret;
    }

    if (starting(p)) {
        double rval = -1;
        double r2val = NADBL;
        const(double) * xvec;
        int t1 = p.dset.t1;
        int t2 = p.dset.t2;
        int pd = 1;

        if (f == GenFunc.f_lrvar && null_node(r)) {
            ; /* OK, second arg is optional */
        } else {
            rval = node_get_scalar(r, p);
        }

        if (l.t == GenSym.num) {
            if (f == GenFunc.f_npv) {
                /* a scalar first arg may be OK for Net Present Value */
                t1 = 0;
                t2 = 0;
                xvec = &l.v.xval;
            } else {
                /* but not otherwise */
                p.err = GretlError.types;
                return null;
            }
        } else if (l.t == GenSym.mat) {
            int n = gretl_vector_get_length(l.v.m);

            if (n == 0) {
                p.err = GretlError.types;
                return null;
            }
            t1 = 0;
            t2 = n - 1;
            xvec = l.v.m.val;
        } else {
            /* got a series on the left */
            pd = p.dset.pd;
            xvec = l.v.xvec;
        }

        if (f == GenFunc.f_lrvar && !null_node(r2)) {
            /* optional third arg */
            r2val = node_get_scalar(r2, p);
        }

        if (!p.err) {
            ret = aux_scalar_node(p);
        }
        if (p.err) {
            return ret;
        }

        switch (f) {
        case GenFunc.f_lrvar:
            ret.v.xval = gretl_long_run_variance(t1, t2, xvec, cast(int) rval, r2val);
            break;
        case GenFunc.f_quantile:
            ret.v.xval = gretl_quantile(t1, t2, xvec, rval, Gretlopt.OPT_NONE, &p.err);
            break;
        case GenFunc.f_npv:
            ret.v.xval = gretl_npv(t1, t2, xvec, rval, pd, &p.err);
            break;
        case GenFunc.f_isconst:
            ret.v.xval = panel_isconst(t1, t2, pd, xvec, cast(int) rval);
            break;
        default:
            break;
        }

    } else {
        ret = aux_any_node(p);
    }

    return ret;
}

Node * isconst_or_dum_node(Node * l, Node * r, Parser * p, int f)
{
    if (f == GenFunc.f_isdummy || null_node(r)) {
        return series_scalar_func(l, f, null, p);
    } else if (l.t == GenSym.mat) {
        node_type_error(f, 1, GenSym.series, l, p);
        return null;
    } else if (!dataset_is_panel(p.dset)) {
        p.err = GretlError.pdwrong;
        return null;
    } else {
        return series_scalar_scalar_func(l, r, null, f, p);
    }
}

/* Series on left, scalar or string on right, as in
   x[23] or somevar["CA"]. We return the selected
   scalar value from the series -- unless the series
   is string-valued, in which case we return the
   string value for the given observation.
*/

Node * series_obs(Node * l, Node * r, Parser * p)
{
    int strval = strvals_node(l, p);
    Node * ret;

    ret = strval ? aux_string_node(p) : aux_scalar_node(p);

    if (ret != null) {
        int t = -1; /* invalid */

        if (r.t == GenSym.str) {
            t = dateton(r.v.str, p.dset);
            if (t < 0) {
                if (dataset_has_markers(p.dset)) {
                    gretl_errmsg_sprintf_replace(_("Invalid observation specifier \"%s\"").ptr,
                                                 r.v.str);
                }
                p.err = GretlError.data;
            }
        } else {
            /* plain integer */
            t = node_get_int(r, p);
            if (!p.err) {
                if (t > 0 && t <= p.dset.n) {
                    t--; /* OK, convert to zero based */
                } else {
                    gretl_errmsg_sprintf(_("Index value %d is out of bounds").ptr, t);
                    p.err = GretlError.data;
                }
            }
        }

        if (!p.err) {
            if (strval) {
                const(char) * s =
                    series_get_string_for_obs(p.dset, l.vnum, t);

                if (s == null) {
                    ret.v.str = gretl_strdup("");
                } else {
                    ret.v.str = gretl_strdup(s);
                }
            } else {
                ret.v.xval = l.v.xvec[t];
            }
        }
    }

    return ret;
}

Node * series_ljung_box(Node * l, Node * r, Parser * p)
{
    Node * ret = aux_scalar_node(p);

    if (ret != null && starting(p)) {
        const(double) * x = l.v.xvec;
        int k = node_get_int(r, p);
        int t1 = p.dset.t1;
        int t2 = p.dset.t2;

        if (!p.err && k <= 0) {
            gretl_errmsg_sprintf(_("Invalid lag order %d").ptr, k);
            p.err = GretlError.data;
        }

        if (!p.err) {
            p.err = series_adjust_sample(x, &t1, &t2);
        }

        if (!p.err) {
            ret.v.xval = ljung_box(k, t1, t2, x, &p.err);
        }
    }

    return ret;
}

Node * series_polyfit(Node * l, Node * r, Parser * p)
{
    Node * ret = aux_series_node(p);

    if (ret != null && starting(p)) {
        const(double) * x = l.v.xvec;
        int order = node_get_int(r, p);

        if (!p.err) {
            p.err = poly_trend(x, ret.v.xvec, p.dset, order);
        }
    }

    return ret;
}

Node * series_lag(Node * l, Node * r, Parser * p)
{
    Node * ret = null;
    const(double) * x = l.v.xvec;
    int k = -(node_get_int(r, p));
    int t1, t2;

    if (!p.err && l.vnum == 0) {
        gretl_errmsg_set(_("The constant cannot be lagged").ptr);
        p.err = GretlError.types;
    }

    if (!p.err) {
        ret = aux_series_node(p);
    }

    if (ret == null) {
        return null;
    }

    t1 = autoreg(p) ? p.obs : p.dset.t1;
    t2 = autoreg(p) ? p.obs : p.dset.t2;

    lag_calc(ret.v.xvec, x, k, t1, t2, GenSym.b_asn, 1.0, p);

    return ret;
}

Node * series_sort_by(Node * l, Node * r, Parser * p)
{
    Node * ret = aux_series_node(p);

    if (ret != null && starting(p)) {
        if (l.t == GenSym.series && r.t == GenSym.series) {
            p.err = gretl_sort_by(l.v.xvec, r.v.xvec, ret.v.xvec, p.dset);
        } else if (l.t != GenSym.series) {
            p.err = e_types(l);
        } else {
            p.err = e_types(r);
        }
    }

    return ret;
}

Node * vector_sort(Node * l, int f, Parser * p)
{
    Node * ret = (l.t == GenSym.series)? aux_series_node(p) :
        aux_matrix_node(p);

    if (ret != null && starting(p)) {
        if (l.t == GenSym.series) {
            p.err = sort_series(l.v.xvec, ret.v.xvec, f, p.dset);
        } else if (l.t == GenSym.num) {
            ret.v.m = gretl_matrix_from_scalar(l.v.xval);
        } else if (gretl_is_null_matrix(l.v.m)) {
            ret.v.m = gretl_null_matrix_new();
        } else if (l.v.m.is_complex) {
            p.err = GretlError.cmplx;
        } else {
            int descending = (f == GenFunc.f_dsort);

            ret.v.m = gretl_vector_sort(l.v.m, descending, &p.err);
        }
    }

    return ret;
}

Node * vector_values(Node * l, int f, Parser * p)
{
    Node * ret = aux_matrix_node(p);

    if (ret != null && starting(p)) {
        const(double) * x = null;
        int n = 0;

        if (l.t == GenSym.num) {
            n = 1;
            x = &l.v.xval;
        } else if (l.t == GenSym.series) {
            n = sample_size(p.dset);
            x = l.v.xvec + p.dset.t1;
        } else if (gretl_is_null_matrix(l.v.m)) {
            ret.v.m = gretl_null_matrix_new();
        } else {
            n = gretl_vector_get_length(l.v.m);
            x = l.v.m.val;
        }

        if (n > 0 && x != null) {
            Gretlopt opt = (f == GenFunc.f_values)? Gretlopt.OPT_S : Gretlopt.OPT_NONE;

            ret.v.m = gretl_matrix_values(x, n, opt, &p.err);
        } else if (ret.v.m == null) {
            p.err = GretlError.data;
        }
    }

    return ret;
}

Node * do_irr(Node * l, Parser * p)
{
    Node * ret = aux_scalar_node(p);

    if (ret != null && starting(p)) {
        const(double) * x = null;
        int pd = 1, n = 0;

        if (l.t == GenSym.num) {
            n = 1;
            x = &l.v.xval;
        } else if (l.t == GenSym.series) {
            n = sample_size(p.dset);
            x = l.v.xvec + p.dset.t1;
            pd = p.dset.pd;
        } else if (!gretl_is_null_matrix(l.v.m)) {
            n = gretl_vector_get_length(l.v.m);
            x = l.v.m.val;
        }

        if (n > 0 && x != null) {
            ret.v.xval = gretl_irr(x, n, pd, &p.err);
        } else {
            p.err = GretlError.data;
        }
    }

    return ret;
}

int deseasonalize(const(double) * x, double * y,
                  const(char) * vname, int f,
                  Node * r, Parser * p)
{
    GretlBundle * b = null;
    int tramo = 0;
    int err = 0;

    if (!null_node(r)) {
        if (r.t == GenSym.str) {
            /* legacy usage */
            if (!strcmp(r.v.str, "T")) {
                tramo = 1;
            } else if (strcmp(r.v.str, "X")) {
                err = GretlError.invarg;
            }
        } else {
            /* current style */
            b = node_get_bundle(r, p);
        }
    }
    if (!err) {
        err = seasonally_adjust_series(x, y, vname, p.dset,
                                       tramo, b, p.prn);
    }

    return err;
}

bool is_panel_stat(int f)
{
    return (f == GenFunc.f_pnobs ||
            f == GenFunc.f_pmin ||
            f == GenFunc.f_pmax ||
            f == GenFunc.f_psum ||
            f == GenFunc.f_pmean ||
            f == GenFunc.f_pxsum ||
            f == GenFunc.f_pxnobs ||
            f == GenFunc.f_pxmean ||
            f == GenFunc.f_psd);
}

/* Functions taking a series as argument and returning a series.  Note
   that the @r node may be null or may contain an auxiliary parameter,
   as follows:

   @f == F_DESEAS: bundle or (backward compat) string
   @f == panel stats function: series
   otherwise: scalar

   In case @f == F_HPFILT, the @o node can be used to pass
   an additional optional argument.
*/

Node * series_series_func(Node * l, Node * r, Node * o,
                          int f, Parser * p)
{
    Node * ret = null;

    if ((f == GenFunc.f_sdiff || f == GenFunc.f_deseas) &&
        !dataset_is_seasonal(p.dset)) {
        p.err = GretlError.pdwrong;
        return null;
    } else if ((f == GenFunc.f_fracdiff || f == GenFunc.f_fraclag ||
                f == GenFunc.f_hpfilt || f == GenFunc.f_tramolin) &&
               !dataset_is_time_series(p.dset)) {
        gretl_errmsg_set(_("This function requires time-series data").ptr);
        p.err = GretlError.data;
        return null;
    }

    /* check the @r argument */
    if ((f == GenFunc.f_fracdiff || f == GenFunc.f_fraclag) && null_node(r)) {
        /* fractional difference needed */
        p.err = GretlError.args;
    } else if (null_node(r)) {
        ; /* not present, OK */
    } else if (f == GenFunc.f_deseas) {
        if (r.t != GenSym.str && r.t != GenSym.bundle && r.t != GenSym.u_addr) {
            p.err = e_types(r);
        }
    } else if (is_panel_stat(f)) {
        if (r.t != GenSym.series) {
            node_type_error(f, 2, GenSym.series, r, p);
        }
    } else if (!scalar_node(r)) {
        node_type_error(f, 2, GenSym.num, r, p);
    }

    if (l.t == GenSym.mat && f == GenFunc.f_boxcox) {
        /* Do all columns of matrix input: this could
           be generalized to some other functions?
        */
        double d = node_get_scalar(r, p);

        if (!p.err) {
            ret = aux_matrix_node(p);
        }
        if (!p.err) {
            ret.v.m = boxcox_matrix(l.v.m, d, &p.err);
        }
        return ret;
    }

    ret = aux_series_node(p);

    if (ret != null) {
        GretlMatrixRaw * tmp = null;
        const(char) * vname = null;
        double parm = NADBL;
        int oneside = 0;
        const(double) * z = null;
        const(double) * x;
        double * y;

        if (l.t == GenSym.mat) {
            cast_to_series(l, f, &tmp, null, null, p);
        } else if (l.vnum > 0) {
            vname = p.dset.varname[l.vnum];
        }
        if (!null_node(r)) {
            if (r.t == GenSym.series) {
                z = r.v.xvec;
            } else if (scalar_node(r)) {
                parm = node_get_scalar(r, p);
            }
        }
        if (p.err) {
            return null;
        }

        x = l.v.xvec;
        y = ret.v.xvec;

        switch (f) {
        case GenFunc.f_hpfilt:
            oneside = node_get_bool(o, p, 0);
            if (!p.err && oneside) {
                p.err = oshp_filter(x, y, p.dset, parm, Gretlopt.OPT_NONE);
            } else if (!p.err) {
                p.err = hp_filter(x, y, p.dset, parm, Gretlopt.OPT_NONE);
            }
            break;
        case GenFunc.f_fracdiff:
            p.err = fracdiff_series(x, y, parm, 1, autoreg(p) ? p.obs : -1, p.dset);
            break;
        case GenFunc.f_fraclag:
            p.err = fracdiff_series(x, y, parm, 0, autoreg(p) ? p.obs : -1, p.dset);
            break;
        case GenFunc.f_boxcox:
            p.err = boxcox_series(x, y, parm, p.dset);
            break;
        case GenFunc.f_diff:
        case GenFunc.f_ldiff:
        case GenFunc.f_sdiff:
            p.err = diff_series(x, y, f, p.dset);
            break;
        case GenFunc.f_odev:
            p.err = orthdev_series(x, y, p.dset);
            break;
        case GenFunc.f_cum:
            p.err = cum_series(x, y, p.dset);
            break;
        case GenFunc.f_deseas:
            p.err = deseasonalize(x, y, vname, f, r, p);
            break;
        case GenFunc.f_tramolin:
            p.err = tramo_linearize_series(x, y, p.dset);
            break;
        case GenFunc.f_interpol:
            p.err = interpolate_series(x, y, p.dset);
            break;
        case GenFunc.f_resample:
            if (!na(parm)) {
                p.err = block_resample_series(x, y, cast(int) parm, p.dset);
            } else {
                p.err = resample_series(x, y, p.dset);
            }
            break;
        case GenFunc.f_pnobs:
        case GenFunc.f_pmin:
        case GenFunc.f_pmax:
        case GenFunc.f_psum:
        case GenFunc.f_pmean:
        case GenFunc.f_pxsum:
        case GenFunc.f_pxnobs:
        case GenFunc.f_pxmean:
        case GenFunc.f_psd:
            p.err = panel_statistic(x, y, p.dset, f, z);
            break;
        case GenFunc.f_ranking:
            p.err = rank_series(x, y, GenFunc.f_sort, p.dset);
            break;
        default:
            break;
        }

        if (l.t == GenSym.mat) {
            l.v.m = tmp;
        }
    }

    return ret;
}

Node * do_panel_shrink(Node * l, int noskip, Parser * p)
{
    Node * ret = aux_matrix_node(p);

    if (ret != null && starting(p)) {
        ret.v.m = panel_shrink(l.v.xvec, noskip, p.dset, &p.err);
    }

    return ret;
}

Node * do_panel_expand(Node * l, Node * r, Parser * p)
{
    Node * ret = aux_series_node(p);

    if (ret != null && starting(p)) {
        Gretlopt opt = Gretlopt.OPT_NONE;

        if (node_get_bool(r, p, 0)) {
            opt = Gretlopt.OPT_X; /* expand by unit, not time */
        }
        if (!p.err) {
            p.err = panel_expand(l.v.m, ret.v.xvec, opt, p.dset);
        }
    }

    return ret;
}

/* pergm function takes series or column vector arg, returns matrix:
   if we come up with more functions on that pattern, the following
   could be extended
*/

Node * pergm_node(Node * l, Node * r, Parser * p)
{
    Node * ret = null;

    if (!null_or_scalar(r)) {
        /* optional 'r' node must be scalar */
        node_type_error(GenFunc.f_pergm, 2, GenSym.num, r, p);
    } else if (l.t == GenSym.mat && gretl_vector_get_length(l.v.m) == 0) {
        /* if 'l' node is not a series, must be a vector */
        node_type_error(GenFunc.f_pergm, 1, GenSym.series, l, p);
    } else {
        ret = aux_matrix_node(p);
    }

    if (!p.err) {
        const(double) * x = null;
        int t1 = 0, t2 = 0;
        int width = -1;

        if (l.t == GenSym.series) {
            x = l.v.xvec;
            t1 = p.dset.t1;
            t2 = p.dset.t2;
        } else if (l.t == GenSym.mat) {
            x = l.v.m.val;
            t1 = 0;
            t2 = gretl_vector_get_length(l.v.m) - 1;
        }

        if (!null_node(r)) {
            width = node_get_int(r, p);
        }

        if (!p.err) {
            ret.v.m = periodogram_matrix(x, t1, t2, width, &p.err);
        }
    }

    return ret;
}

/* the C99 complex functions, declared here with our own complex type
   (Cmplx is ABI-compatible with double complex) */

extern (C) {
    pragma(mangle, "cacos")  Cmplx c_cacos(Cmplx z);
    pragma(mangle, "casin")  Cmplx c_casin(Cmplx z);
    pragma(mangle, "catan")  Cmplx c_catan(Cmplx z);
    pragma(mangle, "ccos")   Cmplx c_ccos(Cmplx z);
    pragma(mangle, "csin")   Cmplx c_csin(Cmplx z);
    pragma(mangle, "ctan")   Cmplx c_ctan(Cmplx z);
    pragma(mangle, "cacosh") Cmplx c_cacosh(Cmplx z);
    pragma(mangle, "casinh") Cmplx c_casinh(Cmplx z);
    pragma(mangle, "catanh") Cmplx c_catanh(Cmplx z);
    pragma(mangle, "ccosh")  Cmplx c_ccosh(Cmplx z);
    pragma(mangle, "csinh")  Cmplx c_csinh(Cmplx z);
    pragma(mangle, "ctanh")  Cmplx c_ctanh(Cmplx z);
    pragma(mangle, "cexp")   Cmplx c_cexp(Cmplx z);
    pragma(mangle, "clog")   Cmplx c_clog(Cmplx z);
    pragma(mangle, "csqrt")  Cmplx c_csqrt(Cmplx z);
    pragma(mangle, "conj")   Cmplx c_conj(Cmplx z);
    pragma(mangle, "cabs")   double c_cabs(Cmplx z);
}

alias CFuncPtr = extern (C) Cmplx function(Cmplx);

CFuncPtr get_complex_counterpart(void * func)
{
    if (func == cast(void *) &acos) return &c_cacos;
    if (func == cast(void *) &asin) return &c_casin;
    if (func == cast(void *) &atan) return &c_catan;
    if (func == cast(void *) &cos)  return &c_ccos;
    if (func == cast(void *) &sin)  return &c_csin;
    if (func == cast(void *) &tan)  return &c_ctan;
    if (func == cast(void *) &acosh) return &c_cacosh;
    if (func == cast(void *) &asinh) return &c_casinh;
    if (func == cast(void *) &atanh) return &c_catanh;
    if (func == cast(void *) &cosh) return &c_ccosh;
    if (func == cast(void *) &sinh) return &c_csinh;
    if (func == cast(void *) &tanh) return &c_ctanh;
    if (func == cast(void *) &exp)  return &c_cexp;
    if (func == cast(void *) &log)  return &c_clog;
    if (func == cast(void *) &sqrt) return &c_csqrt;

    return null;
}

bool cmplx_to_double(int f)
{
    return (f == GenFunc.f_carg || f == GenFunc.f_cmod ||
            f == GenFunc.f_real || f == GenFunc.f_imag ||
            f == GenFunc.f_cquad || f == GenFunc.f_abs);
}

/* application of scalar function to each element of matrix */

Node * apply_matrix_func(Node * t, Node * f, Parser * p)
{
    GretlMatrixRaw * m = null;
    GretlMatrixRaw tmp;
    int ret_complex = 0;
    Node * ret;

    if (t.t == GenSym.num) {
        m = &tmp;
        gretl_matrix_init(m);
        m.rows = m.cols = 1;
        m.val = &t.v.xval;
    } else {
        m = t.v.m;
        ret_complex = m.is_complex;
    }

    if (m.is_complex && cmplx_to_double(f.t)) {
        ret_complex = 0;
    }

    ret = aux_sized_matrix_node(p, m.rows, m.cols, ret_complex);
    if (ret == null) {
        return ret;
    }

    if (m.is_complex) {
        if (f.t == GenFunc.f_abs) {
            apply_cmatrix_dfunc(ret.v.m, m, &c_cabs);
        } else if (!ret_complex) {
            apply_cmatrix_dfunc(ret.v.m, m, cast(CDFunc) f.v.ptr);
        } else if (f.t == GenFunc.f_conj) {
            apply_cmatrix_cfunc(ret.v.m, m, &c_conj);
        } else if (f.t == GenSym.u_neg || f.t == GenSym.u_pos || f.t == GenSym.u_not) {
            apply_cmatrix_unary_op(ret.v.m, m, f.t);
        } else {
            CFuncPtr cfunc = null;

            if (f.v.ptr != null) {
                cfunc = get_complex_counterpart(f.v.ptr);
            }
            if (cfunc == null) {
                /* gatekeeper for complex */
                p.err = function_real_only(f.t);
            } else {
                p.err = apply_cmatrix_cfunc(ret.v.m, m, cfunc);
            }
        }
    } else if (f.t == GenFunc.f_real || f.t == GenFunc.f_imag) {
        if (f.t == GenFunc.f_real) {
            size_t msize = m.rows * m.cols * double.sizeof;

            memcpy(ret.v.m.val, m.val, msize);
        } else {
            int i, n = m.rows * m.cols;

            for (i=0; i<n; i++) {
                ret.v.m.val[i] = 0;
            }
        }
    } else {
        alias DFunc = extern (C) double function(double);
        DFunc dfunc = cast(DFunc) f.v.ptr;
        int i, n = m.rows * m.cols;

        if (dfunc != null) {
            for (i=0; i<n && !p.err; i++) {
                ret.v.m.val[i] = dfunc(m.val[i]);
            }
        } else {
            for (i=0; i<n && !p.err; i++) {
                ret.v.m.val[i] = real_apply_func(m.val[i], f.t, p);
            }
        }
    }

    return ret;
}

GretlMatrixRaw * matrix_from_scalars(Node ** nn, int m,
                                  int nsep, int seppos,
                                  Parser * p)
{
    GretlMatrixRaw * M;
    Node * n;
    int r = nsep + 1;
    int c = (seppos > 0)? seppos : m;
    int nelem = m - nsep;
    double x;
    int i, j, k;

    /* check that all rows are the same length */

    if (nelem != r * c) {
        p.err = GretlError.parse;
    } else if (nsep > 0) {
        k = 0;
        for (i=0; i<m; i++) {
            n = nn[i];
            if (null_node(n)) {
                if (i - k != seppos) {
                    p.err = GretlError.parse;
                    break;
                }
                k = i + 1;
            }
        }
    }

    if (p.err) {
        pprintf(p.prn, _("Matrix specification is not coherent").ptr);
        pputc(p.prn, '\n');
        return null;
    }

    M = gretl_matrix_alloc(r, c);
    if (M == null) {
        p.err = GretlError.alloc;
    } else {
        k = 0;
        for (i=0; i<r && !p.err; i++) {
            for (j=0; j<c; j++) {
                n = nn[k++];
                if (null_node(n)) {
                    n = nn[k++];
                }
                x = node_get_scalar(n, p);
                gretl_matrix_set(M, i, j, x);
            }
        }
    }

    return M;
}

int * full_series_list(const(DatasetRaw) * dset, int * err)
{
    int * list = null;

    if (dset.v < 2) {
        *err = GretlError.data;
        return null;
    }

    list = gretl_consecutive_list_new(1, dset.v - 1);
    if (list == null) {
        *err = GretlError.alloc;
        return null;
    }

    return list;
}

GretlMatrixRaw * real_matrix_from_list(const(int) * list,
                                    const(DatasetRaw) * dset,
                                    Parser * p)
{
    GretlMatrixRaw * M = null;

    if (list != null && list[0] == 0) {
        M = gretl_null_matrix_new();
    } else {
        const(GretlMatrixRaw) * mmask = get_matrix_mask();

        if (mmask != null) {
            M = gretl_matrix_data_subset_special(list, dset,
                                                 mmask, &p.err);
        } else {
            int missop = MMissingCode.ok;

            if (libset_get_bool(SetKey.skip_missing)) {
                missop = MMissingCode.skip;
            }
            M = gretl_matrix_data_subset(list, dset, dset.t1, dset.t2,
                                         missop, &p.err);
        }
    }

    return M;
}

GretlMatrixRaw * matrix_from_list(Node * n, Parser * p)
{
    GretlMatrixRaw * M = null;
    int * list = null;
    int freelist = 0;

    if (n != null) {
        if (n.t == GenSym.list) {
            list = n.v.ivec;
        } else {
            p.err = GretlError.data;
        }
    } else {
        list = full_series_list(p.dset, &p.err);
        freelist = 1;
    }

    if (!p.err) {
        M = real_matrix_from_list(list, p.dset, p);
    }

    if (freelist) {
        free(list);
    }

    return M;
}

void * arg_get_data(Node * n, int ref_, GretlType * type,
                    UserVar ** puv)
{
    void * data = null;

    *puv = n.uv;

    if (n.t == GenSym.series) {
        if (ref_) {
            *type = GretlType.series_ref;
            data = &n.vnum;
        } else if (n.vname != null) {
            /* FIXME conditionality here? */
            *type = GretlType.useries;
            data = &n.vnum;
        } else {
            *type = GretlType.series;
            data = n.v.xvec;
        }
    } else if (n.t == GenSym.num) {
        *type = ref_ ? GretlType.scalar_ref : GretlType.double_;
        data = &n.v.xval;
    } else if (n.t == GenSym.mat) {
        *type = ref_ ? GretlType.matrix_ref : GretlType.matrix;
        data = n.v.m;
    } else if (n.t == GenSym.bundle) {
        *type = ref_ ? GretlType.bundle_ref : GretlType.bundle;
        data = n.v.b;
    } else if (n.t == GenSym.array) {
        *type = ref_ ? GretlType.array_ref : GretlType.array;
        data = n.v.a;
    } else if (n.t == GenSym.str) {
        *type = ref_ ? GretlType.string_ref : GretlType.string_;
        data = n.v.str;
    } else if (n.t == GenSym.list) {
        *type = GretlType.list;
        data = n.v.ivec;
    } else if (n.t == GenSym.sub_addr) {
        data = sub_addr_get_data(n, type, puv);
    } else {
        *type = GretlType.none;
    }

    return data;
}

Node * suitable_ufunc_ret_node(Parser * p,
                               GretlType t)
{
    if (t == GretlType.double_) {
        return aux_scalar_node(p);
    } else if (t == GretlType.series) {
        return aux_empty_series_node(p);
    } else if (t == GretlType.matrix) {
        return aux_matrix_node(p);
    } else if (t == GretlType.list) {
        return aux_list_node(p);
    } else if (t == GretlType.string_) {
        return aux_string_node(p);
    } else if (t == GretlType.bundle) {
        return aux_bundle_node(p);
    } else if (gretl_array_type(t)) {
        return aux_array_node(p);
    } else {
        p.err = GretlError.types;
        return null;
    }
}

bool ok_ufunc_sym(int s)
{
    return (s == GenSym.num || s == GenSym.series || s == GenSym.mat ||
            s == GenSym.list || s == GenSym.u_addr || s == GenSym.dum ||
            s == GenSym.str || s == GenSym.empty || s == GenSym.bundle ||
            s == GenSym.array || s == GenSym.sub_addr);
}

/* evaluate a user-defined function */

Node * eval_ufunc(Node * l, Node * r, Parser * p)
{
    Node * ret = null;
    const(char) * funname = l.vname;
    Ufunc * uf = cast(Ufunc *) l.v.ptr;
    Fncall * fc = null;
    GretlType rtype = cast(GretlType) 0;
    int i, nparam, argc = 0;

    rtype = cast(GretlType) user_func_get_return_type(uf);

    if (!p.err && rtype == GretlType.void_) {
        if (p.targ == GenSym.unk && p.lh.name[0] == '\0' && p.lh.expr == null) {
            /* never reached? */
            p.targ = GenSym.empty;
        } else if (p.targ != GenSym.empty) {
            gretl_errmsg_sprintf(_("The function %s does not return any value").ptr,
                                 funname);
            p.err = GretlError.types;
        }
    }

    if (!p.err) {
        /* get the argument and param counts */
        argc = r.v.bn.n_nodes;
        nparam = fn_n_params(uf);
        if (argc > nparam) {
            gretl_errmsg_sprintf(_("Number of arguments (%d) does not match the number of\nparameters for function %s (%d)").ptr,
                                 argc, funname, nparam);
            p.err = GretlError.data;
        }
    }

    if (p.err) {
        /* no sense in continuing */
        return null;
    }

    if (l == p.tree && (p.flags & GenFlags.catch_)) {
        gretl_warnmsg_set(_("\"catch\" should not be used on calls to user-defined functions").ptr);
    }

    fc = user_func_get_fncall(uf);
    if (fc == null) {
        p.err = GretlError.alloc;
        return null;
    }

    /* check the function argument nodes */

    for (i=0; i<argc && !p.err; i++) {
        Node * arg = r.v.bn.n[i];
        int reftype = 0;

        if (!ok_ufunc_sym(arg.t)) {
            gretl_errmsg_sprintf("%s: invalid argument type %s".ptr, funname,
                                 typestr(arg.t));
            p.err = GretlError.types;
        }
        if (p.err) {
            break;
        }
        if (arg.t == GenSym.u_addr) {
            /* address node: switch to the 'content' sub-node */
            arg = ptr_node_get_referent_node(arg, p);
            if (!p.err) {
                reftype = 1;
            }
        } else if (arg.t == GenSym.dum && arg.v.idnum != DumId.null_) {
            p.err = GretlError.types;
        }
        if (!p.err) {
            /* assemble info and push argument */
            UserVar * uv = null;
            GretlType argt = cast(GretlType) 0;
            void * data;

            data = arg_get_data(arg, reftype, &argt, &uv);
            p.err = push_function_arg(fc, arg.vname, uv, argt, data);
        }
        if (p.err) {
            fprintf(stderr, "%s: error evaluating arg %d\n".ptr, funname, i+1);
        }
    }

    /* try sending args to function */

    if (!p.err) {
        char ** pdescrip = null;
        SeriesTable ** ptab = null;
        void * altp = null;
        void * retp = &altp;

        /* special cases */
        if (rtype == GretlType.void_) {
            retp = null;
        } else if (rtype == GretlType.matrix) {
            if (p.targ == GenSym.unk && p.tree == l) {
                /* target type not specified, and function returns
                   a matrix -> set target type to matrix
                */
                p.targ = GenSym.mat;
            }
        } else if (rtype == GretlType.list) {
            if (p.targ == GenSym.empty && p.tree == l) {
                /* this function offers a list return, but the
                   caller hasn't assigned it and it's not being
                   used as an argument to a further function, so
                   ignore the return value
                */
                retp = null;
            }
        }

        if ((p.flags & GenFlags.ufret) && rtype == GretlType.series) {
            /* arrange to pick up description of generated series, if any */
            pdescrip = &p.lh.label;
            /* and perhaps string values attached to series */
            ptab = &p.lh.stab;
        }

        p.err = gretl_function_exec_full(fc, rtype, p.dset, retp,
                                         pdescrip, ptab, p.prn);

        if (rtype == GretlType.numeric) {
            /* determine which numeric type we actually got */
            rtype = fncall_get_return_type(fc);
        }

        if (!p.err && retp != null) {
            ret = suitable_ufunc_ret_node(p, rtype);
        }

        if (!p.err && ret != null) {
            if (rtype == GretlType.double_) {
                ret.v.xval = *cast(double *) altp;
            } else if (rtype == GretlType.series) {
                if (ret.v.xvec != null) {
                    free(ret.v.xvec);
                }
                ret.v.xvec = cast(double *) altp;
            } else if (rtype == GretlType.matrix) {
                if (is_tmp_node(ret)) {
                    gretl_matrix_free(ret.v.m);
                }
                ret.v.m = cast(GretlMatrixRaw *) altp;
            } else if (rtype == GretlType.list) {
                if (is_tmp_node(ret)) {
                    free(ret.v.ivec);
                }
                if (altp != null) {
                    ret.v.ivec = cast(int *) altp;
                } else {
                    ret.v.ivec = gretl_list_new(0);
                }
            } else if (rtype == GretlType.string_) {
                if (is_tmp_node(ret)) {
                    free(ret.v.str);
                }
                ret.v.str = cast(char *) altp;
            } else if (rtype == GretlType.bundle) {
                if (is_tmp_node(ret)) {
                    gretl_bundle_destroy(ret.v.b);
                }
                ret.t = GenSym.bundle;
                ret.v.b = cast(GretlBundle *) altp;
            } else if (gretl_array_type(rtype)) {
                if (is_tmp_node(ret)) {
                    gretl_array_destroy(ret.v.a);
                }
                ret.t = GenSym.array;
                ret.v.a = cast(GretlArray *) altp;
            }
        }
    }

    return ret;
}

/* USE_RLIB is not defined in this port */

Node * eval_Rfunc(Node * t, Node * r, Parser * p)
{
    gretl_errmsg_set(_("R functions cannot be called via this gretl build").ptr);
    p.err = 1;
    return null;
}

/* Extracting what is nominally a series, stored in the form of a column
   vector, from a bundle on node @l. We assume a series is wanted on
   output unless the result is being directly assigned to a matrix.
   Whether we can truly get it as a series depends on the current state
   of the dataset. Failing that, we'll get it as a vector.
*/

Node * bundled_series_node(Node * l,
                           GretlMatrixRaw * m,
                           int is_virtual,
                           int * is_tmp,
                           Parser * p)
{
    Node * ret = null;
    double * y = null;
    int mt1 = gretl_matrix_get_t1(m);
    int mt2 = gretl_matrix_get_t2(m);
    int n = m.rows;
    int t, s;

    if (p.targ == GenSym.mat && l.parent == p.tree) {
        /* a matrix is wanted */
        ret = aux_matrix_node(p);
        if (!p.err) {
            ret.v.m = m;
            *is_tmp = is_virtual;
        }
        return ret;
    }

    if (n == p.dset.n || n == sample_size(p.dset) ||
        (mt1 >= 0 && mt2 >= mt1)) {
        ret = aux_series_node(p);
        if (!p.err) {
            y = ret.v.xvec;
        }
    }

    if (y != null) {
        /* transcribe to series */
        if (n == p.dset.n) {
            for (t=p.dset.t1; t<=p.dset.t2; t++) {
                y[t] = m.val[t];
            }
        } else if (n == sample_size(p.dset)) {
            for (t=p.dset.t1, s=0; t<=p.dset.t2; t++, s++) {
                y[t] = m.val[s];
            }
        } else {
            /* align using mt1 */
            for (t=mt1; t<mt1 + n && t<=p.dset.t2; t++) {
                if (t >= p.dset.t1) {
                    y[t] = m.val[t-mt1];
                }
            }
        }
        *is_tmp = 1;
    } else {
        /* vector fallback */
        ret = aux_matrix_node(p);
        if (!p.err) {
            ret.v.m = m;
            *is_tmp = is_virtual;
        }
    }

    return ret;
}

void attach_bundled_list(Node * n,
                         int * vals,
                         int * is_tmp,
                         Parser * p)
{
    if (stored_list_ok(vals, p.dset)) {
        /* extract as list */
        n.v.ivec = vals;
        *is_tmp = 0;
    } else {
        /* fallback: extract as row vector */
        n.t = GenSym.mat;
        n.v.m = gretl_list_to_vector(vals, &p.err);
        *is_tmp = 1;
    }
}

Node * virtual_object_node(Node * n,
                           const(char) * key,
                           GretlType type,
                           int * is_tmp,
                           Parser * p)
{
    GretlMatrixRaw * m;
    Node * ret = null;

    m = bundle_get_virtual_object(n.v.b, type, key, p.dset, &p.err);

    if (!p.err) {
        if (type == GretlType.series) {
            ret = bundled_series_node(n, m, 1, is_tmp, p);
            if (!p.err && ret.t == GenSym.series) {
                gretl_matrix_free(m);
            }
        } else {
            ret = aux_matrix_node(p);
            if (!p.err) {
                ret.v.m = m;
                *is_tmp = 1;
            }
        }
    }

    return ret;
}

/* Getting an object from within a bundle: on the left is the bundle
   reference, on the right should be a string: the key to look up to get
   content.
*/

Node * get_bundle_member(Node * l, Node * r, Parser * p)
{
    char * key = r.v.str;
    GretlType type;
    int querying;
    int is_virtual = 0;
    int is_tmp = 0;
    void * val = null;
    Node * ret = null;

    querying = (p.flags & GenFlags.objqry)? 1 : 0;

    if (querying && !gretl_bundle_has_key(l.v.b, key)) {
        ret = aux_empty_node(p);
        ret.flags |= NodeFlags.mut;
        return ret;
    }

    val = gretl_bundle_get_element(l.v.b, key, &type,
                                   &is_virtual, &p.err);
    if (p.err) {
        return ret;
    }

    if (querying) {
        ret = aux_empty_node(p);
        ret.t = cast(short) gen_type_from_gretl_type(type);
        return ret;
    }

    if (type != GretlType.series && !is_virtual) {
        ret = aux_node_for_type(type, p);
    }

    if (is_virtual) {
        ret = virtual_object_node(l, cast(const(char) *) val,
                                  type, &is_tmp, p);
    } else if (gretl_is_scalar_type(type)) {
        ret.v.xval = gretl_bundle_get_scalar(l.v.b, key, null);
    } else if (type == GretlType.string_) {
        ret.v.str = cast(char *) val;
    } else if (type == GretlType.matrix) {
        ret.v.m = cast(GretlMatrixRaw *) val;
    } else if (type == GretlType.series) {
        ret = bundled_series_node(l, cast(GretlMatrixRaw *) val, 0, &is_tmp, p);
    } else if (type == GretlType.bundle) {
        ret.v.b = cast(GretlBundle *) val;
    } else if (type == GretlType.array) {
        ret.v.a = cast(GretlArray *) val;
    } else if (type == GretlType.list) {
        attach_bundled_list(ret, cast(int *) val, &is_tmp, p);
    } else {
        p.err = GretlError.data;
    }

    if (ret != null) {
        if (is_tmp) {
            ret.flags |= NodeFlags.tmp;
        } else {
            ret.flags = cast(ubyte)(ret.flags & ~NodeFlags.tmp);
        }
        ret.flags |= NodeFlags.mut;
    }

    return ret;
}

/* implements the inbundle() function */

Node * test_bundle_key(Node * l, Node * r, Parser * p)
{
    Node * ret = aux_scalar_node(p);

    if (ret != null) {
        GretlBundle * bundle = l.v.b;
        const(char) * key = r.v.str;
        GretlType type;

        type = gretl_bundle_get_member_type(bundle, key, null);
        ret.v.xval = gretl_type_get_order(type);
    }

    return ret;
}

Node * get_bundle_array(Node * n, int f, Parser * p)
{
    Node * ret = aux_array_node(p);

    if (ret != null) {
        if (f == GenFunc.f_getkeys) {
            ret.v.a = gretl_bundle_get_keys(n.v.b, &p.err);
        } else {
            /* HF_JBTERMS */
            alias JFunc = extern (C) GretlArray * function(GretlBundle *, int *);
            JFunc jfunc = cast(JFunc) get_plugin_function("json_bundle_get_terminals");

            if (jfunc == null) {
                p.err = GretlError.fopen;
            } else {
                ret.v.a = jfunc(n.v.b, &p.err);
            }
        }
    }

    return ret;
}

/* USE_CURL is not defined in this port */

Node * curl_bundle_node(Node * n, Parser * p)
{
    Node * ret = aux_scalar_node(p);

    gretl_errmsg_set(_("Internet access not supported").ptr);
    p.err = GretlError.data;

    return ret;
}

Node * function_depth_node(Node * n, Parser * p)
{
    Node * ret = aux_scalar_node(p);

    if (ret != null) {
        ret.v.xval = gretl_function_depth();
    }

    return ret;
}

Node * plot_colors_node(Node * n, Parser * p)
{
    Node * ret = null;
    double f = 0.0;

    if (!null_node(n)) {
        f = node_get_scalar(n, p);
        if (na(f) || f < -1.0 || f > 1.0) {
            p.err = GretlError.invarg;
        }
    }
    if (!p.err) {
        ret = aux_array_node(p);
    }
    if (ret != null) {
        ret.v.a = get_user_colors(f, &p.err);
    }

    return ret;
}

Node * lpsolve_bundle_node(Node * n, Parser * p)
{
    Node * ret = aux_bundle_node(p);

    if (ret != null) {
        alias LpFunc = extern (C) GretlBundle * function(GretlBundle *, Prn *, int *);
        LpFunc lpfunc = cast(LpFunc) get_plugin_function("gretl_lpsolve");

        if (lpfunc == null) {
            p.err = GretlError.fopen;
        } else {
            ret.v.b = lpfunc(n.v.b, p.prn, &p.err);
        }
    }

    return ret;
}

GretlBundle * node_get_bundle(Node * n, Parser * p)
{
    GretlBundle * b = null;

    if (n.t == GenSym.bundle) {
        b = n.v.b;
    } else if (n.t == GenSym.u_addr) {
        n = ptr_node_get_referent_node(n, p);
        if (!p.err) {
            if (n.t != GenSym.bundle) {
                p.err = e_types(n);
            } else {
                b = n.v.b;
            }
        }
    } else {
        p.err = e_types(n);
    }

    return b;
}

Node * svm_driver_node(Node * t, Node * n, Parser * p)
{
    Node * e;
    Node * ret = null;
    int * list = null;
    GretlBundle * bparm = null;
    GretlBundle * bmod = null;
    GretlBundle * bprob = null;
    int i, k = n.v.bn.n_nodes;

    if (k < 2 || k > 4) {
        n_args_error(k, 2, 4, GenFunc.f_svm, p);
    }

    for (i=0; i<k && !p.err; i++) {
        e = n.v.bn.n[i];
        if (i == 0) {
            list = node_get_list(e, p);
        } else if (i == 1) {
            bparm = node_get_bundle(e, p);
        } else if (i == 2) {
            if (!null_node(e)) {
                bmod = node_get_bundle(e, p);
            }
        } else {
            if (!null_node(e)) {
                bprob = node_get_bundle(e, p);
            }
        }
    }

    if (!p.err) {
        ret = aux_series_node(p);
    }

    if (ret != null) {
        alias PFunc = extern (C) int function(const(int) *, GretlBundle *,
                                              GretlBundle *, GretlBundle *,
                                              double *, int *, DatasetRaw *, Prn *);
        int got_yhat = 0;

        PFunc pfunc = cast(PFunc) get_plugin_function("gretl_svm_driver");
        if (pfunc == null) {
            p.err = GretlError.fopen;
        } else {
            p.err = pfunc(list, bparm, bmod, bprob, ret.v.xvec,
                          &got_yhat, p.dset, p.prn);
            if (!p.err && !got_yhat) {
                /* change the return type to scalar NA */
                free(ret.v.xvec);
                ret.t = GenSym.num;
                ret.v.xval = NADBL;
            }
        }
    }

    free(list);

    return ret;
}

GretlBundle * bvar_get_bundle(Node * n, Parser * p)
{
    GretlBundle * b = null;

    if (n.v.idnum == BundleDataIndex.model) {
        b = bundle_from_model(null, p.dset, &p.err);
    } else if (n.v.idnum == BundleDataIndex.system) {
        b = bundle_from_system(null, 0, p.dset, &p.err);
    } else if (n.v.idnum == BundleDataIndex.sysinfo) {
        GretlBundle * tmp = get_sysinfo_bundle(&p.err);

        if (!p.err) {
            b = gretl_bundle_copy(tmp, &p.err);
        }
    } else if (n.v.idnum == RetrievalIndex.result) {
        GretlType type = cast(GretlType) 0;
        void * ptr = get_last_result_data(&type, &p.err);

        if (type == GretlType.bundle) {
            b = cast(GretlBundle *) ptr;
        } else if (!p.err) {
            p.err = GretlError.types;
        }
    } else {
        p.err = GretlError.data;
    }

    return b;
}

Node * dollar_bundle_node(Node * n, Parser * p)
{
    Node * ret = aux_bundle_node(p);

    if (ret != null) {
        ret.v.b = bvar_get_bundle(n, p);
    }

    return ret;
}

Node * type_string_node(Node * n, Parser * p)
{
    Node * ret = aux_string_node(p);

    if (ret != null) {
        int t = node_get_int(n, p);

        if (!p.err) {
            const(char) * s = "null";

            if (t == 1) {
                s = "scalar";
            } else if (t == 2) {
                s = "series";
            } else if (t == 3) {
                s = "matrix";
            } else if (t == 4) {
                s = "string";
            } else if (t == 5) {
                s = "bundle";
            } else if (t == 6) {
                s = "array";
            } else if (t == 7) {
                s = "list";
            }

            ret.v.str = gretl_strdup(s);
            if (ret.v.str == null) {
                p.err = GretlError.alloc;
            }
        }
    }

    return ret;
}

double * scalar_to_series(Node * n, Parser * p)
{
    double * ret = null;
    int t;

    if (p.dset == null || p.dset.n == 0) {
        p.err = GretlError.nodata;
    } else {
        ret = cast(double *) malloc(p.dset.n * (*ret).sizeof);
        if (ret == null) {
            p.err = GretlError.alloc;
        } else {
            for (t=0; t<p.dset.n; t++) {
                if (t >= p.dset.t1 && t <= p.dset.t2) {
                    ret[t] = n.v.xval;
                } else {
                    ret[t] = NADBL;
                }
            }
        }
    }

    return ret;
}

/* We come here only when setting a bundle-member or an element
   of an array -- otherwise we use get_check_return_type().
   Note 2021-08-12: in principle we could allow the case
   where @spec is an array type and @rhs is the singular of
   that type, and support auto-promotion of (e.g.) a string
   to an array of strings. But I'm not sure that's a good idea.
*/

int lhs_type_check(GretlType spec, GretlType rhs, int t)
{
    int err = 0;

    if (spec != 0 && spec != rhs) {
        if (t == GenSym.bundle) {
            gretl_errmsg_sprintf(_("Expected %s but got %s").ptr,
                                 gretl_type_get_name(spec),
                                 gretl_type_get_name(rhs));
        } else {
            gretl_errmsg_sprintf(_("Specified type %s does not match array type %s").ptr,
                                 gretl_type_get_name(spec),
                                 gretl_type_get_name(rhs));
        }
        err = GretlError.types;
    }

    return err;
}

void * get_mod_assign_result(void * lp, GretlType ltype,
                             Node * r, Parser * p)
{
    void * ret = null;
    Node * l;
    Node * op;

    if (p.op == GenSym.inc || p.op == GenSym.dec) {
        /* handle increment/decrement postfix operator */
        if (ltype == GretlType.double_) {
            double x = *cast(double *) lp;

            if (!na(x)) {
                x += (p.op == GenSym.inc)? 1 : -1;
                *cast(double *) lp = x;
            }
            ret = lp;
        } else {
            p.err = GretlError.types;
        }
        return ret; /* handled */
    }

    /* create binary tree to hold left and right */
    l = newempty();
    op = newb2(p.op, l, r);

    /* if that went OK, put the relevant type specifier
       and pointer onto the LHS sub-node
    */
    if (op == null || l == null) {
        p.err = GretlError.alloc;
    } else if (ltype == GretlType.matrix) {
        l.t = GenSym.mat;
        l.v.m = cast(GretlMatrixRaw *) lp;
    } else if (ltype == GretlType.double_) {
        l.t = GenSym.num;
        l.v.xval = *cast(double *) lp;
    } else if (ltype == GretlType.string_) {
        l.t = GenSym.str;
        l.v.str = cast(char *) lp;
    } else if (ltype == GretlType.bundle) {
        l.t = GenSym.bundle;
        l.v.b = cast(GretlBundle *) lp;
    } else if (ltype == GretlType.array) {
        l.t = GenSym.array;
        l.v.a = cast(GretlArray *) lp;
    } else if (ltype == GretlType.list) {
        l.t = GenSym.list;
        l.v.ivec = cast(int *) lp;
        if (p.op == GenSym.b_add) {
            /* reinterpret '+' when appending to a list */
            op.t = GenSym.b_lcat;
        }
    } else if (ltype == GretlType.series) {
        l.t = GenSym.series;
        l.v.xvec = cast(double *) lp;
    } else {
        p.err = GretlError.types;
    }

    if (!p.err) {
        /* FIXME parser state variables? */
        GenFlags saveflags = p.flags;
        int savetarg = p.targ;
        Node * ev;

        p.targ = l.t;
        p.flags = GenFlags.start;
        if (saveflags & GenFlags.listdef) {
            p.flags |= GenFlags.listdef;
        }
        ev = eval(op, p);

        if (!p.err) {
            /* get @ret off node @ev and clean up */
            if (ev.t == GenSym.mat) {
                ret = ev.v.m;
            } else if (ev.t == GenSym.bundle) {
                ret = ev.v.b;
            } else if (ev.t == GenSym.str) {
                ret = ev.v.str;
            } else if (ev.t == GenSym.array) {
                ret = ev.v.a;
            } else if (ev.t == GenSym.list) {
                ret = ev.v.ivec;
            } else if (ev.t == GenSym.series) {
                ret = ev.v.xvec;
            } else if (ev.t == GenSym.num) {
                ret = lp;
                *cast(double *) lp = ev.v.xval;
            } else {
                p.err = e_types(ev);
            }
        }

        p.targ = savetarg;
        p.flags = saveflags;
        free(ev);
    }

    /* thought: if @p is reusable, should we try preserving the
       nodes allocated here?
    */

    /* trash temporary nodes */
    free(op);
    free(l);

    if (ret == null && !p.err) {
        p.err = GretlError.data;
    }

    return ret;
}

/* ".=" : we need a scalar (possibly complex) on the RHS */

int dot_assign_to_matrix(GretlMatrixRaw * m, Parser * p)
{
    Node * n = p.ret;
    int err = 0;

    if (scalar_node(n)) {
        gretl_matrix_fill(m, node_get_scalar(n, p));
    } else if (cscalar_node(n)) {
        err = gretl_cmatrix_fill(m, n.v.m.z[0]);
    } else {
        err = e_types(n);
    }

    return err;
}

bool empty_rhs_ok(int t)
{
    return (t == GretlType.bundle || gretl_is_array_type(cast(GretlType) t));
}

/* Setting an object in a bundle under a given key string. We get here
   only if p->lh.expr is non-NULL.
*/

int set_bundle_value(Node * lhs, Node * rhs, Parser * p)
{
    Node * lh1 = lhs.L;
    Node * lh2 = lhs.R;
    GretlType targ = cast(GretlType) 0;
    GretlType type = cast(GretlType) 0;
    GretlBundle * bundle;
    void * ptr = null;
    char * key = null;
    int donate = 0;
    int err = 0;

    if (lh1.t != GenSym.bundle) {
        return GretlError.data;
    } else if (lh2.t != GenSym.str && lh2.t != GenSym.mspec) {
        return GretlError.data;
    }

    bundle = lh1.v.b;
    key = lh2.t == GenSym.str ? lh2.v.str : lh2.v.mspec.lsel.str;

    if (bundle == null || key == null) {
        return GretlError.data;
    }

    if ((p.flags & GenFlags.priv) && p.op == GenSym.b_asn && null_node(rhs)) {
        /* this is an internal "special" that implements removal
           of a bundle member via the "delete" command
        */
        return gretl_bundle_delete_data(bundle, key);
    }

    if (p.op != GenSym.b_asn) {
        /* We must have an existing bundle member under @key, and
           its type will determine the type of the result of
           inflected assignment.
        */
        GretlType ltype = cast(GretlType) 0;
        void * lp;

        lp = gretl_bundle_get_target(bundle, key, &ltype, &err);
        if (!err) {
            targ = gretl_type_from_gen_type(p.targ);
            err = lhs_type_check(targ, ltype, GenSym.bundle);
        }
        if (p.op == GenSym.b_dotasn) {
            /* accepted only for matrices */
            if (!err) {
                if (ltype == GretlType.matrix) {
                    err = dot_assign_to_matrix(cast(GretlMatrixRaw *) lp, p);
                } else {
                    err = GretlError.types;
                }
            }
            return err; /* handled */
        }
        if (!err) {
            ptr = get_mod_assign_result(lp, ltype, rhs, p);
            err = p.err;
            if (p.op == GenSym.inc || p.op == GenSym.dec) {
                return err; /* handled */
            }
        }
        if (!err) {
            type = ltype;
            if (ptr != lp) {
                donate = 1; /* donate: is this always right? */
            }
        }
        goto push_data;
    }

    /* Note: @targ is the gretl type specified by the caller for the
       bundle member (if any, this need not be supplied), and @type is
       the gretl type of the object arising on the RHS.  It's an error
       if @targ is non-zero and @type does not agree with it -- except
       for the case where @targ is given as "series" and we get a
       suitable matrix on the right. As of 2015-10-03, when we get a
       request to put a series into a bundle we actually put in a
       matrix, which in fact makes it easier to get a series back out
       again.
    */

    if (p.targ == GenSym.array) {
        /* This may be either the type of a pre-existing LHS array, or
           a specific array-type passed to generate().
        */
        targ = p.lh.gtype;
    } else {
        targ = gretl_type_from_gen_type(p.targ);
    }

    if (targ == GretlType.none && null_node(rhs)) {
        /* at this point @targ is indeterminate, but maybe there's
           an existing member to fix its value?
        */
        void * lp = gretl_bundle_get_target(bundle, key, &targ, null);

        if (targ == GretlType.array) {
            targ = gretl_array_get_type(cast(GretlArray *) lp);
        }
    }

    if (!err && targ == GretlType.list) {
        ptr = node_get_list(rhs, p);
        err = p.err;
        if (!err) {
            type = GretlType.list;
            donate = 1;
        }
    } else if (!err) {
        switch (rhs.t) {
        case GenSym.num:
            if (targ == GretlType.series) {
                ptr = scalar_to_series(rhs, p);
                if (p.err) {
                    err = p.err;
                } else {
                    type = GretlType.series;
                    donate = 1;
                }
            } else if (targ == GretlType.matrix) {
                ptr = gretl_matrix_from_scalar(rhs.v.xval);
                type = GretlType.matrix;
                donate = 1;
            } else {
                ptr = &rhs.v.xval;
                type = GretlType.double_;
            }
            break;
        case GenSym.str:
            ptr = rhs.v.str;
            type = GretlType.string_;
            donate = !reusable(p) && is_tmp_node(rhs);
            break;
        case GenSym.mat:
            if (targ == GretlType.double_ && scalar_matrix_node(rhs)) {
                ptr = &rhs.v.m.val[0];
                type = GretlType.double_;
            } else if (targ == GretlType.series) {
                ptr = copy_colvec_as_series(rhs, p);
                if (!p.err) {
                    type = GretlType.series;
                    donate = 1;
                } else {
                    err = GretlError.types;
                }
            } else {
                ptr = rhs.v.m;
                type = GretlType.matrix;
                donate = is_tmp_node(rhs);
            }
            break;
        case GenSym.series:
            ptr = rhs.v.xvec;
            type = GretlType.series;
            donate = !reusable(p) && is_tmp_node(rhs);
            break;
        case GenSym.bundle:
            ptr = rhs.v.b;
            type = GretlType.bundle;
            break;
        case GenSym.array:
            ptr = rhs.v.a;
            /* get more specific type for comparison with
               what the user specified (if anything)
            */
            type = gretl_array_get_type(rhs.v.a);
            donate = !reusable(p) && is_tmp_node(rhs);
            break;
        case GenSym.list:
            if (targ == GretlType.matrix) {
                ptr = gretl_list_to_vector(rhs.v.ivec, &p.err);
                if (!p.err) {
                    type = GretlType.matrix;
                    donate = 1;
                }
            } else {
                ptr = rhs.v.ivec;
                type = GretlType.list;
                donate = is_tmp_node(rhs);
            }
            break;
        case GenSym.empty:
            /* "null" is OK as (re-)initializer for bundle or array */
            if (empty_rhs_ok(targ)) {
                if (targ == GretlType.bundle) {
                    ptr = gretl_bundle_new();
                    if (ptr == null) {
                        err = GretlError.alloc;
                    }
                } else {
                    ptr = gretl_array_new(targ, 0, &err);
                }
                if (!err) {
                    type = targ;
                    donate = 1;
                }
            } else {
                err = GretlError.types;
            }
            break;
        default:
            err = GretlError.types;
            break;
        }
    }

    if (!err) {
        /* check for result type-incompatible with user's spec */
        err = lhs_type_check(targ, type, GenSym.bundle);
    }

 push_data:

    if (!err) {
        if (gretl_is_array_type(type)) {
            /* revert to generic array type for the functions below */
            type = GretlType.array;
        }
        if (donate) {
            /* it's OK to hand over the data pointer */
            if (type == GretlType.series) {
                err = gretl_bundle_donate_series(bundle, key, cast(double *) ptr, p.dset);
            } else {
                err = gretl_bundle_donate_data(bundle, key, ptr, type);
            }
            if (ptr == rhs.v.ptr) {
                rhs.v.ptr = null; /* avoid freeing! */
            }
        } else {
            /* the data must be copied into the bundle */
            if (type == GretlType.series) {
                err = gretl_bundle_set_series(bundle, key, cast(double *) ptr, p.dset);
            } else {
                err = gretl_bundle_set_data(bundle, key, ptr, type);
            }
        }
        if (!err && type == GretlType.matrix) {
            /* for use by genr_get_output_matrix() */
            p.lh.mret = cast(GretlMatrixRaw *) ptr;
        }
    }

    return err;
}

int set_array_value(Node * lhs, Node * rhs, Parser * p)
{
    Node * lh1 = lhs.L;
    Node * lh2 = lhs.R;
    GretlType atype = cast(GretlType) 0;
    GretlType type = cast(GretlType) 0;
    GretlType targ = cast(GretlType) 0;
    GretlArray * array = null;
    void * ptr = null;
    int idx = 0;
    int donate = 0;
    int err = 0;

    if (lh1.t != GenSym.array) {
        return GretlError.types;
    }

    if (lh2.t == GenSym.mspec) {
        /* FIXME allow a range here? */
        idx = mspec_get_array_index(lh2.v.mspec, &err);
        if (err) {
            return err;
        }
    } else {
        idx = cast(int) lh2.v.xval;
    }

    array = lh1.v.a;
    if (array == null) {
        return GretlError.data;
    } else if (idx <= 0 || idx > gretl_array_get_length(array)) {
        gretl_errmsg_sprintf(_("Index value %d is out of bounds").ptr, idx);
        return GretlError.data;
    }

    atype = gretl_array_get_content_type(array);
    targ = gretl_type_from_gen_type(p.targ);
    err = lhs_type_check(targ, atype, GenSym.array);

    idx--; /* convert index to 0-based */

    if (!err && p.op != GenSym.b_asn) {
        GretlType ltype = cast(GretlType) 0;
        void * lp;

        lp = gretl_array_get_element(array, idx, &ltype, &err);
        if (p.op == GenSym.b_dotasn) {
            if (!err) {
                if (ltype == GretlType.matrix) {
                    err = dot_assign_to_matrix(cast(GretlMatrixRaw *) lp, p);
                } else {
                    err = GretlError.types;
                }
            }
            return err; /* handled */
        }
        if (!err) {
            ptr = get_mod_assign_result(lp, ltype, rhs, p);
            err = p.err;
            if (p.op == GenSym.inc || p.op == GenSym.dec) {
                return err; /* handled */
            }
        }
        if (!err) {
            type = ltype;
            if (ptr != lp) {
                donate = 1; /* donate: always right? */
            }
        }
        goto push_data;
    }

    if (!err && atype == GretlType.list) {
        ptr = node_get_list(rhs, p);
        err = p.err;
        if (!err) {
            type = GretlType.list;
            donate = 1;
        }
    } else if (!err) {
        switch (rhs.t) {
        case GenSym.num:
            if (atype == GretlType.matrix) {
                ptr = gretl_matrix_from_scalar(rhs.v.xval);
                type = GretlType.matrix;
                donate = 1;
            }
            break;
        case GenSym.str:
            ptr = rhs.v.str;
            type = GretlType.string_;
            donate = !reusable(p) && is_tmp_node(rhs);
            break;
        case GenSym.mat:
            ptr = rhs.v.m;
            type = GretlType.matrix;
            donate = is_tmp_node(rhs);
            break;
        case GenSym.bundle:
            ptr = rhs.v.b;
            type = GretlType.bundle;
            break;
        case GenSym.list:
            ptr = rhs.v.ivec;
            type = GretlType.list;
            donate = is_tmp_node(rhs);
            break;
        case GenSym.array:
            ptr = rhs.v.a;
            type = GretlType.array;
            donate = is_tmp_node(rhs);
            break;
        default:
            err = GretlError.types;
            break;
        }
    }

    if (!err && type != atype) {
        err = GretlError.types;
    }

 push_data:

    if (!err) {
        if (donate) {
            /* it's OK to hand over the data pointer */
            err = gretl_array_set_element(array, idx, ptr, type, 0);
            if (ptr == rhs.v.ptr) {
                rhs.v.ptr = null; /* gone! */
            }
        } else {
            /* the data must be copied into the array */
            err = gretl_array_set_element(array, idx, ptr, type, 1);
        }
        if (!err && type == GretlType.matrix) {
            /* for use by genr_get_output_matrix() */
            p.lh.mret = cast(GretlMatrixRaw *) ptr;
        }
    }

    return err;
}

/* setting member of list: only straight assignment is accepted */

int set_list_value(Node * lhs, Node * rhs, Parser * p)
{
    Node * lh1 = lhs.L;
    Node * lh2 = lhs.R;
    int * list = null;
    int idx = 0, v = -1;
    int err = 0;

    if (p.op != GenSym.b_asn) {
        gretl_errmsg_sprintf(_("'%s' : not implemented for this type").ptr,
                             get_opstr(p.op));
        return GretlError.types;
    }

    if (lh2.t == GenSym.mspec) {
        idx = mspec_get_array_index(lh2.v.mspec, &err);
        if (err) {
            return err;
        }
    } else {
        idx = cast(int) lh2.v.xval;
    }

    list = lh1.v.ivec;
    if (list == null) {
        return GretlError.data;
    } else if (idx < 1 || idx > list[0]) {
        gretl_errmsg_sprintf(_("Index value %d is out of bounds").ptr, idx);
        return GretlError.data;
    }

    if (rhs.t == GenSym.num) {
        v = node_get_int(rhs, p);
    } else if (rhs.t == GenSym.series) {
        v = rhs.vnum;
    } else {
        p.err = e_types(rhs);
    }

    if (!p.err) {
        if (v < 0 || v >= p.dset.v) {
            gretl_errmsg_set(_("Invalid list element").ptr);
            p.err = GretlError.data;
        }
    }

    if (!p.err) {
        list[idx] = v;
    }

    return p.err;
}

/* setting element of string: only straight assignment accepted */

int set_string_value(Node * lhs, Node * rhs, Parser * p)
{
    Node * lh1 = lhs.L;
    Node * lh2 = lhs.R;
    char * s1 = null;
    const(char) * s2 = null;
    int idx = 0, err = 0;

    if (p.op != GenSym.b_asn) {
        gretl_errmsg_sprintf(_("'%s' : not implemented for this type").ptr,
                             get_opstr(p.op));
        return GretlError.types;
    } else if (rhs.t != GenSym.str) {
        return e_types(rhs);
    }

    if (lh2.t == GenSym.mspec) {
        idx = mspec_get_array_index(lh2.v.mspec, &err);
        if (err) {
            return err;
        }
    } else {
        idx = cast(int) lh2.v.xval;
    }

    /* @s1 may be modified below */
    s1 = lh1.v.str;

    if (s1 == null) {
        return GretlError.data;
    } else if (idx < 1 || idx > g_utf8_strlen(s1, -1)) {
        gretl_errmsg_sprintf(_("Index value %d is out of bounds").ptr, idx);
        return GretlError.data;
    }

    s2 = rhs.v.str;
    if (g_utf8_strlen(s2, -1) != 1) {
        return GretlError.invarg;
    } else if (g_utf8_strlen(s1, -1) == strlen(s1) && strlen(s2) == 1) {
        /* simple: no multibyte characters */
        s1[idx-1] = s2[0];
    } else {
        /* handle the multibyte case */
        char * tmp = gretl_utf8_replace_char(s1, s2, idx - 1);

        if (strlen(tmp) <= strlen(s1)) {
            strcpy(s1, tmp);
            free(tmp);
        } else {
            UserVar * uv = get_user_var_by_data(s1);

            if (uv != null) {
                p.err = user_var_replace_value(uv, tmp,
                                               GretlType.string_);
            } else {
                p.err = GretlError.data;
                free(tmp);
            }
        }
    }

    return p.err;
}

int set_series_obs_value(Node * lhs, Node * rhs, Parser * p)
{
    Node * lh1 = lhs.L;
    Node * lh2 = lhs.R;
    double ** Z = p.dset.Z;
    char * label = null;
    double x = NADBL;
    int op = p.op;
    int v, t;

    if (lh1.t == GenSym.series && (lh2.t == GenSym.num || lh2.t == GenSym.mspec)) {
        ; /* OK */
    } else {
        return GretlError.types;
    }

    v = lh1.vnum;
    if (v <= 0 || v >= p.dset.v) {
        return GretlError.data;
    } else if (object_is_const(null, v)) {
        return overwrite_err(p.dset.varname[v]);
    }

    if (lh2.t == GenSym.mspec) {
        t = mspec_get_series_index(lh2.v.mspec, p);
    } else {
        t = node_get_int(lh2, p);
    }
    if (t < 1 || t > p.dset.n) {
        return GretlError.data;
    }
    /* convert to 0-based */
    t--;

    if (rhs == null) {
        if (p.op == GenSym.inc) {
            x = 1;
            op = GenSym.b_add;
        } else if (p.op == GenSym.dec) {
            x = 1;
            op = GenSym.b_sub;
        } else {
            return GretlError.types; /* ? */
        }
    } else if (rhs.t == GenSym.str) {
        if (is_string_valued(p.dset, v)) {
            label = rhs.v.str;
        } else {
            return GretlError.types;
        }
    } else {
        x = node_get_scalar(rhs, p);
        if (p.err) {
            return p.err;
        }
    }

    if (is_string_valued(p.dset, v)) {
        if (label != null) {
            if (op != GenSym.b_asn) {
                p.err = GretlError.types;
            } else {
                p.err = series_set_string_val(p.dset, v, t, label);
            }
        } else {
            x = xy_calc(Z[v][t], x, op, GenSym.num, p);
            if (!p.err) {
                p.err = string_series_assign_value(p.dset, v, t, x);
            }
        }
    } else {
        Z[v][t] = xy_calc(Z[v][t], x, op, GenSym.num, p);
    }

    if (p.err == 0) {
        /* made a change to an element of a series */
        p.flags |= GenFlags.obsval;
        set_dataset_is_changed(p.dset, 1);
    }

    return p.err;
}

/* Here we're replacing a submatrix, by either straight or inflected
   assignment.

   @lhs must be a binary node holding the target matrix on its L
   branch and a matrix subspec on its R branch.  @rhs must hold the
   replacement value: either a matrix or a scalar (or a series
   standing in for a matrix).
*/

int set_matrix_chunk(Node * lhs, Node * rhs, Parser * p)
{
    Node * lh1 = lhs.L;
    Node * lh2 = lhs.R;
    GretlMatrixRaw * m1;
    GretlMatrixRaw * m2 = null;
    MatrixSubspec * spec;
    double rhs_x = NADBL;
    Cmplx rhs_z = Cmplx(NADBL, 0);
    int rhs_scalar = 0;
    int rhs_cscalar = 0;
    int inflected = 0;
    int free_m2 = 0;

    if (p.op == GenSym.b_hcat || p.op == GenSym.b_vcat) {
        /* can't do these things on a submatrix */
        gretl_errmsg_sprintf(_("The operator '%s' is not valid in this context").ptr,
                             get_opstr(p.op));
        return GretlError.types;
    } else if (lh1.t != GenSym.mat) {
        /* is this ever possible? */
        fprintf(stderr, "set_matrix_chunk: got %s, not matrix!\n".ptr,
                getsymb(lh1.t));
        return GretlError.data;
    }

    /* set up the target */
    m1 = lh1.v.m;
    spec = lh2.v.mspec;
    if (m1 == null || spec == null) {
        return GretlError.data;
    }

    /* check the validity of the subspec we got, and adjust it if need
       be in the light of the dimensions of @m.
    */
    if (!spec.checked) {
        p.err = check_matrix_subspec(spec, m1);
        if (p.err) {
            fprintf(stderr, "set_matrix_chunk: check_matrix_subspec failed\n".ptr);
            return p.err;
        }
    }

    /* Is the assignment straight or inflected?  Note that in this
       context there's no distinction between '=' and '.=' and the
       latter doesn't count as inflected.
    */
    if (p.op != GenSym.b_asn && p.op != GenSym.b_dotasn) {
        inflected = 1;
    }

    if (p.op == GenSym.inc || p.op == GenSym.dec) {
        /* treat as add or subtract */
        rhs_x = 1;
        rhs_z = Cmplx(rhs_x, 0);
        rhs_scalar = 1;
    } else if (scalar_node(rhs)) {
        /* single value (could be 1 x 1 matrix) on RHS */
        rhs_x = (rhs.t == GenSym.num)? rhs.v.xval: rhs.v.m.val[0];
        rhs_z = Cmplx(rhs_x, 0);
        rhs_scalar = 1;
    } else if (cscalar_node(rhs)) {
        if (!m1.is_complex) {
            gretl_errmsg_set("Cannot assign complex values to a real matrix");
            p.err = GretlError.types;
        } else {
            m2 = rhs.v.m;
            rhs_z = rhs.v.m.z[0];
            rhs_cscalar = 1;
        }
    } else if (rhs.t == GenSym.mat) {
        /* not a scalar: get the RHS matrix */
        m2 = rhs.v.m;
        if (m2.is_complex && !m1.is_complex) {
            gretl_errmsg_set("Cannot assign complex values to a real matrix");
            p.err = GretlError.types;
        }
    } else if (rhs.t == GenSym.series) {
        /* legacy: this has long been accepted */
        m2 = series_to_matrix(rhs.v.xvec, p);
        free_m2 = 1; /* flag temporary status of @m2 */
    } else {
        p.err = GretlError.types;
    }

    if (p.err) {
        return p.err;
    }

    if (spec.ltype == SelType.element) {
        /* assignment, plain or inflected, to a single
           element of target matrix.
        */
        int i = mspec_get_element(spec);

        if (rhs_cscalar || (rhs_scalar && m1.is_complex)) {
            if (!inflected) {
                m1.z[i] = rhs_z;
            } else {
                m1.z[i] = c_xy_calc(m1.z[i], rhs_z, p.op, p);
            }
        } else if (rhs_scalar) {
            if (!inflected) {
                m1.val[i] = rhs_x;
            } else {
                m1.val[i] = xy_calc(m1.val[i], rhs_x, p.op, GenSym.mat, p);
            }
        } else {
            /* here the RHS must be 1 x 1 */
            p.err = GretlError.nonconf;
        }
        return p.err; /* we're done */
    }

    if (!inflected) {
        if (rhs_cscalar) {
            return assign_scalar_to_submatrix(m1, m2, 0, spec);
        } else if (rhs_scalar) {
            return assign_scalar_to_submatrix(m1, null, rhs_x, spec);
        } else if (is_sel_dummy(spec.ltype)) {
            return gretl_matrix_set_part(m1, m2, 0, spec.ltype);
        }
    }

    if (inflected) {
        /* Here we're doing '+=' or some such, in which case a new
           submatrix must be calculated using the original
           submatrix @a and the newly generated matrix (or
           scalar value).
        */
        GretlMatrixRaw * a = matrix_get_submatrix(m1, spec, 1, &p.err);

        if (!p.err) {
            if (rhs_scalar || rhs_cscalar) {
                if (a.is_complex) {
                    cmatrix_xy_calc(a, a, rhs_z, 0, p.op, p);
                } else {
                    rmatrix_xy_calc(a, a, rhs_x, 0, p.op, p);
                }
                m2 = a; /* assign computed matrix to m2 */
                free_m2 = 1;
            } else {
                GretlMatrixRaw * b = null;

                p.err = real_matrix_calc(a, m2, p.op, &b);
                gretl_matrix_free(a);
                /* replace RHS m2 with computed result */
                if (free_m2) {
                    /* m2 was temp result of series conversion */
                    gretl_matrix_free(m2);
                }
                m2 = b;
                free_m2 = 1;
            }
            /* we now proceed to matrix_replace_submatrix() */
        }
    }

    if (!p.err) {
        /* Write new submatrix @m2 into place: note that we come here
           directly if none of the special conditions above are
           satisfied -- for example, if the newly generated value
           is a matrix and the task is straight assignment. Also
           check for numerical "breakage" in the replacement
           submatrix.
        */
        p.err = matrix_replace_submatrix(m1, m2, spec);
    }

    if (free_m2) {
        gretl_matrix_free(m2);
    }

    return p.err;
}

GretlMatrixRaw * get_corrgm_matrix(Node * l,
                                Node * m,
                                Node * r,
                                Parser * p)
{
    int xcf = (r.t != GenSym.empty);
    int * list = null;
    GretlMatrixRaw * A = null;
    int k;

    /* ensure we've got an order */
    k = node_get_int(m, p);
    if (p.err) {
        return null;
    }

    /* hook up list if arg1 is list */
    if (l.t == GenSym.list) {
        list = l.v.ivec;
    }

    /* if third node is matrix, must be real col vector */
    if (r.t == GenSym.mat) {
        if (r.v.m.cols != 1) {
            p.err = GretlError.nonconf;
            return null;
        } else if (r.v.m.is_complex) {
            p.err = GretlError.nonconf;
            return null;
        }
    }

    if (!xcf) {
        /* acf/pacf */
        if (l.t == GenSym.series) {
            A = acf_matrix(l.v.xvec, k, p.dset, 0, &p.err);
        } else if (l.t == GenSym.mat) {
            A = multi_acf(l.v.m, null, null, k, &p.err);
        } else {
            /* it must be a list */
            A = multi_acf(null, list, p.dset, k, &p.err);
        }
    } else {
        /* cross-correlogram */
        void * px = null;
        void * py = null;
        int xtype = GenSym.series;

        if (list != null) {
            px = list;
            xtype = GenSym.list;
        } else if (l.t == GenSym.mat) {
            px = l.v.m;
            xtype = GenSym.mat;
        } else {
            px = l.v.xvec;
        }

        py = (r.t == GenSym.mat)? cast(void *) r.v.m : cast(void *) r.v.xvec;

        A = multi_xcf(px, xtype, py, r.t, p.dset, k, &p.err);
    }

    return A;
}

GretlMatrixRaw * get_density_matrix(Node * t, double bws,
                                 int ctrl, Parser * p)
{
    alias KdFunc1 = extern (C) GretlMatrixRaw * function(const(double) *, int, double,
                                                      Gretlopt, int *);
    alias KdFunc2 = extern (C) GretlMatrixRaw * function(const(GretlMatrixRaw) *, double,
                                                      Gretlopt, int *);
    KdFunc1 kdfunc1 = null;
    KdFunc2 kdfunc2 = null;
    Gretlopt opt = ctrl ? Gretlopt.OPT_O : Gretlopt.OPT_NONE;
    GretlMatrixRaw * m = null;
    GretlMatrixRaw * X = null;
    const(double) * x = null;
    int free_X = 0;
    int n = 0;

    if (t.t == GenSym.series) {
        n = sample_size(p.dset);
        x = t.v.xvec + p.dset.t1;
    } else if (t.t == GenSym.list) {
        X = gretl_matrix_data_subset(t.v.ivec, p.dset,
                                     p.dset.t1, p.dset.t2,
                                     MMissingCode.skip, &p.err);
        free_X = 1;
    } else {
        /* matrix */
        if (t.v.m.is_complex) {
            p.err = GretlError.cmplx;
        } else {
            n = gretl_vector_get_length(t.v.m);
            if (n > 0) {
                /* vector */
                x = t.v.m.val;
            } else {
                X = t.v.m;
            }
        }
    }

    if (!p.err) {
        if (X != null) {
            kdfunc2 = cast(KdFunc2) get_plugin_function("multiple_kd_matrix");
        } else if (!p.err) {
            kdfunc1 = cast(KdFunc1) get_plugin_function("kernel_density_matrix");
        }
        if (kdfunc2 == null && kdfunc1 == null) {
            p.err = GretlError.fopen;
        } else if (X != null) {
            m = kdfunc2(X, bws, opt, &p.err);
        } else {
            m = kdfunc1(x, n, bws, opt, &p.err);
        }
    }

    if (free_X) {
        gretl_matrix_free(X);
    }

    return m;
}

int aggregate_discrete_check(const(int) * list, const(DatasetRaw) * dset)
{
    int i, vi;

    for (i=1; i<=list[0]; i++) {
        vi = list[i];
        if (!accept_as_discrete(dset, vi, 0)) {
            gretl_errmsg_sprintf(_("The variable '%s' is not discrete").ptr,
                                 dset.varname[vi]);
            return GretlError.data;
        }
    }

    return 0;
}

GretlMatrixRaw * mshape_scalar(double x, int r, int c, int * err)
{
    GretlMatrixRaw * m = gretl_matrix_alloc(r, c);

    if (m == null) {
        *err = GretlError.alloc;
    } else {
        int i, n = r * c;

        for (i=0; i<n; i++) {
            m.val[i] = x;
        }
    }

    return m;
}

void node_get_int_or_series(int * ip, double ** vecp,
                            Node * n, Parser * p)
{
    if (scalar_node(n)) {
        *ip = node_get_int(n, p);
    } else if (n.t == GenSym.series) {
        *vecp = n.v.xvec;
    } else {
        p.err = e_types(n);
    }
}

/* eval_3args_func: evaluate a built-in function that has three
   arguments. The @post_process flag is a convenience for the
   case where the function in question returns a matrix: it
   centralizes the creation or retrieval of an "aux node" of the
   right type, and attaches the computed matrix @A to it. This
   flag must be set to zero for all cases where the function does
   NOT return a matrix. Conversely, when a function DOES return
   a matrix it should assign this to @A, and leave the aux node
   business to the post-processor.
*/

Node * eval_3args_func(Node * l, Node * m, Node * r,
                       int f, Parser * p)
{
    Node * save_aux = p.aux;
    Node * ret = null;
    GretlMatrixRaw * A = null;
    int post_process = 1;

    if (f == GenFunc.f_mshape) {
        if (l.t != GenSym.mat && l.t != GenSym.num) {
            node_type_error(f, 1, GenSym.mat, l, p);
        } else if (!scalar_node(m)) {
            node_type_error(f, 2, GenSym.num, m, p);
        } else if (!null_node(r) && !scalar_node(r)) {
            node_type_error(f, 3, GenSym.num, r, p);
        } else {
            int n, k2, k1 = node_get_int(m, p);

            if (scalar_node(r)) {
                k2 = node_get_int(r, p);
            } else if (l.t == GenSym.num) {
                k2 = 1;
            } else {
                n = l.v.m.rows * l.v.m.cols;
                if (n % k1 == 0) {
                    k2 = n / k1;
                } else {
                    p.err = GretlError.invarg;
                }
            }
            if (!p.err) {
                if (l.t == GenSym.num) {
                    A = mshape_scalar(l.v.xval, k1, k2, &p.err);
                } else {
                    A = gretl_matrix_shape(l.v.m, k1, k2, &p.err);
                }
            }
        }
    } else if (f == GenFunc.f_trimr) {
        if (l.t != GenSym.mat) {
            node_type_error(f, 1, GenSym.mat, l, p);
        } else if (!scalar_node(m)) {
            node_type_error(f, 2, GenSym.num, m, p);
        } else if (!scalar_node(r)) {
            node_type_error(f, 3, GenSym.num, r, p);
        } else {
            int k1 = node_get_int(m, p);
            int k2 = node_get_int(r, p);

            if (!p.err) {
                A = gretl_matrix_trim_rows(l.v.m, k1, k2, &p.err);
            }
        }
    } else if (f == GenFunc.f_svd) {
        /* note: the complex case is supported */
        GretlMatrixRaw * lm = node_get_matrix(l, p, 0, 1);
        GretlMatrixRaw * U = null;
        GretlMatrixRaw * V = null;

        if (!p.err) {
            if (m.t == GenSym.u_addr) {
                U = ptr_node_get_matrix(m, p);
            } else if (m.t != GenSym.empty) {
                node_type_error(f, 2, GenSym.u_addr, m, p);
            }
        }
        if (!p.err) {
            if (r.t == GenSym.u_addr) {
                V = ptr_node_get_matrix(r, p);
            } else if (r.t != GenSym.empty) {
                node_type_error(f, 3, GenSym.u_addr, r, p);
            }
        }
        if (!p.err) {
            A = user_matrix_SVD(lm, U, V, &p.err);
        }
    } else if (f == GenFunc.f_qr) {
        GretlMatrixRaw * X = node_get_matrix(l, p, 0, 1);
        GretlMatrixRaw * R = null;
        GretlMatrixRaw * P = null;

        if (!p.err) {
            if (m.t == GenSym.u_addr) {
                R = ptr_node_get_matrix(m, p);
            } else if (m.t != GenSym.empty) {
                node_type_error(f, 2, GenSym.u_addr, m, p);
            }
        }
        if (!p.err) {
            if (r.t == GenSym.u_addr) {
                P = ptr_node_get_matrix(r, p);
            } else if (r.t != GenSym.empty) {
                node_type_error(f, 3, GenSym.u_addr, r, p);
            }
        }
        if (!p.err) {
            if (X.is_complex) {
                if (P != null) {
                    A = gretl_cmatrix_QR_pivot_decomp(X, R, P, &p.err);
                } else {
                    A = gretl_cmatrix_QR_decomp(X, R, &p.err);
                }
            } else {
                A = user_matrix_QR_decomp(X, R, P, &p.err);
            }
        }
    } else if (f == GenFunc.f_varsimul) {
        GretlMatrixRaw * m1 = node_get_real_matrix(l, p, 0, 1);
        GretlMatrixRaw * m2 = node_get_real_matrix(m, p, 1, 2);
        GretlMatrixRaw * m3 = node_get_real_matrix(r, p, 2, 3);

        if (!p.err) {
            A = gretl_matrix_varsimul(m1, m2, m3, &p.err);
        }
    } else if (f == GenFunc.f_eigen) {
        GretlMatrixRaw * lm = node_get_matrix(l, p, 0, 1);
        GretlMatrixRaw * v1 = null;
        GretlMatrixRaw * v2 = null;

        if (l.t != GenSym.mat) {
            node_type_error(f, 1, GenSym.mat, l, p);
        } else {
            if (!null_node(m)) {
                v1 = ptr_node_get_matrix(m, p);
            }
            if (!null_node(r)) {
                v2 = ptr_node_get_matrix(r, p);
            }
        }
        if (!p.err) {
            if (lm.is_complex) {
                A = gretl_zgeev(lm, v1, v2, &p.err);
            } else {
                A = gretl_dgeev(lm, v1, v2, &p.err);
            }
        }
    } else if (f == GenFunc.f_schur) {
        GretlMatrixRaw * Z = null;
        GretlMatrixRaw * W = null;

        if (l.t != GenSym.mat || !l.v.m.is_complex) {
            node_type_error(f, 1, GenSym.mat, l, p);
        } else {
            if (!null_node(m)) {
                Z = ptr_node_get_matrix(m, p);
            }
            if (!null_node(r)) {
                W = ptr_node_get_matrix(r, p);
            }
        }
        if (!p.err) {
            A = gretl_zgees(l.v.m, Z, W, &p.err);
        }
    } else if (f == GenFunc.f_corrgm) {
        if (l.t != GenSym.series && l.t != GenSym.mat && l.t != GenSym.list) {
            node_type_error(f, 1, GenSym.series, l, p);
        } else if (!scalar_node(m)) {
            node_type_error(f, 2, GenSym.num, m, p);
        } else if (r.t != GenSym.empty && r.t != GenSym.series && r.t != GenSym.mat) {
            node_type_error(f, 3, GenSym.series, r, p);
        } else {
            A = get_corrgm_matrix(l, m, r, p);
        }
    } else if (f == GenFunc.f_seq) {
        if (!scalar_node(l)) {
            node_type_error(f, 1, GenSym.num, l, p);
        } else if (!scalar_node(m)) {
            node_type_error(f, 2, GenSym.num, m, p);
        } else if (!null_or_scalar(r)) {
            node_type_error(f, 3, GenSym.num, r, p);
        } else {
            double start = node_get_scalar(l, p);
            double end = node_get_scalar(m, p);
            double step = (!null_node(r))? node_get_scalar(r, p) : 1.0;

            if (!p.err) {
                A = gretl_matrix_seq(start, end, step, &p.err);
            }
        }
    } else if (f == GenFunc.f_strncmp) {
        post_process = 0;
        if (l.t != GenSym.str) {
            node_type_error(f, 1, GenSym.str, l, p);
        } else if (m.t != GenSym.str) {
            node_type_error(f, 2, GenSym.str, m, p);
        } else if (!null_or_scalar(r)) {
            node_type_error(f, 3, GenSym.num, r, p);
        } else {
            ret = aux_scalar_node(p);
            if (ret != null) {
                if (!null_node(r)) {
                    int len = node_get_int(r, p);

                    if (!p.err) {
                        ret.v.xval = strncmp(l.v.str, m.v.str, len);
                    }
                } else {
                    ret.v.xval = strcmp(l.v.str, m.v.str);
                }
            }
        }
    } else if (f == GenFunc.f_weekday || f == GenFunc.f_isoweek) {
        int series_mode = 0;
        int basic_mode = 0;

        post_process = 0;

        if (null_node(m) && null_node(r)) {
            /* got a single argument */
            if (l.t == GenSym.series) {
                basic_mode = series_mode = 1;
            } else if (scalar_node(l)) {
                basic_mode = 1;
            } else {
                p.err = e_types(l);
            }
        } else if (l.t == GenSym.series && m.t == GenSym.series && r.t == GenSym.series) {
            series_mode = 1;
        } else if (!(scalar_node(l) && scalar_node(m) && scalar_node(r))) {
            p.err = GretlError.types;
        }

        if (!p.err && series_mode) {
            const(double) * mx = basic_mode ? null : m.v.xvec;
            const(double) * rx = basic_mode ? null : r.v.xvec;

            reset_p_aux(p, save_aux);
            ret = aux_series_node(p);
            if (ret != null && f == GenFunc.f_weekday) {
                p.err = fill_day_of_week_array(ret.v.xvec, l.v.xvec,
                                               mx, rx, p.dset);
            } else if (ret != null) {
                p.err = fill_isoweek_array(ret.v.xvec, l.v.xvec,
                                           mx, rx, p.dset);
            }
        } else if (!p.err) {
            ret = aux_scalar_node(p);
            if (ret != null) {
                int mo, day, yr = node_get_int(l, p);

                if (!p.err && basic_mode) {
                    char[16] tstr;
                    int n;

                    sprintf(tstr.ptr, "%d", yr);
                    n = sscanf(tstr.ptr, "%4d%2d%2d", &yr, &mo, &day);
                    if (n != 3) {
                        p.err = GretlError.invarg;
                    }
                } else if (!p.err) {
                    mo = node_get_int(m, p);
                    day = node_get_int(r, p);
                }
                if (!p.err && f == GenFunc.f_weekday) {
                    int julian = 0;

                    if (yr < 0) {
                        yr = -yr;
                        julian = 1;
                    }
                    ret.v.xval = legacy_day_of_week(yr, mo, day, julian, &p.err);
                } else if (!p.err) {
                    ret.v.xval = iso_week_number(yr, mo, day, &p.err);
                }
            }
        }
    } else if (f == GenFunc.f_dayspan) {
        guint32 ed1 = node_get_guint32(l, p);
        guint32 ed2 = node_get_guint32(m, p);
        int wkdays = node_get_int(r, p);

        post_process = 0;
        if (!p.err) {
            ret = aux_scalar_node(p);
        }
        if (!p.err) {
            ret.v.xval = day_span(ed1, ed2, wkdays, &p.err);
        }
    } else if (f == GenFunc.f_smplspan) {
        post_process = 0;
        if (l.t == GenSym.str && m.t == GenSym.str && r.t == GenSym.num) {
            ret = aux_scalar_node(p);
            if (ret != null) {
                int pd = node_get_int(r, p);

                if (!p.err) {
                    ret.v.xval = sample_span(l.v.str, m.v.str,
                                             pd, &p.err);
                }
            }
        } else {
            p.err = GretlError.types;
        }
    } else if (f == GenFunc.f_kdensity) {
        if (l.t != GenSym.series && l.t != GenSym.list && l.t != GenSym.mat) {
            node_type_error(f, 1, GenSym.series, l, p);
        } else if (m.t != GenSym.num && m.t != GenSym.empty) {
            node_type_error(f, 2, GenSym.num, m, p);
        } else if (r.t != GenSym.num && r.t != GenSym.empty) {
            node_type_error(f, 3, GenSym.num, r, p);
        } else {
            double bws = (m.t != GenSym.empty)? m.v.xval : 1.0;
            int ctrl = (r.t != GenSym.empty)? cast(int) r.v.xval : 0;

            A = get_density_matrix(l, bws, ctrl, p);
        }
    } else if (f == GenFunc.f_monthlen) {
        double * movec = null;
        double * yrvec = null;
        int wkdays, julian = 0;
        int mo = 0, yr = 0;
        int rettype = GenSym.num;

        post_process = 0;
        wkdays = node_get_int(r, p);
        if (!p.err && (wkdays < 5 || wkdays > 7)) {
            p.err = GretlError.invarg;
        }
        if (!p.err) {
            node_get_int_or_series(&mo, &movec, l, p);
            if (!p.err) {
                if (movec != null) {
                    rettype = GenSym.series;
                } else if (mo < 1 || mo > 12) {
                    p.err = GretlError.invarg;
                }
            }
        }
        if (!p.err) {
            node_get_int_or_series(&yr, &yrvec, m, p);
            if (!p.err) {
                if (yrvec != null) {
                    rettype = GenSym.series;
                } else if (yr < 0) {
                    yr = -yr;
                    julian = 1;
                }
            }
        }
        reset_p_aux(p, save_aux);
        if (!p.err && rettype == GenSym.num) {
            ret = aux_scalar_node(p);
            if (ret != null) {
                ret.v.xval = get_days_in_month(mo, yr, wkdays, julian);
            }
        } else if (!p.err) {
            ret = aux_series_node(p);
            if (ret != null) {
                p.err = fill_monthlen_array(ret.v.xvec,
                                            p.dset.t1, p.dset.t2,
                                            wkdays, mo, yr,
                                            movec, yrvec,
                                            julian);
            }
        }
    } else if (f == GenFunc.f_setnote || f == GenFunc.f_brename) {
        post_process = 0;
        if (l.t != GenSym.bundle) {
            node_type_error(f, 1, GenSym.bundle, l, p);
        } else if (m.t != GenSym.str) {
            node_type_error(f, 2, GenSym.str, m, p);
        } else if (r.t != GenSym.str) {
            node_type_error(f, 3, GenSym.str, r, p);
        } else {
            reset_p_aux(p, save_aux);
            ret = aux_scalar_node(p);
            if (!p.err && f == GenFunc.f_setnote) {
                ret.v.xval = gretl_bundle_set_note(l.v.b, m.v.str, r.v.str);
            } else if (!p.err) {
                p.err = gretl_bundle_rekey_data(l.v.b, m.v.str, r.v.str);
                if (!p.err) {
                    ret.v.xval = 0;
                }
            }
        }
    } else if (f == GenFunc.f_bwfilt) {
        GretlMatrixRaw * tmp = null;

        post_process = 0;
        if (l.t != GenSym.series) {
            if (l.t == GenSym.mat) {
                cast_to_series(l, f, &tmp, null, null, p);
            } else {
                node_type_error(f, 1, GenSym.series, l, p);
            }
        } else if (m.t != GenSym.num) {
            node_type_error(f, 2, GenSym.num, m, p);
        } else if (r.t != GenSym.num) {
            node_type_error(f, 3, GenSym.num, r, p);
        } else {
            ret = aux_series_node(p);
            if (!p.err) {
                p.err = butterworth_filter(l.v.xvec, ret.v.xvec, p.dset,
                                           cast(int) m.v.xval, r.v.xval);
            }
        }
        if (tmp != null) {
            l.v.m = tmp;
        }
    } else if (f == GenFunc.f_mlag) {
        GretlMatrixRaw * m1 = node_get_real_matrix(l, p, 0, 1);
        GretlMatrixRaw * m2 = node_get_real_matrix(m, p, 1, 2);

        if (p.err) {
            ; /* skip the rest */
        } else if (r.t != GenSym.num && r.t != GenSym.empty) {
            /* optional scalar */
            node_type_error(f, 3, GenSym.num, r, p);
        } else {
            double missval = (r.t == GenSym.num)? r.v.xval : 0.0;

            A = gretl_matrix_lag(m1, m2, Gretlopt.OPT_L, missval);
        }
    } else if (f == GenFunc.f_lrcovar) {
        GretlMatrixRaw * mc = node_get_real_matrix(l, p, 0, 1);
        int d = 1; /* demean the matrix arg? */

        if (!p.err) {
            d = node_get_bool(r, p, d);
        }
        if (!p.err) {
            A = long_run_covariance(mc, d, &p.err);
        }
    } else if (f == GenFunc.f_eigsolve) {
        GretlMatrixRaw * m1 = node_get_real_matrix(l, p, 0, 1);
        GretlMatrixRaw * m2 = node_get_real_matrix(m, p, 1, 2);

        if (p.err) {
            ; /* skip the rest */
        } else if (r.t != GenSym.empty && r.t != GenSym.u_addr) {
            /* optional matrix-pointer */
            node_type_error(f, 3, GenSym.u_addr, r, p);
        } else {
            GretlMatrixRaw * V = null;

            if (r.t == GenSym.u_addr) {
                V = ptr_node_get_matrix(r, p);
            }
            if (!p.err) {
                A = user_gensymm_eigenvals(m1, m2, V, &p.err);
            }
        }
    } else if (f == GenFunc.f_princomp) {
        if (l.t != GenSym.mat) {
            node_type_error(f, 1, GenSym.mat, l, p);
        } else if (m.t != GenSym.num) {
            node_type_error(f, 2, GenSym.num, m, p);
        } else if (r.t != GenSym.empty && r.t != GenSym.num) {
            /* optional boolean */
            node_type_error(f, 3, GenSym.num, r, p);
        } else {
            int cov = null_node(r) ? 0 : node_get_int(r, p);
            Gretlopt opt = cov ? Gretlopt.OPT_V : Gretlopt.OPT_NONE;
            int k = node_get_int(m, p);

            if (!p.err) {
                A = gretl_matrix_pca(l.v.m, k, opt, &p.err);
            }
        }
    } else if (f == GenFunc.f_halton) {
        if (l.t != GenSym.num) {
            node_type_error(f, 1, GenSym.num, l, p);
        } else if (m.t != GenSym.num) {
            node_type_error(f, 2, GenSym.num, m, p);
        } else if (r.t != GenSym.empty && r.t != GenSym.num) {
            /* optional offset */
            node_type_error(f, 3, GenSym.num, r, p);
        } else {
            int offset = null_node(r) ? 10 : node_get_int(r, p);
            int rows = node_get_int(l, p);
            int cols = node_get_int(m, p);

            if (!p.err) {
                A = halton_matrix(rows, cols, offset, &p.err);
            }
        }
    } else if (f == GenFunc.f_iwishart) {
        if (l.t != GenSym.mat && l.t != GenSym.num) {
            node_type_error(f, 1, GenSym.mat, l, p);
        } else if (!scalar_node(m)) {
            node_type_error(f, 2, GenSym.num, m, p);
        } else if (r.t != GenSym.empty && r.t != GenSym.num) {
            /* optional number of replications */
            node_type_error(f, 3, GenSym.num, r, p);
        } else {
            GretlMatrixRaw * S = node_get_real_matrix(l, p, 0, 0);
            int v = node_get_int(m, p);
            int N = null_node(r) ? 0 : node_get_int(r, p);

            if (!p.err) {
                if (N == 0) {
                    A = inverse_wishart_matrix(S, v, &p.err);
                } else {
                    A = inverse_wishart_sequence(S, v, N, &p.err);
                }
            }
        }
    } else if (f == GenFunc.f_aggrby) {
        DatasetRaw * dset = p.dset;
        DatasetRaw * mdset = null;
        const(char) * fnname = null;
        const(double) * x = null;
        const(double) * y = null;
        int * xlist = null;
        int * ylist = null;
        GretlMatrixRaw * xm = null;
        GretlMatrixRaw * ym = null;
        int matrix_version = 0;
        int data_args = 0;
        int mat_args = 0;

        /* if present, the optional r node must hold a string */
        if (!null_or_string(r)) {
            node_type_error(f, 3, GenSym.str, r, p);
        } else if (!null_node(r)) {
            fnname = r.v.str;
        }
        if (!p.err && !null_node(l)) {
            /* the [optional] l node holds @x */
            if (l.t == GenSym.series) {
                x = l.v.xvec; data_args++;
            } else if (l.t == GenSym.list) {
                xlist = l.v.ivec; data_args++;
            } else if (l.t == GenSym.mat) {
                xm = l.v.m; data_args++;
                mat_args++;
            } else {
                p.err = GretlError.types;
            }
        }
        if (!p.err) {
            /* the [required] m node holds @byvar */
            if (m.t == GenSym.series) {
                y = m.v.xvec; data_args++;
            } else if (m.t == GenSym.list) {
                ylist = m.v.ivec; data_args++;
            } else if (m.t == GenSym.mat) {
                ym = m.v.m; data_args++;
                mat_args++;
            } else {
                p.err = GretlError.types;
            }
        }
        if (!p.err) {
            if (mat_args > 0 && mat_args < data_args) {
                /* matrix and series/list arguments can't be mixed */
                p.err = GretlError.types;
            } else if (mat_args > 0) {
                static if (AGGR_PREFER_MATRIX) {
                    matrix_version = 1;
                } else {
                    /* convert to dataset (with borrowed Z) */
                    dset = mdset = matrix_dset_plus_lists(xm, ym, &xlist,
                                                          &ylist, &p.err);
                }
            }
        }
        if (!p.err && ylist != null) {
            p.err = aggregate_discrete_check(ylist, dset);
        }
        if (!p.err && matrix_version) {
            A = matrix_aggregate(xm, ym, fnname, &p.err);
        } else if (!p.err) {
            A = aggregate_by(x, y, xlist, ylist, fnname, dset, &p.err);
        }
        if (mdset != null) {
            destroy_dataset(mdset);
            free(xlist);
            free(ylist);
        }
    } else if (f == GenFunc.f_mweights) {
        if (!scalar_node(l)) {
            node_type_error(f, 1, GenSym.num, l, p);
        } else if (m.t != GenSym.mat) {
            node_type_error(f, 2, GenSym.mat, m, p);
        } else if (!scalar_node(r) && r.t != GenSym.str) {
            node_type_error(f, 3, GenSym.num, r, p);
        } else {
            int length = node_get_int(l, p);
            int method = node_get_midas_method(r, p);
            GretlMatrixRaw * wm = node_get_real_matrix(m, p, 1, 2);

            if (!p.err) {
                A = midas_weights(length, wm, method, &p.err);
            }
        }
    } else if (f == GenFunc.f_mgradient) {
        if (!scalar_node(l)) {
            node_type_error(f, 1, GenSym.num, l, p);
        } else if (m.t != GenSym.mat) {
            node_type_error(f, 2, GenSym.mat, m, p);
        } else if (!scalar_node(r) && r.t != GenSym.str) {
            node_type_error(f, 3, GenSym.num, r, p);
        } else {
            int length = node_get_int(l, p);
            int method = node_get_midas_method(r, p);
            GretlMatrixRaw * gm = node_get_real_matrix(m, p, 1, 2);

            if (!p.err) {
                A = midas_gradient(length, gm, method, &p.err);
            }
        }
    } else if (f == GenFunc.f_resample) {
        int blocklen = 0, draws = 0;

        if (l.t != GenSym.mat) {
            node_type_error(f, 1, GenSym.mat, l, p);
        }
        if (!p.err && !null_node(m)) {
            blocklen = node_get_int(m, p);
        }
        if (!p.err && !null_node(r)) {
            draws = node_get_int(r, p);
        }
        if (!p.err) {
            if (blocklen != 0) {
                A = gretl_matrix_block_resample(l.v.m, blocklen, draws, &p.err);
            } else {
                A = gretl_matrix_resample(l.v.m, draws, &p.err);
            }
        }
    } else if (f == GenFunc.hf_regls) {
        post_process = 0;
        if (null_node(l) && null_node(m) && null_node(r)) {
            /* doing automatic MPI: no args needed */
            alias RegFunc = extern (C) int function(Prn *);
            RegFunc regfunc = cast(RegFunc) get_plugin_function("regls_xv_mpi");

            if (regfunc == null) {
                p.err = GretlError.fopen;
            } else {
                p.err = regfunc(p.prn);
            }
        } else {
            alias RegFunc = extern (C) int function(const(GretlMatrixRaw) *, const(GretlMatrixRaw) *,
                                                    GretlBundle *, Prn *);
            RegFunc regfunc = cast(RegFunc) get_plugin_function("gretl_regls");

            if (regfunc == null) {
                p.err = GretlError.fopen;
            } else if (l.t != GenSym.mat || m.t != GenSym.mat) {
                p.err = GretlError.types;
            } else {
                if (r.t == GenSym.u_addr) {
                    r = ptr_node_get_referent_node(r, p);
                }
                if (r.t != GenSym.bundle) {
                    p.err = GretlError.types;
                }
            }
            if (!p.err) {
                p.err = regfunc(l.v.m, m.v.m, r.v.b, p.prn);
            }
        }
        if (!p.err) {
            ret = aux_scalar_node(p);
            ret.v.xval = 0;
        }
    } else if (f == GenFunc.f_stack) {
        int length = 0;
        int offset = 0;
        int * list = null;

        post_process = 0;
        ret = aux_empty_series_node(p);
        list = node_get_list(l, p);
        if (null_node(m)) {
            p.err = GretlError.args;
        } else {
            length = node_get_int(m, p);
        }
        if (!p.err && !null_node(r)) {
            offset = node_get_int(r, p);
        }
        if (!p.err) {
            p.err = build_stacked_series(&ret.v.xvec, list, length, offset,
                                         p.dset);
        }
        free(list);
    } else if (f == GenFunc.f_vma) {
        if (l.t != GenSym.mat) {
            /* matrix A, required */
            node_type_error(f, 1, GenSym.mat, l, p);
        } else if (m.t != GenSym.mat && m.t != GenSym.empty) {
            /* matrix C, optional */
            node_type_error(f, 2, GenSym.mat, m, p);
        } else if (r.t != GenSym.num && r.t != GenSym.empty) {
            /* horizon, optional */
            node_type_error(f, 3, GenSym.num, r, p);
        } else {
            GretlMatrixRaw * compan_top = node_get_real_matrix(l, p, 0, 1);
            int horizon = null_node(r) ? 24: node_get_int(r, p);
            int n = compan_top.rows;
            GretlMatrixRaw * C = null;

            if (m.t != GenSym.empty) {
                C = node_get_real_matrix(m, p, 1, 2);
                if (C.rows != n || C.cols !=n) {
                    p.err = GretlError.nonconf;
                }
            }
            if (!p.err) {
                A = vma_rep(compan_top, C, horizon, &p.err);
            }
        }
    } else if (f == GenFunc.f_bcheck) {
        GretlBundle * lb = null;
        GretlBundle * mb = null;
        GretlArray * reqd = null;

        post_process = 0;
        if ((lb = ptr_node_get_bundle(l, p)) == null) {
            node_type_error(f, 1, GenSym.bundle, l, p);
        } else if (!null_node(m) && m.t != GenSym.bundle) {
            node_type_error(f, 2, GenSym.bundle, m, p);
        } else if (!null_node(r) && r.t != GenSym.array) {
            node_type_error(f, 3, GenSym.array, r, p);
        } else {
            ret = aux_scalar_node(p);
        }
        if (!p.err) {
            if (!null_node(m)) {
                mb = m.v.b;
            }
            if (!null_node(r)) {
                reqd = r.v.a;
            }
        }
        if (!p.err) {
            ret.v.xval = gretl_bundle_extract_args(lb, mb, reqd,
                                                   null, p.prn,
                                                   &p.err);
        }
    } else if (f == GenFunc.f_sphcorr) {
        int mode = null_node(m) ? 0 : node_get_int(m, p);
        GretlMatrixRaw * J = null;

        if (!p.err) {
            if (mode < 0 || mode > 2) {
                p.err = GretlError.invarg;
            } else if (l.t != GenSym.mat) {
                node_type_error(f, 1, GenSym.mat, l, p);
            } else if (!null_node(r) && r.t != GenSym.u_addr) {
                node_type_error(f, 3, GenSym.u_addr, r, p);
            } else if (!null_node(r)) {
                if (null_node(m)) {
                    mode = 1;
                }
                J = ptr_node_get_matrix(r, p);
            }
        }
        if (!p.err) {
            if (mode == 0) {
                A = omega_from_R(l.v.m, &p.err);
            } else {
                A = R_from_omega(l.v.m, mode == 2, J, &p.err);
            }
        }
    } else if (f == GenFunc.f_sdc || f == GenFunc.f_cdemean || f == GenFunc.f_stdize) {
        GretlMatrixRaw * X = null;
        int param = f == GenFunc.f_stdize ? 1 : 0;
        int skip_na = 0;

        if (l.t != GenSym.mat) {
            p.err = GretlError.types;
        } else if (f == GenFunc.f_sdc) {
            X = l.v.m;
        } else {
            A = gretl_matrix_copy(l.v.m);
            if (A == null) {
                p.err = GretlError.alloc;
            }
        }
        if (!p.err && !null_node(m)) {
            param = node_get_int(m, p);
        }
        if (!p.err && !null_node(r)) {
            skip_na = node_get_bool(r, p, 0);
        }
        if (!p.err) {
            if (f == GenFunc.f_sdc) {
                A = gretl_matrix_column_sd(X, param, skip_na, &p.err);
            } else if (f == GenFunc.f_cdemean && param == 0) {
                p.err = gretl_matrix_center(A, skip_na);
            } else if (f == GenFunc.f_cdemean) {
                p.err = gretl_matrix_standardize(A, 1, skip_na);
            } else if (f == GenFunc.f_stdize && param == -1) {
                p.err = gretl_matrix_center(A, skip_na);
            } else {
                p.err = gretl_matrix_standardize(A, param, skip_na);
            }
        }
    } else if (f == GenFunc.f_tsolvepd) {
        if (l.t != GenSym.mat) {
            /* matrix a, required */
            node_type_error(f, 1, GenSym.mat, l, p);
        } else if (m.t != GenSym.mat) {
            /* matrix x, required */
            node_type_error(f, 2, GenSym.mat, m, p);
        } else {
            GretlMatrixRaw * a = mat_node_get_real_matrix(l, p);
            GretlMatrixRaw * x = mat_node_get_real_matrix(m, p);
            UserVar * uv = null;

            if (!p.err && !null_node(r)) {
                uv = ptr_node_get_uvar(r, GenSym.num, p);
            }
            if (!p.err) {
                double d;
                double * pdet = (uv != null)? &d : null;

                A = gretl_toeplitz_solve(a, null, x, 1, pdet, &p.err);
                if (!p.err && uv != null) {
                    user_var_set_scalar_value(uv, d);
                }
            }
        }
    } else if (f == GenFunc.f_kmeans) {
        post_process = 0;
        if (l.t != GenSym.mat && l.t != GenSym.list) {
            /* matrix or list, required */
            node_type_error(f, 1, 0, l, p);
        } else if (m.t != GenSym.num && m.t != GenSym.mat) {
            /* int k or matrix c0, required */
            node_type_error(f, 2, 0, m, p);
        } else {
            GretlMatrixRaw * a = null;
            GretlMatrixRaw * c0 = null;
            GretlBundle * opts = null;
            int * alist = null;
            int k = 0;

            if (l.t == GenSym.mat) {
                a = mat_node_get_real_matrix(l, p);
            } else {
                alist = node_get_list(l, p);
            }
            if (m.t == GenSym.num) {
                k = node_get_int(m, p);
            } else {
                c0 = m.v.m;
            }
            if (!p.err && !null_node(r)) {
                /* arg3: optional bundle */
                if (r.t == GenSym.bundle) {
                    opts = r.v.b;
                } else {
                    node_type_error(f, 3, GenSym.bundle, r, p);
                }
            }
            if (!p.err) {
                ret = aux_bundle_node(p);
            }
            if (!p.err) {
                ret.v.b = gretl_kmeans(a, alist, k, c0, opts,
                                       p.dset, p.prn, &p.err);
            }
            free(alist);
        }
    }

    if (post_process) {
        if (!p.err) {
            reset_p_aux(p, save_aux);
            ret = aux_matrix_node(p);
            if (!p.err) {
                ret.v.m = A;
            }
        }
        if (p.err) {
            /* don't leak memory on error */
            gretl_matrix_free(A);
        }
    }

    return ret;
}
