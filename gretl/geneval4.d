/* Port of the fourth part (lines 14334-18374) of lib/src/geneval.c: geneval4.c. See geneval0.d for the notes that
   apply to all five parts of the file. */

module gretl.geneval4;


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
import gretl.adf_kpss : kpss_critvals;
import gretl.calendar : epoch_day_from_julian_ymd, epoch_day_from_ymd, iso_basic_to_extended;
import gretl.dataio : check_varname, date_maj_min, ntolabel;
import gretl.dataset : dataset_add_NA_series, dataset_get_mapfile, dataset_get_structure, get_current_map, get_panel_group_names, rename_series, series_alphabetize_strings, series_destroy_string_table, series_get_string_vals, series_reorder_strings, series_set_string_vals;
import gretl.dataseth : annual_data, calendar_data, dataset_has_panel_time, dataset_is_panel, dataset_is_time_series, decennial_data, quarterly_or_monthly, sample_size;
import gretl.describe : get_last_result_data;
import gretl.geneval1 : aux_list_node, aux_series_node, aux_sized_matrix_node, e_types, free_node, gretl_type_from_gen_type, is_function_call, mat_node_get_real_matrix, n_args_error, no_data_error, node_get_bool, node_get_guint32, node_get_int, node_get_real_matrix, node_get_scalar, node_replace_scalar, ok_list_node, ok_list_node_plus, ptr_node_get_matrix, ptr_node_get_referent_node, ptr_node_get_uvar, reset_p_aux, tmp_matrix_from_series, typestr;
import gretl.geneval2 : list_from_strings_array, list_to_series_func;
import gretl.geneval3 : eval_ufunc, get_bundle_member, matrix_from_list, matrix_from_scalars, pretend_matrix_scalar_func, real_matrix_from_list, series_scalar_func;
import gretl.geneval5 : eval, series_to_matrix;
import gretl.genmain : current_series_index, generate;
import gretl.gmem : g_free;
import gretl.graphingh : gretlRGB;
import gretl.gretl_errors : get_gretl_errno, gretl_errmsg_set, gretl_errmsg_sprintf, gretl_warnmsg_sprintf;
import gretl.gretl_list : ellipsis_list, gretl_consecutive_list_new, gretl_is_midas_list, gretl_list_append_list, gretl_list_new, gretl_lists_join_with_separator, varname_match_list;
import gretl.gretl_matrix : gretl_matrix_alloc, gretl_matrix_commute, gretl_matrix_copy, gretl_matrix_covariogram, gretl_matrix_free, gretl_null_matrix_new, gretl_toeplitz_solve;
import gretl.gretl_model : gretl_model_get_data_element;
import gretl.gretl_prn : Prn, pprintf, pputc, pputs;
import gretl.gretl_utils : get_last_break, get_last_lnl, get_last_pvals_matrix, get_last_pvalue, get_last_test_matrix, get_last_test_statistic, get_last_test_type, gretl_stopwatch, gretl_version_number;
import gretl.gstrfuncs : g_strdup_printf;
import gretl.libset : get_matrix_mask, gretl_iterating, libset_get_bool, libset_get_double, libset_get_int;
import gretl.missingh : NADBL, na;
import gretl.pvalues : hypergeo;
import gretl.strutils : gretl_strdup, strings_array_sort;
import gretl.matrix_extra : gretl_matrix_data_subset, gretl_matrix_from_scalar, gretl_quadrule_matrix_new, gretl_vector_from_series, vector_from_strings;
import gretl.gretl_cmatrix : gretl_cmatrix_build;
import gretl.gretl_bfgs : user_NR;
import gretl.gretl_normal : gretl_GHK;
import gretl.varh : GretlVar;
import gretl.printscan : do_printf, do_sprintf_function, do_sscanf;
import gretl.var : gretl_FEVD_from_bundle, gretl_IRF_from_bundle, gretl_VAR_get_FEVD_matrix;
import gretl.compare : QLR_pval;
import gretl.graphing : gretl_png_font_string, numeric_color_from_string;
import gretl.gretl_array : arglist_validate, gretl_array_append_object, gretl_array_destroy, gretl_array_from_strings, gretl_array_get_length, gretl_array_get_strings, gretl_array_get_type, gretl_array_new, gretl_array_set_type;
import gretl.gretl_bundle : bundle_from_model, gretl_bundle_destroy, gretl_bundle_donate_data, gretl_bundle_donate_series, gretl_bundle_get_array, gretl_bundle_get_keys, gretl_bundle_get_member_type, gretl_bundle_get_private_data, gretl_bundle_get_series, gretl_bundle_get_type, gretl_bundle_has_key, gretl_bundle_new, gretl_bundle_set_data, gretl_bundle_set_series, kalman_bundle_new, type_can_be_bundled;
import gretl.uservar : get_user_var_by_name, get_user_var_of_type_by_name, user_var_get_type_by_name, user_var_replace_value, user_var_set_scalar_value;
import gretl.gretl_typemap : gretl_type_get_name;
import gretl.genfuncs : bkbp_filter, exponential_movavg_series, felogit_rec_loglik, fill_dataset_dates_series, filter_matrix, filter_series, gretl_bessel, gretl_loess, list_from_matrix, make_panel_time_var, midas_multipliers, movavg_series, nadaraya_watson, substitute_values, two_matrices_minmax, two_series_minmax;
import gretl.gretl_color : colormix_array;
import gretl.usermat : user_matrix_GHK, user_matrix_ols, user_matrix_rls;
import gretl.random : gretl_rand_get_seed, gretl_rand_int_minmax;

/* unported dependencies: per-file declarations from the 2026c headers */
extern (C) {
    Node * bncopy(Node * t, int * err);
    int function_lookup(const(char) * s);
    int geoplot_driver(const(char) * fname, GretlBundle * map, int plv, const(double) * plx, const(DatasetRaw) * dset, GretlBundle * opts);
    void * get_genr_function_pointer(int f);
    void * get_last_model(GretlObjType * type);
    int get_loop_renaming();
    GretlMatrixRaw * get_tdisagg_matrix(TdisaggInfo * tdi, DatasetRaw * dset, GretlBundle * b, GretlBundle * r, Prn * prn, int * err);
    Ufunc * get_user_function_by_name(const(char) * name);
    const(char) * getsymb(int t);
    int gretl_function_depth();
    int kalman_bundle_filter(GretlBundle * b, Prn * prn, int * errp);
    GretlMatrixRaw * kalman_bundle_simdata(GretlBundle * b, const(GretlMatrixRaw) * U, Prn * prn, int * err);
    GretlMatrixRaw * kalman_bundle_simulate(GretlBundle * b, const(GretlMatrixRaw) * U, int get_state, Prn * prn, int * err);
    int kalman_bundle_smooth(GretlBundle * b, int dist, Prn * prn);
    GretlMatrixRaw * last_model_get_boot_ci(int cnum, const(DatasetRaw) * dset, int B, double alpha, int method, int studentize, int * err);
    double last_model_get_boot_pval(int cnum, const(DatasetRaw) * dset, int B, int method, int * err);
    GretlMatrixRaw * last_model_get_irf_matrix(int targ, int shock, double alpha, const(DatasetRaw) * dset, int * err);
    GretlMatrixRaw * matrix_chowlin(const(GretlMatrixRaw) * Y, const(GretlMatrixRaw) * X, int s, int * err);
    int mvar_lookup(const(char) * s);
    Node * new_node(int t);
    int parser_ensure_error_buffer(Parser * p);
    void sample_range_get_extrema(const(DatasetRaw) * dset, int * t1, int * t2);
    GretlMatrixRaw * saved_object_build_matrix(const(char) * oname, int idx, const(DatasetRaw) * dset, int * err);
    void * saved_object_get_array(const(char) * oname, int idx, const(DatasetRaw) * dset, int * err);
    GretlType saved_object_get_data_type(const(char) * name, ModelDataIndex idx);
    int * saved_object_get_list(const(char) * oname, int idx, int * err);
    GretlMatrixRaw * saved_object_get_matrix(const(char) * oname, int idx, int * err);
    double saved_object_get_scalar(const(char) * oname, int idx, DatasetRaw * dset, int * err);
    int saved_object_get_series(double * x, const(char) * oname, int idx, const(DatasetRaw) * dset);
    char * saved_object_get_string(const(char) * oname, int idx, const(DatasetRaw) * dset, int * err);
    const(char) * mvarname(int t);
}

Node * geoplot_node(Node * l, Node * m, Node * r, Parser * p)
{
    Node * ret = aux_scalar_node(p);

    if (!p.err) {
        const(char) * mapfile = null;
        GretlBundle * mapbun = null;
        int plv = -1;         /* payload series ID */
        double * plx = null;  /* payload data */
        GretlBundle * opts = null;

        if (null_node(l) || l.t == GenSym.str || l.t == GenSym.bundle) {
            /* [map-fname-or-bundle] [,series] [,options] */
            if (!null_node(l)) {
                if (l.t == GenSym.str) {
                    mapfile = l.v.str;
                } else {
                    mapbun = l.v.b;
                }
            }
            if (!null_node(m)) {
                if (m.t == GenSym.series) {
                    plv = m.vnum;
                    plx = m.v.xvec;
                } else if (m.t == GenSym.bundle) {
                    opts = m.v.b;
                } else {
                    p.err = GretlError.invarg;
                }
            }
            if (!p.err && !null_node(r)) {
                if (opts == null) {
                    opts = r.v.b;
                } else {
                    p.err = GretlError.invarg;
                }
            }
        } else if (l.t == GenSym.series) {
            /* series [,options] (map is implicit) */
            plv = l.vnum;
            plx = l.v.xvec;
            if (!null_node(m)) {
                if (m.t == GenSym.bundle) {
                    opts = m.v.b;
                } else {
                    p.err = GretlError.invarg;
                }
            }
            if (!p.err && !null_node(r)) {
                p.err = GretlError.invarg;
            }
        } else if (null_node(l) && null_node(m) && null_node(r)) {
            ; /* implicit map, no payload, no options */
        } else {
            p.err = GretlError.invarg;
        }

        if (!p.err) {
            int gerr = geoplot_driver(mapfile, mapbun, plv, plx,
                                      p.dset, opts);

            ret.v.xval = gerr;
            p.err = gerr;
        }
    }

    return ret;
}

int scan_to_vector(Node * n, const(char) * fmt,
                   const(char) * arg, int * err)
{
    GretlMatrixRaw * m = null;
    UserVar * uvar = null;
    char ** S = null;
    int ns = 0;
    int nmax = 0;

    uvar = get_user_var_of_type_by_name(arg, GretlType.matrix);
    if (uvar == null) {
        *err = GretlError.invarg;
        return 0;
    }

    if (gretl_array_get_type(n.v.a) != GretlType.strings) {
        *err = e_types(n);
    } else {
        S = gretl_array_get_strings(n.v.a, &ns);
        m = vector_from_strings(S, ns, fmt, &nmax, err);
    }

    if (!*err) {
        *err = user_var_replace_value(uvar, m, GretlType.matrix);
    }

    return nmax;
}

/* Code to handle printf() and sprintf(). The @l node is a
   format string; the @r node holds varargs or may be NULL.
*/

Node * print_func_node(Node * l, Node * r, int f, Parser * p)
{
    Node * ret;

    if (f == GenFunc.f_sprintf) {
        ret = aux_string_node(p);
    } else {
        ret = aux_scalar_node(p);
    }

    if (ret != null) {
        const(char) * fmt = l.v.str;
        const(char) * args = null;
        int n = 0;

        if (!null_node(r)) {
            args = r.v.str;
        }
        if (f == GenFunc.f_sprintf) {
            ret.v.str = do_sprintf_function(fmt, args, p.dset, &p.err);
        } else {
            p.err = do_printf(fmt, args, p.dset, p.prn, &n);
        }
        if (!p.err && f == GenFunc.f_printf) {
            ret.v.xval = n;
        }
    }

    return ret;
}

/* Code to handle sscanf(). The @l node must be a string, or array of
   strings; @m holds a format string; and @r holds varargs.
*/

Node * sscanf_node(Node * l, Node * m, Node * r, Parser * p)
{
    Node * ret = aux_scalar_node(p);

    if (ret != null) {
        const(char) * fmt = m.v.str;
        const(char) * args = r.v.str;
        const(char) * src = null;
        int n = 0;

        if (l.t == GenSym.array) {
            /* scanning array of strings to vector */
            n = scan_to_vector(l, fmt, r.v.str, &p.err);
        } else {
            src = l.v.str;
            p.err = do_sscanf(src, fmt, args, p.dset, &n);
        }
        if (!p.err) {
            ret.v.xval = n;
        }
    }

    return ret;
}

int x_to_period(double x, char c, int * julian, int * err)
{
    if (julian != null && c == 'y') {
        if (x < 0) {
            *julian = 1;
            x = -x;
        } else {
            *julian = 0;
        }
    }

    if (na(x)) {
        /* note: error not flagged here */
        return -1;
    } else if (x < 0 || fabs(x) > INT_MAX) {
        *err = GretlError.invarg;
        return -1;
    } else {
        int k = cast(int) x;
        int ret = cast(int) x;

        if (c == 'y' && k <= 0) {
            ret = -1;
        } else if (c == 'm' && (k < 1 || k > 12)) {
            ret = -1;
        } else if (c == 'd' && (k < 1 || k > 31)) {
            ret = -1;
        }

        if (ret <= 0) {
            fprintf(stderr, "epochday: got %c = %d!\n".ptr, c, k);
            *err = GretlError.invarg;
        }

        return ret;
    }
}

void fill_xymd(double * targ, double x)
{
    int rem;

    targ[0] = floor(x / 10000);
    rem = cast(int) (x - 10000 * targ[0]);
    targ[1] = floor(rem / 100);
    targ[2] = rem - 100 * targ[1];
}

void bad_date_message(int y, int m, int d)
{
    gretl_warnmsg_sprintf("%04d-%02d-%02d: %s".ptr, y, m, d,
                          _("non-existent date").ptr);
}

/* epochday policy: NAs for year, month or day give an NA result;
   non-NA but inherently out-of-bounds values for y, m or d (for
   example, negative y, m > 12, d > 31) will produce an error,
   and otherwise non-existent dates produce NA.
*/

Node * eval_epochday(Node * ny, Node * nm, Node * nd, Parser * p)
{
    Node * ret = null;
    double *[3] x = null;
    double[3] xymd;
    int[3] ymd = -1;
    const(char) * code = "ymd";
    int basic_input = 0;
    int n_series = 0;
    int julian = 0;
    double sval;
    int i;

    if (null_node(nm) && null_node(nd)) {
        /* try for ISO 8601 basic input */
        basic_input = 1;
        if (scalar_node(ny)) {
            sval = node_get_scalar(ny, p);
            if (!p.err) {
                fill_xymd(xymd.ptr, sval);
                for (i=0; i<3 && !p.err; i++) {
                    ymd[i] = x_to_period(xymd[i], code[i], null, &p.err);
                }
            }
        } else if (ny.t == GenSym.series) {
            x[0] = ny.v.xvec;
            n_series = 1;
        } else {
            node_type_error(GenFunc.f_epochday, 1, GenSym.num, ny, p);
        }
    } else {
        Node *[3] nodes = [ny, nm, nd];

        for (i=0; i<3 && !p.err; i++) {
            if (scalar_node(nodes[i])) {
                sval = node_get_scalar(nodes[i], p);
                if (!p.err) {
                    ymd[i] = x_to_period(sval, code[i], &julian, &p.err);
                }
            } else if (nodes[i].t == GenSym.series) {
                x[i] = nodes[i].v.xvec;
                n_series++;
            } else {
                node_type_error(GenFunc.f_epochday, i+1, GenSym.num, nodes[i], p);
            }
        }
    }

    if (!p.err) {
        double edt;
        int t, t1, t2;
        int y = ymd[0];
        int m = ymd[1];
        int d = ymd[2];

        if (n_series > 0) {
            t1 = p.dset.t1;
            t2 = p.dset.t2;
            ret = aux_series_node(p);
        } else {
            t1 = t2 = 0;
            ret = aux_scalar_node(p);
        }

        for (t=t1; t<=t2 && !p.err; t++) {
            if (basic_input) {
                if (x[0] != null) {
                    fill_xymd(xymd.ptr, x[0][t]);
                    y = x_to_period(xymd[0], 'y', null, &p.err);
                    m = x_to_period(xymd[1], 'm', null, &p.err);
                    d = x_to_period(xymd[2], 'd', null, &p.err);
                }
            } else {
                y = (x[0] == null)? y : x_to_period(x[0][t], 'y', &julian, &p.err);
                m = (x[1] == null)? m : x_to_period(x[1][t], 'm', null, &p.err);
                d = (x[2] == null)? d : x_to_period(x[2][t], 'd', null, &p.err);
            }
            if (p.err) {
                break;
            } else if (y < 0 || m < 0 || d < 0) {
                /* got an NA somewhere */
                edt = NADBL;
            } else {
                if (julian) {
                    edt = epoch_day_from_julian_ymd(y, m, d);
                } else {
                    edt = epoch_day_from_ymd(y, m, d);
                }
                if (edt <= 0) {
                    bad_date_message(y, m, d);
                    edt = NADBL;
                }
            }
            if (n_series > 0) {
                ret.v.xvec[t] = edt;
            } else {
                ret.v.xval = edt;
            }
        }
    }

    return ret;
}

/* Bessel function handler: the @r node can be of scalar, series or
   matrix type. But the @m node (order) must hold a scalar.
*/

Node * eval_bessel_func(Node * l, Node * m, Node * r, Parser * p)
{
    char ftype;
    double v;
    Node * ret = null;

    if (!starting(p) && r.t != GenSym.series) {
        return aux_any_node(p);
    }

    ftype = l.v.str[0];
    v = node_get_scalar(m, p);

    if (r.t == GenSym.num) {
        double x = r.v.xval;

        ret = aux_scalar_node(p);
        if (ret != null) {
            ret.v.xval = gretl_bessel(ftype, v, x, &p.err);
        }
    } else if (r.t == GenSym.mat) {
        const(GretlMatrixRaw) * x = r.v.m;
        int i, n = x.rows * x.cols;

        ret = aux_sized_matrix_node(p, x.rows, x.cols, 0);
        if (ret != null) {
            for (i=0; i<n && !p.err; i++) {
                ret.v.m.val[i] = gretl_bessel(ftype, v, x.val[i], &p.err);
            }
        }
    } else if (r.t == GenSym.series) {
        const(double) * x = r.v.xvec;
        int t1 = autoreg(p) ? p.obs : p.dset.t1;
        int t2 = autoreg(p) ? p.obs : p.dset.t2;
        int t;

        ret = aux_series_node(p);
        if (ret != null) {
            for (t=t1; t<=t2 && !p.err; t++) {
                ret.v.xvec[t] = gretl_bessel(ftype, v, x[t], &p.err);
            }
        }
    }

    return ret;
}

/* replace_value: non-interactive search-and-replace for series and
   matrices.  @src holds the series or matrix of which we want a
   modified copy; @n0 holds the value (or vector of values) to be
   replaced; and @n1 holds the replacement value(s). It would be nice
   to extend this to lists.
*/

Node * replace_value(Node * src, Node * n0, Node * n1, Parser * p)
{
    GretlVector * vx0 = null;
    GretlVector * vx1 = null;
    double x0 = 0, x1 = 0;
    int k0 = -1, k1 = -1;
    Node * ret = null;

    if (!starting(p)) {
        return aux_any_node(p);
    }

    /* n0: the original value, to be replaced */
    if (n0.t == GenSym.num) {
        x0 = n0.v.xval;
    } else if (n0.t == GenSym.mat) {
        vx0 = n0.v.m;
        if (gretl_is_null_matrix(vx0)) {
            p.err = GretlError.data;
        } else if ((k0 = gretl_vector_get_length(vx0)) == 0) {
            p.err = GretlError.nonconf;
        }
    } else {
        node_type_error(GenFunc.f_replace, 1, GenSym.num, n0, p);
    }

    if (p.err) {
        return null;
    }

    /* n1: the replacement value */
    if (n1.t == GenSym.num) {
        x1 = n1.v.xval;
    } else if (n1.t == GenSym.mat) {
        vx1 = n1.v.m;
        if (gretl_is_null_matrix(vx1)) {
            p.err = GretlError.data;
        } else if ((k1 = gretl_vector_get_length(vx1)) == 0) {
            p.err = GretlError.nonconf;
        }
    } else {
        node_type_error(GenFunc.f_replace, 2, GenSym.num, n1, p);
    }

    if (!p.err) {
        if (n0.t == GenSym.num && n1.t == GenSym.mat) {
            /* can't replace scalar with vector */
            p.err = GretlError.types;
        } else if (k0 > 0 && k1 > 0 && k0 != k1) {
            /* if they're both vectors, they must be
               the same length */
            p.err = GretlError.nonconf;
        }
    }

    if (!p.err) {
        if (src.t == GenSym.series) {
            ret = aux_series_node(p);
        } else if (src.t == GenSym.mat) {
            ret = aux_matrix_node(p);
        } else {
            node_type_error(GenFunc.f_replace, 3, GenSym.series, src, p);
        }
    }

    if (!p.err) {
        double * px0 = (vx0 != null)? vx0.val : &x0;
        double * px1 = (vx1 != null)? vx1.val : &x1;
        GretlMatrixRaw * m = null;
        const(double) * x = null;
        double * targ = null;
        int n = 0;

        if (k0 < 0) k0 = 1;
        if (k1 < 0) k1 = 1;

        if (src.t == GenSym.series) {
            n = sample_size(p.dset);
            x = src.v.xvec + p.dset.t1;    /* source array */
            targ = ret.v.xvec + p.dset.t1; /* target array */
        } else if (src.t == GenSym.mat) {
            m = src.v.m;
            ret.v.m = gretl_matrix_copy(m);
            if (ret.v.m == null) {
                p.err = GretlError.alloc;
            } else {
                n = m.rows * m.cols;
                x = m.val;           /* source array */
                targ = ret.v.m.val; /* target array */
            }
        }

        if (!p.err) {
            substitute_values(targ, x, n, px0, k0, px1, k1);
        }
    }

    return ret;
}

Node * isoconv_node(Node * t, Node * n, Parser * p)
{
    Node * e;
    Node * ret = null;
    const(double) * x = null;
    double *[3] ymd = null;
    int i, k = n.v.bn.n_nodes;

    if (p.dset == null) {
        p.err = GretlError.nodata;
    } else if (k < 3 || k > 4) {
        n_args_error(k, 3, 4, t.t, p);
    }

    if (p.err) {
        return null;
    }

    for (i=0; i<k && !p.err; i++) {
        e = n.v.bn.n[i];
        if (i == 0) {
            if (e.t == GenSym.series) {
                x = e.v.xvec + p.dset.t1;
            } else {
                node_type_error(t.t, 1, GenSym.series, e, p);
            }
        } else if (i == 3 && null_node(e)) {
            ; /* OK for the last arg to be omitted */
        } else if (e.t != GenSym.u_addr) {
            node_type_error(t.t, i+1, GenSym.u_addr, e, p);
        } else {
            e = ptr_node_get_referent_node(e, p);
            if (e == null || e.t != GenSym.series) {
                node_type_error(t.t, i+1, GenSym.series, e, p);
            } else {
                ymd[i-1] = p.dset.Z[e.vnum] + p.dset.t1;
            }
        }
    }

    if (!p.err) {
        ret = aux_scalar_node(p);
    }

    if (!p.err) {
        int nobs = sample_size(p.dset);

        ret.v.xval = iso_basic_to_extended(x, ymd[0], ymd[1], ymd[2], nobs);
    }

    return ret;
}

/* The arguments here are:

   @A: array to which an element is to be added
   @n: node holding candidate array element
*/

int check_array_element_type(GretlArray * A, Node * n)
{
    GretlType t = gretl_array_get_type(A);
    int ok = 0;

    if (t == GretlType.any) {
        /* The array type is not yet determinate; this may be
           the case when when we're looking at the first element.
           If the type n->t is acceptable we use it to set the
           type of @A.
        */
        t = GretlType.none;
        if (n.t == GenSym.mat || n.t == GenSym.num) {
            t = GretlType.matrices;
        } else if (n.t == GenSym.str) {
            t = GretlType.strings;
        } else if (n.t == GenSym.bundle) {
            t = GretlType.bundles;
        } else if (n.t == GenSym.list) {
            t = GretlType.lists;
        } else if (n.t == GenSym.array) {
            t = GretlType.arrays;
        }
        if (t > GretlType.none) {
            gretl_array_set_type(A, t);
            ok = 1;
        }
    } else {
        /* We're looking for a match between the array type
           and the type of the candidate element.
        */
        if (t == GretlType.matrices) {
            ok = (n.t == GenSym.mat || n.t == GenSym.num);
        } else if (t == GretlType.strings) {
            ok = n.t == GenSym.str;
        } else if (t == GretlType.bundles) {
            ok = n.t == GenSym.bundle;
        } else if (t == GretlType.lists) {
            ok = n.t == GenSym.list;
        } else if (t == GretlType.arrays) {
            ok = n.t == GenSym.array;
        }
    }

    return ok ? 0 : GretlError.types;
}

void node_nullify_ptr(Node * n)
{
    if      (n.t == GenSym.mat)    n.v.m = null;
    else if (n.t == GenSym.str)    n.v.str = null;
    else if (n.t == GenSym.bundle) n.v.b = null;
    else if (n.t == GenSym.list)   n.v.ivec = null;
    else if (n.t == GenSym.array)  n.v.a = null;
    else if (n.t == GenSym.series) n.v.xvec = null;
}

/* supports retrieval of data for candidate array elements
   or bundle members
*/

void * node_get_ptr(Node * n, int f, Parser * p, int * donate)
{
    void * ptr = null;
    int t = n.t;

    /* default to copying the node's data */
    *donate = 0;

    if (f == GenFunc.f_defbundle || f == GenFunc.f_defargs) {
        /* specific to bundles */
        if (t == GenSym.array) {
            ptr = n.v.a;
        } else if (t == GenSym.series) {
            ptr = n.v.xvec;
        } else if (t == GenSym.num) {
            ptr = &n.v.xval;
        }
    }

    if (ptr == null) {
        /* common to arrays and bundles */
        if (t == GenSym.mat) {
            ptr = n.v.m;
        } else if (t == GenSym.str) {
            ptr = n.v.str;
        } else if (t == GenSym.bundle) {
            ptr = n.v.b;
        } else if (t == GenSym.list) {
            ptr = n.v.ivec;
        } else if (t == GenSym.array) {
            ptr = n.v.a;
        }
    }

    if (t == GenSym.num) {
        *donate = 1;
    } else if (!reusable(p) && is_tmp_node(n)) {
        *donate = 1;
        node_nullify_ptr(n);
    }

    return ptr;
}

/* given an FARGS node, detect if the first argument
   is a pointer to bundle */

int bundle_pointer_arg0(Node * t)
{
    Node * n = t.L;

    if (n.v.bn.n_nodes > 0) {
        Node * n0 = n.v.bn.n[0];

        if (n0.t == GenSym.u_addr && ubundle_node(n0.L)) {
            return 1;
        }
    }

    return 0;
}

struct FnArgspec {
    int f;
    int argc_min;
    int argc_max;
}

/* Check the number of arguments, @k, supplied for the multi-arg
   function @f.
*/

int check_argc(int f, int k, Parser * p)
{
    FnArgspec[19] fn_argspecs = [
        FnArgspec(GenFunc.f_bkfilt,    1, 4),
        FnArgspec(GenFunc.f_filter,    1, 5),
        FnArgspec(GenFunc.f_mcovg,     4, 4),
        FnArgspec(GenFunc.f_mrls,      4, 6),
        FnArgspec(GenFunc.f_fols,      3, 6),
        FnArgspec(GenFunc.f_loess,     2, 6),
        FnArgspec(GenFunc.f_ghk,       4, 5),
        FnArgspec(GenFunc.f_quadtab,   1, 4),
        FnArgspec(GenFunc.f_qlrpval,   4, 4),
        FnArgspec(GenFunc.f_bootci,    1, 5),
        FnArgspec(GenFunc.f_bootpval,  1, 3),
        FnArgspec(GenFunc.f_nadarwat,  2, 5),
        FnArgspec(GenFunc.f_hyp2f1,    4, 4),
        FnArgspec(GenFunc.f_chowlin,   2, 3),
        FnArgspec(GenFunc.f_midasmult, 1, 3),
        FnArgspec(GenFunc.f_tdisagg,   3, 5),
        FnArgspec(GenFunc.f_commute,   2, 5),
        FnArgspec(GenFunc.f_toepsolv,  3, 4),
        FnArgspec(GenFunc.f_rgbmix,    3, 4)
    ];
    int argc_min = 2;
    int argc_max = 4;
    int i;

    for (i=0; i<fn_argspecs.length; i++) {
        if (fn_argspecs[i].f == f) {
            argc_min = fn_argspecs[i].argc_min;
            argc_max = fn_argspecs[i].argc_max;
            break;
        }
    }

    if (k < argc_min || k > argc_max) {
        n_args_error(k, argc_min, argc_max, f, p);
    }

    return p.err;
}

int no_ts_data(Parser * p)
{
    if (!dataset_is_time_series(p.dset)) {
        gretl_errmsg_set(_("This function requires time-series data").ptr);
        p.err = GretlError.data;
    }

    return p.err;
}

bool nargs_needs_ts(int f)
{
    return f == GenFunc.f_bkfilt;
}

/* evaluate a built-in function that has more than three arguments */

Node * eval_nargs_func(Node * t, Node * n, Parser * p)
{
    Node * e;
    Node * ret = null;
    int i, k = n.v.bn.n_nodes;

    check_argc(t.t, k, p);
    if (p.err) {
        return null;
    }

    if (nargs_needs_ts(t.t) && no_ts_data(p)) {
        return null;
    }

    if (t.t == GenFunc.f_bkfilt) {
        const(double) * x = null;
        int[3] bk = 0;

        for (i=0; i<k && !p.err; i++) {
            e = n.v.bn.n[i];
            if (i == 0) {
                if (e.t != GenSym.series) {
                    node_type_error(t.t, 1, GenSym.series, e, p);
                } else {
                    x = e.v.xvec;
                }
            } else if (null_node(e)) {
                ; /* NULL arguments are OK */
            } else {
                bk[i-1] = node_get_int(e, p);
            }
        }
        if (!p.err) {
            ret = aux_series_node(p);
        }
        if (!p.err) {
            p.err = bkbp_filter(x, ret.v.xvec, p.dset, bk[0], bk[1], bk[2]);
        }
    } else if (t.t == GenFunc.f_filter) {
        const(double) * x = null; /* series */
        GretlMatrixRaw * X = null;
        GretlMatrixRaw * C = null;
        GretlMatrixRaw * A = null;
        GretlMatrixRaw * x0 = null;
        double y0 = 0;

        for (i=0; i<k && !p.err; i++) {
            e = n.v.bn.n[i];
            if (i == 0) {
                /* the series or matrix to filter */
                if (e.t != GenSym.series && e.t != GenSym.mat) {
                    node_type_error(t.t, i+1, 0, e, p);
                } else if (e.t == GenSym.series && no_ts_data(p)) {
                    ; /* p->err is set */
                } else if (e.t == GenSym.series) {
                    x = e.v.xvec;
                } else {
                    X = e.v.m;
                }
            } else if (i == 1) {
                /* matrix for MA polynomial (but we'll take a scalar) */
                if (e.t != GenSym.mat && e.t != GenSym.num && e.t != GenSym.empty) {
                    node_type_error(t.t, i+1, GenSym.mat, e, p);
                } else if (e.t != GenSym.empty) {
                    C = node_get_real_matrix(e, p, 0, 2);
                }
            } else if (i == 2) {
                /* matrix for AR polynomial (but we'll take a scalar) */
                if (e.t != GenSym.mat && e.t != GenSym.num && e.t != GenSym.empty) {
                    node_type_error(t.t, i+1, GenSym.mat, e, p);
                } else if (e.t != GenSym.empty) {
                    A = node_get_real_matrix(e, p, 1, 3);
                }
            } else if (i == 3) {
                /* initial (optional scalar) value for output series */
                if (e.t != GenSym.empty && !scalar_node(e)) {
                    node_type_error(t.t, i+1, GenSym.num, e, p);
                } else if (e.t != GenSym.empty) {
                    y0 = node_get_scalar(e, p);
                    if (!p.err && na(y0)) {
                        p.err = GretlError.missdata;
                    }
                }
            } else if (i == 4) {
                /* optional pre-sample x vector */
                if (e.t != GenSym.mat && e.t != GenSym.num && e.t != GenSym.empty) {
                    node_type_error(t.t, i+1, GenSym.mat, e, p);
                } else if (e.t != GenSym.empty) {
                    x0 = node_get_real_matrix(e, p, 1, 5);
                }
            }
        }
        if (!p.err) {
            if (X != null) {
                /* matrix output wanted */
                ret = aux_matrix_node(p);
                if (!p.err) {
                    ret.v.m = filter_matrix(X, A, C, y0, x0, &p.err);
                }
            } else if (x != null) {
                /* series output */
                ret = aux_series_node(p);
                if (!p.err) {
                    p.err = filter_series(x, ret.v.xvec, p.dset,
                                          A, C, y0, x0);
                }
            }
        }
    } else if (t.t == GenFunc.f_mcovg) {
        GretlMatrixRaw * X = null;
        GretlVector * u = null;
        GretlVector * w = null;
        int targ, maxlag = 0;

        for (i=0; i<k && !p.err; i++) {
            targ = (i == 3)? GenSym.num : GenSym.mat;
            e = n.v.bn.n[i];
            if ((i == 1 || i == 2) && null_node(e)) {
                ; /* for u or w, NULL is acceptable */
            } else if (e.t != targ) {
                node_type_error(t.t, i+1, targ, e, p);
            } else if (i == 0) {
                X = mat_node_get_real_matrix(e, p);
            } else if (i == 1) {
                u = mat_node_get_real_matrix(e, p);
            } else if (i == 2) {
                w = mat_node_get_real_matrix(e, p);
            } else if (i == 3) {
                maxlag = cast(int) e.v.xval;
            }
        }
        if (!p.err) {
            ret = aux_matrix_node(p);
        }
        if (!p.err) {
            ret.v.m = gretl_matrix_covariogram(X, u, w, maxlag, &p.err);
        }
    } else if (t.t == GenFunc.f_mols || t.t == GenFunc.f_mpols) {
        Gretlopt opt = (t.t == GenFunc.f_mpols)? Gretlopt.OPT_M : Gretlopt.OPT_NONE;
        GretlMatrixRaw *[2] M = null;
        GretlMatrixRaw * U = null;
        GretlMatrixRaw * V = null;
        char[2] freemat = 0;

        for (i=0; i<k && !p.err; i++) {
            e = n.v.bn.n[i];
            if (i < 2) {
                if (e.t == GenSym.series) {
                    M[i] = tmp_matrix_from_series(e, p);
                    freemat[i] = 1;
                } else {
                    M[i] = node_get_real_matrix(e, p, i, i+1);
                }
            } else {
                if (null_node(e)) {
                    ; /* OK */
                } else if (e.t != GenSym.u_addr) {
                    node_type_error(t.t, i+1, GenSym.u_addr, e, p);
                } else if (i == 2) {
                    U = ptr_node_get_matrix(e, p);
                } else {
                    V = ptr_node_get_matrix(e, p);
                }
            }
        }
        if (!p.err) {
            ret = aux_matrix_node(p);
        }
        if (!p.err) {
            ret.v.m = user_matrix_ols(M[0], M[1], null, U, V, null,
                                      opt, &p.err);
        }
        if (freemat[0]) gretl_matrix_free(M[0]);
        if (freemat[1]) gretl_matrix_free(M[1]);
    } else if (t.t == GenFunc.f_mrls) {
        GretlMatrixRaw *[4] M = null;
        GretlMatrixRaw * U = null;
        GretlMatrixRaw * V = null;

        for (i=0; i<k && !p.err; i++) {
            e = n.v.bn.n[i];
            if (i < 4) {
                M[i] = node_get_real_matrix(e, p, i, i+1);
            } else {
                if (null_node(e)) {
                    ; /* OK */
                } else if (e.t != GenSym.u_addr) {
                    node_type_error(t.t, i+1, GenSym.u_addr, e, p);
                } else if (i == 4) {
                    U = ptr_node_get_matrix(e, p);
                } else {
                    V = ptr_node_get_matrix(e, p);
                }
            }
        }
        if (!p.err) {
            ret = aux_matrix_node(p);
        }
        if (!p.err) {
            ret.v.m = user_matrix_rls(M[0], M[1], M[2], M[3], U, V, &p.err);
        }
    } else if (t.t == GenFunc.f_fols) {
        GretlMatrixRaw *[3] M = null;
        GretlMatrixRaw * U = null;
        GretlMatrixRaw * V = null;
        GretlMatrixRaw * A = null;

        for (i=0; i<k && !p.err; i++) {
            e = n.v.bn.n[i];
            if (i < 3) {
                M[i] = node_get_real_matrix(e, p, i, i+1);
            } else {
                if (null_node(e)) {
                    ; /* OK */
                } else if (e.t != GenSym.u_addr) {
                    node_type_error(t.t, i+1, GenSym.u_addr, e, p);
                } else if (i == 3) {
                    U = ptr_node_get_matrix(e, p);
                } else if (i == 4) {
                    V = ptr_node_get_matrix(e, p);
                } else {
                    A = ptr_node_get_matrix(e, p);
                }
            }
        }
        if (!p.err) {
            ret = aux_matrix_node(p);
        }
        if (!p.err) {
            ret.v.m = user_matrix_ols(M[0], M[1], M[2], U, V, A,
                                      Gretlopt.OPT_F, &p.err);
        }
    } else if (t.t == GenFunc.f_nrmax) {
        GretlMatrixRaw * b = null;
        const(char) * sf = null;
        const(char) * sg = null;
        const(char) * sh = null;

        for (i=0; i<k && !p.err; i++) {
            e = n.v.bn.n[i];
            if (i == 0) {
                b = mat_node_get_real_matrix(e, p);
            } else if (i == 1) {
                if (e.t != GenSym.str) {
                    node_type_error(t.t, i+1, GenSym.str, e, p);
                } else {
                    sf = e.v.str;
                }
            } else if (null_node(e)) {
                ; /* OK */
            } else if (e.t != GenSym.str) {
                node_type_error(t.t, i+1, GenSym.str, e, p);
            } else if (i == 2) {
                sg = e.v.str;
            } else {
                sh = e.v.str;
            }
        }
        if (!p.err) {
            if (!gretl_vector_get_length(b)) {
                p.err = GretlError.types;
            } else if (!is_function_call(sf) ||
                       (sg != null && !is_function_call(sg)) ||
                       (sh != null && !is_function_call(sh))) {
                p.err = GretlError.types;
            }
        }
        if (!p.err) {
            ret = aux_scalar_node(p);
        }
        if (!p.err) {
            int minimize = alias_reversed(t) ? 1 : 0;

            ret.v.xval = user_NR(b, sf, sg, sh, p.dset,
                                 minimize, p.prn, &p.err);
        }
    } else if (t.t == GenFunc.f_loess) {
        const(double) * y = null;
        const(double) * x = null;
        double bandwidth = 0.5;
        int poly_order = 1;
        Gretlopt opt = Gretlopt.OPT_NONE;

        for (i=0; i<k && !p.err; i++) {
            e = n.v.bn.n[i];
            if (i < 2) {
                if (e.t != GenSym.series) {
                    node_type_error(t.t, i+1, GenSym.series, e, p);
                } else if (i == 0) {
                    y = e.v.xvec;
                } else {
                    x = e.v.xvec;
                }
            } else if (i == 2 || i == 3) {
                if (e.t != GenSym.num && e.t != GenSym.empty) {
                    node_type_error(t.t, i+1, GenSym.num, e, p);
                } else if (i == 2 && e.t == GenSym.num) {
                    poly_order = node_get_int(e, p);
                } else if (e.t == GenSym.num) {
                    bandwidth = e.v.xval;
                }
            } else {
                if (e.t != GenSym.empty && e.t != GenSym.num) {
                    node_type_error(t.t, i+1, GenSym.num, e, p);
                } else {
                    int ival = node_get_int(e, p);

                    if (!p.err && ival != 0) {
                        if (i == 4) {
                            opt |= Gretlopt.OPT_R;
                        } else {
                            opt |= Gretlopt.OPT_O;
                        }
                    }
                }
            }
        }
        if (!p.err) {
            ret = aux_series_node(p);
            if (ret != null) {
                p.err = gretl_loess(y, x, poly_order, bandwidth,
                                    opt, p.dset, ret.v.xvec);
            }
        }
    } else if (t.t == GenFunc.f_ghk) {
        GretlMatrixRaw *[4] M = null;
        GretlMatrixRaw * dP = null;

        for (i=0; i<k && !p.err; i++) {
            e = n.v.bn.n[i];
            if (i < 4) {
                M[i] = node_get_real_matrix(e, p, i, i+1);
            } else {
                /* the optional last argument */
                if (null_node(e)) {
                    ; /* OK */
                } else if (e.t != GenSym.u_addr) {
                    node_type_error(t.t, i+1, GenSym.u_addr, e, p);
                } else {
                    dP = ptr_node_get_matrix(e, p);
                }
            }
        }
        if (!p.err) {
            ret = aux_matrix_node(p);
        }
        if (!p.err) {
            if (dP == null) {
                ret.v.m = gretl_GHK(M[0], M[1], M[2], M[3], &p.err);
            } else {
                ret.v.m = user_matrix_GHK(M[0], M[1], M[2], M[3],
                                          dP, &p.err);
            }
        }
    } else if (t.t == GenFunc.f_quadtab) {
        int order = -1, method = 1;
        double a = NADBL;
        double b = NADBL;

        for (i=0; i<k && !p.err; i++) {
            e = n.v.bn.n[i];
            if (i == 0) {
                order = node_get_int(e, p);
            } else if (!null_or_scalar(e)) {
                node_type_error(t.t, i+1, GenSym.num, e, p);
            } else if (i == 1) {
                method = node_get_int(e, p);
            } else if (i == 2) {
                a = node_get_scalar(e, p);
            } else {
                b = node_get_scalar(e, p);
            }
        }
        if (!p.err) {
            ret = aux_matrix_node(p);
        }
        if (!p.err) {
            ret.v.m = gretl_quadrule_matrix_new(order, method,
                                                a, b, &p.err);
        }
    } else if (t.t == GenFunc.f_irf) {
        int targ = 0, shock = 0;
        double alpha = 0.0;
        GretlBundle * vb = null;

        for (i=0; i<k && !p.err; i++) {
            e = n.v.bn.n[i];
            if (i == 0) {
                targ = node_get_int(e, p);
            } else if (i == 1) {
                shock = node_get_int(e, p);
            } else if (i == 2) {
                /* optional bootstrap alpha */
                if (e.t != GenSym.empty) {
                    alpha = node_get_scalar(e, p);
                }
            } else {
                /* final optional arg must be bundle */
                if (e.t != GenSym.empty && e.t != GenSym.bundle) {
                    node_type_error(t.t, 4, GenSym.bundle, e, p);
                } else if (e.t == GenSym.bundle) {
                    vb = e.v.b;
                }
            }
        }
        if (!p.err) {
            ret = aux_matrix_node(p);
        }
        if (!p.err) {
            /* convert indices to zero-based */
            targ--;
            shock--;
            if (vb != null) {
                ret.v.m = gretl_IRF_from_bundle(vb, targ, shock, alpha,
                                                p.dset, &p.err);
            } else {
                ret.v.m = last_model_get_irf_matrix(targ, shock, alpha,
                                                    p.dset, &p.err);
            }
        }
    } else if (t.t == GenFunc.f_qlrpval) {
        double X2 = NADBL;
        double p1 = 0, p2 = 0;
        int df = 0;

        for (i=0; i<k && !p.err; i++) {
            e = n.v.bn.n[i];
            if (i == 0) {
                X2 = node_get_scalar(e, p);
            } else if (i == 1) {
                df = node_get_int(e, p);
            } else if (i == 2) {
                p1 = node_get_scalar(e, p);
            } else {
                p2 = node_get_scalar(e, p);
            }
        }
        if (!p.err) {
            ret = aux_scalar_node(p);
        }
        if (!p.err) {
            ret.v.xval = QLR_pval(X2, df, p1, p2);
        }
    } else if (t.t == GenFunc.f_bootci) {
        int cnum = -1, method = 1, B = 0;
        int studentize = 1;
        double alpha = NADBL;

        for (i=0; i<k && !p.err; i++) {
            e = n.v.bn.n[i];
            if (i == 0) {
                cnum = node_get_int(e, p);
            } else if (!null_or_scalar(e)) {
                node_type_error(t.t, i+1, GenSym.num, e, p);
            } else if (i == 1) {
                B = node_get_int(e, p);
            } else if (i == 2) {
                alpha = node_get_scalar(e, p);
            } else if (i == 3) {
                method = node_get_int(e, p);
            } else {
                studentize = node_get_int(e, p);
            }
        }
        if (!p.err) {
            ret = aux_matrix_node(p);
        }
        if (!p.err) {
            ret.v.m = last_model_get_boot_ci(cnum, p.dset, B, alpha, method,
                                             studentize, &p.err);
        }
    } else if (t.t == GenFunc.f_bootpval) {
        int cnum = -1, method = 1, B = 0;

        for (i=0; i<k && !p.err; i++) {
            e = n.v.bn.n[i];
            if (i == 0) {
                cnum = node_get_int(e, p);
            } else if (!null_or_scalar(e)) {
                node_type_error(t.t, i+1, GenSym.num, e, p);
            } else if (i == 1) {
                B = node_get_int(e, p);
            } else {
                method = node_get_int(e, p);
            }
        }
        if (!p.err) {
            ret = aux_scalar_node(p);
        }
        if (!p.err) {
            ret.v.xval = last_model_get_boot_pval(cnum, p.dset, B,
                                                  method, &p.err);
        }
    } else if (t.t == GenFunc.f_movavg) {
        const(double) * x = null;
        double d = 0, y0 = NADBL;
        int len = 0, ctrl = -9999;
        int EMA = 0;

        for (i=0; i<k && !p.err; i++) {
            e = n.v.bn.n[i];
            if (i > 1 && null_node(e)) {
                continue; /* OK */
            } else if (i == 0) {
                if (e.t == GenSym.series) {
                    x = e.v.xvec;
                } else {
                    node_type_error(t.t, i+1, GenSym.series, e, p);
                }
            } else if (i == 1) {
                d = node_get_scalar(e, p);
                if (d < 1.0 && d > 0.0) {
                    EMA = 1;
                } else if (d < 1.0) {
                    p.err = GretlError.invarg;
                } else {
                    len = node_get_int(e, p);
                }
            } else if (i == 2) {
                ctrl = node_get_int(e, p);
            } else {
                y0 = node_get_scalar(e, p);
            }
        }
        if (!p.err) {
            ret = aux_series_node(p);
        }
        if (!p.err) {
            if (ctrl == -9999) {
                /* set the respective defaults */
                ctrl = EMA ? 1 : 0;
            }
            if (EMA) {
                p.err = exponential_movavg_series(x, ret.v.xvec, p.dset,
                                                  d, ctrl, y0);
            } else {
                p.err = movavg_series(x, ret.v.xvec, p.dset, len, ctrl);
            }
        }
    } else if (t.t == GenFunc.f_nadarwat) {
        const(double) * x = null;
        const(double) * y = null;
        int LOO = 0;
        double h = 0;
        double trim = NADBL;

        for (i=0; i<k && !p.err; i++) {
            e = n.v.bn.n[i];
            if (i < 2 && e.t != GenSym.series) {
                node_type_error(t.t, i+1, GenSym.series, e, p);
            }
            if (p.err) {
                break;
            }
            if (i == 0) {
                x = e.v.xvec;
            } else if (i == 1) {
                y = e.v.xvec;
            } else if (i == 2) {
                /* set bandwidth? */
                if (null_node(e)) {
                    ; /* OK: will use the default */
                } else {
                    h = node_get_scalar(e, p);
                    if (h < 0 && k > 3) {
                        gretl_errmsg_sprintf(_("Bandwidth must be non-negative").ptr);
                        p.err = GretlError.invarg;
                        break;
                    } else if (h < 0) {
                        /* it's a legacy thing */
                        LOO = 1;
                        h = -h;
                    }
                }
            } else if (i == 3) {
                /* Leave-one-out */
                LOO = node_get_bool(e, p, 0);
            } else if (i == 4) {
                /* trim spec? (overrides legacy "set" value) */
                trim = node_get_scalar(e, p);
            }
        }
        if (!p.err) {
            ret = aux_series_node(p);
        }
        if (!p.err) {
            if (na(trim)) {
                trim = libset_get_double(SetKey.nadarwat_trim);
            }
            p.err = nadaraya_watson(x, y, h, p.dset, LOO,
                                    trim, ret.v.xvec);
        }
    } else if (t.t == GenFunc.f_hyp2f1) {
        GretlMatrixRaw * x = null;
        double[3] a;

        for (i=0; i<k && !p.err; i++) {
            e = n.v.bn.n[i];
            if (i < 3) {
                a[i] = node_get_scalar(e, p);
            } else {
                x = node_get_real_matrix(e, p, 0, 3);
            }
        }
        if (!p.err) {
            ret = aux_matrix_node(p);
        }
        if (!p.err) {
            ret.v.m = gretl_matrix_alloc(x.rows, x.cols);
            if (ret.v.m == null) {
                p.err = GretlError.alloc;
            } else {
                int nel = x.rows * x.cols;
                double xi, y;

                for (i=0; i<nel; i++) {
                    xi = x.val[i];
                    if (xi < -1.0 || xi > 1.0) {
                        y = NADBL;
                    } else {
                        y = hypergeo(a[0], a[1], a[2], xi);
                    }
                    ret.v.m.val[i] = y;
                }
            }
        }
    } else if (t.t == GenFunc.f_chowlin) {
        GretlMatrixRaw * Y = null;
        GretlMatrixRaw * X = null;
        int fac = 0;

        for (i=0; i<k && !p.err; i++) {
            e = n.v.bn.n[i];
            if (i == 0) {
                /* Y matrix */
                if (e.t == GenSym.mat) {
                    Y = e.v.m;
                } else {
                    p.err = e_types(e);
                }
            } else if (i == 1) {
                /* expansion factor */
                fac = node_get_int(e, p);
            } else if (i == 2) {
                /* optional X matrix  */
                if (e.t == GenSym.mat) {
                    X = e.v.m;
                } else if (!null_node(e)) {
                    p.err = e_types(e);
                }
            }
        }
        if (!p.err) {
            ret = aux_matrix_node(p);
        }
        if (!p.err) {
            ret.v.m = matrix_chowlin(Y, X, fac, &p.err);
        }
    } else if (t.t == GenFunc.f_midasmult) {
        GretlBundle * mb = null;
        int cum = 0;
        int idx = 1;

        for (i=0; i<k && !p.err; i++) {
            e = n.v.bn.n[i];
            if (i == 0) {
                if (e.t == GenSym.bundle) {
                    mb = e.v.b;
                } else {
                    p.err = e_types(e);
                }
            } else if (i == 1) {
                /* cumulate? */
                cum = node_get_bool(e, p, 0);
            } else if (!null_node(e)) {
                /* index of MIDAS term? */
                idx = node_get_int(e, p);
            }
        }
        if (!p.err) {
            ret = aux_matrix_node(p);
        }
        if (!p.err) {
            ret.v.m = midas_multipliers(mb, cum, idx, &p.err);
        }
    } else if (t.t == GenFunc.f_tdisagg) {
        TdisaggInfo tdi;
        GretlBundle * b = null;
        GretlBundle * r = null;

        for (i=0; i<k && !p.err; i++) {
            e = n.v.bn.n[i];
            if (i == 0) {
                /* Y: matrix, series or list */
                if (e.t == GenSym.mat) {
                    tdi.Y = e.v.m;
                } else if (e.t == GenSym.series) {
                    tdi.ynum = e.vnum;
                    tdi.yval = e.v.xvec;
                } else if (e.t == GenSym.list) {
                    tdi.ylist = e.v.ivec;
                } else {
                    p.err = e_types(e);
                }
            } else if (i == 1) {
                /* X: matrix, series, list or null */
                if (e.t == GenSym.mat) {
                    tdi.X = e.v.m;
                } else if (e.t == GenSym.series) {
                    tdi.xnum = e.vnum;
                    tdi.xval = e.v.xvec;
                } else if (e.t == GenSym.list) {
                    tdi.xlist = e.v.ivec;
                    if (gretl_is_midas_list(e.v.ivec, p.dset)) {
                        tdi.xmidas = 1;
                    }
                } else if (!null_node(e)) {
                    p.err = e_types(e);
                }
            } else if (i == 2) {
                /* integer expansion factor */
                tdi.efac = node_get_int(e, p);
            } else if (i == 3) {
                /* optional options bundle */
                if (e.t == GenSym.bundle) {
                    b = e.v.b;
                } else if (!null_node(e)) {
                    p.err = e_types(e);
                }
            } else if (e.t == GenSym.bundle) {
                /* optional retrieval bundle */
                r = e.v.b;
            } else if (!null_node(e)) {
                p.err = e_types(e);
            }
        }
        if (!p.err) {
            ret = aux_matrix_node(p);
        }
        if (!p.err) {
            tdi.targ = gretl_type_from_gen_type(p.targ);
            ret.v.m = get_tdisagg_matrix(&tdi, p.dset, b, r,
                                         p.prn, &p.err);
        }
    } else if (t.t == GenFunc.f_commute) {
        int rowdim = 0, coldim = 0;
        int post = 0, add_id = 0;
        GretlMatrixRaw * A = null;

        for (i=0; i<k && !p.err; i++) {
            e = n.v.bn.n[i];
            if (i == 0) {
                if (e.t == GenSym.mat) {
                    A = e.v.m;
                } else {
                    p.err = e_types(e);
                }
            } else if (i == 1) {
                /* row dimension */
                rowdim = node_get_int(e, p);
                coldim = rowdim;
            } else if (i == 2 && !null_node(e)) {
                /* column dimension */
                coldim = node_get_int(e, p);
            } else if (i == 3) {
                /* postmultiply instead of premultiply? */
                post = node_get_bool(e, p, 0);
            } else if (i == 4) {
                /* add identity matrix? */
                add_id = node_get_bool(e, p, 0);
            }
        }
        if (!p.err) {
            ret = aux_matrix_node(p);
        }
        if (!p.err) {
            ret.v.m = gretl_matrix_commute(A, rowdim, coldim, !post, add_id, &p.err);
        }
    } else if (t.t == GenFunc.f_toepsolv) {
        GretlMatrixRaw * m1 = node_get_real_matrix(n.v.bn.n[0], p, 0, 1);
        GretlMatrixRaw * m2 = node_get_real_matrix(n.v.bn.n[1], p, 1, 2);
        GretlMatrixRaw * m3 = node_get_real_matrix(n.v.bn.n[2], p, 2, 3);
        UserVar * uv = null;

        if (k == 4) {
            uv = ptr_node_get_uvar(n.v.bn.n[3], GenSym.num, p);
        }
        if (!p.err) {
            ret = aux_matrix_node(p);
        }
        if (!p.err) {
            double d;
            double * pdet = (uv != null)? &d : null;

            ret.v.m = gretl_toeplitz_solve(m1, m2, m3, 0, pdet, &p.err);
            if (!p.err && uv != null) {
                user_var_set_scalar_value(uv, d);
            }
        }
    } else if (t.t == GenFunc.f_rgbmix) {
        guint32[2] c = 0;
        double * f = null;
        int do_plot = 0;
        int nf = 0;

        for (i=0; i<k && !p.err; i++) {
            e = n.v.bn.n[i];
            if (i < 2) {
                /* color specs */
                if (e.t == GenSym.num) {
                    c[i] = node_get_guint32(e, p);
                } else if (e.t == GenSym.str) {
                    c[i] = numeric_color_from_string(e.v.str, &p.err);
                } else {
                    p.err = GretlError.types;
                }
            } else if (i == 2) {
                /* mixing spec(s) */
                if (e.t == GenSym.mat) {
                    nf = gretl_vector_get_length(e.v.m);
                    if (nf == 0) {
                        p.err = GretlError.invarg;
                    } else {
                        f = e.v.m.val;
                    }
                } else {
                    p.err = GretlError.types;
                }
            } else {
                /* plot option */
                do_plot = node_get_bool(e, p, 0);
            }
        }
        if (!p.err) {
            ret = aux_array_node(p);
        }
        if (!p.err) {
            ret.v.a = colormix_array(c[0], c[1], f, nf, do_plot, &p.err);
        }
    } else if (t.t == GenFunc.hf_felogitr) {
        GretlMatrixRaw * U = null;
        GretlMatrixRaw * X = null;
        int tt = 0, T = 0;

        for (i=0; i<k && !p.err; i++) {
            e = n.v.bn.n[i];
            if (i == 0) {
                /* number of ones for unit i */
                tt = node_get_int(e, p);
            } else if (i == 1) {
                /* number of observations for unit i */
                T = node_get_int(e, p);
            } else if (i == 2) {
                /* T x 1 index function */
                if (e.t == GenSym.mat) {
                    U = e.v.m;
                } else {
                    p.err = e_types(e);
                }
            } else if (i == 3) {
                /* T x k matrix of covariates */
                if (e.t == GenSym.mat) {
                    X = e.v.m;
                } else {
                    p.err = e_types(e);
                }
            }
        }
        if (!p.err) {
            /* @t < 0 or @T exceeding the size of @U or @X would make
               felogit_rec_loglik() read out of bounds; @t > @T is a
               legitimate case handled internally (returns a zero
               vector) so it is not rejected here. Added from Claude
               2026-08-20.
            */
            if (U == null || X == null || tt < 0 ||
                T > gretl_vector_get_length(U) || T > X.rows) {
                p.err = GretlError.invarg;
            }
        }
        if (!p.err) {
            ret = aux_matrix_node(p);
            ret.v.m = felogit_rec_loglik(tt, T, U, X);
        }
    }

    return ret;
}

GretlType get_type_for_array(Node * t, Parser * p)
{
    if (p.lh.gtype > 0) {
        if (t == p.tree) {
            return p.lh.gtype;
        } else if (p.tree.t == GenSym.query &&
                   (t == p.tree.M || t == p.tree.R)) {
            return p.lh.gtype;
        }
    }

    return GretlType.none;
}

/* evaluate an object definition function (bundle, array, list) */

Node * object_def_node(Node * t, Node * n, Parser * p)
{
    Node * e;
    Node * ret = null;
    GretlType ptype = cast(GretlType) 0;
    int i, f = t.t;
    int k = n.v.bn.n_nodes;

    if ((f == GenFunc.f_defbundle || f == GenFunc.f_defargs) && k % 2 != 0) {
        gretl_errmsg_sprintf(_("%s: requires an even number of arguments").ptr,
                             getsymb(f));
        p.err = GretlError.invarg;
        return null;
    }

    if (f == GenFunc.f_defarray) {
        /* If genr simply defines an array, a type is needed:
           we can get this from p->lh.gtype or from the first
           argument to defarray().
        */
        ptype = get_type_for_array(t, p);
        if (ptype == 0 && k == 0) {
            gretl_errmsg_set(_("defarray: no type was specified").ptr);
            p.err = GretlError.data;
            return null;
        }
    }

    if (f == GenFunc.f_defarray) {
        GretlArray * A;
        int donate = 0;
        void * ptr;

        if (ptype > 0) {
            /* respect the array type specified on input */
            A = gretl_array_new(ptype, 0, &p.err);
        } else {
            /* we'll get the type from the first arg */
            A = gretl_array_new(GretlType.any, 0, &p.err);
        }

        for (i=0; i<k && !p.err; i++) {
            donate = 0;
            e = n.v.bn.n[i];
            p.err = check_array_element_type(A, e);
            if (p.err) {
                break;
            }
            if (e.t == GenSym.num) {
                ptr = gretl_matrix_from_scalar(e.v.xval);
                donate = 1;
            } else {
                ptr = node_get_ptr(e, t.t, p, &donate);
            }
            if (donate) {
                /* copying is not required */
                p.err = gretl_array_append_object(A, ptr, 0);
            } else {
                p.err = gretl_array_append_object(A, ptr, 1);
            }
        }
        if (p.err) {
            gretl_array_destroy(A);
        } else {
            ret = aux_array_node(p);
            if (ret != null) {
                ret.v.a = A;
            }
        }
    } else if (f == GenFunc.f_defbundle || f == GenFunc.f_defargs) {
        GretlBundle * b = gretl_bundle_new();
        GretlType gtype;
        char * key = null;
        void * ptr;
        int donate = 0;

        if (b == null) {
            p.err = GretlError.alloc;
        }
        for (i=0; i<k && !p.err; i++) {
            ptr = null;
            e = n.v.bn.n[i];
            if (i == 0 || i % 2 == 0) {
                /* we need a key string */
                if (e.t == GenSym.str) {
                    key = e.v.str;
                } else {
                    p.err = e_types(e);
                }
            } else {
                /* we need some valid content */
                ptr = node_get_ptr(e, t.t, p, &donate);
                gtype = gretl_type_from_gen_type(e.t);
                if (gtype == GretlType.series) {
                    if (donate) {
                        gretl_bundle_donate_series(b, key, cast(double *) ptr, p.dset);
                    } else {
                        gretl_bundle_set_series(b, key, cast(double *) ptr, p.dset);
                    }
                } else if (type_can_be_bundled(gtype)) {
                    if (donate) {
                        gretl_bundle_donate_data(b, key, ptr, gtype);
                    } else {
                        gretl_bundle_set_data(b, key, ptr, gtype);
                    }
                } else {
                    p.err = e_types(e);
                }
            }
        }
        if (p.err) {
            gretl_bundle_destroy(b);
        } else {
            ret = aux_bundle_node(p);
            if (ret != null) {
                ret.v.b = b;
            }
        }
    } else if (f == GenFunc.f_deflist) {
        int * li;
        int * full_list = gretl_list_new(0);

        for (i=0; i<k && !p.err; i++) {
            li = null;
            e = n.v.bn.n[i];
            if (k == 1 && e.t == GenSym.array) {
                li = list_from_strings_array(e.v.a, p);
            } else if (ok_list_node_plus(e)) {
                li = node_get_list(e, p);
            } else if (e.t == GenSym.mat) {
                li = node_get_list(e, p);
            } else {
                p.err = e_types(e);
            }
            if (!p.err && li[0] > 0) {
                p.err = gretl_list_append_list(&full_list, li);
            }
            free(li);
        }
        if (p.err) {
            free(full_list);
        } else {
            ret = aux_list_node(p);
            if (ret != null) {
                ret.v.ivec = full_list;
            }
        }
    }

    return ret;
}

/* helper for fevalb() in case a bundled series is found */

Node * fevalb_get_bundled_series(GretlBundle * b,
                                 const(char) * key,
                                 Parser * p)
{
    Node * ret = null;
    GretlMatrixRaw * m = null;
    int s;

    s = sample_size(p.dset);
    if (s == 0) {
        no_data_error(p);
        return null;
    }

    m = gretl_bundle_get_series(b, key, null);
    if (m == null) {
        p.err = GretlError.data;
    } else if (m.rows != p.dset.n && m.rows != s) {
        p.err = GretlError.nonconf;
    }

    if (!p.err) {
        int v = current_series_index(p.dset, key);

        if (v < 0) {
            /* we need a new series */
            p.err = dataset_add_NA_series(p.dset, 1);
            if (!p.err) {
                p.err = rename_series(p.dset, p.dset.v - 1,
                                      key, Gretlopt.OPT_NONE);
            }
        }
        if (!p.err) {
            double * dest = p.dset.Z[v] + p.dset.t1;
            double * src = m.val;
            size_t sz = s * (*src).sizeof;

            if (m.rows == p.dset.n) {
                src += p.dset.t1;
            }
            memcpy(dest, src, sz);
            ret = new_node(GenSym.series);
            ret.vnum = v;
            ret.v.xvec = p.dset.Z[v];
            ret.vname = gretl_strdup(key);
        }
    }

    return ret;
}

/* helper function for fevalb() */

Node ** multi_node_from_bundle(GretlBundle * b, int * pargc,
                               Parser * p)
{
    GretlArray * keys = null;
    GretlArray * args = null;
    Node ** nn = null;
    int argc;
    int i;

    keys = gretl_bundle_get_keys(b, &p.err);
    if (keys == null) {
        return null;
    }
    argc = gretl_array_get_length(keys);

    /* check for ordering array */
    if (gretl_bundle_has_key(b, "arglist")) {
        args = gretl_bundle_get_array(b, "arglist", &p.err);
        if (!p.err && !arglist_validate(keys, args)) {
            p.err = GretlError.invarg;
        }
        if (p.err) {
            gretl_array_destroy(keys);
            return null;
        }
        argc -= 1;
    }

    nn = cast(Node **) calloc(argc, (*nn).sizeof);

    if (nn == null) {
        p.err = GretlError.alloc;
    } else {
        Node * save_aux = p.aux;
        char ** S = null;
        GretlType t;
        Node bn;
        Node sn;
        int ns = 0;

        bn.t = GenSym.bundle;
        bn.v.b = b;
        sn.t = GenSym.str;
        p.aux = null;

        if (args != null) {
            S = gretl_array_get_strings(args, &ns);
        } else {
            S = gretl_array_get_strings(keys, &ns);
            strings_array_sort(&S, &ns, Gretlopt.OPT_NONE);
        }

        p.aux = null;
        for (i=0; i<ns; i++) {
            t = gretl_bundle_get_member_type(b, S[i], null);
            if (t == GretlType.series) {
                nn[i] = fevalb_get_bundled_series(b, S[i], p);
            } else {
                sn.v.str = S[i];
                nn[i] = get_bundle_member(&bn, &sn, p);
            }
            if (p.err) {
                break;
            }
        }
        p.aux = save_aux;
    }

    if (!p.err) {
        *pargc = argc;
    } else if (nn != null) {
        for (i=0; i<argc; i++) {
            if (nn[i] != null) {
                free_node(nn[i], p);
            }
        }
        free(nn);
        nn = null;
    }

    gretl_array_destroy(keys);

    return nn;
}

/* In the "traditional" feval() case @l is a multi-node holding the
   function name plus arguments and @r is not referenced.

   In the bundlized fevalb() case, @l holds the name of the function
   to be called and @r holds a bundle of arguments,
*/

Node * eval_feval(int f, Node * l, Node * r, Parser * p)
{
    Node * ret = null;
    Node ** nn = null;
    char * fname = null;
    Ufunc * u = null;
    int argc = 0;
    int fid = 0;

    if (f == GenFunc.f_fevalb) {
        GretlBundle * argb = r.v.b;

        fname = l.v.str;
        nn = multi_node_from_bundle(argb, &argc, p);
    } else {
        int k = l.v.bn.n_nodes;
        Node * e = null;

        if (k < 1) {
            p.err = GretlError.args;
        } else {
            e = l.v.bn.n[0];
            if (e.t != GenSym.str) {
                node_type_error(f, 1, GenSym.str, e, p);
            }
        }
        if (!p.err) {
            argc = k - 1;
            fname = e.v.str;
            nn = l.v.bn.n + 1; /* array of arg nodes */
        }
    }

    if (p.err) {
        return null;
    }

    /* first try for a built-in function */
    fid = function_lookup(fname);
    if (fid != 0) {
        Node tmp;
        Node aux;
        Node mn;
        int np;

        /* temporary function node */
        tmp.t = cast(short) fid;
        if (fid < GenFunc.fp_max) {
            tmp.v.ptr = get_genr_function_pointer(fid);
        }

        np = func1_symb(fid) ? 1 : func2_symb(fid) ? 2 :
            func3_symb(fid) ? 3 : -1;
        aux.t = GenSym.empty;

        if (np > 0) {
            /* known max number of arguments */
            if (argc > np) {
                gretl_errmsg_sprintf("%s: too many arguments".ptr, fname);
                p.err = GretlError.data;
            } else if (np == 1) {
                tmp.L = argc > 0 ? nn[0] : &aux;
            } else if (np == 2) {
                tmp.L = argc > 0 ? nn[0] : &aux;
                tmp.R = argc > 1 ? nn[1] : &aux;
            } else if (np == 3) {
                tmp.L = argc > 0 ? nn[0] : &aux;
                tmp.M = argc > 1 ? nn[1] : &aux;
                tmp.R = argc > 2 ? nn[2] : &aux;
            }
        } else {
            /* multi-arg function */
            mn.t = GenSym.fargs;
            mn.v.bn.n_nodes = argc;
            mn.v.bn.n = nn;
            tmp.L = &mn;
        }
        if (!p.err) {
            ret = eval(&tmp, p);
        }
        if (tmp.aux != null) {
            /* transfer guardianship to the parent of @l */
            if (l.parent.aux != null) {
                free_node(l.parent.aux, p);
            }
            l.parent.aux = tmp.aux;
            tmp.aux = null;
        }
    }

    if (!p.err && fid == 0) {
        /* try for a user function */
        u = get_user_function_by_name(fname);
        if (u != null) {
            Node tmp;
            Node mn;

            tmp.t = GenSym.ufun;
            tmp.vname = fname;
            tmp.v.ptr = u;
            mn.v.bn.n_nodes = argc;
            mn.v.bn.n = nn;
            tmp.R = &mn;
            ret = eval_ufunc(&tmp, &mn, p);
        }
    }

    if (!p.err && fid == 0 && u == null) {
        gretl_errmsg_sprintf(_("%s: function not found").ptr, fname);
        p.err = GretlError.data;
    }

    if (f == GenFunc.f_fevalb) {
        int i;

        for (i=0; i<argc; i++) {
            free_node(nn[i], p);
        }
        free(nn);
    }

    return ret;
}

/* try to get a matrix from @n, even if it's not in fact a
   matrix node as such, provided we can make a matrix out
   of its content
*/

GretlMatrixRaw * node_get_matrix_lenient(Node * n,
                                      int ok,
                                      Parser * p)
{
    GretlMatrixRaw * m = null;

    if (n.t == GenSym.num && (ok == GenSym.num)) {
        m = gretl_matrix_from_scalar(n.v.xval);
    } else if (n.t == GenSym.series && (ok == GenSym.series)) {
        m = gretl_vector_from_series(n.v.xvec, p.dset.t1,
                                     p.dset.t2);
    } else if (n.t == GenSym.list && (ok == GenSym.series)) {
        m = gretl_matrix_data_subset(n.v.ivec, p.dset,
                                     p.dset.t1, p.dset.t2,
                                     MMissingCode.ok, &p.err);
    } else {
        p.err = e_types(n);
    }

    if (!p.err && m == null) {
        p.err = GretlError.alloc;
    }

    return m;
}

GretlBundle * get_kalman_bundle_arg(Node * n, Parser * p)
{
    GretlBundle * b = null;
    Node * e = null;

    if (n.t == GenSym.fargs) {
        /* multi-arguments node */
        e = n.v.bn.n[0];
        e = e.L;
    } else if (n.t == GenSym.u_addr) {
        e = ptr_node_get_referent_node(n, p);
    } else if (n.t == GenSym.bundle) {
        e = n;
    }

    if (e == null || e.t != GenSym.bundle) {
        p.err = e_types(e);
    } else {
        b = e.v.b; /* get the actual bundle */
        if (gretl_bundle_get_type(b) != BundleType.kalman ||
            gretl_bundle_get_private_data(b) == null) {
            p.err = GretlError.types;
            b = null;
        }
    }

    if (p.err) {
        gretl_errmsg_set("Argument 1 must point to a state-space bundle");
    }

    return b;
}

Node * eval_kalman_bundle_func(Node * t, Node * n, Parser * p)
{
    Node * e;
    Node * ret = null;
    int i, k = n.v.bn.n_nodes;

    if (t.t == GenFunc.f_ksetup) {
        GretlMatrixRaw *[5] M = null;
        int[5] copy = 0;
        int dkvar = 0;
        int nmats = 0;

        if (k < 4) {
            n_args_error(k, 4, 4, t.t, p);
        } else if (k >= 5) {
            nmats = 5;
            if (k > 5) {
                dkvar = node_get_bool(n.v.bn.n[5], p, 0);
            }
        } else {
            nmats = 4;
        }
        for (i=0; i<nmats && !p.err; i++) {
            e = n.v.bn.n[i];
            if (e.t == GenSym.mat) {
                M[i] = mat_node_get_real_matrix(e, p);
                if (!p.err) {
                    if (is_tmp_node(e)) {
                        e.v.m = null;
                    } else {
                        copy[i] = 1;
                    }
                }
            } else if (i == 0) {
                /* obsy: accept series or list */
                M[i] = node_get_matrix_lenient(e, GenSym.series, p);
            } else {
                /* system matrices, state variance */
                M[i] = node_get_matrix_lenient(e, GenSym.num, p);
            }
        }
        if (!p.err) {
            GretlBundle * b = kalman_bundle_new(M.ptr, copy.ptr, nmats, dkvar, &p.err);

            if (!p.err) {
                ret = aux_bundle_node(p);
                if (ret != null) {
                    ret.v.b = b;
                }
            }
        }
    } else if (t.t == GenFunc.f_kfilter) {
        GretlBundle * b = get_kalman_bundle_arg(n, p);

        if (!p.err && k != 1) {
            n_args_error(k, 1, 1, t.t, p);
        }
        if (!p.err) {
            ret = aux_scalar_node(p);
        }
        if (!p.err) {
            ret.v.xval = kalman_bundle_filter(b, p.prn, &p.err);
        }
    } else if (t.t == GenFunc.f_kdsmooth) {
        GretlBundle * b = get_kalman_bundle_arg(n, p);
        int dist = 1;

        if (!p.err) {
            if (k == 2) {
                e = n.v.bn.n[1];
                if (node_get_bool(e, p, 0)) {
                    dist = 2;
                }
            } else if (k < 1 || k > 2) {
                n_args_error(k, 1, 2, t.t, p);
            }
        }
        if (!p.err) {
            ret = aux_scalar_node(p);
            if (!p.err) {
                ret.v.xval = kalman_bundle_smooth(b, dist, p.prn);
            }
        }
    } else if (t.t == GenFunc.f_ksmooth) {
        GretlBundle * b = get_kalman_bundle_arg(n, p);

        if (!p.err && k != 1) {
            n_args_error(k, 1, 1, t.t, p);
        }
        if (!p.err) {
            ret = aux_scalar_node(p);
            if (!p.err) {
                ret.v.xval = kalman_bundle_smooth(b, 0, p.prn);
            }
        }
    } else if (t.t == GenFunc.f_ksimul) {
        /* we need a bundle pointer, a matrix,
           and perhaps an optional boolean
        */
        GretlBundle * b = get_kalman_bundle_arg(n, p);
        GretlMatrixRaw * U = null;
        int freeU = 0;
        int get_state = 0;

        if (!p.err && k != 2 && k != 3) {
            n_args_error(k, 2, 3, t.t, p);
        }
        for (i=1; i<k && !p.err; i++) {
            e = n.v.bn.n[i];
            if (i == 1) {
                if (e.t == GenSym.mat) {
                    U = mat_node_get_real_matrix(e, p);
                } else {
                    U = node_get_matrix_lenient(e, GenSym.series, p);
                    if (U != null) {
                        freeU = 1;
                    }
                }
            } else {
                get_state = node_get_int(e, p);
            }
        }
        if (!p.err) {
            ret = aux_matrix_node(p);
        }
        if (!p.err) {
            ret.v.m = kalman_bundle_simulate(b, U, get_state,
                                             p.prn, &p.err);
        }
        if (freeU) {
            gretl_matrix_free(U);
        }
    }

    return ret;
}

Node * kalman_data_node(Node * l, Node * r, Parser * p)
{
    GretlBundle * b = get_kalman_bundle_arg(l, p);
    Node * ret = null;

    if (!p.err) {
        ret = aux_matrix_node(p);
    }
    if (!p.err) {
        ret.v.m = kalman_bundle_simdata(b, r.v.m, p.prn, &p.err);
    }

    return ret;
}

/* Create a matrix using selected series, or a mixture of series and
   lists, or more than one list.  We proceed by setting up a "dummy"
   dataset and constructing a list that indexes into it.  (We can't
   use a regular list, in the general case, since some of the series
   may be temporary variables that are not part of the "real"
   dataset.)
*/

GretlMatrixRaw * assemble_matrix(Node ** nn, int nnodes, Parser * p)
{
    Node * n;
    GretlMatrixRaw * m = null;
    const(int) * list;
    double ** Z = null;
    int * dumlist;
    int i, j, k = 0;

    if (nnodes == 1 && get_matrix_mask() == null) {
        /* take a shortcut if we just got a single series
           and there's no "matrix mask" in place
        */
        n = nn[0];
        if (n.t == GenSym.series) {
            m = series_to_matrix(n.v.xvec, p);
            return m;
        }
    }

    /* how many columns will we need? */
    for (i=0; i<nnodes; i++) {
        n = nn[i];
        if (n.t == GenSym.list) {
            k += n.v.ivec[0];
        } else if (n.t == GenSym.series) {
            k++;
        }
    }

    /* create dummy data array */
    Z = cast(double **) malloc(k * (*Z).sizeof);
    if (Z == null) {
        p.err = GretlError.alloc;
        return null;
    }

    /* and a list associated with Z */
    dumlist = gretl_consecutive_list_new(0, k-1);
    if (dumlist == null) {
        p.err = GretlError.alloc;
        free(Z);
        return null;
    }

    /* attach series pointers to Z */
    k = 0;
    for (i=0; i<nnodes; i++) {
        n = nn[i];
        if (n.t == GenSym.list) {
            list = n.v.ivec;
            for (j=1; j<=list[0]; j++) {
                Z[k++] = p.dset.Z[list[j]];
            }
        } else if (n.t == GenSym.series) {
            Z[k++] = n.v.xvec;
        }
    }

    if (!p.err) {
        DatasetRaw dumset;

        dumset.Z = Z;
        dumset.v = k;
        dumset.n = p.dset.n;
        dumset.t1 = p.dset.t1;
        dumset.t2 = p.dset.t2;

        m = real_matrix_from_list(dumlist, &dumset, p);
    }

    free(dumlist);
    free(Z);

    return m;
}

bool ok_matdef_sym(int s)
{
    return (s == GenSym.num || s == GenSym.series || s == GenSym.empty ||
            s == GenSym.dum || s == GenSym.list || s == GenSym.array);
}

/* composing a matrix from scalars, series or lists */

Node * matrix_def_node(Node * t, Node * n, Parser * p)
{
    GretlMatrixRaw * M = null;
    Node * e;
    Node * ret = null;
    int k = n.v.bn.n_nodes;
    int nnum = 0, nvec = 0;
    int dum = 0, nsep = 0;
    int nlist = 0;
    int seppos = -1;
    int i;

    if (autoreg(p)) {
        fprintf(stderr, "You can't define a matrix in this context\n".ptr);
        p.err = GretlError.types;
        return null;
    }

    for (i=0; i<k && !p.err; i++) {
        e = n.v.bn.n[i];
        if (!ok_matdef_sym(e.t) && !scalar_matrix_node(e)) {
            p.err = e_types(e);
            break;
        }
        if (scalar_node(e)) {
            nnum++;
        } else if (e.t == GenSym.series) {
            nvec++;
        } else if (e.t == GenSym.dum) {
            dum++;
        } else if (e.t == GenSym.list) {
            nlist++;
        } else if (e.t == GenSym.empty) {
            if (nsep == 0) {
                seppos = i;
            }
            nsep++;
        }
        if (dum && k != 1) {
            /* dummy, array must be singleton nodes */
            p.err = GretlError.types;
        } else if ((nvec || nlist) && nnum) {
            /* can't mix series/lists with scalars */
            p.err = GretlError.types;
        } else if ((nvec || nlist) && nsep) {
            /* can't have row separators in a matrix
               composed of series or lists */
            p.err = GretlError.types;
        }
    }

    if (!p.err) {
        Node ** nn = n.v.bn.n;

        if (nvec > 0 || nlist > 1) {
            M = assemble_matrix(nn, k, p);
        } else if (nnum > 0) {
            M = matrix_from_scalars(nn, k, nsep, seppos, p);
        } else if (nlist) {
            M = matrix_from_list(nn[0], p);
        } else if (dum) {
            n = nn[0];
            if (n.v.idnum == DumId.dataset) {
                if (gretl_function_depth() > 0) {
                    gretl_errmsg_set("'dataset' is not recognized as a list within functions");
                    p.err = GretlError.data;
                } else {
                    M = matrix_from_list(null, p);
                }
            } else {
                pprintf(p.prn, "Wrong sort of dummy var\n".ptr);
                p.err = GretlError.types;
            }
        } else {
            /* empty matrix def */
            M = gretl_null_matrix_new();
        }
    }

    if (p.err) {
        if (M != null) {
            gretl_matrix_free(M);
        }
    } else {
        ret = aux_matrix_node(p);
        if (ret != null) {
            ret.v.m = M;
        } else {
            gretl_matrix_free(M);
        }
    }

    return ret;
}

Node * gen_series_from_string(Node * l, Node * r, Parser * p)
{
    Node * ret = null;
    char * line;
    int vnum = -1;
    int err = 0;

    line = g_strdup_printf("%s=%s".ptr, l.v.str, r.v.str);
    err = generate(line, p.dset, GretlType.series,
                   Gretlopt.OPT_U, p.prn);

    if (!err) {
        vnum = current_series_index(p.dset, l.v.str);
    }

    ret = aux_scalar_node(p);
    if (ret != null) {
        ret.v.xval = vnum;
    }

    g_free(line);

    return ret;
}

const(double) * xvec_from_matrix(GretlMatrixRaw * m,
                                 Parser * p,
                                 int * subset,
                                 int * err)
{
    const(double) * ret = null;

    if (gretl_is_null_matrix(m) || m.cols != 1) {
        *err = GretlError.types;
    } else if (m.rows == p.dset.n) {
        ret = m.val;
    } else if (m.rows == sample_size(p.dset)) {
        *subset = 1;
        ret = m.val;
    } else {
        *err = GretlError.types;
    }

    return ret;
}

Node * gen_series_node(Node * l, Node * r, Parser * p)
{
    Node * ret = null;

    if (p.dset == null || p.dset.n == 0) {
        no_data_error(p);
    } else if (l.t == GenSym.str && r.t == GenSym.str) {
        return gen_series_from_string(l, r, p);
    } else if (l.t != GenSym.str || (r.t != GenSym.series && r.t != GenSym.mat && r.t != GenSym.num)) {
        p.err = GretlError.types;
    } else {
        char * vname = l.v.str;
        int v = current_series_index(p.dset, vname);
        const(double) * xvec = null;
        double xt = NADBL;
        int subset = 0;

        if (v == 0) {
            /* can't overwite the constant */
            p.err = GretlError.invarg;
        } else if (v == -1) {
            /* no series called @vname: can we create one? */
            if (user_var_get_type_by_name(vname) != GretlType.none) {
                p.err = GretlError.types;
            } else {
                p.err = check_varname(vname);
            }
        }

        if (!p.err) {
            /* get the data source */
            if (r.t == GenSym.series) {
                xvec = r.v.xvec;
            } else if (r.t == GenSym.mat) {
                xvec = xvec_from_matrix(r.v.m, p, &subset, &p.err);
            } else {
                xt = r.v.xval;
            }
        }

        if (!p.err && v == -1) {
            /* add a new series if required */
            p.err = dataset_add_NA_series(p.dset, 1);
            if (!p.err) {
                v = p.dset.v - 1;
                strcpy(p.dset.varname[v], vname);
            }
        }

        if (!p.err) {
            /* copy data across */
            int t, i = 0;

            for (t=p.dset.t1; t<=p.dset.t2; t++) {
                if (xvec != null) {
                    xt = subset ? xvec[i++] : xvec[t];
                }
                p.dset.Z[v][t] = xt;
            }
        }

        if (!p.err) {
            ret = aux_scalar_node(p);
            if (ret != null) {
                ret.v.xval = v;
            }
        }
    }

    return ret;
}

Node * gen_array_node(Node * t, Node * n, Parser * p)
{
    GretlType ptype;
    Node * ret = null;

    ptype = get_type_for_array(t, p);

    if (!null_or_scalar(n)) {
        p.err = e_types(n);
    } else if (ptype > 0) {
        ; /* OK, we have a type specification */
    } else if (p.targ != GenSym.unk && p.targ != GenSym.array) {
        p.err = GretlError.types;
    } else {
        gretl_errmsg_set(_("array: no type was specified").ptr);
        p.err = GretlError.data;
    }
    if (!p.err) {
        int len = null_node(n) ? 0 : node_get_int(n, p);

        if (!p.err) {
            ret = aux_array_node(p);
            if (!p.err) {
                ret.v.a = gretl_array_new(ptype, len, &p.err);
            }
        }
    }

    return ret;
}

Node * get_series_stringvals(Node * l, Node * r, Parser * p)
{
    Node * ret = aux_array_node(p);

    if (!p.err) {
        int v = l.vnum;

        if (strvals_node(l, p)) {
            int sub = node_get_bool(r, p, 0);
            int n_strs = 0;
            char ** S;

            S = series_get_string_vals(p.dset, v, &n_strs, sub);
            ret.v.a = gretl_array_from_strings(S, n_strs, 1,
                                               &p.err);
        } else {
            ret.v.a = gretl_array_new(GretlType.strings, 0, &p.err);
        }
    }

    return ret;
}

Node * stringify_series(Node * l, Node * r, int f, Parser * p)
{
    Node * ret = aux_scalar_node(p);

    if (ret != null) {
        if (f == GenFunc.f_strvsort && !null_node(r)) {
            ret.v.xval = series_reorder_strings(p.dset, l.vnum, r.v.a);
        } else if (f == GenFunc.f_strvsort) {
            ret.v.xval = series_alphabetize_strings(p.dset, l.vnum);
        } else if (null_node(r)) {
            ret.v.xval = series_destroy_string_table(p.dset, l.vnum);
        } else {
            ret.v.xval = series_set_string_vals(p.dset, l.vnum, r.v.a);
        }
        if (ret.v.xval != 0.0) {
            p.err = GretlError.data;
        }
    }

    return ret;
}

enum Fork {
    l,
    r,
    both,
    none
}

/* Determine whether or not a series is constant in boolean terms,
   i.e. all elements zero, or all non-zero, over the relevant range.
   If so, return FORK_L (all 1) or FORK_R (all 0), otherwise
   return FORK_UNK.
*/

int vec_branch(const(double) * c, Parser * p)
{
    int c1, t, t1, t2;
    int ret;

    t1 = autoreg(p) ? p.obs : p.dset.t1;
    t2 = autoreg(p) ? p.obs : p.dset.t2;

    c1 = (c[t1] != 0.0);
    ret = (c1)? Fork.l : Fork.r;

    for (t=t1; t<=t2; t++) {
        if (!na(c[t])) {
            if ((c1 && c[t] == 0) || (!c1 && c[t] != 0)) {
                ret = Fork.both;
                break;
            }
        }
    }

    return ret;
}

/* Given a series condition in a ternary "?" expression, return the
   evaluated counterpart.  We evaluate both forks and select based on
   the value of the condition at each observation.  We accept only
   scalar (NUM or 1x1 matrix) and series (SERIES) types on input, and
   always produce a SERIES type on output.
*/

Node * query_eval_series(const(double) * c, Node * n, Parser * p)
{
    Node * l = null;
    Node * r = null;
    Node * ret = null;
    Node * save_aux = p.aux;
    double * xvec = null;
    double * yvec = null;
    double x = NADBL, y = NADBL;
    double xt, yt;
    int t, t1, t2;
    int branch;

    branch = vec_branch(c, p);

    if (autoreg(p) || branch != Fork.r) {
        l = eval(n.M, p);
        if (p.err) {
            return null;
        }
        if (l.t == GenSym.series) {
            xvec = l.v.xvec;
        } else if (scalar_node(l)) {
            x = node_get_scalar(l, p);
        } else {
            p.err = e_types(l);
            return null;
        }
    }

    if (autoreg(p) || branch != Fork.l) {
        r = eval(n.R, p);
        if (p.err) {
            return null;
        }
        if (r.t == GenSym.series) {
            yvec = r.v.xvec;
        } else if (scalar_node(r)) {
            y = node_get_scalar(r, p);
        } else {
            p.err = e_types(r);
            return null;
        }
    }

    reset_p_aux(p, save_aux);
    ret = aux_series_node(p);

    t1 = autoreg(p) ? p.obs : p.dset.t1;
    t2 = autoreg(p) ? p.obs : p.dset.t2;

    for (t=t1; t<=t2; t++) {
        if (na(c[t])) {
            ret.v.xvec[t] = NADBL;
        } else {
            xt = (xvec != null)? xvec[t] : x;
            yt = (yvec != null)? yvec[t] : y;
            ret.v.xvec[t] = (c[t] != 0.0)? xt : yt;
        }
    }

    return ret;
}

/* The condition in the ternary query operator is a scalar,
   which has been evaluated to @x, and which must now be
   interpreted as a boolean.
*/

Node * query_eval_scalar(double x, Node * n, Parser * p)
{
    Node * save_aux = p.aux;
    Node * l = null;
    Node * r = null;
    Node * ret = null;
    int branch, gen_t = 0;

    branch = na(x) ? Fork.none : (x != 0 ? Fork.l : Fork.r);

    if (autoreg(p) || branch != Fork.r) {
        l = eval(n.M, p);
        if (p.err) {
            return null;
        } else {
            gen_t = l.t;
        }
    }

    if (autoreg(p) || branch != Fork.l) {
        r = eval(n.R, p);
        if (p.err) {
            return null;
        } else {
            gen_t = r.t;
        }
    }

    if (branch == Fork.none) {
        if (gen_t == GenSym.num) {
            reset_p_aux(p, save_aux);
            ret = aux_scalar_node(p);
            if (ret != null) {
                ret.v.xval = NADBL;
            }
        } else {
            gretl_errmsg_set(_("indeterminate condition for '?' operator").ptr);
            p.err = GretlError.data;
        }
    } else if (branch == Fork.l) {
        ret = l;
    } else if (branch == Fork.r) {
        ret = r;
    }

    return ret;
}

/* The following allows for @n to hold a scalar, real matrix or
   complex matrix, and also allows for the case where @n's payload is
   real-valued but a complex result is required, signalled by @need_z.
*/

int query_term_get_value(Node * n, int i, int j,
                         double * py, Cmplx * pz,
                         int need_z)
{
    if (n.t == GenSym.mat) {
        if (n.v.m == null) {
            fprintf(stderr, "missing matrix in query_term_get_value: n->v.m = %p\n".ptr,
                    cast(void *) n.v.m);
            return GretlError.data;
        } else if (n.v.m.is_complex) {
            *pz = gretl_cmatrix_get(n.v.m, i, j);
        } else {
            *py = gretl_matrix_get(n.v.m, i, j);
            if (need_z) {
                *pz = Cmplx(*py, 0);
            }
        }
    } else {
        *py = n.v.xval;
        if (need_z) {
            *pz = Cmplx(*py, 0);
        }
    }

    return 0;
}

/* the condition in the ternary query operator is a matrix */

Node * query_eval_matrix(GretlMatrixRaw * m, Node * n, Parser * p)
{
    Node * save_aux = p.aux;
    Node * l;
    Node * r;
    Node * ret = null;
    GretlMatrixRaw * mret;
    int lcomplex = 0;
    int rcomplex = 0;

    if (gretl_is_null_matrix(m)) {
        p.err = GretlError.types;
        return null;
    }

    l = eval(n.M, p);
    if (!p.err) {
        r = eval(n.R, p);
    }
    if (p.err) {
        return null;
    }

    if ((l.t != GenSym.num && l.t != GenSym.mat) ||
        (r.t != GenSym.num && r.t != GenSym.mat)) {
        if (gretl_matrix_is_scalar(m)) {
            double x = m.val[0];

            if (na(x)) {
                gretl_errmsg_set(_("indeterminate condition for '?' operator").ptr);
                p.err = GretlError.data;
            } else if (x != 0) {
                ret = l;
            } else {
                ret = r;
            }
        } else {
            p.err = GretlError.types;
        }
        return ret;
    }

    if (l.t == GenSym.mat) {
        if (l.v.m.cols != m.cols || l.v.m.rows != m.rows) {
            p.err = GretlError.nonconf;
            return null;
        } else if (l.v.m.is_complex) {
            lcomplex = 1;
        }
    } else if (r.t == GenSym.mat) {
        if (r.v.m.cols != m.cols || r.v.m.rows != m.rows) {
            p.err = GretlError.nonconf;
            return null;
        } else if (r.v.m.is_complex) {
            rcomplex = 1;
        }
    }

    if (lcomplex || rcomplex) {
        mret = gretl_cmatrix_build(m, null, 0, 0, &p.err);
    } else {
        mret = gretl_matrix_copy(m);
        if (mret == null) {
            p.err = GretlError.alloc;
        }
    }
    if (p.err) {
        return null;
    }

    reset_p_aux(p, save_aux);
    ret = aux_matrix_node(p);

    if (!p.err) {
        int need_z = mret.is_complex;
        Cmplx z = Cmplx(0, 0);
        double x, y = 0;
        int j, i;

        for (j=0; j<m.cols; j++) {
            for (i=0; i<m.rows; i++) {
                x = gretl_matrix_get(m, i, j);
                if (isnan(x)) {
                    if (mret.is_complex) {
                        gretl_cmatrix_set(mret, i, j, Cmplx(x, x));
                    } else {
                        gretl_matrix_set(mret, i, j, x);
                    }
                } else if (x != 0.0) {
                    p.err = query_term_get_value(l, i, j, &y, &z, need_z);
                    if (p.err) {
                        break;
                    }
                    if (mret.is_complex) {
                        gretl_cmatrix_set(mret, i, j, z);
                    } else {
                        gretl_matrix_set(mret, i, j, y);
                    }
                } else {
                    p.err = query_term_get_value(r, i, j, &y, &z, need_z);
                    if (p.err) {
                        break;
                    }
                    if (mret.is_complex) {
                        gretl_cmatrix_set(mret, i, j, z);
                    } else {
                        gretl_matrix_set(mret, i, j, y);
                    }
                }
            }
        }
        ret.v.m = mret;
    }

    return ret;
}

/* Evaluate a ternary "query" expression: C ? X : Y.  The condition C
   may be a scalar, series or matrix. On entry to this function @c will
   hold the evaluated C; the alternates are represented by the M and R
   sub-nodes of @t, not yet evaluated.
*/

Node * eval_query(Node * t, Node * c, Parser * p)
{
    Node * ret = null;

    if (c.t == GenSym.num) {
        ret = query_eval_scalar(c.v.xval, t, p);
    } else if (c.t == GenSym.series) {
        ret = query_eval_series(c.v.xvec, t, p);
    } else if (c.t == GenSym.mat) {
        /* note 2025-03-15: the special handling for 1 x 1 c->v.m,
           defined out below, resulted in breakage elsewhere
        */
        ret = query_eval_matrix(c.v.m, t, p);
    } else {
        /* invalid type for boolean condition */
        p.err = e_types(c);
    }

    return ret;
}

bool dvar_scalar(int i)   { return (i > 0 && i < RetrievalIndex.scalar_max); }
bool dvar_series(int i)   { return (i > RetrievalIndex.scalar_max && i < RetrievalIndex.series_max); }
bool dvar_matrix(int i)   { return (i == RetrievalIndex.now); }
bool dvar_variant1(int i) { return (i == RetrievalIndex.test_stat || i == RetrievalIndex.test_pval); }
bool dvar_variant2(int i) { return (i == RetrievalIndex.result); }

bool no_data(const(DatasetRaw) * p)
{
    return (p == null || p.n == 0);
}

double dvar_get_scalar(int i, const(DatasetRaw) * dset)
{
    switch (i) {
    case RetrievalIndex.nobs:
        return (dset == null) ? NADBL :
            (dset.n == 0 ? 0 : sample_size(dset));
    case RetrievalIndex.nvars:
        return (dset == null)? NADBL : dset.v;
    case RetrievalIndex.pd:
        return (no_data(dset))? NADBL : dset.pd;
    case RetrievalIndex.panel_pd:
        return (dataset_is_panel(dset))? dset.panel_pd : NADBL;
    case RetrievalIndex.t1:
        return (no_data(dset))? NADBL : dset.t1 + 1;
    case RetrievalIndex.t2:
        return (no_data(dset))? NADBL : dset.t2 + 1;
    case RetrievalIndex.tmax:
        if (no_data(dset)) {
            return NADBL;
        } else {
            int tmax;

            sample_range_get_extrema(dset, null, &tmax);
            return tmax + 1;
        }
    case RetrievalIndex.datatype:
        return dataset_get_structure(dset);
    case RetrievalIndex.test_pval:
        return get_last_pvalue();
    case RetrievalIndex.test_stat:
        return get_last_test_statistic();
    case RetrievalIndex.test_lnl:
        return get_last_lnl();
    case RetrievalIndex.test_brk:
        return get_last_break();
    case RetrievalIndex.stopwatch:
        return gretl_stopwatch();
    case RetrievalIndex.windows:
        return 0;
    case RetrievalIndex.version_:
        return gretl_version_number(GRETL_VERSION.ptr);
    case RetrievalIndex.errno:
        return get_gretl_errno();
    case RetrievalIndex.seed:
        return gretl_rand_get_seed();
    case RetrievalIndex.huge:
        return libset_get_double(SetKey.conv_huge);
    case RetrievalIndex.loglevel:
        return libset_get_int(SetKey.loglevel);
    case RetrievalIndex.logstamp:
        return libset_get_bool(SetKey.logstamp);
    default:
        return NADBL;
    }
}

int date_series_ok(const(DatasetRaw) * dset)
{
    if (calendar_data(dset)) {
        return 1;
    } else if (quarterly_or_monthly(dset)) {
        return 1;
    } else if (annual_data(dset) || decennial_data(dset)) {
        return 1;
    } else if (dataset_has_panel_time(dset)) {
        return 1;
    } else {
        return 0;
    }
}

int dvar_get_series(double * x, int i, const(DatasetRaw) * dset)
{
    int tseries = dataset_is_time_series(dset);
    int panel = dataset_is_panel(dset);
    int YMD = calendar_data(dset);
    int t, err = 0;

    if (dset == null || dset.n == 0) {
        return GretlError.nodata;
    }
    if (i == RetrievalIndex.obsmin && dset.pd < 2) {
        return GretlError.pdwrong;
    }
    if (i == RetrievalIndex.obsmic && !YMD) {
        return GretlError.pdwrong;
    }
    if (i == RetrievalIndex.punit && !panel) {
        return GretlError.pdwrong;
    }
    if (i == RetrievalIndex.time && !panel && !tseries) {
        return GretlError.pdwrong;
    }
    if (i == RetrievalIndex.dates && !date_series_ok(dset)) {
        return GretlError.pdwrong;
    }

    if (i == RetrievalIndex.obsmaj) {
        if (dset.pd == 1 && !tseries) {
            i = RetrievalIndex.index;
        } else if (panel) {
            i = RetrievalIndex.punit;
        }
    }

    if (i == RetrievalIndex.index) {
        for (t=0; t<dset.n; t++) {
            x[t] = t + 1;
        }
    } else if (YMD && i != RetrievalIndex.index && i != RetrievalIndex.dates) {
        /* Watch out: we're handling most calendar-data cases
           here, so we have to explicitly exclude cases that
           require different treatment.
        */
        char[12] obs;
        int y, m, d;

        for (t=0; t<dset.n && !err; t++) {
            ntolabel(obs.ptr, t, dset);
            if (sscanf(obs.ptr, YMD_READ_FMT.ptr, &y, &m, &d) != 3) {
                err = GretlError.data;
            } else if (i == RetrievalIndex.obsmaj) {
                x[t] = y;
            } else if (i == RetrievalIndex.obsmin) {
                x[t] = m;
            } else {
                x[t] = d;
            }
        }
    } else if (i == RetrievalIndex.time) {
        if (panel) {
            make_panel_time_var(x, dset);
        } else {
            for (t=0; t<dset.n; t++) {
                x[t] = cast(double) (t + 1);
            }
        }
    } else if (i == RetrievalIndex.punit) {
        for (t=0; t<dset.n; t++) {
            x[t] = t / dset.pd + 1;
        }
    } else if (i == RetrievalIndex.obsmaj) {
        int maj;

        for (t=0; t<dset.n; t++) {
            date_maj_min(t, dset, &maj, null);
            x[t] = maj;
        }
    } else if (i == RetrievalIndex.obsmin) {
        int min;

        for (t=0; t<dset.n; t++) {
            date_maj_min(t, dset, null, &min);
            x[t] = min;
        }
    } else if (i == RetrievalIndex.dates) {
        err = fill_dataset_dates_series(dset, x);
    } else {
        err = GretlError.data;
    }

    return err;
}

GretlMatrixRaw * dvar_get_matrix(int i, int * err)
{
    GretlMatrixRaw * m = null;

    switch (i) {
    case RetrievalIndex.test_stat:
        m = get_last_test_matrix(err);
        break;
    case RetrievalIndex.test_pval:
        m = get_last_pvals_matrix(err);
        break;
    case RetrievalIndex.now:
        m = gretl_matrix_alloc(1, 2);
        if (m == null) {
            *err = GretlError.alloc;
        } else {
            /* package epoch second and ISO 8601 basic date */
            long t = time(null);
            tm tmv;
            int y, mon, d;

            localtime_r(&t, &tmv);
            y = tmv.tm_year + 1900;
            mon = tmv.tm_mon + 1;
            d = tmv.tm_mday;
            m.val[0] = cast(double) t;
            m.val[1] = y * 10000 + mon * 100 + d;
        }
        break;
    default:
        *err = GretlError.data;
        break;
    }

    return m;
}

Node * dollar_var_node(Node * t, Parser * p)
{
    Node * ret = null;

    if (starting(p)) {
        int idx = t.v.idnum;

        if (dvar_scalar(idx)) {
            ret = aux_scalar_node(p);
            if (ret != null) {
                ret.v.xval = dvar_get_scalar(idx, p.dset);
            }
        } else if (dvar_series(idx)) {
            ret = aux_series_node(p);
            if (ret != null) {
                p.err = dvar_get_series(ret.v.xvec, idx, p.dset);
            }
        } else if (dvar_matrix(idx)) {
            ret = aux_matrix_node(p);
            if (ret != null) {
                ret.v.m = dvar_get_matrix(idx, &p.err);
            }
        } else if (idx == RetrievalIndex.pngfont) {
            ret = aux_string_node(p);
            if (ret != null) {
                ret.v.str = gretl_png_font_string();
            }
        } else if (idx == RetrievalIndex.grpnames) {
            ret = aux_array_node(p);
            if (ret != null) {
                ret.v.a = cast(GretlArray *) get_panel_group_names(p.dset, &p.err);
            }
        } else if (idx == RetrievalIndex.mapfile) {
            const(char) * fname = dataset_get_mapfile(p.dset);

            ret = aux_string_node(p);
            if (!p.err) {
                if (fname == null) {
                    ret.v.str = gretl_strdup("");
                } else {
                    ret.v.str = gretl_strdup(fname);
                }
            }
        } else if (idx == RetrievalIndex.map) {
            ret = aux_bundle_node(p);
            if (!p.err) {
                ret.v.b = get_current_map(p.dset, null, &p.err);
            }
        } else if (dvar_variant1(idx)) {
            GretlType type = cast(GretlType) get_last_test_type();

            if (type == GretlType.matrix) {
                ret = aux_matrix_node(p);
                if (ret != null) {
                    ret.v.m = dvar_get_matrix(idx, &p.err);
                }
            } else {
                /* scalar or none */
                ret = aux_scalar_node(p);
                if (ret != null) {
                    ret.v.xval = dvar_get_scalar(idx, p.dset);
                }
            }
        } else if (dvar_variant2(idx)) {
            GretlType type = cast(GretlType) 0;
            void * ptr = get_last_result_data(&type, &p.err);

            if (type == GretlType.matrix) {
                ret = aux_matrix_node(p);
                if (ret != null) {
                    ret.v.m = cast(GretlMatrixRaw *) ptr;
                }
            } else if (type == GretlType.bundle) {
                ret = aux_bundle_node(p);
                if (ret != null) {
                    ret.v.b = cast(GretlBundle *) ptr;
                }
            } else if (!p.err) {
                p.err = GretlError.types;
            }
        }
    } else {
        ret = aux_any_node(p);
    }

    return ret;
}

Node * fevd_node(Node * l, Node * m, Node * r, Parser * p)
{
    Node * ret = aux_matrix_node(p);

    if (ret != null) {
        int targ = node_get_int(l, p);
        int shock = 0;

        if (!p.err && !null_node(m)) {
            shock = node_get_int(m, p);
        }
        if (!p.err && !null_node(r) && r.t != GenSym.bundle) {
            p.err = GretlError.invarg;
        }

        if (!p.err) {
            /* convert @targ, @shock to zero-based */
            targ -= 1;
            shock -= 1;
            if (r.t == GenSym.bundle) {
                ret.v.m = gretl_FEVD_from_bundle(r.v.b, targ, shock,
                                                 p.dset, &p.err);
            } else {
                GretlObjType otype;
                void * ptr;

                ptr = get_last_model(&otype);
                if (ptr == null || otype != GretlObjType.var) {
                    p.err = GretlError.badstat;
                } else {
                    ret.v.m = gretl_VAR_get_FEVD_matrix(cast(GretlVar *) ptr, targ, shock, 0,
                                                        p.dset, &p.err);
                }
            }
        }
    }

    return ret;
}

GretlType object_var_type(int idx, const(char) * oname,
                          int * needs_data)
{
    GretlType vtype = GretlType.none;

    *needs_data = 0;

    if (model_data_scalar(idx)) {
        vtype = GretlType.double_;
    } else if (model_data_series(idx)) {
        vtype = GretlType.series;
    } else if (model_data_matrix(idx)) {
        vtype = GretlType.matrix;
    } else if (model_data_matrix_builder(idx)) {
        vtype = GretlType.matrix;
        *needs_data = 1;
    } else if (model_data_list(idx)) {
        vtype = GretlType.list;
    } else if (model_data_string(idx)) {
        vtype = GretlType.string_;
    } else if (model_data_array(idx)) {
        vtype = GretlType.array;
    } else if (idx == BundleDataIndex.model) {
        vtype = GretlType.bundle;
    } else if (idx == BundleDataIndex.system) {
        vtype = GretlType.bundle;
    }

    if (idx == ModelDataIndex.uhat || idx == ModelDataIndex.yhat || idx == ModelDataIndex.sigma) {
        /* could be a matrix */
        vtype = saved_object_get_data_type(oname, cast(ModelDataIndex) idx);
    }

    return vtype;
}

/* For example, $coeff(sqft); or model1.$coeff(const) */

Node * dollar_str_node(Node * t, Model * pmod, Parser * p)
{
    Node * ret = aux_scalar_node(p);

    if (ret != null && starting(p)) {
        const(char) * str = null;
        Node * l = t.L;
        Node * r = t.R;

        if (r.t == GenSym.str) {
            str = r.v.str;
        } else {
            /* could be element of strings array? */
            Node * e = eval(r, p);

            if (e.t == GenSym.str) {
                str = e.v.str;
            } else {
                p.err = e_types(e);
            }
        }
        if (!p.err) {
            ret.v.xval = gretl_model_get_data_element(pmod, l.v.idnum, str,
                                                      p.dset, &p.err);
        }
        if (p.err && r.t == GenSym.str) {
            pprintf(p.prn, _("'%s': invalid argument for %s()\n").ptr,
                    r.v.str, mvarname(l.v.idnum));
        }
    }

    return ret;
}

/* Retrieve a data item from a model (single equation or system of
   some kind. We're not sure in advance here of the type of the data
   item (scalar, matrix, etc.): we look that up with object_var_type().

   If retrieval is from a named model, the model name will be given
   by the @vname member of @t and the key will be on node @k.
   Otherwise @k will be NULL and @t will already hold the index of
   the required data item.
*/

Node * model_var_node(Node * t, Node * k, Parser * p)
{
    Node * ret = null;

    if (starting(p)) {
        const(char) * oname = null;
        int idx, needs_data = 0;
        GretlType vtype;

        if (t.t == GenSym.dbundle) {
            /* pseudo-bundle: holds name of model */
            oname = t.vname;
            idx = mvar_lookup(k.v.str);
            if (idx == 0) {
                p.err = GretlError.data;
                return null;
            }
        } else {
            /* @t already holds the data-item index */
            idx = t.v.idnum;
        }

        /* determine the type of the requested data */
        vtype = object_var_type(idx, oname, &needs_data);

        if (vtype == GretlType.double_) {
            ret = aux_scalar_node(p);
        } else if (vtype == GretlType.series) {
            ret = aux_series_node(p);
        } else if (vtype == GretlType.list) {
            ret = aux_list_node(p);
        } else if (vtype == GretlType.string_) {
            ret = aux_string_node(p);
        } else if (vtype == GretlType.bundle) {
            ret = aux_bundle_node(p);
        } else if (vtype == GretlType.array) {
            ret = aux_array_node(p);
        } else {
            ret = aux_matrix_node(p);
        }

        if (ret == null) {
            return ret;
        } else if (vtype == GretlType.double_) {
            ret.v.xval = saved_object_get_scalar(oname, idx, p.dset,
                                                 &p.err);
        } else if (vtype == GretlType.series) {
            p.err = saved_object_get_series(ret.v.xvec, oname, idx, p.dset);
        } else if (vtype == GretlType.list) {
            ret.v.ivec = saved_object_get_list(oname, idx, &p.err);
        } else if (vtype == GretlType.string_) {
            ret.v.str = saved_object_get_string(oname, idx, p.dset,
                                                &p.err);
        } else if (vtype == GretlType.bundle) {
            ret.v.b = bundle_from_model(null, p.dset, &p.err);
        } else if (vtype == GretlType.matrix) {
            if (needs_data) {
                ret.v.m = saved_object_build_matrix(oname, idx,
                                                    p.dset, &p.err);
            } else {
                ret.v.m = saved_object_get_matrix(oname, idx, &p.err);
            }
        } else if (vtype == GretlType.array) {
            ret.v.a = cast(GretlArray *) saved_object_get_array(oname, idx, p.dset, &p.err);
        }
    } else {
        ret = aux_any_node(p);
    }

    return ret;
}

Node * wildlist_node(Node * n, Parser * p)
{
    Node * ret = aux_list_node(p);

    if (ret != null && starting(p)) {
        int * list = varname_match_list(p.dset, n.v.str,
                                        &p.err);

        ret.v.ivec = list;
    }

    return ret;
}

Node * ellipsis_list_node(Node * l, Node * r, Parser * p)
{
    Node * ret = aux_list_node(p);

    if (ret != null && starting(p)) {
        int v1 = useries_node(l) ? l.vnum : cast(int) l.v.xval;
        int v2 = useries_node(r) ? r.vnum : cast(int) r.v.xval;

        ret.v.ivec = ellipsis_list(p.dset, v1, v2, &p.err);
    }

    return ret;
}

Node * mat2list_node(Node * l, Node * r, Parser * p)
{
    Node * ret = aux_list_node(p);

    if (ret != null) {
        const(char) * pfx = null_node(r) ? null : r.v.str;

        ret.v.ivec = list_from_matrix(l.v.m, pfx, p.dset, &p.err);
    }

    return ret;
}

/* see if a plain NUM node can be interpreted as holding a
   series ID, in the context of creating a list */

int could_be_series_id(Node * n, Parser * p)
{
    if (n.t == GenSym.num && p.dset != null && p.dset.n > 0) {
        int k = node_get_int(n, p);

        if (!p.err && k >= 0 && k < p.dset.v) {
            return 1;
        }
    }

    return 0;
}

Node * list_join_node(Node * l, Node * r, Parser * p)
{
    Node * ret = aux_list_node(p);

    if (ret != null && starting(p)) {
        int * L1 = node_get_list(l, p);
        int * L2 = node_get_list(r, p);

        if (!p.err) {
            ret.v.ivec = gretl_lists_join_with_separator(L1, L2);
            if (ret.v.ivec == null) {
                p.err = GretlError.alloc;
            }
        }

        free(L1);
        free(L2);
    }

    return ret;
}

Node * two_scalars_func(Node * l, Node * r, int t, Parser * p)
{
    Node * ret = aux_scalar_node(p);

    if (ret != null && starting(p)) {
        double xl = node_get_scalar(l, p);
        double xr = node_get_scalar(r, p);

        if (!na(xl) && !na(xr)) {
            if (t == GenFunc.f_min) {
                ret.v.xval = (xl < xr)? xl : xr;
            } else if (t == GenFunc.f_max) {
                ret.v.xval = (xl > xr)? xl : xr;
            } else if (t == GenFunc.f_randint) {
                int k;

                p.err = gretl_rand_int_minmax(&k, 1, cast(int) xl, cast(int) xr);
                if (!p.err) {
                    ret.v.xval = k;
                }
            }
        }
    }

    return ret;
}

/* Return the main or max value "within" the argument attached
   to the node @n.
*/

Node * within_minmax_node(Node * n, int f, Parser * p)
{
    Node * ret = null;

    if (n.t == GenSym.num) {
        ret = pretend_matrix_scalar_func(n, f, p);
    } else if (n.t == GenSym.series) {
        ret = series_scalar_func(n, f, null, p);
    } else if (n.t == GenSym.mat) {
        if (n.v.m.is_complex) {
            p.err = GretlError.cmplx;
        } else {
            ret = series_scalar_func(n, f, null, p);
        }
    } else if (ok_list_node(n, p)) {
        ret = list_to_series_func(n, f, p);
    } else {
        p.err = GretlError.invarg;
    }

    return ret;
}

/* Return the main or max value "between" the arguments attached
   to the nodes @l and @r.
*/

Node * between_minmax_node(Node * l, Node * r, int f, Parser * p)
{
    Node * ret = null;

    if (scalar_node(l) && scalar_node(r)) {
        ret = two_scalars_func(l, r, f, p);
    } else if (l.t != r.t) {
        p.err = GretlError.types;
    } else if (l.t == GenSym.series) {
        ret = aux_series_node(p);
        if (ret != null) {
            p.err = two_series_minmax(l.v.xvec, r.v.xvec,
                                      ret.v.xvec, f, p.dset);
        }
    } else if (l.t == GenSym.mat) {
        if (l.v.m.is_complex || r.v.m.is_complex) {
            p.err = GretlError.cmplx;
        } else {
            ret = aux_matrix_node(p);
        }
        if (ret != null) {
            ret.v.m = two_matrices_minmax(l.v.m, r.v.m,
                                          f, &p.err);
        }
    } else {
        p.err = GretlError.invarg;
    }

    return ret;
}

Node * kpss_crit_node(Node * l, Node * r, Parser * p)
{
    Node * ret = aux_matrix_node(p);

    if (ret != null && starting(p)) {
        int T = node_get_int(l, p);
        int trend = node_get_int(r, p);

        if (!p.err) {
            ret.v.m = kpss_critvals(T, trend, &p.err);
        }
    }

    return ret;
}

Node * scalar_postfix_node(Node * n, Parser * p)
{
    Node * ret = aux_scalar_node(p);

    if (ret != null && starting(p)) {
        if (n.vname != null && p.lh.name[0] != '\0' &&
            strcmp(n.vname, p.lh.name.ptr) == 0) {
            /* undefined behavior */
            fprintf(stderr, "BAD NUM_P/NUM_M\n".ptr);
            gretl_errmsg_sprintf(_("The result for %s is not well defined").ptr,
                                 n.vname);
            p.err = GretlError.data;
            ret.v.xval = NADBL;
        } else {
            double x = n.v.xval;

            ret.v.xval = x;
            if (n.t == GenSym.num_p) {
                p.err = node_replace_scalar(n, x + 1.0);
            } else {
                p.err = node_replace_scalar(n, x - 1.0);
            }
        }
    }

    return ret;
}

int series_calc_nodes(Node * l, Node * r)
{
    int ret = 0;

    if (l.t == GenSym.series) {
        ret = (r.t == GenSym.series || r.t == GenSym.num || scalar_matrix_node(r));
    } else if (r.t == GenSym.series) {
        ret = scalar_node(l);
    }

    return ret;
}

int cast_series_to_list(Parser * p, Node * n, short f)
{
    if (p.tree.t == GenFunc.f_genseries || p.tree.t == GenSym.ufun) {
        /* any other cases when we shouldn't do this "cast"? */
        return 0;
    } else if (p.targ == GenSym.list && useries_node(n)) {
        return (f == GenFunc.f_log || f == GenFunc.f_diff ||
                f == GenFunc.f_ldiff || f == GenFunc.f_sdiff ||
                f == GenFunc.f_odev);
    } else {
        return 0;
    }
}

Node * lhs_terminal_node(Node * t, Node * l, Node * r,
                         Parser * p)
{
    /* Pass through eval'd @l and @r subnodes, but don't eval
       the parent @t itself */
    Node * ret = aux_parent_node(p);

    ret.t = t.t; /* transcribe type */
    ret.L = l;   /* evaluated left-hand */
    ret.R = r;   /* evaluated right-hand */

    /* prevent double-freeing of children @l and @r */
    ret.flags |= NodeFlags.lht;

    return ret;
}

/* reattach_series: on successive executions of a given
   compiled "genr", the "xvec" pointer recorded on a
   SERIES node will have become invalid if:

   (1) the series has been renamed,
   (2) the dataset has been differently sub-sampled, or
   (3) the series has been deleted (should be impossible).

   In case (2) the ID number of the series should still
   be valid, so it ought to be sufficient to reconnect the
   xvec pointer as in the second branch below. In case (1),
   however, it's necessary to look up the ID number by name
   again (as things stand).

   Note that n->v.xvec will be NULL when this function is
   reached only in case a genr is attached to a loop that
   is saved across function calls -- then the pointer is
   reset to NULL on each call to the function (but not
   on each iteration of the loop itself).
*/

void reattach_series(Node * n, Parser * p)
{
    /* FIXME relax the "1" below? */
    if (true || n.v.xvec == null || get_loop_renaming()) {
        /* do a full reset */
        n.vnum = current_series_index(p.dset, n.vname);
        if (n.vnum < 0) {
            gretl_errmsg_sprintf("'%s': not a series".ptr, n.vname);
            p.err = GretlError.data;
        } else {
            n.v.xvec = p.dset.Z[n.vnum];
        }
    } else if (n.vnum >= 0 && n.vnum < p.dset.v) {
        /* handles the case of the dataset moving due to
           sub-sampling, provided that no series are deleted
           (the name -> ID mapping remains unchanged)
        */
        n.v.xvec = p.dset.Z[n.vnum];
    }
}

/* The message is assembled with D's std.format rather than a glib
   GString (a small departure from the C code, same text). */

void reattach_data_error(Node * n, Parser * p)
{
    import std.format : format;
    import std.string : fromStringz, toStringz;

    string msg = format("'%s': expected %s", fromStringz(n.vname),
                        fromStringz(getsymb(n.t)));

    if (n.uv == null) {
        msg ~= " but name look-up failed";
        p.err = GretlError.data;
    } else {
        string s = fromStringz(gretl_type_get_name(n.uv.type)).idup;

        if (n.uv.ptr == null) {
            msg ~= format(", found %s with NULL data", s);
            p.err = GretlError.data;
        } else {
            msg ~= format(" but found %s", s);
            p.err = GretlError.types;
        }
    }

    gretl_errmsg_set(toStringz(msg));
}

/* Deep debugging to detect any cases where we've failed to NULL-ify
   uservars (n->uv) that may have been relocated in memory between
   visits to geneval. (The UVDEBUG variant is not ported.)
*/

void node_reattach_data(Node * n, Parser * p)
{
    if (n.t == GenSym.series) {
        reattach_series(n, p);
    } else {
        GretlType type = cast(GretlType) 0;
        void * data = null;

        if (n.uv == null || user_qsorting || (n.t == GenSym.list && gretl_iterating())) {
            n.uv = get_user_var_by_name(n.vname);
        }

        if (n.uv != null) {
            data = n.uv.ptr;
            type = n.uv.type;
        }

        if (data == null) {
            reattach_data_error(n, p);
        } else if (uscalar_node(n)) {
            if (type == GretlType.double_) {
                n.v.xval = *cast(double *) data;
            } else if (type == GretlType.matrix) {
                /* allow type-mutation */
                n.t = GenSym.mat;
                n.v.m = cast(GretlMatrixRaw *) data;
            } else {
                reattach_data_error(n, p);
            }
        } else if (n.t == GenSym.mat) {
            if (type == GretlType.matrix) {
                n.v.m = cast(GretlMatrixRaw *) data;
            } else if (type == GretlType.double_) {
                /* allow type-mutation */
                n.t = GenSym.num;
                n.v.xval = *cast(double *) data;
            } else {
                reattach_data_error(n, p);
            }
        } else if (n.t == GenSym.list && type == GretlType.list) {
            n.v.ivec = cast(int *) data;
        } else if (n.t == GenSym.bundle && type == GretlType.bundle) {
            n.v.b = cast(GretlBundle *) data;
        } else if (n.t == GenSym.str && type == GretlType.string_) {
            n.v.str = cast(char *) data;
        } else if (n.t == GenSym.array && type == GretlType.array) {
            n.v.a = cast(GretlArray *) data;
        } else if (n.t == GenSym.osl) {
            n.v.ptr = data;
        } else {
            reattach_data_error(n, p);
        }
    }
}

void node_type_error(int ntype, int argnum, int goodt,
                     Node * bad, Parser * p)
{
    const(char) * nstr;

    if (ntype == 0) {
        p.err = GretlError.types;
        return;
    } else if (bad != null && bad.t == GenSym.undef) {
        p.err = GretlError.unkvar;
        return;
    }

    parser_ensure_error_buffer(p);

    if (ntype == GenSym.lag) {
        nstr = (goodt == GenSym.num)? "lag order" : "lag variable";
    } else {
        nstr = getsymb(ntype);
    }

    if (goodt == 0 && bad != null) {
        if (null_node(bad)) {
            pprintf(p.prn, _("%s: insufficient arguments").ptr, nstr);
        } else if (argnum <= 0) {
            pprintf(p.prn, _("%s: invalid argument type %s").ptr,
                    nstr, typestr(bad.t));
        } else {
            pprintf(p.prn, _("%s: argument %d: invalid type %s").ptr,
                    nstr, argnum, typestr(bad.t));
        }
        pputc(p.prn, '\n');
        p.err = GretlError.types;
        return;
    }

    if (ntype < GenSym.op_max) {
        if (argnum <= 0) {
            pprintf(p.prn, _("%s: operand should be %s").ptr,
                    nstr, typestr(goodt));
        } else {
            pprintf(p.prn, _("%s: operand %d should be %s").ptr,
                    nstr, argnum, typestr(goodt));
        }
    } else {
        if (argnum <= 0) {
            pprintf(p.prn, _("%s: argument should be %s").ptr,
                    nstr, typestr(goodt));
        } else {
            pprintf(p.prn, _("%s: argument %d should be %s").ptr,
                    nstr, argnum, typestr(goodt));
        }
    }

    if (bad != null) {
        pprintf(p.prn, _(", is %s").ptr, typestr(bad.t));
        if (bad.t == GenSym.str) {
            pprintf(p.prn, " (%s)".ptr, bad.v.str);
        }
    }
    pputc(p.prn, '\n');

    if (!strcmp(nstr, "&")) {
        pputs(p.prn, "(for logical AND, please use \"&&\")\n".ptr);
    } else if (!strcmp(nstr, "|")) {
        pputs(p.prn, "(for logical OR, please use \"||\")\n".ptr);
    }

    p.err = GretlError.types;
}

int node_is_true(Node * n, Parser * p)
{
    double x = node_get_scalar(n, p);

    return !na(x) && x != 0.0;
}

int node_is_false(Node * n, Parser * p)
{
    return (node_get_scalar(n, p) == 0.0);
}

/* clean up a node holding evaluated sub-nodes of FARGS */

void destroy_multi(Node * m)
{
    free(m.v.bn.n);
    free(m);
}

Node * evaluate_multi_args(Node * nn, Parser * p)
{
    Node * ret = bncopy(nn, &p.err);
    int i;

    for (i=0; i<nn.v.bn.n_nodes && !p.err; i++) {
        ret.v.bn.n[i] = eval(nn.v.bn.n[i], p);
    }

    return ret;
}

bool ok_nullarg_sym(int s)
{
    return (s == GenSym.fargs || s > GenFunc.fp_max || s == GenSym.u_addr);
}

/* determine if an EMPTY node is acceptable in context */

int empty_ok(Node * n, Parser * p)
{
    if (n.vname == null) {
        /* regular EMPTY node, OK */
        return 1;
    } else if (n.parent != null && ok_nullarg_sym(n.parent.t)) {
        /* null argument pass-through, OK */
        return 1;
    } else {
        p.err = GretlError.types;
        gretl_errmsg_sprintf("'%s': no such object".ptr, n.vname);
        return 0;
    }
}
