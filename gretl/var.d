/* Port of lib/src/var.c: vector autoregressions (VAR) and vector error-correction models (VECM): estimation, lag
   selection, impulse responses, forecast-error variance decompositions, tests, saving to XML and bundles. The
   headers var.h, johansen.h and vartest.h are varh.d, johansenh.d and vartesth.d. */

module gretl.var;


import core.stdc.ctype : isdigit, isspace;
import core.stdc.errno : errno, ERANGE;
import core.stdc.math : ceil, exp, fabs, floor, fmod, isinf, isnan, log, log10, pow, sqrt;
import core.stdc.stdio : fputs, fprintf, printf, puts, sprintf, snprintf, sscanf, stderr;
import core.stdc.stdlib : abs, atoi, calloc, free, getenv, malloc, realloc;
import core.stdc.string : memcpy, memset, strcmp, strcpy, strlen, strchr, strstr, strcat, strncat, strncmp, strspn, strcspn, memmove, strrchr;

import std.algorithm : max, min;

import gretl.gretl_errorsh : GretlError;
import gretl.gretl_matrix;
import gretl.gretl_matrixh;
import gretl.johansenh : JohansenCode, JohansenInfo, jcode, jrank, auto_restr;
import gretl.libgretlh;
import gretl.missingh : NADBL, na;
import gretl.gretl_intl : gretl_pop_c_numeric_locale, gretl_push_c_numeric_locale;
import gretl.gretl_restricth : GretlRestriction;
import gretl.varh;
import gretl.xmlh;
import gretl.vartesth : N_IVALS, VAR_S_DFCORR;

private alias XUC = const(xmlChar) *;

import gretl.calendar : weekday_from_date;
import gretl.compareh : ModelAuxCode;
import gretl.dataio : get_subperiod, ntolabel;
import gretl.dataseth : dataset_is_incomplete_daily, dated_daily_data;
import gretl.describe : gretl_mean, gretl_stddev;
import gretl.estimate : makevcv;
import gretl.genmainh : ModelDataIndex;
import gretl.genparseh : ConstId;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_errors : gretl_errmsg_set, gretl_errmsg_sprintf, gretl_error_clear;
import gretl.gretl_list : gretl_list_copy, gretl_list_duplicates, gretl_list_has_separator, gretl_list_new, gretl_list_purge_const, gretl_list_sort, gretl_list_split_on_separator, gretl_list_to_vector, in_gretl_list;
import gretl.gretl_listh : LISTSEP;
import gretl.gretl_model : clear_model, gretl_model_allocate_storage, gretl_model_array_destroy, gretl_model_array_new, gretl_model_from_XML, gretl_model_get_vcv_element, gretl_model_serialize, gretl_model_smpl_init;
import gretl.gretl_prn : Prn, pprintf, pputc, pputs;
import gretl.gretl_utils : copyvec;
import gretl.gretl_xml : gretl_matrix_serialize, gretl_xml_get_double_array, gretl_xml_get_list, gretl_xml_get_matrix, gretl_xml_get_prop_as_double, gretl_xml_get_prop_as_int, gretl_xml_get_prop_as_string, gretl_xml_put_double, gretl_xml_put_double_array, gretl_xml_put_int, gretl_xml_put_tagged_list;
import gretl.interacth : MAXSAVENAME;
import gretl.libset : libset_get_int;
import gretl.libseth : SetKey;
import gretl.matrix_extra : gretl_matrix_data_subset, gretl_matrix_print_to_prn;
import gretl.matrix_extrah : MMissingCode;
import gretl.objstackh : SavedObjectFlags;
import gretl.plugins : get_plugin_function;
import gretl.pvalues : chisq_cdf_comp;
import gretl.strutils : dot_atof, gretl_strdup, gretl_xml_encode, strings_array_dup, strings_array_new, strings_array_new_with_length;
import gretl.vartest : VAR_LR_lag_test, VAR_do_lagsel, VAR_portmanteau_test, VAR_wald_omit_tests;
import gretl.varprint : gretl_VAR_print;
import gretl.irfboot : VAR_coeff_matrix_from_VECM, irf_bootstrap;
import gretl.gretl_restrict : gretl_restriction_get_options, rset_add_results, rset_get_R_matrix, rset_get_q_matrix, rset_record_LR_result;
import gretl.system : multi_eqn_wald_test;
import gretl.graphing : gretl_VAR_plot_impulse_response;
import gretl.gretl_array : gretl_array_destroy;
import gretl.gretl_bundle : gretl_bundle_donate_data, gretl_bundle_get_bundle, gretl_bundle_get_int, gretl_bundle_get_list, gretl_bundle_get_matrix, gretl_bundle_new, gretl_bundle_set_int, gretl_bundle_set_list, gretl_bundle_set_matrix, gretl_bundle_set_note, gretl_bundle_set_scalar, gretl_bundle_set_string;
import gretl.gretl_commands : gretl_command_word;

/* unported dependencies: per-file declarations from the 2026c headers */
extern (C) {
}


enum VarTask {
    var_estimate = 1,
    lagsel,
    vecm_estimate,
    vecm_ctest
}

int VAR_add_models(GretlVar * var, const(DatasetRaw) * dset)
{
    int n = var.neqns;
    int i, err = 0;

    if (var.models != null) {
        for (i=0; i<n; i++) {
            clear_model(var.models[i]);
        }
    } else {
        var.models = gretl_model_array_new(n);
        if (var.models == null) {
            err = GretlError.alloc;
        }
    }

    if (var.models != null) {
        for (i=0; i<n; i++) {
            gretl_model_smpl_init(var.models[i], dset);
        }
    }

    return err;
}

/* Here we allocate storage for a full-size companion matrix and fill in
   the "boilerplate" rows at the bottom. The initial (VAR coefficient)
   rows are left as zeros.
*/

int VAR_allocate_companion_matrix(GretlVar * var)
{
    int n, dim, err = 0;

    if (var.A != null) {
    /* the job is already done */
        return 0;
    }

    n = var.neqns;
    dim = n * levels_order(var);
    var.A = gretl_zero_matrix_new(dim, dim);

    if (var.A == null) {
        err = GretlError.alloc;
    } else {
    int j, jmax = dim - n;

    for (j=0; j<jmax; j++) {
        gretl_matrix_set(var.A, j+n, j, 1);
    }
    }

    return err;
}

/* Given an m x n matrix @A, with m < n, construct the full n x n
   companion matrix.
*/

GretlMatrixRaw * companionize(const(GretlMatrixRaw) * A, int * err)
{
    GretlMatrixRaw * ret = null;
    int m = A.rows, n = A.cols;
    int d = n - m;

    if (d < 0) {
    *err = GretlError.invarg;
    return null;
    }

    if (d == 0) {
    /* shouldn't get here, but just in case */
    ret = gretl_matrix_copy(A);
    } else {
    double aij;
    int i, j;

    ret = gretl_zero_matrix_new(n, n);
    if (ret != null) {
        for (j=0; j<n; j++) {
        for (i=0; i<m; i++) {
            aij = gretl_matrix_get(A, i, j);
            gretl_matrix_set(ret, i, j, aij);
        }
        if (j < d) {
            gretl_matrix_set(ret, j+m, j, 1);
        }
        }
    }
    }

    if (ret == null) {
    *err = GretlError.alloc;
    }

    return ret;
}

/* Given the full (square) companion matrix, return a version without
   the trailing rows holding (I ~ 0). The inverse operation to
   companionize() above.
*/

GretlMatrixRaw * decompanionize(const(GretlMatrixRaw) * A, int neqns, GretlMatrixMod mod)
{
    GretlMatrixRaw * ret;
    int np = A.cols;

    if (mod == GretlMatrixMod.transpose) {
    ret = gretl_matrix_alloc(np, neqns);
    } else {
    ret = gretl_matrix_alloc(neqns, np);
    }

    if (ret != null) {
    double aij;
    int i, j;

    for (j=0; j<np; j++) {
        for (i=0; i<neqns; i++) {
        aij = gretl_matrix_get(A, i, j);
        if (mod == GretlMatrixMod.transpose) {
            gretl_matrix_set(ret, j, i, aij);
        } else {
            gretl_matrix_set(ret, i, j, aij);
        }
        }
    }
    }

    return ret;
}

int VAR_allocate_cholesky_matrix(GretlVar * var)
{
    int n = var.neqns;
    int err = 0;

    if (var.C != null) {
    /* the job is already done */
        return 0;
    }

    var.C = gretl_zero_matrix_new(n, n);
    if (var.C == null) {
        err = GretlError.alloc;
    }

    return err;
}

int VAR_allocate_residuals_matrix(GretlVar * var)
{
    int err = 0;

    if (var.E != null) {
    /* the job is already done */
        return 0;
    }

    var.E = gretl_zero_matrix_new(var.T, var.neqns);
    if (var.E == null) {
        err = GretlError.alloc;
    }

    return err;
}

int VAR_add_XTX_matrix(GretlVar * var)
{
    int k = var.X.cols;
    int err = 0;

    var.XTX = gretl_zero_matrix_new(k, k);

    if (var.XTX == null) {
        err = GretlError.alloc;
    } else {
        gretl_matrix_multiply_mod(var.X, GretlMatrixMod.transpose,
                                  var.X, GretlMatrixMod.none,
                                  var.XTX, GretlMatrixMod.none);
        err = gretl_invert_matrix(var.XTX);
    }

    return err;
}

int n_restricted_terms(GretlVar * v)
{
    int n = 0;

    if (v.jinfo != null &&
        (v.jinfo.code == JohansenCode.rest_const ||
         v.jinfo.code == JohansenCode.rest_trend)) {
        n++;
    }

    if (v.rlist != null) {
        n += v.rlist[0];
    }

    return n;
}

void gretl_VAR_clear(GretlVar * var)
{
    var.ci = 0;
    var.refcount = 0;
    var.err = 0;
    var.neqns = var.order = 0;
    var.t1 = var.t2 = var.T = var.df = 0;
    var.ifc = var.ncoeff = 0;
    var.detflags = 0;
    var.robust = 0;
    var.LBs = 0;
    var.xcols = 0;

    var.lags = null;
    var.ylist = null;
    var.xlist = null;
    var.rlist = null;

    var.Y = null;
    var.X = null;
    var.B = null;
    var.XTX = null;
    var.A = null;
    var.L = null;
    var.E = null;
    var.C = null;
    var.S = null;
    var.F = null;
    var.V = null;
    var.ord = null;

    var.models = null;
    var.Fvals = null;
    var.Ivals = null;

    var.ll = var.ldet = var.LR = NADBL;
    var.AIC = var.BIC = var.HQC = NADBL;
    var.LB = NADBL;

    var.jinfo = null;
    var.name = null;
}

int VAR_add_regular_seasonals(GretlVar * v, const(DatasetRaw) * dset, int * pk)
{
    int pd1 = dset.pd - 1;
    int i, t, k = *pk;
    double s0, s1;
    int err = 0;

    if (v.ci == GretlCmdIndex.vecm) {
    s1 = 1 - 1.0 / dset.pd;
    s0 = s1 - 1;
    } else {
    s1 = 1;
    s0 = 0;
    }

    if (dataset_is_incomplete_daily(dset)) {
    /* we have to read the date of each observation,
       since there are holes in the time series
    */
    int dow;

    for (t=0; t<v.T; t++) {
        dow = weekday_from_date(dset.S[t + v.t1]);
            dow = dow == 7 ? 0 : dow; /* compatibility */
        for (i=0; i<pd1; i++) {
        gretl_matrix_set(v.X, t, k+i, (dow == i+1)? s1 : s0);
        }
    }
    } else {
    /* we can proceed mostly algorithmically */
    int per, offset = 0;

    if (dated_daily_data(dset)) {
        char[OBSLEN] obs;

        ntolabel(obs.ptr, v.t1, dset);
        per = weekday_from_date(obs.ptr);
            per = per == 7 ? 0 : per; /* compatibility */
            offset = 1;
    } else {
        per = get_subperiod(v.t1, dset, null);
    }
    for (t=0; t<v.T; t++) {
        for (i=0; i<pd1; i++) {
        gretl_matrix_set(v.X, t, k+i, (per == i+offset)? s1 : s0);
        }
        per = (per < pd1)? per + 1 : 0;
    }
    }

    *pk += pd1;

    return err;
}

bool lag_wanted(GretlVar * v, int i)
{
    return (v.lags == null || in_gretl_list(v.lags, i));
}

/* Construct the common X matrix (composed of lags of the core variables
   plus other terms if applicable). This is in common between VARs and
   VECMs.
*/

void VAR_fill_X(GretlVar * v, int p, const(DatasetRaw) * dset)
{
    const(double) * x;
    int i, j, s, t, vi;
    int k = 0; /* X column index */

    /* const first */
    if (v.detflags & DetFlags.const_) {
        s = 0;
        for (t=v.t1; t<=v.t2; t++) {
            gretl_matrix_set(v.X, s++, k, 1.0);
        }
        k++;
    }

    /* add lagged Ys */
    for (i=0; i<v.neqns; i++) {
        vi = v.ylist[i+1];
        for (j=1; j<=p; j++) {
            if (!lag_wanted(v, j)) {
                continue;
            }
            x = dset.Z[vi];
            s = 0;
            for (t=v.t1; t<=v.t2; t++) {
                if (v.ci == GretlCmdIndex.vecm) {
                    gretl_matrix_set(v.X, s++, k, x[t-j] - x[t-j-1]);
                } else {
                    gretl_matrix_set(v.X, s++, k, x[t-j]);
                }
            }
            k++;
        }
    }

    /* add any exogenous vars */
    if (v.xlist != null) {
        for (i=1; i<=v.xlist[0]; i++) {
            vi = v.xlist[i];
            s = 0;
            for (t=v.t1; t<=v.t2; t++) {
                gretl_matrix_set(v.X, s++, k, dset.Z[vi][t]);
            }
            k++;
        }
    }

    /* add seasonals if wanted */
    if (v.detflags & DetFlags.seas) {
    VAR_add_regular_seasonals(v, dset, &k);
    }

    /* add trend if wanted */
    if (v.detflags & DetFlags.trend) {
        s = 0;
        for (t=v.t1; t<=v.t2; t++) {
            gretl_matrix_set(v.X, s++, k, cast(double) (t + 1));
        }
        k++;
    }

    if (v.X != null) {
        gretl_matrix_set_t1(v.X, v.t1);
        gretl_matrix_set_t2(v.X, v.t2);
    }

}

/* Construct the matrix of dependent variables for a plain VAR. */

void VAR_fill_Y(GretlVar * v, const(DatasetRaw) * dset)
{
    int i, vi, s, t;

    for (i=0; i<v.neqns; i++) {
        vi = v.ylist[i+1];
        s = 0;
        for (t=v.t1; t<=v.t2; t++) {
            gretl_matrix_set(v.Y, s++, i, dset.Z[vi][t]);
        }
    }

    if (v.Y != null) {
        gretl_matrix_set_t1(v.Y, v.t1);
        gretl_matrix_set_t2(v.Y, v.t2);
    }

}

/* Construct the combined matrix of dependent variables for a VECM. */

void VECM_fill_Y(GretlVar * v, const(DatasetRaw) * dset, GretlMatrixRaw * Y)
{
    const(double) * yi;
    int i, vi, s, t;
    int k = 0;

    for (i=0; i<v.neqns; i++) {
        vi = v.ylist[i+1];
        yi = dset.Z[vi];
        k = i + v.neqns;
        s = 0;
        for (t=v.t1; t<=v.t2; t++) {
            gretl_matrix_set(Y, s, i, yi[t] - yi[t-1]);
            gretl_matrix_set(Y, s, k, yi[t-1]);
            s++;
        }
    }

    if (auto_restr(v)) {
        int trend = (v.jinfo.code == JohansenCode.rest_trend);

        k++;
        for (t=0; t<v.T; t++) {
            gretl_matrix_set(Y, t, k, trend ? (v.t1 + t) : 1);
        }
    }

    if (v.rlist != null) {
        /* There's room for debate here (?) over whether we should
           use the current value or the first lag of restricted
           exogenous variables. But using the current value agrees
           with Ox. See also VAR_set_sample().
        */
        for (i=0; i<v.rlist[0]; i++) {
            vi = v.rlist[i+1];
            k++;
            s = 0;
            for (t=v.t1; t<=v.t2; t++) {
                gretl_matrix_set(Y, s++, k, dset.Z[vi][t]); /* was t-1 */
            }
        }
    }

}

/* Check that we have no series duplicated in a VAR specification */

int VAR_check_lists(GretlVar * v, const(DatasetRaw) * dset)
{
    int *[3] LL = [v.ylist, v.xlist, v.rlist];
    int dup, err = 0;
    int i, vi;

    /* check for any series included in more than one
       of v->ylist, v->xlist, v->rlist
    */
    if (LL[1] || LL[2]) {
    for (i=1; i<=LL[0][0] && !err; i++) {
        vi = LL[0][i];
        if (LL[1] && in_gretl_list(LL[1], vi)) {
        err = GretlError.invarg;
        } else if (LL[2] && in_gretl_list(LL[2], vi)) {
        err = GretlError.invarg;
        }
    }
    if (!err && LL[1] && LL[2]) {
        for (i=1; i<=LL[1][0] && !err; i++) {
        vi = LL[1][i];
        if (in_gretl_list(LL[2], vi)) {
            err = GretlError.invarg;
        }
        }
    }
    if (err && v.ci == GretlCmdIndex.vecm) {
        gretl_errmsg_set(_("ylist, xlist and rxlist must be disjoint").ptr);
    } else if (err) {
        gretl_errmsg_set(_("ylist and xlist must be disjoint").ptr);
    }
    }

    for (i=0; i<3 && !err; i++) {
    /* check for any within-list duplicates */
    if (LL[i]) {
        dup = gretl_list_duplicates(LL[i], cast(GretlCmdIndex) 0);
        if (dup >= 0) {
        err = GretlError.invarg;
        gretl_errmsg_sprintf(_("duplicate term %s in %s").ptr,
                     dset.varname[dup],
                     i == 0 ? "ylist".ptr :
                     i == 1 ? "xlist".ptr : "rxlist".ptr);
        }
    }
    }

    return err;
}

/* Split the user-supplied list, if need be, and construct the lists of
   endogenous and (possibly) exogenous vars.  Note that deterministic
   terms are handled separately, via option flags.
*/

int VAR_make_lists(GretlVar * v, const(int) * list, const(DatasetRaw) * dset)
{
    int err = 0;

    if (gretl_list_has_separator(list)) {
        err = gretl_list_split_on_separator(list, &v.ylist, &v.xlist);
    } else {
        v.ylist = gretl_list_copy(list);
        if (v.ylist == null) {
            err = GretlError.alloc;
        }
    }

    if (!err && (v.ylist == null || v.ylist[0] < 1)) {
        /* first test for at least 1 endog var */
        err = GretlError.args;
    }

    if (!err && v.ci == GretlCmdIndex.vecm) {
        if (v.xlist != null && gretl_list_has_separator(v.xlist)) {
            int * tmp = null;

            err = gretl_list_split_on_separator(v.xlist, &tmp, &v.rlist);
            if (!err) {
                free(v.xlist);
                v.xlist = tmp;
            }
        }
    }

    if (!err) {
        gretl_list_purge_const(v.ylist, dset);
        if (v.xlist != null) {
            gretl_list_purge_const(v.xlist, dset);
        }
        if (v.rlist != null) {
            gretl_list_purge_const(v.rlist, dset);
        }
    }

    if (!err && (v.ylist == null || v.ylist[0] < 1)) {
        /* re-test after (possibly) losing const */
        err = GretlError.args;
    }

    if (!err) {
    err = VAR_check_lists(v, dset);
    }

    return err;
}

/* Starting from the given sample range, construct the feasible
   estimation range for a VAR or VECM.  Flag an error if there are
   missing values within the sample period.  Also see
   system_adjust_t1t2() in lib/src/system.c (perhaps the two could be
   consolidated at some point).
*/

int VAR_set_sample(GretlVar * v, const(DatasetRaw) * dset)
{
    int i, vi, t, err = 0;

    /* advance t1 if needed */

    for (t=dset.t1; t<=dset.t2; t++) {
        int miss = 0, p, s = t - levels_order(v);

        for (i=1; i<=v.ylist[0] && !miss; i++) {
            vi = v.ylist[i];
            if (na(dset.Z[vi][t]) || s < 0) {
                v.t1 += 1;
                miss = 1;
            }
            for (p=s; p<t && !miss; p++) {
                if (na(dset.Z[vi][p])) {
                    v.t1 += 1;
                    miss = 1;
                }
            }
        }
        if (v.xlist != null && !miss) {
            for (i=1; i<=v.xlist[0] && !miss; i++) {
                vi = v.xlist[i];
                if (na(dset.Z[vi][t])) {
                    v.t1 += 1;
                    miss = 1;
                }
            }
        }

        /* In estimating a VECM we need the level (or first lag?)
           of each restricted exogenous var to serve as dependent
           variable in one of the initial OLS equations. See also
           VECM_fill_Y().
        */
        if (v.rlist != null && !miss) {
            for (i=1; i<=v.rlist[0] && !miss; i++) {
                vi = v.rlist[i];
                if (na(dset.Z[vi][t])) {
                    v.t1 += 1;
                    miss = 1;
                }
            }
        }
        if (!miss) {
            break;
        }
    }

    /* retard t2 if needed */

    for (t=dset.t2; t>=v.t1; t--) {
        int miss = 0;

        for (i=1; i<=v.ylist[0] && !miss; i++) {
            vi = v.ylist[i];
            if (na(dset.Z[vi][t])) {
                v.t2 -= 1;
                miss = 1;
            }
        }
        if (v.xlist != null && !miss) {
            for (i=1; i<=v.xlist[0] && !miss; i++) {
                vi = v.xlist[i];
                if (na(dset.Z[vi][t])) {
                    v.t2 -= 1;
                    miss = 1;
                }
            }
        }
        if (v.rlist != null && !miss) {
            for (i=1; i<=v.rlist[0] && !miss; i++) {
                vi = v.rlist[i];
                if (na(dset.Z[vi][t])) {
                    v.t2 -= 1;
                    miss = 1;
                }
            }
        }
        if (!miss) {
            break;
        }
    }

    /* reject sample in case of internal missing values */

    for (t=v.t1; t<=v.t2 && !err; t++) {
        for (i=1; i<=v.ylist[0] && !err; i++) {
            vi = v.ylist[i];
            if (na(dset.Z[vi][t])) {
                err = GretlError.missdata;
            }
        }
        if (v.xlist != null && !err) {
            for (i=1; i<=v.xlist[0] && !err; i++) {
                vi = v.xlist[i];
                if (na(dset.Z[vi][t])) {
                    err = GretlError.missdata;
                }
            }
        }
        if (v.rlist != null && !err) {
            for (i=1; i<=v.rlist[0] && !err; i++) {
                vi = v.rlist[i];
                if (na(dset.Z[vi][t])) {
                    err = GretlError.missdata;
                }
            }
        }
    }

    return err;
}

/* Account for deterministic terms and check for non-negative degrees of
   freedom.
*/

int VAR_check_df_etc(GretlVar * v, const(DatasetRaw) * dset, Gretlopt opt)
{
    int nl, err = 0;

    v.T = v.t2 - v.t1 + 1;

    nl = var_n_lags(v);
    v.ncoeff = nl * v.neqns;

    if (v.xlist != null) {
        /* user-specified exogenous vars */
        v.ncoeff += v.xlist[0];
    }

    if (v.ci == GretlCmdIndex.vecm) {
        if (!(opt & Gretlopt.OPT_N) && !(opt & Gretlopt.OPT_R)) {
            v.detflags |= DetFlags.const_;
            v.ncoeff += 1;
            v.ifc = 1;
        }
    } else if (!(opt & Gretlopt.OPT_N)) {
        v.detflags |= DetFlags.const_;
        v.ncoeff += 1;
        v.ifc = 1;
    }

    if ((opt & Gretlopt.OPT_D) && dset.pd != 1) {
        v.detflags |= DetFlags.seas;
        v.ncoeff += dset.pd - 1;
    }

    if (opt & Gretlopt.OPT_T) {
        v.detflags |= DetFlags.trend;
        v.ncoeff += 1;
    }

    v.df = v.T - v.ncoeff;

    if (v.df < 0) {
        err = GretlError.df;
    }

    return err;
}

int VAR_add_basic_matrices(GretlVar * v)
{
    int err = 0;

    v.Y = gretl_matrix_alloc(v.T, v.neqns);
    v.E = gretl_matrix_alloc(v.T, v.neqns);

    if (v.Y == null || v.E == null) {
        err = GretlError.alloc;
    }

    if (!err && v.ncoeff > 0) {
        v.X = gretl_matrix_alloc(v.T, v.ncoeff);
        v.B = gretl_matrix_alloc(v.ncoeff, v.neqns);
        if (v.X == null || v.B == null) {
            err = GretlError.alloc;
        } else {
            /* record the number of columns of X */
            v.xcols = v.X.cols;
        }
    }

    return err;
}

void set_to_NA(double * x, int n)
{
    int i;

    for (i=0; i<n; i++) {
        x[i] = NADBL;
    }
}

/* Main function for constructing a new VAR struct, which may be used
   for estimating a VAR, a VECM, or various auxiliary tasks.
*/

GretlVar * gretl_VAR_new(int code, int order, int rank, const(int) * lags, const(int) * list, const(DatasetRaw) * dset, Gretlopt opt, int * errp)
{
    GretlVar * var;
    int ci, err = 0;

    ci = (code >= VarTask.vecm_estimate)? GretlCmdIndex.vecm : GretlCmdIndex.var;

    if (order < 0) {
        gretl_errmsg_sprintf(_("Invalid lag order %d").ptr, order);
        *errp = GretlError.invarg;
        return null;
    }

    var = cast(typeof(var)) malloc((*var).sizeof);
    if (var == null) {
        *errp = GretlError.alloc;
        return null;
    }

    gretl_VAR_clear(var);

    var.ci = ci;
    var.order = order;
    var.lags = gretl_list_copy(lags);

    if (ci == GretlCmdIndex.var) {
        if (opt & Gretlopt.OPT_H) {
            var.robust = VarRobust.hac;
        } else if (opt & Gretlopt.OPT_R) {
            var.robust = VarRobust.hc;
        }
    }

    err = VAR_make_lists(var, list, dset);

    if (!err && rank > var.ylist[0]) {
        gretl_errmsg_sprintf(_("vecm: rank %d is out of bounds").ptr, rank);
        err = GretlError.data;
    }

    if (!err) {
        var.neqns = var.ylist[0];
        var.t1 = dset.t1;
        var.t2 = dset.t2;
    }

    /* FIXME below: some of these allocations are redundant
       for some uses of the VAR struct */

    if (!err) {
        err = VAR_set_sample(var, dset);
    }

    if (!err) {
        err = VAR_check_df_etc(var, dset, opt);
    }

    if (!err) {
        err = VAR_add_basic_matrices(var);
    }

    if (!err && var.ci == GretlCmdIndex.var) {
        err = VAR_allocate_companion_matrix(var);
    }

    if (!err && var.ci == GretlCmdIndex.var) {
        err = VAR_allocate_cholesky_matrix(var);
    }

    if (!err && code != VarTask.lagsel) {
        err = VAR_add_models(var, dset);
    }

    if (!err && code == VarTask.var_estimate) {
        int m = var.neqns * var.neqns + var.neqns;

        var.Fvals = cast(typeof(var.Fvals)) malloc(m * (*var.Fvals).sizeof);
        if (var.Fvals == null) {
            err = 1;
        } else {
            set_to_NA(var.Fvals, m);
        }
    }

    if (!err && code == VarTask.var_estimate) {
        var.Ivals = cast(typeof(var.Ivals)) malloc(N_IVALS * (*var.Ivals).sizeof);
        if (var.Ivals == null) {
            err = 1;
        } else {
            set_to_NA(var.Ivals, N_IVALS);
        }
    }

    if (!err && var.ci == GretlCmdIndex.vecm) {
        var.jinfo = johansen_info_new(var, rank, opt);
        if (var.jinfo == null) {
            err = GretlError.alloc;
        } else if (var.detflags & DetFlags.seas) {
            var.jinfo.seasonals = dset.pd - 1;
        }
    }

    if (!err) {
        if (var.ci == GretlCmdIndex.var) {
            /* this is deferred for a VECM */
            VAR_fill_Y(var, dset);
        }
        VAR_fill_X(var, order, dset);
    }

    if (err) {
        gretl_VAR_free(var);
        var = null;
    }

    *errp = err;

    return var;
}

void johansen_info_free(GretlVar * var)
{
    JohansenInfo *jv = var.jinfo;

    gretl_matrix_free(jv.R0);
    gretl_matrix_free(jv.R1);

    gretl_matrix_free(jv.S00);
    gretl_matrix_free(jv.S11);
    gretl_matrix_free(jv.S01);

    gretl_matrix_free(jv.evals);
    gretl_matrix_free(jv.Beta);
    gretl_matrix_free(jv.Alpha);
    gretl_matrix_free(jv.Gamma);
    gretl_matrix_free(jv.JC);
    gretl_matrix_free(jv.JVC);
    gretl_matrix_free(jv.Cse);
    gretl_matrix_free(jv.Bvar);
    gretl_matrix_free(jv.Bse);
    gretl_matrix_free(jv.Ase);
    gretl_matrix_free(jv.R);
    gretl_matrix_free(jv.q);
    gretl_matrix_free(jv.Ra);
    gretl_matrix_free(jv.qa);

    gretl_matrix_free(jv.YY);
    gretl_matrix_free(jv.RR);
    gretl_matrix_free(jv.BB);
    if (jv.G != null) {
        gretl_array_destroy(jv.G);
    }

    free(jv);
}

void gretl_VAR_free(GretlVar * var)
{
    if (var == null) {
    return;
    }

    var.refcount -= 1;
    if (var.refcount > 0) {
        return;
    }

    free(var.lags);
    free(var.ylist);
    free(var.xlist);
    free(var.rlist);

    gretl_matrix_free(var.Y);
    gretl_matrix_free(var.X);
    gretl_matrix_free(var.B);
    gretl_matrix_free(var.XTX);

    gretl_matrix_free(var.A);
    gretl_matrix_free(var.L);
    gretl_matrix_free(var.E);
    gretl_matrix_free(var.C);
    gretl_matrix_free(var.S);
    gretl_matrix_free(var.F);
    gretl_matrix_free(var.V);
    gretl_matrix_free(var.ord);

    free(var.Fvals);
    free(var.Ivals);
    free(var.name);

    if (var.models != null) {
        gretl_model_array_destroy(var.models, var.neqns);
    }

    if (var.jinfo != null) {
        johansen_info_free(var);
    }

    free(var);

}

int VAR_add_fcast_variance(GretlVar * var, GretlMatrixRaw * F, int n_static)
{
    GretlMatrixRaw * se = null;
    double ftj, vti;
    int k, n = var.neqns;
    int i, j, s;
    int err = 0;

    if (n_static < 0) {
        fprintf(stderr, "n_static = %d\n".ptr, n_static);
        n_static = 0;
    }

    k = F.rows - n_static;

    if (k > 0) {
        se = gretl_VAR_get_fcast_se(var, k);
        if (se == null) {
            err = GretlError.alloc;
            goto bailout;
        }
    }

    for (j=0; j<n; j++) {
        for (s=0; s<F.rows; s++) {
            ftj = gretl_matrix_get(F, s, j);
            if (na(ftj)) {
                gretl_matrix_set(F, s, n + j, NADBL);
            } else {
                i = s - n_static;
                if (i < 0) {
                    vti = sqrt(gretl_matrix_get(var.S, j, j));
                } else {
                    vti = gretl_matrix_get(se, i, j);
                }
                gretl_matrix_set(F, s, n + j, vti);
            }
        }
    }

    if (se != null) {
        gretl_matrix_free(se);
    }

 bailout:

    if (err) {
        for (i=0; i<n; i++) {
            for (s=0; s<F.rows; s++) {
                gretl_matrix_set(F, s, n + i, NADBL);
            }
        }
    }

    return err;
}

/* determine start of dynamic portion of forecast */

int VAR_get_tdyn(GretlVar * var, int t1, int t2, Gretlopt opt)
{
    int td;

    if (opt & Gretlopt.OPT_D) {
        /* force dynamic */
        td = t1;
    } else if (opt & Gretlopt.OPT_S) {
        /* force static */
        td = t2 + 1;
    } else {
        /* out of sample */
        td = var.t2 + 1;
    }

    return td;
}

int VAR_add_forecast(GretlVar * var, int t1, int t2, const(DatasetRaw) * dset, Gretlopt opt)
{
    const(Model) * pmod;
    double fti, xti, xtid;
    int tdyn, nf = t2 - t1 + 1;
    int i, j, k, s, t, p;
    int lag, vj, m = 0;
    int fcols;

    fcols = 2 * var.neqns;

    /* rows = number of forecast periods; cols = 1 to hold forecast
       for each variable, plus 1 to hold variance for each variable
       if forecast is dynamic.
    */

    var.F = gretl_zero_matrix_new(nf, fcols);
    if (var.F == null) {
        return GretlError.alloc;
    }

    if (var.detflags & DetFlags.seas) {
        m = get_subperiod(t1, dset, null);
    }

    /* start of dynamic portion of forecast? */
    tdyn = VAR_get_tdyn(var, t1, t2, opt);

    for (t=t1, s=0; t<=t2; t++, s++) {
        int miss = 0;

        for (i=0; i<var.neqns; i++) {
            pmod = var.models[i];
            fti = 0.0;
            k = 0;

            /* constant, if present */
            if (pmod.ifc) {
                fti += pmod.coeff[k++];
            }

            /* lags of stochastic vars */
            for (j=0; j<var.neqns; j++) {
                vj = var.ylist[j+1];
                for (lag=1; lag<=var.order; lag++) {
                    if (!lag_wanted(var, lag)) {
                        continue;
                    }
                    p = t - lag;
                    xtid = NADBL;
                    if (p < tdyn || s - lag < 0) {
                        /* use actual data if possible */
                        if (p < 0) {
                            xti = NADBL;
                        } else {
                            xti = dset.Z[vj][p];
                        }
                    } else {
                        /* prior forecast value preferred */
                        if (p >= 0) {
                            xtid = dset.Z[vj][p];
                        }
                        xti = gretl_matrix_get(var.F, s-lag, j);
                    }
                    if (!na(xti)) {
                        fti += pmod.coeff[k] * xti;
                    } else if (!na(xtid)) {
                        fti += pmod.coeff[k] * xtid;
                    } else {
                        miss = 1;
                    }
                    k++;
                }
            }

            /* exogenous vars, if any */
            if (!miss && var.xlist != null) {
                for (j=1; j<=var.xlist[0]; j++) {
                    vj = var.xlist[j];
                    xti = dset.Z[vj][t];
                    if (na(xti)) {
                        miss = 1;
                    } else {
                        fti += pmod.coeff[k] * xti;
                    }
                    k++;
                }
            }

            /* other deterministics, if present */
            if (!miss) {
                if (var.detflags & DetFlags.seas) {
                    if (m < dset.pd - 1) {
                        fti += pmod.coeff[k+m];
                    }
                }
                if (var.detflags & DetFlags.trend) {
                    k = pmod.ncoeff - 1;
                    fti += pmod.coeff[k] * (t + 1);
                }
            }

            if (miss) {
                fti = NADBL;
            }

            gretl_matrix_set(var.F, s, i, fti);
        }

        if (var.detflags & DetFlags.seas) {
            m = (m < dset.pd - 1)? m + 1 : 0;
        }
    }

    VAR_add_fcast_variance(var, var.F, tdyn - t1);

    if (var.F != null) {
        gretl_matrix_set_t1(var.F, t1);
        gretl_matrix_set_t2(var.F, t2);
    }

    return 0;
}

int VECM_add_forecast(GretlVar * var, int t1, int t2, const(DatasetRaw) * dset, Gretlopt opt)
{
    GretlMatrixRaw * B = null;
    double s0 = 0, s1 = 1;
    int order = levels_order(var);
    int nexo = (var.xlist != null)? var.xlist[0] : 0;
    int nseas = var.jinfo.seasonals;
    int tdyn, nf = t2 - t1 + 1;
    int i, j, k, vj, s, t;
    int fcols, m = 0;

    fcols = 2 * var.neqns;

    var.F = gretl_zero_matrix_new(nf, fcols);
    if (var.F == null) {
        return GretlError.alloc;
    }

    B = VAR_coeff_matrix_from_VECM(var, 0);
    if (B == null) {
        gretl_matrix_free(var.F);
        var.F = null;
        return GretlError.alloc;
    }

    /* start of dynamic portion of forecast? */
    tdyn = VAR_get_tdyn(var, t1, t2, opt);

    if (nseas > 0) {
        m = get_subperiod(t1, dset, null);
        s1 -= 1.0 / dset.pd;
        s0 = s1 - 1;
    }

    for (t=t1, s=0; t<=t2; t++, s++) {

        for (i=0; i<var.neqns; i++) {
            double bij, xtj, xtjd, fti = 0.0;
            int col = 0;

            /* unrestricted constant, if present */
            if (var.ifc) {
                bij = gretl_matrix_get(B, i, col++);
                fti += bij;
            }

            /* lags of endogenous vars */
            for (j=0; j<var.neqns; j++) {
                vj = var.ylist[j+1];
                for (k=1; k<=order; k++) {
                    /* FIXME gappy lags */
                    if (t - k < 0) {
                        fti = NADBL;
                        break;
                    }
                    bij = gretl_matrix_get(B, i, col++);
                    xtjd = NADBL;
                    if (t < tdyn || s - k < 0) {
                        /* pre-forecast value */
                        xtj = dset.Z[vj][t-k];
                    } else {
                        /* prior forecast value preferred */
                        xtjd = dset.Z[vj][t-k];
                        xtj = gretl_matrix_get(var.F, s-k, j);
                    }
                    if (!na(xtj)) {
                        fti += bij * xtj;
                    } else if (!na(xtjd)) {
                        fti += bij * xtjd;
                    } else {
                        fti = NADBL;
                        break;
                    }
                }
                if (na(fti)) {
                    break;
                }
            }

            if (na(fti)) {
                goto set_fcast;
            }

            /* exogenous vars, if present */
            for (j=0; j<nexo; j++) {
                vj = var.xlist[j+1];
                xtj = dset.Z[vj][t];
                if (na(xtj)) {
                    fti = NADBL;
                } else {
                    bij = gretl_matrix_get(B, i, col++);
                    fti += bij * xtj;
                }
            }

            if (na(fti)) {
                goto set_fcast;
            }

            /* seasonals, if present */
            for (j=0; j<nseas; j++) {
                xtj = (m == j)? s1 : s0;
                bij = gretl_matrix_get(B, i, col++);
                fti += bij * xtj;
            }

            if (jcode(var) == JohansenCode.unrest_trend) {
                /* unrestricted trend */
                bij = gretl_matrix_get(B, i, col++);
                fti += bij * (t + 1);
            } else if (jcode(var) == JohansenCode.rest_const) {
                /* restricted constant */
                fti += gretl_matrix_get(B, i, col++);
            } else if (jcode(var) == JohansenCode.rest_trend) {
                /* restricted trend */
                bij = gretl_matrix_get(B, i, col++);
                fti += bij * t;
            }

            /* restricted exog vars */
            if (var.rlist != null) {
                for (j=0; j<var.rlist[0]; j++) {
                    vj = var.rlist[j+1];
                    xtj = dset.Z[vj][t-1]; /* ?? */
                    if (na(xtj)) {
                        fti = NADBL;
                    } else {
                        bij = gretl_matrix_get(B, i, col++);
                        fti += bij * xtj;
                    }
                }
            }

        set_fcast:
            gretl_matrix_set(var.F, s, i, fti);
        }

        if (nseas > 0) {
            m = (m < dset.pd - 1)? m + 1 : 0;
        }
    }

    gretl_matrix_free(B);

    VAR_add_fcast_variance(var, var.F, tdyn - t1);

    gretl_matrix_set_t1(var.F, t1);
    gretl_matrix_set_t2(var.F, t2);

    return 0;
}

const(GretlMatrixRaw) * gretl_VAR_get_forecast_matrix(GretlVar * var, int t1, int t2, DatasetRaw * dset, Gretlopt opt, int * err)
{
    if (var.F != null) {
        /* there's a forecast attached, but it may not be what we want */
        gretl_matrix_free(var.F);
        var.F = null;
    }

    if (var.ci == GretlCmdIndex.vecm) {
        *err = VECM_add_forecast(var, t1, t2, dset, opt);
    } else {
        *err = VAR_add_forecast(var, t1, t2, dset, opt);
    }

    return var.F;
}

void VAR_dw_rho(Model * pmod)
{
    double ut, u1;
    double xx = 0;
    double ut1 = 0, u11 = 0;
    int t, s;

    for (t=pmod.t1; t<=pmod.t2; t++)  {
        s = t - 1;
        if (s >= 0 && !na(pmod.uhat[s])) {
            ut = pmod.uhat[t];
            u1 = pmod.uhat[s];
            xx += (ut - u1) * (ut - u1);
            ut1 += ut * u1;
            u11 += u1 * u1;
        }
    }

    pmod.dw = xx / pmod.ess;
    pmod.rho = ut1 / u11;
}

/* set basic per-model statistics after estimation has been
   done by matrix methods */

int set_VAR_model_stats(GretlVar * var, int i)
{
    Model * pmod = var.models[i];
    const(double) * y;
    double u, x, SSR = 0, TSS = 0;
    int t;

    /* get the appropriate column of var->Y as a plain array */
    y = var.Y.val + i * var.T;

    pmod.ybar = gretl_mean(0, var.T - 1, y);
    pmod.sdy = gretl_stddev(0, var.T - 1, y);

    for (t=0; t<var.T; t++) {
        u = gretl_matrix_get(var.E, t, i);
        SSR += u * u;
        x = (var.ifc)? (y[t] - pmod.ybar) : y[t];
        TSS += x * x;
        pmod.uhat[t + pmod.t1] = u;
        pmod.yhat[t + pmod.t1] = y[t] - u;
    }

    pmod.ess = SSR;
    pmod.sigma = sqrt(SSR / pmod.dfd);
    pmod.tss = TSS;
    pmod.rsq = 1.0 - SSR / TSS;
    pmod.adjrsq = 1.0 - (SSR / (pmod.dfd)) / (TSS / (pmod.nobs - 1));
    pmod.fstt = ((TSS - SSR) / pmod.dfn) / (SSR / pmod.dfd);

    pmod.lnL = NADBL;

    VAR_dw_rho(pmod);

    return 0;
}

int gretl_VAR_do_error_decomp(const(GretlMatrixRaw) * S, GretlMatrixRaw * C, const(GretlMatrixRaw) * ord)
{
    int g = gretl_matrix_rows(S);
    GretlMatrixRaw * tmp = null;
    double x;
    int i, j, r, c;
    int err = 0;

    /* copy cross-equation covariance matrix */
    tmp = gretl_matrix_copy(S);
    if (tmp == null) {
        err = GretlError.alloc;
    }

    if (ord != null) {
        for (i=0; i<g; i++) {
            r = cast(int) (ord.val[i]);
            for (j=0; j<g; j++) {
                c = cast(int) (ord.val[j]);
                x = gretl_matrix_get(S, r, c);
                gretl_matrix_set(tmp, i, j, x);
            }
        }
    }

    /* lower-triangularize and decompose */
    if (!err) {
        for (i=0; i<g-1; i++) {
            for (j=i+1; j<g; j++) {
                gretl_matrix_set(tmp, i, j, 0.0);
            }
        }
        err = gretl_matrix_cholesky_decomp(tmp);
    }

    /* write the decomposition (lower triangle) into the C matrix */
    if (!err) {
        for (i=0; i<g; i++) {
            r = cast(int) ((ord != null)? ord.val[i] : i);
            for (j=0; j<=i; j++) {
                c = cast(int) ((ord != null)? ord.val[j] : j);
                x = gretl_matrix_get(tmp, i, j);
                gretl_matrix_set(C, r, c, x);
            }
        }
    }

    if (tmp != null) {
        gretl_matrix_free(tmp);
    }

    return err;
}

const(GretlMatrixRaw) * gretl_VAR_get_residual_matrix(GretlVar * var)
{
    if (var.E != null) {
        gretl_matrix_set_t1(var.E, var.t1);
        gretl_matrix_set_t2(var.E, var.t2);
    }

    return var.E;
}

bool samesize(const(GretlMatrixRaw) * p, const(GretlMatrixRaw) * q)
{
    return (p.rows == q.rows && p.cols == q.cols);
}

GretlMatrixRaw * gretl_VAR_get_fitted_matrix(GretlVar * var)
{
    GretlMatrixRaw * Yh = null;

    if (var.Y != null && var.E != null && samesize(var.Y, var.E)) {
        Yh = gretl_matrix_copy(var.Y);
        if (Yh != null) {
            gretl_matrix_subtract_from(Yh, var.E);
            gretl_matrix_set_t1(Yh, var.t1);
            gretl_matrix_set_t2(Yh, var.t2);
        }
    }

    return Yh;
}

GretlMatrixRaw * gretl_VAR_get_vma_matrix(GretlVar * var, const(DatasetRaw) * dset, int * err)
{
    int h = default_VAR_horizon(dset);
    int ar, n = var.neqns;
    int n2 = n * n;
    GretlMatrixRaw * VMA = null;
    GretlMatrixRaw * Tmp1, Tmp2;
    double x;
    int i, j, ii, jj;

    if (var.A == null) {
        /* companion matrix is absent */
        *err = GretlError.badstat;
        return null;
    }

    ar = var.A.rows;

    Tmp1 = gretl_identity_matrix_new(ar);
    Tmp2 = gretl_matrix_alloc(ar, ar);
    if (Tmp1 == null || Tmp2 == null) {
        *err = GretlError.alloc;
        goto bailout;
    }

    VMA = gretl_zero_matrix_new(h, n2);
    if (VMA == null) {
        *err = GretlError.alloc;
        goto bailout;
    }

    /* compose first row of VMA = vec(I(n))' */
    for (i=0; i<n2; i+=n+1) {
        gretl_matrix_set(VMA, 0, i, 1.0);
    }

    for (i=1; i<h; i++) {
        /* Tmp <-  A * Tmp */
        gretl_matrix_multiply(var.A, Tmp1, Tmp2);
        gretl_matrix_copy_values(Tmp1, Tmp2);
        /* VMA |= vec(Tmp[1:n,1:n])' */
        ii = jj = 0;
        for (j=0; j<n2; j++) {
            x = gretl_matrix_get(Tmp1, ii++, jj);
            gretl_matrix_set(VMA, i, j, x);
            if (ii == n) {
                jj++;
                ii = 0;
            }
        }
    }

 bailout:

    gretl_matrix_free(Tmp1);
    gretl_matrix_free(Tmp2);

    return VMA;
}

GretlMatrixRaw * gretl_VAR_get_portmanteau_test(GretlVar * var, int * err)
{
    GretlMatrixRaw * m = null;

    if (na(var.LB)) {
        *err = GretlError.badstat;
    } else {
        m = gretl_matrix_alloc(1, 4);
        if (m == null) {
            *err = GretlError.alloc;
        } else {
            int k = levels_order(var);
            int df = var.neqns * var.neqns * (var.LBs - k);
            double pv = chisq_cdf_comp(df, var.LB);
            char ** S = strings_array_new(4);

            m.val[0] = var.LB;
            m.val[1] = var.LBs;
            m.val[2] = df;
            m.val[3] = pv;
            if (S != null) {
                S[0] = gretl_strdup("LB");
                S[1] = gretl_strdup("maxlag");
                S[2] = gretl_strdup("df");
                S[3] = gretl_strdup("pvalue");
                gretl_matrix_set_colnames(m, S);
            }
            *err = 0;
        }
    }

    return m;
}

int gretl_VAR_get_variable_number(GretlVar * var, int k)
{
    int vnum = 0;

    if (var.models != null && k >= 0 && k < var.neqns) {
        if (var.models[k] != null && var.models[k].list != null) {
            vnum = var.models[k].list[1];
        }
    }

    return vnum;
}

int gretl_VAR_get_n_equations(GretlVar * var)
{
    return (var == null)? 0 : var.neqns;
}

int gretl_VAR_get_t1(GretlVar * var)
{
    return (var == null)? 0 : var.t1;
}

int gretl_VAR_get_t2(GretlVar * var)
{
    return (var == null)? 0 : var.t2;
}

int gretl_var_get_sample(GretlVar * var, int * t1, int * t2)
{
    if (var != null && var.models != null && var.models[0] != null) {
        Model * pmod = var.models[0];

        if (pmod.smpl.t1 >= 0 && pmod.smpl.t2 > pmod.smpl.t1) {
            *t1 = pmod.smpl.t1;
            *t2 = pmod.smpl.t2;
            return 0;
        }
    }

    return GretlError.data;
}

const(Model) * gretl_VAR_get_model(GretlVar * var, int i)
{
    if (var != null && var.models != null && i < var.neqns) {
        return var.models[i];
    } else {
        return null;
    }
}

int periods_from_pd(int pd)
{
    int periods = 10;

    if (pd == 4) {
        /* quarterly: try 5 years */
        periods = 20;
    } else if (pd == 12) {
        /* monthly: two years */
        periods = 24;
    } else if (pd == 7 || pd == 6 || pd == 5) {
        /* daily: three weeks */
        periods = 3 * pd;
    }

    return periods;
}

int default_VAR_horizon(const(DatasetRaw) * dset)
{
    int h = libset_get_int(SetKey.horizon);

    if (h <= 0) {
    if (dset != null) {
        h = periods_from_pd(dset.pd);
    } else {
        h = 20; /* default */
    }
    }

    return h;
}

GretlMatrixRaw * reorder_responses(GretlVar * var, int * err)
{
    GretlMatrixRaw * S, C;
    int i, j, r, c;
    double x;

    S = gretl_matrix_copy(var.S);
    C = gretl_matrix_alloc(var.neqns, var.neqns);

    if (S == null || C == null) {
        gretl_matrix_free(S);
        gretl_matrix_free(C);
        *err = GretlError.alloc;
        return null;
    }

    /* rearrange @S */
    for (i=0; i<var.neqns; i++) {
        r = cast(int) (var.ord.val[i]);
        for (j=0; j<var.neqns; j++) {
            c = cast(int) (var.ord.val[j]);
            x = gretl_matrix_get(var.S, r, c);
            gretl_matrix_set(S, i, j, x);
        }
    }

    /* decompose rearranged @S */
    gretl_matrix_cholesky_decomp(S);

    /* C = rearranged decomposed S */
    for (i=0; i<var.neqns; i++) {
        r = cast(int) (var.ord.val[i]);
        for (j=0; j<var.neqns; j++) {
            c = cast(int) (var.ord.val[j]);
            x = gretl_matrix_get(S, i, j);
            gretl_matrix_set(C, r, c, x);
        }
    }

    gretl_matrix_free(S);

    return C;
}

void copy_north_west(GretlMatrixRaw * targ, const(GretlMatrixRaw) * src, int add)
{
    int i, j, k = src.cols;
    double xij;

    if (!add) {
        gretl_matrix_zero(targ);
    }

    for (j=0; j<k; j++) {
        for (i=0; i<k; i++) {
            xij = gretl_matrix_get(src, i, j);
            if (add) {
                xij += gretl_matrix_get(targ, i, j);
            }
            gretl_matrix_set(targ, i, j, xij);
        }
    }
}

int real_point_responses(const(GretlMatrixRaw) * A, const(GretlMatrixRaw) * C, GretlMatrixRaw * resp, int targ, int shock)
{
    GretlMatrixRaw * rtmp;
    GretlMatrixRaw * ctmp;
    double rij;
    int n = C.rows;
    int np = A.cols;
    int i, j, k, t;

    /* workspace */
    rtmp = gretl_matrix_alloc(np, n);
    ctmp = gretl_zero_matrix_new(np, n);
    if (rtmp == null || ctmp == null) {
    gretl_matrix_free(rtmp);
    gretl_matrix_free(ctmp);
    return GretlError.alloc;
    }

    for (t=0; t<resp.rows; t++) {
        if (t == 0) {
            /* initial estimated responses */
        copy_north_west(rtmp, C, 0);
        } else {
            /* calculate further estimated responses */
            gretl_matrix_multiply(A, rtmp, ctmp);
            gretl_matrix_copy_values(rtmp, ctmp);
        }
        if (resp.cols == 1) {
            resp.val[t] = gretl_matrix_get(rtmp, targ, shock);
        } else if (shock >= 0) {
            /* all targets for one shock */
            for (i=0; i<n; i++) {
                rij = gretl_matrix_get(rtmp, i, shock);
                gretl_matrix_set(resp, t, i, rij);
            }
        } else if (targ >= 0) {
            /* all shocks for one target */
            for (j=0; j<n; j++) {
                rij = gretl_matrix_get(rtmp, targ, j);
                gretl_matrix_set(resp, t, j, rij);
            }
        } else {
            /* all shocks, all targets */
            k = 0;
            for (i=0; i<n; i++) {
                for (j=0; j<n; j++) {
                    rij = gretl_matrix_get(rtmp, i, j);
                    gretl_matrix_set(resp, t, k++, rij);
                }
            }
        }
    }

    gretl_matrix_free(rtmp);
    gretl_matrix_free(ctmp);

    return 0;
}

/* Calculate the VMA representation based on the A and C matrices,
   for a given horizon.
*/

GretlMatrixRaw * vma_rep(GretlMatrixRaw * A, GretlMatrixRaw * C, int horizon, int * err)
{
    GretlMatrixRaw * resp = null;
    GretlMatrixRaw * realC = C;
    GretlMatrixRaw * realA = A;
    int neqns, dim;

    if (horizon <= 0) {
    *err = GretlError.invarg;
    return null;
    }

    neqns = A.rows;
    dim = A.cols;

    if (realC == null) {
    /* manufacture "plain" @C for convenience */
    realC = gretl_identity_matrix_new(neqns);
    if (realC == null) {
        *err = GretlError.alloc;
    }
    }

    if (!*err && dim > neqns) {
    /* the incoming @A is just the top @neqns rows
       of the companion matrix
    */
    realA = companionize(A, err);
    }

    if (!*err) {
    int nresp = neqns * neqns;

    resp = gretl_matrix_alloc(horizon, nresp);

    if (resp == null) {
        *err = GretlError.alloc;
    } else {
        if (neqns == 1) {
        /* "all" means "first" in this case */
        *err = real_point_responses(realA, realC, resp, 0, 0);
        } else {
        *err = real_point_responses(realA, realC, resp, -1, -1);
        }
    }
    }

    if (realC != C) {
    gretl_matrix_free(realC);
    }
    if (realA != A) {
    gretl_matrix_free(realA);
    }
    if (*err && resp != null) {
        gretl_matrix_free(resp);
        resp = null;
    }

    return resp;
}

GretlMatrixRaw * gretl_VAR_get_point_responses(GretlVar * var, int targ, int shock, int periods, int * err)
{
    GretlMatrixRaw * resp = null;
    GretlMatrixRaw * C = var.C;
    int nresp = 1;

    if (shock >= var.neqns) {
        fprintf(stderr, "Shock variable out of bounds\n".ptr);
        *err = GretlError.invarg;
    }
    if (!*err && targ >= var.neqns) {
        fprintf(stderr, "Target variable out of bounds\n".ptr);
        *err = GretlError.invarg;
    }
    if (!*err && periods <= 0) {
        fprintf(stderr, "Invalid number of periods\n".ptr);
        *err = GretlError.invarg;
    }

    if (!*err && var.ord != null) {
        C = reorder_responses(var, err);
    }

    if (*err) {
    return null;
    }

    if (targ < 0 && shock < 0) {
    /* compute all */
    nresp = var.neqns * var.neqns;
    } else if (targ < 0 || shock < 0) {
    nresp = var.neqns;
    }

    resp = gretl_matrix_alloc(periods, nresp);

    if (resp == null) {
        *err = GretlError.alloc;
    } else {
    *err = real_point_responses(var.A, C, resp, targ, shock);
    }

    if (C != var.C) {
        gretl_matrix_free(C);
    }

    if (*err && resp != null) {
        gretl_matrix_free(resp);
        resp = null;
    }

    return resp;
}

/**
 * gretl_VAR_get_impulse_response:
 * @var: pointer to VAR struct.
 * @targ: index of the target or response variable.
 * @shock: index of the source or shock variable.
 * @periods: number of periods over which to compute the response.
 * @alpha: determines confidence level for bootstrap interval.
 * @dset: dataset struct.
 * @err: location to receive error code.
 *
 * Computes the response of @targ to a perturbation of @shock
 * in the context of @var: @targ and @shock are zero-based indices
 * relative to the structure of @var.  For example if @targ = 0 and
 * @shock = 1, we compute the response of the dependent variable in
 * the first VAR equation to a perturbation of the variable that
 * appears as dependent in the second VAR equation.
 *
 * If @alpha is zero, the response matrix returned is a column vector
 * of length @periods, giving the point estimate of the response
 * function.  Otherwise, the response matrix returned has
 * three columns, containing the point estimate, the @alpha / 2
 * quantile and the 1 - @alpha / 2 quantile, where the quantiles
 * are based on bootstrap replication, with resampling of the
 * original residuals with replacement.
 *
 * Returns: matrix containing the estimated impulse responses,
 * with or without a confidence interval.
 */

GretlMatrixRaw * gretl_VAR_get_impulse_response(GretlVar * var, int targ, int shock, int periods, double alpha, const(DatasetRaw) * dset, int * err)
{
    GretlMatrixRaw * point = null;
    GretlMatrixRaw * full = null;
    GretlMatrixRaw * ret = null;

    if (periods == 0) {
        if (dset != null) {
            periods = default_VAR_horizon(dset);
        } else {
            *err = GretlError.data;
            return null;
        }
    }

    if (alpha != 0 && (alpha < 0.01 || alpha > 0.6)) {
        *err = GretlError.data;
    }

    point = gretl_VAR_get_point_responses(var, targ, shock, periods, err);

    if (dset == null || dset.Z == null || alpha == 0.0) {
        /* no data matrix provided, or no alpha given:
           just return point estimate */
        ret = point;
    } else if (point != null) {
        full = irf_bootstrap(var, targ, shock, periods,
                             alpha, point, dset, err);
    gretl_matrix_free(point);
        ret = full;
    }

    return ret;
}

GretlMatrixRaw * gretl_VAR_get_fcast_se(GretlVar * var, int periods)
{
    int k = var.neqns * levels_order(var);
    GretlMatrixRaw * vtmp = null;
    GretlMatrixRaw * vt = null;
    GretlMatrixRaw * se = null;
    double vti;
    int i, t;

    if (periods <= 0) {
        fprintf(stderr, "Invalid number of periods\n".ptr);
        return null;
    }

    se = gretl_zero_matrix_new(periods, var.neqns);
    vt = gretl_matrix_alloc(k, k);
    vtmp = gretl_matrix_alloc(k, k);

    if (se == null || vt == null || vtmp == null) {
        gretl_matrix_free(se);
        gretl_matrix_free(vt);
        gretl_matrix_free(vtmp);
        return null;
    }

    for (t=0; t<periods; t++) {
        if (t == 0) {
            /* initial variance */
            copy_north_west(vt, var.S, 0);
        } else {
            /* calculate further variances */
            gretl_matrix_copy_values(vtmp, vt);
            gretl_matrix_qform(var.A, GretlMatrixMod.none,
                               vtmp, vt, GretlMatrixMod.none);
            copy_north_west(vt, var.S, 1);
        }
        for (i=0; i<var.neqns; i++) {
            vti = gretl_matrix_get(vt, i, i);
            gretl_matrix_set(se, t, i, sqrt(vti));
        }
    }

    gretl_matrix_free(vt);
    gretl_matrix_free(vtmp);

    return se;
}

GretlMatrixRaw * gretl_VAR_get_fcast_decomp(GretlVar * var, int targ, int periods, int * err)
{
    int n = var.neqns;
    int k = n * levels_order(var);
    GretlMatrixBlock * B;
    GretlMatrixRaw * idx, vtmp;
    GretlMatrixRaw * cic, vt;
    GretlMatrixRaw * vd = null;
    GretlMatrixRaw * C = var.C;
    int i, t;

    *err = 0;

    if (targ >= n) {
        fprintf(stderr, "Target variable out of bounds\n".ptr);
        *err = GretlError.data;
    }

    if (!*err && periods <= 0) {
        fprintf(stderr, "Invalid number of periods\n".ptr);
        *err = GretlError.data;
    }

    if (!*err && var.ord != null) {
        C = reorder_responses(var, err);
    }

    if (*err) {
    return null;
    }

    B = gretl_matrix_block_new(&idx, n, n,
                               &cic, n, n,
                               &vt,  k, k,
                               &vtmp, k, k,
                               null);
    if (B == null) {
        *err = GretlError.alloc;
        goto bailout;
    }

    vd = gretl_zero_matrix_new(periods, n + 1);
    if (vd == null) {
        *err = GretlError.alloc;
        goto bailout;
    }

    gretl_matrix_zero(idx);

    for (i=0; i<n && !*err; i++) {
        double vti;

        /* adjust index matrix */
        if (i > 0) {
            gretl_matrix_set(idx, i-1, i-1, 0.0);
        }
        gretl_matrix_set(idx, i, i, 1.0);

        for (t=0; t<periods && !*err; t++) {
            if (t == 0) {
                /* calculate initial variances */
                *err = gretl_matrix_qform(C, GretlMatrixMod.none,
                                          idx, cic, GretlMatrixMod.none);
                copy_north_west(vt, cic, 0);
            } else {
                /* calculate further variances */
                gretl_matrix_copy_values(vtmp, vt);
                *err = gretl_matrix_qform(var.A, GretlMatrixMod.none,
                                          vtmp, vt, GretlMatrixMod.none);
                copy_north_west(vt, cic, 1);
            }
            if (!*err) {
                vti = gretl_matrix_get(vt, targ, targ);
                gretl_matrix_set(vd, t, i, vti);
            }
        }
    }

    for (t=0; t<periods && !*err; t++) {
        double vi, vtot = 0.0;

        for (i=0; i<n; i++) {
            vtot += gretl_matrix_get(vd, t, i);
        }
        /* normalize variance contributions as % shares */
        for (i=0; i<n; i++) {
            vi = gretl_matrix_get(vd, t, i);
            gretl_matrix_set(vd, t, i, 100.0 * vi / vtot);
        }

        gretl_matrix_set(vd, t, var.neqns, sqrt(vtot));
    }

 bailout:

    gretl_matrix_block_destroy(B);

    if (C != var.C) {
        gretl_matrix_free(C);
    }

    if (*err && vd != null) {
        gretl_matrix_free(vd);
        vd = null;
    }

    return vd;
}

GretlMatrixRaw * gretl_VAR_get_FEVD_matrix(GretlVar * var, int targ, int shock, int horizon, const(DatasetRaw) * dset, int * err)
{
    GretlMatrixRaw * vd, V;
    double vjk;
    int h = horizon;
    int n = var.neqns;
    int imin, imax;
    int i, j, k, kk;

    if (h <= 0) {
        h = default_VAR_horizon(dset);
    }

    if (targ < 0) {
        /* do the whole thing */
        k = n * n;
        imin = 0;
        imax = n;
    } else {
        /* got a specific target */
        k = n;
        imin = targ;
        imax = targ + 1;
    }

    V = gretl_matrix_alloc(h, k);
    if (V == null) {
        *err = GretlError.alloc;
        return null;
    }

    kk = 0;
    for (i=imin; i<imax && !*err; i++) {
        vd = gretl_VAR_get_fcast_decomp(var, i, h, err);
        if (!*err) {
            for (k=0; k<n; k++) {
                for (j=0; j<h; j++) {
                    vjk = gretl_matrix_get(vd, j, k);
                    gretl_matrix_set(V, j, kk, vjk / 100.0);
                }
                kk++;
            }
            gretl_matrix_free(vd);
        }
    }

    /* If @shock is specific, not all, we now proceed to carve
       out the columns of @V that are actually wanted. This is
       not as efficient as it might be -- we're computing more
       than we need above -- but it's simple and will do for
       the present.
    */

    if (!*err && shock >= 0) {
        GretlMatrixRaw * V2 = null;

        if (targ >= 0) {
            /* just one column wanted */
            V2 = gretl_column_vector_alloc(h);
            if (V2 != null) {
                for (j=0; j<h; j++) {
                    V2.val[j] = gretl_matrix_get(V, j, shock);
                }
            }
        } else {
            /* n columns wanted */
            V2 = gretl_matrix_alloc(h, n);
            if (V2 != null) {
                k = shock;
                for (j=0; j<n; j++) {
                    for (i=0; i<h; i++) {
                        vjk = gretl_matrix_get(V, i, k);
                        gretl_matrix_set(V2, i, j, vjk);
                    }
                    k += n;
                }
            }
        }
        if (V2 == null) {
            *err = GretlError.alloc;
        } else {
            gretl_matrix_free(V);
            V = V2;
        }
    }

    if (*err && V != null) {
        gretl_matrix_free(V);
        V = null;
    }

    return V;
}

GretlMatrixRaw * gretl_VAR_get_full_FEVD_matrix(GretlVar * var, const(DatasetRaw) * dset, int * err)
{
    return gretl_VAR_get_FEVD_matrix(var, -1, -1, -1, dset, err);
}

int var_max_order(const(int) * list, const(DatasetRaw) * dset)
{
    int T = dset.t2 - dset.t1 + 1;
    int nstoch = 0, ndet = 0;
    int gotsep = 0;
    int order = 1;
    int i;

    for (i=1; i<=list[0]; i++) {
        if (list[i] == LISTSEP) {
            gotsep = 1;
            continue;
        }
        if (!gotsep) {
            nstoch++;
        } else {
            ndet++;
        }
    }

    order = (T - ndet) / nstoch;

    while (order > 0) {
        int t1 = (order > dset.t1)? order : dset.t1;

        T = dset.t2 - t1 + 1;
        if (nstoch * order + ndet > T) {
            order--;
        } else {
            break;
        }
    }

    return order - 1;
}

int VAR_add_roots(GretlVar * var)
{
    int err = 0;

    if (var.A == null) {
        fprintf(stderr, "VAR_add_roots: var->A is missing\n".ptr);
        return GretlError.data;
    }

    /* save eigenvalues of companion form matrix */
    var.L = legacy_general_matrix_eigenvals(var.A, &err);

    if (err) {
        gretl_matrix_free(var.L);
        var.L = null;
    }

    return err;
}

const(GretlMatrixRaw) * gretl_VAR_get_roots(GretlVar * var, int * err)
{
    if (var == null) {
        fprintf(stderr, "gretl_VAR_get_roots: VAR is NULL\n".ptr);
        *err = GretlError.data;
        return null;
    }

    if (var.L == null) {
        /* roots not computed yet */
        *err = VAR_add_roots(var);
    }

    return var.L;
}

/* used in context of LR tests: see vartest.c */

double gretl_VAR_ldet(GretlVar * var, const(GretlMatrixRaw) * E, int * err)
{
    GretlMatrixRaw * S = null;
    double ldet = NADBL;

    S = gretl_matrix_alloc(var.neqns, var.neqns);

    if (S == null) {
        *err = GretlError.alloc;
    } else {
        gretl_matrix_multiply_mod(E, GretlMatrixMod.transpose,
                                  E, GretlMatrixMod.none,
                                  S, GretlMatrixMod.none);
        gretl_matrix_divide_by_scalar(S, var.T);
        ldet = gretl_vcv_log_determinant(S, err);
        gretl_matrix_free(S);
    }

    return ldet;
}

/* identify columns of the VAR X matrix that contain the final_
   lag of an endogenous variable */

int omit_column(GretlVar * var, int nl, int j)
{
    if (var.ifc) {
        return j % nl == 0 && j < 1 + var.neqns * nl;
    } else {
        return (j+1) % nl == 0 && j < var.neqns * nl;
    }
}

/* make and record residuals for LR test on last lag */

GretlMatrixRaw * VAR_short_residuals(GretlVar * var, int * err)
{
    GretlMatrixRaw * X = null;
    GretlMatrixRaw * B = null;
    GretlMatrixRaw * E = null;
    double x;
    /* note: removing one lag from each equation */
    int g = var.ncoeff - var.neqns;
    int j, t, k, nl;

    E = gretl_matrix_alloc(var.T, var.neqns);
    X = gretl_matrix_alloc(var.T, g);
    B = gretl_matrix_alloc(g, var.neqns);

    if (E == null || X == null || B == null) {
        *err = GretlError.alloc;
        goto bailout;
    }

    nl = var_n_lags(var);

    k = 0;
    for (j=0; j<var.ncoeff; j++) {
        /* loop across the cols of var->X */
        if (j > 0 && omit_column(var, nl, j)) {
            continue;
        }
        for (t=0; t<var.T; t++) {
            x = gretl_matrix_get(var.X, t, j);
            gretl_matrix_set(X, t, k, x);
        }
        k++;
    }

    /* put residuals from "short" estimation into E */
    *err = gretl_matrix_multi_ols(var.Y, X, B, E, null);

 bailout:

    gretl_matrix_free(X);
    gretl_matrix_free(B);

    if (*err) {
        gretl_matrix_free(E);
        E = null;
    }

    return E;
}

int VAR_add_stats(GretlVar * var, int code)
{
    GretlMatrixRaw * E1 = null;
    int err = 0;

    if (var.order > 1 && code == VarTask.var_estimate) {
        E1 = VAR_short_residuals(var, &err);
    }

    var.S = gretl_matrix_alloc(var.neqns, var.neqns);
    if (var.S == null) {
        err = GretlError.alloc;
    }

    if (!err) {
        gretl_matrix_multiply_mod(var.E, GretlMatrixMod.transpose,
                                  var.E, GretlMatrixMod.none,
                                  var.S, GretlMatrixMod.none);
        /* for computing log-determinant, don't apply df
           correction (?) */
        gretl_matrix_divide_by_scalar(var.S, var.T);
    }

    if (!err) {
        var.ldet = gretl_vcv_log_determinant(var.S, &err);

    }

    if (!err) {
        int T = var.T;
        int g = var.neqns;
        int p = var.ncoeff;
        int k = g * p;

        var.ll = -(g * T / 2.0) * (ConstId.pi + 1) - (T / 2.0) * var.ldet;
        var.AIC = (-2.0 * var.ll + 2.0 * k) / T;
        var.BIC = (-2.0 * var.ll + k * log(T)) / T;
        var.HQC = (-2.0 * var.ll + 2.0 * k * log(log(T))) / T;
    }

    if (E1 != null) {
        if (!err) {
            VAR_LR_lag_test(var, E1);
        }
        gretl_matrix_free(E1);
    }

    if (!err) {
        VAR_portmanteau_test(var);
    }

    return err;
}

void gretl_VAR_param_names(GretlVar * v, char ** params, const(DatasetRaw) * dset)
{
    char[12] lagstr;
    int i, j, n, k = 0;

    if (v.detflags & DetFlags.const_) {
        strcpy(params[k++], dset.varname[0]);
    }

    for (i=1; i<=v.ylist[0]; i++) {
        for (j=1; j<=v.order; j++) {
            if (!lag_wanted(v, j)) {
                continue;
            }
            sprintf(lagstr.ptr, "_%d".ptr, j);
            n = cast(int) strlen(lagstr.ptr);
            if (v.ci == GretlCmdIndex.vecm) {
                strcpy(params[k], "d_");
                n += 2;
            }
            strncat(params[k], dset.varname[v.ylist[i]],
                    VNAMELEN - n - 1);
            strncat(params[k], lagstr.ptr, n);
            k++;
        }
    }

    if (v.xlist != null) {
        for (i=1; i<=v.xlist[0]; i++) {
            strcpy(params[k++], dset.varname[v.xlist[i]]);
        }
    }

    if (v.detflags & DetFlags.seas) {
        for (i=1; i<dset.pd; i++) {
            sprintf(params[k++], "S%d".ptr, i);
        }
    }

    if (v.detflags & DetFlags.trend) {
        strcpy(params[k++], "time");
    }

    if (v.ci == GretlCmdIndex.vecm) {
        int rank = jrank(v);

        for (i=0; i<rank; i++) {
            sprintf(params[k++], "EC%d".ptr, i+1);
        }
    }
}

/* public because called from irfboot.c */

void VAR_write_A_matrix(GretlVar * v, GretlMatrixMod mod)
{
    int i, ii, j, k, lag;
    int n = v.neqns;
    int dim = n * v.order;
    double bij;

    for (j=0; j<n; j++) {
        k = lag = ii = 0;
        for (i=0; i<dim; i++) {
            if (lag_wanted(v, lag+1)) {
                bij = gretl_matrix_get(v.B, ii + v.ifc, j);
                ii++;
            } else {
                bij = 0;
            }
        if (mod == GretlMatrixMod.transpose) {
        gretl_matrix_set(v.A, n * lag + k, j, bij);
        } else {
        gretl_matrix_set(v.A, j, n * lag + k, bij);
        }
            if (lag < v.order - 1) {
                lag++;
            } else {
                lag = 0;
                k++;
            }
        }
    }
}

int VAR_depvar_name(GretlVar * var, int i, const(char) * yname)
{
    Model * pmod = var.models[i];

    if (var.ci == GretlCmdIndex.var) {
        pmod.depvar = gretl_strdup(yname);
    } else {
        pmod.depvar = cast(typeof(pmod.depvar)) malloc(VNAMELEN);
        if (pmod.depvar != null) {
            strcpy(pmod.depvar, "d_");
            strncat(pmod.depvar, yname, VNAMELEN - 3);
        }
    }

    return (pmod.depvar == null)? GretlError.alloc: 0;
}

/* transcribe the per-equation output from a VAR or VECM into
   MODEL structs, if wanted */

int transcribe_VAR_models(GretlVar * var, const(DatasetRaw) * dset, const(GretlMatrixRaw) * XTX)
{
    Model * pmod;
    char ** params = null;
    int yno, N = dset.n;
    int ecm = (var.ci == GretlCmdIndex.vecm);
    double x;
    int i, j, jmax;
    int err = 0;

    params = strings_array_new_with_length(var.ncoeff, VNAMELEN);
    if (params == null) {
        return GretlError.alloc;
    }

    gretl_VAR_param_names(var, params, dset);

    jmax = (var.B != null)? var.B.rows : 0;

    for (i=0; i<var.neqns && !err; i++) {
        yno = var.ylist[i+1];

        pmod = var.models[i];
        pmod.ID = i + 1;
        pmod.ci = (ecm)? GretlCmdIndex.ols : GretlCmdIndex.var;
        pmod.aux = (ecm)? ModelAuxCode.vecm: ModelAuxCode.var;

        pmod.full_n = N;
        pmod.nobs = var.T;
        pmod.t1 = var.t1;
        pmod.t2 = var.t2;
        pmod.ncoeff = var.ncoeff;
        pmod.ifc = var.ifc;
        pmod.dfn = var.ncoeff - pmod.ifc;
        pmod.dfd = (ecm)? var.df : pmod.nobs - pmod.ncoeff;

        err = gretl_model_allocate_storage(pmod);

        if (!err) {
            VAR_depvar_name(var, i, dset.varname[yno]);
            if (i == 0) {
                pmod.params = params;
            } else {
                pmod.params = strings_array_dup(params, var.ncoeff);
                if (pmod.params == null) {
                    err = GretlError.alloc;
                }
            }
        }

        if (!err) {
            pmod.nparams = var.ncoeff;
            pmod.list = gretl_list_new(1);
            if (pmod.list == null) {
                err = GretlError.alloc;
            }
        }

        if (!err) {
            pmod.list[1] = yno;
            set_VAR_model_stats(var, i);
            for (j=0; j<jmax; j++) {
                pmod.coeff[j] = gretl_matrix_get(var.B, j, i);
                if (XTX != null) {
                    if (XTX.rows <= var.ncoeff) {
                        x = gretl_matrix_get(XTX, j, j);
                        pmod.sderr[j] = pmod.sigma * sqrt(x);
                    } else {
                        int jj = i * var.ncoeff + j;

                        x = gretl_matrix_get(XTX, jj, jj);
                        pmod.sderr[j] = pmod.sigma * sqrt(x);
                    }
                }
            }
        }
    }

    return err;
}

int * lags_from_laglist(const(int) * llist, int * err)
{
    int * lags = null;

    if (llist[0] == 0) {
        *err = GretlError.data;
    } else {
        lags = gretl_list_copy(llist);
        if (lags == null) {
            *err = GretlError.alloc;
        }
    }

    if (lags != null) {
        gretl_list_sort(lags);
        if (lags[1] < 1) {
            *err = GretlError.data;
            free(lags);
            lags = null;
        }
    }

    return lags;
}

/**
 * gretl_VAR:
 * @order: lag order for the VAR.
 * @laglist: specific list of lags, or NULL.
 * @list: specification for the first model in the set.
 * @dset: dataset struct.
 * @opt: if includes %OPT_R, use robust VCV;
 *       if includes %OPT_H, use HAC VCV;
 *       if includes %OPT_I, print impulse responses;
 *       if includes %OPT_F, print forecast variance decompositions;
 *       if includes %OPT_D, add seasonal dummies;
 *       if includes %OPT_M, set min lag for lag-selection.
 *       if includes %OPT_N, do not include a constant.
 *       if includes %OPT_Q, do not show individual regressions.
 *       if includes %OPT_T, include a linear trend.
 *       if includes %OPT_L, test for optimal lag length (only).
 *       if includes %OPT_S, silent (no printing).
 * @prn: gretl printing struct.
 * @err: location to receive error code.
 *
 * Estimate a vector auto-regression (VAR), print and save
 * the results.
 *
 * Returns: pointer to VAR struct, which may be %NULL on error.
 */

GretlVar * real_gretl_VAR(int order, int * laglist, int * list, const(DatasetRaw) * dset, Gretlopt opt, Prn * prn, int * err, GretlMatrixRaw ** pm)
{
    GretlVar * var = null;
    int code = (opt & Gretlopt.OPT_L)? VarTask.lagsel : VarTask.var_estimate;
    int * lags = null;

    if ((opt & Gretlopt.OPT_M) && !(opt & Gretlopt.OPT_L)) {
    gretl_errmsg_sprintf("%s: inapplicable option".ptr, "--minlag".ptr);
    *err = GretlError.badopt;
    return null;
    }

    if (laglist != null) {
        lags = lags_from_laglist(laglist, err);
        if (*err) {
            return null;
        }
    }

    /* allocation and initial set-up */
    var = gretl_VAR_new(code, order, 0, lags, list, dset,
                        opt, err);
    if (var == null) {
        return null;
    }

    /* run the regressions */
    if (getenv("VAR_USE_QR") != null) {
        *err = gretl_matrix_QR_ols(var.Y, var.X,
                                   var.B, var.E,
                                   &var.XTX, null);
    } else {
        /* use Cholesky or QR as needed */
        *err = gretl_matrix_multi_ols(var.Y, var.X,
                                      var.B, var.E,
                                      &var.XTX);
    }

    if (!*err) {
        if (code == VarTask.lagsel) {
            /* doing lag-length selection */
            *err = VAR_add_stats(var, code);
            if (!*err) {
                *err = VAR_do_lagsel(var, dset, opt, prn, pm);
            }
        } else {
            /* regular VAR estimation */
            *err = transcribe_VAR_models(var, dset, null);

            if (!*err) {
                VAR_write_A_matrix(var, GretlMatrixMod.none);
                /* note: the following also sets std. errors */
                *err = VAR_wald_omit_tests(var);
            }

            if (!*err) {
                *err = VAR_add_stats(var, code);
            }

            if (!*err) {
                *err = gretl_VAR_do_error_decomp(var.S, var.C, null);
            }

            if (!*err && prn != null) {
                gretl_VAR_print(var, dset, opt, prn);
            }
        }
    }

    if (code == VarTask.lagsel || (*err && var != null)) {
        gretl_VAR_free(var);
        var = null;
    }

    return var;
}

GretlVar * gretl_VAR(int order, int * laglist, int * list, const(DatasetRaw) * dset, Gretlopt opt, Prn * prn, int * err)
{
    return real_gretl_VAR(order, laglist, list, dset,
              opt, prn, err, null);
}

int gui_VAR_lagsel(int order, int * list, const(DatasetRaw) * dset, Gretlopt opt, Prn * prn, GretlMatrixRaw ** pm)
{
    int err = 0;

    real_gretl_VAR(order, null, list, dset,
           opt, prn, &err, pm);
    return err;
}

void print_johansen_sigmas(const(JohansenInfo) * jv, Prn * prn)
{
    int i, j;

    pprintf(prn, "%s\n\n".ptr, _("Sample variance-covariance matrices for residuals").ptr);

    pprintf(prn, " %s (S00)\n\n".ptr, _("VAR system in first differences").ptr);
    for (i=0; i<jv.S00.rows; i++) {
        for (j=0; j<jv.S00.rows; j++) {
            pprintf(prn, "%#12.5g".ptr, gretl_matrix_get(jv.S00, i, j));
        }
        pputc(prn, '\n');
    }

    pprintf(prn, "\n %s (S11)\n\n".ptr, _("System with levels as dependent variable").ptr);
    for (i=0; i<jv.S11.rows; i++) {
        for (j=0; j<jv.S11.rows; j++) {
            pprintf(prn, "%#12.5g".ptr, gretl_matrix_get(jv.S11, i, j));
        }
        pputc(prn, '\n');
    }

    pprintf(prn, "\n %s (S01)\n\n".ptr, _("Cross-products").ptr);
    for (i=0; i<jv.S01.rows; i++) {
        for (j=0; j<jv.S01.cols; j++) {
            pprintf(prn, "%#12.5g".ptr, gretl_matrix_get(jv.S01, i, j));
        }
        pputc(prn, '\n');
    }

    pputc(prn, '\n');
}

/* called from gretl_restriction_finalize() if the target
   is a VECM */

int gretl_VECM_test(GretlVar * vecm, GretlRestriction * rset, const(DatasetRaw) * dset, Gretlopt opt, Prn * prn)
{
    alias JfunFn = extern (C) int function(GretlVar *, GretlRestriction *, const(DatasetRaw) *, Gretlopt opt, Prn *);
    JfunFn jfun;
    int err = 0;

    if (vecm.jinfo == null || rset == null) {
        return GretlError.data;
    }

    gretl_error_clear();

    jfun = cast(JfunFn) get_plugin_function("vecm_test_restriction");

    if (jfun == null) {
        err = 1;
    } else {
        err = jfun(vecm, rset, dset, opt, prn);
    }

    return err;
}

int var_add_full_vcv(GretlVar * var)
{
    GretlMatrixRaw * V;
    Model * pmod;
    double vij;
    int i, k, nc;
    int bi, bj;
    int vi, vj;
    int err = 0;

    nc = var.neqns * var.ncoeff;
    V = gretl_zero_matrix_new(nc, nc);
    if (V == null) {
        return GretlError.alloc;
    }

    k = var.ncoeff;

    vi = vj = 0;
    for (i=0; i<var.neqns; i++) {
        pmod = var.models[i];
        if (pmod.vcv == null) {
            err = makevcv(pmod, pmod.sigma);
            if (err) {
                break;
            }
        }
        for (bi=0; bi<k; bi++) {
            for (bj=0; bj<k; bj++) {
                vij = gretl_model_get_vcv_element(pmod, bi, bj, k);
                gretl_matrix_set(V, vi+bi, vj+bj, vij);
            }
        }
        vi += k;
        vj += k;
    }

    if (err) {
        gretl_matrix_free(V);
    } else {
        var.V = V;
    }

    return err;
}

/* called from gretl_restriction_finalize() if the target
   is a plain VAR */

int gretl_VAR_test(GretlVar * var, GretlRestriction * rset, const(DatasetRaw) * dset, Gretlopt opt, Prn * prn)
{
    int err = 0;

    if (rset == null) {
        return GretlError.data;
    }

    gretl_error_clear();

    if (var.V == null) {
        err = var_add_full_vcv(var);
    }

    if (!err) {
        const(GretlMatrixRaw) * R = rset_get_R_matrix(rset);
        const(GretlMatrixRaw) * q = rset_get_q_matrix(rset);
        int dfu = var.T - var.ncoeff;
        int r = var.B.rows;
        int c = var.B.cols;

        /* pass var->B as a column vector */
        gretl_matrix_reuse(var.B, r*c, 1);
        err = multi_eqn_wald_test(var.B, var.V, R, q, dfu,
                                  opt, prn);
        gretl_matrix_reuse(var.B, r, c);
    }

    return err;
}

int johansen_test_complete(GretlVar * jvar, const(DatasetRaw) * dset, Gretlopt opt, Prn * prn)
{
    alias JfunFn = extern (C) int function(GretlVar *, const(DatasetRaw) *, Gretlopt, Prn *);
    JfunFn jfun;
    int err = 0;

    gretl_error_clear();

    jfun = cast(JfunFn) get_plugin_function("johansen_coint_test");

    if (jfun == null) {
        err = 1;
    } else {
        err = jfun(jvar, dset, opt, prn);
    }

    return err;
}

int johansen_estimate_complete(GretlVar * jvar, GretlRestriction * rset, const(DatasetRaw) * dset, Prn * prn)
{
    alias JfunFn = extern (C) int function(GretlVar *, GretlRestriction *, const(DatasetRaw) *, Prn *);
    JfunFn jfun;
    int err = 0;

    gretl_error_clear();

    jfun = cast(JfunFn) get_plugin_function("johansen_estimate");

    if (jfun == null) {
        err = 1;
    } else {
        err = jfun(jvar, rset, dset, prn);
    }

    return err;
}

/* N.B. we allow for the possibility that this allocation has
   already been done. */

int allocate_johansen_extra_matrices(GretlVar * v)
{
    if (v.jinfo.R0 != null &&
        v.jinfo.S00 != null &&
        v.jinfo.YY != null) {
        return 0;
    } else {
        int xc = v.X == null ? 0 : v.X.cols;
        int p0 = v.neqns;
        int p1 = p0 + n_restricted_terms(v);
        int p = p0 + p1;

        clear_gretl_matrix_err();

        if (v.jinfo.R0 == null) {
            v.jinfo.R0 = gretl_matrix_alloc(v.T, p0);
            v.jinfo.R1 = gretl_matrix_alloc(v.T, p1);
        }

        if (v.jinfo.S00 == null) {
            v.jinfo.S00 = gretl_matrix_alloc(p0, p0);
            v.jinfo.S11 = gretl_matrix_alloc(p1, p1);
            v.jinfo.S01 = gretl_matrix_alloc(p0, p1);
        }

        if (xc > 0 && v.ncoeff > 0 && v.jinfo.YY == null) {
            v.jinfo.YY = gretl_matrix_alloc(v.T, p);
            v.jinfo.RR = gretl_matrix_alloc(v.T, p);
            v.jinfo.BB = gretl_matrix_alloc(v.X.cols, p);
        }

        return get_gretl_matrix_err();
    }
}

void johansen_fill_S_matrices(GretlVar * v)
{
    gretl_matrix_multiply_mod(v.jinfo.R0, GretlMatrixMod.transpose,
                              v.jinfo.R0, GretlMatrixMod.none,
                              v.jinfo.S00, GretlMatrixMod.none);
    gretl_matrix_multiply_mod(v.jinfo.R1, GretlMatrixMod.transpose,
                              v.jinfo.R1, GretlMatrixMod.none,
                              v.jinfo.S11, GretlMatrixMod.none);
    gretl_matrix_multiply_mod(v.jinfo.R0, GretlMatrixMod.transpose,
                              v.jinfo.R1, GretlMatrixMod.none,
                              v.jinfo.S01, GretlMatrixMod.none);

    gretl_matrix_divide_by_scalar(v.jinfo.S00, v.T);
    gretl_matrix_divide_by_scalar(v.jinfo.S11, v.T);
    gretl_matrix_divide_by_scalar(v.jinfo.S01, v.T);
}

/* We come here if there are no regressors at stage 1 of the Johansen
   procedure. This happens if the associated VAR is of order 1 (which
   becomes order 0 on differencing) and there are no unrestricted
   exogenous variables.  The residuals are then just the values of the
   variables themselves.
*/

int johansen_degenerate_stage_1(GretlVar * v, const(DatasetRaw) * dset)
{
    const(double) ** Z = cast(const(double) **) dset.Z;
    GretlMatrixRaw * R0 = v.jinfo.R0;
    GretlMatrixRaw * R1 = v.jinfo.R1;
    int i, vi, s, t, j = 0;

    for (i=0; i<v.neqns; i++) {
        vi = v.ylist[i+1];
        s = 0;
        for (t=v.t1; t<=v.t2; t++) {
            gretl_matrix_set(R0, s, j, Z[vi][t] - Z[vi][t-1]);
            gretl_matrix_set(R1, s, j, Z[vi][t-1]);
            s++;
        }
        j++;
    }

    if (auto_restr(v)) {
        int trend = (jcode(v) == JohansenCode.rest_trend);

        for (t=0; t<v.T; t++) {
            gretl_matrix_set(R1, t, j, (trend)? v.t1 + t : 1);
        }
        j++;
    }

    if (v.rlist != null) {
        for (i=0; i<v.rlist[0]; i++) {
            vi = v.rlist[i+1];
            s = 0;
            for (t=v.t1; t<=v.t2; t++) {
                gretl_matrix_set(R1, s++, j, Z[vi][t]); /* was t-1 */
            }
            j++;
        }
    }

    johansen_fill_S_matrices(v);

    return 0;
}

void print_stage_1_coeffs(GretlVar * v, GretlMatrixRaw * B, Prn * prn)
{
    GretlMatrixRaw tmp;

    gretl_matrix_init(&tmp);

    tmp.rows = B.rows;
    tmp.cols = v.neqns;
    tmp.val = B.val;

    gretl_matrix_print_to_prn(&tmp, "\nCoefficients, VAR in differences",
                              prn);

    tmp.cols += n_restricted_terms(v);
    tmp.val += v.neqns * tmp.rows;

    gretl_matrix_print_to_prn(&tmp, "Coefficients, eqns in lagged levels",
                              prn);
}

void johansen_partition_residuals(GretlVar * v)
{
    int n = v.neqns * v.T;
    int m = n + n_restricted_terms(v) * v.T;
    GretlMatrixRaw * R = v.jinfo.RR;

    memcpy(v.jinfo.R0.val, R.val, n * double.sizeof);
    memcpy(v.jinfo.R1.val, R.val + n, m * double.sizeof);
}

/* For Johansen analysis: estimate VAR in differences along with the
   other auxiliary regressions required to compute the relevant matrices
   of residuals, for concentration of the log-likelihood.  Then compute
   S00, S11, S01. See for example James Hamilton, "Time Series
   Analysis", section 20.2.
*/

int johansen_stage_1(GretlVar * v, const(DatasetRaw) * dset, Gretlopt opt, Prn * prn)
{
    int err;

    err = allocate_johansen_extra_matrices(v);
    if (err) {
        fprintf(stderr, "allocate_johansen_extra_matrices: err = %d\n".ptr, err);
        return err;
    }

    if (v.ncoeff == 0 || v.jinfo.BB == null) {
        /* nothing to concentrate out */
        if (opt & Gretlopt.OPT_V) {
            pputs(prn, "\nNo initial VAR estimation is required\n\n");
        }
        johansen_degenerate_stage_1(v, dset);
    } else {
        GretlMatrixRaw * Y = v.jinfo.YY;
        GretlMatrixRaw * B = v.jinfo.BB;
        GretlMatrixRaw * R = v.jinfo.RR;

        VECM_fill_Y(v, dset, Y);

        err = gretl_matrix_multi_SVD_ols(Y, v.X, B, R, null);

        if (!err) {
            if (opt & Gretlopt.OPT_V) {
                print_stage_1_coeffs(v, B, prn);
            }
            johansen_partition_residuals(v);
            johansen_fill_S_matrices(v);
        }
    }

    return err;
}

Gretlopt opt_from_jcode(JohansenCode jc)
{
    Gretlopt opt = Gretlopt.OPT_NONE;

    if (jc == JohansenCode.no_const) {
        opt = Gretlopt.OPT_N;
    } else if (jc == JohansenCode.unrest_trend) {
        opt = Gretlopt.OPT_T;
    } else if (jc == JohansenCode.rest_const) {
        opt =  Gretlopt.OPT_R;
    } else if (jc == JohansenCode.rest_trend) {
        opt = Gretlopt.OPT_A;
    }

    return opt;
}

JohansenCode jcode_from_opt(Gretlopt opt)
{
    JohansenCode jc = JohansenCode.unrest_const;

    if (opt & Gretlopt.OPT_N) {
        jc = JohansenCode.no_const;
    } else if (opt & Gretlopt.OPT_T) {
        jc = JohansenCode.unrest_trend;
    } else if (opt & Gretlopt.OPT_R) {
        jc = JohansenCode.rest_const;
    } else if (opt & Gretlopt.OPT_A) {
        jc = JohansenCode.rest_trend;
    }

    return jc;
}

JohansenInfo * johansen_info_new(GretlVar * var, int rank, Gretlopt opt)
{
    JohansenInfo * jv = cast(JohansenInfo *) malloc(JohansenInfo.sizeof);

    if (jv == null) {
        return null;
    }

    jv.ID = 0;
    jv.code = jcode_from_opt(opt);
    jv.rank = rank;
    jv.seasonals = 0;

    jv.R0 = null;
    jv.R1 = null;
    jv.S00 = null;
    jv.S11 = null;
    jv.S01 = null;
    jv.evals = null;
    jv.Beta = null;
    jv.Alpha = null;
    jv.Gamma = null;
    jv.JC = null;
    jv.JVC = null;
    jv.Cse = null;
    jv.Bse = null;
    jv.Ase = null;
    jv.Bvar = null;
    jv.R = null;
    jv.q = null;
    jv.Ra = null;
    jv.qa = null;
    jv.YY = null;
    jv.RR = null;
    jv.BB = null;
    jv.G = null;

    jv.ll0 = jv.prior_ll = NADBL;
    jv.lrdf = jv.prior_df = 0;

    if (rank > 0) {
        jv.Alpha = gretl_zero_matrix_new(var.neqns, rank);
        if (jv.Alpha == null) {
            free(jv);
            jv = null;
        }
    }

    return jv;
}

void coint_test_header(GretlVar * v, const(DatasetRaw) * dset, Prn * prn)
{
    char[OBSLEN] stobs;
    char[OBSLEN] endobs;

    pprintf(prn, "\n%s:\n".ptr, _("Johansen test").ptr);
    pprintf(prn, "%s = %d\n".ptr, _("Number of equations").ptr, v.neqns);
    pprintf(prn, "%s = %d\n".ptr, _("Lag order").ptr, levels_order(v));
    pprintf(prn, "%s: %s - %s (T = %d)\n".ptr, _("Estimation period").ptr,
            ntolabel(stobs.ptr, v.t1, dset),
            ntolabel(endobs.ptr, v.t2, dset), v.T);

}

int jvar_check_allocation(GretlVar * v, const(DatasetRaw) * dset)
{
    int err = VAR_add_models(v, dset);

    if (!err) {
        err = VAR_allocate_companion_matrix(v);
    }
    if (!err) {
        err = VAR_allocate_cholesky_matrix(v);
    }
    if (!err) {
        err = VAR_allocate_residuals_matrix(v);
    }

    return err;
}

/* Driver function for Johansen analysis.  An appropriately initialized
   "jvar" must have been set up already.
*/

int johansen_driver(GretlVar * jvar, GretlRestriction * rset, const(DatasetRaw) * dset, Gretlopt opt, Prn * prn)
{
    int r = jrank(jvar);

    if (r == 0) {
        /* doing cointegration test */
        coint_test_header(jvar, dset, prn);
    }

    jvar.err = johansen_stage_1(jvar, dset, opt, prn);
    if (jvar.err) {
        return jvar.err;
    }

    if (opt & Gretlopt.OPT_V) {
        print_johansen_sigmas(jvar.jinfo, prn);
    }

    if (r > 0) {
        /* estimating VECM, not just doing cointegration test */
        jvar.err = jvar_check_allocation(jvar, dset);
    }

    if (!jvar.err) {
        /* call the johansen plugin to finish the job */
        if (r > 0) {
            jvar.err = johansen_estimate_complete(jvar, rset, dset,
                                                   prn);
        } else {
            jvar.err = johansen_test_complete(jvar, dset, opt, prn);
        }
    }

    return jvar.err;
}

GretlVar * johansen_wrapper(int code, int order, int rank, const(int) * lags, const(int) * list, GretlRestriction * rset, const(DatasetRaw) * dset, Gretlopt opt, Prn * prn, int * err)
{
    GretlVar * jvar;

    jvar = gretl_VAR_new(code, order - 1, rank, lags, list, dset,
                         opt, err);
    if (jvar != null && !jvar.err) {
        *err = jvar.err = johansen_driver(jvar, rset, dset, opt, prn);
    }

    return jvar;
}

/**
 * johansen_test:
 * @order: lag order for test.
 * @list: list of variables to test for cointegration.
 * @dset: dataset struct.
 * @opt: %OPT_A: include constant plus restricted trend; %OPT_D:
 * include centered seasonals; %OPT_N: no constant; %OPT_R:
 * restricted constant; %OPT_T: constant and unrestricted trend
 * (note: default "case" is unrestricted constant).
 * %OPT_V: produce verbose results.
 * @prn: gretl printing struct.
 *
 * Carries out the Johansen test for cointegration.
 *
 * Returns: pointer to struct containing information on
 * the test.
 */

GretlVar * johansen_test(int order, const(int) * list, const(DatasetRaw) * dset, Gretlopt opt, Prn * prn)
{
    int err = 0;

    return johansen_wrapper(VarTask.vecm_ctest, order, 0, null, list, null,
                            dset, opt, prn, &err);
}

/**
 * johansen_test_simple:
 * @order: lag order for test.
 * @list: list of variables to test for cointegration.
 * @dset: dataset struct.
 * @opt: %OPT_A: include constant plus restricted trend; %OPT_D:
 * include centered seasonals; %OPT_N: no constant; %OPT_R:
 * restricted constant; %OPT_T: constant and unrestricted trend
 * (note: default "case" is unrestricted constant);
 * %OPT_V: produce verbose results; %OPT_Q: just print the tests;
 * %OPT_S: don't print anything.
 * @prn: gretl printing struct.
 *
 * Carries out the Johansen test for cointegration and prints the
 * results (but unlike johansen_test(), does not return the
 * allocated results in VAR form).
 *
 * Returns: 0 on success, non-zero code on error.
 */

int johansen_test_simple(int order, const(int) * list, const(DatasetRaw) * dset, Gretlopt opt, Prn * prn)
{
    GretlVar * jvar = null;
    int err = 0;

    jvar = johansen_wrapper(VarTask.vecm_ctest, order, 0, null, list, null,
                            dset, opt, (opt & Gretlopt.OPT_S)? null : prn,
                            &err);

    if (jvar != null) {
        gretl_VAR_free(jvar);
    }

    return err;
}

/**
 * gretl_VECM:
 * @order: lag order.
 * @rank: cointegration rank.
 * @list: list of endogenous variables, possibly plus
 * exogenous variables.
 * @dset: dataset struct.
 * @opt: may include OPT_N ("nc"), OPT_R ("rc"), OPT_A ("crt"),
 * OPT_T ("ct"), OPT_D (include seasonals), OPT_F (show variance
 * decompositions), OPT_I (show impulse responses), OPT_V
 * (verbose operation), OPT_Q (quiet), OPT_S (silent).
 * @prn: gretl printing struct.
 * @err: location to receive error code.
 *
 * Returns: pointer to struct containing a VECM system,
 * or %NULL on failure.
 */

GretlVar * gretl_VECM(int order, int rank, int * list, const(DatasetRaw) * dset, Gretlopt opt, Prn * prn, int * err)
{
    GretlVar * jvar = null;
    int * lags = null;

    if (rank <= 0) {
        gretl_errmsg_sprintf(_("vecm: rank %d is out of bounds").ptr, rank);
        *err = GretlError.data;
        return null;
    }

    jvar = johansen_wrapper(VarTask.vecm_estimate, order, rank, lags, list,
                            null, dset, opt, prn, err);

    if (jvar != null && !jvar.err) {
        gretl_VAR_print(jvar, dset, opt, prn);
    }

    return jvar;
}

int * VAR_list_composite(const(int) * ylist, const(int) * xlist, const(int) * rlist)
{
    int * big = null;
    int i, k, n = ylist[0];

    if (xlist != null && xlist[0] > 0) {
        n += xlist[0] + 1;
    }

    if (rlist != null && rlist[0] > 0) {
        n += rlist[0] + 1;
        if (xlist == null || xlist[0] == 0) {
            /* extra separator needed */
            n++;
        }
    }

    big = gretl_list_new(n);
    if (big == null) {
        return null;
    }

    k = 1;

    for (i=1; i<=ylist[0]; i++) {
        big[k++] = ylist[i];
    }

    if (xlist != null && xlist[0] > 0) {
        big[k++] = LISTSEP;
        for (i=1; i<=xlist[0]; i++) {
            big[k++] = xlist[i];
        }
    }

    if (rlist != null && rlist[0] > 0) {
        if (xlist == null || xlist[0] == 0) {
            /* placeholder for empty xlist */
            big[k++] = LISTSEP;
        }
        big[k++] = LISTSEP;
        for (i=1; i<=rlist[0]; i++) {
            big[k++] = rlist[i];
        }
    }

    return big;
}

int * rebuild_full_VAR_list(GretlVar * var)
{
    int * list;

    if (var.xlist == null && var.rlist == null) {
        list = gretl_list_copy(var.ylist);
    } else {
        list = VAR_list_composite(var.ylist, var.xlist,
                                  var.rlist);
    }

    return list;
}

/**
 * real_gretl_restricted_vecm:
 * @orig: orginal VECM model.
 * @rset: restriction information.
 * @dset: dataset struct.
 * @prn: gretl printing struct.
 * @err: location to receive error code.
 *
 * Returns: pointer to struct containing information on
 * the newly restricted VECM system or %NULL on failure.
 */

GretlVar * real_gretl_restricted_vecm(GretlVar * orig, GretlRestriction * rset, const(DatasetRaw) * dset, Prn * prn, int * err)
{
    GretlVar * jvar = null;
    Gretlopt jopt = Gretlopt.OPT_NONE;
    int * list = null;

    if (orig == null || orig.jinfo == null || rset == null) {
        *err = GretlError.data;
        return null;
    }

    list = rebuild_full_VAR_list(orig);
    if (list == null) {
        *err = GretlError.alloc;
        return null;
    }

    jopt |= opt_from_jcode(orig.jinfo.code);

    if (orig.jinfo.seasonals > 0) {
        jopt |= Gretlopt.OPT_D;
    }

    jvar = johansen_wrapper(VarTask.vecm_estimate, orig.order + 1,
                            orig.jinfo.rank, orig.lags,
                            list, rset, dset,
                            jopt, prn, err);

    if (jvar != null && !jvar.err) {
        Gretlopt ropt, prnopt = Gretlopt.OPT_NONE;
        int df;

        df = jvar.jinfo.lrdf - orig.jinfo.lrdf;

        if (df > 0) {
            double x = 2 * (orig.ll - jvar.ll);
            double pv = chisq_cdf_comp(df, x);

            rset_add_results(rset, x, pv, jvar.ll);
            rset_record_LR_result(rset);
        }

        jvar.jinfo.prior_ll = orig.ll;
        jvar.jinfo.prior_df = orig.jinfo.lrdf;

        ropt = gretl_restriction_get_options(rset);
        if (ropt & Gretlopt.OPT_Q) {
            prnopt = Gretlopt.OPT_Q;
        }

        if (!(ropt & Gretlopt.OPT_S)) {
            /* FIXME OPT_I, OPT_F, impulses and decomp? */
            gretl_VAR_print(jvar, dset, prnopt, prn);
        }
    }

    free(list);

    return jvar;
}

void gretl_VAR_set_name(GretlVar * var, const(char) * name)
{
    if (name == var.name) {
        return;
    }

    if (var.name == null) {
        var.name = cast(typeof(var.name)) malloc(MAXSAVENAME);
    }

    if (var.name != null) {
        *var.name = '\0';
        strncat(var.name, name, MAXSAVENAME - 1);
    }
}

const(char) * gretl_VAR_get_name(GretlVar * var)
{
    return var.name;
}

double * gretl_VAR_get_resid_series(GretlVar * var, int eqnum, int * err)
{
    Model * pmod;
    double * u = null;

    if (var.models == null || eqnum < 0 || eqnum >= var.neqns) {
        *err = GretlError.badstat;
        return null;
    }

    pmod = var.models[eqnum];
    u = copyvec(pmod.uhat, pmod.full_n);

    if (u == null) {
        *err = GretlError.alloc;
    }

    return u;
}

int gretl_VAR_set_ordering(GretlVar * var, GretlMatrixRaw * ord)
{
    gretl_matrix_free(var.ord);
    var.ord = ord;
    return 0;
}

int gretl_VAR_do_irf(GretlVar * var, const(char) * line, const(DatasetRaw) * dset)
{
    int targ = -1, shock = 1;
    int h = 20, boot = 0;
    double alpha = 0.10;
    int err = 0;
    const(char) * s;

    s = strstr(line, "--targ=");
    if (s != null) {
        targ = atoi(s + 7) - 1;
    }

    s = strstr(line, "--shock=");
    if (s != null) {
        shock = atoi(s + 8) - 1;
    }

    s = strstr(line, "--horizon=");
    if (s != null) {
        h = atoi(s + 10);
    }

    s = strstr(line, "--alpha=");
    if (s != null) {
        alpha = dot_atof(s + 8);
    }

    if (strstr(line, "--bootstrap") != null) {
        boot = 1;
    }

    if (targ < 0 || shock < 0 || h <= 0 ||
        alpha < .01 || alpha > 0.5) {
        err = GretlError.invarg;
    } else if (boot) {
        err = gretl_VAR_plot_impulse_response(var, targ, shock,
                                              h, alpha, dset,
                                              Gretlopt.OPT_NONE);
    } else {
        err = gretl_VAR_plot_impulse_response(var, targ, shock,
                                              h, 0, dset,
                                              Gretlopt.OPT_NONE);
    }

    return err;
}

int gretl_VAR_get_highest_variable(GretlVar * var)
{
    int i, vi, vmax = 0;

    if (var.ylist != null) {
        for (i=1; i<=var.ylist[0]; i++) {
            vi = var.ylist[i];
            if (vi > vmax) {
                vmax = vi;
            }
        }
    }

    if (var.xlist != null) {
        for (i=1; i<=var.xlist[0]; i++) {
            vi = var.xlist[i];
            if (vi > vmax) {
                vmax = vi;
            }
        }
    }

    if (var.rlist != null) {
        for (i=1; i<=var.rlist[0]; i++) {
            vi = var.rlist[i];
            if (vi > vmax) {
                vmax = vi;
            }
        }
    }

    return vmax;
}

int gretl_VECM_id(GretlVar * vecm)
{
    __gshared int nvecm = 0;

    if (vecm.jinfo.ID == 0) {
        vecm.jinfo.ID = ++nvecm;
    }

    return vecm.jinfo.ID;
}

int gretl_VECM_n_beta(GretlVar * vecm)
{
    int nb = 0;

    if (vecm.jinfo != null && vecm.jinfo.Beta != null) {
        nb = gretl_matrix_rows(vecm.jinfo.Beta);
    }

    return nb;
}

int gretl_VECM_n_alpha(GretlVar * vecm)
{
    int na = 0;

    if (vecm.jinfo != null && vecm.jinfo.Alpha != null) {
        na = gretl_matrix_rows(vecm.jinfo.Alpha);
    }

    return na;
}

int gretl_VECM_rank(GretlVar * vecm)
{
    int r = 0;

    if (vecm.jinfo != null) {
        r = vecm.jinfo.rank;
    }

    return r;
}

int beta_restricted_VECM(GretlVar * vecm)
{
    if (vecm.jinfo != null && vecm.jinfo.R != null) {
        return 1;
    }

    return 0;
}

int alpha_restricted_VECM(GretlVar * vecm)
{
    if (vecm.jinfo != null && vecm.jinfo.Ra != null) {
        return 1;
    }

    return 0;
}

int restricted_VECM(GretlVar * vecm)
{
    if (vecm.jinfo != null &&
        (vecm.jinfo.R != null || vecm.jinfo.Ra != null)) {
        return 1;
    }

    return 0;
}

const(GretlMatrixRaw) * gretl_VECM_R_matrix(GretlVar * vecm)
{
    return (vecm.jinfo != null)? vecm.jinfo.R : null;
}

const(GretlMatrixRaw) * gretl_VECM_q_matrix(GretlVar * vecm)
{
    return (vecm.jinfo != null)? vecm.jinfo.q : null;
}

const(GretlMatrixRaw) * gretl_VECM_Ra_matrix(GretlVar * vecm)
{
    return (vecm.jinfo != null)? vecm.jinfo.Ra : null;
}

const(GretlMatrixRaw) * gretl_VECM_qa_matrix(GretlVar * vecm)
{
    return (vecm.jinfo != null)? vecm.jinfo.qa : null;
}

double * gretl_VAR_get_series(GretlVar * var, const(DatasetRaw) * dset, int idx, const(char) * key, int * err)
{
    double * x = null;
    const(char) * msel;
    int t, col = 0;

    if (var == null || idx != ModelDataIndex.uhat) { /* FIXME generalize this */
        *err = GretlError.badstat;
        return null;
    }

    msel = strchr(key, '[');
    if (msel == null || sscanf(msel, "[,%d]".ptr, &col) != 1) {
        *err = GretlError.parse;
    } else if (col <= 0 || col > var.neqns) {
        *err = GretlError.data;
    }

    if (!*err) {
        x = cast(typeof(x)) malloc(dset.n * (*x).sizeof);
        if (x == null) {
            *err = GretlError.alloc;
        }
    }

    if (!*err) {
        const(Model) * pmod = var.models[col-1];

        if (pmod == null || pmod.full_n != dset.n) {
            *err = GretlError.data;
            free(x);
            x = null;
        } else {
            for (t=0; t<dset.n; t++) {
                x[t] = pmod.uhat[t];
            }
        }
    }

    return x;
}

/* get coefficients or standard errors from the models
   in the VAR */

GretlMatrixRaw * VAR_matrix_from_models(GretlVar * var, int idx, int * err)
{
    GretlMatrixRaw * m = null;
    const(Model) * pmod;
    double x;
    int i, j;

    if (idx == ModelDataIndex.vecg) {
        if (var.ci != GretlCmdIndex.vecm || var.order == 0) {
            *err = GretlError.badstat;
            return null;
        }
    }

    if (idx == ModelDataIndex.vecg) {
        m = gretl_matrix_alloc(var.neqns, var.neqns * var.order);
    } else {
        m = gretl_matrix_alloc(var.models[0].ncoeff, var.neqns);
    }

    if (m == null) {
        *err = GretlError.alloc;
        return null;
    }

    if (idx == ModelDataIndex.vecg) {
        int k, mcol, cnum;

        for (j=0; j<var.neqns; j++) {
            pmod = var.models[j];
            mcol = 0;
            for (i=0; i<var.order; i++) {
                cnum = pmod.ifc + i;
                for (k=0; k<var.neqns; k++) {
                    x = pmod.coeff[cnum];
                    gretl_matrix_set(m, j, mcol++, x);
                    cnum += var.order;
                }
            }
        }
    } else {
        for (j=0; j<var.neqns; j++) {
            pmod = var.models[j];
            for (i=0; i<pmod.ncoeff; i++) {
                if (idx == ModelDataIndex.coeff) {
                    x = pmod.coeff[i];
                } else {
                    x = pmod.sderr[i];
                }
                gretl_matrix_set(m, i, j, x);
            }
        }
    }

    return m;
}

GretlMatrixRaw * alt_VECM_get_EC_matrix(GretlVar * vecm, const(DatasetRaw) * dset, int * err)
{
    const(GretlMatrixRaw) * B = vecm.jinfo.Beta;
    int r = jrank(vecm);
    GretlMatrixRaw * EC = null;
    double xti, xj;
    int s, t, T;
    int i, j, k;

    if (dset == null || dset.Z == null) {
        *err = GretlError.badstat;
        return null;
    }

    for (i=1; i<=vecm.ylist[0]; i++) {
        if (vecm.ylist[i] >= dset.v) {
            *err = GretlError.data;
            return null;
        }
    }

    T = vecm.t2 - vecm.t1 + 1;

    EC = gretl_matrix_alloc(T, r);
    if (EC == null) {
        *err = GretlError.alloc;
        return null;
    }

    s = 0;

    for (t=vecm.t1; t<=vecm.t2; t++) {
        for (j=0; j<r; j++) {
            xj = 0.0;
            /* beta * X(t-1) */
            k = 0;
            for (i=0; i<vecm.neqns; i++) {
                xti = dset.Z[vecm.ylist[i+1]][t-1];
                if (na(xti)) {
                    xj = NADBL;
                    break;
                }
                xj += xti * gretl_matrix_get(B, k++, j);
            }

            /* restricted const or trend */
            if (auto_restr(vecm) && !na(xj)) {
                xti = gretl_matrix_get(B, k++, j);
                if (jcode(vecm) == JohansenCode.rest_trend) {
                    xti *= t;
                }
                xj += xti;
            }

            /* restricted exog vars */
            if (vecm.rlist != null && !na(xj)) {
                for (i=0; i<vecm.rlist[0]; i++) {
                    xti = dset.Z[vecm.rlist[i+1]][t-1];
                    if (na(xti)) {
                        xj = NADBL;
                        break;
                    }
                    xj += xti * gretl_matrix_get(B, k++, j);
                }
            }

            if (na(xj)) {
                gretl_matrix_set(EC, s, j, NADBL);
            } else {
                gretl_matrix_set(EC, s, j, xj);
            }
        }
        s++;
    }

    gretl_matrix_set_t1(EC, vecm.t1);
    gretl_matrix_set_t2(EC, vecm.t2);

    return EC;
}

GretlMatrixRaw * VECM_get_EC_matrix(GretlVar * v, const(DatasetRaw) * dset, int * err)
{
    GretlMatrixRaw * EC = null;
    double x;
    int rank, k, k0;
    int j, t, T;

    rank = jrank(v);
    if (rank == 0) {
        *err = GretlError.badstat;
        return null;
    }

    if (v.X == null) {
        fprintf(stderr, "VECM_get_EC_matrix: v->X is NULL\n".ptr);
        *err = GretlError.badstat;
        return null;
    } else if (v.X.cols < v.ncoeff) {
        fprintf(stderr, "VECM_get_EC_matrix: v->X is short of cols\n".ptr);
        return alt_VECM_get_EC_matrix(v, dset, err);
    }

    T = v.X.rows;

    EC = gretl_matrix_alloc(T, rank);
    if (EC == null) {
        *err = GretlError.alloc;
        return null;
    }

    k0 = v.ncoeff - rank;

    for (j=0, k=k0; j<rank; j++, k++) {
        for (t=0; t<T; t++) {
            x = gretl_matrix_get(v.X, t, k);
            gretl_matrix_set(EC, t, j, x);
        }
    }

    gretl_matrix_set_t1(EC, v.t1);
    gretl_matrix_set_t2(EC, v.t2);

    return EC;
}

bool vecm_matrix(int i)
{
    return (i == ModelDataIndex.jalpha || i == ModelDataIndex.jbeta || i == ModelDataIndex.jvbeta || i == ModelDataIndex.js00 || i == ModelDataIndex.js11 || i == ModelDataIndex.js01 || i == ModelDataIndex.ec || i == ModelDataIndex.evals);
}

GretlMatrixRaw * gretl_VAR_get_matrix(GretlVar * var, int idx, int * err)
{
    const(GretlMatrixRaw) * src = null;
    GretlMatrixRaw * M = null;
    int copy = 1;

    if (var == null) {
        *err = GretlError.badstat;
        return null;
    }

    if (idx == ModelDataIndex.uhat) {
        src = gretl_VAR_get_residual_matrix(var);
    } else if (idx == ModelDataIndex.yhat) {
        M = gretl_VAR_get_fitted_matrix(var);
        copy = 0;
    } else if (idx == ModelDataIndex.compan) {
        src = var.A;
    } else if (idx == ModelDataIndex.coeff || idx == ModelDataIndex.se || idx == ModelDataIndex.vecg) {
        M = VAR_matrix_from_models(var, idx, err);
        copy = 0;
    } else if (idx == ModelDataIndex.xtxinv) {
        src = var.XTX;
    } else if (idx == ModelDataIndex.sigma) {
        src = var.S;
    } else if (idx == ModelDataIndex.pmanteau) {
        M = gretl_VAR_get_portmanteau_test(var, err);
        copy = 0;
    } else if (vecm_matrix(idx)) {
        if (var.jinfo != null) {
            switch (idx) {
            default:
                break;
            case ModelDataIndex.evals:
                src = var.jinfo.evals;
                break;
            case ModelDataIndex.jalpha:
                src = var.jinfo.Alpha;
                break;
            case ModelDataIndex.jbeta:
                src = var.jinfo.Beta;
                break;
            case ModelDataIndex.jvbeta:
                src = var.jinfo.Bvar;
                break;
            case ModelDataIndex.js00:
                src = var.jinfo.S00;
                break;
            case ModelDataIndex.js11:
                src = var.jinfo.S11;
                break;
            case ModelDataIndex.js01:
                src = var.jinfo.S01;
                break;
            case ModelDataIndex.ec:
                M = VECM_get_EC_matrix(var, null, err);
                copy = 0;
                break;
            }
        }
    }

    if (copy) {
        if (src == null) {
            *err = GretlError.badstat;
        } else {
            M = gretl_matrix_copy(src);
            if (M == null) {
                *err = GretlError.alloc;
            }
        }
    }

    return M;
}

/* retrieve EC (j >= 0 && j < rank) as a full-length series */

double * gretl_VECM_get_EC(GretlVar * vecm, int j, const(DatasetRaw) * dset, int * err)
{
    const(GretlMatrixRaw) * B = vecm.jinfo.Beta;
    int r = jrank(vecm);
    double * x = null;
    double xti;
    int i, k, t, t0;

    if (j < 0 || j >= r) {
        *err = GretlError.data;
        return null;
    }

    for (i=1; i<=vecm.ylist[0]; i++) {
        if (vecm.ylist[i] >= dset.v) {
            *err = GretlError.data;
            return null;
        }
    }

    x = cast(typeof(x)) malloc(dset.n * (*x).sizeof);
    if (x == null) {
        *err = GretlError.alloc;
        return null;
    }

    t0 = (dset.t1 >= 1)? dset.t1 : 1;

    for (t=0; t<dset.n; t++) {
        if (t < t0 || t > dset.t2) {
            x[t] = NADBL;
            continue;
        }
        x[t] = 0.0;

        k = 0;

        /* beta * X(t-1) */
        for (i=0; i<vecm.neqns; i++) {
            xti = dset.Z[vecm.ylist[i+1]][t-1];
            if (na(xti)) {
                x[t] = NADBL;
                break;
            }
            x[t] += xti * gretl_matrix_get(B, k++, j);
        }

        /* restricted const or trend */
        if (auto_restr(vecm) && !na(x[t])) {
            xti = gretl_matrix_get(B, k++, j);
            if (jcode(vecm) == JohansenCode.rest_trend) {
                xti *= t;
            }
            x[t] += xti;
        }

        /* restricted exog vars */
        if (vecm.rlist != null && !na(x[t])) {
            for (i=0; i<vecm.rlist[0]; i++) {
                xti = dset.Z[vecm.rlist[i+1]][t-1];
                if (na(xti)) {
                    x[t] = NADBL;
                    break;
                }
                x[t] += xti * gretl_matrix_get(B, k++, j);
            }
        }
    }

    return x;
}

GretlVar * gretl_VAR_rebuilder_new()
{
    GretlVar * var = cast(GretlVar *) malloc(GretlVar.sizeof);

    if (var == null) {
        return null;
    }

    gretl_VAR_clear(var);

    return var;
}

int VAR_retrieve_jinfo(xmlNodePtr node, xmlDocPtr doc, GretlVar * var)
{
    xmlNodePtr cur = node.children;
    JohansenInfo *jinfo;
    int ID, code, rank;
    int seas;
    int got = 0;
    int err = 0;

    got += gretl_xml_get_prop_as_int(node, "ID", &ID);
    got += gretl_xml_get_prop_as_int(node, "code", &code);
    got += gretl_xml_get_prop_as_int(node, "rank", &rank);
    got += gretl_xml_get_prop_as_int(node, "seasonals", &seas);

    if (got != 4) {
        return GretlError.data;
    }

    jinfo = johansen_info_new(var, rank, Gretlopt.OPT_NONE);
    if (jinfo == null) {
        return GretlError.alloc;
    }

    jinfo.ID = ID;
    jinfo.code = cast(JohansenCode) code;
    jinfo.rank = rank;
    jinfo.seasonals = seas;

    gretl_xml_get_prop_as_double(node, "ll0", &jinfo.ll0);
    gretl_xml_get_prop_as_int(node, "bdf", &jinfo.lrdf);
    gretl_xml_get_prop_as_double(node, "oldll", &jinfo.prior_ll);
    gretl_xml_get_prop_as_int(node, "olddf", &jinfo.prior_df);

    cur = node.children;

    while (cur != null && !err) {
        if (!xmlStrcmp(cur.name, cast(XUC) "gretl-matrix")) {
            char * mname;

            gretl_xml_get_prop_as_string(cur, "name", &mname);
            if (mname == null) {
                err = GretlError.data;
            } else {
                if (!strcmp(mname, "u")) {
                    jinfo.R0 = gretl_xml_get_matrix(cur, doc, &err);
                } else if (!strcmp(mname, "v")) {
                    jinfo.R1 = gretl_xml_get_matrix(cur, doc, &err);
                } else if (!strcmp(mname, "Suu")) {
                    jinfo.S00 = gretl_xml_get_matrix(cur, doc, &err);
                } else if (!strcmp(mname, "Svv")) {
                    jinfo.S11 = gretl_xml_get_matrix(cur, doc, &err);
                } else if (!strcmp(mname, "Suv")) {
                    jinfo.S01 = gretl_xml_get_matrix(cur, doc, &err);
                } else if (!strcmp(mname, "evals")) {
                    jinfo.evals = gretl_xml_get_matrix(cur, doc, &err);
                } else if (!strcmp(mname, "Beta")) {
                    jinfo.Beta = gretl_xml_get_matrix(cur, doc, &err);
                } else if (!strcmp(mname, "Alpha")) {
                    jinfo.Alpha = gretl_xml_get_matrix(cur, doc, &err);
                } else if (!strcmp(mname, "Bvar")) {
                    jinfo.Bvar = gretl_xml_get_matrix(cur, doc, &err);
                } else if (!strcmp(mname, "Bse")) {
                    jinfo.Bse = gretl_xml_get_matrix(cur, doc, &err);
                } else if (!strcmp(mname, "R")) {
                    jinfo.R = gretl_xml_get_matrix(cur, doc, &err);
                } else if (!strcmp(mname, "q")) {
                    jinfo.q = gretl_xml_get_matrix(cur, doc, &err);
                } else if (!strcmp(mname, "Ra")) {
                    jinfo.Ra = gretl_xml_get_matrix(cur, doc, &err);
                } else if (!strcmp(mname, "qa")) {
                    jinfo.qa = gretl_xml_get_matrix(cur, doc, &err);
                }
                free(mname);
            }
        }
        cur = cur.next;
    }

    if (err) {
        johansen_info_free(var);
    } else {
        var.jinfo = jinfo;
    }

    fprintf(stderr, "VAR_retrieve_jinfo: err = %d\n".ptr, err);

    return err;
}

int VAR_retrieve_equations(xmlNodePtr node, xmlDocPtr doc, Model ** models, int neqns, const(DatasetRaw) * dset)
{
    xmlNodePtr cur = node.children;
    int i = 0, err = 0;

    while (cur != null && !err) {
        Model * pmod;

        if (!xmlStrcmp(cur.name, cast(XUC) "gretl-model")) {
            pmod = gretl_model_from_XML(cur, doc, dset, &err);
            if (pmod != null) {
                models[i++] = pmod;
            }
        }
        cur = cur.next;
    }

    fprintf(stderr, "VAR_retrieve_equations: got %d models (neqns = %d), err = %d\n".ptr,
            i, neqns, err);

    return err;
}

/* for backward compatibility */

int make_VAR_global_lists(GretlVar * var)
{
    Model * pmod = var.models[0];
    int * list = pmod.list;
    int p = var.order;
    int n = var.neqns;
    int nx, np = n * p;
    int i, j, ifc;

    if (list == null || list[0] < 3) {
        return GretlError.data;
    }

    /* the _endogenous_ vars start in position 2 or 3 (3 if a constant
       is included), and there are (order * neqns) such terms */

    ifc = (list[2] == 0);
    nx = list[0] - 1 - ifc - np;

    var.ylist = gretl_list_new(n);
    if (var.ylist == null) {
        return GretlError.alloc;
    }

    if (nx > 0) {
        var.xlist = gretl_list_new(nx);
        if (var.xlist == null) {
            free(var.ylist);
            var.ylist = null;
            return GretlError.alloc;
        }
    }

    for (i=0; i<n; i++) {
        pmod = var.models[i];
        var.ylist[i+1] = pmod.list[1];
    }

    j = 2 + ifc + np;
    for (i=0; i<nx; i++) {
        var.xlist[i+1] = list[j++];
    }

    return 0;
}

int rebuild_VAR_matrices(GretlVar * var)
{
    Model * pmod;
    double x;
    int gotA = (var.A != null);
    int gotX = (var.X != null);
    int j, i;
    int err = 0;

    if (var.E == null) {
        err = VAR_allocate_residuals_matrix(var);
    }
    if (!err && var.A == null) {
        err = VAR_allocate_companion_matrix(var);
    }
    if (!err && gotX && var.XTX == null) {
        err = VAR_add_XTX_matrix(var);
    }
    if (!err && var.C == null) {
        err = VAR_allocate_cholesky_matrix(var);
    }
    if (!err && var.B == null) {
        var.B = gretl_matrix_alloc(var.models[0].ncoeff,
                                    var.neqns);
        if (var.B == null) {
            err = GretlError.alloc;
        }
    }

    for (j=0; j<var.neqns && !err; j++) {
        pmod = var.models[j];
        for (i=0; i<pmod.ncoeff; i++) {
            x = pmod.coeff[i];
            gretl_matrix_set(var.B, i, j, x);
        }
        for (i=0; i<var.T; i++) {
            x = pmod.uhat[pmod.t1 + i];
            gretl_matrix_set(var.E, i, j, x);
        }
    }

    if (!err && !gotA) {
        /* note: for VECMs, A should be retrieved from the session
           file, and gotA should be non-zero
        */
        VAR_write_A_matrix(var, GretlMatrixMod.none);
    }

    if (!err && !gotX) {
        fprintf(stderr, "Can't we rebuild VAR->X somehow?\n".ptr);
    }

    return err;
}

GretlVar * gretl_VAR_from_XML(xmlNodePtr node, xmlDocPtr doc, const(DatasetRaw) * dset, int * err)
{
    GretlVar * var;
    Model * pmod;
    xmlNodePtr cur;
    char * vname;
    int i, n, got = 0;

    var = gretl_VAR_rebuilder_new();
    if (var == null) {
        *err = GretlError.alloc;
        return null;
    }

    got += gretl_xml_get_prop_as_int(node, "ecm", &var.ci);
    got += gretl_xml_get_prop_as_int(node, "neqns", &var.neqns);
    got += gretl_xml_get_prop_as_int(node, "order", &var.order);

    if (got < 3) {
        *err = GretlError.data;
        goto bailout;
    }

    var.ci = (var.ci == 0)? GretlCmdIndex.var : GretlCmdIndex.vecm;

    gretl_xml_get_prop_as_string(node, "name", &vname);
    if (vname != null) {
        gretl_VAR_set_name(var, vname);
        free(vname);
    }

    /* these are not show-stoppers */
    gretl_xml_get_prop_as_int(node, "robust", &var.robust);
    gretl_xml_get_prop_as_int(node, "detflags", &var.detflags);
    gretl_xml_get_prop_as_int(node, "LBs", &var.LBs);
    gretl_xml_get_prop_as_double(node, "LB", &var.LB);

    var.models = cast(typeof(var.models)) malloc(var.neqns * (*var.models).sizeof);
    if (var.models == null) {
        *err = GretlError.alloc;
        goto bailout;
    }

    for (i=0; i<var.neqns; i++) {
        var.models[i] = null;
    }

    gretl_push_c_numeric_locale();

    cur = node.children;

    while (cur != null && !*err) {
        if (!xmlStrcmp(cur.name, cast(XUC) "lags")) {
            var.lags = gretl_xml_get_list(cur, doc, err);
        } else if (!xmlStrcmp(cur.name, cast(XUC) "ylist")) {
            var.ylist = gretl_xml_get_list(cur, doc, err);
        } else if (!xmlStrcmp(cur.name, cast(XUC) "xlist")) {
            var.xlist = gretl_xml_get_list(cur, doc, err);
        } else if (!xmlStrcmp(cur.name, cast(XUC) "rlist")) {
            var.rlist = gretl_xml_get_list(cur, doc, err);
        } else if (!xmlStrcmp(cur.name, cast(XUC) "Fvals")) {
            var.Fvals = gretl_xml_get_double_array(cur, doc, &n, err);
        } else if (!xmlStrcmp(cur.name, cast(XUC) "Ivals")) {
            var.Ivals = gretl_xml_get_double_array(cur, doc, &n, err);
        } else if (!xmlStrcmp(cur.name, cast(XUC) "equations")) {
            *err = VAR_retrieve_equations(cur, doc, var.models, var.neqns, dset);
        } else if (!xmlStrcmp(cur.name, cast(XUC) "gretl-johansen")) {
            *err = VAR_retrieve_jinfo(cur, doc, var);
        } else if (!xmlStrcmp(cur.name, cast(XUC) "gretl-matrix")) {
            char * mname;

            gretl_xml_get_prop_as_string(cur, "name", &mname);
            if (mname == null) {
                *err = GretlError.data;
            } else {
                if (!strcmp(mname, "A")) {
                    var.A = gretl_xml_get_matrix(cur, doc, err);
                } else if (!strcmp(mname, "X")) {
                    var.X = gretl_xml_get_matrix(cur, doc, err);
                } else if (!strcmp(mname, "Y")) {
                    var.Y = gretl_xml_get_matrix(cur, doc, err);
                } else if (!strcmp(mname, "ord")) {
                    var.ord = gretl_xml_get_matrix(cur, doc, err);
                }
                free(mname);
            }
        }
        cur = cur.next;
    }

    gretl_pop_c_numeric_locale();

    if (*err) {
        goto bailout;
    }

    pmod = var.models[0];

    /* Note: get "global" info from the first model (equation) */
    var.ncoeff = pmod.ncoeff;
    var.ifc = pmod.ifc;
    var.t1 = pmod.t1;
    var.t2 = pmod.t2;
    var.T = var.t2 - var.t1 + 1;
    var.df = var.T - var.ncoeff;

    if (var.ylist == null) {
        *err = make_VAR_global_lists(var);
    }
    if (!*err) {
        *err = rebuild_VAR_matrices(var);
    }
    if (!*err) {
        *err = VAR_add_stats(var, 0);
    }
    if (!*err) {
        *err = gretl_VAR_do_error_decomp(var.S, var.C, null);
    }

    /* FIXME vecm tests on beta/alpha */

 bailout:

    if (*err) {
        gretl_VAR_free(var);
        var = null;
    }

    return var;
}

void johansen_serialize(JohansenInfo * j, Prn * prn)
{
    pprintf(prn, "<gretl-johansen ID=\"%d\" code=\"%d\" rank=\"%d\" ".ptr,
            j.ID, j.code, j.rank);
    pprintf(prn, "seasonals=\"%d\" ".ptr, j.seasonals);

    if (j.lrdf > 0 && !na(j.ll0)) {
        gretl_xml_put_double("ll0", j.ll0, prn);
        gretl_xml_put_int("bdf", j.lrdf, prn);
    }

    if (j.prior_df > 0 && !na(j.prior_ll)) {
        gretl_xml_put_double("oldll", j.prior_ll, prn);
        gretl_xml_put_int("olddf", j.prior_df, prn);
    }

    pputs(prn, ">\n");

    gretl_matrix_serialize(j.R0, "u", prn);
    gretl_matrix_serialize(j.R1, "v", prn);
    gretl_matrix_serialize(j.S00, "Suu", prn);
    gretl_matrix_serialize(j.S11, "Svv", prn);
    gretl_matrix_serialize(j.S01, "Suv", prn);
    gretl_matrix_serialize(j.evals, "evals", prn);
    gretl_matrix_serialize(j.Beta, "Beta", prn);
    gretl_matrix_serialize(j.Alpha, "Alpha", prn);
    gretl_matrix_serialize(j.Bvar, "Bvar", prn);
    gretl_matrix_serialize(j.Bse, "Bse", prn);
    gretl_matrix_serialize(j.Gamma, "Gamma", prn);
    gretl_matrix_serialize(j.JC, "long_run", prn);
    gretl_matrix_serialize(j.JVC, "cov_long_run", prn);
    gretl_matrix_serialize(j.Cse, "se_long_run", prn);
    gretl_matrix_serialize(j.R, "R", prn);
    gretl_matrix_serialize(j.q, "q", prn);
    gretl_matrix_serialize(j.Ra, "Ra", prn);
    gretl_matrix_serialize(j.qa, "qa", prn);

    pputs(prn, "</gretl-johansen>\n");
}

/* Retrieve enough VECM-related info from @b to carry
   out an IRF bootstrap: this requires more testing.
*/

int retrieve_johansen_basics(GretlVar * var, GretlBundle * b)
{
    int err = 0;

    var.jinfo = cast(typeof(var.jinfo)) calloc(1, (*var.jinfo).sizeof);

    if (var.jinfo == null) {
        err = GretlError.alloc;
    } else {
        int i;
        int[10] e = 0;

        var.jinfo.code = cast(JohansenCode) gretl_bundle_get_int(b, "code", &e[0]);
        var.jinfo.rank = gretl_bundle_get_int(b, "rank", &e[1]);
        var.jinfo.seasonals = gretl_bundle_get_int(b, "seasonals", &e[2]);
        var.jinfo.R0 = gretl_bundle_get_matrix(b, "u", &e[3]);
        var.jinfo.R1 = gretl_bundle_get_matrix(b, "v", &e[4]);
        var.jinfo.S00 = gretl_bundle_get_matrix(b, "Suu", &e[5]);
        var.jinfo.S11 = gretl_bundle_get_matrix(b, "Svv", &e[6]);
        var.jinfo.S01 = gretl_bundle_get_matrix(b, "Suv", &e[7]);
        var.jinfo.Beta = gretl_bundle_get_matrix(b, "Beta", &e[8]);
        var.jinfo.Alpha = gretl_bundle_get_matrix(b, "Alpha", &e[9]);

        for (i=0; i<10; i++) {
            if (e[i]) {
                err = e[i];
                break;
            }
        }

        if (!err) {
            /* restriction matrices: may not be present */
            var.jinfo.R = gretl_bundle_get_matrix(b, "R", null);
            var.jinfo.q = gretl_bundle_get_matrix(b, "q", null);
            var.jinfo.Ra = gretl_bundle_get_matrix(b, "Ra", null);
            var.jinfo.qa = gretl_bundle_get_matrix(b, "qa", null);
        }
    }

    if (err && var.jinfo != null) {
        free(var.jinfo);
        var.jinfo = null;
    }

    return err;
}

/* Load constant plus seasonals, exogenous regressors, trend
   (if present, in each case) from var->X into @R.
*/

int load_exo_terms(GretlVar * var, GretlMatrixRaw * R)
{
    double * dest = R.val + var.T;
    const(double) * src;
    size_t sz = var.T * double.sizeof;
    int k = levels_order(var);
    int n = var.neqns;
    int t, j;

    for (t=0; t<var.T; t++) {
        gretl_matrix_set(R, t, 0, 1.0);
    }

    src = var.X.val + (var.ifc + n * (k-1)) * var.T;

    for (j=1; j<R.cols; j++) {
    memcpy(dest, src, sz);
    dest += var.T;
    src += var.T;
    }

    return 0;
}

double get_x_diff(const(GretlMatrixRaw) * X, const(GretlMatrixRaw) * P, int s, int j)
{
    double[2] x;
    int i;

    for (i=0; i<2; i++) {
    if (s >= 0) {
        x[i] = gretl_matrix_get(X, s, j);
    } else {
        x[i] = gretl_matrix_get(P, P.rows + s, j);
    }
    s--;
    }

    return x[0] - x[1];
}

/* See Johansen's 1995 book: the Sigma matrix is defined on page 188 and
   used in the formula for the C-matrix variance on page 190.

   We follow Paolo Paruolo, 'Asymptotic Inference on the Moving Average
   Impact Matrix in Cointegrated I(1) Var Systems', Econometric Theory
   13(1), Feb. 1997, pp. 79-118, by netting out the effect of any
   deterministic and/or otherwise exogenous terms: see the use of @R
   below.
*/

GretlMatrixRaw * johansen_Sigma(GretlVar * var, const(GretlMatrixRaw) * b, DatasetRaw * dset, GretlMatrixRaw * R, int * err)
{
    GretlMatrixRaw * S = null;
    GretlMatrixRaw * P = null;
    GretlMatrixRaw * X = null;
    GretlMatrixRaw * M = null;
    int p = var.neqns;
    int k = levels_order(var);
    int r = jrank(var);

    if (R != null) {
    /* fill R with seasonals, etc., as needed */
    load_exo_terms(var, R);
    }

    if (k > 1) {
    /* We'll need X values prior to the estimation range,
       to compute the required lagged X-difference(s).
    */
    P = gretl_matrix_data_subset(var.ylist, dset,
                     var.t1-k, var.t1-2,
                     MMissingCode.error,
                     err);
    }

    /* the lagged level of X, for X_{t-1} * \beta */
    X = gretl_matrix_data_subset(var.ylist, dset,
                 var.t1-1, var.t2-1,
                 MMissingCode.error,
                 err);

    if (!err && k < 2) {
    /* the simple case: no X-differences */
    M = gretl_matrix_alloc(X.rows, r);
    gretl_matrix_multiply(X, b, M);
    } else if (!*err) {
    /* lagged differences of X are needed */
    int mc = r + (k-1) * p;
    int lag, j, t, col;
    double dx;

    M = gretl_matrix_alloc(X.rows, mc);
    gretl_matrix_reuse(M, -1, r);
    gretl_matrix_multiply(X, b, M);
    gretl_matrix_reuse(M, -1, mc);

    for (lag=1, col=r; lag<k; lag++) {
        for (j=0; j<p; j++) {
        for (t=0; t<var.T; t++) {
            dx = get_x_diff(X, P, 1 + t - lag, j);
            gretl_matrix_set(M, t, col, dx);
        }
        col++;
        }
    }
    }

    if (!*err && R != null) {
    GretlMatrixRaw * E = gretl_matrix_alloc(var.T, M.cols);

    if (E == null) {
        *err = GretlError.alloc;
    } else {
        *err = gretl_matrix_multi_ols(M, R, null, E, null);
    }
    if (!*err) {
        S = gretl_covariance_matrix(E, 0, 0, err);
    }
    gretl_matrix_free(E);
    } else if (!*err) {
    S = gretl_covariance_matrix(M, 0, 0, err);
    }

    gretl_matrix_free(P);
    gretl_matrix_free(X);
    gretl_matrix_free(M);

    return S;
}

/* See Johansen's 1995 book, around page 190: produce the variance
   matrix of the vec of the "long-run matrix", C, and also, as a
   convenience, a matrix holding the standard errors of the elements
   of C.
*/

int add_johansen_lr_variance(GretlVar * var, DatasetRaw * dset)
{
    JohansenInfo *ji = var.jinfo;
    GretlMatrixBlock * B;
    GretlMatrixRaw * VC = null;
    GretlMatrixRaw * V1 = null;
    GretlMatrixRaw * V2 = null;
    GretlMatrixRaw * Cse = null;
    GretlMatrixRaw * aa = null;
    GretlMatrixRaw * abar = null;
    GretlMatrixRaw * b = null;
    GretlMatrixRaw * tmp = null;
    GretlMatrixRaw * xiT = null;
    GretlMatrixRaw * R = null;
    GretlMatrixRaw * Sigma = null;
    double x;
    int p = var.neqns;
    int k = levels_order(var);
    int nreg = 0;
    int i, j, nc;
    int err = 0;

    /* Will we need regressors in forming Sigma? Yes, if
       we have any unrestricted exogenous terms.
    */
    nreg = ji.seasonals;
    if (var.xlist != null) {
    nreg += var.xlist[0];
    }
    if (ji.code == JohansenCode.unrest_trend) {
    nreg++;
    }
    if (nreg > 0) {
        nreg++; /* for the constant */
    }
    if (nreg > 0) {
    /* Here we just allocate @R; this matrix will be filled
       in johansen_Sigma() if it's non-NULL.
    */
    R = gretl_matrix_alloc(var.T, nreg);
    if (R == null) {
        return GretlError.alloc;
    }
    }

    VC = gretl_matrix_alloc(p*p, p*p);
    if (VC == null) {
    gretl_matrix_free(R);
    return GretlError.alloc;
    }

    if (ji.Beta.rows > p) {
    b = gretl_matrix_alloc(p, ji.Beta.cols);
    for (j=0; j<ji.Beta.cols; j++) {
        for (i=0; i<p; i++) {
        x = gretl_matrix_get(ji.Beta, i, j);
        gretl_matrix_set(b, i, j, x);
        }
    }
    } else {
    b = gretl_matrix_copy(ji.Beta);
    }

    if (k > 1) {
    nc = ji.rank + (k-1) * p;
    } else {
    nc = ji.rank;
    }

    B = gretl_matrix_block_new(&V1, p, p,
                   &aa, ji.rank, ji.rank,
                   &abar, p, ji.rank,
                   &xiT, p, nc,
                   &tmp, p, p,
                   &V2, p, p, null);
    if (B == null) {
    err = GretlError.alloc;
    goto bailout;
    }

    /* V1 = C * vb.Omega * C' */
    err = gretl_matrix_qform(ji.JC, GretlMatrixMod.none, var.S,
                 V1, GretlMatrixMod.none);

    /* abar = alpha * inv(alpha'alpha) */
    gretl_matrix_multiply_mod(ji.Alpha, GretlMatrixMod.transpose,
                  ji.Alpha, GretlMatrixMod.none,
                  aa, GretlMatrixMod.none);
    err = gretl_invert_symmetric_matrix(aa);
    if (err) {
    goto bailout;
    } else {
    gretl_matrix_multiply(ji.Alpha, aa, abar);
    }

    /* xiT = (C'Gamma' - I(p)) * abar */
    gretl_matrix_multiply_mod(ji.JC, GretlMatrixMod.transpose,
                  ji.Gamma, GretlMatrixMod.transpose,
                  tmp, GretlMatrixMod.none);
    for (i=0; i<p; i++) {
    x = gretl_matrix_get(tmp, i, i);
    gretl_matrix_set(tmp, i, i, x - 1.0);
    }
    if (k > 1) {
    gretl_matrix_reuse(xiT, p, ji.rank);
    }
    gretl_matrix_multiply(tmp, abar, xiT);

    if (k > 1) {
    /* append k - 1 instances of C' */
    gretl_matrix_reuse(xiT, p, nc);
    for (i=1, j=ji.rank; i<k && !err; i++) {
        err = gretl_matrix_inscribe_matrix(xiT, ji.JC, 0, j,
                           GretlMatrixMod.transpose);
        j += p;
    }
    }

    Sigma = johansen_Sigma(var, b, dset, R, &err);

    if (!err) {
    /* V2 = xiT * invpd(Sigma) * xiT' */
    gretl_invpd(Sigma);
    gretl_matrix_qform(xiT, GretlMatrixMod.none, Sigma,
               V2, GretlMatrixMod.none);

    /* VC = (V2 ** V1) / T */
    gretl_matrix_kronecker_product(V2, V1, VC);
    gretl_matrix_divide_by_scalar(VC, var.T);
    Cse = gretl_matrix_alloc(p, p);
    for (i=0; i<p*p; i++) {
        Cse.val[i] = sqrt(gretl_matrix_get(VC, i, i));
    }
    }

 bailout:

    gretl_matrix_free(Sigma);
    gretl_matrix_free(R);
    gretl_matrix_free(b);
    gretl_matrix_block_destroy(B);

    if (err) {
    gretl_matrix_free(VC);
    } else {
    ji.JVC = VC;
    ji.Cse = Cse;
    }

    return err;
}

GretlMatrixRaw * johansen_perp(const(GretlMatrixRaw) * a, int * err)
{
    GretlMatrixRaw * tmp = gretl_matrix_alloc(a.cols, a.rows);
    GretlMatrixRaw * ret;

    gretl_matrix_transpose(tmp, a);
    ret = gretl_matrix_right_nullspace(tmp, err);
    gretl_matrix_free(tmp);

    return ret;
}

/* Add the "long-run impact" matrix, which Johansen (1995) calls C
   (Lütkepohl calls it $\Xi$).
*/

void add_johansen_C(JohansenInfo * j)
{
    GretlMatrixRaw * ap, bp;
    GretlMatrixRaw * tmp = null;
    GretlMatrixRaw * C1 = null;
    GretlMatrixRaw * JC = null;
    int p = j.Alpha.rows;
    int err = 0;

    /* make \alpha_{\perp} in ap */
    ap = johansen_perp(j.Alpha, &err);

    /* make \beta_{\perp} in bp */
    if (!err && j.Beta.rows > p) {
        double bik;
        int i, k;

        tmp = gretl_matrix_alloc(p, j.rank);
        for (k=0; k<j.rank; k++) {
            for (i=0; i<p; i++) {
                bik = gretl_matrix_get(j.Beta, i, k);
                gretl_matrix_set(tmp, i, k, bik);
            }
        }
        bp = johansen_perp(tmp, &err);
        gretl_matrix_free(tmp);
    } else if (!err) {
        bp = johansen_perp(j.Beta, &err);
    }

    if (err) {
        return;
    }

    /* JC = bp * inv(ap' * Gamma * bp) * ap' */
    tmp = gretl_matrix_alloc(ap.cols, j.Gamma.cols);
    C1 = gretl_matrix_alloc(ap.cols, bp.cols);
    JC = gretl_matrix_alloc(p, p);
    gretl_matrix_multiply_mod(ap, GretlMatrixMod.transpose,
                              j.Gamma, GretlMatrixMod.none,
                              tmp, GretlMatrixMod.none);
    gretl_matrix_multiply(tmp, bp, C1);
    err = gretl_invert_general_matrix(C1);
    if (!err) {
        gretl_matrix_realloc(tmp, p, C1.cols);
        gretl_matrix_multiply(bp, C1, tmp);
        gretl_matrix_multiply_mod(tmp, GretlMatrixMod.none,
                                  ap, GretlMatrixMod.transpose,
                                  JC, GretlMatrixMod.none);
    }

    if (err) {
        gretl_matrix_free(JC);
    } else {
        j.JC = JC;
    }

    gretl_matrix_free(tmp);
    gretl_matrix_free(C1);
    gretl_matrix_free(ap);
    gretl_matrix_free(bp);
}

GretlBundle * johansen_bundlize(GretlVar * var, DatasetRaw * dset)
{
    GretlBundle * b = gretl_bundle_new();
    JohansenInfo *j = var.jinfo;

    if (b == null) {
        return null;
    }

    /* "j->ID"? */
    gretl_bundle_set_int(b, "code", j.code);
    gretl_bundle_set_int(b, "rank", j.rank);
    gretl_bundle_set_int(b, "seasonals", j.seasonals);

    if (j.lrdf > 0 && !na(j.ll0)) {
        gretl_bundle_set_scalar(b, "ll0", j.ll0);
        gretl_bundle_set_scalar(b, "bdf", j.lrdf);
    }

    gretl_bundle_set_matrix(b, "u", j.R0);
    gretl_bundle_set_matrix(b, "v", j.R1);
    gretl_bundle_set_matrix(b, "Suu", j.S00);
    gretl_bundle_set_matrix(b, "Svv", j.S11);
    gretl_bundle_set_matrix(b, "Suv", j.S01);
    gretl_bundle_set_matrix(b, "evals", j.evals);
    gretl_bundle_set_matrix(b, "Beta", j.Beta);
    gretl_bundle_set_matrix(b, "Alpha", j.Alpha);
    gretl_bundle_set_matrix(b, "Bvar", j.Bvar);
    gretl_bundle_set_matrix(b, "Bse", j.Bse);

    if (j.Gamma != null) {
        gretl_bundle_set_matrix(b, "Gamma", j.Gamma);
    }
    if (j.R != null) {
        gretl_bundle_set_matrix(b, "R", j.R);
    }
    if (j.q != null) {
        gretl_bundle_set_matrix(b, "q", j.q);
    }
    if (j.Ra != null) {
        gretl_bundle_set_matrix(b, "Ra", j.Ra);
    }
    if (j.qa != null) {
        gretl_bundle_set_matrix(b, "qa", j.qa);
    }
    if (j.G != null) {
    gretl_bundle_donate_data(b, "G_array", j.G,
                 GretlType.array);
    j.G = null; /* don't double-free */
    }

    if (j.Alpha != null && j.Beta != null && j.Gamma != null) {
        add_johansen_C(j);
        if (j.JC != null) {
            gretl_bundle_set_matrix(b, "long_run", j.JC);
        add_johansen_lr_variance(var, dset);
        if (j.JVC != null) {
        gretl_bundle_set_matrix(b, "cov_long_run", j.JVC);
        }
        if (j.Cse != null) {
        gretl_bundle_set_matrix(b, "se_long_run", j.Cse);
        }
    }
    }

    return b;
}

int gretl_VAR_serialize(GretlVar * var, SavedObjectFlags flags, Prn * prn)
{
    char * xmlname = null;
    int g = var.neqns;
    int m = g * g + g;
    int i, err = 0;

    if (var.name == null || *var.name == '\0') {
        xmlname = gretl_strdup("none");
    } else {
        xmlname = gretl_xml_encode(var.name);
    }

    pprintf(prn, "<gretl-VAR name=\"%s\" saveflags=\"%d\" ".ptr, xmlname, cast(int) flags);
    free(xmlname);

    pprintf(prn, "ecm=\"%d\" neqns=\"%d\" order=\"%d\" detflags=\"%d\" ".ptr,
            (var.ci == GretlCmdIndex.vecm), var.neqns, var.order, var.detflags);

    if (var.robust) {
        gretl_xml_put_int("robust", var.robust, prn);
    }

    if (var.LBs > 0 && !na(var.LB)) {
        /* Portmanteau test */
        gretl_xml_put_double("LB", var.LB, prn);
        gretl_xml_put_int("LBs", var.LBs, prn);
    }

    pputs(prn, ">\n");

    gretl_xml_put_tagged_list("lags", var.lags, prn);
    gretl_xml_put_tagged_list("ylist", var.ylist, prn);
    gretl_xml_put_tagged_list("xlist", var.xlist, prn);
    gretl_xml_put_tagged_list("rlist", var.rlist, prn);

    gretl_push_c_numeric_locale();

    if (var.Fvals != null) {
        gretl_xml_put_double_array("Fvals", var.Fvals, m, prn);
    }

    if (var.Ivals != null) {
        gretl_xml_put_double_array("Ivals", var.Ivals, N_IVALS, prn);
    }

    if (var.X != null && var.Y != null) {
        /* could be fiddly to reconstruct, needed for IRF bootstrap */
        gretl_matrix_serialize(var.X, "X", prn);
        gretl_matrix_serialize(var.Y, "Y", prn);
    }

    if (var.ord != null) {
        gretl_matrix_serialize(var.ord, "ord", prn);
    }

    if (var.ci == GretlCmdIndex.vecm) {
        /* this is hard to reconstruct for VECMs */
        gretl_matrix_serialize(var.A, "A", prn);
    }

    gretl_pop_c_numeric_locale();

    pputs(prn, "<equations>\n");

    for (i=0; i<var.neqns; i++) {
        gretl_model_serialize(var.models[i], cast(SavedObjectFlags) 0, prn);
    }

    pputs(prn, "</equations>\n");

    if (var.jinfo != null) {
        johansen_serialize(var.jinfo, prn);
    }

    pputs(prn, "</gretl-VAR>\n");

    return err;
}

/* Clean-up of temporary "jinfo" used in IRF bootstrapping. We can't use
   the regular destructor, johansen_info_free(), because the primary
   matrix pointers on var->jinfo are in this case borrowed from a
   $system bundle; but we do need to free any "extra" matrices that have
   been added.
*/

void free_temp_jinfo(GretlVar * var)
{
    /* handle extra matrices */
    gretl_matrix_free(var.jinfo.YY);
    gretl_matrix_free(var.jinfo.RR);
    gretl_matrix_free(var.jinfo.BB);

    free(var.jinfo);
}

/* When destroying a VAR that was reconstituted by reading from a
   $system bundle, we must be careful not to free any elements which
   were just borrowed from the bundle; we set the associated pointers to
   NULL.
*/

void destroy_VAR_from_bundle(GretlVar * var)
{
    GretlMatrixRaw **[9] mm = [&var.A, &var.C, &var.X, &var.Y, &var.XTX, &var.B, &var.S, &var.E, null];
    int i, imin = 0;

    /* var->A is NOT borrowed */
    imin = 1;

    /* nullify borrowed matrices */
    for (i=imin; mm[i] != null; i++) {
    *mm[i] = null;
    }

    /* nullify borrowed lists */
    var.xlist = var.ylist = var.rlist = null;

    if (var.jinfo != null) {
    free_temp_jinfo(var);
    var.jinfo = null;
    }

    gretl_VAR_free(var);
}

/* Reconstitute a VAR or VECM (partially) from the bundle @b, for the
   purposes of generating an FEVD or IRF. How far we need to go in
   rebuilding the GRETL_VAR struct depends on the purpose, and in the
   IRF case, whether bootstrapping is required. For FEVD we basically
   just need the A and C matrices from the bundle.
*/

GretlVar * build_VAR_from_bundle(GretlBundle * b, int irf, int boot, int * err)
{
    GretlVar * var = cast(GretlVar *) malloc(GretlVar.sizeof);
    int i;
    int[7] ierr;

    if (var == null) {
        *err = GretlError.alloc;
        return null;
    }

    gretl_VAR_clear(var);

    if (gretl_bundle_get_int(b, "ecm", null)) {
        var.ci = GretlCmdIndex.vecm;
    } else {
        var.ci = GretlCmdIndex.var;
    }

    /* get integer dimensions */
    var.neqns = gretl_bundle_get_int(b, "neqns", &ierr[0]);
    var.order = gretl_bundle_get_int(b, "order", &ierr[1]);
    var.ncoeff = gretl_bundle_get_int(b, "ncoeff", &ierr[2]);
    var.t1 = gretl_bundle_get_int(b, "t1", &ierr[3]);
    var.t2 = gretl_bundle_get_int(b, "t2", &ierr[4]);
    var.df = gretl_bundle_get_int(b, "df", &ierr[6]);
    var.T  = gretl_bundle_get_int(b, "T", &ierr[5]);
    for (i=0; i<7; i++) {
        if (ierr[i]) {
            *err = ierr[i];
            break;
        }
    }

    if (!*err) {
    /* convert sample range to zero-based */
    var.t1 -= 1;
    var.t2 -= 1;
    var.ifc = gretl_bundle_get_int(b, "ifc", null);
        /* note: borrowing of lists from bundle */
        var.ylist = gretl_bundle_get_list(b, "ylist", null);
        var.xlist = gretl_bundle_get_list(b, "xlist", null);
        var.rlist = gretl_bundle_get_list(b, "rlist", null);
    }

    if (!*err) {
        /* note: borrowing of matrices */
        GretlMatrixRaw **[8] mm = [&var.A, &var.C, &var.X, &var.Y, &var.XTX, &var.B, &var.S, &var.E];
        static const(char) *[8] keys = ["A".ptr, "C".ptr, "X".ptr, "Y".ptr, "xtxinv".ptr, "coeff".ptr, "sigma".ptr, "uhat".ptr];
        int n = 2 + 2*irf + 4*boot;

        for (i=0; i<n && !*err; i++) {
            *mm[i] = gretl_bundle_get_matrix(b, keys[i], err);
        }
    if (!*err) {
        /* matrix @A from the bundle must be converted into
           the full companion matrix
        */
        GretlMatrixRaw * CA;

        if (var.ci == GretlCmdIndex.vecm) {
        /* Backward compatibility for VECM order: before
           March 2022 this was recorded as the order in
           differences, and we need to reconstitute that
           value here -- unless and until we revise the
           internal definition of var->order for VECMs.
        */
        var.order = var.A.cols / var.neqns - 1;
        }
        CA = companionize(var.A, err);
        var.A = CA;
    }
        if (!*err && var.ci == GretlCmdIndex.vecm && boot) {
            /* not fully ready, won't be reached at present? */
            GretlBundle * jb;

            jb = gretl_bundle_get_bundle(b, "vecm_info", err);
            if (!*err) {
                *err = retrieve_johansen_basics(var, jb);
            }
        }
    }

    if (!*err && var.neqns == 0) {
    /* maybe we have a minimal "home-made" bundle? */
    if (var.C != null && var.C.cols > 0) {
        var.neqns = var.C.cols;
    } else {
        *err = GretlError.invarg;
    }
    }

    if (*err) {
    destroy_VAR_from_bundle(var);
        var = null;
    }

    return var;
}

GretlMatrixRaw * gretl_FEVD_from_bundle(GretlBundle * b, int targ, int shock, const(DatasetRaw) * dset, int * err)
{
    GretlVar * var = build_VAR_from_bundle(b, 0, 0, err);
    GretlMatrixRaw * ret = null;

    if (var != null) {
        ret = gretl_VAR_get_FEVD_matrix(var, targ, shock, 0, dset, err);
    destroy_VAR_from_bundle(var);
    }

    return ret;
}

GretlMatrixRaw * gretl_IRF_from_bundle(GretlBundle * b, int targ, int shock, double alpha, const(DatasetRaw) * dset, int * err)
{
    GretlMatrixRaw * ret = null;
    GretlVar * var = null;
    int boot = 0;

    if (alpha != 0 && (alpha < 0.01 || alpha > 0.6)) {
        *err = GretlError.invarg;
        return null;
    } else if (alpha != 0) {
        boot = 1;
    }

    var = build_VAR_from_bundle(b, 1, boot, err);

    if (!*err && var != null) {
        ret = gretl_VAR_get_impulse_response(var, targ, shock,
                                             0, alpha, dset, err);
    }

    if (var != null) {
    destroy_VAR_from_bundle(var);
    }

    return ret;
}

GretlMatrixRaw * make_detflags_matrix(GretlVar * var)
{
    GretlMatrixRaw * d = gretl_zero_matrix_new(1, 3);

    if (var.ifc || (var.detflags & DetFlags.const_)) {
        d.val[0] = 1;
    }
    if (var.detflags & DetFlags.trend) {
        d.val[1] = 1;
    }
    if (var.detflags & DetFlags.seas) {
        d.val[2] = 1;
    }

    return d;
}

int gretl_VAR_bundlize(GretlVar * var, DatasetRaw * dset, GretlBundle * b)
{
    GretlMatrixRaw * d;
    int err = 0;

    gretl_bundle_set_string(b, "command", gretl_command_word(var.ci));

    if (var.name != null) {
        gretl_bundle_set_string(b, "name", var.name);
    }
    gretl_bundle_set_int(b, "ecm", var.ci == GretlCmdIndex.vecm);
    gretl_bundle_set_int(b, "neqns", var.neqns);
    gretl_bundle_set_int(b, "ncoeff", var.ncoeff);
    gretl_bundle_set_int(b, "robust", var.robust);
    gretl_bundle_set_int(b, "t1", var.t1 + 1);
    gretl_bundle_set_int(b, "t2", var.t2 + 1);
    gretl_bundle_set_int(b, "df", var.df);
    gretl_bundle_set_int(b, "T", var.T);
    gretl_bundle_set_int(b, "ifc", var.ifc);
    /* changed 2022-03-31 */
    gretl_bundle_set_int(b, "order", levels_order(var));

    if (var.models != null && var.models[0] != null) {
    gretl_bundle_set_int(b, "sample_t1", var.models[0].smpl.t1 + 1);
    gretl_bundle_set_int(b, "sample_t2", var.models[0].smpl.t2 + 1);
    }

    gretl_bundle_set_scalar(b, "lnl", var.ll);
    gretl_bundle_set_scalar(b, "ldet", var.ldet);

    if (var.LBs > 0 && !na(var.LB)) {
        /* Portmanteau test */
        gretl_bundle_set_scalar(b, "Ljung_Box", var.LB);
        gretl_bundle_set_scalar(b, "LB_order", var.LBs);
    }

    /* deterministic terms: pack in matrix */
    d = make_detflags_matrix(var);
    gretl_bundle_donate_data(b, "detflags", d,
                             GretlType.matrix);

    /* lists: lags, ylist, xlist, rlist */
    if (var.lags != null) {
        GretlMatrixRaw * v = gretl_list_to_vector(var.lags, &err);

        if (!err) {
            gretl_bundle_donate_data(b, "lags", v,
                                     GretlType.matrix);
            gretl_bundle_set_note(b, "lags", "gappy lags vector");
        }
    }
    if (var.ylist != null) {
        gretl_bundle_set_list(b, "ylist", var.ylist);
    }
    if (var.xlist != null) {
        gretl_bundle_set_list(b, "xlist", var.xlist);
    }
    if (var.rlist != null) {
        gretl_bundle_set_list(b, "rlist", var.rlist);
    }

    /* doubles arrays: Fvals, Ivals? */

    if (var.A != null) {
    /* companion matrix */
    if (var.A.rows == var.A.cols) {
        gretl_bundle_set_matrix(b, "A", var.A);
    } else {
        /* A needs trimming */
        GretlMatrixRaw * A = decompanionize(var.A, var.neqns, GretlMatrixMod.none);

        if (A != null) {
        gretl_bundle_donate_data(b, "A", A, GretlType.matrix);
        }
    }
    }
    if (var.C != null) {
        gretl_bundle_set_matrix(b, "C", var.C);
    }
    if (var.B != null) {
        gretl_bundle_set_matrix(b, "coeff", var.B);
    }
    if (var.S != null) {
        gretl_bundle_set_matrix(b, "sigma", var.S);
    }
    if (var.XTX != null) {
        gretl_bundle_set_matrix(b, "xtxinv", var.XTX);
    }
    if (var.E != null) {
        gretl_bundle_set_matrix(b, "uhat", var.E);
    }
    if (var.X != null && var.Y != null) {
        gretl_bundle_set_matrix(b, "X", var.X);
        gretl_bundle_set_matrix(b, "Y", var.Y);
    }
    if (var.ord != null) {
        gretl_bundle_set_matrix(b, "ord", var.ord);
    }

    if (var.jinfo != null) {
        /* VECM specific info */
        GretlBundle * jb = johansen_bundlize(var, dset);

        if (jb != null) {
            err = gretl_bundle_donate_data(b, "vecm_info", jb,
                                           GretlType.bundle);
        }
    }

    return err;
}

