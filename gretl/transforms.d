/* Port of transforms.c */

module gretl.transforms;


import core.stdc.ctype : isdigit;
import core.stdc.math : ceil, fabs, log;
import core.stdc.stdio : fprintf, sprintf, stderr;
import core.stdc.stdlib : abs, free, malloc, qsort, realloc;
import core.stdc.string : strcat, strcmp, strcpy, strlen, strncat, strncmp;

import gretl.dataset : dataset_add_series, panel_sample_size, series_get_lag, series_get_midas_freq,
    series_get_midas_period, series_get_label, series_get_parent_id, series_get_parent_name,
    series_get_stack_level, series_get_string_for_value, series_get_transform,
    series_is_coded, series_is_integer_valued, series_record_label, series_set_discrete,
    series_set_label, series_set_lag, series_set_midas_freq, series_set_midas_period,
    series_set_parent, series_set_transform, series_zero_flags, is_string_valued;
import gretl.dataseth : dataset_is_panel, dataset_is_seasonal, dataset_is_time_series;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_errors : gretl_errmsg_set, gretl_errmsg_sprintf, gretl_warnmsg_set, set_gretl_warning;
import gretl.gretl_errorsh : GretlError, GretlWarningCodes;
import gretl.gretl_list : gretl_is_midas_list, gretl_list_append_list, gretl_list_append_term,
    gretl_list_delete_at_pos, gretl_list_new, gretl_list_set_midas, gretl_null_list,
    in_gretl_list;
import gretl.gretl_matrix : gretl_matrix_QR_decomp, gretl_matrix_alloc, gretl_matrix_free,
    gretl_zero_matrix_new;
import gretl.random : gretl_rand_int_minmax;
import gretl.gretl_matrixh;
import gretl.gstrfuncs : g_strdup_printf;
import gretl.gmem : g_free;
import gretl.libgretlh;
import gretl.matrix_extrah : MMissingCode;
import gretl.missingh : NADBL, na;
import gretl.strutils : gretl_strdup, gretl_utf8_truncate, make_varname_unique, strings_array_add,
    strings_array_free;
import gretl.gretl_utils : accept_as_discrete, count_distinct_values, gretl_compare_doubles, gretl_isdummy, gretl_iszero, rearrange_id_array;
import gretl.genmain : current_series_index, series_index;
import gretl.matrix_extra : gretl_matrix_data_subset;
import gretl.uservar : gretl_is_user_var;
import gretl.genfuncs : orthdev_series, panel_unit_first_obs, standardize_series;

/* ---- declarations of not-yet-ported items used below ----
   Replace each with an import when its defining file is ported. */

extern (C) {
    int gretl_function_depth();
}

/**
 * SECTION:transforms
 * @short_description: standard transformations of series in the dataset
 * @title: Transformations
 * @include: libgretl.h
 *
 * Functions to generate standard transformations (logs, lags, first
 * differences and so on) of series in the dataset.
 */

private enum VarComp {
    identical,
    new_has_missing,
    old_has_missing,
    differ
}

/* "transformation codes" that extend the command index values */

private enum ExtraTransform {
    inverse = GretlCmdIndex.nc + 1,
    resample,
    hfdiff,
    hfldiff
}

/* Note: in this comparison @test is the newly generated series
   and @cand is the previously existing candidate series.
*/

private int check_vals(const(double) * test, const(double) * cand, int n)
{
    int t, ret = VarComp.identical;

    for (t=0; t<n; t++) {
        if (test[t] == cand[t] || (na(test[t]) && na(cand[t]))) {
            /* values the same, or both NA */
            continue;
        }
        /* values differ at obs @t */
        if (na(test[t]) && !na(cand[t])) {
            if (ret == VarComp.identical || ret == VarComp.new_has_missing) {
                ret = VarComp.new_has_missing;
            } else {
                ret = VarComp.differ;
            }
        } else if (!na(test[t]) && na(cand[t])) {
            if (ret == VarComp.identical || ret == VarComp.old_has_missing) {
                ret = VarComp.old_has_missing;
            } else {
                ret = VarComp.differ;
            }
        } else {
            ret = VarComp.differ;
        }
        if (ret == VarComp.differ) {
            break;
        }
    }

    return ret;
}

private int make_transform_varname(char * vname, const(char) * orig, int ci,
                                   int aux, int len)
{
    *vname = '\0';

    if (ci == GretlCmdIndex.diff) {
        strcpy(vname, "d_");
        strncat(vname, orig, len - 2);
    } else if (ci == GretlCmdIndex.ldiff) {
        strcpy(vname, "ld_");
        strncat(vname, orig, len - 3);
    } else if (ci == GretlCmdIndex.sdiff) {
        strcpy(vname, "sd_");
        strncat(vname, orig, len - 3);
    } else if (ci == GretlCmdIndex.orthdev) {
        strcpy(vname, "o_");
        strncat(vname, orig, len - 2);
    } else if (ci == GretlCmdIndex.logs) {
        strcpy(vname, "l_");
        strncat(vname, orig, len - 2);
    } else if (ci == GretlCmdIndex.square) {
        strcpy(vname, "sq_");
        strncat(vname, orig, len - 3);
    } else if (ci == GretlCmdIndex.lags) {
        char[12] ext;

        if (aux >= 0) {
            /* an actual lag (or contemporaneous value) */
            sprintf(ext.ptr, "_%d", aux);
            strncat(vname, orig, len - strlen(ext.ptr));
            strcat(vname, ext.ptr);
        } else {
            /* a lead: in naming the output series we'll
               try to avoid appending a plain digit to
               a series name that already ends in a digit,
               since this may get confusing
            */
            int n1, n2;

            sprintf(ext.ptr, "%d", -aux);
            n1 = cast(int) strlen(orig);
            n2 = cast(int) strlen(ext.ptr);
            if (n1 + n2 + 2 < VNAMELEN) {
                if (isdigit(orig[n1-1])) {
                    /* separate the digits */
                    sprintf(vname, "%s_f%d", orig, -aux);
                } else {
                    sprintf(vname, "%s%d", orig, -aux);
                }
            } else {
                /* Oh, well, can't be helped */
                strncat(vname, orig, len - n2);
                strcat(vname, ext.ptr);
            }
        }
    } else if (ci == GretlCmdIndex.dummify) {
        char[16] ext;

        if (aux < 0) {
            sprintf(ext.ptr, "_m%d", abs(aux));
        } else {
            sprintf(ext.ptr, "_%d", aux);
        }
        strcpy(vname, "D");
        strncat(vname, orig, len - strlen(ext.ptr) - 1);
        strcat(vname, ext.ptr);
    } else if (ci == ExtraTransform.inverse) {
        strcpy(vname, "i_");
        strncat(vname, orig, len - 2);
    } else if (ci == GretlCmdIndex.stdize) {
        if (aux < 0) {
            strcpy(vname, "c_");
        } else {
            strcpy(vname, "s_");
        }
        strncat(vname, orig, len - 2);
    } else if (ci == ExtraTransform.resample) {
        strcpy(vname, "rs_");
        strncat(vname, orig, len - 3);
    }

    return 0;
}

private int make_transform_label(char * label, const(char) * parent,
                                 int ci, int aux)
{
    char * s = null;
    int err = 0;

    if (ci == GretlCmdIndex.diff) {
        s = g_strdup_printf(_("= first difference of %s").ptr, parent);
    } else if (ci == GretlCmdIndex.ldiff) {
        s = g_strdup_printf(_("= log difference of %s").ptr, parent);
    } else if (ci == GretlCmdIndex.sdiff) {
        s = g_strdup_printf(_("= seasonal difference of %s").ptr, parent);
    } else if (ci == GretlCmdIndex.logs) {
        s = g_strdup_printf(_("= log of %s").ptr, parent);
    } else if (ci == GretlCmdIndex.square) {
        s = g_strdup_printf(_("= %s squared").ptr, parent);
    } else if (ci == GretlCmdIndex.lags) {
        /* @aux = lag */
        if (aux >= 0) {
            s = g_strdup_printf("= %s(t - %d)", parent, aux);
        } else {
            s = g_strdup_printf("= %s(t + %d)", parent, -aux);
        }
    } else if (ci == ExtraTransform.inverse) {
        s = g_strdup_printf("= 1/%s", parent);
    } else if (ci == GretlCmdIndex.stdize) {
        /* @aux = dfcorr */
        if (aux < 0) {
            s = g_strdup_printf(_("= centered %s").ptr, parent);
        } else {
            s = g_strdup_printf(_("= standardized %s").ptr, parent);
        }
    } else if (ci == ExtraTransform.resample) {
        s = g_strdup_printf(_("= resampled %s").ptr, parent);
    } else if (ci == ExtraTransform.hfdiff) {
        s = g_strdup_printf(_("= high-frequency difference of %s").ptr, parent);
    } else if (ci == ExtraTransform.hfldiff) {
        s = g_strdup_printf(_("= high-frequency log difference of %s").ptr, parent);
    } else {
        err = 1;
    }

    if (!err) {
        strcpy(label, gretl_utf8_truncate(s, MAXLABEL-1));
        g_free(s);
    }

    return err;
}

/**
 * standard_lag_of:
 * @v: ID number of series to test.
 * @parent: ID of potential parent series.
 * @dset: dataset information.
 *
 * Returns: the lag order of series @v, if it is marked as
 * a lag of @parent, otherwise 0.
 */

int standard_lag_of(int v, int parent, const(DatasetRaw) * dset)
{
    int pv = 0, ret = 0;

    if (dset == null || v <= 0 || v >= dset.v) {
        return 0;
    }

    if (series_get_transform(dset, v) == GretlCmdIndex.lags) {
        pv = series_get_parent_id(dset, v);
        if (pv == parent) {
            ret = series_get_lag(dset, v);
        }
    }

    return ret;
}

/**
 * is_standard_diff:
 * @v: ID number of variable to test.
 * @dset: dataset information.
 * @parent: location to receive ID number of parent variable,
 * or NULL.
 *
 * Returns: 1 if the variable @v is marked as being the first
 * difference of some "parent" variable in the dataset,
 * otherwise 0.
 */

int is_standard_diff(int v, const(DatasetRaw) * dset, int * parent)
{
    int pv = 0, ret = 0;

    if (v <= 0 || v >= dset.v) {
        return 0;
    }

    if (series_get_transform(dset, v) == GretlCmdIndex.diff) {
        pv = series_get_parent_id(dset, v);
        if (pv > 0) {
            if (parent != null) {
                *parent = pv;
            }
            ret = 1;
        }
    }

    return ret;
}

private void make_xp_varname(char * vname, int v1, int v2,
                             const(DatasetRaw) * dset,
                             int len)
{
    int n1 = cast(int) strlen(dset.varname[v1]);
    int n2 = cast(int) strlen(dset.varname[v2]);
    int cut = n1 + n2 + 1 - len;
    int cut1 = 0, cut2 = 0;

    if (cut > 0) {
        cut1 = cut2 = cut/2;
        if (cut % 2) {
            cut2++;
        }
    }

    *vname = '\0';
    strncat(vname, dset.varname[v1], n1 - cut1);
    strcat(vname, "_");
    strncat(vname, dset.varname[v2], n2 - cut2);
}

/* Array into which to write a generated variable, prior to testing
   whether or not the same var already exists.
*/

private __gshared double * testvec_x;
private __gshared int testvec_nx;

private double * testvec(int n)
{
    int t;

    if (n == 0) {
        /* clean up */
        free(testvec_x);
        testvec_x = null;
        testvec_nx = 0;
        return null;
    }

    if (n > testvec_nx) {
        double * y = cast(double *) realloc(testvec_x, n * (*testvec_x).sizeof);

        if (y == null) {
            free(testvec_x);
            testvec_x = null;
            testvec_nx = 0;
            return null;
        }

        testvec_nx = n;
        testvec_x = y;
    }

    for (t=0; t<n; t++) {
        testvec_x[t] = NADBL;
    }

    return testvec_x;
}

/**
 * gretl_transforms_cleanup:
 *
 * Called by libgretl_cleanup(). Frees any memory allocated
 * as workspace for the creation of transformed variables.
 */

void gretl_transforms_cleanup()
{
    testvec(0);
}

/* write lagged values of variable v into xlag */

private int get_lag(int v, int lag, double * xlag,
                    const(DatasetRaw) * dset)
{
    const(double) * x = dset.Z[v];
    int t1 = (lag > 0)? lag : 0;
    int t2 = dset.n - 1;
    int t, s;

    for (t=0; t<dset.n; t++) {
        xlag[t] = NADBL;
    }

    for (t=t1; t<=t2; t++) {
        s = t - lag;
        if (dataset_is_panel(dset) &&
            t / dset.pd != s / dset.pd) {
            continue;
        }
        if (s >= 0 && s < dset.n) {
            xlag[t] = x[s];
        }
    }

    return 0;
}

/* write log of variable v into logvec */

private int get_log(int v, double * logvec, const(DatasetRaw) * dset)
{
    double xx;
    int t, err = 0;

    for (t=dset.t1; t<=dset.t2 && !err; t++) {
        xx = dset.Z[v][t];
        if (na(xx) || xx <= 0.0) {
            logvec[t] = NADBL;
            set_gretl_warning(GretlWarningCodes.genmiss);
        } else {
            logvec[t] = log(xx);
        }
    }

    return err;
}

/* write some sort of difference of variable v into diffvec */

private int get_diff(int v, double * diffvec, int ci,
                     const(DatasetRaw) * dset)
{
    double x0, x1;
    int t, t0, t1;

    t0 = (ci == GretlCmdIndex.sdiff)? dset.pd : 1;
    t1 = (dset.t1 > t0)? dset.t1 : t0;

    for (t=t1; t<=dset.t2; t++) {
        if (dset.structure == TsCodes.stacked_time_series &&
            panel_unit_first_obs(t, dset)) {
            continue;
        }

        x0 = dset.Z[v][t];
        x1 = dset.Z[v][t - t0];

        if (ci == GretlCmdIndex.ldiff) {
            if (!na(x0) && !na(x1) && x0 > 0 && x1 > 0) {
                diffvec[t] = log(x0) - log(x1);
            }
        } else {
            if (!na(x0) && !na(x1)) {
                diffvec[t] = x0 - x1;
            }
        }
    }

    return 0;
}

/* write a matrix containing high-frequency differences of
   the series in @list over the current sample range
*/

private GretlMatrixRaw * get_hf_diffs(const(int) * list,
                                   int ci, double mult,
                                   const(DatasetRaw) * dset,
                                   int * err)
{
    GretlMatrixRaw * dX = null;
    const(double) ** Z = cast(const(double) **) dset.Z;
    double x0, x1, dti;
    int n = list[0];
    int v1 = list[1];
    int i, vi, s, t, T;

    T = dset.t2 - dset.t1 + 1;
    dX = gretl_zero_matrix_new(T, n);

    if (dX == null) {
        *err = GretlError.alloc;
        return null;
    }

    for (t=T-1; t>=0; t--) {
        s = t + dset.t1;
        for (i=0; i<n; i++) {
            vi = list[i+1];
            x0 = Z[vi][s];
            if (i < n-1) {
                x1 = Z[list[i+2]][s];
            } else if (s > 0) {
                x1 = Z[v1][s-1];
            } else {
                x1 = NADBL;
            }
            dti = NADBL;
            if (!na(x0) && !na(x1)) {
                if (ci == GretlCmdIndex.diff) {
                    dti = mult * (x0 - x1);
                } else if (x0 > 0 && x1 > 0) {
                    dti = mult * log(x0/x1);
                }
            }
            gretl_matrix_set(dX, t, i, dti);
        }
    }

    return dX;
}

/* orthogonal deviations */

private int get_orthdev(int v, double * xvec, const(DatasetRaw) * dset)
{
    return orthdev_series(dset.Z[v], xvec, dset);
}

/* write square or cross-product into xvec */

private int get_xpx(int vi, int vj, double * xvec, const(DatasetRaw) * dset)
{
    double xit, xjt;
    int t;

    for (t=dset.t1; t<=dset.t2; t++) {
        xit = dset.Z[vi][t];
        xjt = dset.Z[vj][t];
        if (na(xit) || na(xjt)) {
            xvec[t] = NADBL;
        } else {
            xvec[t] = xit * xjt;
        }
    }

    return 0;
}

/* write reciprocal into xvec */

private int get_inverse(int v, double * xvec, const(DatasetRaw) * dset)
{
    double xt;
    int t;

    for (t=dset.t1; t<=dset.t2; t++) {
        xt = dset.Z[v][t];
        if (na(xt) || xt == 0.0) {
            xvec[t] = NADBL;
        } else {
            xvec[t] = 1.0 / xt;
        }
    }

    return 0;
}

/* write standardized series v into targ */

private int get_std(int v, int dfc, double * targ, const(DatasetRaw) * dset)
{
    const(double) * x = dset.Z[v];

    return standardize_series(x, targ, dfc, dset);
}

/* write resampled series into rsvec */

private int get_resampled(int v, double * rsvec,
                          const(DatasetRaw) * dset,
                          const(int) * z)
{
    const(double) * x = dset.Z[v];
    int t, j, s;

    if (dataset_is_panel(dset)) {
        int n = panel_sample_size(dset);
        int i, T = dset.pd;

        s = dset.t1;
        for (i=0; i<n; i++) {
            j = z[i] * T;
            for (t=0; t<T; t++) {
                rsvec[s++] = x[j + t];
            }
        }
    } else {
        s = 0;
        for (t=dset.t1; t<=dset.t2; t++) {
            j = z[s++];
            rsvec[t] = x[j];
        }
    }

    return 0;
}

/* write dummy for (v == value) into xvec */

private int get_discdum(int v, double val, double * xvec,
                        const(DatasetRaw) * dset)
{
    double xt;
    int t;

    for (t=dset.t1; t<=dset.t2; t++) {
        xt = dset.Z[v][t];
        if (na(xt)) {
            xvec[t] = NADBL;
        } else {
            xvec[t] = (xt == val)? 1.0 : 0.0;
        }
    }

    return 0;
}

private enum TransformResult {
    var_added_ok,
    var_exists_ok,
    var_add_failed,
    varname_duplicate
}

/* Note: new behavior as of 2008-02-11 */

private enum TR_OVERWRITE = true;

private int transform_handle_duplicate(int ci, int lag, int v,
                                       const(double) * x,
                                       const(char) * label,
                                       DatasetRaw * dset,
                                       int parent,
                                       int origv)
{
    const(char) * vlabel = series_get_label(dset, v);
    int ret = TransformResult.varname_duplicate;
    int no_overwrite = 0;
    int t, ok = 0;

    if (vlabel != null && *label && *vlabel) {
        if (!strcmp(label, vlabel)) {
            /* labels identical, so OK? */
            ok = 1;
        } else {
            /* labels differ */
            no_overwrite = 1;
        }
    }

    if (!no_overwrite) {
        int vp = series_get_parent_id(dset, v);

        no_overwrite = parent != vp;
    }

    if (TR_OVERWRITE) {
        if (!ok && v < origv && !no_overwrite) {
            ok = 1;
        }
    }

    if (ok) {
        for (t=0; t<dset.n; t++) {
            dset.Z[v][t] = x[t];
        }
        if (*label != '\0') {
            series_set_label(dset, v, label);
        }
        series_set_transform(dset, v, ci);
        series_set_lag(dset, v, lag);
        series_zero_flags(dset, v);
        ret = TransformResult.var_exists_ok;
    }

    return ret;
}

private int check_add_transform(int ci, int lag, int vnum, const(double) * x,
                                const(char) * vname, const(char) * label,
                                DatasetRaw * dset, int parent, int origv)
{
    int t, ret = TransformResult.var_added_ok;

    if (vnum < dset.v) {
        /* a variable of this name already exists */
        int chk = check_vals(x, dset.Z[vnum], dset.n);

        /* heuristic: we'll assume that if the variables have the same
           name, and differ only in respect of one being more complete
           than the other (having valid values where the other has
           missing values), then we should not create a new variable,
           but rather use the more complete series.

           Watch out for cases where this is not the desired behavior!
        */

        if (chk == VarComp.identical) {
            ret = TransformResult.var_exists_ok;
        } else if (chk == VarComp.new_has_missing) {
            /* is this right? */
            ret = TransformResult.var_exists_ok;
        } else if (chk == VarComp.old_has_missing) {
            for (t=0; t<dset.n; t++) {
                dset.Z[vnum][t] = x[t];
            }
            ret = TransformResult.var_exists_ok;
        } else {
            ret = transform_handle_duplicate(ci, lag, vnum, x, label, dset,
                                             parent, origv);
        }
    } else {
        /* no var of this name, working from scratch */
        if (dataset_add_series(dset, 1)) {
            ret = TransformResult.var_add_failed;
        } else {
            strcpy(dset.varname[vnum], vname);
            series_set_label(dset, vnum, label);
            for (t=0; t<dset.n; t++) {
                dset.Z[vnum][t] = x[t];
            }
        }
    }

    return ret;
}

private enum RESPECT_LAG_LEVEL = true;

private int get_lag_ID(int srcv, int lag, const(DatasetRaw) * dset)
{
    int fd = gretl_function_depth();
    const(char) * parent;
    const(char) * vname = dset.varname[srcv];
    int i;

    for (i=1; i<dset.v; i++) {
        if (RESPECT_LAG_LEVEL) {
            if (series_get_stack_level(dset, i) != fd) {
                continue;
            }
        }
        if (lag == series_get_lag(dset, i)) {
            parent = series_get_parent_name(dset, i);
            if (parent != null && !strcmp(vname, parent)) {
                return i;
            }
        }
    }

    return -1;
}

/* get_transform: create specified transformation of variable v if
   this variable does not already exist.

   The dummies resulting from DUMMIFY are automatically marked as
   discrete.

   origv is the number of variables in the dataset prior to the
   current round of adding.

   Return the ID number of the transformed var, or -1 on error.
*/

private int get_transform(int ci, int v, int aux, double x,
                          DatasetRaw * dset, int startlen, int origv,
                          const(int) * idxvec)
{
    char[VNAMELEN] vname = 0;
    char[MAXLABEL] label = 0;
    int vno = -1, err = 0;
    int len, lag = 0;
    const(char) * srcname;
    double * vx;

    vx = testvec(dset.n);
    if (vx == null) {
        return -1;
    }

    if (ci == GretlCmdIndex.lags) {
        lag = aux;
        err = get_lag(v, lag, vx, dset);
    } else if (ci == GretlCmdIndex.logs) {
        err = get_log(v, vx, dset);
    } else if (ci == GretlCmdIndex.diff || ci == GretlCmdIndex.ldiff ||
               ci == GretlCmdIndex.sdiff) {
        err = get_diff(v, vx, ci, dset);
    } else if (ci == GretlCmdIndex.orthdev) {
        err = get_orthdev(v, vx, dset);
    } else if (ci == GretlCmdIndex.square) {
        /* "aux" = second variable number */
        err = get_xpx(v, aux, vx, dset);
    } else if (ci == GretlCmdIndex.dummify) {
        /* "x" = value for dummy */
        err = get_discdum(v, x, vx, dset);
    } else if (ci == ExtraTransform.inverse) {
        err = get_inverse(v, vx, dset);
    } else if (ci == GretlCmdIndex.stdize) {
        err = get_std(v, aux, vx, dset);
    } else if (ci == ExtraTransform.resample) {
        err = get_resampled(v, vx, dset, idxvec);
    }

    if (err) {
        return -1;
    }

    if (ci == GretlCmdIndex.lags && (vno = get_lag_ID(v, aux, dset)) > 0) {
        /* special case: pre-existing lag */
        err = check_add_transform(ci, lag, vno, vx, vname.ptr, label.ptr, dset, v, origv);
        if (err != TransformResult.var_exists_ok) {
            vno = -1;
        }
        return vno;
    }

    if (ci == GretlCmdIndex.dummify || (ci == GretlCmdIndex.square && v != aux)) {
        /* special cases in respect of labeling */
        char * s = null;

        if (ci == GretlCmdIndex.square && v != aux) {
            s = g_strdup_printf(_("= %s times %s").ptr, dset.varname[v],
                                dset.varname[aux]);
        } else if (is_string_valued(dset, v)) {
            const(char) * sval = series_get_string_for_value(dset, v, x);

            if (sval != null && *sval != '\0') {
                s = g_strdup_printf(_("dummy for %s = '%s'").ptr, dset.varname[v], sval);
            } else {
                s = g_strdup_printf(_("dummy for %s = %g").ptr, dset.varname[v], x);
            }
        } else {
            s = g_strdup_printf(_("dummy for %s = %g").ptr, dset.varname[v], x);
        }
        strcpy(label.ptr, gretl_utf8_truncate(s, MAXLABEL-1));
        g_free(s);
    } else {
        make_transform_label(label.ptr, dset.varname[v], ci, aux);
    }

    srcname = get_mangled_name_by_id(v);
    if (srcname == null) {
        srcname = dset.varname[v];
    }

    for (len=startlen; len<=VNAMELEN; len++) {
        if (vname[0] != '\0' && len == VNAMELEN) {
            /* last resort: hack the name */
            make_varname_unique(vname.ptr, 0, dset);
        } else if (ci == GretlCmdIndex.square && v != aux) {
            make_xp_varname(vname.ptr, v, aux, dset, len);
        } else if (ci == GretlCmdIndex.dummify && series_is_integer_valued(dset, v)) {
            make_transform_varname(vname.ptr, srcname, ci, cast(int) x, len);
        } else {
            make_transform_varname(vname.ptr, srcname, ci, aux, len);
        }
        vno = series_index(dset, vname.ptr);
        err = check_add_transform(ci, lag, vno, vx, vname.ptr, label.ptr, dset, v, origv);
        if (err != TransformResult.var_exists_ok && err != TransformResult.var_added_ok) {
            vno = -1;
        }
        if (err != TransformResult.varname_duplicate) {
            /* either we're OK, or there was a fatal error */
            break;
        }
    }

    if (!err && vno > 0) {
        if (ci == GretlCmdIndex.dummify) {
            series_set_parent(dset, vno, dset.varname[v]);
            series_set_transform(dset, vno, ci);
            series_set_discrete(dset, vno, 1);
        } else if (ci == GretlCmdIndex.lags || ci == GretlCmdIndex.diff ||
                   ci == GretlCmdIndex.logs) {
            series_set_parent(dset, vno, dset.varname[v]);
            series_set_transform(dset, vno, ci);
        }

        if (ci == GretlCmdIndex.lags) {
            series_set_lag(dset, vno, aux);
        } else {
            series_set_lag(dset, vno, 0);
        }
    }

    return vno;
}

/**
 * laggenr:
 * @v: ID number in dataset of source variable.
 * @lag: the order of the lag to create.
 * @dset: dataset struct.
 *
 * Creates the specified lag of variable @v if this variable does
 * not already exist.
 *
 * Returns: the ID number of the lagged variable, or -1 on error.
 */

int laggenr(int v, int lag, DatasetRaw * dset)
{
    int lno;

    if (lag > dset.n || -lag > dset.n) {
        gretl_errmsg_sprintf(_("Invalid lag order %d").ptr, lag);
        lno = -1;
    } else if (v == 0) {
        gretl_errmsg_set(_("The constant cannot be lagged").ptr);
        lno = -1;
    } else if (lag == 0) {
        lno = v;
    } else {
        lno = get_transform(GretlCmdIndex.lags, v, lag, 0.0, dset,
                            VNAMELEN - 3, dset.v, null);
    }

    return lno;
}

/**
 * laggenr_from_to:
 * @v: ID number in dataset of source variable.
 * @fromlag: start of lag/lead range.
 * @tolag: end of lag/lead range.
 * @dset: dataset struct.
 * @err: location to receive error code.
 *
 * Creates the specified lags of variable @v if they do not
 * already exist.
 *
 * Returns: list of lag variables, or NULL or on error.
 */

int * laggenr_from_to(int v, int fromlag, int tolag,
                      DatasetRaw * dset, int * err)
{
    int * llist;
    int i, lv, nlags = -1;
    int p, count_down = 0;

    if (fromlag > tolag) {
        nlags = fromlag - tolag + 1;
        count_down = 1;
    } else {
        nlags = tolag - fromlag + 1;
    }

    if (nlags <= 0) {
        *err = GretlError.data;
        return null;
    }

    llist = gretl_list_new(nlags);
    if (llist == null) {
        *err = GretlError.alloc;
        return null;
    }

    p = fromlag;

    for (i=0; i<nlags; i++) {
        lv = laggenr(v, -p, dset);
        if (lv < 0) {
            *err = GretlError.data;
            free(llist);
            llist = null;
            break;
        }
        llist[i+1] = lv;
        if (count_down) {
            p--;
        } else {
            p++;
        }
    }

    return llist;
}

/**
 * loggenr:
 * @v: ID number in dataset of source variable.
 * @dset: dataset struct.
 *
 * Creates the natural log of variable @v if this variable does
 * not already exist.
 *
 * Returns: the ID number of the log variable, or -1 on error.
 */

int loggenr(int v, DatasetRaw * dset)
{
    return get_transform(GretlCmdIndex.logs, v, 0, 0.0, dset,
                         VNAMELEN - 3, dset.v, null);
}

/**
 * invgenr:
 * @v: ID number in dataset of source variable.
 * @dset: dataset struct.
 *
 * Creates the reciprocal of variable @v if this variable does
 * not already exist.
 *
 * Returns: the ID number of the reciprocal, or -1 on error.
 */

int invgenr(int v, DatasetRaw * dset)
{
    return get_transform(ExtraTransform.inverse, v, 0, 0.0, dset,
                         VNAMELEN - 3, dset.v, null);
}

/**
 * diffgenr:
 * @v: ID number in dataset of source variable.
 * @ci: DIFF (first difference), LDIFF (log difference) or SDIFF
 * (seasonal difference).
 * @dset: dataset struct.
 *
 * Creates the first difference (or log- or seasonal difference,
 * depending on the value of @ci) of variable @v, if the
 * differenced variable does not already exist.
 *
 * Returns: the ID number of the differenced variable, or -1 on error.
 */

int diffgenr(int v, GretlCmdIndex ci, DatasetRaw * dset)
{
    if (ci != GretlCmdIndex.diff && ci != GretlCmdIndex.ldiff &&
        ci != GretlCmdIndex.sdiff) {
        return -1;
    }

    if (ci == GretlCmdIndex.sdiff && !dataset_is_seasonal(dset)) {
        return -1;
    }

    return get_transform(ci, v, 0, 0.0, dset,
                         VNAMELEN - 3, dset.v, null);
}

/**
 * xpxgenr:
 * @vi: ID number in dataset of first source variable.
 * @vj: ID number in dataset of second source variable.
 * @dset: dataset struct.
 *
 * Creates the cross product of variables @vi and @vj if this
 * variable does not already exist.
 *
 * Returns: the ID number of the cross-product variable,
 * or -1 on error.
 */

int xpxgenr(int vi, int vj, DatasetRaw * dset)
{
    if (vi == vj) {
        if (gretl_isdummy(dset.t1, dset.t2, dset.Z[vi])) {
            return -1;
        }
    }

    return get_transform(GretlCmdIndex.square, vi, vj, 0.0, dset,
                         VNAMELEN - 3, dset.v, null);
}

private int get_starting_length(const(int) * list, DatasetRaw * dset, int trim)
{
    int width = VNAMELEN - 3 - trim;
    const(char) * vni;
    const(char) * vnj;
    int len, maxlen = 0;
    int conflict = 0;
    int i, j;

    if (list[0] == 1) {
        return VNAMELEN - 1;
    }

    for (i=1; i<=list[0]; i++) {
        len = cast(int) strlen(dset.varname[list[i]]);
        if (len > maxlen) {
            maxlen = len;
        }
    }

    if (maxlen <= width) {
        /* no problem: generated names will fit in VNAMELEN - 3 chars */
        return VNAMELEN - 3;
    }

    for (len=width; len<=maxlen; len++) {
        conflict = 0;
        for (i=1; i<=list[0] && !conflict; i++) {
            vni = dset.varname[list[i]];
            for (j=i+1; j<=list[0] && !conflict; j++) {
                vnj = dset.varname[list[j]];
                if (!strncmp(vni, vnj, len)) {
                    conflict = 1;
                }
            }
        }
        if (!conflict) {
            break;
        }
    }

    len += trim;

    if (len < VNAMELEN - 3) {
        len = VNAMELEN - 3;
    } else if (len > VNAMELEN - 1) {
        len = VNAMELEN - 1;
    }

    return len;
}

private struct MangledName {
    int v;
    char * s;
}

private __gshared MangledName * mnames;
private __gshared int n_mnames;

private void destroy_mangled_names()
{
    int i;

    for (i=0; i<n_mnames; i++) {
        free(mnames[i].s);
    }
    free(mnames);
    mnames = null;
    n_mnames = 0;
}

private int push_mangled_name(int v, const(char) * s)
{
    MangledName * tmp;

    tmp = cast(MangledName *) realloc(mnames, (n_mnames + 1) * (*tmp).sizeof);

    if (tmp == null) {
        destroy_mangled_names();
        return GretlError.alloc;
    }

    mnames = tmp;
    mnames[n_mnames].v = v;
    mnames[n_mnames].s = gretl_strdup(s);

    n_mnames++;

    return 0;
}

private char * get_mangled_name_by_id(int v)
{
    int i;

    for (i=0; i<n_mnames; i++) {
        if (mnames[i].v == v) {
            return mnames[i].s;
        }
    }

    return null;
}

private int make_mangled_name(int v, const(char) * s, int nc)
{
    char[32] tmp;
    char[32] sfx;
    int n, err;

    if (get_mangled_name_by_id(v)) {
        return 0;
    }

    tmp[0] = sfx[0] = '\0';
    /* we know the bit out here (+nc) must be distinct? */
    strcat(sfx.ptr, s + nc);
    n = cast(int) strlen(sfx.ptr);
    strncat(tmp.ptr, s, nc - n);
    strcat(tmp.ptr, sfx.ptr);

    err = push_mangled_name(v, tmp.ptr);

    return err;
}

private enum XPROD = GretlCmdIndex.nc + 1;

private int transform_preprocess_list(int * list, const(DatasetRaw) * dset, int f)
{
    int maxc = VNAMELEN - 3;
    int longnames = 0;
    char ** S = null;
    int i, v, ok;
    int err = 0;

    if (f == GretlCmdIndex.square || f == GretlCmdIndex.ldiff ||
        f == GretlCmdIndex.sdiff || f == ExtraTransform.resample) {
        /* 3-character prefixes */
        maxc--;
    } else if (f == GretlCmdIndex.dummify) {
        /* "D" plus suffix */
        maxc--;
    }

    for (i=1; i<=list[0]; i++) {
        v = list[i];
        ok = 1;

        if (f == GretlCmdIndex.logs) {
            if (v == 0) { /* 1 -> 0: FIXME?? */
                ok = 1;
            }
            if (gretl_isdummy(dset.t1, dset.t2, dset.Z[v])) {
                ok = 0;
            }
        } else if (f == GretlCmdIndex.square) {
            if (gretl_isdummy(dset.t1, dset.t2, dset.Z[v])) {
                ok = 0;
            }
        } else if (f == GretlCmdIndex.lags) {
            if (v == 0) {
                gretl_errmsg_set(_("The constant cannot be lagged").ptr);
                ok = 0;
            }
        } else if (f == GretlCmdIndex.dummify) {
            ok = v > 0 && accept_as_discrete(dset, v, 0);
        }
        if (!ok) {
            gretl_list_delete_at_pos(list, i--);
        } else if (strlen(dset.varname[v]) > maxc) {
            strings_array_add(&S, &longnames, dset.varname[v]);
        }
    }

    if (list[0] == 0) {
        err = GretlError.types;
    }

    if (longnames > 0) {
        if (!err) {
            int j, herr = 0;

            for (i=0; i<longnames && !herr; i++) {
                for (j=i+1; j<longnames && !herr; j++) {
                    if (!strncmp(S[i], S[j], maxc)) {
                        v = series_index(dset, S[i]);
                        herr = make_mangled_name(v, S[i], maxc);
                        v = series_index(dset, S[j]);
                        herr += make_mangled_name(v, S[j], maxc);
                    }
                }
            }
        }

        strings_array_free(S, longnames);
    }

    return err;
}

/**
 * list_loggenr:
 * @list: on entry, list of variables to process; on exit,
 * holds the ID numbers of the generated variables.
 * @dset: dataset struct.
 *
 * Generates and adds to the data set the natural logs of the
 * variables given in @list.
 *
 * Returns: 0 on success, error code on error.
 */

int list_loggenr(int * list, DatasetRaw * dset)
{
    int origv = dset.v;
    int tnum, i, j, v;
    int startlen;
    int l0 = 0;
    int err;

    err = transform_preprocess_list(list, dset, GretlCmdIndex.logs);
    if (err) {
        return err;
    }

    startlen = get_starting_length(list, dset, 2);

    j = 1;
    for (i=1; i<=list[0]; i++) {
        v = list[i];
        tnum = get_transform(GretlCmdIndex.logs, v, 0, 0.0, dset, startlen,
                             origv, null);
        if (tnum > 0) {
            list[j++] = tnum;
            l0++;
        }
    }

    list[0] = l0;

    destroy_mangled_names();

    return (l0 > 0)? 0 : GretlError.logs;
}

/**
 * list_stdgenr:
 * @list: on entry, list of variables to process; on exit,
 * holds the ID numbers of the generated variables.
 * @dset: dataset struct.
 *
 * Generates and adds to the dataset standardized versions
 * of the series given in @list.
 *
 * Returns: 0 on success, error code on error.
 */

int list_stdgenr(int * list, DatasetRaw * dset, Gretlopt opt)
{
    int origv = dset.v;
    int tnum, i, j, v;
    int startlen;
    int dfc = 1;
    int l0 = 0;
    int err;

    err = transform_preprocess_list(list, dset, GretlCmdIndex.stdize);
    if (err) {
        return err;
    }

    if (opt & Gretlopt.OPT_C) {
        dfc = -1;
    } else if (opt & Gretlopt.OPT_N) {
        dfc = 0;
    }

    startlen = get_starting_length(list, dset, 2);

    j = 1;
    for (i=1; i<=list[0]; i++) {
        v = list[i];
        tnum = get_transform(GretlCmdIndex.stdize, v, dfc, 0.0, dset, startlen,
                             origv, null);
        if (tnum > 0) {
            list[j++] = tnum;
            l0++;
        }
    }

    list[0] = l0;

    destroy_mangled_names();

    return (l0 > 0)? 0 : GretlError.data;
}

private int * make_lags_list(int * list, int order)
{
    int i, v, nl = 0;

    for (i=1; i<=list[0]; i++) {
        v = list[i];
        if (v > 0) {
            nl += order;
        }
    }

    return gretl_list_new(nl);
}

private int lag_wanted(int p, const(GretlMatrixRaw) * v, int n)
{
    int ret = 1;

    if (v != null) {
        int i;

        ret = 0; /* reverse the assumption */
        for (i=0; i<n; i++) {
            if (v.val[i] == p) {
                ret = 1;
                break;
            }
        }
    }

    return ret;
}

private int process_hf_lags_input(int hf_min, int hf_max,
                                  DatasetRaw * dset,
                                  int compfac,
                                  int * lf_min, int * lf_max,
                                  int * skip_first,
                                  int * n_terms)
{
    if (dset.pd != 1 && dset.pd != 12 && dset.pd != 4) {
        return GretlError.pdwrong;
    }

    *lf_min = cast(int) ceil(hf_min / cast(double) compfac);
    *lf_max = cast(int) ceil(hf_max / cast(double) compfac);
    *skip_first = (hf_min - 1) % compfac;
    *n_terms = hf_max - hf_min + 1;

    if (*skip_first < 0) {
        /* handle leads */
        *skip_first = compfac + *skip_first;
    }

    return 0;
}

/**
 * list_laggenr:
 * @plist: on entry, pointer to list of variables to process.  On exit
 * the list holds the ID numbers of the lag variables.
 * @lmin: minimum lag to include (defaults to 1).
 * @lmax: maximum lag to include (or 0 for automatic).
 * @lvec: (alternative to @lmin, @lmax) vector holding lags to generate.
 * @dset: dataset struct.
 * @compfac: compaction factor (for MIDAS lists only, otherwise give 0).
 * @opt: may contain OPT_L to order the list by lag rather than
 * by variable.
 *
 * Generates and adds to the data set lagged values of the
 * variables given in the list pointed to by @plist.
 *
 * Returns: 0 on successful completion, 1 on error.
 */

int list_laggenr(int ** plist, int lmin, int lmax,
                 const(GretlMatrixRaw) * lvec,
                 DatasetRaw * dset, int compfac,
                 Gretlopt opt)
{
    int origv = dset.v;
    int * list = *plist;
    int * laglist = null;
    int l, i, j, v, lv;
    int startlen, l0 = 0;
    int skip_first = 0;
    int veclen = 0;
    int n_terms = 0;
    int err = 0;

    if (compfac < 0) {
        return GretlError.invarg;
    } else if (compfac > 0) {
        /* MIDAS: we want @lmin and @lmax, not @lvec */
        if (lvec != null) {
            return GretlError.invarg;
        }
        if (!gretl_is_midas_list(list, dset)) {
            gretl_warnmsg_set(_("The argument does not seem to be a MIDAS list").ptr);
        }
    }

    if (lvec != null) {
        /* got a vector of lags in @lvec */
        veclen = gretl_vector_get_length(lvec);
        if (veclen == 0) {
            err = GretlError.invarg;
        } else {
            lmin = cast(int) lvec.val[0];
            lmax = cast(int) lvec.val[veclen-1];
            n_terms = veclen;
        }
    } else if (compfac > 0) {
        /* the MIDAS case */
        int hmin = lmin, hmax = lmax;

        if (hmax < hmin) {
            err = GretlError.invarg;
        } else {
            err = process_hf_lags_input(hmin, hmax, dset, compfac, &lmin,
                                        &lmax, &skip_first, &n_terms);
        }
    } else {
        /* non-MIDAS, using order = @lmax, ignoring @lmin */
        if (lmax < 0 || lmax > dset.n) {
            gretl_errmsg_sprintf(_("Invalid lag order %d").ptr, lmax);
            return GretlError.invarg;
        } else if (lmax == 0) {
            lmax = default_lag_order(dset);
        }
        lmin = 1;
        n_terms = lmax;
    }

    if (err) {
        return err;
    }

    /* Note: at this point @lmin and @lmax are defined in terms of
       the frequency of @dset, even when we're handling MIDAS data.
       In the regular case, @n_terms holds the number of lag terms
       wanted per series in the input @list, but in the MIDAS case
       @n_terms is the total number of hf lag terms wanted (since
       the input @list actually represents a single series).
    */

    err = transform_preprocess_list(list, dset, GretlCmdIndex.lags);
    if (err) {
        return err;
    }

    if (compfac) {
        laglist = gretl_list_new(n_terms);
    } else {
        laglist = make_lags_list(list, n_terms);
    }

    if (laglist == null) {
        destroy_mangled_names();
        return GretlError.alloc;
    }

    startlen = get_starting_length(list, dset, (lmax > 9)? 3 : 2);

    j = 1;

    if (compfac > 0) {
        /* MIDAS high-frequency lags, by lags */
        int hfpmax = n_terms + skip_first;
        int mp, hfp = 0;

        for (l=lmin; l<=lmax; l++) {
            for (i=1; i<=list[0]; i++) {
                hfp++;
                if (hfp <= skip_first) {
                    continue;
                } else if (hfp > hfpmax) {
                    break;
                }
                v = list[i];
                lv = get_transform(GretlCmdIndex.lags, v, l, 0.0, dset,
                                   startlen, origv, null);
                if (lv > 0) {
                    mp = series_get_midas_period(dset, v);
                    series_set_midas_period(dset, lv, mp);
                    mp = series_get_midas_freq(dset, v);
                    series_set_midas_freq(dset, lv, mp);
                    laglist[j++] = lv;
                    l0++;
                }
            }
        }
    } else if (opt & Gretlopt.OPT_L) {
        /* order the list by lags */
        for (l=lmin; l<=lmax; l++) {
            if (!lag_wanted(l, lvec, veclen)) {
                continue;
            }
            for (i=1; i<=list[0]; i++) {
                v = list[i];
                lv = get_transform(GretlCmdIndex.lags, v, l, 0.0, dset,
                                   startlen, origv, null);
                if (lv > 0) {
                    laglist[j++] = lv;
                    l0++;
                }
            }
        }
    } else {
        /* order by variable */
        for (i=1; i<=list[0]; i++) {
            v = list[i];
            for (l=lmin; l<=lmax; l++) {
                if (!lag_wanted(l, lvec, veclen)) {
                    continue;
                }
                lv = get_transform(GretlCmdIndex.lags, v, l, 0.0, dset,
                                   startlen, origv, null);
                if (lv > 0) {
                    laglist[j++] = lv;
                    l0++;
                }
            }
        }
    }

    destroy_mangled_names();

    laglist[0] = l0;

    free(*plist);
    *plist = laglist;

    return 0;
}

/**
 * default_lag_order:
 * @dset: data information struct.
 *
 * Returns: the default lag order for generating lags, performing
 * autocorrelation test, and so on.
 */

int default_lag_order(const(DatasetRaw) * dset)
{
    int order = 1;

    if (!dataset_is_panel(dset)) {
        order = (dset.pd < 52)? dset.pd : 14;
    }

    return order;
}

/**
 * list_diffgenr:
 * @list: on entry, list of variables to process; on exit,
 * ID numbers of the generated variables.
 * @ci: must be DIFF, LDIFF or SDIFF.
 * @dset: pointer to dataset struct.
 *
 * Generate differences of the variables in @list, and add them
 * to the data set.  If @ci is DIFF these are ordinary first
 * differences; if @ci is LDIFF they are log differences; and
 * if @ci is SDIFF they are seasonal differences.
 *
 * Returns: 0 on successful completion, non-zero on error.
 */

int list_diffgenr(int * list, GretlCmdIndex ci, DatasetRaw * dset)
{
    int origv = dset.v;
    int i, v, startlen;
    int tnum, l0 = 0;
    int err;

    if (list[0] == 0) {
        return 0;
    }

    if (ci != GretlCmdIndex.diff && ci != GretlCmdIndex.ldiff &&
        ci != GretlCmdIndex.sdiff) {
        return 1;
    }

    if (ci == GretlCmdIndex.sdiff && !dataset_is_seasonal(dset)) {
        return GretlError.pdwrong;
    }

    err = transform_preprocess_list(list, dset, ci);
    if (err) {
        return err;
    }

    startlen = get_starting_length(list, dset, (ci == GretlCmdIndex.diff)? 2 : 3);

    for (i=1; i<=list[0] && !err; i++) {
        v = list[i];
        tnum = get_transform(ci, v, 0, 0.0, dset,
                             startlen, origv, null);
        if (tnum < 0) {
            err = 1;
        } else {
            list[i] = tnum;
            l0++;
        }
    }

    list[0] = l0;

    destroy_mangled_names();

    return err;
}

private int check_hf_difflist(const(int) * list,
                              const(DatasetRaw) * dset,
                              int ci, int * n_add)
{
    char[VNAMELEN] vname;
    int i, v, li, err = 0;

    for (i=1; i<=list[0] && !err; i++) {
        li = list[i];
        make_transform_varname(vname.ptr, dset.varname[li],
                               ci, 0, VNAMELEN - 8);
        v = current_series_index(dset, vname.ptr);
        if (v > 0) {
            *n_add -= 1;
        } else if (gretl_is_user_var(vname.ptr)) {
            gretl_errmsg_sprintf(_("%s: collides with existing object name").ptr,
                                 vname.ptr);
            err = GretlError.types;
        }
    }

    return err;
}

/**
 * hf_list_diffgenr:
 * @list: on entry, midas list of variables to process; on exit,
 * ID numbers of the generated variables.
 * @ci: must be DIFF or LDIFF.
 * @parm: optional scalar multiplier.
 * @dset: pointer to dataset struct.
 *
 * Generate high-frequency differences of the variables in @list
 * and add them to the data set.  If @ci is DIFF these are ordinary
 * first differences; if @ci is LDIFF they are log differences.
 * If @parm is not NA then the new series are all multiplied by
 * the specified value (as in multiplication of log-differences
 * by 100).
 *
 * Returns: 0 on successful completion, non-zero on error.
 */

int hf_list_diffgenr(int * list, GretlCmdIndex ci, double parm, DatasetRaw * dset)
{
    GretlMatrixRaw * dX = null;
    int n_add, hfci = 0;
    int set_midas = 1;
    int err = 0;

    if (list[0] == 0) {
        return 0;
    }

    if (ci == GretlCmdIndex.diff) {
        hfci = ExtraTransform.hfdiff;
    } else if (ci == GretlCmdIndex.ldiff) {
        hfci = ExtraTransform.hfldiff;
    } else {
        return GretlError.invarg;
    }

    if (!dataset_is_time_series(dset)) {
        return GretlError.pdwrong;
    }

    n_add = list[0];
    err = check_hf_difflist(list, dset, ci, &n_add);
    if (err) {
        return err;
    }

    if (!gretl_is_midas_list(list, dset)) {
        gretl_warnmsg_set(_("The argument does not seem to be a MIDAS list").ptr);
        set_midas = 0;
    }

    if (na(parm)) {
        /* set to default */
        parm = 1.0;
    }

    dX = get_hf_diffs(list, ci, parm, dset, &err);

    if (!err) {
        char[VNAMELEN] vname;
        char[MAXLABEL] label;
        const(char) * parent;
        int i, s, t, n = dX.cols;
        int vi, v0 = dset.v;
        int offset = 0;

        err = dataset_add_series(dset, n_add);

        if (!err) {
            for (i=0; i<n; i++) {
                parent = dset.varname[list[i+1]];
                make_transform_varname(vname.ptr, parent, ci, 0,
                                       VNAMELEN - 8);
                vi = current_series_index(dset, vname.ptr);
                if (vi < 0) {
                    vi = v0 + offset;
                    offset++;
                    strcpy(dset.varname[vi], vname.ptr);
                    make_transform_label(label.ptr, parent, hfci, 0);
                    series_record_label(dset, vi, label.ptr);
                }
                s = 0;
                for (t=dset.t1; t<=dset.t2; t++) {
                    dset.Z[vi][t] = gretl_matrix_get(dX, s++, i);
                }
                list[i+1] = vi;
            }
            list[0] = n;
        }
        gretl_matrix_free(dX);
    }

    if (!err && set_midas) {
        gretl_list_set_midas(list, dset);
    }

    return err;
}

/**
 * list_orthdev:
 * @list: list of variables to process.
 * @dset: dataset struct.
 *
 * Generate orthogonal deviations of the variables in @list, and add
 * them to the data set.
 *
 * Returns: 0 on success, error code on error.
 */

int list_orthdev(int * list, DatasetRaw * dset)
{
    int origv = dset.v;
    int i, v, startlen;
    int tnum, l0 = 0;
    int err;

    if (list[0] == 0) {
        return 0;
    }

    if (!dataset_is_panel(dset)) {
        return GretlError.pdwrong;
    }

    err = transform_preprocess_list(list, dset, GretlCmdIndex.orthdev);
    if (err) {
        return err;
    }

    startlen = get_starting_length(list, dset, 2);

    for (i=1; i<=list[0] && !err; i++) {
        v = list[i];
        tnum = get_transform(GretlCmdIndex.orthdev, v, 0, 0.0, dset,
                             startlen, origv, null);
        if (tnum < 0) {
            err = 1;
        } else {
            list[i] = tnum;
            l0++;
        }
    }

    list[0] = l0;

    destroy_mangled_names();

    return err;
}

/**
 * list_xpxgenr:
 * @plist: pointer to list of variables to process.  On exit
 * the list holds the ID numbers of the squares (and possibly
 * cross-products).
 * @dset: dataset struct.
 * @opt: If OPT_O, both squares and cross-products are generated,
 * otherwise only squares.
 *
 * Generates and adds to the data set squares and (if @opt is OPT_O)
 * cross-products of the variables given in the list pointed to
 * by @plist.
 *
 * Returns: 0 on success, error code on error.
 */

int list_xpxgenr(int ** plist, DatasetRaw * dset, Gretlopt opt)
{
    int origv = dset.v;
    int * list = *plist;
    int * xpxlist = null;
    int tnum, i, j, k, vi, vj;
    int startlen, l0;
    int f, err;

    f = (opt & Gretlopt.OPT_O) ? XPROD : GretlCmdIndex.square;

    err = transform_preprocess_list(list, dset, f);
    if (err) {
        return err;
    }

    l0 = list[0];

    if (opt & Gretlopt.OPT_O) {
        int maxterms = (l0 * l0 + l0) / 2;

        xpxlist = gretl_list_new(maxterms);
        if (xpxlist == null) {
            destroy_mangled_names();
            return GretlError.alloc;
        }
    } else {
        xpxlist = list;
    }

    startlen = get_starting_length(list, dset, 3);
    xpxlist[0] = 0;

    k = 1;
    for (i=1; i<=l0; i++) {
        vi = list[i];
        /* don't square dummies, but compute interactions with other series
           if needed
        */
        if (!gretl_isdummy(dset.t1, dset.t2, dset.Z[vi])) {
            tnum = get_transform(GretlCmdIndex.square, vi, vi, 0.0, dset,
                                 startlen, origv, null);
            if (tnum > 0) {
                xpxlist[k++] = tnum;
                xpxlist[0] += 1;
            }
        }

        if (opt & Gretlopt.OPT_O) {
            for (j=i+1; j<=l0; j++) {
                vj = list[j];
                tnum = xpxgenr(vi, vj, dset);
                if (tnum > 0) {
                    xpxlist[k++] = tnum;
                    xpxlist[0] += 1;
                }
            }
        }
    }

    destroy_mangled_names();

    if (opt & Gretlopt.OPT_O) {
        free(*plist);
        *plist = xpxlist;
    }

    return (xpxlist[0] > 0)? 0 : GretlError.squares;
}

int list_resample(int * list, DatasetRaw * dset)
{
    int origv = dset.v;
    int * z = null;
    int tnum, i, j, vi;
    int n, startlen, l0;
    int k1, k2;
    int err;

    err = transform_preprocess_list(list, dset, ExtraTransform.resample);
    if (err) {
        return err;
    }

    if (dataset_is_panel(dset)) {
        k1 = dset.t1 / dset.pd;
        k2 = dset.t2 / dset.pd;
    } else {
        k1 = dset.t1;
        k2 = dset.t2;
    }

    n = k2 - k1 + 1;

    z = cast(int *) malloc(n * (*z).sizeof);
    if (z == null) {
        return GretlError.alloc;
    }

    /* generate uniform drawings from sample range */
    gretl_rand_int_minmax(z, n, k1, k2);

    l0 = list[0];
    startlen = get_starting_length(list, dset, 3);
    list[0] = 0;

    j = 1;
    for (i=1; i<=l0; i++) {
        vi = list[i];
        tnum = get_transform(ExtraTransform.resample, vi, 0, 0.0, dset,
                             startlen, origv, z);
        if (tnum > 0) {
            list[j++] = tnum;
            list[0] += 1;
        }
    }

    free(z);
    destroy_mangled_names();

    return (list[0] > 0)? 0 : GretlError.data;
}

/**
 * list_dropcoll:
 * @list: on entry, list of variables to process; on exit,
 * the original list minus collinear terms.
 * @dset: dataset struct.
 *
 * Drop collinear variables from @list.
 *
 * Returns: 0 on success, error code on error.
 */

int list_dropcoll(int * list, double eps, DatasetRaw * dset)
{
    GretlMatrixRaw * R = null;
    GretlMatrixRaw * X = null;
    const(double) * x;
    int i, T, k;
    int err = 0;

    if (list == null) {
        return GretlError.data;
    } else if (list[0] < 2) {
        /* no-op */
        return 0;
    }

    if (eps < 0) {
        return GretlError.data;
    } else if (na(eps)) {
        eps = R_DIAG_MIN;
    }

    for (i=list[0]; i>0; i--) {
        x = dset.Z[list[i]];
        if (gretl_iszero(dset.t1, dset.t2, x)) {
            err = gretl_list_delete_at_pos(list, i);
        }
    }

    X = gretl_matrix_data_subset(list, dset, dset.t1, dset.t2,
                                 MMissingCode.skip, &err);
    if (err) {
        return err;
    }

    T = X.rows;
    k = X.cols;

    if (T < k) {
        gretl_errmsg_sprintf(_("A minimum of %d observations is required").ptr, k);
        err = GretlError.toofew;
    } else {
        R = gretl_matrix_alloc(k, k);
        if (R == null) {
            err = GretlError.alloc;
        } else {
            err = gretl_matrix_QR_decomp(X, R);
        }
    }

    if (!err) {
        double rii;
        int j = 1;

        for (i=0; i<k && !err; i++, j++) {
            rii = gretl_matrix_get(R, i, i);
            if (fabs(rii) < eps) {
                err = gretl_list_delete_at_pos(list, j--);
            }
        }
    }

    gretl_matrix_free(X);
    gretl_matrix_free(R);

    return err;
}

private int real_list_dumgenr(int ** plist, DatasetRaw * dset,
                              double oddval, Gretlopt opt)
{
    int origv = dset.v;
    int * list = *plist;
    int * tmplist = null;
    double * x = null;
    int i, j, t, n;
    int startlen;
    int err;

    err = transform_preprocess_list(list, dset, GretlCmdIndex.dummify);
    if (err) {
        return err;
    }

    tmplist = gretl_null_list();
    if (tmplist == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    x = cast(double *) malloc(dset.n * (*x).sizeof);
    if (x == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    startlen = get_starting_length(list, dset, 3);

    for (i=1; i<=list[0] && !err; i++) {
        int vi = list[i];
        int nuniq, tnum;
        int jmin, jmax;
        double xt;

        n = 0;
        for (t=dset.t1; t<=dset.t2; t++) {
            xt = dset.Z[vi][t];
            if (!na(xt)) {
                x[n++] = xt;
            }
        }
        if (n == 0) {
            continue;
        }

        qsort(x, n, (*x).sizeof, &gretl_compare_doubles);
        nuniq = count_distinct_values(x, n);
        if (nuniq == 1) {
            continue;
        }

        rearrange_id_array(x, nuniq, n);
        jmin = (opt & Gretlopt.OPT_F)? 1 : 0;
        jmax = (opt & Gretlopt.OPT_L)? nuniq - 1 : nuniq;

        for (j=jmin; j<jmax && !err; j++) {
            if (na(oddval) || x[j] != oddval) {
                tnum = get_transform(GretlCmdIndex.dummify, vi, j+1, x[j], dset,
                                     startlen, origv, null);
                if (tnum > 0) {
                    tmplist = gretl_list_append_term(&tmplist, tnum);
                    if (tmplist == null) {
                        err = GretlError.alloc;
                    }
                } else {
                    err = GretlError.data;
                }
            }
        }
    }

    if (!err && !(opt & Gretlopt.OPT_A) && tmplist[0] == 0) {
        /* don't fail here if doing auto-dummify (OPT_A) */
        gretl_errmsg_set(_("dummify: no suitable variables were found").ptr);
        err = GretlError.data;
    }

    free(x);

 bailout:

    if (!err && tmplist[0] > 0) {
        free(*plist);
        *plist = tmplist;
    } else {
        free(tmplist);
    }

    destroy_mangled_names();

    return err;
}

/**
 * list_dumgenr:
 * @plist: pointer to list of variables to process; on exit
 * the list holds the ID numbers of the generated dummies.
 * @dset: dataset struct.
 * @opt: can include OPT_F to drop the first value, OPT_L to drop
 * the last value.
 *
 * For each of the variables given in the list to which @plist
 * points, generates and adds to the data set k dummy variables
 * coding for the k distinct values of the variable in question.
 * All these variables must have already been marked as discrete.
 * If the OPT_F or OPT_L option is given, either the first or
 * the last value of each variable is taken as the "base" and is
 * not given a dummy encoding (that is, only k - 1 dummies are
 * added for each variable).
 *
 * Returns: 0 on success, error code on error.
*/

int list_dumgenr(int ** plist, DatasetRaw * dset, Gretlopt opt)
{
    return real_list_dumgenr(plist, dset, NADBL, opt);
}

/**
 * dumgenr_with_oddval:
 * @plist: pointer to list of variables to process; on exit
 * the list holds the ID numbers of the generated dummies.
 * @dset: dataset struct.
 * @oddval: value which should be skipped when encoding the
 * input values as dummies.
 *
 * For each of the variables given in the list to which @plist
 * points, generates and adds to the data set k dummy variables
 * coding for the k distinct values of the variable in question.
 * All these variables must have already been marked as discrete.
 * if @oddval is not %NADBL, it is treated as the omitted
 * category and only k - 1 dummies are added for each variable).
 *
 * Returns: 0 on success, error code on error.
*/

int dumgenr_with_oddval(int ** plist, DatasetRaw * dset, double oddval)
{
    return real_list_dumgenr(plist, dset, oddval, Gretlopt.OPT_NONE);
}

/**
 * auto_dummify_list:
 * @plist: pointer to list of variables to process.
 * @dset: dataset struct.
 *
 * Produces a list in which any coded series in the original
 * list are replaced by a set of dummy variables. If the
 * original list contains no such series it is not altered.
 * In generating dummy series, the omitted category is the
 * minimum value of the coding variable.
 *
 * Returns: 0 on success, error code on error.
*/

int auto_dummify_list(int ** plist, DatasetRaw * dset)
{
    int * orig = *plist;
    int * list = null;
    int * dlist = null;
    int i, vi;
    int err = 0;

    /* empty starting point for full list */
    list = gretl_null_list();

    for (i=1; i<=orig[0] && !err; i++) {
        vi = orig[i];
        if (series_is_coded(dset, vi)) {
            /* coded series: split into dummies */
            dlist = gretl_list_new(1);
            dlist[1] = vi;
            err = real_list_dumgenr(&dlist, dset, NADBL, Gretlopt.OPT_F | Gretlopt.OPT_A);
            if (!err) {
                err = gretl_list_append_list(&list, dlist);
            }
            free(dlist);
        } else if (!in_gretl_list(list, vi)) {
            /* non-coded, not already present */
            list = gretl_list_append_term(&list, vi);
            if (list == null) {
                err = GretlError.alloc;
            }
        }
    }

    if (err) {
        free(list);
    } else if (list != *plist) {
        /* replace original list */
        free(*plist);
        *plist = list;
    }

    return err;
}

/**
 * list_makediscrete:
 * @list: list of variables to process.
 * @dset: data information struct.
 * @opt: if OPT_R, reverse the operation.
 *
 * Sets the variables given in @list as discrete, unless
 * opt is OPT_R, in which case the variables are set as
 * continuous.
 *
 * Returns: 0 on success, error code on error.
 */

int list_makediscrete(const(int) * list, DatasetRaw * dset, Gretlopt opt)
{
    int disc = !(opt & Gretlopt.OPT_R);
    int i, v, err = 0;

    for (i=1; i<=list[0]; i++) {
        v = list[i];
        if (v > 0) {
            series_set_discrete(dset, v, disc);
        }
    }

    return err;
}
