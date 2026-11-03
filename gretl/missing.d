/* Port of missing.c */

module gretl.missing;


import core.stdc.stdlib : calloc, free, malloc;
import core.stdc.string : memset;

import gretl.dataset : series_is_hidden;
import gretl.gretl_errors : gretl_errmsg_sprintf;
import gretl.gretl_listh : LISTSEP;
import gretl.gretl_model : gretl_model_get_int;
import gretl.gretl_prn : Prn, pprintf, pputs;
import gretl.libgretlh;
import gretl.missingh : NADBL, na;
import gretl.gretl_errorsh : GretlError;
import gretl.strutils : gretl_strdup;
import gretl.gretl_utils : gretl_double_from_string;

/* ---- declarations of not-yet-ported items used below ----
   Replace each with an import when its defining file is ported. */

extern (C) {
    int object_is_const(const(char) * name, int vnum);
}

/* mechanisms for handling missing observations */

/**
 * model_has_missing_obs:
 * @pmod: pointer to model.
 *
 * Returns: 1 if there are missing observations in the
 * model's sample range, otherwise 0.
 */

int model_has_missing_obs(const(Model) * pmod)
{
    int sample = pmod.t2 - pmod.t1 + 1;

    return (pmod.nobs < sample);
}

private int really_missing(int v, int t, const(double) ** Z, int d)
{
    if (d > 0 && Z[d][t] == 0) {
        /* the obs is dummied out */
        return 0;
    } else {
        return na(Z[v][t]);
    }
}

int model_missing(const(Model) * pmod, int t)
{
    if (pmod.missmask != null) {
        return pmod.missmask[t] == '1';
    } else {
        return 0;
    }
}

/* Construct a mask to code for missing observations within the sample
   range for a model.  The mask is an array of char with '0's
   for OK observations and '1's for missing obs, terminated with
   a NUL byte. This format allows copying of the mask using
   strcpy.
*/

private char * model_missmask(const(int) * list, int t1, int t2,
                              int n, const(double) ** Z, int dwt,
                              int * misscount)
{
    char * mask;
    int i, vi, t;

    mask = cast(char *) malloc(n + 1);
    if (mask == null) {
        return null;
    }

    memset(mask, '0', n);
    mask[n] = 0; /* note NUL-termination */

    if (misscount != null) {
        *misscount = 0;
    }

    for (t=t1; t<=t2; t++) {
        for (i=1; i<=list[0]; i++) {
            vi = list[i];
            if (vi > 0 && vi != LISTSEP) {
                if (really_missing(vi, t, Z, dwt)) {
                    /* FIXME dwt case and nobs? */
                    mask[t] = '1';
                    if (misscount != null) {
                        *misscount += 1;
                    }
                    break;
                }
            }
        }
    }

    return mask;
}

/**
 * model_adjust_sample:
 * @pmod: pointer to gretl model.
 * @n: full length of data array.
 * @Z: data array.
 * @misst: location to receive the first observation with a
 * missing value inside the sample range, or NULL.
 *
 * Drops leading or trailing observations from the sample range
 * initially given by the values in the t1 and t2 members of
 * @pmod, if missing values are found among the variables given
 * in the list member of @pmod.
 *
 * Also checks for missing values within the adjusted sample range.
 * If such values are encountered, either flag an error (if
 * @misst != %NULL), or construct a "missing mask" and attach it
 * to @pmod.
 *
 * Returns: if @misst is not NULL, either the ID number of the
 * variable for which a missing value is first found inside the
 * adjusted sample range or 0 if there is no such variable.  If
 * @misst is NULL, returns 1 if there is an error creating
 * the missing obs mask, otherwise 0.
 */

int model_adjust_sample(Model * pmod, const(DatasetRaw) * dset,
                        int * misst)
{
    const(double) ** Z = cast(const(double) **) dset.Z;
    int i, t, dwt = 0;
    int t1min = pmod.t1;
    int t2max = pmod.t2;
    int vi, missobs;
    int move_ends = 1;
    int ret = 0;

    if (gretl_model_get_int(pmod, "wt_dummy")) {
        /* we have a weight variable which is a 0/1 dummy */
        dwt = pmod.nwt;
    }

    /* advance start of sample range to skip missing obs */
    for (t=t1min; t<t2max; t++) {
        missobs = 0;
        for (i=1; i<=pmod.list[0]; i++) {
            vi = pmod.list[i];
            if (vi > 0 && vi != LISTSEP) {
                if (really_missing(vi, t, Z, dwt)) {
                    missobs = 1;
                    break;
                }
            }
        }
        if (missobs) {
            t1min++;
        } else {
            break;
        }
    }

    /* retard end of sample range to skip missing obs */
    for (t=t2max; t>t1min; t--) {
        missobs = 0;
        for (i=1; i<=pmod.list[0]; i++) {
            vi = pmod.list[i];
            if (vi > 0 && vi != LISTSEP) {
                if (really_missing(vi, t, Z, dwt)) {
                    missobs = 1;
                    break;
                }
            }
        }
        if (missobs) {
            t2max--;
        } else {
            break;
        }
    }

    if (misst != null) {
        /* check for missing values within remaining range and
           flag an error in case any are found */
        for (t=t1min; t<=t2max && !ret; t++) {
            for (i=1; i<=pmod.list[0]; i++) {
                vi = pmod.list[i];
                if (vi > 0 && vi != LISTSEP) {
                    if (really_missing(vi, t, Z, dwt)) {
                        /* identify first missing obs and var */
                        *misst = t + 1;
                        ret = vi;
                        break;
                    }
                }
            }
        }
    } else {
        /* construct a mask for missing values within remaining range?
           we do this only if misst == NULL */
        missobs = 0;
        for (t=t1min; t<=t2max; t++) {
            for (i=1; i<=pmod.list[0]; i++) {
                vi = pmod.list[i];
                if (vi > 0 || vi != LISTSEP) {
                    if (really_missing(vi, t, Z, dwt)) {
                        missobs++;
                        break;
                    }
                }
            }
        }
        if (missobs == t2max - t1min + 1) {
            /* no valid observations */
            pmod.errcode = GretlError.missdata;
            ret = 1;
        } else if (missobs > 0) {
            pmod.missmask = model_missmask(pmod.list, pmod.t1, pmod.t2,
                                           dset.n, Z, dwt, null);
            move_ends = 0;
            if (pmod.missmask == null) {
                pmod.errcode = GretlError.alloc;
                ret = 1;
            }
        }
    }

    if (move_ends) {
        pmod.t1 = t1min;
        pmod.t2 = t2max;
    }

    return ret;
}

/**
 * first_missing_index:
 * @x: array to be checked for missing values.
 * @t1: start of range to check.
 * @t2: end of range to check.
 *
 * Returns: the index of the first missing observation in @x
 * over the sample range @t1 to @t2, or -1 if there is no
 * such observation.
 */

int first_missing_index(const(double) * x, int t1, int t2)
{
    int t;

    for (t=t1; t<=t2; t++) {
        if (na(x[t])) {
            return t;
        }
    }

    return -1;
}

/**
 * series_adjust_sample:
 * @x: series to be checked for missing values.
 * @t1: on entry, initial start of sample range; on exit,
 *      start of sample range adjusted for missing values.
 * @t2: on entry, initial end of sample range; on exit, end
 *      of sample range adjusted for missing values.
 *
 * Adjusts @t1 and @t2 so as to drop any leading or trailing
 * missing observations.
 *
 * Returns: E_MISSDATA if interior missing values were found
 * within the (possibly adjusted) sample range, otherwise 0.
 */

int series_adjust_sample(const(double) * x, int * t1, int * t2)
{
    int t, t1min = *t1, t2max = *t2;
    int err = 0;

    for (t=t1min; t<t2max; t++) {
        if (na(x[t])) t1min++;
        else break;
    }

    for (t=t2max; t>t1min; t--) {
        if (na(x[t])) t2max--;
        else break;
    }

    for (t=t1min; t<=t2max; t++) {
        if (na(x[t])) {
            err = GretlError.missdata;
            break;
        }
    }

    *t1 = t1min;
    *t2 = t2max;

    return err;
}

/**
 * list_adjust_sample:
 * @list: list of variables to be tested for missing values,
 * or %NULL to test all series.
 * @t1: on entry, initial start of sample range; on exit,
 *      start of sample range adjusted for missing values.
 * @t2: on entry, initial end of sample range; on exit, end
 *      of sample range adjusted for missing values.
 * @dset: dataset struct.
 * @nmiss: location to receive number of missing values within
 * (possibly adjusted) sample range.
 *
 * Drops leading or trailing observations from the sample range
 * initially given by the values in @t1 and @t2 if missing values
 * are found for any of the variables given in @list.
 *
 * If @nmiss is non-NULL it receives the number of missing values
 * inside the (possibly reduced) sample range, otherwise it is
 * considered an error if there are any such missing values.
 *
 * Returns: 0 on success or %E_MISSDATA or error.
 */

int list_adjust_sample(const(int) * list, int * t1, int * t2,
                       const(DatasetRaw) * dset, int * nmiss)
{
    int i, t, t1min = *t1, t2max = *t2;
    int k, vi, missing, err = 0;

    if (list != null) {
        k = list[0];
    } else {
        /* check all series */
        k = dset.v - 1;
    }

    /* advance start of sample range to skip missing obs? */
    for (t=t1min; t<t2max; t++) {
        missing = 0;
        for (i=1; i<=k; i++) {
            vi = list == null ? i : list[i];
            if (vi > 0 && vi != LISTSEP) {
                if (na(dset.Z[vi][t])) {
                    missing = 1;
                    break;
                }
            }
        }
        if (missing) {
            t1min++;
        } else {
            break;
        }
    }

    /* retard end of sample range to skip missing obs? */
    for (t=t2max; t>t1min; t--) {
        missing = 0;
        for (i=1; i<=k; i++) {
            vi = list == null ? i : list[i];
            if (vi > 0 && vi != LISTSEP) {
                if (na(dset.Z[vi][t])) {
                    missing = 1;
                    break;
                }
            }
        }
        if (missing) {
            t2max--;
        } else {
            break;
        }
    }

    if (nmiss != null) {
        *nmiss = 0;
    }

    /* check for missing values within remaining range */
    for (t=t1min; t<=t2max && !err; t++) {
        missing = 0;
        for (i=1; i<=k; i++) {
            vi = list == null ? i : list[i];
            if (vi > 0 && vi != LISTSEP) {
                if (na(dset.Z[vi][t])) {
                    if (nmiss == null) {
                        err = GretlError.missdata;
                    } else {
                        *nmiss += 1;
                    }
                    break;
                }
            }
        }
    }

    *t1 = t1min;
    *t2 = t2max;

    return err;
}

/* If there are any missing data within the sample range pmod->t1 to
   pmod->t2, as gauged by reference to pmod->uhat, returns an array of
   char of length pmod->t2 - pmod->t1 + 1 with 1s indicating missing
   observations and 0s otherwise. Otherwise (i.e. there are no missing
   data in the given range) returns NULL.
*/

char * model_sample_mask(const(Model) * pmod)
{
    char * ret = null;
    int doit = 0;
    int t;

    for (t=pmod.t1; t<=pmod.t2; t++) {
        if (na(pmod.uhat[t])) {
            doit = 1;
            break;
        }
    }

    if (doit) {
        int i, n = pmod.t2 - pmod.t1 + 1;

        ret = cast(char *) calloc(n, 1);
        for (t=pmod.t1, i=0; t<=pmod.t2; t++, i++) {
            if (na(pmod.uhat[t])) {
                ret[i] = 1;
            }
        }
    }

    return ret;
}

/* A copy of a given model's missmask, or NULL. If non-NULL this will
   be a NUL-terminated array of char of length equal to the full
   length of the dataset on which the model was estimated, with
   elements '1' for each observation in the sample range t >= pmod->t1
   and t <= pmod->t2 at which the model residual is missing and '0'
   for all other observations.
*/
private __gshared char * refmask;

/* For handling the "omit" command, applied to a model that has
   missing values within the sample range. When the model is
   re-estimated with a reduced set of regressors we must use exactly
   the same sample as the original model, otherwise the comparison of
   the original and reduced models will be invalid.
*/

private char * build_refmask_from_model(const(Model) * pmod)
{
    int t, n = pmod.full_n;
    char * mask = cast(char *) malloc(n + 1);

    if (mask != null) {
        memset(mask, '0', n);
        mask[n] = '\0'; /* NUL-terminate */
        for (t=pmod.t1; t<=pmod.t2; t++) {
            if (na(pmod.uhat[t])) {
                mask[t] = '1';
            }
        }
    }

    return mask;
}

/* Set the "reference" missing observations mask, @refmask,
   based on @pmod.
*/

void set_reference_missmask_from_model(const(Model) * pmod)
{
    if (pmod != null) {
        if (refmask != null) {
            free(refmask);
            refmask = null;
        }
        if (pmod.missmask != null) {
            refmask = gretl_strdup(pmod.missmask);
        } else if (model_has_missing_obs(pmod)) {
            refmask = build_refmask_from_model(pmod);
        }
    }
}

int copy_to_reference_missmask(const(char) * mask)
{
    if (refmask != null) {
        free(refmask);
        refmask = null;
    }

    refmask = gretl_strdup(mask);

    return (refmask == null)? GretlError.alloc : 0;
}

/* attach the reference missing obs mask to a specified model */

int apply_reference_missmask(Model * pmod)
{
    if (refmask != null) {
        pmod.missmask = refmask;
        refmask = null; /* note: apply it only once */
    }

    return 0;
}

int reference_missmask_present()
{
    return (refmask != null);
}

private int real_setmiss(double missval, int varno,
                         DatasetRaw * dset)
{
    int i, t, count = 0;
    int start = 1, end = dset.v;

    if (varno) {
        start = varno;
        end = varno + 1;
    }

    for (i=start; i<end; i++) {
        for (t=0; t<dset.n; t++) {
            if (dset.Z[i][t] == missval) {
                dset.Z[i][t] = NADBL;
                count++;
            }
        }
    }

    return count;
}

/**
 * set_miss:
 * @list: list of variables to process, or an empty list or %NULL
 *        to process all variables.
 * @param: string representation of the numerical value to treat as missing.
 * @dset: dataset struct.
 * @prn: pointer to printing struct.
 *
 * Set to "missing" each observation of each series in @list that
 * has the value represented by @param.
 *
 * Returns: 0 on success, non-zero code on failure.
 */

int set_miss(const(int) * list, const(char) * param,
             DatasetRaw * dset, Prn * prn)
{
    int i, vi, count;
    double missval;
    int err = 0;

    if (param == null) {
        return GretlError.args;
    }

    if (dset == null || dset.n == 0) {
        return GretlError.nodata;
    }

    missval = gretl_double_from_string(param, &err);
    if (err) {
        return err;
    }

    if (list == null || list[0] == 0) {
        count = real_setmiss(missval, 0, dset);
        if (count) {
            pprintf(prn, _("Set %d values to \"missing\"\n").ptr, count);
        } else {
            pputs(prn, _("Didn't find any matching observations\n").ptr);
        }
    } else {
        for (i=1; i<=list[0]; i++) {
            vi = list[i];
            if (vi == 0 || object_is_const(dset.varname[vi], vi)) {
                gretl_errmsg_sprintf(_("The variable %s is read-only").ptr,
                                     dset.varname[vi]);
                err = GretlError.data;
                break;
            }
            count = real_setmiss(missval, vi, dset);
            if (count) {
                pprintf(prn, _("%s: set %d observations to \"missing\"\n").ptr,
                        dset.varname[vi], count);
            } else {
                pprintf(prn, _("%s: Didn't find any matching observations\n").ptr,
                        dset.varname[vi]);
            }
        }
    }

    return err;
}

/**
 * missing_obs_fraction:
 * @dset: dataset struct.
 *
 * Returns: the fraction of the observations in @Z for which
 * all variables have missing values (empty rows).
 */

double missing_obs_fraction(const(DatasetRaw) * dset)
{
    int missrow, totmiss = 0;
    int i, t;

    if (dset.n == 0) {
        return 0.0;
    }

    for (t=0; t<dset.n; t++) {
        missrow = 1;
        for (i=1; i<dset.v; i++) {
            if (!na(dset.Z[i][t])) {
                missrow = 0;
                break;
            }
        }
        totmiss += missrow;
    }

    return cast(double) totmiss / dset.n;
}

/**
 * any_missing_user_values:
 * @dset: dataset struct.
 *
 * Returns: 1 if there are missing values for any non-hidden
 * variables within the current sample range, otherwise 0.
 */

int any_missing_user_values(const(DatasetRaw) * dset)
{
    int i, t;

    if (dset.n == 0) {
        return 0;
    }

    for (i=1; i<dset.v; i++) {
        if (!series_is_hidden(dset, i)) {
            for (t=dset.t1; t<=dset.t2; t++) {
                if (na(dset.Z[i][t])) {
                    return 1;
                }
            }
        }
    }

    return 0;
}

/**
 * model_add_missmask:
 * @pmod: pointer to gretl MODEL.
 * @n: length of mask.
 *
 * Allocates a missing obs mask of length @n and attaches it
 * to @pmod. All elements are initialized to '0' (non-missing).
 *
 * Returns: 0 on success, non-zero if allocation fails.
 */

int model_add_missmask(Model * pmod, int n)
{
    pmod.missmask = cast(char *) malloc((n+1) * (*pmod.missmask).sizeof);
    if (pmod.missmask == null) {
        return GretlError.alloc;
    } else {
        memset(pmod.missmask, '0', n);
        pmod.missmask[n] = 0; /* NUL-termination */
        return 0;
    }
}
