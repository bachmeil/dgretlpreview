/* Port of lib/src/subsample.c: sub-sampling of the dataset (restricting by a boolean condition, dummy, missing
   values, random draws or ranges), panel sample setting, resampling, the stored sub-sample masks of models, and
   the full-dataset bookkeeping. The header subsample.h is subsampleh.d. */

module gretl.subsample;


import core.stdc.ctype : isdigit, isspace;
import core.stdc.errno : errno, ERANGE;
import core.stdc.math : ceil, exp, fabs, floor, fmod, isinf, isnan, log, log10, pow, sqrt;
import core.stdc.stdio : fputs, fprintf, printf, puts, sprintf, snprintf, sscanf, stderr;
import core.stdc.stdlib : abs, atoi, calloc, free, getenv, malloc, qsort, realloc;
import core.stdc.string : memcpy, memset, strcmp, strcpy, strlen, strchr, strstr, strcat, strncat, strncmp, strspn, strcspn, memmove, strrchr;

import std.algorithm : max, min;

import gretl.gretl_errorsh : GretlError;
import gretl.gretl_matrix;
import gretl.gretl_matrixh;
import gretl.johansenh : JohansenCode, JohansenInfo, jcode, jrank, auto_restr;
import gretl.libgretlh;
import gretl.missingh : NADBL, na;
import gretl.gretl_restricth : GretlRestriction;
import gretl.varh;

import gretl.xmlh;


private alias XUC = const(xmlChar) *;

import gretl.calendar : calendar_obs_number, epoch_day_from_ymd_basic, get_epoch_day, weekday_from_date, ymd_extended_from_epoch_day;
import gretl.cmd_privateh : CmdFlags, ExecState;
import gretl.dataio : dateton, get_date_x, ntolabel;
import gretl.dataioh : DataClearCode;
import gretl.dataset : clear_datainfo, copy_varinfo, create_auxiliary_dataset, datainfo_new, dataset_allocate_obs_markers, dataset_clear_sample_record, dataset_destroy_hidden_variables, dataset_shrink_obs_range, destroy_dataset, free_Z, get_resampling_seed, is_string_valued, series_get_label, series_is_hidden, series_recode_strings, shrink_varinfo, start_new_Z;
import gretl.dataseth : DatasetMarkerType, annual_data, calendar_data, custom_time_series, dataset_is_cross_section, dataset_is_panel, dataset_is_time_series, dated_daily_data, sample_size;
import gretl.genmain : current_series_index, generate_scalar, generate_series, series_index;
import gretl.genmainh : ModelDataIndex;
import gretl.glib : GDateWeekday, guint32;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_errors : gretl_errmsg_append, gretl_errmsg_set, gretl_errmsg_sprintf, gretl_error_clear;
import gretl.gretl_list : full_var_list, gretl_list_cmp, gretl_list_new;
import gretl.gretl_model : gretl_is_between_model, gretl_model_get_series, highest_numbered_var_in_model;
import gretl.gretl_panel : obs_index_from_aqm, time_series_from_panel, undo_panel_padding;
import gretl.gretl_prn : Prn, pprintf, pputc, pputs;
import gretl.gretl_utils : simple_set_obs;
import gretl.gretl_xml : gretl_xml_put_tagged_string;
import gretl.libset : gretl_in_gui_mode, gretl_messages_on;
import gretl.missing : list_adjust_sample;
import gretl.objstackh : GretlObjType;
import gretl.options : incompatible_options, print_flags;
import gretl.strutils : gretl_scan_varname, gretl_strdup, integer_string, print_time;
import gretl.subsampleh : RESAMPLED, SUBMASK_SENTINEL;
import gretl.uservar : gretl_is_scalar, gretl_scalar_get_value;
import gretl.genfuncs : fill_permutation_vector, get_t_from_obs_string;

/* unported dependencies: per-file declarations from the 2026c headers */
extern (C) {
    int check_models_for_subsample(char * newmask, int * ndropped);
    int db_set_sample(const(char) * star, const(char) * stop, DatasetRaw * dset);
    void * get_last_model(GretlObjType * type);
    int gretl_function_depth();
    int n_stacked_models();
    void sample_range_get_extrema(const(DatasetRaw) * dset, int * t1, int * t2);
}


enum SubsampleMode {
    none,
    drop_missing,
    drop_empty,
    drop_wkends,
    use_dummy,
    boolean,
    random,
    unknown
}

int recode_strvals(DatasetRaw * dset, Gretlopt opt)
{
    int i, err = 0;

    for (i=1; i<dset.v && !err; i++) {
    if (is_string_valued(dset, i)) {
        err = series_recode_strings(dset, i, opt, null);
    }
    }

    return err;
}

/*
  The purpose of the static pointers below: When the user subsamples
  the current dataset in a non-trivial way -- i.e., by selecting cases
  rather than just moving the starting or ending points of the data
  range -- we create a new_ sub-dataset, and we need to keep the full
  dataset around so that it can be restored later.  The pointer
  @fullset is used to record the address of the full dataset.

  In addition, @peerset keeps track of the location of the DATASET
  struct associated with the backed-up full dataset; by means of this,
  we can know when to free the full dataset and when not to (for
  instance, if we're freeing an auxiliary dataset).
*/

__gshared DatasetRaw * fullset;
__gshared DatasetRaw * peerset;

/* accessors for full dataset, when sub-sampled */

DatasetRaw * fetch_full_dataset()
{
    return fullset;
}

int submask_length(const(char) * s)
{
    int n = 1;

    if (s == null || s == RESAMPLED) {
    n = 0;
    } else {
    while (*s != SUBMASK_SENTINEL) {
        n++;
        s++;
    }
    }

    return n;
}

int get_dataset_submask_size(const(DatasetRaw) * dset)
{
    int n = 0;

    if (dset != null) {
    n = submask_length(dset.submask);
    if (n > 0) n--; /* omit the sentinel */
    }

    return n;
}

int get_model_submask_size(const(Model) * pmod)
{
    int n = 0;

    if (pmod != null) {
    n = submask_length(pmod.submask);
    if (n > 0) n--; /* omit the sentinel */
    }

    return n;
}

void free_subsample_mask(char * s)
{
    if (s != RESAMPLED && s != null) {
    free(s);
    }
}

char * copy_subsample_mask(const(char) * src, int * err)
{
    char * ret = null;

    if (src == RESAMPLED) {
    ret = RESAMPLED;
    } else if (src != null) {
    int n = submask_length(src);

    ret = cast(typeof(ret)) malloc(n * (*ret).sizeof);
    if (ret != null) {
        memcpy(ret, src, n);
    } else {
        *err = GretlError.alloc;
    }
    }

    return ret;
}

char * copy_dataset_submask(const(DatasetRaw) * dset, int * err)
{
    char * mask = null;

    if (complex_subsampled()) {
    mask = copy_subsample_mask(dset.submask, err);
    }

    return mask;
}

int write_model_submask(const(Model) * pmod, Prn * prn)
{
    int ret = 0;

    if (pmod.submask == RESAMPLED) {
    pputs(prn, "<submask length=\"0\"></submask>\n");
    ret = 1;
    } else if (pmod.submask != null) {
    int i, n = submask_length(pmod.submask);

    pprintf(prn, "<submask length=\"%d\">".ptr, n);
    for (i=0; i<n; i++) {
        if (pmod.submask[i] == 'p') {
        pputs(prn, "p ");
        } else {
        pprintf(prn, "%d ".ptr, cast(int) pmod.submask[i]);
        }
    }
    pputs(prn, "</submask>\n");
    ret = 1;
    }

    return ret;
}

int write_dataset_submask(const(DatasetRaw) * dset, Prn * prn)
{
    int ret = 0;

    if (dset.submask == RESAMPLED) {
    uint seed = cast(uint) get_resampling_seed();

    pprintf(prn, "<resample seed=%lu\" n=\"%d\"/>\n".ptr, seed, dset.n);
    ret = 1;
    } else if (complex_subsampled()) {
    int i, n = submask_length(dset.submask);

    pprintf(prn, "<submask length=\"%d\">".ptr, n);
    for (i=0; i<n; i++) {
        if (dset.submask[i] == 'p') {
        pputs(prn, "p ");
        } else {
        pprintf(prn, "%d ".ptr, cast(int) dset.submask[i]);
        }
    }
    pputs(prn, "</submask>\n");

    if (dset.restriction != null) {
        gretl_xml_put_tagged_string("restriction", dset.restriction, prn);
    }

    ret = 1;
    }

    return ret;
}

int submask_cmp(const(char) * m1, const(char) * m2)
{
    if (m1 == null && m2 == null) {
    return 0;
    } else if (m1 == null || m2 == null) {
    return 1;
    }

    if (m1 == RESAMPLED || m2 == RESAMPLED) {
    return m1 != RESAMPLED || m2 != RESAMPLED;
    }

    while (*m1 != SUBMASK_SENTINEL && *m2 != SUBMASK_SENTINEL) {
    if (*m1 != *m2) {
        return 1;
    }
    m1++;
    m2++;
    }

    return 0;
}

/* All values apart from the sentinel are initialized to zero; once
   the mask is actually used, 1s will indicate included observations
   and 0s excluded observations.
*/

char * make_submask(int n)
{
    char * mask = cast(char *) calloc(n + 2, 1);

    if (mask != null) {
    mask[n] = SUBMASK_SENTINEL;
    }

    return mask;
}

void set_dataset_resampled(DatasetRaw * dset, uint seed)
{
    dset.submask = RESAMPLED;
    dset.rseed = seed;
}

int dataset_is_resampled(const(DatasetRaw) * dset)
{
    return (dset != null && dset.submask == RESAMPLED);
}

void maybe_free_full_dataset(DatasetRaw * dset)
{
    if (dset == peerset) {
    if (fullset != null) {
        if (fullset.Z != null) {
        free_Z(fullset);
        }
        clear_datainfo(fullset, DataClearCode.subsample);
        free(fullset);
        fullset = null;
    }
    if (dset.submask != null) {
        free(dset.submask);
        dset.submask = null;
    }
    peerset = null;
    }
}

/* we do this on "restore full sample" */

void relink_to_full_dataset(DatasetRaw * dset)
{

    *dset = *fullset;
    free(fullset);
    fullset = null;
    peerset = null;
}

void sync_dataset_shared_members(const(DatasetRaw) * dset)
{
    if (fullset != null && dset != fullset) {
    fullset.varname = cast(typeof(fullset.varname)) dset.varname;
    fullset.varinfo = cast(typeof(fullset.varinfo)) dset.varinfo;
    }
}

/* Sync malloced elements of @fullset that might have been moved via
   realloc().
*/

void sync_dataset_members(const(DatasetRaw) * dset)
{
    if (fullset.v > dset.v) {
    int i;

    for (i=dset.v; i<fullset.v; i++) {
        free(fullset.Z[i]);
        fullset.Z[i] = null;
    }
    fullset.v = dset.v;
    }

    fullset.varname = cast(typeof(fullset.varname)) dset.varname;
    fullset.varinfo = cast(typeof(fullset.varinfo)) dset.varinfo;
    fullset.descrip = cast(typeof(fullset.descrip)) dset.descrip;
}

/* attach_subsample_to_model:
 * @pmod: model to which subsample should be attached.
 * @dset: pointer to current dataset info.
 *
 * If the dataset is currently subsampled, record the subsample
 * information with the model so that it can be retrieved later.
 *
 * Returns: 0 if the recording is not needed, or on success; non-zero
 * error code failure.
 */

int attach_subsample_to_model(Model * pmod, const(DatasetRaw) * dset)
{
    int err = 0;

    if (fullset != null) {
    /* sync, in case anything has moved */
    sync_dataset_members(dset);

    if (pmod.submask != null) {
        free_subsample_mask(pmod.submask);
    }

    pmod.submask = copy_subsample_mask(dset.submask, &err);
    }

    return err;
}

/* This is called from objstack.c for each saved model, to determine
   whether a proposed permanent shrinkage of the dataset will invalidate
   the model.
*/

int subsample_check_model(Model * pmod, char * mask)
{
    if (submask_cmp(pmod.submask, mask)) {
    return GretlError.data;
    } else {
    return 0;
    }
}

/* Called from objstack.c for each saved model, when the user calls for
   the imposition of a permanent subsample restriction on the
   dataset. We've already checked that the "new" subsample matches that
   on which model was estimated, so all we have to do now is sync by
   removing the model's subsample mask.
*/

int remove_model_subsample_info(Model * pmod)
{
    free_subsample_mask(pmod.submask);
    pmod.submask = null;

    return 0;
}

/* check for at least two observations, all contiguous */

int ts_contig(const(char) * ts, int T)
{
    int contig = 1;
    int stopped = 0;
    int t, n = 0;

    for (t=0; t<T && contig; t++) {
    if (ts[t]) {
        if (stopped) {
        contig = 0;
        } else {
        n++;
        }
    } else if (n > 0) {
        stopped = 1;
    }
    }

    return n >= 2 && contig;
}

/* For forecasting purposes, check to see if a subsample mask represents
   a sub-sampling of a panel dataset in the time dimension only.
*/

int is_panel_time_sample(const(char) * mask, const(DatasetRaw) * dset)
{
    int i, t, s, N, T, ret = 1;
    char * ts = null;

    if (dset == null || !dataset_is_panel(dset)) {
    return 0;
    } else if (submask_length(mask) != dset.n + 1) {
    return 0;
    }

    T = dset.pd;
    N = dset.n / T;
    ts = cast(typeof(ts)) calloc(T, 1);

    s = 0;
    for (i=0; i<N && ret; i++) {
    for (t=0; t<T && ret; t++) {
        if (i == 0) {
        /* first unit sets the pattern */
        if (mask[s]) {
            ts[t] = 1;
        }
        } else if ((mask[s] && !ts[t]) || (ts[t] && !mask[s])) {
        /* subsequent units must match */
        ret = 0;
        }
        s++;
    }
    if (i == 0 && !ts_contig(ts, T)) {
        ret = 0;
    }
    }

    free(ts);

    return ret;
}

/* If series have been added to a resampled dataset, we can't bring
   these back to the "full" dataset, which may have a longer or shorter
   series length, and from which there is no definite mapping by row. So
   we just delete them. In this function we destroy their varnames and
   varinfo structures; the numerical arrays get deleted later.
*/

int resample_sync_dataset(DatasetRaw * dset)
{
    int err = 0;

    if (dset.v > fullset.v) {
    err = shrink_varinfo(dset, fullset.v);
    }

    /* sync */
    fullset.varname = cast(typeof(fullset.varname)) dset.varname;
    fullset.varinfo = cast(typeof(fullset.varinfo)) dset.varinfo;
    fullset.descrip = cast(typeof(fullset.descrip)) dset.descrip;

    return err;
}

/* Apparatus for updating full dataset when restoring full sample
   after sub-sampling.
*/

void update_full_data_values(const(DatasetRaw) * dset)
{
    int i, s, t;

    for (i=1; i<fullset.v && i<dset.v; i++) {
    s = 0;
    for (t=0; t<fullset.n; t++) {
        if (dset.submask[t] == 1) {
        fullset.Z[i][t] = dset.Z[i][s++];
        } else if (dset.submask[t] == 'p') {
        /* skip panel padding (?) */
        s++;
        }
    }
    }
}

int update_case_markers(const(DatasetRaw) * dset)
{
    int err = 0;

    if (fullset.markers == DatasetMarkerType.daily_date_strings) {
    ; /* don't mess with them */
    } else if (dated_daily_data(fullset)) {
    ; /* again, don't mess! */
    } else if (dset.markers && !fullset.markers) {
    dataset_allocate_obs_markers(fullset);
    if (fullset.S == null) {
        err = 1;
    } else {
        int t, subt = 0;

        for (t=0; t<fullset.n; t++) {
        if (dset.submask[t]) {
            strcpy(fullset.S[t], dset.S[subt++]);
        } else {
            sprintf(fullset.S[t], "%d".ptr, t + 1);
        }
        }
    }
    }

    return err;
}

int add_new_vars_to_full(DatasetRaw * dset)
{
    int V1 = dset.v;
    int V0 = fullset.v;
    int N = fullset.n;
    double ** newZ = null;
    int i, t, s;
    int err = 0;

    if (V1 <= V0) {
    return 0;
    }

    if (dset.submask == null) {
    return GretlError.nomerge;
    }

    /* allocate expanded data array */
    newZ = cast(typeof(newZ)) realloc(fullset.Z, V1 * (*fullset.Z).sizeof);

    if (newZ == null) {
    return GretlError.alloc;
    }

    fullset.Z = newZ;

    for (i=V0; i<dset.v && !err; i++) {
    fullset.Z[i] = cast(typeof(fullset.Z[i])) malloc(N * (**fullset.Z).sizeof);
    if (fullset.Z[i] == null) {
        err = GretlError.alloc;
    }
    }

    if (err) {
    return GretlError.alloc;
    }

    for (i=V0; i<dset.v; i++) {
    s = 0;
    for (t=0; t<N; t++) {
        fullset.Z[i][t] = (dset.submask[t])?
        dset.Z[i][s++] : NADBL;
    }
    }

    fullset.v = V1;

    return 0;
}

int sync_data_to_full(DatasetRaw * dset)
{
    int err;

    /* update values for pre-existing series, which may have been
       modified via "genr" etc. */
    update_full_data_values(dset);

    /* if case markers were added when subsampled, carry them back */
    update_case_markers(dset);

    /* delete any newly added hidden vars */
    err = dataset_destroy_hidden_variables(dset, fullset.v);

    /* in case any new vars were added when subsampled, try to merge
       them into the full dataset */
    if (!err) {
    err = add_new_vars_to_full(dset);
    }

    return err;
}

/* Here we make a mask representing the "complete" sample restriction
   currently in force, for use in cumulating restrictions. "Complete"
   means that we take into account both an existing subsampling
   restriction, if any, and possible setting of the t1 and t2 members
   of @dset to exclude certain observation ranges.

   The mask returned, if non-NULL, will be the length of the full
   dataset. It will be NULL (without error) if the current dataset is
   neither subsampled nor range- restricted.
*/

char * make_current_sample_mask(DatasetRaw * dset, int * err)
{
    char * currmask = null;
    int range_set;
    int s, t;

    range_set = (dset.t1 > 0 || dset.t2 < dset.n - 1);

    if (dset.submask == null) {
    /* no pre-existing mask, so not subsampled, but
       we should restrict currmask to observations
       included in the current (t1, t2) range, if
       applicable
    */
    if (range_set) {
        currmask = make_submask(dset.n);
        if (currmask == null) {
        *err = GretlError.alloc;
        } else {
        for (t=dset.t1; t<=dset.t2; t++) {
            currmask[t] = 1;
        }
        }
    }
    } else if (dset.submask == RESAMPLED) {
    gretl_errmsg_set(_("You cannot subsample a resampled dataset").ptr);
    *err = GretlError.data;
    } else {
    /* there's a pre-existing mask: in addition we
       mask out observations outside of the current
       (t1, t2) range, if applicable
    */
    currmask = copy_subsample_mask(dset.submask, err);
    if (!*err && range_set) {
        s = -1;
        for (t=0; t<fullset.n; t++) {
        if (dset.submask[t]) s++;
        if (s < dset.t1 || s > dset.t2) {
            currmask[t] = 0;
        }
        }
    }
    }

    return currmask;
}

/* Deal with the case where sampling has been done simply by
   moving the sample-range endpoints
*/

int restore_full_easy(DatasetRaw * dset, ExecState * state)
{
    int t1min, t2max;

    if (state == null) {
    /* not inside a function */
    t1min = 0;
    t2max = dset.n - 1;
    dataset_clear_sample_record(dset);
    } else {
    /* don't go outside the bounds set on entry to
       a function */
    sample_range_get_extrema(dset, &t1min, &t2max);
    }

    if (dset.t1 != t1min || dset.t2 != t2max) {
    dset.t1 = t1min;
    dset.t2 = t2max;
    }

    return 0;
}

/* restore_full_sample:
 * @dset: dataset struct.
 * @state: structure representing program execution state,
 * or %NULL.
 *
 * Restore the full data range, undoing any sub-sampling that was
 * previously performed (either by shifting the endpoints of the
 * sample range or by selecting cases on some criterion or other).
 * However, if @state is not %NULL, and if it contains a submask,
 * the "full" range is relative to that mask.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int restore_full_sample(DatasetRaw * dset, ExecState * state)
{
    int err = 0;

    if (dset == null) {
    return GretlError.nodata;
    } else if (!complex_subsampled()) {
    return restore_full_easy(dset, state);
    }

    if (dset != peerset) {
    fprintf(stderr, "restore_full_sample: dset is not peerset!\n".ptr);
    return GretlError.data;
    }

    /* Beyond this point we are doing a non-trivial restoration
       of a stored "full" dataset which has previously been
       subsampled, e.g., by some boolean criterion.
    */

    if (dset.submask == RESAMPLED) {
    err = resample_sync_dataset(dset);
    } else {
    if (dset.padmask != null) {
        fprintf(stderr, "restore_full_sample: first undo panel padding\n".ptr);
        err = undo_panel_padding(dset);
    }
    if (!err) {
        sync_dataset_members(dset);
        err = sync_data_to_full(dset);
    }
    }

    if (err == GretlError.nomerge) {
        gretl_errmsg_set(_("Missing sub-sample information; can't merge data\n").ptr);
    }

    if (err) {
    return err;
    }

    /* destroy sub-sampled data array */
    free_Z(dset);
    clear_datainfo(dset, DataClearCode.subsample);

    relink_to_full_dataset(dset);

    if (state != null) {
    /* in this case restoring the "full" sample really means, relative
       to the original state
    */
        int fdepth = gretl_function_depth();
        int masked = state.submask != null;

    if (masked) {
        err = restrict_sample_from_mask(state.submask, dset, Gretlopt.OPT_NONE);
        }
        if (!err && (!masked || fdepth > 0)) {
        int t1min, t2max;

        sample_range_get_extrema(dset, &t1min, &t2max);
        if (dset.t1 < t1min) {
        dset.t1 = t1min;
        }
        if (dset.t2 > t2max) {
        dset.t2 = t2max;
        }
    }
    }

    return err;
}

char * expand_mask(const(char) * tmpmask, const(char) * oldmask, int * err)
{
    int oldlen = submask_length(oldmask) - 1;
    char * expanded = make_submask(oldlen);

    if (expanded == null) {
    *err = GretlError.alloc;
    } else {
    /* map @tmpmask onto the full data range */
    int t, i = 0;

    for (t=0; t<oldlen; t++) {
        if (oldmask[t] == 1 && tmpmask[i] != SUBMASK_SENTINEL) {
                expanded[t] = tmpmask[i++];
        }
    }
    }

    return expanded;
}

int intersect_masks(char ** pmask, const(char) * orig, int * err)
{
    char * mask = *pmask;
    int origlen = submask_length(orig) - 1;
    int masklen = submask_length(mask) - 1;
    char * alt = null;
    int i, sn = 0;

    if (masklen < origlen) {
        alt = expand_mask(mask, orig, err);
        if (*err == 0) {
        /* avoid memory leakage */
        free(mask);
            *pmask = mask = alt;
        }
    }

    for (i=0; i<origlen && !*err; i++) {
    if (mask[i] == 1 && orig[i] == 1) {
        sn++;
    } else {
        mask[i] = 0;
    }
    }

    return sn;
}

int make_weekday_mask(const(DatasetRaw) * dset, char * mask)
{
    if (!dated_daily_data(dset) ||
    dset.markers != 0 || dset.pd < 6) {
    return GretlError.pdwrong;
    } else {
    char[OBSLEN] datestr;
    int t, wd;

    for (t=0; t<dset.n; t++) {
        ntolabel(datestr.ptr, t, dset);
        wd = weekday_from_date(datestr.ptr);
        mask[t] = (wd >= GDateWeekday.monday && wd <= GDateWeekday.friday);
    }

    return 0;
    }
}

/* write into @mask: 0 for observations at which _all_ variables
   (or all variables in @list, if @list is non-NULL) have missing
   values, 1 for all other observations (that have at least one
   non-NA). (Refinement: if no list is given, we ignore the generic
   variables "time" and "index".)
*/

int make_empty_mask(const(int) * list, const(DatasetRaw) * dset, char * mask)
{
    int i, vi, t, vt;

    if (list != null && list[0] > 0) {
    /* check specified list of variables */
    for (t=0; t<dset.n; t++) {
        mask[t] = 0;
        for (i=1; i<=list[0]; i++) {
        vi = list[i];
        if (!na(dset.Z[vi][t])) {
            mask[t] = 1;
            break;
        }
        }
    }
    } else {
    /* check (almost) all variables */
    vt = current_series_index(dset, "time");
    vi = current_series_index(dset, "index");

    for (t=0; t<dset.n; t++) {
        mask[t] = 0;
        for (i=1; i<dset.v; i++) {
        if (!series_is_hidden(dset, i) &&
            i != vt && i != vi &&
            !na(dset.Z[i][t])) {
            mask[t] = 1;
            break;
        }
        }
    }
    }

    return 0;
}

/* write into @mask: 0 for observations at which _any_ variable
   (or any variable in @list, if @list is non-NULL) has a missing
   value, 1 for all other observations.
*/

int make_missing_mask(const(int) * list, const(DatasetRaw) * dset, char * mask)
{
    int i, vi, t;

    if (list != null && list[0] > 0) {
    /* check specified list of variables */
    for (t=0; t<dset.n; t++) {
        mask[t] = 1;
        for (i=1; i<=list[0]; i++) {
        vi = list[i];
        if (na(dset.Z[vi][t])) {
            mask[t] = 0;
            break;
        }
        }
    }
    } else {
    /* check all variables */
    for (t=0; t<dset.n; t++) {
        mask[t] = 1;
        for (i=1; i<dset.v; i++) {
        if (!series_is_hidden(dset, i) && na(dset.Z[i][t])) {
            mask[t] = 0;
            break;
        }
        }
    }
    }

    return 0;
}

int copy_dummy_to_mask(char * mask, const(double) * x, int n)
{
    int t, nsel = 0, err = 0;

    for (t=0; t<n && !err; t++) {
    if (x[t] == 1.0) {
        mask[t] = 1;
        nsel++;
    } else if (!na(x[t]) && x[t] != 0.0) { /* NA? */
        err = 1;
    }
    }

    if (!err && nsel == 0) {
    err = GretlError.zero;
    }

    return err;
}

/* When sub-sampling on some boolean criterion, check to see if we can
   meet the criterion by simply adjusting the endpoints of the sample
   range: life will be simpler if that is so.
*/

int mask_get_t1_t2(const(char) * mask, const(DatasetRaw) * dset, int * pt1, int * pt2)
{
    const(char) * s = mask;
    int contig = 0;
    int nseq = 0; /* number of sequences of 1s */
    int t1 = -1;
    int t2 = -1;
    int t;

    for (t=0; s[t] != SUBMASK_SENTINEL; t++) {
        if (s[t] == 1) {
            if (t == 0 || s[t-1] == 0) {
                nseq++;
            }
            if (t1 < 0) {
                t1 = t;
            }
            if (t > t2) {
                t2 = t;
            }
        }
    }

    /* did we get a single sequence of 1s? */
    contig = (nseq == 1);

    if (contig && dataset_is_panel(dset)) {
    int n = t2 - t1 + 1;

    /* sample must leave a whole number of panel units; moreover,
       to retain "panelness" this number must be greater than 1
    */
    if (t1 % dset.pd != 0 || n % dset.pd != 0) {
        contig = 0;
    } else if (n == dset.pd) {
        contig = 0;
    }
    }

    if (contig) {
        *pt1 = t1;
        *pt2 = t2;
    } else {
        *pt1 = *pt2 = -1;
    }

    return contig;
}

int mask_from_temp_dummy_full(const(char) * s, DatasetRaw * dset, char * mask, Prn * prn, int * new_t1, int * new_t2)
{
    char[MAXLINE] formula;
    double * x;
    int err = 0;

    *formula.ptr = '\0';
    strncat(formula.ptr, s, MAXLINE - 1);

    x = generate_series(formula.ptr, dset, prn, &err);

    if (!err) {
    err = copy_dummy_to_mask(mask, x, dset.n);
    if (err == GretlError.zero) {
        gretl_errmsg_set(_("No observations would be left!").ptr);
    } else if (err) {
        gretl_errmsg_sprintf(_("'%s' is not a dummy variable").ptr, "mask".ptr);
    } else if (new_t1 != null && new_t2 != null) {
            mask_get_t1_t2(mask, dset, new_t1, new_t2);
        }
    }

    free(x);

    return err;
}

int mask_from_temp_dummy(const(char) * s, DatasetRaw * dset, char * mask, Prn * prn)
{
    return mask_from_temp_dummy_full(s, dset, mask, prn, null, null);
}

int mask_from_dummy(const(char) * s, const(DatasetRaw) * dset, char * mask)
{
    double * x = null;
    char[VNAMELEN] dname = 0;
    int free_x = 0;
    int err = 0;

    if (!strcmp(s, "$sample")) {
    /* try fetching the sample series from the last model */
    GretlObjType type = cast(GretlObjType) 0;
    void * ptr = get_last_model(&type);

    if (ptr == null || type != GretlObjType.eqn) {
        err = GretlError.data;
    } else {
        strcpy(dname.ptr, s);
        x = cast(typeof(x)) malloc(dset.n * (*x).sizeof);
        if (x == null) {
        err = GretlError.alloc;
        } else {
        free_x = 1;
        err = gretl_model_get_series(x, cast(Model *) ptr, dset, ModelDataIndex.sample);
        }
    }
    } else {
    int dnum;

    gretl_scan_varname(s, dname.ptr);
    dnum = series_index(dset, dname.ptr);
    if (dnum == dset.v) {
        gretl_errmsg_sprintf(_("Variable '%s' not defined").ptr, dname.ptr);
        err = GretlError.data;
    } else {
        x = cast(double *) dset.Z[dnum];
    }
    }

    if (!err) {
    err = copy_dummy_to_mask(mask, x, dset.n);
    if (err) {
        gretl_errmsg_sprintf(_("'%s' is not a dummy variable").ptr, dname.ptr);
    }
    }

    if (free_x) {
    free(x);
    }

    return err;
}

/* How many observations are selected by @mask relative to @dset?
   Observations where mask == 1 are "genuine" observations; those where
   mask == 'p' are rows of NAs inserted to preserve panel structure.

   As of 2026-08-26 the 'p' cases were not being counted here, leading
   to a crash when opening a gretl session in which a panel dataset was
   subsampled with padding. Now these cases are counted and that crash
   is fixed, but we need to be on the lookout for possible bad side
   effects of this change.
*/

int count_selected_cases(const(char) * mask, const(DatasetRaw) * dset)
{
    int i, n = 0;

    for (i=0; i<dset.n; i++) {
    if (mask[i] == 1 || mask[i] == 'p') {
        n++;
    }
    }

    return n;
}

/* panel data: how many distinct cross-sectional units are included in
   the masked subset of observations? */

int count_panel_units(const(char) * mask, const(DatasetRaw) * dset)
{
    int u, ubak = -1;
    int i, n = 0;

    for (i=0; i<dset.n; i++) {
    if (mask[i]) {
        u = i / dset.pd;
        if (u != ubak) {
        n++;
        ubak = u;
        }
    }
    }

    return n;
}

/* Construct a mask for taking a random sub-sample from @dset: we're
   selecting @subn cases without replacement, and there may or may not
   be an existing mask in place that has to be respected.
*/

int make_random_mask(const(char) * s, const(char) * oldmask, DatasetRaw * dset, char * mask)
{
    GretlMatrixRaw * v = null;
    int oldn, subn;
    int i, j;
    int err = 0;

    /* how many cases are requested? */
    subn = smpl_get_int(s, dset, 0, &err);

    if (subn <= 0 || subn >= dset.n) {
    err = GretlError.invarg;
    } else if (oldmask != null) {
    /* dataset is already sub-sampled */
    oldn = 0;
    for (i=0; i<dset.n; i++) {
        if (oldmask[i]) {
        oldn++;
        }
    }
    if (subn >= oldn) {
        err = GretlError.invarg;
    }
    } else {
    oldn = dset.n;
    }

    if (err && err != GretlError.parse) {
    gretl_errmsg_sprintf(_("Invalid number of cases %d").ptr, subn);
    return err;
    }

    v = gretl_zero_matrix_new(subn, 1);
    err = fill_permutation_vector(v, oldn);
    if (err) {
    gretl_matrix_free(v);
    return err;
    }

    /* At this point @v will contain @subn indices into the sub-array of
       1s in @oldmask, or into the observations in @dset if @oldmask is
       NULL. We need to use these to find the indices of the
       observations that should be retained by @mask. Sorting of @v
       might be worthwhile when @oldmask is non-NULL.
    */

    if (oldmask == null) {
    for (j=0; j<subn; j++) {
        i = cast(int) (v.val[j] - 1);
        mask[i] = 1;
    }
    } else {
    int k = 0;

    for (i=0; i<dset.n; i++) {
        if (oldmask[i] == 1) {
        k++;
        for (j=0; j<subn; j++) {
            if (v.val[j] == cast(double) k) {
            mask[i] = 1;
            break;
            }
        }
        }
    }
    }

    gretl_matrix_free(v);

    return err;
}

int backup_full_dataset(DatasetRaw * dset)
{

    if (fullset == null) {
    fullset = cast(typeof(fullset)) malloc((*fullset).sizeof);
    if (fullset == null) {
        return GretlError.alloc;
    }
    }

    if (dset != null) {
    *fullset = *dset;
    peerset = dset;
    }

    return 0;
}

void destroy_full_dataset(DatasetRaw * dset)
{
    if (fullset != null) {
    if (fullset.varname == dset.varname) {
        fullset.varname = null;
    }
    if (fullset.varinfo == dset.varinfo) {
        fullset.varinfo = null;
    }
    if (fullset.descrip == dset.descrip) {
        fullset.descrip = null;
    }
    destroy_dataset(fullset);
    fullset = null;
    }

    dset.varname = null;
    dset.varinfo = null;
    dset.descrip = null;

    free_Z(dset);
    clear_datainfo(dset, DataClearCode.subsample);

    peerset = null;
}

int complex_subsampled()
{
    return (fullset != null && fullset.Z != null);
}

int get_full_length_n()
{
    return (fullset != null) ? fullset.n : 0;
}

int dataset_is_subsampled(const(DatasetRaw) * dset)
{
    if (fullset != null && fullset.Z != null) {
    return 1;
    } else if (dset.t1 > 0 || dset.t2 < dset.n - 1) {
    return 1;
    } else {
    return 0;
    }
}

void copy_data_to_subsample(DatasetRaw * subset, const(DatasetRaw) * dset, int maxv, const(char) * mask)
{
    int i, t, s;

    /* copy data values */
    for (i=1; i<maxv; i++) {
    s = 0;
    for (t=0; t<dset.n; t++) {
        if (mask == null) {
        subset.Z[i][s++] = dset.Z[i][t];
        } else if (mask[t] == 1) {
        subset.Z[i][s++] = dset.Z[i][t];
        } else if (mask[t] == 'p') {
        /* panel padding */
        subset.Z[i][s++] = NADBL;
        }
    }
    }

    /* copy observation markers, if any */
    if (dset.markers && subset.markers) {
    s = 0;
    for (t=0; t<dset.n; t++) {
        if (mask == null || mask[t] == 1 || mask[t] == 'p') {
        strcpy(subset.S[s++], dset.S[t]);
        }
    }
    }

    /* copy panel time info? */
    if (dataset_is_panel(subset)) {
    if (subset.pd == dset.pd) {
        subset.panel_pd = dset.panel_pd;
        subset.panel_sd0 = dset.panel_sd0;
    }
    }

    if (subset.stobs[0] == '\0' || subset.endobs[0] == '\0') {
    /* impose simple default if not already handled */
    strcpy(subset.stobs.ptr, "1");
    sprintf(subset.endobs.ptr, "%d".ptr, subset.n);
    }
}

int get_restriction_mode(Gretlopt opt, const(char) * panmask)
{
    int mode = SubsampleMode.unknown;

    if (panmask != null) {
    mode = SubsampleMode.boolean;
    } else if (opt & (Gretlopt.OPT_M | Gretlopt.OPT_C)) {
    mode = SubsampleMode.drop_missing;
    } else if (opt & Gretlopt.OPT_R) {
    mode = SubsampleMode.boolean;
    } else if (opt & Gretlopt.OPT_N) {
    mode = SubsampleMode.random;
    } else if (opt & Gretlopt.OPT_O) {
    mode = SubsampleMode.use_dummy;
    } else if (opt & Gretlopt.OPT_A) {
    mode = SubsampleMode.drop_empty;
    } else if (opt & Gretlopt.OPT_W) {
    mode = SubsampleMode.drop_wkends;
    }

    return mode;
}

/* Copy list of panel periods from Ti to T0; return 1 if members of
   the list Ti are consecutive, 0 if they are not.
*/

int copy_periods_list(int * T0, const(int) * Ti)
{
    int j, ret = 1;

    for (j=0; j<=Ti[0] && ret; j++) {
    if (j > 1 && Ti[j] != Ti[j-1] + 1) {
        ret = 0;
    } else {
        T0[j] = Ti[j];
    }
    }

    return ret;
}

/* When sub-sampling panel data on some boolean criterion: see if the
   exclusion of certain rows leaves a balanced panel.  Note that the
   requirement is not simply that each unit has the same number of
   temporal observations -- they must have the _same_ observations,
   and in addition the observations must be temporally contiguous.
*/

int mask_leaves_balanced_panel(char * mask, const(DatasetRaw) * dset)
{
    int T = dset.pd;
    int * T0, Ti;
    int i, u, ubak = -1;
    int ret = 1;

    T0 = gretl_list_new(T);
    Ti = gretl_list_new(T);

    if (T0 == null || Ti == null) {
    free(T0);
    free(Ti);
    return 0;
    }

    T0[0] = Ti[0] = 0;

    for (i=0; i<dset.n && ret; i++) {
    if (mask[i]) {
        u = i / T;
        if (u != ubak) {
        if (Ti[0] > 0) {
            if (T0[0] == 0) {
            /* we haven't made the reference list, T0, yet */
            ret = copy_periods_list(T0, Ti);
            } else if (gretl_list_cmp(Ti, T0)) {
            /* the current list differs from the reference one */
            ret = 0;
            }
        }
        Ti[0] = 1;
        Ti[1] = i % T;
        ubak = u;
        } else {
        Ti[0] += 1;
        Ti[Ti[0]] = i % T;
        }
    }
    }

    if (gretl_list_cmp(Ti, T0)) {
    /* check the last unit */
    ret = 0;
    }

    free(T0);
    free(Ti);

    return ret;
}

int make_panel_submask(char * mask, const(DatasetRaw) * dset, int * err)
{
    int T = dset.pd;
    int N = dset.n / T;
    char * umask, pmask;
    int i, np = 0;

    umask = cast(typeof(umask)) calloc(N + T, 1);
    if (umask == null) {
    *err = GretlError.alloc;
    return 0;
    }

    pmask = umask + N;

    for (i=0; i<dset.n; i++) {
    if (mask[i]) {
        umask[i / T] = 1;
        pmask[i % T] = 1;
    }
    }

    for (i=0; i<dset.n; i++) {
    if (!mask[i]) {
        if (umask[i / T] && pmask[i % T]) {
        mask[i] = 'p'; /* mark as padding row */
        np++;
        }
    }
    }

    free(umask);

    return np;
}

int add_daily_date_strings(char * selected, const(DatasetRaw) * dset, DatasetRaw * subset)
{
    int err, t, i = 0;

    err = dataset_allocate_obs_markers(subset);

    if (!err) {
    subset.markers = DatasetMarkerType.daily_date_strings;
    for (t=0; t<dset.n; t++) {
        if (selected[t]) {
        ntolabel(subset.S[i++], t, dset);
        }
    }
    }

    return err;
}

/* See if sub-sampling dated daily data leaves a useable daily dataset
   (as in dropping weekends and/or holidays).
*/

int try_for_daily_subset(char * selected, const(DatasetRaw) * dset, DatasetRaw * subset, Gretlopt * optp)
{
    char[OBSLEN] datestr;
    guint32 ed, ed0 = 0, edbak = 0;
    int t, wd, ngaps = 0;
    int delta, mon_delta;
    int delta_max = 0;
    int started = 0;
    int any_sat = 0;
    int newpd = 5;
    int err = 0;

    if (dset.pd > 5) {
    /* first pass: look at weekend status of included obs */
    for (t=0; t<dset.n; t++) {
        if (selected[t]) {
        ntolabel(datestr.ptr, t, dset);
        wd = weekday_from_date(datestr.ptr);
        if (wd == GDateWeekday.sunday) {
            /* result must be 7-day */
            newpd = 7;
            break;
        } else if (wd == GDateWeekday.saturday) {
            any_sat = 1;
        }
        }
    }
    if (newpd == 5 && any_sat) {
        /* no Sundays, but got a Saurday: result
           must be 6-day data */
        newpd = 6;
    }
    }

    /* determine expected "delta days" for Mondays */
    mon_delta = (newpd == 7) ? 1 : (newpd == 6)? 2 : 3;

    /* second pass: count the calendar gaps */
    for (t=0; t<dset.n; t++) {
    if (selected[t]) {
        ntolabel(datestr.ptr, t, dset);
        wd = weekday_from_date(datestr.ptr);
        ed = get_epoch_day(datestr.ptr);
        if (started) {
        delta = ed - edbak;
        ngaps += (wd == GDateWeekday.monday)? (delta > mon_delta) : (delta > 1);
        if (delta > delta_max) {
            delta_max = delta;
        }
        } else {
        ed0 = ed;
        started = 1;
        }
        edbak = ed;
    }
    }

    /* Now, we should probably restrict the proportion of missing days
       for this treatment: for exclusion of non-trading days 10
       percent should be generous.  But we'll also require that the
       maximum "daily delta" be less than 10.
    */

    if (delta_max < 10) {
    double keepfrac = subset.n / cast(double) dset.n;
    double keepmin = 0.9 * newpd / cast(double) dset.pd;

    if (keepfrac >= keepmin) {
        if (dset.S == null && ngaps > 0) {
        err = add_daily_date_strings(selected, dset, subset);
        /* flag not to destroy subset date strings */
        *optp |= Gretlopt.OPT_P;
        }
        if (!err) {
        subset.structure = TsCodes.time_series;
        subset.pd = newpd;
        subset.sd0 = cast(double) ed0;
        subset.t2 = subset.n - 1;
        ntolabel(subset.stobs.ptr, 0, subset);
        ntolabel(subset.endobs.ptr, subset.t2, subset);
        }
    }
    }

    return err;
}

/* does the given policy lead to an empty sample, or no change
   in the sample, perchance?
*/

int check_subsample_n(int n, DatasetRaw * dset, char ** pmask, Prn * prn)
{
    if (n == 0) {
    gretl_errmsg_set(_("No observations would be left!").ptr);
    return GretlError.data;
    } else if (n == dset.n) {
    /* not really an error, just a no-op */
    if (gretl_messages_on()) {
        pputs(prn, _("No observations were dropped!").ptr);
        pputc(prn, '\n');
    }
    free(*pmask);
    *pmask = null;
    }

    return 0;
}

bool needs_string_arg(int m)
{
    return (m == SubsampleMode.random || m == SubsampleMode.use_dummy || m == SubsampleMode.boolean);
}

int make_restriction_mask(int mode, const(char) * s, const(int) * list, DatasetRaw * dset, const(char) * oldmask, int replace, char ** pmask, int * pt1, int * pt2, Prn * prn)
{
    DatasetRaw * refset = dset;
    char * mask = null;
    int sn = 0, err = 0;

    if (needs_string_arg(mode) && (s == null || *s == '\0')) {
    return GretlError.args;
    }

    if (oldmask != null && mode == SubsampleMode.random) {
    mask = make_submask(fullset.n);
    } else {
    mask = make_submask(dset.n);
    }
    if (mask == null) {
    return GretlError.alloc;
    }

    if (fullset != null && oldmask == null && replace) {
        refset = fullset;
    }

    /* construct subsample mask in one of several possible ways */

    if (mode == SubsampleMode.drop_missing) {
    err = make_missing_mask(list, refset, mask);
    } else if (mode == SubsampleMode.drop_empty) {
    err = make_empty_mask(list, refset, mask);
    } else if (mode == SubsampleMode.drop_wkends) {
    err = make_weekday_mask(refset, mask);
    } else if (mode == SubsampleMode.random) {
    err = make_random_mask(s, oldmask, refset, mask);
    } else if (mode == SubsampleMode.use_dummy) {
    err = mask_from_dummy(s, refset, mask);
    } else if (mode == SubsampleMode.boolean) {
    err = mask_from_temp_dummy(s, refset, mask, prn);
    } else {
    gretl_errmsg_set(_("Sub-sample command failed mysteriously").ptr);
    err = 1;
    }

    if (err) {
    /* exit now on unrecoverable error */
    free(mask);
    return err;
    } else if (mode != SubsampleMode.random) {
        /* try for contiguous sample */
        if (mask_get_t1_t2(mask, dset, pt1, pt2)) {
            *pmask = mask;
            return 0;
        }
    }

    /* cumulate sample restrictions, if appropriate */
    if (oldmask != null && mode != SubsampleMode.random) {
    sn = intersect_masks(&mask, oldmask, &err);
    } else {
    sn = count_selected_cases(mask, refset);
    }

    if (!err) {
        err = check_subsample_n(sn, refset, &mask, prn);
    }

    if (err) {
    free(mask);
    } else {
    *pmask = mask;
    }

    return err;
}

void panel_time_to_time(DatasetRaw * subset, DatasetRaw * dset)
{
    DatasetRaw tset;

    tset.structure = TsCodes.time_series;
    tset.pd = dset.panel_pd;
    tset.sd0 = dset.panel_sd0;
    ntolabel(tset.stobs.ptr, 0, &tset);
    simple_set_obs(subset, dset.panel_pd, tset.stobs.ptr, Gretlopt.OPT_T);
}

void finalize_panel_subset(DatasetRaw * subset, DatasetRaw * dset, int npad)
{
    int pdp = subset.pd;
    int den = 10;

    while ((pdp = pdp / 10) != 0) {
    den *= 10;
    }

    subset.sd0 = 1.0 + 1.0 / den;
    ntolabel(subset.stobs.ptr, 0, subset);
    ntolabel(subset.endobs.ptr, subset.n - 1, subset);

    if (dset.pangrps != null && npad == 0) {
    /* carry over panel group names from full dataset */
    subset.pangrps = gretl_strdup(dset.pangrps);
    }
}

/* This is also used elsewhere:

   in gui2/session.c, when re-establishing a previously sub-sampled data
   state on reopening a saved session

   in gretl_func.c, on exit from a user function when the dataset was
   sub-sampled on entry to the function (and we need to re-establish the
   original sub-sample)
*/

int restrict_sample_from_mask(char * mask, DatasetRaw * dset, Gretlopt opt)
{
    DatasetRaw * subset;
    Gretlopt zopt = Gretlopt.OPT_R;
    int rebalance = (opt & Gretlopt.OPT_B);
    int err = 0;

    if (dset.auxiliary) {
    fprintf(stderr, "restrict_sample_from_mask: attempting to restrict an auxiliary data\n".ptr);
    return GretlError.data;
    }

    if (mask == RESAMPLED) {
    fprintf(stderr, "restrict_sample_from_mask: got RESAMPLED!\n".ptr);
    return GretlError.data;
    }

    if (rebalance && !dataset_is_panel(dset)) {
    gretl_errmsg_set(_("--balanced option is invalid: the (full) dataset is not a panel").ptr);
    return GretlError.badopt;
    }

    subset = datainfo_new();
    if (subset == null) {
    return GretlError.alloc;
    }

    subset.n = count_selected_cases(mask, dset);
    subset.v = dset.v;

    if (dataset_is_panel(dset)) {
    /* are we able to reconstitute a panel? */
    int nunits = count_panel_units(mask, dset);
    int ok = 0, npad = 0;

    if (nunits > 1 && subset.n > nunits) {
        /* there's some possibility of doing so */
        if (rebalance) {
        /* add padding rows? only if this was requested */
        npad = make_panel_submask(mask, dset, &err);
        if (err) {
            free(subset);
            return err;
        }
        ok = 1;
        } else {
        ok = mask_leaves_balanced_panel(mask, dset);
        }
        if (ok) {
        subset.structure = TsCodes.stacked_time_series;
        subset.n += npad;
        subset.pd = subset.n / nunits;
        finalize_panel_subset(subset, dset, npad);
        }
    } else if (nunits == 1 && subset.n == dset.pd) {
        /* time series for single panel unit */
        if (dset.panel_pd > 0) {
        panel_time_to_time(subset, dset);
        } else {
        subset.structure = TsCodes.special_time_series;
        }
    }
    } else if (dated_daily_data(dset)) {
    /* see if we can preserve daily time series */
    err = try_for_daily_subset(mask, dset, subset, &zopt);
    }

    if (!err) {
    /* set up the sub-sampled dataset */
    err = start_new_Z(subset, zopt);
    }

    if (err) {
    free(subset);
    return err;
    }

    /* link (don't copy) varnames and descriptions, since these are
       not dependent on the series length
    */
    subset.varname = dset.varname;
    subset.varinfo = dset.varinfo;
    subset.descrip = dset.descrip;
    subset.mapfile = dset.mapfile;

    /* set up case markers? */
    if (dset.markers) {
    err = dataset_allocate_obs_markers(subset);
    if (err) {
        free_Z(subset);
        free(subset);
        return GretlError.alloc;
    }
    }

    /* copy across data (and case markers, if any) */
    copy_data_to_subsample(subset, dset, dset.v, mask);

    if (opt & Gretlopt.OPT_T) {
    /* --permanent */
    check_models_for_subsample(mask, null);
    destroy_full_dataset(dset);
    recode_strvals(subset, Gretlopt.OPT_NONE);
    } else {
    err = backup_full_dataset(dset);
    subset.submask = copy_subsample_mask(mask, &err);
    }

    /* switch pointers */
    *dset = *subset;
    free(subset);

    return err;
}

/* Below: we do this "precompute" thing if the dataset is already
   subsampled and the user wants to compound the restriction with a
   boolean restriction. One reason for this is that any "obs"
   references in the new_ restriction may get out of whack if we
   restore the full dataset first.  For example, say the spec is
   "obs!=50" to exclude observation 50: presumably the user means to
   exclude the 50th observation in the current, subsampled dataset,
   which may not be the same as the 50th observation in the full
   dataset.
*/

char * precompute_mask(const(char) * s, const(char) * panmask, const(char) * oldmask, DatasetRaw * dset, int * new_t1, int * new_t2, Prn * prn, int * err)
{
    char * tmpmask = null;
    char * newmask = null;

    if (panmask == null) {
    tmpmask = make_submask(dset.n);
    if (tmpmask == null) {
        *err = GretlError.alloc;
    } else {
        /* fill out mask relative to the current, restricted dataset */
        *err = mask_from_temp_dummy_full(s, dset, tmpmask, prn,
                                             new_t1, new_t2);
    }
    if (!*err) {
        if (fullset != null && *new_t1 < 0) {
        newmask = expand_mask(tmpmask, oldmask, err);
        } else {
        newmask = tmpmask;
        tmpmask = null;
        }
    }
    } else if (fullset != null) {
    newmask = expand_mask(panmask, oldmask, err);
    } else {
    newmask = copy_subsample_mask(panmask, err);
    }

    free(tmpmask);

    return newmask;
}

/* Intended for time series data: trim any missing values at the start
   and end of the current sample range, then check the remaining range
   for missing values and flag an error if any are found. If opt
   contains OPT_T (for permanence) and no error is encountered, we
   shrink the dataset to the contiguous range.
*/

int set_contiguous_sample(const(int) * list, DatasetRaw * dset, Gretlopt opt)
{
    int ct1 = dset.t1;
    int ct2 = dset.t2;
    int err = 0;

    if (list != null && list[0] > 0) {
    err = list_adjust_sample(list, &ct1, &ct2, dset, null);
    } else {
    int * biglist = null;
    int nvars = 0;

    biglist = full_var_list(dset, &nvars);
    if (nvars == 0) {
        ; /* no-op */
    } else if (biglist == null) {
        err = GretlError.alloc;
    } else {
        err = list_adjust_sample(biglist, &ct1, &ct2,
                     dset, null);
        free(biglist);
    }
    }

    if (!err) {
    int changed = ct1 != dset.t1 || ct2 != dset.t2;

    dset.t1 = ct1;
    dset.t2 = ct2;
    if (changed && (opt & Gretlopt.OPT_T)) {
        err = perma_sample(dset, Gretlopt.OPT_T, null, null);
    }
    }

    return err;
}

/* Make a string representing the sample restriction. This is for
   reporting purposes. Note that in some cases the incoming @restr
   string may be NULL, for instance if the no-missing option is chosen
   via the gretl GUI.
*/

int make_restriction_string(DatasetRaw * dset, char * old, const(char) * restr, int mode)
{
    char * s = null;
    int n = 0, err = 0;

    dataset_clear_sample_record(dset);

    if (old != null) {
    n = cast(int) (strlen(old));
    }

    if (mode == SubsampleMode.random) {
    n += strlen("random");
    } else if (mode == SubsampleMode.drop_missing) {
    n += strlen("no-missing");
    } else if (restr != null) {
    n += strlen(restr);
    }

    if (n > 0) {
    s = cast(typeof(s)) malloc(n + 5);
    if (s == null) {
        err = GretlError.alloc;
    }
    }

    if (s != null) {
    *s = '\0';
    if (old != null) {
        strcpy(s, old);
        strcat(s, " && ");
    }
    if (mode == SubsampleMode.random) {
        strcat(s, "random");
    } else if (mode == SubsampleMode.drop_missing) {
        strcat(s, "no-missing");
    } else if (restr != null) {
        strcat(s, restr);
    }
    dset.restriction = s;
    }

    return err;
}

int full_sample(const(DatasetRaw) * dset)
{
    if (complex_subsampled()) {
    return 0;
    } else if (dset.t1 > 0 || dset.t2 < dset.n - 1) {
    return 0;
    } else {
    return 1;
    }
}

int handle_ts_perma_restrict(char * mask, DatasetRaw * dset, Gretlopt opt, int t1)
{
    char[OBSLEN] stobs;
    int pd = dset.pd;
    double sd0;
    int err;

    ntolabel(stobs.ptr, t1, dset);
    sd0 = get_date_x(dset.pd, stobs.ptr);

    err = restrict_sample_from_mask(mask, dset, opt);

    if (!err) {
    /* re-establish time-series characteristics */
    dset.structure = TsCodes.time_series;
    dset.pd = pd;
    dset.sd0 = sd0;
    strcpy(dset.stobs.ptr, stobs.ptr);
    ntolabel(dset.endobs.ptr, dset.n - 1, dset);
    }

    return err;
}

int check_restrict_options(Gretlopt opt, const(char) * param, int * replace)
{
    /* We'll accept the redundant combination of the options --dummy
       (OPT_O) and --restrict (OPT_R), but other than that the options
       --dummy,
       --restrict,
       --no-missing (OPT_M),
       --no-all-missing (OPT_A) and
       --random (OPT_N)
       are all mutually incompatible.
    */
    if (incompatible_options(opt, Gretlopt.OPT_O | Gretlopt.OPT_M | Gretlopt.OPT_A | Gretlopt.OPT_N) ||
    incompatible_options(opt, Gretlopt.OPT_R | Gretlopt.OPT_M | Gretlopt.OPT_A | Gretlopt.OPT_N)) {
    return GretlError.badopt;
    }

    /* The --contiguous option (OPT_C) is compatible only with
       --no-missing (which is implied) and --permanent */
    if (opt & Gretlopt.OPT_C) {
    if ((opt & ~(Gretlopt.OPT_C | Gretlopt.OPT_M | Gretlopt.OPT_T)) != Gretlopt.OPT_NONE) {
        return GretlError.badopt;
    }
    }

    if (opt & (Gretlopt.OPT_O | Gretlopt.OPT_R | Gretlopt.OPT_N)) {
    /* parameter is required */
    if (param == null || *param == '\0') {
        return GretlError.args;
    }
    }

    if (opt & Gretlopt.OPT_P) {
    *replace = 1;
    }

    return 0;
}

int check_permanent_option(Gretlopt opt, DatasetRaw * dset, int * n_models)
{
    if (gretl_function_depth() > 0) {
    /* can't permanently shrink the dataset within a function */
    gretl_errmsg_set(_("The dataset cannot be modified at present").ptr);
    return GretlError.data;
    } else if (!(opt & (Gretlopt.OPT_O | Gretlopt.OPT_M | Gretlopt.OPT_A | Gretlopt.OPT_N | Gretlopt.OPT_R | Gretlopt.OPT_C))) {
    /* we need some kind of restriction flag */
    return GretlError.args;
    } else if (gretl_in_gui_mode() && !(opt & Gretlopt.OPT_F)) {
    /* GUI program, without internal "force" option */
    *n_models = n_stacked_models();
    if (*n_models > 0 && !full_sample(dset)) {
        /* too difficult to recover gracefully on error */
        gretl_errmsg_set(_("The full dataset must be restored before imposing a permanent sample restriction").ptr);
        return GretlError.data;
    }
    }

    return 0;
}

/* "Precomputing" a mask means computing a mask based on a currently
   subsampled dataset (before restoring the full dataset, which is a
   part of the subsampling process). We do this if we're cumulating a
   boolean restriction on top of an existing restriction.

   2023-07-20: the comment used to continue: "but not (because serious
   complications can arise) if the current dataset is a panel".
   It now seems that was over-cautious? Old behavior can be restored
   by defining PAN_DONT_PRECOMP to 1.
*/

int do_precompute(int mode, const(char) * panmask, const(char) * oldmask, DatasetRaw * dset)
{
    if (oldmask == null || mode != SubsampleMode.boolean) {
    /* "precomputing" is not required */
    return 0;
    }
    if (panmask != null && fullset != null) {
    /* @panmask will need expanding to the length of the
       full dataset, which is handled by "precomputing"
    */
    return 1;
    }

    return 1;
}

int handle_contiguous_sample(DatasetRaw * dset, int t1, int t2, char * mask, int permanent, Gretlopt opt)
{
    int err = 0;

    if (opt & (Gretlopt.OPT_M | Gretlopt.OPT_A)) {
        /* sampling to skip missing values: this works relative
           to the current sample range
        */
        if (t1 < dset.t1) {
            t1 = dset.t1;
        }
        if (t2 > dset.t2) {
            t2 = dset.t2;
        }
    }

    if (permanent && dataset_is_time_series(dset)) {
        /* apply the restriction, but then re-establish the
           time-series character of the dataset
        */
        err = handle_ts_perma_restrict(mask, dset, opt, t1);
    } else if (!permanent) {
        /* just move the sample range pointers, avoiding
           the overhead of creating a parallel dataset
        */
        dset.t1 = t1;
        dset.t2 = t2;
        /* but maybe record mask? */
    } else {
        err = restrict_sample_from_mask(mask, dset, opt);
    }

    return err;
}

/* Internal version of restrict_sample(), with the extra argument
   @panmask, to represent a time-based sample specification for panel
   data. The latter comes from the "upstream" function
   panel_time_sample().
*/

int real_restrict_sample(const(char) * param, const(int) * list, const(char) * panmask, DatasetRaw * dset, ExecState * state, Gretlopt opt, Prn * prn, int * n_dropped)
{
    char * oldrestr = null;
    char * oldmask = null;
    char * mask = null;
    int free_oldmask = 0;
    int permanent = 0;
    int n_models = 0;
    int n_orig = 0;
    int replace = 0;
    int new_t1 = -1;
    int new_t2 = -1;
    int contig = 0;
    int mode, err = 0;

    if (dset == null || dset.Z == null) {
    return GretlError.nodata;
    }

    n_orig = dset.n;

    /* general check on incoming options */
    err = check_restrict_options(opt, param, &replace);

    /* take a closer look at the --permanent option */
    if (!err && (opt & Gretlopt.OPT_T)) {
    err = check_permanent_option(opt, dset, &n_models);
    if (!err) {
        permanent = 1;
    }
    }

    if (err) {
    return err;
    }

    gretl_error_clear();

    if (opt & Gretlopt.OPT_C) {
    /* --contiguous */
    return set_contiguous_sample(list, dset, opt);
    }

    mode = get_restriction_mode(opt, panmask);
    if (mode == SubsampleMode.unknown) {
    gretl_errmsg_set("Unrecognized sample command");
    return 1;
    }

    if (!replace) {
    /* not replacing but cumulating any existing restrictions */
    oldmask = make_current_sample_mask(dset, &err);
    if (!err) {
        if (oldmask != null && oldmask != RESAMPLED) {
        free_oldmask = 1;
        }
        if (dset.restriction != null) {
        oldrestr = gretl_strdup(dset.restriction);
        }
    }
    } else if (state != null && state.submask != null) {
    /* subsampling within a function: this necessarily
       cumulates with the incoming sample restriction
       recorded in state->submask
    */
    oldmask = state.submask;
    } else if (state != null) {
    /* FIXME? This is a loophole that potentially allows
       overriding of a (t1, t2) sample range at caller
       level, via "restrict replace" sampling within a
       function.
    */
    ;
    }

    if (err) {
    return err;
    }

    if (!replace && do_precompute(mode, panmask, oldmask, dset)) {
    /* we come here only if cumulating restrictions */
    mask = precompute_mask(param, panmask, oldmask, dset,
                               &new_t1, &new_t2, prn, &err);
        if (new_t1 >= 0 && new_t2 >= 0) {
            contig = 1;
            goto contig_next;
        }
    } else if (panmask != null) {
    /* got a panel time-sample mask */
    mask = copy_subsample_mask(panmask, &err);
    } else {
    /* no @panmask, and not already handled by "precompute" above */
        if (replace) {
            err = restore_full_sample(dset, state);
        }
    err = make_restriction_mask(mode, param, list, dset, oldmask,
                    replace, &mask, &new_t1, &new_t2,
                                    prn);
        if (new_t1 >= 0 && new_t2 >= 0) {
            contig = 1;
            goto contig_next;
        }
    }

    if (!err && !full_sample(dset)) {
    /* restore the full data range, for housekeeping purposes */
    err = restore_full_sample(dset, null);
    }

    if (err) {
    /* clean up and get out */
    free(mask);
    if (free_oldmask) {
        free(oldmask);
    }
    return err;
    }

    if (err && state != null && state.submask != null) {
    /* Should we try to restore the status quo ante on error?
       Too tricky, I think -- best to cancel the 'catch' flag
       if it's present.
    */
    if (state.cmd.flags & CmdFlags.catch_) {
        state.cmd.flags ^= CmdFlags.catch_;
        gretl_errmsg_append("Sorry, cannot 'catch' this error", err);
    }
    }

 contig_next:

    if (!err && n_models > 0) {
    err = check_models_for_subsample(mask, n_dropped);
    }

    if (!err && mask != null) {
    if (!contig && !(opt & Gretlopt.OPT_N)) {
        /* don't bother with this for the --random case */
        contig = mask_get_t1_t2(mask, dset, &new_t1, &new_t2);
    }
    if (contig) {
            err = handle_contiguous_sample(dset, new_t1, new_t2,
                                           mask, permanent, opt);
    } else {
        err = restrict_sample_from_mask(mask, dset, opt);
    }
    }

    free(mask);

    if (free_oldmask) {
    free(oldmask);
    }

    if (!err) {
    if (permanent) {
        dataset_clear_sample_record(dset);
    } else {
        make_restriction_string(dset, oldrestr, param, mode);
    }
    if (n_dropped != null) {
        *n_dropped = n_orig - sample_size(dset);
    }
    }

    free(oldrestr);

    return err;
}

/* restrict_sample:
 * @param: restriction string (or %NULL).
 * @list: list of variables in case of OPT_M (or %NULL).
 * @dset: dataset struct.
 * @state: structure representing program state (or %NULL).
 * @opt: option flags.
 * @prn: printing apparatus.
 * @n_dropped: location to receive count of dropped
 * observations, or NULL.
 *
 * Sub-sample the data set, based on the criterion of skipping all
 * observations with missing data values (OPT_M); or using as mask a
 * specified dummy variable (OPT_O); or masking with a specified
 * boolean condition (OPT_R); or selecting at random (OPT_N).
 *
 * In case OPT_M or OPT_C a @list of variables may be supplied; in
 * cases OPT_O, OPT_R and OPT_N, @param must contain specifics.
 *
 * In case OPT_P is included, the restriction will rePlace any
 * existing sample restriction, otherwise the resulting restriction
 * will be the logical product of the new restriction and any
 * existing restriction.
 *
 * In case the original dataset was a panel and OPT_B was given,
 * we'll pad with missing values if necessary, to try to reconstitute
 * a balanced panel.
 *
 * In case OPT_T is included, the sample restriction will be
 * permanent (the dataset is shrunk to the sample), otherwise the
 * restriction can be undone via the command "smpl full".
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int restrict_sample(const(char) * param, const(int) * list, DatasetRaw * dset, ExecState * state, Gretlopt opt, Prn * prn, int * n_dropped)
{
    return real_restrict_sample(param, list, null, dset, state,
                opt, prn, n_dropped);
}

/* perma_sample:
 * @dset: dataset struct.
 * @opt: option flags: must be OPT_T.
 * @prn: printing apparatus.
 * @n_dropped: location to receive count of dropped models,
 * or NULL.
 *
 * Make the current sub-sampling of the dataset permanent.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int perma_sample(DatasetRaw * dset, Gretlopt opt, Prn * prn, int * n_dropped)
{
    if (gretl_function_depth() > 0) {
    gretl_errmsg_set(_("The dataset cannot be modified at present").ptr);
    return GretlError.data;
    } else if (!dataset_is_subsampled(dset)) {
    pputs(prn, "smpl: nothing to be done\n");
    return 0;
    } else if (dset.submask == RESAMPLED) {
    pputs(prn, "smpl: dataset is resampled\n");
    return GretlError.data;
    } else if (opt != Gretlopt.OPT_T) {
    return GretlError.badopt;
    }

    recode_strvals(dset, Gretlopt.OPT_NONE);

    if (dset.submask == null) {
    /* just trim observations */
    return dataset_shrink_obs_range(dset);
    }

    if (n_dropped != null) {
    int err;

    err = check_models_for_subsample(dset.submask, n_dropped);
    if (err) {
        return err;
    }
    } else {
    check_models_for_subsample(dset.submask, null);
    }

    free(dset.submask);
    dset.submask = null;
    dataset_clear_sample_record(dset);

    if (fullset.varname == dset.varname) {
    fullset.varname = null;
    }
    if (fullset.varinfo == dset.varinfo) {
    fullset.varinfo = null;
    }
    if (fullset.descrip == dset.descrip) {
    fullset.descrip = null;
    }

    destroy_dataset(fullset);
    fullset = peerset = null;

    return 0;
}

enum SmplLimit {
    t1,
    t2
}

int smpl_get_int(const(char) * s, DatasetRaw * dset, int relative, int * err)
{
    int k = -1;

    if (integer_string(s)) {
        double k0 = dset.sd0;

        if (annual_data(dset) && !probably_year(s, dset)) {
            k = atoi(s);
        } else if (dated_daily_data(dset)) {
            k = atoi(s);
        } else if (!relative && k0 > 1.0 && k0 == floor(k0)) {
            /* allow for a starting obs such as 250 */
            k = atoi(s) - cast(int) k0 + 1;
        } else {
            k = atoi(s);
        }
    } else {
    double x;

    if (gretl_is_scalar(s)) {
        x = gretl_scalar_get_value(s, null);
    } else {
        x = generate_scalar(s, dset, err);
    }
    if (!na(x) && x < int.max) {
        k = cast(int) x;
    }
    }

    return k;
}

/* On receiving OPT_D to force interpretation of integer @t
   as representing a date, check that this works, and if so
   return the corresponding observation index. This is
   available only for annual data or calendar data, where
   @t is supplied as ISO 8601 basic, YYYYMMDD.
*/

int verify_integer_date(int t, DatasetRaw * dset, int * err)
{
    int ret = -1;

    if (annual_data(dset)) {
    if (t >= dset.sd0 && t < dset.sd0 + dset.n) {
        ret = cast(int) (t - dset.sd0);
    } else {
        *err = GretlError.pdwrong;
    }
    } else if (calendar_data(dset)) {
    guint32 ed = epoch_day_from_ymd_basic(t);
    char * date = null;

    if (ed <= 0) {
        *err = GretlError.pdwrong;
    } else {
        date = ymd_extended_from_epoch_day(ed, 0, err);
    }
    if (!*err) {
        ret = dateton(date, dset);
    }
    free(date);
    } else {
    *err = GretlError.pdwrong;
    }

    return ret;
}

/* Heuristic for interpreting @s as a year rather than a
   1-based observation index.
*/

int probably_year(const(char) * s, DatasetRaw * dset)
{
    int sval = atoi(s);

    return (annual_data(dset) && dset.sd0 > 1000 &&
        sval > 1000 && sval > dset.n);
}

int got_two_ints(const(char) * s, int * pi, int * pt)
{
    int n = 0;

    if (strchr(s, ':') != null) {
    n = sscanf(s, "%d:%d".ptr, pi, pt);
    } else if (strchr(s, '.') != null) {
    n = sscanf(s, "%d.%d".ptr, pi, pt);
    }

    return n == 2;
}

int get_sample_limit(const(char) * s, DatasetRaw * dset, int code, Gretlopt opt)
{
    int ret = -1;
    int i, t;
    int err = 0;

    if (dataset_is_panel(dset) && got_two_ints(s, &i, &t)) {
    if (code == SmplLimit.t1 && t != 1) {
        gretl_errmsg_sprintf(_("'%s': invalid starting observation").ptr, s);
        err = GretlError.data;
    } else if (code == SmplLimit.t2 && t != dset.pd) {
        gretl_errmsg_sprintf(_("'%s': invalid ending observation").ptr, s);
        err = GretlError.data;
    }
    if (err) {
        return -1;
    }
    }

    if (*s == '-' || *s == '+') {
    /* increment/decrement form */
    int incr = smpl_get_int(s + 1, dset, 1, &err);

    if (!err) {
        if (*s == '-') {
        incr = -incr;
        }
        if (dataset_is_panel(dset)) {
        incr *= dset.pd;
        }
        if (code == SmplLimit.t1) {
        ret = dset.t1 + incr;
        } else {
        ret = dset.t2 + incr;
        }
    }
    } else {
    /* absolute form */
    if (opt & Gretlopt.OPT_D) {
        /* force interpretation of integers as dates */
        int k = smpl_get_int(s, dset, 0, &err);

        if (!err) {
        ret = verify_integer_date(k, dset, &err);
        }
    } else {
        if (!integer_string(s) || probably_year(s, dset)) {
        ret = get_t_from_obs_string(s, dset);
        }
        if (ret < 0) {
        gretl_error_clear();
        ret = smpl_get_int(s, dset, 0, &err);
        /* convert to base 0 */
        if (!err) ret--;
        }
    }
    }

    return ret;
}

/* Catch the case where we're in a function and attempting to
   move t1 or t2 out of the range established on entry to the
   function: we don't want to carry forward any error message
   that implies @t was out of the full data range.
*/

void maybe_clear_range_error(int t, DatasetRaw * dset)
{
    if (t >= 0 && t < dset.n) {
    gretl_error_clear();
    }
}

int out_of_bounds_check(int * t, int which, const(char) * s, DatasetRaw * dset, int tmin, int tmax)
{
    int err = 0;

    if (*t < tmin || *t > tmax) {
        /* Prima facie there's an error, but could we be
           wrong about that? Probably only if @s is an
           integer string and we over-interpreted it. So
           try reading @s as a simple 1-based index.
        */
        err = GretlError.invarg;
        if (integer_string(s)) {
            int k = atoi(s) - 1;

            if (k >= tmin && k <= tmax) {
                *t = k;
                err = 0;
            }
        }
    }

    if (err) {
        maybe_clear_range_error(*t, dset);
        if (which == SmplLimit.t2) {
            gretl_errmsg_set(_("error in new ending obs").ptr);
        } else {
            gretl_errmsg_set(_("error in new starting obs").ptr);
        }
    }

    return err;
}

/**
 * set_sample:
 * @start: start of sample range.
 * @stop: end of sample range.
 * @dset: pointer to dataset struct.
 * @opt: options flag(s).
 *
 * Sets the current sample range for @dset to that specified by
 * @start and @stop, if possible.
 *
 * If @opt contains %OPT_T, make the sub-sampling permanent,
 * or in other words shrink the dataset to the selected range.
 *
 * Returns: 0 on successful completion, non-zero error code
 * otherwise.
 */

int set_sample(const(char) * start, const(char) * stop, DatasetRaw * dset, Gretlopt opt)
{
    Gretlopt opt_in = opt;
    int new_t1 = dset.t1;
    int new_t2 = dset.t2;
    int tmin = 0, tmax = 0;
    int nf, err = 0;

    if (dset == null) {
    return GretlError.nodata;
    }

    /* --quiet (OPT_Q), --permanent (OPT_T) and --dates (OPT_D)
       are the only acceptable options here
    */
    opt_in = cast(Gretlopt)(opt_in & ~Gretlopt.OPT_Q);
    opt_in = cast(Gretlopt)(opt_in & ~Gretlopt.OPT_T);
    opt_in = cast(Gretlopt)(opt_in & ~Gretlopt.OPT_D);
    if (opt_in != Gretlopt.OPT_NONE) {
    return GretlError.badopt;
    }

    gretl_error_clear();

    nf = (start != null) + (stop != null);

    if (nf == 2 && dset.n == 0) {
    /* database special */
    return db_set_sample(start, stop, dset);
    }

    if (nf == 0) {
    /* no-op, just print the current sample */
    return 0;
    }

    sample_range_get_extrema(dset, &tmin, &tmax);

    if (nf == 1) {
    /* implicitly just setting the starting observation? */
    if (start == null) {
        return GretlError.args;
    }
    new_t1 = get_sample_limit(start, dset, SmplLimit.t1, opt);
        err = out_of_bounds_check(&new_t1, SmplLimit.t1, start, dset, tmin, tmax);
    goto finish;
    }

    /* now we're looking at the 2 fields case */
    if (strcmp(start, ";")) {
    new_t1 = get_sample_limit(start, dset, SmplLimit.t1, opt);
        err = out_of_bounds_check(&new_t1, SmplLimit.t1, start, dset, tmin, tmax);
    }
    if (!err && strcmp(stop, ";")) {
    new_t2 = get_sample_limit(stop, dset, SmplLimit.t2, opt);
        err = out_of_bounds_check(&new_t2, SmplLimit.t2, stop, dset, tmin, tmax);
    }

    if (!err && new_t1 > new_t2) {
    gretl_error_clear();
    gretl_errmsg_set(_("Invalid null sample").ptr);
    err = GretlError.data;
    }

 finish:

    if (!err && dataset_is_panel(dset)) {
    int T = dset.pd;

    if (new_t1 % T != 0 || (new_t2 + 1) % T != 0) {
        gretl_errmsg_set("Invalid panel sample range");
        err = GretlError.invarg;
    }
    }

    if (!err) {
    dset.t1 = new_t1;
    dset.t2 = new_t2;
    }

    if (!err && (opt & Gretlopt.OPT_T)) {
    err = perma_sample(dset, opt, null, null);
    }

    return err;
}

int t_from_verified_aqm(const(char) * s, int pd)
{
    if (pd == 1) {
    return atoi(s);
    } else {
    return pd * atoi(s) + atoi(s + 5);
    }
}

bool panel_range_ok(int t1, int t2, int T)
{
    return (t1>=0 && t2<T && t1<=t2);
}

int panel_time_sample(const(char) * start, const(char) * stop, int t1, int t2, Gretlopt opt, DatasetRaw * dset, ExecState * s, Prn * prn)
{
    int T = dset.pd;
    int pd = dset.panel_pd;
    int replace = (opt & Gretlopt.OPT_P);
    int err = 0;

    if (replace && dset == peerset && fullset != null) {
    if (dataset_is_panel(fullset)) {
        restore_full_sample(dset, s);
        T = dset.pd;
        pd = dset.panel_pd;
        opt = cast(Gretlopt)(opt & ~Gretlopt.OPT_P);
    } else {
        gretl_errmsg_set("The full dataset is not a panel");
        return GretlError.data;
    }
    }

    if (panel_range_ok(t1, t2, T)) {
    ; /* alright */
    } else if (pd == 1 || pd == 4 || pd == 12) {
    DatasetRaw tset;
    int t0;

    time_series_from_panel(&tset, dset);
    t0 = t_from_verified_aqm(tset.stobs.ptr, tset.pd);
    t1 = obs_index_from_aqm(start, tset.pd, t0, T);
    t2 = obs_index_from_aqm(stop, tset.pd, t0, T);
    if (!panel_range_ok(t1, t2, T)) {
        err = GretlError.data;
    }
    } else if ((pd == 5 || pd == 6 || pd == 7 || pd == 52) &&
           dset.panel_sd0 > 10000) {
    DatasetRaw tset;

    time_series_from_panel(&tset, dset);
    t1 = calendar_obs_number(start, &tset, 0);
    t2 = calendar_obs_number(stop, &tset, 0);
    if (!panel_range_ok(t1, t2, T)) {
        err = GretlError.data;
    }
    } else {
    err = GretlError.data;
    }

    if (err) {
    gretl_errmsg_sprintf("Invalid sample range %s to %s\n".ptr,
                 start, stop);
    err = GretlError.data;
    } else {
    /* construct a mask with 1s indicating observations
       to include, 0s for those to be skipped
    */
    char * mask = make_submask(dset.n);

    if (mask == null) {
        err = GretlError.alloc;
    } else {
        int N = dset.n / T;
        int i, t, k = 0;

        for (i=0; i<N; i++) {
        for (t=0; t<T; t++) {
            if (t >= t1 && t <= t2) {
            mask[k] = 1;
            }
            k++;
        }
        }
        err = real_restrict_sample(null, null, mask, dset, s,
                       opt, prn, null);
        free(mask);
    }
    if (!err) {
        /* should perhaps be redundant? */
        dset.panel_pd = pd;
        dset.panel_sd0 = get_date_x(pd, start);
    }
    }

    return err;
}

int set_panel_sample(const(char) * start, const(char) * stop, Gretlopt opt, DatasetRaw * dset, ExecState * s, Prn * prn)
{
    Gretlopt testopt = opt;
    int s1 = -1, s2 = -1;
    int err = 0;

    testopt = cast(Gretlopt)(testopt & ~Gretlopt.OPT_U);
    testopt = cast(Gretlopt)(testopt & ~Gretlopt.OPT_X);
    testopt = cast(Gretlopt)(testopt & ~Gretlopt.OPT_Q);
    testopt = cast(Gretlopt)(testopt & ~Gretlopt.OPT_P);

    if (testopt != Gretlopt.OPT_NONE) {
    /* we accept only --unit or --time here, plus --quiet or --replace */
    return GretlError.badopt;
    } else if (incompatible_options(opt, Gretlopt.OPT_U | Gretlopt.OPT_X)) {
    /* cannot supply both --unit and --time in a single command */
    return GretlError.badopt;
    } else if (start == null || stop == null) {
    /* we must have start and stop values */
    return GretlError.parse;
    } else if (!dataset_is_panel(dset)) {
    gretl_errmsg_sprintf(_("%s: inapplicable option").ptr, print_flags(opt, GretlCmdIndex.smpl));
    return GretlError.badopt;
    }

    /* first pass: read start and stop as integer values? */
    if (strchr(start, ':') == null && strchr(start, '-') == null) {
    s1 = smpl_get_int(start, dset, 0, &err);
    }
    if (!err && strchr(stop, ':') == null && strchr(stop, '-') == null) {
    s2 = smpl_get_int(stop, dset, 0, &err);
    }

    if (opt & Gretlopt.OPT_X) {
    /* note that an error from smpl_get_int() above is not
       fatal in the --time case, since we might be looking
       at date strings
    */
    return panel_time_sample(start, stop, s1-1, s2-1, opt,
                 dset, s, prn);
    }

    if (!err) {
    if (s1 < 1) {
        gretl_errmsg_sprintf(_("invalid first unit %d").ptr, s1);
        err = GretlError.data;
    } else if (s2 > dset.n / dset.pd) {
        gretl_errmsg_sprintf(_("invalid last unit %d").ptr, s2);
        err = GretlError.data;
    } else if (s2 < s1) {
        gretl_errmsg_set(_("invalid null sample").ptr);
        err = GretlError.data;
    }
    }

    if (!err) {
    int t1 = (s1 - 1) * dset.pd;
    int t2 = s2 * dset.pd - 1;
    int tmin = 0, tmax = 0;

    sample_range_get_extrema(dset, &tmin, &tmax);
    if (t1 < tmin || t2 > tmax) {
        gretl_errmsg_set("sample range out of bounds");
        err = GretlError.data;
    } else {
        dset.t1 = t1;
        dset.t2 = t2;
    }
    }

    return err;
}

/**
 * count_missing_values:
 * @dset: dataset struct.
 * @opt: use %OPT_V for verbose operation, %OPT_A to
 * examine all data.
 * @prn: printing struct.
 * @err: location to receive error code.
 *
 * Prints a count of missing values (if any) in the current
 * dataset over the currently defined sample range (or the
 * entire data range if %OPT_A is given). If %OPT_V is given
 * this includes a count of missing values at each observation;
 * otherwise it just includes global and per-variable counts.
 *
 * Returns: 0 if no missing values are found (or on error),
 * otherwise the total number of missing values.
 */

int count_missing_values(const(DatasetRaw) * dset, Gretlopt opt, Prn * prn, int * err)
{
    int missval = 0, missobs = 0, totvals = 0, oldmiss = 0;
    int T, t1, t2;
    int * missvec;
    double missfrac;
    int i, t, tmiss;

    if (opt & Gretlopt.OPT_A) {
    t1 = 0;
    t2 = dset.n - 1;
    } else {
    t1 = dset.t1;
    t2 = dset.t2;
    }

    T = t2 - t1 + 1;

    missvec = cast(typeof(missvec)) malloc(dset.v * missvec.sizeof);

    if (missvec == null) {
    *err = GretlError.alloc;
    return 0;
    }

    for (i=0; i<dset.v; i++) {
    missvec[i] = 0;
    }

    for (t=t1; t<=t2; t++) {
    tmiss = 0;
    for (i=1; i<dset.v; i++) {
        if (series_is_hidden(dset, i)) {
        continue;
        }
        if (na(dset.Z[i][t])) {
        if (missvec[i] == 0) {
            missvec[0] += 1;
        }
        missvec[i] += 1;
        missval++;
        }
        totvals++;
    }

    tmiss = missval - oldmiss;

    if (tmiss) {
        missobs++;

        if (opt & Gretlopt.OPT_V) {
        /* verbose print by observation */
        if (dset.markers) {
            pprintf(prn, "%8s %4d %s\n".ptr, dset.S[t], tmiss,
                _("missing values").ptr);
        } else {
            char[OBSLEN] tmp;

            ntolabel(tmp.ptr, t, dset);
            pprintf(prn, "%8s %4d %s\n".ptr, tmp.ptr, tmiss,
                _("missing values").ptr);
        }
        }
    }
    oldmiss = missval;
    }

    missfrac = 100.0 * cast(double) missobs / T;

    pprintf(prn, _("\nNumber of observations (rows) with missing data values = %d (%.2f%%)\n").ptr, missobs, missfrac);

    pprintf(prn, _("Total number of missing data values = %d (%.2f%% of total data values)\n").ptr, missval,
        (100.0 * cast(double) missval / totvals));

    if (missvec[0] > 0) {
    pputc(prn, '\n');
    for (i=1; i<dset.v; i++) {
        if (missvec[i] > 0) {
        missfrac = 100.0 * cast(double) missvec[i] / T;
        pprintf(prn, "%8s: %d %s (%.2f%%); %d %s (%.2f%%)\n".ptr,
            dset.varname[i], missvec[i], _("missing values").ptr,
            missfrac, T - missvec[i], _("valid values").ptr,
            100.0 - missfrac);
        }
    }
    }

    free(missvec);

    return missval;
}

void copy_series_info(DatasetRaw * dest, const(DatasetRaw) * src, int maxv)
{
    int i;

    for (i=1; i<maxv; i++) {
    strcpy(dest.varname[i], src.varname[i]);
    if (src.varinfo != null) {
        copy_varinfo(dest.varinfo[i], src.varinfo[i]);
    }
    }
}

/* Below: For a model that was estimated using a data sample which
   differs from the current sample, reconstitute the dataset
   appropriate to the model.  If the model used a subsample, we use
   the special dummy variable (submask) saved with the model to select
   observations from the full dataset.  Another possibility is that
   the model used the full dataset, which has by now been subsampled.

   We attach this dataset to the model, so that it can be used for,
   e.g., residual or fitted versus actual plots.  The attached data
   set will be destroyed when the model itself is destroyed.

   By default we create a dataset containing series up to the
   highest numbered series associated with the model (regressand
   and regressors).  If OPT_F is given we include all series; if
   OPT_G is given we include only the constant.
*/

int add_dataset_to_model(Model * pmod, const(DatasetRaw) * dset, Gretlopt opt)
{
    const(DatasetRaw) * srcset;
    char * mask = null;
    int maxv, sn = 0;

    if (gretl_is_between_model(pmod)) {
    /* special case: we can't reconstruct the dataset here,
       but it may already be attached to @pmod.
    */
    if (pmod.dataset != null) {
        return 0;
    } else {
        return GretlError.data;
    }
    }

    if (pmod.dataset != null) {
    /* FIXME? */
    destroy_dataset(pmod.dataset);
    pmod.dataset = null;
    }

    if (fullset != null) {
    sync_dataset_members(dset);
    srcset = fullset;
    } else {
    srcset = dset;
    }

    if (pmod.submask == null) {
    /* no subsample info: pmod was estimated on the full dataset,
       so we'll reconstruct the full dataset */
    sn = srcset.n;
    } else if (pmod.submask == RESAMPLED) {
    /* can't reconstruct: resampling doesn't correspond to a
       mask over the original dataset */
    return GretlError.data;
    } else {
    /* pmod was estimated on a subsample, which has to
       be reconstructed */
    int t;

    mask = cast(typeof(mask)) calloc(srcset.n, 1);
    if (mask == null) {
        return GretlError.alloc;
    }
    for (t=0; t<srcset.n; t++) {
        if (pmod.submask[t] == 1 || pmod.submask[t] == 'p') {
        /* preserve the panel-padding marker, if any, so that
           copy_data_to_subsample() fills padding rows with
           NAs rather than genuine source-dataset values */
        mask[t] = pmod.submask[t];
        sn++;
        }
    }
    if (sn == 0) {
        free(mask);
        return 1;
    }
    }

    if (opt & Gretlopt.OPT_F) {
    maxv = srcset.v;
    } else if (opt & Gretlopt.OPT_G) {
    maxv = 1;
    } else {
    maxv = highest_numbered_var_in_model(pmod, dset) + 1;
    }

    /* allocate auxiliary dataset */
    pmod.dataset = create_auxiliary_dataset(maxv, sn, Gretlopt.OPT_NONE);
    if (pmod.dataset == null) {
    free(mask);
    return GretlError.alloc;
    }

    /* copy across info on series */
    copy_series_info(pmod.dataset, srcset, maxv);

    /* copy across data */
    copy_data_to_subsample(pmod.dataset, srcset, maxv, mask);

    /* dataset characteristics such as pd: if we're rebuilding the
       full dataset copy these across; but if we're reconstructing a
       subsampled dataset these features from the full dataset will be
       wrong, and we stay with the simple defaults
    */
    if (pmod.submask == null) {
    pmod.dataset.pd = srcset.pd;
    pmod.dataset.sd0 = srcset.sd0;
    strcpy(pmod.dataset.stobs.ptr, srcset.stobs.ptr);
    strcpy(pmod.dataset.endobs.ptr, srcset.endobs.ptr);
    pmod.dataset.structure = srcset.structure;
    }

    free(mask);

    return 0;
}

int submasks_match(const(DatasetRaw) * dset, const(Model) * pmod)
{
    char * s1 = cast(char *) dset.submask;
    char * s2 = cast(char *) pmod.submask;

    if (s1 == RESAMPLED && s2 == RESAMPLED) {
    return dset.n == pmod.full_n &&
        dset.rseed == pmod.smpl.rseed;
    } else if (s1 == RESAMPLED || s2 == RESAMPLED) {
    return 0;
    } else {
    return submask_cmp(s1, s2) == 0;
    }
}

/* check the subsample mask from a model against @dset to
   see if it may have been estimated on a different
   (subsampled) data set from the current one
*/

int model_sample_problem(const(Model) * pmod, const(DatasetRaw) * dset)
{
    int ret = 1;

    if (pmod.ci == GretlCmdIndex.panel && (pmod.opt & Gretlopt.OPT_B)) {
    ; /* special case: panel "between" model */
    } else if (pmod.submask == null) {
    /* the model has no sub-sampling info recorded */
    if (dset.submask == null) {
        /* data set is not sub-sampled either, OK */
        ret = 0;
    } else {
        fputs("dataset is subsampled, model is not\n", stderr);
        gretl_errmsg_set(_("dataset is subsampled, model is not\n").ptr);
        ret = 1;
    }
    } else {
    /* the model does have sub-sampling info recorded */
    if (dset.submask == null) {
        fputs("model is subsampled, dataset is not\n", stderr);
        gretl_errmsg_set(_("model is subsampled, dataset is not\n").ptr);
        ret = 1;
    } else if (submasks_match(dset, pmod)) {
        /* the subsamples (model and current data set) agree, OK */
        ret = 0;
    } else {
        /* the subsamples differ */
        gretl_errmsg_set(_("model and dataset subsamples not the same\n").ptr);
        ret = 1;
    }
    }

    return ret;
}

int fcast_not_feasible(const(Model) * pmod, const(DatasetRaw) * dset)
{
    int ret = GretlError.data;

    if (pmod.ci == GretlCmdIndex.panel && (pmod.opt & Gretlopt.OPT_B)) {
    ; /* special case: panel "between" model */
    } else if (pmod.submask == null) {
    /* @pmod was estimated on the full dataset */
    if (dset.submask == null) {
        /* @dset not sub-sampled either, OK */
        ret = 0;
    } else if (dataset_is_cross_section(dset)) {
        /* @dset is cross-sectional subset, OK? */
        ret = 0;
    } else {
        gretl_errmsg_set(_("dataset is subsampled, model is not\n").ptr);
    }
    } else {
    /* the model does have sub-sampling info recorded */
    if (dset.submask == null) {
        /* the dataset is not sub-sampled */
        if (dataset_is_cross_section(dset)) {
        ret = 0; /* OK? */
        } else if (is_panel_time_sample(pmod.submask, dset)) {
        ret = 0; /* OK? */
        } else {
        gretl_errmsg_set(_("model is subsampled, dataset is not\n").ptr);
        }
    } else if (submasks_match(dset, pmod)) {
        /* the subsamples (model and current data set) agree, OK */
        ret = 0;
    } else {
        /* the subsamples differ */
        if (dataset_is_cross_section(fullset) &&
        dataset_is_cross_section(dset)) {
        ret = 0; /* OK ? */
        } else if (is_panel_time_sample(pmod.submask, fullset) &&
               is_panel_time_sample(dset.submask, fullset)) {
        /* we should be able to make this work */
        ret = 0; /* new_ 2019-03-23 */
        } else {
        gretl_errmsg_set(_("model and dataset subsamples not the same\n").ptr);
        }
    }
    }

    return ret;
}

int same_dataset(const(Model) * pmod, const(DatasetRaw) * dset)
{
    if (pmod.submask == null && dset.submask == null) {
    return 1;
    } else if (pmod.submask != null && dset.submask != null) {
    return submasks_match(dset, pmod);
    } else {
    return 0;
    }
}

void dataset_type_string(char * str, const(DatasetRaw) * dset)
{
    if (dataset_is_time_series(dset)) {
    strcpy(str, _("time series").ptr);
    } else if (dataset_is_panel(dset)) {
        strcpy(str, _("panel").ptr);
    } else {
        strcpy(str, _("undated").ptr);
    }
}

void pd_string(char * str, const(DatasetRaw) * dset)
{
    if (custom_time_series(dset)) {
    strcpy(str, _("special").ptr);
    } else {
    switch (dset.pd) {
    case 1:
        strcpy(str, _("annual").ptr); break;
    case 4:
        strcpy(str, _("quarterly").ptr); break;
    case 12:
        strcpy(str, _("monthly").ptr); break;
    case 24:
        strcpy(str, _("hourly").ptr); break;
    case 52:
        strcpy(str, _("weekly").ptr); break;
    case 5:
        strcpy(str, _("daily (5 days)").ptr); break;
    case 6:
        strcpy(str, _("daily (6 days)").ptr); break;
    case 7:
        strcpy(str, _("daily (7 days)").ptr); break;
    case 10:
        strcpy(str, _("decennial").ptr); break;
    default:
        strcpy(str, _("unknown").ptr); break;
    }
    }
}

void print_sample_obs(const(DatasetRaw) * dset, Prn * prn)
{
    char[OBSLEN] d1;
    char[OBSLEN] d2;

    ntolabel(d1.ptr, dset.t1, dset);
    ntolabel(d2.ptr, dset.t2, dset);

    pprintf(prn, "%s: %s - %s".ptr, _("Current sample").ptr, d1.ptr, d2.ptr);
    pprintf(prn, " (n = %d)\n".ptr, dset.t2 - dset.t1 + 1);
}

void print_sample_status(const(DatasetRaw) * dset, Prn * prn)
{
    char[128] tmp;

    if (complex_subsampled()) {
    pprintf(prn, "%s\n".ptr, _("Full dataset").ptr);
    dataset_type_string(tmp.ptr, fullset);
    pprintf(prn, "%s: %s\n".ptr, _("Type").ptr, tmp.ptr);
    if (dataset_is_time_series(fullset)) {
        pd_string(tmp.ptr, fullset);
        pprintf(prn, "%s: %s\n".ptr, _("Frequency").ptr, tmp.ptr);
    } else if (dataset_is_panel(fullset)) {
        int nu = fullset.n / fullset.pd;

        pprintf(prn, "%s: %d\n".ptr, _("Number of cross-sectional units").ptr, nu);
        pprintf(prn, "%s: %d\n".ptr, _("Number of time periods").ptr, fullset.pd);
    }
    pprintf(prn, "%s: %s - %s (n = %d)\n".ptr, _("Range").ptr,
        fullset.stobs.ptr, fullset.endobs.ptr, fullset.n);

    pprintf(prn, "\n%s\n".ptr, _("Subsampled data").ptr);
    if (dset.restriction != null) {
        pprintf(prn, "%s: %s\n".ptr, _("restriction").ptr, dset.restriction);
    }
    } else {
    pprintf(prn, "%s\n".ptr, _("Full dataset").ptr);
    }

    dataset_type_string(tmp.ptr, dset);
    pprintf(prn, "%s: %s\n".ptr, _("Type").ptr, tmp.ptr);
    if (dataset_is_time_series(dset)) {
    pd_string(tmp.ptr, dset);
    pprintf(prn, "%s: %s\n".ptr, _("Frequency").ptr, tmp.ptr);
    } else if (dataset_is_panel(dset)) {
    int nu = dset.n / dset.pd;

    pprintf(prn, "%s: %d\n".ptr, _("Number of cross-sectional units").ptr, nu);
    pprintf(prn, "%s: %d\n".ptr, _("Number of time periods").ptr, dset.pd);
    }
    if (dset.t1 == 0 && dset.t2 == dset.n - 1) {
    pprintf(prn, "%s: %s - %s (n = %d)\n".ptr, _("Range").ptr,
        dset.stobs.ptr, dset.endobs.ptr, dset.n);
    } else {
    pprintf(prn, "%s: %s - %s (n = %d)\n".ptr, _("Range").ptr,
        dset.stobs.ptr, dset.endobs.ptr, dset.n);
    print_sample_obs(dset, prn);
    if (dset.restriction != null) {
        pprintf(prn, "(%s: %s)\n".ptr, _("restriction").ptr, dset.restriction);
    }
    }
}

/**
 * data_report:
 * @dset: data information struct.
 * @fname: filename for current datafile, or NULL.
 * @prn: gretl printing struct.
 *
 * Write out a summary of the content of the current data set.
 *
 * Returns: 0 on successful completion, non-zero on error.
 */

int data_report(const(DatasetRaw) * dset, const(char) * fname, Prn * prn)
{
    char[OBSLEN] startobs;
    char[OBSLEN] endobs;
    char[MAXLEN] tmp;
    const(char) * vlabel;
    int i;

    ntolabel(startobs.ptr, 0, dset);
    ntolabel(endobs.ptr, dset.n - 1, dset);

    if (fname != null) {
        /* the original code for this function */
        char[48] tstr;

        sprintf(tmp.ptr, _("Data file %s\nas of").ptr,
                (*fname != '\0')? fname : _("(unsaved)").ptr);
        print_time(tstr.ptr);
        pprintf(prn, "%s %s\n\n".ptr, tmp.ptr, tstr.ptr);

        if (dset.descrip != null && *dset.descrip != '\0') {
            pprintf(prn, "%s:\n\n".ptr, _("Description").ptr);
            pputs(prn, dset.descrip);
            pputs(prn, "\n\n");
        }
    }

    dataset_type_string(tmp.ptr, dset);
    pprintf(prn, "%s: %s\n".ptr, _("Type of data").ptr, tmp.ptr);

    if (dataset_is_time_series(dset)) {
    pd_string(tmp.ptr, dset);
    pprintf(prn, "%s: %s\n".ptr, _("Frequency").ptr, tmp.ptr);
    }

    if (dataset_is_panel(dset)) {
        pprintf(prn, "%s: N = %d, T = %d\n".ptr, _("Panel data").ptr,
                dset.n / dset.pd, dset.pd);
    } else {
        pprintf(prn, "%s: %s - %s (n = %d)\n\n".ptr, _("Range").ptr,
                startobs.ptr, endobs.ptr, dset.n);
    }

    pprintf(prn, "%s:\n\n".ptr, _("Listing of variables").ptr);

    for (i=1; i<dset.v; i++) {
    vlabel = series_get_label(dset, i);
    pprintf(prn, "%*s  %s\n".ptr, VNAMELEN, dset.varname[i],
        vlabel == null ? "".ptr : vlabel);
    }

    return 0;
}

