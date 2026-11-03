/* Port of lib/src/gretl_panel.c: estimation and testing for panel data (pooled OLS with panel diagnostics,
   fixed and random effects, between, WLS/ML groupwise models, panel-robust and clustered covariance matrices,
   Hausman, Breusch-Pagan, Wooldridge and related tests, the panel datasets utilities; the header gretl_panel.h
   has only prototypes).

   Differences from the C:
   - the PDEBUG, CDEBUG and BETWEEN_DEBUG diagnostics (compiled out in the C) and the `#if 0` code are not
     ported;
   - file-internal types are UnitTs, PanelmodT, ClusterInfo and PanIdx (the anonymous struct of the panel
     indexation, now the __gshared variable panidx); the anonymous vcv_ops enum is VcvOps;
   - the small function-like macros are functions; the "for testing purposes" globals are __gshared.
   Dependencies not yet ported are declared extern (C) below, in the style of the other ports. */

module gretl.gretl_panel;


import core.stdc.errno : errno;
import core.stdc.math : exp, fabs, floor, isnan, log, pow, sqrt, ceil;
import core.stdc.stdio : fprintf, sprintf, stderr;
import core.stdc.stdlib : calloc, free, malloc, qsort, realloc;
import core.stdc.string : memcpy, memset, strcmp, strcpy, strlen;

import std.algorithm : max, min;

import gretl.clapack_doubleh;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_f2ch : integer;
import gretl.gretl_matrix;
import gretl.gretl_matrixh;
import gretl.libgretlh;
import gretl.missingh : NADBL, na;

import core.stdc.string : strchr, strncat, strncmp;
import core.stdc.math : sin, asin;
import std.math : PI;
import core.stdc.stdio : sscanf;
import core.stdc.stdlib : atoi;
import core.stdc.string : strspn;
import core.stdc.math : nearbyint;
import gretl.gstrfuncs : g_strndup;
import gretl.calendar : calendar_date_string;
import gretl.dataio : ntolabel;
import gretl.dataset : create_auxiliary_dataset, dataset_add_NA_series, dataset_add_allocated_series, dataset_add_series, dataset_destroy_obs_markers, dataset_drop_last_variables, destroy_dataset, is_string_valued, series_attach_string_table, series_get_string_table, series_get_string_vals, set_panel_groups_name;
import gretl.dataseth : DatasetMarkerType, calendar_data, dataset_is_panel, sample_size;
import gretl.describe : gretl_corr, gretl_mean, gretl_variance;
import gretl.estimate : lsq, makevcv;
import gretl.geneval0 : MIN;
import gretl.genmain : current_series_index, generate, series_greatest_index, series_index;
import gretl.genparseh : ConstId;
import gretl.gmem : g_free;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_errors : errmsg, gretl_errmsg_set, gretl_errmsg_sprintf, gretl_error_clear;
import gretl.gretl_list : gretl_list_add, gretl_list_copy, gretl_list_delete_at_pos, gretl_list_drop, gretl_list_insert_list, gretl_list_new, gretl_list_omit, gretl_list_omit_last, gretl_list_sublist, in_gretl_list;
import gretl.gretl_listh : LISTSEP;
import gretl.gretl_model : ModelTest, clear_model, gretl_model_add_panel_varnames, gretl_model_allocate_param_names, gretl_model_destroy_data_item, gretl_model_get_double, gretl_model_get_int, gretl_model_get_list, gretl_model_init, gretl_model_set_cluster_vcv_info, gretl_model_set_data, gretl_model_set_double, gretl_model_set_hac_order, gretl_model_set_int, gretl_model_set_list_as_data, gretl_model_set_matrix_as_data, gretl_model_set_param_name, gretl_model_set_string_as_data, gretl_model_set_vcv_info, gretl_model_smpl_init, gretl_model_test_print_direct, gretl_model_write_vcv, maybe_add_test_to_model, maybe_suppress_time_dummies, mle_criteria, model_test_new, model_test_set_dfd, model_test_set_dfn, model_test_set_pvalue, model_test_set_teststat, model_test_set_value, set_model_id;
import gretl.gretl_modelh : ModelTestType, PanelVCVType, VCVMajorType;
import gretl.gretl_prn : Prn, pprintf, pputc, pputs;
import gretl.gretl_string_table : SeriesTable, series_table_get_strings, series_table_new;
import gretl.gretl_utils : combine_categories, copyvec, count_distinct_values, doubles_array_new, gretl_compare_doubles, gretl_isconst, ijton, ls_criteria, rearrange_id_array, record_matrix_test_result, record_test_result;
import gretl.gretl_utilsh : floateq;
import gretl.gstrfuncs : g_strdup, g_utf8_strlen;
import gretl.libset : get_hac_lag, gretl_echo_on, libset_get_int, libset_set_int;
import gretl.libseth : PanelRobust, SetKey;
import gretl.matrix_extra : gretl_matrix_read_from_file;
import gretl.options : get_optval_string;
import gretl.pvalues : chisq_cdf_comp, normal_pvalue_2, snedecor_cdf_comp, student_pvalue_2;
import gretl.strutils : gretl_lower, gretl_strdup, gretl_strdup_printf, gretl_string_split, gretl_string_split_quoted, has_suffix, strings_array_dup, strings_array_free, strings_array_new, strings_array_new_with_length;
import gretl.transforms : diffgenr, laggenr;
import gretl.modelprint : print_model_from_matrices, printmodel;
import gretl.subsample : complex_subsampled, fetch_full_dataset;
import gretl.compare : wald_omit_chisq, wald_omit_test;
import gretl.gretl_array : get_array_by_name, gretl_array_destroy, gretl_array_get_strings, gretl_array_new, gretl_array_nullify_elements, gretl_array_set_string;
import gretl.uservar : get_string_by_name, gretl_is_string;
import gretl.genfuncs : gen_panel_dummies, imhof;

/* unported dependencies: per-file declarations from the 2026c headers */
extern (C) {
}
enum VcvOps {
    init,
    subtract
}

struct UnitTs {
    int * t0; /* start of longest contiguous time series, per unit */
    int * T; /* length of such sequence, per unit */
}

struct PanelmodT {
    Gretlopt opt; /* option flags */
    int nunits; /* total cross-sectional units in sample range */
    int effn; /* effective (included) cross-section units */
    int T; /* times-series length of panel */
    int Tmax; /* maximum number of usable obs per unit */
    int Tmin; /* minimum number of usable obs per unit > 0 */
    double Tbar; /* harmonic mean of per-unit time-series lengths */
    int NT; /* total observations used (based on pooled model) */
    int ntdum; /* number of time dummies added */
    int * unit_obs; /* array of number of observations per x-sect unit */
    char * varying; /* array to record properties of pooled-model regressors */
    int * vlist; /* list of time-varying variables from pooled model */
    int balanced; /* 1 if the model dataset is balanced, else 0 */
    int dfH; /* number of coeffs for Hausman test */
    int Fdfn; /* numerator df, F for differing intercepts */
    double Fdfd; /* denominator df, F for differing intercepts */
    double theta; /* quasi-demeaning coefficient */
    double theta_bar; /* mean of theta_i (unbalanced case) */
    double Ffe; /* joint significance of fixed effects */
    double BP; /* Breusch-Pagan test statistic */
    double H; /* Hausman test statistic */
    GretlMatrixRaw * bdiff; /* array of coefficient differences */
    GretlMatrixRaw * Sigma; /* Hausman covariance matrix */
    double s2b; /* residual variance, group means regression */
    double s2e; /* residual variance, fixed-effects regression */
    double s2v; /* estimate of between variance */
    double ubPub; /* for use in unbalanced Swamy-Arora */
    double dw; /* Durbin-Watson, for transfer to random effects */
    double rho; /* for transfer to random effects */
    int * small2big; /* data indexation array */
    int * big2small; /* reverse data indexation array */
    Model * pooled; /* reference model (pooled OLS) */
    Model * realmod; /* fixed or random effects model */
    double * re_uhat; /* "fixed" random-effects residuals */
}

struct PanIdx {
    int n; /* number of cross-sectional units */
    int T; /* number of observations per unit */
    int offset; /* sampling offset into full dataset */
}

__gshared PanIdx panidx;

/* for testing purposes */
__gshared int stata_sa;
__gshared int IGLS;
__gshared char[MAXLEN] glsmat;


/* translate from (i = unit, t = time period for that unit) to
   overall 0-based index into the data set */
int panel_index(int i, int t)
{
    return i * panidx.T + t + panidx.offset;
}

bool panel_missing(PanelmodT * p, int t)
{
    return na(p.pooled.uhat[t]);
}


void panel_index_init(const(DatasetRaw) * dset, int nunits, int T)
{
    panidx.n = nunits;
    panidx.T = T;
    panidx.offset = dset.t1;

}

/* Allocate the indexation arrays that allow us to translate between
   the full-length data array and the possibly shorter array used for
   transformed data (de-meaned or quasi-demeaned).
*/

int allocate_data_finders(PanelmodT * pan, int bign)
{
    int s;

    if (pan.small2big != null) {
        /* already done */
        return 0;
    }

    pan.small2big = cast(typeof(pan.small2big)) malloc(pan.NT * (*pan.small2big).sizeof);
    pan.big2small = cast(typeof(pan.big2small)) malloc(bign * (*pan.big2small).sizeof);

    if (pan.small2big == null || pan.big2small == null) {
        return GretlError.alloc;
    }

    for (s=0; s<bign; s++) {
        pan.big2small[s] = -1;
    }

    return 0;
}

int small_index(PanelmodT * p, int t)
{
    return (p.big2small == null)? t : p.big2small[t];
}
int big_index(PanelmodT * p, int t)
{
    return (p.small2big == null)? t : p.small2big[t];
}

void panelmod_init(PanelmodT * pan)
{
    pan.nunits = 0;
    pan.effn = 0;
    pan.T = 0;
    pan.Tmax = 0;
    pan.Tmin = 0;
    pan.Tbar = 0;
    pan.NT = 0;
    pan.ntdum = 0;
    pan.unit_obs = null;
    pan.varying = null;
    pan.vlist = null;
    pan.opt = Gretlopt.OPT_NONE;

    pan.balanced = 1;
    pan.theta = NADBL;
    pan.theta_bar = NADBL;
    pan.ubPub = NADBL;
    pan.dw = NADBL;
    pan.rho = NADBL;

    pan.Ffe = NADBL;
    pan.Fdfn = 0;
    pan.Fdfd = 0;

    pan.BP = NADBL;
    pan.H = NADBL;
    pan.dfH = 0;

    pan.bdiff = null;
    pan.Sigma = null;

    pan.small2big = null;
    pan.big2small = null;

    pan.pooled = null;
    pan.realmod = null;
    pan.re_uhat = null;
}

void panelmod_free(PanelmodT * pan)
{
    free(pan.unit_obs);
    free(pan.varying);
    free(pan.vlist);

    gretl_matrix_free(pan.bdiff);
    gretl_matrix_free(pan.Sigma);

    free(pan.small2big);
    free(pan.big2small);
    free(pan.re_uhat);

    free(pan.realmod);
}

/* test variable number v against the (possibly reduced) regression
   list which contains only time-varying regressors
*/

int var_is_varying(const(PanelmodT) * pan, int v)
{
    if (v != 0) {
        int i;

        for (i=2; i<=pan.vlist[0]; i++) {
            if (pan.vlist[i] == v) {
                return 1;
            }
        }
    }

    return 0;
}

/* retrieve X'X^{-1} from the pooled or fixed effects model */

GretlMatrixRaw * panel_model_xpxinv(Model * pmod, int * err)
{
    GretlMatrixRaw * X;
    int k = pmod.ncoeff;
    double x;
    int i, j, m;


    if (gretl_model_get_int(pmod, "vcv_xpx")) {
        /* already done */
        gretl_model_set_int(pmod, "vcv_xpx", 0);
    } else if (pmod.vcv != null) {
        double s2 = pmod.sigma * pmod.sigma;
        int nv = (k * k + k) / 2;

        for (i=0; i<nv; i++) {
            pmod.vcv[i] /= s2;
        }
    } else {
        *err = makevcv(pmod, 1.0);
        if (*err) {
            return null;
        }
        gretl_model_set_int(pmod, "vcv_xpx", 1);
    }

    X = gretl_matrix_alloc(k, k);
    if (X == null) {
        *err = GretlError.alloc;
        return null;
    }

    m = 0;
    for (j=0; j<k; j++) {
        for (i=j; i<k; i++) {
            x = pmod.vcv[m++];
            gretl_matrix_set(X, i, j, x);
            gretl_matrix_set(X, j, i, x);
        }
    }


    return X;
}

/* Beck and Katz, as outlined in Greene.  We offer the following
   for pooled OLS and FE.  Greene writes (in effect):

   Var(b) = A^{-1} W A^{-1}

   where A = \sum_{i=1}^n X'_i X_i (called "XX" below)
         W = \sum_{i=1}^n \sum_{j=1}^n \sigma_{ij} X'_{i} X_{j}

   and \sigma_{ij} is estimated as (1/T) \sum_{t=1}^T e_i e_j,
   with the e's being OLS (or FE) residuals.

   In computing this, we need to check that W is positive definite:
   this seems not to be guaranteed for unbalanced panels.
*/

int beck_katz_vcv(Model * pmod, PanelmodT * pan, const(DatasetRaw) * dset, GretlMatrixRaw * XX, GretlMatrixRaw * W, GretlMatrixRaw * V)
{
    GretlMatrixRaw * Xi = null;
    GretlMatrixRaw * Xj = null;
    int T = pan.T;
    int k = pmod.ncoeff;
    int s, si, sj;
    int i, j, p, v, t;
    int err = 0;

    Xi = gretl_matrix_alloc(pan.Tmax, k);
    Xj = gretl_matrix_alloc(pan.Tmax, k);

    if (Xi == null || Xj == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    for (i=0; i<pan.nunits; i++) {
        int Ti = pan.unit_obs[i];
        double sii = 0.0;

        if (Ti == 0) {
            continue;
        }

        Xi = gretl_matrix_reuse(Xi, Ti, k);

        s = 0;
        for (t=0; t<T; t++) {
            si = panel_index(i, t);
            if (!na(pmod.uhat[si])) {
                sii += pmod.uhat[si] * pmod.uhat[si];
                si = small_index(pan, si);
                for (p=0; p<k; p++) {
                    v = pmod.list[p + 2];
                    gretl_matrix_set(Xi, s, p, dset.Z[v][si]);
                }
                s++;
            }
        }

        sii /= Ti;
        sii = sqrt(sii);

        /* "diagonal" component of W */
        gretl_matrix_multiply_by_scalar(Xi, sii);
        gretl_matrix_multiply_mod(Xi, GretlMatrixMod.transpose,
                                  Xi, GretlMatrixMod.none,
                                  W, GretlMatrixMod.cumulate);

        for (j=i+1; j<pan.nunits; j++) {
            int Tij = 0;
            double sij = 0.0;

            if (pan.unit_obs[j] == 0) {
                continue;
            }

            /* count matching observations */
            for (t=0; t<T; t++) {
                si = panel_index(i, t);
                sj = panel_index(j, t);
                if (!na(pmod.uhat[si]) && !na(pmod.uhat[sj])) {
                    sij += pmod.uhat[si] * pmod.uhat[sj];
                    Tij++;
                }
            }

            if (Tij == 0) {
                continue;
            }

            sij /= Tij;
            Xi = gretl_matrix_reuse(Xi, Tij, k);
            Xj = gretl_matrix_reuse(Xj, Tij, k);

            s = 0;
            for (t=0; t<T; t++) {
                si = panel_index(i, t);
                sj = panel_index(j, t);
                if (!na(pmod.uhat[si]) && !na(pmod.uhat[sj])) {
                    si = small_index(pan, si);
                    sj = small_index(pan, sj);
                    for (p=0; p<k; p++) {
                        v = pmod.list[p + 2];
                        gretl_matrix_set(Xi, s, p, dset.Z[v][si]);
                        gretl_matrix_set(Xj, s, p, dset.Z[v][sj]);
                    }
                    s++;
                }
            }

            /* cumulate s_ij * Xi'Xj into W (using V as storage) */
            gretl_matrix_multiply_by_scalar(Xi, sij);
            gretl_matrix_multiply_mod(Xi, GretlMatrixMod.transpose,
                                      Xj, GretlMatrixMod.none,
                                      V, GretlMatrixMod.none);
            gretl_matrix_add_to(W, V);
            /* and take in s_ij * Xj'Xi */
            gretl_matrix_add_transpose_to(W, V);
        }
    }

    /* check that the middle term is p.d. */
    gretl_matrix_copy_values(V, W);
    err = gretl_matrix_cholesky_decomp(V);
    if (err) {
        fprintf(stderr, "beck_katz_vcv: matrix W is not p.d.\n".ptr);
        gretl_model_set_int(pmod, "panel_bk_failed", 1);
        goto bailout;
    }

    /* form V = (Xi'Xi)^{-1} W (Xi'Xi)^{-1} */
    gretl_matrix_qform(XX, GretlMatrixMod.none, W,
                       V, GretlMatrixMod.none);

    gretl_model_set_vcv_info(pmod, VCVMajorType.panel, PanelVCVType.bk);

 bailout:

    gretl_matrix_free(Xi);
    gretl_matrix_free(Xj);

    return err;
}

/* helper function for panel VCV clustered by period */

int * get_panel_tobs(PanelmodT * pan)
{
    int * tobs = cast(int *) calloc(pan.T, (int).sizeof);

    if (tobs != null) {
        int i, t;

        for (i=0; i<pan.nunits; i++) {
            for (t=0; t<pan.T; t++) {
                if (!panel_missing(pan, panel_index(i, t))) {
                    tobs[t] += 1;
                }
            }
        }
    }

    return tobs;
}

void augment_clustered_W(const(GretlMatrixRaw) * e, const(GretlMatrixRaw) * X, GretlMatrixRaw * eX, GretlMatrixRaw * W)
{
    gretl_matrix_multiply_mod(e, GretlMatrixMod.transpose,
                              X, GretlMatrixMod.none,
                              eX, GretlMatrixMod.none);
    gretl_matrix_multiply_mod(eX, GretlMatrixMod.transpose,
                              eX, GretlMatrixMod.none,
                              W, GretlMatrixMod.cumulate);
}

void finalize_clustered_vcv(Model * pmod, PanelmodT * pan, const(GretlMatrixRaw) * XX, const(GretlMatrixRaw) * W, GretlMatrixRaw * V, int nc)
{
    double adj;

    /* form V(b) = (X'X)^{-1} W (X'X)^{-1} */
    gretl_matrix_qform(XX, GretlMatrixMod.none, W,
                       V, GretlMatrixMod.none);

    /* Now, what (if anything) are we going to do in terms of
       degrees-of-freedom adjustment?
    */

    if (pmod.opt & Gretlopt.OPT_N) {
        /* --no-df-corr */
        return;
    }

    if (pmod.ci == GretlCmdIndex.ivreg) {
        /* compatible with gretl 2023b: do we want this? */
        return;
    }

    /* initial adjustment using the number of clusters */
    adj = nc / (nc - 1.0);

    if (pmod.ci == GretlCmdIndex.ivreg || pmod.opt & Gretlopt.OPT_N) {
        /* Just apply the @adj calculated above? This is said to
           be the "asymtotic-like" approach in Stata; see
           https://www.stata.com/meeting/13uk/nichols_crse.pdf,
           under "Finite-Sample Adjustment". Right now this
           is never reached!
        */
        ;
    } else {
        /* Apply the full Cameron-Gelbach-Miller adjustment */
        adj *= (pmod.nobs - 1.0) / (pmod.nobs - pmod.ncoeff);
    }

    gretl_matrix_multiply_by_scalar(V, adj);
}

enum CI_TIME = -2;
enum CI_UNIT = (-3);

struct ClusterInfo {
    const(double) * cz; /* series of all cluster-var values */
    GretlMatrixRaw * cvals; /* sorted vector of unique cluster-var values */
    int[2] dcid; /* IDs of up to two @dset series */
    int[3] pcid; /* IDs of up to three @pset series */
    int[2] nc; /* number(s) of clusters */
    byte target; /* 0 or 1 for first or second cluster var */
    byte pooled; /* flag for the pooled OLS case */
}

bool by_time_and_unit(ClusterInfo * ci)
{
    return ci.dcid[0] < -1 && ci.dcid[1] < -1;
}
bool by_time_or_unit(ClusterInfo * ci)
{
    return ci.dcid[0] < -1 || ci.dcid[1] < -1;
}
bool by_time(ClusterInfo * ci)
{
    return ci.dcid[0] == CI_TIME;
}
bool by_unit(ClusterInfo * ci)
{
    return ci.dcid[0] == CI_UNIT;
}
bool generic(ClusterInfo * ci)
{
    return ci.dcid[0] > 0;
}
bool two_way(ClusterInfo * ci)
{
    return ci.dcid[0] != -1 && ci.dcid[1] != -1;
}
bool is_present(int v)
{
    return v != -1;
}

ClusterInfo * cluster_info_new(PanelmodT * pan)
{
    ClusterInfo * ci = cast(ClusterInfo *) calloc(1, (ClusterInfo).sizeof);

    ci.dcid[0] = ci.dcid[1] = -1;
    ci.target = 0;
    ci.pooled = (pan.opt & Gretlopt.OPT_P)? 1 : 0;

    return ci;
}

void cluster_info_free(ClusterInfo * ci)
{
    if (ci != null) {
        gretl_matrix_free(ci.cvals);
        free(ci);
    }
}

const(char) * cluster_name(int v, const(DatasetRaw) * dset)
{
    if (v > 0) {
    return dset.varname[v];
    } else if (v == CI_TIME) {
    return "time";
    } else {
    return "unit";
    }
}

int time_cluster_vcv(Model * pmod, PanelmodT * pan, const(DatasetRaw) * dset, const(GretlMatrixRaw) * XX, GretlMatrixRaw * W, GretlMatrixRaw * V, ClusterInfo * ci)
{
    GretlVector * et = null;
    GretlMatrixRaw * Xt = null;
    GretlVector * eXt = null;
    int * tobs = null;
    int N = pan.effn;
    int k = pmod.ncoeff;
    int idx = ci.target;
    int i, j, v, s, t;
    int n_c = 0;
    int err = 0;

    tobs = get_panel_tobs(pan);
    if (tobs == null) {
        return GretlError.alloc;
    }

    et  = gretl_column_vector_alloc(N);
    Xt  = gretl_matrix_alloc(N, k);
    eXt = gretl_vector_alloc(k);

    if (et == null || Xt == null || eXt == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    for (t=0; t<pan.T; t++) {
        int Nt = tobs[t]; /* # of units represented in period @t */
        int p = 0;

        if (Nt == 0) {
            continue;
        }

        n_c++;
        et = gretl_matrix_reuse(et, Nt, 1);
        Xt = gretl_matrix_reuse(Xt, Nt, k);

        for (i=0; i<pan.nunits; i++) {
            s = panel_index(i, t);
            if (na(pmod.uhat[s])) {
                continue;
            }
            gretl_vector_set(et, p, pmod.uhat[s]);
            s = small_index(pan, s);
            for (j=0; j<k; j++) {
                v = pmod.list[j+2];
                if (s < 0) {
                    gretl_matrix_set(Xt, p, j, 0);
                } else {
                    gretl_matrix_set(Xt, p, j, dset.Z[v][s]);
                }
            }
            p++;
            if (p == Nt) {
                /* we're done with this period */
                break;
            }
        }
        augment_clustered_W(et, Xt, eXt, W);
    }

    finalize_clustered_vcv(pmod, pan, XX, W, V, pan.Tmax);
    if (!two_way(ci)) {
    /* finish the job */
    gretl_model_set_vcv_info(pmod, VCVMajorType.panel, PanelVCVType.time);
    gretl_model_set_int(pmod, "n_clusters", n_c);
    } else {
    ci.nc[idx] = n_c;
    }

 bailout:

    gretl_matrix_free(et);
    gretl_matrix_free(Xt);
    gretl_matrix_free(eXt);
    free(tobs);

    return err;
}

int panel_cval_count(Model * pmod, PanelmodT * pan, ClusterInfo * ci, int i, const(double) * cvar)
{
    double cvi = ci.cvals.val[i];
    int s, t, cc = 0;

    if (cvar != null) {
    /* the pooled OLS case */
    for (t=pmod.t1; t<=pmod.t2; t++) {
        if (!na(pmod.uhat[t]) && cvar[t] == cvi) {
        cc++;
        }
    }
    } else {
    /* other cases */
    for (s=0; s<pmod.nobs; s++) {
        t = big_index(pan, s);
        if (!na(pmod.uhat[t]) && ci.cz[s] == cvi) {
        cc++;
        }
    }
    }

    return cc;
}

int panel_cval_count_max(Model * pmod, PanelmodT * pan, ClusterInfo * ci, const(double) * cvar)
{
    int n = gretl_vector_get_length(ci.cvals);
    int i, cc, cmax = 0;

    for (i=0; i<n; i++) {
        cc = panel_cval_count(pmod, pan, ci, i, cvar);
        if (cc > cmax) {
            cmax = cc;
        }
    }

    return cmax;
}

int generic_cluster_vcv(Model * pmod, PanelmodT * pan, const(DatasetRaw) * dset, const(GretlMatrixRaw) * XX, GretlMatrixRaw * W, GretlMatrixRaw * V, ClusterInfo * ci)
{
    const(double) * cvar = null;
    GretlVector * ei = null;
    GretlMatrixRaw * Xi = null;
    GretlVector * eXi = null;
    int R, M = ci.cvals.rows;
    int k = pmod.ncoeff;
    int idx = ci.target;
    int i, j, v, t;
    int n_c = 0;
    int err = 0;

    if (ci.pooled) {
    cvar = dset.Z[ci.dcid[idx]];
    }

    R = panel_cval_count_max(pmod, pan, ci, cvar);
    ei  = gretl_column_vector_alloc(R);
    Xi  = gretl_matrix_alloc(R, k);
    eXi = gretl_vector_alloc(k);

    if (ei == null || Xi == null || eXi == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    for (i=0; i<M; i++) {
        double cvi = ci.cvals.val[i];
        int Ni = panel_cval_count(pmod, pan, ci, i, cvar);
        int s, p = 0;

        if (Ni == 0) {
            continue;
        }

        ei = gretl_matrix_reuse(ei, Ni, -1);
        Xi = gretl_matrix_reuse(Xi, Ni, -1);

    if (ci.pooled) {
        /* working with the full dataset */
        for (t=pmod.t1; t<=pmod.t2; t++) {
        if (!na(pmod.uhat[t]) && cvar[t] == cvi) {
            gretl_vector_set(ei, p, pmod.uhat[t]);
            for (j=0; j<k; j++) {
            v = pmod.list[j+2];
            gretl_matrix_set(Xi, p, j, dset.Z[v][t]);
            }
            p++;
        }
        if (p == Ni) {
            break;
        }
        }
    } else {
        /* working with dataset from which NAs have been purged */
        for (s=0; s<pmod.nobs; s++) {
        t = big_index(pan, s);
        if (!na(pmod.uhat[t]) && ci.cz[s] == cvi) {
            gretl_vector_set(ei, p, pmod.uhat[t]);
            for (j=0; j<k; j++) {
            v = pmod.list[j+2];
            gretl_matrix_set(Xi, p, j, dset.Z[v][s]);
            }
            p++;
        }
        if (p == Ni) {
            break;
        }
        }
    }
        augment_clustered_W(ei, Xi, eXi, W);
        n_c++;
    }

    finalize_clustered_vcv(pmod, pan, XX, W, V, n_c);

    if (!two_way(ci)) {
        /* not doing two-way clustering, so finish the job */
    const(char) * cname = dset.varname[ci.pcid[idx]];

        gretl_model_set_cluster_vcv_info(pmod, cname, null);
        gretl_model_set_int(pmod, "n_clusters", n_c);
    } else if (idx < 2) {
    ci.nc[idx] = n_c;
    }

 bailout:

    gretl_matrix_free(ei);
    gretl_matrix_free(Xi);
    gretl_matrix_free(eXi);

    return err;
}

int panel_white_vcv(Model * pmod, PanelmodT * pan, const(DatasetRaw) * dset, const(GretlMatrixRaw) * XX, GretlMatrixRaw * W, GretlMatrixRaw * V, ClusterInfo * ci)
{
    GretlVector * ei = null;
    GretlMatrixRaw * Xi = null;
    GretlVector * eXi = null;
    int k = pmod.ncoeff;
    int i, v, s, t;
    int n_c = 0;
    int err = 0;

    ei  = gretl_column_vector_alloc(1);
    Xi  = gretl_matrix_alloc(1, k);
    eXi = gretl_vector_alloc(k);

    if (ei == null || Xi == null || eXi == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    for (s=0; s<pmod.nobs; s++) {
    t = big_index(pan, s);
    if (!na(pmod.uhat[t])) {
        ei.val[0] = pmod.uhat[t];
        for (i=0; i<k; i++) {
        v = pmod.list[i+2];
        gretl_vector_set(Xi, i, dset.Z[v][s]);
        }
        augment_clustered_W(ei, Xi, eXi, W);
        n_c++;
    }
    }

    finalize_clustered_vcv(pmod, pan, XX, W, V, n_c);
    if (!two_way(ci)) {
    /* should we offer "plain White" (HC0) as an option? */
    gretl_model_set_vcv_info(pmod, VCVMajorType.hc, 0);
    }

 bailout:

    gretl_matrix_free(ei);
    gretl_matrix_free(Xi);
    gretl_matrix_free(eXi);

    return err;
}

/* HAC covariance matrix for pooled, fixed- or random-effects models,
   given "fixed T and large N".  In the case of "large T" a different
   form is needed for robustness in respect of autocorrelation.  See
   Arellano, "Panel Data Econometrics" (Oxford, 2003), pages 18-19.
*/

int unit_cluster_vcv(Model * pmod, PanelmodT * pan, const(DatasetRaw) * dset, const(GretlMatrixRaw) * XX, GretlMatrixRaw * W, GretlMatrixRaw * V, ClusterInfo * ci)
{
    GretlVector * ei = null;
    GretlMatrixRaw * Xi = null;
    GretlVector * eXi = null;
    int T = pan.Tmax;
    int k = pmod.ncoeff;
    int idx = ci.target;
    int i, j, v, s, t;
    int err = 0;

    ei  = gretl_column_vector_alloc(T);
    Xi  = gretl_matrix_alloc(T, k);
    eXi = gretl_vector_alloc(k);

    if (ei == null || Xi == null || eXi == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    for (i=0; i<pan.nunits; i++) {
        int Ti = pan.unit_obs[i];
        int p = 0;

        if (Ti == 0) {
            continue;
        }

        ei = gretl_matrix_reuse(ei, Ti, 1);
        Xi = gretl_matrix_reuse(Xi, Ti, k);

        for (t=0; t<pan.T; t++) {
            s = panel_index(i, t);
            if (na(pmod.uhat[s])) {
                continue;
            }
            gretl_vector_set(ei, p, pmod.uhat[s]);
            s = small_index(pan, s);
            for (j=0; j<k; j++) {
                v = pmod.list[j+2];
                if (s < 0) {
                    gretl_matrix_set(Xi, p, j, 0);
                } else {
                    gretl_matrix_set(Xi, p, j, dset.Z[v][s]);
                }
            }
            if (++p == Ti) {
                /* we're done with this unit */
                break;
            }
        }
        augment_clustered_W(ei, Xi, eXi, W);
    }

    finalize_clustered_vcv(pmod, pan, XX, W, V, pan.effn);

    if (!two_way(ci)) {
        /* finish the job if not doing two-way clustering */
    gretl_model_set_vcv_info(pmod, VCVMajorType.panel, PanelVCVType.hac);
    gretl_model_set_int(pmod, "n_clusters", pan.effn);
    } else {
    ci.nc[idx] = pan.effn;
    }

 bailout:

    gretl_matrix_free(ei);
    gretl_matrix_free(Xi);
    gretl_matrix_free(eXi);

    return err;
}

/* In response to "$time" or "$unit" appearing in the context of the
   --cluster option for a basic panel-data estimator, we construct an
   on-the-fly series that represents the time period -- or if @srcv
   equals CI_UNIT, the cross-sectional unit. We need this only in
   case of two-way clustering with one of $time or $unit combined
   with a named series.
*/

void make_unit_or_period(Model * pmod, PanelmodT * pan, double * y, int srcv)
{
    int do_unit = (srcv == CI_UNIT);
    int i, t, s;

    for (i=0; i<pan.nunits; i++) {
        if (pan.unit_obs[i] == 0) {
        continue;
    }
    for (t=0; t<pan.T; t++) {
        s = panel_index(i, t);
        if (na(pmod.uhat[s])) {
        continue;
        }
        s = small_index(pan, s);
        y[s] = do_unit ? i+1 : t+1;
    }
    }
}

int cluster_missval_error(DatasetRaw * dset, int v, int t)
{
    char[OBSLEN] obs;

    gretl_errmsg_sprintf("cluster variable '%s' missing at observation %s".ptr,
             dset.varname[v], ntolabel(obs.ptr, t, dset));
    return GretlError.missdata;
}

/* The role of this function is to transcribe a clustering variable
   into the panel-estimation dataset, @pset, the destination series ID
   being given by ci->pcid. As for the source of the data, there are
   two cases:

   (1) a named series present in the main dataset, @dset; or

   (2) an automatically generated series, responding to the keyword
   "$time" or "$unit".

   The second case arises only in the case of "mixed" two-way
   clustering (where one of the variables comes from @dset and the
   other is automatic). If the cluster specification contains only
   automatic variables there's no need for transcription.

   On entry to this function ci->target should be set to either 0,
   indicating that we're working on the first (or only) cluster
   variable, or 1 to indicate that we're working on the second.
*/

int transcribe_cluster_var(Model * pmod, PanelmodT * pan, DatasetRaw * pset, DatasetRaw * dset, ClusterInfo * ci)
{
    const(double) * src = null;
    double * dest = null;
    int idx = ci.target;
    int srcv, t, s;
    int err = 0;

    /* source and destination of transcription */
    srcv = ci.dcid[idx];
    dest = pset.Z[ci.pcid[idx]];


    if (srcv == CI_TIME || srcv == CI_UNIT) {
    /* the "automatic" case */
    make_unit_or_period(pmod, pan, dest, srcv);
    return 0;
    }

    src = dset.Z[srcv];

    if (ci.pooled) {
    for (t=pmod.t1; t<=pmod.t2; t++) {
            if (!na(pmod.uhat[t]) && na(src[t])) {
        err = cluster_missval_error(dset, srcv, t);
        break;
            }
        }
    } else if (pset.n == dset.n) {
        /* balanced */
        for (t=0; t<dset.n; t++) {
            if (na(src[t])) {
                err = cluster_missval_error(dset, srcv, t);
                break;
            }
        }
        if (!err) {
            memcpy(dest, src, pset.n * (**pset.Z).sizeof);
        }
    } else {
    /* unbalanced */
        for (s=0; s<pmod.nobs && !err; s++) {
            t = big_index(pan, s);
            if (!na(pmod.uhat[t])) {
                if (na(src[t])) {
                    err = cluster_missval_error(dset, srcv, t);
                } else {
                    dest[s] = src[t];
                }
            }
        }
    }

    if (!err) {
    int nread = pmod.nobs;

        if (ci.cvals != null) {
            gretl_matrix_free(ci.cvals);
        }
    if (ci.pooled) {
        dest += pmod.t1;
        nread = pmod.t2 - pmod.t1 + 1;
    }
        ci.cvals = gretl_matrix_values(dest, nread, Gretlopt.OPT_S, &err);
        ci.cz = dest;
    }


    return err;
}

/* For a clustering variable named by @s, find out if it's a named
   series, or a keyword ("$time", "$unit"), or neither (an error). To
   signal the presence of one of the keywords we assign one of CI_TIME
   or CI_UNIT in place of a regular series ID number.
*/

int get_id_or_special(const(char) * s, const(DatasetRaw) * dset, ClusterInfo * ci)
{
    int * pid = &ci.dcid[ci.target];
    int id = current_series_index(dset, s);
    int err = 0;

    if (id > 0) {
    *pid = id;
    } else if (id == 0) {
    /* "const" will not do! */
    err = GretlError.data;
    } else if (!strcmp(s, "$time")) {
    *pid = CI_TIME;
    } else if (!strcmp(s, "$unit")) {
    *pid = CI_UNIT;
    } else {
    err = GretlError.unkvar;
    }

    return err;
}

int process_cluster_string(const(char) * s, const(DatasetRaw) * dset, ClusterInfo * ci)
{
    const(char) * p = strchr(s, ',');
    int err = 0;


    if (p == null) {
    /* no comma: should have a single name */
    ci.target = 0;
    err = get_id_or_special(s, dset, ci);
    } else {
    /* should have two names */
    char * s1 = null;
    char * s2 = null;

    s1 = g_strndup(s, p - s);
    ci.target = 0;
    err = get_id_or_special(s1, dset, ci);
    if (!err) {
        s2 = g_strdup(p + 1);
        ci.target = 1;
        err = get_id_or_special(s2, dset, ci);
        ci.target = 0;
    }
    g_free(s1);
    g_free(s2);
    }

    return err;
}

/* In this function we're looking for one or two identifiers supplied
   via the string parameter to the --cluster option.  We accept
   "$unit" and "$time" as automatic references to the cross-sectional
   or time dimension of the panel; otherwise names of series in @dset
   are needed.

   If this succeeds, we transcribe the first (or only) clustering
   series into the panel dataset, @pset.  If there's a second argument
   we record its ID as dcid[1] under the @ci struct; it will get
   transcribed later, in two_way_cluster_vcv().
*/

int check_cluster_var(Model * pmod, PanelmodT * pan, DatasetRaw * pset, DatasetRaw * dset, ClusterInfo * ci)
{
    /* note: OPT_L means that we came here from an invocation of
       "ols", so the relevant command index for retrieval of an
       option parameter is OLS and not PANEL
    */
    int opt_ci = (pan.opt & Gretlopt.OPT_L)? GretlCmdIndex.ols : GretlCmdIndex.panel;
    const(char) * s = get_optval_string(cast(GretlCmdIndex) opt_ci, Gretlopt.OPT_C);
    int err = 0;


    if (s == null || *s == '\0') {
        err = GretlError.data;
    } else {
    err = process_cluster_string(s, dset, ci);
    }


    if (err) {
    return err;
    } else if (by_time_and_unit(ci)) {
    ; /* nothing more to be done */
    } else if (by_time_or_unit(ci) && !two_way(ci)) {
    ; /* nothing more to be done */
    } else {
    int[3] adds = 0;
        int v, n_add = 0;

    if (ci.pooled && ci.dcid[0] > 0) {
        /* using an existing series for first var */
        ci.pcid[0] = ci.dcid[0];
    } else {
        adds[0] = 1;
    }
    if (is_present(ci.dcid[1])) {
        if (pset == dset && ci.dcid[1] > 0) {
        /* we'll just need one extra series for the combination */
        ci.pcid[1] = ci.dcid[1];
        adds[2] = 1;
        } else {
        /* we'll need more series plus the combination */
        adds[1] = adds[2] = 1;
        }
    }

    n_add = adds[0] + adds[1] + adds[2];

    if (n_add > 0) {
        err = dataset_add_series(pset, n_add);
        if (err) {
        return err;
        }
        v = pset.v - n_add;
        if (adds[0]) {
        ci.pcid[0] = v;
        strcpy(pset.varname[v++], cluster_name(ci.dcid[0], dset));
        }
        if (adds[1]) {
        ci.pcid[1] = v;
        strcpy(pset.varname[v++], cluster_name(ci.dcid[1], dset));
        }
        if (adds[2]) {
        ci.pcid[2] = v;
        strcpy(pset.varname[v], "combo");
        }
    }

    /* transcribe (if needed) the first cluster var */
    err = transcribe_cluster_var(pmod, pan, pset, dset, ci);
    }

    return err;
}

/* Fix-up, à la Cameron, Gelbach and Miller, for non-positive
   semi-definite covariance matrix in case of non-nested
   two-way clustering: set any negative eigenvalues to zero and
   reconstitute @V.
*/

int maybe_eigenfix(GretlMatrixRaw * V)
{
    GretlMatrixRaw * U = gretl_matrix_copy(V);
    GretlMatrixRaw * lam = null;
    int n = V.rows;
    int i, fixit = 0;
    int err = 0;

    if (U == null) {
        err = GretlError.alloc;
    } else {
        lam = gretl_symmetric_matrix_eigenvals(U, 1, &err);
    }

    for (i=0; i<n && !err; i++) {
        if (lam.val[i] < 0) {
            fixit = 1;
            break;
        }
    }

    if (fixit) {
        GretlMatrixRaw * Tmp = gretl_matrix_copy(U);
        double uij, lvj;
        int j;

        if (Tmp == null) {
            err = GretlError.alloc;
            goto bailout;
        }
        for (j=0; j<n; j++) {
            lvj = lam.val[j];
            for (i=0; i<n; i++) {
                if (lvj > 0) {
                    uij = gretl_matrix_get(Tmp, i, j);
                    gretl_matrix_set(Tmp, i, j, uij * lvj);
                } else {
                    gretl_matrix_set(Tmp, i, j, 0);
                }
            }
        }
        gretl_matrix_multiply_mod(Tmp, GretlMatrixMod.none,
                                  U, GretlMatrixMod.transpose,
                                  V, GretlMatrixMod.none);
        gretl_matrix_free(Tmp);
    }

 bailout:

    gretl_matrix_free(lam);
    gretl_matrix_free(U);

    return err;
}

int two_way_cluster_vcv(Model * pmod, PanelmodT * pan, DatasetRaw * pset, DatasetRaw * dset, const(GretlMatrixRaw) * XX, GretlMatrixRaw * W, GretlMatrixRaw * V, ClusterInfo * ci)
{
    GretlMatrixRaw * Tmp;
    int k = V.rows;
    int * v = ci.pcid.ptr;
    int err = 0;

    Tmp = gretl_matrix_alloc(k, k);
    if (Tmp == null) {
        return GretlError.alloc;
    }

    /* prepare for first step */
     ci.target = 0;

    /* Step 1: calculate @V for the first cluster var */
    if (ci.dcid[0] == CI_TIME) {
    err = time_cluster_vcv(pmod, pan, pset, XX, W, V, ci);
    } else if (ci.dcid[0] == CI_UNIT) {
    err = unit_cluster_vcv(pmod, pan, pset, XX, W, V, ci);
    } else {
    err = generic_cluster_vcv(pmod, pan, pset, XX, W, V, ci);
    }

    /* prepare for second step */
    ci.target = 1;

    if (!err && !by_time_and_unit(ci)) {
        /* Step 2a: transcribe the second cluster var if needed */
        err = transcribe_cluster_var(pmod, pan, pset, dset, ci);
    }

    if (!err) {
        /* Step 2b: calculate variance for the second cluster var */
        gretl_matrix_zero(W);
    if (ci.dcid[1] == CI_TIME) {
        err = time_cluster_vcv(pmod, pan, pset, XX, W, Tmp, ci);
    } else if (ci.dcid[1] == CI_UNIT) {
        err = unit_cluster_vcv(pmod, pan, pset, XX, W, Tmp, ci);
    } else {
        err = generic_cluster_vcv(pmod, pan, pset, XX, W, Tmp, ci);
    }
    /* add second variance to first */
    gretl_matrix_add_to(V, Tmp);
    }

    if (!err && by_time_and_unit(ci)) {
    /* Step 3a: calculate variance for combination (= plain White) */
    gretl_matrix_zero(W);
    err = panel_white_vcv(pmod, pan, pset, XX, W, Tmp, ci);
    goto step4;
    }

    if (!err) {
        /* Step 3b: combine the two cluster vars */
        err = combine_categories(pset, v[0], v[1], v[2]);
    }

    if (!err) {
        /* Step 3c: update ci->cz and ci->cvals, then calculate
           variance for the combination
        */
        gretl_matrix_free(ci.cvals);
        ci.cz = pset.Z[ci.pcid[2]];
        ci.cvals = gretl_matrix_values(ci.cz, pmod.nobs, Gretlopt.OPT_S, &err);
        gretl_matrix_zero(W);
    ci.target = 2;
        err = generic_cluster_vcv(pmod, pan, pset, XX, W, Tmp, ci);
    ci.target = 0;
    }

 step4:

    if (!err) {
        gretl_matrix_subtract_from(V, Tmp);
    err = maybe_eigenfix(V);
    }

    if (!err) {
        int nc = MIN(ci.nc[0], ci.nc[1]);

    if (by_time_and_unit(ci)) {
        gretl_model_set_vcv_info(pmod, VCVMajorType.panel, PanelVCVType.both);
        gretl_model_set_int(pmod, "n_clusters", nc);
    } else {
        gretl_model_set_cluster_vcv_info(pmod,
                         pset.varname[v[0]],
                         pset.varname[v[1]]);
        gretl_model_set_int(pmod, "n_clusters", nc);
    }
    }

    gretl_matrix_free(Tmp);

    return err;
}

/* The method of Driscoll and Kraay in "Consistent Covariance Matrix
   Estimation with Spatially Dependent Panel Data", REStat 80:4, 1998,
   pp. 549-560. See also Daniel Hoechle, "Robust Standard Errors for
   Panel Regressions with Cross-Sectional Dependence", Stata Journal
   7:3, 2007, for discussion of his implementation called xtscc.
*/

int driscoll_kraay_vcv(Model * pmod, PanelmodT * pan, const(DatasetRaw) * dset, const(GretlMatrixRaw) * XX, GretlMatrixRaw * W, GretlMatrixRaw * V)
{
    GretlMatrixBlock * B;
    GretlMatrixRaw * H = null;
    GretlMatrixRaw * ht = null;
    GretlMatrixRaw * htj = null;
    GretlMatrixRaw * Wj = null;
    GretlMatrixRaw * S = null;
    double bw; /* Bartlett weight */
    int * tobs = null;
    int k = pmod.ncoeff;
    int i, j, m, vj, s, t;
    int n_c = 0;
    int err = 0;

    tobs = get_panel_tobs(pan);
    if (tobs == null) {
        return GretlError.alloc;
    }

    B = gretl_matrix_block_new(&H, pan.T, k,
                   &S, k, k,
                   &ht, k, 1,
                   &htj, k, 1,
                   &Wj, k, k,
                   null);
    if (B == null) {
    free(tobs);
    return GretlError.alloc;
    }

    gretl_matrix_zero(H);
    gretl_matrix_zero(S);

    /* Maximum lag for Newey-West. Note: do "set hac_lag nw2"
       for agreement with Stata's xtscc */
    m = get_hac_lag(pan.Tmax);

    /* build the H matrix */
    for (t=0; t<pan.T; t++) {
        int Nt = tobs[t]; /* # of units represented in period @t */
        double eit, htjv, hplus;
        int p = 0;

        if (Nt == 0) {
            continue;
        }
    n_c++;
        for (i=0; i<pan.nunits; i++) {
            s = panel_index(i, t);
            eit = pmod.uhat[s];
            if (na(eit)) {
                continue;
            }
            s = small_index(pan, s);
            if (s >= 0) {
                for (j=0; j<k; j++) {
                    /* cumulate xit*eit into row @t of @H */
                    vj = pmod.list[j+2];
                    hplus = dset.Z[vj][s] * eit;
                    htjv = gretl_matrix_get(H, t, j);
                    gretl_matrix_set(H, t, j, htjv + hplus);
                }
            }
            if (++p == Nt) {
                /* we're done with this period */
                break;
            }
        }
    }

    /* compute initial @S = Omega_0 */
    for (t=0; t<pan.T; t++) {
        for (i=0; i<k; i++) {
            ht.val[i] = gretl_matrix_get(H, t, i);
        }
        gretl_matrix_multiply_mod(ht, GretlMatrixMod.none,
                                  ht, GretlMatrixMod.transpose,
                                  S, GretlMatrixMod.cumulate);
    }

    /* cumulate the weighted cross-lag terms */
    for (j=1; j<=m; j++) {
        gretl_matrix_zero(Wj);
        for (t=j; t<pan.T; t++) {
            for (i=0; i<k; i++) {
                ht.val[i] = gretl_matrix_get(H, t, i);
                htj.val[i] = gretl_matrix_get(H, t-j, i);
            }
            gretl_matrix_multiply_mod(ht, GretlMatrixMod.none,
                                      htj, GretlMatrixMod.transpose,
                                      Wj, GretlMatrixMod.cumulate);
        }
        /* add Barlett weight * (Wj + Wj') to @S */
        bw = 1.0 - j / (m + 1.0);
        gretl_matrix_add_self_transpose(Wj);
        gretl_matrix_multiply_by_scalar(Wj, bw);
        gretl_matrix_add_to(S, Wj);
    }

    finalize_clustered_vcv(pmod, pan, XX, S, V, pan.Tmax);
    gretl_model_set_vcv_info(pmod, VCVMajorType.panel, PanelVCVType.dk);
    gretl_model_set_hac_order(pmod, m);
    gretl_model_set_int(pmod, "DKT", n_c);

    gretl_matrix_block_destroy(B);
    free(tobs);

    return err;
}

/* Common setup for Arellano, Beck-Katz, Driscoll-Kraay and other
   clustered VCV estimators. For fixed and random effects, @pset is
   the special panel dataset comprising only usable observations, and
   @dset is the original dataset; for pooled OLS @pset and @dset are
   one and the same.
*/

int panel_robust_vcv(Model * pmod, PanelmodT * pan, DatasetRaw * pset, DatasetRaw * dset)
{
    ClusterInfo * ci = null;
    GretlMatrixRaw * W = null;
    GretlMatrixRaw * V = null;
    GretlMatrixRaw * XX = null;
    int k = pmod.ncoeff;
    int pan_robust = 0;
    int err = 0;

    XX = panel_model_xpxinv(pmod, &err);
    if (err) {
        fprintf(stderr, "panel_robust_vcv: failed at panel_model_xpxinv\n".ptr);
        goto bailout;
    }

    ci = cluster_info_new(pan);


    if (pan.opt & Gretlopt.OPT_C) {
        err = check_cluster_var(pmod, pan, pset, dset, ci);
        if (err) {
            goto bailout;
        }
    } else {
        pan_robust = libset_get_int(SetKey.panel_robust);
    }

    W = gretl_zero_matrix_new(k, k);
    V = gretl_matrix_alloc(k, k);

    if (W == null || V == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    if (two_way(ci)) {
    err = two_way_cluster_vcv(pmod, pan, pset, dset, XX, W, V, ci);
    } else if (by_time(ci)) {
        err = time_cluster_vcv(pmod, pan, pset, XX, W, V, ci);
    } else if (by_unit(ci)) {
    err = unit_cluster_vcv(pmod, pan, pset, XX, W, V, ci);
    } else if (generic(ci)) {
        err = generic_cluster_vcv(pmod, pan, pset, XX, W, V, ci);
    } else if (pan_robust == PanelRobust.beck_katz) {
        err = beck_katz_vcv(pmod, pan, pset, XX, W, V);
    } else if (pan_robust == PanelRobust.driscoll_kraay) {
        driscoll_kraay_vcv(pmod, pan, pset, XX, W, V);
    } else {
    /* default: pan_robust == PanelRobust.arellano */
        err = unit_cluster_vcv(pmod, pan, pset, XX, W, V, ci);
    }


    if (!err) {
        gretl_model_write_vcv(pmod, V);
    }

 bailout:

    gretl_matrix_free(W);
    gretl_matrix_free(V);
    gretl_matrix_free(XX);
    if (ci != null) {
        cluster_info_free(ci);
    }

    if (err && pmod.vcv != null) {
        free(pmod.vcv);
        pmod.vcv = null;
    }

    return err;
}

void femod_regular_vcv(Model * pmod)
{
    if (pmod.vcv == null) {
        /* estimated via Cholesky: no vcv yet */
        makevcv(pmod, pmod.sigma);
    } else {
        /* estimated via QR: "vcv" = (X'X)^{-1} */
        int i, k = pmod.ncoeff;
        int n = k * (k + 1) / 2;
        double s2 = pmod.sigma * pmod.sigma;

        for (i=0; i<n; i++) {
            pmod.vcv[i] *= s2;
        }
    }
}


int panel_DW_pval_ok(const(Model) * pmod)
{
    if (na(pmod.dw)) {
        return 0;
    } else {
/* c2d: CHECK #if DWPVAL_TESTING */
        return 1;
    }
}

/* See Bhargava, Franzini and Narendranathan, "Serial Correlation and
   the Fixed Effects Model", Review of Economic Studies 49, 1982,
   page 536. Strictly speaking what's calculated here is the marginal
   significance level of the DW stat when considered as a d_L value.
   It's unclear if this method can really be extended to unbalanced
   panels.
*/

double BFN_panel_DW_pvalue(Model * pmod, const(DatasetRaw) * dset, int * err)
{
    GretlMatrixRaw * lam = null;
    double r, pv, lamq, sinarg, pi_2T;
    int T = gretl_model_get_int(pmod, "Tmax");
    int N = gretl_model_get_int(pmod, "n_included_units");
    int nlam, k = pmod.ncoeff;
    int i, q;

    if (pmod.ifc) {
        k--; /* don't include the constant */
    }

    nlam = pmod.nobs - N - k;
    lam = gretl_column_vector_alloc(nlam);
    if (lam == null) {
        *err = GretlError.alloc;
        return NADBL;
    }

    pi_2T = PI / (2.0*T);
    sinarg = sin(pi_2T);
    lamq = 4 * sinarg * sinarg;
    r = pmod.dw;

    q = 1;
    for (i=0; i<nlam; i++) {
        lam.val[i] = lamq - r;
        if ((i+1) % N == 0) {
            q++;
            sinarg = sin(q*pi_2T);
            lamq = 4 * sinarg * sinarg;
        }
    }

    pv = imhof(lam, 0, err);
    if (!*err) {
        if (pv < 0) {
            pv = 0;
        }
        gretl_model_set_double(pmod, "dw_pval", pv);
    }


    gretl_matrix_free(lam);

    return pv;
}

/* Durbin-Watson statistic for the pooled or fixed effects model.

   See Bhargava, Franzini and Narendranathan, "Serial Correlation and
   the Fixed Effects Model", Review of Economic Studies 49, 1982,
   pp. 533-549, and also Baltagi and Wu, "Unequally Spaced Panel Data
   Regressions With AR(1) Disturbances", Econometric Theory, 15, 1999,
   pp. 814-823 for discussion of the unbalanced case.
*/

void panel_dwstat(Model * pmod, PanelmodT * pan)
{
    double dwnum = 0.0;
    double rnum = 0.0;
    double rden = 0.0;
    double ut, u1;
    int in_bounds = 1;
    int i, t, s;

    pmod.dw = pmod.rho = NADBL;

    if (pmod.ess <= 0.0) {
        return;
    }


    for (i=0; i<pan.nunits && in_bounds; i++) {
        int started = 0;

        if (pan.unit_obs[i] == 0) {
            continue;
        }
        for (t=0; t<pan.T; t++) {
            s = panel_index(i, t);
            if (s >= pmod.t2) {
                in_bounds = 0;
                break;
            }
            ut = pmod.uhat[s];
            if (!na(ut)) {
                if (t == 0) {
                    started = 1;
                    continue;
                }
                u1 = pmod.uhat[s-1];
                if (na(u1)) {
                    /* implicitly take u1 as 0 */
                    if (started) {
                        dwnum += ut * ut;
                    }
                } else {
                    dwnum += (ut - u1) * (ut - u1);
                    rnum += ut * u1;
                    rden += u1 * u1;
                }
                started = 1;
            }
        }
    }

    if (dwnum > 0.0) {
        pmod.dw = dwnum / pmod.ess;
    }
    if (na(pmod.rho) && rden > 0.0 && !na(rden)) {
        pmod.rho = rnum / rden;
    }
    if (!na(pmod.rho) && (pmod.rho <= -1.0 || pmod.rho >= 1.0)) {
        pmod.rho = NADBL;
    }

    if (pan.opt & Gretlopt.OPT_U) {
        /* make these stats available for random effects reporting */
        pan.dw = pmod.dw;
        pan.rho = pmod.rho;
    }
}

/* Allocate the arrays needed to perform the Hausman test,
   in its matrix formulation.
*/

int matrix_hausman_allocate(PanelmodT * pan)
{
    int k = pan.vlist[0] - 2;
    int err = 0;

    /* array to hold differences between coefficient estimates */
    pan.bdiff = gretl_vector_alloc(k);
    if (pan.bdiff == null) {
        err = GretlError.alloc;
    } else {
        /* array to hold covariance matrix */
        pan.Sigma = gretl_matrix_alloc(k, k);
        if (pan.Sigma == null) {
            gretl_matrix_free(pan.bdiff);
            pan.bdiff = null;
            err = GretlError.alloc;
        }
    }

    if (!err) {
        pan.dfH = k;
    }

    return err;
}

int * real_varying_list(PanelmodT * pan)
{
    int * vlist = gretl_list_copy(pan.vlist);
    int i;

    if (vlist != null) {
        for (i=2; i<=vlist[0]; i++) {
            if (vlist[i] == 0) {
                gretl_list_delete_at_pos(vlist, i);
                break;
            }
        }
    }

    return vlist;
}

/* Construct a version of the dataset from which the group means are
   subtracted, for the "within" regression.  Nota bene: this auxiliary
   dataset is not necessarily of full length: missing observations
   are skipped.
*/

DatasetRaw * within_groups_dataset(const(DatasetRaw) * dset, PanelmodT * pan)
{
    DatasetRaw * wset = null;
    int * vlist = null;
    int i, j, vj, nv;
    int s, t, bigt;
    int err = 0;

    pan.balanced = 1;

    for (i=0; i<pan.nunits; i++) {
        if (pan.unit_obs[i] > 0) {
            if (pan.unit_obs[i] != pan.Tmax) {
                pan.balanced = 0;
            }
        }
    }

    if (pan.NT < dset.n) {
    /* allow for panel obs book-keeping */
        err = allocate_data_finders(pan, dset.n);
        if (err) {
            return null;
        }
    }

    vlist = real_varying_list(pan);
    if (vlist == null) {
        return null;
    }

    nv = pan.vlist[0];


    wset = create_auxiliary_dataset(nv, pan.NT, Gretlopt.OPT_NONE);
    if (wset == null) {
        free(vlist);
        return null;
    }

    for (j=1; j<=vlist[0]; j++) {
        double xbar, gxbar = 0.0;
        int allzero = 1;

        vj = vlist[j];
        gxbar = 0.0;
        s = 0;


        for (i=0; i<pan.nunits; i++) {
            int Ti = pan.unit_obs[i];
            int grpzero = 1;
            int got = 0;

            if (Ti == 0) {
                continue;
            }

            /* first pass: find the group mean */
            xbar = 0.0;
            for (t=0; t<pan.T; t++) {
                bigt = panel_index(i, t);
                if (!panel_missing(pan, bigt)) {
                    xbar += dset.Z[vj][bigt];
                }
            }

            gxbar += xbar;
            xbar /= Ti;

            /* second pass: calculate de-meaned values */
            got = 0;
            for (t=0; t<pan.T && got<Ti; t++) {
                bigt = panel_index(i, t);
                if (!panel_missing(pan, bigt)) {
                    wset.Z[j][s] = dset.Z[vj][bigt] - xbar;
                    if (wset.Z[j][s] != 0.0) {
                        grpzero = 0;
                    }
                    got++;
                    if (pan.small2big != null) {
                        pan.small2big[s] = bigt;
                        pan.big2small[bigt] = s;
                    }
                    s++;
                }
            }

            if (!grpzero) {
                allzero = 0;
            }
        } /* end loop over units */

        gxbar /= pan.NT;

        if (j == 1 && allzero) {
            /* the dependent variable is not time-varying */
            fprintf(stderr, "fixed effects: dependent var is time-invariant\n".ptr);
        }

        /* wZ = data - group mean + grand mean */
        for (s=0; s<pan.NT; s++) {
            wset.Z[j][s] += gxbar;
        }
    }

    free(vlist);

    return wset;
}

/* Construct a quasi-demeaned version of the dataset so we can apply
   least squares to estimate the random effects model.  This dataset
   is not necessarily of full length.  If we're implementing the
   Hausman test using the regression approach this dataset will also
   include "straight" de-meaned counterparts of the time-varying
   variables.
*/

DatasetRaw * random_effects_dataset(const(DatasetRaw) * dset, const(DatasetRaw) * gset, int * relist, int * hlist, PanelmodT * pan)
{
    DatasetRaw * rset;
    double xbar, theta_i;
    int hreg = (hlist != null);
    int v1 = relist[0];
    int v2 = 0;
    int i, j, k, k2, t;
    int vj, s, bigt, u;
    int err = 0;

    if (hreg) {
        /* Apparatus for regression version of Hausman test:
           note that we shouldn't include time dummies here
           since the de-meaned versions would be perfectly
           collinear with the quasi-demeaned ones.
        */
        int hmax = pan.vlist[0] - pan.ntdum;

        for (i=1; i<hmax; i++) {
            if (pan.vlist[i+1] != 0) {
                hlist[0] = ++v2;
                hlist[i-1] = v1 + i - 2;
            }
        }
    }

    if (pan.NT < dset.n) {
        err = allocate_data_finders(pan, dset.n);
        if (err) {
            return null;
        }
    }


    rset = create_auxiliary_dataset(v1 + v2, pan.NT, Gretlopt.OPT_NONE);
    if (rset == null) {
        return null;
    }

    /* build GLS regression list, and process varnames for regression
       version of Hausman test if wanted
    */

    k = 0;
    k2 = v1 - 1;
    for (j=1; j<=v1; j++) {
        vj = pan.pooled.list[j];
        if (vj == 0) {
            relist[j] = 0;
        } else {
            relist[j] = ++k;
            strcpy(rset.varname[k], dset.varname[vj]);
            if (hreg && k2 < rset.v - 1 && j > 1 && var_is_varying(pan, vj)) {
                /* hausman-related term */
                k2++;
                strcpy(rset.varname[k2], "_");
                strncat(rset.varname[k2], dset.varname[vj],
                        VNAMELEN - 2);
            }
        }
    }

    /* Now create the transformed variables: original data minus theta
       times the appropriate group mean (and in addition, original
       data minus the group mean for time-varying vars, if doing the
       Hausman test by the regression method).
    */

    s = u = 0;
    pan.theta_bar = 0.0;

    for (i=0; i<pan.nunits; i++) {
        int Ti = pan.unit_obs[i];
        int got = 0;

        if (Ti == 0) {
            continue;
        }

        if (Ti != pan.Tmax) {
            theta_i = 1.0 - sqrt(pan.s2e / (Ti * pan.s2v + pan.s2e));
        } else {
            theta_i = pan.theta;
        }

        pan.theta_bar += theta_i;

        for (t=0; t<pan.T && got<Ti; t++) {
            bigt = panel_index(i, t);
            k = 0;
            k2 = v1 - 1;
            if (!panel_missing(pan, bigt)) {
                for (j=0; j<v1; j++) {
                    vj = pan.pooled.list[j+1];
                    if (vj == 0) {
                        rset.Z[0][s] -= theta_i;
                    } else {
                        k++;
                        xbar = (k < gset.v)? gset.Z[k][u] : 1.0 / pan.Tmax;
                        rset.Z[k][s] = dset.Z[vj][bigt] - theta_i * xbar;
                        if (hreg && k2 < rset.v - 1 && var_is_varying(pan, vj)) {
                            /* hausman-related term */
                            rset.Z[++k2][s] = dset.Z[vj][bigt] - xbar;
                        }
                    }
                }
                got++;
                if (pan.small2big != null) {
                    pan.small2big[s] = bigt;
                    pan.big2small[bigt] = s;
                }
                s++;
            }
        }
        u++;
    }

    return rset;
}

/* Construct a "mini-dataset" containing the group means, in
   order to run the group-means or "between" regression.
*/

DatasetRaw * group_means_dataset(PanelmodT * pan, const(DatasetRaw) * dset)
{
    DatasetRaw * gset;
    double x;
    int gn = pan.effn;
    int gv = pan.pooled.list[0];
    int i, j, k;
    int s, t, bigt;

    if (pan.balanced && pan.ntdum > 0) {
        gv -= pan.ntdum;
    }


    gset = create_auxiliary_dataset(gv, gn, Gretlopt.OPT_NONE);
    if (gset == null) {
        return null;
    }

    k = 1;
    for (j=1; j<=gv; j++) {
        int vj = pan.pooled.list[j];

        if (vj == 0) {
            continue;
        }

        if (pan.opt & Gretlopt.OPT_B) {
            /* will save the "between" model: so name the variables */
            strcpy(gset.varname[k], dset.varname[vj]);
        }

        s = 0;
        for (i=0; i<pan.nunits; i++) {
            int Ti = pan.unit_obs[i];

            if (Ti == 0) {
                continue;
            }

            x = 0.0;
            for (t=0; t<pan.T; t++) {
                bigt = panel_index(i, t);
                if (!panel_missing(pan, bigt)) {
                    x += dset.Z[vj][bigt];
                }
            }
            gset.Z[k][s++] = x / Ti;
        }
        k++;
    }


    return gset;
}

/* spruce up the between model and attach it to pan */

int save_between_model(Model * pmod, const(int) * blist, DatasetRaw * gset, PanelmodT * pan)
{
    GretlMatrixRaw * uh, yh;
    int i, err = 0;

    pmod.ci = GretlCmdIndex.panel;
    pmod.opt |= Gretlopt.OPT_B;
    pmod.dw = NADBL;
    gretl_model_add_panel_varnames(pmod, gset, null);

    uh = gretl_column_vector_alloc(pmod.nobs);
    yh = gretl_column_vector_alloc(pmod.nobs);

    if (uh == null || yh == null) {
        err = GretlError.alloc;
    } else {
        for (i=0; i<pmod.nobs; i++) {
            uh.val[i] = pmod.uhat[i];
            yh.val[i] = pmod.yhat[i];
        }
        gretl_model_set_matrix_as_data(pmod, "uhat", uh);
        gretl_model_set_matrix_as_data(pmod, "yhat", yh);
    }


    *pan.realmod = *pmod;

    return err;
}

/* Compute @ubPub as the Ti-weighted sum of the squared
   residuals from the Between model, as per Stata (but
   in disagreement with Baltagi and Chang, 1994).
*/

int alt_compute_ubPub(PanelmodT * pan, Model * bmod)
{
    int i, Ti, t = 0;

    pan.ubPub = 0.0;

    for (i=0; i<pan.nunits; i++) {
        Ti = pan.unit_obs[i];
        if (Ti > 0) {
            pan.ubPub += Ti * bmod.uhat[t] * bmod.uhat[t];
            t++;
        }
    }

    return 0;
}

void adjust_gset_data(PanelmodT * pan, DatasetRaw * gset, int step)
{
    int i, j, Ti, t = 0;
    double adj;

    for (i=0; i<pan.nunits; i++) {
        Ti = pan.unit_obs[i];
        if (Ti > 0) {
            adj = step == 0 ? sqrt(Ti) : 1.0/sqrt(Ti);
            for (j=0; j<gset.v; j++) {
                gset.Z[j][t] *= adj;
            }
            t++;
        }
    }
}

/* Compute @ubPub as the sum of squared residuals from a
   Ti-weighted Between regression, as per Baltagi and Chang,
   1994, and also Baltagi, 2013.
*/

int compute_ubPub(PanelmodT * pan, Model * bmod, int * blist, DatasetRaw * gset)
{
    int err;

    /* multiply all data by sqrt(Ti) */
    adjust_gset_data(pan, gset, 0);
    clear_model(bmod);
    *bmod = lsq(blist, gset, GretlCmdIndex.ols, Gretlopt.OPT_A);
    err = bmod.errcode;
    if (!err) {
        pan.ubPub = bmod.ess;
    }
    /* put the original data back */
    adjust_gset_data(pan, gset, 1);

    return err;
}

/* calculate the group means or "between" regression and its error
   variance */

int between_variance(PanelmodT * pan, DatasetRaw * gset)
{
    Gretlopt bopt;
    Model bmod;
    int * blist;
    int i, j, k;
    int err = 0;

    blist = gretl_list_new(gset.v);
    if (blist == null) {
        return GretlError.alloc;
    }

    j = k = 1;
    for (i=1; i<=gset.v; i++) {
        if (pan.pooled.list[i] == 0) {
            blist[k++] = 0;
        } else {
            blist[k++] = j++;
        }
    }

    bopt = (pan.opt & Gretlopt.OPT_B)? Gretlopt.OPT_NONE : Gretlopt.OPT_A;
    bmod = lsq(blist, gset, GretlCmdIndex.ols, bopt);

    if (bmod.errcode == 0) {
        pan.s2b = bmod.ess / (bmod.nobs - bmod.ncoeff);
    } else {
        err = bmod.errcode;
    }


    if (!err && (pan.opt & Gretlopt.OPT_B)) {
        err = save_between_model(&bmod, blist, gset, pan);
    } else {
        if (!err && !pan.balanced && (pan.opt & Gretlopt.OPT_U) &&
            (pan.opt & Gretlopt.OPT_X) && !(pan.opt & Gretlopt.OPT_E)) {
            /* Prepare for the Baltagi-Chang take on Swamy-Arora
               in the case of an unbalanced panel
            */
            if (stata_sa) {
                err = alt_compute_ubPub(pan, &bmod);
            } else {
                err = compute_ubPub(pan, &bmod, blist, gset);
            }
        }
        clear_model(&bmod);
    }

    free(blist);

    return err;
}

/* op is VCV_SUBTRACT for the random effects model.  With the fixed
   effects model we don't have to worry about excluding vcv entries
   for time-invariant variables, since none of these are included.
*/

int vcv_skip(const(Model) * pmod, int i, const(PanelmodT) * pan, int op)
{
    int skip = 0;

    if (pmod.list[i+2] == 0) {
        /* always skip the constant */
        skip = 1;
    } else if (op == VcvOps.subtract && !pan.varying[i]) {
        /* random effects, time-invariant var */
        skip = 1;
    }

    return skip;
}

/* Fill out the covariance matrix for use with the Hausman test:
   the entries that get transcribed here are only those for
   slopes with respect to time-varying variables.
*/

void vcv_slopes(PanelmodT * pan, const(Model) * pmod, int op)
{
    int idx, i, j;
    int mj, mi = 0;
    int sj, si = 0;
    double x;

    for (i=0; i<pan.dfH; i++) {
        if (vcv_skip(pmod, mi, pan, op)) {
            i--;
            mi++;
            continue;
        }
        mj = mi;
        sj = si;
        for (j=i; j<pan.dfH; j++) {
            if (vcv_skip(pmod, mj, pan, op)) {
                j--;
                mj++;
                continue;
            }

            idx = ijton(mi, mj, pmod.ncoeff);

            if (op == VcvOps.subtract) {
                x = gretl_matrix_get(pan.Sigma, si, sj);
                x -= pmod.vcv[idx];
            } else {
                x = pmod.vcv[idx];
            }

            gretl_matrix_set(pan.Sigma, si, sj, x);
            if (si != sj) {
                gretl_matrix_set(pan.Sigma, sj, si, x);
            }

            mj++;
            sj++;
        }
        mi++;
        si++;
    }
}

/* calculate Hausman test statistic, matrix diff style */

int bXb(PanelmodT * pan)
{
    int err;

    err = gretl_invert_symmetric_matrix(pan.Sigma);

    if (!err) {
        pan.H = gretl_scalar_qform(pan.bdiff, pan.Sigma, &err);
    }

    if (err) {
        pan.H = NADBL;
    }

    return err;
}

void fixed_effects_df_correction(Model * pmod, int k)
{
    double dfcorr = sqrt(cast(double) pmod.dfd / (pmod.dfd - k));
    int i;

    pmod.dfd -= k;
    pmod.dfn += k;

    pmod.sigma *= dfcorr;

    for (i=0; i<pmod.ncoeff; i++) {
        pmod.sderr[i] *= dfcorr;
    }
}

/* used for printing fixed- or random-effects estimates
   in the context of the "panel diagnostics" routine
*/

void simple_print_panel_model(Model * pmod, DatasetRaw * dset, const(int) * xlist, Prn * prn)
{
    int * savelist = gretl_list_copy(pmod.list);
    int i;

    for (i=0; i<pmod.ncoeff; i++) {
        pmod.list[i+2] = xlist[i+2];
    }
    printmodel(pmod, dset, Gretlopt.OPT_P, prn);
    free(pmod.list);
    pmod.list = savelist;
}

void print_re_results(PanelmodT * pan, Model * pmod, DatasetRaw * dset, Prn * prn)
{
    pputs(prn, _("Variance estimators:").ptr);
    pputc(prn, '\n');
    pprintf(prn, _(" between = %g").ptr, pan.s2v);
    pputc(prn, '\n');
    pprintf(prn, _(" within = %g").ptr, pan.s2e);
    pputc(prn, '\n');

    if (pan.balanced || pan.s2v == 0) {
        pprintf(prn, _("theta used for quasi-demeaning = %g").ptr, pan.theta);
        pputc(prn, '\n');
    } else {
        pputs(prn, _("Panel is unbalanced: theta varies across units").ptr);
        pputc(prn, '\n');
    }
    pputc(prn, '\n');

    pputs(prn, _("Random effects estimator\nallows for a unit-specific component to the error term\n").ptr);
    pputc(prn, '\n');

    simple_print_panel_model(pmod, dset, pan.pooled.list, prn);
}

int print_fe_results(PanelmodT * pan, Model * pmod, DatasetRaw * dset, Prn * prn)
{
    int dfn = pan.effn - 1;

    pputs(prn, _("Fixed effects estimator\nallows for differing intercepts by cross-sectional unit\n").ptr);
    pputc(prn, '\n');

    simple_print_panel_model(pmod, dset, pan.vlist, prn);

    pprintf(prn, _("Residual variance: %g/(%d - %d) = %g\n").ptr,
            pmod.ess, pmod.nobs, pan.vlist[0] - 1 + dfn, pan.s2e);
    pputc(prn, '\n');

    if (!na(pan.Ffe)) {
        pprintf(prn, _("Joint significance of differing group means:\n").ptr);
        pprintf(prn, " F(%d, %g) = %g %s %g\n".ptr, pan.Fdfn, pan.Fdfd, pan.Ffe,
                _("with p-value").ptr, snedecor_cdf_comp(pan.Fdfn, pan.Fdfd, pan.Ffe));
        pputs(prn, _("(A low p-value counts against the null hypothesis that the pooled OLS model\nis adequate, in favor of the fixed effects alternative.)\n\n").ptr);
    } else {
        pputc(prn, '\n');
    }

    return 0;
}

int time_dummies_wald_test(PanelmodT * pan, Model * pmod)
{
    GretlMatrixRaw * vcv = null;
    GretlVector * b = null;
    double x;
    int i, j, k, bigk;
    int di, dj;
    int err;

    if (pan.ntdum == 0) {
        return 0;
    }

    k = pan.ntdum;
    bigk = pmod.ncoeff;

    if (pmod.vcv == null) {
        err = makevcv(pmod, pmod.sigma);
        if (err) {
            return err;
        }
    }

    b = gretl_column_vector_alloc(k);
    vcv = gretl_matrix_alloc(k, k);
    if (b == null || vcv == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    di = bigk - k;
    for (i=0; i<k; i++) {
        b.val[i] = pmod.coeff[di++];
    }

    di = bigk - k;
    for (i=0; i<k; i++) {
        dj = bigk - k;
        for (j=0; j<=i; j++) {
            x = pmod.vcv[ijton(di, dj++, bigk)];
            gretl_matrix_set(vcv, i, j, x);
            gretl_matrix_set(vcv, j, i, x);
        }
        di++;
    }

    err = gretl_invert_symmetric_matrix(vcv);
    if (err) {
        goto bailout;
    }

    x = gretl_scalar_qform(b, vcv, &err);
    if (err) {
        fprintf(stderr, "Failed to compute test statistic\n".ptr);
        goto bailout;
    }

    if (!err) {
        ModelTest *test = model_test_new(ModelTestType.panel_timedum);

        if (test != null) {
            model_test_set_teststat(test, TestStats.wald_chisq);
            model_test_set_dfn(test, k);
            model_test_set_value(test, x);
            model_test_set_pvalue(test, chisq_cdf_comp(k, x));
            maybe_add_test_to_model(pmod, test);
        }
    }

 bailout:

    gretl_matrix_free(vcv);
    gretl_vector_free(b);

    return err;
}

void save_fixed_effects_F(PanelmodT * pan, Model * wmod)
{
    int robust = (pan.opt & Gretlopt.OPT_R);
    ModelTest *test;

    if (na(pan.Ffe)) {
        return;
    }

    test = model_test_new(robust ? ModelTestType.panel_welch :
                          ModelTestType.panel_f);

    if (test != null) {
        model_test_set_teststat(test, robust ? TestStats.wf : TestStats.f);
        model_test_set_dfn(test, pan.Fdfn);
        model_test_set_dfd(test, pan.Fdfd);
        model_test_set_value(test, pan.Ffe);
        model_test_set_pvalue(test, snedecor_cdf_comp(pan.Fdfn, pan.Fdfd, pan.Ffe));
        maybe_add_test_to_model(wmod, test);
    }
}

/* "regular" = not robust: sums-of-squares based joint test on
   the fixed effects */

void regular_fixed_effects_F(PanelmodT * pan, Model * wmod)
{
    int k_pooled = pan.pooled.list[0];
    int k_fe = pan.vlist[0];

    pan.Fdfn = pan.effn - 1;
    pan.Fdfd = wmod.dfd;

    if (k_pooled > k_fe) {
        pan.Fdfn -= k_pooled - k_fe;
        if (pan.Fdfn <= 0) {
            pan.Ffe = NADBL;
            return;
        }
    }

    pan.Ffe = (pan.pooled.ess - wmod.ess) * pan.Fdfd /
        (wmod.ess * pan.Fdfn);
    if (na(pan.Ffe)) {
        pan.Ffe = NADBL;
    } else if (pan.Ffe < 0) {
        pan.Ffe = 0;
    }
}

int fix_panelmod_list(Model * targ, PanelmodT * pan)
{
    int i;


    free(targ.list);
    targ.list = gretl_list_copy(pan.pooled.list);
    if (targ.list == null) {
        return GretlError.alloc;
    }

    if (pan.opt & Gretlopt.OPT_F) {
        /* fixed effects: remove any non-varying variables */
        for (i=2; i<=targ.list[0]; i++) {
            if (!in_gretl_list(pan.vlist, targ.list[i])) {
                gretl_list_delete_at_pos(targ.list, i--);
            }
        }
    }


    return 0;
}

/* Compute F or chi-square test for the regular regressors, skipping
   @nskip trailing coefficients (e.g. time dummies).
*/

double panel_overall_test(Model * pmod, PanelmodT * pan, int nskip, Gretlopt opt)
{
    double test = NADBL;
    char * mask = null;
    int i, nm;
    int err = 0;

    nm = pmod.ncoeff - pmod.ifc - nskip;

    if (nm <= 1) {
        return test;
    }

    mask = cast(typeof(mask)) calloc(pmod.ncoeff, 1);
    for (i=0; i<nm; i++) {
    mask[i+pmod.ifc] = 1;
    }

    test = wald_omit_test(null, mask, pmod, opt, &err);
    free(mask);

    return test;
}

/* Correct various model statistics, in the case where we estimated
   the fixed effects or "within" model on an auxiliary dataset
   from which the group means were subtracted.
*/

int fix_within_stats(Model * fmod, PanelmodT * pan)
{
    double wrsq, wfstt = NADBL;
    int wdfn, nc = fmod.ncoeff;
    int err = 0;

    err = fix_panelmod_list(fmod, pan);
    if (err) {
        return err;
    }

    fmod.ybar = pan.pooled.ybar;
    fmod.sdy = pan.pooled.sdy;
    fmod.tss = pan.pooled.tss;
    fmod.ifc = 1;

    wrsq = fmod.rsq;
    if (wrsq < 0.0) {
        wrsq = 0.0;
    } else if (na(wrsq)) {
        wrsq = NADBL;
    }

    if (pan.ntdum > 0) {
        wdfn = fmod.ncoeff - 1 - pan.ntdum;
        if (wdfn > 0) {
            wfstt = panel_overall_test(fmod, pan, pan.ntdum,
                                       Gretlopt.OPT_F);
        }
    } else {
        wdfn = fmod.ncoeff - 1;
        if (wdfn > 0) {
            if (pan.opt & Gretlopt.OPT_R) {
                wfstt = panel_overall_test(fmod, pan, 0, Gretlopt.OPT_F);
            } else {
                wfstt = (wrsq / (1.0 - wrsq)) * (cast(double) fmod.dfd / wdfn);
            }
        }
    }

    if (!na(wfstt) && wfstt >= 0.0) {
        ModelTest *test = model_test_new(ModelTestType.within_f);
        int nclu, wdfd = fmod.dfd;

        if (pan.opt & Gretlopt.OPT_C) {
            nclu = gretl_model_get_int(fmod, "n_clusters");
            fmod.dfd = wdfd = nclu - 1;
        } else if (pan.opt & Gretlopt.OPT_R) {
        nclu = gretl_model_get_int(fmod, "DKT");
        if (nclu > 0) {
        /* Driscoll-Kraay */
        fmod.dfd = wdfd = nclu - 1;
        } else {
        fmod.dfd = wdfd = pan.effn - 1;
        }
        }

        if (test != null) {
            model_test_set_teststat(test, TestStats.f);
            model_test_set_dfn(test, wdfn);
            model_test_set_dfd(test, wdfd);
            model_test_set_value(test, wfstt);
            model_test_set_pvalue(test, snedecor_cdf_comp(wdfn, wdfd, wfstt));
            maybe_add_test_to_model(fmod, test);
        }
    }

    /* note: this member is being borrowed for the "Within R-squared" */
    fmod.adjrsq = wrsq;

    /* LSDV-based statistics (FIXME: can we do this in
       the --robust case?)
    */
    if (pan.opt & Gretlopt.OPT_R) {
        fmod.rsq = 1.0 - (fmod.ess / fmod.tss);
        fmod.fstt = NADBL;
    } else {
        fmod.rsq = 1.0 - (fmod.ess / fmod.tss);
        if (fmod.rsq < 0.0) {
            fmod.rsq = 0.0;
        } else {
            fmod.fstt = (fmod.rsq / (1.0 - fmod.rsq)) *
                (cast(double) fmod.dfd / fmod.dfn);
        }
    }

    fmod.ncoeff = fmod.dfn + 1; /* number of params estimated */
    ls_criteria(fmod);
    fmod.ncoeff = nc;

    return err;
}

/* Fixed-effects model: add the per-unit intercept estimates
   to the model in case the user wants to retrieve them.
   Random-effects model: add estimate of the individual effects.

   By this point the model -- even if it been estimated on a short
   dataset -- should have a full-length residual series.
*/

int panel_model_add_ahat(Model * pmod, const(DatasetRaw) * dset, PanelmodT * pan)
{
    double * ahat = null;
    const(double) * x;
    int i, j, t, bigt;
    int n, err = 0;

    n = pmod.full_n;

    ahat = cast(typeof(ahat)) malloc(n * (*ahat).sizeof);
    if (ahat == null) {
        return GretlError.alloc;
    }

    for (t=0; t<n; t++) {
        ahat[t] = NADBL;
    }

    if (pan.opt & Gretlopt.OPT_F) {
        /* fixed effects */
        for (i=0; i<pan.nunits; i++) {
            if (pan.unit_obs[i] > 0) {
                double a = 0.0;

                /* a = y - Xb, where the 'b' is based on de-meaned data */

                for (t=0; t<pan.T; t++) {
                    bigt = panel_index(i, t);
                    if (!na(pmod.uhat[bigt])) {
                        a += dset.Z[pmod.list[1]][bigt];
                        for (j=1; j<pmod.ncoeff; j++) {
                            x = dset.Z[pmod.list[j+2]];
                            a -= pmod.coeff[j] * x[bigt];
                        }
                    }
                }

                a /= pan.unit_obs[i];

                for (t=0; t<pan.T; t++) {
                    bigt = panel_index(i, t);
                    if (!na(pmod.uhat[bigt])) {
                        ahat[bigt] = a;
                    }
                }
            }
        }
    } else {
        /* random effects */
        double uhbar, frac = 0;

        if (pan.balanced) {
            frac = 1.0 - pan.theta;
            frac = 1.0 - frac * frac;
        }

        for (i=0; i<pan.nunits; i++) {
            if (pan.unit_obs[i] > 0) {
                /* get mean residual */
                uhbar = 0.0;
                for (t=0; t<pan.T; t++) {
                    bigt = panel_index(i, t);
                    if (!na(pmod.uhat[bigt])) {
                        uhbar += pmod.uhat[bigt];
                    }
                }
                uhbar /= pan.unit_obs[i];
                /* ahat = frac * uhbar */
                if (!pan.balanced) {
                    frac = pan.s2v / (pan.s2v + pan.s2e / pan.unit_obs[i]);
                }
                for (t=0; t<pan.T; t++) {
                    bigt = panel_index(i, t);
                    if (!na(pmod.uhat[bigt])) {
                        ahat[bigt] = frac * uhbar;
                    }
                }
            }
        }
    }

    err = gretl_model_set_data(pmod, "ahat", ahat,
                               GretlType.double_array,
                               n * (*ahat).sizeof);

    return err;
}

/* If we're estimating the fixed effects model with the --robust
   flag, we should do a robust version of the joint test on
   the fixed effects. Here we use the algorithm of B. L. Welch,
   "On the Comparison of Several Mean Values: An Alternative
   Approach" (Biometrika 38, 1951, pp. 330-336). The variable
   we're testing for difference of means (by individual) is the
   residual from pooled OLS.
*/

int robust_fixed_effects_F(PanelmodT * pan)
{
    Model * pmod = pan.pooled;
    double * u, w, h, xbar;
    double muhat = 0.0;
    double x, s2, W = 0.0;
    double A, B, sum_h;
    int k = pan.effn;
    int i, j, t, s;
    int Ti, bigt;
    int err = 0;

    u = cast(typeof(u)) malloc(pan.Tmax * (*u).sizeof);
    w = cast(typeof(w)) malloc(3 * k * (*w).sizeof);

    if (u == null || w == null) {
        free(u);
        free(w);
        return GretlError.alloc;
    }

    h = w + k;
    xbar = h + k;

    j = 0;
    for (i=0; i<pan.nunits; i++) {
        Ti = pan.unit_obs[i];
        if (Ti > 1) {
            s = 0;
            for (t=0; t<pan.T; t++) {
                bigt = panel_index(i, t);
                if (!panel_missing(pan, bigt)) {
                    u[s++] = pmod.uhat[bigt];
                }
            }
            xbar[j] = gretl_mean(0, s-1, u);
            s2 = gretl_variance(0, s-1, u);
            w[j] = Ti / s2;
            W += w[j];
            muhat += w[j] * xbar[j];
            j++;
        }
    }

    muhat /= W;
    A = sum_h = 0.0;

    j = 0;
    for (i=0; i<pan.nunits; i++) {
        Ti = pan.unit_obs[i];
        if (Ti > 1) {
            x = 1 - w[j]/W;
            h[j] = (x * x) / (Ti - 1);
            sum_h += h[j];
            x = xbar[j] - muhat;
            A += w[j] * x * x;
            j++;
        }
    }

    A /= (k - 1);
    B = 1.0 + (2.0*(k-2.0)/(k * k - 1.0)) * sum_h;

    pan.Ffe = A / B;
    pan.Fdfn = k - 1;
    pan.Fdfd = (k * k - 1.0)/(3 * sum_h);

    free(u);
    free(w);

    return err;
}

int * real_FE_list(PanelmodT * pan)
{
    int * list = gretl_list_copy(pan.pooled.list);
    int i;

    if (list != null) {
        /* purge any non-time-varying variables */
        for (i=2; i<=list[0]; i++) {
            if (!in_gretl_list(pan.vlist, list[i])) {
                gretl_list_delete_at_pos(list, i--);
            }
        }
    }

    return list;
}

/* For testing purposes: retrieve user-specified values for
   s2v and s2e. These would typically be known good values
   in the context of a simulation.
*/

int read_true_variances(PanelmodT * pan)
{
    GretlMatrixRaw * m;
    int err = 0;

    m = gretl_matrix_read_from_file(glsmat.ptr, 0, &err);
    if (m == null) {
        fprintf(stderr, "read_true_variances: no matrix!\n".ptr);
        if (!err) {
            err = GretlError.data;
        }
    } else {
        pan.s2v = m.val[0];
        pan.s2e = m.val[1];
    }

    return err;
}

/* computation of $\hat{\sigma}^2_v$ a la Nerlove, if wanted */

int nerlove_s2v(Model * pmod, const(DatasetRaw) * dset, PanelmodT * pan)
{
    double a;
    double amean;
    double * ahat;
    double wmean;
    double * wi = null;
    const(double) * x;
    int i, j, t, k, bigt;
    int * list = null;
    int err = 0;

    ahat = cast(typeof(ahat)) malloc(pan.effn * (*ahat).sizeof);
    if (ahat == null) {
        return GretlError.alloc;
    }

    if (pan.opt & Gretlopt.OPT_X) {
        /* --unbalanced */
        wi = cast(typeof(wi)) malloc(pan.effn * (*wi).sizeof);
        if (wi == null) {
            free(ahat);
            return GretlError.alloc;
        }
    }

    list = real_FE_list(pan);
    if (list == null) {
        free(ahat);
        free(wi);
        return GretlError.alloc;
    }

    wmean = amean = 0.0;
    k = 0;

    for (i=0; i<pan.nunits; i++) {
        int Ti = pan.unit_obs[i];

        if (Ti > 0) {
            a = 0.0;
            for (t=0; t<pan.T; t++) {
                bigt = panel_index(i, t);
                if (!na(pan.pooled.uhat[bigt])) {
                    a += dset.Z[list[1]][bigt];
                    for (j=1; j<pmod.ncoeff; j++) {
                        x = dset.Z[list[j+2]];
                        a -= pmod.coeff[j] * x[bigt];
                    }
                }
            }
            a /= Ti;
            ahat[k] = a;
            amean += a;
            if (wi != null) {
                wi[k] = Ti / cast(double) pmod.nobs;
                wmean += wi[k] * a;
            }
            k++;
        }
    }

    pan.s2v = 0.0;

    if (wi != null) {
        for (i=0; i<pan.effn; i++) {
            pan.s2v += wi[i] * (ahat[i] - wmean) * (ahat[i] - wmean);
        }
        pan.s2v /= (pan.effn - 1.0) / cast(double) pan.effn;
        free(wi);
    } else {
        amean /= pan.effn;
        for (i=0; i<pan.effn; i++) {
            pan.s2v += (ahat[i] - amean) * (ahat[i] - amean);
        }
        pan.s2v /= pan.effn - 1;
    }

    free(ahat);
    free(list);

    return err;
}

/* robust RE: we need an extra residuals array */

int re_hatvars_prep(PanelmodT * pan)
{
    int t, n = pan.pooled.full_n;

    pan.re_uhat = cast(typeof(pan.re_uhat)) malloc(n * (*pan.re_uhat).sizeof);

    if (pan.re_uhat == null) {
        return GretlError.alloc;
    } else {
        for (t=0; t<n; t++) {
            pan.re_uhat[t] = NADBL;
        }
        return 0;
    }
}

/* When calculating DW using the fixed-effects residuals,
   in the context of estimation of the random-effects
   specification, we need to expand up the residual series
   temporarily.
*/

int expand_fe_uhat(Model * femod, PanelmodT * pan)
{
    double * uhat = null;
    int NT = pan.pooled.full_n;
    int i, s, t;

    uhat = cast(typeof(uhat)) malloc(NT * (*uhat).sizeof);
    if (uhat == null) {
        return GretlError.alloc;
    }

    for (t=0; t<NT; t++) {
        uhat[t] = NADBL;
    }

    s = 0;
    for (i=0; i<pan.nunits; i++) {
        int ti, Ti = pan.unit_obs[i];

        if (Ti == 0) {
            continue;
        }
        for (ti=0; ti<Ti; ti++) {
            t = big_index(pan, s);
            uhat[t] = femod.uhat[s];
            if (s == 0) {
                femod.t1 = t;
            } else if (s == femod.nobs - 1) {
                femod.t2 = t;
            }
            s++;
        }
    }

    /* replace the uhat array on @femod */
    free(femod.uhat);
    femod.uhat = uhat;

    return 0;
}

/* Fix uhat and yhat in two cases.

   (a) When we estimated fixed effects using a de-meaned dataset we
   need to ensure that the uhat and yhat values get written to the
   right observation slots in relation to the full dataset, and also
   that the yhat values get corrected, putting the means back in.

   (b) When estimating the random effects model we need to compute
   residuals based on the untransformed data (and again, place them
   correctly in relation to the full dataset). However, if we're
   going to produce robust standard errors we also need to preserve
   the GLS residuals.

   The placement issue arises because the special datasets used in
   these cases are not necessarily of full length, since they are
   purged of missing observations.
*/

int fix_panel_hatvars(Model * pmod, PanelmodT * pan, const(double) ** Z)
{
    const(double) * y = null;
    double * yhat = pan.pooled.yhat;
    double * uhat = null;
    int n = pan.pooled.full_n;
    double yht, SSR = 0.0;
    int i, j, s, t;
    int err = 0;

    if (yhat == null) {
        fprintf(stderr, "fix_panel_hatvars: pan->pooled->yhat is NULL\n".ptr);
        return GretlError.data;
    }

    y = Z[pan.pooled.list[1]];

    uhat = cast(typeof(uhat)) malloc(n * (*uhat).sizeof);
    if (uhat == null) {
        return GretlError.alloc;
    }

    if (pan.opt & Gretlopt.OPT_U) {
        /* random effects */
        if (pan.opt & Gretlopt.OPT_R) {
            err = re_hatvars_prep(pan);
        }
        if (err) {
            return err;
        }
    }

    for (t=0; t<n; t++) {
        uhat[t] = NADBL;
    }

    s = 0;
    for (i=0; i<pan.nunits; i++) {
        int ti, Ti = pan.unit_obs[i];

        if (Ti == 0) {
            continue;
        }

        for (ti=0; ti<Ti; ti++) {
            t = big_index(pan, s);
            if (pan.opt & Gretlopt.OPT_U) {
                /* random effects */
                yht = 0.0;
                for (j=0; j<pmod.ncoeff; j++) {
                    yht += pmod.coeff[j] * Z[pan.pooled.list[j+2]][t];
                }
                yhat[t] = yht;
                uhat[t] = y[t] - yht;
                SSR += uhat[t] * uhat[t];
                if (pan.re_uhat != null) {
                    /* store both the "fixed" and the GLS residuals */
                    pan.re_uhat[t] = uhat[t];
                    uhat[t] = pmod.uhat[s];
                }
            } else {
                /* fixed effects */
                uhat[t] = pmod.uhat[s];
                yhat[t] = y[t] - uhat[t];
            }
            if (s == 0) {
                pmod.t1 = t;
            } else if (s == pmod.nobs - 1) {
                pmod.t2 = t;
            }
            s++;
        }
    }

    if (pan.opt & Gretlopt.OPT_U) {
        double r;

        /* Defer rewriting of pmod->ess and pmod->sigma to avoid
           screwing up the covariance matrix estimator and Hausman
           test.
        */
        gretl_model_set_double(pmod, "fixed_SSR", SSR);
        r = gretl_corr(pmod.t1, pmod.t2, y, yhat, null);
        if (!na(r)) {
            gretl_model_set_double(pmod, "corr-rsq", r * r);
        }
    }

    pmod.full_n = n;

    /* replace the uhat and yhat arrays on @pmod */
    free(pmod.uhat);
    pmod.uhat = uhat;
    free(pmod.yhat);
    pmod.yhat = yhat;

    /* NULLify stolen pointer */
    pan.pooled.yhat = null;

    return err;
}

int hausman_move_uhat(Model * pmod, PanelmodT * pan)
{
    int n = pan.pooled.full_n;
    double * uhat;
    int i, s, t, ti;

    uhat = cast(typeof(uhat)) malloc(n * (*uhat).sizeof);
    if (uhat == null) {
        return GretlError.alloc;
    }

    for (t=0; t<n; t++) {
        uhat[t] = NADBL;
    }

    s = 0;
    for (i=0; i<pan.nunits; i++) {
        int Ti = pan.unit_obs[i];

        if (Ti > 0) {
            for (ti=0; ti<Ti; ti++) {
                t = big_index(pan, s);
                uhat[t] = pmod.uhat[s];
                s++;
            }
        }
    }

    free(pmod.uhat);
    pmod.uhat = uhat;

    return 0;
}


/* Estimate the fixed-effects model using a parallel dataset
   with the group means subtracted from all variables.
*/

Model fixed_effects_model(PanelmodT * pan, DatasetRaw * dset, Prn * prn)
{
    Model femod;
    Gretlopt lsqopt = Gretlopt.OPT_A | Gretlopt.OPT_X;
    DatasetRaw * wset = null;
    int * felist = null;
    int save_qr = 0;
    int i;


    gretl_model_init(&femod, dset);

    felist = gretl_list_new(pan.vlist[0]);
    if (felist == null) {
        femod.errcode = GretlError.alloc;
        return femod;
    }

    wset = within_groups_dataset(dset, pan);
    if (wset == null) {
        free(felist);
        femod.errcode = GretlError.alloc;
        return femod;
    }

    felist[1] = 1;
    felist[2] = 0;
    for (i=3; i<=felist[0]; i++) {
        felist[i] = i - 1;
    }

    if (pan.opt & Gretlopt.OPT_F) {
        /* suppress auto-removal of collinear terms */
        lsqopt |= Gretlopt.OPT_Z;
    }
    if (pan.opt & Gretlopt.OPT_N) {
        /* suppress df correction */
        lsqopt |= Gretlopt.OPT_N;
    }

    save_qr = libset_get_int(SetKey.use_qr);
    libset_set_int(SetKey.use_qr, 1);

    femod = lsq(felist, wset, GretlCmdIndex.ols, lsqopt);

    libset_set_int(SetKey.use_qr, save_qr);

    if (femod.errcode) {
        fprintf(stderr, "femod.errcode = %d\n".ptr, femod.errcode);
    } else if ((pan.opt & Gretlopt.OPT_F) && femod.list[0] < felist[0]) {
        femod.errcode = GretlError.singular;
    } else {
        if (!(pan.opt & Gretlopt.OPT_N)) {
            /* we estimated a bunch of group means, and have to
               subtract degrees of freedom */
            fixed_effects_df_correction(&femod, pan.effn - 1);
        }
        if (pan.opt & Gretlopt.OPT_F) {
            /* estimating the FE model in its own right */
            if (pan.opt & Gretlopt.OPT_R) {
                /* we have to do this before the pooled residual
                   array is "stolen" for the fixed-effects model
                */
                robust_fixed_effects_F(pan);
            }
            fix_panel_hatvars(&femod, pan, cast(const(double) **) dset.Z);
            if (pan.opt & Gretlopt.OPT_R) {
                femod.errcode = panel_robust_vcv(&femod, pan, wset, dset);
            } else {
                femod_regular_vcv(&femod);
            }
        } else if (pan.opt & Gretlopt.OPT_E) {
            if (IGLS) {
                read_true_variances(pan);
            } else {
                femod.errcode = nerlove_s2v(&femod, dset, pan);
            }
        }
    }

    destroy_dataset(wset);
    free(felist);

    return femod;
}

/* Construct a gretl list containing the index numbers of the
   cross-sectional units included in the fixed-effects
   regression. This is for the purpose of naming the per-unit
   intercepts.
*/

int * fe_units_list(const(PanelmodT) * pan)
{
    int * ulist = null;
    int i, j, n = 0;

    for (i=0; i<pan.nunits; i++) {
        if (pan.unit_obs[i] > 0) {
            n++;
        }
    }

    ulist = gretl_list_new(n);

    if (ulist != null) {
        j = 1;
        for (i=0; i<pan.nunits; i++) {
            if (pan.unit_obs[i] > 0) {
                ulist[j++] = i + 1;
            }
        }
    }

    return ulist;
}

/* Compose a list referencing all variables that were dropped from the
   final panel model relative to the incoming regression
   specification.  This may include some variables that were dropped
   at the stage of running the baseline pooled model (presumably
   because of perfect collinearity).

   In the case of fixed effects it may include additional variables
   dropped due to the fact that they are time-invariant.
*/

int compose_panel_droplist(Model * pmod, PanelmodT * pan)
{
    int fixed_effects = (pmod.opt & Gretlopt.OPT_F);
    const(int) * pooldrop;
    int * dlist;
    int i, ndrop = 0;

    /* regressors dropped at the stage of estimating the
       initial pooled model */
    pooldrop = gretl_model_get_list(pan.pooled, "droplist");
    if (pooldrop != null) {
        ndrop += pooldrop[0];
    }

    if (fixed_effects) {
        /* regressors dropped because time-invariant */
        ndrop += pan.pooled.list[0] - pan.vlist[0];
    }

    if (ndrop == 0) {
        /* nothing to be done */
        return 0;
    }

    dlist = gretl_list_new(ndrop);
    if (dlist == null) {
        return GretlError.alloc;
    }

    i = 1;

    if (pooldrop != null) {
        for (i=1; i<=pooldrop[0]; i++) {
            dlist[i] = pooldrop[i];
        }
    }

    if (fixed_effects) {
        if (pan.vlist[0] < pan.pooled.list[0]) {
            int j, vj;

            for (j=2; j<=pan.pooled.list[0]; j++) {
                vj = pan.pooled.list[j];
                if (!in_gretl_list(pan.vlist, vj)) {
                    dlist[i++] = vj;
                }
            }
        }
    }

    return gretl_model_set_list_as_data(pmod, "droplist", dlist);
}

/* joint chi-square test on regular regressors */

void add_gls_joint_test(Model * pmod, PanelmodT * pan)
{
    double chisq = NADBL;
    int df = 0;

    if (pan.ntdum > 0) {
        df = pmod.ncoeff - 1 - pan.ntdum;
        if (df > 0) {
            chisq = panel_overall_test(pmod, pan, pan.ntdum,
                                       Gretlopt.OPT_X);
        }
    } else {
        df = pmod.ncoeff - 1;
        if (df > 0) {
            chisq = panel_overall_test(pmod, pan, 0, Gretlopt.OPT_X);
        }
    }

    if (!na(chisq) && chisq >= 0.0) {
        ModelTest *test = model_test_new(ModelTestType.re_wald);

        if (test != null) {
            model_test_set_teststat(test, TestStats.wald_chisq);
            model_test_set_dfn(test, df);
            model_test_set_value(test, chisq);
            model_test_set_pvalue(test, chisq_cdf_comp(df, chisq));
            maybe_add_test_to_model(pmod, test);
        }
    }
}

void add_gls_stats(Model * pmod, PanelmodT * pan)
{
    int nc;

    fix_panelmod_list(pmod, pan);

    pmod.ybar = pan.pooled.ybar;
    pmod.sdy = pan.pooled.sdy;
    pmod.tss = pan.pooled.tss;

    pmod.rsq = NADBL;
    pmod.adjrsq = NADBL;
    pmod.fstt = NADBL;
    pmod.rho = pan.rho;
    pmod.dw = pan.dw;

    /* deferred rewriting of ess and sigma, etc. */
    pmod.ess = gretl_model_get_double(pmod, "fixed_SSR");
    pmod.sigma = sqrt(pmod.ess / (pmod.nobs - (pmod.ncoeff - 1)));
    gretl_model_destroy_data_item(pmod, "fixed_SSR");
    nc = pmod.ncoeff;
    pmod.ncoeff = pmod.dfn + 1;
    ls_criteria(pmod);
    pmod.ncoeff = nc;

    add_gls_joint_test(pmod, pan);
}

void add_panel_obs_info(Model * pmod, PanelmodT * pan)
{
    gretl_model_set_int(pmod, "n_included_units", pan.effn);
    gretl_model_set_int(pmod, "panel_T", pan.T);
    gretl_model_set_int(pmod, "Tmin", pan.Tmin);
    gretl_model_set_int(pmod, "Tmax", pan.Tmax);
    if (pan.Tmax > pan.Tmin && pan.Tbar == 0) {
        calculate_Tbar(pan);
    }
    if (pan.Tbar > 0) {
        gretl_model_set_double(pmod, "Tbar", pan.Tbar);
    }
}

void replace_re_residuals(Model * pmod, PanelmodT * pan)
{
    int t;

    /* replace the GLS residuals with the "fixed up" ones */

    for (t=0; t<pmod.full_n; t++) {
        pmod.uhat[t] = pan.re_uhat[t];
    }
}

/* We use this to "finalize" models estimated via fixed effects
   and random effects */

int save_panel_model(Model * pmod, PanelmodT * pan, const(double) ** Z, const(DatasetRaw) * dset)
{
    int err = 0;

    pmod.ci = GretlCmdIndex.panel;

    if (pan.opt & Gretlopt.OPT_F) {
        /* fixed effects */
        int * ulist;

        pmod.opt |= Gretlopt.OPT_F;
        err = fix_within_stats(pmod, pan);
        if (err) {
            return err;
        }
        ulist = fe_units_list(pan);
        gretl_model_add_panel_varnames(pmod, dset, ulist);
        free(ulist);
        panel_model_add_ahat(pmod, dset, pan);
        save_fixed_effects_F(pan, pmod);
    } else {
        /* random effects */
        pmod.opt |= Gretlopt.OPT_U;
        gretl_model_set_double(pmod, "s2v", pan.s2v);
        gretl_model_set_double(pmod, "s2e", pan.s2e);
        if (pan.balanced) {
            gretl_model_set_double(pmod, "theta", pan.theta);
        } else if (!na(pan.theta_bar)) {
            pan.theta_bar /= pan.effn;
            gretl_model_set_double(pmod, "theta_bar", pan.theta_bar);
        }
        gretl_model_add_panel_varnames(pmod, dset, null);
        if (pan.re_uhat != null) {
            replace_re_residuals(pmod, pan);
        }
        add_gls_stats(pmod, pan);
        panel_model_add_ahat(pmod, dset, pan);
        if (pan.opt & Gretlopt.OPT_E) {
            /* record use of Nerlove transformation */
            pmod.opt |= Gretlopt.OPT_E;
        }
        if ((pan.opt & Gretlopt.OPT_X) && !IGLS) {
            /* record use of special unbalanced ANOVA */
            pmod.opt |= Gretlopt.OPT_X;
            if (!(pan.opt & Gretlopt.OPT_E)) {
                const(char) * meth = stata_sa ? "stata" : "bc";

                gretl_model_set_string_as_data(pmod, "anova_method",
                                               gretl_strdup(meth));
            }
        }
    }

    /* compose list of dropped variables, if any */
    compose_panel_droplist(pmod, pan);

    if (!(pan.opt & Gretlopt.OPT_A)) {
        set_model_id(pmod, pan.opt);
    }

    if (pan.opt & Gretlopt.OPT_F) {
        panel_dwstat(pmod, pan);
    }

    if (pan.opt & Gretlopt.OPT_R) {
        pmod.opt |= Gretlopt.OPT_R;
    }

    time_dummies_wald_test(pan, pmod);
    add_panel_obs_info(pmod, pan);
    *pan.realmod = *pmod;

    return err;
}

/* drive the calculation of the fixed effects regression, print the
   results (if wanted), and compute the "within" error variance */

int within_variance(PanelmodT * pan, DatasetRaw * dset, Prn * prn)
{
    Model femod;
    int i, err = 0;

    femod = fixed_effects_model(pan, dset, prn);

    if (femod.errcode) {
        pputs(prn, _("Error estimating fixed effects model\n").ptr);
        errmsg(femod.errcode, prn);
        err = femod.errcode;
        clear_model(&femod);
    } else {
        int den;

        if (pan.opt & Gretlopt.OPT_E) {
            /* Nerlove */
            den = femod.nobs;
        } else {
            /* as per Greene: nT - n - K */
            den = femod.nobs - pan.effn - (pan.vlist[0] - 2);
        }

        if (den == 0) {
            gretl_errmsg_set(_("Inadequate data for panel estimation").ptr);
            return GretlError.df;
        } else {
            pan.s2e = femod.ess / den;
        }


        if (!(pan.opt & Gretlopt.OPT_R)) {
            regular_fixed_effects_F(pan, &femod);
        }
        if (pan.opt & Gretlopt.OPT_V) {
            print_fe_results(pan, &femod, dset, prn);
        }
        if (pan.bdiff != null && pan.Sigma != null) {
            for (i=1; i<femod.ncoeff; i++) {
                pan.bdiff.val[i-1] = femod.coeff[i];
            }
            femod_regular_vcv(&femod);
            vcv_slopes(pan, &femod, VcvOps.init);
        }
        if (pan.opt & Gretlopt.OPT_F) {
            err = save_panel_model(&femod, pan, cast(const(double) **) dset.Z,
                                   dset);
            if (err) {
                clear_model(&femod);
            }
        } else {
            if (pan.opt & Gretlopt.OPT_U) {
                if (expand_fe_uhat(&femod, pan) == 0) {
                    panel_dwstat(&femod, pan);
                }
            }
            clear_model(&femod);
        }
    }

    return err;
}

void print_hausman_result(PanelmodT * pan, Prn * prn)
{
    if (na(pan.H)) {
        pputs(prn, _("Hausman test matrix is not positive definite!\n").ptr);
    } else {
        pprintf(prn, _("Hausman test statistic:\n H = %g with p-value = prob(chi-square(%d) > %g) = %g\n").ptr,
                pan.H, pan.dfH, pan.H, chisq_cdf_comp(pan.dfH, pan.H));
        pputs(prn, _("(A low p-value counts against the null hypothesis that the random effects\nmodel is consistent, in favor of the fixed effects model.)\n").ptr);
    }
}

void save_hausman_result(PanelmodT * pan)
{
    ModelTest *test;

    if (pan.realmod == null || na(pan.H) || pan.dfH == 0) {
        return;
    }

    test = model_test_new(ModelTestType.panel_hausman);

    if (test != null) {
        model_test_set_teststat(test, TestStats.wald_chisq);
        model_test_set_dfn(test, pan.dfH);
        model_test_set_value(test, pan.H);
        if (na(pan.H)) {
            model_test_set_pvalue(test, NADBL);
        } else {
            model_test_set_pvalue(test, chisq_cdf_comp(pan.dfH, pan.H));
        }
        maybe_add_test_to_model(pan.realmod, test);
    }
}

/* Handle the case where collinear terms were dropped
   when doing the Hausman test via the regression
   method in robust mode.
*/

double robust_hausman_fixup(const(int) * hlist, Model * pmod)
{
    int * wlist = gretl_list_copy(hlist);
    double H = NADBL;

    if (wlist != null) {
        int i;

        for (i=wlist[0]; i>0; i--) {
            if (!in_gretl_list(pmod.list, wlist[i])) {
                gretl_list_delete_at_pos(wlist, i);
            }
        }
        H = wald_omit_chisq(wlist, pmod);
        free(wlist);
    }

    return H;
}

/* Estimate the augmented GLS model for the Hausman test;
   return either the error sum of squares or, in the
   robust case, a Wald chi-square statistic.
*/

double hausman_regression_result(PanelmodT * pan, const(int) * relist, const(int) * hlist, DatasetRaw * rset, DatasetRaw * dset, Prn * prn)
{
    double ret = NADBL;
    int * biglist = null;
    int err = 0;

    biglist = gretl_list_add(relist, hlist, &err);

    if (biglist != null) {
        Model hmod;

        gretl_model_init(&hmod, null);
        hmod = lsq(biglist, rset, GretlCmdIndex.ols, Gretlopt.OPT_A);
        if (hmod.errcode == 0) {
            /* Find the number of additional regressors actually used,
               relative to @relist, allowing for the possibility
               that one or more elements of @biglist were dropped
               in estimation of @hmod due to excessive collinearity.
            */
            pan.dfH = hmod.list[0] - relist[0];
            if (pan.dfH > 0) {
                if (pan.opt & Gretlopt.OPT_R) {
                    /* do robust Wald test */
                    if (hmod.full_n < pan.pooled.full_n) {
                        hausman_move_uhat(&hmod, pan);
                    }
                    panel_robust_vcv(&hmod, pan, rset, dset);
                    if (hmod.vcv != null) {
                        if (pan.dfH < hlist[0]) {
                            ret = robust_hausman_fixup(hlist, &hmod);
                        } else {
                            ret = wald_omit_chisq(hlist, &hmod);
                        }
                    }
                } else {
                    /* just record the unrestricted SSR */
                    ret = hmod.ess;
                }
            }
        }
        clear_model(&hmod);
    }

    free(biglist);

    return ret;
}

/* Computation of s2_v in the manner of Swamy and Arora, for an
   unbalanced panel. See Baltagi, Econometric Analysis of Panel
   Data, 3e, section 9.2.1. */

int unbalanced_SA_s2v(PanelmodT * pan, DatasetRaw * dset)
{
    GretlMatrixBlock * B;
    GretlMatrixRaw * ZmZ = null;
    GretlMatrixRaw * PZ = null;
    GretlMatrixRaw * Z = null;
    GretlMatrixRaw * ZPZ = null;
    GretlMatrixRaw * D2 = null;
    GretlMatrixRaw * trmat = null;
    int k = pan.pooled.ncoeff;
    double zjt, z, tr;
    int i, j, s, t, p;
    int bigt, got;
    int err = 0;

    B = gretl_matrix_block_new(&ZmZ, pan.effn, k,
                               &PZ,  pan.NT, k,
                               &Z,   pan.NT, k,
                               &ZPZ, k, k,
                               &D2,  k, k, null);

    if (B == null) {
        return GretlError.alloc;
    }

    s = p = 0;
    for (i=0; i<pan.nunits; i++) {
        int Ti = pan.unit_obs[i];

        if (Ti == 0) {
            continue;
        }

        for (j=0; j<k; j++) {
            zjt = 0.0;
            got = 0;
            /* cumulate sum of observations for unit */
            for (t=0; t<pan.T && got<Ti; t++) {
                bigt = panel_index(i, t);
                if (!panel_missing(pan, bigt)) {
                    z = dset.Z[pan.pooled.list[j+2]][bigt];
                    zjt += z;
                    gretl_matrix_set(Z, p + got, j, z);
                    got++;
                }
            }
            gretl_matrix_set(ZmZ, s, j, zjt);
            for (t=0; t<Ti; t++) {
                gretl_matrix_set(PZ, p + t, j, zjt / Ti);
            }
        }
        s++;
        p += Ti;
    }

    gretl_matrix_multiply_mod(Z, GretlMatrixMod.transpose,
                              PZ, GretlMatrixMod.none,
                              ZPZ, GretlMatrixMod.none);
    gretl_invert_symmetric_matrix(ZPZ);

    gretl_matrix_multiply_mod(ZmZ, GretlMatrixMod.transpose,
                              ZmZ, GretlMatrixMod.none,
                              D2, GretlMatrixMod.none);

    trmat = gretl_matrix_reuse(PZ, k, k);
    gretl_matrix_multiply(ZPZ, D2, trmat);
    tr = gretl_matrix_trace(trmat);


    pan.s2v = (pan.ubPub - (pan.effn - k) * pan.s2e) / (pan.NT - tr);


    gretl_matrix_block_destroy(B);

    return err;
}

/* Calculate the random effects regression.  Print the results
   here if we're doing the "panel diagnostics" test, otherwise
   save the results.
*/

int random_effects(PanelmodT * pan, DatasetRaw * dset, DatasetRaw * gset, Prn * prn)
{
    DatasetRaw * rset;
    Model remod;
    Gretlopt lsqopt = Gretlopt.OPT_A | Gretlopt.OPT_Z;
    int * relist = null;
    int * hlist = null;
    double hres = NADBL;
    int i, err = 0;

    gretl_model_init(&remod, dset);

    /* FGLS regression list */
    relist = gretl_list_new(pan.pooled.list[0]);
    if (relist == null) {
        return GretlError.alloc;
    }

    if (!(pan.opt & Gretlopt.OPT_M)) {
        /* extra regressors for Hausman test, regression approach */
        hlist = gretl_list_new(pan.vlist[0] - 1);
        if (hlist == null) {
            free(relist);
            return GretlError.alloc;
        }
    }

    /* If OPT_E (--nerlove) was given, we've already calculated
       pan->s2v as the variance of the fixed effects. Otherwise
       we're using the Swamy and Arora method, and pan->s2v still
       needs to be computed.

       Note: for unbalanced panels, theta will vary across the
       units in the final calculation.
    */
    if (!(pan.opt & Gretlopt.OPT_E)) {
        if (IGLS) {
            /* get user-specified values */
            err = read_true_variances(pan);
        } else if (!pan.balanced && !na(pan.ubPub)) {
            err = unbalanced_SA_s2v(pan, dset);
        } else {
            pan.s2v = pan.s2b - pan.s2e / pan.Tbar;
        }
        if (pan.s2v < 0) {
            pan.s2v = 0.0;
        }
    }

    /* theta, the quasi-demeaning coefficient */
    pan.theta = 1.0 - sqrt(pan.s2e / (pan.s2e + pan.Tmax * pan.s2v));


    /* make special transformed dataset, and regression list */
    rset = random_effects_dataset(dset, gset, relist, hlist, pan);
    if (rset == null) {
        free(relist);
        free(hlist);
        return GretlError.alloc;
    }

    if (hlist != null) {
        hres = hausman_regression_result(pan, relist, hlist,
                                         rset, dset, prn);
    }

    /* regular random-effects model */
    remod = lsq(relist, rset, GretlCmdIndex.ols, lsqopt);

    if ((err = remod.errcode) != 0) {
        pputs(prn, _("Error estimating random effects model\n").ptr);
        errmsg(err, prn);
    } else {
        if (pan.bdiff != null) {
            /* matrix-diff Hausman variant, if wanted */
            int vi, k = 0;

            for (i=0; i<remod.ncoeff; i++) {
                vi = pan.pooled.list[i+2];
                if (var_is_varying(pan, vi)) {
                    pan.bdiff.val[k++] -= remod.coeff[i];
                }
            }
        }

        /* note: this call moved to here 2017-10-18 so we
           get computation of robust standard errors right
        */
        fix_panel_hatvars(&remod, pan, cast(const(double) **) dset.Z);

        if (pan.opt & Gretlopt.OPT_R) {
            err = panel_robust_vcv(&remod, pan, rset, dset);
            if (err) {
                goto bailout;
            }
        } else {
            double sigma = remod.sigma; /* or... ? */
            makevcv(&remod, sigma);
        }

        if (pan.opt & Gretlopt.OPT_V) {
            print_re_results(pan, &remod, dset, prn);
        }

        if (!na(hres)) {
            if (pan.opt & Gretlopt.OPT_R) {
                /* @hres is already the (robust) Hausman test statistic */
                pan.H = hres;
            } else {
                /* @hres is the unrestricted ess */
                pan.H = (remod.ess / hres - 1.0) * (remod.nobs);
            }
        } else if (pan.Sigma != null) {
            vcv_slopes(pan, &remod, VcvOps.subtract);
        }
    }


 bailout:

    if (!err && (pan.opt & Gretlopt.OPT_U)) {
        save_panel_model(&remod, pan, cast(const(double) **) dset.Z, rset);
    } else {
        clear_model(&remod);
    }

    destroy_dataset(rset);

    free(relist);
    free(hlist);

    return err;
}

void save_breusch_pagan_result(PanelmodT * pan)
{
    ModelTest *test;

    if (pan.realmod == null || na(pan.BP)) {
        return;
    }

    test = model_test_new(ModelTestType.panel_bp);

    if (test != null) {
        model_test_set_teststat(test, TestStats.wald_chisq);
        model_test_set_dfn(test, 1);
        model_test_set_value(test, pan.BP);
        model_test_set_pvalue(test, chisq_cdf_comp(1, pan.BP));
        maybe_add_test_to_model(pan.realmod, test);
    }
}

/* do the panel Breusch-Pagan test and either print or save
   the results */

int breusch_pagan_LM(PanelmodT * pan, Prn * prn)
{
    double A = 0.0;
    int n = pan.pooled.nobs;
    int print_means = 0;
    int i, t, M = 0;

    if ((pan.opt & Gretlopt.OPT_V) && pan.effn <= 10) {
        print_means = 1;
    }

    if (print_means) {
        pputs(prn, _("Means of pooled OLS residuals for cross-sectional units:").ptr);
        pputs(prn, "\n\n");
    }

    for (i=0; i<pan.nunits; i++) {
        int Ti = pan.unit_obs[i];

        if (Ti > 0) {
            double u, usum = 0.0;

            for (t=0; t<pan.T; t++) {
                u = pan.pooled.uhat[panel_index(i, t)];
                if (!na(u)) {
                    usum += u;
                }
            }
            if (print_means) {
                pprintf(prn, _(" unit %2d: %13.5g\n").ptr, i + 1,
                        usum / cast(double) Ti);
            }
            A += usum * usum;
            M += Ti * Ti;
        }
    }

    pan.BP = (cast(double) n * n /(2.0 * (M - n))) * pow((A / pan.pooled.ess) - 1.0, 2);

    if (pan.opt & Gretlopt.OPT_V) {
        pputc(prn, '\n');
        pprintf(prn, _("Breusch-Pagan test statistic:\n LM = %g with p-value = prob(chi-square(1) > %g) = %g\n").ptr,
                pan.BP, pan.BP, chisq_cdf_comp(1, pan.BP));

        pputs(prn, _("(A low p-value counts against the null hypothesis that the pooled OLS model\nis adequate, in favor of the random effects alternative.)\n\n").ptr);
    }

    return 0;
}

void save_panspec_result(PanelmodT * pan, Prn * prn)
{
    GretlMatrixRaw * tests = gretl_column_vector_alloc(3);
    GretlMatrixRaw * pvals = gretl_column_vector_alloc(3);
    char ** S = strings_array_new(3);

    /* fixed-effects poolability F-test */
    tests.val[0] = pan.Ffe;
    pvals.val[0] = snedecor_cdf_comp(pan.Fdfn, pan.Fdfd, pan.Ffe);
    /* random-effects poolability test */
    tests.val[1] = pan.BP;
    pvals.val[1] = chisq_cdf_comp(1, pan.BP);
    /* Hausman test */
    tests.val[2] = pan.H;
    pvals.val[2] = chisq_cdf_comp(pan.dfH, pan.H);

    if (S != null) {
        /* add row names */
        char ** S2;

        S[0] = gretl_strdup("Poolability (Wald)");
        S[1] = gretl_strdup("Poolability (B-P)");
        S[2] = gretl_strdup("Hausman");
        S2 = strings_array_dup(S, 3);
        gretl_matrix_set_rownames(tests, S);
        if (S2 != null) {
            gretl_matrix_set_rownames(pvals, S2);
        }
    }

    if (pan.opt & Gretlopt.OPT_V) {
        print_hausman_result(pan, prn);
    }
    record_matrix_test_result(tests, pvals);
}

int finalize_hausman_test(PanelmodT * pan, int ci, Prn * prn)
{
    int mdiff = 0, err = 0;

    if (pan.bdiff != null && pan.Sigma != null) {
        /* matrix approach */
        mdiff = 1;
        err = bXb(pan);
    } else if (na(pan.H)) {
        /* regression approach bombed somehow? */
        err = GretlError.data;
    }

    if (ci == GretlCmdIndex.panspec) {
        /* the context is the "panspec" command */
        if (!err || (mdiff && err == GretlError.notpd)) {
            save_panspec_result(pan, prn);
        }
    } else {
        /* the context is random-effects estimation */
        if (mdiff && err == GretlError.notpd) {
            pputs(prn, _("Hausman test matrix is not positive definite").ptr);
            pputc(prn, '\n');
        }
        save_hausman_result(pan);
    }

    return err;
}

/* Based on the residuals from pooled OLS, do some accounting to see
   (a) how many cross-sectional units were actually included (after
   omitting any missing values); (b) how many time-series observations
   were included for each unit; and (c) what was the maximum number
   of time-series observations used.  Return 0 if all goes OK,
   non-zero otherwise.
*/

int panel_obs_accounts(PanelmodT * pan)
{
    int * uobs;
    int i, t, bigt;

    uobs = cast(typeof(uobs)) malloc(pan.nunits * (*uobs).sizeof);
    if (uobs == null) {
        return GretlError.alloc;
    }

    pan.NT = 0;        /* to hold total complete observations */
    pan.effn = 0;      /* to hold number of included units */
    pan.Tmax = 0;      /* to hold max time-series length */
    pan.Tmin = pan.T; /* to hold min time-series length > 0 */

    for (i=0; i<pan.nunits; i++) {
        uobs[i] = 0;
        for (t=0; t<pan.T; t++) {
            bigt = panel_index(i, t);
            if (!panel_missing(pan, bigt)) {
                uobs[i] += 1;
            }
        }
        if (uobs[i] > 0) {
            pan.effn += 1;
            if (uobs[i] > pan.Tmax) {
                pan.Tmax = uobs[i];
            }
            if (uobs[i] < pan.Tmin) {
                pan.Tmin = uobs[i];
            }
            pan.NT += uobs[i];
        }
    }

    for (i=0; i<pan.nunits; i++) {
        if (uobs[i] > 0 && uobs[i] != pan.Tmax) {
            pan.balanced = 0;
            break;
        }
    }

    pan.unit_obs = uobs;

    return 0;
}

/* Construct an array of char to record which parameters in the full
   model correspond to time-varying regressors
*/

int panel_set_varying(PanelmodT * pan, const(Model) * pmod)
{
    int i;

    pan.varying = cast(typeof(pan.varying)) calloc(pmod.ncoeff, 1);
    if (pan.varying == null) {
        return GretlError.alloc;
    }

    for (i=0; i<pmod.ncoeff; i++) {
        if (var_is_varying(pan, pmod.list[i+2])) {
            pan.varying[i] = 1;
        }
    }

    return 0;
}

/* find harmonic mean of the number of time-series observations
   per included group */

void calculate_Tbar(PanelmodT * pan)
{
    int i;

    if (pan.balanced) {
        pan.Tbar = pan.Tmax;
    } else {
        double den = 0.0;

        for (i=0; i<pan.nunits; i++) {
            if (pan.unit_obs[i] > 0) {
                den += 1.0 / pan.unit_obs[i];
            }
        }
        pan.Tbar = pan.effn / den;
    }
}

int panelmod_setup(PanelmodT * pan, Model * pmod, const(DatasetRaw) * dset, int ntdum, Gretlopt opt)
{
    int nobs, dsetT = sample_size(dset);
    int err = 0;

    pan.opt = opt;
    pan.pooled = pmod;

    /* assumes (possibly padded) balanced panel dataset */
    pan.nunits = dsetT / dset.pd;
    pan.T = dset.pd;
    nobs = pan.nunits * pan.T;

    if (nobs < dsetT) {
        fprintf(stderr, "dataset nobs = %d, but panel nobs = %d\n".ptr,
                dsetT, nobs);
    }

    panel_index_init(dset, pan.nunits, pan.T);
    pan.ntdum = ntdum;

    err = panel_obs_accounts(pan);

    if (!err && (pan.opt & (Gretlopt.OPT_U | Gretlopt.OPT_F | Gretlopt.OPT_B))) {
        pan.realmod = cast(typeof(pan.realmod)) malloc((*pan.realmod).sizeof);
        if (pan.realmod == null) {
            err = GretlError.alloc;
        }
    }

    IGLS = stata_sa = 0;

    if (!err && (opt & Gretlopt.OPT_X)) {
        if (!(pan.opt & Gretlopt.OPT_U)) {
            /* --unbalanced --random-effects */
            err = GretlError.badopt;
        } else {
            /* probe optional argument to --unbalanced */
            const(char) * s = get_optval_string(GretlCmdIndex.panel, Gretlopt.OPT_X);

            if (s != null) {
                if (!strcmp(s, "stata")) {
                    stata_sa = 1;
                } else if (!strcmp(s, "bc")) {
                    stata_sa = 0;
                } else if (has_suffix(s, ".mat")) {
                    /* hidden testing option */
                    strcpy(glsmat.ptr, s);
                    IGLS = 1;
                } else {
                    gretl_errmsg_sprintf(_("%s: invalid option argument").ptr, s);
                    err = GretlError.invarg;
                }
            }
        }
    }

    if (err && pan.unit_obs != null) {
        free(pan.unit_obs);
        pan.unit_obs = null;
    }

    return err;
}

/* Called in relation to a model estimated by pooled OLS: test for
   both fixed and random effects. Implements the "panspec" command.
*/

int panel_diagnostics(Model * pmod, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    int nerlove = 0;
    int quiet = (opt & Gretlopt.OPT_Q);
    Gretlopt psopt = opt;
    PanelmodT pan;
    int xdf, err = 0;


    if (pmod.ifc == 0) {
        /* at many points we assume the base regression has an
           intercept included */
        return GretlError.noconst;
    }


    if (opt & Gretlopt.OPT_E) {
        nerlove = 1;
    }

    /* Add OPT_V to make the fixed and random effects functions verbose,
       unless we've been passed the --quiet option.
    */
    if (!quiet) {
        psopt |= Gretlopt.OPT_V;
    }
    panelmod_init(&pan);
    err = panelmod_setup(&pan, pmod, dset, 0, psopt);
    if (err) {
        goto bailout;
    }

    if (pan.effn < pan.nunits) {
        fprintf(stderr, "number of units included = %d\n".ptr, pan.effn);
        if (pan.effn <= 0) {
            return GretlError.data;
        }
    }

    /* figure out which of the original regressors are time-varying */
    err = varying_vars_list(dset, &pan);
    if (err) {
        goto bailout;
    }

    err = panel_set_varying(&pan, pmod);
    if (err) {
        goto bailout;
    }

    calculate_Tbar(&pan);

    /* degrees of freedom relative to # of x-sectional units */
    xdf = pan.effn - pmod.ncoeff;


    if (xdf <= 0 && !nerlove) {
        ; /* can't do Random Effects */
    } else if (pan.opt & Gretlopt.OPT_M) {
        /* matrix version of Hausman test */
        err = matrix_hausman_allocate(&pan);
        if (err) {
            goto bailout;
        }
    }

    if (!quiet) {
        /* header */
        pputc(prn, '\n');
        pprintf(prn, _("Diagnostics: using n = %d cross-sectional units\n").ptr,
                pan.effn);
        pputc(prn, '\n');
    }

    err = within_variance(&pan, dset, prn);
    if (err) {
        goto bailout;
    }

    if (xdf <= 0) {
        pprintf(prn, "Omitting group means regression: insufficient degrees of freedom\n".ptr);
        goto bailout;
    }

    if (xdf > 0 && !na(pan.s2e)) {
        DatasetRaw * gset;

        gset = group_means_dataset(&pan, dset);
        if (gset == null) {
            err = GretlError.alloc;
        } else {
            err = between_variance(&pan, gset);
        }

        if (err) {
            pputs(prn, _("Couldn't estimate group means regression\n").ptr);
            if (err == GretlError.singular) {
                /* treat this as non-fatal */
                err = 0;
            }
        } else {
            random_effects(&pan, dset, gset, prn);
            if (pan.theta > 0) {
                /* this test hardly makes sense if GLS = OLS */
                breusch_pagan_LM(&pan, prn);
            }
            finalize_hausman_test(&pan, GretlCmdIndex.panspec, prn);
        }

        if (gset != null) {
            destroy_dataset(gset);
        }
    }

 bailout:

    panelmod_free(&pan);

    return err;
}

int between_model(PanelmodT * pan, const(DatasetRaw) * dset)
{
    DatasetRaw * gset;
    int err = 0;

    gset = group_means_dataset(pan, dset);

    if (gset == null) {
        err = GretlError.alloc;
    } else {
        err = between_variance(pan, gset);
        if (err) {
            destroy_dataset(gset);
        } else {
            pan.realmod.dataset = gset;
        }
    }

    return err;
}

/* When we automatically add dummy variables to the dataset,
   in estimating a panel model with the --time-dummies option,
   should we delete these when we're done, or leave them
   in the dataset? 2013-11-18: we'll delete them if the dataset
   is subsampled, otherwise leave them.
*/

int process_time_dummies(Model * pmod, const(DatasetRaw) * dset, int v)
{
    int i, vi, n = 0;

    for (i=pmod.list[0]; i>1; i--) {
        if (pmod.list[i] >= v) {
            n++;
        }
    }

    if (n > 0) {
        gretl_model_allocate_param_names(pmod, pmod.ncoeff);

        if (pmod.errcode) {
            return pmod.errcode;
        }

        for (i=2; i<=pmod.list[0]; i++) {
            vi = pmod.list[i];
            gretl_model_set_param_name(pmod, i-2, dset.varname[vi]);
        }

        for (i=pmod.list[0]; i>1; i--) {
            if (pmod.list[i] >= v) {
                gretl_list_delete_at_pos(pmod.list, i);
            }
        }
    }

    /* Record the fact that the model used time dummies, in case
       the user wants to add/omit variables or something.
    */
    pmod.opt |= Gretlopt.OPT_D;

    return 0;
}

int add_time_dummies_to_list(const(int) * list, DatasetRaw * dset, int ** plist)
{
    char[VNAMELEN] dname;
    int * biglist = null;
    int ntdum = dset.pd - 1;
    int i, j, v;
    int err = 0;

    biglist = gretl_list_new(list[0] + ntdum);
    if (biglist == null) {
        return GretlError.alloc;
    }

    for (i=1; i<=list[0]; i++) {
        biglist[i] = list[i];
    }

    j = list[0] + 1;

    for (i=2; i<=dset.pd; i++) {
        sprintf(dname.ptr, "dt_%d".ptr, i);
        v = series_greatest_index(dset, dname.ptr);
        if (v == dset.v) {
            err = GretlError.data;
            break;
        } else if (in_gretl_list(list, v)) {
            /* already present */
            biglist[0] -= 1;
        } else {
            biglist[j++] = v;
        }
    }

    if (err) {
        free(biglist);
    } else {
        *plist = biglist;
    }

    return err;
}

int panel_check_for_const(const(int) * list)
{
    int i;

    for (i=2; i<=list[0]; i++) {
        if (list[i] == 0) {
            return 0; /* OK */
        }
    }

    return GretlError.noconst;
}

int get_ntdum(const(int) * orig, const(int) * new_)
{
    int i, n = 0;

    for (i=2; i<=new_[0]; i++) {
        if (!in_gretl_list(orig, new_[i])) {
            n++;
        }
    }

    return n;
}

int save_pooled_model(Model * pmod, PanelmodT * pan, DatasetRaw * dset)
{
    int err = 0;

    if (pan.opt & Gretlopt.OPT_R) {
        err = panel_robust_vcv(pmod, pan, dset, dset);
    if (!err) {
        pmod.opt |= Gretlopt.OPT_R;
        pmod.fstt = panel_overall_test(pmod, pan, 0, Gretlopt.OPT_F);
        pmod.dfd = pan.effn - 1;
    }
    }

    if (!err) {
    gretl_model_set_int(pmod, "pooled", 1);
    add_panel_obs_info(pmod, pan);
    if (!(pan.opt & Gretlopt.OPT_A)) {
        set_model_id(pmod, pan.opt);
    }
    panel_dwstat(pmod, pan);
    }

    return err;
}

/* fixed effects | random effects | between | pooled */

bool estimator_specified(Gretlopt o)
{
    return (o & (Gretlopt.OPT_F|Gretlopt.OPT_U|Gretlopt.OPT_B|Gretlopt.OPT_P)) != 0;
}

/**
 * real_panel_model:
 * @list: list containing model specification.
 * @dset: dataset struct.
 * @opt: may include %OPT_U for the random effects model;
 * %OPT_R for robust standard errors; %OPT_M to use the
 * matrix-difference variant of the Hausman test (random
 * effects only); %OPT_B for the "between" model; %OPT_P for
 * pooled OLS; and %OPT_D to include time dummies.
 * %OPT_C for clustered standard errors is also accepted.
 * If %OPT_U is given, either of the mutually incompatible options
 * %OPT_E and %OPT_X may be given to inflect the calculation of the
 * variance of the individual effects: %OPT_E means use Nerlove's
 * method, %OPT_X means use the special "exact" method of Swamy
 * and Arora in the case of an unbalanced panel. The default is
 * to use the "generic" Swamy-Arora method.
 * @prn: printing struct.
 *
 * Estimates a panel model, by default using the fixed effects
 * estimator.
 *
 * Returns: a #MODEL struct, containing the estimates.
 */

Model real_panel_model(const(int) * list, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    int orig_v = dset.v;
    Model mod;
    PanelmodT pan;
    Gretlopt pan_opt = opt;
    Gretlopt ols_opt = Gretlopt.OPT_A;
    int * olslist = null;
    int nerlove = 0;
    int ntdum = 0;
    int err = 0;

    gretl_model_init(&mod, dset);
    panelmod_init(&pan);

    if (!estimator_specified(opt)) {
        /* default: add OPT_F to save the fixed effects model */
        pan_opt |= Gretlopt.OPT_F;
    }

    if (opt & Gretlopt.OPT_P) {
        /* doing pooled OLS */
        if (opt & Gretlopt.OPT_N) {
            /* no-df-corr */
            ols_opt |= Gretlopt.OPT_N;
        }
    } else {
        /* all other variants */
        mod.errcode = panel_check_for_const(list);
        if (mod.errcode) {
            return mod;
        }
    }

    if (opt & Gretlopt.OPT_D) {
        /* time dummies wanted */
        err = gen_panel_dummies(dset, Gretlopt.OPT_T, prn);
        if (!err) {
            err = add_time_dummies_to_list(list, dset, &olslist);
        }
    } else {
        olslist = gretl_list_copy(list);
        if (olslist == null) {
            err = GretlError.alloc;
        }
    }

    if (err) {
        goto bailout;
    }

    /* baseline: estimate via pooled OLS */
    mod = lsq(olslist, dset, GretlCmdIndex.ols, ols_opt);
    if (mod.errcode) {
        err = mod.errcode;
        fprintf(stderr, "real_panel_model: error %d in initial OLS\n".ptr,
                mod.errcode);
        goto bailout;
    }


    free(olslist);


    if (pan_opt & Gretlopt.OPT_C) {
        /* cluster implies robust */
        pan_opt |= Gretlopt.OPT_R;
    }
    if (pan_opt & Gretlopt.OPT_N) {
        /* no-df-corr (not working!)  */
        pan_opt |= Gretlopt.OPT_N;
    }

    if ((opt & Gretlopt.OPT_D) && !(opt & Gretlopt.OPT_B)) {
        ntdum = get_ntdum(list, mod.list);
    }

    err = panelmod_setup(&pan, &mod, dset, ntdum, pan_opt);
    if (err) {
        goto bailout;
    }

    if (opt & Gretlopt.OPT_P) {
        err = save_pooled_model(&mod, &pan, dset);
        goto bailout;
    }

    if (opt & Gretlopt.OPT_B) {
        err = between_model(&pan, dset);
        goto bailout;
    }

    /* figure out which of the original regressors are time-varying,
       or unit-varying as the case may be
    */
    err = varying_vars_list(dset, &pan);
    if (err) {
        goto bailout;
    }

    err = panel_set_varying(&pan, &mod);
    if (err) {
        goto bailout;
    }

    if (opt & Gretlopt.OPT_U) {
        nerlove = (opt & Gretlopt.OPT_E)? 1 : 0;
    }

    calculate_Tbar(&pan);

    if (opt & Gretlopt.OPT_U) {
        /* trying to do random effects */
        int xdf = pan.effn - (mod.ncoeff - pan.ntdum);

        if (xdf <= 0 && !nerlove) {
            gretl_errmsg_set(_("Couldn't estimate group means regression").ptr);
            err = mod.errcode = GretlError.df;
            goto bailout;
        } else if (pan.opt & Gretlopt.OPT_M) {
            /* matrix version of Hausman test */
            err = matrix_hausman_allocate(&pan);
            if (err) {
                goto bailout;
            }
        }
    }

    err = within_variance(&pan, dset, prn);
    if (err) {
        goto bailout;
    }

    if ((opt & Gretlopt.OPT_U) && !na(pan.s2e)) {
        DatasetRaw * gset;

        breusch_pagan_LM(&pan, prn);

        gset = group_means_dataset(&pan, dset);
        if (gset == null) {
            err = GretlError.alloc;
        } else {
            err = between_variance(&pan, gset);
        }

        if (err && !nerlove) {
            pputs(prn, _("Couldn't estimate group means regression\n").ptr);
        } else {
            err = random_effects(&pan, dset, gset, prn);
            if (!err) {
                save_breusch_pagan_result(&pan);
                finalize_hausman_test(&pan, GretlCmdIndex.panel, prn);
            }
        }

        if (gset != null) {
            destroy_dataset(gset);
        }
    }

 bailout:

    if (!err) {
        if ((pan_opt & (Gretlopt.OPT_F | Gretlopt.OPT_U)) && mod.missmask != null) {
            /* preserve missing obs mask, if any */
            char * mask = mod.missmask;

            mod.missmask = null;
            clear_model(&mod);
            mod = *pan.realmod;
            mod.missmask = mask;
        } else if (!(opt & Gretlopt.OPT_P)) {
            /* not doing pooled OLS, in which case @mod would
               already hold the correct model
            */
            clear_model(&mod);
            mod = *pan.realmod;
        }
        gretl_model_smpl_init(&mod, dset);
        if (opt & Gretlopt.OPT_D) {
            if (pan.ntdum > 0) {
                gretl_model_set_int(&mod, "ntdum", pan.ntdum);
                maybe_suppress_time_dummies(&mod, pan.ntdum);
            }
            if (complex_subsampled()) {
                process_time_dummies(&mod, dset, orig_v);
            }
        }
    }

    panelmod_free(&pan);

    if (err && mod.errcode == 0) {
        mod.errcode = err;
    }

    if (complex_subsampled()) {
        dataset_drop_last_variables(dset, dset.v - orig_v);
    }

    if (!err) {
    /* equip the panel model with an xlist member */
    int * dropped = gretl_model_get_list(&mod, "droplist");
    int * xlist = gretl_list_sublist(mod.list, 2, -1);
    int lerr = 0;

    free(mod.xlist);
    if (dropped != null) {
        int * tmp = gretl_list_drop(xlist, dropped, &lerr);

        mod.xlist = tmp;
        free(xlist);
    } else {
        mod.xlist = xlist;
    }
    }


    return mod;
}

/* Called from qr_estimate.c in case robust VCV estimation is called
   for in the context of TSLS estimation on panel data.
*/

int panel_tsls_robust_vcv(Model * pmod, DatasetRaw * dset)
{
    PanelmodT pan;
    int err = 0;

    panelmod_init(&pan);

    err = panelmod_setup(&pan, pmod, dset, 0, Gretlopt.OPT_NONE);
    if (!err) {
        err = panel_robust_vcv(pmod, &pan, dset, dset);
    }

    panelmod_free(&pan);

    return err;
}

/* write weights for groupwise weighted least squares into the
   last variable in the dataset */

int write_weights_to_dataset(double * uvar, int nunits, int T, DatasetRaw * dset)
{
    int w = dset.v - 1;
    double wi;
    int i, t;

    for (i=0; i<nunits; i++) {
        if (uvar[i] <= 0.0 || na(uvar[i])) {
            wi = 0.0;
        } else {
            wi = 1.0 / uvar[i]; /* sqrt? */
        }
        for (t=0; t<T; t++) {
            dset.Z[w][panel_index(i, t)] = wi;
        }
    }

    return 0;
}

int allocate_weight_var(DatasetRaw * dset)
{
    if (dataset_add_series(dset, 1)) {
        return GretlError.alloc;
    }

    strcpy(dset.varname[dset.v - 1], "unit_wt");

    return 0;
}

double max_coeff_diff(const(Model) * pmod, const(double) * bvec)
{
    double diff, maxdiff = 0.0;
    int i;

    for (i=0; i<pmod.ncoeff; i++) {
        diff = fabs(pmod.coeff[i] - bvec[i]);
        if (diff > maxdiff) {
            maxdiff = diff;
        }
    }

    return maxdiff;
}

enum S2MINOBS = 1;

/* Wald test for groupwise heteroskedasticity, without assuming
   normality of the errors: see Greene, 4e, p. 598.  Note that the
   computation involves a factor of (1/(T_i - 1)) so we cannot
   include groups that have only one observation.
*/

double wald_hetero_test(const(Model) * pmod, double s2, const(double) * uvar, PanelmodT * pan, int * df)
{
    double x, W = 0.0;
    int i, t, Ti;

    *df = 0;

    for (i=0; i<pan.nunits; i++) {
        double fii = 0.0;

        Ti = pan.unit_obs[i];
        if (Ti < 2) {
            continue;
        }

        for (t=0; t<pan.T; t++) {
            x = pmod.uhat[panel_index(i, t)];
            if (!na(x)) {
                x = x * x - uvar[i];
                fii += x * x;
            }
        }

        if (fii <= 0) {
            W = NADBL;
            break;
        }

        fii *= (1.0 / Ti) * (1.0 / (Ti - 1.0));
        x = uvar[i] - s2;
        W += x * x / fii;
        *df += 1;
    }

    if (*df < 2) {
        W = NADBL;
    }

    return W;
}

/* Likelihood-ratio test for groupwise heteroskedasticity: see Greene,
   4e, pp. 597 and 599.  This test requires that we're able to
   calculate the ML estimates (using iterated WLS).
*/

int ml_hetero_test(Model * pmod, double s2, const(double) * uvar, int nunits, const(int) * unit_obs)
{
    ModelTest *test;
    double x2, s2h = 0.0;
    int i, Ti, df = 0;
    int err = 0;

    for (i=0; i<nunits; i++) {
        Ti = unit_obs[i];
        if (Ti > 0) {
            s2h += Ti * log(uvar[i]);
            df++;
        }
    }

    x2 = pmod.nobs * log(s2) - s2h;
    df--;

    test = model_test_new(ModelTestType.groupwise);
    if (test != null) {
        model_test_set_teststat(test, TestStats.lr);
        model_test_set_dfn(test, df);
        model_test_set_value(test, x2);
        model_test_set_pvalue(test, chisq_cdf_comp(df, x2));
        maybe_add_test_to_model(pmod, test);
    } else {
        err = 1;
    }

    return err;
}

double pooled_ll(const(Model) * pmod)
{
    double n = pmod.nobs;

    return -(n / 2.0) * (1.0 + ConstId.pi - log(n) + log(pmod.ess));
}

/* calculate log-likelihood for panel MLE with groupwise
   heteroskedasticity
*/

void panel_ML_ll(Model * pmod, const(double) * uvar, int nunits, const(int) * unit_obs)
{
    double ll = -(pmod.nobs / 2.0) * ConstId.pi;
    int i, Ti;

    for (i=0; i<nunits; i++) {
        Ti = unit_obs[i];
        if (Ti > 0) {
            ll -= (Ti / 2.0) * (1.0 + log(uvar[i]));
        }
    }

    pmod.lnL = ll;
    mle_criteria(pmod, 0);
}

/* we can't estimate a group-specific variance based on just one
   observation */

int singleton_check(const(int) * unit_obs, int nunits)
{
    int i;

    for (i=0; i<nunits; i++) {
        if (unit_obs[i] == 1) {
            return 1;
        }
    }

    return 0;
}

/* compute per-unit error variances */

void unit_error_variances(double * uvar, const(Model) * pmod, PanelmodT * pan, int * df)
{
    int i, t;
    double uit;

    *df = 0;

    for (i=0; i<pan.nunits; i++) {
        if (pan.unit_obs[i] < S2MINOBS) {
            uvar[i] = NADBL;
        } else {
            uvar[i] = 0.0;
            for (t=0; t<pan.T; t++) {
                uit = pmod.uhat[panel_index(i, t)];
                if (!na(uit)) {
                    uvar[i] += uit * uit;
                }
            }
            uvar[i] /= pan.unit_obs[i];
            *df += 1;
        }
    }
}

void print_unit_variances(PanelmodT * pan, double * uvar, int wald_test, Prn * prn)
{
    int i, Ti;

    pputs(prn, " unit    variance\n");

    for (i=0; i<pan.nunits; i++) {
        Ti = pan.unit_obs[i];
        if (na(uvar[i]) || (wald_test && Ti < 2)) {
            pprintf(prn, "%5d%12s (T = %d)\n".ptr, i+1, "NA".ptr, Ti);
        } else {
            pprintf(prn, "%5d%#12g (T = %d)\n".ptr, i+1, uvar[i], Ti);
        }
    }
}

enum SMALLDIFF = 0.0001;
enum WLS_MAX = 30;

/**
 * panel_wls_by_unit:
 * @list: regression list.
 * @dset: dataset struct.
 * @opt: may include %OPT_I (iterate to ML solution),
 * %OPT_V for verbose operation.
 * @prn: for printing details of iterations (or %NULL).
 *
 * Performs weighted least squares (possibly iterated) on
 * a panel model, allowing for groupwise heteroskedasticity.
 *
 * Returns: the estimates, in a %MODEL.
 */

Model panel_wls_by_unit(const(int) * list, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    Model mdl;
    PanelmodT pan;
    Gretlopt wlsopt = Gretlopt.OPT_A | Gretlopt.OPT_Z;
    double * uvar = null;
    double * bvec = null;
    double s2, diff = 1.0;
    int * wlist = null;
    int orig_v = dset.v;
    int i, iter = 0;

    gretl_error_clear();
    panelmod_init(&pan);

    if (opt & Gretlopt.OPT_I) {
        /* iterating: no degrees-of-freedom correction in WLS */
        wlsopt |= Gretlopt.OPT_N;
    }

    /* baseline pooled model */
    mdl = lsq(list, dset, GretlCmdIndex.ols, Gretlopt.OPT_A);
    if (mdl.errcode) {
        goto bailout;
    }

    mdl.errcode = panelmod_setup(&pan, &mdl, dset, 0, Gretlopt.OPT_NONE);
    if (mdl.errcode) {
        goto bailout;
    }

    uvar = cast(typeof(uvar)) malloc(pan.nunits * (*uvar).sizeof);
    if (uvar == null) {
        free(pan.unit_obs);
        mdl.errcode = GretlError.alloc;
        return mdl;
    }

    if (opt & Gretlopt.OPT_I) {
    /* the --iterate option */
        if (singleton_check(pan.unit_obs, pan.nunits)) {
            pprintf(prn, _("Can't produce ML estimates: some units have only one observation").ptr);
            pputc(prn, '\n');
            opt ^= Gretlopt.OPT_I;
        }
    }

    s2 = mdl.ess / mdl.nobs;

    if ((opt & Gretlopt.OPT_V) && (opt & Gretlopt.OPT_I)) {
        pprintf(prn, _("\nOLS error variance = %g\n").ptr, s2);
        pprintf(prn, "log-likelihood = %g\n".ptr, pooled_ll(&mdl));
    }

    if (allocate_weight_var(dset)) {
        mdl.errcode = GretlError.alloc;
        goto bailout;
    }

    if (opt & Gretlopt.OPT_I) {
        bvec = cast(typeof(bvec)) malloc(mdl.ncoeff * (*bvec).sizeof);
        if (bvec == null) {
            mdl.errcode = GretlError.alloc;
            goto bailout;
        }
    }

    /* allocate and construct WLS regression list */
    wlist = gretl_list_new(mdl.list[0] + 1);
    if (wlist == null) {
        mdl.errcode = GretlError.alloc;
        goto bailout;
    }
    wlist[1] = dset.v - 1; /* weight variable: the last var added */
    for (i=2; i<=wlist[0]; i++) {
        wlist[i] = mdl.list[i-1];
    }

    /* If wanted (and possible) iterate to ML solution; otherwise just do
       one-step FGLS estimation.
    */
    while (diff > SMALLDIFF) {
        int df = 0;

        iter++;
        unit_error_variances(uvar, &mdl, &pan, &df);

        if (opt & Gretlopt.OPT_V) {
            if (opt & Gretlopt.OPT_I) {
                pprintf(prn, "\n*** %s %d ***\n".ptr, _("iteration").ptr,
                        iter);
            } else {
                pputc(prn, '\n');
            }
            print_unit_variances(&pan, uvar, 0, prn);
        }

        write_weights_to_dataset(uvar, pan.nunits, pan.T, dset);

        if (opt & Gretlopt.OPT_I) {
            /* save coefficients for comparison */
            for (i=0; i<mdl.ncoeff; i++) {
                bvec[i] = mdl.coeff[i];
            }
        }

        clear_model(&mdl);
        mdl = lsq(wlist, dset, GretlCmdIndex.wls, wlsopt);
        if (mdl.errcode) {
            break;
        } else if (iter > WLS_MAX) {
            mdl.errcode = GretlError.noconv;
            break;
        }

        if (opt & Gretlopt.OPT_I) {
            diff = max_coeff_diff(&mdl, bvec);
            if ((opt & Gretlopt.OPT_V) && iter == 1) {
                pprintf(prn, _("\nFGLS pooled error variance = %g\n").ptr,
                        mdl.ess / mdl.nobs);
            }
        } else {
            /* one-step FGLS */
            break;
        }
    }

    if (!mdl.errcode) {
        mdl.ci = GretlCmdIndex.panel;
        if (!(opt & Gretlopt.OPT_A)) {
            set_model_id(&mdl, pan.opt);
        }
        gretl_model_set_int(&mdl, "n_included_units", pan.effn);
        mdl.opt |= Gretlopt.OPT_H; /* --unit-weights */
        mdl.nwt = 0;

        if (opt & Gretlopt.OPT_I) {
            int df = 0;

            gretl_model_set_int(&mdl, "iters", iter);
            ml_hetero_test(&mdl, s2, uvar, pan.nunits, pan.unit_obs);
            unit_error_variances(uvar, &mdl, &pan, &df);
            panel_ML_ll(&mdl, uvar, pan.nunits, pan.unit_obs);
            if (opt & Gretlopt.OPT_V) {
                pputc(prn, '\n');
            }
        }
    }

 bailout:

    free(pan.unit_obs);
    free(uvar);
    free(wlist);
    free(bvec);

    if (!(opt & Gretlopt.OPT_V)) {
        dataset_drop_last_variables(dset, dset.v - orig_v);
    }

    return mdl;
}

void print_wald_test(double W, int df, double pval, PanelmodT * pan, double * uvar, double s2, Prn * prn)
{
    pprintf(prn, "\n%s:\n".ptr,
            _("Distribution free Wald test for heteroskedasticity").ptr);
    pprintf(prn, " %s(%d) = %g, ".ptr,  _("Chi-square").ptr, df, W);
    pprintf(prn, "%s = %g\n\n".ptr, _("with p-value").ptr, pval);

    if (pan.nunits <= 30) {
        pprintf(prn, "%s = %g\n\n".ptr, _("Pooled error variance").ptr, s2);
        print_unit_variances(pan, uvar, 1, prn);
    }
}

/**
 * groupwise_hetero_test:
 * @pmod: panel model to be tested.
 * @dset: information on the (panel) data set.
 * @opt: may contain OPT_S to attach the test result
 * to @pmod, OPT_I for silent operation.
 * @prn: for printing details of iterations (or %NULL).
 *
 * Performs a Wald test for the null hypothesis that the
 * error variance is uniform across the units.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int groupwise_hetero_test(Model * pmod, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    PanelmodT pan;
    double * uvar = null;
    double s2, W = NADBL;
    int df, err = 0;

    if (pmod.ci == GretlCmdIndex.ols || (pmod.ci == GretlCmdIndex.panel && (pmod.opt & Gretlopt.OPT_F))) {
        ; /* OK for pooled or fixed effects */
    } else if (0 && pmod.ci == GretlCmdIndex.panel && (pmod.opt & Gretlopt.OPT_U)) {
        ; /* just testing! */
    } else {
        return GretlError.notimp;
    }

    gretl_error_clear();
    panelmod_init(&pan);

    err = panelmod_setup(&pan, pmod, dset, 0, Gretlopt.OPT_NONE);
    if (err) {
        return err;
    }

    uvar = cast(typeof(uvar)) malloc(pan.nunits * (*uvar).sizeof);
    if (uvar == null) {
        free(pan.unit_obs);
        return GretlError.alloc;
    }

    s2 = pmod.ess / pmod.nobs;

    unit_error_variances(uvar, pmod, &pan, &df);

    if (df >= 2) {
        W = wald_hetero_test(pmod, s2, uvar, &pan, &df);
    }

    if (na(W)) {
        err = GretlError.data;
    } else {
        double pval = chisq_cdf_comp(df, W);

        if (!(opt & Gretlopt.OPT_I)) {
            print_wald_test(W, df, pval, &pan, uvar, s2, prn);
        }

        if (opt & Gretlopt.OPT_S) {
            ModelTest *test = model_test_new(ModelTestType.groupwise);

            if (test != null) {
                model_test_set_teststat(test, TestStats.wald_chisq);
                model_test_set_dfn(test, df);
                model_test_set_value(test, W);
                model_test_set_pvalue(test, pval);
                maybe_add_test_to_model(pmod, test);
            }
        }

        record_test_result(W, pval);
    }

    free(pan.unit_obs);
    free(uvar);

    return err;
}

int print_ar_aux_model(Model * pmod, DatasetRaw * dset, int j, Prn * prn)
{
    static immutable string[3] heads = ["First differenced equation (dependent, d_y)",
        "Autoregression of residuals (dependent, uhat)",
        "Auxiliary regression including lagged residual"];
    GretlMatrixRaw * cse;
    GretlArray * S;
    int i, vi, err = 0;

    cse = gretl_matrix_alloc(pmod.ncoeff, 2);
    S = gretl_array_new(GretlType.strings, pmod.ncoeff, &err);

    if (cse == null || S == null) {
        return GretlError.alloc;
    }

    for (i=0; i<pmod.ncoeff; i++) {
        gretl_matrix_set(cse, i, 0, pmod.coeff[i]);
        gretl_matrix_set(cse, i, 1, pmod.sderr[i]);
        vi = pmod.list[i+2];
        gretl_array_set_string(S, i, dset.varname[vi], 0);
    }

    if (j >= 0) {
        pprintf(prn, "%s:\n".ptr, _(heads[j]).ptr);
    }
    print_model_from_matrices(cse, null, S, pmod.dfd, Gretlopt.OPT_NONE, prn);
    pprintf(prn, "  n = %d, R-squared = %.4f\n\n".ptr, pmod.nobs, pmod.rsq);

    gretl_matrix_free(cse);
    gretl_array_nullify_elements(S);
    gretl_array_destroy(S);

    return 0;
}

void finalize_autocorr_test(Model * pmod, ModelTest * test, Gretlopt opt, Prn * prn)
{
    if (!(opt & Gretlopt.OPT_I)) {
        if (opt & Gretlopt.OPT_Q) {
            pputc(prn, '\n');
        }
        gretl_model_test_print_direct(test, 1, prn);
    }
    if (opt & Gretlopt.OPT_S) {
        maybe_add_test_to_model(pmod, test);
    } else {
        free(test);
    }
}

/* See Wooldridge's Econometric Analysis of Cross Section
   and Panel Data (MIT Press, 2002), pages 176-7. Test
   for first order autocorrelation in the context of pooled OLS.
*/

int pooled_autocorr_test(Model * pmod, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    int quiet = (opt & Gretlopt.OPT_Q);
    int * ulist = null;
    int i, vi = 0;
    int err = 0;

    if (pmod.full_n != dset.n) {
        return GretlError.data;
    }

    err = dataset_add_NA_series(dset, 1);
    if (!err) {
        vi = dset.v - 1;
        strcpy(dset.varname[vi], "uhat");
        for (i=0; i<dset.n; i++) {
            dset.Z[vi][i] = pmod.uhat[i];
        }
    }

    if (!err) {
        ulist = gretl_list_new(pmod.ncoeff + 2);
        if (ulist == null) {
            err = GretlError.alloc;
        } else {
            for (i=1; i<=pmod.list[0]; i++) {
                ulist[i] = pmod.list[i];
            }
            /* append lagged residual to list */
            ulist[i] = laggenr(vi, 1, dset);
            if (ulist[i] < 0) {
                err = GretlError.data;
            } else {
                strcpy(dset.varname[ulist[i]], "uhat(-1)");
            }
        }
    }

    if (!err) {
        Gretlopt auxopt = Gretlopt.OPT_P | Gretlopt.OPT_R | Gretlopt.OPT_Q;
        Model tmp;

        /* estimate model including lagged residual */
        tmp = real_panel_model(ulist, dset, auxopt, null);
        err = tmp.errcode;
        if (!err && !quiet) {
            if (gretl_echo_on()) {
                pputc(prn, '\n');
            }
            print_ar_aux_model(&tmp, dset, 2, prn);
        }
        if (!err) {
            double tstat, pval;
            int df = tmp.dfd;

            i = tmp.ncoeff - 1;
            tstat = tmp.coeff[i] / tmp.sderr[i];
            pval = student_pvalue_2(df, tstat);
            record_test_result(tstat, pval);

            if ((opt & Gretlopt.OPT_S) || !(opt & Gretlopt.OPT_I)) {
                /* saving the test onto @pmod, or not silent */
                ModelTest *test = model_test_new(ModelTestType.panel_ar);

                if (test != null) {
                    model_test_set_teststat(test, TestStats.student);
                    model_test_set_value(test, tstat);
                    model_test_set_dfn(test, df);
                    model_test_set_pvalue(test, pval);
                    finalize_autocorr_test(pmod, test, opt, prn);
                }
            }
        }
        clear_model(&tmp);
    }

    free(ulist);

    return err;
}

/* See Wooldridge's Econometric Analysis of Cross Section
   and Panel Data (MIT Press, 2002), pages 282-3. Test
   for first order autocorrelation based on the residuals
   from the first-differenced model.
*/

int wooldridge_autocorr_test(Model * pmod, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    Model tmp;
    Gretlopt tmp_opt;
    int quiet = (opt & Gretlopt.OPT_Q);
    int * dlist = null;
    int i, j, vi;
    int clearit = 0;
    int err = 0;

    dlist = gretl_list_new(pmod.ncoeff + 1);
    if (dlist == null) {
        return GretlError.alloc;
    }

    dlist[0] = 0;
    for (i=1, j=1; i<=pmod.list[0] && !err; i++) {
        vi = pmod.list[i];
        if (vi != 0 && !gretl_isconst(dset.t1, dset.t2, dset.Z[vi])) {
            dlist[j] = diffgenr(vi, GretlCmdIndex.diff, dset);
            if (dlist[j] < 0) {
                err = GretlError.data;
            } else {
                j++;
                dlist[0] += 1;
            }
        } else if (i == 1) {
            gretl_errmsg_set(_("The dependent variable is constant").ptr);
            err = GretlError.data;
        }
    }

    /* aux panel models: pooled, robust, quiet */
    tmp_opt = Gretlopt.OPT_P | Gretlopt.OPT_R | Gretlopt.OPT_Q;

    if (!err) {
        /* estimate model in first-differenced form */
        tmp = real_panel_model(dlist, dset, tmp_opt, null);
        err = tmp.errcode;
        if (!err && !quiet) {
            if (gretl_echo_on()) {
                pputc(prn, '\n');
            }
            print_ar_aux_model(&tmp, dset, 0, prn);
        }
    }

    if (!err) {
        err = dataset_add_allocated_series(dset, tmp.uhat);
        tmp.uhat = null;
        clear_model(&tmp);
    }

    if (!err) {
        vi = dset.v - 1; /* the last series added */
        strcpy(dset.varname[vi], "uhat");
        dlist[0] = 2;
        dlist[1] = vi;
        dlist[2] = laggenr(vi, 1, dset);
        if (dlist[2] < 0) {
            err = GretlError.data;
        } else {
            /* regress residual on its first lag */
            tmp = real_panel_model(dlist, dset, tmp_opt, null);
            err = tmp.errcode;
            if (!err && !quiet) {
                strcpy(dset.varname[dlist[2]], "uhat(-1)");
                print_ar_aux_model(&tmp, dset, 1, prn);
            }
            clearit = 1;
        }
    }

    if (!err) {
        double c, s, F, pval;

        c = tmp.coeff[0] + 0.5;
        s = tmp.sderr[0];
        F = c * c / (s * s);
        pval = snedecor_cdf_comp(1, tmp.dfd, F);
        record_test_result(F, pval);

        if ((opt & Gretlopt.OPT_S) || !(opt & Gretlopt.OPT_I)) {
            /* saving the test onto @pmod, or not silent */
            ModelTest *test = model_test_new(ModelTestType.panel_ar);

            if (test != null) {
                model_test_set_teststat(test, TestStats.f);
                model_test_set_value(test, F);
                model_test_set_dfn(test, 1);
                model_test_set_dfd(test, tmp.dfd);
                model_test_set_pvalue(test, pval);
                finalize_autocorr_test(pmod, test, opt, prn);
            }
        }
    }

    if (clearit) {
        clear_model(&tmp);
    }
    free(dlist);

    return err;
}

int panel_autocorr_test(Model * pmod, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    int save_robust = libset_get_int(SetKey.panel_robust);
    int orig_v = dset.v;
    int err;

    libset_set_int(SetKey.panel_robust, PanelRobust.arellano); /* ? */

    if (pmod.ci == GretlCmdIndex.ols || (pmod.opt & Gretlopt.OPT_P)) {
        err = pooled_autocorr_test(pmod, dset, opt, prn);
    } else {
        err = wooldridge_autocorr_test(pmod, dset, opt, prn);
    }

    libset_set_int(SetKey.panel_robust, save_robust);
    dataset_drop_last_variables(dset, dset.v - orig_v);

    return err;
}

/* test for cross-sectional dependence */

int panel_xdepend_test(Model * pmod, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    const(double) * u;
    double rij, rsum = 0.0;
    double arsum = 0.0;
    double ssx, ssy, sxy;
    double xbar, ybar;
    int N1, N2, T, Tij;
    int i, j, t, si, sj;
    int N = 0, Nr = 0;
    int err = 0;

    if (dset.structure != TsCodes.stacked_time_series) {
        return GretlError.pdwrong;
    } else if (pmod.uhat == null) {
        return GretlError.data;
    }

    T = dset.pd;
    N1 = pmod.t1 / T;
    N2 = pmod.t2 / T;
    u = pmod.uhat;

    for (i=N1; i<N2; i++) {
        int Nj = 0;

        for (j=i+1; j<=N2; j++) {
            xbar = ybar = 0.0;
            Tij = 0;
            for (t=0; t<T; t++) {
                si = i * T + t;
                sj = j * T + t;
                if (!na(u[si]) && !na(u[sj])) {
                    Tij++;
                    xbar += u[si];
                    ybar += u[sj];
                }
            }
            if (Tij >= 2) {
                ssx = ssy = sxy = 0.0;
                xbar /= Tij;
                ybar /= Tij;
                for (t=0; t<T; t++) {
                    si = i * T + t;
                    sj = j * T + t;
                    if (!na(u[si]) && !na(u[sj])) {
                        ssx += (u[si] - xbar) * (u[si] - xbar);
                        ssy += (u[sj] - ybar) * (u[sj] - ybar);
                        sxy += (u[si] - xbar) * (u[sj] - ybar);
                    }
                }
                rij = sxy / sqrt(ssx * ssy);
                rsum += sqrt(Tij) * rij;
                arsum += fabs(rij);
                Nr++;
                Nj++;
            }
        }
        if (Nj > 0) {
            N++;
        }
    }

    if (N == 0) {
        err = GretlError.toofew;
    }

    if (!err) {
        double CD, pval;

        N = N + 1;
        CD = sqrt(2.0 / (N * (N - 1.0))) * rsum;
        pval = normal_pvalue_2(CD);

        if (!(opt & Gretlopt.OPT_I)) {
            pputs(prn, _("Pesaran CD test for cross-sectional dependence").ptr);
            pprintf(prn, "\n%s: z = %f,\n".ptr, _("Test statistic").ptr, CD);
            pprintf(prn, "%s = P(|z| > %g) = %.3g\n".ptr, _("with p-value").ptr,
                    CD, pval);
            pprintf(prn, _("Average absolute correlation = %.3f").ptr, arsum / Nr);
            pputc(prn, '\n');
        }

        if (opt & Gretlopt.OPT_S) {
            ModelTest *test = model_test_new(ModelTestType.xdepend);

            if (test != null) {
                model_test_set_teststat(test, TestStats.z);
                model_test_set_value(test, CD);
                model_test_set_pvalue(test, pval);
                maybe_add_test_to_model(pmod, test);
            }
        }

        record_test_result(CD, pval);
    }

    return err;
}

/**
 * switch_panel_orientation:
 * @dset: pointer to dataset struct.
 *
 * Reorganizes @dset, transforming from stacked cross sections to
 * stacked time series or vice versa.  If the transformation is from
 * stacked time series to stacked cross section, the dataset will
 * no longer be acceptable as a panel for gretl's purposes; it
 * may be useful for export purposes, though.
 *
 * Returns: 0 on successful completion, non-zero on error.
 */

int switch_panel_orientation(DatasetRaw * dset)
{
    int oldmode = dset.structure;
    double * tmp;
    char ** markers = null;
    double pdx;
    int T, n;
    int i, j, s, t;

    if (dset.Z == null) {
        return GretlError.nodata;
    } else if (oldmode != TsCodes.stacked_time_series &&
        oldmode != TsCodes.stacked_cross_section) {
        return GretlError.data;
    }

    tmp = cast(typeof(tmp)) malloc(dset.n * (*tmp).sizeof);
    if (tmp == null) {
        return GretlError.alloc;
    }

    if (oldmode == TsCodes.stacked_cross_section) {
        n = dset.pd;
        T = dset.n / n;
    } else {
        T = dset.pd;
        n = dset.n / T;
    }

    /* copy the data series across in transformed order */
    for (i=1; i<dset.v; i++) {
        for (t=0; t<dset.n; t++) {
            /* transcribe to tmp in original order */
            tmp[t] = dset.Z[i][t];
        }
        s = 0;
        if (oldmode == TsCodes.stacked_cross_section) {
            /* convert to stacked time-series */
            for (j=0; j<n; j++) {
                for (t=0; t<T; t++) {
                    dset.Z[i][s++] = tmp[t * n + j];
                }
            }
        } else {
            /* convert to stacked cross-sections */
            for (t=0; t<T; t++) {
                for (j=0; j<n; j++) {
                    dset.Z[i][s++] = tmp[j * T + t];
                }
            }
        }
    }

    /* rearrange observations markers if relevant */
    if (dset.S != null) {
        markers = strings_array_new_with_length(dset.n, OBSLEN);
        if (markers != null) {
            for (t=0; t<dset.n; t++) {
                strcpy(markers[t], dset.S[t]);
            }
            s = 0;
            if (oldmode == TsCodes.stacked_cross_section) {
                for (j=0; j<n; j++) {
                    for (t=0; t<T; t++) {
                        strcpy(dset.S[s++], markers[t * n + j]);
                    }
                }
            } else {
                for (t=0; t<T; t++) {
                    for (j=0; j<n; j++) {
                        strcpy(dset.S[s++], markers[j * T + t]);
                    }
                }
            }
            strings_array_free(markers, dset.n);
        } else {
            /* should we flag an error? */
            dataset_destroy_obs_markers(dset);
        }
    }

    dset.sd0 = 1.0;
    pdx = 0.1;

    /* change the datainfo setup */
    if (oldmode == TsCodes.stacked_cross_section) {
        dset.structure = TsCodes.stacked_time_series;
        dset.pd = T;
        while (T /= 10) {
            pdx *= 0.1;
        }
    } else {
        dset.structure = TsCodes.stacked_cross_section;
        dset.pd = n;
        while (n /= 10) {
            pdx *= 0.1;
        }
    }

    dset.sd0 += pdx;
    ntolabel(dset.stobs.ptr, 0, dset);
    ntolabel(dset.endobs.ptr, dset.n - 1, dset);

    free(tmp);

    return 0;
}

/**
 * panel_isconst:
 * @t1: starting observation.
 * @t2: ending observation.
 * @pd: panel time-series length.
 * @x: data series to examine.
 * @bygroup: use 1 to check for constancy across groups,
 * 0 for constancy across time.
 *
 * Check whether series @x is constant (either over time or
 * across groups) in a panel dataset. Missing values are
 * ignored.
 *
 * Returns: 1 if the variable is constant, otherwise 0.
 */

int panel_isconst(int t1, int t2, int pd, const(double) * x, int bygroup)
{
    double x0 = NADBL;
    int t, ret = 1;

    if (bygroup) {
        /* check for variation across groups */
        int tref;

        for (t=t1; t<=t2 && ret; t++) {
            if (na(x[t])) {
                continue;
            }
            tref = t - pd;
            while (tref >= t1) {
                /* check last obs for the same period */
                x0 = x[t-pd];
                if (na(x0)) {
                    tref -= pd;
                } else {
                    if (x[t] != x0) {
                        ret = 0;
                    }
                    break;
                }
            }
        }
    } else {
        /* check for time-variation */
        int u, ubak = -1;

        for (t=t1; t<=t2 && ret; t++) {
            if (na(x[t])) {
                continue;
            }
            u = t / pd;
            if (u == ubak) {
                /* same group */
                if (!na(x0) && x[t] != x0) {
                    ret = 0;
                }
            } else {
                /* starting a new group */
                x0 = x[t];
                ubak = u;
            }
        }
    }

    return ret;
}

int varying_vars_list(const(DatasetRaw) * dset, PanelmodT * pan)
{
    int i, j, k, t, bigt;
    int err = 0;

    pan.vlist = gretl_list_new(pan.pooled.list[0]);
    if (pan.vlist == null) {
        return GretlError.alloc;
    }

    pan.vlist[0] = 1;
    pan.vlist[1] = pan.pooled.list[1];

    k = 2;

    for (j=2; j<=pan.pooled.list[0]; j++) {
        int vj = pan.pooled.list[j];
        int varies = 0;
        double x0;

        if (vj == 0) {
            /* const */
            pan.vlist[k++] = 0;
            pan.vlist[0] += 1;
            continue;
        }

        for (i=0; i<pan.nunits && !varies; i++) {
            if (pan.unit_obs[i] == 0) {
                continue;
            }
            x0 = NADBL;
            for (t=0; t<pan.T && !varies; t++) {
                bigt = panel_index(i, t);
                if (panel_missing(pan, bigt)) {
                    continue;
                }
                if (na(x0)) {
                    x0 = dset.Z[vj][bigt];
                } else if (dset.Z[vj][bigt] != x0) {
                    varies = 1;
                }
            }
        }

        if (varies) {
            pan.vlist[k++] = vj;
            pan.vlist[0] += 1;
        } else {
            fprintf(stderr, "Variable %d '%s' is time-invariant\n".ptr,
                    vj, dset.varname[vj]);
        }
    }


    return err;
}

/* apparatus for setting panel structure by reading two variables,
   indexing the cross-sectional unit and the time period
   respectively
*/


struct SPoint {
    int obsnum;
    double val1;
    double val2;
}

struct Sorter {
    int sortvar1;
    int sortvar2;
    SPoint *points;
}

void sorter_fill_points(Sorter * s, const(double) ** Z, int n)
{
    int t;

    for (t=0; t<n; t++) {
        s.points[t].obsnum = t;
        s.points[t].val1 = Z[s.sortvar1][t];
        s.points[t].val2 = Z[s.sortvar2][t];
    }
}

extern (C) int compare_obs(const(void) * a, const(void) * b)
{
    SPoint *pa = cast(SPoint *) a;
    SPoint *pb = cast(SPoint *) b;
    int ret;

    ret = (pa.val1 > pb.val1) - (pa.val1 < pb.val1);
    if (ret == 0) {
        ret = (pa.val2 > pb.val2) - (pa.val2 < pb.val2);
    }

    if (ret == 0) {
        /* mark an error: got a duplicated value */
        pa.obsnum = pb.obsnum = -1;
    }

    return ret;
}

int check_indices(Sorter * s, int n)
{
    int i;

    for (i=0; i<n; i++) {
        if (s.points[i].obsnum < 0) {
            gretl_errmsg_sprintf(_("Error: unit %g, period %g: duplicated observation").ptr,
                                 s.points[i].val1, s.points[i].val2);
            return GretlError.data;
        }
    }

    return 0;
}

int panel_is_sorted(DatasetRaw * dset, int uv, int tv)
{
    const(double) * unit = dset.Z[uv];
    const(double) * period = dset.Z[tv];
    int t, ret = 1;

    for (t=1; t<dset.n && ret; t++) {
        if (unit[t] < unit[t-1]) {
            /* the unit variable must be non-decreasing */
            ret = 0;
        } else if (unit[t] == unit[t-1] && period[t] <= period[t-1]) {
            /* the period variable must be increasing within the
               observations for each unit */
            ret = 0;
        }
    }

    return ret;
}

int panel_data_sort_by(DatasetRaw * dset, int uv, int tv, int * ustrs)
{
    int n = dset.n;
    char ** S = null;
    double * tmp = null;
    int i, t;
    Sorter s;
    int err = 0;

    tmp = cast(typeof(tmp)) malloc(n * (*tmp).sizeof);
    if (tmp == null) {
        return GretlError.alloc;
    }

    s.points = cast(typeof(s.points)) malloc(n * (*s.points).sizeof);
    if (s.points == null) {
        free(tmp);
        return GretlError.alloc;
    }

    if (dset.S != null) {
        S = strings_array_new_with_length(n, OBSLEN);
        if (S == null) {
            free(tmp);
            free(s.points);
            return GretlError.alloc;
        }
    }

    s.sortvar1 = uv;
    s.sortvar2 = tv;

    sorter_fill_points(&s, cast(const(double) **) dset.Z, n);

    qsort(cast(void *) s.points, cast(size_t) n, (*s.points).sizeof,
          &compare_obs);

    err = check_indices(&s, n);
    if (err) {
        goto bailout;
    }

    for (i=1; i<dset.v; i++) {
        for (t=0; t<n; t++) {
            tmp[t] = dset.Z[i][t];
        }
        for (t=0; t<n; t++) {
            dset.Z[i][t] = tmp[s.points[t].obsnum];
        }
    }

    if (S != null) {
        *ustrs = 1;
        for (t=0; t<n; t++) {
            strcpy(S[t], dset.S[t]);
            if (S[t][0] == '\0') {
                *ustrs = 0;
            }
        }
        for (t=0; t<n; t++) {
            strcpy(dset.S[t], S[s.points[t].obsnum]);
        }
        for (t=1; t<n && *ustrs; t++) {
            if (dset.Z[uv][t] == dset.Z[uv][t-1] &&
                strcmp(dset.S[t], dset.S[t-1])) {
                *ustrs = 0;
            }
        }
        strings_array_free(S, n);
    }

 bailout:

    free(s.points);
    free(tmp);

    return err;
}

/* Given the variables coding for panel unit and time period,
   construct parallel arrays of int that code for the same
   information, but are zero-based and consecutive.
*/

int normalize_uid_tid(const(double) * tid, int T, const(DatasetRaw) * dset, int uv, int tv, int ** pnuid, int ** pntid)
{
    int * nuid = null;
    int * ntid = null;
    int ui, tmin;
    int i, t;

    nuid = cast(typeof(nuid)) malloc(dset.n * (*nuid).sizeof);
    ntid = cast(typeof(ntid)) malloc(dset.n * (*ntid).sizeof);
    if (nuid == null || ntid == null) {
        free(nuid);
        free(ntid);
        return GretlError.alloc;
    }

    ui = tmin = 0;

    for (i=0; i<dset.n; i++) {
        if (i > 0 && dset.Z[uv][i] > dset.Z[uv][i-1]) {
            tmin = 0;
            ui++;
        } else if (i > 0) {
            tmin++;
        }
        nuid[i] = ui;
        for (t=tmin; t<T; t++) {
            if (dset.Z[tv][i] == tid[t]) {
                ntid[i] = t;
                break;
            }
        }
    }

    *pnuid = nuid;
    *pntid = ntid;

    return 0;
}

/* Handle the case where a sub-sampled panel dataset has been padded
   with rows containing NAs to recreate a (nominally) balanced panel
   structure. We have to get rid of the padding before we try
   restoring the full data range. Note that while this function
   "shrinks" @dset in the sense of reducing the series length as
   recorded in dset->n, it does not "physically" shrink the array
   members of dset->Z.
*/

int undo_panel_padding(DatasetRaw * dset)
{
    char * mask = dset.padmask;
    double * Zi;
    char ** S = null;
    int padded_n = dset.n; /* current series length */
    int real_n = dset.n; /* target series length */
    int i, t;
    int err = 0;

    for (t=0; t<dset.n; t++) {
        if (mask[t]) {
            real_n--;
        }
    }

    fprintf(stderr, "undo_panel_padding: padded n*T = %d, original dset->n = %d\n".ptr,
            padded_n, real_n);

    if (real_n == padded_n) {
        fprintf(stderr, "strange, couldn't find any padding!\n".ptr);
        return GretlError.data;
    }

    /* temporary holder for shorter series */
    Zi = cast(typeof(Zi)) malloc(real_n * (*Zi).sizeof);

    if (Zi == null) {
        return GretlError.alloc;
    }

    if (dset.S != null) {
        /* holder for shorter obs labels array */
        S = strings_array_new_with_length(real_n, OBSLEN);
    }

    if (!err) {
        /* write non-padding rows from dset->Z (and dset->S, if present)
           into the right places in the reduced-size arrays */
        int s;

        for (i=0; i<dset.v; i++) {
            s = 0;
            for (t=0; t<padded_n; t++) {
                if (!mask[t]) {
                   Zi[s] = dset.Z[i][t];
                   if (i == 0 && S != null) {
                       strcpy(S[s], dset.S[t]);
                   }
                   s++;
                }
            }
            /* copy data back to dset->Z */
            memcpy(dset.Z[i], Zi, real_n * (*Zi).sizeof);
        }

        if (dset.S != null && S != null) {
            strings_array_free(dset.S, padded_n);
            dset.S = S;
        }
    }

    free(Zi);

    dset.n = real_n;
    dset.t2 = dset.n - 1;
    free(dset.padmask);
    dset.padmask = null;

    return err;
}

int pad_panel_dataset(const(double) * uid, int uv, int nunits, const(double) * tid, int tv, int nperiods, DatasetRaw * dset, int ustrs, char * mask)
{
    double ** bigZ = null;
    char ** S = null;
    int * nuid = null;
    int * ntid = null;
    int big_n, n_orig = dset.n;
    int i, s, t;
    int err = 0;

    big_n = nunits * nperiods;

    fprintf(stderr, "pad_panel_dataset: n*T = %d*%d = %d but dset->n = %d\n".ptr,
            nunits, nperiods, big_n, n_orig);

    err = normalize_uid_tid(tid, nperiods, dset, uv, tv,
                            &nuid, &ntid);
    if (err) {
        return err;
    }

    bigZ = doubles_array_new(dset.v, big_n);
    if (bigZ == null) {
        err = GretlError.alloc;
    } else {
        dset.n = big_n;
        dset.t2 = dset.n - 1;
        for (i=0; i<dset.v; i++) {
            for (t=0; t<big_n; t++) {
                bigZ[i][t] = (i == 0)? 1.0 : NADBL;
            }
        }
    }

    if (!err && dset.S != null && ustrs) {
        S = strings_array_new_with_length(dset.n, OBSLEN);
    }

    if (!err) {
        int j, buv = 0, btv = 0;

        if (mask != null) {
            for (t=0; t<big_n; t++) {
                mask[t] = 1;
            }
        }

        /* write rows from original Z into the right places in bigZ */
        for (t=0; t<n_orig; t++) {
            s = nuid[t] * nperiods + ntid[t];
            for (i=1; i<dset.v; i++) {
                bigZ[i][s] = dset.Z[i][t];
            }
            if (S != null) {
                strcpy(S[s], dset.S[t]);
            }
            if (mask != null) {
                /* record (non-)padding in @mask */
                mask[s] = 0;
            }
        }

        /* where are uv and tv in bigZ? */
        j = 1;
        for (i=1; i<dset.v; i++) {
            if (i == uv) {
                buv = j;
            } else if (i == tv) {
                btv = j;
            }
            j++;
        }

        /* complete the index info in slots uv and tv */
        i = j = 0;
        for (t=0; t<dset.n; t++) {
            if (t > 0 && t % nperiods == 0) {
                i++;
                j = 0;
            }
            bigZ[buv][t] = uid[i];
            bigZ[btv][t] = tid[j++];
        }

        /* swap the padded arrays into Z */
        for (i=0; i<dset.v; i++) {
            free(dset.Z[i]);
            dset.Z[i] = bigZ[i];
        }

        free(bigZ);
    }

    free(nuid);
    free(ntid);

    if (!err && dset.S != null) {
        /* expand the obs (unit) marker strings appropriately */
        if (S == null) {
            strings_array_free(dset.S, n_orig);
            dset.S = null;
            dset.markers = DatasetMarkerType.no_markers;
        } else {
            char[OBSLEN] si;
            int j;

            for (i=0; i<nunits; i++) {
                t = i * nperiods;
                for (j=0; j<nperiods; j++) {
                    if (S[t][0] != '\0') {
                        strcpy(si.ptr, S[t]);
                        break;
                    }
                    t++;
                }
                t = i * nperiods;
                for (j=0; j<nperiods; j++) {
                    strcpy(S[t++], si.ptr);
                }
            }
            strings_array_free(dset.S, n_orig);
            dset.S = S;
        }
    }

    return err;
}

int check_index_values(const(double) * x, int n)
{
    int i;

    for (i=1; i<n; i++) {
        if (x[i] < 0) {
            return GretlError.data;
        } else if (na(x[i])) {
            return GretlError.missdata;
        }
    }

    return 0;
}

int uv_tv_from_varnames(const(char) * uname, const(char) * tname, const(DatasetRaw) * dset, int * uv, int * tv)
{
    *uv = current_series_index(dset, uname);
    if (*uv <= 0) {
        return GretlError.data;
    }

    *tv = current_series_index(dset, tname);
    if (*tv <= 0) {
        return GretlError.data;
    }

    return 0;
}


void finalize_panel_datainfo(DatasetRaw * dset, int nperiods)
{
    int pdp = nperiods;
    int den = 10;

    while ((pdp = pdp / 10) != 0) {
        den *= 10;
    }
    dset.structure = TsCodes.stacked_time_series;
    dset.pd = nperiods;
    dset.sd0 = 1.0 + 1.0 / den;
    ntolabel(dset.stobs.ptr, 0, dset);
    ntolabel(dset.endobs.ptr, dset.n - 1, dset);
}

int check_full_dataset()
{
    DatasetRaw * fset = fetch_full_dataset();

    if (!dataset_is_panel(fset)) {
        string msg = "You cannot use the --panel-vars option with the setobs command when\n\n* the dataset is currently sub-sampled and\n* the full dataset is not a panel.\n\nIf you first structure your full dataset as a panel, you can then\ndo what you are trying to do.";

        gretl_errmsg_set(_(msg).ptr);
        return GretlError.data;
    }

    return 0;
}

/**
 * set_panel_structure_from_vars:
 * @uv: index of series uniquely identifying units/groups.
 * @tv: index of series uniquely identifying time periods.
 * @dset: pointer to dataset.
 *
 * Sets the panel structure of @dset based on the information
 * in the series with indices @uv and @tv, if possible.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int set_panel_structure_from_vars(int uv, int tv, DatasetRaw * dset)
{
    int subsampled, presorted;
    double * uid = null;
    double * tid = null;
    char * mask = null;
    int n = dset.n;
    int totmiss = 0;
    int fulln = 0;
    int nunits = 0;
    int nperiods = 0;
    int ustrs = 0;
    int err = 0;

    subsampled = complex_subsampled();

    if (subsampled) {
        err = check_full_dataset();
        if (err) {
            return err;
        }
    }


    /* start by making sorted copies of the unit and time
       variables and counting the unique values of each
    */

    uid = copyvec(dset.Z[uv], n);
    tid = copyvec(dset.Z[tv], n);

    if (uid == null || tid == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    qsort(uid, n, (*uid).sizeof, &gretl_compare_doubles);
    nunits = count_distinct_values(uid, dset.n);

    qsort(tid, n, (*tid).sizeof, &gretl_compare_doubles);
    nperiods = count_distinct_values(tid, dset.n);

    /* check that we have a possible panel structure */

    if (nunits == 1 || nperiods == 1 ||
        nunits == n || nperiods == n ||
        n > nunits * nperiods) {
        fprintf(stderr, "Dataset does not have a panel structure\n".ptr);
        fprintf(stderr, " (nunits=%d, nperiods=%d, n=%d)\n".ptr, nunits, nperiods, n);
        err = GretlError.data;
        goto bailout;
    }

    /* figure how many panel-data rows are implied, and how many
       "padding" rows are needed to complete the structure
    */

    fulln = nunits * nperiods;
    totmiss = fulln - dset.n;


    /* determine if the data rows are already in sort order by unit,
       and by period for each unit */

    presorted = panel_is_sorted(dset, uv, tv);


    if (!presorted) {
        if (subsampled) {
            /* We can't re-order the rows in a subsampled dataset, or the
               data will get scrambled when restoration of the full dataset
               is attempted.
            */
            gretl_errmsg_set(_("Sorry, can't do this with a sub-sampled dataset").ptr);
            err = GretlError.data;
        } else {
            /* sort full dataset by unit and period */
            err = panel_data_sort_by(dset, uv, tv, &ustrs);
        }
    }

    if (!err && totmiss > 0) {
        /* establish a nominally balanced panel */
        mask = cast(typeof(mask)) malloc(fulln);
        if (mask == null) {
            err = GretlError.alloc;
        } else {
            rearrange_id_array(uid, nunits, n);
            rearrange_id_array(tid, nperiods, n);
            err = pad_panel_dataset(uid, uv, nunits,
                                    tid, tv, nperiods,
                                    dset, ustrs,
                                    mask);
        }
    }

    if (!err) {
        finalize_panel_datainfo(dset, nperiods);
    }

 bailout:

    free(uid);
    free(tid);

    if (!err && complex_subsampled()) {
        dset.padmask = mask;
    } else {
        free(mask);
    }

    return err;
}

int set_panel_structure_from_varnames(const(char) * uname, const(char) * tname, DatasetRaw * dset)
{
    int n = dset.n;
    int uv = 0, tv = 0;
    int err = 0;

    err = uv_tv_from_varnames(uname, tname, dset, &uv, &tv);

    if (!err) {
        err = check_index_values(dset.Z[uv], n);
    }

    if (!err) {
        err = check_index_values(dset.Z[tv], n);
    }

    if (!err) {
        err = set_panel_structure_from_vars(uv, tv, dset);
    }

    return err;
}

int group_uniqueness_check(char ** S, int n)
{
    int i, j;

    for (i=0; i<n-1; i++) {
        for (j=i+1; j<n; j++) {
            if (!strcmp(S[i], S[j])) {
                gretl_errmsg_sprintf(_("The string '%s' is given for two or more groups").ptr, S[i]);
                return GretlError.data;
            }
        }
    }

    return 0;
}

char ** group_names_from_array(const(char) * aname, int n_groups, int * err)
{
    GretlArray * A = get_array_by_name(aname);
    char ** S = null;
    int ns = 0;

    if (A == null) {
        *err = GretlError.data;
    } else {
        char ** AS = gretl_array_get_strings(A, &ns);

        if (ns != n_groups) {
            *err = GretlError.data;
        } else {
            S = strings_array_dup(AS, ns);
            if (S == null) {
                *err = GretlError.alloc;
            }
        }
    }

    return S;
}

/* usable_groups_series: to be usable as the basis for panel
   groups strings, a series must have values that are (a)
   constant within-group and (b) unique across groups. Note
   that we only get to this test if the series in question
   has already been found to be string-valued.
*/

int usable_groups_series(DatasetRaw * dset, int v, const(char) * vname)
{
    const(double) * x = dset.Z[v];
    int i, j, g = 0;
    int ok = 1;

    if (na(x[0])) {
        return 0;
    }

    for (i=1; i<dset.n && ok; i++) {
        if (na(x[i])) {
            ok = 0;
        } else if (i % dset.pd == 0) {
            /* starting a new group: x-value must not repeat a
               prior group's value */
            for (j=0; j<=g; j++) {
                if (x[i] == x[j*dset.pd]) {
                    gretl_errmsg_sprintf(_("%s: values are not unique per group").ptr,
                                         vname);
                    ok = 0;
                }
            }
            g++;
        } else if (x[i] != x[i-1]) {
            gretl_errmsg_sprintf(_("%s: is not constant within group").ptr, vname);
            ok = 0;
        }
    }

    return ok;
}

int maybe_use_strval_series(DatasetRaw * dset, const(char) * vname, int ng)
{
    int v = current_series_index(dset, vname);
    SeriesTable * st = null;
    int ns = 0;
    int err = 0;

    if (v <= 0) {
        return GretlError.invarg;
    }

    st = series_get_string_table(dset, v);
    if (st == null) {
        gretl_errmsg_sprintf(_("The series %s is not string-valued").ptr, vname);
        return GretlError.invarg;
    }

    series_table_get_strings(st, &ns);


    if (ns < ng) {
        gretl_errmsg_sprintf(_("The series %s holds %d strings but %d are needed").ptr, vname, ns, ng);
        err = GretlError.invarg;
    } else {
        /* note: don't mess with numerical values, just check them */
        if (usable_groups_series(dset, v, vname)) {
            set_panel_groups_name(dset, vname);
        } else {
            err = GretlError.invarg;
        }
    }

    return err;
}

/* Enable construction of a string-valued series holding a name for
   each panel group/unit. The name of this series is set on the
   'pangrps' member of @dset, and the names are then used in panel
   graphs where appropriate.

   @vname should contain the name of a series, either an existing one
   (which may be overwritten) or a new one to create; and @grpnames
   should be (1) a string literal, or (2) the name of a string variable
   (which in either case should hold N space-separated strings, where N
   is the number of panel groups), or (3) the name of an array variable
   holding N strings.
*/

int set_panel_group_strings(const(char) * vname, const(char) * grpnames, DatasetRaw * dset)
{
    const(char) * namestr = null;
    char ** S = null;
    int ng = dset.n / dset.pd;
    int v, orig_v = dset.v;
    int err = 0;

    if (vname == null || *vname == '\0') {
        return GretlError.data;
    }

    if (grpnames == null || *grpnames == '\0') {
        /* We just got the name of a series? That may be
           OK if it's string-valued and has a suitable
           set of values.
        */
        return maybe_use_strval_series(dset, vname, ng);
    }

    if (strchr(grpnames, ' ')) {
        /* group names as string literal */
        namestr = grpnames;
    } else if (gretl_is_string(grpnames)) {
        /* group names in a single string variable */
        namestr = get_string_by_name(grpnames);
    } else {
        /* try for array of strings */
        S = group_names_from_array(grpnames, ng, &err);
    }

    if (namestr != null) {
        /* we must obtain @S by splitting */
        int ns = 0;

        if (strchr(namestr, '"') != null) {
            S = gretl_string_split_quoted(namestr, &ns, null, &err);
        } else {
            S = gretl_string_split(namestr, &ns, " \n\t");
        }
        if (!err && S == null) {
            err = GretlError.alloc;
        }
        if (!err && ns != ng) {
            /* FIXME subsampled case? */
            fprintf(stderr, "Got %d strings but there are %d groups\n".ptr,
                    ns, ng);
            strings_array_free(S, ns);
            S = null;
            err = GretlError.data;
        }
    }

    if (!err && S != null) {
        err = group_uniqueness_check(S, ng);
    }

    if (!err) {
        v = current_series_index(dset, vname);
        if (v < 0) {
            /* we need to add a series */
            char * gen = gretl_strdup_printf("%s", vname);

            err = generate(gen, dset, GretlType.series, Gretlopt.OPT_Q, null);
            if (!err) {
                v = dset.v - 1;
            }
            free(gen);
        }
    }

    if (!err) {
        SeriesTable * st = series_table_new(S, ng, &err);

        if (!err) {
            int i, g = 0;

            series_attach_string_table(dset, v, st);
            for (i=0; i<dset.n; i++) {
                if (i % dset.pd == 0) {
                    g++;
                }
                dset.Z[v][i] = g;
            }
        }
    }

    if (err) {
        if (S != null) {
            strings_array_free(S, ng);
        }
        if (dset.v > orig_v) {
            dataset_drop_last_variables(dset, dset.v - orig_v);
        }
    } else {
        set_panel_groups_name(dset, vname);
    }

    return err;
}

/* utility functions */

int find_time_var(const(DatasetRaw) * dset)
{
    static immutable string[5] tnames = ["year",
        "Year",
        "period",
        "Period",
    "time"];
    int i, v;

    for (i=0; i<tnames.length; i++) {
        v = series_index(dset, tnames[i].ptr);
        if (v < dset.v) {
            return v;
        }
    }

    return 0;
}

int guess_panel_structure(double ** Z, DatasetRaw * dset)
{
    int ret, v = find_time_var(dset);

    if (v == 0) {
        ret = 0; /* can't guess */
    } else if (floateq(Z[v][0], Z[v][1])) {
        /* "year" or whatever is same for first two obs */
        dset.structure = TsCodes.stacked_cross_section;
        ret = TsCodes.stacked_cross_section;
    } else {
        dset.structure = TsCodes.stacked_time_series;
        ret = TsCodes.stacked_time_series;
    }

    return ret;
}

int balanced_panel(const(DatasetRaw) * dset)
{
    int n = dset.t2 - dset.t1 + 1;
    int ret = 0;

    if (n % dset.pd == 0) {
        char[OBSLEN] unit;
        char[OBSLEN] period;

        if (sscanf(dset.endobs.ptr, "%[^:]:%s".ptr, unit.ptr, period.ptr) == 2) {
            if (atoi(period.ptr) == dset.pd) {
                ret = 1;
            }
        }
    }

    return ret;
}

int may_be_time_name(const(char) * vname)
{
    char[VNAMELEN] test;

    strcpy(test.ptr, vname);
    gretl_lower(test.ptr);

    return strcmp(test.ptr, "year") == 0 ||
        strcmp(test.ptr, "time") == 0 ||
    strcmp(test.ptr, "period") == 0;
}

/* See if the panel dataset contains a variable that plausibly
   represents the year of the observations. We look for series labeled
   "year", "time" or "period" (in a case insensitive comparison) and,
   if found, check that the series has the same sequence of values for
   each individual, and the same increment between successive
   time-series observations.  (The increment does not have to be 1, to
   accommodate, e.g., quinquennial or decennial observations.)

   Return the ID number of a suitable variable, or 0 if there's nothing
   suitable.

   This may be used in setting the x-axis tic marks in a panel time-
   series plot.
*/

int plausible_panel_time_var(const(DatasetRaw) * dset)
{
    int i, t, ret = 0;

    for (i=1; i<dset.v && ret==0; i++) {
        if (may_be_time_name(dset.varname[i])) {
            const(double) * x = dset.Z[i];
            int val0 = cast(int) x[0];
            int incr0 = cast(int) (x[1] - x[0]);
            int ok = 1;

            for (t=0; t<dset.n && ok; t++) {
                if (na(x[t]) || x[t] < 0) {
                    ok = 0;
                } else if (t > 0 && t % dset.pd == 0) {
                    if (x[t] != val0) {
                        ok = 0;
                    }
                } else if (t > 1 && x[t] - x[t-1] != incr0) {
                    ok = 0;
                }
            }
            if (ok) {
                ret = i;
            }
        }
    }

    return ret;
}

/* Tests whether the series with index number @v codes for the
   time dimension of a panel dataset. If so, returns 1 and
   records the minimum value and (constant) increment of the
   index; otherwise returns 0.
*/

int is_panel_time_var(const(DatasetRaw) * dset, int v, int tmax, int * minval, int * incr)
{
    int t, ret = 0;

    *minval = 0;
    *incr = 0;

    if (may_be_time_name(dset.varname[v])) {
        const(double) * x = dset.Z[v];
        int val0 = cast(int) x[0];
        int incr0 = cast(int) x[1] - cast(int) x[0];
        int ok = 1;

        for (t=0; t<tmax && ok; t++) {
            if (na(x[t]) || x[t] < 0 || x[t] != floor(x[t])) {
                ok = 0;
            } else if (t > 0 && t % dset.pd == 0) {
                if (x[t] != val0) {
                    ok = 0;
                }
            } else if (t > 1 && x[t] - x[t-1] != incr0) {
                ok = 0;
            }
        }
        if (ok) {
            *minval = val0;
            *incr = incr0;
            ret = 1;
        }
    }

    return ret;
}

int is_panel_unit_var(const(DatasetRaw) * dset, int v, int tmax)
{
    const(double) * x = dset.Z[v];
    int t, ret = 1;

    if (x[0] != 1) {
        return 0;
    }

    for (t=1; t<tmax && ret; t++) {
        if (na(x[t]) || x[t] < 0 || x[t] != floor(x[t])) {
            ret = 0;
        } else if (t % dset.pd == 0) {
            if (x[t] != x[t-1] + 1) {
                ret = 0;
            }
        } else if (x[t] != x[t-1]) {
            ret = 0;
        }
    }

    return ret;
}

/* FIXME: this does not yet handle the dropping of instruments */

int * dpanel_list_omit(const(Model) * orig, const(int) * drop, int * err)
{
    const(int) * old = orig.list;
    int * new_ = gretl_list_copy(old);
    int sep = 0;
    int i, j;

    if (new_ == null) {
        *err = GretlError.alloc;
        return null;
    }

    for (i=2; i<=new_[0]; i++) {
        if (new_[i] == LISTSEP) {
            sep++;
        }
        if (sep == 1) {
            for (j=1; j<=drop[0]; j++) {
                if (drop[j] == new_[i]) {
                    gretl_list_delete_at_pos(new_, i--);
                }
            }
        }
    }


    return new_;
}

/**
 * panel_list_omit:
 * @orig: list specifying original panel model.
 * @drop: list of variables to be omitted.
 * @err: pointer to receive error code.
 *
 * Creates a new panel regression list, by first reconstructing
 * the regression specification from @orig then deleting from
 * the reconstruction the variables found in @drop.
 *
 * Returns: the new, reduced list or %NULL on error.
 */

int * panel_list_omit(const(Model) * orig, const(int) * drop, int * err)
{
    int * newlist = null;
    int i;

    if (orig.ci == GretlCmdIndex.dpanel) {
        return dpanel_list_omit(orig, drop, err);
    }

    /* sorry, can't drop the constant */
    if (drop != null) {
        int cpos = in_gretl_list(drop, 0);

        if (cpos >= 2) {
            gretl_errmsg_set(_("Panel models must include an intercept").ptr);
            *err = GretlError.data;
            return null;
        }
    }

    if (orig.opt & Gretlopt.OPT_F) {
        int * panlist;

        /* fixed-effects lists have the constant removed,
           so we need to put it back first
        */
        panlist = gretl_list_new(orig.list[0] + 1);
        if (panlist != null) {
            panlist[1] = orig.list[1];
            panlist[2] = 0;
            for (i=3; i<=panlist[0]; i++) {
                panlist[i] = orig.list[i-1];
            }
            if (drop == null) {
                newlist = gretl_list_omit_last(panlist, err);
            } else {
                newlist = gretl_list_omit(panlist, drop, 2, err);
            }
            free(panlist);
        }
    } else if (drop == null) {
        newlist = gretl_list_omit_last(orig.list, err);
    } else {
        newlist = gretl_list_omit(orig.list, drop, 2, err);
    }

    return newlist;
}

/* FIXME doesn't handle adding instruments */

int * dpanel_list_add(const(Model) * orig, const(int) * add, int * err)
{
    const(int) * old = orig.list;
    int * new_ = gretl_list_copy(old);
    int sep = 0, pos = old[0] + 1;
    int i;

    if (new_ == null) {
        *err = GretlError.alloc;
        return null;
    }

    for (i=2; i<=old[0]; i++) {
        if (old[i] == LISTSEP) {
            sep++;
            if (sep == 2) {
                pos = i - 1;
            }
        }
    }

    gretl_list_insert_list(&new_, add, pos);
    if (new_ == null) {
        *err = GretlError.alloc;
        return null;
    }


    return new_;
}

/**
 * panel_list_add:
 * @orig: list specifying original panel model.
 * @add: list of variables to be added.
 * @err: pointer to receive error code.
 *
 * Creates a new panel regression list, by first reconstructing
 * the regression specification from @orig then adding to
 * the reconstruction the variables found in @add.
 *
 * Returns: the new, augmented list or %NULL on error.
 */

int * panel_list_add(const(Model) * orig, const(int) * add, int * err)
{
    if (orig.ci == GretlCmdIndex.dpanel) {
        return dpanel_list_add(orig, add, err);
    } else {
        return gretl_list_add(orig.list, add, err);
    }
}

/* calculate the within and between standard deviations for a given
   variable in a panel data set */

int panel_variance_info(const(double) * x, const(DatasetRaw) * dset, double xbar, double * psw, double * psb)
{
    double sw = 0.0, sb = 0.0;
    double xibar, d;
    int effn, effnT;
    int n, T, nT, Ti;
    int i, t, s;

    if (!dataset_is_panel(dset)) {
        return GretlError.pdwrong;
    }

    nT = dset.t2 - dset.t1 + 1;
    T = dset.pd;
    n = nT / T;

    effn = 0;
    effnT = 0;

    for (i=0; i<n; i++) {
        Ti = 0;
        xibar = 0.0;
        for (t=0; t<T; t++) {
            s = dset.t1 + i * T + t;
            if (!na(x[s])) {
                Ti++;
                xibar += x[s];
            }
        }
        if (Ti > 1) {
            xibar /= Ti;
            for (t=0; t<T; t++) {
                s = dset.t1 + i * T + t;
                if (!na(x[s])) {
                    d = x[s] - xibar;
                    sw += d * d;
                }
            }
        }
        if (Ti > 0) {
            /* is this right for singleton observations? */
            d = xibar - xbar;
            sb += d * d;
            effn++;
            effnT += Ti;
        }
    }

    if (effn > 1) {
        /* between variance: \sum_{i=1}^N (\bar{x}_i - \bar{x})^2
           \times (N-1)^{-1}
        */
        sb /= (effn - 1);
        sb = sqrt(sb);
    } else {
        sb = NADBL;
    }

    if (effnT - effn > 0) {
        /* within variance: \sum_{i=1}^N \sum_{t=1}^T (x_{it} - \bar{x}_i)^2
           \times (\sum_{i=1}^N T_i - N)^{-1}
        */
        sw /= (effnT - effn);
        sw = sqrt(sw);
    } else {
        sw = NADBL;
    }

    *psw = sw;
    *psb = sb;

    return 0;
}

int series_is_group_invariant(const(DatasetRaw) * dset, int v)
{
    const(double) * x = dset.Z[v];
    int T = dset.pd;
    int N = dset.n / T;
    int i, t, s;

    for (i=1; i<N; i++) {
        for (t=0; t<T; t++) {
            s = t + i * T;
            if (x[s] != x[t]) {
                return 0;
            }
        }
    }

    return 1;
}

int panel_padding_rows(const(DatasetRaw) * dset)
{
    int missrow, nmiss = 0;
    int i, t;

    for (t=dset.t1; t<=dset.t2; t++) {
        missrow = 1;
        for (i=1; i<dset.v; i++) {
            if (!na(dset.Z[i][t])) {
                missrow = 0;
                break;
            }
        }
        if (missrow) {
            nmiss++;
        }
    }

    return nmiss;
}

/* Given auxiliary time-series info attached to the panel
   dataset @pdset, transcribe it to the time-series dataset
   @tsset.
*/

int time_series_from_panel(DatasetRaw * tset, const(DatasetRaw) * pset)
{
    int err = 0;

    if (pset.panel_pd == 0) {
        /* the panel time dimension is not set */
        tset.structure = TsCodes.special_time_series;
        tset.pd = 1;
        return 0;
    }

    tset.structure = TsCodes.time_series;
    tset.pd = pset.panel_pd;
    tset.sd0 = pset.panel_sd0;

    if (tset.pd == 1) {
        sprintf(tset.stobs.ptr, "%d".ptr, cast(int) tset.sd0);
    } else if (tset.pd == 4 || tset.pd == 12) {
        double dyr = floor(tset.sd0);
        double dp = tset.sd0 - dyr;
        int yr, p;

        yr = cast(int) dyr;
        dp *= (tset.pd == 4)? 10 : 100;
        p = cast(int) nearbyint(dp);
        if (yr > 0 && yr < 9999 && p > 0) {
            if (tset.pd == 4 && p < 4) {
                sprintf(tset.stobs.ptr, "%d:%d".ptr, yr, p);
            } else if (p < 12) {
                sprintf(tset.stobs.ptr, "%d:%02d".ptr, yr, p);
            }
        } else {
            err = 1;
        }
    } else if (calendar_data(tset)) {
        calendar_date_string(tset.stobs.ptr, 0, tset);
    }

    return err;
}

/* Given an annual, quarterly or monthly date string @s in gretl's
   standard format (YYYY, YYYY:Q or YYYY:MM), convert to a zero-based
   observation index relative to a series of frequency @pd and length
   @n starting at @t0. (We'll also accept YYYY:M for monthly data with
   a single-digit month.)

   @t0 is assumed to be given in a form that makes the calculation
   easy, namely the number of periods since the start of the year 0.
   For annual data this is just the year; for quarterly its 4 * year
   plus initial quarter and for monthly 12 * year plus initial month.

   If the resulting index is out of bounds, -1 is returned.
*/

int obs_index_from_aqm(const(char) * s, int pd, int t0, int n)
{
    const(char) * digits = "0123456789";
    int len = cast(int) strlen(s);
    int len_ok = 0;
    int t = 0;

    if (pd == 1) {
    len_ok = (len == 4); /* e.g. 2025 */
    } else if (pd == 4) {
    len_ok = (len == 6); /* e.g. 2025:1 */
    } else if (pd == 12) {
    len_ok = (len == 6 || len == 7); /* 2025:1 or 2025:01 */
    }

    if (!len_ok || strspn(s, digits) != 4) {
        return -1;
    }

    if (pd == 1) {
        t = atoi(s);
    } else if (s[4] != ':' || strspn(s + 5, digits) != len - 5) {
        t = -1;
    } else {
        int sub = atoi(s + 5);

        if (sub > pd) {
            t = -1;
        } else {
            t = pd * atoi(s) + sub;
        }
    }

    if (t >= 0) {
        t -= t0;
        if (t < 0 || t >= n) {
            t = -1;
        }
    }

    return t;
}

/* Scan @dset for a string-valued series that's suitable for defining
   panel group names: its number of string values must equal the
   number of groups, and its value must be constant within each
   group's time series. In addition, if @maxlen > 0 the strings must
   be no longer than @maxlen UTF-8 characters (we're looking to use
   these strings for plotting and we don't want them too long).

   We return the ID number of the first series to match these
   criteria, or 0 if we can't find any such series.
*/

int usable_group_names_series_id(const(DatasetRaw) * dset, int maxlen)
{
    char ** S;
    int N = dset.n / dset.pd;
    int j, t, s, ns, ok;
    int i, ret = 0;

    for (i=1; i<dset.v; i++) {
        if (!is_string_valued(dset, i)) {
            continue;
        }
        S = series_get_string_vals(dset, i, &ns, 0);
        if (ns < N) {
            continue;
        }
        ok = 1;
        for (j=0, s=0; j<N && ok; j++) {
            if (maxlen > 0 && g_utf8_strlen(S[j], -1) > maxlen) {
                /* too long to be usable for plots? */
                ok = 0;
            }
            for (t=0; t<dset.pd && ok; t++) {
                if (t > 0 && dset.Z[i][s] != dset.Z[i][s-1]) {
                    ok = 0;
                } else if (j > 0 && t == 0 &&
                           dset.Z[i][s] == dset.Z[i][s-1]) {
                    ok = 0;
                }
                s++;
            }
        }
        if (ok) {
            ret = i;
            break;
        }
    }

    return ret;
}
