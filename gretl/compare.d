/* Port of lib/src/compare.c: comparison and diagnostic tests on estimated models: omit/add variables (F, Wald and
   LM forms), Chow and CUSUM tests, RESET, non-linearity, heteroskedasticity, autocorrelation, ARCH and other
   auxiliary regressions. The header compare.h is compareh.d. */

module gretl.compare;


import core.stdc.ctype : isdigit, isspace;
import core.stdc.errno : errno, ERANGE;
import core.stdc.math : ceil, exp, fabs, floor, fmod, isinf, isnan, log, log10, pow, sqrt;
import core.stdc.stdio : FILE, fclose, fopen, fputs, fprintf, printf, puts, sprintf, snprintf, sscanf, stderr;
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

import gretl.gretl_intl : gretl_pop_c_numeric_locale, gretl_push_c_numeric_locale;
import gretl.compareh : ModelAuxCode, SeriesSaveCode;
import gretl.dataio : ntolabel;
import gretl.dataset : dataset_add_series, dataset_drop_last_variables, series_get_lag, series_get_parent_id, series_set_label, set_dataset_is_changed;
import gretl.dataseth : dataset_is_panel, dataset_is_time_series;
import gretl.describe : gretl_minmax, gretl_stddev, ljung_box, set_last_result_data;
import gretl.discrete : count_model, duration_model, heckit_model, logistic_model;
import gretl.estimate : ar1_model, ar_model, arch_model, augment_regression_list, dpd_model, hsk_model, ivreg, lad_model, lsq, makevcv, panel_model;
import gretl.genmain : current_series_index;
import gretl.genmainh : ModelDataIndex;
import gretl.graphingh : GptFlags, PlotType;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_errors : errmsg, gretl_errmsg_set, gretl_errmsg_sprintf, gretl_error_clear;
import gretl.gretl_list : gretl_list_add, gretl_list_append_term, gretl_list_copy, gretl_list_delete_at_pos, gretl_list_diff_new, gretl_list_new, gretl_list_omit, gretl_lists_join_with_separator, in_gretl_list, list_members_replaced;
import gretl.gretl_listh : LISTSEP;
import gretl.gretl_model : ModelTest, arma_model_get_n_arma_coeffs, clear_model, coeff_pval, command_ok_for_model, exact_fit_check, get_h0_string_for_test, get_model_count, gretl_is_regular_panel_model, gretl_model_free, gretl_model_get_cluster_vname, gretl_model_get_data, gretl_model_get_double, gretl_model_get_int, gretl_model_get_list, gretl_model_get_matrix, gretl_model_get_param_number, gretl_model_init, gretl_model_new, gretl_model_set_double, gretl_model_set_int, gretl_model_set_list_as_data, impose_model_smpl, maybe_add_test_to_model, model_count_minus, model_test_new, model_test_set_allocated_param, model_test_set_dfd, model_test_set_dfn, model_test_set_opt, model_test_set_order, model_test_set_param, model_test_set_pvalue, model_test_set_teststat, model_test_set_value;
import gretl.gretl_modelh : COUNT_MODEL, ModelTestType;
import gretl.gretl_panel : BFN_panel_DW_pvalue, balanced_panel, panel_DW_pval_ok, panel_autocorr_test, panel_diagnostics, panel_list_add, panel_list_omit;
import gretl.gretl_prn : Prn, pprintf, pputc, pputs;
import gretl.gretl_utils : gretl_double_from_string, record_QLR_test_result, record_test_result;
import gretl.matrix_extra : gretl_coeff_vector_from_model, gretl_vcv_matrix_from_model;
import gretl.missing : model_has_missing_obs, model_missing, set_reference_missmask_from_model;
import gretl.modelprint : printmodel;
import gretl.options : get_optval_string, incompatible_options, set_optval_double, set_optval_string, transcribe_option_flags;
import gretl.plotspech : GpKeySpec;
import gretl.plugins : get_plugin_function;
import gretl.pvalues : chisq_cdf_comp, snedecor_cdf_comp, student_pvalue_2;
import gretl.strutils : gretl_trunc, make_varname_unique;
import gretl.transforms : laggenr;
import gretl.tsls : ivreg_list_add, ivreg_list_omit;
import gretl.plotspec : print_keypos_string;
import gretl.graphing : finalize_plot_input_file, gnuplot_graph_wanted, open_plot_input_file, reset_effective_plot_ci, set_effective_plot_ci;
import gretl.uservar : get_list_by_name;
import gretl.genfuncs : gretl_plotx;

/* unported dependencies: per-file declarations from the 2026c headers */
extern (C) {
    Model logit_probit(int * list, DatasetRaw * dset, int ci, Gretlopt opt, Prn * prn);
    Model quantreg_driver(const(char) * param, const(int) * list, DatasetRaw * dset, Gretlopt opt, Prn * prn);
    Model tobit_driver(const(int) * list, DatasetRaw * dset, Gretlopt opt, Prn * prn);
}

struct Compare {
    int ci; /* command index: ADD or OMIT */
    int model_id; /* ID of model under test */
    int model_ci; /* estimator code for the model under test */
    int dfn; /* numerator degrees of freedom */
    int dfd; /* denominator degrees of freedom */
    double test; /* test statistic */
    double pval; /* p-value for test statistic */
    int stat; /* code identifying the test statistic */
    double F; /* F test statistic */
    double X2; /* Wald chi-square test statistic */
    double LR; /* likelihood ratio test statistic */
    int score; /* number of info stats showing improvement */
    int robust; /* = 1 when robust vcv is in use, else 0 */
    const(int) * testvars; /* list of variables added or omitted */
};

/* Given a @list of series, check them against the independent variables
   included in @pmod and construct a mask with 1s in positions where
   there is a match, 0s otherwise.  If @list is NULL, match all
   regressors other than the constant.
*/

char * mask_from_test_list(const(int) * list, const(Model) * pmod, int * err)
{
    char * mask;
    int nmask = 0;
    int i, j, vj;
    int nx;

    if (pmod.xlist == null) {
    /* we can't do this */
    gretl_errmsg_set(_("The model for testing lacks an xlist").ptr);
    *err = GretlError.data;
    return null;
    }

    nx = pmod.xlist[0];
    mask = cast(typeof(mask)) calloc(pmod.ncoeff, 1);
    if (mask == null) {
    *err = GretlError.alloc;
    return null;
    }

    if (list != null) {
    /* matching members of @list */
    for (j=1; j<=list[0]; j++) {
        vj = list[j];
        for (i=1; i<=nx; i++) {
        if (pmod.xlist[i] == vj) {
            mask[i-1] = 1;
            nmask++;
        }
        }
    }
    } else {
    /* implicitly matching all but const */
    for (i=1; i<=nx; i++) {
        if (pmod.xlist[i] != 0) {
        mask[i-1] = 1;
        }
    }
    }

    if (list != null && nmask != list[0]) {
    fprintf(stderr, "mask from list: list[0] = %d but nmask = %d\n".ptr,
        list[0], nmask);
    *err = GretlError.data;
    }

    return mask;
}

/* Wald (chi-square and/or F) test for a set of zero restrictions on the
   parameters of a given model, based on the covariance matrix of the
   unrestricted model. Suitable for use where the original model is
   estimated by FGLS or IV.  The zero restrictions can be specified by
   either (a) a @list of series to test, or (b) a @mask with 1s in the
   test positions, 0s otherwise. If neither @list not @mask is non-NULL
   the test is for omission of all regressors except the constant.

   The value returned on success is by default a chi-square statistic,
   but if @opt contains OPT_F it's the F form.
*/

double wald_omit_test(const(int) * list, char * mask, Model * pmod, Gretlopt opt, int * err)
{
    GretlMatrixRaw * C = null;
    GretlVector * b = null;
    int free_mask = 0;
    int df = 0;
    double X2 = NADBL;

    if (mask == null) {
    mask = mask_from_test_list(list, pmod, err);
    if (*err) {
        free(mask);
        return NADBL;
    }
    free_mask = 1;
    }

    if (!*err) {
    C = gretl_vcv_matrix_from_model(pmod, mask, err);
    }

    if (!*err) {
    b = gretl_coeff_vector_from_model(pmod, mask, err);
    }

    if (!*err) {
    /* use "left division" to avoid explicit inversion of @C */
    GretlMatrixRaw * tmp, mx;

    mx = gretl_matrix_alloc(1, 1);
    tmp = gretl_matrix_divide(C, b, GretlMatrixMod.none, err);
    if (tmp != null) {
        gretl_matrix_multiply_mod(b, GretlMatrixMod.transpose,
                      tmp, GretlMatrixMod.none,
                      mx, GretlMatrixMod.none);
        X2 = mx.val[0];
        gretl_matrix_free(tmp);
        gretl_matrix_free(mx);
    }
    df = gretl_vector_get_length(b);
    }

    if (free_mask) {
    free(mask);
    }
    gretl_matrix_free(C);
    gretl_matrix_free(b);

    return (opt & Gretlopt.OPT_F)? X2 / df : X2;
}

/**
 * wald_omit_F:
 * @list: list of variables to omit, or NULL.
 * @pmod: model to be tested.
 *
 * Simple form of Wald F-test for omission of variables.  If @list
 * is non-NULL, do the test for the omission of the variables in
 * @list from the model @pmod.  Otherwise test for omission of
 * all variables in @pmod except for the constant.
 *
 * Returns: Calculated F-value, or #NADBL on failure.
 */

double wald_omit_F(const(int) * list, Model * pmod)
{
    int err = 0;

    return wald_omit_test(list, null, pmod, Gretlopt.OPT_F, &err);
}

/**
 * wald_omit_chisq:
 * @list: list of variables to omit, or NULL.
 * @pmod: model to be tested.
 *
 * Simple form of Wald chi-square for omission of variables.  If @list
 * is non-NULL, do the test for the omission of the variables in
 * @list from the model @pmod.  Otherwise test for omission of
 * all variables in @pmod except for the constant.
 *
 * Returns: Calculated chi-square value, or #NADBL on failure.
 */

double wald_omit_chisq(const(int) * list, Model * pmod)
{
    int err = 0;

    return wald_omit_test(list, null, pmod, Gretlopt.OPT_X, &err);
}

int add_diffvars_to_test(ModelTest * test, const(int) * list, const(DatasetRaw) * dset)
{
    char * vnames;
    int i, len = 0;
    int err = 0;

    for (i=1; i<=list[0]; i++) {
    len += strlen(dset.varname[list[i]]) + 1;
    }

    vnames = cast(typeof(vnames)) malloc(len);

    if (vnames == null) {
    err = 1;
    } else {
    *vnames = '\0';
    for (i=1; i<=list[0]; i++) {
        strcat(vnames, dset.varname[list[i]]);
        if (i < list[0]) {
        strcat(vnames, " ");
        }
    }
    model_test_set_allocated_param(test, vnames);
    }

    return err;
}

void print_add_omit_null(const(int) * list, const(DatasetRaw) * dset, Gretlopt opt, Prn * prn)
{
    int nl = list == null ? 0 : list[0];

    if (nl == 1 && !(opt & Gretlopt.OPT_S)) {
    /* single parameter, not in system context */
    pputs(prn, "\n  ");
    pprintf(prn, _("Null hypothesis: the regression parameter is zero for %s").ptr,
        dset.varname[list[1]]);
    pputc(prn, '\n');
    } else if (nl == 0) {
    /* VAR omit specials */
    if ((opt & Gretlopt.OPT_E) && (opt & Gretlopt.OPT_T)) {
        pprintf(prn, "\n  %s: %s\n".ptr, _("Null hypothesis").ptr,
            _("no seasonal effects or trend").ptr);
    } else if (opt & Gretlopt.OPT_E) {
        pprintf(prn, "\n  %s: %s\n".ptr, _("Null hypothesis").ptr,
            _("no seasonal effects").ptr);
    } else if (opt & Gretlopt.OPT_T) {
        pprintf(prn, "\n  %s: %s\n".ptr, _("Null hypothesis").ptr,
            _("no trend").ptr);
    }
    } else {
    /* all other cases */
    const(char) * vname;
    int i, nc = 0;

    pputs(prn, _("\n  Null hypothesis: the regression parameters are zero for the variables\n").ptr);

    pputs(prn, "    ");
    for (i=1; i<=list[0]; i++) {
        vname = dset.varname[list[i]];
        nc += strlen(vname) + 2;
        pprintf(prn, "%s".ptr, vname);
        if (i < list[0]) {
        if (nc > 60) {
            pputs(prn, ",\n    ");
            nc = 0;
        } else {
            pputs(prn, ", ");
        }
        }
    }
    pputc(prn, '\n');
    if (opt & Gretlopt.OPT_E) {
        /* seasonals */
        pprintf(prn, "    %s\n".ptr, _("seasonal dummies").ptr);
    }
    if (opt & Gretlopt.OPT_T) {
        /* trend */
        pputs(prn, "    time\n");
    }
    }
}

/* add/omit: are we printing a revised model following the
   test results? */

int printing_revised_model(Gretlopt opt)
{
    if (opt & Gretlopt.OPT_Q) {
    /* --quiet: no */
    return 0;
    } else if (opt & (Gretlopt.OPT_L | Gretlopt.OPT_W)) {
    /* --lm or --wald test variants: no */
    return 0;
    } else {
    /* then yes */
    return 1;
    }
}

/* print the results from an add or omit test */

void print_compare(Compare * cmp, const(DatasetRaw) * dset, Gretlopt opt, Prn * prn)
{
    if (opt & Gretlopt.OPT_Q) {
    if (cmp.model_id > 0) {
        pprintf(prn, _("Test on Model %d:").ptr, cmp.model_id);
        pputc(prn, '\n');
    }
    } else if (cmp.model_id > 0) {
    if (opt & Gretlopt.OPT_A) {
        /* --auto-omit: add vertical space */
        pputc(prn, '\n');
    }
    pprintf(prn, _("Test on Model %d:").ptr, cmp.model_id);
    pputc(prn, '\n');
    }

    print_add_omit_null(cmp.testvars, dset, Gretlopt.OPT_NONE, prn);

    if (cmp.stat == TestStats.wald_chisq) {
    pprintf(prn, "  %s: %s(%d) = %g, %s %g\n".ptr,  _("Wald test").ptr,
        _("Chi-square").ptr, cmp.dfn, cmp.test,
        _("p-value").ptr, cmp.pval);
    if (!na(cmp.LR)) {
        /* auxiliary LR test? */
        double pval = chisq_cdf_comp(cmp.dfn, cmp.LR);

        pprintf(prn, "  (%s: %s(%d) = %g, %s %g)\n".ptr,  _("LR test").ptr,
            _("Chi-square").ptr, cmp.dfn, cmp.LR,
            _("p-value").ptr, pval);
    } else if (!na(cmp.F)) {
        /* alternate F-form for Wald? */
        double pval = snedecor_cdf_comp(cmp.dfn, cmp.dfd, cmp.F);

        pprintf(prn, "  (%s: F(%d, %d) = %g, %s %g)\n".ptr,
            _("F-form").ptr, cmp.dfn, cmp.dfd, cmp.F,
            _("p-value").ptr, pval);
    }
    } else if (cmp.stat == TestStats.lr) {
    /* note: unused */
    pprintf(prn, "  %s: %s(%d) = %g, %s %g\n".ptr,  _("LR test").ptr,
        _("Chi-square").ptr, cmp.dfn, cmp.test,
        _("p-value").ptr, cmp.pval);
    } else if (cmp.stat == TestStats.f) {
    pprintf(prn, "  %s: %s(%d, %d) = %g, %s %g\n".ptr, _("Test statistic").ptr,
        (cmp.robust)? _("Robust F").ptr : "F".ptr,
        cmp.dfn, cmp.dfd, cmp.test,
        _("p-value").ptr, cmp.pval);
    } else if (cmp.stat == TestStats.lm) {
    pprintf(prn, "  %s: %s(%d) = %g, %s %g\n".ptr,  _("LM test").ptr,
        _("Chi-square").ptr, cmp.dfn, cmp.test,
        _("p-value").ptr, cmp.pval);
    }

    if (cmp.score >= 0) {
    pputs(prn, "  ");
    if (cmp.ci == GretlCmdIndex.add) {
        pprintf(prn, _("Adding variables improved %d of %d information criteria.\n").ptr, cmp.score, InfoCriterion.max);
    } else {
        pprintf(prn, _("Omitting variables improved %d of %d information criteria.\n").ptr, cmp.score, InfoCriterion.max);
    }
    }

    if (!printing_revised_model(opt)) {
    pputc(prn, '\n');
    }
}

/* try attaching the add or omit test to @pmod */

void add_omit_attach(Compare * cmp, Model * pmod, const(DatasetRaw) * dset)
{
    int tcode = (cmp.ci == GretlCmdIndex.omit)? ModelTestType.omit : ModelTestType.add;
    ModelTest *mtest = model_test_new(cast(ModelTestType) tcode);

    if (mtest != null) {
    model_test_set_teststat(mtest, cast(ubyte) cmp.stat);
    model_test_set_dfn(mtest, cmp.dfn);
    if (cmp.stat == TestStats.f) {
        model_test_set_dfd(mtest, cmp.dfd);
    }
    model_test_set_value(mtest, cmp.test);
    model_test_set_pvalue(mtest, cmp.pval);
    add_diffvars_to_test(mtest, cmp.testvars, dset);
    maybe_add_test_to_model(pmod, mtest);
    }
}

/* add a record of improvement/degeneration of the
   model selection criteria, for non-quiet printing:
   @pmodA is the original model, and @pmodB the
   restricted (OMIT) or augmented (ADD) variant
*/

void maybe_add_info_score(Compare * cmp, const(Model) * pmodA, const(Model) * pmodB)
{
    if (pmodA != null && pmodB != null) {
    int i;

    cmp.score = 0;
    for (i=0; i<InfoCriterion.max; i++) {
        if (na(pmodB.criterion[i]) || na(pmodA.criterion[i])) {
        cmp.score = -1;
        break;
        } else if (pmodB.criterion[i] < pmodA.criterion[i]) {
        cmp.score++;
        }
    }
    }
}

void compare_init(Compare * cmp, Model * pmod, int ci, const(int) * testvars)
{
    cmp.ci = ci;
    cmp.stat = 0;
    cmp.test = cmp.pval = NADBL;
    cmp.F = cmp.X2 = cmp.LR = NADBL;
    cmp.score = -1;
    cmp.robust = 0;
    cmp.dfn = testvars[0];
    cmp.testvars = testvars;
    cmp.model_ci = pmod.ci;
    cmp.model_id = pmod.ID;
}

int add_or_omit_compare(Model * pmodA, Model * pmodB, const(int) * testvars, const(DatasetRaw) * dset, int ci, Gretlopt opt, Prn * prn)
{
    Compare cmp;
    Model * umod, rmod;
    int iv_special = 0;
    int err = 0;

    compare_init(&cmp, pmodA, ci, testvars);

    if (ci == GretlCmdIndex.omit) {
    /* the unrestricted model is the original one, 'A' */
    umod = pmodA;
    rmod = pmodB;
    } else {
    /* add: the unrestricted model is the new one, 'B' */
    umod = pmodB;
    rmod = pmodA;
    }

    if (opt & Gretlopt.OPT_B) {
    /* ivreg special: A and B have different
       instruments, not just different regressors
    */
    iv_special = 1;
    cmp.model_id = -1;
    } else if (opt & Gretlopt.OPT_L) {
    cmp.model_id = -1;
    } else if (umod != null && rmod != null) {
    cmp.dfn = umod.ncoeff - rmod.ncoeff;
    }

    if (pmodA.opt & Gretlopt.OPT_R) {
    cmp.robust = 1;
    }

    if (pmodA.ci == GretlCmdIndex.dpanel) {
    /* FIXME: plus some other cases? */
    opt |= Gretlopt.OPT_X;
    }

    cmp.dfd = umod.dfd;

    if (opt & Gretlopt.OPT_L) {
    /* "add" with LM option: the auxiliary regression
       results are in pmodB */
    cmp.stat = TestStats.lm;
    cmp.test = pmodB.nobs * pmodB.rsq;
    cmp.pval = chisq_cdf_comp(cmp.dfn, cmp.test);
    } else if (cmp.model_ci == GretlCmdIndex.ols && !cmp.robust &&
           rmod != null && umod != null &&
           rmod.nobs == umod.nobs &&
           !(opt & Gretlopt.OPT_X)) {
    /* plain OLS with both the restricted and unrestricted
       models available: base F-test on sums of squared
       residuals
    */
    cmp.stat = TestStats.f;
    cmp.test = ((rmod.ess - umod.ess) / umod.ess) *
        cmp.dfd / cmp.dfn;
    cmp.pval = snedecor_cdf_comp(cmp.dfn, cmp.dfd, cmp.test);
    } else if (opt & Gretlopt.OPT_X) {
    /* chi-square form of Wald test is requested */
    cmp.stat = TestStats.wald_chisq;
    cmp.X2 = wald_omit_test(testvars, null, umod, Gretlopt.OPT_X, &err);
    if (!err) {
        cmp.test = cmp.X2;
        cmp.pval = chisq_cdf_comp(cmp.dfn, cmp.test);
    }
    } else {
    /* F-form of Wald test */
    cmp.stat = TestStats.f;
    cmp.F = wald_omit_test(testvars, null, umod, Gretlopt.OPT_F, &err);
    if (!err) {
        cmp.test = cmp.F;
        cmp.pval = snedecor_cdf_comp(cmp.dfn, cmp.dfd, cmp.test);
    }
    }

    /* auxiliary LR test? */
    if (umod != null && rmod != null && !(opt & Gretlopt.OPT_L)) {
    if (!na(umod.lnL) && !na(rmod.lnL)) {
        cmp.LR = 2.0 * (umod.lnL - rmod.lnL);
    }
    }

    if (!err && (na(cmp.test) || na(cmp.pval))) {
    err = GretlError.data;
    }

    if (!err) {
    record_test_result(cmp.test, cmp.pval);
    if (opt & Gretlopt.OPT_S) {
        /* attach test to model */
        add_omit_attach(&cmp, pmodA, dset);
    }
    if (!(opt & Gretlopt.OPT_I)) {
        /* not --silent */
        if (!(opt & Gretlopt.OPT_Q) && cmp.stat != TestStats.lm &&
        !iv_special) {
        maybe_add_info_score(&cmp, pmodA, pmodB);
        }
        print_compare(&cmp, dset, opt, prn);
    }
    }

    return err;
}

int get_extra_var(const(Model) * pmod)
{
    if (COUNT_MODEL(cast(GretlCmdIndex) pmod.ci)) {
    return gretl_model_get_int(pmod, "offset_var");
    } else if (pmod.ci == GretlCmdIndex.duration) {
    return gretl_model_get_int(pmod, "cens_var");
    } else {
    return 0;
    }
}

/* reconstitute full varlist for WLS, AR, count and
   duation models */

int * full_model_list(const(Model) * pmod, const(int) * inlist)
{
    int * flist = null;

    if (pmod.ci == GretlCmdIndex.ar) {
    /* cobble together arlist and @inlist */
    flist = gretl_lists_join_with_separator(pmod.arinfo.arlist,
                        inlist);
    } else {
    int i, len = inlist[0];

    if (pmod.ci == GretlCmdIndex.wls) {
        /* prepend the weight variable */
        len += 1;
    } else if (COUNT_MODEL(cast(GretlCmdIndex) pmod.ci) || pmod.ci == GretlCmdIndex.duration) {
        /* append list separator and offset or censoring var */
        len += 2;
    }

    flist = gretl_list_new(len);
    if (flist == null) {
        return null;
    }

    if (pmod.ci == GretlCmdIndex.wls) {
        flist[1] = pmod.nwt;
        for (i=1; i<=inlist[0]; i++) {
        flist[i+1] = inlist[i];
        }
    } else if (COUNT_MODEL(cast(GretlCmdIndex) pmod.ci) || pmod.ci == GretlCmdIndex.duration) {
        int extra = get_extra_var(pmod);

        for (i=1; i<=inlist[0]; i++) {
        flist[i] = inlist[i];
        }
        flist[flist[0]-1] = LISTSEP;
        flist[flist[0]] = extra;
    }
    }

    return flist;
}

Gretlopt retrieve_dpanel_opts(const(Model) * pmod)
{
    Gretlopt opt = Gretlopt.OPT_NONE;

    if (pmod.opt & Gretlopt.OPT_D) {
    /* --dpdstyle */
    opt |= Gretlopt.OPT_D;
    }

    if (pmod.opt & Gretlopt.OPT_L) {
    /* --system (include levels) */
    opt |= Gretlopt.OPT_L;
    }

    if (gretl_model_get_int(pmod, "asy")) {
    opt |= Gretlopt.OPT_A;
    }

    if (gretl_model_get_int(pmod, "step") == 2) {
    opt |=Gretlopt.OPT_T;
    }

    return opt;
}

int obs_diff_ok(const(Model) * m_old, const(Model) * m_new)
{
    int tdiff, ndiff = m_new.nobs - m_old.nobs;

    if (m_old.ci == GretlCmdIndex.ar1) {
    return 0;
    }

    if (ndiff > 0) {
    tdiff = (m_new.t2 - m_new.t1) - (m_old.t2 - m_old.t1);
    if (ndiff == tdiff) {
        return 1;
    }
    }

    return 0;
}

/* This function is used for "add" and "omit", when we are estimating
   an augmented or reduced version of the original model. It's also
   used in the special case of calculation of the p-value for the
   Durbin-Watson statistic, which requires re-estimation of the
   original specification. In these cases we need to ensure
   comparability, which means we have to retrieve any relevant options
   from the original model and re-apply them.

   The flag @auto_omit indicates that this function is being called
   in the context of the "omit" command with the --auto option.
*/

Model replicate_estimator(const(Model) * orig, int * list, DatasetRaw * dset, Gretlopt myopt, Prn * prn)
{
    Model rep;
    const(char) * param = null;
    const(int) * laglist = null;
    const(char) * cname = null;
    int * full_list = null;
    char[32] altparm = 0;
    int mc = get_model_count();
    int repci = orig.ci;
    int cv, order = 0;
    int first = 1;

    gretl_model_init(&rep, dset);

    /* recreate options and auxiliary vars, if required */

    transcribe_option_flags(&myopt, orig.opt, Gretlopt.OPT_D | Gretlopt.OPT_J | Gretlopt.OPT_R);

    cname = gretl_model_get_cluster_vname(orig);
    if (cname != null) {
    /* FIXME second cluster-var name? */
    cv = current_series_index(dset, cname);
    if (cv > 0) {
        myopt |= Gretlopt.OPT_C;
        if (orig.ci == GretlCmdIndex.ols && gretl_model_get_int(orig, "pooled")) {
        /* transcribe the cluster spec to "panel", since
           we'll be calling panel_model() below
        */
        set_optval_string(GretlCmdIndex.panel, Gretlopt.OPT_C, cname);
        } else {
        set_optval_string(cast(GretlCmdIndex) orig.ci, Gretlopt.OPT_C, cname);
        }
    }
    }

    if (orig.ci == GretlCmdIndex.ar1) {
    if (orig.opt & Gretlopt.OPT_H) {
        myopt |= Gretlopt.OPT_H;
        if (gretl_model_get_int(orig, "no-corc")) {
        myopt |= Gretlopt.OPT_B;
        }
    } else if (orig.opt & Gretlopt.OPT_P) {
        myopt |= Gretlopt.OPT_P;
    }
    } else if (orig.ci == GretlCmdIndex.wls || orig.ci == GretlCmdIndex.ar || get_extra_var(orig)) {
    full_list = full_model_list(orig, list);
    if (full_list == null) {
        rep.errcode = GretlError.alloc;
    } else {
        list = full_list;
    }
    } else if (orig.ci == GretlCmdIndex.hsk) {
    if (gretl_model_get_int(orig, "no-squares")) {
        myopt |= Gretlopt.OPT_N;
    }
    } else if (orig.ci == GretlCmdIndex.dpanel) {
    param = cast(typeof(param)) gretl_model_get_data(orig, "istr");
    laglist = gretl_model_get_list(orig, "ylags");
    myopt |= retrieve_dpanel_opts(orig);
    } else if (orig.ci == GretlCmdIndex.arch) {
    order = gretl_model_get_int(orig, "arch_order");
    } else if (orig.ci == GretlCmdIndex.logit || orig.ci == GretlCmdIndex.probit) {
    if (orig.opt & Gretlopt.OPT_P) {
        /* p-values, not slopes, selected */
        myopt |= Gretlopt.OPT_P;
    } else if (gretl_model_get_int(orig, "ordered")) {
        myopt |= Gretlopt.OPT_D;
    } else if (gretl_model_get_int(orig, "multinom")) {
        myopt |= Gretlopt.OPT_M;
    } else if (orig.ci == GretlCmdIndex.probit && (orig.opt & Gretlopt.OPT_E)) {
        /* random effects */
        int qp = gretl_model_get_int(orig, "quadpoints");

        myopt |= (Gretlopt.OPT_E | Gretlopt.OPT_G);
        set_optval_double(GretlCmdIndex.probit, Gretlopt.OPT_G, qp);
        }
        if (orig.opt & Gretlopt.OPT_N) {
            myopt |= Gretlopt.OPT_N;
        }
    } else if (orig.ci == GretlCmdIndex.panel) {
    if (gretl_model_get_int(orig, "pooled")) {
        /* pooled OLS */
        myopt |= Gretlopt.OPT_P;
    } else if (orig.opt & Gretlopt.OPT_U) {
        /* random effects */
        myopt |= Gretlopt.OPT_U;
    } else if (orig.opt & Gretlopt.OPT_H) {
        /* unit weights */
        myopt |= Gretlopt.OPT_H;
        if (gretl_model_get_int(orig, "iters")) {
        myopt |= Gretlopt.OPT_I;
        }
    }
    } else if (orig.ci == GretlCmdIndex.lad && gretl_model_get_int(orig, "rq")) {
    double x;

    x = gretl_model_get_double(orig, "tau");
    sprintf(altparm.ptr, "%g".ptr, x);
    if (gretl_model_get_int(orig, "rq_nid")) {
        myopt |= Gretlopt.OPT_R;
    }
    x = gretl_model_get_double(orig, "rq_alpha");
    if (!na(x)) {
        myopt |= Gretlopt.OPT_I;
        set_optval_double(GretlCmdIndex.quantreg, Gretlopt.OPT_I, x);
    }
    } else if (orig.ci == GretlCmdIndex.tobit) {
    double x;

    x = gretl_model_get_double(orig, "llimit");
    if (!na(x)) {
        myopt |= Gretlopt.OPT_L;
        set_optval_double(GretlCmdIndex.tobit, Gretlopt.OPT_L, x);
    }
    x = gretl_model_get_double(orig, "rlimit");
    if (!na(x)) {
        myopt |= Gretlopt.OPT_M;
        set_optval_double(GretlCmdIndex.tobit, Gretlopt.OPT_M, x);
    }
    } else if (orig.ci == GretlCmdIndex.ivreg) {
    transcribe_option_flags(&myopt, orig.opt, Gretlopt.OPT_L | Gretlopt.OPT_G);
    }

    if (rep.errcode) {
    goto bailout;
    }

 try_again:

    switch (orig.ci) {

    case GretlCmdIndex.ar:
    rep = ar_model(list, dset, myopt, null);
    break;
    case GretlCmdIndex.ar1:
    rep = ar1_model(list, dset, myopt, null);
    break;
    case GretlCmdIndex.dpanel:
    rep = dpd_model(list, laglist, param, dset, myopt, prn);
    break;
    case GretlCmdIndex.arch:
    rep = arch_model(list, order, dset, myopt);
    break;
    case GretlCmdIndex.logit:
    case GretlCmdIndex.probit:
    rep = logit_probit(list, dset, orig.ci, myopt, null);
    break;
    case GretlCmdIndex.tobit:
    rep = tobit_driver(list, dset, myopt, null);
    break;
    case GretlCmdIndex.lad:
    if (gretl_model_get_int(orig, "rq")) {
        rep = quantreg_driver(altparm.ptr, list, dset, myopt, null);
    } else {
        rep = lad_model(list, dset, myopt);
    }
    break;
    case GretlCmdIndex.poisson:
    case GretlCmdIndex.negbin:
    rep = count_model(list, orig.ci, dset, myopt, null);
    break;
    case GretlCmdIndex.duration:
    rep = duration_model(list, dset, myopt, null);
    break;
    case GretlCmdIndex.heckit:
    rep = heckit_model(list, dset, myopt, null);
    break;
    case GretlCmdIndex.ivreg:
    rep = ivreg(list, dset, myopt);
    break;
    case GretlCmdIndex.logistic:
    {
        double lmax = gretl_model_get_double(orig, "lmax");

        rep = logistic_model(list, lmax, dset, myopt);
    }
    break;
    case GretlCmdIndex.panel:
    rep = panel_model(list, dset, myopt, prn);
    break;
    case GretlCmdIndex.hsk:
    rep = hsk_model(list, dset, myopt);
    break;
    default:
    /* handles OLS, WLS, etc. */
    if (gretl_model_get_int(orig, "pooled")) {
        myopt |= Gretlopt.OPT_P;
        rep = panel_model(list, dset, myopt, prn);
    } else {
        rep = lsq(list, dset, cast(GretlCmdIndex) repci, myopt);
    }
    break;
    }

    /* check that we got the same sample as the original */
    if (!rep.errcode && rep.nobs != orig.nobs) {
    if (first && obs_diff_ok(orig, &rep)) {
        dset.t1 = orig.t1;
        dset.t2 = orig.t2;
        clear_model(&rep);
        first = 0;
        goto try_again;
    } else {
            gretl_errmsg_sprintf(_("replicate_estimator() failed:\n %d observations used versus %d for the original model").ptr,
                                 rep.nobs, orig.nobs);
        rep.errcode = GretlError.data;
    }
    }

 bailout:

    free(full_list);

    /* if the model count went up for an aux regression,
       bring it back down */
    if ((myopt & Gretlopt.OPT_A) && get_model_count() > mc) {
    model_count_minus(null);
    }

    return rep;
}

int add_residual_to_dataset(Model * pmod, DatasetRaw * dset)
{
    int err = 0;

    if (dataset_add_series(dset, 1)) {
    err = GretlError.alloc;
    } else {
    int t, v = dset.v - 1;

    for (t=0; t<dset.n; t++) {
        dset.Z[v][t] = pmod.uhat[t];
    }

    strcpy(dset.varname[v], "uhat");
    series_set_label(dset, v, _("residual").ptr);
    }

    return err;
}

void nonlin_test_header(int code, Prn * prn)
{
    pputc(prn, '\n');
    if (code == ModelAuxCode.sq) {
    pputs(prn, _("Non-linearity test (squared terms)").ptr);
    } else {
    pputs(prn, _("Non-linearity test (log terms)").ptr);
    }
    pputs(prn, "\n\n");
}

int real_nonlinearity_test(Model * pmod, int * list, DatasetRaw * dset, int aux_code, Gretlopt opt, Prn * prn)
{
    Model aux;
    int err;

    err = add_residual_to_dataset(pmod, dset);
    if (err) {
    return err;
    }

    /* replace the dependent var */
    list[1] = dset.v - 1;

    aux = lsq(list, dset, GretlCmdIndex.ols, Gretlopt.OPT_A);
    if (aux.errcode) {
    err = aux.errcode;
    fprintf(stderr, "auxiliary regression failed\n".ptr);
    } else {
    double pval, trsq = aux.rsq * aux.nobs;
    int df = aux.ncoeff - pmod.ncoeff;

    if (df <= 0) {
        /* couldn't add any regressors */
        err = GretlError.singular;
        goto bailout;
    }

    pval = chisq_cdf_comp(df, trsq);
    aux.aux = aux_code;

    if (!(opt & Gretlopt.OPT_I)) {
        /* OPT_I = --silent */
        if (opt & Gretlopt.OPT_Q) {
        nonlin_test_header(aux_code, prn);
        } else {
        printmodel(&aux, dset, Gretlopt.OPT_NONE, prn);
        pputc(prn, '\n');
        }
        pprintf(prn, "  %s: TR^2 = %g,\n  ".ptr, _("Test statistic").ptr, trsq);
        pprintf(prn, "%s = P(%s(%d) > %g) = %g\n\n".ptr,
            _("with p-value").ptr, _("Chi-square").ptr, df, trsq, pval);
    }

    if (opt & Gretlopt.OPT_S) {
        ModelTest *test;

        test = model_test_new((aux_code == ModelAuxCode.sq)?
                  ModelTestType.squares : ModelTestType.logs);
        if (test != null) {
        model_test_set_teststat(test, TestStats.lm);
        model_test_set_dfn(test, df);
        model_test_set_value(test, trsq);
        model_test_set_pvalue(test, chisq_cdf_comp(df, trsq));
        maybe_add_test_to_model(pmod, test);
        }
    }

    record_test_result(trsq, pval);
    }

 bailout:

    clear_model(&aux);

    return err;
}

/**
 * nonlinearity_test:
 * @pmod: pointer to original model.
 * @dset: dataset struct.
 * @aux: AUX_SQ for squares or AUX_LOG for logs
 * @opt: if contains OPT_S, save test results to model; if
 * contains OPT_I, run silently.
 * @prn: gretl printing struct.
 *
 * Run an auxiliary regression to test @pmod for nonlinearity,
 * via the addition of either squares or logs of the original
 * indepdendent variables.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int nonlinearity_test(Model * pmod, DatasetRaw * dset, ModelAuxCode aux, Gretlopt opt, Prn * prn)
{
    int save_t1 = dset.t1;
    int save_t2 = dset.t2;
    int * tmplist = null;
    const(int) orig_nvar = dset.v;
    int err = 0;

    if (!command_ok_for_model(GretlCmdIndex.add, Gretlopt.OPT_NONE, pmod)) {
    return GretlError.notimp;
    }

    if (pmod.ci == GretlCmdIndex.logistic || pmod.ci == GretlCmdIndex.lad) {
    return GretlError.notimp;
    }

    /* check for changes in original list members */
    err = list_members_replaced(pmod, dset);
    if (err) {
    return err;
    }

    /* re-impose the sample that was in force when the original model
       was estimated */
    impose_model_smpl(pmod, dset);

    /* add squares or logs */
    tmplist = augment_regression_list(pmod.list, aux, dset, &err);
    if (err) {
    return err;
    } else if (tmplist[0] == pmod.list[0]) {
    /* no vars were added */
    if (aux == ModelAuxCode.sq) {
        fprintf(stderr, "gretl: generation of squares failed\n".ptr);
        err = GretlError.squares;
    } else if (aux == ModelAuxCode.log) {
        fprintf(stderr, "gretl: generation of logs failed\n".ptr);
        err = GretlError.logs;
    }
    }

    if (!err) {
    err = real_nonlinearity_test(pmod, tmplist, dset, aux,
                     opt, prn);
    }

    /* trash any extra variables generated (squares, logs) */
    dataset_drop_last_variables(dset, dset.v - orig_nvar);

    /* put back into dset what was there on input */
    dset.t1 = save_t1;
    dset.t2 = save_t2;

    free(tmplist);

    return err;
}

int add_vars_missing(const(Model) * pmod, const(int) * addvars, const(DatasetRaw) * dset)
{
    int i, vi, t;

    for (t=pmod.t1; t<=pmod.t2; t++) {
    if (model_missing(pmod, t) ||
        (pmod.yhat != null && na(pmod.yhat[t]))) {
        continue;
    }
    for (i=1; i<=addvars[0]; i++) {
        vi = addvars[i];
        if (na(dset.Z[vi][t])) {
        fprintf(stderr, "add: obs %d OK in model but missing for series %s\n".ptr,
            t+1, dset.varname[vi]);
        return GretlError.missdata;
        }
    }
    }

    return 0;
}

Model LM_add_test(Model * pmod, DatasetRaw * dset, int * list, Gretlopt opt, Prn * prn)
{
    Model aux;
    int err;

    err = add_residual_to_dataset(pmod, dset);
    if (err) {
    gretl_model_init(&aux, dset);
    aux.errcode = err;
    return aux;
    }

    list[1] = dset.v - 1;

    aux = lsq(list, dset, GretlCmdIndex.ols, Gretlopt.OPT_A | Gretlopt.OPT_Z);

    if (aux.errcode) {
    fprintf(stderr, "auxiliary regression failed\n".ptr);
    } else {
    int df = aux.ncoeff - pmod.ncoeff;

    if (df <= 0) {
        /* couldn't add any regressors */
        aux.errcode = GretlError.singular;
    } else {
        if (!(opt & (Gretlopt.OPT_Q | Gretlopt.OPT_I))) {
        aux.aux = ModelAuxCode.add;
        printmodel(&aux, dset, Gretlopt.OPT_S, prn);
        }
    }
    }

    return aux;
}

/* Process the --auto option to the "add" or "omit" command. We should
   have either the standard abbreviation for one of the Information
   Criteria, or an alpha value, as parameter to the option.
*/

int process_auto_option(int ci, Model * orig, int * crit, double * alpha)
{
    static const(char) *[3] cstrs = ["AIC".ptr, "BIC".ptr, "HQC".ptr];
    const(char) * s = get_optval_string(cast(GretlCmdIndex) ci, Gretlopt.OPT_A);
    int i, err = 0;

    if (ci == GretlCmdIndex.omit && s == null) {
        /* omit: the --auto parameter is optional, and
           defaults to using a p-value of 0.10
        */
        *crit = InfoCriterion.max;
        *alpha = 0.10;
        return 0;
    }

    for (i=0; i<InfoCriterion.max; i++) {
        if (!strcmp(s, cstrs[i])) {
            *crit = i;
        }
    }

    if (*crit < 0) {
        /* not yet determined */
        *alpha = gretl_double_from_string(s, &err);
        if (!err && (*alpha < 0.001 || *alpha > 0.99)) {
            err = GretlError.invarg;
        } else {
            *crit = InfoCriterion.max;
        }
    }

    return err;
}

/**
 * add_test_full:
 * @orig: pointer to original model.
 * @pmod: pointer to receive augmented model.
 * @addvars: list of variables to add to original model.
 * @dset: dataset struct.
 * @opt: can contain OPT_Q (quiet) to suppress printing
 * of the new model, OPT_O to print its covariance matrix,
 * OPT_I for silent operation.
 * @prn: gretl printing struct.
 *
 * Re-estimate a given model after adding the specified
 * variables, and records a joint test on the additional
 * variables.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int add_test_full(Model * orig, Model * pmod, const(int) * addvars, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    Model umod;
    int save_t1 = dset.t1;
    int save_t2 = dset.t2;
    int * biglist = null;
    int * add_order = null;
    const(int) orig_nvar = dset.v;
    double alpha = 0;
    int n_add = 0;
    int crit = -1;
    int err = 0;

    if (orig == null || orig.list == null || addvars == null) {
    return GretlError.data;
    }

    n_add = addvars[0];
    if (n_add == 0) {
    return GretlError.noadd;
    }

    if (incompatible_options(opt, Gretlopt.OPT_A | Gretlopt.OPT_L)) {
        /* --auto and --lm cannot be combined */
        return GretlError.badopt;
    }

    if (!command_ok_for_model(GretlCmdIndex.add, opt, orig)) {
    return GretlError.notimp;
    }

    if ((opt & Gretlopt.OPT_L) && pmod != null) {
    /* --lm option incompatible with "full" variant of add test */
    return GretlError.badopt;
    }

    if (exact_fit_check(orig, prn)) {
    return 0;
    }

    /* check for changes in original list members */
    err = list_members_replaced(orig, dset);
    if (err) {
    return err;
    }

    /* check for NAs in add list relative to model */
    err = add_vars_missing(orig, addvars, dset);
    if (err) {
    return err;
    }

    if (opt & Gretlopt.OPT_A) {
        /* check the --auto option */
        err = process_auto_option(GretlCmdIndex.add, orig, &crit, &alpha);
    } else {
        /* create augmented regression list */
        if (orig.ci == GretlCmdIndex.ivreg) {
            biglist = ivreg_list_add(orig.list, addvars, opt, &err);
        } else if (orig.ci == GretlCmdIndex.dpanel) {
            biglist = panel_list_add(orig, addvars, &err);
        } else {
            biglist = gretl_list_add(orig.list, addvars, &err);
        }
    }

    if (err) {
    return err;
    }

    /* impose the sample range in force when the
       original model was estimated
    */
    impose_model_smpl(orig, dset);

    if (opt & Gretlopt.OPT_A) {
        /* Do stepwise augmentation of the original model */
        alias StepwiseAddFn = extern (C) Model function(Model *, const(int) *, int **, int, double, DatasetRaw *, Gretlopt, Prn *);
        StepwiseAddFn stepwise_add;

        stepwise_add = cast(StepwiseAddFn) get_plugin_function("stepwise_add");
        if (stepwise_add == null) {
            err = GretlError.fopen;
        } else {
            umod = stepwise_add(orig, addvars, &add_order,
                crit, alpha, dset, opt, prn);
        }
    } else if (opt & Gretlopt.OPT_L) {
    /* run an LM test */
    umod = LM_add_test(orig, dset, biglist, opt, prn);
    } else {
    /* Run augmented regression, matching the original estimation
       method; use OPT_Z to suppress the elimination of perfectly
       collinear variables.
    */
    Gretlopt ropt = Gretlopt.OPT_Z;

    if (opt & (Gretlopt.OPT_Q | Gretlopt.OPT_I)) {
        /* not printing @umod, so pass along OPT_Q */
        ropt |= Gretlopt.OPT_Q;
    }
    umod = replicate_estimator(orig, biglist, dset, ropt, prn);
    }

    if (umod.errcode) {
    err = umod.errcode;
    errmsg(err, prn);
    }

    if (!(opt & Gretlopt.OPT_A)) {
        if (!err && umod.ncoeff - orig.ncoeff != n_add) {
            gretl_errmsg_sprintf(_("Failed to add %d variable(s)").ptr, n_add);
            err = GretlError.data;
        }
        if (!err) {
            err = add_or_omit_compare(orig, &umod, addvars,
                                      dset, GretlCmdIndex.add, opt, prn);
        }
    }

    if (err || pmod == null) {
    clear_model(&umod);
    free(add_order);
    } else {
    *pmod = umod;
    if (add_order != null) {
        gretl_model_set_list_as_data(pmod, "add_order", add_order);
    }
    }

    /* put dset back as it was on input */
    dataset_drop_last_variables(dset, dset.v - orig_nvar);
    dset.t1 = save_t1;
    dset.t2 = save_t2;

    free(biglist);

    return err;
}

/**
 * add_test:
 * @pmod: pointer to model to be tested.
 * @addvars: list of variables to test.
 * @dset: dataset struct.
 * @opt: can contain OPT_Q (quiet) to suppress printing
 * of the auxiliary model, OPT_I to suppress all printing
 * of results.
 * @prn: gretl printing struct.
 *
 * Performs an LM test on @pmod for the null hypothesis
 * that the @addvars variables do not contribute
 * significant explanatory power.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int add_test(Model * pmod, const(int) * addvars, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    return add_test_full(pmod, null, addvars, dset, opt, prn);
}

int do_wald_omit(const(int) * list, Model * pmod, const(DatasetRaw) * dset, Gretlopt opt, Prn * prn)
{
    int * test = null;
    int err = 0;

    /* test validity of omissions */

    if (pmod.ci == GretlCmdIndex.ivreg) {
    test = ivreg_list_omit(pmod.list, list, opt, &err);
    } else if (pmod.ci == GretlCmdIndex.panel || pmod.ci == GretlCmdIndex.dpanel) {
    test = panel_list_omit(pmod, list, &err);
    } else {
    test = gretl_list_omit(pmod.list, list, 2, &err);
    }

    if (!err) {
    free(test);
    err = add_or_omit_compare(pmod, null, list, dset, GretlCmdIndex.omit,
                  opt, prn);
    }

    return err;
}

struct OmitInfo {
    Model * orig; /* the original model */
    Model * curr; /* the current model in stepwise elimination */
    Model * tmp; /* temporary workspace model */
    int * list; /* the list, from which regressors may be dropped */
    const(int) * cands;  /* list of candidates for elimination, in case some
                          regressors are to be kept regardless */
    double alpha; /* the max alpha when using the p-value criterion */
    double cval; /* value of the criterion in use */
    byte starting; /* flag indicating the first step */
    byte crit; /* identifier of the criterion in use */
    byte use_pval; /* flag indicating that p-values are in use */
}

/* Check whether coefficient @i corresponds to a variable that is
   removable from the model: this is the case if either (a) the @cands
   list is empty, or (b) coefficient @i is the coefficient on one of the
   variables in @cands.
*/

int coeff_is_removable(const(int) * cands, const(Model) * pmod, DatasetRaw * dset, int i)
{
    int ret = 1;

    if (cands != null && cands[0] > 0) {
    const(char) * vname;
    int j, pj;

        /* We're only "allowed" to drop members of the @cands list,
           so reverse the presumption.
        */
    ret = 0;
    for (j=1; j<=cands[0]; j++) {
        vname = dset.varname[cands[j]];
        pj = gretl_model_get_param_number(pmod, dset, vname);
        if (pj == i) {
        ret = 1;
        break;
        }
    }
    }

    return ret;
}

void print_drop(Model * pmod, OmitInfo * oi, int k, DatasetRaw * dset, Prn * prn)
{
    const(char) *[] cstrs = ["AIC", "BIC", "HQC", "p-value"];

    if (oi.starting) {
        pputc(prn, '\n');
        pprintf(prn, _("Sequential elimination using %s").ptr, cstrs[oi.crit]);
        if (oi.use_pval) {
            pprintf(prn, " (α = %.2f)".ptr, oi.alpha);
        }
        pputs(prn, "\n\n");
    }

    pprintf(prn, _(" Dropping %-16s (%s %#g)\n").ptr,
            dset.varname[oi.list[k+2]],
            cstrs[oi.crit], oi.cval);
}

/* Determine whether @pmod contains a regressor that should be dropped
   according to a specified criterion (either the p-value of its
   coefficient exceeds a specified threshold, or its omission results
   in an improvement in a specified information criterion), If so, remove
   this regressor from oi->list and return 1; otherwise return 0.

   Note that if the list oi->cands is non-empty then the search is
   confined to variables in that list.
*/

int auto_drop_var(OmitInfo * oi, DatasetRaw * dset, int ** orderp, Gretlopt opt, Prn * prn, int * err)
{
    Model * pmod;
    Model * tmp = oi.tmp;
    double tstat, tmin = 4.0;
    int imin, imax;
    int i, k = -1;
    int do_drop = 0;
    int ret = 0;

    pmod = oi.starting ? oi.orig : oi.curr;

    imin = pmod.ifc;
    if (pmod.xlist != null) {
    imax = pmod.xlist[0];
    } else {
    imax = pmod.ncoeff;
    }

    for (i=imin; i<imax; i++) {
    if (coeff_is_removable(oi.cands, pmod, dset, i)) {
        tstat = fabs(pmod.coeff[i] / pmod.sderr[i]);
        if (tstat < tmin) {
        tmin = tstat;
        k = i;
        }
    }
    }

    if (k >= 0) {
        if (oi.use_pval) {
            /* Using the p-value approach: we just have to see if
               regressor with the smallest t-ratio has a p-value
               above the max for retention.
            */
            oi.cval = coeff_pval(pmod.ci, tmin, pmod.dfd);
            do_drop = oi.cval > oi.alpha;
        } else {
            /* Using an info criterion: we have to re-estimate the
               model without the regressor with the smallest t-ratio
               to see if its omission improves the criterion.
            */
            int * ltmp = gretl_list_copy(oi.list);
            double crit0 = pmod.criterion[oi.crit];

            gretl_list_delete_at_pos(ltmp, k + 2);
            if (!oi.starting) {
                clear_model(oi.curr);
            }
            set_reference_missmask_from_model(oi.orig);
            *tmp = replicate_estimator(oi.orig, ltmp, dset, Gretlopt.OPT_A, prn);
            if (tmp.errcode) {
                *err = tmp.errcode;
                clear_model(tmp);
            } else {
                oi.cval = tmp.criterion[oi.crit];
                do_drop = oi.cval < crit0;
            }
            free(ltmp);
        }
    }

    if (do_drop) {
    gretl_list_append_term(orderp, oi.list[k+2]);
        if (!(opt & Gretlopt.OPT_I)) {
            /* not --silent */
            print_drop(pmod, oi, k, dset, prn);
        }
        if (oi.use_pval) {
            /* this is already done in the info criterion case */
            gretl_list_delete_at_pos(oi.list, k + 2);
        } else {
            if (!oi.starting) {
                clear_model(oi.curr);
            }
            /* attach the reduced model that was estimated above */
            *(oi.curr) = *tmp;
        }
        ret = 1;
    }

    oi.starting = 0;

    return ret;
}

void list_copy_values(int * targ, const(int) * src)
{
    int i;

    for (i=0; i<=src[0]; i++) {
    targ[i] = src[i];
    }
}

/* Run a loop in which the least significant variable is dropped from
   the regression list, provided its p-value exceeds some specified
   cutoff (or dropping it improves a specified information criterion).
   If @omitlist is non-empty the routine is confined to members of that
   list.

   FIXME this probably still needs work for estimators other than OLS.
*/

Model auto_omit(Model * orig, const(int) * omitlist, int ** orderp, int crit, double alpha, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    OmitInfo oi;
    memset(&oi, 0, oi.sizeof);
    Model omod;
    int allgone = 0;
    int i, drop;
    int err = 0;

    gretl_model_init(&omod, dset);

    oi.list = gretl_list_copy(orig.list);
    if (oi.list == null) {
    omod.errcode = GretlError.alloc;
    return omod;
    }

    oi.orig = orig;
    oi.curr = &omod;
    oi.cands = omitlist;
    oi.alpha = alpha;
    oi.starting = 1;
    oi.crit = cast(byte) crit;
    oi.use_pval = (crit == InfoCriterion.max);
    oi.tmp = oi.use_pval ? null : gretl_model_new();

    drop = auto_drop_var(&oi, dset, orderp, opt, prn, &err);
    if (err) {
        omod.errcode = err;
    } else if (!drop) {
    /* nothing was dropped: set a "benign" error code */
    err = omod.errcode = GretlError.noomit;
    }

    for (i=0; drop > 0 && !allgone && !err; i++) {
    if (i > 0) {
        set_reference_missmask_from_model(orig);
    }
        if (oi.use_pval) {
            /* in the info criterion case replicate_estimator() is called
               by auto_drop_var()
            */
            omod = replicate_estimator(orig, oi.list, dset, Gretlopt.OPT_A, prn);
            err = omod.errcode;
            if (err) {
                break;
            }
        }
        list_copy_values(oi.list, omod.list);
        drop = auto_drop_var(&oi, dset, orderp, opt, prn, &err);
        if (drop && omod.ncoeff == 1) {
            allgone = 1; /* will break */
        }
        if (oi.use_pval) {
            clear_model(&omod);
        }
    }

    if (err) {
        omod.errcode = err;
    } else if (allgone) {
    pputc(prn, '\n');
    pprintf(prn, _("No coefficient has a p-value less than %g").ptr, alpha);
    pputc(prn, '\n');
    } else {
    /* re-estimate the final_ model without the auxiliary flag */
    Gretlopt ropt = Gretlopt.OPT_NONE;

    if (opt & (Gretlopt.OPT_Q | Gretlopt.OPT_I)) {
        /* not printing @omod, so pass along OPT_Q */
        ropt |= Gretlopt.OPT_Q;
    }
    set_reference_missmask_from_model(orig);
        omod = replicate_estimator(orig, oi.list, dset, ropt, prn);
    }

    if (oi.tmp != null) {
        gretl_model_free(oi.tmp);
    }

    free(oi.list);

    return omod;
}

/* create reduced list for "omit" test on model, based on
   the list of variables to be dropped, @omitvars
*/

int make_short_list(Model * orig, const(int) * omitvars, Gretlopt opt, int ** plist)
{
    int * list = null;
    int err = 0;

    if (omitvars == null || omitvars[0] == 0) {
    return GretlError.parse;
    }

    if (orig.ci == GretlCmdIndex.ivreg) {
    list = ivreg_list_omit(orig.list, omitvars, opt, &err);
    } else if (orig.ci == GretlCmdIndex.panel || orig.ci == GretlCmdIndex.dpanel) {
    list = panel_list_omit(orig, omitvars, &err);
    } else {
    list = gretl_list_omit(orig.list, omitvars, 2, &err);
    }

    if (list != null && list[0] == 1) {
    /* only the dependent variable would be left */
    err = GretlError.novars;
    }

    *plist = list;

    return err;
}

int omit_options_inconsistent(Model * pmod, Gretlopt opt)
{
    if (opt & Gretlopt.OPT_B) {
    /* 2sls: omitting variable as instrument */
    if (opt & (Gretlopt.OPT_W | Gretlopt.OPT_A)) {
        /* can't use Wald method on original VCV,
           and can't do --auto */
        return 1;
    } else if (pmod.ci != GretlCmdIndex.ivreg) {
        return 1;
    }
    }

    if ((opt & Gretlopt.OPT_A) && (opt & Gretlopt.OPT_W)) {
    /* --auto and --wald options incompatible */
    return 1;
    }

    return 0;
}

int const_only(Model * pmod)
{
    return (pmod.list[0] == 2 && pmod.list[2] == 0);
}

int omit_test_precheck(Model * pmod, const(int) * omitvars, Gretlopt opt)
{
    if (pmod == null || pmod.list == null) {
    return GretlError.data;
    } else if (!command_ok_for_model(GretlCmdIndex.omit, Gretlopt.OPT_NONE, pmod)) {
    return GretlError.notimp;
    } else if (omit_options_inconsistent(pmod, opt)) {
    return GretlError.badopt;
    } else if (const_only(pmod)) {
    gretl_errmsg_set(_("There are no variables to omit").ptr);
    return GretlError.data;
    }

    return 0;
}

/* This is about whether we should use the relatively new code in the
   "stepwise" plugin rather than the code here in compare.c.
*/

int use_stepwise_omit(Model * orig, int crit)
{
    if (orig.ci != GretlCmdIndex.ols || (orig.opt & Gretlopt.OPT_R)) {
        /* we only handle OLS with classical std errors */
        return 0;
    } else if (crit > InfoCriterion.hqc) {
        /* we only handle the info criterion approach */
        return 0;
    } else if (getenv("USE_AUTO_OMIT")) {
        /* we allow forcing use of the older code */
        return 0;
    } else {
        /* OK, use stepwise */
        return 1;
    }
}

/**
 * omit_test_full:
 * @orig: pointer to original model.
 * @pmod: pointer to receive new model, with vars omitted.
 * @omitvars: list of variables to omit from original model.
 * @dset: dataset struct.
 * @opt: can contain OPT_Q (quiet) to suppress printing
 * of the new model, OPT_O to print its covariance matrix,
 * OPT_I for silent operation; for OPT_A, see below.
 * @prn: gretl printing struct.
 *
 * Re-estimate a given model after removing the variables
 * specified in @omitvars.  Or if OPT_A is given, proceed
 * sequentially, at each step dropping the least significant
 * variable provided its p-value is above a certain threshold
 * (by default 0.10, two-sided).
 *
 * Returns: 0 on successful completion, error code on error.
 */

int omit_test_full(Model * orig, Model * pmod, const(int) * omitvars, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    Model rmod;
    int save_t1 = dset.t1;
    int save_t2 = dset.t2;
    int * tmplist = null;
    int * omit_order = null;
    double alpha = 0;
    int crit = -1;
    int err;

    err = omit_test_precheck(orig, omitvars, opt);
    if (err) {
    return err;
    }

    /* check that vars to omit have not been redefined */
    if ((err = list_members_replaced(orig, dset)) != 0) {
    return err;
    }

    if (opt & Gretlopt.OPT_A) {
        /* doing auto-omit */
        err = process_auto_option(GretlCmdIndex.omit, orig, &crit, &alpha);
    } else {
    err = make_short_list(orig, omitvars, opt, &tmplist);
    }

    if (err) {
        free(tmplist);
        return err;
    }

    /* impose the sample range used for the original model */
    impose_model_smpl(orig, dset);

    /* set the mask for missing obs within the sample range, based
       on the original model */
    set_reference_missmask_from_model(orig);

    if ((opt & Gretlopt.OPT_A) && use_stepwise_omit(orig, crit)) {
        alias StepwiseOmitFn = extern (C) Model function(Model *, const(int) *, int **, int, double, DatasetRaw *, Gretlopt, Prn *);
        StepwiseOmitFn stepwise_omit;

        stepwise_omit = cast(StepwiseOmitFn) get_plugin_function("stepwise_omit");
        if (stepwise_omit == null) {
            err = GretlError.fopen;
        } else {
            rmod = stepwise_omit(orig, omitvars, &omit_order,
                 crit, alpha, dset, opt, prn);
        }
    } else if (opt & Gretlopt.OPT_A) {
    rmod = auto_omit(orig, omitvars, &omit_order,
             crit, alpha, dset, opt, prn);
    } else {
    Gretlopt ropt = Gretlopt.OPT_NONE;

    if (opt & (Gretlopt.OPT_Q | Gretlopt.OPT_I)) {
        /* not printing @rmod, so pass along OPT_Q */
        ropt |= Gretlopt.OPT_Q;
    }
    rmod = replicate_estimator(orig, tmplist, dset, ropt, prn);
    }

    err = rmod.errcode;

    if (err) {
    if (err == GretlError.noomit && (opt & Gretlopt.OPT_I)) {
        ; /* --silent: keep quiet */
    } else {
        errmsg(err, prn);
    }
    } else {
    Model * newmod = null;
    int minpos = 2;
    int * omitlist = null;

    if (orig.ci == GretlCmdIndex.dpanel) {
        /* skip AR spec, separator, and dep var */
        minpos = 4;
    }

    if (rmod.list == null) {
        /* handle the case where sequential elimination
           ended up dropping all regressors */
        int * minlist = gretl_list_new(1);

        minlist[1] = orig.list[1];
        omitlist = gretl_list_diff_new(orig.list, minlist, minpos);
        free(minlist);
    } else {
        newmod = &rmod;
        omitlist = gretl_list_diff_new(orig.list, rmod.list, minpos);
    }
    if (omitlist != null) {
        err = add_or_omit_compare(orig, newmod, omitlist,
                      dset, GretlCmdIndex.omit, opt, prn);
        free(omitlist);
    }
    }

    if (err || pmod == null) {
    clear_model(&rmod);
    free(omit_order);
    } else {
    *pmod = rmod;
    if (omit_order != null) {
        gretl_model_set_list_as_data(pmod, "omit_order", omit_order);
    }
    }

    /* put back into dset what was there on input */
    dset.t1 = save_t1;
    dset.t2 = save_t2;

    free(tmplist);

    return err;
}

/**
 * omit_test:
 * @pmod: pointer to model to be tested.
 * @omitvars: list of variables to test.
 * @dset: dataset struct.
 * @opt: can contain OPT_Q (quiet) to suppress printing
 * of results.
 * @prn: gretl printing struct.
 *
 * Performs a Wald test on @pmod for the null hypothesis
 * that the @omitvars variables do not contribute
 * explanatory power.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int omit_test(Model * pmod, const(int) * omitvars, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    int err = omit_test_precheck(pmod, omitvars, opt);

    if (!err) {
    err = do_wald_omit(omitvars, pmod, dset, opt, prn);
    }

    return err;
}

/**
 * get_DW_pvalue_for_model:
 * @pmod: model to be tested.
 * @dset: dataset struct.
 * @err: location to receive error code.
 *
 * Computes the p-value for the Durbin-Watson statistic for the
 * given model, using the Imhof method.
 *
 * Returns: the p-value, or #NADBL on error.
 */

double get_DW_pvalue_for_model(Model * pmod, DatasetRaw * dset, int * err)
{
    Model dwmod;
    int save_t1 = dset.t1;
    int save_t2 = dset.t2;
    int * list = null;
    double pv = NADBL;

    /* maybe this has already been done? */
    pv = gretl_model_get_double(pmod, "dw_pval");
    if (!na(pv)) {
    return pv;
    }

    if (pmod.ci == GretlCmdIndex.panel && panel_DW_pval_ok(pmod)) {
    /* Use Bhargava et al approximation */
    return BFN_panel_DW_pvalue(pmod, dset, err);
    }

    if (dset == null || dset.Z == null) {
    *err = GretlError.nodata;
    } else if (pmod == null || pmod.list == null) {
    *err = GretlError.data;
    } else if ((pmod.ci != GretlCmdIndex.ols && pmod.ci != GretlCmdIndex.panel) ||
           na(pmod.dw) || model_has_missing_obs(pmod)) {
    *err = GretlError.badstat;
    } else {
    /* check that relevant vars have not been redefined */
    *err = list_members_replaced(pmod, dset);
    }

    if (!*err) {
    list = gretl_list_copy(pmod.list);
    if (list == null) {
        *err = GretlError.alloc;
    }
    }

    if (*err) {
    return NADBL;
    }

    gretl_model_init(&dwmod, dset);

    /* impose the sample range used for the original model */
    impose_model_smpl(pmod, dset);

    dwmod = replicate_estimator(pmod, list, dset, Gretlopt.OPT_A | Gretlopt.OPT_I, null);
    *err = dwmod.errcode;

    if (!*err) {
    pv = gretl_model_get_double(&dwmod, "dw_pval");
    if (!na(pv)) {
        /* record result on the incoming model */
        gretl_model_set_double(pmod, "dw_pval", pv);
    }
    }

    /* put back into dset what was there on input */
    dset.t1 = save_t1;
    dset.t2 = save_t2;

    clear_model(&dwmod);
    free(list);

    return pv;
}

/**
 * reset_test:
 * @pmod: pointer to model to be tested.
 * @dset: dataset struct.
 * @opt: if contains %OPT_S, save test results to model. %OPT_Q
 * suppresses the printout of the auxiliary regression. %OPT_U and
 * %OPT_C stand for "squares only" and "cubes only", respectively.
 * %OPT_I produces silent operation.
 * @prn: gretl printing struct.
 *
 * Carries out and prints Ramsey's RESET test for model specification.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int reset_test(Model * pmod, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    int * newlist = null;
    Model aux;
    double RF;
    int save_t1 = dset.t1;
    int save_t2 = dset.t2;
    int i, t, orig_v = dset.v;
    int addv, use_square, use_cube;
    int robust = 0;
    const(char) * mode;
    int err = 0;

    if (pmod.ci != GretlCmdIndex.ols) {
    return GretlError.olsonly;
    }

    if ((opt & Gretlopt.OPT_R) || (pmod.opt & Gretlopt.OPT_R)) {
    robust = 1;
    }

    err = incompatible_options(opt, Gretlopt.OPT_C | Gretlopt.OPT_U);

    if (err) {
    return err;
    }

    if (exact_fit_check(pmod, prn)) {
    return 0;
    }

    use_square = !(opt & Gretlopt.OPT_C); /* not cubes-only */
    use_cube = !(opt & Gretlopt.OPT_U);   /* not squares-only */

    gretl_model_init(&aux, dset);

    if (opt & Gretlopt.OPT_U) {
    addv = 1;
    mode = "squares only".ptr;
    } else if (opt & Gretlopt.OPT_C) {
    addv = 1;
    mode = "cubes only".ptr;
    } else {
    addv = 2;
    mode = "squares and cubes".ptr;
    }

    impose_model_smpl(pmod, dset);

    if (pmod.ncoeff + addv >= dset.t2 - dset.t1) {
    /* can't run aux regression */
    err = GretlError.df;
    }

    if (!err) {
    newlist = gretl_list_new(pmod.list[0] + addv);
    if (newlist == null) {
        err = GretlError.alloc;
    }
    }

    if (!err) {
    for (i=1; i<=pmod.list[0]; i++) {
        newlist[i] = pmod.list[i];
    }
    if (dataset_add_series(dset, addv)) {
        err = GretlError.alloc;
    }
    }

    if (!err) {
    /* add scaled yhat^2 and/or yhat^3 to data set */
    int vs = orig_v;
    int vc = (opt & Gretlopt.OPT_C)? orig_v : orig_v + 1;
    int k = pmod.list[0] + 1;
        double s;

        s = gretl_stddev(pmod.t1, pmod.t2, pmod.yhat);

    for (t = pmod.t1; t<=pmod.t2; t++) {
        double x = pmod.yhat[t] / s;

        if (use_square) {
        dset.Z[vs][t] = x * x;
        }
        if (use_cube) {
        dset.Z[vc][t] = x * x * x;
        }
    }

    if (use_square) {
        strcpy(dset.varname[vs], "yhat^2");
        newlist[k++] = vs;
    }
    if (use_cube) {
        strcpy(dset.varname[vc], "yhat^3");
        newlist[k++] = vc;
    }
    }

    if (!err) {
    Gretlopt auxopt = Gretlopt.OPT_A;

    if (robust) {
        auxopt |= Gretlopt.OPT_R;
    }
    aux = lsq(newlist, dset, GretlCmdIndex.ols, auxopt);
    err = aux.errcode;
    if (err) {
        errmsg(aux.errcode, prn);
    }
    }

    if (!err) {
    int silent = (opt & Gretlopt.OPT_I);
    double pval;

    aux.aux = ModelAuxCode.reset;

    if (!silent) {
        if (!(opt & Gretlopt.OPT_Q)) {
        printmodel(&aux, dset, Gretlopt.OPT_NONE, prn);
        } else {
        if (!(opt & Gretlopt.OPT_G)) {
            /* GUI special; see gui2/library.c */
            pputc(prn, '\n');
        }
        if (robust) {
            pputs(prn, _("Robust RESET test for specification").ptr);
        } else {
            pputs(prn, _("RESET test for specification").ptr);
        }
        pprintf(prn, " (%s)\n".ptr, _(mode));
        }
    }

    if (robust) {
        /* FIXME wald_omit_F(omitlist, pmod); */
        int * omitlist = gretl_list_diff_new(aux.list, pmod.list, 2);

        RF = wald_omit_F(omitlist, &aux);
        free(omitlist);
    } else {
        RF = ((pmod.ess - aux.ess) / addv) / (aux.ess / aux.dfd);
    }
    pval = snedecor_cdf_comp(addv, aux.dfd, RF);

    if (!silent) {
        const(char) * h0 = get_h0_string_for_test(ModelTestType.reset);

        pprintf(prn, "%s: %s\n".ptr, _("Null hypothesis").ptr, _(h0));
        pprintf(prn, "%s: F = %f,\n".ptr, _("Test statistic").ptr, RF);
        pprintf(prn, "%s = P(F(%d,%d) > %g) = %.3g\n".ptr, _("with p-value").ptr,
            addv, aux.dfd, RF, pval);
        pputc(prn, '\n');
    }

    if (opt & Gretlopt.OPT_S) {
        ModelTest *test = model_test_new(ModelTestType.reset);
        Gretlopt topt = Gretlopt.OPT_NONE;

        if (test != null) {
        if (robust) {
            topt = Gretlopt.OPT_R;
        }
        if (opt & Gretlopt.OPT_U) {
            topt = Gretlopt.OPT_U;
        } else if (opt & Gretlopt.OPT_C) {
            topt = Gretlopt.OPT_C;
        }
        model_test_set_teststat(test, TestStats.reset);
        model_test_set_dfn(test, addv);
        model_test_set_dfd(test, aux.dfd);
        model_test_set_value(test, RF);
        model_test_set_pvalue(test, pval);
        model_test_set_opt(test, topt);
        maybe_add_test_to_model(pmod, test);
        }
    }

    record_test_result(RF, pval);
    }

    free(newlist);
    dataset_drop_last_variables(dset, addv);
    clear_model(&aux);

    dset.t1 = save_t1;
    dset.t2 = save_t2;

    return err;
}

void bg_test_header(int order, Prn * prn, int ivreg)
{
    if (ivreg) {
    pprintf(prn, "\n%s ".ptr, _("Godfrey (1994) test for").ptr);
    } else {
    pprintf(prn, "\n%s ".ptr, _("Breusch-Godfrey test for").ptr);
    }

    if (order > 1) {
    pprintf(prn, "%s %d\n".ptr, _("autocorrelation up to order").ptr,
        order);
    } else {
    pprintf(prn, "%s\n".ptr, _("first-order autocorrelation").ptr);
    }

    pputc(prn, '\n');
}

double ivreg_autocorr_wald_stat(Model * aux, int order, int * err)
{
    GretlVector * b = gretl_vector_alloc(order);
    GretlMatrixRaw * V1 = gretl_matrix_alloc(order, order);
    GretlVector * WT = gretl_vector_alloc(1);
    GretlMatrixRaw * V0 = null;
    double x = NADBL;
    int i, j, ki, kj;

    if (b == null || V1 == null || WT == null) {
    *err = GretlError.alloc;
    } else {
    V0 = gretl_model_get_matrix(aux, ModelDataIndex.vcv, err);
    }

    if (!*err) {
    ki = aux.ncoeff - order;

    for (i=0; i<order; i++) {
        x = aux.coeff[ki];
        gretl_vector_set(b, i, x);
        kj = ki;
        for (j=i; j<order; j++) {
        x = gretl_matrix_get(V0, ki, kj);
        gretl_matrix_set(V1, i, j, x);
        gretl_matrix_set(V1, j, i, x);
        kj++;
        }
        ki++;
    }
    *err = gretl_invert_symmetric_matrix(V1);
    }

    if (!*err) {
    gretl_matrix_qform(b, GretlMatrixMod.none, V1, WT, GretlMatrixMod.none);
    x = gretl_vector_get(WT, 0) / order;
    }

    gretl_vector_free(WT);
    gretl_vector_free(b);
    gretl_matrix_free(V1);
    gretl_matrix_free(V0);

    return x;
}

/**
 * ivreg_autocorr_test:
 * @pmod: pointer to model to be tested.
 * @order: lag order for test.
 * @dset: dataset struct.
 * @opt: if flags include OPT_S, save test results to model;
 * if OPT_Q, be less verbose.
 * @prn: gretl printing struct.
 *
 * Tests the given IV model for autocorrelation of order equal
 * to the specified value, or equal to the frequency of the data if
 * the supplied @order is zero, as per Godfrey (1994), "Testing for
 * Serial Correlation by Variable Addition in Dynamic Models Estimated
 * by Instrumental Variables", RES. Note that none of the
 * asymptotically equivalent tests given on page 553 is used
 * here. Instead, we estimate the model augmented with lags and then
 * perform a Wald-type test. The resulting chi-square statistic is
 * divided by its degrees of freedom as a finite-sample adjustment and
 * compared to an F distribution.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int ivreg_autocorr_test(Model * pmod, int order, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    int smpl_t1 = dset.t1;
    int smpl_t2 = dset.t2;
    int v = dset.v;
    int * addlist = null;
    int * testlist = null;
    double x = 0, pval = 1;
    Model aux;
    int i, t;
    int err = 0;

    if (dataset_is_panel(dset)) {
    return GretlError.notimp;
    }

    if (model_has_missing_obs(pmod)) {
    return GretlError.missdata;
    }

    /* impose original sample range */
    impose_model_smpl(pmod, dset);

    gretl_model_init(&aux, dset);

    if (order <= 0) {
    order = dset.pd;
    }

    if (pmod.ncoeff + order >= dset.t2 - dset.t1) {
    return GretlError.df;
    }

    addlist = gretl_list_new(order);

    if (addlist == null) {
    err = GretlError.alloc;
    } else {
    err = add_residual_to_dataset(pmod, dset);
    }

    if (!err) {
    /* add lags of residual */
    for (i=1; i<=order && !err; i++) {
        int lnum;

        lnum = laggenr(v, i, dset);

        if (lnum < 0) {
        gretl_errmsg_set(_("lagging uhat failed").ptr);
        err = GretlError.lags;
        } else {
        /* set the first entries to 0 for compatibility with Godfrey (1994)
           and PcGive (not perfect) */
        for (t=smpl_t1; t<smpl_t1+i; t++) {
            dset.Z[lnum][t] = 0;
        }
        addlist[i] = lnum;
        }
    }
    if (!err) {
        /* compose augmented regression list */
        testlist = ivreg_list_add(pmod.list, addlist, Gretlopt.OPT_B, &err);
    }
    }

    if (!err) {
    Gretlopt ivopt = Gretlopt.OPT_A;

    transcribe_option_flags(&ivopt, pmod.opt,
                Gretlopt.OPT_L | Gretlopt.OPT_G | Gretlopt.OPT_R);
    aux = ivreg(testlist, dset, ivopt);
    err = aux.errcode;
    }

    if (!err) {
    x = ivreg_autocorr_wald_stat(&aux, order, &err);
    }

    if (!err) {
    aux.aux = ModelAuxCode.ar;
    gretl_model_set_int(&aux, "BG_order", order);
    pval = snedecor_cdf_comp(order, aux.nobs - pmod.ncoeff - order, x);

    if (opt & Gretlopt.OPT_Q) {
        bg_test_header(order, prn, 1);
    } else {
        printmodel(&aux, dset, Gretlopt.OPT_S, prn);
    }

    pputc(prn, '\n');
    pprintf(prn, "%s: Pseudo-LMF = %f,\n".ptr, _("Test statistic").ptr, x);
    pprintf(prn, "%s = P(F(%d,%d) > %g) = %.3g\n".ptr, _("with p-value").ptr,
        order, aux.nobs - pmod.ncoeff, x, pval);
    pputc(prn, '\n');
    record_test_result(x / order, pval);

    if (opt & Gretlopt.OPT_S) {
        ModelTest *test = model_test_new(ModelTestType.autocorr);

        if (test != null) {
        model_test_set_teststat(test, TestStats.lmf);
        model_test_set_dfn(test, order);
        model_test_set_dfd(test, aux.nobs - pmod.ncoeff);
        model_test_set_order(test, order);
        model_test_set_value(test, x);
        model_test_set_pvalue(test, pval);
        maybe_add_test_to_model(pmod, test);
        }
    }
    }

    free(addlist);
    free(testlist);

    dataset_drop_last_variables(dset, dset.v - v);
    clear_model(&aux);

    /* reset sample as it was */
    dset.t1 = smpl_t1;
    dset.t2 = smpl_t2;

    return err;
}

int lb_autocorr_test(Model * pmod, int order, Gretlopt opt, Prn * prn)
{
    double lb, pval = NADBL;
    int df, k;
    int err = 0;

    k = arma_model_get_n_arma_coeffs(pmod);
    df = order - k;

    if (df <= 0) {
    gretl_errmsg_set(_("Insufficient degrees of freedom for test").ptr);
    return GretlError.data;
    }

    lb = ljung_box(order, pmod.t1, pmod.t2, pmod.uhat, &err);

    if (!err) {
    pval = chisq_cdf_comp(df, lb);
    if (na(pval)) {
        err = GretlError.data;
    }
    }

    if (err) {
    gretl_errmsg_set(_("Error calculating Ljung-Box statistic").ptr);
    } else {
    pputc(prn, '\n');
    pprintf(prn, _("Test for autocorrelation up to order %d").ptr,
        order);
    pputs(prn, "\n\n");
    pprintf(prn, "Ljung-Box Q' = %g,\n".ptr, lb);
    pprintf(prn, "%s = P(%s(%d) > %g) = %#.4g\n".ptr, _("with p-value").ptr,
        _("Chi-square").ptr, df, lb, chisq_cdf_comp(df, lb));
    pputc(prn, '\n');
    record_test_result(lb, pval);
    }

    if (!err && (opt & Gretlopt.OPT_S)) {
    ModelTest *test = model_test_new(ModelTestType.autocorr);

    if (test != null) {
        model_test_set_teststat(test, TestStats.lb_chisq);
        model_test_set_dfn(test, df);
        model_test_set_order(test, order);
        model_test_set_value(test, lb);
        model_test_set_pvalue(test, pval);
        maybe_add_test_to_model(pmod, test);
    }
    }

    return err;
}

/**
 * regular_autocorr_test:
 * @pmod: pointer to model to be tested.
 * @order: lag order for test.
 * @dset: dataset struct.
 * @opt: if flags include OPT_S, save test results to model;
 * if OPT_Q, be less verbose; if OPT_I, be silent.
 * @prn: gretl printing struct.
 *
 * Tests the given model for autocorrelation of order equal to
 * the specified value, or equal to the frequency of the data if
 * the supplied @order is zero. Prints TR^2 and LMF test statistics.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int regular_autocorr_test(Model * pmod, int order, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    int save_t1 = dset.t1;
    int save_t2 = dset.t2;
    int * newlist = null;
    Model aux;
    double RSSxe, RSSx = pmod.ess;
    int i, t, n = dset.n, v = dset.v;
    double trsq, LMF, lb, pval = 1.0;
    int err = 0;

    gretl_model_init(&aux, dset);

    if (order <= 0) {
    order = dset.pd;
    }

    if (pmod.ncoeff + order >= pmod.t2 - pmod.t1) {
    return GretlError.df;
    }

    newlist = gretl_list_new(pmod.list[0] + order);

    if (newlist == null) {
    err = GretlError.alloc;
    } else {
    newlist[0] = pmod.list[0] + order;
    for (i=2; i<=pmod.list[0]; i++) {
        newlist[i] = pmod.list[i];
    }
    if (dataset_add_series(dset, 1 + order)) {
        err = GretlError.alloc;
    }
    }

    if (!err) {
    /* add uhat to data set: substitute zeros for
       pre-sample values */
    for (t=0; t<n; t++) {
        if (t < pmod.t1) {
        dset.Z[v][t] = 0.0;
        } else {
        dset.Z[v][t] = pmod.uhat[t];
        }
    }
    strcpy(dset.varname[v], "uhat");
    series_set_label(dset, v, _("residual").ptr);
    /* then order lags of same */
    for (i=1; i<=order; i++) {
        int s, lv = v + i;
        double ul;

        sprintf(dset.varname[lv], "uhat_%d".ptr, i);
        newlist[pmod.list[0] + i] = lv;
        for (t=0; t<dset.n; t++) {
        s = t - i;
        if (s < 0) {
            dset.Z[lv][t] = 0.0;
        } else {
            ul = dset.Z[v][s];
            dset.Z[lv][t] = (na(ul))? 0.0 : ul;
        }
        }
    }
    }

    /* LMF apparatus: see Kiviet, Review of Economic Studies,
       53/2, 1986, equation (5), p. 245.
    */

    if (!err) {
    /* regression on [X~E], using original sample */
    impose_model_smpl(pmod, dset);
    newlist[1] = v;
    aux = lsq(newlist, dset, GretlCmdIndex.ols, Gretlopt.OPT_A);
    err = aux.errcode;
    if (err) {
        errmsg(err, prn);
    } else {
        RSSxe = aux.ess;
    }
    }

    if (!err) {
    int dfd = aux.nobs - pmod.ncoeff - order;
    int lberr;

    aux.aux = ModelAuxCode.ar;
    gretl_model_set_int(&aux, "BG_order", order);
    trsq = aux.rsq * aux.nobs;
    LMF = ((RSSx - RSSxe) / RSSxe) * dfd / order;
    pval = snedecor_cdf_comp(order, dfd, LMF);

    if (pmod.aux == ModelAuxCode.var || (opt & Gretlopt.OPT_I)) {
        ; /* don't print anything here */
    } else {
        if (opt & Gretlopt.OPT_Q) {
        bg_test_header(order, prn, 0);
        } else {
        printmodel(&aux, dset, Gretlopt.OPT_NONE, prn);
        pputc(prn, '\n');
        }
        pprintf(prn, "%s: LMF = %f,\n".ptr, _("Test statistic").ptr, LMF);
        pprintf(prn, "%s = P(F(%d,%d) > %g) = %.3g\n".ptr, _("with p-value").ptr,
            order, aux.nobs - pmod.ncoeff - order, LMF, pval);

        pprintf(prn, "\n%s: TR^2 = %f,\n".ptr,
            _("Alternative statistic").ptr, trsq);
        pprintf(prn, "%s = P(%s(%d) > %g) = %.3g\n\n".ptr, _("with p-value").ptr,
            _("Chi-square").ptr, order, trsq, chisq_cdf_comp(order, trsq));

        lb = ljung_box(order, pmod.t1, pmod.t2, dset.Z[v], &lberr);
        if (!na(lb)) {
        pprintf(prn, "Ljung-Box Q' = %g,\n".ptr, lb);
        pprintf(prn, "%s = P(%s(%d) > %g) = %.3g\n".ptr, _("with p-value").ptr,
            _("Chi-square").ptr, order, lb, chisq_cdf_comp(order, lb));
        }
        pputc(prn, '\n');
    }

    record_test_result(LMF, pval);

    if (opt & Gretlopt.OPT_S) {
        /* save the test onto @pmod */
        ModelTest *test = model_test_new(ModelTestType.autocorr);

        if (test != null) {
        model_test_set_teststat(test, TestStats.lmf);
        model_test_set_dfn(test, order);
        model_test_set_dfd(test, aux.nobs - pmod.ncoeff - order);
        model_test_set_order(test, order);
        model_test_set_value(test, LMF);
        model_test_set_pvalue(test, pval);
        maybe_add_test_to_model(pmod, test);
        }
    }
    }

    free(newlist);
    dataset_drop_last_variables(dset, dset.v - v);
    clear_model(&aux);

    /* reset sample as it was */
    dset.t1 = save_t1;
    dset.t2 = save_t2;

    return err;
}

/**
 * autocorr_test:
 * @pmod: pointer to model to be tested.
 * @order: lag order for test.
 * @dset: dataset struct.
 * @opt: if flags include OPT_S, save test results to model;
 * if OPT_Q, be less verbose; if OPT_I, be silent.
 * @prn: gretl printing struct.
 *
 * Tests the given model for autocorrelation of order equal to
 * the specified value, or equal to the frequency of the data if
 * the supplied @order is zero.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int autocorr_test(Model * pmod, int order, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    if (pmod.ci == GretlCmdIndex.ivreg) {
    return ivreg_autocorr_test(pmod, order, dset, opt, prn);
    } else if (pmod.ci == GretlCmdIndex.arma) {
    return lb_autocorr_test(pmod, order, opt, prn);
    } else if (gretl_is_regular_panel_model(pmod)) {
    if (dset.pd < 3) {
        /* time series not long enough */
        return GretlError.notimp;
    } else {
        return panel_autocorr_test(pmod, dset, opt, prn);
    }
    }

    if (pmod.ci != GretlCmdIndex.ols && pmod.ci != GretlCmdIndex.var) {
    return GretlError.notimp;
    } else if (model_has_missing_obs(pmod)) {
    return GretlError.missdata;
    }

    return regular_autocorr_test(pmod, order, dset, opt, prn);
}

int chow_active(int split, const(double) * x, int t)
{
    if (x != null) {
    return x[t] == 1.0;
    } else {
    return (t >= split);
    }
}

int * get_break_list(const(Model) * pmod, int ci, int * err)
{
    const(char) * lname = get_optval_string(cast(GretlCmdIndex) ci, Gretlopt.OPT_L);
    int * list = null;

    if (lname == null) {
    *err = GretlError.data;
    } else {
    list = get_list_by_name(lname);
    }

    if (list != null) {
    if (list[0] == 0) {
        *err = GretlError.data;
    } else {
        int i;

        for (i=1; i<=list[0] && !*err; i++) {
        if (!in_gretl_list(pmod.list, list[i])) {
            *err = GretlError.data;
        } else if (list[i] == 0) {
            *err = GretlError.data;
        }
        }
    }
    }

    return list;
}

/* compose list of variables to be used for the Chow test and add
   them to the data set */

int * make_chow_list(const(Model) * pmod, DatasetRaw * dset, int split, int dumv, int ci, Gretlopt opt, int * err)
{
    int * chowlist = null;
    int * brklist = null;
    int l0 = pmod.list[0];
    int ninter = 0, newvars = 0;
    int havedum = (dumv > 0);
    int i, t, v = dset.v;

    if (havedum && in_gretl_list(pmod.list, dumv)) {
    gretl_errmsg_sprintf(_("The model already contains %s").ptr,
                 dset.varname[dumv]);
    *err = GretlError.data;
    return null;
    }

    if (opt & Gretlopt.OPT_L) {
    /* --limit-to */
    brklist = get_break_list(pmod, ci, err);
    if (*err) {
        return null;
    } else {
        ninter = brklist[0];
    }
    } else {
    /* number of interaction terms */
    ninter = pmod.ncoeff - pmod.ifc;
    }

    newvars = ninter + 1 - havedum;

    if (dataset_add_series(dset, newvars)) {
    *err = GretlError.alloc;
    } else {
    chowlist = gretl_list_new(pmod.list[0] + ninter + 1);
    if (chowlist == null) {
        *err = GretlError.alloc;
    }
    }

    if (!*err) {
    const(double) * cdum = null;

    for (i=1; i<=l0; i++) {
        chowlist[i] = pmod.list[i];
    }

    if (dumv > 0) {
        /* we have a user-supplied dummy var */
        cdum = dset.Z[dumv];
    } else {
        /* generate the split variable */
        for (t=0; t<dset.n; t++) {
        dset.Z[v][t] = cast(double) (t >= split);
        }
        strcpy(dset.varname[v], "splitdum");
        series_set_label(dset, v, _("dummy variable for Chow test").ptr);
    }

    chowlist[l0 + 1] = (dumv > 0)? dumv : v;

    /* and the interaction terms */
    for (i=0; i<ninter; i++) {
        int pv, sv = v + i + 1 - havedum;

        if (brklist != null) {
        pv = brklist[i+1];
        } else {
        pv = pmod.list[i + 2 + pmod.ifc];
        }

        for (t=0; t<dset.n; t++) {
        if (model_missing(pmod, t)) {
            dset.Z[sv][t] = NADBL;
        } else if (chow_active(split, cdum, t)) {
            dset.Z[sv][t] = dset.Z[pv][t];
        } else {
            dset.Z[sv][t] = 0.0;
        }
        }

        if (havedum) {
        sprintf(dset.varname[sv], "%.2s_".ptr, dset.varname[dumv]);
        } else {
        strcpy(dset.varname[sv], "sd_");
        }
        strncat(dset.varname[sv], dset.varname[pv], VNAMELEN - 4);
        chowlist[l0 + 2 + i] = sv;
    }
    }

    return chowlist;
}

void write_plot_x_range(const(double) * x, int t1, int t2, FILE * fp)
{
    double xmin0, xmin, xmax0, xmax;
    double xrange;

    gretl_minmax(t1, t2, x, &xmin0, &xmax0);
    xrange = xmax0 - xmin0;
    xmin = xmin0 - xrange * .025;
    if (xmin0 >= 0.0 && xmin < 0.0) {
    xmin = 0.0;
    }
    xmax = xmax0 + xrange * .025;
    fprintf(fp, "set xrange [%.10g:%.10g]\n".ptr, xmin, xmax);
}

/* try to ensure that the critical value, @cv, is not squashed
   against the bottom or top of the plot
*/

void write_plot_y_range(const(double) * y, int T, double cv, FILE * fp)
{
    double ymin, ymax, yrange;
    double f = 0.025;

    gretl_minmax(0, T, y, &ymin, &ymax);
    if (cv < ymin) {
    ymin = cv;
    f = 0.05;
    } else if (cv > ymax) {
    ymax = cv;
    f = 0.05;
    }
    yrange = ymax - ymin;
    ymin -= f * yrange;
    ymax += f * yrange;

    fprintf(fp, "set yrange [%.10g:%.10g]\n".ptr, ymin, ymax);
}

int QLR_graph(const(double) * testvec, int t1, int t2, const(DatasetRaw) * dset, int df, int robust)
{
    const(double) * x = gretl_plotx(dset, Gretlopt.OPT_NONE);
    static const(char) *[2] titles = ["Chow F-test for break".ptr, "Robust Wald test for break".ptr];
    const(char) * title;
    alias QlrCritvalFn = extern (C) double function(int);
    QlrCritvalFn qlr_critval;
    double critval = NADBL;
    GptFlags flags = cast(GptFlags) 0;
    FILE * fp;
    int t;
    int err = 0;

    set_effective_plot_ci(GretlCmdIndex.qlrtest);

    if (dataset_is_time_series(dset)) {
    flags = GptFlags.ts | GptFlags.letterbox;
    }

    fp = open_plot_input_file(PlotType.regular, flags, &err);
    if (err) {
    return err;
    }

    qlr_critval = cast(QlrCritvalFn) get_plugin_function("qlr_critval_15_05");
    if (qlr_critval != null) {
    critval = qlr_critval(df);
    }

    print_keypos_string(GpKeySpec.left | GpKeySpec.top, fp);

    if (robust) {
    title = titles[1];
    } else {
    title = titles[0];
    if (!na(critval)) {
        critval /= df;
    }
    }

    gretl_push_c_numeric_locale();

    write_plot_x_range(x, t1, t2, fp);

    if (!na(critval)) {
    write_plot_y_range(testvec, t2 - t1, critval, fp);
    fprintf(fp, "plot \\\n'-' using 1:2 title \"%s\" w lines , \\\n%g title \"%s\" w lines dt 2\n".ptr,
        _(title), critval, _("5% QLR critical value").ptr);
    } else {
    fprintf(fp, "plot \\\n'-' using 1:2 title \"%s\" w lines\n".ptr, _(title));
    }

    for (t=t1; t<=t2; t++) {
    fprintf(fp, "%g %g\n".ptr, x[t], testvec[t-t1]);
    }
    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    err = finalize_plot_input_file(fp);
    reset_effective_plot_ci();

    return err;
}

void save_QLR_test(Model * pmod, const(char) * datestr, double X2, double pval, int df)
{
    ModelTest *test = model_test_new(ModelTestType.qlr);

    if (test != null) {
    model_test_set_teststat(test, TestStats.sup_wald);
    model_test_set_param(test, datestr);
    model_test_set_value(test, X2);
    model_test_set_pvalue(test, pval);
    model_test_set_dfn(test, df);
    maybe_add_test_to_model(pmod, test);
    }
}

/* for internal use by the "qlrtest" command */

double get_QLR_pval(double test, int df, int k1, int k2, Model * pmod)
{
    alias QlrAsyPvalueFn = extern (C) double function(double, int, double);
    QlrAsyPvalueFn qlr_asy_pvalue;
    double pi_1, pi_2, lam0;
    double pval;

    qlr_asy_pvalue = cast(QlrAsyPvalueFn) get_plugin_function("qlr_asy_pvalue");

    if (qlr_asy_pvalue == null) {
    return NADBL;
    }

    pi_1 = (k1 - pmod.t1 + 1) / cast(double) pmod.nobs;
    pi_2 = (k2 - pmod.t1 + 1) / cast(double) pmod.nobs;
    lam0 = (pi_2*(1.0 - pi_1)) / (pi_1*(1.0 - pi_2));

    pval = qlr_asy_pvalue(test, df, lam0);

    return pval;
}

/* for use by the qlrpval() function */

double QLR_pval(double X2, int df, double p1, double p2)
{
    alias QlrAsyPvalueFn = extern (C) double function(double, int, double);
    QlrAsyPvalueFn qlr_asy_pvalue;
    double lamda, pval;

    if (X2 < 0 || df <= 0 || p1 <= 0 || p2 <= p1 || p2 >= 1) {
    return NADBL;
    }

    qlr_asy_pvalue = cast(QlrAsyPvalueFn) get_plugin_function("qlr_asy_pvalue");

    if (qlr_asy_pvalue == null) {
    return NADBL;
    }

    lamda = (p2*(1.0 - p1)) / (p1*(1.0 - p2));
    pval = qlr_asy_pvalue(X2, df, lamda);

    return pval;
}

void QLR_print_result(double test, double pval, const(char) * datestr, int dfn, int dfd, int robust, Prn * prn)
{
    pputs(prn, _("Quandt likelihood ratio test for structural break at an unknown point,\nwith 15 percent trimming").ptr);
    pputs(prn, ":\n");
    if (robust) {
    pprintf(prn, _("The maximum Wald test = %g occurs at observation %s").ptr, test, datestr);
    } else {
    pprintf(prn, _("The maximum F(%d, %d) = %g occurs at observation %s").ptr, dfn, dfd, test, datestr);
    test *= dfn;
    }
    pputc(prn, '\n');

    if (!na(pval)) {
    pprintf(prn, _("Asymptotic p-value = %.6g for chi-square(%d) = %g").ptr,
        pval, dfn, test);
    pputc(prn, '\n');
    }
    pputc(prn, '\n');
}

void save_chow_test(Model * pmod, const(char) * chowstr, double test, double pval, int dfn, int dfd, Gretlopt opt)
{
    int ttype = (opt & Gretlopt.OPT_D)? ModelTestType.chowdum : ModelTestType.chow;
    ModelTest *mt = model_test_new(cast(ModelTestType) ttype);

    if (mt != null) {
    if (dfd == 0) {
        model_test_set_teststat(mt, TestStats.wald_chisq);
    } else {
        model_test_set_teststat(mt, TestStats.f);
    }
    model_test_set_param(mt, chowstr);
    model_test_set_value(mt, test);
    model_test_set_pvalue(mt, pval);
    model_test_set_dfn(mt, dfn);
    model_test_set_dfd(mt, dfd);
    maybe_add_test_to_model(pmod, mt);
    }
}

/*
 * real_chow_test:
 * @chowparm: sample breakpoint; or ID number of dummy
 * variable if given OPT_D; or ignored if given OPT_T.
 * @pmod: pointer to model to be tested.
 * @dset: dataset information.
 * @ci: either CHOW or QLRTEST.
 * @opt: if flags include OPT_S, save test results to model;
 * if OPT_D included, do the Chow test based on a given dummy
 * variable; if OPT_T included, do the QLR test. If OPT_M,
 * just save the test on the model, don't set $test and
 * $pvalue.
 * @prn: gretl printing struct.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int real_chow_test(int chowparm, Model * pmod, DatasetRaw * dset, int ci, Gretlopt opt, Prn * prn)
{
    int save_t1 = dset.t1;
    int save_t2 = dset.t2;
    int * chowlist = null;
    int origv = dset.v;
    Model chow_mod;
    int QLR = (opt & Gretlopt.OPT_T);
    int dumv = 0, split = 0, smax = 0;
    int err = 0;

    if (pmod.ci != GretlCmdIndex.ols) {
    return GretlError.olsonly;
    }

    if (exact_fit_check(pmod, prn)) {
    return 0;
    }

    /* temporarily impose the sample that was in force when the
       original model was estimated */
    impose_model_smpl(pmod, dset);

    gretl_model_init(&chow_mod, dset);

    if (QLR) {
    /* "15 percent trimming": exactly how this should be
       defined is debatable, but the following agrees
       with the strucchange package for R, when the
       Fstats() function is given "from=0.15, to=0.85".
    */
    split = cast(int) (pmod.t1 + floor(0.15 * pmod.nobs) - 1);
    smax = cast(int) (pmod.t1 + floor(0.85 * pmod.nobs) - 1);
    } else if (opt & Gretlopt.OPT_D) {
    /* Chow, using a predefined dummy */
    dumv = chowparm;
    } else {
    /* Chow, using break observation */
    smax = split = chowparm;
    }

    if (!err) {
    chowlist = make_chow_list(pmod, dset, split, dumv, ci,
                  opt, &err);
    }

    if (err) {
    goto bailout;
    }

    if (QLR) {
    /* Quandt likelihood ratio */
    int robust = (pmod.opt & Gretlopt.OPT_R);
    Gretlopt lsqopt = Gretlopt.OPT_A;
    double test, testmax = 0.0;
    double * testvec = null;
    int * testlist = null;
    int dfn = 0, dfd = 0;
    int nextra = dset.v - origv;
    int tmax = 0;
    int i, t;

    set_effective_plot_ci(GretlCmdIndex.qlrtest);
    if (gnuplot_graph_wanted(PlotType.regular, opt, null)) {
        testvec = cast(typeof(testvec)) malloc((smax - split + 1) * (*testvec).sizeof);
    }
    reset_effective_plot_ci();

    if (robust) {
        lsqopt |= Gretlopt.OPT_R;
        testlist = gretl_list_diff_new(chowlist, pmod.list, 2);
        if (testlist == null) {
        err = GretlError.alloc;
        }
    }

    for (t=split; t<=smax && !err; t++) {
        chow_mod = lsq(chowlist, dset, GretlCmdIndex.ols, lsqopt);
        if (chow_mod.errcode) {
        err = chow_mod.errcode;
        errmsg(err, prn);
        break;
        }
        dfn = chow_mod.ncoeff - pmod.ncoeff;
        if (robust) {
        test = wald_omit_chisq(testlist, &chow_mod);
        } else {
        dfd = chow_mod.dfd;
        test = (pmod.ess - chow_mod.ess) * dfd / (chow_mod.ess * dfn);
        }
        if (test > testmax) {
        tmax = t;
        testmax = test;
        }
        if (testvec != null) {
        testvec[t - split] = test;
        }
        clear_model(&chow_mod);
        for (i=0; i<nextra; i++) {
        dset.Z[origv+i][t] = 0.0;
        }
    }

    if (!err) {
        char[OBSLEN] datestr;
        double pval, X2;

        X2 = robust ? testmax : dfn * testmax;
        pval = get_QLR_pval(X2, dfn, split, smax, pmod);
        ntolabel(datestr.ptr, tmax, dset);

        if (!(opt & Gretlopt.OPT_Q)) {
        QLR_print_result(testmax, pval, datestr.ptr, dfn, dfd,
                 robust, prn);
        }
        if (!(opt & Gretlopt.OPT_M)) {
        record_QLR_test_result(X2, pval, tmax + 1);
        }
        if (opt & (Gretlopt.OPT_S | Gretlopt.OPT_M)) {
        save_QLR_test(pmod, datestr.ptr, X2, pval, dfn);
        }
        if (testvec != null) {
        QLR_graph(testvec, split, smax, dset, dfn, robust);
        }
    }

    free(testvec);
    free(testlist);
    } else {
    /* regular (or robust) Chow test */
    int robust = (pmod.opt & Gretlopt.OPT_R);
    Gretlopt lsqopt = Gretlopt.OPT_A;
    int * testlist = null;

    if (robust) {
        lsqopt |= Gretlopt.OPT_R;
    }
    chow_mod = lsq(chowlist, dset, GretlCmdIndex.ols, lsqopt);

    if (chow_mod.errcode) {
        err = chow_mod.errcode;
        errmsg(err, prn);
    } else if (chow_mod.ncoeff <= pmod.ncoeff) {
        err = chow_mod.errcode = GretlError.data;
        errmsg(err, prn);
    } else {
        int dfd = (robust)? 0 : chow_mod.dfd;
        int dfn = chow_mod.ncoeff - pmod.ncoeff;
        double test, pval = NADBL;

        if (!(opt & Gretlopt.OPT_Q)) {
        chow_mod.aux = ModelAuxCode.chow;
        printmodel(&chow_mod, dset, Gretlopt.OPT_NONE, prn);
        }

        if (robust) {
        testlist = gretl_list_diff_new(chow_mod.list, pmod.list, 2);
        if (testlist == null) {
            err = GretlError.alloc;
            goto bailout;
        } else {
            test = wald_omit_chisq(testlist, &chow_mod);
        }
        if (!na(test)) {
            pval = chisq_cdf_comp(dfn, test);
        }
        } else {
        test = (pmod.ess - chow_mod.ess) * dfd / (chow_mod.ess * dfn);
        if (!na(test)) {
            pval = snedecor_cdf_comp(dfn, dfd, test);
        }
        }

        if (na(test)) {
        pputs(prn, _("Couldn't compute Chow test statistic\n").ptr);
        } else if (na(pval)) {
        pprintf(prn, _("Couldn't compute Chow test p-value (test = %g)\n").ptr, test);
        }

        if (!na(test) && !na(pval)) {
        char[VNAMELEN] chowstr;

        if (opt & Gretlopt.OPT_Q) {
            pputc(prn, '\n');
        }
        if (opt & Gretlopt.OPT_D) {
            strcpy(chowstr.ptr, dset.varname[chowparm]);
            pprintf(prn, _("Chow test for structural difference with respect to %s").ptr,
                chowstr.ptr);
        } else {
            ntolabel(chowstr.ptr, chowparm, dset);
            pprintf(prn, _("Chow test for structural break at observation %s").ptr,
                chowstr.ptr);
        }
        pputc(prn, '\n');

        if (robust) {
            pprintf(prn, "  %s(%d) = %g %s %.4f\n".ptr, _("Chi-square").ptr,
                dfn, test, _("with p-value").ptr, pval);
            pprintf(prn, "  %s: F(%d, %d) = %g %s %.4f\n\n".ptr, _("F-form").ptr,
                dfn, chow_mod.dfd, test / dfn, _("with p-value").ptr,
                snedecor_cdf_comp(dfn, chow_mod.dfd, test / dfn));
        } else {
            pprintf(prn, "  F(%d, %d) = %g %s %.4f\n\n".ptr,
                dfn, dfd, test, _("with p-value").ptr, pval);
        }

        if (opt & Gretlopt.OPT_S) {
            save_chow_test(pmod, chowstr.ptr, test, pval, dfn, dfd, opt);
        }

        record_test_result(test, pval);
        }
    }
    clear_model(&chow_mod);
    free(testlist);
    }

 bailout:

    /* clean up extra variables */
    dataset_drop_last_variables(dset, dset.v - origv);
    free(chowlist);

    dset.t1 = save_t1;
    dset.t2 = save_t2;

    return err;
}

/**
 * chow_test:
 * @splitobs: the 0-based observation number at which to split
 * the sample.
 * @pmod: pointer to model to be tested.
 * @dset: dataset struct.
 * @opt: if flags include OPT_S, save test results to model.
 * @prn: gretl printing struct.
 *
 * Tests the given model for structural stability (Chow test)
 * using the sample break-point given by @splitobs
 * and prints the results to @prn. (The second portion of
 * the sample runs from observation @splitobs to the
 * end of the original sample range.)
 *
 * Returns: 0 on successful completion, error code on error.
 */

int chow_test(int splitobs, Model * pmod, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    int err = 0;

    if (splitobs <= 0 || splitobs >= dset.n) {
    gretl_errmsg_set(_("Invalid sample split for Chow test").ptr);
    err = GretlError.data;
    } else {
    err = real_chow_test(splitobs, pmod, dset, GretlCmdIndex.chow, opt, prn);
    }

    return err;
}

/**
 * chow_test_from_dummy:
 * @splitvar: the ID number of a dummy variable that should
 * be used to divide the sample.
 * @pmod: pointer to model to be tested.
 * @dset: dataset struct.
 * @opt: if flags include OPT_S, save test results to model.
 * @prn: gretl printing struct.
 *
 * Tests the given model for structural stability (Chow test),
 * using the dummy variable specified by @splitvar to divide the
 * original sample range of @pmod into two portions.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int chow_test_from_dummy(int splitvar, Model * pmod, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    int err = 0;

    if (splitvar <= 0 || splitvar >= dset.v) {
    gretl_errmsg_set(_("Invalid sample split for Chow test").ptr);
    err = GretlError.data;
    } else {
    err = real_chow_test(splitvar, pmod, dset, GretlCmdIndex.chow,
                 opt | Gretlopt.OPT_D, prn);
    }

    return err;
}

/**
 * QLR_test:
 * @pmod: pointer to model to be tested.
 * @dset: dataset struct.
 * @opt: if flags include OPT_S, save test results to model.
 * @prn: gretl printing struct.
 *
 * Tests the given model for structural stability via the
 * Quandt Likelihood Ratio test for a structural break at
 * an unknown point in the sample range.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int QLR_test(Model * pmod, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    return real_chow_test(0, pmod, dset, GretlCmdIndex.qlrtest,
              opt | Gretlopt.OPT_T, prn);
}

/* compute v'Mv, for symmetric M */

double vprime_M_v(double * v, double * M, int n)
{
    int i, j, jmin, jmax, k;
    double xx, val = 0.0;

    k = jmin = 0;
    for (i=0; i<n; i++) {
    xx = 0.0;
    for (j=jmin; j<n; j++) {
        xx += v[j] * M[k++];
    }
    val += xx * v[i];
    jmin++;
    }

    jmax = 1;
    for (i=1; i<n; i++) {
    k = i;
    xx = 0.0;
    for (j=0; j<jmax; j++) {
        xx += v[j] * M[k];
        k += n - j - 1;
    }
    val += xx * v[i];
    jmax++;
    }

    return val;
}

/* Compute scaled recursive residuals and cumulate \bar{w} via a
   sequence of OLS regressions. See William Greene, Econometric
   Analysis 5e, pp. 135-136.
*/

int cusum_compute(Model * pmod, double * cresid, int T, int k, double * wbar, DatasetRaw * dset)
{
    Model cmod;
    Gretlopt opt = Gretlopt.OPT_X | Gretlopt.OPT_A;
    double yhat;
    double den;
    double * xt;
    int n = T - k;
    int i, j, t;
    int err = 0;

    /* vector of regressors at t */
    xt = cast(typeof(xt)) malloc(k * (*xt).sizeof);
    if (xt == null) {
    return GretlError.alloc;
    }

    /* set minimal sample based on model to be tested */
    dset.t1 = pmod.t1;
    dset.t2 = pmod.t1 + k - 1;

    for (j=0; j<n && !err; j++) {
    cmod = lsq(pmod.list, dset, GretlCmdIndex.ols, opt);
    if (cmod.errcode) {
        err = cmod.errcode;
    } else {
        /* compute one step ahead prediction error: \hat{y}_t
           based on estimates through t - 1
        */
        t = dset.t2 + 1;
        yhat = 0.0;
        for (i=0; i<cmod.ncoeff; i++) {
        xt[i] = dset.Z[cmod.list[i+2]][t];
        yhat += cmod.coeff[i] * xt[i];
        }
        cresid[j] = dset.Z[pmod.list[1]][t] - yhat;
        cmod.ci = GretlCmdIndex.cusum;
        err = makevcv(&cmod, 1.0); /* (X'X)^{-1} */
        if (!err) {
        /* compute scaled residual: divide by
           \sqrt{1 + x_t'(X_{t-1}'X_{t-1})^{-1} x_t}
        */
        den = vprime_M_v(xt, cmod.vcv, cmod.ncoeff);
        cresid[j] /= sqrt(1.0 + den);
        *wbar += cresid[j];
        /* extend the sample by one observation */
        dset.t2 += 1;
        }
    }
    clear_model(&cmod);
    }

    free(xt);

    return err;
}

bool okfreq(int p)
{
    return (p == 1 || p == 4 || p == 12 || p == 24 || p == 52);
}

int cusum_do_graph(double a, double b, const(double) * W, int t1, int k, int m, DatasetRaw * dset, Gretlopt opt)
{
    FILE * fp = null;
    GptFlags flags = cast(GptFlags) 0;
    const(double) * obs = null;
    double frac = 1.0;
    double x0 = 0.0;
    int j, t, err = 0;

    if (dataset_is_time_series(dset)) {
    flags = GptFlags.ts | GptFlags.letterbox;
    }

    fp = open_plot_input_file(PlotType.cusum, flags, &err);
    if (err) {
    return err;
    }

    if (dataset_is_time_series(dset) && okfreq(dset.pd)) {
    b *= dset.pd;
    frac /= dset.pd;
        obs = gretl_plotx(dset, Gretlopt.OPT_NONE);
    if (obs != null) {
        x0 = obs[t1 + k];
    }
    }

    gretl_push_c_numeric_locale();

    fprintf(fp, "set xlabel '%s'\n".ptr, _("Observation").ptr);
    fputs("set nokey\n", fp);

    if (opt & Gretlopt.OPT_R) {
    fprintf(fp, "set title '%s'\n".ptr,
        /* xgettext:no-c-format */
        _("CUSUMSQ plot with 95% confidence band").ptr);
    fprintf(fp, "plot \\\n%g*(x-%g) title '' w dots lt 2, \\\n".ptr, b, x0 - frac);
    fprintf(fp, "%g+%g*(x-%g) title '' w lines lt 2, \\\n".ptr, -a, b, x0 - frac);
    fprintf(fp, "%g+%g*(x-%g) title '' w lines lt 2, \\\n".ptr, a, b, x0 - frac);
    } else {
    fputs("set xzeroaxis\n", fp);
    fprintf(fp, "set title '%s'\n".ptr,
        /* xgettext:no-c-format */
        _("CUSUM plot with 95% confidence band").ptr);
    fprintf(fp, "plot \\\n%g+%g*(x-%g) title '' w lines lt 2, \\\n".ptr, a, b, x0);
    fprintf(fp, "%g-%g*(x-%g) title '' w lines lt 2, \\\n".ptr, -a, b, x0);
    }

    fputs("'-' using 1:2 w linespoints lt 1\n", fp);

    for (j=0; j<m; j++) {
    t = t1 + k + j;
    if (obs != null) {
        fprintf(fp, "%g %g\n".ptr, obs[t], W[j]);
    } else {
        fprintf(fp, "%d %g\n".ptr, t, W[j]);
    }
    }

    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    return finalize_plot_input_file(fp);
}

void cusum_harvey_collier(double wbar, double sigma, int m, Model * pmod, Gretlopt opt, Prn * prn)
{
    double hct, pval;

    hct = (sqrt(cast(double) m) * wbar) / sigma;
    pval = student_pvalue_2(m - 1, hct);
    pprintf(prn, _("\nHarvey-Collier t(%d) = %g with p-value %.4g\n\n").ptr,
        m - 1, hct, pval);

    if (opt & Gretlopt.OPT_S) {
    ModelTest *test = model_test_new(ModelTestType.cusum);

    if (test != null) {
        model_test_set_teststat(test, TestStats.harvey_collier);
        model_test_set_dfn(test, m - 1);
        model_test_set_value(test, hct);
        model_test_set_pvalue(test, pval);
        maybe_add_test_to_model(pmod, test);
    }
    }

    record_test_result(hct, pval);
}

/**
 * cusum_test:
 * @pmod: pointer to model to be tested.
 * @dset: dataset struct.
 * @opt: if flags include %OPT_S, save results of test to model.
 * @prn: gretl printing struct.
 *
 * Tests the given model for parameter stability via the CUSUM test,
 * or if @opt includes %OPT_R, via the CUSUMSQ test; %OPT_Q makes
 * the test quiet; %OPT_U governs the associated plot, if wanted.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int cusum_test(Model * pmod, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    int save_t1 = dset.t1;
    int save_t2 = dset.t2;
    int T = pmod.nobs;
    int k = pmod.ncoeff;
    char[OBSLEN] cumdate;
    double wbar = 0.0;
    double * cresid = null, W = null;
    int quiet = opt & Gretlopt.OPT_Q;
    int m, i, j;
    int err = 0;

    if (pmod.ci != GretlCmdIndex.ols) {
    return GretlError.olsonly;
    } else if (exact_fit_check(pmod, prn)) {
    return 0;
    } else if (model_has_missing_obs(pmod)) {
    return GretlError.missdata;
    }

    /* the number of forecasts */
    m = T - k;

    cresid = cast(typeof(cresid)) malloc(m * (*cresid).sizeof);
    W = cast(typeof(W)) malloc(m * (*W).sizeof);

    if (cresid == null || W == null) {
    err = GretlError.alloc;
    }

    if (!err) {
    err = cusum_compute(pmod, cresid, T, k, &wbar, dset);
    if (err) {
        errmsg(err, prn);
    }
    }

    if (!err) {
    double a, b, den = 0.0, sigma = 0.0;
    int sig;

    if (opt & Gretlopt.OPT_R) {
        /* squares */
        double a1, a2, a3;
        double n = 0.5 * m - 1;

        pprintf(prn, "\n%s\n\n".ptr, _("CUSUMSQ test for stability of parameters").ptr);

        for (j=0; j<m; j++) {
        den += cresid[j] * cresid[j];
        }

        /* see Edgerton and Wells, Oxford Bulletin of Economics and
           Statistics, 56, 1994, pp. 355-365 */
        a1 = 1.358015 / sqrt(n);
        a2 = -0.6701218 / n;
        a3 = -0.8858694 / pow(n, 1.5);

        /* 0.5 * total height of band */
        a = a1 + a2 + a3;
        /* slope of expectation wrt time */
        b = 1.0 / m;
        if (!quiet) {
        pputs(prn, _("Cumulated sum of squared residuals").ptr);
        }
    } else {
        wbar /= T - k;
        pprintf(prn, "\n%s\n\n".ptr, _("CUSUM test for stability of parameters").ptr);
        pprintf(prn, _("mean of scaled residuals = %g\n").ptr, wbar);

        for (j=0; j<m; j++) {
        sigma += (cresid[j] - wbar) * (cresid[j] - wbar);
        }
        sigma /= T - k - 1;
        sigma = sqrt(sigma);
        pprintf(prn, _("sigmahat                 = %g\n\n").ptr, sigma);

        /* height of confidence band for first prediction */
        a = 0.948 * sqrt(cast(double) m);
        /* slope of confidence band limit wrt time */
        b = 2.0 * a / m;
        if (!quiet) {
        pputs(prn, _("Cumulated sum of scaled residuals").ptr);
        }
    }

    pputc(prn, '\n');
    pputs(prn, /* xgettext:no-c-format */
          _("('*' indicates a value outside of 95% confidence band)").ptr);
    pputs(prn, "\n\n");

    for (j=0; j<m; j++) {
        W[j] = 0.0;
        if (opt & Gretlopt.OPT_R) {
        for (i=0; i<=j; i++) {
            W[j] += cresid[i] * cresid[i] / den;
        }
        sig = fabs(W[j] - (j+1) / cast(double) m) > a;
        } else {
        for (i=0; i<=j; i++) {
            W[j] += cresid[i];
        }
        W[j] /= sigma;
        sig = fabs(W[j]) > a + j * b;
        }
        if (!quiet) {
        ntolabel(cumdate.ptr, pmod.t1 + k + j, dset);
        pprintf(prn, " %s %9.3f %s\n".ptr, cumdate.ptr, W[j], sig? "*".ptr : "".ptr);
        }
    }

    if (!(opt & Gretlopt.OPT_R)) {
        cusum_harvey_collier(wbar, sigma, m, pmod, opt, prn);
    }

    /* plot with 95% confidence band if wanted */
    if (gnuplot_graph_wanted(PlotType.cusum, opt, null)) {
        err = cusum_do_graph(a, b, W, pmod.t1, k, m, dset, opt);
    }
    }

    /* restore original sample range */
    dset.t1 = save_t1;
    dset.t2 = save_t2;

    free(cresid);
    free(W);

    return err;
}

/**
 * comfac_test:
 * @pmod: pointer to original model.
 * @dset: dataset struct.
 * @opt: if contains %OPT_S, save test results to model.
 * @prn: gretl printing struct.
 *
 * If @pmod was estimated via an AR(1) estimator, run an
 * auxiliary regression to test the implied common-factor
 * restriction.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int comfac_test(Model * pmod, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    Model cmod;
    int save_t1 = dset.t1;
    int save_t2 = dset.t2;
    int v = dset.v;
    int * biglist = null;
    int clearit = 0;
    int nadd, i, k, t;
    int err;

    if (pmod.ci != GretlCmdIndex.ar1 || (pmod.opt & Gretlopt.OPT_P)) {
    /* can't handle Prais-Winsten? */
    return GretlError.notimp;
    }

    /* check for changes in original list members */
    err = list_members_replaced(pmod, dset);
    if (err) {
    return err;
    }

    biglist = gretl_list_copy(pmod.list);
    if (biglist == null) {
    return GretlError.alloc;
    }

    nadd = 1 + pmod.ncoeff - pmod.ifc;

    err = dataset_add_series(dset, nadd);
    if (err) {
    free(biglist);
    return err;
    }

    /* add lags of the dependent variable and all regressors: some of
       these may be redundant but we'll just let the redundant terms be
       eliminated automatically via gretl's collinearity checking.
    */

    k = v;
    for (i=1; i<=pmod.list[0]; i++) {
    int src, lag, parent;

    src = pmod.list[i];
    if (src == 0) {
        continue;
    }
    for (t=0; t<dset.n; t++) {
        if (t == 0 || na(dset.Z[src][t-1])) {
        dset.Z[k][t] = NADBL;
        } else {
        dset.Z[k][t] = dset.Z[src][t-1];
        }
    }
    biglist = gretl_list_append_term(&biglist, k);
    if (biglist == null) {
        err = GretlError.alloc;
        break;
    }
    lag = series_get_lag(dset, src);
    parent = series_get_parent_id(dset, src);
    if (lag > 0 && parent > 0) {
        char[16] tmp;

        sprintf(tmp.ptr, "_%d".ptr, lag + 1);
        strcpy(dset.varname[k], dset.varname[parent]);
        gretl_trunc(dset.varname[k], 15 - strlen(tmp.ptr));
        strcat(dset.varname[k], tmp.ptr);
    } else {
        strcpy(dset.varname[k], dset.varname[src]);
        gretl_trunc(dset.varname[k], 13);
        strcat(dset.varname[k], "_1");
    }
    k++;
    }

    if (!err) {
    /* re-impose the sample that was in force when the original model
       was estimated */
    impose_model_smpl(pmod, dset);
    cmod = lsq(biglist, dset, GretlCmdIndex.ols, Gretlopt.OPT_A);
    clearit = 1;
    err = cmod.errcode;
    }

    if (!err) {
    if (cmod.nobs != pmod.nobs || cmod.ess > pmod.ess || cmod.dfd >= pmod.dfd) {
        /* something has gone wrong */
        err = GretlError.data;
    }
    }

    if (!err) {
    /* construct an F-test based on the SSR from the original
       AR(1) model and the SSR from the unrestricted model, cmod.
    */
    int dfd = cmod.dfd;
    int dfn = pmod.dfd - dfd - 1; /* account for rho */
    double SSRr = pmod.ess;
    double SSRu = cmod.ess;
    double Ftest = ((SSRr - SSRu)/dfn) / (SSRu/dfd);
    double pval = snedecor_cdf_comp(dfn, dfd, Ftest);

    if (!(opt & Gretlopt.OPT_I)) {
        if (!(opt & Gretlopt.OPT_Q)) {
        cmod.aux = ModelAuxCode.comfac;
        printmodel(&cmod, dset, Gretlopt.OPT_S, prn);
        pputc(prn, '\n');
        }

        pputs(prn, _("Test of common factor restriction").ptr);
        pputs(prn, "\n\n");

        pprintf(prn, "  %s: %s(%d, %d) = %g, ".ptr, _("Test statistic").ptr,
            "F".ptr, dfn, dfd, Ftest);
        pprintf(prn, _("with p-value = %g\n").ptr, pval);
        pputc(prn, '\n');
    }

    if (opt & Gretlopt.OPT_S) {
        ModelTest *test;

        test = model_test_new(ModelTestType.comfac);
        if (test != null) {
        model_test_set_teststat(test, TestStats.f);
        model_test_set_dfn(test, dfn);
        model_test_set_dfd(test, dfd);
        model_test_set_value(test, Ftest);
        model_test_set_pvalue(test, pval);
        maybe_add_test_to_model(pmod, test);
        }
    }

    record_test_result(Ftest, pval);
    }

    if (clearit) {
    clear_model(&cmod);
    }

    /* delete the added variables and restore the original
       sample range */

    dataset_drop_last_variables(dset, nadd);
    free(biglist);

    dset.t1 = save_t1;
    dset.t2 = save_t2;

    return err;
}

/**
 * panel_specification_test:
 * @pmod: pointer to model to be tested.
 * @dset: dataset struct.
 * @opt: option flags.
 * @prn: gretl printing struct.
 *
 * Tests the given pooled model for fixed and random effects.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int panel_specification_test(Model * pmod, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    if (pmod.ci != GretlCmdIndex.ols || !dataset_is_panel(dset)) {
    pputs(prn, _("This test is only relevant for pooled models\n").ptr);
    return 1;
    }

    if (pmod.ifc == 0) {
    pputs(prn, _("This test requires that the model contains a constant\n").ptr);
    return 1;
    }

    if (!balanced_panel(dset)) { /* ?? */
    pputs(prn, _("Sorry, can't do this test on an unbalanced panel.\nYou need to have the same number of observations\nfor each cross-sectional unit").ptr);
    return 1;
    } else {
    panel_diagnostics(pmod, dset, opt, prn);
    }

    return 0;
}

/* wrapper for backward compatibility */

int panel_hausman_test(Model * pmod, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    return panel_specification_test(pmod, dset, opt, prn);
}

/**
 * add_leverage_values_to_dataset:
 * @dset: dataset struct.
 * @m: matrix containing leverage values.
 * @opt: OPT_O = overwrite series on save.
 * @flags: may include SeriesSaveCode.leverage, SeriesSaveCode.influence,
 * and/or SeriesSaveCode.dffits.
 *
 * Adds to the working dataset one or more series calculated by
 * the gretl test for leverage/influence of data points.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int add_leverage_values_to_dataset(DatasetRaw * dset, GretlMatrixRaw * m, Gretlopt opt, int flags)
{
    int overwrite = (opt & Gretlopt.OPT_O);
    int t1, t2;
    int addvars = 0;

    if (flags & SeriesSaveCode.leverage) addvars++;
    if (flags & SeriesSaveCode.influence) addvars++;
    if (flags & SeriesSaveCode.dffits) addvars++;

    if (addvars == 0) {
    return 0;
    }

    if (dataset_add_series(dset, addvars)) {
    return GretlError.alloc;
    }

    t1 = gretl_matrix_get_t1(m);
    t2 = t1 + gretl_matrix_rows(m);

    /* add leverage? */
    if (flags & SeriesSaveCode.leverage) {
    int t, v = dset.v - addvars;
    int j = 0;

    for (t=0; t<dset.n; t++) {
        if (t < t1 || t >= t2) {
        dset.Z[v][t] = NADBL;
        } else {
        dset.Z[v][t] = gretl_matrix_get(m, j++, 0);
        }
    }
    strcpy(dset.varname[v], "lever");
    if (!overwrite) {
        make_varname_unique(dset.varname[v], v, dset);
    }
    series_set_label(dset, v, "leverage values");
    }

    /* add influence? */
    if (flags & SeriesSaveCode.influence) {
    int t, v = dset.v - (addvars - 1);
    int j = 0;

    for (t=0; t<dset.n; t++) {
        if (t < t1 || t >= t2) {
        dset.Z[v][t] = NADBL;
        } else {
        dset.Z[v][t] = gretl_matrix_get(m, j++, 1);
        }
    }
    strcpy(dset.varname[v], "influ");
    if (!overwrite) {
        make_varname_unique(dset.varname[v], v, dset);
    }
    series_set_label(dset, v, "influence values");
    }

    /* add DFFITS? */
    if (flags & SeriesSaveCode.dffits) {
    int t, v = dset.v - (addvars - 2);
    int j = 0;

    for (t=0; t<dset.n; t++) {
        double s, h;

        if (t < t1 || t >= t2) {
        dset.Z[v][t] = NADBL;
        } else {
        /* s = studentized residuals */
        h = gretl_matrix_get(m, j, 0);
        s = gretl_matrix_get(m, j, 2);
        if (na(h) || na(s)) {
            dset.Z[v][t] = NADBL;
        } else {
            dset.Z[v][t] = s * sqrt(h / (1.0 - h));
        }
        j++;
        }
    }
    strcpy(dset.varname[v], "dffits");
    if (!overwrite) {
        make_varname_unique(dset.varname[v], v, dset);
    }
    series_set_label(dset, v, "DFFITS values");
    }

    set_dataset_is_changed(dset, 1);

    return 0;
}

/**
 * leverage_test:
 * @pmod: pointer to model to be tested.
 * @dset: dataset struct.
 * @opt: if OPT_S, add calculated series to data set; operate
 * silently if OPT_Q. If includes OPT_O, do not adjust save
 * names (overwrite any existing series of the same names).
 * @prn: gretl printing struct.
 *
 * Tests the data used in the given model for points with
 * high leverage and influence on the estimates.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int leverage_test(Model * pmod, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    alias ModelLeverageFn = extern (C) GretlMatrixRaw * function(const(Model) *, const(DatasetRaw) *, Gretlopt, Prn *, int *);
    ModelLeverageFn model_leverage;
    GretlMatrixRaw * m;
    int err = 0;

    if (pmod.ci != GretlCmdIndex.ols) {
    return GretlError.olsonly;
    }

    model_leverage = cast(ModelLeverageFn) get_plugin_function("model_leverage");
    if (model_leverage == null) {
    return 1;
    }

    m = model_leverage(pmod, dset, opt, prn, &err);

    if (!err) {
    /* set the $results accessor */
    set_last_result_data(m, GretlType.matrix);
    if (opt & Gretlopt.OPT_S) {
        /* and respond to the --save option */
        err = add_leverage_values_to_dataset(dset, m, opt,
                         SeriesSaveCode.leverage |
                         SeriesSaveCode.influence |
                         SeriesSaveCode.dffits);
    }
    }

    return err;
}

/**
 * vif_test:
 * @pmod: pointer to model to be tested.
 * @dset: dataset struct.
 * @opt: may contain OPT_Q for quiet operation.
 * @prn: gretl printing struct.
 *
 * Calculates and displays the Variance Inflation Factors for
 * the independent variables in the given model.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int vif_test(Model * pmod, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    alias ComputeVifsFn = extern (C) int function(Model *, DatasetRaw *, Gretlopt, Prn *);
    ComputeVifsFn compute_vifs;

    gretl_error_clear();

    compute_vifs = cast(ComputeVifsFn) get_plugin_function("compute_vifs");
    if (compute_vifs == null) {
    return 1;
    }

    return compute_vifs(pmod, dset, opt, prn);
}

/**
 * bkw_test:
 * @pmod: pointer to model to be tested.
 * @dset: dataset struct.
 * @opt: may contain OPT_Q for quiet operation.
 * @prn: gretl printing struct.
 *
 * Calculates and displays the Belsley-Kuh-Welsch collinearity
 * diagnostics for @pmod.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int bkw_test(Model * pmod, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    alias ComputeBkwFn = extern (C) int function(Model *, DatasetRaw *, Gretlopt, Prn *);
    ComputeBkwFn compute_bkw;

    gretl_error_clear();

    compute_bkw = cast(ComputeBkwFn) get_plugin_function("compute_bkw");
    if (compute_bkw == null) {
    return 1;
    }

    return compute_bkw(pmod, dset, opt, prn);
}

