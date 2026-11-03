/* Port of lib/src/nls.c: estimation of nonlinear models by nonlinear least squares, maximum likelihood (with
   BFGS, Newton or Hessian-based variance) and the Gauss-Newton regression; the parameter/spec machinery shared
   with gmm.c is in nlspech.d, and the enumeration and prototypes of nls.h are nlsh.d. */

module gretl.nls;


import core.stdc.ctype : isdigit, isspace;
import core.stdc.errno : errno, ERANGE;
import core.stdc.float_ : DBL_EPSILON;
import core.stdc.math : ceil, exp, fabs, floor, fmod, isinf, isnan, log, log10, pow, sqrt;
import core.stdc.stdio : FILE, fclose, fopen, fputs, fprintf, printf, puts, sprintf, snprintf, sscanf, stderr;
import core.stdc.stdlib : abs, atoi, calloc, free, getenv, malloc, qsort, realloc;
import core.stdc.string : memcpy, memset, strcmp, strcpy, strlen, strchr, strstr, strcat, strncat, strncmp, strspn, strcspn, strtok, memmove, strrchr;

import std.algorithm : max, min;

import gretl.gretl_errorsh : GretlError;
import gretl.gretl_matrix;
import gretl.gretl_matrixh;
import gretl.johansenh : JohansenCode, JohansenInfo, jcode, jrank, auto_restr;
import gretl.libgretlh;
import gretl.missingh : NADBL, na;
import gretl.varh;

import gretl.xmlh;


private alias XUC = const(xmlChar) *;



import gretl.nlspech : NlFlags, Nlspec, Ocset, Parm;
import gretl.nlsh : CriterionType;
import gretl.minpackh : S_fp6, S_fp8;
import gretl.minpack.fdjac2 : fdjac2_;
import gretl.minpack.lmder : lmder_;
import gretl.minpack.lmdif : lmdif_;
import gretl.compare : QLR_test;
import gretl.dataio : get_date_x, ntolabel;
import gretl.dataset : create_auxiliary_dataset, dataset_add_series, dataset_drop_last_variables, destroy_dataset;
import gretl.dataseth : dataset_is_time_series;
import gretl.describe : gretl_mean, gretl_stddev;
import gretl.estimate : lsq, makevcv;
import gretl.geneval0 : Ufunc;
import gretl.genmain : destroy_genr, execute_genr, generate, genr_compile, genr_get_output_matrix, genr_get_output_type, genr_get_output_varnum, genr_is_autoregressive, genr_no_assign, genr_set_na_check, genr_unset_na_check, series_index;
import gretl.genparseh : Parser;
import gretl.gmem : g_free;
import gretl.gretl_bfgs : BFGS_max, hessian_inverse_from_score, newton_raphson_max, numerical_hessian_inverse, numerical_score_matrix, optimizer_get_matrix_name;
import gretl.gretl_bfgsh : BFGS_CRIT_FUNC, BFGS_GRAD_FUNC, HESS_FUNC;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_errors : errmsg, gretl_errmsg_append, gretl_errmsg_set, gretl_errmsg_sprintf;
import gretl.gretl_list : gretl_consecutive_list_new, gretl_list_copy, gretl_list_diff_new, gretl_list_sublist, in_gretl_list;
import gretl.gretl_model : clear_model, gretl_model_add_OPG_vcv, gretl_model_add_QML_vcv, gretl_model_add_hessian_vcv, gretl_model_allocate_storage, gretl_model_destroy_data_item, gretl_model_get_data, gretl_model_get_int, gretl_model_init, gretl_model_set_double, gretl_model_set_int, gretl_model_set_list_as_data, gretl_model_set_string_as_data, mle_criteria, set_model_id;
import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_get_buffer, gretl_print_new, gretl_print_steal_buffer, pprintf, pputc, pputs;
import gretl.gretl_prnh : PrnType;
import gretl.gretl_restrict : FN_NAMELEN;
import gretl.gretl_utils : copyvec, gretl_iszero, ijton, ls_criteria;
import gretl.gstrfuncs : g_strdup, g_strdup_printf, g_strfreev, g_strsplit;
import gretl.libset : gretl_in_gui_mode, gretl_iteration_pop, gretl_iteration_push, libset_get_double, libset_get_int;
import gretl.libseth : OptimCode, SetKey;
import gretl.nlspech : analytic_mode, matrix_deriv, numeric_mode, scalar_deriv, scalar_loglik, scalar_param, suppress_grad_check;
import gretl.options : get_optval_string, incompatible_options;
import gretl.printout : bufgets, bufgets_finalize, bufgets_init, print_iter_info;
import gretl.strutils : count_fields, equation_get_lhs_and_rhs, gretl_namechar_spn, gretl_strdup, gretl_strndup, gretl_word_strdup, shift_string_left, string_is_blank, strings_array_add, strings_array_free, strings_array_new, strings_array_new_with_length, tailstrip;
import gretl.tsls : ivreg_process_lists, tsls_make_endolist;
import gretl.uservarh : UserVar;
import gretl.gretl_array : get_array_by_name, gretl_array_get_strings, gretl_array_get_type;
import gretl.gretl_bundle : get_bundle_by_name, gretl_bundle_get_matrix, gretl_bundle_get_member_type, gretl_bundle_get_scalar, gretl_bundle_set_scalar;
import gretl.uservar : add_auxiliary_scalar, destroy_private_lists, destroy_private_matrices, destroy_private_scalars, get_string_by_name, get_user_var_by_name, gretl_is_scalar, gretl_scalar_add, gretl_scalar_get_value, gretl_scalar_set_value, private_matrix_add, set_auxiliary_scalars, unset_auxiliary_scalars, user_var_delete_by_name, user_var_get_type, user_var_get_value;
import gretl.genfuncs : resample_series;
import gretl.gretl_commands : gretl_command_number, gretl_command_word;
import gretl.usermat : get_matrix_by_name;

/* unported dependencies: per-file declarations from the 2026c headers */
extern (C) {
    int check_gmm_requirements(Nlspec * spec);
    int function_lookup(const(char) * s);
    Ufunc * get_user_function_by_name(const(char) * name);
    int gmm_add_vcv(Model * pmod, Nlspec * spec);
    int gmm_calculate(Nlspec * s, Prn * prn);
    int gmm_missval_check_etc(Nlspec * s);
    void maybe_add_gmm_residual(Model * pmod, const(Nlspec) * spec, const(DatasetRaw) * dset);
    int nlspec_add_ivreg_oc(Nlspec * s, int lhv, const(int) * rlist, const(double) ** Z);
    int nlspec_add_orthcond(Nlspec * s, const(char) * str, const(DatasetRaw) * dset);
    int nlspec_add_weights(Nlspec * s, const(char) * str);
    void nlspec_print_gmm_info(const(Nlspec) * spec, Prn * prn);
    void oc_set_destroy(Ocset * oc);
    int plausible_genr_start(const(char) * s, const(DatasetRaw) * dset);
}


/* file-scope global variables */

__gshared Nlspec private_spec;

void set_numeric_mode(Nlspec * s)
{
    s.flags &= ~NlFlags.analytical;
}

void set_analytic_mode(Nlspec * s)
{
    s.flags |= NlFlags.analytical;
}

void destroy_genrs_array(Nlspec * s)
{
    int i;

    for (i=0; i<s.ngenrs; i++) {
    destroy_genr(s.genrs[i]);
    }

    free(s.genrs);
    s.genrs = null;
    s.ngenrs = 0;
}

int check_lhs_vec(Nlspec * s)
{
    int v = gretl_vector_get_length(s.lvec);

    if (v != s.nobs) {
    if (v > 0) {
            s.t2 = v;
        s.nobs = v;
    } else {
        fprintf(stderr, "LHS vector should be of length %d, is %d\n".ptr,
            s.nobs, v);
        return 1;
    }
    }

    return 0;
}

GretlMatrixRaw * get_derivative_matrix(Nlspec * s, int i, int * err)
{
    GretlMatrixRaw * m = genr_get_output_matrix(s.params[i].dgenr);

    if (m == null) {
    fprintf(stderr, "get_derivative_matrix: got NULL result\n".ptr);
    *err = GretlError.data;
    }

    return m;
}

int check_derivative_matrix(int i, GretlMatrixRaw * m, Nlspec * s)
{
    int r, c, v;

    if (m == null) {
    fprintf(stderr, "param %d, got NULL matrix derivative\n".ptr, i);
    return 1;
    }

    r = gretl_matrix_rows(m);
    c = gretl_matrix_cols(m);
    v = s.params[i].nc;

    if (c != v || (r != 1 && r != s.nobs)) {
    fprintf(stderr, "matrix deriv for param %d is %d x %d: WRONG\n".ptr,
        i, r, c);
    fprintf(stderr, " should be %d x %d, or %d x %d\n".ptr, s.nobs,
        v, 1, v);
    return 1;
    }

    return 0;
}

int nls_dynamic_check(Nlspec * s, char * formula)
{
    Parser * genr;
    double * y;
    int T = s.dset.n;
    int v = s.dv;
    int err = 0;

    /* back up the dependent variable */
    y = copyvec(s.dset.Z[v], T);
    if (y == null) {
    return GretlError.alloc;
    }

    /* compile the formula for the dependent variable and see
       if it is autoregressive */
    strcpy(formula, s.nlfunc);
    adjust_saved_nlfunc(formula);

    genr = genr_compile(formula, s.dset, GretlType.series,
            Gretlopt.OPT_P | Gretlopt.OPT_N, null, &err);

    if (!err && genr_is_autoregressive(genr)) {
    fprintf(stderr, "nls_dynamic_check: setting NlFlags.autoreg\n".ptr);
    s.flags |= NlFlags.autoreg;
    }

    /* restore the dependent variable */
    memcpy(s.dset.Z[v], y, T * (*y).sizeof);

    destroy_genr(genr);
    free(y);

    return err;
}

int scalar_acceptable(Nlspec * s, int i, const(char) * dername)
{
    if (i < s.naux) {
    /* an auxiliary variable: scalars OK */
    return 1;
    } else if (i == s.naux) {
    /* the criterion value: accept a scalar for loglikelihood */
    return s.ci == GretlCmdIndex.mle;
    } else {
    /* i > s->naux: a derivative */
    return dername != null && gretl_is_scalar(dername);
    }
}

int allocate_generators(Nlspec * s)
{
    int np = analytic_mode(s) ? s.nparam : 0;
    int i, nlf = (s.nlfunc != null);
    int err = 0;

    s.ngenrs = s.naux + nlf + np;
    s.genrs = cast(typeof(s.genrs)) malloc(s.ngenrs * (*s.genrs).sizeof);

    if (s.genrs == null) {
    s.ngenrs = 0;
    err = GretlError.alloc;
    } else {
    for (i=0; i<s.ngenrs; i++) {
        s.genrs[i] = null;
    }
    }

    return err;
}

int unassigned_fncall(const(char) * s)
{
    int n = gretl_namechar_spn(s);
    int ret = 0;

    if (n > 0 && n < FN_NAMELEN && s[n] == '(') {
    char[FN_NAMELEN] word;

    *word.ptr = '\0';
    strncat(word.ptr, s, n);
    if (function_lookup(word.ptr) ||
        get_user_function_by_name(word.ptr)) {
        ret = 1;
    }
    }

    return ret;
}

void genr_setup_error(const(char) * s, int err, int exec)
{
    char * msg;

    if (exec) {
    msg = g_strdup_printf(_("The formula '%s'\n produced an error on execution").ptr, s);
    } else {
    msg = g_strdup_printf(_("The formula '%s'\n produced an error on compilation").ptr, s);
    }

    gretl_errmsg_append(msg, err);
    g_free(msg);
}

/* we "compile" the required equations first, so we can subsequently
   execute the compiled versions for maximum efficiency
*/

int nls_genr_setup(Nlspec * s)
{
    char[MAXLINE] formula;
    GretlMatrixRaw * m;
    int i, j, v;
    int err = 0;

    err = allocate_generators(s);
    if (err) {
    return err;
    }

    /* We now loop across the "generators", setting them up
       and checking them. We hook up any auxiliary genrs
       first, then the criterion function, then the analytical
       derivatives, if any.
    */

    /* initialize the index for derivatives */
    j = -1;

    for (i=0; i<s.ngenrs && !err; i++) {
    Gretlopt genopt = Gretlopt.OPT_P | Gretlopt.OPT_N;
    GretlType gentype = GretlType.any;
    GretlType result = cast(GretlType) 0;
    char * dname = null;

    if (i < s.naux) {
        /* auxiliary genrs */
        strcpy(formula.ptr, s.aux[i]);
        if (unassigned_fncall(formula.ptr)) {
        genopt |= Gretlopt.OPT_O;
        }
    } else if (i == s.naux) {
        /* criterion function */
        if (*s.lhname.ptr != '\0') {
        sprintf(formula.ptr, "%s = %s".ptr, s.lhname.ptr, s.nlfunc);
        } else {
        if (s.ci == GretlCmdIndex.nls && dataset_is_time_series(s.dset)) {
            err = nls_dynamic_check(s, formula.ptr);
        }
        sprintf(formula.ptr, "$nl_y = %s".ptr, s.nlfunc);
        }
    } else {
        /* derivative */
        sprintf(s.params[++j].dname.ptr, "$nl_x%d".ptr, i);
        if (scalar_param(s, j)) {
        sprintf(formula.ptr, "%s = %s".ptr, s.params[j].dname.ptr,
            s.params[j].deriv);
        } else {
        sprintf(formula.ptr, "%s = %s".ptr,  s.params[j].dname.ptr,
            s.params[j].deriv);
        gentype = GretlType.matrix;

        }
        dname = s.params[j].dname.ptr;
    }

    if (!err) {
        s.genrs[i] = genr_compile(formula.ptr, s.dset, gentype,
                       genopt, null, &err);
    }

    if (err) {
        fprintf(stderr, "nls: fail at genr_compile for genrs[%d]\n".ptr, i);
        fprintf(stderr, "> '%s'\n".ptr, formula.ptr);
        genr_setup_error(formula.ptr, err, 0);
        break;
    }

    /* see if the formula actually works, and flush out NAs
       while we're at it
    */
    genr_set_na_check(s.genrs[i]);
    err = execute_genr(s.genrs[i], s.dset, s.prn);
    genr_unset_na_check(s.genrs[i]);

    if (err) {
        genr_setup_error(formula.ptr, err, 1);
        fprintf(stderr, "nls: fail at execute_genr for genrs[%d]\n".ptr, i);
        fprintf(stderr, "> '%s'\n".ptr, formula.ptr);
        break;
    }

    /* skip ahead already? */
    if (genr_no_assign(s.genrs[i])) {
        continue;
    }

    v = 0;
    m = null;

    /* modified 2016-04-03: limit the following test(s) to
       the criterion function and derivatives
    */
    if (i >= s.naux) {
        v = genr_get_output_varnum(s.genrs[i]);
        m = genr_get_output_matrix(s.genrs[i]);
        if (v == 0 && m == null) {
        /* not a series, not a matrix: should be scalar */
        result = genr_get_output_type(s.genrs[i]);
        if (result != GretlType.double_ ||
            !scalar_acceptable(s, i, dname)) {
            gretl_errmsg_sprintf(_("The formula '%s'\n did not produce the required output type").ptr,
                     formula.ptr);
            err = GretlError.types;
            break;
        } else {
            gentype = GretlType.double_;
        }
        }
    }

    if (i == s.naux) {
        /* the criterion function */
        if (m != null) {
        s.lvec = m;
        s.lhtype = GretlType.matrix;
        err = check_lhs_vec(s);
        } else if (v > 0) {
        s.lhv = v;
        s.lhtype = GretlType.series;
        } else if (gentype == GretlType.double_) {
        s.lhtype = GretlType.double_;
        } else {
        err = GretlError.types;
        }
    } else if (j >= 0) {
        /* derivatives */
        if (v > 0) {
        s.params[j].dvtype = GretlType.series;
        } else if (m != null) {
        s.params[j].dvtype = GretlType.matrix;
        } else {
        s.params[j].dvtype = GretlType.double_;
        }

        s.params[j].dnum = v;
        s.params[j].dgenr = s.genrs[i];

        if (m != null || !scalar_param(s, j)) {
        err = check_derivative_matrix(j, m, s);
        }
    }

    }

    if (!err && s.hesscall != null) {
    s.hgen = genr_compile(s.hesscall, s.dset, GretlType.any,
                   Gretlopt.OPT_P | Gretlopt.OPT_O, s.prn, &err);
    if (err) {
        fprintf(stderr, "compilation of s->hgen failed\n".ptr);
    }
    }

    if (err) {
    destroy_genrs_array(s);
    }

    return err;
}

/* If i == 0 we're calculating the function; if i > 0 we're calculating
   a derivative.  Either way, we recalculate any auxiliary variables
   first.
*/

int nls_auto_genr(Nlspec * s, int i)
{
    int j;

    s.generr = 0;

    if (s.genrs == null) {
    gretl_iteration_push();
    s.generr = nls_genr_setup(s);
    gretl_iteration_pop();
    if (s.generr) {
        fprintf(stderr, " nls_genr_setup failed\n".ptr);
    }
    return s.generr;
    }

    for (j=0; j<s.naux; j++) {
    gretl_iteration_push();
    s.generr = execute_genr(s.genrs[j], s.dset, s.prn);
    gretl_iteration_pop();
    if (s.generr) {
        return s.generr;
    }
    }

    if (i == 0 && s.nlfunc == null) {
    /* we're done */
    return s.generr;
    }

    j = s.naux + i;
    s.generr = execute_genr(s.genrs[j], s.dset, s.prn);

    /* make sure we have a correct pointer to matrix deriv */
    if (!s.generr && i > 0 && matrix_deriv(s, i-1)) {
    GretlMatrixRaw * m = genr_get_output_matrix(s.genrs[j]);

    s.generr = check_derivative_matrix(i-1, m, s);
    }

    return s.generr;
}

/* wrappers for the above to enhance comprehensibility below */

int nl_calculate_fvec(Nlspec * s)
{
    return nls_auto_genr(s, 0);
}

int nls_calculate_deriv(Nlspec * s, int i)
{
    return nls_auto_genr(s, i + 1);
}

/* end wrappers */

void nlspec_destroy_arrays(Nlspec * s)
{
    free(s.params);
    s.params = null;
    s.nparam = 0;

    free(s.coeff);
    s.coeff = null;
    s.ncoeff = 0;
}

/* add an array of coefficients to the model specification */

int push_matrix_coeffs(Nlspec * s, GretlMatrixRaw * m, int k)
{
    double * coeff;
    int i, nc = s.ncoeff;

    if (k == 0) {
    return GretlError.data;
    }

    coeff = cast(typeof(coeff)) realloc(s.coeff, (nc + k) * (*coeff).sizeof);
    if (coeff == null) {
    return GretlError.alloc;
    }

    for (i=0; i<k; i++) {
    coeff[nc + i] = m.val[i];
    }

    s.coeff = coeff;
    s.ncoeff = nc + k;

    return 0;
}

/* add a scalar coefficient to the model specification */

int push_scalar_coeff(Nlspec * s, double x)
{
    double * coeff;
    int nc = s.ncoeff;

    coeff = cast(typeof(coeff)) realloc(s.coeff, (nc + 1) * (*coeff).sizeof);
    if (coeff == null) {
    return GretlError.alloc;
    }

    coeff[nc] = x;

    s.coeff = coeff;
    s.ncoeff = nc + 1;

    return 0;
}

double get_param_scalar(Parm * p)
{
    double x = NADBL;

    if (p.bundle != null) {
    int err = 0;

    x = gretl_bundle_get_scalar(p.bundle, p.name.ptr, &err);
    } else {
    x = gretl_scalar_get_value(p.name.ptr, null);
    }

    return x;
}

GretlMatrixRaw * get_param_matrix(Parm * p)
{
    GretlMatrixRaw * m = null;

    if (p.bundle != null) {
    int err = 0;

    m = gretl_bundle_get_matrix(p.bundle, p.name.ptr, &err);
    } else {
    m = get_matrix_by_name(p.name.ptr);
    }

    return m;
}

/* add a parameter to the model specification: this may be
   either a scalar or a matrix
*/

int nlspec_push_param(Nlspec * s, const(char) * name, GretlType type, GretlBundle * bundle, char * deriv)
{
    Parm * params, p;
    int np = s.nparam;
    int err;

    params = cast(typeof(params)) realloc(s.params, (np + 1) * (*params).sizeof);
    if (params == null) {
    return GretlError.alloc;
    }

    p = &params[np];

    p.name[0] = '\0';
    strncat(p.name.ptr, name, VNAMELEN - 1);
    p.type = type;
    p.bundle = bundle;
    p.dvtype = GretlType.none;
    p.deriv = deriv;
    p.dnum = 0;
    p.nc = 1;
    p.dname[0] = '\0';
    p.dgenr = null;
    p.mat = null;

    s.params = params;
    s.nparam = np + 1;

    if (type == GretlType.double_) {
    double x = get_param_scalar(p);

    err = push_scalar_coeff(s, x);
    } else {
    GretlMatrixRaw * m = get_param_matrix(p);
    int k = m.rows * m.cols;

    p.mat = m;
    p.nc = k;
    err = push_matrix_coeffs(s, m, k);
    if (!err) {
        s.nmat += 1;
    }
    }

    return err;
}

int check_matrix_param(const(char) * name, GretlBundle * b)
{
    GretlMatrixRaw * m;
    int err = 0;

    if (b == null) {
    m = get_matrix_by_name(name);
    } else {
    m = gretl_bundle_get_matrix(b, name, &err);
    }

    /* 2024-02-07: allow a matrix parameter which is not a vector,
       deleting the check clause 'gretl_vector_get_length(m) == 0'
    */
    if (err || gretl_is_null_matrix(m)) {
    gretl_errmsg_sprintf(_("'%s': expected a scalar or matrix").ptr, name);
    err = GretlError.types;
    }

    return err;
}

/* Do we have a valid name (the name of a scalar or matrix?):
   if so return 0, else return E_TYPES.
*/

int check_param_name(char ** pname, GretlType * type, GretlBundle ** pbundle)
{
    GretlBundle * b = null;
    char * name = *pname;
    int err = 0;

    if (pbundle != null && strchr(name, '.') != null) {
    /* did we get a bundle member? */
    char ** S = g_strsplit(name, ".", 2);

    b = get_bundle_by_name(S[0]);

    if (b != null) {
        GretlType btype;
        int berr = 0;

        btype = gretl_bundle_get_member_type(b, S[1], &berr);
        if (btype == GretlType.double_) {
        *type = GretlType.double_;
        } else if (btype == GretlType.matrix) {
        err = check_matrix_param(S[1], b);
        if (!err) {
            *type = GretlType.matrix;
        }
        } else {
        err = GretlError.types;
        }
        if (!err) {
        free(*pname);
        *pname = gretl_strdup(S[1]);
        *pbundle = b;
        }
    } else {
        err = GretlError.types;
    }
    g_strfreev(S);
    } else if (gretl_is_scalar(name)) {
    *type = GretlType.double_;
    } else {
    err = check_matrix_param(name, null);
    if (!err) {
        *type = GretlType.matrix;
    }
    }

    return err;
}

int nlspec_add_param_names(Nlspec * spec, const(char) * s)
{
    const(char) * test = null;
    int n, err = 0;

    s += cast(int) strspn(s, " ");

    if (*s == '"') {
    /* inline string literal */
    const(char) * p = strchr(s+1, '"');

    if (p == null) {
        err = GretlError.invarg;
    } else {
        n = cast(int) (p - s - 1);
        if (n > 0) {
        spec.parnames = gretl_strndup(s+1, n);
        } else {
        err = GretlError.invarg;
        }
        test = p + 1;
    }
    } else {
    /* name of string variable? */
    char[VNAMELEN] sname;

    n = gretl_namechar_spn(s);

    if (n == 0 || n >= VNAMELEN) {
        err = GretlError.invarg;
    } else {
        UserVar * uv = null;

        test = s + n;
        *sname.ptr = '\0';
        strncat(sname.ptr, s, n);
        uv = get_user_var_by_name(sname.ptr);
        if (uv == null) {
        err = GretlError.invarg;
        } else if (user_var_get_type(uv) == GretlType.string_) {
        /* got a composite string variable */
        const(char) * names = get_string_by_name(sname.ptr);

        if (names == null || *names == '\0') {
            err = GretlError.invarg;
        } else {
            free(spec.parnames);
            spec.parnames = gretl_strdup(names);
            if (spec.parnames == null) {
            err = GretlError.alloc;
            }
        }
        } else if (user_var_get_type(uv) == GretlType.array) {
        /* got an array of strings */
        GretlArray * a = cast(GretlArray *) user_var_get_value(uv);

        if (gretl_array_get_type(a) == GretlType.strings) {
            spec.parnames = gretl_strdup(sname.ptr);
            spec.flags |= NlFlags.names_array;
        } else {
            err = GretlError.invarg;
        }
        } else {
        err = GretlError.invarg;
        }
    }
    }

    if (!err && test != null && !string_is_blank(test)) {
    gretl_errmsg_set("found trailing junk on command line\n");
    err = GretlError.invarg;
    }

    return err;
}

/* For case where analytical derivatives are not given, the user
   must supply a line like:

   params b0 b1 b2 b3

   specifying the parameters to be estimated.  Here we parse such a
   list and add the parameter info to the spec.  The terms in the list
   must be pre-existing scalars or matrices, either objects in their
   own right or members of a bundle.
*/

int nlspec_add_params_from_line(Nlspec * s, const(char) * str)
{
    int i, nf = count_fields(str, null);
    int err = 0;

    if (s.params != null) {
    gretl_errmsg_set(_("Only one 'params' specification is allowed").ptr);
    return GretlError.data;
    } else if (nf == 0) {
    return GretlError.parse;
    }

    for (i=0; i<nf && !err; i++) {
    Gretlopt opt = Gretlopt.OPT_S | Gretlopt.OPT_D | Gretlopt.OPT_U;
    char * name = gretl_word_strdup(str, &str, opt, &err);
    GretlBundle * bundle = null;
    GretlType type = cast(GretlType) 0;

    if (!err) {
        err = check_param_name(&name, &type, &bundle);
    }

    if (!err) {
        err = nlspec_push_param(s, name, type,
                    bundle, null);
    }

    free(name);
    }

    if (err) {
    nlspec_destroy_arrays(s);
    }

    return err;
}

int nlspec_add_scalar_params(Nlspec * spec, int np, double * vals, char ** names, Gretlopt opt)
{
    int i, err = 0;

    if (spec.params != null || np == 0) {
    return GretlError.data;
    }

    for (i=0; i<np && !err; i++) {
    if (opt & Gretlopt.OPT_A) {
        /* doing internal auxiliary NLS (e.g. arma init) */
        err = add_auxiliary_scalar(names[i], vals[i]);
    } else {
        err = gretl_scalar_add(names[i], vals[i]);
    }
    if (!err) {
        err = nlspec_push_param(spec, names[i], GretlType.double_,
                    null, null);
    }
    }

    if (err) {
    nlspec_destroy_arrays(spec);
    }

    return err;
}

/**
 * nlspec_add_param_list:
 * @spec: nls specification.
 * @np: number of parameters.
 * @vals: array of initial parameter values.
 * @names: array of parameter names.
 *
 * Adds to @spec a list of (scalar) parameters to be estimated.
 * For an example of use see nls_example.c in the gretl extra
 * subdirectory.
 *
 * Returns: 0 on success, non-zero error code on error.
 */

int nlspec_add_param_list(Nlspec * spec, int np, double * vals, char ** names)
{
    return nlspec_add_scalar_params(spec, np, vals, names,
                    Gretlopt.OPT_NONE);

}

int aux_nlspec_add_param_list(Nlspec * spec, int np, double * vals, char ** names)
{
    return nlspec_add_scalar_params(spec, np, vals, names,
                    Gretlopt.OPT_A);
}

/* update the 'external' values of scalars or matrices using
   the values produced by the optimizer.
*/

int update_coeff_values(const(double) * b, Nlspec * s)
{
    int i, j, k = 0;
    int err = 0;

    for (i=0; i<s.nparam; i++) {
    Parm * p = &s.params[i];

    if (scalar_param(s, i)) {
        if (p.bundle != null) {
        err = gretl_bundle_set_scalar(p.bundle, p.name.ptr, b[k++]);
        } else {
        err = gretl_scalar_set_value(p.name.ptr, b[k++]);
        }
    } else {
        GretlMatrixRaw * m = get_param_matrix(p);

        if (m == null) {
        fprintf(stderr, "Couldn't find location for coeff %d\n".ptr, k);
        err = GretlError.data;
        } else {
        if (m != p.mat) {
            fprintf(stderr, "*** coeff_address: by name, '%s' is at %p; stored addr = %p\n".ptr, p.name.ptr,
                cast(void *) m, cast(void *) p.mat);
            p.mat = m;
        }
        for (j=0; j<p.nc; j++) {
            m.val[j] = b[k++];
        }
        }
    }
    }

    return err;
}

int nl_coeff_check(Nlspec * s)
{
    int i;

    for (i=0; i<s.ncoeff; i++) {
    if (na(s.coeff[i])) {
        gretl_errmsg_set("Uninitialized parameter");
        return GretlError.data;
    }
    }

    return 0;
}

/* Adjust starting and ending points of sample if need be, to avoid
   missing values; abort if there are missing values within the
   (possibly reduced) sample range.  For this purpose we generate the
   nls residual, or the loglikelihood in case of MLE.
*/

int nl_missval_check(Nlspec * s, const(DatasetRaw) * dset)
{
    int t1 = s.t1, t2 = s.t2;
    int t, v;
    int err = 0;

    /* calculate the function (NLS residual, MLE likelihood) */
    err = nl_calculate_fvec(s);
    if (err) {
    return err;
    }

    if (s.lvec != null || s.lhtype == GretlType.double_) {
    /* the calculation gives a matrix or scalar */
    goto nl_miss_exit;
    }

    /* ID number of LHS variable */
    v = s.lhv;

    for (t1=s.t1; t1<=s.t2; t1++) {
    if (!na(s.dset.Z[v][t1])) {
        break;
    }
    }

    for (t2=s.t2; t2>=t1; t2--) {
    if (!na(s.dset.Z[v][t2])) {
        break;
    }
    }

    if (t2 - t1 + 1 == 0) {
    fprintf(stderr, "nl_missval_check: no valid data\n".ptr);
    return GretlError.data;
    }

    if (t2 - t1 + 1 < s.ncoeff) {
    return GretlError.df;
    }

    for (t=t1; t<=t2; t++) {
    if (na(s.dset.Z[v][t])) {
        fprintf(stderr, "  after setting t1=%d, t2=%d, got NA for var %d (%s) at obs %d\n".ptr, t1, t2, v,
            dset.varname[v], t);
        return GretlError.missdata;
    }
    }

 nl_miss_exit:

    s.t1 = s.real_t1 = t1;
    s.t2 = s.real_t2 = t2;
    if (s.lhtype != GretlType.matrix) {
    s.nobs = t2 - t1 + 1;
    }

    return err;
}

/* get_mle_ll: callback used by BFGS.  Note that this should return
   NADBL in case of numerical problems; that signals to BFGS to
   try a smaller step length.
*/

extern (C) double get_mle_ll(const(double) * b, void * p)
{
    Nlspec * s = cast(Nlspec *) p;
    double x;
    int t, k, err;

    update_coeff_values(b, s);

    err = nl_calculate_fvec(s);
    if (err) {
    return NADBL;
    }

    if (s.lhtype == GretlType.double_) {
    s.crit = gretl_scalar_get_value(s.lhname.ptr, null);
    return s.crit;
    }

    s.crit = 0.0;

    if (s.lhtype == GretlType.matrix) {
    s.lvec = get_matrix_by_name(s.lhname.ptr);
    if (s.lvec == null) {
        fprintf(stderr, "get_mle_ll: s->lvec is gone!\n".ptr);
        return NADBL;
    }
    }

    if (s.lvec != null) {
    k = gretl_vector_get_length(s.lvec);
    for (t=0; t<k; t++) {
        x = s.lvec.val[t];
        if (na(x)) {
        s.crit = NADBL;
        break;
        }
        s.crit += x;
    }
    } else {
    k = s.lhv;
    for (t=s.t1; t<=s.t2; t++) {
        x = s.dset.Z[k][t];
        if (na(x)) {
        s.crit = NADBL;
        break;
        }
        s.crit += x;
    }
    }

    return s.crit;
}

/* this function is used in the context of the minpack callback, and
   also for checking derivatives in the MLE case
*/

int nl_function_calc(double * f, double * x, void * p)
{
    Nlspec * s = cast(Nlspec *) p;
    const(double) * y;
    int t, err = 0;

    /* calculate function given current parameter estimates */
    err = nl_calculate_fvec(s);
    if (err) {
    return err;
    }

    s.crit = 0.0;

    if (s.lvec != null) {
    y = s.lvec.val;
    } else {
    y = s.dset.Z[s.lhv] + s.t1;
    }

    /* transcribe from vector or series to array @f */
    for (t=0; t<s.nobs; t++) {
    if (na(y[t])) {
        fprintf(stderr, "nl_calculate_fvec: produced NA at obs %d\n".ptr, t+1);
        return 1;
    }
    f[t] = y[t];
    if (s.ci == GretlCmdIndex.mle) {
        s.crit += f[t];
    } else {
        /* sum of squares */
        s.crit += f[t] * f[t];
    }
    }

    s.iters += 1;

    if (s.ci == GretlCmdIndex.nls && (s.opt & Gretlopt.OPT_V)) {
    /* nls verbose output */
    print_iter_info(s.iters, s.crit, CriterionType.ssr, s.ncoeff,
            x, null, 0, s.prn);
    }

    return 0;
}

int get_nls_derivs(int T, double * g, DatasetRaw * gdset, void * p)
{
    Nlspec * spec = cast(Nlspec *) p;
    double * gi;
    double x;
    int j, t, k;
    int err = 0;

    if (g != null) {
    /* coming from &nls_calc, writing to flat array */
    gi = g;
    } else if (gdset != null) {
    /* coming from GNR, writing to a dataset */
    gi = gdset.Z[2];
    } else {
    return 1;
    }

    k = 0;

    for (j=0; j<spec.nparam && !err; j++) {
    if (nls_calculate_deriv(spec, j)) {
        return 1;
    }
    if (matrix_deriv(spec, j)) {
        GretlMatrixRaw * m = get_derivative_matrix(spec, j, &err);
        int i;

        for (i=0; !err && i<m.cols; i++) {
        x = gretl_matrix_get(m, 0, i);
        for (t=0; t<T; t++) {
            if (t > 0 && t < m.rows) {
            x = gretl_matrix_get(m, t, i);
            }
            gi[t] = (spec.ci == GretlCmdIndex.mle)? x : -x;
        }
        if (++k == spec.ncoeff) {
            break;
        } if (g != null) {
            gi += T;
        } else {
            gi = gdset.Z[k+2];
        }
        }
    } else if (scalar_deriv(spec, j)) {
        x = gretl_scalar_get_value(spec.params[j].dname.ptr, null);
        for (t=0; t<T; t++) {
        gi[t] = (spec.ci == GretlCmdIndex.mle)? x : -x;
        }
        if (++k == spec.ncoeff) {
        break;
        } else if (g != null) {
        gi += T;
        } else {
        gi = gdset.Z[k+2];
        }
    } else {
        /* derivative is series */
        int s, v = spec.params[j].dnum;

        if (v == 0) {
        /* FIXME? */
        v = spec.params[j].dnum = spec.dset.v - 1;
        }

        /* transcribe from dataset to array g */
        for (t=0; t<T; t++) {
        s = t + spec.t1;
        x = spec.dset.Z[v][s];
        gi[t] = (spec.ci == GretlCmdIndex.mle)? x : -x;
        }
        if (++k == spec.ncoeff) {
        break;
        } else if (g != null) {
        gi += T;
        } else {
        gi = gdset.Z[k+2];
        }
    }
    }

    return err;
}

/* for use with analytical derivatives, at present only for mle */

extern (C) int get_mle_gradient(double * b, double * g, int n, BFGS_CRIT_FUNC llfunc, void * p)
{
    Nlspec * spec = cast(Nlspec *) p;
    GretlMatrixRaw * m;
    double x;
    int i, j, k, t;
    int err = 0;

    update_coeff_values(b, spec);

    i = 0;

    for (j=0; j<spec.nparam && !err; j++) {
    if (nls_calculate_deriv(spec, j)) {
        return 1;
    }
    if (matrix_deriv(spec, j)) {
        m = get_derivative_matrix(spec, j, &err);
        for (k=0; !err && k<m.cols; k++) {
        g[i] = 0.0;
        for (t=0; t<m.rows; t++) {
            x = gretl_matrix_get(m, t, k);
            if (na(x)) {
            fprintf(stderr, "NA in gradient calculation\n".ptr);
            err = 1;
            } else {
            g[i] += x;
            }
        }
        i++;
        }
    } else if (scalar_deriv(spec, j)) {
        x = gretl_scalar_get_value(spec.params[j].dname.ptr, null);
        if (na(x)) {
        fprintf(stderr, "NA in gradient calculation\n".ptr);
        err = 1;
        } else {
        g[i++] = x;
        }
    } else {
        /* the derivative must be a series */
        int v = spec.params[j].dnum;

        g[i] = 0.0;
        for (t=spec.t1; t<=spec.t2; t++) {
        x = spec.dset.Z[v][t];
        if (na(x)) {
            fprintf(stderr, "NA in gradient calculation\n".ptr);
            err = 1;
        } else {
            g[i] += x;
        }
        }
        i++;
    }
    }

    return err;
}

extern (C) int get_mle_hessian(double * b, GretlMatrixRaw * H, void * p)
{
    Nlspec * spec = cast(Nlspec *) p;
    int k = H.rows;
    int err;

    if (b != null) {
    update_coeff_values(b, spec);
    }

    err = execute_genr(spec.hgen, spec.dset, spec.prn);

    if (!err) {
    GretlMatrixRaw * uH = get_matrix_by_name(spec.hname.ptr);

    if (uH == null) {
        err = GretlError.unkvar;
    } else if (uH.rows != k || uH.cols != k) {
        err = GretlError.nonconf;
    } else {
        gretl_matrix_copy_values(H, uH);
    }
    }

    return err;
}

GretlMatrixRaw * mle_hessian_inverse(Nlspec * spec, int * err)
{
    int k = spec.ncoeff;
    GretlMatrixRaw * H = gretl_matrix_alloc(k, k);

    if (H == null) {
    *err = GretlError.alloc;
    } else {
    *err = get_mle_hessian(null, H, spec);
    }

    if (!*err) {
    *err = gretl_invert_symmetric_matrix(H);
    if (*err) {
        fprintf(stderr, "mle_hessian_inverse: failed\n".ptr);
        gretl_matrix_free(H);
        H = null;
    }
    }

    return H;
}

/* Compute auxiliary statistics and add them to the NLS
   model struct */

void add_stats_to_model(Model * pmod, Nlspec * spec)
{
    int dv = spec.dv;
    int uncentered = 0;
    int k = spec.ncoeff;
    double * y = spec.dset.Z[dv];
    double tss = 0.0;
    double ypy = 0.0;
    double d;
    int s, t;

    if (gretl_model_get_int(pmod, "uncentered")) {
    uncentered = 1;
    }

    pmod.ess = spec.crit;
    pmod.sigma = sqrt(pmod.ess / (pmod.nobs - k));

    pmod.ybar = gretl_mean(pmod.t1, pmod.t2, y);
    pmod.sdy = gretl_stddev(pmod.t1, pmod.t2, y);

    s = (spec.missmask != null)? 0 : pmod.t1;

    for (t=pmod.t1; t<=pmod.t2; t++) {
    d = y[s] - pmod.ybar;
    tss += d * d;
    if (uncentered) {
        ypy += y[s] * y[s];
    }
    s++;
    }

    /* before over-writing the Gauss-Newton R^2, record it:
       it should be very small at convergence
    */
    gretl_model_set_double(pmod, "GNR_Rsquared", pmod.rsq);

    if (uncentered) {
    /* add uncentered R-squared */
    pmod.rsq = 1.0 - pmod.ess / ypy;
    } else if (tss == 0.0) {
    pmod.rsq = pmod.adjrsq = NADBL;
    }

    if (tss > 0) {
    double rsq = 1.0 - pmod.ess / tss;

    if (uncentered) {
        gretl_model_set_double(pmod, "centered-R2", rsq);
        pmod.adjrsq = NADBL;
    } else {
        double num = pmod.ess / (pmod.nobs - k);
        double den = tss / (pmod.nobs - 1);

        pmod.rsq = rsq;
        pmod.adjrsq = 1.0 - num / den;
    }
    }
}

/* returns the per-observation contributions to the log
   likelihood */

extern (C) const(double) * mle_llt_callback(const(double) * b, int i, void * p)
{
    Nlspec * s = cast(Nlspec *) p;
    int err;

    update_coeff_values(b, s);
    err = nl_calculate_fvec(s);

    if (err) {
    return null;
    } else if (s.lhtype == GretlType.matrix) {
    s.lvec = get_matrix_by_name(s.lhname.ptr);
    if (s.lvec == null) {
        fprintf(stderr, "mle_llt_callback: s->lvec is gone!\n".ptr);
        return null;
    } else {
        return s.lvec.val;
    }
    } else {
    return s.dset.Z[s.lhv] + s.t1;
    }
}

GretlMatrixRaw * ml_gradient_matrix(Nlspec * spec, int * err)
{
    GretlMatrixRaw * G = null;
    int k = spec.ncoeff;
    int T = spec.nobs;

    if (numeric_mode(spec)) {
    G = numerical_score_matrix(spec.coeff, T, k, &mle_llt_callback,
                   cast(void *) spec, err);
    } else {
    /* using analytical derivatives */
    GretlMatrixRaw * m;
    double x = 0.0;
    int i, j, v, s, t;

    if (spec.nparam == 1 && matrix_deriv(spec, 0)) {
        m = get_derivative_matrix(spec, 0, err);
        if (!*err) {
        G = gretl_matrix_copy(m);
        if (G == null) {
            *err = GretlError.alloc;
        }
        }
        return G;
    }

    G = gretl_matrix_alloc(T, k);
    if (G == null) {
        *err = GretlError.alloc;
        return null;
    }
    j = 0;
    for (i=0; i<spec.nparam; i++) {
        if (matrix_deriv(spec, i)) {
        m = get_derivative_matrix(spec, i, err);
        if (*err != 0) {
            break;
        }
        for (s=0; s<m.cols; s++) {
            x = gretl_matrix_get(m, 0, s);
            for (t=0; t<T; t++) {
            if (t > 0 && t < m.rows) {
                x = gretl_matrix_get(m, t, s);
            }
            gretl_matrix_set(G, t, j, x);
            }
            j++;
        }
        } else if (scalar_deriv(spec, i)) {
        x = gretl_scalar_get_value(spec.params[i].dname.ptr, null);
        for (t=0; t<T; t++) {
            gretl_matrix_set(G, t, j, x);
        }
        j++;
        } else {
        v = spec.params[i].dnum;
        for (t=0; t<T; t++) {
            x = spec.dset.Z[v][t + spec.t1];
            gretl_matrix_set(G, t, j, x);
        }
        j++;
        }
    }
    }

    return G;
}

Gretlopt ml_robust_specifier(Nlspec * spec)
{
    if (spec.opt & Gretlopt.OPT_C) {
    /* clustered */
    return Gretlopt.OPT_C;
    } else {
    const(char) * s = get_optval_string(GretlCmdIndex.mle, Gretlopt.OPT_R);

    if (s != null && !strcmp(s, "hac")) {
        return Gretlopt.OPT_N; /* Newey-West */
    }
    }

    return Gretlopt.OPT_NONE;
}

int mle_add_vcv(Model * pmod, Nlspec * spec)
{
    int err = 0;

    if (spec.opt & Gretlopt.OPT_A) {
    /* auxiliary model: no VCV */
    int i;

    for (i=0; i<pmod.ncoeff; i++) {
        pmod.sderr[i] = NADBL;
    }
    } else if (spec.opt & Gretlopt.OPT_H) {
    err = gretl_model_add_hessian_vcv(pmod, spec.Hinv);
    } else {
    GretlMatrixRaw * G = ml_gradient_matrix(spec, &err);

    if (!err) {
        if ((spec.opt & (Gretlopt.OPT_R | Gretlopt.OPT_C)) && spec.Hinv != null) {
        /* robust option: QML, possibly clustered or HAC */
        Gretlopt vopt = ml_robust_specifier(spec);

        err = gretl_model_add_QML_vcv(pmod, GretlCmdIndex.mle, spec.Hinv,
                          G, spec.dset, vopt, null);
        } else {
        err = gretl_model_add_OPG_vcv(pmod, G, null);
        }
    }
    gretl_matrix_free(G);
    }

    return err;
}

/* NLS: add coefficient covariance matrix and standard errors
   based on GNR */

int add_full_std_errs_to_model(Model * pmod)
{
    double abst, tstat_max = 0;
    int i, k, err = 0;

    if (pmod.vcv == null) {
    err = makevcv(pmod, pmod.sigma);
    if (err) {
        return err;
    }
    }

    for (i=0; i<pmod.ncoeff; i++) {
    k = ijton(i, i, pmod.ncoeff);
    if (pmod.vcv[k] < 0.0) {
        pmod.sderr[i] = NADBL;
    } else {
        pmod.sderr[i] = sqrt(pmod.vcv[k]);
        abst = fabs(pmod.coeff[i] / pmod.sderr[i]);
        if (abst > tstat_max) {
        tstat_max = abst;
        }
    }
    }

    if (tstat_max > 0) {
    gretl_model_set_double(pmod, "GNR_tmax", tstat_max);
    }

    return 0;
}

/* Experimental: watch out for bad stuff! Here we react
   to the case where there's machine-perfect collinearity
   in the Gauss-Newton Regression and one or more terms
   were dropped. We don't have standard errors for those
   terms but we do have coefficient estimates; we set
   the standard errors to NA. Note that @glist here is
   the full list of regressors passed to the GNR.
*/

int add_partial_std_errs_to_model(Model * pmod, const(int) * glist)
{
    int ndrop = glist[0] - pmod.list[0];
    double * coeff, sderr;
    int k;
    int * dlist;
    int i, j;

    if (ndrop <= 0) {
    /* eh? reinstate the error */
    fprintf(stderr, "no droplist found\n".ptr);
    return GretlError.jacobian;
    }

    dlist = gretl_list_diff_new(glist, pmod.list, 2);

    k = glist[0] - 1;
    coeff = cast(typeof(coeff)) malloc(k * (*coeff).sizeof);
    sderr = cast(typeof(sderr)) malloc(k * (*sderr).sizeof);

    if (coeff == null || sderr == null) {
    free(dlist);
    return GretlError.alloc;
    }

    j = 0;
    for (i=0; i<k; i++) {
    if (in_gretl_list(dlist, i+2)) {
        /* missing, skip it */
        sderr[i] = NADBL;
    } else {
        sderr[i] = pmod.sderr[j++];
    }
    }

    free(pmod.coeff);
    pmod.coeff = coeff;
    free(pmod.sderr);
    pmod.sderr = sderr;
    pmod.ncoeff = k;

    /* clean up stuff we don't want */
    gretl_model_destroy_data_item(pmod, "droplist");
    free(pmod.xpx);
    pmod.xpx = null;
    free(pmod.vcv);
    pmod.vcv = null;
    free(dlist);

    /* and record the breakage */
    gretl_model_set_int(pmod, "broken-vcv", 1);

    fprintf(stderr, "added partial standard errors\n".ptr);

    return 0;
}

int add_GNR_std_errs_to_model(Model * pmod, const(int) * list)
{
    int err;

    if (pmod.errcode == GretlError.jacobian) {
    pmod.errcode = err = add_partial_std_errs_to_model(pmod, list);
    } else {
    pmod.errcode = err = add_full_std_errs_to_model(pmod);
    }

    return err;
}

/* transcribe coefficient estimates into model struct */

void add_coeffs_to_model(Model * pmod, double * coeff)
{
    int i;

    for (i=0; i<pmod.ncoeff; i++) {
    pmod.coeff[i] = coeff[i];
    }
}

/* convert (back) from '%s - (%s)' to '%s = %s' */

char * adjust_saved_nlfunc(char * s)
{
    char * p = strchr(s, '-');

    if (p != null) {
    *p = '=';
    }

    p = strchr(s, '(');
    if (p != null) {
    shift_string_left(p, 1);
    }

    p = strrchr(s, ')');
    if (p != null) {
    *p = '\0';
    }

    return s;
}

 /* Attach additional specification info to make it possible to
    reconstruct the model from the saved state.  We need this
    for the "Modify model" option in the gretl GUI, and may
    also make use of it for bootstrapping NLS models.
 */

int nl_model_add_spec_info(Model * pmod, Nlspec * spec)
{
    const(char) * cmd = gretl_command_word(spec.ci);
    Prn * prn;
    char * buf;
    int i, err = 0;

    prn = gretl_print_new(PrnType.buffer, &err);
    if (err) {
    return err;
    }

    pputs(prn, cmd);

    if (pmod.depvar != null) {
    pprintf(prn, " %s\n".ptr, pmod.depvar);
    } else {
    pputc(prn, '\n');
    }

    if (spec.naux > 0) {
    gretl_model_set_int(pmod, "nl_naux", spec.naux);
    }

    for (i=0; i<spec.naux; i++) {
    pprintf(prn, "%s\n".ptr, spec.aux[i]);
    }

    if (spec.ci == GretlCmdIndex.gmm) {
    /* orthog, weights */
    nlspec_print_gmm_info(spec, prn);
    }

    if (numeric_mode(spec)) {
    pprintf(prn, "params".ptr);
    for (i=0; i<spec.nparam; i++) {
        pprintf(prn, " %s".ptr, spec.params[i].name.ptr);
    }
    pputc(prn, '\n');
    } else {
    for (i=0; i<spec.nparam; i++) {
        pprintf(prn, "deriv %s = %s\n".ptr, spec.params[i].name.ptr,
            spec.params[i].deriv);
    }
    }

    if (spec.hesscall != null) {
    pprintf(prn, "hessian %s\n".ptr, spec.hesscall);
    }

    pprintf(prn, "end %s\n".ptr, cmd);

    buf = gretl_print_steal_buffer(prn);
    gretl_model_set_string_as_data(pmod, "nlinfo", buf);
    gretl_print_destroy(prn);

    return 0;
}

int copy_user_parnames(Model * pmod, Nlspec * spec)
{
    const(char) * sep = ",";
    const(char) * pname;
    char * tmp;
    int i, err = 0;

    if (spec.flags & NlFlags.names_array) {
    /* we got an array of strings */
    GretlArray * a = get_array_by_name(spec.parnames);
    int ns = 0;

    if (a != null) {
        char ** S = gretl_array_get_strings(a, &ns);

        for (i=0; i<pmod.ncoeff && !err; i++) {
        if (i < ns) {
            pmod.params[i] = gretl_strdup(S[i]);
        } else {
            pmod.params[i] = gretl_strdup("unnamed");
        }
        }
    } else {
        err = GretlError.data;
    }
    return err;
    }

    /* copy the user-defined string before applying strtok */
    tmp = gretl_strdup(spec.parnames);
    if (tmp == null) {
    return GretlError.alloc;
    }

    if (strchr(tmp, ',') == null) {
    sep = " ";
    }

    for (i=0; i<pmod.ncoeff && !err; i++) {
    pname = strtok((i == 0)? tmp : null, sep);
    if (pname == null) {
        pmod.params[i] = gretl_strdup("unnamed");
    } else {
        pmod.params[i] = gretl_strdup(pname);
        if (pmod.params[i] == null) {
        err = GretlError.alloc;
        }
    }
    }

    free(tmp);

    return err;
}

/* Called for all of NLS, MLE, GMM: attach to the model struct
   the names of the parameters and some other string info
*/

int add_param_names_to_model(Model * pmod, Nlspec * spec)
{
    char[VNAMELEN] pname;
    int nc = pmod.ncoeff;
    int i, j, k, m, n;
    int err = 0;

    pmod.params = strings_array_new(nc);
    if (pmod.params == null) {
    return GretlError.alloc;
    }

    pmod.nparams = nc;

    if (spec.ci == GretlCmdIndex.nls) {
    pmod.depvar = gretl_strdup(spec.nlfunc);
    if (pmod.depvar != null) {
        adjust_saved_nlfunc(pmod.depvar);
    } else {
        err = GretlError.alloc;
    }
    } else if (spec.nlfunc != null) {
    n = cast(int) (strlen(spec.nlfunc) + strlen(spec.lhname.ptr));
    pmod.depvar = cast(typeof(pmod.depvar)) malloc(n + 4);
    if (pmod.depvar != null) {
        sprintf(pmod.depvar, "%s = %s".ptr, spec.lhname.ptr, spec.nlfunc);
    } else {
        err = GretlError.alloc;
    }
    }

    if (err) {
    free(pmod.params);
    return err;
    }

    if (spec.parnames != null) {
    /* handle the case where the user has given a "parnames"
       string or strings array
    */
    err = copy_user_parnames(pmod, spec);
    } else {
    /* compose automatic parameter names */
    i = 0;
    for (j=0; j<spec.nparam && !err; j++) {
        if (scalar_param(spec, j)) {
        pmod.params[i] = gretl_strdup(spec.params[j].name.ptr);
        if (pmod.params[i] == null) {
            err = GretlError.alloc;
        }
        i++;
        } else {
        m = spec.params[j].nc;
        sprintf(pname.ptr, "%d".ptr, m + 1);
        n = cast(int) (VNAMELEN - strlen(pname.ptr) - 3);
        for (k=0; k<m && !err; k++) {
            sprintf(pname.ptr, "%.*s[%d]".ptr, n, spec.params[j].name.ptr, k + 1);
            pmod.params[i] = gretl_strdup(pname.ptr);
            if (pmod.params[i] == null) {
            err = GretlError.alloc;
            }
            i++;
        }
        }
    }
    }

    if (!err && (spec.ci == GretlCmdIndex.nls || gretl_in_gui_mode())) {
    err = nl_model_add_spec_info(pmod, spec);
    }

    if (err && !pmod.errcode) {
    pmod.errcode = err;
    }

    return err;
}

int add_fit_resid_to_nls_model(Model * pmod, Nlspec * spec, int perfect)
{
    DatasetRaw * dset = spec.dset;
    int T = dset.n;
    int yno = spec.dv;
    double * tmp;
    int s, t;
    int err = 0;

    /* we need full-length arrays for uhat, yhat */

    tmp = cast(typeof(tmp)) realloc(pmod.uhat, T * (*tmp).sizeof);
    if (tmp == null) {
    return GretlError.alloc;
    } else {
    pmod.uhat = tmp;
    }

    tmp = cast(typeof(tmp)) realloc(pmod.yhat, T * (*tmp).sizeof);
    if (tmp == null) {
    return GretlError.alloc;
    } else {
    pmod.yhat = tmp;
    }

    /* OK, we got them, now transcribe */

    s = 0;
    for (t=0; t<T; t++) {
    if (t < pmod.t1 || t > pmod.t2) {
        pmod.uhat[t] = pmod.yhat[t] = NADBL;
    } else if (perfect) {
        pmod.uhat[t] = 0.0;
        pmod.yhat[t] = dset.Z[yno][t];
    } else {
        pmod.uhat[t] = spec.fvec[s];
        pmod.yhat[t] = dset.Z[yno][t] - spec.fvec[s];
        s++;
    }
    }

    if (perfect || (spec.flags & NlFlags.autoreg)) {
    pmod.rho = pmod.dw = NADBL;
    }

    return err;
}

/* this may be used later for generating out-of-sample forecasts --
   see nls_forecast() in forecast.c
*/

int transcribe_nls_function(Model * pmod, const(char) * s)
{
    char * formula;
    int err = 0;

    /* skip "depvar - " */
    s += cast(int) strcspn(s, " ") + 3;

    formula = gretl_strdup(s);
    if (s != null) {
    gretl_model_set_string_as_data(pmod, "nl_regfunc", formula);
    } else {
    err = GretlError.alloc;
    }

    return err;
}

int finalize_nls_model(Model * pmod, Nlspec * spec, int perfect, int * glist)
{
    DatasetRaw * dset = spec.dset;
    int err = 0;

    pmod.ci = spec.ci;
    pmod.t1 = spec.t1;
    pmod.t2 = spec.t2;
    pmod.full_n = dset.n;

    pmod.smpl.t1 = spec.dset.t1;
    pmod.smpl.t2 = spec.dset.t2;

    if (spec.lhtype == GretlType.matrix && spec.nobs != dset.n) {
        gretl_model_set_int(pmod, "matrix_data", 1);
    }

    err = add_GNR_std_errs_to_model(pmod, glist);

    add_stats_to_model(pmod, spec);

    if (!err) {
    add_fit_resid_to_nls_model(pmod, spec, perfect);
    }

    if (!err) {
    err = add_param_names_to_model(pmod, spec);
    }

    if (!err) {
    ls_criteria(pmod);
    pmod.fstt = pmod.chisq = NADBL;
    add_coeffs_to_model(pmod, spec.coeff);
    pmod.list[1] = spec.dv;

    /* set additional data on model to be shipped out */
    gretl_model_set_int(pmod, "iters", spec.iters);
    gretl_model_set_double(pmod, "tol", spec.tol);
    if (spec.nlfunc != null) {
        transcribe_nls_function(pmod, spec.nlfunc);
    }
    if (spec.flags & NlFlags.autoreg) {
        gretl_model_set_int(pmod, "dynamic", 1);
    }
    pmod.opt = spec.opt;
    }

    return err;
}

/* Gauss-Newton regression to calculate standard errors for the NLS
   parameters (see Davidson and MacKinnon).  This model is taken
   as the basis for the model struct returned by the "nls" command.
*/

Model GNR(int * glist, DatasetRaw * gdset, Gretlopt opt, Prn * prn)
{
    Gretlopt lsqopt = Gretlopt.OPT_A;
    Model gnr;

    if (opt & Gretlopt.OPT_R) {
    /* robust variance matrix, if wanted */
    lsqopt |= Gretlopt.OPT_R;
    }

    gnr = lsq(glist, gdset, GretlCmdIndex.ols, lsqopt);

    if (gnr.errcode) {
    pputs(prn, _("In Gauss-Newton Regression:\n").ptr);
    errmsg(gnr.errcode, prn);
    } else if (gnr.list[0] < glist[0]) {
    /* excessive collinearity */
    gnr.errcode = GretlError.jacobian;
    gretl_model_set_int(&gnr, "near-singular", 2);
    }

    return gnr;
}

DatasetRaw * make_GNR_dataset(Nlspec * spec, DatasetRaw * dset, int ** pglist, int * perfect, Prn * prn, int * err)
{
    double * uhat = spec.fvec;
    DatasetRaw * gdset = null;
    int * glist;
    int i, j, t, v;
    int T = spec.nobs;

    if (gretl_iszero(0, T - 1, uhat)) {
    pputs(prn, _("Perfect fit achieved\n").ptr);
    *perfect = 1;
    for (t=0; t<spec.nobs; t++) {
        uhat[t] = 1.0; /* will be adjusted later */
    }
    spec.crit = 0.0;
    }

    /* number of variables = const + depvar + spec->ncoeff
       (derivatives) */
    gdset = create_auxiliary_dataset(2 + spec.ncoeff, T, Gretlopt.OPT_NONE);
    glist = gretl_consecutive_list_new(1, spec.ncoeff + 1);

    if (gdset == null || glist == null) {
    destroy_dataset(gdset);
    free(glist);
    *err = GretlError.alloc;
    return null;
    }

    if (dataset_is_time_series(dset)) {
    gdset.structure = dset.structure;
    gdset.pd = dset.pd;
    ntolabel(gdset.stobs.ptr, dset.t1, dset);
    gdset.sd0 = get_date_x(gdset.pd, gdset.stobs.ptr);
    }

    /* dependent variable (NLS residual): write into
       slot 1 in gdset */
    strcpy(gdset.varname[1], "gnr_y");
    for (t=0; t<T; t++) {
    gdset.Z[1][t] = uhat[t];
    }

    /* independent variables: derivatives wrt NLS params,
       starting at slot 2 in gdset */
    for (i=0; i<spec.ncoeff; i++) {
    sprintf(gdset.varname[i+2], "gnr_x%d".ptr, i + 1);
    }
    if (analytic_mode(spec)) {
    get_nls_derivs(T, null, gdset, spec);
    } else {
    for (i=0; i<spec.ncoeff; i++) {
        v = i + 2;
        j = T * i; /* offset into jac array */
        for (t=0; t<T; t++) {
        gdset.Z[v][t] = spec.J.val[j++];
        }
    }
    }

    *pglist = glist;

    return gdset;
}

/* allocate space to copy info into model struct */

int nl_model_allocate(Model * pmod, Nlspec * spec)
{
    int k = spec.ncoeff;

    if (spec.opt & Gretlopt.OPT_A) {
    /* "auxiliary" model: no variance matrix, but
       we need sderr to prevent crash on printout
    */
    pmod.vcv = null;
    pmod.coeff = cast(typeof(pmod.coeff)) malloc(k * (*pmod.coeff).sizeof);
    pmod.sderr = cast(typeof(pmod.sderr)) malloc(k * (*pmod.sderr).sizeof);
    if (pmod.coeff == null || pmod.sderr == null) {
        pmod.errcode = GretlError.alloc;
    }
    } else {
    int nvc = (k * k + k) / 2;

    pmod.coeff = cast(typeof(pmod.coeff)) malloc(k * (*pmod.coeff).sizeof);
    pmod.sderr = cast(typeof(pmod.sderr)) malloc(k * (*pmod.sderr).sizeof);
    pmod.vcv = cast(typeof(pmod.vcv)) malloc(nvc * (*pmod.vcv).sizeof);

    if (pmod.coeff == null ||
        pmod.sderr == null ||
        pmod.vcv == null) {
        pmod.errcode = GretlError.alloc;
    }
    }

    if (pmod.errcode == 0) {
    pmod.ncoeff = k;
    }

    return pmod.errcode;
}

/* Work up the results of estimation into the form of a gretl
   MODEL struct.  This is used for MLE and GMM; in the case of
   NLS the Gauss-Newton regression provides the basis for the
   final_ model structure.
*/

int make_other_nl_model(Model * pmod, Nlspec * spec, const(DatasetRaw) * dset)
{
    nl_model_allocate(pmod, spec);
    if (pmod.errcode) {
    return pmod.errcode;
    }

    pmod.t1 = spec.t1;
    pmod.t2 = spec.t2;
    pmod.nobs = spec.nobs;

    if (spec.lhtype == GretlType.matrix && spec.nobs != dset.n) {
        gretl_model_set_int(pmod, "matrix_data", 1);
    }

    /* hmm */
    pmod.dfn = pmod.ncoeff;
    pmod.dfd = pmod.nobs - pmod.ncoeff;

    pmod.ci = spec.ci;

    if (spec.ci == GretlCmdIndex.mle) {
    pmod.lnL = spec.crit;
    } else {
    /* GMM */
    pmod.lnL = NADBL;
    }

    mle_criteria(pmod, 0);
    add_coeffs_to_model(pmod, spec.coeff);

    if (spec.ci == GretlCmdIndex.gmm && (spec.opt & Gretlopt.OPT_G)) {
    ; /* IVREG via GMM */
    } else {
    pmod.errcode = add_param_names_to_model(pmod, spec);
    }

    if (!pmod.errcode) {
    if (pmod.ci == GretlCmdIndex.mle) {
        pmod.errcode = mle_add_vcv(pmod, spec);
    } else {
        pmod.errcode = gmm_add_vcv(pmod, spec);
    }
    }

    if (!pmod.errcode && pmod.ci == GretlCmdIndex.gmm) {
    maybe_add_gmm_residual(pmod, spec, dset);
    }

    if (!pmod.errcode) {
    if (spec.flags & NlFlags.newton) {
        gretl_model_set_int(pmod, "iters", spec.fncount);
    } else {
        gretl_model_set_int(pmod, "fncount", spec.fncount);
        gretl_model_set_int(pmod, "grcount", spec.grcount);
        gretl_model_set_double(pmod, "tol", spec.tol);
    }
    }

    /* mask invalid stats */
    pmod.sigma = NADBL;
    pmod.rsq = pmod.adjrsq = NADBL;
    pmod.fstt = pmod.chisq = NADBL;
    pmod.rho = pmod.dw = NADBL;

    return pmod.errcode;
}

int add_nls_coeffs(Model * pmod, Nlspec * spec)
{
    pmod.ncoeff = spec.ncoeff;
    pmod.full_n = 0;

    pmod.errcode = gretl_model_allocate_storage(pmod);

    if (!pmod.errcode) {
    add_coeffs_to_model(pmod, spec.coeff);
    }

    return pmod.errcode;
}

/* free up resources associated with the nlspec struct */

void clear_nlspec(Nlspec * spec)
{
    int i;

    if (spec == null) {
    return;
    }

    free(spec.parnames);
    spec.parnames = null;

    if (spec.params != null) {
    for (i=0; i<spec.nparam; i++) {
        free(spec.params[i].deriv);
    }
    free(spec.params);
    spec.params = null;
    }
    spec.nparam = 0;

    free(spec.fvec);
    spec.fvec = null;

    gretl_matrix_free(spec.J);
    spec.J = null;

    free(spec.coeff);
    spec.coeff = null;
    spec.ncoeff = 0;
    spec.nmat = 0;

    if (spec.aux != null) {
    for (i=0; i<spec.naux; i++) {
        free(spec.aux[i]);
    }
    free(spec.aux);
    spec.aux = null;
    }
    spec.naux = 0;

    if (spec.genrs != null) {
    for (i=0; i<spec.ngenrs; i++) {
        destroy_genr(spec.genrs[i]);
    }
    free(spec.genrs);
    spec.genrs = null;
    }
    spec.ngenrs = 0;
    spec.generr = 0;

    if (spec.hgen != null) {
    destroy_genr(spec.hgen);
    spec.hgen = null;
    }

    if (spec.hesscall != null) {
    free(spec.hesscall);
    spec.hesscall = null;
    }

    free(spec.nlfunc);
    spec.nlfunc = null;

    gretl_matrix_free(spec.Hinv);
    spec.Hinv = null;

    spec.flags = cast(NlFlags) 0;
    spec.opt = Gretlopt.OPT_NONE;

    spec.dv = 0;
    spec.lhtype = GretlType.none;
    spec.lhv = 0;
    spec.lvec = null;

    spec.lhname[0] = '\0';
    spec.hname[0] = '\0';

    spec.iters = 0;
    spec.fncount = 0;
    spec.grcount = 0;

    spec.t1 = spec.t2 = 0;
    spec.real_t1 = spec.real_t2 = 0;
    spec.nobs = 0;

    spec.dset = null;
    spec.prn = null;

    if (spec.oc != null) {
    oc_set_destroy(spec.oc);
    spec.oc = null;
    }

    free(spec.missmask);
    spec.missmask = null;
}

/*
   Next block: functions that interface with minpack.

   The details below may be obscure, but here's the basic idea: The
   minpack functions are passed an array ("fvec") that holds the
   calculated values of the function to be minimized, at a given value
   of the parameters, and also (in the case of analytic derivatives)
   an array holding the Jacobian ("jac").  Minpack is also passed a
   callback function that will recompute the values in these arrays,
   given a revised vector of parameter estimates.

   As minpack does its iterative thing, at each step it invokes the
   callback function, supplying its updated parameter estimates and
   saying via a flag variable ("iflag") whether it wants the function
   itself or the Jacobian re-evaluated.

   The libgretl strategy involves holding all the relevant values (nls
   residual, nls derivatives, and nls parameters) as variables in a
   gretl dataset (Z-array and datainfo-struct pair).  The callback
   function that we supply to minpack first transcribes the revised
   parameter estimates into the dataset, then invokes genr() to
   recalculate the residual and derivatives, then transcribes the
   results back into the fvec and jac arrays.
*/

/* callback for lm_calculate (below) to be used by minpack */

extern (C) int nls_calc(int m, int n, double * x, double * fvec, double * jac, int ldjac, int * iflag, void * p)
{
    Nlspec * s = cast(Nlspec *) p;
    int err;

    /* write current coefficient values into dataset */
    update_coeff_values(x, s);

    if (*iflag == 1) {
    /* calculate function at x, results into fvec */
    err = nl_function_calc(fvec, x, p);
    if (err) {
        fprintf(stderr, "nl_function_calc: err = %d\n".ptr, err);
        *iflag = -1;
    }
    } else if (*iflag == 2) {
    /* calculate jacobian at x, results into jac */
    err = get_nls_derivs(m, jac, null, p);
    if (err) {
        fprintf(stderr, "get_nls_derivs: err = %d\n".ptr, err);
        *iflag = -1;
    }
    }

    return 0;
}

/* Below: copied here from minpack chkder.c, with a view to
   figuring if it could be improved (lessening the chance
   of false alarms for correct derivatives).
*/

int chkder(const(double) * x, double * xp, const(double) * fvec, const(double) * fvecp, const(GretlMatrixRaw) * J, int mode, double * err)
{
    const(double) epsmch = DBL_EPSILON;
    double temp, eps = sqrt(epsmch);
    int i, j;

    if (mode == 1) {
    /* mode 1: construct a neighboring vector, xp */
    for (j=0; j<J.cols; j++) {
        temp = eps * fabs(x[j]);
        if (temp == 0.0) {
        temp = eps;
        }
        xp[j] = x[j] + temp;
    }
    } else {
    /* mode 2: assess validity of gradient */
    const(double) factor = 100;
    double d, epsf = factor * epsmch;
    double epslog = log10(eps);

    for (i=0; i<J.rows; i++) {
        err[i] = 0.0;
    }
    for (j=0; j<J.cols; j++) {
        temp = fabs(x[j]);
        if (temp == 0.0) {
        temp = 1.0;
        }
        for (i=0; i<J.rows; i++) {
        err[i] += temp * gretl_matrix_get(J, i, j);
        }
    }
    for (i=0; i<J.rows; i++) {
        temp = 1.0;
        d = fabs(fvecp[i] - fvec[i]);
        if (fvec[i] != 0.0 && fvecp[i] != 0.0 &&
        d >= epsf * fabs(fvec[i])) {
        d = fabs((fvecp[i] - fvec[i]) / eps - err[i]);
        temp = eps * d / (fabs(fvec[i]) + fabs(fvecp[i]));
        }
        err[i] = 1.0;
        if (temp > epsmch && temp < eps) {
        err[i] = (log10(temp) - epslog) / epslog;
        }
        if (temp >= eps) {
        err[i] = 0.0;
        }
    }
    }

    return 0;
}

/* in case the user supplied analytical derivatives for the
   parameters, check them for sanity */

int check_derivatives(Nlspec * spec, Prn * prn)
{
    double * x = spec.coeff;
    double * fvec = spec.fvec;
    double * jac = spec.J.val;
    int m = spec.nobs;
    int n = spec.ncoeff;
    int ldjac = m;
    int iflag;
    double * xp, xerr, fvecp;
    int i, badcount = 0, zerocount = 0;

    /* note: allocate space for xerr and fvecp too */
    xp = cast(typeof(xp)) malloc((n + m + m) * (*xp).sizeof);
    if (xp == null) {
    return GretlError.alloc;
    }

    xerr = xp + n;
    fvecp = xerr + m;

    /* mode 1: x contains the point of evaluation of the function; on
       output xp is set to a neighboring point. */
    chkder(x, xp, fvec, fvecp, spec.J, 1, xerr);

    /* calculate gradient */
    iflag = 2;
    nls_calc(m, n, x, fvec, jac, ldjac, &iflag, spec);
    if (iflag == -1) goto chkderiv_abort;

    /* calculate function, at neighboring point xp */
    iflag = 1;
    nls_calc(m, n, xp, fvecp, jac, ldjac, &iflag, spec);
    if (iflag == -1) goto chkderiv_abort;

    /* mode 2: on input, fvec must contain the functions, the rows of
       fjac must contain the gradients evaluated at x, and fvecp must
       contain the functions evaluated at xp.  On output, xerr contains
       measures of correctness of the respective gradients.
    */
    chkder(x, xp, fvec, fvecp, spec.J, 2, xerr);

    /* examine "xerr" vector */
    for (i=0; i<m; i++) {
    if (xerr[i] == 0.0) {
        zerocount++;
    } else if (xerr[i] < 0.35) {
        badcount++;
    }
    }

    if (zerocount > 0) {
    gretl_errmsg_set(_("NLS: The supplied derivatives seem to be incorrect").ptr);
    fprintf(stderr, "%d out of %d tests gave zero\n".ptr, zerocount, cast(int) m);
    } else if (badcount > 0) {
    pputs(prn, _("Warning: The supplied derivatives may be incorrect, or perhaps\nthe data are ill-conditioned for this function.\n").ptr);
    pprintf(prn, _("%d out of %d gradients looked suspicious.\n\n").ptr,
        badcount, cast(int) m);
    }

 chkderiv_abort:

    free(xp);

    return (zerocount > m / 4);
}

/* drivers for BFGS code below, first for MLE */

int mle_calculate(Nlspec * s, Prn * prn)
{
    BFGS_GRAD_FUNC gradfunc = null;
    HESS_FUNC hessfunc = null;
    int maxit, use_newton = 0;
    int err = 0;

    if (libset_get_int(SetKey.gretl_optim) == OptimCode.newton) {
    use_newton = 1;
    }

    if (analytic_mode(s) && !suppress_grad_check(s)) {
    err = check_derivatives(s, prn);
    }

    if (!err) {
    if (analytic_mode(s)) {
        gradfunc = &get_mle_gradient;
    }
    if (s.hesscall != null) {
        hessfunc = &get_mle_hessian;
        s.flags |= NlFlags.ahess;
        if (!(s.opt & (Gretlopt.OPT_G | Gretlopt.OPT_R | Gretlopt.OPT_C))) {
        /* default to Hessian unless we got a
           conflicting option flag */
        s.opt |= Gretlopt.OPT_H;
        }
    }
    }

    if (!err && use_newton) {
    double crittol = 1.0e-7;
    double gradtol = 1.0e-7;

    maxit = 100;
    s.flags |= NlFlags.newton;
    err = newton_raphson_max(s.coeff, s.ncoeff, maxit,
                 crittol, gradtol, &s.fncount,
                 CriterionType.loglik, &get_mle_ll,
                 gradfunc, hessfunc, s,
                 s.opt, s.prn);
    } else if (!err) {
    maxit = 500;
    err = BFGS_max(s.coeff, s.ncoeff, maxit, s.tol,
               &s.fncount, &s.grcount,
               &get_mle_ll, CriterionType.loglik, gradfunc, s,
               null, s.opt, s.prn);
    }

    if (!err && (s.opt & (Gretlopt.OPT_H | Gretlopt.OPT_R | Gretlopt.OPT_C))) {
    /* doing Hessian or QML covariance matrix */
    if (hessfunc != null) {
        s.Hinv = mle_hessian_inverse(s, &err);
    } else {
        /* Note 2018-10-05: changed the condition for using
           hessian_inverse_from_score(): we were requiring
           both analytic_mode and that the loglikelihood
           function returns per-observation values (not just
           a scalar). But it seems the latter requirement,
           !scalar_loglik(s), is not really necessary.
        */
        if (analytic_mode(s)) {
        s.Hinv = hessian_inverse_from_score(s.coeff, s.ncoeff,
                             gradfunc, &get_mle_ll,
                             s, &err);
        } else {
        double d = 0.01; /* adjust? */

        s.Hinv = numerical_hessian_inverse(s.coeff, s.ncoeff,
                            &get_mle_ll, s, d, &err);
        }
    }

    if (err && !scalar_loglik(s)) {
        pprintf(prn, _("\nError: Hessian non-negative definite? (err = %d); dropping back to OPG\n").ptr, err);
        s.opt = cast(Gretlopt)(s.opt & ~Gretlopt.OPT_H);
        s.opt = cast(Gretlopt)(s.opt & ~Gretlopt.OPT_R);
        err = 0;
    }
    }

    return err;
}

/* Minpack driver for the case where analytical derivatives have been
   supplied */

int lm_calculate(Nlspec * spec, Prn * prn)
{
    int m = spec.nobs;
    int n = spec.ncoeff;
    int lwa = 5 * n + m; /* work array size */
    int ldjac = m; /* leading dimension of jac array */
    int info = 0;
    int nfev = 0;
    int njev = 0;
    int * ipvt;
    double * wa;
    double factor;
    double ftol, xtol, gtol;
    int mode, maxfev, nprint;
    int err = 0;

    wa = cast(typeof(wa)) malloc(lwa * (*wa).sizeof);
    ipvt = cast(typeof(ipvt)) malloc(n * (*ipvt).sizeof);

    if (wa == null || ipvt == null) {
    err = GretlError.alloc;
    goto nls_cleanup;
    }

    if (!suppress_grad_check(spec)) {
    err = check_derivatives(spec, prn);
    if (err) {
        goto nls_cleanup;
    }
    }

    /* mostly use the lmder1() defaults */
    maxfev = 100 * (n + 1);
    nprint = 0;
    mode = 1;
    ftol = xtol = spec.tol;
    gtol = 0.0;
    factor = 100; /* default */

    if (spec.flags & NlFlags.smallstep) {
    /* try to ensure a shorter step-length */
    factor = 1.0;
    }

    /* call minpack */
    lmder_(&nls_calc, m, n, spec.coeff, spec.fvec, spec.J.val, ldjac,
       ftol, xtol, gtol, maxfev, wa, mode, factor, nprint,
       &info, &nfev, &njev, ipvt, wa + n, wa + 2*n, wa + 3*n,
       wa + 4*n, wa + 5*n, spec);

    switch (info) {
    case -1:
    err = 1;
    break;
    case 0:
    gretl_errmsg_set(_("Invalid NLS specification").ptr);
    err = 1;
    break;
    case 1:
    case 2:
    case 3:
    case 4:
    break;
    case 5:
    case 6:
    case 7:
    gretl_errmsg_sprintf(_("NLS: failed to converge after %d iterations").ptr,
                 spec.iters);
    err = 1;
    break;
    default:
    break;
    }

 nls_cleanup:

    free(wa);
    free(ipvt);

    return err;
}

/* callback for lm_approximate (below) to be used by minpack */

extern (C) int nls_calc_approx(int m, int n, double * x, double * fvec, int * iflag, void * p)
{
    int erru, errc, err;

    /* write current parameter values into dataset Z */
    err = erru = update_coeff_values(x, cast(Nlspec *) p);

    /* calculate function at x, results into fvec */
    if (!err) {
    err = errc = nl_function_calc(fvec, x, p);
    }

    if (err) {
    /* flag error to minpack */
    fprintf(stderr, "nls_calc_approx: got error %d from %s\n".ptr,
        err, (erru)? "update_coeff_values".ptr :
        "nl_function_calc".ptr);
    *iflag = -1;
    }

    return 0;
}

/* Minpack driver for the case where the Jacobian must be approximated
   numerically */

int lm_approximate(Nlspec * spec, Prn * prn)
{
    int m = spec.nobs;
    int n = spec.ncoeff;
    int ldjac = m;
    int info = 0;
    int maxfev = 200 * (n + 1); /* max iterations */
    int mode = 1, nprint = 0, nfev = 0;
    int iflag = 0;
    int * ipvt;
    double gtol = 0.0;
    double epsfcn = 0.0, factor = 100.;
    double * dspace, diag, qtf;
    double * wa1, wa2, wa3, wa4;
    int err = 0;

    dspace = cast(typeof(dspace)) malloc((5 * n + m) * (*dspace).sizeof);
    ipvt = cast(typeof(ipvt)) malloc(n * (*ipvt).sizeof);

    if (dspace == null || ipvt == null) {
    err = GretlError.alloc;
    goto nls_cleanup;
    }

    diag = dspace;
    qtf = diag + n;
    wa1 = qtf + n;
    wa2 = wa1 + n;
    wa3 = wa2 + n;
    wa4 = wa3 + n;

    /* call minpack */
    lmdif_(&nls_calc_approx, m, n, spec.coeff, spec.fvec,
       spec.tol, spec.tol, gtol, maxfev, epsfcn, diag,
       mode, factor, nprint, &info, &nfev, spec.J.val, ldjac,
       ipvt, qtf, wa1, wa2, wa3, wa4, spec);

    spec.iters = nfev;

    switch (cast(int) info) {
    case -1:
    err = 1;
    break;
    case 0:
    gretl_errmsg_set(_("Invalid NLS specification").ptr);
    err = 1;
    break;
    case 1:
    case 2:
    case 3:
    case 4:
    break;
    case 5:
    case 6:
    case 7:
    case 8:
    gretl_errmsg_sprintf(_("NLS: failed to converge after %d iterations").ptr,
                 spec.iters);
    err = 1;
    break;
    default:
    break;
    }

    if (!err) {
    double ess = spec.crit;
    int iters = spec.iters;
    Gretlopt opt = spec.opt;

    spec.opt = Gretlopt.OPT_NONE;

    /* call minpack again */
    fdjac2_(&nls_calc_approx, m, n, 0, spec.coeff, spec.fvec,
        spec.J.val, ldjac, &iflag, epsfcn, wa4, spec);
    spec.crit = ess;
    spec.iters = iters;
    spec.opt = opt;
    }

 nls_cleanup:

    free(dspace);
    free(ipvt);

    return err;
}

/* below: various public functions */

/**
 * nlspec_add_param_with_deriv:
 * @spec: pointer to nls specification.
 * @s: string specifying a derivative with respect to a
 *   parameter of the regression function.
 *
 * Adds an analytical derivative to @spec.  This pointer must
 * have previously been obtained by a call to nlspec_new().
 * The required format for @dstr is "%varname = %formula", where
 * %varname is the name of the (scalar) variable holding the parameter
 * in question, and %formula is an expression, of the sort that
 * is fed to gretl's %genr command, giving the derivative of the
 * regression function in @spec with respect to the parameter.
 * The variable holding the parameter must be already present in
 * the dataset.
 *
 * Returns: 0 on success, non-zero error code on error.
 */

int nlspec_add_param_with_deriv(Nlspec * spec, const(char) * s)
{
    const(char) * p = s;
    char * name = null;
    char * deriv = null;
    GretlBundle * b = null;
    GretlType type = cast(GretlType) 0;
    int err = 0;

    if (spec.ci == GretlCmdIndex.gmm) {
    gretl_errmsg_set(_("Analytical derivatives cannot be used with GMM").ptr);
    return GretlError.data;
    }

    if (!strncmp(p, "deriv ", 6)) {
    /* make starting with "deriv" optional */
    p += 6;
    }

    err = equation_get_lhs_and_rhs(p, &name, &deriv);
    if (err) {
    fprintf(stderr, "parse error in deriv string: '%s'\n".ptr, s);
    return GretlError.parse;
    }

    err = check_param_name(&name, &type, &b);

    if (!err) {
    err = nlspec_push_param(spec, name, type, b, deriv);
    if (err) {
        free(deriv);
        deriv = null;
    }
    }

    free(name);

    if (!err) {
    set_analytic_mode(spec);
    }

    return err;
}

int nlspec_add_hessian(Nlspec * spec, const(char) * hesscall)
{
    int err;

    if (spec.ci != GretlCmdIndex.mle) {
    /* we only do this for MLE */
    return GretlError.types;
    }

    if (*spec.hname.ptr != '\0') {
    /* Hessian already added */
    return GretlError.data;
    }

    err = optimizer_get_matrix_name(hesscall, spec.hname.ptr);

    if (!err) {
    spec.hesscall = gretl_strdup(hesscall);
    if (spec.hesscall == null) {
        err = GretlError.alloc;
    }
    }

    return err;
}

int screen_bad_aux(const(char) * line, const(DatasetRaw) * dset, int * pci)
{
    int n = gretl_namechar_spn(line);
    char[FN_NAMELEN] word;
    int ci, err = GretlError.data;

    *word.ptr = '\0';

    if (n > 0) {
    if (n > FN_NAMELEN - 1) {
        n = FN_NAMELEN - 1;
    }
    strncat(word.ptr, line, n);
    }

    *pci = ci = gretl_command_number(word.ptr);

    if (ci == GretlCmdIndex.genr || ci == GretlCmdIndex.print || ci == GretlCmdIndex.printf || ci == GretlCmdIndex.eval) {
    err = 0;
    } else if (plausible_genr_start(line, dset)) {
    err = 0;
    } else if (function_lookup(word.ptr)) {
    err = 0;
    } else if (get_user_function_by_name(word.ptr)) {
    err = 0;
    } else if (ci > 0) {
    gretl_errmsg_sprintf(_("command '%s' not valid in this context").ptr,
                 word.ptr);
    } else {
    gretl_errmsg_sprintf(_("'%s': not valid in this context").ptr,
                 word.ptr);
    }

    return err;
}

/* convert printf from command to function form, if needed */

char * revise_aux_printf(const(char) * s)
{
    char * tmp = null;

    if (!strncmp(s, "printf ", 7)) {
    tmp = g_strdup_printf("printf(%s)".ptr, s + 7);
    } else {
    tmp = g_strdup(s);
    }

    return tmp;
}

/**
 * nlspec_add_aux:
 * @spec: pointer to nls specification.
 * @s: string specifying an auxiliary command (primarily
 * for use in calculating function or derivatives).
 * @dset: pointer to dataset information.
 *
 * Adds the specification of an auxiliary command to @spec,
 * which pointer must have previously been obtained by a call
 * to nlspec_new().
 *
 * Returns: 0 on success, non-zero error code on error.
 */

int nlspec_add_aux(Nlspec * spec, const(char) * s, const(DatasetRaw) * dset)
{
    int ci, err;

    err = screen_bad_aux(s, dset, &ci);
    if (!err) {
    if (ci == GretlCmdIndex.printf) {
        char * tmp = revise_aux_printf(s);

        if (tmp != null) {
        err = strings_array_add(&spec.aux, &spec.naux, tmp);
        g_free(tmp);
        }
    } else {
        err = strings_array_add(&spec.aux, &spec.naux, s);
    }
    }

    return err;
}

/**
 * nlspec_set_regression_function:
 * @spec: pointer to nls specification.
 * @fnstr: string specifying nonlinear regression function.
 * @dset: information on dataset.
 *
 * Adds the regression function to @spec.  This pointer must
 * have previously been obtained by a call to nlspec_new().
 * The required format for @fnstr is "%varname = %formula", where
 * %varname is the name of the dependent variable and %formula
 * is an expression of the sort that is fed to gretl's %genr command.
 * The dependent variable must be already present in the
 * dataset.
 *
 * Returns: 0 on success, non-zero error code on error.
 */

int nlspec_set_regression_function(Nlspec * spec, const(char) * fnstr, const(DatasetRaw) * dset)
{
    const(char) * p = fnstr;
    char * vname = null;
    char * rhs = null;
    int flen, err = 0;

    if (spec.nlfunc != null) {
    free(spec.nlfunc);
    spec.nlfunc = null;
    }

    spec.dv = 0;

    if (!strncmp(p, "nls ", 4) ||
    !strncmp(p, "mle ", 4) ||
    !strncmp(p, "gmm ", 4)) {
    p += 4;
    } else if (!strncmp(p, "gmm", 3)) {
    p += 3;
    }

    if (spec.ci == GretlCmdIndex.gmm && string_is_blank(p)) {
    /* GMM: we don't insist on a function on the first line */
    return 0;
    }

    if (equation_get_lhs_and_rhs(p, &vname, &rhs)) {
    gretl_errmsg_sprintf(_("parse error in '%s'\n").ptr, fnstr);
    err =  GretlError.parse;
    } else if (spec.ci == GretlCmdIndex.nls) {
    /* the dependent variable must be a series */
    spec.dv = series_index(dset, vname);
    if (spec.dv == dset.v) {
        gretl_errmsg_sprintf(_("'%s' is not a known series").ptr, vname);
        err = GretlError.unkvar;
    }
    } else {
    *spec.lhname.ptr = '\0';
    strncat(spec.lhname.ptr, vname, VNAMELEN - 1);
    }

    if (!err) {
    if (spec.ci == GretlCmdIndex.mle || spec.ci == GretlCmdIndex.gmm) {
        spec.nlfunc = gretl_strdup(rhs);
    } else {
        flen = cast(int) (strlen(vname) + strlen(rhs) + 6);
        spec.nlfunc = cast(typeof(spec.nlfunc)) malloc(flen);
        if (spec.nlfunc != null) {
        /* the equation defining the NLS residual */
        sprintf(spec.nlfunc, "%s - (%s)".ptr, vname, rhs);
        }
    }

    if (spec.nlfunc == null) {
        err = GretlError.alloc;
    }
    }

    free(vname);
    free(rhs);

    return err;
}

/**
 * nlspec_set_t1_t2:
 * @spec: pointer to nls specification.
 * @t1: starting observation.
 * @t2: ending observation.
 *
 * Sets the sample range for estimation of @spec.  This pointer must
 * have previously been obtained by a call to nlspec_new().
 */

void nlspec_set_t1_t2(Nlspec * spec, int t1, int t2)
{
    if (spec != null) {
    spec.t1 = t1;
    spec.t2 = t2;
    spec.nobs = t2 - t1 + 1;
    }
}

bool parnames_line(const(char) * s)
{
    return (!strncmp(s, "param_names ", 12));
}

bool param_line(const(char) * s)
{
    return (!strncmp(s, "deriv ", 6) || !strncmp(s, "params ", 7) || !strncmp(s, "orthog ", 7) || !strncmp(s, "weights ", 8) || !strncmp(s, "hessian ", 8));
}

bool cmd_start(const(char) * s)
{
    return (!strncmp(s, "nls ", 4) || !strncmp(s, "mle ", 4) || !strncmp(s, "gmm ", 4) || !strcmp(s, "gmm"));
}

/**
 * nl_parse_line:
 * @ci: %NLS, %MLE or %GMM.
 * @line: string containing information to be added to a
 * nonlinear model specification.
 * @dset: dataset struct.
 * @prn: gretl printing struct (for warning messages).
 *
 * This function is used to create the specification of a
 * nonlinear model, to be estimated via #nl_model (i.e.
 * via NLS, MLE or GMM).  It can be called several times to
 * build up the required information.  See the manual
 * entries for nls, mle and gmm for details.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int nl_parse_line(int ci, const(char) * line, const(DatasetRaw) * dset, Prn * prn)
{
    Nlspec * s = &private_spec;
    int err = 0;

    s.ci = ci;

    if (cmd_start(line)) {
    if (s.nlfunc != null) {
        clear_nlspec(s);
    }
    err = nlspec_set_regression_function(s, line, dset);
    if (!err) {
        nlspec_set_t1_t2(s, dset.t1, dset.t2);
    }
    } else if (param_line(line)) {
    if (s.nlfunc == null && s.ci != GretlCmdIndex.gmm) {
        gretl_errmsg_set(_("No regression function has been specified").ptr);
        err = GretlError.parse;
    } else {
        if (!strncmp(line, "deriv", 5)) {
        if (numeric_mode(s) && s.params != null) {
            gretl_errmsg_set(_("You cannot supply both a \"params\" line and analytical derivatives").ptr);
            err = GretlError.parse;
        } else {
            err = nlspec_add_param_with_deriv(s, line);
        }
        } else if (!strncmp(line, "params", 6)) {
        if (numeric_mode(s)) {
            err = nlspec_add_params_from_line(s, line + 6);
        } else {
            pprintf(prn, _("Analytical derivatives supplied: \"params\" line will be ignored").ptr);
            pputc(prn, '\n');
        }
        } else if (!strncmp(line, "orthog", 6)) {
        err = nlspec_add_orthcond(s, line + 6, dset);
        } else if (!strncmp(line, "weights", 7)) {
        err = nlspec_add_weights(s, line + 7);
        } else if (!strncmp(line, "hessian", 7)) {
        err = nlspec_add_hessian(s, line + 8);
        }
    }
    } else if (parnames_line(line)) {
    err = nlspec_add_param_names(s, line + 12);
    } else {
    err = nlspec_add_aux(s, line, dset);
    }

    if (err) {
    /* remember to clean up! */
    clear_nlspec(s);
    }

    return err;
}

int nl_set_smallstep()
{
    if (private_spec.nlfunc != null) {
    private_spec.flags |= NlFlags.smallstep;
    return 0;
    } else {
    return GretlError.data;
    }
}

__gshared double default_nls_toler = 0;

/**
 * get_default_nls_toler:
 *
 * Returns: the default value used in the convergence criterion
 * for estimation of models using nonlinear least squares.
 */

double get_default_nls_toler()
{
    if (default_nls_toler == 0.0) {
    default_nls_toler = pow(DBL_EPSILON, .75);
    }

    return default_nls_toler;
}

int check_spec_requirements(Nlspec * spec)
{
    int err = 0;

    if (spec.nparam < 1) {
    gretl_errmsg_set(_("No parameters have been specified").ptr);
    err = 1;
    } else if (spec.ci == GretlCmdIndex.gmm) {
    err = check_gmm_requirements(spec);
    }

    return err;
}

/* make any adjustments that may be needed in case the value
   returned by the user's likelihood function is a scalar
*/

int mle_scalar_check(Nlspec * spec)
{
    int err = 0;

    if (scalar_loglik(spec)) {
    /* without per-observation likelihood values, we can do
       neither OPG nor QML */
    if (spec.opt & Gretlopt.OPT_R) {
        gretl_errmsg_set("Scalar loglikelihood: can't do QML");
        err = GretlError.badopt;
    } else if (!(spec.opt & Gretlopt.OPT_A)) {
        /* ensure that we use the Hessian */
        spec.opt |= Gretlopt.OPT_H;
    }
    if (analytic_mode(spec)) {
        /* analytic mode: suppress gradient check */
        spec.opt |= Gretlopt.OPT_S;
    }
    }

    return err;
}

void nls_run_GNR(Model * pmod, Nlspec * spec, Prn * prn)
{
    DatasetRaw * gdset = null;
    int * glist = null;
    int perfect = 0;
    int err = 0;

    gdset = make_GNR_dataset(spec, spec.dset, &glist, &perfect,
                 prn, &err);

    if (err) {
    pmod.errcode = err;
    } else {
    *pmod = GNR(glist, gdset, spec.opt, prn);
    if (pmod.errcode == 0 || pmod.errcode == GretlError.jacobian) {
        if (pmod.errcode == 0 && (spec.opt & Gretlopt.OPT_B)) {
        QLR_test(pmod, gdset, Gretlopt.OPT_Q | Gretlopt.OPT_M, null);
        }
        pmod.errcode = finalize_nls_model(pmod, spec,
                           perfect, glist);
    }
    destroy_dataset(gdset);
    free(glist);
    }
}

int check_lbfgs_option(Nlspec * spec)
{
    const(char) * s = get_optval_string(GretlCmdIndex.gmm, Gretlopt.OPT_L);
    int err = 0;

    if (s != null) {
        GretlMatrixRaw * bounds = get_matrix_by_name(s);

        if (bounds == null) {
            gretl_errmsg_sprintf(_("Unknown matrix '%s'").ptr, s);
            err = GretlError.unkvar;
        } else {
            spec.bounds = bounds;
        }
    }

    return err;
}

/* static function providing the real content for the two public
   wrapper functions below: does NLS, MLE or GMM */

Model real_nl_model(Nlspec * spec, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    Model nlmod;
    int origv = 0;
    int err = 0;

    if (spec == null) {
    /* we use the static spec composed via nl_parse_line() */
    spec = &private_spec;
    }

    gretl_model_init(&nlmod, dset);

    if (dset != null) {
    origv = dset.v;
    }

    if (spec.ci == GretlCmdIndex.gmm) {
        /* can't have both --iterate and --two-step */
    nlmod.errcode = incompatible_options(opt, Gretlopt.OPT_I | Gretlopt.OPT_T);
        if (!nlmod.errcode && (opt & Gretlopt.OPT_L)) {
            nlmod.errcode = check_lbfgs_option(spec);
        }
    } else if (spec.nlfunc == null) {
    gretl_errmsg_set(_("No function has been specified").ptr);
    nlmod.errcode = GretlError.parse;
    }

    if (nlmod.errcode) {
    goto bailout;
    }

    spec.opt = opt;

    /* make @dset and @prn available via @spec */
    spec.dset = dset;
    spec.prn = prn;

    if (spec.opt & Gretlopt.OPT_N) {
    /* ignore any analytical derivatives */
    set_numeric_mode(spec);
    }

    if (numeric_mode(spec) && spec.nparam == 0 && spec.ci != GretlCmdIndex.gmm) {
    gretl_errmsg_sprintf(_("%s: no parameters were specified").ptr,
                 gretl_command_word(spec.ci));
    nlmod.errcode = GretlError.data;
    goto bailout;
    }

    if (check_spec_requirements(spec)) {
    nlmod.errcode = GretlError.parse;
    goto bailout;
    }

    if (nl_coeff_check(spec)) {
    nlmod.errcode = GretlError.data;
    goto bailout;
    }

    if (spec.ci == GretlCmdIndex.gmm) {
    err = gmm_missval_check_etc(spec);
    } else if (dset != null) {
    err = nl_missval_check(spec, dset);
    }

    if (!err && spec.ci == GretlCmdIndex.mle) {
    err = mle_scalar_check(spec);
    }

    if (err) {
    nlmod.errcode = err;
    goto bailout;
    }

    if (spec.lhtype == GretlType.double_) {
    spec.fvec = null;
    spec.J = null;
    } else {
    /* allocate auxiliary arrays */
    size_t fvec_bytes = spec.nobs * (*spec.fvec).sizeof;

    spec.fvec = cast(typeof(spec.fvec)) malloc(fvec_bytes);
    spec.J = gretl_matrix_alloc(spec.nobs, spec.ncoeff);

    if (spec.fvec == null || spec.J == null) {
        nlmod.errcode = GretlError.alloc;
        goto bailout;
    }

    if (spec.lvec != null) {
        memcpy(spec.fvec, spec.lvec.val, fvec_bytes);
    } else {
        /* not scalar or matrix: must be a series */
        if (dset == null || dset.v == 0) {
        nlmod.errcode = GretlError.nodata;
        goto bailout;
        } else {
        double * src = spec.dset.Z[spec.lhv] + spec.t1;

        memcpy(spec.fvec, src, fvec_bytes);
        }
    }
    }

    /* get tolerance from user setting or default */
    if (spec.ci == GretlCmdIndex.mle || spec.ci == GretlCmdIndex.gmm) {
    spec.tol = libset_get_double(SetKey.bfgs_toler);
    } else {
    spec.tol = libset_get_double(SetKey.nls_toler);
    }

    if (spec.ci != GretlCmdIndex.gmm && !(spec.opt & (Gretlopt.OPT_Q | Gretlopt.OPT_M))) {
    pputs(prn, (numeric_mode(spec))?
          _("Using numerical derivatives\n").ptr :
          _("Using analytical derivatives\n").ptr);
    }

    /* now start the actual calculations */

    if (spec.ci == GretlCmdIndex.mle) {
    err = mle_calculate(spec, prn);
    } else if (spec.ci == GretlCmdIndex.gmm) {
    err = gmm_calculate(spec, prn);
    if (err) {
        fprintf(stderr, "gmm_calculate returned %d\n".ptr, err);
    }
    } else {
    /* NLS: invoke the appropriate minpack driver function */
    gretl_iteration_push();
    if (numeric_mode(spec)) {
        err = lm_approximate(spec, prn);
    } else {
        err = lm_calculate(spec, prn);
        if (err) {
        fprintf(stderr, "lm_calculate returned %d\n".ptr, err);
        }
    }
    gretl_iteration_pop();
    }

    if (!(spec.opt & (Gretlopt.OPT_Q | Gretlopt.OPT_M)) && !(spec.flags & NlFlags.newton)) {
    pprintf(prn, _("Tolerance = %g\n").ptr, spec.tol);
    }

    if (!err) {
    if (spec.ci == GretlCmdIndex.nls) {
        if (spec.opt & Gretlopt.OPT_C) {
        /* coefficients only: don't bother with GNR */
        add_nls_coeffs(&nlmod, spec);
        } else {
        /* use Gauss-Newton Regression for covariance matrix,
           standard errors */
        nls_run_GNR(&nlmod, spec, prn);
        }
    } else {
        /* MLE, GMM */
        make_other_nl_model(&nlmod, spec, dset);
    }
    } else if (nlmod.errcode == 0) {
    /* model error code missing */
    nlmod.errcode = err;
    }

    /* ensure that the canonical parameter values get back
       into external scalars or matrices */
    update_coeff_values(spec.coeff, spec);

 bailout:

    if (spec.nmat > 0 && analytic_mode(spec)) {
    destroy_private_matrices();
    }

    destroy_private_scalars();
    clear_nlspec(spec);
    destroy_private_lists();

    if (dset != null && !dset.auxiliary) {
    dataset_drop_last_variables(dset, dset.v - origv);
    }

    if (!nlmod.errcode && !(opt & Gretlopt.OPT_A)) {
    set_model_id(&nlmod, opt);
    }

    return nlmod;
}

/**
 * nl_model:
 * @dset: dataset struct.
 * @opt: may include %OPT_V for verbose output, %OPT_R
 * for robust covariance matrix.
 * @prn: printing struct.
 *
 * Computes estimates of a model via nonlinear least squares,
 * maximum likelihood, or GMM.  The model must have been specified
 * previously, via calls to the function #nl_parse_line.  Those
 * calls determine, among other things, which estimator will
 * be used.
 *
 * Returns: a model struct containing the parameter estimates
 * and associated statistics.
 */

Model nl_model(DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    return real_nl_model(null, dset, opt, prn);
}

/**
 * model_from_nlspec:
 * @spec: nonlinear model specification.
 * @dset: dataset struct.
 * @opt: may include %OPT_V for verbose output, %OPT_A to
 * treat as an auxiliary model, %OPT_C to produce coefficient
 * estimates only (don't bother with GNR to produce standard
 * errors).
 * @prn: printing struct.
 *
 * Computes estimates of the model specified in @spec.
 * The @spec must first be obtained using nlspec_new(), and
 * initialized using nlspec_set_regression_function().  If analytical
 * derivatives are to be used (which is optional but recommended)
 * these are set using nlspec_add_param_with_deriv().
 *
 * Returns: a model struct containing the parameter estimates
 * and associated statistics.
 */

Model model_from_nlspec(Nlspec * spec, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    return real_nl_model(spec, dset, opt, prn);
}

/**
 * nlspec_new:
 * @ci: %NLS, %MLE or %GMM.
 * @dset: information on dataset.
 *
 * Returns: a pointer to a newly allocated nonlinear model
 * specification, or NULL on failure.
 */

Nlspec * nlspec_new(int ci, const(DatasetRaw) * dset)
{
    Nlspec * spec;

    spec = cast(typeof(spec)) malloc((*spec).sizeof);
    if (spec == null) {
    return null;
    }

    spec.nlfunc = null;

    spec.params = null;
    spec.nparam = 0;

    spec.aux = null;
    spec.naux = 0;

    spec.genrs = null;
    spec.ngenrs = 0;
    spec.generr = 0;

    spec.hgen = null;
    spec.hesscall = null;

    spec.fvec = null;
    spec.J = null;

    spec.coeff = null;
    spec.ncoeff = 0;
    spec.nmat = 0;

    spec.Hinv = null;
    spec.bounds = null;

    spec.ci = ci;
    spec.flags = cast(NlFlags) 0;
    spec.opt = Gretlopt.OPT_NONE;

    spec.dv = 0;
    spec.lhtype = GretlType.none;
    *spec.lhname.ptr = '\0';
    *spec.hname.ptr = '\0';
    spec.parnames = null;
    spec.lhv = 0;
    spec.lvec = null;

    spec.iters = 0;
    spec.fncount = 0;
    spec.grcount = 0;

    spec.t1 = spec.real_t1 = dset.t1;
    spec.t2 = spec.real_t2 = dset.t2;
    spec.nobs = spec.t2 - spec.t1 + 1;

    spec.dset = null;
    spec.prn = null;

    spec.oc = null;
    spec.missmask = null;

    return spec;
}

/**
 * nlspec_destroy:
 * @spec: pointer to nls specification.
 *
 * Frees all resources associated with @spec, and frees the
 * pointer itself.
 */

void nlspec_destroy(Nlspec * spec)
{
    clear_nlspec(spec);
    free(spec);
}

/* below: apparatus for estimating a "tsls" model via GMM */

enum IVREG_RESIDNAME = "_gmm_e";
enum IVREG_WEIGHTNAME = "_gmm_V";

int ivreg_nlfunc_setup(Nlspec * spec, Model * pmod, DatasetRaw * dset)
{
    Prn * prn;
    int v = pmod.list[1];
    int k = pmod.ncoeff;
    const(char) * buf;
    int i, err = 0;

    prn = gretl_print_new(PrnType.buffer, &err);
    if (err) {
    return err;
    }

    pprintf(prn, "gmm %s = %s-(".ptr, IVREG_RESIDNAME.ptr, dset.varname[v]);
    for (i=0; i<k; i++) {
    if (i == 0 && pmod.ifc) {
        pputs(prn, "_b0");
    } else {
        v = pmod.list[i+2];
        pprintf(prn, "+_b%d*%s".ptr, i, dset.varname[v]);
    }
    }
    pputc(prn, ')');

    buf = gretl_print_get_buffer(prn);
    err = nlspec_set_regression_function(spec, buf, dset);
    gretl_print_destroy(prn);

    return err;
}

int ivreg_oc_setup(Nlspec * spec, const(int) * ilist, Model * pmod, DatasetRaw * dset, int * rv)
{
    int v = dset.v;
    int i, err;

    /* add GMM residual */
    err = dataset_add_series(dset, 1);

    if (!err) {
    *rv = v;
    strcpy(dset.varname[v], IVREG_RESIDNAME.ptr);
    for (i=0; i<dset.n; i++) {
        dset.Z[v][i] = pmod.uhat[i];
    }
    err = nlspec_add_ivreg_oc(spec, v, ilist, cast(const(double) **) dset.Z);
    }

    return err;
}

int ivreg_weights_setup(Nlspec * spec, const(int) * ilist, Gretlopt opt)
{
    const(char) * mname = null;
    int k = ilist[0];
    GretlMatrixRaw * V;
    int err;

    if (opt & Gretlopt.OPT_H) {
    mname = get_optval_string(GretlCmdIndex.ivreg, Gretlopt.OPT_H);
    }

    if (mname != null) {
    /* the user requested a specific weights matrix */
    V = get_matrix_by_name(mname);
    if (V == null) {
        err = GretlError.unkvar;
    } else if (V.rows != k || V.cols != k) {
        err = GretlError.nonconf;
    } else {
        err = nlspec_add_weights(spec, mname);
    }
    } else {
    /* use generic weights */
    V = gretl_identity_matrix_new(k);
    if (V == null) {
        return GretlError.alloc;
    }
    err = private_matrix_add(V, IVREG_WEIGHTNAME.ptr);
    if (err) {
        gretl_matrix_free(V);
    } else {
        err = nlspec_add_weights(spec, IVREG_WEIGHTNAME.ptr);
    }
    }

    return err;
}

int ivreg_set_params(Nlspec * spec, Model * pmod)
{
    int np = pmod.ncoeff;
    char ** names;
    int i, err;

    names = strings_array_new_with_length(np, 8);
    if (names == null) {
    return GretlError.alloc;
    }

    for (i=0; i<np; i++) {
    sprintf(names[i], "_b%d".ptr, i);
    }

    err = aux_nlspec_add_param_list(spec, np, pmod.coeff, names);

    strings_array_free(names, np);

    return err;
}

int finalize_ivreg_model(Model * pmod, Model * ols, const(int) * biglist, const(int) * mlist, const(int) * ilist, const(double) ** Z, int rv)
{
    int * endolist;
    int yno = mlist[1];
    int t, err = 0;

    pmod.ci = GretlCmdIndex.ivreg;
    pmod.opt = Gretlopt.OPT_G; /* signal use of GMM */

    /* attach tsls-style list */
    pmod.list = gretl_list_copy(biglist);
    if (pmod.list == null) {
    return GretlError.alloc;
    }

    /* get dep. var. info from OLS model */
    pmod.ybar = ols.ybar;
    pmod.sdy = ols.sdy;

    /* steal uhat */
    if (pmod.uhat != null) {
    free(pmod.uhat);
    }
    pmod.uhat = ols.uhat;
    ols.uhat = null;

    /* and yhat */
    if (pmod.yhat != null) {
    free(pmod.yhat);
    }
    pmod.yhat = ols.yhat;
    ols.yhat = null;

    /* insert residuals and fitted */
    for (t=pmod.t1; t<=pmod.t2; t++) {
    pmod.uhat[t] = Z[rv][t];
    pmod.yhat[t] = Z[yno][t] - Z[rv][t];
    }

    free(pmod.xlist);
    pmod.xlist = gretl_list_sublist(mlist, 2, -1);

    endolist = tsls_make_endolist(mlist, ilist, &err);

    if (endolist != null) {
    gretl_model_set_list_as_data(pmod, "endolist", endolist);
    }

    return err;
}

/* Responds when OPT_L is given to the ivreg() function,
   which lives in estimate.c
*/

Model ivreg_via_gmm(const(int) * list, DatasetRaw * dset, Gretlopt opt)
{
    int orig_v = dset.v;
    int * mlist = null, ilist = null;
    Nlspec * spec = null;
    Model olsmod, model;
    int rv = 0;
    int err = 0;

    gretl_model_init(&olsmod, dset);
    gretl_model_init(&model, dset);

    err = ivreg_process_lists(list, &mlist, &ilist);

    if (!err) {
    spec = nlspec_new(GretlCmdIndex.gmm, dset);
    if (spec == null) {
        err = GretlError.alloc;
    }
    }

    if (!err) {
    /* OLS baseline and check */
    olsmod = lsq(mlist, dset, GretlCmdIndex.ols, Gretlopt.OPT_A | Gretlopt.OPT_M | Gretlopt.OPT_Z);
    err = olsmod.errcode;
    }

    if (!err) {
    /* add the definition of the GMM residual */
    err = ivreg_nlfunc_setup(spec, &olsmod, dset);
    }

    if (!err) {
    /* add the orthogonality conditions */
    err = ivreg_oc_setup(spec, ilist, &olsmod, dset, &rv);
    }

    if (!err) {
    /* add the weights matrix */
    err = ivreg_weights_setup(spec, ilist, opt);
    }

    set_auxiliary_scalars();

    if (!err) {
    /* add the scalar parameters */
    err = ivreg_set_params(spec, &olsmod);
    }

    if (!err) {
    /* call the GMM estimation routine */
    model = model_from_nlspec(spec, dset, opt | Gretlopt.OPT_G, null);
    }

    unset_auxiliary_scalars();

    if (err) {
    model.errcode = err;
    } else {
    /* turn the output model into an "ivreg" type */
    finalize_ivreg_model(&model, &olsmod, list, mlist, ilist,
                 cast(const(double) **) dset.Z, rv);
    }

    nlspec_destroy(spec);
    free(mlist);
    free(ilist);
    clear_model(&olsmod);

    dataset_drop_last_variables(dset, dset.v - orig_v);
    user_var_delete_by_name(IVREG_WEIGHTNAME.ptr, null);

    return model;
}

/* For forecasting purposes ("fcast" command): in general
   it's not enough just to run the genr command that
   creates or revises the dependent variable; we may also
   have to generate some intermediate quantities. In
   nl_model_run_aux_genrs() we retrieve any auxiliary
   genrs from @pmod and try running them.
*/

int nl_model_run_aux_genrs(const(Model) * pmod, DatasetRaw * dset)
{
    char[MAXLEN] line;
    const(char) * buf;
    int i, j, n_aux = 0;
    int err = 0;

    n_aux = gretl_model_get_int(pmod, "nl_naux");
    if (n_aux == 0) {
    /* nothing to be done? */
    return 0;
    }

    buf = cast(typeof(buf)) gretl_model_get_data(pmod, "nlinfo");
    if (buf == null) {
    return GretlError.data;
    }

    bufgets_init(buf);

    i = j = 0;
    while (bufgets(line.ptr, line.sizeof, buf) && !err) {
    if (*line.ptr == '#') {
        continue;
    }
    if (i > 0 && j < n_aux) {
        tailstrip(line.ptr);
        err = generate(line.ptr, dset, GretlType.any, Gretlopt.OPT_P, null);
        if (err) {
        genr_setup_error(line.ptr, err, 1);
        }
        j++;
    }
    i++;
    }

    bufgets_finalize(buf);

    return err;
}

/* apparatus for bootstrapping NLS forecast errors */

/* reconstitute the nlspec based on the information saved in @pmod */

int set_nlspec_from_model(const(Model) * pmod, const(DatasetRaw) * dset)
{
    char[MAXLEN] line;
    const(char) * buf;
    int err = 0;

    buf = cast(typeof(buf)) gretl_model_get_data(pmod, "nlinfo");
    if (buf == null) {
    return GretlError.data;
    }

    bufgets_init(buf);

    while (bufgets(line.ptr, line.sizeof, buf) && !err) {
    tailstrip(line.ptr);
    if (!strcmp(line.ptr, "end nls")) {
        break;
    }
    err = nl_parse_line(GretlCmdIndex.nls, line.ptr, dset, null);
    }

    bufgets_finalize(buf);

    if (err) {
    clear_nlspec(&private_spec);
    } else {
    private_spec.opt = pmod.opt;
    if (private_spec.opt & Gretlopt.OPT_V) {
        private_spec.opt ^= Gretlopt.OPT_V;
    }
    }

    return err;
}

int save_scalar_param(Parm * p, double ** px, int * n)
{
    double x = get_param_scalar(p);
    int err = 0;

    if (na(x)) {
    err = GretlError.data;
    } else {
    double * tmp;
    int k = *n + 1;

    tmp = cast(typeof(tmp)) realloc(*px, k * (*tmp).sizeof);
    if (tmp == null) {
        err = GretlError.alloc;
    } else {
        *px = tmp;
        tmp[*n] = x;
        *n = k;
    }
    }

    return err;
}

int save_matrix_param(Parm * p, GretlMatrixRaw *** pm, int * n)
{
    GretlMatrixRaw * m = get_param_matrix(p);
    int err = 0;

    if (m == null) {
    err = GretlError.data;
    } else {
    GretlMatrixRaw ** tmp;
    int k = *n + 1;

    tmp = cast(typeof(tmp)) realloc(*pm, k * (*tmp).sizeof);
    if (tmp == null) {
        err = GretlError.alloc;
    } else {
        *pm = tmp;
        tmp[*n] = gretl_matrix_copy(m);
        if (tmp[*n] == null) {
        err = GretlError.alloc;
        } else {
        *n = k;
        }
    }
    }

    return err;
}

/* Given an NLS model @pmod, fill out the forecast error series
   @fcerr from @ft1 to @ft2 using the bootstrap, based on
   resampling the original residuals. See also forecast.c,
   from where this function is called.  Here we concentrate
   solely on parameter uncertainty.  At present we do this
   only for static forecasts.
*/

int nls_boot_calc(const(Model) * pmod, DatasetRaw * dset, int ft1, int ft2, double * fcerr)
{
    Nlspec * spec;
    GretlMatrixRaw * fcmat = null;
    GretlMatrixRaw ** msave = null;
    double * orig_y = null;
    double * xsave = null;
    double * resu = null;
    int origv = dset.v;
    int yno = pmod.list[1];
    int iters = 100; /* just testing */
    int nx = 0, nm = 0;
    int i, ix, im, s, t, fT;
    int err = 0;

    /* build the 'private' spec, based on @pmod */
    err = set_nlspec_from_model(pmod, dset);
    if (err) {
    return err;
    }

    spec = &private_spec;

    spec.t1 = spec.real_t1 = pmod.t1;
    spec.t2 = spec.real_t2 = pmod.t2;
    spec.nobs = spec.t2 - spec.t1 + 1;

    /* back up the original parameter values */
    for (i=0; i<spec.nparam && !err; i++) {
    Parm * p = &spec.params[i];

    if (p.type == GretlType.double_) {
        err = save_scalar_param(p, &xsave, &nx);
    } else {
        err = save_matrix_param(p, &msave, &nm);
    }
    }

    if (err) {
    goto bailout;
    }

    /* back up the original dependent variable */
    orig_y = copyvec(dset.Z[yno], dset.n);

    /* allocate storage for resampled residuals */
    resu = cast(typeof(resu)) malloc(dset.n * (*resu).sizeof);

    /* allocate forecast matrix */
    fT = ft2 - ft1 + 1;
    fcmat = gretl_matrix_alloc(iters, fT);

    if (orig_y == null || resu == null || fcmat == null) {
    err = GretlError.alloc;
    goto bailout;
    }

    /* allocate arrays to be passed to minpack */
    spec.fvec = cast(typeof(spec.fvec)) malloc(spec.nobs * (*spec.fvec).sizeof);
    spec.J = gretl_matrix_alloc(spec.nobs, spec.ncoeff);

    if (spec.fvec == null || spec.J == null) {
    err = GretlError.alloc;
    goto bailout;
    }

    spec.dset = dset;

    for (i=0; i<iters && !err; i++) {

    /* resample from pmod->uhat into resu */
    resample_series(pmod.uhat, resu, dset);

    /* construct y^* = \hat{y} + u^* */
    for (t=spec.t1; t<=spec.t2; t++) {
        /* FIXME rescale the residuals */
        dset.Z[yno][t] = pmod.yhat[t] + resu[t];
    }

    for (t=spec.t1, s=0; t<=spec.t2; t++, s++) {
        spec.fvec[s] = resu[s];
    }

    /* re-estimate the model on the resampled data */
    if (numeric_mode(spec)) {
        err = lm_approximate(spec, null);
    } else {
        err = lm_calculate(spec, null);
    }

    if (!err) {
        /* generate and record the forecast (FIXME genr expr?) */
        err = generate(pmod.depvar, dset, GretlType.series,
               Gretlopt.OPT_P, null);
        for (t=ft1, s=0; t<=ft2; t++, s++) {
        gretl_matrix_set(fcmat, i, s, dset.Z[yno][t]);
        }
    }
    }

    if (!err) {
    /* compute and record the standard deviations of the forecasts:
       should we add the square root of the residual variance here?
    */
    GretlMatrixRaw * sd;

    sd = gretl_matrix_column_sd(fcmat, 0, 0, &err);
    if (!err) {
        double cfac = sqrt(cast(double) iters / (iters - 1));

        for (t=ft1, s=0; t<=ft2; t++, s++) {
        fcerr[t] = sd.val[s] * cfac;
        }
        gretl_matrix_free(sd);
    }
    }

    if (spec.nmat > 0 && analytic_mode(spec)) {
    destroy_private_matrices();
    }

 bailout:

    /* restore original params from backup */
    ix = im = 0;
    for (i=0; i<spec.nparam; i++) {
    Parm * p = &spec.params[i];

    if (p.type == GretlType.double_ && ix < nx) {
        if (p.bundle != null) {
        gretl_bundle_set_scalar(p.bundle, p.name.ptr, xsave[ix]);
        } else {
        gretl_scalar_set_value(p.name.ptr, xsave[ix]);
        }
        ix++;
    } else if (im < nm) {
        GretlMatrixRaw * m0 = get_param_matrix(p);
        GretlMatrixRaw * m1 = msave[im];
        int k, len = m0.rows * m0.cols;

        for (k=0; k<len; k++) {
        m0.val[k] = m1.val[k];
        }
        gretl_matrix_free(msave[im]);
        im++;
    }
    }

    free(xsave);
    free(msave);

    clear_nlspec(spec);

    dataset_drop_last_variables(dset, dset.v - origv);

    /* restore original y */
    if (orig_y != null) {
    for (t=0; t<dset.n; t++) {
        dset.Z[yno][t] = orig_y[t];
    }
    }

    free(orig_y);
    free(resu);
    gretl_matrix_free(fcmat);

    return err;
}

