/* Port of lib/src/gretl_restrict.c: sets of linear restrictions on the coefficients of models, systems and VECMs
   (parsing of "restrict" blocks, the R and q matrices, Wald/LR/F tests of the restrictions, restricted VECM
   estimation and nonlinear restrictions through a function). The struct behind the opaque GretlRestriction of
   gretl_restrict.h is declared in gretl_restricth.d. */

module gretl.gretl_restrict;


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
import gretl.varh;

import gretl.xmlh;


private alias XUC = const(xmlChar) *;


import gretl.gretl_restricth : GretlRestriction, RRow;

import gretl.bootstraph : BootMethod;
import core.stdc.config : c_long;

/* gretl_func.h is not ported yet */
enum FN_NAMELEN = 32;

/* the C library's strtol and strtod take a plain char ** end pointer */
extern (C) {
    c_long strtol(const(char) * nptr, char ** endptr, int base);
    double strtod(const(char) * nptr, char ** endptr);
}
import gretl.cmd_privateh : ExecState;
import gretl.geneval0 : Ufunc;
import gretl.genmain : generate_matrix, generate_scalar;
import gretl.genmainh : ModelDataIndex;
import gretl.gretl_bfgs : user_fdjac;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_errors : gretl_errmsg_is_set, gretl_errmsg_set, gretl_errmsg_sprintf, gretl_error_clear;
import gretl.gretl_list : gretl_list_copy, gretl_list_sublist;
import gretl.gretl_listh : LISTSEP;
import gretl.gretl_model : command_ok_for_model, exact_fit_check, gretl_model_allocate_storage, gretl_model_free, gretl_model_get_depvar, gretl_model_get_param_name, gretl_model_get_param_number, gretl_model_new, gretl_model_set_int, gretl_model_smpl_init, gretl_model_write_vcv, set_model_id;
import gretl.gretl_modelh : ASYMPTOTIC_MODEL, NONLIST_MODEL;
import gretl.gretl_prn : Prn, pprintf, pputc, pputs;
import gretl.gretl_utils : ls_criteria, positive_int_from_string, record_LR_test_result, record_test_result;
import gretl.matrix_extra : gretl_coeff_vector_from_model, gretl_matrix_print_to_prn, gretl_vcv_matrix_from_model;
import gretl.modelprint : print_coeffs;
import gretl.objstackh : GretlObjType;
import gretl.options : get_optval_string;
import gretl.printout : get_gretl_digits;
import gretl.pvalues : chisq_cdf_comp, normal_pvalue_2, snedecor_cdf_comp;
import gretl.strutils : gretl_charpos, gretl_namechar_spn, gretl_scan_varname, gretl_strdup, string_is_blank, strings_array_free, strings_array_new_with_length, tailstrip;
import gretl.systemh : EquationSystem, SystemFlags;
import gretl.var : alpha_restricted_VECM, beta_restricted_VECM, gretl_VAR_free, gretl_VAR_test, gretl_VECM_R_matrix, gretl_VECM_Ra_matrix, gretl_VECM_n_alpha, gretl_VECM_n_beta, gretl_VECM_q_matrix, gretl_VECM_qa_matrix, gretl_VECM_rank, gretl_VECM_test, real_gretl_restricted_vecm;
import gretl.system : equation_system_get_matrix, get_anonymous_equation_system, system_get_list, system_get_list_length, system_n_indep_vars, system_set_restriction_matrices, system_wald_test;
import gretl.gretl_array : get_strings_array_by_name, gretl_array_get_strings;
import gretl.uservar : destroy_private_matrices, gretl_is_scalar, gretl_scalar_get_value, private_matrix_add;
import gretl.usermat : get_matrix_by_name;

/* unported dependencies: per-file declarations from the 2026c headers */
extern (C) {
    int bootstrap_analysis(Model * pmod, int p, int B, double alpha, const(DatasetRaw) * dset, Gretlopt opt, Prn * prn);
    int bootstrap_ok(int ci);
    int bootstrap_test_restriction(Model * pmod, GretlMatrixRaw * R, GretlMatrixRaw * q, double test, int g, const(DatasetRaw) * dset, Gretlopt opt, int method, Prn * prn);
    void * get_last_model(GretlObjType * type);
    void * get_model_object_and_type(const(char) * name, GretlObjType * type);
    Ufunc * get_user_function_by_name(const(char) * name);
    void gretl_exec_state_set_model(ExecState * s, Model * pmod);
}
enum R_UNSPEC = -9;

enum VecmSide {
    none = 0,
    b = 1 << 0,
    a = 1 << 1
}

bool targ_specified(const(GretlRestriction) * r, int i, int j, char c)
{
    return (((c=='b' && r.bmulti) || (c=='a' && r.amulti)) && r.rows[i].eq != null && r.rows[i].eq[j] != R_UNSPEC);
}

/* @R is putatively the matrix in R \beta -q = 0.  Check that it makes
   sense in that role, i.e., that RR' is non-singular. */

int check_R_matrix(const(GretlMatrixRaw) * R)
{
    GretlMatrixRaw * m;
    int err = 0;

    m = gretl_matrix_alloc(R.rows, R.rows);

    if (m == null) {
    err = GretlError.alloc;
    } else {
    gretl_matrix_multiply_mod(R, GretlMatrixMod.none,
                  R, GretlMatrixMod.transpose,
                  m, GretlMatrixMod.none);

    err = gretl_invert_general_matrix(m);

    if (err == GretlError.singular) {
        gretl_errmsg_set(_("Matrix inversion failed: restrictions may be inconsistent or redundant").ptr);
    }

    gretl_matrix_free(m);
    }

    return err;
}

/* Here we're adding to an existing restriction on a VECM, either in
   respect of beta (@letter = 'b') or in respect of alpha
   (@letter = 'a').
*/

int augment_vecm_restriction(GretlRestriction * rset, char letter)
{
    GretlVar * vecm = cast(GretlVar *) rset.obj;
    const(GretlMatrixRaw) * R0, q0;
    GretlMatrixRaw * R2 = null, q2 = null;
    int err = 0;

    if (letter == 'b') {
    R0 = gretl_VECM_R_matrix(vecm);
    q0 = gretl_VECM_q_matrix(vecm);
    } else {
    R0 = gretl_VECM_Ra_matrix(vecm);
    q0 = gretl_VECM_qa_matrix(vecm);
    }

    R2 = gretl_matrix_row_concat(R0, rset.R, &err);
    if (err) {
    return err;
    }

    err = check_R_matrix(R2);
    if (err) {
    gretl_matrix_free(R2);
    return err;
    }

    if (q0 == null) {
    q2 = gretl_column_vector_alloc(R2.rows);
    if (q2 == null) {
        err = GretlError.alloc;
    } else {
        int i, n = R0.rows;

        for (i=0; i<R2.rows; i++) {
        q2.val[i] = (i < n)? 0.0 : rset.q.val[i-n];
        }
    }
    } else {
    q2 = gretl_matrix_row_concat(q0, rset.q, &err);
    }

    if (err) {
    gretl_matrix_free(R2);
    } else if (letter == 'b') {
    gretl_matrix_replace(&rset.R, R2);
    gretl_matrix_replace(&rset.q, q2);
    } else {
    gretl_matrix_replace(&rset.Ra, R2);
    gretl_matrix_replace(&rset.qa, q2);
    }

    return err;
}

/* Adding an existing restriction on a VECM, in respect of beta or
   alpha, to a new_ restriction in respect of the other matrix.
*/

int add_old_vecm_restriction(GretlRestriction * rset, char letter)
{
    GretlVar * vecm = cast(GretlVar *) rset.obj;
    int err = 0;

    if (letter == 'b') {
    const(GretlMatrixRaw) * R = gretl_VECM_R_matrix(vecm);
    const(GretlMatrixRaw) * q = gretl_VECM_q_matrix(vecm);

    rset.R = gretl_matrix_copy(R);
    if (rset.R == null) {
        err = GretlError.alloc;
    } else if (q != null) {
        rset.q = gretl_matrix_copy(q);
        if (rset.q == null) {
        err = GretlError.alloc;
        }
    }
    } else {
    const(GretlMatrixRaw) * Ra = gretl_VECM_Ra_matrix(vecm);
    const(GretlMatrixRaw) * qa = gretl_VECM_qa_matrix(vecm);

    rset.Ra = gretl_matrix_copy(Ra);
    if (rset.Ra == null) {
        err = GretlError.alloc;
    } else if (qa != null) {
        rset.qa = gretl_matrix_copy(qa);
        if (rset.qa == null) {
        err = GretlError.alloc;
        }
    }
    }

    return err;
}

/* We set things up here such that the restrictions are
   R_b * vec(beta) = q; R_a * vec(alpha') = 0.
*/

int get_R_vecm_column(const(GretlRestriction) * rset, int i, int j, char letter)
{
    const(RRow) * r = rset.rows[i];
    GretlVar * var = cast(GretlVar *) rset.obj;
    int col = 0;

    if (letter == 'b') {
    col = r.bnum[j];
    if (rset.bmulti) {
        col += r.eq[j] * gretl_VECM_n_beta(var);
    }
    } else if (letter == 'a') {
    col = r.bnum[j];
    if (rset.amulti) {
        col += r.eq[j] * gretl_VECM_n_alpha(var);
    }
    }

    return col;
}

int get_R_var_column(const(GretlRestriction) * rset, int i, int j)
{
    const(RRow) * r = rset.rows[i];
    GretlVar * var = cast(GretlVar *) rset.obj;
    int k, col = r.bnum[j];

    for (k=0; k<r.eq[j]; k++) {
    col += var.ncoeff;
    }

    return col;
}

int get_R_sys_column(const(GretlRestriction) * rset, int i, int j)
{
    const(RRow) * r = rset.rows[i];
    EquationSystem * sys = cast(EquationSystem *) rset.obj;
    int col = r.bnum[j];
    int k, n;

    for (k=0; k<r.eq[j]; k++) {
    n = system_get_list_length(sys, k);
    if (n > 0) {
        col += n - 1;
    }
    }

    return col;
}

double get_restriction_param(const(RRow) * r, int k)
{
    int i;

    for (i=0; i<r.nterms; i++) {
    if (r.bnum[i] == k) {
        return r.mult[i];
    }
    }

    return 0.0;
}

void vecm_cross_error()
{
    gretl_errmsg_set("VECM: beta/alpha cross restrictions are not handled");
}

/* See if we have both beta and alpha terms; check that there
   are no cross beta/alpha restrictions. */

int vecm_ab_check(GretlRestriction * rset)
{
    int atotal = 0, btotal = 0;
    int i, j, err = 0;

    for (i=0; i<rset.g && !err; i++) {
    RRow * r = rset.rows[i];
    int a = 0, b = 0;

    for (j=0; j<r.nterms && !err; j++) {
        if (r.letter[j] == 'a') {
        a++;
        atotal++;
        } else {
        b++;
        btotal++;
        }
        if (a > 0 && b > 0) {
        vecm_cross_error();
        err = GretlError.notimp;
        }
    }
    }

    if (!err) {
    /* incoming default is VecmSide.b */
    if (atotal > 0) {
        rset.vecm |= VecmSide.a;
    }
    if (btotal == 0) {
        rset.vecm &= ~VecmSide.b;
    }
    }

    return err;
}

/* check integrity of vecm restrictions, either beta or alpha */

int vecm_x_check(GretlRestriction * rset, char letter)
{
    int * multi = (letter == 'a')? &rset.amulti : &rset.bmulti;
    int * ki = (letter == 'a')? &rset.ga : &rset.gb;
    int unspec = 0, anyspec = 0;
    int i, j, err = 0;

    for (i=0; i<rset.g && !err; i++) {
    RRow * r = rset.rows[i];

    if (r.letter[0] != letter) {
        continue;
    }

    *ki += 1;

    for (j=0; j<r.nterms && !err; j++) {
        if (r.eq != null) {
        if (r.eq[j] == R_UNSPEC) {
            unspec = 1;
        } else if (r.eq[j] >= 0) {
            anyspec = 1;
        }
        }
        if (anyspec && unspec) {
        err = GretlError.parse;
        }
    }
    }

    if (!err && unspec && *multi) {
    for (i=0; i<rset.g; i++) {
        if (rset.rows[i].letter[0] == letter) {
        free(rset.rows[i].eq);
        rset.rows[i].eq = null;
        }
    }
    *multi = 0;
    }

    return err;
}

/* Check the validity of a set of VECM restrictions */

int vecm_restriction_check(GretlRestriction * rset)
{
    int r = gretl_VECM_rank(cast(GretlVar *) rset.obj);
    int err;

    err = vecm_ab_check(rset);

    if (!err && (rset.vecm & VecmSide.b)) {
    err = vecm_x_check(rset, 'b');
    }

    if (!err && (rset.vecm & VecmSide.a)) {
    err = vecm_x_check(rset, 'a');
    }

    if (!err && (rset.vecm & VecmSide.b)) {
    rset.bcols = gretl_VECM_n_beta(cast(GretlVar *) rset.obj);
    if (rset.bmulti) {
        rset.bcols *= r;
    }
    }

    if (!err && (rset.vecm & VecmSide.a)) {
    rset.acols = gretl_VECM_n_alpha(cast(GretlVar *) rset.obj);
    if (rset.amulti) {
        rset.acols *= r;
    }
    }

    if (!err) {
    if (rset.gb > rset.bcols) {
        gretl_errmsg_set("Too many restrictions");
        err = GretlError.nonconf;
    } else if (rset.ga > rset.acols) {
        fprintf(stderr, "rset->ga = %d, rset->acols = %d\n".ptr,
            rset.ga, rset.acols);
        gretl_errmsg_set("Too many restrictions");
        err = GretlError.nonconf;
    }
    }

    return err;
}

/* If we have been given restrictions in row-wise fashion,
   allocate the required R and q matrices */

int rset_allocate_R_q(GretlRestriction * rset, int nc)
{
    rset.R = gretl_zero_matrix_new(rset.g, nc);
    rset.q = gretl_zero_matrix_new(rset.g, 1);

    if (rset.R == null || rset.q == null) {
    gretl_matrix_free(rset.R);
    gretl_matrix_free(rset.q);
    rset.R = rset.q = null;
    return GretlError.alloc;
    }

    return 0;
}

/* For the single-equation case, allocate and fill out the
   R and q matrices. */

int equation_form_matrices(GretlRestriction * rset)
{
    Model * pmod = cast(Model *) rset.obj;
    int nc = 0;
    int col, i, j;

    if (rset.mask == null) {
    nc = pmod.ncoeff;
    } else {
    for (i=0; i<pmod.ncoeff; i++) {
        if (rset.mask[i]) {
        nc++;
        }
    }
    }

    if (rset_allocate_R_q(rset, nc)) {
    return GretlError.alloc;
    }

    for (i=0; i<rset.g; i++) {
    RRow * r = rset.rows[i];
    double x;

    col = 0;
    for (j=0; j<pmod.ncoeff; j++) {
        if (rset.mask[j]) {
        x = get_restriction_param(r, j);
        gretl_matrix_set(rset.R, i, col++, x);
        }
    }
    gretl_vector_set(rset.q, i, r.rhs);
    }

    return 0;
}

/* For the equation-system case, allocate and fill out the
   R and q matrices. */

int sys_form_matrices(GretlRestriction * rset)
{
    int nc = system_n_indep_vars(cast(EquationSystem *) rset.obj);
    int col, i, j;
    int err = 0;

    if (rset_allocate_R_q(rset, nc)) {
    return GretlError.alloc;
    }

    for (i=0; i<rset.g; i++) {
    RRow * r = rset.rows[i];

    for (j=0; j<r.nterms; j++) {
        col = get_R_sys_column(rset, i, j);
        gretl_matrix_set(rset.R, i, col, r.mult[j]);
    }
    gretl_vector_set(rset.q, i, r.rhs);
    }

    err = check_R_matrix(rset.R);

    return err;
}

/* For the VAR case, allocate and fill out the restriction
   matrices. */

int var_form_matrices(GretlRestriction * rset)
{
    GretlVar * var = cast(GretlVar *) rset.obj;
    int nc = var.neqns * var.ncoeff;
    int col, i, j;
    int err = 0;

    if (rset_allocate_R_q(rset, nc)) {
    return GretlError.alloc;
    }

    for (i=0; i<rset.g; i++) {
    RRow * r = rset.rows[i];

    for (j=0; j<r.nterms; j++) {
        col = get_R_var_column(rset, i, j);
        gretl_matrix_set(rset.R, i, col, r.mult[j]);
    }
    gretl_vector_set(rset.q, i, r.rhs);
    }

    err = check_R_matrix(rset.R);

    return err;
}

/* For the VECM case, allocate and fill out the restriction
   matrices. */

int vecm_form_matrices(GretlRestriction * rset)
{
    int i, j, m, err;

    err = vecm_restriction_check(rset);
    if (err) {
    return err;
    }

    if (rset.vecm & VecmSide.b) {
    rset.R = gretl_zero_matrix_new(rset.gb, rset.bcols);
    rset.q = gretl_zero_matrix_new(rset.gb, 1);
    if (rset.R == null || rset.q == null) {
        return GretlError.alloc;
    }
    }

    if (rset.vecm & VecmSide.a) {
    rset.Ra = gretl_zero_matrix_new(rset.ga, rset.acols);
    rset.qa = gretl_zero_matrix_new(rset.ga, 1);
    if (rset.Ra == null || rset.qa == null) {
        return GretlError.alloc;
    }
    }

    /* construct the restriction matrices, beta first then alpha
       (if both are given) */

    for (m=0; m<2; m++) {
    GretlMatrixRaw * R = (m == 0)? rset.R : rset.Ra;
    GretlMatrixRaw * q = (m == 0)? rset.q : rset.qa;
    char letter = (m == 0)? 'b' : 'a';
    int col, row = 0;

    if ((letter == 'b' && rset.vecm == VecmSide.a) ||
        (letter == 'a' && rset.vecm == VecmSide.b)) {
        continue;
    }

    for (i=0; i<rset.g; i++) {
        RRow * r = rset.rows[i];

        if (r.letter[0] == letter) {
        for (j=0; j<r.nterms; j++) {
            col = get_R_vecm_column(rset, i, j, letter);
            gretl_matrix_set(R, row, col, r.mult[j]);
        }
        gretl_vector_set(q, row, r.rhs);
        row++;
        }
    }

    }

    /* cumulate prior restrictions if needed */

    if (!err && beta_restricted_VECM(cast(GretlVar *) rset.obj)) {
    if (rset.vecm & VecmSide.b) {
        err = augment_vecm_restriction(rset, 'b');
    } else {
        err = add_old_vecm_restriction(rset, 'b');
    }
    }

    if (!err && alpha_restricted_VECM(cast(GretlVar *) rset.obj)) {
    if (rset.vecm & VecmSide.a) {
        err = augment_vecm_restriction(rset, 'a');
    } else {
        err = add_old_vecm_restriction(rset, 'a');
    }
    }

    if (!err && rset.R != null) {
    err = check_R_matrix(rset.R);
    }

    if (!err && rset.Ra != null) {
    err = check_R_matrix(rset.Ra);
    }

    return err;
}

/* We were given R and/or q directly by the user: check them
   for sanity. */

int test_user_matrices(GretlRestriction * rset)
{
    GretlMatrixRaw * R = rset.R;
    GretlMatrixRaw * q = rset.q;
    int err = 0;

    if (R == null || q == null) {
    /* we didn't get both parts */
    return GretlError.data;
    }

    if (R.rows != q.rows || q.cols != 1) {
    /* R and q don't work */
    return GretlError.nonconf;
    }

    if (rset.vecm) {
    /* this is supposed to be a restriction on beta:
       test it as such */
    int nb = gretl_VECM_n_beta(cast(GretlVar *) rset.obj);
    int nbr = nb * gretl_VECM_rank(cast(GretlVar *) rset.obj);

    if (R.cols == nb && R.rows <= nb) {
        /* common restriction on columns of beta */
        rset.bmulti = 0;
    } else if (R.cols == nbr && R.rows <= nbr) {
        ; /* OK, full restriction on beta */
    } else {
        gretl_errmsg_set(_("Invalid restriction on beta").ptr);
        err = GretlError.nonconf;
    }
    if (!err) {
        rset.bcols = R.cols;
    }
    } else {
    /* the non-VECM case */
    if (R.rows > rset.gmax) {
        gretl_errmsg_sprintf(_("Too many restrictions: the maximum number is %d").ptr, rset.gmax);
        err = GretlError.nonconf;
    } else if (R.cols != rset.gmax) {
        gretl_errmsg_sprintf(_("The R matrix must have %d columns").ptr,
                 rset.gmax);
        err = GretlError.nonconf;
    }
    }

    return err;
}

/* Set up the matrices needed for testing a set of restrictions:
   general driver function that calls more specific functions
   depending on the case. */

int restriction_set_form_matrices(GretlRestriction * rset)
{
    int err = 0;

    if (rset.R != null || rset.q != null) {
    err = test_user_matrices(rset);
    } else if (rset.otype == GretlObjType.eqn) {
    err = equation_form_matrices(rset);
    } else if (rset.otype == GretlObjType.var) {
    if (gretl_VECM_rank(cast(GretlVar *) rset.obj) > 0) {
        err = vecm_form_matrices(rset);
    } else {
        err = var_form_matrices(rset);
    }
    } else if (rset.g > 0) {
    /* mutivariate system */
    err = sys_form_matrices(rset);
    }

    return err;
}

/* Make a mask with 1s in positions in the array of coeffs where
   a coeff is referenced in one or more restrictions, 0s otherwise.
   We do this only for single-equation restriction sets, and
   only when the restrictions are given in row-wise form.
*/

int restriction_set_make_mask(GretlRestriction * rset)
{
    Model * pmod;
    RRow * r;
    int i, j;

    if (rset.otype != GretlObjType.eqn || rset.obj == null) {
    return GretlError.data;
    }

    pmod = cast(typeof(pmod)) rset.obj;
    rset.mask = cast(typeof(rset.mask)) calloc(pmod.ncoeff, 1);

    if (rset.mask == null) {
    destroy_restriction_set(rset);
    return GretlError.alloc;
    }

    for (i=0; i<rset.g; i++) {
    r = rset.rows[i];
    for (j=0; j<r.nterms; j++) {
        rset.mask[r.bnum[j]] = 1;
    }
    }

    return 0;
}

int count_ops(const(char) * p)
{
    int n = 0;

    while (*p) {
    if (*p == '+' || *p == '-') n++;
    if (*p == '=') break;
    p++;
    }

    return n;
}

int coeff_index_from_scalar(const(char) * s)
{
    int err = 0;
    double x = gretl_scalar_get_value(s, &err);

    return (err || na(x))? -1 : cast(int) x;
}

/* Try to retrieve the coefficient number in a model to be
   restricted, based on the parameter name, or possibly the
   name of a scalar variable.
*/

int bnum_from_name(GretlRestriction * r, const(DatasetRaw) * dset, const(char) * s)
{
    int k = -1;

    if (dset == null || r.obj == null) {
    ;
    } else if (r.otype == GretlObjType.var) {
    GretlVar * var = cast(GretlVar *) r.obj;
    int i, v;

    for (i=1; i<=var.ylist[0]; i++) {
        v = var.ylist[i];
        if (v > 0 && v < dset.v && !strcmp(dset.varname[v], s)) {
        k = i;
        break;
        }
    }
    } else if (r.otype == GretlObjType.eqn) {
    Model * pmod = cast(Model *) r.obj;

    k = gretl_model_get_param_number(pmod, dset, s);
    if (k >= 0) {
        /* convert to 1-based for compatibility with numbers read
           directly: the index will be converted to 0-base below
        */
        k++;
    }
    }

    if (k < 0) {
    /* try for a scalar variable */
    k = coeff_index_from_scalar(s);
    }
    if (k < 0) {
    gretl_errmsg_sprintf(_("%s: couldn't interpret as coefficient number").ptr, s);
    }

    return k;
}

/* Pick apart strings of the form "b[X]" or "b[X,Y]".  If the ",Y" is
   present the "X" element must be an equation number, and the "Y" may
   be a coefficient number or the name of a variable.  If the ",Y" is
   not present, "X" may be a coefficient number or the name of a
   variable.

   Such terms may be preceded by a multiplier, but they should not
   be followed by anything, so we flag an error if the closing "]"
   is not the end of the string.
*/

int pick_apart(GretlRestriction * r, const(char) * s, int * eq, int * bnum, const(DatasetRaw) * dset)
{
    int parens_ok = 0;
    const(char) * junk;
    char[VNAMELEN] s1 = 0;
    char[VNAMELEN] s2 = 0;
    char * targ = s1.ptr;
    char rdelim;
    int i, j, k;
    int err = 0;

    /* then skip past the opening delimiter */
    rdelim = (*s == '[')? ']' : ')';
    s++;

    *eq = *bnum = -1;

    k = gretl_charpos(rdelim, s);
    if (k <= 0 || k > 2*(VNAMELEN-1)) {
    return GretlError.parse;
    }

    junk = s + k + 1;
    while (isspace(*junk)) junk++;
    if (*junk != '\0') {
    fprintf(stderr, "Got trailing junk '%s'\n".ptr, junk);
    return GretlError.parse;
    }

    j = 0;
    for (i=0; i<k; i++) {
    if (s[i] == ',') {
        targ = s2.ptr;
        j = 0;
    } else if (!isspace(s[i])) {
        if (j == VNAMELEN-1) {
        return GretlError.parse;
        }
        targ[j++] = s[i];
    }
    }

    if (targ == s2.ptr) {
    /* got a comma separator: [eqn,bnum] */
    *eq = positive_int_from_string(s1.ptr);
    if (*eq <= 0) {
        err = GretlError.parse;
    }
    if (isdigit(*s2.ptr)) {
        *bnum = positive_int_from_string(s2.ptr);
    } else {
        *bnum = bnum_from_name(r, dset, s2.ptr);
    }
    } else {
    /* only one field: [bnum] */
    *eq = R_UNSPEC;
    if (isdigit(*s1.ptr)) {
        *bnum = positive_int_from_string(s1.ptr);
    } else {
        parens_ok = 1;
        *bnum = bnum_from_name(r, dset, s1.ptr);
    }
    }

    if (!err && rdelim == ')' && !parens_ok) {
    err = GretlError.parse;
    }

    return err;
}

int bbit_trailing_garbage(const(char) * s)
{
    if (*s) {
    while (isspace(*s)) s++;
    if (*s) {
        if (*s == '/') {
        gretl_errmsg_sprintf(_("%s: division not allowed here").ptr, s-1);
        }
        return 1;
    }
    }

    return 0;
}

/* We got a leading letter ('a' or 'b') and now we parse
   what follows.  The simplest case is just a number, giving
   e.g., "b1".  Then we assess the validity of what we've
   found.
*/

int parse_b_bit(GretlRestriction * r, const(char) * s, int * eq, int * bnum, const(DatasetRaw) * dset)
{
    int err = GretlError.parse;

    if (isdigit(cast(ubyte) *s)) {
    char * test;

    *bnum = cast(int) strtol(s, &test, 10);
    if (bbit_trailing_garbage(test)) {
        return err;
    }
    if (r.otype == GretlObjType.var) {
        *eq = R_UNSPEC;
    }
    err = 0;
    } else if (*s == '[' || *s == '(') {
    err = pick_apart(r, s, eq, bnum, dset);
    }

    if (*bnum < 1) {
    if (!gretl_errmsg_is_set() && *bnum >= 0) {
        gretl_errmsg_sprintf(_("Coefficient number (%d) is out of range").ptr,
                 *bnum);
    }
    err = 1;
    } else {
    *bnum -= 1; /* convert to zero base */
    }

    if (*eq == R_UNSPEC) {
    /* didn't get a specific equation number */
    if (r.otype == GretlObjType.eqn) {
        *eq = 0;
    } else if (r.otype != GretlObjType.var) {
        err = GretlError.parse;
    }
    } else if (*eq == 0) {
    gretl_errmsg_sprintf(_("Equation number (%d) is out of range").ptr,
                 *eq);
    err = 1;
    } else {
    *eq -= 1; /* convert to zero base */
    }

    return err;
}

int r_get_scalar(const(char) * s, double * x)
{
    if (sscanf(s, "%lf".ptr, x)) {
    return 1;
    } else {
    int n;

    s += cast(int) strspn(s, " ");
    n = gretl_namechar_spn(s);

    if (n > 0) {
        char[VNAMELEN] tmp;

        *tmp.ptr = '\0';
        strncat(tmp.ptr, s, n);
        if (gretl_is_scalar(tmp.ptr)) {
        *x = gretl_scalar_get_value(tmp.ptr, null);
        if (!na(*x)) {
            return 1;
        }
        }
    }
    }

    return 0;
}

int get_named_param(const(char) * s, Model * pmod, int * bnum)
{
    if (NONLIST_MODEL(cast(GretlCmdIndex) pmod.ci) && pmod.params != null) {
    char[VNAMELEN] test;
    int i;

    *test.ptr = '\0';
    strncat(test.ptr, s, VNAMELEN - 1);
    tailstrip(test.ptr);

    for (i=0; i<pmod.nparams; i++) {
        if (!strcmp(test.ptr, pmod.params[i])) {
        *bnum = i;
        return 1;
        }
    }
    }

    return 0;
}

int starts_b(const(char) * s, GretlRestriction * r)
{
    if (*s == 'b' || (r.vecm && *s == 'a')) {
    char c = *(s+1);

    return isdigit(c) || c == '[' || c == '(';
    }

    return 0;
}

bool could_be_param(GretlRestriction * r, const(char) * s)
{
    return (r.otype == GretlObjType.eqn && !isdigit(*s));
}

int parse_coeff_chunk(GretlRestriction * r, const(char) * s, double * x, int * eq, int * bnum, char * letter, const(DatasetRaw) * dset)
{
    const(char) * s0 = s;
    int done = 0;
    int err = GretlError.parse;

    *eq = 1;

    while (isspace(cast(ubyte) *s)) s++;

    if (starts_b(s, r)) {
    /* something like "b[" or "b2" */
    *letter = *s;
    s++;
    err = parse_b_bit(r, s, eq, bnum, dset);
    *x = 1.0;
    done = 1;
    } else if (could_be_param(r, s)) {
    /* e.g. "beta" for MLE or GMM */
    done = get_named_param(s, cast(Model *) r.obj, bnum);
    if (done) {
        *letter = 'b';
        *x = 1.0;
        *eq = 0;
        err = 0;
    }
    }

    if (!done && r_get_scalar(s, x)) {
    s += cast(int) strspn(s, " ");
    s += cast(int) strcspn(s, " *");
    s += cast(int) strspn(s, " *");
    if (starts_b(s, r)) {
        *letter = *s;
        s++;
        err = parse_b_bit(r, s, eq, bnum, dset);
    } else if (*s == '\0') {
        /* plain numeric value, saved in *x */
        *bnum = 0;
        *letter = '\0';
        err = 0;
    }
    }

    if (err) {
    gretl_errmsg_sprintf(_("parse error in '%s'\n").ptr, s0);
    }

    return err;
}

void destroy_restriction(RRow * r)
{
    if (r != null) {
    free(r.mult);
    free(r.eq);
    free(r.bnum);
    free(r.letter);
    free(r);
    }
}

void destroy_restriction_set(GretlRestriction * rset)
{
    if (rset == null) {
    return;
    }
    if (rset.rows != null) {
    int i;

    for (i=0; i<rset.g; i++) {
        destroy_restriction(rset.rows[i]);
    }
    free(rset.rows);
    }

    free(rset.mask);

    gretl_matrix_free(rset.R);
    gretl_matrix_free(rset.q);
    gretl_matrix_free(rset.Ra);
    gretl_matrix_free(rset.qa);

    free(rset.rfunc);

    free(rset);
}

RRow * restriction_new(int n, int multi, int vecm)
{
    RRow * r;
    int i;

    r = cast(typeof(r)) malloc((*r).sizeof);
    if (r == null) {
    return null;
    }

    r.mult = null;
    r.eq = null;
    r.bnum = null;
    r.letter = null;

    r.mult = cast(typeof(r.mult)) malloc(n * (*r.mult).sizeof);
    r.bnum = cast(typeof(r.bnum)) malloc(n * (*r.bnum).sizeof);
    if (r.mult == null || r.bnum == null) {
    destroy_restriction(r);
    return null;
    }

    for (i=0; i<n; i++) {
    r.mult[i] = 0.0;
    r.bnum[i] = 0;
    }

    if (multi) {
    r.eq = cast(typeof(r.eq)) malloc(n * (*r.eq).sizeof);
    if (r.eq == null) {
        destroy_restriction(r);
        return null;
    }
    }

    if (vecm) {
    r.letter = cast(typeof(r.letter)) malloc(n * (*r.letter).sizeof);
    if (r.letter == null) {
        destroy_restriction(r);
        return null;
    }
    }

    r.nterms = n;
    r.rhs = 0.0;

    return r;
}

RRow * restriction_set_add_row(GretlRestriction * rset, int n_terms)
{
    RRow ** rows = null;
    int n = rset.g;

    rows = cast(typeof(rows)) realloc(rset.rows, (n + 1) * (*rows).sizeof);
    if (rows == null) {
    return null;
    }

    rset.rows = rows;

    rset.rows[n] = restriction_new(n_terms, rset.bmulti, rset.vecm);
    if (rset.rows[n] == null) {
    return null;
    }

    rset.g += 1;

    return rset.rows[n];
}

void print_mult(double mult, int first, Prn * prn)
{
    if (mult == 1.0) {
    if (!first) pputs(prn, " + ");
    } else if (mult == -1.0) {
    if (first) pputs(prn, "-");
    else pputs(prn, " - ");
    } else if (mult > 0.0) {
    if (first) pprintf(prn, "%g*".ptr, mult);
    else pprintf(prn, " + %g*".ptr, mult);
    } else if (mult < 0.0) {
    if (first) pprintf(prn, "%g*".ptr, mult);
    else pprintf(prn, " - %g*".ptr, fabs(mult));
    }
}

void print_restriction(const(GretlRestriction) * rset, int i, const(DatasetRaw) * dset, Prn * prn)
{
    const(RRow) * r = rset.rows[i];
    char letter;
    char[VNAMELEN] vname;
    int j, k;

    for (j=0; j<r.nterms; j++) {
    letter = (r.letter != null)? r.letter[j] : 'b';
    k = r.bnum[j];
    print_mult(r.mult[j], j == 0, prn);
    if (targ_specified(rset, i, j, letter)) {
        pprintf(prn, "%c[%d,%d]".ptr, letter, r.eq[j] + 1, k + 1);
    } else if (rset.otype == GretlObjType.var) {
        pprintf(prn, "%c[%d]".ptr, letter, k + 1);
    } else {
        Model * pmod = cast(Model *) rset.obj;

        gretl_model_get_param_name(pmod, dset, k, vname.ptr);
        if (NONLIST_MODEL(cast(GretlCmdIndex) pmod.ci)) {
        pputs(prn, vname.ptr);
        } else {
        pprintf(prn, "b[%s]".ptr, vname.ptr);
        }
    }
    }

    pprintf(prn, " = %g\n".ptr, r.rhs);
}

/* Pretty-print the set of restrictions.  We do this only if
   the restrictions were given row-wise by the user.
*/

void print_restriction_set(const(GretlRestriction) * rset, const(DatasetRaw) * dset, Gretlopt opt, Prn * prn)
{
    int i;

    if (prn == null) {
    return;
    }

    if (rset.g > 1) {
    pputs(prn, _("Restriction set").ptr);
    } else {
    pprintf(prn, "%s:".ptr, _("Restriction").ptr);
    }
    pputc(prn, '\n');

    for (i=0; i<rset.g; i++) {
    if (rset.g > 1) {
        pprintf(prn, " %d: ".ptr, i + 1);
    } else {
        pputc(prn, ' ');
    }
    print_restriction(rset, i, dset, prn);
    }

    if (opt & Gretlopt.OPT_V) {
    if (rset.R != null || rset.q != null ||
        rset.Ra != null || rset.qa != null) {
        pputc(prn, '\n');
    }
    if (rset.R != null) {
        gretl_matrix_print_to_prn(rset.R, "R", prn);
    }
    if (rset.q != null) {
        gretl_matrix_print_to_prn(rset.q, "q", prn);
    }
    if (rset.Ra != null) {
        gretl_matrix_print_to_prn(rset.Ra, "R(alpha)", prn);
    }
    if (rset.qa != null) {
        gretl_matrix_print_to_prn(rset.qa, "q(alpha)", prn);
    }
    }
}

void print_restriction_from_matrices(const(GretlMatrixRaw) * R, const(GretlMatrixRaw) * q, char letter, int npar, Prn * prn)
{
    double x;
    int eqn, coeff, started;
    int i, j;

    for (i=0; i<R.rows; i++) {
    started = 0;
    coeff = 1;
    eqn = (R.cols > npar)? 1 : 0;
    for (j=0; j<R.cols; j++) {
        x = gretl_matrix_get(R, i, j);
        if (x != 0.0) {
        if (!started) {
            pputs(prn, "  ");
        }
        if (x == 1.0) {
            if (started) {
            pputs(prn, " + ");
            }
        } else if (x == -1.0) {
            if (started) {
            pputs(prn, " - ");
            } else {
            pputc(prn, '-');
            }
        } else if (x > 0.0) {
            if (started) {
            pprintf(prn, " + %g*".ptr, x);
            } else {
            pprintf(prn, "%g*".ptr, x);
            }
        } else if (x < 0.0) {
            if (started) {
            pprintf(prn, " - %g*".ptr, -x);
            } else {
            pprintf(prn, "%g*".ptr, x);
            }
        }
        if (eqn > 0) {
            pprintf(prn, "%c[%d,%d]".ptr, letter, eqn, coeff);
        } else {
            pprintf(prn, "%c%d".ptr, letter, coeff);
        }
        started = 1;
        }
        if ((j + 1) % npar == 0) {
        eqn++;
        coeff = 1;
        } else {
        coeff++;
        }
    }
    pprintf(prn, " = %g\n".ptr, (q == null)? 0.0 : q.val[i]);
    }
}

int add_term_to_restriction(RRow * r, double mult, int eq, int bnum, char letter, int i)
{
    int j;

    for (j=0; j<i; j++) {
    if (r.letter != null && r.letter[j] != letter) {
        continue;
    }
    if (bnum == r.bnum[j] && (r.eq == null || eq == r.eq[j])) {
        /* additional reference to a previously referenced coeff */
        r.mult[j] += mult;
        r.nterms -= 1;
        return 0;
    }
    }

    r.mult[i] = mult;
    r.bnum[i] = bnum;

    if (r.eq != null) {
    r.eq[i] = eq;
    }

    if (r.letter != null) {
    r.letter[i] = letter;
    }

    return 0;
}

GretlRestriction * restriction_set_new(void * ptr, GretlObjType type, Gretlopt opt)
{
    GretlRestriction * rset;

    rset = cast(typeof(rset)) malloc((*rset).sizeof);
    if (rset == null) return null;

    rset.obj = ptr;
    rset.otype = type;
    rset.opt = opt;

    rset.test = NADBL;
    rset.pval = NADBL;
    rset.lnl = NADBL;
    rset.bsum = NADBL;
    rset.bse = NADBL;

    rset.g = rset.gmax = 0;
    rset.gb = rset.ga = 0;
    rset.R = null;
    rset.q = null;
    rset.Ra = null;
    rset.qa = null;
    rset.mask = null;
    rset.rows = null;

    rset.rfunc = null;

    rset.bmulti = 0;
    rset.amulti = 0;
    rset.bcols = 0;
    rset.acols = 0;
    rset.vecm = 0;
    rset.code = TestStats.none;

    if (rset.otype == GretlObjType.eqn) {
    Model * pmod = cast(Model *) ptr;

    rset.gmax = pmod.ncoeff;
    } else if (rset.otype == GretlObjType.sys) {
    rset.gmax = system_n_indep_vars(cast(EquationSystem *) ptr);
    rset.bmulti = 1;
    } else if (rset.otype == GretlObjType.var) {
    GretlVar * var = cast(GretlVar *) ptr;
    int r = gretl_VECM_rank(var);

    if (r == 0) {
        /* a plain VAR */
        rset.gmax = var.neqns * var.ncoeff;
        rset.bmulti = 1;
    } else {
        if (r > 1) {
        rset.bmulti = 1;
        rset.amulti = 1;
        }
        rset.vecm = VecmSide.b;
        rset.gmax = gretl_VECM_n_beta(var) * r;
        rset.gmax += gretl_VECM_n_alpha(var) * r;
    }
    }

    return rset;
}

GretlRestriction * rset_from_VECM(GretlVar * var, int * err)
{
    const(GretlMatrixRaw) * R = gretl_VECM_R_matrix(var);
    const(GretlMatrixRaw) * q = gretl_VECM_q_matrix(var);
    const(GretlMatrixRaw) * Ra = gretl_VECM_Ra_matrix(var);
    const(GretlMatrixRaw) * qa = gretl_VECM_qa_matrix(var);
    GretlRestriction * rset;

    if (R == null && q == null && Ra == null && qa == null) {
    /* no restriction in place */
    return null;
    }

    rset = cast(typeof(rset)) calloc(1, (*rset).sizeof);
    if (rset == null) {
    *err = GretlError.alloc;
    return null;
    }

    rset.obj = var;
    rset.otype = GretlObjType.var;
    rset.opt = Gretlopt.OPT_B;

    rset.R = cast(GretlMatrixRaw *) R;
    rset.q = cast(GretlMatrixRaw *) q;
    rset.Ra = cast(GretlMatrixRaw *) Ra;
    rset.qa = cast(GretlMatrixRaw *) qa;

    return rset;
}

int list_get_bmax(const(int) * list)
{
    int i, bmax = 0;

    for (i=2; i<=list[0]; i++) {
    if (list[i] == LISTSEP) {
        break;
    }
    bmax++;
    }

    return bmax;
}

/* check that the coefficients referenced in a restriction are
   within bounds, relative to the equation or system that is
   to be restricted */

int bnum_out_of_bounds(const(GretlRestriction) * rset, int eq, int bnum, char letter)
{
    int ret = 1;

    if (rset.otype == GretlObjType.var) {
    GretlVar * var = cast(GretlVar *) rset.obj;
    int r = gretl_VECM_rank(var);

    if (r == 0) {
        if (eq >= var.neqns) {
        gretl_errmsg_sprintf(_("Equation number (%d) is out of range").ptr,
                     eq + 1);
        } else if (bnum >= var.ncoeff) {
        gretl_errmsg_sprintf(_("Coefficient number (%d) is out of range").ptr,
                     bnum + 1);
        } else {
        ret = 0;
        }
    } else {
        /* note: @eq and @bnum are transposed here */
        if (eq >= gretl_VECM_rank(var)) {
        gretl_errmsg_sprintf(_("Coefficient number (%d) is out of range").ptr,
                     eq + 1);
        } else if ((letter == 'b' && bnum >= gretl_VECM_n_beta(var)) ||
               (letter == 'a' && bnum >= gretl_VECM_n_alpha(var))) {
        gretl_errmsg_sprintf(_("Equation number (%d) is out of range").ptr,
                     bnum + 1);
        } else {
        ret = 0;
        }
    }
    } else if (rset.otype == GretlObjType.sys) {
    EquationSystem * sys = cast(EquationSystem *) rset.obj;
    const(int) * list = system_get_list(sys, eq);

    if (list == null) {
        gretl_errmsg_sprintf(_("Equation number (%d) is out of range").ptr,
                 eq + 1);
    } else {
        int bmax = list_get_bmax(list);

        if (bnum >= bmax) {
        gretl_errmsg_sprintf(_("Coefficient number (%d) out of range for equation %d").ptr, bnum + 1, eq + 1);
        } else {
        ret = 0;
        }
    }
    } else {
    Model * pmod = cast(Model *) rset.obj;

    if (eq > 0) {
        gretl_errmsg_sprintf(_("Equation number (%d) is out of range").ptr,
                 eq + 1);
    } else if (bnum >= pmod.ncoeff || bnum < 0) {
        if (!gretl_errmsg_is_set()) {
        gretl_errmsg_sprintf(_("Coefficient number (%d) is out of range").ptr,
                     bnum + 1);
        }
    } else {
        ret = 0;
    }
    }

    return ret;
}

/* Parse a matrix specification for a restriction set.  This should
   take the form "R = <matrix>" or "q = <matrix>", where <matrix> is
   the name of a pre-defined matrix.  On this approach we want exactly
   one R definition and exactly one q definition, and no other sort
   of specification can be given.
*/

int read_matrix_line(const(char) * s, GretlRestriction * rset)
{
    const(GretlMatrixRaw) * m;
    char[VNAMELEN] mname;
    char job = *s;
    int err = 0;

    if (rset.rows != null) {
    /* we already have "row-wise" restrictions */
    return GretlError.parse;
    } else if (job == 'R' && rset.R != null) {
    /* duplicate entry for R matrix */
    return GretlError.parse;
    } else if (job == 'q' && rset.q != null) {
    /* duplicate entry for q matrix */
    return GretlError.parse;
    }

    s++;
    while (isspace(cast(ubyte) *s)) s++;
    if (*s != '=') {
    return GretlError.parse;
    }

    s++;
    while (isspace(cast(ubyte) *s)) s++;
    if (gretl_scan_varname(s, mname.ptr) != 1) {
    return GretlError.parse;
    }

    m = get_matrix_by_name(mname.ptr);
    if (m == null) {
    return GretlError.unkvar;
    }

    if (job == 'R') {
    rset.R = gretl_matrix_copy(m);
    if (rset.R == null) {
        err = GretlError.alloc;
    }
    } else if (job == 'q') {
    rset.q = gretl_matrix_copy(m);
    if (rset.q == null) {
        err = GretlError.alloc;
    }
    }

    return err;
}

/* Try to read the name of a function to be called for testing
   a nonlinear restriction, as in "rfunc = <fname>".  This is
   admissable only if it's not mixed with any other sort
   of restriction.
*/

int read_fncall_line(const(char) * s, GretlRestriction * rset)
{
    char[FN_NAMELEN] fname;

    if (rset.rows != null || rset.R != null || rset.q != null) {
    /* other stuff is in the way! */
    return GretlError.parse;
    }

    s += cast(int) strspn(s, " ");

    if (*s != '=') {
    return GretlError.parse;
    }

    s++;
    s += cast(int) strspn(s, " ");

    if (sscanf(s, "%31s".ptr, fname.ptr) != 1) {
    return GretlError.parse;
    }

    s += strlen(fname.ptr);

    if (!string_is_blank(s)) {
    /* there should be no trailing junk */
    return GretlError.parse;
    }

    rset.rfunc = gretl_strdup(fname.ptr);

    return 0;
}

/* We're going to count operators to figure out how many
   terms we have in a restriction, so we'd better check
   first that we don't have any spurious pile-ups of
   operators (that is, '-', '+' and '*'), none of which
   can be repeated consecutively, or separated only by
   spaces.
*/

bool IS_ROP(int c)
{
    return (c=='-' || c == '+' || c == '*');
}

int restriction_op_check(const(char) * s)
{
    int err = 0;

    while (*s && !err) {
    if (IS_ROP(*s)) {
        s++;
        while (*s && !err) {
        if (*s == ' ') {
            ; /* OK so far */
        } else if (!(IS_ROP(*s))) {
            break;
        } else {
            err = GretlError.parse;
        }
        s++;
        }
    }
    s++;
    }

    return err;
}

int parse_restriction_row(GretlRestriction * rset, const(char) * s, const(DatasetRaw) * dset)
{
    RRow * row;
    int nt, sgn = 1;
    int i, j;
    int err = 0;

    if (rset.R != null || rset.q != null) {
    /* can't mix row-wise spec with matrix spec */
    return GretlError.parse;
    }

    err = restriction_op_check(s);
    if (err) {
    return err;
    }

    if (*s == '+' || *s == '-') {
    sgn = (*s == '+')? 1 : -1;
    s++;
    }

    nt = 1 + count_ops(s);

    row = restriction_set_add_row(rset, nt);
    if (row == null) {
    return GretlError.alloc;
    }

    j = 0;

    for (i=0; i<nt; i++) {
    /* work on the left-hand side */
    char[36] chunk;
    int len, bnum = 1, eq = 1;
    int numeric = 0;
    char letter = 'b';
    double mult;

    len = cast(int) strcspn(s, "+-=");
    if (len >= chunk.sizeof) {
        err = 1;
        break;
    }

    *chunk.ptr = '\0';
    strncat(chunk.ptr, s, len);
    s += len;

    err = parse_coeff_chunk(rset, chunk.ptr, &mult, &eq, &bnum,
                &letter, dset);
    if (err) {
        break;
    } else if (letter == '\0') {
        numeric = 1;
    } else if (bnum_out_of_bounds(rset, eq, bnum, letter)) {
        err = GretlError.data;
        break;
    }

    mult *= sgn;

    if (numeric) {
        /* got a numeric constant */
        row.rhs -= mult;
        row.nterms -= 1;
    } else {
        add_term_to_restriction(row, mult, eq, bnum, letter, j++);
    }

    if (*s == '+') {
        sgn = 1;
        s++;
    } else if (*s == '-') {
        sgn = -1;
        s++;
    }
    }

    if (!err) {
    /* work on the right-hand side */
    double rhs = 0.0;

    s += cast(int) strspn(s, " ");

    if (*s != '=') {
        err = GretlError.parse;
    } else {
        /* we should find here a scalar constant or a
           variable of scalar type */
        rhs = generate_scalar(s+1, null, &err);
        if (!err) {
        if (na(rhs)) {
            err = GretlError.data;
        } else {
            row.rhs += rhs;
        }
        }
    }
    }

    if (!err && row.nterms == 0) {
    gretl_errmsg_set(_("Invalid restriction").ptr);
    err = GretlError.parse;
    }

    return err;
}

int parse_strings_line(const(char) * line, GretlRestriction ** prset, const(DatasetRaw) * dset)
{
    char[VNAMELEN] aname;
    int err = 0;

    if (sscanf(line, "%31s".ptr, aname.ptr) != 1) {
    err = GretlError.parse;
    } else {
    GretlArray * a = get_strings_array_by_name(aname.ptr);

    if (a == null) {
        err = GretlError.invarg;
    } else {
        char ** S;
        int i, n;

        S = gretl_array_get_strings(a, &n);
        for (i=0; i<n && !err; i++) {
        err = real_restriction_set_parse_line(prset, S[i], dset, 0);
        }
    }
    }

    return err;
}

int real_restriction_set_parse_line(GretlRestriction ** prset, const(char) * line, const(DatasetRaw) * dset, int first)
{
    GretlRestriction * rset = *prset;
    const(char) * s = line;
    int err = 0;

    gretl_error_clear();

    /* if we have a function name specified already, nothing else
       should be supplied */
    if (rset.rfunc != null) {
    destroy_restriction_set(rset);
    *prset = null;
    return GretlError.parse;
    }

    if (!strncmp(s, "restrict", 8)) {
    if (strlen(line) == 8) {
        if (first) {
        return 0;
        } else {
        return GretlError.parse;
        }
    }
    s += 8;
    while (isspace(cast(ubyte) *s)) s++;
    }

    if (!strncmp(s, "rfunc", 5)) {
    /* nonlinear function given */
    err = read_fncall_line(s + 5, rset);
    } else if (*s == 'R' || *s == 'q') {
    /* restrictions given in matrix form */
    err = read_matrix_line(s, rset);
    } else if (!strncmp(s, "inject", 6)) {
    /* restrictions given as an array of strings */
    err = parse_strings_line(s + 6, prset, dset);
    } else {
    /* a regular linear restriction */
    err = parse_restriction_row(rset, s, dset);
    }

    if (err && *prset != null) {
    destroy_restriction_set(*prset);
    *prset = null;
    }

    return err;
}

int restriction_set_parse_line(GretlRestriction * rset, const(char) * line, const(DatasetRaw) * dset)
{
    if (rset.g > rset.gmax) {
    gretl_errmsg_sprintf(_("Too many restrictions (maximum is %d)").ptr,
                 rset.gmax);
    destroy_restriction_set(rset);
    return GretlError.data;
    }

    return real_restriction_set_parse_line(&rset, line, dset, 0);
}

/* set-up for a set of restrictions for a VAR (vecm, actually) */

GretlRestriction * var_restriction_set_start(const(char) * line, GretlVar * var)
{
    GretlRestriction * rset;

    rset = restriction_set_new(var, GretlObjType.var, Gretlopt.OPT_NONE);
    if (rset == null) {
    gretl_errmsg_set(_("Out of memory!").ptr);
    return null;
    }

    gretl_error_clear();

    if (real_restriction_set_parse_line(&rset, line, null, 1)) {
    if (!gretl_errmsg_is_set()) {
        gretl_errmsg_sprintf(_("parse error in '%s'\n").ptr, line);
    }
    return null;
    }

    return rset;
}

/* set-up for a set of (possibly) cross-equation restrictions, for a
   system of simultaneous equations */

GretlRestriction * cross_restriction_set_start(const(char) * line, EquationSystem * sys)
{
    GretlRestriction * rset;

    rset = restriction_set_new(sys, GretlObjType.sys, Gretlopt.OPT_NONE);
    if (rset == null) {
    gretl_errmsg_set(_("Out of memory!").ptr);
    return null;
    }

    if (real_restriction_set_parse_line(&rset, line, null, 1)) {
    gretl_errmsg_sprintf(_("parse error in '%s'\n").ptr, line);
    return null;
    }

    return rset;
}

/* set-up for a set of restrictions on a single equation */

GretlRestriction * eqn_restriction_set_start(const(char) * line, Model * pmod, const(DatasetRaw) * dset, Gretlopt opt)
{
    GretlRestriction * rset;

    rset = restriction_set_new(pmod, GretlObjType.eqn, opt);
    if (rset == null) {
    gretl_errmsg_set(_("Out of memory!").ptr);
    return null;
    }

    if (real_restriction_set_parse_line(&rset, line, dset, 1)) {
    gretl_errmsg_sprintf(_("parse error in '%s'\n").ptr, line);
    return null;
    }

    return rset;
}

GretlRestriction * restriction_set_start(const(char) * target, Gretlopt opt, int * err)
{
    GretlRestriction * rset = null;
    GretlObjType type = cast(GretlObjType) 0;
    void * ptr = null;

    if (target != null) {
    /* get pointer to named object */
    if (!strcmp(target, "$system")) {
        ptr = get_anonymous_equation_system();
        if (ptr != null) {
        type = GretlObjType.sys;
        }
    } else {
        ptr = get_model_object_and_type(target, &type);
    }
    if (ptr == null) {
        gretl_errmsg_sprintf("'%s': unrecognized target".ptr, target);
        *err = GretlError.data;
    }
    } else {
    /* get pointer to last-estimated model */
    ptr = get_last_model(&type);
    if (ptr == null) {
        *err = GretlError.data;
    }
    }

    if (!*err) {
    if (type != GretlObjType.eqn && type != GretlObjType.sys &&
        type != GretlObjType.var) {
        *err = GretlError.data;
    } else {
        rset = restriction_set_new(ptr, type, opt);
        if (rset == null) {
        *err = GretlError.alloc;
        }
    }
    }

    return rset;
}

int print_restricted_estimates(Model * pmod, const(DatasetRaw) * dset, GretlMatrixRaw * vb, GretlMatrixRaw * S, double s2, int nc, int rk, Prn * prn)
{
    const(double) * b = vb.val;
    char ** names;
    double v;
    double * se;
    int i;

    names = strings_array_new_with_length(nc, VNAMELEN);
    if (names == null) {
    return GretlError.alloc;
    }

    se = cast(typeof(se)) malloc(nc * (*se).sizeof);
    if (se == null) {
    strings_array_free(names, nc);
    return GretlError.alloc;
    }

    pprintf(prn, "%s:\n\n".ptr, _("Restricted estimates").ptr);

    for (i=0; i<nc; i++) {
    v = gretl_matrix_get(S, i, i);
    if (fabs(v) < 1.0e-16) {
        se[i] = 0.0;
    } else {
        se[i] = sqrt(v);
    }
    gretl_model_get_param_name(pmod, dset, i, names[i]);
    }

    print_coeffs(b, se, cast(const(char) **) names, nc, pmod.dfd + rk,
         pmod.ci, prn);

    pputc(prn, '\n');
    pprintf(prn, "  %s = %.*g\n".ptr, _("Standard error of the regression").ptr,
        get_gretl_digits(), sqrt(s2));

    strings_array_free(names, nc);
    free(se);

    return 0;
}

/* Used when generating restricted estimates or bootstrapping, which
   we do for single-equation OLS/WLS only. */

int rset_expand_R(GretlRestriction * rset, int k)
{
    GretlMatrixRaw * R = null;
    double rij;
    int i, j, jj;

    R = gretl_zero_matrix_new(rset.g, k);
    if (R == null) {
    return GretlError.alloc;
    }

    for (i=0; i<rset.g; i++) {
    jj = 0;
    for (j=0; j<k; j++) {
        if (rset.mask[j]) {
        rij = gretl_matrix_get(rset.R, i, jj);
        gretl_matrix_set(R, i, j, rij);
        jj++;
        }
    }
    }

    gretl_matrix_replace(&rset.R, R);

    return 0;
}

void rmod_check_ifc(Model * rmod, const(GretlMatrixRaw) * y)
{
    double x1 = 0.0, x2 = 0.0;
    double reldiff;
    int t, s = 0;

    for (t=rmod.t1; t<=rmod.t2; t++) {
    if (!na(rmod.yhat[t])) {
        x1 += rmod.yhat[t];
        x2 += y.val[s++];
    }
    }

    reldiff = fabs((x1 - x2) / x2);
    rmod.ifc = reldiff < 9.0e-16;
}

/* Respond to the --full option in case of a restriction on a model
   estimated via OLS: arrange to replace the "last model" with the
   restricted OLS estimates (and to print the model if wanted).
*/

int save_restricted_model(ExecState * state, Model * pmod, const(DatasetRaw) * dset, const(GretlMatrixRaw) * y, const(GretlMatrixRaw) * b, const(GretlMatrixRaw) * S, const(GretlMatrixRaw) * u, double s2, GretlRestriction * rset, Prn * prn)
{
    Model * rmod;
    int i, err = 0;

    if (state == null) {
    return GretlError.data;
    }

    rmod = gretl_model_new();
    if (rmod == null) {
    return GretlError.alloc;
    }

    rmod.ci = GretlCmdIndex.ols;
    rmod.ncoeff = gretl_vector_get_length(b);
    rmod.nobs = gretl_vector_get_length(y);
    rmod.full_n = dset.n;
    rmod.t1 = pmod.t1;
    rmod.t2 = pmod.t2;

    rmod.dfn = pmod.dfn - rset.g;
    rmod.dfd = pmod.dfd + rset.g;
    rmod.ybar = pmod.ybar;
    rmod.sdy = pmod.sdy;
    rmod.tss = pmod.tss;

    gretl_model_smpl_init(rmod, dset);

    err = gretl_model_allocate_storage(rmod);

    if (!err) {
    for (i=0; i<rmod.ncoeff; i++) {
        rmod.coeff[i] = b.val[i];
    }
    gretl_model_set_int(rmod, "restricted", 1);
    err = gretl_model_write_vcv(rmod, S);
    }

    if (!err) {
    rmod.list = gretl_list_copy(pmod.list);
    if (rmod.list == null) {
        err = GretlError.alloc;
    } else {
        rmod.xlist = gretl_list_sublist(rmod.list, 2, -1);
    }
    }

    if (!err) {
    double yt, ut;
    int t, s = 0;

    rmod.sigma = sqrt(s2);
    rmod.ess = 0.0;

    for (t=0; t<dset.n; t++) {
        if (!na(pmod.uhat[t])) {
        yt = gretl_vector_get(y, s);
        ut = gretl_vector_get(u, s);
        rmod.uhat[t] = ut;
        rmod.yhat[t] = yt - ut;
        rmod.ess += ut * ut;
        s++;
        }
    }

    rmod_check_ifc(rmod, y); /* FIXME? */

    if (rmod.ifc) {
        double tmp;

        rmod.rsq = 1.0 - (rmod.ess / rmod.tss);
        if (rmod.rsq < 0.0) {
        rmod.rsq = 0.0;
        }
        if (rmod.dfn == 0) {
        rmod.adjrsq = 0.0;
        rmod.fstt = NADBL;
        } else {
        tmp = rmod.tss * rmod.dfd;
        rmod.adjrsq = 1 - (rmod.ess * (rmod.nobs - 1) / tmp);
        tmp = rmod.dfd / cast(double) rmod.dfn;
        rmod.fstt = tmp * (rmod.tss - rmod.ess) / rmod.ess;
        if (rmod.fstt < 0) {
            rmod.fstt = NADBL;
        }
        }
    }

    rmod.ncoeff -= rset.g;
    ls_criteria(rmod);
    rmod.ncoeff += rset.g;
    }

    if (err) {
    gretl_model_free(rmod);
    } else {
    set_model_id(rmod, Gretlopt.OPT_NONE);
    gretl_exec_state_set_model(state, rmod);
    }

    return err;
}

/* generate full restricted estimates: this function is used only for
   single-equation models estimated via OLS
*/

int do_restricted_estimates(ExecState * state, GretlRestriction * rset, const(DatasetRaw) * dset, Prn * prn)
{
    Model * pmod = cast(Model *) rset.obj;
    GretlMatrixBlock * B;
    GretlMatrixRaw * X, y, b, S;
    GretlMatrixRaw * u = null;
    double s2 = 0.0;
    int T = pmod.nobs;
    int k = pmod.ncoeff;
    int i, s, t;
    int yno, err = 0;

    B = gretl_matrix_block_new(&X, T, k,
                   &y, T, 1,
                   &b, k, 1,
                   &S, k, k,
                   null);

    if (B == null) {
    return GretlError.alloc;
    }

    if (gretl_matrix_cols(rset.R) != k) {
    err = rset_expand_R(rset, k);
    if (err) {
        goto bailout;
    }
    }

    yno = gretl_model_get_depvar(pmod);
    if (pmod.xlist == null) {
    err = GretlError.data;
    goto bailout;
    }

    if (rset.opt & Gretlopt.OPT_F) {
    /* full save wanted */
    u = gretl_matrix_alloc(T, 1);
    if (u == null) {
        err = GretlError.alloc;
        goto bailout;
    }
    }

    s = 0;
    for (t=pmod.t1; t<=pmod.t2; t++) {
    if (na(pmod.uhat[t])) {
        continue;
    }
    gretl_vector_set(y, s, dset.Z[yno][t]);
    for (i=0; i<k; i++) {
        gretl_matrix_set(X, s, i, dset.Z[pmod.xlist[i+1]][t]);
    }
    s++;
    }

    err = gretl_matrix_restricted_ols(y, X, rset.R, rset.q,
                      b, S, u, &s2);

    if (!err) {
    if (rset.opt & Gretlopt.OPT_F) {
        err = save_restricted_model(state, pmod, dset, y, b, S,
                    u, s2, rset, prn);
    } else {
        print_restricted_estimates(pmod, dset, b, S, s2,
                       k, rset.g, prn);
    }
    }

 bailout:

    gretl_matrix_block_destroy(B);
    gretl_matrix_free(u);

    return err;
}

/* print or save restricted estimates, single equation */

void print_or_save_result(ExecState * state, GretlRestriction * rset, const(DatasetRaw) * dset, Prn * prn)
{
    Model * pmod = cast(Model *) rset.obj;
    int robust, asym = 0;
    int quiet = (rset.opt & Gretlopt.OPT_Q);
    int silent = (rset.opt & Gretlopt.OPT_S);

    robust = (pmod.opt & Gretlopt.OPT_R);

    if (ASYMPTOTIC_MODEL(cast(GretlCmdIndex) pmod.ci) || rset.rfunc != null) {
    asym = 1;
    } else if (pmod.ci == GretlCmdIndex.panel && (pmod.opt & Gretlopt.OPT_U)) {
    /* random effects */
    asym = 1;
    }

    if (asym) {
    rset.code = TestStats.wald_chisq;
    rset.pval = chisq_cdf_comp(rset.g, rset.test);
    if (!silent) {
        pprintf(prn, "\n%s: %s(%d) = %g, ".ptr, _("Test statistic").ptr,
            (robust)? _("Robust chi^2").ptr: "chi^2".ptr,
            rset.g, rset.test);
    }
    } else {
    rset.code = TestStats.f;
    rset.test /= rset.g;
    rset.pval = snedecor_cdf_comp(rset.g, pmod.dfd, rset.test);
    if (!silent) {
        pprintf(prn, "\n%s: %s(%d, %d) = %g, ".ptr, _("Test statistic").ptr,
            (robust)? _("Robust F").ptr: "F".ptr,
            rset.g, pmod.dfd, rset.test);
    }
    }

    if (!na(rset.pval) && !silent) {
    pprintf(prn, _("with p-value = %g\n").ptr, rset.pval);
    }

    if (!silent) {
    pputc(prn, '\n');
    }

    if (!(rset.opt & Gretlopt.OPT_C)) {
    record_test_result(rset.test, rset.pval);
    }

    if (pmod != null && pmod.ci == GretlCmdIndex.ols &&
    dset != null && dset.Z != null) {
    if ((rset.opt & Gretlopt.OPT_F) || !(quiet || silent)) {
        /* print and/or save restricted estimates */
        do_restricted_estimates(state, rset, dset, prn);
    }
    }
}

/* execute the test, for a single equation */

int test_restriction_set(GretlRestriction * rset, Prn * prn)
{
    Model * pmod = cast(Model *) rset.obj;
    GretlMatrixRaw * vcv = null;
    GretlVector * b = null;
    GretlVector * br = null;
    GretlMatrixRaw * RvR = null;
    int err, freeRvR = 1;

    gretl_error_clear();

    err = restriction_set_form_matrices(rset);
    if (err) {
    return err;
    }

    err = check_R_matrix(rset.R);

    if (!err) {
    b = gretl_coeff_vector_from_model(pmod, rset.mask, &err);
    }

    if (!err) {
    vcv = gretl_vcv_matrix_from_model(pmod, rset.mask, &err);
    }

    if (err) {
    goto bailout;
    }

    if (rset.g == 0) {
    rset.g = gretl_matrix_rows(rset.q);
    }

    br = gretl_column_vector_alloc(rset.g);
    if (br == null) {
    err = GretlError.alloc;
    goto bailout;
    }

    err = gretl_matrix_multiply(rset.R, b, br);
    if (err) {
    fprintf(stderr, "Failed: gretl_matrix_multiply(R, b, br)\n".ptr);
    goto bailout;
    }

    if (rset.opt & Gretlopt.OPT_C) {
    rset.bsum = br.val[0];
    }

    if (!gretl_is_zero_matrix(rset.q)) {
    err = gretl_matrix_subtract_from(br, rset.q);
    if (err) {
        fprintf(stderr, "Failed: gretl_matrix_subtract_from(br, q)\n".ptr);
        goto bailout;
    }
    }

    if (gretl_is_identity_matrix(rset.R)) {
    RvR = vcv;
    freeRvR = 0;
    } else {
    RvR = gretl_matrix_alloc(rset.R.rows, rset.R.rows);
    if (RvR == null) {
        err = GretlError.alloc;
        goto bailout;
    }
    gretl_matrix_qform(rset.R, GretlMatrixMod.none, vcv,
               RvR, GretlMatrixMod.none);
    if (rset.opt & Gretlopt.OPT_C) {
        rset.bse = sqrt(RvR.val[0]);
    }
    }

    err = gretl_invert_symmetric_matrix(RvR);
    if (err) {
    pputs(prn, _("Matrix inversion failed:\n restrictions may be inconsistent or redundant\n").ptr);
    goto bailout;
    }

    rset.test = gretl_scalar_qform(br, RvR, &err);
    if (err) {
    pputs(prn, _("Failed to compute test statistic\n").ptr);
    goto bailout;
    }

 bailout:

    gretl_matrix_free(vcv);
    gretl_vector_free(b);
    gretl_vector_free(br);

    if (freeRvR) {
    gretl_matrix_free(RvR);
    }

    return err;
}

int bootstrap_zero_restriction(GretlRestriction * rset, Model * pmod, const(DatasetRaw) * dset, int method, Prn * prn)
{
    RRow * r = rset.rows[0];
    Gretlopt bopt = Gretlopt.OPT_P | Gretlopt.OPT_R;
    int B = 0;
    int err = 0;

    if (rset.opt & Gretlopt.OPT_S) {
    /* silent */
    bopt |= Gretlopt.OPT_S;
    }

    if (method == BootMethod.pairs) {
    bopt |= Gretlopt.OPT_X;
    } else if (method == BootMethod.wild) {
    bopt |= Gretlopt.OPT_W;
    } else if (method == BootMethod.parametric) {
    bopt |= Gretlopt.OPT_N;
    }

    gretl_restriction_get_boot_params(&B, &bopt);
    err = bootstrap_analysis(pmod, r.bnum[0], B, 0, dset,
                 bopt, prn);
    return err;
}

/* Are we able to do a bootstrap test of the restrictions on
   the given model?  Return non-zero if not. */

int check_bootstrap(Model * pmod, Prn * prn)
{
    int ok = bootstrap_ok(pmod.ci);

    if (!ok) {
    pputs(prn, "Sorry, the bootstrap option is not supported for this test");
    pputc(prn, '\n');
    }

    return (ok)? 0 : GretlError.notimp;
}

/* Do we have a simple b[i] = 0? */

int is_simple_zero_restriction(GretlRestriction * rset)
{
    if (rset.rows != null && rset.g == 1) {
    return (rset.rows[0].nterms == 1 && rset.rows[0].rhs == 0);
    } else {
    return 0;
    }
}

int get_restriction_boot_method(int * method)
{
    const(char) * s = get_optval_string(GretlCmdIndex.restrict, Gretlopt.OPT_B);
    int err = 0;

    if (s != null) {
    if (!strcmp(s, "pairs")) {
        *method = BootMethod.pairs;
    } else if (!strcmp(s, "wild")) {
        *method = BootMethod.wild;
    } else if (!strcmp(s, "parametric")) {
        *method = BootMethod.parametric;
    } else if (strcmp(s, "residuals")) {
        gretl_errmsg_sprintf(_("field '%s' in command is invalid").ptr, s);
        err = GretlError.data;
    }
    }

    return err;
}

int do_single_equation_test(ExecState * state, GretlRestriction * rset, const(DatasetRaw) * dset, Prn * prn)
{
    int err = 0;

    if (rset.rows != null) {
    err = restriction_set_make_mask(rset);
    if (err) {
        return err;
    }
    }

    if (rset.opt & Gretlopt.OPT_B) {
    /* bootstrapping, if possible */
    Model * pmod = cast(Model *) rset.obj;
    int method = 1;

    err = check_bootstrap(pmod, prn);

    if (!err) {
        err = get_restriction_boot_method(&method);
    }

    if (err) {
        return err;
    }

    if (is_simple_zero_restriction(rset)) {
        err = bootstrap_zero_restriction(rset, pmod, dset,
                         method, prn);
    } else {
        err = test_restriction_set(rset, prn);
        if (!err && rset.R.cols != pmod.ncoeff) {
        err = rset_expand_R(rset, pmod.ncoeff);
        }
        if (!err) {
        rset.test /= rset.g;
        err = bootstrap_test_restriction(pmod, rset.R, rset.q,
                         rset.test, rset.g,
                         dset, rset.opt,
                         method, prn);
        }
    }
    } else {
    err = test_restriction_set(rset, prn);
    if (!err) {
        print_or_save_result(state, rset, dset, prn);
    }
    }

    return err;
}

GretlVar * gretl_restricted_vecm(GretlRestriction * rset, const(DatasetRaw) * dset, Gretlopt opt, Prn * prn, int * err)
{
    GretlVar * jvar = null;

    if (rset == null || rset.otype != GretlObjType.var) {
    *err = GretlError.data;
    return null;
    }

    rset.opt |= opt;

    *err = restriction_set_form_matrices(rset);

    if (rset.rows != null && !(opt & Gretlopt.OPT_S)) {
    print_restriction_set(rset, dset, opt, prn);
    }

    if (!*err) {
    jvar = real_gretl_restricted_vecm(cast(GretlVar *) rset.obj, rset, dset,
                      prn, err);
    }

    destroy_restriction_set(rset);

    if (jvar != null && *err != 0) {
    gretl_VAR_free(jvar);
    jvar = null;
    }

    return jvar;
}

enum RCOEFFNAME = "restr__b";

int nonlinear_wald_test(GretlRestriction * rset, Gretlopt opt, Prn * prn)
{
    GretlMatrixRaw * coeff = null;
    GretlMatrixRaw * vcv = null;
    GretlMatrixRaw * J = null;
    GretlMatrixRaw * bread = null;
    GretlMatrixRaw * ham = null;
    char[64] fncall;
    int published = 0;
    int err = 0;

    if (get_user_function_by_name(rset.rfunc) == null) {
    gretl_errmsg_sprintf(_("The symbol '%s' is undefined\n").ptr, rset.rfunc);
    return GretlError.unkvar;
    }

    if (rset.otype == GretlObjType.eqn) {
    Model * pmod = cast(Model *) rset.obj;

    coeff = gretl_coeff_vector_from_model(pmod, null, &err);
    if (!err) {
        vcv = gretl_vcv_matrix_from_model(pmod, null, &err);
    }
    } else if (rset.otype == GretlObjType.sys) {
    EquationSystem * sys = cast(EquationSystem *) rset.obj;

    coeff = equation_system_get_matrix(sys, ModelDataIndex.coeff, &err);
    if (!err) {
        vcv = equation_system_get_matrix(sys, ModelDataIndex.vcv, &err);
    }
    } else {
    return GretlError.notimp; /* relax this later? */
    }

    if (!err) {
    /* "publish" the coeff matrix temporarily */
    err = private_matrix_add(coeff, RCOEFFNAME);
    if (!err) {
        published = 1;
    }
    }

    if (!err) {
    /* formulate function call string and make the call:
       should get a vector */
    sprintf(fncall.ptr, "%s(%s)".ptr, rset.rfunc, RCOEFFNAME.ptr);
    bread = generate_matrix(fncall.ptr, null, &err);
    }

    if (!err) {
    rset.g = gretl_vector_get_length(bread);
    if (rset.g == 0) {
        err = GretlError.nonconf;
    }
    }

    if (!err) {
    J = user_fdjac(coeff, fncall.ptr, 0.0, null, &err);
    }

    if (!err) {
    ham = gretl_matrix_alloc(J.rows, J.rows);
    if (ham == null) {
        err = GretlError.alloc;
    } else {
        err = gretl_matrix_qform(J, GretlMatrixMod.none, vcv,
                     ham, GretlMatrixMod.none);
    }
    }

    if (!err) {
    err = gretl_invpd(ham);
    }

    if (!err) {
    rset.test = gretl_scalar_qform(bread, ham, &err);
    }

    if (!err) {
    print_or_save_result(null, rset, null, prn);
    }

    gretl_matrix_free(vcv);
    gretl_matrix_free(J);
    gretl_matrix_free(bread);
    gretl_matrix_free(ham);

    if (published) {
    destroy_private_matrices();
    } else {
    gretl_matrix_free(coeff);
    }

    return err;
}

/* Respond to "end restrict": in the case of a single equation, go
   ahead and do the test; for a VECM, hand off to the driver in var.c;
   for non-VAR systems, form the restriction matrices R and q and
   attach these to the equation system.
*/

int gretl_restriction_finalize_full(ExecState * state, GretlRestriction * rset, const(DatasetRaw) * dset, Gretlopt opt, Prn * prn)
{
    int silent, err = 0;

    if (rset == null) {
    return 1;
    }

    rset.opt |= opt;
    silent = (rset.opt & Gretlopt.OPT_S);

    if (rset.otype == GretlObjType.eqn && (opt & Gretlopt.OPT_F)) {
    /* --full option is for OLS only */
    Model * pmod = cast(Model *) rset.obj;

    if (pmod.ci != GretlCmdIndex.ols) {
        return GretlError.badopt;
    }
    }

    if (rset.rfunc != null) {
    err = nonlinear_wald_test(rset, opt, prn);
    destroy_restriction_set(rset);
    return err;
    }

    if (rset.otype != GretlObjType.eqn) {
    err = restriction_set_form_matrices(rset);
    if (err) {
        destroy_restriction_set(rset);
        return err;
    }
    }

    if (!silent && rset.rows != null) {
    print_restriction_set(rset, dset, opt, prn);
    }

    if (rset.otype == GretlObjType.var) {
    if (gretl_VECM_rank(cast(GretlVar *) rset.obj) > 0) {
        err = gretl_VECM_test(cast(GretlVar *) rset.obj, rset, dset, rset.opt, prn);
    } else {
        err = gretl_VAR_test(cast(GretlVar *) rset.obj, rset, dset, rset.opt, prn);
    }
    } else if (rset.otype == GretlObjType.sys) {
    EquationSystem * sys = cast(EquationSystem *) rset.obj;

    if (sys.flags & SystemFlags.robust) {
        /* for now, we'll only do a Wald test */
        err = system_wald_test(sys, rset.R, rset.q, opt, prn);
    } else if (rset.opt & Gretlopt.OPT_W) {
        err = system_wald_test(sys, rset.R, rset.q, opt, prn);
    } else {
        system_set_restriction_matrices(sys, rset.R, rset.q);
        rset.R = null;
        rset.q = null;
    }
    } else {
    /* single-equation model */
    err = do_single_equation_test(state, rset, dset, prn);
    }

    if (!(rset.opt & Gretlopt.OPT_C)) {
    destroy_restriction_set(rset);
    }

    return err;
}

int gretl_restriction_finalize(GretlRestriction * rset, const(DatasetRaw) * dset, Gretlopt opt, Prn * prn)
{
    return gretl_restriction_finalize_full(null, rset, dset,
                       opt, prn);
}

/**
 * gretl_sum_test:
 * @list: list of variables to use.
 * @pmod: pointer to model.
 * @dset: information on the data set.
 * @opt: option flags (OPT_Q gives quiet operation).
 * @prn: gretl printing struct.
 *
 * Calculates the sum of the coefficients, relative to the given model,
 * for the variables given in @list.  Prints this estimate along
 * with its standard error, unless OPT_Q is given.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int gretl_sum_test(const(int) * list, Model * pmod, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    GretlRestriction * r;
    char[MAXLEN] line;
    char[36] bstr;
    int quiet = (opt & Gretlopt.OPT_Q);
    int i, len, err = 0;

    if (list[0] < 2) {
    pprintf(prn, _("Invalid input\n").ptr);
    return GretlError.data;
    }

    if (!command_ok_for_model(GretlCmdIndex.coeffsum, Gretlopt.OPT_NONE, pmod)) {
    return GretlError.notimp;
    }

    if (exact_fit_check(pmod, prn)) {
    return 0;
    }

    r = restriction_set_new(pmod, GretlObjType.eqn, Gretlopt.OPT_Q | Gretlopt.OPT_C);
    if (r == null) {
    return 1;
    }

    *line.ptr = '\0';
    len = 0;

    for (i=1; i<=list[0]; i++) {
    sprintf(bstr.ptr, "b[%s]".ptr, dset.varname[list[i]]);
    len += strlen(bstr.ptr) + 4;
    if (len >= MAXLEN - 1) {
        err = GretlError.parse;
        break;
    }
    strcat(line.ptr, bstr.ptr);
    if (i < list[0]) {
        strcat(line.ptr, " + ");
    } else {
        strcat(line.ptr, " = 0");
    }
    }

    if (!err) {
    err = real_restriction_set_parse_line(&r, line.ptr, dset, 1);
    }

    if (!err) {
    err = gretl_restriction_finalize(r, dset, Gretlopt.OPT_NONE, null);
    }

    if (!err) {
    double test;

    if (!quiet) {
        pprintf(prn, "\n%s: ".ptr, _("Variables").ptr);
        for (i=1; i<=list[0]; i++) {
        pprintf(prn, "%s ".ptr, dset.varname[list[i]]);
        }
        pprintf(prn, "\n   %s = %g\n".ptr, _("Sum of coefficients").ptr, r.bsum);
    }

    test = sqrt(r.test);
    if (r.bsum < 0) {
        test = -test;
    }

    if (r.code == TestStats.f) {
        if (!quiet) {
        pprintf(prn, "   %s = %g\n".ptr, _("Standard error").ptr, r.bse);
        pprintf(prn, "   t(%d) = %g ".ptr, pmod.dfd, test);
        pprintf(prn, _("with p-value = %g\n").ptr, r.pval);
        }
        record_test_result(test, r.pval);
    } else if (r.code == TestStats.wald_chisq) {
        if (!quiet) {
        pprintf(prn, "   %s = %g\n".ptr, _("Standard error").ptr, r.bse);
        r.pval = normal_pvalue_2(test);
        pprintf(prn, "   z = %g ".ptr, test);
        pprintf(prn, _("with p-value = %g\n").ptr, r.pval);
        }
        record_test_result(test, r.pval);
    }

    destroy_restriction_set(r);
    }

    return err;
}

__gshared int restrict_B = 0;
__gshared Gretlopt rboot_opt = Gretlopt.OPT_NONE;

int gretl_restriction_set_boot_params(int B, Gretlopt opt)
{
    int err = 0;

    rboot_opt = opt;

    if (B > 0) {
    restrict_B = B;
    } else {
    err = GretlError.data;
    }

    return err;
}

void gretl_restriction_get_boot_params(int * pB, Gretlopt * popt)
{
    *pB = restrict_B;
    *popt |= rboot_opt;

    /* these are ad hoc values */
    restrict_B = 0;
    rboot_opt = Gretlopt.OPT_NONE;
}

Gretlopt gretl_restriction_get_options(const(GretlRestriction) * rset)
{
    return (rset != null)? rset.opt : Gretlopt.OPT_NONE;
}

GretlObjType gretl_restriction_get_type(const(GretlRestriction) * rset)
{
    return (rset != null)? rset.otype : GretlObjType.null_;
}

const(GretlMatrixRaw) * rset_get_R_matrix(const(GretlRestriction) * rset)
{
    return (rset != null)? rset.R : null;
}

const(GretlMatrixRaw) * rset_get_q_matrix(const(GretlRestriction) * rset)
{
    return (rset != null)? rset.q : null;
}

const(GretlMatrixRaw) * rset_get_Ra_matrix(const(GretlRestriction) * rset)
{
    return (rset != null)? rset.Ra : null;
}

const(GretlMatrixRaw) * rset_get_qa_matrix(const(GretlRestriction) * rset)
{
    return (rset != null)? rset.qa : null;
}

int rset_VECM_bcols(const(GretlRestriction) * rset)
{
    return (rset == null || rset.R == null)? 0 : rset.R.cols;
}

int rset_VECM_acols(const(GretlRestriction) * rset)
{
    return (rset == null || rset.Ra == null)? 0 : rset.Ra.cols;
}

void rset_add_results(GretlRestriction * rset, double test, double pval, double lnl)
{
    rset.test = test;
    rset.pval = pval;
    rset.lnl = lnl;
}

void rset_record_LR_result(GretlRestriction * rset)
{
    record_LR_test_result(rset.test, rset.pval, rset.lnl);
}

