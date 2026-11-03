/* Port of lib/src/system.c: systems of simultaneous equations (SUR, 3SLS, FIML, LIML, OLS, TSLS, WLS): defining a
   system from "system"/"equation"/"identity" lines, list handling and variable categorisation, estimation,
   printing, restrictions, forecasts, saving to session and XML. The types of system.h are systemh.d. */

module gretl.system;


import core.stdc.ctype : isspace;
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



import gretl.systemh : EquationSystem, GretlSystemMethods, IdAtom, Identity, LimlData, Predet, SystemFlags;

import core.stdc.config : c_long;

/* the C library's strtol takes a plain char ** end pointer */
extern (C) c_long strtol(const(char) * nptr, char ** endptr, int base);
import gretl.compareh : ModelAuxCode;
import gretl.dataset : series_get_lag, series_get_parent_name, series_get_stack_level, series_is_listarg;
import gretl.describe : gretl_mean, gretl_stddev, ljung_box, multivariate_normality_test;
import gretl.estimate : array_arch_test;
import gretl.genmain : current_series_index, extract_varname, generate_list, series_index;
import gretl.genmainh : ModelDataIndex;
import gretl.gmem : g_free;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_errors : gretl_errmsg_set, gretl_errmsg_sprintf, gretl_error_clear;
import gretl.gretl_list : gretl_list_append_term, gretl_list_copy, gretl_list_delete_at_pos, gretl_list_duplicates, gretl_list_has_separator, gretl_list_new, gretl_list_separator_position, gretl_list_split_on_separator, gretl_list_sublist, gretl_null_list, in_gretl_list, reglist_check_for_const;
import gretl.gretl_listh : LISTSEP;
import gretl.gretl_model : clear_model, gretl_model_destroy_data_item, gretl_model_get_double, gretl_model_get_int, gretl_model_init, gretl_model_set_double, gretl_model_set_int, set_model_id;
import gretl.gretl_prn : Prn, gretl_prn_newline, pprintf, pputc, pputs, tex_format;
import gretl.gretl_utils : record_test_result;
import gretl.gretl_xml : gretl_matrix_serialize, gretl_xml_get_list, gretl_xml_get_matrix, gretl_xml_get_prop_as_int, gretl_xml_get_prop_as_string, gretl_xml_put_tagged_list;
import gretl.gstrfuncs : g_strdup, g_strdup_printf;
import gretl.interacth : MAXSAVENAME;
import gretl.missing : list_adjust_sample;
import gretl.modelprint : printmodel;
import gretl.objstackh : GretlObjType, SavedObjectFlags;
import gretl.plugins : get_plugin_function;
import gretl.printout : print_contemp_covariance_matrix;
import gretl.pvalues : chisq_cdf_comp, snedecor_cdf_comp;
import gretl.strutils : gretl_lower, gretl_strdup, gretl_strndup, gretl_unquote, gretl_xml_encode;
import gretl.transforms : laggenr;
import gretl.tsls : ivreg_process_lists;
import gretl.gretl_array : get_array_by_name, gretl_array_get_element, gretl_array_get_length, gretl_array_get_type;
import gretl.gretl_bundle : gretl_bundle_donate_data, gretl_bundle_new, gretl_bundle_set_int, gretl_bundle_set_list, gretl_bundle_set_matrix, gretl_bundle_set_scalar, gretl_bundle_set_string;
import gretl.uservar : get_list_by_name;
import gretl.usermat : get_matrix_by_name;

/* unported dependencies: per-file declarations from the 2026c headers */
extern (C) {
    EquationSystem * get_equation_system_by_name(const(char) * sname);
    void * get_last_model(GretlObjType * type);
    int gretl_function_depth();
    void gretl_object_ref(void * ptr, GretlObjType type);
    void gretl_object_unref(void * ptr, GretlObjType type);
    int gretl_stack_object_as(void * ptr, GretlObjType type, const(char) * name);
    void set_as_last_model(void * ptr, GretlObjType type);
}

enum IdentityOp {
    plus,
    minus
}

enum AuxListType {
    endog,
    instr
}

enum SysTestType {
    none,
    lr,
    f,
    notimp
}

enum SysSeries {
    uhat,
    yhat
}

__gshared const(char) *[8] gretl_system_method_strings = ["sur".ptr, "3sls".ptr, "fiml".ptr, "liml".ptr, "ols".ptr, "tsls".ptr, "wls".ptr, null];

__gshared const(char) *[8] gretl_system_short_strings = ["SUR".ptr, "3SLS".ptr, "FIML".ptr, "LIML".ptr, "OLS".ptr, "TSLS".ptr, "WLS".ptr, null];

__gshared const(char) *[8] gretl_system_long_strings = ["Seemingly Unrelated Regressions".ptr, "Three-Stage Least Squares".ptr, "Full Information Maximum Likelihood".ptr, "Limited Information Maximum Likelihood".ptr, "Ordinary Least Squares".ptr, "Two-Stage Least Squares".ptr, "Weighted Least Squares".ptr, null];

const(char) * nosystem = "No system of equations has been defined".ptr;
const(char) * badsystem = "Unrecognized equation system type".ptr;
const(char) * toofew = "An equation system must have at least two equations".ptr;

bool sys_anonymous(const(char) * s)
{
    return (strcmp(s, "$system") == 0);
}

/* the C keeps the anonymous systems in a GList; a D array does the same job */
__gshared EquationSystem *[] sysstack;

void sysstack_remove(EquationSystem * sys)
{
    foreach (i, e; sysstack) {
        if (e is sys) {
            sysstack = sysstack[0 .. i] ~ sysstack[i+1 .. $];
            return;
        }
    }
}

EquationSystem * get_anon_system_at_depth(int fd)
{
    foreach (sys; sysstack) {
        if (sys.fd == fd) {
            return sys;
        }
    }

    return null;
}

EquationSystem * get_anonymous_equation_system()
{
    int fd = gretl_function_depth();

    return get_anon_system_at_depth(fd);
}

/* The use-case for push_anon_system() is when a user wants to define
   a system and then estimate it by more than one method (or set it as
   a target for "restrict" before estimation) but does not care to
   give it a specific name via the "name <- system" mechanism -- or
   we're inside a user-defined function where this naming mechanism is
   not available.

   In this case we save the system "anonymously" instead (under a name
   of "$system").
*/

void push_anon_system(EquationSystem * sys)
{
    EquationSystem * old = get_anon_system_at_depth(sys.fd);

    if (old != null) {
        sysstack_remove(old);
        gretl_object_unref(old, GretlObjType.sys);
    }

    gretl_object_ref(sys, GretlObjType.sys);
    sysstack ~= sys;
}

/* For use when terminating function execution: if an
   anonymous equation system was set up at the given
   level of function execution, trash it.
*/

void delete_anonymous_equation_system(int level)
{
    EquationSystem * sys = get_anon_system_at_depth(level);

    if (sys != null) {
        sysstack_remove(sys);
        gretl_object_unref(sys, GretlObjType.sys);
    }
}

void print_system_equation(const(int) * list, const(DatasetRaw) * dset, Prn * prn)
{
    int i, v;

    pputs(prn, "equation");

    for (i=1; i<=list[0]; i++) {
        v = list[i];
        if (v == LISTSEP) {
            pputs(prn, " ;");
        } else if (v > 0 && v < dset.v) {
            pprintf(prn, " %s".ptr, dset.varname[v]);
        } else {
            pprintf(prn, " %d".ptr, v);
        }
    }

    pputc(prn, '\n');
}

void print_system_identity(const(Identity) * ident, const(DatasetRaw) * dset, Gretlopt opt, Prn * prn)
{
    int i;

    if (opt & Gretlopt.OPT_H) {
        pprintf(prn, "Identity: %s = %s ".ptr,
                dset.varname[ident.depvar],
                dset.varname[ident.atoms[0].varnum]);
    } else {
        pprintf(prn, "identity %s = %s ".ptr,
                dset.varname[ident.depvar],
                dset.varname[ident.atoms[0].varnum]);
    }

    for (i=1; i<ident.n_atoms; i++) {
        pprintf(prn, "%c %s ".ptr, (ident.atoms[i].op == IdentityOp.plus)? '+' : '-',
                dset.varname[ident.atoms[i].varnum]);
    }

    pputc(prn, '\n');
}

int sys_max_predet_lag(EquationSystem * sys)
{
    int i, m = 0;

    for (i=0; i<sys.plist[0]; i++) {
        if (sys.pre_vars[i].lag > m) {
            m = sys.pre_vars[i].lag;
        }
    }

    return m;
}

int get_predet_parent(EquationSystem * sys, int v, int * lag)
{
    int i;

    for (i=0; i<sys.plist[0]; i++) {
        if (v == sys.pre_vars[i].id) {
            *lag = sys.pre_vars[i].lag;
            return sys.pre_vars[i].src;
        }
    }

    *lag = 0;

    return -1;
}

/**
 * print_equation_system_info:
 * @sys: gretl equation system.
 * @dset: dataset information.
 * @opt: use OPT_H for printing in a form designed to appear
 * in the header when results of estimation are printed.
 * @prn: printing struct.
 *
 * Prints details of an equation system to @prn.
 */

void print_equation_system_info(EquationSystem * sys, const(DatasetRaw) * dset, Gretlopt opt, Prn * prn)
{
    int header = (opt & Gretlopt.OPT_H);
    int i, vi, lag;

    if (header && sys.name != null && !sys_anonymous(sys.name)) {
        pprintf(prn, "%s %s\n".ptr, _("Equation system").ptr, sys.name);
    }

    if (!header) {
        for (i=0; i<sys.neqns; i++) {
            print_system_equation(sys.lists[i], dset, prn);
        }
    }

    for (i=0; i<sys.nidents; i++) {
        print_system_identity(sys.idents[i], dset, opt, prn);
    }

    if (sys.ylist != null) {
        pputs(prn, (header)? _("Endogenous variables:").ptr : "endog");
        for (i=1; i<=sys.ylist[0]; i++) {
            vi = sys.ylist[i];
            pprintf(prn, " %s".ptr, dset.varname[vi]);
        }
        pputc(prn, '\n');
    }

    if (header) {
        if (sys.pre_vars != null) {
            pputs(prn, _("Predetermined variables:").ptr);
            for (i=0; i<sys.plist[0]; i++) {
                vi = sys.pre_vars[i].src;
                lag = sys.pre_vars[i].lag;
                pprintf(prn, " %s(-%d)".ptr, dset.varname[vi], lag);
            }
            pputc(prn, '\n');
        }
        if (sys.ilist != null && sys.ilist[0] > sys.plist[0]) {
            pputs(prn, _("Exogenous variables:").ptr);
            for (i=1; i<=sys.ilist[0]; i++) {
                vi = sys.ilist[i];
                if (!in_gretl_list(sys.plist, vi)) {
                    pprintf(prn, " %s".ptr, dset.varname[vi]);
                }
            }
            pputc(prn, '\n');
        }
    } else if (sys.ilist != null) {
        pputs(prn, "instr");
        for (i=1; i<=sys.ilist[0]; i++) {
            vi = sys.ilist[i];
            pprintf(prn, " %s".ptr, dset.varname[vi]);
        }
        pputc(prn, '\n');
    }
}

int system_method_from_string(const(char) * s)
{
    int i = 0;

    while (gretl_system_method_strings[i] != null) {
        if (!strcmp(s, gretl_system_method_strings[i])) {
            return i;
        }
        i++;
    }

    return i;
}

const(char) * system_method_full_string(int method)
{
    if (method >= GretlSystemMethods.sur && method < GretlSystemMethods.max) {
        return gretl_system_long_strings[method];
    } else {
        return null;
    }
}

const(char) * system_method_short_string(int method)
{
    if (method >= GretlSystemMethods.sur && method < GretlSystemMethods.max) {
        return gretl_system_method_strings[method];
    } else {
        return null;
    }
}

EquationSystem * equation_system_new(int method, const(char) * name, int * err)
{
    EquationSystem * sys;

    if (method < 0 && name == null) {
        *err = GretlError.data;
        return null;
    }

    sys = cast(typeof(sys)) malloc((*sys).sizeof);
    if (sys == null) {
        *err = GretlError.alloc;
        return null;
    }

    sys.name = null;

    if (name != null) {
        equation_system_set_name(sys, name);
    }

    sys.refcount = 0;
    sys.fd = gretl_function_depth();
    sys.method = method;

    sys.t1 = sys.t2 = 0;
    sys.smpl_t1 = sys.smpl_t2 = 0;
    sys.T = sys.df = 0;

    sys.neqns = 0;
    sys.nidents = 0;

    sys.q = sys.R = null;

    sys.order = sys.flags = sys.iters = 0;

    sys.ll = sys.llu = sys.ldet = NADBL;
    sys.X2 = NADBL;
    sys.ess = NADBL;
    sys.diag_test = 0.0;
    sys.bdiff = 0.0;

    sys.b = sys.vcv = null;
    sys.E = sys.S = null;
    sys.yhat = null;

    sys.Gamma = null;
    sys.A = sys.B = null;
    sys.F = sys.Sr = null;

    sys.tslists = sys.lists = null;
    sys.ilist = sys.ylist = null;
    sys.plist = sys.xlist = null;
    sys.biglist = null;

    sys.pre_vars = null;
    sys.idents = null;

    sys.models = null;
    sys.ldata = null;

    return sys;
}

void system_clear_restrictions(EquationSystem * sys)
{
    if (sys.R != null) {
        free(sys.R);
        sys.R = null;
    }

    if (sys.q != null) {
        free(sys.q);
        sys.q = null;
    }

    sys.flags &= ~SystemFlags.restrict;
}

void system_clear_results(EquationSystem * sys)
{
    sys.iters = 0;

    sys.t1 = sys.t2 = 0;
    sys.smpl_t1 = sys.smpl_t2 = 0;
    sys.T = sys.df = 0;

    sys.ll = NADBL;
    sys.llu = NADBL;
    sys.ldet = NADBL;
    sys.X2 = NADBL;
    sys.ess = NADBL;

    sys.diag_test = 0.0;
    sys.bdiff = 0.0;

    if (sys.b != null) {
        gretl_matrix_replace(&sys.b, null);
    }

    if (sys.vcv != null) {
        gretl_matrix_replace(&sys.vcv, null);
    }

    if (sys.S != null) {
        gretl_matrix_replace(&sys.S, null);
    }

    if (sys.E != null) {
        gretl_matrix_replace(&sys.E, null);
    }

    if (sys.yhat != null) {
        gretl_matrix_replace(&sys.yhat, null);
    }

    if (sys.Gamma != null) {
        gretl_matrix_replace(&sys.Gamma, null);
    }

    if (sys.B != null) {
        gretl_matrix_replace(&sys.B, null);
    }

    if (sys.A != null) {
        gretl_matrix_replace(&sys.A, null);
    }

    if (sys.Sr != null) {
        gretl_matrix_replace(&sys.Sr, null);
    }

    if (sys.F != null) {
        gretl_matrix_replace(&sys.F, null);
    }

    if (sys.ldata != null) {
        free(sys.ldata.lmin);
        free(sys.ldata.ll);
        free(sys.ldata.idf);
        free(sys.ldata);
        sys.ldata = null;
    }
}

void equation_system_destroy(EquationSystem * sys)
{
    int i;

    if (sys == null || sys.lists == null) {
        return;
    }

    sys.refcount -= 1;
    if (sys.refcount > 0) {
        return;
    }

    for (i=0; i<sys.neqns; i++) {
        free(sys.lists[i]);
    }

    free(sys.lists);
    sys.lists = null;

    for (i=0; i<sys.nidents; i++) {
        destroy_ident(sys.idents[i]);
    }

    free(sys.idents);

    free(sys.ylist);
    free(sys.ilist);
    free(sys.xlist);
    free(sys.plist);
    free(sys.biglist);

    free(sys.pre_vars);

    free(sys.name);

    gretl_matrix_free(sys.R);
    gretl_matrix_free(sys.q);

    system_clear_results(sys);

    free(sys);
}

/* 2020-12-09: When a system is defined using TSLS-style
   regression lists, with instruments specified per
   equation, there were problems estimating the system
   via estimators that don't take such a specification.
   The apparatus below is designed to handle this: we
   remove the instruments immediately prior to estimation,
   then restore them immediately afterwards.
*/

/* Make a copy of the incoming per-equation regression
   lists, with any instruments removed.
*/

int lists_strip_insts(EquationSystem * sys)
{
    int i;
    int * eqlist;

    sys.tslists = cast(typeof(sys.tslists)) malloc(sys.neqns * (*sys.tslists).sizeof);

    for (i=0; i<sys.neqns; i++) {
    sys.tslists[i] = sys.lists[i];
    if (gretl_list_has_separator(sys.lists[i])) {
        gretl_list_split_on_separator(sys.tslists[i], &eqlist, null);
        sys.lists[i] = eqlist;
    } else {
        sys.lists[i] = gretl_list_copy(sys.tslists[i]);
    }
    }

    return 0;
}

/* Restore the original per-equation regression lists,
   in case instruments were removed for estimation by
   OLS, WLS, SUR or FIML.
*/

void lists_restore_insts(EquationSystem * sys)
{
    int i;

    for (i=0; i<sys.neqns; i++) {
    free(sys.lists[i]);
    sys.lists[i] = sys.tslists[i];
    sys.tslists[i] = null;
    }

    free(sys.tslists);
    sys.tslists = null;
}

int sys_rearrange_eqn_lists(EquationSystem * sys, const(DatasetRaw) * dset)
{
    int i, err = 0;

    for (i=0; i<sys.neqns; i++) {
        reglist_check_for_const(sys.lists[i], dset);
    }

    if (sys.method != GretlSystemMethods.tsls &&
        sys.method != GretlSystemMethods.sls3 &&
    sys.method != GretlSystemMethods.liml) {
        /* we can't have ';' in equation lists */
    for (i=0; i<sys.neqns; i++) {
        if (gretl_list_has_separator(sys.lists[i])) {
        lists_strip_insts(sys);
        break;
        }
    }
    }

    return err;
}

/* Form a list from row @i of matrix @m, ignoring trailing
   zero elements. We check that series ID numbers are in
   bounds, and also that the list does not contain
   duplicated elements.
*/

int * matrix_row_to_list(const(GretlMatrixRaw) * m, int i, const(DatasetRaw) * dset, int * err)
{
    int * list = null;
    int j, k, n = m.cols;

    /* figure how many elements to read */
    for (j=m.cols-1; j>=0; j--) {
        if (gretl_matrix_get(m, i, j) == 0) {
            n--;
        } else {
            break;
        }
    }

    if (n == 0) {
        *err = GretlError.data;
        return null;
    }

    /* check for out-of-bounds values */
    for (j=0; j<n; j++) {
        k = cast(int) gretl_matrix_get(m, i, j);
        if (k < 0 || k >= dset.v) {
            *err = GretlError.unkvar;
            return null;
        }
    }

    /* construct the list */
    list = gretl_list_new(n);
    if (list == null) {
        *err = GretlError.alloc;
    } else {
        k = 1;
        for (j=0; j<n; j++) {
            list[k++] = cast(int) gretl_matrix_get(m, i, j);
        }
    }

    if (!*err) {
        k = gretl_list_duplicates(list, GretlCmdIndex.equation);
        if (k >= 0) {
            gretl_errmsg_sprintf(_("variable %d duplicated in the command list.").ptr, k);
            *err = GretlError.data;
            free(list);
            list = null;
        }
    }

    return list;
}

int sys_check_eqn_list(const(int) * list)
{
    int dupv = gretl_list_duplicates(list, GretlCmdIndex.equation);

    if (dupv >= 0) {
        gretl_errmsg_sprintf(_("variable %d duplicated in the command list.").ptr, dupv);
        return GretlError.data;
    } else {
        return 0;
    }
}

/* @LY is a simple list of g regressands; either @LX is
   simple list of common regressors or @AX is an array of
   lists. We get either @LX or @AX as non-NULL, never both.
*/

int add_equations_from_lists(EquationSystem * sys, const(int) * LY, const(int) * LX, GretlArray * AX, const(DatasetRaw) * dset)
{
    int nx, n = sys.neqns;
    int g = LY[0];
    int i, j, err = 0;

    if (AX != null && gretl_array_get_length(AX) != g) {
        err = GretlError.args;
    } else {
        int ** L = cast(int **) realloc(sys.lists, (n + g) * (*sys.lists).sizeof);

        if (L == null) {
            err = GretlError.alloc;
        } else {
            sys.lists = L;
            for (i=0; i<g; i++) {
                sys.lists[n+i] = null;
            }
        }
    }

    if (err) {
        return err;
    }

    if (AX != null) {
        /* handle an array of regressor lists */
        int * rhs, list;

        for (i=0; i<g && !err; i++) {
            list = null;
            rhs = cast(int *) gretl_array_get_element(AX, i, null, &err);
            if (!err) {
                nx = rhs[0];
                list = gretl_list_new(1 + nx);
                if (list == null) {
                    err = GretlError.alloc;
                }
            }
            if (!err) {
                list[1] = LY[i+1];
                for (j=1; j<=nx; j++) {
                    list[j+1] = rhs[j];
                }
                err = sys_check_eqn_list(list);
            }
            if (err) {
                free(list);
            } else {
                sys.lists[n++] = list;
            }
        }
    } else {
        /* handle a list of common regressors */
        int nxl = LX[0];
        int * list;
        int * l0 = null;

        for (i=0; i<g && !err; i++) {
            if (i == 0) {
                l0 = list = gretl_list_new(1 + nxl);
            } else {
                list = gretl_list_copy(l0);
            }
            if (list == null) {
                err = GretlError.alloc;
            } else {
                list[1] = LY[i+1];
                if (i == 0) {
                    for (j=1; j<=nxl; j++) {
                        list[j+1] = LX[j];
                    }
                }
                err = sys_check_eqn_list(list);
            }
            if (err) {
                free(list);
            } else {
                sys.lists[n++] = list;
            }
        }
    }

    if (!err) {
        sys.neqns += g;
    }

    return err;
}

int add_equations_from_matrix(EquationSystem * sys, const(GretlMatrixRaw) * m, const(DatasetRaw) * dset)
{
    int * list;
    int ** L;
    int n = sys.neqns;
    int i, err = 0;

    L = cast(typeof(L)) realloc(sys.lists, (n + m.rows) * (*sys.lists).sizeof);
    if (L == null) {
        return GretlError.alloc;
    }

    sys.lists = L;

    for (i=0; i<m.rows && !err; i++) {
        list = matrix_row_to_list(m, i, dset, &err);
        if (!err) {
            err = sys_check_eqn_list(list);
        }
        if (!err) {
            sys.lists[n++] = list;
        }
    }

    if (!err) {
        sys.neqns += m.rows;
    }

    return err;
}

/**
 * equation_system_append_multi:
 * @sys: initialized equation system.
 * @parm1: the name of a pre-defined matrix or list.
 * @parm2: the name of list or array of lists (or %NULL).
 * @dset: pointer to dataset struct.
 *
 * Adds one or more equations to @sys as follows.
 *
 * If @parm2 is %NULL then @parm1 is taken to be the
 * name of a matrix, and we interpret the rows of the
 * specified matrix as lists. Lists of differing length
 * can be accommodated by padding unused trailing elements
 * of short rows with zeros.
 *
 * If @parm2 is non-%NULL then @parm1 is taken to be the name
 * of a list of regressands, one per equation; @parm2,
 * pertaining to regressors, should be the name of a list
 * if the regressors are in common across the equations
 * or an array of lists (one per equation) otherwise.
 *
 * Returns: 0 on success, non-zero on failure, in which case
 * @sys is destroyed.
 */

int equation_system_append_multi(EquationSystem * sys, const(char) * parm1, const(char) * parm2, const(DatasetRaw) * dset)
{
    int err = 0;

    if (sys == null) {
        gretl_errmsg_set(_(nosystem));
        return GretlError.data;
    } else if (parm1 == null) {
        return GretlError.args;
    }

    if (parm2 != null) {
        /* look for two lists, or list plus array */
        const(int) * LY = get_list_by_name(parm1);
        const(int) * LX = get_list_by_name(parm2);
        GretlArray * AX = null;

        if (LY == null) {
            gretl_errmsg_sprintf(_("'%s': no such list").ptr, parm1);
            err = GretlError.unkvar;
        } else if (LX == null) {
            /* try for an array on the right */
            AX = get_array_by_name(parm2);
            if (gretl_array_get_type(AX) != GretlType.lists) {
                gretl_errmsg_sprintf(_("'%s': not a list or array of lists").ptr,
                                     parm2);
                err = GretlError.unkvar;
            }
        }
        if (!err) {
            err = add_equations_from_lists(sys, LY, LX, AX, dset);
        }
    } else {
        /* look for one matrix */
        const(GretlMatrixRaw) * m = get_matrix_by_name(parm1);

        if (m == null) {
            gretl_errmsg_sprintf(_("'%s': no such matrix").ptr, parm1);
            err = GretlError.unkvar;
        } else if (m.rows == 0 || m.cols == 0) {
            err = GretlError.data;
        } else {
            err = add_equations_from_matrix(sys, m, dset);
        }
    }

    if (err) {
        equation_system_destroy(sys);
    }

    return err;
}

/**
 * equation_system_append:
 * @sys: initialized equation system.
 * @list: list containing dependent variable and regressors.
 *
 * Adds an equation, as represented by @list, to @sys.
 *
 * Returns: 0 on success, non-zero on failure, in which case
 * @sys is destroyed.
 */

int equation_system_append(EquationSystem * sys, const(int) * list)
{
    int err = 0;

    if (sys == null) {
        gretl_errmsg_set(_(nosystem));
        err = GretlError.data;
    } else {
        int n = sys.neqns;
        int ** lists;

        lists = cast(typeof(lists)) realloc(sys.lists, (n + 1) * (*sys.lists).sizeof);
        if (lists == null) {
            err = GretlError.alloc;
        } else {
            sys.lists = lists;
            sys.lists[n] = gretl_list_copy(list);
            if (sys.lists[n] == null) {
                err = GretlError.alloc;
            } else {
                sys.neqns += 1;
            }
        }

        if (err) {
            equation_system_destroy(sys);
        }
    }

    return err;
}

/* retrieve the name -- possibly quoted with embedded spaces -- for
   an equation system */

char * get_system_name_from_line(const(char) * s)
{
    const(char) * p = null;
    char * name = null;
    int pchars = 0;

    if (!strncmp(s, "method", 6)) {
        /* skip "method = whatever", with possible
           spaces around '='
        */
        char c = *(s+6);

        if (c == ' ' || c == '=') {
            p = s + 6;
            p += cast(int) strspn(s, " ");
            if (*p == '=') p++;
            p += cast(int) strspn(s, " ");
            p += cast(int) strcspn(s, " ");
            p += cast(int) strspn(s, " ");
            s = p;
        }
    }

    if (*s == '"') {
        if (*(s + 1) != '\0') s++;
        p = s;
        while (*p && *p != '"') {
            if (!isspace(cast(ubyte) *p)) pchars++;
            p++;
        }
        if (*p != '"') {
            /* no closing quote */
            pchars = 0;
        }
    } else {
        p = s;
        while (*p && !isspace(cast(ubyte) *p)) {
            pchars++;
            p++;
        }
    }

    if (pchars > 0) {
        name = gretl_strndup(s, p - s);
    }

    return name;
}

/* parse a system estimation method out of parameter @s */

int sys_get_estimator(const(char) * s)
{
    const(char) * p = strstr(s, "method");
    int method = -1;

    if (p != null) {
        p += 6;
        p += cast(int) strspn(s, " ");
        if (*p == '=') {
            char[9] mstr;

            p++;
            p += cast(int) strspn(s, " ");
            if (sscanf(p, "%8s".ptr, mstr.ptr) == 1) {
                gretl_lower(mstr.ptr);
                method = system_method_from_string(mstr.ptr);
            }
        }
    }

    return method;
}

int parse_sys_start_param(const(char) * s, int * method, char ** name)
{
    int len = cast(int) (cast(int) strlen(s));
    int err = 0;

    if (len > 7 && !strncmp(s, "method=", 7)) {
        char[9] mstr;

        *mstr.ptr = '\0';
        strncat(mstr.ptr, s + 7, 8);
        gretl_lower(mstr.ptr);
        *method = system_method_from_string(mstr.ptr);
        if (*method == GretlSystemMethods.max) {
            /* invalid method was given */
            gretl_errmsg_set(_(badsystem));
            err = GretlError.data;
        }
    } else if (len > 5 && !strncmp(s, "name=", 5)) {
        *name = gretl_strdup(s + 5);
        gretl_unquote(*name, &err);
    } else {
        err = GretlError.parse;
    }

    return err;
}

/**
 * equation_system_start:
 * @param: may incude an estimation method spec.
 * @name: the name to be given to system, if any, otherwise NULL
 * or an empty string.
 * @opt: may include OPT_I for iterative estimation (will be
 * ignored if the the estimation method does not support it).
 * @err: location to receive error code.
 *
 * Start compiling an equation system. Either @param must contain
 * an estimation method or the system must be given a name.
 * If a method is given, the system will be estimated as soon as its
 * definition is complete. If a name is given, the system definition
 * is saved on a stack and it can subsequently be estimated via various
 * methods. If both a name and an estimation method are given, the system
 * is both estimated and saved.
 *
 * Returns: pointer to a new equation system, or %NULL on error.
 */

EquationSystem * equation_system_start(const(char) * param, const(char) * name, Gretlopt opt, int * err)
{
    EquationSystem * sys = null;
    char * oldname = null;
    int method = -1;

    if (param != null) {
        *err = parse_sys_start_param(param, &method, &oldname);
        if (*err) {
            return null;
        }
    }

    if (name != null && *name != '\0') {
        sys = equation_system_new(method, name, err);
    } else if (oldname != null) {
        sys = equation_system_new(method, oldname, err);
        free(oldname);
    } else {
        /* treat the system as "anonymous" */
        sys = equation_system_new(method, "$system", err);
        if (!*err) {
            push_anon_system(sys);
        }
    }

    if (sys != null) {
        if (opt & Gretlopt.OPT_I) {
            sys.flags |= SystemFlags.iterate;
        }
        if (opt & Gretlopt.OPT_Q) {
            sys.flags |= SystemFlags.quiet;
        }
    }

    return sys;
}

/* determine the degrees of freedom for the unrestricted
   estimation of @sys
*/

int system_get_dfu(EquationSystem * sys)
{
    int dfu = sys.T * sys.neqns; /* total observations */
    int i, pos;

    for (i=0; i<sys.neqns; i++) {
        /* subtract the number of parameters */
        pos = gretl_list_separator_position(sys.lists[i]);
        if (pos > 0) {
            dfu -= pos - 2;
        } else {
            dfu -= sys.lists[i][0] - 1;
        }
    }

    return dfu;
}

int get_eqn_ref(EquationSystem * sys, int j)
{
    int i, pos, nparm = 0;

    for (i=0; i<sys.neqns; i++) {
        pos = gretl_list_separator_position(sys.lists[i]);
        if (pos > 0) {
            nparm += pos - 2;
        } else {
            nparm += sys.lists[i][0] - 1;
        }
        if (nparm > j) {
            return i;
        }
    }

    return -1;
}

/* Determine whether the restrictions on @sys embodied in @R
   reference just one equation, and if so return the degrees
   of freedom for that equation. Otherwise return 0.
*/

int maybe_get_single_equation_dfu(EquationSystem * sys, const(GretlMatrixRaw) * R)
{
    int i, j, eq = -1, eq0 = -1;
    int df = 0;

    for (j=0; j<R.cols; j++) {
        for (i=0; i<R.rows; i++) {
            if (gretl_matrix_get(R, i, j) != 0) {
                eq = get_eqn_ref(sys, j);
                if (eq0 < 0) {
                    eq0 = eq;
                } else if (eq != eq0) {
                    /* multiple equations are referenced */
                    return 0;
                }
            }
        }
    }

    if (eq >= 0) {
        df = sys.T - (sys.lists[eq][0] - 1);
    }

    return df;
}

/* Asymptotic F-test, as in William Greene:

   (1/J) * (Rb-q)' * [R*Var(b)*R']^{-1} * (Rb-q)

   or Chi-square version if dfu = 0
*/

int multi_eqn_wald_test(const(GretlMatrixRaw) * b, const(GretlMatrixRaw) * V, const(GretlMatrixRaw) * R, const(GretlMatrixRaw) * q, int dfu, Gretlopt opt, Prn * prn)
{
    GretlMatrixRaw * Rbq, RvR;
    int Rrows, dfn;
    double test = NADBL;
    int err = 0;

    if (R == null || q == null || b == null || V == null) {
        pputs(prn, "Missing matrix in system Wald test!\n");
        return GretlError.data;
    }

    Rrows = gretl_matrix_rows(R);
    dfn = gretl_matrix_rows(R);

    Rbq = gretl_matrix_alloc(Rrows, 1);
    RvR = gretl_matrix_alloc(Rrows, Rrows);

    if (Rbq == null || RvR == null) {
        err = GretlError.alloc;
    } else {
        gretl_matrix_multiply(R, b, Rbq);
        gretl_matrix_subtract_from(Rbq, q);
        gretl_matrix_qform(R, GretlMatrixMod.none, V,
                           RvR, GretlMatrixMod.none);
        err = gretl_invert_symmetric_matrix(RvR);
    }

    if (!err) {
        test = gretl_scalar_qform(Rbq, RvR, &err);
    }

    if (!err) {
        double pval;

        if (dfu == 0) {
            pval = chisq_cdf_comp(dfn, test);
        } else {
            test /= dfn;
            pval = snedecor_cdf_comp(dfn, dfu, test);
        }

        record_test_result(test, pval);

        if (!(opt & Gretlopt.OPT_Q)) {
            if (dfu == 0) {
                pprintf(prn, "%s:\n".ptr, _("Wald test for the specified restrictions").ptr);
                pprintf(prn, "  %s(%d) = %g [%.4f]\n".ptr, _("Chi-square").ptr,
                        dfn, test, pval);
            } else {
                pprintf(prn, "%s:\n".ptr, _("F test for the specified restrictions").ptr);
                pprintf(prn, "  F(%d,%d) = %g [%.4f]\n".ptr, dfn, dfu, test, pval);
            }
            pputc(prn, '\n');
        }
    }

    gretl_matrix_free(Rbq);
    gretl_matrix_free(RvR);

    return err;
}

int get_wald_dfu(EquationSystem * sys, const(GretlMatrixRaw) * R)
{
    int dfu = 0;

    if (sys.method == GretlSystemMethods.ols || sys.method == GretlSystemMethods.tsls) {
        dfu = maybe_get_single_equation_dfu(sys, R);
    }
    if (dfu == 0) {
        dfu = system_get_dfu(sys);
    }

    return dfu;
}

int system_wald_test(EquationSystem * sys, const(GretlMatrixRaw) * R, const(GretlMatrixRaw) * q, Gretlopt opt, Prn * prn)
{
    const(GretlMatrixRaw) * b = sys.b;
    const(GretlMatrixRaw) * V = sys.vcv;
    int dfu = 0;

    if (b == null || V == null) {
        gretl_errmsg_set("restrict: no estimates are available");
        return GretlError.data;
    }

    if (sys.method == GretlSystemMethods.ols || sys.method == GretlSystemMethods.tsls) {
        /* use the F-form? */
        dfu = maybe_get_single_equation_dfu(sys, R);
    }

    return multi_eqn_wald_test(b, V, R, q, dfu, opt, prn);
}

int system_do_LR_test(EquationSystem * sys, double llu, Gretlopt opt, Prn * prn)
{
    double X2, pval, llr = sys.ll;
    int df = gretl_matrix_rows(sys.R);
    int err = 0;

    if (na(llr) || na(llu) || llr == 0.0 || llu == 0.0 || df <= 0) {
        fputs("bad or missing data in system LR test\n", stderr);
        return GretlError.data;
    }

    X2 = 2.0 * (llu - llr);
    pval = chisq_cdf_comp(df, X2);

    if (!(opt & Gretlopt.OPT_Q)) {
        pprintf(prn, "%s:\n".ptr, _("LR test for the specified restrictions").ptr);
        pprintf(prn, "  %s = %g\n".ptr, _("Restricted log-likelihood").ptr, llr);
        pprintf(prn, "  %s = %g\n".ptr, _("Unrestricted log-likelihood").ptr, llu);
        pprintf(prn, "  %s(%d) = %g [%.4f]\n".ptr, _("Chi-square").ptr, df, X2, pval);
        pputc(prn, '\n');
    }

    return err;
}

int sys_test_type(EquationSystem * sys)
{
    int ret = SysTestType.none;

    if (system_n_restrictions(sys) > 0) {
        if (sys.method == GretlSystemMethods.sur ||
            sys.method == GretlSystemMethods.wls) {
            if (sys.flags & SystemFlags.iterate) {
                ret = SysTestType.lr;
            } else {
                ret = SysTestType.f;
            }
        } else if (sys.method == GretlSystemMethods.ols ||
                   sys.method == GretlSystemMethods.tsls ||
                   sys.method == GretlSystemMethods.sls3) {
            ret = SysTestType.f;
        } else if (sys.method == GretlSystemMethods.liml) {
            /* experimental */
            ret = SysTestType.f;
        } else if (sys.method == GretlSystemMethods.fiml) {
            ret = SysTestType.lr;
        }
    }

    return ret;
}

/* Handle the case where sys->b and sys->vcv have been augmented in
   order to test some restrictions: we now have to trim these matrices
   down to their final_ size. The vector @b contains unrestricted
   coefficient estimates, so its length can be used to gauge the
   true number of coeffs.
*/

int shrink_b_and_vcv(const(GretlMatrixRaw) * b, EquationSystem * sys)
{
    int nc = gretl_vector_get_length(b);
    GretlMatrixRaw * V;
    double x;
    int i, j;

    if (sys.vcv == null) {
        return GretlError.data;
    }

    if (sys.vcv.rows == nc) {
        /* no-op (shouldn't happen) */
        return 0;
    }

    V = gretl_matrix_alloc(nc, nc);
    if (V == null) {
        return GretlError.alloc;
    }

    gretl_matrix_reuse(sys.b, nc, 1);

    for (i=0; i<nc; i++) {
        for (j=0; j<nc; j++) {
            x = gretl_matrix_get(sys.vcv, i, j);
            gretl_matrix_set(V, i, j, x);
        }
    }

    gretl_matrix_replace(&sys.vcv, V);

    return 0;
}

alias SystemEstFn = extern (C) int function(EquationSystem *, DatasetRaw *, Gretlopt, Prn *);

int estimate_with_test(EquationSystem * sys, DatasetRaw * dset, int stest, SystemEstFn system_est, Gretlopt opt, Prn * prn)
{
    GretlMatrixRaw * vcv = null;
    GretlMatrixRaw * b = null;
    double llu = 0.0;
    int err = 0;

    /* estimate the unrestricted system first */
    sys.flags &= ~SystemFlags.restrict;
    err = system_est(sys, dset, opt | Gretlopt.OPT_Q, prn);
    sys.flags ^= SystemFlags.restrict;
    if (err) {
        goto bailout;
    }

    /* grab the data from unrestricted estimation */
    b = sys.b;
    sys.b = null;
    if (stest == SysTestType.lr) {
        llu = sys.ll;
    } else if (stest == SysTestType.f) {
        vcv = sys.vcv;
        sys.vcv = null;
    }

    /* now estimate the restricted system */
    system_clear_results(sys);
    err = system_est(sys, dset, opt, prn);
    if (!err) {
        if (stest == SysTestType.lr) {
            err = system_do_LR_test(sys, llu, opt, prn);
        } else if (stest == SysTestType.f) {
            int dfu = get_wald_dfu(sys, sys.R);

            err = multi_eqn_wald_test(b, vcv, sys.R, sys.q,
                                      dfu, opt, prn);
        }
        shrink_b_and_vcv(b, sys);
    }

 bailout:

    gretl_matrix_free(b);
    gretl_matrix_free(vcv);

    return err;
}

void adjust_sys_flags_for_method(EquationSystem * sys, int method)
{
    char oldflags = cast(char) sys.flags;

    sys.flags = 0;

    /* the robust option is for now OLS-only */
    if (oldflags & SystemFlags.robust) {
        if (sys.method == GretlSystemMethods.ols) {
            sys.flags |= SystemFlags.robust;
        }
    }

    /* the iterate option is only available for WLS, SUR or 3SLS */
    if (oldflags & SystemFlags.iterate) {
        if (sys.method == GretlSystemMethods.wls ||
            sys.method == GretlSystemMethods.sur ||
            sys.method == GretlSystemMethods.sls3) {
            sys.flags |= SystemFlags.iterate;
        }
    }

    /* by default, apply df correction for single-equation methods */
    if (sys.method == GretlSystemMethods.ols ||
        sys.method == GretlSystemMethods.tsls ||
        sys.method == GretlSystemMethods.liml) {
        if (oldflags & SystemFlags.dfcorr) {
            sys.flags |= SystemFlags.dfcorr;
        }
    }

    /* carry forward the GEOMEAN flag */
    if (oldflags & SystemFlags.vcv_geomean) {
        sys.flags |= SystemFlags.vcv_geomean;
    }
}

int set_sys_flags_from_opt(EquationSystem * sys, Gretlopt opt)
{
    sys.flags = 0;

    /* the robust option is for now OLS-only */
    if (opt & Gretlopt.OPT_R) {
        if (sys.method == GretlSystemMethods.ols) {
            sys.flags |= SystemFlags.robust;
        } else {
            gretl_errmsg_sprintf(_("Invalid option '--%s'").ptr, "robust".ptr);
            return GretlError.badopt;
        }
    }

    /* the iterate option is available for WLS, SUR or 3SLS */
    if (opt & Gretlopt.OPT_I) {
        if (sys.method == GretlSystemMethods.wls ||
            sys.method == GretlSystemMethods.sur ||
            sys.method == GretlSystemMethods.sls3) {
            sys.flags |= SystemFlags.iterate;
        }
    }

    /* by default, apply a df correction for single-equation methods */
    if (sys.method == GretlSystemMethods.ols ||
        sys.method == GretlSystemMethods.tsls ||
        sys.method == GretlSystemMethods.liml) {
        if (!(opt & Gretlopt.OPT_N)) {
            sys.flags |= SystemFlags.dfcorr;
        }
    }

    if (opt & Gretlopt.OPT_M) {
        sys.flags |= SystemFlags.vcv_geomean;
    }
    if (opt & Gretlopt.OPT_Q) {
        sys.flags |= SystemFlags.quiet;
    }
    if (opt & Gretlopt.OPT_S) {
        /* estimating single equation via LIML */
        sys.flags |= SystemFlags.liml1;
    }

    return 0;
}

/**
 * equation_system_estimate:
 * @sys: pre-defined equation system.
 * @dset: dataset struct.
 * @opt: may include OPT_V for more verbose operation.
 * @prn: printing struct.
 *
 * Estimate a pre-defined equation system and print the results
 * to @prn.
 *
 * Returns: 0 on success, non-zero on error.
 */

int equation_system_estimate(EquationSystem * sys, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    SystemEstFn system_est;
    int err = 0;

    gretl_error_clear();

    if (sys.xlist == null || sys.biglist == null) {
        /* allow for the possibility that we're looking at a
           system restored from a session file
        */
        err = sys_check_lists(sys, dset);
        if (err) {
            return err;
        }
    }

    if (opt == Gretlopt.OPT_UNSET) {
        opt = Gretlopt.OPT_NONE;
        adjust_sys_flags_for_method(sys, sys.method);
    } else {
        err = set_sys_flags_from_opt(sys, opt);
        if (err) {
            return err;
        }
    }

    /* if restrictions are in place, reset the restricted flag */
    if (sys.R != null && sys.q != null) {
        sys.flags |= SystemFlags.restrict;
    }

    /* in case we're re-estimating */
    system_clear_results(sys);
    err = sys_rearrange_eqn_lists(sys, dset);

    if (!err) {
        system_est = cast(SystemEstFn) get_plugin_function("system_estimate");
        if (system_est == null) {
            err = 1;
        }
    }

    if (!err) {
        int stest = sys_test_type(sys);

        if (stest == SysTestType.notimp) {
            pputs(prn, _("Sorry, command not available for this estimator").ptr);
            pputc(prn, '\n');
            err = 1;
        } else if (stest != SysTestType.none) {
            err = estimate_with_test(sys, dset, stest,
                                     system_est, opt, prn);
        } else {
            err = system_est(sys, dset, opt, prn);
        }
    }

    if (sys.tslists != null) {
    lists_restore_insts(sys);
    }

    if (!err) {
        sys.smpl_t1 = dset.t1;
        sys.smpl_t2 = dset.t2;
        if (!(sys.flags & SystemFlags.liml1)) {
            set_as_last_model(sys, GretlObjType.sys);
        }
    }

    return err;
}

int system_adjust_t1t2(EquationSystem * sys, const(DatasetRaw) * dset)
{
    int err;

    if (sys.biglist == null) {
        fprintf(stderr, "system_adjust_t1t2: no 'biglist' present!\n".ptr);
        return GretlError.data;
    }

    sys.t1 = dset.t1;
    sys.t2 = dset.t2;

    err = list_adjust_sample(sys.biglist, &sys.t1, &sys.t2, dset, null);
    if (!err) {
        sys.T = sys.t2 - sys.t1 + 1;
    }

    return err;
}

/* construct a list containing all the variables referenced in the
   system */

int sys_make_biglist(EquationSystem * sys)
{
    int * biglist;
    const(int) * slist;
    const(Identity) * ident;
    int bign = 0;
    int i, j, vj, k;

    for (i=0; i<sys.neqns; i++) {
        bign += sys.lists[i][0] - gretl_list_has_separator(sys.lists[i]);
    }

    for (i=0; i<sys.nidents; i++) {
        bign += 1 + sys.idents[i].n_atoms;
    }

    if (sys.ilist != null) {
        bign += sys.ilist[0];
    }

    biglist = gretl_list_new(bign);
    if (biglist == null) {
        return GretlError.alloc;
    }

    for (i=1; i<=bign; i++) {
        /* invalidate all elements */
        biglist[i] = -1;
    }

    k = 0;

    /* system equations */
    for (i=0; i<sys.neqns; i++) {
        slist = sys.lists[i];
        for (j=1; j<=slist[0]; j++) {
            vj = slist[j];
            if (vj != LISTSEP && !in_gretl_list(biglist, vj)) {
                biglist[++k] = vj;
            }
        }
    }

    /* system identities */
    for (i=0; i<sys.nidents; i++) {
        ident = sys.idents[i];
        for (j=0; j<=ident.n_atoms; j++) {
            if (j == 0) {
                vj = ident.depvar;
            } else {
                vj = ident.atoms[j-1].varnum;
            }
            if (!in_gretl_list(biglist, vj)) {
                biglist[++k] = vj;
            }
        }
    }

    if (sys.ilist != null) {
        /* system instruments */
        for (j=1; j<=sys.ilist[0]; j++) {
            vj = sys.ilist[j];
            if (!in_gretl_list(biglist, vj)) {
                biglist[++k] = vj;
            }
        }
    }

    biglist[0] = k;

    free(sys.biglist);
    sys.biglist = biglist;

    return 0;
}

int sys_get_lag_src(const(char) * vname, const(DatasetRaw) * dset)
{
    int fd, src = series_index(dset, vname);

    if (src == dset.v && vname != null &&
        (fd = gretl_function_depth()) > 0) {
        int i;

        for (i=1; i<dset.v; i++) {
            if (fd == series_get_stack_level(dset, i) &&
                series_is_listarg(dset, i, null) &&
                !strcmp(dset.varname[i], vname)) {
                src = i;
                break;
            }
        }
    }

    return src;
}

int is_tsls_style_instrument(EquationSystem * sys, int v)
{
    int i, j, pos;

    for (i=0; i<sys.neqns; i++) {
        pos = gretl_list_separator_position(sys.lists[i]);
        if (pos > 0) {
            for (j=pos+1; j<=sys.lists[i][0]; j++) {
                if (sys.lists[i][j] == v) {
                    return 1;
                }
            }
        }
    }

    return 0;
}

/* Here we deal with the case where the user has specified the
   system equations "tsls-style": besides the left-hand side
   variables, we also need to identify as endogenous any
   regressors that do not appear as instruments.
*/

int tsls_style_augment_ylist(EquationSystem * sys, int ** pylist)
{
    const(int) * list;
    int pos, ny = (*pylist)[0];
    int i, j, vj;

    for (i=0; i<sys.neqns; i++) {
        list = sys.lists[i];
        pos = gretl_list_separator_position(list);
        if (pos == 0) {
            continue;
        }
        for (j=2; j<pos; j++) {
            vj = list[j];
            if (vj == 0 || in_gretl_list(*pylist, vj)) {
                continue;
            }
            if (!is_tsls_style_instrument(sys, vj)) {
                if (ny < sys.neqns) {
                    (*pylist)[++ny] = vj;
                    (*pylist)[0] = ny;
                } else {
                    gretl_list_append_term(pylist, vj);
                    if (*pylist == null) {
                        return GretlError.alloc;
                    }
                    ny++;
                }
            }
        }
    }

    return ny >= sys.neqns ? 0 : GretlError.data;
}

/* It's possible to specify instruments per equation, as in tsls,
   with a two-part regression list. But this cannot be mixed with
   any other means of specifying endogeneity/exogeneity of
   variables.
*/

int check_for_tsls_style_lists(EquationSystem * sys, int * err)
{
    int i, tsls_style = 0;

    for (i=0; i<sys.neqns; i++) {
        if (gretl_list_has_separator(sys.lists[i])) {
            tsls_style = 1;
            break;
        }
    }

    if (tsls_style) {
        if (sys.ylist != null || sys.ilist != null) {
            gretl_errmsg_set("You can't mix tsls-style lists with 'endog' or 'instr'");
            *err = GretlError.data;
        } else if (sys.nidents > 0) {
            gretl_errmsg_set("You can't mix identities with tsls-style lists");
            *err = GretlError.data;
        }
    }

    return tsls_style;
}

/* Given a "partial" system, with equations expressed in tsls
   style and with more endogenous variables than equations, shift
   the "excess" endogenous terms into @xplist. They will still
   be instrumented when the equations are estimated, but they
   will be treated "as if" exogenous when it comes to building
   the structural form of the system.
*/

void tsls_style_shift_vars(EquationSystem * sys, int * xplist)
{
    int i, vi, nxp = xplist[0];

    for (i=sys.ylist[0]; i>sys.neqns; i--) {
        vi = sys.ylist[i];
        gretl_list_delete_at_pos(sys.ylist, i);
        xplist[++nxp] = vi;
    }

    xplist[0] = nxp;
}

/* For consistency with the way in which per-equation results
   are organized, we should ensure that if a constant is present
   among the exogenous terms in a system it appears as the
   first element of the system's "xlist".
*/

void sys_xlist_reshuffle_const(int * list)
{
    int i, cpos = 0;

    for (i=1; i<=list[0]; i++) {
        if (list[i] == 0) {
            cpos = i;
            break;
        }
    }

    if (cpos > 1) {
        for (i=cpos; i>1; i--) {
            list[i] = list[i-1];
        }
        list[1] = 0;
    }
}

/* prior to system estimation, get all the required lists of variables
   in order
*/

int sys_check_lists(EquationSystem * sys, const(DatasetRaw) * dset)
{
    const(int) * slist;
    const(char) * vname;
    const(Identity) * ident;
    int user_ylist = (sys.ylist != null);
    int user_ilist = (sys.ilist != null);
    int * ylist = null;
    int * xplist = null;
    int src, lag, nlhs;
    int tsls_style;
    int i, j, k, vj;
    int err = 0;

    tsls_style = check_for_tsls_style_lists(sys, &err);
    if (err) {
        return err;
    }

    /* start an empty list for predetermined variables */

    sys.plist = gretl_null_list();
    if (sys.plist == null) {
        return GretlError.alloc;
    }

    /* Compose a minimal (and possibly incomplete) list of endogenous
       variables, including all vars that appear on the left-hand side
       of structural equations or identities.
    */

    nlhs = sys.neqns + sys.nidents;
    ylist = gretl_list_new(nlhs);
    if (ylist == null) {
        return GretlError.alloc;
    }

    k = 0;

    for (i=0; i<sys.neqns; i++) {
        slist = sys.lists[i];
        vj = slist[1];
        if (!in_gretl_list(ylist, vj)) {
            ylist[++k] = vj;
        }
    }

    if (tsls_style) {
        ylist[0] = k;
        err = tsls_style_augment_ylist(sys, &ylist);
        if (err) {
            goto bailout;
        }
    } else {
        for (i=0; i<sys.nidents; i++) {
            ident = sys.idents[i];
            vj = ident.depvar;
            if (!in_gretl_list(ylist, vj)) {
                ylist[++k] = vj;
            }
        }
        ylist[0] = k;
    }

    /* If the user gave a list of endogenous vars (in sys->ylist), check
       that it contains all the presumably endogenous variables we found
       above and recorded in ylist; otherwise use the ylist as computed
       above.
    */
    if (user_ylist) {
        for (j=1; j<=ylist[0]; j++) {
            vj = ylist[j];
            if (!in_gretl_list(sys.ylist, vj)) {
                gretl_errmsg_sprintf("%s appears on the left-hand side of an equation but is not marked as endogenous".ptr,
                                     dset.varname[vj]);
                err = GretlError.data;
                goto bailout;
            }
        }
    } else {
        sys.ylist = ylist;
        ylist = null;
    }

    /* Now compose a big list of all the variables in the system
       (referenced in structural equations or identities, or specified
       as instruments).
    */
    err = sys_make_biglist(sys);
    if (err) {
        goto bailout;
    }

    /* If the user gave an endogenous list, check that it doesn't
       contain any variables that are not actually present in the
       system.
    */
    if (user_ylist) {
        for (j=1; j<=sys.ylist[0]; j++) {
            vj = sys.ylist[j];
            if (!in_gretl_list(sys.biglist, vj)) {
                gretl_errmsg_sprintf("%s is marked as endogenous but is not present in the system".ptr,
                                     dset.varname[vj]);
                err = GretlError.data;
                goto bailout;
            }
        }
    }

    /* By elimination from biglist, using ylist, compose a maximal
       list of exogenous and predetermined variables, "xplist".
    */

    xplist = gretl_list_new(sys.biglist[0]);
    if (xplist == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    k = 0;
    for (j=1; j<=sys.biglist[0]; j++) {
        vj = sys.biglist[j];
        if (!in_gretl_list(sys.ylist, vj)) {
            xplist[++k] = vj;
        }
    }
    xplist[0] = k;

    /* If the user gave a list of instruments (sys->ilist), check
       that it does not contain anything that is in fact clearly
       endogenous; otherwise use the xplist we computed above.
    */
    if (user_ilist) {
        for (j=1; j<=sys.ilist[0]; j++) {
            vj = sys.ilist[j];
            if (in_gretl_list(sys.ylist, vj)) {
                gretl_errmsg_sprintf("%s is marked as an instrument but is endogenous".ptr,
                                     dset.varname[vj]);
                err = GretlError.data;
                goto bailout;
            }
        }
    } else {
        sys.ilist = gretl_list_copy(xplist);
        if (sys.ilist == null) {
            err = GretlError.alloc;
        }
    }

    if (!err && sys.ylist[0] != nlhs) {
        /* Note: check added 2009-08-10, modified 2012-04-05,
           moved in front of the following "xplist" block
           2015-10-29.
        */
        if (tsls_style && sys.ylist[0] > nlhs) {
            /* from the pov of the structural form, endogenous regressors
               without an equation should be treated "as if" exogenous?
            */
            tsls_style_shift_vars(sys, xplist);
        } else {
            gretl_errmsg_sprintf("Found %d endogenous variables but %d equations".ptr,
                                 sys.ylist[0], nlhs);
            err = GretlError.data;
        }
    }

    /* Now narrow down "xplist", removing variables that are in fact
       lags of endogenous variables (and recording them as such): this
       should leave a list of truly exogenous variables.
    */
    for (j=1; j<=xplist[0] && !err; j++) {
        vj = xplist[j];
        lag = series_get_lag(dset, vj);
        if (lag > 0) {
            vname = series_get_parent_name(dset, vj);
            if (vname != null) {
                src = sys_get_lag_src(vname, dset);
                if (in_gretl_list(sys.ylist, src)) {
                    err = add_predet_to_sys(sys, dset, vj, src, lag);
                    gretl_list_delete_at_pos(xplist, j--);
                }
            }
        }
    }

    if (!err) {
        sys_xlist_reshuffle_const(xplist);
    }

    if (!err) {
        free(sys.xlist);
        sys.xlist = xplist;
        xplist = null;
    }

 bailout:

    free(ylist);
    free(xplist);

    return err;
}

int sys_has_user_name(EquationSystem * sys)
{
    return (sys.name != null && *sys.name != '\0' &&
            !sys_anonymous(sys.name));
}

/**
 * equation_system_finalize:
 * @sys: pre-defined equation system.
 * @dset: dataset struct.
 * @opt: may include OPT_V for verbose operation, OPT_S
 * to permit estimation of a single equation, OPT_R for
 * a robust variance matrix (OLS only), OPT_N for no
 * df correction (OLS and TSLS only).
 * @prn: printing struct.
 *
 * Finalize an equation system, e.g. in response to "end system".
 * If the system has a specified name, we save it on a stack of
 * defined systems.  If it has a specified estimation method,
 * we go ahead and estimate it.  If it has both, we do both.
 *
 * Returns: 0 on success, non-zero on error.  If the system is
 * mal-formed, it is destroyed.
 */

int equation_system_finalize(EquationSystem * sys, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    int mineq = 1;
    int err = 0;

    gretl_error_clear();

    if (sys == null) {
        gretl_errmsg_set(_(nosystem));
        return 1;
    } else if (sys.neqns < mineq) {
        gretl_errmsg_set(_(toofew));
        equation_system_destroy(sys);
        return 1;
    } else if (sys.method >= GretlSystemMethods.max) {
        gretl_errmsg_set(_(badsystem));
        equation_system_destroy(sys);
        return 1;
    }

    err = sys_check_lists(sys, dset);

    if (!err && !(opt & Gretlopt.OPT_S) && sys_has_user_name(sys)) {
        /* Save the system for subsequent estimation: but note that we
           should not do this if given OPT_S, for single-equation
           LIML.
        */
        err = gretl_stack_object_as(sys, GretlObjType.sys, sys.name);
    }

    if (!err && sys.method >= 0) {
        if (sys.flags & SystemFlags.quiet) {
            opt |= Gretlopt.OPT_Q;
        }
        err = equation_system_estimate(sys, dset, opt, prn);
    }

    return err;
}

/* Implement the "estimate" command. The full form of the
   command line is

   estimate <sysname> method=<method>

   where <sysname> is the name of a previously defined
   equation system and <method> is a recognized estimation
   method.

   However, if the "last model" is an equation system, we
   can accept just

   estimate method=<method>

   where the system is not named but is implicitly the
   last model.

   The method=<method> field is optional, provided that the
   (named or last-model) system already has a specified
   estimation method; so the minimal form of the command
   line is simply "estimate".
*/

int estimate_named_system(const(char) * sysname, const(char) * param, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    EquationSystem * sys = null;
    int err = 0;

    if (sysname != null) {
        /* we got a name */
        if (sys_anonymous(sysname)) {
            sys = get_anonymous_equation_system();
        } else {
            sys = get_equation_system_by_name(sysname);
        }
        if (sys == null) {
            gretl_errmsg_sprintf(_("'%s': unrecognized name").ptr, sysname);
            err = GretlError.data;
        }
    } else {
        /* no name given: try "last model"? */
        GretlObjType type;
        void * ptr;

        ptr = get_last_model(&type);
        if (ptr == null || type != GretlObjType.sys) {
            gretl_errmsg_sprintf(_("%s: no system was specified").ptr, "estimate".ptr);
            err = GretlError.data;
        } else {
            sys = cast(EquationSystem *) ptr;
        }
    }

    if (err) {
        /* we haven't found a system to estimate yet */
        sys = get_anonymous_equation_system();
        if (sys != null) {
            gretl_error_clear();
            err = 0;
        }
    }

    if (!err) {
        int method;

        if (param != null) {
            method = sys_get_estimator(param);
        } else {
            method = sys.method;
        }

        if (method < 0 || method >= GretlSystemMethods.max) {
            gretl_errmsg_set("estimate: no valid method was specified");
            err = GretlError.data;
        } else {
            sys.method = method;
            err = equation_system_estimate(sys, dset, opt, prn);
        }
    }

    return err;
}

/* effective list length, allowance made for IVREG-style lists */

int get_real_list_length(const(int) * list)
{
    int i, len = list[0];

    for (i=1; i<=list[0]; i++) {
        if (list[i] == LISTSEP) {
            len = i - 1;
            break;
        }
    }

    return len;
}

/* below: used in formulating restriction matrices for
   an equation system */

int system_get_list_length(EquationSystem * sys, int i)
{
    if (i >= 0 && i < sys.neqns) {
        return get_real_list_length(sys.lists[i]);
    } else {
        return 0;
    }
}

int system_max_indep_vars(EquationSystem * sys)
{
    int i, nvi, nv = 0;

    for (i=0; i<sys.neqns; i++) {
        nvi = get_real_list_length(sys.lists[i]) - 1;
        if (nvi > nv) nv = nvi;
    }

    return nv;
}

int system_n_indep_vars(EquationSystem * sys)
{
    int i, nvi, nv = 0;

    for (i=0; i<sys.neqns; i++) {
        nvi = get_real_list_length(sys.lists[i]) - 1;
        nv += nvi;
    }

    return nv;
}

const(char) * system_short_string(const(Model) * pmod)
{
    int i = gretl_model_get_int(pmod, "method");

    return gretl_system_short_strings[i];
}

void equation_system_set_name(EquationSystem * sys, const(char) * name)
{
    if (name == sys.name) {
        return;
    }

    if (sys.name == null) {
        sys.name = cast(typeof(sys.name)) malloc(MAXSAVENAME);
    }

    if (sys.name != null) {
        *sys.name = '\0';
        strncat(sys.name, name, MAXSAVENAME - 1);
    }
}

int * compose_ivreg_list(EquationSystem * sys, int i)
{
    int * list;
    int j, k1, k2;

    if (i >= sys.neqns) {
        return null;
    }

    k1 = sys.lists[i][0];
    k2 = sys.ilist[0];

    list = gretl_list_new(k1 + k2 + 1);
    if (list == null) {
        return null;
    }

    for (j=1; j<=list[0]; j++) {
        if (j <= k1) {
            list[j] = sys.lists[i][j];
        } else if (j == k1 + 1) {
            list[j] = LISTSEP;
        } else {
            list[j] = sys.ilist[j - (k1 + 1)];
        }
    }

    return list;
}

int system_normality_test(EquationSystem * sys, Gretlopt opt, Prn * prn)
{
    int err = 0;

    if (sys.E == null || sys.S == null) {
        err = 1;
    } else {
        err = multivariate_normality_test(sys.E,
                                          sys.S,
                                          opt,
                                          prn);
    }

    return err;
}

double * system_get_resid_series(EquationSystem * sys, int eqnum, DatasetRaw * dset, int * err)
{
    double * u = null;
    int t;

    if (sys.E == null) {
        *err = GretlError.data;
        return null;
    }

    u = cast(typeof(u)) malloc(dset.n * (*u).sizeof);
    if (u == null) {
        *err = GretlError.alloc;
        return null;
    }

    for (t=0; t<dset.n; t++) {
        if (t < sys.t1 || t > sys.t2) {
            u[t] = NADBL;
        } else {
            u[t] = gretl_matrix_get(sys.E, t - sys.t1, eqnum);
        }
    }

    return u;
}

char * system_get_full_string(EquationSystem * sys)
{
    const(char) * lstr = gretl_system_long_strings[sys.method];
    char * ret;

    if (sys.flags & SystemFlags.iterate) {
    ret = g_strdup_printf(_("iterated %s").ptr, _(lstr));
    } else {
    ret = g_strdup(_(lstr));
    }

    return ret;
}

/* simple accessor functions */

int system_want_df_corr(EquationSystem * sys)
{
    return sys.flags & SystemFlags.dfcorr;
}

int system_n_restrictions(EquationSystem * sys)
{
    int nr = 0;

    if (sys.R != null && (sys.flags & SystemFlags.restrict)) {
        nr = gretl_matrix_rows(sys.R);
    }

    return nr;
}

int * system_get_list(EquationSystem * sys, int i)
{
    if (i >= sys.neqns) return null;

    return sys.lists[i];
}

int system_get_depvar(EquationSystem * sys, int i)
{
    if (i >= sys.neqns) return 0;

    return sys.lists[i][1];
}

int system_get_method(EquationSystem * sys)
{
    return sys.method;
}

int * system_get_endog_vars(EquationSystem * sys)
{
    return sys.ylist;
}

int * system_get_instr_vars(EquationSystem * sys)
{
    return sys.ilist;
}

void system_attach_coeffs(EquationSystem * sys, GretlMatrixRaw * b)
{
    gretl_matrix_replace(&sys.b, b);
}

void system_attach_vcv(EquationSystem * sys, GretlMatrixRaw * vcv)
{
    gretl_matrix_replace(&sys.vcv, vcv);
}

void system_attach_sigma(EquationSystem * sys, GretlMatrixRaw * S)
{
    gretl_matrix_replace(&sys.S, S);
}

void system_attach_uhat(EquationSystem * sys, GretlMatrixRaw * E)
{
    gretl_matrix_replace(&sys.E, E);
}

Model * system_get_model(EquationSystem * sys, int i)
{
    if (sys == null || sys.models == null || i >= sys.neqns) {
        return null;
    }

    return sys.models[i];
}

/* for case of applying df correction to cross-equation
   residual covariance matrix */

int system_vcv_geomean(EquationSystem * sys)
{
    return sys.flags & SystemFlags.vcv_geomean;
}

int sys_eqn_indep_coeffs(EquationSystem * sys, int eq)
{
    int nc;

    if (eq >= sys.neqns) {
        nc = 0;
    } else if (sys.R == null) {
        nc = sys.lists[eq][0] - 1;
    } else {
        int nr = gretl_matrix_rows(sys.R);
        int cols = gretl_matrix_cols(sys.R);
        int cmax, cmin = 0;
        int single, i, j;

        nc = sys.lists[eq][0] - 1;

        for (i=0; i<eq; i++) {
            cmin += sys.lists[i][0] - 1;
        }

        cmax = cmin + nc;

        for (i=0; i<nr; i++) {
            single = 1;
            for (j=0; j<cols; j++) {
                if (j >= cmin && j < cmax) {
                    if (gretl_matrix_get(sys.R, i, j) != 0.0) {
                        single = 0;
                        break;
                    }
                }
            }
            if (single) {
                nc--;
            }
        }
    }

    return nc;
}

double system_vcv_denom(EquationSystem * sys, int i, int j)
{
    double den = sys.T;

    if ((sys.flags & SystemFlags.vcv_geomean) &&
        i < sys.neqns && j < sys.neqns) {
        int ki = sys_eqn_indep_coeffs(sys, i);

        if (j == i) {
            den = sys.T - ki;
        } else {
            int kj = sys_eqn_indep_coeffs(sys, j);

            den = (sys.T - ki) * (sys.T - kj);
            den = sqrt(den);
        }
    }

    return den;
}

/* for system over-identification test */

int system_get_overid_df(EquationSystem * sys)
{
    int gl, i, k = 0;

    gl = sys.neqns * sys.ilist[0];

    for (i=0; i<sys.neqns; i++) {
        k += sys.lists[i][0] - 1;
    }

    return gl - k;
}

/* dealing with identities (FIML, LIML)
   (really also needed for LIML? not clear why...) */

int rhs_var_in_identity(EquationSystem * sys, int lhsvar, int rhsvar)
{
    const(Identity) * ident;
    int i, j;

    for (i=0; i<sys.nidents; i++) {
        ident = sys.idents[i];
        if (ident.depvar == lhsvar) {
            for (j=0; j<ident.n_atoms; j++) {
                if (ident.atoms[j].varnum == rhsvar) {
                    return (ident.atoms[j].op == IdentityOp.plus)? 1 : -1;
                }
            }
        }
    }

    return 0;
}

void destroy_ident(Identity * ident)
{
    free(ident.atoms);
    free(ident);
}

Identity * ident_new(int nv)
{
    Identity * ident;

    ident = cast(typeof(ident)) malloc((*ident).sizeof);
    if (ident == null) return null;

    ident.n_atoms = nv;
    ident.atoms = cast(typeof(ident.atoms)) malloc(nv * (*ident.atoms).sizeof);
    if (ident.atoms == null) {
        free(ident);
        ident = null;
    }

    return ident;
}

Identity * parse_identity(const(char) * str, DatasetRaw * dset, int * err)
{
    Identity * ident;
    char * test;
    const(char) * p = str;
    char[VNAMELEN] vname;
    int lag, gotop = 0;
    int i, v, nv, len;

    /* First pass: count the elements and check everything for validity.
       nv is the number of variables appearing on the right-hand side
       of the identity (nv >= 1).
    */

    i = 0;
    while (*p && !*err) {
        p += cast(int) strspn(p, " ");
        if (i == 0) {
            /* left-hand side variable */
            *err = extract_varname(vname.ptr, p, &len);
            if (!*err) {
                v = current_series_index(dset, vname.ptr);
                if (v < 0) {
                    *err = GretlError.unkvar;
                } else {
                    p += len;
                    p += cast(int) strspn(p, " ");
                    if (*p != '=') {
                        *err = GretlError.parse;
                    } else {
                        p++;
                        i++;
                        gotop = 1;
                    }
                }
            }
        } else if (*p == '+' || *p == '-') {
            gotop = 1;
            p++;
        } else if (i > 1 && !gotop) {
            /* rhs: no + or - given */
            *err = GretlError.parse;
        } else {
            /* right-hand size variable (may be lag) */
            *err = extract_varname(vname.ptr, p, &len);
            if (!*err && gotop && len == 0) {
                /* dangling operator */
                *err = GretlError.parse;
            }
            if (!*err) {
                v = current_series_index(dset, vname.ptr);
                if (v < 0) {
                    *err = GretlError.unkvar;
                } else {
                    p += len;
                    if (*p == '(') {
                        lag = cast(int) strtol(p + 1, &test, 10);
                        if (*test != ')') {
                            *err = GretlError.parse;
                        } else if (lag >= 0) {
                            *err = GretlError.data;
                        } else {
                            v = laggenr(v, -lag, dset);
                            if (v < 0) {
                                *err = GretlError.data;
                            } else {
                                p = test + 1;
                            }
                        }
                    }
                    i++;
                    gotop = 0;
                }
            }
        }
    }

    if (gotop) {
        /* trailing operator */
        *err = GretlError.parse;
    }

    if (!*err) {
        nv = i - 1;
        if (nv < 1) {
            *err = GretlError.parse;
        } else {
            ident = ident_new(nv);
            if (ident == null) {
                *err = GretlError.alloc;
                return null;
            }
        }
    }

    if (*err) {
        return null;
    }

    ident.atoms[0].op = IdentityOp.plus;

    /* Second pass: fill out the identity */

    p = str;
    i = 0;
    while (*p) {
        p += cast(int) strspn(p, " ");
        if (i == 0) {
            extract_varname(vname.ptr, p, &len);
            ident.depvar = series_index(dset, vname.ptr);
            p = strchr(p, '=') + 1;
            i++;
        } else if (*p == '+' || *p == '-') {
            ident.atoms[i-1].op = (*p == '+')? IdentityOp.plus : IdentityOp.minus;
            p++;
        } else {
            extract_varname(vname.ptr, p, &len);
            v = series_index(dset, vname.ptr);
            p += len;
            if (*p == '(') {
                lag = cast(int) strtol(p + 1, &test, 10);
                v = laggenr(v, -lag, dset);
                p = test + 1;
            }
            ident.atoms[i-1].varnum = v;
            i++;
        }
    }

    return ident;
}

/* add information regarding a predetermined regressor */

int add_predet_to_sys(EquationSystem * sys, const(DatasetRaw) * dset, int id, int src, int lag)
{
    int n = sys.plist[0];
    int * test;
    Predet * pre;

    if (id < 0 || src < 0 || id >= dset.v || src >= dset.v) {
        /* something screwy */
        return GretlError.data;
    }

    if (in_gretl_list(sys.plist, id)) {
        /* already present */
        return 0;
    }

    pre = cast(typeof(pre)) realloc(sys.pre_vars, (n + 1) * (*pre).sizeof);
    if (pre == null) {
        return GretlError.alloc;
    }

    sys.pre_vars = pre;
    sys.pre_vars[n].id = id;
    sys.pre_vars[n].src = src;
    sys.pre_vars[n].lag = lag;

    /* add to list of predetermined variable IDs */
    test = gretl_list_append_term(&sys.plist, id);
    if (test == null) {
        return GretlError.alloc;
    }

    return 0;
}

int add_identity_to_sys(EquationSystem * sys, const(char) * line, DatasetRaw * dset)
{
    Identity ** pident;
    Identity * ident;
    int ni = sys.nidents;
    int err = 0;

    ident = parse_identity(line, dset, &err);
    if (ident == null) return err;

    /* connect the identity to the equation system */
    pident = cast(typeof(pident)) realloc(sys.idents, (ni + 1) * (*sys.idents).sizeof);
    if (pident == null) {
        destroy_ident(ident);
        return GretlError.alloc;
    }

    sys.idents = pident;
    sys.idents[ni] = ident;
    sys.nidents += 1;

    return 0;
}

int add_aux_list_to_sys(EquationSystem * sys, const(char) * line, DatasetRaw * dset, int which)
{
    int * list;
    int err = 0;

    if (which == AuxListType.endog) {
        if (sys.ylist != null) {
            gretl_errmsg_set("Only one list of endogenous variables may be given");
            return 1;
        }
    } else if (which == AuxListType.instr) {
        if (sys.ilist != null) {
            gretl_errmsg_set("Only one list of instruments may be given");
            return 1;
        }
    } else {
        return GretlError.data;
    }

    line += cast(int) strspn(line, " ");

    list = generate_list(line, dset, GretlCmdIndex.system, &err);

    if (!err) {
        if (which == AuxListType.endog) {
            sys.ylist = list;
        } else {
            sys.ilist = list;
        }
    }

    return err;
}

/**
 * system_parse_line:
 * @sys: initialized equation system.
 * @line: command line.
 * @dset: dataset struct.
 *
 * Modifies @sys according to the command supplied in @line,
 * which must start with "identity" (and supply an identity
 * to be added to the system, or "endog" (and supply a list
 * of endogenous variables), or "instr" (and supply a list of
 * instrumental variables).
 *
 * Returns: 0 on success, non-zero on failure, in which case
 * @sys is destroyed.
 */

int system_parse_line(EquationSystem * sys, const(char) * line, DatasetRaw * dset)
{
    int err = 0;

    gretl_error_clear();

    line += cast(int) strspn(line, " ");

    if (!strncmp(line, "identity ", 9)) {
        err = add_identity_to_sys(sys, line + 8, dset);
    } else if (!strncmp(line, "endog ", 6)) {
        err = add_aux_list_to_sys(sys, line + 5, dset, AuxListType.endog);
    } else if (!strncmp(line, "instr ", 6)) {
        err = add_aux_list_to_sys(sys, line + 5, dset, AuxListType.instr);
    } else {
        err = GretlError.parse;
    }

    if (err) {
        equation_system_destroy(sys);
    }

    return err;
}

void system_set_restriction_matrices(EquationSystem * sys, GretlMatrixRaw * R, GretlMatrixRaw * q)
{
    system_clear_restrictions(sys);

    sys.R = R;
    sys.q = q;

    sys.flags |= SystemFlags.restrict;
}

double * equation_system_get_series(EquationSystem * sys, const(DatasetRaw) * dset, int idx, const(char) * key, int * err)
{
    const(GretlMatrixRaw) * M = null;
    double * x = null;
    const(char) * msel;
    int t, col = 0;

    if (sys == null || (idx != ModelDataIndex.uhat && idx != ModelDataIndex.yhat)) {
        *err = GretlError.badstat;
        return null;
    }

    msel = strchr(key, '[');
    if (msel == null || sscanf(msel, "[,%d]".ptr, &col) != 1) {
        *err = GretlError.parse;
    } else if (col <= 0 || col > sys.neqns) {
        *err = GretlError.data;
    }

    if (!*err) {
        M = (idx == ModelDataIndex.uhat)? sys.E : sys.yhat;
        if (M == null) {
            *err = GretlError.data;
        }
    }

    if (!*err) {
        x = cast(typeof(x)) malloc(dset.n * (*x).sizeof);
        if (x == null) {
            *err = GretlError.alloc;
        }
    }

    if (!*err) {
        int s = 0;

        col--; /* switch to 0-based */
        for (t=0; t<dset.n; t++) {
            if (t < sys.t1 || t > sys.t2) {
                x[t] = NADBL;
            } else {
                x[t] = gretl_matrix_get(M, s++, col);
            }
        }
    }

    return x;
}

int system_add_yhat_matrix(EquationSystem * sys)
{
    double x, avc;
    int k = 0;
    int i, s, t;

    sys.yhat = gretl_matrix_alloc(sys.T, sys.neqns);
    if (sys.yhat == null) {
        return GretlError.alloc;
    }

    gretl_matrix_set_t1(sys.yhat, sys.t1);
    gretl_matrix_set_t2(sys.yhat, sys.t2);

    for (i=0; i<sys.neqns; i++) {
        s = 0;
        for (t=sys.t1; t<=sys.t2; t++) {
            x = sys.models[i].yhat[t];
            gretl_matrix_set(sys.yhat, s++, i, x);
        }
        k += sys.models[i].ncoeff;
    }

    k -= system_n_restrictions(sys);

    avc = k / cast(double) sys.neqns;
    sys.df = cast(int) (sys.T - floor(avc));

    return 0;
}

GretlMatrixRaw * get_stderr_vec(const(GretlMatrixRaw) * V)
{
    GretlMatrixRaw * se;
    int k = V.rows;

    se = gretl_column_vector_alloc(k);

    if (se != null) {
        double x;
        int i;

        for (i=0; i<k; i++) {
            x = gretl_matrix_get(V, i, i);
            se.val[i] = sqrt(x);
        }
    }

    return se;
}

/* retrieve a copy of a specified matrix from an equation
   system */

GretlMatrixRaw * equation_system_get_matrix(EquationSystem * sys, int idx, int * err)
{
    GretlMatrixRaw * M = null;

    if (sys == null) {
        *err = GretlError.badstat;
        return null;
    }

    switch (idx) {
    case ModelDataIndex.coeff:
        if (sys.b == null) {
            *err = GretlError.badstat;
        } else {
            M = gretl_matrix_copy(sys.b);
        }
        break;
    case ModelDataIndex.uhat:
        M = gretl_matrix_copy(sys.E);
        break;
    case ModelDataIndex.yhat:
        M = gretl_matrix_copy(sys.yhat);
        break;
    case ModelDataIndex.vcv:
    case ModelDataIndex.se:
        if (sys.vcv == null) {
            *err = GretlError.badstat;
        } else if (idx == ModelDataIndex.se) {
            M = get_stderr_vec(sys.vcv);
        } else {
            M = gretl_matrix_copy(sys.vcv);
        }
        break;
    case ModelDataIndex.sigma:
        M = gretl_matrix_copy(sys.S);
        break;
    case ModelDataIndex.sysgam:
        if (sys.Gamma == null) {
            *err = GretlError.badstat;
        } else {
            M = gretl_matrix_copy(sys.Gamma);
        }
        break;
    case ModelDataIndex.sysa:
        if (sys.A == null) {
            *err = GretlError.badstat;
        } else {
            M = gretl_matrix_copy(sys.A);
        }
        break;
    case ModelDataIndex.sysb:
        if (sys.B == null) {
            *err = GretlError.badstat;
        } else {
            M = gretl_matrix_copy(sys.B);
        }
        break;
    default:
        *err = GretlError.badstat;
        break;
    }

    if (M == null && !*err) {
        *err = GretlError.alloc;
    }

    return M;
}

int highest_numbered_var_in_system(EquationSystem * sys, const(DatasetRaw) * dset)
{
    int i, j, v, vmax = 0;

    if (sys.biglist != null) {
        for (j=1; j<=sys.biglist[0]; j++) {
            v = sys.biglist[j];
            if (v > vmax) {
                vmax = v;
            }
        }
    } else {
        /* should not happen */
        for (i=0; i<sys.neqns; i++) {
            for (j=1; j<=sys.lists[i][0]; j++) {
                v = sys.lists[i][j];
                if (v == LISTSEP || v >= dset.v) {
                    /* temporary variables, already gone? */
                    continue;
                }
                if (v > vmax) {
                    vmax = v;
                }
            }
        }
    }

    return vmax;
}

Identity * sys_retrieve_identity(xmlNodePtr node, int * err)
{
    Identity * ident;
    xmlNodePtr cur;
    int n_atoms, depvar;
    int i, got = 0;

    got += gretl_xml_get_prop_as_int(node, "n_atoms", &n_atoms);
    got += gretl_xml_get_prop_as_int(node, "depvar", &depvar);
    if (got < 2) {
        *err = GretlError.data;
        return null;
    }

    ident = ident_new(n_atoms);
    if (ident == null) {
        *err = GretlError.alloc;
        return null;
    }

    ident.depvar = depvar;

    cur = node.children;

    i = 0;
    while (cur != null && !*err) {
        if (!xmlStrcmp(cur.name, cast(XUC) "id_atom")) {
            got = gretl_xml_get_prop_as_int(cur, "op", &ident.atoms[i].op);
            got += gretl_xml_get_prop_as_int(cur, "varnum", &ident.atoms[i].varnum);
            if (got < 2) {
                *err = GretlError.data;
            } else {
                i++;
            }
        }
        cur = cur.next;
    }

    if (!*err && i != n_atoms) {
        *err = GretlError.data;
    }

    if (*err) {
        destroy_ident(ident);
        ident = null;
    }

    return ident;
}

EquationSystem * equation_system_from_XML(xmlNodePtr node, xmlDocPtr doc, int * err)
{
    EquationSystem * sys;
    xmlNodePtr cur;
    char * sname;
    int method = 0;
    int i, j, got = 0;

    got += gretl_xml_get_prop_as_string(node, "name", &sname);
    got += gretl_xml_get_prop_as_int(node, "method", &method);

    if (got < 2) {
        *err = GretlError.data;
        return null;
    }

    sys = equation_system_new(method, sname, err);
    if (*err) {
        return null;
    }

    got = 0;
    got += gretl_xml_get_prop_as_int(node, "n_equations", &sys.neqns);
    got += gretl_xml_get_prop_as_int(node, "nidents", &sys.nidents);
    got += gretl_xml_get_prop_as_int(node, "flags", &sys.flags);

    if (got < 3) {
        *err = GretlError.data;
        goto bailout;
    }

    gretl_xml_get_prop_as_int(node, "order", &sys.order);

    sys.lists = cast(typeof(sys.lists)) malloc(sys.neqns * sys.lists.sizeof);
    if (sys.lists == null) {
        *err = GretlError.alloc;
        goto bailout;
    }

    if (sys.nidents > 0) {
        sys.idents = cast(typeof(sys.idents)) malloc(sys.nidents * sys.idents.sizeof);
        if (sys.idents == null) {
            *err = GretlError.alloc;
            goto bailout;
        }
    }

    cur = node.children;

    i = j = 0;
    while (cur != null && !*err) {
        if (!xmlStrcmp(cur.name, cast(XUC) "eqnlist")) {
            sys.lists[i++] = gretl_xml_get_list(cur, doc, err);
        } else if (!xmlStrcmp(cur.name, cast(XUC) "endog_vars")) {
            sys.ylist = gretl_xml_get_list(cur, doc, err);
        } else if (!xmlStrcmp(cur.name, cast(XUC) "instr_vars")) {
            sys.ilist = gretl_xml_get_list(cur, doc, err);
        } else if (!xmlStrcmp(cur.name, cast(XUC) "identity")) {
            sys.idents[j++] = sys_retrieve_identity(cur, err);
        } else if (!xmlStrcmp(cur.name, cast(XUC) "R")) {
            sys.R = gretl_xml_get_matrix(cur, doc, err);
        } else if (!xmlStrcmp(cur.name, cast(XUC) "q")) {
            sys.q = gretl_xml_get_matrix(cur, doc, err);
        }
        cur = cur.next;
    }

    if (!*err && (i != sys.neqns || j != sys.nidents)) {
        *err = GretlError.data;
    }

    if (*err) {
        equation_system_destroy(sys);
        sys = null;
    }

 bailout:

    return sys;
}

void xml_print_identity(Identity * ident, Prn * prn)
{
    int i;

    pprintf(prn, "<identity n_atoms=\"%d\" depvar=\"%d\">\n".ptr,
            ident.n_atoms, ident.depvar);
    for (i=0; i<ident.n_atoms; i++) {
        pprintf(prn, " <id_atom op=\"%d\" varnum=\"%d\"/>\n".ptr,
                ident.atoms[i].op, ident.atoms[i].varnum);
    }
    pputs(prn, "</identity>\n");
}

int equation_system_serialize(EquationSystem * sys, SavedObjectFlags flags, Prn * prn)
{
    char * xmlname = null;
    int tsls_style = 0;
    int i, err = 0;

    if (sys.name == null || *sys.name == '\0') {
        xmlname = gretl_strdup("none");
    } else {
        xmlname = gretl_xml_encode(sys.name);
    }

    pprintf(prn, "<gretl-equation-system name=\"%s\" saveflags=\"%d\" method=\"%d\" ".ptr,
            xmlname, flags, sys.method);
    free(xmlname);

    pprintf(prn, "n_equations=\"%d\" nidents=\"%d\" flags=\"%d\" order=\"%d\">\n".ptr,
            sys.neqns, sys.nidents, sys.flags, sys.order);

    for (i=0; i<sys.neqns; i++) {
        gretl_xml_put_tagged_list("eqnlist", sys.lists[i], prn);
    }

    for (i=0; i<sys.neqns; i++) {
        if (gretl_list_has_separator(sys.lists[i])) {
            tsls_style = 1;
            break;
        }
    }

    if (!tsls_style) {
        gretl_xml_put_tagged_list("endog_vars", sys.ylist, prn);
        gretl_xml_put_tagged_list("instr_vars", sys.ilist, prn);
    }

    for (i=0; i<sys.nidents; i++) {
        xml_print_identity(sys.idents[i], prn);
    }

    gretl_matrix_serialize(sys.R, "R", prn);
    gretl_matrix_serialize(sys.q, "q", prn);

    pputs(prn, "</gretl-equation-system>\n");

    return err;
}

GretlBundle * bundlize_diag_test(EquationSystem * sys)
{
    GretlBundle * b = null;
    double test = NADBL;
    double pval = NADBL;
    int k = sys.S.rows;
    int df = k * (k - 1) / 2;

    if (sys.method == GretlSystemMethods.sur && sys.iters > 0) {
        if (!na(sys.ldet) && sys.diag_test > 0.0) {
            test = sys.T * (sys.diag_test - sys.ldet);
            pval = chisq_cdf_comp(df, test);

        }
    } else {
        test = sys.diag_test;
        pval = chisq_cdf_comp(df, test);
    }

    if (!na(test) && !na(pval)) {
        b = gretl_bundle_new();
        if (b != null) {
            gretl_bundle_set_scalar(b, "test", test);
            gretl_bundle_set_scalar(b, "pvalue", pval);
            gretl_bundle_set_int(b, "df", df);
        }
    }

    return b;
}

int equation_system_bundlize(EquationSystem * sys, GretlBundle * b)
{
    const(char) * s;
    char[16] lname;
    int tsls_style = 0;
    int i, err = 0;

    gretl_bundle_set_string(b, "command", "system");

    gretl_bundle_set_int(b, "neqns", sys.neqns);
    s = system_method_short_string(sys.method);
    if (s != null) {
        gretl_bundle_set_string(b, "method", s);
    }
    if (sys.flags & SystemFlags.robust) {
        gretl_bundle_set_int(b, "robust", 1);
    }

    /* convert to 1-based for user-space */
    gretl_bundle_set_int(b, "t1", sys.t1 + 1);
    gretl_bundle_set_int(b, "t2", sys.t2 + 1);
    gretl_bundle_set_int(b, "sample_t1", sys.smpl_t1 + 1);
    gretl_bundle_set_int(b, "sample_t2", sys.smpl_t2 + 1);

    gretl_bundle_set_matrix(b, "coeff", sys.b);
    gretl_bundle_set_matrix(b, "vcv",   sys.vcv);
    gretl_bundle_set_matrix(b, "sigma", sys.S);
    gretl_bundle_set_matrix(b, "uhat",  sys.E);
    gretl_bundle_set_matrix(b, "yhat",  sys.yhat);

    if (sys.Gamma != null) {
        gretl_bundle_set_matrix(b, "Gamma", sys.Gamma);
    }
    if (sys.A != null) {
        gretl_bundle_set_matrix(b, "A", sys.A);
    }
    if (sys.B != null) {
        gretl_bundle_set_matrix(b, "B", sys.B);
    }

    /* per-equation model structs? */

    for (i=0; i<sys.neqns; i++) {
        sprintf(lname.ptr, "list%d".ptr, i+1);
        gretl_bundle_set_list(b, lname.ptr, sys.lists[i]);
    }

    for (i=0; i<sys.neqns; i++) {
        if (gretl_list_has_separator(sys.lists[i])) {
            tsls_style = 1;
            break;
        }
    }

    if (!tsls_style) {
        gretl_bundle_set_list(b, "endog_vars", sys.ylist);
        gretl_bundle_set_list(b, "instr_vars", sys.ilist);
    }

    if (sys.xlist != null) {
    gretl_bundle_set_list(b, "xlist", sys.xlist);
    }

    if (sys.plist != null) {
        gretl_bundle_set_list(b, "predet_vars", sys.plist);
    }

    /* restrictions? */
    if (sys.R != null) {
        gretl_bundle_set_matrix(b, "R", sys.R);
    }
    if (sys.q != null) {
        gretl_bundle_set_matrix(b, "q", sys.q);
    }

    if (sys.diag_test > 0) {
        GretlBundle * dt = bundlize_diag_test(sys);

        if (dt != null) {
            gretl_bundle_donate_data(b, "diag_test", dt,
                                     GretlType.bundle);
        }
    }

    return err;
}

int sur_ols_diag(EquationSystem * sys)
{
    double s2, ls2sum = 0.0;
    int i, err = 0;

    for (i=0; i<sys.neqns; i++) {
        s2 = gretl_model_get_double(sys.models[i], "ols_sigma_squared");
        if (na(s2)) {
            err = 1;
            break;
        }
        ls2sum += log(s2);
    }

    if (!err) {
        sys.diag_test = ls2sum;
    }

    return err;
}

void print_system_overid_test(EquationSystem * sys, Prn * prn)
{
    int tex = tex_format(prn);
    int df = system_get_overid_df(sys);
    double pv;

    if (sys.method == GretlSystemMethods.fiml && df > 0) {
        double X2;

        if (na(sys.ll) || na(sys.llu) ||
            sys.ll == 0.0 || sys.llu == 0.0) {
            return;
        }

        X2 = 2.0 * (sys.llu - sys.ll);
        pv = chisq_cdf_comp(df, X2);

        if (tex) {
            pprintf(prn, "%s:\\\\\n".ptr, _("LR over-identification test").ptr);
            if (sys.ll < 0) {
                pprintf(prn, "  %s = $-$%g".ptr, _("Restricted log-likelihood").ptr, -sys.ll);
            } else {
                pprintf(prn, "  %s = %g".ptr, _("Restricted log-likelihood").ptr, sys.ll);
            }
            gretl_prn_newline(prn);
            if (sys.llu < 0) {
                pprintf(prn, "  %s = $-$%g".ptr, _("Unrestricted log-likelihood").ptr, -sys.llu);
            } else {
                pprintf(prn, "  %s = %g".ptr, _("Unrestricted log-likelihood").ptr, sys.llu);
            }
            gretl_prn_newline(prn);
            pprintf(prn, "  $\\chi^2(%d)$ = %g [%.4f]\n".ptr, df, X2, pv);
        } else {
            pprintf(prn, "%s:\n".ptr, _("LR over-identification test").ptr);
            pprintf(prn, "  %s = %g\n".ptr, _("Restricted log-likelihood").ptr, sys.ll);
            pprintf(prn, "  %s = %g\n".ptr, _("Unrestricted log-likelihood").ptr, sys.llu);
            pprintf(prn, "  %s(%d) = %g [%.4f]\n\n".ptr, _("Chi-square").ptr, df, X2, pv);
        }
    } else if ((sys.method == GretlSystemMethods.sls3 ||
                sys.method == GretlSystemMethods.sur) && df > 0) {
        if (na(sys.X2) || sys.X2 <= 0.0) {
            if (!tex) {
                pputs(prn, _("Warning: the Hansen-Sargan over-identification test failed.\nThis probably indicates that the estimation problem is ill-conditioned.\n").ptr);
                pputc(prn, '\n');
            }
            return;
        }

        pv = chisq_cdf_comp(df, sys.X2);

        if (tex) {
            pprintf(prn, "\\noindent %s:\\\\\n".ptr,
                    _("Hansen--Sargan over-identification test").ptr);
            pprintf(prn, "  $\\chi^2(%d)$ = %g [%.4f]\\\\\n".ptr, df, sys.X2, pv);
        } else {
            pprintf(prn, "%s:\n".ptr, _("Hansen-Sargan over-identification test").ptr);
            pprintf(prn, "  %s(%d) = %g [%.4f]\n\n".ptr, _("Chi-square").ptr,
                    df, sys.X2, pv);
        }
    }
}

int system_diag_test(EquationSystem * sys, double * test, double * pval)
{
    int k, df, err = 0;

    if (sys.S == null) {
        return GretlError.badstat;
    }

    k = sys.S.rows;
    df = k * (k - 1) / 2;

    if (sys.method == GretlSystemMethods.sur && sys.iters > 0) {
        /* iterated SUR */
        if (!na(sys.ldet) && sys.diag_test != 0.0) {
            double lr = sys.T * (sys.diag_test - sys.ldet);

            if (test != null) {
                *test = lr;
            }
            if (pval != null) {
                *pval = chisq_cdf_comp(df, lr);
            }
        } else {
            err = GretlError.badstat;
        }
    } else if (sys.diag_test > 0) {
        /* other estimators */
        if (test != null) {
            *test = sys.diag_test;
        }
        if (pval != null) {
            *pval = chisq_cdf_comp(df, sys.diag_test);
        }
    } else {
        err = GretlError.badstat;
    }

    return err;
}

int system_print_sigma(EquationSystem * sys, Prn * prn)
{
    double test, pval;
    int tex = tex_format(prn);
    int k, df;

    if (sys.S == null) {
        return GretlError.data;
    }

    k = sys.S.rows;
    df = k * (k - 1) / 2;

    print_contemp_covariance_matrix(sys.S, sys.ldet, prn);

    if (sys.method == GretlSystemMethods.sur && sys.iters > 0) {
        if (!na(sys.ldet) && sys.diag_test != 0.0) {
            test = sys.T * (sys.diag_test - sys.ldet);
            pval = chisq_cdf_comp(df, test);
            if (tex) {
                pprintf(prn, "%s:\\\\\n".ptr, _("LR test for diagonal covariance matrix").ptr);
                pprintf(prn, "  $\\chi^2(%d)$ = %g [%.4f]".ptr, df, test, pval);
                gretl_prn_newline(prn);
            } else {
                pprintf(prn, "%s:\n".ptr, _("LR test for diagonal covariance matrix").ptr);
                pprintf(prn, "  %s(%d) = %g [%.4f]\n".ptr, _("Chi-square").ptr,
                        df, test, pval);
            }
        }
    } else {
        static const(char) *[4] labels = ["Breusch--Pagan test for diagonal covariance matrix".ptr, "Robust Breusch--Pagan test for diagonal covariance matrix".ptr, "Breusch-Pagan test for diagonal covariance matrix".ptr, "Robust Breusch-Pagan test for diagonal covariance matrix".ptr];
        const(char) * label;

        test = sys.diag_test;
        if (sys.flags & SystemFlags.robust) {
            label = tex ? labels[1] : labels[3];
        } else {
            label = tex ? labels[0] : labels[2];
        }
        if (test > 0) {
            pval = chisq_cdf_comp(df, test);
            if (tex) {
                pprintf(prn, "%s:".ptr, _(label));
                gretl_prn_newline(prn);
                pprintf(prn, "  $\\chi^2(%d)$ = %g [%.4f]".ptr, df, test, pval);
                gretl_prn_newline(prn);
            } else {
                pprintf(prn, "%s:\n".ptr, _(label));
                pprintf(prn, "  %s(%d) = %g [%.4f]\n".ptr, _("Chi-square").ptr,
                        df, test, pval);
            }
        }
    }

    pputc(prn, '\n');

    return 0;
}

enum VarCategory {
    endog,
    exog,
    predet
}

int categorize_variable(int vnum, EquationSystem * sys, int * col, int * lag)
{
    int pos, ret = -1;

    *lag = 0;

    pos = in_gretl_list(sys.ylist, vnum);
    if (pos > 0) {
        ret = VarCategory.endog;
    } else {
        pos = in_gretl_list(sys.xlist, vnum);
        if (pos > 0) {
            ret = VarCategory.exog;
        } else {
            int v = get_predet_parent(sys, vnum, lag);

            pos = in_gretl_list(sys.ylist, v);
            if (pos > 0) {
                ret = VarCategory.predet;
            }
        }
    }

    if (pos > 0) {
        *col = pos - 1;
    }

    return ret;
}

int get_col_and_lag(int vnum, EquationSystem * sys, int * col, int * lag)
{
    int vp, pos, err = 0;

    vp = get_predet_parent(sys, vnum, lag);
    pos = in_gretl_list(sys.ylist, vp);
    if (pos > 0) {
        *col = pos - 1;
    } else {
        *col = 0;
        err = 1;
    }

    return err;
}

/* reduced-form error covariance matrix */

int sys_add_RF_covariance_matrix(EquationSystem * sys, int n)
{
    GretlMatrixRaw * G = null, S = null;
    int err = 0;

    if (!gretl_is_identity_matrix(sys.Gamma)) {
        G = gretl_matrix_copy(sys.Gamma);
        if (G == null) {
            return GretlError.alloc;
        }
    }

    if (n > sys.S.rows) {
        S = gretl_zero_matrix_new(n, n);
        if (S == null) {
            gretl_matrix_free(G);
            return GretlError.alloc;
        }
        gretl_matrix_inscribe_matrix(S, sys.S, 0, 0, GretlMatrixMod.none);
    }

    sys.Sr = gretl_matrix_alloc(n, n);
    if (sys.Sr == null) {
        gretl_matrix_free(G);
        gretl_matrix_free(S);
        return GretlError.alloc;
    }

    if (G != null) {
        err = gretl_SVD_invert_matrix(G);
        if (!err) {
            err = gretl_matrix_qform(G, GretlMatrixMod.none,
                                     (S != null)? S : sys.S,
                                     sys.Sr, GretlMatrixMod.none);
        }
        gretl_matrix_free(G);
    } else {
        gretl_matrix_copy_values(sys.Sr, (S != null)? S : sys.S);
    }

    if (S != null) {
        gretl_matrix_free(S);
    }

    return err;
}

int sys_add_structural_form(EquationSystem * sys, const(DatasetRaw) * dset)
{
    const(int) * ylist = sys.ylist;
    const(int) * xlist = sys.xlist;
    int ne = sys.neqns;
    int ni = sys.nidents;
    int n = ne + ni;
    double x = 0.0;
    int type, col = 0;
    int i, j, vj, lag;
    int err = 0;

    if (ylist[0] > n) {
        gretl_errmsg_set("system: can't add structural form");
        return GretlError.data;
    }

    sys.order = sys_max_predet_lag(sys);

    /* allocate coefficient matrices */

    sys.Gamma = gretl_zero_matrix_new(n, n);
    if (sys.Gamma == null) {
        return GretlError.alloc;
    }

    if (sys.order > 0) {
        sys.A = gretl_zero_matrix_new(n, sys.order * n);
        if (sys.A == null) {
            return GretlError.alloc;
        }
    }

    if (xlist[0] > 0) {
        sys.B = gretl_zero_matrix_new(n, xlist[0]);
        if (sys.B == null) {
            return GretlError.alloc;
        }
    }

    /* process stochastic equations */

    for (i=0; i<sys.neqns && !err; i++) {
        const(Model) * pmod = sys.models[i];
        const(int) * mlist = pmod.list;

        for (j=1; j<=mlist[0] && mlist[j]!=LISTSEP; j++) {
            vj = mlist[j];
            type = categorize_variable(vj, sys, &col, &lag);
            x = (j > 1)? pmod.coeff[j-2] : 1.0;
            if (type == VarCategory.endog) {
                if (j == 1) {
                    /* left-hand side variable */
                    gretl_matrix_set(sys.Gamma, i, col, 1.0);
                } else {
                    gretl_matrix_set(sys.Gamma, i, col, -x);
                }
            } else if (type == VarCategory.exog) {
                gretl_matrix_set(sys.B, i, col, x);
            } else if (type == VarCategory.predet) {
                col += n * (lag - 1);
                gretl_matrix_set(sys.A, i, col, x);
            } else {
                gretl_errmsg_sprintf("structural form:\n couldn't categorize series %s in eqn %d\n".ptr, dset.varname[vj], i+1);
                err = GretlError.data;
            }
        }
    }

    /* process identities */

    for (i=0; i<sys.nidents && !err; i++) {
        const(Identity) * ident = sys.idents[i];

        for (j=0; j<=ident.n_atoms; j++) {
            if (j == 0) {
                vj = ident.depvar;
                x = 1.0;
            } else {
                vj = ident.atoms[j-1].varnum;
                x = (ident.atoms[j-1].op)? -1.0 : 1.0;
            }
            type = categorize_variable(vj, sys, &col, &lag);
            if (type == VarCategory.endog) {
                if (j == 0) {
                    gretl_matrix_set(sys.Gamma, ne+i, col, 1.0);
                } else {
                    gretl_matrix_set(sys.Gamma, ne+i, col, -x);
                }
            } else if (type == VarCategory.exog) {
                gretl_matrix_set(sys.B, ne+i, col, x);
            } else if (type == VarCategory.predet) {
                col += n * (lag - 1);
                gretl_matrix_set(sys.A, ne+i, col, x);
            } else {
                gretl_errmsg_sprintf("structural form:\n couldn't categorize series %s in identity %d\n".ptr, dset.varname[vj], i+1);
                err = GretlError.data;
            }
        }
    }

    if (!err) {
        err = sys_add_RF_covariance_matrix(sys, n);
        if (err) {
            fprintf(stderr, "error %d in sys_add_RF_covariance_matrix\n".ptr, err);
        }
    }

    return err;
}

GretlMatrixRaw * sys_companion_matrix(EquationSystem * sys, int * err)
{
    int m = sys.A.rows;
    int n = sys.A.cols;
    GretlMatrixRaw * C;

    if (m == n) {
        C = sys.A;
    } else {
        C = gretl_zero_matrix_new(n, n);
        if (C == null) {
            *err = GretlError.alloc;
        } else {
            gretl_matrix_inscribe_matrix(C, sys.A, 0, 0,
                                         GretlMatrixMod.none);
            gretl_matrix_inscribe_I(C, m, 0, n - m);
        }
    }

    return C;
}

GretlMatrixRaw * sys_get_fcast_se(EquationSystem * sys, int periods, int * err)
{
    int n = sys.neqns + sys.nidents;
    GretlMatrixRaw * Tmp = null;
    GretlMatrixRaw * V0 = null, Vt = null;
    GretlMatrixRaw * C = null, se = null;
    double vti;
    int i, k, t;

    if (periods <= 0) {
        fprintf(stderr, "Invalid number of periods\n".ptr);
        *err = GretlError.data;
        return null;
    }

    if (sys.A == null) {
        fprintf(stderr, "sys->A is NULL\n".ptr);
        *err = GretlError.data;
        return null;
    }

    C = sys_companion_matrix(sys, err);
    if (*err) {
        return null;
    }

    se = gretl_zero_matrix_new(periods, n);
    if (se == null) {
        *err = GretlError.alloc;
        if (C != sys.A) {
            gretl_matrix_free(C);
        }
        return null;
    }

    k = sys.A.cols;

    Vt = gretl_matrix_alloc(k, k);
    V0 = gretl_zero_matrix_new(k, k);
    Tmp = gretl_matrix_alloc(k, k);

    if (Vt == null || V0 == null || Tmp == null) {
        *err = GretlError.alloc;
        goto bailout;
    }

    for (t=0; t<periods; t++) {
        if (t == 0) {
            /* initial variance */
            gretl_matrix_inscribe_matrix(V0, sys.Sr, 0, 0, GretlMatrixMod.none);
            gretl_matrix_copy_values(Vt, V0);
        } else {
            /* calculate further variances */
            gretl_matrix_copy_values(Tmp, Vt);
            gretl_matrix_qform(C, GretlMatrixMod.none,
                               Tmp, Vt, GretlMatrixMod.none);
            gretl_matrix_add_to(Vt, V0);
        }

        for (i=0; i<n; i++) {
            vti = gretl_matrix_get(Vt, i, i);
            gretl_matrix_set(se, t, i, sqrt(vti));
        }
    }

 bailout:

    gretl_matrix_free(V0);
    gretl_matrix_free(Vt);
    gretl_matrix_free(Tmp);

    if (C != sys.A) {
        gretl_matrix_free(C);
    }

    if (*err) {
        gretl_matrix_free(se);
        se = null;
    }

    return se;
}

int sys_add_fcast_variance(EquationSystem * sys, GretlMatrixRaw * F, int n_static)
{
    GretlMatrixRaw * se = null;
    double ftj, vti;
    int n = sys.neqns + sys.nidents;
    int k = F.rows - n_static;
    int i, j, s;
    int err = 0;

    if (k > 0) {
        se = sys_get_fcast_se(sys, k, &err);
        if (err) {
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
                    if (j < sys.neqns) {
                        vti = sqrt(gretl_matrix_get(sys.Sr, j, j));
                    } else {
                        /* LHS of identity */
                        vti = 0.0;
                    }
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

/* Experimental: generate "fitted values " (in-sample, prior to the
   forecast period) for variables that do not appear on the left-hand
   side of any stochastic equation.  But I'm not quite sure how these
   should be defined.
*/

GretlMatrixRaw * sys_get_fitted_values(EquationSystem * sys, int v, int t1, int t2, const(DatasetRaw) * dset, int * err)
{
    GretlMatrixRaw * F = null;
    GretlMatrixRaw * G = null;
    GretlMatrixRaw * yl = null, x = null;
    GretlMatrixRaw * y = null, yh = null;
    const(int) * ylist = sys.ylist;
    const(int) * xlist = sys.xlist;
    const(int) * plist = sys.plist;
    int n = sys.neqns + sys.nidents;
    double xit;
    int col, T;
    int i, vi, s, t, lag;

    if (sys.Gamma == null) {
        *err = GretlError.data;
        return null;
    }

    if (!gretl_is_identity_matrix(sys.Gamma)) {
        G = gretl_matrix_copy(sys.Gamma);
        if (G == null) {
            *err = GretlError.alloc;
        } else {
            *err = gretl_SVD_invert_matrix(G);
        }
        if (*err) {
            gretl_matrix_free(G);
            return null;
        }
    }

    T = t2 - t1 + 1;

    y = gretl_matrix_alloc(n, 1);
    yh = gretl_matrix_alloc(n, 1);
    F = gretl_zero_matrix_new(T, 1);

    if (y == null || yh == null || F == null) {
        *err = GretlError.alloc;
        goto bailout;
    }

    if (sys.order > 0) {
        yl = gretl_matrix_alloc(sys.order * ylist[0], 1);
        if (yl == null) {
            *err = GretlError.alloc;
            goto bailout;
        }
    }

    if (xlist[0] > 0) {
        x = gretl_matrix_alloc(xlist[0], 1);
        if (x == null) {
            *err = GretlError.alloc;
            goto bailout;
        }
    }

    for (t=t1, s=0; t<=t2; t++, s++) {
        int miss = 0;

        /* lags of endogenous vars */
        if (sys.order > 0) {
            gretl_matrix_zero(yl);
            for (i=1; i<=plist[0] && !miss; i++) {
                vi = plist[i];
                get_col_and_lag(vi, sys, &col, &lag);
                xit = dset.Z[vi][t];
                col += n * (lag - 1);
                if (na(xit)) {
                    miss = 1;
                } else {
                    yl.val[col] = xit;
                }
            }
            if (!miss) {
                gretl_matrix_multiply(sys.A, yl, y);
            }
        } else {
            gretl_matrix_zero(y);
        }

        /* exogenous vars */
        if (xlist[0] > 0 && !miss) {
            for (i=1; i<=xlist[0] && !miss; i++) {
                xit = dset.Z[xlist[i]][t];
                if (na(xit)) {
                    miss = 1;
                } else {
                    x.val[i-1] = xit;
                }
            }
            if (!miss) {
                gretl_matrix_multiply_mod(sys.B, GretlMatrixMod.none,
                                          x, GretlMatrixMod.none,
                                          y, GretlMatrixMod.cumulate);
            }
        }

        if (miss) {
            gretl_vector_set(sys.F, s, NADBL);
        } else {
            /* multiply by Gamma^{-1} */
            if (G != null) {
                gretl_matrix_multiply(G, y, yh);
            } else {
                gretl_matrix_multiply(sys.Gamma, y, yh);
            }
            gretl_vector_set(F, s, yh.val[v]);
        }
    }

 bailout:

    gretl_matrix_free(y);
    gretl_matrix_free(yh);
    gretl_matrix_free(yl);
    gretl_matrix_free(x);
    gretl_matrix_free(G);

    if (*err) {
        gretl_matrix_free(F);
        F = null;
    } else {
        gretl_matrix_set_t1(F, t1);
        gretl_matrix_set_t2(F, t2);
    }

    return F;
}

/* trim the first sys->order rows off sys->F */

int sys_F_trim(EquationSystem * sys)
{
    int err, p = sys.order;
    int r = sys.F.rows - p;
    GretlMatrixRaw * tmp = gretl_matrix_alloc(r, sys.F.cols);

    if (tmp == null) {
        err = GretlError.alloc;
    } else {
        err = gretl_matrix_extract_matrix(tmp, sys.F, p, 0,
                                          GretlMatrixMod.none);
    }

    if (!err) {
        gretl_matrix_free(sys.F);
        sys.F = tmp;
    }

    return err;
}

int sys_add_forecast(EquationSystem * sys, int t1, int t2, const(DatasetRaw) * dset, Gretlopt opt)
{
    GretlMatrixRaw * G = null;
    GretlMatrixRaw * yl = null, x = null;
    GretlMatrixRaw * y = null, yh = null;
    const(int) * ylist = sys.ylist;
    const(int) * xlist = sys.xlist;
    const(int) * plist = sys.plist;
    int n = sys.neqns + sys.nidents;
    double xit, xitd;
    int tdyn, col, T, ncols;
    int i, vi, s, t, lag, s0 = 0;
    int err = 0;

    if (sys.Gamma == null || sys.Gamma.rows != n) {
        fprintf(stderr, "sys_add_forecast: Gamma is broken!\n".ptr);
        return GretlError.data;
    }

    if (!gretl_is_identity_matrix(sys.Gamma)) {
        G = gretl_matrix_copy(sys.Gamma);
        if (G == null) {
            err = GretlError.alloc;
        } else {
            err = gretl_SVD_invert_matrix(G);
        }
        if (err) {
            gretl_matrix_free(G);
            return err;
        }
    }

    /* At which observation (tdyn) does the forecast turn
       dynamic, if at all? */

    if ((opt & Gretlopt.OPT_S) || sys.A == null) {
        /* got the --static option (or no dynamics): never */
        tdyn = t2 + 1;
    } else if (opt & Gretlopt.OPT_D) {
        /* got the --dynamic option: from the start */
        tdyn = t1;
    } else {
        /* by default, for a model with dynamics: just after
           the estimation sample ends */
        tdyn = sys.t2 + 1;
    }

    T = t2 - t1 + 1;
    ncols = 2 * n;

    if (sys.order > 0 && (opt & Gretlopt.OPT_O)) {
        /* out-of-sample option: we'll need some extra rows
           for pre-sample values
        */
        T += sys.order;
    }

    y = gretl_matrix_alloc(n, 1);
    yh = gretl_matrix_alloc(n, 1);
    sys.F = gretl_zero_matrix_new(T, ncols);

    if (y == null || yh == null || sys.F == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    if (sys.order > 0) {
        yl = gretl_matrix_alloc(sys.order * ylist[0], 1);
        if (yl == null) {
            err = GretlError.alloc;
            goto bailout;
        }
        if (opt & Gretlopt.OPT_O) {
            /* out-of-sample */
            int p = sys.order;
            double yti;

            for (t=t1-p, s=0; t<t1; t++, s++) {
                for (i=0; i<n; i++) {
                    yti = dset.Z[ylist[i+1]][t];
                    gretl_matrix_set(sys.F, s, i, yti);
                }
            }
            s0 = p;
        }
    }

    if (xlist[0] > 0) {
        x = gretl_matrix_alloc(xlist[0], 1);
        if (x == null) {
            err = GretlError.alloc;
            goto bailout;
        }
    }

    for (t=t1, s=s0; t<=t2; t++, s++) {
        int miss = 0;

        /* lags of endogenous vars */
        if (sys.order > 0) {
            gretl_matrix_zero(yl);
            for (i=1; i<=plist[0] && !miss; i++) {
                vi = plist[i];
                get_col_and_lag(vi, sys, &col, &lag);
                xitd = NADBL;
                if (t < tdyn || s - lag < 0) {
                    /* pre-forecast value */
                    xit = dset.Z[vi][t];
                } else {
                    /* prior forecast value preferred */
                    if (s - lag >= 0) {
                        xitd = xit = dset.Z[vi][t];
                    }
                    xit = gretl_matrix_get(sys.F, s - lag, col);
                }
                col += n * (lag - 1);
                if (!na(xit)) {
                    yl.val[col] = xit;
                } else if (!na(xitd)) {
                    yl.val[col] = xitd;
                } else {
                    miss = 1;
                }
            }
            if (!miss) {
                gretl_matrix_multiply(sys.A, yl, y);
            }
        } else {
            gretl_matrix_zero(y);
        }

        /* exogenous vars */
        if (xlist[0] > 0 && !miss) {
            for (i=1; i<=xlist[0] && !miss; i++) {
                xit = dset.Z[xlist[i]][t];
                if (na(xit)) {
                    miss = 1;
                } else {
                    x.val[i-1] = xit;
                }
            }
            if (!miss) {
                gretl_matrix_multiply_mod(sys.B, GretlMatrixMod.none,
                                          x, GretlMatrixMod.none,
                                          y, GretlMatrixMod.cumulate);
            }
        }

        if (miss) {
            for (i=0; i<n; i++) {
                gretl_matrix_set(sys.F, s, i, NADBL);
            }
        } else {
            /* multiply by Gamma^{-1} */
            if (G != null) {
                gretl_matrix_multiply(G, y, yh);
            } else {
                gretl_matrix_multiply(sys.Gamma, y, yh);
            }
            for (i=0; i<n; i++) {
                gretl_matrix_set(sys.F, s, i, yh.val[i]);
            }
        }
    }

 bailout:

    gretl_matrix_free(y);
    gretl_matrix_free(yh);
    gretl_matrix_free(yl);
    gretl_matrix_free(x);
    gretl_matrix_free(G);

    if (!err && s0 > 0) {
        err = sys_F_trim(sys);
    }

    if (err) {
        gretl_matrix_free(sys.F);
        sys.F = null;
    } else {
        gretl_matrix_set_t1(sys.F, t1);
        gretl_matrix_set_t2(sys.F, t2);
    }

    if (!err) {
        sys_add_fcast_variance(sys, sys.F, tdyn - t1);
    }

    return err;
}

const(GretlMatrixRaw) * system_get_forecast_matrix(EquationSystem * sys, int t1, int t2, DatasetRaw * dset, Gretlopt opt, int * err)
{
    if (sys.F != null) {
        gretl_matrix_replace(&sys.F, null);
    }

    *err = sys_add_forecast(sys, t1, t2, dset, opt);

    return sys.F;
}

/* attach to the system some extra statistics generated in the course
   of estimation via LIML */

int sys_attach_ldata(EquationSystem * sys)
{
    int i, n = sys.neqns;
    int err = 0;

    sys.ldata = cast(typeof(sys.ldata)) malloc((*sys.ldata).sizeof);

    if (sys.ldata == null) {
        err = GretlError.alloc;
    } else {
        sys.ldata.lmin = null;
        sys.ldata.ll = null;
        sys.ldata.idf = null;

        sys.ldata.lmin = cast(typeof(sys.ldata.lmin)) malloc(n * (*sys.ldata.lmin).sizeof);
        sys.ldata.ll = cast(typeof(sys.ldata.ll)) malloc(n * (*sys.ldata.ll).sizeof);
        sys.ldata.idf = cast(typeof(sys.ldata.idf)) malloc(n * (*sys.ldata.idf).sizeof);

        if (sys.ldata.lmin == null ||
            sys.ldata.ll == null ||
            sys.ldata.idf == null) {
            free(sys.ldata.lmin);
            free(sys.ldata.ll);
            free(sys.ldata.idf);
            free(sys.ldata);
            sys.ldata = null;
            err = GretlError.alloc;
        }
    }

    if (!err) {
        const(Model) * pmod;

        for (i=0; i<n; i++) {
            pmod = sys.models[i];
            sys.ldata.lmin[i] = gretl_model_get_double(pmod, "lmin");
            sys.ldata.ll[i] = pmod.lnL;
            sys.ldata.idf[i] = gretl_model_get_int(pmod, "idf");
        }
    }

    return err;
}

int sys_print_reconstituted_models(EquationSystem * sys, const(DatasetRaw) * dset, Prn * prn)
{
    Model mod;
    const(double) * y;
    double x;
    int print_insts = 0;
    int ifc = 0, nc = 0, ncmax = 0;
    int i, j, t, k = 0;
    int err = 0;

    gretl_model_init(&mod, dset);

    mod.t1 = sys.t1;
    mod.t2 = sys.t2;
    mod.nobs = sys.T;
    mod.aux = ModelAuxCode.sys;
    gretl_model_set_int(&mod, "method", sys.method);

    mod.lnL = NADBL;

    if (sys.method == GretlSystemMethods.tsls ||
        sys.method == GretlSystemMethods.sls3 ||
        sys.method == GretlSystemMethods.fiml ||
        sys.method == GretlSystemMethods.liml) {
        mod.ci = GretlCmdIndex.ivreg;
    } else {
        mod.ci = GretlCmdIndex.ols;
    }

    print_insts = sys.method == GretlSystemMethods.tsls ||
        sys.method == GretlSystemMethods.sls3;

    for (i=0; i<sys.neqns && !err; i++) {
        int * mlist = sys.lists[i];
        int freelist = 0;

        mod.ID = i;

        if (print_insts && !gretl_list_has_separator(mlist)) {
            mod.list = compose_ivreg_list(sys, i);
            if (mod.list == null) {
                err = GretlError.alloc;
            } else {
                freelist = 1;
            }
        } else {
            mod.list = mlist;
        }

        if (!err) {
            ifc = nc = 0;
            for (j=2; j<=mod.list[0]; j++) {
                if (mod.list[j] == 0) {
                    ifc = 1;
                } else if (mod.list[j] == LISTSEP) {
                    break;
                }
                nc++;
            }
        }

        if (!err && nc > ncmax) {
            mod.coeff = cast(typeof(mod.coeff)) realloc(mod.coeff, nc * (*mod.coeff).sizeof);
            mod.sderr = cast(typeof(mod.sderr)) realloc(mod.sderr, nc * (*mod.sderr).sizeof);
            ncmax = nc;
        }

        if (mod.coeff == null || mod.sderr == null) {
            err = GretlError.alloc;
            break;
        }

        mod.ncoeff = nc;
        mod.dfn = nc - ifc;

        if (sys.flags & SystemFlags.dfcorr) {
            gretl_model_set_int(&mod, "dfcorr", 1);
            mod.dfd = mod.nobs - nc;
        } else {
            mod.dfd = mod.nobs;
        }

        y = dset.Z[mod.list[1]];
        mod.ybar = gretl_mean(mod.t1, mod.t2, y);
        mod.sdy = gretl_stddev(mod.t1, mod.t2, y);

        mod.ess = 0.0;
        for (t=0; t<sys.T; t++) {
            x = gretl_matrix_get(sys.E, t, i);
            mod.ess += x * x;
        }

        if (sys.method == GretlSystemMethods.ols ||
            sys.method == GretlSystemMethods.tsls ||
            sys.method == GretlSystemMethods.liml) {
            /* single-equation methods */
            mod.sigma = sqrt(mod.ess / mod.dfd);
        } else {
            mod.sigma = sqrt(gretl_matrix_get(sys.S, i, i));
        }

        for (j=0; j<mod.ncoeff; j++) {
            mod.coeff[j] = sys.b.val[k];
            mod.sderr[j] = sqrt(gretl_matrix_get(sys.vcv, k, k));
            k++;
        }

        if (sys.method == GretlSystemMethods.liml && sys.ldata != null) {
            gretl_model_set_double(&mod, "lmin", sys.ldata.lmin[i]);
            mod.lnL = sys.ldata.ll[i];
            gretl_model_set_int(&mod, "idf", sys.ldata.idf[i]);
        }

        printmodel(&mod, dset, Gretlopt.OPT_NONE, prn);

        if (freelist) {
            free(mod.list);
        }

        mod.list = null;
    }

    clear_model(&mod);

    return err;
}

int gretl_system_print(EquationSystem * sys, const(DatasetRaw) * dset, Gretlopt opt, Prn * prn)
{
    const(char) * name = sys.name;
    char * sysstr = null;
    int tex = tex_format(prn);
    int nr = system_n_restrictions(sys);
    int i;

    if (sys.models != null &&
        sys.method == GretlSystemMethods.liml &&
        sys.ldata == null) {
        sys_attach_ldata(sys);
    }

    if (name != null && sys_anonymous(name)) {
        /* don't print internal reference name */
        name = null;
    }

    if (tex) {
        pputs(prn, "\\begin{center}\n");
    } else {
    pputc(prn, '\n');
    }

    sysstr = system_get_full_string(sys);
    if (name != null) {
    if (tex) {
        pprintf(prn, "%s, %s\\\\\n".ptr, _("Equation system").ptr, name);
    } else {
        pprintf(prn, "%s, %s\n".ptr, _("Equation system").ptr, name);
    }
    pprintf(prn, "%s: %s".ptr, _("Estimator").ptr, sysstr);
    } else {
    pprintf(prn, "%s, %s".ptr, _("Equation system").ptr, sysstr);
    }
    g_free(sysstr);

    if (sys.iters > 0) {
        gretl_prn_newline(prn);
    pprintf(prn, _("Convergence achieved after %d iterations\n").ptr, sys.iters);
        if (sys.method == GretlSystemMethods.sur ||
            sys.method == GretlSystemMethods.fiml) {
            if (tex) {
                gretl_prn_newline(prn);
                pprintf(prn, "%s = ".ptr, _("Log-likelihood").ptr);
                if (sys.ll < 0) {
                    pprintf(prn, "$-$%g".ptr, -sys.ll);
                } else {
                    pprintf(prn, "%g".ptr, sys.ll);
                }
            } else {
                pprintf(prn, "%s = %g\n".ptr, _("Log-likelihood").ptr, sys.ll);
            }
        }
    }

    if (tex) {
        pputs(prn, "\n\\end{center}\n\n");
    } else {
        pputc(prn, '\n');
    }

    if (sys.models != null) {
        for (i=0; i<sys.neqns; i++) {
            if (sys.flags & SystemFlags.dfcorr) {
                gretl_model_set_int(sys.models[i], "dfcorr", 1);
            }
            printmodel(sys.models[i], dset, Gretlopt.OPT_NONE, prn);
        }
    } else {
        sys_print_reconstituted_models(sys, dset, prn);
    }

    system_print_sigma(sys, prn);

    if (nr == 0 && (sys.method == GretlSystemMethods.fiml ||
                    sys.method == GretlSystemMethods.sls3 ||
                    sys.method == GretlSystemMethods.sur)) {
        print_system_overid_test(sys, prn);
    }

    return 0;
}

void ensure_asy_printout(EquationSystem * sys)
{
    int i;

    if (sys.models != null) {
        for (i=0; i<sys.neqns; i++) {
            sys.models[i].ci = GretlCmdIndex.ivreg;
        }
    }
}

int system_save_and_print_results(EquationSystem * sys, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    int nr = system_n_restrictions(sys);
    int err = 0;

    if (sys.E != null) {
        gretl_matrix_set_t1(sys.E, sys.t1);
        gretl_matrix_set_t2(sys.E, sys.t2);
    }

    if (sys.iters > 0 && sys.method == GretlSystemMethods.sur && nr == 0) {
        sur_ols_diag(sys);
    }

    sys.ldet = gretl_vcv_log_determinant(sys.S, &err);

    if (!err) {
        err = system_add_yhat_matrix(sys);
    }

    if (!err) {
        err = sys_add_structural_form(sys, dset);
    }

    if (!(opt & Gretlopt.OPT_Q)) {
        if (sys.method == GretlSystemMethods.fiml) {
            ensure_asy_printout(sys);
        }
        gretl_system_print(sys, dset, opt, prn);
    }

    return err;
}

int system_autocorrelation_test(EquationSystem * sys, int order, Gretlopt opt, Prn * prn)
{
    double * u;
    double lb;
    int i, err = 0;

    for (i=0; i<sys.neqns && !err; i++) {
        if (!(opt & Gretlopt.OPT_Q)) {
            pprintf(prn, "%s %d:\n".ptr, _("Equation").ptr, i + 1);
        }
        u = sys.E.val + (i * sys.T);
        lb = ljung_box(order, 0, sys.T - 1, u, &err);
        if (!err && !((opt & Gretlopt.OPT_Q))) {
            pprintf(prn, "%s: %s(%d) = %g [%.4f]\n\n".ptr,
                    _("Ljung-Box Q'").ptr, _("Chi-square").ptr, order,
                    lb, chisq_cdf_comp(order, lb));
        }
    }

    return err;
}

int system_arch_test(EquationSystem * sys, int order, Gretlopt opt, Prn * prn)
{
    const(double) * u;
    int i, err = 0;

    if (!(opt & Gretlopt.OPT_I)) {
        pputc(prn, '\n');
    }

    for (i=0; i<sys.neqns && !err; i++) {
        if (!(opt & Gretlopt.OPT_I)) {
            pprintf(prn, "%s %d:\n".ptr, _("Equation").ptr, i + 1);
        }
        u = sys.E.val + (i * sys.T);
        err = array_arch_test(u, sys.T, order, opt | Gretlopt.OPT_M, prn);
    }

    return err;
}

int system_supports_method(EquationSystem * sys, int method)
{
    if (sys != null) {
        int i;

        for (i=0; i<sys.neqns; i++) {
            if (gretl_list_has_separator(sys.lists[i])) {
                return method == GretlSystemMethods.tsls ||
                    method == GretlSystemMethods.sls3;
            }
        }
    }

    return 1;
}

void finalize_liml_model(Model * pmod, int * reglist, EquationSystem * sys)
{
    *pmod = *sys.models[0];

    gretl_model_destroy_data_item(pmod, "tslsX");
    gretl_model_destroy_data_item(pmod, "endog");
    gretl_model_destroy_data_item(pmod, "method");
    gretl_model_destroy_data_item(pmod, "liml_y");

    free(sys.models[0]);
    free(sys.models);
    sys.models = null;

    pmod.xlist = gretl_list_sublist(reglist, 2, -1);
    pmod.aux = ModelAuxCode.none;
    pmod.opt |= Gretlopt.OPT_L;
    pmod.rsq = pmod.adjrsq = NADBL;
    pmod.fstt = NADBL;
    set_model_id(pmod, Gretlopt.OPT_NONE);
}

/* implement the --liml option to "tsls" */

Model single_equation_liml(const(int) * list, DatasetRaw * dset, Gretlopt opt)
{
    int * mlist = null, ilist = null;
    EquationSystem * sys = null;
    Model model;
    int err = 0;

    gretl_model_init(&model, dset);

    err = ivreg_process_lists(list, &mlist, &ilist);

    if (!err) {
        sys = equation_system_new(GretlSystemMethods.liml, null, &err);
    }
    if (!err) {
        err = equation_system_append(sys, mlist);
    }
    if (!err) {
        sys.ilist = ilist;
        err = equation_system_finalize(sys, dset, Gretlopt.OPT_S, null);
    }

    if (err) {
        model.errcode = err;
    } else {
        finalize_liml_model(&model, mlist, sys);
    }

    equation_system_destroy(sys);
    free(mlist);

    return model;
}

int gretl_system_get_sample(EquationSystem * sys, int * t1, int * t2)
{
    if (sys != null) {
        *t1 = sys.smpl_t1;
        *t2 = sys.smpl_t2;
        return 0;
    }

    return GretlError.data;
}

