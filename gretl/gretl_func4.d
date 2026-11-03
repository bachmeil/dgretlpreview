/* Port of gretl_func.c, part 4 of 4: see gretl_func0.d. */

module gretl.gretl_func4;

import gretl.glist : g_list_previous;
import gretl.gretl_func_ext : flow_control;

import gretl.interacth : ExecFlags;

import gretl.cmd_privateh : LoopSet;
import gretl.gretl_func_ext : get_loop_renaming;
import gretl.gretl_func_ext : get_gui_callback, gretl_cmd_init, gretl_exec_state_clear, gretl_exec_state_init, gretl_if_state_check, gretl_if_state_false, gretl_if_state_record, gretl_if_state_reset, gretl_loop_exec, is_plotting_command, loop_reset_uvars, maybe_exec_line, md_to_gretl, parse_command_line;

import gretl.cmd_privateh : ExecCallback;
import gretl.dataset : series_get_display_name;
import gretl.dataseth : VarFlags;
import gretl.gretl_typemap : gretl_type_get_name;
import gretl.uservar : user_var_get_value_and_type;
import gretl.versionh : GRETL_VERSION;

import gretl.cmd_privateh : cmd_subst;
import gretl.dataset : dataset_add_series_as, dataset_copy_series_as, dataset_drop_last_variables, dataset_drop_observations, dataset_drop_variable, is_string_valued, series_decrement_stack_level, series_get_label, series_get_stack_level, series_get_string_table, series_is_listarg, series_set_display_name, series_set_label, series_set_stack_level, series_unset_flag;
import gretl.dataseth : dataset_is_panel, dataset_is_time_series, quarterly_or_monthly;
import gretl.gretl_func3 : arg_get_data;
import gretl.geneval5 : genr_reset_uvars;
import gretl.genmain : current_series_index, execute_genr, generate, genr_compile, series_index;
import gretl.glib : gchar;
import gretl.glist : GList, g_list_append, g_list_last;
import gretl.gmem : g_free;
import gretl.gretl_array : get_array_by_name, gretl_array_copy, gretl_array_get_data, gretl_array_get_length, gretl_array_get_type, gretl_array_pull_from_stack, gretl_array_set_data, gretl_array_set_type;
import gretl.gretl_bundle : get_bundle_by_name, gretl_bundle_copy, gretl_bundle_pull_from_stack, gretl_bundle_set_creator;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_errors : errmsg_get_with_default, gretl_errmsg_get, gretl_errmsg_set, gretl_errmsg_sprintf;
import gretl.gretl_func1 : current_function_call, current_function_info, currently_called_function, get_user_function_by_name, gretl_function_depth, gretl_function_recursing, maybe_destroy_fncall, no_scalar_default, set_executing_off, set_listargs_from_call, ufunc_free;
import gretl.gretl_func2 : function_package_free, function_package_unload_full_by_filename, is_pdf_ref;
import gretl.gretl_func3 : add_empty_list, help_text_is_markdown, localize_bundled_lists, localize_list, localize_object_as_shell;
import gretl.gretl_list : gretl_list_append_term, gretl_list_copy, gretl_list_delete_at_pos, in_gretl_list;
import gretl.gretl_listh : LISTSEP;
import gretl.gretl_matrix : gretl_matrix_alloc, gretl_matrix_copy, gretl_matrix_free;
import gretl.gretl_model : allocate_working_model;
import gretl.gretl_paths : get_full_filename, gretl_deltree, gretl_file_exists, gretl_remove;
import gretl.gretl_prn : pprintf, pputc, pputs, print_redirection_level;
import gretl.gretl_string_table : gretl_insert_builtin_string, series_table_copy;
import gretl.gretl_typemap : gretl_is_series_type, gretl_type_get_plain_type, gretl_type_get_ref_type;
import gretl.gretl_typemaph : NUMERIC_TYPE;
import gretl.gretl_utils : copyvec, gretl_int_from_double, gretl_version_number, gretl_version_string;
import gretl.gstrfuncs : g_strdup, g_strdup_printf;
import gretl.libset : gretl_debugging_on, gretl_echo_on, gretl_iterating, gretl_messages_on, pop_program_state, push_program_state, set_gretl_echo, set_gretl_messages;
import gretl.missingh : NADBL, na;
import gretl.printout : bufgets_finalize, bufgets_init, bufspace, safe_bufgets;
import gretl.strutils : gretl_namechar_spn, gretl_strdup, gretl_strndup, has_suffix, string_is_blank, strrslash, switch_ext_new;
import gretl.subsample : complex_subsampled, copy_dataset_submask, restore_full_sample, restrict_sample_from_mask, submask_cmp;
import gretl.system : delete_anonymous_equation_system;
import gretl.usermat : get_matrix_by_name, steal_matrix_by_name;
import gretl.uservar : arg_add_as_shell, copy_as_arg, destroy_user_vars_at_level, get_list_by_name, get_string_by_name, get_user_var_by_name, gretl_is_scalar, gretl_scalar_get_value, set_previous_depth, switch_uservar_hash, user_var_adjust_level, user_var_delete_by_name, user_var_get_type, user_var_get_value, user_var_get_value_by_name, user_var_localize, user_var_set_name;
import gretl.uservarh : UserVar;

import gretl.cmd_privateh : Cmd;

import gretl.cmd_privateh : ExecState, Stmt;
import gretl.genmainh : Generator;
import gretl.gretl_string_table : SeriesTable;
import gretl.objstackh : GretlObjType;

import core.stdc.math : fabs, floor, isnan, log, sqrt, pow, exp;
import core.stdc.stdlib : free, malloc, calloc, realloc, qsort, atoi, atof;
import core.stdc.stdio : fprintf, printf, sprintf, snprintf, stderr, stdout, FILE, fopen, fclose, fgets, fputs;
import core.stdc.string : memcpy, memset, strcmp, strlen, strcpy, strncpy, strcat, strchr, strstr, strncmp, strrchr;

import core.stdc.ctype : isalpha, isspace, isdigit;
import core.stdc.errno : errno;
import core.stdc.stdio : fputc, rewind, sscanf, remove;
import core.stdc.stdlib : strtod, strtol;
import core.stdc.string : memmove, strcspn, strncat, strspn;
version (Posix) import core.sys.posix.sys.stat : stat_t; else import gretl.winport : stat_t;
import gretl.gretl_func0;
import gretl.xmlh : xmlChar;
import gretl.gretl_funch;
import gretl.libgretlh;
import gretl.gretl_matrixh;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_prn : Prn;

private alias XUC = const(xmlChar) *;

int localize_const_object(Fncall * call, int i, FnParam * fp)
{
    FnArg * arg = &call.args[i];
    void *data = arg_get_data(arg, GretlType.none);
    int err = 0;

    if (data == null) {
        err = GretlError.data;
    } else if (arg.uvar == null) {
        /* the const argument is an anonymous object */
        err = arg_add_as_shell(fp.name, arg.type, data);
    } else {
        /* a named object: in the simplest case we just
           adjust the level and name of the object for
           the duration of the function call, but note
           that we can do this only once for any given
           object.
        */
        UserVar * thisvar = arg.uvar;
        int j;
        int done = 0;

        for (j=0; j<i; j++) {
            if (call.args[j].uvar == thisvar) {
                /* we've already localized this one! */
                err = arg_add_as_shell(fp.name, arg.type, data);
                done = 1;
                break;
            }
        }
        if (!done) {
            user_var_adjust_level(thisvar, 1);
            user_var_set_name(thisvar, fp.name);
            arg.shifted = 1;
        }
    }

    return err;
}

int localize_series_ref(Fncall * call, FnArg * arg, FnParam * fp, DatasetRaw * dset)
{
    int level = fn_executing + 1;
    int v = arg.val.idnum;

    series_set_stack_level(dset, v, level);
    strcpy(dset.varname[v], fp.name);

    if (!in_gretl_list(call.ptrvars, v)) {
        gretl_list_append_term(&call.ptrvars, v);
    }

    return 0;
}

int argval_get_int(FnParam * param, double x, int * err)
{
    int ret = gretl_int_from_double(x, err);

    if (*err) {
        gretl_errmsg_sprintf(_("%s: expected an integer but found %g").ptr,
                             param.name, x);
    }

    return ret;
}

int real_add_scalar_arg(FnParam * param, double x)
{
    int err = 0;

    if (na(x)) {
        /* always allow NA for scalar args (?) */
        err = copy_as_arg(param.name, GretlType.double_, &x);
    } else {
        if (param.type == GretlType.bool_) {
            if (x != 0.0) {
                x = 1.0;
            }
        } else if (param.type == GretlType.int_ ||
                   param.type == GretlType.obs) {
            x = argval_get_int(param, x, &err);
        }

        if (!err) {
            if ((!na(param.min) && x < param.min) ||
                (!na(param.max) && x > param.max)) {
                gretl_errmsg_sprintf(_("%s: argument value %g is out of bounds").ptr,
                                     param.name, x);
                err = GretlError.data;
            } else {
                err = copy_as_arg(param.name, GretlType.double_, &x);
            }
        }
    }

    return err;
}

/* Scalar function arguments only: if the arg is not supplied, use the
   default that is found in the function specification, if any.
*/

int add_scalar_arg_default(FnParam * param)
{
    if (default_unset(param)) {
        /* should be impossible here, but... */
        return GretlError.data;
    } else {
        return real_add_scalar_arg(param, param.deflt);
    }
}

/* Handle the case of a scalar parameter for which a 1 x 1 matrix
   was given as argument.
*/

int do_matrix_scalar_cast(FnArg * arg, FnParam * param)
{
    GretlMatrixRaw * m = arg.val.m;

    return real_add_scalar_arg(param, m.val[0]);
}

/* Handle the case of a matrix parameter for which a scalar
   was given as argument.
*/

int do_scalar_matrix_cast(FnArg * arg, FnParam * param)
{
    GretlMatrixRaw *m = gretl_matrix_alloc(1, 1);
    int err = 0;

    if (m == null) {
        err = GretlError.alloc;
    } else {
        m.val[0] = arg.val.x;
        err = copy_as_arg(param.name, GretlType.matrix, m);
        gretl_matrix_free(m);
    }

    return err;
}

void fncall_finalize_listvars(Fncall * call)
{
    int i;
    int v;

    for (i=call.listvars[0]; i>0; i--) {
        v = call.listvars[i];
        if (in_gretl_list(call.ptrvars, v)) {
            gretl_list_delete_at_pos(call.listvars, i);
        }
    }

    if (call.listvars[0] == 0) {
        free(call.listvars);
        call.listvars = null;
    }
}

int upnames_match(FnArg * ai, FnArg * aj)
{
    return ai.upname != null && aj.upname != null &&
        strcmp(ai.upname, aj.upname) == 0;
}

/* Note: a caller cannot be allowed to supply a given variable in
   pointer form for more than one argument slot in a function call
   (although ordinary arguments may be repeated).
*/

int duplicated_pointer_arg_check(Fncall * call)
{
    FnArg * ai;
    FnArg * aj;
    int i;
    int j;
    int err = 0;

    for (i=0; i<call.argc && !err; i++) {
        ai = &call.args[i];
        if (gretl_ref_type(ai.type)) {
            for (j=i+1; j<call.argc && !err; j++) {
                aj = &call.args[j];
                if (gretl_ref_type(aj.type) && upnames_match(ai, aj)) {
                    gretl_errmsg_set(_("Duplicated pointer argument: not allowed").ptr);
                    err = GretlError.data;
                }
            }
        }
    }

    return err;
}

int process_series_arg(Fncall * call, FnParam * fp, FnArg * arg, DatasetRaw * dset)
{
    if (arg.type == GretlType.useries) {
        /* an existing named series */
        if (param_is_const(fp) && arg.val.idnum != 0) {
            /* we can pass it by reference, if the series
	       in question is not the one named "const",
	       with ID number 0
	    */
            return localize_series_ref(call, arg, fp, dset);
        } else {
            /* pass it by value */
            return dataset_copy_series_as(dset, arg.val.idnum, fp.name);
        }
    } else if (arg.val.px != null) {
        /* an on-the-fly series */
        return dataset_add_series_as(dset, arg.val.px, fp.name);
    } else {
        return GretlError.data;
    }
}

/* for matrix, bundle, array, string */

int process_object_arg(Fncall * call, int i, FnParam * fp)
{
    FnArg * arg = &call.args[i];

    if (param_is_const(fp)) {
        /* we can pass it by reference */
        return localize_const_object(call, i, fp);
    } else {
        /* pass it by value */
        return copy_as_arg(fp.name, arg.type,
                           arg_get_data(arg, GretlType.none));
    }
}

/* for *matrix, *bundle, *array, *string */

int process_object_ref_arg(FnArg * arg, FnParam * fp)
{
    if (arg.upname != null) {
        return user_var_localize(arg.upname, fp.name, fp.type);
    } else {
        return localize_object_as_shell(arg, fp);
    }
}

/* Note that if we reach here we've already successfully negotiated
   check_function_args(): now we're actually making the argument
   objects (if any) available within the function.
*/

int allocate_function_args(Fncall * call, DatasetRaw * dset)
{
    Ufunc * fun = call.fun;
    FnArg * arg;
    FnParam * fp;
    int i;
    int err;

    err = duplicated_pointer_arg_check(call);

    for (i=0; i<call.argc && !err; i++) {
        arg = &call.args[i];
        fp = &fun.params[i];

        if (arg.upname != null && object_is_const(arg.upname, -1) &&
            gretl_ref_type(fp.type) && !param_is_const(fp)) {
            const(char) * caller = null;

            current_function_info(&caller, null);
            gretl_errmsg_sprintf(_("%s() tries to pass const argument " ~
		"%s to %s() in pointer form with no\n" ~
		"const guarantee\n").ptr,
                                 caller, fp.name, call.fun.name);
            err = GretlError.invarg;
            break;
        }

        if (arg.type == GretlType.none) {
            if (gretl_scalar_type(fp.type)) {
                err = add_scalar_arg_default(fp);
            } else if (fp.type == GretlType.list) {
                err = add_empty_list(call, fp, dset);
            }
            continue;
        }

        if (fp.type == GretlType.numeric) {
            /* param supports overloading */
            if (gretl_scalar_type(arg.type)) {
                err = real_add_scalar_arg(fp, arg.val.x);
            } else if (gretl_is_series_type(arg.type)) {
                process_series_arg(call, fp, arg, dset);
            } else if (arg.type == GretlType.matrix) {
                err = process_object_arg(call, i, fp);
            } else if (arg.type == GretlType.list) {
                err = localize_list(call, arg, fp, dset);
            }
            continue;
        }

        if (gretl_scalar_type(fp.type)) {
            if (arg.type == GretlType.matrix) {
                err = do_matrix_scalar_cast(arg, fp);
            } else {
                err = real_add_scalar_arg(fp, arg.val.x);
            }
        } else if (fp.type == GretlType.series) {
            err = process_series_arg(call, fp, arg, dset);
        } else if (fp.type == GretlType.matrix &&
                   arg.type == GretlType.double_) {
            err = do_scalar_matrix_cast(arg, fp);
        } else if (fp.type == GretlType.matrix ||
                   fp.type == GretlType.bundle ||
                   fp.type == GretlType.string_ ||
                   gretl_array_type(fp.type)) {
            err = process_object_arg(call, i, fp);
        } else if (fp.type == GretlType.list) {
            err = localize_list(call, arg, fp, dset);
        } else if (fp.type == GretlType.series_ref) {
            err = localize_series_ref(call, arg, fp, dset);
        } else if (gretl_ref_type(fp.type)) {
            err = process_object_ref_arg(arg, fp);
        }
        if (!err && (fp.type == GretlType.bundle ||
                     fp.type == GretlType.bundle_ref)) {
            err = localize_bundled_lists(call, arg, fp, dset);
        }
    }

    /* now for any trailing parameters without matching arguments */

    for (i=call.argc; i<fun.n_params && !err; i++) {
        fp = &fun.params[i];
        if (gretl_scalar_type(fp.type)) {
            err = add_scalar_arg_default(fp);
        } else if (fp.type == GretlType.list) {
            err = add_empty_list(call, fp, dset);
        }
    }

    if (call.listvars != null) {
        if (err) {
            free(call.listvars);
            call.listvars = null;
        } else {
            fncall_finalize_listvars(call);
        }
    }

    if (!err) {
        set_listargs_from_call(call, dset);
    }


    return err;
}

/**
 * check_function_needs:
 * @dset: pointer to dataset info.
 * @dreq: function data requirements flag.
 * @minver: function minimum program version requirement.
 * @pkg: the package to which the function belongs.
 *
 * Checks whether the requirements in @dreq and @minver
 * are jointly satisfied by the current dataset and gretl
 * version.
 *
 * Returns: 0 if all is OK; GretlError.data if the current dataset
 * (or lack thereof) does not satisfy @dreq; 1 otherwise.
 */

extern (C) int check_function_needs(const(DatasetRaw) * dset, DataReq dreq, int minver, Fnpkg * pkg)
{
    static __gshared int thisver = 0;

    if (thisver == 0) {
        thisver = gretl_version_number(GRETL_VERSION);
    }

    if (minver > thisver) {
        char[8] vstr;

        gretl_version_string(vstr.ptr, minver);
        if (pkg != null) {
            gretl_errmsg_sprintf(_("The package %s needs gretl version %s").ptr,
                                 pkg.name, vstr.ptr);
        }
        return 1;
    }

    if ((dset == null || dset.v == 0) && dreq != DataReq.nodata_ok) {
        if (pkg != null) {
            gretl_errmsg_sprintf(_("The package %s needs a dataset in place").ptr,
                                 pkg.name);
        }
        return GretlError.data;
    }

    if (dreq == DataReq.needs_ts && !dataset_is_time_series(dset)) {
        if (pkg != null) {
            gretl_errmsg_sprintf(_("The package %s needs time-series data").ptr,
                                 pkg.name);
        }
        return GretlError.data;
    }

    if (dreq == DataReq.needs_panel && !dataset_is_panel(dset)) {
        if (pkg != null) {
            gretl_errmsg_sprintf(_("The package %s needs panel data").ptr,
                                 pkg.name);
        }
        return GretlError.data;
    }

    if (dreq == DataReq.needs_qm && !quarterly_or_monthly(dset)) {
         if (pkg != null) {
            gretl_errmsg_sprintf(_("The package %s needs quarterly or monthly data").ptr,
                                 pkg.name);
        }
        return GretlError.data;
    }

    if (pkg != null) {
        pkg.prechecked = 1;
    }

    return 0;
}

/**
 * package_version_ok:
 * @minver: function package minimum program version requirement.
 * @reqstr: location to write required version string (should be
 * at least 8 bytes), or null.
 *
 * Returns: 1 if the running version of gretl satisfies the
 * minimum version requirement in @minver; otherwise returns 0,
 * in which case the string representation of @minver is written
 * @reqstr if this is non-null.
 */

extern (C) int package_version_ok(int minver, char * reqstr)
{
    static __gshared int thisver = 0;
    int ret = 0;

    if (thisver == 0) {
        thisver = gretl_version_number(GRETL_VERSION);
    }

    ret = thisver >= minver;

    if (!ret && reqstr != null) {
        gretl_version_string(reqstr, minver);
    }

    return ret;
}

int maybe_check_function_needs(const(DatasetRaw) * dset, Ufunc * fun)
{
    if (fun.pkg == null || fun.pkg.prechecked) {
        return 0;
    }


    return check_function_needs(dset, fun.pkg.dreq,
                                fun.pkg.minver, fun.pkg);
}

/* next block: handling function return values */

int handle_scalar_return(Fncall * call, void * ptr, GretlType rtype)
{
    static __gshared double xret;
    const(char) * vname = call.retname;
    int err = 0;

    if (gretl_is_scalar(vname)) {
        xret = gretl_scalar_get_value(vname, null);
    } else {
        GretlMatrixRaw *m = get_matrix_by_name(vname);

        if (gretl_matrix_is_scalar(m)) {
            xret = m.val[0];
        } else {
            UserVar *uv = get_user_var_by_name(vname);
            GretlType t = user_var_get_type(uv);

            gretl_errmsg_sprintf(_("Function %s did not provide the specified return value\n" ~
		"(expected %s, got %s)").ptr, call.fun.name,
                                 gretl_type_get_name(rtype), gretl_type_get_name(t));
            err = GretlError.types;
            xret = NADBL;
        }
    }

    *cast(double **) ptr = &xret;

    return err;
}

int handle_series_return(const(char) * vname, void * ptr, DatasetRaw * dset, int copy, char ** label, SeriesTable ** stab)
{
    int v = series_index(dset, vname);
    double * x = null;
    int err = 0;

    if (!copy && v == 0) {
        copy = 1;
    }

    if (v >= 0 && v < dset.v) {
        if (copy) {
            x = copyvec(dset.Z[v], dset.n);
            if (x == null) {
                err = GretlError.alloc;
            }
        } else {
            x = dset.Z[v];
            dset.Z[v] = null;
        }
    } else {
        err = GretlError.unkvar;
    }

    *cast(double **) ptr = x;

    if (!err) {
        if (label != null) {
            const(char) * vstr = series_get_label(dset, v);

            if (vstr != null) {
                *label = gretl_strdup(vstr);
            } else {
                *label = gretl_strdup("");
            }
        }
        if (stab != null && is_string_valued(dset, v)) {
            SeriesTable *st = series_get_string_table(dset, v);

            *stab = series_table_copy(st);
        }
    }

    return err;
}

int handle_matrix_return(const(char) * name, void * ptr, int copy)
{
    GretlMatrixRaw * ret = null;
    int err = 0;

    if (copy) {
        GretlMatrixRaw *m = get_matrix_by_name(name);

        if (m != null) {
            ret = gretl_matrix_copy(m);
            if (ret == null) {
                err = GretlError.alloc;
            }
        }
    } else {
        ret = steal_matrix_by_name(name);
    }

    if (ret == null && !err) {
        err = GretlError.unkvar;
    }

    *cast(GretlMatrixRaw **) ptr = ret;

    return err;
}

GretlType fix_numeric_return_type(Fncall * call, DatasetRaw * dset)
{
    GretlType t = GretlType.none;
    void * uv;

    uv = user_var_get_value_and_type(call.retname, &t);
    if (uv != null) {
        call.rtype = t;
    } else {
        int v = current_series_index(dset, call.retname);

        if (v >= 0) {
            t = call.rtype = GretlType.series;
        }
    }

    return t;
}

int handle_bundle_return(Fncall * call, void * ptr, int copy)
{
    const(char) * name = call.retname;
    GretlBundle * ret = null;
    int err = 0;

    if (copy) {
        GretlBundle *b = get_bundle_by_name(name);

        if (b != null) {
            ret = gretl_bundle_copy(b, &err);
        }
    } else {
        ret = gretl_bundle_pull_from_stack(name, &err);
    }

    if (ret != null && call.fun.pkg != null &&
        gretl_function_depth() == 1) {
        if (call.fun.pkg_role != UfunRole.bundle_fcast) {
            gretl_bundle_set_creator(ret, call.fun.pkg.name);
        }
    }

    *cast(GretlBundle **) ptr = ret;

    return err;
}

int handle_array_return(Fncall * call, void * ptr, int copy)
{
    const(char) * name = call.retname;
    GretlArray * ret = null;
    int err = 0;

    if (copy) {
        GretlArray *a = get_array_by_name(name);

        if (a != null) {
            ret = gretl_array_copy(a, &err);
        }
    } else {
        ret = gretl_array_pull_from_stack(name, &err);
    }

    if (ret != null) {
	GretlType atype = gretl_array_get_type(ret);

	if (atype != call.rtype) {
	    if (atype == GretlType.any) {
		gretl_array_set_type(ret, call.rtype);
	    } else {
		err = GretlError.types;
	    }
	}
    }

    *cast(GretlArray **) ptr = ret;

    return err;
}

void replace_caller_series(int targ, int src, DatasetRaw * dset)
{
    const(char) * vlabel = series_get_label(dset, src);
    int t;

    /* replace data values */
    for (t=dset.t1; t<=dset.t2; t++) {
        dset.Z[targ][t] = dset.Z[src][t];
    }

    /* replace variable info? */
    series_set_label(dset, targ, vlabel == null ? "" : vlabel);
    series_set_display_name(dset, targ, series_get_display_name(dset, src));
}

/* Deal with a list, named @lname, that exists at the level of a
   function whose execution is now terminating.  Note that this list
   may be the direct return value of the function (in which case @ret
   will be non-null), or it may have been passed to the function as an
   argument (in which case @arg will be non-null), or it may have been
   constructed on the fly.

   If the list was given as a function argument we simply shunt all
   its members to the caller's stack level.  But if it's the direct
   return value we need to overwrite any variables at caller level
   that have been redefined within the function.  If any series have
   been redefined in that way we also need to adjust the list itself,
   replacing the ID numbers of local variables with the caller-level
   IDs -- all supposing the direct return value is being assigned by
   the caller.

   If the list was a temporary construction (a list parameter was
   required, but a single series was given as argument), we just
   destroy it.
*/

int unlocalize_list(Fncall * call, const(char) * lname, FnArg * arg, void * ret, DatasetRaw * dset)
{
    int *list = get_list_by_name(lname);
    int d = gretl_function_depth();
    int upd = d - 1;
    int i;
    int vi;


    if (list == null) {
        return GretlError.data;
    }

    if (arg == null && ret == null) {
        ; /* no-op: we'll trash the list later */
    } else if (arg == null) {
        /* handle list as direct return value */
        int j;
        int lev;
        int overwrite;
        const(char) * vname;

        for (i=1; i<=list[0]; i++) {
            overwrite = 0;
            vi = list[i];
            vname = dset.varname[vi];
            if (in_gretl_list(call.ptrvars, vi)) {
                /* 2021-06-28: detect a somewhat anomalous case */
                fprintf(stderr, "*** return list '%s' contains series %s, passed in " ~
		"pointer form ***\n", lname, vname);
                continue;
            }
            series_unset_flag(dset, vi, VarFlags.listarg);
            if (vi > 0 && vi < dset.v && series_get_stack_level(dset, vi) == d) {
                for (j=1; j<dset.v; j++) {
                    lev = series_get_stack_level(dset, j);
                    if (lev == upd && !strcmp(dset.varname[j], vname)) {
                        overwrite = 1;
                        break;
                    }
                    if (lev == d && j < vi && series_is_listarg(dset, j, null) &&
                        !strcmp(dset.varname[j], vname)) {
                        overwrite = 1;
                        break;
                    }
                }
                if (overwrite) {
                    replace_caller_series(j, vi, dset);
                    /* replace ID number in list */
                    list[i] = j;
                } else {
                    series_set_stack_level(dset, vi, upd);
                }
            }
        }
    } else {
        /* list was given as argument to function */
        for (i=1; i<=list[0]; i++) {
            vi = list[i];
            if (vi == LISTSEP) {
                continue;
            }
            if (series_is_listarg(dset, vi, null)) {
                series_unset_flag(dset, vi, VarFlags.listarg);
                series_set_stack_level(dset, vi, upd);
            }
        }
        if (arg.type != GretlType.list) {
            /* the list was constructed on the fly */
            user_var_delete_by_name(lname, null);
        }
    }

    return 0;
}

int handle_string_return(const(char) * sname, void * ptr)
{
    const(char) * s = get_string_by_name(sname);
    char * ret = null;
    int err = 0;

    if (s == null) {
        err = GretlError.data;
    } else {
        ret = gretl_strdup(s);
        if (ret == null) {
            err = GretlError.alloc;
        }
    }

    *cast(char **) ptr = ret;

    return err;
}

int check_sub_object_return(FnArg * arg, FnParam * fp)
{
    GretlType type = gretl_type_get_plain_type(arg.type);
    void *orig = arg_get_data(arg, type);
    void *curr = user_var_get_value_by_name(fp.name);
    int err = 0;

    if (curr == null) {
        fprintf(stderr, "sub_object return: value has disappeared!\n");
        err = GretlError.data;
    } else if (curr != orig) {
        /* reattachment is needed */
        GretlType uptype = user_var_get_type(arg.uvar);
        void *updata = user_var_get_value(arg.uvar);

        if (uptype == GretlType.array) {
            GretlArray * a = cast(GretlArray *) updata;
            int i, n = gretl_array_get_length(a);
            int found = 0;

            for (i=0; i<n; i++) {
                if (gretl_array_get_data(a, i) == orig) {
                    gretl_array_set_data(a, i, curr);
                    found = 1;
                }
            }
            if (!found) {
                err = GretlError.data;
            }
        } else if (uptype == GretlType.bundle) {
            fprintf(stderr, "sub_object return: uptype = bundle, not handled\n");
            err = GretlError.data;
        } else {
            /* no other types handled */
            err = GretlError.types;
        }
    }

    return err;
}

int is_pointer_arg(Fncall * call, GretlType rtype)
{
    Ufunc * u = call.fun;

    if (call.retname != null) {
        FnParam * fp;
        int i;

        for (i=0; i<call.argc; i++) {
            fp = &u.params[i];
            if ((fp.type == gretl_type_get_ref_type(rtype) ||
                 (fp.flags & FpFlags.const_)) &&
                strcmp(fp.name, call.retname) == 0) {
                return 1;
            }
        }
    }

    return 0;
}

/* On completing execution of a function, move a series at the current
   level up to the level of the caller. This must be done for series
   that were passed to the function by reference rather than being
   copied -- that is, series passed in explicit pointer form and also
   series "loaned" by the caller via the "const" qualifier.
*/

void push_series_to_caller(Ufunc * u, FnArg * arg, DatasetRaw * dset)
{
    int v = arg.val.idnum;

    if (v == 0) {
	/* this is a no-op */
	return;
    } else if (arg.upname == null) {
	/* an internal error! */
        fprintf(stderr, "ERROR in push_series_to_caller: arg->upname is NULL\n");
        return;
    }

    series_decrement_stack_level(dset, v);
    if (series_get_stack_level(dset, v) < 0) {
        fprintf(stderr, "@@@ After decrement in %s, stack level=%d for %s @@@\n",
                u.name, series_get_stack_level(dset, v), arg.upname);
        /* set_stack_level(dset, v, 0); */
    }
    strcpy(dset.varname[v], arg.upname);
}

void push_object_to_caller(FnArg * arg)
{
    user_var_adjust_level(arg.uvar, -1);
    user_var_set_name(arg.uvar, arg.upname);
}

bool null_return(GretlType t)
{
    return t == GretlType.void_ || t == GretlType.none;
}

bool needs_dataset(GretlType t)
{
    return t == GretlType.series || t == GretlType.list || t == GretlType.series_ref || t == GretlType.lists || t == GretlType.lists_ref;
}

int function_assign_returns(Fncall * call, GretlType rtype, DatasetRaw * dset, void * ret, char ** label, SeriesTable ** stab, Prn * prn, int * perr)
{
    Ufunc * u = call.fun;
    int i;
    int err = 0;


    if (*perr == 0 && !null_return(rtype) && call.retname == null) {
        /* missing return value */
        gretl_errmsg_sprintf(_("Function %s did not provide the specified return value\n" ~
		"(expected %s)").ptr, u.name, gretl_type_get_name(rtype));
        *perr = err = GretlError.unkvar;
    } else if (*perr == 0 && needs_dataset(rtype) && dset == null) {
        /* "can't happen" */
        *perr = err = GretlError.data;
    }

    if (*perr == 0) {
        /* first we work on the value directly returned by the
           function (but only if there's no error)
        */
        int copy = is_pointer_arg(call, rtype);

        if (rtype == GretlType.numeric) {
            rtype = fix_numeric_return_type(call, dset);
        }

        if (rtype == GretlType.double_) {
            err = handle_scalar_return(call, ret, rtype);
        } else if (rtype == GretlType.series) {
            err = handle_series_return(call.retname, ret, dset, copy, label, stab);
        } else if (rtype == GretlType.matrix) {
            err = handle_matrix_return(call.retname, ret, copy);
        } else if (rtype == GretlType.list) {
            /* note: in this case the job is finished in
               stop_fncall(); here we just adjust the info on the
               listed variables so they don't get deleted
            */
            err = unlocalize_list(call, call.retname, null, ret, dset);
        } else if (rtype == GretlType.bundle) {
            err = handle_bundle_return(call, ret, copy);
        } else if (gretl_array_type(rtype)) {
            err = handle_array_return(call, ret, copy);
        } else if (rtype == GretlType.string_) {
            err = handle_string_return(call.retname, ret);
        }

        if (err == GretlError.unkvar) {
            pprintf(prn, _("Function %s did not provide the specified return value\n").ptr,
                    u.name);
        }

        *perr = err;
    }

    /* "Indirect return" values and other pointerized args:
       these should be handled even if the function bombed.
    */

    for (i=0; i<call.argc; i++) {
        FnArg * arg = &call.args[i];
        FnParam * fp = &u.params[i];
        int ierr = 0;

        if (needs_dataset(fp.type) && dset == null) {
            ierr = GretlError.data;
        } else if (gretl_ref_type(fp.type)) {
            if (arg.type == GretlType.series_ref) {
                push_series_to_caller(u, arg, dset);
            } else if (arg.upname != null) {
                push_object_to_caller(arg);
            } else if (arg.uvar != null) {
                ierr = check_sub_object_return(arg, fp);
            } else {
                ; /* pure "shell" object: no-op */
            }
        } else if (arg.type == GretlType.useries &&
                   (fp.type == GretlType.series ||
                    fp.type == GretlType.numeric) &&
                   param_is_const(fp)) {
            push_series_to_caller(u, arg, dset);
        } else if ((fp.type == GretlType.matrix ||
                    fp.type == GretlType.bundle ||
                    fp.type == GretlType.string_ ||
                    fp.type == GretlType.numeric ||
                    gretl_array_type(fp.type))) {
            if (arg.shifted) {
                /* non-pointerized const object argument,
                   which we renamed and shifted
                */
                push_object_to_caller(arg);
            }
        } else if (fp.type == GretlType.list) {
            unlocalize_list(call, fp.name, arg, null, dset);
        }
        if (ierr) {
            *perr = err = ierr;
        }
    }


    return err;
}

/* make a record of the sample information at the time a function is
   called */

void record_obs_info(ObsInfo * oi, DatasetRaw * dset)
{
    oi.added = 0;
    oi.structure = dset.structure;
    oi.pd = dset.pd;
    oi.t1 = dset.t1;
    oi.t2 = dset.t2;
    strcpy(oi.stobs, dset.stobs.ptr);
    oi.panel_pd = dset.panel_pd;
    oi.panel_sd0 = dset.panel_sd0;
}

/* On function exit, restore the observations information
   ("setobs" stuff) that was in force on entry.
*/

void restore_obs_info(ObsInfo * oi, DatasetRaw * dset)
{
    dset.structure = oi.structure;
    dset.pd = oi.pd;
    dset.t1 = oi.t1;
    dset.t2 = oi.t2;
    strcpy(dset.stobs.ptr, oi.stobs);
    dset.panel_pd = oi.panel_pd;
    dset.panel_sd0 = oi.panel_sd0;
}

void push_verbosity(Fncall * call)
{
    if (!is_recursing(call)) {
        if (gretl_messages_on()) {
            call.flags |= FCFlags.prev_msgs;
        } else {
            call.flags &= ~FCFlags.prev_msgs;
        }
        if (gretl_echo_on()) {
            call.flags |= FCFlags.prev_echo;
        } else {
            call.flags &= ~FCFlags.prev_echo;
        }
    }
}

void pop_verbosity(Fncall * call)
{
    if (!is_recursing(call)) {
        set_gretl_messages(call.flags & FCFlags.prev_msgs);
        set_gretl_echo(call.flags & FCFlags.prev_echo);
    }
}

/* do the basic housekeeping that is required when a function exits:
   destroy local variables, restore previous sample info, etc.
*/

int stop_fncall(Fncall * call, GretlType rtype, void * ret, DatasetRaw * dset, Prn * prn, int redir_level)
{
    int i, d = gretl_function_depth();
    int delv;
    int anyerr = 0;
    int err = 0;


    /* below: delete series local to the function, taking care not to
       delete any that have been "promoted" to caller level via their
       inclusion in a returned list; also catch any "listarg" series
       that are no longer attached to a list
    */

    if (dset != null) {
        for (i=call.orig_v, delv=0; i<dset.v; i++) {
            if (series_get_stack_level(dset, i) == d) {
                delv++;
            }
        }
        if (delv > 0) {
            if (delv == dset.v - call.orig_v) {
                /* deleting all added series */
                anyerr = dataset_drop_last_variables(dset, delv);
                if (anyerr && !err) {
                    err = anyerr;
                }
            } else {
                for (i=dset.v-1; i>=call.orig_v; i--) {
                    if (series_get_stack_level(dset, i) == d) {
                        anyerr = dataset_drop_variable(i, dset);
                        if (anyerr && !err) {
                            err = anyerr;
                        }
                    }
                }
            }
        }
        if (call.listvars != null) {
            int vi;

            for (i=1; i<=call.listvars[0]; i++) {
                vi = call.listvars[i];
                series_unset_flag(dset, vi, VarFlags.listarg);
                series_set_stack_level(dset, vi, d - 1);
            }
        }
    }

    /* direct list return: write the possibly revised list to the
       return pointer.  Note that we can't do this earlier, because
       the ID numbers of the variables in the return list may be
       changed due to the deletion of function-local variables.
    */

    if (!err && rtype == GretlType.list && ret != null) {
        int *lret = gretl_list_copy(get_list_by_name(call.retname));

        if (lret != null) {
            *cast(int **) ret = lret;
        } else {
            err = GretlError.alloc;
        }
    }

    /* now we're ready to trash function-local vars */
    anyerr = destroy_user_vars_at_level(d);

    if (anyerr && !err) {
        err = anyerr;
    }

    /* if any anonymous equations system was defined: clean up */
    delete_anonymous_equation_system(d);

    if (call.fun.flags & UfunAttrs.uses_set) {
        pop_program_state();
    } else {
        pop_verbosity(call);
    }

    if (dset != null && dset.v > 0) {
        restore_obs_info(&call.obs, dset);
    }

    set_executing_off(&call, dset, prn);

    if (print_redirection_level(prn) > redir_level) {
        gretl_errmsg_set(_("Incorrect use of 'outfile' in function").ptr);
        err = 1;
    }

    if (call != null && call.retname != null) {
        /* we're surely done with this now: avoid potential leak */
        free(call.retname);
        call.retname = null;
    }

    return err;
}

/* reset any saved uservar addresses in context of loops or
   genrs attached to @u */

void reset_saved_uservars(Ufunc * u, int on_recurse)
{
    Stmt * line;
    int i;


    for (i=0; i<u.n_lines; i++) {
        line = &(u.lines[i]);
        if (line_has_loop(line)) {
            loop_reset_uvars(line.ptr);
        }
    }

    /* 2022-09-11: we shouldn't have to do genr_reset_uvars (below)
       any more, with genrs now attached to Fncall rather than line of
       function?
    */
    if (on_recurse) return;

    if (u.call != null && u.call.lgen != null) {
        for (i=0; i<u.call.n_lgen; i++) {
            genr_reset_uvars(u.call.lgen[i].genr);
        }
    }
}

void set_pkgdir(Fnpkg * pkg)
{
    const(char) * p = strrslash(pkg.fname);

    if (p != null) {
        char *pkgdir = gretl_strndup(pkg.fname, p - pkg.fname);

        gretl_insert_builtin_string("pkgdir", pkgdir);
        free(pkgdir);
    }
}

int start_fncall(Fncall * call, DatasetRaw * dset, Prn * prn)
{
    GList *tmp = g_list_last(callstack);
    Fncall * prevcall;

    while (tmp != null) {
        prevcall = cast(Fncall *) tmp.data;
        if (prevcall.fun == call.fun) {
            set_call_recursing(call);
            reset_saved_uservars(call.fun, 1); /* 2022-08-27 */
            break;
        }
        tmp = tmp.prev;
    }

    set_previous_depth(fn_executing);
    fn_executing++;

    if (call.fun.flags & UfunAttrs.uses_set) {
        push_program_state();
    } else {
        push_verbosity(call);
    }

    callstack = g_list_append(callstack, call);


    switch_uservar_hash(fn_executing);

    if (dset != null && dset.v > 0) {
        record_obs_info(&call.obs, dset);
    }

    if (gretl_debugging_on()) {
        set_gretl_echo(1);
        set_gretl_messages(1);
        pputs(prn, "*** ");
        bufspace(gretl_function_depth(), prn);
        pprintf(prn, "executing function %s\n", call.fun.name);
    } else {
        set_gretl_echo(0);
        set_gretl_messages(0);
    }

    if (fn_executing == 1 && call.fun.pkg != null) {
        set_pkgdir(call.fun.pkg);
    }

    return 0;
}

extern (C) int func_exec_callback(ExecState * s, void * ptr, GretlObjType type)
{
    if (s.cmd.ci == GretlCmdIndex.flush || is_plotting_command(s.cmd)) {
        /* we permit "reach-back" into the GUI for these */
        ExecCallback gc = get_gui_callback();

        if (gc != null) {
            gc(s, null, cast(GretlObjType) 0);
        }
    }

    return 0;
}

double arg_get_double_val(FnArg * arg)
{
    if (gretl_scalar_type(arg.type)) {
        return arg.val.x;
    } else if (arg.type == GretlType.matrix) {
        return arg.val.m.val[0];
    } else {
        return NADBL;
    }
}

int check_function_args(Fncall * call, Prn * prn)
{
    Ufunc * u = call.fun;
    FnArg * arg;
    FnParam * fp;
    double x;
    int i;
    int err = 0;

    for (i=0; i<call.argc && !err; i++) {
        int scalar_check = 1;

        arg = &call.args[i];
        fp = &u.params[i];

        if (param_is_optional(fp) && arg.type == GretlType.none) {
            ; /* this is OK */
        } else if (gretl_scalar_type(fp.type) && arg.type == GretlType.double_) {
            ; /* OK: types match */
        } else if (fp.type == GretlType.series && arg.type == GretlType.useries) {
            ; /* OK: types match */
        } else if (gretl_scalar_type(fp.type) &&
                   arg.type == GretlType.matrix &&
                   gretl_matrix_is_scalar(arg.val.m)) {
            ; /* OK: either types match or we can convert */
        } else if (fp.type == GretlType.matrix && arg.type == GretlType.double_) {
            ; /* OK, we can convert */
        } else if (fp.type == GretlType.list && arg.type == GretlType.useries) {
            ; /* OK, we'll handle it (create singleton list) */
        } else if (fp.type == GretlType.list && arg.type == GretlType.none) {
            ; /* OK (special: "null" was passed as argument) */
        } else if (gretl_scalar_type(fp.type) && arg.type == GretlType.none &&
                   !no_scalar_default(fp)) {
            scalar_check = 0; /* OK: arg missing but we have a default value */
        } else if (fp.type == GretlType.numeric && NUMERIC_TYPE(arg.type)) {
            ; /* OK, for overloaded param */
        } else if (fp.type != arg.type) {
            pprintf(prn, _("%s: argument %d is of the wrong type (is %s, should be %s)\n").ptr,
                    u.name, i + 1, gretl_type_get_name(arg.type),
                    gretl_type_get_name(fp.type));
            err = GretlError.types;
        }

        if (!err && scalar_check && fp.type == GretlType.double_) {
            /* check for out-of-bounds scalar value */
            x = arg_get_double_val(arg);
            if ((!na(fp.min) && x < fp.min) ||
                (!na(fp.max) && x > fp.max)) {
                pprintf(prn, _("%s, argument %d: value %g is out of bounds\n").ptr,
                        u.name, i + 1, x);
                err = GretlError.invarg;
            }
        }
    }

    for (i=call.argc; i<u.n_params && !err; i++) {
        /* do we have defaults for any empty args? */
        fp = &u.params[i];
        if (!param_is_optional(fp) && no_scalar_default(fp)) {
            pprintf(prn, _("%s: not enough arguments\n").ptr, u.name);
            err = GretlError.args;
        }
    }


    return err;
}

const(char) * get_funcerr_message(ExecState * state)
{
    const(char) * p = state.cmd.param;

    if (p != null && strlen(p) == gretl_namechar_spn(p)) {
        const(char) * s = get_string_by_name(p);

        if (s != null) {
            return s;
        }
    } else if (p == null) {
        return "none";
    }

    return p;
}

void set_func_error_message(int err, Ufunc * u, ExecState * state, const(char) * cmdline, Stmt * fline)
{
    if (err == GretlError.funcerr) {
        /* let the function writer set the message? */
        const(char) * msg = get_funcerr_message(state);

        if (msg != null && strcmp(msg, "none")) {
            gretl_errmsg_sprintf(_("Error message from %s():\n %s").ptr,
                                 u.name, msg);
            return;
        }
    }

    if (err == GretlError.stop) {
        ; /* no-op */
    } else {
        /* we'll handle this here */
        const(char) * msg = gretl_errmsg_get();
        int showline = 1;

        if (*msg == '\0') {
            msg = errmsg_get_with_default(err);
            gretl_errmsg_set(msg);
        }

        if (*cmdline == '\0' || strncmp(cmdline, "funcerr(", 8) == 0 ||
            strncmp(cmdline, "errorif(", 8) == 0) {
            showline = 0;
        }

        if (fline == null) {
            /* indicates that we don't have a real line number
               available (the error occurred inside a loop)
            */
            if (showline && err != GretlError.fnest) {
                gretl_errmsg_sprintf(_("*** error within loop in function %s\n> %s").ptr,
                                     u.name, cmdline);
            } else {
                gretl_errmsg_sprintf(_("*** error within loop in function %s\n").ptr,
                                     u.name);
            }
        } else {
            if (showline && err != GretlError.fnest) {
		/* 2024-05-12: last arg here was cmdline */
                gretl_errmsg_sprintf(_("*** error in function %s, line %d\n> %s").ptr,
                                     u.name, fline.idx, fline.s);
            } else {
                gretl_errmsg_sprintf(_("*** error in function %s, line %d\n> %s\n").ptr,
                                     u.name, fline.idx, fline.s);
            }
        }

        if (gretl_function_depth() > 1) {
            GList *tmp = g_list_last(callstack);

            if (tmp != null) {
                tmp = g_list_previous(tmp);
                if (tmp != null) {
                    Fncall * call = cast(Fncall *) tmp.data;

                    gretl_errmsg_sprintf(_(" called by function %s").ptr, call.fun.name);
                }
            }
        }
    }
}

Generator * fnline_get_genr(Fncall * call, Stmt * line)
{
    Generator * genr = null;
    int i;

    for (i=0; i<call.n_lgen; i++) {
        if (call.lgen[i].idx == line.idx) {
            genr = call.lgen[i].genr;
            break;
        }
    }


    return genr;
}

int fnline_set_genr(Fncall * call, Stmt * line, Generator * genr)
{
    LineGen * lgen;
    int n = call.n_lgen;
    int err = 0;

    lgen = cast(LineGen *) realloc(call.lgen, (n + 1) * (*lgen).sizeof);
    if (lgen == null) {
        err = GretlError.alloc;
    } else {
        call.lgen = lgen;
        call.lgen[n].idx = line.idx;
        call.lgen[n].genr = genr;
        call.n_lgen = n + 1;
    }

    return err;
}

int generate_return_value2(Fncall * call, ExecState * state, DatasetRaw * dset, const(char) * s, Stmt * line)
{
    gchar *formula = g_strdup_printf("$retval=%s", s);
    int done = 0;
    int err = 0;

    if (line != null) {
        /* try compilation */
        Generator * genr;

        genr = genr_compile(formula, dset, call.fun.rettype,
                            Gretlopt.OPT_P, state.prn, &err);
        if (!err && genr != null) {
            /* succeeded */
            fnline_set_genr(call, line, genr);
            done = 1;
        } else if (err == GretlError.eqn) {
            /* failed, but not fatally */
            fprintf(stderr, "did NOT attach genr at %p\n", line.ptr);
            line.flags |= LnFlags.nocomp;
            err = 0;
        }
    }
    if (!done && !err) {
        err = generate(formula, dset, call.fun.rettype,
                       Gretlopt.OPT_P, state.prn);
    }
    g_free(formula);

    return err;
}

int generate_return_value1(Fncall * call, ExecState * state, DatasetRaw * dset, const(char) * s)
{
    gchar *formula = g_strdup_printf("$retval=%s", s);
    int err;

    err = generate(formula, dset, call.fun.rettype, Gretlopt.OPT_P, state.prn);
    g_free(formula);

    return err;
}

int handle_return_statement(Fncall * call, ExecState * state, DatasetRaw * dset, int gencomp, Stmt * line)
{
    Ufunc * fun = call.fun;
    const(char) * s = null;
    int err = 0;

    if (line != null) {
        Generator *genr = fnline_get_genr(call, line);

        if (genr != null) {
            /* FIXME: when exactly is this reset required? */
            genr_reset_uvars(genr);
            err = execute_genr(genr, dset, state.prn);
            goto last_step;
        } else {
            err = parse_command_line(state, dset, null);
            if (err) {
                goto last_step;
            }
        }
    }

    s = state.cmd.vstart;


    if (fun.rettype == GretlType.void_) {
        if (s == null || *s == '\0') {
            ; /* plain "return" from void function: OK */
        } else if (line == null) {
            gretl_errmsg_sprintf(_("%s: non-null return value '%s' is not valid").ptr,
                                 fun.name, s);
            err = GretlError.types;
        } else {
            gretl_errmsg_sprintf(_("%s, line %d: non-null return value '%s' is not valid").ptr,
                                 fun.name, line.idx, s);
            err = GretlError.types;
        }
        return err; /* handled */
    }

    if (s == null || *s == '\0') {
        if (line == null) {
            gretl_errmsg_sprintf(_("%s: return value is missing").ptr, fun.name);
        } else {
            gretl_errmsg_sprintf(_("%s, line %d: return value is missing").ptr,
                                 fun.name, line.idx);
        }
        err = GretlError.types;
    } else {
        if (gretl_namechar_spn(s) == strlen(s)) {
            /* returning a named variable, simple */
            call.retname = gretl_strdup(s);
        } else {
            if (line != null && gencomp && !cmd_subst(state.cmd)) {
                /* allow trying compilation */
                err = generate_return_value2(call, state, dset, s, line);
            } else {
                /* play it safe */
                err = generate_return_value1(call, state, dset, s);
            }
        }
    }

 last_step:

    if (err) {
        /* FIXME case of line == null? */
        if (s == null && line != null) {
            s = line.s;
        }
        set_func_error_message(err, fun, state, s, line);
    } else if (call.retname == null) {
        call.retname = gretl_strdup("$retval");
    }

    return err;
}

extern (C) int current_function_size()
{
    Ufunc *u = currently_called_function();

    return (u != null)? u.n_lines : 0;
}

__gshared char * return_line;

/* to be called when a "return" statement occurs within a
   loop that's being called by a function
*/

extern (C) int set_function_should_return(const(char) * line)
{
    if (gretl_function_depth() > 0) {
        g_free(return_line);
        return_line = g_strdup(line);
        return 0;
    } else {
        gretl_errmsg_set(_("return: can only be used in a function").ptr);
        return GretlError.parse;
    }
}

int get_return_line(ExecState * state)
{
    if (return_line == null) {
        return 0;
    } else {
        strcpy(state.line, return_line);
        g_free(return_line);
        return_line = null;
        return 1;
    }
}


/* Under "ONE_WORKSPACE" we use a shared command line (workspace) for
   function calls.
*/

ExecState * make_func_exec_state(Cmd * cmd, DatasetRaw * dset, Prn * prn, int * err)
{
    static __gshared char * shared_line;
    char * line = null;
    int free_line = 0;
    ExecState * state;
    Model * model;

    if (shared_line == null) {
	line = shared_line = cast(char *) calloc(MAXLINE, 1);
    } else {
	line = shared_line;
    }

    if (line == null) {
	*err = GretlError.alloc;
	return null;
    }

    state = cast(ExecState *) calloc(1, (*state).sizeof);
    model = allocate_working_model();

    if (state == null || model == null) {
        *err = GretlError.alloc;
    } else {
        *err = gretl_cmd_init(cmd);
    }

    if (!*err) {
        gretl_exec_state_init(state, ExecFlags.function_exec, line, cmd, model, prn);
        if (dset != null) {
            if (dset.submask != null) {
                state.submask = copy_dataset_submask(dset, err);
            }
            state.padded = dset.padmask != null;
        }
        state.callback = &func_exec_callback;
	state.free_line = cast(ubyte) free_line;
    }

    return state;
}

void restore_dataset_for_caller(ExecState * state, DatasetRaw * dset, int orig_t1, int orig_t2, int orig_n)
{
    if (complex_subsampled()) {
        if (state.submask == null) {
            /* we were not sub-sampled on entry */
            restore_full_sample(dset, null);
        } else if (submask_cmp(state.submask, dset.submask)) {
            /* we were sub-sampled differently on entry */
            Gretlopt opt = state.padded ? Gretlopt.OPT_B : Gretlopt.OPT_NONE;

            restore_full_sample(dset, null);
            restrict_sample_from_mask(state.submask, dset, opt);
        }
    } else if (dset.n > orig_n) {
        /* some observations were added inside the function; note
           that this is not allowed if the dataset is subsampled
           on entry
        */
        dataset_drop_observations(dset, dset.n - orig_n);
    }

    dset.t1 = orig_t1;
    dset.t2 = orig_t2;
}

bool do_if_check(int c)
{
    return c == GretlCmdIndex.if_ || c == GretlCmdIndex.elif || c == GretlCmdIndex.else_ || c == GretlCmdIndex.endif;
}
bool may_have_genr(int c)
{
    return c == GretlCmdIndex.if_ || c == GretlCmdIndex.elif || c == GretlCmdIndex.genr;
}

extern (C) int gretl_function_exec_full(Fncall * call, GretlType rtype, DatasetRaw * dset, void * ret, char ** descrip, SeriesTable ** stab, Prn * prn)
{
    Ufunc * u = call.fun;
    ExecState * state = null;
    Generator * genr = null;
    Cmd cmd;

    memset(&cmd, 0, Cmd.sizeof);
    void * ptr = null;
    int orig_n = 0;
    int orig_t1 = 0;
    int orig_t2 = 0;
    int indent0 = 0;
    int started = 0;
    int redir_level = 0;
    int retline = -1;
    int gencomp = 0;
    int n_saved = 0;
    int i;
    int j;
    int err = 0;


    err = maybe_check_function_needs(dset, u);
    if (err) {
        maybe_destroy_fncall(&call);
        return err;
    }

    /* record incoming dataset dimensions */
    if (dset != null) {
        call.orig_v = dset.v;
        orig_n = dset.n;
        orig_t1 = dset.t1;
        orig_t2 = dset.t2;
    } else {
        call.orig_v = 0;
    }

    /* precaution aganst early error */
    indent0 = gretl_if_state_record();

    err = check_function_args(call, prn);
    if (!err) {
        err = allocate_function_args(call, dset);
    }
    if (err) {
        /* get out before allocating further storage */
        maybe_destroy_fncall(&call);
        return err;
    }

    state = make_func_exec_state(&cmd, dset, prn, &err);


    if (!err) {
        start_fncall(call, dset, prn);
        started = 1;
        redir_level = print_redirection_level(prn);
    }

    /* when should we try to compile genrs, loops? */
    /* add clause to cut out functions that recurse */
    gencomp = gretl_iterating() && !get_loop_renaming() &&
	!function_is_recursive(u);

    /* get function lines in sequence and check, parse, execute */

    for (i=0; i<u.n_lines && !err; i++) {
        Stmt * fline = &u.lines[i];
        int this_gencomp = gencomp && !line_no_comp(fline);

        if (gretl_echo_on()) {
            pprintf(prn, "? %s\n", fline.s);
        }
        if (ignore_line(fline)) {
            continue;
        }
        strcpy(state.line, fline.s); /* needed? */
        call.line = fline;

        if (this_gencomp && may_have_genr(fline.ci)) {
            genr = fnline_get_genr(call, fline);
        } else {
            genr = null;
        }


        /* check and adjust the if-state */
        if (do_if_check(fline.ci)) {
            if (fline.ci == GretlCmdIndex.else_ || fline.ci == GretlCmdIndex.endif) {
                state.cmd.ci = fline.ci;
                flow_control(state, null, null);
            } else if (genr != null) {
                state.cmd.ci = fline.ci;
                flow_control(state, dset, &genr);
            } else {
                ptr = this_gencomp ? &genr : null;
                err = maybe_exec_line(state, dset, ptr);
                if (genr != null) {
                    fnline_set_genr(call, fline, genr);
                } else if (ptr != null) {
                    fline.flags |= LnFlags.nocomp;
                }
            }
            if (err) {
                goto err_next;
            }
            if (gretl_if_state_false() && fline.next > 0) {
                /* skip to next relevant statement */
                i = fline.next - 1;
            }
            continue;
        }

        if (fline.ci == GretlCmdIndex.loop) {
            if (fline.ptr != null) {
                state.loop = cast(LoopSet *) fline.ptr;
                err = gretl_loop_exec(state, dset);
                n_saved++;
            } else {
                /* assemble then execute the loop */
                if (this_gencomp && !gretl_function_recursing()) {
                    ptr = &fline.ptr;
                } else {
                    ptr = null;
                }
                for (j=i; j<=fline.next && !err; j++) {
                    strcpy(state.line, u.lines[j].s);
                    err = maybe_exec_line(state, dset, ptr);
                }
                if (!err) {
                    state.loop = cast(LoopSet *) fline.ptr;
		    err = gretl_loop_exec(state, dset);
                }
            }
            if (err) {
                set_func_error_message(err, u, state, state.line, null);
                break;
            } else if (get_return_line(state)) {
                /* a "return" statement was encountered in the loop */
                err = handle_return_statement(call, state, dset,
                                              this_gencomp, null);
                if (i < u.n_lines) {
                    retline = i;
                }
                break;
            }
            /* skip past the loop lines */
            i = fline.next;
            continue;
        } else if (fline.ci == GretlCmdIndex.funcret) {
            err = handle_return_statement(call, state, dset, this_gencomp, fline);
            if (i < u.n_lines) {
                retline = i;
            }
            break;
        } else if (genr != null) {
            err = execute_genr(genr, dset, prn);
            n_saved++;
        } else {
            ptr = (this_gencomp && may_have_genr(fline.ci))? &genr : null;
            err = maybe_exec_line(state, dset, ptr);
            if (genr != null) {
                fnline_set_genr(call, fline, genr);
            } else if (ptr != null) {
                fline.flags |= LnFlags.nocomp;
            }
        }

    err_next:

        if (err) {
            set_func_error_message(err, u, state, state.line, fline);
            break;
        }
    }


    if (dset != null) {
        /* restore the sample that was in place on entry */
        restore_dataset_for_caller(state, dset, orig_t1, orig_t2, orig_n);
    }

    if (err || (retline >= 0)) {
        /* hit an error, or returned prior to the end of the function */
        gretl_if_state_reset(indent0);
    } else {
        err = gretl_if_state_check(indent0);
    }

    function_assign_returns(call, rtype, dset, ret, descrip, stab, prn, &err);

    if (gencomp || n_saved > 0) {
        reset_saved_uservars(call.fun, 0);
    }

    gretl_exec_state_clear(state);
    free(state);

    if (started) {
        int stoperr = stop_fncall(call, rtype, ret, dset, prn, redir_level);

        if (stoperr && !err) {
            err = stoperr;
        }
    } else {
        maybe_destroy_fncall(&call);
    }


    return err;
}

extern (C) int gretl_function_exec(Fncall * call, GretlType rtype, DatasetRaw * dset, void * ret, Prn * prn)
{
    return gretl_function_exec_full(call, rtype, dset, ret,
                                    null, null, prn);
}

/* look up name of supplied argument based on name of variable
   inside function */

extern (C) char * gretl_func_get_arg_name(const(char) * argvar, int * err)
{
    Fncall *call = current_function_call();
    char * ret = null;

    *err = GretlError.data;

    if (call != null) {
        Ufunc * u = call.fun;
        int i;
        int n = call.argc;

        for (i=0; i<n; i++) {
            if (!strcmp(argvar, u.params[i].name)) {
                *err = 0;
                if (call.args[i].upname != null) {
                    ret = gretl_strdup(call.args[i].upname);
                } else {
                    ret = gretl_strdup("");
                }
                if (ret == null) {
                    *err = GretlError.alloc;
                }
                break;
            }
        }
    }

    return ret;
}

/**
 * object_is_const:
 * @name: name of object.
 * @vnum: ID number (specific to objects of type series).
 *
 * Checks whether the named object currently has 'const' status,
 * by virtue of its being made available as a const argument
 * to a user-defined function. Note that @name is the name by
 * which the object is known within the function, which will
 * likely differ from its name in the caller.
 *
 * Returns: non-zero if the object is const, 0 if it is not.
 */

extern (C) int object_is_const(const(char) * name, int vnum)
{
    Fncall *call = current_function_call();
    int ret = 0;

    if (call != null) {
        if (name != null) {
            FnParam * param;
            int i;

            for (i=0; i<call.fun.n_params; i++) {
                param = &call.fun.params[i];
                if (!strcmp(name, param.name)) {
                    ret = param_is_const(param);
                    break;
                }
            }
        }
        if (!ret && vnum > 0 && vnum < call.orig_v) {
            /* We're looking at a series that is not local to
               the called function, but not yet identified as
               read-only. It probably _should_ be read-only
               unless it was given in pointer form.
               Note: this check added 2018-10-18.0
            */
            if (!in_gretl_list(call.ptrvars, vnum)) {
                ret = 1;
            }
        }
    }

    return ret;
}

/**
 * object_is_function_arg:
 * @name: name of object (e.g. matrix).
 *
 * Checks whether the named object has been made available
 * as a function argument.
 *
 * Returns: non-zero (the gretl type of the object) if the
 * object has arg status, 0 otherwise.
 */

extern (C) int object_is_function_arg(const(char) * name)
{
    Fncall *call = current_function_call();

    if (call != null) {
        FnParam * params = call.fun.params;
        int i;

        for (i=0; i<call.fun.n_params; i++) {
            if (!strcmp(name, params[i].name)) {
                return params[i].type;
            }
        }
    }

    return 0;
}

__gshared int allow_full_data;

extern (C) void allow_full_data_access(int s)
{
    if (s > 0) {
        allow_full_data = 1;
    } else if (s == 0) {
        allow_full_data = 0;
    }
}

/**
 * sample_range_get_extrema:
 * @dset: dataset info.
 * @t1: location to receive earliest possible starting
 * observation, or null.
 * @t2: location to receive latest possible ending observation,
 * or null.
 *
 * Fills out @t1 and/or @t2, making allowance for the possibility
 * that we're currently executing a function, on entry to
 * which the sample range was restricted: within the function,
 * we are not allowed to overstep the bounds set on entry.
 */

extern (C) void sample_range_get_extrema(const(DatasetRaw) * dset, int * t1, int * t2)
{
    int[2] tvals = [0, dset.n - 1];

    if (!allow_full_data) {
        Fncall *call = current_function_call();

        if (call != null) {
            tvals[0] = call.obs.t1;
            tvals[1] = call.obs.t2 + call.obs.added;
            /* FIXME this needs to be smarter */
            if (tvals[1] > dset.n - 1) {
                tvals[1] = dset.n - 1;
            }
        }
    }

    if (t1 != null) {
        *t1 = tvals[0];
    }
    if (t2 != null) {
        *t2 = tvals[1];
    }
}

extern (C) void extend_function_sample_range(int addobs)
{
    Fncall *call = current_function_call();

    if (call != null) {
        call.obs.added += addobs;
    }
}

/**
 * series_is_accessible_in_function:
 * @ID: series ID number (0 = constant).
 *
 * Returns: 1 if the series with ID number @ID is accessible
 * by number at the current level of function execution,
 * otherwise 0.0
 */

extern (C) int series_is_accessible_in_function(int ID, const(DatasetRaw) * dset)
{
    /* This is harder to get right than one might think, and
       for now (as of 2018-10-18) we will not attempt to
       screen access to series in this way. However, we do
       strive to ensure that series are treated as read-only
       when they ought to be. See geneval.c for (potential)
       uses of this check, currently def'd out.
    */
    return 1; /* allow all */
}

/**
 * gretl_functions_cleanup:
 *
 * For internal use: frees all resources associated with
 * user-defined functions and function packages.
 */

extern (C) void gretl_functions_cleanup()
{
    int i;

    for (i=0; i<n_ufuns; i++) {
        ufunc_free(ufuns[i]);
    }

    free(ufuns);
    ufuns = null;
    n_ufuns = 0;

    for (i=0; i<n_pkgs; i++) {
        function_package_free(pkgs[i]);
    }

    free(pkgs);
    pkgs = null;
    n_pkgs = 0;
}

void show_gfn_help(Fnpkg * pkg, int gui_help, int markup, Prn * prn)
{
    const(char) * text;
    int mdconv;

    text = gui_help ? pkg.gui_help : pkg.help;
    mdconv = help_text_is_markdown(pkg, gui_help);

    if (mdconv) {
        /* help is in markdown */
        md_to_gretl(text, prn);
    } else {
        if (markup) {
            pputs(prn, "<@itl=\"Help text\">:\n\n");
            pputs(prn, "<mono>\n");
        } else {
            pputs(prn, "Help text:\n");
        }
        pputs(prn, text);
        if (markup) {
            pputs(prn, "\n</mono>");
        }
        pputs(prn, "\n\n");
    }
}

/* Generate help output for a packaged function, either for
   display on the console, or with markup for display in a
   GtkTextView window in the GUI (opt & Gretlopt.OPT_M).
*/

void real_user_function_help(Ufunc * fun, Gretlopt opt, Prn * prn)
{
    Fnpkg * pkg = fun.pkg;
    int markup = (opt & Gretlopt.OPT_M);
    int i;

    if (markup) {
        pprintf(prn, "<@itl=\"Function\">: %s\n", fun.name);
    } else {
        pprintf(prn, "Function: %s\n", fun.name);
    }

    if (pkg != null) {
        if (markup) {
            pprintf(prn, "<@itl=\"Package\">: %s %s (%s)\n", pkg.name, pkg.version_,
                    pkg.date);
            pprintf(prn, "<@itl=\"Author\">: %s\n", pkg.author? pkg.author : "unknown");
            if (pkg.email != null && *pkg.email != '\0') {
                pprintf(prn, "<@itl=\"Email\">: %s\n", pkg.email);
            }
        } else {
            pprintf(prn, "Package: %s %s (%s)\n", pkg.name, pkg.version_,
                    pkg.date);
            pprintf(prn, "Author: %s\n", pkg.author? pkg.author : "unknown");
            if (pkg.email != null && *pkg.email != '\0') {
                pprintf(prn, "Email:  %s\n", pkg.email);
            }
        }
        pputc(prn, '\n');
    }

    if ((opt & Gretlopt.OPT_G) && pkg != null && pkg.gui_help != null) {
        /* GUI-specific help is preferred, and is available */
        show_gfn_help(pkg, 1, markup, prn);
        return;
    }

    if (markup) {
        pputs(prn, "<@itl=\"Parameters\">: ");
    } else {
        pputs(prn, "Parameters: ");
    }

    if (fun.n_params > 0) {
        pputc(prn, '\n');
        for (i=0; i<fun.n_params; i++) {
            pprintf(prn, " %s (%s",
                    fun.params[i].name, gretl_type_get_name(fun.params[i].type));
            if (fun.params[i].descrip != null) {
                pprintf(prn, ": %s)\n", fun.params[i].descrip);
            } else {
                pputs(prn, ")\n");
            }
        }
        pputc(prn, '\n');
    } else {
        pputs(prn, "none\n\n");
    }

    if (markup) {
        pputs(prn, "<@itl=\"Return value\">: ");
    } else {
        pputs(prn, "Return value: ");
    }

    if (fun.rettype != GretlType.none && fun.rettype != GretlType.void_) {
        pprintf(prn, "%s\n\n", gretl_type_get_name(fun.rettype));
    } else {
        pputs(prn, "none\n\n");
    }

    if (fun.docstr != null) {
	pprintf(prn, "%s\n", fun.docstr);
    }

    if (pkg != null && pkg.help != null) {
        if (is_pdf_ref(pkg.help)) {
            gchar *pdfname = g_strdup(pkg.fname);
            gchar *p = strrchr(pdfname, '.');

            p[0] = '\0';
            strcat(p, ".pdf");
            if (markup) {
                pprintf(prn, "<@itl=\"Documentation\">: <@adb=\"%s\">\n\n", pdfname);
            } else {
                pprintf(prn, "See %s\n\n", pdfname);
            }
            g_free(pdfname);
        } else {
            show_gfn_help(pkg, 0, markup, prn);
        }
    }

    if (pkg != null && pkg.sample != null) {
        if (markup) {
            pputs(prn, "<@itl=\"Sample script\">:\n\n");
            pputs(prn, "<code>\n");
        } else {
            pputs(prn, "Sample script:\n\n");
        }
        pputs(prn, pkg.sample);
        if (markup) {
            pputs(prn, "\n</code>\n");
        } else {
            pprintf(prn, "\n\n");
        }
    }
}

/**
 * user_function_help:
 * @fnname: name of function.
 * @opt: may include Gretlopt.OPT_M for adding markup, Gretlopt.OPT_G
 * for preferring GUI-specific help, if available.
 * @prn: printing struct.
 *
 * Looks for a function named @fnname and prints
 * as much help information as possible.
 *
 * Returns: 0 on success, non-zero if the function is not
 * found.
 */

extern (C) int user_function_help(const(char) * fnname, Gretlopt opt, Prn * prn)
{
    Ufunc *fun = get_user_function_by_name(fnname);
    int err = 0;

    if (fun == null) {
        pprintf(prn, _("\"%s\" is not defined.\n").ptr, fnname);
        err = 1;
    } else {
        real_user_function_help(fun, opt, prn);
    }

    return err;
}

/**
 * function_package_has_PDF_doc:
 * @pkg: function-package pointer.
 * @pdfname: location to receive basename of PDF file,
 * if applicable (or null if this is not wanted).
 *
 * Checks whether @pkg is documented in the form of a PDF file.
 *
 * Returns: 1 if so (in which case the name of that file is
 * returned via @pdfname), 0 otherwise.
 */

extern (C) int function_package_has_PDF_doc(Fnpkg * pkg, char ** pdfname)
{
    int ret = 0;

    if (pkg.help != null && !strncmp(pkg.help, "pdfdoc", 6)) {
        ret = 1;
        if (pdfname != null) {
            *pdfname = switch_ext_new(pkg.fname, "pdf");
            if (*pdfname == null) {
                ret = 0;
            }
        }
    }

    return ret;
}

/**
 * function_package_has_gui_help:
 * @pkg: function-package pointer.
 *
 * Checks whether @pkg has GUI-specific help text.
 *
 * Returns: 1 if so, 0 otherwise.
 */

extern (C) int function_package_has_gui_help(Fnpkg * pkg)
{
    return pkg.gui_help != null && !string_is_blank(pkg.gui_help);
}

extern (C) void function_package_set_editor(Fnpkg * pkg, void * editor)
{
    if (pkg != null) {
        pkg.editor = editor;
    }
}

extern (C) void * function_package_get_editor(Fnpkg * pkg)
{
    return pkg == null ? null : pkg.editor;
}

/**
 * delete_function_package:
 * @gfnname: full path to gfn file.
 *
 * Deletes @gfnname, taking care of deleting its enclosing
 * specific subdirectory and all of its contents if applicable.
 *
 * Returns: 0 on success, non-zero code on error.
 */

extern (C) int delete_function_package(const(char) * gfnpath)
{
    char *p = strrslash(cast(char *) gfnpath);
    char * pkgname = null;
    char * pkgdir = null;
    char * pkgsub = null;
    int err = 0;

    errno = 0;

    if (p != null) {
        /* @pkgname: strip extension from @gfnpath */
        pkgname = g_strdup(p + 1);
        p = strrchr(pkgname, '.');
        if (p != null) {
            p[0] = '\0';
        }
        /* @pkgdir: directory in which the gfn lives */
        pkgdir = g_strdup(gfnpath);
        p = strrslash(pkgdir);
        p[0] = '\0';
        p = strrslash(pkgdir);
        /* @pkgsub: subdir in which the gfn lives */
        if (p != null) {
            pkgsub = g_strdup(p + 1);
        }
    }

    if (pkgname != null && pkgdir != null && pkgsub != null &&
        !strcmp(pkgname, pkgsub)) {
        /* We should delete the tree @pkgdir only if the last
           directory in that path, @pkgsub, compares equal to
           the basename of the package, @pkgname. For example:
           foo(.gfn) lives in <path-to-functions>/foo/.
        */
        err = gretl_deltree(pkgdir);
        if (err) {
            gretl_errmsg_sprintf(_("Couldn't delete %s").ptr, pkgdir);
        }
    } else {
        /* just delete the .gfn file itself */
        err = gretl_remove(gfnpath);
        if (err) {
            gretl_errmsg_sprintf(_("Couldn't delete %s").ptr, gfnpath);
        }
    }

    if (err) {
        fprintf(stderr, "failure in delete_function_package: gfnpath '%s'\n",
                gfnpath);
        fprintf(stderr, " pkgname '%s', pkgdir '%s', pkgsub '%s'\n",
                pkgname, pkgdir, pkgsub);
    }

    g_free(pkgname);
    g_free(pkgdir);
    g_free(pkgsub);

    return err;
}

/**
 * uninstall_function_package:
 * @package: name of package (with or without .gfn extension).
 * @opt: may include Gretlopt.OPT_P to "purge" the package.
 * @prn: gretl printer.
 *
 * Returns: 0 on success, non-zero code on error.
 */

extern (C) int uninstall_function_package(const(char) * package_, Gretlopt opt, Prn * prn)
{
    char * gfnname = null;
    char * pkgname = null;
    char[FILENAME_MAX] fname;
    char * p;
    int err;

    if (strchr(package_, '/') != null || strchr(package_, '\\') != null ||
        strstr(package_, "..") != null) {
        /* reject anything that is not a plain package name: this
           guards against path-traversal via a script-supplied
           "pkg <name> remove/unload" command
        */
        gretl_errmsg_sprintf(_("Invalid package name '%s'").ptr, package_);
        return GretlError.data;
    }

    if (has_suffix(package_, ".gfn")) {
        gfnname = g_strdup(package_);
        pkgname = g_strdup(package_);
        p = strrchr(pkgname, '.');
        p[0] = '\0';
    } else {
        gfnname = g_strdup_printf("%s.gfn", package_);
        pkgname = g_strdup(package_);
    }

    fname[0] = '\0';
    err = get_full_filename(gfnname, fname.ptr, Gretlopt.OPT_I);

    if (!err && !gretl_file_exists(fname.ptr)) {
        gretl_errmsg_sprintf(_("Couldn't find %s").ptr, gfnname);
        err = GretlError.fopen;
    }

    if (!err) {
        function_package_unload_full_by_filename(fname.ptr);
        if (opt & Gretlopt.OPT_P) {
            err = delete_function_package(fname.ptr);
        }
    }

    if (!err && gretl_messages_on()) {
        if (opt & Gretlopt.OPT_P) {
            pprintf(prn, _("Removed %s\n").ptr, pkgname);
        } else {
            pprintf(prn, _("Unloaded %s\n").ptr, pkgname);
        }
    }

    g_free(gfnname);
    g_free(pkgname);

    return err;
}

/* full apparatus for hansl indentation follows */

void get_cmdword(const(char) * s, char * word)
{
    if (!strncmp(s, "catch ", 6)) {
        s += 6;
    }

    if (sscanf(s, "%*s <- %8s", word) != 1) {
        sscanf(s, "%8s", word);
    }
}

bool bare_quote(const(char) * p, const(char) * s)
{
    return *p == '"' && (p - s == 0 || *(p - 1) != '\\');
}
bool starts_ccmt(const(char) * p)
{
    return *p == '/' && *(p + 1) == '*';
}
bool ends_ccmt(const(char) * p)
{
    return *p == '*' && *(p + 1) == '/';
}

void check_for_comment(const(char) * s, int * incomm)
{
    const(char) * p = s;
    int commbak = *incomm;
    int started = 0;
    int quoted = 0;

    p = s + strspn(s, " \t");
    if (*incomm == 2 && ends_ccmt(p)) {
	*incomm = 0;
	p += 2;
	p += strspn(p, " ");
	return;
    }

    p = s;
    while (*p) {
        if (!quoted && !*incomm && *p == '#') {
            break;
        }
        if (!*incomm && bare_quote(p, s)) {
            quoted = !quoted;
        }
        if (!quoted) {
            if (starts_ccmt(p)) {
		started = 1;
                *incomm = 1;
                p += 2;
            } else if (ends_ccmt(p)) {
                *incomm = started ? 0 : 3;
                p += 2;
                p += strspn(p, " ");
            }
        }
        if (*p) {
            p++;
        }
    }

    if (*incomm > 0 && *incomm < 3 && commbak) {
        /* on the second or subsequent line of a C-style comment */
        *incomm = 2;
    }
}

/* determine whether a given line is subject to continuation,
   i.e. ends with backslash or comma, but not in a comment
*/

int line_broken(const(char) * s)
{
    int ret = 0;

    if (*s != '\0') {
        int i, n = cast(int) strlen(s);

        for (i=n-1; i>=0; i--) {
            if (s[i] == '\\' || s[i] == ',') {
                ret = 1;
            } else if (!ret && !isspace(s[i])) {
                break;
            } else if (ret && s[i] == '#') {
                ret = 0;
                break;
            }
        }
    }

    return ret;
}

void strip_trailing_whitespace(char * s)
{
    int i, n = cast(int) strlen(s);

    for (i=n-1; i>=0; i--) {
        if (s[i] == '\n' || s[i] == ' ' || s[i] == '\t') {
            s[i] = '\0';
        } else {
            break;
        }
    }

    strcat(s, "\n");
}

/* determine position of unmatched left parenthesis,
   when applicable, if @s starts a function definition
*/

int left_paren_offset(const(char) * s)
{
    const(char) * p = strchr(s, '(');

    if (p != null && strchr(p, ')') == null) {
        return cast(int) (p - s);
    } else {
        return 0;
    }
}

int wordmatch(const(char) * s, const(char) * test)
{
    int n = cast(int) strlen(test);

    return (!strncmp(s, test, n) && (s[n] == '\0' || isspace(s[n])));
}

extern (C) void adjust_indent(const(char) * s, int * this_indent, int * next_indent)
{
    const(char) *[15] block_starts = [
        "loop".ptr, "if".ptr, "nls".ptr, "mle".ptr, "gmm".ptr, "mpi".ptr, "plot".ptr,
        "function".ptr, "restrict".ptr, "system".ptr, "foreign".ptr,
        "outfile".ptr, "gpbuild".ptr, "gibbs".ptr, null
    ];
    int ti = *next_indent;
    int ni = *next_indent;
    int i;
    int matched = 0;

    if (*s == '\0') {
        *this_indent = *next_indent;
        return;
    }

    /* skip "catch" if present */
    if (!strncmp(s, "catch ", 6)) {
        s += 6;
        s += strspn(s, " ");
    }

    for (i=0; block_starts[i] != null && !matched; i++) {
        if (wordmatch(s, block_starts[i])) {
            matched = 1;
            ni++;
        }
    }

    if (!matched) {
        if (wordmatch(s, "end") ||
            wordmatch(s, "endif") ||
            wordmatch(s, "endloop")) {
            ti--;
            ni--;
        } else if (wordmatch(s, "else") ||
                   wordmatch(s, "elif")) {
            ni = ti;
            ti--;
        }
    }

    *this_indent = ti;
    *next_indent = ni;
}

extern (C) void normalize_hansl(const(char) * buf, int tabwidth, Prn * prn)
{
    char * line;
    char[9] word;
    const(char) * ins;
    size_t llen = 1024;
    int this_indent = 0;
    int next_indent = 0;
    int continuation = 0;
    int incomment = 0;
    int inforeign = 0;
    int lp_pos = 0;
    int lp_zero = 0;
    int nsp;

    line = cast(char *) malloc(llen);
    line[0] = '\0';
    bufgets_init(buf);

    while (safe_bufgets(&line, &llen, buf)) {
        int handled = 0;
	int lbreak = 0;

        strip_trailing_whitespace(line);

        if (string_is_blank(line)) {
            pputc(prn, '\n');
            continue;
        }

        check_for_comment(line, &incomment);
        ins = line + strspn(line, " \t");

        if (!incomment) {
	    lbreak = line_broken(line);
            word[0] = '\0';
            get_cmdword(ins, word.ptr);
            if (!strcmp(word.ptr, "foreign")) {
                inforeign = 1;
            } else if (inforeign) {
                if (!strncmp(ins, "end foreign", 11)) {
                    inforeign = 0;
                } else {
                    /* insert foreign line as is */
                    pputs(prn, line);
                    handled = 1;
                }
            } else {
                if (!strcmp(word.ptr, "function")) {
		    /* record position of left parenthesis */
                    lp_pos = left_paren_offset(ins);
                } else if (lp_pos > 0 && strchr(ins, ')') != null) {
		    /* flag end of function signature */
                    lp_zero = 1;
                }
		adjust_indent(word.ptr, &this_indent, &next_indent);
            }
        }

        if (!handled) {
            nsp = this_indent * tabwidth;
            if (incomment > 1) {
		/* C-style comment continuation */
		nsp += 3;
            } else if (continuation) {
                if (lp_pos > 0) {
		    /* in a function signature */
                    nsp = lp_pos + 1;
                } else {
		    /* small offset for continuation */
                    nsp += 2;
                }
            }
            /* insert required number of spaces */
            bufspace(nsp, prn);
            /* insert line, with leading space pruned */
            pputs(prn, ins);
        }

	/* set parameters for next line */
	this_indent = next_indent;
	continuation = lbreak;
        if (lp_zero) {
	    /* scrub the left-paren record */
            lp_zero = lp_pos = 0;
        }
	if (incomment == 3) {
	    incomment = 0;
	}
    }

    bufgets_finalize(buf);
    free(line);
}
