/* Port of gretl_func.c, part 2 of 4: see gretl_func0.d. */

module gretl.gretl_func2;

import gretl.gretl_prn : pprintf2;

import gretl.xmlh : xmlParseDTD, xmlValidityErrorFunc, xmlValidityWarningFunc;


import gretl.gretl_func_ext : is_gretl_addon;
import gretl.xmlh : xmlDtdPtr, xmlFreeDoc, xmlFreeDtd, xmlFreeValidCtxt, xmlNewValidCtxt, xmlStrcmp, xmlValidCtxtPtr, xmlValidateDtd;

import gretl.gretl_commands : gretl_command_word;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_func1 : fn_param_default;
import gretl.gretl_typemap : gretl_type_get_name;
import gretl.gretl_utils : gretl_version_string;
import gretl.gretl_xml : gretl_xml_get_gfn;
import gretl.versionh : GRETL_VERSION;

import gretl.cmd_privateh : Cmd;

import gretl.cmd_privateh : Cmd;

import gretl.cmd_privateh : Cmd;

import gretl.gerror : GError, g_error_free;
import gretl.gfileutils : g_file_get_contents, g_get_current_dir;
import gretl.glib : gsize;
import gretl.gmem : g_free;
import gretl.gretl_commands : gretl_command_number;
import gretl.gretl_errors : gretl_errmsg_set, gretl_errmsg_sprintf, gretl_warnmsg_sprintf;
import gretl.gretl_func1 : add_allocated_ufunc, fn_param_optional, function_package_alloc, get_user_function_by_name, n_free_functions, pkg_key_get_role, real_function_package_unload, set_function_private, ufunc_free, write_function_xml;
import gretl.gretl_list : gretl_list_new;
import gretl.gretl_paths : gretl_chdir, gretl_fopen, gretl_home, gretl_remove, gretl_stat;
import gretl.gretl_prn : gretl_print_destroy, gretl_print_new_with_filename, pprintf, pputc, pputs;
import gretl.gretl_utils : gretl_version_number;
import gretl.gretl_xml : gretl_xml_get_strings_array, gretl_xml_header, gretl_xml_node_get_trimmed_string, gretl_xml_open_doc_root, gretl_xml_put_string, gretl_xml_put_strings_array, gretl_xml_put_tagged_string, gretl_xml_put_tagged_string_plus;
import gretl.gretl_zip : package_make_zipfile;
import gretl.gstrfuncs : g_strdup, g_strdup_printf;
import gretl.missingh : NADBL, na;
import gretl.options : get_optval_string;
import gretl.strutils : dot_atof, ends_with_backslash, gretl_charsub, gretl_strdup, gretl_string_split, gretl_strndup, has_suffix, path_last_element, string_is_blank, strings_array_add, strings_array_dup, strings_array_free, strings_array_position, strings_array_prepend_uniq, strings_array_sort, strrslash, switch_ext, tailstrip;
import gretl.strutilsh : SLASH;
import gretl.xmlh : xmlDocPtr, xmlNodePtr;

import gretl.garray : GArray;

import core.stdc.math : fabs, floor, isnan, log, sqrt, pow, exp;
import core.stdc.stdlib : free, malloc, calloc, realloc, qsort, atoi, atof;
import core.stdc.stdio : fprintf, printf, sprintf, snprintf, stderr, stdout, FILE, fopen, fclose, fgets, fputs;
import core.stdc.string : memcpy, memset, strcmp, strlen, strcpy, strncpy, strcat, strchr, strstr, strncmp, strrchr;

import core.stdc.stdarg : va_arg, va_end, va_list, va_start;
import core.stdc.ctype : isalpha, isspace, isdigit;
import core.stdc.errno : errno;
import core.stdc.stdio : fputc, rewind, sscanf, remove;
import core.stdc.stdlib : strtod, strtol;
import core.stdc.string : memmove, strcspn, strncat, strspn;
version (Posix) import core.sys.posix.sys.stat : stat_t; else import gretl.winport : stat_t;
import gretl.genlex : function_lookup, install_function_override;
import gretl.gretl_func0;
import gretl.xmlh : xmlChar;
import gretl.gretl_funch;
import gretl.libgretlh;
import gretl.gretl_matrixh;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_prn : Prn;

private alias XUC = const(xmlChar) *;

/* construct a name for @pkg based on its filename member:
   take the basename and knock off ".gfn"
*/

void name_package_from_filename(Fnpkg * pkg)
{
    char *p = strrslash(pkg.fname);
    int n;

    if (p != null) {
        p++;
    } else {
        p = pkg.fname;
    }

    n = cast(int) strlen(p);
    if (has_suffix(p, ".gfn")) {
        n -= 4;
    }

    if (n > FN_NAMELEN - 1) {
        n = FN_NAMELEN - 1;
    }

    *pkg.name = '\0';
    strncat(pkg.name, p, n);

}

/* name lookup for functions to be connected to a package,
   allowing for the possibility that they're already
   connected */

Ufunc * get_uf_array_member(const(char) * name, Fnpkg * pkg)
{
    int i;

    for (i=0; i<n_ufuns; i++) {
        if (ufuns[i].pkg == pkg || ufuns[i].pkg == null) {
            if (!strcmp(name, ufuns[i].name)) {
                return ufuns[i];
            }
        }
    }

    return null;
}

void check_special_comments(Ufunc * fun)
{
    int i;

    for (i=0; i<fun.n_lines; i++) {
        if (strstr(fun.lines[i].s, "## no-print ##")) {
            fun.flags |= UfunAttrs.noprint;
        } else if (strstr(fun.lines[i].s, "## menu-only ##")) {
            fun.flags |= UfunAttrs.menu_only;
        } else if (strstr(fun.lines[i].s, "## must-assign ##")) {
            fun.flags |= UfunAttrs.assign;
        }
    }
}

/* before revising the function members of an existing
   package, detach all its current functions
*/

void package_disconnect_funcs(Fnpkg * pkg)
{
    int i;

    if (pkg.pub != null) {
        for (i=0; i<pkg.n_pub; i++) {
            pkg.pub[i].pkg = null;
        }
        free(pkg.pub);
        pkg.pub = null;
        pkg.n_pub = 0;
    }

    if (pkg.priv != null) {
        for (i=0; i<pkg.n_priv; i++) {
            pkg.priv[i].pkg = null;
            set_function_private(pkg.priv[i], 0);
        }
        free(pkg.priv);
        pkg.priv = null;
        pkg.n_priv = 0;
    }
}

/* Given an array of @n function names in @names, set up a
   corresponding array of pointers to the named functions in @pkg,
   Flag an error if any of the function names are bad, or if
   allocation fails.
*/

int set_uf_array_from_names(Fnpkg * pkg, char ** names, int n, int priv)
{
    Ufunc **uf = null;
    Ufunc * fun;
    int i;
    int j;

    /* check the supplied names */
    for (i=0; i<n; i++) {
        for (j=0; j<i; j++) {
            if (!strcmp(names[j], names[i])) {
                gretl_errmsg_sprintf(_("Duplicated function name '%s'").ptr, names[i]);
                return GretlError.data;
            }
        }
        if (get_uf_array_member(names[i], null) == null) {
            fprintf(stderr, "%s: function not found!\n", names[i]);
            return GretlError.data;
        }
    }

    /* allocate storage */
    if (n > 0) {
        uf = cast(Ufunc* *) malloc(n * (*uf).sizeof);
        if (uf == null) {
            return GretlError.alloc;
        }
    }

    /* connect the specified functions */
    for (i=0; i<n; i++) {
        fun = get_uf_array_member(names[i], null);
        fun.pkg = pkg;
        set_function_private(fun, priv);
        check_special_comments(fun);
        uf[i] = fun;
    }

    /* set up the package info */
    if (priv) {
        pkg.priv = uf;
        pkg.n_priv = n;
    } else {
        pkg.pub = uf;
        pkg.n_pub = n;
    }

    return 0;
}

/**
 * function_package_connect_funcs:
 * @pkg: function package.
 * @pubnames: array of names of public functions.
 * @n_pub: number of strings in @pubnames.
 * @privnames: array of names of private functions (or %null).
 * @n_priv: number of strings in @privnames (may be 0).
 *
 * Looks up the functions named in @pubnames and @privnames
 * and adds pointers to these functions to @pkg, hence marking
 * the functions as belonging to @pkg.
 *
 * Returns: 0 on success, non-zero on error.
 */

extern (C) int function_package_connect_funcs(Fnpkg * pkg, char ** pubnames, int n_pub, char ** privnames, int n_priv)
{
    int err;

    /* clear the decks first */
    package_disconnect_funcs(pkg);

    err = set_uf_array_from_names(pkg, pubnames, n_pub, 0);

    if (!err) {
        err = set_uf_array_from_names(pkg, privnames, n_priv, 1);
    }

    return err;
}

/**
 * function_package_new:
 * @fname: filename for package.
 * @pubnames: array of names of public functions.
 * @n_pub: number of strings in @pubnames.
 * @privnames: array of names of private functions (or %null).
 * @n_priv: number of strings in @privnames (may be 0).
 * @err: location to receive error code.
 *
 * Allocates a new package with filename-member @fname, including
 * the public and private functions named in @pubnames and
 * @privnames.  Note that this function does not cause the
 * package to be written to file; for that, see
 * write_function_package().
 *
 * Returns: pointer to package on success, %null on error.
 */

extern (C) Fnpkg * function_package_new(const(char) * fname, char ** pubnames, int n_pub, char ** privnames, int n_priv, int * err)
{
    Fnpkg * pkg = null;

    if (n_pub <= 0) {
        /* we need at least one public function */
        *err = GretlError.data;
    } else {
        pkg = function_package_alloc(fname);
        if (pkg == null) {
            *err = GretlError.alloc;
        }
    }

    if (*err) {
        return null;
    }

    name_package_from_filename(pkg);

    *err = function_package_connect_funcs(pkg, pubnames, n_pub,
                                          privnames, n_priv);

    if (!*err) {
        *err = function_package_record(pkg);
    }

    if (*err) {
        /* note: this does not free the packaged functions, if any */
        function_package_free(pkg);
        pkg = null;
    }

    return pkg;
}

char * trim_text(char * s)
{
    while (isspace(*s)) s++;
    return tailstrip(s);
}

int package_write_translatable_strings(Fnpkg * pkg, Prn * prn)
{
    FILE * fp;
    char * trname;
    char **S = null;
    int trans = 0;
    int i;
    int n = 0;

    /* open a file to contain the results */
    trname = g_strdup_printf("%s-i18n.c", pkg.name);
    fp = gretl_fopen(trname, "wb");
    if (fp == null) {
        gretl_errmsg_sprintf(_("Couldn't open %s").ptr, trname);
        g_free(trname);
        return GretlError.fopen;
    }

    if (pkg.pub != null) {
        /* pick up any parameter descriptions and/or parameter-value
           enumeration strings
        */
        int j;
        int k;

        for (i=0; i<pkg.n_pub; i++) {
            Ufunc * fun = pkg.pub[i];

            for (j=0; j<fun.n_params; j++) {
                FnParam * param = &fun.params[j];

                if (param.descrip != null) {
                    strings_array_add(&S, &n, param.descrip);
                }
                for (k=0; k<param.nlabels; k++) {
                    strings_array_add(&S, &n, param.labels[k]);
                }
            }
        }
    }

    if (pkg.label != null || S != null) {
        /* we got some relevant content */
        fprintf(fp, "const char *%s_translations[] = {\n".ptr, pkg.name);
        if (pkg.label != null) {
            fprintf(fp, "    N_(\"%s\"),\n", pkg.label);
        }
        if (S != null) {
            strings_array_sort(&S, &n, Gretlopt.OPT_U);
            for (i=0; i<n; i++) {
                fprintf(fp, "    N_(\"%s\")", S[i]);
                if (i < n-1) {
                    fputc(',', fp);
                }
                fputc('\n', fp);
            }
            strings_array_free(S, n);
        }
        fputs("};\n", fp);
        trans = 1;
    }

    fclose(fp);

    if (!trans) {
        gretl_remove(trname);
    } else {
        pprintf(prn, "Wrote translations file %s\n", trname);
    }

    g_free(trname);

    return 0;
}

int package_write_index(Fnpkg * pkg, Prn * inprn)
{
    Prn * prn;
    char * idxname;
    int err = 0;

    idxname = g_strdup_printf("%s.xml", pkg.name);

    prn = gretl_print_new_with_filename(idxname, &err);
    if (prn == null) {
        g_free(idxname);
        return err;
    }

    gretl_xml_header(prn);
    pprintf(prn, "<gretl-addon name=\"%s\"", pkg.name);

    if (pkg.dreq == DataReq.needs_ts) {
        pprintf(prn, " %s=\"true\"".ptr, NEEDS_TS.ptr);
    } else if (pkg.dreq == DataReq.needs_qm) {
        pprintf(prn, " %s=\"true\"".ptr, NEEDS_QM.ptr);
    } else if (pkg.dreq == DataReq.needs_panel) {
        pprintf(prn, " %s=\"true\"".ptr, NEEDS_PANEL.ptr);
    } else if (pkg.dreq == DataReq.nodata_ok) {
        pprintf(prn, " %s=\"true\"".ptr, NO_DATA_OK.ptr);
    }
    if (pkg.modelreq > 0) {
        pprintf(prn, " model-requirement=\"%s\"",
                gretl_command_word(pkg.modelreq));
    }
    if (pkg.minver > 0) {
        char[10] vstr;

        pprintf(prn, " minver=\"%s\"",
                gretl_version_string(vstr.ptr, pkg.minver));
    }
    if (pkg.uses_subdir) {
        pputs(prn, " lives-in-subdir=\"true\"");
    }
    if (pkg.data_access) {
        pputs(prn, " wants-data-access=\"true\"");
    }

    pputs(prn, ">\n");

    gretl_xml_put_tagged_string("author",  pkg.author, prn);
    gretl_xml_put_tagged_string("version", pkg.version_, prn);
    gretl_xml_put_tagged_string("date",    pkg.date, prn);
    gretl_xml_put_tagged_string("description", pkg.descrip, prn);

    if (pkg.help_fname != null) {
	gretl_xml_put_tagged_string("helpfile", pkg.help_fname, prn);
    } else if (pkg.help != null && !strncmp(pkg.help, "pdfdoc:", 7)) {
	gretl_xml_put_tagged_string("helpfile", pkg.help + 7, prn);
    }

    if (pkg.tags != null) {
        gretl_xml_put_tagged_string("tags", pkg.tags, prn);
    }
    if (pkg.label != null) {
        gretl_xml_put_tagged_string("label", pkg.label, prn);
    }
    if (pkg.mpath != null) {
        gretl_xml_put_tagged_string("menu-attachment", pkg.mpath, prn);
    }

    pputs(prn, "</gretl-addon>\n");

    gretl_print_destroy(prn);

    pprintf(inprn, "Wrote index file %s\n", idxname);
    g_free(idxname);

    return 0;
}

/* Write out a function package as XML. If @fp is null, then we write
   this as a package file in its own right, using the filename given
   in the package (the 'standalone' case), otherwise we're writing it
   as a component of the functions listing in a gretl session file,
   which already has an associated open stream.
*/

int real_write_function_package(Fnpkg * pkg, Prn * prn, int mpi)
{
    int standalone = (prn == null);
    int i;
    int err = 0;

    if (standalone) {
        prn = gretl_print_new_with_filename(pkg.fname, &err);
        if (prn == null) {
            return err;
        } else {
            gretl_xml_header(prn);
            pputs(prn, "<gretl-functions>\n");
        }
    }

    pputs(prn, "<gretl-function-package");

    if (pkg.name[0] == '\0') {
        name_package_from_filename(pkg);
    }

    pprintf(prn, " name=\"%s\"", pkg.name);

    if (pkg.dreq == DataReq.needs_ts) {
        pprintf(prn, " %s=\"true\"".ptr, NEEDS_TS.ptr);
    } else if (pkg.dreq == DataReq.needs_qm) {
        pprintf(prn, " %s=\"true\"".ptr, NEEDS_QM.ptr);
    } else if (pkg.dreq == DataReq.needs_panel) {
        pprintf(prn, " %s=\"true\"".ptr, NEEDS_PANEL.ptr);
    } else if (pkg.dreq == DataReq.nodata_ok) {
        pprintf(prn, " %s=\"true\"".ptr, NO_DATA_OK.ptr);
    }

    if (pkg.modelreq > 0) {
        pprintf(prn, " model-requirement=\"%s\"",
                gretl_command_word(pkg.modelreq));
    }
    if (pkg.minver > 0) {
        char[10] vstr;

        pprintf(prn, " minver=\"%s\"", gretl_version_string(vstr.ptr, pkg.minver));
    }
    if (pkg.uses_subdir) {
        pprintf(prn, " lives-in-subdir=\"true\"");
    }
    if (pkg.data_access) {
        pprintf(prn, " wants-data-access=\"true\"");
    }

    pputs(prn, ">\n");

    if (pkg.email != null && *pkg.email != '\0') {
        gretl_xml_put_tagged_string_plus("author", pkg.author,
                                         "email", pkg.email,
                                         prn);
    } else {
        gretl_xml_put_tagged_string("author", pkg.author, prn);
    }
    gretl_xml_put_tagged_string("version", pkg.version_, prn);
    gretl_xml_put_tagged_string("date",    pkg.date, prn);
    gretl_xml_put_tagged_string("description", pkg.descrip, prn);

    if (pkg.tags != null) {
        gretl_xml_put_tagged_string("tags", pkg.tags, prn);
    }
    if (pkg.label != null) {
        gretl_xml_put_tagged_string("label", pkg.label, prn);
    }
    if (pkg.mpath != null) {
        gretl_xml_put_tagged_string("menu-attachment", pkg.mpath, prn);
    }
    if (pkg.help != null && !mpi) {
        if (pkg.help_fname != null) {
            pprintf(prn, "<help filename=\"%s\">\n", pkg.help_fname);
        } else {
            pputs(prn, "<help>\n");
        }
        gretl_xml_put_string(trim_text(pkg.help), prn);
        pputs(prn, "\n</help>\n");
    }
    if (pkg.gui_help != null && !mpi) {
        if (pkg.gui_help_fname != null) {
            pprintf(prn, "<gui-help filename=\"%s\">\n",
                              pkg.gui_help_fname);
        } else {
            pputs(prn, "<gui-help>\n");
        }
        gretl_xml_put_string(trim_text(pkg.gui_help), prn);
        pputs(prn, "\n</gui-help>\n");
    }

    if (pkg.datafiles != null) {
        gretl_xml_put_strings_array("data-files",
                                    cast(const(char) **) pkg.datafiles,
                                    pkg.n_files, prn);
    }
    if (pkg.depends != null) {
        gretl_xml_put_strings_array("depends",
                                    cast(const(char) **) pkg.depends,
                                    pkg.n_depends, prn);
    }
    if (pkg.provider != null) {
        gretl_xml_put_tagged_string("provider", pkg.provider, prn);
    }
    if (pkg.conflict != null) {
        gretl_xml_put_tagged_string("conflict", pkg.conflict, prn);
    }
    if (pkg.Rdeps != null) {
        gretl_xml_put_tagged_string("R-depends", pkg.Rdeps, prn);
    }
    if (pkg.pub != null) {
        for (i=0; i<pkg.n_pub; i++) {
            write_function_xml(pkg.pub[i], prn, mpi);
        }
    }
    if (pkg.priv != null) {
        for (i=0; i<pkg.n_priv; i++) {
            write_function_xml(pkg.priv[i], prn, mpi);
        }
    }
    if (pkg.sample != null && !mpi) {
        if (pkg.sample_fname != null) {
            pprintf(prn, "<sample-script filename=\"%s\">\n",
                              pkg.sample_fname);
        } else {
            pputs(prn, "<sample-script>\n");
        }
        gretl_xml_put_string(trim_text(pkg.sample), prn);
        pputs(prn, "\n</sample-script>\n");
    }

    pputs(prn, "</gretl-function-package>\n");

    if (standalone) {
        pputs(prn, "</gretl-functions>\n");
        gretl_print_destroy(prn);
    }

    return err;
}

/**
 * function_package_write_file:
 * @pkg: function package.
 *
 * Write out @pkg as an XML file, using the filename
 * recorded in the package.
 *
 * Returns: 0 on success, non-zero on error.
 */

extern (C) int function_package_write_file(Fnpkg * pkg)
{
    return real_write_function_package(pkg, null, 0);
}

/* below: apparatus for constructing and saving a gfn function
   package from the command line
*/

int is_unclaimed(const(char) * s, char ** S, int n)
{
    int i;

    for (i=0; i<n; i++) {
        if (strcmp(s, S[i]) == 0) {
            return 0;
        }
    }

    return 1;
}

void strip_cr(char * s)
{
    while (*s) {
        if (*s == 0x0d && *(s+1) == 0x0a) {
            /* CR + LF . LF */
            memmove(s, s+1, strlen(s));
            s++;
        }
        s++;
    }
}

char * pkg_aux_content(const(char) * fname, int * err)
{
    char * ret = null;
    GError *gerr = null;
    gsize len = 0;

    g_file_get_contents(fname, &ret, &len, &gerr);

    if (gerr != null) {
        gretl_errmsg_set(gerr.message);
        g_error_free(gerr);
        *err = GretlError.fopen;
    } else if (strchr(ret, '\r')) {
        strip_cr(ret);
    }

    return ret;
}

int pkg_set_dreq(Fnpkg * pkg, const(char) * s)
{
    int err = 0;

    if (!strcmp(s, NEEDS_TS)) {
        pkg.dreq = DataReq.needs_ts;
    } else if (!strcmp(s, NEEDS_QM)) {
        pkg.dreq = DataReq.needs_qm;
    } else if (!strcmp(s, NEEDS_PANEL)) {
        pkg.dreq = DataReq.needs_panel;
    } else if (!strcmp(s, NO_DATA_OK)) {
        pkg.dreq = DataReq.nodata_ok;
    } else {
        err = GretlError.parse;
    }

    return err;
}

int pkg_set_modelreq(Fnpkg * pkg, const(char) * s)
{
    int ci = gretl_command_number(s);

    if (ci > 0) {
        pkg.modelreq = ci;
        return 0;
    } else {
        return GretlError.parse;
    }
}

int pkg_set_datafiles(Fnpkg * pkg, const(char) * s)
{
    int err = 0;

    if (string_is_blank(s)) {
        err = GretlError.data;
    } else {
        int n = 0;

        pkg.datafiles = gretl_string_split(s, &n, null);
        if (pkg.datafiles == null) {
            pkg.n_files = 0;
            err = GretlError.alloc;
        } else {
            pkg.n_files = n;
        }
    }

    if (!err) {
        pkg.uses_subdir = 1;
    }

    return err;
}

int pkg_set_depends(Fnpkg * pkg, const(char) * s)
{
    int err = 0;

    if (string_is_blank(s)) {
        err = GretlError.data;
    } else {
        int n = 0;

        pkg.depends = gretl_string_split(s, &n, null);
        if (pkg.depends == null) {
            pkg.n_depends = 0;
            err = GretlError.alloc;
        } else {
            pkg.n_depends = n;
        }
    }

    return err;
}

int pkg_boolean_from_string(const(char) * s)
{
    if (!strcmp(s, "true")) {
        return 1;
    } else {
        return 0;
    }
}

int pkg_remove_role(Fnpkg * pkg, int role)
{
    Ufunc * u;
    int i;

    for (i=0; i<pkg.n_priv; i++) {
        u = pkg.priv[i];
        if (u.pkg_role == role) {
            u.pkg_role = UfunRole.none;
            return 0;
        }
    }
    for (i=0; i<pkg.n_pub; i++) {
        u = pkg.pub[i];
        if (u.pkg_role == role) {
            u.pkg_role = UfunRole.none;
            return 0;
        }
    }

    return GretlError.data;
}

int valid_list_maker(Ufunc * u)
{
    if (u.n_params == 0) {
        return 1; /* OK */
    } else if (u.n_params == 1 &&
               u.params[0].type == GretlType.list) {
        return 1; /* OK */
    } else {
        return 0;
    }
}

int valid_plot_checker(Ufunc * u)
{
    return (u.n_params == 1 &&
	    u.params[0].type == GretlType.bundle_ref);
}

bool PRIVATE_ONLY(int r)
{
    return r == UfunRole.gui_precheck || r == UfunRole.plot_precheck || r == UfunRole.list_maker || r == UfunRole.r_setup || r == UfunRole.ui_maker;
}

Ufunc * get_package_member_by_name(const(char) * name, Fnpkg * pkg, int public_ok, int private_ok)
{
    int i;

    if (public_ok) {
        for (i=0; i<pkg.n_pub; i++) {
            if (!strcmp(name, pkg.pub[i].name)) {
                return pkg.pub[i];
            }
        }
    }
    if (private_ok) {
        for (i=0; i<pkg.n_priv; i++) {
            if (!strcmp(name, pkg.priv[i].name)) {
                return pkg.priv[i];
            }
        }
    }

    return null;
}

extern (C) int function_set_package_role(const(char) * name, Fnpkg * pkg, const(char) * attr, Prn * prn)
{
    Ufunc * u = null;
    int role = pkg_key_get_role(attr);
    int err = 0;

    if (name == null) {
        /* removing a role altogether */
        pkg_remove_role(pkg, role);
        return 0;
    } else if (role == UfunRole.none) {
        /* canceling a function's role */
        u = get_package_member_by_name(name, pkg, 1, 1);
        if (u != null) {
            u.pkg_role = role;
            return 0;
        } else {
            return GretlError.data;
        }
    }

    /* Check that the function identified by @name satisfies the
       requirements for the role identified by @attr: if so, hook
       the function up in this role.
    */

    if (PRIVATE_ONLY(role)) {
        /* function roles that must be private */
        u = get_package_member_by_name(name, pkg, 0, 1);
        if (u == null) {
            pprintf(prn, "%s: %s: no such private function\n", attr, name);
            return GretlError.data;
        }
        if (role == UfunRole.gui_precheck) {
            if (u.rettype != GretlType.double_) {
                pprintf(prn, "%s: must return a scalar\n", attr);
                err = GretlError.types;
            } else if (u.n_params > 0) {
                pprintf(prn, "%s: no parameters are allowed\n", attr);
                err = GretlError.types;
            }
	} else if (role == UfunRole.plot_precheck) {
            if (u.rettype != GretlType.matrix) {
                pprintf(prn, "%s: must return a matrix (vector)\n", attr);
                err = GretlError.types;
            } else if (!valid_plot_checker(u)) {
                pprintf(prn, "%s: must have a single bundle parameter\n", attr);
                err = GretlError.types;
            }
        } else if (role == UfunRole.list_maker) {
            if (u.rettype != GretlType.list) {
                pprintf(prn, "%s: must return a list\n", attr);
                err = GretlError.types;
            } else if (!valid_list_maker(u)) {
                pprintf(prn, "%s: must have 0 parameters or a single " ~
		"list parameter\n", attr);
                err = GretlError.types;
            }
        } else if (role == UfunRole.r_setup) {
            if (u.rettype != GretlType.void_) {
                pprintf(prn, "%s: should not return anything\n", attr);
                err = GretlError.types;
            } else if (u.n_params > 0) {
                pprintf(prn, "%s: no parameters are allowed\n", attr);
                err = GretlError.types;
            }
        } else if (role == UfunRole.ui_maker) {
            if (u.rettype != GretlType.bundle) {
                pprintf(prn, "%s: should return a bundle\n", attr);
                err = GretlError.types;
            } else if (u.n_params > 0) {
                pprintf(prn, "%s: no parameters are allowed\n", attr);
                err = GretlError.types;
            }
        }
        if (!err) {
            u.pkg_role = role;
        }
        return err; /* private cases handled */
    }

    /* all other special-role functions may be public */
    u = get_package_member_by_name(name, pkg, 1, 1);
    if (u == null) {
        pprintf(prn, "%s: %s: no such function\n", attr, name);
        return GretlError.data;
    }

    if (role == UfunRole.gui_main) {
        ; /* OK, type doesn't matter */
    } else {
        /* bundle-print, bundle-plot, etc: these cases require
           some further checks
        */
        int fcast = (role == UfunRole.bundle_fcast);
        int j;
        int pmin = fcast ? 2 : 1;
        GretlType pjt;

        if (u.n_params == 0) {
            pprintf(prn, "%s: must take a %s argument\n", attr,
                    gretl_type_get_name(GretlType.bundle_ref));
            err = GretlError.types;
        }
        for (j=0; j<u.n_params && !err; j++) {
            pjt = u.params[j].type;
            if (j == 0 && pjt != GretlType.bundle_ref) {
                pprintf(prn, "%s: first param type must be %s\n",
                        attr, gretl_type_get_name(GretlType.bundle_ref));
                err = GretlError.types;
            } else if (j == 1 && pjt != GretlType.int_) {
                pprintf(prn, "%s: second param type must be %s\n",
                        attr, gretl_type_get_name(GretlType.int_));
                err = GretlError.types;
            } else if (j == 2 && fcast && pjt != GretlType.int_) {
                pprintf(prn, "%s: third param type must be %s\n",
                        attr, gretl_type_get_name(GretlType.int_));
                err = GretlError.types;
            } else if (j > pmin && !fn_param_optional(u, j) &&
                       na(fn_param_default(u, j))) {
                pprintf(prn, "%s: extra params must be optional\n",
                        attr);
                err = GretlError.types;
            }
        }
    }

    if (!err) {
        u.pkg_role = role;
    }

    return err;
}

/* Called from the GUI package editor to check whether
   a given function of name @name can be shown as a
   candidate for a specified GUI-special @role.
*/

extern (C) int function_ok_for_package_role(const(char) * name, int role)
{
    Ufunc * u = null;
    int i;
    int err = 0;

    if (name == null || role == UfunRole.none) {
        return 0;
    }

    for (i=0; i<n_ufuns; i++) {
        if (!strcmp(name, ufuns[i].name)) {
            u = ufuns[i];
            break;
        }
    }

    if (u == null) {
        return 0;
    }

    if (role == UfunRole.gui_precheck) {
        if (u.rettype != GretlType.double_) {
            err = GretlError.types;
        } else if (u.n_params > 0) {
            err = GretlError.types;
        }
        return !err; /* found */
    } else if (role == UfunRole.plot_precheck) {
	if (u.rettype != GretlType.matrix) {
	    err = GretlError.types;
	} else if (!valid_plot_checker(u)) {
	    err = GretlError.types;
	}
	return !err; /* found */
    } else if (role == UfunRole.list_maker) {
        if (u.rettype != GretlType.list) {
            err = GretlError.types;
        } else if (!valid_list_maker(u)) {
            err = GretlError.types;
        }
        return !err; /* found */
    } else if (role == UfunRole.r_setup) {
        if (u.rettype != GretlType.void_) {
            err = GretlError.types;
        } else if (u.n_params > 0) {
            err = GretlError.types;
        }
        return !err; /* found */
    } else if (role == UfunRole.ui_maker) {
        if (u.rettype != GretlType.bundle) {
            err = GretlError.types;
        } else if (u.n_params > 0) {
            err = GretlError.types;
        }
        return !err; /* found */
    }

    if (role == UfunRole.gui_main) {
        ; /* OK, we don't mind what type it is */
    } else {
        /* bundle-print, bundle-plot, etc. */
        if (u.n_params == 0) {
            err = GretlError.types;
        }
        for (i=0; i<u.n_params && !err; i++) {
            if (i == 0 && u.params[i].type != GretlType.bundle_ref) {
                err = GretlError.types;
            } else if (i == 1 && u.params[i].type != GretlType.int_) {
                err = GretlError.types;
            } else if (i > 1 && !fn_param_optional(u, i) &&
                       na(fn_param_default(u, i))) {
                err = GretlError.types;
            }
        }
    }

    return !err;
}

int pkg_set_funcs_attribute(Fnpkg * pkg, const(char) * s, int flag)
{
    char **S;
    int ns = 0;
    int err = 0;

    S = gretl_string_split(s, &ns, null);
    if (ns == 0) {
        err = GretlError.data;
    } else {
        int i;
        int j;
        int match;

        for (i=0; i<ns && !err; i++) {
            match = 0;
            for (j=0; j<pkg.n_pub; j++) {
                if (!strcmp(S[i], pkg.pub[j].name)) {
                    pkg.pub[j].flags |= flag;
                    match = 1;
                    break;
                }
            }
            if (!match) {
                err = GretlError.data;
            }
        }

        strings_array_free(S, ns);
    }

    return err;
}

void function_package_set_auxfile(Fnpkg * pkg, const(char) * id, const(char) * fname)
{
    char * test = null;

    /* Maybe set the source filename for an element
       of a function package read from file, but only
       if it's not the standard, default filename for
       the element in question.
    */

    if (!strcmp(id, "help-fname")) {
        test = g_strdup_printf("%s_help.txt", pkg.name);
    } else if (!strcmp(id, "gui-help-fname")) {
        test = g_strdup_printf("%s_gui_help.txt", pkg.name);
    } else if (!strcmp(id, "sample-fname")) {
        test = g_strdup_printf("%s_sample.inp", pkg.name);
    }

    if (test != null) {
        if (strcmp(fname, test)) {
           function_package_set_properties(pkg, id, fname, null);
        }
        g_free(test);
    }
}

int check_for_pdf_ref(const(char) * s, Fnpkg * pkg, int * err)
{
    int len;
    int ret = 0;

    if (!strncmp(s, "pdfdoc:", 7)) {
        s += 7;
    }

    len = cast(int) strlen(s);
    if (len < 64 && strchr(s, ' ') == null && has_suffix(s, ".pdf")) {
        char *tmp = gretl_strndup(s, len - 4);

        if (strcmp(tmp, pkg.name)) {
            gretl_errmsg_sprintf(_("PDF doc should be called %s.pdf").ptr, pkg.name);
            *err = GretlError.data;
        } else {
            ret = 1;
        }
        free(tmp);
    }

    return ret;
}

int is_pdf_ref(const(char) * s)
{
    if (!strncmp(s, "pdfdoc:", 7)) {
        s += 7;
    }
    return strlen(s) < 64 && strchr(s, ' ') == null &&
        has_suffix(s, ".pdf");
}

int validate_R_depends(const(char) * s)
{
    int err = 0;

    if (strncmp(s, "R ", 2)) {
        err = GretlError.data;
    } else {
        /* skip "R " */
        s += 2;
        s += strcspn(s, " ");
        if (*s == '\0') {
            err = GretlError.data;
        } else {
            /* skip R version */
            s += strcspn(s, " ");
            while (*s && !err) {
                s += strspn(s, " ");
                if (*s == '\0') break;
                /* R package? */
                s += strcspn(s, " ");
                s += strspn(s, " ");
                if (*s == '\0') {
                    /* no version given? */
                    err = GretlError.data;
                } else {
                    /* skip package version */
                    s += strcspn(s, " ");
                }
            }
        }
    }

    if (err) {
        gretl_errmsg_set(_("Invalid R-depends line").ptr);
    }

    return err;
}

int function_package_set_version(Fnpkg * pkg, const(char) * vstr)
{
    int err = 0;

    if (pkg.version_ != null) {
	free(pkg.version_);
	pkg.version_ = null;
    }
    if (!strcmp(vstr, "@VERSION@")) {
        /* auto-versioning for gretl addons */
	if (is_gretl_addon(pkg.name)) {
	    pkg.version_ = gretl_strdup(GRETL_VERSION);
	} else {
	    err = GretlError.data;
	}
    } else {
	pkg.version_ = gretl_strdup(vstr);
    }

    return err;
}

/* Having assembled and checked the function-listing for a new
   package, now retrieve the additional information from the
   spec file (named by @fname, opened as @fp).
*/

int new_package_info_from_spec(Fnpkg * pkg, const(char) * fname, FILE * fp, Gretlopt opt, Prn * prn)
{
    const(char) * okstr;
    const(char) * failed;
    char * currdir = null;
    int quiet = (opt & Gretlopt.OPT_Q);
    char * p;
    char[1024] line;
    int got = 0;
    int err = 0;

    if (opt & Gretlopt.OPT_G) {
        /* GUI use */
        okstr = "<@ok>\n";
        failed = "<@fail>\n";
    } else {
        okstr = "OK\n";
        failed = "failed\n";
    }

    if (strrslash(cast(char *) fname) != null) {
        /* directory change needed */
        char[FILENAME_MAX] dirname;

        strcpy(dirname.ptr, fname);
        p = strrslash(dirname.ptr);
        p[0] = '\0';
        currdir = g_get_current_dir();
        err = gretl_chdir(dirname.ptr);
    }


    while (fgets(line.ptr, line.sizeof, fp) && !err) {
        if (*line.ptr == '#' || string_is_blank(line.ptr)) {
            continue;
        }
        tailstrip(line.ptr);
        p = strchr(line.ptr, '=');
        if (p == null) {
            continue;
        } else {
            p++;
            p += strspn(p, " ");
	    if (!strncmp(line.ptr, "public".ptr, 6)) {
		continue;
	    } else if (!strncmp(line.ptr, "author".ptr, 6)) {
                err = function_package_set_properties(pkg, "author".ptr, p, null);
                if (!err) got++;
            } else if (!strncmp(line.ptr, "email".ptr, 5)) {
                err = function_package_set_properties(pkg, "email".ptr, p, null);
            } else if (!strncmp(line.ptr, "version".ptr, 7)) {
                err = function_package_set_version(pkg, p);
                if (!err) got++;
            } else if (!strncmp(line.ptr, "date".ptr, 4)) {
                err = function_package_set_properties(pkg, "date".ptr, p, null);
                if (!err) got++;
            } else if (!strncmp(line.ptr, "description".ptr, 11)) {
                err = function_package_set_properties(pkg, "description".ptr, p, null);
                if (!err) got++;
            } else if (!strncmp(line.ptr, "tags".ptr, 4)) {
                err = function_package_set_properties(pkg, "tags".ptr, p, null);
            } else if (!strncmp(line.ptr, "label".ptr, 5)) {
                err = function_package_set_properties(pkg, "label".ptr, p, null);
            } else if (!strncmp(line.ptr, "menu-attachment".ptr, 15)) {
                err = function_package_set_properties(pkg, "menu-attachment".ptr, p, null);
            } else if (!strncmp(line.ptr, "provider".ptr, 8)) {
                if (!quiet) {
                    pprintf(prn, "Recording provider %s\n", p);
                }
                err = function_package_set_properties(pkg, "provider".ptr, p, null);
            } else if (!strncmp(line.ptr, "conflict".ptr, 8)) {
                if (!quiet) {
                    pprintf(prn, "Recording conflict %s\n", p);
                }
                err = function_package_set_properties(pkg, "conflict".ptr, p, null);
            } else if (!strncmp(line.ptr, "help".ptr, 4)) {
                char * hstr = null;
                int pdfdoc;

                pdfdoc = check_for_pdf_ref(p, pkg, &err);
                if (pdfdoc) {
                    if (!quiet) {
                        pprintf(prn, "Recording help reference %s\n", p);
                    }
                    hstr = g_strdup_printf("pdfdoc:%s", p);
                } else if (!err) {
                    if (!quiet) {
                        pprintf(prn, "Looking for help text in %s... ", p);
                    }
                    hstr = pkg_aux_content(p, &err);
                    if (err) {
                        pputs(prn, failed);
                    } else {
                        pputs(prn, okstr);
                    }
                }
                if (!err) {
                    err = function_package_set_properties(pkg, "help".ptr, hstr, null);
                    if (!err) {
                        got++;
                        if (!pdfdoc) {
                            function_package_set_auxfile(pkg, "help-fname", p);
                        }
                    }
                }
                g_free(hstr);
            } else if (!strncmp(line.ptr, "gui-help".ptr, 8)) {
                char * ghstr = null;

                if (!quiet) {
                    pprintf(prn, "Looking for GUI help text in %s... ", p);
                }
                ghstr = pkg_aux_content(p, &err);
                if (err) {
                    pputs(prn, failed);
                } else {
                    pputs(prn, okstr);
                    err = function_package_set_properties(pkg, "gui-help".ptr, ghstr, null);
                    if (!err) {
                        function_package_set_auxfile(pkg, "gui-help-fname", p);
                    }
                }
                g_free(ghstr);
            } else if (!strncmp(line.ptr, "sample-script".ptr, 13)) {
                char * script = null;

                if (!quiet) {
                    pprintf(prn, "Looking for sample script in %s... ", p);
                }
                script = pkg_aux_content(p, &err);
                if (err) {
                    pputs(prn, failed);
                } else {
                    pputs(prn, okstr);
                    err = function_package_set_properties(pkg, "sample-script".ptr, script, null);
                    if (!err) {
                        got++;
                        function_package_set_auxfile(pkg, "sample-fname", p);
                    }
                }
                g_free(script);
            } else if (!strncmp(line.ptr, "data-files".ptr, 10)) {
                if (!quiet) {
                    pprintf(prn, "Recording data-file list: %s\n", p);
                }
                err = pkg_set_datafiles(pkg, p);
            } else if (!strncmp(line.ptr, "depends".ptr, 7)) {
                if (!quiet) {
                    pprintf(prn, "Recording dependency list: %s\n", p);
                }
                err = pkg_set_depends(pkg, p);
            } else if (!strncmp(line.ptr, "R-depends".ptr, 9)) {
                err = validate_R_depends(p);
                if (!err) {
                    if (!quiet) {
                        pprintf(prn, "Recording R dependency list: %s\n", p);
                    }
                    err = function_package_set_properties(pkg, "R-depends".ptr, p, null);
                }
            } else if (!strncmp(line.ptr, "data-requirement".ptr, 16)) {
                err = pkg_set_dreq(pkg, p);
            } else if (!strncmp(line.ptr, "model-requirement".ptr, 17)) {
                err = pkg_set_modelreq(pkg, p);
            } else if (!strncmp(line.ptr, "min-version".ptr, 11)) {
                pkg.minver = gretl_version_number(p);
                got++;
            } else if (!strncmp(line.ptr, "lives-in-subdir".ptr, 15)) {
                pkg.uses_subdir = cast(ubyte) pkg_boolean_from_string(p);
            } else if (!strncmp(line.ptr, "wants-data-access".ptr, 17)) {
                pkg.data_access = cast(ubyte) pkg_boolean_from_string(p);
            } else if (!strncmp(line.ptr, "no-print".ptr, 8)) {
                err = pkg_set_funcs_attribute(pkg, p, UfunAttrs.noprint);
            } else if (!strncmp(line.ptr, "menu-only".ptr, 9)) {
                err = pkg_set_funcs_attribute(pkg, p, UfunAttrs.menu_only);
            } else if (!strncmp(line.ptr, "must-assign".ptr, 11)) {
                err = pkg_set_funcs_attribute(pkg, p, UfunAttrs.assign);
            } else {
                const(char) * key;
		int found = 0;
                int i;

                for (i=0; pkg_lookups[i].key != null; i++) {
                    key = pkg_lookups[i].key;
                    if (!strncmp(line.ptr, key, strlen(key))) {
                        err = function_set_package_role(p, pkg, key, prn);
                        if (!err && !quiet) {
                            pprintf(prn, "%s function is %s, %s", key, p, okstr);
                        }
			found = 1;
                        break;
                    }
                }
		if (!found) {
		    fprintf(stderr, "*** unrecognized '%s' ***\n".ptr, line.ptr);
		}
            }
        }
    }

    if (!err && pkg.provider != null) {
        /* in case provider isn't registered as a dependency */
        err = strings_array_prepend_uniq(&pkg.depends, &pkg.n_depends,
                                         pkg.provider);
    }

    if (currdir != null) {
        /* go back where we came from */
        gretl_chdir(currdir);
        g_free(currdir);
    }

    if (!err && got < 7) {
        gretl_errmsg_set(_("Some required information was missing").ptr);
        err = GretlError.data;
    }

    return err;
}

int is_public_funcs_line(const(char) * line, int * offset)
{
    const(char) * s = line;

    if (!strncmp(s, "public", 6)) {
        s += 6;
        s += strspn(s, " ");
        if (*s == '=') {
            *offset = cast(int) (1 + s - line);
            return 1;
        }
    }

    return 0;
}

Fnpkg * new_pkg_from_spec_file(const(char) * gfnname, Gretlopt opt, Prn * prn, int * err)
{
    Fnpkg * pkg = null;
    char[FILENAME_MAX] fname;
    char[4096] line;
    char[1024] cont;
    int quiet = (opt & Gretlopt.OPT_Q);
    FILE * fp;

    if (!has_suffix(gfnname, ".gfn")) {
        gretl_errmsg_set(_("Output must have extension \".gfn\"").ptr);
        *err = GretlError.data;
        return null;
    }

    switch_ext(fname.ptr, gfnname, "spec");
    fp = gretl_fopen(fname.ptr, "r"); /* "rb" ? */

    if (fp == null) {
        *err = GretlError.fopen;
    } else {
        char **pubnames = null;
        char **privnames = null;
        int npub = 0;
        int npriv = 0;
        Ufunc * uf;
        int i;

        if (!quiet) {
            pprintf(prn, _("Found spec file '%s'\n").ptr, fname.ptr);
        }

        /* first pass: gather names of public functions */

        while (fgets(line.ptr, line.sizeof, fp) && !*err) {
            int offset = 0;

            if (is_public_funcs_line(line.ptr, &offset)) {
                while (ends_with_backslash(line.ptr)) {
                    gretl_charsub(line.ptr, '\\', '\0');
                    cont[0] = '\0';
                    if (fgets(cont.ptr, cont.sizeof, fp)) {
                        strcat(line.ptr, cont.ptr);
                    } else {
                        *err = GretlError.parse;
                    }
                }
                if (!*err) {
                    tailstrip(line.ptr);
                    pubnames = gretl_string_split(line.ptr + offset, &npub, null);
                }
                break;
            }
        }

        if (npub == 0) {
            pprintf(prn, _("no specification of public functions was found\n").ptr);
            *err = GretlError.data;
        }

        if (!*err) {
            if (!quiet) {
                pprintf(prn, _("number of public interfaces = %d\n").ptr, npub);
            }
            for (i=0; i<npub && !*err; i++) {
                uf = get_user_function_by_name(pubnames[i]);
                if (!quiet) {
                    pprintf(prn, " %s", pubnames[i]);
                }
                if (uf == null) {
                    if (!quiet) {
                        pputs(prn, _(": *** not found").ptr);
                    }
                    gretl_errmsg_sprintf(_("'%s': no such function").ptr, pubnames[i]);
                    *err = GretlError.data;
                }
                if (!quiet) {
                    pputc(prn, '\n');
                }
            }
        }

        /* note: any other functions currently loaded are assumed to be
           private functions for this package */

        if (!*err) {
            npriv = n_free_functions() - npub;
            if (npriv < 0) {
                npriv = 0;
            }
        }

        if (!*err && npriv > 0) {
            npriv = 0;
            for (i=0; i<n_ufuns && !*err; i++) {
                if (ufuns[i].pkg == null && is_unclaimed(ufuns[i].name, pubnames, npub)) {
                    *err = strings_array_add(&privnames, &npriv, ufuns[i].name);
                }
            }
        }

        if (!quiet && !*err && npriv > 0) {
            pprintf(prn, "number of private functions = %d\n", npriv);
            for (i=0; i<npriv; i++) {
                pprintf(prn, " %s\n", privnames[i]);
            }
        }

        if (!*err) {
            pkg = function_package_new(gfnname, pubnames, npub,
                                       privnames, npriv, err);
        }

        strings_array_free(pubnames, npub);
        strings_array_free(privnames, npriv);

        if (!*err) {
            rewind(fp);
            *err = new_package_info_from_spec(pkg, fname.ptr, fp, opt, prn);
        }

        if (*err && pkg != null) {
            real_function_package_unload(pkg, 0);
            pkg = null;
        }

        fclose(fp);
    }

    return pkg;
}

int cli_validate_package_file(const(char) * fname, Gretlopt opt, Prn * prn)
{
    char[FILENAME_MAX] dtdname;
    xmlDocPtr doc = null;
    xmlDtdPtr dtd = null;
    int err;

    err = gretl_xml_open_doc_root(fname, null, &doc, null);
    if (err) {
        pprintf(prn, "Couldn't parse %s\n", fname);
        return 1;
    }

    dtdname[0] = '\0';

    if (opt & Gretlopt.OPT_D) {
        const(char) * dpath = get_optval_string(GretlCmdIndex.makepkg, Gretlopt.OPT_D);

        if (dpath != null && *dpath != '\0') {
            strcat(dtdname.ptr, dpath);
        }
    } else {
        sprintf(dtdname.ptr, "%sfunctions%cgretlfunc.dtd", gretl_home(), SLASH);
    }

    if (*dtdname.ptr != '\0') {
        dtd = xmlParseDTD(null, cast(const(xmlChar) *) dtdname.ptr);
    }

    if (dtd == null) {
        pputs(prn, "Couldn't open DTD to check package\n");
    } else {
        const(char) * pkgname = path_last_element(fname);
        xmlValidCtxtPtr cvp = xmlNewValidCtxt();

        if (cvp == null) {
            pputs(prn, "Couldn't get an XML validation context\n");
            xmlFreeDtd(dtd);
            xmlFreeDoc(doc);
            return 0;
        }

        cvp.userData = cast(void *) prn;
        cvp.error    = cast(xmlValidityErrorFunc) &pprintf2;
        cvp.warning  = cast(xmlValidityWarningFunc) &pprintf2;

        pprintf(prn, "Checking against %s\n", dtdname.ptr);

        if (!xmlValidateDtd(cvp, doc, dtd)) {
            err = 1;
        } else {
            pprintf(prn, _("%s: validated against DTD OK").ptr, pkgname);
            pputc(prn, '\n');
        }

        xmlFreeValidCtxt(cvp);
        xmlFreeDtd(dtd);
    }

    xmlFreeDoc(doc);

    return err;
}

int get_gfn_info_for_zip(const(char) * fname, int * pdfdoc, char *** datafiles, int * n_datafiles)
{
    xmlDocPtr doc = null;
    xmlNodePtr node = null;
    xmlNodePtr cur;
    int found = 0;
    int err = 0;

    node = gretl_xml_get_gfn(fname, &doc, &err);
    if (err) {
        return err;
    }

    cur = node.children;
    while (cur != null && found < 2) {
        if (!xmlStrcmp(cur.name, cast(XUC) "help".ptr)) {
            char * s = null;

            gretl_xml_node_get_trimmed_string(cur, doc, &s);
            *pdfdoc = is_pdf_ref(s);
            free(s);
            found++;
        } else if (!xmlStrcmp(cur.name, cast(XUC) "data-files".ptr)) {
            *datafiles =
                gretl_xml_get_strings_array(cur, doc, n_datafiles,
                                            0, &err);
            found++;
        } else if (!xmlStrcmp(cur.name, cast(XUC) "gretl-function".ptr)) {
            /* we've overshot */
            found = 2;
        }
        cur = cur.next;
    }

    if (doc != null) {
        xmlFreeDoc(doc);
    }

    return err;
}

int cli_build_zip_package(const(char) * fname, const(char) * gfnname, Fnpkg * pkg, Prn * prn)
{
    char **datafiles = null;
    int n_datafiles = 0;
    int pdfdoc = 0;
    int freeit = 0;
    int err = 0;

    if (pkg != null) {
        datafiles = pkg.datafiles;
        n_datafiles = pkg.n_files;
        pdfdoc = is_pdf_ref(pkg.help);
    } else {
        /* grope in the gfn file */
        err = get_gfn_info_for_zip(gfnname,
                                   &pdfdoc,
                                   &datafiles,
                                   &n_datafiles);
        freeit = 1;
    }

    if (!err) {
        err = package_make_zipfile(gfnname,
                                   pdfdoc,
                                   datafiles,
                                   n_datafiles,
                                   null,
                                   fname,
                                   Gretlopt.OPT_NONE,
                                   prn);
    }

    if (freeit) {
        strings_array_free(datafiles, n_datafiles);
    }

    return err;
}

int should_rebuild_gfn(const(char) * gfnname)
{
    char[FILENAME_MAX] testname;
    stat_t b1, b2, b3;
    int err;

    err = gretl_stat(gfnname, &b1);
    if (err) {
        /* gfn not found: so have to rebuild */
        return 1;
    }

    switch_ext(testname.ptr, gfnname, "inp");
    err = gretl_stat(testname.ptr, &b2);
    if (err) {
        /* no corresponding inp: can't rebuild */
        return 0;
    }

    switch_ext(testname.ptr, gfnname, "spec");
    err = gretl_stat(testname.ptr, &b3);
    if (err) {
        /* no corresponding spec: can't rebuild */
        return 0;
    }

    if (b2.st_mtime > b1.st_mtime ||
        b3.st_mtime > b1.st_mtime) {
        /* inp or spec is newer than gfn */
        return 1;
    }

    return 0;
}

/**
 * create_and_write_function_package:
 * @fname: filename for function package.
 * @opt: may include Gretlopt.OPT_I to write a package-index entry,
 * Gretlopt.OPT_T to write translatable strings (for addons).
 * @prn: printer struct for feedback.
 *
 * Create a package based on the functions currently loaded, and
 * write it out as an XML file. Responds to the makepkg command.
 *
 * Returns: 0 on success, non-zero on error.
 */

extern (C) int create_and_write_function_package(const(char) * fname, Gretlopt opt, Prn * prn)
{
    char[FILENAME_MAX] gfnname;
    Fnpkg * pkg = null;
    int build_gfn = 1;
    int build_zip = 0;
    int err = 0;

    if (has_suffix(fname, ".zip")) {
        /* building a zip package */
        switch_ext(gfnname.ptr, fname, "gfn");
        build_gfn = should_rebuild_gfn(gfnname.ptr);
        build_zip = 1;
    } else {
        /* just building a gfn file */
        strcpy(gfnname.ptr, fname);
    }

    if (build_gfn && n_free_functions() == 0) {
        gretl_errmsg_set(_("No functions are available for packaging at present.").ptr);
        err = GretlError.data;
    } else if (build_gfn) {
        pkg = new_pkg_from_spec_file(gfnname.ptr, opt, prn, &err);
        if (pkg != null) {
            err = function_package_write_file(pkg);
            if (!err) {
                err = cli_validate_package_file(gfnname.ptr, opt, prn);
                /* should we delete @gfnname ? */
            }
            if (!err) {
                if (opt & Gretlopt.OPT_T) {
                    package_write_translatable_strings(pkg, prn);
                }
                if (opt & Gretlopt.OPT_I) {
                    package_write_index(pkg, prn);
                }
            }
        }
    }

    if (!err && build_zip) {
        err = cli_build_zip_package(fname, gfnname.ptr, pkg, prn);
    }

    return err;
}

/**
 * function_package_get_name:
 * @pkg: function package.
 *
 * Returns: the name of the package.
 */

extern (C) const(char) * function_package_get_name(Fnpkg * pkg)
{
    return (pkg != null)? pkg.name : null;
}

/**
 * function_package_get_version:
 * @pkg: function package.
 *
 * Returns: the name of the package.
 */

extern (C) double function_package_get_version(Fnpkg * pkg)
{
    if (pkg == null) {
        return NADBL;
    } else {
        return dot_atof(pkg.version_);
    }
}

int maybe_replace_string_var(char ** svar, const(char) * src)
{
    if (src == null) {
        gretl_errmsg_set(_("string value is missing").ptr);
        return GretlError.data;
    } else {
        free(*svar);
        *svar = gretl_strdup(src);
        return (*svar == null)? GretlError.alloc : 0;
    }
}

/* unlike the case above, here we'll accept null for @src, to wipe
   out an existing string var */

int maybe_replace_optional_string_var(char ** svar, const(char) * src)
{
    if (src == null) {
        free(*svar);
        *svar = null;
        return 0;
    } else {
        free(*svar);
        *svar = gretl_strdup(src);
        return (*svar == null)? GretlError.alloc : 0;
    }
}

/* Called (indirectly) from GUI function packager. Note that we're
   careful not to touch UfunAttrs.private_ on the uf.flags side, since
   this flag is not represented in @attrs.
*/

void pkg_set_gui_attrs(Fnpkg * pkg, const(ubyte) * attrs)
{
    Ufunc * uf;
    int i;
    int r;

    for (r=1; r<UfunRole.gui_precheck; r++) {
        uf = null;
        for (i=0; i<n_ufuns; i++) {
            if (ufuns[i].pkg == pkg && ufuns[i].pkg_role == r) {
                uf = ufuns[i];
                break;
            }
        }
        if (uf != null) {
            if (attrs[r-1] & UfunAttrs.noprint) {
                uf.flags |= UfunAttrs.noprint;
            } else {
                uf.flags &= ~UfunAttrs.noprint;
            }
            if (attrs[r-1] & UfunAttrs.menu_only) {
                uf.flags |= UfunAttrs.menu_only;
            } else {
                uf.flags &= ~UfunAttrs.menu_only;
            }
        }
    }
}

/* Called (indirectly) from GUI function packager. Note that we scrub
   UfunAttrs.private_ from the user function flags passed to the packager:
   this flag is handled by a different mechanism.
*/

void pkg_get_gui_attrs(Fnpkg * pkg, ubyte * attrs)
{
    Ufunc * uf;
    int a;
    int i;
    int r;

    for (r=1; r<UfunRole.gui_precheck; r++) {
        uf = null;
        a = 0;
        for (i=0; i<n_ufuns; i++) {
            if (ufuns[i].pkg == pkg && ufuns[i].pkg_role == r) {
                uf = ufuns[i];
                break;
            }
        }
        if (uf != null) {
            if (uf.flags & UfunAttrs.noprint) {
                a |= UfunAttrs.noprint;
            }
            if (uf.flags & UfunAttrs.menu_only) {
                a |= UfunAttrs.menu_only;
            }
        }
        attrs[r-1] = cast(ubyte) a;
    }
}

char ** pkg_strvar_pointer(Fnpkg * pkg, const(char) * key, int * optional)
{
    *optional = 0;

    if (!strcmp(key, "fname")) {
        return &pkg.fname;
    } else if (!strcmp(key, "author")) {
        return &pkg.author;
    } else if (!strcmp(key, "email")) {
        return &pkg.email;
    } else if (!strcmp(key, "version")) {
        return &pkg.version_;
    } else if (!strcmp(key, "date")) {
        return &pkg.date;
    } else if (!strcmp(key, "description")) {
        return &pkg.descrip;
    } else if (!strcmp(key, "help")) {
        return &pkg.help;
    } else if (!strcmp(key, "sample-script")) {
        return &pkg.sample;
    }

    *optional = 1;

    if (!strcmp(key, "tags")) {
        return &pkg.tags; /* FIXME should be non-optional */
    } else if (!strcmp(key, "label")) {
        return &pkg.label;
    } else if (!strcmp(key, "menu-attachment")) {
        return &pkg.mpath;
    } else if (!strcmp(key, "gui-help")) {
        return &pkg.gui_help;
    } else if (!strcmp(key, "R-depends")) {
        return &pkg.Rdeps;
    } else if (!strcmp(key, "help-fname")) {
        return &pkg.help_fname;
    } else if (!strcmp(key, "gui-help-fname")) {
        return &pkg.gui_help_fname;
    } else if (!strcmp(key, "sample-fname")) {
        return &pkg.sample_fname;
    } else if (!strcmp(key, "provider")) {
        return &pkg.provider;
    } else if (!strcmp(key, "conflict")) {
	return &pkg.conflict;
    }

    return null;
}

/* varargs function for setting the properties of a function
   package: the settings take the form of a null-terminated
   set of (key, value) pairs.
*/

extern (C) int function_package_set_properties(Fnpkg * pkg, ...)
{
    va_list ap;
    const(char) * key;
    char **sptr;
    int optional;
    int err = 0;

    va_start(ap, pkg);

    while (err == 0) {
        key = va_arg!(const(char) *)(ap);
        if (key == null) {
            break;
        }

        sptr = pkg_strvar_pointer(pkg, key, &optional);

        if (sptr != null) {
            const(char) * sval = va_arg!(const(char) *)(ap);

            if (optional) {
                err = maybe_replace_optional_string_var(sptr, sval);
            } else {
                err = maybe_replace_string_var(sptr, sval);
            }
            if (!err && !strcmp(key, "help")) {
                if (!strncmp(sval, "pdfdoc", 6) ||
                    is_pdf_ref(sval)) {
                    pkg.uses_subdir = 1;
                }
            }
        } else if (!strcmp(key, "gui-attrs")) {
            const(ubyte) *np = va_arg!(const(ubyte) *)(ap);

            pkg_set_gui_attrs(pkg, np);
        } else {
            int ival = va_arg!(int)(ap);

            if (!strcmp(key, "data-requirement")) {
                pkg.dreq = cast(DataReq) ival;
            } else if (!strcmp(key, "model-requirement")) {
                pkg.modelreq = ival;
            } else if (!strcmp(key, "min-version")) {
                pkg.minver = ival;
            } else if (!strcmp(key, "lives-in-subdir")) {
                pkg.uses_subdir = (ival != 0);
            } else if (!strcmp(key, "wants-data-access")) {
                pkg.data_access = (ival != 0);
            }
        }
    }

    va_end(ap);

    return err;
}

enum PkgListCode {
    publist,
    guilist,
    privlist,
    assnlist
}

/* From a function package get a list of either its public or its
   private functions, in the form of a simple gretl list.  The
   elements of this list are just the positions of these functions in
   the current recorder array for user-functions.  This is fine if
   the information is to be used right away (as in the GUI function
   call dialog), but it should _not_ be assumed that the identifiers
   will remain valid indefinitely.
*/

int * function_package_get_list(Fnpkg * pkg, int code, int n)
{
    int * list = null;
    int j = 0;

    if (n > 0) {
        list = gretl_list_new(n);
        if (list != null) {
            int i;
            int priv;
            int menu_only;

            for (i=0; i<n_ufuns; i++) {
                if (ufuns[i].pkg == pkg) {
                    priv = function_is_private(ufuns[i]);
                    menu_only = function_is_menu_only(ufuns[i]);
                    if (code == PkgListCode.assnlist && !priv &&
                        function_must_assign(ufuns[i])) {
                        list[++j] = i;
                    }
                    if (code == PkgListCode.privlist && priv) {
                        list[++j] = i;
                    } else if (code == PkgListCode.publist && !priv) {
                        list[++j] = i;
                    } else if (code == PkgListCode.guilist && !priv && !menu_only &&
                               !pkg_aux_role(ufuns[i].pkg_role)) {
                        /* in the GUI list of public functions, don't
                           display post-processing functions
                        */
                        list[++j] = i;
                    }
                }
            }
        }
    }

    if (list != null) {
        if (j == 0) {
            free(list);
            list = null;
        } else {
            list[0] = j;
        }
    }

    return list;
}

char * pkg_get_special_func_name(Fnpkg * pkg, UfunRole role)
{
    int i;

    for (i=0; i<pkg.n_pub; i++) {
        if (pkg.pub[i].pkg_role == role) {
            return g_strdup(pkg.pub[i].name);
        }
    }
    for (i=0; i<pkg.n_priv; i++) {
        if (pkg.priv[i].pkg_role == role) {
            return g_strdup(pkg.priv[i].name);
        }
    }

    return null;
}

int pkg_get_special_func_id(Fnpkg * pkg, UfunRole role)
{
    int i;

    for (i=0; i<n_ufuns; i++) {
        if (ufuns[i].pkg == pkg && ufuns[i].pkg_role == role) {
            return i;
        }
    }

    return -1;
}

void handle_optional_string(char ** ps, const(char) * src)
{
    if (src == null) {
        *ps = null;
    } else {
        *ps = g_strdup(src);
    }
}

/* varargs function for retrieving the properties of a function
   package: the arguments after @pkg take the form of a
   null-terminated set of (key, pointer) pairs; values are written to
   the locations given by the pointers.
*/

extern (C) int function_package_get_properties(Fnpkg * pkg, ...)
{
    va_list ap;
    int npub = 0;
    int npriv = 0;
    int **plist;
    const(char) * key;
    int role;
    void * ptr;
    char **ps;
    int * pi;
    int i;
    int err = 0;

    if (pkg == null) {
        return GretlError.data;
    }

    for (i=0; i<n_ufuns; i++) {
        if (ufuns[i].pkg == pkg) {
            if (function_is_private(ufuns[i])) {
                npriv++;
            } else {
                npub++;
            }
        }
    }

    va_start(ap, pkg);

    for (i=1; ; i++) {
        key = va_arg!(const(char) *)(ap);
        if (key == null) {
            break;
        }
        ptr = va_arg!(void *)(ap);
        if (ptr == null) {
            break;
        }
        if (!strcmp(key, "name")) {
            ps = cast(char **) ptr;
            *ps = g_strdup(pkg.name);
        } else if (!strcmp(key, "author")) {
            ps = cast(char **) ptr;
            *ps = g_strdup(pkg.author);
        } else if (!strcmp(key, "email")) {
            ps = cast(char **) ptr;
            if (pkg.email != null) {
                *ps = g_strdup(pkg.email);
            } else {
                *ps = g_strdup(""); /* ? */
            }
        } else if (!strcmp(key, "date")) {
            ps = cast(char **) ptr;
            *ps = g_strdup(pkg.date);
        } else if (!strcmp(key, "version")) {
            ps = cast(char **) ptr;
            *ps = g_strdup(pkg.version_);
        } else if (!strcmp(key, "description")) {
            ps = cast(char **) ptr;
            *ps = g_strdup(pkg.descrip);
        } else if (!strcmp(key, "help")) {
            ps = cast(char **) ptr;
            *ps = g_strdup(pkg.help);
        } else if (!strcmp(key, "gui-help")) {
            ps = cast(char **) ptr;
            handle_optional_string(ps, pkg.gui_help);
        } else if (!strcmp(key, "R-depends")) {
            ps = cast(char **) ptr;
            handle_optional_string(ps, pkg.Rdeps);
        } else if (!strcmp(key, "sample-script")) {
            ps = cast(char **) ptr;
            *ps = g_strdup(pkg.sample);
        } else if (!strcmp(key, "help-fname")) {
            ps = cast(char **) ptr;
            handle_optional_string(ps, pkg.help_fname);
        } else if (!strcmp(key, "gui-help-fname")) {
            ps = cast(char **) ptr;
            handle_optional_string(ps, pkg.gui_help_fname);
        } else if (!strcmp(key, "sample-fname")) {
            ps = cast(char **) ptr;
            handle_optional_string(ps, pkg.sample_fname);
        } else if (!strcmp(key, "tags")) {
            ps = cast(char **) ptr;
            *ps = g_strdup(pkg.tags);
        } else if (!strcmp(key, "label")) {
            ps = cast(char **) ptr;
            *ps = g_strdup(pkg.label);
        } else if (!strcmp(key, "menu-attachment")) {
            ps = cast(char **) ptr;
            *ps = g_strdup(pkg.mpath);
        } else if (!strcmp(key, "provider")) {
            ps = cast(char **) ptr;
            *ps = g_strdup(pkg.provider);
	} else if (!strcmp(key, "conflict")) {
            ps = cast(char **) ptr;
            *ps = g_strdup(pkg.conflict);
        } else if (!strcmp(key, "data-requirement")) {
            pi = cast(int *) ptr;
            *pi = pkg.dreq;
        } else if (!strcmp(key, "model-requirement")) {
            pi = cast(int *) ptr;
            *pi = pkg.modelreq;
        } else if (!strcmp(key, "min-version")) {
            pi = cast(int *) ptr;
            *pi = pkg.minver;
        } else if (!strcmp(key, "lives-in-subdir")) {
            pi = cast(int *) ptr;
            *pi = pkg.uses_subdir;
        } else if (!strcmp(key, "wants-data-access")) {
            pi = cast(int *) ptr;
            *pi = pkg.data_access;
        } else if (!strcmp(key, "publist")) {
            plist = cast(int **) ptr;
            *plist = function_package_get_list(pkg, PkgListCode.publist, npub);
        } else if (!strcmp(key, "gui-publist")) {
            plist = cast(int **) ptr;
            *plist = function_package_get_list(pkg, PkgListCode.guilist, npub);
        } else if (!strcmp(key, "privlist")) {
            plist = cast(int **) ptr;
            *plist = function_package_get_list(pkg, PkgListCode.privlist, npriv);
        } else if (!strcmp(key, "assnlist")) {
            plist = cast(int **) ptr;
            *plist = function_package_get_list(pkg, PkgListCode.assnlist, npriv);
        } else if (!strcmp(key, "gui-main-id")) {
            pi = cast(int *) ptr;
            *pi = pkg_get_special_func_id(pkg, UfunRole.gui_main);
        } else if (!strcmp(key, "gui-attrs")) {
            ubyte *s = cast(ubyte *) ptr;

            pkg_get_gui_attrs(pkg, s);
        } else if ((role = pkg_key_get_role(key)) != 0) {
            ps = cast(char **) ptr;
            *ps = pkg_get_special_func_name(pkg, cast(UfunRole) role);
        }
    }

    va_end(ap);

    return err;
}

/* don't tamper with return value! */

extern (C) const(char) * function_package_get_string(Fnpkg * pkg, const(char) * id)
{
    if (pkg == null || id == null) {
        return null;
    } else if (!strcmp(id, "fname")) {
        return pkg.fname;
    } else if (!strcmp(id, "help-fname")) {
        return pkg.help_fname;
    } else if (!strcmp(id, "gui-help-fname")) {
        return pkg.gui_help_fname;
    } else if (!strcmp(id, "sample-fname")) {
        return pkg.sample_fname;
    } else if (!strcmp(id, "sample-script")) {
        return pkg.sample;
    } else if (!strcmp(id, "help")) {
        return pkg.help;
    } else if (!strcmp(id, "gui-help")) {
        return pkg.gui_help;
    } else if (!strcmp(id, "R-depends")) {
        return pkg.Rdeps;
    } else {
        return null;
    }
}

extern (C) char ** function_package_get_data_files(Fnpkg * pkg, int * n)
{
    char **S = null;

    *n = 0;
    if (pkg.datafiles != null) {
        S = strings_array_dup(pkg.datafiles, pkg.n_files);
        if (S != null) {
            *n = pkg.n_files;
        }
    }

    return S;
}

extern (C) char ** function_package_get_depends(Fnpkg * pkg, int * n)
{
    char **S = null;

    *n = 0;
    if (pkg.depends != null) {
        S = strings_array_dup(pkg.depends, pkg.n_depends);
        if (S != null) {
            *n = pkg.n_depends;
        }
    }

    return S;
}

extern (C) int function_package_set_data_files(Fnpkg * pkg, char ** S, int n)
{
    int err = 0;

    if (pkg.datafiles != null) {
        strings_array_free(pkg.datafiles, pkg.n_files);
        pkg.datafiles = null;
        pkg.n_files = 0;
    }

    if (n > 0) {
        if (S == null) {
            err = GretlError.data;
        } else {
            pkg.datafiles = strings_array_dup(S, n);
            if (pkg.datafiles == null) {
                err = GretlError.alloc;
            } else {
                pkg.n_files = n;
                pkg.uses_subdir = 1;
            }
        }
    }

    return err;
}

extern (C) int function_package_set_depends(Fnpkg * pkg, char ** S, int n)
{
    int err = 0;

    if (pkg.depends != null) {
        strings_array_free(pkg.depends, pkg.n_depends);
        pkg.depends = null;
        pkg.n_depends = 0;
    }

    if (n > 0) {
        if (S == null) {
            err = GretlError.data;
        } else {
            pkg.depends = strings_array_dup(S, n);
            if (pkg.depends == null) {
                err = GretlError.alloc;
            } else {
                pkg.n_depends = n;
            }
        }
    }

    if (!err && pkg.provider != null) {
        err = strings_array_prepend_uniq(&pkg.depends,
                                         &pkg.n_depends,
                                         pkg.provider);
    }

    return err;
}

/* For use by the GUI gfn editor */

extern (C) int function_package_set_must_assign(Fnpkg * pkg, char ** S, int n)
{
    Ufunc * f;
    int i;
    int pos;

    for (i=0; i<pkg.n_pub; i++) {
        f = pkg.pub[i];
        pos = strings_array_position(S, n, f.name);
        if (pos >= 0) {
            f.flags |= UfunAttrs.assign;
        } else {
            f.flags &= ~UfunAttrs.assign;
        }
    }

    return 0;
}

/* quick check to see if there's a gross problem with a package,
   in the context of considering packing it into a gretl
   session file
*/

int validate_function_package(Fnpkg * pkg)
{
    if (pkg.pub == null || pkg.author == null ||
        pkg.version_ == null || pkg.date == null ||
        pkg.descrip == null) {
        return 0;
    } else if (pkg.name[0] == '\0') {
        return 0;
    }

    return 1;
}

/* For session file use, and also MPI: dump currently defined
   functions into a single XML file. Gretlopt.OPT_U can be passed to
   limit the functions to ones that are _not_ members of a
   function package. Gretlopt.OPT_M is passed to indicate usage under
   MPI.
*/

extern (C) int write_loaded_functions_file(const(char) * fname, Gretlopt opt)
{
    Prn * prn;
    int mpicall;
    int unpackaged;
    int i;
    int nf;
    int err = 0;

    unpackaged = (opt & Gretlopt.OPT_U)? 1 : 0;

    if (n_ufuns == 0) {
        return 0;
    }

    if (unpackaged && !has_suffix(fname, ".xml")) {
        gretl_errmsg_set("Storing functions: the filename suffix must be '.xml'");
        return GretlError.invarg;
    }

    prn = gretl_print_new_with_filename(fname, &err);
    if (prn == null) {
        return err;
    }

    mpicall = (opt & Gretlopt.OPT_M)? 1 : 0;

    gretl_xml_header(prn);
    pputs(prn, "<gretl-functions>\n");


    if (!unpackaged) {
        /* write any loaded function packages */
        for (i=0; i<n_pkgs; i++) {
            if (validate_function_package(pkgs[i])) {
                real_write_function_package(pkgs[i], prn, mpicall);
            }
        }
    }

    /* then any unpackaged functions */
    nf = 0;
    for (i=0; i<n_ufuns; i++) {
        if (ufuns[i].pkg == null) {
            write_function_xml(ufuns[i], prn, mpicall);
            nf++;
        }
    }
    if (unpackaged && nf == 0) {
        gretl_errmsg_set("No unpackaged functiond were found");
        err = GretlError.data;
    }

    pputs(prn, "</gretl-functions>\n");

    gretl_print_destroy(prn);

    if (err) {
        gretl_remove(fname);
    }

    return err;
}

/* De-allocate a function package: this can be done in either of two
   modes.  If 'full' is non-zero then we destroy all functions that
   are children of the given package, otherwise we leave the
   function-children alone (just 'detaching' them from the parent
   package).
*/

void real_function_package_free(Fnpkg * pkg, int full)
{
    if (pkg != null) {
        int i;

        if (full) {
            for (i=0; i<pkg.n_pub; i++) {
                ufunc_free(pkg.pub[i]);
            }
            for (i=0; i<pkg.n_priv; i++) {
                ufunc_free(pkg.priv[i]);
            }
        } else {
            for (i=0; i<n_ufuns; i++) {
                if (ufuns[i].pkg == pkg) {
                    /* remove package info */
                    ufuns[i].pkg = null;
                    set_function_private(ufuns[i], 0);
                }
            }
        }

        if (pkg.datafiles != null && pkg.n_files > 0) {
            strings_array_free(pkg.datafiles, pkg.n_files);
        }
        if (pkg.depends != null && pkg.n_depends > 0) {
            strings_array_free(pkg.depends, pkg.n_depends);
        }

        free(pkg.pub);
        free(pkg.priv);
        free(pkg.fname);
        free(pkg.author);
        free(pkg.email);
        free(pkg.version_);
        free(pkg.date);
        free(pkg.descrip);
        free(pkg.help);
        free(pkg.gui_help);
        free(pkg.Rdeps);
        free(pkg.sample);
        free(pkg.help_fname);
        free(pkg.gui_help_fname);
        free(pkg.sample_fname);
        free(pkg.tags);
        free(pkg.label);
        free(pkg.mpath);
        free(pkg.provider);
	free(pkg.conflict);
        free(pkg);
    }
}

void function_package_free(Fnpkg * pkg)
{
    real_function_package_free(pkg, 0);
}

void function_package_free_full(Fnpkg * pkg)
{
    real_function_package_free(pkg, 1);
}

/* is the package with filename @fname already in memory? */

Fnpkg * get_loaded_pkg_by_filename(const(char) * fname, const(char) ** version_)
{
    int i;

    if (fname == null) {
        return null;
    }

    for (i=0; i<n_pkgs; i++) {
        if (!strcmp(fname, pkgs[i].fname)) {
            if (version_ != null) {
                *version_ = pkgs[i].version_;
            }
            return pkgs[i];
        }
    }

    return null;
}

Fnpkg * get_loaded_pkg_by_name(const(char) * pkgname)
{
    int i;

    if (pkgname == null) {
        return null;
    }

    for (i=0; i<n_pkgs; i++) {
        if (!strcmp(pkgname, pkgs[i].name)) {
            return pkgs[i];
        }
    }

    return null;
}

/**
 * function_package_unload_by_filename:
 * @fname: package filename.
 *
 * Unloads the specified function package from memory, if it
 * is currently loaded.  The functions 'owned' by the package
 * are not destroyed; they become available for inclusion in
 * other packages.
 */

extern (C) void function_package_unload_by_filename(const(char) * fname)
{
    Fnpkg *pkg = get_loaded_pkg_by_filename(fname, null);

    if (pkg != null) {
        real_function_package_unload(pkg, 0);
    }
}

/**
 * function_package_unload_full_by_filename:
 * @fname: package filename.
 *
 * Unloads the specified function package from memory, if it
 * is currently loaded.  The functions 'owned' by the package
 * are also unloaded from memory.
 */

extern (C) void function_package_unload_full_by_filename(const(char) * fname)
{
    Fnpkg *pkg = get_loaded_pkg_by_filename(fname, null);

    if (pkg != null) {
        real_function_package_unload(pkg, 1);
    }
}

/**
 * function_package_unload_full:
 * @pkgname: package name.
 *
 * Unloads the specified function package from memory, if it
 * is currently loaded.  The functions 'owned' by the package
 * are also unloaded from memory.
 *
 * Returns: 1 if the specified package was loaded and is
 * now unloaded, otherwise 0.0
 */

extern (C) int function_package_unload_full(const(char) * pkgname)
{
    Fnpkg *pkg = get_loaded_pkg_by_name(pkgname);

    if (pkg != null) {
        real_function_package_unload(pkg, 1);
        return 1;
    } else {
        return 0;
    }
}

/* append a function package to the recorder array of loaded
   packages */

int function_package_record(Fnpkg * pkg)
{
    Fnpkg **tmp;
    int err = 0;

    tmp = cast(Fnpkg* *) realloc(pkgs, (n_pkgs + 1) * (*tmp).sizeof);

    if (tmp == null) {
        err = GretlError.alloc;
    } else {
        pkgs = tmp;
        pkgs[n_pkgs] = pkg;
        n_pkgs++;
    }

    return err;
}

int broken_package_error(Fnpkg * pkg)
{
    gretl_errmsg_sprintf(_("'%s': package contains " ~
		"duplicated function names").ptr,
                         pkg.name);
    return GretlError.data;
}


void do_override(Ufunc * fun)
{
    gretl_warnmsg_sprintf(_("'%s' is the name of a built-in function").ptr,
                          fun.name);
    install_function_override(fun.name, fun.pkg.name, fun);
    fun.pkg.overrides += 1;
}

/* When loading a private function the only real conflict would be
   with a function of the same name owned by the same package.
   Obviously this shouldn't happen but we'll whack it if it does.
*/

int load_private_function(Fnpkg * pkg, int i, Prn * prn)
{
    Ufunc * fun = pkg.priv[i];
    int j;
    int err = 0;

    for (j=0; j<n_ufuns; j++) {
        if (!strcmp(fun.name, ufuns[j].name)) {
            if (pkg == ufuns[j].pkg) {
                return broken_package_error(pkg);
            }
        }
    }


    if (!err) {
        err = add_allocated_ufunc(fun);
        if (!err && function_lookup(fun.name)) {
            do_override(fun);
        }
    }

    return err;
}

/* When loading a public, packaged function we want to avoid conflicts
   with any non-private functions of the same name.  In case we get a
   conflict with a public member of an already loaded distinct package
   we'll flag an error (otherwise things will get too confusing).
*/

int load_public_function(Fnpkg * pkg, int i)
{
    Ufunc * fun = pkg.pub[i];
    int j;
    int done = 0;
    int err = 0;

    for (j=0; j<n_ufuns; j++) {
        if (!strcmp(fun.name, ufuns[j].name)) {
            if (pkg == ufuns[j].pkg) {
                /* name duplication in single package */
                err = broken_package_error(pkg);
            } else if (ufuns[j].pkg == null) {
                /* conflicting unpackaged function */
                ufunc_free(ufuns[j]);
                ufuns[j] = fun;
                done = 1;
            } else if (!function_is_private(ufuns[j])) {
                /* got a conflicting package */
                gretl_errmsg_sprintf(_("The function %s is already defined " ~
		"by package '%s'").ptr, fun.name,
                                     ufuns[j].pkg.name);
                err = GretlError.data;
                break;
            }
        }
    }

    if (!err && !done && function_lookup(fun.name)) {
        if (!strcmp(pkg.name, "PanelTools".ptr) &&
            !strcmp(fun.name, "pxmean".ptr)) {
            /* give this a pass: just use the built-in pxmean() */
            return 0;
        } else if (strcmp(fun.name, "bkw".ptr)) {
            /* for now, don't throw an error on loading Lee Adkins'
               bkw package */
            gretl_errmsg_sprintf(_("Package loading error:\n" ~
		"%s %s defines a function %s(), " ~
		"which conflicts with a built-in\n" ~
		"function in gretl %s. Please check for an " ~
		"updated version of %s.").ptr,
                                 pkg.name, pkg.version_, fun.name,
                                 GRETL_VERSION.ptr, pkg.name);
            err = GretlError.data;
        }
    }

    if (!err && !done) {
        err = add_allocated_ufunc(fun);
    }

    return err;
}

int pkg_in_stack(const(char) * name, GArray * pstack)
{
    char * s;
    int i;

    for (i=0; i<pstack.len; i++) {
        s = (cast(char **) pstack.data)[i];
        if (!strcmp(name, s)) {
            return 1;
        }
    }

    return 0;
}

bool TS_REQUIRED(DataReq d)
{
    return d == DataReq.needs_ts || d == DataReq.needs_qm;
}

int data_requirements_conflict(Fnpkg * a, Fnpkg * b)
{
    if (a.dreq == DataReq.needs_panel && TS_REQUIRED(b.dreq)) {
        return 1;
    } else if (b.dreq == DataReq.needs_panel && TS_REQUIRED(a.dreq)) {
        return 1;
    }
    return 0;
}
