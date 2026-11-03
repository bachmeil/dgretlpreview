/* Port of gretl_func.c, part 3 of 4: see gretl_func0.d. */

module gretl.gretl_func3;

import gretl.garray : g_array_append_vals;
import gretl.gfileutils : g_path_get_basename;


import gretl.gretl_intl : gretl_pop_c_numeric_locale, gretl_push_c_numeric_locale;
import gretl.gretl_func_ext : check_R_depends, get_command_index, is_gretl_addon, md_to_gretl, simple_markdown_detect, statements_get_structure;
import gretl.xmlh : xmlFreeDoc, xmlStrcmp;

import gretl.gfileutils : GFileTest;
import gretl.gretl_pathsh : PkgType;
import gretl.gretl_prnh : PrnType;
import gretl.strutilsh : SLASHSTR;

import gretl.cmd_privateh : CmdFlags;
import gretl.glist : g_list_next;
import gretl.gretl_func1 : FuncsTask;
import gretl.gretl_paths : gretl_opendir;
import gretl.gretl_typemap : gretl_type_get_name;
import gretl.strutils : dot_atof, gretl_strstrip;

import gretl.cmd_privateh : Cmd;

import gretl.cmd_privateh : Cmd;

import gretl.cmd_privateh : Cmd;

import gretl.cmd_privateh : Cmd;
import gretl.dataset : series_set_stack_level;
import gretl.garray : g_array_free, g_array_new;
import gretl.genmain : gretl_is_series, gretl_reserved_word;
import gretl.gfileutils : GDir, g_dir_close, g_dir_read_name, g_dir_rewind, g_file_test, g_get_current_dir, g_path_is_absolute;
import gretl.glib : gchar;
import gretl.glist : GList, g_list_first, g_list_free, g_list_prepend;
import gretl.gmem : g_free;
import gretl.gretl_array : gretl_array_from_strings;
import gretl.gretl_bundle : gretl_bundle_donate_data, gretl_bundle_get_lists, gretl_bundle_set_int, gretl_bundle_set_scalar, gretl_bundle_set_string;
import gretl.gretl_commands : gretl_command_number;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_errors : gretl_errmsg_prepend, gretl_errmsg_set, gretl_errmsg_sprintf, gretl_error_clear;
import gretl.gretl_func1 : add_ufunc, allocate_params, arg_may_be_optional, changes_state, clear_ufunc_data, current_function_call, free_params_array, function_in_use, function_package_alloc, get_user_function_by_name, gretl_function_depth, gretl_function_print_code, maybe_set_param_const, nullify_upname, push_function_line, read_ufunc_from_xml, real_function_package_unload, return_type_from_string, set_compiling_off, set_compiling_on, ufunc_unload;
import gretl.gretl_func2 : data_requirements_conflict, function_package_free, function_package_free_full, function_package_record, get_loaded_pkg_by_filename, get_loaded_pkg_by_name, is_pdf_ref, load_private_function, load_public_function, pkg_in_stack;
import gretl.gretl_list : gretl_list_append_term, in_gretl_list;
import gretl.gretl_paths : get_full_filename, gretl_build_path, gretl_chdir, gretl_fopen, gretl_function_package_get_path;
import gretl.gretl_prn : gretl_print_destroy, gretl_print_get_buffer, gretl_print_new, gretl_print_steal_buffer, pprintf, pputc, pputs;
import gretl.gretl_typemap : gretl_is_array_type, gretl_type_from_string, gretl_type_get_plain_type;
import gretl.gretl_utils : gretl_version_number, gretl_version_string;
import gretl.gretl_xml : gretl_xml_get_gfn, gretl_xml_get_prop_as_bool, gretl_xml_get_prop_as_string, gretl_xml_get_string, gretl_xml_get_strings_array, gretl_xml_node_get_double, gretl_xml_node_get_trimmed_string, gretl_xml_open_doc_root;
import gretl.gstrfuncs : g_strdup, g_strdup_printf, g_strndup;
import gretl.libset : gretl_messages_on;
import gretl.missingh : NADBL, na;
import gretl.printout : bufgets, bufgets_finalize, bufgets_init;
import gretl.strutils : gretl_namechar_spn, gretl_strdup, gretl_string_split_quoted, gretl_strndup, has_suffix, string_is_blank, strings_array_free, strings_array_new;
import gretl.strutilsh : SLASH;
import gretl.uservar : arg_add_as_shell, copy_as_arg, copy_list_as_arg, get_list_by_name, gretl_is_user_var;

import gretl.cmd_privateh : ExecState;
import gretl.garray : GArray;
import gretl.xmlh : xmlDocPtr, xmlNodePtr;

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
import gretl.genlex : function_lookup;
import gretl.gretl_func0;
import gretl.xmlh : xmlChar;
import gretl.gretl_funch;
import gretl.libgretlh;
import gretl.gretl_matrixh;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_prn : Prn;

private alias XUC = const(xmlChar) *;

int load_gfn_dependencies(Fnpkg * pkg, GArray * pstack, int level)
{
    int err = 0;

    if (pkg.depends != null) {
        char * pkgpath;
        int i;

        fprintf(stderr, "*** load_gfn_dependencies (%d) for %s (level %d) ***\n",
                pkg.n_depends, pkg.name, level);

        for (i=0; i<pkg.n_depends && !err; i++) {
            const(char) * depname = pkg.depends[i];
            Fnpkg * dep = null;

            if (get_function_package_by_name(depname) != null) {
                fprintf(stderr, " %s: already loaded\n", depname);
                ; /* OK, already loaded */
            } else if (pkg_in_stack(depname, pstack)) {
                fprintf(stderr, " %s: found in stack\n", depname);
                ; /* don't go into infinite loop! */
            } else {
                fprintf(stderr, " trying for %s\n", depname);
                pkgpath = gretl_function_package_get_path(depname, PkgType.all);
                if (pkgpath == null) {
                    err = GretlError.data;
                    gretl_errmsg_sprintf(_("%s: dependency %s was not found").ptr,
                                         pkg.name, depname);
                } else {
                    err = load_function_package(pkgpath, Gretlopt.OPT_NONE,
                                                pstack, null, &dep,
                                                level + 1);
                    free(pkgpath);
                    if (!err && level == 0) {
                        err = data_requirements_conflict(pkg, dep);
                    }
                    if (!err) {
                        {
        auto tmp_ = depname;

        g_array_append_vals(pstack, &tmp_, 1);
    }
                    }
                }
            }
        }
    }

    return err;
}

/* A 'real load' is in contrast to just reading some info from a
   package, as is done in various GUI contexts.  We do the real load
   in response to some GUI commands, the "include" command, and also
   when re-opening a gretl session file that contains function
   packages.
*/

int real_load_package(Fnpkg * pkg, GArray * pstack, Prn * prn, int level)
{
    int i;
    int err = 0;


    gretl_error_clear();

    if (pstack != null) {
        err = load_gfn_dependencies(pkg, pstack, level);
    }

    if (!err && pkg.pub != null) {
        for (i=0; i<pkg.n_pub && !err; i++) {
            err = load_public_function(pkg, i);
        }
    }

    if (!err && pkg.priv != null) {
        for (i=0; i<pkg.n_priv && !err; i++) {
            err = load_private_function(pkg, i, prn);
        }
    }

    if (!err && pkg.provider != null) {
        /* check that provider really got loaded */
        if (get_function_package_by_name(pkg.provider) == null) {
            gretl_errmsg_sprintf(_("Provider package %s is not loaded\n").ptr,
                                 pkg.provider);
            err = GretlError.data;
        }
    }

    if (!err) {
        /* add to array of loaded packages */
        err = function_package_record(pkg);
    }

    return err;
}

const(char) * data_needs_string(DataReq dr)
{
    if (dr == DataReq.needs_ts) {
        return "Time-series data".ptr;
    } else if (dr == DataReq.needs_qm) {
        return "Quarterly or monthly data".ptr;
    } else if (dr == DataReq.needs_panel) {
        return "Panel data".ptr;
    } else if (dr == DataReq.nodata_ok) {
        return "none".ptr;
    } else {
        return "some sort of dataset".ptr;
    }
}

void pkg_print_depends(Fnpkg * pkg, Prn * prn)
{
    int rdep = pkg.Rdeps != null;
    int i;
    int n = pkg.n_depends + rdep;

    pputs(prn, "<@itl=\"Dependencies\">: ");
    for (i=0; i<pkg.n_depends; i++) {
        pprintf(prn, "%s%s", pkg.depends[i],
                (i < n-1)? ", ".ptr : ".\n".ptr);
    }
    if (rdep) {
        pprintf(prn, "%s%s\n".ptr, pkg.n_depends ? ",".ptr : ".".ptr,
                pkg.Rdeps);
    }
}

int help_text_is_markdown(Fnpkg * pkg, int gui_help)
{
    const(char) * fname;
    const(char) * text;

    text = gui_help ? pkg.gui_help : pkg.help;
    if (text == null) {
        return 0;
    }

    fname = gui_help ? pkg.gui_help_fname : pkg.help_fname;
    if (fname != null && has_suffix(fname, ".md")) {
        /* A filename matching "*.md" should clinch it. */
        return 1;
    }

    if (pkg.date != null) {
        /* Use of markdown in gfn help text was first supported in
           gretl 2023b (2023-07-21); we assume that packages with
           an earlier last-modified date can't be using it.
        */
        int y;
        int m;
        int d;

        if (sscanf(pkg.date, "%4d-%02d-%02d", &y, &m, &d) == 3) {
            if (10000*y + 100*m + d < 20230721) {
                return 0;
            }
        }
    }

    /* We might have anonymous help text (no filename specified),
       in which case we fall back on a simple heuristic.
    */
    return simple_markdown_detect(text);
}

void print_package_info(Fnpkg * pkg, const(char) * fname, Prn * prn)
{
    char[8] vstr;
    int remote;
    int pdfdoc;
    int md = 0;
    int i;

    remote = (strstr(fname, "dltmp.") != null);

    if (!remote && g_path_is_absolute(fname)) {
        pprintf(prn, "<@itl=\"File\">: %s\n\n", fname);
    }

    if (pkg.name[0] == '\0' || pkg.author == null ||
        pkg.minver <= 0 || pkg.descrip == null ||
        pkg.version_ == null || pkg.date == null ||
        pkg.help == null) {
        pprintf(prn, _("\nBroken package! Basic information is missing\n").ptr);
        return;
    }

    gretl_version_string(vstr.ptr, pkg.minver);
    pdfdoc = is_pdf_ref(pkg.help);

    pprintf(prn, "<@itl=\"Package\">: %s %s (%s)\n", pkg.name, pkg.version_,
            pkg.date);
    pprintf(prn, "<@itl=\"Author\">: %s\n", pkg.author);
    if (pkg.email != null && *pkg.email != '\0') {
        pprintf(prn, "<@itl=\"Email\">: %s\n", pkg.email);
    }
    pprintf(prn, "<@itl=\"Required gretl version\">: %s\n", vstr.ptr);
    pprintf(prn, "<@itl=\"Data requirement\">: %s\n", _(data_needs_string(pkg.dreq)));
    pprintf(prn, "<@itl=\"Description\">: %s\n", gretl_strstrip(pkg.descrip));
    if (pkg.n_depends > 0 || pkg.Rdeps != null) {
        pkg_print_depends(pkg, prn);
    }
    if (pkg.provider != null) {
        pprintf(prn, "<@itl=\"Provider\">: %s\n", pkg.provider);
    }
    if (pkg.conflict != null) {
        pprintf(prn, "<@itl=\"Conflict\">: %s\n", pkg.conflict);
    }

    if (pdfdoc) {
        const(char) * s = strrchr(pkg.help, ':');

        if (remote) {
            pprintf(prn, "<@itl=\"Documentation\">: %s\n\n", s + 1);
        } else {
            gchar *localpdf = g_strdup(fname);
            gchar *p = strrchr(localpdf, '.');

            p[0] = '\0';
            strcat(p, ".pdf");
            pprintf(prn, "<@itl=\"Documentation\">: <@adb=\"%s\">\n\n", localpdf);
            g_free(localpdf);
        }
    } else {
        pputc(prn, '\n');
    }

    if (pkg.pub != null) {
	/* note: pkg.pub is non-null only if the package is loaded */
        if (pkg.n_pub == 1) {
            if (strcmp(pkg.pub[0].name, pkg.name)) {
                pputs(prn, "<@itl=\"Public interface\">: ");
                pputs(prn, "<mono>\n");
                pprintf(prn, "%s()\n", pkg.pub[0].name);
                pputs(prn, "</mono>\n\n");
            }
        } else {
            pputs(prn, "<@itl=\"Public interfaces\">:\n\n");
            pputs(prn, "<mono>\n");
            for (i=0; i<pkg.n_pub; i++) {
                pprintf(prn, "  %s()\n", pkg.pub[i].name);
            }
            pputs(prn, "</mono>\n\n");
        }
    }

    md = help_text_is_markdown(pkg, 0);

    if (!pdfdoc) {
        if (md) {
            md_to_gretl(pkg.help, prn);
        } else {
            pputs(prn, "<@itl=\"Help text\">:\n\n");
            pputs(prn, "<mono>\n");
            pputs(prn, pkg.help);
            pputs(prn, "\n\n");
            pputs(prn, "</mono>\n");
        }
    }

    if (remote && pkg.sample != null) {
        pputs(prn, "<@itl=\"Sample script\">:\n\n");
        pputs(prn, "<code>\n");
        pputs(prn, pkg.sample);
        pputs(prn, "</code>\n");
        pputc(prn, '\n');
    }
}

/* start apparatus for augmenting "pkg query" to cover contents
   of the "examples" directory of a function package, as of
   February/March 2024.0
*/

int resource_wanted(const(char) * s)
{
    return has_suffix(s, ".inp") ||
	has_suffix(s, ".gdt") ||
	has_suffix(s, ".gdtb");
}

int pkg_has_examples(Fnpkg * pkg)
{
    int i;

    for (i=0; i<pkg.n_files; i++) {
        if (!strcmp(pkg.datafiles[i], "examples")) {
            return 1;
        }
    }

    return 0;
}

int pkg_has_datafiles(Fnpkg * pkg)
{
    int i;

    for (i=0; i<pkg.n_files; i++) {
        if (strlen(pkg.datafiles[i]) > 4 &&
	    has_suffix(pkg.datafiles[i], ".gdt")) {
            return 1;
        }
    }

    return 0;
}

/* helper for compare_resources() below */

int same_subdir(const(char) * sa, const(char) * sb)
{
    int n = cast(int) strcspn(sa, SLASHSTR);

    if (strcspn(sb, SLASHSTR) != n) {
	return 0;
    } else {
	return strncmp(sa, sb, n) == 0;
    }
}

/* qsort callback for resource filenames */

extern (C) int compare_resources(const(void) * a, const(void) * b)
{
    const(char) * sa = *cast(const(char) **) a;
    const(char) * sb = *cast(const(char) **) b;
    int ac = strchr(sa, SLASH) != null;
    int bc = strchr(sb, SLASH) != null;

    if (ac + bc == 1) {
	/* sort top-level files before subdir contents */
	return ac ? 1 : -1;
    } else {
	/* sort scripts before data files, per directory */
	int same_dir = (ac + bc == 0) ? 1 : same_subdir(sa, sb);

	if (same_dir) {
	    ac = has_suffix(sa, ".inp");
	    bc = has_suffix(sb, ".inp");
	    if (ac + bc == 1) {
		return ac ? -1 : 1;
	    }
	}
	/* just alphabetize */
	return strcmp(sa, sb);
    }
}

char ** package_resources_info(Fnpkg * pkg, const(char) * fname, char * resdir, int * pn_res)
{
    char **resources = null;
    char * pkg_base = null;
    const(char) * dname;
    const(char) * subdir;
    char[128] tmp;
    char * savedir = null;
    GDir *dir = null;
    GDir *child = null;
    int path_len = 0;
    int n_res = 0;

    path_len = cast(int) (strlen(fname) - strlen(pkg.name) - 4);
    pkg_base = g_strndup(fname, path_len);

    gretl_build_path(resdir, pkg_base, "examples".ptr, null);
    dir = gretl_opendir(resdir);
    if (dir == null) {
	g_free(pkg_base);
	return null;
    }

    /* cd into the examples directory to make life easier */
    savedir = g_get_current_dir();
    gretl_chdir(resdir);

    /* first get a count of relevant files */
    while ((dname = g_dir_read_name(dir)) != null) {
	if (g_file_test(dname, GFileTest.is_dir)) {
	    child = gretl_opendir(dname);
	    if (child != null) {
		while ((dname = g_dir_read_name(child)) != null) {
		    if (resource_wanted(dname)) {
			n_res++;
		    }
		}
		g_dir_close(child);
	    }
	} else if (resource_wanted(dname)) {
	    n_res++;
	}
    }

    if (n_res > 0) {
	/* if there are any relevant files, create an array */
	int i = 0;

	resources = strings_array_new(n_res);
	g_dir_rewind(dir);

	while ((dname = g_dir_read_name(dir)) != null) {
	    if (g_file_test(dname, GFileTest.is_dir)) {
		child = gretl_opendir(dname);
		if (child != null) {
		    subdir = dname;
		    while ((dname = g_dir_read_name(child)) != null) {
			if (resource_wanted(dname)) {
			    sprintf(tmp.ptr, "%s%c%s", subdir, SLASH, dname);
			    resources[i++] = gretl_strdup(tmp.ptr);
			}
		    }
		    g_dir_close(child);
		}
	    } else if (resource_wanted(dname)) {
		resources[i++] = gretl_strdup(dname);
	    }
	}
    }

    if (dir != null) {
        g_dir_close(dir);
    }

    if (savedir != null) {
	/* return to the original working directory */
	gretl_chdir(savedir);
	g_free(savedir);
    }

    *pn_res = n_res;
    g_free(pkg_base);

    if (resources != null) {
	qsort(resources, n_res, (*resources).sizeof, &compare_resources);
    }

    return resources;
}

char ** package_datafiles_info(Fnpkg * pkg, const(char) * fname, char * resdir, int * pn_res)
{
    char **resources = null;
    const(char) * dname;
    GDir *dir = null;
    int path_len = 0;
    int n_res = 0;

    path_len = cast(int) (strlen(fname) - strlen(pkg.name) - 5);
    resdir[0] = '\0';
    strncat(resdir, fname, path_len);

    dir = gretl_opendir(resdir);
    if (dir == null) {
	return null;
    }

    /* first get a count of relevant files */
    while ((dname = g_dir_read_name(dir)) != null) {
	if (resource_wanted(dname)) {
	    n_res++;
	}
    }

    if (n_res > 0) {
	/* if there are any relevant files, create an array */
	int i = 0;

	resources = strings_array_new(n_res);
	g_dir_rewind(dir);
	while ((dname = g_dir_read_name(dir)) != null) {
	    if (resource_wanted(dname)) {
		resources[i++] = gretl_strdup(dname);
	    }
	}
    }

    if (dir != null) {
        g_dir_close(dir);
    }

    *pn_res = n_res;

    if (resources != null) {
	qsort(resources, n_res, (*resources).sizeof, &compare_resources);
    }

    return resources;
}

/* end apparatus for augmenting "pkg query" to "examples" */

void plain_print_package_info(Fnpkg * pkg, const(char) * fname, Prn * prn)
{
    char[8] vstr;
    int has_examples = 0;
    int i;

    if (g_path_is_absolute(fname) && !strstr(fname, "dltmp.")) {
        pprintf(prn, "File: %s\n", fname);
    }
    if (pkg.name[0] == '\0' || pkg.author == null ||
        pkg.minver <= 0 || pkg.descrip == null ||
        pkg.version_ == null || pkg.date == null ||
        pkg.help == null) {
        pprintf(prn, _("Broken package! Basic information is missing\n").ptr);
        return;
    }

    gretl_version_string(vstr.ptr, pkg.minver);

    pprintf(prn, "Package: %s %s (%s)\n", pkg.name, pkg.version_, pkg.date);
    pprintf(prn, "Author: %s\n", pkg.author);
    if (pkg.email != null && *pkg.email != '\0') {
        pprintf(prn, "Email: %s\n", pkg.email);
    }
    pprintf(prn, "Required gretl version: %s\n", vstr.ptr);
    pprintf(prn, "Data requirement: %s\n", _(data_needs_string(pkg.dreq)));
    pprintf(prn, "Description: %s\n", gretl_strstrip(pkg.descrip));
    if (pkg.n_depends > 0) {
        pputs(prn, "Dependencies: ");
        for (i=0; i<pkg.n_depends; i++) {
            pprintf(prn, "%s%s", pkg.depends[i],
                    (i < pkg.n_depends-1)? ", ".ptr : ".\n".ptr);
        }
    }
    if (pkg.provider != null) {
        pprintf(prn, "Provider: %s\n", pkg.provider);
    }
    if (pkg.conflict != null) {
        pprintf(prn, "Conflict: %s\n", pkg.conflict);
    }
    if (pkg.mpath != null) {
	pprintf(prn, "Menu attachment: %s\n", pkg.mpath);
    }

    has_examples = pkg_has_examples(pkg);

    if (has_examples || pkg_has_datafiles(pkg)) {
        char[FILENAME_MAX] resdir;
        char **resources;
        int n_res = 0;

	if (has_examples) {
	    resources = package_resources_info(pkg, fname, resdir.ptr, &n_res);
	} else {
	    resources = package_datafiles_info(pkg, fname, resdir.ptr, &n_res);
	}
        if (resources != null) {
            pprintf(prn, "Resources in %s:\n", resdir.ptr);
            for (i=0; i<n_res; i++) {
                pprintf(prn, "%s\n", resources[i]);
            }
            strings_array_free(resources, n_res);
        }
    }

    pputc(prn, '\n');
}

void real_bundle_package_info(Fnpkg * pkg, const(char) * fname, GretlBundle * b)
{
    int has_examples = 0;
    int err = 0;

    if (g_path_is_absolute(fname) && !strstr(fname, "dltmp.")) {
        gretl_bundle_set_string(b, "file", fname);
    }
    if (pkg.name[0] == '\0' || pkg.author == null ||
        pkg.minver <= 0 || pkg.descrip == null ||
        pkg.version_ == null || pkg.date == null ||
        pkg.help == null) {
        gretl_bundle_set_int(b, "broken", 1);
        return;
    }

    gretl_bundle_set_string(b, "name", pkg.name);
    gretl_bundle_set_string(b, "author", pkg.author);
    /* for backward compatibility, we need a numerical package
       version in the bundle returned by "pkg query"
    */
    if (is_gretl_addon(pkg.name)) {
	gretl_bundle_set_scalar(b, "version",
				gretl_version_number(pkg.version_));
    } else {
	gretl_bundle_set_scalar(b, "version", dot_atof(pkg.version_));
    }
    gretl_bundle_set_string(b, "date", pkg.date);
    if (pkg.email != null && *pkg.email != '\0') {
        gretl_bundle_set_string(b, "email", pkg.email);
    }
    gretl_bundle_set_int(b, "gretl_version", pkg.minver);
    gretl_bundle_set_string(b, "description", gretl_strstrip(pkg.descrip));
    gretl_bundle_set_int(b, "gretl_version", pkg.minver);
    gretl_bundle_set_string(b, "data_requirement", data_needs_string(pkg.dreq));

    if (pkg.n_depends > 0) {
        GretlArray *D = gretl_array_from_strings(pkg.depends,
                                                  pkg.n_depends,
                                                  1, &err);
        gretl_bundle_donate_data(b, "depends", D, GretlType.array);
    }
    if (pkg.provider != null) {
        gretl_bundle_set_string(b, "provider", pkg.provider);
    }
    if (pkg.conflict != null) {
        gretl_bundle_set_string(b, "conflict", pkg.conflict);
    }
    if (pkg.mpath != null) {
        gretl_bundle_set_string(b, "menu_path", pkg.mpath);
    }

    has_examples = pkg_has_examples(pkg);

    if (has_examples || pkg_has_datafiles(pkg)) {
        char[FILENAME_MAX] resdir;
        char **resources;
        int n_res = 0;

	if (has_examples) {
	    resources = package_resources_info(pkg, fname, resdir.ptr, &n_res);
	} else {
	    resources = package_datafiles_info(pkg, fname, resdir.ptr, &n_res);
	}
        if (resources != null) {
            GretlArray *R = gretl_array_from_strings(resources, n_res, 0, &err);

            if (!err) {
		gretl_bundle_set_string(b, "resource_dir", resdir.ptr);
                gretl_bundle_donate_data(b, "resources", R, GretlType.array);
            }
        }
    }
}

void reflow_package_help(const(char) * help, Prn * prn)
{
    char * rem;
    char * prev;
    char[2048] line;
    int i;
    int lmax = 80;

    bufgets_init(help);

    while (bufgets(line.ptr, line.sizeof, help)) {
        if (strlen(line.ptr) <= lmax) {
            pputs(prn, line.ptr);
        } else {
            rem = line.ptr;
            prev = null;
            while (strlen(rem) > lmax) {
                for (i=lmax-1; i>0; i--) {
                    if (rem[i] == ' ') {
                        rem[i] = '\0';
                        pprintf(prn, "%s\n", rem);
                        rem = rem + i + 1;
                        break;
                    }
                }
                if (rem == line.ptr || rem == prev) {
                    /* let's not get into an infinite loop */
                    break;
                }
                prev = rem;
            }
            if (*rem != '\0') {
                pputs(prn, rem);
            }
        }
    }

    bufgets_finalize(help);
    pputs(prn, "\n\n");
}

/* Simple plain-text help output, or converted markdown if applicable
   -- that is, we're in GUI mode and @pkg has help text in
   markdown. Being in GUI mode is signalled by @pbuf being non-null.
*/

void print_package_help(Fnpkg * pkg, const(char) * fname, char ** pbuf, Prn * prn)
{
    pprintf(prn, "%s %s (%s), %s\n", pkg.name, pkg.version_,
            pkg.date, pkg.author);
    pputs(prn, gretl_strstrip(pkg.descrip));
    pputs(prn, "\n\n");

    if (pbuf != null && help_text_is_markdown(pkg, 0)) {
        Prn *myprn = gretl_print_new(PrnType.buffer, null);

         /* FIXME should the next line really be needed? */
        pprintf(myprn, "<@bld=\"%s\">\n\n", pkg.name);
        md_to_gretl(pkg.help, myprn);
        *pbuf = gretl_print_steal_buffer(myprn);
        gretl_print_destroy(myprn);
    } else {
         reflow_package_help(pkg.help, prn);
    }
}

void print_package_code(Fnpkg * pkg, int tabwidth, Prn * prn)
{
    int i;

    pprintf(prn, "# hansl code from package %s %s (%s)\n\n",
            pkg.name, pkg.version_, pkg.date);

    if (pkg.priv != null) {
        pputs(prn, "# private functions\n\n");
        for (i=0; i<pkg.n_priv; i++) {
            gretl_function_print_code(pkg.priv[i], tabwidth, prn);
            pputc(prn, '\n');
        }
    }

    if (pkg.pub != null) {
        pputs(prn, "# public functions\n\n");
        for (i=0; i<pkg.n_pub; i++) {
            gretl_function_print_code(pkg.pub[i], tabwidth, prn);
            pputc(prn, '\n');
        }
    }
}

/* allocate a Fnpkg structure and read from XML file into it */

Fnpkg * real_read_package(xmlDocPtr doc, xmlNodePtr node, const(char) * fname, int get_funcs, int * err)
{
    xmlNodePtr cur;
    Fnpkg * pkg;
    char * tmp = null;


    pkg = function_package_alloc(fname);
    if (pkg == null) {
        *err = GretlError.alloc;
        return null;
    }

    gretl_xml_get_prop_as_string(node, "name", &tmp);

    if (tmp == null) {
        fprintf(stderr, "real_read_package: package has no name\n");
        *err = GretlError.data;
        function_package_free(pkg);
        return null;
    }

    strncat(pkg.name, tmp, FN_NAMELEN - 1);
    free(tmp);

    if (gretl_xml_get_prop_as_bool(node, NEEDS_TS)) {
        pkg.dreq = DataReq.needs_ts;
    } else if (gretl_xml_get_prop_as_bool(node, NEEDS_QM)) {
        pkg.dreq = DataReq.needs_qm;
    } else if (gretl_xml_get_prop_as_bool(node, NEEDS_PANEL)) {
        pkg.dreq = DataReq.needs_panel;
    } else if (gretl_xml_get_prop_as_bool(node, NO_DATA_OK)) {
        pkg.dreq = DataReq.nodata_ok;
    }

    if (gretl_xml_get_prop_as_string(node, "model-requirement", &tmp)) {
        pkg.modelreq = gretl_command_number(tmp);
        free(tmp);
    }

    if (gretl_xml_get_prop_as_string(node, "minver", &tmp)) {
        pkg.minver = gretl_version_number(tmp);
        free(tmp);
    }

    pkg.uses_subdir = cast(ubyte) gretl_xml_get_prop_as_bool(node, "lives-in-subdir");
    pkg.data_access = cast(ubyte) gretl_xml_get_prop_as_bool(node, "wants-data-access");

    cur = node.children;

    while (cur != null) {
        if (!xmlStrcmp(cur.name, cast(XUC) "author".ptr)) {
            gretl_xml_node_get_trimmed_string(cur, doc, &pkg.author);
            gretl_xml_get_prop_as_string(cur, "email", &pkg.email);
        } else if (!xmlStrcmp(cur.name, cast(XUC) "version".ptr)) {
            gretl_xml_node_get_trimmed_string(cur, doc, &pkg.version_);
        } else if (!xmlStrcmp(cur.name, cast(XUC) "date".ptr)) {
            gretl_xml_node_get_trimmed_string(cur, doc, &pkg.date);
        } else if (!xmlStrcmp(cur.name, cast(XUC) "description".ptr)) {
            gretl_xml_node_get_trimmed_string(cur, doc, &pkg.descrip);
        } else if (!xmlStrcmp(cur.name, cast(XUC) "help".ptr)) {
            gretl_xml_node_get_trimmed_string(cur, doc, &pkg.help);
            gretl_xml_get_prop_as_string(cur, "filename", &pkg.help_fname);
        } else if (!xmlStrcmp(cur.name, cast(XUC) "gui-help".ptr)) {
            gretl_xml_node_get_trimmed_string(cur, doc, &pkg.gui_help);
            gretl_xml_get_prop_as_string(cur, "filename", &pkg.gui_help_fname);
        } else if (!xmlStrcmp(cur.name, cast(XUC) "R-depends".ptr)) {
            gretl_xml_node_get_trimmed_string(cur, doc, &pkg.Rdeps);
        } else if (!xmlStrcmp(cur.name, cast(XUC) "sample-script".ptr)) {
            gretl_xml_node_get_trimmed_string(cur, doc, &pkg.sample);
            gretl_xml_get_prop_as_string(cur, "filename", &pkg.sample_fname);
        } else if (!xmlStrcmp(cur.name, cast(XUC) "tags".ptr)) {
            gretl_xml_node_get_trimmed_string(cur, doc, &pkg.tags);
        } else if (!xmlStrcmp(cur.name, cast(XUC) "label".ptr)) {
            gretl_xml_node_get_trimmed_string(cur, doc, &pkg.label);
        } else if (!xmlStrcmp(cur.name, cast(XUC) "menu-attachment".ptr)) {
            gretl_xml_node_get_trimmed_string(cur, doc, &pkg.mpath);
        } else if (!xmlStrcmp(cur.name, cast(XUC) "provider".ptr)) {
            gretl_xml_node_get_trimmed_string(cur, doc, &pkg.provider);
        } else if (!xmlStrcmp(cur.name, cast(XUC) "conflict".ptr)) {
            gretl_xml_node_get_trimmed_string(cur, doc, &pkg.conflict);
        } else if (!xmlStrcmp(cur.name, cast(XUC) "data-files".ptr)) {
            pkg.datafiles =
                gretl_xml_get_strings_array(cur, doc, &pkg.n_files,
                                            0, err);
        } else if (!xmlStrcmp(cur.name, cast(XUC) "depends".ptr)) {
            pkg.depends =
                gretl_xml_get_strings_array(cur, doc, &pkg.n_depends,
                                            0, err);
        }

        cur = cur.next;
    }

    if (get_funcs) {
        cur = node.children;
        while (cur != null && !*err) {
            if (!xmlStrcmp(cur.name, cast(XUC) "gretl-function".ptr)) {
                *err = read_ufunc_from_xml(cur, doc, pkg);
            }
            cur = cur.next;
        }
    }


    return pkg;
}

/* read aggregated file that may contain one or more function packages
   and one or more "loose" functions */

extern (C) int read_session_functions_file(const(char) * fname)
{
    Fnpkg * pkg = null;
    xmlDocPtr doc = null;
    xmlNodePtr node = null;
    xmlNodePtr cur;
    int err = 0;


    err = gretl_xml_open_doc_root(fname, "gretl-functions", &doc, &node);
    if (err) {
        return err;
    }


    /* get any function packages from this file */
    cur = node.children;
    while (cur != null && !err) {
        if (!xmlStrcmp(cur.name, cast(XUC) "gretl-function-package".ptr)) {
            pkg = real_read_package(doc, cur, fname, 1, &err);
            if (!err) {
                err = real_load_package(pkg, null, null, 0);
            }
        }
        cur = cur.next;
    }

    /* then get any unpackaged functions */
    if (!err) {
        cur = node.children;
        while (cur != null && !err) {
            if (!xmlStrcmp(cur.name, cast(XUC) "gretl-function".ptr)) {
                err = read_ufunc_from_xml(cur, doc, null);
            }
            cur = cur.next;
        }
    }

    if (doc != null) {
        xmlFreeDoc(doc);
    }


    return err;
}

/* Parse an XML function package file and return an allocated package
   struct with the functions attached. Note that this function does
   not actually "load" the package (making it available to users) and
   does not check dependencies.
*/

Fnpkg * read_package_file(const(char) * fname, int get_funcs, int * err)
{
    Fnpkg * pkg = null;
    xmlDocPtr doc = null;
    xmlNodePtr node = null;


    node = gretl_xml_get_gfn(fname, &doc, err);
    if (!*err) {
        pkg = real_read_package(doc, node, fname, get_funcs, err);
    }

    if (doc != null) {
        xmlFreeDoc(doc);
    }


    return pkg;
}

/* If @ps is non-null, we retrieve a single string from the package;
   otherwise if @pss is non-null we retrieve an array of strings and
   also write the length of this array into @ns.
*/

void package_peek_attribute(const(char) * fname, const(char) * tag, char ** ps, char *** pss, int * ns)
{
    xmlDocPtr doc = null;
    xmlNodePtr node = null;
    xmlNodePtr cur;
    int err = 0;

    node = gretl_xml_get_gfn(fname, &doc, &err);

    if (node != null) {
        cur = node.children;
        while (cur != null) {
            if (!xmlStrcmp(cur.name, cast(XUC) tag)) {
		if (pss != null) {
		    *pss = gretl_xml_get_strings_array(cur, doc, ns, 0, &err);
		} else {
		    *ps = gretl_xml_get_string(cur, doc);
		}
                break;
            }
            cur = cur.next;
        }
    }

    if (doc != null) {
        xmlFreeDoc(doc);
    }
}

extern (C) char ** package_peek_dependencies(const(char) * fname, int * ndeps)
{
    char **ss = null;

    package_peek_attribute(fname, "depends", null, &ss, ndeps);
    return ss;
}

extern (C) char * package_peek_conflict(const(char) * fname)
{
    char * s = null;

    package_peek_attribute(fname, "conflict", &s, null, null);
    return s;
}

/**
 * function_package_is_loaded:
 * @fname: full path to gfn file.
 * @version: location to receive version info, or null.
 *
 * Returns: 1 if the function package with filename @fname is
 * loaded in memory, otherwise 0.0
 */

extern (C) int function_package_is_loaded(const(char) * fname, const(char) ** version_)
{
    return (get_loaded_pkg_by_filename(fname, version_) != null);
}

/**
 * gfn_is_loaded:
 * @gfnname: basename of gfn file, including suffix.
 *
 * Returns: 1 if the function package with basename @gfnname is
 * loaded in memory, otherwise 0.0
 */

extern (C) int gfn_is_loaded(const(char) * gfnname)
{
    int ret = 0;

    if (strchr(gfnname, SLASH) == null && has_suffix(gfnname, ".gfn")) {
        /* plain basename of gfn file */
        gchar *test = g_strdup(gfnname);
        gchar *p = strrchr(test, '.');

        p[0] = '\0';
        if (get_function_package_by_name(test)) {
            ret = 1;
        }
        g_free(test);
    }

    return ret;
}

int not_mpi_duplicate()
{
    return 0;
}

Fnpkg * check_for_loaded(const(char) * fname, Gretlopt opt)
{
    Fnpkg *pkg = get_loaded_pkg_by_filename(fname, null);

    if (pkg != null) {
        if (opt & Gretlopt.OPT_F) {
            /* force re-reading from file */
            real_function_package_unload(pkg, 1);
            pkg = null;
        }
    }

    return pkg;
}

/**
 * load_function_package:
 * @fname: full path to gfn file.
 * @opt: may include Gretlopt.OPT_F to force loading even when
 * the package is already loaded.
 * @prn: gretl printer.
 * @ppkg: optional pointer to receive Fnpkg pointer.
 *
 * Loads the function package located by @fname into
 * memory, if possible. Supports gretl's "include" command
 * for gfn files.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int load_function_package(const(char) * fname, Gretlopt opt, GArray * pstack, Prn * prn, Fnpkg ** ppkg, int level)
{
    Fnpkg * pkg;
    int err = 0;

    pkg = check_for_loaded(fname, opt);
    if (pkg != null) {
        return 0;
    }

    pkg = read_package_file(fname, 1, &err);

    if (!err) {
        if (ppkg != null) {
            *ppkg = pkg;
        }
        if (pkg.Rdeps != null) {
            err = check_R_depends(pkg.name, pkg.Rdeps, prn);
        }
    }

    if (!err) {
        /* Let's double-check that we don't have a colliding package
           (it would have to be with a different filename).
        */
        Fnpkg * oldpkg;

        oldpkg = get_function_package_by_name(pkg.name);
        if (oldpkg != null) {
            real_function_package_unload(oldpkg, 1);
        }
        err = real_load_package(pkg, pstack, prn, level);
    }

    if (err) {
	gchar *tmp = g_strdup_printf(_("Error loading %s\n").ptr, fname);

	gretl_errmsg_prepend(tmp, err);
	g_free(tmp);
    } else if (pkg != null && prn != null && not_mpi_duplicate()) {
        pprintf(prn, "%s %s, %s (%s)\n", pkg.name, pkg.version_,
                pkg.date, pkg.author);
    }

    return err;
}

/**
 * include_gfn:
 * @fname: full path to gfn file.
 * @opt: may include Gretlopt.OPT_F to force loading even when
 * the package is already loaded.
 * @prn: gretl printer.
 *
 * Loads the function package located by @fname into
 * memory, if possible. Supports gretl's "include" command
 * for gfn files.
 *
 * Returns: 0 on success, non-zero code on error.
 */

extern (C) int include_gfn(const(char) * fname, Gretlopt opt, Prn * prn)
{
    GArray * pstack;
    char * p;
    char * pkgname;
    int err;

    pkgname = g_path_get_basename(fname);
    p = strrchr(pkgname, '.');
    if (p != null) {
        p[0] = '\0';
    }

    pstack = g_array_new(0, 0, (char *).sizeof);
    {
        auto tmp_ = pkgname;

        g_array_append_vals(pstack, &tmp_, 1);
    }
    err = load_function_package(fname, opt, pstack, prn, null, 0);
    g_array_free(pstack, 1);
    g_free(pkgname);

    return err;
}

/**
 * get_function_package_by_filename:
 * @fname: gfn filename.
 * @err: location to receive error code.
 *
 * If the package whose filename is @fname is already loaded,
 * returns the package pointer, otherwise attempts to load the
 * package from file. Similar to include_gfn() but returns
 * the package pointer rather than just a status code.
 *
 * Returns: package-pointer on success, null on error.
 */

extern (C) Fnpkg * get_function_package_by_filename(const(char) * fname, int * err)
{
    Fnpkg * pkg = null;
    int i;
    int myerr = 0;

    for (i=0; i<n_pkgs; i++) {
        if (!strcmp(fname, pkgs[i].fname)) {
            pkg = pkgs[i];
            break;
        }
    }

    if (pkg == null) {
        Prn *prn = gretl_print_new(PrnType.buffer, null);

        myerr = include_gfn(fname, Gretlopt.OPT_NONE, prn);
        if (!myerr) {
            for (i=0; i<n_pkgs; i++) {
                if (!strcmp(fname, pkgs[i].fname)) {
                    pkg = pkgs[i];
                    break;
                }
            }
        } else if (prn != null) {
            const(char) * buf = gretl_print_get_buffer(prn);

            gretl_errmsg_set(buf);
        }
        gretl_print_destroy(prn);
    }

    if (err != null) {
        *err = myerr;
    }

    return pkg;
}

/**
 * get_function_package_by_name:
 * @pkgname: name of function package.
 *
 * Returns: pointer to function package if a package named
 * @pkgname is already in memory, otherwise null.
 */

extern (C) Fnpkg * get_function_package_by_name(const(char) * pkgname)
{
    int i;

    if (has_suffix(pkgname, ".gfn") || has_suffix(pkgname, ".zip")) {
        /* just in case: strip off the extension */
        gchar *tmp = g_strdup(pkgname);
        gchar *p = strrchr(tmp, '.');
        Fnpkg * ret = null;

        p[0] = '\0';
        for (i=0; i<n_pkgs; i++) {
            if (!strcmp(tmp, pkgs[i].name)) {
                ret = pkgs[i];
                break;
            }
        }
        g_free(tmp);
        return ret;
    } else {
        for (i=0; i<n_pkgs; i++) {
            if (!strcmp(pkgname, pkgs[i].name)) {
                return pkgs[i];
            }
        }
        return null;
    }
}

/**
 * get_function_package_path_by_name:
 * @pkgname: name of function package.
 *
 * Returns: the full path to the named function package if it
 * is already in memory, otherwise null.
 */

extern (C) const(char) * get_function_package_path_by_name(const(char) * pkgname)
{
    int i;

    for (i=0; i<n_pkgs; i++) {
        if (!strcmp(pkgname, pkgs[i].name)) {
            return pkgs[i].fname;
        }
    }

    return null;
}

int gfn_version_fail(const(char) * fname, Prn * prn)
{
    FILE *fp = gretl_fopen(fname, "r");
    int err = 0;

    if (fp != null) {
        char[128] line;

        if (fgets(line.ptr, line.sizeof, fp) &&
            strchr(line.ptr, '<') == null &&
            strstr(line.ptr, "requires gretl".ptr) != null) {
            gretl_errmsg_set(gretl_strstrip(line.ptr));
            err = 1;
        }
        fclose(fp);
    }

    return err;
}

char * sample_script_from_xml(const(char) * pkgname, int * err)
{
    char * ret = null;
    char * gfnname = null;
    char[FILENAME_MAX] fullname;
    xmlDocPtr doc = null;
    xmlNodePtr node = null;

    gfnname = g_strdup_printf("%s.gfn", pkgname);
    *err = get_full_filename(gfnname, fullname.ptr, Gretlopt.OPT_I);

    if (!*err) {
        node = gretl_xml_get_gfn(fullname.ptr, &doc, err);
    }

    if (!*err) {
	xmlNodePtr cur = node.children;

	while (cur != null) {
            if (!xmlStrcmp(cur.name, cast(XUC) "sample-script".ptr)) {
                gretl_xml_node_get_trimmed_string(cur, doc, &ret);
                break;
            }
            cur = cur.next;
        }
    }

    g_free(gfnname);
    if (doc != null) {
	xmlFreeDoc(doc);
    }

    return ret;
}

extern (C) int grab_package_sample(const(char) * pkgname, char ** pscript)
{
    int err = 0;

    if (pkgname == null) {
	err = GretlError.data;
    } else {
	Fnpkg *pkg = get_loaded_pkg_by_name(pkgname);
	char * ss = null;

	if (pkg != null) {
	    ss = gretl_strdup(pkg.sample);
	} else {
	    ss = sample_script_from_xml(pkgname, &err);
	}
	*pscript = ss;
	if (ss == null) {
	    gretl_errmsg_sprintf(_("Couldn't find sample script for '%s'").ptr, pkgname);
	    err = GretlError.data;
	}
    }

    return err;
}

/* Retrieve summary info, sample script, or code listing for a
   function package, identified by its filename.  This is called
   (indirectly) from the GUI -- see below for the actual callbacks.
*/

int real_print_gfn_data(const(char) * fname, char ** pbuf, Prn * prn, int tabwidth, int task, GretlBundle * b)
{
    Fnpkg * pkg = null;
    int free_pkg = 0;
    int err = 0;


    if (task == FuncsTask.info && prn != null && strstr(fname, "dltmp.")) {
        if (gfn_version_fail(fname, prn)) {
            return 1;
        }
    }

    pkg = get_loaded_pkg_by_filename(fname, null);


    if (pkg == null) {
        /* package is not loaded, read it now */
        int get_funcs = 1;

        if (task == FuncsTask.help || task == FuncsTask.info ||
            task == FuncsTask.sample || task == FuncsTask.query) {
            get_funcs = 0;
        }
        pkg = read_package_file(fname, get_funcs, &err);
        free_pkg = 1;
    }

    if (!err) {
        if (task == FuncsTask.help) {
            print_package_help(pkg, fname, pbuf, prn);
        } else if (task == FuncsTask.info) {
            print_package_info(pkg, fname, prn);
        } else if (task == FuncsTask.query) {
            if (b != null) {
                real_bundle_package_info(pkg, fname, b);
            } else {
                plain_print_package_info(pkg, fname, prn);
            }
        } else if (task == FuncsTask.sample) {
            pputs(prn, pkg.sample);
        } else if (task == FuncsTask.code) {
            /* called by print_function_package_code() */
            print_package_code(pkg, tabwidth, prn);
        }
        if (free_pkg) {
            function_package_free_full(pkg);
        }
    }


    return err;
}

/* callback used in the GUI function package browser */

extern (C) int print_function_package_info(const(char) * fname, int gui_mode, Prn * prn)
{
    int mode = gui_mode ? FuncsTask.info : FuncsTask.query;

    return real_print_gfn_data(fname, null, prn, 0, mode, null);
}

/* callback used by "pkg" command with query action + --quiet */

extern (C) int bundle_function_package_info(const(char) * fname, GretlBundle * b)
{
    return real_print_gfn_data(fname, null, null, 0, FuncsTask.query, b);
}

/* callback used in the GUI function package browser (gui/datafiles.c) */

extern (C) int print_function_package_code(const(char) * fname, int tabwidth, Prn * prn)
{
    return real_print_gfn_data(fname, null, prn, tabwidth, FuncsTask.code, null);
}

/* callback used in the GUI function package browser */

extern (C) int print_function_package_sample(const(char) * fname, int tabwidth, Prn * prn)
{
    return real_print_gfn_data(fname, null, prn, tabwidth, FuncsTask.sample, null);
}

/* callback used via command line */

extern (C) int print_function_package_help(const(char) * fname, char ** pbuf, Prn * prn)
{
    return real_print_gfn_data(fname, pbuf, prn, 0, FuncsTask.help, null);
}

void maybe_fix_broken_date(char ** pdate)
{
    char * date = *pdate;

    if (strlen(date) < 10 && date[4] == '-') {
        int y;
        int m;
        int d;

        if (sscanf(date, "%d-%d-%d", &y, &m, &d) == 3) {
            char[12] tmp;

            sprintf(tmp.ptr, "%d-%02d-%02d".ptr, y, m, d);
            free(*pdate);
            *pdate = gretl_strdup(tmp.ptr);
        }
    }
}

int get_pdfdoc_status(xmlNodePtr sub, xmlDocPtr doc)
{
    char * tmp = null;
    int ret = 0;

    gretl_xml_node_get_trimmed_string(sub, doc, &tmp);
    if (tmp != null) {
        if (!strncmp(tmp, "pdfdoc", 6) || is_pdf_ref(tmp)) {
            ret = 1;
        }
        free(tmp);
    }

    return ret;
}

/* Read the header from a function package file -- this is used when
   displaying the available packages in the GUI (see gui/datafiles.c).
   We write the description into *pdesc, a string representation of
   version number into *pver, date into *pdate, and author's name into
   *pauthor.

   The *pdfdoc pointer can be used to collect an int value indicating
   whether the package contains PDF documentation (1/0), or it can be
   set to null. None of the other arguments can be null.
*/

extern (C) int get_function_file_header(const(char) * fname, char ** pdesc, char ** pver, char ** pdate, char ** pauthor, int * pdfdoc)
{
    xmlDocPtr doc = null;
    xmlNodePtr node = null;
    xmlNodePtr cur;
    int docdone = 0;
    int ndone = 0;
    int err = 0;

    if (pdesc == null || pver == null || pdate == null || pauthor == null) {
        fprintf(stderr, "get_function_file_header: missing parameter\n");
        return GretlError.data;
    }

    node = gretl_xml_get_gfn(fname, &doc, &err);
    if (err) {
        return err;
    }

    if (pdfdoc == null) {
        /* not wanted, so count it as "done" already */
        docdone = 1;
    } else {
        *pdfdoc = 0;
    }

    cur = node.children;
    while (cur != null) {
        if (!xmlStrcmp(cur.name, cast(XUC) "description".ptr)) {
            gretl_xml_node_get_trimmed_string(cur, doc, pdesc);
            ndone++;
        } else if (!xmlStrcmp(cur.name, cast(XUC) "version".ptr)) {
            gretl_xml_node_get_trimmed_string(cur, doc, pver);
            ndone++;
        } else if (!xmlStrcmp(cur.name, cast(XUC) "date".ptr)) {
            gretl_xml_node_get_trimmed_string(cur, doc, pdate);
            if (*pdate != null) {
                maybe_fix_broken_date(pdate);
            }
            ndone++;
        } else if (!xmlStrcmp(cur.name, cast(XUC) "author".ptr)) {
            gretl_xml_node_get_trimmed_string(cur, doc, pauthor);
            ndone++;
        } else if (!docdone && !xmlStrcmp(cur.name, cast(XUC) "help".ptr)) {
            *pdfdoc = get_pdfdoc_status(cur, doc);
            docdone = 1;
        }
        if (docdone && ndone == 4) {
            break;
        }
        cur = cur.next;
    }

    if (doc != null) {
        xmlFreeDoc(doc);
    }

    if (*pdesc == null) {
        *pdesc = gretl_strdup(_("No description available").ptr);
    }
    if (*pver == null) {
        *pver = gretl_strdup("unknown");
    }
    if (*pdate == null) {
        *pdate = gretl_strdup("unknown");
    }
    if (*pauthor == null) {
        *pauthor = gretl_strdup("unknown");
    }

    if (*pdesc == null || *pver == null ||
        *pdate == null || *pauthor == null) {
        err = GretlError.alloc;
    }

    return err;
}

extern (C) int package_has_menu_attachment(const(char) * fname, char ** pkgname, char ** attach, char ** label)
{
    xmlDocPtr doc = null;
    xmlNodePtr node = null;
    xmlNodePtr cur;
    char * tmp = null;
    int skipit = 0;
    int got_label = 0;
    int got_attach = 0;
    int err = 0;

    node = gretl_xml_get_gfn(fname, &doc, &err);
    if (err) {
        return 0;
    }

    if (pkgname != null) {
        gretl_xml_get_prop_as_string(node, "name", &tmp);
        if (tmp != null && !strcmp(tmp, "ridge")) {
            /* don't give "ridge" a menu attachment */
            free(tmp);
            skipit = 1;
        } else {
            *pkgname = tmp;
        }
    }

    if (skipit) {
        xmlFreeDoc(doc);
        return 0;
    }

    cur = node.children;
    while (cur != null) {
        if (!xmlStrcmp(cur.name, cast(XUC) "menu-attachment".ptr)) {
            gretl_xml_node_get_trimmed_string(cur, doc, &tmp);
            if (tmp != null) {
                got_attach = 1;
                if (attach != null) {
                    *attach = tmp;
                } else {
                    free(tmp);
                }
            }
        } else if (!xmlStrcmp(cur.name, cast(XUC) "label".ptr)) {
            gretl_xml_node_get_trimmed_string(cur, doc, &tmp);
            if (tmp != null) {
                got_label = 1;
                if (label != null) {
                    *label = tmp;
                } else {
                    free(tmp);
                }
            }
        } else if (!xmlStrcmp(cur.name, cast(XUC) "help".ptr)) {
            /* we've overshot */
            break;
        }
        if (got_attach && got_label) {
            break;
        }
        cur = cur.next;
    }

    if (doc != null) {
        xmlFreeDoc(doc);
    }

    return got_attach && got_label;
}

extern (C) int package_needs_zipping(const(char) * fname, int * pdfdoc, char *** datafiles, int * n_files)
{
    xmlDocPtr doc = null;
    xmlNodePtr node = null;
    xmlNodePtr cur;
    char * tmp = null;
    int retmax = 1;
    int ret = 0;
    int err = 0;

    node = gretl_xml_get_gfn(fname, &doc, &err);
    if (err) {
        return 0;
    }

    if (datafiles != null) {
        retmax = 2;
    }

    cur = node.children;
    while (cur != null) {
        if (!xmlStrcmp(cur.name, cast(XUC) "help".ptr)) {
            gretl_xml_node_get_trimmed_string(cur, doc, &tmp);
            if (tmp != null && !strncmp(tmp, "pdfdoc:", 7)) {
                if (pdfdoc != null) {
                    *pdfdoc = 1;
                }
                ret++;
            }
            free(tmp);
        } else if (!xmlStrcmp(cur.name, cast(XUC) "data-files".ptr)) {
            if (datafiles != null) {
                *datafiles = gretl_xml_get_strings_array(cur, doc, n_files,
                                                         0, &err);
            }
            ret++;
        } else if (!xmlStrcmp(cur.name, cast(XUC) "gretl-function".ptr)) {
            /* we've overshot */
            break;
        }
        if (ret == retmax) {
            break;
        }
        cur = cur.next;
    }

    if (doc != null) {
        xmlFreeDoc(doc);
    }

    return ret;
}

extern (C) double gfn_file_get_version(const(char) * fname)
{
    xmlDocPtr doc = null;
    xmlNodePtr node;
    double version_ = NADBL;
    int err = 0;

    node = gretl_xml_get_gfn(fname, &doc, &err);

    if (!err) {
        xmlNodePtr cur = node.children;

        while (cur != null) {
            if (!xmlStrcmp(cur.name, cast(XUC) "version".ptr)) {
                gretl_xml_node_get_double(cur, doc, &version_);
                break;
            }
            cur = cur.next;
        }
    }

    if (doc != null) {
        xmlFreeDoc(doc);
    }

    return version_;
}

/**
 * gretl_is_public_user_function:
 * @name: name to test.
 *
 * Returns: 1 if @name is the name of a user-defined
 * function that is not a private member of a function
 * package, otherwise 0.0
 */

extern (C) int gretl_is_public_user_function(const(char) * name)
{
    Ufunc *fun = get_user_function_by_name(name);

    return (fun != null && !function_is_private(fun));
}

extern (C) void set_current_function_package(Fnpkg * pkg)
{
    current_pkg = pkg;
}

int skip_private_func(Ufunc * ufun)
{
    int skip = 0;

    if (function_is_private(ufun)) {
        /* skip it, unless we're "authorized" to edit it,
           i.e. coming from the gui package editor */
        skip = 1;
        if (ufun.pkg != null && ufun.pkg == current_pkg) {
            skip = 0;
        }
    }

    return skip;
}

int check_function_name(const(char) * name, Ufunc ** pfun, const(DatasetRaw) * dset, Prn * prn)
{
    int i;
    int err = 0;


    if (!isalpha(cast(ubyte) *name)) {
        gretl_errmsg_set(_("Function names must start with a letter").ptr);
        err = 1;
    } else if (gretl_reserved_word(name)) {
        err = 1;
    } else if (function_lookup(name)) {
        gretl_errmsg_sprintf(_("'%s' is the name of a built-in function").ptr, name);
        err = 1;
    } else if (gretl_is_user_var(name)) {
        gretl_errmsg_sprintf(_("'%s' is the name of a variable").ptr, name);
        err = 1;
    } else if (gretl_is_series(name, dset)) {
        gretl_errmsg_sprintf(_("'%s' is the name of a variable").ptr, name);
        err = 1;
    } else {
        /* @name is OK, now check for existing function of the same name */
        for (i=0; i<n_ufuns; i++) {
            if (!skip_private_func(ufuns[i]) && !strcmp(name, ufuns[i].name)) {
                if (ufuns[i].pkg != null && ufuns[i].pkg != current_pkg) {
                    /* don't allow overwriting */
                    gretl_errmsg_sprintf(_("The function %s is already defined " ~
		"by package '%s'").ptr, name,
                                         ufuns[i].pkg.name);
                    err = GretlError.data;
                } else if (pfun != null) {
                    clear_ufunc_data(ufuns[i]);
                    *pfun = ufuns[i];
                } else {
                    ufunc_unload(ufuns[i]);
                }
                break;
            }
        }
    }

    return err;
}

int maybe_delete_function(const(char) * fname, Prn * prn)
{
    Ufunc *fun = get_user_function_by_name(fname);
    int err = 0;

    if (fun == null) {
        ; /* no-op */
    } else if (function_in_use(fun)) {
        gretl_errmsg_sprintf(_("%s: function is in use").ptr, fname);
        err = 1;
    } else if (fun.pkg != null) {
        gretl_errmsg_sprintf(_("%s: function belongs to package").ptr, fname);
        err = 1;
    } else {
        ufunc_unload(fun);
        if (gretl_messages_on()) {
            pprintf(prn, _("Deleted function '%s'\n").ptr, fname);
        }
    }

    return err;
}

/* next: apparatus for parsing function definitions */

int comma_count(const(char) * s)
{
    int quoted = 0;
    int braced = 0;
    int nc = 0;

    while (*s) {
        if (!quoted) {
            if (*s == '{') {
                braced++;
            } else if (*s == '}') {
                braced--;
            }
        }
        if (*s == '"') {
            quoted = !quoted;
        } else if (!quoted && !braced && *s == ',') {
            nc++;
        }
        s++;
    }

    return nc;
}

int check_parm_min_max(FnParam * p, const(char) * name, int * nvals)
{
    int err = 0;

    if (p.type != GretlType.double_ && na(p.deflt)) {
        gretl_errmsg_sprintf(_("%s: parameters of type %s cannot have a " ~
		"default value of NA").ptr, name,
                             gretl_type_get_name(p.type));
        err = GretlError.data;
    }

    if (!err && !na(p.min) && !na(p.max)) {
        if (p.min > p.max) {
            gretl_errmsg_sprintf("%s: min > max", name);
            err = GretlError.data;
        } else if (!na(p.step) && p.step > p.max - p.min) {
            gretl_errmsg_sprintf("%s: step > max - min", name);
            err = GretlError.data;
        }
        if (!err) {
            *nvals = cast(int) p.max - cast(int) p.min + 1;
        }
    }

    if (!err && !na(p.deflt) && !default_unset(p)) {
        if (!na(p.min) && p.deflt < p.min) {
            gretl_errmsg_sprintf(_("%s: default value out of bounds").ptr, name);
            err = GretlError.data;
        } else if (!na(p.max) && p.deflt > p.max) {
            gretl_errmsg_sprintf(_("%s: default value out of bounds").ptr, name);
            err = GretlError.data;
        }
    }

    return err;
}

int colon_count(const(char) * p)
{
    int n = 0;

    while (*p && *p != ']') {
        if (*p == ':') {
            n++;
        }
        p++;
    }

    return n;
}


int read_min_max_deflt(char ** ps, FnParam * param, const(char) * name, int * nvals)
{
    char * p = *ps;
    int err = 0;

    gretl_push_c_numeric_locale();

    errno = 0;

    if (!strcmp(p, "[null]")) {
        gretl_errmsg_set(_("'null' is not a valid default for scalars").ptr);
        err = GretlError.types;
    } else if (!strcmp(p, "[auto]")) {
        /* FIXME */
        err = GretlError.types;
    } else if (param.type == GretlType.bool_) {
        if (!strncmp(p, "[TRUE]", 6)) {
            param.deflt = 1;
        } else if (!strncmp(p, "[FALSE]", 7)) {
            param.deflt = 0;
        } else if (sscanf(p, "[%lf]", &param.deflt) != 1) {
            err = GretlError.parse;
        }
    } else if (!strncmp(p, "[$xlist]", 8)) {
        param.deflt = INT_USE_XLIST;
    } else if (!strncmp(p, "[$mylist]", 9)) {
        param.deflt = INT_USE_MYLIST;
    } else {
        int i, len, nf = colon_count(p) + 1;
        char * test;
        char[32] valstr;
        char * q = p + 1;
        double[4] x = [NADBL, NADBL, UNSET_VALUE, NADBL];

        for (i=0; i<nf && !err; i++) {
            len = cast(int) strcspn(q, ":]");
            if (len > 31) {
                err = GretlError.parse;
            } else if (len > 0) {
                valstr[0] = '\0';
                strncat(valstr.ptr, q, len);
                if (!strcmp(valstr.ptr, "NA".ptr)) {
                    x[i] = NADBL;
                } else {
                    x[i] = strtod(valstr.ptr, &test);
                    if (*test != '\0' || errno) {
                        err = GretlError.parse;
                    }
                }
            }
            q += len + 1;
        }

        if (!err) {
            if (nf == 1) {
                /* we take a single value with no colons as
                   indicating a default */
                param.deflt = x[0];
            } else {
                param.min = x[0];
                param.max = x[1];
                param.deflt = x[2];
                param.step = x[3];
            }
            err = check_parm_min_max(param, name, nvals);
        }
    }

    gretl_pop_c_numeric_locale();

    if (!err) {
        p = strchr(p, ']');
        if (p == null) {
            err = GretlError.parse;
        } else {
            *ps = p + 1;
        }
    }

    return err;
}

int read_param_option(char ** ps, FnParam * param)
{
    int err = 0;


    if (!strncmp(*ps, "[auto]", 6)) {
        set_param_auto(param);
        *ps += 6;
    } else if (!strncmp(*ps, "[null]", 6)) {
        set_param_optional(param);
        *ps += 6;
    } else {
        gretl_errmsg_sprintf(_("got invalid field '%s'").ptr, *ps);
        err = GretlError.parse;
    }

    return err;
}

/* get the descriptive string for a function parameter */

int read_param_descrip(char ** ps, FnParam * param)
{
    int offset = 1;
    char * p = *ps + offset;
    int len = 0;
    int err = GretlError.parse;

    while (*p) {
        if (*p == '"') {
            /* found closing quote */
            err = 0;
            break;
        }
        len++;
        p++;
    }

    if (!err && len > 0) {
        p = *ps + offset;
        param.descrip = gretl_strndup(p, len);
        if (param.descrip == null) {
            err = GretlError.alloc;
        } else {
            *ps = p + len + 1;
        }
    }

    return err;
}

/* Get the value labels for a function parameter. The syntactic
   element we're looking at starts with "{".
*/

int read_param_labels(char ** ps, FnParam * param, const(char) * name, int nvals)
{
    char * p = *ps;
    int offset = 1;
    int len = 0;
    int err = GretlError.parse;

    p += offset;

    while (*p) {
        if (*p == '}') {
            /* OK, found closing brace */
            err = 0;
            break;
        }
        len++;
        p++;
    }

    if (!err && len > 0) {
        char * tmp;

        p = *ps + offset;
        tmp = gretl_strndup(p, len);

        if (tmp == null) {
            err = GretlError.alloc;
        } else {
            *ps = p + len + 1;
            param.labels = gretl_string_split_quoted(tmp, &param.nlabels,
                                                      " ,", &err);
            free(tmp);
            if (!err && param.nlabels != nvals) {
                gretl_errmsg_sprintf(_("%s: found %d values but %d value-labels").ptr,
                                     name, nvals, param.nlabels);
                err = GretlError.data;
            }
        }
    }

    return err;
}

void trash_param_info(char * name, FnParam * param)
{
    free(name);
    free(param.descrip);
    param.descrip = null;
    if (param.nlabels > 0) {
        strings_array_free(param.labels, param.nlabels);
        param.labels = null;
        param.nlabels = 0;
    }
}

int parse_function_param(char * s, FnParam * param, int i)
{
    GretlType type = GretlType.none;
    char[22] tstr = 0;
    char * name;
    int len;
    int nvals = 0;
    int const_flag = 0;
    int err = 0;

    while (isspace(*s)) s++;


    /* pick up the "const" flag if present */
    if (!strncmp(s, "const ", 6)) {
        const_flag = 1;
        s += 6;
        while (isspace(*s)) s++;
    }

    /* get parameter type, required */
    len = gretl_namechar_spn(s);
    if (len > 21) {
        err = GretlError.parse;
    } else {
        strncat(tstr.ptr, s, len);
        s += len;
        while (isspace(*s)) s++;
        if (*s == '*') {
            strcat(tstr.ptr, " *".ptr);
            s++;
        }
        if (*tstr.ptr == '\0') {
            gretl_errmsg_set(_("Expected a type identifier").ptr);
            err = GretlError.parse;
        } else {
            type = gretl_type_from_string(tstr.ptr);
            if (type == 0) {
                gretl_errmsg_sprintf(_("Unrecognized data type '%s'").ptr, tstr.ptr);
                err = GretlError.parse;
            } else if (!ok_function_arg_type(type)) {
                gretl_errmsg_sprintf(_("%s: invalid parameter type").ptr, tstr.ptr);
                err = GretlError.invarg;
            }
        }
    }

    if (err) {
        return err;
    }

    /* get the required parameter name */
    while (isspace(*s)) s++;
    len = gretl_namechar_spn(s);
    if (len == 0) {
        gretl_errmsg_sprintf(_("parameter %d: name is missing").ptr, i + 1);
        err = GretlError.parse;
    } else if (len >= VNAMELEN) {
	gretl_errmsg_set(_("Identifier exceeds the maximum of 31 characters").ptr);
	err = GretlError.parse;
    } else {
        name = gretl_strndup(s, len);
        if (name == null) {
            err = GretlError.alloc;
        } else if (gretl_reserved_word(name)) {
            free(name);
            err = GretlError.data;
        }
    }

    if (err) {
        return err;
    }

    param.type = type;
    if (const_flag) {
        maybe_set_param_const(param);
    }

    s += len;
    s += strspn(s, " ");

    /* now scan for various optional extras: first we may have
       a [min:max:default] specification (scalars only), or a
       [null|auto] spec to allow no argument for "pointer"-type
       parameters (broadly interpreted)
    */

    if (*s == '[') {
        if (gretl_scalar_type(type)) {
            err = read_min_max_deflt(&s, param, name, &nvals);
        } else if (arg_may_be_optional(type)) {
            err = read_param_option(&s, param);
        } else {
            gretl_errmsg_sprintf(_("'%s': error scanning default value").ptr,
                                 name);
            err = GretlError.parse;
        }
        s += strspn(s, " ");
    }

    /* then we may have a double-quoted descriptive string
       for the parameter */

    if (!err && *s == '"') {
        err = read_param_descrip(&s, param);
        s += strspn(s, " ");
    }

    /* and finally we may have a set of value-labels enclosed
       in braces, but this is valid only if we have been able
       to determine a definite number of admissible values
       for the parameter
    */

    if (!err && *s == '{') {
        if (nvals == 0) {
            err = GretlError.parse;
        } else {
            err = read_param_labels(&s, param, name, nvals);
        }
        s += strspn(s, " ");
    }

    if (!err && *s != '\0') {
        /* got trailing unparseable stuff */
        err = GretlError.parse;
    }

    if (!err) {
        param.name = name;
    } else {
        trash_param_info(name, param);
    }


    return err;
}

void arg_tail_strip(char * s)
{
    int i, n = cast(int) strlen(s);

    for (i=n-1; i>=0; i--) {
        if (isspace(s[i])) {
            s[i] = '\0';
        } else if (s[i] == ')') {
            s[i] = '\0';
            break;
        } else {
            break;
        }
    }
}

int parse_function_parameters(const(char) * str, FnParam ** pparams, int * pnp, const(char) * fname, Prn * prn)
{
    FnParam * params = null;
    char * p;
    char * s;
    int i;
    int j;
    int np = 0;
    int err = 0;

    s = gretl_strdup(str);
    if (s == null) {
        return GretlError.alloc;
    }

    if (!strcmp(s, ")")) {
        /* we got a void function "foo()" */
        free(s);
        return 0;
    }

    /* strip trailing ')' and space */
    arg_tail_strip(s);
    np = comma_count(s) + 1;
    if (np == 1 && !strcmp(s, "void")) {
        free(s);
        return 0;
    }


    params = allocate_params(np);
    if (params == null) {
        err = GretlError.alloc;
    }

    if (!err) {
        int quoted = 0;
        int braced = 0;

        p = s;
        while (*p) {
            if (!quoted) {
                if (*p == '{') {
                    braced++;
                } else if (*p == '}') {
                    braced--;
                }
            }
            if (*p == '"') {
                quoted = !quoted;
            } else if (!quoted && !braced && *p == ',') {
                p[0] = '\0';
            }
            p++;
        }
        p = s;
        if (braced != 0) {
            err = GretlError.parse;
        }
        for (i=0; i<np && !err; i++) {
            err = parse_function_param(p, &params[i], i);
            p += strlen(p) + 1;
        }
    }

    free(s);

    for (i=0; i<np && !err; i++) {
        for (j=i+1; j<np && !err; j++) {
            if (!strcmp(params[i].name, params[j].name)) {
                gretl_errmsg_sprintf(_("%s: duplicated parameter name '%s'").ptr,
                                     fname, params[i].name);
                err = GretlError.data;
            }
        }
    }

    if (err) {
        free_params_array(params, np);
    } else {
        *pparams = params;
        *pnp = np;
    }

    return err;
}

int get_two_words(const(char) * s, char * w1, char * w2, int * err)
{
    int nf = 0;

    if (strncmp(s, "function ", 9)) {
        *err = GretlError.parse;
    } else {
        int n;

        *w1 = *w2 = '\0';

        s += 9;
        s += strspn(s, " ");

        /* @w1 should be a return type, except in the case
           of "function foo delete"
        */
        n = cast(int) strcspn(s, " ");

        if (n == 0 || n >= FN_NAMELEN) {
            *err = GretlError.parse;
        } else {
            strncat(w1, s, n);
            nf++;
            s += n;
            s += strspn(s, " ");
            n = cast(int) strcspn(s, " (");
        }

        /* @w2 should generally be the function name */
        if (n > 0 && n < FN_NAMELEN) {
            strncat(w2, s, n);
            nf++;
        } else if (n >= FN_NAMELEN) {
            gretl_errmsg_set(_("Identifier exceeds the maximum of 31 characters").ptr);
            *err = GretlError.parse;
        }
    }

    return nf;
}

/**
 * gretl_start_compiling_function:
 * @line: command line.
 * @prn: printing struct for feedback.
 *
 * Responds to a statement of the form "function ...".  In most
 * cases, embarks on compilation of a function, but this
 * also handles the construction "function foo delete".
 *
 * Returns: 0 on success, non-zero on error.
 */

extern (C) int gretl_start_compiling_function(const(char) * line, const(DatasetRaw) * dset, Prn * prn)
{
    Ufunc * fun = null;
    FnParam * params = null;
    int nf;
    int n_params = 0;
    GretlType rettype = GretlType.none;
    char[FN_NAMELEN] s1;
    char[FN_NAMELEN] s2;
    const(char) * p = null;
    char * name = null;
    int err = 0;

    if (gretl_function_depth() > 0) {
        return GretlError.fnest;
    }

    nf = get_two_words(line, s1.ptr, s2.ptr, &err);

    if (err) {
        return err;
    } else if (nf < 2) {
        gretl_errmsg_set(_("A function definition must have a return type and name").ptr);
        return GretlError.parse;
    }

    if (!strcmp(s2.ptr, "clear".ptr) || !strcmp(s2.ptr, "delete".ptr)) {
        return maybe_delete_function(s1.ptr, prn);
    }

    /* If we didn't get a special such as "function foo delete",
       then @s1.ptr should contain the return type and @s2.ptr the
       name.
    */

    name = s2.ptr;
    rettype = return_type_from_string(s1.ptr, &err);

    if (!err) {
        /* note: this handles a name collision */
        err = check_function_name(name, &fun, dset, prn);
    }

    if (!err) {
        /* now for the args bit */
        p = strchr(line, '(');
        if (p == null || strchr(p, ')') == null) {
            err = GretlError.parse;
        } else {
	    /* skip the left paren and any following spaces */
            p++;
	    p += strspn(p, " \t");
        }
    }

    if (!err) {
        err = parse_function_parameters(p, &params, &n_params,
                                        name, prn);
    }

    if (err) {
        pprintf(prn, "> %s\n", line);
    }

    if (!err && fun == null) {
        fun = add_ufunc(name);
        if (fun == null) {
            free_params_array(params, n_params);
            err = GretlError.alloc;
        }
    }


    if (!err) {
        strcpy(fun.name, name);
        fun.params = params;
        fun.n_params = n_params;
        fun.rettype = rettype;
        current_fdef = fun;
        set_compiling_on();
    } else {
        current_fdef = null;
    }

    return err;
}

void python_check(const(char) * line)
{
    char[8] s1;
    char[16] s2;

    if (sscanf(line, "%7s %15s", s1.ptr, s2.ptr) == 2) {
        if (!strcmp(s1.ptr, "foreign".ptr) && strstr(s2.ptr, "ytho".ptr)) {
            compiling_python = 1;
        } else if (!strcmp(s1.ptr, "end".ptr) && !strcmp(s2.ptr, "foreign".ptr)) {
            compiling_python = 0;
        }
    }
}

int ufunc_get_structure(Ufunc * u)
{
    int ret;
    int recurses = 0;

    ret = statements_get_structure(u.lines, u.n_lines,
				   u.name, &recurses);
    if (recurses) {
	u.flags |= UfunAttrs.recurses;
    }

    return ret;
}

bool NEEDS_IF(int c)
{
    return c == GretlCmdIndex.else_ || c == GretlCmdIndex.elif || c == GretlCmdIndex.endif;
}

bool FLOW_CI(int c)
{
    return c == GretlCmdIndex.if_ || c == GretlCmdIndex.else_ || c == GretlCmdIndex.elif || c == GretlCmdIndex.endif;
}

/**
 * gretl_function_append_line:
 * @s: pointer to execution state.
 *
 * Continuation of definition of user-function.
 *
 * Returns: 0 on success, non-zero on error.
 */

extern (C) int gretl_function_append_line(ExecState * state)
{
    static __gshared int ifdepth;
    static __gshared int last_flow;
    Cmd *cmd = state.cmd;
    char * line = state.line;
    char * origline = null;
    char * s;
    Ufunc * fun = current_fdef;
    int blank = 0;
    int ignore = 0;
    int i;
    int err = 0;

    if (fun == null) {
        fprintf(stderr, "gretl_function_append_line: fun is NULL\n");
        return 1;
    }

    blank = string_is_blank(line);


    /* note: avoid losing comment lines */
    origline = gretl_strdup(line);

    if (!blank) {
	s = line + strspn(line, " \t");
	state.line = s;    /* temporary! */
        err = get_command_index(state, GretlCmdIndex.func, 0);
	state.line = line; /* reset to original */
    }
    if (blank || err) {
        goto next_step;
    }

    /* carry out some basic structural checks */
    if (cmd.ci == GretlCmdIndex.quit) {
        gretl_errmsg_sprintf(_("%s: \"quit\" cannot be used in a function").ptr,
                             fun.name);
        err = GretlError.parse;
        cmd.ci = 0;
    } else if (cmd.flags & CmdFlags.endfun) {
        if (fun.n_lines == 0) {
            gretl_errmsg_sprintf(_("%s: empty function").ptr, fun.name);
            err = 1;
        }
        set_compiling_off();
    } else if (cmd.ci == GretlCmdIndex.set && changes_state(s + 3)) {
        fun.flags |= UfunAttrs.uses_set;
    } else if (cmd.ci == GretlCmdIndex.func) {
        err = GretlError.fnest;
    } else if (FLOW_CI(cmd.ci)) {
        if (cmd.ci == GretlCmdIndex.if_) {
            fun.flags |= UfunAttrs.has_flow;
            ifdepth++;
        } else if (ifdepth == 0) {
            gretl_errmsg_sprintf(_("%s: unbalanced if/else/endif").ptr, fun.name);
            err = GretlError.parse;
        } else if (cmd.ci == GretlCmdIndex.else_ && last_flow == GretlCmdIndex.else_) {
            gretl_errmsg_sprintf(_("%s: unbalanced if/else/endif").ptr, fun.name);
            err = GretlError.parse;
        } else if (cmd.ci == GretlCmdIndex.endif) {
            ifdepth--;
        }
        last_flow = cmd.ci;
    } else if (cmd.ci == GretlCmdIndex.loop) {
        fun.flags |= UfunAttrs.has_flow;
    } else if (cmd.ci < 0) {
        ignore = 1;
    }

  next_step:

    if (err) {
        set_compiling_off();
    }

    if (compiling) {
        if (blank) {
            /* just increment the index (placeholder) */
            fun.line_idx += 1;
        } else {
            /* actually add the line */
            i = fun.n_lines;
            err = push_function_line(fun, origline, cmd.ci, 1, null);
            if (!err) {
                origline = null; /* successfully donated */
                if (ignore) {
                    fun.lines[i].flags |= LnFlags.ignore;
                } else {
                    python_check(line);
                }
            }
        }
    } else {
        /* finished compilation */
        if (!err && ifdepth != 0) {
            gretl_errmsg_sprintf(_("%s: unbalanced if/else/endif").ptr, fun.name);
            err = GretlError.parse;
        }
        /* reset static var */
        ifdepth = 0;
        last_flow = 0;
    }

    if (!err && !compiling && (fun.flags & UfunAttrs.has_flow)) {
        ufunc_get_structure(fun);
    }

    free(origline);
    cmd.flags &= ~CmdFlags.endfun;

    if (err) {
        ufunc_unload(fun);
    }

    return err;
}

/* next block: handling function arguments */

/* Given a list of variables supplied as an argument to a function,
   copy the list under the name assigned by the function and
   make the variables referenced in that list accessible within
   the function. We also handle the case where the given arg
   was not a real list, but we can make a list out of it for the
   purposes of the function.
*/

/* Missing or "null" arg . gives an empty list, rather
   than no list. This is traditional behavior though it's
   (now) inconsistent with what happens with other types
   of argument.
*/

int add_empty_list(Fncall * call, FnParam * fp, DatasetRaw * dset)
{
    int err = 0;

    if (dset == null || dset.n == 0) {
        err = GretlError.nodata;
    } else {
        int * tmp = null;

        copy_list_as_arg(fp.name, tmp, &err);
    }

    return err;
}

void localize_list_members(Fncall * call, int * list, const(char) * lname, DatasetRaw * dset)
{
    int i;
    int vi;
    int level = fn_executing + 1;

    for (i=1; i<=list[0]; i++) {
        vi = list[i];
        if (vi > 0 && vi < dset.v) {
            if (!in_gretl_list(call.listvars, vi)) {
                gretl_list_append_term(&call.listvars, vi);
            }
            series_set_stack_level(dset, vi, level);
        }
    }
}

int localize_list(Fncall * call, FnArg * arg, FnParam * fp, DatasetRaw * dset)
{
    int * list = null;
    int err = 0;

    if (arg.type == GretlType.list) {
        /* actual list arg . copy to function level */
        list = arg.val.list;
        err = copy_as_arg(fp.name, GretlType.list, list);
        call.lists = g_list_prepend(call.lists, fp.name);
    } else if (arg.type == GretlType.useries) {
        /* series arg . becomes a singleton list */
        int[2] tmp = [1, arg.val.idnum];

        list = copy_list_as_arg(fp.name, tmp.ptr, &err);
        nullify_upname(arg);
    } else {
        /* "can't happen" */
        err = GretlError.data;
    }

    if (!err && list == null) {
        err = GretlError.alloc;
    }

    if (!err) {
        localize_list_members(call, list, fp.name, dset);
    }


    return err;
}

/* For the series with index number @ID, find the name of
   the list parameter of the current function call under
   which this series was made accessible to the function.
   We call this only if we know already that series @ID
   was in fact passed to the function via a list.
*/

extern (C) const(char) * series_get_list_parent(int ID)
{
    Fncall *call = current_function_call();
    const(char) * ret = null;

    if (call != null && call.lists != null) {
        GList * L = call.lists;
        int * list;

        while (L != null) {
            list = get_list_by_name(cast(const(char) *) L.data);
            if (list != null && in_gretl_list(list, ID)) {
                ret = cast(const(char) *) L.data;
                break;
            }
            L = L.next;
        }
    }

    return ret;
}

int localize_bundled_lists(Fncall * call, FnArg * arg, FnParam * fp, DatasetRaw * dset)
{
    GretlBundle * b = arg.val.b;
    GList *ll = gretl_bundle_get_lists(b);
    int err = 0;

    if (ll != null) {
        GList *lli = g_list_first(ll);

        while (lli != null) {
            localize_list_members(call, cast(int *) lli.data, null, dset);
            lli = g_list_next(lli);
        }
        g_list_free(ll);
    }

    return err;
}

void * arg_get_data(FnArg * arg, GretlType type)
{
    void * data = null;

    if (type == 0) {
        type = arg.type;
    }

    if (type == GretlType.matrix) {
        data = arg.val.m;
    } else if (type == GretlType.bundle) {
        data = arg.val.b;
    } else if (type == GretlType.string_) {
        data = arg.val.str;
    } else if (gretl_is_array_type(type)) {
        data = arg.val.a;
    }

    return data;
}

/* handle the case where the GUI passed an anonymous
   object in "pointerized" form */

int localize_object_as_shell(FnArg * arg, FnParam * fp)
{
    GretlType type = gretl_type_get_plain_type(arg.type);
    void *data = arg_get_data(arg, type);

    return arg_add_as_shell(fp.name, type, data);
}
