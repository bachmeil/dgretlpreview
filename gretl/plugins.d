/* Port of plugins.c (plugins.h has only prototypes: no header module)

   Plugins are shared objects loaded with dlopen. Note that the stock C plugins resolve their libgretl
   symbols against libgretl itself; they can only be loaded here if the executable provides those
   symbols (see log.md). HAVE_MPFR and WIN32 branches are dropped; HAVE_GMP is not defined, so asking
   for the GMP plugin fails as in a no-GMP build. */

module gretl.plugins;


import core.stdc.stdio : fprintf, stderr;
import core.stdc.string : strcat, strcpy;
version (Posix) import core.sys.posix.dlfcn : RTLD_LAZY, RTLD_NOW, dlclose, dlerror, dlopen, dlsym; else import gretl.winport : RTLD_LAZY, RTLD_NOW, dlclose, dlerror, dlopen, dlsym;

import gretl.ghashtable : GHashTable, g_hash_table_destroy, g_hash_table_insert, g_hash_table_lookup, g_hash_table_new,
    g_str_equal, g_str_hash;
import gretl.glib : gchar, gpointer;
import gretl.gmem : g_free;
import gretl.gretl_errors : gretl_errmsg_set, gretl_errmsg_sprintf;
import gretl.gretl_paths : gretl_plugin_path;
import gretl.gstrfuncs : g_strdup_printf;
import gretl.libgretlh;
import gretl.plugin.gnumeric_import : gnumeric_get_data;
import gretl.plugin.excel_import : xls_get_data;
import gretl.plugin.stata_export : stata_export;
import gretl.plugin.stata_import : dta_get_data;
import gretl.plugin.ods_import : ods_get_data;
import gretl.plugin.arma : arma_model, arma_select, gui_arma_select;
import gretl.plugin.garch : garch_model;
import gretl.plugin.johansen : johansen_boot_round, johansen_coint_test, johansen_estimate, trace_pvalue, vecm_test_restriction;
import gretl.plugin.stats_tables : dw_lookup, get_IPS_critvals, IPS_tbar_moments, IPS_tbar_rho_moments, qlr_asy_pvalue, qlr_critval_15_05, rank_sum_lookup, stock_yogo_lookup;
import gretl.plugin.sysest : system_estimate;
import gretl.plugin.urcdist : dfgls_pvalue, mackinnon_pvalue;
import gretl.plugin.xlsx_import : xlsx_get_data;

enum PluginCode {
    none,
    xls_import,
    xlsx_import,
    gnumeric_import,
    ods_import,
    johansen,
    vif,
    leverage,
    mp_ols,
    pca,
    progress_bar,
    range_mean,
    stats_tables,
    sysest,
    tramo_x12a,
    nistcheck,
    arma,
    arma_x12,
    garch,
    urcdist,
    kernel,
    fractal,
    poisson,
    mailer,
    eviews_import,
    stata_import,
    spss_import,
    sas_import,
    jmulti_import,
    zipfile,
    dpanel,
    heckit,
    odbc,
    quantreg,
    intreg,
    anova,
    duration,
    interpolate,
    iso3166,
    biprobit,
    reprobit,
    panurc,
    json_get,
    xml_get,
    stata_export,
    svm,
    regls,
    geoplot,
    purebin,
    bdstest,
    lpsolve,
    stepwise,
    kmeans
}

struct PluginInfo {
    int pnum;           /* index number of plugin */
    const(char) * pname;  /* name of plugin */
    void * handle;       /* handle obtained via dlopen or similar */
}

struct PluginFunctionInfo {
    const(char) * name;   /* name of function */
    int index;          /* index of the plugin that supplies it */
}

__gshared PluginInfo[53] plugins = [
    PluginInfo(0, null, null),
    PluginInfo(PluginCode.xls_import, "excel_import", null),
    PluginInfo(PluginCode.xlsx_import, "xlsx_import", null),
    PluginInfo(PluginCode.gnumeric_import, "gnumeric_import", null),
    PluginInfo(PluginCode.ods_import, "ods_import", null),
    PluginInfo(PluginCode.johansen, "johansen", null),
    PluginInfo(PluginCode.vif, "vif", null),
    PluginInfo(PluginCode.leverage, "leverage", null),
    PluginInfo(PluginCode.mp_ols, "mp_ols", null),
    PluginInfo(PluginCode.pca, "pca", null),
    PluginInfo(PluginCode.progress_bar, "progress_bar", null),
    PluginInfo(PluginCode.range_mean, "range-mean", null),
    PluginInfo(PluginCode.stats_tables, "stats_tables", null),
    PluginInfo(PluginCode.sysest, "sysest", null),
    PluginInfo(PluginCode.tramo_x12a, "tramo-x12a", null),
    PluginInfo(PluginCode.nistcheck, "nistcheck", null),
    PluginInfo(PluginCode.arma, "arma", null),
    PluginInfo(PluginCode.arma_x12, "arma_x12", null),
    PluginInfo(PluginCode.garch, "garch", null),
    PluginInfo(PluginCode.urcdist, "urcdist", null),
    PluginInfo(PluginCode.kernel, "kernel", null),
    PluginInfo(PluginCode.fractal, "fractals", null),
    PluginInfo(PluginCode.poisson, "poisson", null),
    PluginInfo(PluginCode.mailer, "mailer", null),
    PluginInfo(PluginCode.eviews_import, "eviews_import", null),
    PluginInfo(PluginCode.stata_import, "stata_import", null),
    PluginInfo(PluginCode.spss_import, "spss_import", null),
    PluginInfo(PluginCode.sas_import, "sas_import", null),
    PluginInfo(PluginCode.jmulti_import, "jmulti_import", null),
    PluginInfo(PluginCode.zipfile, "gretlzip", null),
    PluginInfo(PluginCode.dpanel, "dpanel", null),
    PluginInfo(PluginCode.heckit, "heckit", null),
    PluginInfo(PluginCode.odbc, "odbc_import", null),
    PluginInfo(PluginCode.quantreg, "quantreg", null),
    PluginInfo(PluginCode.intreg, "interval", null),
    PluginInfo(PluginCode.anova, "anova", null),
    PluginInfo(PluginCode.duration, "duration", null),
    PluginInfo(PluginCode.interpolate, "interpolate", null),
    PluginInfo(PluginCode.iso3166, "iso3166", null),
    PluginInfo(PluginCode.biprobit, "biprobit", null),
    PluginInfo(PluginCode.reprobit, "reprobit", null),
    PluginInfo(PluginCode.panurc, "panurc", null),
    PluginInfo(PluginCode.json_get, "json_get", null),
    PluginInfo(PluginCode.xml_get, "xml_get", null),
    PluginInfo(PluginCode.stata_export, "stata_export", null),
    PluginInfo(PluginCode.svm, "svm", null),
    PluginInfo(PluginCode.regls, "regls", null),
    PluginInfo(PluginCode.geoplot, "geoplot", null),
    PluginInfo(PluginCode.purebin, "purebin", null),
    PluginInfo(PluginCode.bdstest, "bdstest", null),
    PluginInfo(PluginCode.lpsolve, "lpsolve", null),
    PluginInfo(PluginCode.stepwise, "stepwise", null),
    PluginInfo(PluginCode.kmeans, "kmeans", null)
];

/* C also lists "mp_vector_ln" under HAVE_MPFR (not defined here) */
immutable PluginFunctionInfo[101] plugin_functions = [
    PluginFunctionInfo("xls_get_data", PluginCode.xls_import),
    PluginFunctionInfo("xlsx_get_data", PluginCode.xlsx_import),
    PluginFunctionInfo("gnumeric_get_data", PluginCode.gnumeric_import),
    PluginFunctionInfo("ods_get_data", PluginCode.ods_import),
    PluginFunctionInfo("wf1_get_data", PluginCode.eviews_import),
    PluginFunctionInfo("dta_get_data", PluginCode.stata_import),
    PluginFunctionInfo("sav_get_data", PluginCode.spss_import),
    PluginFunctionInfo("xport_get_data", PluginCode.sas_import),
    PluginFunctionInfo("jmulti_get_data", PluginCode.jmulti_import),
    PluginFunctionInfo("johansen_coint_test", PluginCode.johansen),
    PluginFunctionInfo("johansen_estimate", PluginCode.johansen),
    PluginFunctionInfo("johansen_boot_round", PluginCode.johansen),
    PluginFunctionInfo("vecm_test_restriction", PluginCode.johansen),
    PluginFunctionInfo("trace_pvalue", PluginCode.johansen),
    PluginFunctionInfo("model_leverage", PluginCode.leverage),
    PluginFunctionInfo("leverage_data_dialog", PluginCode.leverage),
    PluginFunctionInfo("compute_vifs", PluginCode.vif),
    PluginFunctionInfo("compute_bkw", PluginCode.vif),
    PluginFunctionInfo("bkw_matrix", PluginCode.vif),
    PluginFunctionInfo("mplsq", PluginCode.mp_ols),
    PluginFunctionInfo("matrix_mp_ols", PluginCode.mp_ols),
    PluginFunctionInfo("mp_vector_raise_to_power", PluginCode.mp_ols),
    PluginFunctionInfo("mp_bw_filter", PluginCode.mp_ols),
    PluginFunctionInfo("mp_midas_weights", PluginCode.mp_ols),
    PluginFunctionInfo("mp_midas_gradient", PluginCode.mp_ols),
    PluginFunctionInfo("pca_from_cmatrix", PluginCode.pca),
    PluginFunctionInfo("show_progress", PluginCode.progress_bar),
    PluginFunctionInfo("range_mean_graph", PluginCode.range_mean),
    PluginFunctionInfo("dw_lookup", PluginCode.stats_tables),
    PluginFunctionInfo("rank_sum_lookup", PluginCode.stats_tables),
    PluginFunctionInfo("stock_yogo_lookup", PluginCode.stats_tables),
    PluginFunctionInfo("get_IPS_critvals", PluginCode.stats_tables),
    PluginFunctionInfo("IPS_tbar_moments", PluginCode.stats_tables),
    PluginFunctionInfo("IPS_tbar_rho_moments", PluginCode.stats_tables),
    PluginFunctionInfo("qlr_asy_pvalue", PluginCode.stats_tables),
    PluginFunctionInfo("qlr_critval_15_05", PluginCode.stats_tables),
    PluginFunctionInfo("system_estimate", PluginCode.sysest),
    PluginFunctionInfo("write_tx_data", PluginCode.tramo_x12a),
    PluginFunctionInfo("exec_tx_script", PluginCode.tramo_x12a),
    PluginFunctionInfo("adjust_series", PluginCode.tramo_x12a),
    PluginFunctionInfo("linearize_series", PluginCode.tramo_x12a),
    PluginFunctionInfo("run_nist_tests", PluginCode.nistcheck),
    PluginFunctionInfo("arma_model", PluginCode.arma),
    PluginFunctionInfo("arma_select", PluginCode.arma),
    PluginFunctionInfo("gui_arma_select", PluginCode.arma),
    PluginFunctionInfo("arma_x12_model", PluginCode.arma_x12),
    PluginFunctionInfo("garch_model", PluginCode.garch),
    PluginFunctionInfo("count_data_estimate", PluginCode.poisson),
    PluginFunctionInfo("heckit_estimate", PluginCode.heckit),
    PluginFunctionInfo("interval_estimate", PluginCode.intreg),
    PluginFunctionInfo("tobit_via_intreg", PluginCode.intreg),
    PluginFunctionInfo("biprobit_estimate", PluginCode.biprobit),
    PluginFunctionInfo("biprobit_adjust_vcv", PluginCode.biprobit),
    PluginFunctionInfo("reprobit_estimate", PluginCode.reprobit),
    PluginFunctionInfo("stepwise_add", PluginCode.stepwise),
    PluginFunctionInfo("stepwise_omit", PluginCode.stepwise),
    PluginFunctionInfo("mackinnon_pvalue", PluginCode.urcdist),
    PluginFunctionInfo("dfgls_pvalue", PluginCode.urcdist),
    PluginFunctionInfo("kernel_density", PluginCode.kernel),
    PluginFunctionInfo("array_kernel_density", PluginCode.kernel),
    PluginFunctionInfo("kernel_density_matrix", PluginCode.kernel),
    PluginFunctionInfo("multiple_kd_matrix", PluginCode.kernel),
    PluginFunctionInfo("kmeans", PluginCode.kmeans),
    PluginFunctionInfo("hurst_exponent", PluginCode.fractal),
    PluginFunctionInfo("email_file", PluginCode.mailer),
    PluginFunctionInfo("gretl_native_make_zipfile", PluginCode.zipfile),
    PluginFunctionInfo("gretl_native_unzip", PluginCode.zipfile),
    PluginFunctionInfo("gretl_native_zip_datafile", PluginCode.zipfile),
    PluginFunctionInfo("dpd_estimate", PluginCode.dpanel),
    PluginFunctionInfo("gretl_odbc_check_dsn", PluginCode.odbc),
    PluginFunctionInfo("gretl_odbc_get_data", PluginCode.odbc),
    PluginFunctionInfo("rq_driver", PluginCode.quantreg),
    PluginFunctionInfo("lad_driver", PluginCode.quantreg),
    PluginFunctionInfo("gretl_anova", PluginCode.anova),
    PluginFunctionInfo("duration_estimate", PluginCode.duration),
    PluginFunctionInfo("time_disaggregate", PluginCode.interpolate),
    PluginFunctionInfo("iso_country", PluginCode.iso3166),
    PluginFunctionInfo("iso_country_array", PluginCode.iso3166),
    PluginFunctionInfo("iso_country_series", PluginCode.iso3166),
    PluginFunctionInfo("real_levin_lin", PluginCode.panurc),
    PluginFunctionInfo("json_get_string", PluginCode.json_get),
    PluginFunctionInfo("json_get_bundle", PluginCode.json_get),
    PluginFunctionInfo("json_bundle_get_terminals", PluginCode.json_get),
    PluginFunctionInfo("bundle_to_json", PluginCode.json_get),
    PluginFunctionInfo("xml_get", PluginCode.xml_get),
    PluginFunctionInfo("stata_export", PluginCode.stata_export),
    PluginFunctionInfo("gretl_svm_driver", PluginCode.svm),
    PluginFunctionInfo("gretl_regls", PluginCode.regls),
    PluginFunctionInfo("regls_xv_mpi", PluginCode.regls),
    PluginFunctionInfo("gretl_glasso", PluginCode.regls),
    PluginFunctionInfo("map_get_data", PluginCode.geoplot),
    PluginFunctionInfo("shp_get_bundle", PluginCode.geoplot),
    PluginFunctionInfo("geoplot", PluginCode.geoplot),
    PluginFunctionInfo("purebin_read_data", PluginCode.purebin),
    PluginFunctionInfo("purebin_write_data", PluginCode.purebin),
    PluginFunctionInfo("purebin_read_subset", PluginCode.purebin),
    PluginFunctionInfo("purebin_read_varnames", PluginCode.purebin),
    PluginFunctionInfo("bdstest", PluginCode.bdstest),
    PluginFunctionInfo("gretl_lpsolve", PluginCode.lpsolve),
    PluginFunctionInfo(null, 0)
];

private GHashTable * gretl_plugin_hash_init()
{
    GHashTable * ht;
    int i;

    ht = g_hash_table_new(&g_str_hash, &g_str_equal);

    /* Record the plugin index of each plugin function
       in a hash table under the key of the function
       name, permitting quick look-up.
    */
    for (i=0; plugin_functions[i].name != null; i++) {
        g_hash_table_insert(ht, cast(gpointer) plugin_functions[i].name,
                            cast(gpointer) cast(size_t) plugin_functions[i].index);
    }

    return ht;
}

private __gshared GHashTable * pht;

private int plugin_index_lookup(const(char) * name)
{
    gpointer ptr;

    if (name == null) {
        /* cleanup signal */
        if (pht != null) {
            g_hash_table_destroy(pht);
            pht = null;
        }
        return 0;
    }

    if (pht == null) {
        /* construct hash table if not already done */
        pht = gretl_plugin_hash_init();
    }

    ptr = g_hash_table_lookup(pht, name);

    return ptr == null ? 0 : cast(int) cast(size_t) ptr;
}

/**
 * gretl_dlopen:
 * @path: full path to the shared object to be opened.
 * @now: on *nix, if non-zero we call dlopen with the flag
 * RTLD_NOW, else we use RTLD_LAZY.
 *
 * Wrapper for opening a shared code object.
 *
 * Returns: handle to the shared object.
 */

void * gretl_dlopen(const(char) * path, int now)
{
    void * handle = null;

    handle = dlopen(path, (now)? RTLD_NOW : RTLD_LAZY);

    if (handle == null) {
        gretl_errmsg_sprintf(_("Failed to load plugin: %s").ptr, path);
        fprintf(stderr, "%s\n", dlerror());
    }

    return handle;
}

/**
 * gretl_dlsym:
 * @handle: handle to shared object; see gretl_dlopen().
 * @name: name of symbol to look up.
 *
 * Wrapper for obtaining a pointer to a named symbol in a shared
 * object represented by @handle.
 *
 * Returns: pointer corresponding to @name, or NULL.
 */

void * gretl_dlsym(void * handle, const(char) * name)
{
    return dlsym(handle, name);
}

private enum PLUGIN_EXT = ".so";

private void * get_plugin_handle_by_index(int i)
{
    void * handle = plugins[i].handle;

    if (handle == null) {
        /* not opened yet */
        char[MAXLEN] pluginpath = 0;

        strcpy(pluginpath.ptr, gretl_plugin_path());
        strcat(pluginpath.ptr, plugins[i].pname);
        strcat(pluginpath.ptr, PLUGIN_EXT.ptr);

        handle = gretl_dlopen(pluginpath.ptr, 0);
        /* store the pointer we got */
        plugins[i].handle = handle;
    }

    return handle;
}

private void * get_function_address(void * handle,
                                    const(char) * name)
{
    void * funp;

    funp = dlsym(handle, name);
    if (funp == null) {
        gchar * munged = g_strdup_printf("_%s", name);

        funp = dlsym(handle, munged);
        if (funp == null) {
            fprintf(stderr, "%s\n", dlerror());
        }
        g_free(munged);
    }

    return funp;
}

void plugins_cleanup()
{
    int i, n = cast(int) plugins.length;

    /* note: plugins[0] is a dummy entry */

    for (i=1; i<n; i++) {
        if (plugins[i].handle != null) {
            close_plugin(plugins[i].handle);
            plugins[i].handle = null;
        }
    }

    /* tear down the plugin look-up hash table */
    plugin_index_lookup(null);
}

/* Plugins that have been ported to D and are part of this program: looked up before any shared object
   is tried. (The C table above stays as it is, so that the remaining plugins still load by dlopen.)
   Add an entry here for each ported importer or other plugin function. */

private void * builtin_plugin_function(const(char) * funcname)
{
    import core.stdc.string : strcmp;

    if (strcmp(funcname, "xlsx_get_data") == 0) {
        return cast(void *) &xlsx_get_data;
    }
    if (strcmp(funcname, "gnumeric_get_data") == 0) {
        return cast(void *) &gnumeric_get_data;
    }
    if (strcmp(funcname, "ods_get_data") == 0) {
        return cast(void *) &ods_get_data;
    }
    if (strcmp(funcname, "xls_get_data") == 0) {
        return cast(void *) &xls_get_data;
    }
    if (strcmp(funcname, "stata_export") == 0) {
        return cast(void *) &stata_export;
    }
    if (strcmp(funcname, "dta_get_data") == 0) {
        return cast(void *) &dta_get_data;
    }
    if (strcmp(funcname, "system_estimate") == 0) {
        return cast(void *) &system_estimate;
    }
    if (strcmp(funcname, "trace_pvalue") == 0) {
        return cast(void *) &trace_pvalue;
    }
    if (strcmp(funcname, "johansen_coint_test") == 0) {
        return cast(void *) &johansen_coint_test;
    }
    if (strcmp(funcname, "johansen_estimate") == 0) {
        return cast(void *) &johansen_estimate;
    }
    if (strcmp(funcname, "johansen_boot_round") == 0) {
        return cast(void *) &johansen_boot_round;
    }
    if (strcmp(funcname, "vecm_test_restriction") == 0) {
        return cast(void *) &vecm_test_restriction;
    }
    if (strcmp(funcname, "arma_model") == 0) {
        return cast(void *) &arma_model;
    }
    if (strcmp(funcname, "arma_select") == 0) {
        return cast(void *) &arma_select;
    }
    if (strcmp(funcname, "gui_arma_select") == 0) {
        return cast(void *) &gui_arma_select;
    }
    if (strcmp(funcname, "garch_model") == 0) {
        return cast(void *) &garch_model;
    }
    if (strcmp(funcname, "dw_lookup") == 0) {
        return cast(void *) &dw_lookup;
    }
    if (strcmp(funcname, "rank_sum_lookup") == 0) {
        return cast(void *) &rank_sum_lookup;
    }
    if (strcmp(funcname, "stock_yogo_lookup") == 0) {
        return cast(void *) &stock_yogo_lookup;
    }
    if (strcmp(funcname, "qlr_asy_pvalue") == 0) {
        return cast(void *) &qlr_asy_pvalue;
    }
    if (strcmp(funcname, "qlr_critval_15_05") == 0) {
        return cast(void *) &qlr_critval_15_05;
    }
    if (strcmp(funcname, "get_IPS_critvals") == 0) {
        return cast(void *) &get_IPS_critvals;
    }
    if (strcmp(funcname, "IPS_tbar_moments") == 0) {
        return cast(void *) &IPS_tbar_moments;
    }
    if (strcmp(funcname, "IPS_tbar_rho_moments") == 0) {
        return cast(void *) &IPS_tbar_rho_moments;
    }
    if (strcmp(funcname, "mackinnon_pvalue") == 0) {
        return cast(void *) &mackinnon_pvalue;
    }
    if (strcmp(funcname, "dfgls_pvalue") == 0) {
        return cast(void *) &dfgls_pvalue;
    }
    return null;
}

/**
 * get_plugin_function:
 * @funcname: name of function to access.
 *
 * Looks up @funcname in gretl's internal plugin table and
 * attempts to open the plugin object file that offers the
 * given function.
 *
 * Returns: function pointer, or NULL on failure.
 */

void * get_plugin_function(const(char) * funcname)
{
    int i;
    void * funp = builtin_plugin_function(funcname);

    if (funp != null) {
        return funp;
    }

    i = plugin_index_lookup(funcname);

    /* C: "#if !HAVE_GMP"; HAVE_GMP is not defined here */
    if (i == PluginCode.mp_ols) {
        gretl_errmsg_set("GMP is not supported in this build");
        return null;
    }

    if (i > 0) {
        void * handle = get_plugin_handle_by_index(i);

        if (handle != null) {
            funp = get_function_address(handle, funcname);
        } else {
            fprintf(stderr, "%s: get_function_address failed\n",
                    funcname);
        }
    } else {
        fprintf(stderr, "%s: plugin_index_lookup failed\n",
                funcname);
    }

    if (funp == null) {
        gretl_errmsg_set(_("Couldn't load plugin function").ptr);
    }

    return funp;
}

/**
 * close_plugin:
 * @handle: pointer obtained via the handle argument to
 * get_plugin_function().
 *
 * Closes a shared plugin object. (C: KEEP_PLUGINS_OPEN is 0.)
 */

void close_plugin(void * handle)
{
    if (handle != null) {
        dlclose(handle);
    }
}
