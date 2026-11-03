/* Port of gretl_model.c */

module gretl.gretl_model;


import core.stdc.ctype : isdigit;
import core.stdc.math : exp, fabs, log, sin, cos, sqrt;
import core.stdc.stdio : fprintf, sprintf, sscanf, stderr;
import core.stdc.stdlib : atoi, calloc, free, malloc, qsort, realloc;
import core.stdc.string : memcpy, strcat, strcmp, strcpy, strlen, strncat, strrchr;
import std.algorithm.comparison : max;

import gretl.compareh : ModelAuxCode;
import gretl.dataseth : dataset_is_panel, dataset_is_time_series;
import gretl.dataset : destroy_dataset;
import gretl.estimate : makevcv;
import gretl.genmainh : ModelDataIndex;
import gretl.gmem;
import gretl.glib;
import gretl.gretl_commandsh;
import gretl.gretl_errorsh;
import gretl.gretl_listh : LISTSEP;
import gretl.gretl_matrix : gretl_invert_symmetric_matrix, gretl_matrix_alloc, gretl_matrix_copy,
    gretl_matrix_free, gretl_matrix_multiply_by_scalar, gretl_matrix_multiply_mod,
    gretl_matrix_print, gretl_matrix_qform, gretl_matrix_set_colnames,
    gretl_matrix_set_rownames, gretl_matrix_set_t1, gretl_matrix_set_t2,
    gretl_matrix_values, gretl_matrix_XTX_new, gretl_matrix_zero, gretl_zero_matrix_new;
import gretl.gretl_matrixh;
import gretl.gretl_modelh;
import gretl.gretl_prn : Prn, csv_format, plain_format, pprintf, pputc, pputs, prn_delim,
    rtf_format, tex_format;
import gretl.gstrfuncs;
import gretl.interacth : MAXSAVENAME;
import gretl.libgretlh;
import gretl.libset : libset_get_bool;
import gretl.libseth : SetKey;
import gretl.missingh;
import gretl.objstackh;
import gretl.pvalues : chisq_cdf_comp, normal_cdf_inverse, normal_pvalue_2,
    student_cdf_inverse, student_pvalue_2, tcrit95;
import gretl.xmlh;
import gretl.gretl_list : gretl_list_copy, gretl_list_copy_from_pos, gretl_list_new, gretl_list_separator_position, in_gretl_list;
import gretl.strutils : gretl_charsub, gretl_strdup, gretl_strdup_printf, gretl_utf8_truncate, gretl_xml_encode, gretl_xml_validate, strings_array_dup, strings_array_free, strings_array_new, strings_array_new_with_length;
import gretl.options : get_optval_string;
import gretl.gretl_errors : gretl_errmsg_set;
import gretl.missing : model_missing;
import gretl.gretl_utils : copyvec, gretl_compare_doubles, gretl_monotonic_time, ijton, record_test_result;
import gretl.qr_estimate : newey_west_OPG;
import gretl.gretl_xml : gretl_matrix_serialize, gretl_xml_get_cmplx_array, gretl_xml_get_double_array, gretl_xml_get_int_array, gretl_xml_get_list, gretl_xml_get_matrix, gretl_xml_get_prop_as_double, gretl_xml_get_prop_as_int, gretl_xml_get_prop_as_opt, gretl_xml_get_prop_as_string, gretl_xml_get_prop_as_uchar, gretl_xml_get_prop_as_uint64, gretl_xml_get_strings_array, gretl_xml_get_submask, gretl_xml_node_get_double, gretl_xml_node_get_int, gretl_xml_node_get_string, gretl_xml_node_get_trimmed_string, gretl_xml_put_double, gretl_xml_put_double_array, gretl_xml_put_string, gretl_xml_put_strings_array, gretl_xml_put_tagged_list;
import gretl.gretl_paths : gretl_remove;
import gretl.plugins : get_plugin_function;
import gretl.gretl_intl : gretl_pop_c_numeric_locale, gretl_push_c_numeric_locale;
import gretl.describe : free_vmatrix, vmatrix_new;
import gretl.genmain : current_series_index;
import gretl.matrix_extra : gretl_matrix_from_scalar, gretl_vcv_matrix_from_model;
import gretl.gretl_cmatrix : gretl_matrix_fft;
import gretl.subsample : copy_subsample_mask, free_subsample_mask, model_sample_problem, write_model_submask;
import gretl.compare : get_DW_pvalue_for_model;
import gretl.gretl_array : gretl_array_copy, gretl_array_destroy, gretl_array_new, gretl_array_serialize, gretl_array_set_string;
import gretl.gretl_bundle : gretl_bundle_donate_data, gretl_bundle_has_key, gretl_bundle_new, gretl_bundle_set_data, gretl_bundle_set_int, gretl_bundle_set_list, gretl_bundle_set_matrix, gretl_bundle_set_scalar, gretl_bundle_set_string;
import gretl.gretl_commands : gretl_command_number, gretl_command_word;

/* ---- declarations of not-yet-ported items used below ----
   Replace each with an import when its defining file is ported. */

extern (C) {
    Model * get_genr_model(GretlObjType * type);
    Model * get_last_model(GretlObjType * type);
    int gretl_function_depth();
    int gretl_model_protect(Model * pmod);
    void remove_model_from_stack_on_exit(Model * pmod);
}

/* ---- end of declarations ---- */

alias DESTFUNC = void function(void *);

struct ModelDataItem {
    char * key;
    void * ptr;
    int type;
    size_t size;
    DESTFUNC destructor;
}

struct ModelTest {
    int type;
    int order;
    char * param;
    ubyte teststat;
    int dfn;
    double dfd;
    double value;
    double pvalue;
    double crit;
    double alpha;
    Gretlopt opt;
}

private enum CSLEN = 64;

private struct CoeffSep {
    char[CSLEN] str;
    int pos;
}

private enum PNAMELEN = 16; /* for parameter names */

private void gretl_test_init(ModelTest * test, ModelTestType ttype)
{
    test.type = ttype;
    test.order = 0;
    test.param = null;
    test.teststat = cast(ubyte) TestStats.none;
    test.dfn = 0;
    test.dfd = 0;
    test.value = test.pvalue = NADBL;
    test.crit = test.alpha = NADBL;
    test.opt = Gretlopt.OPT_NONE;
}

private const(char) * test_type_key(ModelTestType t)
{
    if (t == ModelTestType.add) {
        return "add_test";
    } else if (t == ModelTestType.arch) {
        return "arch_test";
    } else if (t == ModelTestType.autocorr) {
        return "autocorr_test";
    } else if (t == ModelTestType.chow ||
               t == ModelTestType.chowdum) {
        return "chow_test";
    } else if (t == ModelTestType.cusum) {
        return "cusum_test";
    } else if (t == ModelTestType.qlr) {
        return "qlr_test";
    } else if (t == ModelTestType.groupwise) {
        return "grpwise_test";
    } else if (t == ModelTestType.logs) {
        return "logs_test";
    } else if (t == ModelTestType.normal) {
        return "normality_test";
    } else if (t == ModelTestType.omit) {
        return "omit_test";
    } else if (t == ModelTestType.reset) {
        return "reset_test";
    } else if (t == ModelTestType.squares) {
        return "squares_test";
    } else if (t == ModelTestType.whites) {
        return "whites_test";
    } else if (t == ModelTestType.sargan) {
        return "sargan_test";
    } else if (t == ModelTestType.iv_hausman ||
               t == ModelTestType.panel_hausman) {
        return "hausman_test";
    } else if (t == ModelTestType.panel_f ||
               t == ModelTestType.panel_welch) {
        return "fixed_effects_F";
    } else if (t == ModelTestType.panel_bp ||
               t == ModelTestType.bp) {
        return "bp_test";
    } else if (t == ModelTestType.panel_timedum) {
        return "timedum_test";
    } else if (t == ModelTestType.het_1) {
        return "het1_test";
    } else if (t == ModelTestType.comfac) {
        return "comfac_test";
    } else if (t == ModelTestType.indep) {
        return "independence_test";
    } else if (t == ModelTestType.re) {
        return "rho_test";
    } else if (t == ModelTestType.within_f) {
        return "within_F";
    } else if (t == ModelTestType.re_wald) {
        return "re_wald_test";
    } else if (t == ModelTestType.xdepend) {
        return "cross_sectional_dependence_test";
    } else if (t == ModelTestType.panel_ar) {
        return "panel_ar_test";
    } else if (t == ModelTestType.bds) {
        return "BDS_test";
    } else {
        fprintf(stderr, "test_type_key(): type %d has no key!\n", t);
        return null;
    }
}

private void free_model_data_item(ModelDataItem * item)
{
    if (item.destructor != null) {
        (*item.destructor)(item.ptr);
    } else {
        free(item.ptr);
    }
    free(item.key);
    free(item);
}

private ModelDataItem * create_data_item(const(char) * key, void * ptr,
                                         GretlType type, size_t size,
                                         DESTFUNC destructor)
{
    ModelDataItem * item = cast(ModelDataItem *) malloc(ModelDataItem.sizeof);

    if (item != null) {
        item.key = gretl_strdup(key);
        if (item.key == null) {
            free(item);
            item = null;
        } else {
            item.ptr = ptr;
            item.type = type;
            item.size = size;
            item.destructor = destructor;
        }
    }

    return item;
}

private ModelDataItem * replicate_data_item(const(ModelDataItem) * orig)
{
    ModelDataItem * item = cast(ModelDataItem *) malloc(ModelDataItem.sizeof);
    int err = 0;

    if (item != null) {
        item.key = gretl_strdup(orig.key);
        if (item.key == null) {
            free(item);
            item = null;
        }
    }

    if (item != null) {
        if (orig.type == GretlType.matrix) {
            item.ptr = gretl_matrix_copy(cast(const(GretlMatrixRaw) *) orig.ptr);
        } else if (orig.type == GretlType.array) {
            item.ptr = gretl_array_copy(cast(GretlArray *) orig.ptr, &err);
        } else {
            item.ptr = malloc(orig.size);
        }
        if (item.ptr == null) {
            free(item.key);
            free(item);
            item = null;
        }
    }

    if (item != null) {
        if (orig.type != GretlType.matrix &&
            orig.type != GretlType.array) {
            memcpy(item.ptr, orig.ptr, orig.size);
        }
        item.type = orig.type;
        item.size = orig.size;
        item.destructor = orig.destructor;
    }

    return item;
}

/**
 * gretl_model_set_data_with_destructor:
 * @pmod: pointer to #MODEL.
 * @key: key string for data, used in retrieval.
 * @ptr: data-pointer to be attached to model.
 * @type: type of data to set.
 * @size: size of data in bytes.
 * @destructor: pointer to function that should be used to free
 * the data-pointer in question.
 *
 * Attaches data to @pmod: the data can be retrieved later using
 * gretl_model_get_data().  Note that the data are not "physically"
 * copied to the model; simply, @ptr is recorded on the model.
 * This means that the data referenced by the pointer now in
 * effect belong to @pmod.  When @pmod is cleared with clear_model(),
 * @destructor will be invoked with @ptr as its single argument.
 * If a simple "free" is OK for freeing the data, you can use
 * gretl_model_set_data() instead.
 *
 * The @size is needed in case the model is copied with
 * copy_model(), in which case the target of the copying
 * operation receives a newly allocated copy of the data in
 * question.
 *
 * Returns: 0 on success, 1 on failure.
 */

private int gretl_model_set_data_with_destructor(Model * pmod,
                                                 const(char) * key,
                                                 void * ptr,
                                                 GretlType type,
                                                 size_t size,
                                                 DESTFUNC destructor)
{
    ModelDataItem ** items;
    ModelDataItem * item;
    int i, n;

    for (i=0; i<pmod.n_data_items; i++) {
        item = pmod.data_items[i];
        if (!strcmp(key, item.key)) {
            /* there's a pre-existing item of this name */
            if (item.destructor != null) {
                (*item.destructor)(item.ptr);
            } else {
                free(item.ptr);
            }
            item.type = type;
            item.ptr = ptr;
            item.size = size;
            item.destructor = destructor;
            /* handled */
            return 0;
        }
    }

    n = pmod.n_data_items + 1;

    items = cast(ModelDataItem **) realloc(pmod.data_items, n * (ModelDataItem *).sizeof);
    if (items == null) {
        return 1;
    }

    pmod.data_items = items;

    item = create_data_item(key, ptr, type, size, destructor);
    if (item == null) {
        return 1;
    }

    pmod.data_items[n - 1] = item;
    pmod.n_data_items += 1;

    return 0;
}

/**
 * gretl_model_set_data:
 * @pmod: pointer to #MODEL.
 * @key: key string for data, used in retrieval.
 * @ptr: data-pointer to be attached to model.
 * @type: type of the data to set.
 * @size: size of data in bytes.
 *
 * Attaches data to @pmod: the data can be retrieved later using
 * gretl_model_get_data().  Note that the data are not "physically"
 * copied to the model; simply, @ptr is recorded on the model,
 * meaning that @pmod takes ownership of the data. The data
 * pointer will be freed when @pmod is cleared with clear_model().
 * If the data item has structure that requires special treatment
 * on freeing, use a type-specific function instead, for example
 * gretl_model_set_matrix_as_data().
 *
 * The @size is needed just in case the model is copied with
 * copy_model(), in which case the target of the copying
 * operation receives a newly allocated copy of the data in
 * question.
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_model_set_data(Model * pmod, const(char) * key, void * ptr,
                         GretlType type, size_t size)
{
    return gretl_model_set_data_with_destructor(pmod, key, ptr, type,
                                                size, null);
}

private void matrix_free_callback(void * p)
{
    gretl_matrix_free(cast(GretlMatrixRaw *) p);
}

/**
 * gretl_model_set_matrix_as_data:
 * @pmod: pointer to #MODEL.
 * @key: key string, used in retrieval.
 * @m: matrix to attach.
 *
 * Attaches @m to @pmod as data, recoverable via the key @key
 * using gretl_model_get_data().
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_model_set_matrix_as_data(Model * pmod, const(char) * key,
                                   GretlMatrixRaw * m)
{
    return gretl_model_set_data_with_destructor(pmod, key, cast(void *) m,
                                                GretlType.matrix, 0,
                                                &matrix_free_callback);
}

private void array_free_callback(void * p)
{
    gretl_array_destroy(cast(GretlArray *) p);
}

/**
 * gretl_model_set_array_as_data:
 * @pmod: pointer to #MODEL.
 * @key: key string, used in retrieval.
 * @A: array to attach.
 *
 * Attaches @A to @pmod as data, recoverable via the key @key
 * using gretl_model_get_data().
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_model_set_array_as_data(Model * pmod, const(char) * key,
                                  GretlArray * A)
{
    return gretl_model_set_data_with_destructor(pmod, key, cast(void *) A,
                                                GretlType.array, 0,
                                                &array_free_callback);
}

/**
 * gretl_model_set_list_as_data:
 * @pmod: pointer to #MODEL.
 * @key: key string, used in retrieval.
 * @list: list to attach.
 *
 * Attaches @list to @pmod as data, recoverable via the key @key
 * using gretl_model_get_list().  Note that the model takes
 * ownership of the supplied list.
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_model_set_list_as_data(Model * pmod, const(char) * key, int * list)
{
    size_t size = (list[0] + 1) * int.sizeof;

    return gretl_model_set_data_with_destructor(pmod, key, cast(void *) list,
                                                GretlType.list, size,
                                                null);
}

/**
 * gretl_model_set_string_as_data:
 * @pmod: pointer to #MODEL.
 * @key: key string, used in retrieval.
 * @str: string to attach.
 *
 * Attaches @str to @pmod as data, recoverable via the key @key
 * using gretl_model_get_data().
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_model_set_string_as_data(Model * pmod, const(char) * key, char * str)
{
    size_t size = strlen(str) + 1;

    return gretl_model_set_data_with_destructor(pmod, key, cast(void *) str,
                                                GretlType.string_, size,
                                                null);
}

/**
 * gretl_model_set_int:
 * @pmod: pointer to #MODEL.
 * @key: key string, used in retrieval.
 * @val: integer value to set.
 *
 * Records an integer value on a model: the value can be retrieved
 * later using gretl_model_get_int(), using the appropriate @key.
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_model_set_int(Model * pmod, const(char) * key, int val)
{
    int * valp;
    int err = 0;

    /* if value is already set, reset it */
    valp = cast(int *) gretl_model_get_data(pmod, key);
    if (valp != null) {
        *valp = val;
        return 0;
    }

    /* otherwise start from scratch */
    valp = cast(int *) malloc(int.sizeof);
    if (valp == null) {
        err = GretlError.alloc;
    } else {
        *valp = val;
        err = gretl_model_set_data(pmod, key, valp, GretlType.int_,
                                   int.sizeof);
        if (err) {
            free(valp);
        }
    }

    return err;
}

/**
 * gretl_model_set_double:
 * @pmod: pointer to model.
 * @key: key string, used in retrieval.
 * @val: double-precision value to set.
 *
 * Records a floating-point value on @pmod: the value can be
 * retrieved later using gretl_model_get_double() with the
 * appropriate @key.
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_model_set_double(Model * pmod, const(char) * key, double val)
{
    double * valp;
    int err;

    /* if value is already set, reset it */
    valp = cast(double *) gretl_model_get_data(pmod, key);
    if (valp != null) {
        *valp = val;
        return 0;
    }

    valp = cast(double *) malloc(double.sizeof);
    if (valp == null) return 1;

    *valp = val;

    err = gretl_model_set_data(pmod, key, valp, GretlType.double_,
                               double.sizeof);
    if (err) free(valp);

    return err;
}

private VCVInfo * vcv_info_new()
{
    VCVInfo * vi;

    vi = cast(VCVInfo *) malloc(VCVInfo.sizeof);

    if (vi != null) {
        vi.vmaj = vi.vmin = 0;
        vi.order = 0;
        vi.flags = cast(VCVFlags) 0;
        vi.bw = NADBL;
        vi.cv1 = null;
        vi.cv2 = null;
    }

    return vi;
}

private void vcv_info_free(void * data)
{
    VCVInfo * vi = cast(VCVInfo *) data;

    if (vi != null) {
        free(vi.cv1);
        free(vi.cv2);
        free(vi);
    }
}

/**
 * gretl_model_set_full_vcv_info:
 *
 * Returns: 0 on success, 1 on failure.
 */

private int gretl_model_set_full_vcv_info(Model * pmod, int vmaj, int vmin,
                                          int order, int flags, double bw,
                                          const(char) * cv1, const(char) * cv2)
{
    VCVInfo * vi;
    int prev = 0;
    int err = 0;

    vi = cast(VCVInfo *) gretl_model_get_data(pmod, "vcv_info");

    if (vi == null) {
        vi = vcv_info_new();
        if (vi == null) {
            return GretlError.alloc;
        }
    } else {
        prev = 1;
        free(vi.cv1);
        free(vi.cv2);
        vi.cv1 = vi.cv2 = null;
    }

    vi.vmaj = vmaj;
    vi.vmin = vmin;
    vi.order = order;
    vi.flags = cast(VCVFlags) flags;
    vi.bw = bw;

    if (cv1 != null) {
        vi.cv1 = gretl_strdup(cv1);
    }
    if (cv2 != null) {
        vi.cv2 = gretl_strdup(cv2);
    }

    if (!prev) {
        err = gretl_model_set_data_with_destructor(pmod, "vcv_info", vi,
                                                   GretlType.struct_,
                                                   VCVInfo.sizeof,
                                                   &vcv_info_free);
    }

    return err;
}

/**
 * gretl_model_set_vcv_info:
 * @pmod: pointer to model.
 * @vmaj: top-level VCV type.
 * @vmin: variant under @vmaj, if applicable.
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_model_set_vcv_info(Model * pmod, int vmaj, int vmin)
{
    return gretl_model_set_full_vcv_info(pmod, vmaj, vmin,
                                         0, 0, 0, null, null);
}

/**
 * gretl_model_set_cluster_vcv_info:
 * @pmod: pointer to model.
 * @cv1: name of (first) cluster variable.
 * @cv2: name of second cluster variable or NULL.
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_model_set_cluster_vcv_info(Model * pmod,
                                     const(char) * cv1,
                                     const(char) * cv2)
{
    return gretl_model_set_full_vcv_info(pmod, VCVMajorType.cluster,
                                         0, 0, 0, 0,
                                         cv1, cv2);
}

/**
 * gretl_model_set_hac_vcv_info:
 * @pmod: pointer to model.
 * @kern: kernel type.
 * @order: lag order.
 * @flags: bitflags.
 * @bw: QS bandwidth, or #NADBL if not applicable.
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_model_set_hac_vcv_info(Model * pmod, int kern,
                                 int order, int flags,
                                 double bw)
{
    return gretl_model_set_full_vcv_info(pmod, VCVMajorType.hac, kern,
                                         order, flags, bw,
                                         null, null);
}

/**
 * gretl_model_set_hac_order:
 * @pmod: pointer to model.
 * @order: lag order.
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_model_set_hac_order(Model * pmod, int order)
{
    VCVInfo * vi = cast(VCVInfo *) gretl_model_get_data(pmod, "vcv_info");

    if (vi != null) {
        vi.order = order;
        return 0;
    } else {
        return GretlError.data;
    }
}

/**
 * gretl_model_get_vcv_type:
 * @pmod: pointer to model.
 *
 * Returns: index of variance-covariance matrix type.
 */

int gretl_model_get_vcv_type(const(Model) * pmod)
{
    VCVInfo * vi = cast(VCVInfo *) gretl_model_get_data(pmod, "vcv_info");

    if (vi != null) {
        return vi.vmaj;
    } else {
        return 0;
    }
}

/**
 * gretl_model_get_hc_version:
 * @pmod: pointer to model.
 *
 * Returns: the "HC" (HCCME) variant employed in @pmod,
 * or -1 if this does not apply.
 */

int gretl_model_get_hc_version(const(Model) * pmod)
{
    VCVInfo * vi = cast(VCVInfo *) gretl_model_get_data(pmod, "vcv_info");

    if (vi != null && vi.vmaj == VCVMajorType.hc) {
        return vi.vmin;
    } else {
        return -1;
    }
}

/**
 * gretl_model_get_cluster_vname:
 * @pmod: pointer to model.
 *
 * Returns: the name of the (first) clustering variable used for
 * the variance-covariance matrix in @pmod, or NULL if there is no
 * such variable.
 */

const(char) * gretl_model_get_cluster_vname(const(Model) * pmod)
{
    VCVInfo * vi = cast(VCVInfo *) gretl_model_get_data(pmod, "vcv_info");

    if (vi != null && vi.vmaj == VCVMajorType.cluster) {
        return vi.cv1;
    } else {
        return null;
    }
}

/**
 * gretl_model_get_cluster_vname2:
 * @pmod: pointer to model.
 *
 * Returns: the name of the second clustering variable used for
 * the variance-covariance matrix in @pmod, or NULL if there is no
 * such variable.
 */

const(char) * gretl_model_get_cluster_vname2(const(Model) * pmod)
{
    VCVInfo * vi = cast(VCVInfo *) gretl_model_get_data(pmod, "vcv_info");

    if (vi != null && vi.vmaj == VCVMajorType.cluster) {
        return vi.cv2;
    } else {
        return null;
    }
}

/**
 * gretl_model_get_data_full:
 * @pmod: pointer to model.
 * @key: key string.
 * @copied: location to receive flag indicating whether the
 * return value is an allocated copy of the original data.
 * @type: location to receive data type.
 * @sz: location to receive the size of the data.
 *
 * Returns: the data pointer identified by @key, or %NULL on failure.
 * If a non-zero value is written to @copied this indicates that the
 * return value is a copy of the original (and therefore it is the
 * caller's responsibility to free the data when it is no longer
 * required).
 */

void * gretl_model_get_data_full(const(Model) * pmod, const(char) * key,
                                 GretlType * type, int * copied,
                                 size_t * sz)
{
    void * ret = null;
    GretlType itype = cast(GretlType) 0;
    size_t isize = 0;
    int alloced = 0;
    int i, found = 0;

    if (pmod == null) {
        return null;
    }

    for (i=0; i<pmod.n_data_items; i++) {
        if (!strcmp(key, pmod.data_items[i].key)) {
            ret = cast(void *) pmod.data_items[i].ptr;
            itype = cast(GretlType) pmod.data_items[i].type;
            isize = pmod.data_items[i].size;
            found = 1;
            break;
        }
    }

    if (!found && pmod.tests != null) {
        const(ModelTest) * test;
        const(char) * tkey;
        GretlBundle * b;

        for (i=0; i<pmod.ntests; i++) {
            test = &pmod.tests[i];
            tkey = test_type_key(cast(ModelTestType) test.type);
            if (tkey != null && !strcmp(key, tkey)) {
                b = bundlize_test(test);
                if (b != null) {
                    ret = b;
                    itype = GretlType.bundle;
                    alloced = 1;
                }
                break;
            }
        }
    }

    if (ret != null) {
        if (type != null) {
            *type = itype;
        }
        if (sz != null) {
            *sz = isize;
        }
        if (copied != null) {
            *copied = alloced;
        }
    }

    return ret;
}

/**
 * gretl_model_get_data:
 * @pmod: pointer to model.
 * @key: key string.
 *
 * Returns: the data pointer identified by @key, or %NULL on failure.
 */

void * gretl_model_get_data(const(Model) * pmod, const(char) * key)
{
    return gretl_model_get_data_full(pmod, key, null, null, null);
}

/**
 * gretl_model_steal_data:
 * @pmod: pointer to model.
 * @key: key string.
 *
 * Returns: the data pointer identified by @key, or %NULL on failure.
 * The caller takes ownership.
 */

void * gretl_model_steal_data(Model * pmod, const(char) * key)
{
    void * ret = null;
    int i;

    if (pmod == null || key == null) {
        return null;
    }

    for (i=0; i<pmod.n_data_items; i++) {
        if (!strcmp(key, pmod.data_items[i].key)) {
            ret = pmod.data_items[i].ptr;
            discard_model_data_item(pmod, key, 0);
            break;
        }
    }

    return ret;
}

/**
 * gretl_model_get_int:
 * @pmod: pointer to model.
 * @key: key string.
 *
 * Returns: the integer value identified by @key, or 0 on failure.
 */

int gretl_model_get_int(const(Model) * pmod, const(char) * key)
{
    int * valp = null;
    int i;

    if (pmod == null) {
        return 0;
    }

    for (i=0; i<pmod.n_data_items; i++) {
        if (pmod.data_items[i].type != GretlType.int_) {
            continue;
        }
        if (!strcmp(key, pmod.data_items[i].key)) {
            valp = cast(int *) pmod.data_items[i].ptr;
            return *valp;
        }
    }

    return 0;
}

/**
 * gretl_model_get_double:
 * @pmod: pointer to model.
 * @key: key string.
 *
 * Returns: the double-precision value identified by @key, or
 * #NADBL on failure.
 */

double gretl_model_get_double(const(Model) * pmod, const(char) * key)
{
    double * valp = null;
    int i;

    if (pmod == null) {
        return NADBL;
    }

    for (i=0; i<pmod.n_data_items; i++) {
        if (pmod.data_items[i].type != GretlType.double_) {
            continue;
        }
        if (!strcmp(key, pmod.data_items[i].key)) {
            valp = cast(double *) pmod.data_items[i].ptr;
            return *valp;
        }
    }

    return NADBL;
}

/**
 * gretl_model_get_double_default:
 * @pmod: pointer to model.
 * @key: key string.
 * @deflt: default value
 *
 * Returns: the double-precision value identified by @key, or
 * @deflt if there is no such value.
 */

double gretl_model_get_double_default(const(Model) * pmod,
                                      const(char) * key,
                                      double deflt)
{
    double * valp = null;
    double ret = deflt;
    int i;

    if (pmod == null) {
        return NADBL;
    }

    for (i=0; i<pmod.n_data_items; i++) {
        if (pmod.data_items[i].type != GretlType.double_) {
            continue;
        }
        if (!strcmp(key, pmod.data_items[i].key)) {
            valp = cast(double *) pmod.data_items[i].ptr;
            ret = *valp;
            break;
        }
    }

    return ret;
}

/**
 * gretl_model_get_list:
 * @pmod: pointer to model.
 * @key: key string.
 *
 * Returns: the list of integers identified by @key, or
 * %NULL on failure.
 */

int * gretl_model_get_list(const(Model) * pmod, const(char) * key)
{
    int * list = null;
    int i;

    if (pmod == null) {
        return null;
    }

    if (!strcmp(key, "xlist") && pmod.xlist != null) {
        return cast(int *) pmod.xlist;
    }

    for (i=0; i<pmod.n_data_items; i++) {
        if (pmod.data_items[i].type == GretlType.list &&
            !strcmp(key, pmod.data_items[i].key)) {
            list = cast(int *) pmod.data_items[i].ptr;
            break;
        }
    }

    return list;
}

/**
 * gretl_model_set_coeff_separator:
 * @pmod: pointer to model.
 * @s: informative string (or %NULL).
 * @pos: position in the array of coefficients.
 *
 * Arranges for the insertion of the given string (or a blank line if
 * @s is %NULL) at the given position in the array of coefficients,
 * when the model is printed.  The extra line is printed before
 * coefficient @pos, where @pos is interpreted as a zero-based
 * index.
 *
 * Returns: 0 on success, %E_ALLOC on failure.
 */

int gretl_model_set_coeff_separator(Model * pmod, const(char) * s, int pos)
{
    CoeffSep * cs = cast(CoeffSep *) malloc(CoeffSep.sizeof);
    int err = 0;

    if (cs == null) {
        return GretlError.alloc;
    }

    cs.str[0] = '\0';
    if (s != null) {
        strncat(cs.str.ptr, s, CSLEN - 1);
    }
    cs.pos = pos;

    err = gretl_model_set_data(pmod, "coeffsep", cs, GretlType.struct_,
                               CoeffSep.sizeof);
    if (err) {
        free(cs);
    }

    return err;
}

/**
 * gretl_model_get_coeff_separator:
 * @pmod: pointer to model.
 * @ps: location to receive string, if any, or NULL.
 * @ppos: location to receive position, if any, or NULL.
 *
 * Retrieves information that has been set on @pmod regarding the
 * insertion of an extra line when printing the coefficients, if any.
 *
 * Returns: 1 if such information is present, 0 otherwise.
 */

int gretl_model_get_coeff_separator(const(Model) * pmod, char ** ps, int * ppos)
{
    CoeffSep * cs = cast(CoeffSep *) gretl_model_get_data(pmod, "coeffsep");

    if (cs == null) {
        return 0;
    }

    if (ps != null) {
        *ps = gretl_strdup(cs.str.ptr);
    }
    if (ppos != null) {
        *ppos = cs.pos;
    }

    return 1;
}

private void make_cname(const(char) * orig, char * cname)
{
    const(char) * p;

    if (orig == null || *orig == 0) {
        return;
    }

    p = strrchr(orig, '_');

    if (p == null) {
        strcpy(cname, orig);
    } else {
        ubyte c = cast(ubyte) *(p + 1);

        if (isdigit(c)) {
            int lag = atoi(++p);

            sprintf(cname, "ut^2(-%d)", lag);
        }
    }
}

private void plain_ar_coeff_name(char * targ, const(Model) * pmod, int i)
{
    int k = i - pmod.ncoeff;

    if (pmod.arinfo != null && pmod.arinfo.arlist != null &&
        k >= 0 && k < pmod.arinfo.arlist[0]) {
        sprintf(targ, "u_%d", pmod.arinfo.arlist[k+1]);
    } else {
        strcpy(targ, "unknown");
    }
}

/**
 * gretl_model_get_param_name:
 * @pmod: pointer to model.
 * @dset: dataset information.
 * @i: index number for parameter, zero-based, corresponding
 * to position in the %coeff array in @pmod.
 * @targ: string into which to write param name.
 *
 * Writes the appropriate parameter name into @targ, which
 * should be at least #VNAMELEN bytes long.  Usually this is
 * the name of a variable in the dataset, but sometimes it is
 * a special string (e.g. for nonlinear models).
 *
 * Returns: @targ.
 */

char * gretl_model_get_param_name(const(Model) * pmod,
                                  const(DatasetRaw) * dset,
                                  int i, char * targ)
{
    *targ = '\0';

    if (pmod != null) {
        int j = i + 2;
        int k = -1;

        if (pmod.aux == ModelAuxCode.arch) {
            make_cname(dset.varname[pmod.list[j]], targ);
        } else if (pmod.ci == GretlCmdIndex.panel && (pmod.opt & Gretlopt.OPT_H)) {
            /* --unit-weights */
            strcpy(targ, dset.varname[pmod.list[j]]);
        } else if (NONLIST_MODEL(cast(GretlCmdIndex) pmod.ci) ||
                   pmod.ci == GretlCmdIndex.arma || pmod.ci == GretlCmdIndex.panel ||
                   pmod.ci == GretlCmdIndex.dpanel || pmod.ci == GretlCmdIndex.garch ||
                   pmod.ci == GretlCmdIndex.biprobit) {
            k = i;
        } else if (pmod.ci == GretlCmdIndex.mpols && pmod.params != null) {
            k = i;
        } else if ((pmod.ci == GretlCmdIndex.probit || pmod.ci == GretlCmdIndex.logit ||
                    pmod.ci == GretlCmdIndex.heckit) && pmod.params != null) {
            k = i;
        } else if (pmod.ci == GretlCmdIndex.ar && i >= pmod.ncoeff) {
            plain_ar_coeff_name(targ, pmod, i);
        } else if (pmod.ci == GretlCmdIndex.arch && i >= pmod.ncoeff) {
            sprintf(targ, "alpha(%d)", i - pmod.ncoeff);
        } else if (pmod.list == null || j > pmod.list[0] ||
                   pmod.list[j] >= dset.v) {
            k = i;
        } else {
            strcpy(targ, dset.varname[pmod.list[j]]);
        }

        if (k >= 0) {
            if (pmod.params != null) {
                if (g_strlcpy(targ, pmod.params[k], VNAMELEN) >= VNAMELEN) {
                    fprintf(stderr, "gretl_model_get_param_name: "
                            ~ "name too long, truncated\n");
                }
            } else {
                strcpy(targ, "unknown");
            }
        }
    }

    return targ;
}

GretlArray * gretl_model_get_param_names(const(Model) * pmod,
                                         const(DatasetRaw) * dset,
                                         int * err)
{
    GretlArray * names = null;

    if (pmod == null) {
        return null;
    }

    names = gretl_array_new(GretlType.strings,
                            pmod.ncoeff, err);

    if (names != null) {
        char[VNAMELEN] targ;
        int i;

        for (i=0; i<pmod.ncoeff; i++) {
            gretl_model_get_param_name(pmod, dset, i, targ.ptr);
            gretl_array_set_string(names, i, targ.ptr, 1);
        }
    }

    return names;
}

/**
 * gretl_model_get_param_number:
 * @pmod: pointer to model.
 * @dset: dataset information.
 * @s: name of model parameter.
 *
 * Returns: the zero-based index of the coefficient in @pmod
 * corresponding to @s, or -1 if @s is not the name
 * of a parameter.
 */

int gretl_model_get_param_number(const(Model) * pmod, const(DatasetRaw) * dset,
                                 const(char) * s)
{
    char[VNAMELEN] pname;
    char[VNAMELEN] tmp;
    int i, idx = -1;

    if (pmod == null || s == null) {
        return -1;
    }

    if (!strcmp(s, "0")) {
        strcpy(pname.ptr, "const");
    } else {
        strcpy(pname.ptr, s);
    }

    for (i=0; i<pmod.ncoeff; i++) {
        gretl_model_get_param_name(pmod, dset, i, tmp.ptr);
        if (!strcmp(pname.ptr, tmp.ptr)) {
            idx = i;
            break;
        }
    }

    return idx;
}

void free_coeff_intervals(CoeffIntervals * cf)
{
    free(cf.coeff);
    free(cf.maxerr);
    if (cf.names != null) {
        strings_array_free(cf.names, cf.ncoeff);
    }
    free(cf);
}

/**
 * reset_coeff_intervals:
 * @cf: pointer to confidence intervals struct.
 * @alpha: nominal non-coverage, as decimal.
 *
 * Recomputes the intervals in @cf using the given value
 * of @alpha.
 *
 * Returns: 0 on success, non-zero if @alpha is out of
 * bounds.
 */

int reset_coeff_intervals(CoeffIntervals * cf, double alpha)
{
    double se, tbak = cf.t;
    int i;

    if (alpha <= 0 || alpha >= 1) {
        return GretlError.data;
    }

    if (cf.asy) {
        cf.t = normal_cdf_inverse(1 - alpha / 2);
    } else {
        cf.t = student_cdf_inverse(cf.df, 1 - alpha / 2);
    }

    for (i=0; i<cf.ncoeff; i++) {
        if (cf.maxerr[i] > 0) {
            se = cf.maxerr[i] / tbak;
            cf.maxerr[i] = se * cf.t;
        }
    }

    cf.alpha = alpha;

    return 0;
}

/**
 * gretl_model_get_coeff_intervals:
 * @pmod: pointer to gretl model.
 * @dset: dataset information.
 * @opt: TBA.
 *
 * Save the 95 percent confidence intervals for the parameter
 * estimates in @pmod.
 *
 * Returns: pointer to #CONFINT struct containing the results.
 */

CoeffIntervals * gretl_model_get_coeff_intervals(const(Model) * pmod,
                                                 const(DatasetRaw) * dset,
                                                 Gretlopt opt)
{
    CoeffIntervals * cf;
    char[32] pname;
    int nc = pmod.ncoeff;
    const(double) * b = pmod.coeff;
    const(double) * se = pmod.sderr;
    int offset = 0;
    int i, err = 0;

    cf = cast(CoeffIntervals *) malloc(CoeffIntervals.sizeof);
    if (cf == null) {
        return null;
    }

    if ((opt & Gretlopt.OPT_O) && pmod.ifc) {
        /* doing odds ratios: skip intercept */
        offset = 1;
        b++;
        se++;
        nc--;
    }

    cf.ncoeff = nc;
    cf.df = pmod.dfd;
    cf.opt = opt;
    cf.coeff = cf.maxerr = null;
    cf.names = null;

    cf.coeff = cast(double *) malloc(cf.ncoeff * double.sizeof);
    cf.maxerr = cast(double *) malloc(cf.ncoeff * double.sizeof);
    cf.names = strings_array_new(cf.ncoeff);

    if (cf.coeff == null || cf.maxerr == null || cf.names == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    cf.alpha = .05;

    if (ASYMPTOTIC_MODEL(cast(GretlCmdIndex) pmod.ci)) {
        cf.asy = 1;
        cf.t = normal_cdf_inverse(0.975);
    } else {
        cf.asy = 0;
        cf.t = tcrit95(pmod.dfd);
    }

    for (i=0; i<cf.ncoeff && !err; i++) {
        cf.coeff[i] = b[i];
        cf.maxerr[i] = (se[i] > 0)? cf.t * se[i] : 0.0;
        gretl_model_get_param_name(pmod, dset, i+offset, pname.ptr);
        cf.names[i] = gretl_strdup(pname.ptr);
        if (cf.names[i] == null) {
            err = GretlError.alloc;
        }
    }

 bailout:

    if (err) {
        free_coeff_intervals(cf);
        cf = null;
    }

    return cf;
}

GretlMatrixRaw * conf_intervals_matrix(CoeffIntervals * cf)
{
    GretlMatrixRaw * ret;
    const(double) * b = cf.coeff;
    const(double) * me = cf.maxerr;
    double se, lo, hi, pc;
    char ** S = null;
    char[32] head;
    int odds = (cf.opt & Gretlopt.OPT_O) != 0;
    int nrows = cf.ncoeff;
    int ncols = 3;
    int se_col = 0;
    int lo_col = 1;
    int hi_col = 2;
    int i, j;

    if (cf.opt & Gretlopt.OPT_E) {
        se_col = 1;
        lo_col++;
        hi_col++;
        ncols++;
    }

    ret = gretl_matrix_alloc(nrows, ncols);
    gretl_matrix_set_rownames(ret, strings_array_dup(cf.names, nrows));
    S = strings_array_new(ncols);
    pc = 100 * (1 - cf.alpha);

    for (i=0; i<nrows; i++) {
        if (na(me[i])) {
            lo = hi = NADBL;
        } else {
            lo = b[i] - me[i];
            hi = b[i] + me[i];
        }
        for (j=0; j<ncols; j++) {
            if (j == 0) {
                /* coeff or odds ratio */
                gretl_matrix_set(ret, i, j, odds ? exp(b[i]) : b[i]);
                if (i == 0) {
                    S[j] = gretl_strdup(odds? "odds ratio" : "coefficient");
                }
            } else if (j == se_col) {
                /* standard error */
                se = me[i] / cf.t;
                gretl_matrix_set(ret, i, j, odds ? exp(b[i]) * se : se);
                if (i == 0) {
                    S[j] = gretl_strdup("std. error");
                }
            } else if (j == lo_col) {
                /* lower limit */
                gretl_matrix_set(ret, i, j, odds ? exp(lo) : lo);
                if (i == 0) {
                    sprintf(head.ptr, "low%g", pc);
                    S[j] = gretl_strdup(head.ptr);
                }
            } else if (j == hi_col) {
                /* upper limit */
                gretl_matrix_set(ret, i, j, odds ? exp(hi) : hi);
                if (i == 0) {
                    sprintf(head.ptr, "high%g", pc);
                    S[j] = gretl_strdup(head.ptr);
                }
            }
        }
    }

    gretl_matrix_set_colnames(ret, S);

    return ret;
}

int gretl_is_simple_OLS(const(Model) * pmod)
{
    if (pmod.ci != GretlCmdIndex.ols || pmod.list == null) {
        return 0;
    } else if (pmod.list[0] == 3 && pmod.list[2] == 0) {
        return 1;
    } else if (pmod.list[0] == 2 && pmod.list[2] != 0) {
        return 1;
    } else {
        return 0;
    }
}

int gretl_is_arima_model(const(Model) * pmod)
{
    int d = gretl_model_get_int(pmod, "arima_d");
    int D = gretl_model_get_int(pmod, "arima_D");

    return (d > 0 || D > 0);
}

int gretl_is_between_model(const(Model) * pmod)
{
    if (pmod.ci == GretlCmdIndex.panel && (pmod.opt & Gretlopt.OPT_B)) {
        return 1;
    } else {
        return 0;
    }
}

/* "regular" here means pooled OLS, fixed effects or
   random effects
*/

int gretl_is_regular_panel_model(const(Model) * pmod)
{
    if (pmod.ci == GretlCmdIndex.ols && gretl_model_get_int(pmod, "pooled")) {
        return 1;
    } else if (pmod.ci != GretlCmdIndex.panel) {
        return 0;
    } else {
        Gretlopt opt = pmod.opt & (Gretlopt.OPT_P | Gretlopt.OPT_F | Gretlopt.OPT_U);

        return opt != 0;
    }
}

private bool arma_included(const(char) * m, int i)
{
    return (m == null || m[i] == '1');
}

/**
 * arma_model_nonseasonal_AR_order:
 * @pmod: pointer to gretl model.
 *
 * Returns: the non-seasonal autoregressive order of @pmod, or 0 if
 * @pmod is not an ARMA model.
 */

int arma_model_nonseasonal_AR_order(const(Model) * pmod)
{
    int p = 0;

    if (pmod.ci == GretlCmdIndex.arma) {
        p = pmod.list[1];
    }

    return p;
}

/**
 * arma_model_nonseasonal_MA_order:
 * @pmod: pointer to gretl model.
 *
 * Returns: the non-seasonal moving-average order of @pmod, or 0 if
 * @pmod is not an ARMA model.
 */

int arma_model_nonseasonal_MA_order(const(Model) * pmod)
{
    int q = 0;

    if (pmod.ci == GretlCmdIndex.arma) {
        if (gretl_is_arima_model(pmod)) {
            q = pmod.list[3];
        } else {
            q = pmod.list[2];
        }
    }

    return q;
}

/**
 * arma_model_max_AR_lag:
 * @pmod: pointer to gretl model.
 *
 * Returns: the maximum autoregressive lag in @pmod, or 0 if
 * @pmod is not an ARMA model. The maximum AR lag takes into
 * account any differencing (seasonal and/or non-seasonal) in
 * an ARIMA model.
 */

int arma_model_max_AR_lag(const(Model) * pmod)
{
    int pmax = 0;

    if (pmod.ci == GretlCmdIndex.arma) {
        int p = arma_model_nonseasonal_AR_order(pmod);
        int P = gretl_model_get_int(pmod, "arma_P");
        int s = gretl_model_get_int(pmod, "arma_pd");
        int d = gretl_model_get_int(pmod, "arima_d");
        int D = gretl_model_get_int(pmod, "arima_D");

        pmax = p + s * P;
        pmax += d + s * D;
    }

    return pmax;
}

/**
 * arma_model_max_MA_lag:
 * @pmod: pointer to gretl model.
 *
 * Returns: the maximum moving-average lag in @pmod, or 0 if
 * @pmod is not an ARMA model.
 */

int arma_model_max_MA_lag(const(Model) * pmod)
{
    int qmax = 0;

    if (pmod.ci == GretlCmdIndex.arma) {
        int q, Q, pd;

        q = arma_model_nonseasonal_MA_order(pmod);
        Q = gretl_model_get_int(pmod, "arma_Q");

        if (Q == 0) {
            qmax = q;
        } else {
            pd = gretl_model_get_int(pmod, "arma_pd");
            qmax = Q * pd + q;
        }
    }

    return qmax;
}

/* From Box and Jenkins, 1976, pp 506-7, "Program 4": a clever
   algorithm for "unscrambling the coefficients", or in effect
   producing reduced-form AR coefficients that take into account any
   differencing, seasonal and/or non-seasonal.
*/

private int ar_coeff_integrate(double * c0, int d, int D, int s, int pmax)
{
    int pstar = pmax + d + s * D;
    double * c1;
    int i, j, pp;

    c1 = cast(double *) malloc((pstar + 1) * double.sizeof);
    if (c1 == null) {
        return GretlError.alloc;
    }

    for (j=0; j<=pstar; j++) {
        c1[j] = 0.0;
    }

    pp = pmax;

    for (i=0; D>0 && i<D; i++) {
        for (j=0; j<=pstar; j++) {
            if (j < s) {
                c1[j] = c0[j];
            } else if (j <= pp) {
                c1[j] = c0[j] - c0[j-s];
            } else if (j <= pp + s) {
                c1[j] = -c0[j-s];
            }
        }
        pp += s;
        for (j=0; j<=pstar; j++) {
            c0[j] = c1[j];
        }
    }

    for (i=0; d>0 && i<d; i++) {
        for (j=0; j<=pstar; j++) {
            if (j < 1) {
                c1[j] = c0[j];
            } else if (j <= pp) {
                c1[j] = c0[j] - c0[j-1];
            } else if (j <= pp + 1) {
                c1[j] = -c0[j-1];
            }
        }
        pp += 1;
        for (j=0; j<=pstar; j++) {
            c0[j] = c1[j];
        }
    }

    free(c1);

    return 0;
}

private int arma_included_lags(int k, const(char) * mask)
{
    int i, nk = k;

    if (mask != null) {
        nk = 0;
        for (i=0; i<k; i++) {
            if (mask[i] == '1') nk++;
        }
    }

    return nk;
}

private int arma_AR_lags(const(Model) * pmod)
{
    const(char) * pmask = cast(const(char) *) gretl_model_get_data(pmod, "pmask");
    int p = arma_model_nonseasonal_AR_order(pmod);

    return arma_included_lags(p, pmask);
}

private int arma_MA_lags(const(Model) * pmod)
{
    const(char) * qmask = cast(const(char) *) gretl_model_get_data(pmod, "qmask");
    int q = arma_model_nonseasonal_MA_order(pmod);

    return arma_included_lags(q, qmask);
}

/**
 * arma_model_AR_MA_coeffs:
 * @pmod: pointer to gretl model.
 * @phi_star: pointer to receive AR coeff vector.
 * @theta_star: pointer to receive MA coeff vector.
 * @opt: use OPT_I to get the "integrated" variant, which
 * handles differencing.
 *
 * Creates consolidated versions of the AR and MA coefficient vectors
 * from @pmod.  If @pmod includes seasonal ARMA terms, the vectors are
 * suitably expanded to include the interactions between seasonal and
 * non-seasonal terms.  If OPT_I is given and the dependent variable
 * has been differenced, the AR coefficients are integrated to account
 * for the differencing.  These are the \Phi^* and \Theta^* as used by
 * Box and Jenkins for forecasting.
 *
 * Returns: 0 on success, non-zero on error.
 */

int arma_model_AR_MA_coeffs(const(Model) * pmod,
                            GretlVector ** phi_star,
                            GretlVector ** theta_star,
                            Gretlopt opt)
{
    GretlVector * ac = null;
    GretlVector * mc = null;
    int * ainfo;
    int err = 0;

    if (pmod.ci != GretlCmdIndex.arma) {
        err = GretlError.data;
    } else if ((ainfo = cast(int *) gretl_model_get_data(pmod, "ainfo")) == null) {
        fprintf(stderr, "AR_MA_coeffs: no 'ainfo' available!\n");
        err = GretlError.data;
    } else {
        const(double) * phi = null;
        const(double) * Phi = null;
        const(double) * theta = null;
        const(double) * Theta = null;
        const(char) * pmask = cast(const(char) *) gretl_model_get_data(pmod, "pmask");
        const(char) * qmask = cast(const(char) *) gretl_model_get_data(pmod, "qmask");
        int p  = ainfo[1];
        int q  = ainfo[2];
        int P  = ainfo[3];
        int Q  = ainfo[4];
        int np = ainfo[5];
        int nq = ainfo[6];
        int d  = ainfo[7];
        int D  = ainfo[8];
        int s  = ainfo[9];
        int pmax = p + s * P;
        int pstar = pmax + d + s * D;
        int qmax = q + s * Q;
        double x, y;
        int i, j, k, ii;

        if (pstar > 0) {
            /* we have some AR terms */
            ac = gretl_zero_matrix_new(pstar + 1, 1);
            if (ac == null) {
                err = GretlError.alloc;
            }
        }

        if (!err && qmax > 0) {
            /* we have some MA terms */
            mc = gretl_zero_matrix_new(qmax + 1, 1);
            if (mc == null) {
                gretl_vector_free(ac);
                ac = null;
                err = GretlError.alloc;
            }
        }

        if (!err) {
            phi = pmod.coeff + pmod.ifc; /* non-seasonal AR coeffs */
            Phi = phi + np;              /* seasonal AR coeffs */
            theta = Phi + P;             /* non-seasonal MA coeffs */
            Theta = theta + nq;          /* seasonal MA coeffs */
        }

        if (ac != null) {
            for (j=0; j<=P; j++) {
                x = (j == 0)? -1 : Phi[j-1];
                k = 0;
                for (i=0; i<=p; i++) {
                    y = (i == 0)? -1 :
                        (arma_included(pmask, i-1))? phi[k++] : 0;
                    ii = i + s * j;
                    ac.val[ii] -= x * y;
                }
            }
            if ((opt & Gretlopt.OPT_I) && (D > 0 || d > 0)) {
                ar_coeff_integrate(ac.val, d, D, s, pmax);
            }
        }

        if (mc != null) {
            for (j=0; j<=Q; j++) {
                x = (j == 0)? 1 : Theta[j-1];
                k = 0;
                for (i=0; i<=q; i++) {
                    y = (i == 0)? 1 :
                        (arma_included(qmask, i-1))? theta[k++] : 0;
                    ii = i + s * j;
                    mc.val[ii] += x * y;
                }
            }
        }
    }

    if (!err) {
        *phi_star = ac;
        *theta_star = mc;
    }

    return err;
}

/*
  arma_model_spectrum computes

   s2    C(exp(i*omega)) * C(exp(-i*omega))
  ---- * ----------------------------------
  2*pi   A(exp(i*omega)) * A(exp(-i*omega))

  for omega that goes from 0 to pi in T steps. The
  @phi matrix should contain the AR parameters and
  the @theta matrix the MA parameters.
*/

private GretlMatrixRaw * arma_model_spectrum(GretlMatrixRaw * phi,
                                          GretlMatrixRaw * theta,
                                          double s2, int T,
                                          int * err)
{
    GretlVector * ret = null;
    double num, den, nre_t, nim_t, dre_t, dim_t, xt;
    double c, s;
    double ar, ma;
    double scale;
    int nar, nma;
    int i, n, t;

    ret = gretl_matrix_alloc(T, 2);

    if (ret == null) {
        *err = GretlError.alloc;
        return ret;
    }

    nar = phi != null ? phi.rows : 0;
    nma = theta != null ? theta.rows : 0;

    n = max(nar, nma);
    scale = s2/M_2PI;

    for (t=0; t<T; t++) {
        xt = t * M_PI / (T-1);
        ar = nar > 0 ? gretl_vector_get(phi, 0) : 1.0;
        ma = nma > 0 ? gretl_vector_get(theta, 0) : 1.0;
        nre_t = ma * ma;
        dre_t = ar * ar;
        nim_t = dim_t = 0;

        for (i=1; i<n; i++) {
            ar = i < nar ? gretl_vector_get(phi, i) : 0.0;
            ma = i < nma ? gretl_vector_get(theta, i) : 0.0;
            c = cos(i * xt);
            s = sin(i * xt);
            nre_t += c * ma;
            nim_t += s * ma;
            dre_t += c * ar;
            dim_t += s * ar;
        }

        num = nre_t * nre_t + nim_t * nim_t;
        den = dre_t * dre_t + dim_t * dim_t;
        gretl_matrix_set(ret, t, 0, xt);
        gretl_matrix_set(ret, t, 1, scale * num/den);
    }

    return ret;
}

private double arima_diff(const(double) * x, int t,
                          int * delta, int k, int * err)
{
    double dxt = x[t];
    int i, p;

    if (na(dxt)) {
        *err = GretlError.missdata;
    }

    for (i=0; i<k && !*err; i++) {
        if (delta[i] != 0) {
            p = t - i - 1;
            if (p < 0 || na(x[p])) {
                dxt = NADBL;
                *err = GretlError.missdata;
            } else {
                dxt -= delta[i] * x[p];
            }
        }
    }

    return dxt;
}

/* get_arma_yvec: here we reconstruct the dependent variable from an
   ARMA model, differencing as we go in the ARIMA case.  If the model
   includes any regressors besides the constant we have to net out
   their effect in order to produce a periodogram that is comparable
   with the ARMA spectrum. A complication is that in the ARIMAX case
   the regressors may also have to be differenced: this is indicated
   by the presence of the "xdiff" flag on the model.
*/

private GretlVector * get_arma_yvec(const(Model) * pmod,
                                    const(DatasetRaw) * dset,
                                    int * err)
{
    GretlVector * y = null;
    const(double) * beta = null;
    double yt, xti;
    int * xlist = null;
    int * delta = null;
    int yno = gretl_model_get_depvar(pmod);
    int T = pmod.nobs;
    int i, j, t, t1 = pmod.t1;
    int s, d, D, nx = 0;
    int k = 0, xdiff = 0;

    if (yno < 1 || yno >= dset.v) {
        *err = GretlError.data;
        return null;
    }

    xlist = gretl_model_get_x_list(pmod);

    if (xlist != null) {
        for (i=1; i<=xlist[0]; i++) {
            if (xlist[i] >= dset.v) {
                *err = GretlError.data;
                break;
            } else if (xlist[i] != 0) {
                nx++;
            }
        }
        if (!*err && nx > 0) {
            beta = arma_model_get_x_coeffs(pmod);
        }
    }

    if (*err) {
        free(xlist);
        return null;
    }

    d = gretl_model_get_int(pmod, "arima_d");
    D = gretl_model_get_int(pmod, "arima_D");

    if (d > 0 || D > 0) {
        delta = arima_delta_coeffs(d, D, dset.pd);
        if (delta == null) {
            *err = GretlError.alloc;
        } else {
            k = d + dset.pd * D;
            xdiff = gretl_model_get_int(pmod, "xdiff");
        }
    }

    if (!*err) {
        y = gretl_matrix_alloc(T, 1);
        if (y == null) {
            *err = GretlError.alloc;
        }
    }

    if (!*err) {
        for (t=0; t<T && !*err; t++) {
            s = t + t1;
            if (delta != null) {
                yt = arima_diff(dset.Z[yno], s, delta, k, err);
            } else {
                yt = dset.Z[yno][s];
            }
            if (beta != null) {
                /* ARMAX: net out X_t\beta */
                for (i=1, j=0; i<=xlist[0]; i++) {
                    if (xlist[i] != 0) {
                        if (xdiff) {
                            xti = arima_diff(dset.Z[xlist[i]], s,
                                             delta, k, err);
                        } else {
                            xti = dset.Z[xlist[i]][s];
                        }
                        yt -= beta[j++] * xti;
                    }
                }
            }
            gretl_vector_set(y, t, yt);
        }
    }

    free(xlist);
    free(delta);

    if (*err) {
        gretl_matrix_free(y);
        y = null;
    }

    return y;
}

GretlMatrixRaw * arma_spectrum_plot_data(const(Model) * pmod,
                                      const(DatasetRaw) * dset,
                                      int * err)
{
    GretlMatrixRaw * pdata = null;
    GretlMatrixRaw * phi = null;
    GretlMatrixRaw * theta = null;
    GretlMatrixRaw * spec = null;
    GretlMatrixRaw * pergm = null;
    GretlMatrixRaw * y = null;
    int i, T = pmod.nobs;
    int grid = (T-1)/2;
    double s2 = pmod.sigma * pmod.sigma;

    *err = arma_model_AR_MA_coeffs(pmod, &phi, &theta, Gretlopt.OPT_NONE);
    if (*err) {
        return null;
    }

    /* different sign convention from forecasting */
    if (phi != null) {
        gretl_matrix_multiply_by_scalar(phi, -1.0);
    }

    spec = arma_model_spectrum(phi, theta, s2, grid, err);
    if (*err) {
        goto bailout;
    }

    y = get_arma_yvec(pmod, dset, err);

    if (!*err) {
        pergm = gretl_matrix_fft(y, err);
    }

    if (!*err) {
        pdata = gretl_matrix_alloc(grid, 4);
        if (pdata == null) {
            *err = GretlError.alloc;
        } else {
            Cmplx z;

            for (i=0; i<grid; i++) {
                gretl_matrix_set(pdata, i, 0, gretl_matrix_get(spec, i, 0));
                gretl_matrix_set(pdata, i, 1, gretl_matrix_get(spec, i, 1));
                z = gretl_cmatrix_get(pergm, i+1, 0);
                gretl_matrix_set(pdata, i, 2, z.r);
                gretl_matrix_set(pdata, i, 3, z.i);
            }
        }
    }

 bailout:

    gretl_matrix_free(phi);
    gretl_matrix_free(theta);
    gretl_matrix_free(spec);
    gretl_matrix_free(pergm);
    gretl_matrix_free(y);

    return pdata;
}

/**
 * regarma_model_AR_coeffs:
 * @pmod: pointer to gretl model.
 * @phi0: pointer to receive AR coeff vector.
 * @pp: pointer to receive length of @phi0.
 *
 * Creates a consolidated version of the AR coefficients from @pmod.
 * If @pmod includes seasonal AR terms the vector is suitably expanded
 * to include the interactions between seasonal and non-seasonal
 * terms, but it is not integrated with respect to any differencing of
 * the dependent variable.
 *
 * Returns: 0 on success, non-zero on error.
 */

int regarma_model_AR_coeffs(const(Model) * pmod,
                            double ** phi0,
                            int * pp)
{
    const(char) * pmask = cast(const(char) *) gretl_model_get_data(pmod, "pmask");
    int p = arma_model_nonseasonal_AR_order(pmod);
    int np = arma_included_lags(p, pmask);
    int P = gretl_model_get_int(pmod, "arma_P");
    int s = gretl_model_get_int(pmod, "arma_pd");
    const(double) * phi = null;
    const(double) * Phi = null;
    double * ac = null;
    double x, y;
    int i, j, k, ii, pmax;
    int err = 0;

    pmax = p + s * P;

    if (pmax == 0) {
        *pp = 0;
        return 0;
    }

    ac = cast(double *) malloc((pmax + 1) * double.sizeof);
    if (ac == null) {
        return GretlError.alloc;
    }

    phi = pmod.coeff + pmod.ifc;
    Phi = phi + np;

    for (i=0; i<=pmax; i++) {
        ac[i] = 0.0;
    }

    for (j=0; j<=P; j++) {
        x = (j == 0)? -1 : Phi[j-1];
        k = 0;
        for (i=0; i<=p; i++) {
            y = (i == 0)? -1 :
                (arma_included(pmask, i-1))? phi[k++] : 0;
            ii = i + s * j;
            ac[ii] -= x * y;
        }
    }

    *phi0 = ac;
    *pp = pmax;

    return err;
}

/**
 * arima_delta_coeffs:
 * @d: order of non-seasonal differencing (<= 2)
 * @D: order of seasonal differencing (<= 2)
 * @s: seasonal periodicity
 *
 * Returns: array of d + s * D coefficients of the lag operator in
 * the expansion of (1-L)^d * (1-L^s)^D. These are given in the
 * negative; for example, if d = 1 then c[0] = 1.
 */

int * arima_delta_coeffs(int d, int D, int s)
{
    int i, k = d + s * D;
    int * c = cast(int *) malloc(k * int.sizeof);

    if (c == null) {
        return null;
    }

    for (i=0; i<k; i++) {
        c[i] = 0;
    }

    if (d == 1) {
        c[0] = 1;
    } else if (d == 2) {
        c[0] = 2;
        c[1] = -1;
    }

    if (D > 0) {
        c[s-1] += 1;
        if (d > 0) {
            c[s] -= 1;
        }
        if (d == 2) {
            c[s] -= 1;
            c[s+1] += 1;
        }
        if (D == 2) {
            c[s-1] += 1;
            c[2*s-1] -= 1;
            if (d > 0) {
                c[s] -= 1;
                c[2*s] += 1;
            }
            if (d == 2) {
                c[s] -= 1;
                c[2*s] += 1;
                c[s+1] += 1;
                c[2*s+1] -= 1;
            }
        }
    }

    return c;
}

/**
 * arma_model_get_x_coeffs:
 * @pmod: pointer to gretl model.
 *
 * Returns: pointer to the array of coefficients on the exogenous
 * regressors in @pmod, or %NULL if the model is not ARMA or if there
 * are no such regressors.
 */

const(double) * arma_model_get_x_coeffs(const(Model) * pmod)
{
    const(double) * xc = null;

    if (pmod.ci == GretlCmdIndex.arma && gretl_model_get_int(pmod, "armax")) {
        xc = pmod.coeff;
        xc += pmod.ifc;
        xc += arma_AR_lags(pmod);
        xc += arma_MA_lags(pmod);
        xc += gretl_model_get_int(pmod, "arma_P");
        xc += gretl_model_get_int(pmod, "arma_Q");
    }

    return xc;
}

/**
 * arma_model_get_n_arma_coeffs:
 * @pmod: pointer to gretl model.
 *
 * Returns: the sum of the numbers of AR and MA coefficients in
 * @pmod.
 */

int arma_model_get_n_arma_coeffs(const(Model) * pmod)
{
    int npq = 0;

    if (pmod.ci == GretlCmdIndex.arma) {
        npq += arma_AR_lags(pmod);
        npq += arma_MA_lags(pmod);
        npq += gretl_model_get_int(pmod, "arma_P");
        npq += gretl_model_get_int(pmod, "arma_Q");
    }

    return npq;
}

/**
 * gretl_model_ahat_vec:
 * @pmod: pointer to gretl model.
 * @err: location to receive error code.
 *
 * If @pmod was estimated on panel data via fixed or random
 * effects, provides a column vector holding the "ahat" values
 * or estimated individual effects for the individuals
 * included in the active dataset when the model was estimated.
 * Some of these values may be NA if some units were not
 * actually included in estimation.
 *
 * Returns: allocated column vector on success or NULL on failure.
 */

GretlMatrixRaw * gretl_model_ahat_vec(const(Model) * pmod,
                                   int * err)
{
    GretlVector * a_vec = null;
    double * ahat;
    int i, s, t, N, T;

    ahat = cast(double *) gretl_model_get_data(pmod, "ahat");
    T = gretl_model_get_int(pmod, "panel_T");

    if (ahat == null || T == 0) {
        *err = GretlError.badstat;
        return null;
    }

    /* number of units included in model dataset */
    N = pmod.full_n / T;

    a_vec = gretl_column_vector_alloc(N);
    if (a_vec == null) {
        *err = GretlError.alloc;
        return null;
    }

    for (i=0; i<N; i++) {
        a_vec.val[i] = NADBL;
        for (t=0; t<T; t++) {
            s = i * T + t;
            if (!na(ahat[s])) {
                /* pick up first non-missing value */
                a_vec.val[i] = ahat[s];
                break;
            }
        }
    }

    return a_vec;
}

private int arbond_get_depvar(const(Model) * pmod)
{
    int i;

    for (i=1; i<=pmod.list[0]; i++) {
        if (pmod.list[i] == LISTSEP && i < pmod.list[0]) {
            return pmod.list[i+1];
        }
    }

    return 0;
}

private int arma_depvar_pos(const(Model) * pmod)
{
    int sep1 = 0, sep2 = 0;
    int seasonal = 0;
    int arima = 0;
    int i, dvpos;

    for (i=1; i<pmod.list[0]; i++) {
        if (pmod.list[i] == LISTSEP) {
            if (sep1 == 0) {
                sep1 = i;
            } else {
                sep2 = i;
            }
        }
    }

    if (sep2) {
        dvpos = sep2 + 1;
    } else if (sep1) {
        dvpos = sep1 + 1;
    } else {
        if (gretl_model_get_int(pmod, "arma_P") ||
            gretl_model_get_int(pmod, "arima_D") ||
            gretl_model_get_int(pmod, "arma_Q")) {
            seasonal = 1;
        }
        if (gretl_model_get_int(pmod, "arima_d") ||
            gretl_model_get_int(pmod, "arima_D")) {
            arima = 1;
        }
        if (arima) {
            dvpos = (seasonal)? 9 : 5;
        } else {
            dvpos = (seasonal)? 7 : 4;
        }
    }

    /* safety: should be impossible */
    if (dvpos == LISTSEP) {
        fprintf(stderr, "internal error in arma_depvar_pos\n");
        dvpos = 0;
    }

    return dvpos;
}

/**
 * gretl_model_get_depvar:
 * @pmod: pointer to gretl model.
 *
 * Returns: the ID number of the dependent variable in @pmod.
 */

int gretl_model_get_depvar(const(Model) * pmod)
{
    int dv = gretl_model_get_int(pmod, "yno");

    if (dv > 0) {
        return dv;
    }

    if (pmod != null && pmod.list != null) {
        if (pmod.ci == GretlCmdIndex.garch) {
            dv = pmod.list[4];
        } else if (pmod.ci == GretlCmdIndex.arma) {
            dv = pmod.list[arma_depvar_pos(pmod)];
        } else if (pmod.ci == GretlCmdIndex.dpanel) {
            dv = arbond_get_depvar(pmod);
        } else {
            dv = pmod.list[1];
        }
    }

    return dv;
}

/**
 * gretl_model_get_depvar_name:
 * @pmod: pointer to gretl model.
 * @dset: dataset information.
 *
 * Returns: the name of the dependent variable in @pmod.
 */

const(char) * gretl_model_get_depvar_name(const(Model) * pmod,
                                          const(DatasetRaw) * dset)
{
    int dv;

    if (pmod.depvar != null) {
        return pmod.depvar;
    }

    dv = gretl_model_get_int(pmod, "yno");

    if (dv == 0) {
        if (pmod != null && pmod.list != null) {
            if (pmod.ci == GretlCmdIndex.garch) {
                dv = pmod.list[4];
            } else if (pmod.ci == GretlCmdIndex.arma) {
                dv = pmod.list[arma_depvar_pos(pmod)];
            } else if (pmod.ci == GretlCmdIndex.dpanel) {
                dv = arbond_get_depvar(pmod);
            } else {
                dv = pmod.list[1];
            }
        }
    }

    return dset.varname[dv];
}

/**
 * gretl_model_get_x_list:
 * @pmod: pointer to gretl model.
 *
 * Returns: an allocated copy of the list of independent
 * variables included in @pmod, or %NULL on failure.
 */

int * gretl_model_get_x_list(const(Model) * pmod)
{
    if (pmod.xlist != null) {
        return gretl_list_copy(pmod.xlist);
    } else {
        return null;
    }
}

/**
 * gretl_model_get_secondary_list:
 * @pmod: model to examine.
 *
 * Retrieve an allocated copy of the secondary list from
 * @pmod: e.g. the list of instruments in the case of %IVREG, or
 * the selection equation list for %HECKIT.
 *
 * Returns: allocated list or %NULL on error.
 */

int * gretl_model_get_secondary_list(const(Model) * pmod)
{
    int * list = null;
    int pos;

    pos = gretl_list_separator_position(pmod.list);

    if (pos > 0) {
        list = gretl_list_copy_from_pos(pmod.list, pos + 1);
    }

    return list;
}

/**
 * gretl_model_get_y_list:
 * @pmod: model to examine.
 *
 * Retrieve an allocated copy of the list of dependent variables
 * for @pmod: in almost all cases this will have a single
 * element; an exception is biprobit.
 *
 * Returns: allocated list or %NULL on error.
 */

int * gretl_model_get_y_list(const(Model) * pmod)
{
    int * list = null;

    if (pmod.list == null) {
        return null;
    }

    if (pmod.ci == GretlCmdIndex.biprobit) {
        list = gretl_list_new(2);
        if (list != null) {
            list[1] = pmod.list[1];
            list[2] = pmod.list[2];
        }
    } else {
        list = gretl_list_new(1);
        if (list != null) {
            list[1] = gretl_model_get_depvar(pmod);
        }
    }

    return list;
}

/**
 * gretl_model_allocate_storage:
 * @pmod: pointer to model.
 *
 * Allocates space for coefficients and standard errors,
 * residuals and fitted values in @pmod. The sizes of
 * the arrays are based on the %ncoeff and %full_n
 * members of @pmod, which must be set first. The
 * residuals and fitted values are initialized to
 * gretl's missing value.
 *
 * Returns: 0 on success, %E_ALLOC on error.
 */

int gretl_model_allocate_storage(Model * pmod)
{
    int k = pmod.ncoeff;
    int T = pmod.full_n;

    if (k > 0) {
        int i;

        pmod.coeff = cast(double *) malloc(k * double.sizeof);
        if (pmod.coeff == null) {
            return GretlError.alloc;
        }
        pmod.sderr = cast(double *) malloc(k * double.sizeof);
        if (pmod.sderr == null) {
            return GretlError.alloc;
        }
        for (i=0; i<k; i++) {
            pmod.coeff[i] = pmod.sderr[i] = NADBL;
        }
    }

    if (T > 0) {
        int t;

        pmod.uhat = cast(double *) malloc(T * double.sizeof);
        if (pmod.uhat == null) {
            return GretlError.alloc;
        }
        pmod.yhat = cast(double *) malloc(T * double.sizeof);
        if (pmod.yhat == null) {
            return GretlError.alloc;
        }
        for (t=0; t<T; t++) {
            pmod.uhat[t] = pmod.yhat[t] = NADBL;
        }
    }

    return 0;
}

/**
 * gretl_model_new_vcv:
 * @pmod: pointer to model.
 * @nelem: pointer to receive number of elements in
 * the packed array, or %NULL;
 *
 * Allocates space for a packed coefficient covariance matrix
 * in @pmod (if such space is not already allocated).  Sets
 * all entries in the array to zero.
 *
 * Returns: 0 on success, %E_ALLOC on error.
 */

int gretl_model_new_vcv(Model * pmod, int * nelem)
{
    int nv = pmod.ncoeff;
    int nxpx = (nv * nv + nv) / 2;
    int i, err = 0;

    if (pmod.vcv == null) {
        /* not already allocated */
        pmod.vcv = cast(double *) malloc(nxpx * double.sizeof);
        if (pmod.vcv == null) {
            err = GretlError.alloc;
        }
    }

    if (pmod.vcv != null) {
        for (i=0; i<nxpx; i++) {
            pmod.vcv[i] = 0.0;
        }
        if (nelem != null) {
            *nelem = nxpx;
        }
    }

    return err;
}

private double vcv_get_se(double vii, int restricted)
{
    if (restricted) {
        vii = fabs(vii) < 1.0e-17 ? 0.0 : vii;
    }

    return (na(vii) || vii < 0)? NADBL : sqrt(vii);
}

/**
 * gretl_model_write_vcv:
 * @pmod: pointer to model.
 * @V: covariance matrix.
 *
 * Write the covariance matrix @V into the model @pmod, using the
 * special packed format that is required by the MODEL struct,
 * and set the standard errors to the square root of the diagonal
 * elements of this matrix.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_model_write_vcv(Model * pmod, const(GretlMatrixRaw) * V)
{
    int i, j, k, n;
    double x;
    double * tmp;
    int err = 0;

    if (gretl_is_null_matrix(V)) {
        return 0; /* no-op */
    }

    if (V.cols != V.rows) {
        return GretlError.nonconf;
    }

    k = V.rows;
    n = (k * k + k) / 2;

    /* reallocate model vcv in case it's wrongly sized */
    tmp = cast(double *) realloc(pmod.vcv, n * double.sizeof);
    if (tmp == null) {
        err = GretlError.alloc;
    } else {
        pmod.vcv = tmp;
    }

    /* same for standard errors array */
    if (!err) {
        tmp = cast(double *) realloc(pmod.sderr, k * double.sizeof);
        if (tmp == null) {
            err = GretlError.alloc;
        } else {
            pmod.sderr = tmp;
        }
    }

    if (!err) {
        int restricted, idx = 0;

        restricted = gretl_model_get_int(pmod, "restricted");

        for (i=0; i<k; i++) {
            for (j=i; j<k; j++) {
                x = gretl_matrix_get(V, i, j);
                pmod.vcv[idx++] = x;
                if (i == j) {
                    pmod.sderr[i] = vcv_get_se(x, restricted);
                }
            }
        }
    }

    return err;
}

private GretlMatrixRaw * make_cluster_vals(Model * pmod, int cvar,
                                        const(DatasetRaw) * dset,
                                        int * err)
{
    GretlMatrixRaw * cvals = null;
    double ct;
    double * cdata;
    int s, t;

    cdata = cast(double *) malloc(pmod.nobs * double.sizeof);

    if (cdata == null) {
        *err = GretlError.alloc;
    } else {
        s = 0;
        for (t=pmod.t1; t<=pmod.t2 && !*err; t++) {
            if (!model_missing(pmod, t)) {
                ct = dset.Z[cvar][t];
                if (na(ct)) {
                    *err = GretlError.missdata;
                } else {
                    cdata[s++] = ct;
                }
            }
        }
    }

    if (!*err) {
        cvals = gretl_matrix_values(cdata, pmod.nobs, Gretlopt.OPT_S, err);
        if (!*err && gretl_vector_get_length(cvals) < 2) {
            gretl_matrix_free(cvals);
            cvals = null;
            *err = GretlError.data;
        }
    }

    free(cdata);

    return cvals;
}

private int model_make_clustered_GG(Model * pmod, int ci,
                                    const(GretlMatrixRaw) * G,
                                    GretlMatrixRaw * GG,
                                    const(DatasetRaw) * dset,
                                    int * pcvar,
                                    int * pnc)
{
    GretlMatrixRaw * cvals = null;
    GretlMatrixRaw * Gi = null;
    const(double) * cZ = null;
    const(char) * cname;
    int cvar, n_c = 0;
    int k = G.cols;
    int i, j, t;
    int err = 0;

    cname = get_optval_string(cast(GretlCmdIndex) ci, Gretlopt.OPT_C);
    if (cname == null) {
        return GretlError.parse;
    }

    cvar = current_series_index(dset, cname);

    if (cvar < 1 || cvar >= dset.v) {
        err = GretlError.unkvar;
    } else {
        cvals = make_cluster_vals(pmod, cvar, dset, &err);
    }

    if (!err) {
        Gi = gretl_matrix_alloc(1, k);
        if (Gi == null) {
            gretl_matrix_free(cvals);
            err = GretlError.alloc;
        }
    }

    if (err) {
        return err;
    }

    n_c = gretl_vector_get_length(cvals);
    cZ = dset.Z[cvar];

    for (i=0; i<n_c; i++) {
        double cvi = cvals.val[i];
        int s = 0;

        gretl_matrix_zero(Gi);

        for (t=pmod.t1; t<=pmod.t2; t++) {
            if (!model_missing(pmod, t)) {
                if (cZ[t] == cvi) {
                    for (j=0; j<k; j++) {
                        Gi.val[j] += gretl_matrix_get(G, s, j);
                    }
                }
                s++;
            }
        }

        gretl_matrix_multiply_mod(Gi, GretlMatrixMod.transpose,
                                  Gi, GretlMatrixMod.none,
                                  GG, GretlMatrixMod.cumulate);
    }

    if (!err) {
        *pcvar = cvar;
        *pnc = n_c;
    }

    gretl_matrix_free(cvals);
    gretl_matrix_free(Gi);

    return err;
}

private alias BiprobitAdjFn = extern (C) void function(Model *, GretlMatrixRaw *);

private int do_biprobit_adjustment(Model * pmod,
                                   GretlMatrixRaw * V)
{
    BiprobitAdjFn biprobit_adj;

    biprobit_adj = cast(BiprobitAdjFn) get_plugin_function("biprobit_adjust_vcv");

    if (biprobit_adj == null) {
        return GretlError.fopen;
    } else {
        (*biprobit_adj)(pmod, V);
        return 0;
    }
}

/**
 * gretl_model_add_QML_vcv:
 * @pmod: pointer to model.
 * @ci: command index for model.
 * @H: inverse of the (negative) Hessian, k x k.
 * @G: score matrix, T x k.
 * @dset: pointer to dataset (can be NULL if not doing
 * clustering).
 * @opt: may include OPT_C for cluster-robust variant.
 * @pV: location to receive the full covariance matrix, or NULL.
 *
 * Write a QML covariance matrix into the model @pmod, and set
 * the standard errors to the square root of the diagonal
 * elements of this matrix. The @ci argument, specifying the
 * estimator for @pmod, is required only if @opt includes
 * OPT_C; otherwise it is ignored.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_model_add_QML_vcv(Model * pmod, int ci,
                            const(GretlMatrixRaw) * H,
                            const(GretlMatrixRaw) * G,
                            const(DatasetRaw) * dset,
                            Gretlopt opt,
                            GretlMatrixRaw ** pV)
{
    GretlMatrixRaw * GG = null;
    GretlMatrixRaw * V = null;
    int cvar = 0, n_c = 0;
    int k = H.rows;
    int err = 0;

    V = gretl_matrix_alloc(k, k);

    if (V == null) {
        return GretlError.alloc;
    }

    if (!err) {
        if (opt & Gretlopt.OPT_C) {
            /* clustered */
            GG = gretl_zero_matrix_new(k, k);
            if (GG == null) {
                err = GretlError.alloc;
            } else {
                err = model_make_clustered_GG(pmod, ci, G, GG,
                                              dset, &cvar, &n_c);
            }
        } else if (opt & Gretlopt.OPT_N) {
            /* Newey-West */
            GG = newey_west_OPG(G, &err);
        } else {
            /* regular QML using OPG */
            GG = gretl_matrix_XTX_new(G);
            if (GG == null) {
                err = GretlError.alloc;
            }
        }
    }

    if (!err) {
        err = gretl_matrix_qform(H, GretlMatrixMod.none, GG,
                                 V, GretlMatrixMod.none);
        if (!err && (opt & Gretlopt.OPT_C)) {
            /* clustering: use stata-style df adjustment */
            double dfc = n_c / (n_c - 1.0);

            gretl_matrix_multiply_by_scalar(V, dfc);
        }
    }

    if (!err) {
        if (ci == GretlCmdIndex.biprobit) {
            do_biprobit_adjustment(pmod, V);
        }
        err = gretl_model_write_vcv(pmod, V);
    }

    if (!err) {
        if (opt & Gretlopt.OPT_C) {
            const(char) * cname = dset.varname[cvar];

            gretl_model_set_int(pmod, "n_clusters", n_c);
            gretl_model_set_cluster_vcv_info(pmod, cname, null);
            pmod.opt |= Gretlopt.OPT_C;
        } else {
            MLVCVType vt = (opt & Gretlopt.OPT_N)? MLVCVType.hac : MLVCVType.qml;

            gretl_model_set_vcv_info(pmod, VCVMajorType.ml, vt);
        }
        pmod.opt |= Gretlopt.OPT_R;
    }

    gretl_matrix_free(GG);

    if (pV != null) {
        *pV = V;
    } else {
        gretl_matrix_free(V);
    }

    return err;
}

/**
 * gretl_model_add_hessian_vcv:
 * @pmod: pointer to model.
 * @H: inverse of the (negative) Hessian.
 *
 * Write @H into the model @pmod as its covariance matrix, and
 * set the standard errors to the square roots of the diagonal
 * elements of this matrix.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_model_add_hessian_vcv(Model * pmod,
                                const(GretlMatrixRaw) * H)
{
    int err = gretl_model_write_vcv(pmod, H);

    if (!err) {
        gretl_model_set_vcv_info(pmod, VCVMajorType.ml, MLVCVType.hessian);
    }

    return err;
}

/**
 * gretl_model_add_OPG_vcv:
 * @pmod: pointer to model.
 * @G: T x k gradient matrix.
 * @pV: location to receive the full covariance matrix, or NULL.
 *
 * Compute (G'G)^{-1}, write this into @pmod as its covariance matrix,
 * and set the standard errors to the square roots of the diagonal
 * elements of this matrix.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_model_add_OPG_vcv(Model * pmod,
                            const(GretlMatrixRaw) * G,
                            GretlMatrixRaw ** pV)
{
    GretlMatrixRaw * GG = gretl_matrix_XTX_new(G);
    int err = 0;

    if (GG == null) {
        return GretlError.alloc;
    }

    err = gretl_invert_symmetric_matrix(GG);

    if (!err) {
        err = gretl_model_write_vcv(pmod, GG);
        if (!err) {
            gretl_model_set_vcv_info(pmod, VCVMajorType.ml, MLVCVType.op);
        }
    } else {
        /* diagnostic: recreate and print G'G */
        gretl_matrix_multiply_mod(G, GretlMatrixMod.transpose,
                                  G, GretlMatrixMod.none,
                                  GG, GretlMatrixMod.none);
        gretl_matrix_print(GG, "non-p.d. G'G");
    }

    if (pV != null) {
        *pV = GG;
    } else {
        gretl_matrix_free(GG);
    }

    return err;
}

private double * copy_vcv_subset(const(Model) * pmod)
{
    double * V;
    int nc = pmod.ncoeff;
    int k = pmod.list[0] - 1;
    int n = k * (k + 1) / 2;
    int i, j;

    V = cast(double *) malloc(n * double.sizeof);
    if (V == null) {
        return null;
    }

    for (i=0; i<k; i++) {
        for (j=0; j<=i; j++) {
            V[ijton(i, j, k)] = pmod.vcv[ijton(i, j, nc)];
        }
    }

    return V;
}

/**
 * gretl_model_get_vcv:
 * @pmod: pointer to model.
 * @dset: dataset information.
 *
 * Supplies the caller with a copy of the variance-covariance
 * matrix for the parameter estimates in @pmod, in a format
 * suitable for printing.  See also free_vcv().  To get the
 * covariance matrix as a gretl_matrix, see
 * gretl_vcv_matrix_from_model().
 *
 * Returns: #VMatrix struct or %NULL on error.
 */

VMatrix * gretl_model_get_vcv(Model * pmod, const(DatasetRaw) * dset)
{
    char[VNAMELEN] varname;
    int i, k = pmod.ncoeff;
    int special = 0;
    VMatrix * vcv;

    vcv = vmatrix_new();

    if (vcv == null) {
        return null;
    }

    /* special for fixed effects panel model: strip out the
       vcv elements for per-unit dummies, if present
    */
    if (pmod.ci == GretlCmdIndex.panel) {
        int k2 = pmod.list[0] - 1;

        if (k > k2) {
            k = k2;
            special = 1;
        }
    }

    vcv.names = strings_array_new(k);
    if (vcv.names == null) {
        free(vcv);
        return null;
    }

    for (i=0; i<k; i++) {
        gretl_model_get_param_name(pmod, dset, i, varname.ptr);
        vcv.names[i] = gretl_strdup(varname.ptr);
        if (vcv.names[i] == null) {
            free_vmatrix(vcv);
            return null;
        }
    }

    if (pmod.vcv == null && makevcv(pmod, pmod.sigma)) {
        free_vmatrix(vcv);
        return null;
    }

    if (special) {
        /* copy subset of vcv */
        vcv.vec = copy_vcv_subset(pmod);
    } else {
        /* copy full vcv */
        vcv.vec = copyvec(pmod.vcv, k * (k + 1)  / 2);
    }

    if (vcv.vec == null) {
        free_vmatrix(vcv);
        return null;
    }

    vcv.ci = pmod.ci;
    vcv.dim = k;
    vcv.t1 = pmod.t1;
    vcv.t2 = pmod.t2;

    return vcv;
}

/**
 * gretl_model_get_vcv_element:
 * @pmod: pointer to model.
 * @i: row (0-based).
 * @j: column(0-based).
 * @np: the number of parameters represented in @pmod's
 * covariance matrix.
 *
 * Returns: the (@i, @j) element of @pmod's covariance matrix,
 * or #NADBL on failure.
 */

double gretl_model_get_vcv_element(const(Model) * pmod,
                                   int i, int j,
                                   int np)
{
    if (pmod.vcv == null) {
        return NADBL;
    } else {
        int k = ijton(i, j, np);

        return pmod.vcv[k];
    }
}

/**
 * gretl_model_write_coeffs:
 * @pmod: pointer to model.
 * @b: array of coefficients.
 * @k: number of elements in @b.
 *
 * Write the coefficients @b into the model @pmod, whose
 * coefficient array is resized appropriately if need be.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_model_write_coeffs(Model * pmod, double * b, int k)
{
    size_t sz = k * double.sizeof;
    int err = 0;

    if (pmod.coeff == null || pmod.ncoeff != k) {
        double * tmp = cast(double *) realloc(pmod.coeff, sz);

        if (tmp == null) {
            err = GretlError.alloc;
        } else {
            pmod.coeff = tmp;
        }
    }

    if (!err) {
        memcpy(pmod.coeff, b, sz);
        pmod.ncoeff = k;
    }

    return err;
}

/**
 * impose_model_smpl:
 * @pmod: pointer to model.
 * @dset: dataset information.
 *
 * Sets on @dset the sample range (starting and ending
 * observations) that was in effect when @pmod was estimated.
 * This is not always the same as the data range over which
 * @pmod was actually estimated (e.g. in case of
 * autoregressive models, where observations are dropped
 * to allow for lags).
 */

void impose_model_smpl(const(Model) * pmod, DatasetRaw * dset)
{
    dset.t1 = pmod.smpl.t1;
    dset.t2 = pmod.smpl.t2;
}

/* HAVE_X12A is taken to be defined: see log.md */

private void maybe_delete_x12_file(const(Model) * pmod)
{
    char * fname = null;

    fname = cast(char *) gretl_model_get_data(pmod, "x12a_output");
    if (fname != null) {
        gretl_remove(fname);
    }
}

private void destroy_all_data_items(Model * pmod)
{
    int i;

    if (pmod.n_data_items == 0) {
        return;
    }

    maybe_delete_x12_file(pmod);

    for (i=0; i<pmod.n_data_items; i++) {
        free_model_data_item(pmod.data_items[i]);
    }

    free(pmod.data_items);
    pmod.data_items = null;
}

private int discard_model_data_item(Model * pmod, const(char) * key,
                                    int free_data)
{
    ModelDataItem * junk = null;
    int i, targ = 0;
    int err = 0;

    for (i=0; i<pmod.n_data_items; i++) {
        if (!strcmp(key, pmod.data_items[i].key)) {
            junk = pmod.data_items[i];
            targ = i;
            break;
        }
    }

    if (junk == null) {
        err = 1;
    } else {
        int new_n_items = pmod.n_data_items - 1;

        if (new_n_items == 0) {
            free(pmod.data_items);
            pmod.data_items = null;
        } else {
            ModelDataItem ** items;

            for (i=targ; i<new_n_items; i++) {
                pmod.data_items[i] = pmod.data_items[i+1];
            }
            items = cast(ModelDataItem **) realloc(pmod.data_items,
                                                   new_n_items * (ModelDataItem *).sizeof);
            if (items != null) {
                pmod.data_items = items;
            }
        }

        pmod.n_data_items = new_n_items;

        if (free_data) {
            /* deep free the data item */
            free_model_data_item(junk);
        } else {
            /* just free the item, not the actual data */
            free(junk.key);
            free(junk);
        }
    }

    return err;
}

/**
 * gretl_model_destroy_data_item:
 * @pmod: pointer to model.
 * @key: key string.
 *
 * Looks up the data pointer, attached to @pmod, that is
 * identified by @key, and if a pointer is found, frees
 * it (or applies the destructor function that was set for
 * the item, if any) and removes it from the model's list of
 * data items.  If you want to remove the item from the
 * model's list without freeing the underlying data pointer,
 * use gretl_model_detach_data_item().
 *
 * Returns: 0 on success, 1 on failure (pointer not found).
 */

int gretl_model_destroy_data_item(Model * pmod, const(char) * key)
{
    return discard_model_data_item(pmod, key, 1);
}

/**
 * gretl_model_detach_data_item:
 * @pmod: pointer to model.
 * @key: key string.
 *
 * Look up the data item, attached to @pmod, that is
 * identified by @key, and if an item is found remove
 * it from the model's list of such items.  The data
 * pointer associated with @key is not touched.  If you
 * want the underlying resources associated with @key to be
 * freed, use gretl_model_destroy_data_item().
 *
 * Returns: 0 on success, 1 on failure (key not found).
 */

int gretl_model_detach_data_item(Model * pmod, const(char) * key)
{
    return discard_model_data_item(pmod, key, 0);
}

/**
 * gretl_model_set_auxiliary:
 * @pmod: pointer to model.
 * @aux: code indicating a model's function in an
 * auxiliary role (typically, in relation to a hypothesis
 * test on another model).
 *
 * Sets an auxiliary code on @pmod, which may be relevant
 * for how the model is printed.
 */

void gretl_model_set_auxiliary(Model * pmod, ModelAuxCode aux)
{
    pmod.aux = aux;
}

/**
 * gretl_model_add_arinfo:
 * @pmod: pointer to model.
 * @nterms: number of autoregressive coefficients.
 *
 * Performs initial setup for structure to hold info
 * on autoregressive coefficients.
 *
 * Returns: 0 on success, 1 on error.
 */

int gretl_model_add_arinfo(Model * pmod, int nterms)
{
    int i;

    pmod.arinfo = cast(ArInfo *) malloc(ArInfo.sizeof);
    if (pmod.arinfo == null) {
        return 1;
    }

    if (nterms == 0) {
        pmod.arinfo.arlist = null;
        pmod.arinfo.rho = null;
        pmod.arinfo.sderr = null;
        return 0;
    }

    pmod.arinfo.arlist = gretl_list_new(nterms);
    if (pmod.arinfo.arlist == null) {
        free(pmod.arinfo);
        pmod.arinfo = null;
        return 1;
    }

    pmod.arinfo.rho = cast(double *) malloc(nterms * double.sizeof);
    if (pmod.arinfo.rho == null) {
        free(pmod.arinfo.arlist);
        free(pmod.arinfo);
        pmod.arinfo = null;
        return 1;
    }

    pmod.arinfo.sderr = cast(double *) malloc(nterms * double.sizeof);
    if (pmod.arinfo.sderr == null) {
        free(pmod.arinfo.arlist);
        free(pmod.arinfo.rho);
        free(pmod.arinfo);
        pmod.arinfo = null;
        return 1;
    }

    for (i=0; i<nterms; i++) {
        pmod.arinfo.sderr[i] = pmod.arinfo.rho[i] = NADBL;
    }

    return 0;
}

private void model_stats_init(Model * pmod)
{
    int i;

    pmod.ess = pmod.tss = NADBL;
    pmod.sigma = NADBL;
    pmod.rsq = pmod.adjrsq = NADBL;
    pmod.fstt = pmod.chisq = NADBL;
    pmod.lnL = NADBL;
    pmod.ybar = pmod.sdy = NADBL;
    pmod.dw = pmod.rho = NADBL;

    for (i=0; i<InfoCriterion.max; i++) {
        pmod.criterion[i] = NADBL;
    }
}

private void gretl_model_init_pointers(Model * pmod)
{
    pmod.list = null;
    pmod.xlist = null;
    pmod.submask = null;
    pmod.missmask = null;
    pmod.coeff = null;
    pmod.sderr = null;
    pmod.yhat = null;
    pmod.uhat = null;
    pmod.xpx = null;
    pmod.vcv = null;
    pmod.arinfo = null;
    pmod.name = null;
    pmod.depvar = null;
    pmod.params = null;
    pmod.tests = null;
    pmod.dataset = null;
    pmod.data_items = null;
}

/**
 * gretl_model_init:
 * @pmod: pointer to model.
 * @dset: pointer to dataset.
 *
 * Initializes a gretl #MODEL, including setting its pointer members
 * to %NULL. This initialization should be done if the caller has
 * declared a #MODEL struct directly, rather than obtaining a pointer to
 * #MODEL using gretl_model_new() (in which case the initialization is
 * done automatically).
 */

void gretl_model_init(Model * pmod, const(DatasetRaw) * dset)
{
    if (pmod == null) return;

    pmod.ID = 0;
    pmod.refcount = 0;
    pmod.ci = 0;
    pmod.opt = Gretlopt.OPT_NONE;
    pmod.full_n = 0;
    pmod.t1 = pmod.t2 = 0;
    pmod.nobs = 0;

    if (dset != null) {
        pmod.smpl.t1 = dset.t1;
        pmod.smpl.t2 = dset.t2;
        pmod.smpl.rseed = dset.rseed;
    } else {
        pmod.smpl.t1 = 0;
        pmod.smpl.t2 = 0;
        pmod.smpl.rseed = 0;
    }

    pmod.ncoeff = 0;
    pmod.ntests = 0;
    pmod.nparams = 0;
    pmod.errcode = 0;
    pmod.ifc = 0;
    pmod.nwt = 0;
    pmod.aux = ModelAuxCode.none;
    pmod.esttime = 0;

    model_stats_init(pmod);

    gretl_model_init_pointers(pmod);
    pmod.n_data_items = 0;
}

/**
 * gretl_model_smpl_init:
 * @pmod: pointer to model.
 * @dset: dataset information.
 *
 * Records the start and end of the current sample range in
 * the model @pmod, which may be necessary for future reference
 * if a hypothesis test is to be performed.  Note that this
 * sample range may not be the same as the data range over
 * which the model is actually estimated (for example, in the
 * case of autoregressive models where observations have to
 * be dropped to allow for lags).
 */

void gretl_model_smpl_init(Model * pmod, const(DatasetRaw) * dset)
{
    pmod.smpl.t1 = dset.t1;
    pmod.smpl.t2 = dset.t2;
    pmod.smpl.rseed = dset.rseed;
}

/**
 * gretl_model_new:
 *
 * Allocates memory for a gretl #MODEL struct and initializes the struct,
 * using gretl_model_init().
 *
 * Returns: pointer to model (or %NULL if allocation fails).
 */

Model * gretl_model_new()
{
    Model * pmod = cast(Model *) malloc(Model.sizeof);

    if (pmod != null) {
        gretl_model_init(pmod, null);
    }

    return pmod;
}

/**
 * gretl_model_array_new:
 * @n: number of models in array.
 *
 * Allocates memory for an array of @n gretl #MODEL structs and
 * initializes each model using gretl_model_init().
 *
 * Returns: pointer to models array (or %NULL if allocation fails).
 */

Model ** gretl_model_array_new(int n)
{
    Model ** models;
    int i, j;

    models = cast(Model **) malloc(n * (Model *).sizeof);

    if (models != null) {
        for (i=0; i<n; i++) {
            models[i] = gretl_model_new();
            if (models[i] == null) {
                for (j=0; j<i; j++) {
                    free(models[i]);
                }
                free(models);
                models = null;
                break;
            }
        }
    }

    return models;
}

/**
 * allocate_working_model:
 *
 * Allocates memory for gretl #MODEL struct and initializes it.
 * The model is "protected" against deletion.
 *
 * Returns: pointer to model (or %NULL if allocation fails).
 */

Model * allocate_working_model()
{
    Model * model = gretl_model_new();

    if (model != null) {
        int err = gretl_model_protect(model);

        if (err) {
            gretl_model_free(model);
            model = null;
        }
    }

    return model;
}

/**
 * gretl_model_array_destroy:
 * @models: array of gretl models.
 * @n: number of models in array.
 *
 * Frees all resources associated with an array of models, which
 * should have been obtained via gretl_model_array_new().
 */

void gretl_model_array_destroy(Model ** models, int n)
{
    int i;

    if (models != null) {
        for (i=0; i<n; i++) {
            clear_model(models[i]);
            free(models[i]);
        }
        free(models);
    }
}

void destroy_working_model(Model * model)
{
    if (model != null) {
        gretl_model_free_on_exit(model);
    }
}

private void clear_ar_info(Model * pmod)
{
    if (pmod.arinfo.arlist) {
        free(pmod.arinfo.arlist);
    }
    if (pmod.arinfo.rho) {
        free(pmod.arinfo.rho);
    }
    if (pmod.arinfo.sderr) {
        free(pmod.arinfo.sderr);
    }

    free(pmod.arinfo);
    pmod.arinfo = null;
}

/**
 * gretl_model_destroy_tests:
 * @pmod: pointer to model.
 *
 * Clears any hypothesis test structs that have been attached
 * to @pmod.
 */

void gretl_model_destroy_tests(Model * pmod)
{
    if (pmod != null && pmod.ntests > 0) {
        int i;

        for (i=0; i<pmod.ntests; i++) {
            if (pmod.tests[i].param != null) {
                free(pmod.tests[i].param);
            }
        }
        free(pmod.tests);
        pmod.tests = null;
        pmod.ntests = 0;
    }
}

/**
 * clear_model:
 * @pmod: pointer to model.
 *
 * Clears a gretl #MODEL, freeing all allocated storage and setting
 * pointer members to %NULL.  Also frees any data pointers attached
 * via gretl_model_set_data().  The model pointer itself is not
 * freed, so this function may be called on a #MODEL which has been
 * declared directly by the caller; in that case the caller should
 * pass the address of the #MODEL in question.
 */

void clear_model(Model * pmod)
{
    if (pmod != null) {
        if (pmod.list != null) free(pmod.list);
        if (pmod.xlist != null) free(pmod.xlist);
        if (pmod.missmask != null) free(pmod.missmask);
        if (pmod.coeff != null) free(pmod.coeff);
        if (pmod.sderr != null) free(pmod.sderr);
        if (pmod.yhat != null) free(pmod.yhat);
        if (pmod.uhat != null) free(pmod.uhat);
        if (pmod.xpx != null) free(pmod.xpx);
        if (pmod.vcv != null) free(pmod.vcv);
        if (pmod.name != null) free(pmod.name);
        if (pmod.depvar != null) free(pmod.depvar);

        if (pmod.submask != null) {
            free_subsample_mask(pmod.submask);
        }

        if (pmod.arinfo != null) {
            clear_ar_info(pmod);
        }
        if (pmod.params != null) {
            strings_array_free(pmod.params, pmod.nparams);
        }

        destroy_dataset(pmod.dataset);
        gretl_model_destroy_tests(pmod);
        destroy_all_data_items(pmod);
    }

    /* this may be redundant */
    gretl_model_init(pmod, null);
}

void clear_model_xpx(Model * pmod)
{
    if (pmod.xpx != null) {
        free(pmod.xpx);
        pmod.xpx = null;
    }
}

/**
 * gretl_model_free:
 * @pmod: pointer to #MODEL.
 *
 * Free allocated content of @pmod then the pointer itself.
 */

void gretl_model_free(Model * pmod)
{
    if (pmod != null) {
        clear_model(pmod);
        free(pmod);
    }
}

/**
 * gretl_model_free_on_exit:
 * @pmod: pointer to #MODEL.
 *
 * Free allocated content of @pmod then the pointer itself,
 * without regard to the model's reference count.
 */

void gretl_model_free_on_exit(Model * pmod)
{
    if (pmod != null) {
        remove_model_from_stack_on_exit(pmod);
        clear_model(pmod);
        free(pmod);
    }
}

private void copy_test(ModelTest * targ, const(ModelTest) * src)
{
    targ.type = src.type;

    if (src.param != null && *src.param != '\0') {
        targ.param = gretl_strdup(src.param);
    } else {
        targ.param = null;
    }

    targ.teststat = src.teststat;
    targ.dfn = src.dfn;
    targ.dfd = src.dfd;
    targ.order = src.order;
    targ.value = src.value;
    targ.pvalue = src.pvalue;
    targ.crit = src.crit;
    targ.alpha = src.alpha;
    targ.opt = src.opt;
}

private int copy_model_tests(Model * targ, const(Model) * src)
{
    int i, n = src.ntests;

    if (n <= 0 || src.tests == null) {
        return 0;
    }

    targ.tests = cast(ModelTest *) malloc(n * ModelTest.sizeof);
    if (targ.tests == null) return 1;

    for (i=0; i<n; i++) {
        copy_test(&targ.tests[i], &src.tests[i]);
    }

    return 0;
}

private int real_add_test_to_model(Model * pmod, ModelTest * test)
{
    ModelTest * tests = null;
    int nt = pmod.ntests;
    int err = 0;

    tests = cast(ModelTest *) realloc(pmod.tests, (nt + 1) * ModelTest.sizeof);

    if (tests == null) {
        err = GretlError.alloc;
    } else {
        pmod.tests = tests;
        pmod.ntests += 1;
        copy_test(&pmod.tests[nt], test);
    }

    return err;
}

private int attach_model_params_from_xml(xmlNodePtr node, xmlDocPtr doc,
                                         Model * pmod)
{
    char ** S;
    int slop = (pmod.ci == GretlCmdIndex.garch);
    int np = 0;
    int err = 0;

    S = gretl_xml_get_strings_array(node, doc, &np, slop, &err);

    if (!err) {
        pmod.params = S;
        pmod.nparams = np;
    }

    return err;
}

int attach_model_tests_from_xml(Model * pmod, xmlNodePtr node)
{
    ModelTest test;
    xmlNodePtr cur = node.children;
    int got, err = 0;

    while (cur != null && !err) {
        gretl_test_init(&test, cast(ModelTestType) 0);
        got = 0;
        got += gretl_xml_get_prop_as_int(cur, "type", &test.type);
        got += gretl_xml_get_prop_as_uchar(cur, "teststat", &test.teststat);
        got += gretl_xml_get_prop_as_int(cur, "dfn", &test.dfn);
        got += gretl_xml_get_prop_as_double(cur, "dfd", &test.dfd);
        got += gretl_xml_get_prop_as_int(cur, "order", &test.order);
        got += gretl_xml_get_prop_as_double(cur, "value", &test.value);
        got += gretl_xml_get_prop_as_double(cur, "pvalue", &test.pvalue);
        got += gretl_xml_get_prop_as_double(cur, "crit", &test.crit);
        got += gretl_xml_get_prop_as_double(cur, "alpha", &test.alpha);
        if (got < 7) {
            err = GretlError.data;
        } else {
            gretl_xml_get_prop_as_opt(cur, "opt", &test.opt);
            gretl_xml_get_prop_as_string(cur, "param", &test.param);
            err = real_add_test_to_model(pmod, &test);
            free(test.param); /* copied by real_add_test_to_model */
        }
        cur = cur.next;
    }

    return err;
}

private void serialize_test(const(ModelTest) * src, Prn * prn)
{
    pprintf(prn, "<test type=\"%d\" ", src.type);

    if (src.param != null && *src.param != '\0') {
        pprintf(prn, "param=\"%s\" ", src.param);
    }

    pprintf(prn, "teststat=\"%d\" ", cast(int) src.teststat);
    pprintf(prn, "dfn=\"%d\" ", src.dfn);
    pprintf(prn, "dfd=\"%g\" ", src.dfd);
    pprintf(prn, "order=\"%d\" ", src.order);
    pprintf(prn, "value=\"%.15g\" ", src.value);
    pprintf(prn, "pvalue=\"%.15g\" ", src.pvalue);

    if (!na(src.crit)) {
        pprintf(prn, "crit=\"%g\" ", src.crit);
        pprintf(prn, "alpha=\"%g\" ", src.alpha);
    }

    if (src.opt != Gretlopt.OPT_NONE) {
        pprintf(prn, "opt=\"%u\" ", cast(guint32) src.opt);
    }

    pputs(prn, "/>\n");
}

private GretlBundle * bundlize_test(const(ModelTest) * src)
{
    GretlBundle * b = gretl_bundle_new();

    if (b == null) {
        return null;
    }

    if (src.param != null && *src.param != '\0') {
        gretl_bundle_set_string(b, "param", src.param);
    }

    if (src.dfn > 0) {
        gretl_bundle_set_scalar(b, "dfn", cast(double) src.dfn);
    }
    if (src.dfd > 0) {
        gretl_bundle_set_scalar(b, "dfd", src.dfd);
    }
    if (src.order > 0) {
        gretl_bundle_set_scalar(b, "order", src.order);
    }
    if (!na(src.value)) {
        gretl_bundle_set_scalar(b, "test", src.value);
    }
    if (!na(src.pvalue)) {
        gretl_bundle_set_scalar(b, "pvalue", src.pvalue);
    }
    if (!na(src.crit)) {
        gretl_bundle_set_scalar(b, "crit", src.crit);
        gretl_bundle_set_scalar(b, "alpha", src.alpha);
    }

    return b;
}

private int serialize_model_tests(const(Model) * pmod, Prn * prn)
{
    int i, n = pmod.ntests;

    if (n <= 0 || pmod.tests == null) {
        return 0;
    }

    pprintf(prn, "<tests count=\"%d\">\n", pmod.ntests);

    for (i=0; i<n; i++) {
        serialize_test(&pmod.tests[i], prn);
    }

    pputs(prn, "</tests>\n");

    return 0;
}

private int testvals_differ(double x, double y)
{
    double eq_tol = 1.0e-10;
    double reldiff;

    if (x == 0.0) {
        reldiff = fabs(y);
    } else if (y == 0.0) {
        reldiff = fabs(x);
    } else if (x > y) {
        reldiff = fabs((x - y) / y);
    } else {
        reldiff = fabs((y - x) / x);
    }

    return reldiff > eq_tol;
}

private int test_params_differ(const(char) * p, const(char) * s)
{
    int ret = 0;

    if (s == null && p != null) {
        ret = 1;
    } else if (p == null && s != null) {
        ret = 1;
    } else if (p != null && s != null) {
        ret = strcmp(p, s);
    }

    return ret;
}

private int model_tests_differ(ModelTest * mt1, ModelTest * mt2)
{
    int ret = 0;

    if (mt1.type != mt2.type) {
        ret = 1;
    } else if (mt1.order != mt2.order) {
        ret = 1;
    } else if (mt1.teststat != mt2.teststat) {
        ret = 1;
    } else if (test_params_differ(mt1.param, mt2.param)) {
        ret = 1;
    } else if (testvals_differ(mt1.value, mt2.value)) {
        ret = 1;
    }

    return ret;
}

/**
 * model_test_free:
 * @test: object to free.
 */

void model_test_free(ModelTest * test)
{
    free(test.param);
    free(test);
}

/**
 * model_test_new:
 * @ttype: type of test to add.
 *
 * Returns: new #ModelTest pointer, or %NULL on failure.
 */

ModelTest * model_test_new(ModelTestType ttype)
{
    ModelTest * test = cast(ModelTest *) malloc(ModelTest.sizeof);

    if (test != null) {
        gretl_test_init(test, ttype);
    }

    return test;
}

/**
 * maybe_add_test_to_model:
 * @pmod: pointer to model.
 * @test: model test to be added.
 *
 * Adds a #ModelTest to @pmod, if the test in question has
 * not already been performed and recorded.  Note that this
 * function takes care of freeing @test.
 *
 * Returns: 1 if the test was added, otherwise 0.
 */

int maybe_add_test_to_model(Model * pmod, ModelTest * test)
{
    int i, done = 0, add = 0;

    if (test == null || test.teststat == TestStats.none) {
        return 0;
    }

    for (i=0; i<pmod.ntests; i++) {
        if (!model_tests_differ(test, &pmod.tests[i])) {
            done = 1;
        }
    }

    if (!done) {
        int n = pmod.ntests + 1;
        ModelTest * tests;

        tests = cast(ModelTest *) realloc(pmod.tests, n * ModelTest.sizeof);
        if (tests != null) {
            pmod.tests = tests;
            copy_test(&pmod.tests[n-1], test);
            pmod.ntests += 1;
            add = 1;
        }
    }

    free(test.param);
    free(test);

    return add;
}

int gretl_model_get_normality_test(const(Model) * pmod, Prn * prn)
{
    const(ModelTest) * test = null;
    int i, err = 0;

    for (i=0; i<pmod.ntests; i++) {
        if (pmod.tests[i].type == ModelTestType.normal) {
            test = &pmod.tests[i];
            break;
        }
    }

    if (test == null) {
        err = GretlError.badstat;
    } else {
        record_test_result(test.value, test.pvalue);
        gretl_model_test_print(pmod, i, prn);
    }

    return err;
}

void model_test_set_teststat(ModelTest * test, ubyte ts)
{
    test.teststat = ts;
}

void model_test_set_order(ModelTest * test, int order)
{
    test.order = order;
}

void model_test_set_dfn(ModelTest * test, int df)
{
    test.dfn = df;
}

void model_test_set_dfd(ModelTest * test, double df)
{
    test.dfd = df;
}

void model_test_set_value(ModelTest * test, double val)
{
    test.value = val;
}

void model_test_set_pvalue(ModelTest * test, double pval)
{
    test.pvalue = pval;
}

void model_test_set_opt(ModelTest * test, Gretlopt opt)
{
    test.opt = opt;
}

void model_test_set_crit_and_alpha(ModelTest * test,
                                   double crit,
                                   double alpha)
{
    test.crit = crit;
    test.alpha = alpha;
}

void model_test_set_param(ModelTest * test, const(char) * s)
{
    test.param = gretl_strdup(s);
}

void model_test_set_allocated_param(ModelTest * test, char * s)
{
    test.param = s;
}

private struct TestStrings {
    ModelTestType ID;
    string descrip;
    string H0;
}

private immutable TestStrings[] tstrings = [
    TestStrings(ModelTestType.add,
                N_("Test for addition of variables"),
                N_("parameters are zero for the variables")),
    TestStrings(ModelTestType.arch,
                N_("Test for ARCH of order %s"),
                N_("no ARCH effect is present")),
    TestStrings(ModelTestType.autocorr,
                N_("LM test for autocorrelation up to order %s"),
                N_("no autocorrelation")),
    TestStrings(ModelTestType.chow,
                N_("Chow test for structural break at observation %s"),
                N_("no structural break")),
    TestStrings(ModelTestType.chowdum,
                N_("Chow test for structural difference with respect to %s"),
                N_("no structural difference")),
    TestStrings(ModelTestType.cusum,
                N_("CUSUM test for parameter stability"),
                N_("no change in parameters")),
    TestStrings(ModelTestType.qlr,
                N_("QLR test for structural break"),
                N_("no structural break")),
    TestStrings(ModelTestType.groupwise,
                N_("Likelihood ratio test for groupwise heteroskedasticity"),
                N_("the units have a common error variance")),
    TestStrings(ModelTestType.logs,
                N_("Non-linearity test (logs)"),
                N_("relationship is linear")),
    TestStrings(ModelTestType.normal,
                N_("Test for normality of residual"),
                N_("error is normally distributed")),
    TestStrings(ModelTestType.omit,
                N_("Test for omission of variables"),
                N_("parameters are zero for the variables")),
    TestStrings(ModelTestType.reset,
                N_("RESET test for specification"),
                N_("specification is adequate")),
    TestStrings(ModelTestType.squares,
                N_("Non-linearity test (squares)"),
                N_("relationship is linear")),
    TestStrings(ModelTestType.whites,
                N_("White's test for heteroskedasticity"),
                N_("heteroskedasticity not present")),
    TestStrings(ModelTestType.bp,
                N_("Breusch-Pagan test for heteroskedasticity"),
                N_("heteroskedasticity not present")),
    TestStrings(ModelTestType.sargan,
                N_("Sargan over-identification test"),
                N_("all instruments are valid")),
    TestStrings(ModelTestType.iv_hausman,
                N_("Hausman test"),
                N_("OLS estimates are consistent")),
    TestStrings(ModelTestType.panel_hausman,
                N_("Hausman test"),
                N_("GLS estimates are consistent")),
    TestStrings(ModelTestType.panel_f,
                N_("Test for differing group intercepts"),
                N_("The groups have a common intercept")),
    TestStrings(ModelTestType.panel_bp,
                N_("Breusch-Pagan test"),
                N_("Variance of the unit-specific error = 0")),
    TestStrings(ModelTestType.panel_timedum,
                N_("Wald joint test on time dummies"),
                N_("No time effects")),
    TestStrings(ModelTestType.het_1,
                N_("Pesaran-Taylor test for heteroskedasticity"),
                N_("heteroskedasticity not present")),
    TestStrings(ModelTestType.comfac,
                N_("Test of common factor restriction"),
                N_("restriction is acceptable")),
    TestStrings(ModelTestType.indep,
                N_("Test of independence"),
                N_("rho = 0")),
    TestStrings(ModelTestType.re,
                N_("LR test for rho = 0"),
                null),
    TestStrings(ModelTestType.within_f,
                N_("Joint test on named regressors"),
                null),
    TestStrings(ModelTestType.re_wald,
                N_("Joint test on named regressors"),
                null),
    TestStrings(ModelTestType.panel_welch,
                N_("Robust test for differing group intercepts"),
                N_("The groups have a common intercept")),
    TestStrings(ModelTestType.xdepend,
                N_("Pesaran CD test for cross-sectional dependence"),
                N_("No cross-sectional dependence")),
    TestStrings(ModelTestType.panel_ar,
                N_("Wooldridge test for autocorrelation in panel data"),
                N_("No first-order autocorrelation (rho = -0.5)")),
    TestStrings(ModelTestType.bds,
                N_("BDS test for nonlinearity"),
                null),
    TestStrings(ModelTestType.max, null, null)
];

private string test_get_descrip(const(ModelTest) * test)
{
    string dfree =
        N_("Distribution free Wald test for heteroskedasticity");
    string lr_indep =
        N_("Likelihood ratio test of independence");

    if (test.type == ModelTestType.groupwise &&
        test.teststat == TestStats.wald_chisq) {
        return dfree;
    } else if (test.type == ModelTestType.indep &&
               test.teststat == TestStats.lr) {
        return lr_indep;
    } else {
        int i;

        for (i=0; tstrings[i].ID < ModelTestType.max; i++) {
            if (test.type == tstrings[i].ID) {
                return tstrings[i].descrip;
            }
        }
    }

    return null;
}

private string test_get_opt_string(const(ModelTest) * test)
{
    if (test.type == ModelTestType.reset) {
        if (test.opt & Gretlopt.OPT_U) {
            return N_("squares only");
        } else if (test.opt & Gretlopt.OPT_C) {
            return N_("cubes only");
        }
    } else if (test.type == ModelTestType.whites) {
        if (test.opt & Gretlopt.OPT_X) {
            return N_("squares only");
        }
    } else if (test.type == ModelTestType.bp) {
        if (test.opt & Gretlopt.OPT_R) {
            return N_("robust variant");
        }
    }

    return null;
}

private void print_test_opt(const(ModelTest) * test, Prn * prn)
{
    string optstr = test_get_opt_string(test);

    if (optstr != null) {
        pprintf(prn, " (%s)", _(optstr).ptr);
    }
}

private int gretl_test_print_heading(const(ModelTest) * test, Prn * prn)
{
    string descrip;
    const(char) * param = null;
    char[16] ordstr;

    descrip = test_get_descrip(test);
    if (descrip == null) {
        return 1;
    }

    if (test.order > 0) {
        sprintf(ordstr.ptr, "%d", test.order);
        param = ordstr.ptr;
    } else if (test.type == ModelTestType.chow ||
               test.type == ModelTestType.chowdum) {
        param = test.param;
    }

    if (param != null) {
        if (plain_format(prn)) {
            pprintf(prn, _(descrip).ptr, param);
        } else {
            pprintf(prn, _(descrip).ptr, param);
        }
    } else {
        if (plain_format(prn)) {
            pputs(prn, _(descrip).ptr);
        } else {
            pputs(prn, _(descrip).ptr);
        }
    }

    if (test.opt) {
        print_test_opt(test, prn);
    }

    return 0;
}

private int print_add_omit_varnames(const(char) * s, Prn * prn)
{
    string[3] endings = [
        "\n    ",
        "\\\\\n \\qquad ",
        "\\par\n    "
    ];
    string ending = null;

    if (s == null || *s == '\0') {
        return 1;
    }

    if (plain_format(prn)) {
        ending = endings[0];
    } else if (tex_format(prn)) {
        ending = endings[1];
    } else if (rtf_format(prn)) {
        ending = endings[2];
    } else {
        return 1;
    }

    pputs(prn, ending.ptr);

    while (*s) {
        if (*s == ' ') {
            pputs(prn, ending.ptr);
        } else {
            pputc(prn, *s);
        }
        s++;
    }

    return 0;
}

const(char) * get_h0_string_for_test(ModelTestType ttype)
{
    string ret = null;
    int i;

    for (i=0; tstrings[i].ID < ModelTestType.max; i++) {
        if (ttype == tstrings[i].ID) {
            ret = tstrings[i].H0;
        }
    }

    return ret.ptr;
}

private void gretl_test_print_h_0(const(ModelTest) * test, int heading,
                                  Prn * prn)
{
    string H0 = null;
    int i;

    if (test.type == ModelTestType.panel_ar &&
        test.teststat == TestStats.student) {
        H0 = N_("No first-order autocorrelation (rho = 0)");
    } else {
        for (i=0; tstrings[i].ID < ModelTestType.max; i++) {
            if (test.type == tstrings[i].ID) {
                H0 = tstrings[i].H0;
            }
        }
    }

    if (test.type == ModelTestType.indep &&
        test.teststat == TestStats.wald_chisq &&
        H0 != null && !strcmp(H0.ptr, "rho = 0")) {
        /* adjust for QML biprobit */
        H0 = N_("atanh(rho) = 0");
    }

    if (H0 == null && test.type != ModelTestType.within_f &&
        test.type != ModelTestType.re_wald) {
        return;
    }

    if (heading) {
        pputs(prn, " -");
        if (tex_format(prn)) {
            pputc(prn, '-');
        }
    }

    if (H0 == null) {
        return;
    }

    if (plain_format(prn)) {
        pprintf(prn, "\n  %s: ", _("Null hypothesis").ptr);
        pputs(prn, _(H0).ptr);
    } else if (tex_format(prn)) {
        pprintf(prn, "\\\\\n\\quad %s: ", _("Null hypothesis").ptr);
        if (!strcmp(H0.ptr, "rho = 0")) {
            pputs(prn, "$\\rho = 0$");
        } else {
            pputs(prn, _(H0).ptr);
        }
    } else if (rtf_format(prn)) {
        pprintf(prn, "\\par\n %s: ", _("Null hypothesis").ptr);
        pputs(prn, _(H0).ptr);
    }

    if (test.type == ModelTestType.add || test.type == ModelTestType.omit) {
        print_add_omit_varnames(test.param, prn);
    }
}

private gchar * get_test_stat_string(const(ModelTest) * test, Prn * prn)
{
    int tex = tex_format(prn);
    gchar * ret = null;

    switch (test.teststat) {
    case TestStats.lm:
        ret = g_strdup_printf("LM = %g", test.value);
        break;
    case TestStats.f:
    case TestStats.reset:
        if (tex) {
            ret = g_strdup_printf("$F(%d, %g)$ = %g", test.dfn, test.dfd, test.value);
        } else {
            ret = g_strdup_printf("F(%d, %g) = %g", test.dfn, test.dfd, test.value);
        }
        break;
    case TestStats.sup_wald:
        if (tex) {
            ret = g_strdup_printf("max $\\chi^2(%d)$ = %g (%s)", test.dfn,
                                  test.value, test.param);
        } else {
            ret = g_strdup_printf(_("chi-square(%d) = %g at observation %s").ptr,
                                  test.dfn, test.value, test.param);
        }
        break;
    case TestStats.lmf:
        ret = g_strdup_printf("LMF = %g", test.value);
        break;
    case TestStats.wf:
        if (tex) {
            ret = g_strdup_printf("Welch $F(%d, %.1f)$ = %g", test.dfn, test.dfd, test.value);
        } else {
            ret = g_strdup_printf("Welch F(%d, %.1f) = %g", test.dfn, test.dfd, test.value);
        }
        break;
    case TestStats.harvey_collier:
        if (tex) {
            ret = g_strdup_printf("Harvey--Collier $t(%d)$ = %g", test.dfn, test.value);
        } else {
            ret = g_strdup_printf("Harvey-Collier t(%d) = %g", test.dfn, test.value);
        }
        break;
    case TestStats.normal_chisq:
        if (tex) {
            ret = g_strdup_printf("$\\chi^2(2)$ = %g", test.value);
        } else {
            ret = g_strdup_printf("%s(2) = %g", _("Chi-square").ptr, test.value);
        }
        break;
    case TestStats.lr:
    case TestStats.wald_chisq:
    case TestStats.lb_chisq:
        if (tex) {
            if (na(test.value)) {
                ret = g_strdup_printf("$\\chi^2(%d)$ = NA (%s)", test.dfn, _("failed").ptr);
            } else {
                ret = g_strdup_printf("$\\chi^2(%d)$ = %g", test.dfn, test.value);
            }
        } else {
            if (na(test.value)) {
                ret = g_strdup_printf("%s(%d) = NA (%s)", _("Chi-square").ptr, test.dfn,
                                      _("failed").ptr);
            } else {
                ret = g_strdup_printf("%s(%d) = %g", _("Chi-square").ptr, test.dfn, test.value);
            }
        }
        break;
    case TestStats.z:
        if (tex) {
            ret = g_strdup_printf("$z$ = %g", test.value);
        } else {
            ret = g_strdup_printf("z = %g", test.value);
        }
        break;
    case TestStats.student:
        if (tex) {
            ret = g_strdup_printf("$t$(%d) = %g", test.dfn, test.value);
        } else {
            ret = g_strdup_printf("t(%d) = %g", test.dfn, test.value);
        }
        break;
    default:
        break;
    }

    return ret;
}

private gchar * get_test_pval_string(const(ModelTest) * test, Prn * prn)
{
    int tex = tex_format(prn);
    gchar * ret = null;

    switch (test.teststat) {
    case TestStats.lm:
        if (tex) {
            ret = g_strdup_printf("$P$($\\chi^2(%d) >$ %g) = %g",
                                  test.dfn, test.value, test.pvalue);
        } else {
            ret = g_strdup_printf("P(%s(%d) > %g) = %g", _("Chi-square").ptr,
                                  test.dfn, test.value, test.pvalue);
        }
        break;
    case TestStats.f:
    case TestStats.reset:
        if (tex) {
            ret = g_strdup_printf("$P$($F(%d, %g) >$ %g) = %g",
                                  test.dfn, test.dfd, test.value, test.pvalue);
        } else {
            ret = g_strdup_printf("P(F(%d, %g) > %g) = %g",
                                  test.dfn, test.dfd, test.value, test.pvalue);
        }
        break;
    case TestStats.wf:
        if (tex) {
            ret = g_strdup_printf("$P$($F(%d, %.1f) >$ %g) = %g",
                                  test.dfn, test.dfd, test.value, test.pvalue);
        } else {
            ret = g_strdup_printf("P(F(%d, %.1f) > %g) = %g",
                                  test.dfn, test.dfd, test.value, test.pvalue);
        }
        break;
    case TestStats.lmf:
        if (tex) {
            ret = g_strdup_printf("$P$($F(%d, %g) >$ %g) = %g",
                                  test.dfn, test.dfd, test.value, test.pvalue);
        } else {
            ret = g_strdup_printf("P(F(%d, %g) > %g) = %g",
                                  test.dfn, test.dfd, test.value, test.pvalue);
        }
        break;
    case TestStats.harvey_collier:
        if (tex) {
            ret = g_strdup_printf("$P$($t_{%d} >$ %g) = %g",
                                  test.dfn, test.value, test.pvalue);
        } else {
            ret = g_strdup_printf("P(t(%d) > %g) = %g",
                                  test.dfn, test.value, test.pvalue);
        }
        break;
    case TestStats.student:
        if (tex) {
            ret = g_strdup_printf("$P$($|t| >$ %g) = %g",
                                  fabs(test.value), test.pvalue);
        } else {
            ret = g_strdup_printf("P(|t| > %g) = %g",
                                  fabs(test.value), test.pvalue);
        }
        break;
    case TestStats.normal_chisq:
    case TestStats.lr:
    case TestStats.wald_chisq:
    case TestStats.sup_wald:
    case TestStats.z:
        if (na(test.value)) {
            ; /* leave @ret NULL */
        } else if (na(test.pvalue)) {
            ret = g_strdup("NA");
        } else {
            ret = g_strdup_printf("%g", test.pvalue);
        }
        break;
    default:
        break;
    }

    return ret;
}

private void csv_print_test(const(ModelTest) * src, Prn * prn)
{
    string descrip = null;
    char c = prn_delim(prn);

    descrip = test_get_descrip(src);

    if (descrip != null) {
        string optstr = test_get_opt_string(src);

        if (optstr != null) {
            pprintf(prn, "\"%s, %s\"\n", descrip.ptr, optstr.ptr);
        } else {
            pprintf(prn, "\"%s\"\n", descrip.ptr);
        }
    }

    if (src.param != null && *src.param != '\0') {
        pprintf(prn, "\"%s\"%c\"%s\"\n", _("parameter").ptr, c, src.param);
    }

    if (src.dfn > 0 && src.dfd > 0) {
        pprintf(prn, "\"%s\"%c%d\n", _("dfn").ptr, c, src.dfn);
        pprintf(prn, "\"%s\"%c%g\n", _("dfd").ptr, c, src.dfd);
    } else if (src.dfn > 0) {
        pprintf(prn, "\"%s\"%c%d\n", _("df").ptr, c, src.dfn);
    }
    if (src.order) {
        pprintf(prn, "\"%s\"%c%d\n", _("lag order").ptr, c, src.order);
    }
    pprintf(prn, "\"%s\"%c%.15g\n", _("test statistic").ptr, c, src.value);
    if (!na(src.pvalue)) {
        pprintf(prn, "\"%s\"%c%.15g\n", _("p-value").ptr, c, src.pvalue);
    }

    if (!na(src.crit)) {
        double a = 100 * src.alpha;
        gchar * buf;

        buf = g_strdup_printf(_("%g percent critical value").ptr, a);
        pprintf(prn, "\"%s\"%c%g\n", buf, c, src.crit);
        g_free(buf);
    }

    pputc(prn, '\n');
}

private bool asy_test(int t)
{
    return (t == TestStats.wald_chisq ||
            t == TestStats.z);
}

private bool asy_pval(const(ModelTest) * test)
{
    return (test.teststat == TestStats.sup_wald ||
            (test.opt & Gretlopt.OPT_A));
}

private bool boot_pval(const(ModelTest) * test)
{
    return (test.opt & Gretlopt.OPT_B) != 0;
}

void gretl_model_test_print_direct(const(ModelTest) * test, int heading, Prn * prn)
{
    string tstr;
    gchar * buf = null;

    if (rtf_format(prn)) {
        pputs(prn, "\\par \\ql ");
    }

    if (heading) {
        gretl_test_print_heading(test, prn);
    }

    gretl_test_print_h_0(test, heading, prn);

    tstr = asy_test(test.teststat)? N_("Asymptotic test statistic") :
        N_("Test statistic");

    buf = get_test_stat_string(test, prn);

    if (plain_format(prn)) {
        pprintf(prn, "\n  %s: %s\n", _(tstr).ptr, buf);
    } else if (tex_format(prn)) {
        pprintf(prn, "\\\\\n\\quad %s: %s\\\\\n", _(tstr).ptr, buf);
    } else if (rtf_format(prn)) {
        pprintf(prn, "\\par\n %s: %s\\par\n", _(tstr).ptr, buf);
    }

    g_free(buf);
    buf = get_test_pval_string(test, prn);

    if (buf != null) {
        string pvstr =
            asy_pval(test) ? N_("with asymptotic p-value") :
            boot_pval(test) ? N_("with bootstrap p-value") :
            N_("with p-value");

        if (plain_format(prn)) {
            pprintf(prn, "  %s = %s\n\n", _(pvstr).ptr, buf);
        } else if (tex_format(prn)) {
            pprintf(prn, "\\quad %s = %s\\\\\n", _(pvstr).ptr, buf);
        } else if (rtf_format(prn)) {
            pprintf(prn, " %s = %s\\par\n\n", _(pvstr).ptr, buf);
        }
    } else if (!na(test.crit) && !na(test.alpha)) {
        double a = test.alpha * 100.0;

        if (plain_format(prn)) {
            buf = g_strdup_printf(_("%g percent critical value").ptr, a);
            pprintf(prn, "  (%s = %.2f)\n\n", buf, test.crit);
        } else if (tex_format(prn)) {
            buf = g_strdup_printf(_("%g percent critical value").ptr, a);
            pprintf(prn, "\\quad (%s = %.2f)\\\\\n", buf, test.crit);
        } else if (rtf_format(prn)) {
            buf = g_strdup_printf(_("%g percent critical value").ptr, a);
            pprintf(prn, " (%s = %.2f)\\par\n\n", buf, test.crit);
        }
    } else {
        pputc(prn, '\n');
    }

    g_free(buf);
}

void gretl_model_test_print(const(Model) * pmod, int i, Prn * prn)
{
    if (i >= 0 && i < pmod.ntests) {
        if (csv_format(prn)) {
            csv_print_test(&pmod.tests[i], prn);
        } else {
            gretl_model_test_print_direct(&pmod.tests[i], 1, prn);
        }
    }
}

void gretl_model_print_last_test(const(Model) * pmod, Prn * prn)
{
    gretl_model_test_print(pmod, pmod.ntests - 1, prn);
}

private ArInfo * copy_ar_info(const(ArInfo) * src)
{
    ArInfo * targ;

    targ = cast(ArInfo *) malloc(ArInfo.sizeof);
    if (targ == null) return null;

    if (src.arlist != null) {
        int m = src.arlist[0];

        if ((targ.rho = copyvec(src.rho, m)) == null) {
            free(targ);
            return null;
        }

        if ((targ.sderr = copyvec(src.sderr, m)) == null) {
            free(targ.rho);
            free(targ);
            return null;
        }

        targ.arlist = gretl_list_copy(src.arlist);
        if (targ.arlist == null) {
            free(targ.rho);
            free(targ.sderr);
            free(targ);
            return null;
        }
    }

    return targ;
}

private int copy_model_params(Model * targ, const(Model) * src)
{
    int err = 0;

    if (src.nparams > 0) {
        targ.params = strings_array_dup(src.params, src.nparams);

        if (targ.params == null) {
            err = GretlError.alloc;
        } else {
            targ.nparams = src.nparams;
        }
    }

    return err;
}

private void serialize_coeff_sep(ModelDataItem * item, Prn * prn)
{
    CoeffSep * cs = cast(CoeffSep *) item.ptr;

    pprintf(prn, " pos=\"%d\"", cs.pos);
    if (cs.str[0] != '\0') {
        pputs(prn, " string=\"");
        gretl_xml_put_string(cs.str.ptr, prn);
        pputc(prn, '"');
    }
    pputs(prn, "/>\n");
}

private void serialize_vcv_info(ModelDataItem * item, Prn * prn)
{
    VCVInfo * vi = cast(VCVInfo *) item.ptr;

    pprintf(prn, " vmaj=\"%d\"", vi.vmaj);
    pprintf(prn, " vmin=\"%d\"", vi.vmin);
    if (vi.order > 0) {
        pprintf(prn, " order=\"%d\"", vi.order);
    }
    if (vi.flags > 0) {
        pprintf(prn, " flags=\"%d\"", cast(int) vi.flags);
    }
    if (vi.bw > 0 && !na(vi.bw)) {
        pprintf(prn, " bw=\"%.12g\"", vi.bw);
    }
    if (vi.cv1 != null) {
        pprintf(prn, " cv1=\"%s\"", vi.cv1);
    }
    if (vi.cv2 != null) {
        pprintf(prn, " cv2=\"%s\"", vi.cv2);
    }
    pputs(prn, "/>\n");
}

/* FIXME updating and placement of these functions */

private struct TypeMapper {
    int type;
    string name;
    string compat;
}

private immutable TypeMapper[10] mapper = [
    TypeMapper(GretlType.int_,          "int",         "1"),
    TypeMapper(GretlType.list,          "list",        "2"),
    TypeMapper(GretlType.double_,       "double",      "3"),
    TypeMapper(GretlType.int_array,     "intarray",    "4"),
    TypeMapper(GretlType.double_array,  "doublearray", "5"),
    TypeMapper(GretlType.string_,       "string",      "6"),
    TypeMapper(GretlType.cmplx_array,   "cmplxarray",  "8"),
    TypeMapper(GretlType.struct_,       "struct" ,     "9"),
    TypeMapper(GretlType.matrix,        "matrix" ,    "10"),
    TypeMapper(GretlType.none,          "none",        "0")
];

private int type_from_type_string(const(char) * s)
{
    int i;

    if (isdigit(*s)) {
        /* backward compatibility */
        for (i=0; mapper[i].type != GretlType.none; i++) {
            if (*s == mapper[i].compat[0]) {
                return mapper[i].type;
            }
        }
        if (*s == '7') {
            /* note: was "chararray" */
            return GretlType.string_;
        }
    } else {
        for (i=0; mapper[i].type != GretlType.none; i++) {
            if (!strcmp(s, mapper[i].name.ptr)) {
                return mapper[i].type;
            }
        }
        if (!strcmp(s, "chararray")) {
            return GretlType.string_;
        }
    }

    return GretlType.none;
}

private const(char) * gretl_type_name(int t)
{
    int i;

    for (i=0; mapper[i].type != GretlType.none; i++) {
        if (t == mapper[i].type) {
            return mapper[i].name.ptr;
        }
    }

    return "unknown";
}

private void serialize_model_data_items(const(Model) * pmod, Prn * prn)
{
    ModelDataItem * item;
    int i, j, nelem;

    pprintf(prn, "<data-items count=\"%d\">\n", pmod.n_data_items);

    for (i=0; i<pmod.n_data_items; i++) {
        item = cast(ModelDataItem *) pmod.data_items[i];
        nelem = 0;

        pprintf(prn, "<data-item key=\"%s\" type=\"%s\"",
                item.key, gretl_type_name(item.type));

        if (!strcmp(item.key, "coeffsep")) {
            serialize_coeff_sep(item, prn);
            continue;
        } else if (!strcmp(item.key, "vcv_info")) {
            serialize_vcv_info(item, prn);
            continue;
        }

        if (item.type == GretlType.int_array) {
            nelem = cast(int) (item.size / int.sizeof);
        } else if (item.type == GretlType.double_array) {
            nelem = cast(int) (item.size / double.sizeof);
        } else if (item.type == GretlType.cmplx_array) {
            nelem = cast(int) (item.size / Cmplx.sizeof);
        }

        if (nelem > 0) {
            pprintf(prn, " count=\"%d\">\n", nelem);
        } else if (item.type == GretlType.string_) {
            pputc(prn, '>');
        } else {
            pputs(prn, ">\n");
        }

        if (item.type == GretlType.int_) {
            pprintf(prn, "%d", *cast(int *) item.ptr);
        } else if (item.type == GretlType.double_) {
            double x = *cast(double *) item.ptr;

            if (na(x)) {
                pputs(prn, "NA");
            } else {
                pprintf(prn, "%.15g", x);
            }
        } else if (item.type == GretlType.int_array) {
            int * vals = cast(int *) item.ptr;

            for (j=0; j<nelem; j++) {
                pprintf(prn, "%d ", vals[j]);
            }
        } else if (item.type == GretlType.double_array) {
            double * vals = cast(double *) item.ptr;

            for (j=0; j<nelem; j++) {
                if (na(vals[j])) {
                    pputs(prn, "NA ");
                } else {
                    pprintf(prn, "%.15g ", vals[j]);
                }
            }
        } else if (item.type == GretlType.cmplx_array) {
            Cmplx * vals = cast(Cmplx *) item.ptr;
            double x;

            for (j=0; j<nelem; j++) {
                x = vals[j].r;
                if (na(x)) {
                    pputs(prn, "NA ");
                } else {
                    pprintf(prn, "%.15g ", x);
                }
                x = vals[j].i;
                if (na(x)) {
                    pputs(prn, "NA ");
                } else {
                    pprintf(prn, "%.15g ", x);
                }
            }
        } else if (item.type == GretlType.list) {
            int * list = cast(int *) item.ptr;

            for (j=0; j<=list[0]; j++) {
                pprintf(prn, "%d ", list[j]);
            }
        } else if (item.type == GretlType.string_) {
            gretl_xml_put_string(cast(const(char) *) item.ptr, prn);
        } else if (item.type == GretlType.matrix) {
            GretlMatrixRaw * m = cast(GretlMatrixRaw *) item.ptr;

            gretl_matrix_serialize(m, null, prn);
        } else if (item.type == GretlType.array) {
            GretlArray * A = cast(GretlArray *) item.ptr;

            gretl_array_serialize(A, prn);
        } else {
            ; /* no-op: not handled */
        }

        pputs(prn, "</data-item>\n");
    }

    pputs(prn, "</data-items>\n");
}

private int copy_model_data_items(Model * targ, const(Model) * src)
{
    int i, n = src.n_data_items;
    int err = 0;

    targ.data_items = cast(ModelDataItem **) malloc(n * (ModelDataItem *).sizeof);
    if (targ.data_items == null) {
        return 1;
    }

    for (i=0; i<n; i++) {
        targ.data_items[i] = null;
    }

    for (i=0; i<n; i++) {
        targ.data_items[i] = replicate_data_item(src.data_items[i]);
        if (targ.data_items[i] == null) {
            err = 1;
            break;
        }
    }

    if (err) {
        for (i=0; i<n; i++) {
            free(targ.data_items[i]);
        }
        free(targ.data_items);
        targ.data_items = null;
        targ.n_data_items = 0;
    } else {
        targ.n_data_items = n;
    }

    return err;
}

/* ensure we get a valid identifier for the bundle key */

private char * item_key_to_bundle_key(char * targ,
                                      const(char) * src)
{
    strcpy(targ, src);
    return gretl_charsub(targ, '-', '_');
}

private GretlMatrixRaw * my_matrix_from_list(ModelDataItem * item)
{
    int * list = cast(int *) item.ptr;
    GretlMatrixRaw * m;
    int i;

    m = gretl_matrix_alloc(1, list[0]);

    if (m != null) {
        for (i=0; i<list[0]; i++) {
            m.val[i] = list[i+1];
        }
    }

    return m;
}

private GretlMatrixRaw * matrix_from_cmplx(ModelDataItem * item)
{
    int nr = cast(int) (item.size / Cmplx.sizeof);
    Cmplx * c = cast(Cmplx *) item.ptr;
    GretlMatrixRaw * m;
    int i;

    m = gretl_matrix_alloc(nr, 2);

    if (m != null) {
        for (i=0; i<nr; i++) {
            gretl_matrix_set(m, i, 0, c[i].r);
            gretl_matrix_set(m, i, 1, c[i].i);
        }
    }

    return m;
}

int bundlize_model_data_items(const(Model) * pmod, GretlBundle * b)
{
    ModelDataItem * item;
    char[VNAMELEN] bkey;
    GretlMatrixRaw * m;
    double xval;
    int i, ival;
    int err = 0;

    for (i=0; i<pmod.n_data_items && !err; i++) {
        item = cast(ModelDataItem *) pmod.data_items[i];
        item_key_to_bundle_key(bkey.ptr, item.key);
        if (gretl_bundle_has_key(b, bkey.ptr)) {
            continue; /* item already present */
        } else if (item.type == GretlType.int_) {
            ival = *cast(int *) item.ptr;
            err = gretl_bundle_set_scalar(b, bkey.ptr, ival);
        } else if (item.type == GretlType.double_) {
            xval = *cast(double *) item.ptr;
            err = gretl_bundle_set_scalar(b, bkey.ptr, xval);
        } else if (item.type == GretlType.matrix) {
            err = gretl_bundle_set_matrix(b, bkey.ptr, cast(const(GretlMatrixRaw) *) item.ptr);
        } else if (item.type == GretlType.list) {
            if (pmod.ci == GretlCmdIndex.midasreg && !strcmp(bkey.ptr, "seplist")) {
                ; /* internal use only, skip it */
            } else if (pmod.ci == GretlCmdIndex.arma && !strcmp(bkey.ptr, "ainfo")) {
                m = my_matrix_from_list(item);
                if (m != null) {
                    err = gretl_bundle_donate_data(b, bkey.ptr, m,
                                                   GretlType.matrix);
                }
            } else {
                err = gretl_bundle_set_list(b, bkey.ptr, cast(const(int) *) item.ptr);
            }
        } else if (item.type == GretlType.array) {
            gretl_bundle_set_data(b, bkey.ptr, item.ptr, cast(GretlType) item.type);
        } else if (item.type == GretlType.cmplx_array) {
            m = matrix_from_cmplx(item);
            if (m != null) {
                err = gretl_bundle_donate_data(b, bkey.ptr, m,
                                               GretlType.matrix);
            }
        }
    }

    if (!err && pmod.tests != null) {
        const(ModelTest) * test;
        const(char) * tkey;
        GretlBundle * bt;

        for (i=0; i<pmod.ntests && !err; i++) {
            test = &pmod.tests[i];
            tkey = test_type_key(cast(ModelTestType) test.type);
            if (tkey != null && !gretl_bundle_has_key(b, tkey)) {
                bt = bundlize_test(test);
                if (bt != null) {
                    err = gretl_bundle_donate_data(b, tkey, bt,
                                                   GretlType.bundle);
                }
            }
        }
    }

    if (!err) {
        /* estimation start and end, 1-based */
        gretl_bundle_set_int(b, "t1", pmod.t1 + 1);
        gretl_bundle_set_int(b, "t2", pmod.t2 + 1);
        /* contemporaneous sample range, 1-based */
        gretl_bundle_set_int(b, "sample_t1", pmod.smpl.t1 + 1);
        gretl_bundle_set_int(b, "sample_t2", pmod.smpl.t2 + 1);
    }

    if (pmod.ci == GretlCmdIndex.panel && (pmod.opt & Gretlopt.OPT_F) &&
        !gretl_bundle_has_key(b, "within_rsq")) {
        /* fixed effects: add within R-squared */
        gretl_bundle_set_scalar(b, "within_rsq", pmod.adjrsq);
    }

    if (!na(pmod.dw) && !gretl_bundle_has_key(b, "dw")) {
        /* add Durbin-Watson statistic if available */
        gretl_bundle_set_scalar(b, "dw", pmod.dw);
    }

    gretl_bundle_set_int(b, "asymptotic", ASYMPTOTIC_MODEL(cast(GretlCmdIndex) pmod.ci));

    return err;
}

void display_model_data_items(const(Model) * pmod)
{
    ModelDataItem * item;
    int i, n = pmod.n_data_items;

    fprintf(stderr, "model has %d data items\n", n);

    for (i=0; i<n; i++) {
        item = cast(ModelDataItem *) pmod.data_items[i];
        fprintf(stderr, "%d '%s': ", i, item.key);
        if (item.type == GretlType.int_) {
            fprintf(stderr, "%d\n", *cast(int *) item.ptr);
        } else if (item.type == GretlType.double_) {
            fprintf(stderr, "%.15g\n", *cast(double *) item.ptr);
        } else {
            fprintf(stderr, "%p\n", cast(void *) item.ptr);
        }
    }
}

private int copy_model(Model * targ, Model * src)
{
    int k = src.ncoeff;
    int m = k * (k + 1) / 2;
    int err = 0;

    /* monolithic copy of structure */
    *targ = *src;

    /* temporarily zero various count members just in case
       copying of allocated info fails */
    targ.ntests = 0;
    targ.nparams = 0;
    targ.n_data_items = 0;

    /* now work on pointer members */
    gretl_model_init_pointers(targ);

    if (src.coeff != null &&
        (targ.coeff = copyvec(src.coeff, src.ncoeff)) == null) {
        return 1;
    }
    if (src.sderr != null &&
        (targ.sderr = copyvec(src.sderr, src.ncoeff)) == null) {
        return 1;
    }
    if (src.uhat != null &&
        (targ.uhat = copyvec(src.uhat, src.full_n)) == null) {
        return 1;
    }
    if (src.yhat != null &&
        (targ.yhat = copyvec(src.yhat, src.full_n)) == null) {
        return 1;
    }
    if (src.submask != null &&
        (targ.submask = copy_subsample_mask(src.submask, &err)) == null) {
        return 1;
    }
    if (src.missmask != null &&
        (targ.missmask = gretl_strdup(src.missmask)) == null) {
        return 1;
    }
    if (src.xpx != null && (targ.xpx = copyvec(src.xpx, m)) == null) {
        return 1;
    }
    if (src.vcv != null && (targ.vcv = copyvec(src.vcv, m)) == null) {
        return 1;
    }
    if (src.arinfo != null &&
        (targ.arinfo = copy_ar_info(src.arinfo)) == null) {
        return 1;
    }
    if (src.ntests > 0 && src.tests != null) {
        copy_model_tests(targ, src);
        if (targ.tests == null) {
            return 1;
        }
        targ.ntests = src.ntests;
    }
    if (src.name != null) {
        targ.name = gretl_strdup(src.name);
        if (targ.name == null) {
            return 1;
        }
    }
    if (src.depvar != null) {
        targ.depvar = gretl_strdup(src.depvar);
        if (targ.depvar == null) {
            return 1;
        }
    }
    if (src.nparams > 0 && src.params != null) {
        err = copy_model_params(targ, src);
        if (err) {
            return err;
        }
    }
    if (src.n_data_items > 0) {
        copy_model_data_items(targ, src);
        if (targ.data_items == null) {
            return 1;
        }
        targ.n_data_items = src.n_data_items;
    }
    if (src.list != null &&
        (targ.list = gretl_list_copy(src.list)) == null) {
        return 1;
    }
    if (src.xlist != null &&
        (targ.xlist = gretl_list_copy(src.xlist)) == null) {
        return 1;
    }
    if (gretl_is_between_model(src) && src.dataset != null) {
        /* special for "between" panel model: transfer the
           group-means dataset to @targ
        */
        targ.dataset = src.dataset;
        src.dataset = null;
    }

    return 0;
}

int gretl_model_serialize(const(Model) * pmod, SavedObjectFlags flags,
                          Prn * prn)
{
    char * xmlname = null;
    int k = pmod.ncoeff;
    int m = k * (k + 1) / 2;
    int err = 0;

    if (pmod.name == null || *pmod.name == '\0') {
        xmlname = gretl_strdup("none");
    } else {
        xmlname = gretl_xml_encode(pmod.name);
    }

    pprintf(prn, "<gretl-model ID=\"%d\" name=\"%s\" saveflags=\"%d\" ",
            pmod.ID, xmlname, cast(int) flags);
    free(xmlname);

    if (pmod.depvar != null) {
        /* In an nls model pmod->depvar may contain '<' or '>' */
        if (gretl_xml_validate(pmod.depvar)) {
            pprintf(prn, "depvar=\"%s\" ", pmod.depvar);
        } else {
            char * tmp = gretl_xml_encode(pmod.depvar);

            pprintf(prn, "depvar=\"%s\" ", tmp);
            free(tmp);
        }
    }

    pprintf(prn, "t1=\"%d\" t2=\"%d\" nobs=\"%d\" ",
            pmod.t1, pmod.t2, pmod.nobs);
    pprintf(prn, "full_n=\"%d\" ncoeff=\"%d\" dfn=\"%d\" dfd=\"%d\" ",
            pmod.full_n, pmod.ncoeff, pmod.dfn, pmod.dfd);
    pprintf(prn, "ifc=\"%d\" ci=\"%s\" opt=\"%u\" nwt=\"%d\" aux=\"%d\" ",
            pmod.ifc, gretl_command_word(pmod.ci), cast(uint) pmod.opt,
            pmod.nwt, pmod.aux);

    gretl_push_c_numeric_locale();

    gretl_xml_put_double("ess", pmod.ess, prn);
    gretl_xml_put_double("tss", pmod.tss, prn);
    gretl_xml_put_double("sigma", pmod.sigma, prn);
    gretl_xml_put_double("rsq", pmod.rsq, prn);
    gretl_xml_put_double("adjrsq", pmod.adjrsq, prn);
    gretl_xml_put_double("fstt", pmod.fstt, prn);
    gretl_xml_put_double("chi-square", pmod.chisq, prn);
    gretl_xml_put_double("lnL", pmod.lnL, prn);
    gretl_xml_put_double("ybar", pmod.ybar, prn);
    gretl_xml_put_double("sdy", pmod.sdy, prn);

    gretl_xml_put_double("crit0", pmod.criterion[0], prn);
    gretl_xml_put_double("crit1", pmod.criterion[1], prn);
    gretl_xml_put_double("crit2", pmod.criterion[2], prn);

    gretl_xml_put_double("dw", pmod.dw, prn);
    gretl_xml_put_double("rho", pmod.rho, prn);

    pputs(prn, ">\n");

    if (pmod.smpl.rseed > 0) {
        pprintf(prn, "<sample t1=\"%d\" t2=\"%d\" rseed=\"%lu\"/>\n",
                pmod.smpl.t1, pmod.smpl.t2, pmod.smpl.rseed);
    } else {
        pprintf(prn, "<sample t1=\"%d\" t2=\"%d\"/>\n",
                pmod.smpl.t1, pmod.smpl.t2);
    }

    gretl_xml_put_double_array("coeff", pmod.coeff, k, prn);
    gretl_xml_put_double_array("sderr", pmod.sderr, k, prn);

    if (pmod.uhat != null) {
        gretl_xml_put_double_array("uhat", pmod.uhat, pmod.full_n, prn);
    }

    if (pmod.yhat != null) {
        gretl_xml_put_double_array("yhat", pmod.yhat, pmod.full_n, prn);
    }

    if (pmod.submask != null) {
        write_model_submask(pmod, prn);
    }

    if (pmod.missmask != null) {
        pputs(prn, "<missmask>");
        pputs(prn, pmod.missmask);
        pputs(prn, "</missmask>\n");
    }

    if (pmod.xpx != null) {
        gretl_xml_put_double_array("xpx", pmod.xpx, m, prn);
    }
    if (pmod.vcv != null) {
        gretl_xml_put_double_array("vcv", pmod.vcv, m, prn);
    }

    if (pmod.arinfo != null && pmod.arinfo.arlist != null) {
        int r = pmod.arinfo.arlist[0];

        pputs(prn, "<arinfo>\n");
        gretl_xml_put_tagged_list("arlist", pmod.arinfo.arlist, prn);
        gretl_xml_put_double_array("rho", pmod.arinfo.rho, r, prn);
        gretl_xml_put_double_array("sderr", pmod.arinfo.sderr, r, prn);
        pputs(prn, "</arinfo>\n");
    }

    if (pmod.ntests > 0 && pmod.tests != null) {
        serialize_model_tests(pmod, prn);
    }

    if (pmod.nparams > 0 && pmod.params != null) {
        gretl_xml_put_strings_array("params",
                                    cast(const(char) **) pmod.params,
                                    pmod.nparams, prn);
    }

    if (pmod.list != null) {
        gretl_xml_put_tagged_list("list", pmod.list, prn);
    }
    if (pmod.xlist != null) {
        gretl_xml_put_tagged_list("xlist", pmod.xlist, prn);
    }

    if (pmod.n_data_items > 0) {
        serialize_model_data_items(pmod, prn);
    }

    pputs(prn, "</gretl-model>\n");

    /* note: the DATASET element of pmod is not
       handled here */

    gretl_pop_c_numeric_locale();

    return err;
}

/* next block: functions for reconstituting a model from
   its XML representation */

private int model_submask_from_xml(xmlNodePtr node, xmlDocPtr doc,
                                   Model * pmod)
{
    char * mask;
    int err;

    err = gretl_xml_get_submask(node, doc, &mask);
    if (!err) {
        pmod.submask = mask;
    }

    return err;
}

private int retrieve_model_coeff_separator(xmlNodePtr cur, Model * pmod)
{
    CoeffSep * cs = cast(CoeffSep *) malloc(CoeffSep.sizeof);
    char * tmp = null;
    int err = 0;

    if (cs != null) {
        cs.str[0] = '\0';
        gretl_xml_get_prop_as_int(cur, "pos", &cs.pos);
        gretl_xml_get_prop_as_string(cur, "string", &tmp);
        if (tmp != null) {
            strcpy(cs.str.ptr, tmp);
            free(tmp);
        }
        err = gretl_model_set_data(pmod, "coeffsep", cs,
                                   GretlType.struct_,
                                   CoeffSep.sizeof);
    }

    return err;
}

private int retrieve_model_vcv_info(xmlNodePtr cur, Model * pmod)
{
    VCVInfo * vi = vcv_info_new();
    int err = 0;

    if (vi != null) {
        char * s = null;
        int ival;
        double x;

        gretl_xml_get_prop_as_int(cur, "vmaj", &vi.vmaj);
        gretl_xml_get_prop_as_int(cur, "vmin", &vi.vmin);
        if (gretl_xml_get_prop_as_int(cur, "order", &ival)) {
            vi.order = ival;
        }
        if (gretl_xml_get_prop_as_int(cur, "flags", &ival)) {
            vi.flags = cast(VCVFlags) ival;
        }
        if (gretl_xml_get_prop_as_double(cur, "bw", &x)) {
            vi.bw = x;
        }
        if (gretl_xml_get_prop_as_string(cur, "cv1", &s)) {
            vi.cv1 = s;
        }
        if (gretl_xml_get_prop_as_string(cur, "cv2", &s)) {
            vi.cv2 = s;
        }
        err = gretl_model_set_data_with_destructor(pmod, "vcv_info", vi,
                                                   GretlType.struct_,
                                                   VCVInfo.sizeof,
                                                   &vcv_info_free);
    }

    return err;
}

/* backward compatibility: transcribe old-style "gretl_set_int"
   settings to bits in pmod->opt, where applicable
*/

private int gretl_model_set_int_compat(Model * pmod,
                                       const(char) * key,
                                       int ival)
{
    Gretlopt opt = pmod.opt;
    int err = 0;

    if (ival == 0) {
        return 0;
    }

    if (pmod.ci == GretlCmdIndex.panel) {
        if (!strcmp(key, "between")) {
            pmod.opt |= Gretlopt.OPT_B;
        } else if (!strcmp(key, "fixed-effects")) {
            pmod.opt |= Gretlopt.OPT_F;
        } else if (!strcmp(key, "random-effects")) {
            pmod.opt |= Gretlopt.OPT_U;
        } else if (!strcmp(key, "unit-weights")) {
            pmod.opt |= Gretlopt.OPT_H;
        }
    } else if (pmod.ci == GretlCmdIndex.ar1) {
        if (!strcmp(key, "hilu")) {
            pmod.opt |= Gretlopt.OPT_H;
        } else if (!strcmp(key, "pwe")) {
            pmod.opt |= Gretlopt.OPT_P;
        }
    } else if (pmod.ci == GretlCmdIndex.gmm) {
        if (!strcmp(key, "iterated")) {
            pmod.opt |= Gretlopt.OPT_I;
        } else if (!strcmp(key, "two-step")) {
            pmod.opt |= Gretlopt.OPT_T;
        }
    } else if (pmod.ci == GretlCmdIndex.heckit) {
        if (!strcmp(key, "two-step")) {
            pmod.opt |= Gretlopt.OPT_T;
        }
    } else if (pmod.ci == GretlCmdIndex.logit || pmod.ci == GretlCmdIndex.probit) {
        if (!strcmp(key, "show-pvals")) {
            pmod.opt |= Gretlopt.OPT_P;
        }
    }

    if (!strcmp(key, "robust")) {
        pmod.opt |= Gretlopt.OPT_R;
    } else if (!strcmp(key, "time-dummies")) {
        pmod.opt |= Gretlopt.OPT_D;
    } else if (!strcmp(key, "no-df-corr")) {
        pmod.opt |= Gretlopt.OPT_N;
    }

    if (pmod.opt == opt) {
        /* not handled above */
        err = gretl_model_set_int(pmod, key, ival);
    }

    return err;
}

/* backward compat: convert ad hoc settings to VCVInfo */

private void compat_compose_vcv_info(Model * pmod)
{
    VCVInfo vi = VCVInfo(0, 0, 0, cast(VCVFlags) 0, NADBL, null, null);
    int ival;

    if (gretl_model_get_int(pmod, "hc")) {
        vi.vmaj = VCVMajorType.hc;
        vi.vmin = gretl_model_get_int(pmod, "hc_version");
    } else if ((ival = gretl_model_get_int(pmod, "ml_vcv")) != 0) {
        vi.vmaj = VCVMajorType.ml;
        vi.vmin = ival;
    } else if (gretl_model_get_int(pmod, "panel_hac")) {
        vi.vmaj = VCVMajorType.panel;
        vi.vmin = PanelVCVType.hac;
    } else if (gretl_model_get_int(pmod, "panel_bk")) {
        vi.vmaj = VCVMajorType.panel;
        vi.vmin = PanelVCVType.bk;
    } else if (gretl_model_get_int(pmod, "using_hac") ||
               gretl_model_get_int(pmod, "hac_kernel") ||
               gretl_model_get_int(pmod, "hac_lag")) {
        vi.vmaj = VCVMajorType.hac;
        vi.vmin = gretl_model_get_int(pmod, "hac_kernel");
        vi.order = gretl_model_get_int(pmod, "hac_lag");
        vi.flags = cast(VCVFlags) gretl_model_get_int(pmod, "hac_prewhiten");
        if (vi.vmin == HACKernel.qs) {
            vi.bw = gretl_model_get_double(pmod, "qs_bandwidth");
        }
    } else if (pmod.ci == GretlCmdIndex.lad && gretl_model_get_int(pmod, "rq")) {
        double a = gretl_model_get_double(pmod, "rq_alpha");

        if (na(a)) {
            /* doing VCV, not confidence intervals */
            vi.vmaj = VCVMajorType.rq;
            vi.vmin = (gretl_model_get_int(pmod, "rq_nid"))?
                RQVCVType.nid : RQVCVType.asy;
        }
    }

    if (vi.vmaj > 0) {
        VCVInfo * pvi = vcv_info_new();

        if (pvi != null) {
            *pvi = vi;
            gretl_model_set_data_with_destructor(pmod, "vcv_info", pvi,
                                                 GretlType.struct_,
                                                 VCVInfo.sizeof,
                                                 &vcv_info_free);
        }
    }
}

/* backward compat: convert old list separator */

private bool old_listsep(int v)
{
    return (v == 999 || v == 9999);
}

private void maybe_convert_listsep(int * list, const(DatasetRaw) * dset)
{
    int i;

    for (i=1; i<=list[0]; i++) {
        if (old_listsep(list[i]) && list[i] >= dset.v) {
            list[i] = LISTSEP;
        }
    }
}

private int model_data_items_from_xml(xmlNodePtr node, xmlDocPtr doc,
                                      Model * pmod, int fixopt)
{
    xmlNodePtr cur;
    int n_items;
    int got_vcv = 0;
    int err = 0;

    if (!gretl_xml_get_prop_as_int(node, "count", &n_items)) {
        return 1;
    }

    cur = node.children;

    while (cur != null && !err) {
        char * key = null;
        char * typestr = null;
        int ignore_err = 0;
        int t, nelem = 0;

        if (!gretl_xml_get_prop_as_string(cur, "type", &typestr) ||
            !gretl_xml_get_prop_as_string(cur, "key", &key)) {
            err = 1;
            break;
        }

        t = type_from_type_string(typestr);

        if (!strcmp(key, "roots")) {
            /* could be a backward-compat problem */
            ignore_err = 1;
        }

        if (!strcmp(key, "coeffsep") || !strcmp(key, "10")) {
            /* special, with backward compatibility */
            err = retrieve_model_coeff_separator(cur, pmod);
        } else if (!strcmp(key, "vcv_info")) {
            err = retrieve_model_vcv_info(cur, pmod);
            if (!err) {
                got_vcv = 1;
            }
        } else if (t == GretlType.int_) {
            int ival;

            if (!gretl_xml_node_get_int(cur, doc, &ival)) {
                err = 1;
            } else if (fixopt) {
                err = gretl_model_set_int_compat(pmod, key, ival);
            } else {
                err = gretl_model_set_int(pmod, key, ival);
            }
        } else if (t == GretlType.double_) {
            double xval;

            if (!gretl_xml_node_get_double(cur, doc, &xval)) {
                err = 1;
            } else {
                err = gretl_model_set_double(pmod, key, xval);
            }
        } else if (t == GretlType.list) {
            int * list = gretl_xml_get_list(cur, doc, &err);

            if (!err && list != null) {
                err = gretl_model_set_list_as_data(pmod, key, list);
            }
        } else if (t == GretlType.string_) {
            char * s;
            int gots;

            if (!strcmp(key, "pmask") || !strcmp(key, "qmask")) {
                /* these masks must not have leading white space */
                gots = gretl_xml_node_get_trimmed_string(cur, doc, &s);
            } else {
                gots = gretl_xml_node_get_string(cur, doc, &s);
            }
            if (!gots) {
                err = 1;
            } else {
                err = gretl_model_set_string_as_data(pmod, key, s);
            }
        } else if (t == GretlType.int_array) {
            int * ivals = gretl_xml_get_int_array(cur, doc, &nelem, &err);

            if (nelem > 0) {
                err = gretl_model_set_data(pmod, key, ivals, cast(GretlType) t,
                                           nelem * int.sizeof);
            }
        } else if (t == GretlType.double_array) {
            double * xvals = gretl_xml_get_double_array(cur, doc, &nelem, &err);

            if (nelem > 0) {
                err = gretl_model_set_data(pmod, key, xvals, cast(GretlType) t,
                                           nelem * double.sizeof);
            }
        } else if (t == GretlType.cmplx_array) {
            Cmplx * cvals = gretl_xml_get_cmplx_array(cur, doc, &nelem, &err);

            if (nelem > 0) {
                err = gretl_model_set_data(pmod, key, cvals, cast(GretlType) t,
                                           nelem * Cmplx.sizeof);
            }
        } else if (t == GretlType.matrix) {
            xmlNodePtr child = cur.children;
            GretlMatrixRaw * m;

            if (child == null) {
                err = GretlError.data;
            } else {
                m = gretl_xml_get_matrix(child, doc, &err);
                if (!err) {
                    err = gretl_model_set_matrix_as_data(pmod, key, m);
                }
            }
        }

        if (err && ignore_err) {
            fprintf(stderr, "ignoring error reloading %s\n", key);
            err = 0;
        }

        xmlFree(key);
        xmlFree(typestr);

        cur = cur.next;
    }

    if (!err && !got_vcv) {
        compat_compose_vcv_info(pmod);
    }

    return err;
}

private const(xmlChar) * XUC(const(char) * s)
{
    return cast(const(xmlChar) *) s;
}

private int arinfo_from_xml(xmlNodePtr node, xmlDocPtr doc,
                            Model * pmod)
{
    xmlNodePtr cur;
    int n, err = 0;

    if (gretl_model_add_arinfo(pmod, 0)) {
        return 1;
    }

    cur = node.children;

    while (cur != null && !err) {
        if (!xmlStrcmp(cur.name, XUC("arlist"))) {
            pmod.arinfo.arlist = gretl_xml_get_list(cur, doc, &err);
        } else if (!xmlStrcmp(cur.name, XUC("rho"))) {
            pmod.arinfo.rho = gretl_xml_get_double_array(cur, doc, &n, &err);
        } else if (!xmlStrcmp(cur.name, XUC("sderr"))) {
            pmod.arinfo.sderr = gretl_xml_get_double_array(cur, doc, &n, &err);
        }
        cur = cur.next;
    }

    return err;
}

/**
 * gretl_model_from_XML:
 * @node: XML node from which to read.
 * @doc: pointer to XML document.
 * @dset: dataset information.
 * @err: location to receive error code.
 *
 * Reads info on a gretl model from the given XML @node
 * and @doc, and reconstitutes the model in memory.
 *
 * Returns: allocated model, or %NULL on failure.
 */

Model * gretl_model_from_XML(xmlNodePtr node, xmlDocPtr doc,
                             const(DatasetRaw) * dset,
                             int * err)
{
    Model * pmod;
    char * buf = null;
    char * modname = null;
    xmlNodePtr cur;
    int fixopt = 0;
    int n, got = 0;

    pmod = gretl_model_new();
    if (pmod == null) {
        *err = GretlError.alloc;
        return null;
    }

    got += gretl_xml_get_prop_as_int(node, "ID", &pmod.ID);
    got += gretl_xml_get_prop_as_int(node, "t1", &pmod.t1);
    got += gretl_xml_get_prop_as_int(node, "t2", &pmod.t2);
    got += gretl_xml_get_prop_as_int(node, "nobs", &pmod.nobs);
    got += gretl_xml_get_prop_as_int(node, "full_n", &pmod.full_n);
    got += gretl_xml_get_prop_as_int(node, "ncoeff", &pmod.ncoeff);
    got += gretl_xml_get_prop_as_int(node, "dfn", &pmod.dfn);
    got += gretl_xml_get_prop_as_int(node, "dfd", &pmod.dfd);
    got += gretl_xml_get_prop_as_int(node, "ifc", &pmod.ifc);

    got += gretl_xml_get_prop_as_string(node, "ci", &buf);

    if (got < 10) {
        fprintf(stderr, "gretl_model_from_XML: got(1) = %d (expected 12)\n", got);
        if (buf != null) {
            free(buf);
        }
        *err = GretlError.data;
        goto bailout;
    }

    gretl_xml_get_prop_as_int(node, "nwt", &pmod.nwt);
    gretl_xml_get_prop_as_int(node, "aux", &pmod.aux);

    if (gretl_xml_get_prop_as_int(node, "opt", &n)) {
        pmod.opt = cast(Gretlopt) cast(guint32) n;
    } else {
        /* backward compat may be needed */
        fixopt = 1;
    }

    if (isdigit(*buf)) {
        /* backward compatibility: the model command is given
           as a number, not a string */
        pmod.ci = atoi(buf);
    } else {
        pmod.ci = gretl_command_number(buf);
    }

    free(buf);

    gretl_push_c_numeric_locale();

    gretl_xml_get_prop_as_double(node, "ess", &pmod.ess);
    gretl_xml_get_prop_as_double(node, "tss", &pmod.tss);
    gretl_xml_get_prop_as_double(node, "sigma", &pmod.sigma);
    gretl_xml_get_prop_as_double(node, "rsq", &pmod.rsq);
    gretl_xml_get_prop_as_double(node, "adjrsq", &pmod.adjrsq);
    gretl_xml_get_prop_as_double(node, "fstt", &pmod.fstt);
    gretl_xml_get_prop_as_double(node, "chi-square", &pmod.chisq);
    gretl_xml_get_prop_as_double(node, "lnL", &pmod.lnL);
    gretl_xml_get_prop_as_double(node, "ybar", &pmod.ybar);
    gretl_xml_get_prop_as_double(node, "sdy", &pmod.sdy);

    gretl_xml_get_prop_as_double(node, "crit0", &pmod.criterion[0]);
    gretl_xml_get_prop_as_double(node, "crit1", &pmod.criterion[1]);
    gretl_xml_get_prop_as_double(node, "crit2", &pmod.criterion[2]);

    gretl_xml_get_prop_as_double(node, "dw", &pmod.dw);
    gretl_xml_get_prop_as_double(node, "rho", &pmod.rho);

    gretl_xml_get_prop_as_string(node, "name", &modname);
    gretl_xml_get_prop_as_string(node, "depvar", &pmod.depvar);

    cur = node.children;

    while (cur != null) {
        if (!xmlStrcmp(cur.name, XUC("sample"))) {
            gretl_xml_get_prop_as_int(cur, "t1", &pmod.smpl.t1);
            gretl_xml_get_prop_as_int(cur, "t2", &pmod.smpl.t2);
            gretl_xml_get_prop_as_uint64(cur, "rseed", &pmod.smpl.rseed);
        } else if (!xmlStrcmp(cur.name, XUC("coeff"))) {
            pmod.coeff = gretl_xml_get_double_array(cur, doc, &n, err);
        } else if (!xmlStrcmp(cur.name, XUC("sderr"))) {
            pmod.sderr = gretl_xml_get_double_array(cur, doc, &n, err);
        } else if (!xmlStrcmp(cur.name, XUC("uhat"))) {
            pmod.uhat = gretl_xml_get_double_array(cur, doc, &n, err);
        } else if (!xmlStrcmp(cur.name, XUC("yhat"))) {
            pmod.yhat = gretl_xml_get_double_array(cur, doc, &n, err);
        } else if (!xmlStrcmp(cur.name, XUC("xpx"))) {
            /* note: allow this one to fail silently */
            pmod.xpx = gretl_xml_get_double_array(cur, doc, &n, null);
        } else if (!xmlStrcmp(cur.name, XUC("vcv"))) {
            pmod.vcv = gretl_xml_get_double_array(cur, doc, &n, err);
        } else if (!xmlStrcmp(cur.name, XUC("list"))) {
            pmod.list = gretl_xml_get_list(cur, doc, err);
        } else if (!xmlStrcmp(cur.name, XUC("xlist"))) {
            pmod.xlist = gretl_xml_get_list(cur, doc, err);
        } else if (!xmlStrcmp(cur.name, XUC("tests"))) {
            *err = attach_model_tests_from_xml(pmod, cur);
        } else if (!xmlStrcmp(cur.name, XUC("params"))) {
            *err = attach_model_params_from_xml(cur, doc, pmod);
        } else if (!xmlStrcmp(cur.name, XUC("submask"))) {
            fprintf(stderr, "getting submask...\n");
            *err = model_submask_from_xml(cur, doc, pmod);
            fprintf(stderr, "got it\n");
        } else if (!xmlStrcmp(cur.name, XUC("missmask"))) {
            if (!gretl_xml_node_get_string(cur, doc, &pmod.missmask)) {
                *err = 1;
            }
        } else if (!xmlStrcmp(cur.name, XUC("arinfo"))) {
            *err = arinfo_from_xml(cur, doc, pmod);
        } else if (!xmlStrcmp(cur.name, XUC("data-items"))) {
            *err = model_data_items_from_xml(cur, doc, pmod, fixopt);
        }

        if (*err) {
            fprintf(stderr, "gretl_model_from_XML: block 3: err = %d reading '%s'\n",
                    *err, cast(const(char) *) cur.name);
            break;
        }

        cur = cur.next;
    }

    gretl_pop_c_numeric_locale();

    /* backward compatibility */
    if (!*err && pmod.nparams > pmod.ncoeff && pmod.depvar == null) {
        int i, np = pmod.nparams - 1;

        pmod.depvar = pmod.params[0];
        for (i=0; i<np; i++) {
            pmod.params[i] = pmod.params[i+1];
        }
        pmod.params[np] = null;
        pmod.nparams = np;
    }

    if (!*err && pmod.opt != 0) {
        /* compat for some options whose flags have been changed */
        if (pmod.opt & Gretlopt.OPT_E) {
            if (pmod.ci == GretlCmdIndex.garch) {
                /* standardized residuals */
                pmod.opt &= ~Gretlopt.OPT_E;
                pmod.opt |= Gretlopt.OPT_Z;
            }
        }
        if (pmod.opt & Gretlopt.OPT_W) {
            if (pmod.ci == GretlCmdIndex.panel || pmod.ci == GretlCmdIndex.ivreg) {
                /* weights */
                pmod.opt &= ~Gretlopt.OPT_W;
                pmod.opt |= Gretlopt.OPT_H;
            } else if (pmod.ci == GretlCmdIndex.duration) {
                /* Weibull */
                pmod.opt &= ~Gretlopt.OPT_W;
                pmod.opt |= Gretlopt.OPT_B;
            }
        }
    }

    if (!*err && pmod.list != null) {
        maybe_convert_listsep(pmod.list, dset);
    }

    if (modname != null) {
        gretl_model_set_name(pmod, modname);
        free(modname);
    }

 bailout:

    if (*err) {
        if (pmod != null) {
            gretl_model_free(pmod);
            pmod = null;
        }
    }

    return pmod;
}

/**
 * gretl_model_copy:
 * @pmod: pointer to #MODEL to copy.
 *
 * Does a deep copy of @pmod: allocates a new #MODEL pointer
 * which has its own allocated copies of all the pointer
 * members of @pmod.  The only feature of @pmod that is
 * not duplicated is the reference count, which is set
 * to zero in the copy.
 *
 * Returns: the copied model, or %NULL on failure.
 */

Model * gretl_model_copy(Model * pmod)
{
    Model * new_ = cast(Model *) malloc(Model.sizeof);

    if (new_ != null) {
        int err = copy_model(new_, pmod);

        if (err) {
            clear_model(new_);
            free(new_);
            new_ = null;
        } else {
            /* nota bene */
            new_.refcount = 0;
        }
    }

    return new_;
}

/**
 * swap_models:
 * @targ: pointer to target #MODEL.
 * @src: pointer to source #MODEL.
 *
 * Swaps the content of the two model pointers.
 */

void swap_models(Model * targ, Model * src)
{
    Model tmp = *targ;

    *targ = *src;
    *src = tmp;
}

private bool normality_test(int ci, Gretlopt opt)
{
    return (ci == GretlCmdIndex.modtest && (opt & Gretlopt.OPT_N));
}

/**
 * command_ok_for_model:
 * @test_ci: index of command to be tested.
 * @opt: option for command to be tested.
 * @pmod: the gretl model in question.
 *
 * Check to see if the model-related command in question is
 * meaningful and acceptable in the context of the estimator
 * associated with @pmod. Note, though, that this function may
 * give a "false positive": to be quite sure, we may need to
 * know more about the model (e.g. specific options used).
 * See also model_test_ok().
 *
 * Returns: 1 if the command seems OK, otherwise 0.
 */

int command_ok_for_model(int test_ci, Gretlopt opt,
                         const(Model) * pmod)
{
    int between = 0;
    int regular_panel = 0;
    int mci = pmod.ci;
    int ok = 1;

    if (mci == GretlCmdIndex.midasreg) {
        /* treat as a case of NLS */
        mci = GretlCmdIndex.nls;
    }

    if (mci == GretlCmdIndex.nls && test_ci == GretlCmdIndex.fcast) {
        return 1;
    }

    if (mci == GretlCmdIndex.biprobit) {
        /* can we support anything else? */
        return test_ci == GretlCmdIndex.restrict;
    }

    if (test_ci == GretlCmdIndex.bkw) {
        /* most models should be OK */
        return pmod.ncoeff > 1 && (pmod.vcv != null || pmod.xpx != null);
    }

    if (NONLIST_MODEL(cast(GretlCmdIndex) mci)) {
        return (test_ci == GretlCmdIndex.restrict ||
                test_ci == GretlCmdIndex.tabprint ||
                (mci != GretlCmdIndex.mle && normality_test(test_ci, opt)));
    }

    between = gretl_is_between_model(pmod);
    regular_panel = gretl_is_regular_panel_model(pmod);

    switch (test_ci) {
    case GretlCmdIndex.add:
        if (mci == GretlCmdIndex.arma || mci == GretlCmdIndex.garch ||
            mci == GretlCmdIndex.heckit || mci == GretlCmdIndex.intreg) {
            ok = 0;
        } else if (mci == GretlCmdIndex.panel && (pmod.opt & Gretlopt.OPT_B)) {
            /* the panel "between" model */
            ok = 0;
        } else if (opt & (Gretlopt.OPT_L | Gretlopt.OPT_A)) {
            /* --lm and --auto variants: OLS only */
            ok = (mci == GretlCmdIndex.ols);
        }
        break;
    case GretlCmdIndex.omit:
        if (mci == GretlCmdIndex.arma || mci == GretlCmdIndex.garch ||
            mci == GretlCmdIndex.intreg || mci == GretlCmdIndex.dpanel) {
            /* FIXME INTREG should be workable? */
            ok = 0;
        } else if (between) {
            ok = 0;
        }
        break;
    case GretlCmdIndex.vif:
        if (mci == GretlCmdIndex.ivreg || mci == GretlCmdIndex.arma || mci == GretlCmdIndex.garch ||
            mci == GretlCmdIndex.panel || mci == GretlCmdIndex.dpanel) {
            ok = 0;
        }
        break;
    case GretlCmdIndex.eqnprint:
        if (mci == GretlCmdIndex.arma || mci == GretlCmdIndex.dpanel ||
            mci == GretlCmdIndex.heckit || mci == GretlCmdIndex.intreg) {
            ok = 0;
        }
        break;
    case GretlCmdIndex.modtest:
        if (opt & Gretlopt.OPT_H) {
            /* ARCH */
            ok = (mci != GretlCmdIndex.arch && mci != GretlCmdIndex.garch);
        } else if (opt & Gretlopt.OPT_C) {
            /* common factor restriction */
            ok = (mci == GretlCmdIndex.ar1);
        } else if (opt & Gretlopt.OPT_D) {
            /* x-sectional dependence */
            ok = !between;
        } else if (opt & Gretlopt.OPT_N) {
            /* normality */
            if (mci == GretlCmdIndex.logit || mci == GretlCmdIndex.heckit || mci == GretlCmdIndex.duration) {
                /* POISSON, NEGBIN? */
                ok = 0;
            }
        } else if (mci != GretlCmdIndex.ols) {
            if (mci == GretlCmdIndex.ivreg && (opt & (Gretlopt.OPT_A | Gretlopt.OPT_W))) {
                /* Autocorr. and H'sked. supported for IVREG */
                ok = 1;
            } else if (mci == GretlCmdIndex.arma && (opt & Gretlopt.OPT_A)) {
                ok = 1;
            } else if (mci == GretlCmdIndex.panel && (opt & Gretlopt.OPT_P)) {
                ok = 1;
            } else if (regular_panel && (opt & Gretlopt.OPT_A)) {
                ok = 1;
            } else {
                ok = 0;
            }
        }
        break;
    case GretlCmdIndex.chow:
    case GretlCmdIndex.cusum:
    case GretlCmdIndex.qlrtest:
    case GretlCmdIndex.leverage:
    case GretlCmdIndex.reset:
    case GretlCmdIndex.panspec:
        /* OLS-only tests */
        ok = (mci == GretlCmdIndex.ols);
        break;
    case GretlCmdIndex.restrict:
        if (mci == GretlCmdIndex.lad || mci == GretlCmdIndex.quantreg) {
            ok = 0;
        }
        break;
    default:
        break;
    }

    return ok;
}

/**
 * model_test_ok:
 * @ci: index of a model test command.
 * @opt: option associated with test command, if any.
 * @pmod: the model to be tested.
 * @dset: dataset information.
 *
 * A more rigorous superset of command_ok_for_model().  Use
 * this function if the extra information is available.
 *
 * Returns: 1 if the test command @ci (with possible option
 * @opt) is acceptable in the context of the model @pmod, and
 * the dataset described by @dset, otherwise 0.
 */

int model_test_ok(int ci, Gretlopt opt, const(Model) * pmod,
                  const(DatasetRaw) * dset)
{
    /* preliminary screening */
    int ok = command_ok_for_model(ci, opt, pmod);

    if (!ok) {
        return 0;
    }

    /* for now we'll treat MIDASREG as a case of NLS */
    if (ci == GretlCmdIndex.midasreg) {
        ci = GretlCmdIndex.nls;
    }

    if (ok && pmod.missmask != null) {
        /* can't do these with embedded missing obs? */
        if (gretl_is_regular_panel_model(pmod)) {
            ; /* OK? */
        } else if (ci == GretlCmdIndex.cusum || ci == GretlCmdIndex.bds ||
            (ci == GretlCmdIndex.modtest && (opt & (Gretlopt.OPT_A | Gretlopt.OPT_H)))) {
            ok = 0;
        }
    }

    if (ok && pmod.ncoeff == 1) {
        /* only one coefficient was estimated */
        if (ci == GretlCmdIndex.coeffsum) {
            ok = 0;
        } else if (pmod.ifc && ci == GretlCmdIndex.modtest) {
            /* const only: rule out squares, logs, h'sked */
            if (opt & (Gretlopt.OPT_W | Gretlopt.OPT_B | Gretlopt.OPT_S | Gretlopt.OPT_L)) {
                ok = 0;
            }
        }
    }

    if (ci == GretlCmdIndex.modtest && (opt & Gretlopt.OPT_A) &&
        gretl_is_regular_panel_model(pmod)) {
        return 1;
    }

    if (ok && !dataset_is_time_series(dset)) {
        /* rule out time-series-only tests */
        if (ci == GretlCmdIndex.cusum || ci == GretlCmdIndex.qlrtest || ci == GretlCmdIndex.bds ||
            (ci == GretlCmdIndex.modtest && (opt & (Gretlopt.OPT_H | Gretlopt.OPT_A)))) {
            ok = 0;
        }
    }

    if (ok && !dataset_is_panel(dset)) {
        /* rule out panel-only tests */
        if (ci == GretlCmdIndex.panspec) {
            ok = 0;
        } else if (ci == GretlCmdIndex.modtest && (opt & (Gretlopt.OPT_P | Gretlopt.OPT_D))) {
            ok = 0;
        }
    }

    if (ok && pmod.ncoeff - pmod.ifc <= 1 && ci == GretlCmdIndex.vif) {
        /* needs at least two independent variables */
        ok = 0;
    }

    if (ok && ci == GretlCmdIndex.modtest && (opt & Gretlopt.OPT_C)) {
        /* common factor test */
        if (pmod.opt & Gretlopt.OPT_P) {
            /* ?? check what this means! */
            ok = 0;
        }
    }

    if (ok && ci == GretlCmdIndex.restrict) {
        if (gretl_model_get_int(pmod, "null-model")) {
            ok = 0;
        }
    }

    return ok;
}

private __gshared int gretl_model_count = 0;

int get_model_count()
{
    return gretl_model_count;
}

void set_model_count(int c)
{
    gretl_model_count = c;
}

int model_count_plus()
{
    return ++gretl_model_count;
}

void model_count_minus(Model * pmod)
{
    if (pmod == null) {
        --gretl_model_count;
    } else if (pmod.ID > 0) {
        gretl_model_count = pmod.ID - 1;
        pmod.ID = 0;
    }
}

void set_model_id(Model * pmod, Gretlopt opt)
{
    /* always record model estimation time */
    pmod.esttime = gretl_monotonic_time();

    /* Experimental: limit setting of sequential model
       number to "visible" models estimated in main
       script or session (not in functions).
    */
    if (opt & Gretlopt.OPT_A) {
        /* An auxiliary model? Likely, but OPT_A has
           special meaning for a few estimators
        */
        if (pmod.ci != GretlCmdIndex.dpanel && pmod.ci != GretlCmdIndex.garch &&
            pmod.ci != GretlCmdIndex.arma) {
            return;
        }
    }

    if ((opt & Gretlopt.OPT_Q) && !(opt & Gretlopt.OPT_W)) {
        /* --quiet and not --window, so "invisible" */
        return;
    } else if (gretl_function_depth() > 0) {
        /* model was estimated inside a function */
        return;
    }

    if (pmod.errcode == 0) {
        pmod.ID = ++gretl_model_count;
    }
}

void model_list_to_string(int * list, char * buf)
{
    char[8] numstr;
    int i;

    for (i=1; i<=list[0]; i++) {
        if (list[i] == LISTSEP) {
            strcat(buf, "; ");
        } else {
            sprintf(numstr.ptr, "%d ", list[i]);
            strcat(buf, numstr.ptr);
        }
    }
}

int highest_numbered_var_in_model(const(Model) * pmod,
                                  const(DatasetRaw) * dset)
{
    const(char) * vname;
    int i, v, vmax = 0;
    int gotsep = 0;

    if (pmod.ci == GretlCmdIndex.mle || pmod.ci == GretlCmdIndex.gmm || pmod.list == null) {
        return 0;
    }

    for (i=1; i<=pmod.list[0]; i++) {
        v = pmod.list[i];
        if (v == LISTSEP) {
            gotsep = 1;
            continue;
        }
        if (v >= dset.v) {
            /* temporary variables, already gone? */
            continue;
        }
        if ((pmod.ci == GretlCmdIndex.arma || pmod.ci == GretlCmdIndex.garch) && !gotsep) {
            /* real vars start after LISTSEP */
            continue;
        }
        if (v > vmax) {
            vmax = v;
        }
        if (pmod.ci == GretlCmdIndex.nls) {
            /* only the dependent var can be tested */
            break;
        }
    }

    /* clustered standard errors? */
    vname = gretl_model_get_cluster_vname(pmod);
    if (vname != null) {
        v = current_series_index(dset, vname);
        if (v > vmax) vmax = v;
    }
    vname = gretl_model_get_cluster_vname2(pmod);
    if (vname != null) {
        v = current_series_index(dset, vname);
        if (v > vmax) vmax = v;
    }

    /* auxiliary variables for some model types */

    if (pmod.ci == GretlCmdIndex.wls) {
        if (pmod.nwt > vmax) vmax = pmod.nwt;
    } else if (pmod.ci == GretlCmdIndex.intreg) {
        v = gretl_model_get_int(pmod, "lovar");
        if (v > vmax) vmax = v;
        v = gretl_model_get_int(pmod, "hivar");
        if (v > vmax) vmax = v;
    } else if (COUNT_MODEL(cast(GretlCmdIndex) pmod.ci)) {
        v = gretl_model_get_int(pmod, "offset_var");
        if (v > vmax) vmax = v;
    } else if (pmod.ci == GretlCmdIndex.duration) {
        v = gretl_model_get_int(pmod, "cens_var");
        if (v > vmax) vmax = v;
    }

    return vmax;
}

int mle_criteria(Model * pmod, int addk)
{
    int err = 0;

    if (na(pmod.lnL)) {
        pmod.criterion[InfoCriterion.aic] = NADBL;
        pmod.criterion[InfoCriterion.bic] = NADBL;
        pmod.criterion[InfoCriterion.hqc] = NADBL;
        err = 1;
    } else {
        int k = pmod.ncoeff + addk;
        int n = pmod.nobs;

        pmod.criterion[InfoCriterion.aic] = -2.0 * pmod.lnL + 2.0 * k;
        pmod.criterion[InfoCriterion.bic] = -2.0 * pmod.lnL + k * log(cast(double) n);
        pmod.criterion[InfoCriterion.hqc] = -2.0 * pmod.lnL + 2 * k * log(log(cast(double) n));
    }

    return err;
}

int model_use_zscore(const(Model) * pmod)
{
    if (pmod == null) {
        /* modprint */
        return 1;
    } else if (gretl_model_get_int(pmod, "dfcorr")) {
        /* override ASYMPTOTIC_MODEL if need be */
        return 0;
    } else if (pmod.ci == GretlCmdIndex.ols && (pmod.opt & Gretlopt.OPT_N)) {
        /* OLS with --no-df-corr */
        return 1;
    } else if (ASYMPTOTIC_MODEL(cast(GretlCmdIndex) pmod.ci)) {
        return 1;
    } else if (pmod.ci == GretlCmdIndex.panel && (pmod.opt & (Gretlopt.OPT_U | Gretlopt.OPT_N))) {
        return 1;
    } else if ((pmod.opt & Gretlopt.OPT_R) && libset_get_bool(SetKey.robust_z)) {
        return 1;
    } else if (gretl_model_get_int(pmod, "asy")) {
        return 1;
    } else {
        return 0;
    }
}

double coeff_pval(int ci, double x, int df)
{
    double p = NADBL;

    if (!na(x)) {
        if (df == 0 || ASYMPTOTIC_MODEL(cast(GretlCmdIndex) ci)) {
            p = normal_pvalue_2(x);
        } else {
            p = student_pvalue_2(df, x);
        }
    }

    return p;
}

double model_coeff_pval(const(Model) * pmod, double x)
{
    double p = NADBL;

    if (!na(x)) {
        if (model_use_zscore(pmod)) {
            p = normal_pvalue_2(x);
        } else if (pmod.ci == GretlCmdIndex.ar) {
            /* backward compat (?) */
            int k = pmod.arinfo.arlist[0];
            int dfd = pmod.dfd;

            if (k > 1) {
                dfd += pmod.ncoeff - k;
            }
            p = student_pvalue_2(dfd, x);
        } else {
            p = student_pvalue_2(pmod.dfd, x);
        }
    }

    return p;
}

/**
 * gretl_model_allocate_param_names:
 * @pmod: pointer to target model.
 * @k: number of strings to allocate.
 *
 * Allocate an array of @k strings to hold the names given to
 * the associated  coefficients, in a model where these strings
 * are not simply given by the names of the independent variables.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_model_allocate_param_names(Model * pmod, int k)
{
    pmod.params = strings_array_new_with_length(k, PNAMELEN);

    if (pmod.params == null) {
        pmod.errcode = GretlError.alloc;
    }

    if (!pmod.errcode) {
        pmod.nparams = k;
    }

    return pmod.errcode;
}

/**
 * gretl_model_set_param_name:
 * @pmod: pointer to target model.
 * @i: 0-based index of value to set.
 * @name: string value to set.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_model_set_param_name(Model * pmod, int i, const(char) * name)
{
    if (pmod.params == null || i < 0 || i >= pmod.nparams ||
        name == null) {
        return GretlError.data;
    } else {
        pmod.params[i][0] = '\0';
        if (strlen(name) >= PNAMELEN) {
            strncat(pmod.params[i], name, PNAMELEN - 2);
            strcat(pmod.params[i], "~");
        } else {
            strncat(pmod.params[i], name, PNAMELEN - 1);
        }
        return 0;
    }
}

private int count_coeffs(int p, const(char) * pmask,
                         int q, const(char) * qmask,
                         int P, int Q, int r, int ifc)
{
    int i, n = P + Q + r + ifc;

    for (i=0; i<p; i++) {
        if (arma_included(pmask, i)) {
            n++;
        }
    }

    for (i=0; i<q; i++) {
        if (arma_included(qmask, i)) {
            n++;
        }
    }

    return n;
}

/**
 * gretl_model_add_arma_varnames:
 * @pmod: pointer to target model.
 * @dset: dataset information.
 * @yno: ID number of dependent variable.
 * @p: non-seasonal (max) AR order.
 * @q: non-seasonal (max) MA order.
 * @pmask: for masking out specific AR lags.
 * @qmask: for masking out specific MA lags.
 * @P: seasonal AR order.
 * @Q: seasonal MA order.
 * @r: number of exogenous regressors (excluding the constant).
 *
 * Composes a set of names to be given to the regressors in an
 * ARMA model.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_model_add_arma_varnames(Model * pmod, const(DatasetRaw) * dset,
                                  int yno, int p, int q,
                                  const(char) * pmask, const(char) * qmask,
                                  int P, int Q,
                                  int r)
{
    int nullmod = 0;
    int nc, xstart;
    int i, j;

    nc = count_coeffs(p, pmask, q, qmask, P, Q, r, pmod.ifc);

    if (pmod.depvar != null) {
        free(pmod.depvar);
    }

    pmod.depvar = gretl_strdup(dset.varname[yno]);
    if (pmod.depvar == null) {
        pmod.errcode = GretlError.alloc;
        return 1;
    }

    if (pmod.nparams > 0 && pmod.params != null) {
        for (i=0; i<pmod.nparams; i++) {
            free(pmod.params[i]);
        }
        free(pmod.params);
    }

    nullmod = gretl_model_get_int(pmod, "null-model");

    if (nc == 0 && nullmod) {
        /* special case of null model */
        nc = 1;
    }

    pmod.params = strings_array_new_with_length(nc, VNAMELEN);
    if (pmod.params == null) {
        free(pmod.depvar);
        pmod.depvar = null;
        pmod.errcode = GretlError.alloc;
        return 1;
    }

    pmod.nparams = nc;

    if (pmod.ifc || nullmod) {
        strcpy(pmod.params[0], dset.varname[0]);
        j = 1;
    } else {
        j = 0;
    }

    for (i=0; i<p; i++) {
        if (arma_included(pmask, i)) {
            sprintf(pmod.params[j++], "phi_%d", i + 1);
        }
    }

    for (i=0; i<P; i++) {
        sprintf(pmod.params[j++], "Phi_%d", i + 1);
    }

    for (i=0; i<q; i++) {
        if (arma_included(qmask, i)) {
            sprintf(pmod.params[j++], "theta_%d", i + 1);
        }
    }

    for (i=0; i<Q; i++) {
        sprintf(pmod.params[j++], "Theta_%d", i + 1);
    }

    xstart = arma_depvar_pos(pmod) + 1;

    for (i=0; i<r; i++) {
        strcpy(pmod.params[j++], dset.varname[pmod.list[xstart+i]]);
    }

    return 0;
}

/**
 * gretl_model_add_panel_varnames:
 * @pmod: pointer to target model.
 * @dset: dataset information.
 * @ulist: list of index numbers of cross-sectional
 * units included in the model.
 *
 * Composes a set of names to be given to the regressors in an
 * panel model.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_model_add_panel_varnames(Model * pmod, const(DatasetRaw) * dset,
                                   const(int) * ulist)
{
    int np = pmod.ncoeff;
    int i, j, v;

    pmod.depvar = gretl_strdup(dset.varname[pmod.list[1]]);
    if (pmod.depvar == null) {
        pmod.errcode = GretlError.alloc;
        return 1;
    }

    pmod.params = strings_array_new_with_length(np, VNAMELEN);
    if (pmod.params == null) {
        pmod.errcode = GretlError.alloc;
        return 1;
    }

    pmod.nparams = np;

    j = 1;
    for (i=0; i<np; i++) {
        v = pmod.list[i+2];
        if (v < dset.v) {
            strcpy(pmod.params[i], dset.varname[v]);
        } else if (ulist != null) {
            sprintf(pmod.params[i], "ahat_%d", ulist[j++]);
        } else {
            sprintf(pmod.params[i], "ahat_%d", j++);
        }
    }

    return 0;
}

/**
 * gretl_model_add_allocated_varnames:
 * @pmod: pointer to target model.
 * @vnames: array of names of independent variables.
 *
 * Attaches an allocated set of variable names to be used
 * when printing model results, for use in special cases
 * where we can't just reference names from the list of
 * regressors attached to the model.  The number of strings
 * must match the number of coefficients, given by the
 * %ncoeff member of @pmod.
 *
 * Note that @pmod "takes charge" of the array @vnames:
 * this will be freed when the model is destroyed.
 */

void gretl_model_add_allocated_varnames(Model * pmod, char ** vnames)
{
    pmod.nparams = pmod.ncoeff;
    pmod.params = vnames;
}

/**
 * gretl_model_add_y_median:
 * @pmod: pointer to target model.
 * @y: array containing the dependent variable.
 *
 * Calculates the median of @y using the valid observations
 * with the model's sample range and attaches the median
 * to the model as data under the key %ymedian.
 *
 * Returns: 0 on success or error code on error.
 */

int gretl_model_add_y_median(Model * pmod, const(double) * y)
{
    int T = pmod.t2 - pmod.t1 + 1;
    double * sy;
    double m;
    int t, n, ok, n2p;

    sy = cast(double *) malloc(T * double.sizeof);

    if (sy == null) {
        return GretlError.alloc;
    }

    n = 0;
    for (t=pmod.t1; t<=pmod.t2; t++) {
        if (pmod.uhat != null) {
            ok = !na(pmod.uhat[t]);
        } else {
            ok = !model_missing(pmod, t);
        }
        if (ok) {
            sy[n++] = y[t];
        }
    }

    if (n == 0) {
        free(sy);
        return GretlError.data;
    }

    qsort(sy, n, double.sizeof, &gretl_compare_doubles);

    n2p = (T = n / 2) + 1;
    m = (n % 2)? sy[n2p - 1] : 0.5 * (sy[T - 1] + sy[n2p - 1]);

    gretl_model_set_double(pmod, "ymedian", m);

    free(sy);

    return 0;
}

int gretl_model_add_normality_test(Model * pmod, double X2)
{
    ModelTest * test = model_test_new(ModelTestType.normal);
    int err = 0;

    if (test != null) {
        model_test_set_teststat(test, cast(ubyte) TestStats.normal_chisq);
        model_test_set_dfn(test, 2);
        model_test_set_value(test, X2);
        model_test_set_pvalue(test, chisq_cdf_comp(2, X2));
        maybe_add_test_to_model(pmod, test);
    } else {
        err = GretlError.alloc;
    }

    return err;
}

ModelTest * gretl_model_get_test(Model * pmod, ModelTestType ttype)
{
    int i;

    for (i=0; i<pmod.ntests; i++) {
        if (pmod.tests[i].type == ttype) {
            return &pmod.tests[i];
        }
    }

    return null;
}

private int xneq(double x, double y)
{
    double reldiff;

    if (x == 0.0) {
        reldiff = fabs(y);
    } else if (y == 0.0) {
        reldiff = fabs(x);
    } else if (x > y) {
        reldiff = fabs((x - y) / y);
    } else {
        reldiff = fabs((y - x) / x);
    }

    return reldiff > 1.5e-12;
}

/* try to tell if an OLS model with two independent variables is
   actually a quadratic model (x_2 = x_1^2). */

private int model_is_quadratic(const(Model) * pmod, const(int) * xlist,
                               const(DatasetRaw) * dset)
{
    const(double) * x1 = dset.Z[xlist[2]];
    const(double) * x2 = dset.Z[xlist[3]];
    int t;

    for (t=pmod.t1; t<=pmod.t2; t++) {
        if (!na(x1[t]) && xneq(x2[t], x1[t] * x1[t])) {
            return 0;
        }
    }

    return 1;
}

/* we'll be a bit conservative here so as not to make
   fools of ourselves */

private bool fitted_formula_ok(int c)
{
    return (c == GretlCmdIndex.ols || c == GretlCmdIndex.wls ||
            c == GretlCmdIndex.hsk || c == GretlCmdIndex.ivreg ||
            c == GretlCmdIndex.lad || c == GretlCmdIndex.logistic);
}

/**
 * gretl_model_get_fitted_formula:
 * @pmod: pointer to target model.
 * @xvar: ID number of variable that _may_ be "x" in the model.
 * @dset: dataset information.
 *
 * If @pmod is a simple linear, quadratic or logistic model,
 * and if @xvar is in fact the "x" variable from the model,
 * returns a string representing the formula for generating the
 * fitted values as a function of x.  This formula may be used
 * in the context of a fitted versus actual plot.
 *
 * Returns: formula for fitted values, or %NULL if this is
 * not available.
 */

char * gretl_model_get_fitted_formula(const(Model) * pmod, int xvar,
                                      const(DatasetRaw) * dset)
{
    const(DatasetRaw) * mset;
    int * xlist = null;
    char * ret = null;

    if (xvar == 0 || pmod.ncoeff > 3 || !fitted_formula_ok(pmod.ci)) {
        return null;
    }

    xlist = gretl_model_get_x_list(pmod);
    if (xlist == null) {
        return null;
    }

    if (pmod.dataset != null) {
        mset = pmod.dataset;
    } else {
        mset = dset;
    }

    gretl_push_c_numeric_locale();

    if (pmod.ci == GretlCmdIndex.logistic) {
        if (pmod.ifc && pmod.ncoeff == 2 && xvar == pmod.list[3]) {
            double lmax = gretl_model_get_double(pmod, "lmax");

            if (!na(lmax)) {
                ret = gretl_strdup_printf("yformula: %g/(1.0+exp(-(%g%s%g*x)))",
                                          lmax, pmod.coeff[0],
                                          (pmod.coeff[1] >= 0)? "+".ptr : "".ptr,
                                          pmod.coeff[1]);
            }
        }
    } else if (!pmod.ifc && pmod.ncoeff == 1 && xvar == xlist[1]) {
        ret = gretl_strdup_printf("yformula: %g*x", pmod.coeff[0]);
    } else if (pmod.ifc && pmod.ncoeff == 2 && xvar == xlist[2]) {
        ret = gretl_strdup_printf("yformula: %g%s%g*x", pmod.coeff[0],
                                  (pmod.coeff[1] >= 0)? "+".ptr : "".ptr,
                                  pmod.coeff[1]);
    } else if (pmod.ifc && pmod.ncoeff == 3 && xvar == xlist[2]) {
        if (model_is_quadratic(pmod, xlist, mset)) {
            ret = gretl_strdup_printf("yformula: %g%s%g*x%s%g*x**2", pmod.coeff[0],
                                      (pmod.coeff[1] >= 0)? "+".ptr : "".ptr,
                                      pmod.coeff[1],
                                      (pmod.coeff[2] >= 0)? "+".ptr : "".ptr,
                                      pmod.coeff[2]);
        }
    }

    gretl_pop_c_numeric_locale();
    free(xlist);

    return ret;
}

/**
 * gretl_model_set_name:
 * @pmod: pointer to target model.
 * @name: the name to give the model.
 *
 * Sets the name of the given model; this is used in
 * printing the model and in displaying it in the
 * gretl GUI. Note that a model's name must be no more
 * than #MAXSAVENAME bytes in length, including the
 * terminating NUL byte; @name is truncated if it is
 * too long.
 */

void gretl_model_set_name(Model * pmod, const(char) * name)
{
    if (name == pmod.name) {
        return;
    }

    if (pmod.name == null) {
        pmod.name = cast(char *) calloc(1, MAXSAVENAME);
    }

    if (pmod.name != null) {
        *pmod.name = '\0';
        if (strlen(name) >= MAXSAVENAME) {
            char * tmp = g_strdup(name);

            strcpy(pmod.name, gretl_utf8_truncate(tmp, MAXSAVENAME-1));
            g_free(tmp);
        } else {
            strcpy(pmod.name, name);
        }
    }
}

/**
 * gretl_model_get_name:
 * @pmod: pointer to gretl model.
 *
 * Returns: the name that has been set for @pmod, if any.
 * Note that the value returned may be %NULL.
 */

const(char) * gretl_model_get_name(const(Model) * pmod)
{
    if (pmod != null) {
        return pmod.name;
    }

    return null;
}

/**
 * gretl_model_get_scalar:
 * @pmod: pointer to target model.
 * @idx: index for the scalar value that is wanted.
 * @dset: dataset struct.
 * @err: location to receive error code (required).
 *
 * Retrieves a specified scalar statistic from @pmod:
 * @idx must be less than %M_SCALAR_MAX.
 *
 * Returns: the requested statistic, or #NADBL on failure,
 * in which case @err will contain a non-zero error code.
 */

double gretl_model_get_scalar(Model * pmod, ModelDataIndex idx,
                              DatasetRaw * dset, int * err)
{
    double x = NADBL;

    if (pmod == null) {
        *err = GretlError.badstat;
        return x;
    }

    if (idx == ModelDataIndex.gmmcrit && pmod.ci != GretlCmdIndex.gmm) {
        *err = GretlError.badstat;
        return x;
    }

    switch (idx) {
    case ModelDataIndex.ess:
    case ModelDataIndex.gmmcrit:
        x = pmod.ess;
        break;
    case ModelDataIndex.rsq:
        x = pmod.rsq;
        break;
    case ModelDataIndex.lnl:
        x = pmod.lnL;
        break;
    case ModelDataIndex.aic:
        x = pmod.criterion[InfoCriterion.aic];
        break;
    case ModelDataIndex.bic:
        x = pmod.criterion[InfoCriterion.bic];
        break;
    case ModelDataIndex.hqc:
        x = pmod.criterion[InfoCriterion.hqc];
        break;
    case ModelDataIndex.sigma:
        x = pmod.sigma;
        break;
    case ModelDataIndex.trsq:
        if (!na(pmod.rsq)) {
            x = pmod.nobs * pmod.rsq;
        }
        break;
    case ModelDataIndex.df:
        x = cast(double) pmod.dfd;
        break;
    case ModelDataIndex.ncoeff:
        x = cast(double) pmod.ncoeff; /* is this always available? */
        break;
    case ModelDataIndex.t:
        x = cast(double) pmod.nobs;
        break;
    case ModelDataIndex.dw:
        x = pmod.dw;
        break;
    case ModelDataIndex.dwpval:
        x = get_DW_pvalue_for_model(pmod, dset, err);
        break;
    case ModelDataIndex.fstt:
        x = pmod.fstt;
        break;
    case ModelDataIndex.chisq:
        x = pmod.chisq;
        break;
    default:
        break;
    }

    if (idx != ModelDataIndex.dwpval && na(x)) {
        *err = GretlError.badstat;
    }

    return x;
}

/**
 * gretl_model_get_series:
 * @x: series to fill (must be of length dset->n).
 * @pmod: pointer to target model.
 * @dset: dataset information.
 * @idx: index for the series that is wanted.
 *
 * Retrieves from @pmod a copy of a specified series (for
 * example, regression residuals); @idx must be greater than
 * %M_SCALAR_MAX and less than %M_SERIES_MAX.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_model_get_series(double * x, Model * pmod,
                           const(DatasetRaw) * dset,
                           ModelDataIndex idx)
{
    const(double) * src = null;
    int t;

    if (pmod.t2 - pmod.t1 + 1 > dset.n ||
        model_sample_problem(pmod, dset)) {
        gretl_errmsg_set(
               ((idx == ModelDataIndex.uhat)?
                _("Can't retrieve uhat: data set has changed") :
                (idx == ModelDataIndex.yhat)?
                _("Can't retrieve yhat: data set has changed") :
                (idx == ModelDataIndex.h)?
                _("Can't retrieve ht: data set has changed") :
                _("Can't retrieve series: data set has changed")).ptr);
        return GretlError.badstat;
    }

    if (pmod.ci == GretlCmdIndex.biprobit &&
        (idx == ModelDataIndex.uhat || idx == ModelDataIndex.yhat)) {
        return GretlError.badstat;
    }

    if (idx == ModelDataIndex.uhat) {
        src = pmod.uhat;
    } else if (idx == ModelDataIndex.yhat) {
        src = pmod.yhat;
    } else if (idx == ModelDataIndex.llt) {
        src = cast(const(double) *) gretl_model_get_data(pmod, "llt");
    } else if (idx == ModelDataIndex.ahat) {
        src = cast(const(double) *) gretl_model_get_data(pmod, "ahat");
    } else if (idx == ModelDataIndex.h) {
        src = cast(const(double) *) gretl_model_get_data(pmod, "garch_h");
    }

    if (src == null && idx != ModelDataIndex.sample) {
        return GretlError.badstat;
    }

    /* allow for internal "just testing" usage */
    if (x == null) {
        return 0;
    }

    if (idx == ModelDataIndex.sample) {
        for (t=0; t<dset.n; t++) {
            if (t < pmod.t1 || t > pmod.t2) {
                x[t] = 0.0;
            } else {
                x[t] = model_missing(pmod, t)? 0 : 1;
            }
        }
    } else {
        for (t=0; t<dset.n; t++) {
            if (t < pmod.t1 || t > pmod.t2) {
                x[t] = NADBL;
            } else {
                x[t] = src[t];
            }
        }
    }

    return 0;
}

private GretlMatrixRaw * tauvec_result(const(Model) * pmod,
                                    const(GretlVector) * tau,
                                    int idx)
{
    GretlMatrixRaw * S = cast(GretlMatrixRaw *) gretl_model_get_data(pmod, "rq_sequence");
    GretlMatrixRaw * ret = null;

    if (S != null) {
        double xij;
        int src = idx == ModelDataIndex.coeff ? 0 : 1;
        int c = gretl_vector_get_length(tau);
        int r = pmod.ncoeff;
        int i, j, t;

        ret = gretl_matrix_alloc(r, c);

        t = 0;
        for (i=0; i<r; i++) {
            t = i * c;
            for (j=0; j<c; j++) {
                xij = gretl_matrix_get(S, t++, src);
                gretl_matrix_set(ret, i, j, xij);
            }
        }
    }

    return ret;
}

private GretlMatrixRaw * model_get_estvec(const(Model) * pmod, int idx, int * err)
{
    const(double) * src;
    GretlMatrixRaw * v = null;
    int i;

    if ((v = cast(GretlMatrixRaw *) gretl_model_get_data(pmod, "rq_tauvec")) != null) {
        GretlMatrixRaw * ret = tauvec_result(pmod, v, idx);

        if (ret == null) {
            *err = GretlError.badstat;
        }
        return ret;
    }

    src = (idx == ModelDataIndex.coeff)? pmod.coeff : pmod.sderr;

    if (src == null) {
        *err = GretlError.badstat;
        return null;
    }

    v = gretl_column_vector_alloc(pmod.ncoeff);
    if (v == null) {
        *err = GretlError.alloc;
    } else {
        for (i=0; i<pmod.ncoeff; i++) {
            gretl_vector_set(v, i, src[i]);
        }
    }

    return v;
}

private GretlMatrixRaw * model_get_hatvec(const(Model) * pmod, int idx, int * err)
{
    GretlMatrixRaw * v = null;
    const(double) * src = null;
    int t;

    if (idx == ModelDataIndex.uhat) {
        src = pmod.uhat;
    } else if (idx == ModelDataIndex.yhat) {
        src = pmod.yhat;
    } else if (idx == ModelDataIndex.llt) {
        src = cast(const(double) *) gretl_model_get_data(pmod, "llt");
    }

    if (src == null) {
        *err = GretlError.badstat;
        return null;
    }

    /* FIXME: we reject NAs when creating a gretl_matrix -- but should
       we just substitute NaNs?
    */
    for (t=pmod.t1; t<=pmod.t2; t++) {
        if (na(src[t])) {
            *err = GretlError.missdata;
            break;
        }
    }

    if (!*err) {
        v = gretl_column_vector_alloc(pmod.t2 - pmod.t1 + 1);
        if (v == null) {
            *err = GretlError.alloc;
        } else {
            for (t=pmod.t1; t<=pmod.t2; t++) {
                gretl_vector_set(v, t - pmod.t1, src[t]);
            }
        }
    }

    return v;
}

private GretlMatrixRaw * model_get_special_vec(const(Model) * pmod, ModelDataIndex idx, int * err)
{
    GretlMatrixRaw * v = null;
    double * mdata = null;
    int t;

    if (idx == ModelDataIndex.ahat) {
        mdata = cast(double *) gretl_model_get_data(pmod, "ahat");
    } else if (idx == ModelDataIndex.h) {
        mdata = cast(double *) gretl_model_get_data(pmod, "garch_h");
    }

    if (mdata == null) {
        fprintf(stderr, "model_get_special_vec: mdata is NULL\n");
        *err = GretlError.badstat;
    }

    v = gretl_column_vector_alloc(pmod.t2 - pmod.t1 + 1);

    if (v == null) {
        *err = GretlError.alloc;
    } else {
        gretl_matrix_set_t1(v, pmod.t1);
        gretl_matrix_set_t2(v, pmod.t2);
        for (t=pmod.t1; t<=pmod.t2; t++) {
            /* FIXME: is this indexation always right? */
            gretl_vector_set(v, t - pmod.t1, mdata[t]);
        }
    }

    return v;
}

private GretlMatrixRaw * model_get_rhovec(const(Model) * pmod, int * err)
{
    GretlMatrixRaw * r = null;

    if (pmod.ci == GretlCmdIndex.ar1) {
        double x = gretl_model_get_double(pmod, "rho_gls");

        r = gretl_matrix_from_scalar(x);
    } else if (pmod.ci != GretlCmdIndex.ar) {
        r = gretl_matrix_from_scalar(pmod.rho);
    } else if (pmod.arinfo == null ||
               pmod.arinfo.arlist == null ||
               pmod.arinfo.rho == null) {
        *err = GretlError.invarg;
    } else {
        int i, l0 = pmod.arinfo.arlist[0];
        int lmax = pmod.arinfo.arlist[l0];

        r = gretl_vector_alloc(lmax);
        if (r != null) {
            gretl_matrix_zero(r);
            for (i=1; i<=l0; i++) {
                gretl_vector_set(r, pmod.arinfo.arlist[i] - 1,
                                 pmod.arinfo.rho[i-1]);
            }
        }
    }

    return r;
}

private GretlMatrixRaw * model_get_intervals_matrix(const(Model) * pmod, int * err)
{
    GretlMatrixRaw * m;
    GretlMatrixRaw * ci;

    ci = cast(GretlMatrixRaw *) gretl_model_get_data(pmod, "coeff_intervals");
    if (ci == null) {
        *err = GretlError.badstat;
        return null;
    }

    m = gretl_matrix_copy(ci);

    if (m == null) {
        *err = GretlError.alloc;
    }

    return m;
}

private GretlMatrixRaw * model_get_special_test(const(Model) * pmod,
                                             int type, int * err)
{
    GretlMatrixRaw * r = null;
    int i, found = 0;

    for (i=0; i<pmod.ntests; i++) {
        found = (pmod.tests[i].type == type);
        if (found) {
            r = gretl_vector_alloc(3);
            if (r == null) {
                *err = GretlError.alloc;
                return null;
            }
            r.val[0] = pmod.tests[i].value;
            r.val[1] = pmod.tests[i].dfn;
            r.val[2] = pmod.tests[i].pvalue;
            break;
        }
    }

    if (!found) {
        *err = GretlError.badstat;
    }

    return r;
}

/**
 * gretl_model_get_matrix:
 * @pmod: pointer to target model.
 * @idx: index for the matrix that is wanted.
 * @err: location to receive error code (required).
 *
 * Retrieves from @pmod a copy of a specified matrix (for
 * example, regression residuals); @idx must be greater than
 * %M_SERIES_MAX and less than %M_MATRIX_MAX.
 *
 * Returns: the allocated matrix, or %NULL on failure,
 * in which case @err will contain a non-zero error code.
 */

GretlMatrixRaw * gretl_model_get_matrix(Model * pmod, ModelDataIndex idx,
                                     int * err)
{
    GretlMatrixRaw * M = null;

    if (pmod == null) {
        *err = GretlError.badstat;
        return null;
    }

    if ((idx == ModelDataIndex.uhat || idx == ModelDataIndex.yhat) &&
        (pmod.ci == GretlCmdIndex.biprobit || gretl_is_between_model(pmod))) {
        /* special: uhat and yhat are matrices */
        const(char) * utag = pmod.ci == GretlCmdIndex.biprobit ? "bp_uhat" : "uhat";
        const(char) * ytag = pmod.ci == GretlCmdIndex.biprobit ? "bp_yhat" : "yhat";
        GretlMatrixRaw * P;

        if (idx == ModelDataIndex.uhat) {
            P = cast(GretlMatrixRaw *) gretl_model_get_data(pmod, utag);
        } else {
            P = cast(GretlMatrixRaw *) gretl_model_get_data(pmod, ytag);
        }

        if (P == null) {
            *err = GretlError.badstat;
        } else {
            M = gretl_matrix_copy(P);
            if (M == null) {
                *err = GretlError.alloc;
            }
        }
        return M;
    }

    switch (idx) {
    case ModelDataIndex.uhat:
    case ModelDataIndex.yhat:
    case ModelDataIndex.llt:
        M = model_get_hatvec(pmod, idx, err);
        break;
    case ModelDataIndex.coeff:
    case ModelDataIndex.se:
        M = model_get_estvec(pmod, idx, err);
        break;
    case ModelDataIndex.coeff_ci:
        M = model_get_intervals_matrix(pmod, err);
        break;
    case ModelDataIndex.vcv:
        M = gretl_vcv_matrix_from_model(pmod, null, err);
        break;
    case ModelDataIndex.ahat:
        if (gretl_model_get_data(pmod, "ahat") == null) {
            *err = GretlError.badstat;
        } else {
            M = model_get_special_vec(pmod, ModelDataIndex.ahat, err);
        }
        break;
    case ModelDataIndex.h:
    case ModelDataIndex.sigma:
        if (gretl_model_get_data(pmod, "garch_h") == null) {
            *err = GretlError.badstat;
        } else {
            M = model_get_special_vec(pmod, ModelDataIndex.h, err);
        }
        break;
    case ModelDataIndex.rho:
        M = model_get_rhovec(pmod, err);
        break;
    case ModelDataIndex.hausman:
        if (pmod.ci == GretlCmdIndex.ivreg) {
            M = model_get_special_test(pmod, ModelTestType.iv_hausman, err);
        } else if (pmod.ci == GretlCmdIndex.panel) {
            M = model_get_special_test(pmod, ModelTestType.panel_hausman, err);
        } else {
            *err = GretlError.badstat;
        }
        break;
    case ModelDataIndex.sargan:
        M = model_get_special_test(pmod, ModelTestType.sargan, err);
        break;
    case ModelDataIndex.ehat:
        M = cast(GretlMatrixRaw *) gretl_model_get_data(pmod, "ehat");
        if (M == null) {
            *err = GretlError.badstat;
        } else {
            M = gretl_matrix_copy(M);
            if (M == null) {
                *err = GretlError.alloc;
            }
        }
        break;
    case ModelDataIndex.oddsratios:
        M = cast(GretlMatrixRaw *) gretl_model_get_data(pmod, "oddsratios");
        if (M == null) {
            *err = GretlError.badstat;
        } else {
            M = gretl_matrix_copy(M);
            if (M == null) {
                *err = GretlError.alloc;
            }
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

private double get_vcv_element(Model * pmod, const(char) * s,
                               const(DatasetRaw) * dset)
{
    char[VNAMELEN] v1str;
    char[VNAMELEN] v2str;
    char[16] fmt;
    int v1 = 0, v2 = 0;
    int i, j, k, gotit;
    double ret = NADBL;

    if (pmod == null || pmod.list == null) {
        return NADBL;
    }

    sprintf(fmt.ptr, "%%%d[^,],%%%ds", VNAMELEN-1, VNAMELEN-1);

    if (sscanf(s, fmt.ptr, v1str.ptr, v2str.ptr) != 2) {
        return NADBL;
    }

    v1 = gretl_model_get_param_number(pmod, dset, v1str.ptr);
    v2 = gretl_model_get_param_number(pmod, dset, v2str.ptr);

    if (v1 < 0 || v2 < 0) {
        return NADBL;
    }

    /* make model vcv matrix if need be */
    if (pmod.vcv == null && makevcv(pmod, pmod.sigma)) {
        return NADBL;
    }

    /* now find the right entry */
    if (v1 > v2) {
        k = v1;
        v1 = v2;
        v2 = k;
    }

    gotit = 0;
    k = 0;
    for (i=0; i<pmod.ncoeff && !gotit; i++) {
        for (j=0; j<pmod.ncoeff; j++) {
            if (j < i) {
                continue;
            }
            if (i == v1 && j == v2) {
                ret = pmod.vcv[k];
                gotit = 1;
                break;
            }
            k++;
        }
    }

    return ret;
}

/* retrieve a specific element from one of the arrays of data
   on a model */

double gretl_model_get_data_element(Model * pmod, int idx, const(char) * s,
                                    const(DatasetRaw) * dset, int * err)
{
    GretlObjType type;
    double x = NADBL;
    int vi = 0;

    if (pmod == null) {
        pmod = get_genr_model(&type);
        if (pmod == null || type != GretlObjType.eqn) {
            pmod = get_last_model(&type);
            if (pmod == null || type != GretlObjType.eqn) {
                *err = GretlError.invarg;
                return x;
            }
        }
    }

    if (idx == ModelDataIndex.rho) {
        int k = atoi(s);

        if (k == 1 && pmod.ci == GretlCmdIndex.ar1) {
            x = gretl_model_get_double(pmod, "rho_gls");
        } else if (k == 1 && pmod.ci != GretlCmdIndex.ar) {
            x = pmod.rho;
        } else if (pmod.arinfo == null ||
                   pmod.arinfo.arlist == null ||
                   pmod.arinfo.rho == null) {
            *err = GretlError.invarg;
        } else {
            vi = in_gretl_list(pmod.arinfo.arlist, k);
            if (vi > 0) {
                x = pmod.arinfo.rho[vi-1];
            } else {
                *err = GretlError.invarg;
            }
        }
    } else if (idx == ModelDataIndex.vcv) {
        x = get_vcv_element(pmod, s, dset);
        if (na(x)) {
            *err = GretlError.invarg;
        }
    } else if (idx == ModelDataIndex.coeff || idx == ModelDataIndex.se) {
        vi = gretl_model_get_param_number(pmod, dset, s);
        if (vi < 0) {
            *err = GretlError.invarg;
        } else {
            if (idx == ModelDataIndex.coeff && pmod.coeff != null) {
                x = pmod.coeff[vi];
            } else if (idx == ModelDataIndex.se && pmod.sderr != null) {
                x = pmod.sderr[vi];
            } else {
                *err = GretlError.invarg;
            }
        }
    }

    return x;
}

int exact_fit_check(const(Model) * pmod, Prn * prn)
{
    if (pmod.rsq == 1.0) {
        pputs(prn, _("The model exhibits an exact linear fit").ptr);
        pputc(prn, '\n');
        return 1;
    }

    return 0;
}

void maybe_suppress_time_dummies(Model * pmod, int ndum)
{
    const(char) * s = get_optval_string(cast(GretlCmdIndex) pmod.ci, Gretlopt.OPT_D);

    if (s != null && !strcmp(s, "noprint")) {
        gretl_model_set_int(pmod, "skipdums", ndum);
    }
}
