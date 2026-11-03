/* Port of lib/src/gretl_bundle.c: gretl bundles (string-keyed hash tables of scalars, strings, matrices, series,
   lists, arrays and other bundles): creation, member access, copying, merging, comparing, printing, XML and
   binary (de)serialisation, the "sysinfo" bundle, function-package helpers and the bundle-valued functions of
   the genr language. The type GretlBundle (struct gretl_bundle_ in the C, opaque in gretl_bundle.h) is declared in
   libgretlh.d, and the glib hash table and list are ghashtable.d and glist.d. */

module gretl.gretl_bundle;


import core.stdc.ctype : isdigit, isspace;
import core.stdc.float_ : DBL_EPSILON;
import core.stdc.math : round, nearbyint, atan2, cos, sin, ceil, exp, fabs, floor, fmod, isinf, isnan, log, log10, pow, sqrt;
import core.stdc.stdio : fread, fwrite, FILE, fclose, fopen, fgets, fputc, fputs, rewind, fprintf, printf, puts, sprintf, snprintf, sscanf, stderr;
import core.stdc.stdlib : abs, atof, atoi, calloc, free, getenv, malloc, qsort, realloc;
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




import core.stdc.config : c_long;

/* the C library's strtol takes a plain char ** end pointer */
import gretl.gretl_intl : gretl_pop_c_numeric_locale, gretl_push_c_numeric_locale;
import gretl.gstrfuncs : g_strdup;
import gretl.ghashtable;
import gretl.glist;
import gretl.kalmanh : Kalman;
import gretl.gstring;
import gretl.gretl_prnh : PrnType;
version (Posix) import core.sys.posix.unistd : gethostname; else import gretl.winport : gethostname;

/* g_get_host_name(): the name of the host, in a static buffer */
const(char) * g_get_host_name()
{
    __gshared char[256] hostname = 0;

    gethostname(hostname.ptr, hostname.length - 1);
    return hostname.ptr;
}
import gretl.gretl_paths : gretl_dotdir;
import gretl.gfileutils : g_file_get_contents;
import gretl.genmain : generate_matrix;
import gretl.genmainh : ModelDataIndex, RetrievalIndex;
import gretl.gerror : GError, g_error_free;
import gretl.glib : gboolean, gint, gpointer, gsize, guint;
import gretl.gmem : g_free;
import gretl.graphing : gnuplot_version_string;
import gretl.gretl_array : gretl_array_copy, gretl_array_deserialize, gretl_array_destroy, gretl_array_get_bundle, gretl_array_get_content_type, gretl_array_get_data, gretl_array_get_length, gretl_array_get_strings, gretl_array_get_type, gretl_array_new, gretl_array_serialize, gretl_array_set_string, gretl_arrays_are_equal;
import gretl.gretl_bundleh : BundleType, BundledItem;
import gretl.gretl_errors : gretl_errmsg_set, gretl_errmsg_sprintf;
import gretl.gretl_list : gretl_list_cmp, gretl_list_copy;
import gretl.gretl_model : bundlize_model_data_items, gretl_model_get_depvar_name, gretl_model_get_double, gretl_model_get_matrix, gretl_model_get_param_names, gretl_model_get_scalar, gretl_model_get_series, gretl_model_get_y_list;
import gretl.gretl_paths : foreign_info, gretl_build_path, gretl_maybe_prepend_dir, gretl_maybe_switch_dir, gretl_remove;
import gretl.gretl_prn : Prn, gretl_gzip_print_new, gretl_print_destroy, gretl_print_new, gretl_print_new_with_filename, pprintf, pputc, pputs;
import gretl.gretl_utils : avx_support, blas_variant_string, get_blas_details, get_cpu_details, get_stack_size, gretl_int_from_double, gretl_uint32_from_double, gretl_uint64_from_double;
import gretl.gretl_xml : gretl_list_serialize, gretl_matrix_serialize, gretl_xml_get_double_array, gretl_xml_get_list, gretl_xml_get_matrix, gretl_xml_get_prop_as_bool, gretl_xml_get_prop_as_int, gretl_xml_get_type_property, gretl_xml_header, gretl_xml_node_get_double, gretl_xml_node_get_int, gretl_xml_node_get_trimmed_string, gretl_xml_node_get_uint32, gretl_xml_node_get_uint64, gretl_xml_open_doc_root, gretl_xml_put_string;
import gretl.gstrfuncs : g_ascii_strcasecmp;
import gretl.libset : gretl_in_gui_mode;
import gretl.matrix_extra : gretl_matrix_from_scalar, gretl_matrix_write_constructor;
import gretl.objstackh : GretlObjType;
import gretl.plugins : get_plugin_function;
import gretl.printout : bufspace;
import gretl.strutils : gretl_namechar_spn, gretl_strdup, gretl_strndup, has_suffix, strings_array_donate, strings_array_free, strings_array_new, strings_array_position;
import gretl.system : equation_system_bundlize;
import gretl.systemh : EquationSystem;
import gretl.uservarh : UserVar;
import gretl.var : gretl_VAR_bundlize;
import gretl.uservar : destroy_private_uvars, get_user_var_by_data, get_user_var_of_type_by_name, private_bundle_add_as_shell, private_matrix_add_as_shell, private_scalar_add, user_var_add, user_var_get_name, user_var_get_value, user_var_replace_value, user_var_steal_value, user_var_unstack_value;
import gretl.gretl_typemap : gretl_is_scalar_type, gretl_type_get_name, gretl_type_get_order;
import gretl.gretl_commands : gretl_command_word;

/* unported dependencies: per-file declarations from the 2026c headers */
extern (C) {
    int genr_function_word(const(char) * s);
    void * get_genr_model(GretlObjType * type);
    void * get_last_model(GretlObjType * type);
    int gretl_max_mpi_processes();
    int gretl_n_physical_cores();
    int gretl_n_processors();
    int is_user_function(const(char) * name);
    GretlBundle * kalman_bundle_copy(const(GretlBundle) * src, int * err);
    int kalman_bundle_n_members(GretlBundle * b);
    GretlBundle * kalman_deserialize(void * p1, void * p2, int * err);
    void kalman_free(Kalman * K);
    Kalman * kalman_new_minimal(GretlMatrixRaw ** M, int * copy, int nmat, int dkvar, int * err);
    int kalman_serialize(void * kptr, Prn * prn);
    int maybe_delete_kalman_element(void * kptr, const(char) * key, int * err);
    void * maybe_retrieve_kalman_element(void * kptr, const(char) * key, GretlType * type, int * private_, int * err);
    int maybe_set_kalman_element(void * kptr, const(char) * key, void * vptr, GretlType vtype, int copy, int * err);
    int memory_stats(double * vals);
    const(char) * mvarname(int t);
    int print_kalman_bundle_info(void * kptr, Prn * prn);
    void set_bundle_pkg(const(char) * pkgname);
}

/**
 * gretl_bundle:
 *
 * An opaque type; use the relevant accessor functions.
 */

__gshared GretlBundle * sysinfo_bundle = null;

int gretl_bundle_is_stacked(GretlBundle * b)
{
    UserVar * u = get_user_var_by_data(b);

    return u != null;
}

/* gets the number of keys in the bundle's hash table */

int gretl_bundle_get_n_keys(GretlBundle * b)
{
    int n_items = 0;

    if (b != null && b.ht != null) {
        n_items = g_hash_table_size(b.ht);
    }

    return n_items;
}

/* gets total number of members including any "special"
   contents outside of the hash table */

int gretl_bundle_get_n_members(GretlBundle * b)
{
    int nmemb = 0;

    if (b != null) {
        if (b.type == BundleType.kalman) {
            nmemb += kalman_bundle_n_members(b);
        }
        if (b.ht != null) {
            nmemb += g_hash_table_size(b.ht);
        }
    }

    return nmemb;
}

void maybe_append_list(gpointer key, gpointer value, gpointer p)
{
    BundledItem * item = cast(BundledItem *) value;
    GList **plist = cast(GList **) p;

    if (item.type == GretlType.list) {
        *plist = g_list_append(*plist, item.data);
    }
}

/* If bundle @b contains any members of type list, compose
   and return a GList that contains all of them.
*/

GList * gretl_bundle_get_lists(GretlBundle * b)
{
    GList *list = null;

    g_hash_table_foreach(b.ht, &maybe_append_list, &list);

    return list;
}

int gretl_bundle_has_content(GretlBundle * b)
{
    int ret = 0;

    if (b != null && b.ht != null &&
        (b.type == BundleType.kalman ||
         g_hash_table_size(b.ht) > 0)) {
        ret = 1;
    }

    return ret;
}

int type_can_be_bundled(GretlType type)
{
    if (type == GretlType.int_ ||
        type == GretlType.uint32 ||
        type == GretlType.uint64 ||
        type == GretlType.bool_) {
        /* for the current purpose */
        type = GretlType.double_;
    }

    return (type == GretlType.double_ ||
            type == GretlType.string_ ||
            type == GretlType.matrix ||
            type == GretlType.series ||
            type == GretlType.bundle ||
            type == GretlType.array ||
            type == GretlType.list);
}

bool bundled_scalar(GretlType t)
{
    return (t == GretlType.double_ || t == GretlType.int_ || t == GretlType.uint32 || t == GretlType.uint64);
}

bool is_virtual_object(const(char) * s)
{
    return (s != null && strlen(s) > 15 && (strncmp(s, "virtual_series:", 15) == 0 || strncmp(s, "virtual_matrix:", 15) == 0));
}

GretlType virtual_type(const(char) * s)
{
    if (!strncmp(s + 8, "series", 6)) {
    return GretlType.series;
    } else if (!strncmp(s + 8, "matrix", 6)) {
    return GretlType.matrix;
    } else {
    return GretlType.none;
    }
}

const(char) * get_gen_str(const(char) * s)
{
    return s + 15;
}

int bundled_item_copy_in_data(BundledItem * item, void * ptr)
{
    int err = 0;

    if (item.is_virtual) {
    item.data = gretl_strdup(cast(char *) ptr);
    return 0;
    }

    switch (item.type) {
    case GretlType.double_:
        item.data = cast(typeof(item.data)) malloc(double.sizeof);
        if (item.data != null) {
            double * dp = cast(double *) item.data;

            *dp = *cast(double *) ptr;
        }
        break;
    case GretlType.int_:
        item.data = cast(typeof(item.data)) malloc(int.sizeof);
        if (item.data != null) {
            int * ip = cast(int *) item.data;

            *ip = *cast(int *) ptr;
        }
        break;
    case GretlType.uint32:
        item.data = cast(typeof(item.data)) malloc(uint.sizeof);
        if (item.data != null) {
            uint *up = cast(uint *) item.data;

            *up = *cast(uint *) ptr;
        }
        break;
    case GretlType.uint64:
        item.data = cast(typeof(item.data)) malloc(ulong.sizeof);
        if (item.data != null) {
            ulong *up = cast(ulong *) item.data;

            *up = *cast(ulong *) ptr;
        }
        break;
    case GretlType.string_:
        item.data = gretl_strdup(cast(char *) ptr);
        break;
    case GretlType.matrix:
    case GretlType.series:
    item.data = gretl_matrix_copy(cast(GretlMatrixRaw *) ptr);
        break;
    case GretlType.list:
        item.data = gretl_list_copy(cast(const(int) *) ptr);
        break;
    case GretlType.bundle:
        item.data = gretl_bundle_copy(cast(GretlBundle *) ptr, &err);
        break;
    case GretlType.array:
        item.data = gretl_array_copy(cast(GretlArray *) ptr, &err);
        break;
    default:
        err = GretlError.types;
        break;
    }

    if (!err && item.data == null) {
        err = GretlError.alloc;
    }

    return err;
}

/* allocate and fill out a 'value' (type plus data pointer) that will
   be inserted into a bundle's hash table */

BundledItem * bundled_item_new(GretlType type, void * ptr, int copy, const(char) * note, int * err)
{
    BundledItem * item = cast(BundledItem *) malloc(BundledItem.sizeof);

    if (item == null) {
        *err = GretlError.alloc;
        return null;
    }

    item.type = type;
    item.key = null;
    item.note = null;
    item.is_virtual = 0;

    if (item.type == GretlType.string_ && is_virtual_object(cast(const(char) *) ptr)) {
    item.is_virtual = 1;
    item.type = virtual_type(cast(const(char) *) ptr);
    }

    if (!copy && !bundled_scalar(item.type)) {
        item.data = ptr;
    } else {
        *err = bundled_item_copy_in_data(item, ptr);
        if (*err) {
            free(item);
            item = null;
        }
    }

    if (item != null && note != null) {
        item.note = gretl_strdup(note);
    }

    return item;
}

void bundled_item_free_data(BundledItem * item)
{
    if (item.is_virtual) {
    free(item.data);
    } else if (item.type == GretlType.matrix ||
           item.type == GretlType.series) {
        gretl_matrix_free(cast(GretlMatrixRaw *) item.data);
    } else if (item.type == GretlType.bundle) {
        gretl_bundle_destroy(cast(GretlBundle *) item.data);
    } else if (item.type == GretlType.array) {
        gretl_array_destroy(cast(GretlArray *) item.data);
    } else {
        free(item.data);
    }
}

/* Note: we come here only if the replacement type is the same as the
   original type, apart from the case of inter-conversion of the various
   scalar types.
*/

int bundled_item_replace_data(BundledItem * item, void * ptr, GretlType src_t, int copy)
{
    int err = 0;

    if (ptr == item.data) {
        return 0;
    }

    if (bundled_scalar(item.type)) {
        /* storage is already available */
        if (item.type == GretlType.double_) {
            double * dp = cast(double *) item.data;

            if (src_t == GretlType.double_) {
                *dp = *cast(double *) ptr;
            } else if (src_t == GretlType.int_) {
                *dp = cast(double) *cast(int *) ptr;
            } else {
                *dp = cast(double) *cast(uint *) ptr;
            }
        } else if (item.type == GretlType.int_) {
            int * ip = cast(int *) item.data;

            if (src_t == GretlType.int_) {
                *ip = *cast(int *) ptr;
            } else if (src_t == GretlType.double_) {
                *ip = gretl_int_from_double(*cast(double *) ptr, &err);
            } else {
                *ip = cast(int) *cast(uint *) ptr;
            }
        } else if (item.type == GretlType.uint32) {
            uint *up = cast(uint *) item.data;

            if (src_t == GretlType.uint32) {
                *up = *cast(uint *) ptr;
            } else if (src_t == GretlType.double_) {
                *up = gretl_uint32_from_double(*cast(double *) ptr, &err);
            } else {
                *up = cast(uint) *cast(int *) ptr;
            }
        } else if (item.type == GretlType.uint64) {
            ulong *up = cast(ulong *) item.data;

            if (src_t == GretlType.uint64) {
                *up = *cast(ulong *) ptr;
            } else if (src_t == GretlType.double_) {
                *up = gretl_uint64_from_double(*cast(double *) ptr, &err);
            } else {
                *up = cast(ulong) *cast(int *) ptr;
            }
        }
    } else {
        /* free then copy or donate */
        bundled_item_free_data(item);
        if (copy) {
            err = bundled_item_copy_in_data(item, ptr);
        } else {
            item.data = ptr;
        }
    }

    if (!err && item.data == null) {
        err = GretlError.alloc;
    }

    if (item.note != null) {
        free(item.note);
        item.note = null;
    }

    return err;
}

/* callback invoked when a bundle's hash table is destroyed */

void bundle_item_destroy(gpointer data)
{
    BundledItem * item = cast(BundledItem *) data;

    bundled_item_free_data(item);
    g_free(item.key);
    free(item.note);
    free(item);
}

/**
 * gretl_bundle_destroy:
 * @bundle: bundle to destroy.
 *
 * Frees all contents of @bundle as well as the pointer itself.
 */

void gretl_bundle_destroy(GretlBundle * bundle)
{
    if (bundle != null) {
        if (bundle.ht != null) {
            g_hash_table_destroy(bundle.ht);
        }
        free(bundle.creator);
        if (bundle.type == BundleType.kalman) {
            kalman_free(cast(Kalman *) bundle.data);
        }
        free(bundle);
    }
}

/**
 * gretl_bundle_void_content:
 * @bundle: target bundle.
 *
 * Frees all contents of @bundle.
 */

void gretl_bundle_void_content(GretlBundle * bundle)
{
    if (bundle == null) {
        return;
    }

    if (bundle.creator != null) {
        free(bundle.creator);
        bundle.creator = null;
    }

    if (bundle.ht != null && g_hash_table_size(bundle.ht) > 0) {
        g_hash_table_destroy(bundle.ht);
        bundle.ht = g_hash_table_new_full(&g_str_hash, &g_str_equal,
                                           null, &bundle_item_destroy);
    }

    if (bundle.type == BundleType.kalman) {
        kalman_free(cast(Kalman *) bundle.data);
        bundle.data = null;
        bundle.type = BundleType.plain;
    }
}

/**
 * gretl_bundle_new:
 *
 * Returns: a newly allocated, empty gretl bundle.
 */

GretlBundle * gretl_bundle_new()
{
    GretlBundle * b = cast(GretlBundle *) malloc(GretlBundle.sizeof);

    if (b != null) {
        b.type = BundleType.plain;
        b.ht = g_hash_table_new_full(&g_str_hash, &g_str_equal,
                                      null, &bundle_item_destroy);
        b.creator = null;
        b.data = null;
    }

    return b;
}

/* Determine whether @name is the name of a saved bundle. */

int gretl_is_bundle(const(char) * name)
{
    return get_user_var_of_type_by_name(name, GretlType.bundle) != null;
}

/**
 * get_bundle_by_name:
 * @name: the name to look up.
 *
 * Returns: pointer to a saved bundle, if found, else NULL.
 */

GretlBundle * get_bundle_by_name(const(char) * name)
{
    GretlBundle * b = null;

    if (name != null && *name != '\0') {
        UserVar * u = get_user_var_of_type_by_name(name, GretlType.bundle);

        if (u != null) {
            b = cast(GretlBundle *) user_var_get_value(u);
        }
    }

    return b;
}

int gretl_bundle_has_data(GretlBundle * b, const(char) * key)
{
    gpointer p = g_hash_table_lookup(b.ht, key);

    return (p != null);
}

/* Supports gretl_bundle_get_data(), gretl_bundle_get_element()
   and gretl_bundle_get_target().
*/

void * real_bundle_get_data(GretlBundle * bundle, const(char) * key, GretlType * type, int * is_virtual, int get_target, int * err)
{
    void * ret = null;
    int private_ = 0;
    int myerr = 0;

    if (bundle == null) {
        myerr = GretlError.data;
        goto finish;
    }

    if (!get_target && bundle.type == BundleType.kalman) {
        ret = maybe_retrieve_kalman_element(bundle.data, key,
                                            type, &private_,
                                            &myerr);
    }

    if (!myerr && ret == null && !private_) {
        /* try for a regular bundle member */
        gpointer p = g_hash_table_lookup(bundle.ht, key);

        if (p != null) {
            BundledItem * item = cast(BundledItem *) p;

        ret = item.data;
        if (type != null) {
        *type = item.type;
        }
        if (is_virtual != null) {
        *is_virtual = item.is_virtual;
        }
        } else {
            if (err != null) {
                gretl_errmsg_sprintf("\"%s\": %s".ptr, key, _("no such item").ptr);
            }
            myerr = GretlError.data;
        }
    }

 finish:

    if (err != null) {
        *err = myerr;
    }

    return ret;
}

/**
 * gretl_bundle_get_element:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @type: location to receive data type, or NULL.
 * @is_virtual: location to receive "is_virtual" flag.
 * @err: location to receive error code, or NULL.
 *
 * Returns: the item pointer associated with @key in the
 * specified @bundle, or NULL if there is no such item.
 * If @err is non-NULL, its content is set to a non-zero
 * value if @bundle contains no item with key @key. If
 * the intent is simply to determine whether @bundle contains
 * an item under the specified @key, @err should generally be
 * left NULL.
 */

void * gretl_bundle_get_element(GretlBundle * bundle, const(char) * key, GretlType * type, int * is_virtual, int * err)
{
    return real_bundle_get_data(bundle, key, type,
                                is_virtual, 0, err);
}

/**
 * gretl_bundle_get_data:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @type: location to receive data type, or NULL.
 * @err: location to receive error code, or NULL.
 *
 * Returns: the item pointer associated with @key in the
 * specified @bundle, or NULL if there is no such item.
 * If @err is non-NULL, its content is set to a non-zero
 * value if @bundle contains no item with key @key. If
 * the intent is simply to determine whether @bundle contains
 * an item under the specified @key, @err should generally be
 * left NULL.
 *
 * Note that the value returned is the actual data pointer from
 * within the bundle, not a copy of the data; so the pointer
 * must not be freed, and in general its content should not
 * be modified.
 */

void * gretl_bundle_get_data(GretlBundle * bundle, const(char) * key, GretlType * type, int * err)
{
    return real_bundle_get_data(bundle, key, type,
                                null, 0, err);
}

/**
 * gretl_bundle_get_data_full:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @type: location to receive data type, or NULL.
 * @is_virtual: location to receive is_virtual status, or NULL.
 * @err: location to receive error code, or NULL.
 *
 * Returns: the item pointer associated with @key in the
 * specified @bundle, or NULL if there is no such item.
 * If @err is non-NULL, its content is set to a non-zero
 * value if @bundle contains no item with key @key. If
 * the intent is simply to determine whether @bundle contains
 * an item under the specified @key, @err should generally be
 * left NULL.
 *
 * Note that the value returned is the actual data pointer from
 * within the bundle, not a copy of the data; so the pointer
 * must not be freed, and in general its content should not
 * be modified.
 */

void * gretl_bundle_get_data_full(GretlBundle * bundle, const(char) * key, GretlType * type, int * is_virtual, int * err)
{
    return real_bundle_get_data(bundle, key, type,
                                is_virtual, 0, err);
}

/* As gretl_bundle_get_data() except that the caller is looking
   for a "target" object within @bundle with a view to modifying
   it. Called only from geneval.c, in response to, for example,
   b.X = Y.
*/

void * gretl_bundle_get_target(GretlBundle * bundle, const(char) * key, GretlType * type, int * err)
{
    /* FIXME virtual member? */
    return real_bundle_get_data(bundle, key, type,
                                null, 1, err);
}

/**
 * gretl_bundle_steal_data:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @type: location to receive data type, or NULL.
 * @err: location to receive error code, or NULL.
 *
 * Works like gretl_bundle_get_data() except that the data
 * pointer in question is removed from @bundle before it is
 * returned to the caller; so in effect the caller assumes
 * ownership of the item.
 *
 * Returns: the item pointer associated with @key in the
 * specified @bundle, or NULL if there is no such item.
 */

void * gretl_bundle_steal_data(GretlBundle * bundle, const(char) * key, GretlType * type, int * err)
{
    void * ret = null;

    if (bundle == null) {
        if (err != null) {
            *err = GretlError.data;
        }
    } else {
        gpointer p = g_hash_table_lookup(bundle.ht, key);

        if (p != null) {
            GList *keys = g_hash_table_get_keys(bundle.ht);
            BundledItem * item = cast(BundledItem *) p;
            char * keycpy = null;

        if (item.is_virtual) {
        if (err != null) {
            *err = GretlError.invarg;
        }
        return null;
        }

            ret = item.data;
            if (type != null) {
                *type = item.type;
            }
            while (keys) {
                if (!strcmp(cast(const(char) *) keys.data, key)) {
                    keycpy = cast(char *) keys.data;
                    break;
                } else if (keys.next) {
                    keys = keys.next;
                } else {
                    break;
                }
            }
            g_hash_table_steal(bundle.ht, key);
        /* g_hash_table_steal() frees neither the data nor
           the associated key: here we free the key manually,
           having saved a pointer to it.
        */
            g_free(keycpy);
            g_list_free(keys);
            free(item);
        } else if (err != null) {
            gretl_errmsg_sprintf("\"%s\": %s".ptr, key, _("no such item").ptr);
            *err = GretlError.data;
        }
    }

    return ret;
}

void * gretl_bundle_get_private_data(GretlBundle * bundle)
{
    return bundle.data;
}

BundleType gretl_bundle_get_type(GretlBundle * bundle)
{
    return bundle.type;
}

/**
 * gretl_bundle_get_member_type:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @err:location to receive error code, or NULL.
 *
 * Returns: the data type associated with @key in the
 * specified @bundle, or 0 on failure. Set the @err
 * argument to NULL if you do not want an error flagged
 * in case there's no such key in @bundle.
 */

GretlType gretl_bundle_get_member_type(GretlBundle * bundle, const(char) * key, int * err)
{
    GretlType ret = GretlType.none;
    int private_ = 0;
    int myerr = 0;

    if (bundle == null) {
        myerr = GretlError.data;
    } else if (bundle.type == BundleType.kalman) {
        maybe_retrieve_kalman_element(bundle.data, key,
                                      &ret, &private_,
                                      &myerr);
    }

    if (!myerr && ret == GretlType.none && !private_) {
        /* try for a regular bundle member */
        gpointer p = g_hash_table_lookup(bundle.ht, key);

        if (p != null) {
            BundledItem * item = cast(BundledItem *) p;

            ret = item.type;
        } else if (err != null) {
            gretl_errmsg_sprintf("\"%s\": %s".ptr, key, _("no such item").ptr);
            myerr = GretlError.data;
        }
    }

    if (err != null) {
        *err = myerr;
    }

    return ret;
}

/**
 * gretl_bundle_has_key:
 * @bundle: bundle to access.
 * @key: name of key to test.
 *
 * Returns: 1 if there is an item under the given @key in the
 * specified @bundle, 0 otherwise.
 */

int gretl_bundle_has_key(GretlBundle * bundle, const(char) * key)
{
    int ret = 0;

    if (bundle != null && key != null) {
        gpointer p = g_hash_table_lookup(bundle.ht, key);

        ret = (p != null);
    }

    return ret;
}

/**
 * gretl_bundle_get_matrix:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @err: location to receive error code.
 *
 * Returns: the matrix associated with @key in the
 * specified @bundle, if any; otherwise NULL.
 */

GretlMatrixRaw * gretl_bundle_get_matrix(GretlBundle * bundle, const(char) * key, int * err)
{
    GretlMatrixRaw * m = null;
    GretlType type;
    void * ptr;
    int myerr = 0;

    ptr = gretl_bundle_get_data(bundle, key, &type, &myerr);

    if (ptr == null) {
    myerr = GretlError.invarg;
    } else if (type != GretlType.matrix) {
    myerr = GretlError.types;
    } else {
    m = cast(GretlMatrixRaw *) ptr;
    }

    if (err != null) {
        *err = myerr;
    }

    return m;
}

/**
 * gretl_bundle_get_array:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @err: location to receive error code.
 *
 * Returns: the array associated with @key in the
 * specified @bundle, if any; otherwise NULL.
 */

GretlArray * gretl_bundle_get_array(GretlBundle * bundle, const(char) * key, int * err)
{
    GretlArray * a = null;
    GretlType type;
    void * ptr;
    int myerr = 0;

    ptr = gretl_bundle_get_data(bundle, key, &type, &myerr);

    if (ptr == null) {
    myerr = GretlError.invarg;
    } else if (type != GretlType.array) {
        myerr = GretlError.types;
    } else {
    a = cast(GretlArray *) ptr;
    }

    if (err != null) {
        *err = myerr;
    }

    return a;
}

/**
 * gretl_bundle_get_bundle:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @err: location to receive error code.
 *
 * Returns: the bundle associated with @key in the
 * specified @bundle, if any; otherwise NULL.
 */

GretlBundle * gretl_bundle_get_bundle(GretlBundle * bundle, const(char) * key, int * err)
{
    GretlBundle * b = null;
    GretlType type;
    void * ptr;
    int myerr = 0;

    ptr = gretl_bundle_get_data(bundle, key, &type, &myerr);

    if (ptr == null) {
    myerr = GretlError.invarg;
    } else if (type != GretlType.bundle) {
        myerr = GretlError.types;
    } else {
    b = cast(GretlBundle *) ptr;
    }

    if (err != null) {
        *err = myerr;
    }

    return b;
}

/**
 * gretl_bundle_get_series:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @err: location to receive error code.
 *
 * Returns: the series associated with @key in the
 * specified @bundle, if any, in the form of a column
 * vector; otherwise NULL.
 */

GretlMatrixRaw * gretl_bundle_get_series(GretlBundle * bundle, const(char) * key, int * err)
{
    GretlMatrixRaw * m = null;
    GretlType type;
    int is_virtual = 0;
    void * ptr;
    int myerr = 0;

    ptr = gretl_bundle_get_data_full(bundle, key, &type,
                     &is_virtual, &myerr);

    if (ptr == null || is_virtual) {
    myerr = GretlError.invarg;
    } else if (type != GretlType.series) {
        myerr = GretlError.types;
    } else {
        m = cast(GretlMatrixRaw *) ptr;
    }

    if (err != null) {
        *err = myerr;
    }

    return m;
}

/**
 * gretl_bundle_get_list:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @err: location to receive error code.
 *
 * Returns: the list associated with @key in the
 * specified @bundle, if any; otherwise NULL.
 */

int * gretl_bundle_get_list(GretlBundle * bundle, const(char) * key, int * err)
{
    int * list = null;
    GretlType type;
    void * ptr;
    int myerr = 0;

    ptr = gretl_bundle_get_data(bundle, key, &type, &myerr);

    if (ptr == null) {
    myerr = GretlError.invarg;
    } else if (type != GretlType.list) {
        myerr = GretlError.types;
    } else {
    list = cast(int *) ptr;
    }

    if (err != null) {
        *err = myerr;
    }

    return list;
}

/**
 * gretl_bundle_get_scalar:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @err: location to receive error code, or NULL.
 *
 * Returns: the scalar value associated with @key in the
 * specified @bundle, if any; otherwise #NADBL.
 */

double gretl_bundle_get_scalar(GretlBundle * bundle, const(char) * key, int * err)
{
    double x = NADBL;
    GretlType type;
    void * ptr;
    int myerr = 0;

    ptr = gretl_bundle_get_data(bundle, key, &type, &myerr);

    if (ptr == null) {
        myerr = GretlError.invarg;
    } else if (type != GretlType.double_ &&
               type != GretlType.int_ &&
               type != GretlType.uint32 &&
               type != GretlType.uint64 &&
               type != GretlType.matrix) {
        myerr = GretlError.types;
    }

    if (myerr == 0) {
        if (type == GretlType.double_) {
            x = *cast(double *) ptr;
        } else if (type == GretlType.int_) {
            int k = *cast(int *) ptr;

            x = cast(double) k;
        } else if (type == GretlType.uint32) {
            uint u = *cast(uint *) ptr;

            x = cast(double) u;
        } else if (type == GretlType.uint64) {
            ulong u = *cast(ulong *) ptr;

            x = cast(double) u;
        } else {
            /* must be a matrix */
            GretlMatrixRaw * m = cast(GretlMatrixRaw *) ptr;

            if (gretl_matrix_is_scalar(m)) {
                x = m.val[0];
            } else {
                myerr = GretlError.types;
            }
        }
    }

    if (err != null) {
        *err = myerr;
    }

    return x;
}

/**
 * gretl_bundle_get_int:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @err: location to receive error code.
 *
 * Returns: the integer value associated with @key in the
 * specified @bundle, if any; otherwise 0.
 */

int gretl_bundle_get_int(GretlBundle * bundle, const(char) * key, int * err)
{
    int i = 0;
    GretlType type;
    void * ptr;
    int myerr = 0;

    ptr = gretl_bundle_get_data(bundle, key, &type, &myerr);

    if (ptr == null) {
    myerr = GretlError.invarg;
    } else {
        if (type == GretlType.int_) {
            int * pi = cast(int *) ptr;

            i = *pi;
        } else if (type == GretlType.uint32) {
            uint u;
            uint * pu = cast(uint *) ptr;

            u = *pu;
            if (u <= int.max) {
                i = cast(int) u;
            } else {
                myerr = GretlError.types;
            }
        } else if (type == GretlType.uint64) {
            ulong u;
            ulong * pu = cast(ulong *) ptr;

            u = *pu;
            if (u <= int.max) {
                i = cast(int) u;
            } else {
                myerr = GretlError.types;
            }
        } else if (type == GretlType.double_) {
            double x;
            double * px = cast(double *) ptr;

            x = *px;
            if (x >= int.min && x <= int.max && x == floor(x)) {
                i = cast(int) x;
            } else {
                myerr = GretlError.types;
            }
        } else {
            myerr = GretlError.types;
        }
    }

    if (err != null) {
        *err = myerr;
    }

    return i;
}

/**
 * gretl_bundle_get_int_deflt:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @deflt: default integer value.
 *
 * Returns: the integer value associated with @key in the
 * specified @bundle, if any; otherwise @deflt.
 */

int gretl_bundle_get_int_deflt(GretlBundle * bundle, const(char) * key, int deflt)
{
    int val = deflt;
    GretlType type;
    void * ptr;

    ptr = gretl_bundle_get_data(bundle, key, &type, null);

    if (ptr != null) {
        if (type == GretlType.int_) {
            int * pi = cast(int *) ptr;

            val = *pi;
        } else if (type == GretlType.uint32) {
            uint *pu = cast(uint *) ptr;

            val = cast(int) *pu;
        } else if (type == GretlType.uint64) {
            ulong *pu = cast(ulong *) ptr;

            val = cast(int) *pu;
        } else if (type == GretlType.double_) {
            double * px = cast(double *) ptr;

            val = cast(int) *px;
        }
    }

    return val;
}

/**
 * gretl_bundle_get_bool:
 * @bundle: bundle to access.
 * @key: name of key to access, if present.
 * @deflt: default boolean value.
 *
 * Returns: the boolean value associated with @key in the
 * specified @bundle, if any; otherwise @deflt.
 */

int gretl_bundle_get_bool(GretlBundle * bundle, const(char) * key, int deflt)
{
    int val = deflt;
    GretlType type;
    void * ptr;

    ptr = gretl_bundle_get_data(bundle, key, &type, null);

    if (ptr != null) {
        if (type == GretlType.int_) {
            int * pi = cast(int *) ptr;

            val = (*pi != 0);
        } else if (type == GretlType.uint32) {
            uint *pu = cast(uint *) ptr;

            val = (*pu != 0);
        } else if (type == GretlType.uint32) {
            ulong *pu = cast(ulong *) ptr;

            val = (*pu != 0);
        } else if (type == GretlType.double_) {
            double * px = cast(double *) ptr;

            val = (*px != 0);
        }
    }

    return val;
}

/**
 * gretl_bundle_get_uint32:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @err: location to receive error code.
 *
 * Returns: the 32-bit unsigned integer value associated with @key
 * in the specified @bundle, if any; otherwise 0.
 */

uint gretl_bundle_get_uint32(GretlBundle * bundle, const(char) * key, int * err)
{
    uint u = 0;
    GretlType type;
    void * ptr;
    int myerr = 0;

    ptr = gretl_bundle_get_data(bundle, key, &type, &myerr);

    if (ptr == null) {
    myerr = GretlError.invarg;
    } else {
        if (type == GretlType.uint32) {
            uint *pu = cast(uint *) ptr;

            u = *pu;
        } else if (type == GretlType.int_) {
            int i;
            int * pi = cast(int *) ptr;

            i = *pi;
            if (i >= 0) {
                u = cast(uint) i;
            } else {
                myerr = GretlError.types;
            }
        } else if (type == GretlType.double_) {
            double x;
            double * px = cast(double *) ptr;

            x = *px;
            if (x >= 0 && x <= uint.max && x == floor(x)) {
                u = cast(uint) x;
            } else {
                myerr = GretlError.types;
            }
        } else {
            myerr = GretlError.types;
        }
    }

    if (err != null) {
        *err = myerr;
    }

    return u;
}

/**
 * gretl_bundle_get_uint64:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @err: location to receive error code.
 *
 * Returns: the 64-bit unsigned integer value associated with @key
 * in the specified @bundle, if any; otherwise 0.
 */

ulong gretl_bundle_get_uint64(GretlBundle * bundle, const(char) * key, int * err)
{
    ulong u = 0;
    GretlType type;
    void * ptr;
    int myerr = 0;

    ptr = gretl_bundle_get_data(bundle, key, &type, &myerr);

    if (ptr == null) {
    myerr = GretlError.invarg;
    } else {
        if (type == GretlType.uint64) {
            ulong *pu = cast(ulong *) ptr;

            u = *pu;
        } else if (type == GretlType.uint32) {
            uint *pu = cast(uint *) ptr;

            u = cast(ulong) *pu;
        } else if (type == GretlType.int_) {
            int i;
            int * pi = cast(int *) ptr;

            i = *pi;
            if (i >= 0) {
                u = cast(ulong) i;
            } else {
                myerr = GretlError.types;
            }
        } else if (type == GretlType.double_) {
            double x;
            double * px = cast(double *) ptr;

            x = *px;
            if (x >= 0 && x <= cast(double) ulong.max && x == floor(x)) {
                u = cast(ulong) x;
            } else {
                myerr = GretlError.types;
            }
        } else {
            myerr = GretlError.types;
        }
    }

    if (err != null) {
        *err = myerr;
    }

    return u;
}

/**
 * gretl_bundle_get_string:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @err: location to receive error code.
 *
 * Returns: the string value associated with @key in the
 * specified @bundle, if any; otherwise NULL.
 */

const(char) * gretl_bundle_get_string(GretlBundle * bundle, const(char) * key, int * err)
{
    const(char) * s = null;
    GretlType type;
    int is_virtual = 0;
    void * ptr;
    int myerr = 0;

    ptr = gretl_bundle_get_data_full(bundle, key, &type,
                     &is_virtual, &myerr);

    if (is_virtual) {
    s = get_gen_str(cast(const(char) *) ptr);
    } else if (type == GretlType.string_) {
    s = cast(const(char) *) ptr;
    } else {
    myerr = GretlError.types;
    }

    if (err != null) {
        *err = myerr;
    }

    return s;
}

/**
 * gretl_bundle_get_strings:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @ns: location to receive the number of strings, if any.
 *
 * Returns: the array of strings associated with @key in the
 * specified @bundle, if any; otherwise NULL.
 */

const(char) ** gretl_bundle_get_strings(GretlBundle * bundle, const(char) * key, int * ns)
{
    const(char) ** S = null;
    GretlType type = cast(GretlType) 0;
    void * ptr;

    ptr = gretl_bundle_get_data(bundle, key, &type, null);

    if (ptr != null && type == GretlType.array) {
        GretlArray * a = cast(GretlArray *) ptr;

        S = cast(const(char) **) gretl_array_get_strings(a, ns);
    }

    return S;
}

/**
 * gretl_bundle_get_note:
 * @bundle: bundle to access.
 * @key: name of key to access.
 *
 * Returns: the note associated with @key in the
 * specified @bundle, if any; otherwise NULL.
 */

const(char) * gretl_bundle_get_note(GretlBundle * bundle, const(char) * key)
{
    const(char) * ret = null;

    if (bundle != null) {
        gpointer p = g_hash_table_lookup(bundle.ht, key);

        if (p != null) {
            BundledItem * item = cast(BundledItem *) p;

            ret = item.note;
        }
    }

    return ret;
}

/**
 * gretl_bundle_get_creator:
 * @bundle: bundle to access.
 *
 * Returns: the name of the package that created @bundle, if any,
 * otherwise NULL.
 */

const(char) * gretl_bundle_get_creator(GretlBundle * bundle)
{
    return (bundle != null)? bundle.creator : null;
}

/**
 * gretl_bundle_get_content:
 * @bundle: bundle to access.
 *
 * Returns: the content of @bundle, which is in fact
 * a GHashTable object.
 */

void * gretl_bundle_get_content(GretlBundle * bundle)
{
    return (bundle == null)? null : cast(void *) bundle.ht;
}

/**
 * gretl_bundles_swap_content:
 * @b1: first bundle.
 * @b2: second bundle.
 *
 * Exchanges the entire contents of the two bundles.
 *
 * Returns: 0 on success, error code on failure.
 */

int gretl_bundles_swap_content(GretlBundle * b1, GretlBundle * b2)
{
    if (b1 == null || b1.type != BundleType.plain ||
        b2 == null || b2.type != BundleType.plain) {
        return GretlError.data;
    } else {
        void * tmp = b1.ht;

        b1.ht = b2.ht;
        b2.ht = cast(GHashTable *) tmp;
        return 0;
    }
}

int real_bundle_set_data(GretlBundle * b, const(char) * key, void * ptr, GretlType type, int copy, const(char) * note)
{
    int err = 0, done = 0;

    if (key == null || key[0] == '\0') {
        gretl_errmsg_sprintf("real_bundle_set_data: missing key string".ptr);
        return GretlError.data;
    }

    if (b.type == BundleType.kalman) {
        done = maybe_set_kalman_element(b.data, key,
                                        ptr, type, copy,
                                        &err);
    }

    if (!done && !err) {
        BundledItem * item = cast(BundledItem *) g_hash_table_lookup(b.ht, key);
        int replace = 0;

        if (item != null) {
            replace = 1;
            if (item.type == type && !item.is_virtual) {
                /* we can take a shortcut */
                return bundled_item_replace_data(item, ptr, type, copy);
            }
        }

        item = bundled_item_new(type, ptr, copy, note, &err);
        if (!err) {
            item.key = g_strdup(key);
            if (replace) {
                g_hash_table_replace(b.ht, item.key, item);
            } else {
                g_hash_table_insert(b.ht, item.key, item);
            }
        }
    }

    return err;
}

/**
 * gretl_bundle_donate_data:
 * @bundle: target bundle.
 * @key: name of key to create or replace.
 * @ptr: data pointer.
 * @type: type of data.
 *
 * Sets the data type and pointer to be associated with @key in
 * the bundle given by @name. If @key is already present in
 * the bundle's hash table the original value is replaced
 * and destroyed. The value of @ptr is transcribed into the
 * bundle, which therefore "takes ownership" of the data;
 * compare gretl_bundle_set_data().
 *
 * Returns: 0 on success, error code on error.
 */

int gretl_bundle_donate_data(GretlBundle * bundle, const(char) * key, void * ptr, GretlType type)
{
    if (key != null && ptr == null) {
        gretl_errmsg_sprintf("'%s': got NULL data value".ptr, key);
        return GretlError.data;
    }

    return real_bundle_set_data(bundle, key, ptr, type, 0, null);
}

/**
 * gretl_bundle_set_data:
 * @bundle: target bundle.
 * @key: name of key to create or replace.
 * @ptr: data pointer.
 * @type: the gretl type of the data.
 *
 * Sets the data type and pointer to be associated with @key in
 * the bundle given by @name. If @key is already present in
 * the bundle's hash table the original value is replaced
 * and destroyed. The content of @ptr is copied into the
 * bundle; compare gretl_bundle_donate_data().
 *
 * Returns: 0 on success, error code on error.
 */

int gretl_bundle_set_data(GretlBundle * bundle, const(char) * key, void * ptr, GretlType type)
{
    int err;

    if (type == GretlType.series) {
    fprintf(stderr, "gretl_bundle_set_data: got GRETL_TYPE_SERIES\n".ptr);
    fprintf(stderr, " use gretl_bundle_set_series() instead\n".ptr);
    return GretlError.types;
    }

    if (bundle == null) {
        err = GretlError.unkvar;
    } else {
        err = real_bundle_set_data(bundle, key, ptr, type, 1, null);
    }

    return err;
}

/**
 * gretl_bundle_set_string:
 * @bundle: target bundle.
 * @key: name of key to create or replace.
 * @str: the string to set.
 *
 * Sets @str as a member of @bundle under the name @key.
 * If @key is already present in the bundle the original
 * value is replaced and destroyed.
 *
 * Returns: 0 on success, error code on error.
 */

int gretl_bundle_set_string(GretlBundle * bundle, const(char) * key, const(char) * str)
{
    return gretl_bundle_set_data(bundle, key, cast(void *) str, GretlType.string_);
}

/**
 * gretl_bundle_set_scalar:
 * @bundle: target bundle.
 * @key: name of key to create or replace.
 * @val: the value to set.
 *
 * Sets @val as a member of @bundle under the name @key.
 * If @key is already present in the bundle the original
 * value is replaced and destroyed.
 *
 * Returns: 0 on success, error code on error.
 */

int gretl_bundle_set_scalar(GretlBundle * bundle, const(char) * key, double val)
{
    return gretl_bundle_set_data(bundle, key, &val, GretlType.double_);
}

/**
 * gretl_bundle_set_int:
 * @bundle: target bundle.
 * @key: name of key to create or replace.
 * @val: the integer value to set.
 *
 * Sets @val as a member of @bundle under the name @key.
 * If @key is already present in the bundle the original
 * value is replaced and destroyed.
 *
 * Returns: 0 on success, error code on error.
 */

int gretl_bundle_set_int(GretlBundle * bundle, const(char) * key, int val)
{
    return gretl_bundle_set_data(bundle, key, &val, GretlType.int_);
}

/**
 * gretl_bundle_set_uint32:
 * @bundle: target bundle.
 * @key: name of key to create or replace.
 * @val: the 32-bit unsigned integer value to set.
 *
 * Sets @val as a member of @bundle under the name @key.
 * If @key is already present in the bundle the original
 * value is replaced and destroyed.
 *
 * Returns: 0 on success, error code on error.
 */

int gretl_bundle_set_uint32(GretlBundle * bundle, const(char) * key, uint val)
{
    return gretl_bundle_set_data(bundle, key, &val, GretlType.uint32);
}

/**
 * gretl_bundle_set_uint64:
 * @bundle: target bundle.
 * @key: name of key to create or replace.
 * @val: the 64-bit unsigned integer value to set.
 *
 * Sets @val as a member of @bundle under the name @key.
 * If @key is already present in the bundle the original
 * value is replaced and destroyed.
 *
 * Returns: 0 on success, error code on error.
 */

int gretl_bundle_set_uint64(GretlBundle * bundle, const(char) * key, ulong val)
{
    return gretl_bundle_set_data(bundle, key, &val, GretlType.uint64);
}

/**
 * gretl_bundle_set_series:
 * @bundle: target bundle.
 * @key: name of key to create or replace.
 * @x: array of doubles.
 * @dset: the current dataset.
 *
 * Copies @x as a member of @bundle under the name @key.
 * If @key is already present in the bundle the original
 * value is replaced and destroyed.
 *
 * Returns: 0 on success, error code on error.
 */

int gretl_bundle_set_series(GretlBundle * bundle, const(char) * key, const(double) * x, const(DatasetRaw) * dset)
{
    GretlMatrixRaw * m;
    int len = dset.n;

    /* FIXME: do we want to make this vector the full length of the
       dataset, or restrict it to the current sample length and set t1,
       t2 on it?
    */
    m = gretl_matrix_alloc(len, 1);
    memcpy(m.val, x, len * (*x).sizeof);

    return real_bundle_set_data(bundle, key, m,
                GretlType.series,
                0, null);
}

/**
 * gretl_bundle_donate_series:
 * @bundle: target bundle.
 * @key: name of key to create or replace.
 * @x: array of doubles.
 * @dset: the current dataset.
 *
 * Sets @x as a member of @bundle under the name @key:
 * the bundle takes ownership of the array.
 * If @key is already present in the bundle the original
 * value is replaced and destroyed.
 *
 * Returns: 0 on success, error code on error.
 */

int gretl_bundle_donate_series(GretlBundle * bundle, const(char) * key, double * x, const(DatasetRaw) * dset)
{
    GretlMatrixRaw * m = null;

    if (x == null) {
    gretl_errmsg_set("bundled series: no data supplied");
    return GretlError.data;
    }

    m = cast(typeof(m)) malloc((*m).sizeof);
    gretl_matrix_init_full(m, dset.n, 1, x);

    return real_bundle_set_data(bundle, key, m,
                GretlType.series,
                0, null);
}

/**
 * gretl_bundle_set_list:
 * @bundle: target bundle.
 * @key: name of key to create or replace.
 * @list: gretl list.
 *
 * Sets @list as a member of @bundle under the name @key.
 * If @key is already present in the bundle the original
 * value is replaced and destroyed.
 *
 * Returns: 0 on success, error code on error.
 */

int gretl_bundle_set_list(GretlBundle * bundle, const(char) * key, const(int) * list)
{
    return gretl_bundle_set_data(bundle, key, cast(void *) list,
                                 GretlType.list);
}

/**
 * gretl_bundle_set_matrix:
 * @bundle: target bundle.
 * @key: name of key to create or replace.
 * @m: gretl matrix.
 *
 * Sets @m as a member of @bundle under the name @key.
 * If @key is already present in the bundle the original
 * value is replaced and destroyed.
 *
 * Returns: 0 on success, error code on error.
 */

int gretl_bundle_set_matrix(GretlBundle * bundle, const(char) * key, const(GretlMatrixRaw) * m)
{
    return gretl_bundle_set_data(bundle, key, cast(void *) m,
                                 GretlType.matrix);
}

/**
 * gretl_bundle_delete_data:
 * @bundle: target bundle.
 * @key: name of key to delete.
 *
 * Deletes the data item under @key from @bundle, if
 * such an item is present.
 */

int gretl_bundle_delete_data(GretlBundle * bundle, const(char) * key)
{
    int done = 0;
    int err = 0;

    if (bundle == null) {
        return GretlError.data;
    }

    if (bundle.type == BundleType.kalman) {
        done = maybe_delete_kalman_element(bundle.data, key, &err);
    if (err) {
        return err;
    }
    }

    if (!done) {
    gpointer p = g_hash_table_lookup(bundle.ht, key);

    if (p == null) {
        gretl_errmsg_sprintf("\"%s\": no such item".ptr, key);
        err = GretlError.unkvar;
    } else {
        done = g_hash_table_remove(bundle.ht, key);
        if (!done) {
        err = GretlError.invarg;
        }
        }
    }

    return err;
}

/**
 * gretl_bundle_rekey_data:
 * @bundle: target bundle.
 * @oldkey: name of key to access.
 * @newkey: revised key.
 *
 * If @bundle contains an item under the key @oldkey, changes
 * the string by which the item is identified to @newkey.
 *
 * Returns: 0 on success, error code on failure.
 */

int gretl_bundle_rekey_data(GretlBundle * bundle, const(char) * oldkey, const(char) * newkey)
{
    BundledItem * item = null;
    int err = 0;

    if (strcmp(oldkey, newkey) == 0) {
        return 0;
    }

    if (bundle != null) {
        item = cast(BundledItem *) g_hash_table_lookup(bundle.ht, oldkey);
    }

    if (item == null) {
        err = GretlError.data;
    } else {
        g_hash_table_steal(bundle.ht, oldkey);
        g_free(item.key);
        item.key = g_strdup(newkey);
        g_hash_table_insert(bundle.ht, item.key, item);
    }

    return err;
}

/**
 * gretl_bundle_set_note:
 * @bundle: target bundle.
 * @key: name of key to access.
 * @note: note to add.
 *
 * Adds a descriptive note to the item under @key in @bundle.
 * If a note is already present it is replaced by the new one.
 *
 * Returns: 0 on success, error code on error.
 */

int gretl_bundle_set_note(GretlBundle * bundle, const(char) * key, const(char) * note)
{
    int err = 0;

    if (bundle == null) {
        err = GretlError.unkvar;
    } else {
        gpointer p = g_hash_table_lookup(bundle.ht, key);

        if (p == null) {
            err = GretlError.data;
        } else {
            BundledItem * item = cast(BundledItem *) p;

            free(item.note);
            item.note = gretl_strdup(note);
        }
    }

    return err;
}

/**
 * gretl_bundle_set_creator:
 * @bundle: target bundle.
 * @name: name of function package that built the bundle.
 *
 * Sets the "creator" attribute of @bundle. This is called
 * automatically when a bundle is returned to top-level
 * userspace by a packaged function.
 *
 * Returns: 0 on success, error code on error.
 */

int gretl_bundle_set_creator(GretlBundle * bundle, const(char) * name)
{
    int err = 0;

    if (bundle == null) {
        err = GretlError.data;
    } else {
        free(bundle.creator);
        if (name == null) {
            bundle.creator = null;
        } else {
            bundle.creator = gretl_strdup(name);
            if (bundle.creator == null) {
                err = GretlError.alloc;
            }
        }
    }

    return err;
}

/* replicate on a target bundle a BundledItem from some other
   other bundle, provided the target bundle does not already
   have a BundledItem under the same key
*/

void copy_new_bundled_item(gpointer key, gpointer value, gpointer p)
{
    BundledItem * item = cast(BundledItem *) value;
    GretlBundle * targ = cast(GretlBundle *) p;

    if (!gretl_bundle_has_data(targ, cast(const(char) *) key)) {
        real_bundle_set_data(targ, cast(const(char) *) key,
                             item.data, item.type,
                             1, item.note);
    }
}

/* replicate on a target bundle a BundledItem from some other
   bundle */

void copy_bundled_item(gpointer key, gpointer value, gpointer p)
{
    BundledItem * item = cast(BundledItem *) value;
    GretlBundle * targ = cast(GretlBundle *) p;

    real_bundle_set_data(targ, cast(const(char) *) key,
                         item.data, item.type,
                         1, item.note);
}

/* Create a new bundle as the union of two existing bundles:
   we first copy bundle1 in its entirety, then append any elements
   of bundle2 whose keys that are not already present in the
   copy-target. In case bundle1 and bundle2 share any keys, the
   value copied to the target is therefore that from bundle1.
*/

GretlBundle * gretl_bundle_union(GretlBundle * bundle1, GretlBundle * bundle2, int * err)
{
    GretlBundle * b = null;

    if (bundle2.type == BundleType.kalman) {
        gretl_errmsg_set(_("bundle union: the right-hand operand cannot be a kalman bundle").ptr);
        *err = GretlError.data;
    } else {
        b = gretl_bundle_copy(bundle1, err);
    }

    if (!*err) {
        g_hash_table_foreach(bundle2.ht, &copy_new_bundled_item, b);
    }

    return b;
}

/* argument-bundle checking apparatus */

struct Bchecker {
    GretlBundle * b;
    char ** bounds_keys;
    GretlMatrixRaw * bounds_mat;
    int n_bounds;
    char ** ignore;
    int n_ignores;
    int * ret;
    int * err;
    Prn * prn;
}

void bchecker_free(Bchecker * bc)
{
    if (bc.n_bounds > 0) {
        gretl_matrix_free(bc.bounds_mat);
        strings_array_free(bc.bounds_keys, bc.n_bounds);
    }
    free(bc);
}

Bchecker * bchecker_new(GretlBundle * b, int * ret, int * err, GretlArray * ignore, Prn * prn)
{
    Bchecker * bc = cast(Bchecker *) calloc(1, Bchecker.sizeof);

    if (bc == null) {
        *err = GretlError.alloc;
        return null;
    }

    bc.b = b;
    bc.ret = ret;
    bc.err = err;
    bc.prn = prn;

    if (gretl_bundle_has_key(b, "bounds")) {
        GretlBundle * bb;

        bb = gretl_bundle_get_bundle(b, "bounds", err);
        if (!*err) {
            *err = unpack_bounds_info(bc, bb);
        }
    }

    if (!*err && ignore != null) {
        char ** S = null;
        int ns = 0;

        /* note: borrowing here */
        S = gretl_array_get_strings(ignore, &ns);
        if (S != null) {
            bc.ignore = S;
            bc.n_ignores = ns;
        }
    }

    return bc;
}

/* Check for mismatched types in @targ vs @src, allowing for
   convertibility between scalar types.
*/

int bundled_types_mismatch(BundledItem * targ, BundledItem * src)
{
    int err = 0;

    if (src.type == targ.type) {
        return 0; /* OK */
    } else if (gretl_is_scalar_type(targ.type) &&
               gretl_is_scalar_type(src.type)) {
        /* may be inter-convertible: the main potential
           problem, diagnosed here, lies in conversion
           from double to integer types
        */
        if (src.type == GretlType.double_) {
            double * px = cast(double *) src.data;

            if (targ.type == GretlType.bool_) {
                err = na(*px);
            } else if (targ.type == GretlType.int_) {
                gretl_int_from_double(*px, &err);
            } else if (targ.type == GretlType.uint32) {
                gretl_uint32_from_double(*px, &err);
            } else if (targ.type == GretlType.uint64) {
                gretl_uint64_from_double(*px, &err);
            }
        }
    } else {
        /* other mismatches will not work */
        err = 1;
    }

    return err;
}

int should_ignore(const(char) * key, Bchecker * bc)
{
    int i;

    for (i=0; i<bc.n_ignores; i++) {
        if (!strcmp(key, bc.ignore[i])) {
            return 1;
        }
    }

    return 0;
}

double bundled_item_get_scalar(BundledItem * item)
{
    if (item.type == GretlType.double_) {
        return *cast(double *) item.data;
    } else if (item.type == GretlType.int_) {
        return cast(double) *cast(int *) item.data;
    } else {
        return NADBL;
    }
}

/* check @src, a presumably scalar bundled input value,
   against any relevant bounds specification in @bc
*/

int bcheck_bounds(const(char) * key, BundledItem * src, Bchecker * bc)
{
    double lo, hi, val;
    const(char) * s;
    int i;

    for (i=0; i<bc.n_bounds; i++) {
        s = bc.bounds_keys[i];
        if (!strcmp(s, src.key)) {
            /* found an associated limits item */
            if (!gretl_is_scalar_type(src.type)) {
                pprintf(bc.prn, _("%s: should be scalar\n").ptr, src.key);
                return 1;
            }
            val = bundled_item_get_scalar(src);
            if (na(val)) {
                pprintf(bc.prn, _("%s: invalid (missing) value\n").ptr, src.key);
                return 1;
            }
            lo = gretl_matrix_get(bc.bounds_mat, i, 0);
            hi = gretl_matrix_get(bc.bounds_mat, i, 1);
            if (!na(lo) && val < lo) {
                pprintf(bc.prn, _("%s: value %g is out of bounds (lower limit %g)\n").ptr,
                        src.key, val, lo);
                return 1;
            } else if (!na(hi) && val > hi) {
                pprintf(bc.prn, _("%s: value %g is out of bounds (upper limit %g)\n").ptr,
                        src.key, val, hi);
                return 1;
            }
            break;
        }
    }

    return 0;
}

/* Given a sub-bundle named 'bounds', unpack it to form an
   array of n >= 1 key-strings and an associated n x 2 matrix
   holding lower and upper bounds.
*/

int unpack_bounds_info(Bchecker * bc, GretlBundle * bounds)
{
    GretlMatrixRaw * bi, bb = null;
    char ** S = null;
    int i, ns = 0;
    int err = 0;

    S = gretl_bundle_get_keys_raw(bounds, &ns);
    if (S == null) {
        err = 1;
    } else {
        bb = gretl_matrix_alloc(ns, 2);
        if (bb == null) {
            err = 1;
        }
    }

    for (i=0; i<ns && !err; i++) {
        bi = gretl_bundle_get_matrix(bounds, S[i], &err);
        if (!err && gretl_vector_get_length(bi) != 2) {
            err = GretlError.invarg;
        }
        if (!err) {
            gretl_matrix_set(bb, i, 0, bi.val[0]);
            gretl_matrix_set(bb, i, 1, bi.val[1]);
        }
    }

    if (!err) {
        bc.bounds_keys = S;
        bc.bounds_mat = bb;
        bc.n_bounds = ns;
    }

    return err;
}

void check_bundled_item(gpointer key, gpointer value, gpointer p)
{
    BundledItem * targ, src = cast(BundledItem *) value;
    Bchecker * bc = cast(Bchecker *) p;

    if (*bc.ret || *bc.err) {
        /* don't waste time if we already hit an error */
        return;
    }

    if (bc.ignore != null && should_ignore(cast(const(char) *) key, bc)) {
        return;
    }

    if (bc.n_bounds > 0) {
        *bc.ret = bcheck_bounds(cast(const(char) *) key, src, bc);
        if (*bc.ret) {
            return;
        }
    }

    /* look up @key (from input) in the template bundle */
    targ = cast(typeof(targ)) g_hash_table_lookup(bc.b.ht, cast(const(char) *) key);

    if (targ == null) {
        /* extraneous key in input */
        pprintf(bc.prn, "bcheck: unrecognized key '%s'\n".ptr, key);
        *bc.ret = 2;
    } else if (targ.type == GretlType.matrix &&
               src.type == GretlType.double_) {
        double x = *cast(double *) src.data;
        GretlMatrixRaw * m = gretl_matrix_from_scalar(x);

        *bc.err = bundled_item_replace_data(targ, m, GretlType.matrix, 0);
    } else if (targ.type == GretlType.double_ &&
               src.type == GretlType.matrix) {
        GretlMatrixRaw * m = cast(GretlMatrixRaw *) src.data;

        if (gretl_matrix_is_scalar(m)) {
            *bc.err = bundled_item_replace_data(targ, &m.val[0],
                                                 GretlType.double_, 0);
        } else {
            *bc.ret = 3;
        }
    } else if (bundled_types_mismatch(targ, src)) {
        *bc.ret = 3;
    } else {
        /* transcribe input value -> template */
        *bc.err = bundled_item_replace_data(targ, src.data, src.type, 1);
        if (*bc.err) {
            pprintf(bc.prn, "bcheck: failed to copy '%s'\n".ptr, key);
        }
    }

    if (*bc.ret == 3) {
        pprintf(bc.prn, "bcheck: '%s' should be %s, is %s\n".ptr, key,
                gretl_type_get_name(targ.type),
                gretl_type_get_name(src.type));
    }
}

/**
 * gretl_bundle_extract_args:
 * @defaults: bundle containing keys for all supported
 * inputs, both optional and required (if any).
 * @input: bundle supplied by caller.
 * @reqd: array of strings identifying required keys, if any
 * (or NULL).
 * @ignore: array of strings identifying keys that are not
 * associated with options but whose presence should not
 * provoke an error.
 * @prn: pointer to printing struct for display of error
 * messages.
 * @err:location to receive error code.
 *
 * This function checks @input against @defaults. It flags an
 * error (a) if a required element is missing, (b) if @input
 * contains an unrecognized key, or (c) if the type of any element
 * in @input fails to match the type of the corresponding element
 * in @defaults. If there's no error the content of @defaults
 * is updated from @input; that is, default values are replaced by
 * values selected by the caller.
 *
 * Note the the @err pointer receives a non-zero value only
 * if this function fails, which is distinct from the case
 * where it successfully diagnoses an error in @input.
 *
 * Returns: 0 if @input is deemed valid in light of @defaults,
 * non-zero otherwise.
 */

int gretl_bundle_extract_args(GretlBundle * defaults, GretlBundle * input, GretlArray * reqd, GretlArray * ignore, Prn * prn, int * err)
{
    int ret = 0;

    if (reqd != null && gretl_array_get_type(reqd) != GretlType.strings) {
        *err = GretlError.types;
        return *err;
    }

    if (input == null) {
        if (reqd == null) {
            /* OK, it's a no-op */
            return 0;
        } else {
            /* null @input is not acceptable given @reqd */
            pprintf(prn, "bcheck: required input bundle is missing\n".ptr);
            return 1;
        }
    }

    if (reqd != null) {
        int i, n = gretl_array_get_length(reqd);
        const(char) * key;

        for (i=0; i<n; i++) {
            key = cast(const(char) *) gretl_array_get_data(reqd, i);
            if (!gretl_bundle_has_key(input, key)) {
                pprintf(prn, "bcheck: required key '%s' is missing\n".ptr, key);
                ret = 1;
                break;
            }
        }
    }

    if (ret == 0) {
        Bchecker * bc = bchecker_new(defaults, &ret, err, ignore, prn);

        if (bc != null) {
            g_hash_table_foreach(input.ht, &check_bundled_item, bc);
            bchecker_free(bc);
        }
    }

    return ret;
}

int gretl_bundle_append(GretlBundle * bundle1, GretlBundle * bundle2)
{
    if (bundle2.type == BundleType.kalman) {
        gretl_errmsg_set(_("bundle union: the right-hand operand cannot be a kalman bundle").ptr);
        return GretlError.data;
    } else {
        g_hash_table_foreach(bundle2.ht, &copy_new_bundled_item, bundle1);
        return 0;
    }
}

/**
 * gretl_bundle_copy:
 * @bundle: gretl bundle to be copied.
 * @err: location to receive error code.
 *
 * Returns: a "deep copy" of @bundle (all the items in @bundle
 * are themselves copied), or NULL on failure.
 */

GretlBundle * gretl_bundle_copy(GretlBundle * bundle, int * err)
{
    GretlBundle * bcpy = null;

    if (bundle == null) {
        *err = GretlError.data;
    } else {
        if (bundle.type == BundleType.kalman) {
            bcpy = kalman_bundle_copy(bundle, err);
        } else {
            bcpy = gretl_bundle_new();
            if (bcpy == null) {
                *err = GretlError.alloc;
            }
        }
        if (!*err) {
            g_hash_table_foreach(bundle.ht, &copy_bundled_item, bcpy);
        }
    }

    return bcpy;
}

/**
 * gretl_bundle_copy_as:
 * @name: name of source bundle.
 * @copyname: name for copy.
 *
 * Look for a saved bundle specified by @name, and if found,
 * make a full copy and save it under @copyname. This is
 * called from geneval.c on completion of assignment to a
 * bundle named @copyname, where the returned value on the
 * right-hand side is a pre-existing saved bundle.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_bundle_copy_as(const(char) * name, const(char) * copyname)
{
    GretlBundle * b0, b1 = null;
    UserVar * u;
    int prev = 0;
    int err = 0;

    if (!strcmp(name, "$sysinfo")) {
        b0 = sysinfo_bundle;
    } else {
        u = get_user_var_of_type_by_name(name, GretlType.bundle);
        if (u == null) {
            return GretlError.data;
        } else {
            b0 = cast(GretlBundle *) user_var_get_value(u);
        }
    }

    u = get_user_var_of_type_by_name(copyname, GretlType.bundle);

    if (u != null) {
        b1 = cast(GretlBundle *) user_var_steal_value(u);
        if (b1 != null) {
            gretl_bundle_destroy(b1);
        }
        prev = 1;
    }

    b1 = gretl_bundle_copy(b0, &err);

    if (!err) {
        if (prev) {
            err = user_var_replace_value(u, b1, GretlType.bundle);
        } else {
            err = user_var_add(copyname, GretlType.bundle, b1);
        }
    }

    return err;
}

struct BItemPrinter {
    Prn * prn;
    int indent;
    int tree;
}

void real_print_bundled_item(BundledItem * item, BItemPrinter * bip)
{
    GretlType t = item.type;
    const(char) * kstr = item.key;
    Prn * prn = bip.prn;
    int indent;

    if (bip.tree && t == GretlType.bundle) {
        return;
    }

    indent = 2 + 2 * bip.indent;
    bufspace(indent, prn);

    if (item.is_virtual) {
    const(char) * ts = gretl_type_get_name(item.type);
    char * s = cast(char *) item.data;

    if (strlen(s) > 15) {
        pprintf(prn, "%s (virtual %s: %s)".ptr, kstr, ts, s + 15);
    } else {
        pprintf(prn, "%s (virtual %s)".ptr, kstr, ts);
    }
    } else if (t == GretlType.double_) {
        double x = *cast(double *) item.data;

        if (na(x)) {
            pprintf(prn, "%s = NA".ptr, kstr);
        } else {
            pprintf(prn, "%s = %g".ptr, kstr, x);
        }
    } else if (t == GretlType.int_) {
        int i = *cast(int *) item.data;

        pprintf(prn, "%s = %d".ptr, kstr, i);
    } else if (t == GretlType.uint32) {
        uint u = *cast(uint *) item.data;

        pprintf(prn, "%s = %u".ptr, kstr, u);
    } else if (t == GretlType.uint64) {
        ulong u = *cast(ulong *) item.data;

        pprintf(prn, "%s = %lu".ptr, kstr, u);
    } else if (t == GretlType.string_) {
        char * s = cast(char *) item.data;
        int n = cast(int) (cast(int) strlen(s));

        if (n + indent < 70) {
            if (strchr(s, '"') != null) {
                pprintf(prn, "%s = '%s'".ptr, kstr, s);
            } else {
                pprintf(prn, "%s = \"%s\"".ptr, kstr, s);
            }
        } else {
            pprintf(prn, "%s (%s, %d bytes)".ptr, kstr,
                    gretl_type_get_name(t), n);
        }
    } else if (t == GretlType.bundle) {
        pprintf(prn, "%s (%s)".ptr, kstr, gretl_type_get_name(t));
    } else if (t == GretlType.matrix) {
        GretlMatrixRaw * m = cast(GretlMatrixRaw *) item.data;

        if (m.rows == 1 && m.cols == 1) {
            pprintf(prn, "%s = %g".ptr, kstr, m.val[0]);
        } else {
            pprintf(prn, "%s (%s: %d x %d)".ptr, kstr,
                    gretl_type_get_name(t),
                    m.rows, m.cols);
        }
    } else if (t == GretlType.series) {
    GretlMatrixRaw * m = cast(GretlMatrixRaw *) item.data;

        pprintf(prn, "%s (%s: length %d)".ptr, kstr,
                gretl_type_get_name(t), m.rows);
    } else if (t == GretlType.list) {
        pprintf(prn, "%s (%s)".ptr, kstr, gretl_type_get_name(t));
    } else if (t == GretlType.array) {
        GretlArray * a = cast(GretlArray *) item.data;
        int n = gretl_array_get_length(a);

        t = gretl_array_get_type(a);
        pprintf(prn, "%s = array of %s, length %d".ptr, kstr,
                gretl_type_get_name(t), n);
    }

    if (item.note != null) {
        pprintf(prn, " %s\n".ptr, item.note);
    } else {
        pputc(prn, '\n');
    }
}

/* wrapper for use with GList foreach */

void print_bundled_item(gpointer value, gpointer p)
{
    BundledItem * item = cast(BundledItem *) value;
    BItemPrinter * bip = cast(BItemPrinter *) p;

    real_print_bundled_item(item, bip);
}

void bundle_header(const(char) * name, const(char) * creator, int empty, Prn * prn)
{
    if (name != null) {
        if (empty) {
            pprintf(prn, "bundle %s: empty\n".ptr, name);
        } else if (creator != null) {
            pprintf(prn, "bundle %s, created by %s:\n".ptr, name, creator);
        } else {
            pprintf(prn, "bundle %s:\n".ptr, name);
        }
    } else {
        if (empty) {
            pputs(prn, "bundle: empty\n");
        } else if (creator != null) {
            pprintf(prn, "bundle created by %s:\n".ptr, creator);
        } else {
            pputs(prn, "bundle:\n");
        }
    }
}

int real_bundle_print(GretlBundle * bundle, int indent, int tree, Prn * prn)
{
    BItemPrinter bip = BItemPrinter(prn, indent, tree);
    GList *L;

    if (bundle == null) {
        return GretlError.data;
    }

    L = gretl_bundle_get_sorted_items(bundle);

    if (indent > 0) {
        /* child, when printing tree */
        int n_items = g_hash_table_size(bundle.ht);

        if (bundle.type == BundleType.plain && n_items == 0) {
            pputs(prn, "empty\n");
        } else if (bundle.type == BundleType.kalman) {
            print_kalman_bundle_info(bundle.data, prn);
            if (n_items > 0) {
                pputs(prn, "\nOther content\n");
                g_list_foreach(L, &print_bundled_item, &bip);
            }
        } else if (n_items > 0) {
            g_list_foreach(L, &print_bundled_item, &bip);
        }
    } else {
        int n_items = g_hash_table_size(bundle.ht);
        UserVar * u = get_user_var_by_data(bundle);
        const(char) * name = null;

        if (u != null) {
            name = user_var_get_name(u);
        }
        if (bundle.type == BundleType.plain && n_items == 0) {
            bundle_header(name, null, 1, prn);
        } else {
            bundle_header(name, bundle.creator, 0, prn);
            if (bundle.type == BundleType.kalman) {
                print_kalman_bundle_info(bundle.data, prn);
                if (n_items > 0) {
                    pputs(prn, "\nOther content\n");
                    g_list_foreach(L, &print_bundled_item, &bip);
                }
            } else if (n_items > 0) {
                g_list_foreach(L, &print_bundled_item, &bip);
            }
        }
    }

    g_list_free(L);

    return 0;
}

int is_bundles_array(BundledItem * item)
{
    if (item.type == GretlType.array) {
        return gretl_array_get_content_type(cast(GretlArray *) item.data) == GretlType.bundle;
    }

    return 0;
}

int print_bundle_tree(GretlBundle * b, int level, Prn * prn)
{
    GretlArray * K = gretl_bundle_get_keys(b, null);
    BItemPrinter bip = BItemPrinter(prn, level, 1);
    BundledItem * item;
    const(char) * key;
    GretlArray * A;
    GretlBundle * bj;
    char ** keys;
    int nb, nj, nmemb;
    int i, j, n = 0;

    keys = gretl_array_get_strings(K, &n);
    if (n == 0) {
        bufspace(2*(level+1), prn);
        pputs(prn, "empty\n");
    }

    for (i=0; i<n; i++) {
        key = keys[i];
        item = cast(BundledItem *) g_hash_table_lookup(b.ht, key);
        if (item.type == GretlType.bundle) {
            nmemb = gretl_bundle_get_n_members(cast(GretlBundle *) item.data);
            bufspace(2*(level+1), prn);
            if (nmemb == 0) {
                pprintf(prn, "%s = bundle (empty)\n".ptr, key);
            } else {
                pprintf(prn, "%s = bundle (%d members)\n".ptr, key, nmemb);
                print_bundle_tree(cast(GretlBundle *) item.data, level + 1, prn);
            }
        } else if (is_bundles_array(item)) {
            A = cast(GretlArray *) item.data;
            nb = gretl_array_get_length(A);
            bufspace(2*(level+1), prn);
            pprintf(prn, "%s = array of bundles (%d elements)\n".ptr, key, nb);
            for (j=0; j<nb; j++) {
                bj = gretl_array_get_bundle(A, j);
                nj = gretl_bundle_get_n_members(bj);
                bufspace(2*(level+2), prn);
                pprintf(prn, "%s[%d] (%d members):\n".ptr, key, j + 1, nj);
                print_bundle_tree(bj, level + 2, prn);
            }
        } else {
            real_print_bundled_item(item, &bip);
        }
    }

    gretl_array_destroy(K);

    return 0;
}

/**
 * gretl_bundle_print:
 * @bundle: gretl bundle.
 * @prn: gretl printer.
 *
 * Prints to @prn a list of the keys defined in @bundle, along
 * with descriptive notes, if any.
 *
 * Returns: 0 on success, non-zero code on failure.
 */

int gretl_bundle_print(GretlBundle * bundle, Prn * prn)
{
    int err = real_bundle_print(bundle, 0, 0, prn);

    if (!err) {
        pputc(prn, '\n');
    }

    return err;
}

/* helper function for gretl_bundle_write_constructor() */

void write_item_constructor(gpointer value, gpointer p)
{
    BundledItem * item = cast(BundledItem *) value;
    const(char) * kstr = item.key;
    GretlType t = item.type;
    GString * gs = cast(GString *) p;

    if (item.is_virtual) {
        char * s = cast(char *) item.data;

    g_string_append_printf(gs, "%s=\"%s\",".ptr, kstr, s);
    } else if (t == GretlType.double_) {
        double x = *cast(double *) item.data;

        if (na(x)) {
            g_string_append_printf(gs, "%s=NA,".ptr, kstr);
        } else {
            g_string_append_printf(gs, "%s=%g,".ptr, kstr, x);
        }
    } else if (t == GretlType.int_) {
        int i = *cast(int *) item.data;

        g_string_append_printf(gs, "%s=%d,".ptr, kstr, i);
    } else if (t == GretlType.uint32) {
        uint u = *cast(uint *) item.data;

        g_string_append_printf(gs, "%s=%u,".ptr, kstr, u);
    } else if (t == GretlType.uint64) {
        ulong u = *cast(ulong *) item.data;

        g_string_append_printf(gs, "%s=%lu".ptr, kstr, u);
    } else if (t == GretlType.string_) {
        char * s = cast(char *) item.data;

    g_string_append_printf(gs, "%s=\"%s\",".ptr, kstr, s);
    } else if (t == GretlType.matrix) {
    GretlMatrixRaw * m = cast(GretlMatrixRaw *) item.data;
    char * ms = gretl_matrix_write_constructor(m);

    g_string_append_printf(gs, "%s=%s,".ptr, kstr, ms);
    g_free(ms);
    } else {
    fprintf(stderr, "write bundle constructor: skipping unhandled member '%s'\n".ptr, kstr);
    }
}

/**
 * gretl_bundle_write_constructor:
 * @bundle: gretl bundle.
 *
 * For internal use. This function is designed to create a
 * string that can be used in a hansl command or function
 * call to create a copy of the @bundle argument. In that
 * context, only rather simple bundle elements are supported:
 * scalars, strings and small matrices. Large matrices, as well
 * as members of type bundle, array, series or list, are
 * ignored.
 *
 * Returns: a GLib-allocated string; use g_free() to destroy it.
 */

char * gretl_bundle_write_constructor(GretlBundle * bundle)
{
    int n_items = g_hash_table_size(bundle.ht);
    char * ret = null;

    if (n_items == 0) {
    return g_strdup("_()");
    } else {
    GList *L = g_hash_table_get_values(bundle.ht);
    gsize sz = 4 + 8 * n_items;
    GString * gs;

    gs = g_string_sized_new(sz);
    g_string_append(gs, "_(");
    gretl_push_c_numeric_locale();
    g_list_foreach(L, &write_item_constructor, gs);
    gretl_pop_c_numeric_locale();
    sz = cast(int) (strlen(gs.str));
    if (gs.str[sz-1] == ',') {
        gs.str[sz-1] = ')';
    } else {
        g_string_append_c(gs, ')');
    }
    g_list_free(L);
    ret = g_string_free(gs, 0);
    }

    return ret;
}

/**
 * gretl_bundle_debug_print:
 * @bundle: gretl bundle.
 * @msg: extra string to print, or NULL.
 *
 * Prints to stderr a list of the keys defined in @bundle, along
 * with descriptive notes, if any. If @msg is non-NULL it is
 * printed first.
 */

void gretl_bundle_debug_print(GretlBundle * bundle, const(char) * msg)
{
    Prn * prn = gretl_print_new(PrnType.stderr, null);

    if (msg != null) {
        pputs(prn, msg);
        pputc(prn, '\n');
    }
    real_bundle_print(bundle, 0, 0, prn);
    pputc(prn, '\n');
    gretl_print_destroy(prn);
}

int gretl_bundle_print_tree(GretlBundle * bundle, Prn * prn)
{
    int err;

    if (bundle == null) {
        err = GretlError.data;
    } else {
        UserVar * u = get_user_var_by_data(bundle);

        if (u != null) {
            pprintf(prn, "bundle %s\n".ptr, user_var_get_name(u));
        } else {
            pputs(prn, "anonymous bundle\n");
        }
        err = print_bundle_tree(bundle, 0, prn);
        if (!err) {
            pputc(prn, '\n');
        }
    }

    return err;
}

/* Called from gretl_func.c on return, to remove
   a given bundle from the stack of named bundles in
   preparation for handing it over to the caller,
   who will take ownership of it.
*/

GretlBundle * gretl_bundle_pull_from_stack(const(char) * name, int * err)
{
    GretlBundle * b = null;
    UserVar * u;

    u = get_user_var_of_type_by_name(name, GretlType.bundle);

    if (u != null) {
        b = cast(GretlBundle *) user_var_unstack_value(u);
    }

    if (b == null) {
        *err = GretlError.data;
    }

    return b;
}

/* serialize a gretl bundled item as XML */

int xml_put_bundled_item(gpointer keyp, gpointer value, gpointer p)
{
    const(char) * key = cast(const(char) *) keyp;
    BundledItem * item = cast(BundledItem *) value;
    size_t length;
    Prn * prn = cast(Prn *) p;
    int err = 0;

    if (item.type == GretlType.string_) {
        if (item.data == null) {
            fprintf(stderr, "bundle -> XML: skipping NULL string %s\n".ptr, key);
            return 0;
        }
    }

    pprintf(prn, "<bundled-item key=\"%s\" type=\"%s\"".ptr, key,
            gretl_type_get_name(item.type));

    if (item.is_virtual) {
    pprintf(prn, " is_virtual=\"1\"".ptr);
    }
    if (item.note != null) {
        pprintf(prn, " note=\"%s\"".ptr, item.note);
    }

    if (item.type == GretlType.string_ ||
    gretl_is_scalar_type(item.type) ||
    item.is_virtual) {
        pputc(prn, '>');
    } else {
        pputs(prn, ">\n");
    }

    if (item.type == GretlType.string_ || item.is_virtual) {
    length = cast(int) (strlen(cast(char *) item.data));
    if (length > 1000000000) {
        gretl_errmsg_sprintf("string %s: exceeds maximum string length".ptr, key);
        err = GretlError.data;
    } else {
        gretl_xml_put_string(cast(char *) item.data, prn);
    }
    } else if (item.type == GretlType.double_) {
        double x = *cast(double *) item.data;

        if (na(x)) {
            pputs(prn, "NA");
        } else {
            pprintf(prn, "%.16g".ptr, x);
        }
    } else if (item.type == GretlType.int_) {
        int i = *cast(int *) item.data;

        pprintf(prn, "%d".ptr, i);
    } else if (item.type == GretlType.uint32) {
        uint u = *cast(uint *) item.data;

        pprintf(prn, "%u".ptr, u);
    } else if (item.type == GretlType.uint64) {
        ulong u = *cast(ulong *) item.data;

        pprintf(prn, "%lu".ptr, u);
    } else if (item.type == GretlType.matrix ||
           item.type == GretlType.series) {
        GretlMatrixRaw * m = cast(GretlMatrixRaw *) item.data;

    err = gretl_matrix_serialize(m, null, prn);
    if (err) {
        gretl_errmsg_sprintf("matrix %s: exceeds maximum XML size".ptr, key);
    }
    } else if (item.type == GretlType.bundle) {
        GretlBundle * b = cast(GretlBundle *) item.data;

        err = gretl_bundle_serialize(b, null, prn);
    } else if (item.type == GretlType.array) {
        GretlArray * a = cast(GretlArray *) item.data;

        gretl_array_serialize(a, prn);
    } else if (item.type == GretlType.list) {
        int * list = cast(int *) item.data;

        gretl_list_serialize(list, null, prn);
    } else {
        fprintf(stderr, "bundle -> XML: skipping %s\n".ptr, key);
    }

    pputs(prn, "</bundled-item>\n");

    return err;
}

int gretl_bundle_serialize(GretlBundle * b, const(char) * name, Prn * prn)
{
    int err = 0;

    pputs(prn, "<gretl-bundle");
    if (name != null) {
        pprintf(prn, " name=\"%s\"".ptr, name);
    }
    if (b.creator != null && *b.creator != '\0') {
        pprintf(prn, " creator=\"%s\"".ptr, b.creator);
    }
    if (b.type == BundleType.kalman) {
        pputs(prn, " type=\"kalman\"");
    }
    pputs(prn, ">\n");

    if (b.type == BundleType.kalman) {
        kalman_serialize(b.data, prn);
    }

    if (b.ht != null) {
    GHashTableIter iter;
    gpointer keyp, value;

    g_hash_table_iter_init(&iter, b.ht);
    while (g_hash_table_iter_next(&iter, &keyp, &value) && !err) {
        err = xml_put_bundled_item(keyp, value, prn);
    }
    }

    pputs(prn, "</gretl-bundle>\n");

    return err;
}

int load_bundled_items(GretlBundle * b, xmlNodePtr cur, xmlDocPtr doc)
{
    GretlType type;
    int is_virtual;
    int size;
    char * key;
    int err = 0;

    while (cur != null && !err) {
        if (!xmlStrcmp(cur.name, cast(XUC) "bundled-item")) {
        size = 0;
            key = cast(char *) xmlGetProp(cur, cast(XUC) "key");
            type = gretl_xml_get_type_property(cur);
        is_virtual = gretl_xml_get_prop_as_bool(cur, "is_virtual");
        gretl_xml_get_prop_as_int(cur, "size", &size); /* legacy */
            if (key == null || type == 0) {
                err = GretlError.data;
            } else {
        if (type == GretlType.string_ || is_virtual) {
                    char * s;

                    if (!gretl_xml_node_get_trimmed_string(cur, doc, &s)) {
                        err = GretlError.data;
                    } else {
                        err = gretl_bundle_donate_data(b, key, s, GretlType.string_);
                    }
                } else if (type == GretlType.double_) {
                    double x;

                    if (!gretl_xml_node_get_double(cur, doc, &x)) {
                        err = GretlError.data;
                    } else {
                        err = gretl_bundle_set_data(b, key, &x, type);
                    }
                } else if (type == GretlType.int_) {
                    int i;

                    if (!gretl_xml_node_get_int(cur, doc, &i)) {
                        err = GretlError.data;
                    } else {
                        err = gretl_bundle_set_data(b, key, &i, type);
                    }
                } else if (type == GretlType.uint32) {
                    uint u;

                    if (!gretl_xml_node_get_uint32(cur, doc, &u)) {
                        err = GretlError.data;
                    } else {
                        err = gretl_bundle_set_data(b, key, &u, type);
                    }
                } else if (type == GretlType.uint64) {
                    ulong u;

                    if (!gretl_xml_node_get_uint64(cur, doc, &u)) {
                        err = GretlError.data;
                    } else {
                        err = gretl_bundle_set_data(b, key, &u, type);
                    }
                } else if (type == GretlType.series && size > 0) {
            /* legacy representation of bundled series */
                    double * xvec = gretl_xml_get_double_array(cur, doc, &size, &err);

                    if (!err) {
            DatasetRaw dset;

            dset.n = size;
            dset.t1 = 0;
            dset.t2 = dset.n - 1;
            gretl_bundle_donate_series(b, key, xvec, &dset);
                    }
                } else if (type == GretlType.matrix ||
               type == GretlType.series) {
                    xmlNodePtr child = cur.children;
                    GretlMatrixRaw * m;

                    if (child == null) {
                        err = GretlError.data;
                    } else {
                        m = gretl_xml_get_matrix(child, doc, &err);
                        if (!err) {
                            err = gretl_bundle_donate_data(b, key, m, type);
                        }
                    }
                } else if (type == GretlType.bundle) {
                    xmlNodePtr child = cur.children;
                    GretlBundle * baby;

                    if (child == null) {
                        err = GretlError.data;
                    } else {
                        baby = gretl_bundle_deserialize(child, doc, &err);
                        if (!err) {
                            err = gretl_bundle_donate_data(b, key, baby, type);
                        }
                    }
                } else if (type == GretlType.array) {
                    xmlNodePtr child = cur.children;
                    GretlArray * a;

                    if (child == null) {
                        err = GretlError.data;
                    } else {
                        a = gretl_array_deserialize(child, doc, &err);
                        if (!err) {
                            err = gretl_bundle_donate_data(b, key, a, type);
                        }
                    }
                } else if (type == GretlType.list) {
                    xmlNodePtr child = cur.children;
                    int * list;

                    if (child == null) {
                        err = GretlError.data;
                    } else {
                        list = gretl_xml_get_list(child, doc, &err);
                        if (!err) {
                            err = gretl_bundle_donate_data(b, key, list, type);
                        }
                    }
                } else {
                    fprintf(stderr, "bundle: ignoring unhandled type %d\n".ptr, type);
                }

                if (!err) {
                    char * note = cast(char *) xmlGetProp(cur, cast(XUC) "note");

                    if (note != null) {
                        gretl_bundle_set_note(b, key, note);
                        xmlFree(note);
                    }
                }

                xmlFree(key);
            }
        }
        cur = cur.next;
    }

    return err;
}

/* For internal use only: @p1 should be of type xmlNodePtr and @p2
   should be an xmlDocPtr. We suppress the actual pointer types in the
   signature so that it's possible for a module to include
   gretl_bundle.h without including the full libxml headers.
*/

GretlBundle * gretl_bundle_deserialize(void * p1, void * p2, int * err)
{
    xmlNodePtr cur, node = cast(xmlNode *) p1;
    xmlDocPtr doc = cast(xmlDoc *) p2;
    GretlBundle * b = null;
    char * btype = null;
    char * creator = null;

    btype = cast(char *) xmlGetProp(node, cast(XUC) "type");
    creator = cast(char *) xmlGetProp(node, cast(XUC) "creator");

    cur = node.children;

    if (btype != null && !strcmp(btype, "kalman")) {
        b = kalman_deserialize(cur, doc, err);
    } else {
        b = gretl_bundle_new();
        if (b == null) {
            *err = GretlError.alloc;
        } else if (creator != null) {
            b.creator = creator;
            creator = null;
        }
    }

    free(btype);
    free(creator);

    if (b != null) {
        *err = load_bundled_items(b, cur, doc);
        if (*err) {
            fprintf(stderr, "deserialize bundle: bundle is broken (err = %d)\n".ptr, *err);
            gretl_bundle_destroy(b);
            b = null;
        }
    }

    return b;
}

int call_bundle_to_json(GretlBundle * b, const(char) * fname, int control, char ** buf)
{
    alias JfuncFn = extern (C) int function(GretlBundle *, const(char) *, Gretlopt, char **);
    JfuncFn jfunc;

    jfunc = cast(JfuncFn) get_plugin_function("bundle_to_json");
    if (jfunc == null) {
        return GretlError.fopen;
    } else {
        Gretlopt opt = Gretlopt.OPT_NONE;

        if (control & 2) {
            opt |= Gretlopt.OPT_A;
        }
        if (control & 4) {
            opt |= Gretlopt.OPT_P;
        }
        return jfunc(b, fname, opt, buf);
    }
}

int gretl_bundle_write_to_file(GretlBundle * b, const(char) * fname, int control)
{
    char[FILENAME_MAX] fullname;
    Prn * prn = null;
    int err = 0;

    if (2 * strlen(fname) >= fullname.sizeof) {
        gretl_errmsg_sprintf("Filename too long, %d bytes".ptr,
                             cast(int) strlen(fname));
        return GretlError.data;
    }

    if (control & 1) {
        gretl_build_path(fullname.ptr, gretl_dotdir(), fname, null);
    } else {
        strcpy(fullname.ptr, fname);
        gretl_maybe_switch_dir(fullname.ptr);
    }

    if (has_suffix(fname, ".json") || has_suffix(fname, ".geojson")) {
        return call_bundle_to_json(b, fullname.ptr, control, null);
    }

    if (has_suffix(fname, ".gz")) {
        prn = gretl_gzip_print_new(fullname.ptr, -1, &err);
    } else {
        prn = gretl_print_new_with_filename(fullname.ptr, &err);
    }

    if (prn != null) {
        gretl_push_c_numeric_locale();
        gretl_xml_header(prn);
        err = gretl_bundle_serialize(b, null, prn);
        gretl_print_destroy(prn);
        gretl_pop_c_numeric_locale();
    if (err) {
        gretl_remove(fullname.ptr);
    }
    }

    return err;
}

char * gretl_bundle_write_to_json_buffer(GretlBundle * b, int * err)
{
    char * buf = null;

    *err = call_bundle_to_json(b, null, 0, &buf);
    return buf;
}

GretlBundle * read_json_bundle(const(char) * fname, int * err)
{
    alias JfuncFn = extern (C) GretlBundle * function(const(char) *, const(char) *, GretlArray *, int *);
    JfuncFn jfunc;
    GretlBundle * b = null;
    GError *gerr = null;
    char * JSON = null;
    gsize len = 0;
    gboolean ok;

    ok = g_file_get_contents(fname, &JSON, &len, &gerr);

    if (ok) {
        jfunc = cast(JfuncFn) get_plugin_function("json_get_bundle");
        if (jfunc == null) {
            *err = GretlError.fopen;
        } else {
            b = jfunc(JSON, null, null, err);
        }
        g_free(JSON);
    } else if (gerr != null) {
        *err = GretlError.fopen;
        gretl_errmsg_set(gerr.message);
        g_error_free(gerr);
    }

    return b;
}

GretlBundle * read_shapefile_bundle(const(char) * fname, int * err)
{
    alias BfuncFn = extern (C) GretlBundle * function(const(char) *, int *);
    BfuncFn bfunc;
    GretlBundle * b = null;

    bfunc = cast(BfuncFn) get_plugin_function("shp_get_bundle");
    if (bfunc == null) {
        *err = GretlError.fopen;
    } else {
        b = bfunc(fname, err);
    }

    return b;
}

GretlBundle * gretl_bundle_read_from_file(const(char) * fname, int from_dotdir, int * err)
{
    char[FILENAME_MAX] fullname;
    xmlDocPtr doc = null;
    xmlNodePtr cur = null;
    GretlBundle * b = null;

    if (2 * strlen(fname) >= fullname.sizeof) {
    gretl_errmsg_sprintf("Filename too long, %d bytes".ptr,
                 cast(int) strlen(fname));
    *err = GretlError.data;
    return null;
    }

    if (from_dotdir) {
        gretl_build_path(fullname.ptr, gretl_dotdir(), fname, null);
    } else {
        strcpy(fullname.ptr, fname);
        gretl_maybe_prepend_dir(fullname.ptr);
    }

    if (has_suffix(fname, ".json") || has_suffix(fname, ".geojson")) {
        b = read_json_bundle(fullname.ptr, err);
    } else if (has_suffix(fname, ".shp")) {
        b = read_shapefile_bundle(fullname.ptr, err);
    } else {
        *err = gretl_xml_open_doc_root(fullname.ptr, "gretl-bundle", &doc, &cur);
        if (!*err) {
            gretl_push_c_numeric_locale();
            b = gretl_bundle_deserialize(cur, doc, err);
            gretl_pop_c_numeric_locale();
            xmlFreeDoc(doc);
        }
    }

    if (*err && b != null) {
        gretl_bundle_destroy(b);
        b = null;
    }

    return b;
}

GretlBundle * gretl_bundle_read_from_buffer(const(char) * buf, int len, int * err)
{
    xmlDocPtr doc = null;
    GretlBundle * b = null;

    doc = xmlReadMemory(buf, len, null, null,
                        RetrievalIndex.huge | XML_PARSE_NOBLANKS);

    if (doc == null) {
        gretl_errmsg_set(_("xmlReadMemory failed").ptr);
        *err = 1;
    } else {
        xmlNodePtr cur = xmlDocGetRootElement(doc);

        if (cur == null) {
            gretl_errmsg_set(_("xmlDocGetRootElement failed").ptr);
            *err = 1;
        } else {
            gretl_push_c_numeric_locale();
            b = gretl_bundle_deserialize(cur, doc, err);
            gretl_pop_c_numeric_locale();
        }
        xmlFreeDoc(doc);
    }

    if (*err && b != null) {
        gretl_bundle_destroy(b);
        b = null;
    }

    return b;
}

/* get the key strings from @b in the form of a gretl_array */

GretlArray * gretl_bundle_get_keys(GretlBundle * b, int * err)
{
    GretlArray * A = null;
    int myerr = 0;

    if (b == null || b.ht == null) {
        myerr = GretlError.data;
    } else {
        GList *keys = g_hash_table_get_keys(b.ht);
        guint n;

        if (keys != null && (n = g_list_length(keys)) > 0) {
            A = gretl_array_new(GretlType.strings, n, &myerr);
            if (!myerr) {
                GList *L = g_list_first(keys);
                int i = 0;

                while (L != null) {
                    gretl_array_set_string(A, i++, cast(char *) L.data, 1);
                    L = L.next;
                }
            }
        } else {
            A = gretl_array_new(GretlType.strings, 0, &myerr);
        }
        if (keys != null) {
            g_list_free(keys);
        }
    }

    if (err != null) {
        *err = myerr;
    }

    return A;
}

/* Get the key strings from @b in the form of a "raw" array
   of C type char *, terminated by a NULL pointer.
*/

char ** gretl_bundle_get_keys_raw(GretlBundle * b, int * ns)
{
    char ** S = null;

    *ns = 0;

    if (b != null && b.ht != null) {
        GList *keys = g_hash_table_get_keys(b.ht);
        guint n;

        if (keys != null && (n = g_list_length(keys)) > 0) {
            S = strings_array_new(n + 1);
            if (S != null) {
                GList *L = g_list_first(keys);
                int i = 0;

                while (L != null) {
                    S[i++] = gretl_strdup(cast(const(char) *) L.data);
                    L = g_list_next(L);
                }
                *ns = n;
                S[i] = null; /* terminator */
            }
        }
        if (keys != null) {
            g_list_free(keys);
        }
    }

    return S;
}

/* apparatus for supporting virtual bundle members */

enum T_UNSET = -1;
enum T_INVALID = -2;

struct Vsinfo {
    int t1;
    int t2;
    int nids;
    char ** ids;
}

/* Check a bundled object for inclusion in the private namespace for
   generating a virtual series. We're looking for objects whose names
   appear in the array vi->ids of length vi->nids; these objects will be
   either matrices or scalars.
*/

void check_bundled_object(gpointer listitem, gpointer p)
{
    BundledItem * bi = cast(BundledItem *) listitem;
    Vsinfo *vi = cast(Vsinfo *) p;

    if (strings_array_position(vi.ids, vi.nids, bi.key) >= 0) {
    if (bi.type == GretlType.matrix) {
        /* found a matrix named in vi->ids */
        GretlMatrixRaw * m = cast(GretlMatrixRaw *) bi.data;

        if (gretl_matrix_is_dated(m)) {
        int t1 = gretl_matrix_get_t1(m);
        int t2 = gretl_matrix_get_t2(m);

        if (vi.t1 == T_UNSET) {
            vi.t1 = t1;
            vi.t2 = t2;
        }
        }
        private_matrix_add_as_shell(m, bi.key);
    } else if (bi.type == GretlType.double_) {
        /* found a scalar named in vi->ids */
        double val = *cast(double *) bi.data;

        private_scalar_add(val, bi.key);
    } else if (bi.type == GretlType.bundle) {
        /* found a sub-bundle in vi->ids */
        GretlBundle * b = cast(GretlBundle *) bi.data;

        private_bundle_add_as_shell(b, bi.key);
    }
    }
}

int set_bundle_namespace(GretlBundle * b, Vsinfo * vi)
{
    GList *list;

    list = g_hash_table_get_values(b.ht);
    g_list_foreach(list, &check_bundled_object, vi);
    g_list_free(list);

    return 0;
}

/* Crawl the expression @s looking for sub-strings that are the
   identifiers of gretl objects. Return an array of such identifiers.
*/

char ** ids_in_expr(const(char) * s, int * ns)
{
    char ** S = null;
    char * id;
    int n;

    while (*s) {
    n = gretl_namechar_spn(s);
    if (n > 0 && s[n] == '(') {
        s += n + 1;
    } else if (n > 0) {
        id = gretl_strndup(s, n);
        /* If @id is not already in @S and looks like a variable id,
           append it to @S via "donation". Otherwise free it.
        */
        if (strings_array_position(S, *ns, id) < 0 &&
        !(genr_function_word(id) || is_user_function(id))) {
        strings_array_donate(&S, ns, id);
        } else {
        free(id);
        }
        s += n;
    } else {
        s++;
    }
    }

    return S;
}

GretlMatrixRaw * bundle_get_virtual_object(GretlBundle * b, GretlType type, const(char) * s, DatasetRaw * dset, int * err)
{
    GretlMatrixRaw * ret = null;
    Vsinfo vi = Vsinfo(T_UNSET, T_UNSET, 0, null);
    const(char) * gen_str = get_gen_str(s);
    const(char) * creator;

    vi.ids = ids_in_expr(gen_str, &vi.nids);
    if (vi.nids == 0) {
    *err = GretlError.invarg;
    goto bailout;
    }

    set_bundle_namespace(b, &vi);
    creator = gretl_bundle_get_creator(b);
    set_bundle_pkg(creator);
    ret = generate_matrix(gen_str, dset, err);
    set_bundle_pkg(null);

    if (type == GretlType.series && *err == 0) {
    if (ret.cols > 1) {
        *err = GretlError.types;
    } else if (vi.t1 >= 0 && vi.t2 >= vi.t1 &&
           ret.rows == vi.t2 - vi.t1 + 1) {
        /* heuristic: could perhaps be wrong? */
        gretl_matrix_set_t1(ret, vi.t1);
        gretl_matrix_set_t2(ret, vi.t2);
    }
    }

    destroy_private_uvars();
    strings_array_free(vi.ids, vi.nids);

 bailout:

    if (*err) {
    gretl_errmsg_sprintf("'%s': invalid virtual %s specifier".ptr,
                 gretl_type_get_name(type), s);
    }

    return ret;
}

/* end of apparatus for supporting virtual series in bundles */

int make_sysinfo_bundle()
{
    GretlBundle * b, fb;
    char * s1, s2, s3 = null;
    int ival = 0;

    b = sysinfo_bundle = gretl_bundle_new();
    if (b == null) {
    return GretlError.alloc;
    }

    gretl_bundle_set_scalar(b, "mpi", cast(double) ival);
    ival = gretl_max_mpi_processes();
    gretl_bundle_set_scalar(b, "mpimax", cast(double) ival);
    ival = gretl_n_processors();
    gretl_bundle_set_scalar(b, "nproc", cast(double) ival);
    ival = gretl_n_physical_cores();
    gretl_bundle_set_scalar(b, "ncores", cast(double) ival);
    ival = get_stack_size();
    gretl_bundle_set_scalar(b, "stack_size", cast(double) ival);
    ival = gretl_in_gui_mode();
    gretl_bundle_set_scalar(b, "gui_mode", cast(double) ival);
    ival = 0;
    gretl_bundle_set_scalar(b, "omp", cast(double) ival);
    ival = (void *).sizeof == 8 ? 64 : 32;
    gretl_bundle_set_scalar(b, "wordlen", cast(double) ival);
    gretl_bundle_set_string(b, "os", "linux");
    gretl_bundle_set_string(b, "hostname", g_get_host_name());
    gretl_bundle_set_string(b, "blas", blas_variant_string());
    gretl_bundle_set_string(b, "gnuplot", gnuplot_version_string());
    if (get_blas_details(&s1, &s2, &s3)) {
    gretl_bundle_set_string(b, "blascore", s1);
    gretl_bundle_set_string(b, "blas_parallel", s2);
    if (s3 != null) {
        gretl_bundle_set_string(b, "blas_version", s3);
    }
    }
    /* allocated string */
    gretl_bundle_donate_data(b, "cpuid", get_cpu_details(),
                 GretlType.string_);
    ival = avx_support();
    gretl_bundle_set_scalar(b, "avx", cast(double) ival);
    /* information pertaining to 'foreign' programs */
    fb = foreign_info();
    if (fb != null) {
    gretl_bundle_donate_data(b, "foreign", fb,
                 GretlType.bundle);
    }

    return 0;
}

GretlBundle * get_sysinfo_bundle(int * err)
{
    if (sysinfo_bundle == null) {
    *err = make_sysinfo_bundle();
    }

    if (sysinfo_bundle != null) {
    /* Memory usage should be the only element that
       might need updating, and we won't stop the
       show if this fails.
    */
    GretlMatrixRaw * memvals = gretl_matrix_alloc(1, 2);

    if (memvals != null) {
        char ** S = cast(char **) malloc(2 * (char *).sizeof);

        memory_stats(memvals.val);
        S[0] = gretl_strdup("MBtotal");
        S[1] = gretl_strdup("MBfree");
        gretl_matrix_set_colnames(memvals, S);
        gretl_bundle_donate_data(sysinfo_bundle,
                     "mem", memvals,
                     GretlType.matrix);
    }
    }

    return sysinfo_bundle;
}

void * sysinfo_bundle_get_data(const(char) * key, GretlType * type, int * err)
{
    GretlBundle * b = get_sysinfo_bundle(err);
    void * ret = null;

    if (b != null) {
        ret = gretl_bundle_get_data(b, key, type, err);
    }

    return ret;
}

/* For a single-equation model, create a bundle containing all the data
   available via $-accessors.
*/

GretlBundle * bundle_from_model(Model * pmod, DatasetRaw * dset, int * err)
{
    GretlBundle * b = null;
    GretlMatrixRaw * m;
    GretlArray * a;
    double * x;
    double val;
    int * list;
    const(char) * s;
    const(char) * key;
    int i, t, berr;

    if (pmod == null) {
        GretlObjType type = cast(GretlObjType) 0;
        void * p = get_genr_model(&type);

        if (p == null || type != GretlObjType.eqn) {
            p = get_last_model(&type);
            if (p == null || type != GretlObjType.eqn) {
                gretl_errmsg_sprintf(_("%s: no data available").ptr, "$model".ptr);
                *err = GretlError.data;
                return null;
            }
        }
        pmod = cast(Model *) p;
    }

    x = cast(typeof(x)) malloc(dset.n * (*x).sizeof);
    if (x == null) {
        *err = GretlError.alloc;
        return null;
    }

    b = gretl_bundle_new();
    if (b == null) {
        free(x);
        *err = GretlError.alloc;
        return null;
    }

    for (i=ModelDataIndex.ess; i<ModelDataIndex.scalar_max && !*err; i++) {
        /* model scalars */
        berr = 0;
        if (i == ModelDataIndex.dwpval) {
            /* Durbin-Watson p-value: don't include this unless
               it has already been computed and attached to the
               model, since it may involve heavy computation.
            */
            val = gretl_model_get_double(pmod, "dwpval");
            if (!na(val)) {
                *err = gretl_bundle_set_scalar(b, "dwpval", val);
            }
        } else {
            val = gretl_model_get_scalar(pmod, cast(ModelDataIndex) i, dset, &berr);
            if (!berr) {
                key = mvarname(i) + 1;
                *err = gretl_bundle_set_scalar(b, key, val);
            }
        }
    }

    for (i=ModelDataIndex.scalar_max+1; i<ModelDataIndex.series_max && !*err; i++) {
        /* model series */
        for (t=0; t<dset.n; t++) {
            x[t] = NADBL;
        }
        berr = gretl_model_get_series(x, pmod, dset, cast(ModelDataIndex) i);
        if (!berr) {
            key = mvarname(i) + 1;
            *err = gretl_bundle_set_series(b, key, x, dset);
        }
    }

    for (i=ModelDataIndex.series_max+1; i<ModelDataIndex.matrix_max && !*err; i++) {
        /* model matrices */
        berr = 0;
        m = gretl_model_get_matrix(pmod, cast(ModelDataIndex) i, &berr);
        if (!berr) {
            key = mvarname(i) + 1;
            *err = gretl_bundle_donate_data(b, key, m,
                                            GretlType.matrix);
        }
    }

    for (i=ModelDataIndex.mbuild_max+1; i<ModelDataIndex.list_max && !*err; i++) {
        /* model lists */
        list = null;
        if (i == ModelDataIndex.xlist) {
            list = gretl_list_copy(pmod.xlist);
        } else if (i == ModelDataIndex.ylist) {
            list = gretl_model_get_y_list(pmod);
        }
        if (list != null) {
            key = mvarname(i) + 1;
            *err = gretl_bundle_donate_data(b, key, list,
                                            GretlType.list);
        }
    }

    for (i=ModelDataIndex.list_max+1; i<ModelDataIndex.parnames && !*err; i++) {
        /* model strings */
        s = null;
        if (i == ModelDataIndex.depvar) {
            s = gretl_model_get_depvar_name(pmod, dset);
        } else if (i == ModelDataIndex.command) {
            s = gretl_command_word(pmod.ci);
        }
        if (s != null && *s != '\0') {
            key = mvarname(i) + 1;
            *err = gretl_bundle_set_string(b, key, s);
        }
    }

    for (i=ModelDataIndex.parnames; i<ModelDataIndex.max && !*err; i++) {
        /* model arrays */
        a = null;
        if (i == ModelDataIndex.parnames) {
            a = gretl_model_get_param_names(pmod, dset, err);
        }
        if (a != null) {
            key = mvarname(i) + 1;
            *err = gretl_bundle_donate_data(b, key, a,
                                            GretlType.array);
        }
    }

    if (!*err) {
        /* miscellaneous data items */
        *err = bundlize_model_data_items(pmod, b);
    }

    free(x);

    /* don't return a broken bundle */
    if (*err && b != null) {
        gretl_bundle_destroy(b);
        b = null;
    }

    return b;
}

/* For an estimated system of some sort, create a bundle containing
   relevant data.
*/

GretlBundle * bundle_from_system(void * ptr, int type, DatasetRaw * dset, int * err)
{
    GretlVar * var = null;
    GretlObjType otype = cast(GretlObjType) type;
    EquationSystem * sys = null;
    GretlBundle * b = null;

    if (ptr == null) {
        ptr = get_genr_model(&otype);
        if (ptr == null || otype == GretlObjType.eqn) {
            ptr = get_last_model(&otype);
        }
    }

    if (ptr == null) {
        gretl_errmsg_sprintf(_("%s: no data available").ptr, "$system".ptr);
        *err = GretlError.data;
    } else if (otype == GretlObjType.var) {
        var = cast(GretlVar *) ptr;
    } else if (otype == GretlObjType.sys) {
        sys = cast(EquationSystem *) ptr;
    } else {
        gretl_errmsg_sprintf(_("%s: no data available").ptr, "$system".ptr);
        *err = GretlError.data;
    }

    if (!*err) {
        b = gretl_bundle_new();
        if (b == null) {
            *err = GretlError.alloc;
            return null;
        }
    }

    if (var != null) {
        *err = gretl_VAR_bundlize(var, dset, b);
    } else if (sys != null) {
        *err = equation_system_bundlize(sys, b);
    }

    /* don't return a broken bundle */
    if (*err && b != null) {
        gretl_bundle_destroy(b);
        b = null;
    }

    return b;
}

GretlBundle * kalman_bundle_new(GretlMatrixRaw ** M, int * copy, int nmat, int dkvar, int * err)
{
    GretlBundle * b = gretl_bundle_new();

    if (b == null) {
        *err = GretlError.alloc;
    } else {
        b.type = BundleType.kalman;
        b.data = kalman_new_minimal(M, copy, nmat, dkvar, err);
    }

    /* don't return a broken bundle */
    if (*err && b != null) {
        gretl_bundle_destroy(b);
        b = null;
    }

    return b;
}

gint sort_bundled_items(const(void) * a, const(void) * b)
{
    const(BundledItem) * ia = cast(const(BundledItem) *) a;
    const(BundledItem) * ib = cast(const(BundledItem) *) b;
    int ta = gretl_type_get_order(ia.type);
    int tb = gretl_type_get_order(ib.type);
    int ret = ta - tb;

    if (ret == 0) {
        ret = g_ascii_strcasecmp(ia.key, ib.key);
    }

    return ret;
}

GList * gretl_bundle_get_sorted_items(GretlBundle * b)
{
    GList *blist;

    blist = g_hash_table_get_values(b.ht);
    blist = g_list_sort(blist, &sort_bundled_items);

    return blist;
}

int gretl_bundles_are_equal(GretlBundle * b1, GretlBundle * b2)
{
    GList *k1 = g_hash_table_get_keys(b1.ht);
    void * v1, v2;
    const(char) * s1;
    GretlType t1, t2;
    int virt1, virt2;
    int err;
    int ret = 1;

    while (k1 != null && ret) {
    s1 = cast(const(char) *) k1.data;
    v1 = gretl_bundle_get_data_full(b1, s1, &t1, &virt1, &err);
    if (!gretl_bundle_has_key(b2, s1)) {
        ret = 0;
    } else {
        v2 = gretl_bundle_get_data_full(b2, s1, &t2, &virt2, &err);
        if (t2 != t1 || virt1 != virt2) {
        ret = 0;
        } else if (t1 == GretlType.string_ || virt1) {
        if (strcmp(cast(const(char) *) v1, cast(const(char) *) v2)) {
            ret = 0;
        }
        } else if (gretl_is_scalar_type(t1)) {
        double d1 = gretl_bundle_get_scalar(b1, s1, null);
        double d2 = gretl_bundle_get_scalar(b2, s1, null);
        int nas = na(d1) + na(d2);

        if (d1 != d2 && nas != 2) {
            ret = 0;
        }
        } else if (t1 == GretlType.matrix) {
        int eq = gretl_matrices_are_equal(cast(const(GretlMatrixRaw) *) v1, cast(const(GretlMatrixRaw) *) v2, 0, &err);

        if (eq != 1) {
            ret = 0;
        }
        } else if (t1 == GretlType.bundle) {
        if (!gretl_bundles_are_equal(cast(GretlBundle *) v1, cast(GretlBundle *) v2)) {
            ret = 0;
        }
        } else if (t1 == GretlType.array) {
        if (!gretl_arrays_are_equal(cast(const(GretlArray) *) v1, cast(const(GretlArray) *) v2)) {
            ret = 0;
        }
        } else {
        /* must be GRETL_TYPE_LIST */
        if (gretl_list_cmp(cast(const(int) *) v1, cast(const(int) *) v2)) {
            ret = 0;
        }
        }
    }
    k1 = k1.next;
    }

    if (ret) {
    GList *k2 = g_hash_table_get_keys(b2.ht);

    while (k2 != null && ret) {
        if (!gretl_bundle_has_key(b1, cast(const(char) *) k2.data)) {
        ret = 0;
        }
        k2 = k2.next;
    }
    g_list_free(k2);
    }

    g_list_free(k1);

    return ret;
}

void gretl_bundle_cleanup()
{
    if (sysinfo_bundle != null) {
        gretl_bundle_destroy(sysinfo_bundle);
        sysinfo_bundle = null;
    }
}

