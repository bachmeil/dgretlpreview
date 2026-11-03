/* Port of lib/src/uservar.c: the table of user-defined variables (scalars, matrices, strings, lists, bundles,
   arrays): adding, deleting, replacing and looking up by name or by data pointer, scoping by function depth (two
   hash tables), private and shell variables, scalar fast paths, copying between levels, and (de)serialisation of
   the stored variables. The type UserVar and the callback type are in uservarh.d (uservar_priv.h is folded into
   it). The header uservar.h has only prototypes. */

module gretl.uservar;


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
import gretl.ghashtable;
import gretl.glist;
import gretl.uservarh : UVFlags, UserVar, UvarAction, USER_VAR_FUNC;
import gretl.geneval1 : get_const_by_name;
import gretl.glib : gpointer;
import gretl.gretl_array : gretl_array_copy, gretl_array_destroy, gretl_array_get_type, gretl_array_new, gretl_array_print;
import gretl.gretl_bundle : get_bundle_by_name, gretl_bundle_copy, gretl_bundle_deserialize, gretl_bundle_destroy, gretl_bundle_get_scalar, gretl_bundle_new, gretl_bundle_print, gretl_bundle_print_tree, gretl_bundle_serialize, gretl_bundle_set_creator;
import gretl.gretl_cmatrix : gretl_cmatrix_print;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_errors : gretl_errmsg_sprintf;
import gretl.gretl_list : gretl_is_midas_list, gretl_list_append_list, gretl_list_copy, gretl_list_delete_at_pos, gretl_list_drop, gretl_list_new, gretl_list_print, gretl_null_list, in_gretl_list;
import gretl.gretl_paths : gretl_fopen;
import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_new_with_filename, pprintf, pputc, pputs;
import gretl.gretl_string_table : get_built_in_string_by_name;
import gretl.gretl_xml : gretl_list_serialize, gretl_matrix_serialize, gretl_xml_get_list, gretl_xml_get_matrix, gretl_xml_get_prop_as_string, gretl_xml_header, gretl_xml_open_doc_root;
import gretl.libset : gretl_messages_on;
import gretl.matrix_extra : gretl_matrix_print_to_prn;
import gretl.options : get_optval_string;
import gretl.printout : bufgets_finalize, list_series;
import gretl.strutils : gretl_strdup;
import gretl.strutilsh : SLASH;
import gretl.gretl_typemap : gretl_is_array_ref_type, gretl_type_from_string, gretl_type_get_name, gretl_type_get_plain_type;
import gretl.genfuncs : list_ok_dollar_vars;

/* unported dependencies: per-file declarations from the 2026c headers */
extern (C) {
    int const_lookup(const(char) * s);
    int current_function_size();
    int gretl_function_depth();
    int scalar_is_read_only_index(const(char) * name);
}


enum LEVEL_AUTO = -1;
enum LEV_PRIVATE = -1;

__gshared UserVar ** uvars = null;
__gshared int n_vars = 0;
__gshared int n_alloc = 0;
__gshared int scalar_imin = 0;

/* callback for the benefit of the edit scalars window
   in the gretl GUI */

__gshared void function() scalar_edit_callback = null;

/* callback for adding or deleting icons representing
   things in the GUI session window */

__gshared USER_VAR_FUNC user_var_callback = null;

enum UV_CHUNK = 32;

bool var_is_private(const(UserVar) * u)
{
    return (u.flags & UVFlags.private_) || u.name[0] == '$' || u.name[0] == '_';
}
bool var_is_shell(const(UserVar) * u)
{
    return (u.flags & UVFlags.shell) != 0;
}

double * na_ptr()
{
    double * px = cast(double *) malloc(double.sizeof);

    if (px != null) {
        *px = NADBL;
    }

    return px;
}

UserVar * user_var_new(const(char) * name, GretlType type, void * value, int * err)
{
    UserVar * u;

    if (type == GretlType.none) {
        *err = GretlError.data;
        fputs("user_var_new: type = GRETL_TYPE_NONE\n", stderr);
        return null;
    }

    u = cast(typeof(u)) malloc((*u).sizeof);

    if (u == null) {
        *err = GretlError.alloc;
    } else {
        u.type = type;
        u.level = gretl_function_depth();
        u.flags = (u.level == 0)? UVFlags.main : cast(UVFlags) 0;
        *u.name.ptr = '\0';
        strncat(u.name.ptr, name, VNAMELEN - 1);
        u.ptr = null;

        if (type == GretlType.matrix) {
            GretlMatrixRaw * m = cast(GretlMatrixRaw *) value;

            if (m == null) {
                u.ptr = gretl_null_matrix_new();
            } else if (get_user_var_by_data(m) != null) {
                /* this check should be redundant? */
                u.ptr = gretl_matrix_copy(m);
            } else {
                u.ptr = value;
            }
        } else if (type == GretlType.bundle) {
            if (value == null) {
                u.ptr = gretl_bundle_new();
            } else {
                u.ptr = value;
            }
        } else if (type == GretlType.string_) {
            if (value == null) {
                u.ptr = gretl_strdup("");
            } else {
                u.ptr = value;
            }
        } else if (type == GretlType.list) {
            if (value == null) {
                u.ptr = gretl_null_list();
            } else {
                u.ptr = value;
            }
        } else if (type == GretlType.double_) {
            if (value == null) {
                u.ptr = na_ptr();
            } else {
                u.ptr = value;
            }
        } else if (gretl_array_type(type) || type == GretlType.any) {
            if (value == null) {
                u.ptr = gretl_array_new(type, 0, err);
            } else {
                u.ptr = value;
            }
            u.type = GretlType.array;
        } else {
            fprintf(stderr, "user_var_new error, type=%d (%s)\n".ptr, type,
                    gretl_type_get_name(type));
            *err = GretlError.data;
        }
    }

    if (u.ptr == null) {
        if (!*err) {
            *err = GretlError.alloc;
        }
        free(u);
        u = null;
    }

    return u;
}

void uvar_free_value(UserVar * u)
{
    if (u.ptr == null) {
        return;
    } else if (u.type == GretlType.matrix) {
        gretl_matrix_free(cast(GretlMatrixRaw *) u.ptr);
    } else if (u.type == GretlType.bundle) {
        gretl_bundle_destroy(cast(GretlBundle *) u.ptr);
    } else if (u.type == GretlType.string_) {
        bufgets_finalize(cast(const(char) *) u.ptr);
        free(u.ptr);
    } else if (u.type == GretlType.array) {
        gretl_array_destroy(cast(GretlArray *) u.ptr);
    } else {
        /* scalar, list */
        free(u.ptr);
    }
}

__gshared GHashTable * uvh0 = null;       /* for use at "main" exec level */
__gshared GHashTable * uvh1 = null;       /* for use within functions */
__gshared GHashTable * uvars_hash = null; /* pointer to one or other of the above */
__gshared int previous_d = -1;    /* record of previous "function depth" */

void set_previous_depth(int d)
{
    previous_d = d;
}

int get_previous_depth()
{
    return previous_d;
}

void switch_uservar_hash(int level)
{

    if (level == 0) {
        uvars_hash = uvh0;
        if (uvh1 != null) {
            g_hash_table_remove_all(uvh1);
        }
    } else {
        uvars_hash = uvh1;
    }
}

void uvar_hash_destroy()
{
    if (uvh0 != null) {
        g_hash_table_destroy(uvh0);
        uvh0 = null;
    }

    if (uvh1 != null) {
        g_hash_table_destroy(uvh1);
        uvh1 = null;
    }

    /* also NULL the convenience pointer */
    uvars_hash = null;

    set_previous_depth(-1);
}

void user_var_destroy(UserVar * u)
{

    if (uvars_hash != null) {
        g_hash_table_remove(uvars_hash, u.name.ptr);
    }

    if (!var_is_shell(u)) {
        uvar_free_value(u);
    }

    free(u);
}

int resize_uvar_stack(int n)
{
    int err = 0;

    if (n > n_alloc) {
        int n_new = n_alloc + UV_CHUNK;
        UserVar ** tmp;

        tmp = cast(typeof(tmp)) realloc(uvars, n_new * (*tmp).sizeof);
        if (tmp == null) {
            err = GretlError.alloc;
        } else {
            uvars = tmp;
            n_alloc = n_new;
        }
    }

    return err;
}

void set_nvars(int n, const(char) * caller)
{
    n_vars = n;
}

int bname_is_temp(const(char) * name)
{
    return !strncmp(name, "btmp___", 7) && isdigit(name[7]);
}

int real_user_var_add(const(char) * name, GretlType type, void * value, Gretlopt opt, UserVar ** pu)
{
    UserVar * u;
    int err = 0;

    u = user_var_new(name, type, value, &err);

    if (u == null) {
        fprintf(stderr, "real_user_var_add: name='%s', value=%p, u=NULL\n".ptr,
                name, value);
    if (pu != null) {
        *pu = null;
    }
        return err ? err : GretlError.data;
    }

    /* We use OPT_P for a private variable, OPT_A when adding as a
       function argument, OPT_S when adding as a "shell" variable,
       OPT_C when auto-casting a 1 x 1 matrix result to a scalar.
    */

    if (!err) {
        err = resize_uvar_stack(n_vars + 1);
        if (!err) {
            if (opt & Gretlopt.OPT_P) {
                u.flags = UVFlags.private_;
            }
        if (opt & Gretlopt.OPT_S) {
                u.flags |= UVFlags.shell;
            }
            if (opt & Gretlopt.OPT_A) {
                u.flags &= ~UVFlags.main;
                u.level += 1;
            }
            if (opt & Gretlopt.OPT_C) {
                u.flags |= UVFlags.nodecl;
            }
            uvars[n_vars] = u;
            set_nvars(n_vars + 1, "user_var_add");
        }
    }

    if (!err && user_var_callback != null && u.level == 0 &&
        !(opt & (Gretlopt.OPT_P | Gretlopt.OPT_S)) && *name != '$' &&
        (type == GretlType.matrix || type == GretlType.bundle) &&
        !(type == GretlType.bundle && bname_is_temp(name))) {
        return user_var_callback(name, type, UvarAction.add);
    }

    if (pu != null) {
    *pu = u;
    }

    return err;
}

/**
 * user_var_add:
 * @name: name to give the variable.
 * @type: the type of the variable.
 * @value: pointer to value for variable.
 *
 * Adds a new user-variable with the given characteristics.
 * Note that the user-variable takes ownership of the
 * supplied @value; this should be copied first if need be.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int user_var_add(const(char) * name, GretlType type, void * value)
{
    return real_user_var_add(name, type, value, Gretlopt.OPT_NONE, null);
}

/**
 * alt_user_var_add:
 * @name: name to give the variable.
 * @type: the type of the variable.
 * @value: pointer to value for variable.
 *
 * Adds a new UserVar with the given characteristics.
 * Note that the UserVar takes ownership of the supplied
 * @value; this should be copied first if need be.
 *
 * Returns: pointer to the new UserVar, or NULL on failure.
 */

UserVar * alt_user_var_add(const(char) * name, GretlType type, void * value)
{
    UserVar * uv = null;

    real_user_var_add(name, type, value, Gretlopt.OPT_NONE, &uv);
    return uv;
}

int private_matrix_add(GretlMatrixRaw * M, const(char) * name)
{
    return real_user_var_add(name, GretlType.matrix,
                 M, Gretlopt.OPT_P, null);
}

int private_matrix_add_as_shell(GretlMatrixRaw * M, const(char) * name)
{
    return real_user_var_add(name, GretlType.matrix,
                 M, Gretlopt.OPT_P | Gretlopt.OPT_S, null);
}

int private_bundle_add_as_shell(GretlBundle * b, const(char) * name)
{
    return real_user_var_add(name, GretlType.bundle,
                 b, Gretlopt.OPT_P | Gretlopt.OPT_S, null);
}

int private_scalar_add(double val, const(char) * name)
{
    double * px = cast(double *) malloc(double.sizeof);
    int err;

    if (px == null) {
        err = GretlError.alloc;
    } else {
        *px = val;
        err = real_user_var_add(name, GretlType.double_,
                                px, Gretlopt.OPT_P, null);
    }

    return err;
}

/**
 * user_var_delete_by_name:
 * @name: name of the variable to delete.
 * @prn: pointer to gretl printer, or NULL.
 *
 * Deletes the specified user-variable.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int user_var_delete_by_name(const(char) * name, Prn * prn)
{
    GretlType type = cast(GretlType) 0;
    int level = gretl_function_depth();
    UserVar * targ = null;
    int i, j, k = 0;
    int err = 0;

    for (i=0; i<n_vars; i++) {
        if (uvars[i].level == level && !strcmp(uvars[i].name.ptr, name)) {
            targ = uvars[i];
            k = i;
            break;
        }
    }

    if (targ == null) {
        return GretlError.unkvar;
    }

    if (level > 0 && (targ.flags & UVFlags.main)) {
        gretl_errmsg_sprintf(_("%s: cannot be deleted here").ptr, targ.name.ptr);
        return GretlError.data;
    }

    if (user_var_callback != null && level == 0 &&
        !var_is_private(targ) &&
        (targ.type == GretlType.matrix ||
         targ.type == GretlType.bundle)) {
        /* run this deletion through the GUI program to ensure
           that things stay in sync
        */
        return user_var_callback(name, targ.type, UvarAction.delete_);
    }

    type = targ.type;
    user_var_destroy(targ);
    for (j=k; j<n_vars-1; j++) {
        uvars[j] = uvars[j+1];
    }
    resize_uvar_stack(n_vars - 1);
    set_nvars(n_vars - 1, "user_var_delete_by_name");

    if (prn != null && gretl_messages_on()) {
        pprintf(prn, _("Deleted %s").ptr, name);
        pputc(prn, '\n');
    }
    if (level == 0 && type == GretlType.double_ &&
        scalar_edit_callback != null) {
        scalar_edit_callback();
    }

    return err;
}

int user_var_delete(UserVar * uvar)
{
    int i, j, err = GretlError.unkvar;

    for (i=0; i<n_vars; i++) {
        if (uvar == uvars[i]) {
            user_var_destroy(uvars[i]);
            for (j=i; j<n_vars-1; j++) {
                uvars[j] = uvars[j+1];
            }
            set_nvars(n_vars - 1, "user_var_delete");
            err = 0;
            break;
        }
    }

    return err;
}

/* Try to guess whether the currently-called function is big enough
   (number of lines of code) to make it worthwhile to construct a hash
   table for uservars at its level of execution, namely @uvh1, given
   that we'll have to empty the table on exit from the function.

   The number-of-lines threshold here is obviously kinda arbitrary;
   some systematic experimentation might be useful.
*/

int use_uvh1()
{
    /* maybe add: && !gretl_function_recursing() */
    return current_function_size() > 40;
}

UserVar * get_user_var_of_type_by_name(const(char) * name, GretlType type)
{
    int prev_d = get_previous_depth();
    int d = gretl_function_depth();
    int i, imin = 0;
    UserVar * u = null;

    if (name == null || *name == '\0') {
        return null;
    }

    if (type == GretlType.double_) {
        /* support "auxiliary scalars" mechanism */
        imin = scalar_imin;
    }

    if (d != prev_d) {
        if (d == 0) {
            /* we're now at "main" level */
            if (uvh0 == null) {
                uvh0 = g_hash_table_new(&g_str_hash, &g_str_equal);
            }
            if (uvh1 != null) {
                g_hash_table_remove_all(uvh1);
            }
            uvars_hash = uvh0;
        } else if (!use_uvh1()) {
            /* exec'ing a function, hash table not wanted */
            if (prev_d > 0 && uvh1 != null) {
                g_hash_table_remove_all(uvh1);
            }
            uvars_hash = null;
        } else {
            /* exec'ing a function, hash table wanted */
            if (uvh1 == null) {
                uvh1 = g_hash_table_new(&g_str_hash, &g_str_equal);
            } else if (prev_d > 0 && uvh1 != null) {
                g_hash_table_remove_all(uvh1);
            }
            uvars_hash = uvh1;
        }
        set_previous_depth(d);
    }

    if (uvars_hash != null) {
        /* first resort: try a hash look-up */
        u = cast(typeof(u)) g_hash_table_lookup(uvars_hash, name);
        /* but verify type, if specified */
        if (u != null && type != GretlType.any && u.type != type) {
            u = null;
        }
    }

    if (u == null) {
        /* "On demand" hashing: if we're successful in looking
           up a variable in the traditional manner, then
           insert it into the uservars hash table.
        */
        for (i=imin; i<n_vars; i++) {
            if (uvars[i].level == d &&
                (type == GretlType.any || uvars[i].type == type) &&
                !strcmp(uvars[i].name.ptr, name)) {
                u = uvars[i];
                if (uvars_hash != null) {
                    g_hash_table_insert(uvars_hash, u.name.ptr, u);
                }
                break;
            }
        }
    }

    return u;
}

UserVar * get_user_var_by_name(const(char) * name)
{
    return get_user_var_of_type_by_name(name, GretlType.any);
}

GretlType user_var_get_type_by_name(const(char) * name)
{
    UserVar * u;

    u = get_user_var_of_type_by_name(name, GretlType.any);

    return u == null ? GretlType.none : u.type;
}

GretlType user_var_get_specific_type(const(char) * name)
{
    UserVar * u;

    u = get_user_var_of_type_by_name(name, GretlType.any);

    if (u != null && u.type == GretlType.array) {
    return gretl_array_get_type(cast(GretlArray *) u.ptr);
    } else {
    return u == null ? GretlType.none : u.type;
    }
}

void * user_var_get_value_and_type(const(char) * name, GretlType * type)
{
    void * ret = null;
    UserVar * u;

    u = get_user_var_of_type_by_name(name, GretlType.any);

    if (u != null) {
        ret = u.ptr;
        *type = u.type;
    } else {
        *type = GretlType.none;
    }

    return ret;
}

const(char) * uservar_name_complete(const(char) * s)
{
    const(char) * ret = null;

    if (uvars_hash != null) {
        GList *hk = g_hash_table_get_keys(uvars_hash);
        int n = cast(int) (cast(int) strlen(s));

        while (hk != null) {
            if (!strncmp(cast(const(char) *) hk.data, s, n)) {
                ret = cast(const(char) *) hk.data;
                break;
            }
            hk = hk.next;
        }
        g_list_free(hk);
    }

    return ret;
}

int gretl_is_user_var(const(char) * name)
{
    return get_user_var_by_name(name) != null;
}

UserVar * get_user_var_by_data(const(void) * data)
{
    int i, d = gretl_function_depth();

    if (data == null) {
        return null;
    }

    for (i=0; i<n_vars; i++) {
        if (uvars[i] != null && uvars[i].level == d &&
            uvars[i].ptr == data) {
            return uvars[i];
        }
    }

    return null;
}

const(char) * user_var_get_name(UserVar * uvar)
{
    return uvar == null ? null : uvar.name.ptr;
}

const(char) * user_var_get_name_by_data(const(void) * data)
{
    UserVar * u = get_user_var_by_data(data);

    return u == null ? null : u.name.ptr;
}

int user_var_get_level(UserVar * uvar)
{
    return (uvar == null)? -1 : uvar.level;
}

int user_var_get_flags(UserVar * uvar)
{
    return (uvar == null)? 0 : cast(int) uvar.flags;
}

int user_var_set_flag(UserVar * uvar, UVFlags flag)
{
    if (uvar != null) {
        uvar.flags |= flag;
        return 0;
    } else {
        return GretlError.invarg;
    }
}

int user_var_unset_flag(UserVar * uvar, UVFlags flag)
{
    if (uvar != null) {
        uvar.flags &= ~flag;
        return 0;
    } else {
        return GretlError.invarg;
    }
}

void user_var_privatize_by_name(const(char) * name)
{
    UserVar * u = get_user_var_by_name(name);

    if (u != null) {
        u.flags |= UVFlags.private_;
    }
}

void * user_var_get_value(UserVar * uvar)
{
    return (uvar == null)? null : uvar.ptr;
}

GretlType user_var_get_type(UserVar * uvar)
{
    return (uvar == null)? cast(GretlType) 0 : uvar.type;
}

void * user_var_get_value_by_name(const(char) * name)
{
    UserVar * u = get_user_var_by_name(name);

    return (u == null)? null : u.ptr;
}

/* special for scalars since user_var_get_value returns
   a pointer */

double user_var_get_scalar_value(UserVar * uvar)
{
    if (uvar != null && uvar.type == GretlType.double_) {
        return *cast(double *) uvar.ptr;
    } else {
        return NADBL;
    }
}

int user_var_set_scalar_value(UserVar * uvar, double x)
{
    if (uvar != null && uvar.type == GretlType.double_) {
        *cast(double *) uvar.ptr = x;
        return 0;
    } else {
        return GretlError.data;
    }
}

int user_var_adjust_level(UserVar * uvar, int adj)
{
    if (uvar == null) {
        return GretlError.unkvar;
    } else {
        uvar.level += adj;
        return 0;
    }
}

/* Note: the following should be called only from internal
   contexts in which we know that the attempted renaming
   is not broken (e.g. trying to assign to @uvar a name
   that is already taken by some other object).
*/

int user_var_set_name(UserVar * uvar, const(char) * name)
{
    int err = 0;

    if (uvar == null) {
        err = GretlError.data;
    } else {
        *uvar.name.ptr = '\0';
        strncat(uvar.name.ptr, name, VNAMELEN - 1);
    }

    return err;
}

/**
 * user_var_localize:
 * @origname: name of variable at caller level.
 * @localname: name to be used within function.
 * @type: type of the variable.
 *
 * On entry to a function, renames the named variable (provided
 * as an argument) and sets its level so that is is accessible
 * within the function.
 *
 * Returns: 0 on success, non-zero on error.
 */

int user_var_localize(const(char) * origname, const(char) * localname, GretlType type)
{
    UserVar * u;
    int err = 0;

    if (gretl_is_array_ref_type(type)) {
        type = GretlType.array;
    } else {
        type = gretl_type_get_plain_type(type);
    }

    u = get_user_var_of_type_by_name(origname, type);

    if (u == null) {
        err = GretlError.data;
    } else {
        user_var_set_name(u, localname);
        u.level += 1;
    }

    return err;
}

int user_var_count_for_type(GretlType type)
{
    int i, n = 0;

    for (i=0; i<n_vars; i++) {
        if (uvars[i].type == type) {
            n++;
        }
    }

    return n;
}

int n_user_matrices()
{
    return user_var_count_for_type(GretlType.matrix);
}

int n_user_scalars()
{
    return user_var_count_for_type(GretlType.double_);
}

int n_user_lists()
{
    return user_var_count_for_type(GretlType.list);
}

int n_user_bundles()
{
    return user_var_count_for_type(GretlType.bundle);
}

/**
 * user_var_replace_value:
 * @uvar: user variable.
 * @value: the new value to place as the value or @uvar.
 * @type: the type of the replacement value.
 *
 * Replaces the value of @uvar; the existing value is
 * freed first.
 *
 * Returns: 0 on success, non-zero on error.
 */

int user_var_replace_value(UserVar * uvar, void * value, GretlType type)
{
    int err = 0;

    if (uvar == null) {
        err = GretlError.unkvar;
    } else if (value != uvar.ptr && (uvar.flags & UVFlags.norepl)) {
        gretl_errmsg_sprintf(_("The variable %s is read-only").ptr, uvar.name.ptr);
        err = GretlError.data;
    } else if (type != uvar.type) {
        err = GretlError.types; /* assume the worst */
        if (uvar.type == GretlType.array && uvar.ptr != null) {
            /* but we might be OK */
            if (type == gretl_array_get_type(cast(GretlArray *) uvar.ptr)) {
                err = 0;
            }
        }
        if (err) {
            fputs("*** user_var_replace_value: type mismatch ***\n", stderr);
            fprintf(stderr, " (expected %s but got %s)\n".ptr,
                    gretl_type_get_name(uvar.type), gretl_type_get_name(type));
        }
    }

    if (!err && value != uvar.ptr) {
        if (uvar.ptr != null) {
            uvar_free_value(uvar);
        }
        uvar.ptr = value;
    }

    return err;
}

/**
 * user_var_set_pointer:
 * @uvar: user variable.
 * @ptr: the new pointer to set on @uvar.
 *
 * Unlike @user_var_replace_value, this function does not free
 * the existing data pointer on @uvar. It should therefore be
 * used only when another pointer to the prior data exists,
 * otherwise a memory leak would result. Plus there's no type
 * checking so use with care.
 *
 * Returns: 0 on success, error code if @uvar is NULL.
 */

int user_var_set_pointer(UserVar * uvar, void * ptr)
{
    if (uvar != null) {
    uvar.ptr = ptr;
    return 0;
    } else {
    return GretlError.data;
    }
}

char * user_string_resize(const(char) * name, size_t len, int * err)
{
    UserVar * u;

    u = get_user_var_of_type_by_name(name, GretlType.string_);

    if (u == null) {
        *err = GretlError.invarg;
        return null;
    } else {
        char * orig = cast(char *) u.ptr;

        if (orig == null || len > strlen(orig) + 1) {
            char * tmp = cast(char *) realloc(u.ptr, len);

            if (tmp == null) {
                *err = GretlError.alloc;
            } else {
                u.ptr = tmp;
            }
        }
    }

    return cast(char *) u.ptr;
}

char * user_string_reset(const(char) * name, const(char) * repl, int * err)
{
    UserVar * u;

    u = get_user_var_of_type_by_name(name, GretlType.string_);

    if (u == null) {
        *err = GretlError.invarg;
        return null;
    } else {
        free(u.ptr);
        if (repl == null) {
            u.ptr = gretl_strdup("");
        } else {
            u.ptr = gretl_strdup(repl);
        }
        return cast(char *) u.ptr;
    }
}

int check_array_type_compat(GretlType type, UserVar * u)
{
    int err = 0;

    if (u.type != GretlType.array) {
        err = GretlError.types;
    } else {
        /* we also need a more specific check here */
        if (type != gretl_array_get_type(cast(GretlArray *) u.ptr)) {
            err = GretlError.types;
        }
    }

    return err;
}

int user_var_add_or_replace(const(char) * name, GretlType type, void * value)
{
    UserVar * u = get_user_var_by_name(name);
    int err = 0;

    if (u != null) {
        if (gretl_array_type(type)) {
            err = check_array_type_compat(type, u);
        } else if (u.type != type) {
            err = GretlError.types;
        }
        if (!err) {
            err = user_var_replace_value(u, value, type);
        }
    } else {
        err = real_user_var_add(name, type, value, Gretlopt.OPT_NONE, null);
    }

    return err;
}

void * user_var_steal_value(UserVar * uvar)
{
    void * ret = null;

    if (uvar != null) {
        ret = uvar.ptr;
        uvar.ptr = null;
    }

    return ret;
}

/* FIXME: are both the above and the below necessary? */

void * user_var_unstack_value(UserVar * uvar)
{
    void * ret = null;
    int i, j;

    for (i=0; i<n_vars; i++) {
        if (uvar == uvars[i]) {
            ret = uvar.ptr;
            uvars[i].ptr = null;
            user_var_destroy(uvars[i]);
            for (j=i; j<n_vars-1; j++) {
                uvars[j] = uvars[j+1];
            }
            set_nvars(n_vars - 1, "user_var_unstack_value");
            break;
        }
    }

    return ret;
}

int user_matrix_replace_matrix_by_name(const(char) * name, GretlMatrixRaw * m)
{
    UserVar * u = get_user_var_by_name(name);

    if (u != null) {
        return user_var_replace_value(u, m, GretlType.matrix);
    } else {
        return GretlError.data;
    }
}

GList * user_var_names_for_type(GretlType type)
{
    GList *list = null;
    int i;

    for (i=0; i<n_vars; i++) {
        if (uvars[i].type == type) {
            list = g_list_append(list, cast(gpointer) uvars[i].name.ptr);
        }
    }

    return list;
}

GList * user_var_list_for_type(GretlType type)
{
    GList *list = null;
    int i;

    for (i=0; i<n_vars; i++) {
        if (uvars[i].type == type) {
            list = g_list_append(list, cast(gpointer) uvars[i]);
        }
    }

    return list;
}

/**
 * set_user_var_callback:
 * @callback: function function to put in place.
 *
 * Sets the callback function to be invoked when a user-defined
 * matrix is added to or removed from the stack of saved objects.
 * Intended for synchronizing the GUI program with the saved object
 * state.
 */

void set_user_var_callback(USER_VAR_FUNC callback)
{
    user_var_callback = callback;
}

void set_scalar_edit_callback(void function() callback)
{
    scalar_edit_callback = callback;
}

/**
 * arg_add_as_shell:
 * @name: the name to be given to the "shell" variable.
 * @type: the type of the variable.
 * @value: the value pointer
 *
 * The value in question is added to the stack of named
 * variables under the name @name with the shell flag
 * set. This is used (a) when an anonymous matrix is given as
 * a %const argument to a user-defined function and (b) when
 * an anonymous bundle is given as the argument corresponding
 * to a bundle-pointer parameter. The @value becomes
 * accessible by @name within the function, but is protected
 * from destruction on exit from the function.
 *
 * Returns: 0 on success, non-zero on error.
 */

int arg_add_as_shell(const(char) * name, GretlType type, void * value)
{
    return real_user_var_add(name, type, value, Gretlopt.OPT_S | Gretlopt.OPT_A, null);
}

/**
 * copy_matrix_as:
 * @m: the original matrix.
 * @newname: the name to be given to the copy.
 * @fnarg: 0 for regular use.
 *
 * A copy of matrix @m is added to the stack of saved matrices
 * under the name @newname.
 *
 * The @fnarg argument should be non-zero only if this function
 * is used to handle the case where a matrix is given as the argument
 * to a user-defined function.
 *
 * Returns: 0 on success, non-zero on error.
 */

int copy_matrix_as(const(GretlMatrixRaw) * m, const(char) * newname, int fnarg)
{
    GretlMatrixRaw * m2 = gretl_matrix_copy(m);
    int err = 0;

    if (m2 == null) {
        err = GretlError.alloc;
    } else {
        Gretlopt opt = fnarg ? Gretlopt.OPT_A : Gretlopt.OPT_NONE;

        err = real_user_var_add(newname, GretlType.matrix,
                m2, opt, null);
    }

    return err;
}

int copy_as_arg(const(char) * param_name, GretlType type, void * value)
{
    void * copyval = null;
    GretlType cpytype = type;
    int err = 0;

    if (type == GretlType.matrix) {
        GretlMatrixRaw * mcpy = gretl_matrix_copy(cast(GretlMatrixRaw *) value);

        if (mcpy == null) {
            err = GretlError.alloc;
        } else {
            copyval = mcpy;
        }
    } else if (type == GretlType.list) {
        int * lcpy = gretl_list_copy(cast(int *) value);

        if (lcpy == null) {
            err = GretlError.alloc;
        } else {
            copyval = lcpy;
        }
    } else if (type == GretlType.string_) {
        char * scpy = gretl_strdup(cast(char *) value);

        if (scpy == null) {
            err = GretlError.alloc;
        } else {
            copyval = scpy;
        }
    } else if (type == GretlType.double_) {
        double * px = cast(double *) malloc(double.sizeof);

        if (px == null) {
            err = GretlError.alloc;
        } else {
            *px = *cast(double *) value;
            copyval = px;
        }
    } else if (type == GretlType.bundle) {
        GretlBundle * bcpy = gretl_bundle_copy(cast(GretlBundle *) value, &err);

        if (!err) {
            copyval = bcpy;
        }
    } else if (gretl_array_type(type)) {
        GretlArray * acpy = gretl_array_copy(cast(const(GretlArray) *) value, &err);

        if (!err) {
            copyval = acpy;
            cpytype = gretl_array_get_type(acpy);
        }
    }

    if (!err) {
        err = real_user_var_add(param_name, cpytype,
                copyval, Gretlopt.OPT_A, null);
    }

    return err;
}

int * copy_list_as_arg(const(char) * param_name, int * list, int * err)
{
    int * ret = null;

    *err = copy_as_arg(param_name, GretlType.list, list);
    if (!*err) {
        ret = cast(typeof(ret)) uvars[n_vars-1].ptr;
    }

    return ret;
}

void destroy_user_vars()
{
    int i, j;

    for (i=0; i<n_vars; i++) {
        if (uvars[i] == null) {
            break;
        }
        user_var_destroy(uvars[i]);
        for (j=i; j<n_vars-1; j++) {
            uvars[j] = uvars[j+1];
        }
        uvars[n_vars-1] = null;
        i--;
    }

    if (uvh0 != null || uvh1 != null) {
        uvar_hash_destroy();
    }

    set_nvars(0, "destroy_user_vars");

    free(uvars);
    uvars = null;
    n_alloc = 0;
}

int uvar_levels_match(UserVar * u, int level)
{
    int ret = 0;

    if (u.level == level) {
        ret = 1;
    } else if (level == LEV_PRIVATE && var_is_private(u)) {
        ret = 1;
    }

    return ret;
}

int real_destroy_user_vars_at_level(int level, int type, int imin)
{
    int i, j, nv = imin;
    int err = 0;

    for (i=imin; i<n_vars; i++) {
        if (uvars[i] == null) {
            break;
        }
        if (type > 0 && uvars[i].type != type) {
            /* preserve this variable */
            nv++;
            continue;
        }
        if (uvar_levels_match(uvars[i], level)) {
            user_var_destroy(uvars[i]);
            /* shuffle the remainder down one place */
            for (j=i; j<n_vars-1; j++) {
                uvars[j] = uvars[j+1];
            }
            uvars[n_vars-1] = null;
            i--;
        } else {
            /* preserving */
            nv++;
        }
    }

    set_nvars(nv, "real_destroy_user_vars_at_level");

    return err;
}

int destroy_user_vars_via_callback(int type)
{
    UserVar ** delvars = null;
    int i, j, n = 0;
    int err = 0;

    for (i=0; i<n_vars; i++) {
        if (uvars[i].level == 0 && uvars[i].type == type) {
            n++;
        }
    }

    if (n == 0) {
        return 0;
    }

    delvars = cast(typeof(delvars)) malloc(n * (*delvars).sizeof);
    if (delvars == null) {
        return GretlError.alloc;
    }

    j = 0;
    for (i=0; i<n_vars; i++) {
        if (uvars[i].level == 0 && uvars[i].type == type) {
            delvars[j++] = uvars[i];
        }
    }

    for (j=0; j<n && !err; j++) {
        err = user_var_callback(delvars[j].name.ptr,
                                   delvars[j].type,
                                   UvarAction.delete_);
    }

    free(delvars);

    return err;
}

/**
 * destroy_user_vars_at_level:
 * @level: stack level of function execution.
 *
 * Destroys and removes from the stack of user matrices all
 * matrices that were created at the given @level.  This is
 * part of the cleanup that is performed when a user-defined
 * function terminates.
 *
 * Returns: 0 on success, non-zero on error.
 */

int destroy_user_vars_at_level(int level)
{
    return real_destroy_user_vars_at_level(level, 0, 0);
}

int destroy_private_uvars()
{
    return real_destroy_user_vars_at_level(LEV_PRIVATE, 0, 0);
}

int destroy_private_matrices()
{
    return real_destroy_user_vars_at_level(LEV_PRIVATE,
                                           GretlType.matrix,
                                           0);
}

int destroy_private_lists()
{
    return real_destroy_user_vars_at_level(LEV_PRIVATE,
                                           GretlType.list,
                                           0);
}

int delete_user_vars_of_type(GretlType type, Prn * prn)
{
    int err = 0;

    if (type == GretlType.matrix ||
        type == GretlType.bundle ||
        type == GretlType.array  ||
        type == GretlType.string_ ||
        type == GretlType.double_ ||
        type == GretlType.list) {
        int level = gretl_function_depth();

        if (level == 0 && user_var_callback != null &&
            (type == GretlType.matrix || type == GretlType.bundle)) {
            err = destroy_user_vars_via_callback(type);
        } else {
            err = real_destroy_user_vars_at_level(level, type, 0);
        }

        if (!err && gretl_messages_on()) {
            pprintf(prn, _("Deleted all variables of type %s\n").ptr,
                    gretl_type_get_name(type));
        }
    } else {
        err = GretlError.types;
    }

    return err;
}

/**
 * destroy_private_scalars:
 *
 * Gets rid of private or "internal" scalars whose
 * names begin with '$'.
 */

void destroy_private_scalars()
{
    real_destroy_user_vars_at_level(LEV_PRIVATE,
                                    GretlType.double_,
                                    0);
}

char * temp_name_for_bundle()
{
    char[VNAMELEN] tmpname;
    int i, nb = 0;

    for (i=0; i<n_vars; i++) {
        if (uvars[i].type == GretlType.bundle) {
            nb++;
        }
    }

    sprintf(tmpname.ptr, "btmp___%d".ptr, nb);
    return gretl_strdup(tmpname.ptr);
}

int xml_put_user_matrix(UserVar * u, Prn * prn)
{
    int err = 0;

    if (u != null && u.ptr != null) {
        err = gretl_matrix_serialize(cast(const(GretlMatrixRaw) *) u.ptr, u.name.ptr, prn);
    }

    return err;
}

void write_scalar_value(double x, const(char) * fmt, Prn * prn)
{
    if (na(x)) {
        pprintf(prn, "%g\n".ptr, x);
    } else {
        pprintf(prn, fmt, x);
    }
}

void serialize_scalar_value(double x, Prn * prn)
{
    if (na(x)) {
        pprintf(prn, "%g".ptr, x);
    } else {
        pprintf(prn, "%.16g".ptr, x);
    }
}

/**
 * print_scalars:
 * @prn: pointer to gretl printing struct.
 *
 * Prints names and values of any saved scalars.
 */

void print_scalars(Prn * prn)
{
    double x;
    int level = gretl_function_depth();
    int len, ns = 0, maxlen = 0;
    int i;

    for (i=0; i<n_vars; i++) {
        if (uvars[i].type == GretlType.double_ &&
            uvars[i].level == level) {
            len = cast(int) (strlen(uvars[i].name.ptr));
            if (len > maxlen) {
                maxlen = len;
            }
            ns++;
        }
    }

    if (ns == 0) {
        pprintf(prn, "%s\n".ptr, _("none").ptr);
        return;
    }

    pputc(prn, '\n');

    for (i=0; i<n_vars; i++) {
        if (uvars[i].type == GretlType.double_ &&
            uvars[i].level == level) {
            x = *cast(double *) uvars[i].ptr;
            pprintf(prn, " %*s = ".ptr, maxlen, uvars[i].name.ptr);
            write_scalar_value(x, "%.16g\n", prn);
        }
    }

    pputc(prn, '\n');
}

void print_scalar_by_name(const(char) * name, Prn * prn)
{
    UserVar * u;

    u = get_user_var_of_type_by_name(name, GretlType.double_);

    if (u != null) {
        double x = *cast(double *) u.ptr;

        pprintf(prn, "\n%15s = ".ptr, u.name.ptr);
        write_scalar_value(x, "% #.8g\n", prn);
    }
}

/* "auxiliary scalars": this apparatus is used when we want to do
   "private" NLS estimation (e.g. in ARMA initialization).  It ensures
   that the scalar NLS parameters don't collide with the public scalar
   namespace. FIXME.
*/

void set_auxiliary_scalars()
{
    scalar_imin = n_vars;
}

void unset_auxiliary_scalars()
{
    real_destroy_user_vars_at_level(0, GretlType.double_, scalar_imin);
    scalar_imin = 0;
}

int real_scalar_add(const(char) * name, double val, Gretlopt opt)
{
    UserVar * u = get_user_var_by_name(name);
    int level = gretl_function_depth();
    int err = 0;

    if (u != null) {
        if (u.type == GretlType.double_) {
            *cast(double *) u.ptr = val;
        } else {
            err = GretlError.types;
        }
        return err;
    } else {
        double * px = cast(double *) malloc(double.sizeof);

        if (px == null) {
            err = GretlError.alloc;
        } else {
            *px = val;
            err = real_user_var_add(name, GretlType.double_,
                                    px, opt, null);
        }

        if (!err && level == 0 && scalar_edit_callback != null) {
            scalar_edit_callback();
        }
    }

    return err;
}

int gretl_scalar_add(const(char) * name, double val)
{
    return real_scalar_add(name, val, Gretlopt.OPT_NONE);
}

int gretl_scalar_add_mutable(const(char) * name, double val)
{
    return real_scalar_add(name, val, Gretlopt.OPT_C);
}

int gretl_scalar_convert_to_matrix(UserVar * u)
{
    GretlMatrixRaw * m = null;

    if (u == null) {
        return GretlError.unkvar;
    } else if (u.type != GretlType.double_) {
        return GretlError.types;
    }

    m = gretl_matrix_alloc(1, 1);
    if (m == null) {
        return GretlError.alloc;
    }

    m.val[0] = *cast(double *) u.ptr;
    free(u.ptr);
    u.ptr = m;
    u.type = GretlType.matrix;

    if (gretl_function_depth() == 0) {
        if (scalar_edit_callback != null) {
            scalar_edit_callback();
        }
        if (user_var_callback != null) {
            user_var_callback(u.name.ptr, GretlType.matrix, UvarAction.add);
        }
    }

    return 0;
}

int add_auxiliary_scalar(const(char) * name, double val)
{
    double * px = cast(double *) malloc(double.sizeof);
    int err;

    /* Note that unlike gretl_scalar_add() above, this function
       adds a new_ scalar variable unconditionally; it never
       modifies the value of an existing scalar of the same
       name.
    */

    if (px == null) {
        err = GretlError.alloc;
    } else {
        *px = val;
        err = real_user_var_add(name, GretlType.double_,
                                px, Gretlopt.OPT_NONE, null);
    }

    return err;
}

int gretl_scalar_set_value(const(char) * name, double val)
{
    UserVar * u;
    int err = 0;

    u = get_user_var_of_type_by_name(name, GretlType.double_);

    if (u == null) {
        gretl_errmsg_sprintf(_("%s: no such scalar").ptr, name);
        err = GretlError.data;
    } else if (scalar_is_read_only_index(name)) {
        err = GretlError.data;
        gretl_errmsg_sprintf(_("The variable %s is currently read-only").ptr, name);
    } else {
        *cast(double *) u.ptr = val;

        if (scalar_edit_callback != null) {
            scalar_edit_callback();
        }
    }

    return err;
}

/* get the value from a user variable of scalar type */

double gretl_scalar_get_value(const(char) * name, int * err)
{
    UserVar * u;
    double ret = NADBL;

    u = get_user_var_of_type_by_name(name, GretlType.double_);

    if (u != null) {
        ret = *cast(double *) u.ptr;
    } else {
        ret = get_const_by_name(name, err);
    }

    return ret;
}

double maybe_get_bundled_scalar(const(char) * name, int * err)
{
    const(char) * p = strchr(name, '.');
    GretlBundle * b = null;
    char[VNAMELEN] bname;
    char[VNAMELEN] key;
    double x = NADBL;

    *bname.ptr = '\0';
    strncat(bname.ptr, name, p - name);
    b = get_bundle_by_name(bname.ptr);

    if (b == null) {
        *err = GretlError.invarg;
    } else {
        *key.ptr = '\0';
        strncat(key.ptr, p + 1, VNAMELEN - 1);
        x = gretl_bundle_get_scalar(b, key.ptr, err);
    }

    return x;
}

/* more "permissive" than gretl_scalar_get_value(): allows
   for @name being the identifier for a 1 x 1 matrix, or
   bundle.member
*/

double get_scalar_value_by_name(const(char) * name, int * err)
{
    double ret = NADBL;
    UserVar * u;

    if (strchr(name, '.')) {
        ret = maybe_get_bundled_scalar(name, err);
        goto bailout;
    }

    u = get_user_var_by_name(name);

    if (u != null) {
        if (u.type == GretlType.double_) {
            ret = *cast(double *) u.ptr;
        } else if (u.type == GretlType.matrix) {
            GretlMatrixRaw * m = cast(GretlMatrixRaw *) u.ptr;

            if (gretl_matrix_is_scalar(m)) {
                ret = m.val[0];
            } else {
                *err = GretlError.types;
            }
        } else {
            *err = GretlError.types;
        }
    } else {
        ret = get_const_by_name(name, err);
    }

 bailout:

    if (*err) {
        gretl_errmsg_sprintf(_("'%s': not a scalar").ptr, name);
    }

    return ret;
}

int gretl_is_scalar(const(char) * name)
{
    int ret = 0;

    if (get_user_var_of_type_by_name(name, GretlType.double_) != null) {
        ret = 1;
    }

    if (!ret) {
        ret = const_lookup(name);
    }

    return ret;
}

/**
 * get_string_by_name:
 * @name: the name of the string variable to access.
 *
 * Returns: the value of string variable @name, or %NULL
 * if there is no such variable. Note that this is the
 * actual string value, not a copy thereof, compare
 * copy_string_by_name().
 */

char * get_string_by_name(const(char) * name)
{
    UserVar * u = null;

    if (name == null) {
        return null;
    }

    u = get_user_var_of_type_by_name(name, GretlType.string_);
    if (u != null) {
        return cast(char *) u.ptr;
    } else {
        return get_built_in_string_by_name(name);
    }
}

/**
 * copy_string_by_name:
 * @name: the name of the string variable to access.
 * @err: location to receive error code.
 *
 * Returns: a copy of the value of string variable @name,
 * or %NULL on failure.
 */

char * copy_string_by_name(const(char) * name, int * err)
{
    UserVar * u;
    const(char) * s;
    char * ret = null;

    u = get_user_var_of_type_by_name(name, GretlType.string_);

    if (u != null) {
        s = cast(const(char) *) u.ptr;
    } else {
        s = get_built_in_string_by_name(name);
    }

    if (s == null) {
        *err = GretlError.data;
    } else {
        ret = gretl_strdup(s);
        if (ret == null) {
            *err = GretlError.alloc;
        }
    }

    return ret;
}

/**
 * gretl_is_string:
 * @name: name to test.
 *
 * Returns: 1 if @name is the name of a currently defined
 * string variable, otherwise 0.
 */

int gretl_is_string(const(char) * name)
{
    if (*name == '@' && *(name + 1) != '@') {
        name++;
    }

    if (get_user_var_of_type_by_name(name, GretlType.string_) != null) {
        return 1;
    } else if (get_built_in_string_by_name(name) != null) {
        return 1;
    } else {
        return 0;
    }
}

int is_user_string(const(char) * name)
{
    if (*name == '@' && *(name + 1) != '@') {
        name++;
    }

    if (get_user_var_of_type_by_name(name, GretlType.string_) != null) {
        return 1;
    } else {
        return 0;
    }
}

int max_varno_in_saved_lists()
{
    int * list;
    int i, j, vmax = 0;

    for (i=0; i<n_vars; i++) {
        if (uvars[i].type == GretlType.list) {
            list = cast(typeof(list)) uvars[i].ptr;
            if (list != null) {
                for (j=1; j<=list[0]; j++) {
                    if (list[j] > vmax) {
                        vmax = list[j];
                    }
                }
            }
        }
    }

    return vmax;
}

int var_is_deleted(const(int) * dlist, int dmin, int i)
{
    int v = dmin + i - 1;

    if (dlist != null) {
        return in_gretl_list(dlist, v);
    } else {
        return (v >= dmin);
    }
}

/**
 * gretl_lists_revise:
 * @dlist: list of variables to be deleted (or NULL).
 * @dmin: lowest ID number of deleted var (referenced only
 * if @dlist is NULL).
 *
 * Goes through any saved lists, adjusting the ID numbers
 * they contain to reflect the deletion from the dataset of
 * certain variables: those referenced in @dlist, if given,
 * or if @dlist is NULL, those variables with IDs greater
 * than or equal to @dmin.
 *
 * Returns: 0 on success, non-zero code on failure.
 */

int gretl_lists_revise(const(int) * dlist, int dmin)
{
    int * list, maplist;
    int lmax = 0;
    int i, j, k;

    if (dlist != null) {
        /* determine lowest deleted ID */
        dmin = dlist[1];
        for (i=2; i<=dlist[0]; i++) {
            if (dlist[i] > 0 && dlist[i] < dmin) {
                dmin = dlist[i];
            }
        }
    }

    /* find highest ID ref'd in any saved list */
    for (j=0; j<n_vars; j++) {
        if (uvars[j].type == GretlType.list) {
            list = cast(typeof(list)) uvars[j].ptr;
            if (list != null) {
                for (i=1; i<=list[0]; i++) {
                    if (list[i] > lmax) {
                        lmax = list[i];
                    }
                }
            }
        }
    }

    if (lmax < dmin) {
        /* nothing to be done */
        return 0;
    }

    /* make mapping from old to new IDs */

    maplist = gretl_list_new(lmax - dmin + 1);
    if (maplist == null) {
        return GretlError.alloc;
    }

    j = dmin;

    for (i=1; i<=maplist[0]; i++) {
        if (var_is_deleted(dlist, dmin, i)) {
            maplist[i] = -1;
        } else {
            maplist[i] = j++;
        }
    }

    /* use mapping to revise saved lists */
    for (j=0; j<n_vars; j++) {
        if (uvars[j].type == GretlType.list) {
            list = cast(typeof(list)) uvars[j].ptr;
            if (list != null) {
                for (i=list[0]; i>0; i--) {
                    k = list[i] - dmin + 1;
                    if (k >= 1) {
                        if (maplist[k] == -1) {
                            gretl_list_delete_at_pos(list, i);
                        } else {
                            list[i] = maplist[k];
                        }
                    }
                }
            }
        }
    }

    free(maplist);

    return 0;
}

/**
 * gretl_lists_cleanup:
 *
 * Frees all resources associated with the internal
 * apparatus for saving and retrieving named lists.
 */

void gretl_lists_cleanup()
{
    real_destroy_user_vars_at_level(0,
                                    GretlType.list,
                                    0);
}

/* below: serialization of user vars to XML, plus de-serialization
   from XML -- for use in GUI session mechanism
*/

void write_user_scalars(Prn * prn)
{
    double x;
    int i;

    for (i=0; i<n_vars; i++) {
        if (uvars[i].type == GretlType.double_) {
            x = *cast(double *) uvars[i].ptr;
            pprintf(prn, " <gretl-scalar name=\"%s\" value=\"".ptr, uvars[i].name.ptr);
            serialize_scalar_value(x, prn);
            pputs(prn, "\"/>\n");
        }
    }
}

void write_user_matrices(Prn * prn)
{
    int i;

    for (i=0; i<n_vars; i++) {
        if (uvars[i].type == GretlType.matrix) {
            xml_put_user_matrix(uvars[i], prn);
        }
    }
}

void write_user_lists(Prn * prn)
{
    int i;

    for (i=0; i<n_vars; i++) {
        if (uvars[i].type == GretlType.list) {
            gretl_list_serialize(cast(const(int) *) uvars[i].ptr,
                                 uvars[i].name.ptr,
                                 prn);
        }
    }
}

void write_user_bundles(Prn * prn)
{
    int i, err = 0;

    for (i=0; i<n_vars && !err; i++) {
        if (uvars[i].type == GretlType.bundle) {
            err = gretl_bundle_serialize(cast(GretlBundle *) uvars[i].ptr,
                     uvars[i].name.ptr,
                     prn);
        }
    }

    // return err; NOT YET
}

int read_user_scalars(xmlDocPtr doc, xmlNodePtr cur)
{
    char * name, val;
    double x;
    int n, err = 0;

    cur = cur.children;

    gretl_push_c_numeric_locale();

    while (cur != null && !err) {
        if (!xmlStrcmp(cur.name, cast(XUC) "gretl-scalar")) {
            name = cast(char *) xmlGetProp(cur, cast(XUC) "name");
            val = cast(char *) xmlGetProp(cur, cast(XUC) "value");
            if (name == null || val == null) {
                err = 1;
            } else {
                n = sscanf(val, "%lf".ptr, &x);
                if (n < 1) {
                    x = NADBL;
                }
                err = gretl_scalar_add(name, x);
            }
            free(name);
            free(val);
        }
        cur = cur.next;
    }

    gretl_pop_c_numeric_locale();

    return err;
}

int read_user_matrices(xmlDocPtr doc, xmlNodePtr cur)
{
    GretlMatrixRaw * m;
    char * name;
    int err = 0;

    cur = cur.children;

    while (cur != null && !err) {
        if (!xmlStrcmp(cur.name, cast(XUC) "gretl-matrix")) {
            name = cast(char *) xmlGetProp(cur, cast(XUC) "name");
            if (name == null) {
                err = 1;
            } else {
                m = gretl_xml_get_matrix(cur, doc, &err);
                if (m != null) {
                    err = user_var_add(name, GretlType.matrix, m);
                }
                free(name);
            }
        }
        cur = cur.next;
    }

    return err;
}

int read_user_lists(xmlDocPtr doc, xmlNodePtr cur)
{
    int * list;
    char * name;
    int err = 0;

    cur = cur.children;

    while (cur != null && !err) {
        if (!xmlStrcmp(cur.name, cast(XUC) "list")) {
            if (!gretl_xml_get_prop_as_string(cur, "name", &name)) {
                err = GretlError.data;
            } else {
                list = gretl_xml_get_list(cur, doc, &err);
                if (!err) {
                    err = user_var_add(name, GretlType.list, list);
                }
                free(name);
            }
        }
        cur = cur.next;
    }

    return err;
}

int read_user_bundles(xmlDocPtr doc, xmlNodePtr cur)
{
    int err = 0;

    gretl_push_c_numeric_locale();

    cur = cur.children;

    while (cur != null && !err) {
        if (!xmlStrcmp(cur.name, cast(XUC) "gretl-bundle")) {
            char * name = cast(char *) xmlGetProp(cur, cast(XUC) "name");

            if (name == null) {
                err = 1;
            } else {
                char * creator = null;
                GretlBundle * b;

                b = gretl_bundle_deserialize(cur, doc, &err);
                if (!err) {
                    creator = cast(char *) xmlGetProp(cur, cast(XUC) "creator");
                    gretl_bundle_set_creator(b, creator);
                    err = user_var_add(name, GretlType.bundle, b);
                }
                free(name);
                free(creator);
            }
        }
        cur = cur.next;
    }

    gretl_pop_c_numeric_locale();

    return err;
}

alias XmlWriteFunc = void function(Prn *);
alias XmlReadFunc = int function(xmlDocPtr, xmlNodePtr);

struct UvarFile {
    GretlType type;
    const(char) * typestr;
    XmlWriteFunc write_func;
    XmlReadFunc read_func;
}

__gshared UvarFile[4] uvar_files = [

    UvarFile(GretlType.double_, "scalars".ptr, &write_user_scalars, &read_user_scalars),
    UvarFile(GretlType.matrix, "matrices".ptr, &write_user_matrices, &read_user_matrices),
    UvarFile(GretlType.list, "lists".ptr, &write_user_lists, &read_user_lists),
    UvarFile(GretlType.bundle, "bundles".ptr, &write_user_bundles, &read_user_bundles)
];

int serialize_user_vars(const(char) * dirname)
{
    GretlType type;
    const(char) * typestr;
    alias WriteFuncFn = void function(Prn *);
    WriteFuncFn write_func;
    char[MAXLEN] path;
    Prn * prn;
    int i, n, ni;
    int err = 0;

    n = cast(int) uvar_files.length;

    gretl_push_c_numeric_locale();

    for (i=0; i<n; i++) {
        type = uvar_files[i].type;
        ni = user_var_count_for_type(type);
        if (ni > 0) {
            int errp = 0;

            typestr = uvar_files[i].typestr;
            sprintf(path.ptr, "%s%c%s.xml".ptr, dirname, SLASH, typestr);
            write_func = uvar_files[i].write_func;
            prn = gretl_print_new_with_filename(path.ptr, &errp);
            if (prn == null) {
                err++;
                continue;
            }
            gretl_xml_header(prn);
            pprintf(prn, "<gretl-%s count=\"%d\">\n".ptr, typestr, ni);
            write_func(prn);
            pprintf(prn, "</gretl-%s>\n".ptr, typestr);
            gretl_print_destroy(prn);
        }
    }

    gretl_pop_c_numeric_locale();

    if (err > 0) {
        fprintf(stderr, "Failed writing %d UserVar files\n".ptr, err);
        err = GretlError.fopen;
    }

    return err;
}

int deserialize_user_vars(const(char) * dirname)
{
    xmlDocPtr doc = null;
    xmlNodePtr cur = null;
    const(char) * typestr;
    alias ReadFuncFn = int function(xmlDocPtr, xmlNodePtr);
    ReadFuncFn read_func;
    char[16] root_name;
    char[MAXLEN] path;
    FILE * fp;
    int i, n;
    int n_failed = 0;
    int err = 0;

    n = cast(int) uvar_files.length;

    for (i=0; i<n; i++) {
        int err_i = 0;

        typestr = uvar_files[i].typestr;
        sprintf(path.ptr, "%s%c%s.xml".ptr, dirname, SLASH, typestr);

        fp = gretl_fopen(path.ptr, "r");
        if (fp == null) {
            /* OK, no user-vars of this type */
            continue;
        }
        fclose(fp);
        sprintf(root_name.ptr, "gretl-%s".ptr, typestr);
        err_i = gretl_xml_open_doc_root(path.ptr, root_name.ptr, &doc, &cur);
        if (!err_i) {
            read_func = uvar_files[i].read_func;
            err_i = read_func(doc, cur);
        }
        if (doc != null) {
            xmlFreeDoc(doc);
            doc = null;
        }
        if (err_i) {
            n_failed++;
            if (!err) {
                err = err_i;
            }
        }
    }

    if (n_failed > 0) {
        fprintf(stderr, "Failed reading %d UserVar files\n".ptr, n_failed);
    }

    return err;
}

int print_user_var_by_name(const(char) * name, const(DatasetRaw) * dset, Gretlopt opt, Prn * prn)
{
    UserVar * u = get_user_var_by_name(name);
    int err = 0;

    if (u == null || u.ptr == null) {
        return GretlError.data;
    }

    if (u.type == GretlType.double_) {
        print_scalar_by_name(name, prn);
    } else if (u.type == GretlType.matrix) {
        GretlMatrixRaw * tmp = cast(GretlMatrixRaw *) u.ptr;

        if (tmp.is_complex || opt & Gretlopt.OPT_C) {
            err = gretl_cmatrix_print(tmp, name, prn);
    } else {
            gretl_matrix_print_to_prn(tmp, name, prn);
        }
    } else if (u.type == GretlType.bundle) {
        if (opt & Gretlopt.OPT_Y) {
            gretl_bundle_print_tree(cast(GretlBundle *) u.ptr, prn);
        } else {
            gretl_bundle_print(cast(GretlBundle *) u.ptr, prn);
        }
    } else if (u.type == GretlType.array) {
        gretl_array_print(cast(GretlArray *) u.ptr, prn);
    } else if (u.type == GretlType.list) {
        gretl_list_print(cast(const(int) *) u.ptr, dset, prn);
    } else if (u.type == GretlType.string_) {
        pputs(prn, cast(char *) u.ptr);
        pputc(prn, '\n');
    }

    return err;
}

int uvar_type_match(UserVar * u, GretlType t)
{
    if (u.type == t) {
        return 1;
    } else if (u.type == GretlType.array &&
               gretl_array_type(t)) {
        return t == gretl_array_get_type(cast(GretlArray *) u.ptr);
    } else {
        return 0;
    }
}

int list_user_vars_of_type(const(DatasetRaw) * dset, Prn * prn)
{
    const(char) * typename;
    GretlType t;

    typename = get_optval_string(GretlCmdIndex.varlist, Gretlopt.OPT_T);
    if (typename == null) {
        return GretlError.invarg;
    }

    if (!strcmp(typename, "accessor")) {
        list_ok_dollar_vars(cast(DatasetRaw *) dset, prn);
        return 0;
    }

    t = gretl_type_from_string(typename);
    if (t == GretlType.none) {
        return GretlError.invarg;
    }

    if (t == GretlType.series) {
        list_series(dset, Gretlopt.OPT_NONE, prn);
    } else if (t == GretlType.double_) {
        print_scalars(prn);
    } else if (t == GretlType.list ||
               t == GretlType.matrix ||
               t == GretlType.bundle ||
               t == GretlType.array ||
               t == GretlType.string_ ||
               gretl_array_type(t)) {
        int i, n = 0;

        pprintf(prn, _("variables of type %s:").ptr, typename);
        for (i=0; i<n_vars; i++) {
            if (uvar_type_match(uvars[i], t)) {
                if (n == 0) {
                    pputc(prn, '\n');
                }
                if (uvars[i].name[0] == '\0') {
                    pputs(prn, _("  (unnamed)\n").ptr);
                } else if (t == GretlType.array) {
                    GretlType at = gretl_array_get_type(cast(GretlArray *) uvars[i].ptr);

                    pprintf(prn, "  %s (%s)\n".ptr, uvars[i].name.ptr,
                            gretl_type_get_name(at));
                } else {
                    pprintf(prn, "  %s\n".ptr, uvars[i].name.ptr);
                }
                n++;
            }
        }
        if (n == 0) {
            pprintf(prn, " %s\n".ptr, _("none").ptr);
        }
        pputc(prn, '\n');
    } else {
        return GretlError.invarg;
    }

    return 0;
}

int leads_midas_list(int ID, const(DatasetRaw) * dset, char * listname)
{
    int level = gretl_function_depth();
    int * list;
    int i, ret = 0;

    for (i=0; i<n_vars && !ret; i++) {
        if (uvars[i].type == GretlType.list &&
            uvars[i].level == level) {
            list = cast(typeof(list)) uvars[i].ptr;
            if (list[0] > 2 && list[1] == ID) {
                ret = gretl_is_midas_list(list, dset);
                if (ret && listname != null) {
                    strcpy(listname, uvars[i].name.ptr);
                }
            }
        }
    }

    return ret;
}

int in_midas_list(int ID, const(DatasetRaw) * dset, char * listname)
{
    int level = gretl_function_depth();
    int * list;
    int i, ret = 0;

    for (i=0; i<n_vars && !ret; i++) {
        if (uvars[i].type == GretlType.list &&
            uvars[i].level == level) {
            list = cast(typeof(list)) uvars[i].ptr;
            if (list[0] > 2 && in_gretl_list(list, ID)) {
                ret = gretl_is_midas_list(list, dset);
                if (ret && listname != null) {
                    strcpy(listname, uvars[i].name.ptr);
                }
            }
        }
    }

    return ret;
}

const(char) * get_listname_by_consecutive_content(int l0, int l1)
{
    int level = gretl_function_depth();
    const(char) * ret = null;
    int i;
    int j;
    int * list;

    for (i=0; i<n_vars; i++) {
        if (uvars[i].type == GretlType.list &&
            uvars[i].level == level) {
            list = cast(typeof(list)) uvars[i].ptr;
            if (list[0] == l0 && list[1] == l1) {
                int found = 1;

                for (j=2; j<=l0; j++) {
                    if (list[j] != list[j-1] + 1) {
                        found = 0;
                        break;
                    }
                }
                if (found) {
                    return uvars[i].name.ptr;
                }
            }
        }
    }

    return ret;
}

/* Dropping terms from the list @targ inside a function:
   this is tricky with regard to the auto-generated
   "time" variable, which may be a member of a list that
   was passed as an argument yet not "visible" by name
   within the function. Here we attempt to fix this by
   transcribing the ID number for "time" from the caller's
   namespace into the @drop list -- if the latter is
   trying to drop this variable.
*/

void check_auto_time_var(const(int) * targ, int * drop, const(DatasetRaw) * dset)
{
    int i, vi, tnum = 0;

    for (i=1; i<=targ[0]; i++) {
        vi = targ[i];
        if (!strcmp(dset.varname[vi], "time")) {
            tnum = vi;
            break;
        }
    }

    if (tnum > 0) {
        for (i=drop[0]; i>0; i--) {
            vi = drop[i];
            if (!strcmp(dset.varname[vi], "time")) {
                drop[i] = tnum;
            }
        }
    }
}

/* functions called from geneval.c when "editing" a named list */

int user_list_append(UserVar * uvar, const(int) * add)
{
    int err = 0;

    if (uvar == null || user_var_get_type(uvar) != GretlType.list) {
        err = GretlError.data;
    } else {
        const(int) * list = cast(const(int) *) user_var_get_value(uvar);
        int * tmp = gretl_list_copy(list);

        if (tmp == null) {
            err = GretlError.alloc;
        } else {
            err = gretl_list_append_list(&tmp, add);
            if (!err) {
                user_var_replace_value(uvar, tmp, GretlType.list);
            }
        }
    }

    return err;
}

int user_list_subtract(UserVar * uvar, int * sub, const(DatasetRaw) * dset)
{
    int err = 0;

    if (uvar == null || user_var_get_type(uvar) != GretlType.list) {
        err = GretlError.data;
    } else {
        const(int) * list = cast(const(int) *) user_var_get_value(uvar);
        int * tmp;

        if (gretl_function_depth() > 0) {
            check_auto_time_var(list, sub, dset);
        }
        tmp = gretl_list_drop(list, sub, &err);
        if (!err) {
            user_var_replace_value(uvar, tmp, GretlType.list);
        }
    }

    return err;
}

int user_list_replace(UserVar * uvar, const(int) * src)
{
    int err = 0;

    if (uvar == null || user_var_get_type(uvar) != GretlType.list) {
        err = GretlError.data;
    } else {
        int * tmp = gretl_list_copy(src);

        if (tmp == null) {
            err = GretlError.alloc;
        } else {
            user_var_replace_value(uvar, tmp, GretlType.list);
        }
    }

    return err;
}

/**
 * remember_list:
 * @list: array of integers, the first element being a count
 * of the following elements.
 * @name: name to be given to the list.
 * @prn: printing struct.
 *
 * Adds a copy of @list to the stack of saved lists and associates
 * it with @name, unless there is already a list with the given
 * name in which case the original list is replaced.  A status
 * message is printed to @prn.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int remember_list(const(int) * list, const(char) * name, Prn * prn)
{
    int * lcpy = gretl_list_copy(list);
    int err = 0;

    if (lcpy == null) {
        err = (list == null)? GretlError.data : GretlError.alloc;
    } else {
        UserVar * orig;

        orig = get_user_var_of_type_by_name(name, GretlType.list);

        if (orig != null) {
            /* replace existing list of same name */
            user_var_replace_value(orig, lcpy, GretlType.list);
            if (prn != null && gretl_messages_on()) {
                pprintf(prn, _("Replaced list '%s'\n").ptr, name);
            }
        } else {
            err = user_var_add(name, GretlType.list, lcpy);
            if (!err && prn != null && gretl_messages_on()) {
                pprintf(prn, _("Added list '%s'\n").ptr, name);
            }
        }
    }

    return err;
}

/**
 * get_list_by_name:
 * @name: the name of the list to be found.
 *
 * Looks up @name in the stack of saved variables, at the current level
 * of function execution, and retrieves the associated list.
 *
 * Returns: the list, or NULL if the lookup fails.
 */

int * get_list_by_name(const(char) * name)
{
    UserVar * u;
    int * ret = null;

    u = get_user_var_of_type_by_name(name, GretlType.list);

    if (u != null) {
        ret = cast(int *) user_var_get_value(u);
    }

    return ret;
}

/**
 * gretl_is_list:
 * @name: the name to test.
 *
 * Returns: 1 if @name is the name of a saved list, 0
 * otherwise.
 */

int gretl_is_list(const(char) * name)
{
    return get_user_var_of_type_by_name(name, GretlType.list) != null;
}

