/* Port of lib/src/gretl_typemap.c: relations between the gretl data types: array (plural) and singular types,
   reference ("pointer") and plain types, type names and their parsing, the "gen" type keywords, the ordering of
   types used by typeof(), scalar/series/array predicates and the name-conflict message. The header
   gretl_typemap.h is gretl_typemaph.d (the NUMERIC_TYPE macro); its prototypes belong to this module.

   Differences from the C: the "gen types" lookup table (a glib hash table) is the string-keyed table of
   ghashtable.d, and its int values travel as pointers as in the C; the translation-marked messages use the
   gettext stand-in. */

module gretl.gretl_typemap;


import core.stdc.string : strcmp, strncmp;

import gretl.ghashtable;
import gretl.glib : gchar, gpointer;
import gretl.gstrfuncs : g_strdup_printf;
import gretl.libgretlh;

private struct TypeMap {
    GretlType std;     /* a "standard" type: non-array, non-reference */
    GretlType stdref;  /* the "reference" form of the standard type, if any */
    GretlType plural;  /* the associated array type, if any */
    GretlType plref;   /* reference form of array type, if any */
}

private immutable TypeMap[7] gretl_type_map = [
    TypeMap(GretlType.matrix, GretlType.matrix_ref, GretlType.matrices, GretlType.matrices_ref),
    TypeMap(GretlType.bundle, GretlType.bundle_ref, GretlType.bundles, GretlType.bundles_ref),
    TypeMap(GretlType.string_, GretlType.string_ref, GretlType.strings, GretlType.strings_ref),
    TypeMap(GretlType.array, GretlType.array_ref, GretlType.arrays, GretlType.arrays_ref),
    TypeMap(GretlType.list, GretlType.list_ref, GretlType.lists, GretlType.lists_ref),
    TypeMap(GretlType.series, GretlType.series_ref, GretlType.none, GretlType.none),
    TypeMap(GretlType.double_, GretlType.scalar_ref, GretlType.none, GretlType.none)
];

private enum N_TYPES = gretl_type_map.length;

GretlType gretl_type_get_plural(GretlType type)
{
    if (type == 0) return GretlType.none;

    foreach (i; 0 .. N_TYPES) {
        if (type == gretl_type_map[i].std) {
            return gretl_type_map[i].plural;
        }
    }

    return GretlType.none;
}

GretlType gretl_type_get_singular(GretlType type)
{
    if (type == 0) return GretlType.none;

    foreach (i; 0 .. N_TYPES) {
        if (type == gretl_type_map[i].plural) {
            return gretl_type_map[i].std;
        }
    }

    return GretlType.none;
}

GretlType gretl_type_get_ref_type(GretlType type)
{
    if (type == 0) return GretlType.none;

    foreach (i; 0 .. N_TYPES) {
        if (type == gretl_type_map[i].std) {
            return gretl_type_map[i].stdref;
        }
        if (type == gretl_type_map[i].plural) {
            return gretl_type_map[i].plref;
        }
        if (type == gretl_type_map[i].stdref || type == gretl_type_map[i].plref) {
            /* allow pass-through */
            return type;
        }
    }

    return GretlType.none;
}

GretlType gretl_type_get_plain_type(GretlType type)
{
    if (type == 0) return GretlType.none;

    foreach (i; 0 .. N_TYPES) {
        if (type == gretl_type_map[i].stdref) {
            return gretl_type_map[i].std;
        }
        if (type == gretl_type_map[i].plref) {
            return gretl_type_map[i].plural;
        }
        if (type == gretl_type_map[i].std || type == gretl_type_map[i].plural) {
            /* allow pass-through */
            return type;
        }
    }

    return GretlType.none;
}

int gretl_is_array_type(GretlType type)
{
    if (type == 0) return 0;

    foreach (i; 0 .. N_TYPES) {
        if (type == gretl_type_map[i].plural) {
            return 1;
        }
    }

    return 0;
}

int gretl_is_array_ref_type(GretlType type)
{
    if (type == 0) return 0;

    foreach (i; 0 .. N_TYPES) {
        if (type == gretl_type_map[i].plref) {
            return 1;
        }
    }

    return 0;
}

int gretl_is_arrayable_type(GretlType type)
{
    return gretl_type_get_plural(type) > 0;
}

/**
 * gretl_type_get_name:
 * @type: a gretl type.
 *
 * Returns: the name of @type, or "invalid" on failure.
 */

const(char) * gretl_type_get_name(GretlType type)
{
    switch (type) {
    case GretlType.bool_:       return "bool";
    case GretlType.int_:        return "int";
    case GretlType.uint32:      return "unsigned32";
    case GretlType.uint64:      return "unsigned64";
    case GretlType.obs:         return "obs";
    case GretlType.double_:     return "scalar";
    case GretlType.series:      return "series";
    case GretlType.useries:     return "series";
    case GretlType.matrix:      return "matrix";
    case GretlType.list:        return "list";
    case GretlType.bundle:      return "bundle";
    case GretlType.array:       return "array";
    case GretlType.string_:     return "string";
    case GretlType.scalar_ref:  return "scalar *";
    case GretlType.series_ref:  return "series *";
    case GretlType.matrix_ref:  return "matrix *";
    case GretlType.bundle_ref:  return "bundle *";
    case GretlType.array_ref:   return "array *";
    case GretlType.string_ref:  return "string *";
    case GretlType.list_ref:    return "list *";

    case GretlType.strings:      return "strings";
    case GretlType.matrices:     return "matrices";
    case GretlType.bundles:      return "bundles";
    case GretlType.lists:        return "lists";
    case GretlType.arrays:       return "arrays";

    case GretlType.strings_ref:  return "strings *";
    case GretlType.matrices_ref: return "matrices *";
    case GretlType.bundles_ref:  return "bundles *";
    case GretlType.lists_ref:    return "lists *";
    case GretlType.arrays_ref:   return "arrays *";

    case GretlType.date:         return "date"; /* ODBC special */

    case GretlType.void_:       return "void";
    case GretlType.numeric:     return "numeric";
    case GretlType.none:        return "null";
    case GretlType.any:         return "any";
    default:
        return "invalid";
    }
}

private bool is_ref_suffix(const(char) * p)
{
    return !strcmp(p, " *") || !strcmp(p, "ref");
}

GretlType gretl_type_from_string(const(char) * s)
{
    const(char) * p;

    if (!strncmp(s, "matrix", 6)) {
        p = s + 6;
        if (*p == '\0') {
            return GretlType.matrix;
        } else if (is_ref_suffix(p)) {
            return GretlType.matrix_ref;
        }
    } else if (!strncmp(s, "series", 6)) {
        p = s + 6;
        if (*p == '\0') {
            return GretlType.series;
        } else if (is_ref_suffix(p)) {
            return GretlType.series_ref;
        }
    } else if (!strncmp(s, "scalar", 6)) {
        p = s + 6;
        if (*p == '\0') {
            return GretlType.double_;
        } else if (is_ref_suffix(p)) {
            return GretlType.scalar_ref;
        }
    } else if (!strncmp(s, "bundle", 6)) {
        p = s + 6;
        if (*p == 's') {
            p++;
            if (*p == '\0') {
                return GretlType.bundles;
            } else if (is_ref_suffix(p)) {
                return GretlType.bundles_ref;
            }
        } else if (*p == '\0') {
            return GretlType.bundle;
        } else if (is_ref_suffix(p)) {
            return GretlType.bundle_ref;
        }
    } else if (!strncmp(s, "array", 5)) {
        p = s + 5;
        if (*p == 's') {
            p++;
            if (*p == '\0') {
                return GretlType.arrays;
            } else if (is_ref_suffix(p)) {
                return GretlType.arrays_ref;
            }
        } else if (*p == '\0') {
            return GretlType.array;
        } else if (is_ref_suffix(p)) {
            return GretlType.array_ref;
        }
    } else if (!strncmp(s, "string", 6)) {
        p = s + 6;
        if (*p == 's') {
            p++;
            if (*p == '\0') {
                return GretlType.strings;
            } else if (is_ref_suffix(p)) {
                return GretlType.strings_ref;
            }
        } else if (*p == '\0') {
            return GretlType.string_;
        } else if (is_ref_suffix(p)) {
            return GretlType.string_ref;
        }
    } else if (!strncmp(s, "matrices", 8)) {
        p = s + 8;
        if (*p == '\0') {
            return GretlType.matrices;
        } else if (is_ref_suffix(p)) {
            return GretlType.matrices_ref;
        }
    } else if (!strncmp(s, "list", 4)) {
        p = s + 4;
        if (*p == 's') {
            p++;
            if (*p == '\0') {
                return GretlType.lists;
            } else if (is_ref_suffix(p)) {
                return GretlType.lists_ref;
            }
        } else if (*p == '\0') {
            return GretlType.list;
        } else if (is_ref_suffix(p)) {
            return GretlType.list_ref;
        }
    } else if (!strncmp(s, "numeric", 7)) {
        return GretlType.numeric;
    } else {
        /* aliases */
        if (!strcmp(s, "bool"))       return GretlType.bool_;
        if (!strcmp(s, "boolean"))    return GretlType.bool_;
        if (!strcmp(s, "int"))        return GretlType.int_;
        if (!strcmp(s, "unsigned32")) return GretlType.uint32;
        if (!strcmp(s, "unsigned64")) return GretlType.uint64;
        if (!strcmp(s, "obs"))        return GretlType.obs;
    }

    return GretlType.none;
}

private struct Lookup {
    const(char) * word;
    GretlType type;
}

private __gshared const(Lookup)[11] gentypes = [
    Lookup("series".ptr,   GretlType.series),
    Lookup("scalar".ptr,   GretlType.double_),
    Lookup("matrix".ptr,   GretlType.matrix),
    Lookup("string".ptr,   GretlType.string_),
    Lookup("list".ptr,     GretlType.list),
    Lookup("bundle".ptr,   GretlType.bundle),
    Lookup("strings".ptr,  GretlType.strings),
    Lookup("matrices".ptr, GretlType.matrices),
    Lookup("bundles".ptr,  GretlType.bundles),
    Lookup("lists".ptr,    GretlType.lists),
    Lookup("arrays".ptr,   GretlType.arrays)
];

GretlType gretl_get_gen_type(const(char) * s)
{
    __gshared GHashTable * ht = null;
    gpointer p;
    GretlType t = GretlType.none;

    if (s == null) {
        /* the clean-up signal */
        if (ht != null) {
            g_hash_table_destroy(ht);
            ht = null;
        }
        return GretlType.none;
    }

    if (ht == null) {
        ht = g_hash_table_new(&g_str_hash, &g_str_equal);
        foreach (i; 0 .. gentypes.length) {
            g_hash_table_insert(ht, cast(gpointer) gentypes[i].word, cast(gpointer) cast(size_t) gentypes[i].type);
        }
    }

    p = g_hash_table_lookup(ht, s);
    if (p != null) {
        t = cast(GretlType) cast(int) cast(size_t) p;
    }

    return t;
}

/* Note: this must agree with the doc for the typeof() function */

int gretl_type_get_order(GretlType type)
{
    if (gretl_scalar_type(type)) {
        return 1;
    } else if (type == GretlType.series) {
        return 2;
    } else if (type == GretlType.matrix) {
        return 3;
    } else if (type == GretlType.string_) {
        return 4;
    } else if (type == GretlType.bundle) {
        return 5;
    } else if (type == GretlType.array) {
        return 6;
    } else if (type == GretlType.list) {
        return 7;
    } else {
        return 0;
    }
}

int gretl_is_scalar_type(GretlType type)
{
    return type == GretlType.bool_ ||
        type == GretlType.int_ ||
        type == GretlType.uint32 ||
        type == GretlType.uint64 ||
        type == GretlType.double_;
}

int gretl_is_series_type(GretlType type)
{
    return type == GretlType.series ||
        type == GretlType.useries;
}

void gretl_typemap_cleanup()
{
    gretl_get_gen_type(null);
}

gchar * name_conflict_message(const(char) * name, GretlType type)
{
    gchar * s = null;

    if (type == GretlType.series) {
        s = g_strdup_printf(_("A series named %s already exists").ptr, name);
    } else if (type == GretlType.matrix) {
        s = g_strdup_printf(_("A matrix named %s already exists").ptr, name);
    } else if (type == GretlType.double_) {
        s = g_strdup_printf(_("A scalar named %s already exists").ptr, name);
    } else if (type == GretlType.list) {
        s = g_strdup_printf(_("A list named %s already exists").ptr, name);
    } else if (type == GretlType.string_) {
        s = g_strdup_printf(_("A string named %s already exists").ptr, name);
    } else if (type == GretlType.bundle) {
        s = g_strdup_printf(_("A bundle named %s already exists").ptr, name);
    } else if (type == GretlType.array) {
        s = g_strdup_printf(_("An array named %s already exists").ptr, name);
    }

    return s;
}

/* Used in gui/fncall.c, where (until it becomes time to initiate an
   actual function call) we don't make a hard distinction between
   plain and pointer types.
*/

int gretl_type_mismatch(GretlType t1, GretlType t2)
{
    if (t1 == t2) {
        return 0;
    } else {
        t1 = gretl_type_get_plain_type(t1);
        t2 = gretl_type_get_plain_type(t2);
        if (t1 == t2) {
            return 0;
        }
    }

    return 1;
}
