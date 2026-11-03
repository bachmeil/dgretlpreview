/*
   Replacement for the glib GHashTable, backed by a D associative array.

   Why this file exists: the C code uses glib's GHashTable. Rather than bind
   to glib, this module provides functions with the glib names and
   signatures, implemented with a built-in D associative array (AA), so
   that the ported call sites read like the C.

   Limitations (deliberate, to keep this small):

   - Only string-keyed tables are supported, i.e. g_hash_table_new() must be
     given g_str_hash and g_str_equal. Anything else fails an assert. The
     keys are NUL-terminated C strings.
   - Only the functions needed so far are provided (new, insert, lookup,
     destroy). Add others here (remove, foreach, size, ...) as ported files
     need them.
   - Keys are NOT copied, as with glib: a key must outlive the table.
     (libset.d uses string literals, which do.)
   - The table is GC-allocated, so this is a GC dependency. If a betterc
     build is ever wanted, replace the body of this module with a
     hand-written hash table, keeping the same function signatures.
   - A table shared between threads must be declared __gshared by its user,
     and is not synchronised, as with glib.

   Extended for gretl_bundle.d: g_hash_table_new_full (key/value destroy notifies, called as glib calls them),
   replace, remove, steal, size, foreach, get_keys, get_values (the latter two build GLists, see glist.d) and the
   GHashTableIter functions (a snapshot of the keys is iterated). The AA is keyed by a copy of the key string,
   so the caller's key memory may be freed by key_destroy without affecting the table; the original key pointer
   is kept in the entry and is what get_keys/foreach hand back. Iteration order is that of the D AA
   (unspecified, as in glib).

   Used by: libset.d, gretl_bundle.d.
*/

module gretl.ghashtable;


import core.stdc.string : strcmp;
import std.string : fromStringz;

import gretl.glib : gboolean, gconstpointer, gpointer, guint;
import gretl.glist : GList, g_list_append;

alias GHashFunc = guint function(gconstpointer key);
alias GEqualFunc = gboolean function(gconstpointer a, gconstpointer b);
alias GDestroyNotify = void function(gpointer data);
alias GHFunc = void function(gpointer key, gpointer value, gpointer user_data);

private struct HEntry {
    gpointer key;
    gpointer value;
}

struct GHashTable {
    HEntry[string] map;
    GDestroyNotify key_destroy;
    GDestroyNotify value_destroy;
}

struct GHashTableIter {
    GHashTable * ht;
    string[] keys;
    size_t pos;
}

/* Tags for g_hash_table_new(); the AA does the hashing and comparison. */

guint g_str_hash(gconstpointer v)
{
    /* not used by the AA; present so call sites stay as in the C */
    guint h = 5381;
    for (const(ubyte) * p = cast(const(ubyte) *) v; *p != 0; p++) {
        h = h * 33 + *p;
    }
    return h;
}

gboolean g_str_equal(gconstpointer v1, gconstpointer v2)
{
    return strcmp(cast(const(char) *) v1, cast(const(char) *) v2) == 0;
}

GHashTable * g_hash_table_new_full(GHashFunc hash_func, GEqualFunc key_equal_func,
                                   GDestroyNotify key_destroy_func, GDestroyNotify value_destroy_func)
{
    assert(hash_func == &g_str_hash && key_equal_func == &g_str_equal,
           "ghashtable.d supports only string keys");

    GHashTable * ht = new GHashTable;

    ht.key_destroy = key_destroy_func;
    ht.value_destroy = value_destroy_func;
    return ht;
}

GHashTable * g_hash_table_new(GHashFunc hash_func, GEqualFunc key_equal_func)
{
    return g_hash_table_new_full(hash_func, key_equal_func, null, null);
}

private string keystr(gconstpointer key)
{
    return fromStringz(cast(const(char) *) key).idup;
}

/* glib returns TRUE if the key is new. For an existing key the old key is kept and the passed key is freed (with
   key_destroy), the old value is freed and replaced. */

gboolean g_hash_table_insert(GHashTable * hash_table, gpointer key, gpointer value)
{
    string k = keystr(key);
    HEntry * e = k in hash_table.map;

    if (e is null) {
        hash_table.map[k] = HEntry(key, value);
        return true;
    }
    if (hash_table.key_destroy !is null && e.key != key) {
        hash_table.key_destroy(key);
    }
    if (hash_table.value_destroy !is null && e.value != value) {
        hash_table.value_destroy(e.value);
    }
    e.value = value;
    return false;
}

/* like insert, but the old key (not the new one) is freed */

gboolean g_hash_table_replace(GHashTable * hash_table, gpointer key, gpointer value)
{
    string k = keystr(key);
    HEntry * e = k in hash_table.map;

    if (e is null) {
        hash_table.map[k] = HEntry(key, value);
        return true;
    }
    if (hash_table.key_destroy !is null && e.key != key) {
        hash_table.key_destroy(e.key);
    }
    if (hash_table.value_destroy !is null && e.value != value) {
        hash_table.value_destroy(e.value);
    }
    e.key = key;
    e.value = value;
    return false;
}

/* returns null if the key is not found */

gpointer g_hash_table_lookup(GHashTable * hash_table, gconstpointer key)
{
    HEntry * e = keystr(key) in hash_table.map;

    return e !is null ? e.value : null;
}

gboolean g_hash_table_remove(GHashTable * hash_table, gconstpointer key)
{
    string k = keystr(key);
    HEntry * e = k in hash_table.map;

    if (e is null) {
        return false;
    }
    HEntry old = *e;

    hash_table.map.remove(k);
    if (hash_table.key_destroy !is null) {
        hash_table.key_destroy(old.key);
    }
    if (hash_table.value_destroy !is null) {
        hash_table.value_destroy(old.value);
    }
    return true;
}

/* removes the entry without calling the destroy notifies */

gboolean g_hash_table_steal(GHashTable * hash_table, gconstpointer key)
{
    string k = keystr(key);

    if ((k in hash_table.map) is null) {
        return false;
    }
    hash_table.map.remove(k);
    return true;
}

guint g_hash_table_size(GHashTable * hash_table)
{
    return cast(guint) hash_table.map.length;
}

void g_hash_table_foreach(GHashTable * hash_table, GHFunc func, gpointer user_data)
{
    foreach (k; hash_table.map.keys) {
        HEntry * e = k in hash_table.map;

        if (e !is null) {
            func(e.key, e.value, user_data);
        }
    }
}

GList * g_hash_table_get_keys(GHashTable * hash_table)
{
    GList * l = null;

    foreach (k, e; hash_table.map) {
        l = g_list_append(l, e.key);
    }
    return l;
}

GList * g_hash_table_get_values(GHashTable * hash_table)
{
    GList * l = null;

    foreach (k, e; hash_table.map) {
        l = g_list_append(l, e.value);
    }
    return l;
}

void g_hash_table_iter_init(GHashTableIter * iter, GHashTable * hash_table)
{
    iter.ht = hash_table;
    iter.keys = hash_table.map.keys;
    iter.pos = 0;
}

gboolean g_hash_table_iter_next(GHashTableIter * iter, gpointer * key, gpointer * value)
{
    while (iter.pos < iter.keys.length) {
        HEntry * e = iter.keys[iter.pos++] in iter.ht.map;

        if (e !is null) {
            if (key !is null) {
                *key = e.key;
            }
            if (value !is null) {
                *value = e.value;
            }
            return true;
        }
    }
    return false;
}

/* removes every entry, calling the destroy notifies */

void g_hash_table_remove_all(GHashTable * hash_table)
{
    g_hash_table_destroy(hash_table);
}

/* calls the destroy notifies on all entries, then empties the table */

void g_hash_table_destroy(GHashTable * hash_table)
{
    foreach (k, e; hash_table.map) {
        if (hash_table.key_destroy !is null) {
            hash_table.key_destroy(e.key);
        }
        if (hash_table.value_destroy !is null) {
            hash_table.value_destroy(e.value);
        }
    }
    hash_table.map = null;
}

unittest
{
    GHashTable * ht = g_hash_table_new(&g_str_hash, &g_str_equal);
    int a, b;

    assert(g_hash_table_insert(ht, cast(gpointer) "alpha".ptr, &a));
    assert(g_hash_table_insert(ht, cast(gpointer) "beta".ptr, &b));
    assert(!g_hash_table_insert(ht, cast(gpointer) "beta".ptr, &a));
    assert(g_hash_table_lookup(ht, "alpha".ptr) == &a);
    assert(g_hash_table_lookup(ht, "beta".ptr) == &a);
    assert(g_hash_table_lookup(ht, "gamma".ptr) is null);
    assert(g_hash_table_size(ht) == 2);
    assert(g_hash_table_remove(ht, "alpha".ptr) && g_hash_table_size(ht) == 1);
    g_hash_table_destroy(ht);
    assert(g_hash_table_lookup(ht, "beta".ptr) is null);
}
