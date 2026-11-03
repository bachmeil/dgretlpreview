/*
   Replacement for the glib GPtrArray, backed by a D dynamic array.

   Why this file exists: the C code uses glib's GPtrArray (a growable array
   of pointers). Rather than bind to glib, this module provides the glib
   names and the public fields (pdata, len) so that ported call sites read
   like the C, e.g. array.pdata[i].

   Limitations (deliberate, to keep this small):

   - Only the functions needed so far are provided (new, add, free). Add
     others here (remove_index, index, set_size, free_func, ...) as ported
     files need them.
   - pdata and len are refreshed after every add. A pdata pointer obtained
     earlier may be invalidated by g_ptr_array_add, as with glib.
   - Storage is GC-allocated, so this is a GC dependency. If a betterc build
     is ever wanted, replace the body of this module with a hand-written
     growable array, keeping the same names and fields.
   - The array is not synchronised, as with glib.

   Used by: libset.d.
*/

module gretl.gptrarray;


import gretl.glib : gboolean, gpointer, guint;

struct GPtrArray {
    gpointer * pdata;      /* as in glib: the elements */
    guint len;             /* as in glib: the number of elements */
    private gpointer[] store;
}

GPtrArray * g_ptr_array_new()
{
    return new GPtrArray;
}

void g_ptr_array_add(GPtrArray * array, gpointer data)
{
    array.store ~= data;
    array.pdata = array.store.ptr;
    array.len = cast(guint) array.store.length;
}

/* glib frees the array, and the element segment too if free_seg is true;
   otherwise it returns the segment. Here the GC owns the memory, so the
   array is just emptied (and the segment returned if it is wanted). */

gpointer * g_ptr_array_free(GPtrArray * array, gboolean free_seg)
{
    gpointer * seg = free_seg ? null : array.pdata;

    array.store = null;
    array.pdata = null;
    array.len = 0;
    return seg;
}

unittest
{
    GPtrArray * a = g_ptr_array_new();
    int x, y;

    assert(a.len == 0);
    g_ptr_array_add(a, &x);
    g_ptr_array_add(a, &y);
    assert(a.len == 2 && a.pdata[0] == &x && a.pdata[1] == &y);
    g_ptr_array_free(a, 1);
    assert(a.len == 0 && a.pdata is null);
}
