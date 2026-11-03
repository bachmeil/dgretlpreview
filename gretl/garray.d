/*
   Replacement for the glib GArray (a growable array of fixed-size
   elements), backed by a D dynamic array of bytes.

   Why this file exists: so that the port does not depend on glib. The
   glib names and the public fields (data, len) are kept so that ported
   call sites read like the C, e.g. g_array_index(a, T, i) becomes
   (cast(T *) a.data)[i].

   Limitations (deliberate, to keep this small):

   - Only the functions needed so far are provided (new, append_vals,
     remove_index, free). Add others here (g_array_sized_new, insert_vals,
     set_size, ...) as ported files need them.
   - The zero_terminated and clear_ arguments to g_array_new are accepted
     but only zero_terminated is honoured (one extra zeroed element is
     kept after the data). Elements are always zero-initialised by the
     GC, which is a superset of clear_.
   - data and len are refreshed after every change. A data pointer
     obtained earlier may be invalidated by g_array_append_vals, as with
     glib.
   - Storage is GC-allocated, so this is a GC dependency. If a betterc
     build is ever wanted, replace the body of this module with a
     hand-written growable array, keeping the same names and fields.
   - The array is not synchronised, as with glib.

   Used by: gretl_prn.d.
*/

module gretl.garray;


import core.stdc.string : memmove, memcpy;

import gretl.glib : gboolean, gchar, gconstpointer, guint;

struct GArray {
    gchar * data;          /* as in glib: the elements */
    guint len;             /* as in glib: the number of elements */
    private ubyte[] store;
    private size_t elt_size;
    private bool zero_terminated;
}

private void garray_sync(GArray * a)
{
    a.data = cast(gchar *) a.store.ptr;
    a.len = cast(guint) (a.store.length / a.elt_size);
    if (a.zero_terminated) {
        a.len--;
    }
}

GArray * g_array_new(gboolean zero_terminated, gboolean clear_, guint element_size)
{
    GArray * a = new GArray;

    a.elt_size = element_size;
    a.zero_terminated = zero_terminated != 0;
    if (a.zero_terminated) {
        a.store = new ubyte[element_size];
    }
    garray_sync(a);
    return a;
}

GArray * g_array_append_vals(GArray * array, gconstpointer data, guint len)
{
    size_t esz = array.elt_size;
    size_t n = array.len;
    size_t nbytes = len * esz;

    /* keep the terminating zero element, if any, at the end */
    array.store.length = (n + len + (array.zero_terminated ? 1 : 0)) * esz;
    memcpy(array.store.ptr + n * esz, data, nbytes);
    if (array.zero_terminated) {
        array.store[(n + len) * esz .. (n + len + 1) * esz] = 0;
    }
    garray_sync(array);
    return array;
}

GArray * g_array_remove_index(GArray * array, guint index)
{
    size_t esz = array.elt_size;
    size_t n = array.len;

    assert(index < n);
    memmove(array.store.ptr + index * esz, array.store.ptr + (index + 1) * esz,
            (n - index - 1 + (array.zero_terminated ? 1 : 0)) * esz);
    array.store.length = array.store.length - esz;
    garray_sync(array);
    return array;
}

/* glib frees the array, and the element segment too if free_segment is
   true; otherwise it returns the segment. Here the GC owns the memory, so
   the array is just emptied (and the segment returned if it is wanted). */

gchar * g_array_free(GArray * array, gboolean free_segment)
{
    gchar * seg = free_segment ? null : array.data;

    array.store = null;
    array.data = null;
    array.len = 0;
    return seg;
}

unittest
{
    foreach (zt; [0, 1]) {
        GArray * a = g_array_new(zt, 0, cast(guint) int.sizeof);
        int[3] v = [10, 20, 30];

        assert(a.len == 0);
        foreach (i; 0 .. 3) {
            g_array_append_vals(a, &v[i], 1);
        }
        assert(a.len == 3);
        assert((cast(int *) a.data)[0 .. 3] == [10, 20, 30]);
        if (zt) {
            assert((cast(int *) a.data)[3] == 0);
        }
        g_array_remove_index(a, 1);
        assert(a.len == 2);
        assert((cast(int *) a.data)[0 .. 2] == [10, 30]);
        if (zt) {
            assert((cast(int *) a.data)[2] == 0);
        }
        g_array_remove_index(a, 1);
        g_array_remove_index(a, 0);
        assert(a.len == 0);
        g_array_free(a, 1);
        assert(a.data is null);
    }
}
