/*
   Replacement for glib's g_free, so that the port does not depend on glib.

   glib's g_free is free() under the hood (and accepts null), and the
   strings returned by the g_strdup family (see gstrfuncs.d) are allocated
   with malloc, so memory can be passed freely between g_free and the
   C library's free(). This matches the port's current "malloc/free now,
   maybe the GC later" approach.

   Only g_free is provided so far; add g_malloc, g_realloc, ... as ported
   files need them.

   Used by: gretl_prn.d, adf_kpss.d.
*/

module gretl.gmem;


import core.stdc.stdlib : free;

import gretl.glib : gpointer;

void g_free(gpointer mem)
{
    free(mem);
}
