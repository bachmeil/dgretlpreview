/* Port of varinfo_priv.h */

module gretl.varinfo_privh;


/* struct definition: basically private but needs to be shared between
   lib/src/dataset.c and plugin/purebin.c

   Warning (from the C header): this struct must not be modified, on
   pain of breaking the currently favoured means of saving gretl data in
   binary form (see purebin.c in the plugin directory). The "purebin"
   reader works on the assumption that the size and membership of this
   struct are known.

   The G_OS_WIN32 32-bit variant VARINFO64 is not ported.
*/

import gretl.glib : gint64;
import gretl.dataseth : SeriesTable, VarFlags;
import gretl.libgretlh : MAXDISP, VNAMELEN;

struct VarInfo {
    char * label;
    char[MAXDISP] display_name;
    char[VNAMELEN] parent;
    VarFlags flags;
    char compact_method;
    gint64 mtime;
    short transform;    /* note: command index of transform */
    short lag;
    short stack_level;
    short midas_period;
    char midas_freq;
    short orig_pd;
    SeriesTable * st;
}
