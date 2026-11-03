/* D aliases for the glib basic types */

module gretl.glib;


alias gchar = char;
alias guchar = ubyte;
alias gshort = short;
alias gushort = ushort;
alias gint = int;
alias guint = uint;
alias glong = long;
alias gulong = ulong;

alias gint8 = byte;
alias guint8 = ubyte;
alias gint16 = short;
alias guint16 = ushort;
alias gint32 = int;
alias guint32 = uint;
alias gint64 = long;
alias guint64 = ulong;

alias gboolean = int;
alias gfloat = float;
alias gdouble = double;

alias gsize = size_t;
alias gssize = ptrdiff_t;

alias gpointer = void *;
alias gconstpointer = const(void) *;

/* Version macros */
enum GLIB_MAJOR_VERSION = 2;
enum GLIB_MINOR_VERSION = 88;

int G_ENCODE_VERSION(int major, int minor) {
    return (major << 16 | minor << 8);
}

enum GLIB_VERSION_2_28 = G_ENCODE_VERSION(2, 28);
/* Minor version is even, so this is the current minor version */
enum GLIB_VERSION_CUR_STABLE = G_ENCODE_VERSION(GLIB_MAJOR_VERSION, GLIB_MINOR_VERSION);

/* GDateWeekday (from glib's gdate.h), as used by dataset.c */
enum GDateWeekday {
    bad_weekday,
    monday,
    tuesday,
    wednesday,
    thursday,
    friday,
    saturday,
    sunday
}

/* glib's pointer/integer conversion macros */
gpointer GINT_TO_POINTER(int i)
{
    return cast(gpointer) cast(size_t) i;
}

int GPOINTER_TO_INT(gpointer p)
{
    return cast(int) cast(size_t) p;
}
