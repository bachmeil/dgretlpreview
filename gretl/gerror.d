/*
   Replacement for glib's GError (the part of it the port uses): a plain
   struct holding an error code and a message, allocated with malloc.

   - g_error_free() releases the struct and its message.
   - g_set_error_literal_msg() (not a glib function) is the helper the
     replacement modules use to create one: it copies @msg into a new GError
     stored in *error, if @error is non-null.

   glib's domains and the g_set_error() machinery are not provided; the code
   is always 1.

   Used by: gspawn.d, gfileutils.d, gstrfuncs.d, dataio.d.
*/

module gretl.gerror;


import core.stdc.stdlib : free, malloc;
import core.stdc.string : memcpy;

import gretl.glib : gchar;

struct GError {
    int code;
    gchar * message;
}

void g_error_free(GError * error)
{
    if (error != null) {
        free(error.message);
        free(error);
    }
}

void g_set_error_literal_msg(GError ** error, string msg)
{
    if (error != null) {
        GError * e = cast(GError *) malloc(GError.sizeof);

        e.code = 1;
        e.message = cast(gchar *) malloc(msg.length + 1);
        memcpy(e.message, msg.ptr, msg.length);
        e.message[msg.length] = 0;
        *error = e;
    }
}
