/* The two functions of lib/src/interact.c that the plotting code needs: check_stringvar_name(), which validates
   (and if wanted creates) the string variable that receives plot output, and set_plot_produced(), which records
   that a plot was written so that the command-execution code can announce it.

   interact.c itself is the command interpreter and GUI glue and is not ported; if it ever is, these two
   functions move back into it. The "plot_ok" flag of the C is exposed here through plot_was_produced() and
   clear_plot_produced() for whatever replaces the interpreter. */

module gretl.interact_extra;


import core.stdc.string : strcmp;

import gretl.dataio : check_identifier;
import gretl.gretl_errors : gretl_errmsg_sprintf;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_typemap : gretl_type_get_name;
import gretl.gretl_utils : gretl_type_from_name;
import gretl.libgretlh;
import gretl.uservar : user_var_add;

int check_stringvar_name(const(char) * name, int allow_new, const(DatasetRaw) * dset)
{
    GretlType t = gretl_type_from_name(name, dset);
    int err = 0;

    if (t != GretlType.none && t != GretlType.string_) {
        gretl_errmsg_sprintf(_("'%s' is of type %s").ptr, name, gretl_type_get_name(t));
        err = GretlError.types;
    } else if (t == GretlType.none) {
        if (allow_new) {
            /* create the variable if possible */
            err = check_identifier(name);
            if (!err) {
                err = user_var_add(name, GretlType.string_, null);
            }
        } else {
            /* otherwise require that the variable already exists */
            gretl_errmsg_sprintf(_("'%s' : not a string variable").ptr, name);
            err = GretlError.data;
        }
    }

    return err;
}

private __gshared int plot_ok;

void set_plot_produced()
{
    plot_ok = 1;
}

/* accessors for the flag, which the C interpreter reads and clears directly */

int plot_was_produced()
{
    return plot_ok;
}

void clear_plot_produced()
{
    plot_ok = 0;
}
