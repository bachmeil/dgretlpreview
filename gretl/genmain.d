/* Port of lib/src/genmain.c: driver module for 'genr' and related commands (the header genmain.h is
   genmainh.d).

   The syntax-tree evaluator (geneval1.d - geneval5.d, which supply realgen(), gen_save_or_print(), gen_cleanup()
   and node_get_list()) is ported; the lexer and parser proper (genlex.c, genparse.c) are not yet, so the
   remaining parser entry points are declared extern (C) below.

   Differences from the C:
   - the GDEBUG diagnostics (compiled out in the C) are dropped;
   - the C's untyped "targtype" values (the GenSym members and, for array targets, a GretlType) are ints;
   - genflags values are the GenFlags enum; realgen() takes the flags as an int.
   See bugs.md #90-#93 for suspected upstream problems found while porting. */

module gretl.genmain;


import core.stdc.stdio : fprintf, sprintf, sscanf, stderr;
import core.stdc.stdlib : free, malloc;
import core.stdc.string : strcat, strchr, strcmp, strcpy, strlen, strncat, strncmp, strspn;
import core.stdc.ctype : isdigit;

import gretl.dataset : get_current_dataset, series_attach_string_table, series_delete_metadata, series_get_stack_level,
    series_is_listarg, series_set_flag, series_set_label, series_set_lag, series_set_mtime, series_set_parent,
    series_set_transform, set_dataset_is_changed;
import gretl.dataseth : VarFlags;
import gretl.geneval2 : node_get_list;
import gretl.geneval5 : gen_cleanup, gen_save_or_print, realgen;
import gretl.genmainh : ModelDataIndex;
import gretl.genparseh;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_errors : gretl_errmsg_sprintf;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_matrix : gretl_matrix_alloc, gretl_matrix_copy;
import gretl.gretl_matrixh : GretlMatrixRaw, gretl_is_null_matrix, gretl_matrix_is_scalar;
import gretl.gretl_prn : Prn, pprintf, pputc, pputs;
import gretl.gretl_typemaph : NUMERIC_TYPE;
import gretl.gretl_utils : copyvec, gretl_int_from_double, positive_int_from_string;
import gretl.libgretlh;
import gretl.libset : gretl_messages_on;
import gretl.missingh : NADBL, na;
import gretl.strutils : gretl_namechar_spn, gretl_strdup, string_is_blank;
import gretl.printscan : generate_obs_markers;
import gretl.uservar : get_list_by_name, gretl_scalar_get_value;
import gretl.genfuncs : gen_panel_dummies, gen_seasonal_dummies, gen_time, gen_unit, gen_wkday;
import gretl.gretl_commands : gretl_command_number;
import gretl.usermat : get_matrix_by_name;

/* unported dependencies: per-file declarations from the 2026c headers */
extern (C) {
    int gretl_function_depth();
    int function_lookup(const(char) * s);
}

private bool setting_obsval(const(Parser) * p)
{
    return (p.flags & GenFlags.obsval) != 0;
}

private void write_scalar_message(const(Parser) * p, Prn * prn)
{
    double x = gretl_scalar_get_value(p.lh.name.ptr, null);

    if (p.lh.t == GenSym.num) {
        pprintf(prn, _("Replaced scalar %s").ptr, p.lh.name.ptr);
    } else {
        pprintf(prn, _("Generated scalar %s").ptr, p.lh.name.ptr);
    }

    if (na(x)) {
        pputs(prn, " = NA");
    } else {
        pprintf(prn, " = %g", x);
    }
}

/* Note: p->lh.name may be empty but it's never NULL, while
   the vname member of a NODE may be NULL.
*/

private void gen_write_message(const(Parser) * p, int oldv, Prn * prn)
{
    const(char) * name = p.lh.name.ptr;
    int targ = p.targ;
    int vnum = p.lh.vnum;
    int t = p.lh.t;

    if (p.lhres != null) {
        /* compound LHS object */
        Node * lhs = cast(Node *) p.lhres;
        Node * lh1 = lhs.L;

        if (lh1.vname != null) {
            /* not an "embedded" object */
            name = lh1.vname;
            t = targ = lh1.t;
            if (t == GenSym.series) {
                vnum = lh1.vnum;
                targ = GenSym.num;
            }
        } else {
            return;
        }
        if (t != GenSym.series && t != GenSym.list && t != GenSym.mat) {
            /* we'll not print a message for modification
               of container types (bundle, array)
            */
            return;
        }
    }

    if (targ == GenSym.num) {
        if (setting_obsval(p)) {
            /* setting specific observation in series */
            pprintf(prn, _("Modified series %s (ID %d)").ptr,
                    name, vnum);
        } else {
            write_scalar_message(p, prn);
        }
    } else if (targ == GenSym.series) {
        if (vnum < oldv) {
            pprintf(prn, _("Replaced series %s (ID %d)").ptr, name, vnum);
        } else {
            pprintf(prn, _("Generated series %s (ID %d)").ptr, name, vnum);
        }
    } else if (targ == GenSym.mat) {
        GretlMatrixRaw * m = get_matrix_by_name(name);

        if (p.lhres != null) {
            pprintf(prn, _("Modified matrix %s").ptr, name);
        } else if (t == GenSym.mat) {
            pprintf(prn, _("Replaced matrix %s").ptr, name);
        } else {
            pprintf(prn, _("Generated matrix %s").ptr, name);
        }
        if (gretl_matrix_is_scalar(m)) {
            pprintf(prn, " = {%g}", m.val[0]);
        }
    } else if (targ == GenSym.list) {
        if (p.lhres != null) {
            pprintf(prn, _("Modified list %s").ptr, name);
        } else if (t == GenSym.list) {
            pprintf(prn, _("Replaced list %s").ptr, name);
        } else {
            pprintf(prn, _("Generated list %s").ptr, name);
        }
    } else if (targ == GenSym.str) {
        if (p.lhres != null) {
            pprintf(prn, _("Modified string %s").ptr, name);
        } else if (t == GenSym.str) {
            pprintf(prn, _("Replaced string %s").ptr, name);
        } else {
            pprintf(prn, _("Generated string %s").ptr, name);
        }
    } else {
        return;
    }

    pputc(prn, '\n');
}

private int maybe_record_lag_info(Parser * p)
{
    const(char) * s = p.input;
    int n = cast(int) strlen(p.lh.name.ptr);
    char[VNAMELEN] vname = 0;
    char[16] fmt = 0;
    int lag;

    if (!strncmp(s, "genr ", 5)) {
        s += 5;
    } else if (!strncmp(s, "series ", 7)) {
        s += 7;
    }

    s += strspn(s, " ");

    if (!strncmp(s, p.lh.name.ptr, n)) {
        s += n;
        s += strspn(s, " ");
        if (*s == '=') s++;
        s += strspn(s, " ");
    }

    /* ensure we don't preserve stale metadata */
    if (p.lh.vnum > 0) {
        series_delete_metadata(p.dset, p.lh.vnum);
    }

    sprintf(fmt.ptr, "%%%d[^ ()](%%d)", VNAMELEN-1);

    if (sscanf(s, fmt.ptr, vname.ptr, &lag) == 2) {
        s = strchr(s, ')');
        if (s != null && string_is_blank(s + 1)) {
            int pv = series_index(p.dset, vname.ptr);

            if (pv < p.dset.v && function_lookup(vname.ptr)) {
                /* rule out the case of a series name shadowing
                   a built-in function
                */
                pv = -1;
            }
            if (pv > 0 && pv < p.dset.v) {
                series_set_parent(p.dset, p.lh.vnum, p.dset.varname[pv]);
                series_set_transform(p.dset, p.lh.vnum, GretlCmdIndex.lags);
                series_set_lag(p.dset, p.lh.vnum, -lag);
            }
        }
    }

    return 0;
}

private int get_genseries_target(Parser * p)
{
    if (p.tree.t == GenFunc.f_genseries && (p.flags & GenFlags.ufret)) {
        return cast(int) p.ret.v.xval;
    } else {
        return 0;
    }
}

private void series_write_metadata(Parser * p, int oldv)
{
    const(char) * src = null;
    int vnum = 0;

    if (p.targ == GenSym.series) {
        vnum = p.lh.vnum;
    } else if (p.tree != null) {
        vnum = get_genseries_target(p);
    }

    if (vnum == 0) {
        /* this is relevant only for series */
        return;
    }

    if (p.lh.expr != null) {
        /* don't touch the metadata if we generated a single
           observation in a series
        */
        return;
    }

    if (p.targ == GenSym.series && p.input != null) {
        maybe_record_lag_info(p);
    }

    if (vnum < oldv && p.targ == GenSym.series) {
        series_set_mtime(p.dset, vnum);
    }

    if (p.lh.label != null && (p.flags & GenFlags.ufret)) {
        src = p.lh.label;
    } else if (p.rhs != null && strcmp(p.rhs, "NA")) {
        src = p.rhs;
    }

    if (src != null && *src != '\0') {
        if (strlen(src) > MAXLABEL - 1) {
            /* truncate if necessary */
            char[MAXLABEL] tmp = 0;

            tmp[0] = '\0';
            strncat(tmp.ptr, src, MAXLABEL - 4);
            strcat(tmp.ptr, "...");
            series_set_label(p.dset, vnum, tmp.ptr);
        } else {
            series_set_label(p.dset, vnum, src);
        }
        if (src == p.rhs) {
            /* in case the label is a formula */
            series_set_flag(p.dset, vnum, VarFlags.generated);
        }
    }

    if (p.lh.stab != null) {
        series_attach_string_table(p.dset, vnum, p.lh.stab);
        p.lh.stab = null;
    }
}

private immutable string[] reswords = [
    /* constants */
    "const",
    "NA",
    "null",
    "empty",
    "obs", /* not exactly a constant, but hey */
    /* types */
    "scalar",
    "series",
    "matrix",
    "string",
    "list",
    "bundle",
    "array",
    "void",
    /* control flow */
    "for",
    /* debugging instructions, etc. */
    "continue",
    "next",
    "to"
];

/**
 * gretl_reserved_word:
 * @str: string to be tested.
 *
 * Returns: non-zero if @str is a reserved word that cannot
 * figure as the name of a user-defined variable, otherwise 0.
 */

int gretl_reserved_word(const(char) * str)
{
    int n = cast(int) reswords.length;
    int i, ret = gretl_command_number(str);

    /* the names of built-in functions are deliberately
       not reserved */

    for (i=0; i<n && !ret; i++) {
        if (!strcmp(str, reswords[i].ptr)) {
            ret = GretlError.invarg;
        }
    }

    if (ret) {
        gretl_errmsg_sprintf(_("'%s' is a reserved word").ptr, str);
    }

    return ret;
}

/**
 * extract_varname:
 * @targ: target string into which to write name.
 * @src: source string.
 * @len: location to receive the length of the extracted portion.
 *
 * Writes up to #VNAMELEN - 1 characters from @s into @vname.
 *
 * Returns: 0 on success, non-zero if the number of valid varname
 * characters in @s is greater than #VNAMELEN - 1.
 */

int extract_varname(char * targ, const(char) * src, int * len)
{
    int err = 0;

    *targ = '\0';
    *len = gretl_namechar_spn(src);

    if (*len >= VNAMELEN) {
        /* too long to be a valid variable name */
        err = GretlError.unkvar;
    } else {
        strncat(targ, src, *len);
    }

    return err;
}

private int try_for_listvar(const(DatasetRaw) * dset, const(char) * s)
{
    char[VNAMELEN] vname = 0;
    char[VNAMELEN] lname = 0;
    char[16] fmt = 0;

    sprintf(fmt.ptr, "%%%d[^.].%%%ds", VNAMELEN-1, VNAMELEN-1);

    if (sscanf(s, fmt.ptr, lname.ptr, vname.ptr) == 2) {
        int * list = get_list_by_name(lname.ptr);

        if (list != null) {
            int i, vi;

            for (i=1; i<=list[0]; i++) {
                vi = list[i];
                if (!strcmp(vname.ptr, dset.varname[vi])) {
                    return vi;
                }
            }
        }
    }

    return dset.v;
}

/* The condition for recognizing a series by name. If we're inside a
   user function (@fd > 0) it must exist at the current level of
   function execution and its tenure at that level must not just be
   the result of its being a member of a list that was passed as an
   argument. Otherwise it's just a name comparison.
*/

private int series_matched(const(char) * s, int fd,
                           const(DatasetRaw) * dset, int i)
{
    if (fd == 0) {
        return strcmp(dset.varname[i], s) == 0;
    } else {
        return fd == series_get_stack_level(dset, i) &&
            !series_is_listarg(dset, i, null) &&
            strcmp(dset.varname[i], s) == 0;
    }
}

/**
 * real_series_index:
 * @dset: data information struct.
 * @varname: name of variable to test.
 * @greatest: search down from highest ID.
 *
 * Returns: the ID number of the variable whose name is given,
 * or the next available ID number if there is no variable of
 * that name.
 */

private int real_series_index(const(DatasetRaw) * dset, const(char) * s,
                              int greatest)
{
    int ret = -1;

    if (dset != null) {
        int i, fd;

        ret = dset.v; /* initialize to "next" series ID */

        if (s == null || *s == '\0' || isdigit(*s)) {
            goto bailout;
        }
        if (strcmp(s, "const") == 0) {
            ret = 0;
            goto bailout;
        }
        if (strchr(s, '.') != null) {
            ret = try_for_listvar(dset, s);
            goto bailout;
        }

        fd = gretl_function_depth();

        if (greatest) {
            for (i=dset.v-1; i>0; i--) {
                if (series_matched(s, fd, dset, i)) {
                    ret = i;
                    break;
                }
            }
        } else {
            for (i=1; i<dset.v; i++) {
                if (series_matched(s, fd, dset, i)) {
                    ret = i;
                    break;
                }
            }
        }
    }

 bailout:

    return ret;
}

/**
 * series_index:
 * @dset: data information struct.
 * @varname: name of variable to test.
 *
 * Returns: the ID number of the variable whose name is given,
 * or the next available ID number if there is no variable of
 * that name.
 */

int series_index(const(DatasetRaw) * dset, const(char) * varname)
{
    return real_series_index(dset, varname, 0);
}

/**
 * series_greatest_index:
 * @dset: data information struct.
 * @varname: name of variable to test.
 *
 * Returns: the ID number of the variable whose name is given,
 * or the next available ID number if there is no variable of
 * that name. In contrast to series_index() this variant searches
 * down from the greatest current series ID.
 */

int series_greatest_index(const(DatasetRaw) * dset, const(char) * varname)
{
    return real_series_index(dset, varname, 1);
}

/**
 * caller_series_index:
 * @dset: data information struct.
 * @varname: name of series to test.
 *
 * Returns: the ID number of the series whose name is given,
 * at the caller's level of function execution, or -1 if there's
 * no matching series.
 */

int caller_series_index(const(DatasetRaw) * dset, const(char) * vname)
{
    int i, targ = gretl_function_depth() - 1;

    if (dset != null && targ >= 0) {
        for (i=0; i<dset.v; i++) {
            if (series_get_stack_level(dset, i) == targ &&
                !strcmp(dset.varname[i], vname)) {
                return i;
            }
        }
    }

    return -1;
}

int current_series_index(const(DatasetRaw) * dset, const(char) * vname)
{
    int v = -1;

    if (dset != null && dset.v > 0 &&
        vname != null && *vname != '\0') {
        v = real_series_index(dset, vname, 0);
        if (v >= dset.v) {
            v = -1;
        }
    }

    return v;
}

int gretl_is_series(const(char) * name, const(DatasetRaw) * dset)
{
    if (dset == null) {
        dset = get_current_dataset();
    }

    if (dset == null) {
        return 0;
    } else {
        int v = real_series_index(dset, name, 0);

        return (v >= 0 && v < dset.v);
    }
}

int genr_special_word(const(char) * s)
{
    if (!strcmp(s, "dummy") ||
        !strcmp(s, "cdummy") ||
        !strcmp(s, "timedum") ||
        !strcmp(s, "unitdum") ||
        !strcmp(s, "time") ||
        !strcmp(s, "index") ||
        !strcmp(s, "unit") ||
        !strcmp(s, "weekday")) {
        return 1;
    } else if (!strncmp(s, "dummy:", 6) && *(s+6) != '\0') {
        return 1;
    } else if (!strncmp(s, "cdummy:", 7) && *(s+7) != '\0') {
        return 1;
    } else {
        return 0;
    }
}

__gshared private GretlType genr_last_type;

GretlType genr_get_last_output_type()
{
    return genr_last_type;
}

private int gen_special(const(char) * s, const(char) * line,
                        DatasetRaw * dset, Prn * prn, Parser * p)
{
    string msg = null;
    int orig_v = dset.v;
    int write_label = 0;
    int vnum = -1;
    int err = 0;

    if (dset == null || dset.n == 0) {
        return GretlError.nodata;
    }

    if (!strcmp(s, "markers")) {
        return generate_obs_markers(line, dset);
    } else if (!strcmp(s, "dummy")) {
        err = gen_seasonal_dummies(dset, 0, 0);
        if (!err) {
            msg = N_("Periodic dummy variables generated.\n");
        }
    } else if (!strcmp(s, "cdummy")) {
        err = gen_seasonal_dummies(dset, 0, 1);
        if (!err) {
            msg = N_("Centered periodic dummy variables generated.\n");
        }
    } else if (!strcmp(s, "timedum")) {
        err = gen_panel_dummies(dset, Gretlopt.OPT_T, prn);
        if (!err) {
            msg = N_("Panel dummy variables generated.\n");
        }
    } else if (!strcmp(s, "unitdum")) {
        err = gen_panel_dummies(dset, Gretlopt.OPT_NONE, prn);
        if (!err) {
            msg = N_("Panel dummy variables generated.\n");
        }
    } else if (!strcmp(s, "time")) {
        err = gen_time(dset, 1, &vnum);
        write_label = 1;
    } else if (!strcmp(s, "index")) {
        err = gen_time(dset, 0, &vnum);
        write_label = 1;
    } else if (!strcmp(s, "unit")) {
        err = gen_unit(dset, &vnum);
        write_label = 1;
    } else if (!strcmp(s, "weekday")) {
        err = gen_wkday(dset, &vnum);
        write_label = 1;
    } else if (!strncmp(s, "dummy:", 6)) {
        int ref_ = positive_int_from_string(s + 6);

        if (ref_ > 0 && ref_ <= dset.pd) {
            err = gen_seasonal_dummies(dset, ref_, 0);
        } else {
            err = GretlError.invarg;
        }
        if (!err) {
            msg = N_("Periodic dummy variables generated.\n");
        }
    } else if (!strncmp(s, "cdummy:", 7)) {
        int ref_ = positive_int_from_string(s + 7);

        if (ref_ > 0 && ref_ <= dset.pd) {
            err = gen_seasonal_dummies(dset, ref_, 1);
        } else {
            err = GretlError.invarg;
        }
        if (!err) {
            msg = N_("Centered periodic dummy variables generated.\n");
        }
    }

    if (msg != null && gretl_messages_on()) {
        pputs(prn, _(msg).ptr);
    }

    if (!err && write_label) {
        strcpy(p.lh.name.ptr, s);
        p.lh.vnum = vnum;
        p.dset = dset;
        p.targ = GenSym.series;
        p.flags = cast(GenFlags) 0;
        p.lhres = null;
        p.err = 0;
        p.prn = prn;
        if (prn != null && gretl_messages_on()) {
            gen_write_message(p, orig_v, prn);
        }
    }

    if (dset.v > orig_v) {
        set_dataset_is_changed(dset, 1);
        genr_last_type = GretlType.series;
    }

    return err;
}

private int is_genr_special(const(char) * s, char * vspec, const(char) ** rem)
{
    if (strncmp(s, "genr ", 5)) {
        /* all specials begin with "genr " */
        return 0;
    }

    s += 5;
    while (*s == ' ') s++;

    if (genr_special_word(s)) {
        if (vspec != null) {
            *vspec = '\0';
            strncat(vspec, s, VNAMELEN - 1);
        }
        if (rem != null) {
            *rem = s;
        }
        return 1;
    }

    if (!strncmp(s, "markers", 7) && strchr(s, '=')) {
        if (vspec != null) {
            strcpy(vspec, "markers");
        }
        if (rem != null) {
            s = strchr(s, '=') + 1;
            while (*s == ' ') s++;
            *rem = s;
        }
        return 1;
    }

    return 0;
}

private int maybe_unassigned_fncall(const(char) * s)
{
    return s[strlen(s)-1] == ')';
}

private bool gen_silent(GenFlags f)
{
    return (f & (GenFlags.discard | GenFlags.priv | GenFlags.decl)) != 0;
}

int generate(const(char) * line, DatasetRaw * dset,
             GretlType gtype, Gretlopt opt,
             Prn * prn)
{
    char[VNAMELEN] vname = 0;
    const(char) * subline = null;
    GretlType out_t;
    int oldv, targtype = GenSym.unk;
    GenFlags flags = cast(GenFlags) 0;
    Parser p;

    if (line == null) {
        return GretlError.args;
    }

    if (gtype == GretlType.none) {
        flags |= GenFlags.discard;
    } else if (gtype == GretlType.double_) {
        targtype = GenSym.num;
    } else if (gtype == GretlType.series) {
        targtype = GenSym.series;
    } else if (gtype == GretlType.matrix) {
        targtype = GenSym.mat;
    } else if (gtype == GretlType.string_) {
        targtype = GenSym.str;
    } else if (gtype == GretlType.bundle) {
        targtype = GenSym.bundle;
    } else if (gtype == GretlType.list) {
        targtype = GenSym.list;
    } else if (gtype == GretlType.bool_) {
        targtype = GenSym.num;
        flags |= GenFlags.anon;
    } else if (gretl_array_type(gtype)) {
        targtype = gtype;
    }

    if (opt & Gretlopt.OPT_P) {
        /* internal use of generate() */
        flags |= GenFlags.priv;
    }
    if (opt & Gretlopt.OPT_Q) {
        flags |= GenFlags.quiet;
    }
    if (opt & Gretlopt.OPT_C) {
        flags |= GenFlags.catch_;
    }
    if (opt & Gretlopt.OPT_U) {
        flags |= GenFlags.ufret;
    }
    if (opt & Gretlopt.OPT_O) {
        /* special for function call, no assignment */
        if (maybe_unassigned_fncall(line)) {
            targtype = GenSym.empty;
            flags |= GenFlags.void_;
        } else {
            return GretlError.parse;
        }
    }

    oldv = (dset != null)? dset.v : 0;

    if (is_genr_special(line, vname.ptr, &subline)) {
        return gen_special(vname.ptr, subline, dset, prn, &p);
    }

    realgen(line, &p, dset, prn, cast(int) flags, targtype);

    if (!p.err && targtype != GenSym.empty) {
        gen_save_or_print(&p, prn);
        if (!p.err && !gen_silent(p.flags)) {
            series_write_metadata(&p, oldv);
            if (gretl_messages_on() && prn != null && !(opt & Gretlopt.OPT_Q)) {
                gen_write_message(&p, oldv, prn);
            }
        }
    }

    genr_last_type = out_t = genr_get_output_type(&p);
    if (out_t == GretlType.series || out_t == GretlType.list) {
        set_dataset_is_changed(dset, 1);
    }
    if (!p.err && gtype == GretlType.numeric && !NUMERIC_TYPE(out_t)) {
        p.err = GretlError.types;
    } else if (p.err == 1) {
        /* a fairly good guess? */
        p.err = GretlError.parse;
    }

    gen_cleanup(&p);

    return p.err;
}

/* Get a pointer to a sub-object of type @gtype (an array or bundle),
   answering to the specification in @spec. This gets called in the
   context of a "foreach" loop if we can't find a "top-level" array
   or bundle as data source; we then see if it could be an embedded
   object.
*/

void * genr_get_pointer(const(char) * spec, GretlType gtype, int * errp)
{
    GenFlags flags = GenFlags.priv | GenFlags.quiet;
    void * ret = null;
    Parser p;
    int err;

    err = realgen(spec, &p, null, null, cast(int) flags, GenSym.unk);

    if (!err) {
        if (p.tree.t != GenSym.bmemb && p.tree.t != GenSym.osl) {
            /* we didn't get a sub-object */
            err = GretlError.types;
        } else if (gtype == GretlType.array && p.ret.t == GenSym.array) {
            ret = p.ret.v.a;
            p.ret.v.a = null;
        } else if (gtype == GretlType.bundle && p.ret.t == GenSym.bundle) {
            ret = p.ret.v.b;
            p.ret.v.b = null;
        } else {
            err = GretlError.types;
        }
    }

    gen_cleanup(&p);

    if (errp != null) {
        *errp = err;
    }

    return ret;
}

/* retrieve a scalar result directly */

private double generate_scalar_full(const(char) * s, DatasetRaw * dset,
                                    Prn * prn, int * err)
{
    Parser p;
    double x = NADBL;

    *err = realgen(s, &p, dset, prn, cast(int) (GenFlags.priv | GenFlags.anon), GenSym.num);

    if (!*err) {
        if (p.ret.t == GenSym.mat) {
            GretlMatrixRaw * m = p.ret.v.m;

            if (gretl_matrix_is_scalar(m)) {
                x = p.ret.v.m.val[0];
            } else if (!gretl_is_null_matrix(m)) {
                fprintf(stderr, "generate_scalar: got %d x %d matrix\n",
                        m.rows, m.cols);
                *err = GretlError.types;
            }
        } else if (p.ret.t == GenSym.num) {
            x = p.ret.v.xval;
        } else {
            *err = GretlError.types;
        }
    } else if (*err == 1) {
        *err = GretlError.parse;
    }

    gen_cleanup(&p);

    return x;
}

double generate_scalar(const(char) * s, DatasetRaw * dset, int * err)
{
    return generate_scalar_full(s, dset, null, err);
}

/* retrieve a boolean result directly: called only from
   flow_control.c
*/

double generate_boolean(const(char) * s, DatasetRaw * dset, Prn * prn, int * err)
{
    double x = generate_scalar_full(s, dset, prn, err);

    return (*err || na(x)) ? NADBL : cast(double) (x != 0.0);
}

/* retrieve an integer result directly */

int generate_int(const(char) * s, DatasetRaw * dset, int * err)
{
    double x = generate_scalar_full(s, dset, null, err);
    int ret = -1;

    if (!*err) {
        ret = gretl_int_from_double(x, err);
    }

    return ret;
}

/* Execute statement @s, bypassing the command tokenizer,
   when we know that it contains a call to a function whose
   return value (if any) we do not want to assign. Prime
   example: the bundle-print function for a function package.
*/

int generate_void(const(char) * s, DatasetRaw * dset, Prn * prn)
{
    Parser p;
    int err;

    err = realgen(s, &p, dset, prn, cast(int) (GenFlags.priv | GenFlags.void_), GenSym.empty);
    gen_cleanup(&p);

    return err;
}

/* retrieve a series result directly */

double * generate_series(const(char) * s, DatasetRaw * dset, Prn * prn,
                         int * err)
{
    Parser p;
    double * x = null;

    *err = realgen(s, &p, dset, prn, cast(int) (GenFlags.priv | GenFlags.anon), GenSym.series);

    if (!*err) {
        Node * n = p.ret;

        if (n.t == GenSym.series) {
            if (n.flags & NodeFlags.tmp) {
                /* steal the generated series */
                x = n.v.xvec;
                n.v.xvec = null;
            } else {
                x = copyvec(n.v.xvec, p.dset.n);
            }
        } else {
            *err = GretlError.types;
        }
    } else if (*err == 1) {
        *err = GretlError.parse;
    }

    gen_cleanup(&p);

    return x;
}

/* retrieve a matrix result directly */

GretlMatrixRaw * generate_matrix(const(char) * s, DatasetRaw * dset,
                              int * err)
{
    GretlMatrixRaw * m = null;
    Parser p;

    *err = realgen(s, &p, dset, null, cast(int) (GenFlags.priv | GenFlags.anon), GenSym.mat);

    if (!*err) {
        Node * n = p.ret;

        if (n.t == GenSym.mat) {
            if (n.flags & NodeFlags.tmp) {
                /* steal the generated matrix */
                m = n.v.m;
                n.v.m = null;
            } else {
                m = gretl_matrix_copy(n.v.m);
                if (m == null) {
                    *err = GretlError.alloc;
                }
            }
        } else if (n.t == GenSym.num) {
            if (na(n.v.xval)) {
                *err = GretlError.nan;
            } else {
                m = gretl_matrix_alloc(1, 1);
                if (m == null) {
                    *err = GretlError.alloc;
                } else {
                    m.val[0] = n.v.xval;
                }
            }
        } else {
            *err = GretlError.types;
        }
    } else if (*err == 1) {
        *err = GretlError.parse;
    }

    gen_cleanup(&p);

    return m;
}

/* retrieve a string result directly */

char * generate_string(const(char) * s, DatasetRaw * dset, int * err)
{
    Parser p;
    char * ret = null;

    *err = realgen(s, &p, dset, null, cast(int) (GenFlags.priv | GenFlags.anon), GenSym.str);

    if (!*err) {
        Node * n = p.ret;

        if (n.t == GenSym.str) {
            if (n.flags & NodeFlags.tmp) {
                /* steal the generated string */
                ret = n.v.str;
                n.v.str = null;
            } else {
                ret = gretl_strdup(n.v.str);
            }
        } else {
            *err = GretlError.types;
        }
    } else if (*err == 1) {
        *err = GretlError.parse;
    }

    gen_cleanup(&p);

    return ret;
}

/* retrieve a list result directly */

int * generate_list(const(char) * s, DatasetRaw * dset, int ci, int * err)
{
    int * ret = null;
    GenFlags flags = GenFlags.priv | GenFlags.anon;
    Parser p;

    if (dset == null) {
        *err = GretlError.nodata;
        return null;
    }

    if (ci == GretlCmdIndex.print) {
        flags |= GenFlags.prnlist;
    }

    *err = realgen(s, &p, dset, null, cast(int) flags, GenSym.list);

    if (!*err) {
        ret = node_get_list(p.ret, &p);
        *err = p.err;
    }

    gen_cleanup(&p);

    return ret;
}

/* retrieve a bundle result directly */

GretlBundle * generate_bundle(const(char) * s, DatasetRaw * dset, int * err)
{
    GretlBundle * ret = null;
    GenFlags flags = GenFlags.priv | GenFlags.anon;
    Parser p;

    *err = realgen(s, &p, dset, null, cast(int) flags, GenSym.bundle);

    if (!*err) {
        Node * n = p.ret;

        if (n.t == GenSym.bundle) {
            ret = n.v.b;
            n.v.b = null;
        }
    }

    gen_cleanup(&p);

    return ret;
}

void * generate_gretl_object(const(char) * s, DatasetRaw * dset,
                             GretlType * type, int * err)
{
    void * ret = null;
    GenFlags flags = GenFlags.priv | GenFlags.anon;
    Parser p;

    *type = GretlType.none;
    *err = realgen(s, &p, dset, null, cast(int) flags, GenSym.unk);

    if (!*err) {
        Node * n = p.ret;

        if (n.t == GenSym.num) {
            *type = GretlType.double_;
            ret = &n.v.xval;
        } else if (n.t == GenSym.mat) {
            *type = GretlType.matrix;
            ret = n.v.m;
            n.v.m = null;
        } else if (n.t == GenSym.str) {
            *type = GretlType.string_;
            ret = n.v.str;
            n.v.str = null;
        } else if (n.t == GenSym.bundle) {
            *type = GretlType.bundle;
            ret = n.v.b;
            n.v.b = null;
        } else if (n.t == GenSym.array) {
            *type = GretlType.array;
            ret = n.v.a;
            n.v.a = null;
        } else if (n.t == GenSym.list) {
            *type = GretlType.list;
            ret = n.v.ivec;
            n.v.ivec = null;
        } else if (n.t == GenSym.series) {
            *type = GretlType.series;
            ret = n.v.xvec;
            n.v.xvec = null;
        }
    }

    gen_cleanup(&p);

    return ret;
}

/* create a parsed tree that can be evaluated later,
   probably multiple times */

Parser * genr_compile(const(char) * s, DatasetRaw * dset,
                      GretlType gtype, Gretlopt opt,
                      Prn * prn, int * err)
{
    Parser * p;
    GenFlags flags = GenFlags.compile;
    int targtype = GenSym.unk;

    if (is_genr_special(s, null, null)) {
        *err = GretlError.eqn;
        return null;
    }

    p = cast(Parser *) malloc((*p).sizeof);
    if (p == null) {
        *err = GretlError.alloc;
        return null;
    }

    if (gtype == GretlType.none) {
        flags |= GenFlags.discard;
    } else if (gtype == GretlType.double_) {
        targtype = GenSym.num;
    } else if (gtype == GretlType.series) {
        targtype = GenSym.series;
    } else if (gtype == GretlType.matrix) {
        targtype = GenSym.mat;
    } else if (gtype == GretlType.string_) {
        targtype = GenSym.str;
    } else if (gtype == GretlType.bundle) {
        targtype = GenSym.bundle;
    } else if (gtype == GretlType.list) {
        targtype = GenSym.list;
    } else if (gtype == GretlType.bool_) {
        targtype = GenSym.num;
        flags |= GenFlags.anon;
    } else if (gretl_array_type(gtype)) {
        targtype = gtype;
    }

    if (opt & Gretlopt.OPT_P) {
        /* internal use of generate() */
        flags |= GenFlags.priv;
    }
    if (opt & Gretlopt.OPT_O) {
        /* special for function call, no assignment */
        targtype = GenSym.empty;
        flags |= GenFlags.void_;
    }
    if (opt & Gretlopt.OPT_N) {
        /* "no exec": compile but don't run */
        flags |= GenFlags.noexec;
    }
    if (opt & Gretlopt.OPT_A) {
        /* anonymous: no assignment to named variable */
        flags |= GenFlags.anon;
    }

    *err = realgen(s, p, dset, prn, cast(int) flags, targtype);

    if (*err == 0 && p != null &&
        !(opt & Gretlopt.OPT_N) && p.targ != GenSym.empty) {
        gen_save_or_print(p, prn);
        if (p.err) {
            *err = p.err;
        }
    }

    if (*err) {
        destroy_genr(p);
        p = null;
    }

    return p;
}

/* run a previously compiled generator */

int execute_genr(Parser * p, DatasetRaw * dset, Prn * prn)
{
    int oldv = (dset != null)? dset.v : 0;

    realgen(null, p, dset, prn, cast(int) GenFlags.exec, GenSym.unk);

    if (!p.err && p.targ != GenSym.empty) {
        gen_save_or_print(p, prn);
        if (dset != null && !(p.flags & GenFlags.priv)) {
            series_write_metadata(p, oldv);
        }
    }

    if (p.err) {
        gen_cleanup(p);
    }

    return p.err;
}

double evaluate_scalar_genr(Parser * p, DatasetRaw * dset,
                            Prn * prn, int * err)
{
    int flags = cast(int) (GenFlags.exec | GenFlags.priv | GenFlags.anon);
    double x = NADBL;

    *err = realgen(null, p, dset, null, flags, GenSym.num);

    if (!*err) {
        if (p.ret.t == GenSym.mat) {
            GretlMatrixRaw * m = p.ret.v.m;

            if (gretl_matrix_is_scalar(m)) {
                x = p.ret.v.m.val[0];
            } else if (!gretl_is_null_matrix(m)) {
                fprintf(stderr, "evaluate_if_cond: got %d x %d matrix\n",
                        m.rows, m.cols);
                *err = GretlError.types;
            }
        } else if (p.ret.t == GenSym.num) {
            x = p.ret.v.xval;
        } else {
            *err = GretlError.types;
        }
    } else if (*err == 1) {
        *err = GretlError.parse;
    }

    gen_cleanup(p);

    return x;
}

double evaluate_if_cond(Parser * p, DatasetRaw * dset, Prn * prn, int * err)
{
    double x = evaluate_scalar_genr(p, dset, prn, err);

    if (!*err && na(x)) {
        return NADBL;
    } else {
        return *err ? x : cast(double) (x != 0.0);
    }
}

/* destroy a previously compiled generator */

void destroy_genr(Parser * p)
{
    if (p != null) {
        p.flags = (p.flags & GenFlags.altinp) ? GenFlags.altinp : cast(GenFlags) 0;
        gen_cleanup(p);
        free(p);
    }
}

GretlType genr_get_output_type(const(Parser) * p)
{
    GretlType t = GretlType.none;

    if (!p.err) {
        if (p.targ == GenSym.num) {
            t = GretlType.double_;
        } else if (p.targ == GenSym.series) {
            t = GretlType.series;
        } else if (p.targ == GenSym.mat) {
            t = GretlType.matrix;
        } else if (p.targ == GenSym.list) {
            t = GretlType.list;
        } else if (p.targ == GenSym.bundle) {
            t = GretlType.bundle;
        } else if (p.targ == GenSym.array) {
            t = GretlType.array;
        }
    }

    return t;
}

int genr_get_output_varnum(const(Parser) * p)
{
    return p.lh.vnum;
}

GretlMatrixRaw * genr_get_output_matrix(Parser * p)
{
    GretlMatrixRaw * m = p.lh.mret;

    if (p.targ != GenSym.mat) {
        /* matrix under bundle or array?
           nullify pointer in case the type changes
        */
        p.lh.mret = null;
    }

    return m;
}

double genr_get_output_scalar(const(Parser) * p)
{
    if (p.targ == GenSym.num) {
        return gretl_scalar_get_value(p.lh.name.ptr, null);
    } else {
        return NADBL;
    }
}

int genr_no_assign(const(Parser) * p)
{
    return (p.flags & (GenFlags.discard | GenFlags.void_)) != 0;
}

int genr_is_autoregressive(const(Parser) * p)
{
    return (p.flags & GenFlags.autoreg) != 0;
}

void genr_set_na_check(Parser * p)
{
    p.flags |= GenFlags.natest;
}

void genr_unset_na_check(Parser * p)
{
    p.flags &= ~GenFlags.natest;
}
