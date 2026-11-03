/* Port of lib/src/printscan.c: the printf, sprintf and sscanf commands/functions of the scripting language. */

module gretl.printscan;


import core.stdc.ctype : isdigit, isspace;
import core.stdc.errno : errno, ERANGE;
import core.stdc.math : exp, fabs, floor, isinf, isnan, log, pow, sqrt;
import core.stdc.stdio : fprintf, sprintf, snprintf, sscanf, stderr;
import core.stdc.stdlib : abs, atoi, calloc, free, malloc, realloc;
import core.stdc.string : memcpy, memset, strcmp, strcpy, strlen, strchr, strstr, strcat, strncat, strncmp, strspn, strcspn;

import std.algorithm : max, min;

import gretl.gretl_errorsh : GretlError;
import gretl.gretl_matrix;
import gretl.gretl_matrixh;
import gretl.libgretlh;
import gretl.missingh : NADBL, na;


import core.stdc.config : c_long, c_ulong;
import gretl.gretl_intl : gretl_pop_c_numeric_locale, gretl_push_c_numeric_locale;

/* the C library's strtol and friends take const char * and a plain char ** end pointer */
extern (C) {
    c_long strtol(const(char) * nptr, char ** endptr, int base);
    c_ulong strtoul(const(char) * nptr, char ** endptr, int base);
    double strtod(const(char) * nptr, char ** endptr);
}
import gretl.dataset : dataset_allocate_obs_markers, dataset_drop_last_variables, is_string_valued, series_get_string_for_obs;
import gretl.genmain : current_series_index, generate, generate_scalar, generate_string, genr_get_last_output_type;
import gretl.gmem : g_free;
import gretl.gretl_errors : gretl_errmsg_set, gretl_errmsg_sprintf, gretl_error_clear;
import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_get_buffer, gretl_print_new, gretl_print_reset_buffer, gretl_print_steal_buffer, pprintf, pputc, pputs;
import gretl.gretl_prnh : PrnType;
import gretl.gretl_utils : gretl_int_from_string;
import gretl.gstrfuncs : g_strdup_printf, g_strndup, g_utf8_offset_to_pointer, g_utf8_strlen;
import gretl.matrix_extra : gretl_matrix_print_with_format;
import gretl.strutils : dot_atof, gretl_strdup, gretl_strdup_printf, gretl_strndup, numeric_string, printf_escape, shift_string_left;
import gretl.uservarh : UserVar;
import gretl.printout : max_obs_marker_length, print_obs_marker;
import gretl.uservar : get_string_by_name, get_user_var_of_type_by_name, gretl_is_scalar, gretl_scalar_get_value, user_matrix_replace_matrix_by_name, user_var_delete_by_name, user_var_replace_value;
import gretl.genfuncs : gen_time;
import gretl.usermat : get_matrix_by_name;

/* unported dependencies: per-file declarations from the 2026c headers */
extern (C) {
}

void printf_series(int v, const(DatasetRaw) * dset, const(char) * fmt, int wid, int prec, int wstar, int pstar, Prn * prn)
{
    int obslen = max_obs_marker_length(dset);
    double x;
    int t;

    for (t=dset.t1; t<=dset.t2; t++) {
    print_obs_marker(t, dset, obslen, prn);
    x = dset.Z[v][t];
    if (!na(x)) {
        x *= (1 + 0x1p-52);
        if (wstar && pstar) {
        pprintf(prn, fmt, wid, prec, x);
        } else if (wstar || pstar) {
        wid = (wstar)? wid : prec;
        pprintf(prn, fmt, wid, x);
        } else {
        pprintf(prn, fmt, x);
        }
    }
    pputc(prn, '\n');
    }
}

char * printf_get_string(char * s, DatasetRaw * dset, int t, int * err)
{
    char * ret = null;

    if (!strncmp(s, "marker", 6) && dset != null && dset.S != null) {
    /* special, not in regular genr */
    const(char) * p = dset.S[t];
    const(char) * q = s + 6;
    int len = cast(int) strlen(p);
    int offset = 0;

    while (isspace(*q)) q++;

    if (*q == '+') {
        offset = atoi(++q);
        if (offset < 0) {
        offset = 0;
        }
        len -= offset;
    }
    if (len >= 0) {
        ret = gretl_strndup(p + offset, len);
        if (ret == null) {
        *err = GretlError.alloc;
        }
    }
    } else {
    ret = generate_string(s, dset, err);
    if (*err && dset != null && t >= 0) {
        /* using a string-valued series? */
        int v = current_series_index(dset, s);
        const(char) * str;

        if (v > 0 && is_string_valued(dset, v)) {
        str = series_get_string_for_obs(dset, v, t);
        if (str != null) {
            ret = gretl_strdup(str);
            *err = 0;
        }
        }
    }
    }

    if (ret == null && !*err) {
    *err = GretlError.data;
    }

    return ret;
}

/* here we're trying to get a field width or precision
   specifier */

int printf_get_int(const(char) * s, DatasetRaw * dset, int * err)
{
    int ret = 0;

    ret = gretl_int_from_string(s, err);

    if (*err) {
    *err = 0;
    gretl_error_clear();
    /* we'll allow a bit of slop with regard to
       integerhood here */
    ret = cast(int) generate_scalar(s, dset, err);
    }

    if (!*err && abs(ret) > 255) {
    *err = GretlError.data;
    }

    return ret;
}

/* gen_arg_val: when this function is called we have already
   determined that the argument string @s is not associated with a
   string format specifier, is not a literal numerical value, and is
   not the name of an existing scalar, matrix or series.  So now we
   attempt to generate a value (of type unknown).
*/

int gen_arg_val(const(char) * s, DatasetRaw * dset, double * px, GretlMatrixRaw ** pm, int * pv, int * scalar)
{
    const(char) * name = "$ptmp";
    char * genstr;
    int err = 0;

    genstr = gretl_strdup_printf("%s=%s", name, s);
    if (genstr == null) {
    return GretlError.alloc;
    }

    err = generate(genstr, dset, GretlType.any, Gretlopt.OPT_P, null);

    if (!err) {
    int type = genr_get_last_output_type();

    if (type == GretlType.double_) {
        *scalar = 1;
        *px = gretl_scalar_get_value(name, null);
        user_var_delete_by_name(name, null);
    } else if (type == GretlType.matrix) {
        GretlMatrixRaw * m = get_matrix_by_name(name);

        *pm = gretl_matrix_copy(m);
        if (*pm == null) {
        err = GretlError.alloc;
        }
        user_var_delete_by_name(name, null);
    } else if (type == GretlType.series) {
        *pv = current_series_index(dset, name);
    } else {
        err = GretlError.types; /* FIXME cleanup? */
    }
    }

    free(genstr);

    return err;
}

/* dup argv (@s) up to the next free comma */

char * get_next_arg(const(char) * s, int * len, int * err)
{
    const(char) * p;
    char * arg = null;
    int par = 0;
    int br = 0;
    int cbr = 0;
    int quoted = 0;
    int n = 0;

    if (s == null) {
    *err = GretlError.parse;
    return null;
    }

    *len = cast(int) strspn(s, ", ");
    s += *len;
    p = s;

    while (*p) {
    if (*p == '"' && (p == s || *(p-1) != '\\')) {
        quoted = !quoted;
    }
    if (!quoted) {
        if (*p == '(') par++;
        else if (*p == ')') par--;
        else if (*p == '[') br++;
        else if (*p == ']') br--;
        else if (*p == '{') cbr++;
        else if (*p == '}') cbr--;
    }
    if (!quoted && !par && !br && !cbr && *p == ',') {
        break;
    }
    p++;
    }

    n = cast(int) (p - s);
    *len += n;

    if (n > 0) {
    arg = gretl_strndup(s, n);
    if (arg == null) {
        *err = GretlError.alloc;
    }
    } else {
    *err = GretlError.args;
    }

    return arg;
}

/* conversion characters that may be preceded by 'l' */
bool l_ok(int c)
{
    return (c == 'd' || c == 'u' || c == 'x');
}

/* conversion characters that may be preceded by 'v' */
bool v_ok(int c)
{
    return (c == 'f' || c == 'g' || c == 'e' || c == 'E');
}

/* dup format string up to the end of the current conversion:
   e.g. "%10.4f", "%6g", "%.8g", %3s", "%.*s", "%#12.4g"
*/

char * get_printf_format_chunk(const(char) * s, char * conv, int * len, int * wstar, int * pstar, int * err)
{
    const(char) * cnvchars = "eEfgGduxs";
    const(char) * numchars = "0123456789";
    char * chunk = null;
    const(char) * p = s;
    int n;

    p++; /* move past % */

    /* '#', if present, must come first */
    if (*p == '#') {
    p++;
    }

    /* zero padding? */
    if (*p == '0') {
    p++;
    }

    /* left justification? */
    if (*p == '-') {
    p++;
    }

    /* leading space before positive number? */
    if (*p == ' ') {
    p++;
    }

    /* always print sign? */
    if (*p == '+') {
    p++;
    }

    /* optional width? */
    n = cast(int) strspn(p, numchars);
    if (n == 0 && *p == '*') {
    /* variable version */
    *wstar = 1;
    p++;
    } else {
    p += n;
    }

    /* optional dot plus precision? */
    if (*p == '.') {
    p++;
    n = cast(int) strspn(p, numchars);
    if (n == 0 && *p == '*') {
        *pstar = 1;
        p++;
    } else {
        p += n;
    }
    }

    /* now we should have a conversion character */
    if (*p == '\0' || strchr(cnvchars, *p) == null) {
    if (*p == 'l' && l_ok(*(p+1))) {
        conv[0] = *p;
        conv[1] = *(p+1);
        p++;
    } else if (*p == 'v') {
        if (v_ok(*(p+1))) {
        /* allow "vf", "vg" for matrices */
        conv[0] = *p;
        conv[1] = *(p+1);
        p++;
        } else {
        /* allow plain 'v' -> "vg" */
        conv[0] = *p;
        }
    } else {
        fprintf(stderr, "bad conversion '%c'\n".ptr, *p);
        *err = GretlError.parse;
        return null;
    }
    }

    if (!conv[0]) {
    conv[0] = *p;
    }
    p++;
    *len = n = cast(int) (p - s);

    if (n > 0) {
    chunk = gretl_strndup(s, n);
    if (chunk == null) {
        *err = GretlError.alloc;
    }
    } else {
    *err = GretlError.parse;
    }

    return chunk;
}

int printf_as_int(double x, int fc, int lflag, char * fmt, int wstar, int pstar, int wid, int prec, Prn * prn)
{
    if ((fc == 'u' || fc == 'x') && lflag) {
    if (x < 0 || x > cast(double) ulong.max) {
        return GretlError.invarg;
    }
    if (wstar && pstar) {
        pprintf(prn, fmt, wid, prec, cast(ulong) x);
    } else if (wstar || pstar) {
        wid = (wstar)? wid : prec;
        pprintf(prn, fmt, wid, cast(ulong) x);
    } else {
        pprintf(prn, fmt, cast(ulong) x);
    }
    } else if (fc == 'u' || fc == 'x') {
    if (x < 0 || x > cast(double) uint.max) {
        return GretlError.invarg;
    }
    if (wstar && pstar) {
        pprintf(prn, fmt, wid, prec, cast(uint) x);
    } else if (wstar || pstar) {
        wid = (wstar)? wid : prec;
        pprintf(prn, fmt, wid, cast(uint) x);
    } else {
        pprintf(prn, fmt, cast(uint) x);
    }
    } else if (fc == 'd' && lflag) {
    if (x < cast(double) long.min || x > cast(double) long.max) {
        return GretlError.invarg;
    }
    if (wstar && pstar) {
        pprintf(prn, fmt, wid, prec, cast(long) x);
    } else if (wstar || pstar) {
        wid = (wstar)? wid : prec;
        pprintf(prn, fmt, wid, cast(long) x);
    } else {
        pprintf(prn, fmt, cast(long) x);
    }
    } else {
    if (x < cast(double) int.min || x > cast(double) int.max) {
        return GretlError.invarg;
    }
    if (wstar && pstar) {
        pprintf(prn, fmt, wid, prec, cast(int) x);
    } else if (wstar || pstar) {
        wid = (wstar)? wid : prec;
        pprintf(prn, fmt, wid, cast(int) x);
    } else {
        pprintf(prn, fmt, cast(int) x);
    }
    }

    return 0;
}

/* Allow for UTF-8 strings where the number of characters
   is not equal to the number of bytes: we'll adjust a
   printf specification such as "%10s" so that the result
   comes out 10 characters wide, rather than 10 bytes.
*/

void maybe_adjust_string_format(char * str, int * wstar, int * pstar, int * wid, int * prec, char ** fmt)
{
    int adj = cast(int) (strlen(str) - g_utf8_strlen(str, -1));

    if (adj > 0) {
    const(char) * s = *fmt;
    const(char) * p;
    char * tmp = null;
    int uwid = 0, uprec = 0;

    if (*wstar) {
        /* got width in variable form */
        uwid = *wid;
    } else {
        /* maybe got inline width? */
        if (isdigit(s[1]) || s[1] == '-') {
        uwid = atoi(s + 1);
        }
    }

    if (*pstar) {
        /* got "precision" in variable form */
        p = g_utf8_offset_to_pointer(str, *prec);
        uprec = cast(int) (p - str);
    } else if ((p = strchr(s, '.')) != null) {
        /* got inline precision? */
        uprec = atoi(p + 1);
        p = g_utf8_offset_to_pointer(str, uprec);
        uprec = cast(int) (p - str);
    }

    if (uwid != 0) {
        if (uprec > 0 && uprec < strlen(str)) {
        char * trunc = g_strndup(str, uprec);

        adj = cast(int) (strlen(trunc) - g_utf8_strlen(trunc, -1));
        g_free(trunc);
        }
        uwid = (uwid < 0)? uwid - adj : uwid + adj;
    }

    if (uwid != 0 && uprec > 0) {
        tmp = g_strdup_printf("%%%d.%ds".ptr, uwid, uprec);
    } else if (uwid != 0) {
        tmp = g_strdup_printf("%%%ds".ptr, uwid);
    } else if (uprec > 0) {
        tmp = g_strdup_printf("%%.%ds".ptr, uprec);
    }

    if (tmp != null) {
        /* replace the incoming format */
        free(*fmt);
        *fmt = gretl_strdup(tmp);
        g_free(tmp);
        *wstar = *pstar = 0;
    }
    }
}

/* Extract the next conversion spec from *pfmt, find and evaluate the
   corresponding elements in *pargs, and print the result. Note
   that @t should be negative (and therefore ignored) unless
   we're doing the "genr markers" special thing.
*/

int print_arg(const(char) ** pfmt, const(char) ** pargs, DatasetRaw * dset, int t, Prn * prn)
{
    const(char) * intconv = "dxu";
    char * fmt = null;
    char * arg = null;
    char * str = null;
    char[3] conv = 0;
    GretlMatrixRaw * m = null;
    int free_m = 0;
    int series_v = -1;
    int free_v = 0;
    double x = NADBL;
    int flen = 0, alen = 0;
    int wstar = 0, pstar = 0;
    int wid = 0, prec = 0;
    int got_scalar = 0;
    int fc = 0;
    int err = 0;

    /* select current conversion format */
    fmt = get_printf_format_chunk(*pfmt, conv.ptr, &flen, &wstar, &pstar, &err);
    if (err) {
    return err;
    }

    *pfmt += flen;

    if (conv[0] == 'l' && l_ok(conv[1])) {
    /* e.g. "ld", "lx" */
    fc = conv[1];
    } else {
    fc = conv[0];
    }

    if (wstar) {
    /* evaluate field width specifier */
    arg = get_next_arg(*pargs, &alen, &err);
    if (!err) {
        x = printf_get_int(arg, dset, &err);
    }
    free(arg);
    if (err) {
        goto bailout;
    }
    *pargs += alen;
    wid = cast(int) x;
    }

    if (pstar) {
    /* evaluate precision specifier */
    arg = get_next_arg(*pargs, &alen, &err);
    if (!err) {
        x = printf_get_int(arg, dset, &err);
    }
    free(arg);
    if (err) {
        goto bailout;
    }
    *pargs += alen;
    prec = cast(int) x;
    }

    /* get next substantive arg */
    arg = get_next_arg(*pargs, &alen, &err);

    if (!err && !strcmp(arg, "time") &&
    current_series_index(dset, arg) < 0) {
    /* pre-generate "time" if needed */
    err = gen_time(dset, 1, null);
    }

    if (!err) {
    int v = -1;

    if (fc == 's') {
        /* must be printing a string */
        str = printf_get_string(arg, dset, t, &err);
    } else if (numeric_string(arg)) {
        /* printing a scalar value */
        x = dot_atof(arg);
        got_scalar = 1;
    } else if (gretl_is_scalar(arg)) {
        /* printing a named scalar */
        x = gretl_scalar_get_value(arg, null);
        got_scalar = 1;
    } else if ((m = get_matrix_by_name(arg)) != null) {
        /* printing a named matrix */
        if (t >= 0) {
        err = GretlError.invarg;
        }
    } else if ((v = current_series_index(dset, arg)) >= 0) {
        /* printing a named series, unless t >= 0 */
        if (t < 0) {
        series_v = v;
        } else {
        x = dset.Z[v][t];
        got_scalar = 1;
        }
    } else if (t >= 0) {
        /* doing "genr markers" */
        err = GretlError.invarg;
    } else {
        /* try treating as generated value */
        err = gen_arg_val(arg, dset, &x, &m, &v, &got_scalar);
        if (m != null) {
        free_m = 1;
        } else if (v > 0) {
        series_v = v;
        free_v = 1;
        }
    }
    *pargs += alen;
    }

    free(arg);

    if (!err && fc == 'v' && m == null) {
    /* 'v' for variant is for matrices only */
    err = GretlError.parse;
    }

    if (err) {
    goto bailout;
    }

    if (m != null && m.rows == 1 && m.cols == 1) {
    /* print 1 x 1 matrix as scalar? */
    if (!isnan(m.val[0])) {
        x = m.val[0];
        got_scalar = 1;
    }
    }

    if (got_scalar && na(x)) {
    fc = fmt[flen - 1] = 's';
    str = gretl_strdup("NA");
    if (conv[0] == 'l') {
        /* strip the 'l' length modifier: "%ld" -> "%s", so we
           don't form an invalid "%ls" (wide-string) conversion
        */
        fmt[flen - 2] = 's';
        fmt[flen - 1] = '\0';
    }
    }

    /* do the actual printing */

    if (m != null && !got_scalar) {
    /* printing a matrix */
    if (!wstar) wid = -1;
    if (!pstar) prec = -1;
    gretl_matrix_print_with_format(m, fmt, wid, prec, prn);
    } else if (series_v >= 0) {
    /* printing a series */
    printf_series(series_v, dset, fmt, wid, prec,
              wstar, pstar, prn);
    } else if (fc == 's') {
    /* printing a string */
    if (!strcmp(fmt+1, "s")) {
        pputs(prn, str);
    } else {
        maybe_adjust_string_format(str, &wstar, &pstar,
                       &wid, &prec, &fmt);
        if (wstar && pstar) {
        pprintf(prn, fmt, wid, prec, str);
        } else if (wstar || pstar) {
        wid = (wstar)? wid : prec;
        pprintf(prn, fmt, wid, str);
        } else {
        pprintf(prn, fmt, str);
        }
    }
    } else if (strchr(intconv, fc)) {
    /* printing a scalar as int */
    int lflag = (conv[0] == 'l');

    err = printf_as_int(x, fc, lflag, fmt, wstar, pstar,
                wid, prec, prn);
    } else if (got_scalar) {
    /* printing a (non-missing) scalar value */
    if (!isnan(x) && !isinf(x)) {
        char * s;

        if (pstar && prec >= 15) {
        ; /* don't mess with it */
        } else if (((s = strchr(fmt, '.')) !is null) && isdigit(*(s+1))) {
        if (atoi(s+1) < 15) {
            x *= (1 + 0x1p-52);
        }
        } else {
        x *= (1 + 0x1p-52);
        }
    }
    if (wstar && pstar) {
        pprintf(prn, fmt, wid, prec, x);
    } else if (wstar || pstar) {
        wid = (wstar)? wid : prec;
        pprintf(prn, fmt, wid, x);
    } else {
        pprintf(prn, fmt, x);
    }
    }

 bailout:

    if (err) {
    pputc(prn, '\n');
    }

    free(fmt);
    free(str);

    if (free_m) {
    gretl_matrix_free(m);
    } else if (free_v) {
    dataset_drop_last_variables(dset, 1);
    }

    return err;
}

/* supports both regular printf and generation of obs markers */

int real_do_printf(const(char) * format, const(char) * args, DatasetRaw * dset, Prn * inprn, int * nchars, int t)
{
    Prn * prn = null;
    int err = 0;

    if (inprn == null) {
    return 0;
    }

    gretl_error_clear();

    /* We'll buffer the output locally in case there's an
       error part-way through the printing.
    */
    prn = gretl_print_new(PrnType.buffer, &err);
    if (err) {
        return err;
    }

    if (format != null) {
    const(char) * p = format;
    const(char) * q = args;
    char c;

    while (*p && !err) {
        if (*p == '%' && *(p+1) == '%') {
        pputc(prn, '%');
        p += 2;
        } else if (*p == '%') {
        err = print_arg(&p, &q, dset, t, prn);
        } else if (*p == '\\') {
        if (*(p+1) == '\0') {
            /* trailing backslash: treat as literal, and don't
               step past the string terminator */
            pputc(prn, *p);
            p++;
        } else {
            c = printf_escape(*(p+1), &err);
            if (c == 0) {
            /* tolerate slop! */
            pputc(prn, *p);
            pputc(prn, *(p+1));
            } else {
            pputc(prn, c);
            }
            p += 2;
        }
        } else {
        pputc(prn, *p);
        p++;
        }
    }

    if (!err && q != null && *q != '\0') {
        char * s = g_strdup_printf(_("unprocessed argument(s): '%s'").ptr, q);

        gretl_errmsg_sprintf("printf: %s".ptr, s);
        g_free(s);
        err = GretlError.parse;
    }
    }

    if (err) {
    pputc(inprn, '\n');
    } else {
    const(char) * buf = gretl_print_get_buffer(prn);

    if (buf != null) {
        if (nchars != null) {
        *nchars = cast(int) strlen(buf);
        }
        pputs(inprn, buf);
    }
    }

    gretl_print_destroy(prn);

    return err;
}

/**
 * do_printf:
 * @format: format string.
 * @args: list of arguments as string.
 * @dset: dataset struct.
 * @prn: printing struct.
 * @nchars: location to receive number of characters printed.
 *
 * Implements a somewhat limited version of C's printf()
 * for use in gretl scripts.
 *
 * Returns: 0 on success, non-zero on error.
 */

int do_printf(const(char) * format, const(char) * args, DatasetRaw * dset, Prn * prn, int * nchars)
{
    return real_do_printf(format, args, dset, prn, nchars, -1);
}

/**
 * do_sprintf_function:
 * @format: format string.
 * @args: list of arguments as string.
 * @dset: dataset struct.
 * @err: location to receive error code.
 *
 * Implements the hansl-function form of sprintf.
 *
 * Returns: constructed string, or NULL on error.
 */

char * do_sprintf_function(const(char) * format, const(char) * args, DatasetRaw * dset, int * err)
{
    const(char) * p = format;
    const(char) * q = args;
    char * buf = null;
    char c;
    Prn * prn;

    if (format == null || *format == '\0') {
    gretl_errmsg_set("sprintf: format string is missing");
    *err = GretlError.data;
    return null;
    }

    prn = gretl_print_new(PrnType.buffer, err);
    if (*err) {
    fprintf(stderr, "sprintf: creation of buffer failed\n".ptr);
    return null;
    }

    while (*p && !*err) {
    if (*p == '%' && *(p+1) == '%') {
        pputc(prn, '%');
        p += 2;
    } else if (*p == '%') {
        *err = print_arg(&p, &q, dset, -1, prn);
    } else if (*p == '\\') {
        if (*(p+1) == '\0') {
        /* trailing backslash: treat as literal, and don't
           step past the string terminator */
        pputc(prn, *p);
        p++;
        } else {
        c = printf_escape(*(p+1), err);
        if (c == 0) {
            /* tolerate slop! */
            pputc(prn, *p);
            pputc(prn, *(p+1));
        } else {
            pputc(prn, c);
        }
        p += 2;
        }
    } else {
        pputc(prn, *p);
        p++;
    }
    }

    if (!*err && q != null && *q != '\0') {
    char * s = g_strdup_printf(_("unprocessed argument(s): '%s'").ptr, q);

    gretl_errmsg_sprintf("sprintf: %s".ptr, s);
    g_free(s);
    *err = GretlError.parse;
    }

    if (!*err) {
    buf = gretl_print_steal_buffer(prn);
    }

    gretl_print_destroy(prn);

    return buf;
}

/* below: sscanf apparatus */

struct BracketScan {
    int reverse;
    char * chrs;
};

int bscan_escape(char * s)
{
    int i, n = cast(int) strlen(s);

    for (i=0; i<n; i++) {
    if (s[i] == '\\' && (i == 0 || s[i-1] != '\\')) {
        if (s[i+1] == 'n') {
        s[i] = '\n';
        shift_string_left(s + i + 1, 1);
        i++;
        } else if (s[i+1] == 't') {
        s[i] = '\t';
        shift_string_left(s + i + 1, 1);
        i++;
        }
    }
    }

    return 0;
}

int parse_BracketScanner(const(char) ** pfmt, BracketScan * bscan)
{
    const(char) * p, s = *pfmt;
    int k = 0, n = 0, closed = 0;

    s++;
    k++;

    if (*s == '^') {
    bscan.reverse = 1;
    s++;
    k++;
    }
    if (*s == ']') {
    s++;
    k++;
    }

    p = s;
    while (*p) {
    if (*p == ']') {
        closed = 1;
        break;
    }
    p++;
    n++;
    }

    if (!closed || n == 0) {
    return GretlError.parse;
    }

    bscan.chrs = gretl_strndup(s, n);
    if (bscan.chrs == null) {
    return GretlError.alloc;
    }

    bscan_escape(bscan.chrs);

    *pfmt += (k + n + 1);

    return 0;
}

enum WIDTH_UNSPEC = -1;

char * get_scanf_format_chunk(const(char) * s, int * fc, int * len, int * width, int * supp, BracketScan * bscan, int * err)
{
    const(char) * cnvchars = "dfgsmx";
    const(char) * numchars = "0123456789";
    char * chunk = null;
    const(char) * p = s;
    int n;

    p++; /* % */

    if (*p == '*') {
    /* suppress assignment */
    *supp = 1;
    p++;
    }

    /* optional max width? */
    n = cast(int) strspn(p, numchars);
    if (n > 0) {
    *width = atoi(p);
    p += n;
    } else {
    *width = WIDTH_UNSPEC;
    }

    /* now we should have a conversion spec */
    if (*p == '[') {
    *err = parse_BracketScanner(&p, bscan);
    if (!*err) {
        *fc = 's';
    }
    } else {
    if (*p == 'l' && *(p+1) == 'f') {
        p++;
    }
    if (*p == '\0' || strchr(cnvchars, *p) == null) {
        fprintf(stderr, "bad conversion '%c'\n".ptr, *p);
        *err = GretlError.parse;
        return null;
    }
    *fc = *p++;
    }

    *len = n = cast(int) (p - s);

    if (n > 0) {
    chunk = gretl_strndup(s, n);
    if (chunk == null) {
        *err = GretlError.alloc;
    }
    } else {
    *err = GretlError.parse;
    }

    return chunk;
}

int scan_string(char * targ, const(char) ** psrc, int width, BracketScan * bscan, int * ns)
{
    int n, err = 0;

    *psrc += cast(int) strspn(*psrc, " \t\n");

    if (bscan.chrs != null) {
    if (bscan.reverse) {
        n = cast(int) strcspn(*psrc, bscan.chrs);
    } else {
        n = cast(int) strspn(*psrc, bscan.chrs);
    }
    } else {
    n = cast(int) strcspn(*psrc, " \t\n");
    }

    if (width != WIDTH_UNSPEC && n > width) {
    n = width;
    }

    if (targ != null) {
    UserVar * uvar;

    uvar = get_user_var_of_type_by_name(targ, GretlType.string_);

    if (uvar == null) {
        gretl_errmsg_sprintf(_("%s: not a string variable").ptr, targ);
        err = GretlError.data;
    } else if (n > 0) {
        char * conv = gretl_strndup(*psrc, n);

        if (conv == null) {
        err = GretlError.alloc;
        } else {
        user_var_replace_value(uvar, conv, GretlType.string_);
        *psrc += n;
        }
        if (!err) {
        *ns += 1;
        }
    }
    } else {
    *psrc += n;
    }

    return err;
}

int scan_scalar(char * targ, const(char) ** psrc, int fc, int width, DatasetRaw * dset, int * ns)
{
    char * endp = null;
    long k = 0;
    ulong u = 0;
    double x = 0;
    int err = 0;

    errno = 0;

    if (fc == 'd') {
    k = strtol(*psrc, &endp, 10);
    } else if (fc == 'x') {
    u = strtoul(*psrc, &endp, 16);
    } else {
    x = strtod(*psrc, &endp);
    }

    if (width != WIDTH_UNSPEC) {
    int eaten = cast(int) (endp - *psrc);

    if (eaten > width) {
        char * tmp = gretl_strndup(*psrc, width);

        if (tmp == null) {
        return GretlError.alloc;
        }
        if (fc == 'd') {
        k = strtol(tmp, &endp, 10);
        } else if (fc == 'x') {
        u = strtoul(tmp, &endp, 16);
        } else {
        x = strtod(tmp, &endp);
        }
        endp = cast(char *) *psrc + width;
        free(tmp);
    }
    }

    /* "genr" operates in the C locale regardless */
    gretl_push_c_numeric_locale();

    if (!errno && endp != *psrc && targ != null) {
    char * genline;

    if (fc == 'd') {
        genline = gretl_strdup_printf("%s=%ld", targ, k);
    } else if (fc == 'x') {
        genline = gretl_strdup_printf("%s=%lu", targ, u);
    } else {
        genline = gretl_strdup_printf("%s=%.16g", targ, x);
    }
    /* we use GRETL_TYPE_ANY below to allow for the possibility
       that we're assigning to, e.g., an element of a matrix
    */
    err = generate(genline, dset, GretlType.any, Gretlopt.OPT_P, null);
    if (!err) {
        *ns += 1;
    }
    free(genline);
    }

    /* match the push above */
    gretl_pop_c_numeric_locale();

    *psrc = endp;

    return err;
}

int scan_matrix(char * targ, const(char) ** psrc, int rows, int * ns)
{
    const(char) * src = *psrc;
    char * endp = null;
    GretlMatrixRaw * m = null;
    double x;
    int r = 0, c = 0, c0 = 0;
    int i, j, err = 0;

    errno = 0;

    while (*src) {
    strtod(src, &endp);
    if (endp == src || errno == ERANGE) {
        break;
    }
    c++;
    src = endp;
    src += cast(int) strspn(src, " \t");
    if (*src == '\r' && *(src+1) == '\n') {
        src++;
    }
    if (*src == '\n' || *src == '\r' || *src == '\0') {
        if (r == 0) {
        c0 = c;
        } else if (c != c0) {
        break;
        }
        r++;
        if (r == rows) {
        break;
        }
        c = 0;
    }
    }

    if (targ == null) {
    *psrc = src;
    return 0;
    }

    if (get_matrix_by_name(targ) == null) {
    return GretlError.unkvar;
    }

    if (r > 0 && c0 > 0) {
    m = gretl_matrix_alloc(r, c0);
    if (m == null) {
        return GretlError.alloc;
    }

    src = *psrc;
    for (i=0; i<r; i++) {
        for (j=0; j<c0; j++) {
        x = strtod(src, &endp);
        gretl_matrix_set(m, i, j, x);
        src = endp;
        src += cast(int) strspn(src, " \t\n\r");
        }
    }

    err = user_matrix_replace_matrix_by_name(targ, m);
    if (err) {
        gretl_matrix_free(m);
    } else {
        /* number of elements scanned */
        *ns += m.rows * m.cols;
    }
    }

    *psrc = src;

    return err;
}

int scan_arg(const(char) ** psrc, const(char) ** pfmt, const(char) ** pargs, DatasetRaw * dset, int * ns)
{
    char * fmt = null;
    char * arg = null;
    char * str = null;
    char * targ = null;
    BracketScan bscan;
    int flen = 0, alen = 0;
    int wid = 0, supp = 0;
    int fc = 0;
    int err = 0;

    bscan.reverse = 0;
    bscan.chrs = null;

    /* select current conversion format */
    fmt = get_scanf_format_chunk(*pfmt, &fc, &flen, &wid, &supp,
                 &bscan, &err);
    if (err) {
    return err;
    }

    *pfmt += flen;

    if (!supp) {
    /* get next substantive arg */
    arg = get_next_arg(*pargs, &alen, &err);
    if (!err) {
        targ = arg;
        if (*targ == '&') {
        targ++;
        }
        if (*targ == '\0') {
        err = GretlError.parse;
        } else {
        *pargs += alen;
        }
    }
    }

    /* do the actual scanning */

    if (!err && fc == 's') {
    err = scan_string(targ, psrc, wid, &bscan, ns);
    } else if (!err && fc == 'm') {
    err = scan_matrix(targ, psrc, wid, ns);
    } else if (!err) {
    err = scan_scalar(targ, psrc, fc, wid, dset, ns);
    }

    free(fmt);
    free(arg);
    free(str);
    free(bscan.chrs);

    return err;
}

/**
 * do_sscanf:
 * @src: the source string.
 * @format: the format string.
 * @args: the arguments, as composite string.
 * @dset: dataset struct.
 * @n_items: location to receive the number of scanned items.
 *
 * Implements a somewhat simplified version of C's sscanf()
 * for use in hansl scripts.
 *
 * Returns: 0 on success, non-zero on error.
 */

int do_sscanf(const(char) * src, const(char) * format, const(char) * args, DatasetRaw * dset, int * n_items)
{
    const(char) * p = src;
    const(char) * q = format;
    const(char) * r = args;
    int nscan = 0;
    int err = 0;

    gretl_error_clear();

    while (*p && *q && !err) {
    if (isspace(*q)) {
        while (isspace(*p)) p++;
        while (isspace(*q)) q++;
    }
    if (*p == '\t' && *q == '\\' && *(q+1) == 't') {
        p++;
        q += 2;
    } else if (*p == '\n' && *q == '\\' && *(q+1) == 'n') {
        p++;
        q += 2;
    } else if (*p == *q) {
        p++;
        q++;
    } else if (*q == '%') {
        err = scan_arg(&p, &q, &r, dset, &nscan);
    } else {
        break;
    }
    }

    if (!err) {
    *n_items = nscan;
    }

    return err;
}

/* Apparatus to support the command form of printf: the @quoted
   argument indicates whether @parm was in quotes on input.
*/

int do_printf_command(const(char) * parm, const(char) * args, DatasetRaw * dset, Prn * prn, int quoted)
{
    const(char) * format;

    if (quoted) {
    format = parm;
    } else {
    /* @parm should be a string variable */
    format = get_string_by_name(parm);
    }

    if (format == null) {
    return GretlError.invarg;
    } else {
    return do_printf(format, args, dset, prn, null);
    }
}

/* The originating command is of form:

     genr markers = format.ptr, arg1,...

   We're assuming that we're just getting the part after
   the equals sign here, in the variable @s.
*/

int generate_obs_markers(const(char) * s, DatasetRaw * dset)
{
    char[16] format = 0;
    char[32] args = 0;
    Prn * prn = null;
    int n, t, err = 0;

    n = sscanf(s, "\"%15[^\"]\", %31[^\r\n]".ptr, format.ptr, args.ptr);
    if (n != 2) {
    return GretlError.parse;
    }

    prn = gretl_print_new(PrnType.buffer, &err);
    if (err) {
    return err;
    }

    if (dset.S == null) {
    err = dataset_allocate_obs_markers(dset);
    }

    if (!err) {
    const(char) * buf;

    for (t=0; t<dset.n && !err; t++) {
        gretl_print_reset_buffer(prn);
        err = real_do_printf(format.ptr, args.ptr, dset, prn, null, t);
        if (!err) {
        buf = gretl_print_get_buffer(prn);
        dset.S[t][0] = '\0';
        strncat(dset.S[t], buf, OBSLEN - 1);
        }
    }
    }

    gretl_print_destroy(prn);

    return err;
}

