/* Port of gretl_errors.c */

module gretl.gretl_errors;


import core.stdc.errno : errno;
import core.stdc.stdarg : va_end, va_list, va_start;
import core.stdc.stdio : fprintf, snprintf, vsnprintf, stderr;
import core.stdc.string : strcat, strcmp, strlen, strncat, strstr, strerror;

import gretl.gmem : g_free;
import gretl.gretl_errorsh : GretlError, GretlWarningCodes;
import gretl.gretl_prn : Prn, pprintf, print_redirection_filename, print_redirection_level;
import gretl.gstrfuncs : g_strdup;
import gretl.libgretlh;
import gretl.libset : gretl_in_batch_mode, gretl_warnings_on;
import gretl.strutils : gretl_strdup;

/* ---- declarations of not-yet-ported items used below ----
   Replace each with an import when its defining file is ported. */

extern (C) {
    void current_function_info(const(char) ** fname, const(char) ** pname);
    int gretl_compiling_function();
    int gretl_compiling_loop();
    int gretl_function_depth();
}

private enum ERRLEN = 2048;

private __gshared char[ERRLEN] gretl_errmsg = 0;
private __gshared char[ERRLEN] gretl_warnmsg = 0;

private __gshared int gretl_errno;
private __gshared int gretl_warnnum;

private immutable string[] gretl_error_messages = [
    null,  
    null,  
    N_("Data error"),  /* E_DATA = 2 */
    N_("Exact or near collinearity encountered"),  /* E_SINGULAR */
    N_("Insufficient degrees of freedom for regression"),  /* E_DF */
    N_("Dependent variable is all zeros, aborting regression"),  /* E_ZERO */
    N_("Total sum of squares was not positive"),  /* E_TSS */
    N_("Sum of squared residuals negative!"),  /* E_ESS */
    N_("Sorry, command not available for this estimator"),  /* E_NOTIMP */
    N_("Unspecified error -- FIXME"),  /* E_UNSPEC */
    N_("This command won't work with the current periodicity"),  /* E_PDWRONG */
    N_("Error attempting to open file"),  /* E_FOPEN */
    N_("Out of memory!"),  /* E_ALLOC */
    N_("No formula supplied in genr"),  /* E_EQN */
    N_("Unknown variable name in command"),  /* E_UNKVAR */
    N_("Command has insufficient arguments"),  /* E_ARGS */
    N_("This command is implemented only for OLS models"),  /* E_OLSONLY */
    N_("Invalid argument"),  /* E_INVARG */
    N_("Syntax error"),  /* E_PARSE */
    N_("No independent variables left after omissions"),  /* E_NOVARS */
    N_("No independent variables were omitted"),  /* E_NOOMIT */
    N_("No new independent variables were added"),  /* E_NOADD */
    N_("One or more \"added\" vars were already present"),  /* E_ADDDUP */
    N_("Error generating logarithms"),  /* E_LOGS */
    N_("Error generating squares"),  /* E_SQUARES */
    N_("Error generating lagged variables"),  /* E_LAGS */
    N_("Attempting to take square root of negative number"),  /* E_SQRT */
    N_("Excessive exponent in genr formula"),  /* E_HIGH */
    N_("Need valid starting and ending observations"),  /* E_OBS */
    N_("You must include a constant in this sort of model"),  /* E_NOCONST */
    N_("The statistic you requested is not available"),  /* E_BADSTAT */
    N_("Missing sub-sample information; can't merge data"),  /* E_NOMERGE */
    N_("The convergence criterion was not met"),  /* E_NOCONV */
    N_("The operation was canceled"),  /* E_CANCEL */
    N_("Missing values encountered"),  /* E_MISSDATA */
    N_("Not a Number in calculation"),  /* E_NAN */
    N_("Matrices not conformable for operation"),  /* E_NONCONF */
    N_("Data types not conformable for operation"),  /* E_TYPES */
    N_("Incompatible options"),  /* E_BADOPT */
    N_("The restrictions do not identify the parameters"),  /* E_NOIDENT */
    N_("External command failed"),  /* E_EXTERNAL */
    N_("Maximum length of command line (65536 bytes) exceeded"),  /* E_TOOLONG */
    N_("No dataset is in place"),  /* E_NODATA */
    N_("Matrix is not positive definite"),  /* E_NOTPD */
    N_("Failed to calculate Jacobian"),  /* E_JACOBIAN */
    N_("Insufficient observations for this operation"),  /* E_TOOFEW */
    N_("You cannot define a function within a function"),  /* E_FNEST */
    N_("Out of bounds error"),  /* E_BOUNDS */
    N_("Error executing function"),  /* E_FUNCERR */
    N_("Execution aborted by request"),  /* E_STOP */
    N_("'catch' cannot be used in this context"),  /* E_BADCATCH */
    N_("complex arguments/operands not supported"),  /* E_CMPLX */
    N_("mixed complex/real arguments not supported"),  /* E_MIXED */
    N_("Dependencies not met"),  /* E_DEPENDS */
    N_("Matrix argument contains non-finite values"),  /* E_NONFIN */
    N_("This command requires time-series data"),  /* E_TSONLY */
    null,  /* E_DB_DUP */
    null,  /* E_OK */
    null,  /* E_MAX */
];

private immutable string[] gretl_warning_messages = [
    null,
    N_("gradient is not close to zero"),                 /* W_GRADIENT */
    N_("generated missing values"),                      /* W_GENMISS */
    N_("generated non-finite values"),                   /* W_GENNAN */
    null                                                 /* W_MAX */
];

private const(char) * look_up_errmsg(int err)
{
    if (err > 0 && err < GretlError.max) {
        return gretl_error_messages[err].ptr;
    } else if (err == 0) {
        return "";
    } else {
        fprintf(stderr, "look_up_errmsg: out of bounds code %d\n", err);
        return "missing error message!";
    }
}

private const(char) * look_up_warnmsg(int w)
{
    if (w > 0 && w < GretlWarningCodes.max) {
        return gretl_warning_messages[w].ptr;
    } else {
        fprintf(stderr, "look_up_warnmsg: out of bounds code %d\n", w);
        return "missing warning message!";
    }
}

private __gshared int error_printed;
private __gshared int alarm_set;

/**
 * errmsg_get_with_default:
 * @err: gretl error code (see #gretl_error_codes).
 *
 * Returns: a specific error message if available,
 * otherwise a generic error message corresponding to the
 * given @err.
 */

const(char) * errmsg_get_with_default(int err)
{
    const(char) * ret = "";

    if ((err > 0 && err < GretlError.max) || err == -1) {
        if (gretl_errmsg[0] != '\0') {
            ret = gretl_errmsg.ptr;
        } else if (err > 0) {
            const(char) * deflt = look_up_errmsg(err);

            if (deflt != null) {
                /* (the C code passes @deflt through gettext) */
                ret = deflt;
            }
        }
    }

    return ret;
}

/**
 * gretl_warnmsg_get:
 *
 * Returns: the current gretl warning message, or %NULL if no
 * warning is currently in place.
 */

const(char) * gretl_warnmsg_get()
{
    const(char) * ret = null;

    if (gretl_warnnum) {
        if (gretl_warnmsg[0] != '\0') {
            ret = gretl_warnmsg.ptr;
            /* note; can't zero the message here! */
        } else {
            const(char) * deflt = look_up_warnmsg(gretl_warnnum);

            if (deflt != null) {
                ret = deflt;
            }
        }
        gretl_warnnum = 0;
    }

    return ret;
}

/**
 * errmsg:
 * @err: gretl error code (see #error_codes).
 * @prn: gretl printing struct.
 *
 * Prints to @prn a specific error message if available,
 * otherwise a generic error message corresponding to the
 * given @err.
 */

void errmsg(int err, Prn * prn)
{
    if (!error_printed && prn != null) {
        const(char) * msg = errmsg_get_with_default(err);

        if (print_redirection_level(prn) > 0) {
            /* FIXME can we get this message to appear at the
               "top level" of @prn? */
            const(char) * fname = print_redirection_filename(prn);

            if (fname != null) {
                fprintf(stderr, "error when 'outfile' (%s) active\n %s\n",
                        fname, msg);
            } else {
                fprintf(stderr, "error when 'outfile' active\n %s\n", msg);
            }
        }
        pprintf(prn, "%s\n", msg);
        error_printed = 1;
    }
}

private void print_function_info(Prn * prn)
{
    if (gretl_function_depth() > 0) {
        const(char) * fname = null;
        const(char) * pname = null;

        current_function_info(&fname, &pname);
        if (fname != null && pname != null) {
            pprintf(prn, "%s %s (%s %s):\n", _("In regard to function").ptr,
                    fname, _("package").ptr, pname);
        } else if (fname != null) {
            pprintf(prn, "%s %s:\n", _("In regard to function").ptr, fname);
        }
    }
}

/**
 * warnmsg:
 * @prn: gretl printing struct.
 *
 * If a gretl warning is set, prints a message to @prn
 * and zeros the warning signal.
 */

void warnmsg(Prn * prn)
{
    if (prn == null || gretl_warnnum == 0) {
        return;
    }

    if (!gretl_warnings_on()) {
        gretl_warnmsg[0] = '\0';
        gretl_warnnum = 0;
        return;
    }

    if (gretl_warnmsg[0] != '\0') {
        print_function_info(prn);
        pprintf(prn, "%s: %s\n", _("Warning").ptr, gretl_warnmsg.ptr);
        gretl_warnmsg[0] = '\0';
    } else {
        const(char) * s = look_up_warnmsg(gretl_warnnum);

        print_function_info(prn);
        pprintf(prn, "%s: %s\n", _("Warning").ptr, s);
    }

    gretl_warnnum = 0;
}

/**
 * gretl_errmsg_get:
 *
 * Returns: a specific error message if available,
 * otherwise an empty string.
 */

const(char) * gretl_errmsg_get()
{
    return gretl_errmsg.ptr;
}

/**
 * gretl_errmsg_set:
 * @str: an error message.
 *
 * If %gretl_errmsg is currently blank, copy the given string into
 * the message space; or if the error message is not blank but
 * sufficient space remains, append @str to the message.
 */

void gretl_errmsg_set(const(char) * str)
{
    if (alarm_set && gretl_errmsg[0] != '\0') {
        /* leave the current error message in place */
        return;
    }

    if (gretl_errmsg[0] == '\0') {
        strncat(gretl_errmsg.ptr, str, ERRLEN - 1);
    } else if (strcmp(gretl_errmsg.ptr, str)) {
        /* should we do the following? */
        int n = cast(int) strlen(gretl_errmsg.ptr);
        int m = cast(int) strlen(str);

        if (n + m + 2 < ERRLEN) {
            strcat(gretl_errmsg.ptr, "\n");
            strcat(gretl_errmsg.ptr, str);
        }
    }
}

/**
 * gretl_errmsg_append:
 * @str: an error message.
 * @err: the current error state, if any.
 *
 * Add @str to the current gretl error message, starting a
 * new line, if space permits.
 */

void gretl_errmsg_append(const(char) * str, int err)
{
    const(char) * s = null;
    char * tmp = null;

    if (gretl_errmsg[0]) {
        s = tmp = g_strdup(gretl_errmsg.ptr);
    } else if (err > 0 && err < GretlError.max) {
        s = look_up_errmsg(err);
    }
    if (s == null) {
        snprintf(gretl_errmsg.ptr, ERRLEN, "%s", str);
    } else {
        snprintf(gretl_errmsg.ptr, ERRLEN, "%s\n%s", s, str);
    }

    g_free(tmp);
}

/**
 * gretl_errmsg_prepend:
 * @str: an error message.
 * @err: the current error state, if any.
 *
 * Prepend @str to the current gretl error message, ending with a
 * newline, if space permits.
 */

void gretl_errmsg_prepend(const(char) * str, int err)
{
    const(char) * s = null;
    char * tmp = null;

    if (gretl_errmsg[0]) {
        s = tmp = g_strdup(gretl_errmsg.ptr);
    } else if (err > 0 && err < GretlError.max) {
        s = look_up_errmsg(err);
    }
    if (s == null) {
        snprintf(gretl_errmsg.ptr, ERRLEN, "%s", str);
    } else if (str[strlen(str)-1] != '\n') {
        snprintf(gretl_errmsg.ptr, ERRLEN, "%s\n%s", str, s);
    } else {
        snprintf(gretl_errmsg.ptr, ERRLEN, "%s%s", str, s);
    }

    g_free(tmp);
}

/**
 * gretl_errmsg_ensure:
 * @str: an error message.
 *
 * If %gretl_errmsg is currently blank, copy the given string into
 * the message space.
 */

void gretl_errmsg_ensure(const(char) * str)
{
    if (gretl_errmsg[0] == '\0') {
        strncat(gretl_errmsg.ptr, str, ERRLEN - 1);
    }
}

/**
 * gretl_warnmsg_set:
 * @str: a warning message.
 *
 * Copy the given string into the warning message space.
 */

void gretl_warnmsg_set(const(char) * str)
{
    gretl_warnmsg[0] = '\0';
    strncat(gretl_warnmsg.ptr, str, ERRLEN - 1);
    gretl_warnnum = GretlWarningCodes.max;
}

/**
 * gretl_errmsg_sprintf:
 * @fmt: format string.
 * @...: arguments, as to sprintf.
 *
 * Append a formatted message to the current gretl
 * error message.
 */

extern (C) void gretl_errmsg_sprintf(const(char) * fmt, ...)
{
    if (gretl_errmsg[0] == '\0') {
        va_list ap;

        va_start(ap, fmt);
        vsnprintf(gretl_errmsg.ptr, ERRLEN, fmt, ap);
        va_end(ap);
    } else if (strstr(gretl_errmsg.ptr, "*** error in fun") &&
               strstr(fmt, "*** error in fun")) {
        /* don't print more than one "error in function"
           message, as this gets confusing
        */
        ;
    } else {
        /* find the number of characters left */
        int len0 = cast(int) strlen(gretl_errmsg.ptr);
        int n = ERRLEN - len0 - 2;

        if (n > 31) {
            char[ERRLEN] tmp;
            va_list ap;

            tmp[0] = '\0';
            va_start(ap, fmt);
            vsnprintf(tmp.ptr, n, fmt, ap);
            va_end(ap);

            if (gretl_errmsg[len0 - 1] != '\n') {
                strcat(gretl_errmsg.ptr, "\n");
            }
            strcat(gretl_errmsg.ptr, tmp.ptr);
        }
    }
}

extern (C) void gretl_errmsg_sprintf_replace(const(char) * fmt, ...)
{
    va_list ap;

    gretl_errmsg[0] = '\0';
    va_start(ap, fmt);
    vsnprintf(gretl_errmsg.ptr, ERRLEN, fmt, ap);
    va_end(ap);
}

/**
 * gretl_warnmsg_sprintf:
 * @fmt: format string.
 * @...: arguments, as to sprintf.
 *
 * Write a formatted message to the current gretl
 * warning message space.
 */

extern (C) void gretl_warnmsg_sprintf(const(char) * fmt, ...)
{
    va_list ap;

    gretl_warnmsg[0] = '\0';

    va_start(ap, fmt);
    vsnprintf(gretl_warnmsg.ptr, ERRLEN, fmt, ap);
    va_end(ap);

    gretl_warnnum = GretlWarningCodes.max;
}

char * gretl_strerror(int errnum)
{
    return strerror(errnum);
}

/**
 * gretl_errmsg_set_from_errno:
 * @s: string to prepend to error message, or %NULL.
 *
 * If %gretl_errmsg is currently blank, copy the string
 * returned by %strerror into the message space; or if the
 * error message is not blank but sufficient space remains,
 * append the new error info to the message.
 */

void gretl_errmsg_set_from_errno(const(char) * s, int errnum)
{
    char * msg = null;

    if (errnum) {
        msg = gretl_strerror(errnum);
        errno = 0;
    }

    if (msg != null) {
        if (s != null) {
            gretl_errmsg_sprintf("%s: %s", s, msg);
        } else {
            gretl_errmsg_set(msg);
        }
    }
}

/**
 * gretl_error_clear:
 *
 * Blank out any previously recorded error message.
 */

int gretl_error_clear()
{
    if (!alarm_set) {
        gretl_errmsg[0] = '\0';
    }
    error_printed = 0;
    errno = 0;

    return 0;
}

/**
 * gretl_errmsg_is_set:
 *
 * Returns: 1 if the gretl error message is currently
 * set (not blank), otherwise 0.
 */

int gretl_errmsg_is_set()
{
    return (gretl_errmsg[0] != '\0');
}

/**
 * maybe_save_gretl_errmsg:
 *
 * Returns: An allocated copy of the current gretl error
 * message, or NULL if @err is zero or the current message
 * is blank.
 */

char * maybe_save_gretl_errmsg(int err)
{
    if (err && gretl_errmsg[0] != '\0') {
        return gretl_strdup(gretl_errmsg.ptr);
    } else {
        return null;
    }
}

/* setting the "alarm" prevents gretl_errmsg from being
   overwritten
*/

void set_gretl_alarm(int val)
{
    alarm_set = val;
}

private void real_set_gretl_errno(int err)
{
    gretl_errno = err;
}

void set_gretl_errno(int err)
{
    real_set_gretl_errno(err);
}

void set_gretl_warning(int w)
{
    gretl_warnnum = w;
}

int get_gretl_errno()
{
    int err = gretl_errno;

    real_set_gretl_errno(0);
    return err;
}

int check_gretl_errno()
{
    return gretl_errno;
}

int check_gretl_warning()
{
    return gretl_warnnum;
}

int gretl_error_is_fatal()
{
    if (gretl_compiling_function()) {
        return 1;
    } else if (gretl_compiling_loop()) {
        return 1;
    } else {
        return gretl_in_batch_mode();
    }
}

int invalid_field_error(const(char) * s)
{
    gretl_errmsg_sprintf(_("field '%s' in command is invalid").ptr, s);
    return GretlError.data;
}
