/* Port of lib/src/varprint.c: printing of VAR and VECM results (coefficient tables, information criteria, impulse
   responses, forecast-error variance decompositions, Johansen output). The header varprint.h has only prototypes. */

module gretl.varprint;


import core.stdc.ctype : isdigit, isspace;
import core.stdc.errno : errno, ERANGE;
import core.stdc.math : ceil, exp, fabs, floor, fmod, isinf, isnan, log, log10, pow, sqrt;
import core.stdc.stdio : fputs, fprintf, printf, puts, sprintf, snprintf, sscanf, stderr;
import core.stdc.stdlib : abs, atoi, calloc, free, getenv, malloc, realloc;
import core.stdc.string : memcpy, memset, strcmp, strcpy, strlen, strchr, strstr, strcat, strncat, strncmp, strspn, strcspn, memmove, strrchr;

import std.algorithm : max, min;

import gretl.gretl_errorsh : GretlError;
import gretl.gretl_matrix;
import gretl.gretl_matrixh;
import gretl.johansenh : JohansenCode, JohansenInfo, jcode, jrank, auto_restr;
import gretl.libgretlh;
import gretl.missingh : NADBL, na;
import gretl.gretl_restricth : GretlRestriction;
import gretl.varh;
import gretl.xmlh;


private alias XUC = const(xmlChar) *;

import gretl.dataio : ntolabel;
import gretl.gmem : g_free;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_prn : Prn, gretl_print_toggle_doc_flag, gretl_prn_newline, plain_format, pprintf, pputc, pputs, rtf_format, tex_format;
import gretl.gstrfuncs : g_strdup_printf, g_utf8_strlen;
import gretl.modelprint : printmodel;
import gretl.printout : bufspace, max_namelen_in_list, maybe_trim_varname, print_contemp_covariance_matrix;
import gretl.printouth : SHORTLEN;
import gretl.pvalues : chisq_cdf_comp, snedecor_cdf_comp;
import gretl.var : copy_north_west, default_VAR_horizon, gretl_VAR_get_fcast_decomp, gretl_VECM_n_alpha, gretl_VECM_n_beta, reorder_responses;
import gretl.gretl_restrict : print_restriction_from_matrices;

/* unported dependencies: per-file declarations from the 2026c headers */
extern (C) {
    void tex_print_VAR_ll_stats(GretlVar * var, Prn * prn);
    void tex_print_VECM_coint_eqns(GretlVar * vecm, const(DatasetRaw) * dset, Prn * prn);
    void tex_print_VECM_omega(GretlVar * vecm, const(DatasetRaw) * dset, Prn * prn);
    void tex_print_double(double x, Prn * prn);
    char * tex_escape(char * targ, const(char) * src);
}

/**
 * gretl_VAR_print_sigma:
 * @var: pointer to gretl VAR structure.
 * @prn: printing object.
 *
 * Prints to @prn the contemporaneous (cross-equation) variance
 * matrix for @var.
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_VAR_print_sigma(GretlVar * var, Prn * prn)
{
    int err = 0;

    if (var.S == null) {
    err = 1;
    } else {
    print_contemp_covariance_matrix(var.S, var.ldet, prn);
    }

    return err;
}

/* printing of impulse responses and variance decompositions */

enum IRF_ROW_MAX = 4;
enum IRF_WIDTH = 13;

enum VDC_ROW_MAX = 5;
enum VDC_WIDTH = 11;

enum VarPrint {
    irf,
    vdc
}

void VAR_RTF_row_spec(int ncols, Prn * prn)
{
    int lcol = 800, colwid = 1600;
    int i, cellx = lcol;

    pputs(prn, "{\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262");
    for (i=0; i<ncols; i++) {
    pprintf(prn, "\\cellx%d".ptr, cellx);
    cellx += colwid;
    }
    pputc(prn, '\n');
}

void VAR_info_header_block(int code, int v, int block, const(DatasetRaw) * dset, Prn * prn)
{
    int tex = tex_format(prn);
    int rtf = rtf_format(prn);
    char[48] vname;

    if (tex) {
    pputs(prn, "\\vspace{1em}\n\n");
    if (code == VarPrint.irf) {
        pprintf(prn, _("Responses to a one-standard error shock in %s").ptr,
            tex_escape(vname.ptr, dset.varname[v]));
    } else {
        pprintf(prn, _("Decomposition of variance for %s").ptr,
            tex_escape(vname.ptr, dset.varname[v]));
    }
    if (block == 0) {
        pputs(prn, "\n\n");
    } else {
        pprintf(prn, " (%s)\n\n".ptr, _("continued").ptr);
    }
    pprintf(prn, "\\vspace{1em}\n\n\\begin{longtable}{%s}\n".ptr,
        (code == VarPrint.irf)? "rrrrr".ptr : "rrrrrr".ptr);
    } else if (rtf) {
    pputs(prn, "\\par\n\n");
    if (code == VarPrint.irf) {
        pprintf(prn, _("Responses to a one-standard error shock in %s").ptr,
            dset.varname[v]);
    } else {
        pprintf(prn, _("Decomposition of variance for %s").ptr,
            dset.varname[v]);
    }
    if (block == 0) {
        pputs(prn, "\\par\n\n");
    } else {
        pprintf(prn, " (%s)\\par\n\n".ptr, _("continued").ptr);
    }
    /* FIXME */
    VAR_RTF_row_spec((code == VarPrint.irf)? IRF_ROW_MAX : VDC_ROW_MAX, prn);
    } else {
    if (code == VarPrint.irf) {
        pprintf(prn, _("Responses to a one-standard error shock in %s").ptr,
            dset.varname[v]);
    } else {
        pprintf(prn, _("Decomposition of variance for %s").ptr,
            dset.varname[v]);
    }
    if (block == 0) {
        pputs(prn, "\n\n");
    } else {
        pprintf(prn, " (%s)\n\n".ptr, _("continued").ptr);
    }
    }

    /* first column: period number header */
    if (tex) {
    pprintf(prn, "%s & ".ptr, _("period").ptr);
    } else if (rtf) {
    pprintf(prn, "\\intbl \\qc %s\\cell ".ptr, _("period").ptr);
    } else {
    pputs(prn, _("period").ptr);
    }
}

void VAR_info_print_vname(int i, int v, int endrow, int width, const(DatasetRaw) * dset, Prn * prn)
{
    int tex = tex_format(prn);
    int rtf = rtf_format(prn);
    char[32] vname;

    if (tex) {
    pprintf(prn, " %s ".ptr, tex_escape(vname.ptr, dset.varname[v]));
    if (endrow) {
        pputs(prn, "\\\\");
    } else {
        pputs(prn, "& ");
    }
    } else if (rtf) {
    pprintf(prn, "\\qc %s\\cell".ptr, dset.varname[v]);
    if (endrow) {
        pputs(prn, " \\intbl \\row");
    }
    } else {
    pprintf(prn, "%*s".ptr, width, dset.varname[v]);
    }
}

void VAR_info_print_period(int t, Prn * prn)
{
    if (tex_format(prn)) {
    pprintf(prn, "%d & ".ptr, t);
    } else if (rtf_format(prn)) {
    pprintf(prn, "\\intbl \\qc %d\\cell ".ptr, t);
    } else {
    pprintf(prn, " %3d  ".ptr, t);
    }
}

void VAR_info_end_row(Prn * prn)
{
    if (tex_format(prn)) {
    pputs(prn, "\\\\\n");
    } else if (rtf_format(prn)) {
    pputs(prn, "\\intbl \\row\n");
    } else {
    pputc(prn, '\n');
    }
}

void VAR_info_end_table(Prn * prn)
{
    if (tex_format(prn)) {
    pputs(prn, "\\end{longtable}\n\n");
    } else if (rtf_format(prn)) {
    pputs(prn, "}\n");
    } else {
    pputc(prn, '\n');
    }
}

int varprint_namelen(GretlVar * var, const(DatasetRaw) * dset, int rmax, int block)
{
    int len, maxlen = 0;
    int i, k, v;

    for (i=0; i<rmax; i++) {
    k = rmax * block + i - 1;
    if (k < 0) {
        continue;
    }
    if (k >= var.neqns) {
        break;
    }
    v = var.ylist[k+1];
    len = cast(int) (strlen(dset.varname[v]));
    if (len > maxlen) {
        maxlen = len;
    }
    }

    return maxlen;
}

/**
 * gretl_VAR_print_impulse_response:
 * @var: pointer to VAR struct.
 * @shock: index number of the "shock" variable.
 * @periods: number of periods over which to print response.
 * @dset: dataset information.
 * @prn: gretl printing object.
 *
 * Prints to @prn the estimated responses of the endogenous
 * variables in @var to a one-standard deviation shock in
 * the specified variable: @shock is a zero-based index into
 * the equations of the VAR so for example if @shock = 1,
 * the responses are to a shock in the second endogenous
 * variable in the VAR specification.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_VAR_print_impulse_response(GretlVar * var, int shock, int periods, const(DatasetRaw) * dset, Prn * prn)
{
    GretlMatrixRaw * rtmp = null, ctmp = null;
    GretlMatrixRaw * C = var.C;
    int rows = var.neqns * levels_order(var);
    int block, blockmax;
    int tex = tex_format(prn);
    int rtf = rtf_format(prn);
    int i, t, vsrc, err = 0;

    if (prn == null) {
    return 0;
    }

    if (shock >= var.neqns) {
    fprintf(stderr, "Shock variable out of bounds\n".ptr);
    return 1;
    }

    if (var.ord != null) {
    C = reorder_responses(var, &err);
    if (err) {
        return err;
    }
    }

    rtmp = gretl_matrix_alloc(rows, var.neqns);
    ctmp = gretl_matrix_alloc(rows, var.neqns);

    if (rtmp == null || ctmp == null) {
    err = GretlError.alloc;
    goto bailout;
    }

    vsrc = var.ylist[shock + 1];

    blockmax = var.neqns / IRF_ROW_MAX;
    if (var.neqns % IRF_ROW_MAX) {
    blockmax++;
    }

    for (block=0; block<blockmax && !err; block++) {
    int k, vtarg, endrow;
    int namelen, width;
    double r;

    VAR_info_header_block(VarPrint.irf, vsrc, block, dset, prn);

    namelen = varprint_namelen(var, dset, IRF_ROW_MAX, block);
    width = (namelen < IRF_WIDTH - 1)? IRF_WIDTH : namelen + 1;

    for (i=0; i<IRF_ROW_MAX; i++) {
        k = IRF_ROW_MAX * block + i;
        if (k >= var.neqns) {
        break;
        }
        vtarg = var.ylist[k+1];
        endrow = !(i < IRF_ROW_MAX - 1 && k < var.neqns - 1);
        VAR_info_print_vname(i, vtarg, endrow, width, dset, prn);
    }

    if (tex || rtf) {
        pputc(prn, '\n');
    } else {
        pputs(prn, "\n\n");
    }

    for (t=0; t<periods && !err; t++) {
        VAR_info_print_period(t + 1, prn);
        if (t == 0) {
        /* calculate initial estimated responses */
        copy_north_west(rtmp, C, 0);
        } else {
        /* calculate further estimated responses */
        err = gretl_matrix_multiply(var.A, rtmp, ctmp);
        gretl_matrix_copy_values(rtmp, ctmp);
        }

        if (err) break;

        for (i=0; i<IRF_ROW_MAX; i++) {
        k = IRF_ROW_MAX * block + i;
        if (k >= var.neqns) {
            break;
        }
        r = gretl_matrix_get(rtmp, k, shock);
        if (tex) {
            tex_print_double(r, prn);
            if (i < IRF_ROW_MAX - 1 && k < var.neqns - 1) {
            pputs(prn, " & ");
            }
        } else if (rtf) {
            pprintf(prn, "\\qc %.5g\\cell ".ptr, r);
        } else {
            if (i == 0) pputc(prn, ' ');
            pprintf(prn, "%#*.5g ".ptr, width - 1, r);
        }
        }

        VAR_info_end_row(prn);
    }

    VAR_info_end_table(prn);
    }

 bailout:

    if (rtmp != null) gretl_matrix_free(rtmp);
    if (ctmp != null) gretl_matrix_free(ctmp);

    if (C != var.C) {
    gretl_matrix_free(C);
    }

    return err;
}

int gretl_VAR_print_all_impulse_responses(GretlVar * var, const(DatasetRaw) * dset, int horizon, Prn * prn)
{
    int i, err = 0;

    if (horizon <= 0) {
    horizon = default_VAR_horizon(dset);
    }

    if (rtf_format(prn)) {
    pputs(prn, "{\\rtf1\\par\n\\qc ");
    }

    for (i=0; i<var.neqns && !err; i++) {
    err = gretl_VAR_print_impulse_response(var, i, horizon, dset,
                           prn);
    }

    if (rtf_format(prn)) {
    pputs(prn, "}\n");
    }

    return err;
}

/**
 * gretl_VAR_print_fcast_decomp:
 * @var: pointer to VAR struct.
 * @targ:
 * @periods: number of periods over which to print decomposition.
 * @dset: dataset information.
 * @prn: gretl printing struct.
 *
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_VAR_print_fcast_decomp(GretlVar * var, int targ, int periods, const(DatasetRaw) * dset, Prn * prn)
{
    int vtarg;
    GretlMatrixRaw * vd = null;
    int block, blockmax;
    int tex = tex_format(prn);
    int rtf = rtf_format(prn);
    int i, t, err = 0;

    if (prn == null) {
    return 0;
    }

    if (targ >= var.neqns) {
    fprintf(stderr, "Target variable out of bounds\n".ptr);
    return 1;
    }

    vd = gretl_VAR_get_fcast_decomp(var, targ, periods, &err);
    if (err) {
    pprintf(prn, _("Forecast decomposition failed\n").ptr);
    return err;
    }

    vtarg = var.ylist[targ + 1];

    blockmax = (var.neqns + 1) / VDC_ROW_MAX;
    if ((var.neqns + 1) % VDC_ROW_MAX) {
    blockmax++;
    }

    for (block=0; block<blockmax; block++) {
    int k, vsrc, endrow;
    int namelen, width;
    double r;

    VAR_info_header_block(VarPrint.vdc, vtarg, block, dset, prn);

    namelen = varprint_namelen(var, dset, VDC_ROW_MAX, block);
    width = (namelen < VDC_WIDTH - 1)? VDC_WIDTH : namelen + 1;

    for (i=0; i<VDC_ROW_MAX; i++) {
        k = VDC_ROW_MAX * block + i - 1;
        if (k < 0) {
        if (tex) {
            pprintf(prn, " %s & ".ptr, _("std. error").ptr);
        } else if (rtf) {
            pprintf(prn, " \\qc %s\\cell ".ptr, _("std. error").ptr);
        } else {
            pprintf(prn, " %14s".ptr, _("std. error").ptr);
        }
        continue;
        }
        if (k >= var.neqns) {
        break;
        }
        vsrc = var.ylist[k+1];
        endrow = !(i < VDC_ROW_MAX - 1 && k < var.neqns - 1);
        VAR_info_print_vname(i, vsrc, endrow, width, dset, prn);
    }

    if (tex || rtf) {
        pputc(prn, '\n');
    } else {
        pputs(prn, "\n\n");
    }

    for (t=0; t<periods && !err; t++) {
        VAR_info_print_period(t + 1, prn);
        for (i=0; i<VDC_ROW_MAX; i++) {
        k = VDC_ROW_MAX * block + i - 1;
        if (k < 0) {
            /* standard error column */
            r = gretl_matrix_get(vd, t, var.neqns);
            if (tex) {
            pprintf(prn, "%g & ".ptr, r);
            } else if (rtf) {
            pprintf(prn, "\\qc %g\\cell".ptr, r);
            } else {
            pprintf(prn, " %14g ".ptr, r);
            }
            continue;
        }
        if (k >= var.neqns) {
            break;
        }
        r = gretl_matrix_get(vd, t, k);
        if (tex) {
            pprintf(prn, "$%.4f$".ptr, r);
            if (i < VDC_ROW_MAX - 1 && k < var.neqns - 1) {
            pputs(prn, " & ");
            }
        } else if (rtf) {
            pprintf(prn, "\\qc %.4f\\cell".ptr, r);
        } else {
            pprintf(prn, "%*.4f ".ptr, width - 1, r);
        }
        }

        VAR_info_end_row(prn);
    }

    VAR_info_end_table(prn);
    }

    if (vd != null) {
    gretl_matrix_free(vd);
    }

    return err;
}

int gretl_VAR_print_all_fcast_decomps(GretlVar * var, const(DatasetRaw) * dset, int horizon, Prn * prn)
{
    int i, err = 0;

    if (horizon <= 0) {
    horizon = default_VAR_horizon(dset);
    }

    if (rtf_format(prn)) {
    pputs(prn, "{\\rtf1\\par\n\\qc ");
    }

    for (i=0; i<var.neqns && !err; i++) {
    err = gretl_VAR_print_fcast_decomp(var, i, horizon, dset,
                       prn);
    }

    if (rtf_format(prn)) {
    pputs(prn, "}\n");
    }

    return err;
}

void print_Johansen_test_case(JohansenCode jcode, Prn * prn)
{
    static const(char) *[5] jcase = ["Case 1: No constant".ptr, "Case 2: Restricted constant".ptr, "Case 3: Unrestricted constant".ptr, "Case 4: Restricted trend, unrestricted constant".ptr, "Case 5: Unrestricted trend and constant".ptr];

    if (jcode <= JohansenCode.unrest_trend) {
    if (plain_format(prn) || tex_format(prn)) {
        pputs(prn, _(jcase[jcode]));
    } else {
        pputs(prn, _(jcase[jcode]));
    }
    }
}

char * vecm_beta_varname(char * vname, GretlVar * v, const(DatasetRaw) * dset, int i)
{
    const(char) * src = "";

    if (i < v.neqns) {
    src = dset.varname[v.ylist[i+1]];
    } else if (auto_restr(v) && i == v.neqns) {
    src = (jcode(v) == JohansenCode.rest_const)? "const" : "trend";
    } else if (v.rlist != null) {
    int k = i - v.ylist[0] - auto_restr(v) + 1;

    src = dset.varname[v.rlist[k]];
    }

    maybe_trim_varname(vname, src);

    return vname;
}

int max_beta_namelen(GretlVar * v, const(DatasetRaw) * dset)
{
    int r = gretl_matrix_rows(v.jinfo.Beta);
    char[32] s;
    int i, ni, n = 0;

    for (i=0; i<r; i++) {
    vecm_beta_varname(s.ptr, v, dset, i);
    ni = cast(int) (strlen(s.ptr));
    if (ni > n) {
        n = ni;
    }
    }

    return n;
}

void print_VECM_coint_eqns(GretlVar * jvar, const(DatasetRaw) * dset, Prn * prn)
{
    JohansenInfo *jv = jvar.jinfo;
    int rtf = rtf_format(prn);
    char[16] namefmt;
    char[16] s;
    char[32] vname;
    int rows = gretl_matrix_rows(jv.Beta);
    int nwid;
    int i, j;
    double x;

    pprintf(prn, "beta (%s".ptr, _("cointegrating vectors").ptr);
    if (jv.Bse != null) {
    pprintf(prn, ", %s)".ptr, _("standard errors in parentheses").ptr);
    } else {
    pputc(prn, ')');
    }

    gretl_prn_newline(prn);
    gretl_prn_newline(prn);

    nwid = max_beta_namelen(jvar, dset) + 1;
    sprintf(namefmt.ptr, "%%-%ds".ptr, nwid);

    for (i=0; i<rows; i++) {
    vecm_beta_varname(vname.ptr, jvar, dset, i);
    if (rtf) {
        pputs(prn, vname.ptr);
    } else {
        pprintf(prn, namefmt.ptr, vname.ptr);
    }

    /* coefficients */
    for (j=0; j<jv.rank; j++) {
        x = gretl_matrix_get(jv.Beta, i, j);
        if (rtf) {
        pprintf(prn, "\t%#.5g ".ptr, x);
        } else {
        pprintf(prn, "%#12.5g ".ptr, x);
        }
    }
    gretl_prn_newline(prn);

    if (jv.Bse != null) {
        /* standard errors */
        if (rtf) {
        pputs(prn, "\t");
        } else {
        bufspace(nwid + 1, prn);
        }
        for (j=0; j<jv.rank; j++) {
        x = gretl_matrix_get(jv.Bse, i, j);
        sprintf(s.ptr, "(%#.5g)".ptr, x);
        if (rtf) {
            pprintf(prn, "\t%s".ptr, s.ptr);
        } else {
            pprintf(prn, "%12s ".ptr, s.ptr);
        }
        }
        gretl_prn_newline(prn);
    }
    }

    gretl_prn_newline(prn);

    rows = gretl_matrix_rows(jv.Alpha);

    pprintf(prn, "alpha (%s".ptr, _("adjustment vectors").ptr);
    if (jv.Ase != null) {
    pprintf(prn, ", %s)".ptr, _("standard errors in parentheses").ptr);
    } else {
    pputc(prn, ')');
    }

    gretl_prn_newline(prn);
    gretl_prn_newline(prn);

    for (i=0; i<rows; i++) {
    const(char) * src = dset.varname[jvar.ylist[i+1]];

    if (rtf) {
        pputs(prn, src);
    } else {
        maybe_trim_varname(vname.ptr, src);
        pprintf(prn, namefmt.ptr, vname.ptr);
    }

    for (j=0; j<jv.rank; j++) {
        x = gretl_matrix_get(jv.Alpha, i, j);
        if (rtf) {
        pprintf(prn, "\t%#.5g ".ptr, x);
        } else {
        pprintf(prn, "%#12.5g ".ptr, x);
        }
    }
    gretl_prn_newline(prn);

    if (jv.Ase != null) {
        if (rtf) {
        pputs(prn, "\t");
        } else {
        bufspace(nwid + 1, prn);
        }
        for (j=0; j<jv.rank; j++) {
        x = gretl_matrix_get(jv.Ase, i, j);
        sprintf(s.ptr, "(%#.5g)".ptr, x);
        if (rtf) {
            pprintf(prn, "\t%s".ptr, s.ptr);
        } else {
            pprintf(prn, "%12s ".ptr, s.ptr);
        }
        }
        gretl_prn_newline(prn);
    }
    }

    gretl_prn_newline(prn);
}

void print_VECM_omega(GretlVar * jvar, const(DatasetRaw) * dset, Prn * prn)
{
    int rtf = rtf_format(prn);
    int * list = jvar.ylist;
    const(char) * src;
    char[32] vname = 0;
    int vwidth = 13;
    int w0, wi = 12;
    int i, j;

    pprintf(prn, "%s:\n".ptr, _("Cross-equation covariance matrix").ptr);
    gretl_prn_newline(prn);

    w0 = max_namelen_in_list(list, dset) + 1;
    if (w0 > vwidth) {
    vwidth = w0;
    }

    /* top row: names of Y variables */

    for (i=0; i<jvar.neqns; i++) {
    src = dset.varname[list[i+1]];
    if (plain_format(prn)) {
        maybe_trim_varname(vname.ptr, src);
    }
    if (i == 0) {
        if (rtf) {
        pprintf(prn, "\t\t%s".ptr, src);
        } else {
        wi = cast(int) (strlen(vname.ptr));
        if (wi < vwidth) wi = vwidth;
        pprintf(prn, "%*s".ptr, w0 + wi, vname.ptr);
        }
    } else {
        if (rtf) {
        pprintf(prn, "\t%s".ptr, src);
        } else {
        wi = cast(int) (strlen(vname.ptr) + 1);
        if (wi < vwidth) wi = vwidth;
        pprintf(prn, "%*s".ptr, wi, vname.ptr);
        }
    }
    }

    gretl_prn_newline(prn);

    /* subsequent rows: Y name plus values */

    for (i=0; i<jvar.neqns; i++) {
    src = dset.varname[list[i+1]];
    if (plain_format(prn)) {
        maybe_trim_varname(vname.ptr, src);
    }
    if (rtf) {
        pputs(prn, src);
        if (strlen(src) < 8) {
        pputc(prn, '\t');
        }
    } else {
        pprintf(prn, "%-*s ".ptr, w0, vname.ptr);
    }
    for (j=0; j<jvar.neqns; j++) {
        if (rtf) {
        pprintf(prn, "\t%#.5g".ptr, gretl_matrix_get(jvar.S, i, j));
        } else {
        src = dset.varname[list[j+1]];
        wi = cast(int) (strlen(src));
        if (wi >= SHORTLEN) {
            wi = SHORTLEN - 1;
        }
        if (wi < vwidth - 1) wi = vwidth - 1;
        pprintf(prn, "%#*.5g ".ptr, wi, gretl_matrix_get(jvar.S, i, j));
        }
    }
    gretl_prn_newline(prn);
    }

    gretl_prn_newline(prn);

    pprintf(prn, "%s = %g".ptr, _("determinant").ptr, exp(jvar.ldet));

    gretl_prn_newline(prn);
}

/* FIXME TeX and RTF */

void vecm_print_restrictions(GretlVar * vecm, Prn * prn)
{
    if (vecm.jinfo.R != null) {
    pputs(prn, "\n\n");
    pputs(prn, _("Restrictions on beta:").ptr);
    pputc(prn, '\n');
    print_restriction_from_matrices(vecm.jinfo.R, vecm.jinfo.q,
                    'b', gretl_VECM_n_beta(vecm),
                    prn);
    pputc(prn, '\n');
    }

    if (vecm.jinfo.Ra != null) {
    if (vecm.jinfo.R == null) {
        pputs(prn, "\n\n");
    }
    pputs(prn, _("Restrictions on alpha:").ptr);
    pputc(prn, '\n');
    print_restriction_from_matrices(vecm.jinfo.Ra, vecm.jinfo.qa,
                    'a', gretl_VECM_n_alpha(vecm),
                    prn);
    pputc(prn, '\n');
    }
}

void print_LR_stat(double x, int df, Prn * prn)
{
    double pv = chisq_cdf_comp(df, x);

    if (na(pv)) {
    return;
    }

    if (tex_format(prn)) {
    pprintf(prn, "$2 (l_u - l_r) = %g$".ptr, x);
    gretl_prn_newline(prn);
    pprintf(prn, "$P(\\chi^2_{%d} > %g) = %g$".ptr, df, x, pv);
    } else if (rtf_format(prn)) {
    pprintf(prn, "2 * (lu - lr) = %g".ptr, x);
    gretl_prn_newline(prn);
    pprintf(prn, "P(%s(%d) > %g) = %g".ptr, _("Chi-square").ptr, df, x, pv);
    } else {
    pprintf(prn, "2 * (lu - lr) = %g".ptr, x);
    gretl_prn_newline(prn);
    pprintf(prn, "P(%s(%d) > %g) = %g".ptr, _("Chi-square").ptr, df, x, pv);
    }
}

enum LrKind {
    absolute,
    relative
}

void vecm_print_LR_test(GretlVar * vecm, Prn * prn, int code)
{
    double ll0, x;
    int df;

    if (code == LrKind.relative) {
    ll0 = vecm.jinfo.prior_ll;
    df = vecm.jinfo.lrdf - vecm.jinfo.prior_df;
    gretl_prn_newline(prn);
    pputs(prn, _("Relative to prior restriction").ptr);
    pputc(prn, ':');
    gretl_prn_newline(prn);
    } else {
    ll0 = vecm.jinfo.ll0;
    df = vecm.jinfo.lrdf;
    }

    x = 2.0 * (ll0 - vecm.ll);

    if (tex_format(prn)) {
    pprintf(prn, _("Unrestricted loglikelihood $(l_u) = %.8g$").ptr, ll0);
    gretl_prn_newline(prn);
    pprintf(prn, _("Restricted loglikelihood $(l_r) = %.8g$").ptr, vecm.ll);
    } else {
    pprintf(prn, _("Unrestricted loglikelihood (lu) = %.8g").ptr, ll0);
    gretl_prn_newline(prn);
    pprintf(prn, _("Restricted loglikelihood (lr) = %.8g").ptr, vecm.ll);
    }

    gretl_prn_newline(prn);
    print_LR_stat(x, df, prn);
    gretl_prn_newline(prn);
}

void print_vecm_header_info(GretlVar * vecm, int * lldone, Prn * prn)
{
    JohansenInfo *J = vecm.jinfo;
    gretl_prn_newline(prn);

    if (vecm.jinfo == null) {
    return;
    }

    pprintf(prn, "%s = %d".ptr, _("Cointegration rank").ptr, jrank(vecm));
    gretl_prn_newline(prn);
    print_Johansen_test_case(jcode(vecm), prn);

    if (J.R != null || J.Ra != null) {
    vecm_print_restrictions(vecm, prn);
    if (!na(J.ll0) && J.lrdf > 0) {
        vecm_print_LR_test(vecm, prn, LrKind.absolute);
        *lldone = 1;
    }
    if (!na(J.prior_ll) && J.prior_df > 0) {
        vecm_print_LR_test(vecm, prn, LrKind.relative);
        *lldone = 1;
    }
    if (!*lldone) {
        pputc(prn, '\n');
    }
    } else {
    pputc(prn, '\n');
    }
}

void VAR_print_LB_stat(GretlVar * var, Prn * prn)
{
    int k = levels_order(var);
    int df = var.neqns * var.neqns * (var.LBs - k);
    double pv = chisq_cdf_comp(df, var.LB);

    if (df <= 0) {
    return;
    }

    if (tex_format(prn)) {
    pprintf(prn, "\\noindent\n%s: LB(%d) = %g, %s = %d [%.4f]\\par\n".ptr,
        _("Portmanteau test").ptr, var.LBs, var.LB,
        _("df").ptr, df, pv);
    } else if (rtf_format(prn)) {
    pprintf(prn, "%s: LB(%d) = %g, %s = %d [%.4f]\\par\n".ptr,
        _("Portmanteau test").ptr, var.LBs, var.LB,
        _("df").ptr, df, pv);
    } else {
    pprintf(prn, "%s: LB(%d) = %g, %s = %d [%.4f]\n".ptr,
        _("Portmanteau test").ptr, var.LBs, var.LB,
        _("df").ptr, df, pv);
    }
}

int Ivals_ok(GretlVar * var)
{
    return var.Ivals != null && !na(var.Ivals[0])
    && !na(var.Ivals[1]) && !na(var.Ivals[2]);
}

int max_Ftest_label_len(GretlVar * var, const(DatasetRaw) * dset, int maxlag)
{
    char * tmp = null;
    int len, len1, len2;
    int maxnamelen = 0;
    int i, v;

    for (i=0; i<var.neqns; i++) {
    v = var.models[i].list[1];
    len = cast(int) (strlen(dset.varname[v]));
    if (len > maxnamelen) {
        maxnamelen = len;
    }
    }

    if (maxnamelen >= SHORTLEN) {
    maxnamelen = SHORTLEN - 1;
    }

    tmp = g_strdup_printf(_("All lags of %s").ptr, "x".ptr);
    len1 = cast(int) (g_utf8_strlen(tmp, -1) + maxnamelen - 1);
    g_free(tmp);

    tmp = g_strdup_printf(_("All vars, lag %d").ptr, maxlag);
    len2 = cast(int) (g_utf8_strlen(tmp, -1));
    g_free(tmp);

    len = (len1 > len2)? len1 : len2;

    return (len > 25)? (len + 1) : 25;
}

/**
 * gretl_VAR_print:
 * @var: pointer to VAR struct.
 * @dset: dataset information.
 * @opt: if includes %OPT_I, include impulse responses; if
 * includes %OPT_F, include forecast variance decompositions;
 * if includes %OPT_Q, don't print individual regressions.
 * @prn: pointer to printing struct.
 *
 * Prints the models in @var, along with relevant F-tests and
 * possibly impulse responses and variance decompositions.
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_VAR_print(GretlVar * var, const(DatasetRaw) * dset, Gretlopt opt, Prn * prn)
{
    char[OBSLEN] startdate;
    char[OBSLEN] enddate;
    char * label = null;
    int vecm = (var.ci == GretlCmdIndex.vecm);
    int dfd = var.models[0].dfd;
    int tex = tex_format(prn);
    int rtf = rtf_format(prn);
    int quiet = (opt & Gretlopt.OPT_Q);
    int nlags, maxlag, nextlag;
    int llen, fwidth = 0;
    int lldone = 0;
    double pv;
    int i, j, k, v;

    if (prn == null || (opt & Gretlopt.OPT_S)) {
    return 0;
    }

    nlags = var_n_lags(var);
    maxlag = var_max_lag(var);

    if (var.lags != null) {
    nextlag = var.lags[nlags - 1];
    } else {
    nextlag = maxlag - 1;
    }

    ntolabel(startdate.ptr, var.t1, dset);
    ntolabel(enddate.ptr, var.t2, dset);

    if (rtf) {
    pputs(prn, "{\\rtf1\\par\n\\qc ");
    }

    if (vecm) {
    label = g_strdup_printf(_("VECM system, lag order %d").ptr, levels_order(var));
    } else {
    label = g_strdup_printf(_("VAR system, lag order %d").ptr, var.order);
    }

    if (tex) {
    pputs(prn, "\\begin{center}");
    pprintf(prn, "\n%s\\\\\n".ptr, label);
    pprintf(prn, _("%s estimates, observations %s--%s ($T=%d$)").ptr,
        (vecm)? _("Maximum likelihood").ptr : _("OLS").ptr, startdate.ptr, enddate.ptr, var.T);
    if (vecm) {
        print_vecm_header_info(var, &lldone, prn);
    }
    pputs(prn, "\n\\end{center}\n");
    } else if (rtf) {
    gretl_print_toggle_doc_flag(prn);
    pprintf(prn, "\n%s\\par\n".ptr, label);
    pprintf(prn, _("%s estimates, observations %s-%s (T = %d)").ptr,
        (vecm)? _("Maximum likelihood").ptr : _("OLS").ptr, startdate.ptr, enddate.ptr, var.T);
    if (vecm) {
        print_vecm_header_info(var, &lldone, prn);
    }
    pputs(prn, "\\par\n\n");
    } else {
    pprintf(prn, "\n%s\n".ptr, label);
    pprintf(prn, _("%s estimates, observations %s-%s (T = %d)").ptr,
        (vecm)? _("Maximum likelihood").ptr : _("OLS").ptr, startdate.ptr, enddate.ptr, var.T);
    if (vecm) {
        print_vecm_header_info(var, &lldone, prn);
    }
    pputc(prn, '\n');
    }

    g_free(label);
    label = null;

    if (vecm) {
    if (tex_format(prn)) {
        tex_print_VECM_coint_eqns(var, dset, prn);
    } else {
        print_VECM_coint_eqns(var, dset, prn);
    }
    }

    if (tex) {
    tex_print_VAR_ll_stats(var, prn);
    } else if (rtf) {
    if (!lldone) {
        pprintf(prn, "%s = %.8g\\par\n".ptr, _("Log-likelihood").ptr, var.ll);
    }
    pprintf(prn, "%s = %.8g\\par\n".ptr, _("Determinant of covariance matrix").ptr,
        exp(var.ldet));
    pprintf(prn, "%s = %.4f\\par\n".ptr, _("AIC").ptr, var.AIC);
    pprintf(prn, "%s = %.4f\\par\n".ptr, _("BIC").ptr, var.BIC);
    pprintf(prn, "%s = %.4f\\par\n".ptr, _("HQC").ptr, var.HQC);
    } else {
    if (!lldone) {
        pprintf(prn, "%s = %.8g\n".ptr, _("Log-likelihood").ptr, var.ll);
    }
    pprintf(prn, "%s = %.8g\n".ptr, _("Determinant of covariance matrix").ptr, exp(var.ldet));
    pprintf(prn, "%s = %.4f\n".ptr, _("AIC").ptr, var.AIC);
    pprintf(prn, "%s = %.4f\n".ptr, _("BIC").ptr, var.BIC);
    pprintf(prn, "%s = %.4f\n".ptr, _("HQC").ptr, var.HQC);
    }

    if (var.LBs > 0 && !na(var.LB)) {
    VAR_print_LB_stat(var, prn);
    }

    if (vecm && !quiet) {
    pputc(prn, '\n');
    }

    k = 0;

    for (i=0; i<var.neqns; i++) {
    double Fval;
    char[24] Fstr;

    if (!quiet) {
        printmodel(var.models[i], dset, Gretlopt.OPT_NONE, prn);
    } else {
        if (var.ci != GretlCmdIndex.vecm) {
            v = var.models[i].list[1];
        if (tex) {
            pputs(prn, "\n\\begin{center}\n");
            pprintf(prn, "%s\\\\[1em]\n".ptr, _("Equation for ").ptr);
            pprintf(prn, "%s\\\n".ptr, dset.varname[v]);
            pputs(prn, "\n\\end{center}\n");
        } else if (rtf) {
            pprintf(prn, "\\par\n%s".ptr, _("Equation for ").ptr);
            pprintf(prn, "%s:\\par\n\n".ptr, dset.varname[v]);
        } else {
            pprintf(prn, "\n%s".ptr, _("Equation for ").ptr);
            pprintf(prn, "%s:\n".ptr, dset.varname[v]);
        }
        }
    }

    if (vecm) {
        continue;
    }

    if (var.order == 0) {
        goto skip_lag_tests;
    }

    if (tex) {
        pputs(prn, "\n\\begin{center}\n");
        pprintf(prn, "%s\\\\[1em]\n".ptr, _("F-tests of zero restrictions").ptr);
        pputs(prn, "\\begin{tabular}{lll}\n");
    } else if (rtf) {
        pprintf(prn, "%s:\\par\n\n".ptr, _("F-tests of zero restrictions").ptr);
    } else {
        pprintf(prn, "%s:\n\n".ptr, _("F-tests of zero restrictions").ptr);
    }

    if (!tex) {
        fwidth = max_Ftest_label_len(var, dset, maxlag);
    }

    for (j=0; j<var.neqns; j++) {
        Fval = var.Fvals[k];
        if (!na(Fval)) {
        char * lagstr = null;
        const(char) * vname;

        pv = snedecor_cdf_comp(nlags, dfd, Fval);
        v = (var.models[j]).list[1];
        vname = dset.varname[v];
        if (tex) {
            char[32] tname;

            pprintf(prn, _("All lags of %s").ptr, tex_escape(tname.ptr, vname));
            pputs(prn, " & ");
            pprintf(prn, "$F(%d, %d) = %g$ & ".ptr, nlags, dfd, Fval);
            pprintf(prn, "[%.4f]\\\\\n".ptr, pv);
        } else if (rtf) {
            lagstr = g_strdup_printf(_("All lags of %s").ptr, vname);
            llen = cast(int) (strlen(lagstr));
            pputs(prn, lagstr);
            bufspace(fwidth - llen, prn);
            pprintf(prn, "F(%d, %d) = %8.5g ".ptr, nlags, dfd, Fval);
            pprintf(prn, "[%.4f]\\par\n".ptr, pv);
        } else {
            char[SHORTLEN] tmp;

            maybe_trim_varname(tmp.ptr, vname);
            lagstr = g_strdup_printf(_("All lags of %s").ptr, tmp.ptr);
            llen = cast(int) (g_utf8_strlen(lagstr, -1));
            pputs(prn, lagstr);
            bufspace(fwidth - llen, prn);
            sprintf(Fstr.ptr, "F(%d, %d)".ptr, nlags, dfd);
            pprintf(prn, "%12s = %#8.5g [%.4f]\n".ptr, Fstr.ptr, Fval, pv);
        }
        g_free(lagstr);
        }
        k++;
    }

    skip_lag_tests:

    if (var.order > 0) {
        Fval = var.Fvals[k];
        if (!na(Fval)) {
        char * lagstr = null;

        pv = snedecor_cdf_comp(var.neqns, dfd, Fval);
        if (tex) {
            pprintf(prn, _("All vars, lag %d").ptr, maxlag);
            pputs(prn, " & ");
            pprintf(prn, "$F(%d, %d) = %g$ & ".ptr, var.neqns, dfd, Fval);
            pprintf(prn, "[%.4f]\\\\\n".ptr, pv);
        } else if (rtf) {
            lagstr = g_strdup_printf(_("All vars, lag %d").ptr, maxlag);
            llen = cast(int) (strlen(lagstr));
            pputs(prn, lagstr);
            bufspace(fwidth - llen, prn);
            pprintf(prn, "F(%d, %d) = %8.5g ".ptr, var.neqns, dfd, Fval);
            pprintf(prn, "[%.4f]\\par\n".ptr, pv);
        } else {
            lagstr = g_strdup_printf(_("All vars, lag %d").ptr, maxlag);
            llen = cast(int) (g_utf8_strlen(lagstr, -1));
            pputs(prn, lagstr);
            bufspace(fwidth - llen, prn);
            sprintf(Fstr.ptr, "F(%d, %d)".ptr, var.neqns, dfd);
            pprintf(prn, "%12s = %#8.5g [%.4f]\n".ptr, Fstr.ptr, Fval, pv);
        }
        g_free(lagstr);
        }
        k++;
    }

    if (tex) {
        pputs(prn, "\\end{tabular}\n\\end{center}\n\n\\clearpage\n\n");
    } else if (rtf) {
        pputs(prn, "\\par\\n\n");
    }
    }

    /* global LR test on max lag */
    if (!na(var.LR)) {
    char * h0str, h1str;
    int df = var.neqns * var.neqns;

    pputc(prn, '\n');

    h0str = g_strdup_printf(_("the longest lag is %d").ptr, nextlag);
    h1str = g_strdup_printf(_("the longest lag is %d").ptr, maxlag);

    pv = chisq_cdf_comp(df, var.LR);

    if (tex) {
        pprintf(prn, "\\noindent %s ---\\par\n".ptr, _("For the system as a whole").ptr);
        pprintf(prn, "%s: %s\\par\n".ptr, _("Null hypothesis").ptr, h0str);
        pprintf(prn, "%s: %s\\par\n".ptr, _("Alternative hypothesis").ptr, h1str);
        pprintf(prn, "%s: $\\chi^2_{%d}$ = %.3f [%.4f]\\par\n".ptr,
            _("Likelihood ratio test").ptr, df, var.LR, pv);
    } else if (rtf) {
        pprintf(prn, "\\par %s\n".ptr, _("For the system as a whole").ptr);
        pprintf(prn, "\\par %s: %s\n".ptr, _("Null hypothesis").ptr, h0str);
        pprintf(prn, "\\par %s: %s\n".ptr, _("Alternative hypothesis").ptr, h1str);
        pprintf(prn, "\\par %s: %s(%d) = %g [%.4f]\n".ptr, _("Likelihood ratio test").ptr,
            _("Chi-square").ptr, df, var.LR, pv);
    } else {
        int ordlen = (var.order > 10)? 2 : 1;

        pprintf(prn, "%s:\n\n".ptr, _("For the system as a whole").ptr);
        pprintf(prn, "  %s: %s\n".ptr, _("Null hypothesis").ptr, h0str);
        pprintf(prn, "  %s: %s\n".ptr, _("Alternative hypothesis").ptr, h1str);
        pprintf(prn, "  %s: %s(%d) = %g [%.4f]\n".ptr, _("Likelihood ratio test").ptr,
            _("Chi-square").ptr, df, var.LR, pv);
        if (Ivals_ok(var)) {
        /* Info criteria comparison */
        pprintf(prn, "\n  %s:\n".ptr, _("Comparison of information criteria").ptr);
        pputs(prn, "  ");
        pprintf(prn, _("Lag order %*d").ptr, ordlen, var.order);
        pprintf(prn, ": AIC = %#.6g, BIC = %#.6g, HQC = %#.6g\n".ptr,
            var.AIC, var.BIC, var.HQC);
        pputs(prn, "  ");
        pprintf(prn, _("Lag order %*d").ptr, ordlen, var.order - 1);
        pprintf(prn, ": AIC = %#.6g, BIC = %#.6g, HQC = %#.6g\n".ptr,
            var.Ivals[0], var.Ivals[1], var.Ivals[2]);
        }
    }
    g_free(h0str);
    g_free(h1str);
    }

    if (vecm) {
    if (tex_format(prn)) {
        tex_print_VECM_omega(var, dset, prn);
    } else {
        print_VECM_omega(var, dset, prn);
        pputc(prn, '\n');
    }
    } else if (var.order > 0) {
    pputc(prn, '\n');
    }

    if (opt & Gretlopt.OPT_I) {
    gretl_VAR_print_all_impulse_responses(var, dset, 0, prn);
    }

    if (opt & Gretlopt.OPT_F) {
    gretl_VAR_print_all_fcast_decomps(var, dset, 0, prn);
    }

    if (rtf) {
    pputs(prn, "}\n");
    }

    return 0;
}

