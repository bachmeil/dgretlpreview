/* Port of lib/src/modelprint.c: printing of estimated models in plain text, TeX, RTF and CSV formats (coefficient
   tables, the "middle table" of model statistics, tests and information criteria, per-model-type extras; the
   header modelprint.h is modelprinth.d).

   Differences from the C:
   - the anonymous enums of the file are the named enums ChiTest, AltStats, MinusType, MiddleType and MiddleKey, and
     the file-internal structs are Middletab and Printval;
   - the RTF/TeX/text row-format macros are string enums, the model-type test macros (ordered_model, tsls_model, ...)
     are functions;
   - the arrays of N_() strings are arrays of C string pointers (the translation function is the identity here);
   - the `#if 0` alternatives are not ported.
   Dependencies not yet ported are declared extern (C) below, in the style of the other ports. */

module gretl.modelprint;


import core.stdc.errno : errno;
import core.stdc.math : exp, fabs, floor, isnan, log, pow, sqrt, ceil, log10, atan2;
import core.stdc.stdio : fprintf, sscanf, sprintf, snprintf, stderr;
import core.stdc.stdlib : calloc, free, malloc, qsort, realloc;
import core.stdc.string : memcpy, memset, strcmp, strcpy, strlen, strchr, strstr, strcat, strncat, strrchr, memmove;

import std.algorithm : max, min;

import gretl.gretl_errorsh : GretlError;
import gretl.gretl_matrix;
import gretl.gretl_matrixh;
import gretl.libgretlh;
import gretl.missingh : NADBL, na;
import gretl.gretl_intl : get_local_decpoint;
import gretl.gstrfuncs : g_unichar_iswide, g_utf8_get_char, g_utf8_next_char, g_utf8_strlen, g_utf8_validate;
import gretl.modelprinth;
import gretl.printouth : GRETL_MP_DIGITS, SHORTLEN;
import gretl.systemh : GretlSystemMethods;
import gretl.compareh : ModelAuxCode;
import gretl.dataio : ntolabel;
import gretl.dataset : is_string_valued, series_get_string_for_value;
import gretl.dataseth : dataset_is_cross_section, dataset_is_panel;
import gretl.gmem : g_free;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_intl : east_asian_locale;
import gretl.gretl_list : gretl_list_separator_position, in_gretl_list;
import gretl.gretl_model : arma_model_nonseasonal_AR_order, arma_model_nonseasonal_MA_order, coeff_pval, gretl_is_between_model, gretl_model_get_coeff_separator, gretl_model_get_data, gretl_model_get_depvar, gretl_model_get_depvar_name, gretl_model_get_double, gretl_model_get_int, gretl_model_get_list, gretl_model_get_param_name, gretl_model_set_double, gretl_model_test_print, model_coeff_pval, model_use_zscore;
import gretl.gretl_modelh : ASYMPTOTIC_MODEL, ArmaFlags, CoeffIntervals, HACKernel, MLVCVType, NONLIST_MODEL, POOLED_MODEL, COUNT_MODEL, PanelVCVType, VCVInfo, VCVMajorType;
import gretl.gretl_paths : gretl_x12_is_x13;
import gretl.gretl_prn : Prn, csv_format, gretl_print_has_minus, gretl_print_set_delim, gretl_print_set_has_minus, gretl_prn_newline, plain_format, pprintf, pputc, pputs, prn_delim, rtf_doc_format, rtf_format, tex_doc_format, tex_format;
import gretl.gstrfuncs : g_strdup, g_strdup_printf;
import gretl.libseth : HACopts;
import gretl.missing : model_has_missing_obs;
import gretl.plugins : get_plugin_function;
import gretl.pvalues : chisq_cdf_comp, normal_cdf, normal_pvalue_2, snedecor_cdf_comp, student_pvalue_2;
import gretl.strutils : gretl_strdup, string_is_blank, strings_array_free, strings_array_new_with_length;
import gretl.printout : bufspace, get_gretl_digits, gretl_fix_exponent, gretl_sprint_fullwidth_double, maybe_trim_varname, outcovmx, print_centered;
import gretl.system : system_short_string;
import gretl.gretl_array : gretl_array_get_strings;

/* unported dependencies: per-file declarations from the 2026c headers */
extern (C) {
    int doing_nls();
    char * get_midas_term_line(const(Model) * pmod, int i);
    void gretl_tex_preamble(Prn * prn, int fmt);
    void make_tex_coeff_name(const(Model) * pmod, const(DatasetRaw) * dset, int i, char * name);
    void tex_coeff_table_end(Prn * prn);
    int tex_coeff_table_start(const(char) ** cols, Gretlopt opt, Prn * prn);
    char * tex_escape(char * targ, const(char) * src);
    char * tex_escape_bounded(char * targ, size_t targ_len, const(char) * src);
    char * tex_escape_special(char * targ, size_t targ_len, const(char) * src);
    char * tex_modify_exponent(char * s);
    void tex_print_coeff(const(ModelCoeff) * mc, Prn * prn);
    char * tex_rl_double(double x, char * s);
    char * tex_rl_float(double x, char * s, int dig);
    char * tex_sprint_double(double x, char * s);
    char * tex_sprint_double_digits(double x, char * s, int dig);
}

enum RTFTAB = "\\par \\ql \\tab ";

int XDIGITS(const(Model) * m)
{
    return ((m.ci == GretlCmdIndex.mpols)? GRETL_MP_DIGITS : get_gretl_digits());
}

int FDIGITS(const(Model) * m)
{
    return ((m.ci == GretlCmdIndex.mpols)? GRETL_MP_DIGITS : 5);
}

bool ordered_model(const(Model) * m)
{
    return ((m.ci == GretlCmdIndex.logit || m.ci == GretlCmdIndex.probit) && gretl_model_get_int(m, "ordered"));
}

bool binary_model(const(Model) * m)
{
    return ((m.ci == GretlCmdIndex.logit || m.ci == GretlCmdIndex.probit) && gretl_model_get_int(m, "binary"));
}

bool multinomial_model(const(Model) * m)
{
    return (m.ci == GretlCmdIndex.logit && gretl_model_get_int(m, "multinom"));
}

bool logit_probit_model(const(Model) * m)
{
    return (m.ci == GretlCmdIndex.logit || m.ci == GretlCmdIndex.probit);
}

bool re_probit_model(const(Model) * m)
{
    return (m.ci == GretlCmdIndex.probit && (m.opt & Gretlopt.OPT_E));
}

bool liml_equation(const(Model) * m)
{
    return (gretl_model_get_int(m, "method") == GretlSystemMethods.liml);
}

bool tsls_model(const(Model) * m)
{
    return (m.ci == GretlCmdIndex.ivreg && !(m.opt & Gretlopt.OPT_L) && !(m.opt & Gretlopt.OPT_G) && !m.aux);
}

bool liml_model(const(Model) * m)
{
    return (m.ci == GretlCmdIndex.ivreg && (m.opt & Gretlopt.OPT_L));
}

bool gmm_model(const(Model) * m)
{
    return (m.ci == GretlCmdIndex.gmm || (m.ci == GretlCmdIndex.ivreg && (m.opt & Gretlopt.OPT_G)));
}

bool intreg_model(const(Model) * m)
{
    return (m.ci == GretlCmdIndex.intreg || m.ci == GretlCmdIndex.tobit);
}

bool hessian_maybe_fishy(const(Model) * m)
{
    return (m.ci == GretlCmdIndex.arma || (m.ci == GretlCmdIndex.probit && (m.opt & Gretlopt.OPT_E)));
}

void model_coeff_init(ModelCoeff * mc)
{
    mc.b = NADBL;
    mc.se = NADBL;
    mc.tval = NADBL;
    mc.pval = NADBL;
    mc.slope = NADBL;
    mc.lo = mc.hi = NADBL;
    mc.show_tval = 1;
    mc.show_pval = 1;
    mc.df_pval = 0;
    mc.multi = 0;
    mc.name[0] = '\0';
}

int glyph_count(const(char) * str)
{
    int ng = 0;

    for (const(char) *p = str; *p; p = g_utf8_next_char(p)) {
        uint u = g_utf8_get_char(p);

        if (g_unichar_iswide(u)) {
            ng += 2;
        } else if (u > 0) {
            ng++;
        }
    }

    return ng;
}

int char_len(const(char) * s)
{
    static int ea = -1;

    if (ea < 0) {
    ea = east_asian_locale();
    }

    if (g_utf8_validate(s, -1, null)) {
    return ea ? glyph_count(s) : cast(int) g_utf8_strlen(s, -1);
    } else {
    return cast(int) strlen(s);
    }
}

void plain_print_double(char * s, int d, double x, Prn * prn)
{
    if (na(x)) {
    strcpy(s, "NA");
    } else if (x < 0 && gretl_print_has_minus(prn)) {
    char[32] tmp;

    *s = '\0';
    strcat(s, "−"); /* U+2212: minus */
    sprintf(tmp.ptr, "%.*g".ptr, d, -x);
    strcat(s, tmp.ptr);
    } else {
    sprintf(s, "%.*g".ptr, d, x);
    }
}

/* for use when printing a user-defined model */

void print_model_stats_table(const(double) * stats, const(char) ** names, int ns, Gretlopt opt, Prn * prn)
{
    int oneline = opt & Gretlopt.OPT_I;
    char[32] tmp1;
    char[32] tmp2;
    char pt = cast(char) get_local_decpoint();
    char delim;
    int i, d;

    if (plain_format(prn)) {
    pputc(prn, '\n');
    } else if (tex_format(prn)) {
    pputs(prn, "\\medskip\n\n");
    pprintf(prn, "\\begin{tabular}{lr@{%c}l}\n".ptr, pt);
    }

    d = get_gretl_digits();
    delim = (pt == ',')? ';' : ',';

    for (i=0; i<ns; i++) {
    if (plain_format(prn)) {
        plain_print_double(tmp1.ptr, d, stats[i], prn);
        if (oneline && ns > 1) {
        if (i == 0) {
            pprintf(prn, "  %s = %s".ptr, names[i], tmp1.ptr);
        } else if (i == ns-1) {
            pprintf(prn, "%c %s = %s\n".ptr, delim, names[i], tmp1.ptr);
        } else {
            pprintf(prn, "%c %s = %s".ptr, delim, names[i], tmp1.ptr);
        }
        } else {
        pprintf(prn, "  %s = %s\n".ptr, names[i], tmp1.ptr);
        }
    } else if (tex_format(prn)) {
        tex_escape_special(tmp1.ptr, tmp1.sizeof, names[i]);
        tex_rl_double(stats[i], tmp2.ptr);
        pprintf(prn, "%s & %s \\\\\n".ptr, tmp1.ptr, tmp2.ptr);
    } else if (rtf_format(prn)) {
        char * ename = rtf_escape_new(names[i]);

        if (na(stats[i])) {
        pprintf(prn, (RTFTAB ~ "%s = NA\n").ptr, ename);
        } else {
        pprintf(prn, (RTFTAB ~ "%s = %g\n").ptr, ename, stats[i]);
        }
        free(ename);
    } else if (csv_format(prn)) {
        char * ename = csv_escape_new(names[i]);

        pprintf(prn, "\"%s\"%c%.15g\n".ptr, ename, prn_delim(prn), stats[i]);
        free(ename);
    }
    }

    if (tex_format(prn)) {
    pputs(prn, "\\end{tabular}");
    }
}

void ensure_vsep(Prn * prn)
{
    if (tex_format(prn)) {
    pputs(prn, "\n\\vspace{1ex}\n");
    } else if (rtf_format(prn)) {
    pputs(prn, "\\par\n");
    }
}

void garch_variance_line(const(Model) * pmod, Prn * prn)
{
    const(char) * varstr = "Unconditional error variance".ptr;
    double v = pmod.sigma * pmod.sigma;
    int d = get_gretl_digits();

    ensure_vsep(prn);

    if (plain_format(prn)) {
    double LR = gretl_model_get_double(pmod, "garch_LR");
    int LRdf = gretl_model_get_int(pmod, "garch_LR_df");

    pprintf(prn, "%s = %.*g\n".ptr, _(varstr), d, v);
    if (pmod.opt & Gretlopt.OPT_Z) {
        pprintf(prn, "%s\n".ptr, _("The residuals are standardized").ptr);
    }
    if (LR >= 0 && LRdf > 0) {
        pprintf(prn, "%s:\n".ptr, _("Likelihood ratio test for (G)ARCH terms").ptr);
        pprintf(prn, "  %s(%d) = %g [%g]\n".ptr,
            _("Chi-square").ptr, LRdf, LR, chisq_cdf_comp(LRdf, LR));
    }
    pputc(prn, '\n');
    } else if (rtf_format(prn)) {
    pprintf(prn, (RTFTAB ~ "%s = %g\n").ptr, _(varstr), v);
    } else if (tex_format(prn)) {
    pprintf(prn, "%s = %g \\\\\n".ptr, _(varstr), v);
    } else if (csv_format(prn)) {
    pprintf(prn, "\"%s\"%c%.15g\n".ptr, _(varstr), prn_delim(prn), v);
    }
}

void print_intreg_info(const(Model) * pmod, const(DatasetRaw) * dset, Prn * prn)
{
    static const(char) *[5] nstrs = ["Left-unbounded observations".ptr, "Right-unbounded observations".ptr, "Bounded observations".ptr, "Point observations".ptr, "Pseudo-point observations".ptr];
    char * lstr = null, rstr = null;
    int nl = gretl_model_get_int(pmod, "n_left");
    int nr = gretl_model_get_int(pmod, "n_right");
    int nb = -1, np = -1, nfp = -1;
    int digits = get_gretl_digits();
    double llim = 0, rlim = NADBL;
    double se_sigma;

    if (pmod.ci == GretlCmdIndex.intreg) {
    nb = gretl_model_get_int(pmod, "n_both");
    np = gretl_model_get_int(pmod, "n_point");
    nfp = gretl_model_get_int(pmod, "n_fpoint");
    } else {
    nstrs[0] = "Left-censored observations".ptr;
    nstrs[1] = "Right-censored observations".ptr;
    if (pmod.opt & Gretlopt.OPT_L) {
        llim = gretl_model_get_double(pmod, "llimit");
        if (tex_format(prn)) {
        lstr = g_strdup_printf(" (%s $\\le$ %g)".ptr,
                       dset.varname[pmod.list[1]],
                       llim);
        } else {
        lstr = g_strdup_printf(" (%s <= %g)".ptr,
                       dset.varname[pmod.list[1]],
                       llim);
        }
    }
    if ((pmod.opt & Gretlopt.OPT_M) && nr > 0) {
        rlim = gretl_model_get_double(pmod, "rlimit");
        if (!na(rlim) && tex_format(prn)) {
        rstr = g_strdup_printf(" (%s $\\ge$ %g)".ptr,
                       dset.varname[pmod.list[1]],
                       rlim);
        } else if (!na(rlim)) {
        rstr = g_strdup_printf(" (%s >= %g)".ptr,
                       dset.varname[pmod.list[1]],
                       rlim);
        }
    }
    }

    ensure_vsep(prn);

    se_sigma = gretl_model_get_double(pmod, "se_sigma");

    if (plain_format(prn)) {
    pprintf(prn, "%s = %.*g".ptr, _("sigma").ptr, digits, pmod.sigma);
    if (!na(se_sigma)) {
        pprintf(prn, " (%g)".ptr, se_sigma);
    }
    pputc(prn, '\n');
    pprintf(prn, "%s: %d%s\n".ptr, _(nstrs[0]), nl, lstr == null ? "".ptr : lstr);
    pprintf(prn, "%s: %d%s\n".ptr, _(nstrs[1]), nr, rstr == null ? "".ptr : rstr);
    if (nb >= 0 && np >= 0) {
        pprintf(prn, "%s: %d\n".ptr, _(nstrs[2]), nb);
        pprintf(prn, "%s: %d\n".ptr, _(nstrs[3]), np);
    }
    if (nfp > 0) {
        pprintf(prn, "%s: %d\n".ptr, _(nstrs[4]), nfp);
    }
    pputc(prn, '\n');
    } else if (rtf_format(prn)) {
    pprintf(prn, (RTFTAB ~ "%s = %g").ptr, _("sigma").ptr, pmod.sigma);
    if (!na(se_sigma)) {
        pprintf(prn, " (%g)".ptr, se_sigma);
    }
    pputc(prn, '\n');
    pprintf(prn, (RTFTAB ~ "%s: %d%s\n").ptr, _(nstrs[0]), nl, lstr == null ? "".ptr : lstr);
    pprintf(prn, (RTFTAB ~ "%s: %d%s\n").ptr, _(nstrs[1]), nr, rstr == null ? "".ptr : rstr);
    if (nb >= 0 && np >= 0) {
        pprintf(prn, (RTFTAB ~ "%s: %d\n").ptr, _(nstrs[2]), nb);
        pprintf(prn, (RTFTAB ~ "%s: %d\n").ptr, _(nstrs[3]), np);
    }
    if (nfp > 0) {
        pprintf(prn, (RTFTAB ~ "%s: %d\n").ptr, _(nstrs[4]), nfp);
    }
    } else if (tex_format(prn)) {
    pprintf(prn, "$\\hat{\\sigma}$ = %g".ptr, pmod.sigma);
    if (!na(se_sigma)) {
        pprintf(prn, " (%g)".ptr, se_sigma);
    }
    pputs(prn, " \\\\\n");
    pprintf(prn, "%s: %d%s \\\\\n".ptr, _(nstrs[0]), nl, lstr == null ? "".ptr : lstr);
    pprintf(prn, "%s: %d%s \\\\\n".ptr, _(nstrs[1]), nr, rstr == null ? "".ptr : rstr);
    if (nb >= 0 && np >= 0) {
        pprintf(prn, "%s: %d \\\\\n".ptr, _(nstrs[2]), nb);
        pprintf(prn, "%s: %d \\\\\n".ptr, _(nstrs[3]), np);
    }
    if (nfp > 0) {
        pprintf(prn, "%s: %d \\\\\n".ptr, _(nstrs[4]), nfp);
    }
    } else if (csv_format(prn)) {
    int d = prn_delim(prn);

    pprintf(prn, "%s%c%.15g\n".ptr, _("sigma").ptr, d, pmod.sigma);
    pprintf(prn, "\"%s\"%c%d\n".ptr, _(nstrs[0]), d, nl);
    pprintf(prn, "\"%s\"%c%d\n".ptr, _(nstrs[1]), d, nr);
    if (nb >= 0 && np >= 0) {
        pprintf(prn, "\"%s\"%c%d\n".ptr, _(nstrs[2]), d, nb);
        pprintf(prn, "\"%s\"%c%d\n".ptr, _(nstrs[3]), d, np);
    }
    if (nfp > 0) {
        pprintf(prn, "\"%s\"%c%d\n".ptr, _(nstrs[4]), d, nfp);
    }
    }

    g_free(lstr);
    g_free(rstr);
}

void rsqline(const(Model) * pmod, Prn * prn)
{
    if (!na(pmod.rsq) && plain_format(prn)) {
    pprintf(prn, "  %s = %f\n".ptr, _("Unadjusted R-squared").ptr, pmod.rsq);
    }
}

void ssrline(const(Model) * pmod, Prn * prn)
{
    if (!na(pmod.ess) && plain_format(prn)) {
    pprintf(prn, "  %s = %.*g\n".ptr, _("Sum of squared residuals").ptr,
        XDIGITS(pmod), pmod.ess);
    }
}

void rssline(const(Model) * pmod, Prn * prn)
{
    if (!na(pmod.ess) && !na(pmod.tss) && plain_format(prn)) {
    pprintf(prn, "  %s = %.*g\n".ptr, _("Explained sum of squares").ptr,
        XDIGITS(pmod), pmod.tss - pmod.ess);
    }
}

void print_liml_equation_data(const(Model) * pmod, Prn * prn)
{
    double lmin = gretl_model_get_double(pmod, "lmin");
    int idf = gretl_model_get_int(pmod, "idf");

    if (!na(lmin)) {
    ensure_vsep(prn);
    if (!pmod.aux) {
        if (tex_format(prn)) {
        pprintf(prn, "%s = %g\\\\\n".ptr, _("Smallest eigenvalue").ptr, lmin);
        } else {
        pprintf(prn, "%s = %g\n".ptr, _("Smallest eigenvalue").ptr, lmin);
        }
    }
    if (idf > 0) {
        double X2 = pmod.nobs * log(lmin);
        double pv = chisq_cdf_comp(idf, X2);

        if (tex_format(prn)) {
        pprintf(prn, "%s: $\\chi^2(%d)$ = %g [%.4f] \\\\\n".ptr,
            _("LR over-identification test").ptr, idf, X2, pv);
        } else {
        pprintf(prn, "%s: ".ptr, _("LR over-identification test").ptr);
        pprintf(prn, "%s(%d) = %g [%.4f]\n\n".ptr, _("Chi-square").ptr,
            idf, X2, pv);
        }
    } else if (idf == 0) {
        pprintf(prn, "%s\n\n".ptr, _("Equation is just identified").ptr);
    }
    }
}

void print_panel_AR_test(double z, int order, Prn * prn)
{
    pprintf(prn, _("Test for AR(%d) errors:").ptr, order);

    if (na(z)) {
    if (tex_format(prn)) {
        pputs(prn, " & $z$ = NA");
    } else {
        pputs(prn, " z = NA");
    }
    } else {
    double pv = normal_pvalue_2(z);

    if (tex_format(prn)) {
        char[32] numstr;

        tex_sprint_double_digits(z, numstr.ptr, 4);
        pprintf(prn, " & $z$ = %s [%.4f]".ptr, numstr.ptr, pv);
    } else {
        pprintf(prn, " z = %g [%.4f]".ptr, z, pv);
    }
    }

    gretl_prn_newline(prn);
}

enum ChiTest {
    dp_sargan,
    dp_hansen,
    dp_wald,
    dp_wald_time,
    j_test,
    overdisp
}

void print_model_chi2_test(const(Model) * pmod, double x, int j, Prn * prn)
{
    static const(char) *[6] strs = ["Sargan over-identification test".ptr, "Hansen over-identification test".ptr, "Wald (joint) test".ptr, "Wald (time dummies)".ptr, "J test".ptr, "Overdispersion test".ptr];
    static const(char) *[6] texstrs = ["Sargan test".ptr, "Hansen test".ptr, "Wald (joint) test".ptr, "Wald (time dummies)".ptr, "J test".ptr, "Overdispersion test".ptr];
    double pv;
    int df;

    if (j == ChiTest.dp_sargan) {
    df = gretl_model_get_int(pmod, "sargan_df");
    } else if (j == ChiTest.dp_hansen) {
    df = gretl_model_get_int(pmod, "hansen_df");
    } else if (j == ChiTest.dp_wald) {
    df = gretl_model_get_int(pmod, "wald_df");
    } else if (j == ChiTest.dp_wald_time) {
    df = gretl_model_get_int(pmod, "wald_time_df");
    } else if (j == ChiTest.j_test) {
    df = gretl_model_get_int(pmod, "J_df");
    } else if (j == ChiTest.overdisp) {
    df = 1;
    } else {
    return;
    }

    if (na(x)) {
    if (j == ChiTest.dp_sargan || j == ChiTest.dp_hansen || j == ChiTest.dp_wald) {
        if (df < 0) {
        df = 0;
        }
        if (tex_format(prn)) {
        pprintf(prn, "%s: ".ptr, _(texstrs[j]));
        pprintf(prn, " $\\chi^2(%d)$ = NA".ptr, df);
        } else if (plain_format(prn)) {
        if (gmm_model(pmod)) {
            pputs(prn, "  ");
        }
        pprintf(prn, "%s: ".ptr, _(strs[j]));
        pprintf(prn, "%s(%d) = NA".ptr, _("Chi-square").ptr, df);
        } else {
        pprintf(prn, "%s: ".ptr, _(strs[j]));
        pprintf(prn, "%s(%d) = NA".ptr, _("Chi-square").ptr, df);
        }
        gretl_prn_newline(prn);
    }
    return;
    }

    pv = chisq_cdf_comp(df, x);
    if (na(pv)) {
    return;
    }

    if (tex_format(prn)) {
    pprintf(prn, "%s: ".ptr, _(texstrs[j]));
    pprintf(prn, " $\\chi^2(%d)$ = %g [%.4f]".ptr, df, x, pv);
    } else if (plain_format(prn)) {
    if (gmm_model(pmod)) {
        pputs(prn, "  ");
    }
    pprintf(prn, "%s: ".ptr, _(strs[j]));
    pprintf(prn, "%s(%d) = %g [%.4f]".ptr, _("Chi-square").ptr, df, x, pv);
    } else {
    pprintf(prn, "%s: ".ptr, _(strs[j]));
    pprintf(prn, "%s(%d) = %g [%.4f]".ptr, _("Chi-square").ptr, df, x, pv);
    }

    gretl_prn_newline(prn);
}

int GMM_crit_line(const(Model) * pmod, Prn * prn)
{
    double Q = pmod.ess;
    double TQ = pmod.ess * pmod.nobs;

    if (plain_format(prn)) {
    pprintf(prn, "  %s: Q = %.*g (TQ = %.*g)\n".ptr, _("GMM criterion").ptr,
        XDIGITS(pmod), Q, XDIGITS(pmod), TQ);
    } else if (rtf_format(prn)) {
    pprintf(prn, (RTFTAB ~ "%s: Q = %g (TQ = %g)\n").ptr, _("GMM criterion").ptr,
        Q, TQ);
    } else if (tex_format(prn)) {
    char[32] x1;
    char[32] x2;

    tex_sprint_double(Q, x1.ptr);
    tex_sprint_double(TQ, x2.ptr);
    pprintf(prn, "%s, $Q$ = %s ($TQ$ = %s)\\\\\n".ptr,
        _("GMM criterion").ptr, x1.ptr, x2.ptr);
    } else if (csv_format(prn)) {
    pprintf(prn, "\"%s\"%c%.15g\n".ptr, _("GMM criterion").ptr,
        prn_delim(prn), Q);
    }

    return 0;
}

void print_GMM_stats(const(Model) * pmod, Prn * prn)
{
    double x;

    ensure_vsep(prn);

    GMM_crit_line(pmod, prn);
    x = gretl_model_get_double(pmod, "J_test");
    if (!na(x)) {
    print_model_chi2_test(pmod, x, ChiTest.j_test, prn);
    }

    if (!tex_format(prn)) {
    gretl_prn_newline(prn);
    }
}

void print_DPD_stats(const(Model) * pmod, Prn * prn)
{
    double x, y;
    int k;

    /* 2021-01-28: mask stats not calculated when showing
       the 1-step model in the context of 2-step estimation
       with the --verbose flag
    */
    if (pmod.ID < 0) {
    return;
    }

    ensure_vsep(prn);

    if (tex_format(prn)) {
    pputs(prn, "\\begin{tabular}{ll}\n");
    }

    k = gretl_model_get_int(pmod, "ninst");
    if (k > 0) {
    pprintf(prn, _("Number of instruments = %d").ptr, k);
    gretl_prn_newline(prn);
    }

    x = gretl_model_get_double(pmod, "AR1");
    y = gretl_model_get_double(pmod, "AR2");
    if (!na(x) && !na(y)) {
    print_panel_AR_test(x, 1, prn);
    print_panel_AR_test(y, 2, prn);
    }

    x = gretl_model_get_double(pmod, "sargan");
    if (!na(x)) {
    print_model_chi2_test(pmod, x, ChiTest.dp_sargan, prn);
    }
    x = gretl_model_get_double(pmod, "hansen");
    if (!na(x)) {
    print_model_chi2_test(pmod, x, ChiTest.dp_hansen, prn);
    }
    x = gretl_model_get_double(pmod, "wald");
    if (!na(x)) {
    print_model_chi2_test(pmod, x, ChiTest.dp_wald, prn);
    }
    x = gretl_model_get_double(pmod, "wald_time");
    if (!na(x)) {
    print_model_chi2_test(pmod, x, ChiTest.dp_wald_time, prn);
    }

    if (tex_format(prn)) {
    pputs(prn, "\\end{tabular}\n");
    } else {
    gretl_prn_newline(prn);
    }
}

void print_overdisp_test(const(Model) * pmod, Prn * prn)
{
    double x = gretl_model_get_double(pmod, "overdisp");

    if (!na(x)) {
    ensure_vsep(prn);
    print_model_chi2_test(pmod, x, ChiTest.overdisp, prn);
    if (!tex_format(prn)) {
        gretl_prn_newline(prn);
    }
    }
}

void print_duration_alpha(const(Model) * pmod, Prn * prn)
{
    if ((pmod.opt & Gretlopt.OPT_B) && plain_format(prn)) {
    /* Weibull */
    double a = 1.0 / pmod.coeff[pmod.ncoeff - 1];
    double sa = a * a * pmod.sderr[pmod.ncoeff - 1];

    ensure_vsep(prn);
    pprintf(prn, "1/sigma = %g (%g)\n".ptr, a, sa);
    if (!tex_format(prn)) {
        gretl_prn_newline(prn);
    }
    }
}

void print_probit_rho(const(Model) * pmod, Prn * prn)
{
    ensure_vsep(prn);

    if (re_probit_model(pmod)) {
    if (tex_format(prn)) {
        pprintf(prn, "$\\hat{\\sigma}_u$ = %.5f\n".ptr, pmod.sigma);
    } else {
        pprintf(prn, "sigma_u = %g\n".ptr, pmod.sigma);
    }
    }

    if (tex_format(prn)) {
    double r = pmod.rho;

    if (r < 0) {
        pprintf(prn, "$\\hat{\\rho}$ = $-$%.5f\n".ptr, fabs(r));
    } else {
        pprintf(prn, "$\\hat{\\rho}$ = %.5f\n".ptr, r);
    }
    } else {
    pprintf(prn, "rho = %g\n".ptr, pmod.rho);
    }

    gretl_prn_newline(prn);
}

void print_GNR_info(const(Model) * pmod, Prn * prn)
{
    double R2 = gretl_model_get_double(pmod, "GNR_Rsquared");
    double tmax = gretl_model_get_double(pmod, "GNR_tmax");
    int rd = gretl_model_get_int(pmod, "near-singular");
    const(char) * msg;

    if (na(R2) && na(tmax) && rd == 0) {
    return;
    }

    ensure_vsep(prn);

    if (!na(R2) && !na(tmax)) {
    if (tex_format(prn)) {
        pprintf(prn, "GNR: $R^2$ = %g, max $|t|$ = %g".ptr, R2, tmax);
    } else {
        pprintf(prn, "GNR: R-squared = %g, max |t| = %g".ptr, R2, tmax);
    }
    gretl_prn_newline(prn);

    if (R2 > 1.0e-8 || tmax > 1.0e-4) {
        msg = "Warning: convergence is questionable".ptr;
    } else {
        msg = "Convergence seems to be reasonably complete".ptr;
    }
    pputs(prn, _(msg));
    gretl_prn_newline(prn);
    } else {
    pputs(prn, _("GNR: got incomplete results").ptr);
    gretl_prn_newline(prn);
    }

    if (rd > 0) {
    msg = (rd > 1)? "Warning: Jacobian is rank-deficient".ptr :
        "Warning: Jacobian near to rank-deficiency".ptr;
    pputs(prn, _(msg));
    gretl_prn_newline(prn);
    }

    if (!tex_format(prn)) {
    gretl_prn_newline(prn);
    }
}

void maybe_print_lad_warning(const(Model) * pmod, Prn * prn)
{
    if (gretl_model_get_int(pmod, "nonunique")) {
    pputs(prn, _("Warning: solution is probably not unique").ptr);
    pputc(prn, '\n');
    }
}

void maybe_print_hessian_warning(const(Model) * pmod, Prn * prn)
{
    if (gretl_model_get_int(pmod, "hess-error")) {
    pputs(prn, _("Warning: couldn't compute numerical Hessian").ptr);
    pputc(prn, '\n');
    }
    if (gretl_model_get_int(pmod, "non-pd-hess")) {
    pputs(prn, _("Warning: non-pd Hessian (but still nonsingular)").ptr);
    pputc(prn, '\n');
    }
}

void panel_variance_lines(const(Model) * pmod, Prn * prn)
{
    double s2v = gretl_model_get_double(pmod, "s2v");
    double s2e = gretl_model_get_double(pmod, "s2e");
    double theta = gretl_model_get_double(pmod, "theta");
    double theta_bar = gretl_model_get_double(pmod, "theta_bar");
    double rsq = gretl_model_get_double(pmod, "corr-rsq");

    if (na(s2v) || na(s2e)) {
    return;
    }

    ensure_vsep(prn);

    if (plain_format(prn)) {
    pprintf(prn, "%s = %g\n".ptr, _("'Between' variance").ptr, s2v);
    pprintf(prn, "%s = %g\n".ptr, _("'Within' variance").ptr, s2e);
    if (!na(theta)) {
        pprintf(prn, "%s = %g\n".ptr, _("theta used for quasi-demeaning").ptr, theta);
    } else if (!na(theta_bar)) {
        pprintf(prn, "%s = %g\n".ptr, _("mean theta").ptr, theta_bar);
    }
    if (!na(rsq)) {
        pprintf(prn, "corr(y,yhat)^2 = %g\n".ptr, rsq);
    }
    pputc(prn, '\n');
    } else if (tex_format(prn)) {
    char[32] xstr;

    tex_sprint_double(s2v, xstr.ptr);
    pprintf(prn, "$\\hat{\\sigma}^2_v$ = %s \\\\\n".ptr, xstr.ptr);
    tex_sprint_double(s2e, xstr.ptr);
    pprintf(prn, "$\\hat{\\sigma}^2_{\\varepsilon}$ = %s \\\\\n".ptr, xstr.ptr);
    if (!na(theta)) {
        tex_sprint_double(theta, xstr.ptr);
        pprintf(prn, "$\\theta$ = %s \\\\\n".ptr, xstr.ptr);
    } else if (!na(theta_bar)) {
        tex_sprint_double(theta_bar, xstr.ptr);
        pprintf(prn, "$\\bar{\\theta}$ = %s \\\\\n".ptr, xstr.ptr);
    }
    } else if (rtf_format(prn)) {
    pprintf(prn, (RTFTAB ~ "%s = %g").ptr, _("'Between' variance").ptr, s2v);
    pprintf(prn, (RTFTAB ~ "%s = %g").ptr, _("'Within' variance").ptr, s2e);
    if (!na(theta)) {
        pprintf(prn, (RTFTAB ~ "%s = %g").ptr, _("theta used for quasi-demeaning").ptr, theta);
    } else if (!na(theta_bar)) {
        pprintf(prn, (RTFTAB ~ "%s = %g").ptr, _("mean theta").ptr, theta_bar);
    }
    } else if (csv_format(prn)) {
    char d = prn_delim(prn);

    pprintf(prn, "\"%s\"%c%.15g\n".ptr, _("'Between' variance").ptr, d, s2v);
    pprintf(prn, "\"%s\"%c%.15g\n".ptr, _("'Within' variance").ptr, d, s2e);
    if (!na(theta)) {
        pprintf(prn, "\"%s\"%c%.15g\n".ptr, _("theta used for quasi-demeaning").ptr,
            d, theta);
    } else if (!na(theta_bar)) {
        pprintf(prn, "\"%s\"%c%.15g\n".ptr, _("mean theta").ptr,
            d, theta_bar);
    }
    }
}

double durbins_h(Model * pmod, int * err)
{
    int ldv = gretl_model_get_int(pmod, "ldepvar");
    double se = pmod.sderr[ldv - 2];
    int T = pmod.nobs;
    double h = NADBL;

    if (pmod.ess <= 0.0 || na(se) || na(pmod.rho)) {
    /* the model is flaky */
    *err = GretlError.badstat;
    } else if ((T * se * se) >= 1.0) {
    /* h is undefined */
    *err = GretlError.sqrt;
    } else {
    h = pmod.rho * sqrt(T / (1 - T * se * se));
    gretl_model_set_double(pmod, "durbin_h", h);
    }

    return h;
}

int least_significant_coeff(const(Model) * pmod)
{
    double x, tmin = 3.2;
    int i, imax, k = 0;

    if (pmod.list == null ||
    gretl_list_separator_position(pmod.list) > 0) {
    return 0;
    }

    /* Don't go beyond the last _listed_ regressor: some
       models may have "extra" coefficients that should be
       ignored here.
    */
    imax = pmod.list[0] - 1;

    for (i=pmod.ifc; i<imax; i++) {
    if (pmod.sderr[i] > 0) {
        x = fabs(pmod.coeff[i] / pmod.sderr[i]);
        if (x < tmin) {
        tmin = x;
        k = i;
        }
    }
    }

    if (tmin < 3.0) {
    x = model_coeff_pval(pmod, tmin);
    if (!na(x) && x > .10) {
        return pmod.list[k+2];
    }
    }

    return 0;
}

void pval_max_line(const(Model) * pmod, const(DatasetRaw) * dset, Prn * prn)
{
    int k = pmod.ncoeff - pmod.ifc;

    if (k < 3) return;

    k = least_significant_coeff(pmod);

    if (k > 0 && k < dset.v) {
    if (pmod.ifc) {
        pprintf(prn, _("Excluding the constant, p-value was highest for variable %d (%s)").ptr, k, dset.varname[k]);
    } else {
        pprintf(prn, _("P-value was highest for variable %d (%s)").ptr,
            k, dset.varname[k]);
    }
    pputs(prn, "\n\n");
    }
}

void print_aux_string(const(Model) * pmod, Prn * prn)
{
    int aux = pmod.aux;
    int plain = plain_format(prn);
    int tex = tex_format(prn);
    int csv = csv_format(prn);
    int close = 1;

    if (plain || tex) {
    pputc(prn, '\n');
    } else if (csv) {
    pputc(prn, '"');
    }

    if (aux == ModelAuxCode.sq) {
    pputs(prn, _("Auxiliary regression for non-linearity test (squared terms)").ptr);
    } else if (aux == ModelAuxCode.log) {
    pputs(prn, _("Auxiliary regression for non-linearity test (log terms)").ptr);
    } else if (aux == ModelAuxCode.add) {
    pputs(prn, _("Auxiliary regression for added variables").ptr);
    } else if (aux == ModelAuxCode.white) {
    pputs(prn, _("White's test for heteroskedasticity").ptr);
    if (pmod.opt & Gretlopt.OPT_X) {
        pprintf(prn, " (%s)".ptr, _("squares only").ptr);
    }
    } else if (aux == ModelAuxCode.bp) {
    pputs(prn, _("Breusch-Pagan test for heteroskedasticity").ptr);
    } else if (aux == ModelAuxCode.het_1) {
    pputs(prn, _("Pesaran-Taylor test for heteroskedasticity").ptr);
    } else if (aux == ModelAuxCode.chow) {
    pputs(prn, _("Augmented regression for Chow test").ptr);
    } else if (aux == ModelAuxCode.coint) {
    if (tex_format(prn)) {
        pputs(prn, _("Cointegrating regression -- ").ptr);
    } else {
        pputs(prn, _("Cointegrating regression - ").ptr);
    }
    } else if (aux == ModelAuxCode.adf) {
    if (tex_format(prn)) {
        pputs(prn, _("Augmented Dickey--Fuller regression").ptr);
    } else {
        pputs(prn, _("Augmented Dickey-Fuller regression").ptr);
    }
    } else if (aux == ModelAuxCode.df) {
    if (tex_format(prn)) {
        pputs(prn, _("Dickey--Fuller regression").ptr);
    } else {
        pputs(prn, _("Dickey-Fuller regression").ptr);
    }
    } else if (aux == ModelAuxCode.kpss) {
    pputs(prn, _("KPSS regression").ptr);
    } else if (aux == ModelAuxCode.reset) {
    pputs(prn, _("Auxiliary regression for RESET specification test").ptr);
    } else if (aux == ModelAuxCode.groupwise) {
    pputs(prn, _("Groupwise heteroskedasticity").ptr);
    } else if (aux == ModelAuxCode.comfac) {
    pputs(prn, _("Augmented regression for common factor test").ptr);
    } else {
    close = 0;
    }

    if (close) {
    if (plain || tex) {
        pputc(prn, '\n');
    } else if (csv) {
        pputs(prn, "\"\n");
    } else { /* RTF */
        pputs(prn, "\\par\n");
    }
    }
}

const(char) * simple_estimator_string(int ci, Prn * prn)
{
    if (ci == GretlCmdIndex.ols || ci == GretlCmdIndex.var) return "OLS".ptr;
    else if (ci == GretlCmdIndex.wls)  return "WLS".ptr;
    else if (ci == GretlCmdIndex.arch) return "WLS (ARCH)".ptr;
    else if (ci == GretlCmdIndex.hsk)  return "Heteroskedasticity-corrected".ptr;
    else if (ci == GretlCmdIndex.ar)   return "AR".ptr;
    else if (ci == GretlCmdIndex.lad)  return "LAD".ptr;
    else if (ci == GretlCmdIndex.quantreg) return "Quantile regression".ptr;
    else if (ci == GretlCmdIndex.mpols) return "High-Precision OLS".ptr;
    else if (ci == GretlCmdIndex.probit) return "Probit".ptr;
    else if (ci == GretlCmdIndex.logit)  return "Logit".ptr;
    else if (ci == GretlCmdIndex.tobit)  return "Tobit".ptr;
    else if (ci == GretlCmdIndex.heckit) return "Heckit".ptr;
    else if (ci == GretlCmdIndex.poisson) return "Poisson".ptr;
    else if (ci == GretlCmdIndex.negbin) return "Negative Binomial".ptr;
    else if (ci == GretlCmdIndex.duration) return "Duration".ptr;
    else if (ci == GretlCmdIndex.nls) return "NLS".ptr;
    else if (ci == GretlCmdIndex.mle) return "ML".ptr;
    else if (ci == GretlCmdIndex.gmm) return "GMM".ptr;
    else if (ci == GretlCmdIndex.logistic) return "Logistic".ptr;
    else if (ci == GretlCmdIndex.garch) return "GARCH".ptr;
    else if (ci == GretlCmdIndex.intreg) return "Interval estimates".ptr;
    else if (ci == GretlCmdIndex.dpanel) return "Dynamic panel".ptr;
    else if (ci == GretlCmdIndex.biprobit) return "Bivariate probit".ptr;
    else if (ci == GretlCmdIndex.midasreg) return "MIDAS".ptr;
    else return "";
}

const(char) * estimator_string(const(Model) * pmod, Prn * prn)
{
    if (pmod.ci == GretlCmdIndex.ar1) {
    if (pmod.opt & Gretlopt.OPT_H) {
        if (tex_format(prn)) return "Hildreth--Lu".ptr;
        else return "Hildreth-Lu".ptr;
    } else if (pmod.opt & Gretlopt.OPT_P) {
        if (tex_format(prn)) return "Prais--Winsten".ptr;
        else return "Prais-Winsten".ptr;
    } else {
        if (tex_format(prn)) return "Cochrane--Orcutt".ptr;
        else return "Cochrane-Orcutt".ptr;
    }
    } else if (pmod.ci == GretlCmdIndex.arma) {
    if (gretl_model_get_int(pmod, "armax")) {
        return "ARMAX".ptr;
    } else if (gretl_model_get_int(pmod, "arima_d") ||
           gretl_model_get_int(pmod, "arima_D")) {
        return "ARIMA".ptr;
    } else {
        return "ARMA".ptr;
    }
    } else if (POOLED_MODEL(pmod)) {
    return "Pooled OLS".ptr;
    } else if (pmod.ci == GretlCmdIndex.panel) {
    if (pmod.opt & Gretlopt.OPT_F) {
        return "Fixed-effects".ptr;
    } else if (pmod.opt & Gretlopt.OPT_U) {
        return "Random-effects (GLS)".ptr;
    } else if (pmod.opt & Gretlopt.OPT_H) {
        if (gretl_model_get_int(pmod, "iters")) {
        return "Maximum Likelihood".ptr;
        } else {
        return "WLS".ptr;
        }
    } else {
        return "Between-groups".ptr;
    }
    } else if (pmod.ci == GretlCmdIndex.dpanel) {
    if (gretl_model_get_int(pmod, "step") == 2) {
        return "2-step dynamic panel".ptr;
    } else {
        return "1-step dynamic panel".ptr;
    }
    } else if (gmm_model(pmod)) {
    if (pmod.opt & Gretlopt.OPT_T) {
        return "2-step GMM".ptr;
    } else if (pmod.opt & Gretlopt.OPT_I) {
        return "Iterated GMM".ptr;
    } else if (gretl_model_get_int(pmod, "step") == 2) {
        return "2-step GMM".ptr;
    } else if (gretl_model_get_int(pmod, "step") > 2) {
        return "Iterated GMM".ptr;
    } else {
        return "1-step GMM".ptr;
    }
    } else if (pmod.ci == GretlCmdIndex.logit) {
    if (gretl_model_get_int(pmod, "ordered")) {
        return "Ordered Logit".ptr;
    } else if (gretl_model_get_int(pmod, "multinom")) {
        return "Multinomial Logit".ptr;
    } else {
        return "Logit".ptr;
    }
    } else if (pmod.ci == GretlCmdIndex.probit) {
    if (pmod.opt & Gretlopt.OPT_E) {
        return "Random-effects probit".ptr;
    } else if (gretl_model_get_int(pmod, "ordered")) {
        return "Ordered Probit".ptr;
    } else {
        return "Probit".ptr;
    }
    } else if (pmod.ci == GretlCmdIndex.logistic) {
    if (pmod.opt & Gretlopt.OPT_F) {
        return "Fixed-effects logistic".ptr;
    } else {
        return "Logistic".ptr;
    }
    } else if (pmod.ci == GretlCmdIndex.heckit) {
    if (pmod.opt & Gretlopt.OPT_T) {
        return "Two-step Heckit".ptr;
    } else {
        return "ML Heckit".ptr;
    }
    } else if (pmod.ci == GretlCmdIndex.lad || pmod.ci == GretlCmdIndex.quantreg) {
    if (gretl_model_get_int(pmod, "rq")) {
        return "Quantile estimates".ptr;
    } else {
        return "LAD".ptr;
    }
    } else if (pmod.ci == GretlCmdIndex.ivreg) {
    if (pmod.opt & Gretlopt.OPT_L) {
        return "LIML".ptr;
    } else {
        return "TSLS".ptr;
    }
    } else if (pmod.ci == GretlCmdIndex.negbin) {
    if (pmod.opt & Gretlopt.OPT_M) {
        return "Negative Binomial 1".ptr;
    } else {
        return "Negative Binomial".ptr;
    }
    } else if (pmod.ci == GretlCmdIndex.duration) {
    if (pmod.opt & Gretlopt.OPT_E) {
        return "Duration (exponential)".ptr;
    } else if (pmod.opt & Gretlopt.OPT_L) {
        return "Duration (log-logistic)".ptr;
    } else if (pmod.opt & Gretlopt.OPT_Z) {
        return "Duration (log-normal)".ptr;
    } else {
        return "Duration (Weibull)".ptr;
    }
    } else if (pmod.ci == GretlCmdIndex.ols &&
           gretl_model_get_int(pmod, "restricted")) {
    return "Restricted OLS".ptr;
    } else if (pmod.ci == GretlCmdIndex.midasreg) {
    if (gretl_model_get_int(pmod, "umidas")) {
        return "MIDAS (OLS)".ptr;
    } else {
        return "MIDAS (NLS)".ptr;
    }
    } else {
    return simple_estimator_string(pmod.ci, prn);
    }
}

int any_tests(const(Model) * pmod)
{
    if (pmod.ntests > 0) {
    return 1;
    }

    if (pmod.ci == GretlCmdIndex.ivreg && gretl_model_get_int(pmod, "stage1-dfn")) {
    return 1;
    }

    return 0;
}

void get_stock_yogo_critvals(int n, int K2, Gretlopt opt, GretlMatrixRaw ** p1, GretlMatrixRaw ** p2)
{
    alias LookupFn = extern (C) GretlMatrixRaw * function(int, int, int);
    LookupFn lookup;

    lookup = cast(LookupFn) get_plugin_function("stock_yogo_lookup");

    if (lookup != null) {
    if (opt & Gretlopt.OPT_L) {
        /* LIML test size */
        *p2 = lookup(n, K2, 3);
    } else {
        /* TSLS relative bias, test size */
        *p1 = lookup(n, K2, 1);
        *p2 = lookup(n, K2, 2);
    }
    }
}

void plain_print_sy_vals(GretlMatrixRaw * v, double g, int k, Gretlopt opt, Prn * prn)
{
    int i, gpos = -1;
    double x;

    pputs(prn, "  ");

    if (k == 1) {
    pputs(prn, _("Critical values for TSLS bias relative to OLS:\n").ptr);
    } else if (opt & Gretlopt.OPT_L) {
    /* xgettext:no-c-format */
    pputs(prn, _("Critical values for desired LIML maximal size, when running\n  tests at a nominal 5% significance level:\n").ptr);
    } else {
    /* xgettext:no-c-format */
    pputs(prn, _("Critical values for desired TSLS maximal size, when running\n  tests at a nominal 5% significance level:\n").ptr);
    }

    pprintf(prn, "\n%9s".ptr, (k == 1)? _("bias").ptr : _("size").ptr);
    for (i=0; i<4; i++) {
    pprintf(prn, "%8g%%".ptr, 100 * gretl_matrix_get(v, 0, i));
    }
    pprintf(prn, "\n%9s".ptr, _("value").ptr);
    for (i=0; i<4; i++) {
    x = gretl_matrix_get(v, 1, i);
    if (gpos < 0 && g > x) {
        gpos = i;
    }
    pprintf(prn, "%9.2f".ptr, x);
    }
    pputs(prn, "\n\n  ");

    if (gpos == 0) {
    x = gretl_matrix_get(v, 0, 0);
    if (k == 1) {
        pprintf(prn, _("Relative bias is probably less than %g%%").ptr, 100 * x);
    } else {
        pprintf(prn, _("Maximal size is probably less than %g%%").ptr, 100 * x);
    }
    } else {
    if (gpos < 0) {
        gpos = 3;
    } else {
        gpos--;
    }
    x = gretl_matrix_get(v, 0, gpos);
    if (k == 1) {
        pprintf(prn, _("Relative bias may exceed %g%%").ptr, 100 * x);
    } else {
        pprintf(prn, _("Maximal size may exceed %g%%").ptr, 100 * x);
    }
    }

    pputs(prn, "\n\n");
}

void maybe_print_weak_insts_test(const(Model) * pmod, Prn * prn)
{
    const(char) * head = "Weak instrument test".ptr;
    double F = gretl_model_get_double(pmod, "stage1-F");
    double g = gretl_model_get_double(pmod, "gmin");
    int digits = get_gretl_digits();
    int got_critvals = 0;
    int dfn = 0, dfd = 0;

    if (na(F) && na(g)) {
    return;
    }

    if (plain_format(prn)) {
    pprintf(prn, "%s - \n".ptr, _(head));
    } else if (tex_format(prn)) {
    pprintf(prn, "%s -- \\\\\n".ptr, _(head));
    } else if (rtf_format(prn)) {
    pprintf(prn, "%s - \\par\n".ptr, _(head));
    }

    if (!na(F)) {
    /* got first-stage F-test (single endogenous regressor) */
    dfn = gretl_model_get_int(pmod, "stage1-dfn");
    dfd = gretl_model_get_int(pmod, "stage1-dfd");

    if (plain_format(prn)) {
        pprintf(prn, "  %s (%d, %d) = %.*g\n".ptr, _("First-stage F-statistic").ptr,
            dfn, dfd, digits, F);
    } else if (tex_format(prn)) {
        char[32] x1str;

        tex_sprint_double(F, x1str.ptr);
        pprintf(prn, "\\quad First-stage $F(%d, %d)$ = %s \\\\\n".ptr, dfn, dfd, x1str.ptr);
    } else if (rtf_format(prn)) {
        pprintf(prn, "  %s (%d, %d) = %g\n".ptr, _("First-stage F-statistic").ptr,
            dfn, dfd, F);
    }
    } else {
    /* got minimum eigenvalue test statistic */
    if (plain_format(prn)) {
        pprintf(prn, "  %s = %.*g\n".ptr, _("Cragg-Donald minimum eigenvalue").ptr,
            digits, g);
    } else if (tex_format(prn)) {
        char[32] x1str;

        tex_sprint_double(g, x1str.ptr);
        pprintf(prn, "\\quad  %s = %s \\\\\n".ptr, _("Cragg--Donald minimum eigenvalue").ptr,
            x1str.ptr);
    } else if (rtf_format(prn)) {
        pprintf(prn, "  %s = %g\n".ptr, _("Cragg-Donald minimum eigenvalue").ptr, g);
    }
    }

    if (!na(g)) {
    /* print Stock-Yogo critical values, if available */
    GretlMatrixRaw * bvals = null;
    GretlMatrixRaw * svals = null;
    int n, K2;

    if (na(F)) {
        n = gretl_model_get_int(pmod, "n");
        K2 = gretl_model_get_int(pmod, "K2");
    } else {
        n = 1;
        K2 = dfn;
    }

    get_stock_yogo_critvals(n, K2, pmod.opt, &bvals, &svals);

    if (bvals != null) {
        got_critvals = 1;
        if (plain_format(prn)) {
        plain_print_sy_vals(bvals, g, 1, pmod.opt, prn);
        }
        gretl_matrix_free(bvals);
    }

    if (svals != null) {
        got_critvals = 1;
        if (plain_format(prn)) {
        plain_print_sy_vals(svals, g, 2, pmod.opt, prn);
        }
        gretl_matrix_free(svals);
    }
    }

    if (!na(F) && !got_critvals && plain_format(prn)) {
    pprintf(prn, "  %s\n\n".ptr, _("A value < 10 may indicate weak instruments").ptr);
    }
}

void print_model_tests(const(Model) * pmod, Prn * prn)
{
    int i;

    if (tex_format(prn)) {
    pputs(prn, "\\vspace{1em}\n\\begin{raggedright}\n");
    for (i=0; i<pmod.ntests; i++) {
        if (i > 0) {
        pputs(prn, "\\vspace{1ex}\n");
        }
        gretl_model_test_print(pmod, i, prn);
    }
    if (pmod.ntests > 0) {
        pputs(prn, "\\vspace{1ex}\n");
    }
    maybe_print_weak_insts_test(pmod, prn);
    pputs(prn, "\\end{raggedright}\n");
    } else {
    for (i=0; i<pmod.ntests; i++) {
        gretl_model_test_print(pmod, i, prn);
    }
    maybe_print_weak_insts_test(pmod, prn);
    }
}

int print_ivreg_instruments(const(Model) * pmod, const(DatasetRaw) * dset, Prn * prn)
{
    static const(char) *[2] labels = ["Instrumented".ptr, "Instruments".ptr];
    const(int) * list;
    char[VNAMELEN] vname;
    int tex = tex_format(prn);
    int i, j, imin, jmin, vi;

    jmin = (pmod.aux)? 1 : 0;

    /* below: first print list of endogenous variables, then
       list of instruments */

    for (j=jmin; j<2; j++) {
    int ccount, nv = 0;

    if (j == 0) {
        list = gretl_model_get_list(pmod, "endolist");
        imin = (list == null)? -1 : 1;
    } else {
        /* j = 1: instruments */
        list = pmod.list;
        imin = gretl_list_separator_position(list);
        if (imin > 0) {
        imin++;
        } else {
        imin = -1;
        }
    }

    if (imin < 0) {
        continue;
    }

    ccount = pprintf(prn, "%s: ".ptr, _(labels[j]));

    for (i=imin; i<=list[0]; i++) {
        vi = list[i];
        if (vi >= dset.v) {
        fprintf(stderr, "print_ivreg_instruments: bad varnum %d\n".ptr, vi);
        continue;
        }
        if (tex) {
        tex_escape_bounded(vname.ptr, vname.sizeof, dset.varname[vi]);
        } else {
        strcpy(vname.ptr, dset.varname[vi]);
        }
        ccount += pprintf(prn, "%s ".ptr, vname.ptr);
        if (ccount >= 64) {
        if (tex) {
            pputs(prn, "\\\\\n");
        } else if (rtf_format(prn)) {
            pputs(prn, "\\par\n");
        } else {
            pputs(prn, "\n  ");
        }
        ccount = 0;
        }
        nv++;
    }

    if (nv == 0) {
        pputs(prn, "none??");
    }

    gretl_prn_newline(prn);
    }

    return 0;
}

void dpd_asy_vcv_line(Prn * prn)
{
    if (csv_format(prn)) {
    pprintf(prn, "\"%s\"".ptr, _("Asymptotic standard errors").ptr);
    } else {
    pputs(prn, _("Asymptotic standard errors").ptr);
    }

    pputc(prn, '\n');
}

void panel_vcv_line(const(VCVInfo) * vi, Prn * prn)
{
    if (vi.vmin == PanelVCVType.hac) {
    if (csv_format(prn)) {
        pprintf(prn, "\"%s\"".ptr, _("Standard errors clustered by unit").ptr);
    } else {
        pputs(prn, _("Standard errors clustered by unit").ptr);
    }
    pputc(prn, '\n');
    } else if (vi.vmin == PanelVCVType.bk) {
    if (csv_format(prn)) {
        pprintf(prn, "\"%s\"".ptr, _("Beck-Katz standard errors").ptr);
    } else if (tex_format(prn)) {
        pputs(prn, _("Beck--Katz standard errors").ptr);
    } else {
        pputs(prn, _("Beck-Katz standard errors").ptr);
    }
    pputc(prn, '\n');
    } else if (vi.vmin == PanelVCVType.time) {
    if (csv_format(prn)) {
        pprintf(prn, "\"%s\"".ptr, _("Standard errors clustered by time").ptr);
    } else {
        pputs(prn, _("Standard errors clustered by time").ptr);
    }
    pputc(prn, '\n');
    } else if (vi.vmin == PanelVCVType.dk) {
    if (csv_format(prn)) {
        pputc(prn, '"');
        pprintf(prn, _("Driscoll-Kraay standard errors, bandwidth %d").ptr, vi.order);
        pputc(prn, '"');
    } else if (tex_format(prn)) {
        pprintf(prn, _("Driscoll--Kraay standard errors, bandwidth %d").ptr, vi.order);
    } else {
        pprintf(prn, _("Driscoll-Kraay standard errors, bandwidth %d").ptr, vi.order);
    }
    pputc(prn, '\n');
    } else if (vi.vmin == PanelVCVType.both) {
    if (csv_format(prn)) {
        pprintf(prn, "\"%s\"".ptr, _("Standard errors clustered by unit and time").ptr);
    } else {
        pputs(prn, _("Standard errors clustered by unit and time").ptr);
    }
    pputc(prn, '\n');
    }
}

void cluster_vcv_line(const(Model) * pmod, const(VCVInfo) * vi, const(DatasetRaw) * dset, Prn * prn)
{
    char * cstr;

    if (vi.cv1 != null && vi.cv2 != null) {
    cstr = g_strdup_printf(_("Standard errors clustered by %s and %s").ptr,
                   vi.cv1, vi.cv2);
    } else if (vi.cv1 != null) {
    int n_c = gretl_model_get_int(pmod, "n_clusters");

    cstr = g_strdup_printf(_("Standard errors clustered by %d values of %s").ptr,
                   n_c, vi.cv1);
    } else {
    cstr = g_strdup(_("Clustered standard errors").ptr);
    }

    if (csv_format(prn)) {
    pprintf(prn, "\"%s\"".ptr, cstr);
    } else {
    pputs(prn, cstr);
    }
    pputc(prn, '\n');

    g_free(cstr);
}

void beck_katz_failed_line(Prn * prn)
{
    if (plain_format(prn)) {
    pputs(prn, _("Could not compute Beck-Katz standard errors").ptr);
    pputc(prn, '\n');
    }
}

void nls_robust_failed_line(Prn * prn)
{
    if (plain_format(prn)) {
    pputs(prn, _("Could not compute robust standard errors").ptr);
    pputc(prn, '\n');
    }
}

void hac_vcv_line(const(Model) * pmod, const(VCVInfo) * vi, Prn * prn)
{
    static const(char) *[3] kstrs = ["Bartlett kernel".ptr, "Parzen kernel".ptr, "QS kernel".ptr];
    int missvals = 0;

    if (vi.vmin == HACKernel.qs) {
    pprintf(prn, _("HAC standard errors, bandwidth %.2f").ptr, vi.bw);
    } else {
    pprintf(prn, _("HAC standard errors, bandwidth %d").ptr, vi.order);
    }

    pprintf(prn, ", %s".ptr, _(kstrs[vi.vmin]));
    if (vi.flags) {
    pprintf(prn, ", %s".ptr, _("prewhitened").ptr);
    }
    pputc(prn, '\n');

    missvals = gretl_model_get_int(pmod, "hac_missvals");
    if (missvals) {
    pprintf(prn, _("Observations not contiguous: %s method used").ptr,
        missvals == HACopts.es ? "ES".ptr : "AM".ptr);
    pputc(prn, '\n');
    }
}

void hc_vcv_line(const(VCVInfo) * vi, Prn * prn)
{
    int hcv = vi.vmin;
    int jack = 0;

    if (vi.vmin == 4) {
    jack = 1;
    hcv--;
    }

    pprintf(prn, "%s, %s%sHC%d%s".ptr,
        _("Heteroskedasticity-robust standard errors").ptr,
        (jack)? "".ptr : _("variant").ptr,
        (jack)? "".ptr : " ".ptr,
        hcv, (jack)? "a (jackknife)".ptr : "".ptr);

    if (rtf_format(prn)) {
    pputs(prn, "\\par\n");
    } else {
    pputc(prn, '\n');
    }
}

void ml_vcv_line(const(VCVInfo) * vi, Prn * prn)
{
    int tex = tex_format(prn);
    const(char) * s = null;

    switch (vi.vmin) {
    case MLVCVType.hessian:
    s = "Standard errors based on Hessian".ptr;
    break;
    case MLVCVType.im:
    s = "Standard errors based on Information Matrix".ptr;
    break;
    case MLVCVType.op:
    s = "Standard errors based on Outer Products matrix".ptr;
    break;
    case MLVCVType.qml:
    s = "QML standard errors".ptr;
    break;
    case MLVCVType.hac:
    s = "HAC standard errors".ptr;
    break;
    case MLVCVType.bw:
    if (tex) {
        s = "Bollerslev--Wooldridge standard errors".ptr;
    } else {
        s = "Bollerslev-Wooldridge standard errors".ptr;
    }
    break;
    case MLVCVType.vcvmax:
    s = "Warning: could not compute standard errors".ptr;
    break;
    default:
    break;
    }

    if (s != null) {
    if (csv_format(prn)) {
        pprintf(prn, "\"%s\"\n".ptr, _(s));
    } else {
        pprintf(prn, "%s\n".ptr, _(s));
    }
    }
}

void rq_vcv_line(const(Model) * pmod, Prn * prn)
{
    int robust = gretl_model_get_int(pmod, "rq_nid");
    double a = gretl_model_get_double(pmod, "rq_alpha");
    int free_s = 0;
    char * s;

    if (!na(a)) {
    if (robust) {
        s = g_strdup_printf("With robust %g percent confidence intervals".ptr,
                100 * (1 - a));
    } else {
        s = g_strdup_printf("With %g percent confidence intervals".ptr,
                100 * (1 - a));
    }
    free_s = 1;
    } else if (robust) {
    s = cast(char *) "Robust (sandwich) standard errors".ptr;
    } else {
    s = cast(char *) "Asymptotic standard errors assuming IID errors".ptr;
    }

    if (csv_format(prn)) {
    pprintf(prn, "\"%s\"".ptr, _(s));
    } else {
    pprintf(prn, "%s".ptr, _(s));
    }

    gretl_prn_newline(prn);

    if (free_s) {
    g_free(s);
    }
}

void tex_vecm_depvar_name(char * s, const(char) * vname)
{
    char[14] tmp;
    int gotit = 0;

    if (sscanf(vname, "d_%13s".ptr, tmp.ptr)) {
    char[24] myvar;

    tex_escape_bounded(myvar.ptr, myvar.sizeof, tmp.ptr);
    sprintf(s, "$\\Delta$%s".ptr, myvar.ptr);
    gotit = 1;
    }

    if (!gotit) {
    tex_escape(s, vname);
    }
}

void tex_dpd_depvar_name(char * s, const(char) * vname)
{
    char[32] vnesc;

    tex_escape_bounded(vnesc.ptr, vnesc.sizeof, vname);
    sprintf(s, "$\\Delta$%s".ptr, vnesc.ptr);
}

void print_model_vcv_info(const(Model) * pmod, const(DatasetRaw) * dset, Prn * prn)
{
    VCVInfo *vi = null;

    if (pmod.ci == GretlCmdIndex.quantreg) {
    rq_vcv_line(pmod, prn);
    } else if (gretl_model_get_int(pmod, "panel_bk_failed")) {
    beck_katz_failed_line(prn);
    } else if (pmod.ci == GretlCmdIndex.dpanel && gretl_model_get_int(pmod, "asy")) {
    dpd_asy_vcv_line(prn);
    } else if ((pmod.ci == GretlCmdIndex.nls || pmod.ci == GretlCmdIndex.midasreg) &&
           gretl_model_get_int(pmod, "non-robust")) {
    nls_robust_failed_line(prn);
    } else {
    vi = cast(typeof(vi)) gretl_model_get_data(pmod, "vcv_info");
    }

    if (vi != null) {
    switch (vi.vmaj) {
    case VCVMajorType.hc:
        hc_vcv_line(vi, prn);
        break;
    case VCVMajorType.hac:
        hac_vcv_line(pmod, vi, prn);
        break;
    case VCVMajorType.ml:
        ml_vcv_line(vi, prn);
        break;
    case VCVMajorType.panel:
        panel_vcv_line(vi, prn);
        break;
    case VCVMajorType.cluster:
        cluster_vcv_line(pmod, vi, dset, prn);
        break;
    default:
        break;
    }
    }
}

void print_extra_list(const(char) * tag, const(int) * list, const(DatasetRaw) * dset, Prn * prn)
{
    int i, v, len;

    len = pputs(prn, _(tag));

    for (i=1; i<=list[0]; i++) {
    v = list[i];
    if (v < dset.v) {
        len += pprintf(prn, " %s".ptr, dset.varname[v]);
    } else {
        len += pprintf(prn, " %d".ptr, v);
    }
    if (len > 68 && i < list[0]) {
        pputc(prn, '\n');
        len = 0;
    }
    }

    pputc(prn, '\n');
}

void print_model_zerolist(const(Model) * pmod, const(DatasetRaw) * dset, Prn * prn)
{
    const(int) * zlist = gretl_model_get_list(pmod, "zerolist");
    const(char) * tag = "Omitted because all values were zero:".ptr;

    if (gretl_is_between_model(pmod)) {
    return;
    }

    print_extra_list(tag, zlist, dset, prn);
}

void print_model_droplist(const(Model) * pmod, const(DatasetRaw) * dset, Prn * prn)
{
    const(int) * dlist = gretl_model_get_list(pmod, "droplist");
    const(char) * tag = "Omitted due to exact collinearity:".ptr;

    if (gretl_is_between_model(pmod) && pmod.dataset == null) {
    /* any droplist references will probably be wrong,
       and may be out of bounds */
    return;
    } else if (dlist[0] > 10) {
    /* let's not spew a ton of series names here */
    pputs(prn, _(tag));
    pprintf(prn, _(" %d series").ptr, dlist[0]);
    pputc(prn, '\n');
    return;
    }

    print_extra_list(tag, dlist, dset, prn);
}

void print_ivreg_droplist(const(Model) * pmod, const(DatasetRaw) * dset, Prn * prn)
{
    const(int) * dlist = gretl_model_get_list(pmod, "inst_droplist");
    int i, v;

    pputs(prn, _("Redundant instruments:").ptr);
    for (i=1; i<=dlist[0]; i++) {
    v = dlist[i];
    if (v < dset.v) {
        pprintf(prn, " %s".ptr, dset.varname[v]);
    } else {
        pprintf(prn, " %d".ptr, v);
    }
    }
    pputc(prn, '\n');
}

void print_arma_depvar(const(Model) * pmod, const(DatasetRaw) * dset, Prn * prn)
{
    int tex = tex_format(prn);
    int yno = gretl_model_get_depvar(pmod);
    int d = gretl_model_get_int(pmod, "arima_d");
    int D = gretl_model_get_int(pmod, "arima_D");
    char[64] vname;

    *vname.ptr = 0;

    if (tex) {
    char[32] tmp;

    if (d > 0 || D > 0) {
        strcat(vname.ptr, "$");
    }
    if (d == 1) {
        strcat(vname.ptr, "(1-L)");
    } else if (d == 2) {
        strcat(vname.ptr, "(1-L)^2");
    }
    if (D == 1) {
        strcat(vname.ptr, "(1-L^s)");
    } else if (D == 2) {
        strcat(vname.ptr, "(1-L^s)^2");
    }
    if (d > 0 || D > 0) {
        strcat(vname.ptr, "$");
    }
    tex_escape_bounded(tmp.ptr, tmp.sizeof, dset.varname[yno]);
    strcat(vname.ptr, tmp.ptr);
    } else {
    if (d == 1) {
        strcat(vname.ptr, "(1-L)");
    } else if (d == 2) {
        strcat(vname.ptr, "(1-L)^2");
    }
    if (D == 1) {
        strcat(vname.ptr, "(1-Ls)");
    } else if (D == 2) {
        strcat(vname.ptr, "(1-Ls)^2");
    }
    if (d > 0 || D > 0) {
        strcat(vname.ptr, " ");
    }
    strcat(vname.ptr, dset.varname[yno]);
    }

    pprintf(prn, "%s: %s".ptr, _("Dependent variable").ptr, vname.ptr);
}

void arma_extra_info(const(Model) * pmod, Prn * prn)
{
    int acode = gretl_model_get_int(pmod, "arma_flags");

    if (acode & ArmaFlags.x12a) {
    if (gretl_x12_is_x13()) {
        pputs(prn, _("Estimated using X-13-ARIMA").ptr);
    } else {
        pputs(prn, _("Estimated using X-12-ARIMA").ptr);
    }
    pputs(prn, " (");
    pputs(prn, (acode & ArmaFlags.exact)? _("exact ML").ptr : _("conditional ML").ptr);
    pputs(prn, ")\n");
    } else if (acode & ArmaFlags.ols) {
    if (gretl_model_get_int(pmod, "null-model")) {
        return;
    }
    pputs(prn, _("Estimated using least squares").ptr);
    pputs(prn, " (");
    pputs(prn, _("= MLE").ptr);
    pputs(prn, ")\n");
    } else if (acode & ArmaFlags.exact) {
    int asnum = gretl_model_get_int(pmod, "as_algo");

    if (asnum > 0) {
        pprintf(prn, _("Estimated using AS %d").ptr, asnum);
    } else {
        pputs(prn, _("Estimated using Kalman filter").ptr);
    }
    pputs(prn, " (");
    pputs(prn, _("exact ML").ptr);
    pputs(prn, ")\n");
    } else if (acode & ArmaFlags.ls) {
    pputs(prn, _("Estimated using least squares").ptr);
    pputs(prn, " (");
    pputs(prn, _("conditional ML").ptr);
    pputs(prn, ")\n");
    } else {
    pputs(prn, _("Estimated using BHHH method").ptr);
    pputs(prn, " (");
    pputs(prn, _("conditional ML").ptr);
    pputs(prn, ")\n");
    }
}

void godfrey_test_string(int ci, int order, Prn * prn)
{
    pputc(prn, '\n');

    if (ci == GretlCmdIndex.ivreg) {
    if (order > 1) {
        pprintf(prn, _("Godfrey (1994) test for autocorrelation up to order %d").ptr,
            order);
    } else {
        pputs(prn, _("Godfrey (1994) test for first-order autocorrelation").ptr);
    }
    } else {
    if (order > 1) {
        pprintf(prn, _("Breusch-Godfrey test for autocorrelation up to order %d").ptr,
            order);
    } else {
        pputs(prn, _("Breusch-Godfrey test for first-order autocorrelation").ptr);
    }
    }

    pputc(prn, '\n');
}

/* The selection variable should be the first variable following
   the list separator */

const(char) * heckit_selvar_name(const(Model) * pmod, const(DatasetRaw) * dset)
{
    const(int) * list = pmod.list;
    int pos = gretl_list_separator_position(list);

    if (pos > 0 && pos < list[0] && list[pos+1] < dset.v) {
    return dset.varname[list[pos+1]];
    } else {
    return null;
    }
}

void print_intreg_depvar(const(Model) * pmod, const(DatasetRaw) * dset, Prn * prn)
{
    int lov = gretl_model_get_int(pmod, "lovar");
    int hiv = gretl_model_get_int(pmod, "hivar");

    if (lov < dset.v && hiv < dset.v) {
    pprintf(prn, "%s: %s".ptr, _("Lower limit").ptr, dset.varname[lov]);
    pprintf(prn, ", %s: %s".ptr, _("Upper limit").ptr, dset.varname[hiv]);
    }
}

void maybe_print_T(const(Model) * pmod, const(DatasetRaw) * dset, const(char) * start, Prn * prn)
{
    if (pmod.ci == GretlCmdIndex.heckit) {
    return;
    } else if (!strcmp(start, "1") && !model_has_missing_obs(pmod)) {
    return;
    } else {
    int xsect = dataset_is_cross_section(dset);

    if (model_has_missing_obs(pmod) || !xsect || strcmp(start, "1")) {
        static const(char) *[2] nstrs = ["n".ptr, "T".ptr];
        const(char) * nstr = xsect ? nstrs[0] : nstrs[1];

        if (tex_format(prn)) {
        pprintf(prn, " ($%s$ = %d)".ptr, _(nstr), pmod.nobs);
        } else {
        pprintf(prn, " (%s = %d)".ptr, _(nstr), pmod.nobs);
        }
    }
    }
}

void make_obs_sep(char * targ, const(char) * obs, int tex)
{
    if (tex) {
    strcpy(targ, "--");
    } else if (strchr(obs, '-') != null) {
    strcpy(targ, ":");
    } else {
    strcpy(targ, "-");
    }
}

void maybe_show_midas_method(const(Model) * pmod, Prn * prn)
{
    if (gretl_model_get_int(pmod, "umidas") == 0) {
    gretl_prn_newline(prn);
    if (gretl_model_get_int(pmod, "GSS")) {
        pputs(prn, _("Using line search with conditional OLS").ptr);
    } else if (gretl_model_get_int(pmod, "BFGS")) {
        pputs(prn, _("Using L-BFGS-B with conditional OLS").ptr);
    } else if (tex_format(prn)) {
        pputs(prn, _("Using Levenberg--Marquardt algorithm").ptr);
    } else {
        pputs(prn, _("Using Levenberg-Marquardt algorithm").ptr);
    }
    }
}

void maybe_show_SA_method(const(Model) * pmod, Prn * prn)
{
    const(char) * s = cast(const(char) *) gretl_model_get_data(pmod, "anova_method");

    if (s != null) {
    gretl_prn_newline(prn);
    if (!strcmp(s, "stata")) {
        pputs(prn, _("Using Swamy-Arora as per Stata").ptr);
    } else {
        pputs(prn, _("Using Swamy-Arora as per Baltagi-Chang").ptr);
    }
    }
}

void print_model_heading(const(Model) * pmod, const(DatasetRaw) * dset, Gretlopt opt, Prn * prn)
{
    char[OBSLEN] startdate;
    char[OBSLEN] enddate;
    char[48] vname;
    char[4] datesep;
    int t1 = pmod.t1, t2 = pmod.t2;
    int tex = tex_format(prn);
    int csv = csv_format(prn);
    int dvnl = 1;
    int order = 0;

    /* flag indicating whether we've printed an opening quote,
       in the case of CSV output */
    int csv_iniquote = 0;

    if (pmod.aux != ModelAuxCode.var && pmod.aux != ModelAuxCode.vecm) {
    ntolabel(startdate.ptr, t1, dset);
    ntolabel(enddate.ptr, t2, dset);
    make_obs_sep(datesep.ptr, startdate.ptr, tex);
    }

    switch (pmod.aux) {
    case ModelAuxCode.sq:
    case ModelAuxCode.log:
    case ModelAuxCode.add:
    case ModelAuxCode.white:
    case ModelAuxCode.bp:
    case ModelAuxCode.het_1:
    case ModelAuxCode.chow:
    case ModelAuxCode.coint:
    case ModelAuxCode.adf:
    case ModelAuxCode.df:
    case ModelAuxCode.kpss:
    case ModelAuxCode.reset:
    case ModelAuxCode.groupwise:
    case ModelAuxCode.comfac:
    case ModelAuxCode.biprob:
    print_aux_string(pmod, prn);
    break;
    case ModelAuxCode.ar:
    order = gretl_model_get_int(pmod, "BG_order");
    if (order > 0) {
        godfrey_test_string(pmod.ci, order, prn);
    } else {
        pputc(prn, '\n');
    }
    break;
    case ModelAuxCode.arch:
    order = gretl_model_get_int(pmod, "arch_order");
    pputc(prn, '\n');
    pprintf(prn, _("Test for ARCH of order %d").ptr, order);
    pputc(prn, '\n');
    break;
    case ModelAuxCode.sys:
    pprintf(prn, "%s %d: ".ptr, _("Equation").ptr, pmod.ID + 1);
    break;
    case ModelAuxCode.var:
    pprintf(prn, "\n%s %d: ".ptr, _("Equation").ptr, pmod.ID);
    break;
    case ModelAuxCode.vecm:
    pprintf(prn, "%s %d: ".ptr, _("Equation").ptr, pmod.ID);
    break;
    case ModelAuxCode.aux:
    pputc(prn, '\n');
    break;
    default:
    /* 2017-03-12: first test was 'pmod->ID < 0' */
    if (pmod.ID <= 0 || (opt & Gretlopt.OPT_S)) {
        if (!csv) {
        pputc(prn, '\n');
        }
    } else if (pmod.name) {
        char[32] modname;

        if (tex) {
        tex_escape_bounded(modname.ptr, modname.sizeof, pmod.name);
        } else {
        strcpy(modname.ptr, pmod.name);
        }

        if (csv) {
        pprintf(prn, "\"%s:\"\n".ptr, modname.ptr);
        csv_iniquote = 1;
        } else if (strlen(pmod.name) > 8) {
        pprintf(prn, "\n%s:\n".ptr, modname.ptr);
        } else {
        pprintf(prn, "\n%s: ".ptr, modname.ptr);
        }
    } else {
        if (csv) {
        pprintf(prn, "\"%s %d: ".ptr, _("Model").ptr, pmod.ID);
        csv_iniquote = 1;
        } else {
        pprintf(prn, "\n%s %d: ".ptr, _("Model").ptr, pmod.ID);
        }
    }
    break;
    }

    if (csv && !csv_iniquote) {
    pputc(prn, '"');
    }

    if (pmod.aux == ModelAuxCode.var || pmod.aux == ModelAuxCode.vecm) {
    ;
    } else if (pmod.aux == ModelAuxCode.sys) {
    pprintf(prn, _("%s, using observations %s%s%s").ptr,
        _(system_short_string(pmod)),
        startdate.ptr, datesep.ptr, enddate.ptr);
    maybe_print_T(pmod, dset, startdate.ptr, prn);
    } else if (!dataset_is_panel(dset)) {
    int mc, Tmax = pmod.t2 - pmod.t1 + 1;
    const(char) * estr = estimator_string(pmod, prn);
    const(char) * fmt;

    if (dset.v == 0 || gretl_model_get_int(pmod, "matrix_data")) {
        fmt = "%s, using %d observations".ptr;
        pprintf(prn, _(fmt), _(estr), pmod.nobs);
    } else if (char_len(estr) > 32) {
        fmt = "%s, obs. %s%s%s".ptr;
        pprintf(prn, _(fmt), _(estr), startdate.ptr,
            (tex)? "--".ptr : "-".ptr, enddate.ptr);
        maybe_print_T(pmod, dset, startdate.ptr, prn);
    } else {
        fmt = "%s, using observations %s%s%s".ptr;
        pprintf(prn, _(fmt), _(estr), startdate.ptr,
            datesep.ptr, enddate.ptr);
        maybe_print_T(pmod, dset, startdate.ptr, prn);
    }

    if (pmod.ci == GretlCmdIndex.midasreg) {
        maybe_show_midas_method(pmod, prn);
    }

    if (gretl_model_get_int(pmod, "binary_obs_dropped")) {
        /* the relevant info should be printed already */
        mc = 0;
    } else if (pmod.ci == GretlCmdIndex.heckit) {
        mc = Tmax - gretl_model_get_int(pmod, "totobs");
    } else {
        mc = Tmax - pmod.nobs;
        if (pmod.ci == GretlCmdIndex.wls) {
        mc -= gretl_model_get_int(pmod, "wt_zeros");
        }
    }

        if (gretl_model_get_int(pmod, "JWJ_fail")) {
            gretl_prn_newline(prn);
            pprintf(prn, _("J'WJ is not positive definite, std. errors not available").ptr);
        }
    if (mc > 0) {
        gretl_prn_newline(prn);
        pprintf(prn, "%s: %d".ptr, _("Missing or incomplete observations dropped").ptr,
            mc);
    }
    } else {
    /* panel data */
    int effn = gretl_model_get_int(pmod, "n_included_units");
    int Tmin = gretl_model_get_int(pmod, "Tmin");
    int Tmax = gretl_model_get_int(pmod, "Tmax");

    pprintf(prn, _("%s, using %d observations").ptr,
        _(estimator_string(pmod, prn)), pmod.nobs);

    if (pmod.opt & Gretlopt.OPT_U) {
        /* random effects */
        if (pmod.opt & Gretlopt.OPT_E) {
        gretl_prn_newline(prn);
        if (pmod.opt & Gretlopt.OPT_X) {
            pputs(prn, _("Using weighted Nerlove transformation").ptr);
        } else {
            pputs(prn, _("Using Nerlove's transformation").ptr);
        }
        } else if (pmod.opt & Gretlopt.OPT_X) {
        /* Swamy-Arora special */
        maybe_show_SA_method(pmod, prn);
        }
    }

    if (effn > 0) {
        gretl_prn_newline(prn);
        pprintf(prn, _("Included %d cross-sectional units").ptr, effn);
    }
    if (Tmin > 0 && Tmax > 0) {
        gretl_prn_newline(prn);
        if (Tmin == Tmax) {
        pprintf(prn, _("Time-series length = %d").ptr, Tmin);
        } else {
        pprintf(prn, _("Time-series length: minimum %d, maximum %d").ptr,
            Tmin, Tmax);
        }
    }
    if (pmod.ci == GretlCmdIndex.dpanel) {
        if (pmod.opt & Gretlopt.OPT_L) {
        gretl_prn_newline(prn);
        pputs(prn, _("Including equations in levels").ptr);
        }
        if (pmod.opt & Gretlopt.OPT_X) {
        gretl_prn_newline(prn);
        pputs(prn, _("H-matrix as per Ox/DPD").ptr);
        }
    }
    }

    if (csv) pputc(prn, '"');

    if (pmod.aux != ModelAuxCode.var && pmod.aux != ModelAuxCode.vecm) {
    gretl_prn_newline(prn);
    }

    if (pmod.ci == GretlCmdIndex.arma && plain_format(prn)) {
    arma_extra_info(pmod, prn);
    }

    if (csv) pputc(prn, '"');

    /* special formulations for dependent variable in various cases */
    if (pmod.aux == ModelAuxCode.sq || pmod.aux == ModelAuxCode.log || pmod.aux == ModelAuxCode.add) {
    pprintf(prn, "%s: %s".ptr, _("Dependent variable").ptr,
        (tex)? "$\\hat{u}$".ptr : "uhat".ptr);
    } else if (pmod.aux == ModelAuxCode.white || pmod.aux == ModelAuxCode.het_1) {
    pprintf(prn, "%s: %s".ptr, _("Dependent variable").ptr,
        (tex)? "$\\hat{u}^2$".ptr : "uhat^2".ptr);
    } else if (pmod.aux == ModelAuxCode.bp) {
    const(char) * fmt;

    if (pmod.opt & Gretlopt.OPT_R) {
        fmt = "scaled %s (Koenker robust variant)".ptr;
    } else {
        fmt = "scaled %s".ptr;
    }
    pprintf(prn, "%s: ".ptr, _("Dependent variable").ptr);
    pprintf(prn, _(fmt), (tex)? "$\\hat{u}^2$".ptr : "uhat^2".ptr);
    } else if (pmod.aux == ModelAuxCode.arch) {
    pprintf(prn, "%s: %s".ptr, _("Dependent variable").ptr,
        (tex)? "$u_t^2$".ptr : "ut^2".ptr);
    } else if (pmod.ci == GretlCmdIndex.nls || pmod.ci == GretlCmdIndex.mle || pmod.ci == GretlCmdIndex.gmm) {
    if (pmod.depvar != null) {
        if (tex) {
        pprintf(prn, "\\verb!%s!".ptr, pmod.depvar);
        } else {
        pputs(prn, pmod.depvar);
        }
    } else {
        dvnl = 0;
    }
    } else if (pmod.ci == GretlCmdIndex.arma) {
    print_arma_depvar(pmod, dset, prn);
    } else if (pmod.ci == GretlCmdIndex.intreg) {
    print_intreg_depvar(pmod, dset, prn);
    } else if (pmod.ci == GretlCmdIndex.biprobit) {
    dvnl = 0;
    } else {
    const(char) * dvname = gretl_model_get_depvar_name(pmod, dset);

    if (tex) {
        if (pmod.aux == ModelAuxCode.vecm) {
        tex_vecm_depvar_name(vname.ptr, dvname);
        } else if (pmod.ci == GretlCmdIndex.dpanel) {
        tex_dpd_depvar_name(vname.ptr, dvname);
        } else {
        tex_escape_bounded(vname.ptr, vname.sizeof, dvname);
        }
    }

    if (pmod.aux == ModelAuxCode.var || pmod.aux == ModelAuxCode.vecm) {
        pputs(prn, (tex)? vname.ptr : dvname);
    } else {
        pprintf(prn, "%s: %s".ptr, _("Dependent variable").ptr,
            (tex)? vname.ptr : dvname);
    }
    }

    if (csv) pputc(prn, '"');

    if (dvnl) {
    gretl_prn_newline(prn);
    }

    /* supplementary strings below the estimator and sample info */

    if (pmod.ci == GretlCmdIndex.ivreg) {
    /* list of instruments for IV estimation */
    int method = gretl_model_get_int(pmod, "method");

    if (method != GretlSystemMethods.fiml && method != GretlSystemMethods.liml) {
        print_ivreg_instruments(pmod, dset, prn);
    }
    } else if (pmod.ci == GretlCmdIndex.lad || pmod.ci == GretlCmdIndex.quantreg) {
    /* tau for quantile regression */
    double tau = gretl_model_get_double(pmod, "tau");

    if (!na(tau)) {
        if (tex) {
        pprintf(prn, "$\\tau$ = %g".ptr, tau);
        } else {
        pprintf(prn, "tau = %g".ptr, tau);
        }
        gretl_prn_newline(prn);
    }
    } else if (pmod.ci == GretlCmdIndex.heckit) {
    /* selection variable for Heckit */
    const(char) * selvar = heckit_selvar_name(pmod, dset);

    if (selvar != null) {
        if (csv) pputc(prn, '"');
        if (tex) {
        tex_escape_bounded(vname.ptr, vname.sizeof, selvar);
        }
        pprintf(prn, "%s: %s".ptr, _("Selection variable").ptr,
            (tex)? vname.ptr : selvar);
        if (csv) pputc(prn, '"');
        pputc(prn, '\n');
    }
    } else if (pmod.ci == GretlCmdIndex.ar1) {
    double rho = gretl_model_get_double(pmod, "rho_gls");

    if (!na(rho)) {
        if (tex) {
        pprintf(prn, "$\\rho$ = %g".ptr, rho);
        } else {
        pprintf(prn, "rho = %g".ptr, rho);
        }
        gretl_prn_newline(prn);
    }
    } else if (pmod.ci == GretlCmdIndex.probit && (pmod.opt & Gretlopt.OPT_E)) {
    int qp = gretl_model_get_int(pmod, "quadpoints");

    if (qp > 0) {
        pprintf(prn, _("Using %d quadrature points").ptr, qp);
        gretl_prn_newline(prn);
    }
    } else if (pmod.ci == GretlCmdIndex.hsk && (pmod.opt & Gretlopt.OPT_N)) {
    pputs(prn, _("Without squared terms in variance equation").ptr);
    gretl_prn_newline(prn);
    }

    /* VCV variants */
    print_model_vcv_info(pmod, dset, prn);

    if (pmod.ci == GretlCmdIndex.panel && (pmod.opt & Gretlopt.OPT_H) && !pmod.aux) {
    /* WLS on panel data */
    if (tex) {
        pputs(prn, "\\\\\n");
    }
    if (gretl_model_get_int(pmod, "iters")) {
        pprintf(prn, _("Allowing for groupwise heteroskedasticity").ptr);
    } else {
        pprintf(prn, _("Weights based on per-unit error variances").ptr);
    }
    pputc(prn, '\n');
    } else if (pmod.ci == GretlCmdIndex.wls && !pmod.aux) {
    /* weight variable for WLS */
    if (tex) {
        tex_escape_bounded(vname.ptr, vname.sizeof, dset.varname[pmod.nwt]);
    }
    if (csv) pputc(prn, '"');
    pprintf(prn, "%s: %s".ptr, _("Variable used as weight").ptr,
        (tex)? vname.ptr : dset.varname[pmod.nwt]);
    if (csv) pputc(prn, '"');
    pputc(prn, '\n');
    } else if (pmod.ci == GretlCmdIndex.arch) {
    /* weight variable for ARCH */
    if (csv) pputc(prn, '"');
    pprintf(prn, "%s: %s".ptr, _("Variable used as weight").ptr,
        (tex)? "$1/\\hat{\\sigma}_t$".ptr : "1/sigma".ptr);
    if (csv) pputc(prn, '"');
    pputc(prn, '\n');
    } else if (pmod.ci == GretlCmdIndex.ar1) {
    /* rhohat for AR1 (TeX only) */
    if (tex) {
        pprintf(prn, "$\\hat{\\rho}$ = %g\n".ptr,
            gretl_model_get_double(pmod, "rho_gls".ptr));
    }
    } else if (pmod.ci == GretlCmdIndex.logistic) {
    /* y-hat formula for logistic regression */
    if (tex) {
        pprintf(prn, "$\\hat{y} \\approx E\\left(%g / (1 + e^{-X\\hat{\\beta}})\\right)$\n".ptr,
            gretl_model_get_double(pmod, "lmax".ptr));
    } else {
        pprintf(prn, "yhat =~ E(%g / (1 + exp(-X*b)))\n".ptr,
            gretl_model_get_double(pmod, "lmax".ptr));
    }
    }

    /* IVREG: message about redundant instruments */
    if (plain_format(prn) && pmod.ci == GretlCmdIndex.ivreg &&
    gretl_model_get_list(pmod, "inst_droplist") != null) {
    print_ivreg_droplist(pmod, dset, prn);
    }

    /* messages about collinear and/or zero regressors */
    if (plain_format(prn)) {
    if (gretl_model_get_list(pmod, "zerolist") != null) {
        print_model_zerolist(pmod, dset, prn);
    }
    if (gretl_model_get_list(pmod, "droplist") != null) {
        print_model_droplist(pmod, dset, prn);
    }
    }

    if (plain_format(prn) && (pmod.ci == GretlCmdIndex.lad || pmod.ci == GretlCmdIndex.quantreg)) {
    maybe_print_lad_warning(pmod, prn);
    }

    if (plain_format(prn) && hessian_maybe_fishy(pmod)) {
    maybe_print_hessian_warning(pmod, prn);
    }

    if (pmod.missmask == null && gretl_model_get_int(pmod, "wt_dummy")) {
    /* FIXME alt formats */
    pprintf(prn, "%s %d\n".ptr,
        _("Weight var is a dummy variable, effective obs =").ptr,
        pmod.nobs);
    }

    if (rtf_format(prn)) {
    pputs(prn, "\\par\n");
    } else {
    pputc(prn, '\n');
    }
}

void model_format_start(Prn * prn)
{
    if (tex_format(prn)) {
    if (tex_doc_format(prn)) {
        gretl_tex_preamble(prn, 0);
    }
    pputs(prn, "\\begin{center}\n");
    } else if (rtf_format(prn)) {
    gretl_print_set_has_minus(prn);
    if (rtf_doc_format(prn)) {
        pputs(prn, "{\\rtf1\\par\n\\qc ");
    } else {
        pputs(prn, "\\par\n\\qc ");
    }
    }
}

enum RTF_MULTICOL = "\\trowd \\trqc \\trgaph30\\trleft-30\\trrh262\\cellx8000\n\\intbl";

enum RTF_COEFF_ROW = "\\trowd \\trqc \\trgaph30\\trleft-30\\trrh262\\cellx1900\\cellx3300\\cellx4700\\cellx6100\\cellx7500\\cellx8000\n\\intbl";

enum RTF_BINARY_ROW = "\\trowd \\trqc \\trgaph30\\trleft-30\\trrh262\\cellx1900\\cellx3300\\cellx4700\\cellx6100\\cellx8000\n\\intbl";

enum RTF_INTVL_ROW = "\\trowd \\trqc \\trgaph30\\trleft-30\\trrh262\\cellx1900\\cellx3300\\cellx4700\\cellx6100\n\\intbl";

enum RTF_ROOT_ROW = "\\trowd \\trqc \\trgaph30\\trleft-30\\trrh262\\cellx500\\cellx1500\\cellx2900\\cellx4300\\cellx5700\\cellx7100\n\\intbl";

/* below: this is used when we're doing something other than a plain
   text print of a model. When using TeX, returns the number of
   columns in the table, otherwise just returns 0.
*/

int alt_print_coeff_table_start(const(Model) * pmod, int ci, Prn * prn)
{
    const(char) * tlabel;
    int use_param = 0;
    int slopes = 0;
    int intervals = 0;
    int seqcols = 0;
    int mp = 0, ret = 0;

    if (model_use_zscore(pmod)) {
    tlabel = (tex_format(prn))? "$z$".ptr : "z".ptr;
    } else {
    tlabel = (tex_format(prn))? "$t$-ratio".ptr : "t-ratio".ptr;
    }

    if (pmod != null) {
    GretlMatrixRaw * m;

    use_param = NONLIST_MODEL(cast(GretlCmdIndex) pmod.ci);
    slopes = binary_model(pmod) && !(pmod.opt & Gretlopt.OPT_P);
    intervals = (gretl_model_get_data(pmod, "coeff_intervals") != null);
    m = cast(typeof(m)) gretl_model_get_data(pmod, "rq_sequence");
    seqcols = gretl_matrix_cols(m);
    mp = (pmod.ci == GretlCmdIndex.mpols);
    }

    if (csv_format(prn)) {
    char d = prn_delim(prn);

    if (mp) {
        pprintf(prn, "%c\"%s\"%c\"%s\"\n".ptr,
            d, _("coefficient").ptr, d, _("std. error").ptr);
    } else if (slopes) {
        pprintf(prn, "%c\"%s\"%c\"%s\"%c\"%s\"%c\"%s\"\n".ptr,
            d, _("coefficient").ptr, d, _("std. error").ptr,
            d, _(tlabel), d, _("slope at mean").ptr);
    } else if (use_param) {
        pprintf(prn, "%c\"%s\"%c\"%s\"%c\"%s\"%c\"%s\"\n".ptr,
            d, _("estimate").ptr, d, _("std. error").ptr,
            d, _(tlabel), d, _("p-value").ptr);
    } else if (intervals) {
        pprintf(prn, "%c\"%s\"%c\"%s\"%c\"%s\"\n".ptr,
            d, _("coefficient").ptr, d, _("lower").ptr,
            d, _("upper").ptr);
    } else if (seqcols == 3) {
        pprintf(prn, "%c\"%s\"%c\"%s\"%c\"%s\"%c\"%s\"\n".ptr,
            d, "tau".ptr, d, _("coefficient").ptr, d, _("lower").ptr,
            d, _("upper").ptr);
    } else if (seqcols == 2) {
        pprintf(prn, "%c\"%s\"%c\"%s\"%c\"%s\"\n".ptr,
            d, _("coefficient").ptr, d, _("std. error").ptr,
            d, _(tlabel));
    } else {
        pprintf(prn, "%c\"%s\"%c\"%s\"%c\"%s\"%c\"%s\"\n".ptr,
            d, _("coefficient").ptr, d, _("std. error").ptr,
            d, _(tlabel), d, _("p-value").ptr);
    }
    } else {
    const(char) *[6] cols = null;
    int i;

    cols[0] = " ";
    cols[1] = (use_param)? "Estimate".ptr : "Coefficient".ptr;

    if (intervals) {
        cols[2] = "Lower".ptr;
        cols[3] = "Upper".ptr;
    } else if (seqcols == 3) {
        cols[2] = "tau".ptr;
        cols[3] = "Lower".ptr;
        cols[4] = "Upper".ptr;
    } else {
        cols[2] = (tex_format(prn))? "Std.\\ Error".ptr : "Std. Error".ptr;
        if (!mp) {
        cols[3] = tlabel;
        cols[4] = (slopes)? "Slope".ptr : "p-value".ptr;
        }
    }

    if (tex_format(prn)) {
        Gretlopt tabopt = (ci == GretlCmdIndex.modprint)? Gretlopt.OPT_U : Gretlopt.OPT_NONE;

        if (slopes) {
        tabopt |= Gretlopt.OPT_B; /* "binary" */
        }
        if (mp) {
        tabopt |= Gretlopt.OPT_M; /* multiple precision */
        }
        ret = tex_coeff_table_start(cols.ptr, tabopt, prn);
    } else if (rtf_format(prn)) {
        pputc(prn, '{');
        if (intervals) {
        pputs(prn, RTF_INTVL_ROW.ptr);
        } else if (slopes) {
        pputs(prn, RTF_BINARY_ROW.ptr);
        } else {
        pputs(prn, RTF_COEFF_ROW.ptr);
        }
        for (i=0; cols[i] != null; i++) {
        if (slopes && i == 4) {
            pprintf(prn, " \\qc {\\i %s{\\super *}}\\cell".ptr, _(cols[i]));
        } else {
            pprintf(prn, " \\qc {\\i %s}\\cell".ptr, _(cols[i]));
        }
        }
        if (!slopes && !intervals) {
        pputs(prn, " \\ql \\cell");
        }
        pputs(prn, " \\intbl \\row\n");
    }
    }

    return ret;
}

void print_coeff_table_end(const(Model) * pmod, Prn * prn)
{
    if (plain_format(prn) || csv_format(prn)) {
    pputc(prn, '\n');
    } else if (tex_format(prn)) {
    tex_coeff_table_end(prn);
    } else if (rtf_format(prn)) {
    pputs(prn, "}\n\n");
    }

    if (pmod == null) {
    /* no MODEL is given when we're called via the
       "modprint" command */
    return;
    }

    /* Maybe print "near-singularity" message? Not in the
       case of NLS or MIDAS regression, when this will be
       handled in conjunction with printing GNR info.
    */
    if (pmod.ci == GretlCmdIndex.nls || pmod.ci == GretlCmdIndex.midasreg) {
    ;
    } else if (plain_format(prn) && gretl_model_get_int(pmod, "near-singular")) {
    const(char) * msg = "Warning: data matrix close to singularity!".ptr;

    pprintf(prn, "%s\n\n".ptr, _(msg));
    }
}

void model_format_end(Prn * prn)
{
    if (tex_format(prn)) {
    pputs(prn, "\n\\end{center}\n");
    if (tex_doc_format(prn)) {
        pputs(prn, "\n\\end{document}\n");
    }
    } else if (rtf_doc_format(prn)) {
    pputs(prn, "\n}\n");
    }
}

void addconst_message(const(Model) * pmod, Prn * prn)
{
    if (gretl_model_get_int(pmod, "addconst")) {
    pprintf(prn, "\n\n%s\n".ptr,
        _("WARNING: The constant was present among the regressors but not among the\ninstruments, so it has been automatically added to the instrument list.\nThis behavior may change in future versions, so you may want to adjust your\nscripts accordingly.\n").ptr);
    }
}

enum AltStats {
    orig_stats,
    wtd_stats,
    rhodiff_stats,
    transform_stats
}

void alternate_stats_message(int i, Prn * prn)
{
    static const(char) *[4] msg = ["Statistics based on the original data".ptr, "Statistics based on the weighted data".ptr, "Statistics based on the rho-differenced data".ptr, "Statistics based on the transformed data".ptr];

    if (plain_format(prn)) {
    pprintf(prn, "%s:\n\n".ptr, _(msg[i]));
    } else if (tex_format(prn)) {
    pprintf(prn, "\\vspace{1em}%s:\n\n".ptr, _(msg[i]));
    } else if (csv_format(prn)) {
    pprintf(prn, "\"%s\"\n".ptr, _(msg[i]));
    } else {
    /* RTF */
    pprintf(prn, "\\par \\qc\n%s:\n\n".ptr, _(msg[i]));
    }
}

void print_whites_results(const(Model) * pmod, Prn * prn)
{
    double X = pmod.rsq * pmod.nobs;
    int df = pmod.ncoeff - 1;
    double pv = chisq_cdf_comp(df, X);

    if (plain_format(prn)) {
    pprintf(prn, "\n%s: TR^2 = %f,\n".ptr, _("Test statistic").ptr, X);
    pprintf(prn, "%s = P(%s(%d) > %f) = %f\n\n".ptr,
        _("with p-value").ptr, _("Chi-square").ptr, df, X, pv);
    } else if (rtf_format(prn)) { /* FIXME */
    pprintf(prn, "\\par \\ql\n%s: TR{\\super 2} = %f,\n".ptr, _("Test statistic").ptr,
        X);
    pprintf(prn, "%s = P(%s(%d) > %f) = %f\n\n".ptr,
        _("with p-value").ptr, _("Chi-square").ptr, df, X, pv);
    } else if (tex_format(prn)) {
    pprintf(prn, "\n%s: $TR^2$ = %f,\n".ptr, _("Test statistic").ptr, X);
    pprintf(prn, "%s = $P$($\\chi^2(%d)$ > %f) = %f\n\n".ptr,
        _("with p-value").ptr, df, X, pv);
    }
}

void print_bp_results(const(Model) * pmod, Prn * prn)
{
    double pv, X = gretl_model_get_double(pmod, "BPLM");
    int df = pmod.ncoeff - 1;

    if (na(X)) {
    return;
    }

    pv = chisq_cdf_comp(df, X);

    if (plain_format(prn)) {
    pprintf(prn, "\n%s: LM = %f,\n".ptr, _("Test statistic").ptr, X);
    pprintf(prn, "%s = P(%s(%d) > %f) = %f\n\n".ptr,
        _("with p-value").ptr, _("Chi-square").ptr, df, X, pv);
    } else if (rtf_format(prn)) { /* FIXME */
    pprintf(prn, "\\par \\ql\n%s: LM = %f,\n".ptr, _("Test statistic").ptr,
        X);
    pprintf(prn, "%s = P(%s(%d) > %f) = %f\n\n".ptr,
        _("with p-value").ptr, _("Chi-square").ptr, df, X, pv);
    } else if (tex_format(prn)) {
    pprintf(prn, "\n%s: LM = %f,\n".ptr, _("Test statistic").ptr, X);
    pprintf(prn, "%s = $P$($\\chi^2(%d)$ > %f) = %f\n\n".ptr,
        _("with p-value").ptr, df, X, pv);
    }
}

void print_HET_1_results(const(Model) * pmod, Prn * prn)
{
    double z = fabs(pmod.coeff[1]) / pmod.sderr[1];
    double pv = 2.0 * (1 - normal_cdf(z));

    if (plain_format(prn)) {
    pprintf(prn, "\n%s: HET_1 = |%f| / %f = %f,\n".ptr, _("Test statistic").ptr,
        pmod.coeff[1], pmod.sderr[1], z);
    pprintf(prn, "%s = 2 * P(z > %f) = %.3g\n\n".ptr,
        _("with p-value").ptr, z, pv);
    } else if (rtf_format(prn)) { /* FIXME */
    pprintf(prn, "\\par \\ql\n%s: HET_1 = %f,\n".ptr, _("Test statistic").ptr, z);
    pprintf(prn, "%s = 2 * P(z > %f) = %.3g\n\n".ptr,
        _("with p-value").ptr, z, pv);
    } else if (tex_format(prn)) {
    pprintf(prn, "\n%s: \verb|HET_1| = %f,\n".ptr, _("Test statistic").ptr, z);
    pprintf(prn, "%s = $2 \times P$($z$ > %f) = %f\n\n".ptr,
        _("with p-value").ptr, z, pv);
    }
}

void maybe_print_jll(const(Model) * pmod, int lldig, Prn * prn)
{
    double jll = gretl_model_get_double(pmod, "jll");

    if (!na(jll)) {
    const(char) * lp;
    char * jllstr;
    char[32] xstr;

    lp = cast(typeof(lp)) gretl_model_get_data(pmod, "log-parent");
    jllstr = g_strdup_printf(_("Log-likelihood for %s").ptr, lp);
    if (lldig > 0) {
        plain_print_double(xstr.ptr, lldig, jll, prn);
        pprintf(prn, "  (%s = %s)\n".ptr, jllstr, xstr.ptr);
    } else {
        plain_print_double(xstr.ptr, XDIGITS(pmod), jll, prn);
        pprintf(prn, "%s = %s\n\n".ptr, jllstr, xstr.ptr);
    }
    g_free(jllstr);
    }
}

bool fixed_effects_model(const(Model) * m)
{
    return (m.ci == GretlCmdIndex.panel && (m.opt & Gretlopt.OPT_F));
}

bool random_effects_model(const(Model) * m)
{
    return (m.ci == GretlCmdIndex.panel && (m.opt & Gretlopt.OPT_U));
}

bool weighted_model(const(Model) * m)
{
    return (m.ci == GretlCmdIndex.hsk || m.ci == GretlCmdIndex.arch || (m.ci == GretlCmdIndex.wls && !gretl_model_get_int(m, "wt_dummy")) || (m.ci == GretlCmdIndex.panel && (m.opt & Gretlopt.OPT_H)));
}

bool panel_ML_model(const(Model) * m)
{
    return (m.ci == GretlCmdIndex.panel && (m.opt & Gretlopt.OPT_H) && gretl_model_get_int(m, "iters"));
}

bool non_weighted_panel(const(Model) * m)
{
    return (m.ci == GretlCmdIndex.panel && !(m.opt & Gretlopt.OPT_H));
}

enum MID_STATS = 14;

struct Middletab {
    const(char) ** key; /* stats strings */
    double * val; /* stats values */
    int minus; /* variant of minus sign in use */
    int mlen; /* max length of translated string */
    int[MID_STATS] ipos; /* is i-th stat integer? */
    int nls; /* translation on? (0/1) */
    int d; /* CSV field delimiter */
    int multi; /* Using multiple precision? (0/1) */
    char[48] txt_fmt; /* format for plain text output */
}

void set_mtab_string_width(Middletab * mt)
{
    int len, badkey = 0, maxlen = 0;
    int i, j;

    for (i=0, j=0; i<7; i++, j+=2) {
    if (!na(mt.val[j])) {
        len = cast(int) g_utf8_strlen(_(mt.key[j]), -1);
        if (len > maxlen) {
        maxlen = len;
        badkey = j;
        }
        len = cast(int) g_utf8_strlen(_(mt.key[j+1]), -1);
        if (len > maxlen) {
        maxlen = len;
        badkey = j+1;
        }
    }
    }

    if (maxlen > 22) {
    fprintf(stderr, "Can't make compact model stats table -- the max\nlength translated string is %d chars, should be < 23\n".ptr,
        maxlen);
    fprintf(stderr, "offending string: '%s' ->\n '%s'\n".ptr,
        mt.key[badkey], _(mt.key[badkey]));
    }

    mt.mlen = maxlen;
}

void middletab_prepare_format(Middletab * mt, int j)
{
    const(char) * s1 = mt.key[j];
    const(char) * s2 = mt.key[j+1];
    int d1 = 0, d2 = 0;
    int len = mt.mlen;

    if (mt.nls) {
    /* calculate bytes minus glyphs */
    d1 = cast(int) (strlen(_(s1)) - g_utf8_strlen(_(s1), -1));
    d2 = cast(int) (strlen(_(s2)) - g_utf8_strlen(_(s2), -1));
    }

    if (mt.multi) {
    len = (len < 20)? 20 : len + 1;
    sprintf(mt.txt_fmt.ptr, "  %%-%ds  %%s\n  %%-%ds  %%s\n".ptr, len + d1, len + d2);
    } else if (len < 20) {
    sprintf(mt.txt_fmt.ptr, "%%-%ds%%s   %%-%ds%%s\n".ptr, 20 + d1, 20 + d2);
    } else if (len < 23) {
    sprintf(mt.txt_fmt.ptr, "%%-%ds%%s   %%-%ds%%s\n".ptr, len + d1 + 1, len + d2 + 1);
    } else {
    sprintf(mt.txt_fmt.ptr, "  %%-%ds  %%s\n  %%-%ds  %%s\n".ptr, len + d1, len + d2);
    }
}

char * print_csv(char * s, double x)
{
    if (na(x)) {
    strcpy(s, "NA");
    } else {
    sprintf(s, "%.15g".ptr, x);
    }

    return s;
}

enum MinusType {
    hyphen,
    utf,
    tex
}

enum RSQ_POS = 4;

void mtab_numstart(char * s, double x, int d7, int minus)
{
    if (x < 0) {
    if (minus == MinusType.utf) {
        strcpy(s, (d7)? " −" : "−"); /* U+2212: minus */
    } else if (minus == MinusType.tex) {
        strcpy(s, "$-$");
    } else {
        strcpy(s, (d7)? " -" : "-"); /* ASCII: use hyphen */
    }
    } else {
    strcpy(s, " ");
    }
}

bool seven_digits(double x)
{
    return x > 999999 && x < 1.0e8;
}

/* Try to pack as much as we can of a given number into a fixed width
   of 8 characters (a leading minus, if needed, is a ninth).
*/

char * print_eight(char * s, Middletab * mt, int i)
{
    double ax, x = mt.val[i];
    int d7 = (x < 0 && seven_digits(-x));
    char[16] tmp;

    if (mt.ipos[i]) {
    sprintf(s, "%9d".ptr, cast(int) x);
    return s;
    }

    if (na(x)) {
    sprintf(s, "%9s".ptr, "NA".ptr);
    return s;
    }

    if (i == RSQ_POS || i == RSQ_POS + 1) {
    /* R-squared: don't use scientific notation */
    sprintf(s, "%9.6f".ptr, x);
    return s;
    }

    mtab_numstart(s, x, d7, mt.minus);

    ax = fabs(x);

    if (d7) {
    sprintf(tmp.ptr, "%.0f".ptr, ax);
    } else if (ax < 0.00001 || ax > 99999999) {
    char * p;

    sprintf(tmp.ptr, "%#.3g".ptr, ax);
    p = strrchr(tmp.ptr, (ax < 1)? '-' : '+');
    if (p == null) {
        sprintf(tmp.ptr, "%8.6f".ptr, ax);
    } else if (strlen(p) == 4) {
        if (*(p+1) == '0') {
        memmove(p+1, p+2, 3);
        } else {
        sprintf(tmp.ptr, "%#.2g".ptr, ax);
        }
    }
    } else if (ax < 10.0) {
    sprintf(tmp.ptr, "%8.6f".ptr, ax);
    } else {
    double lx = log10(ax);
    int ldig = cast(int) (ceil(lx) + (lx == floor(lx)));

    ldig = (ldig > 7)? 7 : ldig;
    sprintf(tmp.ptr, "%8.*f".ptr, 8 - ldig - 1, ax);
    }

    strcat(s, tmp.ptr);

    if (mt.minus == MinusType.tex && strchr(s, 'e') != null) {
    tex_modify_exponent(s);
    }

    return s;
}

char * print_fifteen(char * s, double x, int minus)
{
    if (na(x)) {
    strcpy(s, " NA");
    } else if (minus == MinusType.tex) {
    char * p;

    if (x < 0) {
        sprintf(s, "$-$%.15E".ptr, -x);
    } else {
        sprintf(s, "%.15E".ptr, x);
    }

    if ((p = strstr(s, "E-")) != null) {
        char[8] tmp;

        sprintf(tmp.ptr, "E--%s".ptr, p + 2);
        strcpy(p, tmp.ptr);
    }
    } else if (minus == MinusType.utf) {
    if (x < 0) {
        sprintf(s, "−%.15E".ptr, -x);
    } else {
        sprintf(s, "% .15E".ptr, x);
    }
    } else {
    sprintf(s, "% .15E".ptr, x);
    }

    return s;
}

enum RTF_MT_ROW = "\\trowd \\trqc \\trgaph30\\trleft-30\\trrh262\\cellx2500\\cellx3800\\cellx4200\\cellx6700\\cellx8000\n";

enum RTF_MULTI_ROW = "\\trowd \\trqc \\trgaph30\\trleft-30\\trrh262\\cellx2500\\cellx6000\n";

enum RTF_MT_FMT = "\\intbl\\ql %s\\cell\\qr %s\\cell \\qc \\cell\\ql %s\\cell\\qr %s\\cell\\intbl\\row\n";

enum RTF_MULTI_FMT = "\\intbl\\ql %s\\cell\\qr %s\\cell\\intbl\\row\n";

enum TXT_MT_FMT = "%-20s%s   %-20s%s\n";

enum MiddleType {
    regular,
    transform,
    orig
}

void middle_table_row(Middletab * mt, int j, Prn * prn)
{
    const(char) * s1 = mt.key[j];
    const(char) * s2 = mt.key[j+1];
    char[48] x1;
    char[48] x2;
    int k = j + 1;

    if (tex_format(prn)) {
    if (mt.multi) {
        pprintf(prn, "%s & %s \\\\\n%s & %s \\\\\n".ptr,
            _(s1), print_fifteen(x1.ptr, mt.val[j], mt.minus),
            _(s2), print_fifteen(x2.ptr, mt.val[k], mt.minus));
    } else {
        pprintf(prn, "%s & %s & %s & %s \\\\\n".ptr,
            _(s1), print_eight(x1.ptr, mt, j),
            _(s2), print_eight(x2.ptr, mt, k));
    }
    } else if (rtf_format(prn)) {
    if (mt.multi) {
        pputs(prn, RTF_MULTI_ROW.ptr);
        pprintf(prn, RTF_MULTI_FMT.ptr, _(s1),
            print_fifteen(x1.ptr, mt.val[j], mt.minus));
        pputs(prn, RTF_MULTI_ROW.ptr);
        pprintf(prn, RTF_MULTI_FMT.ptr, _(s2),
            print_fifteen(x1.ptr, mt.val[k], mt.minus));
    } else {
        pputs(prn, RTF_MT_ROW.ptr);
        pprintf(prn, RTF_MT_FMT.ptr,
            _(s1), print_eight(x1.ptr, mt, j),
            _(s2), print_eight(x2.ptr, mt, k));
    }
    } else if (csv_format(prn)) {
    pprintf(prn, "\"%s\"%c%s%c\"%s\"%c%s\n".ptr,
        _(s1), mt.d, print_csv(x1.ptr, mt.val[j]), mt.d,
        _(s2), mt.d, print_csv(x2.ptr, mt.val[k]));
    } else {
    if (mt.nls || mt.multi) {
        middletab_prepare_format(mt, j);
    }
    if (mt.multi) {
        pprintf(prn, mt.txt_fmt.ptr,
            _(s1), print_fifteen(x1.ptr, mt.val[j], mt.minus),
            _(s2), print_fifteen(x2.ptr, mt.val[k], mt.minus));
    } else {
        pprintf(prn, mt.txt_fmt.ptr,
            _(s1), print_eight(x1.ptr, mt, j),
            _(s2), print_eight(x2.ptr, mt, k));
    }
    }
}

int string_is_translated(const(char) * s)
{
    if (strcmp(s, "Hannan-Quinn") && strcmp(s, "rho")) {
    return strcmp(s, _(s));
    } else {
    /* give the benefit of the doubt */
    return 1;
    }
}

/* If we haven't found a translation for a new-style short string,
   use the corresponding old-style long string instead.
*/

void maybe_remedy_translations(const(char) ** S, int n)
{
    static const(char) *[8] old_key = ["Mean of dependent variable".ptr, "Standard deviation of dep. var.".ptr, "Sum of squared residuals".ptr, "Standard error of the regression".ptr, "Unadjusted R-squared".ptr, "Akaike information criterion".ptr, "Schwarz Bayesian criterion".ptr, "Hannan-Quinn Information Criterion".ptr];
    int i;

    for (i=0; i<n; i++) {
    if (old_key[i] != null && !string_is_translated(S[i])) {
        /* new-style string is not translated */
        if (strcmp(old_key[i], _(old_key[i]))) {
        /* but the old-style one is, so we'll use it */
        S[i] = old_key[i];
        }
    }
    }
}

/* order of keys for model stats in middle table */

enum MiddleKey {
    ybar,
    sdy,
    ssr,
    ser,
    rsq,
    r22,
    fst,
    pvf,
    lnl,
    aka,
    scz,
    hq,
    rho,
    dw
}

/* print the block of statistics that appears beneath of the
   table of coefficients, standard errors, etc.
*/

void print_middle_table(Model * pmod, Prn * prn, int code)
{
    const char *note = /* TRANSLATORS: please do not translate literally: this is for your use in describing locale-specific abbreviations. It is OK to leave it untranslated, in which case it will not be printed. */ "note on model statistics abbreviations here".ptr;
    int rtf = rtf_format(prn);
    int tex = tex_format(prn);
    int csv = csv_format(prn);
    char * teststr = null;
    static const(char) *[14] key = ["Mean dependent var".ptr, "S.D. dependent var".ptr, "Sum squared resid".ptr, "S.E. of regression".ptr, "R-squared".ptr, "Adjusted R-squared".ptr, "F-statistic".ptr, "P-value(F)".ptr, "Log-likelihood".ptr, "Akaike criterion".ptr, "Schwarz criterion".ptr, "Hannan-Quinn".ptr, "rho".ptr, "Durbin-Watson".ptr];
    double[MID_STATS] val = [
    pmod.ybar,
    pmod.sdy,
    pmod.ess,
    pmod.sigma,
    pmod.rsq,
    pmod.adjrsq,
    pmod.fstt,
    NADBL,
    pmod.lnL,
    pmod.criterion[InfoCriterion.aic],
    pmod.criterion[InfoCriterion.bic],
    pmod.criterion[InfoCriterion.hqc],
    pmod.rho,
    pmod.dw
    ];
    Middletab mtab;
    int i, j;

    mtab.mlen = 0;
    mtab.d = 0;
    mtab.minus = MinusType.hyphen;
    for (i=0; i<MID_STATS; i++) {
    mtab.ipos[i] = 0;
    }
    mtab.nls = 0;
    mtab.multi = (pmod.ci == GretlCmdIndex.mpols);

    mtab.nls = doing_nls();
    if (mtab.nls) {
    maybe_remedy_translations(key.ptr, MID_STATS);
    }

    if (tex) {
    /* some special strings for TeX output */
    mtab.minus = MinusType.tex;
    key[MiddleKey.rsq] = "$R^2$";
    key[MiddleKey.r22] = "Adjusted $R^2$".ptr;
    key[MiddleKey.pvf] = "P-value($F$)".ptr;
    key[MiddleKey.hq]  = "Hannan--Quinn";
    key[MiddleKey.rho] = "$\\hat{\\rho}$";
    key[MiddleKey.dw]  = "Durbin--Watson";
    } else if (gretl_print_has_minus(prn)) {
    /* print a 'real' minus sign? */
    mtab.minus = MinusType.utf;
    }

    if (pmod.aux == ModelAuxCode.vecm || pmod.aux == ModelAuxCode.coint) {
    /* VECM equation or Engle-Granger test: suppress F-test */
    val[MiddleKey.fst] = val[MiddleKey.pvf] = NADBL;
    } else if (!na(pmod.fstt) && pmod.dfd > 0) {
    /* format F-stat and get its p-value */
    int nc = gretl_model_get_int(pmod, "n_clusters");
    int dfd = nc > 1 ? nc - 1 : pmod.dfd;

    if (tex) {
        teststr = g_strdup_printf("$F(%d, %d)$".ptr, pmod.dfn, dfd);
    } else if ((pmod.ci == GretlCmdIndex.panel || pmod.ci == GretlCmdIndex.logistic) &&
           (pmod.opt & Gretlopt.OPT_F)) {
        teststr = g_strdup_printf(_("LSDV F(%d, %d)").ptr, pmod.dfn, dfd);
    } else {
        teststr = g_strdup_printf("F(%d, %d)".ptr, pmod.dfn, dfd);
    }
    key[MiddleKey.fst] = teststr;
    val[MiddleKey.pvf] = snedecor_cdf_comp(pmod.dfn, dfd, pmod.fstt);
    } else if (!na(pmod.chisq)) {
    /* alternative: chi-square and its p-value */
    teststr = g_strdup_printf("%s(%d)".ptr, _("Chi-square").ptr, pmod.dfn);
    key[MiddleKey.fst] = teststr;
    val[MiddleKey.fst] = pmod.chisq;
    key[MiddleKey.pvf] = "p-value".ptr;
    val[MiddleKey.pvf] = chisq_cdf_comp(pmod.dfn, val[6]);
    }

    /* special variants of R-squared */
    if (gretl_model_get_int(pmod, "uncentered")) {
    key[MiddleKey.rsq] = (tex)? "Uncentered $R^2$".ptr :
        "Uncentered R-squared".ptr; /* 22: */
    key[MiddleKey.r22] = (tex)? "Centered $R^2$".ptr :
        "Centered R-squared".ptr;  /* 22: */
    val[MiddleKey.r22] = gretl_model_get_double(pmod, "centered-R2");
    } else if (COUNT_MODEL(cast(GretlCmdIndex) pmod.ci)) {
    key[MiddleKey.rsq] = (tex)? "McFadden $R^2$".ptr :
        "McFadden R-squared".ptr;  /* 22: McFadden pseudo-R^2 */
    } else if (binary_model(pmod) || ordered_model(pmod)) {
    if (pmod.opt & Gretlopt.OPT_S) {
        key[MiddleKey.rsq] = (tex)? "Estrella $R^2$".ptr :
        "Estrella R-squared".ptr;  /* 22: Estrella pseudo-R^2 */
    } else {
        key[MiddleKey.rsq] = (tex)? "McFadden $R^2$".ptr :
        "McFadden R-squared".ptr;  /* 22: McFadden pseudo-R^2 */
    }
    } else if ((pmod.ci == GretlCmdIndex.panel || pmod.ci == GretlCmdIndex.logistic) &&
           (pmod.opt & Gretlopt.OPT_F)) {
    /* 22: panel, fixed effects */
    key[MiddleKey.rsq] = (tex)? "LSDV $R^2$".ptr : "LSDV R-squared".ptr;
    key[MiddleKey.r22] = (tex)? "Within $R^2$".ptr : "Within R-squared".ptr;
    }

    if (pmod.ci == GretlCmdIndex.dpanel) {
    for (i=0; i<MID_STATS; i++) {
        if (i < MiddleKey.ssr || i > MiddleKey.ser) {
        val[i] = NADBL;
        }
    }
    } else if (pmod.ci == GretlCmdIndex.ivreg && (pmod.opt & Gretlopt.OPT_G)) {
    /* IVREG via GMM */
    for (i=MiddleKey.ssr; i<MID_STATS; i++) {
        val[i] = NADBL;
    }
    } else if (pmod.aux == ModelAuxCode.sys) {
    /* only dep. var. stats, SSR and SE, unless LIML */
    if (liml_equation(pmod) && !gretl_model_get_int(pmod, "restricted")) {
        for (i=4; i<MID_STATS; i++) {
        if (i < MiddleKey.lnl || i > MiddleKey.aka) {
            val[i] = NADBL;
        }
        }
        key[MiddleKey.aka] = "Smallest eigenvalue".ptr; /* 22: */
        val[MiddleKey.aka] = gretl_model_get_double(pmod, "lmin");
    } else {
        for (i=MiddleKey.fst; i<MID_STATS; i++) {
        val[i] = NADBL;
        }
    }
    } else if (pmod.ci == GretlCmdIndex.arma) {
    key[MiddleKey.ssr] = "Mean of innovations".ptr; /* 22: Mean of ARMA innovations */
    val[MiddleKey.ssr] = gretl_model_get_double(pmod, "mean_error");
    key[MiddleKey.ser] = "S.D. of innovations".ptr; /* 22: Std. dev. of ARMA innovations */
    } else if (pmod.ci == GretlCmdIndex.lad || pmod.ci == GretlCmdIndex.quantreg) {
    key[MiddleKey.ybar] = "Median depend. var".ptr;  /* 22: Median of dependent variable */
    val[MiddleKey.ybar] = gretl_model_get_double(pmod, "ymedian");
    key[MiddleKey.ssr] = "Sum absolute resid".ptr;  /* 22: Sum of absolute residuals */
    val[MiddleKey.ssr] = gretl_model_get_double(pmod, "ladsum");
    key[MiddleKey.ser] = "Sum squared resid".ptr;
    val[MiddleKey.ser] = pmod.ess;
    for (i=MiddleKey.rsq; i<MID_STATS; i++) {
        if (i < MiddleKey.lnl || i > MiddleKey.hq) {
        val[i] = NADBL;
        }
    }
    } else if (logit_probit_model(pmod)) {
    val[MiddleKey.ssr] = val[MiddleKey.ser] = NADBL;
    val[MiddleKey.fst] = val[MiddleKey.pvf] = NADBL;
    val[MiddleKey.rho] = val[MiddleKey.dw] = NADBL;
    } else if (pmod.ci == GretlCmdIndex.biprobit) {
    val[MiddleKey.ssr] = val[MiddleKey.ser] = NADBL;
    val[MiddleKey.fst] = val[MiddleKey.pvf] = NADBL;
    val[MiddleKey.rho] = val[MiddleKey.dw] = NADBL;
    } else if (pmod.ci == GretlCmdIndex.heckit) {
    key[MiddleKey.ssr] = (tex)? "$\\hat{\\sigma}$" : "sigma".ptr;
    val[MiddleKey.ssr] = pmod.sigma;
    key[MiddleKey.ser] = (tex)? "$\\hat{\\rho}$" : "rho".ptr;
    val[MiddleKey.ser] = pmod.rho;
    for (i=MiddleKey.rsq; i<MID_STATS; i++) {
        /* may add an R-squared (item 4)? */
        if (i < MiddleKey.lnl || i > MiddleKey.hq) {
        val[i] = NADBL;
        }
    }
    } else if (panel_ML_model(pmod) || pmod.ci == GretlCmdIndex.garch) {
    for (i=MiddleKey.ssr; i<MID_STATS; i++) {
        if (i < MiddleKey.lnl || i > MiddleKey.hq) {
        val[i] = NADBL;
        }
    }
    } else if (pmod.ci == GretlCmdIndex.mle || ordered_model(pmod)) {
    for (i=0; i<MID_STATS; i++) {
        if (i < MiddleKey.lnl || i > MiddleKey.hq) {
        val[i] = NADBL;
        }
    }
    } else if (pmod.ci != GretlCmdIndex.var && pmod.aux != ModelAuxCode.vecm &&
           !na(pmod.rho) && gretl_model_get_int(pmod, "ldepvar")) {
        int err = 0;
        double h = durbins_h(pmod, &err);

    key[MiddleKey.dw] = (tex)? "Durbin's $h$".ptr : "Durbin's h".ptr;
    val[MiddleKey.dw] = err ? NADBL : h;
    } else if (intreg_model(pmod)) {
    for (i=0; i<MID_STATS; i++) {
        if (i < MiddleKey.fst || i > MiddleKey.hq) {
        val[i] = NADBL;
        }
    }
    }

    if (code == MiddleType.transform) {
    /* transformed data: don't print mean, s.d. of dep. var. */
    val[MiddleKey.ybar] = val[MiddleKey.sdy] = NADBL;
    } else if (code == MiddleType.orig) {
    /* print a limited range of stats */
    val[MiddleKey.ssr] = gretl_model_get_double(pmod, "ess_orig");
    val[MiddleKey.ser] = gretl_model_get_double(pmod, "sigma_orig");
    for (i=MiddleKey.rsq; i<MID_STATS; i++) {
        val[i] = NADBL;
    }
    if (pmod.ci == GretlCmdIndex.logistic) {
        val[MiddleKey.lnl] = gretl_model_get_double(pmod, "jll");
        val[MiddleKey.aka] = gretl_model_get_double(pmod, "jaic");
    }
    }

    /* start the table */
    if (tex) {
    pputs(prn, "\\vspace{1ex}\n");
    if (mtab.multi) {
        pputs(prn, "\\begin{tabular}{lr}\n");
    } else {
        pputs(prn, "\\begin{tabular}{lrlr}\n");
    }
    } else if (rtf) {
    pputs(prn, "\\par\n{");
    } else if (csv) {
    mtab.d = prn_delim(prn);
    }

    mtab.key = key.ptr;
    mtab.val = val.ptr;
    strcpy(mtab.txt_fmt.ptr, TXT_MT_FMT.ptr);

    if (plain_format(prn) && mtab.nls) {
    set_mtab_string_width(&mtab);
    }

    /* print the various statistics */
    for (i=0, j=0; i<7; i++, j+=2) {
    if (!na(val[j])) {
        middle_table_row(&mtab, j, prn);
    }
    }

    if (tex) {
    pputs(prn, "\\end{tabular}\n\n");
    } else if (rtf) {
    pputs(prn, "}\n\n"); /* close RTF table */
    }

    if (plain_format(prn) && strcmp(note, _(note)) &&
    !string_is_blank(_(note))) {
    pputs(prn, _(note));
    pputc(prn, '\n');
    }

    if (!tex && !rtf) {
    pputc(prn, '\n');
    }

    g_free(teststr);
}

void print_model_iter_info(const(Model) * pmod, Prn * prn)
{
    int iters = gretl_model_get_int(pmod, "iters");

    if (iters > 0) {
    /* pputc(prn, '\n'); */
    pprintf(prn, _("Convergence achieved after %d iterations\n").ptr, iters);
    } else {
    int fncount = gretl_model_get_int(pmod, "fncount");
    int grcount = gretl_model_get_int(pmod, "grcount");

    if (fncount > 0) {
        pputc(prn, '\n');
        pprintf(prn, _("Function evaluations: %d\n").ptr, fncount);
        pprintf(prn, _("Evaluations of gradient: %d\n").ptr, grcount);
    }
    }
}

void aux_print_info_criteria(const(Model) * pmod, Prn * prn)
{
    static const(char) *[3] istrs = ["AIC".ptr, "BIC".ptr, "HQC".ptr];
    int i, n = 0;

    for (i=0; i<InfoCriterion.max; i++) {
    if (!na(pmod.criterion[i])) {
        if (n > 0) {
        pputc(prn, ' ');
        }
        pprintf(prn, "  %s: %g".ptr, _(istrs[i]), pmod.criterion[i]);
        n++;
    }
    }

    if (n > 0) {
    if (gretl_model_get_int(pmod, "eg-resids")) {
        pputc(prn, '\n');
    } else {
        pputs(prn, "\n\n");
    }
    }
}

void set_csv_delim(Prn * prn)
{
    char[4] test;

    sprintf(test.ptr, "%.1f".ptr, 1.0);

    if (test[1] == ',') {
    gretl_print_set_delim(prn, '\t');
    } else {
    gretl_print_set_delim(prn, ',');
    }
}

/**
 * printmodel:
 * @pmod: pointer to gretl model.
 * @dset: data information struct.
 * @opt: may contain %OPT_O to print covariance matrix, %OPT_S
 * to get a "simple" print (just coefficients and standard
 * errors), %OPT_P for printing model in the context of panel
 * diagnostics only.
 * @prn: gretl printing struct.
 *
 * Print to @prn the estimates in @pmod plus associated statistics.
 *
 * Returns: 0 on success, 1 if some of the values to print were %NaN.
 */

int printmodel(Model * pmod, const(DatasetRaw) * dset, Gretlopt opt, Prn * prn)
{
    int gotnan = 0;

    if (prn == null || (opt & Gretlopt.OPT_Q)) {
    return 0;
    }

    if (csv_format(prn)) {
    set_csv_delim(prn);
    }

    if (plain_format(prn)) {
    print_model_iter_info(pmod, prn);
    }

    if (gretl_is_between_model(pmod) && pmod.dataset != null) {
    dset = pmod.dataset;
    }

    model_format_start(prn);

    if (!(opt & Gretlopt.OPT_P)) {
    print_model_heading(pmod, dset, opt, prn);
    }

    if (plain_format(prn)) {
    gotnan = plain_print_coefficients(pmod, dset, prn);
    } else {
    gotnan = alt_print_coefficients(pmod, dset, prn);
    }

    print_coeff_table_end(pmod, prn);

    if (pmod.aux == ModelAuxCode.df || pmod.aux == ModelAuxCode.adf ||
    pmod.aux == ModelAuxCode.kpss) {
    if (!gretl_model_get_int(pmod, "dfgls")) {
        aux_print_info_criteria(pmod, prn);
    }
    goto close_format;
    }

    if (pmod.aux == ModelAuxCode.arch || pmod.aux == ModelAuxCode.reset) {
    goto close_format;
    } else if (opt & Gretlopt.OPT_P) {
    goto close_format;
    }

    /* other auxiliary regressions */

    if (pmod.aux == ModelAuxCode.het_1) {
    rsqline(pmod, prn);
    print_HET_1_results(pmod, prn);
    goto close_format;
    } else if (pmod.aux == ModelAuxCode.white) {
    rsqline(pmod, prn);
    print_whites_results(pmod, prn);
    goto close_format;
    } else if (pmod.aux == ModelAuxCode.bp) {
    rssline(pmod, prn);
    print_bp_results(pmod, prn);
    goto close_format;
    } else if (pmod.aux == ModelAuxCode.sq ||
           pmod.aux == ModelAuxCode.log ||
           pmod.aux == ModelAuxCode.ar ||
           pmod.aux == ModelAuxCode.add) {
    rsqline(pmod, prn);
    goto close_format;
    } else if (pmod.aux == ModelAuxCode.comfac) {
    ssrline(pmod, prn);
    goto close_format;
    }

    if ((pmod.ci == GretlCmdIndex.ols || pmod.ci == GretlCmdIndex.mpols) && (opt & Gretlopt.OPT_S)) {
    /* --simple-print */
    if (!na(pmod.ess) && !na(pmod.rsq) && plain_format(prn)) {
        int uc = gretl_model_get_int(pmod, "uncentered");
            const(char) * s2 = uc ? _("Uncentered R-squared").ptr : _("R-squared").ptr;

            if (pmod.ci == GretlCmdIndex.mpols) {
                const(char) * s1 = _("Sum squared resid").ptr;
                int n1 = cast(int) strlen(s1);
                int n2 = cast(int) strlen(s2);
                int n = 1 + max(n1, n2);
                char[32] tmp;

                pprintf(prn, "  %-*s %s\n".ptr, n, s1,
                        print_fifteen(tmp.ptr, pmod.ess, 0));
                pprintf(prn, "  %-*s %s\n\n".ptr, n, s2,
                        print_fifteen(tmp.ptr, pmod.rsq, 0));
            } else {
                pprintf(prn, "%s = %g, %s = %f\n\n".ptr, _("SSR").ptr, pmod.ess,
                        s2, pmod.rsq);
            }
    }
    goto close_format;
    }

    if (weighted_model(pmod)) {
    alternate_stats_message(AltStats.wtd_stats, prn);
    print_middle_table(pmod, prn, MiddleType.transform);
    alternate_stats_message(AltStats.orig_stats, prn);
    print_middle_table(pmod, prn, MiddleType.orig);
    goto pval_max;
    }

    if (pmod.ci == GretlCmdIndex.logistic) {
    alternate_stats_message(AltStats.transform_stats, prn);
    print_middle_table(pmod, prn, MiddleType.transform);
    alternate_stats_message(AltStats.orig_stats, prn);
    print_middle_table(pmod, prn, MiddleType.orig);
    goto pval_max;
    }

    if (pmod.ci == GretlCmdIndex.ar || pmod.ci == GretlCmdIndex.ar1) {
    alternate_stats_message(AltStats.rhodiff_stats, prn);
    print_middle_table(pmod, prn, MiddleType.transform);
    alternate_stats_message(AltStats.orig_stats, prn);
    print_middle_table(pmod, prn, MiddleType.orig);
    goto pval_max;
    }

    /* print table of model stats */
    if (pmod.ci != GretlCmdIndex.gmm) {
    print_middle_table(pmod, prn, MiddleType.regular);
    }

    /* random effects probit */
    if (re_probit_model(pmod)) {
    print_probit_rho(pmod, prn);
    }

    /* additional stats/info for some cases */
    if ((pmod.aux == ModelAuxCode.sys && liml_equation(pmod)) ||
    liml_model(pmod)) {
    print_liml_equation_data(pmod, prn);
    } else if (pmod.ci == GretlCmdIndex.arma) {
    print_arma_roots(pmod, prn);
    } else if (pmod.ci == GretlCmdIndex.garch) {
    garch_variance_line(pmod, prn);
    } else if (pmod.ci == GretlCmdIndex.heckit) {
    print_heckit_stats(pmod, prn);
    } else if (random_effects_model(pmod)) {
    panel_variance_lines(pmod, prn);
    } else if (gmm_model(pmod)) {
    print_GMM_stats(pmod, prn);
    } else if (pmod.ci == GretlCmdIndex.dpanel) {
    print_DPD_stats(pmod, prn);
    } else if (logit_probit_model(pmod)) {
    if (!pmod.aux) {
        logit_probit_stats(pmod, prn);
    }
    } else if (tsls_model(pmod) && plain_format(prn)) {
    addconst_message(pmod, prn);
    } else if (intreg_model(pmod)) {
    print_intreg_info(pmod, dset, prn);
    } else if (pmod.ci == GretlCmdIndex.poisson) {
    print_overdisp_test(pmod, prn);
    } else if (pmod.ci == GretlCmdIndex.duration) {
    print_duration_alpha(pmod, prn);
    } else if (pmod.ci == GretlCmdIndex.nls ||
           (pmod.ci == GretlCmdIndex.midasreg &&
        !gretl_model_get_int(pmod, "umidas"))) {
    print_GNR_info(pmod, prn);
    }

    /* FIXME alternate R^2 measures (within, centered) */

    if (plain_format(prn) && !pmod.aux) {
    maybe_print_jll(pmod, 0, prn);
    }

 pval_max:

    if (plain_format(prn) && pmod.ci != GretlCmdIndex.mle && pmod.ci != GretlCmdIndex.panel &&
    pmod.ci != GretlCmdIndex.arma && pmod.ci != GretlCmdIndex.nls && pmod.ci != GretlCmdIndex.gmm &&
    pmod.ci != GretlCmdIndex.lad && pmod.ci != GretlCmdIndex.heckit && pmod.ci != GretlCmdIndex.midasreg &&
    pmod.ci != GretlCmdIndex.dpanel && pmod.ci != GretlCmdIndex.garch && pmod.ci != GretlCmdIndex.duration &&
    !ordered_model(pmod) && !multinomial_model(pmod) &&
    !COUNT_MODEL(cast(GretlCmdIndex) pmod.ci) && !intreg_model(pmod) &&
    pmod.ci != GretlCmdIndex.biprobit && pmod.ci != GretlCmdIndex.quantreg && !pmod.aux) {
    pval_max_line(pmod, dset, prn);
    }

 close_format:

    if (opt & Gretlopt.OPT_V) {
    ols_print_anova(pmod, prn);
    }

    if (opt & Gretlopt.OPT_O) {
    outcovmx(pmod, dset, prn);
    }

    if (any_tests(pmod) && !(opt & Gretlopt.OPT_S)) {
    print_model_tests(pmod, prn);
    }

    model_format_end(prn);

    if (gotnan) {
    pmod.errcode = GretlError.nan;
    }

    return gotnan;
}

void print_pval_str(double pval, char * str)
{
    if (pval < .0001) {
    sprintf(str, "<%.4f".ptr, 0.0001);
    } else {
    sprintf(str, "%.4f".ptr, pval);
    }
}

/* used only for "alternate" model-output formats */

int prepare_model_coeff(const(Model) * pmod, const(DatasetRaw) * dset, int i, int adfnum, ModelCoeff * mc, Prn * prn)
{
    int gotnan = 0;

    model_coeff_init(mc);

    mc.show_pval = !binary_model(pmod) || (pmod.opt & Gretlopt.OPT_P);

    if (tex_format(prn)) {
    make_tex_coeff_name(pmod, dset, i, mc.name.ptr);
    } else {
    gretl_model_get_param_name(pmod, dset, i, mc.name.ptr);
    }

    if (na(pmod.coeff[i])) {
    gotnan = 1;
    } else {
    mc.b = pmod.coeff[i];
    }

    if (!na(pmod.sderr[i])) {
    mc.se = pmod.sderr[i];
    }

    if (pmod.ci == GretlCmdIndex.duration && i > 0 && !strcmp(mc.name.ptr, "sigma")) {
    mc.tval = NADBL;
    mc.pval = NADBL;
    mc.show_tval = 0;
    mc.show_pval = 0;
    return gotnan;
    }

    if (!na(mc.b) && !na(mc.se) && mc.se > 0.0) {
    mc.tval = mc.b / mc.se;
    if (na(mc.tval)) {
        mc.tval = NADBL;
    }
    }

    if (mc.show_pval && !na(mc.tval)) {
    if (i == adfnum) {
        mc.pval = gretl_model_get_double(pmod, "dfpval");
        mc.df_pval = 1;
    } else {
        mc.pval = model_coeff_pval(pmod, mc.tval);
    }
    }

    if (!gotnan && !mc.show_pval &&
    binary_model(pmod) &&
    pmod.list[i+2] != 0) {
    double * slopes = cast(double *) gretl_model_get_data(pmod, "slopes");

    if (slopes != null) {
        mc.slope = slopes[i];
    }
    }

    if (pmod.ci == GretlCmdIndex.mpols) {
    mc.multi = 1;
    }

    return gotnan;
}

void rtf_print_double(double xx, Prn * prn)
{
    char[32] numstr;

    xx = screen_zero(xx);

    if (xx < 0 && gretl_print_has_minus(prn)) {
    sprintf(numstr.ptr, "%#.*g".ptr, get_gretl_digits(), fabs(xx));
    gretl_fix_exponent(numstr.ptr);
    pprintf(prn, " \\qc −%s\\cell".ptr, numstr.ptr);
    } else {
    sprintf(numstr.ptr, "%#.*g".ptr, get_gretl_digits(), xx);
    gretl_fix_exponent(numstr.ptr);
    pprintf(prn, " \\qc %s\\cell".ptr, numstr.ptr);
    }
}

/* Escape a string for safe embedding in generated RTF text: double
   any literal backslash or brace, which are RTF's own control-word/
   group-delimiter characters. Without this, a coefficient or
   statistic name containing e.g. "\field{\*\fldinst ...}" would be
   emitted as live RTF markup rather than literal text. Returns a
   newly allocated string; the caller must free() it.
*/

char * rtf_escape_new(const(char) * src)
{
    const(char) * s;
    char * ret, p;
    int len = 0;

    if (src == null) {
    return gretl_strdup("");
    }

    for (s=src; *s; s++) {
    if (*s == '\\' || *s == '{' || *s == '}') {
        len++;
    }
    len++;
    }

    ret = cast(typeof(ret)) calloc(1, len + 1);
    if (ret == null) {
    return null;
    }

    p = ret;
    for (s=src; *s; s++) {
    if (*s == '\\' || *s == '{' || *s == '}') {
        *p++ = '\\';
    }
    *p++ = *s;
    }

    return ret;
}

/* Escape a string for safe embedding in a quoted CSV field: double
   any embedded double-quote (RFC 4180), and if the value would
   start with a character a spreadsheet application treats as a
   formula trigger ('=','+','-','@', tab, or CR), prefix it with a
   benign leading apostrophe to neutralize that. Returns a newly
   allocated string; the caller must free() it.
*/

char * csv_escape_new(const(char) * src)
{
    const(char) * s;
    char * ret, p;
    int len = 0;
    int prefix = 0;

    if (src == null) {
    return gretl_strdup("");
    }

    if (*src == '=' || *src == '+' || *src == '-' || *src == '@' ||
    *src == '\t' || *src == '\r') {
    prefix = 1;
    len++;
    }

    for (s=src; *s; s++) {
    if (*s == '"') {
        len++;
    }
    len++;
    }

    ret = cast(typeof(ret)) calloc(1, len + 1);
    if (ret == null) {
    return null;
    }

    p = ret;
    if (prefix) {
    *p++ = '\'';
    }
    for (s=src; *s; s++) {
    if (*s == '"') {
        *p++ = '"';
    }
    *p++ = *s;
    }

    return ret;
}

void rtf_print_coeff(const(ModelCoeff) * mc, Prn * prn)
{
    char * ename;

    if (!na(mc.lo)) {
    pputs(prn, RTF_INTVL_ROW.ptr);
    } else {
    pputs(prn, RTF_COEFF_ROW.ptr);
    }

    ename = rtf_escape_new(mc.name.ptr);
    pprintf(prn, "\\ql %s\\cell".ptr, ename);
    free(ename);

    if (na(mc.b)) {
    pprintf(prn, " \\qc %s\\cell".ptr, _("undefined").ptr);
    } else {
    rtf_print_double(mc.b, prn);
    }

    if (!na(mc.lo) && !na(mc.hi)) {
    rtf_print_double(mc.lo, prn);
    rtf_print_double(mc.hi, prn);
    goto rtf_finish;
    }

    if (na(mc.se)) {
    pprintf(prn, " \\qc %s\\cell".ptr, _("undefined").ptr);
    pprintf(prn, " \\qc %s\\cell".ptr, _("undefined").ptr);
    pprintf(prn, " \\qc %s\\cell".ptr, _("undefined").ptr);
    goto rtf_finish;
    }

    rtf_print_double(mc.se, prn);

    if (!na(mc.tval)) {
    if (mc.tval < 0 && gretl_print_has_minus(prn)) {
        pprintf(prn, " \\qc −%#.4g\\cell".ptr, -mc.tval);
    } else {
        pprintf(prn, " \\qc %#.4g\\cell".ptr, mc.tval);
    }
    } else if (!mc.show_tval) {
    pputs(prn, " \\qc \\cell");
    } else {
    pprintf(prn, " \\qc %s\\cell".ptr, _("undefined").ptr);
    }

    if (!na(mc.slope)) {
    rtf_print_double(mc.slope, prn);
    } else if (mc.show_pval) {
    if (na(mc.pval)) {
        if (mc.df_pval) {
        pprintf(prn, " \\qc %s\\cell".ptr, _("unknown").ptr);
        } else {
        pprintf(prn, " \\qc %s\\cell".ptr, _("undefined").ptr);
        }
    } else {
        char[16] pvalstr;

        print_pval_str(mc.pval, pvalstr.ptr);
        pprintf(prn, " \\qc %s\\cell".ptr, pvalstr.ptr);

        if (mc.pval < 0.01) {
        pputs(prn, " \\ql ***\\cell");
        } else if (mc.pval < 0.05) {
        pputs(prn, " \\ql **\\cell");
        } else if (mc.pval < 0.10) {
        pputs(prn, " \\ql *\\cell");
        } else {
        pputs(prn, " \\ql \\cell");
        }
    }
    } else {
    pputs(prn, " \\qc \\cell \\ql \\cell");
    }

 rtf_finish:

    pputs(prn, " \\intbl \\row\n");
}

void csv_print_coeff(const(ModelCoeff) * mc, Prn * prn)
{
    char * ename = csv_escape_new(mc.name.ptr);
    char d = prn_delim(prn);

    pprintf(prn, "\"%s\"".ptr, ename);
    free(ename);

    if (na(mc.b)) {
    pprintf(prn, "%c\"%s\"".ptr, d, _("undefined").ptr);
    } else {
    pprintf(prn, "%c%.15g".ptr, d, mc.b);
    }

    if (!na(mc.lo) && !na(mc.hi)) {
    pprintf(prn, "%c%.15g".ptr, d, mc.lo);
    pprintf(prn, "%c%.15g\n".ptr, d, mc.hi);
    return;
    }

    /* get out if std error is undefined */
    if (na(mc.se)) {
    pprintf(prn, "%c\"%s\"\n".ptr, d, _("undefined").ptr);
    return;
    }

    pprintf(prn, "%c%.15g".ptr, d, mc.se);

    if (!na(mc.tval)) {
    pprintf(prn, "%c%.15g".ptr, d, mc.tval);
    } else if (!mc.show_tval) {
    pputs(prn, "%c\n");
    } else {
    pprintf(prn, "%c\"%s\"\n".ptr, d, _("undefined").ptr);
    }

    if (!na(mc.slope)) {
    /* slope for binary models */
    pprintf(prn, "%c%.15g".ptr, d, mc.slope);
    } else if (mc.show_pval) {
    if (na(mc.pval)) {
        if (mc.df_pval) {
        pprintf(prn, "%c\"%s\"\n".ptr, d, _("unknown").ptr);
        } else {
        pprintf(prn, "%c\"%s\"\n".ptr, d, _("undefined").ptr);
        }
    } else {
        pprintf(prn, "%c%.15g".ptr, d, mc.pval);
    }
    }

    pputc(prn, '\n');
}

void alt_print_coeff(const(ModelCoeff) * mc, Prn * prn)
{
    if (rtf_format(prn)) {
    rtf_print_coeff(mc, prn);
    } else if (tex_format(prn)) {
    tex_print_coeff(mc, prn);
    } else if (csv_format(prn)) {
    csv_print_coeff(mc, prn);
    }
}

char * tex_handle_underscores(const(char) * s, int n)
{
    char * ret = cast(char *) malloc(strlen(s) + 1 + n);
    int i, j = 0;

    if (ret != null) {

    for (i=0; s[i] != '\0'; i++) {
        if (s[i] == '_' && (i == 0 || s[i-1] != '\\')) {
        ret[j++] = '\\';
        ret[j++] = '_';
        } else {
        ret[j++] = s[i];
        }
    }

    ret[j] = '\0';
    }

    return ret;
}

int has_raw_underscores(const(char) * s)
{
    int i, n = 0;

    for (i=0; s[i] != '\0'; i++) {
    if (s[i] == '_' && (i == 0 || s[i-1] != '\\')) {
        n++;
    }
    }

    return n;
}

void print_coeff_separator(const(char) * s, int width, Prn * prn)
{
    int havestr = (s != null && *s != '\0');
    int vspace = (width > 0);

    if (plain_format(prn)) {
    if (havestr) {
        if (vspace) {
        pputs(prn, "\n  ");
        } else {
        pputs(prn, "  ");
        }
        if (width > 0) {
        print_centered(_(s), width, prn);
        } else {
        pputs(prn, _(s));
        }
        pputc(prn, '\n');
    }
    if (vspace) {
        pputc(prn, '\n');
    }
    } else if (tex_format(prn)) {
    if (havestr) {
        int n = has_raw_underscores(s);
        char * repl = null;

        if (n > 0) {
        repl = tex_handle_underscores(s, n);
        }
        pputs(prn, "\\\\ [-8pt]\n");
        pprintf(prn, "\\multicolumn{%d}{c}{%s} \\\\[1ex]\n".ptr, width,
            repl != null ? _(repl) : _(s));
        free(repl);
    } else {
        pputs(prn, "\\\\ \n");
    }
    } else if (rtf_format(prn)) {
    pputs(prn, RTF_MULTICOL.ptr);
    if (havestr) {
        pprintf(prn, "\\qc %s".ptr, _(s));
    }
    pputs(prn, "\\cell\\intbl\\row\n");
    } else if (csv_format(prn)) {
    if (havestr) {
        pprintf(prn, "\n\"%s\"\n".ptr, _(s));
    } else {
        pputc(prn, '\n');
    }
    }
}

void print_coeff_left_string(const(char) * s, Prn * prn)
{
    if (plain_format(prn)) {
    pprintf(prn, " %s:\n".ptr, s);
    } else if (tex_format(prn)) {
    char[48] tmp;

    tex_escape_bounded(tmp.ptr, tmp.sizeof, s);
    pputs(prn, "\\\\ [-8pt]\n");
    pprintf(prn, "%s \\\\[1ex]\n".ptr, tmp.ptr);
    } else if (rtf_format(prn)) {
    pputs(prn, RTF_MULTICOL.ptr);
    pprintf(prn, "\\ql %s".ptr, s);
    pputs(prn, "\\cell\\intbl\\row\n");
    } else if (csv_format(prn)) {
    pprintf(prn, "\n\"%s\"\n".ptr, s);
    }
}

int print_rq_sequence(const(Model) * pmod, const(DatasetRaw) * dset, Prn * prn)
{
    GretlVector * tauvec = cast(GretlVector *) gretl_model_get_data(pmod, "rq_tauvec");
    GretlMatrixRaw * B = cast(GretlMatrixRaw *) gretl_model_get_data(pmod, "rq_sequence");
    double tau, bi, se = NADBL;
    double blo = 0, bhi = 0;
    char[16] test;
    int n, bcols, taulen = 5;
    int namelen = 8;
    int offset;
    int ntau, i, j, k = 0;
    static const(char) *[6] headings = ["tau".ptr, "coefficient".ptr, "lower".ptr, "upper".ptr, "std. error".ptr, "t-ratio".ptr];
    const(char) * head;

    if (tauvec == null || B == null) {
    return GretlError.data;
    }

    ntau = gretl_vector_get_length(tauvec);
    bcols = gretl_matrix_cols(B);

    for (i=2; i<=pmod.list[0]; i++) {
    n = char_len(dset.varname[pmod.list[i]]);
    if (n > namelen) {
        namelen = n;
    }
    }

    for (j=0; j<ntau; j++) {
    sprintf(test.ptr, "%g".ptr, gretl_vector_get(tauvec, j));
    n = cast(int) strlen(test.ptr);
    if (n > taulen) {
        taulen = n;
    }
    }

    /* headings row */

    bufspace(namelen + 4, prn);
    pprintf(prn, " %-*s".ptr, taulen, _(headings[0]));
    for (i=1; i<4; i++) {
    j = (i == 1)? 1 : (bcols == 3)? i : i + 2;
    head = _(headings[j]);
    n = char_len(head);
    offset = (i == 1)? 13 - n : 12 - n;
    bufspace(offset, prn);
    pputs(prn, head);
    pputc(prn, ' ');
    }
    pputc(prn, '\n');

    /* separator row */
    pputs(prn, "  ");
    n = namelen + 2 + taulen + 1 + 13 * 3;
    for (i=0; i<n; i++) {
    pputc(prn, '-');
    }
    pputc(prn, '\n');

    for (i=0; i<pmod.ncoeff; i++) {
    pprintf(prn, "  %-*s  ".ptr, namelen, dset.varname[pmod.list[i+2]]);
    for (j=0; j<ntau; j++) {
        tau = gretl_vector_get(tauvec, j);
        bi  = gretl_matrix_get(B, k, 0);
        if (bcols == 3) {
        blo = gretl_matrix_get(B, k, 1);
        bhi = gretl_matrix_get(B, k, 2);
        } else {
        se = gretl_matrix_get(B, k, 1);
        }
        if (j > 0) {
        bufspace(namelen + 4, prn);
        }
        if (bcols == 3) {
        pprintf(prn, "%.*f  %#12.6g %#12.6g %#12.6g\n".ptr, taulen - 2,
            tau, bi, blo, bhi);
        } else {
        pprintf(prn, "%.*f  %#12.6g %#12.6g %#12.6g\n".ptr, taulen - 2,
            tau, bi, se, bi / se);
        }
        k++;
    }
    if (i < pmod.ncoeff - 1) {
        pputc(prn, '\n');
    }
    }

    return 0;
}

struct Printval {
    char[36] s;
    int lw, rw;
    double x;
}

void get_number_dims(Printval * v, int * lmax, int * rmax)
{
    char[36] tmp;
    const(char) * s = v.s.ptr;
    char * p;
    int i, n;

    while (*s == ' ') s++;
    n = cast(int) strlen(s);
    for (i=n-1; i>0 && s[i] == ' '; i--) {
    n--;
    }

    *tmp.ptr = '\0';
    strncat(tmp.ptr, s, n);

    p = strchr(tmp.ptr, '.');
    if (p == null) {
    p = strchr(tmp.ptr, ',');
    }

    if (p == null) {
    v.lw = char_len(tmp.ptr);
    v.rw = 0;
    } else {
    n = char_len(tmp.ptr);
    *p = '\0';
    v.lw = char_len(tmp.ptr);
    v.rw = n - v.lw;
    }

    if (v.lw > *lmax) *lmax = v.lw;
    if (v.rw > *rmax) *rmax = v.rw;
}

void print_padded_head(const(char) * s, int w, Prn * prn)
{
    int hlen = char_len(s);
    int offset = (w - hlen) / 2;
    int pad = w - hlen - offset;

    bufspace(offset, prn);
    pputs(prn, s);
    bufspace(pad, prn);
}

void print_padded_value(Printval * val, int w, int lmax, int addoff, Prn * prn)
{
    int offset = lmax - val.lw + addoff;
    int pad = w - (lmax + val.rw + addoff);

    bufspace(offset, prn);
    pputs(prn, val.s.ptr);
    bufspace(pad, prn);
}

Printval ** allocate_printvals(int n, int m)
{
    Printval ** vals;
    int i, j;

    vals = cast(typeof(vals)) malloc(n * (*vals).sizeof);

    for (i=0; i<n && vals != null; i++) {
    vals[i] = cast(typeof(vals[i])) malloc(m * (**vals).sizeof);
    if (vals[i] == null) {
        for (j=0; j<i; j++) {
        free(vals[j]);
        }
        free(vals);
        vals = null;
    }
    }

    return vals;
}

int get_ar_data(const(Model) * pmod, double ** pb, double ** pse, int * pk)
{
    double * b, se;
    int i, k, err = 0;

    if (pmod.arinfo == null ||
    pmod.arinfo.arlist == null ||
    pmod.arinfo.rho == null ||
    pmod.arinfo.sderr == null) {
    return GretlError.data;
    }

    k = pmod.arinfo.arlist[0];

    b = cast(typeof(b)) malloc((pmod.ncoeff + k) * (*b).sizeof);
    se = cast(typeof(se)) malloc((pmod.ncoeff + k) * (*se).sizeof);

    if (b == null || se == null) {
    free(b);
    free(se);
    return GretlError.alloc;
    }

    for (i=0; i<pmod.ncoeff; i++) {
    b[i] = pmod.coeff[i];
    se[i] = pmod.sderr[i];
    }

    for (i=0; i<k; i++) {
    b[pmod.ncoeff + i] = pmod.arinfo.rho[i];
    se[pmod.ncoeff + i] = pmod.arinfo.sderr[i];
    }

    *pb = b;
    *pse = se;
    *pk = k;

    return err;
}

int get_arch_data(const(Model) * pmod, double ** pb, double ** pse, int * pk)
{
    double * ab = cast(double *) gretl_model_get_data(pmod, "arch_coeff");
    double * ase = cast(double *) gretl_model_get_data(pmod, "arch_sderr");
    int order = gretl_model_get_int(pmod, "arch_order");
    double * b, se;
    int err = 0;

    if (ab != null && ase != null && order > 0) {
    int i, k = order + 1;

    b = cast(typeof(b)) malloc((pmod.ncoeff + k) * (*b).sizeof);
    se = cast(typeof(se)) malloc((pmod.ncoeff + k) * (*se).sizeof);

    if (b == null || se == null) {
        free(b);
        free(se);
        return GretlError.alloc;
    }

    for (i=0; i<pmod.ncoeff; i++) {
        b[i] = pmod.coeff[i];
        se[i] = pmod.sderr[i];
    }

    for (i=0; i<k; i++) {
        b[pmod.ncoeff + i] = ab[i];
        se[pmod.ncoeff + i] = ase[i];
    }

    *pb = b;
    *pse = se;
    *pk = k;
    } else {
    err = GretlError.data;
    }

    return err;
}

void figure_colsep(int namelen, int ncols, int * w, int * colsep)
{
    int j, n = namelen + ncols * *colsep;

    for (j=0; j<ncols; j++) {
    n += w[j];
    }
    n = 68 - (n + 2);
    if (n > ncols) {
    *colsep += 1;
    }
}

int print_sep_row(int namelen, int ncols, int * w, int colsep, Prn * prn)
{
    int j, n;

    pputc(prn, '\n');
    bufspace(2, prn);
    n = namelen + ncols * colsep;
    for (j=0; j<ncols; j++) {
    n += w[j];
    }
    for (j=0; j<n; j++) {
    pputc(prn, '-');
    }
    pputc(prn, '\n');

    return n;
}

void put_asts(double pv, Prn * prn)
{
    if (pv < 0.01) {
    pputs(prn, " ***");
    } else if (pv < 0.05) {
    pputs(prn, " **");
    } else if (pv < 0.10) {
    pputs(prn, " *");
    }
}

int alt_print_aux_coeffs(const(double) * b, const(double) * se, const(char) ** names, int nc, int df, int ci, Prn * prn)
{
    ModelCoeff mc;
    int i;

    for (i=0; i<nc; i++) {
    if (na(b[i])) {
        return GretlError.nan;
    }
    }

    alt_print_coeff_table_start(null, ci, prn);

    model_coeff_init(&mc);

    /* print row values */
    for (i=0; i<nc; i++) {
    double sei = isnan(se[i]) ? NADBL : se[i];

    mc.b = b[i];
    mc.se = sei;
    if (na(sei) || sei <= 0) {
        mc.tval = mc.pval = NADBL;
    } else {
        mc.tval = b[i] / sei;
        mc.pval = coeff_pval(ci, mc.tval, df);
    }
    if (tex_format(prn)) {
        tex_escape_special(mc.name.ptr, mc.name.sizeof, names[i]);
    } else {
        *mc.name.ptr = '\0';
        strncat(mc.name.ptr, names[i], MC_NAMELEN - 1);
    }
    alt_print_coeff(&mc, prn);
    }

    print_coeff_table_end(null, prn);

    return 0;
}

int plain_print_aux_coeffs(const(double) * b, const(double) * se, const(char) ** names, int nc, int df, int ci, Prn * prn)
{
    Printval ** vals;
    Printval * vij;
    static const(char) *[4] headings = ["coefficient".ptr, "std. error".ptr, "t-ratio".ptr, "p-value".ptr];
    const(char) * head;
    int[4] lmax = 0;
    int[4] rmax = 0;
    int[4] w;
    int[4] addoff = 0;
    int hlen;
    double tval, pval;
    int namelen = 0;
    int colsep = 2;
    int n, d, i, j;
    int err = 0;

    vals = allocate_printvals(nc, 4);
    if (vals == null) {
    return GretlError.alloc;
    }

    if (df == 0 || ASYMPTOTIC_MODEL(cast(GretlCmdIndex) ci)) {
    headings[2] = "z".ptr;
    }

    for (i=0; i<nc; i++) {
    double sei = isnan(se[i]) ? NADBL : se[i];

    if (na(b[i])) {
        err = GretlError.nan;
        goto bailout;
    }
    n = char_len(names[i]);
    if (n > namelen) {
        namelen = n;
    }
    if (na(sei) || sei <= 0) {
        tval = pval = NADBL;
    } else {
        tval = b[i] / sei;
        pval = coeff_pval(ci, tval, df);
    }
    for (j=0; j<4; j++) {
        if (j < 2) {
        /* coeff, standard error */
        d = get_gretl_digits();
        vals[i][j].x = (j == 0)? b[i] : sei;
        } else if (j == 2) {
        /* t-ratio */
        d = 4;
        vals[i][j].x = tval;
        } else {
        /* p-value */
        d = -4;
        vals[i][j].x = pval;
        }
        gretl_sprint_fullwidth_double(vals[i][j].x, d, vals[i][j].s.ptr, prn);
        get_number_dims(&vals[i][j], &lmax[j], &rmax[j]);
    }
    }

    if (namelen < 8) {
    namelen = 8;
    }

    /* figure appropriate column separation */
    for (j=0; j<4; j++) {
    w[j] = lmax[j] + rmax[j];
    head = _(headings[j]);
    hlen = char_len(head);
    if (hlen > w[j]) {
        addoff[j] = (hlen - w[j]) / 2;
        w[j] = hlen;
    }
    }
    figure_colsep(namelen, 4, w.ptr, &colsep);

    /* print headings.ptr */
    bufspace(namelen + 2 + colsep, prn);
    for (j=0; j<4; j++) {
    print_padded_head(_(headings[j]), w[j], prn);
    if (j < 3) {
        bufspace(colsep, prn);
    }
    }

    /* separator row */
    print_sep_row(namelen, 4, w.ptr, colsep, prn);

    /* print row values */
    for (i=0; i<nc; i++) {
    pprintf(prn, "  %-*s".ptr, namelen, names[i]);
    bufspace(colsep, prn);
    for (j=0; j<4; j++) {
        vij = &vals[i][j];
        print_padded_value(vij, w[j], lmax[j], addoff[j], prn);
        if (j < 3) {
        bufspace(colsep, prn);
        } else if (!na(vij.x)) {
        put_asts(vij.x, prn);
        }
    }
    pputc(prn, '\n');
    }

 bailout:

    for (i=0; i<nc; i++) {
    free(vals[i]);
    }
    free(vals);

    return err;
}

int print_coeff_intervals(const(CoeffIntervals) * cf, Prn * prn)
{
    const(double) * b = cf.coeff;
    const(double) * me = cf.maxerr;
    char ** names = cast(char **) cf.names;
    Printval ** vals;
    Printval * vij;
    static const(char) *[3] headings0 = ["coefficient".ptr, "low".ptr, "high".ptr];
    static const(char) *[4] headings1 = ["coefficient".ptr, "std. error".ptr, "low".ptr, "high".ptr];
    const(char) ** headings = headings0.ptr;
    const(char) * head;
    int nfields = (cf.opt & Gretlopt.OPT_E)? 4 : 3;
    int[4] lmax;
    int[4] rmax;
    int[4] w;
    int[4] addoff;
    int hlen, odds = 0;
    double se, lo, hi, tail;
    int se_col = 0;
    int lo_col = 1;
    int hi_col = 2;
    int namelen = 0;
    int colsep = 2;
    int n, d, i, j;
    int err = 0;

    if (cf.opt & Gretlopt.OPT_E) {
    headings = headings1.ptr;
    se_col = 1;
    lo_col++;
    hi_col++;
    }

    if (cf.opt & Gretlopt.OPT_O) {
    headings[0] = "odds ratio".ptr;
    odds = 1;
    }

    vals = allocate_printvals(cf.ncoeff, nfields);
    if (vals == null) {
    return GretlError.alloc;
    }

    for (i=0; i<nfields; i++) {
    lmax[i] = rmax[i] = w[i] = addoff[i] = 0;
    }

    d = get_gretl_digits();
    tail = cf.alpha / 2;

    /* xgettext:no-c-format */
    pprintf(prn, _("%g%% confidence intervals").ptr, 100 * (1 - cf.alpha));
    if (cf.asy) {
    pprintf(prn, "\nz(%g) = %.4f\n\n".ptr, tail, cf.t);
    } else {
    pprintf(prn, "\nt(%d, %g) = %.3f\n\n".ptr, cf.df, tail, cf.t);
    }

    for (i=0; i<cf.ncoeff; i++) {
    if (na(b[i])) {
        err = GretlError.nan;
        goto bailout;
    }
    n = char_len(names[i]);
    if (n > namelen) {
        namelen = n;
    }
    if (na(me[i])) {
        lo = hi = NADBL;
    } else {
        lo = b[i] - me[i];
        hi = b[i] + me[i];
    }
    for (j=0; j<nfields; j++) {
        if (j == 0) {
        /* coeff or odds ratio */
        vals[i][j].x = odds ? exp(b[i]) : b[i];
        } else if (j == se_col) {
        /* standard error */
        se = me[i] / cf.t;
        vals[i][j].x = odds ? exp(b[i]) * se : se;
        } else if (j == lo_col) {
        /* lower limit */
        vals[i][j].x = odds ? exp(lo) : lo;
        } else if (j == hi_col) {
        /* upper limit */
        vals[i][j].x = odds ? exp(hi) : hi;
        }
        gretl_sprint_fullwidth_double(vals[i][j].x, d, vals[i][j].s.ptr, prn);
        get_number_dims(&vals[i][j], &lmax[j], &rmax[j]);
    }
    }

    if (namelen < 8) {
    namelen = 8;
    }

    /* figure appropriate column separation */
    for (j=0; j<nfields; j++) {
    w[j] = lmax[j] + rmax[j];
    head = _(headings[j]);
    hlen = char_len(head);
    if (hlen > w[j]) {
        addoff[j] = (hlen - w[j]) / 2;
        w[j] = hlen;
    }
    }
    figure_colsep(namelen, 3, w.ptr, &colsep);

    /* print headings */
    bufspace(namelen + 2 + colsep, prn);
    for (j=0; j<nfields; j++) {
    print_padded_head(_(headings[j]), w[j], prn);
    if (j < nfields-1) {
        bufspace(colsep, prn);
    }
    }

    /* separator row */
    print_sep_row(namelen, nfields, w.ptr, colsep, prn);

    /* print row values */
    for (i=0; i<cf.ncoeff; i++) {
    pprintf(prn, "  %-*s".ptr, namelen, names[i]);
    bufspace(colsep, prn);
    for (j=0; j<nfields; j++) {
        vij = &vals[i][j];
        print_padded_value(vij, w[j], lmax[j], addoff[j], prn);
        if (j < nfields-1) {
        bufspace(colsep, prn);
        }
    }
    pputc(prn, '\n');
    }

 bailout:

    for (i=0; i<cf.ncoeff; i++) {
    free(vals[i]);
    }
    free(vals);

    return err;
}

/* Called by external functions that want to print a coefficient
   table that does not reside in a MODEL.
*/

int print_coeffs(const(double) * b, const(double) * se, const(char) ** names, int nc, int df, int ci, Prn * prn)
{
    if (plain_format(prn)) {
    return plain_print_aux_coeffs(b, se, names, nc, df, ci, prn);
    } else {
    return alt_print_aux_coeffs(b, se, names, nc, df, ci, prn);
    }
}

int plain_print_mp_coeffs(const(Model) * pmod, const(DatasetRaw) * dset, Prn * prn)
{
    Printval ** vals;
    static const(char) *[2] headings = ["coefficient".ptr, "std. error".ptr];
    const(double) * b = pmod.coeff;
    const(double) * se = pmod.sderr;
    char ** names = null;
    const(char) * head;
    int[2] lmax = 0;
    int[2] rmax = 0;
    int[2] w;
    int[2] addoff = 0;
    int hlen;
    int n, nc = pmod.ncoeff;
    int namelen = 0;
    int colsep = 2;
    int i, j, minus;
    int err = 0;

    vals = allocate_printvals(nc, 2);
    if (vals == null) {
    return GretlError.alloc;
    }

    names = strings_array_new_with_length(nc, 32);
    if (names == null) {
    err = GretlError.alloc;
    goto bailout;
    }

    minus = (gretl_print_has_minus(prn))? MinusType.utf : MinusType.hyphen;

    for (i=0; i<nc; i++) {
    if (na(b[i])) {
        err = GretlError.nan;
        goto bailout;
    }
    gretl_model_get_param_name(pmod, dset, i, names[i]);
    n = char_len(names[i]);
    if (n > namelen) {
        namelen = n;
    }
    for (j=0; j<2; j++) {
        vals[i][j].x = (j==0)? b[i] : se[i];
        print_fifteen(vals[i][j].s.ptr, vals[i][j].x, minus);
        vals[i][j].lw = lmax[j] = 2;
        vals[i][j].rw = rmax[j] = 19;
    }
    }

    if (namelen < 8) {
    namelen = 8;
    }

    /* figure appropriate column separation */

    for (j=0; j<2; j++) {
    w[j] = lmax[j] + rmax[j];
    head = _(headings[j]);
    hlen = char_len(head);
    if (hlen > w[j]) {
        addoff[j] = (hlen - w[j]) / 2;
        w[j] = hlen;
    }
    }

    figure_colsep(namelen, 2, w.ptr, &colsep);

    /* print headings.ptr */

    bufspace(namelen + 2 + colsep, prn);
    for (j=0; j<2; j++) {
    print_padded_head(_(headings[j]), w[j], prn);
    if (j < 3) {
        bufspace(colsep, prn);
    }
    }

    /* separator row */
    print_sep_row(namelen, 2, w.ptr, colsep, prn);

    /* print row values */

    for (i=0; i<nc; i++) {
    pprintf(prn, "  %-*s".ptr, namelen, names[i]);
    bufspace(colsep, prn);
    for (j=0; j<2; j++) {
        print_padded_value(&vals[i][j], w[j], lmax[j], addoff[j], prn);
        if (j == 0) {
        bufspace(colsep, prn);
        }
    }
    pputc(prn, '\n');
    }

 bailout:

    for (i=0; i<nc; i++) {
    free(vals[i]);
    }
    free(vals);

    strings_array_free(names, nc);

    return err;
}

const(char) * get_col_heading(const(char) ** S, int j, int slopes, int intervals)
{
    if (j == 3 && slopes) {
    return _(S[j+1]);
    } else if (intervals && (j == 1 || j == 2)) {
    return _(S[j+4]);
    } else {
    return _(S[j]);
    }
}

void print_count_offset(const(Model) * pmod, const(DatasetRaw) * dset, Printval * val, int namelen, int colsep, int w, int lmax, int addoff, Prn * prn)
{
    int offvar = gretl_model_get_int(pmod, "offset_var");

    if (offvar > 0) {
    char[VNAMELEN + 6] name; /* "log(" + varname + ")" + NUL */
    int n;

    snprintf(name.ptr, name.sizeof, "log(%s)", dset.varname[offvar]);
    n = cast(int) strlen(name.ptr);
    pprintf(prn, "\n  %-*s".ptr, namelen, name.ptr);
    if (n > namelen) {
        colsep -= n - namelen;
    }
    if (colsep > 0) {
        bufspace(colsep, prn);
    }
    sprintf(val.s.ptr, "%.1f".ptr, 1.0);
    val.lw = val.rw = 1;
    print_padded_value(val, w, lmax, addoff, prn);
    pputc(prn, '\n');
    }
}

void print_ar_sum(const(Model) * pmod, Prn * prn)
{
    if (pmod.arinfo != null && pmod.arinfo.arlist[0] > 1) {
    double arsum = 0.0;
    int i;

    for (i=0; i<pmod.arinfo.arlist[0]; i++) {
        arsum += pmod.arinfo.rho[i];
    }
    pprintf(prn, "\n%s = %#g\n".ptr, _("Sum of AR coefficients").ptr, arsum);
    }
}

void mn_logit_coeffsep(char * sep, size_t sep_len, const(Model) * pmod, const(DatasetRaw) * dset, int i)
{
    const(char) * vname = gretl_model_get_depvar_name(pmod, dset);
    int v = pmod.list[1];

    if (is_string_valued(dset, v)) {
        /* the string-value label has no length limit anywhere it's
           stored (series_table_add_string(), gretl_string_table.c),
           so bound this write explicitly rather than trusting it
        */
        const(char) * value = series_get_string_for_value(dset, v, i+1);

    snprintf(sep, sep_len, "%s = %s", vname, value);
    } else {
        const(GretlMatrixRaw) * y = cast(const(GretlMatrixRaw) *) gretl_model_get_data(pmod, "yvals");
        int val = cast(int) ((y != null)? y.val[i] : i);

    snprintf(sep, sep_len, "%s = %d", vname, val);
    }
}

int separator_wanted(int i, int seppos, char ** sepstr, const(Model) * pmod)
{
    int ret = 0;

    if (i == seppos) {
        /* the straightforward criterion */
        ret = 1;
    } else if (pmod.ci == GretlCmdIndex.midasreg && sepstr != null) {
        const(int) * seplist = gretl_model_get_list(pmod, "seplist");
        int j = 0;

        if (seplist != null && (j = in_gretl_list(seplist, i)) > 0) {
            char * mtl = get_midas_term_line(pmod, j-1);

            if (*sepstr != null) {
                /* avoid memory leak */
                free(*sepstr);
            }
            *sepstr = mtl;
            if (*sepstr != null) {
                ret = 1;
            }
        }
    } else if (pmod.ci == GretlCmdIndex.biprobit && i == pmod.ncoeff - 1) {
        /* the last coefficient is biprobit rho */
        ret = 1;
    }

    return ret;
}

int plain_print_coeffs(const(Model) * pmod, const(DatasetRaw) * dset, Prn * prn)
{
    Printval ** vals;
    Printval * vij;
    static const(char) *[7] headings = ["coefficient".ptr, "std. error".ptr, "t-ratio".ptr, "p-value".ptr, "slope".ptr, "lower".ptr, "upper".ptr];
    const(double) * b = pmod.coeff;
    const(double) * se = pmod.sderr;
    double * slopes = null;
    char ** names = null;
    const(char) * head;
    char * sepstr = null;
    double * xb = null;
    double * xse = null;
    int seppos = -1, cblock = 0;
    int coeff_digits;
    int[4] lmax = 0;
    int[4] rmax = 0;
    int[4] w;
    int[4] addoff = 0;
    int hlen;
    double tval, pval = 0.0;
    int n, d, nc = pmod.ncoeff;
    int show_slope, adfnum = -1;
    int intervals = 0;
    int dotlen, namelen = 0;
    int colsep = 2;
    int ncols = 4;
    int i, j, k;
    int err = 0;

    if (model_use_zscore(pmod)) {
        headings[2] = "z".ptr;
    }

    if (pmod.ci == GretlCmdIndex.ar || pmod.ci == GretlCmdIndex.arch) {
        k = 0;
        if (pmod.ci == GretlCmdIndex.ar) {
            err = get_ar_data(pmod, &xb, &xse, &k);
        } else {
            err = get_arch_data(pmod, &xb, &xse, &k);
        }
        if (!err) {
            b = xb;
            se = xse;
            nc += k;
        }
    } else if (pmod.ci == GretlCmdIndex.lad || pmod.ci == GretlCmdIndex.quantreg) {
        GretlMatrixRaw * m = cast(GretlMatrixRaw *) gretl_model_get_data(pmod, "coeff_intervals");

        if (m != null) {
            se = m.val;
            intervals = 1;
            ncols = 3;
        }
    } else if (NONLIST_MODEL(cast(GretlCmdIndex) pmod.ci)) {
        headings[0] = "estimate".ptr;
    }

    if (err) {
        return err;
    }

    nc -= gretl_model_get_int(pmod, "skipdums");

    vals = allocate_printvals(nc, ncols);
    if (vals == null) {
        return GretlError.alloc;
    }

    names = strings_array_new_with_length(nc, 32);
    if (names == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    show_slope = binary_model(pmod) && !(pmod.opt & Gretlopt.OPT_P);
    if (show_slope) {
        slopes = cast(typeof(slopes)) gretl_model_get_data(pmod, "slopes");
    }

    if (pmod.aux == ModelAuxCode.adf || pmod.aux == ModelAuxCode.df) {
        adfnum = gretl_model_get_int(pmod, "dfnum");
    }

    gretl_model_get_coeff_separator(pmod, &sepstr, &seppos);
    if (seppos == -1) {
        if (pmod.ci == GretlCmdIndex.garch && pmod.list[0] > 4) {
            seppos = pmod.list[0] - 4;
        } else if (pmod.ci == GretlCmdIndex.ar || pmod.ci == GretlCmdIndex.arch) {
            seppos = pmod.ncoeff;
        } else if (multinomial_model(pmod)) {
            cblock = gretl_model_get_int(pmod, "cblock");
        }
    }

    coeff_digits = get_gretl_digits();

    for (i=0; i<nc; i++) {
        if (na(b[i])) {
            err = GretlError.nan;
            goto bailout;
        }
        gretl_model_get_param_name(pmod, dset, i, names[i]);
        n = char_len(names[i]);
        if (n > namelen) {
            namelen = n;
        }
        if (na(se[i]) || se[i] <= 0) {
            tval = pval = NADBL;
        } else {
            tval = b[i] / se[i];
            if (slopes != null) {
                /* last column: actually slope at mean */
                pval = (pmod.list[i+2] == 0)? 0 : slopes[i];
            } else if (i == adfnum) {
                /* special Dickey-Fuller p-value */
                pval = gretl_model_get_double(pmod, "dfpval");
            } else if (!intervals) {
                /* regular p-value */
                pval = model_coeff_pval(pmod, tval);
            }
        }
        for (j=0; j<ncols; j++) {
            if (j < 2) {
                /* coeff, standard error or lower c.i. limit */
                d = coeff_digits;
                vals[i][j].x = (j==0)? b[i] : se[i];
            } else if (j == 2) {
                /* t-ratio or upper c.i. limit */
                if (intervals) {
                    d = coeff_digits;
                    vals[i][j].x = se[i + nc];
                } else {
                    d = 4;
                    vals[i][j].x = tval;
                }
            } else {
                /* p-value or slope */
                d = (show_slope)? coeff_digits : -4;
                vals[i][j].x = pval;
            }
            if (show_slope && j == 3 && pmod.list[i+2] == 0) {
                /* don't show 'slope' for constant */
                vals[i][j].s[0] = '\0';
                vals[i][j].lw = vals[i][j].rw = 0;
            } else if (pmod.ci == GretlCmdIndex.duration && j > 1 &&
                       !strcmp(names[i], "sigma")) {
                /* suppress result for H0: sigma = 0 */
                vals[i][j].x = NADBL;
                vals[i][j].s[0] = '\0';
                vals[i][j].lw = vals[i][j].rw = 0;
            } else {
                gretl_sprint_fullwidth_double(vals[i][j].x, d, vals[i][j].s.ptr, prn);
                get_number_dims(&vals[i][j], &lmax[j], &rmax[j]);
            }
        }
    }

    if (namelen < 8) {
        namelen = 8;
    } else if (namelen > SHORTLEN) {
        namelen = SHORTLEN;
    }

    /* figure appropriate column separation */

    for (j=0; j<ncols; j++) {
        w[j] = lmax[j] + rmax[j];
        head = get_col_heading(headings.ptr, j, show_slope, intervals);
        hlen = char_len(head);
        if (hlen > w[j]) {
            addoff[j] = (hlen - w[j]) / 2;
            w[j] = hlen;
        }
    }

    figure_colsep(namelen, ncols, w.ptr, &colsep);

    /* print headings.ptr */

    bufspace(namelen + 2 + colsep, prn);
    for (j=0; j<ncols; j++) {
        head = get_col_heading(headings.ptr, j, show_slope, intervals);
        print_padded_head(head, w[j], prn);
        if (j < ncols - 1) {
            bufspace(colsep, prn);
        }
    }

    /* separator row */
    dotlen = print_sep_row(namelen, ncols, w.ptr, colsep, prn);

    /* biprobit special: name of first dependent variable */
    if (pmod.ci == GretlCmdIndex.biprobit) {
        print_coeff_left_string(gretl_model_get_depvar_name(pmod, dset),
                                prn);
    }

    /* print row values */

    k = 0;
    for (i=0; i<nc; i++) {
        char[SHORTLEN] tmp;
        int namepad;

        if (separator_wanted(i, seppos, &sepstr, pmod)) {
            if (pmod.ci == GretlCmdIndex.biprobit) {
                pputc(prn, '\n');
                if (i == seppos) {
                    print_coeff_left_string(sepstr, prn);
                }
            } else {
                print_coeff_separator(sepstr, dotlen, prn);
            }
            free(sepstr);
            sepstr = null;
        } else if (cblock > 0 && i % cblock == 0) {
            char[32] mnlsep;

            mn_logit_coeffsep(mnlsep.ptr, mnlsep.sizeof, pmod, dset, ++k);
        print_coeff_separator(mnlsep.ptr, 0, prn);
    }
    maybe_trim_varname(tmp.ptr, names[i]);
    namepad = cast(int) (namelen + strlen(tmp.ptr) - g_utf8_strlen(tmp.ptr, -1));
    pprintf(prn, "  %-*s".ptr, namepad, tmp.ptr);
    bufspace(colsep, prn);
    for (j=0; j<ncols; j++) {
        vij = &vals[i][j];
        print_padded_value(vij, w[j], lmax[j], addoff[j], prn);
        if (j < ncols - 1) {
        bufspace(colsep, prn);
        } else if (!show_slope && !na(vij.x)) {
        put_asts(vij.x, prn);
        }
    }
    pputc(prn, '\n');
    }

    if (pmod.ci == GretlCmdIndex.ar) {
    print_ar_sum(pmod, prn);
    } else if (COUNT_MODEL(cast(GretlCmdIndex) pmod.ci)) {
    print_count_offset(pmod, dset, &vals[0][0], namelen,
               colsep, w[0], lmax[0], addoff[0],
               prn);
    }

 bailout:

    for (i=0; i<nc; i++) {
    free(vals[i]);
    }
    free(vals);

    strings_array_free(names, nc);
    free(xb);
    free(xse);

    return err;
}

int plain_print_coefficients(const(Model) * pmod, const(DatasetRaw) * dset, Prn * prn)
{
    if (pmod.ncoeff == 0) {
    return 0;
    } else if (gretl_model_get_data(pmod, "rq_sequence") != null) {
    return print_rq_sequence(pmod, dset, prn);
    } else if (pmod.ci == GretlCmdIndex.mpols) {
    return plain_print_mp_coeffs(pmod, dset, prn);
    } else {
    return plain_print_coeffs(pmod, dset, prn);
    }
}

void alt_print_arch_terms(const(Model) * pmod, Prn * prn)
{
    double * a = cast(double *) gretl_model_get_data(pmod, "arch_coeff");
    double * se = cast(double *) gretl_model_get_data(pmod, "arch_sderr");
    int order = gretl_model_get_int(pmod, "arch_order");

    if (a != null && se != null && order > 0) {
    ModelCoeff mc;
    int i;

    gretl_prn_newline(prn);

    for (i=0; i<=order; i++) {
        model_coeff_init(&mc);
        mc.b = a[i];
        mc.se = se[i];
        mc.tval = a[i] / se[i];
        mc.pval = student_pvalue_2(pmod.nobs - (order + 1), mc.tval);

        if (tex_format(prn)) {
        sprintf(mc.name.ptr, "$\\alpha_%d$".ptr, i);
        } else {
        sprintf(mc.name.ptr, "alpha(%d)".ptr, i);
        }

        alt_print_coeff(&mc, prn);
    }
    }
}

void alt_print_count_offset(const(Model) * pmod, const(DatasetRaw) * dset, Prn * prn)
{
    int offvar = gretl_model_get_int(pmod, "offset_var");

    if (offvar > 0) {
    char[VNAMELEN + 6] name; /* "log(" + varname + ")" + NUL */

    snprintf(name.ptr, name.sizeof, "log(%s)", dset.varname[offvar]);
    if (plain_format(prn)) {
        pprintf(prn, "\n  %-13s         1.0\n".ptr, name.ptr);
    } else if (rtf_format(prn)) {
        pputs(prn, RTF_COEFF_ROW.ptr);
        pprintf(prn, "\\ql %s\\cell\\qc 1.0\\cell".ptr, name.ptr);
        pputs(prn, "\\qc \\cell\\qc \\cell \\qc \\cell \\intbl \\row\n");
    } else if (tex_format(prn)) {
        char[48] tmp;

        tex_escape_bounded(tmp.ptr, tmp.sizeof, name.ptr);
        pprintf(prn, "{\\rm %s} & \\multicolumn{1}{c}{1.0} \\\\\n".ptr, tmp.ptr);
    }
    }
}

int alt_print_coefficients(const(Model) * pmod, const(DatasetRaw) * dset, Prn * prn)
{
    GretlMatrixRaw * intervals = null;
    char * sepstr = null;
    int seppos = -1;
    ModelCoeff mc;
    int adfnum = -1;
    int nc = pmod.ncoeff;
    int cols, gotnan = 0;
    int i, err = 0;

    if (gretl_model_get_data(pmod, "rq_sequence") != null) {
    pputs(prn, "Sorry, not implemented yet!\n");
    return 1;
    }

    cols = alt_print_coeff_table_start(pmod, pmod.ci, prn);

    if (pmod.ci == GretlCmdIndex.panel) {
    nc = pmod.list[0] - 1;
    }

    gretl_model_get_coeff_separator(pmod, &sepstr, &seppos);
    if (seppos == -1 && pmod.ci == GretlCmdIndex.garch && pmod.list[0] > 4) {
    seppos = pmod.list[0] - 4;
    }

    if (pmod.ci == GretlCmdIndex.lad || pmod.ci == GretlCmdIndex.quantreg) {
    intervals = cast(typeof(intervals)) gretl_model_get_data(pmod, "coeff_intervals");
    }

    if (pmod.aux == ModelAuxCode.df || pmod.aux == ModelAuxCode.adf) {
    adfnum = gretl_model_get_int(pmod, "dfnum");
    }

    if (pmod.ci == GretlCmdIndex.biprobit) {
    print_coeff_left_string(gretl_model_get_depvar_name(pmod, dset),
                prn);
    }

    for (i=0; i<nc; i++) {
    err = prepare_model_coeff(pmod, dset, i, adfnum, &mc, prn);
    if (err) gotnan = 1;
    if (separator_wanted(i, seppos, &sepstr, pmod)) {
        if (pmod.ci == GretlCmdIndex.biprobit) {
        print_coeff_left_string(sepstr, prn);
        } else {
        print_coeff_separator(sepstr, cols, prn);
        }
        free(sepstr);
        sepstr = null;
    }
    if (intervals != null) {
        mc.lo = gretl_matrix_get(intervals, i, 0);
        mc.hi = gretl_matrix_get(intervals, i, 1);
        mc.show_pval = 0;
    }
    alt_print_coeff(&mc, prn);
    }

    if (pmod.ci == GretlCmdIndex.ar) {
    alt_print_rho_terms(pmod, prn);
    } else if (pmod.ci == GretlCmdIndex.arch) {
    alt_print_arch_terms(pmod, prn);
    } else if (COUNT_MODEL(cast(GretlCmdIndex) pmod.ci)) {
    alt_print_count_offset(pmod, dset, prn);
    }

    return gotnan;
}

void print_rho(const(ArInfo) * arinfo, int c, int dfd, Prn * prn)
{
    double tval = arinfo.rho[c] / arinfo.sderr[c];
    char[32] ustr;

    if (tex_format(prn)) {
    char[32] s1;
    char[32] s2;
    char[32] s3;
    char[32] s4;

    tex_rl_double(arinfo.rho[c], s1.ptr);
    tex_rl_double(arinfo.sderr[c], s2.ptr);
    tex_rl_float(tval, s3.ptr, 4);
    tex_rl_float(student_pvalue_2(dfd, tval), s4.ptr, 4);

    sprintf(ustr.ptr, "$\\hat{u}_{t-%d}$".ptr, arinfo.arlist[c+1]);

    pprintf(prn, "%s &\n  %s &\n    %s &\n      %s &\n        %s \\\\\n".ptr,
        ustr.ptr, s1.ptr, s2.ptr, s3.ptr, s4.ptr);
    } else if (rtf_format(prn)) {
    char[16] pvalstr;
    double pval;

    pputs(prn, RTF_COEFF_ROW.ptr);
    pprintf(prn, "\\ql u(-%d)\\cell".ptr, arinfo.arlist[c+1]);
    rtf_print_double(arinfo.rho[c], prn);
    rtf_print_double(arinfo.sderr[c], prn);
    pprintf(prn, " \\qc %.4f\\cell".ptr, tval);
    pval = student_pvalue_2(dfd, tval);
    print_pval_str(pval, pvalstr.ptr);
    pprintf(prn, " \\qc %s\\cell".ptr, pvalstr.ptr);
    if (pval < 0.01) {
        pputs(prn, " \\ql ***\\cell");
    } else if (pval < 0.05) {
        pputs(prn, " \\ql **\\cell");
    } else if (pval < 0.10) {
        pputs(prn, " \\ql *\\cell");
    } else {
        pputs(prn, " \\ql \\cell");
    }
    pputs(prn, " \\intbl \\row\n");
    }
}

void alt_print_rho_terms(const(Model) * pmod, Prn * prn)
{
    int i, dfd;

    if (pmod.arinfo == null ||
    pmod.arinfo.arlist == null ||
    pmod.arinfo.rho == null ||
    pmod.arinfo.sderr == null) {
    return;
    }

    if (pmod.arinfo.arlist[0] > 1) {
    dfd = pmod.dfd + (pmod.ncoeff - pmod.arinfo.arlist[0]);
    } else {
    dfd = pmod.dfd;
    }

    for (i=1; i<=pmod.arinfo.arlist[0]; i++) {
    print_rho(pmod.arinfo, i - 1, dfd, prn);
    }
}

const(char) * roots_hdr = "                        Real  Imaginary    Modulus  Frequency".ptr;
const(char) * root_fmt = "%8s%3d%17.4f%11.4f%11.4f%11.4f\n";
const(char) * roots_sep = "  -----------------------------------------------------------";

void print_root(double rx, double ix, double mod, double fr, int i, int hline, Prn * prn)
{
    if (plain_format(prn)) {
    pprintf(prn, root_fmt, _("Root").ptr, i, rx, ix, mod, fr);
    } else if (tex_format(prn)) {
    pprintf(prn, "& %s & %d & $%.4f$ & $%.4f$ & $%.4f$ & $%.4f$ \\\\ ".ptr,
        _("Root").ptr, i, rx, ix, mod, fr);
    if (hline) {
        pputs(prn, "\\hline\n");
    } else {
        pputc(prn, '\n');
    }
    } else if (rtf_format(prn)) {
    pputs(prn, RTF_ROOT_ROW.ptr);
    pprintf(prn, "\\ql \\cell \\ql %s %d \\cell \\qr %.4f\\cell \\qr %.4f\\cell \\qr %.4f\\cell \\qr %.4f\\cell \\intbl \\row\n".ptr,
        _("Root").ptr, i, rx, ix, mod, fr);
    }
}

void root_start(const(char) * tag, Prn * prn)
{
    if (plain_format(prn)) {
    pprintf(prn, "  %s\n".ptr, _(tag));
    } else if (tex_format(prn)) {
    pprintf(prn, "%s \\\\ \n".ptr, _(tag));
    } else if (rtf_format(prn)) {
    pputs(prn, RTF_ROOT_ROW.ptr);
    pprintf(prn, "\\ql %s\\cell\\ql \\cell\\ql \\cell\\ql \\cell\\ql \\cell\\ql\\cell \\intbl \\row\n".ptr, _(tag));
    }
}

void mod_and_freq(Cmplx * c, double * mod, double * frq)
{
    if (c.i != 0) {
    *mod = sqrt(c.r * c.r + c.i * c.i);
    } else {
    *mod = fabs(c.r);
    }
    *frq = atan2(c.i, c.r) / M_2PI;
}

void print_arma_roots(const(Model) * pmod, Prn * prn)
{
    Cmplx * roots = cast(Cmplx *) gretl_model_get_data(pmod, "roots");

    if (roots != null) {
    int p = arma_model_nonseasonal_AR_order(pmod);
    int q = arma_model_nonseasonal_MA_order(pmod);
    int P = gretl_model_get_int(pmod, "arma_P");
    int Q = gretl_model_get_int(pmod, "arma_Q");
    int i, k, hline;
    double mod, fr;

    if (plain_format(prn)) {
        pprintf(prn, "%s\n%s\n".ptr, _(roots_hdr), roots_sep);
    } else if (tex_format(prn)) {
        pputs(prn, "\n\\vspace{1em}\n\n");
        pputs(prn, "\\begin{tabular}{llrrrrr}\n");
        pprintf(prn, "& & & %s & %s & %s & %s \\\\ \\hline\n".ptr,
            _("Real").ptr, _("Imaginary").ptr, _("Modulus").ptr, _("Frequency").ptr);
    } else if (rtf_format(prn)) {
        pputs(prn, ("\n\\par\n{" ~ RTF_ROOT_ROW).ptr);
        pprintf(prn, "\\qr \\cell \\qc \\cell \\qc {\\i %s}\\cell \\qc {\\i %s}\\cell \\qc {\\i %s}\\cell \\qc {\\i %s}\\cell \\intbl \\row\n".ptr,
            _("Real").ptr, _("Imaginary").ptr, _("Modulus").ptr, _("Frequency").ptr);
    }

    if (p > 0) {
        k = 1;
        root_start("AR".ptr, prn);
        for (i=0; i<p; i++) {
        mod_and_freq(&roots[i], &mod, &fr);
        if (i == p - 1 && q == 0 && P == 0 && Q == 0) {
            hline = 1;
        } else {
            hline = 0;
        }
        print_root(roots[i].r, roots[i].i, mod, fr, k++, hline, prn);
        }
    }

    if (P > 0) {
        k = 1;
        root_start("AR (seasonal)".ptr, prn);
        for (i=p; i<p+P; i++) {
        mod_and_freq(&roots[i], &mod, &fr);
        if (i == p + P - 1 && q == 0 && Q == 0) {
            hline = 1;
        } else {
            hline = 0;
        }
        print_root(roots[i].r, roots[i].i, mod, fr, k++, hline, prn);
        }
    }

    if (q > 0) {
        k = 1;
        root_start("MA".ptr, prn);
        for (i=p+P; i<p+P+q; i++) {
        mod_and_freq(&roots[i], &mod, &fr);
        if (i == p + P + q - 1 && Q == 0) {
            hline = 1;
        } else {
            hline = 0;
        }
        print_root(roots[i].r, roots[i].i, mod, fr, k++, hline, prn);
        }
    }

    if (Q > 0) {
        k = 1;
        root_start("MA (seasonal)".ptr, prn);
        for (i=p+P+q; i<p+P+q+Q; i++) {
        mod_and_freq(&roots[i], &mod, &fr);
        if (i == p + P + q + Q - 1) {
            hline = 1;
        } else {
            hline = 0;
        }
        print_root(roots[i].r, roots[i].i, mod, fr, k++, hline, prn);
        }
    }

    if (plain_format(prn)) {
        pprintf(prn, "%s\n\n".ptr, roots_sep);
    } else if (tex_format(prn)) {
        pputs(prn, "\\end{tabular}\n");
    } else if (rtf_format(prn)) {
        pputs(prn, "}\n");
    }
    }
}

void print_heckit_stats(const(Model) * pmod, Prn * prn)
{
    int totobs = gretl_model_get_int(pmod, "totobs");
    int cenobs = totobs - pmod.nobs;
    double cenpc = (100.0 * cenobs) / totobs;

    ensure_vsep(prn);

    if (plain_format(prn)) {
    pprintf(prn, "%s: %d\n".ptr, _("Total observations").ptr, totobs);
    pprintf(prn, "%s: %d (%.1f%%)\n".ptr, _("Censored observations").ptr,
        cenobs, cenpc);
    pputc(prn, '\n');
    } else if (rtf_format(prn)) {
    pprintf(prn, (RTFTAB ~ "%s: %d\n").ptr, _("Total observations").ptr, totobs);
    pprintf(prn, (RTFTAB ~ "%s: %d (%.1f%%)\n").ptr, _("Censored observations").ptr,
        cenobs, cenpc);
    } else if (tex_format(prn)) {
    pprintf(prn, "%s: %d \\\\\n".ptr, _("Total observations").ptr, totobs);
    pprintf(prn, "%s: %d (%.1f\\%%) \\\\\n".ptr, _("Censored observations").ptr,
        cenobs, cenpc);
    }
}

void plain_print_act_pred(const(int) * ap, Prn * prn)
{
    int leftlen;
    int numwidth = 1;
    int i, bign = 0;

    for (i=0; i<4; i++) {
    if (ap[i] > bign) {
        bign = ap[i];
    }
    }

    while (bign /= 10) {
    numwidth++;
    }

    leftlen = char_len(_("Actual").ptr) + 3;

    bufspace(leftlen + 2, prn);
    pputs(prn, _("Predicted").ptr);
    pputc(prn, '\n');
    bufspace(leftlen + 3, prn);
    pprintf(prn, "%*d   %*d\n".ptr, numwidth, 0, numwidth, 1);
    bufspace(2, prn);
    pputs(prn, _("Actual").ptr);
    pprintf(prn, " 0  %*d   %*d\n".ptr, numwidth, ap[0], numwidth, ap[1]);
    bufspace(leftlen, prn);
    pprintf(prn, "1  %*d   %*d\n".ptr, numwidth, ap[2], numwidth, ap[3]);
    pputc(prn, '\n');
}

int limdep_df(const(Model) * pmod)
{
    int df = gretl_model_get_int(pmod, "lr_df");

    if (df == 0) {
    df = pmod.dfn;
    }

    return df;
}

void logit_probit_stats(const(Model) * pmod, Prn * prn)
{
    const(int) * act_pred = null;
    int binary, slopes, correct = 0;
    double pc_correct;
    int df = 0;

    if ((pmod.opt & Gretlopt.OPT_M) || gretl_model_get_int(pmod, "ordered")) {
    /* ordered logit/probit or multinomial logit */
    binary = slopes = 0;
    } else {
    binary = 1;
    slopes = !(pmod.opt & Gretlopt.OPT_P);
    }

    /* for overall likelihood ratio test */
    if (!na(pmod.chisq)) {
    df = limdep_df(pmod);
    }

    if (binary) {
    act_pred = cast(typeof(act_pred)) gretl_model_get_data(pmod, "discrete_act_pred");
    if (act_pred != null) {
        correct = act_pred[0] + act_pred[3];
    }
    } else {
    correct = gretl_model_get_int(pmod, "correct");
    }

    pc_correct = 100 * cast(double) correct / pmod.nobs;

    ensure_vsep(prn);

    if (plain_format(prn)) {
    if (correct > 0) {
        pprintf(prn, "%s = %d (%.1f%%)\n".ptr,
            _("Number of cases 'correctly predicted'").ptr,
            correct, pc_correct);
    }
    if (binary) {
        double fXb = gretl_model_get_double(pmod, "fXb");

        if (!na(fXb)) {
        pprintf(prn, "f(beta'x) %s = %.3f\n".ptr, _("at mean of independent vars").ptr,
            fXb);
        }
    }
    if (df) {
        pprintf(prn, "%s: %s(%d) = %g [%.4f]\n".ptr,
            _("Likelihood ratio test").ptr, _("Chi-square").ptr,
            df, pmod.chisq, chisq_cdf_comp(df, pmod.chisq));
    }
    pputc(prn, '\n');
    if (act_pred != null) {
        plain_print_act_pred(act_pred, prn);
    }
    } else if (rtf_format(prn)) {
    pputc(prn, '\n');
    if (slopes) {
        pprintf(prn, "\\par {\\super *}%s\n".ptr, _("Evaluated at the mean").ptr);
    }
    if (correct > 0) {
        pprintf(prn, "\\par %s = %d (%.1f%%)\n".ptr,
            _("Number of cases 'correctly predicted'").ptr,
            correct, pc_correct);
    }
    if (binary) {
        pprintf(prn, "\\par f(beta'x) %s = %.3f\n".ptr, _("at mean of independent vars").ptr,
            pmod.sdy);
    }
    if (df) {
        pprintf(prn, "\\par %s: %s(%d) = %g [%.4f]\n".ptr,
            _("Likelihood ratio test").ptr, _("Chi-square").ptr,
            df, pmod.chisq, chisq_cdf_comp(df, pmod.chisq));
    }
    pputc(prn, '\n');
    } else if (tex_format(prn)) {
    if (slopes) {
        pprintf(prn, "\\begin{center}\n$^*$%s\n\\end{center}\n".ptr,
            _("Evaluated at the mean").ptr);
    }
    if (correct > 0 || df) {
        pputs(prn, "\\vspace{1em}\n\\begin{raggedright}\n");
        if (correct > 0) {
        pprintf(prn, "%s = %d (%.1f %s)\\\\\n".ptr,
            _("Number of cases `correctly predicted'").ptr,
            correct, pc_correct, _("percent").ptr);
        }
        if (df) {
        pprintf(prn, "%s: $\\chi^2(%d)$ = %.3f [%.4f]\\\\\n".ptr,
            _("Likelihood ratio test").ptr,
            df, pmod.chisq, chisq_cdf_comp(df, pmod.chisq));
        }
        pputs(prn, "\\end{raggedright}\n");
    }
    }
}

int ols_print_anova(const(Model) * pmod, Prn * prn)
{
    double mst, msr, mse, rss;
    int n, c1, c2, c3;

    if (pmod.ci != GretlCmdIndex.ols || !pmod.ifc ||
    na(pmod.ess) || na(pmod.tss)) {
    return GretlError.notimp;
    }

    pprintf(prn, "%s:\n\n".ptr, _("Analysis of Variance").ptr);

    if (pmod.dfn == 0) {
    /* degenerate model: const only */
    rss = 0.0;
    } else {
    rss = pmod.tss - pmod.ess;
    }

    c1 = cast(int) g_utf8_strlen(_("Sum of squares").ptr, -1);
    c2 = cast(int) g_utf8_strlen(_("df").ptr, -1);
    c3 = cast(int) g_utf8_strlen(_("Mean square").ptr, -1);

    c1 = (c1 < 35)? 35 : c1;
    c2 = (c2 > 8)? c2 + 1 : (c2 < 8)? 8 : c2;
    c3 = (c3 > 16)? c3 + 1 : (c3 < 16)? 16 : c3;

    /* header strings are right-aligned */
    n = cast(int) g_utf8_strlen(_("Sum of squares").ptr, -1);
    bufspace(c1 - n, prn);
    pputs(prn, _("Sum of squares").ptr);
    n = cast(int) g_utf8_strlen(_("df").ptr, -1);
    bufspace(c2 + 1 - n, prn);
    pputs(prn, _("df").ptr);
    n = cast(int) g_utf8_strlen(_("Mean square").ptr, -1);
    bufspace(c3 + 1 - n, prn);
    pputs(prn, _("Mean square").ptr);
    pputs(prn, "\n\n");
    c1 = 16;

    /* Mean Square, regression */
    msr = rss / pmod.dfn;
    /* string left-aligned with initial offset of 2 */
    n = cast(int) g_utf8_strlen(_("Regression").ptr, -1);
    bufspace(2, prn);
    pputs(prn, _("Regression").ptr);
    bufspace(16 - n, prn);
    if (pmod.dfn == 0) {
    pprintf(prn, " %*g %*d %*s\n".ptr, c1, rss, c2, pmod.dfn, c3, _("undefined").ptr);
    } else {
    pprintf(prn, " %*g %*d %*g\n".ptr, c1, rss, c2, pmod.dfn, c3, msr);
    }

    /* Mean Square, errors */
    mse = pmod.ess / pmod.dfd;
    /* string left-aligned with initial offset of 2 */
    n = cast(int) g_utf8_strlen(_("Residual").ptr, -1);
    bufspace(2, prn);
    pputs(prn, _("Residual").ptr);
    bufspace(16 - n, prn);
    pprintf(prn, " %*g %*d %*g\n".ptr, c1, pmod.ess, c2, pmod.dfd, c3, mse);

    /* Mean Square, total */
    mst = pmod.tss / (pmod.nobs - 1);
    /* string left-aligned with initial offset of 2 */
    n = cast(int) g_utf8_strlen(_("Total").ptr, -1);
    bufspace(2, prn);
    pputs(prn, _("Total").ptr);
    bufspace(16 - n, prn);
    pprintf(prn, " %*g %*d %*g\n".ptr, c1, pmod.tss, c2, pmod.nobs - 1, c3, mst);

    pprintf(prn, "\n  R^2 = %g / %g = %.6f\n".ptr, rss, pmod.tss, rss / pmod.tss);

    if (pmod.dfn == 0) {
    pprintf(prn, "  F(%d, %d) %s\n\n".ptr, pmod.dfn, pmod.dfd, _("undefined").ptr);
    return 0;
    }

    if (pmod.ess == 0 || rss == 0.0) {
    pprintf(prn, "  F(%d, %d) = %g / %g (%s)\n\n".ptr, pmod.dfn, pmod.dfd,
        msr, mse, _("undefined").ptr);
    } else {
    double F = msr / mse;
    double pv = snedecor_cdf_comp(pmod.dfn, pmod.dfd, F);

    pprintf(prn, "  F(%d, %d) = %g / %g = %g ".ptr,
        pmod.dfn, pmod.dfd, msr, mse, F);
    if (pv < .0001) {
        pprintf(prn, "[%s %.3g]\n\n".ptr, _("p-value").ptr, pv);
    } else if (!na(pv)) {
        pprintf(prn, "[%s %.4f]\n\n".ptr, _("p-value").ptr, pv);
    }
    }

    return 0;
}

/**
 * print_model_from_matrices:
 * @cs: k x 2 matrix containing coefficients and standard errors.
 * @adds: matrix containing an additional p statistics, or NULL.
 * @names: array of strings containing all required names.
 * @df: degrees of freedom, or 0 for asymptotic.
 * @opt: may contain OPT_I fpr printing the extra stats inline.
 * @prn: gretl printer.
 *
 * Prints to @prn the coefficient table and optional additional statistics
 * for a model estimated "by hand". Mainly useful for user-written functions.
 *
 * The number of string in the array @names must be >= k + p, where k is
 * the number of coefficients and p the number of additional statistics
 * given in @adds.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int print_model_from_matrices(const(GretlMatrixRaw) * cs, const(GretlMatrixRaw) * adds, GretlArray * names, int df, Gretlopt opt, Prn * prn)
{
    int k = gretl_matrix_rows(cs);
    int p = gretl_vector_get_length(adds);
    const(double) * b = cs.val;
    const(double) * se = b + k;
    char ** S;
    int ns = 0;

    S = gretl_array_get_strings(names, &ns);
    if (S == null || ns < k + p) {
    return GretlError.nonconf;
    }

    if (plain_format(prn)) {
    /* newline here is useless for TeX and makes RTF choke */
    pputc(prn, '\n');
    } else if (csv_format(prn)) {
    set_csv_delim(prn);
    }

    model_format_start(prn);

    print_coeffs(b, se, cast(const(char) **) S, k, df, GretlCmdIndex.modprint, prn);

    if (p > 0) {
    print_model_stats_table(adds.val, cast(const(char) **) S + k,
                p, opt, prn);
    }

    if (plain_format(prn)) {
    pputc(prn, '\n');
    }

    model_format_end(prn);

    return 0;
}

Gretlopt get_printmodel_opt(const(Model) * pmod, Gretlopt opt)
{
    Gretlopt ret = Gretlopt.OPT_NONE;

    /* Screen out any irrelevant and possible confusing
       options given with an estimation command, leaving
       only options that are intended for the printing
       routine.
    */

    /* first handle cases where we should not print */
    if (opt & Gretlopt.OPT_Q) {
    return Gretlopt.OPT_Q;
    } else if (pmod.ci == GretlCmdIndex.mle && (opt & Gretlopt.OPT_A) && !(opt & Gretlopt.OPT_V)) {
    /* mle with the --auxiliary option */
    return Gretlopt.OPT_Q;
    }

    if (opt & Gretlopt.OPT_O) {
    ret |= Gretlopt.OPT_O; /* show covariance matrix */
    }
    if (opt & Gretlopt.OPT_S) {
    ret |= Gretlopt.OPT_S; /* "simple": reduced output */
    }
    if (pmod.ci == GretlCmdIndex.ols && (opt & Gretlopt.OPT_V)) {
    ret |= Gretlopt.OPT_V; /* anova (OLS only) */
    }

    return ret;
}
