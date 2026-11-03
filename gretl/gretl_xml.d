/* Port of gretl_xml.c */

module gretl.gretl_xml;


import core.stdc.config : c_long;
import core.stdc.ctype : isspace;
import core.stdc.errno : ERANGE, errno;
import core.stdc.float_ : DBL_MAX;
import core.stdc.stdint : uint32_t, uint64_t;
import core.stdc.stdio : FILE, FILENAME_MAX, SEEK_CUR, fclose, fprintf, fputs, fread, fscanf, fseek,
    perror, remove, sprintf, sscanf, stderr;
import core.stdc.stdlib : atof, atoi, calloc, free, malloc, realloc, strtod, strtol;
import core.stdc.string : memcpy, strchr, strcmp, strcpy, strcspn, strlen, strncat, strncmp, strrchr,
    strspn;
version (Posix) import core.sys.posix.sys.stat : stat_t; else import gretl.winport : stat_t;
import etc.c.zlib : gzFile, gzclose, gzprintf, gzputs;

import gretl.dataseth;
import gretl.glib : gchar;
import gretl.gmem : g_free;
import gretl.gstrfuncs : g_strdup_printf;
import gretl.swap_bytesh : reverse_double;
import gretl.subsampleh : RESAMPLED, SUBMASK_SENTINEL;
import gretl.gretl_commandsh;
import gretl.gretl_errors : gretl_errmsg_ensure, gretl_errmsg_set, gretl_errmsg_sprintf;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_list : gretl_list_new, gretl_null_list;
import gretl.gretl_listh : LISTSEP;
import gretl.gretl_matrixh;
import gretl.gretl_prn : Prn, pprintf, pputs;
import gretl.libgretlh;
import gretl.missingh : NADBL, na;
import gretl.strutils : colonize_obs, dot_atof, gretl_charsub, gretl_strdup, gretl_string_split, gretl_strndup, gretl_xml_encode, gretl_xml_encode_to_buf, gretl_xml_validate, has_suffix, strings_array_free, strings_array_new, strrslash, switch_ext_new;
import gretl.gretl_matrix : gretl_matrix_alloc, gretl_matrix_free, gretl_matrix_get_colnames, gretl_matrix_get_rownames,
    gretl_matrix_get_t1, gretl_matrix_get_t2, gretl_matrix_is_dated, gretl_matrix_set_colnames,
    gretl_matrix_set_complex_full, gretl_matrix_set_rownames, gretl_matrix_set_t1, gretl_matrix_set_t2,
    gretl_null_matrix_new;
import gretl.gretl_utils : doubles_array_new, gretl_int_from_string;
import gretl.gretl_errors : gretl_error_clear, gretl_errmsg_is_set, gretl_warnmsg_set, gretl_warnmsg_sprintf;
import gretl.gretl_list : in_gretl_list, gretl_list_copy, gretl_list_delete_at_pos;
import gretl.options : get_compression_option;
import gretl.libset : gretl_in_gui_mode;
import gretl.gretl_prn : gretl_gzip_print_new, gretl_print_destroy, gretl_print_new_with_filename, pputc;
import gretl.dataset : dataset_allocate_obs_markers, dataset_allocate_varnames, dataset_drop_last_variables, destroy_dataset, datainfo_new, is_string_valued, panel_group_names_ok, series_attach_string_table, series_destroy_string_table, series_get_compact_method, series_get_display_name, series_get_label, series_get_lag, series_get_midas_freq, series_get_midas_period, series_get_orig_pd, series_get_parent_name, series_get_string_table, series_get_string_vals, series_get_transform, series_is_coded, series_is_discrete, series_is_midas_anchor, series_recode_strings, series_set_compact_method, series_set_display_name, series_set_flag, series_set_label, series_set_lag, series_set_midas_anchor, series_set_midas_freq, series_set_midas_period, series_set_orig_pd, series_set_parent, series_set_transform, dataset_add_series;
import gretl.dataio : check_varname, dateton, is_gzipped, merge_or_replace_data, get_merge_opts, ntolabel;
import gretl.xmlh;
import gretl.gretl_paths : gretl_build_path, gretl_deltree, gretl_dotdir, gretl_fopen, gretl_gzopen, gretl_make_dotpath, gretl_mkdir, gretl_mktemp, gretl_stat;
import gretl.plugins : get_plugin_function;
import gretl.gretl_zip : gretl_unzip_into;
import gretl.gretl_intl : gretl_pop_c_numeric_locale, gretl_push_c_numeric_locale;
import gretl.genmain : series_index;
import gretl.calendar : get_epoch_day;
import gretl.gretl_string_table : series_table_new;
import gretl.gretl_panel : panel_padding_rows, set_panel_structure_from_vars, switch_panel_orientation;
import gretl.subsample : dataset_is_subsampled;
import gretl.gretl_typemap : gretl_type_from_string;
import gretl.gretl_commands : gretl_command_number, gretl_command_word;

private enum GRETLDATA_VERSION = "1.5";
private enum GRETLDATA_COMPAT = "1.3"; /* status quo as of gretl 2018b */

private alias XUC = const(xmlChar) *;
private alias ProgressFn = extern (C) int function(double, double, int);

private enum G_LITTLE_ENDIAN = 1234;
private enum G_BIG_ENDIAN = 4321;
version (LittleEndian) private enum G_BYTE_ORDER = G_LITTLE_ENDIAN;
else private enum G_BYTE_ORDER = G_BIG_ENDIAN;

/* unported dependencies: per-file declarations from the 2026c headers */
extern (C) {
    int include_gfn(const(char) * fname, Gretlopt opt, Prn * prn);
    int read_session_functions_file(const(char) * fname);
}

private int read_file_fail()
{
    const(xmlError) * err = xmlGetLastError();

    if (err != null) {
        const(char) * msg = err.message;

        fprintf(stderr, "xmlError: domain %d\n", err.domain);
        fprintf(stderr, "message: %s\n", msg != null ? msg : "null");
        if (err.file != null) {
            fprintf(stderr, "file: %s\n", err.file);
            fprintf(stderr, "line: %d\n", err.line);
        }
    } else {
        fprintf(stderr, "xmlGetLastError() gave NULL\n");
    }

    return GretlError.fopen;
}

private xmlDocPtr read_gzipped_xml(const(char) * fname, int options)
{
    import etc.c.zlib : gzclose, gzopen, gzread;

    ubyte[] buf;
    ubyte[65536] chunk;
    int n;

    auto fz = gzopen(fname, "rb");
    if (fz is null) {
        return null;
    }
    while ((n = gzread(fz, chunk.ptr, cast(uint) chunk.length)) > 0) {
        buf ~= chunk[0 .. n];
    }
    gzclose(fz);
    if (buf.length == 0) {
        return null;
    }
    return xmlReadMemory(cast(const(char) *) buf.ptr, cast(int) buf.length, fname, null, options);
}

int gretl_xml_open_doc_root(const(char) * fname,
                            const(char) * rootname,
                            xmlDocPtr * pdoc,
                            xmlNodePtr * pnode)
{
    xmlDocPtr doc;
    xmlNodePtr node = null;
    int options;
    int err = 0;

    *pdoc = null;
    if (pnode != null) {
        *pnode = null;
    }

    options = XML_PARSE_HUGE | XML_PARSE_NOBLANKS;

    doc = xmlReadFile(fname, null, options);
    /* some libxml2 builds (vcpkg, recent Ubuntu) do not decompress gzipped files themselves */
    if (doc == null && is_gzipped(fname)) {
        doc = read_gzipped_xml(fname, options);
    }
    if (doc == null) {
        gretl_errmsg_sprintf(_("xmlReadFile failed on %s").ptr, fname);
        err = read_file_fail();
    }

    if (!err && pnode != null) {
        node = xmlDocGetRootElement(doc);
        if (node == null) {
            gretl_errmsg_sprintf(_("%s: empty document").ptr, fname);
            xmlFreeDoc(doc);
            err = 1;
        }
    }

    if (!err && node != null && rootname != null) {
        if (xmlStrcmp(node.name, cast(XUC) rootname)) {
            gretl_errmsg_sprintf(_("File of the wrong type, root node not %s").ptr,
                                 rootname);
            fprintf(stderr, "Unexpected root node '%s'\n", cast(char *) node.name);
            xmlFreeDoc(doc);
            err = 1;
        }
    }

    if (!err) {
        *pdoc = doc;
        if (pnode != null) {
            *pnode = node;
        }
    }

    return err;
}

xmlNodePtr gretl_xml_get_gfn(const(char) * fname,
                             xmlDocPtr * pdoc,
                             int * err)
{
    xmlNodePtr ret = null;
    xmlNodePtr node = null;

    *err = gretl_xml_open_doc_root(fname, "gretl-functions", pdoc, &node);

    if (!*err) {
        xmlNodePtr cur = node.children;

        while (cur != null) {
            if (!xmlStrcmp(cur.name, cast(XUC) "gretl-function-package")) {
                ret = cur;
                break;
            }
            cur = cur.next;
        }
        if (ret == null) {
            *err = GretlError.data;
        }
    }

    return ret;
}

int gretl_xml_read_buffer(const(char) * buf,
                          const(char) * rootname,
                          xmlDocPtr * pdoc,
                          xmlNodePtr * pnode)
{
    xmlDocPtr doc;
    xmlNodePtr node = null;
    int err = 0;

    *pdoc = null;
    if (pnode != null) {
        *pnode = null;
    }

    doc = xmlReadMemory(buf, cast(int) strlen(buf), "include.xml", null,
                        XML_PARSE_HUGE | XML_PARSE_NOBLANKS);
    if (doc == null) {
        gretl_errmsg_set(_("xmlReadMemory failed").ptr);
        err = 1;
    }

    if (!err && pnode != null) {
        node = xmlDocGetRootElement(doc);
        if (node == null) {
            gretl_errmsg_set(_("xmlReadMemory: empty document").ptr);
            xmlFreeDoc(doc);
            err = 1;
        }
    }

    if (!err && node != null && rootname != null) {
        if (xmlStrcmp(node.name, cast(XUC) rootname)) {
            gretl_errmsg_sprintf(_("File of the wrong type, root node not %s").ptr,
                                 rootname);
            fprintf(stderr, "Unexpected root node '%s'\n", cast(char *) node.name);
            xmlFreeDoc(doc);
            err = 1;
        }
    }

    if (!err) {
        *pdoc = doc;
        if (pnode != null) {
            *pnode = node;
        }
    }

    return err;
}

private const(char) * compact_method_to_string(int method)
{
    if (method == CompactMethod.sum) return "COMPACT_SUM";
    else if (method == CompactMethod.avg) return "COMPACT_AVG";
    else if (method == CompactMethod.sop) return "COMPACT_SOP";
    else if (method == CompactMethod.eop) return "COMPACT_EOP";
    else return "COMPACT_UNSET";
}

private int compact_string_to_int(const(char) * str)
{
    if (!strcmp(str, "COMPACT_SUM")) return CompactMethod.sum;
    else if (!strcmp(str, "COMPACT_AVG")) return CompactMethod.avg;
    else if (!strcmp(str, "COMPACT_SOP")) return CompactMethod.sop;
    else if (!strcmp(str, "COMPACT_EOP")) return CompactMethod.eop;
    else return CompactMethod.unset;
}

/* given a full filename in @src, write to @dest a "simple"
   counterpart without leading path or extension
*/

private char * simple_fname(char * dest, const(char) * src)
{
    const(char) * s;
    char * p;

    if (src == null) {
        strcpy(dest, "gdttmp");
        return dest;
    }

    s = strrslash(cast(char *) src);

    /* take last part of src filename */
    if (s != null) {
        strcpy(dest, s + 1);
    } else {
        strcpy(dest, src);
    }

    /* trash any extension */
    p = strrchr(dest, '.');
    if (p != null && strlen(dest) > 3) {
        *p = '\0';
    }

    return dest;
}

private const(char) * data_structure_string(int s)
{
    switch (s) {
    case TsCodes.time_series:
    case TsCodes.special_time_series:
        return "time-series";
    case TsCodes.stacked_time_series:
        return "stacked-time-series";
    case TsCodes.stacked_cross_section:
        return "stacked-cross-section";
    default:
        return "cross-section";
    }
}

/* Retrieve the index number of the i^{th} series to be stored. This
   function allows for the case of @list = NULL (saving all series)
   as well as storing of a subset of series specified by @list.
*/

private int savenum(const(int) * list, int i)
{
    if (list != null) {
        return list[i];
    } else {
        return i;
    }
}

/**
 * gretl_xml_put_int:
 * @tag: name to give value.
 * @i: value to put (as attribute)
 * @fp: file to which to write.
 *
 * Writes to @fp a string of the form "\%s=\%d".
 */

void gretl_xml_put_int(const(char) * tag, int i, Prn * prn)
{
    pprintf(prn, "%s=\"%d\" ", tag, i);
}

void gretl_xml_put_uint32(const(char) * tag, uint32_t u, Prn * prn)
{
    pprintf(prn, "%s=\"%u\" ", tag, u);
}

void gretl_xml_put_uint64(const(char) * tag, uint64_t u, Prn * prn)
{
    pprintf(prn, "%s=\"%llu\" ", tag, u);
}

/**
 * gretl_xml_put_double:
 * @tag: name to give value.
 * @x: value to put (as attribute)
 * @fp: file to which to write.
 *
 * Writes to @fp a string of the form "\%s=\%.17g" if the value of
 * @x is valid, otherwise "\%s=NA".
 */

void gretl_xml_put_double(const(char) * tag, double x, Prn * prn)
{
    if (na(x)) {
        pprintf(prn, "%s=\"NA\" ", tag);
    } else {
        pprintf(prn, "%s=\"%.17g\" ", tag, x);
    }
}

void gretl_xml_put_double_array(const(char) * tag, const(double) * x, int n,
                                Prn * prn)
{
    int i;

    pprintf(prn, "<%s count=\"%d\">\n", tag, n);
    for (i=0; i<n; i++) {
        if (na(x[i])) {
            pputs(prn, "NA ");
        } else {
            pprintf(prn, "%.17g ", x[i]);
        }
    }
    pprintf(prn, "</%s>\n", tag);
}

/**
 * gretl_xml_put_string:
 * @str: string to put.
 * @fp: file to which to write.
 *
 * Write @str to @fp.  If @str needs to have XML-special
 * characters escaped, this will be done automatically.
 * If @str is NULL, this is considered a no-op.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int gretl_xml_put_string(const(char) * str, Prn * prn)
{
    int err = 0;

    if (str == null) {
        return 0;
    }

    if (gretl_xml_validate(str)) {
        pputs(prn, str);
    } else {
        char * xstr = gretl_xml_encode(str);

        if (xstr != null) {
            pputs(prn, xstr);
            free(xstr);
        } else {
            err = GretlError.alloc;
        }
    }

    return err;
}

void gretl_xml_put_strings_array(const(char) * tag, const(char) ** strs,
                                 int n, Prn * prn)
{
    int i;

    if (n < 0) {
        /* assume NULL termination */
        const(char) ** S = strs;

        n = 0;
        while (*S != null) {
            n++;
            S++;
        }
    }

    pprintf(prn, "<%s count=\"%d\">\n", tag, n);
    for (i=0; i<n; i++) {
        gretl_xml_put_string(strs[i], prn);
        pputs(prn, " ");
    }
    pprintf(prn, "</%s>\n", tag);
}

void gretl_xml_put_strings_array_quoted(const(char) * tag,
                                        const(char) ** strs, int n,
                                        Prn * prn)
{
    int i;

    pprintf(prn, "<%s count=\"%d\">\n", tag, n);
    for (i=0; i<n; i++) {
        pputs(prn, "\"");
        gretl_xml_put_string(strs[i], prn);
        pputs(prn, "\" ");
    }
    pprintf(prn, "</%s>\n", tag);
}

void gretl_xml_put_tagged_string(const(char) * tag, const(char) * str,
                                 Prn * prn)
{
    pprintf(prn, "<%s>", tag);
    gretl_xml_put_string(str, prn);
    pprintf(prn, "</%s>\n", tag);
}

/**
 * gretl_xml_put_tagged_string_plus:
 * @tag: name to give string.
 * @str: string to put.
 * @attrib: name of attribute.
 * @attval: value of attribute.
 * @fp: file to which to write.
 *
 * Similar to gretl_xml_put_tagged_string(), but allows
 * (in fact, requires) an attribute name and value along
 * with the element name and value.
 */

void gretl_xml_put_tagged_string_plus(const(char) * tag,
                                      const(char) * str,
                                      const(char) * attrib,
                                      const(char) * attval,
                                      Prn * prn)
{
    pprintf(prn, "<%s %s=\"", tag, attrib);
    gretl_xml_put_string(attval, prn);
    pputs(prn, "\">");
    gretl_xml_put_string(str, prn);
    pprintf(prn, "</%s>\n", tag);
}

void gretl_list_serialize(const(int) * list, const(char) * name,
                          Prn * prn)
{
    int i;

    if (list == null) {
        return;
    }

    if (name == null) {
        pputs(prn, "<list>\n");
    } else {
        pprintf(prn, "<list name=\"%s\">\n", name);
    }

    for (i=0; i<=list[0]; i++) {
        pprintf(prn, "%d ", list[i]);
    }
    pputs(prn, "</list>\n");
}

void gretl_xml_put_tagged_list(const(char) * tag, const(int) * list,
                               Prn * prn)
{
    int i;

    if (list == null) {
        return;
    }

    pprintf(prn, "<%s>\n", tag);
    for (i=0; i<=list[0]; i++) {
        if (list[i] == LISTSEP) {
            pputs(prn, "; ");
        } else {
            pprintf(prn, "%d ", list[i]);
        }
    }
    pprintf(prn, "</%s>\n", tag);
}

/**
 * gretl_matrix_serialize:
 * @m: matrix to be written.
 * @name: name for matrix, or NULL.
 * @fp: stream to which to print.
 *
 * Returns: 0 on success, non-zero if the size of the
 * serialization exceeds XML_MAX_HUGE_LENGTH.
 */

int gretl_matrix_serialize(const(GretlMatrixRaw) * m,
                           const(char) * name,
                           Prn * prn)
{
    int is_complex = 0;
    const(char) ** S;
    size_t len = 0;
    double x;
    int i, j;
    int err = 0;

    if (m == null) {
        return 0;
    }

    if (name == null) {
        pprintf(prn, "<gretl-matrix rows=\"%d\" cols=\"%d\"",
                m.rows, m.cols);
    } else {
        pprintf(prn, "<gretl-matrix name=\"%s\" rows=\"%d\" cols=\"%d\"",
                name, m.rows, m.cols);
    }

    if (m.is_complex) {
        is_complex = 1;
        pputs(prn, " complex=\"1\"");
    }

    if (gretl_matrix_is_dated(m)) {
        int mt1 = gretl_matrix_get_t1(m);
        int mt2 = gretl_matrix_get_t2(m);

        pprintf(prn, " t1=\"%d\" t2=\"%d\"", mt1, mt2);
    }

    S = gretl_matrix_get_colnames(m);
    if (S != null) {
        pputs(prn, " colnames=\"");
        for (j=0; j<m.cols; j++) {
            pputs(prn, S[j]);
            pputs(prn, (j < m.cols - 1)? " " : "\"");
        }
    }

    S = gretl_matrix_get_rownames(m);
    if (S != null) {
        pputs(prn, " rownames=\"");
        for (j=0; j<m.rows; j++) {
            pputs(prn, S[j]);
            pputs(prn, (j < m.rows - 1)? " " : "\"");
        }
    }

    pputs(prn, ">\n");

    if (is_complex) {
        gretl_matrix_set_complex_full(cast(GretlMatrixRaw *) m, 0);
    }

    for (i=0; i<m.rows; i++) {
        for (j=0; j<m.cols; j++) {
            x = gretl_matrix_get(m, i, j);
            len += pprintf(prn, "%.17g ", x);
        }
        len += pputs(prn, "\n");
    }

    if (is_complex) {
        gretl_matrix_set_complex_full(cast(GretlMatrixRaw *) m, 1);
    }

    pputs(prn, "</gretl-matrix>\n");

    if (len > cast(size_t) 1000000000) {
        if (name != null) {
            gretl_errmsg_sprintf("matrix %s: exceeds maximum XML size".ptr, name);
        }
        err = GretlError.data;
    }

    return err;
}

/**
 * gretl_finite_scalar_serialize:
 * @x: numerical value.
 * @name: name of object.
 * @fp: stream for writing.
 *
 * Note: not suitable for writing non-finite values.
 */

void gretl_finite_scalar_serialize(double x, const(char) * name, Prn * prn)
{
    pprintf(prn, "<scalar name=\"%s\" value=\"%.17g\"/>\n", name, x);
}

void gretl_string_serialize(const(char) * s, const(char) * name, Prn * prn)
{
    pprintf(prn, "<string name=\"%s\" value=\"%s\"/>\n", name, s);
}

/**
 * gretl_xml_get_prop_as_int:
 * @node: XML node pointer.
 * @tag: name by which integer property is known.
 * @i: location to write int value.
 *
 * Returns: 1 if an int is found and read successfully, 0
 * otherwise.
 */

int gretl_xml_get_prop_as_int(xmlNodePtr node, const(char) * tag,
                              int * i)
{
    xmlChar * tmp = xmlGetProp(node, cast(XUC) tag);
    int ret = 0;

    if (tmp != null) {
        ret = sscanf(cast(const(char) *) tmp, "%d", i);
        free(tmp);
    }

    return ret;
}

int gretl_xml_get_prop_as_opt(xmlNodePtr node, const(char) * tag,
                              Gretlopt * o)
{
    xmlChar * tmp = xmlGetProp(node, cast(XUC) tag);
    int ret = 0;

    if (tmp != null) {
        int i;

        ret = sscanf(cast(const(char) *) tmp, "%d", &i);
        *o = cast(Gretlopt) i;
        free(tmp);
    }

    return ret;
}

uint32_t gretl_xml_get_prop_as_uint32(xmlNodePtr node,
                                      const(char) * tag,
                                      uint32_t * u)
{
    xmlChar * tmp = xmlGetProp(node, cast(XUC) tag);
    int ret = 0;

    if (tmp != null) {
        ret = sscanf(cast(const(char) *) tmp, "%u", u);
        free(tmp);
    }

    return ret;
}

uint64_t gretl_xml_get_prop_as_uint64(xmlNodePtr node,
                                      const(char) * tag,
                                      uint64_t * u)
{
    xmlChar * tmp = xmlGetProp(node, cast(XUC) tag);
    int ret = 0;

    if (tmp != null) {
        ret = sscanf(cast(const(char) *) tmp, "%llu", u);
        free(tmp);
    }

    return ret;
}

int gretl_xml_get_prop_as_uchar(xmlNodePtr node, const(char) * tag,
                                ubyte * u)
{
    xmlChar * tmp = xmlGetProp(node, cast(XUC) tag);
    int ret = 0;

    if (tmp != null) {
        *u = cast(ubyte) atoi(cast(const(char) *) tmp);
        free(tmp);
        ret = 1;
    }

    return ret;
}

int gretl_xml_get_prop_as_double(xmlNodePtr node, const(char) * tag,
                                 double * x)
{
    char * p, s = cast(char *) xmlGetProp(node, cast(XUC) tag);
    int ret = 0;

    *x = NADBL;

    if (s != null) {
        p = s;
        p += strspn(p, " \r\n");
        if (strncmp(p, "NA", 2)) {
            *x = atof(p);
        }
        free(s);
        ret = 1;
    }

    return ret;
}

int gretl_xml_get_prop_as_string(xmlNodePtr node, const(char) * tag,
                                 char ** pstr)
{
    xmlChar * tmp = xmlGetProp(node, cast(XUC) tag);
    int ret = 0;

    if (tmp != null) {
        *pstr = cast(char *) tmp;
        ret = 1;
    }

    return ret;
}

GretlType gretl_xml_get_type_property(xmlNodePtr node)
{
    char * tmp = cast(char *) xmlGetProp(node, cast(XUC) "type");
    GretlType type = cast(GretlType) 0;

    if (tmp != null) {
        if (!strncmp(tmp, "gretl-", 6)) {
            /* allow, e.g. "gretl-matrix" for "matrix" */
            type = gretl_type_from_string(tmp + 6);
        } else {
            type = gretl_type_from_string(tmp);
        }
        free(tmp);
    }

    return type;
}

int gretl_xml_get_prop_as_bool(xmlNodePtr node, const(char) * tag)
{
    xmlChar * tmp = xmlGetProp(node, cast(XUC) tag);
    int ret = 0;

    if (tmp != null) {
        if (!strcmp(cast(char *) tmp, "true") ||
            !strcmp(cast(char *) tmp, "1")) {
            ret = 1;
        }
        free(tmp);
    }

    return ret;
}

int gretl_xml_node_get_int(xmlNodePtr node, xmlDocPtr doc, int * i)
{
    xmlChar * tmp;
    int ret = 0;

    tmp = xmlNodeListGetString(doc, node.children, 1);

    if (tmp != null) {
        *i = atoi(cast(const(char) *) tmp);
        free(tmp);
        ret = 1;
    }

    return ret;
}

int gretl_xml_node_get_uint32(xmlNodePtr node, xmlDocPtr doc,
                              uint32_t * u)
{
    xmlChar * tmp;
    int ret = 0;

    tmp = xmlNodeListGetString(doc, node.children, 1);

    if (tmp != null) {
        int n = sscanf(cast(const(char) *) tmp, "%u", u);

        free(tmp);
        ret = (n == 1);
    }

    return ret;
}

int gretl_xml_node_get_uint64(xmlNodePtr node, xmlDocPtr doc,
                              uint64_t * u)
{
    xmlChar * tmp;
    int ret = 0;

    tmp = xmlNodeListGetString(doc, node.children, 1);

    if (tmp != null) {
        int n = sscanf(cast(const(char) *) tmp, "%llu", u);

        free(tmp);
        ret = (n == 1);
    }

    return ret;
}

int gretl_xml_node_get_double(xmlNodePtr node, xmlDocPtr doc,
                              double * x)
{
    char * s, p;
    int ret = 0;

    s = cast(char *) xmlNodeListGetString(doc, node.children, 1);

    if (s != null) {
        p = s;
        p += strspn(p, " \r\n");
        if (!strncmp(p, "NA", 2)) {
            *x = NADBL;
        } else {
            *x = atof(p);
        }
        free(s);
        ret = 1;
    }

    return ret;
}

int gretl_xml_node_get_string(xmlNodePtr node, xmlDocPtr doc,
                              char ** pstr)
{
    xmlChar * tmp;
    int ret = 0;

    tmp = xmlNodeListGetString(doc, node.children, 1);

    if (tmp != null) {
        *pstr = cast(char *) tmp;
        ret = 1;
    }

    return ret;
}

char * gretl_xml_get_string(xmlNodePtr node, xmlDocPtr doc)
{
    return cast(char *) xmlNodeListGetString(doc, node.children, 1);
}

/**
 * gretl_xml_node_get_trimmed_string:
 * @node: XML node pointer.
 * @doc: XML document pointer.
 * @pstr: location to receive string.
 *
 * Reads a string from @node and trims both leading and trailing
 * white space.
 *
 * Returns: 1 if a string is found and read successfully, 0
 * otherwise.
 */

int gretl_xml_node_get_trimmed_string(xmlNodePtr node, xmlDocPtr doc,
                                      char ** pstr)
{
    char * tmp;
    char * s;
    int i, len, ret = 0;

    tmp = cast(char *) xmlNodeListGetString(doc, node.children, 1);

    if (tmp != null) {
        s = tmp;
        s += strspn(s, " \t\n\r");
        len = cast(int) strlen(s);
        for (i=len-1; i>=0; i--) {
            if (s[i] == ' ' || s[i] == '\t' ||
                s[i] == '\r' || s[i] == '\n') {
                len--;
            } else {
                break;
            }
        }
        *pstr = gretl_strndup(s, len);
        if (*pstr != null) {
            ret = 1;
        }
        free(tmp);
    } else {
        /* empty string */
        *pstr = gretl_strdup("");
        if (*pstr != null) {
            ret = 1;
        }
    }

    return ret;
}

int * gretl_xml_get_list(xmlNodePtr node, xmlDocPtr doc, int * err)
{
    xmlChar * tmp;
    const(char) * p;
    int * list = null;
    int i, n;

    tmp = xmlNodeListGetString(doc, node.children, 1);

    if (tmp == null) {
        *err = GretlError.data;
    } else {
        p = cast(const(char) *) tmp;
        p += strspn(p, " \r\n"); /* skip space (get to first value) */
        if (sscanf(p, "%d", &n) != 1) {
            *err = GretlError.data;
        } else if (n == 0) {
            list = gretl_null_list();
            free(tmp);
            return list;
        } else if (n < 0) {
            *err = GretlError.data;
        } else {
            p += strcspn(p, " \r\n"); /* skip non-space (get beyond value) */
            list = gretl_list_new(n);
            if (list == null) {
                *err = GretlError.alloc;
            }
        }

        if (list != null && !*err) {
            for (i=1; i<=n && !*err; i++) {
                p += strspn(p, " \r\n"); /* skip space (get to next value) */
                if (*p == ';') {
                    list[i] = LISTSEP;
                } else if (sscanf(p, "%d", &list[i]) != 1) {
                    *err = GretlError.data;
                }
                p += strcspn(p, " \r\n"); /* skip non-space (get beyond value) */
            }
        }

        free(tmp);
    }

    if (list != null && *err) {
        free(list);
        list = null;
    }

    return list;
}

int gretl_xml_child_get_string(xmlNodePtr node, xmlDocPtr doc,
                               const(char) * name, char ** pstr)
{
    xmlNodePtr cur;
    xmlChar * tmp;
    int ret = 0;

    *pstr = null;

    cur = node.children;

    while (cur != null) {
        if (!xmlStrcmp(cur.name, cast(XUC) name)) {
            tmp = xmlNodeListGetString(doc, cur.children, 1);
            if (tmp != null) {
                *pstr = cast(char *) tmp;
                ret = 1;
            }
            break;
        }
        cur = cur.next;
    }

    return ret;
}

private bool SMALLVAL(double x)
{
    return x > -1e-40 && x < 1e-40;
}

/* Sanity ceiling on the element count read from a <count>/<size> XML
   attribute in a data/session/function-package file. This bounds the
   allocation requested so that a crafted file cannot force a
   multi-gigabyte (or, on ILP32 builds, size_t-overflowing) malloc
   just by declaring an implausible count.
*/
private enum XML_ARRAY_MAXLEN = 100000000;

private void * gretl_xml_get_array(xmlNodePtr node, xmlDocPtr doc,
                                   GretlType type,
                                   int * nelem, int * err)
{
    xmlChar * tmp = xmlGetProp(node, cast(XUC) "count");
    const(char) * test = null;
    int * ivals = null;
    double * xvals = null;
    Cmplx * cvals = null;
    void * ptr = null;
    int nread = 0;
    c_long ln;
    int i, n = 0;

    *nelem = 0;

    if (tmp == null) {
        tmp = xmlGetProp(node, cast(XUC) "size");
    }
    if (tmp == null) {
        fprintf(stderr, "gretl_xml_get_array: didn't find count\n");
        *err = GretlError.data;
        return null;
    }

    errno = 0;
    ln = strtol(cast(const(char) *) tmp, &test, 10);
    if (errno || test == cast(const(char) *) tmp || ln <= 0 ||
        ln > XML_ARRAY_MAXLEN) {
        n = 0;
    } else {
        n = cast(int) ln;
    }

    free(tmp);

    if (n <= 0) {
        return null;
    }

    if (type == GretlType.int_array) {
        ivals = cast(int *) malloc(n * (*ivals).sizeof);
        ptr = ivals;
    } else if (type == GretlType.double_array) {
        xvals = cast(double *) malloc(n * (*xvals).sizeof);
        ptr = xvals;
    } else if (type == GretlType.cmplx_array) {
        cvals = cast(Cmplx *) malloc(n * (*cvals).sizeof);
        ptr = cvals;
    }

    if (ptr == null) {
        *err = GretlError.alloc;
        return null;
    }

    tmp = xmlNodeListGetString(doc, node.children, 1);

    if (tmp == null) {
        *err = GretlError.data;
    } else {
        const(char) * s = cast(const(char) *) tmp;
        const(char) * test2; /* C: inner "test" shadows the outer one */

        errno = 0;

        if (type == GretlType.double_array) {
            double x;

            for (i=0; i<n && !*err && *s; i++) {
                while (isspace(*s)) s++;
                x = strtod(s, &test2);
                if (!strncmp(test2, "NA", 2)) {
                    x = NADBL;
                    s = test2 + 2;
                } else {
                    s = test2;
                    if (*s != '\0' && !isspace(*s)) {
                        *err = GretlError.data;
                    } else if (errno) {
                        perror(null);
                        if (!SMALLVAL(x)) {
                            x = NADBL;
                        }
                        errno = 0;
                    }
                }
                xvals[i] = x;
                nread++;
            }
        } else if (type == GretlType.int_array) {
            c_long kl;

            for (i=0; i<n && !*err && *s; i++) {
                while (isspace(*s)) s++;
                kl = strtol(s, &test2, 10);
                if (errno) {
                    *err = GretlError.data;
                } else if (*test2 != '\0' && !isspace(*test2)) {
                    *err = GretlError.data;
                } else {
                    s = test2;
                    ivals[i] = cast(int) kl;
                    nread++;
                }
            }
        } else if (type == GretlType.cmplx_array) {
            double x;
            int n2 = n * 2;
            int rval = 1;

            for (i=0; i<n2 && !*err && *s; i++) {
                while (isspace(*s)) s++;
                x = strtod(s, &test2);
                if (!strncmp(test2, "NA", 2)) {
                    x = NADBL;
                    test2 += 2;
                } else if (errno) {
                    if (SMALLVAL(x)) {
                        errno = 0;
                    } else {
                        perror(null);
                        *err = GretlError.data;
                    }
                } else if (*test2 != '\0' && !isspace(*test2)) {
                    *err = GretlError.data;
                }
                if (!*err) {
                    s = test2;
                    if (rval) {
                        cvals[nread].r = x;
                        rval = 0;
                    } else {
                        cvals[nread].i = x;
                        rval = 1;
                        nread++;
                    }
                }
            }
        }

        free(tmp);

        if (nread < n) {
            fprintf(stderr, "expected %d items in array, but got %d\n", n, nread);
            *err = GretlError.data;
        }
    }

    if (ptr != null && *err) {
        free(ptr);
        ptr = null;
    }

    if (!*err) {
        *nelem = n;
    }

    return ptr;
}

int * gretl_xml_get_int_array(xmlNodePtr node, xmlDocPtr doc,
                              int * nelem, int * err)
{
    return cast(int *) gretl_xml_get_array(node, doc, GretlType.int_array,
                                           nelem, err);
}

double * gretl_xml_get_double_array(xmlNodePtr node, xmlDocPtr doc,
                                    int * nelem, int * err)
{
    int myerr = 0;

    if (err == null) {
        err = &myerr;
    }

    return cast(double *) gretl_xml_get_array(node, doc, GretlType.double_array,
                                              nelem, err);
}

Cmplx * gretl_xml_get_cmplx_array(xmlNodePtr node, xmlDocPtr doc,
                                  int * nelem, int * err)
{
    return cast(Cmplx *) gretl_xml_get_array(node, doc, GretlType.cmplx_array,
                                             nelem, err);
}

private char * chunk_strdup(const(char) * src, const(char) ** ptr, int * err)
{
    char * targ = null;

    if (*src == '\0') {
        *ptr = src;
    } else {
        const(char) * p;
        int len = 0;

        src += strspn(src, " \n");
        p = src;

        if (*src == '"') {
            p = ++src;
            while (*src && *src != '"') {
                len++;
                src++;
            }
            if (*src == '"') {
                src++;
            }
        } else {
            while (*src && !isspace(*src)) {
                len++;
                src++;
            }
        }

        if (ptr != null) {
            *ptr = src;
        }

        if (len > 0) {
            targ = gretl_strndup(p, len);
            if (targ == null) {
                *err = GretlError.alloc;
            }
        }
    }

    if (targ == null && !*err) {
        *err = GretlError.data;
    }

    return targ;
}

/**
 * gretl_xml_get_strings_array:
 * @node: XML node pointer.
 * @doc: XML document pointer.
 * @nelem: location to receive number of elements in array.
 * @slop: if non-zero, allow the number of strings to fall
 * short of the recorded string count by one.
 * @err: location to receive error code.
 *
 * Returns: allocated array of strings read from @node, or
 * %NULL on failure.
 */

char ** gretl_xml_get_strings_array(xmlNodePtr node, xmlDocPtr doc,
                                    int * nelem, int slop, int * err)
{
    xmlChar * tmp = xmlGetProp(node, cast(XUC) "count");
    char ** S = null;
    const(char) * p;
    int i, n = 0;

    if (tmp == null) {
        *err = GretlError.data;
        return null;
    }

    n = atoi(cast(const(char) *) tmp);
    free(tmp);

    if (n > 0) {
        S = strings_array_new(n);
        if (S == null) {
            *err = GretlError.alloc;
        } else {
            tmp = xmlNodeListGetString(doc, node.children, 1);
            if (tmp == null) {
                fprintf(stderr, "xmlNodeListGetString failed\n");
                *err = GretlError.data;
            } else {
                p = cast(const(char) *) tmp;
                for (i=0; i<n && !*err; i++) {
                    S[i] = chunk_strdup(p, &p, err);
                    if (*err == GretlError.data) {
                        if (i == n - 1 && slop) {
                            /* somebody got off by one somewhere? */
                            *err = 0;
                            n--;
                        } else {
                            /* treat as non-fatal? */
                            fprintf(stderr, "repairing missing string value!\n");
                            S[i] = gretl_strdup("empty string");
                            *err = 0;
                        }
                    }
                }
                free(tmp);
            }
        }
    }

    if (S != null && *err) {
        strings_array_free(S, n);
        S = null;
    }

    if (!*err) {
        *nelem = n;
    }

    return S;
}

int gretl_xml_child_get_strings_array(xmlNodePtr node, xmlDocPtr doc,
                                      const(char) * name, char *** pstrs,
                                      int * nstrs)
{
    xmlNodePtr cur = node.children;
    int ret = 0;

    while (cur != null) {
        if (!xmlStrcmp(cur.name, cast(XUC) name)) {
            int err = 0;

            *pstrs = gretl_xml_get_strings_array(cur, doc, nstrs, 0, &err);
            ret = !err;
            break;
        }
        cur = cur.next;
    }

    return ret;
}

private int get_matrix_values_via_file(GretlMatrixRaw * m, const(char) * s)
{
    gchar * fname;
    FILE * fp;
    int err = 0;

    fname = gretl_make_dotpath("matrix.xml.XXXXXX");
    if (fname == null) {
        return GretlError.alloc;
    }

    fp = gretl_mktemp(fname, "wb");
    if (fp == null) {
        g_free(fname);
        return GretlError.fopen;
    }

    fputs(s, fp);
    fclose(fp);
    fp = gretl_fopen(fname, "r");

    if (fp == null) {
        err = GretlError.fopen;
    } else {
        double x;
        int i, j;

        for (i=0; i<m.rows && !err; i++) {
            for (j=0; j<m.cols && !err; j++) {
                if (fscanf(fp, "%lf", &x) != 1) {
                    err = GretlError.data;
                } else {
                    gretl_matrix_set(m, i, j, x);
                }
            }
        }

        fclose(fp);
    }

    remove(fname);
    g_free(fname);

    return err;
}

private int maybe_add_matrix_labels(GretlMatrixRaw * m,
                                    const(char) * s,
                                    int byrow)
{
    int n, err = 0;

    n = (byrow)? m.rows : m.cols;

    if (s != null && *s != '\0') {
        char ** S;
        int ns;

        S = gretl_string_split(s, &ns, " \n\t");
        if (S == null) {
            err = GretlError.alloc;
        } else if (ns != n) {
            err = GretlError.nonconf;
            strings_array_free(S, ns);
        } else if (byrow) {
            gretl_matrix_set_rownames(m, S);
        } else {
            gretl_matrix_set_colnames(m, S);
        }
    }

    return err;
}

/**
 * gretl_xml_get_matrix:
 * @node: XML node pointer.
 * @doc: XML document pointer.
 * @err: location to receive error code.
 *
 * Returns: allocated gretl matrix read from @node, or %NULL
 * on failure.
 */

GretlMatrixRaw * gretl_xml_get_matrix(xmlNodePtr node,
                                   xmlDocPtr doc,
                                   int * err)
{
    GretlMatrixRaw * m = null;
    char * names = null;
    xmlChar * tmp = null;
    const(char) * p;
    double x;
    int rows = 0, cols = 0;
    int t1 = -1, t2 = -1;
    int is_complex = 0;
    int i, j;

    if (!gretl_xml_get_prop_as_int(node, "rows", &rows) ||
        !gretl_xml_get_prop_as_int(node, "cols", &cols)) {
        *err = GretlError.data;
        return null;
    }

    if (gretl_xml_get_prop_as_bool(node, "complex")) {
        is_complex = 1;
        rows *= 2; /* 8-byte indexing */
    }

    if (rows == 0 || cols == 0) {
        /* allow case of empty matrix */
        if (rows == 0 && cols == 0) {
            m = gretl_null_matrix_new();
        } else {
            m = gretl_matrix_alloc(rows, cols);
        }
        if (m == null) {
            *err = GretlError.alloc;
        }
        return m;
    } else if (rows <= 0 || cols <= 0) {
        *err = GretlError.data;
        return null;
    }

    gretl_xml_get_prop_as_int(node, "t1", &t1);
    gretl_xml_get_prop_as_int(node, "t2", &t2);

    m = gretl_matrix_alloc(rows, cols);
    if (m == null) {
        *err = GretlError.alloc;
        return null;
    }

    names = cast(char *) xmlGetProp(node, cast(XUC) "colnames");
    if (names != null) {
        maybe_add_matrix_labels(m, names, 0);
    }

    names = cast(char *) xmlGetProp(node, cast(XUC) "rownames");
    if (names != null) {
        maybe_add_matrix_labels(m, names, 1);
    }

    tmp = xmlNodeListGetString(doc, node.children, 1);
    if (tmp == null) {
        gretl_matrix_free(m);
        *err = GretlError.data;
        return null;
    }

    p = cast(const(char) *) tmp;

    gretl_push_c_numeric_locale();

    if (rows * cols > 5000) {
        /* it's relatively slow to crawl along the string holding
           many matrix elements using sscanf plus str* functions
        */
        *err = get_matrix_values_via_file(m, p);
    } else {
        p += strspn(p, " \r\n");
        for (i=0; i<rows && !*err; i++) {
            for (j=0; j<cols && !*err; j++) {
                if (sscanf(p, "%lf", &x) != 1) {
                    *err = GretlError.data;
                } else {
                    gretl_matrix_set(m, i, j, x);
                }
                p += strspn(p, " \r\n");
                p += strcspn(p, " \r\n");
            }
        }
    }

    gretl_pop_c_numeric_locale();

    free(tmp);

    if (*err) {
        gretl_matrix_free(m);
        m = null;
    } else {
        if (is_complex) {
            gretl_matrix_set_complex_full(m, 1);
        }
        if (t1 >= 0 && t2 >= t1) {
            gretl_matrix_set_t1(m, t1);
            gretl_matrix_set_t2(m, t2);
        }
    }

    return m;
}

/**
 * gretl_xml_get_submask:
 * @node: XML node pointer.
 * @doc: XML document pointer.
 * @pmask: location to receive allocated mask.
 *
 * In XML a mask is written as an array of numerals, '0' or '1',
 * or 'p' to represent panel padding, each separated by single spaces.
 * What's written to @pmask is an array of coniguous bytes, each with
 * value 0, 1 or 'p', terminated by 0x0f (SUBMASK_SENTINEL, see subsample.h).
 *
 * Returns: 0 on success, non-zero on failure.
 */

int gretl_xml_get_submask(xmlNodePtr node, xmlDocPtr doc, char ** pmask)
{
    char * mask = null;
    int i, len;
    int err = 0;

    if (!gretl_xml_get_prop_as_int(node, "length", &len)) {
        return 1;
    }

    if (len == 0) {
        *pmask = RESAMPLED;
        return 0;
    } else if (len < 0) {
        return 1;
    }

    mask = cast(char *) calloc(len, 1);

    if (mask == null) {
        err = 1;
    } else {
        xmlChar * tmp = xmlNodeListGetString(doc, node.children, 1);

        if (tmp == null) {
            err = 1;
        } else {
            char c;
            char * s = cast(char *) tmp;
            size_t rem = strlen(s);

            for (i=0; i<len; i++) {
                if (rem <= 0) {
                    /* declared "length" exceeds the actual number of
                       tokens present in the text content: bail out
                       rather than reading past the string buffer
                    */
                    fprintf(stderr, "Invalid mask!\n");
                    err = 1;
                    break;
                }
                if (*s == 'p') {
                    c = mask[i] = 'p';
                } else {
                    c = mask[i] = cast(char) atoi(s);
                }
                if (c != 0 && c != 1 && c != 'p' && c != SUBMASK_SENTINEL) {
                    fprintf(stderr, "Invalid mask!\n");
                    err = 1;
                    break;
                }
                s += 2; /* skip to next int */
                rem = (rem > 2)? (rem - 2) : 0;
            }
            if (!err && mask[len-1] != SUBMASK_SENTINEL) {
                fprintf(stderr, "Invalid mask!\n");
                err = 1;
            }
            free(tmp);
        }
        if (err) {
            free(mask);
            mask = null;
        }
    }

    if (!err) {
        *pmask = mask;
    }

    return err;
}

void gretl_xml_header(Prn * prn)
{
    pputs(prn, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
}

/**
 * gretl_matrix_write_as_gdt:
 * @fname: name of file to write.
 * @X: matrix, variable in columns.
 * @varnames: column names.
 * @labels: descriptive labels for the variables, or %NULL.
 *
 * Write out a .gdt data file containing the elements of
 * of the given matrix.
 *
 * Returns: 0 on successful completion, non-zero on error.
 */

int gretl_matrix_write_as_gdt(const(char) * fname,
                              const(GretlMatrixRaw) * X,
                              const(char) ** varnames,
                              const(char) ** labels)
{
    gzFile fz = null;
    char[MAXLEN] datname = 0;
    char * xmlbuf = null;
    ProgressFn show_progress = null;
    double msize = 0;
    int T = X.rows;
    int k = X.cols;
    int in_c_locale = 0;
    int i, t, err = 0;

    fz = gretl_gzopen(fname, "wb");

    if (fz == null) {
        gretl_errmsg_sprintf(_("Couldn't open %s for writing").ptr, fname);
        return 1;
    }

    msize = T * k * double.sizeof;

    if (msize > 100000) {
        fprintf(stderr, "Writing %.0f Kbytes of data\n", msize / 1024);
    } else {
        msize = 0;
    }

    if (msize > 0) {
        show_progress = cast(typeof(show_progress)) get_plugin_function("show_progress");
        if (show_progress == null) {
            msize = 0;
        } else {
            show_progress(0, msize, ProgressFlags.save_init);
        }
    }

    simple_fname(datname.ptr, fname);
    xmlbuf = gretl_xml_encode(datname.ptr);
    if (xmlbuf == null) {
        err = 1;
        goto cleanup;
    }

    gzprintf(fz, ("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
             ~ "<!DOCTYPE gretldata SYSTEM \"gretldata.dtd\">\n\n"
             ~ "<gretldata version=\"%s\" name=\"%s\" frequency=\"1\" "
             ~ "startobs=\"1\" endobs=\"%d\" type=\"cross-section\">\n").ptr,
             GRETLDATA_VERSION.ptr, datname.ptr, T);

    free(xmlbuf);

    gretl_push_c_numeric_locale();
    in_c_locale = 1;

    gzprintf(fz, "<variables count=\"%d\">\n", k);

    for (i=0; i<k; i++) {
        gzprintf(fz, "<variable name=\"%s\"", varnames[i]);
        if (labels != null && labels[i] != null) {
            gzprintf(fz, "\n label=\"%s\"", labels[i]);
        }
        gzputs(fz, "\n/>\n");
    }

    gzputs(fz, "</variables>\n");

    gzprintf(fz, "<observations count=\"%d\" labels=\"false\">\n", T);

    for (t=0; t<T; t++) {
        gzputs(fz, "<obs>");
        for (i=0; i<k; i++) {
            gzprintf(fz, "%.12g ", gretl_matrix_get(X, t, i));
        }
        gzputs(fz, "</obs>\n");
        if (msize > 0 && t && (t % 50 == 0)) {
            show_progress(50, T, ProgressFlags.none);
        }
    }

    gzputs(fz, "</observations>\n</gretldata>\n");

 cleanup:

    if (in_c_locale) {
        gretl_pop_c_numeric_locale();
    }

    if (msize > 0) {
        show_progress(0, T, ProgressFlags.finish);
    }

    gzclose(fz);

    return err;
}

private int string_table_count(const(DatasetRaw) * dset,
                               const(int) * list,
                               int nvars)
{
    int i, v, n = 0;

    for (i=1; i<=nvars; i++) {
        v = savenum(list, i);
        if (is_string_valued(dset, v)) {
            n++;
        }
    }

    return n;
}

private void maybe_print_panel_info(const(DatasetRaw) * dset,
                                    int skip_padding,
                                    Prn * prn)
{
    int names = panel_group_names_ok(dset);
    int pd = dset.panel_pd;
    double sd0 = dset.panel_sd0;
    int times = pd > 0 && sd0 > 0.0;

    if (names || times || skip_padding) {
        pputs(prn, "<panel-info\n");
        if (names) {
            pprintf(prn, " group-names=\"%s\"\n", dset.pangrps);
        }
        if (times) {
            pprintf(prn, " time-frequency=\"%d\"\n", pd);
            pprintf(prn, " time-start=\"%.10g\"\n", sd0);
        }
        if (skip_padding) {
            pputs(prn, " skip-padding=\"1\"\n");
        }
        pputs(prn, "/>\n");
    }
}

private int row_is_padding(const(DatasetRaw) * dset, int t, int vmax)
{
    int i;

    for (i=1; i<vmax; i++) {
        if (!na(dset.Z[i][t])) {
            return 0;
        }
    }

    return 1;
}

private Prn * open_gdt_write_stream(const(char) * fname,
                                    Gretlopt opt)
{
    int err = 0;
    Prn * prn;

    if (opt & Gretlopt.OPT_Z) {
        /* using zlib compression */
        int gzlevel = get_compression_option(GretlCmdIndex.store);

        prn = gretl_gzip_print_new(fname, gzlevel, &err);
    } else {
        /* no compression */
        prn = gretl_print_new_with_filename(fname, &err);
    }

    return prn;
}

private enum BIN_HDRLEN = 24;

private void gdt_swap_endianness(DatasetRaw * dset)
{
    int i, t;

    for (i=1; i<dset.v; i++) {
        for (t=0; t<dset.n; t++) {
            reverse_double(dset.Z[i][t]);
        }
    }
}

private int read_binary_header(FILE * fp, int order)
{
    char[BIN_HDRLEN] hdr = 0;
    uint chk;
    int err = 0;

    chk = cast(uint) fread(hdr.ptr, 1, BIN_HDRLEN, fp);

    if (chk != BIN_HDRLEN) {
        err = GretlError.data;
    } else {
        int bin_order = 0;

        if (strncmp(hdr.ptr, "gretl-bin:", 10)) {
            err = GretlError.data;
        } else if (!strcmp(hdr.ptr + 10, "little-endian")) {
            bin_order = G_LITTLE_ENDIAN;
        } else if (!strcmp(hdr.ptr + 10, "big-endian")) {
            bin_order = G_BIG_ENDIAN;
        } else {
            err = GretlError.data;
        }
        if (!err && bin_order != order) {
            err = GretlError.data;
        }
    }

    if (err) {
        gretl_errmsg_set("Error reading binary data file");
    }

    return err;
}

private void na_convert(double * x, int n)
{
    int i;

    for (i=0; i<n; i++) {
        if (x[i] == DBL_MAX) {
            x[i] = NADBL;
        }
    }
}

private int read_binary_data(const(char) * fname,
                             DatasetRaw * dset,
                             int order,
                             double gdtversion,
                             int fullv,
                             const(int) * vlist)
{
    char * bname;
    FILE * fp;
    int err = 0;

    bname = switch_ext_new(fname, "bin");
    fp = gretl_fopen(bname, "rb");

    if (fp == null) {
        err = GretlError.fopen;
    } else {
        int T = dset.n;
        c_long offset = cast(c_long) (T * double.sizeof);
        size_t got;
        int i, k = 1;

        err = read_binary_header(fp, order);

        for (i=1; i<fullv && !err; i++) {
            if (vlist == null || in_gretl_list(vlist, i)) {
                got = fread(dset.Z[k], double.sizeof, T, fp);
                if (got != T) {
                    err = GretlError.data;
                }
                if (!err && gdtversion < 1.4) {
                    /* we need to convert old-style NAs */
                    na_convert(dset.Z[k], dset.n);
                }
                k++;
            } else {
                fseek(fp, offset, SEEK_CUR);
            }
        }
        fclose(fp);
    }

    free(bname);

    if (!err && order != G_BYTE_ORDER) {
        gdt_swap_endianness(dset);
    }

    return err;
}

/* Here we're testing a series to see if it can be represented in full
   precision (discarding any artifacts) using the format "%.15g". The
   criterion is that each value must have a least one trailing zero
   when printed using "% .14e". We don't necessarily check every
   element in the series; assuming a reasonable degree of homogeneity
   in the data it ought to be enough if 100 members pass the test.
*/

private int p15_OK(const(DatasetRaw) * dset, int v)
{
    const(double) * x = dset.Z[v];
    char[32] s = 0;
    int t, i, n_ok = 0;
    int ret = 1;

    for (t=dset.t1; t<=dset.t2; t++) {
        if (!na(x[t])) {
            sprintf(s.ptr, "% .14e", x[t]);
            for (i=16; i>13; i--) {
                if (s[i] != '0') {
                    break;
                }
            }
            if (i == 16) {
                ret = 0;
                break;
            }
            n_ok++;
        }
        if (n_ok > 100) {
            break;
        }
    }

    return ret;
}

/* apparatus from trimming string-values on gdt save,
   and restoration afterwards
*/

private struct StrvalSaver {
    int * list;          /* list of string-valued series */
    double * X;          /* matrix of numeric codes */
    SeriesTable ** st;   /* array of code -> string mappings */
}

private void strval_saver_destroy(StrvalSaver * ss)
{
    if (ss != null) {
        free(ss.list);
        free(ss.X);
        free(ss.st);
        free(ss);
    }
}

/* @dset: the dataset we're storing.
   @nvars: the number of series to be stored.
   @nsv: the number of string-valued series.
   @list: list of series to be stored, or NULL for all.
*/

private StrvalSaver * strval_saver_setup(DatasetRaw * dset,
                                         int nvars, int nsv,
                                         int * list, int * err)
{
    StrvalSaver * ss;
    int n_changed = 0;
    size_t sz;
    int i, n;

    ss = cast(StrvalSaver *) calloc(1, StrvalSaver.sizeof);
    if (ss == null) {
        *err = GretlError.alloc;
        return null;
    }

    n = sample_size(dset);
    ss.list = gretl_list_new(nsv);
    ss.X = cast(double *) malloc(n * nsv * double.sizeof);
    ss.st = cast(SeriesTable **) calloc(nsv, (SeriesTable *).sizeof);
    sz = n * double.sizeof;

    if (ss.list == null || ss.X == null || ss.st == null) {
        *err = GretlError.alloc;
        strval_saver_destroy(ss);
        ss = null;
    } else {
        double * xj = ss.X;
        int changed = 0;
        int v, j = 0;

        for (i=1; i<=nvars && !*err; i++) {
            v = savenum(list, i);
            if (is_string_valued(dset, v)) {
                /* add @v to the strvars list */
                ss.list[j+1] = v;
                /* and get a pointer to its table */
                ss.st[j] = series_get_string_table(dset, v);
                /* copy its code values into column @j of @X */
                memcpy(xj, dset.Z[v] + dset.t1, sz);
                *err = series_recode_strings(dset, v, Gretlopt.OPT_P, &changed);
                if (changed) {
                    n_changed++;
                } else {
                    ss.st[j] = null;
                }
                j++;
                xj += n;
            }
        }
    }

    if (ss != null && n_changed == 0) {
        /* it turns out we have no need for @ss */
        strval_saver_destroy(ss);
        ss = null;
    }

    return ss;
}

private void strval_saver_restore(DatasetRaw * dset,
                                  StrvalSaver * ss)
{
    double * xi = ss.X;
    int n = sample_size(dset);
    size_t sz = n * double.sizeof;
    int i, v;

    for (i=1; i<=ss.list[0]; i++) {
        v = ss.list[i];
        if (ss.st[i-1] != null) {
            series_destroy_string_table(dset, v);
            series_attach_string_table(dset, v, ss.st[i-1]);
            memcpy(dset.Z[v] + dset.t1, xi, sz);
        }
        xi += n;
    }

    strval_saver_destroy(ss);
}

private int gdt_write_observations(const(DatasetRaw) * dset,
                                   const(int) * storelist,
                                   int skip_padding,
                                   int full_range,
                                   ProgressFn show_progress,
                                   Prn * prn)
{
    char[32] buf = 0;
    char[32] numstr = 0;
    char * p15 = null;
    int gdt_digits = 17;
    int i, v, t, nvars;
    int T, t1, t2;
    int err = 0;

    nvars = storelist != null ? storelist[0] : dset.v - 1;

    if (full_range) {
        t1 = 0;
        t2 = dset.n - 1;
    } else {
        t1 = dset.t1;
        t2 = dset.t2;
    }
    T = t2 - t1 + 1;

    p15 = cast(char *) calloc(nvars, 1);
    if (p15 != null) {
        for (i=0; i<nvars; i++) {
            v = savenum(storelist, i+1);
            p15[i] = cast(char) p15_OK(dset, v);
        }
    }

    pputs(prn, "<observations ");
    pprintf(prn, "count=\"%d\" labels=\"%s\"", T,
            (dset.S != null)? "true".ptr : "false".ptr);
    pputs(prn, ">\n");

    for (t=t1; t<=t2; t++) {
        if (skip_padding && row_is_padding(dset, t, dset.v)) {
            continue;
        }
        pputs(prn, "<obs");
        if (dset.S != null) {
            err = gretl_xml_encode_to_buf(buf.ptr, dset.S[t], buf.sizeof);
            if (!err) {
                pprintf(prn, " label=\"%s\"", buf.ptr);
            }
        }
        pputs(prn, ">");
        for (i=1; i<=nvars; i++) {
            v = savenum(storelist, i);
            if (na(dset.Z[v][t])) {
                strcpy(numstr.ptr, "NA ");
            } else if (p15 == null || p15[i-1] == 0) {
                /* use full default precision if required */
                sprintf(numstr.ptr, "%.*g ", gdt_digits, dset.Z[v][t]);
            } else {
                sprintf(numstr.ptr, "%.15g ", dset.Z[v][t]);
            }
            pputs(prn, numstr.ptr);
        }
        if (skip_padding) {
            int unit = 1 + t / dset.pd;
            int time = t % dset.pd + 1;

            sprintf(numstr.ptr, "%d %d ", unit, time);
            pputs(prn, numstr.ptr);
        }
        pputs(prn, "</obs>\n");

        if (show_progress != null && t && ((t - dset.t1) % 50 == 0)) {
            show_progress(50, sample_size(dset), ProgressFlags.none);
        }
    }

    pputs(prn, "</observations>\n");

    free(p15);

    return 0;
}

private int gdt_write_series_info(const(DatasetRaw) * dset,
                                  const(int) * storelist,
                                  int skip_padding,
                                  Prn * prn)
{
    int nvars = storelist != null ? storelist[0] : dset.v - 1;
    char[1024] xmlbuf = 0;
    int i, v, err;

    if (skip_padding) {
        pprintf(prn, "<variables count=\"%d\">\n", nvars + 2);
    } else {
        pprintf(prn, "<variables count=\"%d\">\n", nvars);
    }

    for (i=1; i<=nvars; i++) {
        const(char) * vstr;
        int vprop, mpd;

        v = savenum(storelist, i);
        gretl_xml_encode_to_buf(xmlbuf.ptr, dset.varname[v], xmlbuf.sizeof);
        pprintf(prn, "<variable name=\"%s\"", xmlbuf.ptr);

        if (dset.varinfo == null) {
            pputs(prn, "\n/>\n");
            continue;
        }

        /* string properties */
        vstr = series_get_label(dset, v);
        if (vstr != null && *vstr != '\0') {
            err = gretl_xml_encode_to_buf(xmlbuf.ptr, vstr, xmlbuf.sizeof);
            if (!err) {
                pprintf(prn, "\n label=\"%s\"", xmlbuf.ptr);
            }
        }
        vstr = series_get_display_name(dset, v);
        if (vstr != null && *vstr != '\0') {
            err = gretl_xml_encode_to_buf(xmlbuf.ptr, vstr, xmlbuf.sizeof);
            if (!err) {
                pprintf(prn, "\n displayname=\"%s\"", xmlbuf.ptr);
            }
        }
        vstr = series_get_parent_name(dset, v);
        if (vstr != null) {
            err = gretl_xml_encode_to_buf(xmlbuf.ptr, vstr, xmlbuf.sizeof);
            if (!err) {
                pprintf(prn, "\n parent=\"%s\"", xmlbuf.ptr);
            }
        }

        /* integer properties */
        vprop = series_get_transform(dset, v);
        if (vprop != 0) {
            const(char) * tr = gretl_command_word(vprop);

            pprintf(prn, "\n transform=\"%s\"", tr);
        }
        vprop = series_get_lag(dset, v);
        if (vprop != 0) {
            pprintf(prn, "\n lag=\"%d\"", vprop);
        }
        vprop = series_get_compact_method(dset, v);
        if (vprop != CompactMethod.unset) {
            const(char) * meth = compact_method_to_string(vprop);

            pprintf(prn, "\n compact-method=\"%s\"", meth);
        }

        /* boolean properties */
        if (series_is_discrete(dset, v)) {
            pputs(prn, "\n discrete=\"true\"");
        }
        if (series_is_coded(dset, v)) {
            pputs(prn, "\n coded=\"true\"");
        }
        if (series_is_midas_anchor(dset, v)) {
            pputs(prn, "\n hf-anchor=\"true\"");
        }
        if ((mpd = series_get_midas_period(dset, v)) > 0) {
            pprintf(prn, "\n midas_period=\"%d\"", mpd);
        }
        if ((mpd = series_get_midas_freq(dset, v)) > 0) {
            pprintf(prn, "\n midas_freq=\"%d\"", mpd);
        }
        if ((mpd = series_get_orig_pd(dset, v)) > 0) {
            pprintf(prn, "\n orig_pd=\"%d\"", mpd);
        }

        pputs(prn, "\n/>\n");
    }

    if (skip_padding) {
        pputs(prn, "<variable name=\"unit__\"\n/>\n");
        pputs(prn, "<variable name=\"time__\"\n/>\n");
    }

    pputs(prn, "</variables>\n");

    return 0;
}

private int gdt_write_string_tables(const(DatasetRaw) * dset,
                                    const(int) * storelist,
                                    int ntabs, Prn * prn)
{
    int nvars = storelist != null ? storelist[0] : dset.v - 1;
    char[64] xmlbuf = 0;
    char * sbuf;
    char ** strs;
    int i, j, v, n_strs;

    pprintf(prn, "<string-tables count=\"%d\">\n", ntabs);

    for (i=1; i<=nvars; i++) {
        v = savenum(storelist, i);
        if (!is_string_valued(dset, v)) {
            continue;
        }
        strs = series_get_string_vals(dset, v, &n_strs, 0);
        gretl_xml_encode_to_buf(xmlbuf.ptr, dset.varname[v], xmlbuf.sizeof);
        pprintf(prn, "<valstrings owner=\"%s\" count=\"%d\">", xmlbuf.ptr, n_strs);
        for (j=0; j<n_strs; j++) {
            sbuf = gretl_xml_encode(strs[j]);
            if (sbuf == null || *sbuf == '\0') {
                fprintf(stderr, "string values for var %d: string %d is empty\n",
                        i, j);
                pputs(prn, "\"empty string\" ");
            } else {
                pprintf(prn, "\"%s\" ", sbuf);
            }
            free(sbuf);
        }
        pputs(prn, "</valstrings>\n");
    }

    pputs(prn, "</string-tables>\n");

    return 0;
}

/* end apparatus for trimming string-values */

private int real_write_gdt(const(char) * fname,
                           Prn * inprn,
                           const(int) * storelist,
                           const(DatasetRaw) * dset,
                           Gretlopt opt,
                           int progress)
{
    Prn * prn = null;
    int T = sample_size(dset);
    int dslength = T;
    char[OBSLEN] startdate = 0;
    char[OBSLEN] enddate = 0;
    char[MAXLEN] datname = 0;
    char[32] freqstr = 0;
    char[1024] xmlbuf = 0;
    ProgressFn show_progress = null;
    double dsize = 0;
    int nvars;
    int skip_padding = 0;
    int full_range = 0;
    int uerr = 0;
    int err = 0;

    /* how many series are we storing? */
    nvars = storelist != null ? storelist[0] : dset.v - 1;
    if (nvars < 0) {
        gretl_errmsg_set("Invalid number of series to save");
        return GretlError.data;
    }

    if (inprn != null) {
        /* for clipboard */
        prn = inprn;
    } else {
        const(char) * path = fname;
        gchar * fullname = null;

        if (!has_suffix(fname, ".gdt")) {
            /* force use of .gdt extension for native XML data */
            path = fullname = g_strdup_printf("%s.gdt", fname);
        }
        prn = open_gdt_write_stream(path, opt);
        g_free(fullname);
    }

    if (prn == null) {
        return GretlError.fopen;
    }

    dsize = T * nvars * double.sizeof;

    if (dsize > 100000) {
        fprintf(stderr, "Writing %.0f Kbytes of data\n", dsize / 1024);
    } else if (progress) {
        /* suppress progress bar for smaller data */
        progress = 0;
    }

    if (progress) {
        show_progress = cast(typeof(show_progress)) get_plugin_function("show_progress");
        if (show_progress == null) {
            progress = 0;
        } else {
            show_progress(0, dsize, ProgressFlags.save_init);
        }
    }

    if (opt & Gretlopt.OPT_F) {
        dslength = dset.n;
        ntolabel(startdate.ptr, 0, dset);
        ntolabel(enddate.ptr, dset.n - 1, dset);
    } else {
        ntolabel(startdate.ptr, dset.t1, dset);
        ntolabel(enddate.ptr, dset.t2, dset);
    }

    simple_fname(datname.ptr, fname);
    uerr = gretl_xml_encode_to_buf(xmlbuf.ptr, datname.ptr, xmlbuf.sizeof);
    if (uerr) {
        strcpy(xmlbuf.ptr, "unknown");
    }

    if (custom_time_series(dset)) {
        sprintf(freqstr.ptr, "special:%d", dset.pd);
    } else {
        sprintf(freqstr.ptr, "%d", dset.pd);
    }

    pprintf(prn, ("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
            ~ "<!DOCTYPE gretldata SYSTEM \"gretldata.dtd\">\n\n"
            ~ "<gretldata version=\"%s\" name=\"%s\" frequency=\"%s\" "
            ~ "n=\"%d\" startobs=\"%s\" endobs=\"%s\" ").ptr, GRETLDATA_VERSION.ptr,
            xmlbuf.ptr, freqstr.ptr, dslength, startdate.ptr, enddate.ptr);

    pprintf(prn, "type=\"%s\"", data_structure_string(dset.structure));

    if (dset.rseed > 0) {
        /* record resampling info */
        pprintf(prn, " rseed=\"%llu\"", dset.rseed);
    }
    if (dset.mapfile != null) {
        /* record map link */
        uerr = gretl_xml_encode_to_buf(xmlbuf.ptr, dset.mapfile, xmlbuf.sizeof);
        if (uerr == 0) {
            pprintf(prn, " mapfile=\"%s\"", xmlbuf.ptr);
        }
    }
    if (opt & Gretlopt.OPT_F) {
        pputs(prn, " mpi-transfer=\"true\"");
    }

    pputs(prn, ">\n");

    if (dataset_is_panel(dset) && !dataset_has_markers(dset) &&
        nvars == dset.v - 1 && dsize > 1024 * 1024 * 10) {
        /* we have more than 10 MB of panel data */
        int padrows = panel_padding_rows(dset);

        if (padrows > 0.4 * dset.n) {
            fprintf(stderr, "skip-padding: dropping %d rows\n", padrows);
            skip_padding = 1;
            T -= padrows;
        }
    }

    /* deal with description, if any */
    if (dset.descrip != null) {
        char * dbuf = gretl_xml_encode(dset.descrip);

        if (dbuf != null) {
            pputs(prn, "<description>");
            pputs(prn, dbuf);
            pputs(prn, "</description>\n");
            free(dbuf);
        }
    }

    if (opt & Gretlopt.OPT_F) {
        pprintf(prn, "<sample t1=\"%d\" t2=\"%d\"/>\n", dset.t1, dset.t2);
        full_range = 1;
    }

    gretl_push_c_numeric_locale();

    /* write listing of series */
    gdt_write_series_info(dset, storelist, skip_padding, prn);

    /* write listing of observations */
    if (nvars > 0) {
        gdt_write_observations(dset, storelist, skip_padding,
                               full_range, show_progress, prn);
    }

    /* maybe write string tables */
    if (dset.varinfo != null) {
        int ntabs = string_table_count(dset, storelist, nvars);

        if (ntabs > 0) {
            gdt_write_string_tables(dset, storelist, ntabs, prn);
        }
    }

    if (dataset_is_panel(dset)) {
        maybe_print_panel_info(dset, skip_padding, prn);
    }

    gretl_pop_c_numeric_locale();

    pputs(prn, "</gretldata>\n");

    if (progress) {
        show_progress(0, dset.t2 - dset.t1 + 1, ProgressFlags.finish);
    }
    if (prn != inprn) {
        gretl_print_destroy(prn);
    }

    return err;
}

private int write_purebin(const(char) * fname, const(int) * list,
                          const(DatasetRaw) * dset, Gretlopt opt)
{
    extern (C) int function(const(char) *, const(int) *, const(DatasetRaw) *,
                            Gretlopt) writer;
    int err = 0;

    writer = cast(typeof(writer)) get_plugin_function("purebin_write_data");

    if (writer == null) {
        err = 1;
    } else {
        err = writer(fname, list, dset, opt);
    }

    return err;
}

private int gdt_store_prepare(const(DatasetRaw) * dset,
                              const(int) * inlist,
                              Gretlopt opt,
                              int ** plist,
                              StrvalSaver ** pss)
{
    int * list = null;
    int nvars = dset.v - 1;
    int ntabs = 0;
    int err = 0;

    /* which series are we saving? */
    if (inlist != null) {
        int cpos = in_gretl_list(inlist, 0);

        list = gretl_list_copy(inlist);
        if (cpos > 0) {
            gretl_list_delete_at_pos(list, cpos);
        }
        nvars = list[0];
        *plist = list;
    }

    /* handling of string-valued series when sub-sampled */
    if (dataset_is_subsampled(dset) && !(opt & Gretlopt.OPT_P) &&
        dset.varinfo != null) {
        ntabs = string_table_count(dset, list, nvars);
        if (ntabs > 0) {
            /* set up for trimming strvals */
            StrvalSaver * ss;

            ss = strval_saver_setup(cast(DatasetRaw *) dset, nvars,
                                    ntabs, list, &err);
            *pss = ss;
        }
    }

    if (err && list != null) {
        free(list);
        *plist = null;
    }

    return err;
}

/**
 * gretl_write_gdt:
 * @fname: name of file to write.
 * @list: list of variables to write (or %NULL to write all).
 * @dset: dataset struct.
 * @opt: if %OPT_Z write gzipped data, else uncompressed; if
 * %OPT_F write full data and include sample range data.
 * @progress: may be 1 when called from gui to display progress
 * bar in case of a large data write; generally should be 0.
 *
 * Write out in xml a data file containing the values of the given set
 * of variables.
 *
 * Returns: 0 on successful completion, non-zero on error.
 */

int gretl_write_gdt(const(char) * fname, const(int) * list,
                    const(DatasetRaw) * dset, Gretlopt opt,
                    int progress)
{
    int gdtb = has_suffix(fname, ".gdtb");
    StrvalSaver * ss = null;
    int * storelist = null;
    int err;

    err = gdt_store_prepare(dset, list, opt, &storelist, &ss);
    if (err) {
        return err;
    }

    if (gdtb) {
        /* binary format for gretl >= 2020b */
        err = write_purebin(fname, storelist, dset, opt);
    } else {
        /* plain gdt file */
        err = real_write_gdt(fname, null, storelist, dset, opt, progress);
    }

    if (ss != null) {
        strval_saver_restore(cast(DatasetRaw *) dset, ss);
    }
    free(storelist);

    return err;
}

/* Called from gui/guiprint.c when for copying a dataset
   to the clipboard.
*/

int gretl_write_gdt_to_prn(Prn * prn, const(int) * list,
                           const(DatasetRaw) * dset)
{
    StrvalSaver * ss = null;
    int * storelist = null;
    int err;

    err = gdt_store_prepare(dset, list, Gretlopt.OPT_NONE, &storelist, &ss);
    if (err) {
        return err;
    }

    err = real_write_gdt(null, prn, storelist, dset, Gretlopt.OPT_NONE, 0);

    if (ss != null) {
        strval_saver_restore(cast(DatasetRaw *) dset, ss);
    }
    free(storelist);

    return err;
}

private void transcribe_string(char * targ, const(char) * src, int maxlen)
{
    *targ = '\0';
    strncat(targ, src, maxlen - 1);
}

private int process_sample(xmlNodePtr node, int * sample,
                           DatasetRaw * dset)
{
    xmlChar * st1 = xmlGetProp(node, cast(XUC) "t1");
    xmlChar * st2 = xmlGetProp(node, cast(XUC) "t2");
    int err = 0;

    if (st1 == null || st2 == null) {
        err = GretlError.data;
    } else {
        int t1 = gretl_int_from_string(cast(const(char) *) st1, &err);
        int t2 = gretl_int_from_string(cast(const(char) *) st2, &err);

        if (!err && (t1 < 0 || t2 >= dset.n || t2 < t1)) {
            err = GretlError.data;
        }
        if (!err) {
            sample[0] = t1;
            sample[1] = t2;
        }
    }

    free(st1);
    free(st2);

    return err;
}

/* Note: if @probe is non-zero, this means that we're really just
   scraping series names from the data file, and so we should not
   start allocating memory for dset->Z based on the number of
   series we find.
*/

private int process_varlist(xmlNodePtr node, DatasetRaw * dset, int probe)
{
    xmlNodePtr cur;
    xmlChar * tmp = xmlGetProp(node, cast(XUC) "count");
    int i, nv = 0, err = 0;

    if (tmp != null) {
        if (sscanf(cast(char *) tmp, "%d", &nv) == 1) {
            dset.v = nv + 1;
        } else {
            gretl_errmsg_set(_("Failed to parse count of variables").ptr);
            err = GretlError.data;
        }
        if (!err && dataset_allocate_varnames(dset)) {
            err = GretlError.alloc;
        }
        if (!err && !probe) {
            dset.Z = doubles_array_new(dset.v, 0);
            if (dset.Z == null) {
                err = GretlError.alloc;
            }
        }
        free(tmp);
    } else {
        gretl_errmsg_set(_("Got no variables").ptr);
        err = GretlError.data;
    }

    if (err) {
        return err;
    } else if (nv == 0) {
        return 0;
    }

    /* now get individual variable info: names and labels */
    cur = node.children;
    while (cur && xmlIsBlankNode(cur)) {
        cur = cur.next;
    }

    if (cur == null) {
        gretl_errmsg_set(_("Got no variables").ptr);
        return GretlError.data;
    }

    i = 1;
    while (cur != null) {
        if (!xmlStrcmp(cur.name, cast(XUC) "variable")) {
            if (i >= dset.v) {
                /* there are more <variable> elements than specified
                   by @count: stop before writing out of bounds
                */
                gretl_errmsg_set(_("Number of variables does not match declaration").ptr);
                err = GretlError.data;
                break;
            }
            tmp = xmlGetProp(cur, cast(XUC) "name");
            if (tmp != null) {
                if (strcmp(cast(char *) tmp, "catch") && strcmp(cast(char *) tmp, "const")) {
                    /* temporary hack: allow 'catch', even though it's
                       officially reserved, on account of Ramanathan
                       data file data6-2.gdt
                    */
                    err = check_varname(cast(const(char) *) tmp);
                }
                if (!err) {
                    transcribe_string(dset.varname[i], cast(char *) tmp, VNAMELEN);
                }
                free(tmp);
            } else {
                gretl_errmsg_sprintf(_("Variable %d has no name").ptr, i);
                err = GretlError.data;
            }
            if (err) {
                return err;
            }
            tmp = xmlGetProp(cur, cast(XUC) "label");
            if (tmp != null) {
                series_set_label(dset, i, cast(char *) tmp);
                free(tmp);
            }
            tmp = xmlGetProp(cur, cast(XUC) "displayname");
            if (tmp != null) {
                series_set_display_name(dset, i, cast(char *) tmp);
                free(tmp);
            }
            tmp = xmlGetProp(cur, cast(XUC) "parent");
            if (tmp != null) {
                series_set_parent(dset, i, cast(char *) tmp);
                free(tmp);
            }
            tmp = xmlGetProp(cur, cast(XUC) "transform");
            if (tmp != null) {
                int ci = gretl_command_number(cast(char *) tmp);

                series_set_transform(dset, i, ci);
                free(tmp);
            }
            tmp = xmlGetProp(cur, cast(XUC) "lag");
            if (tmp != null) {
                series_set_lag(dset, i, atoi(cast(char *) tmp));
                free(tmp);
            }
            tmp = xmlGetProp(cur, cast(XUC) "compact-method");
            if (tmp != null) {
                series_set_compact_method(dset, i, compact_string_to_int(cast(char *) tmp));
                free(tmp);
            }
            tmp = xmlGetProp(cur, cast(XUC) "discrete");
            if (tmp != null) {
                if (!strcmp(cast(char *) tmp, "true")) {
                    series_set_flag(dset, i, VarFlags.discrete);
                }
                free(tmp);
            }
            tmp = xmlGetProp(cur, cast(XUC) "coded");
            if (tmp != null) {
                if (!strcmp(cast(char *) tmp, "true")) {
                    series_set_flag(dset, i, VarFlags.coded);
                }
                free(tmp);
            }
            tmp = xmlGetProp(cur, cast(XUC) "hf-anchor");
            if (tmp != null) {
                if (!strcmp(cast(char *) tmp, "true")) {
                    series_set_midas_anchor(dset, i);
                }
                free(tmp);
            }
            tmp = xmlGetProp(cur, cast(XUC) "midas_period");
            if (tmp != null) {
                int mpd = atoi(cast(char *) tmp);

                if (mpd > 0) {
                    series_set_midas_period(dset, i, mpd);
                }
                free(tmp);
            }
            tmp = xmlGetProp(cur, cast(XUC) "midas_freq");
            if (tmp != null) {
                int mpc = atoi(cast(char *) tmp);

                if (mpc > 0) {
                    series_set_midas_freq(dset, i, mpc);
                }
                free(tmp);
            }
            tmp = xmlGetProp(cur, cast(XUC) "orig_pd");
            if (tmp != null) {
                int opd = atoi(cast(char *) tmp);

                if (opd > 0) {
                    series_set_orig_pd(dset, i, opd);
                }
                free(tmp);
            }
            i++;
        }
        cur = cur.next;
    }

    if (i != dset.v) {
        gretl_errmsg_set(_("Number of variables does not match declaration").ptr);
        err = GretlError.data;
    }

    return err;
}

private int process_varlist_subset(xmlNodePtr node, DatasetRaw * dset,
                                   int * fullv, const(int) * vlist)
{
    xmlNodePtr cur;
    xmlChar * tmp = xmlGetProp(node, cast(XUC) "count");
    int nv = vlist[0];
    int i, k, err = 0;

    if (tmp != null) {
        if (sscanf(cast(char *) tmp, "%d", fullv) != 1) {
            err = GretlError.data;
        }
        free(tmp);
    } else {
        err = GretlError.data;
    }

    if (err) {
        return err;
    }

    *fullv += 1; /* allow for const */

    cur = node.children;
    while (cur && xmlIsBlankNode(cur)) {
        cur = cur.next;
    }

    if (cur == null) {
        gretl_errmsg_set(_("Got no variables").ptr);
        return GretlError.data;
    }

    /* allocate the dataset content */

    dset.v = nv + 1;
    if (!err && dataset_allocate_varnames(dset)) {
        err = GretlError.alloc;
    }
    if (!err) {
        dset.Z = cast(double **) malloc(dset.v * (*dset.Z).sizeof);
        if (dset.Z == null) {
            err = GretlError.alloc;
        }
    }

    if (err) {
        return err;
    }

    /* actually read the info */

    k = i = 1;

    while (cur != null) {
        if (!xmlStrcmp(cur.name, cast(XUC) "variable")) {
            int wanted = 0;

            tmp = xmlGetProp(cur, cast(XUC) "name");
            if (tmp != null) {
                if (in_gretl_list(vlist, i)) {
                    wanted = 1;
                    transcribe_string(dset.varname[k], cast(char *) tmp, VNAMELEN);
                }
                free(tmp);
            } else {
                return GretlError.data;
            }

            if (!wanted) {
                i++;
                cur = cur.next;
                continue;
            }

            tmp = xmlGetProp(cur, cast(XUC) "label");
            if (tmp != null) {
                series_set_label(dset, k, cast(char *) tmp);
                free(tmp);
            }
            tmp = xmlGetProp(cur, cast(XUC) "displayname");
            if (tmp != null) {
                series_set_display_name(dset, k, cast(char *) tmp);
                free(tmp);
            }
            tmp = xmlGetProp(cur, cast(XUC) "parent");
            if (tmp != null) {
                series_set_parent(dset, k, cast(char *) tmp);
                free(tmp);
            }
            tmp = xmlGetProp(cur, cast(XUC) "transform");
            if (tmp != null) {
                int ci = gretl_command_number(cast(char *) tmp);

                series_set_transform(dset, k, ci);
                free(tmp);
            }
            tmp = xmlGetProp(cur, cast(XUC) "lag");
            if (tmp != null) {
                series_set_lag(dset, k, atoi(cast(char *) tmp));
                free(tmp);
            }
            tmp = xmlGetProp(cur, cast(XUC) "compact-method");
            if (tmp != null) {
                series_set_compact_method(dset, k, compact_string_to_int(cast(char *) tmp));
                free(tmp);
            }
            tmp = xmlGetProp(cur, cast(XUC) "discrete");
            if (tmp != null) {
                if (!strcmp(cast(char *) tmp, "true")) {
                    series_set_flag(dset, k, VarFlags.discrete);
                }
                free(tmp);
            }
            tmp = xmlGetProp(cur, cast(XUC) "coded");
            if (tmp != null) {
                if (!strcmp(cast(char *) tmp, "true")) {
                    series_set_flag(dset, k, VarFlags.coded);
                }
                free(tmp);
            }
            tmp = xmlGetProp(cur, cast(XUC) "midas_period");
            if (tmp != null) {
                int mpd = atoi(cast(char *) tmp);

                if (mpd > 0) {
                    series_set_midas_period(dset, k, mpd);
                }
                free(tmp);
            }
            tmp = xmlGetProp(cur, cast(XUC) "midas_freq");
            if (tmp != null) {
                int mpf = atoi(cast(char *) tmp);

                if (mpf > 0) {
                    series_set_midas_freq(dset, k, mpf);
                }
                free(tmp);
            }
            tmp = xmlGetProp(cur, cast(XUC) "hf-anchor");
            if (tmp != null) {
                if (!strcmp(cast(char *) tmp, "true")) {
                    series_set_midas_anchor(dset, k);
                }
                free(tmp);
            }

            i++;
            k++;
        }

        cur = cur.next;
    }

    return err;
}

private void set_underflow_warning(int n, const(char) * fname)
{
    gchar * msg;

    fprintf(stderr, "underflow in %s\n", fname);
    msg = g_strdup_printf(_("Data file contains %d subnormal values").ptr, n);
    gretl_warnmsg_set(msg);
    g_free(msg);
}

/* Read the values for all (or selected) variables at
   observation @t */

private int process_values(DatasetRaw * dset,
                           int t, char * s,
                           int fullv,
                           const(int) * vlist,
                           int * n_uflow)
{
    int try_again = 1;
    char * s0 = s;
    char * alt = null;
    char * test;
    double x;
    int i, k;
    int err = 0;

 reread:

    gretl_error_clear();

    for (i=1, k=1; i<fullv && !err; i++) {
        while (isspace(*s)) s++;
        if (vlist != null && !in_gretl_list(vlist, i)) {
            s += strcspn(s, " \t\r\n");
        } else {
            x = strtod(s, &test);
            if (errno == ERANGE && SMALLVAL(x)) {
                errno = 0; /* underflow, treat as OK? */
                fprintf(stderr, "warning, underflow: %g for series %d (%s) at obs %d\n",
                        x, i, dset.varname[i], t + 1);
                *n_uflow += 1;
                s = test;
            } else if (errno) {
                fprintf(stderr, "%s: %d: bad data\n", __FILE__.ptr, __LINE__);
                perror(null);
                err = GretlError.data;
            } else if (!strncmp(test, "NA", 2)) {
                x = NADBL;
                s = test + 2;
            } else if (*test != '\0' && !isspace(*test)) {
                err = 1;
            } else {
                s = test;
            }
            if (!err && t < dset.n) {
                dset.Z[k++][t] = x;
            }
        }
    }

    s += strspn(s, " \t");
    if (try_again && strchr(s, ',')) {
        /* remedial read in case of decimal commas */
        alt = gretl_strdup(s0);
        gretl_charsub(alt, ',', '.');
        s = alt;
        err = 0;
        try_again = 0;
        goto reread;
    }

    if (err && !gretl_errmsg_is_set()) {
        gretl_errmsg_sprintf(_("Failed to parse data values at obs %d").ptr, t+1);
    }

    free(alt);

    return err;
}

private int read_observations(xmlDocPtr doc, xmlNodePtr node,
                              DatasetRaw * dset, double dsize,
                              int binary, double gdtversion,
                              const(char) * fname)
{
    xmlNodePtr cur;
    xmlChar * tmp;
    int n, i, t;
    ProgressFn show_progress = null;
    int progbar = 0;
    int n_uflow = 0;
    int err = 0;

    tmp = xmlGetProp(node, cast(XUC) "count");
    if (tmp == null) {
        return GretlError.data;
    }

    if (sscanf(cast(char *) tmp, "%d", &n) == 1) {
        if (dset.n > 0 && n != dset.n) {
            gretl_errmsg_set(_("Inconsistent record of number of observations").ptr);
            err = GretlError.data;
        } else {
            dset.n = n;
        }
    } else {
        gretl_errmsg_set(_("Failed to parse number of observations").ptr);
        err = GretlError.data;
    }

    free(tmp);
    if (err) {
        return err;
    }

    if (dsize > 100000 && !binary) {
        show_progress = cast(typeof(show_progress)) get_plugin_function("show_progress");
        if (show_progress != null) {
            progbar = 1;
        }
    }

    tmp = xmlGetProp(node, cast(XUC) "labels");
    if (tmp) {
        if (!strcmp(cast(char *) tmp, "true")) {
            if (dataset_allocate_obs_markers(dset)) {
                return GretlError.alloc;
            }
        }
        free(tmp);
    } else {
        return GretlError.data;
    }

    if (dset.endobs[0] == '\0') {
        sprintf(dset.endobs.ptr, "%d", dset.n);
    }

    dset.t2 = dset.n - 1;

    for (i=0; i<dset.v; i++) {
        dset.Z[i] = cast(double *) malloc(dset.n * double.sizeof);
        if (dset.Z[i] == null) {
            return GretlError.alloc;
        }
    }

    for (t=0; t<dset.n; t++) {
        dset.Z[0][t] = 1.0;
    }

    if (binary) {
        err = read_binary_data(fname, dset, binary, gdtversion,
                               dset.v, null);
        if (!dset.markers) {
            goto bailout;
        }
    }

    /* now get individual obs info: labels and values */
    cur = node.children;
    while (cur && xmlIsBlankNode(cur)) {
        cur = cur.next;
    }

    if (cur == null) {
        gretl_errmsg_set(_("Got no observations\n").ptr);
        return GretlError.data;
    }

    if (progbar) {
        show_progress(0, dsize, ProgressFlags.load_init);
    }

    t = 0;
    while (cur != null) {
        if (!xmlStrcmp(cur.name, cast(XUC) "obs")) {
            if (dset.markers) {
                tmp = xmlGetProp(cur, cast(XUC) "label");
                if (tmp) {
                    transcribe_string(dset.S[t], cast(char *) tmp, OBSLEN);
                    free(tmp);
                } else {
                    gretl_errmsg_sprintf(_("Case marker missing at obs %d").ptr, t+1);
                    return GretlError.data;
                }
            }
            if (!binary) {
                tmp = xmlNodeListGetRawString(doc, cur.children, 1);
                if (tmp) {
                    err = process_values(dset, t, cast(char *) tmp, dset.v, null, &n_uflow);
                    free(tmp);
                } else if (dset.v > 1) {
                    gretl_errmsg_sprintf(_("Values missing at observation %d").ptr, t+1);
                    err = GretlError.data;
                }
            }
            t++;
        }

        if (err) {
            break;
        }

        cur = cur.next;

        if (cur != null && t == dset.n) {
            /* got too many observations */
            t = dset.n + 1;
            break;
        }

        if (progbar && t > 0 && t % 50 == 0) {
            show_progress(50, dset.n, ProgressFlags.none);
        }
    }

 bailout:

    if (progbar) {
        show_progress(0, dset.n, ProgressFlags.finish);
    }

    if (!err && t != dset.n) {
        gretl_errmsg_set(_("Number of observations does not match declaration").ptr);
        err = GretlError.data;
    }

    if (!err && n_uflow > 0) {
        set_underflow_warning(n_uflow, fname);
    }

    return err;
}

private int read_observations_subset(xmlDocPtr doc,
                                     xmlNodePtr node,
                                     DatasetRaw * dset,
                                     int binary,
                                     double gdtversion,
                                     const(char) * fname,
                                     int fullv,
                                     const(int) * vlist,
                                     Gretlopt opt)
{
    xmlNodePtr cur;
    xmlChar * tmp;
    int n, i, t;
    int n_uflow = 0;
    int err = 0;

    tmp = xmlGetProp(node, cast(XUC) "count");
    if (tmp == null) {
        return GretlError.data;
    }

    if (sscanf(cast(char *) tmp, "%d", &n) == 1) {
        if (dset.n > 0 && n != dset.n) {
            gretl_errmsg_set(_("Inconsistent record of number of observations").ptr);
            err = GretlError.data;
        } else {
            dset.n = n;
        }
    } else {
        gretl_errmsg_set(_("Failed to parse number of observations").ptr);
        err = GretlError.data;
    }

    free(tmp);
    if (err) {
        return err;
    }

    if (opt & Gretlopt.OPT_M) {
        tmp = xmlGetProp(node, cast(XUC) "labels");
        if (tmp) {
            if (!strcmp(cast(char *) tmp, "true")) {
                if (dataset_allocate_obs_markers(dset)) {
                    return GretlError.alloc;
                }
            }
            free(tmp);
        } else {
            return GretlError.data;
        }
    }

    if (dset.endobs[0] == '\0') {
        sprintf(dset.endobs.ptr, "%d", dset.n);
    }

    dset.t2 = dset.n - 1;

    for (i=0; i<dset.v; i++) {
        dset.Z[i] = cast(double *) malloc(dset.n * double.sizeof);
        if (dset.Z[i] == null) {
            return GretlError.alloc;
        }
    }

    for (t=0; t<dset.n; t++) {
        dset.Z[0][t] = 1.0;
    }

    if (binary) {
        err = read_binary_data(fname, dset, binary, gdtversion,
                               fullv, vlist);
        if (!dset.markers) {
            goto bailout;
        }
    }

    /* now get individual obs info: labels and values */
    cur = node.children;
    while (cur && xmlIsBlankNode(cur)) {
        cur = cur.next;
    }

    if (cur == null) {
        gretl_errmsg_set(_("Got no observations\n").ptr);
        return GretlError.data;
    }

    t = 0;
    while (cur != null) {
        if (!xmlStrcmp(cur.name, cast(XUC) "obs")) {
            if (dset.markers) {
                tmp = xmlGetProp(cur, cast(XUC) "label");
                if (tmp) {
                    transcribe_string(dset.S[t], cast(char *) tmp, OBSLEN);
                    free(tmp);
                } else {
                    gretl_errmsg_sprintf(_("Case marker missing at obs %d").ptr, t+1);
                    return GretlError.data;
                }
            }
            if (!binary) {
                tmp = xmlNodeListGetRawString(doc, cur.children, 1);
                if (tmp) {
                    err = process_values(dset, t, cast(char *) tmp, fullv, vlist, &n_uflow);
                    free(tmp);
                } else if (dset.v > 1) {
                    gretl_errmsg_sprintf(_("Values missing at observation %d").ptr, t+1);
                    err = GretlError.data;
                }
            }
            t++;
        }

        if (err) {
            break;
        }

        cur = cur.next;

        if (cur != null && t == dset.n) {
            /* got too many observations */
            t = dset.n + 1;
            goto bailout;
        }
    }

 bailout:

    if (!err && t != dset.n) {
        gretl_errmsg_set(_("Number of observations does not match declaration").ptr);
        err = GretlError.data;
    }

    if (!err && n_uflow > 0) {
        set_underflow_warning(n_uflow, fname);
    }

    return err;
}

private int owner_id(const(DatasetRaw) * dset, const(char) * s)
{
    int i;

    for (i=1; i<dset.v; i++) {
        if (!strcmp(s, dset.varname[i])) {
            return i;
        }
    }

    return -1;
}

private int process_string_tables(xmlDocPtr doc,
                                  xmlNodePtr node,
                                  DatasetRaw * dset,
                                  int subset)
{
    xmlNodePtr cur = null;
    xmlChar * tmp;
    int ntabs = 0;
    int err = 0;

    tmp = xmlGetProp(node, cast(XUC) "count");

    if (tmp == null) {
        err = GretlError.data;
    } else {
        if (sscanf(cast(char *) tmp, "%d", &ntabs) != 1) {
            err = GretlError.data;
        }
        free(tmp);
    }

    if (!err) {
        cur = node.children;
        while (cur && xmlIsBlankNode(cur)) {
            cur = cur.next;
        }
        if (cur == null) {
            err = GretlError.data;
        }
    }

    if (err) {
        fprintf(stderr, "process_string_tables: returning err = %d\n", err);
        return err;
    }

    while (cur != null && !err) {
        if (!xmlStrcmp(cur.name, cast(XUC) "valstrings")) {
            xmlChar * owner = xmlGetProp(cur, cast(XUC) "owner");
            SeriesTable * st;
            char ** strs = null;
            int v = 0, n_strs = 0;

            if (owner == null) {
                err = GretlError.data;
            } else {
                v = owner_id(dset, cast(const(char) *) owner);
                if (v <= 0 && !subset) {
                    fprintf(stderr, "process_string_tables: invalid owner_id\n");
                    err = GretlError.data;
                }
            }
            if (v > 0) {
                strs = gretl_xml_get_strings_array(cur, doc, &n_strs,
                                                   0, &err);
                if (err) {
                    fprintf(stderr, "process_string_tables: get_strings_array "
                            ~ "gave error %d\n", err);
                } else {
                    st = series_table_new(strs, n_strs, &err);
                    if (err) {
                        strings_array_free(strs, n_strs);
                    } else {
                        series_attach_string_table(dset, v, st);
                    }
                }
            }
            free(owner);
        }
        if (!err) {
            cur = cur.next;
        }
    }

    return err;
}

private int process_panel_info(xmlNodePtr cur, DatasetRaw * dset,
                               int * repad)
{
    xmlChar * tmp;
    double sd0 = 0.0;
    int pd = 0;

    tmp = xmlGetProp(cur, cast(XUC) "group-names");
    if (tmp != null) {
        dset.pangrps = cast(char *) tmp;
    }

    tmp = xmlGetProp(cur, cast(XUC) "time-frequency");
    if (tmp != null) {
        pd = atoi(cast(const(char) *) tmp);
        free(tmp);
    }

    tmp = xmlGetProp(cur, cast(XUC) "time-start");
    if (tmp != null) {
        sd0 = atof(cast(const(char) *) tmp);
        free(tmp);
    }

    tmp = xmlGetProp(cur, cast(XUC) "skip-padding");
    if (tmp != null) {
        *repad = 1;
        free(tmp);
    }

    if (pd > 0 && sd0 > 0.0) {
        dset.panel_pd = pd;
        dset.panel_sd0 = sd0;
    }

    return 0;
}

private double get_gdt_version(xmlNodePtr node)
{
    xmlChar * tmp = xmlGetProp(node, cast(XUC) "version");
    double v = 1.0;

    if (tmp != null) {
        v = dot_atof(cast(char *) tmp);
        free(tmp);
    }

    return v;
}

private int xml_get_data_structure(xmlNodePtr node, int * dattype)
{
    xmlChar * tmp = xmlGetProp(node, cast(XUC) "type");
    int err = 0;

    if (tmp == null) {
        gretl_errmsg_set(_("Required attribute 'type' is missing from data file").ptr);
        err = 1;
    } else {
        if (!strcmp(cast(char *) tmp, "cross-section")) {
            *dattype = TsCodes.cross_section;
        } else if (!strcmp(cast(char *) tmp, "time-series")) {
            *dattype = TsCodes.time_series;
        } else if (!strcmp(cast(char *) tmp, "stacked-time-series")) {
            *dattype = TsCodes.stacked_time_series;
        } else if (!strcmp(cast(char *) tmp, "stacked-cross-section")) {
            *dattype = TsCodes.stacked_cross_section;
        } else {
            gretl_errmsg_set(_("Unrecognized type attribute for data file").ptr);
            err = 1;
        }
        free(tmp);
    }

    return err;
}

private int xml_get_data_frequency(xmlNodePtr node, int * pd, int * dattype)
{
    xmlChar * tmp = xmlGetProp(node, cast(XUC) "frequency");
    int err = 0;

    *pd = 1;

    if (tmp != null) {
        if (!strncmp(cast(char *) tmp, "special", 7)) {
            *dattype = TsCodes.special_time_series;
            if (sscanf(cast(char *) tmp + 7, ":%d", pd) == 1) {
                fprintf(stderr, "custom time series, frequency %d\n", *pd);
            } else {
                fprintf(stderr, "custom time series, using frequency 1\n");
            }
        } else if (sscanf(cast(char *) tmp, "%d", pd) != 1) {
            gretl_errmsg_set(_("Failed to parse data frequency").ptr);
            err = 1;
        }
        free(tmp);
    }

    return err;
}

private int xml_get_nobs(xmlNodePtr node, int * n)
{
    xmlChar * tmp = xmlGetProp(node, cast(XUC) "n");
    int err = 0;

    if (tmp != null) {
        *n = gretl_int_from_string(cast(const(char) *) tmp, &err);
        if (*n <= 0) {
            gretl_errmsg_sprintf("%d: invalid number of observations".ptr, *n);
            err = GretlError.data;
        }
        free(tmp);
    } else {
        *n = 0;
    }

    return err;
}

private int likely_calendar(const(char) * s)
{
    return strchr(s, '-') || strchr(s, '/');
}

private int xml_get_startobs(xmlNodePtr node, double * sd0, char * stobs,
                             int caldata)
{
    xmlChar * tmp = xmlGetProp(node, cast(XUC) "startobs");
    int err = 0;

    if (tmp != null) {
        char[OBSLEN] obstr = 0;
        int y, m, d;

        obstr[0] = '\0';
        strncat(obstr.ptr, cast(char *) tmp, OBSLEN - 1);
        gretl_charsub(obstr.ptr, ':', '.');

        if (sscanf(obstr.ptr, "%d/%d/%d", &y, &m, &d) == 3) {
            /* handle legacy gdt dates */
            gretl_charsub(obstr.ptr, '/', '-');
        }

        if (likely_calendar(obstr.ptr) && caldata) {
            uint32_t ed = get_epoch_day(obstr.ptr);

            if (ed <= 0) {
                err = 1;
            } else {
                *sd0 = ed;
            }
        } else {
            double x;

            if (sscanf(obstr.ptr, "%lf", &x) != 1) {
                err = 1;
            } else {
                *sd0 = x;
            }
        }

        if (err) {
            gretl_errmsg_set(_("Failed to parse startobs").ptr);
        } else {
            stobs[0] = '\0';
            strncat(stobs, obstr.ptr, OBSLEN - 1);
            colonize_obs(stobs);
        }

        free(tmp);
    }

    return err;
}

private int xml_get_endobs(xmlNodePtr node, char * endobs, int caldata)
{
    xmlChar * tmp = xmlGetProp(node, cast(XUC) "endobs");
    int err = 0;

    if (tmp != null) {
        char[OBSLEN] obstr = 0;

        obstr[0] = '\0';
        strncat(obstr.ptr, cast(char *) tmp, OBSLEN - 1);
        gretl_charsub(obstr.ptr, ':', '.');

        if (caldata) {
            int y, m, d;
            uint32_t ed;

            if (sscanf(obstr.ptr, "%d/%d/%d", &y, &m, &d) == 3) {
                /* handle legacy gdt file */
                gretl_charsub(obstr.ptr, '/', '-');
            }

            ed = get_epoch_day(obstr.ptr);
            if (ed <= 0) {
                err = 1;
            }
        } else {
            double x;

            if (sscanf(obstr.ptr, "%lf", &x) != 1) {
                err = 1;
            }
        }

        if (err) {
            gretl_errmsg_set(_("Failed to parse endobs").ptr);
        } else {
            endobs[0] = '\0';
            strncat(endobs, obstr.ptr, OBSLEN - 1);
            colonize_obs(endobs);
        }

        free(tmp);
    }

    return err;
}

private int gdt_binary_order(xmlNodePtr node)
{
    xmlChar * tmp = xmlGetProp(node, cast(XUC) "binary");
    int ret = 0;

    if (tmp != null) {
        if (!strcmp(cast(char *) tmp, "little-endian")) {
            ret = G_LITTLE_ENDIAN;
        } else if (!strcmp(cast(char *) tmp, "big-endian")) {
            ret = G_BIG_ENDIAN;
        }
        free(tmp);
    }

    return ret;
}

private int lag_from_label(int v, const(DatasetRaw) * dset, int * lag)
{
    const(char) * test = series_get_label(dset, v);
    char pm = 0;
    char[20] fmt = 0;
    char[VNAMELEN] vname = 0;
    int pv = 0;

    if (test == null) {
        return 0;
    }

    sprintf(fmt.ptr, "= %%%d[^(](t %%c %%d)", VNAMELEN - 1);

    if (sscanf(test, fmt.ptr, vname.ptr, &pm, lag) == 3) {
        pv = series_index(dset, vname.ptr);
        pv = (pv < dset.v)? pv : 0;
    }

    return pv;
}

private int dummy_child_from_label(int v, const(DatasetRaw) * dset)
{
    const(char) * test = series_get_label(dset, v);
    char[VNAMELEN] vname = 0;
    double val;
    int pv = 0;

    if (test == null) {
        return 0;
    }

    if (sscanf(test, _("dummy for %s = %lf").ptr, vname.ptr, &val) == 2 ||
        sscanf(test, "dummy for %s = %lf", vname.ptr, &val) == 2) {
        pv = series_index(dset, vname.ptr);
        pv = (pv < dset.v)? pv : 0;
    }

    return pv;
}

private void record_transform_info(DatasetRaw * dset, double gdtver)
{
    int i, p, pv;

    for (i=1; i<dset.v; i++) {
        if (series_get_transform(dset, i) == GretlCmdIndex.lags) {
            /* already handled */
            continue;
        }
        pv = lag_from_label(i, dset, &p);
        if (pv > 0) {
            series_set_parent(dset, i, dset.varname[pv]);
            series_set_transform(dset, i, GretlCmdIndex.lags);
            series_set_lag(dset, i, p);
        } else if (gdtver < 1.1) {
            pv = dummy_child_from_label(i, dset);
            if (pv > 0) {
                series_set_parent(dset, i, dset.varname[pv]);
                series_set_transform(dset, i, GretlCmdIndex.dummify);
            }
        }
    }
}

private void data_read_message(const(char) * fname, DatasetRaw * dset, Prn * prn)
{
    pprintf(prn, _("\nRead datafile %s\n").ptr, fname);
    pprintf(prn, (_("periodicity: %d, maxobs: %d\n")
                  ~ _("observations range: %s to %s\n")).ptr,
            (custom_time_series(dset))? 1 : dset.pd,
            dset.n, dset.stobs.ptr, dset.endobs.ptr);
    pputc(prn, '\n');
}

private c_long get_filesize(const(char) * fname)
{
    stat_t buf;
    int err;

    err = gretl_stat(fname, &buf);

    return (err)? -1 : cast(c_long) buf.st_size;
}

private int remedy_empty_data(DatasetRaw * dset)
{
    int add_idx = gretl_in_gui_mode();
    int i, err = 0;

    if (dset.Z[0] == null) {
        dset.Z[0] = cast(double *) malloc(dset.n * double.sizeof);
        if (dset.Z[0] == null) {
            return GretlError.alloc;
        }
        for (i=0; i<dset.n; i++) {
            dset.Z[0][i] = 1.0;
        }
    }

    if (add_idx) {
        err = dataset_add_series(dset, 1);
    }

    if (!err) {
        int j;

        dset.t1 = 0;
        dset.t2 = dset.n - 1;
        if (add_idx) {
            strcpy(dset.varname[1], "index");
            series_set_label(dset, 1, _("index variable").ptr);
            for (j=0; j<dset.n; j++) {
                dset.Z[1][j] = cast(double) (j + 1);
            }
        }
    }

    return err;
}

private void check_for_daily_date_strings(DatasetRaw * dset)
{
    int m, d, n;
    int y1 = 0, y2 = 0;
    int oldfmt = 0;

    n = sscanf(dset.S[0], YMD_READ_FMT.ptr, &y1, &m, &d);
    if (n != 3) {
        oldfmt = 1;
        n = sscanf(dset.S[0], "%d/%d/%d", &y1, &m, &d);
    }

    if (n == 3) {
        int k = dset.n - 1;

        if (oldfmt) {
            n = sscanf(dset.S[k], "%d/%d/%d", &y2, &m, &d);
        } else {
            n = sscanf(dset.S[k], YMD_READ_FMT.ptr, &y2, &m, &d);
        }
    }

    if (n == 3 && y2 >= y1) {
        dset.markers = cast(char) DatasetMarkerType.daily_date_strings;
    }
}

private int replace_panel_padding(DatasetRaw * dset)
{
    int uv = dset.v - 2;
    int tv = dset.v - 1;
    int err = 0;

    if (!strcmp(dset.varname[uv], "unit__") &&
        !strcmp(dset.varname[tv], "time__")) {
        err = set_panel_structure_from_vars(uv, tv, dset);
        if (!err) {
            dataset_drop_last_variables(dset, 2);
        }
    }

    return err;
}

private void future_datafile_warning(double v1, double v2)
{
    const(char) * fmt =
        (_("The version of this datafile (%g) is higher than that fully\n")
         ~ _("supported by this build of gretl (%g). Some features may not be\n")
         ~ _("correctly recognized.\n")).ptr;

    gretl_warnmsg_sprintf(fmt, v1, v2);
}

private int real_read_gdt(const(char) * fname, const(char) * srcname,
                          DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    DatasetRaw * tmpset;
    xmlDocPtr doc = null;
    xmlNodePtr cur;
    int gotvars = 0, err = 0;
    int caldata = 0, repad = 0;
    int[2] sample = 0;
    int readsmpl = 0;
    double gdtversion = 1.0;
    double myversion;
    int in_c_locale = 0;
    int gz, binary = 0;
    c_long fsz;

    gretl_error_clear();
    fsz = get_filesize(fname);
    gz = is_gzipped(fname);

    if (fsz < 0) {
        return GretlError.fopen;
    } else if (fsz > 100000) {
        fprintf(stderr, "%s %.0f Kbytes %s...\n",
                gz ? "Uncompressing".ptr : "Reading".ptr,
                cast(double) fsz / 1024.0, "of data".ptr);
    }

    tmpset = datainfo_new();
    if (tmpset == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    err = gretl_xml_open_doc_root(fname, "gretldata", &doc, &cur);
    if (err) {
        goto bailout;
    }

    gdtversion = get_gdt_version(cur);
    myversion = dot_atof(GRETLDATA_VERSION.ptr);
    if (gdtversion > myversion) {
        if (gdtversion == 1.5 && myversion == 1.4) {
            ; /* this shouldn't be a problem */
        } else {
            future_datafile_warning(gdtversion, myversion);
        }
    }

    /* optional */
    readsmpl = gretl_xml_get_prop_as_bool(cur, "mpi-transfer");

    /* optional */
    gretl_xml_get_prop_as_uint64(cur, "rseed", &tmpset.rseed);

    /* optional */
    gretl_xml_get_prop_as_string(cur, "mapfile", &tmpset.mapfile);

    /* required only if the <observations> element is absent */
    err = xml_get_nobs(cur, &tmpset.n);
    if (err) {
        goto bailout;
    }

    /* set some required datainfo parameters */

    err = xml_get_data_structure(cur, &tmpset.structure);
    if (err) {
        goto bailout;
    }

    err = xml_get_data_frequency(cur, &tmpset.pd, &tmpset.structure);
    if (err) {
        goto bailout;
    }

    gretl_push_c_numeric_locale();
    in_c_locale = 1;

    strcpy(tmpset.stobs.ptr, "1");
    caldata = dataset_is_daily(tmpset) || dataset_is_weekly(tmpset);

    err = xml_get_startobs(cur, &tmpset.sd0, tmpset.stobs.ptr, caldata);
    if (err) {
        goto bailout;
    }

    tmpset.endobs[0] = '\0';
    caldata = calendar_data(tmpset);

    err = xml_get_endobs(cur, tmpset.endobs.ptr, caldata);
    if (err) {
        goto bailout;
    }

    binary = gdt_binary_order(cur);

    /* Now walk the tree */
    cur = cur.children;
    while (cur != null && !err) {
        if (!xmlStrcmp(cur.name, cast(XUC) "description")) {
            tmpset.descrip = cast(char *)
                xmlNodeListGetString(doc, cur.children, 1);
        } else if (readsmpl && !xmlStrcmp(cur.name, cast(XUC) "sample")) {
            err = process_sample(cur, sample.ptr, tmpset);
            if (err) {
                fprintf(stderr, "error processing sample range\n");
            } else {
                readsmpl = 2;
            }
        } else if (!xmlStrcmp(cur.name, cast(XUC) "variables")) {
            err = process_varlist(cur, tmpset, 0);
            if (err) {
                fprintf(stderr, "error processing varlist\n");
            } else {
                gotvars = 1;
            }
        } else if (!xmlStrcmp(cur.name, cast(XUC) "observations")) {
            if (!gotvars) {
                gretl_errmsg_set(_("Variables information is missing").ptr);
                err = 1;
            } else {
                double dsize = (opt & Gretlopt.OPT_B)? cast(double) fsz : 0;

                err = read_observations(doc, cur, tmpset, dsize,
                                        binary, gdtversion, fname);
                if (err) {
                    fprintf(stderr, "error %d in read_observations\n", err);
                }
            }
        } else if (!xmlStrcmp(cur.name, cast(XUC) "string-tables")) {
            if (!gotvars) {
                gretl_errmsg_set(_("Variables information is missing").ptr);
                err = GretlError.data;
            } else {
                err = process_string_tables(doc, cur, tmpset, 0);
                if (err) {
                    fprintf(stderr, "error %d processing string tables\n", err);
                }
            }
        } else if (!xmlStrcmp(cur.name, cast(XUC) "panel-info")) {
            if (!gotvars) {
                gretl_errmsg_set(_("Variables information is missing").ptr);
                err = GretlError.data;
            } else {
                err = process_panel_info(cur, tmpset, &repad);
                if (err) {
                    fprintf(stderr, "error %d processing panel info\n", err);
                }
            }
        }
        if (!err) {
            cur = cur.next;
        }
    }

    if (!err && !gotvars) {
        gretl_errmsg_set(_("Variables information is missing").ptr);
        err = 1;
    }

    if (!err && caldata && tmpset.S != null) {
        check_for_daily_date_strings(tmpset);
    }

    if (!err && repad) {
        err = replace_panel_padding(tmpset);
    }

    if (!err) {
        if (srcname == null) {
            srcname = fname;
        }
        data_read_message(srcname, tmpset, prn);
        err = merge_or_replace_data(dset, &tmpset,
                                    get_merge_opts(opt), prn);
    }

 bailout:

    if (in_c_locale) {
        gretl_pop_c_numeric_locale();
    }

    if (doc != null) {
        xmlFreeDoc(doc);
    }

    if (dset != null) {
        /* pre-process stacked cross-sectional panels: put into canonical
           stacked time series form
        */
        if (!err && dset.structure == TsCodes.stacked_cross_section) {
            err = switch_panel_orientation(dset);
        }
        if (!err && dset.v == 1) {
            err = remedy_empty_data(dset);
        }
        if (!err && gdtversion < 1.2) {
            record_transform_info(dset, gdtversion);
        }
    }

    if (err && tmpset != null) {
        destroy_dataset(tmpset);
    }

    if (err && err != GretlError.alloc) {
        /* ensure we don't just show "unspecified error" */
        gchar * msg;

        if (err == GretlError.data) {
            msg = g_strdup_printf(_("'%s': malformed gdt file").ptr, fname);
        } else {
            msg = g_strdup_printf(_("Couldn't read '%s'").ptr, fname);
        }

        gretl_errmsg_ensure(msg);
        g_free(msg);
    }

    if (!err && readsmpl == 2) {
        /* the gdt file contained a sample-range spec */
        dset.t1 = sample[0];
        dset.t2 = sample[1];
    }

    return err;
}

private int real_read_gdt_subset(const(char) * fname,
                                 DatasetRaw * dset,
                                 int * vlist,
                                 Gretlopt opt)
{
    DatasetRaw * tmpset;
    xmlDocPtr doc = null;
    xmlNodePtr cur;
    double gdtversion = 1.0;
    int gotvars = 0, gotobs = 0;
    int caldata = 0;
    int in_c_locale = 0;
    int binary = 0;
    int fullv = 0;
    int err = 0;

    gretl_error_clear();

    tmpset = datainfo_new();
    if (tmpset == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    err = gretl_xml_open_doc_root(fname, "gretldata", &doc, &cur);
    if (err) {
        goto bailout;
    }

    /* set some datainfo parameters */

    err = xml_get_nobs(cur, &tmpset.n);
    if (err) {
        goto bailout;
    }

    err = xml_get_data_structure(cur, &tmpset.structure);
    if (err) {
        goto bailout;
    }

    err = xml_get_data_frequency(cur, &tmpset.pd, &tmpset.structure);
    if (err) {
        goto bailout;
    }

    gdtversion = get_gdt_version(cur);

    gretl_push_c_numeric_locale();
    in_c_locale = 1;

    strcpy(tmpset.stobs.ptr, "1");
    caldata = dataset_is_daily(tmpset) || dataset_is_weekly(tmpset);

    err = xml_get_startobs(cur, &tmpset.sd0, tmpset.stobs.ptr, caldata);
    if (err) {
        goto bailout;
    }

    tmpset.endobs[0] = '\0';
    caldata = calendar_data(tmpset);

    err = xml_get_endobs(cur, tmpset.endobs.ptr, caldata);
    if (err) {
        goto bailout;
    }

    binary = gdt_binary_order(cur);

    /* Now walk the tree */
    cur = cur.children;
    while (cur != null && !err) {
        if (!xmlStrcmp(cur.name, cast(XUC) "variables")) {
            err = process_varlist_subset(cur, tmpset, &fullv, vlist);
            if (!err) {
                gotvars = 1;
            }
        } else if (!xmlStrcmp(cur.name, cast(XUC) "observations")) {
            if (!gotvars) {
                gretl_errmsg_set(_("Variables information is missing").ptr);
                err = GretlError.data;
            } else {
                err = read_observations_subset(doc, cur, tmpset,
                                               binary, gdtversion,
                                               fname, fullv, vlist,
                                               opt);
            }
            if (!err) {
                gotobs = 1;
            }
        } else if (!xmlStrcmp(cur.name, cast(XUC) "string-tables")) {
            if (!gotvars) {
                gretl_errmsg_set(_("Variables information is missing").ptr);
                err = 1;
            } else {
                err = process_string_tables(doc, cur, tmpset, 1);
            }
        }
        if (!err) {
            cur = cur.next;
        }
    }

    if (!err && !gotvars) {
        gretl_errmsg_set(_("Variables information is missing").ptr);
        err = 1;
    }

    if (!err && !gotobs) {
        gretl_errmsg_set(_("No observations were found").ptr);
        err = 1;
    }

    if (!err && caldata && tmpset.S != null) {
        check_for_daily_date_strings(tmpset);
    }

 bailout:

    if (in_c_locale) {
        gretl_pop_c_numeric_locale();
    }

    if (doc != null) {
        xmlFreeDoc(doc);
    }

    if (!err) {
        *dset = *tmpset;
        free(tmpset);
    } else {
        destroy_dataset(tmpset);
    }

    return err;
}

private int real_read_gdt_varnames(const(char) * fname,
                                   char *** vnames,
                                   int * nvars)
{
    DatasetRaw * tmpset;
    xmlDocPtr doc = null;
    xmlNodePtr cur;
    int gotvars = 0;
    int caldata = 0;
    int in_c_locale = 0;
    int err = 0;

    gretl_error_clear();

    tmpset = datainfo_new();
    if (tmpset == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    err = gretl_xml_open_doc_root(fname, "gretldata", &doc, &cur);
    if (err) {
        goto bailout;
    }

    /* set some datainfo parameters */

    err = xml_get_data_structure(cur, &tmpset.structure);
    if (err) {
        goto bailout;
    }

    err = xml_get_data_frequency(cur, &tmpset.pd, &tmpset.structure);
    if (err) {
        goto bailout;
    }

    gretl_push_c_numeric_locale();
    in_c_locale = 1;

    strcpy(tmpset.stobs.ptr, "1");
    caldata = dataset_is_daily(tmpset) || dataset_is_weekly(tmpset);

    err = xml_get_startobs(cur, &tmpset.sd0, tmpset.stobs.ptr, caldata);
    if (err) {
        goto bailout;
    }

    tmpset.endobs[0] = '\0';
    caldata = calendar_data(tmpset);

    err = xml_get_endobs(cur, tmpset.endobs.ptr, caldata);
    if (err) {
        goto bailout;
    }

    /* Now walk the tree */
    cur = cur.children;
    while (cur != null && !err) {
        if (!xmlStrcmp(cur.name, cast(XUC) "variables")) {
            err = process_varlist(cur, tmpset, 1);
            if (!err) {
                gotvars = 1;
            }
            break;
        }
        cur = cur.next;
    }

    if (!err && !gotvars) {
        gretl_errmsg_set(_("Variables information is missing").ptr);
        err = 1;
    }

 bailout:

    if (in_c_locale) {
        gretl_pop_c_numeric_locale();
    }

    if (doc != null) {
        xmlFreeDoc(doc);
    }

    if (!err) {
        *vnames = tmpset.varname;
        *nvars = tmpset.v;
        tmpset.varname = null;
    }

    destroy_dataset(tmpset);

    return err;
}

private int read_gbin(const(char) * fname, DatasetRaw * dset,
                      Gretlopt opt, Prn * prn)
{
    extern (C) int function(const(char) *, DatasetRaw *, Gretlopt, Prn * prn) reader;
    int err = 0;

    reader = cast(typeof(reader)) get_plugin_function("purebin_read_data");

    if (reader == null) {
        err = 1;
    } else {
        err = reader(fname, dset, opt, prn);
    }

    return err;
}

private int read_gbin_subset(const(char) * fname, DatasetRaw * dset,
                             int * vlist, Gretlopt opt)
{
    extern (C) int function(const(char) *, DatasetRaw *, int *, Gretlopt) reader;
    int err = 0;

    reader = cast(typeof(reader)) get_plugin_function("purebin_read_subset");

    if (reader == null) {
        err = 1;
    } else {
        err = reader(fname, dset, vlist, opt);
    }

    return err;
}

private int read_gbin_varnames(const(char) * fname, char *** pS,
                               int * pns)
{
    extern (C) int function(const(char) *, char ***, int *) reader;
    int err = 0;

    reader = cast(typeof(reader)) get_plugin_function("purebin_read_varnames");

    if (reader == null) {
        err = 1;
    } else {
        err = reader(fname, pS, pns);
    }

    return err;
}

private int is_purebin_file(const(char) * fname)
{
    FILE * fp = gretl_fopen(fname, "rb");
    int ret = 0;

    if (fp != null) {
        char[16] buf = 0;

        if (fread(buf.ptr, 1, 14, fp) == 14 &&
            !strcmp(buf.ptr, "gretl-purebin")) {
            ret = 1;
        }
        fclose(fp);
    }

    return ret;
}

/**
 * gretl_read_gdt:
 * @fname: name of file to open for reading.
 * @dset: dataset struct.
 * @opt: use OPT_B to display gui progress bar; may also
 * use OPT_T when appending to panel data (see the "append"
 * command in the gretl manual). Otherwise use OPT_NONE.
 * @prn: where any messages should be written.
 *
 * Read data from native file into gretl's workspace.
 *
 * Returns: 0 on successful completion, non-zero otherwise.
 */

int gretl_read_gdt(const(char) * fname, DatasetRaw * dset,
                   Gretlopt opt, Prn * prn)
{
    int gdtb = has_suffix(fname, ".gdtb");

    if (gdtb && is_purebin_file(fname)) {
        return read_gbin(fname, dset, opt, prn);
    } else if (gdtb) {
        /* zipfile with gdt + binary */
        gchar * zdir;
        int id = -1;
        int err;

        if (id >= 0) {
            zdir = g_strdup_printf("%stmp%d-unzip", gretl_dotdir(), id);
        } else {
            zdir = g_strdup_printf("%stmp-unzip", gretl_dotdir());
        }
        err = gretl_mkdir(zdir);

        if (!err) {
            err = gretl_unzip_into(fname, zdir);
            if (err) {
                gretl_errmsg_ensure("Problem opening data file");
            } else {
                char[FILENAME_MAX] xmlfile = 0;

                gretl_build_path(xmlfile.ptr, zdir, "data.xml".ptr, null);
                err = real_read_gdt(xmlfile.ptr, fname, dset, opt, prn);
            }
            gretl_deltree(zdir);
        }

        g_free(zdir);
        return err;
    } else {
        /* plain XML file */
        return real_read_gdt(fname, null, dset, opt, prn);
    }
}

/**
 * gretl_read_gdt_subset:
 * @fname: name of file to open for reading.
 * @dset: dataset struct.
 * @vlist: list of series to extract.
 * @opt: may include OPT_M to retrieve the observation
 * markers associated with the data, if any.
 *
 * Read specified series from native file into @dset,
 * which should be "empty" on input.
 *
 * Returns: 0 on successful completion, non-zero otherwise.
 */

int gretl_read_gdt_subset(const(char) * fname, DatasetRaw * dset,
                          int * vlist, Gretlopt opt)
{
    int gdtb = has_suffix(fname, ".gdtb");
    int err = 0;

    if (gdtb && is_purebin_file(fname)) {
        err = read_gbin_subset(fname, dset, vlist, opt);
    } else if (gdtb) {
        /* legacy: zipfile with gdt + binary */
        gchar * zdir;
        int inner_err; /* C: "int err;" here shadows the outer err (see bugs.md) */

        zdir = g_strdup_printf("%stmp-unzip", gretl_dotdir());
        inner_err = gretl_mkdir(zdir);

        if (!inner_err) {
            inner_err = gretl_unzip_into(fname, zdir);
            if (inner_err) {
                gretl_errmsg_ensure("Problem opening data file");
            } else {
                char[FILENAME_MAX] xmlfile = 0;

                gretl_build_path(xmlfile.ptr, zdir, "data.xml".ptr, null);
                inner_err = real_read_gdt_subset(xmlfile.ptr, dset, vlist, opt);
            }
            gretl_deltree(zdir);
        }

        g_free(zdir);
    } else {
        /* plain XML gdt file */
        err = real_read_gdt_subset(fname, dset, vlist, opt);
    }

    return err;
}

/**
 * gretl_read_native_varnames:
 * @fname: name of native gretl data file to open for reading.
 * @vnames: location to receive array of series names.
 * @nvars: location to receive the number of series.
 *
 * Read the array of series names from the specified file.
 *
 * Returns: 0 on successful completion, non-zero otherwise.
 */

int gretl_read_native_varnames(const(char) * fname,
                               char *** vnames,
                               int * nvars)
{
    int gdtb = has_suffix(fname, ".gdtb");
    int err = 0;

    if (gdtb && is_purebin_file(fname)) {
        err = read_gbin_varnames(fname, vnames, nvars);
    } else if (gdtb) {
        /* zipfile with gdt + binary */
        gchar * zdir;
        int inner_err; /* C: "int err;" here shadows the outer err (see bugs.md) */

        zdir = g_strdup_printf("%stmp-unzip", gretl_dotdir());
        inner_err = gretl_mkdir(zdir);

        if (!inner_err) {
            inner_err = gretl_unzip_into(fname, zdir);
            if (inner_err) {
                gretl_errmsg_ensure("Problem opening data file");
            } else {
                char[FILENAME_MAX] xmlfile = 0;

                gretl_build_path(xmlfile.ptr, zdir, "data.xml".ptr, null);
                inner_err = real_read_gdt_varnames(xmlfile.ptr, vnames,
                                                   nvars);
            }
            gretl_deltree(zdir);
        }

        g_free(zdir);
    } else {
        /* plain XML file */
        err = real_read_gdt_varnames(fname, vnames, nvars);
    }

    return err;
}

/**
 * gretl_get_gdt_description:
 * @fname: name of file to try.
 *
 * Read data description for gretl native data file.
 *
 * Returns: buffer containing description, or NULL on failure.
 */

char * gretl_get_gdt_description(const(char) * fname, int * err)
{
    xmlDocPtr doc;
    xmlNodePtr cur;
    xmlChar * buf = null;

    gretl_error_clear();

    if (!has_suffix(fname, ".gdt")) {
        gretl_errmsg_set("gretl_get_gdt_description() is only for gdt files");
        *err = GretlError.invarg;
        return null;
    }

    *err = gretl_xml_open_doc_root(fname, "gretldata", &doc, &cur);
    if (*err) {
        return null;
    }

    cur = cur.children;
    while (cur != null) {
        if (!xmlStrcmp(cur.name, cast(XUC) "description")) {
            buf = xmlNodeListGetString(doc, cur.children, 1);
            break;
        }
        cur = cur.next;
    }

    if (buf == null) {
        gretl_errmsg_set("No description was found");
        *err = GretlError.data;
    }

    xmlFreeDoc(doc);

    return cast(char *) buf;
}

private char * gretl_xml_get_doc_type(const(char) * fname, int * err)
{
    xmlDocPtr doc;
    xmlNodePtr node;
    int options;
    char * ret = null;

    options = XML_PARSE_HUGE;

    doc = xmlReadFile(fname, null, options);

    if (doc == null) {
        gretl_errmsg_sprintf(_("xmlReadFile failed on %s").ptr, fname);
        *err = read_file_fail();
    } else {
        node = xmlDocGetRootElement(doc);
        if (node == null) {
            gretl_errmsg_sprintf(_("%s: empty document").ptr, fname);
            *err = GretlError.data;
        } else {
            ret = gretl_strdup(cast(char *) node.name);
            if (ret == null) {
                *err = GretlError.alloc;
            }
        }
    }

    if (doc != null) {
        xmlFreeDoc(doc);
    }

    return ret;
}

/* This is called in response to the "include" command in the CLI
   program, the GUI program, and in interact.c, if we detect that the
   named file is XML.
*/

int load_XML_functions_file(const(char) * fname, Gretlopt opt, Prn * prn)
{
    char * rootname = null;
    int err = 0;

    rootname = gretl_xml_get_doc_type(fname, &err);
    if (err) {
        return err;
    }

    if (!strcmp(rootname, "gretl-functions")) {
        if (has_suffix(fname, ".gfn")) {
            err = include_gfn(fname, opt, prn);
        } else {
            err = read_session_functions_file(fname);
        }
    } else {
        err = GretlError.data;
    }

    free(rootname);

    return err;
}

void gretl_xml_init()
{
    xmlInitParser();
}

void gretl_xml_cleanup()
{
    xmlCleanupParser();
}
