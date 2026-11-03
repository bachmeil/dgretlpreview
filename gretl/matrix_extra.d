/* Port of lib/src/matrix_extra.c: conversions between gretl_matrix and other data types, matrix printing, reading
   and writing matrices from/to file, masking/cutting of rows and columns, quadrature rules and a few more
   matrix operations (the header matrix_extra.h is matrix_extrah.d).

   Differences from the C:
   - the Windows branches (win32_fscan_nonfinite, win32_xna_out) and the big-endian byte-swapping branches (swap_bytes.h)
     are not ported: the code assumes a little-endian host;
   - gretl_matrix_write_constructor() builds its string with D's std.array/std.format instead of a glib GString;
   - the `sign` macro of the quadrature code is a function, MAX/MIN come from std.algorithm.
   Dependencies not yet ported are declared extern (C) below, in the style of the other ports. */

module gretl.matrix_extra;


import core.stdc.errno : errno;
import core.stdc.math : ceil, fabs, floor, log, log10, sqrt, tgamma;
import core.stdc.stdio : EOF, FILE, SEEK_SET, fclose, fgetc, fprintf, fputc, fputs, fread, fscanf, fseek, ftell, fwrite,
    snprintf, sprintf, sscanf, stderr, ungetc;
import core.stdc.stdlib : calloc, free, malloc, strtod, strtol;
import core.stdc.string : memcpy, strcat, strchr, strcmp, strcpy, strlen, strncat, strncpy, strstr;
import core.stdc.ctype : isalpha, isdigit;

import std.algorithm : max, min;
import std.math : PI, SQRT2;

import gretl.dataio : ntolabel, write_data;
import gretl.dataset : create_auxiliary_dataset, destroy_dataset, series_record_display_name;
import gretl.describe : gretl_mean;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_errors : gretl_errmsg_set, gretl_errmsg_sprintf;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_intl : gretl_pop_c_numeric_locale, gretl_push_c_numeric_locale;
import gretl.gretl_matrix;
import gretl.gretl_matrixh;
import gretl.estimate : makevcv;
import gretl.gretl_paths : gretl_build_path, gretl_dotdir, gretl_fopen, gretl_maybe_prepend_dir, gretl_maybe_switch_dir,
    gretl_remove;
import gretl.gretl_prn : Prn, gretl_gzip_print_new, gretl_print_destroy, gretl_print_new, gretl_print_new_with_filename, pprintf,
    pputc, pputs;
import gretl.gretl_prnh : PrnType;
import gretl.gretl_utils : doubles_array_new, gretl_gunzip, gretl_ispositive, ijton;
import gretl.graphingh : FitType;
import gretl.libgretlh;
import gretl.libset : get_csv_na_write_string, get_data_export_delimiter, libset_get_bool;
import gretl.libseth : SetKey;
import gretl.missing : list_adjust_sample;
import gretl.matrix_extrah : MMissingCode, QuadMethod;
import gretl.missingh : NADBL, na;
import gretl.options : get_optval_string;
import gretl.plugins : get_plugin_function;
import gretl.strutils : gretl_delchar, gretl_normalize_varname, gretl_strdup, has_suffix, numeric_string, strings_array_free,
    strings_array_new;

/* unported dependencies: per-file declarations from the 2026c headers */
extern (C) {
}

import std.array : appender;
import std.format : format;
import gretl.gstrfuncs : g_strdup, g_utf8_strlen;
import gretl.gretl_cmatrix : gretl_cmatrix_print_range, gretl_cmatrix_printf;
import gretl.printout : bufspace, max_obs_marker_length;
import gretl.usermat : get_matrix_by_name;

/**
 * gretl_vector_from_array:
 * @x: pointer to array of elements.
 * @n: number of elements.
 * @mod: modifier flag: either %GRETL_MOD_NONE, or %GRETL_MOD_SQUARE
 * to use the squares of the elements of @x.
 *
 * Returns: pointer to a newly allocated gretl_vector containing
 * the elements of x (or their squares), or NULL on failure.
 * Missing values in @x are skipped.
 */

GretlVector * gretl_vector_from_array(const(double) * x, int n, GretlMatrixMod mod)
{
    GretlMatrixRaw * v;
    double xi;

    v = gretl_column_vector_alloc(n);

    if (v != null) {
        int i = 0, j = 0;

        while (j < n) {
            xi = x[i++];
            if (!na(xi)) {
                if (mod == GretlMatrixMod.square) {
                    v.val[j] = xi * xi;
                } else {
                    v.val[j] = xi;
                }
                j++;
            }
        }
    }

    return v;
}

/**
 * gretl_vector_from_series:
 * @x: series from data array.
 * @t1: starting observation.
 * @t2: ending observation.
 *
 * Returns: a newly allocated gretl_vector containing the values
 * of the given data series for the given range, or NULL on failure.
 * Any missing values in the input array are preserved as NaNs in
 * the output.
 */

GretlVector * gretl_vector_from_series(const(double) * x, int t1, int t2)
{
    GretlMatrixRaw * v = null;
    int n = t2 - t1 + 1;

    if (n > 0) {
        size_t sz = n * (*x).sizeof;

        v = gretl_column_vector_alloc(n);
        if (v != null) {
            memcpy(v.val, x + t1, sz);
            gretl_matrix_set_t1(v, t1);
            gretl_matrix_set_t2(v, t2);
        }
    }

    return v;
}

/**
 * gretl_matrix_from_2d_array:
 * @X: two-dimensional array of doubles.
 * @rows: number of rows in target matrix.
 * @cols: number of columns in target matrix.
 *
 * Returns: allocated gretl_matrix, the elements of which are set to
 * the values in @X, or NULL on allocation failure.
 */

GretlMatrixRaw * gretl_matrix_from_2d_array(const(double) ** X, int rows, int cols)
{
    GretlMatrixRaw * m;
    int i, j, k = 0;

    m = gretl_matrix_alloc(rows, cols);

    if (m != null) {
        for (j=0; j<cols; j++) {
            for (i=0; i<rows; i++) {
                m.val[k++] = X[j][i];
            }
        }
    }

    return m;
}

/**
 * gretl_matrix_from_scalar:
 * @x: scalar to be "promoted".
 *
 * Returns: allocated 1 x 1 gretl_matrix, the single element
 * of which is set to @x, or NULL on allocation failure.
 */

GretlMatrixRaw * gretl_matrix_from_scalar(double x)
{
    GretlMatrixRaw * m = gretl_matrix_alloc(1, 1);

    m.val[0] = x;

    return m;
}

int count_selection(const(char) * s, int n)
{
    int i, c = 0;

    for (i=0; i<n; i++) {
        if (s[i] != 0) c++;
    }

    return c;
}

/**
 * gretl_vcv_matrix_from_model:
 * @pmod: pointer to model
 * @select: char array indicating which rows and columns to select
 * (or NULL for the full matrix).
 * @err: location to receive error code.
 *
 * Produces all or part of the covariance matrix for @pmod
 * in the form of a gretl_matrix.  Storage is allocated, to be freed
 * by the caller.  If @select is not NULL, it should be an array
 * with non-zero elements in positions corresponding to the
 * desired rows (and columns), and zero elements otherwise.
 *
 * Returns: the covariance matrix, or NULL on error.
 */

GretlMatrixRaw * gretl_vcv_matrix_from_model(Model * pmod, const(char) * select, int * err)
{
    GretlMatrixRaw * vcv;
    int i, j, idx, nc;
    int ii, jj;
    int k = pmod.ncoeff;

    /* first ensure the model actually has a vcv */
    *err = makevcv(pmod, pmod.sigma);
    if (*err) {
        fprintf(stderr, "gretl_vcv_matrix_from_model: no $vcv!\n".ptr);
        return null;
    }

    if (select == null) {
        nc = k;
    } else {
        nc = count_selection(select, k);
    }

    if (nc == 0) {
        *err = GretlError.data;
        return null;
    }

    vcv = gretl_matrix_alloc(nc, nc);
    if (vcv == null) {
        *err = GretlError.alloc;
        return null;
    }

    ii = 0;
    for (i=0; i<k; i++) {
        if (select != null && !select[i]) {
            continue;
        }
        jj = 0;
        for (j=0; j<=i; j++) {
            if (select != null && !select[j]) {
                continue;
            }
            idx = ijton(i, j, pmod.ncoeff);
            gretl_matrix_set(vcv, ii, jj, pmod.vcv[idx]);
            if (jj != ii) {
                gretl_matrix_set(vcv, jj, ii, pmod.vcv[idx]);
            }
            jj++;
        }
        ii++;
    }

    return vcv;
}

/**
 * gretl_coeff_vector_from_model:
 * @pmod: pointer to model
 * @select: char array indicating which rows to select
 * (or NULL for the full vector).
 * @err: location to receive error code.
 *
 * Produces all or part of the coefficient vector for @pmod
 * in the form of a gretl column vector.  Storage is allocated, to be freed
 * by the caller.  If @select is non-NULL, it should be an array
 * with non-zero elements in positions corresponding to the
 * desired rows and zero elements otherwise.
 *
 * Returns: the coefficient vector, or NULL on error.
 */

GretlVector * gretl_coeff_vector_from_model(const(Model) * pmod, const(char) * select, int * err)
{
    GretlVector * b;
    int i, j, nc;
    int k = pmod.ncoeff;

    if (select == null) {
        nc = k;
    } else {
        nc = count_selection(select, k);
    }

    if (nc == 0) {
        *err = GretlError.data;
        return null;
    }

    b = gretl_column_vector_alloc(nc);
    if (b == null) {
        *err = GretlError.alloc;
        return null;
    }

    j = 0;
    for (i=0; i<k; i++) {
        if (select != null && !select[i]) {
            continue;
        }
        b.val[j++] = pmod.coeff[i];
    }

    return b;
}

/**
 * gretl_covariance_matrix_from_varlist:
 * @list: list of variables by ID number.
 * @dset: dataset struct.
 * @means: pointer to pick up vector of means, or NULL to discard.
 * @errp: pointer to receive non-zero error code in case of
 * failure, or NULL.
 *
 * Returns: the variance-covariance matrix of the listed variables
 * (over the currently defined data sample), or NULL in case of
 * failure.
 */

GretlMatrixRaw * gretl_covariance_matrix_from_varlist(const(int) * list,
                                                   const(DatasetRaw) * dset,
                                                   GretlMatrixRaw ** means,
                                                   int * errp)
{
    GretlMatrixRaw * V;
    GretlVector * xbar;
    int k = list[0];
    int t, i, j, vi, vj;
    double vv, x, y;
    int n, err = 0;

    V = gretl_matrix_alloc(k, k);
    xbar = gretl_vector_alloc(k);

    if (V == null || xbar == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    for (i=0; i<k && !err; i++) {
        xbar.val[i] = gretl_mean(dset.t1, dset.t2, dset.Z[list[i+1]]);
        if (na(xbar.val[i])) {
            err = GretlError.data;
        }
    }

    for (i=0; i<k && !err; i++) {
        vi = list[i+1];
        for (j=i; j<k; j++) {
            vj = list[j+1];
            vv = 0.0;
            n = 0;
            for (t=dset.t1; t<=dset.t2; t++) {
                x = dset.Z[vi][t];
                y = dset.Z[vj][t];
                if (na(x) || na(y)) {
                    continue;
                }
                vv += (x - xbar.val[i]) * (y - xbar.val[j]);
                n++;
            }
            if (n < 2) {
                err = GretlError.data;
                vv = NADBL;
            } else {
                vv /= (n - 1); /* or plain n? */
            }
            gretl_matrix_set(V, i, j, vv);
            gretl_matrix_set(V, j, i, vv);
        }
    }

 bailout:

    if (means != null && !err) {
        *means = xbar;
    } else {
        gretl_vector_free(xbar);
    }

    if (err) {
        gretl_matrix_free(V);
        V = null;
        if (errp != null) {
            *errp = err;
        }
    }

    return V;
}

/**
 * gretl_matrix_row_to_array:
 * @m: source matrix.
 * @i: the row from which values should be copied.
 * @x: array of doubles large enough to hold a row from @m.
 *
 * Copies the values from row @i of matrix @m into the array
 * @x, which should already be allocated to the correct size.
 *
 * Returns: 0 on sucess, non-zero if the row is out of bounds.
 */

int gretl_matrix_row_to_array(const(GretlMatrixRaw) * m, int i, double * x)
{
    int j, err = 0;

    if (i < 0 || i >= gretl_matrix_rows(m)) {
        err = GretlError.data;
    } else {
        for (j=0; j<m.cols; j++) {
            x[j] = gretl_matrix_get(m, i, j);
        }
    }

    return err;
}

/**
 * gretl_matrix_get_columns:
 * @m: source matrix.
 * @err: location to receive error code.
 *
 * Constructs a two-dimensional array in which each sub-array
 * is a pointer to a column of @m. The content of these arrays
 * belongs to @m and must not be freed; only the returned
 * pointer itself should be freed.
 *
 * Returns: the constructed array on success or NULL on failure.
 */

double ** gretl_matrix_get_columns(const(GretlMatrixRaw) * m, int * err)
{
    double ** X = null;

    if (gretl_is_null_matrix(m)) {
        *err = GretlError.data;
    } else {
        double * val = cast(double *) m.val;
        int j;

        X = doubles_array_new(m.cols, 0);

        if (X == null) {
            *err = GretlError.alloc;
        } else {
            for (j=0; j<m.cols; j++) {
                X[j] = val;
                val += m.rows;
            }
        }
    }

    return X;
}

int get_mask_count(const(char) * mask, int n)
{
    int i, k = 0;

    for (i=0; i<n; i++) {
        if (mask[i]) k++;
    }

    return k;
}

void add_dataset_colnames(GretlMatrixRaw * M, const(int) * list, const(DatasetRaw) * dset)
{
    int i, vi, nv;
    char ** S = null;
    int err = 0;

    if (list != null) {
        nv = list[0];
    } else {
        /* all series except "const" */
        nv = dset.v - 1;
    }

    /* we won't treat errors here as fatal */

    if (nv != M.cols) {
        /* "can't happen" */
        return;
    }

    S = strings_array_new(nv);
    if (S == null) {
        return;
    }

    for (i=1; i<=nv; i++) {
        vi = list == null ? i : list[i];
        S[i-1] = gretl_strdup(dset.varname[vi]);
        if (S[i-1] == null) {
            strings_array_free(S, nv);
            return;
        }
    }

    err = gretl_matrix_set_colnames(M, S);
    if (err) {
        strings_array_free(S, nv);
    }
}

GretlMatrixRaw * real_gretl_matrix_data_subset(const(int) * list,
                                            const(DatasetRaw) * dset,
                                            int t1, int t2,
                                            const(char) * mask,
                                            int op, int * err)
{
    GretlMatrixRaw * M = null;
    int T, Tmax = t2 - t1 + 1;
    int do_obs_info = 1;
    int k, mt1 = 0, mt2 = 0;
    int skip;
    int j, vj, s, t;

    if (list != null) {
        k = list[0];
    } else {
        k = dset.v - 1;
    }

    if (t1 < 0 && t2 < 0) {
        /* use full dataset */
        t1 = 0;
        t2 = dset.n - 1;
        T = Tmax = dset.n;
        do_obs_info = 0;
    }

    if (k <= 0 || Tmax <= 0 || t2 >= dset.n) {
        *err = GretlError.data;
        return null;
    }

    *err = 0;
    T = Tmax;

    if (mask != null) {
        T -= get_mask_count(mask, Tmax);
    } else if (op == MMissingCode.trim) {
        *err = list_adjust_sample(list, &t1, &t2, dset, null);
        if (*err) {
            return null;
        } else {
            Tmax = T = t2 - t1 + 1;
        }
    } else if (op == MMissingCode.skip || op == MMissingCode.error) {
        for (t=t1; t<=t2 && !*err; t++) {
            for (j=1; j<=k; j++) {
                vj = list == null ? j : list[j];
                if (na(dset.Z[vj][t])) {
                    if (op == MMissingCode.skip) {
                        T--;
                    } else {
                        *err = GretlError.missdata;
                        return null;
                    }
                    break;
                }
            }
        }
    }

    if (T < 0) {
        /* "can't happen" */
        *err = GretlError.data;
    } else {
        M = gretl_matrix_alloc(T, k);
        if (M == null) {
            *err = GretlError.alloc;
        }
    }

    if (*err || T == 0) {
        return M;
    }

    if (do_obs_info) {
        mt1 = t1; mt2 = t2;
    }

    s = 0;
    for (t=t1; t<=t2; t++) {
        skip = 0;
        if (mask != null) {
            skip = mask[t - t1];
        } else if (op == MMissingCode.skip) {
            for (j=1; j<=k; j++) {
                vj = list == null ? j : list[j];
                if (na(dset.Z[vj][t])) {
                    skip = 1;
                    break;
                }
            }
        }
        if (!skip) {
            for (j=1; j<=k; j++) {
                vj = list == null ? j : list[j];
                gretl_matrix_set(M, s, j-1, dset.Z[vj][t]);
            }
            if (s == 0) {
                mt1 = t;
            } else if (s == T - 1) {
                mt2 = t;
                break;
            }
            s++;
        }
    }

    if (!*err && (mask != null || op == MMissingCode.ok)) {
        int n = k * T;

        for (j=0; j<n && !*err; j++) {
            if (na(M.val[j])) {
                if (op != MMissingCode.ok) {
                    *err = GretlError.missdata;
                }
            }
        }
    }

    if (*err) {
        gretl_matrix_free(M);
        M = null;
    } else {
        if (do_obs_info && T == mt2 - mt1 + 1) {
            gretl_matrix_set_t1(M, mt1);
            gretl_matrix_set_t2(M, mt2);
        }
        if (dset.varname != null) {
            add_dataset_colnames(M, list, dset);
        }
    }

    return M;
}

/**
 * gretl_matrix_data_subset:
 * @list: list of series to process, or %NULL to include
 * all series except "const".
 * @dset: pointer to dataset struct.
 * @t1: starting observation.
 * @t2: ending observation.
 * @missop: how to handle missing observations.
 * @err: location to receive error code.
 *
 * Creates a matrix holding the subset of variables from
 * @dset specified by @list, over the sample range @t1 to @t2,
 * inclusive. Series are in columns. The @missop flag can
 * be %M_MISSING_OK to indicate that it's OK to include
 * missing values in the matrix, %M_MISSING_ERROR (it's an
 * error if any missing values are found), or %M_MISSING_SKIP
 * (observations with any missing values are omitted from the
 * matrix).
 *
 * Returns: allocated matrix or NULL on failure.
 */

GretlMatrixRaw * gretl_matrix_data_subset(const(int) * list,
                                       const(DatasetRaw) * dset,
                                       int t1, int t2,
                                       int missop,
                                       int * err)
{
    return real_gretl_matrix_data_subset(list, dset, t1, t2, null,
                                         missop, err);
}

/**
 * gretl_matrix_data_subset_masked:
 * @list: list of variable to process.
 * @dset: dataset struct.
 * @t1: starting observation.
 * @t2: ending observation.
 * @mask: missing observations mask, or NULL.
 * @err: location to receive error code.
 *
 * Creates a gretl matrix holding the subset of variables from
 * @dset specified by @list, over the sample range @t1 to @t2,
 * inclusive.  Variables are in columns.  @mask should be an
 * array of char of length (@t2 - @t1 + 1) with 1s in the positions
 * of observations to exclude from the subset and zeros elsewhere.
 * This apparatus can be used to exclude missing observations.
 *
 * Returns: allocated matrix or NULL on failure.
 */

GretlMatrixRaw * gretl_matrix_data_subset_masked(const(int) * list,
                                              const(DatasetRaw) * dset,
                                              int t1, int t2,
                                              const(char) * mask,
                                              int * err)
{
    return real_gretl_matrix_data_subset(list, dset, t1, t2, mask,
                                         MMissingCode.error, err);
}

/* count the number of selected rows in the current
   sample range */

int mmask_row_count(const(GretlMatrixRaw) * mask, const(DatasetRaw) * dset)
{
    int t, m = 0;

    for (t=dset.t1; t<=dset.t2; t++) {
        if (mask.val[t] != 0) {
            m++;
        }
    }

    return m;
}

/**
 * gretl_matrix_data_subset_special:
 * @list: list of variables to process.
 * @dset: dataset struct.
 * @mmask: matrix holding desired observations mask.
 * @err: location to receive error code.
 *
 * Creates a gretl matrix holding the subset of variables from
 * @dset specified by @list, using the observations (data rows)
 * selected by non-zero elements in @mmask. This is designed
 * to support the gretl "libset" variable %matrix_mask.
 *
 * The length of the vector @mmask must equal the number of
 * observations in the dataset, the member %n of @dset.
 *
 * Returns: allocated matrix or NULL on failure.
 */

GretlMatrixRaw * gretl_matrix_data_subset_special(const(int) * list,
                                               const(DatasetRaw) * dset,
                                               const(GretlMatrixRaw) * mmask,
                                               int * err)
{
    int n = gretl_vector_get_length(mmask);
    GretlMatrixRaw * X = null;

    if (list == null || n != dset.n) {
        *err = GretlError.data;
    } else if (list[0] == 0) {
        X = gretl_null_matrix_new();
    } else {
        int T = mmask_row_count(mmask, dset);
        int k = list[0];

        if (T == 0) {
            X = gretl_null_matrix_new();
        } else {
            X = gretl_matrix_alloc(T, k);
        }

        if (X != null && T > 0) {
            const(double) * xi;
            int i, s, t;

            for (i=0; i<k; i++) {
                xi = dset.Z[list[i+1]];
                s = 0;
                for (t=dset.t1; t<=dset.t2; t++) {
                    if (mmask.val[t] != 0) {
                        if (s == 0) {
                            gretl_matrix_set_t1(X, t);
                        } else if (s == T - 1) {
                            gretl_matrix_set_t2(X, t);
                        }
                        gretl_matrix_set(X, s++, i, xi[t]);
                    }
                }
            }
        }
    }

    if (X == null && *err == 0) {
        *err = GretlError.alloc;
    }

    return X;
}

void check_matrix_varnames(DatasetRaw * dset)
{
    int i, j, err = 0;

    for (i=1; i<dset.v; i++) {
        for (j=1; j<i; j++) {
            if (!strcmp(dset.varname[i], dset.varname[j])) {
                err = 1;
                break;
            }
        }
    }

    if (err) {
        for (i=1; i<dset.v; i++) {
            sprintf(dset.varname[i], "v%d".ptr, i);
        }
    }
}

/**
 * gretl_dataset_from_matrix:
 * @m: source matrix.
 * @list: list of columns (1-based) to include, or NULL to
 * include all columns.
 * @opt: may include OPT_B to attempt "borrowing" of data;
 * OPT_N to use plain numbers as variable names; OPT_R to
 * use row-names as observation markers, if present; OPT_S
 * when called from the "store" command.
 * @err: location to receive error code.
 *
 * Creates a gretl dataset from matrix @m, either using the
 * columns specified in @list or using all columns if @list
 * is NULL.
 *
 * Returns: pointer to new dataset information struct on success,
 * or NULL on failure.
 */

DatasetRaw * gretl_dataset_from_matrix(const(GretlMatrixRaw) * m,
                                    const(int) * list,
                                    Gretlopt opt,
                                    int * err)
{
    DatasetRaw * dset = null;
    const(char) ** cnames = null;
    const(char) ** rnames = null;
    int use_rownames;
    int i, t, col, nv, T;

    if (gretl_is_null_matrix(m)) {
        *err = GretlError.data;
        return null;
    }

    use_rownames = (opt & (Gretlopt.OPT_R | Gretlopt.OPT_S)) != 0;

    T = m.rows;
    nv = m.cols;

    if (list != null) {
        for (i=1; i<=list[0]; i++) {
            col = list[i];
            if (col < 1 || col > nv) {
                gretl_errmsg_sprintf("Variable number %d is out of bounds".ptr, col);
                *err = GretlError.data;
                break;
            }
        }
        nv = list[0];
    }

    if (!*err) {
        Gretlopt dsopt = Gretlopt.OPT_NONE;

        if (opt & Gretlopt.OPT_B) {
            dsopt |= Gretlopt.OPT_B;
        }
        if (use_rownames) {
            rnames = gretl_matrix_get_rownames(m);
            if (rnames != null) {
                dsopt |= Gretlopt.OPT_M;
            }
        }
        dset = create_auxiliary_dataset(nv + 1, T, dsopt);
        if (dset == null) {
            *err = GretlError.alloc;
        }
    }

    if (*err) {
        return null;
    }

    cnames = gretl_matrix_get_colnames(m);

    for (i=1; i<=nv; i++) {
        double * src;

        col = (list != null)? list[i] - 1 : i - 1;
        src = cast(double *) m.val + T * col;
        if (opt & Gretlopt.OPT_B) {
            /* "borrowing" */
            dset.Z[i] = src;
        } else {
            /* copying */
            memcpy(dset.Z[i], src, T * (*src).sizeof);
        }
        if (cnames != null) {
            if (gretl_normalize_varname(dset.varname[i], cnames[col], 0, i)) {
                series_record_display_name(dset, i, cnames[col]);
            }
        } else if (opt & Gretlopt.OPT_N) {
            sprintf(dset.varname[i], "%d".ptr, col + 1);
        } else if (opt & Gretlopt.OPT_R) {
            sprintf(dset.varname[i], "v%d".ptr, i);
        } else {
            sprintf(dset.varname[i], "col%d".ptr, col + 1);
        }
    }

    if (rnames != null) {
        for (t=0; t<T; t++) {
            dset.S[t][0] = '\0';
            strncat(dset.S[t], rnames[t], OBSLEN-1);
        }
    }

    if (cnames != null && (opt & Gretlopt.OPT_S)) {
        /* we must have valid varnames for "store" */
        check_matrix_varnames(dset);
    }

    return dset;
}

/**
 * gretl_dataset_from_matrices:
 * @m1: first source matrix.
 * @m2: second source matrix.
 * @err: location to receive error code.
 *
 * Creates a gretl dataset from the side-by-side concatenation
 * of @m1 and @m2. Storage for the values is "borrowed" from the
 * columns of the matrices. Column names are used to label the
 * series, if present; otherwise the labels are "v1", "v2" and so on.
 *
 * The two matrices must have the same number of rows. Rownames on
 * either matrix are disregarded.
 *
 * Returns: pointer to new dataset on success, or NULL on failure.
 */

DatasetRaw * gretl_dataset_from_matrices(const(GretlMatrixRaw) * m1,
                                      const(GretlMatrixRaw) * m2,
                                      int * err)
{
    DatasetRaw * dset = null;
    const(char) ** cnames1 = null;
    const(char) ** cnames2 = null;
    const(char) ** cnames = null;
    const(GretlMatrixRaw) * m;
    int i, col, nv, T;

    if (gretl_is_null_matrix(m1) ||
        gretl_is_null_matrix(m2)) {
        *err = GretlError.data;
        return null;
    } else if (m1.rows != m2.rows) {
        *err = GretlError.nonconf;
        return null;
    }

    T = m1.rows;
    nv = m1.cols + m2.cols;

    dset = create_auxiliary_dataset(nv + 1, T, Gretlopt.OPT_B);
    if (dset == null) {
        *err = GretlError.alloc;
        return null;
    }

    cnames1 = gretl_matrix_get_colnames(m1);
    cnames2 = gretl_matrix_get_colnames(m2);
    m = m1;
    cnames = cnames1;
    col = 0;

    for (i=1; i<=nv; i++) {
        dset.Z[i] = cast(double *) m.val + T * col;
        if (cnames != null) {
            if (gretl_normalize_varname(dset.varname[i], cnames[col], 0, i)) {
                series_record_display_name(dset, i, cnames[col]);
            }
        } else {
            sprintf(dset.varname[i], "v%d".ptr, i);
        }
        if (m == m1 && i == m1.cols) {
            m = m2;
            cnames = cnames2;
            col = 0;
        } else {
            col++;
        }
    }

    return dset;
}

int write_matrix_as_dataset(const(char) * fname, Gretlopt opt, Prn * prn)
{
    const(char) * mname;
    GretlMatrixRaw * m;
    DatasetRaw * mdset;
    int err = 0;

    mname = get_optval_string(GretlCmdIndex.store, Gretlopt.OPT_A);
    m = get_matrix_by_name(mname);
    if (m == null) {
        return GretlError.data;
    }

    mdset = gretl_dataset_from_matrix(m, null, Gretlopt.OPT_S, &err);
    if (!err) {
        opt = cast(Gretlopt)(opt & ~Gretlopt.OPT_A);
        err = write_data(fname, null, mdset, opt, prn);
    }

    destroy_dataset(mdset);

    return err;
}

/**
 * gretl_plotfit_matrices:
 * @yvar: the y variable.
 * @xvar: the x variable.
 * @fit: type of fit sought.
 * @t1: starting observation.
 * @t2: ending observation.
 * @py: location to receive y vector.
 * @pX: location to receive X matrix.
 *
 * Creates a vector y and matrix X based on the input @yvar,
 * @xvar and @fit, using the given sample range.  An observation
 * is skipped if either @yvar or @xvar is missing at that
 * observation.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_plotfit_matrices(const(double) * yvar, const(double) * xvar,
                           FitType fit, int t1, int t2,
                           GretlMatrixRaw ** py, GretlMatrixRaw ** pX)
{
    GretlMatrixRaw * y = null;
    GretlMatrixRaw * X = null;
    char * mask = null;
    double xt;
    int T = t2 - t1 + 1;
    int n = 0;
    int i, j, k, s, t;
    int err = 0;

    if (T <= 0) {
        return GretlError.data;
    }

    if (fit == FitType.loglin && !gretl_ispositive(t1, t2, yvar, 1)) {
        gretl_errmsg_set(_("Non-positive values encountered").ptr);
        return GretlError.data;
    } else if (fit == FitType.linlog && !gretl_ispositive(t1, t2, xvar, 1)) {
        gretl_errmsg_set(_("Non-positive values encountered").ptr);
        return GretlError.data;
    }

    mask = cast(char *) calloc(T, 1);
    if (mask == null) {
        return GretlError.alloc;
    }

    for (s=0; s<T; s++) {
        t = s + t1;
        if (na(yvar[t]) || (xvar != null && na(xvar[t]))) {
            mask[s] = 1;
        } else {
            n++;
        }
    }

    if (n == 0) {
        free(mask);
        return GretlError.missdata;
    }

    if (fit == FitType.cubic) {
        k = 4;
    } else if (fit == FitType.quadratic) {
        k = 3;
    } else if (fit == FitType.loess) {
        k = 1;
    } else {
        k = 2;
    }

    y = gretl_column_vector_alloc(n);
    X = gretl_matrix_alloc(n, k);
    if (y == null || X == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    i = 0;
    for (s=0; s<T; s++) {
        t = s + t1;
        if (!mask[s]) {
            j = 0;
            if (fit == FitType.loglin) {
                y.val[i] = log(yvar[t]);
            } else {
                y.val[i] = yvar[t];
            }
            if (fit != FitType.loess) {
                gretl_matrix_set(X, i, j++, 1.0);
            }
            xt = (xvar != null)? xvar[t] : s;
            if (fit == FitType.inverse) {
                gretl_matrix_set(X, i, j++, 1.0 / xt);
            } else if (fit == FitType.linlog) {
                gretl_matrix_set(X, i, j++, log(xt));
            } else {
                gretl_matrix_set(X, i, j++, xt);
            }
            if (fit == FitType.quadratic || fit == FitType.cubic) {
                gretl_matrix_set(X, i, j++, xt * xt);
            }
            if (fit == FitType.cubic) {
                gretl_matrix_set(X, i, j, xt * xt * xt);
            }
            i++;
        }
    }

 bailout:

    free(mask);

    if (err) {
        gretl_matrix_free(y);
        gretl_matrix_free(X);
    } else {
        *py = y;
        *pX = X;
    }

    return err;
}

int skip_matrix_comment(FILE * fp, int * rows, int * cols,
                        int * is_complex, int * err)
{
    int c, ret = 0;

    c = fgetc(fp);

    if (c == '#') {
        char * p;
        char[16] test = 0;
        int n, i = 0;

        ret = 1;
        while (c != '\n' && c != EOF && c != -1) {
            c = fgetc(fp);
            if (i < 15) {
                test[i] = cast(char) c;
            }
            i++;
        }
        if ((p = strstr(test.ptr, "rows:")) != null) {
            if (sscanf(p + 5, "%d", &n) == 1) {
                *rows = n;
            }
        } else if ((p = strstr(test.ptr, "columns:")) != null) {
            if (sscanf(p + 8, "%d", &n) == 1) {
                *cols = n;
            }
        } else if ((p = strstr(test.ptr, "complex:")) != null) {
            if (sscanf(p + 8, "%d", &n) == 1) {
                *is_complex = (n == 1);
            }
        }
    } else {
        ungetc(c, fp);
    }

    if (c == EOF || c == -1) {
        fprintf(stderr, "reached premature end of file\n".ptr);
        *err = GretlError.data;
    }

    return ret;
}

GretlMatrixRaw * read_binary_matrix_file(FILE * fp, int * err)
{
    GretlMatrixRaw * A = null;
    char[20] header = 0;
    int is_complex = 0;
    int[2] dim;

    if (fread(header.ptr, 1, 19, fp) < 19) {
        *err = GretlError.data;
    } else if (fread(dim.ptr, (*dim.ptr).sizeof, 2, fp) < 2) {
        *err = GretlError.data;
    } else if (dim[0] <= 0 || dim[1] <= 0) {
        *err = GretlError.data;
    }

    if (!*err) {
        if (!strcmp(header.ptr, "gretl_binary_matrix")) {
            ; /* OK */
        } else if (!strcmp(header.ptr, "gretl_binar_cmatrix")) {
            is_complex = 1;
        } else {
            *err = GretlError.data;
        }
    }

    if (!*err) {
        A = gretl_matrix_alloc(dim[0], dim[1]);
        if (A == null) {
            *err = GretlError.alloc;
        }
    }

    if (!*err) {
        size_t n = dim[0] * dim[1];

        if (fread(A.val, (*A.val).sizeof, n, fp) < n) {
            *err = GretlError.data;
            gretl_matrix_free(A);
            A = null;
        } else if (is_complex) {
            gretl_matrix_set_complex_full(A, 1);
        }
    }

    fclose(fp);

    return A;
}

/* In reading matrices, accept "NA" and "." (Stata) for NaN? */

double unix_scan_NA(FILE * fp, int * err)
{
    char[3] test = 0;

    if (fscanf(fp, "%2s", test.ptr) < 1) {
        *err = GretlError.data;
        return 0;
    } else if (!strcmp(test.ptr, "NA") || !strcmp(test.ptr, ".")) {
        return NADBL;
    } else {
        gretl_errmsg_sprintf(_("got invalid field '%s'").ptr, test.ptr);
        *err = GretlError.data;
        return 0;
    }
}

/**
 * gretl_matrix_read_from_file:
 * @fname: name of file.
 * @import: non-zero means we're importing via dotdir.
 * @err: location to receive error code.
 *
 * Reads a matrix from a file by the name @fname; if it's
 * a text file the column separator must be space or tab. It is
 * assumed that the dimensions of the matrix (number of rows and
 * columns) are found on the first line of the file, so no
 * heuristics are necessary. In case of error @err is filled
 * appropriately.
 *
 * Returns: The matrix as read from file, or NULL.
 */

GretlMatrixRaw * gretl_matrix_read_from_file(const(char) * fname, int import_, int * err)
{
    char[FILENAME_MAX] fullname;
    char * unzname = null;
    int r = -1, c = -1, n = 0;
    int gz, bin = 0;
    GretlMatrixRaw * A = null;
    int is_complex = 0;
    FILE * fp = null;

    gz = has_suffix(fname, ".gz");

    if (!gz) {
        bin = has_suffix(fname, ".bin");
    }

    if (import_) {
        gretl_build_path(fullname.ptr, gretl_dotdir(), fname, null);
    } else {
        strcpy(fullname.ptr, fname);
        gretl_maybe_prepend_dir(fullname.ptr);
    }

    if (gz) {
        char[FILENAME_MAX] tmp;

        gretl_build_path(tmp.ptr, gretl_dotdir(), "_unztmp_.mat".ptr, null);
        *err = gretl_gunzip(fullname.ptr, tmp.ptr);
        if (!*err) {
            fp = gretl_fopen(tmp.ptr, "rb");
            unzname = g_strdup(tmp.ptr);
        }
    } else {
        fp = gretl_fopen(fullname.ptr, "rb");
    }

    if (fp == null) {
        *err = GretlError.fopen;
        free(unzname);
        return null;
    }

    if (bin) {
        return read_binary_matrix_file(fp, err);
    }

    /* 'Eat' any leading comment lines starting with '#',
       but while we go, scan for Matlab-style rows/columns
       specification.
    */
    while (!*err && skip_matrix_comment(fp, &r, &c, &is_complex, err)) {
        ;
    }

    if (!*err) {
        if (r >= 0 && c >= 0) {
            /* we got dimensions from the preamble */
            n = 2;
        } else {
            n = fscanf(fp, "%d %ds\n".ptr, &r, &c);
        }
        if (n < 2 || r < 0 || c < 0) {
            fprintf(stderr, "error reading rows, cols (r=%d, c=%d)\n".ptr,
                    r, c);
            *err = GretlError.data;
        } else {
            A = gretl_matrix_alloc(r, c);
            if (A == null) {
                *err = GretlError.alloc;
            }
        }
    }

    if (!*err) {
        import core.stdc.config : c_long;
        c_long pos = 0;
        double x;
        int i, j;

        gretl_push_c_numeric_locale();

        for (i=0; i<r && !*err; i++) {
            for (j=0; j<c && !*err; j++) {
                pos = ftell(fp);
                n = fscanf(fp, "%lf".ptr, &x);
                if (n != 1) {
                    fseek(fp, pos, SEEK_SET);
                    x = unix_scan_NA(fp, err);
                }
                if (*err) {
                    fprintf(stderr, "error reading row %d, column %d\n".ptr, i+1, j+1);
                } else {
                    gretl_matrix_set(A, i, j, x);
                }
            }
        }

        gretl_pop_c_numeric_locale();
    }

    fclose(fp);

    if (unzname != null) {
        gretl_remove(unzname);
        free(unzname);
    }

    if (!*err && is_complex) {
        gretl_matrix_set_complex_full(A, 1);
    }

    if (*err && A != null) {
        gretl_matrix_free(A);
        A = null;
    }

    return A;
}

int matrix_to_csv(const(GretlMatrixRaw) * A, const(char) * fname)
{
    const(char) * na_str;
    char delim;
    FILE * fp;
    double aij;
    int i, j;

    fp = gretl_fopen(fname, "wb");
    if (fp == null) {
        return GretlError.fopen;
    }

    na_str = get_csv_na_write_string();
    delim = get_data_export_delimiter();
    if (delim == '\t') {
        /* we can't allow this! */
        delim = ',';
    }

    gretl_push_c_numeric_locale();

    for (i=0; i<A.rows; i++) {
        for (j=0; j<A.cols; j++) {
            aij = gretl_matrix_get(A, i, j);
            if (na(aij)) {
                fputs(na_str, fp);
            } else {
                fprintf(fp, "%.17g".ptr, aij);
            }
            if (j < A.cols - 1) {
                fputc(delim, fp);
            } else {
                fputc('\n', fp);
            }
        }
    }

    gretl_pop_c_numeric_locale();
    fclose(fp);

    return 0;
}

/**
 * gretl_matrix_write_to_file:
 * @A: matrix to write.
 * @fname: name of file.
 * @use_dotdir: non-zero means we're exporting via dotdir.
 *
 * Writes the matrix @A to a file by the name @fname; if it's a text
 * file, the column separator is tab. The number of rows and columns
 * are written on the first line of the file.
 *
 * Returns: 0 on successful completion, non-zero code on error.
 */

int gretl_matrix_write_to_file(GretlMatrixRaw * A, const(char) * fname, int use_dotdir)
{
    char[FILENAME_MAX] targ;
    int r, c, i, j;
    int is_complex = 0;
    Prn * prn = null;
    FILE * fp = null;
    int csv = 0;
    int gz, bin = 0;
    int err = 0;

    gz = has_suffix(fname, ".gz");

    if (!gz) {
        bin = has_suffix(fname, ".bin");
    }

    if (!gz && !bin) {
        csv = has_suffix(fname, ".csv");
        if (csv && A.is_complex) {
            return GretlError.cmplx;
        }
    }

    if (use_dotdir) {
        gretl_build_path(targ.ptr, gretl_dotdir(), fname, null);
    } else {
        fname = gretl_maybe_switch_dir(fname);
        strcpy(targ.ptr, fname);
    }

    if (csv) {
        return matrix_to_csv(A, targ.ptr);
    }

    if (bin) {
        fp = gretl_fopen(targ.ptr, "wb");
    } else if (gz) {
        prn = gretl_gzip_print_new(targ.ptr, -1, &err);
    } else {
        prn = gretl_print_new_with_filename(targ.ptr, &err);
    }

    if (fp == null && prn == null) {
        return GretlError.fopen;
    }

    if (A.is_complex) {
        /* reset to "pretend real" for writing */
        is_complex = 1;
        gretl_matrix_set_complex_full(A, 0);
    }

    r = A.rows;
    c = A.cols;

    if (bin) {
        const(char) * header = is_complex ? "gretl_binar_cmatrix" :
            "gretl_binary_matrix";
        int[2] dim = [r, c];
        size_t n = r * c;

        fwrite(header, 1, strlen(header), fp);
        fwrite(dim.ptr, (*dim.ptr).sizeof, 2, fp);
        fwrite(A.val, (*A.val).sizeof, n, fp);
        fclose(fp);
    } else {
        /* !bin: textual representation */
        int format_g = libset_get_bool(SetKey.mwrite_g);
        char pad, d = '\t';
        double x;

        if (is_complex) {
            pprintf(prn, "# rows: %d\n".ptr, r);
            pprintf(prn, "# columns: %d\n".ptr, c);
            pprintf(prn, "# complex: %d\n".ptr, 1);
        } else {
            pprintf(prn, "%d%c%d\n".ptr, r, d, c);
        }

        gretl_push_c_numeric_locale();

        for (i=0; i<r; i++) {
            for (j=0; j<c; j++) {
                pad = (j == c-1)? '\n' : d;
                x = gretl_matrix_get(A, i, j);
                if (format_g) {
                    pprintf(prn, "%g".ptr, x);
                } else {
                    pprintf(prn, "%26.18E".ptr, x);
                }
                pputc(prn, pad);
            }
        }

        gretl_pop_c_numeric_locale();
        gretl_print_destroy(prn);
    }

    if (is_complex) {
        /* reset A's original status */
        gretl_matrix_set_complex_full(A, 1);
    }

    return err;
}

void make_numstr(char * s, double x)
{
    if (x == -0.0) {
        x = 0.0;
    }

    if (x != x) {
        strcpy(s, "nan");
        return;
    }

    sprintf(s, "%#.5g".ptr, x);

    /* remove surplus, or add deficient, zeros */

    if (strstr(s, ".00000")) {
        s[strlen(s) - 1] = 0;
    } else {
        char * p = s;
        int n = 0;

        while (*p) {
            if (isdigit(*p)) {
                n++;
            } else if (isalpha(*p)) {
                n = 0;
                break;
            }
            p++;
        }
        if (n > 0 && n < 5) {
            strcat(s, "0");
        }
    }
}

int max_numchars(const(GretlMatrixRaw) * m)
{
    char[24] s;
    int i, n = m.rows * m.cols;
    int c, cmax = 0;

    for (i=0; i<n && cmax<6; i++) {
        sprintf(s.ptr, "%g".ptr, m.val[i]);
        c = cast(int) strlen(s.ptr);
        if (c > cmax) {
            cmax = c;
        }
    }

    return cmax;
}

int max_label_length(const(char) ** names, int n)
{
    int i, len, maxlen = 0;

    for (i=0; i<n; i++) {
        if (names[i] != null) {
            len = cast(int) g_utf8_strlen(names[i], -1);
            if (len > maxlen) {
                maxlen = len;
            }
        }
    }

    return maxlen;
}

void real_matrix_print_to_prn(const(GretlMatrixRaw) * m,
                              const(char) * msg,
                              int plain,
                              const(char) ** colnames,
                              const(DatasetRaw) * dset,
                              int imin, int imax,
                              Prn * prn)
{
    const(char) ** rownames = null;
    char[32] numstr;
    char[OBSLEN] obs;
    double x;
    int strwidth = 12;
    int rnamelen = 0;
    int dated = 0;
    int cmax = 0, cpad = 2;
    int mt1 = 0, mt2 = 0;
    int i, j;

    if (prn == null) {
        return;
    }

    if (m == null) {
        if (msg != null && *msg != '\0') {
            pprintf(prn, _("%s: matrix is NULL").ptr, msg);
        } else {
            pputs(prn, _("matrix is NULL").ptr);
        }
        pputc(prn, '\n');
        return;
    } else if (m.rows == 0 || m.cols == 0) {
        if (msg != null && *msg != '\0') {
            pprintf(prn, _("%s: matrix is empty (%d x %d)").ptr,
                    msg, m.rows, m.cols);
        } else {
            pprintf(prn, _("matrix is empty (%d x %d)").ptr,
                    m.rows, m.cols);
        }
        pputc(prn, '\n');
        return;
    }

    dated = gretl_matrix_is_dated(m);
    if (dated) {
        mt1 = gretl_matrix_get_t1(m);
        mt2 = gretl_matrix_get_t2(m);
    }

    /* @plain means skip the header stuff */

    if (msg != null && *msg != '\0' && !plain) {
        pprintf(prn, "%s (%d x %d)".ptr, msg, m.rows, m.cols);
        if (imin > 0 && imax > 0) {
            pprintf(prn, " rows %d to %d\n\n".ptr, imin+1, imax);
        } else if (dated && dset == null) {
            pprintf(prn, " [t1 = %d, t2 = %d]\n\n".ptr, mt1 + 1, mt2 + 1);
        } else {
            pputs(prn, "\n\n");
        }
    }

    cmax = max_numchars(m);
    if (cmax > 5) {
        cmax = 0;
    }

    if (colnames == null) {
        colnames = gretl_matrix_get_colnames(m);
    }
    if (colnames != null) {
        int n, maxlen = 0;

        for (i=0; i<m.cols; i++) {
            n = cast(int) strlen(colnames[i]);
            if (n+1 > maxlen) {
                maxlen = n+1;
            }
        }
        if (maxlen < strwidth) {
            strwidth = maxlen < 9 ? 9 : maxlen;
        }
    }

    rownames = gretl_matrix_get_rownames(m);

    if (rownames != null) {
        rnamelen = max_label_length(rownames, m.rows);
    } else if (dated && dset != null) {
        rnamelen = max_obs_marker_length(dset);
    }

    if (colnames != null) {
        int numeric_names = 1;

        if (rnamelen > 0) {
            bufspace(rnamelen + 1, prn);
        }
        for (j=0; j<m.cols; j++) {
            pprintf(prn, "%*.*s ".ptr, strwidth, strwidth, colnames[j]);
            if (!numeric_string(colnames[j])) {
                numeric_names = 0;
            }
        }
        pputc(prn, '\n');
        if (numeric_names) {
            pputc(prn, '\n');
        }
    }

    if (imin < 0) imin = 0;
    if (imax < 0) imax = m.rows;

    for (i=imin; i<imax; i++) {
        if (rnamelen > 0) {
            if (rownames != null) {
                pprintf(prn, "%*s ".ptr, rnamelen, rownames[i]);
            } else {
                ntolabel(obs.ptr, i + mt1, dset);
                pprintf(prn, "%*s ".ptr, rnamelen, obs.ptr);
            }
        }
        for (j=0; j<m.cols; j++) {
            x = gretl_matrix_get(m, i, j);
            if (cmax > 0) {
                if (colnames != null) {
                    bufspace(strwidth - cmax - cpad, prn);
                }
                if (x != x) {
                    strcpy(numstr.ptr, "nan");
                } else {
                    sprintf(numstr.ptr, "%g".ptr, x);
                }
                pprintf(prn, "%*s ".ptr, cmax + cpad, numstr.ptr);
            } else {
                make_numstr(numstr.ptr, x);
                pprintf(prn, "%*s ".ptr, strwidth, numstr.ptr);
            }
        }
        pputc(prn, '\n');
    }

    pputc(prn, '\n');
}

/**
 * gretl_matrix_print_to_prn:
 * @m: matrix to print.
 * @msg: accompanying message text (or NULL if no message is wanted).
 * @prn: pointer to gretl printing struct.
 *
 * Prints the matrix @m to @prn.
 */

void gretl_matrix_print_to_prn(const(GretlMatrixRaw) * m, const(char) * msg, Prn * prn)
{
    real_matrix_print_to_prn(m, msg, 0, null, null,
                             -1, -1, prn);
}

/**
 * gretl_matrix_print_range:
 * @m: matrix to print.
 * @rmin: first row to print.
 * @prn: last row to print.
 *
 * Prints the specified rows of matrix @m to @prn.
 */

void gretl_matrix_print_range(const(GretlMatrixRaw) * m,
                              const(char) * msg,
                              int rmin, int rmax,
                              Prn * prn)
{
    if (m.is_complex) {
        gretl_cmatrix_print_range(m, msg, rmin, rmax, prn);
    } else {
        real_matrix_print_to_prn(m, msg, 0, null, null,
                                 rmin, rmax, prn);
    }
}

/**
 * gretl_matrix_print_with_col_heads:
 * @m: T x k matrix to print.
 * @title: accompanying title (or NULL if no title is wanted).
 * @heads: array of k strings to identify the columns.
 * @prn: pointer to gretl printing struct.
 *
 * Prints the matrix @m to @prn, with column headings given
 * by @heads.
 */

void gretl_matrix_print_with_col_heads(const(GretlMatrixRaw) * m,
                                       const(char) * title,
                                       const(char) ** heads,
                                       const(DatasetRaw) * dset,
                                       Prn * prn)
{
    real_matrix_print_to_prn(m, title, 0, heads, dset,
                             -1, -1, prn);
}

void maybe_print_col_heads(const(GretlMatrixRaw) * m,
                           const(char) * fmt,
                           int wid, int prec,
                           int icast, int llen,
                           Prn * prn)
{
    const(char) ** heads;

    heads = gretl_matrix_get_colnames(m);

    if (heads != null) {
        int numeric = 1;
        const(char) * s;
        char[32] wtest;
        double x;
        int j, n, d;

        x = gretl_matrix_get(m, 0, 0);

        if (icast) {
            if (wid >= 0 && prec >= 0) {
                snprintf(wtest.ptr, 32, fmt, wid, prec, cast(int) x);
            } else if (wid >= 0 || prec >= 0) {
                n = (wid >= 0)? wid : prec;
                snprintf(wtest.ptr, 32, fmt, n, cast(int) x);
            } else {
                snprintf(wtest.ptr, 32, fmt, cast(int) x);
            }
        } else {
            if (wid >= 0 && prec >= 0) {
                snprintf(wtest.ptr, 32, fmt, wid, prec, x);
            } else if (wid >= 0 || prec >= 0) {
                n = (wid >= 0)? wid : prec;
                snprintf(wtest.ptr, 32, fmt, n, x);
            } else {
                snprintf(wtest.ptr, 32, fmt, x);
            }
        }

        if (llen > 0) {
            bufspace(llen + 1, prn);
        }

        n = cast(int) strlen(wtest.ptr);
        for (j=0; j<m.cols; j++) {
            s = heads[j];
            d = cast(int) (strlen(s) - g_utf8_strlen(s, -1));
            pprintf(prn, "%*s".ptr, n + d, s);
            if (!numeric_string(s)) {
                numeric = 0;
            }
        }
        pputc(prn, '\n');
        if (numeric) {
            pputc(prn, '\n');
        }
    }
}

/* Return an array of char holding 1s for columns in
   @m that are all integers, 0s otherwise. Except
   that if we find that either none or all the columns
   if @m are composed of integers we scratch the
   array and signal the result via the @intcast pointer.
*/

char * get_intcols(const(GretlMatrixRaw) * m, int * intcast)
{
    char * ret = cast(char *) malloc(m.cols);

    if (ret != null) {
        double x;
        int i, j, n = 0;

        for (j=0; j<m.cols; j++) {
            ret[j] = 1;
            for (i=0; i<m.rows; i++) {
                x = gretl_matrix_get(m, i, j);
                if (x != floor(x)) {
                    ret[j] = 0;
                    break;
                }
            }
            n += ret[j];
        }
        if (n == 0 || n == m.cols) {
            /* no, or all, integer columns */
            free(ret);
            ret = null;
            *intcast = n > 0;
        }
    }

    return ret;
}

/* Take a format which is designed for floating-point
   values and convert it for use with integers: preserve
   the width but discard the precision, and change the
   trailing 'f', 'g' or whatever to 'd'.
*/

void revise_integer_format(char * ifmt)
{
    char * p = strchr(ifmt, '.');

    if (p != null) {
        *p = '\0';
        strcat(p, "d");
    } else {
        ifmt[strlen(ifmt)-1] = 'd';
    }
}

/**
 * gretl_matrix_print_with_format:
 * @m: matrix to print.
 * @fmt: a printf-type format string.
 * @wid: an integer width, or 0.
 * @prec: an integer precision, or 0.
 * @prn: pointer to gretl printing struct.
 *
 * Prints the matrix @m to @prn, with the elements of @m
 * formatted as specified by @fmt, @wid and @prec.
 * The arguments @wid and/or @prec are required only if
 * @fmt contains one or more placeholders ("*") to be filled out.
 * For example, if @fmt is "\%14.6f" then neither @wid nor
 * @prec is needed and both should be set to 0; if
 * @fmt is "\%*.6f" then a @wid value is needed but @prec
 * should be 0.
 */

void gretl_matrix_print_with_format(const(GretlMatrixRaw) * m,
                                    const(char) * fmt,
                                    int wid, int prec,
                                    Prn * prn)
{
    if (prn == null) {
        return;
    }

    if (gretl_is_null_matrix(m) || fmt == null || *fmt == '\0') {
        real_matrix_print_to_prn(m, null, 1, null, null, -1, -1, prn);
    } else if (m.is_complex) {
        gretl_cmatrix_printf(m, fmt, prn);
    } else {
        const(char) ** rownames = null;
        char * ifmt = null;
        char * xfmt = null;
        char * intcols = null;
        int llen = 0, intcast = 0;
        int cpos = cast(int) strlen(fmt) - 1;
        int variant = 0;
        char c = fmt[cpos];
        int i, j, k;
        double x;

        xfmt = gretl_strdup(fmt);

        if (c == 'd' || c == 'u' || c == 'x' || c == 'l') {
            intcast = 1;
        } else if (fmt[cpos-1] == 'v') {
            /* new-style "variant" column format */
            intcols = get_intcols(m, &intcast);
            gretl_delchar('v', xfmt);
            variant = 1;
        } else if (c == 'v') {
            /* old-style "variant" column format */
            intcols = get_intcols(m, &intcast);
            xfmt[cpos] = 'g';
            variant = 1;
        }

        if (intcast || intcols) {
            ifmt = gretl_strdup(xfmt);
            if (variant) {
                revise_integer_format(ifmt);
            }
        }

        rownames = gretl_matrix_get_rownames(m);
        if (rownames != null) {
            llen = max_label_length(rownames, m.rows);
        }

        maybe_print_col_heads(m, xfmt, wid, prec, intcast, llen, prn);

        for (i=0; i<m.rows; i++) {
            if (rownames != null) {
                pprintf(prn, "%*s ".ptr, llen, rownames[i]);
            }
            for (j=0; j<m.cols; j++) {
                x = gretl_matrix_get(m, i, j);
                if (intcols != null) {
                    intcast = intcols[j];
                }
                if (intcast) {
                    if (wid > 0 && prec > 0) {
                        pprintf(prn, ifmt, wid, prec, cast(int) x);
                    } else if (wid > 0 || prec > 0) {
                        k = (wid > 0)? wid : prec;
                        pprintf(prn, ifmt, k, cast(int) x);
                    } else {
                        pprintf(prn, ifmt, cast(int) x);
                    }
                } else {
                    if (wid > 0 && prec > 0) {
                        pprintf(prn, xfmt, wid, prec, x);
                    } else if (wid > 0 || prec > 0) {
                        k = (wid > 0)? wid : prec;
                        pprintf(prn, xfmt, k, x);
                    } else {
                        pprintf(prn, xfmt, x);
                    }
                }
            }
            pputc(prn, '\n');
        }

        free(intcols);
        free(ifmt);
        free(xfmt);
    }
}

/* gretl_matrix_write_constructor() is intended for internal use in
   gretl_bundle.c, for creating a string representation of a matrix
   contained in a bundle when the bundle is passed as argument to
   gretl_bundle_write_constructor(). In that context it will not be
   appropriate to show the constructor for a large matrix in extenso,
   hence the size limitation here. However it's not clear what exactly
   we should do with an oversize matrix.
*/

char * gretl_matrix_write_constructor(const(GretlMatrixRaw) * m)
{
    int n;

    if (gretl_is_null_matrix(m)) {
        return g_strdup("{}");
    }

    n = m.rows * m.cols;

    if (n > 16) {
        /* should we do something else here? */
        return g_strdup("{}");
    } else {
        auto gs = appender!string();
        double mij;
        int i, j;

        gs.put('{');
        for (i=0; i<m.rows; i++) {
            for (j=0; j<m.cols; j++) {
                mij = gretl_matrix_get(m, i, j);
                if (na(mij)) {
                    gs.put("NA");
                } else {
                    gs.put(format("%g", mij));
                }
                if (j < m.cols - 1) {
                    gs.put(',');
                } else if (i < m.rows - 1) {
                    gs.put(';');
                }
            }
        }
        gs.put('}');
        gs.put('\0');
        return g_strdup(gs.data.ptr);
    }
}

/**
 * debug_print_matrix:
 * @m: matrix to print.
 * @msg: accompanying message text (or NULL if no message is wanted).
 *
 * Prints the matrix @m to stderr with high precision, along
 * with the address of the matrix struct.
 */

void debug_print_matrix(const(GretlMatrixRaw) * m, const(char) * msg)
{
    char[64] full = 0;

    if (msg != null) {
        strncpy(full.ptr, msg, 32);
        sprintf(full.ptr + strlen(full.ptr), " (%p)".ptr, cast(void *) m);
    } else {
        sprintf(full.ptr, " (%p)".ptr, cast(void *) m);
    }

    if (m != null) {
        int i, n = m.rows * m.cols;
        int d = cast(int) ceil(log10(cast(double) n));

        fprintf(stderr, "%s\n".ptr, full.ptr);
        for (i=0; i<n; i++) {
            fprintf(stderr, "val[%0*d] = % .10E\n".ptr, d, i, m.val[i]);
        }
    } else {
        Prn * prn;
        int err = 0;

        prn = gretl_print_new(PrnType.stderr, &err);
        if (!err) {
            gretl_matrix_print_to_prn(m, full.ptr, prn);
            gretl_print_destroy(prn);
        }
    }
}

int count_unmasked_elements(const(char) * mask, int n)
{
    int i, c = 0;

    for (i=0; i<n; i++) {
        c += (mask[i] == 0);
    }

    return c;
}

/**
 * gretl_matrix_cut_rows:
 * @m: matrix to process.
 * @mask: character array of length equal to the rows of @m,
 * with 1s indicating rows to be cut, 0s for rows to be
 * retained.
 *
 * In-place reduction of @m based on @mask: the masked rows
 * are cut out of @m.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_matrix_cut_rows(GretlMatrixRaw * m, const(char) * mask)
{
    int i, j, k, n;
    double x;

    if (m == null || mask == null) {
        return GretlError.data;
    }

    n = count_unmasked_elements(mask, m.rows);

    for (j=0; j<m.cols; j++) {
        k = 0;
        for (i=0; i<m.rows; i++) {
            if (!mask[i]) {
                x = gretl_matrix_get(m, i, j);
                m.val[j * n + k] = x;
                k++;
            }
        }
    }

    m.rows = n;

    return 0;
}

/**
 * gretl_matrix_cut_cols:
 * @m: matrix to process.
 * @mask: character array of length equal to the cols of @m,
 * with 1s indicating cols to be cut, 0s for cols to be
 * retained.
 *
 * In-place reduction of @m based on @mask: the masked cols
 * are cut out of @m.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_matrix_cut_cols(GretlMatrixRaw * m, const(char) * mask)
{
    int i, j, k, n;
    double x;

    if (m == null || mask == null) {
        return GretlError.data;
    }

    n = count_unmasked_elements(mask, m.cols);

    k = 0;
    for (j=0; j<m.cols; j++) {
        if (!mask[j]) {
            for (i=0; i<m.rows; i++) {
                x = gretl_matrix_get(m, i, j);
                m.val[k++] = x;
            }
        }
    }

    m.cols = n;

    return 0;
}

/**
 * gretl_matrix_cut_rows_cols:
 * @m: square matrix to process.
 * @mask: character array of length equal to the dimension
 * of @m, with 1s indicating rows and columns to be cut, 0s
 * for rows/columns to be retained.
 *
 * In-place reduction of @m based on @mask: the masked rows
 * and columns are cut out of @m. (The data array within @m
 * is not "physically" resized.)
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_matrix_cut_rows_cols(GretlMatrixRaw * m, const(char) * mask)
{
    GretlMatrixRaw * tmp;
    double x;
    int i, j, k, l, n;

    if (m == null || mask == null) {
        return GretlError.data;
    } else if (m.rows != m.cols) {
        return GretlError.nonconf;
    }

    n = count_unmasked_elements(mask, m.rows);
    if (n == 0) {
        gretl_matrix_reuse(m, 0, 0);
        return 0;
    }

    /* create smaller temporary matrix */
    tmp = gretl_matrix_alloc(n, n);
    if (tmp == null) {
        return GretlError.alloc;
    }

    /* copy unmasked values into temp matrix */
    k = 0;
    for (i=0; i<m.rows; i++) {
        if (!mask[i]) {
            l = 0;
            for (j=0; j<m.cols; j++) {
                if (!mask[j]) {
                    x = gretl_matrix_get(m, i, j);
                    gretl_matrix_set(tmp, k, l++, x);
                }
            }
            k++;
        }
    }

    /* redimension the original matrix, copy the values back,
       and free the temp matrix */
    gretl_matrix_reuse(m, n, n);
    gretl_matrix_copy_values(m, tmp);
    gretl_matrix_free(tmp);

    return 0;
}

/**
 * gretl_matrix_zero_row_mask:
 * @m: matrix to process.
 * @err: location to receive error code.
 *
 * Checks matrix @m for rows that are all zero.  If there are
 * any such rows, constructs a mask of length equal to the
 * number of rows in @m, with 1s indicating zero rows, 0s
 * elsewhere.  If there are no such rows, returns NULL.
 *
 * E_ALLOC is written to @err in case a mask should have
 * been constructed but allocation failed.
 *
 * Returns: allocated mask or NULL.
 */

char * gretl_matrix_zero_row_mask(const(GretlMatrixRaw) * m, int * err)
{
    char * mask = null;
    int row0, any0 = 0;
    int i, j;

    mask = cast(char *) calloc(m.rows, 1);
    if (mask == null) {
        *err = GretlError.alloc;
        return null;
    }

    for (i=0; i<m.rows; i++) {
        row0 = 1;
        for (j=0; j<m.cols; j++) {
            if (gretl_matrix_get(m, i, j) != 0.0) {
                row0 = 0;
                break;
            }
        }
        if (row0) {
            mask[i] = 1;
            any0 = 1;
        }
    }

    if (!any0) {
        free(mask);
        mask = null;
    }

    return mask;
}

/**
 * gretl_matrix_zero_col_mask:
 * @m: matrix to process.
 * @err: location to receive error code.
 *
 * Checks matrix @m for columns that are all zero.  If there are
 * any such columns, constructs a mask of length equal to the
 * number of columns in @m, with 1s indicating zero columns, 0s
 * elsewhere.  If there are no such columns, returns NULL.
 *
 * E_ALLOC is written to @err in case a mask should have
 * been constructed but allocation failed.
 *
 * Returns: allocated mask or NULL.
 */

char * gretl_matrix_zero_col_mask(const(GretlMatrixRaw) * m, int * err)
{
    char * mask = null;
    int col0, any0 = 0;
    int i, j;

    mask = cast(char *) calloc(m.cols, 1);
    if (mask == null) {
        *err = GretlError.alloc;
        return null;
    }

    for (j=0; j<m.cols; j++) {
        col0 = 1;
        for (i=0; i<m.rows; i++) {
            if (gretl_matrix_get(m, i, j) != 0.0) {
                col0 = 0;
                break;
            }
        }
        if (col0) {
            mask[j] = 1;
            any0 = 1;
        }
    }

    if (!any0) {
        free(mask);
        mask = null;
    }

    return mask;
}

/**
 * gretl_matrix_zero_diag_mask:
 * @m: matrix to process.
 * @err: location to receive error code.
 *
 * Checks square matrix @m for diagonal elements that are zero.
 * If there are any such, constructs a mask of length equal to
 * the number of rows (and columns) of @m, with 1s indicating
 * the zero diagonal entries. If there are no zero diagonal
 * elements, returns NULL.
 *
 * E_ALLOC is written to @err in case a mask should have
 * been constructed but allocation failed.
 *
 * Returns: allocated mask or NULL.
 */

char * gretl_matrix_zero_diag_mask(const(GretlMatrixRaw) * m, int * err)
{
    char * mask = null;
    int i, trim = 0;

    if (gretl_is_null_matrix(m)) {
        return null;
    }

    if (m.rows != m.cols) {
        *err = GretlError.nonconf;
        return null;
    }

    for (i=0; i<m.rows; i++) {
        if (gretl_matrix_get(m, i, i) == 0.0) {
            trim = 1;
            break;
        }
    }

    if (trim) {
        mask = cast(char *) calloc(m.rows, 1);
        if (mask == null) {
            *err = GretlError.alloc;
        } else {
            for (i=0; i<m.rows; i++) {
                if (gretl_matrix_get(m, i, i) == 0.0) {
                    mask[i] = 1;
                }
            }
        }
    }

    return mask;
}

/**
 * gretl_matrix_rank_mask:
 * @m: matrix to process.
 * @err: location to receive error code.
 *
 * Performs a QR decomposition of matrix @m and uses this
 * to assess the rank of @m.  If @m is not of full rank,
 * constructs a mask of length equal to the numbers of
 * columns in @m, with 1s in positions corresponding
 * to diagonal elements of R that are effectively 0, and
 * 0s elsewhere.  If @m is of full column rank, NULL is
 * returned.
 *
 * E_ALLOC is written to @err in case a mask should have
 * been constructed but allocation failed.
 *
 * Returns: allocated mask or NULL.
 */

char * gretl_matrix_rank_mask(const(GretlMatrixRaw) * m, int * err)
{
    GretlMatrixRaw * Q = null;
    GretlMatrixRaw * R = null;
    char * mask = null;
    double Ri;
    int fullrank = 1;
    int i, n = m.cols;

    Q = gretl_matrix_copy(m);
    if (Q == null) {
        *err = GretlError.alloc;
        return null;
    }

    R = gretl_matrix_alloc(n, n);
    if (R == null) {
        gretl_matrix_free(Q);
        *err = GretlError.alloc;
        return null;
    }

    *err = gretl_matrix_QR_decomp(Q, R);

    if (!*err) {
        mask = cast(char *) calloc(n, 1);
        if (mask == null) {
            *err = GretlError.alloc;
        }
    }

    if (!*err) {
        for (i=0; i<n; i++) {
            Ri = gretl_matrix_get(R, i, i);
            if (fabs(Ri) < R_DIAG_MIN) {
                mask[i] = 1;
                fullrank = 0;
            }
        }
    }

    if (*err || fullrank) {
        free(mask);
        mask = null;
    }

    gretl_matrix_free(Q);
    gretl_matrix_free(R);

    return mask;
}

/**
 * gretl_matrix_mp_ols:
 * @y: dependent variable vector.
 * @X: matrix of independent variables.
 * @b: vector to hold coefficient estimates.
 * @vcv: matrix to hold the covariance matrix of the coefficients,
 * or NULL if this is not needed.
 * @uhat: vector to hold the regression residuals, or NULL if
 * these are not needed.
 * @s2: pointer to receive residual variance, or NULL.  Note:
 * if @s2 is NULL, the "vcv" estimate will be plain (X'X)^{-1}.
 *
 * Computes OLS estimates using Cholesky factorization, via
 * the GMP multiple-precision library, and puts the
 * coefficient estimates in @b.  Optionally, calculates the
 * covariance matrix in @vcv and the residuals in @uhat.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int gretl_matrix_mp_ols(const(GretlVector) * y, const(GretlMatrixRaw) * X,
                        GretlVector * b, GretlMatrixRaw * vcv,
                        GretlVector * uhat, double * s2)
{
    alias MpOls = extern (C) int function(const(GretlVector) *,
                                          const(GretlMatrixRaw) *,
                                          GretlVector *,
                                          GretlMatrixRaw *,
                                          GretlVector *,
                                          double *);
    MpOls matrix_mp_ols = cast(MpOls) get_plugin_function("matrix_mp_ols");

    if (matrix_mp_ols == null) {
        return 1;
    }

    return matrix_mp_ols(y, X, b, vcv, uhat, s2);
}

/* following: quadrature sources based on John Burkhart's GPL'd
   code at http://people.sc.fsu.edu/~jburkardt
*/

double sign(double x)
{
    return x < 0.0 ? -1.0 : 1.0;
}

/* Diagonalize a Jacobi (symmetric tridiagonal) matrix. On entry, @d
   contains the diagonal and @e the sub-diagonal. On exit, @d contains
   quadrature nodes or abscissae and @z the square roots of the
   corresponding weights; @e is used as workspace.

   This is a C adaptation of subroutine IMTQLX, by Sylvan Elhay,
   Jaroslav Kautsky and John Burkardt. See
   http://people.sc.fsu.edu/~jburkardt/f_src/quadrature_weights/
   specifically, qw_golub_welsch.f90
*/

int diag_jacobi(int n, double * d, double * e, double * z)
{
    double b, c, f, g;
    double p, r, s;
    int j, k, l, m = 0;
    int maxiter = 30;
    int i, ii, mml;
    int err = 0;

    if (n == 1) {
        return 0;
    }

    e[n-1] = 0.0;

    errno = 0;

    for (l=1; l<=n; l++) {
        j = 0;
        for ( ; ; ) {
            for (m=l; m<=n; m++) {
                if (m == n) {
                    break;
                }
                p = fabs(d[m-1]) + fabs(d[m]);
                if (fabs(e[m-1]) <= double.epsilon * p) {
                    break;
                }
            }
            p = d[l-1];
            if (m == l) {
                break;
            }
            if (j >= maxiter) {
                fprintf(stderr, "diag_jacobi: iteration limit exceeded\n".ptr);
                return GretlError.noconv;
            }
            j++;
            g = (d[l] - p) / (2.0 * e[l-1]);
            r =  sqrt(g * g + 1.0);
            g = d[m-1] - p + e[l-1] / (g + fabs(r) * sign(g));
            s = 1.0;
            c = 1.0;
            p = 0.0;
            mml = m - l;

            for (ii=1; ii<=mml; ii++) {
                i = m - ii;
                f = s * e[i-1];
                b = c * e[i-1];
                if (fabs(g) <= fabs(f)) {
                    c = g / f;
                    r =  sqrt(c * c + 1.0);
                    e[i] = f * r;
                    s = 1.0 / r;
                    c = c * s;
                } else {
                    s = f / g;
                    r =  sqrt(s * s + 1.0);
                    e[i] = g * r;
                    c = 1.0 / r;
                    s = s * c;
                }
                g = d[i] - p;
                r = (d[i-1] - g) * s + 2.0 * c * b;
                p = s * r;
                d[i] = g + p;
                g = c * r - b;
                f = z[i];
                z[i] = s * z[i-1] + c * f;
                z[i-1] = c * z[i-1] - s * f;
            }

            d[l-1] = d[l-1] - p;
            e[l-1] = g;
            e[m-1] = 0.0;
        }
    }

    /* sorting */

    for (ii=2; ii<=m; ii++) {
        i = ii - 1;
        k = i;
        p = d[i-1];
        for (j=ii; j<=n; j++) {
            if (d[j-1] < p) {
                k = j;
                p = d[j-1];
            }
        }
        if (k != i) {
            d[k-1] = d[i-1];
            d[i-1] = p;
            p = z[i-1];
            z[i-1] = z[k-1];
            z[k-1] = p;
        }
    }

    if (errno) {
        /* some calculation went badly */
        err = GretlError.nan;
        errno = 0;
    }

    return err;
}

/* Construct the Jacobi matrix for a quadrature rule: @d
   holds the diagonal, @e the subdiagonal.
*/

double make_jacobi(int n, int kind, double * d, double * e)
{
    double z0 = 0.0;
    int i;

    if (kind == QuadMethod.legendre) {
        double xi, xj;

        z0 = 2.0;
        for (i=1; i<=n; i++) {
            xi = i;
            xj = 2.0 * i;
            e[i-1] = sqrt(xi * xi / (xj * xj - 1.0));
        }
    } else if (kind == QuadMethod.laguerre) {
        z0 = 1.0; /* tgamma(1.0) */
        for (i=1; i<=n; i++) {
            d[i-1] = 2.0 * i - 1.0;
            e[i-1] = i;
        }
    } else if (kind == QuadMethod.ghermite) {
        z0 = tgamma(0.5);
        for (i=1; i<=n; i++) {
            e[i-1] = sqrt(i / 2.0);
        }
    }

    return z0;
}

/* compute basic Gauss quadrature formula */

int gauss_quad_basic(int n, int kind, double * x, double * w)
{
    double z0;
    double * tmp;
    int i, err;

    tmp = cast(double *) malloc(n * (*tmp).sizeof);
    if (tmp == null) {
        return GretlError.alloc;
    }

    z0 = make_jacobi(n, kind, x, tmp);
    w[0] = sqrt(z0);
    err = diag_jacobi(n, x, tmp, w);

    if (!err) {
        for (i=0; i<n; i++) {
            w[i] = w[i] * w[i];
        }
    }

    free(tmp);

    return err;
}

/* Scale Legendre quadrature rule to the interval (a, b) */

int legendre_scale(int n, double * x, double * w,
                   double a, double b)
{
    double shft, slp;
    int i;

    if (na(a) || na(b)) {
        return GretlError.invarg;
    }

    if (fabs(b - a) <= double.epsilon) {
        fprintf(stderr, "legendre: |b - a| too small\n".ptr);
        return GretlError.data;
    }

    shft = (a + b) / 2.0;
    slp = (b - a) / 2.0;

    for (i=0; i<n; i++) {
        x[i] = shft + slp * x[i];
        w[i] = w[i] * slp;
    }

    return 0;
}

int hermite_scale(int n, double * x, double * w,
                  double mu, double sigma)
{
    double rtpi = sqrt(PI);
    int i;

    if (na(mu) || na(sigma) || sigma <= 0.0) {
        return GretlError.invarg;
    }

    for (i=0; i<n; i++) {
        x[i] = mu + SQRT2 * sigma * x[i];
        w[i] /= rtpi;
    }

    return 0;
}

/**
 * gretl_quadrule_matrix_new:
 * @n: the order (i.e. the number of abscissae and weights).
 * @method: should be one of the values in #QuadMethod.
 * @a: optional parameter (see below).
 * @b: optional parameter (see below).
 * @err: location to receive error code.
 *
 * Calculates a quadrature "rule" (i.e. a set of abscissae or
 * nodes and associated weights) for use in numerical integration.
 * The three supported methods are Gauss-Hermite, Gauss-Legendre
 * and Gauss-Laguerre.
 *
 * If the optional parameters are not to be used, both should
 * be given the value NADBL.
 *
 * If the method is QUAD_GHERMITE, the default weights are W(x)
 * = exp(-x^2), but @a and @b can be used to apply a normal
 * scaling with mean @a and standard deviation @b. The limits
 * of integration are in principle minus infinity and plus
 * infinity.
 *
 * If the method is QUAD_LEGENDRE the weights are W(x) = 1 and
 * the default limits of integration are -1 and 1, but @a and @b
 * can be used to adjust the lower and upper limit respectively.
 *
 * In the QUAD_LAGUERRE case W(x) = exp(-x) and the limits
 * of integration are 0 and plus infinity; @a and @b are
 * ignored.
 *
 * Returns: an @n x 2 matrix containing the abscissae in the first
 * column and the weights in the second, or NULL on failure.
 */

GretlMatrixRaw * gretl_quadrule_matrix_new(int n, int method,
                                        double a, double b,
                                        int * err)
{
    GretlMatrixRaw * m;

    if (method < QuadMethod.ghermite || method >= QuadMethod.invalid) {
        *err = GretlError.data;
        return null;
    }

    if (n < 0) {
        *err = GretlError.data;
        return null;
    } else if (n == 0) {
        return gretl_null_matrix_new();
    }

    m = gretl_zero_matrix_new(n, 2);

    if (m == null) {
        *err = GretlError.alloc;
    } else {
        double * x = m.val;
        double * w = x + n;

        *err = gauss_quad_basic(n, method, x, w);

        if (!*err) {
            if (method == QuadMethod.legendre) {
                if (na(a) && na(b)) {
                    ; /* default */
                } else if (a == -1.0 && b == 1.0) {
                    ; /* default */
                } else {
                    *err = legendre_scale(n, x, w, a, b);
                }
            } else if (method == QuadMethod.ghermite) {
                if (!na(a) || !na(b)) {
                    *err = hermite_scale(n, x, w, a, b);
                }
            } else if (method == QuadMethod.laguerre) {
                ; /* a and b are ignored */
            }
        }
    }

    if (*err && m != null) {
        gretl_matrix_free(m);
        m = null;
    }

    return m;
}

GretlMatrixRaw * gretl_gauss_hermite_matrix_new(int n, int * err)
{
    return gretl_quadrule_matrix_new(n, QuadMethod.ghermite, 0, 1, err);
}

/**
 * gretl_matrix_2d_convolution:
 * @A: first matrix.
 * @B: second matrix.
 * @err: location to receive error code.
 *
 * Computes the 2-dimensional convolution of the matrices
 * @A and @B. The code borrows from Octave's conv2.m.
 *
 * Returns: newly allocated convolution matrix, or %NULL on
 * failure.
 */

GretlMatrixRaw * gretl_matrix_2d_convolution(const(GretlMatrixRaw) * A,
                                          const(GretlMatrixRaw) * B,
                                          int * err)
{
    GretlMatrixRaw * C;
    double sum;
    int ci, cj, bj, aj, bi, ai;
    int Am = A.rows;
    int An = A.cols;
    int Bm = B.rows;
    int Bn = B.cols;
    int Cm = Am + Bm - 1;
    int Cn = An + Bn - 1;
    int edgM = Bm - 1;
    int edgN = Bn - 1;

    C = gretl_matrix_alloc(Cm, Cn);

    for (ci=0; ci < Cm; ci++) {
        for (cj=0; cj < Cn; cj++) {
            sum = 0;
            for (bj = Bn - 1 - max(0, edgN-cj), aj = max(0, cj-edgN);
                 bj >= 0 && aj < An; bj--, aj++) {
                bi = Bm - 1 - max(0, edgM-ci);
                ai = max(0, ci-edgM);
                for ( ; bi >= 0 && ai < Am; bi--, ai++) {
                    sum += gretl_matrix_get(A, ai, aj) * gretl_matrix_get(B, bi, bj);
                }
            }
            gretl_matrix_set(C, ci, cj, sum);
        }
    }

    return C;
}

GretlMatrixRaw * gretl_matrix_which(const(GretlMatrixRaw) * cond, int * err)
{
    GretlMatrixRaw * ret = null;
    double cij;
    int r = cond.rows;
    int c = cond.cols;
    int n = r * c;
    int m = 0;
    int i, j, k;

    *err = 0;

    for (i=0; i<n; i++) {
        m += cond.val[i] != 0;
    }
    if (m == 0) {
        return gretl_null_matrix_new();
    }

    if (min(r, c) > 1) {
        ret = gretl_matrix_alloc(m, 2);
    } else {
        ret = gretl_matrix_alloc(m, 1);
    }
    if (ret == null) {
        *err = GretlError.alloc;
        return null;
    }

    if (min(r, c) > 1) {
        k = 0;
        for (j=0; j<c; j++) {
            for (i=0; i<r; i++) {
                cij = gretl_matrix_get(cond, i, j);
                if (cij != 0) {
                    cij = r == 1 ? j+1 : i+1;
                    gretl_matrix_set(ret, k, 0, cij);
                    if (c > 1) {
                        gretl_matrix_set(ret, k, 1, j+1);
                    }
                    k++;
                    if (k == m) {
                        goto done;
                    }
                }
            }
        }
    } else {
        k = 0;
        for (i=0; i<max(r, c); i++) {
            if (cond.val[i] != 0) {
                ret.val[k++] = i+1;
                if (k == m) {
                    goto done;
                }
            }
        }
    }

 done:

    return ret;
}

/* scan format string: should be "%<num>m" or just "%m" */

int check_matrix_scan_format(const(char) * fmt, int ns, int * n)
{
    int err = 0;

    if (*fmt == '%') {
        const(char) * p;
        long len = strtol(++fmt, &p, 10);

        if (p == fmt) {
            /* no "max rows" given */
            *n = ns;
        } else if (len < 0) {
            err = GretlError.invarg;
        } else if (len > ns) {
            *n = ns;
        } else {
            *n = cast(int) len;
        }
        if (!err && strcmp(p, "m")) {
            err = GretlError.invarg;
        }
    } else {
        err = GretlError.invarg;
    }

    return err;
}

GretlMatrixRaw * vector_from_strings(char ** S, int ns,
                                  const(char) * fmt,
                                  int * nvals,
                                  int * err)
{
    GretlMatrixRaw * m = null;
    int n = 0;

    *nvals = 0;
    *err = check_matrix_scan_format(fmt, ns, &n);

    if (!*err) {
        if (n == 0) {
            m = gretl_null_matrix_new();
        } else {
            m = gretl_column_vector_alloc(n);
        }
        if (m == null) {
            *err = GretlError.alloc;
        }
    }

    if (!*err && n > 0) {
        char * s;
        char * p;
        double x;
        int i;

        gretl_push_c_numeric_locale();
        for (i=0; i<n; i++) {
            s = S[i];
            if (s == null || *s == '\0') {
                m.val[i] = NADBL;
            } else {
                x = strtod(s, &p);
                if (*p == '\0') {
                    m.val[i] = x;
                    *nvals += 1;
                } else {
                    m.val[i] = NADBL;
                }
            }
        }
        gretl_pop_c_numeric_locale();
    }

    return m;
}
