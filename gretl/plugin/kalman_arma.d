/* Port of plugin/kalman_arma.c: exact ML for ARMA using the gretl-native Kalman filter (module gretl.kalman).

   As of January 2022 gretl needs this for one case that is not handled by the faster dedicated ARMA codes AS 154 and AS 197,
   namely ARIMA (with non-zero integration order) plus within-sample missing values.

   In the C this file is #included by arma.c and so can call the static functions of arma.c; here they are imported from gretl.plugin.arma
   (form_arma_y_vector, form_arma_X_matrix, arma_QML_vcv, arma_OPG_vcv, arma_use_opg, arima_ydiff_only, extern (C) there). The callbacks
   kalman_arma_ll and kalman_arma_llt_callback, and kalman_arma itself, are public because arma.c uses them. */

module gretl.plugin.kalman_arma;

import core.stdc.math : sqrt;
import core.stdc.stdio : fprintf, stderr;
import core.stdc.stdlib : free, malloc;

import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_bfgs : BFGS_defaults, BFGS_max, numerical_hessian_inverse;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_matrix;
import gretl.gretl_matrixh;
import gretl.gretl_model : arima_delta_coeffs, gretl_model_allocate_storage, gretl_model_set_int, gretl_model_set_matrix_as_data, gretl_model_set_vcv_info, gretl_model_write_vcv;
import gretl.gretl_modelh : ArmaFlags, MLVCVType, VCVMajorType;
import gretl.gretl_prn : Prn, pputs;
import gretl.gretl_utils : copyvec, doubles_array_free, doubles_array_new;
import gretl.kalman;
import gretl.libgretlh;
import gretl.libset : libset_get_bool, libset_set_bool;
import gretl.libseth : SetKey;
import gretl.missingh : NADBL, na;
import gretl.nlsh : CriterionType;
import gretl.plugin.arma_priv;
import gretl.plugin.arma : arima_ydiff_only, arma_OPG_vcv, arma_model_add_roots, maybe_correct_MA, write_arma_model_stats, arma_QML_vcv, arma_use_opg, form_arma_X_matrix, form_arma_y_vector;

enum ARMA_DEBUG = 0;
enum KALMAN_ALL = 999;


struct KHelper {
    GretlMatrixBlock * Bk;
    GretlMatrixRaw * a;
    GretlMatrixRaw * P;
    GretlMatrixRaw * T;
    GretlMatrixRaw * B;
    GretlMatrixRaw * Z;
    GretlMatrixRaw * Q;
    GretlMatrixRaw * V;
    GretlMatrixRaw * avar;

    GretlMatrixRaw * avar2;
    GretlMatrixRaw * vQ;

    GretlMatrixRaw * T_; /* used only for ARIMA via levels */
    GretlMatrixRaw * Q_; /* ditto */
    GretlMatrixRaw * P_; /* ditto */

    ArmaInfo *kainfo;
};

private __gshared int kalman_do_ma_check = 1;

private void kalman_helper_free(KHelper * kh)
{
    if (kh != null) {
        gretl_matrix_block_destroy(kh.Bk);
        gretl_matrix_free(kh.avar2);
        gretl_matrix_free(kh.vQ);
        gretl_matrix_free(kh.T_);
        gretl_matrix_free(kh.Q_);
        gretl_matrix_free(kh.P_);
        free(kh);
    }
}

private KHelper * kalman_helper_new(ArmaInfo * ainfo, int r, int k)
{
    KHelper * kh;
    int r0;
    int r2;
    int err = 0;

    kh = cast(KHelper *) malloc((*kh).sizeof);
    if (kh == null) {
        return null;
    }

    r0 = ainfo.r0;
    r2 = r0 * r0;

    kh.avar2 = kh.vQ = null;
    kh.T_ = kh.Q_ = kh.P_ = null;

    kh.Bk = gretl_matrix_block_new(&kh.a, r, 1,
				    &kh.P, r, r,
				    &kh.T, r, r,
				    &kh.B, k, 1,
				    &kh.Z, r, 1,
				    &kh.Q, r, r,
				    &kh.V, ainfo.fullT, 1,
				    &kh.avar, r2, r2,
				    null);

    if (kh.Bk == null) {
        err = GretlError.alloc;
    } else if (arma_using_vech(ainfo)) {
        int m = r0 * (r0 + 1) / 2;

        kh.avar2 = gretl_matrix_alloc(m, m);
        kh.vQ = gretl_column_vector_alloc(m);
        if (kh.avar2 == null || kh.vQ == null) {
            err = GretlError.alloc;
        }
    } else {
        kh.vQ = gretl_column_vector_alloc(r2);
        if (kh.vQ == null) {
            err = GretlError.alloc;
        }
    }

    if (!err && arima_levels(ainfo)) {
        kh.T_ = gretl_matrix_alloc(r0, r0);
        kh.Q_ = gretl_matrix_alloc(r0, r0);
        kh.P_ = gretl_matrix_alloc(r0, r0);
        if (kh.T_ == null || kh.Q_ == null || kh.P_ == null) {
            err = GretlError.alloc;
        }
    }

    if (err) {
        kalman_helper_free(kh);
        kh = null;
    } else {
        kh.kainfo = ainfo;
    }

    return kh;
}

/* Get the dimension of the state-space representation: note
   that this is augmented if we're estimating an ARIMA
   model using the levels formulation in order to handle
   missing values -- see Harvey and Pierse, "Estimating
   Missing Observations in Economic Time Series", JASA 1984.
*/

private int ainfo_get_state_size(ArmaInfo * ainfo)
{
    int plen = ainfo.p + ainfo.pd * ainfo.P;
    int qlen = ainfo.q + ainfo.pd * ainfo.Q;
    int r = (plen > qlen + 1)? plen : qlen + 1;

    ainfo.r0 = r;

    if (arima_levels(ainfo)) {
        r += ainfo.d + ainfo.pd * ainfo.D;
    }

    return r;
}

private int allocate_ac_mc(ArmaInfo * ainfo)
{
    int m = (ainfo.P > 0) + (ainfo.Q > 0);
    int err = 0;

    if (m > 0) {
        double * ac = null;
        double * mc = null;
        int n;
        int i = 0;

        ainfo.aux = doubles_array_new(m, 0);
        if (ainfo.aux == null) {
            return GretlError.alloc;
        }

        if (ainfo.P > 0) {
            n = 1 + ainfo.p + ainfo.pd * ainfo.P;
            ac = cast(double *) malloc(n * (*ac).sizeof);
            if (ac == null) {
                err = GretlError.alloc;
            } else {
                ainfo.aux[i++] = ac;
            }
        }

        if (!err && ainfo.Q > 0) {
            n = 1 + ainfo.q + ainfo.pd * ainfo.Q;
            mc = cast(double *) malloc(n * (*mc).sizeof);
            if (mc == null) {
                err = GretlError.alloc;
            } else {
                ainfo.aux[i++] = mc;
            }
        }

        if (err) {
            doubles_array_free(ainfo.aux, m);
        } else {
            ainfo.n_aux = m;
        }
    }

    return err;
}

private void write_big_phi(const(double) * phi, const(double) * Phi, ArmaInfo * ainfo, GretlMatrixRaw * F)
{
    int pmax = ainfo.p + ainfo.pd * ainfo.P;
    double * ac = ainfo.aux[0];
    double x;
    double y;
    int i;
    int j;
    int k;
    int ii;

    for (i=0; i<=pmax; i++) {
        ac[i] = 0.0;
    }

    for (j=-1; j<ainfo.P; j++) {
        x = (j < 0)? -1 : Phi[j];
        k = 0;
        for (i=-1; i<ainfo.p; i++) {
            if (i < 0) {
                y = -1;
            } else if (AR_included(ainfo, i)) {
                y = phi[k++];
            } else {
                y = 0.0;
            }
            ii = (j+1) * ainfo.pd + (i+1);
            ac[ii] -= x * y;
        }
    }

    for (i=0; i<pmax; i++) {
        gretl_matrix_set(F, 0, i, ac[i+1]);
    }
}

private void write_big_theta(const(double) * theta, const(double) * Theta, ArmaInfo * ainfo, GretlMatrixRaw * H, GretlMatrixRaw * F)
{
    int qmax = ainfo.q + ainfo.pd * ainfo.Q;
    int i = (ainfo.P > 0)? 1 : 0;
    double * mc = ainfo.aux[i];
    double x;
    double y;
    int j;
    int k;
    int ii;

    for (i=0; i<=qmax; i++) {
        mc[i] = 0.0;
    }

    for (j=-1; j<ainfo.Q; j++) {
        x = (j < 0)? 1 : Theta[j];
        k = 0;
        for (i=-1; i<ainfo.q; i++) {
            if (i < 0) {
                y = 1;
            } else if (MA_included(ainfo, i)) {
                y = theta[k++];
            } else {
                y = 0.0;
            }
            ii = (j+1) * ainfo.pd + (i+1);
            mc[ii] += x * y;
        }
    }

    for (i=1; i<=qmax; i++) {
        if (H != null) {
            H.val[i] = mc[i];
        } else {
            gretl_matrix_set(F, ainfo.r0, i, mc[i]);
        }
    }
}

private void condense_row(GretlMatrixRaw * targ, const(GretlMatrixRaw) * src, int targrow, int srcrow, int n)
{
    double x;
    int i;
    int j;
    int k;
    int g;
    int targcol = 0;

    for (j=0; j<n; j++) {
        for (i=j; i<n; i++) {
            k = j * n + i;
            g = (k % n) * n + k / n;
            x = gretl_matrix_get(src, srcrow, k);
            if (g != k) {
                x += gretl_matrix_get(src, srcrow, g);
            }
            gretl_matrix_set(targ, targrow, targcol++, x);
        }
    }
}

private void condense_state_vcv(GretlMatrixRaw * targ, const(GretlMatrixRaw) * src, int n)
{
    int posr = 0;
    int posc = 0;
    int i;
    int j;

    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            if (j >= i) {
                condense_row(targ, src, posr++, posc, n);
            }
            posc++;
        }
    }
}

private int kalman_matrices_init(ArmaInfo * ainfo, KHelper * kh, const(double) * y)
{
    int r0 = ainfo.r0;
    int r = kh.T.rows;

    gretl_matrix_zero(kh.B);
    gretl_matrix_zero(kh.a);
    gretl_matrix_zero(kh.P);
    gretl_matrix_zero(kh.T);
    gretl_matrix_inscribe_I(kh.T, 1, 0, r0 - 1);

    gretl_matrix_zero(kh.Q);
    gretl_matrix_set(kh.Q, 0, 0, 1.0);

    gretl_matrix_zero(kh.Z);
    gretl_vector_set(kh.Z, 0, 1.0);

    if (arima_levels(ainfo)) {
        /* write additional constant elements of T, Z and a */
        int d = ainfo.d;
        int D = ainfo.D;
        int s = ainfo.pd;
        int i;
        int k = d + s * D;
        int *c = arima_delta_coeffs(d, D, s);
        double y0;

        if (c == null) {
            return GretlError.alloc;
        }
        for (i=0; i<k; i++) {
            gretl_matrix_set(kh.T, r0, r0 + i, c[i]);
        }
        gretl_matrix_set(kh.T, r0, 0, 1.0);
        if (r - r0 > 1) {
            gretl_matrix_inscribe_I(kh.T, r0 + 1, r0, k - 1);
        }
        for (i=0; i<k; i++) {
            gretl_vector_set(kh.Z, r0 + i, c[i]);
            /* lagged data */
            y0 = y[ainfo.t1 - 1 - i];
            if (ainfo.yscale != 1.0 && !na(y0)) {
                y0 -= ainfo.yshift;
                y0 *= ainfo.yscale;
            }
            gretl_vector_set(kh.a, r0 + i, y0);
        }
        free(c);

static if (ARMA_DEBUG) {
        gretl_matrix_print(kh.a, "a0 (arima via levels)");
}
        /* initialize the plain-arma "shadow" matrices */
        gretl_matrix_zero(kh.T_);
        gretl_matrix_inscribe_I(kh.T_, 1, 0, r0 - 1);
        gretl_matrix_zero(kh.Q_);
        gretl_matrix_set(kh.Q_, 0, 0, 1.0);
        gretl_matrix_zero(kh.P_);
    } else if (ainfo.np == 0 && ainfo.P == 0) {
        /* initialize P to identity matrix */
        gretl_matrix_inscribe_I(kh.P, 0, 0, kh.P.rows);
    }

    return 0;
}

private int write_kalman_matrices(KHelper * kh, const(double) * b, int idx)
{
    ArmaInfo * ainfo = kh.kainfo;
    const(double) * phi = b + ainfo.ifc;
    const(double) * Phi = phi + ainfo.np;
    const(double) * theta = Phi + ainfo.P;
    const(double) * Theta = theta + ainfo.nq;
    const(double) * beta = Theta + ainfo.Q;
    double mu = (ainfo.ifc)? b[0] : 0.0;
    int rewrite_B = 0;
    int rewrite_T = 0;
    int rewrite_Z = 0;
    int i;
    int k;
    int err = 0;

    if (idx == KALMAN_ALL) {
        rewrite_B = rewrite_T = rewrite_Z = 1;
    } else {
        /* called in context of calculating score, for OPG matrix */
        int pmax = ainfo.ifc + ainfo.np + ainfo.P;
        int tmax = pmax + ainfo.nq + ainfo.Q;

        if (ainfo.ifc && idx == 0) {
            rewrite_B = 1;
        } else if (idx >= ainfo.ifc && idx < pmax) {
            rewrite_T = 1;
        } else if (idx >= ainfo.ifc && idx < tmax) {
            rewrite_Z = 1;
        } else {
            rewrite_B = 1;
        }
    }

    /* revise for pure MA model */
    if (ainfo.np == 0 && ainfo.P == 0 && !arima_levels(ainfo)) {
        rewrite_T = 0;
    }
    /* and for case of no constant or other regressors */
    if (ainfo.ifc == 0 && ainfo.nexo == 0) {
        rewrite_B = 0;
    }

    /* See Hamilton, Time Series Analysis, ch 13, p. 375 */

    if (rewrite_B) {
        /* const and coeffs on exogenous vars */
        gretl_vector_set(kh.B, 0, mu);
        for (i=0; i<ainfo.nexo; i++) {
            gretl_vector_set(kh.B, i + 1, beta[i]);
        }
    }

    if (rewrite_Z) {
        /* form the Z' vector using theta and/or Theta */
        if (ainfo.Q > 0) {
            write_big_theta(theta, Theta, ainfo, kh.Z, null);
        } else {
            k = 0;
            for (i=0; i<ainfo.q; i++) {
                if (MA_included(ainfo, i)) {
                    gretl_vector_set(kh.Z, i+1, theta[k++]);
                } else {
                    gretl_vector_set(kh.Z, i+1, 0.0);
                }
            }
        }
    }

    if (rewrite_T) {
        /* form the T matrix using phi and/or Phi */
        GretlMatrixRaw * T = (kh.T_ != null)? kh.T_ : kh.T;
        GretlMatrixRaw * Q = (kh.Q_ != null)? kh.Q_ : kh.Q;
        GretlMatrixRaw * P = (kh.P_ != null)? kh.P_ : kh.P;

        if (ainfo.P > 0) {
            write_big_phi(phi, Phi, ainfo, T);
        } else {
            k = 0;
            for (i=0; i<ainfo.p; i++) {
                if (AR_included(ainfo, i)) {
                    gretl_matrix_set(T, 0, i, phi[k++]);
                } else {
                    gretl_matrix_set(T, 0, i, 0.0);
                }
            }
        }

        if (arima_levels(ainfo)) {
            /* the full T matrix incorporates \theta */
            if (ainfo.Q > 0) {
                write_big_theta(theta, Theta, ainfo, null, kh.T);
            } else {
                k = 0;
                for (i=0; i<ainfo.q; i++) {
                    if (MA_included(ainfo, i)) {
                        gretl_matrix_set(kh.T, ainfo.r0, i+1, theta[k++]);
                    } else {
                        gretl_matrix_set(kh.T, ainfo.r0, i+1, 0.0);
                    }
                }
            }
        }

        /* form $P_{1|0}$ (MSE) matrix, as per Hamilton, ch 13, p. 378. */

        gretl_matrix_kronecker_product(T, T, kh.avar);
        gretl_matrix_I_minus(kh.avar);
        if (arma_using_vech(ainfo)) {
            condense_state_vcv(kh.avar2, kh.avar, gretl_matrix_rows(T));
            gretl_matrix_vectorize_h(kh.vQ, Q);
            err = gretl_LU_solve(kh.avar2, kh.vQ);
            if (!err) {
                gretl_matrix_unvectorize_h(P, kh.vQ);
            }
        } else {
            gretl_matrix_vectorize(kh.vQ, Q);
            err = gretl_LU_solve(kh.avar, kh.vQ);
            if (!err) {
                gretl_matrix_unvectorize(P, kh.vQ);
            }
        }
    }

    if (arima_levels(ainfo)) {
        /* complete the job on T, Q, P */
        gretl_matrix_inscribe_matrix(kh.T, kh.T_, 0, 0, GretlMatrixMod.none);
        gretl_matrix_inscribe_matrix(kh.Q, kh.Q_, 0, 0, GretlMatrixMod.none);
        gretl_matrix_inscribe_matrix(kh.P, kh.P_, 0, 0, GretlMatrixMod.none);
    }

    return err;
}

private int rewrite_kalman_matrices(Kalman * K, const(double) * b, int i)
{
    KHelper * kh = cast(KHelper *) kalman_get_data(K);
    int err = write_kalman_matrices(kh, b, i);

    if (!err) {
        kalman_set_initial_state_vector(K, kh.a);
        kalman_set_initial_MSE_matrix(K, kh.P);
    }

    return err;
}

/* used only in obtaining the OPG, if wanted */

extern (C) const(double) * kalman_arma_llt_callback(const(double) * b, int i, void * data)
{
    Kalman * K = cast(Kalman *) data;
    KHelper * kh = cast(KHelper *) kalman_get_data(K);
    int err;

    rewrite_kalman_matrices(K, b, i);
    err = kfilter_standard(K, null);

    return (err)? null : kh.V.val;
}

extern (C) double kalman_arma_ll(const(double) * b, void * data)
{
    Kalman * K = cast(Kalman *) data;
    KHelper * kh = cast(KHelper *) kalman_get_data(K);
    ArmaInfo * ainfo = kh.kainfo;
    int offset = ainfo.ifc + ainfo.np + ainfo.P;
    double * theta = cast(double *) b + offset;
    double * Theta = theta + ainfo.nq;
    double ll = NADBL;
    int err = 0;

    if (kalman_do_ma_check && maybe_correct_MA(ainfo, theta, Theta)) {
        pputs(kalman_get_printer(K), _("MA estimate(s) out of bounds\n").ptr);
        return NADBL;
    }

    err = rewrite_kalman_matrices(K, b, KALMAN_ALL);

    if (!err) {
        err = kfilter_standard(K, null);
        ll = kalman_get_loglik(K);
    }

    return ll;
}

private int kalman_arma_finish(Model * pmod, ArmaInfo * ainfo, const(DatasetRaw) * dset, Kalman * K, double * b, Gretlopt opt, Prn * prn)
{
    KHelper * kh = cast(KHelper *) kalman_get_data(K);
    int do_opg = arma_use_opg(opt);
    int i;
    int t;
    int k = ainfo.nc;
    int QML = (opt & Gretlopt.OPT_R);
    double s2;
    int err;

    pmod.t1 = ainfo.t1;
    pmod.t2 = ainfo.t2;
    pmod.nobs = ainfo.T;
    pmod.ncoeff = ainfo.nc;
    pmod.full_n = dset.n;

    /* in the Kalman case the basic model struct is empty, so we
       have to allocate for coefficients, residuals and so on
    */
    err = gretl_model_allocate_storage(pmod);
    if (err) {
        return err;
    }

    for (i=0; i<k; i++) {
        pmod.coeff[i] = b[i];
    }

    i = 0;
    for (t=pmod.t1; t<=pmod.t2; t++) {
        pmod.uhat[t] = gretl_vector_get(kh.V, i++);
    }

    s2 = kalman_get_arma_variance(K);
    pmod.sigma = sqrt(s2);
    pmod.lnL = kalman_get_loglik(K);

    if (!do_opg) {
        /* covariance matrix based on Hessian (perhaps QML) */
        GretlMatrixRaw * Hinv;
        double d = 0.0; /* adjust? */

        kalman_do_ma_check = 0;
        Hinv = numerical_hessian_inverse(b, ainfo.nc, &kalman_arma_ll,
                                         K, d, &err);
        kalman_do_ma_check = 1;
        if (!err) {
            if (QML) {
                err = arma_QML_vcv(pmod, Hinv, K, 0, b, s2, k, ainfo.T, prn);
            } else {
                err = gretl_model_write_vcv(pmod, Hinv);
                if (!err) {
                    gretl_model_set_vcv_info(pmod, VCVMajorType.ml, MLVCVType.hessian);
                }
            }
        } else if (!(opt & Gretlopt.OPT_H)) {
            /* fallback when Hessian not explicitly requested */
            err = 0;
            do_opg = 1;
            gretl_model_set_int(pmod, "hess-error", 1);
        }
        gretl_matrix_free(Hinv);
    }

    if (do_opg) {
        err = arma_OPG_vcv(pmod, K, 0, b, s2, k, ainfo.T, prn);
        if (!err) {
            gretl_model_set_vcv_info(pmod, VCVMajorType.ml, MLVCVType.op);
            pmod.opt |= Gretlopt.OPT_G;
        }
    }

    if (!err) {
        write_arma_model_stats(pmod, ainfo, dset);
        arma_model_add_roots(pmod, ainfo, b);
        gretl_model_set_int(pmod, "arma_flags", ArmaFlags.exact);
        if (arma_lbfgs(ainfo)) {
            pmod.opt |= Gretlopt.OPT_L;
        }
        if (arima_ydiff_only(ainfo)) {
            pmod.opt |= Gretlopt.OPT_Y;
        }
    }

    return err;
}

private int kalman_undo_y_scaling(ArmaInfo * ainfo, GretlMatrixRaw * y, double * b, Kalman * K)
{
    double *beta = b + ainfo.ifc + ainfo.np + ainfo.P +
        ainfo.nq + ainfo.Q;
    int i;
    int t;
    int T = ainfo.t2 - ainfo.t1 + 1;
    int err = 0;

    if (ainfo.ifc) {
        b[0] /= ainfo.yscale;
        b[0] += ainfo.yshift;
    }

    for (i=0; i<ainfo.nexo; i++) {
        beta[i] /= ainfo.yscale;
    }

    i = ainfo.t1;
    for (t=0; t<T; t++) {
        y.val[t] /= ainfo.yscale;
        y.val[t] += ainfo.yshift;
    }

    if (na(kalman_arma_ll(b, K))) {
        err = 1;
    }

    return err;
}

private void free_arma_X_matrix(ArmaInfo * ainfo, GretlMatrixRaw * X)
{
    if (X == ainfo.dX) {
        gretl_matrix_free(ainfo.dX);
        ainfo.dX = null;
    } else {
        gretl_matrix_free(X);
    }
}

private int add_smoothed_y(Kalman * K, Model * pmod, ArmaInfo * ainfo)
{
    GretlMatrixRaw * S;
    int err = 0;

    S = kalman_smooth(K, Gretlopt.OPT_NONE, null, &err);

static if (0) {
    fprintf(stderr, "HERE ainfo.r0 = %d\n", ainfo.r0);
    fprintf(stderr, " ainfo.t1 = %d, ainfo.t2 = %d\n", ainfo.t1, ainfo.t2);
}

    if (S != null) {
	GretlMatrixRaw * m = gretl_matrix_alloc(S.rows, 1);
	int t;

	if (m == null) {
	    err = GretlError.alloc;
	} else {
	    for (t=0; t<S.rows; t++) {
		m.val[t] = gretl_matrix_get(S, t+1, ainfo.r0);
	    }
	    m.val[S.rows-1] = NADBL;
	    gretl_matrix_set_t1(m, ainfo.t1);
	    gretl_matrix_set_t2(m, ainfo.t2);
	    gretl_model_set_matrix_as_data(pmod, "smstate", m);
	}
	gretl_matrix_free(S);
    }

    return err;
}

int kalman_arma(const(double) * coeff, const(DatasetRaw) * dset, ArmaInfo * ainfo, Model * pmod, Gretlopt opt)
{
    Kalman * K = null;
    KHelper * kh = null;
    GretlMatrixRaw * y = null;
    GretlMatrixRaw * X = null;
    int r;
    int k = 1 + ainfo.nexo; /* number of exog vars plus space for const */
    double * b;
    int err = 0;

    b = copyvec(coeff, ainfo.nc);
    if (b == null) {
        return GretlError.alloc;
    }

    y = form_arma_y_vector(ainfo, &err);

    if (!err && ainfo.nexo > 0) {
        if (ainfo.dX != null) {
            X = ainfo.dX;
        } else {
            X = form_arma_X_matrix(ainfo, dset, &err);
        }
    }

    if (!err) {
        err = allocate_ac_mc(ainfo);
    }

    if (err) {
        goto bailout;
    }

    r = ainfo_get_state_size(ainfo);

    /* when should we use vech apparatus? */
    if (r > 4) {
        set_arma_use_vech(ainfo);
    }

    kh = kalman_helper_new(ainfo, r, k);
    if (kh == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    kalman_matrices_init(ainfo, kh, dset.Z[ainfo.yno]);

    K = kalman_new(kh.a, kh.P, kh.T, kh.B, kh.Z, kh.Q,
                   null, y, X, null, kh.V, &err);

    if (err) {
        fprintf(stderr, "kalman_new(): err = %d\n", err);
    } else {
        int save_lbfgs = libset_get_bool(SetKey.use_lbfgs);
        double toler;
        int maxit;

        kalman_attach_printer(K, ainfo.prn);
        kalman_attach_data(K, kh);
	kalman_set_arma_ll(K);

        BFGS_defaults(&maxit, &toler, GretlCmdIndex.arma);

        if (save_lbfgs) {
            ainfo.pflags |= PrivFlags.lbfgs;
        } else if (opt & Gretlopt.OPT_L) {
            libset_set_bool(SetKey.use_lbfgs, 1);
            ainfo.pflags |= PrivFlags.lbfgs;
        }

        err = BFGS_max(b, ainfo.nc, maxit, toler,
                       &ainfo.fncount, &ainfo.grcount,
                       &kalman_arma_ll, CriterionType.loglik,
                       null, K, null, opt | Gretlopt.OPT_A,
                       ainfo.prn);

        if (save_lbfgs == 0 && (opt & Gretlopt.OPT_L)) {
            libset_set_bool(SetKey.use_lbfgs, 0);
        }

        if (err) {
            fprintf(stderr, "kalman_arma: optimizer returned %d\n", err);
        } else if (ainfo.yscale == 1.0 && ainfo.r0 > 0) {
	    /* experimental, could handle non-unit yscale with some work? */
	    add_smoothed_y(K, pmod, ainfo);
	}
    }

    if (!err && ainfo.yscale != 1.0) {
        kalman_undo_y_scaling(ainfo, y, b, K);
    }

    if (!err) {
	gretl_model_set_int(pmod, "fncount", ainfo.fncount);
	gretl_model_set_int(pmod, "grcount", ainfo.grcount);
        err = kalman_arma_finish(pmod, ainfo, dset, K, b,
                                 opt, ainfo.prn);
    }

 bailout:

    if (err) {
        pmod.errcode = err;
    }

    kalman_free(K);
    kalman_helper_free(kh);

    gretl_matrix_free(y);
    free_arma_X_matrix(ainfo, X);
    free(b);

    return err;
}
