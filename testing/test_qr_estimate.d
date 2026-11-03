/* Tests for qr_estimate.d: the QR (plain and pivoting) and Cholesky
   regression routines, robust covariance estimators, and the HAC helpers.
   Besides checking that the calls run, we check that the different
   estimators agree on the coefficients (which they must), and that the
   robust standard errors are finite. */

import core.stdc.math : fabs, sqrt;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;
import core.stdc.string : strcpy;

import gretl.dataset : create_new_dataset;
import gretl.estimate : lsq;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_list : gretl_list_new;
import gretl.gretl_matrix : gretl_matrix_alloc, gretl_matrix_free;
import gretl.gretl_matrixh;
import gretl.libgretlh;
import gretl.libset : libset_get_int, libset_set_bool, libset_set_int;
import gretl.libseth : SetKey;
import gretl.missingh : na;
import gretl.options : set_optval_string;
import gretl.qr_estimate;

int n_tests, n_failed;

void check(bool ok, string what)
{
    n_tests++;
    if (!ok) n_failed++;
    printf("%s %.*s\n", ok ? "ok:  ".ptr : "FAIL:".ptr, cast(int) what.length, what.ptr);
}

ulong lcg_state = 4242;

double urand()
{
    lcg_state = lcg_state * 6364136223846793005UL + 1442695040888963407UL;
    return ((lcg_state >> 11) / cast(double) (1UL << 53)) - 0.5;
}

/* 0 const, 1 y, 2 x1, 3 x2, 4 cluster id, 5 x3 = x1 + x2 (collinear) */

DatasetRaw * make_dataset(int T)
{
    DatasetRaw * dset = create_new_dataset(6, T, 0);

    foreach (i, nm; ["const", "y", "x1", "x2", "grp", "x3"]) {
        strcpy(dset.varname[i], (nm ~ "\0").ptr);
    }
    dset.pd = 1;
    dset.structure = TsCodes.time_series;
    dset.t1 = 0;
    dset.t2 = T - 1;
    foreach (t; 0 .. T) {
        double x1 = urand(), x2 = urand();

        dset.Z[2][t] = x1;
        dset.Z[3][t] = x2;
        dset.Z[1][t] = 1.0 + 2.0 * x1 - 1.5 * x2 + 0.5 * urand() * (1 + fabs(x1));
        dset.Z[4][t] = t % 5;
        dset.Z[5][t] = x1 + x2;
    }
    return dset;
}

int * mklist(int[] v)
{
    int * l = gretl_list_new(cast(int) v.length);

    foreach (i, x; v) l[i+1] = x;
    return l;
}

bool close_to(double a, double b, double tol = 1e-8)
{
    return !na(a) && !na(b) && fabs(a - b) <= tol * (1 + fabs(b));
}

void main()
{
    setvbuf(stdout, null, _IOLBF, 4096);
    DatasetRaw * dset = make_dataset(60);
    int err;

    /* baseline: Cholesky (the default) */
    Model mc = lsq(mklist([1, 0, 2, 3]), dset, GretlCmdIndex.ols, Gretlopt.OPT_NONE);

    check(mc.errcode == 0 && mc.ncoeff == 3, "lsq (Cholesky) estimates the model");
    printf("Cholesky coeffs: %g %g %g\n", mc.coeff[0], mc.coeff[1], mc.coeff[2]);

    /* plain QR (forced with the use_qr setting) */
    libset_set_int(SetKey.use_qr, 1);
    Model m1 = lsq(mklist([1, 0, 2, 3]), dset, GretlCmdIndex.ols, Gretlopt.OPT_NONE);

    check(m1.errcode == 0, "lsq with use_qr=1 (plain QR)");
    foreach (i; 0 .. 3) {
        check(close_to(m1.coeff[i], mc.coeff[i]) && close_to(m1.sderr[i], mc.sderr[i]),
              "QR and Cholesky agree on coefficient and standard error");
    }
    check(close_to(m1.ess, mc.ess) && close_to(m1.rsq, mc.rsq), "QR and Cholesky agree on ESS and R-squared");

    /* QR with column pivoting */
    libset_set_int(SetKey.use_qr, 2);
    Model m2 = lsq(mklist([1, 0, 2, 3]), dset, GretlCmdIndex.ols, Gretlopt.OPT_NONE);

    check(m2.errcode == 0 && close_to(m2.coeff[1], mc.coeff[1]), "pivoting QR agrees with Cholesky");

    /* collinear regressor: plain and pivoting QR should both drop one term */
    libset_set_int(SetKey.use_qr, 1);
    Model m3 = lsq(mklist([1, 0, 2, 3, 5]), dset, GretlCmdIndex.ols, Gretlopt.OPT_NONE);

    printf("collinear (plain QR): errcode %d, ncoeff %d\n", m3.errcode, m3.ncoeff);
    check(m3.errcode == 0 && m3.ncoeff == 3, "collinear regressor dropped (plain QR)");
    libset_set_int(SetKey.use_qr, 2);
    Model m4 = lsq(mklist([1, 0, 2, 3, 5]), dset, GretlCmdIndex.ols, Gretlopt.OPT_NONE);

    printf("collinear (pivoting QR): errcode %d, ncoeff %d\n", m4.errcode, m4.ncoeff);
    check(m4.errcode == 0 && m4.ncoeff == 3, "collinear regressor dropped (pivoting QR)");
    libset_set_int(SetKey.use_qr, 0);

    /* robust covariance estimators (all of these go through the QR code).
       For time-series data lsq() uses HAC unless force_hc is set. */
    libset_set_bool(SetKey.force_hc, 1);
    double[4] se1;

    foreach (v; 0 .. 4) {
        libset_set_int(SetKey.hc_version, v);
        Model mr = lsq(mklist([1, 0, 2, 3]), dset, GretlCmdIndex.ols, Gretlopt.OPT_R);

        printf("HC%d: errcode %d, se = %g %g %g\n", v, mr.errcode, mr.sderr[0], mr.sderr[1], mr.sderr[2]);
        check(mr.errcode == 0 && mr.sderr[1] > 0 && close_to(mr.coeff[1], mc.coeff[1]),
              "robust (HC) standard errors computed, coefficients unchanged");
        se1[v] = mr.sderr[1];
    }
    check(close_to(se1[1], se1[0] * sqrt(60.0 / (60.0 - 3.0))), "HC1 = HC0 scaled by T/(T-k)");
    check(se1[2] > se1[0] && se1[3] > se1[2], "HC2 and HC3 inflate the standard error");
    libset_set_int(SetKey.hc_version, 1);
    libset_set_bool(SetKey.force_hc, 0);

    foreach (kern; 0 .. 3) {
        foreach (pw; 0 .. 2) {
            libset_set_int(SetKey.hac_kernel, kern);
            libset_set_bool(SetKey.prewhiten, pw);
            Model mh = lsq(mklist([1, 0, 2, 3]), dset, GretlCmdIndex.ols, Gretlopt.OPT_R);

            printf("HAC kernel %d prewhiten %d: errcode %d, se = %g %g %g\n", kern, pw,
                   mh.errcode, mh.sderr[0], mh.sderr[1], mh.sderr[2]);
            check(mh.errcode == 0 && mh.sderr[1] > 0, "HAC standard errors computed");
        }
    }
    libset_set_int(SetKey.hac_kernel, 0);
    libset_set_bool(SetKey.prewhiten, 0);

    /* clustered standard errors (clustering variable chosen through the option store) */
    set_optval_string(GretlCmdIndex.ols, Gretlopt.OPT_C, "grp");
    Model mcl = lsq(mklist([1, 0, 2, 3]), dset, GretlCmdIndex.ols, Gretlopt.OPT_R | Gretlopt.OPT_C);

    printf("clustered: errcode %d, se = %g %g %g\n", mcl.errcode, mcl.sderr[0], mcl.sderr[1], mcl.sderr[2]);
    check(mcl.errcode == 0 && mcl.sderr[1] > 0 && mcl.dfd == 4, "clustered standard errors (5 clusters)");

    /* Durbin-Watson via the QR path (OPT_I with OPT_R) is exercised by the ADF/EG tests */

    /* HAC helpers */
    check(close_to(hac_weight(0, 4, 2), 1.0 - 2.0 / 5.0), "hac_weight (Bartlett)");
    check(hac_weight(1, 4, 1) > hac_weight(1, 4, 3), "hac_weight (Parzen) decreases with lag");
    check(!na(qs_hac_weight(3.0, 2)), "qs_hac_weight returns a number");

    GretlMatrixRaw * H = gretl_matrix_alloc(40, 2);

    foreach (t; 0 .. 40) {
        gretl_matrix_set(H, t, 0, 1.0);
        gretl_matrix_set(H, t, 1, urand());
    }
    int m; double bt;
    err = newey_west_bandwidth(H, null, 0, 0, &m, &bt);
    printf("newey_west_bandwidth: err %d, m %d, bt %g\n", err, m, bt);
    check(err == 0 && m >= 0, "newey_west_bandwidth");

    GretlMatrixRaw * XOX = HAC_XOX(H, null, null, 0, &err);

    check(err == 0 && XOX !is null && XOX.rows == 2, "HAC_XOX");
    gretl_matrix_free(XOX);
    GretlMatrixRaw * lr = long_run_covariance(H, 1, &err);

    check(err == 0 && lr !is null, "long_run_covariance");
    gretl_matrix_free(lr);
    GretlMatrixRaw * og = newey_west_OPG(H, &err);

    check(err == 0 && og !is null, "newey_west_OPG");
    gretl_matrix_free(og);
    gretl_matrix_free(H);

    printf("\n%d checks, %d failed\n", n_tests, n_failed);
    if (n_failed) {
        throw new Exception("some tests failed");
    }
}
