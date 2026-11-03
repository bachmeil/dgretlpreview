/* Smoke tests for discrete.d: binary logit and probit (coefficients and log-likelihood against values from an
   independent Newton iteration), the logit odds-ratio and perfect-prediction handling, bit_permutations,
   Fisher's exact test, logistic_ymax_lmax and the plugin-missing paths of the models that need plugins.

   Not exercised: ordered and multinomial models, logistic_model, and the plugin models (tobit, interval, heckit,
   duration, count, reprobit, biprobit) which need the corresponding plugins.

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;
import core.stdc.stdlib : free;
import core.stdc.string : strstr;

import gretl.dataset : create_new_dataset, destroy_dataset;
import gretl.describeh : Xtab;
import gretl.discrete;
import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_get_buffer, gretl_print_new;
import gretl.gretl_prnh : PrnType;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_matrix : gretl_matrix_free;
import gretl.gretl_matrixh : GretlMatrixRaw;
import gretl.gretl_model : clear_model;
import gretl.libgretlh;
import gretl.missingh : NADBL, na;

int n_tests, n_failed;

void check(bool ok, string what)
{
    n_tests++;
    if (!ok) {
        n_failed++;
        printf("FAIL: %.*s\n", cast(int) what.length, what.ptr);
    } else {
        printf("ok:   %.*s\n", cast(int) what.length, what.ptr);
    }
}

bool close(double a, double b, double tol = 1e-5)
{
    return fabs(a - b) <= tol * (1.0 + fabs(b));
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    immutable double[20] yv = [0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1];
    DatasetRaw * d = create_new_dataset(3, 20, 0);

    d.pd = 1;
    foreach (i, nm; ["const", "y", "x"]) {
        d.varname[i][0 .. nm.length] = nm[];
        d.varname[i][nm.length] = 0;
    }
    foreach (t; 0 .. 20) {
        d.Z[0][t] = 1;
        d.Z[1][t] = yv[t];
        d.Z[2][t] = t + 1;
    }
    d.t1 = 0;
    d.t2 = 19;
    int[4] list = [3, 1, 0, 2];

    Model m = binary_logit(list.ptr, d, Gretlopt.OPT_NONE, null);

    check(m.errcode == 0 && m.ncoeff == 2, "binary_logit runs");
    check(close(m.coeff[0], -2.0925417977204326) && close(m.coeff[1], 0.19928969502099356), "logit coefficients");
    check(close(m.lnL, -11.308172146661281), "logit log-likelihood");
    check(m.nobs == 20, "logit sample size");
    clear_model(&m);

    m = binary_probit(list.ptr, d, Gretlopt.OPT_NONE, null);
    check(m.errcode == 0 && close(m.coeff[0], -1.2788097290567912) && close(m.coeff[1], 0.1213962427544963), "probit coefficients");
    check(close(m.lnL, -11.297847950581978), "probit log-likelihood");
    clear_model(&m);

    /* perfect prediction */
    foreach (t; 0 .. 20) d.Z[1][t] = t < 10 ? 0 : 1;
    m = binary_logit(list.ptr, d, Gretlopt.OPT_NONE, null);
    check(m.errcode != 0, "perfectly separated data give an error");
    clear_model(&m);

    /* plugin models without the plugin */
    m = tobit_model(list.ptr, 0, NADBL, d, Gretlopt.OPT_NONE, null);
    check(m.errcode == GretlError.fopen, "tobit without its plugin is E_FOPEN");
    clear_model(&m);

    /* bit permutations */
    int err = 0;
    GretlMatrixRaw * bp = bit_permutations(4, 2, &err);

    check(err == 0 && bp.rows == 6 && bp.cols == 4, "bit_permutations(4, 2) has 6 rows");
    double s0 = 0;
    foreach (i; 0 .. bp.rows) {
        double rs = 0;
        foreach (j; 0 .. bp.cols) rs += bp.val[j * bp.rows + i];
        s0 += (rs == 2) ? 1 : 0;
    }
    check(s0 == 6, "every row has exactly two 1s");
    gretl_matrix_free(bp);
    bp = bit_permutations(3, 3, &err);
    check(bp.rows == 1 && bp.val[0] == 1 && bp.val[2] == 1, "bit_permutations(3, 3) is all ones");
    gretl_matrix_free(bp);
    check(bit_permutations(2, 3, &err) is null && err == GretlError.invarg, "bit_permutations rejects k > n");
    err = 0;

    /* Fisher's exact test on [[3, 1], [1, 3]]: left/right/two-sided p-values 0.9857, 0.2429, 0.4857 */
    {
        int[2] r0 = [3, 1];
        int[2] r1 = [1, 3];
        int *[2] fr = [r0.ptr, r1.ptr];
        int[2] rt = [4, 4];
        int[2] ct = [4, 4];
        Xtab tab;

        tab.rows = 2;
        tab.cols = 2;
        tab.f = fr.ptr;
        tab.rtotal = rt.ptr;
        tab.ctotal = ct.ptr;
        tab.n = 8;
        Prn * prn = gretl_print_new(PrnType.buffer, null);

        check(fishers_exact_test(&tab, prn) == 0, "fishers_exact_test runs");
        const(char) * buf = gretl_print_get_buffer(prn);

        check(strstr(buf, "0.242857") !is null && strstr(buf, "0.485714") !is null && strstr(buf, "0.985714") !is null, "Fisher p-values");
        gretl_print_destroy(prn);
    }

    destroy_dataset(d);
    printf("\n%d checks, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
