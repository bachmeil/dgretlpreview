/* Port of kalman.c: the Kalman filter and smoother (the "kalman" bundle of the hansl language): the filter in its legacy,
   sequential (a la KFAS) and de Jong (exact diffuse) forms, state and disturbance smoothing, simulation, the user-level bundle
   interface and its XML serialization, and the ARMA log-likelihood filter used by arma.c.

   The enumerations of kalman.c become KFlags, KVarType, KSmoType and KCodeType (members lower-cased, prefixes dropped), the two
   anonymous enumerations of matrix symbols KMat and KOut, the macros that inspect a Kalman struct are small functions, and
   the structs stepinfo_, seq_info_, dejong_info_ and kalman_ are StepInfo, SeqInfo, DejongInfo and Kalman. The type Kalman
   here is the complete struct (in C it is private to this file); gretl.kalmanh declares it opaque for the modules that only
   pass pointers. The functions that other modules declare extern (C) (kalman_bundle_filter and so on) are extern (C) here, so
   that those declarations bind to them.

   The file ends in C with the inclusion of kalman_legacy.c (the smoothers for the standard filter, anderson_moore_smooth and
   legacy_koopman_smooth); its contents are at the end of this file. */

module gretl.kalman;

import gretl.gretl_matrix;
import gretl.gretl_matrixh;
import gretl.libgretlh;

/* the user-function interface (gretl_func.c, not ported): used to run the user's functions that update time-varying matrices */
extern (C) {
    Ufunc * get_user_function_by_name(const(char) * name);
    Fncall * fncall_new(Ufunc * fun, int preserve);
    int gretl_function_exec(Fncall * call, int rtype, DatasetRaw * dset, void * ret, Prn * prn);
    int push_anon_function_arg(Fncall * fc, GretlType type, void * value);
    char ** gretl_function_retrieve_code(Ufunc * u, int * n);
    const(char) * fn_param_name(const(Ufunc) * u, int i);
}

import gretl.libgretlh : _;

import gretl.geneval0 : Fncall;
import gretl.gretl_bfgs : numgrad_in_progress;
import gretl.gretl_bundle : gretl_bundle_donate_data, gretl_bundle_get_matrix, gretl_bundle_get_private_data, gretl_bundle_set_matrix, kalman_bundle_new;
import gretl.gretl_errors : gretl_errmsg_set, gretl_errmsg_sprintf;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_matrix : clear_gretl_matrix_err, get_gretl_matrix_err, gretl_LU_solve, gretl_identity_matrix_new, gretl_invert_symmetric_matrix2, gretl_invert_triangular_matrix, gretl_is_identity_matrix, gretl_is_zero_matrix, gretl_matrix_I_minus, gretl_matrix_add_transpose_to, gretl_matrix_alloc, gretl_matrix_block_destroy, gretl_matrix_block_new, gretl_matrix_cholesky_decomp, gretl_matrix_copy, gretl_matrix_copy_values, gretl_matrix_copy_values_shaped, gretl_matrix_extract_matrix, gretl_matrix_free, gretl_matrix_inscribe_I, gretl_matrix_inscribe_matrix, gretl_matrix_is_symmetric, gretl_matrix_moore_penrose, gretl_matrix_multiply, gretl_matrix_multiply_by_scalar, gretl_matrix_multiply_mod, gretl_matrix_print, gretl_matrix_psd_root, gretl_matrix_qform, gretl_matrix_realloc, gretl_matrix_reuse, gretl_matrix_transpose, gretl_matrix_zero, gretl_null_matrix_new, gretl_scalar_qform, gretl_square_matrix_transpose, gretl_zero_matrix_new;
import gretl.gretl_matrixh : gretl_column_vector_alloc, gretl_matrix_cols, gretl_matrix_get, gretl_matrix_rows, gretl_matrix_set, gretl_vector_get, gretl_vector_get_length, gretl_vector_set;
import gretl.gretl_prn : pprintf, pputs;
import gretl.gretl_typemap : gretl_type_get_name;
import gretl.gretl_xml : gretl_finite_scalar_serialize, gretl_matrix_serialize, gretl_string_serialize, gretl_xml_get_matrix, gretl_xml_get_prop_as_double, gretl_xml_get_prop_as_string;
import gretl.libgretlh : LN_2_PI, VNAMELEN;
import gretl.libset : gretl_iteration_depth;
import gretl.matrix_extra : gretl_matrix_from_scalar;
import gretl.missingh : NADBL, na;
import gretl.strutils : gretl_strdup, strings_array_add, strings_array_new;

import gretl.geneval0 : Ufunc;
import gretl.gretl_matrix : GretlMatrixBlock;
import gretl.gretl_matrixh : GretlMatrixMod, GretlMatrixRaw;
import gretl.gretl_prn : Prn;
import gretl.libgretlh : GretlBundle, GretlType, Gretlopt;

import gretl.xmlh;
private alias XUC = const(xmlChar) *;

import core.stdc.math : ceil, exp, fabs, floor, fmod, isinf, isnan, log, pow, sqrt;
import core.stdc.stdio : fprintf, fputs, printf, puts, snprintf, sprintf, sscanf, stderr;
import core.stdc.stdlib : abs, atoi, calloc, free, getenv, malloc, realloc;
import core.stdc.string : memcpy, memmove, memset, strchr, strcmp, strcpy, strlen, strncmp, strstr;

enum K_TINY = 1.0e-12;
enum P0_ZERO = 1;
enum SUPPORT_LEGACY = 1;

private __gshared int kdebug = 0;
private __gshared int trace = 0;
private __gshared int identify = 0;

enum KFlags {
    bundle  = 1 << 0,  /* user-defined filter, inside bundle */
    diffuse = 1 << 1,  /* using diffuse P_{1|0} */
    forward = 1 << 2,  /* running forward filtering pass */
    sim     = 1 << 3,  /* running simulation */
    check   = 1 << 4,  /* checking user-defined matrices */
    arma_ll = 1 << 5,  /* filtering for ARMA estimation */
    extra   = 1 << 6   /* recording extra info (att, Ptt) */
}

/* How exactly are the disturbance variances represented? */
enum KVarType {
    std_var, /* standard, basic variance representation: VS, VY */
    dj_var,  /* as per DeJong, supporting cross correlation */
    dk_var   /* as per Durbin and Koopman, with selection */
}

/* On filtering pass, what are we preparing for? */
enum KSmoType {
    none = 0, /* not preparing for smoothing */
    state,    /* preparing for state smoothing */
    dist      /* preparing for disturbance smoothing */
}

/* Which filtering and smoothing code are we using? */
enum KCodeType {
    legacy, /* legacy Kalman code */
    sequen, /* sequential, a la KFAS */
    dejong  /* de Jong and Lin approach */
}


/* recorder for time-varying T and/or Z */
struct StepInfo {
    GretlMatrixRaw * T; /* (r * r) x N */
    GretlMatrixRaw * ZT; /* (r * n) x N */
}


/* holder for data specific to sequential processing */
struct SeqInfo {
    GretlMatrixRaw * g; /* obs disturbance variance vector */
    GretlMatrixRaw * Z; /* untransposed Z matrix */
    GretlMatrixRaw * Zi; /* row i of Z */
    GretlMatrixRaw * Kti; /* column i of Kt */
    GretlMatrixRaw * Pki; /* P at step i */
    GretlMatrixRaw * Kki; /* K∞ at i */
    GretlMatrixRaw * m1; /* workspace */
    GretlMatrixRaw * mm; /* workspace */
    GretlMatrixRaw * Linv; /* for diagonalization */
    GretlMatrixRaw * Finf; /* recorder for smoothing */
    GretlMatrixRaw * Kinf; /* ditto */
    GretlMatrixRaw * y0; /* original data (when diagonalizing) */
    GretlMatrixRaw * V; /* sequential counterpart of regular V */
    GretlMatrixRaw * F; /* sequential counterpart of regular F */
    int free_g;
    int free_Z;
    int Nd;
}


/* holder for data specific to de Jong (exact) representation */
struct DejongInfo {
    GretlMatrixRaw * Vt; /* for augmented recursion */
    GretlMatrixRaw * At; /* ditto */
    GretlMatrixRaw * Vv; /* ditto */
    GretlMatrixRaw * Qt; /* ditto */
    GretlMatrixRaw * S; /* north-west block of Qt */
    GretlMatrixRaw * s; /* vector embedded in Qt */
    GretlMatrixRaw * A; /* recorder for smoother */
    GretlMatrixRaw * V; /* ditto */
    GretlMatrixRaw * Am; /* ditto (redundant?) */
    double qm; /* qt for t = m */
    double ldS; /* log determinant of S-inverse */
    int Nd; /* size of diffuse recorders */
}

struct Kalman {
    int flags; /* bitflags for recording various options */
    KVarType vartype; /* KVarType.std_var, KVarType.dj_var or KVarType.dk_var */
    KCodeType code; /* algorithms used: KCodeType.legacy, KCodeType.sequen or KCodeType.dejong */
    KSmoType smo_prep; /* preparing for smoothing (of type)? */
    int exact; /* using an exact initial diffuse method? */
    int user_code; /* code selected by user (sequential or dejong) */

    int r; /* rows of a = number of elements in state */
    int n; /* columns of y = number of observables */
    int k; /* columns of B = number of exogenous vars in obs eqn */
    int p; /* length of combined disturbance vector */
    int q; /* length of state disturbance vector */
    int N; /* rows of y = number of observations */
    int okN; /* N - number of missing observations */
    int t; /* current time step, when filtering */
    int d; /* (diffuse) time-step at which standard iterations start */
    int j; /* (sequential) observable at which standard iters start */
    int TI; /* flag for K.T is identity matrix */

    int ifc; /* boolean: obs equation includes an implicit constant? */
    int dj_initted; /* flag for basic de Jong initialization done */

    double qsum; /* \sum_{t=1}^N q_t = v_t' F_t^{-1} v_t */
    double loglik; /* log-likelihood */
    double s2; /* = qsum / k */

    /* continuously updated matrices */
    GretlMatrixRaw * a0; /* r x 1: state vector, before updating */
    GretlMatrixRaw * a1; /* r x 1: state vector, after updating */
    GretlMatrixRaw * P0; /* r x r: MSE matrix, before updating */
    GretlMatrixRaw * P1; /* r x r: MSE matrix, after updating */
    GretlMatrixRaw * vt; /* n x 1: one-step forecast error(s), time t */

    /* for the exact diffuse sequential case */
    GretlMatrixRaw * Pk0; /* P∞,t */
    GretlMatrixRaw * PK; /* P∞, all t */

    /* for the de Jong variant */
    GretlMatrixRaw * Lt; /* r x r */
    GretlMatrixRaw * Jt; /* r x p */
    GretlMatrixRaw * rr; /* r x r */
    GretlMatrixRaw * L; /* recorder */
    GretlMatrixRaw * J; /* recorder */
    GretlMatrixRaw * VS0; /* baseline for detecting change */
    GretlMatrixRaw * VY0; /* ditto */

    /* input data matrices: note that the order matters for various
       functions, including matrix_is_varying()
    */
    GretlMatrixRaw * T; /* r x r: state transition matrix */
    GretlMatrixRaw * BT; /* k x n: B', coeffs on exogenous vars, obs eqn */
    GretlMatrixRaw * ZT; /* r x n: Z', coeffs on state variables, obs eqn */
    GretlMatrixRaw * VS; /* r x r: contemp covariance matrix, state eqn */
    GretlMatrixRaw * VY; /* n x n: contemp covariance matrix, obs eqn */
    GretlMatrixRaw * mu; /* r x 1: constant term in state transition */
    GretlMatrixRaw * y; /* N x n: dependent variable vector (or matrix) */
    GretlMatrixRaw * x; /* N x k: independent variables matrix */
    GretlMatrixRaw * aini; /* r x 1: a_0 */
    GretlMatrixRaw * Pini; /* r x r: P_0 */

    /* user inputs for cross-correlated disturbances */
    GretlMatrixRaw * H; /* r x p: state var factor */
    GretlMatrixRaw * G; /* n x p: obs var factor */
    /* and the cross-matrix itself */
    GretlMatrixRaw * HG; /* r x n */

    /* OR: Durbin-Koopman R and Q (as in statevar = RQR') */
    GretlMatrixRaw * Q; /* q x q, symmetric */
    GretlMatrixRaw * R; /* r x q, selection matrix */
    GretlMatrixRaw * QRT; /* Q*R' */

    /* apparatus for registering time-variation of matrices */
    char * matcall;
    char * varying;

    /* export matrices */
    GretlMatrixRaw * V; /* N x n: forecast errors, all time-steps */
    GretlMatrixRaw * F; /* N x n(n+1)/2: MSE for observables, all time-steps */
    GretlMatrixRaw * A; /* N x r: state vector, all time-steps */
    GretlMatrixRaw * P; /* N x r(r+1)/2: MSE for state, all time-steps */
    GretlMatrixRaw * K; /* N x rn: gain matrix, all time-steps */
    GretlMatrixRaw * LL; /* N x 1: loglikelihood, all time-steps */
    GretlMatrixRaw * U; /* N x r+n: smoothed disturbances */
    GretlMatrixRaw * Vsd; /* std. dev. of smoothed disturbance */

    /* recorder for F-inverse, all time-steps */
    GretlMatrixRaw * Finv;

    /* struct needed only when smoothing in the time-varying case */
    StepInfo * step;

    /* struct needed for sequential approach */
    SeqInfo * sinfo;

    /* struct needed for de Jong exact approach */
    DejongInfo * djinfo;

    /* workspace matrices */
    GretlMatrixBlock * Blk; /* holder for the following */
    GretlMatrixRaw * PZ;
    GretlMatrixRaw * Ft;
    GretlMatrixRaw * iFt;
    GretlMatrixRaw * Kt;
    GretlMatrixRaw * Mt;
    GretlMatrixRaw * Ct;

    GretlBundle * b; /* the bundle of which this struct is a member */
    void * data; /* handle for attaching additional info */
    Prn * prn; /* verbose printer */
}


enum KOut {
    v,
    f,
    a,
    big_p,
    k,
    ll,
}

/* max number of time-varying matrices: T, B, Z, VS, VY, mu */
enum K_N_MATCALLS = 6;

private void set_kalman_running(Kalman * K) { K.flags |= KFlags.forward; }
private int kalman_is_running(const(Kalman) * K) { return K.flags & KFlags.forward; }
private int kalman_simulating(const(Kalman) * K) { return K.flags & KFlags.sim; }
private int kalman_checking(const(Kalman) * K) { return K.flags & KFlags.check; }
private int kalman_diffuse(const(Kalman) * K) { return (K.flags & KFlags.diffuse)? 1 : 0; }
private int kalman_arma_ll(const(Kalman) * K) { return K.flags & KFlags.arma_ll; }
private int kalman_extra(const(Kalman) * K) { return K.flags & KFlags.extra; }
private int kalman_is_bundle(const(Kalman) * K) { return (K.flags & KFlags.bundle)? 1 : 0; }

private bool kalman_sequential(const(Kalman) * K) { return K.code == KCodeType.sequen; }
private bool kalman_dejong(const(Kalman) * K) { return K.code == KCodeType.dejong; }

private bool kalman_djvar(const(Kalman) * K) { return K.vartype == KVarType.dj_var; }
private bool kalman_dkvar(const(Kalman) * K) { return K.vartype == KVarType.dk_var; }

private bool filter_is_varying(const(Kalman) * K) { return K.matcall != null; }

int is_kalman_bundle(GretlBundle * b)
{
    Kalman *K = cast(Kalman *) gretl_bundle_get_private_data(b);

    return K != null;
}

/* symbolic identifiers for input matrices: note that potentially
   time-varying matrices must appear first in the enumeration, and
   the order must match the order of the "input data matrices" in
   the Kalman struct (above).
*/

enum KMat {
    t = 0,
    bt,
    zt,
    vs,
    vy,
    m,
    y,
    x,
    a,
    p,
    r,
    mmax /* sentinel */
}

/* variant of gretl_matrix_copy_values for us when we already
   know that the matrices are non-null and conformable
*/

private void fast_copy_values(GretlMatrixRaw * B, const(GretlMatrixRaw) * A)
{
    memcpy(B.val, A.val, B.rows * B.cols * double.sizeof);
}

private void fast_write_I(GretlMatrixRaw * A)
{
    int i;
    int j;
    int k = 0;

    for (j=0; j<A.cols; j++) {
        for (i=0; i<A.rows; i++) {
            A.val[k++] = (i == j)? 1 : 0;
        }
    }
}

private int matrix_changed(const(GretlMatrixRaw) * A, const(GretlMatrixRaw) * B)
{
    int i;
    int n = A.rows * A.cols;

    for (i=0; i<n; i++) {
	if (A.val[i] != B.val[i]) {
	    return 1;
	}
    }

    return 0;
}

private void set_kalman_stopped(Kalman * K)
{
    K.flags &= ~KFlags.forward;
    K.t = 0;
}

private void kalman_check_env()
{
    char *s1 = getenv("KALMAN_DEBUG");
    char *s2 = getenv("KALMAN_TRACE");

    if (s1 != null) kdebug = atoi(s1);
    if (s2 != null) trace  = atoi(s2);
}

private void free_stepinfo(Kalman * K)
{
    if (K.step != null) {
        gretl_matrix_free(K.step.T);
        gretl_matrix_free(K.step.ZT);
        free(K.step);
        K.step = null;
    }
}

private void free_seq_info(Kalman * K)
{
    if (K.sinfo != null) {
        if (K.sinfo.free_g) {
            gretl_matrix_free(K.sinfo.g);
        }
        if (K.sinfo.free_Z) {
            gretl_matrix_free(K.sinfo.Z);
        }
        gretl_matrix_free(K.sinfo.Zi);
        gretl_matrix_free(K.sinfo.Kti);
        gretl_matrix_free(K.sinfo.Pki);
        gretl_matrix_free(K.sinfo.Kki);
        gretl_matrix_free(K.sinfo.m1);
        gretl_matrix_free(K.sinfo.mm);
        gretl_matrix_free(K.sinfo.Linv);
        gretl_matrix_free(K.sinfo.Finf);
        gretl_matrix_free(K.sinfo.Kinf);
        gretl_matrix_free(K.sinfo.y0);
        gretl_matrix_free(K.sinfo.V);
        gretl_matrix_free(K.sinfo.F);

        free(K.sinfo);
        K.sinfo = null;
    }
}

private void free_dejong_info(Kalman * K)
{
    if (K.djinfo != null) {
        gretl_matrix_free(K.djinfo.Vt);
        gretl_matrix_free(K.djinfo.At);
        gretl_matrix_free(K.djinfo.Vv);
        gretl_matrix_free(K.djinfo.Qt);
        gretl_matrix_free(K.djinfo.S);
        gretl_matrix_free(K.djinfo.s);

        gretl_matrix_free(K.djinfo.A);
        gretl_matrix_free(K.djinfo.V);
        gretl_matrix_free(K.djinfo.Am);

        free(K.djinfo);
        K.djinfo = null;
    }
}

extern (C) void kalman_free(Kalman * K)
{
    if (K == null) {
        return;
    }

    /* internally allocated matrices */
    gretl_matrix_free(K.a0);
    gretl_matrix_free(K.a1);
    gretl_matrix_free(K.P0);
    gretl_matrix_free(K.P1);
    gretl_matrix_free(K.vt);
    gretl_matrix_free(K.Pk0);
    gretl_matrix_free(K.L);
    gretl_matrix_free(K.J);
    gretl_matrix_free(K.Finv);

    /* internally allocated workspace */
    gretl_matrix_block_destroy(K.Blk);

    if (kalman_is_bundle(K)) {
        GretlMatrixRaw ** [] mptr = [
            &K.T, &K.BT, &K.ZT, &K.VS, &K.VY,
            &K.mu, &K.y, &K.x, &K.aini, &K.Pini,
            &K.R
        ];
        int i;

        for (i=0; i<KMat.mmax; i++) {
            gretl_matrix_free(*mptr[i]);
        }

        /* @K also owns these "export" matrices */
        gretl_matrix_free(K.V);
        gretl_matrix_free(K.F);
        gretl_matrix_free(K.A);
        gretl_matrix_free(K.P);
        gretl_matrix_free(K.K);
        gretl_matrix_free(K.LL);
        gretl_matrix_free(K.U);
        gretl_matrix_free(K.Vsd);
    }

    /* time-variation info */
    free(K.matcall);
    free(K.varying);

    if (!kalman_sequential(K)) {
        gretl_matrix_free(K.H);
        gretl_matrix_free(K.G);
        gretl_matrix_free(K.HG);
        gretl_matrix_free(K.Jt);
	gretl_matrix_free(K.VS0);
	gretl_matrix_free(K.VY0);
    }

    if (kalman_dkvar(K)) {
        /* Durbin-Koopman errors info */
        gretl_matrix_free(K.Q);
        gretl_matrix_free(K.QRT);
    }

    if (K.step != null) {
        free_stepinfo(K);
    }
    if (K.sinfo != null) {
        free_seq_info(K);
    }
    if (K.djinfo != null) {
        free_dejong_info(K);
    }

    free(K);
}

private Kalman * kalman_new_empty(int flags)
{
    Kalman * K = cast(Kalman *) malloc(Kalman.sizeof);

    if (K != null) {
        K.exact = 0;
        K.user_code = KCodeType.legacy;
	K.code = KCodeType.legacy;
        K.vartype = KVarType.std_var;
        K.aini = K.Pini = null;
        K.a0 = K.a1 = null;
        K.P0 = K.P1 = null;
        K.Pk0 = null;
        K.PK = null;
        K.vt = null;
        K.Blk = null;
        K.T = K.BT = K.ZT = null;
        K.VS = K.VY = null;
        K.H = K.G = K.HG = K.Jt = null;  /* de Jong variance */
        K.Q = K.R = K.QRT = null; /* Durbin-Koopman variance */
        K.V = K.F = K.A = K.K = K.LL = K.P = null;
        K.Finv = null;
	K.VS0 = K.VY0 = null;
        K.L = K.J = null;
        K.y = K.x = null;
        K.mu = null;
        K.U = null;
        K.Vsd = null;
        K.matcall = null;
        K.varying = null;
        K.step = null;
        K.sinfo = null;
        K.djinfo = null;
        K.flags = flags;
        K.t = 0;
        K.prn = null;
        K.data = null;
        K.b = null;
        K.d = 0;
        K.j = 0;
	K.smo_prep = KSmoType.none;
	K.dj_initted = 0;
    }

    return K;
}

private void kalman_set_dimensions(Kalman * K)
{
    K.r = gretl_matrix_rows(K.T);  /* T.rows defines r */
    K.k = gretl_matrix_rows(K.BT); /* BT.rows defines k */
    K.n = gretl_matrix_cols(K.y);  /* y.cols defines n */

    if (!kalman_simulating(K)) {
        /* y.rows defines N, except when simulating */
        K.N = gretl_matrix_rows(K.y);
    }

    K.okN = K.N;

    /* K.p is non-zero only under cross-correlation; in that case the
       matrix given as 'Q' in Kalman set-up in fact represents H (as
       in v_t = H \varepsilon_t) and it must be r x p, where p is the
       number of elements in the "combined" disturbance vector
       \varepsilon_t.
    */
    K.p = (K.H != null)? gretl_matrix_cols(K.H): 0;

    /* K.q is non-zero only when using the Durbin-Koopman
       representation of the variance of the state disturbance as
       R*Q*R'
    */
    K.q = (K.Q != null)? gretl_matrix_rows(K.Q): 0;
}

private int missing_matrix_error(const(char) * name)
{
    if (name == null) {
        gretl_errmsg_set(_("kalman: a required matrix is missing").ptr);
    } else {
        gretl_errmsg_sprintf(_("kalman: required matrix %s is missing").ptr,
                             name);
    }
    return GretlError.data;
}

private int check_matrix_dims(Kalman * K, const(GretlMatrixRaw) * m, int i)
{
    int r = 0, c = 0, symm = (i == KMat.vs || i == KMat.vy);
    int err = 0;

    if (i == KMat.t || i == KMat.vs || i == KMat.p) {
        r = c = K.r;
    } else if (i == KMat.bt) {
        r = K.k;
        c = K.n;
    } else if (i == KMat.zt) {
        r = K.r;
        c = K.n;
    } else if (i == KMat.vy)  {
        r = c = K.n;
    } else if (i == KMat.a || i == KMat.m) {
        r = K.r;
        c = 1;
    }

    if (m.rows != r || m.cols != c) {
        gretl_errmsg_sprintf(_("kalman: %s is %d x %d, should be %d x %d\n").ptr,
                             kalman_matrix_name(i), m.rows, m.cols, r, c);
        err = GretlError.nonconf;
    } else if (symm && !gretl_matrix_is_symmetric(m)) {
        gretl_errmsg_sprintf(_("kalman: %s is not symmetric\n").ptr,
                kalman_matrix_name(i));
        err = GretlError.nonconf;
    }

    return err;
}

private int maybe_resize_export_matrix(Kalman * K, GretlMatrixRaw * m, int i)
{
    int rows = K.N;
    int cols = 0;
    int err = 0;

    if (i == KOut.v) {
        cols = K.n;
    } else if (i == KOut.f) {
        /* vech(F) per row */
        cols = (K.n * K.n + K.n) / 2;
    } else if (i == KOut.a) {
        cols = K.r;
    } else if (i == KOut.big_p) {
	/* vech(P) per row */
	cols = (K.r * K.r + K.r) / 2;
    } else if (i == KOut.ll) {
        cols = 1;
    } else if (i == KOut.k) {
        cols = K.r * K.n;
    } else {
        err = GretlError.data;
    }

    if (!err && (m.rows != rows || m.cols != cols)) {
        if (trace) {
            printf("resize export matrix %d from %d x %d to %d x %d\n",
                   i, m.rows, m.cols, rows, cols);
        }
        err = gretl_matrix_realloc(m, rows, cols);
        if (!err) {
            gretl_matrix_zero(m);
        }
    }

    return err;
}

private int kalman_check_dimensions(Kalman * K)
{
    int err = 0;

    if (trace) {
        printf("kalman_check_dimensions()\n");
    }

    if (K.r < 1 || K.n < 1 || K.N < 2) {
        /* the state and observation vectors must have at least one
           element, and there must be at least two observations
        */
        err = GretlError.data;
    }

    /* T is mandatory, should be r x r */
    if (!err) {
        err = check_matrix_dims(K, K.T, KMat.t);
    }

    /* Z' is mandatory, should be r x n */
    if (!err) {
        err = check_matrix_dims(K, K.ZT, KMat.zt);
    }

    if (K.H != null) {
        /* de Jong variance representation */
        if (K.H.rows != K.r || K.H.cols != K.p) {
            fprintf(stderr, "H is %d x %d, p=%d, r=%d, n=%d\n",
                    K.H.rows, K.H.cols, K.p, K.r, K.n);
            err = GretlError.nonconf;
        }
        if (!err && K.G != null) {
            if (K.G.rows != K.n || K.G.cols != K.p) {
                fprintf(stderr, "G is %d x %d, p=%d, r=%d, n=%d\n",
                        K.G.rows, K.G.cols, K.p, K.r, K.n);
                err = GretlError.nonconf;
            }
        }
    } else if (K.R != null) {
        /* state disturbances as per Durbin-Koopman */
        if (K.Q == null) {
            err = missing_matrix_error("statevar");
        } else if (K.R.rows != K.r || K.R.cols != K.Q.rows) {
            fprintf(stderr, "R is %d x %d, Q is %d x %d, r=%d\n",
                    K.R.rows, K.R.cols, K.Q.rows, K.Q.cols,
                    K.r);
            err = GretlError.nonconf;
        }
    } else {
        /* VS is mandatory, should be r x r and symmetric */
        if (!err) {
            err = check_matrix_dims(K, K.VS, KMat.vs);
        }
        /* VY should be n x n and symmetric, if present */
        if (!err && K.VY != null) {
            err = check_matrix_dims(K, K.VY, KMat.vy);
        }
    }

    /* initial a should be r x 1, if present */
    if (!err && K.aini != null) {
        err = check_matrix_dims(K, K.aini, KMat.a);
    }

    /* initial P should be r x r, if present */
    if (!err && K.Pini != null) {
        err = check_matrix_dims(K, K.Pini, KMat.p);
    }

    /* B' should be k x n, if present (BT.rows defines k) */
    if (!err) {
        if (K.BT != null) {
            err = check_matrix_dims(K, K.BT, KMat.bt);
        } else if (K.x != null) {
            /* B is null => can't have a non-null x */
            err = GretlError.nonconf;
        }
    }

    /* mu should be r x 1, if present */
    if (!err && K.mu != null) {
        err = check_matrix_dims(K, K.mu, KMat.m);
    }

    if (err) {
        goto bailout;
    }

    K.ifc = 0;

    /* x should have N rows to match y; and it should have either k or k - 1
       columns (the latter case indicating an implicit const) */
    if (K.x != null) {
        if (K.x.rows < K.N) {
            fprintf(stderr, "kalman: %s has %d rows, should have %d\n",
                    kalman_matrix_name(KMat.x), K.x.rows, K.N);
            return GretlError.nonconf;
        } else if (K.x.cols != K.k && K.x.cols != K.k - 1) {
            fprintf(stderr, "kalman: %s has %d columns, should have %d or %d\n",
                    kalman_matrix_name(KMat.x), K.x.cols, K.k, K.k - 1);
            return GretlError.nonconf;
        } else if (K.x.cols == K.k - 1) {
            /* register the implicit const */
            K.ifc = 1;
        }
    } else if (K.k == 1) {
        /* B has one row but there's no x => implicit const */
        K.ifc = 1;
    } else if (K.k > 1) {
        /* B has more than one row but there's no x => error */
        return missing_matrix_error("obsxmat");
    }

    /* Below we have the "export" matrices for providing results. If
       these are present but not sized correctly we'll try to fix them
       up -- but note that they are not used in a simulation run.
    */

    if (kalman_simulating(K)) {
        return err;
    }

    /* big V should be N x n */
    if (K.V != null) {
        err = maybe_resize_export_matrix(K, K.V, KOut.v);
    }

    /* big F should be N x (n*n+n) / 2 */
    if (!err && K.F != null) {
        err = maybe_resize_export_matrix(K, K.F, KOut.f);
    }

    /* big A should be N x r */
    if (!err && K.A != null) {
        err = maybe_resize_export_matrix(K, K.A, KOut.a);
    }

    /* big P should be N x nr, or N x r*r */
    if (!err && K.P != null) {
        err = maybe_resize_export_matrix(K, K.P, KOut.big_p);
    }

    /* K (gain) should be N x (r*n) */
    if (!err && K.K != null) {
        err = maybe_resize_export_matrix(K, K.K, KOut.k);
    }

    /* LL should be N x 1 */
    if (!err && K.LL != null) {
        err = maybe_resize_export_matrix(K, K.LL, KOut.ll);
    }

 bailout:

    if (err) {
        fprintf(stderr, "kalman_check_dimensions: err = %d\n", err);
    }

    return err;
}

private int kalman_init(Kalman * K)
{
    int err = 0;

    if (trace) {
        printf("kalman_init(): aini %p, Pini %p\n",
               cast(void *) K.aini, cast(void *) K.Pini);
    }

    K.qsum = NADBL;
    K.loglik = NADBL;
    K.s2 = NADBL;

    clear_gretl_matrix_err();

    if (K.aini != null) {
        K.a0 = gretl_matrix_copy(K.aini);
        K.a1 = gretl_matrix_copy(K.aini);
    } else {
        K.a0 = gretl_zero_matrix_new(K.r, 1);
        K.a1 = gretl_zero_matrix_new(K.r, 1);
    }

    if (K.Pini != null) {
        K.P0 = gretl_matrix_copy(K.Pini);
        K.P1 = gretl_matrix_copy(K.Pini);
    } else {
        K.P0 = gretl_zero_matrix_new(K.r, K.r);
        K.P1 = gretl_zero_matrix_new(K.r, K.r);
    }

    /* forecast error vector, per observation */
    K.vt = gretl_matrix_alloc(K.n, 1);

    err = get_gretl_matrix_err();
    if (err) {
        return err;
    }

    K.Blk = gretl_matrix_block_new(&K.PZ,  K.r, K.n, /* P*Z */
                                    &K.Ft,  K.n, K.n, /* (Z'*P*Z + R)^{-1} */
                                    &K.iFt, K.n, K.n, /* Ft-inverse */
                                    &K.Kt,  K.r, K.n, /* gain at t */
                                    &K.Mt,  K.r, K.n, /* intermediate term */
                                    &K.Ct,  K.r, K.r, /* intermediate term */
                                    &K.Lt,  K.r, K.r, /* ditto */
                                    &K.rr,  K.r, K.r, /* ditto */
                                    null);

    if (K.Blk == null) {
        err = GretlError.alloc;
    }

    if (!err && !kalman_is_bundle(K)) {
        /* in the "user-bundle" case we do this later */
        err = set_initial_statevar(K);
    }

    return err;
}

/* The following section includes functions that support the "plain C"
   Kalman API, as opposed to the bundle-based interface that subserves
   userland state-space functionality. We use this API in gretl's arma
   plugin in the special case of ARIMA with non-zero order of
   integration along with missing values. It's also possible that
   third-party users of libgretl might wish to use it.
*/

/**
 * kalman_new:
 * @a: r x 1 initial state vector.
 * @P: r x r initial precision matrix.
 * @T: r x r state transition matrix.
 * @BT: n x k matrix of coefficients on exogenous variables in the
 * observation equation, transposed.
 * @ZT: n x r matrix of coefficients on the state variables in the
 * observation equation, transposed.
 * @VS: r x r contemporaneous covariance matrix for the errors in the
 * state equation.
 * @VY: n x n contemporaneous covariance matrix for the errors in the
 * observation equation (or null if this is not applicable).
 * @y: T x n matrix of observable variable(s).
 * @x: T x k matrix of exogenous variable(s).  May be null if there
 * are no exogenous variables, or if there's only a constant.
 * @mu: r x 1 vector of constants in the state transition, or null.
 * @V: T x n matrix in which to record forecast errors (or null if
 * this is not required).
 * @err: location to receive error code.
 *
 * Allocates and initializes a Kalman struct, which can subsequently
 * be used for filtering and smoothing.
 *
 * Returns: pointer to allocated struct, or null on failure, in
 * which case @err will receive a non-zero code.
 */

Kalman * kalman_new(GretlMatrixRaw * a, GretlMatrixRaw * P, GretlMatrixRaw * T, GretlMatrixRaw * BT, GretlMatrixRaw * ZT, GretlMatrixRaw * VS, GretlMatrixRaw * VY, GretlMatrixRaw * y, GretlMatrixRaw * x, GretlMatrixRaw * mu, GretlMatrixRaw * V, int * err)
{
    Kalman * K;

    *err = 0;

    kalman_check_env();

    if (trace) {
        printf("kalman_new()\n");
    }

    if (y == null || T == null || ZT == null || VS == null) {
        fprintf(stderr, "kalman_new: y=%p, T=%p, ZT=%p, VS=%p\n",
                cast(void *) y, cast(void *) T, cast(void *) ZT, cast(void *) VS);
        *err = missing_matrix_error(null);
        return null;
    }

    K = kalman_new_empty(0);
    if (K == null) {
        *err = GretlError.alloc;
        return null;
    }

    /* use pointers for input matrices, don't copy */
    K.T = T;
    K.BT = BT;
    K.ZT = ZT;
    K.VS = VS;
    K.VY = VY;
    K.y = y;
    K.x = x;
    K.aini = a;
    K.Pini = P;
    K.mu = mu;

    /* output, but again use external pointer */
    K.V = V;

    kalman_set_dimensions(K);

    *err = kalman_check_dimensions(K);
    if (*err) {
        free(K);
        return null;
    }

    *err = kalman_init(K);

    if (*err) {
        kalman_free(K);
        K = null;
    } else {
        gretl_matrix_zero(K.vt);
    }

    return K;
}

/**
 * kalman_get_loglik:
 * @K: pointer to Kalman struct.
 *
 * Retrieves the log-likelhood calculated via a filtering run.
 *
 * Returns: ll value, or #NADBL on failure.
 */

double kalman_get_loglik(const(Kalman) * K)
{
    return K.loglik;
}

/**
 * kalman_get_arma_variance:
 * @K: pointer to Kalman struct.
 *
 * Retrieves the estimated variance for an ARMA model
 * estimated using the Kalman filter.
 *
 * Returns: sigma-squared value, or #NADBL on failure.
 */

double kalman_get_arma_variance(const(Kalman) * K)
{
    if (na(K.qsum)) {
        return NADBL;
    } else {
        return K.qsum / K.okN;
    }
}

/**
 * kalman_set_initial_state_vector:
 * @K: pointer to Kalman struct.
 * @a: vector of values to set.
 *
 * Resets the initial value of the state vector in a Kalman
 * struct, using the values from @a.  See also kalman_new().
 *
 * Returns: 0 on success, non-zero on error.
 */

int kalman_set_initial_state_vector(Kalman * K, const(GretlMatrixRaw) * a)
{
    return gretl_matrix_copy_values(K.a0, a);
}

/**
 * kalman_set_initial_MSE_matrix:
 * @K: pointer to Kalman struct.
 * @P: matrix of values to set.
 *
 * Resets the initial value of the MSE matrix in a Kalman
 * struct, using the values from @P.  See also kalman_new().
 *
 * Returns: 0 on success, non-zero on error.
 */

int kalman_set_initial_MSE_matrix(Kalman * K, const(GretlMatrixRaw) * P)
{
    return gretl_matrix_copy_values(K.P0, P);
}

void kalman_attach_data(Kalman * K, void * data)
{
    if (K != null) {
        K.data = data;
    }
}

void * kalman_get_data(const(Kalman) * K)
{
    return (K != null)? cast(void *) K.data : null;
}

void kalman_attach_printer(Kalman * K, Prn * prn)
{
    if (K != null) {
        K.prn = prn;
    }
}

Prn * kalman_get_printer(const(Kalman) * K)
{
    return (K != null)? cast(Prn *) K.prn : null;
}

void kalman_set_arma_ll(Kalman * K)
{
    K.flags |= KFlags.arma_ll;
}

/* end of functions dedicated to non-bundle C API */

private int statemat_out_of_bounds(Kalman * K)
{
    GretlMatrixRaw * evals;
    double r;
    double c;
    double x;
    int i;
    int err = 0;

    evals = legacy_general_matrix_eigenvals(K.T, &err);

    for (i=0; i<evals.rows && !err; i++) {
        r = gretl_matrix_get(evals, i, 0);
        c = gretl_matrix_get(evals, i, 1);
        x = sqrt(r*r + c*c);
        if (x >= 1.0) {
            fprintf(stderr, "T: modulus of eigenvalue %d = %g\n", i+1, x);
            err = GretlError.singular;
        }
    }

    gretl_matrix_free(evals);

    return err;
}

enum kappa = 1.0e7; /* 1.0e6? */

private int diffuse_Pini(Kalman * K)
{
    if (trace) {
        printf("diffuse_Pini(), exact = %d\n", K.exact);
    }

    if (K.exact) {
        if (kalman_dejong(K)) {
            fast_copy_values(K.P0, K.VS); /* H*H' */
        } else {
            /* sequential case */
            gretl_matrix_zero(K.P0);
            /* initialize P∞ */
            if (K.Pk0 == null) {
                K.Pk0 = gretl_identity_matrix_new(K.r);
                if (K.Pk0 == null) {
                    return GretlError.alloc;
                }
            } else {
                fast_write_I(K.Pk0);
            }
        }
    } else {
        /* old-style, using big kappa */
        int i;

        gretl_matrix_zero(K.P0);
        for (i=0; i<K.r; i++) {
            gretl_matrix_set(K.P0, i, i, kappa);
        }
    }

    return 0;
}

private int hamilton_Pini(Kalman * K)
{
    GretlMatrixRaw * Svar;
    GretlMatrixRaw * vQ;
    int r2 = K.r * K.r;
    int err = 0;

    if (trace) {
        printf("hamilton_Pini()\n");
    }

    Svar = gretl_matrix_alloc(r2, r2);
    vQ = gretl_column_vector_alloc(r2);

    if (Svar == null || vQ == null) {
        gretl_matrix_free(Svar);
        gretl_matrix_free(vQ);
        return GretlError.alloc;
    }

    gretl_matrix_kronecker_product(K.T, K.T, Svar);
    gretl_matrix_I_minus(Svar);
    gretl_matrix_vectorize(vQ, K.VS);

    err = gretl_LU_solve(Svar, vQ);
    if (err) {
        /* failed: are some of the eigenvalues out of bounds? */
        err = statemat_out_of_bounds(K);
        if (err == GretlError.singular) {
            err = diffuse_Pini(K);
            K.flags |= KFlags.diffuse;
        }
	if (kdebug) {
	    printf("Hamilton-style P0 failed, set KFlags.diffuse\n");
	}
    } else {
	if (kdebug) {
	    printf("Hamilton-style P0 succeeded\n");
	}
        gretl_matrix_unvectorize(K.P0, vQ);
    }

    gretl_matrix_free(Svar);
    gretl_matrix_free(vQ);

    return err;
}

/* If the user has not given an initial value for P_{1|0}, compute
   this automatically as per Hamilton, ch 13, p. 378.  This works only
   if the eigenvalues of K.T lie inside the unit circle.  Failing
   that, or if the --diffuse option is given for the user Kalman
   filter, we apply a diffuse initialization.
*/

private int set_initial_statevar(Kalman * K)
{
    if (K.Pini != null) {
        fast_copy_values(K.P0, K.Pini);
        return 0;
    } else if (K.flags & KFlags.diffuse) {
        return diffuse_Pini(K);
    } else {
        return hamilton_Pini(K);
    }
}

/* Write the vech of @src into row @t of @targ */

private void record_to_vech(GretlMatrixRaw * targ, const(GretlMatrixRaw) * src, int n, int t)
{
    int i;
    int j;
    int k = 0;
    double x;

    for (j=0; j<n; j++) {
        for (i=j; i<n; i++) {
            x = gretl_matrix_get(src, i, j);
            gretl_matrix_set(targ, t, k++, x);
        }
    }
}

/* Write the vec of @src into row @t of @targ */

private void record_to_vec(GretlMatrixRaw * targ, const(GretlMatrixRaw) * src, int t)
{
    int i;

    for (i=0; i<targ.cols; i++) {
        gretl_matrix_set(targ, t, i, src.val[i]);
    }
}

/* Write the vec of @src into column @t of @targ */

private void record_to_col(GretlMatrixRaw * targ, const(GretlMatrixRaw) * src, int t)
{
    size_t sz = targ.rows * double.sizeof;

    memcpy(targ.val + t * targ.rows, src.val, sz);
}

/* copy from vector @src into row @t of @targ */

private void record_to_row(GretlMatrixRaw * targ, const(GretlMatrixRaw) * src, int t)
{
    int i;

    for (i=0; i<targ.cols; i++) {
        gretl_matrix_set(targ, t, i, src.val[i]);
    }
}

/* copy from vector @src into row @t of @targ,
   starting at column offset @j in @targ */

private void record_to_row_offset(GretlMatrixRaw * targ, const(GretlMatrixRaw) * src, int t, int j)
{
    int i, n = gretl_vector_get_length(src);

    for (i=0; i<n; i++) {
        gretl_matrix_set(targ, t, i+j, src.val[i]);
    }
}

private void set_row_to_value(GretlMatrixRaw * targ, int t, double x)
{
    int i;

    for (i=0; i<targ.cols; i++) {
        gretl_matrix_set(targ, t, i, x);
    }
}

private int check_DK_dims(const(GretlMatrixRaw) * Q, const(GretlMatrixRaw) * R)
{
    int m = Q.rows; /* Q: m x m */
    int r = R.rows; /* R: r x m, with m <= r */

    if (Q.cols != m) {
	gretl_errmsg_set("Q must be square");
	return GretlError.invarg;
    } else if (R.cols != m) {
	gretl_errmsg_set("cols(R) must equal rows(Q)");
	return GretlError.invarg;
    } else if (m > r) {
	gretl_errmsg_set("rows(Q) must not exceed rows(R))");
	return GretlError.invarg;
    } else {
	return 0;
    }
}

/* Supports the hansl function ksetup() for creating a Kalman bundle */

extern (C) Kalman * kalman_new_minimal(GretlMatrixRaw * * M, int * copy, int nmat, int dkvar, int * err)
{
    GretlMatrixRaw ** [5] targ;
    Kalman * K;
    int i;

    if (trace) {
        printf("kalman_new_minimal (dkvar = %d)\n", dkvar);
    }

    *err = 0;

    if (M[0] == null || M[1] == null || M[2] == null || M[3] == null) {
        fprintf(stderr, "kalman_new_minimal: nmat=%d, y=%p, Z=%p, T=%p, Q=%p\n",
                nmat, cast(void *) M[0], cast(void *) M[1], cast(void *) M[2], cast(void *) M[3]);
        *err = missing_matrix_error(null);
        return null;
    }

    K = kalman_new_empty(KFlags.bundle);
    if (K == null) {
        *err = GretlError.alloc;
        return null;
    }

    kalman_check_env();

    targ[0] = &K.y;
    targ[1] = &K.ZT;
    targ[2] = &K.T;

    if (nmat == 5) {
        if (dkvar) {
	    *err = check_DK_dims(M[3], M[4]);
	    if (*err) {
		goto bailout;
	    }
            K.vartype = KVarType.dk_var;
            targ[3] = &K.Q;
            targ[4] = &K.R;
        } else {
            K.vartype = KVarType.dj_var;
            targ[3] = &K.H;
            targ[4] = &K.G;
        }
    } else {
        targ[3] = &K.VS; /* statevar */
    }

    for (i=0; i<nmat; i++) {
        if (copy[i]) {
            *targ[i] = gretl_matrix_copy(M[i]);
        } else {
            *targ[i] = M[i];
        }
    }

    kalman_set_dimensions(K);

    if (K.vartype != KVarType.std_var) {
        *err = kalman_revise_variance(K);
    }
    if (!*err) {
        *err = kalman_check_dimensions(K);
    }
    if (!*err) {
        *err = kalman_init(K);
    }

 bailout:

    if (*err) {
        kalman_free(K);
        K = null;
    } else {
        gretl_matrix_zero(K.vt);
    }

    return K;
}

private int matrix_is_varying(Kalman * K, int i)
{
    if (K.matcall != null) {
        if (K.varying == null) {
            check_for_matrix_updates(K, null);
        }
        if (K.varying != null) {
            return K.varying[i];
        }
    }

    return 0;
}

enum KUpdate {
    initial, /* initialization of matrices */
    step     /* refreshing matrices per time-step */
}

/* Ensure that we have suitable matrices into which to write
   VS, VY and HG'
*/

private int ensure_cross_covariance_matrices(Kalman * K)
{
    int[] r = [K.H.rows, K.G.rows, K.H.rows];
    int[] c = [K.H.rows, K.G.rows, K.G.rows];
    GretlMatrixRaw ** [] V = [&K.VS, &K.VY, &K.HG];
    GretlMatrixRaw ** targ;
    int i;
    int err = 0;

    for (i=0; i<3 && !err; i++) {
        targ = V[i];
        if (*targ == null) {
            *targ = gretl_matrix_alloc(r[i], c[i]);
        } else {
            GretlMatrixRaw * Vi = *targ;

            if (Vi.rows != r[i] || Vi.cols != c[i]) {
                gretl_matrix_realloc(Vi, r[i], c[i]);
            }
        }
        if (*targ == null) {
            err = GretlError.alloc;
        }
    }

    return err;
}

private int ensure_DK_covariance_matrices(Kalman * K)
{
    int[] r = [K.r, K.q];
    int[] c = [K.r, K.r];
    GretlMatrixRaw ** [] V = [&K.VS, &K.QRT];
    GretlMatrixRaw ** targ;
    int i;
    int err = 0;

    for (i=0; i<2 && !err; i++) {
        targ = V[i];
        if (*targ == null) {
            *targ = gretl_matrix_alloc(r[i], c[i]);
        } else {
            GretlMatrixRaw * Vi = *targ;

            if (Vi.rows != r[i] || Vi.cols != c[i]) {
                gretl_matrix_realloc(Vi, r[i], c[i]);
            }
        }
        if (*targ == null) {
            err = GretlError.alloc;
        }
    }

    return err;
}

/* After reading H and G from the user, either at (re-)initialization or
   at a given time-step when either of these matrices is time-varying,
   record the user input in K.H and K.G and form the 'real' VS and VY.
   But note that in the time-step case it may be that only one of VS, VY
   needs to be treated in this way (if only one is time-varying, only
   one will have been redefined via a function call).
*/

private int kalman_update_crossinfo(Kalman * K, int mode)
{
    int err = 0;

    /* Note that H and G may be needed as such for simulation */

    if (mode == KUpdate.initial) {
        err = ensure_cross_covariance_matrices(K);
        if (err) {
            return err;
        }
    }

    if (mode == KUpdate.initial || matrix_is_varying(K, KMat.vs)) {
        /* (re)create VS using modified H */
        err = gretl_matrix_multiply_mod(K.H, GretlMatrixMod.none,
                                        K.H, GretlMatrixMod.transpose,
                                        K.VS, GretlMatrixMod.none);
    }

    if (!err && (mode == KUpdate.initial || matrix_is_varying(K, KMat.vy))) {
        /* (re)create VY using modified G */
        err = gretl_matrix_multiply_mod(K.G, GretlMatrixMod.none,
                                        K.G, GretlMatrixMod.transpose,
                                        K.VY, GretlMatrixMod.none);
    }

    if (!err && (mode == KUpdate.initial || matrix_is_varying(K, KMat.vs) ||
                 matrix_is_varying(K, KMat.vy))) {
        /* (re)create HG' using modified H and/or G */
        err = gretl_matrix_multiply_mod(K.H, GretlMatrixMod.none,
                                        K.G, GretlMatrixMod.transpose,
                                        K.HG, GretlMatrixMod.none);
    }

    return err;
}

private int kalman_update_dkvar(Kalman * K, int mode)
{
    int err = 0;

    if (mode == KUpdate.initial) {
        err = ensure_DK_covariance_matrices(K);
        if (err) {
            return err;
        }
    }

    if (mode == KUpdate.initial || matrix_is_varying(K, KMat.vs)) {
        /* (re)create VS and QR' using modified Q */
        gretl_matrix_qform(K.R, GretlMatrixMod.none, K.Q,
                           K.VS, GretlMatrixMod.none);
        err = gretl_matrix_multiply_mod(K.Q, GretlMatrixMod.none,
                                        K.R, GretlMatrixMod.transpose,
                                        K.QRT, GretlMatrixMod.none);
     }

    return err;
}

/* kalman_revise_DJ_variance: the user has actually given H in place
   of VS and G in place of VY; so we have to form VS, VY, and HG'
   (cross-correlated disturbances).

   This function is called in the course of initial set-up of a
   filter, and also when the Kalman matrices are being re-checked at
   the start of filtering, smoothing or simulation.
*/

private int kalman_revise_DJ_variance(Kalman * K)
{
    int err = 0;

    if (K.H == null || K.G == null) {
        fprintf(stderr, "K.H %p, K.G %p\n", cast(void *) K.H, cast(void *) K.G);
        return missing_matrix_error("'statevar' or 'obsvar'");
    }

    err = kalman_update_crossinfo(K, KUpdate.initial);

    if (err) {
        fprintf(stderr, "kalman_revise_cross_variance: err = %d\n", err);
    }

    return err;
}

private int kalman_revise_DK_variance(Kalman * K)
{
    int err = 0;

    if (K.Q == null || K.R == null) {
        return missing_matrix_error("'statevar' or 'statsel'");
    }

    err = kalman_update_dkvar(K, KUpdate.initial);

    if (err) {
        fprintf(stderr, "kalman_revise_DK_variance: err = %d\n", err);
    }

    return err;
}

private int kalman_revise_variance(Kalman * K)
{
    if (K.vartype == KVarType.dj_var) {
        return kalman_revise_DJ_variance(K);
    } else {
        return kalman_revise_DK_variance(K);
    }
}

private void kalman_print_state(Kalman * K)
{
   int j;

    if (K.t > 4) return;

    fprintf(stderr, "t = %d:\n", K.t);

    for (j=0; j<K.n; j++) {
        fprintf(stderr, "y[%d] = %0.10g, err[%d] = %0.10g\n", j,
                gretl_matrix_get(K.y, K.t, j),
                j, gretl_vector_get(K.vt, j));
    }

    gretl_matrix_print(K.Kt, "K.Kt");
    gretl_matrix_print(K.a0, "K.a0");
    gretl_matrix_print(K.P0, "K.P0");
}

/* On filtering: record the state and/or its variance, as
   wanted.
*/

private int kalman_record_state(Kalman * K)
{
    int err = 0;

    if (K.A != null) {
        record_to_row(K.A, K.a0, K.t);
    }
    if (K.P != null) {
        record_to_vech(K.P, K.P0, K.r, K.t);
    }
    if (K.exact && K.djinfo != null && K.djinfo.A != null) {
	/* record "augmented" quantities */
        record_to_col(K.djinfo.A, K.djinfo.At, K.t);
        record_to_col(K.djinfo.V, K.djinfo.Vt, K.t);
    }

    return err;
}

private double get_xjt(Kalman * K, int j)
{
    if (K.ifc) {
        return (j == 0)? 1.0 : gretl_matrix_get(K.x, K.t, j-1);
    } else {
        return gretl_matrix_get(K.x, K.t, j);
    }
}

/* Read from the appropriate row of x (N x k) and multiply by B' to
   form B'x_t.  Note: the flag K.ifc is used to indicate that the
   observation equation has an implicit constant, with an entry in
   the B matrix (the first) but no explicit entry in the x matrix.

   The case where x is null and B an n-vector (implicit constant)
   is also handled.

   There's no need to store B'x_t: we either want to add this to,
   or subtract it from, the n-vector @targ (the operator being
   indicated by the @mod argument).

   Returns: the number of missing values encountered.
*/

private int kalman_do_Bx(Kalman * K, GretlMatrixRaw * targ, GretlMatrixMod mod)
{
    double xjt;
    double bji;
    double bxi;
    int i;
    int j;
    int missvals = 0;

    for (i=0; i<K.n; i++) {
        if (K.x == null) {
            /* the implicit constant case */
            bxi = K.BT.val[i];
        } else {
            bxi = 0;
            for (j=0; j<K.k; j++) {
                xjt = get_xjt(K, j);
                bji = gretl_matrix_get(K.BT, j, i);
                if (bji != 0) {
                    /* here we'll implicitly take 0 * NA as 0 */
                    if (na(xjt)) {
                        missvals++;
                    } else {
                        bxi += xjt * bji;
                    }
                }
            }
        }
        if (mod == GretlMatrixMod.decrement) {
            targ.val[i] -= bxi;
        } else {
            targ.val[i] += bxi;
        }
    }

    return missvals;
}

/* version of the above for use with the sequential approach */

private double kalman_Bx_uni(Kalman * K, int i)
{
    double xjt;
    double bji;
    double bxi;
    int j;

    if (K.x == null) {
        bxi = K.BT.val[i];
    } else {
        bxi = 0;
        for (j=0; j<K.k; j++) {
            xjt = get_xjt(K, j);
            bji = gretl_matrix_get(K.BT, j, i);
            if (bji != 0 && na(xjt)) {
                return NADBL;
            } else {
                bxi += xjt * bji;
            }
        }
    }

    return bxi;
}

/* Compute the one-step ahead forecast error:

   v_t = y_t - B_t*x_t - Z_t*a_t

   Returns the number of valid observables at time step t,
   or 0 if the forecast error cannot be computed due to
   missing values.
*/

private int compute_forecast_error(Kalman * K)
{
    int i;

    /* initialize v_t to y_t */
    for (i=0; i<K.n; i++) {
        K.vt.val[i] = gretl_matrix_get(K.y, K.t, i);
    }

    if (K.BT != null) {
        /* subtract effect of exogenous terms, if any */
        kalman_do_Bx(K, K.vt, GretlMatrixMod.decrement);
    }

    /* check for any missing values in v_t */
    for (i=0; i<K.n; i++) {
        if (na(K.vt.val[i])) {
            return 0;
        }
    }

    /* subtract contribution from state */
    gretl_matrix_multiply_mod(K.ZT,  GretlMatrixMod.transpose,
                              K.a0, GretlMatrixMod.none,
                              K.vt, GretlMatrixMod.decrement);

    return K.n;
}

/* Given a unified function to update one or more of the potentially
   time-varying matrices, try to figure out which matrix or matrices
   are actually modified by this function.
*/

private int check_for_matrix_updates(Kalman * K, Ufunc * uf)
{
    char **lines;
    int nlines = 0;

    if (K.varying != null) {
        free(K.varying);
        K.varying = null;
    }

    if (uf == null) {
        uf = get_user_function_by_name(K.matcall);
        if (uf == null) {
            gretl_errmsg_sprintf(_("Couldn't find function '%s'").ptr, K.matcall);
            return GretlError.data;
        }
    }

    K.varying = cast(char *) calloc(K_N_MATCALLS, 1);

    lines = gretl_function_retrieve_code(uf, &nlines);

    if (lines != null) {
        const char *bname = fn_param_name(uf, 0);
        char[VNAMELEN+1] test;
        const(char) * s;
        int n = cast(int) strlen(bname) + 1;
        int i;
        int j;

        sprintf(test.ptr, "%s.", bname);
        for (i=0; i<nlines; i++) {
            if (!strncmp(lines[i], test.ptr, n)) {
                for (j=KMat.t; j<=KMat.m; j++) {
                    s = kalman_matrix_name(j);
                    if (!strncmp(lines[i] + n, s, strlen(s))) {
			if (kdebug) {
			    fprintf(stderr, "matrix %s is varying\n", s);
			}
                        K.varying[j] = 1;
                        break;
                    }
                }
            }
        }
        free(lines);
    }

    return 0;
}

/* Function to update any time-varying matrices, for use with a Kalman
   bundle. Bypasses the regular "genr" apparatus, passing the attached
   bundle directly to the given user function after it has been found
   by name.
*/

private int kalman_update_matrices(Kalman * K, Prn * prn)
{
    Ufunc * uf;
    Fncall * fc;
    int err = 0;

    uf = get_user_function_by_name(K.matcall);
    if (uf == null) {
        gretl_errmsg_sprintf(_("Couldn't find function '%s'").ptr, K.matcall);
        return GretlError.data;
    }

    if (K.varying == null) {
        check_for_matrix_updates(K, uf);
    }

    fc = fncall_new(uf, 0);
    err = push_anon_function_arg(fc, GretlType.bundle_ref, K.b);

    if (!err) {
        err = gretl_function_exec(fc, GretlType.none, null, null, prn);
    }
    if (err) {
        fprintf(stderr, "kalman_update_matrices: call='%s', err=%d\n",
                K.matcall, err);
    }

    return err;
}

/* On getting a new K.ZT in case of time-variation, redo any
   transformations that may be needed.
*/

private int adjust_sequential_Z(Kalman * K)
{
    int err = 0;

    if (K.r > 1 || K.n > 1) {
        if (kdebug > 1) {
            fprintf(stderr, "doing adjust_sequential_Z\n");
        }
        /* transposition needed */
        if (K.sinfo.Z == null) {
            fprintf(stderr, "K.sinfo.Z not allocated!\n");
            err = GretlError.data;
        } else {
            err = gretl_matrix_transpose(K.sinfo.Z, K.ZT);
        }
        if (!err && K.sinfo.Linv != null) {
            /* diagonalization needed (FIXME conditionality) */
            gretl_matrix_multiply_mod(K.sinfo.Linv, GretlMatrixMod.none,
                                      K.ZT, GretlMatrixMod.transpose,
                                      K.sinfo.Z, GretlMatrixMod.none);
        }
    }

    return err;
}

/* If we have any time-varying coefficient matrices, refresh these for
   the current time step. This is called on a forward filtering pass.
*/

private int kalman_refresh_matrices(Kalman * K, Prn * prn)
{
    GretlMatrixRaw ** [] mptr = [
        &K.T, &K.BT, &K.ZT, &K.VS, &K.VY, &K.mu
    ];
    int cross_update = 0;
    int dkvar_update = 0;
    int i;
    int err = 0;

    if (kalman_djvar(K)) {
        mptr[3] = &K.H;
        mptr[4] = &K.G;
    } else if (kalman_dkvar(K)) {
        mptr[3] = &K.Q;
        mptr[4] = &K.R;
    }

    if (K.matcall != null) {
        err = kalman_update_matrices(K, prn);
    }

    for (i=0; i<K_N_MATCALLS && !err; i++) {
        if (matrix_is_varying(K, i)) {
            if (kalman_djvar(K) && (i == KMat.vs || i == KMat.vy)) {
                /* handle revised H and/or G */
                cross_update = 1;
            } else if (kalman_dkvar(K) && i == KMat.vs) {
                /* handle revised Q and/or R */
                dkvar_update = 1;
            } else {
                err = check_matrix_dims(K, *mptr[i], i);
            }
            if (err) {
                fprintf(stderr, "kalman_refresh_matrices: err = %d at t = %d\n",
                        err, K.t);
            } else if (i == KMat.zt && kalman_sequential(K)) {
                err = adjust_sequential_Z(K);
            } else if (i == KMat.t) {
                K.TI = gretl_is_identity_matrix(K.T);
            }
        }
    }

    if (!err && K.step != null) {
        /* keep a record of T and/or Z' at the given time step */
        if (K.step.T != null) {
            record_to_col(K.step.T, K.T, K.t);
        }
        if (K.step.ZT != null) {
            if (K.sinfo != null) {
                record_to_col(K.step.ZT, K.sinfo.Z, K.t);
            } else {
                record_to_col(K.step.ZT, K.ZT, K.t);
            }
        }
    }

    if (!err) {
        if (cross_update) {
            /* cross-correlated case */
            err = kalman_update_crossinfo(K, KUpdate.step);
        } else if (dkvar_update) {
            /* Durbin-Koopman case */
            err = kalman_update_dkvar(K, KUpdate.step);
        }
    }

    return err;
}

/* Handle missing values in y_t or x_t, for the case of univariate y_t
   or when all y_t values are missing. This is called from
   kfilter_standard() and kfilter_dejong(), but not when the
   sequential method is being used (kfilter_sequential).
*/

private void handle_missing_obs(Kalman * K)
{
    /* state update: a1 = T*a0 */
    gretl_matrix_multiply(K.T, K.a0, K.a1);
    /* handle stconst if present */
    if (K.mu != null) {
        gretl_matrix_add_to(K.a1, K.mu);
    }
    fast_copy_values(K.a0, K.a1);

    /* var(state) update: P1 = T*P0*T' + VS (C = 0) */
    fast_copy_values(K.P1, K.VS);
    gretl_matrix_qform(K.T, GretlMatrixMod.none, K.P0,
                       K.P1, GretlMatrixMod.cumulate);
    fast_copy_values(K.P0, K.P1);

    /* record stuff if wanted */
    if (K.F != null) {
        if (K.smo_prep) {
            set_row_to_value(K.F, K.t, 0.0);
        } else {
            set_row_to_value(K.F, K.t, NADBL); /* 0 ? */
        }
    }
    if (K.LL != null) {
        gretl_vector_set(K.LL, K.t, 0.0); /* ? */
    }
    if (K.K != null) {
        set_row_to_value(K.K, K.t, 0.0);
    }
    if (K.V != null) {
        if (K.smo_prep) {
            set_row_to_value(K.V, K.t, 0.0);
        } else {
            set_row_to_value(K.V, K.t, NADBL); /* 0 ? */
        }
    }
    if (K.L != null) {
	record_to_col(K.L, K.T, K.t);
    }
    if (K.J != null) {
	record_to_col(K.J, K.H, K.t);
    }
    /* more to be done if we're in the de Jong diffuse phase? */
}

private int ensure_Finv(Kalman * K)
{
    GretlMatrixRaw * Fi = K.Finv;
    int c = (K.n*K.n + K.n) / 2;

    if (Fi != null && (Fi.rows != K.N || Fi.cols != c)) {
        gretl_matrix_free(Fi);
    }

    K.Finv = gretl_matrix_alloc(K.N, c);

    return K.Finv == null ? GretlError.alloc : 0;
}

/**
 * kfilter_standard:
 * @K: pointer to Kalman struct.
 * @prn: printing apparatus (or null).
 *
 * Generates a series of one-step ahead forecasts for y, based on
 * information in the Kalman struct @K.
 *
 * Returns: 0 on success, non-zero on error.
 */

int kfilter_standard(Kalman * K, Prn * prn)
{
    double ll0 = K.n * LN_2_PI;
    double sumldet = 0;
    int nt;
    int err = 0;

    if (trace) {
        printf("kfilter_standard()\n");
    }
    if (kdebug > 1) {
        fprintf(stderr, "\n*** kfilter_standard: N=%d, n=%d ***\n",
                K.N, K.n);
    }

    if (K.smo_prep) {
        err = ensure_Finv(K);
        if (err) {
            return err;
        }
    }

    K.qsum = K.loglik = 0.0;
    K.s2 = NADBL;
    K.okN = K.N;
    set_kalman_running(K);

    for (K.t = 0; K.t < K.N && !err; K.t += 1) {
        double llt = NADBL;
        double ldet = 0;
        double qt = 0;

        if (K.A != null || K.P != null) {
            kalman_record_state(K);
        }

        if (filter_is_varying(K)) {
            /* we have time-varying coefficients */
            err = kalman_refresh_matrices(K, prn);
            if (err) {
                K.loglik = NADBL;
                break;
            }
        }

        if (kdebug > 1) {
            kalman_print_state(K);
        }

        /* calculate v_t, checking for missing values */
        nt = compute_forecast_error(K);
        if (nt == 0) {
            /* skip this observation */
            K.okN -= 1;
            handle_missing_obs(K);
            continue;
        }

        /* calculate F_t = ZPZ' [+ VY] */
        if (K.VY != null) {
            fast_copy_values(K.Ft, K.VY);
            gretl_matrix_qform(K.ZT, GretlMatrixMod.transpose, K.P0,
                               K.Ft, GretlMatrixMod.cumulate);
        } else {
            gretl_matrix_qform(K.ZT, GretlMatrixMod.transpose, K.P0,
                               K.Ft, GretlMatrixMod.none);
        }

        /* calculate M_t = TPZ' [+ HG'] */
        gretl_matrix_multiply(K.P0, K.ZT, K.PZ);
        gretl_matrix_multiply(K.T, K.PZ, K.Mt);
        if (K.HG != null) {
            gretl_matrix_add_to(K.Mt, K.HG);
        }

        /* standard Kalman procedure */
        fast_copy_values(K.iFt, K.Ft);
        err = gretl_invert_symmetric_matrix2(K.iFt, &ldet);
        if (err) {
            fprintf(stderr, "kfilter_standard: failed to invert Ft\n");
        } else {
            qt = gretl_scalar_qform(K.vt, K.iFt, &err);
        }

        if (K.F != null) {
            /* we're recording F_t for all t */
            if (K.smo_prep) {
                /* record inverse */
                record_to_vech(K.Finv, K.iFt, nt, K.t);
            } else {
                /* record F_t itself */
                record_to_vech(K.F, K.Ft, nt, K.t);
            }
        }

        /* determine and record loglikelihood */
        if (err) {
            K.loglik = NADBL;
            break;
        } else {
            llt = -0.5 * (ll0 + ldet + qt);
            if (na(llt)) {
                K.loglik = NADBL;
                break;
            }
            K.loglik += llt;
            K.qsum += qt;
            sumldet += ldet;
        }
        if (K.LL != null) {
            gretl_vector_set(K.LL, K.t, llt);
        }

        /* Calculate gain K_t = M_t F_t^{-1}, and C matrix */
        gretl_matrix_multiply(K.Mt, K.iFt, K.Kt);
        gretl_matrix_multiply_mod(K.Kt, GretlMatrixMod.none,
                                  K.Mt, GretlMatrixMod.transpose,
                                  K.Ct, GretlMatrixMod.none);
        if (!err && K.K != null) {
            /* record the gain */
            record_to_vec(K.K, K.Kt, K.t);
        }

        /* update state: a1 = T a0 + K_t v_t [+ mu] */
        gretl_matrix_multiply(K.T, K.a0, K.a1);
        if (K.mu != null) {
            gretl_matrix_add_to(K.a1, K.mu);
        }
        gretl_matrix_multiply_mod(K.Kt, GretlMatrixMod.none,
                                  K.vt, GretlMatrixMod.none,
                                  K.a1, GretlMatrixMod.cumulate);
        fast_copy_values(K.a0, K.a1);

        /* update var(state): P1 = TPT' + VS - C */
        fast_copy_values(K.P1, K.VS);
        gretl_matrix_qform(K.T, GretlMatrixMod.none,
                           K.P0, K.P1, GretlMatrixMod.cumulate);
        gretl_matrix_subtract_from(K.P1, K.Ct);
        fast_copy_values(K.P0, K.P1);

        /* record forecast errors if wanted */
        if (!err && K.V != null) {
            record_to_row(K.V, K.vt, K.t);
        }
    }

    set_kalman_stopped(K);

    if (na(K.loglik)) {
        err = GretlError.nan;
    } else if (kalman_arma_ll(K)) {
        double ll1 = 1.0 + LN_2_PI + log(K.qsum / K.okN);

        K.loglik = -0.5 * (K.okN * ll1 + sumldet);
    } else {
	int d = kalman_diffuse(K) ? K.r : 0;

        K.s2 = K.qsum / (K.n * K.okN - d);
    }

    if (kdebug) {
        fprintf(stderr, "kfilter_standard: err=%d, ll=%#0.8g\n",
		err, K.loglik);
    }

    return err;
}

struct K_input_mat {
    int sym;
    const(char) * name;
}

/* mapping to names used in setting elements of Kalman bundle,
   either via ksetup() or manually
*/

private immutable K_input_mat[] K_input_mats = [
    K_input_mat(KMat.y, "obsy"),
    K_input_mat(KMat.zt, "obsymat"),
    K_input_mat(KMat.x, "obsx"),
    K_input_mat(KMat.bt, "obsxmat"),
    K_input_mat(KMat.vy, "obsvar"),
    K_input_mat(KMat.t, "statemat"),
    K_input_mat(KMat.vs, "statevar"),
    K_input_mat(KMat.m, "stconst"),
    K_input_mat(KMat.a, "inistate"),
    K_input_mat(KMat.p, "inivar"),
    K_input_mat(KMat.r, "statesel")
];

/* indices of input matrices that cannot be set via ksetup() */

private immutable int[] extra_input_mats = [
    KMat.bt, /* obsxmat */
    KMat.vy, /* obsvar */
    KMat.m,  /* stconst */
    KMat.x,  /* obsx */
    KMat.a,  /* inistate */
    KMat.p,  /* inivar */
    KMat.r   /* statesel */
];

private immutable int n_extra_mats = cast(int) extra_input_mats.length;

private int obsy_check(Kalman * K)
{
    if (K.y == null) {
        return missing_matrix_error("obsy");
    } else if (K.y.rows != K.N || K.y.cols != K.n) {
        fprintf(stderr, "obsy_check: K.y should be %d x %d, is %d x %d\n",
                K.N, K.n, gretl_matrix_rows(K.y),
                gretl_matrix_cols(K.y));
        return GretlError.nonconf;
    } else {
        return 0;
    }
}

private int kalman_bundle_recheck_matrices(Kalman * K, Prn * prn)
{
    int err = 0;

    K.flags |= KFlags.check;

    if (filter_is_varying(K)) {
        err = kalman_update_matrices(K, prn);
    }

    K.flags ^= KFlags.check;

    if (!err && (K.ZT == null || K.T == null || K.VS == null)) {
        fprintf(stderr, "kalman_bundle_kalman_recheck_matrices: Z=%p, T=%p, Q=%p\n",
                K.ZT, K.T, K.VS);
        err = missing_matrix_error(null);
    }

    if (err) {
        return err;
    }

    /* maybe redundant? */
    kalman_set_dimensions(K);

    if (gretl_matrix_rows(K.T) != K.r ||
        gretl_matrix_rows(K.BT) != K.k) {
        err = GretlError.nonconf;
    } else if (!kalman_simulating(K)) {
        err = obsy_check(K);
    }

    if (!err && K.vartype != KVarType.std_var) {
        err = kalman_revise_variance(K);
    }

    if (!err) {
        err = kalman_check_dimensions(K);
    }

    if (!err) {
        if (K.aini != null) {
            gretl_matrix_copy_values(K.a0, K.aini);
        } else {
            gretl_matrix_zero(K.a0);
        }
        err = set_initial_statevar(K);
    }

    return err;
}

private const(char) * kalman_matrix_name(int sym)
{
    int i;

    for (i=0; i<KMat.mmax; i++) {
        if (K_input_mats[i].sym == sym) {
            return K_input_mats[i].name;
        }
    }

    /* failed */
    return "matrix";
}

private int kalman_ensure_output_matrices(Kalman * K)
{
    int err = 0;

    if (trace) {
        printf("kalman_ensure_output_matrices()\n");
    }

    if (K.V == null) {
        K.V = gretl_null_matrix_new();
    }
    if (K.F == null) {
        K.F = gretl_null_matrix_new();
    }
    if (K.A == null) {
        K.A = gretl_null_matrix_new();
    }
    if (K.P == null) {
        K.P = gretl_null_matrix_new();
    }
    if (K.K == null) {
        K.K = gretl_null_matrix_new();
    }
    if (K.LL == null) {
        K.LL = gretl_null_matrix_new();
    }
    if (K.V == null || K.F == null || K.A == null ||
        K.P == null || K.K == null || K.LL == null) {
        err = GretlError.alloc;
    }

    return err;
}

private int ensure_sequential_info(Kalman * K)
{
    if (K.sinfo == null) {
        return kalman_add_seq_info(K);
    } else {
        return 0;
    }
}

private int ensure_dejong_info(Kalman * K)
{
    if (K.djinfo == null) {
        return kalman_add_dejong_info(K);
    } else {
        return 0;
    }
}

/* Kalman_filter() is called (only) by kalman_bundle_filter(), which
   implements the userspace function kfilter(). It performs various
   checks and (conditional) allocations then hands off to one of
   kfilter_standard(), kfilter_sequential() or kfilter_dejong() to do
   the actual filtering.
*/

private int kalman_filter(Kalman * K, Prn * prn, int * errp)
{
    int err = kalman_ensure_output_matrices(K);

    if (trace) {
        printf("kalman_filter()\n");
    }

    if (kalman_diffuse(K) && K.exact) {
        if (!kalman_dejong(K)) {
            K.code = KCodeType.sequen;
        }
    }

    if (!err) {
        gretl_matrix_zero(K.vt);
        /* includes setting of K.P0 */
        err = kalman_bundle_recheck_matrices(K, prn);
    }

    if (!err) {
        if (kalman_sequential(K)) {
            err = ensure_sequential_info(K);
        } else if (K.exact) {
            /* only if diffuse? */
            err = ensure_dejong_info(K);
        }
    }

    if (K.exact && kalman_sequential(K)) {
        /* in case we're re-running the filter */
        fast_write_I(K.Pk0);
    }

    if (!err) {
        if (kalman_sequential(K)) {
            err = kfilter_sequential(K, prn);
        } else if (kalman_dejong(K)) {
	    err = kfilter_dejong(K, prn);
	} else {
	    err = kfilter_standard(K, prn);
	}
    }

    if (err != GretlError.nan) {
        *errp = err;
    } else {
        /* we'll flag GretlError.nan with a return value of 1 but
           won't count it as a 'true' error
        */
        *errp = 0;
    }

    return err;
}

/* Implements the userland kfilter() function */

extern (C) int kalman_bundle_filter(GretlBundle * b, Prn * prn, int * errp)
{
    Kalman *K = cast(Kalman *) gretl_bundle_get_private_data(b);

    if (trace) {
        printf("kalman_bundle_filter()\n");
    }

    K.code = cast(KCodeType) K.user_code;
    K.b = b; /* attach bundle pointer */

    return kalman_filter(K, prn, errp);
}

/* write row @t of matrix @src into matrix @targ */

private void load_from_row(GretlMatrixRaw * targ, const(GretlMatrixRaw) * src, int t)
{
    int j;

    for (j=0; j<src.cols; j++) {
        targ.val[j] = gretl_matrix_get(src, t, j);
    }
}

/* Copy row @t from @src into @targ; or add row @t of @src to @targ;
   or subtract row @t of @src from @targ.  We allow the possibility
   that the length of vector @targ is less than the number of columns
   in @src, but not the converse.
*/

private int vector_from_row_mod(GretlMatrixRaw * targ, const(GretlMatrixRaw) * src, int t, GretlMatrixMod mod)
{
    int i, n = gretl_vector_get_length(targ);
    double x;

    if (n > src.cols) {
        fprintf(stderr, "vector_from_row: targ length = %d, but src " ~
                "has %d columns\n", n, src.cols);
        return GretlError.nonconf;
    }

    for (i=0; i<n; i++) {
        x = gretl_matrix_get(src, t, i);
        if (mod == GretlMatrixMod.cumulate) {
            targ.val[i] += x;
        } else if (mod == GretlMatrixMod.decrement) {
            targ.val[i] -= x;
        } else {
            targ.val[i] = x;
        }
    }

    return 0;
}

/* Similar to vector_from_row_offset() except that a column offset,
   @j, is supported for the reading of a row from @src, and we don't
   support the @mod option.
*/

private int vector_from_row_offset(GretlMatrixRaw * targ, const(GretlMatrixRaw) * src, int t, int j)
{
    int i, n = gretl_vector_get_length(targ);
    double x;

    if (n > src.cols - j) {
        return GretlError.nonconf;
    }

    for (i=0; i<n; i++) {
        x = gretl_matrix_get(src, t, i + j);
        targ.val[i] = x;
    }

    return 0;
}

/* Row @t of @src represents the vech of an n x n matrix: extract the
   row and apply the inverse operation of vech to reconstitute the
   matrix in @targ.
*/

private void load_from_vech(GretlMatrixRaw * targ, const(GretlMatrixRaw) * src, int n, int t)
{
    int i;
    int j;
    int k = 0;
    double x;

    for (j=0; j<n; j++) {
        for (i=j; i<n; i++) {
            x = gretl_matrix_get(src, t, k++);
            gretl_matrix_set(targ, i, j, x);
            if (i != j) {
                gretl_matrix_set(targ, j, i, x);
            }
        }
    }
}

/* Row @t of @src represents the vec of a certain matrix: extract the
   row and reconstitute the matrix in @targ.
*/

private int load_from_vec(GretlMatrixRaw * targ, const(GretlMatrixRaw) * src, int t)
{
    int i;
    int k = targ.rows * targ.cols;

    for (i=0; i<k; i++) {
        targ.val[i] = gretl_matrix_get(src, t, i);
    }

    return 0;
}

/* Column @t of @src represents the vec of a certain matrix: extract the
   column and reconstitute the matrix in @targ.
*/

private void load_from_col(GretlMatrixRaw * targ, const(GretlMatrixRaw) * src, int t)
{
    size_t sz = src.rows * double.sizeof;

    memcpy(targ.val, src.val + t * src.rows, sz);
}

private int retrieve_Tt(Kalman * K)
{
    if (K.step == null || K.step.T == null) {
        return GretlError.data;
    } else {
        load_from_col(cast(GretlMatrixRaw *) K.T, K.step.T, K.t);
        return 0;
    }
}

private int retrieve_Zt(Kalman * K)
{
    if (K.step == null || K.step.ZT == null) {
        return GretlError.data;
    } else {
        if (K.sinfo != null) {
            load_from_col(cast(GretlMatrixRaw *) K.sinfo.Z, K.step.ZT, K.t);
        } else {
            load_from_col(cast(GretlMatrixRaw *) K.ZT, K.step.ZT, K.t);
        }
        return 0;
    }
}

private int load_filter_data(Kalman * K, int no_collapse, int * err)
{
    int missing = 0;

    /* FIXME maybe we should have a byte array of length K.N
       to identify missing observations (constructed on the
       filtering step if needed). This could speed things
       up below,
    */

    /* load the forecast error into K.vt */
    load_from_row(K.vt, K.V, K.t);

    /* state: get T_t and/or Z_t if need be */
    if (matrix_is_varying(K, KMat.t)) {
        *err = retrieve_Tt(K);
    }
    if (!*err && matrix_is_varying(K, KMat.zt)) {
        *err = retrieve_Zt(K);
    }
    if (*err) {
        return 0;
    }

    /* load the state and its MSE */
    load_from_row(K.a0, K.A, K.t);
    if (no_collapse) {
        /* FIXME? looks weird but seems to help */
        if (K.t == 0) {
            gretl_matrix_zero(K.P0);
        } else {
            load_from_vech(K.P0, K.P, K.r, K.t-1);
        }
    } else {
        load_from_vech(K.P0, K.P, K.r, K.t);
    }

    /* load the gain and F^{-1} */
    load_from_vec(K.Kt, K.K, K.t);
    load_from_vech(K.iFt, K.Finv, K.n, K.t);

    /* heuristic for missing observable (see note above) */
    missing = gretl_is_zero_matrix(K.Kt);

    return missing ? 0 : K.n;
}

/* kalman_ensure_stepinfo() is called in two cases:

  (a) If we're doing sequential filtering, K.ZT is time-varying, and
   K.n > 1, we need a record of K.ZT per time step for building a
   proper F matrix (pevar).

  (b) If we're doing smoothing for a system that has time-varying
   coefficients in K.T or K.ZT we record the vec of these matrices
   for each time-step on the forward pass.
*/

private int kalman_ensure_stepinfo(Kalman * K, int smoothing)
{
    int alloc_T = 0;
    int alloc_Z = 0;
    int err = 0;

    if (K.step != null) {
        /* already started */
        if (smoothing && matrix_is_varying(K, KMat.t)) {
            if (K.step.T != null &&
                K.step.T.rows == K.r * K.r &&
                K.step.T.cols == K.N) {
                ; /* T already OK */
            } else {
                gretl_matrix_free(K.step.T);
                alloc_T = 1;
            }
        }
        if (matrix_is_varying(K, KMat.zt)) {
            if (K.step.ZT != null &&
                K.step.ZT.rows == K.r * K.n &&
                K.step.ZT.cols == K.N) {
                ; /* ZT already OK */
            } else {
                gretl_matrix_free(K.step.ZT);
                alloc_Z = 1;
            }
        }
    } else {
        /* starting from scratch */
        K.step = cast(StepInfo *) malloc(StepInfo.sizeof);
        if (K.step == null) {
            return GretlError.alloc;
        }
        K.step.T = K.step.ZT = null;
        alloc_T = smoothing && matrix_is_varying(K, KMat.t);
        alloc_Z = matrix_is_varying(K, KMat.zt);
    }

    if (alloc_T) {
        K.step.T = gretl_matrix_alloc(K.r * K.r, K.N);
        if (K.step.T == null) {
            err = GretlError.alloc;
        }
    }
    if (!err && alloc_Z) {
        K.step.ZT = gretl_matrix_alloc(K.r * K.n, K.N);
        if (K.step.ZT == null) {
            err = GretlError.alloc;
        }
    }

    if (err) {
        free_stepinfo(K);
    }

    return err;
}

private int kalman_add_seq_info(Kalman * K)
{
    int err = 0;

    K.sinfo = cast(SeqInfo *) calloc(1, SeqInfo.sizeof);
    if (K.sinfo == null) {
        return GretlError.alloc;
    }

    K.sinfo.Zi  = gretl_matrix_alloc(1, K.r);
    K.sinfo.Kti = gretl_matrix_alloc(K.r, 1);
    K.sinfo.Pki = gretl_matrix_alloc(K.r, K.r);
    K.sinfo.Kki = gretl_matrix_alloc(K.r, 1);
    K.sinfo.m1  = gretl_matrix_alloc(K.r, 1);
    K.sinfo.mm  = gretl_matrix_alloc(K.r, K.r);

    if (K.n > 1) {
        K.sinfo.V = gretl_matrix_alloc(K.N, K.n);
        K.sinfo.F = gretl_matrix_alloc(K.N, K.n);
    } else {
        K.sinfo.V = null;
        K.sinfo.F = null;
    }

    K.sinfo.free_Z = 0;
    K.sinfo.free_g = 0;
    K.sinfo.Linv = null;
    K.sinfo.y0 = null;

    if (K.r > 1 || K.n > 1) {
        /* we'll need untransposed Z */
        K.sinfo.Z = gretl_zero_matrix_new(K.n, K.r);
        K.sinfo.free_Z = 1;
    } else {
        K.sinfo.Z = K.ZT;
    }

    if (K.VY != null) {
        if (K.n > 1) {
            K.sinfo.g = gretl_zero_matrix_new(K.n, 1);
            K.sinfo.free_g = 1;
            /* we may need this too */
            K.sinfo.Linv = gretl_matrix_copy(K.VY);
        } else {
            K.sinfo.g = K.VY;
        }
    } else {
        K.sinfo.g = null;
    }

    K.sinfo.Finf = null;
    K.sinfo.Kinf = null;
    K.sinfo.Nd = 0;

    /* FIXME proper error checking above */

    return err;
}

private int kalman_add_dejong_info(Kalman * K)
{
    int err = 0;

    K.djinfo = cast(DejongInfo *) calloc(1, DejongInfo.sizeof);
    if (K.djinfo == null) {
        return GretlError.alloc;
    }

    K.djinfo.Vt = gretl_zero_matrix_new(K.n, K.r);
    K.djinfo.At = gretl_zero_matrix_new(K.r, K.r);
    K.djinfo.Vv = gretl_matrix_alloc(K.n, K.r + 1);
    K.djinfo.Qt = gretl_zero_matrix_new(K.r + 1, K.r + 1);
    K.djinfo.S  = gretl_matrix_alloc(K.r, K.r);
    K.djinfo.s  = gretl_matrix_alloc(K.r, 1);
    K.djinfo.A  = null;
    K.djinfo.V  = null;
    K.djinfo.Am = null;
    K.djinfo.qm = 0;

    return err;
}

private int real_kalman_smooth(Kalman * K, int dist, Prn * prn)
{
    int DKstyle = 0;
    int err = 0;

    if (trace) {
        printf("real_kalman_smooth(), K.code %d, dist %d\n", K.code, dist);
    }

    if (dist == 2) {
	DKstyle = 1;
    }

    err = kalman_ensure_output_matrices(K);
    if (err) {
        return err;
    }

    if (matrix_is_varying(K, KMat.t) || matrix_is_varying(K, KMat.zt)) {
        /* add recorder for T_t and/or Z_t */
        err = kalman_ensure_stepinfo(K, 1);
        if (err) {
            goto bailout;
        }
    }

    if (K.exact && kalman_sequential(K)) {
        int Nd = matrix_is_varying(K, KMat.zt) ? K.N : 4 * K.r;
        int rows = K.r * K.r;

        K.PK = gretl_zero_matrix_new(rows, Nd);
    }

    if (!err) {
        err = kalman_bundle_recheck_matrices(K, prn);
    }

    if (!err) {
        if (kalman_sequential(K)) {
            err = ensure_sequential_info(K);
        } else if (K.exact) {
            err = ensure_dejong_info(K);
        }
    }

    if (!err) {
        /* prior forward pass */
	K.smo_prep = dist ? KSmoType.dist : KSmoType.state;
        if (kalman_sequential(K)) {
            err = kfilter_sequential(K, null);
        } else if (kalman_dejong(K)) {
            err = kfilter_dejong(K, null);
        } else {
	    err = kfilter_standard(K, null);
	}
        K.smo_prep = KSmoType.none;
    }

    K.t = 0;

    if (!err) {
        if (kalman_sequential(K)) {
            err = ksmooth_sequential(K, dist);
	} else if (kalman_dejong(K)) {
            if (dist) {
                err = dist_smooth_dejong(K, DKstyle);
	    } else {
		err = state_smooth_dejong(K);
	    }
	} else {
            if (dist) {
		err = legacy_koopman_smooth(K, DKstyle);
	    } else {
		err = anderson_moore_smooth(K);
	    }
        }
    }

    /* free "special case" storage */
    if (K.PK != null) {
        gretl_matrix_free(K.PK);
        K.PK = null;
    }

 bailout:

    /* trash the "StepInfo" storage */
    free_stepinfo(K);

    return err;
}

/* For use with userspace bundle-based API: implements
   the ksmooth() and kdsmooth() functions. The @dist parameter
   works as follows:

   0 . state smoothing
   1 . disturbance smoothing
   2 . disturbance smoothing with "DKstyle" variance
*/

extern (C) int kalman_bundle_smooth(GretlBundle * b, int dist, Prn * prn)
{
    Kalman *K = cast(Kalman *) gretl_bundle_get_private_data(b);

    if (trace) {
        printf("kalman_bundle_smooth(), dist = %d\n", dist);
    }

    K.b = b; /* attach bundle pointer */

    return real_kalman_smooth(K, dist, prn);
}

/* For use with "plain C" API. Maybe support disturbance smoothing via
   @opt? At present only state smoothing is supported.
*/

GretlMatrixRaw * kalman_smooth(Kalman * K, Gretlopt opt, Prn * prn, int * err)
{
    GretlMatrixRaw * S = null;

    if (trace) {
        printf("kalman_smooth()\n");
    }

    *err = real_kalman_smooth(K, 0, prn);

    if (K.A != null && K.P != null) {
        int r = K.A.rows;
        int c = K.A.cols;
        size_t Asize = r * c * double.sizeof;
        size_t Psize = 0;

        if (opt & Gretlopt.OPT_M) {
            /* include MSE of state */
            c += K.P.cols;
            Psize = r * K.P.cols * double.sizeof;
        }
        if (r > 0 && c > 0) {
            S = gretl_matrix_alloc(r, c);
            if (S == null) {
                *err = GretlError.alloc;
            }
        } else {
            *err = GretlError.data;
        }
        if (S != null) {
            memcpy(S.val, K.A.val, Asize);
            if (opt & Gretlopt.OPT_M) {
                memcpy(S.val + r * K.A.cols, K.P.val, Psize);
            }
        }
    }

    return S;
}

/* note: it's OK for @S to be null (if the simulated state is not
   wanted), so watch out for that!
*/

private int kalman_simulate(Kalman * K, const(GretlMatrixRaw) * U, const(GretlMatrixRaw) * Sim0, GretlMatrixRaw * Y, GretlMatrixRaw * S, Prn * prn)
{
    GretlMatrixRaw * yt = null;
    GretlMatrixRaw * et = null;
    int obs_offset = 0;
    int obsdist = 0;
    int tmin = 0;
    int err = 0;

    if (Sim0 != null) {
        /* Sim0 contains the state for t = 1 */
        err = gretl_matrix_copy_values(K.a0, Sim0);
        if (err) {
            return err;
        }
    }

    yt = gretl_zero_matrix_new(K.n, 1);
    if (yt == null) {
        return GretlError.alloc;
    }

    if (K.p > 0) {
        et = gretl_matrix_alloc(K.p, 1);
        if (et == null) {
            gretl_matrix_free(yt);
            return GretlError.alloc;
        }
    }

    if (Y.cols == K.r + K.n) {
        /* combined (state, obs) in @Y */
        S = Y;
        obs_offset = K.r;
    }

    if (K.p == 0 && K.VY != null) {
        /* we want to read observation disturbances */
        obsdist = 1;
    }

    for (K.t = tmin; K.t < K.N && !err; K.t += 1) {
        if (filter_is_varying(K)) {
            err = kalman_refresh_matrices(K, prn);
            if (err) {
                break;
            }
        }

        /* y_t = B'*x_t + Z'*a_t + w_t */
        gretl_matrix_multiply_mod(K.ZT, GretlMatrixMod.transpose,
                                  K.a0, GretlMatrixMod.none,
                                  yt, GretlMatrixMod.none);
        if (K.BT != null) {
            /* handle missing values? */
            kalman_do_Bx(K, yt, GretlMatrixMod.cumulate);
        }
        if (K.p > 0) {
            /* G \varepsilon_t */
            load_from_row(et, U, K.t);
            gretl_matrix_multiply(K.G, et, K.vt);
        } else if (obsdist) {
            vector_from_row_offset(K.vt, U, K.t, K.r);
        }
        gretl_matrix_add_to(yt, K.vt);

        /* record the t-dated observables */
        record_to_row_offset(Y, yt, K.t, obs_offset);

        /* record the t-dated state? */
        if (S != null && tmin == 0) {
            record_to_row_offset(S, K.a0, K.t, 0);
        }

        /* a_{t+1} = T*a_t + v_t */
        gretl_matrix_multiply(K.T, K.a0, K.a1);
        if (K.p > 0) {
            /* H \varepsilon_t */
            gretl_matrix_multiply_mod(K.H, GretlMatrixMod.none,
                                      et, GretlMatrixMod.none,
                                      K.a1, GretlMatrixMod.cumulate);
        } else {
            vector_from_row_mod(K.a1, U, K.t, GretlMatrixMod.cumulate);
        }

        if (K.mu != null) {
            gretl_matrix_add_to(K.a1, K.mu);
        }

        /* record the (t+1)-dated state? */
        if (S != null && tmin == 1) {
            record_to_row_offset(S, K.a1, K.t, 0);
        }

        fast_copy_values(K.a0, K.a1);
    }

    gretl_matrix_free(yt);
    gretl_matrix_free(et);

    return err;
}

private int check_simul_inputs(Kalman * K, const(GretlMatrixRaw) * U, const(GretlMatrixRaw) * Sim0, const(GretlMatrixRaw) * SimX, Prn * prn)
{
    int err = 0;

    if (U == null) {
        err = missing_matrix_error("U");
    } else {
        int ncols;

        if (K.p > 0) {
            /* cross-correlated */
            ncols = K.p;
        } else {
            ncols = K.VY == null ? K.r : K.r + K.n;
        }

        if (U.cols != ncols) {
            pprintf(prn, _("U should have %d columns but has %d\n").ptr,
                    ncols, U.cols);
            err = GretlError.nonconf;
        }
    }

    if (!err && Sim0 != null) {
        if (Sim0.rows != K.r || Sim0.cols != 1) {
            pprintf(prn, _("simstart should be %d x %d, is %d x %d\n").ptr,
                    K.r, 1, Sim0.rows, Sim0.cols);
        }
    }

    if (!err && K.x != null) {
        /* do we have enough "obsx" data? */
        const(GretlMatrixRaw) * X = SimX != null ? SimX : K.x;

        if (X.rows < U.rows) {
            pprintf(prn, _("obsx should have %d rows but has %d\n").ptr,
                    U.rows, X.rows);
            err = GretlError.nonconf;
        }
    }

    return err;
}

extern (C) GretlMatrixRaw * kalman_bundle_simulate(GretlBundle * b, const(GretlMatrixRaw) * U, int get_state, Prn * prn, int * err)
{
    Kalman *K = cast(Kalman *) gretl_bundle_get_private_data(b);
    const(GretlMatrixRaw) * Sim0 = null;
    const(GretlMatrixRaw) * SimX = null;
    GretlMatrixRaw * Ret = null;
    GretlMatrixRaw * savex = null;
    int saveN;

    if (K == null) {
        *err = GretlError.data;
        return null;
    }

    /* try accessing auxiliary info from the bundle */
    Sim0 = gretl_bundle_get_matrix(b, "simstart", null);
    if (K.x != null) {
        SimX = gretl_bundle_get_matrix(b, "simx", null);
    }

    *err = check_simul_inputs(K, U, Sim0, SimX, prn);
    if (*err) {
        return null;
    }

    K.b = b; /* attach bundle pointer */

    saveN = K.N;
    savex = K.x;

    /* we let U temporarily define the sample length */
    K.N = U.rows;

    /* and we allow temporary replacement of K.x */
    if (SimX != null) {
        K.x = cast(GretlMatrixRaw *) SimX;
    }

    /* set state */
    K.flags |= KFlags.sim;

    /* now, are the other needed matrices in place? */
    *err = kalman_bundle_recheck_matrices(K, prn);

    /* matrix to hold simulated observables, and state
       if wanted */
    if (!*err) {
        int ncols = get_state ? (K.r + K.n) : K.n;

        Ret = gretl_matrix_alloc(K.N, ncols);
        if (Ret == null) {
            *err = GretlError.alloc;
        }
    }

    if (!*err) {
        *err = kalman_simulate(K, U, Sim0, Ret, null, prn);
    }

    if (*err) {
        gretl_matrix_free(Ret);
        Ret = null;
    }

    /* restore state */
    K.flags &= ~KFlags.sim;
    K.N = saveN;
    K.x = savex;

    return Ret;
}

private int matrix_is_diagonal(const(GretlMatrixRaw) * m)
{
    int i;
    int j;

    for (j=0; j<m.cols; j++) {
        for (i=0; i<m.rows; i++) {
            if (i != j && gretl_matrix_get(m, i, j) != 0.0) {
                return 0;
            }
        }
    }

    return 1;
}

private int simdata_refresh_QR(Kalman * K, Prn * prn)
{
    int err = 0;

    if (matrix_is_varying(K, KMat.vs) || matrix_is_varying(K, KMat.vy)) {
        err = kalman_update_matrices(K, prn);
    }

    return err;
}

/* Return a matrix in which the standard normal variates
   in @U are scaled according to K.VS, and K.VY if present.
*/

extern (C) GretlMatrixRaw * kalman_bundle_simdata(GretlBundle * b, const(GretlMatrixRaw) * U, Prn * prn, int * err)
{
    Kalman *K = cast(Kalman *) gretl_bundle_get_private_data(b);
    GretlMatrixRaw * E = null;

    if (K == null || U == null) {
        *err = GretlError.data;
        return null;
    }

    if (K.p > 0) {
        if (U.cols != K.p) {
            *err = GretlError.data;
            return null;
        } else {
            /* nothing to be done, really */
            E = gretl_matrix_copy(U);
        }
    } else {
        int n = K.VY == null ? 0 : K.n;
        int t;
        int j;
        int rn = K.r + n;
        int T = U.rows;
        int varying = 0;
        double vjj;
        double utj;

        if (U.cols != rn) {
            *err = GretlError.data;
            return null;
        }

        E = gretl_matrix_alloc(U.rows, rn);
        if (E == null) {
            *err = GretlError.alloc;
            return null;
        }

        if (matrix_is_varying(K, KMat.vs) || matrix_is_varying(K, KMat.vy)) {
            varying = 1;
        }

        K.b = b;
        set_kalman_running(K);

        if (matrix_is_diagonal(K.VS) &&
            (K.VY == null || matrix_is_diagonal(K.VY))) {
            for (t=0; t<T && !*err; t++) {
                if (varying) {
                    K.t = t;
                    *err = simdata_refresh_QR(K, prn);
                }
                for (j=0; j<rn && !*err; j++) {
                    if (j < K.r) {
                        vjj = gretl_matrix_get(K.VS, j, j);
                    } else {
                        vjj = gretl_matrix_get(K.VY, j-K.r, j-K.r);
                    }
                    utj = gretl_matrix_get(U, t, j);
                    gretl_matrix_set(E, t, j, sqrt(vjj) * utj);
                }
            }
        } else {
            GretlMatrixRaw * V = gretl_zero_matrix_new(rn, rn);
            GretlMatrixRaw * Ut = null;
            GretlMatrixRaw * Et = null;

            if (V == null) {
                *err = GretlError.alloc;
                goto bailout;
            }

            if (varying) {
                Ut = gretl_matrix_alloc(1, rn);
                Et = gretl_matrix_alloc(1, rn);
                if (Ut == null || Et == null) {
                    gretl_matrix_free(V);
                    *err = GretlError.alloc;
                    goto bailout;
                }
            }

            if (varying) {
                for (t=0; t<T && !*err; t++) {
                    K.t = t;
                    *err = simdata_refresh_QR(K, prn);
                    if (!*err) {
                        gretl_matrix_inscribe_matrix(V, K.VS, 0, 0,
                                                     GretlMatrixMod.none);
                        if (n > 0) {
                            gretl_matrix_inscribe_matrix(V, K.VY, K.r, K.r,
                                                         GretlMatrixMod.none);
                        }
                        *err = gretl_matrix_psd_root(V, 0);
                        if (*err) {
                            gretl_errmsg_set(_("Failed to compute factor of Omega_t").ptr);
                        } else {
                            load_from_row(Ut, U, t);
                            gretl_matrix_multiply_mod(Ut, GretlMatrixMod.none,
                                                      V, GretlMatrixMod.transpose,
                                                      Et, GretlMatrixMod.none);
                            record_to_row(E, Et, t);
                        }
                    }
                }

                gretl_matrix_free(Ut);
                gretl_matrix_free(Et);
            } else {
                gretl_matrix_inscribe_matrix(V, K.VS, 0, 0,
                                             GretlMatrixMod.none);
                if (n > 0) {
                    gretl_matrix_inscribe_matrix(V, K.VY, K.r, K.r,
                                                 GretlMatrixMod.none);
                }
                *err = gretl_matrix_psd_root(V, 0);
                if (*err) {
                    gretl_errmsg_set(_("Failed to compute factor of Omega").ptr);
                } else {
                    gretl_matrix_multiply_mod(U, GretlMatrixMod.none,
                                              V, GretlMatrixMod.transpose,
                                              E, GretlMatrixMod.none);
                }
            }

            gretl_matrix_free(V);
        }
    }

 bailout:

    set_kalman_stopped(K);

    if (E == null && !*err) {
        *err = GretlError.alloc;
    } else if (E != null && *err) {
        gretl_matrix_free(E);
        E = null;
    }

    return E;
}

/*
   below: functions to support the "wrapping" of a Kalman
   struct in a gretl bundle
*/

private int check_replacement_dims(const(GretlMatrixRaw) * orig, const(GretlMatrixRaw) * repl, int id)
{
    int err = 0;

    if (id == KMat.zt || id == KMat.vy || id == KMat.t || id == KMat.vs ||
        id == KMat.bt || id == KMat.m || id == KMat.a || id == KMat.p || id == KMat.r) {
        if (repl.rows != orig.rows || repl.cols != orig.cols) {
            err = GretlError.data;
        }
    } else if (id == KMat.y || id == KMat.x) {
        if (repl.cols != orig.cols) {
            err = GretlError.data;
        }
    }

    if (err) {
        gretl_errmsg_set(_("You cannot resize a state-space system matrix").ptr);
    }

    return err;
}

/* On input @targ is a pointer to the matrix to be replaced, if
   there's already a matrix corresponding to @id in place; @src is
   the new matrix, to be copied in if @copy is non-zero or
   otherwise just attached. On output @targ points to the new
   matrix.
*/

private int add_or_replace_k_matrix(Kalman * K, GretlMatrixRaw ** targ, GretlMatrixRaw * src, int id, int copy)
{
    int err = 0;

    if (*targ != src) {
        if (*targ != null) {
            err = check_replacement_dims(*targ, src, id);
            if (err) {
                fprintf(stderr, " bad kalman matrix replacement!\n");
                return err;
            }
            fast_copy_values(*targ, src);
            if (!copy) {
                gretl_matrix_free(src);
            }
        } else if (copy) {
            *targ = gretl_matrix_copy(src);
            if (*targ == null) {
                err = GretlError.alloc;
            }
        } else {
            *targ = src;
        }
    }

    return err;
}

private GretlMatrixRaw ** get_input_matrix_target_by_id(Kalman * K, int i)
{
    GretlMatrixRaw ** targ = null;

    if (i == KMat.t) {
        targ = &K.T;
    } else if (i == KMat.bt) {
        targ = &K.BT;
    } else if (i == KMat.zt) {
        targ = &K.ZT;
    } else if (i == KMat.vs) {
        /* variance of state */
        if (kalman_djvar(K)) {
            targ = &K.H;
        } else if (kalman_dkvar(K)) {
            targ = &K.Q;
        } else {
            targ = &K.VS;
        }
    } else if (i == KMat.vy) {
        /* variance of observable */
        if (kalman_djvar(K)) {
            targ = &K.G;
        } else {
            targ = &K.VY;
        }
    } else if (i == KMat.m) {
        targ = &K.mu;
    } else if (i == KMat.y) {
        targ = &K.y;
    } else if (i == KMat.x) {
        targ = &K.x;
    } else if (i == KMat.a) {
        targ = &K.aini;
    } else if (i == KMat.p) {
        targ = &K.Pini;
    } else if (i == KMat.r) {
        targ = &K.R;
    }

    return targ;
}

/* Try attaching a matrix to a Kalman bundle: similar to
   kalman_bundle_set_matrix() below, but allowing for the possibility
   that the @data input of type @vtype has to be converted first.
*/

private int kalman_bundle_try_set_matrix(Kalman *K, void *data, GretlType vtype, int id, int copy)
{
    GretlMatrixRaw ** targ;
    int err = 0;

    /* determine the location for this matrix */
    targ = get_input_matrix_target_by_id(K, id);

    if (targ == null) {
        err = GretlError.data;
    } else {
        GretlMatrixRaw * m;

        if (vtype == GretlType.matrix) {
            m = cast(GretlMatrixRaw *) data;
            err = add_or_replace_k_matrix(K, targ, m, id, copy);
        } else if (vtype == GretlType.double_) {
            m = gretl_matrix_from_scalar(*cast(double *) data);
            if (m == null) {
                err = GretlError.alloc;
            } else {
                err = add_or_replace_k_matrix(K, targ, m, id, 0);
            }
        } else {
            err = GretlError.types;
        }
    }

    return err;
}

/* Called by kalman_deserialize() when reconstructing a Kalman bundle
   from XML, and also by kalman_bundle_copy() when duplicating such a
   bundle. In these contexts we know we have a gretl_matrix on input.
*/

private int kalman_bundle_set_matrix(Kalman *K, GretlMatrixRaw * m, int i, int copy)
{
    GretlMatrixRaw ** targ;
    int err = 0;

    targ = get_input_matrix_target_by_id(K, i);

    if (targ == null) {
        err = GretlError.data;
    } else {
        err = add_or_replace_k_matrix(K, targ, m, i, copy);
    }

    return err;
}

private immutable(char) * [9] kalman_output_matrix_names = [
    "prederr",
    "pevar",
    "state",
    "stvar",
    "gain",
    "llt",
    "smdist",
    "smdisterr",
    "uhat"
];

enum K_N_OUTPUTS = 9;

private GretlMatrixRaw ** kalman_output_matrix(Kalman * K, const(char) * key)
{
    GretlMatrixRaw ** pm = null;

    if (!strcmp(key, "prederr")) {
        pm = &K.V;
    } else if (!strcmp(key, "pevar")) {
        pm = &K.F;
    } else if (!strcmp(key, "state")) {
        pm = &K.A;
    } else if (!strcmp(key, "stvar")) {
        pm = &K.P;
    } else if (!strcmp(key, "gain")) {
        pm = &K.K;
    } else if (!strcmp(key, "llt")) {
        pm = &K.LL;
    } else if (!strcmp(key, "smdist")) {
        pm = &K.U;
    } else if (!strcmp(key, "smdisterr")) {
        pm = &K.Vsd;
    } else if (!strcmp(key, "uhat")) {
        pm = &K.vt;
    }

    return pm;
}

private int is_on_demand(const(char) * key)
{
    return (strcmp(key, "smdist") == 0 ||
            strcmp(key, "smdisterr") == 0);
}

private int output_matrix_slot(const(char) * s)
{
    int i;

    for (i=0; i<K_N_OUTPUTS; i++) {
        if (!strcmp(s, kalman_output_matrix_names[i])) {
            return i;
        }
    }

    return -1;
}

enum K_N_SCALARS = 15;

enum KScalar {
    t = 0,
    diffuse,
    exact,
    cross,
    dkvar,
    seq,
    dej,
    s2,
    lnl,
    r,
    n,
    N,
    p,
    d,
    j
}

private immutable(char) * [K_N_SCALARS] kalman_output_scalar_names = [
    "t",
    "diffuse",
    "exact",
    "cross",
    "dkvar",
    "sequential",
    "dejong",
    "s2",
    "lnl",
    "r",
    "n",
    "N",
    "p",
    "d",
    "j"
];

private double * kalman_output_scalar(Kalman * K, const(char) * key)
{
    /* static storage for on-the-fly scalars */
    static __gshared double[K_N_SCALARS] retval = 0;
    int i;
    int idx = -1;

    for (i=0; i<K_N_SCALARS; i++) {
        if (!strcmp(key, kalman_output_scalar_names[i])) {
            idx = i;
            break;
        }
    }

    if (idx < 0 && !strcmp(key, "T")) {
        /* backward compatibility */
        idx = KScalar.N;
    }

    if (idx < 0) {
        return null;
    }

    switch (idx) {
    case KScalar.t:
        if (kalman_is_running(K)) {
            retval[idx] = K.t + 1;
        } else {
            retval[idx] = kalman_checking(K) ? 1 : 0;
        }
        break;
    case KScalar.diffuse:
        retval[idx] = kalman_diffuse(K) + K.exact;
        break;
    case KScalar.exact:
        retval[idx] = K.exact;
        break;
    case KScalar.cross:
        retval[idx] = (K.vartype == KVarType.dj_var);
        break;
    case KScalar.dkvar:
        retval[idx] = (K.vartype == KVarType.dk_var);
        break;
    case KScalar.seq:
        retval[idx] = kalman_sequential(K);
        break;
    case KScalar.dej:
        retval[idx] = kalman_dejong(K);
        break;
    case KScalar.s2:
        retval[idx] = K.s2;
        break;
    case KScalar.lnl:
        retval[idx] = K.loglik;
        break;
    case KScalar.r:
        retval[idx] = K.r;
        break;
    case KScalar.n:
        retval[idx] = K.n;
        break;
    case KScalar.N:
        retval[idx] = K.N;
        break;
    case KScalar.p:
        retval[idx] = K.p;
        break;
    case KScalar.d:
        retval[idx] = K.d;
        break;
    case KScalar.j:
        retval[idx] = K.j;
        break;
    default:
        break;
    }

    return &retval[idx];
}

private const(GretlMatrixRaw) * k_input_matrix_by_id(Kalman * K, int i)
{
    const(GretlMatrixRaw) * m = null;

    if (i == KMat.t) {
        m = K.T;
    } else if (i == KMat.bt) {
        m = K.BT;
    } else if (i == KMat.zt) {
        m = K.ZT;
    } else if (i == KMat.vs) {
        if (kalman_djvar(K)) {
            m = K.H;
        } else if (kalman_dkvar(K)) {
            m = K.Q;
        } else {
            m = K.VS;
        }
    } else if (i == KMat.vy) {
        if (kalman_djvar(K)) {
            m = K.G;
        } else {
            m = K.VY;
        }
    } else if (i == KMat.m) {
        m = K.mu;
    } else if (i == KMat.y) {
        m = K.y;
    } else if (i == KMat.x) {
        m = K.x;
    } else if (i == KMat.a) {
        m = K.aini;
    } else if (i == KMat.p) {
        m = K.Pini;
    } else if (i == KMat.r) {
        m = K.R;
    }

    return m;
}

private int input_matrix_slot(const(char) * s)
{
    int i;

    for (i=0; i<KMat.mmax; i++) {
        if (!strcmp(s, K_input_mats[i].name)) {
            return K_input_mats[i].sym;
        }
    }

    return -1;
};

private GretlType kalman_extra_type(const(char) * key)
{
    if (!strcmp(key, "simstart") || !strcmp(key, "simx")) {
        return GretlType.matrix;
    } else {
        return GretlType.none;
    }
}

/* Respond to "diffuse" setting: 0 = no, 1 = traditional, 2 = exact */

private int kalman_set_diffuse(Kalman * K, int d)
{
    int err = 0;

    if (trace) {
        printf("kalman_set_diffuse: d = %d\n", d);
    }
    if (d != 0 && d != 1 && d != 2) {
        return GretlError.invarg;
    } else if (d > 0) {
	K.flags |= KFlags.diffuse;
	if (d == 2) {
	    /* exact initial diffuse mode */
	    K.exact = 1;
	    if (K.code == KCodeType.legacy) {
		/* the legacy code doesn't handle "exact" */
		K.code = kalman_djvar(K) ? KCodeType.dejong : KCodeType.sequen;
	    }
	} else {
	    /* d == 1: traditional "kappa" mode */
	    K.exact = 0;
	}
        err = diffuse_Pini(K);
    } else {
        /* d == 0 */
        K.exact = 0;
        K.flags &= ~KFlags.diffuse;
        if (K.Pk0 != null) {
            gretl_matrix_free(K.Pk0);
            K.Pk0 = null;
        }
    }

    return err;
}

/* Respond to user's setting of "sequential" or "dejong": @code tells
   us which term we got, and @s tells us whether the value was 0 or 1.
*/

private int kalman_set_code(Kalman * K, int code, int s)
{
    if (s && code == KCodeType.sequen && K.vartype == KVarType.dj_var) {
        gretl_errmsg_set(_("kalman: the 'sequential' setting is not compatible with\n" ~
                         "cross-correlated disturbances").ptr);
        return GretlError.invarg;
    } else if (s && (code != K.user_code)) {
        /* select the specified optional code */
        K.user_code = code;
    } else if (!s && (code == K.user_code)) {
        /* deselect the specified optional code */
        K.user_code = KCodeType.legacy;
    }

    /* is this clause actually in the right place? */
    if (K.user_code == KCodeType.sequen && kalman_diffuse(K) && K.Pk0 == null) {
        K.Pk0 = gretl_identity_matrix_new(K.r);
        if (K.Pk0 == null) {
            return GretlError.alloc;
        }
    }

    if (trace) {
        printf("kalman_set_code: code %d and s %d . K.code %d\n",
               code, s, K.code);
    }

    return 0;
}

/* Called by real_bundle_set_data() in gretl_bundle.c.  The return
   value indicates whether the putative setting was handled (1) or not
   (0). Not being handled here is not necessarily an error.

   The @copy flag here is inherited from the specific caller of
   real_bundle_set_data(): @copy = 1 if the caller was
   gretl_bundle_set_data(), 0 if it was gretl_bundle_donate_data().
   Either way the Kalman struct takes ownership.
*/

extern (C) int maybe_set_kalman_element(void * kptr, const(char) * key, void * vptr, GretlType vtype, int copy, int * err)
{
    GretlType targtype;
    Kalman * K = cast(Kalman *) kptr;
    int fncall = 0;
    int i;
    int id = -1;
    int done = 0;

    if (K == null) {
        *err = GretlError.data;
        return 0;
    }

    if (kdebug) {
	fprintf(stderr, "maybe_set_kalman_element: key '%s', type %s\n",
		key, gretl_type_get_name(vtype));
    }

    /* Check for optional "extra" Kalman matrices that
       live outside of the Kalman struct itself.
    */
    targtype = kalman_extra_type(key);
    if (targtype != GretlType.none) {
        if (vtype != targtype) {
            *err = GretlError.types;
        }
        return 0;
    }

    if (!strcmp(key, "diffuse") || !strcmp(key, "sequential") ||
        !strcmp(key, "dejong") || !strcmp(key, "extra")) {
        /* scalar config settings */
        if (vtype == GretlType.double_) {
            double v = *cast(double *) vptr;

            if (!strcmp(key, "diffuse")) {
                *err = kalman_set_diffuse(K, cast(int) v);
            } else if (!strcmp(key, "sequential")) {
                *err = kalman_set_code(K, KCodeType.sequen, cast(int) v);
	    } else if (!strcmp(key, "dejong")) {
		*err = kalman_set_code(K, KCodeType.dejong, cast(int) v);
            } else {
                /* note: just for KFAS comparison */
                K.flags |= KFlags.extra;
            }
        } else {
            *err = GretlError.types;
        }
        return 1; /* done, error or no */
    }

    if (!strcmp(key, "timevar_call")) {
        /* try for a function call specifier (string) */
        if (vtype == GretlType.string_) {
            fncall = 1;
        } else {
            *err = GretlError.types;
        }
    } else {
        /* try for a matrix specifier */
        i = input_matrix_slot(key);
        if (i >= 0) {
            if (vtype == GretlType.matrix ||
                vtype == GretlType.double_) {
                id = i;
            } else {
		gretl_errmsg_sprintf(_("%s: expected matrix but got %s").ptr,
				     key, gretl_type_get_name(vtype));
                *err = GretlError.types;
            }
        }
    }

    if (*err) {
        return 0;
    } else if (fncall) {
        if (copy) {
            K.matcall = gretl_strdup(cast(char *) vptr);
        } else {
            K.matcall = cast(char *) vptr;
        }
        /* re-evaluate what's actually varying */
        *err = check_for_matrix_updates(K, null);
        if (!*err) {
            done = 1;
        }
    } else if (id < 0) {
        if (kalman_output_matrix(K, key) != null ||
            kalman_output_scalar(K, key) != null) {
            *err = GretlError.data;
            gretl_errmsg_sprintf(_("The member %s is read-only").ptr, key);
        }
    } else {
        *err = kalman_bundle_try_set_matrix(K, vptr, vtype, id, copy);
        done = (*err == 0);
    }

    return done;
}

extern (C) int maybe_delete_kalman_element(void * kptr, const(char) * key, int * err)
{
    Kalman * K = cast(Kalman *) kptr;
    GretlMatrixRaw ** pm;
    int done = 0;

    if (K == null) {
        return 0;
    }

    if (kalman_output_scalar(K, key) != null ||
        input_matrix_slot(key) >= 0 || !strcmp(key, "uhat")) {
        /* note: the matrix under the key "uhat" is part of
           the internal Kalman apparatus
	*/
        gretl_errmsg_sprintf(_("%s: cannot be deleted").ptr, key);
        *err = GretlError.data;
    } else if ((pm = kalman_output_matrix(K, key)) != null) {
        /* OK to delete a user-output matrix? */
        gretl_matrix_free(*pm);
        *pm = null;
    } else if (!strcmp(key, "timevar_call")) {
        /* OK to delete time-variation call */
        if (K.matcall != null) {
            free(K.matcall);
            K.matcall = null;
            free(K.varying);
            K.varying = null;
            done = 1;
        } else {
            *err = GretlError.invarg;
        }
    }

    return done;
}

/* Take the smoothed disturbances recorded in K.b.etahat and
   K.b.epshat (if present) and write them into K.U (aka "smdist") in
   the format documented in the Guide.
*/

private GretlMatrixRaw * fill_smdist(Kalman * K)
{
    GretlMatrixRaw * U1 = gretl_bundle_get_matrix(K.b, "etahat", null);
    GretlMatrixRaw * U2 = gretl_bundle_get_matrix(K.b, "epshat", null);
    GretlMatrixRaw * U = K.U;

    if (U1 != null) {
        int j;
        int t;
        int q = U1.cols;
        int n = (U2 != null)? U2.cols : 0;
        double utj;

        if (U == null || U.rows != K.N || U.cols != q + n) {
            /* we need new storage */
            gretl_matrix_free(U);
            K.U = U = gretl_matrix_alloc(K.N, q + n);
            if (U == null) {
                return null;
            }
        }
        for (t=0; t<K.N; t++) {
            for (j=0; j<q; j++) {
                utj = gretl_matrix_get(U1, t, j);
                gretl_matrix_set(U, t, j, utj);
            }
            for (j=0; j<n; j++) {
                utj = gretl_matrix_get(U2, t, j);
                gretl_matrix_set(U, t, j+q, utj);
            }
        }
    }

    return U;
}

/* Take the representations of the variance of the smoothed
   disturbances recorded in K.b.veta and K.b.veps (if present) and
   write their square roots into K.Vsd (aka "smdisterr") in the
   format documented in the Guide.
*/

private GretlMatrixRaw * fill_smdisterr(Kalman * K)
{
    GretlMatrixRaw * V1 = gretl_bundle_get_matrix(K.b, "veta", null);
    GretlMatrixRaw * V2 = gretl_bundle_get_matrix(K.b, "veps", null);
    GretlMatrixRaw * V = K.Vsd;

    if (V1 != null) {
        int q = kalman_dkvar(K) ? K.q : K.r;
        int n = (V2 != null)? K.n : 0;
        int c2 = (V2 != null)? V2.cols : 0;
        int c1 = V1.cols;
        int diag1 = 0;
        int diag2 = 0;
        int j;
        int k;
        int t;
        double vtj;

        diag1 = (q > 1 && c1 == q*q);
        diag2 = (n > 1 && c2 == n*n);

        if (V == null || V.rows != K.N || V.cols != q + n) {
            /* we need new storage */
            gretl_matrix_free(V);
            K.Vsd = V = gretl_zero_matrix_new(K.N, q + n);
            if (V == null) {
                return null;
            }
        }
        k = 0;
        for (j=0; j<q; j++) {
            for (t=0; t<K.N; t++) {
                vtj = gretl_matrix_get(V1, t, diag1 ? j*(q+1) : j);
                gretl_matrix_set(V, t, k, sqrt(vtj));
            }
            k++;
        }
        for (j=0; j<n; j++) {
            for (t=0; t<K.N; t++) {
                vtj = gretl_matrix_get(V2, t, diag2 ? j*(n+1) : j);
                gretl_matrix_set(V, t, k, sqrt(vtj));
            }
            k++;
        }
    }

    return V;
}

/* When the sequential method is being used but the observable is
   multivariate, transform K.V and K.F for compatibility with the
   standard multivariate output. This corresponds to the function
   mvInnovations in KFAS. Added March, 2025.
*/

private int mv_transform(Kalman * K)
{
    GretlMatrixRaw * Vt = gretl_zero_matrix_new(K.n, 1);
    GretlMatrixRaw * PM = gretl_zero_matrix_new(K.n, K.r);
    GretlMatrixRaw * At = gretl_zero_matrix_new(K.r, 1);
    GretlMatrixRaw * Pt = gretl_zero_matrix_new(K.r, K.r);
    GretlMatrixRaw * Z = K.sinfo.Z;
    GretlMatrixRaw * Ft = K.Ft;
    int t;

    if (Vt == null || PM == null) {
        return GretlError.alloc;
    }

    for (t=0; t<K.N; t++) {
        if (matrix_is_varying(K, KMat.zt)) {
            K.t = t;
            retrieve_Zt(K);
        }
        load_from_row(At, K.A, t);
        load_from_row(Vt, K.y, t);
        gretl_matrix_multiply_mod(Z,  GretlMatrixMod.none,
                                  At, GretlMatrixMod.none,
                                  Vt, GretlMatrixMod.decrement);
        record_to_row(K.V, Vt, t);
        load_from_vech(Pt, K.P, K.r, t);
        gretl_matrix_multiply(Z, Pt, PM);
        gretl_matrix_copy_values(Ft, K.VY);
        gretl_matrix_multiply_mod(PM, GretlMatrixMod.none,
                                  Z,  GretlMatrixMod.transpose,
                                  Ft, GretlMatrixMod.cumulate);
        record_to_vech(K.F, Ft, K.n, t);
    }

    gretl_matrix_free(Vt);
    gretl_matrix_free(PM);
    gretl_matrix_free(At);
    gretl_matrix_free(Pt);

    return 0;
}

private GretlMatrixRaw * construct_kalman_matrix(Kalman * K, const(char) * key)
{
    GretlMatrixRaw * m = null;

    if (!strcmp(key, "smdist")) {
        m = fill_smdist(K);
    } else if (!strcmp(key, "smdisterr")) {
        m = fill_smdisterr(K);
    }

    return m;
}

/* maybe_retrieve_kalman_element():

   This function handles retrieval of elements from a Kalman bundle by
   @key. Such elements may be of type matrix, scalar or string.

   Most retrievable items will already be available (or not available)
   for access, but a couple are 'constructables' that will be built on
   demand.

   We set *private_ to 1 if the given @key refers to a member of the
   Kalman struct inside the bundle, as opposed to a regular bundle
   member.
*/

extern (C) void * maybe_retrieve_kalman_element(void * kptr, const(char) * key, GretlType * type, int * private_, int * err)
{
    GretlMatrixRaw ** pm = null;
    Kalman * K = cast(Kalman *) kptr;
    void * ret = null;
    int i;
    int id = -1;

    if (K == null) {
        *err = GretlError.data;
        return null;
    }

    if (kdebug) {
	fprintf(stderr, "maybe_retrieve_kalman_element: key '%s', K.b %p\n",
                key, cast(void *) K.b);
    }

    /* initialize to nothing doing */
    *type = GretlType.none;

    if (!strcmp(key, "timevar_call")) {
        /* try for time-variation function call */
        *private_ = 1;
        if (K.matcall != null) {
            ret = K.matcall;
            *type = GretlType.string_;
        }
        return ret; /* handled */
    }

    /* try for an Kalman input matrix specifier */
    for (i=0; i<KMat.mmax; i++) {
        if (!strcmp(key, K_input_mats[i].name)) {
            *private_ = 1;
            id = K_input_mats[i].sym;
            ret = cast(GretlMatrixRaw *) k_input_matrix_by_id(K, id);
            if (ret != null) {
                *type = GretlType.matrix;
            }
            return ret; /* handled */
        }
    }

    /* try for an output matrix */
    pm = kalman_output_matrix(K, key);
    if (pm != null) {
        *private_ = 1;
        if (is_on_demand(key)) {
            ret = construct_kalman_matrix(K, key);
        } else if (*pm != null) {
            ret = *pm;
        }
        if (ret != null) {
            *type = GretlType.matrix;
        }
        return ret; /* handled */
    }

    /* no match yet: try scalar member */
    ret = kalman_output_scalar(K, key);
    if (ret != null) {
        *private_ = 1;
        *type = GretlType.double_;
    }

    return ret;
}

private int output_matrix_count(Kalman * K)
{
    GretlMatrixRaw ** pm;
    int i;
    int n = 0;

    for (i=0; i<K_N_OUTPUTS; i++) {
        pm = kalman_output_matrix(K, kalman_output_matrix_names[i]);
        n += (pm != null && *pm != null);
    }

    return n;
}

extern (C) int print_kalman_bundle_info(void * kptr, Prn * prn)
{
    Kalman * K = cast(Kalman *) kptr;
    int err = 0;

    if (K == null) {
        pputs(prn, _("Kalman struct: empty\n").ptr);
        err = GretlError.data;
    } else {
        const(GretlMatrixRaw) * m;
        GretlMatrixRaw ** pm;
        double * px;
        const(char) * name;
        int i;
        int id;

        pputs(prn, _("\nKalman input matrices\n").ptr);

        for (i=0; i<KMat.mmax; i++) {
            id = K_input_mats[i].sym;
            m = k_input_matrix_by_id(K, id);
            if (m != null) {
                pprintf(prn, " %s: ", K_input_mats[i].name);
                pprintf(prn, "%d x %d\n", m.rows, m.cols);
            }
        }

        if (output_matrix_count(K) > 0) {
            pputs(prn, _("\nKalman output matrices\n").ptr);
            for (i=0; i<K_N_OUTPUTS; i++) {
                name = kalman_output_matrix_names[i];
                pm = kalman_output_matrix(K, name);
                if (pm != null && *pm != null) {
                    m = *pm;
                    pprintf(prn, " %s: ", name);
                    pprintf(prn, "%d x %d\n", m.rows, m.cols);
                }
            }
        }

        pputs(prn, _("\nKalman scalars\n").ptr);

        for (i=0; i<K_N_SCALARS; i++) {
            name = kalman_output_scalar_names[i];
            pprintf(prn, " %s: ", name);
            px = kalman_output_scalar(K, name);
            if (px == null || na(*px)) {
                pputs(prn, "NA\n");
            } else {
                pprintf(prn, "%g\n", *px);
            }
        }

        if (K.matcall != null) {
            pputs(prn, "\nKalman strings\n");
            pprintf(prn, " timevar_call: %s\n", K.matcall);
        }
    }

    return err;
}

/* For use in context of a Kalman bundle: serialize the information in
   the Kalman struct to XML
*/

extern (C) int kalman_serialize(void * kptr, Prn * prn)
{
    Kalman * K = cast(Kalman *) kptr;
    const(GretlMatrixRaw) * m;
    GretlMatrixRaw ** pm;
    double * px;
    const(char) * name;
    int i;
    int err = 0;

    if (K == null) {
        fputs("kalman_serialize: got NULL\n", stderr);
        return GretlError.data;
    }

    pputs(prn, "<gretl-kalman>\n");

    for (i=0; i<KMat.mmax; i++) {
        m = k_input_matrix_by_id(K, K_input_mats[i].sym);
        if (m != null) {
            gretl_matrix_serialize(m, K_input_mats[i].name, prn);
        }
    }

    for (i=0; i<K_N_OUTPUTS; i++) {
        name = kalman_output_matrix_names[i];
        pm = kalman_output_matrix(K, name);
        if (pm != null && *pm != null) {
            gretl_matrix_serialize(*pm, name, prn);
        }
    }

    for (i=0; i<=KScalar.lnl; i++) {
        name = kalman_output_scalar_names[i];
        px = kalman_output_scalar(K, name);
        if (px != null && !na(*px)) {
            gretl_finite_scalar_serialize(*px, name, prn);
        }
    }

    if (K.matcall != null) {
        gretl_string_serialize(K.matcall, "timevar_call", prn);
    }

    pputs(prn, "</gretl-kalman>\n");

    return err;
}

private int required_matrix_slot(const(char) * s)
{
    if (!strcmp(s, "obsy"))     return 0;
    if (!strcmp(s, "obsymat"))  return 1;
    if (!strcmp(s, "statemat")) return 2;
    if (!strcmp(s, "statevar")) return 3;
    if (!strcmp(s, "obsvar"))   return 4;
    return -1;
};

/* For use in context of a Kalman bundle: deserialize the Kalman
   struct from XML
*/

extern (C) GretlBundle * kalman_deserialize(void * p1, void * p2, int * err)
{
    xmlNodePtr cur;
    xmlNodePtr node = cast(xmlNodePtr) p1;
    xmlDocPtr doc = cast(xmlDocPtr) p2;
    GretlMatrixRaw *[5] Mreq = null;
    GretlMatrixRaw *[KMat.mmax] Mopt = null;
    GretlMatrixRaw *[K_N_OUTPUTS] Mout = null;
    char * tvcall = null;
    double s2 = NADBL;
    double lnl = NADBL;
    int[5] copy = 0;
    int i;
    int nmats = 0;
    int Kflags = 0;
    int vtype = 0;
    GretlMatrixRaw * m;
    double x;
    char * key;
    char * strv;
    GretlBundle * b = null;

    while (node != null && !*err) {
        if (!xmlStrcmp(node.name, cast(XUC) "gretl-kalman")) {
            cur = node.children;
            while (cur != null && !*err) {
                key = cast(char *) xmlGetProp(cur, cast(XUC) "name");
                if (!xmlStrcmp(cur.name, cast(XUC) "gretl-matrix")) {
                    /* pick up Kalman matrices */
                    m = gretl_xml_get_matrix(cur, doc, err);
                    if ((i = required_matrix_slot(key)) >= 0) {
                        nmats++;
                        Mreq[i] = m;
                    } else if ((i = input_matrix_slot(key)) >= 0) {
                        Mopt[i] = m;
                    } else if ((i = output_matrix_slot(key)) >= 0) {
                        Mout[i] = m;
                    }
                } else if (!xmlStrcmp(cur.name, cast(XUC) "scalar")) {
                    /* pick up Kalman scalars */
                    if (gretl_xml_get_prop_as_double(cur, "value", &x)) {
                        if (!strcmp(key, "diffuse") && x > 0) {
                            Kflags |= KFlags.diffuse;
                        } else if (!strcmp(key, "cross") && x > 0) {
                            vtype = KVarType.dj_var;
                        } else if (!strcmp(key, "dkvar") && x > 0) {
                            vtype = KVarType.dk_var;
                        } else if (!strcmp(key, "s2")) {
                            s2 = x;
                        } else if (!strcmp(key, "lnl")) {
                            lnl = x;
                        }
                    }
                } else if (!xmlStrcmp(cur.name, cast(XUC) "string")) {
                    /* pick up Kalman strings */
                    if (!strcmp(key, "timevar_call") &&
                        gretl_xml_get_prop_as_string(cur, "value", &strv)) {
                            tvcall = strv;
                    }
                }
                free(key);
                cur = cur.next;
            }
            break;
        }
        node = node.next;
    }

    if (vtype == KVarType.dk_var) {
        if (Mreq[KMat.vs] == null || Mopt[KMat.r] == null) {
            *err = GretlError.data;
            goto bailout;
        } else {
            Mopt[KMat.vy] = Mreq[4];
            Mreq[4] = Mopt[KMat.r];
            Mopt[KMat.r] = null;
        }
    } else if (vtype == KVarType.std_var && nmats == 5) {
        /* drop obsvar from initialization */
        Mopt[KMat.vy] = Mreq[4];
        Mreq[4] = null;
        nmats--;
    }

    if ((vtype > 0 && nmats != 5) || (vtype == 0 && nmats != 4)) {
        *err = GretlError.data;
    } else {
        b = kalman_bundle_new(Mreq.ptr, copy.ptr, nmats, (vtype == KVarType.dk_var)? 1 : 0, err);
        if (b != null) {
            Kalman *K = cast(Kalman *) gretl_bundle_get_private_data(b);
            GretlMatrixRaw ** pm;
            const(char) * name;

            K.flags = Kflags;
            K.vartype = cast(KVarType) vtype;
            K.s2 = s2;
            K.loglik = lnl;

            for (i=0; i<KMat.mmax; i++) {
                if (Mopt[i] != null) {
                    kalman_bundle_set_matrix(K, Mopt[i], i, 0);
                }
            }
            for (i=0; i<K_N_OUTPUTS; i++) {
                if (Mout[i] != null) {
                    name = kalman_output_matrix_names[i];
                    pm = kalman_output_matrix(K, name);
                    *pm = Mout[i];
                }
            }
            K.matcall = tvcall;
        }
    }

 bailout:

    if (*err) {
        /* clean up */
        for (i=0; i<5; i++) {
            gretl_matrix_free(Mreq[i]);
        }
        for (i=0; i<KMat.mmax; i++) {
            gretl_matrix_free(Mopt[i]);
        }
        for (i=0; i<K_N_OUTPUTS; i++) {
            gretl_matrix_free(Mout[i]);
        }
        free(tvcall);
    }

    return b;
}

/* Called from gretl_bundle.c to meet the case where the user calls
   for a Kalman bundle to be copied: here we create a new Kalman
   struct and copy across the required elements (since they are
   not regular bundle members).
*/

extern (C) GretlBundle * kalman_bundle_copy(const(GretlBundle) * src, int * err)
{
    Kalman * K;
    Kalman * Knew;
    GretlBundle * b = null;
    GretlMatrixRaw *[5] M = null;
    int[5] copy = [1, 1, 1, 1, 1];
    GretlMatrixRaw * m;
    GretlMatrixRaw ** pm;
    GretlMatrixRaw ** pm1;
    const(char) * name;
    int dkopt = 0;
    int i;
    int id;
    int k = 4;

    K = cast(Kalman *) gretl_bundle_get_private_data(cast(GretlBundle *) src);

    if (K == null) {
        *err = GretlError.data;
        return null;
    }

    /* set pointers to required Kalman matrices */
    M[0] = K.y;
    M[1] = K.ZT;
    M[2] = K.T;

    /* variants dependent on presence/absence of cross-correlation */
    if (kalman_djvar(K)) {
        M[3] = K.H;
        M[4] = K.G;
        k = 5;
    } else if (kalman_dkvar(K)) {
        M[3] = K.Q;
        M[4] = K.R;
        k = 5;
        dkopt = 1;
    } else {
        M[3] = K.VS;
    }

    b = kalman_bundle_new(M.ptr, copy.ptr, k, dkopt, err);

    if (*err) {
        return b;
    }

    Knew = cast(Kalman *) gretl_bundle_get_private_data(b);
    Knew.flags = K.flags;

    /* add any "extra" matrices, beyond the required ones */
    for (i=0; i<n_extra_mats && !*err; i++) {
        id = extra_input_mats[i];
        m = cast(GretlMatrixRaw *) k_input_matrix_by_id(K, id);
        if (m != null) {
            *err = kalman_bundle_set_matrix(Knew, m, id, 1);
        }
    }

    for (i=0; i<K_N_OUTPUTS && !*err; i++) {
        name = kalman_output_matrix_names[i];
        pm = kalman_output_matrix(K, name);
        if (pm != null && *pm != null) {
            pm1 = kalman_output_matrix(Knew, name);
            *pm1 = gretl_matrix_copy(*pm);
        }
    }

    Knew.ifc = K.ifc;
    Knew.s2 = K.s2;
    Knew.loglik = K.loglik;
    Knew.vartype = K.vartype;

    if (K.flags & KFlags.diffuse) {
        Knew.flags |= KFlags.diffuse;
    }

    if (K.matcall != null) {
        Knew.matcall = gretl_strdup(K.matcall);
    }

    return b;
}

/* for use in constructing GUI bundle save menu */

char ** kalman_bundle_get_matrix_names(Kalman * K, int * ns)
{
    char **S = null;
    GretlMatrixRaw ** pm;
    const(char) * name;
    int i;
    int id;
    int err = 0;

    *ns = 0;

    for (i=0; i<KMat.mmax && !err; i++) {
        id = K_input_mats[i].sym;
        if (k_input_matrix_by_id(K, id) != null) {
            err = strings_array_add(&S, ns, K_input_mats[i].name);
        }
    }

    for (i=0; i<K_N_OUTPUTS && !err; i++) {
        name = kalman_output_matrix_names[i];
        pm = kalman_output_matrix(K, name);
        if (pm != null && *pm != null) {
            err = strings_array_add(&S, ns, name);
        }
    }

    return S;
}

/* also for use in constructing GUI bundle save menu */

char ** kalman_bundle_get_scalar_names(Kalman * K, int * ns)
{
    char **S;

    *ns = K_N_SCALARS - na(K.s2) - na(K.loglik);
    S = strings_array_new(*ns);

    if (S != null) {
        int i;
        int j = 0;

        for (i=0; i<K_N_SCALARS; i++) {
            if (i == KScalar.s2 && na(K.s2)) {
                continue;
            } else if (i == KScalar.lnl && na(K.loglik)) {
                continue;
            } else {
                S[j++] = gretl_strdup(kalman_output_scalar_names[i]);
            }
        }
    }

    return S;
}

/* to support the nelem() function for Kalman bundles */

extern (C) int kalman_bundle_n_members(GretlBundle * b)
{
    Kalman *K = cast(Kalman *) gretl_bundle_get_private_data(b);
    int n = 0;

    if (K != null) {
        int i;
        int id;

        n = K_N_SCALARS;

        for (i=0; i<KMat.mmax; i++) {
            id = K_input_mats[i].sym;
            n += (k_input_matrix_by_id(K, id) != null);
        }

        n += output_matrix_count(K);
    }

    return n;
}

/* Start of sequential filtering and smoothing code. Based on Durbin
   and Koopman, Time Series Analysis by State Space Methods, 2e
   (Oxford, 2012), section 6.4 in particular. The C implementation
   here owes a lot to Jouni Helske's fortran code in his KFAS package
   for R (see kfilter2.f90). Some helper functions come first.
*/

private void dj_from_Finf(const(GretlMatrixRaw) * Finf, int * d, int * j, double tiny)
{
    int cmax = Finf.cols - 1;
    int rmax = Finf.rows - 1;
    int i;
    int k;

    *d = *j = 0;

    for (k=cmax; k >= 0 && *d == 0; k--) {
        for (i=rmax; i>=0; i--) {
            if (gretl_matrix_get(Finf, i, k) > tiny) {
                *d = k + 1;
                *j = i + 1;
                break;
            }
        }
    }
}

private int bundle_add_matrix(GretlBundle * b, const(char) * key, GretlMatrixRaw * m)
{
    return gretl_bundle_donate_data(b, key, m, GretlType.matrix);
}

/* simple implementation assuming no error check required */

private double dotprod(const(GretlMatrixRaw) * a, const(GretlMatrixRaw) * b)
{
    int i, n = gretl_vector_get_length(a);
    double dp = 0;

    for (i=0; i<n; i++) {
        dp += a.val[i] * b.val[i];
    }

    return dp;
}

private void increment_state(GretlMatrixRaw * a, const(GretlMatrixRaw) * b, GretlMatrixRaw * tmp, double x)
{
    fast_copy_values(tmp, b);
    gretl_matrix_multiply_by_scalar(tmp, x);
    gretl_matrix_add_to(a, tmp);
}

private void decrement_state_var(GretlMatrixRaw * P, const(GretlMatrixRaw) * K, GretlMatrixRaw * tmp, double x)
{
    gretl_matrix_multiply_mod(K, GretlMatrixMod.none,
                              K, GretlMatrixMod.transpose,
                              tmp, GretlMatrixMod.none);
    gretl_matrix_multiply_by_scalar(tmp, x);
    gretl_matrix_subtract_from(P, tmp);
}

private void state_cross_update(GretlMatrixRaw * Pti, const(GretlMatrixRaw) * Kti, const(GretlMatrixRaw) * Kki, GretlMatrixRaw * tmp, double Fti, double Fkinv)
{
    gretl_matrix_multiply_mod(Kki, GretlMatrixMod.none,
                              Kki, GretlMatrixMod.transpose,
                              tmp, GretlMatrixMod.none);
    gretl_matrix_multiply_by_scalar(tmp, Fkinv * Fti * Fkinv);
    gretl_matrix_add_to(Pti, tmp);
    gretl_matrix_multiply_mod(Kti, GretlMatrixMod.none,
                              Kki, GretlMatrixMod.transpose,
                              tmp, GretlMatrixMod.none);
    gretl_matrix_multiply_mod(Kki, GretlMatrixMod.none,
                              Kti, GretlMatrixMod.transpose,
                              tmp, GretlMatrixMod.cumulate);
    gretl_matrix_multiply_by_scalar(tmp, Fkinv);
    gretl_matrix_subtract_from(Pti, tmp);
}

private int ensure_DK_smo_recorders(Kalman * K, int Nd)
{
    if (K.sinfo.Finf == null) {
        K.sinfo.Finf = gretl_zero_matrix_new(K.n, Nd);
        K.sinfo.Kinf = gretl_zero_matrix_new(K.n * K.r, Nd);
        K.sinfo.Nd = Nd;
    }

    return 0;
}

private int kalman_ldl(Kalman * K)
{
    GretlMatrixRaw * L = K.sinfo.Linv;
    GretlMatrixRaw * g = K.sinfo.g;
    double gj;
    double lij;
    int i;
    int j;
    int err;

    fast_copy_values(L, K.VY);
    err = gretl_matrix_cholesky_decomp(L);

    if (!err) {
        for (j=0; j<K.n; j++) {
            gj = gretl_matrix_get(L, j, j);
            g.val[j] = gj * gj;
            for (i=j; i<K.n; i++) {
                lij = gretl_matrix_get(L, i, j);
                gretl_matrix_set(L, i, j, lij / gj);
            }
        }
        /* record L so we can undo the transformation */
        gretl_bundle_set_matrix(K.b, "L", L);
        err = gretl_invert_triangular_matrix(L, 'L');
    }

    return err;
}

private int kalman_diagonalize(Kalman * K)
{
    GretlMatrixRaw * tmp = gretl_matrix_alloc(K.n, 1);
    GretlMatrixRaw * yt = gretl_matrix_alloc(K.n, 1);
    int t;
    int err;

    err = kalman_ldl(K);
    if (err) {
        gretl_matrix_free(tmp);
        gretl_matrix_free(yt);
        return err;
    }

    /* adjust Z */
    gretl_matrix_multiply_mod(K.sinfo.Linv, GretlMatrixMod.none,
                              K.ZT, GretlMatrixMod.transpose,
                              K.sinfo.Z, GretlMatrixMod.none);

    if (K.sinfo.y0 == null) {
        K.sinfo.y0 = gretl_matrix_copy(K.y);
        if (K.sinfo.y0 == null) {
            err = GretlError.alloc;
        }
    }

    if (!err) {
        /* adjust the observable */
        for (t=0; t<K.N; t++) {
            load_from_row(tmp, K.sinfo.y0, t);
            gretl_matrix_multiply(K.sinfo.Linv, tmp, yt);
            record_to_vec(K.y, yt, t);
        }
    }

    gretl_matrix_free(tmp);
    gretl_matrix_free(yt);

    return err;
}

private int handle_sequential_transforms(Kalman * K)
{
    int err = 0;

    err = gretl_matrix_transpose(K.sinfo.Z, K.ZT);

    if (K.n > 1 && K.VY != null) {
        /* obsvar needs transformation */
        int i;

        if (matrix_is_diagonal(K.VY)) {
            for (i=0; i<K.n; i++) {
                K.sinfo.g.val[i] = gretl_matrix_get(K.VY, i, i);
            }
        } else {
            if (kdebug) {
                fprintf(stderr, "kfilter: diagonalizing\n");
            }
            err = kalman_diagonalize(K);
        }
    }

    return err;
}

private double uni_get_vti(Kalman * K, const(GretlMatrixRaw) * Zi, const(GretlMatrixRaw) * ati, int i)
{
    double vti = gretl_matrix_get(K.y, K.t, i);

    if (!na(vti)) {
        if (K.BT != null) {
            vti -= kalman_Bx_uni(K, i);
        }
        if (!na(vti)) {
            vti -= dotprod(Zi, ati);
        }
    }

    return vti;
}

private double tiny_value(const(GretlMatrixRaw) * Z)
{
    int i;
    int n = Z.rows * Z.cols;
    double x;
    double zmin = 1.0e+30;
    double macheps = 2.20e-16;

    for (i=0; i<n; i++) {
        x = fabs(Z.val[i]);
        if (x > 0 && x < zmin) {
            zmin = x;
        }
    }

    x = sqrt(macheps) * zmin * zmin;
    return x;
}

private int kfilter_sequential(Kalman * K, Prn * prn)
{
    SeqInfo * si = K.sinfo;
    GretlMatrixRaw * Kkt = K.Mt;
    GretlMatrixRaw * at = K.a0;
    GretlMatrixRaw * Pt = K.P0;
    GretlMatrixRaw * ati = K.a1;
    GretlMatrixRaw * Pti = K.P1;
    GretlMatrixRaw * Pk = K.Pk0;
    GretlMatrixRaw * V = K.V;
    GretlMatrixRaw * F = K.F;

    GretlMatrixRaw * Z = si.Z;
    GretlMatrixRaw * g = si.g;
    GretlMatrixRaw * Zi = si.Zi;
    GretlMatrixRaw * Kti = si.Kti;
    GretlMatrixRaw * m1 = si.m1;
    GretlMatrixRaw * mm = si.mm;
    GretlMatrixRaw * Pki = null;
    GretlMatrixRaw * Kki = null;

    /* used only for KFAS comparison */
    GretlMatrixRaw * att = null;
    GretlMatrixRaw * Ptt = null;

    double vti;
    double Fti;
    double llct;
    double Ftinv;
    double Fkinv;
    double Fki = 0;
    double qt;
    double ldt;
    double k_tiny = 0;

    int exact_diffuse = kalman_diffuse(K) && K.exact;
    int rankPk = exact_diffuse ? K.r : 0;
    int d = K.exact ? 0 : -1;
    int j = K.exact ? 0 : -1;
    int Nd = 0;
    int t;
    int i;
    int err = 0;

    if (identify || trace) {
        printf("kfilter_sequential(), diffuse = %d, exact = %d, smo_prep = %d\n",
	       kalman_diffuse(K), K.exact, K.smo_prep);
    }

    if (K.exact) {
        Pki = si.Pki;
        Kki = si.Kki;
    }

    if (K.n > 1) {
        if (matrix_is_varying(K, KMat.zt)) {
            err = kalman_ensure_stepinfo(K, 0);
            if (err) {
                return err;
            }
        }
        /* select the sequential-special V and F matrices */
        if (trace) {
            printf("  select V, F from SeqInfo\n");
        }
        V = si.V;
        F = si.F;
    }

    if (K.n > 1 || K.r > 1) {
        handle_sequential_transforms(K);
    }

    if (K.smo_prep && K.exact) {
        Nd = matrix_is_varying(K, KMat.zt) ? K.N : 4*K.r;
        ensure_DK_smo_recorders(K, Nd);
    }

    if (kalman_extra(K)) {
        att = gretl_matrix_alloc(K.N, K.r);
        Ptt = gretl_matrix_alloc(K.N, K.r * K.r);
    }

    K.TI = gretl_is_identity_matrix(K.T);
    K.okN = K.n * K.N;
    K.qsum = 0;
    K.loglik = 0;
    set_kalman_running(K);

    if (kdebug) {
        fprintf(stderr, "\n*** kfilter_sequential: r=%d, n=%d, exact=%d ***\n",
                K.r, K.n, K.exact);
    }

    if (K.exact) {
        k_tiny = tiny_value(Z);
    }

    for (t=0; t<K.N; t++) {
        K.t = t;
        record_to_vec(K.A, at, t);
        record_to_vech(K.P, Pt, K.r, t);
        fast_copy_values(ati, at);
        fast_copy_values(Pti, Pt);
        if (Pki != null) {
            fast_copy_values(Pki, Pk);
        }
        if (filter_is_varying(K)) {
            /* we have time-varying matrices */
            err = kalman_refresh_matrices(K, prn);
            if (err) {
                K.loglik = NADBL;
                break;
            }
        }
        qt = ldt = llct = 0;
        if (K.smo_prep && d == 0 && t < Nd) {
            record_to_col(K.PK, Pk, t);
        }
        for (i=0; i<K.n; i++) {
            if (kdebug && t < 2 && K.exact) {
                fprintf(stderr, "t,i = %d,%d, rankPk = %d\n", t, i, rankPk);
            }
            load_from_row(Zi, Z, i);
            gretl_matrix_multiply_mod(Pti, GretlMatrixMod.none,
                                      Zi, GretlMatrixMod.transpose,
                                      Kti, GretlMatrixMod.none);
            Fti = dotprod(Zi, Kti);
            if (g != null) {
                Fti += g.val[i];
            }
            gretl_matrix_set(F, t, i, Fti);
            record_to_col(K.Kt, Kti, i);
            if (d == 0) {
                /* still initial */
                gretl_matrix_multiply_mod(Pki, GretlMatrixMod.none,
                                          Zi, GretlMatrixMod.transpose,
                                          Kki, GretlMatrixMod.none);
                Fki = dotprod(Zi, Kki);
                if (K.smo_prep && t < Nd) {
                    gretl_matrix_set(si.Finf, i, t, Fki);
                    record_to_col(Kkt, Kki, i);
                }
            }

            vti = uni_get_vti(K, Zi, ati, i);
            gretl_matrix_set(V, t, i, vti);
            if (na(vti)) {
                K.okN -= 1;
                continue;
            }

            if (rankPk > 0) {
                if (Fki > k_tiny) {
                    K.okN -= 1;
                    Fkinv = 1.0 / Fki;
                    ldt += log(Fki);
                    increment_state(ati, Kki, m1, Fkinv * vti);
                    state_cross_update(Pti, Kti, Kki, mm, Fti, Fkinv);
                    decrement_state_var(Pki, Kki, mm, Fkinv);
                    --rankPk;
                } else if (Fti > 0) {
                    Ftinv = 1.0 / Fti;
                    qt += vti * vti * Ftinv;
                    ldt += log(Fti);
                    llct += LN_2_PI;
                    increment_state(ati, Kti, m1, Ftinv * vti);
                    decrement_state_var(Pti, Kti, mm, Ftinv);
                }
            } else {
                /* rankPk = 0 */
                if (j == 0) {
                    j = (K.n == 1)? 1 : i;
                    if (kdebug) {
                        fprintf(stderr, "*** filter: got j = %d (i=%d)\n", j, i);
                    }
                }
                Ftinv = 1/Fti;
                qt += vti * vti * Ftinv;
                ldt += log(Fti);
                llct += LN_2_PI;
                increment_state(ati, Kti, m1, Ftinv * vti);
                decrement_state_var(Pti, Kti, mm, Ftinv);
            }
        } /* observables at time t */

	if (K.K != null) {
	    /* record the gain */
	    record_to_vec(K.K, K.Kt, t);
	}
        if (K.smo_prep && d == 0 && t < Nd) {
	    /* record K∞ */
	    record_to_col(si.Kinf, Kkt, t);
        }
        if (j > 0 && d == 0) {
            d = t + 1; /* note: 1-based */
            if (kdebug) {
                fprintf(stderr, "*** filter: got d = %d\n", d);
            }
        }

        /* update for t+1 */
        if (K.TI) {
            /* shortcut in case T = I_m */
            fast_copy_values(at, ati);
            fast_copy_values(Pt, Pti);
            gretl_matrix_add_to(Pt, K.VS);
            if (Pki != null) {
                fast_copy_values(Pk, Pki);
            }
        } else {
            gretl_matrix_multiply(K.T, ati, at);
            fast_copy_values(Pt, K.VS);
            gretl_matrix_qform(K.T, GretlMatrixMod.none, Pti,
                               Pt, GretlMatrixMod.cumulate);
            if (Pki != null) {
                gretl_matrix_qform(K.T, GretlMatrixMod.none, Pki,
                                   Pk, GretlMatrixMod.none);
            }
        }
        if (K.mu != null) {
            gretl_matrix_add_to(at, K.mu);
        }

        if (K.LL != null) {
            K.LL.val[t] = -0.5 * (llct + ldt + qt);
            K.loglik += K.LL.val[t];
        } else {
            K.loglik -= 0.5 * (llct + ldt + qt);
        }
        K.qsum += qt;
        if (att != null) {
            record_to_vec(att, ati, t);
        }
        if (Ptt != null) {
            record_to_vec(Ptt, Pti, t);
        }
    } /* end of time loop */

    if (0) { /* Hmm, conditionality? */
        K.d = K.N;
        K.j = K.n;
    } else if (K.smo_prep && K.exact) {
        dj_from_Finf(si.Finf, &d, &j, k_tiny);
    }

    if (kalman_diffuse(K) && K.exact) {
        K.d = d > 0 ? d : 0;
        K.j = j > 0 ? j : 0;
    }

    set_kalman_stopped(K);

    K.s2 = K.qsum / K.okN;

    if (kdebug) {
        fprintf(stderr, "*** after filtering: d=%d, j=%d ***\n", d, j);
        fprintf(stderr, "sequential: loglik = %#0.8g\n", K.loglik);
    }

    if (kalman_extra(K)) {
        if (att != null) {
            bundle_add_matrix(K.b, "att", att);
        }
        if (Ptt != null) {
            bundle_add_matrix(K.b, "Ptt", Ptt);
        }
    }

    if (K.n > 1) {
        if (trace) {
            printf( "  call mv_transform\n");
        }
        mv_transform(K);
    }

    return 0;
}

/* sequential smoothing follows */

/* cumulant matrices */

struct Cumulants {
    GretlMatrixRaw * r0;
    GretlMatrixRaw * N0;
    GretlMatrixRaw * L0;
    GretlMatrixRaw * r1;
    GretlMatrixRaw * N1;
    GretlMatrixRaw * N2;
    GretlMatrixRaw * Nt;
    GretlMatrixRaw * mm;
    GretlMatrixRaw * m1;
};

private int allocate_cumulants(Cumulants * c, Kalman * K)
{
    int r = K.r;

    c.r0 = gretl_zero_matrix_new(r, 1);
    c.N0 = gretl_zero_matrix_new(r, r);
    c.L0 = gretl_zero_matrix_new(r, r);
    c.Nt = gretl_zero_matrix_new(r, r);
    c.mm = gretl_matrix_alloc(r, r);
    c.m1 = gretl_matrix_alloc(r, 1);
    if (K.exact) {
        c.r1 = gretl_zero_matrix_new(r, 1);
        c.N1 = gretl_zero_matrix_new(r, r);
        c.N2 = gretl_zero_matrix_new(r, r);
    }

    return 0;
}

private void clear_cumulants(Cumulants * c, Kalman * K)
{
    gretl_matrix_free(c.r0);
    gretl_matrix_free(c.N0);
    gretl_matrix_free(c.L0);
    gretl_matrix_free(c.Nt);
    gretl_matrix_free(c.mm);
    gretl_matrix_free(c.m1);
    if (K.exact) {
        gretl_matrix_free(c.r1);
        gretl_matrix_free(c.N1);
        gretl_matrix_free(c.N2);
    }
}

/* complex smoothing iteration when F∞ > 0 */

private void fkpos(double fkinv, double Fti, double vti, const(GretlMatrixRaw) * Kki, const(GretlMatrixRaw) * Kti, const(GretlMatrixRaw) * Zti, Cumulants * c)
{
    GretlMatrixRaw * Linf;
    GretlMatrixRaw * Lmid;
    GretlMatrixRaw * mm2;
    int m = Zti.cols;

    mm2 = gretl_matrix_alloc(m, m);

    /* Linf = I(m) - Kki * Zti * fkinv */
    Linf = gretl_identity_matrix_new(m);
    gretl_matrix_multiply(Kki, Zti, c.mm);
    gretl_matrix_multiply_by_scalar(c.mm, fkinv);
    gretl_matrix_subtract_from(Linf, c.mm);

    /* Lmid = Kki * Fti * fkinv - Kti */
    Lmid = gretl_matrix_copy(Kki);
    gretl_matrix_multiply_by_scalar(Lmid, Fti * fkinv);
    gretl_matrix_subtract_from(Lmid, Kti);

    /* L0 = Lmid * Zti * fkinv */
    gretl_matrix_multiply(Lmid, Zti, c.L0);
    gretl_matrix_multiply_by_scalar(c.L0, fkinv);

    /* r1 = Linf' * r1 + L0' * r0 */
    fast_copy_values(c.m1, c.r1);
    gretl_matrix_multiply_mod(Linf, GretlMatrixMod.transpose,
                              c.m1, GretlMatrixMod.none,
                              c.r1, GretlMatrixMod.none);
    gretl_matrix_multiply_mod(c.L0, GretlMatrixMod.transpose,
                              c.r0, GretlMatrixMod.none,
                              c.r1, GretlMatrixMod.cumulate);

    /* r1 += vti * fkinv * Zti' */
    gretl_matrix_copy_values_shaped(c.m1, Zti);
    gretl_matrix_multiply_by_scalar(c.m1, vti * fkinv);
    gretl_matrix_add_to(c.r1, c.m1);

    /* r0 = Linf' * r0 */
    gretl_matrix_multiply_mod(Linf, GretlMatrixMod.transpose,
                              c.r0, GretlMatrixMod.none,
                              c.m1, GretlMatrixMod.none);
    fast_copy_values(c.r0, c.m1);

    /* N2 = Linf' * N2 * Linf */
    gretl_matrix_multiply_mod(Linf, GretlMatrixMod.transpose,
                              c.N2, GretlMatrixMod.none,
                              c.mm, GretlMatrixMod.none);
    gretl_matrix_multiply(c.mm, Linf, c.N2);

    /* N2 -= Fti * fkinv * fkinv * Zti' * Zti */
    gretl_matrix_multiply_mod(Zti, GretlMatrixMod.transpose,
                              Zti, GretlMatrixMod.none,
                              c.mm, GretlMatrixMod.none);
    gretl_matrix_multiply_by_scalar(c.mm, Fti * fkinv * fkinv);
    gretl_matrix_subtract_from(c.N2, c.mm);

    /* N2 += L0' * N0 * L0 */
    gretl_matrix_multiply_mod(c.L0, GretlMatrixMod.transpose,
                              c.N0, GretlMatrixMod.none,
                              c.mm, GretlMatrixMod.none);
    gretl_matrix_multiply(c.mm, c.L0, mm2);
    gretl_matrix_add_to(c.N2, mm2);

    /* mm = Linf' * N1 * L0 */
    gretl_matrix_multiply_mod(Linf, GretlMatrixMod.transpose,
                              c.N1, GretlMatrixMod.none,
                              mm2, GretlMatrixMod.none);
    gretl_matrix_multiply(mm2, c.L0, c.mm);

    /* N2 += mm + mm' */
    gretl_matrix_add_to(c.N2, c.mm);
    gretl_matrix_add_transpose_to(c.N2, c.mm);

    /* N1 = Linf' * N1 * Linf */
    gretl_matrix_multiply(mm2, Linf, c.N1);

    /* N1 += fkinv * Zti' * Zti */
    gretl_matrix_multiply_mod(Zti, GretlMatrixMod.transpose,
                              Zti, GretlMatrixMod.none,
                              c.mm, GretlMatrixMod.none);
    gretl_matrix_multiply_by_scalar(c.mm, fkinv);
    gretl_matrix_add_to(c.N1, c.mm);

    /* N1 += L0' * N0 * Linf */
    gretl_matrix_multiply(c.N0, Linf, c.mm);
    gretl_matrix_multiply_mod(c.L0, GretlMatrixMod.transpose,
                              c.mm, GretlMatrixMod.none,
                              c.N1, GretlMatrixMod.cumulate);

    /* N0 = Linf' * N0 * Linf */
    gretl_matrix_multiply_mod(Linf, GretlMatrixMod.transpose,
                              c.N0, GretlMatrixMod.none,
                              c.mm, GretlMatrixMod.none);
    gretl_matrix_multiply(c.mm, Linf, c.N0);

    gretl_matrix_free(mm2);
    gretl_matrix_free(Linf);
    gretl_matrix_free(Lmid);
}

/* relatively simple smoothing iteration when F∞ = 0 */

private void fkzero(double ftinv, double vti, const(GretlMatrixRaw) * Kti, const(GretlMatrixRaw) * Zti, Cumulants * c)
{
    /* L0 = I(m) - ftinv * Kti * Zti */
    fast_write_I(c.L0);
    gretl_matrix_multiply(Kti, Zti, c.mm);
    gretl_matrix_multiply_by_scalar(c.mm, ftinv);
    gretl_matrix_subtract_from(c.L0, c.mm);

    /* r0 = L0' * r0 + vti * ftinv * Zti' */
    gretl_matrix_multiply_mod(c.L0, GretlMatrixMod.transpose,
                              c.r0, GretlMatrixMod.none,
                              c.m1, GretlMatrixMod.none);
    fast_copy_values(c.r0, c.m1);
    gretl_matrix_copy_values_shaped(c.m1, Zti);
    gretl_matrix_multiply_by_scalar(c.m1, vti * ftinv);
    gretl_matrix_add_to(c.r0, c.m1);

    if (c.r1 != null) {
        /* r1 = L0' * r1 */
        fast_copy_values(c.m1, c.r1);
        gretl_matrix_multiply_mod(c.L0, GretlMatrixMod.transpose,
                                  c.m1, GretlMatrixMod.none,
                                  c.r1, GretlMatrixMod.none);
    }

    /* N0 = L0' * N0 * L0 + ftinv * Z[i,]' * Z[i,] */
    gretl_matrix_multiply_mod(c.L0, GretlMatrixMod.transpose,
                              c.N0, GretlMatrixMod.none,
                              c.mm, GretlMatrixMod.none);
    gretl_matrix_multiply(c.mm, c.L0, c.N0);
    gretl_matrix_multiply_mod(Zti, GretlMatrixMod.transpose,
                              Zti, GretlMatrixMod.none,
                              c.mm, GretlMatrixMod.none);
    gretl_matrix_multiply_by_scalar(c.mm, ftinv);
    gretl_matrix_add_to(c.N0, c.mm);

    if (c.N1 != null && c.N2 != null) {
        /* N1 = N1 * L0 */
        gretl_matrix_multiply(c.N1, c.L0, c.mm);
        fast_copy_values(c.N1, c.mm);
        /* N2 = N2 * L0 */
        gretl_matrix_multiply(c.N2, c.L0, c.mm);
        fast_copy_values(c.N2, c.mm);
    }
}

/* smoothing update of state and its variance */

private int dagger_calc(Cumulants * c, GretlMatrixRaw * Pt, GretlMatrixRaw * Pk, Kalman * K)
{
    GretlMatrixBlock * B;
    GretlMatrixRaw * rdag;
    GretlMatrixRaw * Ndag;
    GretlMatrixRaw * Pdag;
    GretlMatrixRaw * at = c.m1;
    GretlMatrixRaw * tmp;
    int m = c.N0.rows;
    int m2 = m * 2;
    int mm = m * m;
    size_t rsize = m * double.sizeof;
    size_t Nsize = mm * double.sizeof;

    /* workspace */
    B = gretl_matrix_block_new(&tmp, m, m2, &rdag, m2, 1,
                               &Ndag, m2, m2, &Pdag, m, m2,
                               null);
    if (B == null) {
        return GretlError.alloc;
    }

    /* rdag = r0 | r1 */
    memcpy(rdag.val, c.r0.val, rsize);
    memcpy(rdag.val + m, c.r1.val, rsize);

    /* Ndag = (N0 ~ N1') | (N1 ~ N2) */
    gretl_matrix_inscribe_matrix(Ndag, c.N0, 0, 0, GretlMatrixMod.none);
    gretl_matrix_inscribe_matrix(Ndag, c.N1, 0, m, GretlMatrixMod.transpose);
    gretl_matrix_inscribe_matrix(Ndag, c.N1, m, 0, GretlMatrixMod.none);
    gretl_matrix_inscribe_matrix(Ndag, c.N2, m, m, GretlMatrixMod.none);

    load_from_row(at, K.A, K.t);
    load_from_vech(Pt, K.P, K.r, K.t);
    load_from_col(Pk, K.PK, K.t);

    /* Pdag = Pt ~ Pk */
    memcpy(Pdag.val, Pt.val, Nsize);
    memcpy(Pdag.val + mm, Pk.val, Nsize);

    if (kdebug) {
        fprintf(stderr, "dagger calc, t=%d\n", K.t);
        gretl_matrix_print(Pt, "Pt");
        gretl_matrix_print(Pdag, "Pdag");
        gretl_matrix_print(at, "at");
    }

    /* ahat[t,] = (at + Pdag * rdag)' */
    gretl_matrix_multiply_mod(Pdag, GretlMatrixMod.none,
                              rdag, GretlMatrixMod.none,
                              at, GretlMatrixMod.cumulate);
    record_to_vec(K.A, at, K.t);

    /* Vhat[t,] = vec(Pt - Pdag * Ndag * Pdag')' */
    gretl_matrix_multiply(Pdag, Ndag, tmp);
    gretl_matrix_multiply_mod(tmp, GretlMatrixMod.none,
                              Pdag, GretlMatrixMod.transpose,
                              Pt, GretlMatrixMod.decrement);
    record_to_vech(K.P, Pt, K.r, K.t);

    gretl_matrix_block_destroy(B);

    return 0;
}

private void eps_smooth_Ft(double ftinv, double gi, double vti, const(GretlMatrixRaw) * Kti, Cumulants * c, GretlMatrixRaw * epshat, GretlMatrixRaw * veps, int t, int i, int dist)
{
    double x = dotprod(Kti, c.r0);

    gretl_matrix_set(epshat, t, i, gi * ftinv * (vti - x));

    gretl_matrix_multiply(c.N0, Kti, c.m1);
    gretl_matrix_multiply_by_scalar(c.m1, ftinv * ftinv);
    if (dist == 2) {
        x = gi - gi*gi * (ftinv + dotprod(Kti, c.m1));
    } else {
        x = gi*gi * (ftinv + dotprod(Kti, c.m1));
    }
    gretl_matrix_set(veps, t, i, x);
}

private void eps_smooth_Fk(double fkinv, double gi, const(GretlMatrixRaw) * Kki, Cumulants * c, GretlMatrixRaw * epshat, GretlMatrixRaw * veps, int t, int i, int dist)
{
    double x = dotprod(Kki, c.r0);

    gretl_matrix_set(epshat, t, i, -gi * x * fkinv);

    gretl_matrix_multiply(c.N0, Kki, c.m1);
    if (dist == 2) {
        x = gi - gi*gi * fkinv * fkinv * dotprod(Kki, c.m1);
    } else {
        x = gi*gi * fkinv * fkinv * dotprod(Kki, c.m1);
    }
    gretl_matrix_set(veps, t, i, x);
}

private void eta_smooth(const(GretlMatrixRaw) * Q, const(GretlMatrixRaw) * QRT, Cumulants * c, GretlMatrixRaw * etahat, GretlMatrixRaw * veta, int t, int dist)
{
    GretlMatrixRaw * tmp;
    int q = Q.rows;
    int m = QRT.cols;

    tmp = gretl_matrix_reuse(c.mm, q, q);

    if (dist == 2) {
        /* Q - QRT * Nt * (QRT)' */
        fast_copy_values(tmp, Q);
        gretl_matrix_qform(QRT, GretlMatrixMod.none, c.Nt,
                           tmp, GretlMatrixMod.decrement);
    } else {
        /* QRT * Nt * (QRT)' */
        gretl_matrix_qform(QRT, GretlMatrixMod.none, c.Nt,
                           tmp, GretlMatrixMod.none);
    }

    record_to_vec(veta, tmp, t);

    if (t > 0) {
        /* QRT * r */
        gretl_matrix_reuse(tmp, -1, 1);
        gretl_matrix_multiply(QRT, c.r0, tmp);
        record_to_vec(etahat, tmp, t-1);
    }

    gretl_matrix_reuse(c.mm, m, m);
}

private void state_smooth(const(GretlMatrixRaw) * at, const(GretlMatrixRaw) * Pt, Cumulants * c, Kalman * K)
{
    /* Pt - Pt * N0 * Pt' */
    fast_copy_values(c.mm, Pt);
    gretl_matrix_qform(Pt, GretlMatrixMod.none, c.N0,
                       c.mm, GretlMatrixMod.decrement);
    record_to_vech(K.P, c.mm, K.r, K.t);

    /* at + Pt * r0 */
    gretl_matrix_multiply(Pt, c.r0, c.m1);
    gretl_matrix_add_to(c.m1, at);
    record_to_vec(K.A, c.m1, K.t);
}

private void regular_backdate(Cumulants * c, GretlMatrixRaw * T)
{
    /* N0 = T' * N0 * T */
    gretl_matrix_multiply_mod(T, GretlMatrixMod.transpose,
                              c.N0, GretlMatrixMod.none,
                              c.mm, GretlMatrixMod.none);
    gretl_matrix_multiply(c.mm, T, c.N0);

    /* r0 = T' * r0 */
    gretl_matrix_multiply_mod(T, GretlMatrixMod.transpose,
                              c.r0, GretlMatrixMod.none,
                              c.m1, GretlMatrixMod.none);
    fast_copy_values(c.r0, c.m1);
}

private void diffuse_backdate(Cumulants * c, GretlMatrixRaw * T)
{
    /* N0 = T' * N0 * T */
    gretl_matrix_multiply_mod(T, GretlMatrixMod.transpose,
                              c.N0, GretlMatrixMod.none,
                              c.mm, GretlMatrixMod.none);
    gretl_matrix_multiply(c.mm, T, c.N0);

    /* N1 = T' * N1 * T */
    gretl_matrix_multiply_mod(T, GretlMatrixMod.transpose,
                              c.N1, GretlMatrixMod.none,
                              c.mm, GretlMatrixMod.none);
    gretl_matrix_multiply(c.mm, T, c.N1);

    /* N2 = T' * N2 * T */
    gretl_matrix_multiply_mod(T, GretlMatrixMod.transpose,
                              c.N2, GretlMatrixMod.none,
                              c.mm, GretlMatrixMod.none);
    gretl_matrix_multiply(c.mm, T, c.N2);

    /* r0 = T' * r0 */
    gretl_matrix_multiply_mod(T, GretlMatrixMod.transpose,
                              c.r0, GretlMatrixMod.none,
                              c.m1, GretlMatrixMod.none);
    fast_copy_values(c.r0, c.m1);
    /* r1 = T' * r1 */
    gretl_matrix_multiply_mod(T, GretlMatrixMod.transpose,
                              c.r1, GretlMatrixMod.none,
                              c.m1, GretlMatrixMod.none);
    fast_copy_values(c.r1, c.m1);
}

private void kalman_bundle_add_epshat(Kalman * K, GretlMatrixRaw * epshat)
{
    GretlMatrixRaw * L;

    /* add epshat "as is" */
    bundle_add_matrix(K.b, "epshat", epshat);

    L = gretl_bundle_get_matrix(K.b, "L", null);

    if (L != null) {
        /* The triangular matrix L will be present only when the
           variance matrix of eps is not diagonal. In that case
           we'll try adding estimates of the untransformed
           disturbances.
        */
        GretlMatrixRaw * real_epshat;
        int T = epshat.rows;
        int k = epshat.cols;

        real_epshat = gretl_matrix_alloc(T, k);

        if (real_epshat != null) {
            double eti;
            double lij;
            double etj;
            int i;
            int j;
            int jmax;
            int t;

            /* real_epshat[t,] = (L * epshat[t,]')' */
            for (t=0; t<T; t++) {
                jmax = 1;
                for (i=0; i<k; i++) {
                    eti = 0.0;
                    for (j=0; j<jmax; j++) {
                        lij = gretl_matrix_get(L, i, j);
                        etj = gretl_matrix_get(epshat, t, j);
                        eti += lij * etj;
                    }
                    gretl_matrix_set(real_epshat, t, i, eti);
                    jmax++;
                }
            }
            bundle_add_matrix(K.b, "real_epshat", real_epshat);
        }
    }
}

private int state_dist_setup(Kalman * K, GretlMatrixRaw ** Q, GretlMatrixRaw ** QRT, GretlMatrixRaw ** etahat, GretlMatrixRaw ** veta)
{
    int q = K.r;

    if (K.Q != null) {
        q = K.Q.rows;
        *Q = K.Q;
        *QRT = K.QRT;
    } else {
        *Q = K.VS;
        *QRT = K.VS;
    }

    *etahat = gretl_zero_matrix_new(K.N, q);
    *veta = gretl_zero_matrix_new(K.N, q*q);

    return 1;
}

private int ksmooth_sequential(Kalman * K, int dist)
{
    GretlMatrixBlock * B = null;
    SeqInfo * si = K.sinfo;
    Cumulants c;
    GretlMatrixRaw * Q = null;
    GretlMatrixRaw * QRT = null;
    GretlMatrixRaw * V = null;
    GretlMatrixRaw * F = null;
    GretlMatrixRaw * g = null;

    int i;
    int t;
    int N = K.N;
    int d = K.d;
    int j = K.j;
    int eps_smo = 0;
    int eta_smo = 0;
    int err = 0;

    /* matrices to be returned in bundle */
    GretlMatrixRaw * epshat = null;
    GretlMatrixRaw * veps = null;
    GretlMatrixRaw * etahat = null;
    GretlMatrixRaw * veta = null;

    /* block-workspace matrices */
    GretlMatrixRaw * Pt;
    GretlMatrixRaw * Pk;
    GretlMatrixRaw * vt;
    GretlMatrixRaw * Ft;
    GretlMatrixRaw * Kt;
    GretlMatrixRaw * at;
    GretlMatrixRaw * Fk;
    GretlMatrixRaw * Kk;
    GretlMatrixRaw * Zti;
    GretlMatrixRaw * Kti;
    GretlMatrixRaw * Kki;

    /* workspace scalars */
    double vti;
    double Fti;
    double Fki;
    double ftinv;
    double fkinv;
    double k_tiny = 0;

    /* convenience pointers */
    V = K.n > 1 ? si.V : K.V;
    F = K.n > 1 ? si.F : K.F;
    g  = si.g;

    if (trace) {
        printf("ksmooth_sequential(), dist = %d, d=%d, j=%d\n",
               dist, K.d, K.j);
    }
    if (kdebug) {
        fprintf(stderr, "ksmooth_sequential: dist=%d, d=%d, j=%d\n",
                dist, K.d, K.j);
    }

    if (K.exact) {
        B = gretl_matrix_block_new(&Kt, K.r, K.n,
				   &Pt, K.r, K.r,
                                   &Ft, K.n, 1,
				   &vt, K.n, 1,
                                   &at, K.r, 1,
				   &Fk, K.n, 1,
                                   &Kk, K.r, K.n,
				   &Pk, K.r, K.r,
                                   &Zti, 1, K.r,
				   &Kti, K.r, 1,
                                   &Kki, K.r, 1, null);
    } else {
        B = gretl_matrix_block_new(&Kt, K.r, K.n,
				   &Pt, K.r, K.r,
                                   &Ft, K.n, 1,
				   &vt, K.n, 1,
                                   &at, K.r, 1,
				   &Zti, 1, K.r,
                                   &Kti, K.r, 1, null);
    }

    allocate_cumulants(&c, K);

    if (dist) {
        if (g != null) {
            /* smoothing of obs disturbances wanted */
            epshat = gretl_zero_matrix_new(N, K.n);
            veps = gretl_zero_matrix_new(N, K.n);
            eps_smo = 1;
        }
        /* smoothing of state disturbances wanted */
        eta_smo = state_dist_setup(K, &Q, &QRT, &etahat, &veta);
    }

    if (K.exact) {
	k_tiny = tiny_value(si.Z); /* ? */
    }

    if (kdebug) {
        fprintf(stderr, "\n*** sequential smoother: d=%d, j=%d, dist=%d ***\n",
                d, j, dist);
    }

    for (t=K.N-1; t>=d && !err; t--) {
        K.t = t;
        if (matrix_is_varying(K, KMat.t)) {
            err = retrieve_Tt(K);
        }
        if (!err && matrix_is_varying(K, KMat.zt)) {
            err = retrieve_Zt(K);
        }
        load_from_row(vt, V, t);
        load_from_row(Ft, F, t);
        load_from_row(Kt, K.K, t);
        load_from_row(at, K.A, t);
        load_from_vech(Pt, K.P, K.r, K.t);

        /* loop across observables */
        for (i=K.n-1; i>=0; i--) {
            vti = vt.val[i];
            Fti = Ft.val[i];
            if (na(vti) || Fti == 0) {
                if (eps_smo) {
                    gretl_matrix_set(veps, t, i, g.val[i]);
                }
                continue;
            }
            ftinv = 1.0 / Fti;
            load_from_row(Zti, si.Z, i);
            load_from_col(Kti, Kt, i);
            if (eps_smo) {
                eps_smooth_Ft(ftinv, g.val[i], vti, Kti, &c,
                              epshat, veps, t, i, dist);
            }
            fkzero(ftinv, vti, Kti, Zti, &c);
        }

        if (eta_smo) {
            eta_smooth(Q, QRT, &c, etahat, veta, t, dist);
        }

        /* smoothed state and its variance */
        state_smooth(at, Pt, &c, K);

        /* regular backdate for t-1 */
        if (t > 0) {
            fast_copy_values(c.Nt, c.N0);
            if (!K.TI) {
                regular_backdate(&c, K.T);
            }
        }
    } /* end first time loop */

    if (d > 0) {
        K.t = t = d - 1; /* first time-step in diffuse phase */
        if (matrix_is_varying(K, KMat.t)) {
            err = retrieve_Tt(K);
        }
        if (!err && matrix_is_varying(K, KMat.zt)) {
            err = retrieve_Zt(K);
        }
        load_from_row(Kt, K.K, t);
        load_from_row(Ft, F, t);
        load_from_row(vt, V, t);

        /* countdown to observable @j */
        for (i=K.n-1; i>=j; i--) {
            vti = vt.val[i];
            Fti = Ft.val[i];
            if (na(vti) || Fti == 0) {
                if (eps_smo) {
                    gretl_matrix_set(veps, t, i, g.val[i]);
                }
                continue;
            }
            ftinv = 1.0 / Fti;
            load_from_row(Zti, si.Z, i);
            load_from_col(Kti, Kt, i);
            if (eps_smo) {
                eps_smooth_Ft(ftinv, g.val[i], vti, Kti, &c,
                              epshat, veps, t, i, dist);
            }
            fkzero(ftinv, vti, Kti, Zti, &c);
        }

        /* retrieve the diffuse values */
        load_from_col(Fk, si.Finf, t);
        load_from_col(Kk, si.Kinf, t);
        load_from_col(Pk, K.PK, t);

        /* further countdown to first observable */
        for (i=j-1; i>=0; i--) {
            Fki = Fk.val[i];
            Fti = Ft.val[i];
            vti = vt.val[i];
            if (kdebug) {
                fprintf(stderr, "t,i = %d,%d, i-countdown (2) j to 0, " ~
                        "Fki = %g\n", t, i, Fki);
            }
            if (na(vti)) {
                if (eps_smo) {
                    gretl_matrix_set(veps, t, i, g.val[i]);
                }
                continue;
            }
            if (Fki > k_tiny) {
                fkinv = 1.0 / Fki;
                load_from_row(Zti, si.Z, i);
                load_from_col(Kki, Kk, i);
                load_from_col(Kti, Kt, i);
                if (eps_smo) {
                    eps_smooth_Fk(fkinv, g.val[i], Kki, &c,
                                  epshat, veps, t, i, dist);
                }
                fkpos(fkinv, Fti, vti, Kki, Kti, Zti, &c);
            } else if (Fti > 0) {
                ftinv = 1.0 / Fti;
                load_from_row(Zti, si.Z, i);
                load_from_col(Kti, Kt, i);
                if (eps_smo) {
                    eps_smooth_Ft(ftinv, g.val[i], vti, Kti, &c,
                                  epshat, veps, t, i, dist);
                }
                fkzero(ftinv, vti, Kti, Zti, &c);
            }
        }

        if (eta_smo) {
            eta_smooth(Q, QRT, &c, etahat, veta, t, dist);
        }
        if (kdebug) {
            fprintf(stderr, "dagger calc, call 1 (t=%d)\n", t);
        }
        dagger_calc(&c, Pt, Pk, K);

        if (t > 0) {
            if (kdebug) {
                fprintf(stderr, "call diffuse_backdate\n");
            }
            fast_copy_values(c.Nt, c.N0);
            if (!K.TI) {
                diffuse_backdate(&c, K.T);
            }
        }

        for (t=d-2; t>=0 && !err; t--) {
            K.t = t;
            // printf("t,i = %d,%d, t-countdown d-2 to 0\n", t, i);
            if (matrix_is_varying(K, KMat.t)) {
                err = retrieve_Tt(K);
            }
            if (!err && matrix_is_varying(K, KMat.zt)) {
                err = retrieve_Zt(K);
            }
            load_from_row(vt, V, t);
            load_from_row(Ft, F, t);
            load_from_row(Kt, K.K, t);
            load_from_col(Fk, si.Finf, t);
            load_from_col(Kk, si.Kinf, t);

            for (i=K.n-1; i>=0; i--) {
                if (kdebug) {
                    fprintf(stderr, "t,i = %d,%d, i-countdown p to 0\n", t, i);
                }
                vti = vt.val[i];
                Fti = Ft.val[i];
                Fki = Fk.val[i];
                if (na(vti)) {
                    if (eps_smo) {
                        gretl_matrix_set(veps, t, i, g.val[i]);
                    }
                    continue;
                }
                if (Fki > k_tiny) {
                    fkinv = 1.0 / Fki;
                    load_from_row(Zti, si.Z, i);
                    load_from_col(Kki, Kk, i);
                    load_from_col(Kti, Kt, i);
                    if (eps_smo) {
                        eps_smooth_Fk(fkinv, g.val[i], Kki, &c,
                                      epshat, veps, t, i, dist);
                    }
                    fkpos(fkinv, Fti, vti, Kki, Kti, Zti, &c);
                } else if (Fti > 0) {
                    ftinv = 1.0 / Fti;
                    load_from_row(Zti, si.Z, i);
                    load_from_col(Kti, Kt, i);
                    if (eps_smo) {
                        eps_smooth_Ft(ftinv, g.val[i], vti, Kti, &c,
                                      epshat, veps, t, i, dist);
                    }
                    fkzero(ftinv, vti, Kti, Zti, &c);
                }
            }

            if (kdebug) {
                fprintf(stderr, "dagger calc, call 2 (t=%d)\n", t);
            }
            dagger_calc(&c, Pt, Pk, K);

            if (eta_smo) {
                eta_smooth(Q, QRT, &c, etahat, veta, t, dist);
            }

            if (t > 0) {
                fast_copy_values(c.Nt, c.N0);
                if (!K.TI) {
                    diffuse_backdate(&c, K.T);
                }
            }
        } /* end t=d-2..0 */
    }

    /* add smoothed quantities to @b */
    if (eps_smo) {
        kalman_bundle_add_epshat(K, epshat);
        bundle_add_matrix(K.b, "veps",  veps);
    }
    if (eta_smo) {
        bundle_add_matrix(K.b, "etahat", etahat);
        bundle_add_matrix(K.b, "veta", veta);
    }

    gretl_matrix_block_destroy(B);
    clear_cumulants(&c, K);

    return 0;
}

/* end of sequential smoother functions */

/* Start de Jong filtering and smoothing functions.  These are based
   on the 2003 paper by de Jong and Lin, "Smoothing with an unknown
   initial condition" (Journal of Time Series Analysis 24:2). See also
   de Jong's original paper on this topic, "The diffuse Kalman filter"
   (Annals of Statistics, June 1991).
*/

private int ensure_extra_dj_smo_recorders(Kalman * K, int Nd)
{
    if (K.exact && K.djinfo.A == null) {
        K.djinfo.A = gretl_zero_matrix_new(K.r * K.r, Nd);
        K.djinfo.V = gretl_zero_matrix_new(K.n * K.r, Nd);
        if (K.djinfo.A == null || K.djinfo.V == null) {
            return GretlError.alloc;
        }
    }

    K.djinfo.Nd = Nd;

    return 0;
}

private int ensure_basic_dj_smo_recorder(Kalman * K)
{
    int err = 0;

    if (K.L == null) {
        K.L = gretl_zero_matrix_new(K.r * K.r, K.N);
        if (K.L == null) {
            err = GretlError.alloc;
        }
    }
    if (!err && K.J == null && K.smo_prep == KSmoType.dist) {
        K.J = gretl_zero_matrix_new(K.r * K.p, K.N);
        if (K.J == null) {
            err = GretlError.alloc;
        }
    }

    return err;
}

private int ensure_dj_variance_matrices(Kalman * K)
{
    GretlMatrixRaw * tmp;
    int err = 0;

    if (K.p == 0) {
        K.p = K.r;
        if (K.VY != null) {
            K.p += K.n;
        }
    }

    if (K.H == null) {
        /* create H from K.VS */
        tmp = gretl_matrix_copy(K.VS);
        if (tmp == null) {
            err = GretlError.alloc;
        } else {
            err = gretl_matrix_psd_root(tmp, 0);
        }
        if (!err) {
	    int offset = K.VY != null ? K.n : 0;

            K.H = gretl_zero_matrix_new(K.r, K.p);
            err = gretl_matrix_inscribe_matrix(K.H, tmp, 0, offset, GretlMatrixMod.none);
        }
	gretl_matrix_free(tmp);
    }

    if (!err && K.G == null && K.VY != null) {
        /* create G from K.VY */
        tmp = gretl_matrix_copy(K.VY);
        if (tmp == null) {
            err = GretlError.alloc;
        } else {
            err = gretl_matrix_psd_root(tmp, 0);
        }
        if (!err) {
            K.G = gretl_zero_matrix_new(K.n, K.p);
            err = gretl_matrix_inscribe_matrix(K.G, tmp, 0, 0, GretlMatrixMod.none);
        }
	gretl_matrix_free(tmp);
    }

    if (kalman_djvar(K) && !err && K.HG == null && K.H != null && K.G != null) {
        tmp = gretl_matrix_alloc(K.r, K.n);
        if (tmp == null) {
            err = GretlError.alloc;
        } else {
            gretl_matrix_multiply_mod(K.H, GretlMatrixMod.none,
                                      K.G, GretlMatrixMod.transpose,
                                      tmp, GretlMatrixMod.none);
            if (gretl_is_zero_matrix(tmp)) {
                gretl_matrix_free(tmp);
            } else {
                K.HG = tmp;
            }
        }
    }

    if (!err && K.Jt == null) {
        K.Jt = gretl_zero_matrix_new(K.r, K.p);
    }

    if (!err && K.VS0 == null && !kalman_djvar(K)) {
	/* establish baseline for detecting changes */
	K.VS0 = gretl_matrix_copy(K.VS);
	if (K.VY != null) {
	    K.VY0 = gretl_matrix_copy(K.VY);
	}
    }

    return err;
}

/* The following function assumes that K.H (and K.G, if present) are
   already allocated at appropriate sizes. We come here only if the
   caller has not used the de Jong representation of the disturbance
   variances, in which case BFGS will be tweaking K.VS rather than
   K.H, so that K.H must be updated accordingly. And similarly for
   K.VY and K.G if applicable.
*/

private int update_dj_variance_matrices(Kalman * K)
{
    int err = 0;

    if (matrix_changed(K.VS, K.VS0)) {
	/* update H from new K.VS */
	fast_copy_values(K.VS0, K.VS);
	err = gretl_matrix_psd_root(K.VS0, 0);
	if (!err) {
	    int offset = K.VY != null ? K.n : 0;

	    err = gretl_matrix_inscribe_matrix(K.H, K.VS0, 0, offset, GretlMatrixMod.none);
	}
	fast_copy_values(K.VS0, K.VS);
    }

    if (err || K.VY == null) {
	return err;
    }

    if (matrix_changed(K.VY, K.VY0)) {
        /* update G from new K.VY */
	fast_copy_values(K.VY0, K.VY);
	err = gretl_matrix_psd_root(K.VY0, 0);
        if (!err) {
            err = gretl_matrix_inscribe_matrix(K.G, K.VY0, 0, 0, GretlMatrixMod.none);
        }
	fast_copy_values(K.VY0, K.VY);
    }

    return err;
}

private int dejong_diffuse_filter_step(Kalman * K, double * pqt, double * pldt)
{
    DejongInfo * dj = K.djinfo;
    GretlMatrixRaw * lam;
    double lmin = 0;
    int err = 0;

    /* Vt = -Zt * At */
    gretl_matrix_zero(dj.Vt);
    gretl_matrix_multiply_mod(K.ZT, GretlMatrixMod.transpose,
                              dj.At, GretlMatrixMod.none,
                              dj.Vt, GretlMatrixMod.decrement);
    if (dj.V != null) {
        record_to_col(dj.V, dj.Vt, K.t);
    }

    /* Vv = Vt ~ vt */
    gretl_matrix_inscribe_matrix(dj.Vv, dj.Vt,
                                 0, 0, GretlMatrixMod.none);
    gretl_matrix_inscribe_matrix(dj.Vv, K.vt,
                                 0, K.r, GretlMatrixMod.none);

    /* Q[t+1] : Qt += Vv' * iFt * Vv */
    gretl_matrix_qform(dj.Vv, GretlMatrixMod.transpose,
                       K.iFt, dj.Qt, GretlMatrixMod.cumulate);
    if (kdebug > 1) {
	gretl_matrix_print(dj.Qt, "Qt");
    }

    /* A[t+1]: At = T*At + Kt*Vt */
    gretl_matrix_multiply(K.Kt, dj.Vt, K.rr);
    if (K.TI) {
        gretl_matrix_add_to(dj.At, K.rr);
    } else {
        gretl_matrix_multiply_mod(K.T, GretlMatrixMod.none,
                                  dj.At, GretlMatrixMod.none,
                                  K.rr, GretlMatrixMod.cumulate);
        fast_copy_values(dj.At, K.rr);
    }

    /* Sinv = Qt[1:r,1:r] */
    gretl_matrix_extract_matrix(dj.S, dj.Qt, 0, 0, GretlMatrixMod.none);
    if (K.r == 1) {
        lmin = dj.S.val[0];
    } else {
        lam = gretl_symmetric_matrix_eigenvals(dj.S, 0, &err);
        lmin = lam.val[0];
        gretl_matrix_free(lam);
    }
    if (kdebug > 1) {
	fprintf(stderr, "dejong_diffuse_filter_step: t=%d, lmin=%g\n", K.t, lmin);
    }

    if (na(lmin)) {
	err = GretlError.nan;
    } else if (lmin > 1.0e-7) {
        /* handle the collapse */
        K.d = K.t + 1;
	if (kdebug) {
	    fprintf(stderr, "*** determined m = %d ***\n", K.d);
	}

        /* Sinv = Qt[1:r,1:r] */
        gretl_matrix_extract_matrix(dj.S, dj.Qt, 0, 0, GretlMatrixMod.none);

        /* -s = Qt[1:r,r+1] */
        gretl_matrix_extract_matrix(dj.s, dj.Qt, 0, K.r, GretlMatrixMod.none);
        gretl_matrix_multiply_by_scalar(dj.s, -1.0);

        /* S = invpd(Sinv) */
        gretl_invert_symmetric_matrix2(dj.S, &dj.ldS);
        *pldt = dj.ldS;

        /* at += At*S*s */
        gretl_matrix_multiply(dj.At, dj.S, K.rr);
        gretl_matrix_multiply_mod(K.rr, GretlMatrixMod.none,
                                  dj.s, GretlMatrixMod.none,
                                  K.a0, GretlMatrixMod.cumulate);

        /* Pt += At*S*At' */
        gretl_matrix_qform(dj.At, GretlMatrixMod.none,
                           dj.S, K.P0, GretlMatrixMod.cumulate);

        /* qm = Qt[r+1,r+1] - s'*S*s */
        dj.qm = gretl_matrix_get(dj.Qt, K.r, K.r);
        dj.qm -= gretl_scalar_qform(dj.s, dj.S, &err); /* FIXME? */
        *pqt = dj.qm;

        /* record */
	if (dj.Am == null) {
	    dj.Am = gretl_matrix_copy(dj.At);
	} else {
	    fast_copy_values(dj.Am, dj.At);
	}
    }

    return err;
}

/* Carry out some calculations that are required when the diffuse
   initial recursion never collapsed to the regular Kalman
   filter. @sum_preldet is the sum of log determinants of F_t, t = 1 to
   N, all of which observations were in the diffuse phase.
*/

private int handle_no_collapse_case(Kalman * K, double nl2pi, double sum_preldet)
{
    DejongInfo * dj = K.djinfo;
    double cterm = K.okN * nl2pi;
    double qm;
    double qadj;
    int err = 0;

    qm = gretl_matrix_get(dj.Qt, K.r, K.r);
    gretl_matrix_extract_matrix(dj.s, dj.Qt, 0, K.r, GretlMatrixMod.none);
    gretl_matrix_multiply_by_scalar(dj.s, -1.0);
    err = gretl_matrix_moore_penrose(dj.S, 1.0e-7);
    if (!err) {
        qadj = gretl_scalar_qform(dj.s, dj.S, &err);
    }
    if (!err) {
        K.loglik = -0.5 * (sum_preldet + qm - qadj + cterm);
    }
    if (dj.Am != null) {
        gretl_matrix_free(dj.Am);
	dj.Am = null;
    }

    return err;
}

private int kfilter_dejong(Kalman * K, Prn * prn)
{
    double nl2pi = K.n * LN_2_PI;
    double sumldet = 0;
    double sum_preldet = 0;
    double llt;
    double llct;
    double ldt;
    double qt;
    double pldt;
    int exact_diffuse;
    int loglik_1991 = 0;
    int nt;
    int Nd = K.N;
    int err = 0;

    exact_diffuse = kalman_diffuse(K) && K.exact;

    if (identify || trace) {
	printf("kfilter_dejong() exact_diffuse = %d\n", exact_diffuse);
    }

    if (K.dj_initted == 0) {
	err = ensure_dj_variance_matrices(K);
	if (err) {
	    return err;
	}
	K.dj_initted = 1;
    }

    if (K.smo_prep) {
        err = ensure_Finv(K);
        if (!err) {
            ensure_basic_dj_smo_recorder(K);
        }
        if (!err && K.exact) {
            Nd = K.N; // ? matrix_is_varying(K, KMat.zt) ? K.N : 4*K.r;
            err = ensure_extra_dj_smo_recorders(K, Nd);
        }
        if (err) {
            return err;
        }
    }

    if (exact_diffuse) {
	/* (re-)initialize djinfo */
        fast_write_I(K.djinfo.At);
	gretl_matrix_zero(K.djinfo.Qt);
	gretl_matrix_zero(K.djinfo.Vt);
    }

    if (kdebug) {
        fprintf(stderr, "\n*** kfilter_dejong: diffuse %d, exact %d, smo %d, djvar %d\n",
                kalman_diffuse(K) ? 1 : 0, K.exact, K.smo_prep ? 1 : 0, kalman_djvar(K));
    }

    if (exact_diffuse) {
        K.d = K.N;
    } else {
        K.d = -1;
    }
    K.qsum = K.loglik = 0.0;
    K.s2 = NADBL;
    K.TI = gretl_is_identity_matrix(K.T);
    K.okN = K.N;

    if (!kalman_djvar(K) && (gretl_iteration_depth() > 0 || numgrad_in_progress())) {
	/* guard against "back door" changes in dist. variance parameters under mle */
        err = update_dj_variance_matrices(K);
    }

    set_kalman_running(K);

    for (K.t = 0; K.t < K.N && !err; K.t += 1) {
        llct = ldt = qt = 0;
        llt = NADBL;

        if (K.A != null || K.P != null) {
            kalman_record_state(K);
        }
        if (filter_is_varying(K)) {
            /* we have time-varying coefficients */
            err = kalman_refresh_matrices(K, prn);
            if (err) {
                K.loglik = NADBL;
                break;
            }
        }
        if (kdebug > 1) {
            kalman_print_state(K);
        }
        /* calculate v_t, checking for missing values */
        nt = compute_forecast_error(K);
        if (nt == 0) {
            /* skip this observation */
            K.okN -= 1;
            handle_missing_obs(K);
	    if (kdebug) {
		fprintf(stderr, "dejong filter: NA at t=%d\n", K.t);
	    }
            continue;
        }
        /* and record forecast errors if wanted */
        if (K.V != null) {
            record_to_row(K.V, K.vt, K.t);
        }
        /* calculate F_t = ZPZ' [+ VY] (where VY = GG') */
        if (K.VY != null) {
            fast_copy_values(K.Ft, K.VY);
            gretl_matrix_qform(K.ZT, GretlMatrixMod.transpose, K.P0,
                               K.Ft, GretlMatrixMod.cumulate);
        } else {
            gretl_matrix_qform(K.ZT, GretlMatrixMod.transpose, K.P0,
                               K.Ft, GretlMatrixMod.none);
        }
        /* calculate M_t = TPZ' [+ HG'] */
        gretl_matrix_multiply(K.P0, K.ZT, K.PZ);
        if (K.TI) {
            fast_copy_values(K.Mt, K.PZ);
        } else {
            gretl_matrix_multiply(K.T, K.PZ, K.Mt);
        }
        if (K.HG != null) {
            gretl_matrix_add_to(K.Mt, K.HG);
        }
        /* calculate F_t and its inverse */
        fast_copy_values(K.iFt, K.Ft);
        if (K.t < K.d && !loglik_1991) {
	    /* conditionality regarding the log-determinant? */
	    err = gretl_invert_symmetric_matrix2(K.iFt, &pldt);
            sum_preldet += pldt;
	} else {
	    err = gretl_invert_symmetric_matrix2(K.iFt, &ldt);
	}
        if (err) {
            fprintf(stderr, "kfilter_dejong: failed to invert Ft\n");
	    break;
        }

        if (K.F != null) {
            /* we're recording F_t for all t */
            if (K.smo_prep) {
                /* record inverse */
                record_to_vech(K.Finv, K.iFt, nt, K.t);
            } else {
                /* record F_t itself */
                record_to_vech(K.F, K.Ft, nt, K.t);
            }
        }

        /* gain: K_t = M_t F_t^{-1} */
        gretl_matrix_multiply(K.Mt, K.iFt, K.Kt);
        /* L_t = T_t - K_t Z_t */
        fast_copy_values(K.Lt, K.T);
        gretl_matrix_multiply_mod(K.Kt, GretlMatrixMod.none,
                                  K.ZT, GretlMatrixMod.transpose,
                                  K.Lt, GretlMatrixMod.decrement);

        /* J_t = H_t - K_t * G_t */
        fast_copy_values(K.Jt, K.H);
        if (K.G != null) {
            gretl_matrix_multiply_mod(K.Kt, GretlMatrixMod.none,
                                      K.G,  GretlMatrixMod.none,
                                      K.Jt, GretlMatrixMod.decrement);
        }

        /* record matrices for the smoother */
        if (K.L != null) {
            record_to_col(K.L, K.Lt, K.t);
        }
        if (K.J != null) {
            record_to_col(K.J, K.Jt, K.t);
        }
        if (K.K != null) {
            record_to_vec(K.K, K.Kt, K.t);
        }

        /* update state: a1 = T a0 + K_t v_t [+ mu] */
        gretl_matrix_multiply(K.T, K.a0, K.a1);
        if (K.mu != null) {
            gretl_matrix_add_to(K.a1, K.mu);
        }
        gretl_matrix_multiply_mod(K.Kt, GretlMatrixMod.none,
                                  K.vt, GretlMatrixMod.none,
                                  K.a1, GretlMatrixMod.cumulate);
        fast_copy_values(K.a0, K.a1);

        /* update var(state): P1 = TPL' + HJ' */
        gretl_matrix_multiply_mod(K.H, GretlMatrixMod.none,
                                  K.Jt, GretlMatrixMod.transpose,
                                  K.P1, GretlMatrixMod.none);
        gretl_matrix_multiply_mod(K.P0, GretlMatrixMod.none,
                                  K.Lt, GretlMatrixMod.transpose,
                                  K.rr, GretlMatrixMod.none);
        if (K.TI) {
            gretl_matrix_add_to(K.P1, K.rr);
        } else {
            gretl_matrix_multiply_mod(K.T,  GretlMatrixMod.none,
                                      K.rr, GretlMatrixMod.none,
                                      K.P1, GretlMatrixMod.cumulate);
        }
        fast_copy_values(K.P0, K.P1);

        if (exact_diffuse && K.t < K.d) {
	    /* still in diffuse phase */
            err = dejong_diffuse_filter_step(K, &qt, &ldt);
        } else {
            /* standard Kalman phase */
            qt = gretl_scalar_qform(K.vt, K.iFt, &err);
            llct = nl2pi;
        }

        /* determine and record loglikelihood */
        if (err) {
            K.loglik = NADBL;
            break;
        } else {
            llt = -0.5 * (llct + ldt + qt);
            if (na(llt)) {
                K.loglik = NADBL;
                break;
            }
            K.loglik += llt;
            K.qsum += qt;
            sumldet += ldt;
        }
        if (K.LL != null) {
            gretl_vector_set(K.LL, K.t, llt);
        }
    }

    set_kalman_stopped(K);

    if (na(K.loglik)) {
        err = GretlError.nan;
    } else if (exact_diffuse) {
        if (K.d == K.N) {
            err = handle_no_collapse_case(K, nl2pi, sum_preldet);
        } else if (loglik_1991) {
            /* Deal with extra terms set out in De Jong (1991)?
               Not activated at present.
            */
            double ll_adj = K.djinfo.qm + K.djinfo.ldS;

            //fprintf(stderr, "HERE qm = %g, ldS = %g, ll_adj = %g\n",
            //        K.djinfo.qm, K.djinfo.ldS, ll_adj);
            K.loglik -= 0.5 * ll_adj;
        } else {
            /* FIXME: the following needs more testing */
            K.s2 = K.qsum / (K.n * K.okN - K.r);
        }
    } else if (kalman_arma_ll(K)) {
        double ll1 = 1.0 + LN_2_PI + log(K.qsum / K.okN);

        K.loglik = -0.5 * (K.okN * ll1 + sumldet);
    } else {
        int d = (kalman_diffuse(K) && !K.exact)? K.r : 0;

        /* as per kfilter_standard() code: do we want this? */
        K.s2 = K.qsum / (K.n * K.okN - d);
    }

    if (kdebug) {
        fprintf(stderr, "kfilter_dejong: err=%d, ll=%#0.8g, d=%d\n",
                err, K.loglik, K.d);
    }

    return err;
}

/* r_{t-1} = Z_t' F_t^{-1} v_t + L_t' r_t */

private void rt_recursion(Kalman * K, GretlMatrixRaw * n1, GretlMatrixRaw * r0, GretlMatrixRaw * r1)
{
    gretl_matrix_multiply(K.iFt, K.vt, n1);
    if (K.t == K.N - 1) {
        /* simplification for the last observation */
        gretl_matrix_multiply(K.ZT, n1, r0);
    } else {
        gretl_matrix_multiply(K.ZT, n1, r1);
        gretl_matrix_multiply_mod(K.Lt, GretlMatrixMod.transpose,
                                  r0, GretlMatrixMod.none,
                                  r1, GretlMatrixMod.cumulate);
        fast_copy_values(r0, r1);
    }
}

/* N_{t-1} = Z_t' F_t^{-1} Z_t + L_t' N_t L_t */

private void Nt_recursion(Kalman * K, GretlMatrixRaw * N0, GretlMatrixRaw * N1)
{
    if (K.t == K.N - 1) {
        /* simplification for the last observation */
        gretl_matrix_qform(K.ZT, GretlMatrixMod.none,
                           K.iFt, N0, GretlMatrixMod.none);
     } else {
        gretl_matrix_qform(K.ZT, GretlMatrixMod.none,
                           K.iFt, N1, GretlMatrixMod.none);
        load_from_col(K.Lt, K.L, K.t);
        gretl_matrix_qform(K.Lt, GretlMatrixMod.transpose,
                           N0, N1, GretlMatrixMod.cumulate);
        fast_copy_values(N0, N1);
    }
}

private int record_dhat_and_Psi(Kalman * K, GretlMatrixRaw * dhat, GretlMatrixRaw * Psi, GretlMatrixRaw * rt, GretlMatrixRaw * Nt)
{
    DejongInfo * dj = K.djinfo;
    GretlMatrixRaw * tmp;
    int err = 0;

    /* dhat = S * (s + Am' * rt) */
    tmp = gretl_matrix_copy(dj.s);
    gretl_matrix_multiply_mod(dj.Am, GretlMatrixMod.transpose,
                              rt, GretlMatrixMod.none,
                              tmp, GretlMatrixMod.cumulate);
    gretl_matrix_multiply(dj.S, tmp, dhat);
    gretl_matrix_free(tmp);

    /* Psi = S - S * Am' * Nt * Am * S */
    fast_copy_values(Psi, dj.S);
    gretl_matrix_qform(dj.Am, GretlMatrixMod.transpose,
                       Nt, K.rr, GretlMatrixMod.none);
    gretl_matrix_qform(dj.S, GretlMatrixMod.none, K.rr,
                       Psi, GretlMatrixMod.decrement);

    return err;
}

private int dejong_diffuse_smoother_step(Kalman * K, GretlMatrixRaw * rt, GretlMatrixRaw * Rt, GretlMatrixRaw * Bt, GretlMatrixRaw * Nt, GretlMatrixRaw * dhat, GretlMatrixRaw * Psi, int no_collapse)
{
    DejongInfo * dj = K.djinfo;
    GretlMatrixRaw * tmp1 = null;
    GretlMatrixRaw * tmp2 = null;
    GretlMatrixRaw * Rr = null;
    int err = 0;

    load_from_col(dj.Vt, dj.V, K.t);
    load_from_col(dj.At, dj.A, K.t);
    load_from_col(K.Lt, K.L, K.t);

    tmp1 = gretl_matrix_alloc(K.r, K.r + 1);
    tmp2 = gretl_matrix_alloc(K.n, K.r + 1);

    /* Rr = Rt ~ rt */
    Rr = gretl_matrix_alloc(K.r, K.r + 1);
    gretl_matrix_inscribe_matrix(Rr, Rt, 0, 0, GretlMatrixMod.none);
    gretl_matrix_inscribe_matrix(Rr, rt, 0, K.r, GretlMatrixMod.none);

    /* Vv = Vt ~ vt */
    gretl_matrix_inscribe_matrix(dj.Vv, dj.Vt, 0, 0, GretlMatrixMod.none);
    gretl_matrix_inscribe_matrix(dj.Vv, K.vt, 0, K.r, GretlMatrixMod.none);

    /* Rr = Z'*iFt*Vv + Lt'*Rr */
    gretl_matrix_multiply_mod(K.Lt, GretlMatrixMod.transpose,
                              Rr, GretlMatrixMod.none,
                              tmp1, GretlMatrixMod.none);
    fast_copy_values(Rr, tmp1);
    gretl_matrix_multiply(K.iFt, dj.Vv, tmp2);
    gretl_matrix_multiply_mod(K.ZT, GretlMatrixMod.none,
                              tmp2, GretlMatrixMod.none,
                              Rr, GretlMatrixMod.cumulate);

    /* Rt = Rr[,1:r] ; rt = Rr[,r+1] */
    gretl_matrix_extract_matrix(Rt, Rr, 0, 0, GretlMatrixMod.none);
    gretl_matrix_extract_matrix(rt, Rr, 0, K.r, GretlMatrixMod.none);

    /* Bt = At + Pt * Rt */
    fast_copy_values(Bt, dj.At);
    gretl_matrix_multiply_mod(K.P0, GretlMatrixMod.none,
                              Rt, GretlMatrixMod.none,
                              Bt, GretlMatrixMod.cumulate);

    /* aht = at + Pt * rt */
    gretl_matrix_multiply_mod(K.P0, GretlMatrixMod.none,
                              rt, GretlMatrixMod.none,
                              K.a0, GretlMatrixMod.cumulate);

    /* aht += Bt * dhat */
    gretl_matrix_multiply_mod(Bt, GretlMatrixMod.none,
                              dhat, GretlMatrixMod.none,
                              K.a0, GretlMatrixMod.cumulate);

    /* vht = Pt - Pt * Nt * Pt */
    fast_copy_values(K.P1, K.P0);
    gretl_matrix_qform(K.P0, GretlMatrixMod.transpose,
                       Nt, K.P1, GretlMatrixMod.decrement);

    if (no_collapse) {
        /* vht += Bt * Psi * Bt' */
        gretl_matrix_qform(Bt, GretlMatrixMod.none,
                           Psi, K.P1, GretlMatrixMod.cumulate);
    } else {
        GretlMatrixRaw * BSC = gretl_matrix_reuse(tmp1, K.r, K.r);

        /* Ct = Pt * (Rt + Nt * At) */
        gretl_matrix_multiply(Nt, dj.At, K.rr);
        gretl_matrix_add_to(K.rr, Rt);
        gretl_matrix_multiply(K.P0, K.rr, K.Ct);

        /* BSC = Bt * S * Ct' */
        gretl_matrix_multiply_mod(dj.S, GretlMatrixMod.none,
                                  K.Ct, GretlMatrixMod.transpose,
                                  K.rr, GretlMatrixMod.none);
        gretl_matrix_multiply(Bt, K.rr, BSC);

        /* vht += Bt * Psi * Bt' - BSC - BSC' */
        gretl_matrix_qform(Bt, GretlMatrixMod.none,
                           Psi, K.rr, GretlMatrixMod.none);
        gretl_matrix_subtract_from(K.rr, BSC);
        gretl_square_matrix_transpose(BSC);
        gretl_matrix_subtract_from(K.rr, BSC);
        gretl_matrix_add_to(K.P1, K.rr);
    }

    record_to_row(K.A, K.a0, K.t);
    record_to_vech(K.P, K.P1, K.r, K.t);

    gretl_matrix_free(tmp1);
    gretl_matrix_free(tmp2);
    gretl_matrix_free(Rr); /* FIXME make reusable? */

    return err;
}

private int state_smooth_dejong(Kalman * K)
{
    GretlMatrixBlock * B;
    GretlMatrixBlock * B2 = null;
    GretlMatrixRaw * r0;
    GretlMatrixRaw * r1;
    GretlMatrixRaw * N0;
    GretlMatrixRaw * N1;
    GretlMatrixRaw * n1;
    GretlMatrixRaw * Bt;
    GretlMatrixRaw * Rt;
    GretlMatrixRaw * dhat;
    GretlMatrixRaw * Psi;
    int no_collapse = 0;
    int nt = K.n;
    int t;
    int err = 0;

    no_collapse = (K.d == K.N);

    if (trace) {
        printf("state_smooth_dejong()\n");
    }
    if (kdebug) {
        fprintf(stderr, "state_smooth_dejong: no_collapse = %d\n", no_collapse);
    }

    B = gretl_matrix_block_new(&r0,  K.r, 1,
                               &r1,  K.r, 1,
                               &N0,  K.r, K.r,
                               &N1,  K.r, K.r,
                               &n1,  K.n, 1,
                               null);
    if (B == null) {
        err = GretlError.alloc;
    }

    if (!err && K.exact) {
        B2 = gretl_matrix_block_new(&Bt,   K.r, K.r,
                                    &Rt,   K.r, K.r,
                                    &dhat, K.r, 1,
                                    &Psi,  K.r, K.r,
                                    null);
        if (B2 == null) {
            err = GretlError.alloc;
        }
    }

    if (err) {
        return err;
    }

    gretl_matrix_zero(r0);
    gretl_matrix_zero(N0);
    if (B2 != null) {
        gretl_matrix_zero(Rt);
    }

    if (no_collapse) {
        gretl_matrix_multiply(K.djinfo.S, K.djinfo.s, dhat);
        gretl_matrix_copy_values(Psi, K.djinfo.S);
    }

    for (t=K.N-1; t>=0 && !err; t--) {
        K.t = t;
	nt = load_filter_data(K, no_collapse, &err);
        if (err) {
            break;
        } else if (kdebug && nt == 0) {
            fprintf(stderr, "state_smooth_dejong: nt=0 at t=%d\n", K.t);
        }

        /* compute N_{t-1} */
	Nt_recursion(K, N0, N1);

        if (t >= K.d) {
            /* regular recursion: compute r_{t-1} */
	    rt_recursion(K, n1, r0, r1);
            /* a_{t|T} = a_{t|t-1} + P_{t|t-1} r_{t-1} */
            fast_copy_values(K.a1, K.a0);
            gretl_matrix_multiply_mod(K.P0, GretlMatrixMod.none,
                                      r0, GretlMatrixMod.none,
                                      K.a1, GretlMatrixMod.cumulate);
            record_to_row(K.A, K.a1, t);
            /* P_{t|T} = P_{t|t-1} - P_{t|t-1} N_{t-1} P_{t|t-1} */
            fast_copy_values(K.P1, K.P0);
            gretl_matrix_qform(K.P0, GretlMatrixMod.none, N0,
                               K.P1, GretlMatrixMod.decrement);
            record_to_vech(K.P, K.P1, K.r, t);
            if (t == K.d) {
                record_dhat_and_Psi(K, dhat, Psi, r1, N1);
            }
        } else {
            /* t < K.d: diffuse phase (FIXME nt == 0?) */
            dejong_diffuse_smoother_step(K, r1, Rt, Bt, N1,
                                         dhat, Psi, no_collapse);
        }
    }

    gretl_matrix_block_destroy(B2);
    gretl_matrix_block_destroy(B);

    return err;
}

private int diffuse_dist_smooth_step(Kalman * K, GretlMatrixRaw * Rt, GretlMatrixRaw * Bt, GretlMatrixRaw * Ct, GretlMatrixRaw * Nt, GretlMatrixRaw * dhat, GretlMatrixRaw * Psi, GretlMatrixRaw * nr, GretlMatrixRaw * nu_t, GretlMatrixRaw * Vnu_t, int DKstyle, int no_collapse)
{
    DejongInfo * dj = K.djinfo;
    GretlMatrixRaw * BSC = null;
    GretlMatrixRaw * pp = null;
    GretlMatrixRaw * rp = null;
    int err = 0;

    if (!no_collapse) {
        BSC = gretl_matrix_alloc(K.p, K.p);
        rp = gretl_matrix_alloc(K.r, K.p);
    }
    pp = gretl_matrix_alloc(K.p, K.p);

    /* Bt = G' * iFt * Vt + Jt' * Rt */
    gretl_matrix_multiply_mod(K.Jt, GretlMatrixMod.transpose,
                              Rt, GretlMatrixMod.none,
                              Bt, GretlMatrixMod.none);
    if (K.G != null) {
        gretl_matrix_multiply(K.iFt, dj.Vt, nr);
        gretl_matrix_multiply_mod(K.G, GretlMatrixMod.transpose,
                                  nr, GretlMatrixMod.none,
                                  Bt, GretlMatrixMod.cumulate);
    }

    if (!no_collapse) {
        /* Ct = Jt' * (Rt + Nt * A_{t+1}) */
        gretl_matrix_multiply(Nt, dj.At, K.rr);
        gretl_matrix_add_to(K.rr, Rt);
        gretl_matrix_multiply_mod(K.Jt, GretlMatrixMod.transpose,
                                  K.rr, GretlMatrixMod.none,
                                  Ct, GretlMatrixMod.none);
    }

    /* nu_t += Bt * dhat */
    gretl_matrix_multiply_mod(Bt, GretlMatrixMod.none,
                              dhat, GretlMatrixMod.none,
                              nu_t, GretlMatrixMod.cumulate);

    if (!no_collapse) {
        /* BSC = Bt * S * Ct' */
        gretl_matrix_multiply_mod(dj.S, GretlMatrixMod.none,
                                  Ct, GretlMatrixMod.transpose,
                                  rp, GretlMatrixMod.none);
        gretl_matrix_multiply(Bt, rp, BSC);
    }

    /* Vnu_adj = Bt * Psi * Bt' - BSC - BSC' */
    gretl_matrix_qform(Bt, GretlMatrixMod.none, Psi,
                       pp, GretlMatrixMod.none);
    if (!no_collapse) {
        gretl_matrix_subtract_from(pp, BSC);
        gretl_square_matrix_transpose(BSC);
        gretl_matrix_subtract_from(pp, BSC);
    }

    if (DKstyle) {
        gretl_matrix_add_to(Vnu_t, pp);
    } else {
        gretl_matrix_subtract_from(Vnu_t, pp);
    }

    if (!no_collapse) {
        gretl_matrix_free(BSC);
        gretl_matrix_free(rp);
    }
    gretl_matrix_free(pp);

    return err;
}

private int decompose_nu(Kalman * K, GretlMatrixRaw * etahat, GretlMatrixRaw * veta, GretlMatrixRaw * epshat, GretlMatrixRaw * veps, GretlMatrixRaw * nu, GretlMatrixRaw * vnu, GretlMatrixRaw * H_nu, GretlMatrixRaw * G_nu, GretlMatrixRaw * nu_t, GretlMatrixRaw * vnu_t)
{
    GretlMatrixRaw * vrr;
    GretlMatrixRaw * vnn = null;
    int t;

    vrr = gretl_matrix_alloc(K.r, K.r);
    if (K.G != null) {
        vnn = gretl_matrix_alloc(K.n, K.n);
    }

    for (t=0; t<K.N; t++) {
        /* retrieve t-dated @nu and @vnu */
        load_from_col(nu_t, nu, t);
        load_from_col(vnu_t, vnu, t);
        /* work on etahat, veta */
        gretl_matrix_multiply(K.H, nu_t, H_nu);
        record_to_row(etahat, H_nu, t);
        gretl_matrix_qform(K.H, GretlMatrixMod.none,
                           vnu_t, vrr, GretlMatrixMod.none);
        record_to_vec(veta, vrr, t);
        /* if applicable, work on epshat, veps */
        if (K.G != null) {
            gretl_matrix_multiply(K.G, nu_t, G_nu);
            record_to_row(epshat, G_nu, t);
            gretl_matrix_qform(K.G, GretlMatrixMod.none,
                               vnu_t, vnn, GretlMatrixMod.none);
            record_to_vec(veps, vnn, t);
        }
    }

    gretl_matrix_free(vrr);
    gretl_matrix_free(vnn);

    return 0;
}

private int dist_smooth_dejong(Kalman * K, int DKstyle)
{
    GretlMatrixBlock * B;
    GretlMatrixBlock * B2 = null;
    GretlMatrixRaw * nu_t;
    GretlMatrixRaw * vnu_t;
    GretlMatrixRaw * r0;
    GretlMatrixRaw * r1;
    GretlMatrixRaw * N0;
    GretlMatrixRaw * N1;
    GretlMatrixRaw * n1;
    GretlMatrixRaw * Bt;
    GretlMatrixRaw * Ct;
    GretlMatrixRaw * Rt;
    GretlMatrixRaw * nr;
    GretlMatrixRaw * dhat;
    GretlMatrixRaw * Psi;
    GretlMatrixRaw * nu = null;
    GretlMatrixRaw * vnu = null;
    GretlMatrixRaw * H_nu;
    GretlMatrixRaw * G_nu = null;
    GretlMatrixRaw * etahat;
    GretlMatrixRaw * veta;
    GretlMatrixRaw * epshat = null;
    GretlMatrixRaw * veps = null;
    DejongInfo * dj = K.djinfo;
    int no_collapse;
    int nt = K.n;
    int t;
    int err = 0;

    no_collapse = (K.d == K.N);

    if (trace) {
        fprintf(stderr, "dist_smooth_dejong(), DKstyle %d, K.d = %d\n",
                DKstyle, K.d);
        gretl_matrix_print(K.G, "K.G");
        gretl_matrix_print(K.H, "K.H");
    }

    /* the "global" disturbances */
    nu = gretl_zero_matrix_new(K.p, K.N);
    vnu = gretl_zero_matrix_new(K.p * K.p, K.N);

    etahat = gretl_matrix_alloc(K.N, K.r);
    veta = gretl_matrix_alloc(K.N, K.r * K.r);
    if (K.G != null) {
        epshat = gretl_matrix_alloc(K.N, K.n);
        veps = gretl_matrix_alloc(K.N, K.n * K.n);
        G_nu = gretl_matrix_alloc(K.n, 1);
    }

    B = gretl_matrix_block_new(&r0,  K.r, 1,
                               &r1,  K.r, 1,
                               &N0,  K.r, K.r,
                               &N1,  K.r, K.r,
                               &n1,  K.n, 1,
                               &nr,  K.n, K.r,
                               &nu_t,  K.p, 1,
                               &vnu_t, K.p, K.p,
                               &H_nu,  K.r, 1,
                               null);
    if (B == null) {
        err = GretlError.alloc;
    }

    if (!err && K.exact) {
        B2 = gretl_matrix_block_new(&Bt,   K.p, K.r,
                                    &Ct,   K.p, K.r,
                                    &Rt,   K.r, K.r,
                                    &dhat, K.r, 1,
                                    &Psi,  K.r, K.r,
                                    null);
        if (B2 == null) {
            err = GretlError.alloc;
        }
    }

    if (err) {
        return err;
    }

    gretl_matrix_zero(r0);
    gretl_matrix_zero(N0);
    if (B2 != null) {
        gretl_matrix_zero(Rt);
    }

    for (t=K.N-1; t>=0 && !err; t--) {
        K.t = t;
        nt = load_filter_data(K, 0, &err);
        if (err) {
            break;
        } else if (kdebug && nt == 0) {
	    fprintf(stderr, "dist_smooth_dejong: nt=0 at t=%d\n", K.t);
        }

	load_from_col(K.Jt, K.J, K.t);
	load_from_col(K.Lt, K.L, K.t);

	/* nu_t = G' * iFt * vt + Jt' * rt */
	gretl_matrix_multiply_mod(K.Jt, GretlMatrixMod.transpose,
				  r0, GretlMatrixMod.none,
				  nu_t, GretlMatrixMod.none);
	if (K.G != null) {
	    gretl_matrix_multiply(K.iFt, K.vt, n1);
	    gretl_matrix_multiply_mod(K.G, GretlMatrixMod.transpose,
				      n1, GretlMatrixMod.none,
				      nu_t, GretlMatrixMod.cumulate);
	}

        if (DKstyle) {
            /* vnu = I(p) - (G' * iFt * G + Jt' * Nt * Jt) */
            gretl_matrix_inscribe_I(vnu_t, 0, 0, K.p);
            gretl_matrix_qform(K.Jt, GretlMatrixMod.transpose,
                               N0, vnu_t, GretlMatrixMod.decrement); /* N1? */
            if (K.G != null) {
                gretl_matrix_qform(K.G, GretlMatrixMod.transpose,
                                   K.iFt, vnu_t, GretlMatrixMod.decrement);
            }
        } else {
	    /* vnu = G' * iFt * G + Jt' * Nt * Jt */
	    gretl_matrix_qform(K.Jt, GretlMatrixMod.transpose, N0,
			       vnu_t, GretlMatrixMod.none); /* N1? */
	    if (K.G != null) {
		gretl_matrix_qform(K.G, GretlMatrixMod.transpose, K.iFt,
				   vnu_t, GretlMatrixMod.cumulate);
	    }
        }

        if (t == K.d) {
            /* rm1 = Z' * iFt * vt + Lt' * rt */
            rt_recursion(K, n1, r0, r1);
            /* Nm1 = Z' * iFt * Z + Lt' * Nt * Lt */
            Nt_recursion(K, N0, N1);
            record_dhat_and_Psi(K, dhat, Psi, r1, N1); /* r0? N0? */
        } else if (t < K.d) {
            load_from_col(dj.Vt, dj.V, t);
            if (t == K.d - 1) {
                fast_copy_values(dj.At, dj.Am);
            } else {
                load_from_col(dj.At, dj.A, t);
            }
            diffuse_dist_smooth_step(K, Rt, Bt, Ct, N1, dhat, Psi,
                                     nr, nu_t, vnu_t, DKstyle,
                                     no_collapse);
        }

        /* record t-dated results */
        record_to_col(nu, nu_t, t);
        record_to_col(vnu, vnu_t, t);

        if (t > 1) {
            /* t-1 dated quantities */
            if (t == K.d) {
                fast_copy_values(r0, r1);
                fast_copy_values(N0, N1);
            } else {
                /* rt = Z' * iFt * vt + Lt' * rt */
                rt_recursion(K, n1, r0, r1);
                /* Nt = Z' * iFt * Z + Lt' * Nt * Lt */
                Nt_recursion(K, N0, N1);
            }
            if (t < K.d) {
		/* Rt = Z' * iFt * Vt + Lt' * Rt */
		gretl_matrix_multiply_mod(K.Lt, GretlMatrixMod.transpose,
					  Rt, GretlMatrixMod.none,
					  K.rr, GretlMatrixMod.none);
		gretl_matrix_multiply(K.iFt, dj.Vt, nr);
		gretl_matrix_multiply_mod(K.ZT, GretlMatrixMod.none,
					  nr, GretlMatrixMod.none,
					  K.rr, GretlMatrixMod.cumulate);
                fast_copy_values(Rt, K.rr);
            }
        }
    }

    decompose_nu(K, etahat, veta, epshat, veps, nu, vnu,
                 H_nu, G_nu, nu_t, vnu_t);

    /* add results to K.b */
    bundle_add_matrix(K.b, "etahat", etahat);
    bundle_add_matrix(K.b, "veta", veta);
    if (epshat != null) {
        bundle_add_matrix(K.b, "epshat", epshat);
        bundle_add_matrix(K.b, "veps",  veps);
    }

    gretl_matrix_free(nu);
    gretl_matrix_free(vnu);
    gretl_matrix_free(G_nu);
    gretl_matrix_block_destroy(B);
    gretl_matrix_block_destroy(B2);

    return err;
}

/* ---- the contents of kalman_legacy.c (included at the end of kalman.c under SUPPORT_LEGACY) ---- */


/* for use in Durbin-Koopman disturbance smoothing */
struct DkDistInfo {
    GretlMatrixRaw * etahat;
    GretlMatrixRaw * veta;
    GretlMatrixRaw * epshat;
    GretlMatrixRaw * veps;
    GretlMatrixRaw * D;
    GretlMatrixRaw * vut;
    GretlMatrixRaw * vwt;
};

private void clear_dist_info(DkDistInfo * di, int err)
{
    /* workspace matrices: free unconditionally */
    gretl_matrix_free(di.vut);
    gretl_matrix_free(di.vwt);
    gretl_matrix_free(di.D);
    if (err) {
        /* On success (err == 0), these matrices are added to
           the Kalman bundle, and so must not be freed here.
        */
        gretl_matrix_free(di.etahat);
        gretl_matrix_free(di.veta);
        gretl_matrix_free(di.epshat);
        gretl_matrix_free(di.veps);
    }
}

/* For disturbance smoothing: ensure we have on hand matrices that are
   correctly sized to hold estimates of the variance of the
   disturbance(s) in the state and (if applicable) observation
   equations.
*/

private int allocate_dist_info(Kalman * K, DkDistInfo * di)
{
    int n = K.VY == null ? 0 : K.n;
    int err = 0;

    /* null everything in case of error */
    di.vut = di.vwt = null;
    di.veta = di.veps = null;
    di.etahat = di.epshat = null;
    di.D = null;

    /* FIXME what about the DK_VAR case? */

    /* recorder matrices for all time steps */
    di.etahat = gretl_matrix_alloc(K.N, K.r);
    di.veta = gretl_matrix_alloc(K.N, K.r * K.r);
    if (di.etahat == null || di.veta == null) {
        err = GretlError.alloc;
    }
    if (!err && n > 0) {
        di.epshat = gretl_matrix_alloc(K.N, n);
        di.veps = gretl_matrix_alloc(K.N, n);
        di.D = gretl_matrix_alloc(n, n);
        if (di.epshat == null ||
            di.veps == null ||
            di.D == null) {
            err = GretlError.alloc;
        }
    }

    /* step-t square matrices */
    if (!err) {
        di.vwt = gretl_matrix_alloc(K.r, K.r);
        if (di.vwt == null) {
            err = GretlError.alloc;
        }
    }
    if (!err && n > 0) {
        di.vut = gretl_matrix_alloc(n, n);
        if (di.vut == null) {
            err = GretlError.alloc;
        }
    }

    if (err) {
        /* destroy everything */
        clear_dist_info(di, 1);
    }

    return err;
}

/* Write the diagonal elements of the square matrix @src
   into row @t of @targ.
*/

private void record_to_diag(GretlMatrixRaw * targ, const(GretlMatrixRaw) * src, int t)
{
    int i;
    int n = src.rows;
    double xti;

    for (i=0; i<n; i++) {
        xti = gretl_matrix_get(src, i, i);
        gretl_matrix_set(targ, t, i, xti);
    }
}

private int ksmooth_refresh_matrices(Kalman * K, Prn * prn)
{
    GretlMatrixRaw ** [2] mptr = [&K.VS, &K.VY];
    int[2] idx = [KMat.vs, KMat.vy];
    int cross_update = 0;
    int dkvar_update = 0;
    int i;
    int j;
    int err = 0;

    if (kalman_djvar(K)) {
        mptr[0] = &K.H;
        mptr[1] = &K.G;
    }

    if (K.matcall != null) {
        err = kalman_update_matrices(K, prn);
    }

    for (i=0; i<2 && !err; i++) {
        j = idx[i];
        if (matrix_is_varying(K, j)) {
            if (kalman_djvar(K) && (j == KMat.vs || j == KMat.vy)) {
                /* handle revised H and/or G */
                cross_update = 1;
            } else if (kalman_dkvar(K) && j == KMat.vs) {
                /* handle revised Q and/or R */
                dkvar_update = 1;
            } else {
                err = check_matrix_dims(K, *mptr[i], j);
            }
            if (err) {
                fprintf(stderr, "ksmooth_refresh_matrices: err = %d at t = %d\n",
                        err, K.t);
            }
        }
    }

    if (!err) {
        if (cross_update) {
            /* cross-correlated case */
            err = kalman_update_crossinfo(K, KUpdate.step);
        } else if (dkvar_update) {
            /* Durbin-Koopman case */
            err = kalman_update_dkvar(K, KUpdate.step);
        }
    }

    return err;
}

enum {
    SM_TYPE_NONE, /* not doing smoothing */
    SM_STATE_STD, /* regular state smoother */
    SM_DIST_BKWD, /* disturbance smoother, backward pass */
    SM_DIST_FRWD, /* disturbance smoother, forward pass */
};

/* For use with smoothing: load what we need for step @t, from the
   record that was kept on the prior forecasting pass. What we need
   from the forward pass depends on what exactly we're computing,
   which is conveyed by @smtype. We must handle time-variation in
   any matrices for which that is relevant to the smoother.
*/

private int legacy_load_filter_data(Kalman * K, int smtype)
{
    int err = 0;

    /* load the forecast error */
    load_from_row(K.vt, K.V, K.t);

    if (smtype == SM_DIST_BKWD) {
        /* disturbances */
        if (filter_is_varying(K)) {
            ksmooth_refresh_matrices(K, null);
        }
    } else {
        /* state: get T_t and/or Z_t if need be */
        if (matrix_is_varying(K, KMat.t)) {
            err = retrieve_Tt(K);
        }
        if (!err && matrix_is_varying(K, KMat.zt)) {
            err = retrieve_Zt(K);
        }
        if (err) {
            return err;
        }
    }

    if (smtype < SM_DIST_BKWD) {
        /* load the state and its MSE */
        load_from_row(K.a0, K.A, K.t);
        load_from_vech(K.P0, K.P, K.r, K.t);
    }

    if (smtype == SM_STATE_STD || smtype == SM_DIST_BKWD) {
        /* load the gain and F^{-1} */
        load_from_vec(K.Kt, K.K, K.t);
        load_from_vech(K.iFt, K.Finv, K.n, K.t);
    }

    return err;
}

/* Calculate the variance of the smoothed disturbances for the
   cross-correlated case. See Koopman, Shephard and Doornik (1998),
   page 19, var(\varepsilon_t|Y_n).
*/

private int legacy_combined_dist_variance(Kalman * K, DkDistInfo * di, GretlMatrixRaw * Nt, GretlMatrixBlock * BX, int DKstyle)
{
    GretlMatrixRaw * DG;
    GretlMatrixRaw * KN;
    GretlMatrixRaw * VP;
    GretlMatrixRaw * NH;
    GretlMatrixRaw * NK;

    DG  = gretl_matrix_block_get_matrix(BX, 0);
    KN  = gretl_matrix_block_get_matrix(BX, 1);
    VP  = gretl_matrix_block_get_matrix(BX, 2);
    NH  = gretl_matrix_block_get_matrix(BX, 3);

    /* First chunk in Koopman's notation:
       G_t'(D_t*G_t - K_t'*N_t*H_t)
    */
    KN = gretl_matrix_reuse(KN, K.n, K.r);
    gretl_matrix_multiply(di.D, K.G, DG);
    gretl_matrix_multiply_mod(K.Kt, GretlMatrixMod.transpose,
                              Nt, GretlMatrixMod.none,
                              KN, GretlMatrixMod.none);
    gretl_matrix_multiply_mod(KN, GretlMatrixMod.none,
                              K.H, GretlMatrixMod.none,
                              DG, GretlMatrixMod.decrement);
    gretl_matrix_multiply_mod(K.G, GretlMatrixMod.transpose,
                              DG, GretlMatrixMod.none,
                              VP, GretlMatrixMod.none);

    /* Second chunk, to be added to the above
       H_t'(N_t*H_t - N_t*K_t*G_t)
    */
    NK = gretl_matrix_reuse(KN, K.r, K.n);
    gretl_matrix_multiply(Nt, K.H, NH);
    gretl_matrix_multiply_mod(Nt, GretlMatrixMod.transpose,
                              K.Kt, GretlMatrixMod.none,
                              NK, GretlMatrixMod.none);
    gretl_matrix_multiply_mod(NK, GretlMatrixMod.none,
                              K.G, GretlMatrixMod.none,
                              NH, GretlMatrixMod.decrement);
    gretl_matrix_multiply_mod(K.H, GretlMatrixMod.transpose,
                              NH, GretlMatrixMod.none,
                              VP, GretlMatrixMod.cumulate);

    if (DKstyle) {
        /* VP = I_p - VP */
        double vii;
        int i;

        gretl_matrix_multiply_by_scalar(VP, -1.0);
        for (i=0; i<K.p; i++) {
            vii = gretl_matrix_get(VP, i, i);
            gretl_matrix_set(VP, i, i, 1.0 + vii);
        }
    }

    /* VP (p x p) holds the variance of \epsilon_t
       conditional on Y_n: now form the per-equation
       disturbance variance matrices, @vwt and @vut,
       for this time-step.
    */
    gretl_matrix_qform(K.H, GretlMatrixMod.none, VP,
                       di.vwt, GretlMatrixMod.none);
    gretl_matrix_qform(K.G, GretlMatrixMod.none, VP,
                       di.vut, GretlMatrixMod.none);

    return 0;
}

private int legacy_dist_variance(Kalman * K, DkDistInfo * di, GretlMatrixRaw * Nt, GretlMatrixBlock * BX, int DKstyle)
{
    int err = 0;

    if (di.D != null) {
	/* needed only in presence of obs disturbance */
	/* D_t = F_t^{-1} + K_t' N_t K_t */
	fast_copy_values(di.D, K.iFt);
	if (K.t < K.N - 1) {
	    gretl_matrix_qform(K.Kt, GretlMatrixMod.transpose,
			       Nt, di.D, GretlMatrixMod.cumulate);
	}
    }

    if (K.p > 0) {
        /* cross-correlated disturbance variance */
        err = legacy_combined_dist_variance(K, di, Nt, BX, DKstyle);
        if (!err) {
            record_to_diag(di.veta, di.vwt, K.t);
            record_to_diag(di.veps, di.vut, K.t);
        }
    } else {
	/* variance of state disturbance */
	if (DKstyle) {
	    /* HH' - HH' N_t HH' */
	    fast_copy_values(di.vwt, K.VS);
	    gretl_matrix_qform(K.VS, GretlMatrixMod.transpose,
			       Nt, di.vwt, GretlMatrixMod.decrement);
	} else {
	    /* HH' N_t HH' */
	    gretl_matrix_qform(K.VS, GretlMatrixMod.transpose,
			       Nt, di.vwt, GretlMatrixMod.none);
	}
	record_to_diag(di.veta, di.vwt, K.t);

        /* variance of obs disturbance, if applicable */
	if (K.VY != null) {
	    if (DKstyle) {
		/* GG' - GG D_t GG' */
		fast_copy_values(di.vut, K.VY);
		gretl_matrix_qform(K.VY, GretlMatrixMod.transpose,
				   di.D, di.vut, GretlMatrixMod.decrement);
	    } else {
		/* GG' D_t GG' */
		gretl_matrix_qform(K.VY, GretlMatrixMod.transpose,
				   di.D, di.vut, GretlMatrixMod.none);
	    }
	    record_to_diag(di.veps, di.vut, K.t);
        }
    }

    return err;
}

/* This iteration is in common between the state smoother
   (Anderson-Moore) and the disturbance smoother (Koopman).
*/

private void LrN_iteration(Kalman * K, GretlMatrixRaw * L, GretlMatrixRaw * n1, GretlMatrixRaw * r0, GretlMatrixRaw * r1, GretlMatrixRaw * N0, GretlMatrixRaw * N1)
{
    if (K.t < K.N - 1) {
	/* L_t = T_t - K_t Z_t */
	fast_copy_values(L, K.T);
	gretl_matrix_multiply_mod(K.Kt, GretlMatrixMod.none,
				  K.ZT, GretlMatrixMod.transpose,
				  L, GretlMatrixMod.decrement);
    }

    /* r_{t-1} = Z_t' F_t^{-1} v_t + L_t' r_t */
    gretl_matrix_multiply(K.iFt, K.vt, n1);
    if (K.t == K.N - 1) {
	gretl_matrix_multiply(K.ZT, n1, r0);
    } else {
	gretl_matrix_multiply(K.ZT, n1, r1);
	gretl_matrix_multiply_mod(L, GretlMatrixMod.transpose,
				  r0, GretlMatrixMod.none,
				  r1, GretlMatrixMod.cumulate);
	fast_copy_values(r0, r1);
    }

    /* N_{t-1} = Z_t' F_t^{-1} Z_t + L_t' N_t L_t */
    if (K.t == K.N - 1) {
	gretl_matrix_qform(K.ZT, GretlMatrixMod.none,
			   K.iFt, N0, GretlMatrixMod.none);
    } else {
	gretl_matrix_qform(K.ZT, GretlMatrixMod.none,
			   K.iFt, N1, GretlMatrixMod.none);
	gretl_matrix_qform(L, GretlMatrixMod.transpose,
			   N0, N1, GretlMatrixMod.cumulate);
	fast_copy_values(N0, N1);
    }
}

/* Initial smoothed state: a + P*r0 */

private void koopman_calc_a0(Kalman * K, GretlMatrixRaw * r0)
{
    if (K.Pini != null) {
        gretl_matrix_multiply(K.Pini, r0, K.a0);
    } else {
        set_initial_statevar(K);
        gretl_matrix_multiply(K.P0, r0, K.a0);
    }
    if (K.aini != null) {
        gretl_matrix_add_to(K.a0, K.aini);
    }
    record_to_row(K.A, K.a0, 0);
}

/* Disturbance smoothing -- see Koopman, Shephard and Doornik
   (SsfPack doc), section 4.4.
*/

private int legacy_koopman_smooth(Kalman * K, int DKstyle)
{
    DkDistInfo di;
    GretlMatrixBlock * B;
    GretlMatrixBlock * BX = null;
    GretlMatrixRaw * u;
    GretlMatrixRaw * L;
    GretlMatrixRaw * R;
    GretlMatrixRaw * r0;
    GretlMatrixRaw * r1;
    GretlMatrixRaw * N0;
    GretlMatrixRaw * N1;
    GretlMatrixRaw * n1;
    GretlMatrixRaw * tr;
    GretlMatrixRaw * DG = null;
    GretlMatrixRaw * KN = null;
    GretlMatrixRaw * RZS = null;
    GretlMatrixRaw * NH = null;
    GretlMatrixRaw * Ut = null;
    int ft_min = 0;
    int t;
    int err = 0;

    if (trace) {
        printf("koopman_smooth(), legacy, DKstyle = %d\n", DKstyle);
    }

    B = gretl_matrix_block_new(&u,  K.n, 1,
                               &L,  K.r, K.r,
                               &R,  K.N, K.r,
                               &r0, K.r, 1,
                               &r1, K.r, 1,
                               &N0, K.r, K.r,
                               &N1, K.r, K.r,
                               &n1, K.n, 1,
                               &tr, K.r, 1,
                               null);

    if (B == null) {
        return GretlError.alloc;
    }

    /* for variance of smoothed disturbances */
    err = allocate_dist_info(K, &di);

    if (K.p > 0) {
	/* cross-correlated disturbances */
        BX = gretl_matrix_block_new(&DG,  K.n, K.p,
                                    &KN,  K.n, K.r,
                                    &RZS, K.p, K.p,
                                    &NH,  K.r, K.p,
                                    &Ut,  K.p, 1,
                                    null);
        if (BX == null) {
            err = GretlError.alloc;
        }
    }

    if (err) {
        goto bailout;
    }

    gretl_matrix_zero(r0);
    gretl_matrix_zero(N0);

    /* The backward recursion */

    for (t=K.N-1; t>=0 && !err; t--) {
	K.t = t;
	err = legacy_load_filter_data(K, SM_DIST_BKWD);
        if (err) {
            break;
        }

        /* u_t = F_t^{-1} v_t - K_t' r_t */
        gretl_matrix_multiply(K.iFt, K.vt, u);
        if (t < K.N - 1) {
            gretl_matrix_multiply_mod(K.Kt, GretlMatrixMod.transpose,
                                      r0, GretlMatrixMod.none,
                                      u, GretlMatrixMod.decrement);
        }
        /* Store u_t values in K.V: these are needed in
           the forward pass to compute the smoothed
           disturbances. Also save r_t in R.
        */
        record_to_row(K.V, u, t);
	record_to_row(R, r0, t);

	/* compute variance of disturbances */
	err = legacy_dist_variance(K, &di, N0, BX, DKstyle);
	if (err) {
	    break;
	}

	/* compute r_{t-1}, N_{t-1} */
	LrN_iteration(K, L, n1, r0, r1, N0, N1);

	if (t == 0) {
	    /* compute initial smoothed state */
	    koopman_calc_a0(K, r0);
        }
    }

    /* Forward iteration for smoothed disturbances, all time steps,
       plus smoothed state from t = 1 onward.
    */
    for (t=ft_min; t<K.N; t++) {
	K.t = t;
        err = legacy_load_filter_data(K, SM_DIST_FRWD);
        if (err) {
            break;
        }

	/* state disturbance */
        load_from_row(r0, R, t);
        if (K.p > 0) {
            gretl_matrix_multiply_mod(K.H, GretlMatrixMod.transpose,
                                      r0, GretlMatrixMod.none,
                                      Ut, GretlMatrixMod.none);
            load_from_row(K.vt, K.V, t);
            gretl_matrix_multiply_mod(K.G, GretlMatrixMod.transpose,
                                      K.vt, GretlMatrixMod.none,
                                      Ut, GretlMatrixMod.cumulate);
            gretl_matrix_multiply(K.H, Ut, r1);
        } else {
            gretl_matrix_multiply(K.VS, r0, r1);
        }
        record_to_row(R, r1, t);
        record_to_row(di.etahat, r1, t);

	/* observation disturbance */
	if (K.VY != null) {
	    if (K.p > 0) {
		gretl_matrix_multiply(K.G, Ut, n1);
	    } else {
		gretl_matrix_multiply(K.VY, K.vt, n1);
	    }
            record_to_row(di.epshat, n1, t);
	}

	/* state: a_{t+1} = T a_t + w_t (or + H*eps_t) */
        if (t > 0) {
            load_from_row(K.a0, K.A, t-1);
            gretl_matrix_multiply(K.T, K.a0, K.a1);
            vector_from_row_mod(K.a1, R, t-1, GretlMatrixMod.cumulate);
            if (K.mu != null) {
                gretl_matrix_add_to(K.a1, K.mu);
            }
            record_to_row(K.A, K.a1, t);
        }
    }

    if (!err) {
        if (di.etahat != null) {
            bundle_add_matrix(K.b, "etahat",  di.etahat);
            bundle_add_matrix(K.b, "veta", di.veta);
        }
        if (di.epshat != null) {
            bundle_add_matrix(K.b, "epshat",  di.epshat);
            bundle_add_matrix(K.b, "veps", di.veps);
        }
    }

 bailout:

    gretl_matrix_block_destroy(B);
    gretl_matrix_block_destroy(BX);
    clear_dist_info(&di, err);

    return err;
}

/* Anderson-Moore state smoothing */

private int anderson_moore_smooth(Kalman * K)
{
    GretlMatrixBlock * B;
    GretlMatrixRaw * r0;
    GretlMatrixRaw * r1;
    GretlMatrixRaw * N0;
    GretlMatrixRaw * N1;
    GretlMatrixRaw * n1;
    GretlMatrixRaw * L;
    int t;
    int err = 0;

    if (trace) {
        printf("anderson_moore_smooth()\n");
        if (K.K == null || K.Finv == null) {
            printf(" K.K %p, K.Finv %p\n", cast(void *) K.K,
                   cast(void *) K.Finv);
        } else {
            printf(" K.K is %d x %d, K.Finv is %d x %d\n",
                   K.K.rows, K.K.cols, K.Finv.rows,
                   K.Finv.cols);
        }
    }

    B = gretl_matrix_block_new(&r0,  K.r, 1,
                               &r1,  K.r, 1,
                               &N0,  K.r, K.r,
                               &N1,  K.r, K.r,
                               &n1,  K.n, 1,
                               &L,   K.r, K.r,
                               null);
    if (B == null) {
        return GretlError.alloc;
    }

    gretl_matrix_zero(r0);
    gretl_matrix_zero(N0);

    for (t=K.N-1; t>=0 && !err; t--) {
	K.t = t;
        err = legacy_load_filter_data(K, SM_STATE_STD);
        if (err) {
            break;
        }

	/* compute r_{t-1}, N_{t-1} */
	LrN_iteration(K, L, n1, r0, r1, N0, N1);

        /* a_{t|T} = a_{t|t-1} + P_{t|t-1} r_{t-1} */
        fast_copy_values(K.a1, K.a0);
        gretl_matrix_multiply_mod(K.P0, GretlMatrixMod.none,
                                  r0, GretlMatrixMod.none,
                                  K.a1, GretlMatrixMod.cumulate);
        record_to_row(K.A, K.a1, t);

        /* P_{t|T} = P_{t|t-1} - P_{t|t-1} N_{t-1} P_{t|t-1} */
        fast_copy_values(K.P1, K.P0);
        gretl_matrix_qform(K.P0, GretlMatrixMod.none, N0,
                           K.P1, GretlMatrixMod.decrement);
        record_to_vech(K.P, K.P1, K.r, t);
    }

    gretl_matrix_block_destroy(B);

    return err;
}
