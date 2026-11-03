/* Port of libset.h */

module gretl.libseth;


import gretl.libset : libset_set_double;

enum VECMnorm {
    phillips,
    diag,
    first,
    none,
    max
}

enum OptimCode {
    auto_,
    bfgs,
    newton,
    max
}

enum HACopts {
    refuse,
    es,
    am
}

enum PanelRobust {
    arellano,
    beck_katz,
    driscoll_kraay
}

enum SetKey {
    use_cwd = 1 << 0,  /* store: use current dir as default */
    echo_on = 1 << 1,  /* echoing commands or not */
    msgs_on = 1 << 2,  /* emitting non-error messages or not */
    force_decpoint = 1 << 3,  /* override locale decimal character */
    use_svd = 1 << 4,  /* SVD decomposition is matrix OLS default */
    prewhiten = 1 << 5,  /* HAC pre-whitening? */
    force_hc = 1 << 6,  /* don't use HAC for time series */
    use_lbfgs = 1 << 7,  /* prefer LBFGS to BFGS? */
    shell_ok = 1 << 8,  /* "shell" facility is approved? */
    warnings = 1 << 9,  /* print numerical warning messages */
    skip_missing = 1 << 10,  /* skip NAs when building matrix from series */
    bfgs_rstep = 1 << 11,  /* use Richardson in BFGS numerical gradient */
    robust_z = 1 << 12,  /* use z- not t-score with HCCM/HAC */
    mwrite_g = 1 << 13,  /* use %g format with mwrite() */
    mpi_use_smt = 1 << 14,  /* MPI: use hyperthreads by default */
    state_flag_max = 1 << 15,  /* separator */
    /* state small int (but non-boolean) vars */
    gretl_optim,
    vecm_norm,
    garch_vcv,
    garch_alt_vcv,
    arma_vcv,
    wildboot_dist,
    fdjac_qual,
    use_qr,
    max_verbose,
    hc_version,
    panel_robust,
    hac_kernel,
    hac_lag,
    user_hac_lag,
    lbfgs_mem,
    quantile_type,
    state_small_int_max,  /* separator: start state int vars */
    horizon,
    bootrep,
    loop_maxiter,
    bfgs_maxiter,
    bfgs_verbskip,
    boot_iters,
    bhhh_maxiter,
    rq_maxiter,
    gmm_maxiter,
    state_int_max,  /* separator: end state int vars */
    conv_huge,
    nls_toler,
    bfgs_toler,
    bfgs_maxgrad,
    bhhh_toler,
    qs_bandwidth,
    nadarwat_trim,
    state_float_max,  /* separator: end state floats */
    csv_write_na,
    csv_read_na,
    initvals,
    initcurv,
    matmask,
    state_vars_max,  /* separator */
    /* non-state vars follow */
    gretl_debug,
    gretl_assert,
    datacols,
    plot_collect,
    r_functions,
    r_lib,
    loglevel,
    logstamp,
    csv_digits,
    hac_missvals,
    ns_small_int_max,  /* separator */
    gmp_bits,
    ns_max,  /* separator */
    blas_mnk_min,
    omp_n_threads,
    simd_k_max,
    simd_mn_min,
    ns_int_max,  /* separator */
    seed,
    csv_delim,
    stopwatch,
    verbose,
    sv_workdir,
    sv_logfile,
    graph_theme,
    disp_digits,
    tex_plot_opts,
    setvar_max  /* sentinel */
}

alias ShowActivityFunc = void function();

int set_nls_toler(double x)
{
    return libset_set_double(SetKey.nls_toler, x);
}
