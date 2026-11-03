/* Port of gretl_commands.h */

module gretl.gretl_commandsh;


enum GretlCmdIndex {
    add = 1,
    adf,
    anova,
    append,
    ar,
    ar1,
    arch,
    arma,
    bds,
    biprobit,
    bkw,
    break_,
    bxplot,
    chow,
    clear,
    coeffsum,
    coint,
    coint2,
    continue_,
    corr,
    corrgm,
    cusum,
    data,
    datamod,
    deleet,
    diff,
    difftest,
    discrete,
    dpanel,
    dummify,
    duration,
    elif,
    else_,
    end,
    endif,
    endloop,
    eqnprint,
    equation,
    estimate,
    eval,
    fcast,
    flush,
    foreign,
    fractint,
    freq,
    func,
    funcerr,
    garch,
    genr,
    gibbs,
    gmm,
    gnuplot,
    gpbuild,
    graphpg,
    gridplot,
    heckit,
    help,
    hfplot,
    hsk,
    hurst,
    if_,
    include,
    info,
    intreg,
    join,
    kdplot,
    kpss,
    labels,
    lad,
    lags,
    ldiff,
    leverage,
    levinlin,
    logistic,
    logit,
    logs,
    loop,
    mahal,
    makepkg,
    markers,
    meantest,
    midasreg,
    mle,
    modeltab,
    modprint,
    modtest,
    mpi,
    mpols,
    negbin,
    nls,
    normtest,
    nulldata,
    ols,
    omit,
    open,
    orthdev,
    outfile,
    panel,
    panplot,
    panspec,
    pca,
    pergm,
    plot,
    poisson,
    print,
    printf,
    probit,
    pval,
    quantreg,
    qlrtest,
    qqplot,
    quit,
    rename,
    reset,
    restrict,
    rmplot,
    run,
    runs,
    scatters,
    sdiff,
    set,
    setinfo,
    setobs,
    setopt,
    setmiss,
    shell,
    smpl,
    spearman,
    square,
    stdize,
    store,
    summary,
    system,
    tabprint,
    textplot,
    tobit,
    ivreg,
    var,
    varlist,
    vartest,
    vecm,
    vif,
    wls,
    xcorrgm,
    xtab,
    funcret,
    catch_,
    pkg,
    nc
}

bool TEXTSAVE_OK(GretlCmdIndex c) {
    return (c == GretlCmdIndex.adf ||
            c == GretlCmdIndex.bds ||
            c == GretlCmdIndex.bkw ||
            c == GretlCmdIndex.chow ||
            c == GretlCmdIndex.coeffsum ||
            c == GretlCmdIndex.coint ||
            c == GretlCmdIndex.coint2 ||
            c == GretlCmdIndex.corr ||
            c == GretlCmdIndex.corrgm ||
            c == GretlCmdIndex.cusum ||
            c == GretlCmdIndex.fcast ||
            c == GretlCmdIndex.freq ||
            c == GretlCmdIndex.panspec ||
            c == GretlCmdIndex.kpss ||
            c == GretlCmdIndex.leverage ||
            c == GretlCmdIndex.levinlin ||
            c == GretlCmdIndex.modtest ||
            c == GretlCmdIndex.meantest ||
            c == GretlCmdIndex.mpols ||
            c == GretlCmdIndex.normtest ||
            c == GretlCmdIndex.pca ||
            c == GretlCmdIndex.pergm ||
            c == GretlCmdIndex.print ||
            c == GretlCmdIndex.printf ||
            c == GretlCmdIndex.pval ||
            c == GretlCmdIndex.qlrtest ||
            c == GretlCmdIndex.reset ||
            c == GretlCmdIndex.rmplot ||
            c == GretlCmdIndex.runs ||
            c == GretlCmdIndex.spearman ||
            c == GretlCmdIndex.summary ||
            c == GretlCmdIndex.textplot ||
            c == GretlCmdIndex.vif ||
            c == GretlCmdIndex.vartest ||
            c == GretlCmdIndex.xcorrgm ||
            c == GretlCmdIndex.xtab);
}

bool NEEDS_MODEL_CHECK(GretlCmdIndex c) {
    return (c == GretlCmdIndex.add ||
            c == GretlCmdIndex.bkw ||
            c == GretlCmdIndex.chow ||
            c == GretlCmdIndex.coeffsum ||
            c == GretlCmdIndex.cusum ||
            c == GretlCmdIndex.eqnprint ||
            c == GretlCmdIndex.fcast ||
            c == GretlCmdIndex.panspec ||
            c == GretlCmdIndex.leverage ||
            c == GretlCmdIndex.modtest ||
            c == GretlCmdIndex.omit ||
            c == GretlCmdIndex.qlrtest ||
            c == GretlCmdIndex.reset ||
            c == GretlCmdIndex.tabprint ||
            c == GretlCmdIndex.vif);
}

bool GRAPHING_COMMAND(GretlCmdIndex c) {
    return (c == GretlCmdIndex.gnuplot ||
            c == GretlCmdIndex.plot ||
            c == GretlCmdIndex.bxplot ||
            c == GretlCmdIndex.scatters ||
            c == GretlCmdIndex.corrgm ||
            c == GretlCmdIndex.xcorrgm ||
            c == GretlCmdIndex.pergm ||
            c == GretlCmdIndex.rmplot ||
            c == GretlCmdIndex.hfplot ||
            c == GretlCmdIndex.hurst ||
            c == GretlCmdIndex.leverage ||
            c == GretlCmdIndex.qlrtest ||
            c == GretlCmdIndex.qqplot ||
            c == GretlCmdIndex.kdplot ||
            c == GretlCmdIndex.freq ||
            c == GretlCmdIndex.fcast ||
            c == GretlCmdIndex.corr ||
            c == GretlCmdIndex.cusum ||
            c == GretlCmdIndex.panplot ||
            c == GretlCmdIndex.gridplot);
}
