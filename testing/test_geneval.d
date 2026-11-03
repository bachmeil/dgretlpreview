/* Smoke tests for geneval1.d - geneval5.d (the "genr" syntax-tree evaluator). The lexer and parser (genlex.c,
   gensyntax.c) are not ported, so there is no way yet to evaluate a formula given as text: instead the syntax trees
   are built by hand, as the parser would build them, and handed to eval() or to the individual node functions.

   Build and run with testing/run_tests.sh. */

import core.stdc.math : sqrt;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;
import core.stdc.stdlib : calloc, free;
import core.stdc.string : strcmp;

import gretl.dataset : create_new_dataset, destroy_dataset;
import gretl.geneval0;
import gretl.geneval1 : get_aux_node, free_tree, xy_calc;
import gretl.geneval2 : escape_strsplit_sep, variant_string_func;
import gretl.geneval3 : get_bundle_member;
import gretl.geneval4 : dvar_get_scalar, eval_epochday, node_type_error;
import gretl.geneval5 : eval, get_op, get_opstr;
import gretl.genmainh : RetrievalIndex;
import gretl.genparseh;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_list : gretl_list_new;
import gretl.gretl_matrix : gretl_matrix_alloc, gretl_matrix_free;
import gretl.gretl_matrixh : GretlMatrixRaw;
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

/* ---- helpers to build syntax trees by hand ---- */

Node * mk(int t)
{
    Node * n = cast(Node *) calloc(1, Node.sizeof);

    n.t = cast(short) t;
    return n;
}

Node * num(double x)
{
    Node * n = mk(GenSym.num);

    n.v.xval = x;
    return n;
}

Node * str(string s)
{
    Node * n = mk(GenSym.str);

    n.v.str = cast(char *) (s ~ "\0").ptr;
    return n;
}

Node * series(DatasetRaw * d, int v)
{
    Node * n = mk(GenSym.series);

    n.vnum = v;
    n.v.xvec = d.Z[v];
    return n;
}

Node * mat(int r, int c, double[] vals...)
{
    Node * n = mk(GenSym.mat);

    n.v.m = gretl_matrix_alloc(r, c);
    n.v.m.val[0 .. vals.length] = vals[];
    return n;
}

Node * bin(int op, Node * l, Node * r)
{
    Node * n = mk(op);

    n.L = l;
    n.R = r;
    return n;
}

Node * fn1(int f, Node * l)
{
    Node * n = mk(f);

    n.L = l;
    n.R = mk(GenSym.empty);
    return n;
}

Parser newparser(DatasetRaw * d)
{
    Parser p;

    p.dset = d;
    p.flags = GenFlags.start;
    p.targ = GenSym.unk;
    return p;
}

/* evaluate @tree, run @test on the result, then release the tree */

void run(Node * tree, Parser * p, void delegate(Node *) test)
{
    p.err = 0;
    p.tree = tree;
    Node * r = eval(tree, p);

    test(r);
    free_tree(tree, p, FreeCode.tree);
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    DatasetRaw * d = create_new_dataset(3, 6, 0);

    d.pd = 1;
    d.structure = TsCodes.cross_section;
    foreach (i, nm; ["const", "x", "y"]) {
        d.varname[i][0 .. nm.length] = nm[];
        d.varname[i][nm.length] = 0;
    }
    foreach (t; 0 .. 6) {
        d.Z[0][t] = 1.0;
        d.Z[1][t] = t + 1.0;     /* 1 .. 6 */
        d.Z[2][t] = 2.0 * t;     /* 0, 2 .. 10 */
    }
    Parser p = newparser(d);

    /* scalar arithmetic and comparisons */
    run(bin(GenSym.b_add, num(2), num(3.5)), &p, (Node * r) {
        check(p.err == 0 && r.t == GenSym.num && r.v.xval == 5.5, "scalar 2 + 3.5");
    });
    run(bin(GenSym.b_mul, num(4), bin(GenSym.b_sub, num(10), num(7))), &p, (Node * r) {
        check(r.t == GenSym.num && r.v.xval == 12, "scalar 4 * (10 - 7), nested tree");
    });
    run(bin(GenSym.b_pow, num(2), num(10)), &p, (Node * r) {
        check(r.v.xval == 1024, "scalar 2 ^ 10");
    });
    run(bin(GenSym.b_gt, num(3), num(2)), &p, (Node * r) {
        check(r.v.xval == 1, "scalar 3 > 2");
    });
    run(bin(GenSym.b_and, num(0), num(1)), &p, (Node * r) {
        check(r.v.xval == 0, "logical AND short-circuits on a false left side");
    });
    run(bin(GenSym.b_or, num(1), num(0)), &p, (Node * r) {
        check(r.v.xval == 1, "logical OR short-circuits on a true left side");
    });
    run(bin(GenSym.b_add, num(NADBL), num(1)), &p, (Node * r) {
        check(na(r.v.xval), "NA propagates through scalar addition");
    });

    /* the query operator */
    {
        Node * q = mk(GenSym.query);

        q.L = num(1);
        q.M = num(10);
        q.R = num(20);
        run(q, &p, (Node * r) {
            check(r.t == GenSym.num && r.v.xval == 10, "1 ? 10 : 20");
        });
        q = mk(GenSym.query);
        q.L = num(0);
        q.M = num(10);
        q.R = num(20);
        run(q, &p, (Node * r) {
            check(r.v.xval == 20, "0 ? 10 : 20");
        });
    }

    /* series */
    run(bin(GenSym.b_add, series(d, 1), num(10)), &p, (Node * r) {
        check(r.t == GenSym.series && r.v.xvec[0] == 11 && r.v.xvec[5] == 16, "series + scalar");
    });
    run(bin(GenSym.b_mul, series(d, 1), series(d, 2)), &p, (Node * r) {
        check(r.t == GenSym.series && r.v.xvec[3] == 4 * 6 && r.v.xvec[0] == 0, "series * series");
    });
    run(bin(GenSym.b_gt, series(d, 1), num(3)), &p, (Node * r) {
        check(r.t == GenSym.series && r.v.xvec[2] == 0 && r.v.xvec[3] == 1, "series > scalar is a 0/1 series");
    });

    /* functions of one argument, applied to a scalar, a series and a matrix */
    {
        Node * f = fn1(GenFunc.f_sqrt, num(16));

        f.v.ptr = cast(void *) &sqrt;
        run(f, &p, (Node * r) {
            check(r.t == GenSym.num && r.v.xval == 4, "sqrt of a scalar");
        });
        f = fn1(GenFunc.f_sqrt, series(d, 1));
        f.v.ptr = cast(void *) &sqrt;
        run(f, &p, (Node * r) {
            check(r.t == GenSym.series && r.v.xvec[3] == 2 && r.v.xvec[0] == 1, "sqrt of a series");
        });
        f = fn1(GenFunc.f_sqrt, mat(1, 3, 1, 4, 9));
        f.v.ptr = cast(void *) &sqrt;
        run(f, &p, (Node * r) {
            check(r.t == GenSym.mat && r.v.m.val[2] == 3, "sqrt of a matrix");
        });
    }
    run(fn1(GenSym.u_neg, num(5)), &p, (Node * r) {
        check(r.v.xval == -5, "unary minus");
    });

    /* series -> scalar functions */
    run(fn1(GenFunc.f_sum, series(d, 1)), &p, (Node * r) {
        check(r.t == GenSym.num && r.v.xval == 21, "sum(x)");
    });
    run(fn1(GenFunc.f_mean, series(d, 1)), &p, (Node * r) {
        check(r.v.xval == 3.5, "mean(x)");
    });
    run(fn1(GenFunc.f_max, series(d, 2)), &p, (Node * r) {
        check(r.v.xval == 10, "max(y) of a single argument");
    });
    run(bin(GenFunc.f_max, num(3), num(9)), &p, (Node * r) {
        check(r.v.xval == 9, "max(3, 9) of two arguments");
    });

    /* matrices */
    run(bin(GenSym.b_mul, mat(2, 2, 1, 2, 3, 4), mat(2, 2, 5, 6, 7, 8)), &p, (Node * r) {
        /* column-major storage: A = [1 3; 2 4], B = [5 7; 6 8] */
        check(r.t == GenSym.mat && r.v.m.val[0] == 23 && r.v.m.val[1] == 34 && r.v.m.val[2] == 31 && r.v.m.val[3] == 46,
              "matrix * matrix");
    });
    run(bin(GenSym.b_add, mat(1, 3, 1, 2, 3), num(10)), &p, (Node * r) {
        check(r.t == GenSym.mat && r.v.m.val[0] == 11 && r.v.m.val[2] == 13, "matrix + scalar");
    });
    run(bin(GenSym.b_dotmult, mat(1, 3, 1, 2, 3), mat(1, 3, 4, 5, 6)), &p, (Node * r) {
        check(r.t == GenSym.mat && r.v.m.val[1] == 10, "matrix .* matrix");
    });
    run(bin(GenSym.b_add, mat(2, 2, 1, 2, 3, 4), mat(1, 3, 1, 2, 3)), &p, (Node * r) {
        check(p.err == GretlError.nonconf, "nonconformable matrix addition is an error");
    });
    run(fn1(GenFunc.f_rows, mat(3, 2, 1, 2, 3, 4, 5, 6)), &p, (Node * r) {
        check(r.t == GenSym.num && r.v.xval == 3, "rows()");
    });

    /* strings */
    run(bin(GenSym.b_eq, str("abc"), str("abc")), &p, (Node * r) {
        check(r.t == GenSym.num && r.v.xval == 1, "string == string");
    });
    run(bin(GenSym.b_neq, str("abc"), str("abd")), &p, (Node * r) {
        check(r.v.xval == 1, "string != string");
    });
    {
        Node * l = str("Hello World");
        Node * f = mk(GenFunc.f_tolower);

        f.L = l;
        Node * r = variant_string_func(l, null, null, GenFunc.f_tolower, &p);

        check(r !is null && strcmp(r.v.str, "hello world") == 0, "tolower() on a string");
        p.err = 0;
    }
    {
        char * s = escape_strsplit_sep("a\\tb\\n");

        check(s[0] == 'a' && s[1] == '\t' && s[2] == 'b' && s[3] == '\n' && s[4] == 0, "escape_strsplit_sep");
        free(s);
    }


    /* more series functions: lag, sort, correlation, quantile */
    run(bin(GenSym.lag, series(d, 1), num(-1)), &p, (Node * r) {
        check(r.t == GenSym.series && na(r.v.xvec[0]) && r.v.xvec[3] == 3, "x(-1)");
    });
    {
        Node * f = fn1(GenFunc.f_dsort, mat(1, 4, 3, 1, 4, 2));

        run(f, &p, (Node * r) {
            check(r.t == GenSym.mat && r.v.m.val[0] == 4 && r.v.m.val[3] == 1, "dsort() of a vector");
        });
    }
    run(bin(GenFunc.f_cor, series(d, 1), series(d, 2)), &p, (Node * r) {
        check(r.t == GenSym.num && r.v.xval > 0.9999 && r.v.xval < 1.0001, "corr(x, y)");
    });
    run(bin(GenFunc.f_quantile, series(d, 1), num(0.5)), &p, (Node * r) {
        check(r.t == GenSym.num && r.v.xval == 3.5, "quantile(x, 0.5)");
    });

    /* matrix functions */
    run(fn1(GenFunc.f_transp, mat(2, 3, 1, 2, 3, 4, 5, 6)), &p, (Node * r) {
        check(r.t == GenSym.mat && r.v.m.rows == 3 && r.v.m.cols == 2 && r.v.m.val[1] == 3, "transp()");
    });
    run(fn1(GenFunc.f_inv, mat(2, 2, 2, 0, 0, 4)), &p, (Node * r) {
        check(r.t == GenSym.mat && r.v.m.val[0] == 0.5 && r.v.m.val[3] == 0.25, "inv() of a diagonal matrix");
    });
    run(bin(GenFunc.f_zeros, num(2), num(3)), &p, (Node * r) {
        check(r.t == GenSym.mat && r.v.m.rows == 2 && r.v.m.cols == 3 && r.v.m.val[5] == 0, "zeros(2, 3)");
    });
    run(fn1(GenFunc.f_cols, mat(2, 3, 1, 2, 3, 4, 5, 6)), &p, (Node * r) {
        check(r.t == GenSym.num && r.v.xval == 3, "cols()");
    });

    /* three-argument functions */
    {
        Node * f = mk(GenFunc.f_seq);

        f.L = num(1);
        f.M = num(9);
        f.R = num(4);
        run(f, &p, (Node * r) {
            check(r.t == GenSym.mat && r.v.m.val[0] == 1 && r.v.m.val[1] == 5 && r.v.m.val[2] == 9, "seq(1, 9, 4)");
        });
    }

    /* more string functions */
    {
        Node * f = mk(GenFunc.f_instring);

        f.L = str("hello world");
        f.M = str("o w");
        f.R = mk(GenSym.empty);
        run(f, &p, (Node * r) {
            check(r.t == GenSym.num && r.v.xval == 1, "instring()");
        });
    }
    {
        Node * f = mk(GenFunc.f_strlen);

        f.L = str("h\u00e9llo");
        f.R = mk(GenSym.empty);
        run(f, &p, (Node * r) {
            check(r.t == GenSym.num && r.v.xval == 5, "strlen() counts characters, not bytes");
        });
    }

    /* the type-error reporting */
    {
        Node * bad = str("x");

        p.err = 0;
        node_type_error(GenFunc.f_sqrt, 1, GenSym.num, bad, &p);
        check(p.err == GretlError.types, "node_type_error sets E_TYPES");
        p.err = 0;
    }

    /* the assignment-operator parser */
    {
        char[3] a = "+=";
        char[3] b = "=x";
        char[3] c = "++";

        check(get_op(a.ptr) == GenSym.b_add && get_op(b.ptr) == GenSym.b_asn && get_op(c.ptr) == GenSym.inc, "get_op");
        check(strcmp(get_opstr(GenSym.b_mul), "*=") == 0 && strcmp(get_opstr(GenSym.b_asn), "=") == 0, "get_opstr");
    }

    /* xy_calc on plain numbers, and the dataset-variable accessors */
    check(xy_calc(6, 3, GenSym.b_div, GenSym.num, &p) == 2 && xy_calc(6, 3, GenSym.b_sub, GenSym.num, &p) == 3, "xy_calc");
    check(dvar_get_scalar(RetrievalIndex.nobs, d) == 6 && dvar_get_scalar(RetrievalIndex.nvars, d) == 3
          && dvar_get_scalar(RetrievalIndex.pd, d) == 1, "dvar_get_scalar: $nobs, $nvars, $pd");
    check(na(dvar_get_scalar(RetrievalIndex.nobs, null)), "dvar_get_scalar with no dataset is NA");

    /* epoch days */
    {
        Node * y = num(1970);
        Node * m = num(1);
        Node * dd = num(2);
        Node * r = eval_epochday(y, m, dd, &p);

        check(p.err == 0 && r !is null && r.t == GenSym.num && r.v.xval == epoch_day_1970_01_02(), "epochday(1970,1,2)");
        p.err = 0;
        y.v.xval = 2024;
        m.v.xval = 13;
        eval_epochday(y, m, dd, &p);
        check(p.err == GretlError.invarg, "epochday rejects month 13");
        p.err = 0;
    }

    destroy_dataset(d);
    printf("\n%d checks, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}

/* 1970-01-01 is day 719163 of the proleptic Gregorian count used by g_date, so 1970-01-02 is 719164 */
double epoch_day_1970_01_02() { return 719164; }
