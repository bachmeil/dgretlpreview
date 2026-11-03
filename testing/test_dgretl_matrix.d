/* Tests for dgretl.matrix: copy by default (m[a .. b, c .. d]), explicit reference views (sub, rows, cols), strides,
   views of views, lifetime of the data, bounds errors.

   Build and run with testing/run_tests.sh. */

import core.stdc.stdio : printf;

import std.exception : collectException;
import std.conv : to;
import std.math : abs, cos, exp, fabs, isNaN, log, sin, sqrt;

import dgretl;

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

GretlMatrix make(int r, int c)
{
    auto m = GretlMatrix(gretl_matrix_alloc(r, c));
    foreach (j; 0 .. c)
        foreach (i; 0 .. r)
            m[i, j] = 10.0 * i + j;
    return m;
}

void clobber(GretlMatrix x)
{
    x[0, 0] = -1;
}

/* a view that outlives every other handle on the matrix */
Submatrix lonelyView()
{
    auto m = make(4, 4);
    return m.sub[1 .. 3, 1 .. 3];
}

bool close(double a, double b)
{
    return a - b < 1e-9 && b - a < 1e-9;
}

int main()
{
    auto m = make(5, 4);
    check(m.nrow == 5 && m.ncol == 4, "nrow and ncol");

    /* copy by default */
    auto c = m[1 .. 3, 2 .. 4];
    check(c.nrow == 2 && c.ncol == 2 && c[0, 0] == 12 && c[1, 1] == 23, "m[a..b, c..d] has the right size and elements");
    clobber(c);
    c[1, 1] = 99;
    check(m[1, 2] == 12 && m[2, 3] == 23, "changing the copy leaves the original alone");
    auto all = m[0 .. $, 0 .. $];
    check(all.nrow == 5 && all.ncol == 4 && all[4, 3] == 43, "$ in a slice");
    check(m.dup[3, 1] == 31, "dup");

    /* explicit references */
    auto s = m.sub[1 .. 4, 1 .. 3];
    check(s.nrow == 3 && s.ncol == 2 && s.stride == 5 && s[0, 0] == 11 && s[2, 1] == 32, "sub: size, stride, elements");
    check(!s.isContiguous, "a block is not contiguous");
    s[0, 1] = 777;
    check(m[1, 2] == 777, "writing through a sub changes the matrix");
    m[3, 1] = -5;
    check(s[2, 0] == -5, "a sub sees later changes of the matrix");

    auto r = m.rows[2 .. 4];
    check(r.nrow == 2 && r.ncol == 4 && r[0, 3] == 23 && r[1, 0] == 30, "rows");
    r[1, 0] = 300;
    check(m[3, 0] == 300, "writing through rows");

    auto cc = m.cols[1 .. 3];
    check(cc.nrow == 5 && cc.ncol == 2 && cc.isContiguous && cc[4, 1] == 42, "cols are contiguous");
    check(m.cols[2 .. $].ncol == 2, "$ in cols");

    /* views of views, through UFCS functions of both types */
    auto ss = s.sub[1 .. 3, 0 .. 1];
    check(ss.nrow == 2 && ss.ncol == 1 && ss[0, 0] == 21 && ss[1, 0] == -5, "sub of a sub");
    auto rr = r.cols[1 .. 3];
    check(rr.nrow == 2 && rr[0, 0] == 21 && rr[1, 1] == 32, "cols of rows");
    check(s.rows[1 .. 2][0, 0] == 21, "rows of a sub");
    auto sc = s[0 .. 2, 0 .. 1];
    sc[0, 0] = 0;
    check(m[1, 1] == 11 && sc.nrow == 2, "slicing a Submatrix with [a..b, c..d] copies too");

    /* the data live as long as a view does */
    auto v = lonelyView();
    check(v[0, 0] == 11 && v[1, 1] == 22, "a view keeps the data alive after the matrix is gone");
    auto v2 = v;
    v = Submatrix.init;
    check(v2[1, 0] == 21, "copies of a view share the data");

    /* errors */
    check(collectException!GretlException(m.sub[0 .. 6, 0 .. 1]) !is null, "sub: rows out of range");
    check(collectException!GretlException(m.cols[3 .. 5]) !is null, "cols: out of range");
    check(collectException!GretlException(m.rows[3 .. 2]) !is null, "rows: reversed range");
    check(collectException!GretlException(s[3, 0]) !is null, "Submatrix index out of range");
    check(collectException!GretlException(m[5, 0]) !is null, "GretlMatrix index out of range");

    /* iteration */
    auto it = make(3, 2);
    double tot = 0;
    size_t nr;
    foreach (row; it.byRow) { nr++; tot += row[0, 0] + row[0, 1]; check(row.nrow == 1 && row.ncol == 2, "byRow yields 1 x ncol views"); }
    check(nr == 3 && tot == 0 + 1 + 10 + 11 + 20 + 21, "byRow visits every row once");
    size_t nc;
    foreach (col; it.byColumn) { nc++; check(col.nrow == 3 && col.ncol == 1, "byColumn yields nrow x 1 views"); col[1, 0] = -1; }
    check(nc == 2 && it[1, 0] == -1 && it[1, 1] == -1 && it[0, 0] == 0, "writing through a byColumn element changes the matrix");
    check(it.byRow.length == 3 && it.byColumn[1][2, 0] == 21 && it.byRow[2][0, 1] == 21, "length and indexing");
    auto blk = it.sub[1 .. 3, 0 .. 2];
    size_t nb;
    foreach (row; blk.byRow) nb += row.ncol;
    check(nb == 4 && blk.byRow.front[0, 1] == -1 && blk.byColumn.back[1, 0] == 21, "byRow and byColumn of a Submatrix");
    auto br = it.byRow;
    br.popFront();
    br.popBack();
    check(br.length == 1 && br.front[0, 0] == -1, "popFront and popBack");
    check(make(0, 3).byRow.empty && make(3, 0).byColumn.empty, "empty matrices");

    /* Row and Col as ranges of doubles */
    auto g = make(3, 4);
    double rs = 0;
    foreach (val; Row(g, 1)) rs += val;
    check(rs == 10 + 11 + 12 + 13, "foreach over Row");
    double cs = 0;
    foreach (val; Col(g, 2)) cs += val;
    check(cs == 2 + 12 + 22, "foreach over Col");
    size_t isum;
    double wsum = 0;
    foreach (jj, val; Row(g, 2)) { isum += jj; wsum += jj * val; }
    check(isum == 6 && wsum == 21 + 2 * 22 + 3 * 23, "foreach with an index");
    foreach (ref val; Col(g, 1)) val *= 2;
    check(g[0, 1] == 2 && g[1, 1] == 22 && g[2, 1] == 42 && g[1, 0] == 10, "foreach ref changes the matrix");
    check(Row(g, 0).length == 4 && Col(g, 0).length == 3 && Row(g, 1)[3] == 13 && Col(g, 3)[2] == 23, "length and indexing");
    Row(g, 0)[0] = 55;
    check(g[0, 0] == 55, "indexed assignment");
    check(Row(g, 1).dup == [10.0, 22.0, 12.0, 13.0], "dup gives a double[]");
    auto rg = Row(g, 1);
    rg.popFront();
    rg.popBack();
    check(rg.length == 2 && rg.front == 22 && rg.back == 12, "range primitives");
    int count;
    foreach (val; Row(g, 0)) { if (++count == 2) break; }
    check(count == 2, "break in foreach");
    auto sv = g.sub[1 .. 3, 1 .. 4];
    check(Row(sv, 1).dup == [42.0, 22.0, 23.0], "Row of a Submatrix uses its stride");
    check(Col(sv, 0).dup == [22.0, 42.0], "Col of a Submatrix");
    check(collectException!GretlException(Row(g, 3)) !is null && collectException!GretlException(Col(g, 4)) !is null, "Row/Col out of range");
    import std.algorithm : map, sum;
    check(Row(g, 2).map!(x => x * 2).sum == 2 * (20 + 42 + 22 + 23), "works with std.algorithm");

    /* el: all the elements */
    auto e = make(3, 4);
    double esum = 0;
    size_t ecount;
    foreach (val; e.el) { esum += val; ecount++; }
    check(ecount == 12 && esum == (0 + 10 + 20) * 4 + (0 + 1 + 2 + 3) * 3, "foreach over el visits every element once");
    bool idxOk = true;
    foreach (p, val; e.el) idxOk = idxOk && val == 10.0 * p.row + p.col;
    check(idxOk, "el.row and el.col index the element");
    foreach (ref val; e.el) val += 1;
    check(e[0, 0] == 1 && e[2, 3] == 24, "foreach ref over el changes the matrix");
    foreach (p, ref val; e.el) if (p.row == p.col) val = 0;
    check(e[1, 1] == 0 && e[2, 2] == 0 && e[0, 1] == 2, "assigning by position");
    auto eb = e.sub[1 .. 3, 1 .. 4];
    size_t nb2;
    bool okb = true;
    foreach (p, val; eb.el) { nb2++; okb = okb && val == e[p.row + 1, p.col + 1]; }
    check(nb2 == 6 && okb && eb.el.length == 6, "el of a Submatrix uses its stride and its own indices");
    int stop;
    foreach (val; e.el) { if (++stop == 5) break; }
    check(stop == 5, "break in el");
    size_t none;
    foreach (val; make(0, 3).el) none++;
    check(none == 0, "el of an empty matrix");

    /* operators */
    auto a = make(2, 3);                 // [0 1 2; 10 11 12]
    auto b = make(2, 3);
    foreach (p, ref val; b.el) val = 1 + p.row + 2 * p.col;   // [1 3 5; 2 4 6]
    auto sm = a + b;
    check(sm.nrow == 2 && sm.ncol == 3 && sm[0, 0] == 1 && sm[1, 2] == 18, "m + mm");
    check((a - b)[1, 1] == 7 && (a - b)[0, 2] == -3, "m - mm");
    check((a + 1)[1, 2] == 13 && (1 + a)[0, 0] == 1 && (a - 1)[0, 0] == -1 && (1 - a)[1, 0] == -9, "m +- a and a +- m");
    check((a * 2)[1, 1] == 22 && (2 * a)[1, 1] == 22 && (a / 2)[1, 2] == 6, "m * a, a * m, m / a");
    check((12 / b)[0, 0] == 12 && (12 / b)[1, 2] == 2, "a / m divides a by each element");
    auto sv1 = a.sub[0 .. 2, 1 .. 3] + b.cols[0 .. 2];
    check(sv1.nrow == 2 && sv1[0, 0] == 1 + 1 && sv1[1, 1] == 12 + 4, "Submatrix + Submatrix with strides");
    check((a.sub[1 .. 2, 0 .. 3] * 2)[0, 2] == 24 && (3 - a.rows[1 .. 2])[0, 0] == -7, "Submatrix and scalar");
    check((a + b.sub[0 .. 2, 0 .. 3])[1, 0] == 12, "GretlMatrix + Submatrix");
    check(collectException!GretlException(a + b.sub[0 .. 2, 0 .. 2]) !is null, "size mismatch in + throws");

    /* matrix multiplication: [1 2; 3 4] * [5 6; 7 8] = [19 22; 43 50] */
    auto p1 = make(2, 2), p2 = make(2, 2);
    p1[0, 0] = 1; p1[0, 1] = 2; p1[1, 0] = 3; p1[1, 1] = 4;
    p2[0, 0] = 5; p2[0, 1] = 6; p2[1, 0] = 7; p2[1, 1] = 8;
    auto pr = p1 * p2;
    check(pr[0, 0] == 19 && pr[0, 1] == 22 && pr[1, 0] == 43 && pr[1, 1] == 50, "m * mm (matrix product)");
    auto big = make(5, 4), big2 = make(4, 3);
    auto bp = big * big2;
    auto vp = big.sub[0 .. 5, 0 .. 4] * big2.sub[0 .. 4, 0 .. 3];
    bool same = bp.nrow == 5 && bp.ncol == 3 && vp.nrow == 5 && vp.ncol == 3;
    foreach (i; 0 .. 5) foreach (j; 0 .. 3) {
        double ref_ = 0;
        foreach (k; 0 .. 4) ref_ += big[i, k] * big2[k, j];
        same = same && (bp[i, j] - ref_) < 1e-9 && (ref_ - bp[i, j]) < 1e-9 && vp[i, j] == ref_;
    }
    check(same, "matrix product agrees for whole matrices (BLAS) and for views (loop)");
    check((p1.sub[0 .. 2, 1 .. 2] * p2.sub[0 .. 1, 0 .. 2]).nrow == 2 && (p1.sub[0 .. 2, 1 .. 2] * p2.sub[0 .. 1, 0 .. 2])[1, 1] == 24, "product of views with a stride");
    check(collectException!GretlException(big * big) !is null, "size mismatch in * throws");

    /* element by element */
    auto hm = a.el * b;
    check(hm[1, 2] == 12 * 6 && hm[0, 1] == 3, "m.el * mm");
    check((a * b.el)[1, 0] == 20 && (a.el * b.el)[1, 1] == 44, "m * mm.el and m.el * mm.el");
    check((b.el / b)[1, 1] == 1 && (a / b.el)[1, 2] == 2 && (a.el / b.el)[0, 2] == 2.0 / 5, "element-wise division");
    check((a.el * 2)[1, 1] == 22 && (2 * a.el)[1, 1] == 22 && (12 / b.el)[1, 2] == 2 && (a.el / 2)[0, 2] == 1, "elements with a scalar");
    check((a.sub[0 .. 2, 1 .. 3].el * b.sub[0 .. 2, 0 .. 2])[1, 1] == 12 * 4, "element-wise on views");
    check(collectException!GretlException(a.el * b.sub[0 .. 2, 0 .. 2]) !is null, "size mismatch in .el * throws");
    check((a.el + b)[1, 2] == 18 && (a + b.el)[0, 0] == 1 && (a.el - b.el)[1, 1] == 7 && (a.el + 1)[0, 0] == 1 && (1 - a.el)[1, 0] == -9,
          ".el with + and - is the same as without");

    /* right division: x * mm = m. Square: [1 2; 3 4] / [5 6; 7 8] = [3 -2; 2 -1] */
    auto rd = p1 / p2;
    check(rd.nrow == 2 && rd.ncol == 2 && close(rd[0, 0], 3) && close(rd[0, 1], -2) && close(rd[1, 0], 2) && close(rd[1, 1], -1), "m / mm for a square mm");
    auto back = rd * p2;
    check(close(back[0, 0], 1) && close(back[0, 1], 2) && close(back[1, 0], 3) && close(back[1, 1], 4), "(m / mm) * mm gives m back");

    /* OLS: y a 1 x n row, X' a p x n matrix, y / X' = (X'X)^-1 X'y as a row. y = 1 + 2 x exactly, with noise-free data */
    auto Xt = make(2, 5), yrow = make(1, 5);
    foreach (t; 0 .. 5) { Xt[0, t] = 1; Xt[1, t] = t; yrow[0, t] = 1 + 2.0 * t; }
    auto beta = yrow / Xt;
    check(beta.nrow == 1 && beta.ncol == 2 && close(beta[0, 0], 1) && close(beta[0, 1], 2), "y / X' is the OLS coefficient row (exact data)");
    yrow[0, 2] = yrow[0, 2] + 1.0;                    // add an error: compare with the normal equations
    auto b2 = yrow / Xt;
    double sxx = 0, sx = 0, sy = 0, sxy = 0;
    foreach (t; 0 .. 5) { sx += t; sy += yrow[0, t]; sxx += t * t; sxy += t * yrow[0, t]; }
    double slope = (5 * sxy - sx * sy) / (5 * sxx - sx * sx), icpt = (sy - slope * sx) / 5;
    check(close(b2[0, 1], slope) && close(b2[0, 0], icpt), "y / X' matches the OLS formulas with an error term");
    check((yrow.sub[0 .. 1, 0 .. 5] / Xt.sub[0 .. 2, 0 .. 5]).ncol == 2 && close((yrow.sub[0 .. 1, 0 .. 5] / Xt.sub[0 .. 2, 0 .. 5])[0, 1], slope), "m / mm on views");
    /* solve(X, y) for data in columns */
    auto Xc = make(5, 2), yc = make(5, 1);
    foreach (t; 0 .. 5) { Xc[t, 0] = 1; Xc[t, 1] = t; yc[t, 0] = yrow[0, t]; }
    auto bs = solve(Xc, yc);
    check(bs.nrow == 2 && bs.ncol == 1 && close(bs[0, 0], icpt) && close(bs[1, 0], slope), "solve(X, y) is the OLS coefficient column");
    check(close(bs[0, 0], b2[0, 0]) && close(bs[1, 0], b2[0, 1]), "solve(X, y) agrees with y' / X'");
    auto bs2 = solve(Xc.sub[0 .. 5, 0 .. 2], yc.rows[0 .. 5]);
    check(close(bs2[1, 0], slope), "solve on views");
    auto Ym = make(5, 2);
    foreach (t; 0 .. 5) { Ym[t, 0] = yc[t, 0]; Ym[t, 1] = 3 - 0.5 * t; }
    auto bm = solve(Xc, Ym);
    check(bm.nrow == 2 && bm.ncol == 2 && close(bm[0, 1], 3) && close(bm[1, 1], -0.5) && close(bm[1, 0], slope), "solve with several columns in y");
    auto sq = solve(p2, p1);
    check(close((p2 * sq)[0, 0], 1) && close((p2 * sq)[1, 1], 4), "solve with a square X gives X \\ y");
    check(collectException!GretlException(solve(Xc, make(4, 1))) !is null, "solve with different row counts throws");
    check(collectException!GretlException(solve(make(1, 2), make(1, 1))) !is null, "underdetermined solve throws");
    /* t, inv, kron */
    auto tr = t(a);                       // a is 2 x 3
    check(tr.nrow == 3 && tr.ncol == 2 && tr[2, 1] == 12 && tr[1, 0] == 1 && t(tr)[1, 2] == 12, "t(m)");
    check(t(a.sub[0 .. 2, 1 .. 3])[1, 1] == 12 && a.t[0, 1] == 10, "t of a view and as a method");
    auto pi = inv(p1);                    // [1 2; 3 4]^-1 = [-2 1; 1.5 -0.5]
    check(close(pi[0, 0], -2) && close(pi[0, 1], 1) && close(pi[1, 0], 1.5) && close(pi[1, 1], -0.5), "inv(m)");
    auto id = p1 * p1.inv;
    check(close(id[0, 0], 1) && close(id[1, 1], 1) && close(id[0, 1], 0) && close(id[1, 0], 0), "m * inv(m) is the identity");
    check(close(inv(p1.sub[0 .. 2, 0 .. 2])[1, 0], 1.5), "inv of a view");
    check(collectException!GretlException(inv(a)) !is null, "inv of a non-square matrix throws");
    auto sing = make(2, 2);
    sing[0, 0] = 1; sing[0, 1] = 2; sing[1, 0] = 2; sing[1, 1] = 4;
    check(collectException!GretlException(inv(sing)) !is null, "inv of a singular matrix throws");
    auto kp = kron(p1, p2);               // [1 2; 3 4] (x) [5 6; 7 8]
    check(kp.nrow == 4 && kp.ncol == 4 && kp[0, 0] == 5 && kp[0, 3] == 12 && kp[1, 2] == 14 && kp[2, 1] == 18 && kp[3, 3] == 32, "kron(a, b)");
    auto kr = a.kron(b.sub[0 .. 1, 0 .. 2]);
    check(kr.nrow == 2 && kr.ncol == 6 && kr[1, 4] == 12 && kr[0, 5] == 2 * 3, "a.kron(b) with a view, sizes");
    check(kron(make(0, 2), p1).nrow == 0 && kron(p1, make(2, 0)).ncol == 0, "kron with an empty matrix");

    /* powers */
    auto sq2 = p1.el ^^ 2;
    check(sq2[0, 0] == 1 && sq2[0, 1] == 4 && sq2[1, 0] == 9 && sq2[1, 1] == 16, "m.el ^^ 2 squares each element");
    check((2 ^^ p1.el)[1, 1] == 16 && (p1.el ^^ p1.el)[1, 0] == 27, "a ^^ m.el and m.el ^^ mm.el");
    check((p1.el ^^ 0.5)[1, 1] == 2 && (p1.el ^^ p2.el)[0, 1] == 64, "m.el ^^ a, m.el ^^ mm.el");
    check(!__traits(compiles, p1 ^^ 2) && !__traits(compiles, 2 ^^ p1), "m ^^ n and a ^^ m are refused");
    auto pw = power(p1, 2), pm = p1 * p1;
    check(pw[0, 0] == pm[0, 0] && pw[0, 1] == pm[0, 1] && pw[1, 0] == pm[1, 0] && pw[1, 1] == pm[1, 1], "power(m, 2) is m * m");
    auto p5 = power(p1, 5), pm5 = p1 * p1 * p1 * p1 * p1;
    check(close(p5[0, 0], pm5[0, 0]) && close(p5[1, 1], pm5[1, 1]) && close(power(p1, 1)[1, 0], 3), "power(m, 5), power(m, 1)");
    check(power(p1, 0)[0, 0] == 1 && power(p1, 0)[0, 1] == 0 && close(power(p1, -1)[0, 0], -2) && close(power(p1, -2)[1, 1], inv(p1 * p1)[1, 1]), "power(m, 0) and negative powers");
    check(close(power(p1.sub[0 .. 2, 0 .. 2], 3)[0, 0], power(p1, 3)[0, 0]), "power of a view");
    check(collectException!GretlException(power(a, 2)) !is null, "power of a non-square matrix throws");

    /* chol: [4 2; 2 3] = U'U, U = [2 1; 0 sqrt(2)] */
    auto spd = make(2, 2);
    spd[0, 0] = 4; spd[0, 1] = 2; spd[1, 0] = 2; spd[1, 1] = 3;
    auto U = chol(spd);
    check(close(U[0, 0], 2) && close(U[0, 1], 1) && U[1, 0] == 0 && close(U[1, 1], 1.4142135623730951), "chol(m) is the upper triangle");
    auto L = spd.chol(Triangle.lower);
    check(close(L[0, 0], 2) && close(L[1, 0], 1) && L[0, 1] == 0 && close(L[1, 1], 1.4142135623730951), "chol(m, Triangle.lower) is the lower triangle");
    auto ut = t(U) * U, ll = L * t(L);
    check(close(ut[0, 1], 2) && close(ut[1, 1], 3) && close(ll[0, 0], 4) && close(ll[1, 1], 3), "U'U and L L' give m back");
    check(close(chol(spd.sub[0 .. 2, 0 .. 2])[1, 1], 1.4142135623730951), "chol of a view");
    check(collectException!GretlException(chol(p1 - make(2, 2) * 100)) !is null, "chol of a matrix that is not positive definite throws");
    check(collectException!GretlException(chol(a)) !is null, "chol of a non-square matrix throws");

    /* diag */
    auto dm = make(3, 3);
    diag(dm) = [1.1, 3.4, 1.7];
    check(dm[0, 0] == 1.1 && dm[1, 1] == 3.4 && dm[2, 2] == 1.7 && dm[0, 1] == 1, "diag(m) = [..] sets the diagonal");
    check(diag(dm).length == 3 && diag(dm)[1] == 3.4 && diag(dm).dup == [1.1, 3.4, 1.7], "diag(m) reads");
    diag(dm) = Row(dm, 1);          // row 1 is [10, 3.4, 12], overlapping the diagonal
    check(dm[0, 0] == 10 && dm[1, 1] == 3.4 && dm[2, 2] == 12, "diag(m) = Row(m, 1), with overlap");
    diag(dm) = 0;
    check(dm[0, 0] == 0 && dm[1, 1] == 0 && dm[2, 2] == 0 && dm[2, 1] == 21, "diag(m) = scalar");
    foreach (ref val; diag(dm)) val += 7;
    check(dm[1, 1] == 7 && dm[0, 1] == 1, "foreach ref over diag");
    auto dn = make(2, 4);
    check(diag(dn).length == 2 && diag(dn).dup == [0.0, 11.0] && diag(make(4, 2)).dup == [0.0, 11.0], "diag of a non-square matrix");
    check(diag(dm.sub[1 .. 3, 1 .. 3]).dup == [7.0, 7.0] && diag(dm.sub[1 .. 3, 0 .. 3]).dup == [10.0, 21.0], "diag of a view");
    diag(dm.sub[1 .. 3, 0 .. 3]) = [5.0, 6.0];
    check(dm[1, 0] == 5 && dm[2, 1] == 6, "assigning to the diagonal of a view");
    Row(dm, 0) = [9.0, 8.0, 7.0];
    check(dm[0, 0] == 9 && dm[0, 2] == 7, "Row(m, i) = [..] assigns too");
    check(collectException!GretlException(diag(dm) = [1.0, 2.0]) !is null, "assigning the wrong number of values throws");
    check(diag(make(0, 3)).length == 0, "diag of an empty matrix");
    auto identity = make(2, 2);
    identity = identity * 0;
    diag(identity) = 1;
    check(identity[0, 0] == 1 && identity[0, 1] == 0 && identity[1, 1] == 1, "identity matrix idiom");

    /* eigen */
    auto s2 = make(2, 2);
    s2[0, 0] = 2; s2[0, 1] = 1; s2[1, 0] = 1; s2[1, 1] = 2;
    auto ev = eigen(s2);
    check(ev.values.nrow == 2 && ev.values.ncol == 1 && close(ev.values[0, 0], 1) && close(ev.values[1, 0], 3), "eigen: values, ascending");
    check(ev.vectors.nrow == 2 && ev.vectors.ncol == 2 && close(abs(ev.vectors[0, 0]), 0.7071067811865476) && close(ev.vectors[0, 0] * ev.vectors[1, 0], -0.5), "eigen: first eigenvector is (1, -1)/sqrt 2 up to sign");
    auto s4 = make(4, 4);
    foreach (i; 0 .. 4) foreach (j; 0 .. 4) s4[i, j] = 1.0 / (1 + i + j) + (i == j ? i : 0);   // a symmetric matrix
    auto e4 = eigen(s4);
    auto vv = e4.vectors;
    auto resid = s4 * vv;                   // s4 v = v diag(values)
    bool eok = true, ascending = true, orth = true;
    foreach (i; 0 .. 4) foreach (j; 0 .. 4) eok = eok && close(resid[i, j], vv[i, j] * e4.values[j, 0]);
    foreach (j; 1 .. 4) ascending = ascending && e4.values[j - 1, 0] <= e4.values[j, 0];
    auto vtv = t(vv) * vv;
    foreach (i; 0 .. 4) foreach (j; 0 .. 4) orth = orth && close(vtv[i, j], i == j ? 1 : 0);
    check(eok && ascending && orth, "eigen: m v = v diag(values), ascending, orthonormal (4 x 4)");
    double trSum = 0, sumev = 0;
    foreach (j; 0 .. 4) { trSum += s4[j, j]; sumev += e4.values[j, 0]; }
    check(close(trSum, sumev), "the eigenvalues add up to the trace");
    auto novec = eigen(s4, false);
    check(close(novec.values[3, 0], e4.values[3, 0]) && novec.vectors.nrow == 0, "eigen without vectors");
    check(close(eigen(s4.sub[0 .. 4, 0 .. 4]).values[0, 0], e4.values[0, 0]) && close(s4.eigen.values[1, 0], e4.values[1, 0]), "eigen of a view, as a method");
    check(close(s4[0, 0], 1) && close(s4[3, 3], 3 + 1.0 / 7), "eigen leaves the matrix alone");
    check(collectException!GretlException(eigen(p1)) !is null, "eigen of a non-symmetric matrix throws");
    check(collectException!GretlException(eigen(a)) !is null, "eigen of a non-square matrix throws");
    check(eigen(make(0, 0)).values.nrow == 0, "eigen of an empty matrix");

    /* eigenGeneral */
    auto rot = make(2, 2);                      // [0 -1; 1 0]: eigenvalues +-i
    rot[0, 0] = 0; rot[0, 1] = -1; rot[1, 0] = 1; rot[1, 1] = 0;
    auto eg = eigenGeneral(rot);
    check(eg.values.nrow == 2 && eg.values.ncol == 2 && close(eg.values[0, 0], 0) && close(abs(eg.values[0, 1]), 1) && close(eg.values[0, 1], -eg.values[1, 1]) && !eg.isReal,
          "eigenGeneral: complex conjugate pair for a rotation");
    /* m v = lambda v with v = vr + i vi and lambda = a + i b:  m vr = a vr - b vi,  m vi = a vi + b vr */
    bool rok = true;
    auto mvr = rot * eg.vectorsRe, mvi = rot * eg.vectorsIm;
    foreach (j; 0 .. 2) foreach (i; 0 .. 2) {
        const la = eg.values[j, 0], lb = eg.values[j, 1];
        rok = rok && close(mvr[i, j], la * eg.vectorsRe[i, j] - lb * eg.vectorsIm[i, j]) && close(mvi[i, j], la * eg.vectorsIm[i, j] + lb * eg.vectorsRe[i, j]);
    }
    check(rok, "eigenGeneral: m v = lambda v for complex eigenpairs");
    auto tri = make(3, 3);                      // upper triangular, eigenvalues on the diagonal
    foreach (i; 0 .. 3) foreach (j; 0 .. 3) tri[i, j] = i > j ? 0 : (i == j ? 2.0 + 3 * i : 1.0);
    auto et = eigenGeneral(tri);
    double sumRe = 0, nComplex = 0;
    foreach (i; 0 .. 3) { sumRe += et.values[i, 0]; nComplex += et.values[i, 1] == 0 ? 0 : 1; }
    check(et.isReal && close(sumRe, 2 + 5 + 8) && nComplex == 0, "eigenGeneral: real eigenvalues of a triangular matrix");
    bool tok = true;
    auto mvt = tri * et.vectorsRe;
    foreach (j; 0 .. 3) foreach (i; 0 .. 3) tok = tok && close(mvt[i, j], et.values[j, 0] * et.vectorsRe[i, j]) && et.vectorsIm[i, j] == 0;
    check(tok, "eigenGeneral: real eigenvectors");
    auto e1 = eigenGeneral(p1);                 // [1 2; 3 4]: (5 +- sqrt 33)/2
    const lo = (5 - 5.744562646538029) / 2, hi = (5 + 5.744562646538029) / 2;
    check(e1.isReal && ((close(e1.values[0, 0], lo) && close(e1.values[1, 0], hi)) || (close(e1.values[0, 0], hi) && close(e1.values[1, 0], lo))), "eigenGeneral of [1 2; 3 4]");
    auto one = eigenGeneral(make(1, 1));
    check(one.values.nrow == 1 && one.vectorsRe[0, 0] == 1 && one.isReal, "eigenGeneral of a 1 x 1 matrix");
    check(eigenGeneral(tri, false).vectorsRe.nrow == 0 && eigenGeneral(tri, false).values.nrow == 3, "eigenGeneral without vectors");
    check(close(eigenGeneral(tri.sub[0 .. 3, 0 .. 3]).values[0, 0], et.values[0, 0]) && close(tri.eigenGeneral.values[2, 0], et.values[2, 0]), "eigenGeneral of a view, as a method");
    check(collectException!GretlException(eigenGeneral(a)) !is null, "eigenGeneral of a non-square matrix throws");
    check(eigenGeneral(make(0, 0)).values.nrow == 0, "eigenGeneral of an empty matrix");

    /* svd */
    auto dg = make(2, 2);
    dg[0, 0] = 3; dg[0, 1] = 0; dg[1, 0] = 0; dg[1, 1] = 4;
    auto sd = svd(dg);
    check(sd.s.nrow == 2 && sd.s.ncol == 1 && close(sd.s[0, 0], 4) && close(sd.s[1, 0], 3), "svd: singular values of diag(3, 4), largest first");
    foreach (shape; [[3, 3], [5, 3], [3, 5], [6, 1], [1, 4]]) {
        const nrw = shape[0], ncl = shape[1], k = nrw < ncl ? nrw : ncl;
        auto am = make(nrw, ncl);
        foreach (i; 0 .. nrw) foreach (j; 0 .. ncl) am[i, j] = 1.0 / (1 + i + 2 * j) + (i == j ? 1 : 0) - 0.1 * i * j;
        auto svr = svd(am);
        bool shapeOk = svr.u.nrow == nrw && svr.u.ncol == k && svr.v.nrow == ncl && svr.v.ncol == k && svr.s.nrow == k;
        auto us = make(nrw, k);
        foreach (i; 0 .. nrw) foreach (j; 0 .. k) us[i, j] = svr.u[i, j] * svr.s[j, 0];
        auto rec = us * t(svr.v);
        bool recOk = true, descending = true, orthU = true, orthV = true;
        foreach (i; 0 .. nrw) foreach (j; 0 .. ncl) recOk = recOk && close(rec[i, j], am[i, j]);
        foreach (j; 1 .. k) descending = descending && svr.s[j - 1, 0] >= svr.s[j, 0];
        auto utu = t(svr.u) * svr.u, vtv2 = t(svr.v) * svr.v;
        foreach (i; 0 .. k) foreach (j; 0 .. k) { orthU = orthU && close(utu[i, j], i == j ? 1 : 0); orthV = orthV && close(vtv2[i, j], i == j ? 1 : 0); }
        check(shapeOk && recOk && descending && orthU && orthV, "svd: shapes, reconstruction, order, orthonormality for " ~ to!string(nrw) ~ " x " ~ to!string(ncl));
    }
    auto fsv = svd(make(5, 3), true, true);
    auto fu = t(fsv.u) * fsv.u;
    check(fsv.u.nrow == 5 && fsv.u.ncol == 5 && fsv.v.nrow == 3 && fsv.v.ncol == 3 && close(fu[4, 4], 1) && close(fu[0, 4], 0) && fsv.s.nrow == 3, "svd full: u is 5 x 5, v is 3 x 3");
    auto fsw = svd(make(3, 5), true, true);
    check(fsw.u.nrow == 3 && fsw.u.ncol == 3 && fsw.v.nrow == 5 && fsw.v.ncol == 5, "svd full of a wide matrix");
    auto ns = svd(make(5, 3), false);
    check(ns.u.nrow == 0 && ns.v.nrow == 0 && ns.s.nrow == 3 && close(ns.s[0, 0], svd(make(5, 3)).s[0, 0]), "svd without vectors");
    /* accuracy for a tall, nearly rank-deficient matrix: smallest singular value is about 1e-9 / sqrt 2 */
    auto nr3 = make(3, 2);
    nr3[0, 0] = 1; nr3[0, 1] = 1; nr3[1, 0] = 0; nr3[1, 1] = 1e-9; nr3[2, 0] = 0; nr3[2, 1] = 0;
    auto sn = svd(nr3, false);
    check(close(sn.s[0, 0], 1.4142135623730951) && abs(sn.s[1, 0] - 7.0710678118654752e-10) < 1e-12, "svd of a tall, nearly singular matrix keeps the small singular value");
    check(svd(p1.sub[0 .. 2, 0 .. 2]).s.nrow == 2 && close(p1.svd.s[0, 0], 5.464985704219043), "svd of a view, as a method");
    check(svd(make(0, 3)).s.nrow == 0, "svd of an empty matrix");

    /* qr */
    foreach (shape; [[3, 3], [6, 3], [5, 1]]) {
        const qnr = shape[0], qnc = shape[1];
        auto qa = make(qnr, qnc);
        foreach (i; 0 .. qnr) foreach (j; 0 .. qnc) qa[i, j] = 1.0 / (1 + i + 2 * j) + (i == j ? 1 : 0) - 0.1 * i * j;
        auto dec = qr(qa);
        bool qshape = dec.q.nrow == qnr && dec.q.ncol == qnc && dec.r.nrow == qnc && dec.r.ncol == qnc;
        auto prod = dec.q * dec.r;
        auto qtq = t(dec.q) * dec.q;
        bool rec = true, orthq = true, upper = true;
        foreach (i; 0 .. qnr) foreach (j; 0 .. qnc) rec = rec && close(prod[i, j], qa[i, j]);
        foreach (i; 0 .. qnc) foreach (j; 0 .. qnc) { orthq = orthq && close(qtq[i, j], i == j ? 1 : 0); if (i > j) upper = upper && dec.r[i, j] == 0; }
        check(qshape && rec && orthq && upper, "qr: shapes, q r = m, orthonormal q, upper triangular r for " ~ to!string(qnr) ~ " x " ~ to!string(qnc));
    }
    auto q2 = qr(p1);                                  // [1 2; 3 4]: |r00| = sqrt 10, |r11| = 2 / sqrt 10 * ... = determinant / |r00|
    check(close(abs(q2.r[0, 0]), 3.1622776601683795) && close(abs(q2.r[0, 0] * q2.r[1, 1]), 2), "qr: |r00| and |det| of [1 2; 3 4]");
    check(close(qr(p1.sub[0 .. 2, 0 .. 2]).r[1, 1], q2.r[1, 1]) && close(p1.qr.r[0, 0], q2.r[0, 0]), "qr of a view, as a method");
    auto input = make(3, 2);
    auto before = input[2, 1];
    qr(input);
    check(input[2, 1] == before, "qr leaves the matrix alone");
    check(collectException!GretlException(qr(make(2, 3))) !is null, "qr of a wide matrix throws");
    check(qr(make(0, 0)).q.nrow == 0, "qr of an empty matrix");

    /* lu */
    foreach (shape; [[3, 3], [5, 3], [3, 5], [1, 4]]) {
        const lnr = shape[0], lnc = shape[1], lk = lnr < lnc ? lnr : lnc;
        auto la = make(lnr, lnc);
        foreach (i; 0 .. lnr) foreach (j; 0 .. lnc) la[i, j] = (i + 2 * j) % 5 - 1.5 + (i == j ? 0.25 : 0) + 0.1 * i * j;
        auto f = lu(la);
        bool lshape = f.l.nrow == lnr && f.l.ncol == lk && f.u.nrow == lk && f.u.ncol == lnc && f.p.nrow == lnr && f.perm.length == lnr;
        auto pa = f.p * la, lu_ = f.l * f.u;
        bool prod = true, unit = true, trian = true, rows = true;
        foreach (i; 0 .. lnr) foreach (j; 0 .. lnc) prod = prod && close(pa[i, j], lu_[i, j]);
        foreach (i; 0 .. lnr) foreach (j; 0 .. lk) { if (i == j) unit = unit && f.l[i, j] == 1; if (j > i) trian = trian && f.l[i, j] == 0; }
        foreach (i; 0 .. lk) foreach (j; 0 .. lnc) if (i > j) trian = trian && f.u[i, j] == 0;
        foreach (i; 0 .. lnr) foreach (j; 0 .. lnc) rows = rows && pa[i, j] == la[f.perm[i], j];
        check(lshape && prod && unit && trian && rows && !f.singular, "lu: p m = l u, triangular factors and perm for " ~ to!string(lnr) ~ " x " ~ to!string(lnc));
    }
    auto lf = lu(p1);                                  // [1 2; 3 4]: pivots on the 3
    check(lf.perm == [cast(size_t) 1, 0] && close(lf.u[0, 0], 3) && close(lf.u[1, 1], 2.0 / 3) && close(lf.l[1, 0], 1.0 / 3), "lu of [1 2; 3 4]: the pivot, l and u");
    auto sg = make(2, 2);
    sg[0, 0] = 1; sg[0, 1] = 2; sg[1, 0] = 2; sg[1, 1] = 4;
    check(lu(sg).singular && !lu(p1).singular, "lu flags a singular matrix");
    check(close(lu(p1.sub[0 .. 2, 0 .. 2]).u[1, 1], 2.0 / 3) && close(p1.lu.l[1, 0], 1.0 / 3), "lu of a view, as a method");
    check(lu(make(0, 3)).perm.length == 0 && lu(make(0, 3)).u.ncol == 3, "lu of an empty matrix");

    /* schur */
    foreach (sz; [2, 3, 5]) {
        auto sa = make(sz, sz);
        foreach (i; 0 .. sz) foreach (j; 0 .. sz) sa[i, j] = ((i * 3 + j * 5) % 7) - 2.5 + (i == j ? 0.5 : 0) + (j == i + 1 ? 2 : 0) - (i == j + 1 ? 2 : 0);
        auto sch = schur(sa);
        auto back2 = sch.z * sch.t * t(sch.z);
        auto ztz = t(sch.z) * sch.z;
        bool sback = true, sorth = true, quasi = true;
        foreach (i; 0 .. sz) foreach (j; 0 .. sz) { sback = sback && close(back2[i, j], sa[i, j]); sorth = sorth && close(ztz[i, j], i == j ? 1 : 0); }
        foreach (i; 2 .. sz) foreach (j; 0 .. i - 1) quasi = quasi && sch.t[i, j] == 0;      // below the first subdiagonal
        double trT = 0, trA = 0, sumReal = 0;
        foreach (i; 0 .. sz) { trT += sch.t[i, i]; trA += sa[i, i]; sumReal += sch.values[i, 0]; }
        check(sback && sorth && quasi && close(trT, trA) && close(sumReal, trA) && sch.values.ncol == 2, "schur: z t z' = m, orthogonal z, quasi-triangular t, eigenvalues, size " ~ to!string(sz));
    }
    auto sr = schur(rot);                              // [0 -1; 1 0]: eigenvalues +-i, one 2 x 2 block
    check(close(abs(sr.values[0, 1]), 1) && close(sr.values[0, 0], 0), "schur: a rotation keeps its complex pair as a 2 x 2 block");
    auto sup = schur(tri);                             // already upper triangular
    check(close(sup.t[0, 0] + sup.t[1, 1] + sup.t[2, 2], 15) && sup.t[1, 0] == 0 && sup.t[2, 0] == 0 && sup.t[2, 1] == 0, "schur of a triangular matrix is triangular");
    check(schur(tri, false).z.nrow == 0 && close(schur(tri, false).t[0, 0], sup.t[0, 0]), "schur without z");
    check(close(schur(p1.sub[0 .. 2, 0 .. 2]).t[0, 0] + schur(p1).t[1, 1], 5) && close(p1.schur.t[0, 1] * 0, 0), "schur of a view, as a method");
    check(collectException!GretlException(schur(a)) !is null, "schur of a non-square matrix throws");
    check(schur(make(0, 0)).t.nrow == 0, "schur of an empty matrix");

    /* pinv */
    auto pp = pinv(p1), ip = inv(p1);
    check(pp.nrow == 2 && pp.ncol == 2 && close(pp[0, 0], ip[0, 0]) && close(pp[0, 1], ip[0, 1]) && close(pp[1, 0], ip[1, 0]) && close(pp[1, 1], ip[1, 1]), "pinv of a nonsingular matrix is inv");
    auto rd2 = make(2, 2);                                // rank 1: [1 2; 2 4], pinv = m / 25
    rd2[0, 0] = 1; rd2[0, 1] = 2; rd2[1, 0] = 2; rd2[1, 1] = 4;
    auto pr2 = pinv(rd2);
    check(close(pr2[0, 0], 1.0 / 25) && close(pr2[0, 1], 2.0 / 25) && close(pr2[1, 0], 2.0 / 25) && close(pr2[1, 1], 4.0 / 25), "pinv of a rank-1 matrix");
    foreach (shape; [[3, 2], [2, 3], [5, 3], [1, 4], [4, 1]]) {
        const pnr = shape[0], pnc = shape[1];
        auto pa2 = make(pnr, pnc);
        foreach (i; 0 .. pnr) foreach (j; 0 .. pnc) pa2[i, j] = 1.0 / (1 + i + 2 * j) + (i == j ? 1 : 0) - 0.1 * i * j;
        auto pi2 = pinv(pa2);
        auto apa = pa2 * pi2 * pa2, pap = pi2 * pa2 * pi2, ap = pa2 * pi2, pai = pi2 * pa2;
        bool pshape = pi2.nrow == pnc && pi2.ncol == pnr, mp1 = true, mp2 = true, sym1 = true, sym2 = true;
        foreach (i; 0 .. pnr) foreach (j; 0 .. pnc) mp1 = mp1 && close(apa[i, j], pa2[i, j]);
        foreach (i; 0 .. pnc) foreach (j; 0 .. pnr) mp2 = mp2 && close(pap[i, j], pi2[i, j]);
        foreach (i; 0 .. pnr) foreach (j; 0 .. pnr) sym1 = sym1 && close(ap[i, j], ap[j, i]);
        foreach (i; 0 .. pnc) foreach (j; 0 .. pnc) sym2 = sym2 && close(pai[i, j], pai[j, i]);
        check(pshape && mp1 && mp2 && sym1 && sym2, "pinv: the four Moore-Penrose conditions for " ~ to!string(pnr) ~ " x " ~ to!string(pnc));
    }
    auto tallm = make(4, 2);
    foreach (i; 0 .. 4) { tallm[i, 0] = 1; tallm[i, 1] = i; }
    auto viaN = inv(t(tallm) * tallm) * t(tallm), viaP = pinv(tallm);
    bool tallOk = true;
    foreach (i; 0 .. 2) foreach (j; 0 .. 4) tallOk = tallOk && close(viaN[i, j], viaP[i, j]);
    check(tallOk, "pinv of a tall full-rank matrix is inv(m'm) m'");
    auto zm = make(2, 3);
    zm = zm * 0;
    check(pinv(zm).nrow == 3 && pinv(zm).ncol == 2 && pinv(zm)[2, 1] == 0, "pinv of the zero matrix is zero");
    auto nearly = make(2, 2);
    nearly[0, 0] = 1; nearly[0, 1] = 0; nearly[1, 0] = 0; nearly[1, 1] = 1e-8;
    check(abs(pinv(nearly)[1, 1] - 1e8) < 1, "pinv keeps a small singular value by default");
    check(pinv(nearly, 1e-6)[1, 1] == 0, "a larger tol drops it");
    check(close(pinv(p1.sub[0 .. 2, 0 .. 2])[0, 0], ip[0, 0]) && close(p1.pinv[1, 1], ip[1, 1]), "pinv of a view, as a method");
    check(pinv(make(0, 3)).nrow == 3 && pinv(make(0, 3)).ncol == 0, "pinv of an empty matrix");

    /* det, logdet, trace */
    check(close(det(p1), -2) && close(p1.det, -2), "det of a 2 x 2 matrix");
    auto m33 = make(3, 3);
    foreach (i; 0 .. 3) foreach (j; 0 .. 3) m33[i, j] = [2.0, 1, 0, 1, 3, 1, 0, 1, 4][j * 3 + i];   // symmetric tridiagonal, det = 18
    check(close(det(m33), 18), "det of a 3 x 3 matrix");
    auto m5 = make(5, 5);
    foreach (i; 0 .. 5) foreach (j; 0 .. 5) m5[i, j] = (i == j ? 3.0 : 0) + (i + j) % 3 * 0.5 - (i > j ? 0.25 : 0);
    auto d5 = det(m5);
    auto lu5 = lu(m5);
    double d5b = 1;
    foreach (i; 0 .. 5) d5b *= lu5.u[i, i];
    int swaps;
    foreach (i; 0 .. 5) { foreach (j; i + 1 .. 5) if (lu5.perm[j] < lu5.perm[i]) swaps++; }
    if (swaps % 2) d5b = -d5b;
    check(close(d5, d5b), "det of a 5 x 5 matrix agrees with the LU factors");
    check(det(sing) == 0 && det(make(0, 0)) == 1, "det of a singular matrix is 0, of a 0 x 0 matrix 1");
    check(close(det(p1 * p2), det(p1) * det(p2)) && close(det(t(m5)), d5), "det of a product, det of the transpose");
    check(close(det(p1.sub[0 .. 2, 0 .. 2]), -2), "det of a view");
    check(collectException!GretlException(det(a)) !is null, "det of a non-square matrix throws");
    check(close(logdet(m33), log(18.0)) && close(m33.logdet, log(18.0)), "logdet is the log of the determinant");
    check(close(logdet(m5, true), log(fabs(d5))) && (d5 > 0 ? close(logdet(m5), log(d5)) : isNaN(logdet(m5))), "logdet of a 5 x 5 matrix");
    check(isNaN(logdet(p1)) && close(logdet(p1, true), log(2.0)), "a negative determinant: NaN, or the log of the absolute value");
    check(logdet(sing) == -double.infinity && logdet(sing, true) == -double.infinity, "logdet of a singular matrix is -infinity");
    auto bigm = make(400, 400);
    foreach (i; 0 .. 400) foreach (j; 0 .. 400) bigm[i, j] = i == j ? 100.0 : 0.0;
    check(close(logdet(bigm), 400 * log(100.0)) && det(bigm) == double.infinity, "logdet stays finite where the determinant overflows");
    check(logdet(make(0, 0)) == 0 && collectException!GretlException(logdet(a)) !is null, "logdet of a 0 x 0 matrix, of a non-square matrix");
    check(trace(p1) == 5 && p2.trace == 13 && trace(p1.sub[0 .. 2, 0 .. 2]) == 5 && trace(make(0, 0)) == 0, "trace");
    check(trace(m5 * 2) == 35 && close(trace(p1 * p2), trace(p2 * p1)) && trace(m33.sub[1 .. 3, 1 .. 3]) == 7, "trace of a scaled matrix and of products and views");
    check(collectException!GretlException(trace(a)) !is null, "trace of a non-square matrix throws");

    /* expm */
    auto ez = expm(make(3, 3) * 0);
    bool eid = true;
    foreach (i; 0 .. 3) foreach (j; 0 .. 3) eid = eid && close(ez[i, j], i == j ? 1 : 0);
    check(eid, "expm of the zero matrix is the identity");
    auto edg = make(2, 2);
    edg[0, 0] = 1; edg[0, 1] = 0; edg[1, 0] = 0; edg[1, 1] = -2;
    auto ed = expm(edg);
    check(close(ed[0, 0], exp(1.0)) && close(ed[1, 1], exp(-2.0)) && close(ed[0, 1], 0) && close(ed[1, 0], 0), "expm of a diagonal matrix");
    auto nil = make(2, 2);
    nil[0, 0] = 0; nil[0, 1] = 1; nil[1, 0] = 0; nil[1, 1] = 0;
    auto en = expm(nil);
    check(close(en[0, 0], 1) && close(en[0, 1], 1) && close(en[1, 0], 0) && close(en[1, 1], 1), "expm of a nilpotent matrix is I + m");
    auto th = make(2, 2);
    th[0, 0] = 0; th[0, 1] = -0.7; th[1, 0] = 0.7; th[1, 1] = 0;
    auto er = expm(th);
    check(close(er[0, 0], cos(0.7)) && close(er[0, 1], -sin(0.7)) && close(er[1, 0], sin(0.7)) && close(er[1, 1], cos(0.7)), "expm of a skew matrix is a rotation");
    auto ei = expm(p1) * expm(p1 * -1);
    check(close(ei[0, 0], 1) && close(ei[1, 1], 1) && abs(ei[0, 1]) < 1e-8 && abs(ei[1, 0]) < 1e-8, "expm(m) * expm(-m) is the identity");
    check(close(det(expm(m33)) / exp(trace(m33)), 1) && close(expm(m33 * 2)[0, 1] / (expm(m33) * expm(m33))[0, 1], 1), "det(expm(m)) = exp(trace(m)) and expm(2m) = expm(m)^2");
    check(close(expm(p1.sub[0 .. 2, 0 .. 2])[0, 0], expm(p1)[0, 0]) && close(p1.expm[1, 1], expm(p1)[1, 1]), "expm of a view, as a method");
    check(collectException!GretlException(expm(a)) !is null && expm(make(0, 0)).nrow == 0, "expm of a non-square matrix throws, of a 0 x 0 matrix is empty");

    /* rank and eye */
    check(rank(p1) == 2 && rank(rd2) == 1 && rank(make(3, 3) * 0) == 0 && rank(make(0, 3)) == 0, "rank of full, rank-1, zero and empty matrices");
    check(rank(make(5, 3)) == 2 && rank(make(3, 5)) == 2 && p1.rank == 2, "rank of the (row + column) test matrices, as a method");
    auto lowrank = make(6, 4);                          // columns 0 and 1 are independent, 2 = 0 + 1, 3 = 2 * 0
    foreach (i; 0 .. 6) { lowrank[i, 0] = 1.0 + i; lowrank[i, 1] = (i * i) % 5 - 1.5; lowrank[i, 2] = lowrank[i, 0] + lowrank[i, 1]; lowrank[i, 3] = 2 * lowrank[i, 0]; }
    check(rank(lowrank) == 2 && rank(t(lowrank)) == 2, "rank of a 6 x 4 matrix of rank 2 and its transpose");
    auto elong = make(9, 2);                            // the library's rank gives 1 here (bugs.md #111); the singular values are 4.2 and 2e-9
    foreach (i; 0 .. 9) { elong[i, 0] = 1; elong[i, 1] = 1 + (i % 2 ? 1e-9 : -1e-9); }
    check(rank(elong) == 2 && rank(elong, 1e-6) == 1, "rank of an elongated nearly rank-1 matrix, and with a larger tol");
    check(rank(p1.sub[0 .. 2, 0 .. 2]) == 2 && rank(m5) == 5, "rank of a view and of the 5 x 5 matrix");
    auto id3 = eye(3);
    bool ide = id3.nrow == 3 && id3.ncol == 3;
    foreach (i; 0 .. 3) foreach (j; 0 .. 3) ide = ide && id3[i, j] == (i == j ? 1.0 : 0.0);
    check(ide, "eye(3) is the identity");
    auto e23 = eye(2, 4);
    check(e23.nrow == 2 && e23.ncol == 4 && e23[0, 0] == 1 && e23[1, 1] == 1 && e23[1, 2] == 0 && e23[0, 3] == 0 && eye(4, 2)[1, 1] == 1 && eye(4, 2)[3, 1] == 0, "eye(nr, nc)");
    check(eye(0).nrow == 0 && eye(0, 3).ncol == 3 && close(det(eye(5)), 1) && rank(eye(4)) == 4, "eye of size 0, det and rank of eye");
    auto emul = p1 * eye(2), emul2 = eye(2) * p1;
    check(emul[1, 0] == 3 && emul[0, 1] == 2 && emul2[1, 1] == 4 && (eye(2) * 3)[1, 1] == 3, "eye is the identity for the product");

    /* norm */
    check(close(norm(p1), sqrt(30.0)) && close(p1.norm(Norm.frobenius), sqrt(30.0)), "frobenius norm of [1 2; 3 4]");
    check(norm(p1, Norm.one) == 6 && norm(p1, Norm.infinity) == 7 && norm(p1, Norm.max) == 4, "one, infinity and max norms of [1 2; 3 4]");
    check(close(norm(p1, Norm.spectral), 5.464985704219043) && close(norm(p1, Norm.spectral), svd(p1).s[0, 0]), "spectral norm is the largest singular value");
    auto vcol = make(3, 1);
    vcol[0, 0] = 3; vcol[1, 0] = -4; vcol[2, 0] = 12;
    check(close(norm(vcol), 13) && close(norm(vcol, Norm.spectral), 13) && norm(vcol, Norm.one) == 19 && norm(vcol, Norm.infinity) == 12 && norm(t(vcol), Norm.infinity) == 19, "norms of a vector");
    auto tiny = vcol * 1e-200, huge = vcol * 1e200;
    check(close(norm(tiny) / 1e-200, 13) && close(norm(huge) / 1e200, 13), "frobenius does not underflow or overflow");
    check(norm(p1.sub[0 .. 2, 0 .. 1], Norm.one) == 4 && norm(p1.sub[0 .. 2, 0 .. 1], Norm.infinity) == 3 && close(norm(p1.sub[1 .. 2, 0 .. 2]), 5), "norms of views");
    check(norm(make(0, 3)) == 0 && norm(make(2, 2) * 0) == 0 && norm(make(2, 2) * 0, Norm.spectral) == 0, "norm of an empty and a zero matrix");
    auto withNan = make(2, 2);
    withNan[0, 1] = double.nan;
    check(isNaN(norm(withNan)) && isNaN(norm(withNan, Norm.one)) && isNaN(norm(withNan, Norm.infinity)) && isNaN(norm(withNan, Norm.max)) && isNaN(norm(withNan, Norm.spectral)), "norm of a matrix with a NaN is NaN");
    check(norm(p1 * p2, Norm.spectral) <= norm(p1, Norm.spectral) * norm(p2, Norm.spectral) + 1e-9 && norm(p1 * p2, Norm.one) <= norm(p1, Norm.one) * norm(p2, Norm.one), "the norms are submultiplicative");

    /* cond */
    auto sp1 = svd(p1).s;
    check(close(cond(p1), sp1[0, 0] / sp1[1, 0]) && close(p1.cond, 14.933034373659268), "cond of [1 2; 3 4] (2-norm)");
    check(close(cond(p1, Norm.one), 6 * 3.5) && close(cond(p1, Norm.infinity), 7 * 3.0) && close(cond(p1, Norm.frobenius), sqrt(30.0) * sqrt(7.5)), "cond in the one, infinity and Frobenius norms");
    check(cond(eye(4)) == 1 && close(cond(edg), 2.0 / 1) && cond(eye(3), Norm.one) == 1, "cond of the identity, of diag(1, -2)");
    check(cond(rd2) == double.infinity || cond(rd2) > 1e15, "cond of a singular matrix (2-norm) is infinite or huge");
    check(cond(sing, Norm.one) == double.infinity && cond(sing, Norm.infinity) == double.infinity && cond(sing, Norm.frobenius) == double.infinity, "cond of a singular matrix in the other norms is infinity");
    check(cond(make(0, 0) , Norm.one) == 0 && cond(make(0, 3)) == 0, "cond of an empty matrix is 0");
    auto condTall = make(3, 2);
    condTall[0, 0] = 1; condTall[0, 1] = 0; condTall[1, 0] = 0; condTall[1, 1] = 4; condTall[2, 0] = 0; condTall[2, 1] = 0;
    check(close(cond(condTall), 4) && close(cond(t(condTall)), 4), "cond of a non-square matrix is the ratio of the singular values");
    check(close(cond(p1.sub[0 .. 2, 0 .. 2]), cond(p1)) && close(cond(p1.sub[0 .. 2, 0 .. 2], Norm.one), cond(p1, Norm.one)), "cond of a view");
    check(collectException!GretlException(cond(condTall, Norm.one)) !is null && collectException!GretlException(cond(p1, Norm.max)) !is null, "cond: a non-square matrix with the one norm, and Norm.max, throw");
    auto hilb = make(8, 8);
    foreach (i; 0 .. 8) foreach (j; 0 .. 8) hilb[i, j] = 1.0 / (i + j + 1);
    check(cond(hilb) > 1e9 && cond(hilb) < 1e11, "cond of the 8 x 8 Hilbert matrix is about 1.5e10");

    /* sqrtm */
    auto dd = make(2, 2);
    dd[0, 0] = 4; dd[0, 1] = 0; dd[1, 0] = 0; dd[1, 1] = 9;
    auto rtD = sqrtm(dd);
    check(close(rtD[0, 0], 2) && close(rtD[1, 1], 3) && close(rtD[0, 1], 0) && close(rtD[1, 0], 0), "sqrtm of diag(4, 9)");
    auto rtS = sqrtm(spd), rtSS = rtS * rtS;                        // spd = [4 2; 2 3]
    check(close(rtSS[0, 0], 4) && close(rtSS[0, 1], 2) && close(rtSS[1, 0], 2) && close(rtSS[1, 1], 3) && close(rtS[0, 1], rtS[1, 0]), "sqrtm of a symmetric positive definite matrix squares back and is symmetric");
    auto ones = make(2, 2);
    ones[0, 0] = 1; ones[0, 1] = 1; ones[1, 0] = 1; ones[1, 1] = 1;
    auto rtO = sqrtm(ones);
    check(close(rtO[0, 0], 0.7071067811865476) && close(rtO[0, 1], 0.7071067811865476) && close(rtO[1, 0], 0.7071067811865476) && close(rtO[1, 1], 0.7071067811865476), "sqrtm of a singular positive semidefinite matrix");
    check(sqrtm(make(2, 2) * 0)[1, 1] == 0 && close(sqrtm(eye(3))[2, 2], 1) && close(sqrtm(eye(3))[0, 2], 0), "sqrtm of the zero matrix and the identity");
    auto ns3 = make(3, 3);                                      // not symmetric, positive eigenvalues
    foreach (i; 0 .. 3) foreach (j; 0 .. 3) ns3[i, j] = [4.0, 0, 1, 1, 9, 0, 0, 2, 16][j * 3 + i];
    auto rtN = sqrtm(ns3), rtNN = rtN * rtN;
    bool nsq = true;
    foreach (i; 0 .. 3) foreach (j; 0 .. 3) nsq = nsq && close(rtNN[i, j], ns3[i, j]);
    check(nsq, "sqrtm of a non-symmetric matrix squares back (Denman-Beavers)");
    auto upper2 = make(2, 2);
    upper2[0, 0] = 4; upper2[0, 1] = 1; upper2[1, 0] = 0; upper2[1, 1] = 9;
    auto rtU = sqrtm(upper2);
    check(close(rtU[0, 0], 2) && close(rtU[1, 1], 3) && close(rtU[1, 0], 0) && close(rtU[0, 1], 0.2), "sqrtm of an upper triangular matrix is upper triangular (the (0,1) entry is 1 / (2 + 3))");
    check(collectException!GretlException(sqrtm(p1)) !is null, "no real root: [1 2; 3 4] has a negative eigenvalue");
    auto negdef = make(2, 2);
    negdef[0, 0] = -1; negdef[0, 1] = 0; negdef[1, 0] = 0; negdef[1, 1] = -2;
    check(collectException!GretlException(sqrtm(negdef)) !is null, "a negative definite symmetric matrix throws");
    check(close(sqrtm(spd.sub[0 .. 2, 0 .. 2])[0, 0], rtS[0, 0]) && close(spd.sqrtm[1, 1], rtS[1, 1]), "sqrtm of a view, as a method");
    check(collectException!GretlException(sqrtm(a)) !is null && sqrtm(make(0, 0)).nrow == 0, "sqrtm of a non-square matrix throws, of a 0 x 0 matrix is empty");

    /* logm */
    auto lg0 = logm(eye(3));
    bool lz = true;
    foreach (i; 0 .. 3) foreach (j; 0 .. 3) lz = lz && abs(lg0[i, j]) < 1e-15;
    check(lz, "logm of the identity is zero");
    auto lde = make(2, 2);
    lde[0, 0] = exp(1.0); lde[0, 1] = 0; lde[1, 0] = 0; lde[1, 1] = exp(2.0);
    check(close(logm(lde)[0, 0], 1) && close(logm(lde)[1, 1], 2) && close(logm(lde)[0, 1], 0), "logm of diag(e, e^2)");
    auto lgS = logm(spd), eS = expm(lgS);
    check(close(eS[0, 0], 4) && close(eS[0, 1], 2) && close(eS[1, 1], 3) && close(lgS[0, 1], lgS[1, 0]), "expm(logm(m)) = m for a symmetric positive definite matrix");
    auto lgN = logm(ns3), eN = expm(lgN);
    bool lnok = true;
    foreach (i; 0 .. 3) foreach (j; 0 .. 3) lnok = lnok && abs(eN[i, j] - ns3[i, j]) < 1e-8;
    check(lnok, "expm(logm(m)) = m for a non-symmetric matrix");
    auto lgU = logm(upper2);                                  // [4 1; 0 9]: log 4, log 9 on the diagonal, (log 9 - log 4) / 5 above it
    check(close(lgU[0, 0], log(4.0)) && close(lgU[1, 1], log(9.0)) && close(lgU[1, 0], 0) && close(lgU[0, 1], (log(9.0) - log(4.0)) / 5), "logm of an upper triangular matrix");
    auto small = p1 * 0.1, lsmall = logm(expm(small));
    check(close(lsmall[0, 0], small[0, 0]) && close(lsmall[0, 1], small[0, 1]) && close(lsmall[1, 0], small[1, 0]) && close(lsmall[1, 1], small[1, 1]), "logm(expm(m)) = m for a small matrix");
    auto lrot = logm(expm(th));                               // a rotation by 0.7: complex eigenvalues, real logarithm
    check(close(lrot[0, 1], -0.7) && close(lrot[1, 0], 0.7) && abs(lrot[0, 0]) < 1e-10, "logm of a rotation is the skew matrix");
    auto lbig = make(2, 2);
    lbig[0, 0] = 1000; lbig[0, 1] = 3; lbig[1, 0] = 0; lbig[1, 1] = 0.001;
    check(close(logm(lbig)[0, 0], log(1000.0)) && close(logm(lbig)[1, 1], log(0.001)), "logm of a matrix far from the identity");
    check(collectException!GretlException(logm(p1)) !is null, "logm of a matrix with a negative eigenvalue throws");
    check(collectException!GretlException(logm(rd2)) !is null && collectException!GretlException(logm(negdef)) !is null, "logm of a singular or negative definite symmetric matrix throws");
    check(close(logm(spd.sub[0 .. 2, 0 .. 2])[0, 0], lgS[0, 0]) && close(spd.logm[1, 1], lgS[1, 1]), "logm of a view, as a method");
    check(collectException!GretlException(logm(a)) !is null && logm(make(0, 0)).nrow == 0, "logm of a non-square matrix throws, of a 0 x 0 matrix is empty");

    /* sort and order */
    auto so = make(5, 3);          // columns: key (with ties), second key, a label
    double[5] key = [2, 1, 2, 1, 3], key2 = [9, 8, 7, 8, 1];
    foreach (i; 0 .. 5) { so[i, 0] = key[i]; so[i, 1] = key2[i]; so[i, 2] = 100 + i; }
    check(order(so, 0) == [cast(size_t) 1, 3, 0, 2, 4], "order by the first column, ties keep their original order (stable)");
    check(so.order(0, 1) == [cast(size_t) 1, 3, 2, 0, 4], "order by the first column, then the second");
    check(so.order(1) == [cast(size_t) 4, 2, 1, 3, 0], "order by the second column");
    auto sorted = sort(so, 0, 1);
    bool sortedOk = sorted.nrow == 5 && sorted.ncol == 3;
    foreach (rw; 0 .. 5) sortedOk = sortedOk && sorted[rw, 2] == 100 + [1, 3, 2, 0, 4][rw] && sorted[rw, 0] == key[[1, 3, 2, 0, 4][rw]];
    check(sortedOk, "sort moves whole rows");
    check(so[0, 0] == 2 && so[4, 0] == 3 && so[2, 2] == 102, "sort leaves the matrix alone");
    check(so.sort(2)[4, 2] == 104 && so.sort(2)[0, 2] == 100, "sort by a column that is already inOrder");
    auto tied = make(4, 2);
    foreach (i; 0 .. 4) { tied[i, 0] = 1; tied[i, 1] = 10 - i; }
    check(order(tied, 0) == [cast(size_t) 0, 1, 2, 3] && order(tied, 0, 1) == [cast(size_t) 3, 2, 1, 0], "all rows tied: unchanged; then broken by the second column");
    auto nanm = make(4, 1);
    nanm[0, 0] = 3; nanm[1, 0] = double.nan; nanm[2, 0] = 1; nanm[3, 0] = double.nan;
    check(order(nanm, 0) == [cast(size_t) 2, 0, 1, 3] && isNaN(sort(nanm, 0)[3, 0]) && sort(nanm, 0)[0, 0] == 1, "NaN sorts last");
    auto vsub = so.sub[1 .. 5, 0 .. 2];                             // rows 1..4, keys 1 2 1 3 / 8 7 8 1
    check(order(vsub, 0) == [cast(size_t) 0, 2, 1, 3] && sort(vsub, 0, 1)[1, 1] == 8 && vsub.sort(1)[0, 1] == 1, "order and sort of a view");
    check(order(so.cols[0 .. 1], 0).length == 5 && so.rows[0 .. 0].order(0).length == 0 && sort(make(0, 3), 1).nrow == 0, "single column and empty cases");
    check(collectException!GretlException(order(so, 3)) !is null && collectException!GretlException(sort(so, 0, 5)) !is null, "a column out of range throws");
    check(!__traits(compiles, order(so)), "at least one column is required");
    auto big1 = make(200, 2);
    foreach (i; 0 .. 200) { big1[i, 0] = (i * 37) % 11; big1[i, 1] = i; }
    auto bigSorted = sort(big1, 0, 1);
    bool inOrder = true;
    foreach (rw; 1 .. 200) inOrder = inOrder && (bigSorted[rw - 1, 0] < bigSorted[rw, 0] || (bigSorted[rw - 1, 0] == bigSorted[rw, 0] && bigSorted[rw - 1, 1] < bigSorted[rw, 1]));
    check(inOrder, "sort of 200 rows is in lexicographic order");

    /* cumsum and cumprod */
    auto cumS = cumsum(a);                                  // a = [0 1 2; 10 11 12]
    check(cumS.nrow == 2 && cumS.ncol == 3 && cumS[0, 0] == 0 && cumS[1, 0] == 10 && cumS[1, 1] == 12 && cumS[1, 2] == 14 && cumS[0, 2] == 2, "cumsum down the columns");
    auto cumSR = cumsum(a, Along.row);
    check(cumSR[0, 0] == 0 && cumSR[0, 1] == 1 && cumSR[0, 2] == 3 && cumSR[1, 0] == 10 && cumSR[1, 1] == 21 && cumSR[1, 2] == 33, "cumsum across the rows");
    auto cumM = make(3, 3);
    foreach (i; 0 .. 3) foreach (j; 0 .. 3) cumM[i, j] = 1 + i + j;                 // [1 2 3; 2 3 4; 3 4 5]
    auto cumP = cumprod(cumM), cumPR = cumprod(cumM, Along.row);
    check(cumP[0, 0] == 1 && cumP[1, 0] == 2 && cumP[2, 0] == 6 && cumP[2, 1] == 2 * 3 * 4 && cumP[2, 2] == 3 * 4 * 5, "cumprod down the columns");
    check(cumPR[0, 2] == 6 && cumPR[1, 1] == 6 && cumPR[1, 2] == 24 && cumPR[2, 2] == 60 && cumPR[2, 0] == 3, "cumprod across the rows");
    auto vec5 = make(5, 1);
    foreach (i; 0 .. 5) vec5[i, 0] = i + 1;
    check(cumsum(vec5)[4, 0] == 15 && cumprod(vec5)[4, 0] == 120 && cumsum(t(vec5), Along.row)[0, 4] == 15 && cumsum(t(vec5))[0, 4] == 5, "vectors in both directions");
    check(cumsum(a)[1, 2] == a.cols[2 .. 3][0, 0] + a.cols[2 .. 3][1, 0] && cumsum(a.sub[0 .. 2, 1 .. 3], Along.row)[1, 1] == 11 + 12 && cumprod(a.sub[1 .. 2, 0 .. 3])[0, 2] == 12 && a.cumsum[1, 0] == 10, "cumsum and cumprod of views, as methods");
    auto cumN = make(3, 2);
    foreach (i; 0 .. 3) foreach (j; 0 .. 2) cumN[i, j] = 1;
    cumN[1, 0] = double.nan;
    check(cumsum(cumN)[0, 0] == 1 && isNaN(cumsum(cumN)[1, 0]) && isNaN(cumsum(cumN)[2, 0]) && cumsum(cumN)[2, 1] == 3 && isNaN(cumprod(cumN, Along.row)[1, 1]), "a NaN spreads along its column (or row)");
    auto csBefore = a[1, 2];
    cumsum(a);
    check(a[1, 2] == csBefore && cumsum(make(0, 3)).nrow == 0 && cumprod(make(3, 0)).ncol == 0, "the original is unchanged, empty matrices");
    auto lastCol = cumsum(big1.cols[1 .. 2]);
    check(lastCol[199, 0] == 199 * 200 / 2, "cumsum of a long column");

    /* cbind, rbind, ~ and ~= */
    auto cb = cbind(p1, p2);                                // [1 2 5 6; 3 4 7 8]
    check(cb.nrow == 2 && cb.ncol == 4 && cb[0, 0] == 1 && cb[0, 2] == 5 && cb[1, 3] == 8 && cb[1, 1] == 4, "cbind of two matrices");
    auto cb3 = cbind(p1, p2, p1);
    check(cb3.ncol == 6 && cb3[1, 5] == 4 && cb3[0, 4] == 1, "cbind of three");
    auto cbs = cbind(1, p1);
    check(cbs.ncol == 3 && cbs[0, 0] == 1 && cbs[1, 0] == 1 && cbs[0, 1] == 1 && cbs[1, 2] == 4, "cbind(1, m) adds a column of ones");
    auto cbv = cbind(p1, [10.0, 20.0], 7);
    check(cbv.ncol == 4 && cbv[0, 2] == 10 && cbv[1, 2] == 20 && cbv[0, 3] == 7 && cbv[1, 3] == 7, "cbind with an array and a number");
    const(double)[] cser = [1.5, 2.5];
    check(cbind(cser, cser).ncol == 2 && cbind(cser, [3, 4])[1, 1] == 4, "cbind of const(double)[] series and an int array");
    check(cbind(1, 2, 3).nrow == 1 && cbind(1, 2, 3).ncol == 3 && cbind([1.0, 2, 3], [4.0, 5, 6]).nrow == 3, "cbind of numbers is one row; of vectors, columns");
    auto cbview = cbind(a.cols[0 .. 1], a.sub[0 .. 2, 2 .. 3]);
    check(cbview.ncol == 2 && cbview[1, 0] == 10 && cbview[1, 1] == 12, "cbind of views");
    auto rb = rbind(p1, p2);                                // [1 2; 3 4; 5 6; 7 8]
    check(rb.nrow == 4 && rb.ncol == 2 && rb[0, 0] == 1 && rb[2, 0] == 5 && rb[3, 1] == 8 && rb[1, 1] == 4, "rbind of two matrices");
    auto rbv = rbind(p1, [9.0, 10.0], 0);
    check(rbv.nrow == 4 && rbv[2, 0] == 9 && rbv[2, 1] == 10 && rbv[3, 0] == 0 && rbv[3, 1] == 0, "rbind with an array (a row) and a number");
    check(rbind(1, 2, 3).ncol == 1 && rbind(1, 2, 3).nrow == 3 && rbind([1, 2], [3, 4]).nrow == 2 && rbind([1, 2], [3, 4])[1, 0] == 3, "rbind of numbers is one column; of arrays, rows");
    check(rbind(a.rows[0 .. 1], a.rows[1 .. 2])[1, 2] == 12, "rbind of views");
    GretlMatrix acc;                                        // an empty start, grown in a loop
    foreach (i; 0 .. 3) acc = rbind(acc, [cast(double) i, 10.0 * i]);
    check(acc.nrow == 3 && acc.ncol == 2 && acc[2, 1] == 20 && acc[1, 0] == 1, "rbind onto an empty GretlMatrix.init");
    check(rbind(GretlMatrix.init, GretlMatrix.init).nrow == 0 && cbind(GretlMatrix.init, p1).ncol == 2, "0 x 0 matrices are left out");
    check(collectException!GretlException(cbind(p1, big)) !is null && collectException!GretlException(rbind(p1, make(2, 3))) !is null, "cbind and rbind with sizes that do not fit throw");
    check(collectException!GretlException(cbind(p1, [1.0, 2, 3])) !is null && collectException!GretlException(rbind(p1, [1.0])) !is null, "an array of the wrong length throws");
    auto tl = p1 ~ p2;
    check(tl.nrow == 2 && tl.ncol == 4 && tl[0, 2] == 5, "m ~ mm");
    auto tl2 = 1 ~ p1 ~ [9.0, 8.0];
    check(tl2.ncol == 4 && tl2[0, 0] == 1 && tl2[1, 3] == 8 && tl2[1, 2] == 4, "a ~ m ~ array");
    check((a.sub[0 .. 2, 0 .. 2] ~ b.cols[0 .. 1]).ncol == 3 && ([5.0, 6.0] ~ p1)[1, 0] == 6, "views with ~, an array on the left");
    auto grow = make(2, 1);
    grow[0, 0] = 1; grow[1, 0] = 2;
    grow ~= [3.0, 4.0];
    grow ~= 9;
    grow ~= p1;
    check(grow.nrow == 2 && grow.ncol == 5 && grow[0, 1] == 3 && grow[1, 1] == 4 && grow[1, 2] == 9 && grow[0, 3] == 1 && grow[1, 4] == 4 && grow[1, 0] == 2, "m ~= adds columns");
    auto alias1 = grow, view1 = grow.sub[0 .. 2, 0 .. 2];
    grow ~= [100.0, 200.0];
    check(grow.ncol == 6 && grow[1, 5] == 200 && alias1.ncol == 5 && view1[1, 1] == 4 && alias1[1, 4] == 4, "m ~= leaves other handles and views untouched");
    view1[0, 0] = -1;
    check(alias1[0, 0] == -1 && grow[0, 0] == 1, "after ~=, the old handle and the grown one no longer share data");
    grow ~= grow;
    check(grow.ncol == 12 && grow[1, 11] == 200 && grow[1, 6] == 2, "m ~= m");
    GretlMatrix fresh;
    fresh ~= [1.0, 2.0];
    fresh ~= [3.0, 4.0];
    check(fresh.nrow == 2 && fresh.ncol == 2 && fresh[1, 1] == 4, "~= onto an empty GretlMatrix.init");
    check(collectException!GretlException(grow ~= [1.0, 2.0, 3.0]) !is null && grow.ncol == 12, "~= with the wrong length throws and changes nothing");
    auto many = make(3, 0 + 1);
    foreach (i; 0 .. 3) many[i, 0] = i;
    foreach (k; 0 .. 500) many ~= many.cols[0 .. 1];            // a view of itself: copies each time, still correct
    check(many.ncol == 501 && many[2, 500] == 2, "many appends");
    auto loopm = make(3, 1);
    foreach (k; 0 .. 2000) loopm ~= [1.0, 2.0, 3.0];
    check(loopm.ncol == 2001 && loopm[2, 2000] == 3 && loopm[0, 1000] == 1, "2000 appends in place");
    auto rawBefore = loopm.ptr;
    loopm ~= [7.0, 8.0, 9.0];
    check(loopm.ptr is rawBefore && loopm.ncol == 2002, "~= grows the same matrix when nothing else refers to it");
    auto holder = loopm;
    loopm ~= [7.0, 8.0, 9.0];
    check(loopm.ptr !is rawBefore && holder.ptr is rawBefore && holder.ncol == 2002 && loopm.ncol == 2003, "~= makes a new matrix when another handle exists");

    check(collectException!GretlException(p1 / big) !is null, "m / mm with different column counts throws");
    check(collectException!GretlException(Xt.sub[0 .. 2, 0 .. 1] / Xt.sub[0 .. 2, 0 .. 1]) !is null, "underdetermined m / mm throws");
    check((a * 1.0 + 0).nrow == 2 && (make(0, 3) + make(0, 3)).nrow == 0, "empty matrices");

    /* an empty view is fine */
    check(m.sub[2 .. 2, 0 .. 4].nrow == 0, "empty view");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
