/* Port of gensyntax.c: the parser module for 'genr' and related commands: the recursive-descent parser that builds the expression tree from the
   tokens of genlex.d (`expr`, over the levels expr0 .. expr4, term, factor, powterm and base), the node constructors (new_node, newempty, newdbl,
   newb2, newb3, newbn, bncopy, newref, newstr), and the argument parsers for function calls, strings, bundles, matrix definitions and slices.

   The C ids (F_PRINTF, P_COM, G_RPR, EOT, OP_MAX, P_SLICING, R_..., ...) are the members of the named enums of genparseh.d / genmainh.d, replaced by
   a table derived from them; the anonymous enum of StrArg.RIGHT_STR / StrArg.MID_STR is StrArg; the table of argument records is a __gshared array of ArgRecord.
   The function-like macros are functions. The debug output (SDEBUG / MDEBUG, off in the C) is not in the port. new_node, newempty, newdbl, newb2,
   bncopy and expr are extern (C) and public: other modules declared them extern (C) while this file was missing. */

module gretl.gensyntax;

import gretl.gretl_func1 : get_user_function_by_name;
import gretl.geneval1 : free_tree, getsymb;
import gretl.geneval5 : parser_advance, parser_char_index, parser_getc, parser_print_input, parser_ungetc;
import gretl.genlex : context_error, getsymb_full, lex, function_lookup, obs_node, set_doing_genseries, set_parsing_query;
import gretl.gstrfuncs : g_strchomp, g_strchug;
import gretl.gretl_array : get_array_by_name;
import gretl.gretl_bundle : gretl_is_bundle;
import gretl.gretl_errors : gretl_errmsg_set, gretl_errmsg_sprintf;
import gretl.gretl_prn : pprintf, pputc, pputs;
import gretl.gretl_string_table : get_built_in_string_by_name;
import gretl.strutils : gretl_namechar_spn, gretl_strdup, gretl_strndup, tailstrip;
import gretl.uservar : get_string_by_name, get_user_var_of_type_by_name;
import gretl.uservarh : UserVar;

import core.stdc.math : fabs, floor, isnan, log, sqrt, pow, exp;
import core.stdc.stdlib : free, malloc, calloc, realloc, atoi, atof, strtol, strtod;
import core.stdc.stdio : fprintf, printf, sprintf, snprintf, sscanf, stderr;
import core.stdc.string : strncat, strrchr, memcpy, memset, strcmp, strlen, strcpy, strncpy, strcat, strchr, strstr, strncmp, strspn, strcspn;
import core.stdc.ctype : isalpha, isdigit, isspace, isalnum, isupper, islower, isprint;

import gretl.genparseh;
import gretl.genmainh;
import gretl.libgretlh;
import gretl.gretl_matrixh;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_prn : Prn;



/* Parser module for 'genr' and related commands */




extern (C) Node * new_node(int t)
{
    Node * n = cast(Node *) calloc(1, Node.sizeof);


    if (n != null) {
	n.t = cast(short) t;
	n.vnum = NO_VNUM;
	n.L = n.M = n.R = null;
	n.parent = null;
    }

    return n;
}

extern (C) Node * newempty()
{
    Node *n = new_node(GenSym.empty);

    if (n != null) {
	n.v.idnum = 0;
    }

    return n;
}

private Node * newref(Parser * p, int t)
{
    Node *n = new_node(t);

    if (n != null) {
	if (t == GenSym.series) {
	    n.vnum = p.idnum;
	    n.v.xvec = p.dset.Z[n.vnum];
	    n.vname = p.idstr;
	} else if (t == GenSym.num || t == GenSym.num_p || t == GenSym.num_m) {
	    UserVar * u = cast(UserVar *) p.data;

	    n.vname = p.idstr;
	    n.v.xval = *cast(double *) u.ptr;
	    n.uv = u;
	} else if (t == GenSym.mat || t == GenSym.list || t == GenSym.bundle ||
		   t == GenSym.array || t == GenSym.str) {
	    UserVar * u = cast(UserVar *) p.data;

	    n.vname = p.idstr;
	    n.v.ptr = u.ptr;
	    n.uv = u;
	} else if (t == GenSym.osl) {
	    UserVar * u = cast(UserVar *) p.data;

	    n.vname = p.idstr;
	    n.v.ptr = u.ptr;
	    n.uv = u;
	    n.flags |= NodeFlags.lht;
	} else if (t == GenSym.ufun) {
	    n.vname = p.idstr; /* function name */
	    n.v.ptr = p.data;  /* pointer to function */
        } else if (t == GenSym.rfun) {
            n.vname = p.idstr; /* function name */
	} else if (t == GenSym.dbundle) {
	    n.v.idnum = p.idnum;
	} else if (t == GenSym.mmemb) {
	    /* for case of named model */
	    n.vname = p.idstr;
	    n.v.idnum = -1;
	    n.t = GenSym.dbundle; /* roll into DBUNDLE case */
	} else if (t == GenSym.undef) {
	    n.vname = p.idstr;
	} else if (t == GenSym.uobj || t == GenSym.wlist) {
	    n.v.str = p.idstr;
	} else {
	    n.v.idnum = p.idnum;
	}
    }

    return n;
}

/* node storing an anonymous string value */

private Node * newstr(char * s)
{
    Node *n = null;

    if (s == null) {
	fprintf(stderr, "newstr: input is null\n");
    } else {
	n = new_node(GenSym.str);
	if (n != null) {
	    n.v.str = s;
	    n.flags = NodeFlags.tmp;
	}
    }

    return n;
}

/* node storing an anonymous floating point value */

extern (C) Node * newdbl(double x)
{
    Node *n = new_node(GenSym.num);

    if (n != null) {
	n.v.xval = x;
    }

    return n;
}

/* node for binary operator or two-argument function */

extern (C) Node * newb2(int t, Node * l, Node * r)
{
    Node *n = new_node(t);

    if (n != null) {
	n.L = l;
	n.R = r;
    }

    return n;
}

/* node for ternary operator */

private Node * newb3(int t, Node * l)
{
    Node *n = new_node(t);

    if (n != null) {
	n.L = l;
    }

    return n;
}

/* node for unknown number of subnodes */

private Node * newbn(int t)
{
    Node *n = new_node(t);

    if (n != null) {
	n.v.bn.n_nodes = 0;
	n.v.bn.n = null;
    }

    return n;
}

extern (C) Node * bncopy(Node * t, int * err)
{
    Node *n = new_node(t.t);

    if (n != null) {
	int i;
	int k = t.v.bn.n_nodes;
	Node ** nn = cast(Node **) malloc(k * (Node *).sizeof);

	if (nn == null) {
	    *err = GretlError.alloc;
	} else {
	    n.v.bn.n_nodes = k;
	    n.v.bn.n = nn;
	    for (i=0; i<k; i++) {
		n.v.bn.n[i] = t.v.bn.n[i];
	    }
	}
    }

    return n;
}

private int push_bn_node(Node * t, Node * n)
{
    Node **nn;
    int k = t.v.bn.n_nodes;

    nn = cast(Node* *) realloc(t.v.bn.n, (k + 1) * (*nn).sizeof);
    if (nn == null) {
	return GretlError.alloc;
    }

    t.v.bn.n = nn;
    t.v.bn.n[k] = n;
    t.v.bn.n_nodes += 1;
    n.parent = t;


    return 0;
}

private void expected_symbol_error(int c, Parser * p, int badc)
{
    parser_print_input(p);

    if (badc) {
	gretl_errmsg_sprintf(_("Expected '%c' but found '%c'\n").ptr, c,
			     badc);
    } else {
	const char *found = getsymb_full(p.sym, p);

	if (found == null || *found == '\0') {
	    gretl_errmsg_sprintf(_("Expected '%c' but formula ended\n").ptr, c);
	} else {
	    gretl_errmsg_sprintf(_("Expected '%c' but found '%s'\n").ptr, c,
				 found);
	}
    }

    p.err = GretlError.parse;
}

private void unmatched_symbol_error(int c, Parser * p)
{
    if (c == '(' || c == ')') {
	fprintf(stderr, "*** gensyntax: unmatched_symbol_error '%c'\n", c);
    } else {
	fprintf(stderr, "*** gensyntax: unmatched_symbol_error (%c)\n", c);
    }
    parser_print_input(p);
    pprintf(p.prn, _("Unmatched '%c'\n").ptr, c);
    p.err = GretlError.parse;
}

private Node * base(Parser * p, Node * up)
{
    Node *t = null;

    if (p.err) {
	return null;
    }


    switch (p.sym) {
    case GenSym.cnum:
	t = newdbl(p.xval);
	lex(p);
	break;
    case GenSym.cstr:
	t = newstr(p.idstr);
	lex(p);
	break;
    case GenSym.num:
    case GenSym.series:
    case GenSym.mat:
    case GenSym.bundle:
    case GenSym.array:
    case GenSym.str:
    case GenSym.uobj:
    case GenSym.con:
    case GenSym.dvar:
    case GenSym.mvar:
    case GenSym.list:
    case GenSym.wlist:
    case GenSym.dbundle:
    case GenSym.num_p:
    case GenSym.num_m:
    case GenSym.osl:
    case GenSym.undef:
	t = newref(p, p.sym);
	lex(p);
	break;
    case GenSym.dum:
	if (p.idnum == DumId.null_ || p.idnum == DumId.empty) {
	    t = newempty();
	} else {
	    t = newref(p, p.sym);
	}
	lex(p);
	break;
    case GenSym.empty:
    case GenSym.nullarg:
	t = newempty();
	if (p.sym == GenSym.nullarg) {
	    t.vname = p.idstr;
	}
	lex(p);
	break;
    case GenSym.b_range:
    case GenSym.p_dot:
	lex(p);
	t = expr(p);
	break;
    case GenSym.b_join:
	/* list joiner with empty LHS */
	t = newempty();
	break;
    case GenSym.g_lpr: /* left paren '(' */
	lex(p);
	t = expr(p);
	if (p.sym == GenSym.g_rpr) {
	    lex(p);
	} else if (p.err == 0) {
	    expected_symbol_error(')', p, 0);
	}
	break;
    case GenSym.g_lbr: /* left bracket '[' */
	if (up == null) {
	    goto deferr;
	}
	if (up.t == GenSym.obs) {
	    t = obs_node(p);
	}
	if (p.sym == GenSym.g_rbr) {
	    lex(p);
	} else if (p.err == 0) {
	    expected_symbol_error(']', p, 0);
	}
	break;
    default:
    deferr:
	context_error(0, p, "base");
	break;
    }


    return t;
}

/* Special for the unary '&' operator: the operand must
   be a named series or "user_var", with an additional
   restriction on the type.
*/

private Node * u_addr_base(Parser * p)
{
    Node *t = base(p, null);

    if (t != null && !p.err) {
	switch (t.t) {
	case GenSym.series:
	    if (t.v.xvec == null) {
		p.err = GretlError.types;
	    }
	    break;
	case GenSym.num:
	case GenSym.mat:
	case GenSym.bundle:
	case GenSym.array:
	case GenSym.str:
	case GenSym.list:
	    if (t.uv == null) {
		p.err = GretlError.types;
	    }
	    break;
	case GenSym.osl:
	case GenSym.undef:
	    break;
	case GenSym.empty:
	    if (t.vname == null) {
		p.err = GretlError.types;
	    }
	    break;
	default:
	    p.err = GretlError.types;
	}

	if (p.err) {
	    pputs(p.prn, _("Wrong type of operand for unary '&'").ptr);
	    pputc(p.prn, '\n');
	}
    }

    return t;
}

private void unwrap_string_arg(Parser * p)
{
    int n = cast(int) strlen(p.idstr);

    if (p.idstr[n-1] == '"') {
	p.idstr[n-1] = '\0';
    } else {
	unmatched_symbol_error('"', p);
    }
}

bool varargs_func(int f)
{
    return f == GenFunc.f_printf || f == GenFunc.f_sprintf || f == GenFunc.f_sscanf;
}

bool testc(int c) { return c == '(' || c == '.' || c == '['; }

/* Grab a string argument. Note: we have a mechanism in genlex.c for
   retrieving arguments that take the form of quoted string literals
   or names of string variables.  The special use of this function is
   to grab a literal string without requiring the user to wrap it in
   quotes; we use it only where we know the only acceptable argument
   is a string.  This function assumes that the argument in question
   is the last (or only) argument to a function: we look for a closing
   parenthesis and flag an error if we don't find one.

   Depending on the context we may or may not want to "consume" the
   trailing right paren that follows the string of interest; that's
   controlled by the second argument.

   Added 2013-08-25: The above is all very well, but it breaks the
   case where a function that returns a string is given as an
   argument, rather than a plain string variable or string literal --
   the function call was getting passed as a string literal. A block
   is now inserted to test for this case.
*/

private Node * get_final_string_arg(Parser * p, Node * t, int sym, int eat_last)
{
    const(char) * src = null;
    int n;
    int wrapped = 0;
    int strvar = 0;

    while (p.ch == ' ') {
	parser_getc(p);
    }

    if (!varargs_func(sym) && sym != GenFunc.f_exists) {
	/* Check for a nested function call (2013-08-25) or
	   bundle/array member (2015-09-25). Further fixes
	   applied July 2018.0
	*/
	src = p.point - 1;
	n = gretl_namechar_spn(src);
	if (n > 0 && n < VNAMELEN && testc(src[n])) {
	    Node *ret = null;
	    char[VNAMELEN] tmp;
	    char c = src[n];

	    tmp[0] = '\0';
	    strncat(tmp.ptr, src, n);
	    src = null;
	    if (c == '(') {
		if (function_lookup(tmp.ptr) || get_user_function_by_name(tmp.ptr)) {
		    /* watch out: this is VERY finicky! */
		    /* was: ret = base(p, t); */
		    lex(p);
		    ret = expr(p);
		}
	    } else if (gretl_is_bundle(tmp.ptr)) {
		lex(p);
		ret = expr(p);
	    } else if (c == '[' && get_array_by_name(tmp.ptr)) {
		lex(p);
		ret = expr(p);
	    }
	    if (ret != null) {
		if (eat_last) {
		    /* consume trailing right paren */
		    lex(p);
		}
		return ret;
	    }
	}
    }

    if (p.ch == '"') {
	wrapped = 1;
	if (!varargs_func(sym)) {
	    parser_getc(p);
	}
    }

    if (p.ch == ')') {
	if (wrapped) {
	    p.err = GretlError.parse;
	    return null;
	}
	/* handle empty arg */
	p.idstr = gretl_strdup("");
    } else {
	char p0 = '(', p1 = ')';
	int i;
	int paren = 1;
	int quoted = wrapped;
	int close = -1;
	int started = 0;
	const(char) * s = p.point;

	if (p.ch == p0) paren++;

	/* find length of string to closing paren */
	i = 0;
	while (*s) {
	    if (!quoted) {
		if (*s == p0) {
		    paren++;
		} else if (*s == p1) {
		    paren--;
		}
	    }
	    if (paren == 0) {
		close = i;
		break;
	    }
	    if (*s == '"') {
		quoted = !quoted;
	    }
	    s++;
	    i++;
	}

	if (close < 0) {
	    unmatched_symbol_error(p0, p);
	    return null;
	}

	for (i=0; i<=close; i++) {
	    if (!started && !isspace(p.ch)) {
		src = p.point - 1;
		p.idstr = gretl_strndup(src, close - i + 1);
		started = 1;
	    }
	    parser_getc(p);
	}
    }

    if (p.idstr == null) {
	p.err = GretlError.alloc;
	return null;
    }

    tailstrip(p.idstr);
    if (eat_last) {
	/* consume trailing right paren */
	parser_getc(p);
    }
    lex(p);

    if (!varargs_func(sym)) {
	if (wrapped) {
	    unwrap_string_arg(p);
	} else if (sym != GenFunc.f_isdiscr &&
		   sym != GenFunc.f_typeof &&
		   sym != GenFunc.f_exists) {
	    /* not quoted: give priority to string variables
	       unless we need the _names_ of string variables
	       rather then their content
	    */
	    p.data = get_user_var_of_type_by_name(p.idstr,
						   GretlType.string_);
	    if (p.data != null) {
		strvar = 1;
	    } else {
		char *s = get_built_in_string_by_name(p.idstr);

		if (s != null) {
		    free(p.idstr);
		    p.idstr = gretl_strdup(s);
		}
	    }
	}
    }


    if (p.err) {
	return null;
    } else if (strvar == 1) {
	return newref(p, GenSym.str);
    } else {
	return newstr(p.idstr);
    }
}

/* Here we're allowing for the possibility that an "fncall"
   argument is given in the form of a string variable. If the
   relevant portion of the statement contains any quotation
   marks that cannot be the case.
*/

private void try_for_string_var(Parser * p, int close, int anyquote)
{
    int err = 0;

    if (!anyquote) {
	char *sname = gretl_strndup(p.point - 1, close + 1);
	char *sval = get_string_by_name(sname);

	if (sval != null) {
	    p.idstr = gretl_strdup(sval);
	} else {
	    err = 1;
	}
	free(sname);
    } else {
	err = 1;
    }

    if (err) {
	unmatched_symbol_error(')', p);
    }
}

enum StrArg {
    RIGHT_STR  = 1 << 0,
    MID_STR    = 1 << 1
}

/* get_literal_string_arg() is mostly used to retrieve a string
   that defines a function call, e.g. in BFGSmax(). It can also
   be used to trap a non-quoted string that functions as an
   "enum" value, as in the "vcat" flag for mpireduce().
   In the former case @opt will be 0, in the latter StrArg.MID_STR.
*/

private Node * get_literal_string_arg(Parser * p, int opt)
{
    int close = 0;

    while (p.ch == ' ') {
	parser_getc(p);
    }

    if (p.ch == '"') {
	/* handle the case where the argument is quoted */
	close = parser_char_index(p, '"');
	if (close < 0) {
	    unmatched_symbol_error('"', p);
	    return null;
	}
	p.idstr = gretl_strndup(p.point, close);
	close++;
    } else if (opt == StrArg.MID_STR) {
	/* handle the "enum"-type unquoted case */
	int i1 = parser_char_index(p, ',');
	int i2 = parser_char_index(p, ')');

	if (i2 < 0) {
	    unmatched_symbol_error('(', p);
	    return null;
	}
	close = (i1 >= 0 && i1 < i2)? i1 : i2;
	p.idstr = gretl_strndup(p.point - 1, close + 1);
    } else {
	/* handle the function-call case: the terminator of
	   the argument will be a bare comma (unquoted, not
	   in parentheses) or a right paren which matches
	   the paren that opens the argument list of
	   the function call.
	*/
	const(char) * s = p.point;
	int gotparen = 0;
	int quoted = 0;
	int anyquote = 0;
	int paren = 0;
	int i = 0;

	if (p.ch == '(') {
	    /* an "anonymous function"? */
	    gotparen = paren = 1;
	}

	while (*s) {
	    if (*s == '"') {
		anyquote = 1;
		quoted = !quoted;
	    }
	    if (!quoted) {
		if (*s == '(') {
		    gotparen++;
		    paren++;
		} else if (*s == ')') {
		    if (paren == 0) {
			paren--;
			close = i;
			break;
		    }
		    paren--;
		}
		if (paren == 0) {
		    if (gotparen) {
			/* include right paren */
			close = i + 1;
			break;
		    } else if (*s == ',') {
			/* don't include comma */
			close = i;
			break;
		    }
		}
	    }
	    s++;
	    i++;
	}

	if (paren > 0) {
	    unmatched_symbol_error('(', p);
	} else if (paren < 0) {
	    try_for_string_var(p, close, anyquote);
	} else if (quoted) {
	    unmatched_symbol_error('"', p);
	} else {
	    p.idstr = gretl_strndup(p.point - 1, close + 1);
	}
    }

    if (p.idstr == null) {
	if (!p.err) {
	    p.err = GretlError.alloc;
	}
    } else {
	tailstrip(p.idstr);
	parser_advance(p, close);
	lex(p);
    }

    return (p.idstr == null)? null : newstr(p.idstr);
}

private Node * get_bundle_member_name(Parser * p, int dollarize)
{
    Node *ret = null;


    if (p.ch == '[') {
	/* using bracketed key notation */
	parser_getc(p); /* eat opening '[' */
	if (!p.err) {
	    lex(p); /* get next character */
	    ret = expr(p);
	}
	if (!p.err) {
	    if (p.sym != GenSym.g_rbr) {
		unmatched_symbol_error('[', p);
	    } else {
		lex(p); /* eat closing ']' */
	    }
	}
    } else if (p.ch == '.') {
	/* using bundle dot notation */
	int i;
	int n;

	if (dollarize && *p.point == '$') {
	    /* allow leading '$' for objects under named models */
	    n = 1 + gretl_namechar_spn(p.point + 1);
	} else {
	    /* otherwise should be standard valid identifier */
	    n = gretl_namechar_spn(p.point);
	}
	if (n == 0 || n >= VNAMELEN) {
	    p.err = GretlError.parse;
	} else {
	    p.idstr = gretl_strndup(p.point, n);
	    if (p.idstr == null) {
		p.err = GretlError.alloc;
	    } else {
		for (i=0; i<=n; i++) {
		    parser_getc(p);
		}
		lex(p);
		ret = newstr(p.idstr);
	    }
	}
    } else {
	p.err = GretlError.parse;
    }

    return ret;
}

/* Find the name of the (putative) series to the right of '.'
   in an expression of the form <listname>.<series-name>,
   and return it in a string node.
*/

private Node * series_name_node(Parser * p)
{
    int n = gretl_namechar_spn(p.point);
    Node *ret = null;

    if (n == 0 || n >= VNAMELEN) {
	p.err = GretlError.parse;
    } else {
	p.idstr = gretl_strndup(p.point, n);
	if (p.idstr == null) {
	    p.err = GretlError.alloc;
	} else {
	    int i;

	    for (i=0; i<=n; i++) {
		parser_getc(p);
	    }
	    lex(p);
	    ret = newstr(p.idstr);
	}
    }

    return ret;
}

private void get_matrix_def(Node * t, Parser * p)
{
    Node *n;
    char cexp = 0;


    if (p.sym == GenSym.g_lcb) {
	lex(p);
	if (!p.err && p.sym == GenSym.g_rcb) {
	    /* empty matrix def, OK */
	    lex(p);
	    return;
	}
	while (p.ch != 0 && !p.err) {
	    n = expr(p);
	    if (p.err) {
		break;
	    } else {
		p.err = push_bn_node(t, n);
	    }
	    if (p.sym == GenSym.p_com) {
		/* comma: on to the next column */
		lex(p);
	    } else if (p.sym == GenSym.p_semi) {
		/* semicolon: on to the next row */
		n = newempty();
		p.err = push_bn_node(t, n);
		lex(p);
	    } else if (p.sym == GenSym.g_rcb) {
		/* right curly bracket: reached the end */
		break;
	    } else {
		/* something that doesn't belong here */
		cexp = '}';
		break;
	    }
	}
    } else {
	cexp = '{';
    }


    if (!p.err && cexp == 0) {
	if (p.sym == GenSym.g_rcb) {
	    lex(p);
	} else {
	    unmatched_symbol_error('{', p);
	}
    }

    if (cexp && !p.err) {
	expected_symbol_error(cexp, p, 0);
    }
}

void set_slice_on(Parser * p) { p.flags |= GenFlags.slicing; }
void set_slice_off(Parser * p) { p.flags &= ~GenFlags.slicing; }

void set_lag_parse_on(Parser * p) { p.flags |= GenFlags.lagprse; }
void set_lag_parse_off(Parser * p) { p.flags &= ~GenFlags.lagprse; }

private void get_slice_parts(Node * t, Parser * p)
{
    int slice_upsym = p.upsym;
    char cexp = 0;


    set_slice_on(p);

    if (p.sym == GenSym.g_lbr) {
	lex(p);
	if (p.sym == GenSym.p_com) {
	    /* empty row spec, OK */
	    t.L = newempty();
	} else {
	    t.L = expr(p);
	}
	if (p.sym == GenSym.p_col) {
	    /* second part of colon-separated range */
	    t.L = newb2(GenSym.subsl, t.L, null);
	    lex(p);
	    if (p.sym == GenSym.p_com || p.sym == GenSym.g_rbr) {
		/* reached end: second part implicitly empty */
		t.L.R = newempty();
	    } else {
		t.L.R = expr(p);
	    }
	}
	if (p.sym == GenSym.g_rbr) {
	    /* no comma, no second arg string: may be OK */
	    t.R = null;
	    lex(p);
	    set_slice_off(p);
	    return;
	}
	if (p.sym == GenSym.p_com) {
	    p.upsym = slice_upsym; /* 2021-09-19 */
	    lex(p);
	    if (p.sym == GenSym.g_rbr) {
		/* empty column spec, OK */
		t.R = newempty();
	    } else {
		t.R = expr(p);
	    }
	    if (p.sym == GenSym.p_col) {
		/* second part of colon-separated range */
		t.R = newb2(GenSym.subsl, t.R, null);
		lex(p);
		if (p.sym == GenSym.g_rbr) {
		    /* reached end: second part implicitly empty */
		    t.R.R = newempty();
		} else {
		    t.R.R = expr(p);
		}
	    }
	    if (p.sym == GenSym.g_rbr) {
		lex(p);
	    } else {
		cexp = ']';
	    }
	} else {
	    cexp = ',';
	}
    } else {
	cexp = '[';
    }

    if (cexp && !p.err) {
	expected_symbol_error(cexp, p, 0);
    }

    set_slice_off(p);
}

/* contexts where empty arguments are never acceptable */
bool no_empty(int f) { return f == GenFunc.f_defbundle || f == GenFunc.f_defarray || f == GenFunc.f_deflist; }

private void attach_child(Node * parent, Node * child, int f, int np, int i, Parser * p)
{
    if (p.err) {
	return;
    }


    /* catch erroneous cases */
    if (np > 0 && i == np) {
	gretl_errmsg_sprintf("%s: %s", getsymb_full(f, p),
			     _("too many arguments").ptr);
	p.err = GretlError.invarg;
    } else if (child.t == GenSym.empty && no_empty(f)) {
	p.err = GretlError.parse;
    }

    if (p.err) {
	free_tree(child, p, 0);
	return;
    }

    if (np == 1) {
	/* 1-place node */
	parent.L = child;
	child.parent = parent;
    } else if (np == 2) {
	/* 2-place node */
	if (i == 0) {
	    parent.L = child;
	} else {
	    parent.R = child;
	}
	child.parent = parent;
    } else if (np == 3) {
	/* 3-place node */
	if (i == 0) {
	    parent.L = child;
	} else if (i == 1) {
	    parent.M = child;
	} else {
	    parent.R = child;
	}
	child.parent = parent;
    } else if (child != null) {
	/* n-place node */
	p.err = push_bn_node(parent, child);
	child.parent = parent;
    }
}

private void pad_parent(Node * parent, int np, int i, Parser * p)
{
    int j;

    for (j=i; j<np; j++) {
	attach_child(parent, newempty(), 0, np, j, p);
    }
}

struct ArgRecord {
    int f;
    int[4] vec;
}

/* Here we record the positions of function-call arguments to
   built-in functions that take such arguments. These require
   special handling: they should be passed forward as string
   literals even if they are not wrapped in double-quotes.
*/

__gshared ArgRecord[] fncall_argrec = [
    ArgRecord(GenFunc.f_bfgsmax, [0, 1, 1, 0]),
    ArgRecord(GenFunc.f_nrmax, [0, 1, 1, 1]),
    ArgRecord(GenFunc.f_fdjac, [0, 1, 0, 0]),
    ArgRecord(GenFunc.f_simann, [0, 1, 0, 0]),
    ArgRecord(GenFunc.f_bfgscmax, [0, 0, 1, 1]),
    ArgRecord(GenFunc.f_nmmax, [0, 1, 0, 0]),
    ArgRecord(GenFunc.f_gssmax, [0, 1, 0, 0]),
    ArgRecord(GenFunc.f_numhess, [0, 1, 0, 0]),
    ArgRecord(GenFunc.f_fzero, [1, 0, 0, 0]),
];

private const(int) * get_callargs(int f)
{
    int i, n = cast(int) fncall_argrec.length;

    for (i=0; i<n; i++) {
	if (f == fncall_argrec[i].f) {
	    return fncall_argrec[i].vec.ptr;
	}
    }

    if (fncall_func(f)) {
	/* inconsistent with info in genparse.h */
	fprintf(stderr, "ERROR: function %s has an 'fncall' arg" ~
			" but no argrecord\n", getsymb(f));
    }

    return null;
}

private int next_arg_is_string(int i, const(int) * callargs, int k, int opt)
{
    if (callargs && callargs[i]) {
	return 1;
    }
    if ((opt & StrArg.MID_STR) && i > 0 && i < k-1) {
	return 1;
    }
    if ((opt & StrArg.RIGHT_STR) && i == k-1) {
	return 1;
    }
    return 0;
}

/* Get up to @np comma-separated arguments (possibly optional).
   However, if np < 0 this is a signal to get as many arguments as we
   can find (the number unknown in advance).
*/

private void get_args(Node * t, Parser * p, int f, int np, int opt, int * next)
{
    Node *child;
    const(int) * callargs;
    char cexp = 0;
    int i = 0;


    if (p.sym != GenSym.g_lpr) {
	expected_symbol_error('(', p, 0);
	return;
    }

    if (f == GenFunc.f_stack) {
	/* for stack(), the first arg should be a list */
	p.flags |= GenFlags.listdef;
    }

    callargs = get_callargs(f);
    if (callargs == null || callargs[0] == 0) {
	/* nothing special about the first arg */
	lex(p);
    }

    if (f == GenFunc.f_stack) {
	/* revert to regular parsing */
	p.flags &= ~GenFlags.listdef;
    }

    while (p.sym != GenSym.g_rpr && p.ch != '\0' && !p.err) {
	/* get the next argument */
	if (p.sym == GenSym.p_com) {
	    /* implies an empty argument slot */
	    child = newempty();
	} else if (i < 4 && callargs && callargs[i]) {
	    /* a function-call argument: don't insist on quotation */
	    child = get_literal_string_arg(p, 0);
	} else if (i > 0 && i < np - 1 && (opt & StrArg.MID_STR)) {
	    child = get_literal_string_arg(p, opt);
	} else if (i == np - 1 && (opt & StrArg.RIGHT_STR)) {
	    child = get_final_string_arg(p, t, f, 0);
	} else {
	    child = expr(p);
	}

	if (!p.err) {
	    attach_child(t, child, f, np, i++, p);
	}

	if (p.err || p.sym == GenSym.g_rpr) {
	    break;
	} else if (p.sym == GenSym.p_com) {
	    /* turn off flag for accepting string as first arg */
	    p.flags &= ~GenFlags.getstr;
	    if (next_arg_is_string(i, callargs, np, opt)) {
		/* the next arg needs special handling */
		p.sym = GenSym.empty;
	    } else {
		lex(p);
		if (p.sym == GenSym.g_rpr) {
		    /* trailing empty argument slot */
		    attach_child(t, newempty(), f, np, i++, p);
		}
	    }
	} else {
	    /* either arg-separating comma or closing paren was expected */
	    cexp = (i < np)? ',' : ')';
	    break;
	}
    }

    if (!p.err) {
	if (cexp) {
	    expected_symbol_error(cexp, p, 0);
	} else if (p.sym == GenSym.g_rpr) {
	    if (i < np) {
		pad_parent(t, np, i, p);
	    }
	    lex(p);
	    if (p.sym == GenSym.g_lbr) {
		*next = '[';
	    }
	} else {
	    expected_symbol_error(')', p, 0);
	}
    }

}

private void get_assertion(Node * t, Parser * p)
{
    const char *s = strrchr(p.point, ')');
    char * str;
    int n = cast(int) (s - p.point + 1);
    int next = 0;

    str = gretl_strndup(p.point - 1, n);

    get_args(t, p, GenFunc.f_assert, 1, 0, &next);
    if (!p.err) {
	str = g_strchomp(g_strchug(str));
	t.R = newstr(str);
    }
}

/* For defining a bundle via _().ptr: get one or more comma-separated
   terms: each must take the form key=value, or just key if the
   key and the name of a pre-defined object are one and the same.
*/

private void get_bundle_pairs(Node * t, Parser * p, int * next)
{
    Node *child;
    const(char) * src;
    int n;
    int j;
    int i = 0;


    /* allow for an empty definition */
    while (p.ch == ' ') parser_getc(p);
    if (p.ch == ')') {
	parser_getc(p);
	lex(p);
	return;
    }

    while (p.ch && !p.err) {
	/* first get an unquoted key */
	while (p.ch == ' ') parser_getc(p);
	src = p.point - 1;
	n = gretl_namechar_spn(src);
	if (n == 0) {
	    p.err = GretlError.parse;
	    break;
	}
	for (j=0; j<n; j++) {
	    parser_getc(p);
	}
	p.idstr = gretl_strndup(src, n);
	child = newstr(p.idstr);
	attach_child(t, child, 0, -1, i++, p);
	while (p.ch == ' ') parser_getc(p);
	if (p.ch == '=') {
	    /* parse the following expression */
	    parser_getc(p);
	    while (p.ch == ' ') parser_getc(p);
	    lex(p);
	    child = expr(p);
	} else {
	    /* backtrack to get named object */
	    parser_advance(p, cast(int) -(p.point - src));
	    lex(p);
	    child = base(p, null);
	}
	if (!p.err) {
	    attach_child(t, child, 0, -1, i++, p);
	}
	if (p.sym == GenSym.p_com) {
	    ; /* OK */
	} else if (p.sym == GenSym.g_rpr) {
	    lex(p);
	    break;
	} else {
	    gretl_errmsg_set(_("Expected ',' or ')'").ptr);
	    p.err = GretlError.parse;
	}
    }
}

private Node * powterm(Parser * p, Node * l)
{
    int sym = p.sym;
    int opt = Gretlopt.OPT_NONE;
    int next = 0;
    Node *t = null;

    if (p.err) {
	return null;
    }


    if (l != null) {
	/* powterm recursion: swallowing prior node @l */
	if (sym == GenSym.bmemb || sym == GenSym.dbmemb) {
	    t = newb2(sym, l, null);
	    if (t != null) {
		parser_ungetc(p);
		t.R = get_bundle_member_name(p, 0);
	    }
	} else if (sym == GenSym.g_lbr) {
	    /* "OSL": we're being somewhat agnostic here regarding
	       the type of object of which we're taking a slice;
	       That will be sorted out at eval() time.
	    */
	    t = newb2(GenSym.osl, l, null);
	    if (t != null) {
		t.R = new_node(GenSym.slraw);
		if (t.R != null) {
		    get_slice_parts(t.R, p);
		}
	    }
	} else if (sym == GenSym.g_lpr) {
	    /* (bundle member or listvar) plus lag spec? */
	    if (l.t == GenSym.bmemb || l.t == GenSym.listvar) {
		set_lag_parse_on(p);
		t = newb2(GenSym.lag, l, null);
		if (t != null) {
		    /* 2024-05-02: call base() instead of expr(): this
		       fixes the bug reported in msg17192.html at
		       https://www.mail-archive.com/gretl-users@gretlml.univpm.it/
		    */
		    t.R = base(p, null);
		}
		set_lag_parse_off(p);
	    } else {
		context_error('(', p, "powterm");
	    }
	}
	goto maybe_recurse;
    }

    if (string_last_func(sym)) {
	opt |= StrArg.RIGHT_STR;
    }
    if (string_mid_func(sym)) {
	opt |= StrArg.MID_STR;
    }

    if (sym == GenSym.u_addr) {
        t = new_node(sym);
        if (t != null) {
            lex(p);
	    t.L = u_addr_base(p);
        }
    } else if (func2_symb(sym)) {
	int unset = 0;

	if (sym == GenFunc.f_genseries) {
	    set_doing_genseries(1);
	    unset = 1;
	}
	t = new_node(sym);
	if (t != null) {
	    lex(p);
	    if (sym == GenFunc.f_assert) {
		get_assertion(t, p);
	    } else {
		get_args(t, p, sym, 2, opt, &next);
	    }
	}
	if (unset) {
	    set_doing_genseries(0);
	}
    } else if (func3_symb(sym)) {
	t = newb3(sym, null);
	if (t != null) {
	    lex(p);
	    if (str0_func(sym)) {
		p.flags |= GenFlags.getstr;
	    }
	    get_args(t, p, sym, 3, opt, &next);
	}
    } else if (str0_func(sym)) {
	t = new_node(sym);
	if (t != null) {
	    lex(p);
	    t.L = newbn(GenSym.fargs);
	    if (t != null) {
		p.flags |= GenFlags.getstr;
		get_args(t.L, p, sym, -1, opt, &next);
	    }
	}
    } else if (string_arg_func(sym)) {
	t = new_node(sym);
	if (t != null) {
	    lex(p);
	    t.L = get_final_string_arg(p, t, sym, 1);
	}
    } else if (func1_symb(sym)) {
	if (undef_arg_ok(sym)) {
	    set_parsing_query(1);
	}
	t = new_node(sym);
	if (t != null) {
	    if (sym < GenFunc.fp_max) {
		t.v.ptr = p.data; /* attach function pointer */
	    }
	    lex(p);
	    get_args(t, p, sym, 1, opt, &next);
	}
	if (undef_arg_ok(sym)) {
	    set_parsing_query(0);
	}
    } else if (sym == GenSym.lag || sym == GenSym.obs) {
	if (sym == GenSym.lag) {
	    set_lag_parse_on(p);
	}
	t = new_node(sym);
	if (t != null) {
	    t.L = newref(p, p.upsym);
	    lex(p);
	    t.R = base(p, t);
	}
	if (sym == GenSym.lag) {
	    set_lag_parse_off(p);
	}
    } else if (sym == GenSym.osl) {
	t = new_node(sym);
	if (t != null) {
	    t.L = newref(p, p.upsym);
	    t.R = new_node(GenSym.slraw);
	    if (t.R != null) {
		lex(p);
		get_slice_parts(t.R, p);
	    }
	}
    } else if (sym == GenSym.dmstr) {
	t = new_node(sym);
	if (t != null) {
	    t.L = newref(p, GenSym.mvar);
	    lex(p);
	    t.R = get_final_string_arg(p, t, sym, 1);
	}
    } else if (sym == GenSym.bmemb) {
	if (p.data == null) {
	    p.err = GretlError.parse;
	    return null;
	}
	t = new_node(sym);
	if (t != null) {
	    t.L = newref(p, GenSym.bundle);
	    t.R = get_bundle_member_name(p, 0);
	}
    } else if (sym == GenSym.dbmemb) {
	t = new_node(sym);
	if (t != null) {
	    t.L = newref(p, GenSym.dbundle);
	    t.R = get_bundle_member_name(p, 0);
	}
    } else if (sym == GenSym.mmemb) {
	t = new_node(GenSym.dbmemb);
	if (t != null) {
	    t.L = newref(p, GenSym.mmemb);
	    t.R = get_bundle_member_name(p, 1);
	}
    } else if (sym == GenSym.listvar) {
	t = new_node(GenSym.listvar);
	if (t != null) {
	    t.L = newref(p, GenSym.list);
	    t.R = series_name_node(p);
	}
    } else if (sym == GenSym.g_lpr) {
	/* parenthesized expression */
	t = base(p, null);
    } else if (sym == GenSym.g_lcb) {
	/* explicit matrix definition */
	t = new_node(GenFunc.f_defmat);
	if (t != null) {
	    t.L = newbn(GenSym.fargs);
	    if (t.L != null) {
		get_matrix_def(t.L, p);
	    }
	}
    } else if (sym == GenSym.ufun || sym == GenSym.rfun) {
        t = newref(p, sym);
        if (t != null) {
            lex(p);
            t.R = newbn(GenSym.fargs);
	    if (t.R != null) {
		get_args(t.R, p, sym, -1, opt, &next);
	    }
        }
    } else if (sym == GenFunc.f_defargs) {
	t = new_node(sym);
	if (t != null) {
	    lex(p);
	    t.L = newbn(GenSym.fargs);
	    if (t.L != null) {
		get_bundle_pairs(t.L, p, &next);
	    }
	}
    } else if (funcn_symb(sym)) {
	t = new_node(sym);
	if (t != null) {
	    lex(p);
	    t.L = newbn(GenSym.fargs);
	    if (t.L != null) {
		int k = -1;

		if (sym == GenFunc.f_nrmax ||
		    sym == GenFunc.f_bfgscmax ||
		    sym == GenFunc.f_movavg) {
		    k = 4;
		}
		get_args(t.L, p, sym, k, opt, &next);
	    }
	}
    } else {
	t = base(p, null);
    }

    if (als_func(sym) && (p.flags & GenFlags.aliased)) {
	/* transfer flag to newly created node */
	t.flags |= NodeFlags.als;
	p.flags ^= GenFlags.aliased;
    }

 maybe_recurse:

    if (!p.err) {
	if (t != null && next == '[') {
	    /* support func(args)[slice], etc. */
	    t = newb2(GenSym.osl, t, null);
	    if (t != null) {
		t.R = new_node(GenSym.slraw);
		if (t.R != null) {
		    get_slice_parts(t.R, p);
		}
	    }
	} else if (t != null) {
	    if (p.sym == GenSym.bmemb || p.sym == GenSym.dbmemb) {
		/* these types can recurse */
		t = powterm(p, t);
	    } else if (p.sym == GenSym.g_lbr || p.sym == GenSym.g_lpr) {
		t = powterm(p, t);
	    }
	}
    }


    return t;
}

private int got_hex_val(Parser * p)
{
    int c = *p.point;

    return p.ch == '0' && (c == 'x' || c == 'X');
}

/* Test for whether the ' symbol must represent the unary
   transposition operator rather than binary transpose-multiply,
   based on the following symbol, @t.
*/

bool must_be_unary(int t)
{
    return t == GenSym.eot || t < GenSym.op_max || t == GenSym.p_com || t == GenSym.g_rpr || t == GenSym.g_rbr || t == GenSym.g_rcb || t == GenSym.p_col;
}

private Node * factor(Parser * p)
{
    int sym = p.sym == GenSym.b_sub ? GenSym.u_neg :
	p.sym == GenSym.b_add ? GenSym.u_pos : p.sym;
    Node *t;


    if (sym == GenSym.u_neg && got_hex_val(p)) {
	gretl_errmsg_set("hexadecimal values must be unsigned");
	p.err = GretlError.types;
    }
    if (p.err) {
	return null;
    }

    if (unary_op(sym) && sym != GenSym.u_addr) {
        t = new_node(sym);
        if (t != null) {
            lex(p);
	    t.L = factor(p);
        }
    } else {
	t = powterm(p, null);
	if (t != null) {
	    int upsym = p.sym;

	    while (!p.err && (p.sym == GenSym.b_pow ||
			       p.sym == GenSym.b_dotpow ||
			       p.sym == GenSym.b_trmul)) {
		t = newb2(p.sym, t, null);
		if (t != null) {
		    lex(p);
		    if (upsym == GenSym.b_trmul && must_be_unary(p.sym)) {
			/* dummy RHS for unary transpose */
			t.R = newempty();
		    } else {
			t.R = factor(p);
		    }
		}
	    }
	}
    }


    return t;
}

private Node * term(Parser * p)
{
    Node *t;


    if (p.err || (t = factor(p)) == null) {
	return null;
    }

    while (!p.err && (p.sym == GenSym.b_mul || p.sym == GenSym.b_div ||
		       p.sym == GenSym.b_ldiv || p.sym == GenSym.b_mod ||
		       p.sym == GenSym.b_dotmult || p.sym == GenSym.b_dotdiv ||
		       p.sym == GenSym.b_kron)) {
	t = newb2(p.sym, t, null);
	if (t != null) {
	    lex(p);
	    t.R = factor(p);
	}
    }


    return t;
}

private Node * expr4(Parser * p)
{
    Node *t;


    if (p.err || (t = term(p)) == null) {
	return null;
    }

    while (!p.err && (p.sym == GenSym.b_add || p.sym == GenSym.b_sub ||
		       p.sym == GenSym.b_dotadd || p.sym == GenSym.b_dotsub ||
		       p.sym == GenSym.b_hcat || p.sym == GenSym.b_vcat ||
		       p.sym == GenSym.b_lcat || p.sym == GenSym.b_range ||
		       p.sym == GenSym.b_ellip || p.sym == GenSym.b_join)) {
	t = newb2(p.sym, t, null);
	if (t != null) {
	    lex(p);
	    t.R = term(p);
	}
    }


    return t;
}

private Node * expr3(Parser * p)
{
    Node *t;


    if (p.err || (t = expr4(p)) == null) {
	return null;
    }

    while (!p.err && (p.sym == GenSym.b_gt || p.sym == GenSym.b_lt ||
		       p.sym == GenSym.b_dotgt || p.sym == GenSym.b_dotlt ||
		       p.sym == GenSym.b_gte || p.sym == GenSym.b_lte ||
		       p.sym == GenSym.b_dotgte || p.sym == GenSym.b_dotlte)) {
	t = newb2(p.sym, t, null);
	if (t != null) {
	    lex(p);
	    t.R = expr4(p);
	}
    }


    return t;
}

private Node * expr2(Parser * p)
{
    Node *t;


    if (p.err || (t = expr3(p)) == null) {
	return null;
    }

    while (!p.err && (p.sym == GenSym.b_eq || p.sym == GenSym.b_neq ||
		       p.sym == GenSym.b_doteq || p.sym == GenSym.b_dotneq)) {
	t = newb2(p.sym, t, null);
	if (t != null) {
	    lex(p);
	    t.R = expr3(p);
	}
    }


    return t;
}

private Node * expr1(Parser * p)
{
    Node *t;


    if (p.err || (t = expr2(p)) == null) {
	return null;
    }

    while (!p.err && p.sym == GenSym.b_and) {
	t = newb2(p.sym, t, null);
	if (t != null) {
	    p.flags |= GenFlags.and;
	    lex(p);
	    t.R = expr2(p);
	    p.flags &= ~GenFlags.and;
	}
    }


    return t;
}

private Node * expr0(Parser * p)
{
    Node *t;


    if (p.err || (t = expr1(p)) == null) {
	return null;
    }

    while (!p.err && p.sym == GenSym.b_or) {
	t = newb2(p.sym, t, null);
	if (t != null) {
	    p.flags |= GenFlags.or;
	    lex(p);
	    t.R = expr1(p);
	    p.flags &= ~GenFlags.or;
	}
    }


    return t;
}

extern (C) Node * expr(Parser * p)
{
    Node *t;


    if (p.err || (t = expr0(p)) == null) {
	return null;
    }

    while (!p.err && p.sym == GenSym.query) {
	t = newb3(p.sym, t);
	if (t != null) {
	    set_parsing_query(1);
	    lex(p);
	    if (!p.err) {
		t.M = expr(p);
		if (p.sym == GenSym.p_col) {
		    lex(p);
		    if (!p.err) {
			t.R = expr(p);
		    }
		} else {
		    expected_symbol_error(':', p, 0);
		}
	    }
	    set_parsing_query(0);
	}
    }


    return t;
}
