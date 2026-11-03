/* Port of lib/src/gretl_color.c: mixing of two ARGB colors (the colormix() function), splitting a packed color
   into its channels, and the gnuplot display of the mixes. */

module gretl.gretl_color;


import core.stdc.math : nearbyint;
import core.stdc.stdio : FILE, fprintf, fputc, fputs, sprintf;
import core.stdc.stdlib : free, malloc;
import core.stdc.string : memset;

import gretl.gretl_array : gretl_array_new, gretl_array_set_string;
import gretl.graphing : finalize_plot_input_file, open_plot_input_file, print_rgb_hash;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.graphingh : GptFlags, PlotType, gretlRGB;
import gretl.libgretlh;
import gretl.libset : gretl_in_gui_mode;
import gretl.options : set_optval_string;

/* unported dependencies: per-file declarations from the 2026c headers */
extern (C) {
    void manufacture_gui_callback(int ci);
}

struct ARGB {
    ubyte a; /* alpha channel [0,255] */
    ubyte r; /* red channel [0,255] */
    ubyte g; /* green channel [0,255] */
    ubyte b; /* blue channel [0,255] */
}

union Channels {
    uint u32;
    ubyte[4] u8;
}

/* rgb color mixing: (1-f) * c1 + f * c2 */

private ARGB argb_mix(ARGB * c1, ARGB * c2, double f)
{
    ARGB ret;
    double x;

    memset(&ret, 0, ret.sizeof);

    x = (1-f) * c1.a + f * c2.a;
    ret.a = cast(ubyte) nearbyint(x);
    x = (1-f) * c1.r + f * c2.r;
    ret.r = cast(ubyte) nearbyint(x);
    x = (1-f) * c1.g + f * c2.g;
    ret.g = cast(ubyte) nearbyint(x);
    x = (1-f) * c1.b + f * c2.b;
    ret.b = cast(ubyte) nearbyint(x);

    return ret;
}

private ARGB ARGB_from_guint32(uint u)
{
    Channels ch;
    ARGB output;

    ch.u32 = u;
    output.a = ch.u8[3];
    output.r = ch.u8[2];
    output.g = ch.u8[1];
    output.b = ch.u8[0];

    return output;
}

void decompose_argb(uint u, ubyte * a, ubyte * r, ubyte * g, ubyte * b)
{
    Channels ch;

    ch.u32 = u;
    *a = ch.u8[3];
    *r = ch.u8[2];
    *g = ch.u8[1];
    *b = ch.u8[0];
}

GretlArray * colormix_array(gretlRGB c0, gretlRGB c1, const(double) * f, int nf, int do_plot, int * err)
{
    GretlArray * ret = null;
    uint * uvec = null;
    char[12] color = 0;
    ARGB argb0, argb1, mix;
    uint u;
    int i;

    ret = gretl_array_new(GretlType.strings, nf, err);
    if (*err) {
        return null;
    }

    if (do_plot) {
        uvec = cast(uint *) malloc(nf * (*uvec).sizeof);
    }

    argb0 = ARGB_from_guint32(c0);
    argb1 = ARGB_from_guint32(c1);

    for (i=0; i<nf; i++) {
        mix = argb_mix(&argb0, &argb1, f[i]);
        if (mix.a > 0) {
            u = (cast(uint) mix.a << 24) | (cast(uint) mix.r << 16) | (cast(uint) mix.g << 8) | mix.b;
            sprintf(color.ptr, "0x%08x", u);
        } else {
            u = (cast(uint) mix.r << 16) | (cast(uint) mix.g << 8) | mix.b;
            sprintf(color.ptr, "0x%06x", u);
        }
        gretl_array_set_string(ret, i, color.ptr, 1);
        if (uvec != null) {
            uvec[i] = u;
        }
    }

    if (*err == 0 && uvec != null) {
        show_mixed_colors(c0, c1, uvec, f, nf);
    }

    free(uvec);

    return ret;
}

private int show_mixed_colors(gretlRGB c0, gretlRGB c1, uint * mix, const(double) * f, int n)
{
    FILE * fp;
    char[10][3] cs;
    int i, j;
    int err = 0;

    set_optval_string(GretlCmdIndex.gnuplot, Gretlopt.OPT_U, "display");
    fp = open_plot_input_file(PlotType.regular, cast(GptFlags) 0, &err);
    if (err) {
        return err;
    }

    print_rgb_hash(cs[0].ptr, c0);
    print_rgb_hash(cs[1].ptr, c1);
    print_rgb_hash(cs[2].ptr, mix[0]);

    fprintf(fp, "set title \"mix %s and %s via RGB\\n", cs[0].ptr, cs[1].ptr);
    if (n == 1) {
        fprintf(fp, "f = %g\"\n", f[0]);
        fprintf(fp, "set xtics (\"%s\" 1, \"%s\" 2, \"%s\" 3)\n", cs[0].ptr, cs[1].ptr, cs[2].ptr);
    } else {
        fprintf(fp, "f = ");
        for (i=0; i<n; i++) {
            fprintf(fp, "%g", f[i]);
            fputs(i < n-1 ? ", " : "\"\n", fp);
        }
        fputs("set xtics (", fp);
        for (i=0; i<n; i++) {
            print_rgb_hash(cs[2].ptr, mix[i]);
            fprintf(fp, "\"%s\" %d", cs[2].ptr, i+1);
            fputs(i < n-1 ? ", " : ")\n", fp);
        }
    }
    fputs("set border 0\n", fp);
    fputs("set nokey\n", fp);
    fputs("set noytics\n", fp);
    fputs("set boxwidth 0.9 relative\n", fp);
    fputs("set style fill solid 1.0\n", fp);
    if (n == 1) {
        /* show c0, c1, mix */
        fprintf(fp, "set linetype 1 lc rgb \"%s\"\n", cs[0].ptr);
        fprintf(fp, "set linetype 2 lc rgb \"%s\"\n", cs[1].ptr);
        fprintf(fp, "set linetype 3 lc rgb \"%s\"\n", cs[2].ptr);
        n = 3;
    } else {
        /* just show the mixes */
        for (i=0; i<n; i++) {
            print_rgb_hash(cs[2].ptr, mix[i]);
            fprintf(fp, "set linetype %d lc rgb \"%s\"\n", i+1, cs[2].ptr);
        }
    }
    fputs("$data << EOD\n", fp);
    for (i=0; i<n; i++) {
        fprintf(fp, "%d ", i+1);
        for (j=0; j<n; j++) {
            fprintf(fp, "%d", j == i ? 1 : 0);
            fputc(j < n-1 ? ' ' : '\n', fp);
        }
    }
    fputs("EOD\n", fp);
    fputs("plot \\\n", fp);
    for (i=0; i<n; i++) {
        fprintf(fp, " '$data' using 1:%d w boxes", i+2);
        fputs(i < n-1 ? ", \\\n" : "\n", fp);
    }

    err = finalize_plot_input_file(fp);
    if (!err && gretl_in_gui_mode()) {
        manufacture_gui_callback(GretlCmdIndex.gnuplot);
    }

    return err;
}
