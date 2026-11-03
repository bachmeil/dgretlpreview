/* Tests for gretl_commands.d: the command-word table agrees with the GretlCmdIndex enum in both directions, aliases,
   the help-index hack for "tsplots" and command completion.

   Build and run with testing/run_tests.sh. */

import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;
import core.stdc.string : strcmp;

import gretl.gretl_commands;
import gretl.gretl_commandsh : GretlCmdIndex;

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

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    bool roundtrip = true;
    foreach (int i; 1 .. cast(int) GretlCmdIndex.nc) {
        const(char) * w = gretl_command_word(i);
        if (w is null || *w == '\0' || gretl_command_number(w) != i) {
            printf("  mismatch at %d\n", i);
            roundtrip = false;
        }
    }
    check(roundtrip, "every command word maps back to its enum index");

    check(strcmp(gretl_command_word(GretlCmdIndex.ols), "ols") == 0, "word of ols");
    check(strcmp(gretl_command_word(GretlCmdIndex.gridplot), "gridplot") == 0, "word of gridplot");
    check(strcmp(gretl_command_word(GretlCmdIndex.ivreg), "tsls") == 0, "ivreg is spelled tsls");
    check(*gretl_command_word(0) == '\0' && *gretl_command_word(GretlCmdIndex.nc) == '\0' && *gretl_command_word(-3) == '\0', "out-of-range indices give the empty string");

    check(gretl_command_number("ols") == GretlCmdIndex.ols && gretl_command_number("nonesuch") == 0, "lookup of known and unknown words");
    check(gretl_command_number("series") == GretlCmdIndex.genr && gretl_command_number("hausman") == GretlCmdIndex.panspec, "aliases resolve");
    check(word_is_genr_alias("matrix") == 1 && word_is_genr_alias("ols") == 0 && word_is_genr_alias("arma") == 0, "word_is_genr_alias");

    check(gretl_help_index("ols") == GretlCmdIndex.ols && gretl_help_index("tsplots") == GretlCmdIndex.nc && gretl_help_index("zzz") == 0, "help index and the tsplots special case");

    const(char) * c = gretl_command_complete("scat");
    check(c !is null && strcmp(c, "scatters") == 0, "completion of a prefix");
    check(gretl_command_complete("zzz") is null, "completion with no match");
    c = gretl_command_complete_next("s", GretlCmdIndex.scatters + 1);
    check(c !is null && strcmp(c, "sdiff") == 0, "completion continues from an index");

    gretl_command_hash_cleanup();

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
