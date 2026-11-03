/*
   Replacements for the parts of glib's GDate that calendar.c uses, written in plain D (no glib):
   g_date_clear, g_date_valid_dmy, g_date_valid_julian, g_date_set_dmy, g_date_set_julian,
   g_date_get_julian/_weekday/_day/_month/_year/_day_of_year, g_date_is_leap_year, g_date_new_julian,
   g_date_free, g_date_get_iso8601_week_of_year and g_date_strftime.

   As in glib, a "Julian day" here is the day number counted from 1 for 1 January, AD 1, on the
   proleptic Gregorian calendar (not an astronomers' Julian day number); the weekday numbering is
   Monday = 1 to Sunday = 7 (G_DATE_MONDAY..G_DATE_SUNDAY); valid years are 1 to 65535 (G_MAXUINT16).
   Conversions use the days-from-civil algorithm. g_date_strftime() calls the C library's strftime()
   on a struct tm built from the date.

   Used by: calendar.d.
*/

module gretl.gdate;


import core.stdc.time : strftime, tm;

import gretl.glib : gboolean, gsize;

enum G_DATE_BAD_JULIAN = 0;
enum G_DATE_MONDAY = 1;
enum G_DATE_SATURDAY = 6;
enum G_DATE_SUNDAY = 7;

struct GDate {
    uint julian_days;
    bool julian;        /* is the julian_days member valid? */
}

/* Hinnant's days_from_civil/civil_from_days: days relative to 1970-01-01 */

private long days_from_civil(long y, uint m, uint d)
{
    y -= m <= 2;
    long era = (y >= 0 ? y : y - 399) / 400;
    uint yoe = cast(uint) (y - era * 400);
    uint doy = (153 * (m > 2 ? m - 3 : m + 9) + 2) / 5 + d - 1;
    uint doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;

    return era * 146097 + cast(long) doe - 719468;
}

private void civil_from_days(long z, out long y, out uint m, out uint d)
{
    z += 719468;
    long era = (z >= 0 ? z : z - 146096) / 146097;
    uint doe = cast(uint) (z - era * 146097);
    uint yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;

    y = cast(long) yoe + era * 400;
    uint doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    uint mp = (5 * doy + 2) / 153;

    d = doy - (153 * mp + 2) / 5 + 1;
    m = mp < 10 ? mp + 3 : mp - 9;
    y += m <= 2;
}

/* 1 for 0001-01-01 */
private enum EPOCH_OFFSET = 719163;

private bool leap(int y)
{
    return (y % 4 == 0 && y % 100 != 0) || y % 400 == 0;
}

private int month_days(int y, int m)
{
    static immutable int[13] md = [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];

    return (m == 2 && leap(y)) ? 29 : md[m];
}

void g_date_clear(GDate * date, uint n_dates)
{
    foreach (i; 0 .. n_dates) {
        date[i] = GDate.init;
    }
}

gboolean g_date_valid_dmy(ubyte day, int month, int year)
{
    return day >= 1 && month >= 1 && month <= 12 && year >= 1 && year <= 65535 &&
        day <= month_days(year, month);
}

gboolean g_date_valid_julian(uint julian_date)
{
    return julian_date != G_DATE_BAD_JULIAN;
}

void g_date_set_dmy(GDate * date, ubyte day, int month, int year)
{
    if (!g_date_valid_dmy(day, month, year)) {
        *date = GDate.init;
        return;
    }
    date.julian_days = cast(uint) (days_from_civil(year, month, day) + EPOCH_OFFSET);
    date.julian = true;
}

void g_date_set_julian(GDate * date, uint julian_date)
{
    if (!g_date_valid_julian(julian_date)) {
        *date = GDate.init;
        return;
    }
    date.julian_days = julian_date;
    date.julian = true;
}

uint g_date_get_julian(const(GDate) * date)
{
    return date.julian ? date.julian_days : G_DATE_BAD_JULIAN;
}

private void dmy_of(const(GDate) * date, out int y, out int m, out int d)
{
    long yy;
    uint mm, dd;

    civil_from_days(cast(long) date.julian_days - EPOCH_OFFSET, yy, mm, dd);
    y = cast(int) yy;
    m = cast(int) mm;
    d = cast(int) dd;
}

int g_date_get_weekday(const(GDate) * date)
{
    if (!date.julian) {
        return 0;
    }
    return cast(int) ((date.julian_days - 1) % 7) + 1;
}

ubyte g_date_get_day(const(GDate) * date)
{
    int y, m, d;

    if (!date.julian) {
        return 0;
    }
    dmy_of(date, y, m, d);
    return cast(ubyte) d;
}

int g_date_get_month(const(GDate) * date)
{
    int y, m, d;

    if (!date.julian) {
        return 0;
    }
    dmy_of(date, y, m, d);
    return m;
}

ushort g_date_get_year(const(GDate) * date)
{
    int y, m, d;

    if (!date.julian) {
        return 0;
    }
    dmy_of(date, y, m, d);
    return cast(ushort) y;
}

uint g_date_get_day_of_year(const(GDate) * date)
{
    int y, m, d;

    if (!date.julian) {
        return 0;
    }
    dmy_of(date, y, m, d);
    return date.julian_days - cast(uint) (days_from_civil(y, 1, 1) + EPOCH_OFFSET) + 1;
}

gboolean g_date_is_leap_year(ushort year)
{
    return leap(year);
}

GDate * g_date_new_julian(uint julian_day)
{
    import core.stdc.stdlib : malloc;

    GDate * d = cast(GDate *) malloc(GDate.sizeof);

    *d = GDate.init;
    g_date_set_julian(d, julian_day);
    return d;
}

void g_date_free(GDate * date)
{
    import core.stdc.stdlib : free;

    free(date);
}

private int weeks_in_year(int y)
{
    int p = (y + y / 4 - y / 100 + y / 400) % 7;
    int q = ((y - 1) + (y - 1) / 4 - (y - 1) / 100 + (y - 1) / 400) % 7;

    return (p == 4 || q == 3) ? 53 : 52;
}

uint g_date_get_iso8601_week_of_year(const(GDate) * date)
{
    int y, m, d;

    if (!date.julian) {
        return 0;
    }
    dmy_of(date, y, m, d);
    int wd = g_date_get_weekday(date);
    int w = (cast(int) g_date_get_day_of_year(date) - wd + 10) / 7;

    if (w < 1) {
        return weeks_in_year(y - 1);
    } else if (w > weeks_in_year(y)) {
        return 1;
    }
    return w;
}

gsize g_date_strftime(char * s, gsize slen, const(char) * format, const(GDate) * date)
{
    tm t;
    int y, m, d;

    if (!date.julian || slen == 0 || format[0] == '\0') {
        return 0;
    }
    dmy_of(date, y, m, d);
    t.tm_year = y - 1900;
    t.tm_mon = m - 1;
    t.tm_mday = d;
    t.tm_wday = g_date_get_weekday(date) % 7;
    t.tm_yday = cast(int) g_date_get_day_of_year(date) - 1;
    t.tm_isdst = -1;
    return strftime(s, slen, format, &t);
}

unittest {
    import core.stdc.string : strcmp;

    GDate g;

    g_date_clear(&g, 1);
    g_date_set_dmy(&g, 1, 1, 1);
    assert(g_date_get_julian(&g) == 1 && g_date_get_weekday(&g) == 1);   /* a Monday */
    g_date_set_dmy(&g, 1, 1, 1970);
    assert(g_date_get_julian(&g) == 719163 && g_date_get_weekday(&g) == 4);   /* a Thursday */
    g_date_set_dmy(&g, 29, 2, 2024);
    assert(g_date_get_day(&g) == 29 && g_date_get_month(&g) == 2 && g_date_get_year(&g) == 2024);
    assert(g_date_get_day_of_year(&g) == 60);
    g_date_set_julian(&g, g_date_get_julian(&g) + 1);
    assert(g_date_get_day(&g) == 1 && g_date_get_month(&g) == 3);
    assert(!g_date_valid_dmy(29, 2, 2023) && g_date_valid_dmy(29, 2, 2000) && !g_date_valid_dmy(29, 2, 1900));
    assert(!g_date_valid_dmy(1, 13, 2000) && !g_date_valid_dmy(0, 1, 2000) && !g_date_valid_dmy(1, 1, 0));
    assert(!g_date_valid_julian(0));

    /* round trip over a long range */
    foreach (j; 1 .. 800_000) {
        GDate a, b;
        int y, m, d;

        g_date_set_julian(&a, j);
        dmy_of(&a, y, m, d);
        g_date_set_dmy(&b, cast(ubyte) d, m, y);
        assert(g_date_get_julian(&b) == j);
        if (j > 800_000 - 3) break;
    }

    /* ISO weeks: 2021-01-03 is week 53 of 2020; 2018-12-31 is week 1 of 2019; 2024-06-17 is week 25 */
    g_date_set_dmy(&g, 3, 1, 2021);
    assert(g_date_get_iso8601_week_of_year(&g) == 53);
    g_date_set_dmy(&g, 31, 12, 2018);
    assert(g_date_get_iso8601_week_of_year(&g) == 1);
    g_date_set_dmy(&g, 17, 6, 2024);
    assert(g_date_get_iso8601_week_of_year(&g) == 25);

    char[32] buf = 0;

    g_date_set_dmy(&g, 5, 3, 2024);
    assert(g_date_strftime(buf.ptr, buf.length, "%Y-%m-%d %a", &g) == 14 && strcmp(buf.ptr, "2024-03-05 Tue") == 0);
}
