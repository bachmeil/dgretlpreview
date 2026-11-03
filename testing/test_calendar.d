/* Smoke tests for calendar.d: epoch days and their round trips (Gregorian and Julian calendars),
   weekdays, daily dataset calendars with 5/6/7-day weeks, Excel serial dates, ISO weeks, Easter,
   strftime/strptime wrappers.

   Build and run with testing/run_tests.sh. */

import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;
import core.stdc.stdlib : free;
import core.stdc.string : strcmp, strcpy;

import gretl.calendar;
import gretl.dataset : create_new_dataset, destroy_dataset;
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

void test_epoch()
{
    check(epoch_day_from_ymd(1, 1, 1) == 1 && epoch_day_from_ymd(1970, 1, 1) == 719163, "epoch_day_from_ymd");
    check(epoch_day_from_ymd(2023, 2, 29) == 0 && epoch_day_from_ymd(2024, 2, 29) != 0, "invalid dates give 0");
    check(get_epoch_day("2024-03-05") == epoch_day_from_ymd(2024, 3, 5), "get_epoch_day (ISO)");
    check(get_epoch_day("24-03-05") == epoch_day_from_ymd(2024, 3, 5) && get_epoch_day("20240305") == epoch_day_from_ymd(2024, 3, 5) &&
          get_epoch_day("2024/03/05") == epoch_day_from_ymd(2024, 3, 5), "get_epoch_day (2-digit year, basic, slashes)");
    check(get_epoch_day("2024-13-05") == 0 && get_epoch_day("garbage") == 0, "get_epoch_day rejects bad input");

    int y, m, d;
    uint ed = epoch_day_from_ymd(1999, 12, 31);

    check(ymd_bits_from_epoch_day(ed, &y, &m, &d) == 0 && y == 1999 && m == 12 && d == 31, "ymd_bits_from_epoch_day");
    check(weekday_from_epoch_day(ed) == 5 && weekday_from_date("2024-03-05") == 2, "weekdays (Fri, Tue)");

    int err = 0;
    char * s = ymd_extended_from_epoch_day(ed, 0, &err);

    check(err == 0 && strcmp(s, "1999-12-31") == 0, "ymd_extended_from_epoch_day");
    free(s);
    check(ymd_basic_from_epoch_day(ed, 0, &err) == 19991231.0 && epoch_day_from_ymd_basic(19991231) == ed, "basic format round trip");

    /* the Julian calendar: 1582-10-05 (Julian) is 1582-10-15 (Gregorian) */
    uint ej = epoch_day_from_julian_ymd(1582, 10, 5);

    check(ej == epoch_day_from_ymd(1582, 10, 15), "Julian/Gregorian correspondence");
    check(julian_ymd_bits_from_epoch_day(ej, &y, &m, &d) == 0 && y == 1582 && m == 10 && d == 5, "Julian bits round trip");
    check(nearby_epoch_day(2024, 3, 3, 5) == epoch_day_from_ymd(2024, 3, 4), "nearby_epoch_day skips a Sunday for 5-day weeks");
    check(nearby_epoch_day(2024, 3, 2, 5) == epoch_day_from_ymd(2024, 3, 4), "... and a Saturday");
}

void test_dataset_calendar()
{
    DatasetRaw * dset = create_new_dataset(2, 12, 0);

    dset.pd = 5;
    dset.structure = TsCodes.time_series;
    dset.sd0 = epoch_day_from_ymd(2024, 3, 4);   /* a Monday */
    strcpy(dset.stobs.ptr, "2024-03-04");

    char[16] buf = 0;

    check(calendar_date_string(buf.ptr, 0, dset) == 0 && strcmp(buf.ptr, "2024-03-04") == 0, "calendar_date_string t=0");
    check(calendar_date_string(buf.ptr, 5, dset) == 0 && strcmp(buf.ptr, "2024-03-11") == 0, "5-day week: obs 5 is the next Monday");
    check(calendar_obs_number("2024-03-11", dset, 0) == 5 && calendar_obs_number("2024-03-04", dset, 0) == 0, "calendar_obs_number (5-day)");
    check(calendar_obs_number("2024-03-09", dset, 0) == -1, "a Saturday is not in a 5-day calendar");
    check(epoch_day_from_t(7, dset) == epoch_day_from_ymd(2024, 3, 13), "epoch_day_from_t");

    dset.pd = 7;
    check(calendar_obs_number("2024-03-11", dset, 0) == 7 && calendar_date_string(buf.ptr, 7, dset) == 0 &&
          strcmp(buf.ptr, "2024-03-11") == 0, "7-day calendar");
    dset.pd = 52;
    check(calendar_obs_number("2024-03-25", dset, 0) == 3, "weekly calendar");
    destroy_dataset(dset);

    check(get_days_in_month(2, 2024, 7, 0) == 29 && get_days_in_month(2, 2023, 7, 0) == 28, "get_days_in_month");
    check(get_days_in_month(3, 2024, 5, 0) == 21 && get_days_in_month(3, 2024, 6, 0) == 26, "relevant days in March 2024 (5- and 6-day weeks)");
    check(day_starts_month(1, 4, 2024, 7, null) == 1 && day_starts_month(1, 6, 2024, 5, null) == 0 && day_starts_month(3, 6, 2024, 5, null) == 1,
          "day_starts_month (June 1, 2024 is a Saturday)");
    check(day_ends_month(31, 3, 2024, 7) == 1 && day_ends_month(29, 3, 2024, 5) == 1, "day_ends_month");
    check(month_day_index(2024, 3, 11, 5) == 7 && days_in_month_before(2024, 3, 11, 5) == 6, "month_day_index / days_in_month_before");
    char[16] ds = 0;
    check(date_to_daily_index("2024-03-11", 5) == 6 && daily_index_to_date(ds.ptr, 2024, 3, 6, 5) == 0 && strcmp(ds.ptr, "2024-03-11") == 0,
          "date_to_daily_index / daily_index_to_date");
    int err = 0;

    check(day_span(epoch_day_from_ymd(2024, 3, 1), epoch_day_from_ymd(2024, 3, 31), 5, &err) == 21, "day_span (5-day)");
}

void test_misc()
{
    char[16] buf = 0;

    MS_excel_date_string(buf.ptr, 45000, 0, 0);
    check(strcmp(buf.ptr, "2023-03-15") == 0, "MS_excel_date_string (serial 45000)");
    MS_excel_date_string(buf.ptr, 59, 0, 0);
    check(strcmp(buf.ptr, "1900-02-28") == 0, "... 59 is 1900-02-28");
    MS_excel_date_string(buf.ptr, 45000, 4, 0);
    check(strcmp(buf.ptr, "2023:1") == 0, "... quarterly");
    MS_excel_date_string(buf.ptr, 45000, 12, 0);
    check(strcmp(buf.ptr, "2023:03") == 0, "... monthly");

    check(iso_week_from_date("2021-01-03") == 53 && iso_week_from_date("2018-12-31") == 1, "ISO week numbers");
    check(easterdate(2024) == 3.31 && easterdate(2025) == 4.2 && easterdate(2000) == 4.23, "easterdate");
    check(get_dec_date("2024-07-01") > 2024.49 && get_dec_date("2024-07-01") < 2024.51 && na(get_dec_date("nonsense")), "get_dec_date");
    check(legacy_day_of_week(2024, 3, 3, 0, null) == 0.0 && legacy_day_of_week(2024, 3, 5, 0, null) == 2.0, "legacy_day_of_week (Sunday = 0)");

    double[3] b = [20240305, NADBL, 20230229];
    double[3] yy, mm, dd;

    iso_basic_to_extended(b.ptr, yy.ptr, mm.ptr, dd.ptr, 3);
    check(yy[0] == 2024 && mm[0] == 3 && dd[0] == 5 && na(yy[1]) && na(yy[2]), "iso_basic_to_extended (valid, NA, invalid date)");

    double[3] ml, mov = [1, 2, 4];
    check(fill_monthlen_array(ml.ptr, 0, 2, 7, 0, 2024, mov.ptr, null, 0) == 0 && ml[0] == 31 && ml[1] == 31 + 0 - 2 && ml[2] == 30, "fill_monthlen_array");

    char[40] out_ = 0;

    check(gretl_strfdate(out_.ptr, 40, "%d/%m/%Y", epoch_day_from_ymd(2024, 3, 5)) == 10 && strcmp(out_.ptr, "05/03/2024") == 0, "gretl_strfdate");
    check(gretl_alt_strfdate(out_.ptr, 40, 1, epoch_day_from_julian_ymd(1500, 2, 29)) == 10 && strcmp(out_.ptr, "1500-02-29") == 0, "gretl_alt_strfdate (Julian leap day)");

    /* Unix time 1700000000 = 2023-11-14 22:13:20 UTC */
    check(gretl_strftime(out_.ptr, 40, "%Y-%m-%d %H:%M:%S", 1700000000, 0.0) > 0 && strcmp(out_.ptr, "2023-11-14 22:13:20") == 0, "gretl_strftime (UTC)");
    check(gretl_strftime(out_.ptr, 40, "%H:%M", 1700000000, 3600.0) > 0 && strcmp(out_.ptr, "23:13") == 0, "gretl_strftime (UTC+1)");
    check(gretl_strftime(out_.ptr, 40, "8601", 1700000000, 0.0) > 0 && strcmp(out_.ptr, "2023-11-14T22:13:20Z") == 0, "gretl_strftime (ISO 8601, UTC)");
    check(gretl_strftime(out_.ptr, 40, "8601", 1700000000, 5400.0) > 0 && strcmp(out_.ptr, "2023-11-14T23:43:20+01:30") == 0, "gretl_strftime (ISO 8601, +01:30)");

    double dt;
    char * rem = gretl_strptime("2023-11-14 22:13:20+01:00", "%Y-%m-%d %H:%M:%S%z", &dt);

    check(rem !is null && dt == 1700000000.0 - 3600, "gretl_strptime with a time-zone offset");
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);
    test_epoch();
    test_dataset_calendar();
    test_misc();
    printf("\n%d checks, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
