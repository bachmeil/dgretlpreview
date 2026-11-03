/* Port of dataset.h */

module gretl.dataseth;


import core.stdc.string : strchr;
import gretl.libgretlh;

/* opaque: holds the string values of a string-valued series */
public import gretl.gretl_string_table : SeriesTable;

enum DatasetMarkerType {
    no_markers = 0,
    regular_markers,
    daily_date_strings
}

enum VarFlags {
    discrete   = 1 << 0,
    hidden     = 1 << 1,
    generated  = 1 << 2,
    listarg    = 1 << 3,
    timecol    = 1 << 4,
    hfanchor   = 1 << 5,
    coded      = 1 << 6
}

enum DatasetOp {
    none,
    addobs,
    compact,
    expand,
    transpose,
    sortby,
    dsortby,
    resample,
    clear,
    renumber,
    insobs,
    pad_daily,
    unpad_daily
}

enum DataCopyFlag {
    copy_values,
    grab_values
}

/* Symbolic codes for various methods of compacting data
   series (i.e. converting from a higher to a lower
   frequency). wday is applicable only when
   converting from daily to weekly frequency. */
enum CompactMethod {
    unset,   /* no data compaction method is set */
    sum,     /* take sum of higher frequency data */
    avg,     /* take mean of higher frequency data */
    sop,     /* use start-of-period value */
    eop,     /* use end-of-period value */
    wday,    /* use a specified day of the week */
    spread,  /* spread out into multiple series */
    max      /* sentinel value */
}

bool dataset_is_cross_section(const(DatasetRaw) * p) {
    return (p !is null && p.structure == TsCodes.cross_section);
}

bool dataset_is_time_series(const(DatasetRaw) * p) {
    return (p !is null &&
            (p.structure == TsCodes.time_series ||
             p.structure == TsCodes.special_time_series));
}

bool dataset_is_seasonal(const(DatasetRaw) * p) {
    return (p !is null &&
            (p.structure == TsCodes.time_series ||
             p.structure == TsCodes.special_time_series) &&
            p.pd > 1);
}

bool custom_time_series(const(DatasetRaw) * p) {
    return (p !is null && p.structure == TsCodes.special_time_series);
}

bool dataset_is_daily(const(DatasetRaw) * p) {
    return (p !is null &&
            p.structure == TsCodes.time_series &&
            (p.pd == 5 ||
             p.pd == 6 ||
             p.pd == 7));
}

bool dataset_is_incomplete_daily(const(DatasetRaw) * p) {
    return (p !is null &&
            p.structure == TsCodes.time_series &&
            (p.pd == 5 ||
             p.pd == 6 ||
             p.pd == 7) &&
            p.markers == DatasetMarkerType.daily_date_strings);
}

bool dataset_is_weekly(const(DatasetRaw) * p) {
    return (p !is null && p.structure == TsCodes.time_series && p.pd == 52);
}

bool dataset_is_hourly(const(DatasetRaw) * p) {
    return (p !is null && p.structure == TsCodes.time_series && p.pd == 24);
}

bool dataset_is_decennial(const(DatasetRaw) * p) {
    return (p !is null && p.structure == TsCodes.time_series && p.pd == 10);
}

bool dated_daily_data(const(DatasetRaw) * p) {
    return (p !is null &&
            p.structure == TsCodes.time_series &&
            (p.pd == 5 ||
             p.pd == 6 ||
             p.pd == 7) &&
            p.sd0 > 100000);
}

bool undated_daily_data(const(DatasetRaw) * p) {
    return (p !is null &&
            p.structure == TsCodes.time_series &&
            (p.pd == 5 ||
             p.pd == 6 ||
             p.pd == 7) &&
            p.sd0 == 1);
}

bool dated_weekly_data(const(DatasetRaw) * p) {
    return (p !is null && p.structure == TsCodes.time_series && p.pd == 52 && p.sd0 > 100000);
}

bool calendar_data(const(DatasetRaw) * p) {
    return (p !is null &&
            p.structure == TsCodes.time_series &&
            (p.pd == 5 ||
             p.pd == 6 ||
             p.pd == 7 ||
             p.pd == 52) &&
            (p.sd0 > 100000 ||
             (strchr(p.stobs.ptr, '-') !is null)));
}

bool quarterly_or_monthly(const(DatasetRaw) * p) {
    return (p !is null && p.structure == TsCodes.time_series && (p.pd == 4 || p.pd == 12));
}

bool annual_data(const(DatasetRaw) * p) {
    return (p !is null && p.structure == TsCodes.time_series && p.pd == 1);
}

bool decennial_data(const(DatasetRaw) * p) {
    return (p !is null && p.structure == TsCodes.time_series && p.pd == 10 && p.sd0 > 1000);
}

bool dataset_is_panel(const(DatasetRaw) * p) {
    return (p !is null && p.structure == TsCodes.stacked_time_series);
}

bool dataset_is_seasonal_panel(const(DatasetRaw) * p) {
    return (p !is null && p.structure == TsCodes.stacked_time_series && p.panel_pd > 1);
}

bool dataset_has_markers(const(DatasetRaw) * p) {
    return (p !is null && p.markers != 0 && p.S !is null);
}

bool dataset_has_panel_time(const(DatasetRaw) * p) {
    return (p !is null &&
            p.structure == TsCodes.stacked_time_series &&
            p.panel_pd > 0 &&
            p.panel_sd0 > 0.0);
}

int sample_size(const(DatasetRaw) * p) {
    return (p is null) ? 0 : (p.t2 - p.t1 + 1);
}

double dset_get_data(DatasetRaw * d, int i, int t) {
    return d.Z[i][t];
}

void dset_set_data(DatasetRaw * d, int i, int t, double x) {
    d.Z[i][t] = x;
}
