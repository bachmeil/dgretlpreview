/* Port of calendar.h */

module gretl.calendarh;


/*
 * FOUR_DIGIT_YEAR:
 * @y: year given to only two digits.
 *
 * Returns: a guess at the 4-digit year intended when
 * two digits were provided. But really, who knows?
 */
int FOUR_DIGIT_YEAR(int y)
{
    return (y < 50) ? y + 2000 : y + 1900;
}
