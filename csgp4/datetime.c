#include "datetime.h"

#include <stdint.h>
#include <sys/time.h>
#include <math.h>

#include "norad_in.h"

#define to_radians(degrees)(degrees * pi/180.0)

double absolute_days_yd(int year, double days)
{
    long int previousYear = year - 1;
    /*
     * + days in previous years ignoring leap days
     * + Julian leap days before this year
     * - minus prior century years
     * + plus prior years divisible by 400 days
     */
    long int daysSoFar = 365L * previousYear
                        + previousYear / 4L
                        - previousYear / 100L
                        + previousYear / 400L;

    return ((double)daysSoFar + days - 1.0);
}


long get_system_ticks(void)
{
    struct timeval tv;
    gettimeofday(&tv, 0);
    const long int system_ticks = UNIX_EPOCH + tv.tv_sec*1000000 + tv.tv_usec;
    return system_ticks;
}

long int get_ticks_from_yd(int year, double days)
{
    double absolute_days = absolute_days_yd(year, days);
    long int ticks = (long int)(absolute_days * TICKS_PER_DAY);
    return ticks;
}

double to_julian(long int ticks)
{
    return (((double)ticks/TICKS_PER_DAY) + 1721425.5);
}

double to_greenwich_sidereal_time(long int ticks)
{
    // julian date of previous midnight
    double jd = to_julian(ticks);
    double jd0 = floor(jd + 0.5) - 0.5;
    // julian centuries since epoch
    double t   = (jd0 - 2451545.0) / 36525.0;
    double jdf = jd - jd0;

    double gt  = 24110.54841 + t * (8640184.812866 + t * (0.093104 - t * 6.2E-6));
    gt  += jdf * 1.00273790935 * 86400.0;

    // 360.0 / 86400.0 = 1.0 / 240.0
    double result = to_radians(gt/240.0) - twopi * floor(to_radians(gt/240.0) / twopi);
    return result;
}

double to_j2000(long int ticks)
{
    return (to_julian(ticks) - 2415020.0);
}

double tsince(long t1, long t2)
{
    long diff = t2 - t1;
    double result = (double)diff / TICKS_PER_MINUTE;
    return result;
}
