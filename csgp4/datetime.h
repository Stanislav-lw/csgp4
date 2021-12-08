#ifndef DATETIME_H
#define DATETIME_H

#define TICKS_PER_DAY 86400000000L
#define TICKS_PER_HOUR 3600000000L
#define TICKS_PER_MINUTE 60000000L
#define TICKS_PER_SECOND 1000000L
#define TICKS_PER_MILLISECOND 1000L
#define TICKS_PER_MICROSECOND 1L
#define UNIX_EPOCH 62135596800000000L
#define MAX_VALUE_TICKS 315537897599999999L
#define GREGORIAN_START 49916304000000000L // 1582-Oct-15

#ifdef __cplusplus
extern "C" {
#endif
double tsince(long int t1, long int t2);
long int get_system_ticks(void);
double absolute_days_yd(int year, double days);
long int get_ticks_from_yd(int year, double days);
double to_julian(long int ticks);
double to_greenwich_sidereal_time(long int ticks);
double to_j2000(long int ticks);
#ifdef __cplusplus
}                       /* end of 'extern "C"' section */
#endif

#endif /* DATETIME_H */
