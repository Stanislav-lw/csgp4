/* Copyright (C) 2018, Project Pluto.  See LICENSE.  */

#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

#include "norad.h"
#include "norad_in.h"
#include "datetime.h"

#define minutes_per_day_squared (minutes_per_day * minutes_per_day)
#define minutes_per_day_cubed (minutes_per_day * minutes_per_day_squared)

/* TLEs have four angles on line 2,  given in the form DDD.DDDD.  This
   can be parsed more quickly as an integer,  then cast to double and
   converted to radians,  all in one step.    */

static int get_angle(const char *buff)
{
   int rval = 0;
   while (*buff == ' ') {
      buff++;
   }
   while (*buff != ' ') {
      if (*buff != '.') {
         rval = rval * 10 + (int)(*buff - '0');
      }
      buff++;
   }
   return rval;
}

/* Converts the quasi scientific notation of the "Motion Dot Dot/6" or
"BSTAR" field to double.  The input will always be of the form

sdddddSe

   ....where s is blank or + or -;  ddddd is a five-digit mantissa;
S is + or - or blank;  and e is a single-digit exponent.  A decimal
point is assumed before the five-digit mantissa.  */

static double sci(const char *string)
{
   double rval = 0.;
   if (string[1] != ' ') {
       const int ival = atoi(string);
       if (ival) {
           rval = (double)ival * 1.e-5;
           if(string[7] != '0') {
               int exponent = string[7] - '0';
               if( string[6] == '-') {
                   while(exponent--) {
                       rval *= .1;
                   }
               } else {
                   while(exponent--) {
                       rval *= 10.;
                   }
               }
           }
       }
   }
   return rval;
}


/* Does a checksum modulo 10 on the given line.  Digits = their
value, '-' = 1, all other chars = 0.  Returns 0 if ok, a negative
value if it's definitely not a TLE line,  positive if it's all OK
except the checksum.  This last was added because people sometimes
want to use TLEs without worrying about the checksum. */

static int tle_checksum(const char *buff)
{
    int rval = 0;
    int count = 69;
    if ((*buff != '1' && *buff != '2') || buff[1] != ' ') {
        return -1;
    }
    while (--count) {
        if (*buff > '0' && *buff <= '9') {
            rval += *buff - '0';
        } else if ( *buff == '-') {
            rval++;
        }
        if (*buff < ' ' || *buff > 'z') {           /* invalid character */
            return -2;
        }
        buff++;
    }
    rval -= *buff++ - '0';
    if (*buff > ' ') {                 /* line unterminated */
        rval = -3;
    } else {
        rval %= 10;
        if (rval < 0) {
            rval += 10;
        }
    }
    return (rval % 10);
}

static int mutant_dehex(char ichar)
{
    int rval;
    if (ichar <= '9' && ichar >= '0') {
        rval = ichar - '0';
    } else if (ichar >= 'A' && ichar <= 'Z') {
        rval = ichar + 10 - 'A';
    } else {
        rval = -1;
    }
    return rval;
}

/* The "standard" SDP4 model fails badly for very high-flying satellites.
As a non-standard extension,  I'm simply storing state vectors for them,
using the following somewhat odd scheme :

1 40391U 15007B   15091.99922241 sxxxxxxxx syyyyyyyy szzzzzzzzH  9997
2 49391 [valid range, accuracy]  saaaaaaaa sbbbbbbbb scccccccc    0 8

   Epoch,  int'l & NORAD IDs are stored in the standard manner.  The
'ephemeris type' is H (rather than the otherwise universal 0).  The
xyz position and vx, vy, vz velocity are stored as 8-digit signed
base-36 integers,  hence a range of +/- 36^8 = about +/- 2.82x10^12.

  x, y, z are in meters,  and hence cover a range +/- 18.9 AU.
vx, vy, vz are in 10^-4 m/s,  range +/- 94% c.  */

static double get_high_value(const char *iptr)
{
   int64_t rval = 0;
   if(*iptr == '+' || *iptr == '-') {
       int i, digit;
       for( i = 1; i < 9; i++) {
           digit = mutant_dehex(iptr[i]);
           rval = rval * (int64_t)36LL + (int64_t)digit;
       }
       if (*iptr == '-') {
           rval = -rval;
       }
   }
   return ((double)rval);
}

/* Traditionally,  NORAD numbers were stored as five digits.  In 2020, new
detectors threatened to go past 100K objects;  the 'Alpha-5' scheme allows
the first byte to be replaced by an uppercase letter,  with I and O
skipped.  That gets us to 339999 :

https://www.space-track.org/documentation#tle-alpha5

   Note that Alpha-5 is referred to as a "stopgap".  Near the bottom of
the above link,  "space-track.org encourages users to switch to... XML,
KVN,  or JSON",  (partly) because these will handle nine-digit catalog
numbers.

   I have implemented the following (unofficial,  my proposal) scheme to
go beyond the Alpha-5 limit of 340000 possible numbers.  To do so,  we
use all 34^5 = 45435424 possible combinations;  i.e.,  each of the five
characters can be a digit or a letter.  For lack of a better name,  call
it 'Super-5'.  It does not get us to a full nine digits,  is only
supported by me,  and is not as easy to parse visually as Alpha-5.

   We also could use lowercase letters and some others to get 10^9
combinations within five bytes with backward compatibility.  (Which may
eventually be needed;  it sounds as if Space-Track may make some use of
the full nine-digit range.)

   d = digit, L = letter,  x = either.

(1) ddddd = 'traditional' scheme provides 100000 combinations;
         Numbers 0 to 99999

(2) Ldddd = Alpha-5 scheme adds 240000
         Numbers 100000 to 339999;     A0000 to Z9999

(3) xxxxL = 34^4*24      = 32072064 more  (start of 'Super-5' range)
         Numbers 340000 to 32412063;   0000A to ZZZZZ

(4) xxxLd = 34^3*24*10   =  9432960 more
         Numbers 32412064 to 41845023; 000A0 to ZZZZ9

(5) xxLdd = 34^2*24*100  =  2774400 more
         Numbers 41845024 to 44619423; 00A00 to ZZZ99

(6) xLddd = 34*24*1000   =   816000 more  (end of 'Super-5' range)
         Numbers 44619424 to 45435423; 0A000 to ZZ999     */

static int base34_to_int(char c)
{
    int offset;
    if (c == ' ') {
        return 0;
    }
    if (c >= '0' && c <= '9') {
        offset = '0';
    } else if (c >= 'A' && c <= 'H') {
        offset = 'A' - 10;
    } else if (c >= 'J' && c <= 'N') {
        offset = 'J' - 18;
    } else if (c >= 'P' && c <= 'Z') {
        offset = 'P' - 23;
    } else {
        return -1;
    }
    return (c - offset);
}

static int get_norad_number(const char *buff)
{
   size_t i;
   int digits[5], rval = 0;
   for(i = 0; i < 5; i++) {
       digits[i] = base34_to_int(buff[i]);
       if(digits[i] == -1) {       /* not a valid number */
           return 0;
       }
   }
   if (digits[4] > 9) {      /* case (3): last char is uppercase */
      rval = 340000 + (digits[4] - 10) + 24 * (digits[3]
               + digits[2] * 34 + digits[1] * 34 * 34 + digits[0] * 34 * 34 * 34);
   } else if (digits[3] > 9) {    /* case (4) above */
      rval = 340000 + 32072064 + digits[4] + (digits[3] - 10) * 10
            + 10 * 24 * (digits[2] + digits[1] * 34 + digits[0] * 34 * 34);
   } else if (digits[2] > 9) {    /* case (5) above */
      rval = 340000 + 32072064 + 9432960 + digits[4]
            + digits[3] * 10 + (digits[2] - 10) * 100
            + 2400 * (digits[1] + digits[0] * 34);
   } else if (digits[1] > 9) {   /* case (6) above */
      rval = 340000 + 32072064 + 9432960 + 2774400 + digits[4]
            + digits[3] * 10 + digits[2] * 100 + (digits[1] - 10) * 1000
            + digits[0] * 24000;
   } else {       /* last four digits are 0-9;  'standard' NORAD desig */
      rval = digits[0] * 10000 + atoi(buff + 1);
   }
   return rval;
}

static double get_eight_places(const char *ptr)
{
   int lval = atoi(ptr);
   int rval = atoi(ptr + 4);
   return ((double)lval + (double)rval * 1e-8);
}

/* parse_elements returns:
   0 if the elements are parsed without error;
   1 if they're OK except the first line has a checksum error;
   2 if they're OK except the second line has a checksum error;
   3 if they're OK except both lines have checksum errors;
   a negative value if the lines aren't at all parseable */

int parse_elements(const char *line1, const char *line2, tle_t *tle)
{
    int rval = 0;
    int checksum_problem = 0;
    if(*line1 != '1' || *line2 != '2') {
        rval = -4;
    } else {
        rval = tle_checksum(line1);
        if (rval > 0) {
            checksum_problem = 1;  /* there's a checksum problem,  but it's */
            rval = 0;              /* not fatal; continue processing the TLE */
        }
    }
    if (rval) {
        rval -= 100;
    } else {
        rval = tle_checksum(line2);
        if(rval > 0) {
            checksum_problem |= 2;  /* there's a checksum problem,  but it's */
            rval = 0;               /* not fatal; continue processing the TLE */
        }
    }
   if (!rval) {
      char tbuff[13];
      int year = line1[19] - '0';
      if (line1[18] >= '0') {
         year += (line1[18] - '0') * 10;
      }
      if (year < 57) {          /* cycle around Y2K */
         year += 2000;
      } else {
         year += 1900;
      }
      double days = get_eight_places(line1 + 20);
      tle->epoch = get_ticks_from_yd(year, days);
      tle->norad_number = get_norad_number(line1 + 2);
      memcpy(tbuff, line1 + 64, 4);
      tbuff[4] = '\0';
      tle->bulletin_number = atoi(tbuff);
      tle->classification = line1[7];       /* almost always 'U' */
      memcpy(tle->intl_designation, line1 + 9, 8);
      tle->intl_designation[8] = '\0';
      memcpy(tbuff, line2 + 63, 5);
      tbuff[5] = '\0';
      tle->revolution_number = atoi(tbuff);
      tle->ephemeris_type = line1[62];
      if (tle->ephemeris_type == 'H') {
         double *state_vect = &tle->xincl;
         for(int i = 0; i < 3; i++) {
            state_vect[i] = get_high_value(line1 + 33 + i * 10);
            state_vect[i + 3] = get_high_value(line2 + 33 + i * 10) * 1e-4;
         }
         return 0;
      }
      int mean_anomaly = get_angle(line2 + 43);
      int xnodeo_a = get_angle(line2 + 17);
      int omegao_a = get_angle(line2 + 34);
      int inclination = get_angle(line2 + 8);
      tle->xmo = (double)mean_anomaly * (pi / 180e+4);
      tle->xnodeo = (double)xnodeo_a * (pi / 180e+4);
      tle->omegao = (double)omegao_a * (pi / 180e+4);
      tle->xincl = (double)inclination * (pi / 180e+4);
      tle->eo = atoi(line2 + 26) * 1.e-7;
      /* Make sure mean motion is null-terminated, since rev. no. may immediately follow. */
      memcpy(tbuff, line2 + 51, 12);
      tbuff[12] = '\0';
      /* Input mean motion,  derivative of mean motion and second
         deriv of mean motion,  are all in revolutions and days.
         Convert them here to radians and minutes:                 */
      tle->xno = get_eight_places(tbuff) * twopi / minutes_per_day;
      tle->bstar = sci(line1 + 53) * ae;
   }
   return (rval ? rval : checksum_problem);
}
