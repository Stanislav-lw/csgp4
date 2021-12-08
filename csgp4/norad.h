/* Copyright (C) 2018, Project Pluto.  See LICENSE.  */
#ifndef NORAD_H
#define NORAD_H

/* Two-line-element satellite orbital data */
typedef struct {
    long int epoch;
    double xincl, xnodeo, eo, omegao, xmo, xno, bstar;
    int norad_number, bulletin_number, revolution_number;
    char classification;
    char ephemeris_type;
    char intl_designation[9];
} tle_t;

/* Table of abbreviations */
/* xnodeo - ascending node,
   omegao - argument perigee,
   eo     - eccentricity,
   xincl  - inclination,
   aodp   - recovered semi major axis,
   xnodp  - recovered mean motion,
   xmo    - mean anomaly (radians),
   xno    - mean motion  (radians/minute),
   epoch,                (microseconds) */


#define N_SGP4_PARAMS         24

/* SGP4 can return zero,  or any of the following error/warning codes.
   The 'warnings' result in a mathematically reasonable value being returned,
   and perigee within the earth is completely reasonable for an object that's
   just left the earth or is about to hit it.  The 'errors' mean that no
   reasonable position/velocity was determined. */

#define SXPX_SUCCESS                       0
#define SXPX_WARN_ORBIT_WITHIN_EARTH      -1
#define SXPX_WARN_PERIGEE_WITHIN_EARTH    -2
#define SXPX_ERR_NEARLY_PARABOLIC         -3
#define SXPX_ERR_NEGATIVE_MAJOR_AXIS      -4
#define SXPX_ERR_NEGATIVE_XN              -5
#define SXPX_ERR_CONVERGENCE_FAIL         -6

/* The function parse_elements can return zero, or any of the following error/warning codes.
   The 'warning' results can be returned only if one or both lines has a checksum error.
   If function return a negative value the lines aren't at all parseable */

#define TLE_SUCCESS          0
#define TLE_WARN_FL_CHECKSUM 1
#define TLE_WARN_SL_CHECKSUM 2
#define TLE_WARN_BL_CHECKSUM 3

#ifdef __cplusplus
extern "C" {
#endif
void SGP4_init(double *params, const tle_t *tle);
int SGP4(double tsince, const tle_t *tle, const double *params, double *pos, double *vel);
int parse_elements(const char *line1, const char *line2, tle_t *tle);
#ifdef __cplusplus
}                       /* end of 'extern "C"' section */
#endif

#endif /* NORAD_H */
