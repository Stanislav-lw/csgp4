/* Copyright (C) 2018, Project Pluto.  See LICENSE.  */

#ifndef NORAD_IN_H
#define NORAD_IN_H

/* Common "internal" arguments between deep-space functions;  users of
   the satellite routines shouldn't need to bother with any of this */

typedef struct
{
  double
  /* Common between SGP4 and SDP4: */
  aodp, cosio, sinio, omgdot, xmdot, xnodot, xnodp,
  /* Used by dpinit part of Deep() */
  eosq, betao, cosio2, sing, betao2;
} deep_arg_t;

typedef struct
{
   double coef, coef1, tsi, s4, unused_a3ovk2, eta;
} init_t;

/* Table of constant values */
#define pi                  3.141592653589793238462643383279502884197
#define twopi               (pi*2.0)
#define e6a                 1.0E-6
#define two_thirds          (2.0 / 3.0)
#define xj3                 -2.53881E-6
#define minus_xj3           2.53881E-6
#define earth_radius_in_km  6378.135
#ifndef minutes_per_day
   #define minutes_per_day  1440.0
#endif
#define ae                  1.0
#define xj2                 1.082616e-3
#define ck2                 (0.5 * xj2 * ae * ae)

#define xj4                 (-1.65597e-6)
#define ck4                 (-0.375 * xj4 * ae * ae * ae * ae)
#define s_const             (ae * (1.0 + 78.0 / earth_radius_in_km))
#define qoms2t              1.880279159015270643865e-9
#define xke                 0.0743669161331734132

#define a3ovk2              (minus_xj3/ck2*ae*ae*ae)

#endif /* NORAD_IN_H */
