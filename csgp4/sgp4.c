/* Copyright (C) 2018, Project Pluto.  See LICENSE.  */

#include <math.h>
#include <stdio.h>
#include "norad.h"
#include "norad_in.h"

#define c2             params[0]
#define c1             params[1]
#define c4             params[2]
#define xnodcf         params[3]
#define t2cof          params[4]
#define p_aodp         params[5]
#define p_cosio        params[6]
#define p_sinio        params[7]
#define p_omgdot       params[8]
#define p_xmdot        params[9]
#define p_xnodot       params[10]
#define p_xnodp        params[11]
#define c5             params[12]
#define d2             params[13]
#define d3             params[14]
#define d4             params[15]
#define delmo          params[16]
#define p_eta          params[17]
#define omgcof         params[18]
#define sinmo          params[19]
#define t3cof          params[20]
#define t4cof          params[21]
#define t5cof          params[22]
#define xmcof          params[23]
#define simple_flag *((int *)(params + 24))

#define MINIMAL_E    1.e-4
#define ECC_EPS      1.e-6     /* Too low for computing further drops. */
#define MAX_KEPLER_ITER 10

void sxpx_common_init(double *params, const tle_t *tle, init_t *init, deep_arg_t *deep_arg);
void sxpall_common_init(const tle_t *tle, deep_arg_t *deep_arg);

void SGP4_init(double *params, const tle_t *tle)
{
    deep_arg_t deep_arg;
    init_t init;

    sxpx_common_init(params, tle, &init, &deep_arg);
    p_aodp = deep_arg.aodp;
    p_cosio = deep_arg.cosio;
    p_sinio = deep_arg.sinio;
    p_omgdot = deep_arg.omgdot;
    p_xmdot = deep_arg.xmdot;
    p_xnodot = deep_arg.xnodot;
    p_xnodp = deep_arg.xnodp;
    p_eta = deep_arg.aodp * tle->eo * init.tsi;

    const double eeta = tle->eo * p_eta;
    /* For perigee less than 220 kilometers, the "simple" flag is set
       and the equations are truncated to linear variation in sqrt a
       and quadratic variation in mean anomaly.  Also, the c3 term,
       the delta omega term, and the delta m term are dropped. */
    simple_flag = ((p_aodp*(1.0 - tle->eo)/ae) < (220.0/earth_radius_in_km+ae));
    if (!simple_flag) {
        const double c1sq = c1 * c1;
        double temp;

        simple_flag = 0.0;
        delmo = 1.0 + p_eta * cos(tle->xmo);
        delmo *= delmo * delmo;
        d2 = 4.0 * p_aodp * init.tsi * c1sq;
        temp = d2 * init.tsi * c1/3;
        d3 = (17.0 * p_aodp + init.s4) * temp;
        d4 = 0.5 * temp * p_aodp * init.tsi * (221.0 * p_aodp + 31.0 * init.s4) * c1;
        t3cof = d2 + 2.0 * c1sq;
        t4cof = 0.25 * (3.0 * d3 + c1 * (12.0 * d2 + 10.0 * c1sq));
        t5cof = 0.2* (3.0 * d4 + 12.0 * c1 * d3 + 6.0 * d2 * d2 + 15.0 * c1sq * (2.0 * d2 + c1sq));
        sinmo = sin(tle->xmo);
        if (tle->eo < MINIMAL_E) {
            omgcof = 0.0;
            xmcof = 0.0;
        } else {
            const double c3 = init.coef * init.tsi * a3ovk2 * p_xnodp * ae * p_sinio/tle->eo;
            xmcof = -two_thirds * init.coef * tle->bstar * ae/eeta;
            omgcof = tle->bstar * c3 * cos(tle->omegao);
        }
    } /* End of if (isFlagClear(SIMPLE_FLAG)) */
    const double etasq = p_eta * p_eta;
    c5 = 2.0 * init.coef1 * p_aodp * deep_arg.betao2 * (1.0 + 2.75 * (etasq + eeta) + eeta * etasq);
}

void sxpx_common_init(double *params, const tle_t *tle, init_t *init, deep_arg_t *deep_arg)
{
    sxpall_common_init(tle, deep_arg);
    const double x3thm1 = 3.0 * deep_arg->cosio2 - 1.0;
    /* For perigee below 156 km, the values
       of s and qoms2t are altered.         */
    init->s4 = s_const;
    double qoms24 = qoms2t;
    const double perige = (deep_arg->aodp * (1.0 - tle->eo) - ae) * earth_radius_in_km;
    if (perige < 156.0) {
        if (perige <= 98.0) {
            init->s4 = 20.0;
        } else {
            init->s4 = perige - 78.0;
        }
        const double temp_val = (120.0 - init->s4) * ae/earth_radius_in_km;
        const double temp_val_squared = temp_val * temp_val;
        qoms24 = temp_val_squared * temp_val_squared;
        init->s4 = init->s4/earth_radius_in_km + ae;
    }  /* End of if(perige <= 156) */
    const double pinv = 1.0/(deep_arg->aodp * deep_arg->betao2);
    const double pinvsq = pinv * pinv;
    init->tsi = 1.0/(deep_arg->aodp - init->s4);
    init->eta = deep_arg->aodp * tle->eo * init->tsi;
    const double etasq = init->eta * init->eta;
    const double eeta = tle->eo * init->eta;
    const double psisq = fabs(1.0 - etasq);
    const double tsi_squared = init->tsi * init->tsi;
    init->coef = qoms24 * tsi_squared * tsi_squared;
    init->coef1 = init->coef/pow(psisq, 3.5);
    c2 = init->coef1 * deep_arg->xnodp * (deep_arg->aodp * (1.0 + 1.5 * etasq + eeta * (4.0 + etasq))
                                          +0.75*ck2*init->tsi/psisq*x3thm1*(8+3*etasq*(8+etasq)));
    c1 = tle->bstar * c2;
    deep_arg->sinio = sin(tle->xincl);
    c4 = 2.0*deep_arg->xnodp*init->coef1*deep_arg->aodp*deep_arg->betao2*
            (init->eta*(2.0+0.5*etasq)+tle->eo*(0.5+2.0*etasq)-2.0*ck2*init->tsi /
             (deep_arg->aodp*psisq)*(-3.0*x3thm1*(1.0-2.0*eeta+etasq*
                                                 (1.5-0.5*eeta))+0.75*(1.0-deep_arg->cosio2) *(2.0*etasq-eeta*(1.0+etasq))*
                                     cos(2.0*tle->omegao))
             );
    const double cosio4 = deep_arg->cosio2 * deep_arg->cosio2;
    const double temp1 = 3.*ck2*pinvsq*deep_arg->xnodp;
    const double temp2 = temp1 * ck2 * pinvsq;
    const double temp3 = 1.25 * ck4 * pinvsq * pinvsq * deep_arg->xnodp;
    deep_arg->xmdot = deep_arg->xnodp
                      + temp1 * deep_arg->betao * x3thm1 / 2.0
                      + temp2 * deep_arg->betao * (13.0 - 78.0 * deep_arg->cosio2 + 137.0 * cosio4)/16.0;
    deep_arg->omgdot = -temp1 * (1.0 - 5.0 * deep_arg->cosio2)/2.0
                       + temp2 * (7.0 - 114.0 * deep_arg->cosio2 + 395.0 * cosio4)/16.0
                       + temp3 * (3.0 - 36.0 * deep_arg->cosio2 + 49.0 * cosio4);
    const double xhdot1 = -temp1 * deep_arg->cosio;
    deep_arg->xnodot = xhdot1 + (temp2 * (4.0 - 19.0 * deep_arg->cosio2) / 2.0
                                 + 2.0 * temp3 * (3.0 - 7.0 * deep_arg->cosio2)) * deep_arg->cosio;
    xnodcf = 3.5 * deep_arg->betao2 * xhdot1 * c1;
    t2cof = 1.5 * c1;
}

void sxpall_common_init(const tle_t *tle, deep_arg_t *deep_arg)
{
    const double a1 = pow(xke/tle->xno, two_thirds);  /* in Earth radii */
    /* Recover original mean motion (xnodp) and
       semimajor axis (aodp) from input elements. */
    deep_arg->cosio = cos(tle->xincl);
    deep_arg->cosio2 = deep_arg->cosio * deep_arg->cosio;
    deep_arg->eosq = tle->eo * tle->eo;
    deep_arg->betao2 = 1.0 - deep_arg->eosq;
    deep_arg->betao = sqrt(deep_arg->betao2);
    const double tval = 1.5 * ck2 * (3.0 * deep_arg->cosio2 - 1.0) / (deep_arg->betao * deep_arg->betao2);
    const double del1 = tval / (a1 * a1);
    const double ao = a1 * (1.0 - del1 * (1.0/3.0 + del1 * (1.0 + 134.0/81.0 * del1)));
    const double delo = tval/(ao * ao);
    deep_arg->xnodp = tle->xno/(1.0 + delo);   /* in radians/minute */
    deep_arg->aodp = ao/(1.0 - delo);
} /* End of SGP4() initialization */



int sxpx_posn_vel(const double xnode, const double a, const double ecc,
                  const double cosio, const double sinio,
                  const double xincl, const double omega,
                  const double xl, double *pos, double *vel);
double centralize_angle(double ival);

int SGP4(double tsince, const tle_t *tle, const double *params, double *pos, double *vel)
{
    /* Update for secular gravity and atmospheric drag. */
    const double xmdf = tle->xmo + p_xmdot * tsince;
    const double omgadf = tle->omegao + p_omgdot * tsince;
    const double xnoddf = tle->xnodeo + p_xnodot * tsince;
    const double tsq = tsince * tsince;
    const double xnode = xnoddf + xnodcf * tsq;
    double omega = omgadf;
    double xmp = xmdf;
    double tempa = 1.0 - c1 * tsince;
    double tempe = tle->bstar * c4 * tsince;
    double templ = t2cof * tsq;
    if (!simple_flag) {
        const double delomg = omgcof * tsince;
        double delm = 1.0 + p_eta * cos(xmdf);

        delm = xmcof * (delm * delm * delm - delmo);
        const double temp = delomg + delm;
        xmp = xmdf + temp;
        omega = omgadf - temp;
        const double tcube = tsq * tsince;
        const double tfour = tsince * tcube;
        tempa = tempa - d2 * tsq - d3 * tcube - d4 * tfour;
        tempe = tempe + tle->bstar * c5 * (sin(xmp) - sinmo);
        templ = templ + t3cof * tcube + tfour * (t4cof + tsince * t5cof);
    }; /* End of if (isFlagClear(SIMPLE_FLAG)) */

    double e = tle->eo - tempe;
    /* A highly arbitrary lower limit on e,  of 1e-6: */
    if (e < ECC_EPS) {
        e = ECC_EPS;
    }
    double a = p_aodp * tempa * tempa;
    if (tempa < 0.0) {      /* force negative a,  to indicate error condition */
        a = -a;
    }
    const double xl = xmp+omega+xnode+p_xnodp*templ;
    return (sxpx_posn_vel(xnode, a, e, p_cosio, p_sinio, tle->xincl, omega, xl, pos, vel));
}

int sxpx_posn_vel(const double xnode, const double a, const double ecc,
                  const double cosio, const double sinio,
                  const double xincl, const double omega,
                  const double xl, double *pos, double *vel)
{
    /* Long period periodics */
    const double axn = ecc*cos(omega);
    const double xlcof = 0.125 * a3ovk2 * sinio * (3.0 + 5.0 * cosio)/(1.0 + cosio);
    const double aycof = 0.25 * a3ovk2 * sinio;
    const double xll = 1.0/(a * (1.0-ecc*ecc)) * xlcof*axn;
    const double aynl = 1.0/(a * (1.0-ecc*ecc)) * aycof;
    const double xlt = xl + xll;
    const double ayn = ecc * sin(omega) + aynl;
    const double elsq = axn * axn + ayn * ayn;
    const double capu = centralize_angle(xlt - xnode);
    const double chicken_factor_on_eccentricity = 1.e-6;

    /* Dundee changes:  items dependent on cosio get recomputed: */
    const double cosio_squared = cosio * cosio;
    const double x3thm1 = 3.0 * cosio_squared - 1.0;
    const double sinio2 = 1.0 - cosio_squared;
    const double x7thm1 = 7.0 * cosio_squared - 1.0;

    /* Extremely decayed satellites can end up "orbiting" within the
       earth.  Eventually,  the semimajor axis becomes zero,
       then negative.  In that case,  or if the orbit is near
       to parabolic,  we zero the posn/vel and quit.  If the
       object has a perigee or apogee indicating a crash,  we
       just flag it. */
    int rval = 0;
    if (a < 0.0) {
        rval = SXPX_ERR_NEGATIVE_MAJOR_AXIS;
    }
    if (elsq > 1.0 - chicken_factor_on_eccentricity) {
        rval = SXPX_ERR_NEARLY_PARABOLIC;
    }
    for (int i = 0; i < 3; i++) {
        pos[i] = 0.0;
        if (vel) {
            vel[i] = 0.0;
        }
    }
    if (rval) {
        return rval;
    }
    if (a * (1.0 - ecc) < 1.0 && a * (1.0 + ecc) < 1.0) {   /* entirely within earth */
        rval = SXPX_WARN_ORBIT_WITHIN_EARTH;     /* remember, e can be negative */
    }
    if (a * (1.0 - ecc) < 1.0 || a * (1.0 + ecc) < 1.0) {  /* perigee within earth */
        rval = SXPX_WARN_PERIGEE_WITHIN_EARTH;
    }
    /* Solve Kepler's' Equation */
    double ecosE, esinE;
    double sinEPW, cosEPW;
    double epw = capu;
    int kepler_iter = 0;
    for(; kepler_iter < MAX_KEPLER_ITER; kepler_iter++) {
        const double newton_raphson_epsilon = 1e-12;
        double f, fdot, delta_epw;
        int do_second_order_newton_raphson = 1;
        sinEPW = sin(epw);
        cosEPW = cos(epw);
        ecosE = axn * cosEPW + ayn * sinEPW;
        esinE = axn * sinEPW - ayn * cosEPW;
        f = capu - epw + esinE;
        if (fabs(f) < newton_raphson_epsilon) {
            break;
        }
        fdot = 1.0 - ecosE;
        delta_epw = f / fdot;
        if (!kepler_iter) {
            const double max_newton_raphson = 1.25 * fabs(ecc);
            do_second_order_newton_raphson = 0;
            if (delta_epw > max_newton_raphson) {
                delta_epw = max_newton_raphson;
            } else if (delta_epw < -max_newton_raphson) {
                delta_epw = -max_newton_raphson;
            } else {
                do_second_order_newton_raphson = 1;
            }
        }
        if (do_second_order_newton_raphson) {
            delta_epw = f / (fdot + 0.5 * esinE * delta_epw);
        }
        epw += delta_epw;
    }
    if (kepler_iter == MAX_KEPLER_ITER) {
        return (SXPX_ERR_CONVERGENCE_FAIL);
    }

    /* Short period preliminary quantities */
    double temp = 1 - elsq;
    const double pl = a * temp;
    double r = a*(1.0 - ecosE);
    double temp2 = a / r;
    const double betal = sqrt(temp);
    temp = esinE/(1 + betal);
    const double cosu = temp2 * (cosEPW - axn + ayn * temp);
    const double sinu = temp2 * (sinEPW - ayn - axn * temp);
    const double u = atan2(sinu, cosu);
    const double sin2u = 2 * sinu * cosu;
    const double cos2u = 2 * cosu * cosu - 1;
    const double temp1 = ck2 / pl;
    temp2 = temp1 / pl;
    /* Update for short periodics */
    const double rk = r * (1 - 1.5 * temp2 * betal * x3thm1) + 0.5 * temp1 * sinio2 * cos2u;
    const double uk = u - 0.25 * temp2 * x7thm1 * sin2u;
    const double xnodek = xnode + 1.5 * temp2 * cosio * sin2u;
    const double xinck = xincl + 1.5 * temp2 * cosio * sinio * cos2u;
    /* Orientation vectors */
    const double sinuk = sin(uk);
    const double cosuk = cos(uk);
    const double sinik = sin(xinck);
    const double cosik = cos(xinck);
    const double sinnok = sin(xnodek);
    const double cosnok = cos(xnodek);
    const double xmx = -sinnok * cosik;
    const double xmy = cosnok * cosik;
    const double ux = xmx * sinuk + cosnok * cosuk;
    const double uy = xmy * sinuk + sinnok * cosuk;
    const double uz = sinik * sinuk;
    /* Position and velocity */
    pos[0] = rk * ux * earth_radius_in_km;
    pos[1] = rk * uy * earth_radius_in_km;
    pos[2] = rk * uz * earth_radius_in_km;
    if (vel) {
        const double rdot = xke * sqrt(a) * esinE/r;
        const double rfdot = xke * sqrt(pl)/r;
        const double xn = xke/(a * sqrt(a));
        const double rdotk = rdot - xn * temp1 * sinio2 * sin2u;
        const double rfdotk = rfdot + xn * temp1 * (sinio2 * cos2u + 1.5 * x3thm1);
        const double vx = xmx * cosuk - cosnok * sinuk;
        const double vy = xmy * cosuk - sinnok * sinuk;
        const double vz = sinik * cosuk;

        vel[0] = (rdotk * ux + rfdotk * vx) * earth_radius_in_km;
        vel[1] = (rdotk * uy + rfdotk * vy) * earth_radius_in_km;
        vel[2] = (rdotk * uz + rfdotk * vz) * earth_radius_in_km;
    }
    return rval;
}

double centralize_angle(double ival)
{
   double rval = fmod(ival, twopi);
   if (rval > pi) {
       rval -= twopi;
   } else if (rval < - pi) {
       rval += twopi;
   }
   return rval;
} /* End of SGP4 */
