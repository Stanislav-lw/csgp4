#include "afuncs.h"

#include <math.h>

#include "norad_in.h"

#define kF (1.0 / 298.26)

static double c_mod(double x, double y)
{
    if (y == 0.0) {
        return x;
    }
    return (x - y * floor(x / y));
}

static double wrap_neg_pos_pI(const double a)
{
    return c_mod(a + pi, twopi) - pi;
}


static double arctan(double sinx, double cosx)
{
    if (cosx == 0.0) {
        if (sinx > 0.0) {
            return (pi/2.0);
        } else {
            return (3.0 * pi/2.0);
        }
    } else {
        if (cosx > 0.0) {
            return atan(sinx/cosx);
        } else {
            return (pi + atan(sinx/cosx));
        }
    }
}



double get_magnitude(double x, double y, double z)
{
    const double magnitude = sqrt(x*x+y*y+z*z);
    return magnitude;
}

double get_latitude(double x, double y, double z)
{
    const double r = sqrt((x * x) + (y * y));
    const double lat = arctan(z, r);
    return lat;
}

double get_longitude(double x, double y, double gst)
{
    const double theta = arctan(y, x);
    const double lon = wrap_neg_pos_pI(theta - gst);
    return lon;
}

double get_altitude(double x, double y, double z)
{
    static const double e2 = kF * (2.0 - kF);
    const double r = sqrt((x * x) + (y * y));
    double lat = arctan(z, r);
    double phi = 0.0;
    double c = 0.0;
    int cnt = 0;
    do {
        phi = lat;
        const double sinphi = sin(phi);
        c = 1.0 / sqrt(1.0 - e2 * sinphi * sinphi);
        lat = arctan(z + earth_radius_in_km * c * e2 * sinphi, r);
        cnt++;
    } while (fabs(lat - phi) >= 1e-10 && cnt < 10);
    const double alt = r / cos(lat) - earth_radius_in_km * c;
    return alt;
}

coordinat_geodetic to_geodetic(const double *position, double gst)
{
    static const double e2 = kF * (2.0 - kF);
    const double theta = arctan(position[1], position[0]);
    const double lon = wrap_neg_pos_pI(theta - gst);
    const double r = sqrt((position[0] * position[0]) + (position[1] * position[1]));
    double lat = arctan(position[2], r);
    double phi = 0.0;
    double c = 0.0;
    int cnt = 0;
    do {
        phi = lat;
        const double sinphi = sin(phi);
        c = 1.0 / sqrt(1.0 - e2 * sinphi * sinphi);
        lat = arctan(position[2] + earth_radius_in_km * c * e2 * sinphi, r);
        cnt++;
    } while (fabs(lat - phi) >= 1e-10 && cnt < 10);
    const double alt = r / cos(lat) - earth_radius_in_km * c;
    coordinat_geodetic result;
    result.altitude = alt;
    result.longitude = lon;
    result.latitude = lat;
    return result;
}
