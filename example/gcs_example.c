#include <stdio.h>
#include <signal.h>

#include "csgp4/afuncs.h"
#include "csgp4/norad.h"
#include "csgp4/datetime.h"

#define pi 3.141592653589793238462643383279502884197

double to_radiant(double degrees)
{
    return degrees*(pi/180.);
}

double to_degress(double radiant)
{
    return radiant*(180./pi);
}



volatile sig_atomic_t stop = 0;

void inthand(int signum)
{
    stop = signum;
}



int main(int argc, char** argv)
{
    if (argc <= 1) {
        return -1;
    }
    signal(SIGINT, inthand);
    char* filename = argv[1];
    FILE *file = fopen(filename, "r");
    char line1[100], line2[100];
    tle_t tle;
    if (file) {
        fgets(line1, sizeof(line1), file);
        fgets(line2, sizeof(line2), file);
        int parse_code = parse_elements(line1, line2, &tle);
        if (parse_code == 0) {
            double sgp_sat_params[N_SGP4_PARAMS];
            SGP4_init(sgp_sat_params, &tle);
            while (!stop) {
                double state_vector[6];
                long int current_ticks = get_system_ticks();
                SGP4(tsince(tle.epoch, current_ticks), &tle, sgp_sat_params, state_vector, state_vector+3);
                double gst = to_greenwich_sidereal_time(current_ticks);
                double magnitude = get_magnitude(state_vector[0], state_vector[1], state_vector[2]);
                coordinat_geodetic current_geodetic = to_geodetic(state_vector, gst);
                double longitude = to_degress(current_geodetic.longitude);
                double latitude = to_degress(current_geodetic.latitude);
                double altitude = current_geodetic.altitude;
                printf("Geodetic posistion: latitude %f, longitude %f, altitude %f, magnitude %f\n",
                       latitude, longitude, altitude, magnitude);
            }
        }
        fclose(file);
        return parse_code;
    } else {
        printf( "Couldn't open input TLE file\n");
        return -1;
    }


}
