#ifndef AFUNCS
#define AFUNCS

typedef struct {
    double longitude;
    double latitude;
    double altitude;
} coordinat_geodetic;

#ifdef __cplusplus
extern "C" {
#endif
double get_magnitude(double x, double y, double z);
double get_latitude(double x, double y, double z);
double get_longitude(double x, double y, double gst);
double get_altitude(double x, double y, double z);
coordinat_geodetic to_geodetic(const double* position, double gst);
#ifdef __cplusplus
}                       /* end of 'extern "C"' section */
#endif

#endif  /* AFUNCS */
