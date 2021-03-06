/**
 * \author Sylvain Marsat, University of Maryland - NASA GSFC
 *
 * \brief C header for the conversion between time scales.
 *
 */

#ifndef _TIMECONVERSION_H
#define _TIMECONVERSION_H

#define _XOPEN_SOURCE 500

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <unistd.h>
#include <getopt.h>
#include <stdbool.h>
#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_spline.h>

#include "constants.h"

/*************************/
/****** Prototypes ******/

/* Function computing the gmst angle from the gps time */
double gmst_angle_from_gpstime(const double gpstime); /* gpstime in seconds */

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _TIMECONVERSION_H */
