/**
 * \author Sylvain Marsat, University of Maryland - NASA GSFC
 *
 * \brief C code headers for the geometric coefficients entering the response for LISA-like detectors.
 *
 */

#ifndef _LISAGEOMETRY_H
#define _LISAGEOMETRY_H

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
#include <gsl/gsl_complex.h>

#include "constants.h"

/**************************************************/
/**************** Prototypes **********************/

/* Function cardinal sine */
double sinc(const double x);

/* Function to compute, given a value of a sky position and polarization, all the complicated time-independent trigonometric coefficients entering the response */
void SetCoeffsG(const double lambda, const double beta, const double psi);

/* Functions evaluating the G_AB functions, combining the two polarization with the spherical harmonics factors */
double complex G21mode(const double f, const double t, const double complex Yfactorplus, const double complex Yfactorcross);
double complex G12mode(const double f, const double t, const double complex Yfactorplus, const double complex Yfactorcross);
double complex G32mode(const double f, const double t, const double complex Yfactorplus, const double complex Yfactorcross);
double complex G23mode(const double f, const double t, const double complex Yfactorplus, const double complex Yfactorcross);
double complex G13mode(const double f, const double t, const double complex Yfactorplus, const double complex Yfactorcross);
double complex G31mode(const double f, const double t, const double complex Yfactorplus, const double complex Yfactorcross);

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LISAGEOMETRY_H */
