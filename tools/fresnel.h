/**
 * \author Sylvain Marsat, University of Maryland - NASA GSFC
 *
 * \brief C header for functions computing coeffcients of Not-A-Knot and Quadratic splines in matrix form.
 *
 */

#ifndef _FRESNEL_H
#define _FRESNEL_H

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
#include "struct.h"


#if defined(__cplusplus)
#define complex _Complex
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

double complex ComputeInt(
  gsl_matrix* splinecoeffsAreal,         /*  */
  gsl_matrix* splinecoeffsAimag,         /*  */
  gsl_matrix* splinecoeffsphase);        /*  */


#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _FRESNEL_H */
