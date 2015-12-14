/**
 * \author Sylvain Marsat, University of Maryland - NASA GSFC
 *
 * \brief C header for the instrumental noise for LISA-type detectors.
 *
 * Formulas taken from Królak&al gr-qc/0401108 (c.f. section III).
 *
 */

#ifndef _LISANOISE_H
#define _LISANOISE_H

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
#include "LISAgeometry.h"


/************************************************************************/
/****** Global variables storing min and max f for the noise PSD  *******/

#define __LISASimFD_Noise_fLow 1.e-5
#define __LISASimFD_Noise_fHigh 1. 

/**************************************************************/
/****** Prototypes: functions evaluating the noise PSD  *******/

/* Function returning the relevant noise function, given a set of TDI observables and a channel */
RealFunction* NoiseFunction(const TDItag tditag, const int nchan);

/* Noise Sn for TDI observables - factors have been scaled out both in the response and the noise */
double SnXYZ(double f);
double Snalphabetagamma(double f);
double SnAXYZ(double f);
double SnEXYZ(double f);
double SnTXYZ(double f);
double SnAalphabetagamma(double f);
double SnEalphabetagamma(double f);
double SnTalphabetagamma(double f);

/* Function returning the relevant noise function, given a set of TDI observables and a channel */
double (*NoiseFunction(const TDItag tditag, const int chan))(double);

/* The noise functions themselves */
/* double NoiseSnA(const double f); */
/* double NoiseSnE(const double f); */
/* double NoiseSnT(const double f); */

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LISANOISE_H */
