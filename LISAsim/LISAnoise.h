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


#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/************************************************************************/
/****** Global variables storing min and max f for the noise PSD  *******/

/* Defines bounds in frequency beyond which we don't trust the instrument model anymore - all waveforms will be cut to this range */
/* Here extended range - allows for instance to taper the FD signal for f>1Hz */
#define __LISASimFD_Noise_fLow 1.e-6
#define __LISASimFD_Noise_fHigh 5.
/* Original, more conservative bounds */
//#define __LISASimFD_Noise_fLow 1.e-5
//#define __LISASimFD_Noise_fHigh 1.

/**************************************************************/
/****** Prototypes: functions evaluating the noise PSD  *******/

/* Function returning the relevant noise function, given a set of TDI observables and a channel */
RealFunctionPtr NoiseFunction(const TDItag tditag, const int nchan);

/* Noise Sn for TDI observables - factors have been scaled out both in the response and the noise */
double SnXYZ(const LISAconstellation *variant, double f);
double Snalphabetagamma(const LISAconstellation *variant, double f);
double SnAXYZ(const LISAconstellation *variant, double f);
double SnEXYZ(const LISAconstellation *variant, double f);
double SnTXYZ(const LISAconstellation *variant, double f);
double SnAalphabetagamma(const LISAconstellation *variant, double f);
double SnEalphabetagamma(const LISAconstellation *variant, double f);
double SnTalphabetagamma(const LISAconstellation *variant, double f);

/* Noise functions for AET(XYZ) without rescaling */
double SnAXYZNoRescaling(const LISAconstellation *variant, double f);
double SnEXYZNoRescaling(const LISAconstellation *variant, double f);
double SnTXYZNoRescaling(const LISAconstellation *variant, double f);

/* Function returning the relevant noise function, given a set of TDI observables and a channel */
/* double (*NoiseFunction(const TDItag tditag, const int chan))(double); */

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
