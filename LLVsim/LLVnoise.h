/**
 * \author Sylvain Marsat, University of Maryland - NASA GSFC
 *
 * \brief C header for the initialization of the instrumental noise for LIGO/VIRGO detectors.
 *
 *
 */

#ifndef _LLVNOISE_H
#define _LLVNOISE_H

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


/**************************************************************************/
/****** Prototypes: functions loading and evaluating the noise PSD  *******/

/* Function parsing the environment variable $LLV_NOISE_DATA_PATH and trying to run LLVSimFD_Noise_Init in each */
void LLVSimFD_Noise_Init_ParsePath(void);
/* Function loading the noise data from a directory */
int LLVSimFD_Noise_Init(const char dir[]);

/* The noise functions themselves */
double NoiseSnLHO(const double f);
double NoiseSnLLO(const double f);
double NoiseSnVIRGO(const double f);

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LLVNOISE_H */
