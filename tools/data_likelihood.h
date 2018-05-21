/**
 * \author John Baker NASA - GSFC
 *
 * \brief C header  for the implementing likelihoods with Data provided on a fixed frequency domain grid.
 *
 *
 */

#ifndef _DATA_LIKELIHOOD_H
#define _DATA_LIKELIHOOD_H

#define _XOPEN_SOURCE 500

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

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
#include "splinecoeffs.h"

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif


void WhitenData(
  ReImUniformFrequencySeries** whitened_data, /* Output */
  ReImUniformFrequencySeries* data,           /* Output */
  ObjectFunction * Snoise,                    /* Noise function */
  double fLow,                                /* if fLow/fHigh > 0 restrict whitented data to withing that range */  
  double fHigh);

double ReIntProdCAmpPhaseSpline(
  CAmpPhaseSpline* splines,                    //input "model"
  ReImUniformFrequencySeries* Dfreqseries,  //input defines ReIm frequency series to be multiplied with spline evals and summed
  int imin,  // min/max  index-range of the freqseries to use include in the sum 
  int imax);

double FDSinglemodeDataOverlap(
  CAmpPhaseFrequencySeries* model,     
  ReImUniformFrequencySeries* datachan);

double FDDataOverlap(
  ListmodesCAmpPhaseFrequencySeries *listmodelchan, /* Waveform channel channel, list of modes in amplitude/phase form */
  ReImUniformFrequencySeries *datachan);                  /* Data channel channel */


#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _DATA_LIKELIHOOD_H */



