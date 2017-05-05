/**
 * \author Sylvain Marsat, University of Maryland - NASA GSFC
 *
 * \brief C code headers for the implementation of the Fourier domain response for LISA-like detectors
 *
 */

#ifndef _LISAFDRESPONSE_H
#define _LISAFDRESPONSE_H

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
#include "waveform.h"
#include "EOBNRv2HMROMstruct.h"
#include "LISAgeometry.h"


#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/**************************************************/
/**************** Prototypes **********************/

/* Core function processing a signal (in the form of a list of modes) through the Fourier-domain LISA response, for given values of the inclination, position in the sky and polarization angle - here simplified version for just the y_21 observable */
int LISASimFDResponse21(
  LISAconstellation *variant,                                 /* Provides specifics on the variant of LISA */
  struct tagListmodesCAmpPhaseFrequencySeries **list,  /* Input/Output: list of modes in Frequency-domain amplitude and phase form as produced by the ROM, and output after FD response processing */
  const double inclination,                                   /* Inclination of the source */
  const double lambda,                                        /* First angle for the position in the sky */
  const double beta,                                          /* Second angle for the position in the sky */
  const double psi);                                          /* Polarization angle */

//WARNING: tRef is ignored for now in the response - i.e. set to 0
/* Core function processing a signal (in the form of a list of modes) through the Fourier-domain LISA response, for given values of the inclination, position in the sky and polarization angle */
int LISASimFDResponsey12(
  LISAconstellation *variant,                                 /* Provides specifics on the variant of LISA */
  struct tagListmodesCAmpPhaseFrequencySeries **list,      /* Input: list of modes in Frequency-domain amplitude and phase form as produced by the ROM */
  struct tagListmodesCAmpPhaseFrequencySeries **listy12,   /* Output: list of contribution of each mode in Frequency-domain amplitude and phase form, in the y12 observable */
  const double tRef,                                          /* Time at coalescence */
  const double lambda,                                        /* First angle for the position in the sky */
  const double beta,                                          /* Second angle for the position in the sky */
  const double inclination,                                   /* Inclination of the source */
  const double psi);                                          /* Polarization angle */

//WARNING: tRef is ignored for now in the response - i.e. set to 0
/* Core function processing a signal (in the form of a list of modes) through the Fourier-domain LISA response, for given values of the inclination, position in the sky and polarization angle */
int LISASimFDResponseTDI3Chan(
  int tagtRefatLISA,                                          /* 0 to measure Tref from SSB arrival, 1 at LISA guiding center */
  LISAconstellation *variant,                                 /* Provides specifics on the variant of LISA */
  struct tagListmodesCAmpPhaseFrequencySeries **list,      /* Input: list of modes in Frequency-domain amplitude and phase form as produced by the ROM */
  struct tagListmodesCAmpPhaseFrequencySeries **listTDI1,  /* Output: list of contribution of each mode in Frequency-domain amplitude and phase form, in the TDI channel 1 */
  struct tagListmodesCAmpPhaseFrequencySeries **listTDI2,  /* Output: list of contribution of each mode in Frequency-domain amplitude and phase form, in the TDI channel 2 */
  struct tagListmodesCAmpPhaseFrequencySeries **listTDI3,  /* Output: list of contribution of each mode in Frequency-domain amplitude and phase form, in the TDI channel 3 */
  const double tRef,                                          /* Time at coalescence */
  const double lambda,                                        /* First angle for the position in the sky */
  const double beta,                                          /* Second angle for the position in the sky */
  const double inclination,                                   /* Inclination of the source */
  const double psi,                                           /* Polarization angle */
  const double maxf,                                          /* Maximal frequency to consider - used to ignore hard-to-resolve response at f>1Hz - NOTE: for now, no recomputation of the boundary, so when not resampling can lose a bit of support between the last frequency point covered and maxf */
  const TDItag tditag);                                       /* Selector for the set of TDI observables */

// int LISASimFDResponseTDI1Chan(
//   struct tagListmodesCAmpPhaseFrequencySeries **list,      /* Input: list of modes in Frequency-domain amplitude and phase form as produced by the ROM */
//   struct tagListmodesCAmpPhaseFrequencySeries **listTDI,   /* Output: list of contribution of each mode in Frequency-domain amplitude and phase form, in the TDI channel 1 */
//   const double tRef,                                          /* Time at coalescence */
//   const double lambda,                                        /* First angle for the position in the sky */
//   const double beta,                                          /* Second angle for the position in the sky */
//   const double inclination,                                   /* Inclination of the source */
//   const double psi,                                           /* Polarization angle */
//   const double maxf,                                          /* Maximal frequency to consider - used to ignore hard-to-resolve response at f>1Hz - NOTE: for now, no recomputation of the boundary, so when not resampling can lose a bit of support between the last frequency point covered and maxf */
//   const TDItag tditag);                                       /* Selector for the set of TDI observables */

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LISAFDRESPONSE_H */
