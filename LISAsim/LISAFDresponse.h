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
#include "EOBNRv2HMROMstruct.h"
#include "LISAgeometry.h"


/**************************************************/
/**************** Prototypes **********************/

/* Core function processing a signal (in the form of a list of modes) through the Fourier-domain LISA response, for given values of the inclination, position in the sky and polarization angle - here simplified version for just the y_21 observable */
int LISASimFDResponse21(
  struct tagListmodesCAmpPhaseFrequencySeries **list,  /* Input/Output: list of modes in Frequency-domain amplitude and phase form as produced by the ROM, and output after FD response processing */
  const double inclination,                                   /* Inclination of the source */
  const double lambda,                                        /* First angle for the position in the sky */
  const double beta,                                          /* Second angle for the position in the sky */
  const double psi);                                          /* Polarization angle */

//WARNING: tRef is ignored for now in the response - i.e. set to 0
/* Core function processing a signal (in the form of a list of modes) through the Fourier-domain LISA response, for given values of the inclination, position in the sky and polarization angle */
int LISASimFDResponseTDIAET(
  struct tagListmodesCAmpPhaseFrequencySeries **list,  /* Input: list of modes in Frequency-domain amplitude and phase form as produced by the ROM */
  struct tagListmodesCAmpPhaseFrequencySeries **listA,  /* Output: list of contribution of each mode in Frequency-domain amplitude and phase form, in the second-generation TDI observable A */
  struct tagListmodesCAmpPhaseFrequencySeries **listE,  /* Output: list of contribution of each mode in Frequency-domain amplitude and phase form, in the second-generation TDI observable E */
  struct tagListmodesCAmpPhaseFrequencySeries **listT,  /* Output: list of contribution of each mode in Frequency-domain amplitude and phase form, in the second-generation TDI observable T */
  const double tRef,                                          /* Time at coalescence */
  const double lambda,                                        /* First angle for the position in the sky */
  const double beta,                                          /* Second angle for the position in the sky */
  const double inclination,                                   /* Inclination of the source */
  const double psi);                                          /* Polarization angle */ 

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LISAFDRESPONSE_H */
