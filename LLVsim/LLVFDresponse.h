/**
 * \author Sylvain Marsat, University of Maryland - NASA GSFC
 *
 * \brief C code headers for the implementation of the Fourier domain response for LIGO-VIRGO detectors
 *
 *
 */

#ifndef _LLVFDRESPONSE_H
#define _LLVFDRESPONSE_H

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
#include "LLVgeometry.h"
#include "struct.h"
#include "waveform.h"
#include "timeconversion.h"


/**************************************************/
/**************** Prototypes **********************/

/* Function to convert string input network string to Networktag */
Networktag ParseNetworktag(char* string);

/* Core function processing a signal (in the form of a list of modes) through the Fourier-domain LLV response (for a given detector), for given values of the inclination, position in the sky and polarization angle */
int LLVSimFDResponse(
  struct tagListmodesCAmpPhaseFrequencySeries **listhlm,  /* Input: list of modes in Frequency-domain amplitude and phase form as produced by the ROM */
  struct tagListmodesCAmpPhaseFrequencySeries **lists,    /* Output: list of contribution of each mode in the detector signal, in Frequency-domain amplitude and phase form, for the given detector and sky position */
  const double gpstime,                                   /* GPS time (s) when the signal at coalescence reaches geocenter */
  const double ra,                                        /* Position in the sky: J2000.0 right ascension (rad) */
  const double dec,                                       /* Position in the sky: J2000.0 declination (rad) */
  const double inclination,                               /* Inclination of the source (rad) */
  const double psi,                                       /* Polarization angle (rad) */
  const Detectortag tag);                                 /* Tag identifying the detector */

  /* Core function processing a signal (in the form of a list of modes) through the Fourier-domain LLV response, for given values of the inclination, position in the sky and polarization angle */
  int LLVSimFDResponse3Det(
    struct tagListmodesCAmpPhaseFrequencySeries **list,      /* Input: list of modes in Frequency-domain amplitude and phase form as produced by the ROM */
    struct tagListmodesCAmpPhaseFrequencySeries **listDet1,  /* Output: list of contribution of each mode in Frequency-domain amplitude and phase form, in the detector 1 */
    struct tagListmodesCAmpPhaseFrequencySeries **listDet2,  /* Output: list of contribution of each mode in Frequency-domain amplitude and phase form, in the detector 2 */
    struct tagListmodesCAmpPhaseFrequencySeries **listDet3,  /* Output: list of contribution of each mode in Frequency-domain amplitude and phase form, in the detector 3 */
    const double gpstime,                                    /* GPS time (s) when the signal at coalescence reaches geocenter */
    const double ra,                                            /* First angle for the position in the sky */
    const double dec,                                           /* Second angle for the position in the sky */
    const double inclination,                                   /* Inclination of the source */
    const double psi,                                           /* Polarization angle */
    const Networktag tag);                               /* Selector for the detector network */

/* Function setting the response matrix of a given detector, in cartesian coordinates */
void SetMatrixD(
  gsl_matrix* D,                       /* Output: matrix of the detector response Dij */
  const Detectortag tag);              /* Tag identifying the detector */
/* Function setting the position of a detector, in cartesian coordinates */
void SetVectorXd(
  gsl_vector* Xd,                      /* Output: position vector of the detector */
  const Detectortag tag);              /* Tag identifying the detector */

/* Function setting the cartesian coordinates of the wave frame vectors (X,Y,Z), given the position in the sky and polarization */
void SetVectorsXYZ(
  gsl_vector* X,                       /* Output: cartesian vector of the wave frame unit vector X */
  gsl_vector* Y,                       /* Output: cartesian vector of the wave frame unit vector Y */
  gsl_vector* Z,                       /* Output: cartesian vector of the wave frame unit vector Z */
  const double theta,                  /* First angle for the position in the sky (Earth-based spherical angle) */
  const double phi,                    /* Second angle for the position in the sky (Earth-based spherical angle) */
  const double psi);                   /* Polarization angle */

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LLVFDRESPONSE_H */
