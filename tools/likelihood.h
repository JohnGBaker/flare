/**
 * \author Sylvain Marsat, University of Maryland - NASA GSFC
 *
 * \brief C header for the computation of the Fourier-domain overlaps, likelihoods.
 *
 *
 */

#ifndef _LIKELIHOOD_H
#define _LIKELIHOOD_H

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
#include "wip.h"


/****** Prototypes: utilities *******/

/* Newtonian estimate of the relation Mf(deltat/M) (for the 22 mode) - gives the starting geometric frequency for a given mass ratio and a given geometric duration of the observations */
double NewtonianfoftGeom(
  const double q,                      /* Mass ratio m1/m2 */
  const double t);                     /* Duration of the observations in geometric units (t/M) */

/* Function to evaluate a Noise function  */
void EvaluateNoise(
  gsl_vector* noisevalues,                         /* Output: vector of the noise values */
  gsl_vector* freq,                                /* Input: vector of frequencies on which to evaluate */
  double (*Snoise)(double),                        /* Noise function */
  double fLow,                                     /* Lower bound of the frequency window for the detector */
  double fHigh);                                   /* Upper bound of the frequency window for the detector */

/* Function building a frequency vector with logarithmic sampling */
void SetLogFrequencies(
  gsl_vector* freqvector,    /* Output pointer to gsl_vector, already allocated */
  const double fmin,         /* Lower bound of the frequency interval */
  const double fmax,         /* Upper bound of the frequency interval */
  const int nbpts);          /* Number of points */

/* Function determining logarithmically spaced frequencies from a waveform given as a list of modes, a minimal frequency and a number of points */
void ListmodesSetLogFrequencies(
  struct tagListmodesCAmpPhaseFrequencySeries *list,     /* Waveform, list of modes in amplitude/phase form */
  double fLow,                                           /* Additional lower frequency limit */
  double fHigh,                                          /* Additional upper frequency limit */
  int nbpts,                                             /* Number of frequency samples */
  gsl_vector* freqvector);                               /* Output: vector of frequencies, already allocated with nbpts */

/****** Prototypes: overlaps and likelihood computation *******/

/***************************** Functions for overlaps using linear integration ******************************/

/* Function computing the overlap (h1|h2) between two given modes, for a given noise function - uses simple trapeze integration on logarithmically sampled frequencies  */
double FDSinglemodeLogLinearOverlap(
  struct tagCAmpPhaseFrequencySeries *freqseries1, /* First mode h1, in amplitude/phase form */
  struct tagCAmpPhaseFrequencySeries *freqseries2, /* Second mode h2, in amplitude/phase form */
  double (*Snoise)(double),                        /* Noise function */
  double fLow,                                     /* Lower bound of the frequency window for the detector */
  double fHigh);                                   /* Upper bound of the frequency window for the detector */

/* Function computing the overlap (h1|h2) between two waveforms given as lists of mode contributions (factos sYlm already included), for a given noise function - uses simple trapeze integration on logarithmically sampled frequencies - generates the frequency series in Re/Im form by summing the mode contributions first, then computes the overlap */
double FDListmodesLogLinearOverlap(
  struct tagListmodesCAmpPhaseFrequencySeries *list1,    /* First waveform, list of modes in amplitude/phase form */
  struct tagListmodesCAmpPhaseFrequencySeries *list2,    /* First waveform, list of modes in amplitude/phase form */
  double (*Snoise)(double),                              /* Noise function */
  double fLow,                                           /* Lower bound of the frequency window for the detector */
  double fHigh,                                          /* Upper bound of the frequency window for the detector */
  double fstartobs1,                                     /* Starting frequency for the 22 mode of wf 1 - as determined from a limited duration of the observation - set to 0 to ignore */
  double fstartobs2);                                    /* Starting frequency for the 22 mode of wf 2 - as determined from a limited duration of the observation - set to 0 to ignore */

/* Function computing the overlap (h1|h2) between h1 given in Re/Im form and h2 given as a list of mode contributions (factos sYlm already included), for a given vector of noise values - uses simple trapeze integration - generates the frequency series of h2 in Re/Im form by summing the mode contributions first, then computes the overlap */
double FDOverlapReImvsListmodesCAmpPhase(
  struct tagReImFrequencySeries *freqseries1,       /* First waveform, in Re/Im form */
  struct tagListmodesCAmpPhaseFrequencySeries *list2,    /* Second waveform, list of modes in amplitude/phase form */
  gsl_vector* noisevalues,                               /* Vector for the noise values on the freq of h1 */
  double fstartobs2);                                    /* Starting frequency for the 22 mode of wf 2 - as determined from a limited duration of the observation - set to 0 to ignore */


/* Function computing the overlap (h1|h2) between two waveforms given as Re/Im frequency series (common freq values), for a given vector of noise values - uses simple trapeze integration */
double FDOverlapReImvsReIm(
  struct tagReImFrequencySeries *h1,  /* First waveform, frequency series in Re/Im form */
  struct tagReImFrequencySeries *h2,  /* Second waveform, frequency series in Re/Im form */
  gsl_vector* noisevalues);           /* Vector for the noise values on common freq of the freqseries */

/***************************** Functions for overlaps using amplitude/phase (wip) ******************************/

/* Function computing the overlap (h1|h2) between two given modes in amplitude/phase form, for a given noise function - uses the amplitude/phase representation (wip) */
double FDSinglemodeWIPOverlap(
  struct tagCAmpPhaseFrequencySeries *freqseries1, /* First mode h1, in amplitude/phase form */
  struct tagCAmpPhaseFrequencySeries *freqseries2, /* Second mode h2, in amplitude/phase form */
  double (*Snoise)(double),                        /* Noise function */
  double fLow,                                     /* Lower bound of the frequency window for the detector */
  double fHigh);                                   /* Upper bound of the frequency window for the detector */

/* Function computing the overlap (h1|h2) between two waveforms given as list of modes, for a given noise function - no cross-products between modes are taken into account - two additional parameters for the starting 22-mode frequencies (then properly scaled for the other modes) for a limited duration of the observations */
double FDListmodesWIPOverlap(
  struct tagListmodesCAmpPhaseFrequencySeries *listh1, /* First waveform, list of modes in amplitude/phase form */
  struct tagListmodesCAmpPhaseFrequencySeries *listh2, /* Second waveform, list of modes in amplitude/phase form */
  double (*Snoise)(double),                            /* Noise function */
  double fLow,                                         /* Lower bound of the frequency window for the detector */
  double fHigh,                                        /* Upper bound of the frequency window for the detector */
  double fstartobs1,                                   /* Starting frequency for the 22 mode of wf 1 - as determined from a limited duration of the observation - set to 0 to ignore */
  double fstartobs2);                                  /* Starting frequency for the 22 mode of wf 2 - as determined from a limited duration of the observation - set to 0 to ignore */

/************** Functions for overlap/likelihood allowing to switch between wip and loglinear integration *****************/

/* Wrapping of FDListmodesWIPOverlap or FDListmodesLogLinearOverlap according to tagint */
double FDListmodesOverlap(
  struct tagListmodesCAmpPhaseFrequencySeries *listh1, /* First mode h1, list of modes in amplitude/phase form */
  struct tagListmodesCAmpPhaseFrequencySeries *listh2, /* Second mode h2, list of modes in amplitude/phase form */
  double (*Snoise)(double),                            /* Noise function */
  double fLow,                                         /* Lower bound of the frequency window for the detector */
  double fHigh,                                        /* Upper bound of the frequency window for the detector */
  double fstartobs1,                                   /* Starting frequency for the 22 mode of wf 1 - as determined from a limited duration of the observation - set to 0 to ignore */
  double fstartobs2,                                   /* Starting frequency for the 22 mode of wf 2 - as determined from a limited duration of the observation - set to 0 to ignore */
  int tagint);                                         /* Tag choosing the integrator: 0 for wip, 1 for log linear integration */

/* Function computing the log likelihood (h|s) - 1/2 (h|h) - 1/2 (s|s), with s the signal, h the template, and where we keep the constant term (s|s) - passed to the function as a parameter - all the cross-products between modes are taken into account - two additionals parameters for the starting 22-mode frequencies (then properly scaled for the other modes) for a limited duration of the observations */
double FDLogLikelihood(
  struct tagListmodesCAmpPhaseFrequencySeries *lists,  /* Input: list of modes for the signal s, in Frequency-domain amplitude and phase form */
  struct tagListmodesCAmpPhaseFrequencySeries *listh,  /* Input: list of modes for the template, in Frequency-domain amplitude and phase form */
  double (*Snoise)(double),                            /* Noise function */
  double fLow,                                         /* Lower bound of the frequency window for the detector */
  double fHigh,                                        /* Upper bound of the frequency window for the detector */
  double ss,                                           /* Inner product (s|s), constant to be computed elsewhere and passed as an argument */
  double hh,                                           /* Inner product (h|h), constant to be computed elsewhere and passed as an argument */
  double fstartobss,                                   /* Starting frequency for the 22 mode of s - as determined from a limited duration of the observation - set to 0 to ignore */
  double fstartobsh,                                   /* Starting frequency for the 22 mode of h - as determined from a limited duration of the observation - set to 0 to ignore */
  int tagint);                                         /* Tag choosing the integrator: 0 for wip, 1 for log linear integration */

/************** Functions for overlap/likelihood specific to loglinear integration *****************/

/* Function computing the log likelihood -1/2(h-s|h-s), with s the signal, h the template (both given as frequency series in Re/Im form) - h, s assumed to be given on the same set of frequencies - same for the vector of noise values - fstartobs for s, h has already been taken into account */
double FDLogLikelihoodReIm(
  struct tagReImFrequencySeries *s,    /* First waveform (injection), frequency series in Re/Im form */
  struct tagReImFrequencySeries *h,    /* Second waveform (template), frequency series in Re/Im form */
  gsl_vector* noisevalues);            /* Vector for the noise values on common freq of the freqseries */

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LIKELIHOOD_H */
