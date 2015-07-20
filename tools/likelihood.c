/**
 * \author Sylvain Marsat, University of Maryland - NASA GSFC
 *
 * \brief C code for the implementation of the Fourier-domain overlaps, likelihoods.
 *
 */


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
#include "likelihood.h"

#include "wip.h"

#include <time.h> /* for testing */


/***********************************************************/
/********* Core functions to compute overlaps **************/

/* Function computing the overlap (h1|h2) between two given modes, for a given noise function */
double FDSinglemodeOverlap(
  struct tagCAmpPhaseFrequencySeries *freqseries1, /* First mode h1, in amplitude/phase form */
  struct tagCAmpPhaseFrequencySeries *freqseries2, /* Second mode h2, in amplitude/phase form */
  double (*Snoise)(double),                        /* Noise function */
  double fLow,                                     /* Lower bound of the frequency window for the detector */
  double fHigh)                                    /* Upper bound of the frequency window for the detector */
{
  /* Should add some error checking */
  CAmpPhaseFrequencySeries* h1 = freqseries1;
  CAmpPhaseFrequencySeries* h2 = freqseries2;

  double *f1 = h1->freq->data;
  int n1 = h1->freq->size;
  double *h1Ar = h1->amp_real->data;
  double *h1Ai = h1->amp_imag->data;
  double *h1p = h1->phase->data;

  double *f2 = h2->freq->data;
  int n2 = h2->freq->size;
  double *h2Ar = h2->amp_real->data;
  double *h2Ai = h2->amp_imag->data;
  double *h2p = h2->phase->data;

  /* fLow or fHigh <= 0 means use intersection of signal domains  */
  //NOTE: factor 4 was apparently missing
  double overlap = 4.*wip_phase(f1, n1, f2, n2, h1Ar, h1Ai, h1p, h2Ar, h2Ai, h2p, Snoise, 1.0, fLow, fHigh);
  return overlap;
}

/* Function computing the overlap (h1|h2) between two waveforms given as list of modes, for a given noise function - all the cross-products between modes are taken into account */
double FDListmodesOverlap(
  struct tagListmodesCAmpPhaseFrequencySeries *listh1, /* First mode h1, list of modes in amplitude/phase form */
  struct tagListmodesCAmpPhaseFrequencySeries *listh2, /* Second mode h2, list of modes in amplitude/phase form */
  double (*Snoise)(double),                            /* Noise function */
  double fLow,                                         /* Lower bound of the frequency window for the detector */
  double fHigh)                                        /* Upper bound of the frequency window for the detector */
{
  double overlap = 0;

  /* Main loop over the modes - goes through all the modes present */
  //printf("FDListmodesOverlap:\n");
  ListmodesCAmpPhaseFrequencySeries* listelementh1 = listh1;
  while(listelementh1) {
    ListmodesCAmpPhaseFrequencySeries* listelementh2 = listh2;
    while(listelementh2) {
      //printf("(%d%d,%d%d): %g\n", listelementh1->l, listelementh1->m, listelementh2->l, listelementh2->m, FDSinglemodeOverlap(listelementh1->freqseries, listelementh2->freqseries, Snoise, fLow, fHigh));
      overlap += FDSinglemodeOverlap(listelementh1->freqseries, listelementh2->freqseries, Snoise, fLow, fHigh);
      listelementh2 = listelementh2->next;
    }
    listelementh1 = listelementh1->next;
  }
  return overlap;
}

/* Function computing the overlap (h1|h2) between two waveforms given as list of modes, for a given noise function - no cross-products between modes are taken into account */
double FDListmodesOverlapNoCrossTerms(
  struct tagListmodesCAmpPhaseFrequencySeries *listh1, /* First mode h1, list of modes in amplitude/phase form */
  struct tagListmodesCAmpPhaseFrequencySeries *listh2, /* Second mode h2, list of modes in amplitude/phase form */
  double (*Snoise)(double),                            /* Noise function */
  double fLow,                                         /* Lower bound of the frequency window for the detector */
  double fHigh)                                        /* Upper bound of the frequency window for the detector */
{
  double overlap = 0;

  /* Main loop over the modes - goes through all the modes present */
  ListmodesCAmpPhaseFrequencySeries* listelementh1 = listh1;
  while(listelementh1) {
    ListmodesCAmpPhaseFrequencySeries* listelementh2 = ListmodesCAmpPhaseFrequencySeries_GetMode(listh2, listelementh1->l, listelementh1->m);
    overlap += FDSinglemodeOverlap(listelementh1->freqseries, listelementh2->freqseries, Snoise, fLow, fHigh);
    listelementh1 = listelementh1->next;
  }
  return overlap;
}

/* Function computing the log likelihood (h|s) - 1/2 (h|h) - 1/2 (s|s), with s the signal, h the template, and where we keep the constant term (s|s) - passed to the function as a parameter - all the cross-products between modes are taken into account */
double FDLogLikelihood(
  struct tagListmodesCAmpPhaseFrequencySeries *lists,  /* Input: list of modes for the signal s, in Frequency-domain amplitude and phase form */
  struct tagListmodesCAmpPhaseFrequencySeries *listh,  /* Input: list of modes for the template, in Frequency-domain amplitude and phase form */
  double (*Snoise)(double),                            /* Noise function */
  double fLow,                                         /* Lower bound of the frequency window for the detector */
  double fHigh,                                        /* Upper bound of the frequency window for the detector */
  double ss,                                           /* Inner product (s|s), constant to be computed elsewhere and passed as an argument */
  double hh)                                          /* Inner product (h|h), constant to be computed elsewhere and passed as an argument */
{
  return FDListmodesOverlap(lists, listh, Snoise, fLow, fHigh) - 1./2 * hh - 1./2 * ss;
}

/* Function computing the overlap (h1|h2) between two waveforms given as list of modes, for a given noise function - no cross-products between modes are taken into account - adds an additional parameter for the starting 22-mode frequency (then properly scaled for the other modes) for a limited duration of the observations */
double LISAFDListmodesOverlap(
  struct tagListmodesCAmpPhaseFrequencySeries *listh1, /* First mode h1, list of modes in amplitude/phase form */
  struct tagListmodesCAmpPhaseFrequencySeries *listh2, /* Second mode h2, list of modes in amplitude/phase form */
  double (*Snoise)(double),                            /* Noise function */
  double fLow,                                         /* Lower bound of the frequency window for the detector */
  double fHigh,                                        /* Upper bound of the frequency window for the detector */
  double fstartobs)                                    /* Starting frequency for the 22 modes - as determined from a limited duration of the observation */
{
  double overlap = 0;

  /* Main loop over the modes - goes through all the modes present */
  //printf("FDListmodesOverlap:\n");
  ListmodesCAmpPhaseFrequencySeries* listelementh1 = listh1;
  while(listelementh1) {
    ListmodesCAmpPhaseFrequencySeries* listelementh2 = listh2;
    while(listelementh2) {
      //printf("(%d%d,%d%d): %g\n", listelementh1->l, listelementh1->m, listelementh2->l, listelementh2->m, FDSinglemodeOverlap(listelementh1->freqseries, listelementh2->freqseries, Snoise, fLow, fHigh));
      double fcutLow = fmax(fLow, fmax(((double) listelementh1->m)/2. * fstartobs, ((double) listelementh2->m)/2. * fstartobs) );
      overlap += FDSinglemodeOverlap(listelementh1->freqseries, listelementh2->freqseries, Snoise, fcutLow, fHigh);
      listelementh2 = listelementh2->next;
    }
    listelementh1 = listelementh1->next;
  }
  return overlap;
}

/* Function computing the log likelihood (h|s) - 1/2 (h|h) - 1/2 (s|s), with s the signal, h the template, and where we keep the constant term (s|s) - passed to the function as a parameter - all the cross-products between modes are taken into account - adds an additional parameter for the starting 22-mode frequency (then properly scaled for the other modes) for a limited duration of the observations */
double LISAFDLogLikelihood(
  struct tagListmodesCAmpPhaseFrequencySeries *lists,  /* Input: list of modes for the signal s, in Frequency-domain amplitude and phase form */
  struct tagListmodesCAmpPhaseFrequencySeries *listh,  /* Input: list of modes for the template, in Frequency-domain amplitude and phase form */
  double (*Snoise)(double),                            /* Noise function */
  double fLow,                                         /* Lower bound of the frequency window for the detector */
  double fHigh,                                        /* Upper bound of the frequency window for the detector */
  double ss,                                           /* Inner product (s|s), constant to be computed elsewhere and passed as an argument */
  double hh,                                           /* Inner product (h|h), constant to be computed elsewhere and passed as an argument */
  double fstartobs)                                    /* Starting frequency for the 22 modes - as determined from a limited duration of the observation */
{
  return LISAFDListmodesOverlap(lists, listh, Snoise, fLow, fHigh, fstartobs) - 1./2 * hh - 1./2 * ss;
}

/* Newtonian estimate of the relation Mf(deltat/M) (for the 22 mode) - gives the starting geometric frequency for a given mass ratio and a given geometric duration of the observations */
double NewtonianfoftGeom(
  const double q,                      /* Mass ratio m1/m2 */
  const double t)                      /* Duration of the observations in geometric units (t/M) */
{
  double nu = q/(1.+q)/(1.+q);
  return 1./PI * pow(256*nu/5. * t, -3./8);
}
