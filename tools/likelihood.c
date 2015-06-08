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
  double (*Snoise)(double) ) /* Noise function */
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
  /* Here we add the complex conjugate apparently missing, copying the data to do it in a non-destructive way -  would be better to implement directly in wip_phase */  
  //JGB: This shouldn't be necessary.  Before 6/8/15 there was a bug making phase sign inconsistent.)
  gsl_vector* amp_imag2 = gsl_vector_alloc(n2);
  gsl_vector* phase2 = gsl_vector_alloc(n2);
  gsl_vector_memcpy(amp_imag2, h2->amp_imag);
  gsl_vector_memcpy(phase2, h2->phase);
  //JGB:So I'm commenting these out...
  //gsl_vector_scale(amp_imag2, -1);
  //gsl_vector_scale(phase2, -1);
  double *h2Ai = amp_imag2->data;
  double *h2p = phase2->data;
  //JGB: Set these somehwere else where it makes sense.  f_min or f_max <= 0 means use intersection of signal domains.
  double f_min=10.0;
  double f_max=-1.0;
  double overlap = wip_phase(f1, n1, f2, n2, h1Ar, h1Ai, h1p, h2Ar, h2Ai, h2p, Snoise, 1.0,f_min,f_max);
  return overlap;
}

/* Function computing the overlap (h1|h2) between two waveforms given as list of modes, for a given noise function - all the cross-products between modes are taken into account */
double FDListmodesOverlap(
  struct tagListmodesCAmpPhaseFrequencySeries *listh1, /* First mode h1, list of modes in amplitude/phase form */
  struct tagListmodesCAmpPhaseFrequencySeries *listh2, /* Second mode h2, list of modes in amplitude/phase form */
  double (*Snoise)(double)) /* Noise function */
{
  double overlap = 0;

  /* Main loop over the modes - goes through all the modes present */
  ListmodesCAmpPhaseFrequencySeries* listelementh1 = listh1;
  while(listelementh1) {
    ListmodesCAmpPhaseFrequencySeries* listelementh2 = listh2;
    while(listelementh2) {
      overlap += FDSinglemodeOverlap(listelementh1->freqseries, listelementh2->freqseries, Snoise);
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
  double (*Snoise)(double)) /* Noise function */
{
  double overlap = 0;

  /* Main loop over the modes - goes through all the modes present */
  ListmodesCAmpPhaseFrequencySeries* listelementh1 = listh1;
  while(listelementh1) {
    ListmodesCAmpPhaseFrequencySeries* listelementh2 = ListmodesCAmpPhaseFrequencySeries_GetMode(listh2, listelementh1->l, listelementh1->m);
    overlap += FDSinglemodeOverlap(listelementh1->freqseries, listelementh2->freqseries, Snoise);
    listelementh1 = listelementh1->next;
  }
  return overlap;
}

/* Function computing the log likelihood (h|s) - 1/2 (h|h) - 1/2 (s|s), with s the signal, h the template, and where we keep the constant term (s|s) - passed to the function as a parameter - all the cross-products between modes are taken into account */
double FDLogLikelihood(
  struct tagListmodesCAmpPhaseFrequencySeries *lists,  /* Input: list of modes for the signal s, in Frequency-domain amplitude and phase form */
  struct tagListmodesCAmpPhaseFrequencySeries *listh,  /* Input: list of modes for the template, in Frequency-domain amplitude and phase form */
  double (*Snoise)(double),                            /* Noise function */
  double ss,                                           /* Inner product (s|s), constant to be computed elsewhere and passed as an argument */
  double hh )                                         /* Inner product (h|h), constant to be computed elsewhere and passed as an argument */
{
  return FDListmodesOverlap(lists, listh, Snoise) - 1./2 * hh - 1/2 * ss;
}
