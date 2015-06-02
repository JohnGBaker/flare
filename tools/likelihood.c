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
  gsl_vector* amp_imag2 = gsl_vector_alloc(n2);
  gsl_vector* phase2 = gsl_vector_alloc(n2);
  gsl_vector_memcpy(amp_imag2, h2->amp_imag);
  gsl_vector_memcpy(phase2, h2->phase);
  gsl_vector_scale(amp_imag2, -1);
  gsl_vector_scale(phase2, -1);
  double *h2Ai = amp_imag2->data;
  double *h2p = phase2->data;

  double overlap = wip_phase(f1, n1, f2, n2, h1Ar, h1Ai, h1p, h2Ar, h2Ai, h2p, Snoise, 1.0);
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

/* Function computing the log likelihood (h|s) - 1/2 (h|h), with s the signal, h the template, and where we discarded the constant term (s|s) - all the cross-products between modes are taken into account */
double FDLogLikelihood(
  struct tagListmodesCAmpPhaseFrequencySeries *lists,  /* Input: list of modes for the signal s, in Frequency-domain amplitude and phase form */
  struct tagListmodesCAmpPhaseFrequencySeries *listh,  /* Input: list of modes for the template, in Frequency-domain amplitude and phase form */
  double (*Snoise)(double)) /* Noise function */
{
  return FDListmodesOverlap(lists, listh, Snoise) - 1./2 * FDListmodesOverlap(listh, listh, Snoise);
}
