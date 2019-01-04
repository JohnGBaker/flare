#ifndef _LIKELIHOODKURZ_H
#define _LIKELIHOODKURZ_H

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
#include "Faddeeva.h"


#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif


/* Function to evaluate a Noise function  */
void EvaluateNoise(
  gsl_vector* noisevalues,                         /* Output: vector of the noise values */
  gsl_vector* freq,                                /* Input: vector of frequencies on which to evaluate */
  ObjectFunction * Snoise,                         /* Noise function */
  double fLow,                                     /* Lower bound of the frequency window for the detector */
  double fHigh);                                   /* Upper bound of the frequency window for the detector */

/* Function building a frequency vector with logarithmic sampling */
void SetLogFrequencies(
  gsl_vector* freqvector,    /* Output pointer to gsl_vector, already allocated */
  const double fmin,         /* Lower bound of the frequency interval */
  const double fmax,         /* Upper bound of the frequency interval */
  const int nbpts);          /* Number of points */
/* Function building a frequency vector with linear sampling */
void SetLinearFrequencies(
  gsl_vector* freqvector,    /* Output pointer to gsl_vector, already allocated */
  const double fmin,         /* Lower bound of the frequency interval */
  const double fmax,         /* Upper bound of the frequency interval */
  const int nbpts);          /* Number of points */



  /* Function computing the overlap (h1|h2) between two given modes, for a given noise function - uses simple trapeze integration on logarithmically sampled frequencies  */
double FDSinglemodeLogLinearOverlap(
  struct tagCAmpPhaseFrequencySeries *freqseries1, /* First mode h1, in amplitude/phase form */
  struct tagCAmpPhaseFrequencySeries *freqseries2, /* Second mode h2, in amplitude/phase form */
  ObjectFunction * Snoise,                         /* Noise function */
  double fLow,                                     /* Lower bound of the frequency window for the detector */
  double fHigh);                                   /* Upper bound of the frequency window for the detector */


/* Function computing the overlap (h1|h2) between two waveforms given as Re/Im frequency series (common freq values), for a given vector of noise values - uses simple trapeze integration */
double FDOverlapReImvsReIm(
  struct tagReImFrequencySeries *h1,  /* First waveform, frequency series in Re/Im form */
  struct tagReImFrequencySeries *h2,  /* Second waveform, frequency series in Re/Im form */
  gsl_vector* noisevalues);           /* Vector for the noise values on common freq of the freqseries */


/* Function computing the overlap (h1|h2) between two given modes in amplitude/phase form, one being already interpolated, for a given noise function - uses the amplitude/phase representation (Fresnel) */
double FDSinglemodeFresnelOverlap(
  struct tagCAmpPhaseFrequencySeries *freqseries1, /* First mode h1, in amplitude/phase form */
  struct tagCAmpPhaseSpline *splines2,             /* Second mode h2, already interpolated in matrix form */
  ObjectFunction * Snoise,                  /* Noise function */
  double fLow,                                     /* Lower bound of the frequency window for the detector */
  double fHigh);                                    /* Upper bound of the frequency window for the detector */


/* Function computing the overlap (h1|h2) between two given modes in amplitude/phase form for each non-correlated channel 1,2,3, one being already interpolated, for a given noise function - uses the amplitude/phase representation (Fresnel) */
double FDSinglemodeFresnelOverlap3Chan(
  struct tagCAmpPhaseFrequencySeries *freqseries1chan1, /* First mode h1 for channel 1, in amplitude/phase form */
  struct tagCAmpPhaseFrequencySeries *freqseries1chan2, /* First mode h1 for channel 2, in amplitude/phase form */
  struct tagCAmpPhaseFrequencySeries *freqseries1chan3, /* First mode h1 for channel 3, in amplitude/phase form */
  struct tagCAmpPhaseSpline *splines2chan1,             /* Second mode h2 for channel 1, already interpolated in matrix form */
  struct tagCAmpPhaseSpline *splines2chan2,             /* Second mode h2 for channel 2, already interpolated in matrix form */
  struct tagCAmpPhaseSpline *splines2chan3,             /* Second mode h2 for channel 3, already interpolated in matrix form */
  ObjectFunction * Snoisechan1,                  /* Noise function */
  ObjectFunction * Snoisechan2,                  /* Noise function */
  ObjectFunction * Snoisechan3,                  /* Noise function */
  double fLow,                                      /* Lower bound of the frequency window for the detector */
  double fHigh);                                     /* Upper bound of the frequency window for the detector */

  /* Function computing the integrand values, combining three non-correlated channels */
int ComputeIntegrandValues3Chan(
  CAmpPhaseFrequencySeries** integrand,     /* Output: values of the integrand on common frequencies (initialized in the function) */
  CAmpPhaseFrequencySeries* freqseries1chan1,    /* Input: frequency series for wf 1, channel 1 */
  CAmpPhaseFrequencySeries* freqseries1chan2,    /* Input: frequency series for wf 1, channel 2 */
  CAmpPhaseFrequencySeries* freqseries1chan3,    /* Input: frequency series for wf 1, channel 3 */
  CAmpPhaseSpline* splines2chan1,                /* Input: splines in matrix form for wf 2, channel 1 */
  CAmpPhaseSpline* splines2chan2,                /* Input: splines in matrix form for wf 2, channel 2 */
  CAmpPhaseSpline* splines2chan3,                /* Input: splines in matrix form for wf 2, channel 3 */
  ObjectFunction * Snoise1,                /* Noise function */
  ObjectFunction * Snoise2,                /* Noise function */
  ObjectFunction * Snoise3,                /* Noise function */
  double fLow,                              /* Lower bound of the frequency - 0 to ignore */
  double fHigh);                             /* Upper bound of the frequency - 0 to ignore */

/* Function computing the integrand values */
void ComputeIntegrandValues(
  CAmpPhaseFrequencySeries** integrand,     /* Output: values of the integrand on common frequencies (initialized in the function) */
  CAmpPhaseFrequencySeries* freqseries1,    /* Input: frequency series for wf 1 */
  CAmpPhaseSpline* splines2,                /* Input: splines in matrix form for wf 2 */
  ObjectFunction * Snoise,           /* Noise function */
  double fLow,                              /* Lower bound of the frequency - 0 to ignore */
  double fHigh);                            /* Upper bound of the frequency - 0 to ignore */


/* Function computing the overlap (h1|h2) between two given modes in amplitude/phase form, for a given noise function - uses the amplitude/phase representation (wip) */
double FDSinglemodeWIPOverlap(
  struct tagCAmpPhaseFrequencySeries *freqseries1, /* First mode h1, in amplitude/phase form */
  struct tagCAmpPhaseFrequencySeries *freqseries2, /* Second mode h2, in amplitude/phase form */
  ObjectFunction * Snoise,                         /* Noise function */
  double fLow,                                     /* Lower bound of the frequency window for the detector */
  double fHigh);                                   /* Upper bound of the frequency window for the detector */


double complex wip_phase (double *f1,
  int n1,
  double *f2,
  int n2,
  double *s1Ar,
  double *s1Ai,
  double  *s1p,
  double *s2Ar,
  double*s2Ai,
  double *s2p,
  ObjectFunction *Snoise,
  // double (*Snoise)(double),
  double scalefactor,
  double min_f,
  double max_f);

double complex wip_adaptive_phase (double *f1, int n1, double *f2, int n2, double *s1Ar, double *s1Ai, double  *s1p, double *s2Ar, double*s2Ai, double *s2p, double (*Snoise)(double), double scalefactor, int downsample,  double errtol, double min_f, double max_f);

void spline_construct(const double xs[],const double ys[],double zs[],int n);
double spline_int(double xp, const double xs[],const double ys[],const double zs[],int n);
void spline_intd3(double xp,
  const double xs[],
  const double ys[],
  const double zs[],
  int n,
  double *q,
  double *dq1,
  double *dq2,
  double *dq3);

// void spline_intd3_vec(const double xp[],int nxp, const double xs[],const double ys[],const double zs[],int n,double q[], double dq1[], double dq2[], double dq3[]);
int spline_findix(double xp, const double xs[],int n);

double complex cexpm1i(double zi);

///This is just the one-step part of the integration routine abstracted
double complex compute_int(
  double fleft,
  double fright,
  double *fs,
  int nf,
  double *Ars,
  double *Arz,
  double *Ais,
  double *Aiz,
  double *dphis,
  double *dphiz);


//Interpolation of integrand on a new f-grid
///If min_f and max_f are nonpositive then the integrand is defined over the intersection of the f1 and f2 domains.
///If either of these is greater than zero, then that value is used for the repsective bound.
void interpolate_ip_integrand(double *f1,
  double *s1Ar,
  double *s1Ai,
  double  *s1p,
  int n1,
  double *f2,
  double *s2Ar,
  double*s2Ai,
  double *s2p,
  int n2,
  ObjectFunction *Snoise,
  // double (*Snoise)(double),
  double scalefactor,
  double *fnew,
  double *Ars,
  double *Ais,
  double *dphis,
  int *nnew,
  double min_f,
  double max_f);

double SnXProposal(double f);
double SnTProposal(double f);
double SnEProposal(double f);
double SnAProposal(double f);


#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LIKELIHOOD_H */
