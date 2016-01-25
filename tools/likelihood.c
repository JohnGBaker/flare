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
#include "splinecoeffs.h"
#include "fresnel.h"
#include "likelihood.h"

#include "wip.h"

#include <time.h> /* for testing */


/* Number of points to be used in linear integration - hardcoded for now */
#define nbptsintdefault 32768 /* Default number of points to use for linear overlaps */


/********************************* Utilities ****************************************/

/* Newtonian estimate of the relation Mf(deltat/M) (for the 22 mode) - gives the starting geometric frequency for a given mass ratio and a given geometric duration of the observations */
double NewtonianfoftGeom(
  const double q,                      /* Mass ratio m1/m2 */
  const double t)                      /* Duration of the observations in geometric units (t/M) */
{
  double nu = q/(1.+q)/(1.+q);
  return 1./PI * pow(256*nu/5. * t, -3./8);
}

/* Function to evaluate a Noise function  */
void EvaluateNoise(
  gsl_vector* noisevalues,                         /* Output: vector of the noise values */
  gsl_vector* freq,                                /* Input: vector of frequencies on which to evaluate */
  double (*Snoise)(double),                        /* Noise function */
  double fLow,                                     /* Lower bound of the frequency window for the detector */
  double fHigh)                                    /* Upper bound of the frequency window for the detector */
{
  int nbpts = (int) freq->size;

  /* Checking the length */
  if( freq->size != noisevalues->size) {
    printf("Error: incompatible sizes in EvaluateNoise.\n");
    exit(1);
  }
  
  /* Checking the boundary frequencies */
  if( gsl_vector_get(freq, 0) < fLow || gsl_vector_get(freq, nbpts-1) > fHigh ) {
    printf("Error: incompatible frequency range in EvaluateNoise.\n");
    exit(1);
  }

  for(int i=0; i<nbpts; i++) {
    gsl_vector_set(noisevalues, i, Snoise(gsl_vector_get(freq, i)));
  }
}

/* Function building a frequency vector with linear or logarithmic sampling */
void SetLinearFrequencies(
  gsl_vector* freqvector,    /* Output pointer to gsl_vector, already allocated */
  const double fmin,         /* Lower bound of the frequency interval */
  const double fmax,         /* Upper bound of the frequency interval */
  const int nbpts)           /* Number of points */
{  
  /* Vector of frequencies with logarithmic spacing */
  double stepf = (fmax-fmin)/(nbpts-1.);
  for(int i=0; i<nbpts; i++) {
    gsl_vector_set(freqvector, i, fmin + i * stepf);
  }
}
void SetLogFrequencies(
  gsl_vector* freqvector,    /* Output pointer to gsl_vector, already allocated */
  const double fmin,         /* Lower bound of the frequency interval */
  const double fmax,         /* Upper bound of the frequency interval */
  const int nbpts)           /* Number of points */
{  
  /* Vector of frequencies with logarithmic spacing */
  double lnratio = log(fmax/fmin);
  double lnfmin = log(fmin);
  for(int i=0; i<nbpts; i++) {
    gsl_vector_set(freqvector, i, exp(lnfmin + (double) i/(nbpts-1.) * lnratio));
  }
}

/* Function determining logarithmically spaced frequencies from a waveform given as a list of modes, a minimal frequency and a number of points */
void ListmodesSetFrequencies(
  struct tagListmodesCAmpPhaseFrequencySeries *list,     /* Waveform, list of modes in amplitude/phase form */
  double fLow,                                           /* Additional lower frequency limit */
  double fHigh,                                          /* Additional upper frequency limit */
  int nbpts,                                             /* Number of frequency samples */
  int tagsampling,                                       /* Tag for linear (0) or logarithmic (1) sampling */
  gsl_vector* freqvector)                                /* Output: vector of frequencies, already allocated with nbpts */
{
  /* Checking the length */
  if((int) freqvector->size != nbpts) {
    printf("Error: incompatible sizes in ListmodesSetFrequencies.\n");
    exit(1);
  }
  
  /* Determining the frequency interval - from the lowest frequency of the 22 mode to the highest frequency covered by at least one mode */
  double minf, maxf;
  ListmodesCAmpPhaseFrequencySeries* listelementmode22 = ListmodesCAmpPhaseFrequencySeries_GetMode(list, 2, 2);
  minf = gsl_vector_get(listelementmode22->freqseries->freq, 0);
  maxf = gsl_vector_get(list->freqseries->freq, (int) list->freqseries->freq->size - 1);
  ListmodesCAmpPhaseFrequencySeries* listelement = list;
  while(listelement) {
    maxf = fmax(maxf, gsl_vector_get(listelement->freqseries->freq, (int) listelement->freqseries->freq->size - 1));
    listelement = listelement->next;
  }

  /* Checking that the waveform covers the fLow */
  if(fLow < minf) printf("Warning: in ListmodesSetFrequencies, fLow not covered by the wave data.\n");

  /* Actual boundaries of the frequency vector - cuts what is below fLow and above fHigh */
  minf = fmax(minf, fLow);
  maxf = fmin(maxf, fHigh);

  /* Setting values of the vector of frequencies, linear or logarithmic sampling */
  if(tagsampling==0) SetLinearFrequencies(freqvector, minf, maxf, nbpts);
  else if(tagsampling==1) SetLogFrequencies(freqvector, minf, maxf, nbpts);
  else {
    printf("Error: incorrect tagsampling in ListmodesSetFrequencies.");
    exit(1);
  }
}

/* Function computing a simple trapeze integration for real data */
static double TrapezeIntegrate(const gsl_vector* x, const gsl_vector* y)
{
  if(x->size!=y->size) {
    printf("Error: trying to apply TrapezeIntegrate on vectors of different lengths.\n");
    exit(1);
  }
  int N = ((int) x->size) - 1;
  double result = 0.;
  for(int i=0; i<N; i++){
    result += (gsl_vector_get(x, i+1) - gsl_vector_get(x, i)) * (gsl_vector_get(y, i) + gsl_vector_get(y, i+1))/2.;
  }
  return result;
}

/***************************** Functions for overlaps using linear integration ******************************/

/* Function computing the overlap (h1|h2) between two given modes, for a given noise function - uses simple trapeze integration on logarithmically sampled frequencies  */
double FDSinglemodeLogLinearOverlap(
  struct tagCAmpPhaseFrequencySeries *freqseries1, /* First mode h1, in amplitude/phase form */
  struct tagCAmpPhaseFrequencySeries *freqseries2, /* Second mode h2, in amplitude/phase form */
  double (*Snoise)(double),                        /* Noise function */
  double fLow,                                     /* Lower bound of the frequency window for the detector */
  double fHigh)                                    /* Upper bound of the frequency window for the detector */
{
  double res;
  int size1 = (int) freqseries1->freq->size;
  int size2 = (int) freqseries2->freq->size;
  int nbpts = nbptsintdefault;

  /* Minimal, maximal frequencies */
  double fmin1 = gsl_vector_get(freqseries1->freq, 0);
  double fmin2 = gsl_vector_get(freqseries2->freq, 0);
  double fmax1 = gsl_vector_get(freqseries1->freq, size1 - 1);
  double fmax2 = gsl_vector_get(freqseries2->freq, size2 - 1);
  double fmin0 = fmax(fLow, fmax(fmin1, fmin2));
  double fmax0 = fmin(fHigh, fmin(fmax1, fmax2));

  /* Vector of frequencies with logarithmic spacing */
  gsl_vector* freqvector = gsl_vector_alloc(nbpts);
  SetLogFrequencies(freqvector, fmin0, fmax0, nbpts);

  /* Initializing the splines */
  /* Note: since this must also apply to mode contribution after processing, real and imaginary parts of the amplitude are present - but since they should differ for LLV detectors by a constant factor, they can be interpolated */
  gsl_interp_accel* accel_amp1real = gsl_interp_accel_alloc();
  gsl_interp_accel* accel_amp2real = gsl_interp_accel_alloc();
  gsl_interp_accel* accel_amp1imag = gsl_interp_accel_alloc();
  gsl_interp_accel* accel_amp2imag = gsl_interp_accel_alloc();
  gsl_interp_accel* accel_phase1 = gsl_interp_accel_alloc();
  gsl_interp_accel* accel_phase2 = gsl_interp_accel_alloc();
  gsl_spline* amp1real = gsl_spline_alloc(gsl_interp_cspline, size1);
  gsl_spline* amp2real = gsl_spline_alloc(gsl_interp_cspline, size2);
  gsl_spline* amp1imag = gsl_spline_alloc(gsl_interp_cspline, size1);
  gsl_spline* amp2imag = gsl_spline_alloc(gsl_interp_cspline, size2);
  gsl_spline* phase1 = gsl_spline_alloc(gsl_interp_cspline, size1);
  gsl_spline* phase2 = gsl_spline_alloc(gsl_interp_cspline, size2);
  gsl_vector* valuesvector = gsl_vector_alloc(nbpts);
  gsl_spline_init(amp1real, gsl_vector_const_ptr(freqseries1->freq,0), gsl_vector_const_ptr(freqseries1->amp_real,0), size1);
  gsl_spline_init(amp1imag, gsl_vector_const_ptr(freqseries1->freq,0), gsl_vector_const_ptr(freqseries1->amp_imag,0), size1);
  gsl_spline_init(amp2real, gsl_vector_const_ptr(freqseries2->freq,0), gsl_vector_const_ptr(freqseries2->amp_real,0), size2);
  gsl_spline_init(amp2imag, gsl_vector_const_ptr(freqseries2->freq,0), gsl_vector_const_ptr(freqseries2->amp_imag,0), size2);
  gsl_spline_init(phase1, gsl_vector_const_ptr(freqseries1->freq,0), gsl_vector_const_ptr(freqseries1->phase,0), size1);
  gsl_spline_init(phase2, gsl_vector_const_ptr(freqseries2->freq,0), gsl_vector_const_ptr(freqseries2->phase,0), size2);

  /* Main loop - vector of values to be evaluated */
  double f, phi1, phi2, Sn;
  double complex A1;
  double complex A2;
  double* freqvectordata = freqvector->data;
  double* valuesvectordata = valuesvector->data;
  int i=0;
  for(i=0; i<nbpts; i++){
    if(i==0) {f = fmax(freqvectordata[i], fmin0);}
    else if(i==nbpts-1) {f = fmin(freqvectordata[i], fmax0);}
    else {f = freqvectordata[i];}
    A1 = gsl_spline_eval(amp1real, f, accel_amp1real) + I*gsl_spline_eval(amp1imag, f, accel_amp1imag);
    A2 = gsl_spline_eval(amp2real, f, accel_amp2real) + I*gsl_spline_eval(amp2imag, f, accel_amp2imag);
    phi1 = gsl_spline_eval(phase1, f, accel_phase1);
    phi2 = gsl_spline_eval(phase2, f, accel_phase2);
    Sn = Snoise(f);
    valuesvectordata[i] = 4.*creal(A1*conj(A2)*cexp(I*(phi1-phi2))/Sn);
  }

  /* Trapeze integration */
  res = TrapezeIntegrate(freqvector, valuesvector);

  /* Clean up */
  gsl_vector_free(freqvector);
  gsl_vector_free(valuesvector);
  gsl_interp_accel_free(accel_amp1real);
  gsl_interp_accel_free(accel_amp2real);
  gsl_interp_accel_free(accel_amp1imag);
  gsl_interp_accel_free(accel_amp2imag);
  gsl_interp_accel_free(accel_phase1);
  gsl_interp_accel_free(accel_phase2);
  gsl_spline_free(amp1real);
  gsl_spline_free(amp2real);
  gsl_spline_free(amp1imag);
  gsl_spline_free(amp2imag);
  gsl_spline_free(phase1);
  gsl_spline_free(phase2);

  return res;
}

/* Function computing the overlap (h1|h2) between two waveforms given as lists of mode contributions (factos sYlm already included), for a given noise function - uses simple trapeze integration on logarithmically sampled frequencies - generates the frequency series in Re/Im form by summing the mode contributions first, then computes the overlap */
double FDListmodesLogLinearOverlap(
  struct tagListmodesCAmpPhaseFrequencySeries *list1,    /* First waveform, list of modes in amplitude/phase form */
  struct tagListmodesCAmpPhaseFrequencySeries *list2,    /* Second waveform, list of modes in amplitude/phase form */
  double (*Snoise)(double),                              /* Noise function */
  double fLow,                                           /* Lower bound of the frequency window for the detector */
  double fHigh,                                          /* Upper bound of the frequency window for the detector */
  double fstartobs1,                                     /* Starting frequency for the 22 mode of wf 1 - as determined from a limited duration of the observation - set to 0 to ignore */
  double fstartobs2)                                     /* Starting frequency for the 22 mode of wf 2 - as determined from a limited duration of the observation - set to 0 to ignore */
{
  /* Number of points to use in the trapeze integration */
  int nbpts = nbptsintdefault;

  /* Determining the frequency interval - from the lowest frequency of the 22 mode to the highest frequency covered by at least one mode */
  double minf1, maxf1, minf2, maxf2, minf, maxf;
  ListmodesCAmpPhaseFrequencySeries* listelement1mode22 = ListmodesCAmpPhaseFrequencySeries_GetMode(list1, 2, 2);
  ListmodesCAmpPhaseFrequencySeries* listelement2mode22 = ListmodesCAmpPhaseFrequencySeries_GetMode(list2, 2, 2);
  minf1 = fmax(gsl_vector_get(listelement1mode22->freqseries->freq, 0), fstartobs1);
  minf2 = fmax(gsl_vector_get(listelement2mode22->freqseries->freq, 0), fstartobs2);
  maxf1 = gsl_vector_get(list1->freqseries->freq, (int) list1->freqseries->freq->size - 1);
  maxf2 = gsl_vector_get(list2->freqseries->freq, (int) list2->freqseries->freq->size - 1);
  ListmodesCAmpPhaseFrequencySeries* listelement1 = list1;
  ListmodesCAmpPhaseFrequencySeries* listelement2 = list2;
  while(listelement1) {
    maxf1 = fmax(maxf1, gsl_vector_get(listelement1->freqseries->freq, (int) listelement1->freqseries->freq->size - 1));
    listelement1 = listelement1->next;
  }
  while(listelement2) {
    maxf2 = fmax(maxf2, gsl_vector_get(listelement2->freqseries->freq, (int) listelement2->freqseries->freq->size - 2));
    listelement2 = listelement2->next;
  }
  /* Actual boundaries to be used in the overlap - also cuts what is not covered by the noise data */
  minf = fmax(fmax(minf1, minf2), fLow);
  maxf = fmin(fmin(maxf1, maxf2), fHigh);

  /* Vector of frequencies used for the overlap */
  gsl_vector* freqoverlap = gsl_vector_alloc(nbpts);
  SetLogFrequencies(freqoverlap, minf, maxf, nbpts);
  
  /* Evaluating each frequency series by interpolating and summing the mode contributions */
  ReImFrequencySeries* freqseries1 = NULL;
  ReImFrequencySeries_Init(&freqseries1, nbpts);
  ReImFrequencySeries_SumListmodesCAmpPhaseFrequencySeries(freqseries1, list1, freqoverlap, fstartobs1);
  ReImFrequencySeries* freqseries2 = NULL;
  ReImFrequencySeries_Init(&freqseries2, nbpts);
  ReImFrequencySeries_SumListmodesCAmpPhaseFrequencySeries(freqseries2, list2, freqoverlap, fstartobs2);

  /* Compute the integrand */
  gsl_vector* valuesoverlap = gsl_vector_alloc(nbpts);
  double* hreal1data = freqseries1->h_real->data;
  double* himag1data = freqseries1->h_imag->data;
  double* hreal2data = freqseries2->h_real->data;
  double* himag2data = freqseries2->h_imag->data;
  double* freqdata = freqoverlap->data;
  for(int i=0; i<nbpts; i++) {
    gsl_vector_set(valuesoverlap, i, 4.*creal( (hreal1data[i] + I*himag1data[i]) * (hreal2data[i] - I*himag2data[i]) / Snoise(freqdata[i])));
  }
  
  /* Final trapeze integration */
  double overlap = TrapezeIntegrate(freqoverlap, valuesoverlap);

  /* Clean up */
  ReImFrequencySeries_Cleanup(freqseries1);
  ReImFrequencySeries_Cleanup(freqseries2);
  gsl_vector_free(freqoverlap);
  gsl_vector_free(valuesoverlap);

  return overlap;
}

/* Function computing the overlap (h1|h2) between h1 given in Re/Im form and h2 given as a list of mode contributions (factors sYlm already included), for a given vector of noise values - uses simple trapeze integration - generates the frequency series of h2 in Re/Im form by summing the mode contributions first, then computes the overlap */
double FDOverlapReImvsListmodesCAmpPhase(
  struct tagReImFrequencySeries *freqseries1,            /* First waveform, in Re/Im form */
  struct tagListmodesCAmpPhaseFrequencySeries *list2,    /* Second waveform, list of modes in amplitude/phase form */
  gsl_vector* noisevalues,                               /* Vector for the noise values on the freq of h1 */
  double fstartobs2)                                     /* Starting frequency for the 22 mode of wf 2 - as determined from a limited duration of the observation - set to 0 to ignore */
{
  /* Check the lengths */
  if(freqseries1->freq->size != noisevalues->size) {
    printf("Error: inconsistent lengths in FDOverlapReImvsListmodesCAmpPhase.\n");
    exit(1);
  }
  
  /* Frequencies used for the overlap */
  int nbpts = (int) freqseries1->freq->size;
  gsl_vector* freqoverlap = freqseries1->freq;
  
  /* Evaluating frequency series 2 by interpolating and summing the mode contributions */
  ReImFrequencySeries* freqseries2 = NULL;
  ReImFrequencySeries_Init(&freqseries2, nbpts);
  ReImFrequencySeries_SumListmodesCAmpPhaseFrequencySeries(freqseries2, list2, freqoverlap, fstartobs2);

  /* Compute the integrand */
  gsl_vector* valuesoverlap = gsl_vector_alloc(nbpts);
  double* hreal1data = freqseries1->h_real->data;
  double* himag1data = freqseries1->h_imag->data;
  double* hreal2data = freqseries2->h_real->data;
  double* himag2data = freqseries2->h_imag->data;
  double* noisedata = noisevalues->data;
  for(int i=0; i<nbpts; i++) {
    gsl_vector_set(valuesoverlap, i, 4.*creal( (hreal1data[i] + I*himag1data[i]) * (hreal2data[i] - I*himag2data[i]) / noisedata[i]));
  }
  
  /* Final trapeze integration */
  double overlap = TrapezeIntegrate(freqoverlap, valuesoverlap);

  /* Clean up */
  ReImFrequencySeries_Cleanup(freqseries2);
  gsl_vector_free(valuesoverlap);

  return overlap;
}

/* Function computing the overlap (h1|h2) between two waveforms given as Re/Im frequency series (common freq values assumed), for a given vector of noise values - uses simple trapeze integration */
double FDOverlapReImvsReIm(
  struct tagReImFrequencySeries *freqseries1,  /* First waveform, frequency series in Re/Im form */
  struct tagReImFrequencySeries *freqseries2,  /* Second waveform, frequency series in Re/Im form */
  gsl_vector* noisevalues)                     /* Vector for the noise values on common freq of the freqseries */
{
  /* Check the lengths */
  if(freqseries1->freq->size != noisevalues->size || freqseries2->freq->size != noisevalues->size) {
    printf("Error: inconsistent lengths in FDOverlapReImvsReIm.\n");
    exit(1);
  }

  /* Frequency vector - assuming they match beyond their mere lengths */
  gsl_vector* freqoverlap = freqseries1->freq;
  int nbpts = (int) freqoverlap->size;

  /* Compute the integrand */
  gsl_vector* valuesoverlap = gsl_vector_alloc((int) freqoverlap->size);
  double* hreal1data = freqseries1->h_real->data;
  double* himag1data = freqseries1->h_imag->data;
  double* hreal2data = freqseries2->h_real->data;
  double* himag2data = freqseries2->h_imag->data;
  double* noisedata = noisevalues->data;
  for(int i=0; i<nbpts; i++) {
    gsl_vector_set(valuesoverlap, i, 4.*creal( (hreal1data[i] + I*himag1data[i]) * (hreal2data[i] - I*himag2data[i]) / noisedata[i]));
  }

  //
  Write_Text_Vector("/Users/marsat/src/flare/test/testcompareSNRtoFFT", "testfreqoverlap.txt", freqoverlap);
  Write_Text_Vector("/Users/marsat/src/flare/test/testcompareSNRtoFFT", "testvaluesoverlap.txt", valuesoverlap);
  
  /* Final trapeze integration */
  double overlap = TrapezeIntegrate(freqoverlap, valuesoverlap);

  /* Clean up */
  gsl_vector_free(valuesoverlap);

  return overlap;
}

/***************************** Functions for overlaps using amplitude/phase (wip) ******************************/

/* Function computing the overlap (h1|h2) between two given modes in amplitude/phase form, for a given noise function - uses the amplitude/phase representation (wip) */
double FDSinglemodeWIPOverlap(
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
  /* NOTE: factor 4 was previously missing */
  double overlap = 4.*wip_phase(f1, n1, f2, n2, h1Ar, h1Ai, h1p, h2Ar, h2Ai, h2p, Snoise, 1.0, fLow, fHigh);
  return overlap;
}

/* Function computing the overlap (h1|h2) between two waveforms given as list of modes, for a given noise function - two additional parameters for the starting 22-mode frequencies (then properly scaled for the other modes) for a limited duration of the observations */
double FDListmodesWIPOverlap(
  struct tagListmodesCAmpPhaseFrequencySeries *listh1, /* First waveform, list of modes in amplitude/phase form */
  struct tagListmodesCAmpPhaseFrequencySeries *listh2, /* Second waveform, list of modes in amplitude/phase form */
  double (*Snoise)(double),                            /* Noise function */
  double fLow,                                         /* Lower bound of the frequency window for the detector */
  double fHigh,                                        /* Upper bound of the frequency window for the detector */
  double fstartobs1,                                   /* Starting frequency for the 22 mode of wf 1 - as determined from a limited duration of the observation - set to 0 to ignore */
  double fstartobs2)                                   /* Starting frequency for the 22 mode of wf 2 - as determined from a limited duration of the observation - set to 0 to ignore */
{
  double overlap = 0;

  /* Main loop over the modes - goes through all the modes present */
  ListmodesCAmpPhaseFrequencySeries* listelementh1 = listh1;
  while(listelementh1) {
    ListmodesCAmpPhaseFrequencySeries* listelementh2 = listh2;
    while(listelementh2) {
      /* Scaling fstartobs1/2 with the appropriate factor of m (for the 21 mode we use m=2) - setting fmin in the overlap accordingly */
      int mmax1 = max(2, listelementh1->m);
      int mmax2 = max(2, listelementh2->m);
      double fcutLow = fmax(fLow, fmax(((double) mmax1)/2. * fstartobs1, ((double) mmax2)/2. * fstartobs2));
      overlap += FDSinglemodeWIPOverlap(listelementh1->freqseries, listelementh2->freqseries, Snoise, fcutLow, fHigh);

      listelementh2 = listelementh2->next;
    }
    listelementh1 = listelementh1->next;
  }
  return overlap;
}

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
  int tagint)                                          /* Tag choosing the integrator: 0 for wip, 1 for log linear integration */
{
  double overlap;
  if(tagint==0) {
    overlap = FDListmodesWIPOverlap(listh1, listh2, Snoise, fLow, fHigh, fstartobs1, fstartobs2);
  }
  else if(tagint==1) {
    overlap = FDListmodesLogLinearOverlap(listh1, listh2, Snoise, fLow, fHigh, fstartobs1, fstartobs2);
  }
  return overlap;
}

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
  int tagint)                                          /* Tag choosing the integrator: 0 for wip, 1 for log linear integration */
{
  double lnL;
  if(tagint==0) {
    lnL = FDListmodesWIPOverlap(lists, listh, Snoise, fLow, fHigh, fstartobss, fstartobsh) - 1./2 * hh - 1./2 * ss;
  }
  else if(tagint==1) {
    lnL = FDListmodesLogLinearOverlap(lists, listh, Snoise, fLow, fHigh, fstartobss, fstartobsh) - 1./2 * hh - 1./2 * ss;
  }
  return lnL;
}

/************** Functions for overlap/likelihood specific to loglinear integration *****************/

/* Function computing the log likelihood -1/2(h-s|h-s), with s the signal, h the template (both given as frequency series in Re/Im form) - h, s assumed to be given on the same set of frequencies - same for the vector of noise values - fstartobs for s, h has already been taken into account */
double FDLogLikelihoodReIm(
  struct tagReImFrequencySeries *s,    /* First waveform (injection), frequency series in Re/Im form */
  struct tagReImFrequencySeries *h,    /* Second waveform (template), frequency series in Re/Im form */
  gsl_vector* noisevalues)             /* Vector for the noise values on common freq of the freqseries */
{
  /* Check the lengths */
  if(s->freq->size != noisevalues->size || h->freq->size != noisevalues->size) {
    printf("Error: inconsistent lengths in FDLogLikelihoodReIm.\n");
    exit(1);
  }

  /* Taking the difference between the frequency series: diff = h-s */
  int nbpts = (int) s->freq->size;
  ReImFrequencySeries* diff = NULL;
  ReImFrequencySeries_Init(&diff, nbpts);
  gsl_vector_memcpy(diff->freq, s->freq);
  gsl_vector_memcpy(diff->h_real, h->h_real);
  gsl_vector_memcpy(diff->h_imag, h->h_imag);
  gsl_vector_sub(diff->h_real, s->h_real);
  gsl_vector_sub(diff->h_imag, s->h_imag);

  /* Likelihood lnL = -1/2(h-s|h-s) */
  double lnL;
  lnL = -1./2 * FDOverlapReImvsReIm(diff, diff, noisevalues);

  /* Clean up */
  ReImFrequencySeries_Cleanup(diff);

  return lnL;
}

/***************************** Functions for overlaps using amplitude/phase (Fresnel) ******************************/

/* Note: for the spines in matrix form, the first column contains the x values, so the coeffs start at 1 */
static double EvalCubic(
  gsl_vector* coeffs,  /**/
  double eps,          /**/
  double eps2,         /**/
  double eps3)         /**/
{
  double p0 = gsl_vector_get(coeffs, 1);
  double p1 = gsl_vector_get(coeffs, 2);
  double p2 = gsl_vector_get(coeffs, 3);
  double p3 = gsl_vector_get(coeffs, 4);
  return p0 + p1*eps + p2*eps2 + p3*eps3;
}
static double EvalQuad(
  gsl_vector* coeffs,  /**/
  double eps,          /**/
  double eps2)         /**/
{
  double p0 = gsl_vector_get(coeffs, 1);
  double p1 = gsl_vector_get(coeffs, 2);
  double p2 = gsl_vector_get(coeffs, 3);
  return p0 + p1*eps + p2*eps2;
}

/* Quadratic Legendre approximation to compute values at minf and maxf when they do not fall on the grid of a freqseries, using the two first/last intervals */
static double EstimateBoundaryLegendreQuad(
  gsl_vector* vectx,   /**/
  gsl_vector* vecty,   /**/
  int j,               /**/
  double xvalue)       /**/
{
  double x0 = 0;
  double x1 = gsl_vector_get(vectx, j+1) - gsl_vector_get(vectx, j);
  double x2 = gsl_vector_get(vectx, j+2) - gsl_vector_get(vectx, j);
  double x = xvalue - gsl_vector_get(vectx, j);
  double y0 = gsl_vector_get(vecty, j);
  double y1 = gsl_vector_get(vecty, j+1);
  double y2 = gsl_vector_get(vecty, j+2);
  if(!(x>=x0 && x<=x2)) {
    printf("Error: value out of bounds in EstimateBoundaryLegendreQuad.\n");
    //
    //printf("%.16e, %.16e, %.16e, %.16e:\n", xvalue, gsl_vector_get(vectx, j), gsl_vector_get(vectx, j+1), gsl_vector_get(vectx, j+2));
    exit(1);
  }
  return y0*(x-x1)*(x-x2)/(x0-x1)/(x0-x2) + y1*(x-x0)*(x-x2)/(x1-x0)/(x1-x2) + y2*(x-x0)*(x-x1)/(x2-x0)/(x2-x1);
}

/* Function computing the integrand values */
void ComputeIntegrandValues(
  CAmpPhaseFrequencySeries** integrand,     /* Output: values of the integrand on common frequencies (initialized in the function) */
  CAmpPhaseFrequencySeries* freqseries1,    /* Input: frequency series for wf 1 */
  CAmpPhaseSpline* splines2,                /* Input: splines in matrix form for wf 2 */
  double (*Snoise)(double),                 /* Noise function */
  double fLow,                              /* Lower bound of the frequency - 0 to ignore */
  double fHigh)                             /* Upper bound of the frequency - 0 to ignore */
{
  //
gsl_set_error_handler(&Err_Handler);

  /* Determining the boundaries of indices */
  gsl_vector* freq1 = freqseries1->freq;
  int imin1 = 0;
  int imax1 = freq1->size - 1;
  double* f1 = freq1->data;
  double f2min = gsl_matrix_get(splines2->quadspline_phase, 0, 0);
  double f2max = gsl_matrix_get(splines2->quadspline_phase, splines2->quadspline_phase->size1 - 1, 0);
  if((fLow>0 && (f1[imax1]<=fLow || f2max<=fLow)) || (fHigh>0 && (f1[imin1]>=fHigh || f2min>=fHigh))) {
    printf("Error: range of frequencies incompatible with fLow, fHigh in IntegrandValues.\n");
    exit(1);
  }
  /* If starting outside, move the ends of the frequency series to be just outside the final minf and maxf */
  double minf = fmax(f1[imin1], f2min);
  double maxf = fmin(f1[imax1], f2max);
  if(fLow>0) {minf = fmax(fLow, minf);}
  if(fHigh>0) {maxf = fmin(fHigh,maxf);}
  while(f1[imin1+1]<=minf) imin1++;
  while(f1[imax1-1]>=maxf) imax1--;
  /* Estimate locally values for freqseries1 at the boundaries */
  double areal1minf = EstimateBoundaryLegendreQuad(freq1, freqseries1->amp_real, imin1, minf);
  double aimag1minf = EstimateBoundaryLegendreQuad(freq1, freqseries1->amp_imag, imin1, minf);
  double phi1minf = EstimateBoundaryLegendreQuad(freq1, freqseries1->phase, imin1, minf);
  double areal1maxf = EstimateBoundaryLegendreQuad(freq1, freqseries1->amp_real, imax1-2, maxf); /* Note the imax1-2 */
  double aimag1maxf = EstimateBoundaryLegendreQuad(freq1, freqseries1->amp_imag, imax1-2, maxf); /* Note the imax1-2 */
  double phi1maxf = EstimateBoundaryLegendreQuad(freq1, freqseries1->phase, imax1-2, maxf); /* Note the imax1-2 */

  /* Initializing output structure */
  int nbpts = imax1 + 1 - imin1;
  CAmpPhaseFrequencySeries_Init(integrand, nbpts);

  /* Loop computing integrand values */
  gsl_vector* freq = (*integrand)->freq;
  gsl_vector* ampreal = (*integrand)->amp_real;
  gsl_vector* ampimag = (*integrand)->amp_imag;
  gsl_vector* phase = (*integrand)->phase;
  double f, eps, eps2, eps3, ampreal1, ampimag1, phase1, ampreal2, ampimag2, phase2, invSn;
  double complex camp;
  double* areal1 = freqseries1->amp_real->data;
  double* aimag1 = freqseries1->amp_imag->data;
  double* phi1 = freqseries1->phase->data;
  gsl_matrix* splineAreal2 = splines2->spline_amp_real;
  gsl_matrix* splineAimag2 = splines2->spline_amp_imag;
  gsl_matrix* quadsplinephase2 = splines2->quadspline_phase;
  int i2 = 0; int j = 0;
  for(int i=imin1; i<=imax1; i++) {
    /* Distinguish the case where we are at minf or maxf */
    if(i==imin1) {
      f = minf;
      ampreal1 = areal1minf;
      ampimag1 = aimag1minf;
      phase1 = phi1minf;
    }
    else if(i==imax1) {
      f = maxf;
      ampreal1 = areal1maxf;
      ampimag1 = aimag1maxf;
      phase1 = phi1maxf;
    }
    else {
      f = gsl_vector_get(freq1, i);
      ampreal1 = areal1[i];
      ampimag1 = aimag1[i];
      phase1 = phi1[i];
    }
    /* Adjust the index in the spline if necessary and compute */
    while(gsl_matrix_get(splines2->quadspline_phase, i2+1, 0)<f) i2++;
    eps = f - gsl_matrix_get(splines2->quadspline_phase, i2, 0);
    eps2 = eps*eps;
    eps3 = eps2*eps;
    gsl_vector_view coeffsampreal2 = gsl_matrix_row(splineAreal2, i2);
    gsl_vector_view coeffsampimag2 = gsl_matrix_row(splineAimag2, i2);
    gsl_vector_view coeffsphase2 = gsl_matrix_row(quadsplinephase2, i2);
    ampreal2 = EvalCubic(&coeffsampreal2.vector, eps, eps2, eps3);
    ampimag2 = EvalCubic(&coeffsampimag2.vector, eps, eps2, eps3);
    phase2 = EvalQuad(&coeffsphase2.vector, eps, eps2);
    invSn = 1./Snoise(f);
    camp = invSn * (ampreal1 + I*ampimag1) * (ampreal2 - I*ampimag2);
    gsl_vector_set(freq, j, f);
    gsl_vector_set(ampreal, j, creal(camp));
    gsl_vector_set(ampimag, j, cimag(camp));
    gsl_vector_set(phase, j, phase1 - phase2);
    j++;
  }
}

/* Function computing the integrand values, combining three non-correlated channels */
void ComputeIntegrandValues3Chan(
  CAmpPhaseFrequencySeries** integrand,     /* Output: values of the integrand on common frequencies (initialized in the function) */
  CAmpPhaseFrequencySeries* freqseries1chan1,    /* Input: frequency series for wf 1, channel 1 */
  CAmpPhaseFrequencySeries* freqseries1chan2,    /* Input: frequency series for wf 1, channel 2 */
  CAmpPhaseFrequencySeries* freqseries1chan3,    /* Input: frequency series for wf 1, channel 3 */
  CAmpPhaseSpline* splines2chan1,                /* Input: splines in matrix form for wf 2, channel 1 */
  CAmpPhaseSpline* splines2chan2,                /* Input: splines in matrix form for wf 2, channel 2 */
  CAmpPhaseSpline* splines2chan3,                /* Input: splines in matrix form for wf 2, channel 3 */
  double (*Snoise1)(double),                /* Noise function for channel 1 */
  double (*Snoise2)(double),                /* Noise function for channel 2 */
  double (*Snoise3)(double),                /* Noise function for channel 3 */
  double fLow,                              /* Lower bound of the frequency - 0 to ignore */
  double fHigh)                             /* Upper bound of the frequency - 0 to ignore */
{
  //
  gsl_set_error_handler(&Err_Handler);

  /* Determining the boundaries of indices - frequency vectors assumed to be the same for channels 1,2,3 */
  gsl_vector* freq1 = freqseries1chan1->freq;
  int imin1 = 0;
  int imax1 = freq1->size - 1;
  double* f1 = freq1->data;
  double f2min = gsl_matrix_get(splines2chan1->quadspline_phase, 0, 0);
  double f2max = gsl_matrix_get(splines2chan1->quadspline_phase, splines2chan1->quadspline_phase->size1 - 1, 0);
  if((fLow>0 && (f1[imax1]<=fLow || f2max<=fLow)) || (fHigh>0 && (f1[imin1]>=fHigh || f2min>=fHigh))) {
    printf("Error: range of frequencies incompatible with fLow, fHigh in IntegrandValues.\n");
    exit(1);
  }
  /* If starting outside, move the ends of the frequency series to be just outside the final minf and maxf */
  double minf = fmax(f1[imin1], f2min);
  double maxf = fmin(f1[imax1], f2max);
  if(fLow>0) {minf = fmax(fLow, minf);}
  if(fHigh>0) {maxf = fmin(fHigh,maxf);}
  while(f1[imin1+1]<=minf) imin1++;
  while(f1[imax1-1]>=maxf) imax1--;
  /* Estimate locally values for freqseries1 at the boundaries - phase vectors assumed to be the same for channels 1,2,3 */
  double areal1chan1minf = EstimateBoundaryLegendreQuad(freq1, freqseries1chan1->amp_real, imin1, minf);
  double aimag1chan1minf = EstimateBoundaryLegendreQuad(freq1, freqseries1chan1->amp_imag, imin1, minf);
  double areal1chan2minf = EstimateBoundaryLegendreQuad(freq1, freqseries1chan2->amp_real, imin1, minf);
  double aimag1chan2minf = EstimateBoundaryLegendreQuad(freq1, freqseries1chan2->amp_imag, imin1, minf);
  double areal1chan3minf = EstimateBoundaryLegendreQuad(freq1, freqseries1chan3->amp_real, imin1, minf);
  double aimag1chan3minf = EstimateBoundaryLegendreQuad(freq1, freqseries1chan3->amp_imag, imin1, minf);
  double phi1minf = EstimateBoundaryLegendreQuad(freq1, freqseries1chan1->phase, imin1, minf);
  double areal1chan1maxf = EstimateBoundaryLegendreQuad(freq1, freqseries1chan1->amp_real, imax1-2, maxf); /* Note the imax1-2 */
  double aimag1chan1maxf = EstimateBoundaryLegendreQuad(freq1, freqseries1chan1->amp_imag, imax1-2, maxf); /* Note the imax1-2 */
  double areal1chan2maxf = EstimateBoundaryLegendreQuad(freq1, freqseries1chan2->amp_real, imax1-2, maxf); /* Note the imax1-2 */
  double aimag1chan2maxf = EstimateBoundaryLegendreQuad(freq1, freqseries1chan2->amp_imag, imax1-2, maxf); /* Note the imax1-2 */
  double areal1chan3maxf = EstimateBoundaryLegendreQuad(freq1, freqseries1chan3->amp_real, imax1-2, maxf); /* Note the imax1-2 */
  double aimag1chan3maxf = EstimateBoundaryLegendreQuad(freq1, freqseries1chan3->amp_imag, imax1-2, maxf); /* Note the imax1-2 */
  double phi1maxf = EstimateBoundaryLegendreQuad(freq1, freqseries1chan1->phase, imax1-2, maxf); /* Note the imax1-2 */

  /* Initializing output structure */
  int nbpts = imax1 + 1 - imin1;
  CAmpPhaseFrequencySeries_Init(integrand, nbpts);

  /* Loop computing integrand values - phases are the same for chan1, chan2 and chan3 */
  gsl_vector* freq = (*integrand)->freq;
  gsl_vector* ampreal = (*integrand)->amp_real;
  gsl_vector* ampimag = (*integrand)->amp_imag;
  gsl_vector* phase = (*integrand)->phase;
  double f, eps, eps2, eps3, ampreal1chan1, ampimag1chan1, ampreal1chan2, ampimag1chan2, ampreal1chan3, ampimag1chan3, phase1, ampreal2chan1, ampimag2chan1, ampreal2chan2, ampimag2chan2, ampreal2chan3, ampimag2chan3, phase2, invSnchan1, invSnchan2, invSnchan3;
  double complex camp;
  double* areal1chan1 = freqseries1chan1->amp_real->data;
  double* aimag1chan1 = freqseries1chan1->amp_imag->data;
  double* areal1chan2 = freqseries1chan2->amp_real->data;
  double* aimag1chan2 = freqseries1chan2->amp_imag->data;
  double* areal1chan3 = freqseries1chan3->amp_real->data;
  double* aimag1chan3 = freqseries1chan3->amp_imag->data;
  double* phi1 = freqseries1chan1->phase->data;
  gsl_matrix* splinechan1real2chan1 = splines2chan1->spline_amp_real;
  gsl_matrix* splinechan1imag2chan1 = splines2chan1->spline_amp_imag;
  gsl_matrix* splinechan1real2chan2 = splines2chan2->spline_amp_real;
  gsl_matrix* splinechan1imag2chan2 = splines2chan2->spline_amp_imag;
  gsl_matrix* splinechan1real2chan3 = splines2chan3->spline_amp_real;
  gsl_matrix* splinechan1imag2chan3 = splines2chan3->spline_amp_imag;
  gsl_matrix* quadsplinephase2 = splines2chan1->quadspline_phase;
  int i2 = 0; int j = 0;
  for(int i=imin1; i<=imax1; i++) {
    /* Distinguish the case where we are at minf or maxf */
    if(i==imin1) {
      f = minf;
      ampreal1chan1 = areal1chan1minf;
      ampimag1chan1 = aimag1chan1minf;
      ampreal1chan2 = areal1chan2minf;
      ampimag1chan2 = aimag1chan2minf;
      ampreal1chan3 = areal1chan3minf;
      ampimag1chan3 = aimag1chan3minf;
      phase1 = phi1minf;
    }
    else if(i==imax1) {
      f = maxf;
      ampreal1chan1 = areal1chan1maxf;
      ampimag1chan1 = aimag1chan1maxf;
      ampreal1chan2 = areal1chan2maxf;
      ampimag1chan2 = aimag1chan2maxf;
      ampreal1chan3 = areal1chan3maxf;
      ampimag1chan3 = aimag1chan3maxf;
      phase1 = phi1maxf;
    }
    else {
      f = gsl_vector_get(freq1, i);
      ampreal1chan1 = areal1chan1[i];
      ampimag1chan1 = aimag1chan1[i];
      ampreal1chan2 = areal1chan2[i];
      ampimag1chan2 = aimag1chan2[i];
      ampreal1chan3 = areal1chan3[i];
      ampimag1chan3 = aimag1chan3[i];
      phase1 = phi1[i];
    }
    /* Adjust the index in the spline if necessary and compute */
    while(gsl_matrix_get(splines2chan1->quadspline_phase, i2+1, 0)<f) i2++;
    eps = f - gsl_matrix_get(splines2chan1->quadspline_phase, i2, 0);
    eps2 = eps*eps;
    eps3 = eps2*eps;
    gsl_vector_view coeffsampreal2chan1 = gsl_matrix_row(splinechan1real2chan1, i2);
    gsl_vector_view coeffsampimag2chan1 = gsl_matrix_row(splinechan1imag2chan1, i2);
    gsl_vector_view coeffsampreal2chan2 = gsl_matrix_row(splinechan1real2chan2, i2);
    gsl_vector_view coeffsampimag2chan2 = gsl_matrix_row(splinechan1imag2chan2, i2);
    gsl_vector_view coeffsampreal2chan3 = gsl_matrix_row(splinechan1real2chan3, i2);
    gsl_vector_view coeffsampimag2chan3 = gsl_matrix_row(splinechan1imag2chan3, i2);
    gsl_vector_view coeffsphase2 = gsl_matrix_row(quadsplinephase2, i2);
    ampreal2chan1 = EvalCubic(&coeffsampreal2chan1.vector, eps, eps2, eps3);
    ampimag2chan1 = EvalCubic(&coeffsampimag2chan1.vector, eps, eps2, eps3);
    ampreal2chan2 = EvalCubic(&coeffsampreal2chan2.vector, eps, eps2, eps3);
    ampimag2chan2 = EvalCubic(&coeffsampimag2chan2.vector, eps, eps2, eps3);
    ampreal2chan3 = EvalCubic(&coeffsampreal2chan3.vector, eps, eps2, eps3);
    ampimag2chan3 = EvalCubic(&coeffsampimag2chan3.vector, eps, eps2, eps3);
    phase2 = EvalQuad(&coeffsphase2.vector, eps, eps2);
    invSnchan1 = 1./Snoise1(f);
    invSnchan2 = 1./Snoise2(f);
    invSnchan3 = 1./Snoise3(f);
    camp = invSnchan1 * (ampreal1chan1 + I*ampimag1chan1) * (ampreal2chan1 - I*ampimag2chan1) + invSnchan2 * (ampreal1chan2 + I*ampimag1chan2) * (ampreal2chan2 - I*ampimag2chan2) + invSnchan3 * (ampreal1chan3 + I*ampimag1chan3) * (ampreal2chan3 - I*ampimag2chan3);
    gsl_vector_set(freq, j, f);
    gsl_vector_set(ampreal, j, creal(camp));
    gsl_vector_set(ampimag, j, cimag(camp));
    gsl_vector_set(phase, j, phase1 - phase2);
    j++;
  }
}

/* Function computing the overlap (h1|h2) between two given modes in amplitude/phase form, one being already interpolated, for a given noise function - uses the amplitude/phase representation (Fresnel) */
double FDSinglemodeFresnelOverlap(
  struct tagCAmpPhaseFrequencySeries *freqseries1, /* First mode h1, in amplitude/phase form */
  struct tagCAmpPhaseSpline *splines2,             /* Second mode h2, already interpolated in matrix form */
  double (*Snoise)(double),                        /* Noise function */
  double fLow,                                     /* Lower bound of the frequency window for the detector */
  double fHigh)                                    /* Upper bound of the frequency window for the detector */
{
  /* Computing the integrand values, on the frequency grid of h1 */
  CAmpPhaseFrequencySeries* integrand = NULL;
  ComputeIntegrandValues(&integrand, freqseries1, splines2, Snoise, fLow, fHigh);

  /* Rescaling the integrand */
  double scaling = 10./gsl_vector_get(integrand->freq, integrand->freq->size-1);
  gsl_vector_scale(integrand->freq, scaling);
  gsl_vector_scale(integrand->amp_real, 1./scaling);
  gsl_vector_scale(integrand->amp_imag, 1./scaling);

  /* Interpolating the integrand */
  CAmpPhaseSpline* integrandspline = NULL;
  BuildSplineCoeffs(&integrandspline, integrand);

  //
  /* printf("Inside FDSinglemodeFresnelOverlap\n"); */
  /* Write_Text_Matrix("/Users/marsat/src/flare/test/testlisaoverlap/temp4", "intAsplineampreal.txt", integrandspline->spline_amp_real); */
  /* Write_Text_Matrix("/Users/marsat/src/flare/test/testlisaoverlap/temp4", "intAsplineampimag.txt", integrandspline->spline_amp_imag); */
  /* Write_Text_Matrix("/Users/marsat/src/flare/test/testlisaoverlap/temp4", "intAsplinephase.txt", integrandspline->quadspline_phase); */

  /* Computing the integral - including here the factor 4 and the real part */
  double overlap = 4.*creal(ComputeInt(integrandspline->spline_amp_real, integrandspline->spline_amp_imag, integrandspline->quadspline_phase));

  /* Clean up */
  CAmpPhaseSpline_Cleanup(integrandspline);
  CAmpPhaseFrequencySeries_Cleanup(integrand);

  return overlap;
}

/* Function computing the overlap (h1|h2) between two given modes in amplitude/phase form for each non-correlated channel 1,2,3, one being already interpolated, for a given noise function - uses the amplitude/phase representation (Fresnel) */
double FDSinglemodeFresnelOverlap3Chan(
  struct tagCAmpPhaseFrequencySeries *freqseries1chan1, /* First mode h1 for channel 1, in amplitude/phase form */
  struct tagCAmpPhaseFrequencySeries *freqseries1chan2, /* First mode h1 for channel 2, in amplitude/phase form */
  struct tagCAmpPhaseFrequencySeries *freqseries1chan3, /* First mode h1 for channel 3, in amplitude/phase form */
  struct tagCAmpPhaseSpline *splines2chan1,             /* Second mode h2 for channel 1, already interpolated in matrix form */
  struct tagCAmpPhaseSpline *splines2chan2,             /* Second mode h2 for channel 2, already interpolated in matrix form */
  struct tagCAmpPhaseSpline *splines2chan3,             /* Second mode h2 for channel 3, already interpolated in matrix form */
  double (*Snoisechan1)(double),                        /* Noise function for channel 1 */
  double (*Snoisechan2)(double),                        /* Noise function for channel 2 */
  double (*Snoisechan3)(double),                        /* Noise function for channel 3 */
  double fLow,                                      /* Lower bound of the frequency window for the detector */
  double fHigh)                                     /* Upper bound of the frequency window for the detector */
{
  /* Computing the integrand values, on the frequency grid of h1 */
  CAmpPhaseFrequencySeries* integrand = NULL;

  /////////
  //ComputeIntegrandValues(&integrand, freqseries1A, splines2A, SnoiseA, fLow, fHigh);
  //ComputeIntegrandValues(&integrand, freqseries1chan2, splines2chan2, Snoisechan2, fLow, fHigh);
  //ComputeIntegrandValues(&integrand, freqseries1T, splines2T, SnoiseT, fLow, fHigh);
  ComputeIntegrandValues3Chan(&integrand, freqseries1chan1, freqseries1chan2, freqseries1chan3, splines2chan1, splines2chan2, splines2chan3, Snoisechan1, Snoisechan2, Snoisechan3, fLow, fHigh);

  /* Rescaling the integrand */
  double scaling = 10./gsl_vector_get(integrand->freq, integrand->freq->size-1);
  gsl_vector_scale(integrand->freq, scaling);
  gsl_vector_scale(integrand->amp_real, 1./scaling);
  gsl_vector_scale(integrand->amp_imag, 1./scaling);

  /* Interpolating the integrand */
  CAmpPhaseSpline* integrandspline = NULL;
  BuildSplineCoeffs(&integrandspline, integrand);

  /* Computing the integral - including here the factor 4 and the real part */
  double overlap = 4.*creal(ComputeInt(integrandspline->spline_amp_real, integrandspline->spline_amp_imag, integrandspline->quadspline_phase));

  //
  //printf("In FDSinglemodeFresnelOverlapAET\n");

  /* Clean up */
  CAmpPhaseSpline_Cleanup(integrandspline);
  CAmpPhaseFrequencySeries_Cleanup(integrand);

  return overlap;
}


/* Function computing the overlap (h1|h2) between two waveforms given as list of modes, one being already interpolated, for a given noise function - two additional parameters for the starting 22-mode frequencies (then properly scaled for the other modes) for a limited duration of the observations */
double FDListmodesFresnelOverlap(
  struct tagListmodesCAmpPhaseFrequencySeries *listh1, /* First waveform, list of modes in amplitude/phase form */
  struct tagListmodesCAmpPhaseSpline *listsplines2,    /* Second waveform, list of modes already interpolated in matrix form */
  double (*Snoise)(double),                            /* Noise function */
  double fLow,                                         /* Lower bound of the frequency window for the detector */
  double fHigh,                                        /* Upper bound of the frequency window for the detector */
  double fstartobs1,                                   /* Starting frequency for the 22 mode of wf 1 - as determined from a limited duration of the observation - set to 0 to ignore */
  double fstartobs2)                                   /* Starting frequency for the 22 mode of wf 2 - as determined from a limited duration of the observation - set to 0 to ignore */
{
  double overlap = 0;

  /* Main loop over the modes - goes through all the modes present */
  ListmodesCAmpPhaseFrequencySeries* listelementh1 = listh1;
  while(listelementh1) {
    ListmodesCAmpPhaseSpline* listelementsplines2 = listsplines2;
    while(listelementsplines2) {
      /* Scaling fstartobs1/2 with the appropriate factor of m (for the 21 mode we use m=2) - setting fmin in the overlap accordingly */
      int mmax1 = max(2, listelementh1->m);
      int mmax2 = max(2, listelementsplines2->m);
      double fcutLow = fmax(fLow, fmax(((double) mmax1)/2. * fstartobs1, ((double) mmax2)/2. * fstartobs2));
      overlap += FDSinglemodeFresnelOverlap(listelementh1->freqseries, listelementsplines2->splines, Snoise, fcutLow, fHigh);

      listelementsplines2 = listelementsplines2->next;
    }
    listelementh1 = listelementh1->next;
  }
  return overlap;
}

/* Function computing the overlap (h1|h2) between two waveforms given as list of modes for each non-correlated channel 1,2,3, one being already interpolated, for a given noise function - two additional parameters for the starting 22-mode frequencies (then properly scaled for the other modes) for a limited duration of the observations */
double FDListmodesFresnelOverlap3Chan(
  struct tagListmodesCAmpPhaseFrequencySeries *listh1chan1, /* First waveform channel channel 1, list of modes in amplitude/phase form */
  struct tagListmodesCAmpPhaseFrequencySeries *listh1chan2, /* First waveform channel channel 2, list of modes in amplitude/phase form */
  struct tagListmodesCAmpPhaseFrequencySeries *listh1chan3, /* First waveform channel channel 3, list of modes in amplitude/phase form */
  struct tagListmodesCAmpPhaseSpline *listsplines2chan1,    /* Second waveform channel channel 1, list of modes already interpolated in matrix form */
  struct tagListmodesCAmpPhaseSpline *listsplines2chan2,    /* Second waveform channel channel 2, list of modes already interpolated in matrix form */
  struct tagListmodesCAmpPhaseSpline *listsplines2chan3,    /* Second waveform channel channel 3, list of modes already interpolated in matrix form */
  double (*Snoise1)(double),                            /* Noise function for channel 1 */
  double (*Snoise2)(double),                            /* Noise function for channel 2 */
  double (*Snoise3)(double),                            /* Noise function for channel 3 */
  double fLow,                                          /* Lower bound of the frequency window for the detector */
  double fHigh,                                         /* Upper bound of the frequency window for the detector */
  double fstartobs1,                                    /* Starting frequency for the 22 mode of wf 1 - as determined from a limited duration of the observation - set to 0 to ignore */
  double fstartobs2)                                    /* Starting frequency for the 22 mode of wf 2 - as determined from a limited duration of the observation - set to 0 to ignore */
{
  double overlap = 0;

  /* Main loop over the modes - goes through all the modes present, the same for all three channels 1,2,3 */
  ListmodesCAmpPhaseFrequencySeries* listelementh1chan1 = listh1chan1;
  while(listelementh1chan1) { /* We use the structure for channel 1 to loop through modes */
    ListmodesCAmpPhaseFrequencySeries* listelementh1chan2 = ListmodesCAmpPhaseFrequencySeries_GetMode(listh1chan2, listelementh1chan1->l, listelementh1chan1->m);
    ListmodesCAmpPhaseFrequencySeries* listelementh1chan3 = ListmodesCAmpPhaseFrequencySeries_GetMode(listh1chan3, listelementh1chan1->l, listelementh1chan1->m);
    ListmodesCAmpPhaseSpline* listelementsplines2chan1 = listsplines2chan1;
    while(listelementsplines2chan1) { /* We use the structure for channel 1 to loop through modes */
      ListmodesCAmpPhaseSpline* listelementsplines2chan2 = ListmodesCAmpPhaseSpline_GetMode(listsplines2chan2, listelementsplines2chan1->l, listelementsplines2chan1->m);
      ListmodesCAmpPhaseSpline* listelementsplines2chan3 = ListmodesCAmpPhaseSpline_GetMode(listsplines2chan3, listelementsplines2chan1->l, listelementsplines2chan1->m);
      /* Scaling fstartobs1/2 with the appropriate factor of m (for the 21 mode we use m=2) - setting fmin in the overlap accordingly */
      int mmax1 = max(2, listelementh1chan1->m);
      int mmax2 = max(2, listelementsplines2chan1->m);
      double fcutLow = fmax(fLow, fmax(((double) mmax1)/2. * fstartobs1, ((double) mmax2)/2. * fstartobs2));
      double overlapmode = FDSinglemodeFresnelOverlap3Chan(listelementh1chan1->freqseries, listelementh1chan2->freqseries, listelementh1chan3->freqseries, listelementsplines2chan1->splines, listelementsplines2chan2->splines, listelementsplines2chan3->splines, Snoise1, Snoise2, Snoise3, fcutLow, fHigh);
      overlap += overlapmode;
      listelementsplines2chan1 = listelementsplines2chan1->next;
    }
    listelementh1chan1 = listelementh1chan1->next;
  }

  return overlap;
}
