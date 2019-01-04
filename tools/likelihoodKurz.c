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
#include "likelihoodKurz.h"
#include <time.h> /* for testing */


/* Number of points to be used in linear integration - hardcoded for now */
#define nbptsintdefault 32768 /* Default number of points to use for linear overlaps */

#define false 0
#define true 1

int spline_natural=false;
int spline_iold=0;
int verbose=0;

int wip_relative=1;
double wip_mindf=1e-8;
double wip_errtol=1.0e-8;
int wip_uselogf=0;
int wip_count;
int wip_adapt_verbose=0;



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

// To avoid loss of precision at small z, cexpm1i(zi) implements cos(zi)-1 + i*sin(zi)
// The real analogue of this function is part of the standard library, and the C99 standard
// reserves the name, but it is not generally implemented
// inline double complex cexpm1i(double zi){
double complex cexpm1i(double zi){
  double rr;//rr=cos(x)-1=2*sin(x/2)**2
  double ri=sin(zi);
  double sinhalf=sin(zi/2.0);
  rr=-2.0*sinhalf*sinhalf;
  return rr+I*ri;
};



/* Function to evaluate a Noise function  */
void EvaluateNoise(
  gsl_vector* noisevalues,                         /* Output: vector of the noise values */
  gsl_vector* freq,                                /* Input: vector of frequencies on which to evaluate */
  ObjectFunction * Snoise,                  /* Noise function */
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
  if( gsl_vector_get(freq, 0) - fLow < -1e-15 || gsl_vector_get(freq, nbpts-1) > fHigh + 1e-15 ) {
    printf("Error: incompatible frequency range in EvaluateNoise.\n");
    printf("freq[0]=%g vs fLow=%g, freq[max]=%g vs fHigh=%g\n",gsl_vector_get(freq, 0), fLow, gsl_vector_get(freq, nbpts-1), fHigh);
    printf(" %i, %i\n",gsl_vector_get(freq, 0) < fLow , gsl_vector_get(freq, nbpts-1) > fHigh);
    printf("%g\n",gsl_vector_get(freq, 0) - fLow);
    exit(1);
  }

  for(int i=0; i<nbpts; i++) {
    gsl_vector_set(noisevalues, i, ObjectFunctionCall(Snoise,gsl_vector_get(freq, i)));
  }
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
    exit(1);
  }
  return y0*(x-x1)*(x-x2)/(x0-x1)/(x0-x2) + y1*(x-x0)*(x-x2)/(x1-x0)/(x1-x2) + y2*(x-x0)*(x-x1)/(x2-x0)/(x2-x1);
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



/***************************** Functions for overlaps using linear integration ******************************/

/* Function computing the overlap (h1|h2) between two given modes, for a given noise function - uses simple trapeze integration on logarithmically sampled frequencies  */
double FDSinglemodeLogLinearOverlap(
  struct tagCAmpPhaseFrequencySeries *freqseries1, /* First mode h1, in amplitude/phase form */
  struct tagCAmpPhaseFrequencySeries *freqseries2, /* Second mode h2, in amplitude/phase form */
  ObjectFunction * Snoise,                  /* Noise function */
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
  // SetLogFrequencies(freqvector, fmin0, fmax0, nbpts);
  SetLinearFrequencies(freqvector, fmin0, fmax0, nbpts);

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
  FILE *fp;
  fp = fopen("TestOlap1.dat", "w");
  double Ast1, Ast2;
  for (i=0; i<size1; i++){
    f = gsl_vector_get(freqseries1->freq, i);
    Ast1 =  gsl_vector_get(freqseries1->amp_real, i);
    Ast2 = gsl_vector_get(freqseries1->amp_imag, i);
    Sn = ObjectFunctionCall(Snoise,f);
    // fprintf(fp, "%g   %g  \n", f, 4.0*creal(A1*conj(A1))/Sn );
    fprintf(fp, "%g   %g  \n", f, 4.0*(Ast1*Ast1 + Ast2*Ast2) );
  }
  for(i=0; i<nbpts; i++){
    if(i==0) {f = fmax(freqvectordata[i], fmin0);}
    else if(i==nbpts-1) {f = fmin(freqvectordata[i], fmax0);}
    else {f = freqvectordata[i];}
    A1 = gsl_spline_eval(amp1real, f, accel_amp1real) + I*gsl_spline_eval(amp1imag, f, accel_amp1imag);
    A2 = gsl_spline_eval(amp2real, f, accel_amp2real) + I*gsl_spline_eval(amp2imag, f, accel_amp2imag);
    phi1 = gsl_spline_eval(phase1, f, accel_phase1);
    phi2 = gsl_spline_eval(phase2, f, accel_phase2);
    Sn = ObjectFunctionCall(Snoise,f);
    // valuesvectordata[i] = 4.*creal(A1*conj(A2)*cexp(I*(phi1-phi2))/Sn);
    valuesvectordata[i] = 4.*creal(A1*conj(A2)*cexp(I*(phi1-phi2)))/Sn;
    // fprintf(fp, "%g   %g  \n", f, valuesvectordata[i]);
    // fprintf(fp, "%g   %g  %g \n", f, creal(A2), phi1-phi2);
  }
  fclose(fp);

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

  /* Final trapeze integration */
  double overlap = TrapezeIntegrate(freqoverlap, valuesoverlap);

  /* Clean up */
  gsl_vector_free(valuesoverlap);

  return overlap;
}


/* Function computing the overlap (h1|h2) between two given modes in amplitude/phase form, one being already interpolated, for a given noise function - uses the amplitude/phase representation (Fresnel) */
double FDSinglemodeFresnelOverlap(
  struct tagCAmpPhaseFrequencySeries *freqseries1, /* First mode h1, in amplitude/phase form */
  struct tagCAmpPhaseSpline *splines2,             /* Second mode h2, already interpolated in matrix form */
  ObjectFunction * Snoise,                  /* Noise function */
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
  ObjectFunction * Snoisechan1,                  /* Noise function */
  ObjectFunction * Snoisechan2,                  /* Noise function */
  ObjectFunction * Snoisechan3,                  /* Noise function */
  double fLow,                                      /* Lower bound of the frequency window for the detector */
  double fHigh)                                     /* Upper bound of the frequency window for the detector */
{
  /* Computing the integrand values, on the frequency grid of h1 */
  CAmpPhaseFrequencySeries* integrand = NULL;
  if(0>ComputeIntegrandValues3Chan(&integrand, freqseries1chan1, freqseries1chan2, freqseries1chan3, splines2chan1, splines2chan2, splines2chan3, Snoisechan1, Snoisechan2, Snoisechan3, fLow, fHigh))return 0;//if allowed freq range does not exist, return 0 for overlap
  // ComputeIntegrandValues(&integrand, freqseries1chan1, splines2chan1, Snoisechan1, fLow, fHigh);
  // ComputeIntegrandValues(&integrand, freqseries1chan3, splines2chan3, Snoisechan3, fLow, fHigh);

  /* Rescaling the integrand */
  double scaling = 10./gsl_vector_get(integrand->freq, integrand->freq->size-1);

//TEST
//scaling = 1.;

  gsl_vector_scale(integrand->freq, scaling);
  gsl_vector_scale(integrand->amp_real, 1./scaling);
  gsl_vector_scale(integrand->amp_imag, 1./scaling);

  //dump
  /*
  for(int ii=0;ii<integrand->freq->size;ii++)
    printf("ii=%i, f=%g, integrand = ( %g, %g, %g )\n",ii,gsl_vector_get(integrand->freq,ii),gsl_vector_get(integrand->amp_real,ii),gsl_vector_get(integrand->amp_imag,ii),gsl_vector_get(integrand->phase,ii));
  */

  /* Interpolating the integrand */
  CAmpPhaseSpline* integrandspline = NULL;
  BuildSplineCoeffs(&integrandspline, integrand);

  /* Computing the integral - including here the factor 4 and the real part */
  double overlap = 4.*creal(ComputeInt(integrandspline->spline_amp_real, integrandspline->spline_amp_imag, integrandspline->quadspline_phase));

  /* Clean up */
  CAmpPhaseSpline_Cleanup(integrandspline);
  CAmpPhaseFrequencySeries_Cleanup(integrand);

  return overlap;
}




/* Function computing the integrand values */
void ComputeIntegrandValues(
  CAmpPhaseFrequencySeries** integrand,     /* Output: values of the integrand on common frequencies (initialized in the function) */
  CAmpPhaseFrequencySeries* freqseries1,    /* Input: frequency series for wf 1 */
  CAmpPhaseSpline* splines2,                /* Input: splines in matrix form for wf 2 */
  ObjectFunction * Snoise,           /* Noise function */
  double fLow,                              /* Lower bound of the frequency - 0 to ignore */
  double fHigh)                             /* Upper bound of the frequency - 0 to ignore */
{
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
    printf("need one of {%g, %g} <= %g and one of {%g, %g}>=%g\n",f1[imax1],f2max,fLow,f1[imin1],f2min,fHigh);
    exit(1);
  }
  /* If starting outside, move the ends of the frequency series to be just outside the final minf and maxf */
  double minf = fmax(f1[imin1], f2min);
  double maxf = fmin(f1[imax1], f2max);
  if(fLow>0) {minf = fmax(fLow, minf);}
  if(fHigh>0) {maxf = fmin(fHigh, maxf);}
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
    invSn = 1./ObjectFunctionCall(Snoise,f);
    camp = invSn * (ampreal1 + I*ampimag1) * (ampreal2 - I*ampimag2);
    gsl_vector_set(freq, j, f);
    gsl_vector_set(ampreal, j, creal(camp));
    gsl_vector_set(ampimag, j, cimag(camp));
    gsl_vector_set(phase, j, phase1 - phase2);
    j++;
  }
}




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
  double fHigh)                             /* Upper bound of the frequency - 0 to ignore */
{
  gsl_set_error_handler(&Err_Handler);

  /* Determining the boundaries of indices - frequency vectors assumed to be the same for channels 1,2,3 */
  gsl_vector* freq1 = freqseries1chan1->freq;
  int imin1 = 0;
  int imax1 = freq1->size - 1;
  double* f1 = freq1->data;
  double f2min = gsl_matrix_get(splines2chan1->quadspline_phase, 0, 0);
  double f2max = gsl_matrix_get(splines2chan1->quadspline_phase, splines2chan1->quadspline_phase->size1 - 1, 0);
  if((fLow>0 && (f1[imax1]<=fLow || f2max<=fLow)) || (fHigh>0 && (f1[imin1]>=fHigh || f2min>=fHigh))) {
    //printf("Error: range of frequencies incompatible with fLow, fHigh in IntegrandValues.\n");
    //printf("need both {%g, %g} > %g and both {%g, %g} < %g\n",f1[imax1],f2max,fLow,f1[imin1],f2min,fHigh);
    return -1;
  }
  /* If starting outside, move the ends of the frequency series to be just outside the final minf and maxf */
  double minf = fmax(f1[imin1], f2min);
  double maxf = fmin(f1[imax1], f2max);
  if(fLow>0) {minf = fmax(fLow, minf);}
  if(fHigh>0) {maxf = fmin(fHigh, maxf);}
  while(f1[imin1+1]<=minf) imin1++;
  while(f1[imax1-1]>=maxf) imax1--;
  //printf("imin=%i, imax=%i\n",imin1,imax1);
  int nbpts = imax1 + 1 - imin1;
  //printf("nbpts=%i\n",nbpts);
  if(nbpts<4)return -1;
  /* Estimate locally values for freqseries1 at the boundaries - phase vectors assumed to be the same for channels 1,2,3 - this is still true now that the response-processed phase includes the signal phase + R-delay phase, which is the same for all channels */
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
    invSnchan1 = 1./ObjectFunctionCall(Snoise1,f);
    invSnchan2 = 1./ObjectFunctionCall(Snoise2,f);
    invSnchan3 = 1./ObjectFunctionCall(Snoise3,f);

    camp = invSnchan1 * (ampreal1chan1 + I*ampimag1chan1) * (ampreal2chan1 - I*ampimag2chan1) + invSnchan2 * (ampreal1chan2 + I*ampimag1chan2) * (ampreal2chan2 - I*ampimag2chan2) + invSnchan3 * (ampreal1chan3 + I*ampimag1chan3) * (ampreal2chan3 - I*ampimag2chan3);
    //dump
    /*
    printf("j=%i, im=%g\n a11=(%g,%g),  a12=(%g,%g),  a13=(%g,%g)\n a11=(%g,%g),  a12=(%g,%g),  a13=(%g,%g)\n",j,cimag(camp),
	   ampreal1chan1,ampimag1chan1,ampreal1chan2,ampimag1chan2,ampreal1chan3,ampimag1chan3,
	   ampreal2chan1,ampimag2chan1,ampreal2chan2,ampimag2chan2,ampreal2chan3,ampimag2chan3);
    printf("in1=%g, in2=%g, in3=%g\n",invSnchan1,invSnchan2,invSnchan3);
    */

    gsl_vector_set(freq, j, f);
    gsl_vector_set(ampreal, j, creal(camp));
    gsl_vector_set(ampimag, j, cimag(camp));
    gsl_vector_set(phase, j, phase1 - phase2);
    j++;
  }
  return 0;
}







/* Function computing the overlap (h1|h2) between two given modes in amplitude/phase form, for a given noise function - uses the amplitude/phase representation (wip) */
double FDSinglemodeWIPOverlap(
  struct tagCAmpPhaseFrequencySeries *freqseries1, /* First mode h1, in amplitude/phase form */
  struct tagCAmpPhaseFrequencySeries *freqseries2, /* Second mode h2, in amplitude/phase form */
  ObjectFunction * Snoise,                  /* Noise function */
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

  printf("SB wip: let's try %9.5f, %3.3e \n", 1.e-3, ObjectFunctionCall(Snoise, 1.e-3) );

  /* fLow or fHigh <= 0 means use intersection of signal domains  */
  /* NOTE: factor 4 was previously missing */
  double overlap = 4.*wip_phase(f1, n1, f2, n2, h1Ar, h1Ai, h1p, h2Ar, h2Ai, h2p, Snoise, 1.0, fLow, fHigh);
  return overlap;
}


//Solves tridiagonal for the matched second derivatives z of the 3rd order interpoling polynomials at the sample points (knots).
//The condition that 1st derivatives are continuous at the knots becomes:
//  dx[i-1] z[i-1] + 2(dx[i]+dx[i-1])z[i] + dx[i]z[i+1] = r[i]
//with:
//  dx[i] = x[i+1]-x[i]
// r[i>0] = 6(dy[i]/D[i]-dy[i-1]/D[i])
//  dy[i] = y[i+1]-y[i]
//
//Solution part 1:
//  For 0<i<n-1 defines a tridiagonal matrix  with diagonals b[i],upperdiags c[i] and lower diags a[i]
//  with numbering starting from i=1 in each case
//  All rows are identical but first and last which depend on BCs
//
//  The iterative solution algorithm by LU decomposition is:
//  (with shorthand: bb = b[i] - a[i-1]*v[i-1])
//  v[1] = c[1]/b[1]
//  s[1] = r[1]/b[1]
//  v[i] = c[i]/bb
//  s[i] = (r[i] - a[i-1]*s[i-1])/bb
//
//  In the bulk: 1<i<n-2 the abc values come directly from the eqn at top:
//  a[i]=c[i]=x[i+1]-x[i]
//  b[i]=2*(x[i+1]-x[i-1])
//
//  The boundary conditions for the first and last row come from some condition which
//  eliminates z[0] and z[n-1] in the i=1 and i=n-2 versions of the continuity condition at top
//  This yields the values for b[1],c[1],b[n-2],a[n-3] in the first and last rows.
//
void spline_construct(const double xs[],const double ys[],double zs[],int n){
  //xs and ys are passed in zs are passed out
  //they all must be same length n>2
  double v[n],s[n];
  double dxm,dx,dydxm,dydx,b,am1;
  int i=1;

  //if(verbose)printf("n=%i, natural=%i\n",n,spline_natural);
  //First step initialization:
  dxm=xs[1]-xs[0];
  dx=xs[2]-xs[1];
  dydxm= (ys[1]-ys[0]) / dxm;
  dydx = (ys[2]-ys[1]) / dx;
  b= 2.0*(dx+dxm);
  if(spline_natural){//natural spline
    //these depend on the first line matrix elements
    v[1] = dx/b;
    s[1] = 6.0*( dydx - dydxm )/b;
  } else {//not-a-knot condition (effectively fourth deriv vanishing)
    v[1] = (dx-dxm) / (b-dxm);
    s[1] =  12.0*( dydx - dydxm )/b*dx/(b-dxm);
  }
  if(verbose) printf( "i,dx,dydx,v,s,znum %i %g %g %g %g\n" ,i,dx,dydx,v[i],s[i]);


  //forward loop
  for(i=2;i<n-2;i++){
    dxm=dx;
    dx=xs[i+1]-xs[i];
    dydxm=dydx;
    dydx = (ys[i+1]-ys[i]) / dx;
    b = 2.0*(dx+dxm) - dxm*v[i-1];
    v[i] = dx / b;
    s[i] = (6.0*( dydx - dydxm ) - s[i-1]*dxm)/ b ;
    if(verbose) printf( "i,dx,dydx,v,s,znum %i %g %g %g %g %g\n" ,i,dx,dydx,v[i],s[i],( dydx - dydxm )/dx);
  }

  //Reverse loop inititalization depends on BC:
  i=n-2;
  dxm=dx;
  dx=xs[i+1]-xs[i];
  dydxm=dydx;
  dydx = (ys[i+1]-ys[i]) / dx;
  if(spline_natural){//natural spline
    //these depend on the last line matrix elements
    //everything is identical to the bulk in this case
    am1=dxm;
    b = 2.0*(dx+dxm) - dx*v[i-1];
  } else { //not-a-knot
    am1=(1-dx/dxm)*(dxm+dx);
    b = (2.0+dx/dxm)*(dx+dxm) - am1*v[i-1];
  }
  v[i] = dx / b;
  s[i] =  (6.0*( dydx - dydxm ) - am1*s[i-1])/b;
  if(verbose) printf( "i,dx,dydx,v,s,znum %i %g %g %g %g %g\n",i,dx,dydx,v[i],s[i],( dydx - dydxm )/dx);


  //reverse loop
  zs[n-1]=0;  //this is temporary, to init the loop, not the final value
  for(i=n-2;i>0;i--){
    zs[i]= s[i] - v[i]*zs[i+1];
    if(verbose)printf( "i,z: %i %g\n",i,zs[i]);
  }

  if(spline_natural){
    zs[0]=0;
    zs[n-1]=0;
  } else {
    zs[0]=zs[1]-(xs[1]-xs[0])/(xs[2]-xs[1])*(zs[2]-zs[1]);
    zs[n-1]=zs[n-2]-(xs[n-2]-xs[n-1])/(xs[n-3]-xs[n-2])*(zs[n-3]-zs[n-2]);
  }
  if(verbose)printf( "i,z: 0 %g\n",zs[0]);
  if(1+verbose){
    for(i=1;i<n-1;i++){
      if(verbose)printf ("i,z, test: %i %g %g\n",i,zs[i],(xs[i]-xs[i-1])*zs[i-1]+2*(xs[i+1]-xs[i-1])*zs[i]+(xs[i+1]-xs[i])*zs[i+1]-6*((ys[i+1]-ys[i])/(xs[i+1]-xs[i])-(ys[i-1]-ys[i])/(xs[i-1]-xs[i])));
    }
  }
}


double spline_int(double xp, const double xs[],const double ys[],const double zs[],int n){
  double dx,t,ct,tct,ctmt,zL,zR,dA,Ax,yL,yR;
  int i;
  i=spline_findix(xp,xs,n);
  //printf("spline_int: %g < xp = %g < %g, i=%i\n",xs[0],xp,xs[n-1], i);
  if(i<0){
    printf("spline_int: xp is out of range.\n");
    exit(1);
  }
  dx = xs[i+1]-xs[i];
  t  = (xp-xs[i])/dx;
  if(t<0 || t>1){//do we allow extrapolation?
    printf("spline_int: Extrapolating? t=%7.4f, xp=%7.4f, xs[%i]=%7.4f,dx=%7.4f\n",t,xp,i,xs[i],dx);
    //int verbosetmp=verbose;verbose=1;
    //i=spline_findix(xp,xs,n);
    //printf("got i=%i\n",i);
    //verbose=verbosetmp;
  }
  ct = 1 - t ;
  tct = t*ct;
  ctmt=ct-t;
  zL=zs[i];
  zR=zs[i+1];
  dA = (zR-zL)/6.0;
  Ax = (zR+zL)/4.0 -ctmt*dA/2.0;
  yL=ys[i];
  yR=ys[i+1];
  //results:
  double q   = yL*ct + yR*t - tct*dx*dx*Ax;
  return q;
}


//Bisection search for interval xi array containing xp
//for internal use.
int spline_findix(double xp, const double xs[],int n){
  int done=false;
  int iL,iR,ic;
  double xc;
  //double xL,xR;
  //set bracketing limits
  iL=0;
  iR=n-1;
  ic=spline_iold;
  if(ic>n-2)ic=n-2;
  if(ic<1)ic=1;
  if(verbose){
    double xL=xs[iL];
    double xR=xs[iR];
    xc=xs[ic];
    printf( "Entering search for %g < %g < %g with x[%i < %i < %i ]= %g\n",xL,xp,xR,iL,ic,iR,xc);
    //for(int i=0;i<n;i++)printf(" xs[%i]=%7.4f\n",i,xs[i]);
    printf(" %7.4f<=xs[i]<=%7.4f\n",xs[0],xs[n-1]);
  }

  //first we try near the last location, hoping to get lucky
  if(xs[ic]<=xp){
    iL=ic;
    if(xs[ic+1]>=xp){
      iL=ic;
      done=true;
    } else if(ic+1<iR && xs[ic+2]>=xp){
      iL=ic+1;
      done=true;
    } else {
      iL=ic+2; //quick search failed
    }
  } else if(xs[ic-1]<=xp) { //try one more on the left
    iL=ic-1;
    done=true;
  } else {
    iR=ic-1;//quick search failed
  }
  if(!done){
    //bisection search
    //xR=xs[iR];
    //xL=xs[iL];
    while(iR-iL>1){
      ic=(iR+iL)/2;
      xc=xs[ic];
      if(xp>=xc){
	iL=ic;
	//xL=xc;
      }
      if(xp<xc){
	iR=ic;
	//xR=xc;
      }
    }
  }
  spline_iold=iL;
  return iL;
}

//The interpolating function on the [i,i+1) interval is then:
//  q(x) = y[i]*ct + y[i+1]*t - t*ct*dx^2*Ax
//with
//  dx = x[i+1]-x[i]
//   t = (x-x[i])/dx
//  ct = 1-t
//  zL = z[i]
//  zR = z[i+1]
//  dA = (zR-zL)/6
//  Ax = (zR+zL)/4 - (ct-t)*dA/2
//and the derivatives are:
//  q'(x) =  dy[i]/dx - dx*( (ct-t)*Ax + t*ct*dA)
// q''(x) =  2( Ax - (ct-t)*dA )
//q'''(x) =  6dA/dx =(zR-zL)/dx
// Arguments:
// xp    = evaluation point
// xs,ys = must be identical to xs that was provided to spline_construct
// zs    = output from spline_construct
// n     = array lengths
// q,dq1,dq2,dq3  = return pointers for the value and its derivatives.
void spline_intd3(double xp, const double xs[],const double ys[],const double zs[],int n,double *q, double *dq1, double *dq2, double *dq3){
  double dx,t,ct,tct,ctmt,zL,zR,dA,Ax,yL,yR;
  int i;
  i=spline_findix(xp,xs,n);
  //printf("spline_intd3: %g < xp = %g < %g, i=%i\n",xs[0],xp,xs[n-1], i);
  if(i<0){
    printf("spline_intd3: xp is out of range.\n");
    exit(1);
  }
  dx = xs[i+1]-xs[i];
  t  = (xp-xs[i])/dx;
  if(t<0 || t>1){//do we allow extrapolation?
    //printf("spline_intd3 t=%7.4f\n",t);
    printf("spline_intd3: Extrapolating? t=%7.4f, xp=%7.4f, xs[%i]=%7.4f,dx=%7.4f\n",t,xp,i,xs[i],dx);
  }
  ct = 1 - t ;
  tct = t*ct;
  ctmt=ct-t;
  zL=zs[i];
  zR=zs[i+1];
  dA = (zR-zL)/6.0;
  Ax = (zR+zL)/4.0 -ctmt*dA/2.0;
  yL=ys[i];
  yR=ys[i+1];
  //results:
  *q   = yL*ct + yR*t - tct*dx*dx*Ax;
  *dq1 = (yR-yL)/dx - dx * ( ctmt*Ax + tct*dA );
  *dq2 = 2 * ( Ax - ctmt*dA );
  *dq3 = (zR-zL)/dx;
  /*if(verbose){
    printf ("results for x: %g %g %g\n",xs[i],xp,xs[i+1]);
    printf("  q: %g %g %g\n", yL,*q,yR);
    printf("dq1: %g \n", (yR-yL)/dx);
    printf("qq2: %g %g %g\n", zL,*dq2,zR);
    }*/
}



//Interpolation of integrand on a new f-grid
///If min_f and max_f are nonpositive then the integrand is defined over the intersection of the f1 and f2 domains.
///If either of these is greater than zero, then that value is used for the repsective bound.
// void interpolate_ip_integrand(double *f1, double *s1Ar, double *s1Ai, double  *s1p, int n1, double *f2, double *s2Ar, double*s2Ai, double *s2p, int n2,  double (*Snoise)(double), double scalefactor, double *fnew, double *Ars,  double *Ais, double *dphis, int *nnew, double min_f, double max_f){
void interpolate_ip_integrand(double *f1, double *s1Ar, double *s1Ai, double  *s1p, int n1, double *f2, double *s2Ar, double*s2Ai, double *s2p, int n2,  ObjectFunction *Snoise, double scalefactor, double *fnew, double *Ars,  double *Ais, double *dphis, int *nnew, double min_f, double max_f){
  //#returns fnew, integrand
  //#f1,f2 should be arrays of ordered positive values
  //#s1,s2,Sn should be a function

  //Construct splines from the input data:
  double s1Arz[n1],s1Aiz[n1],s1pz[n1];
  double s2Arz[n2],s2Aiz[n2],s2pz[n2];
  spline_construct(f1,s1Ar,s1Arz,n1);
  spline_construct(f1,s1Ai,s1Aiz,n1);
  spline_construct(f1,s1p, s1pz, n1);
  spline_construct(f2,s2Ar,s2Arz,n2);
  spline_construct(f2,s2Ai,s2Aiz,n2);
  spline_construct(f2,s2p, s2pz, n2);

  //#Construct new grid based on spacing of the old grid spacings with continuity in df
  //#(could make this smoother by going higher order, but why).
  //#We define the new scaling to have appropriate limits, when either grids spacing is much smaller,
  //#or when they are equal.
  //#With scale=1, defining df12=df1+df2, mu=df1*df2/df12, want:
  //#  dfnew ~ min(df1,fd2) mu        when mu is small
  //#  dfnew ~ df12/2                 when mu is maximal, mu ~ df12/4
  //#  so
  //#  dfnew = mu*(1+4*mu/df12)

  double f_max=fmin(f1[n1-1],f2[n2-1]);
  if(max_f>0){
    if(f_max>=max_f)f_max=max_f;
    //SM: Removed this warning statement: waveforms have a natural higher frequency cutoff which is the last frequency generated by the ROM - max_f will typically be the 8kHz cutoff for Sn, so ok if f_max<max_f
    //else{
    //  printf("wip:interpolate_ip_integrand: Specified upper range bound max_f is beyond data range.\n");
    //}
  }
  double f_min=fmax(f1[0],f2[0]);
  if(min_f>0){
    if(f_min<=min_f)f_min=min_f;
    else{
      //SM: Here we keep the warning, as it indicates that we are really leaving out a physical part of the waveform (presumably because we called the ROM with too low of a total mass)
      printf("wip:interpolate_ip_integrand: Specified lower range bound min_f is beyond data range.\n");
      //TESTING
      printf("f1[0]: %.8e | f2[0]: %.8e | min_f: %.8e \n", f1[0], f2[0], min_f);
    }
  }
  double f=f_min;
  int NnewMax=*nnew;
  int i1=0,i2=0,inew=0;
  printf("Here 1 \n");
  while(1){//skip the first increment
    double df=0;
    if(inew>0){
      while(f1[i1]<f && i1<n1-2) i1+=1;//#must keep i<nf-1
      while(f2[i2]<f && i2<n2-2) i2+=1;
      double df1= f1[i1+1]-f1[i1];
      double df2= f2[i2+1]-f2[i2];
      double df12=df1+df2;
      double mu=df1*df2/df12;
      df=mu*(1+4.0*mu/df12)/scalefactor;
      // printf(" %i: df1=%7.4f ->df2=%7.4f, df = %7.4f\n",inew,df1,df2, df);
    }
    f+=df;
    if(f_max-f<wip_mindf)f=f_max;
    if(inew>=NnewMax){
      //printf(" inew==%i>=%i\n",inew,NnewMax);
      printf("wip:Error, out of space in array.\n");
      exit(1);
    }
    // printf("inew = %d, df = %9.6f, f = %9.6f \n", inew, df, f);
    fnew[inew]=f;
    //printf(" %i: df=%7.4f ->f=%7.4f\n",inew,df,f);
    //#now we compute the integrand on the new grid
    double A1r=spline_int(f,f1,s1Ar,s1Arz,n1);
    double A1i=spline_int(f,f1,s1Ai,s1Aiz,n1);
    // printf("A1: re = %9.6e, im = %9.6e  \n", A1r, A1i);
    double p1=spline_int(f,f1,s1p,s1pz,n1);
    double A2r=spline_int(f,f2,s2Ar,s2Arz,n2);
    double A2i=spline_int(f,f2,s2Ai,s2Aiz,n2);
    // printf("A2: re = %9.6e, im = %9.6e  \n", A2r, A2i);
    double p2=spline_int(f,f2,s2p,s2pz,n2);
    // printf("before Noise: ,  f = %9.6e, p2 = %9.6e  \n", f, p2);
    double Sn= ObjectFunctionCall(Snoise, f);
    // printf("Noise: ,  f = %9.6f, Sn = %9.6e  \n", f, Sn);
    //Ars[inew]=(A1r*A2r-A1i*A2i)/Sn;
    //Ais[inew]=(A1r*A2i+A1i*A2r)/Sn;
    Ars[inew]=(A1r*A2r+A1i*A2i)/Sn;//previously hadn't taken CC of S2
    Ais[inew]=(-A1r*A2i+A1i*A2r)/Sn;
    if(wip_uselogf){
      double ef=exp(f);
      Ars[inew]*=ef;
      Ais[inew]*=ef;
    }
    dphis[inew]=-p2+p1;
    // printf("As: re = %9.6e, im = %9.6e  dph = %9.6e \n", Ars[inew], Ais[inew], dphis[inew]);
    if(f>=f_max)break;
    inew++;
  }
  printf("Here 2 \n");
  fnew[inew]=f_max;
  *nnew=inew+1;
};


//#Now perform the integration...
double complex wip_phase (double *f1, int n1, double *f2, int n2, double *s1Ar, double *s1Ai, double  *s1p, double *s2Ar, double*s2Ai, double *s2p, ObjectFunction *Snoise, double scalefactor, double min_f, double max_f){
// double complex wip_phase (double *f1, int n1, double *f2, int n2, double *s1Ar, double *s1Ai, double  *s1p, double *s2Ar, double*s2Ai, double *s2p, double (*Snoise)(double), double scalefactor, double min_f, double max_f){
  int i;
  double f1x[n1],f2x[n2];
  if(wip_uselogf){
    for(i=0; i<n1;i++){
      double f=f1[i];
      if(f<=0)printf("wip_phase:Trouble using log grid: f1=%g\n",f);
      f1x[i]=log(f);
    }
    for(i=0; i<n2;i++){
      double f=f2[i];
      if(f<=0)printf("wip_phase:Trouble using log grid: f2=%g\n",f);
      f2x[i]=log(f);
    }
  } else {
    for(i=0; i<n1;i++)f1x[i]=f1[i];
    for(i=0; i<n2;i++)f2x[i]=f2[i];
  }

  //#integrate to new grid
  int NfMax=(n1+n2)*scalefactor;//should be adequately large.
  double fs[NfMax],Ars[NfMax],Ais[NfMax],Arz[NfMax],Aiz[NfMax],dphis[NfMax],dphiz[NfMax];//(does this work in C?)
  int nf=NfMax;
  printf("SB: let's try %9.5f, %9.6e \n", 1.e-3, ObjectFunctionCall(Snoise, 1.e-3) );
  printf("SB: interpolating....\n");
  interpolate_ip_integrand(f1x, s1Ar, s1Ai, s1p, n1, f2x, s2Ar, s2Ai, s2p, n2, Snoise, scalefactor, fs, Ars, Ais, dphis, &nf, min_f, max_f);
  printf("SB: interpolation is done\n");
  wip_count=nf;//just for testing/reference
  spline_construct(fs,Ars,Arz,nf);
  spline_construct(fs,Ais,Aiz,nf);
  spline_construct(fs,dphis,dphiz,nf);
  printf("SB: splines are constructed\n");

  //printf( "fs: [");for(i=0;i<nf-1;i++)printf(" %g,",fs[i]);printf(" %g ]\n",fs[nf-1]);
  //printf( "Ars: [");for(i=0;i<nf-1;i++)printf(" %g,",Ars[i]);printf(" %g ]\n",Ars[nf-1]);
  //printf( "Ais: [");for(i=0;i<nf-1;i++)printf(" %g,",Ais[i]);printf(" %g ]\n",Ais[nf-1]);
  //printf( "dphis: [");for(i=0;i<nf-1;i++)printf(" %g,",dphis[i]);printf(" %g ]\n",dphis[nf-1]);


  double complex intsum=0;
  //intd3vec(f,fs,Ars,Azs,nf)
  const double complex sqrti = cexp(0.25*I*PI);
  const double sqrtpi = sqrt(PI);

  for(i=0;i<nf-1;i++){//loop over freq intervals
    //#quadratic phase integration over the current interval (fs[0],fs[2]):
    //#might also try a linear-phase version for comparison...
    double f=fs[i];
    double eps=fs[i+1]-f;
    double eps2=eps*eps;
    //printf("[ %.8g -- %.8g ]\n",f,f+eps);

    //#amp coeffs
    double a0r,a0i,a1r,a1i,a2r,a2i,a3r,a3i;
    spline_intd3(f,fs,Ars,Arz,nf,&a0r,&a1r,&a2r,&a3r);
    spline_intd3(f,fs,Ais,Aiz,nf,&a0i,&a1i,&a2i,&a3i);
    double complex a0=a0r+I*a0i;
    double complex a1=a1r+I*a1i;
    double complex a2=a2r+I*a2i;
    double complex a3=a3r+I*a3i;

    a2=a2/2.0;
    a3=a3/6.0;
    printf("i=%i: As = %g+%gi ,  %g+%gi ,  %g+%gi , %g+%gi \n",i, a0r,a0i, a1r,a1i, a2r/2,a2i/2, a3r/6,a3i/6);
    double complex ampscale=1;
    if(wip_relative){
      ampscale=fabs(a0r)+fabs(a0i)+1.0e-50;
      a0=a0/ampscale;
      a1=a1/ampscale;
      a2=a2/ampscale;
      a3=a3/ampscale;
    }
    //#phase coefficients
    double p0,p1,p2,p3;
    spline_intd3(f,fs,dphis,dphiz,nf,&p0,&p1,&p2,&p3);
    p2=p2/2.0;
    p3=p3/6.0;

    //#Inm   = Integrate(x^n Em(x), {x,0,eps}) / E0
    //#Em(x) = Exp(p1*x+...pm*x^m)
    double p1e=p1*eps;
    double p2e2=p2*eps2;
    double p3e3=p3*eps*eps2;
    double complex E0=cexp(I*p0);
    double complex E2m1=cexpm1i(p1e+p2e2);
    double complex E2=E2m1+1.0;
    //#const phase terms
    double I00, I10, I20, I30;
    I00=eps;
    I10=eps2/2.0;
    I20=eps2*eps/3.0;
    I30=eps2*eps2/4.0;
    //#linear phase terms (if needed), and quadratic phase terms
    double complex I01, I11, I21, I31;
    double complex I02, I12, I22, I32;
    int useapproxquad=p2e2*p2e2*p2e2*p2<wip_errtol;
    printf("useapprquad = %d \n", useapproxquad);
    if(useapproxquad){
      double p1e2=p1e*p1e;
      if(p1e2*p1e2<wip_errtol){// #small p1e approx with errs order p1e^4
      	double complex ip1=I*p1;
      	I01=I00+ip1*(I10+ip1/2.0*(I20+ip1/3.0*I30));
      	I11=I10+ip1*(I20+ip1/2.0*I30);
      	I21=I20+ip1*I30;
      	I31=I30;
      } else {
      	double complex iop1=I/p1;
      	//A space for the cexpm1 function is reserved in the C99 standard, though it need not be implemented.
      	//If it exists it can avoid precision issues with small values of the argument.
      	double complex E1m1=cexpm1i(p1e);
      	double complex E1=1.0+E1m1;
      	I01=-iop1*E1m1;
      	I11=iop1*(I01-I00*E1);
      	I21=2.0*iop1*(I11-I10*E1);
      	I31=3.0*iop1*(I21-I20*E1);
      }
        //#small p2e approx with errs order p2e^2
        double complex ip2=I*p2;
        I02=I01+ip2*I21;
        I12=I11+ip2*I31;
        I22=I21;
        I32=I31;
    } else {
      double complex io2p2=0.5*I/p2;
      double complex s=csqrt(p2);
      double complex z0=p1/s/2.0;
      double complex w0=Faddeeva_w(sqrti*z0,wip_errtol);
      double complex wp=Faddeeva_w(sqrti*( eps*s+z0),wip_errtol);
      double complex ip1=I*p1;
      I02=-0.5/s*sqrti*sqrtpi*( E2*wp - w0 );
      I12= -io2p2*(E2m1-ip1*I02);
      I22=-io2p2*(I00*E2-I02-ip1*I12);
      I32=-io2p2*(2.0*(I10*E2-I12)-ip1*I22);
    }
    printf("i=%i: Ix2s = %g+%gi ,  %g+%gi ,  %g+%gi , %g+%gi \n",i,creal(I02),cimag(I02),creal(I12),cimag(I12),creal(I22),cimag(I22),creal(I32),cimag(I32));
    //#cubic phase terms
    //#small p3e3 approx is our only option (we can keep more terms... above)
    double complex I03, I13, I23, I33;
    double complex ip3=I*p3;
    I03=I02+ip3*I32;
    I13=I12;
    I23=I22;
    I33=I32;
    printf("i=%i: In3s = %g+%gi ,  %g+%gi ,  %g+%gi , %g+%gi \n",i,creal(I03),cimag(I03),creal(I13),cimag(I13),creal(I23),cimag(I23),creal(I33),cimag(I33));
    //#finish
    //printf("i=%i cofac=%g+%gi\n",i,creal(ampscale*E0),cimag(ampscale*E0));
    //printf("i=%i: ampscale=%g+%gi, E0=%g+%gi\n",i,creal(ampscale),cimag(ampscale),creal(E0),cimag(E0));
    //printf("i=%i: As = %g+%gi ,  %g+%gi ,  %g+%gi , %g+%gi \n",i,creal(a0),cimag(a0),creal(a1),cimag(a1),creal(a2),cimag(a2),creal(a3),cimag(a3));
    double complex dint=ampscale*E0*(a0*I03 + a1*I13 + a2*I23 + a3*I33);
    // #step
    intsum += dint;
    printf("i=%i: intsum,dint= %g+%gi, %g+%gi\n",i,creal(intsum),cimag(intsum),creal(dint),cimag(dint));
  }//#end of loop over freq intervals
  printf("SB: loop is done intsum = %f \n", intsum);

  return intsum;
};


///This is just the one-step part of the integration routine abstracted
double complex compute_int(double fleft, double fright, double *fs, int nf, double *Ars, double *Arz, double *Ais, double *Aiz, double *dphis, double *dphiz){
  double f=fleft;
  double eps=fright-fleft;
  double eps2=eps*eps;

  if(wip_adapt_verbose)printf("\n[ %.8g -- %.8g ] Delta=%.4g\n",fleft,fright,fright-fleft);
  //intd3vec(f,fs,Ars,Azs,nf)
  const double complex sqrti = cexp(0.25*I*PI);
  const double sqrtpi = sqrt(PI);

  //#amp coeffs
  double a0r,a0i,a1r,a1i,a2r,a2i,a3r,a3i;
  spline_intd3(f,fs,Ars,Arz,nf,&a0r,&a1r,&a2r,&a3r);
  spline_intd3(f,fs,Ais,Aiz,nf,&a0i,&a1i,&a2i,&a3i);
  double complex a0=a0r+I*a0i;
  double complex a1=a1r+I*a1i;
  double complex a2=a2r+I*a2i;
  double complex a3=a3r+I*a3i;

  a2=a2/2.0;
  a3=a3/6.0;
  if(wip_adapt_verbose)printf(" As = %g+%gi ,  %g+%gi ,  %g+%gi , %g+%gi \n",a0r,a0i,a1r,a1i,a2r/2,a2i/2,a3r/6,a3i/6);
  double complex ampscale=1;
  if(wip_relative){
    ampscale=fabs(a0r)+fabs(a0i)+1.0e-50;
    a0=a0/ampscale;
    a1=a1/ampscale;
    a2=a2/ampscale;
    a3=a3/ampscale;
  }
  //#phase coefficients
  double p0,p1,p2,p3;
  spline_intd3(f,fs,dphis,dphiz,nf,&p0,&p1,&p2,&p3);
  p2=p2/2.0;
  p3=p3/6.0;

  //#Inm   = Integrate(x^n Em(x), {x,0,eps}) / E0
  //#Em(x) = Exp(p1*x+...pm*x^m)
  double p1e=p1*eps;
  double p2e2=p2*eps2;
  double p3e3=p3*eps*eps2;
  double complex E0=cexp(I*p0);
  double complex E2m1=cexpm1i(p1e+p2e2);
  double complex E2=E2m1+1.0;
  //#const phase terms
  double I00, I10, I20, I30;
  I00=eps;
  I10=eps2/2.0;
  I20=eps2*eps/3.0;
  I30=eps2*eps2/4.0;
  if(wip_adapt_verbose)printf(" Ix0s = %g+%gi ,  %g+%gi ,  %g+%gi , %g+%gi \n",creal(I00),cimag(I00),creal(I10),cimag(I10),creal(I20),cimag(I20),creal(I30),cimag(I30));
  //#linear phase terms (if needed), and quadratic phase terms
  double complex I01, I11, I21, I31;
  double complex I02, I12, I22, I32;
  int useapproxquad=p2e2*p2e2*p2e2*p2<wip_errtol;
  if(useapproxquad){
    double p1e2=p1e*p1e;
    if(p1e2*p1e2<wip_errtol){// #small p1e approx with errs order p1e^4
      double complex ip1=I*p1;
      I01=I00+ip1*(I10+ip1/2.0*(I20+ip1/3.0*I30));
      I11=I10+ip1*(I20+ip1/2.0*I30);
      I21=I20+ip1*I30;
      I31=I30;
    } else {
      double complex iop1=I/p1;
      //A space for the cexpm1 function is reserved in the C99 standard, though it need not be implemented.
      //If it exists it can avoid precision issues with small values of the argument.
      double complex E1m1=cexpm1i(p1e);
      double complex E1=1.0+E1m1;
      I01=-iop1*E1m1;
      I11=iop1*(I01-I00*E1);
      I21=2.0*iop1*(I11-I10*E1);
      I31=3.0*iop1*(I21-I20*E1);
    }
    if(wip_adapt_verbose)printf(" Ix1s = %g+%gi ,  %g+%gi ,  %g+%gi , %g+%gi \n",creal(I02),cimag(I01),creal(I11),cimag(I11),creal(I21),cimag(I21),creal(I31),cimag(I31));
    //#small p2e approx with errs order p2e^2
    double complex ip2=I*p2;
    I02=I01+ip2*I21;
    I12=I11+ip2*I31;
    I22=I21;
    I32=I31;
  } else {
    double complex io2p2=0.5*I/p2;
    double complex s=csqrt(p2);
    double complex z0=p1/s/2.0;
    double complex w0=Faddeeva_w(sqrti*z0,wip_errtol);
    double complex wp=Faddeeva_w(sqrti*( eps*s+z0),wip_errtol);
    double complex ip1=I*p1;
    I02=-0.5/s*sqrti*sqrtpi*( E2*wp - w0 );
    I12= -io2p2*(E2m1-ip1*I02);
    I22=-io2p2*(I00*E2-I02-ip1*I12);
    I32=-io2p2*(2.0*(I10*E2-I12)-ip1*I22);
  }
  if(wip_adapt_verbose)printf(" Ix2s = %g+%gi ,  %g+%gi ,  %g+%gi , %g+%gi \n",creal(I02),cimag(I02),creal(I12),cimag(I12),creal(I22),cimag(I22),creal(I32),cimag(I32));
  //#cubic phase terms
  //#small p3e3 approx is our only option (we can keep more terms... above)
  double complex I03, I13, I23, I33;
  double complex ip3=I*p3;
  I03=I02+ip3*I32;
  I13=I12;
  I23=I22;
  I33=I32;
  if(wip_adapt_verbose){
    printf(" In3s = %g+%gi ,  %g+%gi ,  %g+%gi , %g+%gi \n",creal(I03),cimag(I03),creal(I13),cimag(I13),creal(I23),cimag(I23),creal(I33),cimag(I33));
  //#finish
    printf(" cofac=%g+%gi\n",creal(ampscale*E0),cimag(ampscale*E0));
    printf(" ampscale=%g+%gi, E0=%g+%gi\n",creal(ampscale),cimag(ampscale),creal(E0),cimag(E0));
    printf(" As = %g+%gi ,  %g+%gi ,  %g+%gi , %g+%gi \n",creal(a0),cimag(a0),creal(a1),cimag(a1),creal(a2),cimag(a2),creal(a3),cimag(a3));
  }
  double complex dint=ampscale*E0*(a0*I03 + a1*I13 + a2*I23 + a3*I33);
  // #step
  return dint;
};




//#Now perform the integration adaptively.
//This version is the same as above, but implementing an experimental adaptive approach.
//In this approach we test and refine the grid as needed to reach a target accuracy
double complex wip_adaptive_phase (double *f1, int n1, double *f2, int n2, double *s1Ar, double *s1Ai, double  *s1p, double *s2Ar, double*s2Ai, double *s2p, double (*Snoise)(double), double scalefactor, int downsample, double errtol, double min_f, double max_f){
  const int maxlev=20;
  int i;
  double f1x[n1],f2x[n2];

  if(wip_uselogf){
    for(i=0; i<n1;i++){
      double f=f1[i];
      if(f<=0)printf("wip_adaptive_phase:Trouble using log grid: f1=%g\n",f);
      f1x[i]=log(f);
    }
    for(i=0; i<n2;i++){
      double f=f2[i];
      if(f<=0)printf("wip_adaptive_phase:Trouble using log grid: f2=%g\n",f);
      f2x[i]=log(f);
    }
  } else {
    for(i=0; i<n1;i++)f1x[i]=f1[i];
    for(i=0; i<n2;i++)f2x[i]=f2[i];
  }

  //#integrate to new grid
  int NfMax=(n1+n2)*scalefactor;//should be adequately large.
  double fs[NfMax],Ars[NfMax],Ais[NfMax],Arz[NfMax],Aiz[NfMax],dphis[NfMax],dphiz[NfMax];//(does this work in C?)
  int nf=NfMax;
  interpolate_ip_integrand(f1x, s1Ar, s1Ai, s1p, n1, f2x, s2Ar, s2Ai, s2p, n2, Snoise, scalefactor, fs, Ars, Ais, dphis, &nf, min_f, max_f);
  spline_construct(fs,Ars,Arz,nf);
  spline_construct(fs,Ais,Aiz,nf);
  spline_construct(fs,dphis,dphiz,nf);

  //printf( "fs: [");for(i=0;i<nf-1;i++)printf(" %g,",fs[i]);printf(" %g ]\n",fs[nf-1]);
  //printf( "Ars: [");for(i=0;i<nf-1;i++)printf(" %g,",Ars[i]);printf(" %g ]\n",Ars[nf-1]);
  //printf( "Ais: [");for(i=0;i<nf-1;i++)printf(" %g,",Ais[i]);printf(" %g ]\n",Ais[nf-1]);
  //printf( "dphis: [");for(i=0;i<nf-1;i++)printf(" %g,",dphis[i]);printf(" %g ]\n",dphis[nf-1]);

  int ncount=0;
  double complex intsum=0;
  double dfs[NfMax];
  int ndf=(nf-2)/downsample+2;
  //perform downsampling
  for(i=0;i<nf-1;i+=downsample)dfs[i/downsample]=fs[i];
  dfs[ndf-1]=fs[nf-1];
  //printf("fs = {%.4g,%.4g,%.4g,...,%.4g,%.4g,%.4g}\n",fs[0],fs[1],fs[2],fs[nf-3],fs[nf-2],fs[nf-1]);
  //printf("dfs = {%.4g,%.4g,%.4g,...,%.4g,%.4g,%.4g}\n",dfs[0],dfs[1],dfs[2],dfs[ndf-3],dfs[ndf-2],dfs[ndf-1]);
  double basetol=errtol/ndf;
  //printf("ndf=%i, wip_adaptTOL=%0.8g, base_tol=%0.8g,\n",ndf,wip_adaptTOL,basetol);
  for(i=0;i<ndf-1;i+=1){//loop over freq intervals
    int ilevel=0;
    double fleft,fright[maxlev+1];
    double complex right_result[maxlev+1];
    fleft=dfs[i];fright[0]=dfs[i+1];
    int on_right=0;
    int have_trial=0;
    double complex trial=0;
    while(ilevel>=0){
      if(!have_trial){//First time through at this level, we compute the left-half integral.
	if(wip_adapt_verbose)printf("\nTrial:\n");
	trial=compute_int(fleft,fright[ilevel],fs, nf, Ars, Arz, Ais, Aiz, dphis, dphiz);
	have_trial=1;
	ncount++;
	if(wip_adapt_verbose)printf("computed: lev=%2i, trial=%.8g+i%.8g\n",ilevel, creal(trial), cimag(trial));
      }
      int pass=0;
      double complex oldtrial=trial;
      double complex left_result;
      if(wip_adapt_verbose)printf("set: lev=%2i, oldtrial=%.8g+i%.8g\n",ilevel, creal(oldtrial), cimag(oldtrial));
      double fcent=(fleft+fright[ilevel])/2.0;

      if(ilevel>=maxlev){
	printf("wip_adaptive_phase:Reached maxlevel without passing tolerance test.\n");
	pass=1;
      } else {
	int k;
	if(wip_adapt_verbose)for(k=0;k<=ilevel;k++)printf("%i: [%.8g,%.8g]\n",k,fleft,fright[k]);

	if(wip_adapt_verbose)printf("\n  left:\n");
	left_result=compute_int(fleft,fcent,fs, nf, Ars, Arz, Ais, Aiz, dphis, dphiz);
	if(wip_adapt_verbose)printf("\n  right:\n");
	right_result[ilevel]=compute_int(fcent,fright[ilevel],fs, nf, Ars, Arz, Ais, Aiz, dphis, dphiz);
	trial=left_result+right_result[ilevel];
	double complex diff=trial-oldtrial;
	ncount+=2;
	double tol=basetol/(2<<ilevel);
	pass=cabs(diff)<tol;
	if(wip_adapt_verbose){
	  printf("check: lev=%2i, oldtrial=%.8g+i%.8g\n",ilevel, creal(oldtrial), cimag(oldtrial));
	  printf("testing: lev=%2i, fleft=%.8g, fright=%.8g\n",ilevel,fleft,fright[ilevel]);
	  printf("testing: lev=%2i, left_result=%.8g+i%.8g, right_result=%.8g+i%.8g\n",ilevel,creal(left_result),cimag(left_result),creal(right_result[ilevel]),cimag(right_result[ilevel]));
	  printf("testing: lev=%2i, trial=%.8g+i%.8g,  oldtrial=%.8g+i%.8g\n",ilevel,creal(trial),cimag(trial),creal(oldtrial),cimag(oldtrial));
	  printf("testing: lev=%2i, 2^lev=%i, diff=%.8g, tol=%.8g\n",ilevel,2<<ilevel,cabs(diff),tol);
	}
      }
      if(!pass){
	//We fail the test and should continue refining
	fright[ilevel+1]=fcent;
	trial=left_result;
	ilevel++;
	if(wip_adapt_verbose)printf("failing: lev=%2i, fleft=%.8g, fright=%.8g\n",ilevel,fleft,fright[ilevel]);
      } else {
	//We pass the test and thus we finish with this level, and advance the marker.
	intsum+=trial;
	fleft=fright[ilevel];
	fright[ilevel]=fright[ilevel-1];
	trial=right_result[ilevel-1];
	while(ilevel>0&&fright[ilevel-1]<=fleft)ilevel--;
	if(ilevel>0){
	  trial=right_result[ilevel-1];
	  fright[ilevel]=fright[ilevel-1];
	} else {
	  ilevel=-1;
	  have_trial=0;
	}
	if(wip_adapt_verbose)printf("passing:   -->lev=%2i, fleft=%.8g, fright=%.8g\n",ilevel,fleft,fright[ilevel]);
      }
    }
  }
  wip_count=ncount;//just for testing/reference
  return intsum;
}

// Noise

const double ConstL = 2.5e9;

/* Proof mass and optical noises - f in Hz */
/* LISA SRD copied from the LISA Data Challenge pipeline */
static double SpmLISAProposal(const double f) {
  /* Acceleration noise */
  double noise_Sa_a = 9.e-30; /* m^2/sec^4 /Hz */
  /* In acceleration */
  double Sa_a = noise_Sa_a * (1.0 + pow(0.4e-3/f, 2)) * (1.0 + pow((f/8e-3), 4));
  /* In displacement */
  double Sa_d = Sa_a * pow(2.*PI*f, -4);
  /* In relative frequency unit */
  double Sa_nu = Sa_d * pow(2.*PI*f/C_SI, 2);
  double Spm = Sa_nu;
  return Spm;
}

static double SopLISAProposal(const double f) {
  /* Optical Metrology System noise */
  double noise_Soms_d = pow((10e-12), 2); /* m^2/Hz This is proposal value  */
  //double noise_Soms_d = pow((15e-12), 2); /* m^2/Hz This is science requirement value */
  /* In displacement */
  double Soms_d = noise_Soms_d * (1. + pow(2.e-3/f, 4));
  /* In relative frequency unit */
  double Soms_nu = Soms_d * pow(2.*PI*f/C_SI, 2);
  double Sop = Soms_nu;
  return Sop;
}


/* NOTE that Sylvain's notations are factor 4 smaller than LDC for A, E, T */
/* Noise functions for AET(XYZ) without rescaling */
/* Scaling by 2*sin2pifL^2 put back */
double SnAProposal(double f) {
  double twopifL = 2.*PI*ConstL/C_SI*f;
  double c2 = cos(twopifL);
  double c4 = cos(2*twopifL);
  double s2 = sin(twopifL);
  double Spm = SpmLISAProposal(f);
  double Sop = SopLISAProposal(f);
  return 8.*s2*s2 * (2.*(3. + 2.*c2 + c4)*Spm + (2. + c2)*Sop);
}
/* Scaling by 2*sin2pifL^2 put back */
double SnEProposal(double f) {
  double twopifL = 2.*PI*ConstL/C_SI*f;
  double c2 = cos(twopifL);
  double c4 = cos(2*twopifL);
  double s2 = sin(twopifL);
  double Spm = SpmLISAProposal(f);
  double Sop = SopLISAProposal(f);
  return 8.*s2*s2 * (2*(3. + 2*c2 + c4)*Spm + (2 + c2)*Sop);
}
/* Scaling by 8*sin2pifL^2*sinpifL^2 put back*/
double SnTProposal(double f) {
  double pifL = PI*ConstL/C_SI*f;
  double s1 = sin(pifL);
  double s2 = sin(2*pifL);
  double Spm = SpmLISAProposal(f);
  double Sop = SopLISAProposal(f);
  return 32.*s1*s1*s2*s2 * (4*s1*s1*Spm + Sop);
}

/* Rescaled by 4*sin2pifL^2 */
double SnXProposal(double f) {
  double twopifL = 2.*PI*ConstL/C_SI*f;
  double c2 = cos(twopifL);
  double s2 = sin(twopifL);
  double Spm = SpmLISAProposal(f);
  double Sop = SopLISAProposal(f);
  return 16.*s2*s2*( 2.*(1. + c2*c2)*Spm + Sop );
}
