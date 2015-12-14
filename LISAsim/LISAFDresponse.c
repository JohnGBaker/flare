/**
 * \author Sylvain Marsat, University of Maryland - NASA GSFC
 *
 * \brief C code for the implementation of the Fourier domain response for LISA-like detectors.
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
#include "EOBNRv2HMROMstruct.h"
#include "LISAgeometry.h"
#include "LISAFDresponse.h"


/***************************************/
/********* Core functions **************/

/* Core function processing a signal (in the form of a list of modes) through the Fourier-domain LISA response, for given values of the inclination, position in the sky and polarization angle - here simplified version for just the y_21 observable */
int LISASimFDResponse21(
  struct tagListmodesCAmpPhaseFrequencySeries **list,  /* Input/Output: list of modes in Frequency-domain amplitude and phase form as produced by the ROM, and output after FD response processing */
  const double inclination,                                   /* Inclination of the source */
  const double lambda,                                        /* First angle for the position in the sky */
  const double beta,                                          /* Second angle for the position in the sky */
  const double psi)                                           /* Polarization angle */
{
  /* Computing the complicated trigonometric coefficients */
  //clock_t begsetcoeffs = clock();
  SetCoeffsG(lambda, beta, psi);
  //clock_t endsetcoeffs = clock();
  //printf("Set Coeffs time: %g s\n", (double)(endsetcoeffs - begsetcoeffs) / CLOCKS_PER_SEC);

  /* Main loop over the modes - goes through all the modes present */
  ListmodesCAmpPhaseFrequencySeries* listelement = *list;
  while(listelement) {

    /* Definitions: l,m, frequency series and length */
    int l = listelement->l;
    int m = listelement->m;
    CAmpPhaseFrequencySeries* freqseries = listelement->freqseries;
    gsl_vector* freq = freqseries->freq;
    gsl_vector* amp_real = freqseries->amp_real;
    gsl_vector* amp_imag = freqseries->amp_imag;
    gsl_vector* phase = freqseries->phase;
    int len = (int) freq->size;
    double f, tf, bphi;
    double complex camp;

    /* Computing the Ylm combined factors for plus and cross for this mode */
    /* Capital Phi is set to 0 by convention */
    double complex Yfactorplus;
    double complex Yfactorcross;
    if (!(l%2)) {
      Yfactorplus = 1/2 * (SpinWeightedSphericalHarmonic(inclination, 0., -2, l, m) + conj(SpinWeightedSphericalHarmonic(inclination, 0., -2, l, -m)));
      Yfactorcross = I/2 * (SpinWeightedSphericalHarmonic(inclination, 0., -2, l, m) - conj(SpinWeightedSphericalHarmonic(inclination, 0., -2, l, -m)));
    }
    else {
      Yfactorplus = 1/2 * (SpinWeightedSphericalHarmonic(inclination, 0., -2, l, m) - conj(SpinWeightedSphericalHarmonic(inclination, 0., -2, l, -m)));
      Yfactorcross = I/2 * (SpinWeightedSphericalHarmonic(inclination, 0., -2, l, m) + conj(SpinWeightedSphericalHarmonic(inclination, 0., -2, l, -m)));
    }

    /* First step of the processing: orbital delay */
    /* Initializing spline for the phase */
    gsl_spline* spline_phi = gsl_spline_alloc(gsl_interp_cspline, len);
    gsl_interp_accel* accel_phi = gsl_interp_accel_alloc();
    gsl_spline_init(spline_phi, gsl_vector_const_ptr(freq, 0), gsl_vector_const_ptr(phase, 0), len);
    /* Vector keeping track of the Bessel Phase correction - to be used in the next step */
    gsl_vector* besselphi = gsl_vector_alloc(len);

    /* Loop over the frequencies - computing the correction due to the orbital delay */
    for(int j=0; j<len; j++) {
      f = gsl_vector_get(freq, j);
      tf =  gsl_spline_eval_deriv(spline_phi, f, accel_phi)/(2*PI);
      bphi = -2*PI*f*R_SI/C_SI*cos(beta) * cos( Omega_SI*tf - lambda );
      camp = gsl_vector_get(amp_real, j) * cexp(I*bphi); /* Amplitude is real before applying this first delay */
      gsl_vector_set(amp_real, j, creal(camp));
      gsl_vector_set(amp_imag, j, cimag(camp));
      gsl_vector_set(besselphi, j, bphi);
    }
    
    /* Second step of the processing: constellation delay/modulation */
    /* Initializing spline for the bessel phase */
    gsl_spline* spline_besselphi = gsl_spline_alloc(gsl_interp_cspline, len);
    gsl_interp_accel* accel_besselphi = gsl_interp_accel_alloc();
    gsl_spline_init(spline_besselphi, gsl_vector_const_ptr(freq, 0), gsl_vector_const_ptr(besselphi, 0), len);
    /* Loop over the frequencies - computing the correction due to the orbital delay */
    for(int j=0; j<len; j++) {
      f = gsl_vector_get(freq, j);
      tf =  (gsl_spline_eval_deriv(spline_phi, f, accel_phi) + gsl_spline_eval_deriv(spline_besselphi, f, accel_besselphi))/(2*PI);
      camp = G21mode(f, tf, Yfactorplus, Yfactorcross) * (gsl_vector_get(amp_real, j) + I * gsl_vector_get(amp_imag, j));
      /**/
      gsl_vector_set(amp_real, j, creal(camp));
      gsl_vector_set(amp_imag, j, cimag(camp));
    } 

    listelement = listelement->next;

    /* Clean up */
    gsl_spline_free(spline_phi);
    gsl_interp_accel_free(accel_phi);
    gsl_vector_free(besselphi);
    gsl_spline_free(spline_besselphi);
    gsl_interp_accel_free(accel_besselphi);
  }

  return SUCCESS;
}

//WARNING: tRef is ignored for now in the response - i.e. set to 0
/* Core function processing a signal (in the form of a list of modes) through the Fourier-domain LISA response, for given values of the inclination, position in the sky and polarization angle */
int LISASimFDResponseTDI3Chan(
  struct tagListmodesCAmpPhaseFrequencySeries **list,      /* Input: list of modes in Frequency-domain amplitude and phase form as produced by the ROM */
  struct tagListmodesCAmpPhaseFrequencySeries **listTDI1,  /* Output: list of contribution of each mode in Frequency-domain amplitude and phase form, in the TDI channel 1 */
  struct tagListmodesCAmpPhaseFrequencySeries **listTDI2,  /* Output: list of contribution of each mode in Frequency-domain amplitude and phase form, in the TDI channel 2 */
  struct tagListmodesCAmpPhaseFrequencySeries **listTDI3,  /* Output: list of contribution of each mode in Frequency-domain amplitude and phase form, in the TDI channel 3 */
  const double tRef,                                          /* Time at coalescence */
  const double lambda,                                        /* First angle for the position in the sky */
  const double beta,                                          /* Second angle for the position in the sky */
  const double inclination,                                   /* Inclination of the source */
  const double psi,                                           /* Polarization angle */
  const TDItag tditag)                                        /* Selector for the set of TDI observables */
{
  /* Computing the complicated trigonometric coefficients */
  //clock_t begsetcoeffs = clock();
  SetCoeffsG(lambda, beta, psi);
  //clock_t endsetcoeffs = clock();
  //printf("Set Coeffs time: %g s\n", (double)(endsetcoeffs - begsetcoeffs) / CLOCKS_PER_SEC);

  /* Main loop over the modes - goes through all the modes present, stopping when encountering NULL */
  ListmodesCAmpPhaseFrequencySeries* listelement = *list;
  while(listelement) {

    /* Definitions: l,m, frequency series and length */
    int l = listelement->l;
    int m = listelement->m;
    CAmpPhaseFrequencySeries* freqseries = listelement->freqseries;
    gsl_vector* freq = freqseries->freq;
    gsl_vector* amp_real = freqseries->amp_real;
    gsl_vector* amp_imag = freqseries->amp_imag;
    gsl_vector* phase = freqseries->phase;
    int len = (int) freq->size;
    double f, tf, bphi;
    double complex g21mode = 0.;
    double complex g12mode = 0.;
    double complex g32mode = 0.;
    double complex g23mode = 0.;
    double complex g13mode = 0.;
    double complex g31mode = 0.;
    double complex camp;
    double complex camp1;
    double complex camp2;
    double complex camp3;
    double complex factor1 = 0.;
    double complex factor2 = 0.;
    double complex factor3 = 0.;

    /* Computing the Ylm combined factors for plus and cross for this mode */
    /* Capital Phi is set to 0 by convention */
    double complex Yfactorplus;
    double complex Yfactorcross;
    if (!(l%2)) {
      //
      //printf("l even: %d\n", l%2);
      Yfactorplus = 1./2 * (SpinWeightedSphericalHarmonic(inclination, 0., -2, l, m) + conj(SpinWeightedSphericalHarmonic(inclination, 0., -2, l, -m)));
      Yfactorcross = I/2 * (SpinWeightedSphericalHarmonic(inclination, 0., -2, l, m) - conj(SpinWeightedSphericalHarmonic(inclination, 0., -2, l, -m)));
    }
    else {
      //
      //printf("l odd\n");
      Yfactorplus = 1./2 * (SpinWeightedSphericalHarmonic(inclination, 0., -2, l, m) - conj(SpinWeightedSphericalHarmonic(inclination, 0., -2, l, -m)));
      Yfactorcross = I/2 * (SpinWeightedSphericalHarmonic(inclination, 0., -2, l, m) + conj(SpinWeightedSphericalHarmonic(inclination, 0., -2, l, -m)));
    }
    //
/*     printf("Ylm: %g + I*%g\n", creal(SpinWeightedSphericalHarmonic(inclination, 0., -2, l, m)), cimag(SpinWeightedSphericalHarmonic(inclination, 0., -2, l, m))); */
/* printf("Ylminusm: %g + I*%g\n", creal(SpinWeightedSphericalHarmonic(inclination, 0., -2, l, -m)), cimag(SpinWeightedSphericalHarmonic(inclination, 0., -2, l, -m))); */
/*     printf("Yfactorplus : %g + I*%g\n", creal(Yfactorplus), cimag(Yfactorplus)); */
/*     printf("Yfactorcross : %g + I*%g\n", creal(Yfactorcross), cimag(Yfactorcross)); */

    /* First step of the processing: orbital delay */
    /* Initializing spline for the phase */
    gsl_spline* spline_phi = gsl_spline_alloc(gsl_interp_cspline, len);
    gsl_interp_accel* accel_phi = gsl_interp_accel_alloc();
    gsl_spline_init(spline_phi, gsl_vector_const_ptr(freq, 0), gsl_vector_const_ptr(phase, 0), len);
    /* Intermediate vectors for the complex amplitude of h0 and for the  Bessel Phase correction */
    gsl_vector* h0tilde_ampreal = gsl_vector_alloc(len);
    gsl_vector* h0tilde_ampimag = gsl_vector_alloc(len);
    gsl_vector* besselphi = gsl_vector_alloc(len);

    /* Loop over the frequencies - computing the correction due to the orbital delay */
    //clock_t tbegorbital = clock();
    for(int j=0; j<len; j++) {
      f = gsl_vector_get(freq, j);
      tf = gsl_spline_eval_deriv(spline_phi, f, accel_phi)/(2*PI);
      bphi = -2*PI*f*R_SI/C_SI*cos(beta) * cos( Omega_SI*tf - lambda );
      camp = gsl_vector_get(amp_real, j) * cexp(I*bphi); /* Amplitude is real before applying this first delay */
      gsl_vector_set(h0tilde_ampreal, j, creal(camp));
      gsl_vector_set(h0tilde_ampimag, j, cimag(camp));
      gsl_vector_set(besselphi, j, bphi);
    }
    //clock_t tendorbital = clock();
    //printf("Set orbital delay time: %g s\n", (double)(tendorbital - tbegorbital) / CLOCKS_PER_SEC);
    
    /* Second step of the processing: constellation delay/modulation */
    /* Initializing frequency series structure for this mode, for each of the TDI observables */
    CAmpPhaseFrequencySeries *modefreqseries1 = NULL;
    CAmpPhaseFrequencySeries *modefreqseries2 = NULL;
    CAmpPhaseFrequencySeries *modefreqseries3 = NULL;
    CAmpPhaseFrequencySeries_Init(&modefreqseries1, len);
    CAmpPhaseFrequencySeries_Init(&modefreqseries2, len);
    CAmpPhaseFrequencySeries_Init(&modefreqseries3, len);
    gsl_vector* freq1 = modefreqseries1->freq;
    gsl_vector* amp_real1 = modefreqseries1->amp_real;
    gsl_vector* amp_imag1 = modefreqseries1->amp_imag;
    gsl_vector* phase1 = modefreqseries1->phase;
    gsl_vector* freq2 = modefreqseries2->freq;
    gsl_vector* amp_real2 = modefreqseries2->amp_real;
    gsl_vector* amp_imag2 = modefreqseries2->amp_imag;
    gsl_vector* phase2 = modefreqseries2->phase;
    gsl_vector* freq3 = modefreqseries3->freq;
    gsl_vector* amp_real3 = modefreqseries3->amp_real;
    gsl_vector* amp_imag3 = modefreqseries3->amp_imag;
    gsl_vector* phase3 = modefreqseries3->phase;
    /* Initializing spline for the bessel phase */
    gsl_spline* spline_besselphi = gsl_spline_alloc(gsl_interp_cspline, len);
    gsl_interp_accel* accel_besselphi = gsl_interp_accel_alloc();
    gsl_spline_init(spline_besselphi, gsl_vector_const_ptr(freq, 0), gsl_vector_const_ptr(besselphi, 0), len);

    /* Loop over the frequencies - computing the correction due to the constellation delay/modulation */
    //clock_t tbegcontesllation = clock();
    //double timingcumulativeGABmode = 0;
    for(int j=0; j<len; j++) {
      f = gsl_vector_get(freq, j);
      tf = (gsl_spline_eval_deriv(spline_phi, f, accel_phi) + gsl_spline_eval_deriv(spline_besselphi, f, accel_besselphi))/(2*PI) ;
      //clock_t tbegGAB = clock();
      EvaluateGABmode(&g12mode, &g21mode, &g23mode, &g32mode, &g31mode, &g13mode, f, tf, Yfactorplus, Yfactorcross);
      //clock_t tendGAB = clock();
      //timingcumulativeGABmode += (double) (tendGAB-tbegGAB) /CLOCKS_PER_SEC;
      /**/
      EvaluateTDIfactor3Chan(&factor1, &factor2, &factor3, g12mode, g21mode, g23mode, g32mode, g31mode, g13mode, f, tditag);
      double complex amphOtilde = gsl_vector_get(h0tilde_ampreal, j) + I * gsl_vector_get(h0tilde_ampimag, j);
      camp1 = factor1 * amphOtilde;
      camp2 = factor2 * amphOtilde;
      camp3 = factor3 * amphOtilde;
      /**/
      gsl_vector_set(amp_real1, j, creal(camp1));
      gsl_vector_set(amp_imag1, j, cimag(camp1));
      gsl_vector_set(amp_real2, j, creal(camp2));
      gsl_vector_set(amp_imag2, j, cimag(camp2));
      gsl_vector_set(amp_real3, j, creal(camp3));
      gsl_vector_set(amp_imag3, j, cimag(camp3));
    }
    //clock_t tendcontesllation = clock();
  //printf("Set constellation time: %g s\n", (double)(tendcontesllation - tbegcontesllation) / CLOCKS_PER_SEC);
  //printf("GAB cumulated time: %g s\n", timingcumulativeGABmode);

    /* Copying the vectors of frequencies and phases */
    gsl_vector_memcpy(freq1, freq);
    gsl_vector_memcpy(freq2, freq);
    gsl_vector_memcpy(freq3, freq);
    gsl_vector_memcpy(phase1, phase);
    gsl_vector_memcpy(phase2, phase);
    gsl_vector_memcpy(phase3, phase);
    
    /* Append the modes to the ouput list-of-modes structures */
    *listTDI1 = ListmodesCAmpPhaseFrequencySeries_AddModeNoCopy(*listTDI1, modefreqseries1, l, m);
    *listTDI2 = ListmodesCAmpPhaseFrequencySeries_AddModeNoCopy(*listTDI2, modefreqseries2, l, m);
    *listTDI3 = ListmodesCAmpPhaseFrequencySeries_AddModeNoCopy(*listTDI3, modefreqseries3, l, m);

    /* Going to the next mode in the list */
    listelement = listelement->next;

    /* Clean up */
    gsl_spline_free(spline_phi);
    gsl_interp_accel_free(accel_phi);
    gsl_vector_free(h0tilde_ampreal);
    gsl_vector_free(h0tilde_ampimag);
    gsl_vector_free(besselphi);
    gsl_spline_free(spline_besselphi);
    gsl_interp_accel_free(accel_besselphi);
  }

  return SUCCESS;
}
