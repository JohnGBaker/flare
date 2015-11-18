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
int LISASimFDResponseTDIAET(
  struct tagListmodesCAmpPhaseFrequencySeries **list,  /* Input: list of modes in Frequency-domain amplitude and phase form as produced by the ROM */
  struct tagListmodesCAmpPhaseFrequencySeries **listA,  /* Output: list of contribution of each mode in Frequency-domain amplitude and phase form, in the second-generation TDI observable A */
  struct tagListmodesCAmpPhaseFrequencySeries **listE,  /* Output: list of contribution of each mode in Frequency-domain amplitude and phase form, in the second-generation TDI observable E */
  struct tagListmodesCAmpPhaseFrequencySeries **listT,  /* Output: list of contribution of each mode in Frequency-domain amplitude and phase form, in the second-generation TDI observable T */
  const double tRef,                                          /* Time at coalescence */
  const double lambda,                                        /* First angle for the position in the sky */
  const double beta,                                          /* Second angle for the position in the sky */
  const double inclination,                                   /* Inclination of the source */
  const double psi)                                           /* Polarization angle */
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
    double complex g21mode = 0;
    double complex g12mode = 0;
    double complex g32mode = 0;
    double complex g23mode = 0;
    double complex g13mode = 0;
    double complex g31mode = 0;
    double complex camp;
    double complex campA;
    double complex campE;
    double complex campT;

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
    CAmpPhaseFrequencySeries *modefreqseriesA = NULL;
    CAmpPhaseFrequencySeries *modefreqseriesE = NULL;
    CAmpPhaseFrequencySeries *modefreqseriesT = NULL;
    CAmpPhaseFrequencySeries_Init(&modefreqseriesA, len);
    CAmpPhaseFrequencySeries_Init(&modefreqseriesE, len);
    CAmpPhaseFrequencySeries_Init(&modefreqseriesT, len);
    gsl_vector* freqA = modefreqseriesA->freq;
    gsl_vector* amp_realA = modefreqseriesA->amp_real;
    gsl_vector* amp_imagA = modefreqseriesA->amp_imag;
    gsl_vector* phaseA = modefreqseriesA->phase;
    gsl_vector* freqE = modefreqseriesE->freq;
    gsl_vector* amp_realE = modefreqseriesE->amp_real;
    gsl_vector* amp_imagE = modefreqseriesE->amp_imag;
    gsl_vector* phaseE = modefreqseriesE->phase;
    gsl_vector* freqT = modefreqseriesT->freq;
    gsl_vector* amp_realT = modefreqseriesT->amp_real;
    gsl_vector* amp_imagT = modefreqseriesT->amp_imag;
    gsl_vector* phaseT = modefreqseriesT->phase;
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
      /* g21mode = G21mode(f, tf, Yfactorplus, Yfactorcross); */
      /* g12mode = G12mode(f, tf, Yfactorplus, Yfactorcross); */
      /* g32mode = G32mode(f, tf, Yfactorplus, Yfactorcross); */
      /* g23mode = G23mode(f, tf, Yfactorplus, Yfactorcross); */
      /* g13mode = G13mode(f, tf, Yfactorplus, Yfactorcross); */
      /* g31mode = G31mode(f, tf, Yfactorplus, Yfactorcross); */
      EvaluateGABmode(&g12mode, &g21mode, &g23mode, &g32mode, &g31mode, &g13mode, f, tf, Yfactorplus, Yfactorcross);
      //clock_t tendGAB = clock();
      //timingcumulativeGABmode += (double) (tendGAB-tbegGAB) /CLOCKS_PER_SEC;
      /**/
      //Previous version (presumably wrong, to be checked more)
      /* double sin3pifL = sin(3*PI*f*L_SI/C_SI); */
      /* double complex exp3ipifL = cexp(3*I*PI*f*L_SI/C_SI); */
      /* double complex exp2ipifL = cexp(2*I*PI*f*L_SI/C_SI); */
      /* double complex exp4ipifL = cexp(4*I*PI*f*L_SI/C_SI); */
      /* double complex commonfac = -2*I*sin3pifL*exp3ipifL; */
      /* double complex amphOtilde = (gsl_vector_get(amp_real, j) + I * gsl_vector_get(amp_imag, j)); */
      /**/
      /* campA = commonfac/sqrt(2) * ( (g13mode+g31mode)*(exp4ipifL-1.) + (g21mode+g23mode)*(1.-exp2ipifL) + (g32mode+g12mode)*(exp2ipifL-exp4ipifL) ) * amphOtilde; */
      /* campE = commonfac/sqrt(6) * ( (g23mode-g21mode)*(1.+exp2ipifL-2.*exp4ipifL) + (g31mode-g13mode)*(1.-2.*exp2ipifL+exp4ipifL) + (g12mode-g32mode)*(-2.+exp2ipifL+exp4ipifL) ) * amphOtilde; */
      /* campT = commonfac/sqrt(3) * (g31mode-g13mode+g12mode-g21mode+g23mode-g32mode) * (1.+exp2ipifL+exp4ipifL) * amphOtilde; */
      /**/
      double complex expipifL = cexp(I*PI*f*L_SI/C_SI);
      double complex exp2ipifL = expipifL*expipifL;
      double complex exp4ipifL = exp2ipifL*exp2ipifL;
      double complex commonfac = -2*I*exp4ipifL;
      double complex amphOtilde = gsl_vector_get(h0tilde_ampreal, j) + I * gsl_vector_get(h0tilde_ampimag, j);
      /**/
      campA = commonfac/sqrt(2) * ( (g21mode+g23mode)*(exp2ipifL+1.) - exp2ipifL*(g32mode+g12mode) + (g13mode+g31mode) ) * amphOtilde;
      campE = commonfac/sqrt(6) * ( (g31mode-g13mode)*(2*exp2ipifL+1.) + (g32mode-g12mode)*(exp2ipifL+2.) + (g21mode-g23mode)*(exp2ipifL-1.) ) * amphOtilde;
      campT = commonfac/sqrt(3) * expipifL * (g13mode-g31mode+g21mode-g12mode+g32mode-g23mode) * amphOtilde;
      /**/
      gsl_vector_set(amp_realA, j, creal(campA));
      gsl_vector_set(amp_imagA, j, cimag(campA));
      gsl_vector_set(amp_realE, j, creal(campE));
      gsl_vector_set(amp_imagE, j, cimag(campE));
      gsl_vector_set(amp_realT, j, creal(campT));
      gsl_vector_set(amp_imagT, j, cimag(campT));
    }
    //clock_t tendcontesllation = clock();
  //printf("Set constellation time: %g s\n", (double)(tendcontesllation - tbegcontesllation) / CLOCKS_PER_SEC);
  //printf("GAB cumulated time: %g s\n", timingcumulativeGABmode);

    /* Copying the vectors of frequencies and phases */
    gsl_vector_memcpy(freqA, freq);
    gsl_vector_memcpy(freqE, freq);
    gsl_vector_memcpy(freqT, freq);
    gsl_vector_memcpy(phaseA, phase);
    gsl_vector_memcpy(phaseE, phase);
    gsl_vector_memcpy(phaseT, phase);
    
    /* Append the modes to the ouput list-of-modes structures */
    *listA = ListmodesCAmpPhaseFrequencySeries_AddModeNoCopy(*listA, modefreqseriesA, l, m);
    *listE = ListmodesCAmpPhaseFrequencySeries_AddModeNoCopy(*listE, modefreqseriesE, l, m);
    *listT = ListmodesCAmpPhaseFrequencySeries_AddModeNoCopy(*listT, modefreqseriesT, l, m);

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
