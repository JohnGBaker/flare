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
/* Older version of th FD response: two stages, first for the orbital delay (Bessel phase) and second for the constellation delay/modulation */
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
int LISASimFDResponsey12(
  struct tagListmodesCAmpPhaseFrequencySeries **list,      /* Input: list of modes in Frequency-domain amplitude and phase form as produced by the ROM */
  struct tagListmodesCAmpPhaseFrequencySeries **listy12,   /* Output: list of contribution of each mode in Frequency-domain amplitude and phase form, in the y12 observable */
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
    double f, tf;
    double complex g12mode = 0.;
    double complex g21mode = 0.;
    double complex g23mode = 0.;
    double complex g32mode = 0.;
    double complex g31mode = 0.;
    double complex g13mode = 0.;
    double complex camp_y;

    /* Computing the Ylm combined factors for plus and cross for this mode */
    /* Capital Phi is set to 0 by convention */
    double complex Yfactorplus;
    double complex Yfactorcross;
    if (!(l%2)) {
      Yfactorplus = 1./2 * (SpinWeightedSphericalHarmonic(inclination, 0., -2, l, m) + conj(SpinWeightedSphericalHarmonic(inclination, 0., -2, l, -m)));
      Yfactorcross = I/2 * (SpinWeightedSphericalHarmonic(inclination, 0., -2, l, m) - conj(SpinWeightedSphericalHarmonic(inclination, 0., -2, l, -m)));
    }
    else {
      Yfactorplus = 1./2 * (SpinWeightedSphericalHarmonic(inclination, 0., -2, l, m) - conj(SpinWeightedSphericalHarmonic(inclination, 0., -2, l, -m)));
      Yfactorcross = I/2 * (SpinWeightedSphericalHarmonic(inclination, 0., -2, l, m) + conj(SpinWeightedSphericalHarmonic(inclination, 0., -2, l, -m)));
    }

    /* Initializing spline for the phase */
    gsl_spline* spline_phi = gsl_spline_alloc(gsl_interp_cspline, len);
    gsl_interp_accel* accel_phi = gsl_interp_accel_alloc();
    gsl_spline_init(spline_phi, gsl_vector_const_ptr(freq, 0), gsl_vector_const_ptr(phase, 0), len);

    /* Orbital delay, constellation delay/modulation */
    /* Initializing frequency series structure for this mode, for each of the TDI observables */
    CAmpPhaseFrequencySeries *modefreqseries = NULL;
    CAmpPhaseFrequencySeries_Init(&modefreqseries, len);
    gsl_vector* freq_y = modefreqseries->freq;
    gsl_vector* amp_real_y = modefreqseries->amp_real;
    gsl_vector* amp_imag_y = modefreqseries->amp_imag;
    gsl_vector* phase_y = modefreqseries->phase;
    /* Loop over the frequencies - computing the correction due to the constellation delay/modulation */
    //clock_t tbegcontesllation = clock();
    //double timingcumulativeGABmode = 0;
    for(int j=0; j<len; j++) {
      f = gsl_vector_get(freq, j);
      tf = (gsl_spline_eval_deriv(spline_phi, f, accel_phi))/(2*PI);
      //clock_t tbegGAB = clock();
      EvaluateGABmode(&g12mode, &g21mode, &g23mode, &g32mode, &g31mode, &g13mode, f, tf, Yfactorplus, Yfactorcross, 1);
      //clock_t tendGAB = clock();
      //timingcumulativeGABmode += (double) (tendGAB-tbegGAB) /CLOCKS_PER_SEC;
      /**/
      camp_y = (gsl_vector_get(amp_real, j) + I * gsl_vector_get(amp_imag, j)) * g12mode;
      /**/
      gsl_vector_set(amp_real_y, j, creal(camp_y));
      gsl_vector_set(amp_imag_y, j, cimag(camp_y));
    }
    //clock_t tendcontesllation = clock();
  //printf("Set constellation time: %g s\n", (double)(tendcontesllation - tbegcontesllation) / CLOCKS_PER_SEC);
  //printf("GAB cumulated time: %g s\n", timingcumulativeGABmode);

    /* Copying the vectors of frequencies and phases */
    gsl_vector_memcpy(freq_y, freq);
    gsl_vector_memcpy(phase_y, phase);

    /* Append the modes to the ouput list-of-modes structures */
    *listy12 = ListmodesCAmpPhaseFrequencySeries_AddModeNoCopy(*listy12, modefreqseries, l, m);

    /* Going to the next mode in the list */
    listelement = listelement->next;

    /* Clean up */
    gsl_spline_free(spline_phi);
    gsl_interp_accel_free(accel_phi);
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
  const double tRef,                                       /* Time at coalescence */
  const double lambda,                                     /* First angle for the position in the sky */
  const double beta,                                       /* Second angle for the position in the sky */
  const double inclination,                                /* Inclination of the source */
  const double psi,                                        /* Polarization angle */
  const double maxf,                                       /* Maximal frequency to consider - used to ignore hard-to-resolve response at f>1Hz - NOTE: for now, no recomputation of the boundary, so when not resampling can lose a bit of support between the last frequency point covered and maxf */
  const TDItag tditag)                                     /* Selector for the set of TDI observables */
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
    //
    //if(l==2&&m==1) printf("no resampling 21\n");
    //      if(l==2&&m==1) {for(int i=0; i<len; i++) printf("%d, %g, %g, %g, %g\n", i, gsl_vector_get(freq, i), gsl_vector_get(amp_real, i), gsl_vector_get(amp_imag, i), gsl_vector_get(phase, i));};

    double f, tf;
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
      Yfactorplus = 1./2 * (SpinWeightedSphericalHarmonic(inclination, 0., -2, l, m) + conj(SpinWeightedSphericalHarmonic(inclination, 0., -2, l, -m)));
      Yfactorcross = I/2 * (SpinWeightedSphericalHarmonic(inclination, 0., -2, l, m) - conj(SpinWeightedSphericalHarmonic(inclination, 0., -2, l, -m)));
    }
    else {
      Yfactorplus = 1./2 * (SpinWeightedSphericalHarmonic(inclination, 0., -2, l, m) - conj(SpinWeightedSphericalHarmonic(inclination, 0., -2, l, -m)));
      Yfactorcross = I/2 * (SpinWeightedSphericalHarmonic(inclination, 0., -2, l, m) + conj(SpinWeightedSphericalHarmonic(inclination, 0., -2, l, -m)));
    }

    /* Orbital delay, constellation delay/modulation */

    /* Initializing spline for the phase */
    gsl_spline* spline_phi = gsl_spline_alloc(gsl_interp_cspline, len);
    gsl_interp_accel* accel_phi = gsl_interp_accel_alloc();
    gsl_spline_init(spline_phi, gsl_vector_const_ptr(freq, 0), gsl_vector_const_ptr(phase, 0), len);

    /* Determine frequencies to use */
    /* Because we will need the interpolation on the Re/Im amplitude to resolve the structure of the L-response at high fequencies, we add points beyond 0.01 Hz, with a linear sampling */
    /* WARNING : It seemed 600 points between 0.1 and 3 Hz (deltaf=0.005) should give interpolation errors below 1e-4 - considering a simple sin(2 pi f L) */
    /* But first test with TDIA show the sampling should be reduced 10-fold - possibly large increase in cost */
    /* We use here deltaf=0.0005 until the cause is better understood */
    /* NOTE: the structure in the response due to the R-delay gives much larger interpolation errors - here we assume the R-delay term is now treated as a phase */
    int resampled = 0; /* Keeps track of wether or not we resampled and allocated new resources we need to free */
    double maxfsignal = gsl_vector_get(freq, len-1);
    double fHigh = fmin(maxf, maxfsignal);
    double deltaflineartarget = 0.0005;
    int imaxlogsampling = len-1; /* last index to be covered with original sampling */
    while((gsl_vector_get(freq, imaxlogsampling)>fHigh) && imaxlogsampling>0) imaxlogsampling--;
    gsl_vector* freq_resample = NULL; /*  */
    gsl_vector* amp_real_resample = NULL;
    gsl_vector* amp_imag_resample = NULL;
    gsl_vector* phase_resample = NULL;
    int len_resample;
    //
    //printf("%d\n", imaxlogsampling);
    //printf("%g, %g, %g\n", maxf, maxfsignal, fHigh);
    /* We keep original logarithmic sampling below 0.1Hz */
    if((fHigh>0.1) && ((fHigh - gsl_vector_get(freq, imaxlogsampling))/deltaflineartarget)>len-1-imaxlogsampling ) { /* condition to check if the linear sampling will add points - if not, do nothing */
      resampled = 1;
      /* Number of pts in resampled part */
      int nbfreqlinear = ceil((fHigh - gsl_vector_get(freq, imaxlogsampling))/deltaflineartarget);
      double deltaflinear = (fHigh - gsl_vector_get(freq, imaxlogsampling))/(nbfreqlinear + 1);
      /* Initialize new vectors */
      len_resample = imaxlogsampling + 1 + nbfreqlinear;
      freq_resample = gsl_vector_alloc(len_resample);
      amp_real_resample = gsl_vector_alloc(len_resample);
      amp_imag_resample = gsl_vector_alloc(len_resample);
      phase_resample = gsl_vector_alloc(len_resample);
      /* Build interpolation for original amp_real, amp_imag and phase */
      /* NOTE: we could use spline_phi here, written this way for clarity */
      gsl_spline* spline_amp_real = gsl_spline_alloc(gsl_interp_cspline, len);
      gsl_spline* spline_amp_imag = gsl_spline_alloc(gsl_interp_cspline, len);
      gsl_spline* spline_phase = gsl_spline_alloc(gsl_interp_cspline, len);
      gsl_interp_accel* accel_amp_real = gsl_interp_accel_alloc();
      gsl_interp_accel* accel_amp_imag = gsl_interp_accel_alloc();
      gsl_interp_accel* accel_phase = gsl_interp_accel_alloc();
      gsl_spline_init(spline_amp_real, gsl_vector_const_ptr(freq, 0), gsl_vector_const_ptr(amp_real, 0), len);
      gsl_spline_init(spline_amp_imag, gsl_vector_const_ptr(freq, 0), gsl_vector_const_ptr(amp_imag, 0), len);
      gsl_spline_init(spline_phase, gsl_vector_const_ptr(freq, 0), gsl_vector_const_ptr(phase, 0), len);
      /* Set resampled frequencies and values */
      for(int j=0; j<=imaxlogsampling; j++) {
        gsl_vector_set(freq_resample, j, gsl_vector_get(freq, j));
        gsl_vector_set(amp_real_resample, j, gsl_vector_get(amp_real, j));
        gsl_vector_set(amp_imag_resample, j, gsl_vector_get(amp_imag, j));
        gsl_vector_set(phase_resample, j, gsl_vector_get(phase, j));
      }
      double fimax = gsl_vector_get(freq, imaxlogsampling);
      for(int j=imaxlogsampling+1; j<len_resample; j++) {
        f = fimax + (j-imaxlogsampling) * deltaflinear;
        gsl_vector_set(freq_resample, j, f);
        gsl_vector_set(amp_real_resample, j, gsl_spline_eval(spline_amp_real, f, accel_amp_real));
        gsl_vector_set(amp_imag_resample, j, gsl_spline_eval(spline_amp_imag, f, accel_amp_imag));
        gsl_vector_set(phase_resample, j, gsl_spline_eval(spline_phase, f, accel_phase));
      }
      /* Free interpolation functions */
      gsl_spline_free(spline_amp_real);
      gsl_spline_free(spline_amp_imag);
      gsl_spline_free(spline_phase);
      gsl_interp_accel_free(accel_amp_real);
      gsl_interp_accel_free(accel_amp_imag);
      gsl_interp_accel_free(accel_phase);
    }
    else { /* If no resampling, use the values we had as input */
      resampled = 0;
      len_resample = imaxlogsampling + 1; /* If maxf < maxfsignal, we will cut the signal above maxf by simply adjusting the range of indices included (NOTE: without recomputing the exact boundary, so some support is lost due to discretization) */
      freq_resample = freq;
      amp_real_resample = amp_real;
      amp_imag_resample = amp_imag;
      phase_resample = phase;
    }

    /* Initializing frequency series structure for this mode, for each of the TDI observables */
    CAmpPhaseFrequencySeries *modefreqseries1 = NULL;
    CAmpPhaseFrequencySeries *modefreqseries2 = NULL;
    CAmpPhaseFrequencySeries *modefreqseries3 = NULL;
    CAmpPhaseFrequencySeries_Init(&modefreqseries1, len_resample);
    CAmpPhaseFrequencySeries_Init(&modefreqseries2, len_resample);
    CAmpPhaseFrequencySeries_Init(&modefreqseries3, len_resample);
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

    /* Loop over the frequencies */
    //clock_t tbegcontesllation = clock();
    //double timingcumulativeGABmode = 0;
    //
    //printf("loop\n");
      //
      //printf("len_resample: %d\n", len_resample);
    for(int j=0; j<len_resample; j++) {

        //
        //printf("j: %d\n", j);
      f = gsl_vector_get(freq_resample, j);
      tf = (gsl_spline_eval_deriv(spline_phi, f, accel_phi))/(2*PI);
      //clock_t tbegGAB = clock();
      EvaluateGABmode(&g12mode, &g21mode, &g23mode, &g32mode, &g31mode, &g13mode, f, tf, Yfactorplus, Yfactorcross, 0); /* does not include the R-delay term */
      //clock_t tendGAB = clock();
      //timingcumulativeGABmode += (double) (tendGAB-tbegGAB) /CLOCKS_PER_SEC;
      /**/
      EvaluateTDIfactor3Chan(&factor1, &factor2, &factor3, g12mode, g21mode, g23mode, g32mode, g31mode, g13mode, f, tditag);
      double complex amphtilde = gsl_vector_get(amp_real_resample, j) + I * gsl_vector_get(amp_imag_resample, j);
      camp1 = factor1 * amphtilde;
      camp2 = factor2 * amphtilde;
      camp3 = factor3 * amphtilde;
      /* Phase term due to the R-delay, including correction to first order */
      double phaseRdelay = -2*PI*R_SI/C_SI*f*cos(beta)*cos(Omega_SI*tf - lambda) * (1 + R_SI/C_SI*Omega_SI*sin(Omega_SI*tf - lambda));
      double phasewithRdelay = gsl_vector_get(phase_resample, j) + phaseRdelay;

    //
    //if(l==2&&m==1) printf("j, gsl_vector_get(phase_resample, j): %d, %g\n", j, gsl_vector_get(phase_resample, j));
      /**/
      gsl_vector_set(amp_real1, j, creal(camp1));
      gsl_vector_set(amp_imag1, j, cimag(camp1));
      gsl_vector_set(amp_real2, j, creal(camp2));
      gsl_vector_set(amp_imag2, j, cimag(camp2));
      gsl_vector_set(amp_real3, j, creal(camp3));
      gsl_vector_set(amp_imag3, j, cimag(camp3));
      gsl_vector_set(phase1, j, phasewithRdelay);
      gsl_vector_set(phase2, j, phasewithRdelay);
      gsl_vector_set(phase3, j, phasewithRdelay);
    }
    //clock_t tendconstellation = clock();
  //printf("Set constellation time: %g s\n", (double)(tendconstellation - tbegconstellation) / CLOCKS_PER_SEC);
  //printf("GAB cumulated time: %g s\n", timingcumulativeGABmode);

    /* Copying the vectors of frequencies - we have to allow for the case where it has been shortened */
    gsl_vector_view freq_resample_subview = gsl_vector_subvector(freq_resample, 0, len_resample);
    gsl_vector_memcpy(freq1, &freq_resample_subview.vector);
    gsl_vector_memcpy(freq2, &freq_resample_subview.vector);
    gsl_vector_memcpy(freq3, &freq_resample_subview.vector);

    /* Append the modes to the ouput list-of-modes structures */
    *listTDI1 = ListmodesCAmpPhaseFrequencySeries_AddModeNoCopy(*listTDI1, modefreqseries1, l, m);
    *listTDI2 = ListmodesCAmpPhaseFrequencySeries_AddModeNoCopy(*listTDI2, modefreqseries2, l, m);
    *listTDI3 = ListmodesCAmpPhaseFrequencySeries_AddModeNoCopy(*listTDI3, modefreqseries3, l, m);

    /* Going to the next mode in the list */
    listelement = listelement->next;

    /* Clean up */
    gsl_spline_free(spline_phi);
    gsl_interp_accel_free(accel_phi);
    /* If we used resampling, then we need to free the additional resources that were allocated */
    if(resampled) {
      gsl_vector_free(freq_resample);
      gsl_vector_free(amp_real_resample);
      gsl_vector_free(amp_imag_resample);
      gsl_vector_free(phase_resample);
    }
  }

  //
  //printf("end fdresponse\n");

  return SUCCESS;
}
