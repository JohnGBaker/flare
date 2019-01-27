/**
 * \author Sylvain Marsat, University of Maryland - NASA GSFC
 *
 * \brief C header for functions manipulating waveforms.
 *
 */

#ifndef _WAVEFORM_H
#define _WAVEFORM_H

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
#include "EOBNRv2HMROM.h"


#if defined(__cplusplus)
extern "C" {
#define complex _Complex
#elif 0
} /* so that editors will match preceding brace */
#endif

/* NOTE: uses the list of modes of EOBNRv2HMROM (listmode), to be extended when more waveform models are added */

/***************** Function estimating frequency corresponding to a given time to coalescence ****************/

/* Functions computing relations between chirp mass, eta, m1, m2*/
double Mchirpofm1m2(
  const double m1,   /* Mass 1 */
  const double m2);  /* Mass 2 */
double etaofm1m2(
  const double m1,   /* Mass 1 */
  const double m2);  /* Mass 2 */
double m1ofMchirpeta(
  const double Mchirp,   /* Chirp mass */
  const double eta);     /* Symmetric mass ratio */
double m2ofMchirpeta(
  const double Mchirp,   /* Chirp mass */
  const double eta);     /* Symmetric mass ratio */

/* Newtonian estimate of the relation Mf(deltat/M) (for the 22 mode) - gives the starting geometric frequency for a given mass ratio and a given geometric duration of the observations */
double NewtonianfoftGeom(const double q, const double t); /* t here is t/M */
/* Newtonian estimate of the relation f(deltat) (for the 22 mode) - gives the starting frequency in Hz for a given mass ratio and a given geometric duration of the observations */
double Newtonianfoft(const double m1, const double m2, const double t); /* t here is in years, m1-m2 in solar masses */
/* Newtonian estimate of the relation f(deltat) (for the 22 mode freq) - gives the starting geometric frequency for a given time to merger and chirp mass - output in Hz */
double Newtonianfoftchirp(
  const double mchirp,                 /* Chirp mass (solar masses) */
  const double t);                     /* Time in years */
/* Newtonian estimate of the relation deltat(f) (for the 22 mode freq) - gives the time to merger from a starting frequency for a given chirp mass - output in years */
double Newtoniantoffchirp(
  const double mchirp,                 /* Chirp mass (solar masses) */
  const double t);                     /* Freq in Hz */

/***************** Function estimating time to coalescence and min/max frequency ****************/

/* Functions reading from a list of modes the minimal and maximal frequencies */
double ListmodesCAmpPhaseFrequencySeries_maxf(ListmodesCAmpPhaseFrequencySeries* listhlm);
double ListmodesCAmpPhaseFrequencySeries_minf(ListmodesCAmpPhaseFrequencySeries* listhlm);

/* Function estimating initial time from Psi22, according to tf_SPA = -1/(2pi)dPsi/df */
double EstimateInitialTime(ListmodesCAmpPhaseFrequencySeries* listhlm, double fLow);

/***************** Functions to manipulate ReImFrequencySeries structure ****************/

void ReImFrequencySeries_AddCAmpPhaseFrequencySeries(
  struct tagReImFrequencySeries* freqseriesReIm,              /* Output Re/Im frequency series */
  struct tagCAmpPhaseFrequencySeries* freqseriesCAmpPhase,    /* Input CAmp/Phase frequency series, to be interpolated and added to the output */
  double fLow,                                                /* Minimal frequency - set to 0 to ignore */
  double fHigh,                                               /* Maximal frequency - set to 0 to ignore */
  double fstartobsmode);                                      /* Starting frequency in case of limited duration of observations- assumed to have been scaled with the proper factor m/2 for this mode - set to 0 to ignore */
/* Function evaluating a ReImFrequencySeries by interpolating wach mode of a ListmodesCAmpPhaseFrequencySeries and summing them, given a set of frequencies */
void ReImFrequencySeries_SumListmodesCAmpPhaseFrequencySeries(
  struct tagReImFrequencySeries* freqseriesReIm,                    /* Output Re/Im frequency series - already initialized */
  struct tagListmodesCAmpPhaseFrequencySeries* listmodesCAmpPhase,  /* Input CAmp/Phase frequency series, to be interpolated */
  gsl_vector* freq,                                                 /* Input set of frequencies on which evaluating */
  double fLow,                                                      /* Minimal frequency - set to 0 to ignore */
  double fHigh,                                                     /* Maximal frequency - set to 0 to ignore */
  double fstartobs);                                                /* For limited duration of observation, starting frequency for the 22 mode - set to 0 to ignore */
/* Helper function to add a mode to hplus, hcross in Fourier domain
 * - copies the function XLALSimAddMode, which was done only for TD structures */
int FDAddMode(
  ReImFrequencySeries* hptilde,       /* Output: frequency series for hplus */
  ReImFrequencySeries* hctilde,       /* Output: frequency series for hcross */
  ReImFrequencySeries* hlmtilde,      /* Input: frequency series for the mode hlm */
  double theta,                       /* First angle for position in the sky of observer */
  double phi,                         /* Second angle for position in the sky of observer  */
  int l,                              /* First mode number l */
  int m,                              /* Second mode number m */
  int sym);                           /* If 1, assume planar symmetry and add also mode l,-m. Do not if set to 0. */
/* Function evaluating the FD frequency series for hplus, hcross from the modes hlm */
int GeneratehphcFDReImFrequencySeries(
  ReImFrequencySeries** hptilde,                 /* Output: frequency series for hplus */
  ReImFrequencySeries** hctilde,                 /* Output: frequency series for hcross */
  ListmodesCAmpPhaseFrequencySeries* listhlm,    /* Input: frequency series for the mode hlm  */
  double fLow,                                   /* Minimal frequency - set to 0 to ignore */
  double fHigh,                                  /* Maximal frequency - set to 0 to ignore */
  double fstartobs,                              /* For limited duration of observation, starting frequency for the 22 mode - set to 0 to ignore */
  double deltaf,                                 /* Frequency step */
  int nbpt,                                      /* Number of points of output - if 0, determined from deltaF and maximal frequency in input */
  int nbmode,                                    /* Number of modes to add */
  double theta,                                  /* First angle for position in the sky of observer */
  double phi,                                    /* Second angle for position in the sky of observer */
  int sym);                                      /* If 1, assume planar symmetry and add also mode l,-m. Do not if set to 0. */
/* Function evaluating the FD frequency series by summing mode contributions from each hlm */
int GenerateFDReImFrequencySeries(
  ReImFrequencySeries** freqseries,              /* Output: frequency series */
  ListmodesCAmpPhaseFrequencySeries* listhlm,    /* Input: FD modes hlm in the form AmpReal/AmpIm/Phase  */
  double fLow,                                   /* Minimal frequency - set to 0 to ignore */
  double fHigh,                                  /* Maximal frequency - set to 0 to ignore */
  double fstartobs,                              /* For limited duration of observation, starting frequency for the 22 mode - set to 0 to ignore */
  double deltaf,                                 /* Frequency step */
  int nbpt);                                     /* Number of points of output - if 0, determined from deltaf and maximal frequency in input */
/* Function evaluating the FD frequency series for a single mode contribution (l,m) */
int GenerateFDReImFrequencySeriesSingleMode(
  ReImFrequencySeries** freqseries,              /* Output: frequency series */
  ListmodesCAmpPhaseFrequencySeries* listhlm,    /* Input: FD modes hlm in the form AmpReal/AmpIm/Phase  */
  double fLow,                                   /* Minimal frequency - set to 0 to ignore */
  double fHigh,                                  /* Maximal frequency - set to 0 to ignore */
  double fstartobs,                              /* For limited duration of observation, starting frequency for the 22 mode - set to 0 to ignore */
  double deltaf,                                 /* Frequency step */
  int nbpt,                                      /* Number of points of output - if 0, determined from deltaf and maximal frequency in input */
  int l,                                         /* Mode index l */
  int m);                                        /* Mode index m */
/* Function to restrict a frequency series (typically output of a FFT) to a given frequency range */
int RestrictFDReImFrequencySeries(
  ReImFrequencySeries** freqseriesout,           /* Output: truncated frequency series */
  ReImFrequencySeries* freqseriesin,             /* Input: frequency series */
  double fLow,                                   /* Minimal frequency */
  double fHigh);                                 /* Maximal frequency */
/* Function for getting a phase mod 2pi (between -pi and pi) */
double mod2pi(double phase);
/* Function for getting a phase mod pi (between 0 and pi, e.g. polarization) */
double modpi(double phase);
/* Function to unwrap the phase mod 2pi  - acts directly on the gsl_vector representing the phase */
int UnwrapPhase(
  gsl_vector*  phaseout,   /* Output: unwrapped phase vector - already allocated */
  gsl_vector*  phasein);   /* Input: phase vector */
/* Function to convert a time series from Re/Im form to Amp/Phase form - unwrapping the phase */
int ReImTimeSeries_ToAmpPhase(
  AmpPhaseTimeSeries** timeseriesout,             /* Output: Amp/Phase time series */
  ReImTimeSeries* timeseriesin);                  /* Input: Re/Im time series */
/* Function to convert a time series from Amp/Phase form to Re/Im form */
int AmpPhaseTimeSeries_ToReIm(
  ReImTimeSeries** timeseriesout,                 /* Output: Re/Im time series */
  AmpPhaseTimeSeries* timeseriesin);              /* Input: Amp/Phase time series */
/* Function to compute a linear resampling at high frequencies to enforce a maximal deltaf */
/* NOTE: Assumes input frequencies are logarithmic (except maybe first interval) to evaluate when to resample */
int SetMaxdeltafResampledFrequencies(
  gsl_vector** freqr,              /* Output: resampled frequencies */
  gsl_vector* freq,                /* Input: original frequencies */
  const double maxf,               /* Input: maximal frequency - set to 0. to ignore */
  const double deltaf);            /* Input: maximal deltaf aimed for - 0.002Hz appropriate for LISA */
/* Function to compute a linear-in-time resampling at low frequencies to enforce a maximal deltat */
int SetMaxdeltatResampledFrequencies(
  gsl_vector** freqr,              /* Output: resampled frequencies */
  gsl_vector* freq,                /* Input: original frequencies */
  const double deltat,             /* Input: maximal deltat aimed for - fraction of a year, 1/24 (half month) appropriate for 1e-4 interpolation errors */
  const double mchirp,             /* Input: chirp mass, used for approximate t-f correspondence */
  const int m);                    /* Input: chirp mass, used for approximate t-f correspondence */
/* Function to resample a CAmp/Phase frequency series on the specified frequencies */
int CAmpPhaseFrequencySeries_Resample(
  CAmpPhaseFrequencySeries** freqseriesout,         /* Output: CAmp/Phase freq series */
  CAmpPhaseFrequencySeries* freqseriesin,           /* Input: CAmp/Phase freq series */
  gsl_vector* freqr);                                 /* Input: freq vector to resample on */

/***************** Spin weighted spherical harmonics ****************/

/* Additional function reproducing XLALSpinWeightedSphericalHarmonic */
double complex SpinWeightedSphericalHarmonic(double theta, double phi, int s, int l, int m); /* Currently only supports s=-2, l=2,3,4,5 modes */

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _WAVEFORM_H */
