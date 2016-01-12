/**
 * \author Sylvain Marsat, University of Maryland - NASA GSFC
 *
 * \brief C code headers for the geometric coefficients entering the response for LISA-like detectors.
 *
 */

#ifndef _LISAGEOMETRY_H
#define _LISAGEOMETRY_H

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

/*****************************************************/
/**************** TDI variables **********************/

typedef enum TDItag {
  TDIXYZ,
  TDIalphabetagamma,
  TDIAETXYZ,
  TDIAETalphabetagamma,
  TDIX,
  TDIalpha,
  TDIAXYZ,
  TDIEXYZ,
  TDITXYZ,
  TDIAalphabetagamma,
  TDIEalphabetagamma,
  TDITalphabetagamma
} TDItag;

/**************************************************/
/**************** Prototypes **********************/

/* Function to convert string input TDI string to TDItag */
TDItag ParseTDItag(char* string);
/* Function returning the number of channels for a TDItag */
int nbchanTDI(TDItag tditag);

/* Function cardinal sine */
double sinc(const double x);

/* Function to compute, given a value of a sky position and polarization, all the complicated time-independent trigonometric coefficients entering the response */
void SetCoeffsG(const double lambda, const double beta, const double psi);

/* Functions evaluating the G_AB functions, combining the two polarization with the spherical harmonics factors */
double complex G21mode(const double f, const double t, const double complex Yfactorplus, const double complex Yfactorcross);
double complex G12mode(const double f, const double t, const double complex Yfactorplus, const double complex Yfactorcross);
double complex G32mode(const double f, const double t, const double complex Yfactorplus, const double complex Yfactorcross);
double complex G23mode(const double f, const double t, const double complex Yfactorplus, const double complex Yfactorcross);
double complex G13mode(const double f, const double t, const double complex Yfactorplus, const double complex Yfactorcross);
double complex G31mode(const double f, const double t, const double complex Yfactorplus, const double complex Yfactorcross);
int EvaluateGABmode(
  double complex* G12,                     /* Output for G12 */
  double complex* G21,                     /* Output for G21 */
  double complex* G23,                     /* Output for G23 */
  double complex* G32,                     /* Output for G32 */
  double complex* G31,                     /* Output for G31 */
  double complex* G13,                     /* Output for G13 */
  const double f,                          /* Frequency */
  const double t,                          /* Time */
  const double complex Yfactorplus,        /* Spin-weighted spherical harmonic factor for plus */
  const double complex Yfactorcross);      /* Spin-weighted spherical harmonic factor for cross */

/* Functions evaluating the Fourier-domain factors (combinations of the GAB's) for TDI observables */
/* Note: in case only one channel is considered, amplitudes for channels 2 and 3 are simply set to 0 */
/* (allows minimal changes from the old structure that assumed KTV A,E,T - but probably not optimal) */
int EvaluateTDIfactor3Chan(
  double complex* factor1,                       /* Output for factor for TDI channel 1 */
  double complex* factor2,                       /* Output for factor for TDI channel 2 */
  double complex* factor3,                       /* Output for factor for TDI channel 3 */
  const double complex G12,                      /* Input for G12 */
  const double complex G21,                      /* Input for G21 */
  const double complex G23,                      /* Input for G23 */
  const double complex G32,                      /* Input for G32 */
  const double complex G31,                      /* Input for G31 */
  const double complex G13,                      /* Input for G13 */
  const double f,                                /* Frequency */
  const TDItag tditag);                          /* Selector for the TDI observables */
/* int EvaluateTDIfactor1Chan( */
/*   double complex* factor,                       /\* Output for factor for TDI channel *\/ */
/*   const double complex G12,                      /\* Input for G12 *\/ */
/*   const double complex G21,                      /\* Input for G21 *\/ */
/*   const double complex G23,                      /\* Input for G23 *\/ */
/*   const double complex G32,                      /\* Input for G32 *\/ */
/*   const double complex G31,                      /\* Input for G31 *\/ */
/*   const double complex G13,                      /\* Input for G13 *\/ */
/*   const double f,                                /\* Frequency *\/ */
/*   const TDItag tditag);                          /\* Selector for the TDI observable *\/ */

/* Functions for the response in time domain */
double y12TD(
  gsl_spline* splinehp,                    /* Input spline for TD hplus */
  gsl_spline* splinehc,                    /* Input spline for TD hcross */
  gsl_interp_accel* accelhp,               /* Accelerator for hp spline */
  gsl_interp_accel* accelhc,               /* Accelerator for hc spline */
  const double t);                         /* Time */
double y21TD(
  gsl_spline* splinehp,                    /* Input spline for TD hplus */
  gsl_spline* splinehc,                    /* Input spline for TD hcross */
  gsl_interp_accel* accelhp,               /* Accelerator for hp spline */
  gsl_interp_accel* accelhc,               /* Accelerator for hc spline */
  const double t);                         /* Time */
double y23TD(
  gsl_spline* splinehp,                    /* Input spline for TD hplus */
  gsl_spline* splinehc,                    /* Input spline for TD hcross */
  gsl_interp_accel* accelhp,               /* Accelerator for hp spline */
  gsl_interp_accel* accelhc,               /* Accelerator for hc spline */
  const double t);                         /* Time */
double y32TD(
  gsl_spline* splinehp,                    /* Input spline for TD hplus */
  gsl_spline* splinehc,                    /* Input spline for TD hcross */
  gsl_interp_accel* accelhp,               /* Accelerator for hp spline */
  gsl_interp_accel* accelhc,               /* Accelerator for hc spline */
  const double t);                         /* Time */
double y31TD(
  gsl_spline* splinehp,                    /* Input spline for TD hplus */
  gsl_spline* splinehc,                    /* Input spline for TD hcross */
  gsl_interp_accel* accelhp,               /* Accelerator for hp spline */
  gsl_interp_accel* accelhc,               /* Accelerator for hc spline */
  const double t);                         /* Time */
double y13TD(
  gsl_spline* splinehp,                    /* Input spline for TD hplus */
  gsl_spline* splinehc,                    /* Input spline for TD hcross */
  gsl_interp_accel* accelhp,               /* Accelerator for hp spline */
  gsl_interp_accel* accelhc,               /* Accelerator for hc spline */
  const double t);                         /* Time */
int EvaluateTDIXYZTD(
  double* TDIX,                            /* Output: value of TDI observable X */
  double* TDIY,                            /* Output: value of TDI observable Y */
  double* TDIZ,                            /* Output: value of TDI observable Z */
  gsl_spline* splinehp,                    /* Input spline for TD hplus */
  gsl_spline* splinehc,                    /* Input spline for TD hcross */
  gsl_interp_accel* accelhp,               /* Accelerator for hp spline */
  gsl_interp_accel* accelhc,               /* Accelerator for hc spline */
  const double t);                         /* Time */
int EvaluateTDIAETXYZTD(
  double* TDIA,                            /* Output: value of TDI observable X */
  double* TDIE,                            /* Output: value of TDI observable Y */
  double* TDIT,                            /* Output: value of TDI observable Z */
  gsl_spline* splinehp,                    /* Input spline for TD hplus */
  gsl_spline* splinehc,                    /* Input spline for TD hcross */
  gsl_interp_accel* accelhp,               /* Accelerator for hp spline */
  gsl_interp_accel* accelhc,               /* Accelerator for hc spline */
  const double t);                         /* Time */

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LISAGEOMETRY_H */
