/**
 * \author Sylvain Marsat, University of Maryland - NASA GSFC
 *
 * \brief C header for functions windowing and computing FFT/IFFT of time/frequency series.
 *
 */

#ifndef _FFT_H
#define _FFT_H

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

#include <fftw3.h> /* Note: when included AFTER complex.h, fftw_complex type defaults to the native double complex */

#include "constants.h"
#include "struct.h"

/* Window functions */
double WindowFunction(double x, double xi, double xf, double deltaxi, double deltaxf);
double WindowFunctionLeft(double x, double xf, double deltaxf);
double WindowFunctionRight(double x, double xi, double deltaxi);

/* FFT of time series */
/* Note: FFT uses flipped convention (i.e. h(f) = int e^(+2ipift)h(t)) */
int FFTTimeSeries(
  ReImFrequencySeries** freqseries,   /* Output: frequency series */
  RealTimeSeries* timeseries,         /* Input: real time series */
  double twindowbeg,                  /* Extent of the window at beginning (starts at the first point) */
  double twindowend,                  /* Extent of the window at the end (end at the last point) */
  int nzeropad);                      /* For 0-padding: length will be (upper power of 2)*2^nzeropad */

/* IFFT of frequency series */
/* Note: assumes frequency series is FT of real data */
/* Note: FFT uses flipped convention (i.e. h(f) = int e^(+2ipift)h(t)) */
int IFFTFrequencySeriesReal(
  RealTimeSeries** timeseries,       /* Output: real time series*/
  ReImFrequencySeries* freqseries,   /* Input: complex frequency series, assumed to be the FT of a real time series */
  double f1windowbeg,                /* Start of window at the beginning */
  double f2windowbeg,                /* End of window at the beginning */
  double f1windowend,                /* Start of window at the end */
  double f2windowend,                /* End of window at the end */
  int nzeropad);                     /* For 0-padding: length will be (upper power of 2)*2^nzeropad */

/* IFFT of frequency series */
/* Note: assumes frequency series is FT of complex data - produces complex output */
/* Note: FFT uses flipped convention (i.e. h(f) = int e^(+2ipift)h(t)) */
int IFFTFrequencySeries(
  ReImTimeSeries** timeseries,       /* Output: complex time series */
  ReImFrequencySeries* freqseries,   /* Input: complex frequency series, assumed to be the FT of a real time series */
  double f1windowbeg,                /* Start of window at the beginning */
  double f2windowbeg,                /* End of window at the beginning */
  double f1windowend,                /* Start of window at the end */
  double f2windowend,                /* End of window at the end */
  int nzeropad);                     /* For 0-padding: length will be (upper power of 2)*2^nzeropad */

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _FFT_H */
