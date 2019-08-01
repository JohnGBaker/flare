/**
 * \author Sylvain Marsat, University of Maryland - NASA GSFC
 *
 * \brief C header for functions windowing and computing FFT/IFFT of time/frequency series.
 *
 */

#ifndef __COMPUTELISASNR_H__
#define __COMPUTELISASNR_H__ 1

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#define _XOPEN_SOURCE 500

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
#include "waveform.h"
#include "fft.h"
#include "LISAgeometry.h"
#include "LISAFDresponse.h"
#include "LISAutils.h"

/********************************** Structures ******************************************/

/* Parameters for the generation of a LISA waveform (in the form of a list of modes) */
typedef struct tagComputeLISASNRparams {
  double tRef;               /* reference time (s) - GPS time at the frequency representing coalescence */
  double phiRef;             /* reference phase (rad) - phase at the frequency representing coalescence (or at fRef if specified) */
  double fRef;               /* reference frequency at which phiRef is set (Hz, default 0 which is interpreted as Mf=0.14) */
  double m1;                 /* mass of companion 1 (solar masses, default 2e6) */
  double m2;                 /* mass of companion 2 (solar masses, default 1e6) */
  double distance;           /* distance of source (Mpc, default 1e3) */
  double inclination;        /* inclination of source (rad, default pi/3) */
  double lambda;             /* first angle for the position in the sky (rad, default 0) */
  double beta;               /* second angle for the position in the sky (rad, default 0) */
  double polarization;       /* polarization angle (rad, default 0) */

  int nbmode;                /* number of modes to generate (starting with 22) - defaults to 5 (all modes) */
  double deltatobs;          /* Observation duration (years, default=2) */
  double minf;               /* Minimal frequency (Hz, default=0) - when set to 0, use the lowest frequency where the detector noise model is trusted __LISASimFD_Noise_fLow (set somewhat arbitrarily)*/
  double maxf;               /* Maximal frequency (Hz, default=0) - when set to 0, use the highest frequency where the detector noise model is trusted __LISASimFD_Noise_fHigh (set somewhat arbitrarily)*/
  int tagextpn;              /* Tag to allow PN extension of the waveform at low frequencies */
  double Mfmatch;            /* When PN extension allowed, geometric matching frequency: will use ROM above this value. If <=0, use ROM down to the lowest covered frequency */
  int tagtdi;                /* Tag selecting the desired TDI observables */
  int tagint;                /* Tag choosing the integrator: 0 for Fresnel (default), 1 for linear integration */
  int nbptsoverlap;          /* Number of points to use in loglinear overlaps (default 32768) */
  LISAconstellation *variant;  /* A structure defining the LISA constellation features */
  int frozenLISA;            /* Freeze the orbital configuration (default 0) */
  double tfrozenLISA;        /* Time in s at which to freeze LISA (default 0.) */
  ResponseApproxtag responseapprox;    /* Approximation in the GAB and orb response - choices are full (full response, default), lowfL (keep orbital delay frequency-dependence but simplify constellation response) and lowf (simplify constellation and orbital response) - WARNING : at the moment noises are not consistent, and TDI combinations from the GAB are unchanged */
  int delaycorrection;       /* Include the first-order ddot delay correction in phaseRdelay (default 1) - NOTE: treated separately from frozenLISA, while strictly speaking ddot should be zero for a frozen LISA */
  int fromtditdfile;         /* Option for loading time series for TDI observables and FFTing */
  int nlinesinfile;          /* Number of lines of input file */
  char indir[256];           /* Input directory */
  char infile[256];          /* Input file */
  int loadparamsfile;        /* Option to load physical parameters from file and to output result to file (default 0) */
  int nlinesparams;          /* Number of lines in params file */
  char paramsdir[256];       /* Directory for the input/output file */
  char paramsfile[256];      /* Input file with the parameters */
  char outputfile[256];      /* Output file */
} ComputeLISASNRparams;


#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _COMPUTELISASNR_H */
