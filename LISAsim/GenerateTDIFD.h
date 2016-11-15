/**
 * \author Sylvain Marsat, University of Maryland - NASA GSFC
 *
 * \brief C header for functions windowing and computing FFT/IFFT of time/frequency series.
 *
 */

#ifndef __GENERATETDIFD_H__
#define __GENERATETDIFD_H__ 1

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
#include "LISAnoise.h"

/********************************** Structures ******************************************/

typedef enum GenTDIFDtag {
  TDIFD,
  TDIhlm,
  h22FFT,
  yslrFFT
} GenTDIFDtag;

/* Parameters for the generation of a LISA waveform (in the form of a list of modes) */
typedef struct tagGenTDIFDparams {
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
  double minf;               /* Minimal frequency (Hz, default=0) - when set to 0, use the lowest frequency where the detector noise model is trusted __LISASimFD_Noise_fLow (set somewhat arbitrarily)*/
  double maxf;               /* Maximal frequency (Hz, default=0) - when set to 0, use the highest frequency where the detector noise model is trusted __LISASimFD_Noise_fHigh (set somewhat arbitrarily)*/
  double deltatobs;          /* Observation duration (years, default=2) */
  int tagextpn;              /* Tag to allow PN extension of the waveform at low frequencies (default=1) */
  double Mfmatch;            /* When PN extension allowed, geometric matching frequency: will use ROM above this value. If <=0, use ROM down to the lowest covered frequency (default=0.) */
  double deltaf;             /* When generating frequency series from the mode contributions, deltaf for the output (0 to set automatically at 1/2*1/(2T)) */
  double twindowbeg;         /* When generating frequency series from file by FFT, twindowbeg (0 to set automatically at 0.05*duration) */
  double twindowend;         /* When generating frequency series from file by FFT, twindowend (0 to set automatically at 0.01*duration) */
  int tagtdi;                /* Tag selecting the desired output format */
  int taggenwave;            /* Tag selecting the desired output format */
  int restorescaledfactor;   /* If 1, restore the factors that were scaled out of TDI observables */
  int fromtditdfile;         /* Tag for loading time series for TDI observables and FFTing */
  int nsamplesinfile;        /* Number of lines of input file */
  int binaryin;              /* Tag for loading the data in gsl binary form instead of text (default false) */
  int binaryout;             /* Tag for outputting the data in gsl binary form instead of text (default false) */
  char indir[256];           /* Path for the input directory */
  char infile[256];          /* Path for the input file */
  char outdir[256];          /* Path for the output directory */
  char outfileprefix[256];   /* Path for the output file */
} GenTDIFDparams;


#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _GENERATETDIFD_H */
