/**
 * \author Sylvain Marsat, University of Maryland - NASA GSFC
 *
 * \brief C header for functions windowing and computing FFT/IFFT of time/frequency series.
 *
 */

#ifndef __GENERATELLVFD_H__
#define __GENERATELLVFD_H__ 1

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
#include "LLVgeometry.h"
#include "LLVFDresponse.h"

/********************************** Structures ******************************************/

typedef enum GenLLVFDtag {
  LLVFD,
  LLVhlm
} GenLLVFDtag;

/* Parameters for the generation of a LLV waveform (in the form of a list of modes) */
typedef struct tagGenLLVFDparams {
  double tRef;               /* reference time (s) - GPS time at the frequency representing coalescence */
  double phiRef;             /* reference phase (rad) - phase at the frequency representing coalescence (or at fRef if specified) */
  double fRef;               /* reference frequency at which phiRef is set (Hz, default 0 which is interpreted as Mf=0.14) */
  double m1;                 /* mass of companion 1 (solar masses, default 2e6) */
  double m2;                 /* mass of companion 2 (solar masses, default 1e6) */
  double distance;           /* distance of source (Mpc, default 1e3) */
  double inclination;        /* inclination of source (rad, default pi/3) */
  double ra;                 /* first angle for the position in the sky (rad, default 0) */
  double dec;                /* second angle for the position in the sky (rad, default 0) */
  double polarization;       /* polarization angle (rad, default 0) */

  int nbmode;                /* number of modes to generate (starting with 22) - defaults to 5 (all modes) */
  double minf;               /* Minimal frequency, ignore if 0 (Hz, default=0) - will use first frequency covered by the ROM if higher */
  double maxf;               /* Maximal frequency, ignore if 0 (Hz, default=0) - will use last frequency covered by the ROM if lower */
  double deltaf;             /* When generating frequency series from the mode contributions, deltaf for the output (0 to set automatically at 1/2*1/(2T)) */
  int setphiRefatfRef;       /* Flag for adjusting the FD phase at phiRef at the given fRef, which depends also on tRef - if false, treat phiRef simply as an orbital phase shift (minus an observer phase shift) (default=1) */
  int tagnetwork;            /* Tag selecting the desired output format */
  int taggenwave;            /* Tag selecting the desired output format */
  int fromLLVtdfile;         /* Tag for loading time series for LLV detectors and FFTing */
  int nsamplesinfile;        /* Number of lines of input file */
  int binaryin;              /* Tag for loading the data in gsl binary form instead of text (default false) */
  int binaryout;             /* Tag for outputting the data in gsl binary form instead of text (default false) */
  char indir[256];           /* Path for the input directory */
  char infile[256];          /* Path for the input file */
  char outdir[256];          /* Path for the output directory */
  char outfileprefix[256];   /* Path for the output file */
} GenLLVFDparams;


#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _GENERATELLVFD_H */
