/**
 * \author Sylvain Marsat, University of Maryland - NASA GSFC
 *
 * \brief C header for functions windowing and computing FFT/IFFT of time/frequency series.
 *
 */

#ifndef __GENERATEWAVEFORM_H__
#define __GENERATEWAVEFORM_H__ 1

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
#include "EOBNRv2HMROMstruct.h"
#include "EOBNRv2HMROM.h"
#include "fft.h"

/********************************** Structures ******************************************/

typedef enum GenWavetag {
  hlm,
  h22TD,
  hphcFD,
  hphcTD
} GenWavetag;

/* Parameters for the generation of a LISA waveform (in the form of a list of modes) */
typedef struct tagGenWaveParams {
  double tRef;               /* reference time (s) - GPS time at the frequency representing coalescence */
  double phiRef;             /* reference phase (rad) - phase at the frequency representing coalescence (or at fRef if specified) */
  double fRef;               /* reference frequency at which phiRef is set (Hz, default 0 which is interpreted as Mf=0.14) */
  double m1;                 /* mass of companion 1 (solar masses, default 2e6) */
  double m2;                 /* mass of companion 2 (solar masses, default 1e6) */
  double distance;           /* distance of source (Mpc, default 1e3) */
  double inclination;        /* inclination of source (rad, default pi/3) */
  double minf;               /* Minimal frequency, ignore if 0 (Hz, default=0) - will use first frequency covered by the ROM if higher */
  double maxf;               /* Maximal frequency, ignore if 0 (Hz, default=0) - will use last frequency covered by the ROM if lower */
  double deltatobs;          /* Observation duration (years, default=2) */
  int tagextpn;              /* Tag to allow PN extension of the waveform at low frequencies */
  double Mfmatch;            /* When PN extension allowed, geometric matching frequency: will use ROM above this value. If <=0, use ROM down to the lowest covered frequency */
  int setphiRefatfRef;       /* Flag for adjusting the FD phase at phiRef at the given fRef, which depends also on tRef - if false, treat phiRef simply as an orbital phase shift (minus an observer phase shift) (default=1) */
  int nbmode;                /* Number of modes to generate (starting with 22) - defaults to 1 (22 only) */
  int taggenwave;            /* Tag selecting the desired output format */
  double f1windowbeg;        /* If generating h22TD/hphcTD, start frequency for windowing at the beginning - set to 0 to ignore and use max(fstartobs, fLowROM, minf), where fLowROM is either the lowest frequency covered by the ROM or simply minf if PN extension is used (Hz, default=0) */
  double f2windowbeg;        /* If generating h22TD/hphcTD, stop frequency for windowing at the beginning - set to 0 to ignore and use 1.1*f1windowbeg (Hz, default=0) */
  double f1windowend;        /* If generating h22TD/hphcTD, start frequency for windowing at the end - set to 0 to ignore and use 0.995*f2windowend (Hz, default=0) */
  double f2windowend;        /* If generating h22TD/hphcTD, stop frequency for windowing at the end - set to 0 to ignore and use min(maxf, fHighROM), where fHighROM is the highest frequency covered by the ROM (Hz, default=0) */
  int tagh22fromfile;        /* Tag choosing wether to load h22 FD downsampled Amp/Phase from file (default 0) */
  int nsamplesinfile;        /* Number of lines of inputs file */
  int binaryin;              /* Tag for loading the data in gsl binary form instead of text (default false) */
  char indir[256];           /* Input directory */
  char infile[256];          /* Input file name */
  int binaryout;             /* Tag for outputting the data in gsl binary form instead of text (default 0) */
  char outdir[256];          /* Path for the output directory */
  char outfile[256];         /* Path for the output file */
} GenWaveParams;

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _GENERATEWAVEFORM_H */
