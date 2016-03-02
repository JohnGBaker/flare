/**
 * \author Sylvain Marsat, University of Maryland - NASA GSFC
 *
 * \brief C header for functions windowing and computing FFT/IFFT of time/frequency series.
 *
 */

#ifndef __GENERATETDITD_H__
#define __GENERATETDITD_H__ 1

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
#include "LISAgeometry.h"

/********************************** Structures ******************************************/

typedef enum GenWavetag {
  hphcTD,
  hphcFD,
  hlm
} GenWavetag;

/* Parameters for the generation of a LISA waveform (in the form of a list of modes) */
typedef struct tagGenTDITDparams {
  double lambda;             /* first angle for the position in the sky (rad, default 0) */
  double beta;               /* second angle for the position in the sky (rad, default 0) */
  double polarization;       /* polarization angle (rad, default 0) */
  double fLow;               /* Minimal frequency (Hz) - when set to 0 (default), use the first frequency covered by the ROM */
  int tagtdi;                /* Tag selecting the desired output format */
  int nlinesinfile;          /* Number of lines of input file */
  char indir[256];           /* Path for the input directory */
  char infile[256];          /* Path for the input file */
  char outdir[256];          /* Path for the output directory */
  char outfile[256];         /* Path for the output file */
} GenTDITDparams;


#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _GENERATETDITD_H */