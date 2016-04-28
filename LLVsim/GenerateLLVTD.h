/**
 * \author Sylvain Marsat, University of Maryland - NASA GSFC
 *
 * \brief C header for functions windowing and computing FFT/IFFT of time/frequency series.
 *
 */

#ifndef __GENERATELLVTD_H__
#define __GENERATELLVTD_H__ 1

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
typedef struct tagGenLLVTDparams {
  double lambda;             /* first angle for the position in the sky (rad, default 0) */
  double beta;               /* second angle for the position in the sky (rad, default 0) */
  double polarization;       /* polarization angle (rad, default 0) */
  int tagnetwork;            /* Tag selecting the desired detector network */
  int nsamplesinfile;        /* Number of lines of input file */
  int binaryin;              /* Tag for loading the data in gsl binary form instead of text (default false) */
  int binaryout;             /* Tag for outputting the data in gsl binary form instead of text (default false) */
  char indir[256];           /* Path for the input directory */
  char infile[256];          /* Path for the input file */
  char outdir[256];          /* Path for the output directory */
  char outfile[256];         /* Path for the output file */
} GenLLVTDparams;


#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _GENERATELLVTD_H */
