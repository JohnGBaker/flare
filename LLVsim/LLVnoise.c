/**
 * \author Sylvain Marsat, University of Maryland - NASA GSFC
 *
 * \brief C code for the initialization of the instrumental noise for LIGO/VIRGO detectors.
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
#include "LLVnoise.h"


/******************************************************************************/
/****** Global variables storing interpolating splines for the noise PSD ******/

gsl_spline* __LLVSimFD_LHONoiseSpline_init = NULL; /* for initialization only */
gsl_spline** const __LLVSimFD_LHONoiseSpline = &__LLVSimFD_LHONoiseSpline_init;
gsl_spline* __LLVSimFD_LLONoiseSpline_init = NULL; /* for initialization only */
gsl_spline** const __LLVSimFD_LLONoiseSpline = &__LLVSimFD_LLONoiseSpline_init;
gsl_spline* __LLVSimFD_VIRGONoiseSpline_init = NULL; /* for initialization only */
gsl_spline** const __LLVSimFD_VIRGONoiseSpline = &__LLVSimFD_VIRGONoiseSpline_init;
gsl_interp_accel* __LLVSimFD_LHONoiseAccel_init = NULL; /* for initialization only */
gsl_interp_accel** const __LLVSimFD_LHONoiseAccel = &__LLVSimFD_LHONoiseAccel_init;
gsl_interp_accel* __LLVSimFD_LLONoiseAccel_init = NULL; /* for initialization only */
gsl_interp_accel** const __LLVSimFD_LLONoiseAccel = &__LLVSimFD_LLONoiseAccel_init;
gsl_interp_accel* __LLVSimFD_VIRGONoiseAccel_init = NULL; /* for initialization only */
gsl_interp_accel** const __LLVSimFD_VIRGONoiseAccel = &__LLVSimFD_VIRGONoiseAccel_init;
int __LLVSimFD_Noise_setup = FAILURE;

/* The number of points in the noise data files - required as the Read_Vector function needs a gsl vector already initialized to the right length */
#define noisedata_pts 3000

/**************************************************************/
/****** Functions loading and evaluating the noise PSD  *******/

/* Function parsing the environment variable $LLV_NOISE_DATA_PATH and trying to run LLVSimFD_Noise_Init in each */
void LLVSimFD_Noise_Init_ParsePath(void)
{
  if (!__LLVSimFD_Noise_setup) return(SUCCESS);

  int ret = FAILURE;
  char *envpath = NULL;
  char path[32768];
  char *brkt, *word;
  envpath = getenv("LLV_NOISE_DATA_PATH");
  if(!envpath) return(1);
  strncpy(path, envpath, sizeof(path));

  for(word=strtok_r(path,":",&brkt); word; word=strtok_r(NULL,":",&brkt))
  {
    ret = LLVSimFD_Noise_Init_File(word);
    if(ret == SUCCESS) break;
  }
  if(ret!=SUCCESS) {
    printf("Error: unable to find LLVSimFD noise data files in $LLV_NOISE_DATA_PATH\n");
    exit(FAILURE);
  }
  __LLVSimFD_Noise_setup = ret;
  return(ret);
}
/* Function loading the noise data from a directory */
int LLVSimFD_Noise_Init(const char dir[]) {
  if(!__LLVSimFD_Noise_setup) {
    printf("Error: LLVSimFD noise was already set up!");
    exit(1);
  }

  /* Loading noise data in gsl_vectors */
  int ret = SUCCESS;
  gsl_matrix* noise_LHO = gsl_matrix_alloc(noisedata_pts, 2);
  gsl_matrix* noise_LLO = gsl_matrix_alloc(noisedata_pts, 2);
  gsl_matrix* noise_VIRGO = gsl_matrix_alloc(noisedata_pts, 2);
  char* file_LIGO = malloc(strlen(dir)+64);
  char* file_VIRGO = malloc(strlen(dir)+64);
  sprintf(file_LIGO, "%s", "LIGO-P1200087-v18-aLIGO_DESIGN.txt");
  sprintf(file_VIRGO, "%s", "LIGO-P1200087-v18-AdV_DESIGN.txt");
  ret |= Read_Matrix(dir, file_LIGO, noise_LHO);
  ret |= Read_Matrix(dir, file_LIGO, noise_LLO);
  ret |= Read_Matrix(dir, file_VIRGO, noise_VIRGO);

  /* Linear interpolation of the data, and setting the gsl_spline structures */
  if(ret==SUCCESS) {
    gsl_vector* noise_LHO_freq = gsl_vector_alloc(noisedata_pts);
    gsl_vector* noise_LLO_freq = gsl_vector_alloc(noisedata_pts);
    gsl_vector* noise_VIRGO_freq = gsl_vector_alloc(noisedata_pts);
    gsl_vector* noise_LHO_data = gsl_vector_alloc(noisedata_pts);
    gsl_vector* noise_LLO_data = gsl_vector_alloc(noisedata_pts);
    gsl_vector* noise_VIRGO_data = gsl_vector_alloc(noisedata_pts);
    *__LLVSimFD_LHONoiseSpline = gsl_spline_alloc(gsl_interp_linear, noisedata_pts);
    *__LLVSimFD_LLONoiseSpline = gsl_spline_alloc(gsl_interp_linear, noisedata_pts);
    *__LLVSimFD_VIRGONoiseSpline = gsl_spline_alloc(gsl_interp_linear, noisedata_pts);
    gsl_spline_init(*__LLVSimFD_LHONoiseSpline, gsl_vector_const_ptr(noise_LHO_freq, 0), gsl_vector_const_ptr(noise_LHO_data, 0), noisedata_pts);
    gsl_spline_init(*__LLVSimFD_LLONoiseSpline, gsl_vector_const_ptr(noise_LLO_freq, 0), gsl_vector_const_ptr(noise_LLO_data, 0), noisedata_pts);
    gsl_spline_init(*__LLVSimFD_VIRGONoiseSpline, gsl_vector_const_ptr(noise_VIRGO_freq, 0), gsl_vector_const_ptr(noise_VIRGO_data, 0), noisedata_pts);
  }

  /* Cleaning and output */
  free(file_LIGO);
  free(file_VIRGO);
  return(ret);
}

/* The noise functions themselves */
double NoiseSnLHO(const double f) {
  if(__LLVSimFD_Noise_setup==FAILURE) {
    printf("Error: noise interpolation has not been set up\n");
    exit(1);
  }
  double sqrtSn = gsl_spline_eval(*__LLVSimFD_LHONoiseSpline, f, *__LLVSimFD_LHONoiseAccel);
  return sqrtSn * sqrtSn;
}
double NoiseSnLLO(const double f) {
  if(__LLVSimFD_Noise_setup==FAILURE) {
    printf("Error: noise interpolation has not been set up\n");
    exit(1);
  }
  double sqrtSn = gsl_spline_eval(*__LLVSimFD_LLONoiseSpline, f, *__LLVSimFD_LLONoiseAccel);
  return sqrtSn * sqrtSn;
}
double NoiseSnVIRGO(const double f) {
  if(__LLVSimFD_Noise_setup==FAILURE) {
    printf("Error: noise interpolation has not been set up\n");
    exit(1);
  }
  double sqrtSn = gsl_spline_eval(*__LLVSimFD_VIRGONoiseSpline, f, *__LLVSimFD_VIRGONoiseAccel);
  return sqrtSn * sqrtSn;
}
