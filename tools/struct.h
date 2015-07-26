/**
 * \author Sylvain Marsat, University of Maryland - NASA GSFC
 *
 * \brief C header for structures representing a waveform as a list of modes in amplitude/phase form.
 *
 */

#ifndef _STRUCT_H
#define _STRUCT_H

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


/***************************************************/
/*************** Type definitions ******************/

/* Complex frequency series in amplitude and phase representation (for one mode) */
typedef struct tagCAmpPhaseFrequencySeries
{
  gsl_vector* freq;
  gsl_vector* amp_real; /* We authorize complex amplitudes - will be used for the LISA response */
  gsl_vector* amp_imag; /* We authorize complex amplitudes - will be used for the LISA response */
  gsl_vector* phase;
} CAmpPhaseFrequencySeries;

/* Complex frequency series in real/imaginary part representation (for one mode, or their sum) */
typedef struct tagReImFrequencySeries
{
  gsl_vector* freq;
  gsl_vector* h_real;
  gsl_vector* h_imag;
} ReImFrequencySeries;

/* List structure, for a list of modes, each in amplitude and phase form */
typedef struct tagListmodesCAmpPhaseFrequencySeries
{
  CAmpPhaseFrequencySeries*                      freqseries; /* The frequencies series with amplitude and phase */
  int                                            l; /* Node mode l  */
  int                                            m; /* Node submode m  */
  struct tagListmodesCAmpPhaseFrequencySeries*    next; /* Next pointer */
} ListmodesCAmpPhaseFrequencySeries;

/**************************************************************/
/* Function computing the max between two int */
inline int max ( int a, int b ) { return a > b ? a : b; }

/**************************************************************/
/************** GSL error handling and I/O ********************/

/* GSL error handler */
void Err_Handler(const char *reason, const char *file, int line, int gsl_errno);

/* Functions to read/write data from files */
int Read_Vector(const char dir[], const char fname[], gsl_vector *v);
int Read_Matrix(const char dir[], const char fname[], gsl_matrix *m);
int Read_Text_Vector(const char dir[], const char fname[], gsl_vector *v);
int Read_Text_Matrix(const char dir[], const char fname[], gsl_matrix *m);
int Write_Text_Vector(const char dir[], const char fname[], gsl_vector *v);
int Write_Text_Matrix(const char dir[], const char fname[], gsl_matrix *m);

/**********************************************************/
/**************** Internal functions **********************/

/* Functions for list manipulations */
ListmodesCAmpPhaseFrequencySeries* ListmodesCAmpPhaseFrequencySeries_AddModeNoCopy(
	   ListmodesCAmpPhaseFrequencySeries* appended,  /* List structure to prepend to */
	   CAmpPhaseFrequencySeries* freqseries,  /* data to contain */
	   int l, /*< major mode number */
	   int m  /*< minor mode number */
);
ListmodesCAmpPhaseFrequencySeries* ListmodesCAmpPhaseFrequencySeries_GetMode( 
	   ListmodesCAmpPhaseFrequencySeries* const list,  /* List structure to get a particular mode from */
	   int l, /*< major mode number */
	   int m  /*< minor mode number */    
);
void ListmodesCAmpPhaseFrequencySeries_Destroy( 
	   ListmodesCAmpPhaseFrequencySeries* list  /* List structure to destroy; notice that the data is destroyed too */
);

/* Functions to initialize and clean up data structure */
void CAmpPhaseFrequencySeries_Init(
	 CAmpPhaseFrequencySeries **freqseries, /* double pointer for initialization */
	 const int n );                         /* length of the frequency series */
void CAmpPhaseFrequencySeries_Cleanup(CAmpPhaseFrequencySeries *freqseries);
void ReImFrequencySeries_Init(
	 ReImFrequencySeries **freqseries,      /* double pointer for initialization */
	 const int n );                         /* length of the frequency series */
void ReImFrequencySeries_Cleanup(ReImFrequencySeries *freqseries);

/***************** Functions to manipulate ReImFrequencySeries structure ****************/

void ReImFrequencySeries_AddCAmpPhaseFrequencySeries(
  struct tagReImFrequencySeries* freqseriesReIm,              /* Output Re/Im frequency series */
  struct tagCAmpPhaseFrequencySeries* freqseriesCAmpPhase,    /* Input CAmp/Phase frequency series, to be interpolated and added to the output */
  double fstartobsmode);                                      /* Starting frequency in case of limited duration of observations- assumed to have been scaled with the proper factor m/2 for this mode - set to 0 to ignore */
/* Function evaluating a ReImFrequencySeries by interpolating wach mode of a ListmodesCAmpPhaseFrequencySeries and summing them, given a set of frequencies */
void ReImFrequencySeries_SumListmodesCAmpPhaseFrequencySeries(
  struct tagReImFrequencySeries* freqseriesReIm,                    /* Output Re/Im frequency series - already initialized */
  struct tagListmodesCAmpPhaseFrequencySeries* listmodesCAmpPhase,  /* Input CAmp/Phase frequency series, to be interpolated */
  gsl_vector* freq,                                                 /* Input set of frequencies on which evaluating */
  double fstartobs);                                                /* Starting frequency in case of limited duration of observations - set to 0 to ignore */

/***************** Other structure functions ****************/
/* Additional function reproducing XLALSpinWeightedSphericalHarmonic */
double complex SpinWeightedSphericalHarmonic(double theta, double phi, int s, int l, int m); /* Currently only supports s=-2, l=2,3,4,5 modes */

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _STRUCT_H */
