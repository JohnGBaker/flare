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


#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/***************************************************/
/*************** Type definitions ******************/

/* Type for real functions */
typedef double (*RealFunctionPtr)(double);
typedef double (*RealObjectFunctionPtr)(const void *, double);

/* Type for real functions that reference an object */
typedef struct tagObjectFunction
{
  const void * object;
  RealObjectFunctionPtr function;
} ObjectFunction ;
double ObjectFunctionCall(const ObjectFunction*,double);

/* Signal Framing structure */
typedef struct {
  double fmin; //(Hz) minimum relevant freq. (for efficiency)
  double fmax; //(Hz) maximum relevant freq. (for efficiency)
  double tmin; //(s) minimum time extent, relative to waveform reference time
  double tmax; //(s) maximum time extent, relative to waveform reference time
}SignalFraming;
  
/* Complex frequency series in amplitude and phase representation (for one mode) */
typedef struct tagCAmpPhaseFrequencySeries
{
  gsl_vector* freq;
  gsl_vector* amp_real; /* We authorize complex amplitudes - will be used for the LISA response */
  gsl_vector* amp_imag; /* We authorize complex amplitudes - will be used for the LISA response */
  gsl_vector* phase;
} CAmpPhaseFrequencySeries;
/* GSL splines for complex amplitude and phase representation (for one mode) */
typedef struct tagCAmpPhaseGSLSpline
{
  gsl_vector* freq;
  gsl_spline* spline_amp_real; /* We authorize complex amplitudes - will be used for the LISA response */
  gsl_spline* spline_amp_imag; /* We authorize complex amplitudes - will be used for the LISA response */
  gsl_spline* spline_phase;
  gsl_interp_accel* accel_amp_real;
  gsl_interp_accel* accel_amp_imag;
  gsl_interp_accel* accel_phase;
} CAmpPhaseGSLSpline;
/* Splines in matrix form for complex amplitude and phase representation (for one mode) */
typedef struct tagCAmpPhaseSpline
{
  gsl_matrix* spline_amp_real; /* We authorize complex amplitudes - will be used for the LISA response */
  gsl_matrix* spline_amp_imag; /* We authorize complex amplitudes - will be used for the LISA response */
  gsl_matrix* quadspline_phase;
} CAmpPhaseSpline;

/* Complex frequency series in real/imaginary part representation (for one mode, or their sum) */
typedef struct tagReImFrequencySeries
{
  gsl_vector* freq;
  gsl_vector* h_real;
  gsl_vector* h_imag;
} ReImFrequencySeries;

/* Uniform grid complex frequency series in real/imaginary part representation (for one mode, or their sum) */
typedef struct tagReImUniformFrequencySeries
{
  int N;
  double fmin;
  double df;
  gsl_vector* h_real;
  gsl_vector* h_imag;
} ReImUniformFrequencySeries;

/* Complex frequency series in real/imaginary part representation (for one mode, or their sum) */
/* NOTE: for now, exact duplicata of ReImFrequencySeries - differentiated for readability of the code */
typedef struct tagReImTimeSeries
{
  gsl_vector* times;
  gsl_vector* h_real;
  gsl_vector* h_imag;
} ReImTimeSeries;

/* Complex frequency series in amplitude/phase representation (representing one mode) */
/* NOTE: for now, exact duplicata of ReImFrequencySeries - differentiated for readability of the code */
typedef struct tagAmpPhaseTimeSeries
{
  gsl_vector* times;
  gsl_vector* h_amp;
  gsl_vector* h_phase;
} AmpPhaseTimeSeries;

/* Real time series */
/* NOTE: could change the h to something more general, like values - also used for TD 22 amplitude, for instance  */
typedef struct tagRealTimeSeries
{
  gsl_vector* times;
  gsl_vector* h;
} RealTimeSeries;

/* List structure, for a list of modes, each in amplitude and phase form */
typedef struct tagListmodesCAmpPhaseFrequencySeries
{
  CAmpPhaseFrequencySeries*                      freqseries; /* The frequencies series with amplitude and phase */
  int                                            l; /* Node mode l  */
  int                                            m; /* Node submode m  */
  struct tagListmodesCAmpPhaseFrequencySeries*    next; /* Next pointer */
} ListmodesCAmpPhaseFrequencySeries;

/* List structure, for a list of modes, each with interpolated splines in amplitude and phase form */
typedef struct tagListmodesCAmpPhaseSpline
{
  CAmpPhaseSpline*                       splines; /* The frequencies series with amplitude and phase */
  int                                    l;       /* Node mode l  */
  int                                    m;       /* Node submode m  */
  struct tagListmodesCAmpPhaseSpline*    next;    /* Next pointer */
} ListmodesCAmpPhaseSpline;

/**************************************************************/
/* Functions computing the max and min between two int */
int max (int a, int b);
int min (int a, int b);

/**************************************************************/
/************** GSL error handling and I/O ********************/

/* GSL error handler */
void Err_Handler(const char *reason, const char *file, int line, int gsl_errno);

/* Functions to read/write data from files */
int Read_Vector(const char dir[], const char fname[], gsl_vector *v);
int Read_Matrix(const char dir[], const char fname[], gsl_matrix *m);
int Read_Text_Vector(const char dir[], const char fname[], gsl_vector *v);
int Read_Text_Matrix(const char dir[], const char fname[], gsl_matrix *m);
int Write_Vector(const char dir[], const char fname[], gsl_vector *v);
int Write_Matrix(const char dir[], const char fname[], gsl_matrix *m);
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
ListmodesCAmpPhaseSpline* ListmodesCAmpPhaseSpline_AddModeNoCopy(
	   ListmodesCAmpPhaseSpline* appended,  /* List structure to prepend to */
	   CAmpPhaseSpline* freqseries,  /* data to contain */
	   int l, /*< major mode number */
	   int m  /*< minor mode number */
);
ListmodesCAmpPhaseSpline* ListmodesCAmpPhaseSpline_GetMode(
	   ListmodesCAmpPhaseSpline* const list,  /* List structure to get a particular mode from */
	   int l, /*< major mode number */
	   int m  /*< minor mode number */
);
void ListmodesCAmpPhaseSpline_Destroy(
	   ListmodesCAmpPhaseSpline* list  /* List structure to destroy; notice that the data is destroyed too */
);

/* Functions to initialize and clean up data structure */
void CAmpPhaseFrequencySeries_Init(
	 CAmpPhaseFrequencySeries** freqseries, /* double pointer for initialization */
	 const int n );                         /* length of the frequency series */
void CAmpPhaseFrequencySeries_Cleanup(CAmpPhaseFrequencySeries* freqseries);
void CAmpPhaseSpline_Init(
	 CAmpPhaseSpline** splines,             /* double pointer for initialization */
	 const int n );                         /* length of the frequency series setting the splines */
void CAmpPhaseSpline_Cleanup(CAmpPhaseSpline* splines);
void CAmpPhaseGSLSpline_Init(
	 CAmpPhaseGSLSpline** splines,          /* double pointer for initialization */
	 const int n );                         /* length of the frequency series setting the splines */
void CAmpPhaseGSLSpline_Cleanup(CAmpPhaseGSLSpline* splines);
void ReImFrequencySeries_Init(
	 ReImFrequencySeries** freqseries,      /* double pointer for initialization */
	 const int n );                         /* length of the frequency series */
void ReImFrequencySeries_Cleanup(ReImFrequencySeries* freqseries);
void ReImUniformFrequencySeries_Init(
	 ReImUniformFrequencySeries** freqseries,      /* double pointer for initialization */
	 const int n );                                /* length of the frequency series */
void ReImUniformFrequencySeries_Cleanup(ReImUniformFrequencySeries* freqseries);
ReImUniformFrequencySeries * ReImFrequencySeries_ConvertToUniform(ReImFrequencySeries *oldfreqseries, int extrap);
double  Get_UniformFrequency(const ReImUniformFrequencySeries* freqseries, const int index);
void ReImTimeSeries_Init(
	 ReImTimeSeries** timeseries,           /* double pointer for initialization */
	 const int n );                         /* length of the time series */
void ReImTimeSeries_Cleanup(ReImTimeSeries* timeseries);
void AmpPhaseTimeSeries_Init(
	 AmpPhaseTimeSeries** timeseries,       /* double pointer for initialization */
	 const int n );                         /* length of the time series */
void AmpPhaseTimeSeries_Cleanup(AmpPhaseTimeSeries* timeseries);
void RealTimeSeries_Init(
	 RealTimeSeries** timeseries,      /* double pointer for initialization */
	 const int n );                /* length of the time series */
void RealTimeSeries_Cleanup(RealTimeSeries* timeseries);

/***********************************************************************/
/**************** I/O functions for internal structures ****************/

/* Note: at the moment, requires external input for the number of lines in the data */
int Read_RealTimeSeries(RealTimeSeries** timeseries, const char dir[], const char file[], const int nblines, const int binary);
int Read_AmpPhaseTimeSeries(AmpPhaseTimeSeries** timeseries, const char dir[], const char file[], const int nblines, const int binary);
int Read_ReImTimeSeries(ReImTimeSeries** timeseries, const char dir[], const char file[], const int nblines, const int binary);
int Read_ReImFrequencySeries(ReImFrequencySeries** freqseries, const char dir[], const char file[], const int nblines, const int binary);
int Read_ReImUniformFrequencySeries(ReImUniformFrequencySeries** freqseries, const char dir[], const char file[], const int nblines, const int binary);
int Write_ReImFrequencySeries(const char dir[], const char file[], ReImFrequencySeries* freqseries, const int binary);
int Write_ReImUniformFrequencySeries(const char dir[], const char file[], ReImUniformFrequencySeries* freqseries, const int binary);
int Write_RealTimeSeries(const char dir[], const char file[], RealTimeSeries* timeseries, int binary);
int Write_AmpPhaseTimeSeries(const char dir[], const char file[], AmpPhaseTimeSeries* timeseries, int binary);
int Write_ReImTimeSeries(const char dir[], const char file[], ReImTimeSeries* timeseries, int binary);


#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _STRUCT_H */
