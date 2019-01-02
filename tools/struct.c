/**
 * \author Sylvain Marsat, University of Maryland - NASA GSFC
 *
 * \brief C code for structures representing a waveform as a list of modes in amplitude/phase form.
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
#include "struct.h"


/* Call a function that references an object */
double ObjectFunctionCall(const ObjectFunction* this,double f){return this->function(this->object,f);};

/**************************************************************/
/* Functions computing the max and min between two int */
int max (int a, int b) { return a > b ? a : b; }
int min (int a, int b) { return a < b ? a : b; }

/************** GSL error handling and I/O ********************/

/* GSL error handler */
void Err_Handler(const char *reason, const char *file, int line, int gsl_errno) {
  printf("gsl: %s:%d: %s - %d\n", file, line, reason, gsl_errno);
  exit(1);
}

/* Functions to read binary data from files */
int Read_Vector(const char dir[], const char fname[], gsl_vector *v) {
  char *path=malloc(strlen(dir)+64);
  sprintf(path,"%s/%s", dir, fname);
  FILE *f = fopen(path, "rb");
  if (!f) {
    fprintf(stderr, "Error reading data from %s\n", path);
    free(path);
    return(FAILURE);
  }
  int ret = gsl_vector_fread(f, v);
  if (ret != 0) {
    fprintf(stderr, "Error reading data from %s.\n",path);
    free(path);
    return(FAILURE);
  }
  fclose(f);
  free(path);
  return(SUCCESS);
}
int Read_Matrix(const char dir[], const char fname[], gsl_matrix *m) {
  char *path=malloc(strlen(dir)+64);
  sprintf(path,"%s/%s", dir, fname);
  FILE *f = fopen(path, "rb");
  if (!f) {
    fprintf(stderr, "Error reading data from %s\n", path);
    free(path);
    return(FAILURE);
  }
  int ret = gsl_matrix_fread(f, m);
  if (ret != 0) {
    fprintf(stderr, "Error reading data from %s\n", path);
    free(path);
    return(FAILURE);
  }
  fclose(f);
  free(path);
  return(SUCCESS);
}
/* Functions to read text data from files */
int Read_Text_Vector(const char dir[], const char fname[], gsl_vector *v) {
  char *path=malloc(strlen(dir)+64);
  sprintf(path,"%s/%s", dir, fname);
  FILE *f = fopen(path, "rb");
  if (!f) {
    fprintf(stderr, "Error reading data from %s\n", path);
    free(path);
    return(FAILURE);
  }
  int ret = gsl_vector_fscanf(f, v);
  if (ret != 0) {
    fprintf(stderr, "Error reading data from %s.\n",path);
    free(path);
    return(FAILURE);
  }
  fclose(f);
  free(path);
  return(SUCCESS);
}
int Read_Text_Matrix(const char dir[], const char fname[], gsl_matrix *m) {
  char *path=malloc(strlen(dir)+64);
  sprintf(path,"%s/%s", dir, fname);
  FILE *f = fopen(path, "rb");
  if (!f) {
    fprintf(stderr, "Error reading data from %s\n", path);
    free(path);
    return(FAILURE);
  }
  int ret = gsl_matrix_fscanf(f, m);
  if (ret != 0) {
    fprintf(stderr, "Error reading data from %s.\n",path);
    free(path);
    return(FAILURE);
  }
  fclose(f);
  free(path);
  return(SUCCESS);
}
/* Functions to write data to files */
int Write_Vector(const char dir[], const char fname[], gsl_vector *v) {
  char *path=malloc(strlen(dir)+64);
  sprintf(path,"%s/%s", dir, fname);
  FILE *f = fopen(path, "w");
  if (!f) {
    fprintf(stderr, "Error writing data to %s\n", path);
    free(path);
    return(FAILURE);
  }
  int ret = gsl_vector_fwrite(f, v);
  if (ret != 0) {
    fprintf(stderr, "Error writing data to %s\n",path);
    free(path);
    return(FAILURE);
  }
  fclose(f);
  free(path);
  return(SUCCESS);
}
int Write_Matrix(const char dir[], const char fname[], gsl_matrix *m) {
  char *path=malloc(strlen(dir)+64);

  sprintf(path,"%s/%s", dir, fname);
  FILE *f = fopen(path, "w");
  if (!f) {
    fprintf(stderr, "Error writing data to %s\n", path);
    free(path);
    return(FAILURE);
  }
  int ret = gsl_matrix_fwrite(f, m);
  if (ret != 0) {
    fprintf(stderr, "Error writing data to %s\n", path);
    free(path);
    return(FAILURE);
  }
  fclose(f);
  free(path);
  return(SUCCESS);
}
/* Functions to write text data to files */
int Write_Text_Vector(const char dir[], const char fname[], gsl_vector *v) {
  char *path=malloc(strlen(dir)+64);
  sprintf(path,"%s/%s", dir, fname);
  FILE *f = fopen(path, "w");
  if (!f) {
    fprintf(stderr, "Error writing data to %s\n", path);
    free(path);
    return(FAILURE);
  }
  int ret = gsl_vector_fprintf(f, v, "%.16e");
  if (ret != 0) {
    fprintf(stderr, "Error writing data to %s\n",path);
    free(path);
    return(FAILURE);
  }
  fclose(f);
  free(path);
  return(SUCCESS);
}
int Write_Text_Matrix(const char dir[], const char fname[], gsl_matrix *m) {
  char *path=malloc(strlen(dir)+64);
  int ret = 0;

  sprintf(path,"%s/%s", dir, fname);
  FILE *f = fopen(path, "w");
  if (!f) {
    fprintf(stderr, "Error writing data to %s\n",path);
    free(path);
    return(FAILURE);
  }
  int N = (int) m->size1;
  int M = (int) m->size2;
  for(int i=0; i<N; i++){
    for(int j=0; j<M; j++){
      ret |= (fprintf(f, "%.16e ", gsl_matrix_get(m, i, j)) < 0);
    }
    if(i < N-1) ret |= (fprintf(f, "\n") < 0);
  }
  if (ret != 0) {
    fprintf(stderr, "Error writing data to %s\n",path);
    free(path);
    return(FAILURE);
  }
  fclose(f);
  free(path);
  return(SUCCESS);
}

/******** Functions to initialize and clean up CAmpPhaseFrequencySeries structure ********/
void CAmpPhaseFrequencySeries_Init(CAmpPhaseFrequencySeries **freqseries, const int n) {
  if(!freqseries) exit(1);
  /* Create storage for structures */
  if(!*freqseries) *freqseries=malloc(sizeof(CAmpPhaseFrequencySeries));
  else
  {
    CAmpPhaseFrequencySeries_Cleanup(*freqseries);
    *freqseries=malloc(sizeof(CAmpPhaseFrequencySeries));
  }
  gsl_set_error_handler(&Err_Handler);
  (*freqseries)->freq = gsl_vector_alloc(n);
  (*freqseries)->amp_real = gsl_vector_alloc(n);
  (*freqseries)->amp_imag = gsl_vector_alloc(n);
  (*freqseries)->phase = gsl_vector_alloc(n);
}
void CAmpPhaseFrequencySeries_Cleanup(CAmpPhaseFrequencySeries *freqseries) {
  if(freqseries->freq) gsl_vector_free(freqseries->freq);
  if(freqseries->amp_real) gsl_vector_free(freqseries->amp_real);
  if(freqseries->amp_imag) gsl_vector_free(freqseries->amp_imag);
  if(freqseries->phase) gsl_vector_free(freqseries->phase);
  free(freqseries);
}

/******** Functions to initialize and clean up CAmpPhaseSpline structure ********/
void CAmpPhaseSpline_Init(CAmpPhaseSpline **splines, const int n) {
  // if(splines != NULL) exit(1);
  /* Create storage for structures */

  /* FIXME!!! Dangerous thing, need "try-except" here otherwise doesn't work with python  */

  *splines=malloc(sizeof(CAmpPhaseSpline));
  // if(*splines != NULL){
  //   printf("step-0 \n");
  //   *splines=malloc(sizeof(CAmpPhaseSpline));
  //   printf("step1 \n");
  // }else
  // {
  //   CAmpPhaseSpline_Cleanup(*splines);
  // }
  gsl_set_error_handler(&Err_Handler);
  // printf("step1.0 \n");
  (*splines)->spline_amp_real = gsl_matrix_alloc(n, 5);
  // printf("step 2.1 \n");
  (*splines)->spline_amp_imag = gsl_matrix_alloc(n, 5);
  // printf("step 2.2 \n");
  (*splines)->quadspline_phase = gsl_matrix_alloc(n, 4);
  // printf("step 2.3 \n");
}
void CAmpPhaseSpline_Cleanup(CAmpPhaseSpline *splines) {
  if(splines->spline_amp_real) gsl_matrix_free(splines->spline_amp_real);
  if(splines->spline_amp_imag) gsl_matrix_free(splines->spline_amp_imag);
  if(splines->quadspline_phase) gsl_matrix_free(splines->quadspline_phase);
  free(splines);
}

/******** Functions to initialize and clean up CAmpPhaseGSLSpline structure ********/
void CAmpPhaseGSLSpline_Init(CAmpPhaseGSLSpline **splines, const int n) {
  if(!splines) exit(1);
  /* Create storage for structures */
  if(!*splines) *splines=malloc(sizeof(CAmpPhaseGSLSpline));
  else
  {
    CAmpPhaseGSLSpline_Cleanup(*splines);
  }
  gsl_set_error_handler(&Err_Handler);
  /* Note: for freq we won't copy the vector but simply copy the pointer to the existing one, so we don't allocate anything here */
  (*splines)->spline_amp_real = gsl_spline_alloc(gsl_interp_cspline, n);
  (*splines)->spline_amp_imag = gsl_spline_alloc(gsl_interp_cspline, n);
  (*splines)->spline_phase = gsl_spline_alloc(gsl_interp_cspline, n);
  (*splines)->accel_amp_real = gsl_interp_accel_alloc();
  (*splines)->accel_amp_imag = gsl_interp_accel_alloc();
  (*splines)->accel_phase = gsl_interp_accel_alloc();
}
void CAmpPhaseGSLSpline_Cleanup(CAmpPhaseGSLSpline *splines) {
  if(splines->spline_amp_real) gsl_spline_free(splines->spline_amp_real);
  if(splines->spline_amp_imag) gsl_spline_free(splines->spline_amp_imag);
  if(splines->spline_phase) gsl_spline_free(splines->spline_phase);
  if(splines->accel_amp_real) gsl_interp_accel_free(splines->accel_amp_real);
  if(splines->accel_amp_imag) gsl_interp_accel_free(splines->accel_amp_imag);
  if(splines->accel_phase) gsl_interp_accel_free(splines->accel_phase);
  free(splines);
}

/******** Functions to initialize and clean up ReImFrequencySeries structure ********/
void ReImFrequencySeries_Init(ReImFrequencySeries **freqseries, const int n) {
  if(!freqseries) exit(1);
  /* Create storage for structures */
  if(!*freqseries) *freqseries=malloc(sizeof(ReImFrequencySeries));
  else
  {
    ReImFrequencySeries_Cleanup(*freqseries);
  }
  gsl_set_error_handler(&Err_Handler);
  (*freqseries)->freq = gsl_vector_alloc(n);
  (*freqseries)->h_real = gsl_vector_alloc(n);
  (*freqseries)->h_imag = gsl_vector_alloc(n);
}
void ReImFrequencySeries_Cleanup(ReImFrequencySeries *freqseries) {
  if(freqseries->freq) gsl_vector_free(freqseries->freq);
  if(freqseries->h_real) gsl_vector_free(freqseries->h_real);
  if(freqseries->h_imag) gsl_vector_free(freqseries->h_imag);
  free(freqseries);
}

/******** Functions to initialize and clean up ReImTimeSeries structure ********/
void ReImTimeSeries_Init(ReImTimeSeries **timeseries, const int n) {
  if(!timeseries) exit(1);
  /* Create storage for structures */
  if(!*timeseries) *timeseries=malloc(sizeof(ReImTimeSeries));
  else
  {
    ReImTimeSeries_Cleanup(*timeseries);
  }
  gsl_set_error_handler(&Err_Handler);
  (*timeseries)->times = gsl_vector_alloc(n);
  (*timeseries)->h_real = gsl_vector_alloc(n);
  (*timeseries)->h_imag = gsl_vector_alloc(n);
}
void ReImTimeSeries_Cleanup(ReImTimeSeries *timeseries) {
  if(timeseries->times) gsl_vector_free(timeseries->times);
  if(timeseries->h_real) gsl_vector_free(timeseries->h_real);
  if(timeseries->h_imag) gsl_vector_free(timeseries->h_imag);
  free(timeseries);
}
/******** Functions to initialize and clean up AmpPhaseTimeSeries structure ********/
void AmpPhaseTimeSeries_Init(AmpPhaseTimeSeries **timeseries, const int n) {
  if(!timeseries) exit(1);
  /* Create storage for structures */
  if(!*timeseries) *timeseries=malloc(sizeof(AmpPhaseTimeSeries));
  else
  {
    AmpPhaseTimeSeries_Cleanup(*timeseries);
  }
  gsl_set_error_handler(&Err_Handler);
  (*timeseries)->times = gsl_vector_alloc(n);
  (*timeseries)->h_amp = gsl_vector_alloc(n);
  (*timeseries)->h_phase = gsl_vector_alloc(n);
}
void AmpPhaseTimeSeries_Cleanup(AmpPhaseTimeSeries *timeseries) {
  if(timeseries->times) gsl_vector_free(timeseries->times);
  if(timeseries->h_amp) gsl_vector_free(timeseries->h_amp);
  if(timeseries->h_phase) gsl_vector_free(timeseries->h_phase);
  free(timeseries);
}

/******** Functions to initialize and clean up RealTimeSeries structure ********/
void RealTimeSeries_Init(RealTimeSeries **timeseries, const int n) {
  if(!timeseries) exit(1);
  /* Create storage for structures */
  if(!*timeseries) *timeseries=malloc(sizeof(RealTimeSeries));
  else
  {
    RealTimeSeries_Cleanup(*timeseries);
  }
  gsl_set_error_handler(&Err_Handler);
  (*timeseries)->times = gsl_vector_alloc(n);
  (*timeseries)->h = gsl_vector_alloc(n);
}
void RealTimeSeries_Cleanup(RealTimeSeries *timeseries) {
  if(timeseries->times) gsl_vector_free(timeseries->times);
  if(timeseries->h) gsl_vector_free(timeseries->h);
  free(timeseries);
}

/***************** Functions for the ListmodesCAmpPhaseFrequencySeries structure ****************/
ListmodesCAmpPhaseFrequencySeries* ListmodesCAmpPhaseFrequencySeries_AddModeNoCopy(
	   ListmodesCAmpPhaseFrequencySeries* appended,  /* List structure to prepend to */
	   CAmpPhaseFrequencySeries* freqseries,  /* data to contain */
	   int l, /* major mode number */
	   int m  /* minor mode number */)
{
    ListmodesCAmpPhaseFrequencySeries* list;
    /* Check if the node with this mode already exists */
    list = appended;
    while( list ){
      if( l == list->l && m == list->m ){
	break;
      }
      list = list->next;
    }
    if( list ){ /* We don't allow for the case where the mode already exists in the list*/
      printf("Error: Tried to add an already existing mode to a ListmodesCAmpPhaseFrequencySeries ");
      return(NULL);
    } else { /* In that case, we do NOT COPY the input interpolated data, which therefore can't be
		used anywhere else; this will be acceptable as these operations will only be done
		when interpolating the initialization data */
      list = malloc( sizeof(ListmodesCAmpPhaseFrequencySeries) );
    }
    list->l = l;
    list->m = m;
    if( freqseries ){
      list->freqseries = freqseries;
    } else {
      list->freqseries = NULL;
    }
    if( appended ){
      list->next = appended;
    } else {
        list->next = NULL;
    }
    return list;
}
/* Get the element of a ListmodesCAmpPhaseFrequencySeries with a given index */
ListmodesCAmpPhaseFrequencySeries* ListmodesCAmpPhaseFrequencySeries_GetMode(
	   ListmodesCAmpPhaseFrequencySeries* const list,  /* List structure to get a particular mode from */
	   int l, /*< major mode number */
	   int m  /*< minor mode number */ )
{
    if( !list ) return NULL;

    ListmodesCAmpPhaseFrequencySeries *itr = list;
    while( itr->l != l || itr->m != m ){
        itr = itr->next;
        if( !itr ) return NULL;
    }
    return itr; /* The element returned is itself a pointer to a ListmodesCAmpPhaseFrequencySeries */
}
void ListmodesCAmpPhaseFrequencySeries_Destroy(
	   ListmodesCAmpPhaseFrequencySeries* list  /* List structure to destroy; notice that the data is destroyed too */
)
{
  ListmodesCAmpPhaseFrequencySeries* pop;
  while( (pop = list) ){
    if( pop->freqseries ){ /* Destroying the CAmpPhaseFrequencySeries data */
      CAmpPhaseFrequencySeries_Cleanup( pop->freqseries );
    }
    /* Notice that the mode indices l and m are not freed, like in SphHarmTimeSeries struct indices l and m */
    list = pop->next;
    free( pop );
  }
}

/***************** Functions for the ListmodesCAmpPhaseSpline structure ****************/
ListmodesCAmpPhaseSpline* ListmodesCAmpPhaseSpline_AddModeNoCopy(
	   ListmodesCAmpPhaseSpline* appended,  /* List structure to prepend to */
	   CAmpPhaseSpline* splines,  /* data to contain */
	   int l, /* major mode number */
	   int m  /* minor mode number */)
{
    ListmodesCAmpPhaseSpline* list;
    /* Check if the node with this mode already exists */
    list = appended;
    while( list ){
      if( l == list->l && m == list->m ){
	break;
      }
      list = list->next;
    }
    if( list ){ /* We don't allow for the case where the mode already exists in the list*/
      printf("Error: Tried to add an already existing mode to a ListmodesCAmpPhaseSpline ");
      return(NULL);
    } else { /* In that case, we do NOT COPY the input interpolated data, which therefore can't be
		used anywhere else; this will be acceptable as these operations will only be done
		when interpolating the initialization data */
      list = malloc( sizeof(ListmodesCAmpPhaseSpline) );
    }
    list->l = l;
    list->m = m;
    if( splines ){
      list->splines = splines;
    } else {
      list->splines = NULL;
    }
    if( appended ){
      list->next = appended;
    } else {
        list->next = NULL;
    }
    return list;
}
/* Get the element of a ListmodesCAmpPhaseSpline with a given index */
ListmodesCAmpPhaseSpline* ListmodesCAmpPhaseSpline_GetMode(
	   ListmodesCAmpPhaseSpline* const list,  /* List structure to get a particular mode from */
	   int l, /*< major mode number */
	   int m  /*< minor mode number */ )
{
    if( !list ) return NULL;

    ListmodesCAmpPhaseSpline *itr = list;
    while( itr->l != l || itr->m != m ){
        itr = itr->next;
        if( !itr ) return NULL;
    }
    return itr; /* The element returned is itself a pointer to a ListmodesCAmpPhaseSpline */
}
void ListmodesCAmpPhaseSpline_Destroy(
	   ListmodesCAmpPhaseSpline* list  /* List structure to destroy; notice that the data is destroyed too */
)
{
  ListmodesCAmpPhaseSpline* pop;
  while( (pop = list) ){
    if( pop->splines ){ /* Destroying the CAmpPhaseSpline data */
      CAmpPhaseSpline_Cleanup( pop->splines );
    }
    /* Notice that the mode indices l and m are not freed, like in SphHarmTimeSeries struct indices l and m */
    list = pop->next;
    free( pop );
  }
}

/***********************************************************************/
/**************** I/O functions for internal structures ****************/

/* Read waveform Real time series */
int Read_RealTimeSeries(RealTimeSeries** timeseries, const char dir[], const char file[], const int nblines, const int binary)
{
  /* Initalize and read input */
  int ret;
  gsl_matrix* inmatrix =  gsl_matrix_alloc(nblines, 2);
  if(!binary) ret = Read_Text_Matrix(dir, file, inmatrix);
  else ret = Read_Matrix(dir, file, inmatrix);

  /* Initialize structures */
  RealTimeSeries_Init(timeseries, nblines);

  /* Set values */
  gsl_vector_view timesview = gsl_matrix_column(inmatrix, 0);
  gsl_vector_view hview = gsl_matrix_column(inmatrix, 1);
  gsl_vector_memcpy((*timeseries)->times, &timesview.vector);
  gsl_vector_memcpy((*timeseries)->h, &hview.vector);

  /* Clean up */
  gsl_matrix_free(inmatrix);

  return ret;
}

/* Read waveform Amp/Phase time series */
int Read_AmpPhaseTimeSeries(AmpPhaseTimeSeries** timeseries, const char dir[], const char file[], const int nblines, const int binary)
{
  /* Initalize and read input */
  int ret;
  gsl_matrix* inmatrix =  gsl_matrix_alloc(nblines, 3);
  if(!binary) ret = Read_Text_Matrix(dir, file, inmatrix);
  else ret = Read_Matrix(dir, file, inmatrix);

  /* Initialize structures */
  AmpPhaseTimeSeries_Init(timeseries, nblines);

  /* Set values */
  gsl_vector_view timesview = gsl_matrix_column(inmatrix, 0);
  gsl_vector_view hampview = gsl_matrix_column(inmatrix, 1);
  gsl_vector_view hphaseview = gsl_matrix_column(inmatrix, 2);
  gsl_vector_memcpy((*timeseries)->times, &timesview.vector);
  gsl_vector_memcpy((*timeseries)->h_amp, &hampview.vector);
  gsl_vector_memcpy((*timeseries)->h_phase, &hphaseview.vector);

  /* Clean up */
  gsl_matrix_free(inmatrix);

  return ret;
}

/* Read waveform Re/Im time series */
int Read_ReImTimeSeries(ReImTimeSeries** timeseries, const char dir[], const char file[], const int nblines, const int binary)
{
  /* Initalize and read input */
  int ret;
  gsl_matrix* inmatrix =  gsl_matrix_alloc(nblines, 3);
  if(!binary) ret = Read_Text_Matrix(dir, file, inmatrix);
  else ret = Read_Matrix(dir, file, inmatrix);

  /* Initialize structures */
  ReImTimeSeries_Init(timeseries, nblines);

  /* Set values */
  gsl_vector_view timesview = gsl_matrix_column(inmatrix, 0);
  gsl_vector_view hrealview = gsl_matrix_column(inmatrix, 1);
  gsl_vector_view himagview = gsl_matrix_column(inmatrix, 2);
  gsl_vector_memcpy((*timeseries)->times, &timesview.vector);
  gsl_vector_memcpy((*timeseries)->h_real, &hrealview.vector);
  gsl_vector_memcpy((*timeseries)->h_imag, &himagview.vector);

  /* Clean up */
  gsl_matrix_free(inmatrix);

  return ret;
}

/* Output Re/Im frequency series */
int Write_ReImFrequencySeries(const char dir[], const char file[], ReImFrequencySeries* freqseries, const int binary)
{
  /* Initialize output */
  /* Note: assumes hplus, hcross have same length as expected */
  int nbfreq = freqseries->freq->size;
  gsl_matrix* outmatrix = gsl_matrix_alloc(nbfreq, 3);

  /* Set output matrix */
  gsl_matrix_set_col(outmatrix, 0, freqseries->freq);
  gsl_matrix_set_col(outmatrix, 1, freqseries->h_real);
  gsl_matrix_set_col(outmatrix, 2, freqseries->h_imag);

  /* Output */
  int ret;
  if (!binary) ret = Write_Text_Matrix(dir, file, outmatrix);
  else ret = Write_Matrix(dir, file, outmatrix);

  return ret;
}

/* Output real time series */
int Write_RealTimeSeries(const char dir[], const char file[], RealTimeSeries* timeseries, int binary)
{
  /* Initialize output */
  int nbtimes = timeseries->times->size;
  gsl_matrix* outmatrix = gsl_matrix_alloc(nbtimes, 2);

  /* Set data */
  gsl_matrix_set_col(outmatrix, 0, timeseries->times);
  gsl_matrix_set_col(outmatrix, 1, timeseries->h);

  /* Output */
  int ret;
  if(!binary) ret = Write_Text_Matrix(dir, file, outmatrix);
  else ret = Write_Matrix(dir, file, outmatrix);

  return ret;
}

/* Output Amp/Phase time series */
int Write_AmpPhaseTimeSeries(const char dir[], const char file[], AmpPhaseTimeSeries* timeseries, int binary)
{
  /* Initialize output */
  int nbtimes = timeseries->times->size;
  gsl_matrix* outmatrix = gsl_matrix_alloc(nbtimes, 3);

  /* Set data */
  gsl_matrix_set_col(outmatrix, 0, timeseries->times);
  gsl_matrix_set_col(outmatrix, 1, timeseries->h_amp);
  gsl_matrix_set_col(outmatrix, 2, timeseries->h_phase);

  /* Output */
  int ret;
  if(!binary) ret = Write_Text_Matrix(dir, file, outmatrix);
  else ret = Write_Matrix(dir, file, outmatrix);

  return ret;
}

/* Output Re/Im time series */
int Write_ReImTimeSeries(const char dir[], const char file[], ReImTimeSeries* timeseries, int binary)
{
  /* Initialize output */
  int nbtimes = timeseries->times->size;
  gsl_matrix* outmatrix = gsl_matrix_alloc(nbtimes, 3);

  /* Set data */
  gsl_matrix_set_col(outmatrix, 0, timeseries->times);
  gsl_matrix_set_col(outmatrix, 1, timeseries->h_real);
  gsl_matrix_set_col(outmatrix, 2, timeseries->h_imag);

  /* Output */
  int ret;
  if(!binary) ret = Write_Text_Matrix(dir, file, outmatrix);
  else ret = Write_Matrix(dir, file, outmatrix);

  return ret;
}
