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


/************** GSL error handling and I/O ********************/

/* GSL error handler */
void Err_Handler(const char *reason, const char *file, int line, int gsl_errno) {
  printf("gsl: %s:%d: %s - %d\n", file, line, reason, gsl_errno);
}

/* Functions to read data from files */
int Read_Vector(const char dir[], const char fname[], gsl_vector *v) {
  char *path=malloc(strlen(dir)+64);

  sprintf(path,"%s/%s", dir, fname);
  FILE *f = fopen(path, "rb");
  if (!f) {
      return(FAILURE);
  }
  int ret = gsl_vector_fread(f, v);
  if (ret != 0) {
      fprintf(stderr, "Error reading data from %s.\n",path);
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
      return(FAILURE);
  }
  int ret = gsl_matrix_fread(f, m);
  if (ret != 0) {
      fprintf(stderr, "Error reading data from %s.\n", path);
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

/***************** Other structure functions ****************/
/* Function reproducing XLALSpinWeightedSphericalHarmonic
 * - Currently only supports s=-2, l=2,3,4,5 modes */
double complex SpinWeightedSphericalHarmonic(double theta, double phi, int s, int l, int m)
{
  static const char *func = "SpinWeightedSphericalHarmonic";
  double fac;
  double complex ans;
 
  /* sanity checks ... */
  if ( l < abs(s) ) {
    printf("Error - %s: Invalid mode s=%d, l=%d, m=%d - require |s| <= l\n", func, s, l, m );
    exit(1);
  }
  if ( l < abs(m) ) {
    printf("Error - %s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", func, s, l, m );
    exit(1);
  }
 
  if ( s == -2 ) {
    if ( l == 2 ) {
      switch ( m ) {
      case -2:
	fac = sqrt( 5.0 / ( 64.0 * PI ) ) * ( 1.0 - cos( theta ))*( 1.0 - cos( theta ));
	break;
      case -1:
	fac = sqrt( 5.0 / ( 16.0 * PI ) ) * sin( theta )*( 1.0 - cos( theta ));
	break;
         
      case 0:
	fac = sqrt( 15.0 / ( 32.0 * PI ) ) * sin( theta )*sin( theta );
	break;
         
      case 1:
	fac = sqrt( 5.0 / ( 16.0 * PI ) ) * sin( theta )*( 1.0 + cos( theta ));
	break;
         
      case 2:
	fac = sqrt( 5.0 / ( 64.0 * PI ) ) * ( 1.0 + cos( theta ))*( 1.0 + cos( theta ));
	break;     
      default:
	printf("Error - %s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", func, s, l, m );
	exit(1);
	break;
      } /*  switch (m) */
    }  /* l==2*/
    else if ( l == 3 ) {
      switch ( m ) {
      case -3:
	fac = sqrt(21.0/(2.0*PI))*cos(theta/2.0)*pow(sin(theta/2.0),5.0);
	break;
      case -2:
	fac = sqrt(7.0/4.0*PI)*(2.0 + 3.0*cos(theta))*pow(sin(theta/2.0),4.0);
	break;
      case -1:
	fac = sqrt(35.0/(2.0*PI))*(sin(theta) + 4.0*sin(2.0*theta) - 3.0*sin(3.0*theta))/32.0;
	break;
      case 0:
	fac = (sqrt(105.0/(2.0*PI))*cos(theta)*pow(sin(theta),2.0))/4.0;
	break;
      case 1:
	fac = -sqrt(35.0/(2.0*PI))*(sin(theta) - 4.0*sin(2.0*theta) - 3.0*sin(3.0*theta))/32.0;
	break;
           
      case 2:
	fac = sqrt(7.0/PI)*pow(cos(theta/2.0),4.0)*(-2.0 + 3.0*cos(theta))/2.0;
	break;     
           
      case 3:
	fac = -sqrt(21.0/(2.0*PI))*pow(cos(theta/2.0),5.0)*sin(theta/2.0);
	break;     
           
      default:
	printf("Error - %s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", func, s, l, m );
	exit(1);
	break;
      } 
    }   /* l==3 */ 
    else if ( l == 4 ) {
      switch ( m ) {
      case -4:
	fac = 3.0*sqrt(7.0/PI)*pow(cos(theta/2.0),2.0)*pow(sin(theta/2.0),6.0);
	break;
      case -3:
	fac = 3.0*sqrt(7.0/(2.0*PI))*cos(theta/2.0)*(1.0 + 2.0*cos(theta))*pow(sin(theta/2.0),5.0);
	break;
         
      case -2:
	fac = (3.0*(9.0 + 14.0*cos(theta) + 7.0*cos(2.0*theta))*pow(sin(theta/2.0),4.0))/(4.0*sqrt(PI));
	break;
      case -1:
	fac = (3.0*(3.0*sin(theta) + 2.0*sin(2.0*theta) + 7.0*sin(3.0*theta) - 7.0*sin(4.0*theta)))/(32.0*sqrt(2.0*PI));
	break;
      case 0:
	fac = (3.0*sqrt(5.0/(2.0*PI))*(5.0 + 7.0*cos(2.0*theta))*pow(sin(theta),2.0))/16.0;
	break;
      case 1:
	fac = (3.0*(3.0*sin(theta) - 2.0*sin(2.0*theta) + 7.0*sin(3.0*theta) + 7.0*sin(4.0*theta)))/(32.0*sqrt(2.0*PI));
	break;
      case 2:
	fac = (3.0*pow(cos(theta/2.0),4.0)*(9.0 - 14.0*cos(theta) + 7.0*cos(2.0*theta)))/(4.0*sqrt(PI));
	break;     
      case 3:
	fac = -3.0*sqrt(7.0/(2.0*PI))*pow(cos(theta/2.0),5.0)*(-1.0 + 2.0*cos(theta))*sin(theta/2.0);
	break;     
      case 4:
	fac = 3.0*sqrt(7.0/PI)*pow(cos(theta/2.0),6.0)*pow(sin(theta/2.0),2.0);
	break;     
      default:
	printf("Error - %s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", func, s, l, m );
	exit(1);
	break;
      }
    }    /* l==4 */
    else if ( l == 5 ) {
      switch ( m ) {
      case -5:
	fac = sqrt(330.0/PI)*pow(cos(theta/2.0),3.0)*pow(sin(theta/2.0),7.0);
	break;
      case -4:
	fac = sqrt(33.0/PI)*pow(cos(theta/2.0),2.0)*(2.0 + 5.0*cos(theta))*pow(sin(theta/2.0),6.0);
	break;
      case -3:
	fac = (sqrt(33.0/(2.0*PI))*cos(theta/2.0)*(17.0 + 24.0*cos(theta) + 15.0*cos(2.0*theta))*pow(sin(theta/2.0),5.0))/4.0;
	break;
      case -2:
	fac = (sqrt(11.0/PI)*(32.0 + 57.0*cos(theta) + 36.0*cos(2.0*theta) + 15.0*cos(3.0*theta))*pow(sin(theta/2.0),4.0))/8.0;
	break;
      case -1:
	fac = (sqrt(77.0/PI)*(2.0*sin(theta) + 8.0*sin(2.0*theta) + 3.0*sin(3.0*theta) + 12.0*sin(4.0*theta) - 15.0*sin(5.0*theta)))/256.0;
	break;
      case 0:
	fac = (sqrt(1155.0/(2.0*PI))*(5.0*cos(theta) + 3.0*cos(3.0*theta))*pow(sin(theta),2.0))/32.0;
	break;
      case 1:
	fac = sqrt(77.0/PI)*(-2.0*sin(theta) + 8.0*sin(2.0*theta) - 3.0*sin(3.0*theta) + 12.0*sin(4.0*theta) + 15.0*sin(5.0*theta))/256.0;
	break;
      case 2:
	fac = sqrt(11.0/PI)*pow(cos(theta/2.0),4.0)*(-32.0 + 57.0*cos(theta) - 36.0*cos(2.0*theta) + 15.0*cos(3.0*theta))/8.0;
	break;     
      case 3:
	fac = -sqrt(33.0/(2.0*PI))*pow(cos(theta/2.0),5.0)*(17.0 - 24.0*cos(theta) + 15.0*cos(2.0*theta))*sin(theta/2.0)/4.0;
	break;     
      case 4:
	fac = sqrt(33.0/PI)*pow(cos(theta/2.0),6.0)*(-2.0 + 5.0*cos(theta))*pow(sin(theta/2.0),2.0);
	break;     
      case 5:
	fac = -sqrt(330.0/PI)*pow(cos(theta/2.0),7.0)*pow(sin(theta/2.0),3.0);
	break;     
      default:
	printf("Error - %s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", func, s, l, m );
	exit(1);
	break;
      }
    }  /* l==5 */
    else {
      printf("Error - %s: Unsupported mode l=%d (only l in [2,5] implemented)\n", func, s);
      exit(1);
    }
  }
  else {
    printf("Error - %s: Unsupported mode s=%d (only s=-2 implemented)\n", func, s);
    exit(1);
  }
  if (m)
    ans = fac*cexp(I*m*phi);
  else
    ans = fac;
  return ans;
}
