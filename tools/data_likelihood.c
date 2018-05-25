/**
 * \author John Baker - NASA GSFC
 * \autho Sylvain Marsat - AIE
 *
 * \brief C code for the implementing likelihoods with Data provided on a fixed frequency domain grid.
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
#include "splinecoeffs.h"
#include "fresnel.h"
#include "likelihood.h"
#include "data_likelihood.h"

#include <time.h> /* for testing */

//Flags
int IntProdStyle=0;
int NratioCut = 4;

//Data whitening
//This function pre-multiplies the data with the inverse of the noise power

/* Function computing the integrand values */
void WhitenData(
  ReImUniformFrequencySeries** whitened_data, /* Output */
  ReImUniformFrequencySeries* data,           /* Input */
  ObjectFunction * Snoise,                    /* Noise function */
  double fLow,                                /* if fLow/fHigh > 0 restrict whitented data to withing that range */  
  double fHigh)
{
  gsl_set_error_handler(&Err_Handler);
  /* Determining the boundaries of indices */
  int imin=0;
  int imax=data->N-1;
  if( fLow > 0 ){
    while(Get_UniformFrequency(data,imin)<fLow){
      imin++;
      if(imin+1>imax){
	printf("Error: fLow is above the range of the data\n");
	exit(1);
      }
    }
  }
  if( fHigh > 0 ){
    while(Get_UniformFrequency(data,imax)>fHigh){
      imax--;
      if(imax-1<imin+1){
	printf("Error: fHigh is below the range of the (trimmed) data\n");
	exit(1);
      }
    }
  }
  
  /* Initializing output structure */
  int nbpts = imax - imin + 1;
  ReImUniformFrequencySeries_Init(whitened_data, nbpts);
  (*whitened_data)->fmin = Get_UniformFrequency(data,imin);
  (*whitened_data)->df = data->df;
  for(int i=imin; i<=imax; i++) {
    double f = Get_UniformFrequency(data,i);
    //printf("i=%i f=%g Snoise=%p\n",i,f,Snoise);
    double invSn = 1./ObjectFunctionCall(Snoise,f);
    gsl_vector_set((*whitened_data)->h_real, i-imin, gsl_vector_get(data->h_real,i)*invSn );
    gsl_vector_set((*whitened_data)->h_imag, i-imin, gsl_vector_get(data->h_imag,i)*invSn );
  }
};



// Functions for computing the integral 

//Perform explicit spline evaluation for inner-product computation of
//freqseries D with CAmpPhase spline set.  The grid of the D freqseries is summed for a simple integral approximant
double ReIntProdCAmpPhaseSpline(
  CAmpPhaseSpline* splines,                    //input "model"
  ReImUniformFrequencySeries* Dfreqseries,  //input defines ReIm uniform frequency series to be multiplied with spline evals and summed
  int imin,  // min/max  index-range of the freqseries to use include in the sum 
  int imax)
{
  int ispline=0;
  int nspline=splines->quadspline_phase->size1;
  double sum=0;
  double df = Dfreqseries->df;
  //printf("%i < i <%i \n",imin,imax);
  for(int i=imin;i<=imax;i++){
    double f=Get_UniformFrequency(Dfreqseries,i);
    double Dr=gsl_vector_get( Dfreqseries->h_real,i);
    double Di=gsl_vector_get( Dfreqseries->h_imag,i);

    //printf(" f:%g < %g < %g\n",gsl_matrix_get(splines->quadspline_phase, 0, 0),f,gsl_matrix_get(splines->quadspline_phase, splines->quadspline_phase->size1-1, 0));

    /* Adjust the index in the spline if necessary and compute */
    while(ispline<nspline-1&&gsl_matrix_get(splines->quadspline_phase, ispline+1, 0)<f){
      ispline++;
    }
    double eps = f - gsl_matrix_get(splines->quadspline_phase, ispline, 0);
    double eps2 = eps*eps;
    double eps3 = eps2*eps;
    gsl_vector_view coeffsampreal = gsl_matrix_row(splines->spline_amp_real, ispline);
    gsl_vector_view coeffsampimag = gsl_matrix_row(splines->spline_amp_imag, ispline);
    gsl_vector_view coeffsphase = gsl_matrix_row(splines->quadspline_phase, ispline);
    double Ar = EvalCubic(&coeffsampreal.vector, eps, eps2, eps3);
    double Ai = EvalCubic(&coeffsampimag.vector, eps, eps2, eps3);
    double Ph = EvalQuad(&coeffsphase.vector, eps, eps2);
    
    double c = cos(Ph); 
    double s = sin(Ph);
    double modelr = Ar*c-Ai*s;
    double modeli = Ar*s+Ai*c;
    double delta=df;
    if(i==imin||i+1==imax)
      delta=df/2;
    //sum += ( Dr*modelr+Di*modeli + I*(Dr*modeli-Di*modelr) ) * delta ; //We never need the image part and can save a few cycles
    double dsum = ( Dr*modelr+Di*modeli ) * delta ; 
    sum += dsum;
    //if(i>5000 &&i<5005)printf("AmPh %i %g %g\n",i,delta, ( Dr*modelr+Di*modeli ));
  }
  return sum;
};

//Perform explicit spline evaluation for inner-product computation of
//freqseries D with CAmpPhase spline set.  The grid of the D freqseries is summed for a simple integral approximant
double ReIntProdCAmpPhaseSplineB(
  CAmpPhaseSpline* splines,                    //input "model"
  ReImUniformFrequencySeries* Dfreqseries,  //input defines ReIm uniform frequency series to be multiplied with spline evals and summed
  int imin,  // min/max  index-range of the freqseries to use include in the sum 
  int imax)
{
  int ispline=0;
  int nspline=splines->quadspline_phase->size1;
  double sum=0;
  double df = Dfreqseries->df;
  double f0 = Dfreqseries->fmin;
  //printf("%i < i <%i \n",imin,imax);
  double fleft=gsl_matrix_get(splines->quadspline_phase, ispline, 0);
  double fright=gsl_matrix_get(splines->quadspline_phase, ispline+1, 0);
  int reset=1,direct=0;
  if((fright-fleft)/df > NratioCut )direct=1;
  complex double E,F,G,dG;
  
  gsl_vector_view coeffsampreal,coeffsampimag,coeffsphase;
    for(int i=imin;i<=imax;i++){
    double f=f0+i*df;
    double Dr=gsl_vector_get( Dfreqseries->h_real,i);
    double Di=gsl_vector_get( Dfreqseries->h_imag,i);

    /* Adjust the index in the spline if necessary and reset */
    while( ispline<nspline-2 && fright<f ){
      ispline++;
      reset=1;
      fleft=fright;
      //printf("ispline=%i size=%i\n",ispline,splines->quadspline_phase->size1);
      fright=gsl_matrix_get(splines->quadspline_phase, ispline+1, 0);
      if((fright-fleft)/df > NratioCut )direct=1;
      else direct=0;
    }
    if(reset){
      reset=0;
      coeffsampreal = gsl_matrix_row(splines->spline_amp_real, ispline);
      coeffsampimag = gsl_matrix_row(splines->spline_amp_imag, ispline);
      if(!direct)
	coeffsphase = gsl_matrix_row(splines->quadspline_phase, ispline);
      else {
	double eps0 = f - fleft;
	double p0 = gsl_matrix_get(splines->quadspline_phase, ispline,1);
	double p1 = gsl_matrix_get(splines->quadspline_phase, ispline,2);
	double p2 = gsl_matrix_get(splines->quadspline_phase, ispline,3);
	double phi0 = p0 + eps0 * ( p1 + eps0 * p2 );
	E = cexp(I*phi0);
	F = cexp(I*(p1+2*p2*eps0)*df);
	G = cexp(I*p2*df*df);
	dG = G*G;
      }
    }
    double eps = f - fleft;
    double eps2 = eps*eps;
    double eps3 = eps2*eps;
    double Ar = EvalCubic(&coeffsampreal.vector, eps, eps2, eps3);
    double Ai = EvalCubic(&coeffsampimag.vector, eps, eps2, eps3);
    double c,s;
    if(!direct){
      double Ph = EvalQuad(&coeffsphase.vector, eps, eps2);
      c = cos(Ph); 
      s = sin(Ph);
    } else {
      c = creal(E);
      s = cimag(E);
      E *= F*G;
      G *= dG;
    }
    double modelr = Ar*c-Ai*s;
    double modeli = Ar*s+Ai*c;
    double delta=df;
    double prod = Dr*modelr+Di*modeli;
    
    if(i==imin||i+1==imax)
      delta=df/2;
    sum += prod * delta;
  }
  return sum;
};


/* Function to compute an unweighted inner product betwee a model freq series and a data freq series
   using the data series as a the gird for the integral. */
double FDSinglemodeDataOverlap(
  CAmpPhaseFrequencySeries* model,     
  ReImUniformFrequencySeries* datachan)
{
  gsl_set_error_handler(&Err_Handler);

  /* Determining the boundaries of indices */
  /* i runs over the spline elements, k runs over the data grid */
  gsl_vector* freq1 = model->freq;
  int imin1 = 0;
  int imax1 = freq1->size - 1;
  double* f1 = freq1->data;
  double modelfmin = f1[0];
  double modelfmax = f1[ freq1->size - 1 ];
  double datafmin = Get_UniformFrequency( datachan, 0 );
  double datafmax = Get_UniformFrequency( datachan, datachan->N - 1 );
  double df = datachan->df;
  double kmin=0,kmax=datachan->N-1;
  //set limits to reference the data just within model limits
  if(modelfmin>datafmin)
    kmin=(modelfmin-datafmin)/df+1;
  if(modelfmax<datafmax)
    kmax=(modelfmax-datafmin)/df;

  /* Interpolating to complete the integrand, and summing */
  CAmpPhaseSpline* modelspline = NULL;
  BuildSplineCoeffs(&modelspline, model);

  double overlap;
  if(IntProdStyle==1)
    overlap = 4.0*ReIntProdCAmpPhaseSplineB(modelspline,datachan,kmin,kmax);   
  else
    overlap = 4.0*ReIntProdCAmpPhaseSpline(modelspline,datachan,kmin,kmax);   

  CAmpPhaseSpline_Cleanup(modelspline);
  
  return overlap;
}

/* Function computing the overlap (data|model) between tabulared any uniformly space Fourier domain data and a waveform given as list of modes for each channel 1,2,3.  No handling yet for gap where the signal is in band before the beginning of the data */
double FDCAmpPhaseModelDataOverlap(
  ListmodesCAmpPhaseFrequencySeries *listmodelchan, /* Waveform channel channel, list of modes in amplitude/phase form */
  ReImUniformFrequencySeries *datachan)                  /* Data channel channel */
{
  double overlap = 0;

  /* Main loop over the modes - goes through all the modes present, the same for all three channels 1,2,3 */
  ListmodesCAmpPhaseFrequencySeries* listelementmodelchan = listmodelchan;
  while(listelementmodelchan) { 
    //printf("computing overlap for element with max f = %g\n",gsl_vector_get(listelementmodelchan->freqseries->freq,listelementmodelchan->freqseries->freq->size-1));
    double overlapmode = FDSinglemodeDataOverlap(listelementmodelchan->freqseries, datachan);
    overlap += overlapmode;

    listelementmodelchan = listelementmodelchan->next;
  }

  return overlap;
}

/* Function computing the overlap (data|model) between tabulared any uniformly space Fourier domain data and a waveform given as list of modes for each channel 1,2,3.  No handling yet for gap where the signal is in band before the beginning of the data */
double FDReImModelDataOverlap(
  ReImFrequencySeries *modelchan, /* Waveform model channel, ReIm form */
  ReImUniformFrequencySeries *datachan)                  /* Data channel */
{
  double overlap = 0;

  int N=datachan->N;
  double df=datachan->df;
  printf("N=%i, df=%g\n",N,df);
  /* We assume comensurate freq grids */
  for(int i=0;i<N;i++){
    //test freqs
    if(fabs(gsl_vector_get(modelchan->freq,i)-Get_UniformFrequency(datachan,i))>1e-6)
      printf("incompatible freqs at i=%i: %g %g\n",i,gsl_vector_get(modelchan->freq,i),Get_UniformFrequency(datachan,i));
    double delta=df;
    if(i==0||i+1==N)delta=df/2.0;
    double dsnr2=0;
    dsnr2 += gsl_vector_get(modelchan->h_real,i)*gsl_vector_get(datachan->h_real,i);
    dsnr2 += gsl_vector_get(modelchan->h_imag,i)*gsl_vector_get(datachan->h_imag,i);
    overlap += 4.0*delta*dsnr2;
  }
  return overlap;
}





