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
int IntProdStyle=1;
int NratioCut = 4;
int half_edges = 0;
int in_development=0;
int use_logarithmic=0;

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
//freqseries D with CAmpPhase spline set.  The grid of the D freqseries is summed for a simple integral approximant.
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
  printf(" IntProd: %g < f < %g\n",Get_UniformFrequency(Dfreqseries,imin),Get_UniformFrequency(Dfreqseries,imax));
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
    if(half_edges){ //HACK: do we interpret the minimum freq as a "center" or and "edge"; affects above as well
      if(i==imin||i+1==imax)  
	delta=df/2;
    }
    //sum += ( Dr*modelr+Di*modeli + I*(Dr*modeli-Di*modelr) ) * delta ; //We never need the imaginary part and can save a few cycles
    double dsum = ( Dr*modelr+Di*modeli ) * delta ; 
    sum += dsum;
    //if(i>5000 &&i<5005)printf("AmPh %i %g %g\n",i,delta, ( Dr*modelr+Di*modeli ));
  }
  return sum;
};

//Perform explicit spline evaluation for inner-product computation of
//freqseries D with CAmpPhase spline set.  
//Differences of C from plain version (Doesn't include "B" features):
//  -Switched integration range specification from index-based to freq-based
double ReIntProdCAmpPhaseSplineC(
  CAmpPhaseSpline* splines,                    //input "model"
  ReImUniformFrequencySeries* Dfreqseries,  //input defines ReIm uniform frequency series to be multiplied with spline evals and summed
  double mfmin,  // min/max  index-range of the freqseries to use include in the sum 
  double mfmax)
{
  int ispline=0;
  int nspline=splines->quadspline_phase->size1;
  double sum=0;
  double df = Dfreqseries->df;
  double dfleft=0,dfright=0;
  int imin=0,imax=Dfreqseries->N-1;
  //set limits to reference the data just within model limits
  double datafmin = Get_UniformFrequency( Dfreqseries, imin );
  double datafmax = Get_UniformFrequency( Dfreqseries, imax );
  if(mfmin<datafmin-df/2)mfmin=datafmin-df/2;//restrict to data limits (cell centered)
  if(mfmax>datafmax+df/2)mfmax=datafmax+df/2;//restrict to data limits (cell centered)
  if(mfmin>datafmin)
    imin=(mfmin-datafmin)/df+1;
  dfleft=Get_UniformFrequency(Dfreqseries,imin)-mfmin; //Account for the edge
  if(mfmax<datafmax)
    imax=(mfmax-datafmin)/df;
  dfright=mfmax-Get_UniformFrequency(Dfreqseries,imax); //Account for the edge
  //printf("%i < i <%i \n",imin,imax);
  printf(" IntProd: %g < f < %g\n",Get_UniformFrequency(Dfreqseries,imin),Get_UniformFrequency(Dfreqseries,imax));
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
    if(half_edges){ //This is consistent with dfleft==dfright==0
      if(i==imin||i+1==imax)  
	delta=df/2;
    } else {
      if(i==imin){  //This is an order df fix for the edges 
	//printf("df/left/right: %g, %g, %g\n",df, dfleft, dfright);
	delta=df/2+dfleft;
      }
      else if(i+1==imax)  
	delta=df/2+dfright;
    }
    //sum += ( Dr*modelr+Di*modeli + I*(Dr*modeli-Di*modelr) ) * delta ; //We never need the imaginary part and can save a few cycles
    double dsum = ( Dr*modelr+Di*modeli ) * delta ; 
    sum += dsum;
    //if(i>5000 &&i<5005)printf("AmPh %i %g %g\n",i,delta, ( Dr*modelr+Di*modeli ));
  }
  return sum;
};


//Perform explicit spline evaluation for inner-product computation of
//freqseries D with CAmpPhase spline set.
//Differences in B vs A:
//  -Added fast computation of exp(i*phase) fac (when sufficient samples in spline interval)
//  -Also includes diagnostic for estimating the effective time-width of FD spline intervals
double ReIntProdCAmpPhaseSplineB(
  CAmpPhaseSpline* splines,                    //input "model"
  ReImUniformFrequencySeries* Dfreqseries,  //input defines ReIm uniform frequency series to be multiplied with spline evals and summed
  int imin,  // min/max  index-range of the freqseries to use include in the sum 
  int imax)
{
  int ispline=0;
  int nspline=splines->quadspline_phase->size1;
  double sum=0;
  double dTmin,dTmax;
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
    
    if(half_edges){ //HACK: do we interpret the minimum freq as a "center" or and "edge"; affects above as well
      if(i==imin||i+1==imax)
	delta=df/2;
    }
    sum += prod * delta;
    
    /*diagnostic about time*/
    if(i==imin||i==imax){
      double p1=gsl_matrix_get(splines->quadspline_phase, ispline,2);
      double p2=gsl_matrix_get(splines->quadspline_phase, ispline,3);
      double dT=-(p1+p2*eps)/2.0/M_PI;
      if(i==imin)dTmin=dT;
      else dTmax=dT;
    }
  }
  printf(" %g < dT < %g\n",dTmin, dTmax);
    
  return sum;
};

//Perform explicit spline evaluation for inner-product computation of
//freqseries D with CAmpPhase spline set. 
//Difference from plain version:
//   -Adds C features: freq-based limits
//   -Adds B features: fast exp(i*ph)
double ReIntProdCAmpPhaseSplineD(
  CAmpPhaseSpline* splines,                    //input "model"
  ReImUniformFrequencySeries* Dfreqseries,  //input defines ReIm uniform frequency series to be multiplied with spline evals and summed
  double mfmin,  // min/max  index-range of the freqseries to use include in the sum 
  double mfmax)
{
  int ispline=0;
  int nspline=splines->quadspline_phase->size1;
  double sum=0;
  double dTmin,dTmax;
  double df = Dfreqseries->df;
  double f0 = Dfreqseries->fmin;
  //printf("%i < i <%i \n",imin,imax);
  double dfleft=0,dfright=0;
  int imin=0,imax=Dfreqseries->N-1;
  //set limits to reference the data just within model limits
  double datafmin = Get_UniformFrequency( Dfreqseries, imin );
  double datafmax = Get_UniformFrequency( Dfreqseries, imax );
  printf("%g < dataf < %g  : %i <= i <= %i  ",datafmin,datafmax,imin,imax);
  if(mfmin<datafmin-df/2)mfmin=datafmin-df/2;//restrict to data limits (cell centered)
  if(mfmax>datafmax+df/2)mfmax=datafmax+df/2;//restrict to data limits (cell centered)
  if(mfmin>datafmin)
    imin=(mfmin-datafmin)/df+1;
  dfleft=Get_UniformFrequency(Dfreqseries,imin)-mfmin; //Account for the edge
  if(mfmax<datafmax)
    imax=(mfmax-datafmin)/df;
  dfright=mfmax-Get_UniformFrequency(Dfreqseries,imax); //Account for the edge
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
    
    if(0&&half_edges){ //This is consistent with dfleft==dfright==0
      if(i==imin||i+1==imax)
	delta=df/2;
    } else {
      if(i==imin){ //This is an order df fix for the edges 
	//printf("df/left/right: %g, %g, %g\n",df, dfleft, dfright);
	delta=df/2+dfleft;
      } else if(i+1==imax) {  
	delta=df/2+dfright;
      }
    }
    sum += prod * delta;
    
    /*diagnostic about time*/
    if(i==imin||i==imax){
      double p1=gsl_matrix_get(splines->quadspline_phase, ispline,2);
      double p2=gsl_matrix_get(splines->quadspline_phase, ispline,3);
      double dT=-(p1+p2*eps)/2.0/M_PI;
      if(i==imin)dTmin=dT;
      else dTmax=dT;
    }
  }
  printf(" %g < dT < %g\n",dTmin, dTmax);
    
  return sum;
};

//Perform explicit spline evaluation for inner-product computation of
//freqseries Data with CAmpPhase spline set. 
//Difference from D:
//   -Performs downsampling of data, before computing as appropriate
//   -Has option for being passed a stack of data where the downsampling has already been performed.
//   -Special handling of case where amplitude varies faster than phase.
double ReIntProdCAmpPhaseSplineE(
  CAmpPhaseSpline* splines,                 //input spline "model " to be multiplied with data and summed
  ReImUniformFrequencySeries** datastack ,  //input stack of nstack ReIm uniform frequency series progressively decimated x2
  int nstack,
  double mfmin,  // min/max  index-range of the freqseries to use include in the sum 
  double mfmax)
{
  ReImUniformFrequencySeries* Dfreqseries=datastack[0];
  int ispline=0;
  int nspline=splines->quadspline_phase->size1;
  double sum=0;
  double dTmin,dTmax;
  double df = Dfreqseries->df;
  double f0 = Dfreqseries->fmin;
  //printf("%i < i <%i \n",imin,imax);
  double dfleft=0,dfright=0;
  int imin=0,imax=Dfreqseries->N-1;
  //set limits to reference the data just within model limits
  double datafmin = Get_UniformFrequency( Dfreqseries, imin );
  double datafmax = Get_UniformFrequency( Dfreqseries, imax );

  //Verify stacked data
  //printf("nstack=%i\n",nstack);
  if(in_development){
    for(int i=1;i<nstack;i++){
      ReImUniformFrequencySeries* fs=datastack[i];
      double tol=1e-12;
      if(fabs(fs->df-pow(2,i)*df)>tol)printf("ReIntProdCAmpPhaseSplineE: datastack[%i].df=%g != %g\n",i,fs->fmin,f0/pow(2,i));
      if(fabs(fs->fmin-fs->df/2.0-f0+df/2.0)>tol)printf("ReIntProdCAmpPhaseSplineE: datastack[%i].fmin=%g != %g\n",i,fs->fmin,f0+fs->df/2.0-df/2.0);
      //else printf("datum %i OK\n",i); 
    }
  }
  
  //printf("%g < dataf < %g  : %i <= i <= %i  ",datafmin,datafmax,imin,imax);
  if(mfmin<datafmin-df/2)mfmin=datafmin-df/2;//restrict to data limits (cell centered)
  if(mfmax>datafmax+df/2)mfmax=datafmax+df/2;//restrict to data limits (cell centered)
  if(mfmin>datafmin)
    imin=(mfmin-datafmin)/df+1;
  dfleft=Get_UniformFrequency(Dfreqseries,imin)-mfmin; //Account for the edge
  if(mfmax<datafmax)
    imax=(mfmax-datafmin)/df;
  dfright=mfmax-Get_UniformFrequency(Dfreqseries,imax); //Account for the edge
  double fleft=gsl_matrix_get(splines->quadspline_phase, ispline, 0);
  double fright=gsl_matrix_get(splines->quadspline_phase, ispline+1, 0);
  double fend=gsl_matrix_get(splines->quadspline_phase, nspline-1, 0);
  int reset=1,direct=0;
  if((fright-fleft)/df > NratioCut )direct=1;
  ///As we step through the integration, we will average down the data as needed for an effective approximation
  int ndown=1; //downsample factor
  int nextndown=ndown;
  int indown=0;
  int iklev=0,kstep=1,ikmark=0;
  double ndowndeltaf=0; //An offset to f0 for the downsampled data
  double ndf=df;        //Effective step size
  double ndf2=df*df;      
  
  complex double E,F,G,dG;
  int print=0;
  gsl_vector_view coeffsampreal,coeffsampimag,coeffsphase;
  //for(int i=imin;i<=imax;i++){
  int count=0;
  for(int i=imin;i<=imax;i+=ndown){
    //compute the next sample freq value
    double f=f0+i*df+ndowndeltaf;

    if(ndown!=nextndown){
      printf("ndown!=nextndwon\n");
      ndown=nextndown;
      iklev=min(indown,nstack-1);
      kstep=1<<(indown-iklev);
      ikmark=i>>iklev;
      ndowndeltaf=(ndown-1)*0.5*df;
      ndf=df*ndown;
      ndf2=ndf*ndf;
      //printf("i,nextndown,ikmark,iklev,kstep,ndown: %i %i %i %i %i %i\n",i, nextndown,ikmark,iklev,kstep,ndown);
      if((fright-fleft)/ndf > NratioCut )direct=1;
      else direct=0;
      reset=1;
    }
    /* Adjust the index in the spline if necessary as well as ndown and reset */
    ///If this step takes us beyond the right edge of the current spline segment range then we need to
    ///shift to the next spline segment and thus also to "reset" the interpolant curve integration.
    ///In that case, we also recompute the appropriate sampling rate.
    while( ispline<nspline-2 && fright<f ){
      ispline++;
      reset=1;
      fleft=fright;
      fright=gsl_matrix_get(splines->quadspline_phase, ispline+1, 0);
      //adjust ndown:
      double p1=gsl_matrix_get(splines->quadspline_phase, ispline,2);
      double p2=gsl_matrix_get(splines->quadspline_phase, ispline,3);
      double epstol=0.002;
      ///We apply the following criterion for determinine the downsampling rate:
      ///We want: ndf*(-dphidf) < epstol and Dfspline/df < epstol/Nfac with Nfac=6
      ///Could also (or instead of latter cond.) constrain abs((ndf)^2(-d2phidf2)) < 2*epstol
      const int Nfac=1;
      int ndownmax=epstol/df/fmax(abs(-p1),sqrt(abs(p2)));
      //int ndownmax=fmin( epstol/df/fmax(abs(-p1),sqrt(abs(p2))) , (fright-fleft)/df/Nfac);
      //printf("ndownmax=max: %i, %i, %i\n",(int)(epstol/df/abs(p1)),(int)(epstol/df/sqrt(abs(p2))) , (int)((fright-fleft)/df/Nfac));
      int nmax=(fright-fleft)/df;
      int ndownmaxtmp= max(1,ndownmax);
      indown=0;
      while(ndownmaxtmp>>=1)indown++;
      ndown=1<<indown;
      while(ndown*df>fend-fright && ndown>1){ndown/=2;indown--;}
      if(!use_logarithmic&&nstack>1){
	int nextiklev=min(indown,nstack-1);
	int nextkstep=1<<(indown-nextiklev);
	nextndown=ndown;
	int fac=1<<nextiklev;
	//printf("fac,i%%fac:%i,%i\n",fac,i%fac);
	ndown=(fac-i%fac);
	if((1<<iklev)*kstep!=ndown){
	  int ndowntmp=ndown;
	  indown=0;
	  while(ndowntmp>>=1)indown++;
	  if(ndown==1<<indown){//ndown is a power of two
	    iklev=min(indown,nstack-1);
	    ikmark=i>>iklev;
	    kstep=1<<(indown-nextiklev);
	  } else {
	    iklev=0;
	    kstep=ndown;
	  }
	}
	if(i%fac == 0 ){
	  iklev=nextiklev;
	  ikmark=i>>iklev;
	  kstep=nextkstep;
	  ndown=nextndown;
	}
	//printf("nextndown,ikmark,iklev,kstep,ndown: %i %i %i %i %i\n",nextndown,ikmark,iklev,kstep,ndown);
      } else { 
	nextndown=ndown;
      }
      ndowndeltaf=(ndown-1)*0.5*df;
      ndf=df*ndown;
      ndf2=ndf*ndf;
      if((fright-fleft)/ndf > NratioCut )direct=1;
      else direct=0;
    }
    if(reset){
      reset=0;
      coeffsampreal = gsl_matrix_row(splines->spline_amp_real, ispline);
      coeffsampimag = gsl_matrix_row(splines->spline_amp_imag, ispline);
      if(!direct)
	coeffsphase = gsl_matrix_row(splines->quadspline_phase, ispline);
      else {
	f=f0+i*df+ndowndeltaf;//ndowndeltaf may have changed since entering the loop
	//printf("f:%g+%g->%g\n",f0+i*df,ndowndeltaf,f);
	double eps0 = f - fleft;
	double p0 = gsl_matrix_get(splines->quadspline_phase, ispline,1);
	double p1 = gsl_matrix_get(splines->quadspline_phase, ispline,2);
	double p2 = gsl_matrix_get(splines->quadspline_phase, ispline,3);
	double phi0 = p0 + eps0 * ( p1 + eps0 * p2 );
	//printf("p1,p2,eps0,ndf=%g,%g,%g,%g\n",p1,p2,eps0,ndf);
	E = cexp(I*phi0);
	F = cexp(I*(p1+2*p2*eps0)*ndf);
	G = cexp(I*p2*ndf2);
	dG = G*G;
	//printf("phi0,E,F,G=%g,%g+%gi,%g+%gi,%g+%gi\n",phi0,creal(E),cimag(E),creal(F),cimag(F),creal(G),cimag(G));
      }
    }
    double eps = f - fleft;  
    double eps2 = eps*eps;
    double eps3 = eps2*eps;
    double Ar = EvalCubic(&coeffsampreal.vector, eps, eps2, eps3);
    double Ai = EvalCubic(&coeffsampimag.vector, eps, eps2, eps3);
    if(true){
      //To handle cases where the amplitude varies faster than the phase
      //Compute the average ampl over df rather than just the central value
      //A(eps+deps) = A(eps) + A''(eps)/2*deps^2
      //int_{-df/2}^{+df/2}{A(eps+deps)} = 2*A(eps)*(df/2) + 2*A''(eps)(df/2)^3/6
      //                                 = A(eps)*df + A''(eps)*df^3/24
      //avg = A(eps) + A''(eps)*df^2/24
      //A''(eps) = 2*p2 + 6*p3*eps)
      //avg = A(eps) + (p2/12 +p3*eps/4)*df^2
      Ar+=(gsl_vector_get(&coeffsampreal.vector, 2)/12.0+gsl_vector_get(&coeffsampreal.vector, 3)/4.0)*ndf2;
      Ai+=(gsl_vector_get(&coeffsampimag.vector, 2)/12.0+gsl_vector_get(&coeffsampimag.vector, 3)/4.0)*ndf2;
    }      
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
    //double Dr=gsl_vector_get( Dfreqseries->h_real,i);
    //double Di=gsl_vector_get( Dfreqseries->h_imag,i);
    double Dr=0,Di=0;

    if(nstack==1 || (use_logarithmic && ndown<17) || iklev<2){
      //linear method
      for(int k=0;k<ndown;k++){//Could use pre-decimated data here for faster result
	Dr+=gsl_vector_get( Dfreqseries->h_real,i+k);
	Di+=gsl_vector_get( Dfreqseries->h_imag,i+k);
      }
    } else {
      if(use_logarithmic){
	//logarithmic method:
	//Consider an example
	// 70=1000110 .. 86=1010110 - 32 =00110 
	//................. 70..86 -> null            0/0   
	//._. . . . . . ._. 35..43 -> 35..36, 42..43  1/1
	//  .   .   .___.   18..21 -> 20..21,  null   1/1
	//  ._______.        9..10 ->  9..10 ,  null  0/0  
	//          .          ref is 80            100/101
	//
	unsigned int ileft=i,iright=i+ndown;
	unsigned int il=ileft,ir=iright;
	//compute iref;
	//printf("\nileft/right: %i %i\n",ileft,iright);
	int ikmax=0;
	for(unsigned int ik=0;ik<nstack;ik++){
	  //here we progressively knock off the least significant bits until the lower/upper limits differ by no more than 1
	  //(or a larger power of 2 cut-off) or until we run out of stack levels.
	  //printf("ik,il,ir: %i %i %i\n",ik,il,ir);
	  ikmax=ik;
	  unsigned int nextir=ir>>1;
	  unsigned int nextil=((il+1)>>1);
	  if( nextil >= nextir || ik==nstack-1 )break;
	  ir=nextir;
	  il=nextil;
	}
	//printf("ikmax=%i\n",ikmax);
	//First add up whatever we have at the base level ikmax
	int fac=1<<ikmax;
	//printf("ikmax,il,ir: %i %i %i\n",ikmax,il,ir);
	for(int k=il;k<ir;k++){
	  //printf(" +stack[%i][%i] = %g %g\n",ikmax,k,gsl_vector_get( datastack[ikmax]->h_real,k),gsl_vector_get( datastack[ikmax]->h_imag,k));
	  Dr+=gsl_vector_get( datastack[ikmax]->h_real,k);
	  Di+=gsl_vector_get( datastack[ikmax]->h_imag,k);
	}
	Dr*=fac;
	Di*=fac;
	//next, got through the finer levels on the left and right
	il<<=(ikmax);
	ir<<=(ikmax);
	//printf("ikmax,il,ir: %i %i %i\n",ikmax,il,ir);
	unsigned int ilbits=il-ileft, irbits=iright-ir;
	//printf("ilbits=%i irbits=%i\n",ilbits,irbits);
	fac=1;
	for(unsigned int ik=0;irbits|ilbits>0 && ik<ikmax;ik++){
	  //printf("ik,ilbit,irbit,fac: %i %i %i %i\n",ik,ilbits & 1, irbits & 1, fac);
	  if( ilbits & 1 ){
	    unsigned int ind=il-ilbits;
	    //printf(" +stack[%i][%i] on left\n",ik,ind);
	    //printf(" +stack[%i][%i] on left = %g %g\n",ik,ind,gsl_vector_get( datastack[ik]->h_real,ind),gsl_vector_get( datastack[ik]->h_imag,ind));
	    Dr+=gsl_vector_get( datastack[ik]->h_real,ind)*fac;
	    Di+=gsl_vector_get( datastack[ik]->h_imag,ind)*fac;
	  }
	  if( irbits & 1 ){
	    unsigned int ind=irbits+ir-1;
	    //printf(" +stack[%i][%i] on right\n",ik,ind);
	    //printf(" +stack[%i][%i] on right = %g %g\n",ik,ind,gsl_vector_get( datastack[ik]->h_real,ind),gsl_vector_get( datastack[ik]->h_imag,ind));
	    Dr+=gsl_vector_get( datastack[ik]->h_real,ind)*fac;
	    Di+=gsl_vector_get( datastack[ik]->h_imag,ind)*fac;
	  }
	  fac*=2;
	  il>>=1;
	  ir>>=1;
	  ilbits>>=1;
	  irbits>>=1;
	}
      } else { //not use_logarithmic and ndown=nextndown
	for(int ik=0;ik<kstep;ik++){
	  //printf(" +stack[%i][%i] = %g %g\n",iklev,ikmark+ik,gsl_vector_get( datastack[iklev]->h_real,ik),gsl_vector_get( datastack[iklev]->h_imag,ik));
	  Dr+=gsl_vector_get( datastack[iklev]->h_real,ikmark+ik);
	  Di+=gsl_vector_get( datastack[iklev]->h_imag,ikmark+ik);
	}
	ikmark+=kstep;
	int fac=1<<iklev;
	Dr*=fac;
	Di*=fac;
	//printf("ndown,kstep,fac, Dr,Di: %i, %i %i %g %g\n",ndown,kstep,fac ,Dr,Di);

      }
      if(0){
	//This is just temporart for checking (and overwriting the above results!!)
	printf("Dr,Di: %g %g\n",Dr,Di);
	double saveDr=Dr;
	Dr=0;Di=0;
	for(int k=0;k<ndown;k++){
	  Dr+=gsl_vector_get( Dfreqseries->h_real,i+k);
	  Di+=gsl_vector_get( Dfreqseries->h_imag,i+k);
	}
	printf("Dr,Di: %g %g\n",Dr,Di);
	printf("relerr = %g\n",1.0-saveDr/Dr);
	if(fabs(1.0-saveDr/Dr)>1e-7){
	  printf("disagree\n");
	  for(int k=0;k<ndown;k++){
	    for(int ij=0;ij<indown;ij++){
	      int j=1<<ij;
	      if((i+k)%j==0)printf(" %2i %6i %12.8g %12.6g", ij, (i+k)>>ij, gsl_vector_get( datastack[ij]->h_real,(i+k)>>ij),
				   gsl_vector_get( datastack[ij]->h_imag,(i+k)>>ij));
	    }
	    printf("\n");
	  }
	}
	
      }

      
    }
    
    double prod = Dr*modelr+Di*modeli;
    if(0&half_edges){ //This is consistent with dfleft==dfright==0
      if(i==imin||i+1==imax)
	delta=df/2;
    } else {
      if(i==imin){ //This is an order df fix for the edges 
	//printf("df/left/right: %g, %g, %g\n",df, dfleft, dfright);
	delta=df/2+dfleft;
      } else if(i+1==imax) {  
	delta=df/2+dfright;
      }
    }
    count++;
    sum += prod * delta;
    if(print>0){
      printf("i=%i\n",i);
      printf("D=%g+%gi\n",Dr,Di);
      printf("f,df,prod,dsum,sum = %g, %g, %g, %g, %g\n",f,ndf,prod,prod*delta,sum);
      print--;
    }
    
    /*diagnostic about time*/
    if(i==imin||i==imax){
      double p1=gsl_matrix_get(splines->quadspline_phase, ispline,2);
      double p2=gsl_matrix_get(splines->quadspline_phase, ispline,3);
      double dT=-(p1+p2*eps)/2.0/M_PI;
      if(i==imin)dTmin=dT;
      else dTmax=dT;
    }
  }
  printf("Ncount=%6i %g < dT < %g\n",count,dTmin, dTmax);
    
  return sum;
};

/* Function to compute an unweighted inner product betwee a model freq series and a data freq series
   using the data series as a the gird for the integral. */
double FDSinglemodeDataOverlap(
  CAmpPhaseFrequencySeries* model,     
  ReImUniformFrequencySeries** datastack,
  int nstack,
  double bandcut,                                   /*target SNR value that we'd like to assure the excluded region does not attain*/
  ObjectFunction * bandSnoise)                      /* Noise function to use for the band cut calculation*/
{
  gsl_set_error_handler(&Err_Handler);

  //First estimate the relevant band
  gsl_vector* mfreq = model->freq;
  double* mf = mfreq->data;
  double modelfmin = mf[0];
  double modelfmax = mf[ mfreq->size - 1 ];
  //printf("model range %g < f < %g\n",modelfmin,modelfmax);
  if(bandcut>0)FDSinglemodeEstimateBand( model, bandSnoise, bandcut, &modelfmin, &modelfmax);
  printf("%g < f < %g ",modelfmin,modelfmax);

  /* Determining the boundaries of indices */
  /* i runs over the spline elements, k runs over t`he data grid */
  double datafmin = Get_UniformFrequency( datastack[0], 0 );
  double datafmax = Get_UniformFrequency( datastack[0], datastack[0]->N - 1 );
  //printf("data range %g < f < %g\n",datafmin,datafmax);
  double df = datastack[0]->df;
  double kmin=0,kmax=datastack[0]->N-1;
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
    if(half_edges)
      //overlap = 4.0*ReIntProdCAmpPhaseSplineB(modelspline,(*datachan),kmin,kmax);   
      overlap = 4.0*ReIntProdCAmpPhaseSplineD(modelspline,datastack[0],modelfmin,modelfmax);   
    else
      overlap = 4.0*ReIntProdCAmpPhaseSplineE(modelspline,datastack,nstack,modelfmin,modelfmax);   
  else
    //overlap = 4.0*ReIntProdCAmpPhaseSpline(modelspline,datachan,kmin,kmax);   
    overlap = 4.0*ReIntProdCAmpPhaseSplineC(modelspline,datastack[0],modelfmin,modelfmax);   

  CAmpPhaseSpline_Cleanup(modelspline);
  
  return overlap;
}

/* Function computing the overlap (data|model) between tabulared any uniformly space Fourier domain data and a waveform given as list of modes for each channel 1,2,3
 This version will also estimate the relevant band of the signal and use that information to accelerate the computation.*/
double FDCAmpPhaseModelDataOverlap(
  ListmodesCAmpPhaseFrequencySeries *listmodelchan, /* Waveform channel channel, list of modes in amplitude/phase form */
  ReImUniformFrequencySeries **datastack,           /* Data channel channel */
  int nstack,                                       /* number of x2 decimated data series.*/ 
  double bandcut,                                   /*target SNR value that we'd like to assure the excluded region does not attain*/
  ObjectFunction * bandSnoise)                      /* Noise function to use for the band cut calculation*/
  {
  double overlap = 0;

  /* Main loop over the modes - goes through all the modes present, the same for all three channels 1,2,3 */
  ListmodesCAmpPhaseFrequencySeries* listelementmodelchan = listmodelchan;
  while(listelementmodelchan) { 
    printf("mode (%i,%i)\n",listelementmodelchan->l,listelementmodelchan->m);
    //Probably need this fstartobs implementation here, but fstartobs itself might need to be passed in
    //double fstartobs = Newtonianfoft(params->m1, params->m2, globalparams->deltatobs);
    //int mmax1 = max(2, listelementh1chan1->m);
    //double fcutLow = fmax(fLow, ((double) mmax1)/2. * fstartobs); 
    double overlapmode = FDSinglemodeDataOverlap(listelementmodelchan->freqseries, datastack, nstack, bandcut, bandSnoise);
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






