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
int half_edges = 0;

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
//freqseries D with CAmpPhase spline set.  The grid of the D freqseries is summed for a simple integral approximant
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
//freqseries D with CAmpPhase spline set.  The grid of the D freqseries is summed for a simple integral approximant
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
//freqseries Data with CAmpPhase spline set.  The grid of the D freqseries is summed for a simple integral approximant
//In this case we try to decimate the
double ReIntProdCAmpPhaseSplineE(
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
  double fend=gsl_matrix_get(splines->quadspline_phase, nspline-1, 0);
  int reset=1,direct=0;
  if((fright-fleft)/df > NratioCut )direct=1;
  int ndown=1; //downsample factor
  double ndowndeltaf=0,ndf=df;
  
  complex double E,F,G,dG;
  int print=0;
  gsl_vector_view coeffsampreal,coeffsampimag,coeffsphase;
  //for(int i=imin;i<=imax;i++){
  int count=0;
  for(int i=imin;i<=imax;i+=ndown){
    //double f=f0+i*df;
    //double f=f0+(i+(ndown-1)*0.5)*df;ndowndeltaf=(ndown-1)*0.5*df
    double f=f0+i*df+ndowndeltaf;
    /* Adjust the index in the spline if necessary as well as ndown and reset */
    //if(i>13214 && i<13220)printf("i=%i:%g<%g<%g, nddf=%g\n",i,fleft,f,fright,ndowndeltaf);
      
    while( ispline<nspline-2 && fright<f ){
      ispline++;
      reset=1;
      fleft=fright;
      //printf("ispline=%i size=%i\n",ispline,splines->quadspline_phase->size1);
      fright=gsl_matrix_get(splines->quadspline_phase, ispline+1, 0);
      //adjust ndown:
      double p1=gsl_matrix_get(splines->quadspline_phase, ispline,2);
      double epstol=0.01;
      int ndownmax=epstol/(df*fmax(-p1,6/(fright-fleft)));
      int nmax=(fright-fleft)/df;
      while(2*ndown<=ndownmax)ndown*=2;
      while(ndown*df>fend-fright && ndown>1)ndown/=2;
      printf("%i: i=%i, %g<f<%g, di=(%i)->%i, p1=%g, sum=%g\n",ispline,i,fleft,fright,ndownmax,ndown,nmax,p1,sum);
      //print=ndown;
      //ndown=1;
      ndowndeltaf=(ndown-1)*0.5*df;
      ndf=df*ndown;
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
	G = cexp(I*p2*ndf*ndf);
	dG = G*G;
	//printf("phi0,E,F,G=%g,%g+%gi,%g+%gi,%g+%gi\n",phi0,creal(E),cimag(E),creal(F),cimag(F),creal(G),cimag(G));
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
    //double Dr=gsl_vector_get( Dfreqseries->h_real,i);
    //double Di=gsl_vector_get( Dfreqseries->h_imag,i);
    double Dr=0,Di=0;
    for(int k=0;k<ndown;k++){//Could use pre-decimated data here for faster result
      Dr+=gsl_vector_get( Dfreqseries->h_real,i+k);
      Di+=gsl_vector_get( Dfreqseries->h_imag,i+k);
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
  ReImUniformFrequencySeries* datachan,
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
  //printf("--> range %g < f < %g\n",modelfmin,modelfmax);

  /* Determining the boundaries of indices */
  /* i runs over the spline elements, k runs over the data grid */
  double datafmin = Get_UniformFrequency( datachan, 0 );
  double datafmax = Get_UniformFrequency( datachan, datachan->N - 1 );
  //printf("data range %g < f < %g\n",datafmin,datafmax);
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
    if(half_edges)
      //overlap = 4.0*ReIntProdCAmpPhaseSplineB(modelspline,datachan,kmin,kmax);   
      overlap = 4.0*ReIntProdCAmpPhaseSplineD(modelspline,datachan,modelfmin,modelfmax);   
    else
      overlap = 4.0*ReIntProdCAmpPhaseSplineE(modelspline,datachan,modelfmin,modelfmax);   
  else
    //overlap = 4.0*ReIntProdCAmpPhaseSpline(modelspline,datachan,kmin,kmax);   
    overlap = 4.0*ReIntProdCAmpPhaseSplineC(modelspline,datachan,modelfmin,modelfmax);   

  CAmpPhaseSpline_Cleanup(modelspline);
  
  return overlap;
}

/* Function computing the overlap (data|model) between tabulared any uniformly space Fourier domain data and a waveform given as list of modes for each channel 1,2,3
 This version will also estimate the relevant band of the signal and use that information to accelerate the computation.*/
double FDCAmpPhaseModelDataOverlap(
  ListmodesCAmpPhaseFrequencySeries *listmodelchan, /* Waveform channel channel, list of modes in amplitude/phase form */
  ReImUniformFrequencySeries *datachan,             /* Data channel channel */
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
    double overlapmode = FDSinglemodeDataOverlap(listelementmodelchan->freqseries, datachan, bandcut, bandSnoise);
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






