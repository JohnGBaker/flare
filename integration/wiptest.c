//test of wip 
//Written by John G Baker at gsfc.nasa.gov
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <time.h>
#include "wip.h"

double testfn_I_p(double x, double c2,double c3,double x0){
  double xo=x-x0;
  return xo*xo*(c2+xo*c3);
};

double testfn_I_A(double x, double sig,double xc){
  return erf((x-xc)/sig);
};

double complex testfn_dI_A(double x,double c2,double c3,double sig,double x0,double xc){
  const double twoortpi=2.0/sqrt(M_PI);
  double xo=x-x0;
  double dp=xo*(2.0*c2+3*xo*c3);
  double xnorm=(x-xc)/sig;
  //double expfac=exp(-xnorm*xnorm);
  //double complex wfac=I*dp*Faddeeva_w(I*xnorm,1.0e-11)
  //double complex dI=I*dp + ( twoortpi/sig - wfac ) * expfac 
  double complex dI=twoortpi/sig*exp(-xnorm*xnorm) + I*dp*erf(xnorm);
  return dI;
};

double Swhitenoise(double x){
  return 1;
};

extern int wip_count;

int main (){
  //Define the test function
  double c2,c3,x0,xc,a0,sig;
  c2=-1.41;
  //c3=+0.0199;
  c3=+0.0203;
  c2*=5;
  c3*=5;
  xc=3.47;
  x0=-20;
  a0=10.0;
  sig=1.51e2;

  //Defin the f grid
  double fmax, fmin;
  fmax=12.0;
  fmin=1e-4;
  int N0=20,N=N0;

  //wip_uselogf=1

  //#numerical test function check
  if(1){
    srand(time(NULL));
    double x=rand()/((double) RAND_MAX);
    double dx=0.01;
    int twotolevel=1;
    double eps=dx;
    int level;
    for(level=0;level<10;level++){	
      eps=dx/twotolevel;
      double complex Ival = a0*cexp(I*testfn_I_p(x, c2, c3, x0))*testfn_I_A(x, sig, xc);
      double complex dIval = a0*cexp(I*testfn_I_p(x, c2, c3, x0))*testfn_dI_A(x, c2, c3, sig, x0, xc);
      double complex Ivalp = a0*cexp(I*testfn_I_p(x+eps, c2, c3, x0))*testfn_I_A(x+eps, sig, xc);
      double complex Ivalm = a0*cexp(I*testfn_I_p(x-eps, c2, c3, x0))*testfn_I_A(x-eps, sig, xc);
      //printf( "Ivalp=%g+%gi\n",creal(Ivalp),cimag(Ivalp));
      //printf( "Ivalm=%g+%gi\n",creal(Ivalm),cimag(Ivalm));
      //printf( "Ivalp-Ivalm=%g+%gi\n",creal(Ivalp-Ivalm),cimag(Ivalp-Ivalm));
      double complex dIest=(Ivalp-Ivalm)/(eps*2.0);
      //print "level ",level,": Ip=",Ivalp,"\t(",dIest,"-",dIval,"\t)=\t",dIest-dIval,"->\t",(dIest-dIval)*4**level 
      printf("level %i: Ip=%6.5g+%6.5gi\t( %6.5g+%6.5gi - %6.5g+%6.5gi )=\t%6.5g+%6.5gi \t-> %6.5g+%6.5gi\n",level,creal(Ivalp),cimag(Ivalp),creal(dIest),cimag(dIest),creal(dIval),cimag(dIval),creal(dIest-dIval),cimag(dIest-dIval),creal((dIest-dIval)*twotolevel*twotolevel),cimag((dIest-dIval)*twotolevel*twotolevel));
      twotolevel*=2.0;
    }
  }

  double complex rex=a0*(cexp(I*testfn_I_p(fmax, c2, c3, x0))*testfn_I_A(fmax, sig, xc)-cexp(I*testfn_I_p(fmin, c2, c3, x0))*testfn_I_A(fmin, sig, xc));
  rex=conj(rex);
  printf("rex=%g + %gi\n",creal(rex),cimag(rex));

  if(1){
    N=N0;
    printf("\nInterpQuad(vectored) IP (N->%i):\n",N);
    int i;
    double ff[N],s1ar[N],s2ar[N],s1ai[N],s2ai[N],s1p[N],s2p[N];
    for(i=0;i<N;i++){
      ff[i] = i/(N-1.0)*(fmax-fmin)+fmin;
      s1ar[i] = 1.0;
      s1ai[i] = 0.0;
      s1p[i] = 0.0;
      double complex s2a=a0*testfn_dI_A(ff[i], c2, c3, sig, x0, xc);
      s2ar[i] = creal(s2a);
      s2ai[i] = cimag(s2a);
      s2p[i] = testfn_I_p(ff[i], c2, c3, x0);
    }
    double start=((double)clock())/CLOCKS_PER_SEC;
    double complex rs0=wip_phase(ff,N,ff,N,s1ar,s1ai,s1p,s2ar,s2ai,s2p,Swhitenoise,1.0,-1,-1);  
    double end=((double)clock())/CLOCKS_PER_SEC;
    double complex err=rs0-rex;
    printf( "ran %i steps.\n",wip_count);
    printf( "result= %g + %gi  err= %g + %gi\n",creal(rs0),cimag(rs0),creal(err),cimag(err));
    printf( "total time= %g   time per eval= %g\n",end-start,(end-start)/wip_count);
    printf( "----------------------------------------------------\n");
  }
  
  if(1){
    N=N0*4;
    printf("\nInterpQuad(vectored) IP (N->%i):\n",N);
    int i;
    double ff[N],s1ar[N],s2ar[N],s1ai[N],s2ai[N],s1p[N],s2p[N];
    for(i=0;i<N;i++){
      ff[i] = i/(N-1.0)*(fmax-fmin)+fmin;
      s1ar[i] = 1.0;
      s1ai[i] = 0.0;
      s1p[i] = 0.0;
      double complex s2a=a0*testfn_dI_A(ff[i], c2, c3, sig, x0, xc);
      s2ar[i] = creal(s2a);
      s2ai[i] = cimag(s2a);
      s2p[i] = testfn_I_p(ff[i], c2, c3, x0);
    }
    double start=((double)clock())/CLOCKS_PER_SEC;
    double complex rs0=wip_phase(ff,N,ff,N,s1ar,s1ai,s1p,s2ar,s2ai,s2p,Swhitenoise,1.0,-1.0,-1.0);  
    double end=((double)clock())/CLOCKS_PER_SEC;
    double complex err=rs0-rex;
    printf( "ran %i steps.\n",wip_count);
    printf( "result= %g + %gi  err= %g + %gi\n",creal(rs0),cimag(rs0),creal(err),cimag(err));
    printf( "total time= %g   time per eval= %g\n",end-start,(end-start)/wip_count);
    printf( "----------------------------------------------------\n");
  }

  if(1){
    N=N0*16;
    printf("\nInterpQuad(vectored) IP (N->%i):\n",N);
    int i;
    double ff[N],s1ar[N],s2ar[N],s1ai[N],s2ai[N],s1p[N],s2p[N];
    for(i=0;i<N;i++){
      ff[i] = i/(N-1.0)*(fmax-fmin)+fmin;
      s1ar[i] = 1.0;
      s1ai[i] = 0.0;
      s1p[i] = 0.0;
      double complex s2a=a0*testfn_dI_A(ff[i], c2, c3, sig, x0, xc);
      s2ar[i] = creal(s2a);
      s2ai[i] = cimag(s2a);
      s2p[i] = testfn_I_p(ff[i], c2, c3, x0);
    }
    double start=((double)clock())/CLOCKS_PER_SEC;
    double complex rs0=wip_phase(ff,N,ff,N,s1ar,s1ai,s1p,s2ar,s2ai,s2p,Swhitenoise,1.0,-1.0,-1.0);  
    double end=((double)clock())/CLOCKS_PER_SEC;
    double complex err=rs0-rex;
    printf( "ran %i steps.\n",wip_count);
    printf( "result= %g + %gi  err= %g + %gi\n",creal(rs0),cimag(rs0),creal(err),cimag(err));
    printf( "total time= %g   time per eval= %g\n",end-start,(end-start)/wip_count);
    printf( "----------------------------------------------------\n");
  }

  if(1){
    N=N0*64;
    printf("\nInterpQuad(vectored) IP (N->%i):\n",N);
    int i;
    double ff[N],s1ar[N],s2ar[N],s1ai[N],s2ai[N],s1p[N],s2p[N];
    for(i=0;i<N;i++){
      ff[i] = i/(N-1.0)*(fmax-fmin)+fmin;
      s1ar[i] = 1.0;
      s1ai[i] = 0.0;
      s1p[i] = 0.0;
      double complex s2a=a0*testfn_dI_A(ff[i], c2, c3, sig, x0, xc);
      s2ar[i] = creal(s2a);
      s2ai[i] = cimag(s2a);
      s2p[i] = testfn_I_p(ff[i], c2, c3, x0);
    }
    double start=((double)clock())/CLOCKS_PER_SEC;
    double complex rs0=wip_phase(ff,N,ff,N,s1ar,s1ai,s1p,s2ar,s2ai,s2p,Swhitenoise,1.0,-1.0,-1.0);  
    double end=((double)clock())/CLOCKS_PER_SEC;
    double complex err=rs0-rex;
    printf( "ran %i steps.\n",wip_count);
    printf( "result= %g + %gi  err= %g + %gi\n",creal(rs0),cimag(rs0),creal(err),cimag(err));
    printf( "total time= %g   time per eval= %g\n",end-start,(end-start)/wip_count);
    printf( "----------------------------------------------------\n");
  }

  if(1){
    N=N0*256;
    printf("\nInterpQuad(vectored) IP (N->%i):\n",N);
    int i;
    double ff[N],s1ar[N],s2ar[N],s1ai[N],s2ai[N],s1p[N],s2p[N];
    for(i=0;i<N;i++){
      ff[i] = i/(N-1.0)*(fmax-fmin)+fmin;
      s1ar[i] = 1.0;
      s1ai[i] = 0.0;
      s1p[i] = 0.0;
      double complex s2a=a0*testfn_dI_A(ff[i], c2, c3, sig, x0, xc);
      s2ar[i] = creal(s2a);
      s2ai[i] = cimag(s2a);
      s2p[i] = testfn_I_p(ff[i], c2, c3, x0);
    }
    double start=((double)clock())/CLOCKS_PER_SEC;
    double complex rs0=wip_phase(ff,N,ff,N,s1ar,s1ai,s1p,s2ar,s2ai,s2p,Swhitenoise,1.0,-1.0,-1.0);  
    double end=((double)clock())/CLOCKS_PER_SEC;
    double complex err=rs0-rex;
    printf( "ran %i steps.\n",wip_count);
    printf( "result= %g + %gi  err= %g + %gi\n",creal(rs0),cimag(rs0),creal(err),cimag(err));
    printf( "total time= %g   time per eval= %g\n",end-start,(end-start)/wip_count);
    printf( "----------------------------------------------------\n");
  }

}
