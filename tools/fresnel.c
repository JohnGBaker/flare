#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_complex.h>

#include "constants.h"
#include "struct.h"
#include "fresnel.h"


#define acc 1e-8

/* Used for previous version of branching between cases 1a, 1b and 2 in ComputeInt */
/* static const double betathreshold = 1./10; */
/* static const double invbetathreshold = 10.; */
/* static const double lambdathreshold = 0.01; */

/* Used for branching between cases 1a, 1b and 2 in ComputeInt */
static const double factorbranchingcondition = 0.2154434690031884;
static const double case1athreshold = 0.1;

static const double coeffdeltaC[13] = {-0.1,-0.0462962962963,-0.0230769230769,-0.0136554621849,-0.00899470899471,-0.00636363636364,-0.00473664266768,-0.00366161616162,-0.00291467938527,-0.00237483953787,-0.0019721019721,-0.00166370896185,-0.00142235123367};
static const double coeffdeltaS[13] = {-0.0714285714286,-0.0318181818182,-0.0174603174603,-0.0109649122807,-0.00750988142292,-0.00546058879392,-0.00414746543779,-0.00325630252101,-0.00262408157145,-0.00215946843854,-0.00180809015222,-0.00153594771242,-0.0013209013209};

static const double fcoeff[12] = {0.3989422275318915, 1.839999368141345e-7, -0.2992429124371606, 0.0029447692036381, 2.48395059327971, 3.897790465813753, -140.2709114924737, 952.371347498219, -3660.703989590287, 8649.34061360465, -11740.93693284928, 7032.8915074629};
static const double gcoeff[12] = {0, 0.1994718004197629, -0.000125952146250442, -0.7386823670326993, -0.353224701618227, 19.45013452321292, -97.8846433217499, 221.1470704754588, -32.29876222016862, -1102.56036216157, 2554.184335427029, -1962.422566478277};

static double invn[12] = {1., 0.5, 0.3333333333333333, 0.25, 0.2, 0.1666666666666667, 0.1428571428571428, 0.125, 0.1111111111111111, 0.1, 0.0909090909090909, 0.0833333333333333};

/* Function computing cexp(I*x)-1 without loss of accuracy*/
static double complex cexpm1i(double x){
  double rr;//rr=cos(x)-1=2*sin(x/2)**2
  double ri=sin(x);
  double sinhalf=sin(0.5*x);
  rr=-2.0*sinhalf*sinhalf;
  return rr+I*ri;
};

/* Function computing the Fresnel integral for any real (positive) number, with 1e-8 accuracy */
/* Reference used: Computation of Fresnel Integrals. II by K. Mielenz, [J. Res. Natl. Inst. Stand. Technol. 105, 589 (2000)] */
static double complex FresnelE(const double x) {
  double complex res;

  /* Check argument */
  if(x<0) return -FresnelE(-x);

  /* If x<2, Taylor expansion */
  else if(x<2) {
    double C = 0;
    double S = 0;
    double x3 = x*x*x;
    double x4 = x3*x;
    int n = 0;
    double deltaC = x;
    double deltaS = x3/3.;
    while(true) {
      C += deltaC;
      S += deltaS;
      //deltaC = -(4*n+1.)/(4*n+5.)/(2*n+1.)/(2*n+2.) * x4 * deltaC;
      //deltaS = -(4*n+3.)/(4*n+7.)/(2*n+2.)/(2*n+3.) * x4 * deltaS;
      deltaC = coeffdeltaC[n] * x4 * deltaC;
      deltaS = coeffdeltaS[n] * x4 * deltaS;
      n++;
      if( (fabs(deltaC)<acc && fabs(deltaS)<acc) || n>12) break;
    }
    res = C + I*S;
  }

  /* If x>2, Mielenz-Boersma formula using tabulated coeffs */
  else {
    double xinv = 1./x;
    double xinv2 = xinv*xinv;
    double y[12] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    y[0] = xinv;
    for(int i=1; i<12; i++) y[i] = y[i-1]*xinv2;
    double f = 0;
    double g = 0;
    for(int i=0; i<12; i++) {
      f += fcoeff[i] * y[i];
      g += gcoeff[i] * y[i];
    }
    res = sqrt(PI/2) * ( (1+I)/2 - (g+I*f) * cexp(I*x*x));
  }

  return res;
}

static double complex ComputeIntCase1a(
  gsl_vector* coeffsAreal,         /* */
  gsl_vector* coeffsAimag,         /* */
  double p0,                       /* */
  double p1,                       /* */
  double p2,                       /* */  
  double eps)                      /* */
{
  /* Setting up the coefficients and the required accuracy */
  double complex coeffs[12] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  for(int i=0; i<=3; i++) {
    coeffs[i] = gsl_vector_get(coeffsAreal, i+1) + I*gsl_vector_get(coeffsAimag, i+1);
  }
  int maxj0 = 3;
  double eps2 = eps*eps;
  double maxA = fabs(creal(coeffs[0])) + fabs(cimag(coeffs[0]));
  double acctol = 1.7e-5 / (maxA * eps);

  /* Polynomial given by the expansion of cexp(I*p1*x) */
  double complex poly1[8] = {0.,0.,0.,0.,0.,0.,0.,0.};
  double p1abseps = fabs(p1)*eps;
  double complex term1 = 1.; double term1absmax = 1.;
  int maxj1 = 0;
  for(int i=0; i<8; i++) {
    poly1[i] = term1;
    term1absmax = term1absmax * p1abseps * invn[i]; /* Note: invn[i] = 1/(i+1) */ 
    if(term1absmax<acctol) break;
    term1 = term1 * I*p1 * invn[i];
    maxj1++;
  }

  /* Polynomial given by the expansion of cexp(I*p2*x2) */
  double complex poly2[11] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  double p2abseps2 = fabs(p2)*eps*eps;
  double complex term2 = 1.; double term2absmax = 1.;
  int maxj2 = 0;
  for(int i=0; i<4; i++) {
    poly2[2*i] = term2;
    term2absmax = term2absmax * p2abseps2 * invn[i]; /* Note: invn[i] = 1/(i+1) */ 
    if(term2absmax<acctol) break;
    term2 = term2 * I*p2 * invn[i];
    maxj2 += 2;
  }

  /* Multiplying in place the polynomials - we cut according to the fact that we allow for only 12 coeffs up to x^11 (assuming the rest would be negligible anyway), and we know the constant in poly1 and poly2 is always 1. */
  for(int i=min(11, maxj0+maxj1); i>=0; i--) {
    for(int j=max(1, i-maxj0); j<=min(i, maxj1); j++) {
      coeffs[i] += coeffs[i-j] * poly1[j]; 
    }
  }
  maxj0 = min(11, maxj0 + maxj1);
  for(int i=min(11, maxj0+maxj2); i>=0; i--) {
    for(int j=max(1, i-maxj0); j<=min(i, maxj2); j++) {
      coeffs[i] += coeffs[i-j] * poly2[j];
    }
  }
  maxj0 = min(11, maxj0 + maxj2);
  
  /* Computing the integral itself */
  double epstab[maxj0+2];
  epstab[0] = 1.;
  for(int i=1; i<=maxj0+1; i++) {epstab[i] = eps * epstab[i-1];}
  double complex res = 0.;
  for(int i=0; i<=maxj0; i++) {
    res += coeffs[i] * epstab[i+1] * invn[i]; /* Note: invn[i] = 1/(i+1) */
  }
  res = cexp(I*p0) * res;
  return res;
}

static double complex ComputeIntCase1b(
  gsl_vector* coeffsAreal,         /* */
  gsl_vector* coeffsAimag,         /* */
  double p0,                       /* */
  double p1,                       /* */
  double p2,                       /* */  
  double eps)                      /* */
{
  /* Setting up the coefficients and the required accuracy */
  double complex coeffs[12] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  for(int i=0; i<=3; i++) {
    coeffs[i] = gsl_vector_get(coeffsAreal, i+1) + I*gsl_vector_get(coeffsAimag, i+1);
  }
  int maxj0 = 3;
  double eps2 = eps*eps;
  double maxA = fabs(creal(coeffs[0])) + fabs(cimag(coeffs[0]));
  double acctol = 1.7e-5 / maxA;

  /* Polynomial given by the expansion of cexp(I*p2*x2) */
  double complex poly2[11] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  double p2abseps2 = fabs(p2)*eps*eps;
  double complex term2 = 1.; double term2absmax = 1.;
  int maxj2 = 0;
  for(int i=0; i<4; i++) {
    poly2[2*i] = term2;
    term2absmax = term2absmax * p2abseps2 * invn[i]; /* Note: invn[i] = 1/(i+1) */ 
    if(term2absmax<acctol) break;
    term2 = term2 * I*p2 * invn[i];
    maxj2 += 2;
  }

  /* Multiplying in place the polynomials - we cut according to the fact that we allow for only 12 coeffs up to x^11 (assuming the rest would be negligible anyway), and we know the constant in poly2 is always 1. */
  for(int i=min(11, maxj0+maxj2); i>=0; i--) {
    for(int j=max(1, i-maxj0); j<=min(i, maxj2); j++) {
      coeffs[i] += coeffs[i-j] * poly2[j];
    }
  }
  maxj0 = min(11, maxj0 + maxj2);  

  /* Go down the recursion relation */
  double epstab[maxj0];
  epstab[0] = 1.;
  for(int i=1; i<maxj0; i++) {epstab[i] = eps * epstab[i-1];}
  double invp1 = 1./p1;
  double complex coeffE1 = 0.; double complex coeffE1minus1 = 0;
  for(int i=maxj0; i>=1; i--) {
    coeffE1 += -I*invp1 * epstab[i] * coeffs[i];
    coeffs[i-1] += I*invp1 * i * coeffs[i];
  }
  coeffE1minus1 = -I*invp1 * coeffs[0];

  /* Result */
  double complex cexpm1ip1eps = cexpm1i(p1*eps);
  double complex res = cexp(I*p0) * (coeffE1*(cexpm1ip1eps+1.) + coeffE1minus1*cexpm1ip1eps);
  return res;
}

static double complex ComputeIntCase2(
  gsl_vector* coeffsAreal,         /* */
  gsl_vector* coeffsAimag,         /* */
  double p0,                       /* */
  double p1,                       /* */
  double p2,                       /* */  
  double eps,                      /* */
  double sqrtp2)                   /* */
{
  /* Amplitude coeffs and powers of eps */
  double complex coeffsA[4];
  for(int i=0; i<=3; i++) coeffsA[i] = gsl_vector_get(coeffsAreal, i+1) + I*gsl_vector_get(coeffsAimag, i+1);
  double poweps[3];
  poweps[0] = 1; poweps[1] = eps; poweps[2] = eps*eps;

  double complex coeffE2 = 0; double complex coefffresnel = 0; double complex constant = 0;
  double inv2p2 = 1. / (2.*p2);
  /* n=3 and n=2 */
  for(int i=3; i>=2; i--) {
    coeffE2 += -I*inv2p2 * poweps[i-1] * coeffsA[i];
    coeffsA[i-2] += I*inv2p2 * (i-1) * coeffsA[i];
    coeffsA[i-1] += -p1 * inv2p2 * coeffsA[i];
  }
  /* n=1 */
  coeffE2 += -I*inv2p2 * coeffsA[1];
  constant = I*inv2p2 * coeffsA[1];
  coefffresnel += -p1 * inv2p2 * coeffsA[1];
  /* n=0 */
  coefffresnel += coeffsA[0];

  /* Fresnel integral and complex exponential */
  double sqrtp2eps = sqrtp2 * eps;
  double complex E2 = cexp(I * (p1*eps + p2*eps*eps));
  int conjug = 0;
  if(p2<0) {conjug = 1; p1 = -p1; p2 = -p2;}
  double invsqrtp2 = 1./sqrtp2;
  double halfratio = 0.5 * p1 * invsqrtp2;
  double complex fresnel = cexp(-I*halfratio*halfratio) * invsqrtp2 * (FresnelE(halfratio + sqrtp2eps) - FresnelE(halfratio));
  if(conjug) fresnel = conj(fresnel);

  /* Result */
  double complex res = cexp(I*p0) * (constant + coeffE2*E2 + coefffresnel*fresnel);
  return res;
}


double complex ComputeInt(
  gsl_matrix* splinecoeffsAreal,         /*  */
  gsl_matrix* splinecoeffsAimag,         /*  */
  gsl_matrix* quadsplinecoeffsphase)     /*  */
{
  double complex res = 0.;
  /* Number of points - i.e. nb of intervals + 1 */
  /* Assumes that the dimensions match and that the frequency vectors are the same between the different splinecoeffs */
  int nbpts = (int) quadsplinecoeffsphase->size1;

clock_t tbegtotal = clock();
  for(int j=0; j<nbpts-1; j++) {
    double eps = gsl_matrix_get(quadsplinecoeffsphase, j+1, 0) - gsl_matrix_get(quadsplinecoeffsphase, j, 0);
    double p0 = gsl_matrix_get(quadsplinecoeffsphase, j, 1);
    double p1 = gsl_matrix_get(quadsplinecoeffsphase, j, 2);
    double p2 = gsl_matrix_get(quadsplinecoeffsphase, j, 3);

    double sqrtp2 = sqrt(fabs(p2)); double absp1 = fabs(p1);

    gsl_vector_view viewAreal = gsl_matrix_row(splinecoeffsAreal, j);
    gsl_vector_view viewAimag = gsl_matrix_row(splinecoeffsAimag, j);

    if(absp1*eps<case1athreshold) {
      res += ComputeIntCase1a(&viewAreal.vector, &viewAimag.vector, p0, p1, p2, eps);
    }
    else {
      if(fabs(p2)*eps*eps < factorbranchingcondition*pow(absp1*eps, 4./3)) {
	res += ComputeIntCase1b(&viewAreal.vector, &viewAimag.vector, p0, p1, p2, eps);
      }
      else {
	res += ComputeIntCase2(&viewAreal.vector, &viewAimag.vector, p0, p1, p2, eps, sqrtp2);
      }
    }
  }

  return res;
}
