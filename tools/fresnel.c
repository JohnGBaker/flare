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

/* Used for previous version of branching between cases 1a, 1b and 2 in ComputeInt */
/* static const double factorbranchingcondition = 0.2154434690031884; */
/* static const double case1athreshold = 0.1; */
/* Used for case 0 - hoping to cover previous instability */
/* static const double case0threshold = 0.05; */

/* Used for branching between cases 1a, 1b, 2, 3 and 4 in ComputeInt */
/* FIXME: another round of refinement would be welcome, case 3 has some cases with 1e-4 error */
static const double p2threshold = 0.02;
static const double p1threshold1 = 0.1;
static const double p1threshold2 = 100.;
static const double p2p1slopethreshold = 0.2;

static const double coeffdeltaC[13] = {-0.1,-0.0462962962963,-0.0230769230769,-0.0136554621849,-0.00899470899471,-0.00636363636364,-0.00473664266768,-0.00366161616162,-0.00291467938527,-0.00237483953787,-0.0019721019721,-0.00166370896185,-0.00142235123367};
static const double coeffdeltaS[13] = {-0.0714285714286,-0.0318181818182,-0.0174603174603,-0.0109649122807,-0.00750988142292,-0.00546058879392,-0.00414746543779,-0.00325630252101,-0.00262408157145,-0.00215946843854,-0.00180809015222,-0.00153594771242,-0.0013209013209};

static const double fcoeff[12] = {0.3989422275318915, 1.839999368141345e-7, -0.2992429124371606, 0.0029447692036381, 2.48395059327971, 3.897790465813753, -140.2709114924737, 952.371347498219, -3660.703989590287, 8649.34061360465, -11740.93693284928, 7032.8915074629};
static const double gcoeff[12] = {0, 0.1994718004197629, -0.000125952146250442, -0.7386823670326993, -0.353224701618227, 19.45013452321292, -97.8846433217499, 221.1470704754588, -32.29876222016862, -1102.56036216157, 2554.184335427029, -1962.422566478277};

static double invn[12] = {1., 0.5, 0.3333333333333333, 0.25, 0.2, 0.1666666666666667, 0.1428571428571428, 0.125, 0.1111111111111111, 0.1, 0.0909090909090909, 0.0833333333333333};

/* Use for Chebyshev, case 3 */
static double ChebyshevNodes[6] = {0.965925826289068,0.7071067811865475,0.2588190451025207,-0.2588190451025207,-0.7071067811865475,-0.965925826289068};
static double ChebyshevWeights[6][6] = {{0.3333333333333333,0.3333333333333333,0.3333333333333333,0.3333333333333333,0.3333333333333333,0.3333333333333333},{0.3219752754296894,0.2357022603955158,0.0862730150341736,-0.0862730150341736,-0.2357022603955158,-0.3219752754296894},{0.2886751345948129,0.,-0.2886751345948129,-0.2886751345948129,0.,0.2886751345948129},{0.2357022603955158,-0.2357022603955158,-0.2357022603955158,0.2357022603955158,0.2357022603955158,-0.2357022603955158},{0.1666666666666667,-0.3333333333333333,0.1666666666666667,0.1666666666666667,-0.3333333333333333,0.1666666666666667},{0.0862730150341736,-0.2357022603955158,0.3219752754296894,-0.3219752754296894,0.2357022603955158,-0.0862730150341736}};

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

/* In fact equivalent to Case 2 */
//static double complex ComputeIntCase0(
//   gsl_vector* coeffsAreal,         /* */
//   gsl_vector* coeffsAimag,         /* */
//   double p0,                       /* */
//   double p1,                       /* */
//   double p2,                       /* */
//   double eps)                      /* */
// {
//   /* Get poly coeffs and rescale by epsilon */
//   double complex coeffs[4] = {0.,0.,0.,0.};
//   double epspow[4] = {0.,0.,0.,0.};
//   epspow[0] = 1.; epspow[1] = eps; epspow[2] = eps*eps; epspow[3] = eps * epspow[2];
//   for(int i=0; i<=3; i++) {
//     coeffs[i] = epspow[i] * (gsl_vector_get(coeffsAreal, i+1) + I*gsl_vector_get(coeffsAimag, i+1));
//   }
//   double p0r = p0;
//   double p1r = p1*epspow[1];
//   double p2r = p2*epspow[2];
//
//   /* In case p2<0, conjugation */
//   int conjug = 0;
//   if(p2r<0) {
//     conjug=1;
//     p2r = -p2r;
//     p1r = -p1r;
//     p0r = -p0r;
//     for(int i=0; i<=3; i++) coeffs[i] = conj(coeffs[i]) ;
//   }
//
//   //
//   //printf("case0\n");
//
//   /* Prefactor and bounds after change of variable */
//   double sqrtp2r = sqrt(p2r);
//   double complex factor = 1./sqrtp2r*cexp(I*p0r);
//   double a = p1r/(2*sqrtp2r);
//   double b = p1r/(2*sqrtp2r) + sqrtp2r;
//   /* Coeffs in polynomial changed by shift */
//   double alpha = 1./sqrtp2r;
//   double beta = -p1r/(2*p2r);
//   double complex a0 = coeffs[0]; double complex a1 = coeffs[1]; double complex a2 = coeffs[2]; double complex a3 = coeffs[3];
//   double alpha2 = alpha*alpha; double beta2 = beta*beta;
//   double complex a0p = a0 + a1*beta + a2*beta2 + a3*beta2*beta;
//   double complex a1p = a1*alpha + 2*a2*alpha*beta + 3*a3*alpha*beta2;
//   double complex a2p = a2*alpha2 + 3*a3*alpha2*beta;
//   double complex a3p = a3*alpha2*alpha;
//
//   /* Result */
//   double complex res = 0.;
//   res = factor*((a0p + I/2*a2p) * (FresnelE(b) - FresnelE(a)) + (1./2*(a3p - I*a1p) - I/2*a2p*b - I/2*a3p*b*b) * cexp(I*b*b) - (1./2*(a3p - I*a1p) - I/2*a2p*a - I/2*a3p*a*a) * cexp(I*a*a));
//   //
//   //printf("%g, %g, %g \n", creal(factor), creal((a0p + I/2*a2p) * (FresnelE(b) - FresnelE(a))), creal(+ (1./2*(a3p - I*a1p) - I/2*a2p*b - I/2*a3p*b*b) * cexp(I*b*b) - (1./2*(a3p - I*a1p) - I/2*a2p*a - I/2*a3p*a*a) * cexp(I*a*a)));
//   //exit(0);
//
//   /* Remembering conjugation */
//   if(conjug) res = conj(res);
//   return res;
// }

double complex ComputeIntCase1a(
  const double complex* coeffsA,         /* */
  const double p1,                       /* */
  const double p2,                       /* */
  const double scale)                    /* */
{
  /* Setting up the coefficients and the required accuracy */
  double complex coeffs[12] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  for(int i=0; i<=3; i++) coeffs[i] = coeffsA[i];
  int maxj0 = 3;
  double acctol = 1.e-5/scale;

  /* Polynomial given by the expansion of cexp(I*p1*x) */
  double complex poly1[8] = {0.,0.,0.,0.,0.,0.,0.,0.};
  double p1abs = fabs(p1);
  double complex term1 = 1.; double term1absmax = 1.;
  int maxj1 = 0;
  for(int i=0; i<8; i++) {
    poly1[i] = term1;
    term1absmax = term1absmax * p1abs * invn[i]; /* Note: invn[i] = 1/(i+1) */
    if(term1absmax<acctol) break;
    term1 = term1 * I*p1 * invn[i];
    maxj1++;
  }

  /* Polynomial given by the expansion of cexp(I*p2*x2) */
  double complex poly2[11] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  double p2abs = fabs(p2);
  double complex term2 = 1.; double term2absmax = 1.;
  int maxj2 = 0;
  for(int i=0; i<4; i++) {
    poly2[2*i] = term2;
    term2absmax = term2absmax * p2abs * invn[i]; /* Note: invn[i] = 1/(i+1) */
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
  double complex res = 0.;
  for(int i=0; i<=maxj0; i++) {
    res += coeffs[i] * invn[i]; /* Note: invn[i] = 1/(i+1) */
  }
  return res;
}

double complex ComputeIntCase1b(
  const double complex* coeffsA,         /* */
  const double p1,                       /* */
  const double p2,                       /* */
  const double scale)                    /* */
{
  /* Setting up the coefficients and the required accuracy */
  double complex coeffs[12] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  for(int i=0; i<=3; i++) coeffs[i] = coeffsA[i];
  int maxj0 = 3;
  double acctol = 1.e-5/scale;

  /* Polynomial given by the expansion of cexp(I*p2*x2) */
  double complex poly2[11] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  double p2abs = fabs(p2);
  double complex term2 = 1.; double term2absmax = 1.;
  int maxj2 = 0;
  for(int i=0; i<4; i++) {
    poly2[2*i] = term2;
    term2absmax = term2absmax * p2abs * invn[i]; /* Note: invn[i] = 1/(i+1) */
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
  double invp1 = 1./p1;
  double complex coeffE1 = 0.; double complex coeffE1minus1 = 0;
  for(int i=maxj0; i>=1; i--) {
    coeffE1 += -I*invp1 * coeffs[i];
    coeffs[i-1] += I*invp1 * i * coeffs[i];
  }
  coeffE1minus1 = -I*invp1 * coeffs[0];
  /* Result */
  double complex cexpm1ip1 = cexpm1i(p1);
  double complex res = (coeffE1*(cexpm1ip1+1.) + coeffE1minus1*cexpm1ip1);
  return res;
}

double complex ComputeIntCase2(
  const double complex* coeffsA,         /* */
  const double p1,                       /* */
  const double p2)                       /* */
{
  /* Setting up the coefficients */
  double complex coeffs[4] = {0.,0.,0.,0.};
  for(int i=0; i<=3; i++) coeffs[i] = coeffsA[i];
  double complex coeffE2 = 0; double complex coefffresnel = 0; double complex constant = 0;
  /* Conjugate if required */
  int conjug = 0; double p1b = p1; double p2b = p2;
  if(p2<0) { conjug = 1; p1b = -p1; p2b = -p2; for(int i=0; i<=3; i++) coeffs[i] = conj(coeffs[i]); }

  /* n=3 and n=2 */
  double inv2p2 = 1. / (2.*p2b);
  for(int i=3; i>=2; i--) {
    coeffE2 += -I*inv2p2 * coeffs[i];
    coeffs[i-2] += I*inv2p2 * (i-1) * coeffs[i];
    coeffs[i-1] += -p1b * inv2p2 * coeffs[i];
  }
  /* n=1 */
  coeffE2 += -I*inv2p2 * coeffs[1];
  constant = I*inv2p2 * coeffs[1];
  coefffresnel += -p1b * inv2p2 * coeffs[1];
  /* n=0 */
  coefffresnel += coeffs[0];

  /* Fresnel integral and complex exponential */
  double complex E2 = cexp(I * (p1b + p2b));
  double sqrtp2 = sqrt(p2b);
  double invsqrtp2 = 1./sqrtp2;
  double halfratio = 0.5 * p1b * invsqrtp2;
  double complex fresnel = cexp(-I*halfratio*halfratio) * invsqrtp2 * (FresnelE(halfratio + sqrtp2) - FresnelE(halfratio));

  /* Result */
  double complex res = (constant + coeffE2*E2 + coefffresnel*fresnel);
  if(conjug) res = conj(res); /* Remember possible conjugation */
  return res;
}

double complex ComputeIntCase3(
  const double complex* coeffsA,         /* */
  const double p1,                       /* */
  const double p2)                       /* */
{
  /* Prepare change of variables */
  double r = p2/p1;

  /* Compute integrand F on nodes */
  double Fvalues[6] = {0.,0.,0.,0.,0.,0.};
  double xk; double y; double sqrtterm; double factor; double x;
  for(int i=0; i<6; i++) {
    xk = ChebyshevNodes[i];
    y = (1.+xk)/2.;
    sqrtterm = sqrt(1. + 4.*r*(1.+r)*y);
    factor = (1.+r)/sqrtterm;
    x = (2.*(1.+r)*y)/(1. + sqrtterm);
    Fvalues[i] = coeffsA[3];
    for(int j=2; j>=0; j--) Fvalues[i] = coeffsA[j] + x*Fvalues[i];
    Fvalues[i] = factor * Fvalues[i];
  }

  /* Compute Chebyshev coefficients ci */
  double complex c[6] = {0.,0.,0.,0.,0.,0.};
  for(int i=0; i<6; i++) {
    for(int j=0; j<6; j++) c[i] += ChebyshevWeights[i][j] * Fvalues[j];
  }

  /* Coefficients of polynomial in 1/(I*p), p=p1*(1+r) */
  double p = p1*(1.+r);
  double complex coeffszcos[6] = {0.,0.,0.,0.,0.,0.};
  double complex coeffszsin[6] = {0.,0.,0.,0.,0.,0.};
  coeffszcos[0] = c[1] + c[3] + c[5];
  coeffszcos[1] = -8.*c[2] - 32.*c[4];
  coeffszcos[2] = 96.*c[3] + 800.*c[5];
  coeffszcos[3] = -1536.*c[4];
  coeffszcos[4] = 30720.*c[5];
  coeffszcos[5] = 0.;
  coeffszsin[0] = 0.5*c[0] + c[2] + c[4];
  coeffszsin[1] = -2.*c[1] - 18.*c[3] - 50.*c[5];
  coeffszsin[2] = 16.*c[2] + 320.*c[4];
  coeffszsin[3] = -192.*c[3] -6720.*c[5];
  coeffszsin[4] = 3072.*c[4];
  coeffszsin[5] = -61440.*c[5];

  /* Coefficients of cos(p/2) and sin(p/2) */
  double complex z = 1./(I*p);
  double complex coeffcos = coeffszcos[5];
  for(int i=4; i>=0; i--) coeffcos = coeffszcos[i] + z*coeffcos;
  double complex coeffsin = coeffszsin[5];
  for(int i=4; i>=0; i--) coeffsin = coeffszsin[i] + z*coeffsin;
  coeffcos = z * coeffcos;
  coeffsin = z * coeffsin;

  /* Result */
  double pov2 = p/2.;
  return cexp(I*pov2) * (2*cos(pov2)*coeffcos + 2*I*sin(pov2)*coeffsin);
}

double complex ComputeIntCase4(
  const double complex* coeffsA,         /* */
  const double p1,                       /* */
  const double p2)                       /* */
{
  double complex a0 = coeffsA[0];
  double complex a1 = coeffsA[1];
  double complex a2 = coeffsA[2];
  double complex a3 = coeffsA[3];
  double r = 2.*p2/p1;
  double inv1plusr = 1./(1.+r); double inv1plusr2 = inv1plusr*inv1plusr;
  double inv1plusr4 = inv1plusr2*inv1plusr2;
  double invp1 = 1./p1; double invp12 = invp1*invp1; double invp13 = invp12*invp1;

  double complex term1 = cexp(I*(p1 + p2))*(-I*invp1*(a0 + a1 + a2 + a3)*inv1plusr + invp12*(a1 + 2.*a2 + 3.*a3)*inv1plusr2 + I*invp13*(2.*a2 + 6.*a3 + r*(3.*a3 - a1))*inv1plusr4);
  double complex term2 = - (-I*invp1*a0 + invp12*a1 + I*invp13*(2.*a2 - a1*r));
  return term1 + term2;
}

double complex ComputeInt(
  gsl_matrix* splinecoeffsAreal,         /*  */
  gsl_matrix* splinecoeffsAimag,         /*  */
  gsl_matrix* quadsplinecoeffsphase)     /*  */
{
  double complex res = 0.;
  double complex resint = 0.;
  /* Number of points - i.e. nb of intervals + 1 */
  /* Assumes that the dimensions match and that the frequency vectors are the same between the different splinecoeffs */
  int nbpts = (int) quadsplinecoeffsphase->size1;

  //clock_t tbegtotal = clock();
  for(int j=0; j<nbpts-1; j++) {
    double eps = gsl_matrix_get(quadsplinecoeffsphase, j+1, 0) - gsl_matrix_get(quadsplinecoeffsphase, j, 0);
    double epspow[4];
    epspow[0] = 1.;
    epspow[1] = eps;
    epspow[2] = eps*eps;
    epspow[3] = eps*epspow[2];
    double p0 = gsl_matrix_get(quadsplinecoeffsphase, j, 1);

    /* Rescale the interval to [0,1] and scale out the constant amplitude term */
    double p1 = gsl_matrix_get(quadsplinecoeffsphase, j, 2) * eps;
    double p2 = gsl_matrix_get(quadsplinecoeffsphase, j, 3) * epspow[2];
    double complex A0 = gsl_matrix_get(splinecoeffsAreal, j, 1) + I*gsl_matrix_get(splinecoeffsAimag, j, 1);
    double complex A0inv = 1./A0;
    double A0abs = cabs(A0); /* Used for accuracy requirement of adaptive Taylor expansions Ia and Ib */
    double complex coeffsA[4] = {0.,0.,0.,0.};
    coeffsA[0] = 1.;
    for(int i=1; i<=3; i++) coeffsA[i] = epspow[i] * A0inv * (gsl_matrix_get(splinecoeffsAreal, j, i+1) + I*gsl_matrix_get(splinecoeffsAimag, j, i+1));

    /* Factor scaled out */
    double complex factor = eps * A0 * cexp(I*p0);

    double absp1 = fabs(p1); double absp2 = fabs(p2);
    if(absp1<p1threshold1 && absp2<p2threshold) {
      resint = factor * ComputeIntCase1a(coeffsA, p1, p2, A0abs);
    }
    else if(p1threshold1<=absp1 && absp1<p1threshold2 && absp2<p2threshold) {
      resint = factor * ComputeIntCase1b(coeffsA, p1, p2, A0abs);
    }
    else if(absp2>=p2p1slopethreshold*absp1 && absp2>=p2threshold) {
      resint = factor * ComputeIntCase2(coeffsA, p1, p2);
    }
    else if(absp2<p2p1slopethreshold*absp1 && absp2>=p2threshold && absp1<p1threshold2) {
      resint = factor * ComputeIntCase3(coeffsA, p1, p2);
    }
    else {
      resint = factor * ComputeIntCase4(coeffsA, p1, p2);
    }
    res += resint;
  }

  return res;
}
