/**
 * \author Sylvain Marsat, University of Maryland - NASA GSFC
 *
 * \brief C code for the geometric coefficients entering the response for LISA-like detectors.
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
#include "waveform.h"
#include "LISAgeometry.h"

#include <time.h> /* for testing */

/****************************************************************/
/********* Coefficients for the geometric response **************/

/* External storage for cos, sin and coefficients */
static double coeffn1Hn1crossconst, coeffn1Hn1plusconst, coeffn2Hn2crossconst, coeffn2Hn2plusconst, coeffn3Hn3crossconst, coeffn3Hn3plusconst;
static double coeffn1Hn1pluscos[4];
static double coeffn1Hn1plussin[4];
static double coeffn2Hn2pluscos[4];
static double coeffn2Hn2plussin[4];
static double coeffn3Hn3pluscos[4];
static double coeffn3Hn3plussin[4];
static double coeffn1Hn1crosscos[4];
static double coeffn1Hn1crosssin[4];
static double coeffn2Hn2crosscos[4];
static double coeffn2Hn2crosssin[4];
static double coeffn3Hn3crosscos[4];
static double coeffn3Hn3crosssin[4];
static double coeffkn1const, coeffkn2const, coeffkn3const, coeffkp1plusp2const, coeffkp2plusp3const, coeffkp3plusp1const, coeffkp1const, coeffkp2const, coeffkp3const, coeffkRconst;
static double coeffkn1cos[2];
static double coeffkn1sin[2];
static double coeffkn2cos[2];
static double coeffkn2sin[2];
static double coeffkn3cos[2];
static double coeffkn3sin[2];
static double coeffkp1plusp2cos[2];
static double coeffkp1plusp2sin[2];
static double coeffkp2plusp3cos[2];
static double coeffkp2plusp3sin[2];
static double coeffkp3plusp1cos[2];
static double coeffkp3plusp1sin[2];
static double coeffkp1cos[2];
static double coeffkp1sin[2];
static double coeffkp2cos[2];
static double coeffkp2sin[2];
static double coeffkp3cos[2];
static double coeffkp3sin[2];
static double coeffkRcos[2];
static double coeffkRsin[2];

static double cosarray[4];
static double sinarray[4];

#pragma omp threadprivate(coeffn1Hn1crossconst, coeffn1Hn1plusconst, coeffn2Hn2crossconst, coeffn2Hn2plusconst, coeffn3Hn3crossconst, coeffn3Hn3plusconst)
#pragma omp threadprivate(coeffn1Hn1pluscos,coeffn1Hn1plussin,coeffn2Hn2pluscos,coeffn2Hn2plussin,coeffn3Hn3pluscos,coeffn3Hn3plussin)
#pragma omp threadprivate(coeffn1Hn1crosscos,coeffn1Hn1crosssin,coeffn2Hn2crosscos,coeffn2Hn2crosssin,coeffn3Hn3crosscos,coeffn3Hn3crosssin)
#pragma omp threadprivate(coeffkn1const, coeffkn2const, coeffkn3const, coeffkp1plusp2const, coeffkp2plusp3const, coeffkp3plusp1const, coeffkp1const, coeffkp2const, coeffkp3const, coeffkRconst)
#pragma omp threadprivate(coeffkn1cos,coeffkn1sin,coeffkn2cos,coeffkn2sin,coeffkn3cos,coeffkn3sin)
#pragma omp threadprivate(coeffkp1cos,coeffkp1sin,coeffkp2cos,coeffkp2sin,coeffkp3cos,coeffkp3sin)
#pragma omp threadprivate(coeffkp1plusp2cos,coeffkp1plusp2sin,coeffkp2plusp3cos,coeffkp2plusp3sin,coeffkp3plusp1cos,coeffkp3plusp1sin)
#pragma omp threadprivate(coeffkRcos,coeffkRsin,cosarray,sinarray)

/*************************************************************/
/********* Functions for the geometric response **************/

/* Function to convert string input TDI string to TDItag */
TDItag ParseTDItag(char* string) {
  TDItag tag;
  if(strcmp(string, "delayO")==0) tag = delayO;
  else if(strcmp(string, "y12L")==0) tag = y12L;
  else if(strcmp(string, "y12")==0) tag = y12;
  else if(strcmp(string, "TDIXYZ")==0) tag = TDIXYZ;
  else if(strcmp(string, "TDIalphabetagamma")==0) tag = TDIalphabetagamma;
  else if(strcmp(string, "TDIAETXYZ")==0) tag = TDIAETXYZ;
  else if(strcmp(string, "TDIAETalphabetagamma")==0) tag = TDIAETalphabetagamma;
  else if(strcmp(string, "TDIX")==0) tag = TDIX;
  else if(strcmp(string, "TDIalpha")==0) tag = TDIalpha;
  else if(strcmp(string, "TDIAXYZ")==0) tag = TDIAXYZ;
  else if(strcmp(string, "TDIEXYZ")==0) tag = TDIEXYZ;
  else if(strcmp(string, "TDITXYZ")==0) tag = TDITXYZ;
  else if(strcmp(string, "TDIAalphabetagamma")==0) tag = TDIAalphabetagamma;
  else if(strcmp(string, "TDIEalphabetagamma")==0) tag = TDIEalphabetagamma;
  else if(strcmp(string, "TDITalphabetagamma")==0) tag = TDITalphabetagamma;
  else {
    printf("Error in ParseTDItag: string not recognized.\n");
    exit(1);
  }
  return tag;
}

/* Function to convert string input ResponseApprox to tag */
ResponseApproxtag ParseResponseApproxtag(char* string) {
  ResponseApproxtag tag;
  if(strcmp(string, "full")==0) tag = full;
  else if(strcmp(string, "lowfL")==0) tag = lowfL;
  else if(strcmp(string, "lowf")==0) tag = lowf;
  else {
    printf("Error in ParseResponseApproxtag: string not recognized.\n");
    exit(1);
  }
  return tag;
}

/* Compute Solar System Barycenter time tSSB from retarded time at the center of the LISA constellation tL */
double tSSBfromtL(const double tL, const double lambda, const double beta) {
  return tL + R_SI/C_SI*cos(beta)*cos(Omega_SI*tL - lambda) - 1./2*Omega_SI*pow(R_SI/C_SI*cos(beta), 2)*sin(2.*(Omega_SI*tL - lambda));
}
double tLfromttSSB(const double tSSB, const double lambda, const double beta) {
  return tSSB - R_SI/C_SI*cos(beta)*cos(Omega_SI*tSSB - lambda);
}

/* Function cardinal sine */
double sinc(const double x) {
  if (x==0)
    return 1;
  else return sin(x)/x;
}

/* Function to compute, given a value of a sky position and polarization, all the complicated time-independent trigonometric coefficients entering the response */
void SetCoeffsG(const double lambda, const double beta, const double psi) {
  /* Precomputing cosines and sines */
  double coslambda = cos(lambda);
  double sinlambda = sin(lambda);
  double cosbeta = cos(beta);
  double sinbeta = sin(beta);
  double cospsi = cos(psi);
  double sinpsi = sin(psi);

  /* Projection coefficients for hplus in n3.H.n3 */
  /**/
  coeffn3Hn3plusconst = 1./128 * (-4*cospsi*cospsi + 4*sinpsi*sinpsi -27*coslambda*coslambda*cospsi*cospsi -27*sinlambda*sinlambda*sinpsi*sinpsi -4*cosbeta*cosbeta*cospsi*cospsi -4*sinbeta*sinbeta*sinpsi*sinpsi + 4*cosbeta*cosbeta*sinpsi*sinpsi + 4*cospsi*cospsi*sinbeta*sinbeta + 27*coslambda*coslambda*sinpsi*sinpsi + 27*cospsi*cospsi*sinlambda*sinlambda -9*cosbeta*cosbeta*coslambda*coslambda*sinpsi*sinpsi -9*cosbeta*cosbeta*cospsi*cospsi*sinlambda*sinlambda -9*coslambda*coslambda*cospsi*cospsi*sinbeta*sinbeta -9*sinbeta*sinbeta*sinlambda*sinlambda*sinpsi*sinpsi + 9*cosbeta*cosbeta*coslambda*coslambda*cospsi*cospsi + 9*cosbeta*cosbeta*sinlambda*sinlambda*sinpsi*sinpsi + 9*coslambda*coslambda*sinbeta*sinbeta*sinpsi*sinpsi + 9*cospsi*cospsi*sinbeta*sinbeta*sinlambda*sinlambda -54*sqrt3*coslambda*sinlambda*sinpsi*sinpsi + 54*sqrt3*coslambda*cospsi*cospsi*sinlambda -144*coslambda*cospsi*sinbeta*sinlambda*sinpsi -72*sqrt3*coslambda*coslambda*cospsi*sinbeta*sinpsi -18*sqrt3*coslambda*sinbeta*sinbeta*sinlambda*sinpsi*sinpsi -18*sqrt3*cosbeta*cosbeta*coslambda*cospsi*cospsi*sinlambda + 18*sqrt3*coslambda*cospsi*cospsi*sinbeta*sinbeta*sinlambda + 18*sqrt3*cosbeta*cosbeta*coslambda*sinlambda*sinpsi*sinpsi + 72*sqrt3*cospsi*sinbeta*sinlambda*sinlambda*sinpsi);
  /**/
  coeffn3Hn3pluscos[0] = 1./16*cosbeta * (-9*cospsi*cospsi*sinbeta*sinlambda + 9*sinbeta*sinlambda*sinpsi*sinpsi + 18*coslambda*cospsi*sinpsi -7*sqrt3*coslambda*sinbeta*sinpsi*sinpsi + 7*sqrt3*coslambda*cospsi*cospsi*sinbeta + 14*sqrt3*cospsi*sinlambda*sinpsi);
  /**/
  coeffn3Hn3pluscos[1] = -3./64 * (-3*sinpsi*sinpsi + 3*cospsi*cospsi -6*coslambda*coslambda*cospsi*cospsi -6*sinlambda*sinlambda*sinpsi*sinpsi -3*cosbeta*cosbeta*sinpsi*sinpsi -3*cospsi*cospsi*sinbeta*sinbeta + 3*cosbeta*cosbeta*cospsi*cospsi + 3*sinbeta*sinbeta*sinpsi*sinpsi + 6*coslambda*coslambda*sinpsi*sinpsi + 6*cospsi*cospsi*sinlambda*sinlambda -2*cosbeta*cosbeta*coslambda*coslambda*sinpsi*sinpsi -2*cosbeta*cosbeta*cospsi*cospsi*sinlambda*sinlambda -2*coslambda*coslambda*cospsi*cospsi*sinbeta*sinbeta -2*sinbeta*sinbeta*sinlambda*sinlambda*sinpsi*sinpsi + 2*cosbeta*cosbeta*coslambda*coslambda*cospsi*cospsi + 2*cosbeta*cosbeta*sinlambda*sinlambda*sinpsi*sinpsi + 2*coslambda*coslambda*sinbeta*sinbeta*sinpsi*sinpsi + 2*cospsi*cospsi*sinbeta*sinbeta*sinlambda*sinlambda -32*coslambda*cospsi*sinbeta*sinlambda*sinpsi);
  /**/
  coeffn3Hn3pluscos[2] = -1./16*cosbeta * (-6*coslambda*cospsi*sinpsi -3*sinbeta*sinlambda*sinpsi*sinpsi + 3*cospsi*cospsi*sinbeta*sinlambda + sqrt3*coslambda*cospsi*cospsi*sinbeta -sqrt3*coslambda*sinbeta*sinpsi*sinpsi + 2*sqrt3*cospsi*sinlambda*sinpsi);
  /**/
  coeffn3Hn3pluscos[3] = 1./128 * (-3*coslambda*coslambda*cospsi*cospsi -3*sinlambda*sinlambda*sinpsi*sinpsi + 3*coslambda*coslambda*sinpsi*sinpsi + 3*cospsi*cospsi*sinlambda*sinlambda + cosbeta*cosbeta*coslambda*coslambda*cospsi*cospsi + cosbeta*cosbeta*sinlambda*sinlambda*sinpsi*sinpsi + coslambda*coslambda*sinbeta*sinbeta*sinpsi*sinpsi + cospsi*cospsi*sinbeta*sinbeta*sinlambda*sinlambda -cosbeta*cosbeta*coslambda*coslambda*sinpsi*sinpsi -cosbeta*cosbeta*cospsi*cospsi*sinlambda*sinlambda -coslambda*coslambda*cospsi*cospsi*sinbeta*sinbeta -sinbeta*sinbeta*sinlambda*sinlambda*sinpsi*sinpsi -6*sqrt3*coslambda*cospsi*cospsi*sinlambda + 6*sqrt3*coslambda*sinlambda*sinpsi*sinpsi -16*coslambda*cospsi*sinbeta*sinlambda*sinpsi -8*sqrt3*cospsi*sinbeta*sinlambda*sinlambda*sinpsi -2*sqrt3*coslambda*cospsi*cospsi*sinbeta*sinbeta*sinlambda -2*sqrt3*cosbeta*cosbeta*coslambda*sinlambda*sinpsi*sinpsi + 2*sqrt3*coslambda*sinbeta*sinbeta*sinlambda*sinpsi*sinpsi + 2*sqrt3*cosbeta*cosbeta*coslambda*cospsi*cospsi*sinlambda + 8*sqrt3*coslambda*coslambda*cospsi*sinbeta*sinpsi);
  /**/
  coeffn3Hn3plussin[0] = -1./16*cosbeta * (-9*coslambda*sinbeta*sinpsi*sinpsi + 9*coslambda*cospsi*cospsi*sinbeta + 18*cospsi*sinlambda*sinpsi + sqrt3*sinbeta*sinlambda*sinpsi*sinpsi -sqrt3*cospsi*cospsi*sinbeta*sinlambda + 2*sqrt3*coslambda*cospsi*sinpsi);
  /**/
  coeffn3Hn3plussin[1] = 3./64 * (-3*sqrt3*sinpsi*sinpsi + 3*sqrt3*cospsi*cospsi -12*coslambda*sinlambda*sinpsi*sinpsi -3*sqrt3*cosbeta*cosbeta*sinpsi*sinpsi -3*sqrt3*cospsi*cospsi*sinbeta*sinbeta + 3*sqrt3*cosbeta*cosbeta*cospsi*cospsi + 3*sqrt3*sinbeta*sinbeta*sinpsi*sinpsi + 12*coslambda*cospsi*cospsi*sinlambda -16*coslambda*coslambda*cospsi*sinbeta*sinpsi -4*coslambda*sinbeta*sinbeta*sinlambda*sinpsi*sinpsi -4*cosbeta*cosbeta*coslambda*cospsi*cospsi*sinlambda + 4*coslambda*cospsi*cospsi*sinbeta*sinbeta*sinlambda + 4*cosbeta*cosbeta*coslambda*sinlambda*sinpsi*sinpsi + 16*cospsi*sinbeta*sinlambda*sinlambda*sinpsi);
  /**/
  coeffn3Hn3plussin[2] = 1./16*cosbeta * (-3*coslambda*sinbeta*sinpsi*sinpsi + 3*coslambda*cospsi*cospsi*sinbeta + 6*cospsi*sinlambda*sinpsi + sqrt3*sinbeta*sinlambda*sinpsi*sinpsi -sqrt3*cospsi*cospsi*sinbeta*sinlambda + 2*sqrt3*coslambda*cospsi*sinpsi);
  /**/
  coeffn3Hn3plussin[3] = 1./128 * (-6*coslambda*cospsi*cospsi*sinlambda -3*sqrt3*coslambda*coslambda*sinpsi*sinpsi -3*sqrt3*cospsi*cospsi*sinlambda*sinlambda + 3*sqrt3*coslambda*coslambda*cospsi*cospsi + 3*sqrt3*sinlambda*sinlambda*sinpsi*sinpsi + 6*coslambda*sinlambda*sinpsi*sinpsi + sqrt3*cosbeta*cosbeta*coslambda*coslambda*sinpsi*sinpsi + sqrt3*cosbeta*cosbeta*cospsi*cospsi*sinlambda*sinlambda + sqrt3*coslambda*coslambda*cospsi*cospsi*sinbeta*sinbeta + sqrt3*sinbeta*sinbeta*sinlambda*sinlambda*sinpsi*sinpsi -8*cospsi*sinbeta*sinlambda*sinlambda*sinpsi -2*coslambda*cospsi*cospsi*sinbeta*sinbeta*sinlambda -2*cosbeta*cosbeta*coslambda*sinlambda*sinpsi*sinpsi -sqrt3*cosbeta*cosbeta*coslambda*coslambda*cospsi*cospsi -sqrt3*cosbeta*cosbeta*sinlambda*sinlambda*sinpsi*sinpsi -sqrt3*coslambda*coslambda*sinbeta*sinbeta*sinpsi*sinpsi -sqrt3*cospsi*cospsi*sinbeta*sinbeta*sinlambda*sinlambda + 2*coslambda*sinbeta*sinbeta*sinlambda*sinpsi*sinpsi + 2*cosbeta*cosbeta*coslambda*cospsi*cospsi*sinlambda + 8*coslambda*coslambda*cospsi*sinbeta*sinpsi + 16*sqrt3*coslambda*cospsi*sinbeta*sinlambda*sinpsi);

  /* Projection coefficients for hcross in n3.H.n3 */
  /**/
  coeffn3Hn3crossconst = 1./64 * (4*cospsi*sinpsi -27*cospsi*sinlambda*sinlambda*sinpsi -4*cospsi*sinbeta*sinbeta*sinpsi + 4*cosbeta*cosbeta*cospsi*sinpsi + 27*coslambda*coslambda*cospsi*sinpsi -36*coslambda*cospsi*cospsi*sinbeta*sinlambda -18*sqrt3*coslambda*coslambda*cospsi*cospsi*sinbeta -18*sqrt3*sinbeta*sinlambda*sinlambda*sinpsi*sinpsi -9*cospsi*sinbeta*sinbeta*sinlambda*sinlambda*sinpsi -9*cosbeta*cosbeta*coslambda*coslambda*cospsi*sinpsi + 9*cosbeta*cosbeta*cospsi*sinlambda*sinlambda*sinpsi + 9*coslambda*coslambda*cospsi*sinbeta*sinbeta*sinpsi + 18*sqrt3*coslambda*coslambda*sinbeta*sinpsi*sinpsi + 18*sqrt3*cospsi*cospsi*sinbeta*sinlambda*sinlambda + 36*coslambda*sinbeta*sinlambda*sinpsi*sinpsi -54*sqrt3*coslambda*cospsi*sinlambda*sinpsi -18*sqrt3*coslambda*cospsi*sinbeta*sinbeta*sinlambda*sinpsi + 18*sqrt3*cosbeta*cosbeta*coslambda*cospsi*sinlambda*sinpsi);
  /**/
  coeffn3Hn3crosscos[0] = 1./16*cosbeta * (-9*coslambda*sinpsi*sinpsi + 9*coslambda*cospsi*cospsi -7*sqrt3*sinlambda*sinpsi*sinpsi + 7*sqrt3*cospsi*cospsi*sinlambda + 18*cospsi*sinbeta*sinlambda*sinpsi -14*sqrt3*coslambda*cospsi*sinbeta*sinpsi);
  /**/
  coeffn3Hn3crosscos[1] = -3./32 * (-3*cospsi*sinpsi -6*cospsi*sinlambda*sinlambda*sinpsi -3*cosbeta*cosbeta*cospsi*sinpsi + 3*cospsi*sinbeta*sinbeta*sinpsi + 6*coslambda*coslambda*cospsi*sinpsi -8*coslambda*cospsi*cospsi*sinbeta*sinlambda -2*cospsi*sinbeta*sinbeta*sinlambda*sinlambda*sinpsi -2*cosbeta*cosbeta*coslambda*coslambda*cospsi*sinpsi + 2*cosbeta*cosbeta*cospsi*sinlambda*sinlambda*sinpsi + 2*coslambda*coslambda*cospsi*sinbeta*sinbeta*sinpsi + 8*coslambda*sinbeta*sinlambda*sinpsi*sinpsi);
  /**/
  coeffn3Hn3crosscos[2] = 1./16*cosbeta * (-3*coslambda*sinpsi*sinpsi + 3*coslambda*cospsi*cospsi + sqrt3*sinlambda*sinpsi*sinpsi -sqrt3*cospsi*cospsi*sinlambda + 6*cospsi*sinbeta*sinlambda*sinpsi + 2*sqrt3*coslambda*cospsi*sinbeta*sinpsi);
  /**/
  coeffn3Hn3crosscos[3] = 1./64 * (-3*cospsi*sinlambda*sinlambda*sinpsi + 3*coslambda*coslambda*cospsi*sinpsi + cosbeta*cosbeta*cospsi*sinlambda*sinlambda*sinpsi + coslambda*coslambda*cospsi*sinbeta*sinbeta*sinpsi -4*coslambda*cospsi*cospsi*sinbeta*sinlambda -2*sqrt3*coslambda*coslambda*sinbeta*sinpsi*sinpsi -2*sqrt3*cospsi*cospsi*sinbeta*sinlambda*sinlambda -cospsi*sinbeta*sinbeta*sinlambda*sinlambda*sinpsi -cosbeta*cosbeta*coslambda*coslambda*cospsi*sinpsi + 2*sqrt3*coslambda*coslambda*cospsi*cospsi*sinbeta + 2*sqrt3*sinbeta*sinlambda*sinlambda*sinpsi*sinpsi + 4*coslambda*sinbeta*sinlambda*sinpsi*sinpsi + 6*sqrt3*coslambda*cospsi*sinlambda*sinpsi -2*sqrt3*cosbeta*cosbeta*coslambda*cospsi*sinlambda*sinpsi + 2*sqrt3*coslambda*cospsi*sinbeta*sinbeta*sinlambda*sinpsi);
  /**/
  coeffn3Hn3crosssin[0] = -1./16*cosbeta * (-9*sinlambda*sinpsi*sinpsi + 9*cospsi*cospsi*sinlambda + sqrt3*coslambda*cospsi*cospsi -sqrt3*coslambda*sinpsi*sinpsi -18*coslambda*cospsi*sinbeta*sinpsi + 2*sqrt3*cospsi*sinbeta*sinlambda*sinpsi);
  /**/
  coeffn3Hn3crosssin[1] = -3./32 * (-4*coslambda*coslambda*sinbeta*sinpsi*sinpsi -4*cospsi*cospsi*sinbeta*sinlambda*sinlambda + 3*sqrt3*cospsi*sinpsi + 4*coslambda*coslambda*cospsi*cospsi*sinbeta + 4*sinbeta*sinlambda*sinlambda*sinpsi*sinpsi -3*sqrt3*cospsi*sinbeta*sinbeta*sinpsi + 3*sqrt3*cosbeta*cosbeta*cospsi*sinpsi + 12*coslambda*cospsi*sinlambda*sinpsi -4*cosbeta*cosbeta*coslambda*cospsi*sinlambda*sinpsi + 4*coslambda*cospsi*sinbeta*sinbeta*sinlambda*sinpsi);
  /**/
  coeffn3Hn3crosssin[2] = 1./16*cosbeta * (-3*sinlambda*sinpsi*sinpsi + 3*cospsi*cospsi*sinlambda + sqrt3*coslambda*cospsi*cospsi -sqrt3*coslambda*sinpsi*sinpsi -6*coslambda*cospsi*sinbeta*sinpsi + 2*sqrt3*cospsi*sinbeta*sinlambda*sinpsi);
  /**/
  coeffn3Hn3crosssin[3] = 1./64 * (-2*coslambda*coslambda*sinbeta*sinpsi*sinpsi -2*cospsi*cospsi*sinbeta*sinlambda*sinlambda + 2*coslambda*coslambda*cospsi*cospsi*sinbeta + 2*sinbeta*sinlambda*sinlambda*sinpsi*sinpsi -3*sqrt3*coslambda*coslambda*cospsi*sinpsi + 3*sqrt3*cospsi*sinlambda*sinlambda*sinpsi + 6*coslambda*cospsi*sinlambda*sinpsi + sqrt3*cospsi*sinbeta*sinbeta*sinlambda*sinlambda*sinpsi + sqrt3*cosbeta*cosbeta*coslambda*coslambda*cospsi*sinpsi -4*sqrt3*coslambda*sinbeta*sinlambda*sinpsi*sinpsi -2*cosbeta*cosbeta*coslambda*cospsi*sinlambda*sinpsi -sqrt3*cosbeta*cosbeta*cospsi*sinlambda*sinlambda*sinpsi -sqrt3*coslambda*coslambda*cospsi*sinbeta*sinbeta*sinpsi + 2*coslambda*cospsi*sinbeta*sinbeta*sinlambda*sinpsi + 4*sqrt3*coslambda*cospsi*cospsi*sinbeta*sinlambda);

  /* Projection coefficients for hplus in n2.H.n2 */
  /**/
  coeffn2Hn2plusconst = 1./128 * (-4*cospsi*cospsi + 4*sinpsi*sinpsi -27*coslambda*coslambda*cospsi*cospsi -27*sinlambda*sinlambda*sinpsi*sinpsi -4*cosbeta*cosbeta*cospsi*cospsi -4*sinbeta*sinbeta*sinpsi*sinpsi + 4*cosbeta*cosbeta*sinpsi*sinpsi + 4*cospsi*cospsi*sinbeta*sinbeta + 27*coslambda*coslambda*sinpsi*sinpsi + 27*cospsi*cospsi*sinlambda*sinlambda -9*cosbeta*cosbeta*coslambda*coslambda*sinpsi*sinpsi -9*cosbeta*cosbeta*cospsi*cospsi*sinlambda*sinlambda -9*coslambda*coslambda*cospsi*cospsi*sinbeta*sinbeta -9*sinbeta*sinbeta*sinlambda*sinlambda*sinpsi*sinpsi + 9*cosbeta*cosbeta*coslambda*coslambda*cospsi*cospsi + 9*cosbeta*cosbeta*sinlambda*sinlambda*sinpsi*sinpsi + 9*coslambda*coslambda*sinbeta*sinbeta*sinpsi*sinpsi + 9*cospsi*cospsi*sinbeta*sinbeta*sinlambda*sinlambda -54*sqrt3*coslambda*cospsi*cospsi*sinlambda + 54*sqrt3*coslambda*sinlambda*sinpsi*sinpsi -144*coslambda*cospsi*sinbeta*sinlambda*sinpsi -72*sqrt3*cospsi*sinbeta*sinlambda*sinlambda*sinpsi -18*sqrt3*coslambda*cospsi*cospsi*sinbeta*sinbeta*sinlambda -18*sqrt3*cosbeta*cosbeta*coslambda*sinlambda*sinpsi*sinpsi + 18*sqrt3*coslambda*sinbeta*sinbeta*sinlambda*sinpsi*sinpsi + 18*sqrt3*cosbeta*cosbeta*coslambda*cospsi*cospsi*sinlambda + 72*sqrt3*coslambda*coslambda*cospsi*sinbeta*sinpsi);
  /**/
  coeffn2Hn2pluscos[0] = 1./16*cosbeta * (-18*coslambda*cospsi*sinpsi -9*sinbeta*sinlambda*sinpsi*sinpsi + 9*cospsi*cospsi*sinbeta*sinlambda -7*sqrt3*coslambda*sinbeta*sinpsi*sinpsi + 7*sqrt3*coslambda*cospsi*cospsi*sinbeta + 14*sqrt3*cospsi*sinlambda*sinpsi);
  /**/
  coeffn2Hn2pluscos[1] = -3./64 * (-3*sinpsi*sinpsi + 3*cospsi*cospsi -6*coslambda*coslambda*cospsi*cospsi -6*sinlambda*sinlambda*sinpsi*sinpsi -3*cosbeta*cosbeta*sinpsi*sinpsi -3*cospsi*cospsi*sinbeta*sinbeta + 3*cosbeta*cosbeta*cospsi*cospsi + 3*sinbeta*sinbeta*sinpsi*sinpsi + 6*coslambda*coslambda*sinpsi*sinpsi + 6*cospsi*cospsi*sinlambda*sinlambda -2*cosbeta*cosbeta*coslambda*coslambda*sinpsi*sinpsi -2*cosbeta*cosbeta*cospsi*cospsi*sinlambda*sinlambda -2*coslambda*coslambda*cospsi*cospsi*sinbeta*sinbeta -2*sinbeta*sinbeta*sinlambda*sinlambda*sinpsi*sinpsi + 2*cosbeta*cosbeta*coslambda*coslambda*cospsi*cospsi + 2*cosbeta*cosbeta*sinlambda*sinlambda*sinpsi*sinpsi + 2*coslambda*coslambda*sinbeta*sinbeta*sinpsi*sinpsi + 2*cospsi*cospsi*sinbeta*sinbeta*sinlambda*sinlambda -32*coslambda*cospsi*sinbeta*sinlambda*sinpsi);
  /**/
  coeffn2Hn2pluscos[2] = -1./16*cosbeta * (-3*cospsi*cospsi*sinbeta*sinlambda + 3*sinbeta*sinlambda*sinpsi*sinpsi + 6*coslambda*cospsi*sinpsi + sqrt3*coslambda*cospsi*cospsi*sinbeta -sqrt3*coslambda*sinbeta*sinpsi*sinpsi + 2*sqrt3*cospsi*sinlambda*sinpsi);
  /**/
  coeffn2Hn2pluscos[3] = 1./128 * (-3*coslambda*coslambda*cospsi*cospsi -3*sinlambda*sinlambda*sinpsi*sinpsi + 3*coslambda*coslambda*sinpsi*sinpsi + 3*cospsi*cospsi*sinlambda*sinlambda + cosbeta*cosbeta*coslambda*coslambda*cospsi*cospsi + cosbeta*cosbeta*sinlambda*sinlambda*sinpsi*sinpsi + coslambda*coslambda*sinbeta*sinbeta*sinpsi*sinpsi + cospsi*cospsi*sinbeta*sinbeta*sinlambda*sinlambda -cosbeta*cosbeta*coslambda*coslambda*sinpsi*sinpsi -cosbeta*cosbeta*cospsi*cospsi*sinlambda*sinlambda -coslambda*coslambda*cospsi*cospsi*sinbeta*sinbeta -sinbeta*sinbeta*sinlambda*sinlambda*sinpsi*sinpsi -6*sqrt3*coslambda*sinlambda*sinpsi*sinpsi + 6*sqrt3*coslambda*cospsi*cospsi*sinlambda -16*coslambda*cospsi*sinbeta*sinlambda*sinpsi -8*sqrt3*coslambda*coslambda*cospsi*sinbeta*sinpsi -2*sqrt3*coslambda*sinbeta*sinbeta*sinlambda*sinpsi*sinpsi -2*sqrt3*cosbeta*cosbeta*coslambda*cospsi*cospsi*sinlambda + 2*sqrt3*coslambda*cospsi*cospsi*sinbeta*sinbeta*sinlambda + 2*sqrt3*cosbeta*cosbeta*coslambda*sinlambda*sinpsi*sinpsi + 8*sqrt3*cospsi*sinbeta*sinlambda*sinlambda*sinpsi);
  /**/
  coeffn2Hn2plussin[0] = 1./16*cosbeta * (-9*coslambda*sinbeta*sinpsi*sinpsi + 9*coslambda*cospsi*cospsi*sinbeta + 18*cospsi*sinlambda*sinpsi + sqrt3*cospsi*cospsi*sinbeta*sinlambda -2*sqrt3*coslambda*cospsi*sinpsi -sqrt3*sinbeta*sinlambda*sinpsi*sinpsi);
  /**/
  coeffn2Hn2plussin[1] = -3./64 * (-3*sqrt3*sinpsi*sinpsi + 3*sqrt3*cospsi*cospsi -12*coslambda*cospsi*cospsi*sinlambda -3*sqrt3*cosbeta*cosbeta*sinpsi*sinpsi -3*sqrt3*cospsi*cospsi*sinbeta*sinbeta + 3*sqrt3*cosbeta*cosbeta*cospsi*cospsi + 3*sqrt3*sinbeta*sinbeta*sinpsi*sinpsi + 12*coslambda*sinlambda*sinpsi*sinpsi -16*cospsi*sinbeta*sinlambda*sinlambda*sinpsi -4*coslambda*cospsi*cospsi*sinbeta*sinbeta*sinlambda -4*cosbeta*cosbeta*coslambda*sinlambda*sinpsi*sinpsi + 4*coslambda*sinbeta*sinbeta*sinlambda*sinpsi*sinpsi + 4*cosbeta*cosbeta*coslambda*cospsi*cospsi*sinlambda + 16*coslambda*coslambda*cospsi*sinbeta*sinpsi);
  /**/
  coeffn2Hn2plussin[2] = -1./16*cosbeta * (-3*coslambda*sinbeta*sinpsi*sinpsi + 3*coslambda*cospsi*cospsi*sinbeta + 6*cospsi*sinlambda*sinpsi + sqrt3*cospsi*cospsi*sinbeta*sinlambda -2*sqrt3*coslambda*cospsi*sinpsi -sqrt3*sinbeta*sinlambda*sinpsi*sinpsi);
  /**/
  coeffn2Hn2plussin[3] = 1./128 * (-6*coslambda*cospsi*cospsi*sinlambda -3*sqrt3*coslambda*coslambda*cospsi*cospsi -3*sqrt3*sinlambda*sinlambda*sinpsi*sinpsi + 3*sqrt3*coslambda*coslambda*sinpsi*sinpsi + 3*sqrt3*cospsi*cospsi*sinlambda*sinlambda + 6*coslambda*sinlambda*sinpsi*sinpsi + sqrt3*cosbeta*cosbeta*coslambda*coslambda*cospsi*cospsi + sqrt3*cosbeta*cosbeta*sinlambda*sinlambda*sinpsi*sinpsi + sqrt3*coslambda*coslambda*sinbeta*sinbeta*sinpsi*sinpsi + sqrt3*cospsi*cospsi*sinbeta*sinbeta*sinlambda*sinlambda -8*cospsi*sinbeta*sinlambda*sinlambda*sinpsi -2*coslambda*cospsi*cospsi*sinbeta*sinbeta*sinlambda -2*cosbeta*cosbeta*coslambda*sinlambda*sinpsi*sinpsi -sqrt3*cosbeta*cosbeta*coslambda*coslambda*sinpsi*sinpsi -sqrt3*cosbeta*cosbeta*cospsi*cospsi*sinlambda*sinlambda -sqrt3*coslambda*coslambda*cospsi*cospsi*sinbeta*sinbeta -sqrt3*sinbeta*sinbeta*sinlambda*sinlambda*sinpsi*sinpsi + 2*coslambda*sinbeta*sinbeta*sinlambda*sinpsi*sinpsi + 2*cosbeta*cosbeta*coslambda*cospsi*cospsi*sinlambda + 8*coslambda*coslambda*cospsi*sinbeta*sinpsi -16*sqrt3*coslambda*cospsi*sinbeta*sinlambda*sinpsi);

  /* Projection coefficients for hcross in n2.H.n2 */
  /**/
  coeffn2Hn2crossconst = 1./64 * (4*cospsi*sinpsi -27*cospsi*sinlambda*sinlambda*sinpsi -4*cospsi*sinbeta*sinbeta*sinpsi + 4*cosbeta*cosbeta*cospsi*sinpsi + 27*coslambda*coslambda*cospsi*sinpsi -36*coslambda*cospsi*cospsi*sinbeta*sinlambda -18*sqrt3*coslambda*coslambda*sinbeta*sinpsi*sinpsi -18*sqrt3*cospsi*cospsi*sinbeta*sinlambda*sinlambda -9*cospsi*sinbeta*sinbeta*sinlambda*sinlambda*sinpsi -9*cosbeta*cosbeta*coslambda*coslambda*cospsi*sinpsi + 9*cosbeta*cosbeta*cospsi*sinlambda*sinlambda*sinpsi + 9*coslambda*coslambda*cospsi*sinbeta*sinbeta*sinpsi + 18*sqrt3*coslambda*coslambda*cospsi*cospsi*sinbeta + 18*sqrt3*sinbeta*sinlambda*sinlambda*sinpsi*sinpsi + 36*coslambda*sinbeta*sinlambda*sinpsi*sinpsi + 54*sqrt3*coslambda*cospsi*sinlambda*sinpsi -18*sqrt3*cosbeta*cosbeta*coslambda*cospsi*sinlambda*sinpsi + 18*sqrt3*coslambda*cospsi*sinbeta*sinbeta*sinlambda*sinpsi);
  /**/
  coeffn2Hn2crosscos[0] = -1./16*cosbeta * (-9*coslambda*sinpsi*sinpsi + 9*coslambda*cospsi*cospsi -7*sqrt3*cospsi*cospsi*sinlambda + 7*sqrt3*sinlambda*sinpsi*sinpsi + 18*cospsi*sinbeta*sinlambda*sinpsi + 14*sqrt3*coslambda*cospsi*sinbeta*sinpsi);
  /**/
  coeffn2Hn2crosscos[1] = -3./32 * (-3*cospsi*sinpsi -6*cospsi*sinlambda*sinlambda*sinpsi -3*cosbeta*cosbeta*cospsi*sinpsi + 3*cospsi*sinbeta*sinbeta*sinpsi + 6*coslambda*coslambda*cospsi*sinpsi -8*coslambda*cospsi*cospsi*sinbeta*sinlambda -2*cospsi*sinbeta*sinbeta*sinlambda*sinlambda*sinpsi -2*cosbeta*cosbeta*coslambda*coslambda*cospsi*sinpsi + 2*cosbeta*cosbeta*cospsi*sinlambda*sinlambda*sinpsi + 2*coslambda*coslambda*cospsi*sinbeta*sinbeta*sinpsi + 8*coslambda*sinbeta*sinlambda*sinpsi*sinpsi);
  /**/
  coeffn2Hn2crosscos[2] = -1./16*cosbeta * (-3*coslambda*sinpsi*sinpsi + 3*coslambda*cospsi*cospsi + sqrt3*cospsi*cospsi*sinlambda -sqrt3*sinlambda*sinpsi*sinpsi + 6*cospsi*sinbeta*sinlambda*sinpsi -2*sqrt3*coslambda*cospsi*sinbeta*sinpsi);
  /**/
  coeffn2Hn2crosscos[3] = 1./64 * (-3*cospsi*sinlambda*sinlambda*sinpsi + 3*coslambda*coslambda*cospsi*sinpsi + cosbeta*cosbeta*cospsi*sinlambda*sinlambda*sinpsi + coslambda*coslambda*cospsi*sinbeta*sinbeta*sinpsi -4*coslambda*cospsi*cospsi*sinbeta*sinlambda -2*sqrt3*coslambda*coslambda*cospsi*cospsi*sinbeta -2*sqrt3*sinbeta*sinlambda*sinlambda*sinpsi*sinpsi -cospsi*sinbeta*sinbeta*sinlambda*sinlambda*sinpsi -cosbeta*cosbeta*coslambda*coslambda*cospsi*sinpsi + 2*sqrt3*coslambda*coslambda*sinbeta*sinpsi*sinpsi + 2*sqrt3*cospsi*cospsi*sinbeta*sinlambda*sinlambda + 4*coslambda*sinbeta*sinlambda*sinpsi*sinpsi -6*sqrt3*coslambda*cospsi*sinlambda*sinpsi -2*sqrt3*coslambda*cospsi*sinbeta*sinbeta*sinlambda*sinpsi + 2*sqrt3*cosbeta*cosbeta*coslambda*cospsi*sinlambda*sinpsi);
  /**/
  coeffn2Hn2crosssin[0] = -1./16*cosbeta * (-9*cospsi*cospsi*sinlambda + 9*sinlambda*sinpsi*sinpsi + sqrt3*coslambda*cospsi*cospsi -sqrt3*coslambda*sinpsi*sinpsi + 18*coslambda*cospsi*sinbeta*sinpsi + 2*sqrt3*cospsi*sinbeta*sinlambda*sinpsi);
  /**/
  coeffn2Hn2crosssin[1] = -3./32 * (-4*coslambda*coslambda*sinbeta*sinpsi*sinpsi -4*cospsi*cospsi*sinbeta*sinlambda*sinlambda -3*sqrt3*cospsi*sinpsi + 4*coslambda*coslambda*cospsi*cospsi*sinbeta + 4*sinbeta*sinlambda*sinlambda*sinpsi*sinpsi -3*sqrt3*cosbeta*cosbeta*cospsi*sinpsi + 3*sqrt3*cospsi*sinbeta*sinbeta*sinpsi + 12*coslambda*cospsi*sinlambda*sinpsi -4*cosbeta*cosbeta*coslambda*cospsi*sinlambda*sinpsi + 4*coslambda*cospsi*sinbeta*sinbeta*sinlambda*sinpsi);
  /**/
  coeffn2Hn2crosssin[2] = 1./16*cosbeta * (-3*cospsi*cospsi*sinlambda + 3*sinlambda*sinpsi*sinpsi + sqrt3*coslambda*cospsi*cospsi -sqrt3*coslambda*sinpsi*sinpsi + 6*coslambda*cospsi*sinbeta*sinpsi + 2*sqrt3*cospsi*sinbeta*sinlambda*sinpsi);
  /**/
  coeffn2Hn2crosssin[3] = 1./64 * (-2*coslambda*coslambda*sinbeta*sinpsi*sinpsi -2*cospsi*cospsi*sinbeta*sinlambda*sinlambda + 2*coslambda*coslambda*cospsi*cospsi*sinbeta + 2*sinbeta*sinlambda*sinlambda*sinpsi*sinpsi -3*sqrt3*cospsi*sinlambda*sinlambda*sinpsi + 3*sqrt3*coslambda*coslambda*cospsi*sinpsi + 6*coslambda*cospsi*sinlambda*sinpsi + sqrt3*cosbeta*cosbeta*cospsi*sinlambda*sinlambda*sinpsi + sqrt3*coslambda*coslambda*cospsi*sinbeta*sinbeta*sinpsi -4*sqrt3*coslambda*cospsi*cospsi*sinbeta*sinlambda -2*cosbeta*cosbeta*coslambda*cospsi*sinlambda*sinpsi -sqrt3*cospsi*sinbeta*sinbeta*sinlambda*sinlambda*sinpsi -sqrt3*cosbeta*cosbeta*coslambda*coslambda*cospsi*sinpsi + 2*coslambda*cospsi*sinbeta*sinbeta*sinlambda*sinpsi + 4*sqrt3*coslambda*sinbeta*sinlambda*sinpsi*sinpsi);

  /* Projection coefficients for hplus in n1.H.n1 */
  /**/
  coeffn1Hn1plusconst = 1./64 * (-2*cospsi*cospsi + 2*sinpsi*sinpsi -27*coslambda*coslambda*sinpsi*sinpsi -27*cospsi*cospsi*sinlambda*sinlambda -2*cosbeta*cosbeta*cospsi*cospsi -2*sinbeta*sinbeta*sinpsi*sinpsi + 2*cosbeta*cosbeta*sinpsi*sinpsi + 2*cospsi*cospsi*sinbeta*sinbeta + 27*coslambda*coslambda*cospsi*cospsi + 27*sinlambda*sinlambda*sinpsi*sinpsi -9*cosbeta*cosbeta*coslambda*coslambda*cospsi*cospsi -9*cosbeta*cosbeta*sinlambda*sinlambda*sinpsi*sinpsi -9*coslambda*coslambda*sinbeta*sinbeta*sinpsi*sinpsi -9*cospsi*cospsi*sinbeta*sinbeta*sinlambda*sinlambda + 9*cosbeta*cosbeta*coslambda*coslambda*sinpsi*sinpsi + 9*cosbeta*cosbeta*cospsi*cospsi*sinlambda*sinlambda + 9*coslambda*coslambda*cospsi*cospsi*sinbeta*sinbeta + 9*sinbeta*sinbeta*sinlambda*sinlambda*sinpsi*sinpsi + 144*coslambda*cospsi*sinbeta*sinlambda*sinpsi);
  /**/
  coeffn1Hn1pluscos[0] = -1./8*sqrt3*cosbeta * (coslambda*cospsi*cospsi*sinbeta -coslambda*sinbeta*sinpsi*sinpsi + 2*cospsi*sinlambda*sinpsi);
  /**/
  coeffn1Hn1pluscos[1] = -3./32 * (-3*cospsi*cospsi + 3*sinpsi*sinpsi -3*cosbeta*cosbeta*cospsi*cospsi -3*coslambda*coslambda*cospsi*cospsi -3*sinbeta*sinbeta*sinpsi*sinpsi -3*sinlambda*sinlambda*sinpsi*sinpsi + 3*cosbeta*cosbeta*sinpsi*sinpsi + 3*coslambda*coslambda*sinpsi*sinpsi + 3*cospsi*cospsi*sinbeta*sinbeta + 3*cospsi*cospsi*sinlambda*sinlambda + cosbeta*cosbeta*coslambda*coslambda*cospsi*cospsi + cosbeta*cosbeta*sinlambda*sinlambda*sinpsi*sinpsi + coslambda*coslambda*sinbeta*sinbeta*sinpsi*sinpsi + cospsi*cospsi*sinbeta*sinbeta*sinlambda*sinlambda -cosbeta*cosbeta*coslambda*coslambda*sinpsi*sinpsi -cosbeta*cosbeta*cospsi*cospsi*sinlambda*sinlambda -coslambda*coslambda*cospsi*cospsi*sinbeta*sinbeta -sinbeta*sinbeta*sinlambda*sinlambda*sinpsi*sinpsi -16*coslambda*cospsi*sinbeta*sinlambda*sinpsi);
  /**/
  coeffn1Hn1pluscos[2] = 1./8*sqrt3*cosbeta * (coslambda*cospsi*cospsi*sinbeta -coslambda*sinbeta*sinpsi*sinpsi + 2*cospsi*sinlambda*sinpsi);
  /**/
  coeffn1Hn1pluscos[3] = 1./64 * (-3*coslambda*coslambda*sinpsi*sinpsi -3*cospsi*cospsi*sinlambda*sinlambda + 3*coslambda*coslambda*cospsi*cospsi + 3*sinlambda*sinlambda*sinpsi*sinpsi + cosbeta*cosbeta*coslambda*coslambda*sinpsi*sinpsi + cosbeta*cosbeta*cospsi*cospsi*sinlambda*sinlambda + coslambda*coslambda*cospsi*cospsi*sinbeta*sinbeta + sinbeta*sinbeta*sinlambda*sinlambda*sinpsi*sinpsi -cosbeta*cosbeta*coslambda*coslambda*cospsi*cospsi -cosbeta*cosbeta*sinlambda*sinlambda*sinpsi*sinpsi -coslambda*coslambda*sinbeta*sinbeta*sinpsi*sinpsi -cospsi*cospsi*sinbeta*sinbeta*sinlambda*sinlambda + 16*coslambda*cospsi*sinbeta*sinlambda*sinpsi);
  /**/
  coeffn1Hn1plussin[0] = 5./8*sqrt3*cosbeta * (cospsi*cospsi*sinbeta*sinlambda -2*coslambda*cospsi*sinpsi -sinbeta*sinlambda*sinpsi*sinpsi);
  /**/
  coeffn1Hn1plussin[1] = -3./16 * (-3*coslambda*cospsi*cospsi*sinlambda + 3*coslambda*sinlambda*sinpsi*sinpsi + coslambda*sinbeta*sinbeta*sinlambda*sinpsi*sinpsi + cosbeta*cosbeta*coslambda*cospsi*cospsi*sinlambda -4*cospsi*sinbeta*sinlambda*sinlambda*sinpsi -coslambda*cospsi*cospsi*sinbeta*sinbeta*sinlambda -cosbeta*cosbeta*coslambda*sinlambda*sinpsi*sinpsi + 4*coslambda*coslambda*cospsi*sinbeta*sinpsi);
  /**/
  coeffn1Hn1plussin[2] = 1./8*sqrt3*cosbeta * (cospsi*cospsi*sinbeta*sinlambda -2*coslambda*cospsi*sinpsi -sinbeta*sinlambda*sinpsi*sinpsi);
  /**/
  coeffn1Hn1plussin[3] = 1./32 * (-3*coslambda*sinlambda*sinpsi*sinpsi + 3*coslambda*cospsi*cospsi*sinlambda + coslambda*cospsi*cospsi*sinbeta*sinbeta*sinlambda + cosbeta*cosbeta*coslambda*sinlambda*sinpsi*sinpsi -4*coslambda*coslambda*cospsi*sinbeta*sinpsi -coslambda*sinbeta*sinbeta*sinlambda*sinpsi*sinpsi -cosbeta*cosbeta*coslambda*cospsi*cospsi*sinlambda + 4*cospsi*sinbeta*sinlambda*sinlambda*sinpsi);

  /* Projection coefficients for hcross in n1.H.n1 */
  /**/
  coeffn1Hn1crossconst = 1./32 * (2*cospsi*sinpsi -27*coslambda*coslambda*cospsi*sinpsi -2*cospsi*sinbeta*sinbeta*sinpsi + 2*cosbeta*cosbeta*cospsi*sinpsi + 27*cospsi*sinlambda*sinlambda*sinpsi -36*coslambda*sinbeta*sinlambda*sinpsi*sinpsi -9*cosbeta*cosbeta*cospsi*sinlambda*sinlambda*sinpsi -9*coslambda*coslambda*cospsi*sinbeta*sinbeta*sinpsi + 9*cospsi*sinbeta*sinbeta*sinlambda*sinlambda*sinpsi + 9*cosbeta*cosbeta*coslambda*coslambda*cospsi*sinpsi + 36*coslambda*cospsi*cospsi*sinbeta*sinlambda);
  /**/
  coeffn1Hn1crosscos[0] = -1./8*sqrt3*cosbeta * (cospsi*cospsi*sinlambda -sinlambda*sinpsi*sinpsi -2*coslambda*cospsi*sinbeta*sinpsi);
  /**/
  coeffn1Hn1crosscos[1] = -3./16 * (3*cospsi*sinpsi -3*cospsi*sinbeta*sinbeta*sinpsi -3*cospsi*sinlambda*sinlambda*sinpsi + 3*cosbeta*cosbeta*cospsi*sinpsi + 3*coslambda*coslambda*cospsi*sinpsi + cosbeta*cosbeta*cospsi*sinlambda*sinlambda*sinpsi + coslambda*coslambda*cospsi*sinbeta*sinbeta*sinpsi -4*coslambda*cospsi*cospsi*sinbeta*sinlambda -cospsi*sinbeta*sinbeta*sinlambda*sinlambda*sinpsi -cosbeta*cosbeta*coslambda*coslambda*cospsi*sinpsi + 4*coslambda*sinbeta*sinlambda*sinpsi*sinpsi);
  /**/
  coeffn1Hn1crosscos[2] = 1./8*sqrt3*cosbeta * (cospsi*cospsi*sinlambda -sinlambda*sinpsi*sinpsi -2*coslambda*cospsi*sinbeta*sinpsi);
  /**/
  coeffn1Hn1crosscos[3] = 1./32 * (-3*coslambda*coslambda*cospsi*sinpsi + 3*cospsi*sinlambda*sinlambda*sinpsi + cospsi*sinbeta*sinbeta*sinlambda*sinlambda*sinpsi + cosbeta*cosbeta*coslambda*coslambda*cospsi*sinpsi -4*coslambda*sinbeta*sinlambda*sinpsi*sinpsi -cosbeta*cosbeta*cospsi*sinlambda*sinlambda*sinpsi -coslambda*coslambda*cospsi*sinbeta*sinbeta*sinpsi + 4*coslambda*cospsi*cospsi*sinbeta*sinlambda);
  /**/
  coeffn1Hn1crosssin[0] = -5./8*sqrt3*cosbeta * (coslambda*cospsi*cospsi -coslambda*sinpsi*sinpsi + 2*cospsi*sinbeta*sinlambda*sinpsi);
  /**/
  coeffn1Hn1crosssin[1] = -3./8 * (coslambda*coslambda*cospsi*cospsi*sinbeta + sinbeta*sinlambda*sinlambda*sinpsi*sinpsi -coslambda*coslambda*sinbeta*sinpsi*sinpsi -cospsi*cospsi*sinbeta*sinlambda*sinlambda + 3*coslambda*cospsi*sinlambda*sinpsi + coslambda*cospsi*sinbeta*sinbeta*sinlambda*sinpsi -cosbeta*cosbeta*coslambda*cospsi*sinlambda*sinpsi);
  /**/
  coeffn1Hn1crosssin[2] = -1./8*sqrt3*cosbeta * (coslambda*cospsi*cospsi -coslambda*sinpsi*sinpsi + 2*cospsi*sinbeta*sinlambda*sinpsi);
  /**/
  coeffn1Hn1crosssin[3] = 1./16 * (coslambda*coslambda*sinbeta*sinpsi*sinpsi + cospsi*cospsi*sinbeta*sinlambda*sinlambda -coslambda*coslambda*cospsi*cospsi*sinbeta -sinbeta*sinlambda*sinlambda*sinpsi*sinpsi -3*coslambda*cospsi*sinlambda*sinpsi + cosbeta*cosbeta*coslambda*cospsi*sinlambda*sinpsi -coslambda*cospsi*sinbeta*sinbeta*sinlambda*sinpsi);

  /* Coefficients in k.n3 */
  /**/
  coeffkn3const = 3./8*cosbeta * (sinlambda -sqrt3*coslambda);
  /**/
  coeffkn3cos[0] = 3./4 * (-sinbeta);
  /**/
  coeffkn3cos[1] = -1./8*cosbeta * (-sinlambda -sqrt3*coslambda);
  /**/
  coeffkn3sin[0] = -1./4*sqrt3 * (-sinbeta);
  /**/
  coeffkn3sin[1] = 1./8*cosbeta * (-coslambda + sqrt3*sinlambda);

  /* Coefficients in k.n2 */
  /**/
  coeffkn2const = -3./8*cosbeta * (-sinlambda -sqrt3*coslambda);
  /**/
  coeffkn2cos[0] = -3./4 * (-sinbeta);
  /**/
  coeffkn2cos[1] = 1./8*cosbeta * (sinlambda -sqrt3*coslambda);
  /**/
  coeffkn2sin[0] = -1./4*sqrt3 * (-sinbeta);
  /**/
  coeffkn2sin[1] = 1./8*cosbeta * (-coslambda -sqrt3*sinlambda);

  /* Coefficients in k.n1 */
  /**/
  coeffkn1const = 3./4*cosbeta * (-sinlambda);
  /**/
  coeffkn1cos[0] =  0. ;
  /**/
  coeffkn1cos[1] = 1./4*cosbeta * (-sinlambda);
  /**/
  coeffkn1sin[0] = 1./2*sqrt3 * (-sinbeta);
  /**/
  coeffkn1sin[1] = -1./4*cosbeta * (-coslambda);

  /* Coefficients in k.(p1+p2) */
  /**/
  coeffkp1plusp2const = -1./8*cosbeta * (-3*sinlambda -sqrt3*coslambda);
  /**/
  coeffkp1plusp2cos[0] = -1./4 * (-sinbeta);
  /**/
  coeffkp1plusp2cos[1] = 1./24*cosbeta * (3*sinlambda -sqrt3*coslambda);
  /**/
  coeffkp1plusp2sin[0] = -1./4*sqrt3 * (-sinbeta);
  /**/
  coeffkp1plusp2sin[1] = 1./24*cosbeta * (-3*coslambda -sqrt3*sinlambda);

  /* Coefficients in k.(p2+p3) */
  /**/
  coeffkp2plusp3const = 1./4*sqrt3*cosbeta * (-coslambda);
  /**/
  coeffkp2plusp3cos[0] = 1./2 * (-sinbeta);
  /**/
  coeffkp2plusp3cos[1] = -1./4/sqrt3 * (-cosbeta*coslambda);
  /**/
  coeffkp2plusp3sin[0] =  0. ;
  /**/
  coeffkp2plusp3sin[1] = -1./4/sqrt3 * (-cosbeta*sinlambda);

  /* Coefficients in k.(p3+p1) */
  /**/
  coeffkp3plusp1const = -1./8*cosbeta * (3*sinlambda -sqrt3*coslambda);
  /**/
  coeffkp3plusp1cos[0] = -1./4 * (-sinbeta);
  /**/
  coeffkp3plusp1cos[1] = 1./24*cosbeta * (-3*sinlambda -sqrt3*coslambda);
  /**/
  coeffkp3plusp1sin[0] = 1./4*sqrt3 * (-sinbeta);
  /**/
  coeffkp3plusp1sin[1] = -1./24*cosbeta * (-3*coslambda + sqrt3*sinlambda);

  /* Coefficients in k.p1 */
  /**/
  coeffkp1const = -1./4*sqrt3 * (-cosbeta*coslambda);
  /**/
  coeffkp1cos[0] = -1./2 * (-sinbeta);
  /**/
  coeffkp1cos[1] = 1./(4*sqrt3) * (-cosbeta*coslambda);
  /**/
  coeffkp1sin[0] =  0. ;
  /**/
  coeffkp1sin[1] = 1./(4*sqrt3) * (-cosbeta*sinlambda);

  /* Coefficients in k.p2 */
  /**/
  coeffkp2const = 1./8*cosbeta * (3*sinlambda -sqrt3*coslambda);
  /**/
  coeffkp2cos[0] = 1./4 * (-sinbeta);
  /**/
  coeffkp2cos[1] = -1./24*cosbeta * (-3*sinlambda -sqrt3*coslambda);
  /**/
  coeffkp2sin[0] = -1./4*sqrt3 * (-sinbeta);
  /**/
  coeffkp2sin[1] = 1./24*cosbeta * (-3*coslambda + sqrt3*sinlambda);

  /* Coefficients in k.p3 */
  /**/
  coeffkp3const = 1./8*cosbeta * (-3*sinlambda -sqrt3*coslambda);
  /**/
  coeffkp3cos[0] = 1./4 * (-sinbeta);
  /**/
  coeffkp3cos[1] = -1./24*cosbeta * (3*sinlambda -sqrt3*coslambda);
  /**/
  coeffkp3sin[0] = 1./4*sqrt3 * (-sinbeta);
  /**/
  coeffkp3sin[1] = -1./24*cosbeta * (-3*coslambda -sqrt3*sinlambda);

  /* Coefficients in k.R */
  /**/
  coeffkRconst = 0.;
  coeffkRcos[0] = 1. * (-cosbeta*coslambda);
  coeffkRsin[0] = 1. * (-cosbeta*sinlambda);
  coeffkRcos[1] = 0.;
  coeffkRsin[1] = 0.;

}

/*********************** Fourier-domain response ************************/

/* Individual functions GABmode: older version, does not include the orbital delay (was treated separately as Bessel phase) */
/* Collective function EvaluateGABmode: orbital delay included */
/* Conventions changed: now MLDC conventions */

/* Function evaluating G21, combining the two polarization with the spherical harmonics factors */
double complex G21mode(const double f, const double t, const double complex Yfactorplus, const double complex Yfactorcross) {

  for(int j=0; j<4; j++) {
    cosarray[j] = cos((j+1) * Omega_SI * t);
    sinarray[j] = sin((j+1) * Omega_SI * t);
  }
  double n3Pn3plus = coeffn3Hn3plusconst;
  double n3Pn3cross = coeffn3Hn3crossconst;
  for(int j=0; j<4; j++) {
    n3Pn3plus += cosarray[j] * coeffn3Hn3pluscos[j] + sinarray[j] * coeffn3Hn3plussin[j];
    n3Pn3cross += cosarray[j] * coeffn3Hn3crosscos[j] + sinarray[j] * coeffn3Hn3crosssin[j];
  }
  double kn3 = coeffkn3const;
  double kp1plusp2 = coeffkp1plusp2const;
  for(int j=0; j<2; j++) {
    kn3 += cosarray[j] * coeffkn3cos[j] + sinarray[j] * coeffkn3sin[j];
    kp1plusp2 += cosarray[j] * coeffkp1plusp2cos[j] + sinarray[j] * coeffkp1plusp2sin[j];
  }
  return I*PI*f*L_SI/C_SI * (n3Pn3plus*Yfactorplus + n3Pn3cross*Yfactorcross) * sinc( PI*f*L_SI/C_SI * (1.+kn3)) * cexp( I*PI*f*L_SI/C_SI * (1.+kp1plusp2) );
}
/* Function evaluating G12, combining the two polarization with the spherical harmonics factors */
double complex G12mode(const double f, const double t, const double complex Yfactorplus, const double complex Yfactorcross) {

  for(int j=0; j<4; j++) {
    cosarray[j] = cos((j+1) * Omega_SI * t);
    sinarray[j] = sin((j+1) * Omega_SI * t);
  }
  double n3Pn3plus = coeffn3Hn3plusconst;
  double n3Pn3cross = coeffn3Hn3crossconst;
  for(int j=0; j<4; j++) {
    n3Pn3plus += cosarray[j] * coeffn3Hn3pluscos[j] + sinarray[j] * coeffn3Hn3plussin[j];
    n3Pn3cross += cosarray[j] * coeffn3Hn3crosscos[j] + sinarray[j] * coeffn3Hn3crosssin[j];
  }
  double kn3 = coeffkn3const;
  double kp1plusp2 = coeffkp1plusp2const;
  for(int j=0; j<2; j++) {
    kn3 += cosarray[j] * coeffkn3cos[j] + sinarray[j] * coeffkn3sin[j];
    kp1plusp2 += cosarray[j] * coeffkp1plusp2cos[j] + sinarray[j] * coeffkp1plusp2sin[j];
  }

  return I*PI*f*L_SI/C_SI * (n3Pn3plus*Yfactorplus + n3Pn3cross*Yfactorcross) * sinc( PI*f*L_SI/C_SI * (1.-kn3)) * cexp( I*PI*f*L_SI/C_SI * (1.+kp1plusp2) );
}
/* Function evaluating G32, combining the two polarization with the spherical harmonics factors */
double complex G32mode(const double f, const double t, const double complex Yfactorplus, const double complex Yfactorcross) {

  for(int j=0; j<4; j++) {
    cosarray[j] = cos((j+1) * Omega_SI * t);
    sinarray[j] = sin((j+1) * Omega_SI * t);
  }
  double n1Pn1plus = coeffn1Hn1plusconst;
  double n1Pn1cross = coeffn1Hn1crossconst;
  for(int j=0; j<4; j++) {
    n1Pn1plus += cosarray[j] * coeffn1Hn1pluscos[j] + sinarray[j] * coeffn1Hn1plussin[j];
    n1Pn1cross += cosarray[j] * coeffn1Hn1crosscos[j] + sinarray[j] * coeffn1Hn1crosssin[j];
  }
  double kn1 = coeffkn1const;
  double kp2plusp3 = coeffkp2plusp3const;
  for(int j=0; j<2; j++) {
    kn1 += cosarray[j] * coeffkn1cos[j] + sinarray[j] * coeffkn1sin[j];
    kp2plusp3 += cosarray[j] * coeffkp2plusp3cos[j] + sinarray[j] * coeffkp2plusp3sin[j];
  }

  return I*PI*f*L_SI/C_SI * (n1Pn1plus*Yfactorplus + n1Pn1cross*Yfactorcross) * sinc( PI*f*L_SI/C_SI * (1.+kn1)) * cexp( I*PI*f*L_SI/C_SI * (1.+kp2plusp3) );
}
/* Function evaluating G23, combining the two polarization with the spherical harmonics factors */
double complex G23mode(const double f, const double t, const double complex Yfactorplus, const double complex Yfactorcross) {

  for(int j=0; j<4; j++) {
    cosarray[j] = cos((j+1) * Omega_SI * t);
    sinarray[j] = sin((j+1) * Omega_SI * t);
  }
  double n1Pn1plus = coeffn1Hn1plusconst;
  double n1Pn1cross = coeffn1Hn1crossconst;
  for(int j=0; j<4; j++) {
    n1Pn1plus += cosarray[j] * coeffn1Hn1pluscos[j] + sinarray[j] * coeffn1Hn1plussin[j];
    n1Pn1cross += cosarray[j] * coeffn1Hn1crosscos[j] + sinarray[j] * coeffn1Hn1crosssin[j];
  }
  double kn1 = coeffkn1const;
  double kp2plusp3 = coeffkp2plusp3const;
  for(int j=0; j<2; j++) {
    kn1 += cosarray[j] * coeffkn1cos[j] + sinarray[j] * coeffkn1sin[j];
    kp2plusp3 += cosarray[j] * coeffkp2plusp3cos[j] + sinarray[j] * coeffkp2plusp3sin[j];
  }

  return I*PI*f*L_SI/C_SI * (n1Pn1plus*Yfactorplus + n1Pn1cross*Yfactorcross) * sinc( PI*f*L_SI/C_SI * (1.-kn1)) * cexp( I*PI*f*L_SI/C_SI * (1.+kp2plusp3) );
}
/* Function evaluating G13, combining the two polarization with the spherical harmonics factors */
double complex G13mode(const double f, const double t, const double complex Yfactorplus, const double complex Yfactorcross) {

  for(int j=0; j<4; j++) {
    cosarray[j] = cos((j+1) * Omega_SI * t);
    sinarray[j] = sin((j+1) * Omega_SI * t);
  }
  double n2Pn2plus = coeffn2Hn2plusconst;
  double n2Pn2cross = coeffn2Hn2crossconst;
  for(int j=0; j<4; j++) {
    n2Pn2plus += cosarray[j] * coeffn2Hn2pluscos[j] + sinarray[j] * coeffn2Hn2plussin[j];
    n2Pn2cross += cosarray[j] * coeffn2Hn2crosscos[j] + sinarray[j] * coeffn2Hn2crosssin[j];
  }
  double kn2 = coeffkn2const;
  double kp3plusp1 = coeffkp3plusp1const;
  for(int j=0; j<2; j++) {
    kn2 += cosarray[j] * coeffkn2cos[j] + sinarray[j] * coeffkn2sin[j];
    kp3plusp1 += cosarray[j] * coeffkp3plusp1cos[j] + sinarray[j] * coeffkp3plusp1sin[j];
  }

  return I*PI*f*L_SI/C_SI * (n2Pn2plus*Yfactorplus + n2Pn2cross*Yfactorcross) * sinc( PI*f*L_SI/C_SI * (1.+kn2)) * cexp( I*PI*f*L_SI/C_SI * (1.+kp3plusp1) );
}
/* Function evaluating G31, combining the two polarization with the spherical harmonics factors */
double complex G31mode(const double f, const double t, const double complex Yfactorplus, const double complex Yfactorcross) {

  for(int j=0; j<4; j++) {
    cosarray[j] = cos((j+1) * Omega_SI * t);
    sinarray[j] = sin((j+1) * Omega_SI * t);
  }
  double n2Pn2plus = coeffn2Hn2plusconst;
  double n2Pn2cross = coeffn2Hn2crossconst;
  for(int j=0; j<4; j++) {
    n2Pn2plus += cosarray[j] * coeffn2Hn2pluscos[j] + sinarray[j] * coeffn2Hn2plussin[j];
    n2Pn2cross += cosarray[j] * coeffn2Hn2crosscos[j] + sinarray[j] * coeffn2Hn2crosssin[j];
  }
  double kn2 = coeffkn2const;
  double kp3plusp1 = coeffkp3plusp1const;
  for(int j=0; j<2; j++) {
    kn2 += cosarray[j] * coeffkn2cos[j] + sinarray[j] * coeffkn2sin[j];
    kp3plusp1 += cosarray[j] * coeffkp3plusp1cos[j] + sinarray[j] * coeffkp3plusp1sin[j];
  }

  return I*PI*f*L_SI/C_SI * (n2Pn2plus*Yfactorplus + n2Pn2cross*Yfactorcross) * sinc( PI*f*L_SI/C_SI * (1.-kn2)) * cexp( I*PI*f*L_SI/C_SI * (1.+kp3plusp1) );
}

/* Function evaluating all coefficients G12, G21, G23, G32, G31, G13, combining the two polarization with the spherical harmonics factors */
/* Note: includes orbital delay */
int EvaluateGABmode(
  double complex* G12,                     /* Output for G12 */
  double complex* G21,                     /* Output for G21 */
  double complex* G23,                     /* Output for G23 */
  double complex* G32,                     /* Output for G32 */
  double complex* G31,                     /* Output for G31 */
  double complex* G13,                     /* Output for G13 */
  const double f,                          /* Frequency */
  const double t,                          /* Time */
  const double complex Yfactorplus,        /* Spin-weighted spherical harmonic factor for plus */
  const double complex Yfactorcross,       /* Spin-weighted spherical harmonic factor for cross */
  const int tagdelayR,                     /* Tag: when 1, include the phase term of the R-delay */
  const ResponseApproxtag responseapprox)  /* Tag to select possible low-f approximation level in FD response */
{
  /* Precompute array of sine/cosine */
  for(int j=0; j<4; j++) {
    cosarray[j] = cos((j+1) * Omega_SI * t);
    sinarray[j] = sin((j+1) * Omega_SI * t);
  }
  /* Scalar products with k */
  double n1Pn1plus = coeffn1Hn1plusconst;
  double n1Pn1cross = coeffn1Hn1crossconst;
  double n2Pn2plus = coeffn2Hn2plusconst;
  double n2Pn2cross = coeffn2Hn2crossconst;
  double n3Pn3plus = coeffn3Hn3plusconst;
  double n3Pn3cross = coeffn3Hn3crossconst;
  for(int j=0; j<4; j++) {
    n1Pn1plus += cosarray[j] * coeffn1Hn1pluscos[j] + sinarray[j] * coeffn1Hn1plussin[j];
    n1Pn1cross += cosarray[j] * coeffn1Hn1crosscos[j] + sinarray[j] * coeffn1Hn1crosssin[j];
    n2Pn2plus += cosarray[j] * coeffn2Hn2pluscos[j] + sinarray[j] * coeffn2Hn2plussin[j];
    n2Pn2cross += cosarray[j] * coeffn2Hn2crosscos[j] + sinarray[j] * coeffn2Hn2crosssin[j];
    n3Pn3plus += cosarray[j] * coeffn3Hn3pluscos[j] + sinarray[j] * coeffn3Hn3plussin[j];
    n3Pn3cross += cosarray[j] * coeffn3Hn3crosscos[j] + sinarray[j] * coeffn3Hn3crosssin[j];
  }
  /* Scalar products with k */
  double kn1 = coeffkn1const;
  double kn2 = coeffkn2const;
  double kn3 = coeffkn3const;
  double kp1plusp2 = coeffkp1plusp2const;
  double kp2plusp3 = coeffkp2plusp3const;
  double kp3plusp1 = coeffkp3plusp1const;
  double kR = coeffkRconst;
  for(int j=0; j<2; j++) {
    kn1 += cosarray[j] * coeffkn1cos[j] + sinarray[j] * coeffkn1sin[j];
    kn2 += cosarray[j] * coeffkn2cos[j] + sinarray[j] * coeffkn2sin[j];
    kn3 += cosarray[j] * coeffkn3cos[j] + sinarray[j] * coeffkn3sin[j];
    kp1plusp2 += cosarray[j] * coeffkp1plusp2cos[j] + sinarray[j] * coeffkp1plusp2sin[j];
    kp2plusp3 += cosarray[j] * coeffkp2plusp3cos[j] + sinarray[j] * coeffkp2plusp3sin[j];
    kp3plusp1 += cosarray[j] * coeffkp3plusp1cos[j] + sinarray[j] * coeffkp3plusp1sin[j];
    kR += cosarray[j] * coeffkRcos[j] + sinarray[j] * coeffkRsin[j];
  }
  /* Common factors */
  double complex factn1Pn1 = n1Pn1plus*Yfactorplus + n1Pn1cross*Yfactorcross;
  double complex factn2Pn2 = n2Pn2plus*Yfactorplus + n2Pn2cross*Yfactorcross;
  double complex factn3Pn3 = n3Pn3plus*Yfactorplus + n3Pn3cross*Yfactorcross;
  double prefactor = PI*f*L_SI/C_SI;
  double prefactorR = 2*PI*f*R_SI/C_SI;
  double complex factorcexp12 = cexp(I*prefactor * (1.+kp1plusp2));
  double complex factorcexp23 = cexp(I*prefactor * (1.+kp2plusp3));
  double complex factorcexp31 = cexp(I*prefactor * (1.+kp3plusp1));
  double factorsinc12 = sinc( prefactor * (1.-kn3));
  double factorsinc21 = sinc( prefactor * (1.+kn3));
  double factorsinc23 = sinc( prefactor * (1.-kn1));
  double factorsinc32 = sinc( prefactor * (1.+kn1));
  double factorsinc31 = sinc( prefactor * (1.-kn2));
  double factorsinc13 = sinc( prefactor * (1.+kn2));
  /* The tag tagdelayR allows to choose to include or not the R-delay phase term (here leading order) */
  double complex factorcexpkR;
  if(tagdelayR) factorcexpkR = cexp(I*prefactorR * kR);
  else factorcexpkR = 1.;

  /* Take into account level of approximation in for low-f response - choices are full, lowfL or lowf */
  if(responseapprox==lowf) {
    factorcexpkR = 1.;
  }
  if((responseapprox==lowfL)||(responseapprox==lowf)) {
    factorsinc12 = 1.;
    factorsinc21 = 1.;
    factorsinc23 = 1.;
    factorsinc32 = 1.;
    factorsinc31 = 1.;
    factorsinc13 = 1.;
    factorcexp12 = 1.;
    factorcexp23 = 1.;
    factorcexp31 = 1.;
  }

  /* Output result */
  *G12 = I*prefactor * factorcexpkR * factn3Pn3 * factorsinc12 * factorcexp12;
  *G21 = I*prefactor * factorcexpkR * factn3Pn3 * factorsinc21 * factorcexp12;
  *G23 = I*prefactor * factorcexpkR * factn1Pn1 * factorsinc23 * factorcexp23;
  *G32 = I*prefactor * factorcexpkR * factn1Pn1 * factorsinc32 * factorcexp23;
  *G31 = I*prefactor * factorcexpkR * factn2Pn2 * factorsinc31 * factorcexp31;
  *G13 = I*prefactor * factorcexpkR * factn2Pn2 * factorsinc13 * factorcexp31;

  return SUCCESS;
}

/*********************** Fourier-domain TDI factors ************************/

/* Functions evaluating the Fourier-domain factors (combinations of the GAB's) for TDI observables */
/* NOTE: factors have been scaled out, in parallel of what is done for the noise function */
/* Note: in case only one channel is considered, amplitudes for channels 2 and 3 are simply set to 0 */
/* (allows minimal changes from the old structure that assumed KTV A,E,T - but probably not optimal) */
int EvaluateTDIfactor3Chan(
  double complex* factor1,                       /* Output for factor for TDI channel 1 */
  double complex* factor2,                       /* Output for factor for TDI channel 2 */
  double complex* factor3,                       /* Output for factor for TDI channel 3 */
  const double complex G12,                      /* Input for G12 */
  const double complex G21,                      /* Input for G21 */
  const double complex G23,                      /* Input for G23 */
  const double complex G32,                      /* Input for G32 */
  const double complex G31,                      /* Input for G31 */
  const double complex G13,                      /* Input for G13 */
  const double f,                                /* Frequency */
  const TDItag tditag)                           /* Selector for the TDI observables */
{
  /* Notation: x=pifL, z=e^2ix*/
  double x = PI*f*L_SI/C_SI;
  double complex z = cexp(2*I*x);
  double sin2x = sin(2*x);
  switch(tditag) {
    /* For testing purposes: basic yAB observable - no factor */
  case y12:
    *factor1 = G12;
    *factor2 = 0.;
    *factor3 = 0.;
    break;
    /* For testing purposes: basic yABL observable - no factor, same as for yAB */
  case y12L:
    *factor1 = G12;
    *factor2 = 0.;
    *factor3 = 0.;
    break;
    /* First-generation rescaled TDI aet from X,Y,Z */
    /* With x=pifL, factors scaled out: A,E I*sqrt2*sin2x*e2ix - T 2*sqrt2*sin2x*sinx*e3ix */
  case TDIAETXYZ:
    *factor1 = 0.5 * ( (1.+z)*(G31+G13) - G23 - z*G32 - G21 - z*G12 );
    *factor2 = 0.5*invsqrt3 * ( (1.-z)*(G13-G31) + (2.+z)*(G12-G32) + (1.+2*z)*(G21-G23) );
    *factor3 = invsqrt6 * ( G21-G12 + G32-G23 + G13-G31);
    break;
    /* First-generation rescaled TDI aet from alpha, beta, gamma */
    /* With x=pifL, factors scaled out: A,E -I*2sqrt2*sinx*eix - T sinx/(sin3x*eix) */
  case TDIAETalphabetagamma:
    *factor1 = 0.5 * (G13+G31 + z*(G12+G32) - (1.+z)*(G21+G13));
    *factor2 = 0.5*invsqrt3 * ((2.+z)*(G12-G32) + (1.+z)*(G21-G23) + (1.+2*z)*(G13-G31));
    *factor3 = invsqrt3 * (G21-G12 + G32-G23 + G13-G31);
    break;
    /* First-generation TDI XYZ */
    /* With x=pifL, factor scaled out: 2I*sin2x*e2ix */
  case TDIXYZ:
    *factor1 = G21 + z*G12 - G31 - z*G13;
    *factor2 = G32 + z*G23 - G12 - z*G21;
    *factor3 = G13 + z*G31 - G23 - z*G32;
    break;
    /* First-generation TDI alpha beta gamma */
  case TDIalphabetagamma:
    *factor1 = G21-G31 + z*(G13-G12) + z*z*(G32-G23);
    *factor2 = G32-G12 + z*(G21-G23) + z*z*(G13-G31);
    *factor3 = G13-G23 + z*(G32-G31) + z*z*(G21-G12);
    break;
  /* First-generation TDI XYZ */
  case TDIX:
    *factor1 = G21 + z*G12 - G31 - z*G13;
    *factor2 = 0.;
    *factor3 = 0.;
    break;
  /* First-generation TDI alpha beta gamma */
  case TDIalpha:
    *factor1 = G21-G31 + z*(G13-G12) + z*z*(G32-G23);
    *factor2 = 0.;
    *factor3 = 0.;
    break;
  /* First-generation rescaled TDI aet from X,Y,Z */
  /* With x=pifL, factors scaled out: A,E I*sqrt2*sin2x*eix - T 2*sqrt2*sin2x*sinx*e2ix */
  case TDIAXYZ:
    *factor1 = 0.5 * ( (1.+z)*(G31+G13) - G23 - z*G32 - G21 - z*G12 );
    *factor2 = 0.;
    *factor3 = 0.;
    break;
  case TDIEXYZ:
    *factor1 = 0.5*invsqrt3 * ( (1.-z)*(G13-G31) + (2.+z)*(G12-G32) + (1.+2*z)*(G12-G23) );
    *factor2 = 0.;
    *factor3 = 0.;
    break;
  case TDITXYZ:
    *factor1 = invsqrt6 * ( G21-G12 + G32-G23 + G13-G31);
    *factor2 = 0.;
    *factor3 = 0.;
    break;
  /* First-generation rescaled TDI aet from alpha, beta, gamma */
  /* With x=pifL, factors scaled out: A,E -I*2sqrt2*sinx*eix - T sinx/(sin3x*eix) */
  case TDIAalphabetagamma:
    *factor1 = 0.5 * (G13+G31 + z*(G12+G32) - (1.+z)*(G21+G13));
    *factor2 = 0.;
    *factor3 = 0.;
    break;
  case TDIEalphabetagamma:
    *factor1 = 0.5*invsqrt3 * ((2.+z)*(G12-G32) + (1.+z)*(G21-G23) + (1.+2*z)*(G13-G31));
    *factor2 = 0.;
    *factor3 = 0.;
    break;
  case TDITalphabetagamma:
    *factor1 = invsqrt3 * (G21-G12 + G32-G23 + G13-G31);
    *factor2 = 0.;
    *factor3 = 0.;
    break;
  default:
    printf("Error in EvaluateTDIfactor3Chan: tditag not recognized.\n");
    exit(1);
  }
  return SUCCESS;
}

/* Function evaluating the Fourier-domain factors that have been scaled out of TDI observables */
/* The factors scaled out, parallel what is done for the noise functions */
/* Note: in case only one channel is considered, factors for channels 2 and 3 are simply set to 0 */
int ScaledTDIfactor3Chan(
  double complex* factor1,                       /* Output for factor for TDI factor 1 */
  double complex* factor2,                       /* Output for factor for TDI factor 2 */
  double complex* factor3,                       /* Output for factor for TDI factor 3 */
  const double f,                                /* Frequency */
  const TDItag tditag)                           /* Selector for the TDI observables */
{
  /* Notation: x=pifL */
  double x = PI*f*L_SI/C_SI;
  switch(tditag) {
    /* First-generation rescaled TDI aet from X,Y,Z */
  case TDIAETXYZ:
    *factor1 = I*sqrt(2)*sin(2*x)*cexp(2*I*x);
    *factor2 = I*sqrt(2)*sin(2*x)*cexp(2*I*x);
    *factor3 = 2*sqrt(2)*sin(x)*sin(2*x)*cexp(3*I*x);
    break;
    /* First-generation rescaled TDI aet from alpha, beta, gamma */
  case TDIAETalphabetagamma:
    *factor1 = -I*2*sqrt(2)*sin(x)*cexp(I*x);
    *factor2 = -I*2*sqrt(2)*sin(x)*cexp(I*x);
    *factor3 = sin(3*x)/sin(x)*cexp(I*x);
    break;
    /* First-generation TDI XYZ */
  case TDIXYZ:
    *factor1 = 2*I*sin(2*x)*cexp(2*I*x);
    *factor2 = 2*I*sin(2*x)*cexp(2*I*x);
    *factor3 = 2*I*sin(2*x)*cexp(2*I*x);
    break;
    /* First-generation TDI alpha beta gamma */
  case TDIalphabetagamma:
    *factor1 = 1.;
    *factor2 = 1.;
    *factor3 = 1.;
    break;
  /* First-generation TDI XYZ */
  case TDIX:
    *factor1 = 2*I*sin(2*x)*cexp(2*I*x);
    *factor2 = 0.;
    *factor3 = 0.;
    break;
  /* First-generation TDI alpha beta gamma */
  case TDIalpha:
    *factor1 = 1.;
    *factor2 = 0.;
    *factor3 = 0.;
    break;
  /* First-generation rescaled TDI aet from X,Y,Z */
  /* With x=pifL, factors scaled out: A,E I*sqrt2*sin2x*eix - T 2*sqrt2*sin2x*sinx*e2ix */
  case TDIAXYZ:
    *factor1 = I*sqrt(2)*sin(2*x)*cexp(2*I*x);
    *factor2 = 0.;
    *factor3 = 0.;
    break;
  case TDIEXYZ:
    *factor1 = I*sqrt(2)*sin(2*x)*cexp(2*I*x);
    *factor2 = 0.;
    *factor3 = 0.;
    break;
  case TDITXYZ:
    *factor1 = 2*sqrt(2)*sin(x)*sin(2*x)*cexp(3*I*x);
    *factor2 = 0.;
    *factor3 = 0.;
    break;
  /* First-generation rescaled TDI aet from alpha, beta, gamma */
  /* With x=pifL, factors scaled out: A,E -I*2sqrt2*sinx*eix - T sinx/(sin3x*eix) */
  case TDIAalphabetagamma:
    *factor1 = -I*2*sqrt(2)*sin(x)*cexp(I*x);
    *factor2 = 0.;
    *factor3 = 0.;
    break;
  case TDIEalphabetagamma:
    *factor1 = -I*2*sqrt(2)*sin(x)*cexp(I*x);
    *factor2 = 0.;
    *factor3 = 0.;
    break;
  case TDITalphabetagamma:
    *factor1 = sin(3*x)/sin(x)*cexp(I*x);
    *factor2 = 0.;
    *factor3 = 0.;
    break;
  default:
    printf("Error in EvaluateTDIfactor3Chan: tditag not recognized.\n");
    exit(1);
  }
  return SUCCESS;
}

/* Function restoring the factor that have been scaled out of the TDI observables */
/* NOTE: the operation is made in-place, and the input is overwritten */
int RestoreInPlaceScaledFactorTDI(
  ListmodesCAmpPhaseFrequencySeries* listtdi,     /* Output/Input: list of mode contributions to TDI observable */
  TDItag tditag,                                  /* Tag selecting the TDI observable */
  int nchannel)                                   /* TDI channel number */
{
  double complex factor1 = 0;
  double complex factor2 = 0;
  double complex factor3 = 0;
  double complex factor;
  double complex camp;
  ListmodesCAmpPhaseFrequencySeries* listelement = listtdi;
  /* Going throug the list of modes */
  while(listelement) {
    gsl_vector* freq = listelement->freqseries->freq;
    gsl_vector* ampreal = listelement->freqseries->amp_real;
    gsl_vector* ampimag = listelement->freqseries->amp_imag;
    for(int i=0; i<freq->size; i++) {
      ScaledTDIfactor3Chan(&factor1, &factor2, &factor3, gsl_vector_get(freq, i), tditag);
      switch(nchannel) {
      case 1: factor = factor1; break;
      case 2: factor = factor2; break;
      case 3: factor = factor3; break;
      }
      camp = factor * (gsl_vector_get(ampreal, i) + I*gsl_vector_get(ampimag, i));
      gsl_vector_set(ampreal, i, creal(camp));
      gsl_vector_set(ampimag, i, cimag(camp));
    }
    listelement = listelement->next;
  }
  return SUCCESS;
}


/* Functions evaluating the Fourier-domain factors (combinations of the GAB's) for TDI observables */
/* int EvaluateTDIfactor1Chan( */
/*   double complex* factor,                        /\* Output for factor for TDI channel *\/ */
/*   const double complex G12,                      /\* Input for G12 *\/ */
/*   const double complex G21,                      /\* Input for G21 *\/ */
/*   const double complex G23,                      /\* Input for G23 *\/ */
/*   const double complex G32,                      /\* Input for G32 *\/ */
/*   const double complex G31,                      /\* Input for G31 *\/ */
/*   const double complex G13,                      /\* Input for G13 *\/ */
/*   const double f,                                /\* Frequency *\/ */
/*   const TDItag tditag)                           /\* Selector for the TDI observables *\/ */
/* { */
/*   /\* Notation: x=pifL, z = e^2ix*\/ */
/*   double x = PI*f*L_SI/C_SI; */
/*   double complex z = cexp(2*I*x); */
/*   double sin2x = sin(2*x); */
/*   double complex commonfac; */
/*   switch(tditag) { */
/*   /\* First-generation TDI XYZ *\/ */
/*   case TDIX: { */
/*     commonfac = 2*I*z*sin2x; */
/*     *factor = commonfac * (G21 + z*G12 - G31 - z*G13); } */
/*   case TDIY: { */
/*     commonfac = 2*I*z*sin2x; */
/*     *factor = commonfac * (G32 + z*G23 - G12 - z*G21); } */
/*   case TDIZ: { */
/*     commonfac = 2*I*z*sin2x; */
/*     *factor = commonfac * (G13 + z*G31 - G23 - z*G32); } */
/*   /\* First-generation TDI alpha beta gamma *\/ */
/*   case TDIalpha: { */
/*     *factor = G21-G31 + z*(G13-G12) + z*z*(G32-G23); } */
/*   case TDIbeta: { */
/*     *factor = G32-G12 + z*(G21-G23) + z*z*(G13-G31); } */
/*   case TDIgamma: { */
/*     *factor = G13-G23 + z*(G32-G31) + z*z*(G21-G12); } */
/*   /\* First-generation rescaled TDI aet from X,Y,Z *\/ */
/*   /\* With x=pifL, factors scaled out: A,E I*sqrt2*sin2x*eix - T 2*sqrt2*sin2x*sinx*e2ix *\/ */
/*   case TDIAXYZ: { */
/*     *factor = 0.5 * ( (1.+z)*(G31+G13) - G23 - z*G32 - G21 - z*G12 ); } */
/*   case TDIEXYZ: { */
/*     *factor = 0.5*invsqrt3 * ( (1.-z)*(G13-G31) + (2.+z)*(G12-G32) + (1.+2*z)*(G12-G23) ); } */
/*   case TDITXYZ: { */
/*     *factor = invsqrt6 * ( G21-G12 + G32-G23 + G13-G31); } */
/*   /\* First-generation rescaled TDI aet from alpha, beta, gamma *\/ */
/*   /\* With x=pifL, factors scaled out: A,E -I*2sqrt2*sinx*eix - T sinx/(sin3x*eix) *\/ */
/*   case TDIAalphabetagamma: { */
/*     *factor = 0.5 * (G13+G31 + z*(G12+G32) - (1.+z)*(G21+G13)); } */
/*   case TDIEalphabetagamma: { */
/*     *factor = 0.5*invsqrt3 * ((2+z)*(G12-G32) + (1+z)*(G21-G23) + (1.+2*z)*(G13-G31)); } */
/*   case TDITalphabetagamma: { */
/*     *factor = invsqrt3 * (G21-G12 + G32-G23 + G13-G31); } */
/*   default: { */
/*     printf("Error in EvaluateTDIfactor3Chan: tditag not recognized."); */
/*     exit(1); } */
/*   } */
/* } */

/*********************** Time-domain response ************************/

/* Processing single mode in amp/phase form through orbital time delay */
double hOTDAmpPhase(
  double* amp,                             /* Output: amplitude */
  double* phase,                           /* Output: phase */
  gsl_spline* splineamp,                   /* Input spline for TD mode amplitude */
  gsl_spline* splinephase,                 /* Input spline for TD mode phase */
  gsl_interp_accel* accelamp,              /* Accelerator for amp spline */
  gsl_interp_accel* accelphase,            /* Accelerator for phase spline */
  const double t)                          /* Time */
{
  /* Precompute array of sine/cosine */
  for(int j=0; j<4; j++) {
    cosarray[j] = cos((j+1) * Omega_SI * t);
    sinarray[j] = sin((j+1) * Omega_SI * t);
  }
  /* Scalar product k.R */
  double kR = coeffkRconst;
  for(int j=0; j<2; j++) {
    kR += cosarray[j] * coeffkRcos[j] + sinarray[j] * coeffkRsin[j];
  }
  /* Common factor and delay */
  double delay = -(kR*R_SI)/C_SI;

  /* Output result */
  *amp = gsl_spline_eval(splineamp, t+delay, accelamp);
  *phase = gsl_spline_eval(splinephase, t+delay, accelphase);
}

/* Functions evaluating yAB observables in time domain */
/* Note: includes both h22 and h2m2 contributions, assuming planar orbits so that h2-2 = h22* */
double y12LTDfromh22AmpPhase(
  gsl_spline* splineamp,                   /* Input spline for h22 TD amp */
  gsl_spline* splinephase,                 /* Input spline for h22 TD phase */
  gsl_interp_accel* accelamp,              /* Accelerator for amp spline */
  gsl_interp_accel* accelphase,            /* Accelerator for phase spline */
  double complex Y22,                      /* Y22 factor needed to convert h22 to hplus, hcross */
  double complex Y2m2,                     /* Y2-2 factor needed to convert h2-2 to hplus, hcross */
  const double t)                          /* Time */
{
  /* Precompute array of sine/cosine */
  for(int j=0; j<4; j++) {
    cosarray[j] = cos((j+1) * Omega_SI * t);
    sinarray[j] = sin((j+1) * Omega_SI * t);
  }
  /* Scalar products with k */
  double n3Pn3plus = coeffn3Hn3plusconst;
  double n3Pn3cross = coeffn3Hn3crossconst;
  for(int j=0; j<4; j++) {
    n3Pn3plus += cosarray[j] * coeffn3Hn3pluscos[j] + sinarray[j] * coeffn3Hn3plussin[j];
    n3Pn3cross += cosarray[j] * coeffn3Hn3crosscos[j] + sinarray[j] * coeffn3Hn3crosssin[j];
  }
  /* Scalar products with k */
  double kn3 = coeffkn3const;
  double kp1 = coeffkp1const;
  double kp2 = coeffkp2const;
  for(int j=0; j<2; j++) {
    kn3 += cosarray[j] * coeffkn3cos[j] + sinarray[j] * coeffkn3sin[j];
    kp1 += cosarray[j] * coeffkp1cos[j] + sinarray[j] * coeffkp1sin[j];
    kp2 += cosarray[j] * coeffkp2cos[j] + sinarray[j] * coeffkp2sin[j];
  }
  /* Common factor and delay */
  double factorp = (1./(1.-kn3)) * 0.5*n3Pn3plus;
  double factorc = (1./(1.-kn3)) * 0.5*n3Pn3cross;
  double firstdelay = -((kp1 + 1)*L_SI)/C_SI;
  double seconddelay = -(kp2*L_SI)/C_SI;

  /* Values of Y22*h22 + Y2-2*h2-2 at 1 and 2 with delays, and hplus, hcross */
  /* Note: includes both h22 and h2m2 contributions, assuming planar orbits so that h2-2 = h22* */
  double A22at1 = gsl_spline_eval(splineamp, t+firstdelay, accelamp);
  double phi22at1 = gsl_spline_eval(splinephase, t+firstdelay, accelphase);
  double A22at2 = gsl_spline_eval(splineamp, t+seconddelay, accelamp);
  double phi22at2 = gsl_spline_eval(splinephase, t+seconddelay, accelphase);
  double complex Y22h22at1 = Y22 * A22at1 * cexp(I*phi22at1);
  double complex Y22h22at2 = Y22 * A22at2 * cexp(I*phi22at2);
  double complex Y2m2h2m2at1 = Y2m2 * A22at1 * cexp(-I*phi22at1);
  double complex Y2m2h2m2at2 = Y2m2 * A22at2 * cexp(-I*phi22at2);
  double hp1 = creal(Y22h22at1 + Y2m2h2m2at1);
  double hc1 = -cimag(Y22h22at1 + Y2m2h2m2at1);
  double hp2 = creal(Y22h22at2 + Y2m2h2m2at2);
  double hc2 = -cimag(Y22h22at2 + Y2m2h2m2at2);

  /* Result */
  double y12 = factorp*(hp1 - hp2) + factorc*(hc1 - hc2);
  return y12;
}

/* Functions evaluating yAB observables in time domain */
double y12TD(
  gsl_spline* splinehp,                    /* Input spline for TD hplus */
  gsl_spline* splinehc,                    /* Input spline for TD hcross */
  gsl_interp_accel* accelhp,               /* Accelerator for hp spline */
  gsl_interp_accel* accelhc,               /* Accelerator for hc spline */
  const double t)                          /* Time */
{
  /* Precompute array of sine/cosine */
  for(int j=0; j<4; j++) {
    cosarray[j] = cos((j+1) * Omega_SI * t);
    sinarray[j] = sin((j+1) * Omega_SI * t);
  }
  /* Scalar products with k */
  double n3Pn3plus = coeffn3Hn3plusconst;
  double n3Pn3cross = coeffn3Hn3crossconst;
  for(int j=0; j<4; j++) {
    n3Pn3plus += cosarray[j] * coeffn3Hn3pluscos[j] + sinarray[j] * coeffn3Hn3plussin[j];
    n3Pn3cross += cosarray[j] * coeffn3Hn3crosscos[j] + sinarray[j] * coeffn3Hn3crosssin[j];
  }
  /* Scalar products with k */
  double kn3 = coeffkn3const;
  double kp1 = coeffkp1const;
  double kp2 = coeffkp2const;
  double kR = coeffkRconst;
  for(int j=0; j<2; j++) {
    kn3 += cosarray[j] * coeffkn3cos[j] + sinarray[j] * coeffkn3sin[j];
    kp1 += cosarray[j] * coeffkp1cos[j] + sinarray[j] * coeffkp1sin[j];
    kp2 += cosarray[j] * coeffkp2cos[j] + sinarray[j] * coeffkp2sin[j];
    kR += cosarray[j] * coeffkRcos[j] + sinarray[j] * coeffkRsin[j];
  }
  /* Common factor and delay */
  double factorp = (1./(1.-kn3)) * 0.5*n3Pn3plus;
  double factorc = (1./(1.-kn3)) * 0.5*n3Pn3cross;
  double firstdelay = -(kR*R_SI + (kp1 + 1)*L_SI)/C_SI;
  double seconddelay = -(kR*R_SI + kp2*L_SI)/C_SI;

  /* Result */
  double y12 = factorp*(gsl_spline_eval(splinehp, t+firstdelay, accelhp) - gsl_spline_eval(splinehp, t+seconddelay, accelhp)) + factorc*(gsl_spline_eval(splinehc, t+firstdelay, accelhc) - gsl_spline_eval(splinehc, t+seconddelay, accelhc));
  return y12;
}
double y21TD(
  gsl_spline* splinehp,                    /* Input spline for TD hplus */
  gsl_spline* splinehc,                    /* Input spline for TD hcross */
  gsl_interp_accel* accelhp,               /* Accelerator for hp spline */
  gsl_interp_accel* accelhc,               /* Accelerator for hc spline */
  const double t)                          /* Time */
{
  /* Precompute array of sine/cosine */
  for(int j=0; j<4; j++) {
    cosarray[j] = cos((j+1) * Omega_SI * t);
    sinarray[j] = sin((j+1) * Omega_SI * t);
  }
  /* Scalar products with k */
  double n3Pn3plus = coeffn3Hn3plusconst;
  double n3Pn3cross = coeffn3Hn3crossconst;
  for(int j=0; j<4; j++) {
    n3Pn3plus += cosarray[j] * coeffn3Hn3pluscos[j] + sinarray[j] * coeffn3Hn3plussin[j];
    n3Pn3cross += cosarray[j] * coeffn3Hn3crosscos[j] + sinarray[j] * coeffn3Hn3crosssin[j];
  }
  /* Scalar products with k */
  double kn3 = coeffkn3const;
  double kp1 = coeffkp1const;
  double kp2 = coeffkp2const;
  double kR = coeffkRconst;
  for(int j=0; j<2; j++) {
    kn3 += cosarray[j] * coeffkn3cos[j] + sinarray[j] * coeffkn3sin[j];
    kp1 += cosarray[j] * coeffkp1cos[j] + sinarray[j] * coeffkp1sin[j];
    kp2 += cosarray[j] * coeffkp2cos[j] + sinarray[j] * coeffkp2sin[j];
    kR += cosarray[j] * coeffkRcos[j] + sinarray[j] * coeffkRsin[j];
  }
  /* Common factor and delay */
  double factorp = (1./(1.+kn3)) * 0.5*n3Pn3plus;
  double factorc = (1./(1.+kn3)) * 0.5*n3Pn3cross;
  double firstdelay = -(kR*R_SI + (kp2 + 1)*L_SI)/C_SI;
  double seconddelay = -(kR*R_SI + kp1*L_SI)/C_SI;

  /* Result */
  double y21 = factorp*(gsl_spline_eval(splinehp, t+firstdelay, accelhp) - gsl_spline_eval(splinehp, t+seconddelay, accelhp)) + factorc*(gsl_spline_eval(splinehc, t+firstdelay, accelhc) - gsl_spline_eval(splinehc, t+seconddelay, accelhc));
  return y21;
}
double y23TD(
  gsl_spline* splinehp,                    /* Input spline for TD hplus */
  gsl_spline* splinehc,                    /* Input spline for TD hcross */
  gsl_interp_accel* accelhp,               /* Accelerator for hp spline */
  gsl_interp_accel* accelhc,               /* Accelerator for hc spline */
  const double t)                          /* Time */
{
  /* Precompute array of sine/cosine */
  for(int j=0; j<4; j++) {
    cosarray[j] = cos((j+1) * Omega_SI * t);
    sinarray[j] = sin((j+1) * Omega_SI * t);
  }
  /* Scalar products with k */
  double n1Pn1plus = coeffn1Hn1plusconst;
  double n1Pn1cross = coeffn1Hn1crossconst;
  for(int j=0; j<4; j++) {
    n1Pn1plus += cosarray[j] * coeffn1Hn1pluscos[j] + sinarray[j] * coeffn1Hn1plussin[j];
    n1Pn1cross += cosarray[j] * coeffn1Hn1crosscos[j] + sinarray[j] * coeffn1Hn1crosssin[j];
  }
  /* Scalar products with k */
  double kn1 = coeffkn1const;
  double kp2 = coeffkp2const;
  double kp3 = coeffkp3const;
  double kR = coeffkRconst;
  for(int j=0; j<2; j++) {
    kn1 += cosarray[j] * coeffkn1cos[j] + sinarray[j] * coeffkn1sin[j];
    kp2 += cosarray[j] * coeffkp2cos[j] + sinarray[j] * coeffkp2sin[j];
    kp3 += cosarray[j] * coeffkp3cos[j] + sinarray[j] * coeffkp3sin[j];
    kR += cosarray[j] * coeffkRcos[j] + sinarray[j] * coeffkRsin[j];
  }
  /* Common factor and delay */
  double factorp = (1./(1.-kn1)) * 0.5*n1Pn1plus;
  double factorc = (1./(1.-kn1)) * 0.5*n1Pn1cross;
  double firstdelay = -(kR*R_SI + (kp2 + 1)*L_SI)/C_SI;
  double seconddelay = -(kR*R_SI + kp3*L_SI)/C_SI;

  /* Result */
  double y23 = factorp*(gsl_spline_eval(splinehp, t+firstdelay, accelhp) - gsl_spline_eval(splinehp, t+seconddelay, accelhp)) + factorc*(gsl_spline_eval(splinehc, t+firstdelay, accelhc) - gsl_spline_eval(splinehc, t+seconddelay, accelhc));
  return y23;
}
double y32TD(
  gsl_spline* splinehp,                    /* Input spline for TD hplus */
  gsl_spline* splinehc,                    /* Input spline for TD hcross */
  gsl_interp_accel* accelhp,               /* Accelerator for hp spline */
  gsl_interp_accel* accelhc,               /* Accelerator for hc spline */
  const double t)                          /* Time */
{
  /* Precompute array of sine/cosine */
  for(int j=0; j<4; j++) {
    cosarray[j] = cos((j+1) * Omega_SI * t);
    sinarray[j] = sin((j+1) * Omega_SI * t);
  }
  /* Scalar products with k */
  double n1Pn1plus = coeffn1Hn1plusconst;
  double n1Pn1cross = coeffn1Hn1crossconst;
  for(int j=0; j<4; j++) {
    n1Pn1plus += cosarray[j] * coeffn1Hn1pluscos[j] + sinarray[j] * coeffn1Hn1plussin[j];
    n1Pn1cross += cosarray[j] * coeffn1Hn1crosscos[j] + sinarray[j] * coeffn1Hn1crosssin[j];
  }
  /* Scalar products with k */
  double kn1 = coeffkn1const;
  double kp2 = coeffkp2const;
  double kp3 = coeffkp3const;
  double kR = coeffkRconst;
  for(int j=0; j<2; j++) {
    kn1 += cosarray[j] * coeffkn1cos[j] + sinarray[j] * coeffkn1sin[j];
    kp2 += cosarray[j] * coeffkp2cos[j] + sinarray[j] * coeffkp2sin[j];
    kp3 += cosarray[j] * coeffkp3cos[j] + sinarray[j] * coeffkp3sin[j];
    kR += cosarray[j] * coeffkRcos[j] + sinarray[j] * coeffkRsin[j];
  }
  /* Common factor and delay */
  double factorp = (1./(1.+kn1)) * 0.5*n1Pn1plus;
  double factorc = (1./(1.+kn1)) * 0.5*n1Pn1cross;
  double firstdelay = -(kR*R_SI + (kp3 + 1)*L_SI)/C_SI;
  double seconddelay = -(kR*R_SI + kp2*L_SI)/C_SI;

  /* Result */
  double y32 = factorp*(gsl_spline_eval(splinehp, t+firstdelay, accelhp) - gsl_spline_eval(splinehp, t+seconddelay, accelhp)) + factorc*(gsl_spline_eval(splinehc, t+firstdelay, accelhc) - gsl_spline_eval(splinehc, t+seconddelay, accelhc));
  return y32;
}
double y31TD(
  gsl_spline* splinehp,                    /* Input spline for TD hplus */
  gsl_spline* splinehc,                    /* Input spline for TD hcross */
  gsl_interp_accel* accelhp,               /* Accelerator for hp spline */
  gsl_interp_accel* accelhc,               /* Accelerator for hc spline */
  const double t)                          /* Time */
{
  /* Precompute array of sine/cosine */
  for(int j=0; j<4; j++) {
    cosarray[j] = cos((j+1) * Omega_SI * t);
    sinarray[j] = sin((j+1) * Omega_SI * t);
  }
  /* Scalar products with k */
  double n2Pn2plus = coeffn2Hn2plusconst;
  double n2Pn2cross = coeffn2Hn2crossconst;
  for(int j=0; j<4; j++) {
    n2Pn2plus += cosarray[j] * coeffn2Hn2pluscos[j] + sinarray[j] * coeffn2Hn2plussin[j];
    n2Pn2cross += cosarray[j] * coeffn2Hn2crosscos[j] + sinarray[j] * coeffn2Hn2crosssin[j];
  }
  /* Scalar products with k */
  double kn2 = coeffkn2const;
  double kp3 = coeffkp3const;
  double kp1 = coeffkp1const;
  double kR = coeffkRconst;
  for(int j=0; j<2; j++) {
    kn2 += cosarray[j] * coeffkn2cos[j] + sinarray[j] * coeffkn2sin[j];
    kp3 += cosarray[j] * coeffkp3cos[j] + sinarray[j] * coeffkp3sin[j];
    kp1 += cosarray[j] * coeffkp1cos[j] + sinarray[j] * coeffkp1sin[j];
    kR += cosarray[j] * coeffkRcos[j] + sinarray[j] * coeffkRsin[j];
  }
  /* Common factor and delay */
  double factorp = (1./(1.-kn2)) * 0.5*n2Pn2plus;
  double factorc = (1./(1.-kn2)) * 0.5*n2Pn2cross;
  double firstdelay = -(kR*R_SI + (kp3 + 1)*L_SI)/C_SI;
  double seconddelay = -(kR*R_SI + kp1*L_SI)/C_SI;

  /* Result */
  double y31 = factorp*(gsl_spline_eval(splinehp, t+firstdelay, accelhp) - gsl_spline_eval(splinehp, t+seconddelay, accelhp)) + factorc*(gsl_spline_eval(splinehc, t+firstdelay, accelhc) - gsl_spline_eval(splinehc, t+seconddelay, accelhc));
  return y31;
}
double y13TD(
  gsl_spline* splinehp,                    /* Input spline for TD hplus */
  gsl_spline* splinehc,                    /* Input spline for TD hcross */
  gsl_interp_accel* accelhp,               /* Accelerator for hp spline */
  gsl_interp_accel* accelhc,               /* Accelerator for hc spline */
  const double t)                          /* Time */
{
  /* Precompute array of sine/cosine */
  for(int j=0; j<4; j++) {
    cosarray[j] = cos((j+1) * Omega_SI * t);
    sinarray[j] = sin((j+1) * Omega_SI * t);
  }
  /* Scalar products with k */
  double n2Pn2plus = coeffn2Hn2plusconst;
  double n2Pn2cross = coeffn2Hn2crossconst;
  for(int j=0; j<4; j++) {
    n2Pn2plus += cosarray[j] * coeffn2Hn2pluscos[j] + sinarray[j] * coeffn2Hn2plussin[j];
    n2Pn2cross += cosarray[j] * coeffn2Hn2crosscos[j] + sinarray[j] * coeffn2Hn2crosssin[j];
  }
  /* Scalar products with k */
  double kn2 = coeffkn2const;
  double kp3 = coeffkp3const;
  double kp1 = coeffkp1const;
  double kR = coeffkRconst;
  for(int j=0; j<2; j++) {
    kn2 += cosarray[j] * coeffkn2cos[j] + sinarray[j] * coeffkn2sin[j];
    kp3 += cosarray[j] * coeffkp3cos[j] + sinarray[j] * coeffkp3sin[j];
    kp1 += cosarray[j] * coeffkp1cos[j] + sinarray[j] * coeffkp1sin[j];
    kR += cosarray[j] * coeffkRcos[j] + sinarray[j] * coeffkRsin[j];
  }
  /* Common factor and delay */
  double factorp = (1./(1.+kn2)) * 0.5*n2Pn2plus;
  double factorc = (1./(1.+kn2)) * 0.5*n2Pn2cross;
  double firstdelay = -(kR*R_SI + (kp1 + 1)*L_SI)/C_SI;
  double seconddelay = -(kR*R_SI + kp3*L_SI)/C_SI;

  /* Result */
  double y13 = factorp*(gsl_spline_eval(splinehp, t+firstdelay, accelhp) - gsl_spline_eval(splinehp, t+seconddelay, accelhp)) + factorc*(gsl_spline_eval(splinehc, t+firstdelay, accelhc) - gsl_spline_eval(splinehc, t+seconddelay, accelhc));
  return y13;
}

/**/
int EvaluateTDIXYZTD(
  double* TDIX,                            /* Output: value of TDI observable X */
  double* TDIY,                            /* Output: value of TDI observable Y */
  double* TDIZ,                            /* Output: value of TDI observable Z */
  gsl_spline* splinehp,                    /* Input spline for TD hplus */
  gsl_spline* splinehc,                    /* Input spline for TD hcross */
  gsl_interp_accel* accelhp,               /* Accelerator for hp spline */
  gsl_interp_accel* accelhc,               /* Accelerator for hc spline */
  const double t)                          /* Time */
{
  double armdelay = L_SI/C_SI;
  double X = (y31TD(splinehp, splinehc, accelhp, accelhc, t) + y13TD(splinehp, splinehc, accelhp, accelhc, t - armdelay)) + (y21TD(splinehp, splinehc, accelhp, accelhc, t - 2*armdelay) + y12TD(splinehp, splinehc, accelhp, accelhc, t - 3*armdelay)) - (y21TD(splinehp, splinehc, accelhp, accelhc, t) + y12TD(splinehp, splinehc, accelhp, accelhc, t - armdelay)) - (y31TD(splinehp, splinehc, accelhp, accelhc, t - 2*armdelay) + y13TD(splinehp, splinehc, accelhp, accelhc, t - 3*armdelay));
  double Y = (y12TD(splinehp, splinehc, accelhp, accelhc, t) + y21TD(splinehp, splinehc, accelhp, accelhc, t - armdelay)) + (y32TD(splinehp, splinehc, accelhp, accelhc, t - 2*armdelay) + y23TD(splinehp, splinehc, accelhp, accelhc, t - 3*armdelay)) - (y32TD(splinehp, splinehc, accelhp, accelhc, t) + y23TD(splinehp, splinehc, accelhp, accelhc, t - armdelay)) - (y12TD (splinehp, splinehc, accelhp, accelhc, t - 2*armdelay) + y21TD(splinehp, splinehc, accelhp, accelhc, t - 3*armdelay));
  double Z = (y23TD(splinehp, splinehc, accelhp, accelhc, t) + y32TD(splinehp, splinehc, accelhp, accelhc, t - armdelay)) + (y13TD(splinehp, splinehc, accelhp, accelhc, t - 2*armdelay) + y31TD(splinehp, splinehc, accelhp, accelhc, t - 3*armdelay)) - (y13TD(splinehp, splinehc, accelhp, accelhc, t) + y31TD(splinehp, splinehc, accelhp, accelhc, t - armdelay)) - (y23TD(splinehp, splinehc, accelhp, accelhc, t - 2*armdelay) + y32TD(splinehp, splinehc, accelhp, accelhc, t - 3*armdelay));

  /* Output */
  *TDIX = X;
  *TDIY = Y;
  *TDIZ = Z;

  return SUCCESS;
}

/**/
int EvaluateTDIAETXYZTD(
  double* TDIA,                            /* Output: value of TDI observable X */
  double* TDIE,                            /* Output: value of TDI observable Y */
  double* TDIT,                            /* Output: value of TDI observable Z */
  gsl_spline* splinehp,                    /* Input spline for TD hplus */
  gsl_spline* splinehc,                    /* Input spline for TD hcross */
  gsl_interp_accel* accelhp,               /* Accelerator for hp spline */
  gsl_interp_accel* accelhc,               /* Accelerator for hc spline */
  const double t)                          /* Time */
{
  double armdelay = L_SI/C_SI;
  double X = (y31TD(splinehp, splinehc, accelhp, accelhc, t) + y13TD(splinehp, splinehc, accelhp, accelhc, t - armdelay)) + (y21TD(splinehp, splinehc, accelhp, accelhc, t - 2*armdelay) + y12TD(splinehp, splinehc, accelhp, accelhc, t - 3*armdelay)) - (y21TD(splinehp, splinehc, accelhp, accelhc, t) + y12TD(splinehp, splinehc, accelhp, accelhc, t - armdelay)) - (y31TD(splinehp, splinehc, accelhp, accelhc, t - 2*armdelay) + y13TD(splinehp, splinehc, accelhp, accelhc, t - 3*armdelay));
  double Y = (y12TD(splinehp, splinehc, accelhp, accelhc, t) + y21TD(splinehp, splinehc, accelhp, accelhc, t - armdelay)) + (y32TD(splinehp, splinehc, accelhp, accelhc, t - 2*armdelay) + y23TD(splinehp, splinehc, accelhp, accelhc, t - 3*armdelay)) - (y32TD(splinehp, splinehc, accelhp, accelhc, t) + y23TD(splinehp, splinehc, accelhp, accelhc, t - armdelay)) - (y12TD (splinehp, splinehc, accelhp, accelhc, t - 2*armdelay) + y21TD(splinehp, splinehc, accelhp, accelhc, t - 3*armdelay));
  double Z = (y23TD(splinehp, splinehc, accelhp, accelhc, t) + y32TD(splinehp, splinehc, accelhp, accelhc, t - armdelay)) + (y13TD(splinehp, splinehc, accelhp, accelhc, t - 2*armdelay) + y31TD(splinehp, splinehc, accelhp, accelhc, t - 3*armdelay)) - (y13TD(splinehp, splinehc, accelhp, accelhc, t) + y31TD(splinehp, splinehc, accelhp, accelhc, t - armdelay)) - (y23TD(splinehp, splinehc, accelhp, accelhc, t - 2*armdelay) + y32TD(splinehp, splinehc, accelhp, accelhc, t - 3*armdelay));

  /* Output */
  *TDIA = 1./(2*sqrt(2)) * (Z-X);
  *TDIE = 1./(2*sqrt(6)) * (X-2*Y+Z);
  *TDIT = 1./(2*sqrt(3)) * (X+Y+Z);

  return SUCCESS;
}

/**/
int GenerateTDITD3Chanhphc(
  RealTimeSeries** TDI1,                   /* Output: real time series for TDI channel 1 */
  RealTimeSeries** TDI2,                   /* Output: real time series for TDI channel 2 */
  RealTimeSeries** TDI3,                   /* Output: real time series for TDI channel 3 */
  gsl_spline* splinehp,                    /* Input spline for TD hplus */
  gsl_spline* splinehc,                    /* Input spline for TD hcross */
  gsl_interp_accel* accelhp,               /* Accelerator for hp spline */
  gsl_interp_accel* accelhc,               /* Accelerator for hc spline */
  gsl_vector* times,                       /* Vector of times to evaluate */
  int nbptmargin,                          /* Margin set to 0 on both side to avoid problems with delays out of the domain */
  TDItag tditag)                           /* Tag selecting the TDI observables */
{
  /* Initialize output */
  int nbpt = times->size;
  RealTimeSeries_Init(TDI1, nbpt);
  RealTimeSeries_Init(TDI2, nbpt);
  RealTimeSeries_Init(TDI3, nbpt);
  gsl_vector_memcpy((*TDI1)->times, times);
  gsl_vector_memcpy((*TDI2)->times, times);
  gsl_vector_memcpy((*TDI3)->times, times);
  gsl_vector_set_zero((*TDI1)->h);
  gsl_vector_set_zero((*TDI2)->h);
  gsl_vector_set_zero((*TDI3)->h);

  /* Loop over time samples - we take a margin to avoid problems with the domain */
  double t;
  double* tval = times->data;
  double* tdi1 = (*TDI1)->h->data;
  double* tdi2 = (*TDI2)->h->data;
  double* tdi3 = (*TDI3)->h->data;
  double tdi1val = 0, tdi2val = 0, tdi3val = 0;

  /* For testing purposes: basic observable yAB */
  if(tditag==y12) {
    for(int i=nbptmargin; i<nbpt-nbptmargin; i++) {
      t = tval[i];
      tdi1[i] = y12TD(splinehp, splinehc, accelhp, accelhc, t);
      tdi2[i] = 0.;
      tdi3[i] = 0.;
    }
  }
  else if(tditag==TDIXYZ) {
    for(int i=nbptmargin; i<nbpt-nbptmargin; i++) {
      t = tval[i];
      EvaluateTDIXYZTD(&tdi1val, &tdi2val, &tdi3val, splinehp, splinehc, accelhp, accelhc, t);
      tdi1[i] = tdi1val;
      tdi2[i] = tdi2val;
      tdi3[i] = tdi3val;
    }
  }
  else if(tditag==TDIAETXYZ) {
    for(int i=nbptmargin; i<nbpt-nbptmargin; i++) {
      t = tval[i];
      EvaluateTDIAETXYZTD(&tdi1val, &tdi2val, &tdi3val, splinehp, splinehc, accelhp, accelhc, t);
      tdi1[i] = tdi1val;
      tdi2[i] = tdi2val;
      tdi3[i] = tdi3val;
    }
  }
  else {
    printf("Error: in GenerateTDITD3Chan, TDI tag not recognized.\n");
  }

  return SUCCESS;
}

/* Generate hO orbital-delayed for one mode contribution from amp, phase */
int Generateh22TDO(
  AmpPhaseTimeSeries** h22tdO,             /* Output: amp/phase time series for h22TDO */
  gsl_spline* splineamp,                   /* Input spline for TD mode amplitude */
  gsl_spline* splinephase,                 /* Input spline for TD mode phase */
  gsl_interp_accel* accelamp,              /* Accelerator for amp spline */
  gsl_interp_accel* accelphase,            /* Accelerator for phase spline */
  gsl_vector* times,                       /* Vector of times to evaluate */
  int nbptmargin)                          /* Margin set to 0 on both side to avoid problems with delays out of the domain */
{
  /* Initialize output */
  int nbpt = times->size;
  AmpPhaseTimeSeries_Init(h22tdO, nbpt);
  gsl_vector_memcpy((*h22tdO)->times, times);
  gsl_vector_set_zero((*h22tdO)->h_amp);
  gsl_vector_set_zero((*h22tdO)->h_phase);

  /* Loop over time samples - we take a margin to avoid problems with the domain */
  double t;
  double* tval = times->data;
  double* amp = (*h22tdO)->h_amp->data;
  double* phase = (*h22tdO)->h_phase->data;

  /* Loop over time samples */
  for(int i=nbptmargin; i<nbpt-nbptmargin; i++) {
    t = tval[i];
    hOTDAmpPhase(&(amp[i]), &(phase[i]), splineamp, splinephase, accelamp, accelphase, t);
  }

  return SUCCESS;
}

/* Generate y12L from orbital-delayed h22 in amp/phase form */
/* Note: includes both h22 and h2m2 contributions, assuming planar orbits so that h2-2 = h22* */
/* BEWARE: this ignores the fact that processing through orbital delay breaks the h2-2 = h22* symmetry */
int Generatey12LTD(
  RealTimeSeries** y12Ltd,                 /* Output: real time series for y12L */
  gsl_spline* splineamp,                   /* Input spline for h22 TD amplitude */
  gsl_spline* splinephase,                 /* Input spline for h22 TD phase */
  gsl_interp_accel* accelamp,              /* Accelerator for h22 amp spline */
  gsl_interp_accel* accelphase,            /* Accelerator for h22 phase spline */
  gsl_vector* times,                       /* Vector of times to evaluate */
  double Theta,                            /* Inclination */
  double Phi,                              /* Phase */
  int nbptmargin)                          /* Margin set to 0 on both side to avoid problems with delays out of the domain */
{
  /* Initialize output */
  int nbpt = times->size;
  RealTimeSeries_Init(y12Ltd, nbpt);
  gsl_vector_memcpy((*y12Ltd)->times, times);
  gsl_vector_set_zero((*y12Ltd)->h);

  /* Spin-weighted spherical harmonic Y22 and Y2-2 */
  double complex Y22 = SpinWeightedSphericalHarmonic(Theta, Phi, -2, 2, 2);
  double complex Y2m2 = SpinWeightedSphericalHarmonic(Theta, Phi, -2, 2, -2);

  /* Loop over time samples - we take a margin to avoid problems with the domain */
  double t;
  double* tval = times->data;
  double* y12val = (*y12Ltd)->h->data;

  /* Loop over time samples */
  for(int i=nbptmargin; i<nbpt-nbptmargin; i++) {
    t = tval[i];
    y12val[i] = y12LTDfromh22AmpPhase(splineamp, splinephase, accelamp, accelphase, Y22, Y2m2, t);
  }

  return SUCCESS;
}
