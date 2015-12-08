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

/*************************************************************/
/********* Functions for the geometric response **************/

/* Function cardinal sine */
double sinc(const double x) {
  if (x==0)
    return 1;
  else return sin(x)/x;
}

/* Function to compute, given a value of a sky position and polarization, all the complicated time-independent trigonometric coefficients entering the response */
void SetCoeffsG(const double lambda, const double beta, const double psi) {
  /* Precomputing cosines and sines */
  double sqrt3 = sqrt(3);
  double coslambda = cos(lambda);
  double sinlambda = sin(lambda);
  double cosbeta = cos(beta);
  double sinbeta = sin(beta);
  double cospsi = cos(psi);
  double sinpsi = sin(psi);
  //printf("cos(lambda), sin(lambda), cos(beta), sin(beta), cos(psi), sin(psi): %g, %g, %g, %g, %g, %g\n", coslambda, sinlambda, cosbeta, sinbeta, cospsi, sinpsi);

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
  coeffkn3const = 3./8*cosbeta * (-sinlambda + sqrt3*coslambda);
  /**/
  coeffkn3cos[0] = 3./4 * (sinbeta);
  /**/
  coeffkn3cos[1] = -1./8*cosbeta * (sinlambda + sqrt3*coslambda);
  /**/
  coeffkn3sin[0] = -1./4*sqrt3 * (sinbeta);
  /**/
  coeffkn3sin[1] = 1./8*cosbeta * (coslambda -sqrt3*sinlambda);

  /* Coefficients in k.n2 */
  /**/
  coeffkn2const = -3./8*cosbeta * (sinlambda + sqrt3*coslambda);
  /**/
  coeffkn2cos[0] = -3./4 * (sinbeta);
  /**/
  coeffkn2cos[1] = 1./8*cosbeta * (-sinlambda + sqrt3*coslambda);
  /**/
  coeffkn2sin[0] = -1./4*sqrt3 * (sinbeta);
  /**/
  coeffkn2sin[1] = 1./8*cosbeta * (coslambda + sqrt3*sinlambda);

  /* Coefficients in k.n1 */
  /**/
  coeffkn1const = 3./4*cosbeta * (sinlambda);
  /**/
  coeffkn1cos[0] =  0. ;
  /**/
  coeffkn1cos[1] = 1./4*cosbeta * (sinlambda);
  /**/
  coeffkn1sin[0] = 1./2*sqrt3 * (sinbeta);
  /**/
  coeffkn1sin[1] = -1./4*cosbeta * (coslambda);

  /* Coefficients in k.(p1+p2) */
  /**/
  coeffkp1plusp2const = -1./8*cosbeta * (3*sinlambda + sqrt3*coslambda);
  /**/
  coeffkp1plusp2cos[0] = -1./4 * (sinbeta);
  /**/
  coeffkp1plusp2cos[1] = 1./24*cosbeta * (-3*sinlambda + sqrt3*coslambda);
  /**/
  coeffkp1plusp2sin[0] = -1./4*sqrt3 * (sinbeta);
  /**/
  coeffkp1plusp2sin[1] = 1./24*cosbeta * (3*coslambda + sqrt3*sinlambda);

  /* Coefficients in k.(p2+p3) */
  /**/
  coeffkp2plusp3const = 1./4*sqrt3*cosbeta * (coslambda);
  /**/
  coeffkp2plusp3cos[0] = 1./2 * (sinbeta);
  /**/
  coeffkp2plusp3cos[1] = -1./4/sqrt3 * (cosbeta*coslambda);
  /**/
  coeffkp2plusp3sin[0] =  0. ;
  /**/
  coeffkp2plusp3sin[1] = -1./4/sqrt3 * (cosbeta*sinlambda);

  /* Coefficients in k.(p3+p1) */
  /**/
  coeffkp3plusp1const = -1./8*cosbeta * (-3*sinlambda + sqrt3*coslambda);
  /**/
  coeffkp3plusp1cos[0] = -1./4 * (sinbeta);
  /**/
  coeffkp3plusp1cos[1] = 1./24*cosbeta * (3*sinlambda + sqrt3*coslambda);
  /**/
  coeffkp3plusp1sin[0] = 1./4*sqrt3 * (sinbeta);
  /**/
  coeffkp3plusp1sin[1] = -1./24*cosbeta * (3*coslambda -sqrt3*sinlambda);

  /* Coefficients in k.p1 */
  /**/
  coeffkp1const = -1./4*sqrt3 * (cosbeta*coslambda);
  /**/
  coeffkp1cos[0] = -1./2 * (sinbeta);
  /**/
  coeffkp1cos[1] = 1./(4*sqrt3) * (cosbeta*coslambda);
  /**/
  coeffkp1sin[0] =  0. ;
  /**/
  coeffkp1sin[1] = 1./(4*sqrt3) * (cosbeta*sinlambda);

  /* Coefficients in k.p2 */
  /**/
  coeffkp2const = 1./8*cosbeta * (-3*sinlambda + sqrt3*coslambda);
  /**/
  coeffkp2cos[0] = 1./4 * (sinbeta);
  /**/
  coeffkp2cos[1] = -1./24*cosbeta * (3*sinlambda + sqrt3*coslambda);
  /**/
  coeffkp2sin[0] = -1./4*sqrt3 * (sinbeta);
  /**/
  coeffkp2sin[1] = 1./24*cosbeta * (3*coslambda -sqrt3*sinlambda);

  /* Coefficients in k.p3 */
  /**/
  coeffkp3const = 1./8*cosbeta * (3*sinlambda + sqrt3*coslambda);
  /**/
  coeffkp3cos[0] = 1./4 * (sinbeta);
  /**/
  coeffkp3cos[1] = -1./24*cosbeta * (-3*sinlambda + sqrt3*coslambda);
  /**/
  coeffkp3sin[0] = 1./4*sqrt3 * (sinbeta);
  /**/
  coeffkp3sin[1] = -1./24*cosbeta * (3*coslambda + sqrt3*sinlambda);

  /* Coefficients in k.R */
  /**/
  coeffkRconst = 0.;
  coeffkRcos[0] = 1. * (cosbeta*coslambda);
  coeffkRsin[0] = 1. * (cosbeta*sinlambda);
  coeffkRcos[1] = 0.;
  coeffkRsin[1] = 0.;


}

/*********************** Fourier-domain response ************************/

/* Function evaluating G21, combining the two polarization with the spherical harmonics factors */
double complex G21mode(const double f, const double t, const double complex Yfactorplus, const double complex Yfactorcross) {
  
  for(int j=0; j<4; j++) {
    cosarray[j] = cos((j+1) * Omega_SI * t);
    sinarray[j] = sin((j+1) * Omega_SI * t);
  }
  double n3Pn3plus = coeffn3Hn3plusconst;
  double n3Pn3cross = coeffn3Hn3crossconst;
  //printf("n3Pn3cross: %g\n", n3Pn3cross);
  for(int j=0; j<4; j++) {
    n3Pn3plus += cosarray[j] * coeffn3Hn3pluscos[j] + sinarray[j] * coeffn3Hn3plussin[j];
    n3Pn3cross += cosarray[j] * coeffn3Hn3crosscos[j] + sinarray[j] * coeffn3Hn3crosssin[j];
  }
  //printf("n3Pn3cross: %g\n", n3Pn3cross);
  double kn3 = coeffkn3const;
  double kp1plusp2 = coeffkp1plusp2const;
  for(int j=0; j<2; j++) {
    kn3 += cosarray[j] * coeffkn3cos[j] + sinarray[j] * coeffkn3sin[j];
    kp1plusp2 += cosarray[j] * coeffkp1plusp2cos[j] + sinarray[j] * coeffkp1plusp2sin[j];
  }
  //printf("Yfactorplus: %g+i*%g\n", creal(Yfactorplus), cimag(Yfactorplus));
  //printf("Yfactorcross: %g+i*%g\n", creal(Yfactorcross), cimag(Yfactorcross));
  //printf("%g+i*%g\n", creal((n3Pn3plus*Yfactorplus + n3Pn3cross*Yfactorcross)), cimag((n3Pn3plus*Yfactorplus + n3Pn3cross*Yfactorcross)));
  //printf("%g+i*%g\n", creal( sinc( PI*f*L_SI/C_SI * (1.+kn3))), cimag( sinc( PI*f*L_SI/C_SI * (1.+kn3))));
  //printf("%g+i*%g\n", creal(cexp( I*PI*f*L_SI/C_SI * (1.-kp1plusp2) )), cimag(cexp( I*PI*f*L_SI/C_SI * (1.-kp1plusp2) )));
  //printf("%g+i*%g\n", creal(I*PI*f*L_SI/C_SI), cimag(I*PI*f*L_SI/C_SI));
  //printf("result G21: %g+i*%g\n", creal(I*PI*f*L_SI/C_SI * (n3Pn3plus*Yfactorplus + n3Pn3cross*Yfactorcross) * sinc( PI*f*L_SI/C_SI * (1.+kn3)) * cexp( I*PI*f*L_SI/C_SI * (1.-kp1plusp2) )), cimag(I*PI*f*L_SI/C_SI * (n3Pn3plus*Yfactorplus + n3Pn3cross*Yfactorcross) * sinc( PI*f*L_SI/C_SI * (1.+kn3)) * cexp( I*PI*f*L_SI/C_SI * (1.-kp1plusp2) )));
  return I*PI*f*L_SI/C_SI * (n3Pn3plus*Yfactorplus + n3Pn3cross*Yfactorcross) * sinc( PI*f*L_SI/C_SI * (1.+kn3)) * cexp( I*PI*f*L_SI/C_SI * (1.-kp1plusp2) ); 
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

  return I*PI*f*L_SI/C_SI * (n3Pn3plus*Yfactorplus + n3Pn3cross*Yfactorcross) * sinc( PI*f*L_SI/C_SI * (1.-kn3)) * cexp( I*PI*f*L_SI/C_SI * (1.-kp1plusp2) ); 
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

  return I*PI*f*L_SI/C_SI * (n1Pn1plus*Yfactorplus + n1Pn1cross*Yfactorcross) * sinc( PI*f*L_SI/C_SI * (1.+kn1)) * cexp( I*PI*f*L_SI/C_SI * (1.-kp2plusp3) ); 
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

  return I*PI*f*L_SI/C_SI * (n1Pn1plus*Yfactorplus + n1Pn1cross*Yfactorcross) * sinc( PI*f*L_SI/C_SI * (1.-kn1)) * cexp( I*PI*f*L_SI/C_SI * (1.-kp2plusp3) ); 
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

  return I*PI*f*L_SI/C_SI * (n2Pn2plus*Yfactorplus + n2Pn2cross*Yfactorcross) * sinc( PI*f*L_SI/C_SI * (1.+kn2)) * cexp( I*PI*f*L_SI/C_SI * (1.-kp3plusp1) ); 
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

  return I*PI*f*L_SI/C_SI * (n2Pn2plus*Yfactorplus + n2Pn2cross*Yfactorcross) * sinc( PI*f*L_SI/C_SI * (1.-kn2)) * cexp( I*PI*f*L_SI/C_SI * (1.-kp3plusp1) ); 
}

/* Function evaluating all coefficients G12, G21, G23, G32, G31, G13, combining the two polarization with the spherical harmonics factors */
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
  const double complex Yfactorcross)       /* Spin-weighted spherical harmonic factor for cross */
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
  for(int j=0; j<2; j++) {
    kn1 += cosarray[j] * coeffkn1cos[j] + sinarray[j] * coeffkn1sin[j];
    kn2 += cosarray[j] * coeffkn2cos[j] + sinarray[j] * coeffkn2sin[j];
    kn3 += cosarray[j] * coeffkn3cos[j] + sinarray[j] * coeffkn3sin[j];
    kp1plusp2 += cosarray[j] * coeffkp1plusp2cos[j] + sinarray[j] * coeffkp1plusp2sin[j];
    kp2plusp3 += cosarray[j] * coeffkp2plusp3cos[j] + sinarray[j] * coeffkp2plusp3sin[j];
    kp3plusp1 += cosarray[j] * coeffkp3plusp1cos[j] + sinarray[j] * coeffkp3plusp1sin[j];
  }
  /* Common factors */
  double complex factn1Pn1 = n1Pn1plus*Yfactorplus + n1Pn1cross*Yfactorcross;
  double complex factn2Pn2 = n2Pn2plus*Yfactorplus + n2Pn2cross*Yfactorcross;
  double complex factn3Pn3 = n3Pn3plus*Yfactorplus + n3Pn3cross*Yfactorcross;
  double prefactor = PI*f*L_SI/C_SI;
  double complex factorcexp12 = cexp(I*prefactor * (1.-kp1plusp2));
  double complex factorcexp23 = cexp(I*prefactor * (1.-kp2plusp3));
  double complex factorcexp31 = cexp(I*prefactor * (1.-kp3plusp1));
  /* Output result */
  *G12 = I*prefactor * factn3Pn3 * sinc( prefactor * (1.-kn3)) * factorcexp12;
  *G21 = I*prefactor * factn3Pn3 * sinc( prefactor * (1.+kn3)) * factorcexp12;
  *G23 = I*prefactor * factn1Pn1 * sinc( prefactor * (1.-kn1)) * factorcexp23;
  *G32 = I*prefactor * factn1Pn1 * sinc( prefactor * (1.+kn1)) * factorcexp23;
  *G31 = I*prefactor * factn2Pn2 * sinc( prefactor * (1.-kn2)) * factorcexp31;
  *G13 = I*prefactor * factn2Pn2 * sinc( prefactor * (1.+kn2)) * factorcexp31;

  return SUCCESS;
}

/*********************** Time-domain response ************************/

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
  double factorp = (1./(1.+kn3)) * 0.5*n3Pn3plus;
  double factorc = (1./(1.+kn3)) * 0.5*n3Pn3cross;
  double firstdelay = (kR*R_SI + (kp1 - 1)*L_SI)/C_SI;
  double seconddelay = (kR*R_SI + kp2*L_SI)/C_SI;

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
  double factorp = (1./(1.-kn3)) * 0.5*n3Pn3plus;
  double factorc = (1./(1.-kn3)) * 0.5*n3Pn3cross;
  double firstdelay = (kR*R_SI + (kp2 - 1)*L_SI)/C_SI;
  double seconddelay = (kR*R_SI + kp1*L_SI)/C_SI;

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
  double factorp = (1./(1.+kn1)) * 0.5*n1Pn1plus;
  double factorc = (1./(1.+kn1)) * 0.5*n1Pn1cross;
  double firstdelay = (kR*R_SI + (kp2 - 1)*L_SI)/C_SI;
  double seconddelay = (kR*R_SI + kp3*L_SI)/C_SI;

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
  double factorp = (1./(1.-kn1)) * 0.5*n1Pn1plus;
  double factorc = (1./(1.-kn1)) * 0.5*n1Pn1cross;
  double firstdelay = (kR*R_SI + (kp3 - 1)*L_SI)/C_SI;
  double seconddelay = (kR*R_SI + kp2*L_SI)/C_SI;

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
  double factorp = (1./(1.+kn2)) * 0.5*n2Pn2plus;
  double factorc = (1./(1.+kn2)) * 0.5*n2Pn2cross;
  double firstdelay = (kR*R_SI + (kp3 - 1)*L_SI)/C_SI;
  double seconddelay = (kR*R_SI + kp1*L_SI)/C_SI;

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
  double factorp = (1./(1.-kn2)) * 0.5*n2Pn2plus;
  double factorc = (1./(1.-kn2)) * 0.5*n2Pn2cross;
  double firstdelay = (kR*R_SI + (kp1 - 1)*L_SI)/C_SI;
  double seconddelay = (kR*R_SI + kp3*L_SI)/C_SI;

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
