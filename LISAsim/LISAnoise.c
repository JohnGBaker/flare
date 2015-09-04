/**
 * \author Sylvain Marsat, University of Maryland - NASA GSFC
 *
 * \brief C code for the instrumental noise for LISA-type detectors.
 *
 * Formulas taken from Królak&al gr-qc/0401108 (c.f. section III).
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
#include "LISAnoise.h"


/* static double tflight_SI = L_SI/C_SI; */
/* static double twopitflight_SI = 2.*PI*L_SI/C_SI; */

/* Proof mass and optic noises - f in Hz */
static double Spm(const double f) {
  return 2.54e-48 /f/f;
}
static double Sop(const double f) {
  return 1.76e-37 *f*f;
}

/* The noise functions themselves */
double NoiseSnA(const double f) {
  double twopifL = 2.*PI*L_SI/C_SI*f;
  double sinhalf = sin(1./2*twopifL);
  double sin3half = sin(3./2*twopifL);
  double cos1 = cos(twopifL);
  double cos2 = cos(2*twopifL);
  return 32*sinhalf*sinhalf*sin3half*sin3half*( (6 + 4*cos1 + 2*cos2)*Spm(f) + (2 + cos1)*Sop(f) );
}
double NoiseSnE(const double f) {
  double twopifL = 2.*PI*L_SI/C_SI*f;
  double sinhalf = sin(1./2*twopifL);
  double sin3half = sin(3./2*twopifL);
  double cos1 = cos(twopifL);
  double cos2 = cos(2*twopifL);
  return 32*sinhalf*sinhalf*sin3half*sin3half*( (6 + 4*cos1 + 2*cos2)*Spm(f) + (2 + cos1)*Sop(f) );
}
double NoiseSnT(const double f) {
  double twopifL = 2.*PI*L_SI/C_SI*f;
  double sinhalf = sin(1./2*twopifL);
  double sin3half = sin(3./2*twopifL);
  double cos1 = cos(twopifL);
  return 8*(1+2*cos1)*(1+2*cos1)*sin3half*sin3half*( 4*sinhalf*sinhalf*Spm(f) + Sop(f) );
}
