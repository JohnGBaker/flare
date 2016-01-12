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
/* Taken from (4) in McWilliams&al_0911 */
static double Spm(const double f) {
  double invf2 = 1./(f*f);
  return 2.5e-48 * invf2 * sqrt(1. + 1e-8*invf2);
}
static double Sop(const double f) {
  return 1.8e-37 *f*f;
}

/* Noise Sn for TDI observables - factors have been scaled out both in the response and the noise */
/* Rescaled by 4*sin2pifL^2 */
double SnXYZ(double f) {
  double twopifL = 2.*PI*L_SI/C_SI*f;
  double c2 = cos(twopifL);
  return 4*( 2*(1. + c2*c2)*Spm(f) + Sop(f) );
}
/* No rescaling */
double Snalphabetagamma(double f) {
  double pifL = PI*L_SI/C_SI*f;
  double s1 = sin(pifL);
  double s3 = sin(3*pifL);
  return 2*( (4*s3*s3 + 8*s1*s1)*Spm(f) + 3*Sop(f) );
}
/* Rescaled by 2*sin2pifL^2 */
double SnAXYZ(double f) {
  double twopifL = 2.*PI*L_SI/C_SI*f;
  double c2 = cos(twopifL);
  double c4 = cos(2*twopifL);
  return 2*(3. + 2*c2 + c4)*Spm(f) + (2 + c2)*Sop(f);
}
/* Rescaled by 2*sin2pifL^2 */
double SnEXYZ(double f) {
  double twopifL = 2.*PI*L_SI/C_SI*f;
  double c2 = cos(twopifL);
  double c4 = cos(2*twopifL);
  return 2*(3. + 2*c2 + c4)*Spm(f) + (2 + c2)*Sop(f);
}
/* Rescaled by 8*sin2pifL^2*sinpifL^2 */
double SnTXYZ(double f) {
  double pifL = PI*L_SI/C_SI*f;
  double s1 = sin(pifL);
  return 4*s1*s1*Spm(f) + Sop(f);
}
/* Rescaled by 8*sin2pifL^2 */
double SnAalphabetagamma(double f) {
  double twopifL = 2.*PI*L_SI/C_SI*f;
  double c2 = cos(twopifL);
  double c4 = cos(2*twopifL);
  return 2*(3. + 2*c2 + c4)*Spm(f) + (2 + c2)*Sop(f);
}
/* Rescaled by 8*sin2pifL^2 */
double SnEalphabetagamma(double f) {
  double twopifL = 2.*PI*L_SI/C_SI*f;
  double c2 = cos(twopifL);
  double c4 = cos(2*twopifL);
  return 2*(3. + 2*c2 + c4)*Spm(f) + (2 + c2)*Sop(f);
}
/* Rescaled by sin3pifL^2/sinpifL^2 */
double SnTalphabetagamma(double f) {
  double pifL = PI*L_SI/C_SI*f;
  double s1 = sin(pifL);
  return 8*s1*s1*Spm(f) + 2*Sop(f);
}

/* The noise functions themselves
/* Note - we factored out and cancelled the factors of the type sin(n pi f L) */
/* double NoiseSnA(const double f) { */
/*   double twopifL = 2.*PI*L_SI/C_SI*f; */
/*   double cos1 = cos(twopifL); */
/*   double cos2 = cos(2*twopifL); */
/*   return 32*( (6 + 4*cos1 + 2*cos2)*Spm(f) + (2 + cos1)*Sop(f) ); */
/* } */
/* double NoiseSnE(const double f) { */
/*   double twopifL = 2.*PI*L_SI/C_SI*f; */
/*   double cos1 = cos(twopifL); */
/*   double cos2 = cos(2*twopifL); */
/*   return 32*( (6 + 4*cos1 + 2*cos2)*Spm(f) + (2 + cos1)*Sop(f) ); */
/* } */
/* double NoiseSnT(const double f) { */
/*   double twopifL = 2.*PI*L_SI/C_SI*f; */
/*   double sinhalf = sin(0.5*twopifL); */
/*   double cos1 = cos(twopifL); */
/*   return 8*( 4*sinhalf*sinhalf*Spm(f) + Sop(f) ); */
/* } */

/* Function returning the relevant noise function, given a set of TDI observables and a channel */
RealFunctionPtr NoiseFunction(const TDItag tditag, const int nchan)
{
  RealFunctionPtr ptr = NULL;
  switch(tditag) {
  case TDIXYZ: {
    switch(nchan) {
    case 1: ptr = &SnXYZ; break;
    case 2: ptr = &SnXYZ; break;
    case 3: ptr = &SnXYZ; break;
    }
    break;
  }
  case TDIalphabetagamma: {
    switch(nchan) {
    case 1: ptr = &Snalphabetagamma; break;
    case 2: ptr = &Snalphabetagamma; break;
    case 3: ptr = &Snalphabetagamma; break;
    }
    break;
  }
  case TDIAETXYZ: {
    switch(nchan) {
    case 1: ptr = &SnAXYZ; break;
    case 2: ptr = &SnEXYZ; break;
    case 3: ptr = &SnTXYZ; break;
    }
    break;
  }
  case TDIAETalphabetagamma: {
    switch(nchan) {
    case 1: ptr = &SnAalphabetagamma; break;
    case 2: ptr = &SnEalphabetagamma; break;
    case 3: ptr = &SnTalphabetagamma; break;
    }
    break;
  }
  case TDIX: {
    switch(nchan) {
    case 1: ptr = &SnXYZ; break;
    }
    break;
  }
  case TDIalpha: {
    switch(nchan) {
    case 1: ptr = &Snalphabetagamma; break;
    }
    break;
  }
  case TDIAXYZ: {
    switch(nchan) {
    case 1: ptr = &SnAXYZ; break;
    }
  }
  case TDIEXYZ: {
    switch(nchan) {
    case 1: ptr = &SnEXYZ; break;
    }
    break;
  }
  case TDITXYZ: {
    switch(nchan) {
    case 1: ptr = &SnTXYZ; break;
    }
    break;
  }
  case TDIAalphabetagamma: {
    switch(nchan) {
    case 1: ptr = &SnAalphabetagamma; break;
    }
    break;
  }
  case TDIEalphabetagamma: {
    switch(nchan) {
    case 1: ptr = &SnEalphabetagamma; break;
    }
    break;
  }
  case TDITalphabetagamma: {
    switch(nchan) {
    case 1: ptr = &SnTalphabetagamma; break;
    }
    break;
  }
  }
  if(ptr==NULL) {
    printf("Error in NoiseFunction: incorrect argument.\n");
    exit(1);
  }
  return ptr;
}

//Previous version - we had put a noise floor to mitigate cancellation lines
/* double NoiseSnA(const double f) { */
/*   double twopifL = 2.*PI*L_SI/C_SI*f; */
/*   double sinhalf = sin(0.5*twopifL); */
/*   double sin3half = sin(1.5*twopifL); */
/*   double cos1 = cos(twopifL); */
/*   double cos2 = cos(2*twopifL); */
/*   double res = 32*sinhalf*sinhalf*sin3half*sin3half*( (6 + 4*cos1 + 2*cos2)*Spm(f) + (2 + cos1)*Sop(f) ); */
/*   return fmax(res, 1e-46); */
/* } */
/* double NoiseSnE(const double f) { */
/*   double twopifL = 2.*PI*L_SI/C_SI*f; */
/*   double sinhalf = sin(0.5*twopifL); */
/*   double sin3half = sin(1.5*twopifL); */
/*   double cos1 = cos(twopifL); */
/*   double cos2 = cos(2*twopifL); */
/*   double res = 32*sinhalf*sinhalf*sin3half*sin3half*( (6 + 4*cos1 + 2*cos2)*Spm(f) + (2 + cos1)*Sop(f) ); */
/*   return fmax(res, 1e-46); */
/* } */
/* double NoiseSnT(const double f) { */
/*   double twopifL = 2.*PI*L_SI/C_SI*f; */
/*   double sinhalf = sin(0.5*twopifL); */
/*   double sin3half = sin(1.5*twopifL); */
/*   double cos1 = cos(twopifL); */
/*   double res = 8*(1+2*cos1)*(1+2*cos1)*sin3half*sin3half*( 4*sinhalf*sinhalf*Spm(f) + Sop(f) ); */
/*   return fmax(res, 1e-46); */
/* } */
