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
/* WARNING: Taken from (4) in McWilliams&al_0911, but there was a typo there */
/* WARNING: the term (1. + 1e-8*invf2) in Spm should not come with a sqrt */
static double SpmLISA2010(const double f) {
  double invf2 = 1./(f*f);
  //return 2.5e-48 * invf2 * sqrt(1. + 1e-8*invf2);
  //const double Daccel=3.0e-15; //acceleration noise in m/s^2/sqrt(Hz)
  const double Daccel=3.0e-15; //scaled off L3LISA-v1 to for equal-SNR PE experiment
  const double SaccelFF=Daccel*Daccel/4.0/PI/PI/C_SI/C_SI; //f^-2 coeff for fractional-freq noise PSD from accel noise; yields 2.54e-48 from 3e-15;
  double invf8=invf2*invf2*invf2*invf2;
  //Here we add an eyeball approximation based on 4yrs integration with L3LISAReferenceMission looking at a private comm from Neil Cornish 2016.11.12
  double WDWDnoise=5000.0/sqrt(1e-21*invf8 + invf2 + 3e28/invf8)*SaccelFF*invf2;
  return SaccelFF * invf2 * (1. + 1e-8*invf2) + WDWDnoise;
}
static double SopLISA2010(const double f) {
  const double Dop=2.0e-11; //Optical path noise in m/rtHz (Standard LISA)
  const double SopFF=Dop*Dop*4.0*PI*PI/C_SI/C_SI; //f^2 coeff for OP frac-freq noise PSD.  Yields 1.76e-37 for Dop=2e-11.
  return SopFF * f * f;
}

/* Proof mass and optical noises - f in Hz */
/* L3 Reference Mission, from Petiteau LISA-CST-TN-0001 */
static double SpmLISA2017(const double f) {
  double invf2 = 1./(f*f);
  //double invf4=invf2*invf2;
  //double invf8=invf4*invf4;
  //double invf10=invf8*invf2;
  const double twopi2=4.0*PI*PI;
  double ddtsq=twopi2/invf2; //time derivative factor
  const double C2=1.0*C_SI*C_SI; //veloc to doppler
  const double Daccel_white=3.0e-15; //acceleration noise in m/s^2/sqrt(Hz)
  const double Daccel_white2=Daccel_white*Daccel_white;
  const double Dloc=1.7e-12; //local IFO noise in m/sqrt(Hz)
  const double Dloc2=Dloc*Dloc;
  double Saccel_white=Daccel_white2/ddtsq; //PM vel noise PSD (white accel part)
  //double Saccel_red=Saccel_white*(1.0 + 2.12576e-44*invf10 + 3.6e-7*invf2); //reddening factor from Petiteau Eq 1
  double Saccel_red=Saccel_white*(1.0 + 36.0*(pow(3e-5/f,10) + 1e-8*invf2)); //reddening factor from Petiteau Eq 1
  //Saccel_red*=4.0;//Hack to decrease low-freq sens by fac of 2.
  double Sloc=Dloc2*ddtsq/4.0;//Factor of 1/4.0 is in Petiteau eq 2
  double S4yrWDWD=5.16e-27*exp(-pow(f,1.2)*2.9e3)*pow(f,(-7./3.))*0.5*(1.0 + tanh(-(f-2.0e-3)*1.9e3))*ddtsq;//Stas' fit for 4yr noise (converted from Sens curve to position noise by multipyling by 3*L^2/80) which looks comparable to my fit), then converted to velocity noise
  double Spm_vel = ( Saccel_red + Sloc + S4yrWDWD );
  return Spm_vel / C2;//finally convert from velocity noise to fractional-frequency doppler noise.
}
static double SopLISA2017(const double f) {
  //double invf2 = 1./(f*f);
  const double twopi2=4.0*PI*PI;
  double ddtsq=twopi2*f*f; //time derivative factor
  const double C2=C_SI*C_SI; //veloc to doppler
  const double Dloc=1.7e-12; //local IFO noise in m/sqrt(Hz)
  const double Dsci=8.9e-12; //science IFO noise in m/sqrt(Hz)
  const double Dmisc=2.0e-12; //misc. optical path noise in m/sqrt(Hz)
  const double Dop2=Dsci*Dsci+Dloc*Dloc+Dmisc*Dmisc;
  double Sop=Dop2*ddtsq/C2; //f^2 coeff for OP frac-freq noise PSD.  Yields 1.76e-37 for Dop=2e-11.
  return Sop;
}

/* Proof mass and optical noises - f in Hz */
/* LISA Proposal, copied from the LISA Data Challenge pipeline */
static double SpmLISAProposal(const double f) {
  /* Acceleration noise */
  double noise_Sa_a = 9.e-30; /* m^2/sec^4 /Hz */
  /* In acceleration */
  double Sa_a = noise_Sa_a * (1.0 + pow(0.4e-3/f, 2)) * (1.0 + pow((f/8e-3), 4));
  /* In displacement */
  double Sa_d = Sa_a * pow(2.*PI*f, -4);
  /* In relative frequency unit */
  double Sa_nu = Sa_d * pow(2.*PI*f/C_SI, 2);
  double Spm = Sa_nu;
  return Spm;
}
static double SopLISAProposal(const double f) {
  /* Optical Metrology System noise */
  double noise_Soms_d = pow((10e-12), 2); /* m^2/Hz */
  /* In displacement */
  double Soms_d = noise_Soms_d * (1. + pow(2.e-3/f, 4));
  /* In relative frequency unit */
  double Soms_nu = Soms_d * pow(2.*PI*f/C_SI, 2);
  double Sop = Soms_nu;
  return Sop;
}

/* Compute proof mass and optical noises, for a given choice of noise - f in Hz */
static void ComputeLISASpmSop(double* Spm, double* Sop, const double f, LISANoiseType noise) {
  if(noise==LISAProposalnoise) {
    *Spm = SpmLISAProposal(f);
    *Sop = SopLISAProposal(f);
  }
  else if(noise==LISA2017noise) {
    *Spm = SpmLISA2017(f);
    *Sop = SopLISA2017(f);
  }
  else if(noise==LISA2010noise) {
    *Spm = SpmLISA2010(f);
    *Sop = SopLISA2010(f);
  }
  else {
    printf("Error in ComputeLISASpmSop: LISANoiseType not recognized.\n");
    exit(1);
  }
}

/* Noise Sn for TDI observables - factors have been scaled out both in the response and the noise */
/* Rescaled by 4*sin2pifL^2 */
double SnXYZ(const LISAconstellation *variant, double f) {
  double twopifL = 2.*PI*variant->ConstL/C_SI*f;
  double c2 = cos(twopifL);
  double Spm = 0., Sop = 0.;
  ComputeLISASpmSop(&Spm, &Sop, f, variant->noise);
  return 4*( 2*(1. + c2*c2)*Spm + Sop );
}
/* No rescaling */
double Snalphabetagamma(const LISAconstellation *variant, double f) {
  double pifL = PI*variant->ConstL/C_SI*f;
  double s1 = sin(pifL);
  double s3 = sin(3*pifL);
  double Spm = 0., Sop = 0.;
  ComputeLISASpmSop(&Spm, &Sop, f, variant->noise);
  return 2*( (4*s3*s3 + 8*s1*s1)*Spm + 3*Sop );
}
/* Rescaled by 2*sin2pifL^2 */
double SnAXYZ(const LISAconstellation *variant, double f) {
  double twopifL = 2.*PI*variant->ConstL/C_SI*f;
  double c2 = cos(twopifL);
  double c4 = cos(2*twopifL);
  double Spm = 0., Sop = 0.;
  ComputeLISASpmSop(&Spm, &Sop, f, variant->noise);
  return 2*(3. + 2*c2 + c4)*Spm + (2 + c2)*Sop;
}
/* Rescaled by 2*sin2pifL^2 */
double SnEXYZ(const LISAconstellation *variant, double f) {
  double twopifL = 2.*PI*variant->ConstL/C_SI*f;
  double c2 = cos(twopifL);
  double c4 = cos(2*twopifL);
  double Spm = 0., Sop = 0.;
  ComputeLISASpmSop(&Spm, &Sop, f, variant->noise);
  return 2*(3. + 2*c2 + c4)*Spm + (2 + c2)*Sop;
}
/* Rescaled by 8*sin2pifL^2*sinpifL^2 */
double SnTXYZ(const LISAconstellation *variant, double f) {
  double pifL = PI*variant->ConstL/C_SI*f;
  double s1 = sin(pifL);
  double Spm = 0., Sop = 0.;
  ComputeLISASpmSop(&Spm, &Sop, f, variant->noise);
  return 4*s1*s1*Spm + Sop;
}
/* Rescaled by 8*sin2pifL^2 */
double SnAalphabetagamma(const LISAconstellation *variant, double f) {
  double twopifL = 2.*PI*variant->ConstL/C_SI*f;
  double c2 = cos(twopifL);
  double c4 = cos(2*twopifL);
  double Spm = 0., Sop = 0.;
  ComputeLISASpmSop(&Spm, &Sop, f, variant->noise);
  return 2*(3. + 2*c2 + c4)*Spm + (2 + c2)*Sop;
}
/* Rescaled by 8*sin2pifL^2 */
double SnEalphabetagamma(const LISAconstellation *variant, double f) {
  double twopifL = 2.*PI*variant->ConstL/C_SI*f;
  double c2 = cos(twopifL);
  double c4 = cos(2*twopifL);
  double Spm = 0., Sop = 0.;
  ComputeLISASpmSop(&Spm, &Sop, f, variant->noise);
  return 2*(3. + 2*c2 + c4)*Spm + (2 + c2)*Sop;
}
/* Rescaled by sin3pifL^2/sinpifL^2 */
double SnTalphabetagamma(const LISAconstellation *variant, double f) {
  double pifL = PI*variant->ConstL/C_SI*f;
  double s1 = sin(pifL);
  double Spm = 0., Sop = 0.;
  ComputeLISASpmSop(&Spm, &Sop, f, variant->noise);
  return 8*s1*s1*Spm + 2*Sop;
}

/* Noise functions for AET(XYZ) without rescaling */
/* Scaling by 2*sin2pifL^2 put back */
double SnAXYZNoRescaling(const LISAconstellation *variant, double f) {
  double twopifL = 2.*PI*variant->ConstL/C_SI*f;
  double c2 = cos(twopifL);
  double c4 = cos(2*twopifL);
  double s2 = sin(twopifL);
  double Spm = 0., Sop = 0.;
  ComputeLISASpmSop(&Spm, &Sop, f, variant->noise);
  return 2*s2*s2 * (2*(3. + 2*c2 + c4)*Spm + (2 + c2)*Sop);
}
/* Scaling by 2*sin2pifL^2 put back */
double SnEXYZNoRescaling(const LISAconstellation *variant, double f) {
  double twopifL = 2.*PI*variant->ConstL/C_SI*f;
  double c2 = cos(twopifL);
  double c4 = cos(2*twopifL);
  double s2 = sin(twopifL);
  double Spm = 0., Sop = 0.;
  ComputeLISASpmSop(&Spm, &Sop, f, variant->noise);
  return 2*s2*s2 * (2*(3. + 2*c2 + c4)*Spm + (2 + c2)*Sop);
}
/* Scaling by 8*sin2pifL^2*sinpifL^2 put back*/
double SnTXYZNoRescaling(const LISAconstellation *variant, double f) {
  double pifL = PI*variant->ConstL/C_SI*f;
  double s1 = sin(pifL);
  double s2 = sin(2*pifL);
  double Spm = 0., Sop = 0.;
  ComputeLISASpmSop(&Spm, &Sop, f, variant->noise);
  return 8*s1*s1*s2*s2 * (4*s1*s1*Spm + Sop);
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
ObjectFunction NoiseFunction(const LISAconstellation *variant, const TDItag tditag, const int nchan)
{
  ObjectFunction fn;
  switch(tditag) {
  case TDIXYZ:
  case TDIX: {
    switch(nchan) {
    case 1: fn = (ObjectFunction){variant,(RealObjectFunctionPtr)SnXYZ}; break;
    case 2: fn = (ObjectFunction){variant,(RealObjectFunctionPtr)SnXYZ}; break;
    case 3: fn = (ObjectFunction){variant,(RealObjectFunctionPtr)SnXYZ}; break;
    }
    break;
  }
  case TDIalphabetagamma:
  case TDIalpha: {
    switch(nchan) {
    case 1: fn = (ObjectFunction){variant,(RealObjectFunctionPtr)Snalphabetagamma}; break;
    case 2: fn = (ObjectFunction){variant,(RealObjectFunctionPtr)Snalphabetagamma}; break;
    case 3: fn = (ObjectFunction){variant,(RealObjectFunctionPtr)Snalphabetagamma}; break;
    }
    break;
  }
  case TDIAETXYZ:
  case TDIAXYZ:
  case TDIEXYZ:
  case TDITXYZ: {
    switch(nchan) {
    case 1: fn = (ObjectFunction){variant,(RealObjectFunctionPtr)SnAXYZ}; break;
    case 2: fn = (ObjectFunction){variant,(RealObjectFunctionPtr)SnEXYZ}; break;
    case 3: fn = (ObjectFunction){variant,(RealObjectFunctionPtr)SnTXYZ}; break;
    }
    break;
  }
  case TDIAETalphabetagamma:
  case TDIAalphabetagamma:
  case TDIEalphabetagamma:
  case TDITalphabetagamma: {
    switch(nchan) {
    case 1: fn = (ObjectFunction){variant,(RealObjectFunctionPtr)SnAalphabetagamma}; break;
    case 2: fn = (ObjectFunction){variant,(RealObjectFunctionPtr)SnEalphabetagamma}; break;
    case 3: fn = (ObjectFunction){variant,(RealObjectFunctionPtr)SnTalphabetagamma}; break;
    }
    break;
  }
  }
  if(fn.object==NULL) {
    printf("Error in NoiseFunction: incorrect argument.\n");
    exit(1);
  }
  return fn;
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
