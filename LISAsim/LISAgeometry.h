/**
 * \author Sylvain Marsat, University of Maryland - NASA GSFC
 *
 * \brief C code headers for the geometric coefficients entering the response for LISA-like detectors.
 *
 */

#ifndef _LISAGEOMETRY_H
#define _LISAGEOMETRY_H

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

#if defined(__cplusplus)
#define complex _Complex
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/*****************************************************/
/**************** TDI variables **********************/

typedef enum TDItag {
  delayO,                  /* Orbital delay */
  y12L,                    /* Constellation-only response y12L */
  y12,                     /* Complete response y12 (includes orbital delay) */
  TDIXYZ,
  TDIalphabetagamma,
  TDIAETXYZ,
  TDIAETalphabetagamma,
  TDIX,
  TDIalpha,
  TDIAXYZ,
  TDIEXYZ,
  TDITXYZ,
  TDIAalphabetagamma,
  TDIEalphabetagamma,
  TDITalphabetagamma
} TDItag;

/* Enumerator to choose what level of low-f approximation to do in the response */
typedef enum ResponseApproxtag {
  full,
  lowfL,
  lowf
} ResponseApproxtag;

/**************************************************/
/************** LISAconstellation *****************/

typedef enum LISANoiseType{
  LISA2010noise,
  LISA2017noise,
  LISAProposalnoise
} LISANoiseType;

typedef struct tagLISAconstellation {
  double OrbitOmega;
  double OrbitPhi0;
  double OrbitR;
  double ConstOmega;
  double ConstPhi0;
  double ConstL;
  LISANoiseType noise;
} LISAconstellation;
extern LISAconstellation LISAProposal;
extern LISAconstellation LISA2017;
extern LISAconstellation LISA2010;
extern LISAconstellation slowOrbitLISA;
extern LISAconstellation tinyOrbitLISA;
extern LISAconstellation fastOrbitLISA;
extern LISAconstellation bigOrbitLISA;


/**************************************************/
/**************** Prototypes **********************/

/* Function to convert string input TDI string to TDItag */
TDItag ParseTDItag(char* string);

/* Function to convert string input ResponseApprox to tag */
ResponseApproxtag ParseResponseApproxtag(char* string);

/* Function cardinal sine */
double sinc(const double x);

/* Compute Solar System Barycenter time tSSB from retarded time at the center of the LISA constellation tL */
/* NOTE: depends on the sky position given in SSB parameters */
double tSSBfromLframe(const LISAconstellation *variant, const double tL, const double lambdaSSB, const double betaSSB);
/* Compute retarded time at the center of the LISA constellation tL from Solar System Barycenter time tSSB */
double tLfromSSBframe(const LISAconstellation *variant, const double tSSB, const double lambdaSSB, const double betaSSB);
/* Convert L-frame params to SSB-frame params */
/* NOTE: no transformation of the phase -- approximant-dependence with e.g. EOBNRv2HMROM setting phiRef at fRef, and freedom in definition */
int ConvertLframeParamsToSSBframe(
  double* tSSB,
  double* lambdaSSB,
  double* betaSSB,
  double* polSSB,
  const double tL,
  const double lambdaL,
  const double betaL,
  const double psiL,
  const LISAconstellation *variant);
/* Convert SSB-frame params to L-frame params */
/* NOTE: no transformation of the phase -- approximant-dependence with e.g. EOBNRv2HMROM setting phiRef at fRef, and freedom in definition */
int ConvertSSBframeParamsToLframe(
  double* tL,
  double* lambdaL,
  double* betaL,
  double* polL,
  const double tSSB,
  const double lambdaSSB,
  const double betaSSB,
  const double psiSSB,
  const LISAconstellation *variant);

/* Function to compute, given a value of a sky position and polarization, all the complicated time-independent trigonometric coefficients entering the response */
void SetCoeffsG(const double lambda, const double beta, const double psi);

/* Functions evaluating the G_AB functions, combining the two polarization with the spherical harmonics factors */
double complex G21mode(const LISAconstellation *LISAvariant, const double f, const double t, const double complex Yfactorplus, const double complex Yfactorcross);
double complex G12mode(const LISAconstellation *LISAvariant, const double f, const double t, const double complex Yfactorplus, const double complex Yfactorcross);
double complex G32mode(const LISAconstellation *LISAvariant, const double f, const double t, const double complex Yfactorplus, const double complex Yfactorcross);
double complex G23mode(const LISAconstellation *LISAvariant, const double f, const double t, const double complex Yfactorplus, const double complex Yfactorcross);
double complex G13mode(const LISAconstellation *LISAvariant, const double f, const double t, const double complex Yfactorplus, const double complex Yfactorcross);
double complex G31mode(const LISAconstellation *LISAvariant, const double f, const double t, const double complex Yfactorplus, const double complex Yfactorcross);
int EvaluateGABmode(
  const LISAconstellation *LISAvariant,
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
  const ResponseApproxtag responseapprox); /* Tag to select possible low-f approximation level in FD response */

/* Functions evaluating the Fourier-domain factors (combinations of the GAB's) for TDI observables */
/* NOTE: factors have been scaled out, in parallel of what is done for the noise function */
/* Note: in case only one channel is considered, amplitudes for channels 2 and 3 are simply set to 0 */
/* (allows minimal changes from the old structure that assumed KTV A,E,T - but probably not optimal) */
int EvaluateTDIfactor3Chan(
  const LISAconstellation *variant,    /* Description of LISA variant */
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
  const TDItag tditag,                           /* Selector for the TDI observables */
  const ResponseApproxtag responseapprox);       /* Tag to select possible low-f approximation level in FD response */
/* int EvaluateTDIfactor1Chan( */
/*   double complex* factor,                       /\* Output for factor for TDI channel *\/ */
/*   const double complex G12,                      /\* Input for G12 *\/ */
/*   const double complex G21,                      /\* Input for G21 *\/ */
/*   const double complex G23,                      /\* Input for G23 *\/ */
/*   const double complex G32,                      /\* Input for G32 *\/ */
/*   const double complex G31,                      /\* Input for G31 *\/ */
/*   const double complex G13,                      /\* Input for G13 *\/ */
/*   const double f,                                /\* Frequency *\/ */
/*   const TDItag tditag);                          /\* Selector for the TDI observable *\/ */
/* Function evaluating the Fourier-domain factors that have been scaled out of TDI observables */
/* The factors scaled out, parallel what is done for the noise functions */
/* Note: in case only one channel is considered, factors for channels 2 and 3 are simply set to 0 */
int ScaledTDIfactor3Chan(
  const LISAconstellation *variant,    /* Description of LISA variant */
  double complex* factor1,                       /* Output for factor for TDI factor 1 */
  double complex* factor2,                       /* Output for factor for TDI factor 2 */
  double complex* factor3,                       /* Output for factor for TDI factor 3 */
  const double f,                                /* Frequency */
  const TDItag tditag);                          /* Selector for the TDI observables */
/* Function restoring the factor that have been scaled out of the TDI observables */
/* NOTE: the operation is made in-place, and the input is overwritten */
int RestoreInPlaceScaledFactorTDI(
  const LISAconstellation *variant,    /* Description of LISA variant */
  ListmodesCAmpPhaseFrequencySeries* listtdi,     /* Output/Input: list of mode contributions to TDI observable */
  TDItag tditag,                                  /* Tag selecting the TDI observable */
  int nchannel);                                  /* TDI channel number */

/* Functions for the response in time domain */

/* Basic yslr observables (including orbital delay) from hplus, hcross */
double y12TD(
  const LISAconstellation *variant,    /* Description of LISA variant */
  gsl_spline* splinehp,                    /* Input spline for TD hplus */
  gsl_spline* splinehc,                    /* Input spline for TD hcross */
  gsl_interp_accel* accelhp,               /* Accelerator for hp spline */
  gsl_interp_accel* accelhc,               /* Accelerator for hc spline */
  const double t);                         /* Time */
double y21TD(
  const LISAconstellation *variant,    /* Description of LISA variant */
  gsl_spline* splinehp,                    /* Input spline for TD hplus */
  gsl_spline* splinehc,                    /* Input spline for TD hcross */
  gsl_interp_accel* accelhp,               /* Accelerator for hp spline */
  gsl_interp_accel* accelhc,               /* Accelerator for hc spline */
  const double t);                         /* Time */
double y23TD(
  const LISAconstellation *variant,    /* Description of LISA variant */
  gsl_spline* splinehp,                    /* Input spline for TD hplus */
  gsl_spline* splinehc,                    /* Input spline for TD hcross */
  gsl_interp_accel* accelhp,               /* Accelerator for hp spline */
  gsl_interp_accel* accelhc,               /* Accelerator for hc spline */
  const double t);                         /* Time */
double y32TD(
  const LISAconstellation *variant,    /* Description of LISA variant */
  gsl_spline* splinehp,                    /* Input spline for TD hplus */
  gsl_spline* splinehc,                    /* Input spline for TD hcross */
  gsl_interp_accel* accelhp,               /* Accelerator for hp spline */
  gsl_interp_accel* accelhc,               /* Accelerator for hc spline */
  const double t);                         /* Time */
double y31TD(
  const LISAconstellation *variant,    /* Description of LISA variant */
  gsl_spline* splinehp,                    /* Input spline for TD hplus */
  gsl_spline* splinehc,                    /* Input spline for TD hcross */
  gsl_interp_accel* accelhp,               /* Accelerator for hp spline */
  gsl_interp_accel* accelhc,               /* Accelerator for hc spline */
  const double t);                         /* Time */
double y13TD(
  const LISAconstellation *variant,    /* Description of LISA variant */
  gsl_spline* splinehp,                    /* Input spline for TD hplus */
  gsl_spline* splinehc,                    /* Input spline for TD hcross */
  gsl_interp_accel* accelhp,               /* Accelerator for hp spline */
  gsl_interp_accel* accelhc,               /* Accelerator for hc spline */
  const double t);                         /* Time */
/* TDI observables (including orbital delay) from hplus, hcross */
int EvaluateTDIXYZTDhphc(
  const LISAconstellation *variant,    /* Description of LISA variant */
  double* TDIX,                            /* Output: value of TDI observable X */
  double* TDIY,                            /* Output: value of TDI observable Y */
  double* TDIZ,                            /* Output: value of TDI observable Z */
  gsl_spline* splinehp,                    /* Input spline for TD hplus */
  gsl_spline* splinehc,                    /* Input spline for TD hcross */
  gsl_interp_accel* accelhp,               /* Accelerator for hp spline */
  gsl_interp_accel* accelhc,               /* Accelerator for hc spline */
  const double t);                         /* Time */
int EvaluateTDIAETXYZTDhphc(
  const LISAconstellation *variant,    /* Description of LISA variant */
  double* TDIA,                            /* Output: value of TDI observable X */
  double* TDIE,                            /* Output: value of TDI observable Y */
  double* TDIT,                            /* Output: value of TDI observable Z */
  gsl_spline* splinehp,                    /* Input spline for TD hplus */
  gsl_spline* splinehc,                    /* Input spline for TD hcross */
  gsl_interp_accel* accelhp,               /* Accelerator for hp spline */
  gsl_interp_accel* accelhc,               /* Accelerator for hc spline */
  const double t);                         /* Time */

/* Generate hO orbital-delayed for one mode contribution from amp, phase */
int Generateh22TDO(
  const LISAconstellation *variant,    /* Description of LISA variant */
  AmpPhaseTimeSeries** h22tdO,             /* Output: amp/phase time series for h22TDO */
  gsl_spline* splineamp,                   /* Input spline for TD mode amplitude */
  gsl_spline* splinephase,                 /* Input spline for TD mode phase */
  gsl_interp_accel* accelamp,              /* Accelerator for amp spline */
  gsl_interp_accel* accelphase,            /* Accelerator for phase spline */
  gsl_vector* times,                       /* Vector of times to evaluate */
  int nbptmargin);                         /* Margin set to 0 on both side to avoid problems with delays out of the domain */
/* Generate y12L from orbital-delayed h22 in amp/phase form */
int Generatey12LTD(
  const LISAconstellation *variant,    /* Description of LISA variant */
  RealTimeSeries** y12Ltd,                 /* Output: real time series for y12L */
  gsl_spline* splineamp,                   /* Input spline for TD mode amplitude */
  gsl_spline* splinephase,                 /* Input spline for TD mode phase */
  gsl_interp_accel* accelamp,              /* Accelerator for amp spline */
  gsl_interp_accel* accelphase,            /* Accelerator for phase spline */
  gsl_vector* times,                       /* Vector of times to evaluate */
  double Theta,                            /* Inclination */
  double Phi,                              /* Phase */
  int nbptmargin);                         /* Margin set to 0 on both side to avoid problems with delays out of the domain */


/* Generate TDI observables (including orbital delay) for one mode contritbution from amp, phase */
/* NOTE: developed for testing purposes, so far supports only dO and y12L */
int GenerateTDITD3Chanhlm(
  AmpPhaseTimeSeries** hlm,                /* Output: real time series for TDI channel 1 */
  RealTimeSeries** TDI2,                   /* Output: real time series for TDI channel 2 */
  RealTimeSeries** TDI3,                   /* Output: real time series for TDI channel 3 */
  gsl_spline* splineamp,                   /* Input spline for TD mode amplitude */
  gsl_spline* splinephase,                 /* Input spline for TD mode phase */
  gsl_interp_accel* accelamp,              /* Accelerator for amp spline */
  gsl_interp_accel* accelphase,            /* Accelerator for phase spline */
  gsl_vector* times,                       /* Vector of times to evaluate */
  int nbptsmargin,                         /* Margin set to 0 on both side to avoid problems with delays out of the domain */
  double theta,                            /* Inclination angle - used to convert hlm to hplus, hcross for y12L - ignored for dO */
  double phi,                              /* Observer phase - used to convert hlm to hplus, hcross for y12L - ignored for dO */
  TDItag tditag);                          /* Tag selecting the TDI observables */


/* Generate TDI observables (including orbital delay) for one mode contritbution from hplus, hcross */
int GenerateTDITD3Chanhphc(
  const LISAconstellation *variant,    /* Description of LISA variant */
  RealTimeSeries** TDI1,                   /* Output: real time series for TDI channel 1 */
  RealTimeSeries** TDI2,                   /* Output: real time series for TDI channel 2 */
  RealTimeSeries** TDI3,                   /* Output: real time series for TDI channel 3 */
  gsl_spline* splinehp,                    /* Input spline for TD hplus */
  gsl_spline* splinehc,                    /* Input spline for TD hcross */
  gsl_interp_accel* accelhp,               /* Accelerator for hp spline */
  gsl_interp_accel* accelhc,               /* Accelerator for hc spline */
  gsl_vector* times,                       /* Vector of times to evaluate */
  int nbptsmargin,                         /* Margin set to 0 on both side to avoid problems with delays out of the domain */
  TDItag tditag);                          /* Tag selecting the TDI observables */

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LISAGEOMETRY_H */
