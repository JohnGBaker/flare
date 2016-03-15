#ifndef __LISAUTILS_H__
#define __LISAUTILS_H__ 1

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_cdf.h>

#include "constants.h"
#include "struct.h"
#include "EOBNRv2HMROMstruct.h"
#include "EOBNRv2HMROM.h"
#include "waveform.h"
#include "wip.h"
#include "splinecoeffs.h"
#include "fresnel.h"
#include "likelihood.h"
#include "LISAFDresponse.h"
#include "LISAnoise.h"

#if defined(__cplusplus)
extern "C" {
#define complex _Complex
#elif 0
} /* so that editors will match preceding brace */
#endif

/***************** Structure definitions *****************/

/* Parameters for the generation of a LISA waveform (in the form of a list of modes) */
typedef struct tagLISAParams {
  double tRef;               /* reference time (s) - GPS time at the frequency representing coalescence */
  double phiRef;             /* reference phase (rad) - phase at the frequency representing coalescence (or at fRef if specified) */
  double m1;                 /* mass of companion 1 (solar masses, default 2e6) */
  double m2;                 /* mass of companion 2 (solar masses, default 1e6) */
  double distance;           /* distance of source (Mpc, default 1e3) */
  double lambda;             /* first angle for the position in the sky (rad, default 0) */
  double beta;               /* second angle for the position in the sky (rad, default 0) */
  double inclination;        /* inclination of L relative to line of sight (rad, default PI/3) */
  double polarization;       /* polarization angle (rad, default 0) */
  int nbmode;                /* number of modes to generate (starting with 22) - defaults to 5 (all modes) */
} LISAParams;

/* Global parameters for the waveform generation and overlap computation */
typedef struct tagLISAGlobalParams {
  double fRef;               /* reference frequency (Hz, default 0 which is interpreted as Mf=0.14) */
  double deltatobs;          /* max duration of observation (years, default 2) - the start of the signals might be cut in time instead of cut in frequency */
  double minf;               /* Minimal frequency (Hz) - when set to 0 (default), use the first frequency covered by the noise data of the detector */
  int nbmodeinj;             /* number of modes to include in the injection (starting with 22) - defaults to 5 (all modes) */
  int nbmodetemp;            /* number of modes to include in the templates (starting with 22) - defaults to 5 (all modes) */
  int tagint;                /* Tag choosing the integrator: 0 for wip (default), 1 for linear integration */
  int tagtdi;                /* Tag choosing the TDI variables to use */
  int nbptsoverlap;          /* Number of points to use in loglinear overlaps (default 32768) */
} LISAGlobalParams;

typedef struct tagLISASignalCAmpPhase
{
  struct tagListmodesCAmpPhaseFrequencySeries* TDI1Signal;   /* Signal in the TDI channel 1, in the form of a list of the contribution of each mode */
  struct tagListmodesCAmpPhaseFrequencySeries* TDI2Signal;   /* Signal in the TDI channel 2, in the form of a list of the contribution of each mode */
  struct tagListmodesCAmpPhaseFrequencySeries* TDI3Signal;   /* Signal in the TDI channel 3, in the form of a list of the contribution of each mode */
  double TDI123hh;                                           /* Combined Inner product (h|h) for TDI channels 123 */
} LISASignalCAmpPhase;

typedef struct tagLISAInjectionCAmpPhase
{
  struct tagListmodesCAmpPhaseSpline* TDI1Splines;   /* Signal in the TDI channel 1, in the form of a list of splines for the contribution of each mode */
  struct tagListmodesCAmpPhaseSpline* TDI2Splines;   /* Signal in the TDI channel 2, in the form of a list of splines for the contribution of each mode */
  struct tagListmodesCAmpPhaseSpline* TDI3Splines;   /* Signal in the TDI channel 3, in the form of a list of splines for the contribution of each mode */
  double TDI123ss;                                   /* Combined Inner product (s|s) for TDI channels 123 */
} LISAInjectionCAmpPhase;

typedef struct tagLISASignalReIm /* We don't store the SNRs here, as we will use -1/2(h-s|h-s) for the likelihood */
{
  struct tagReImFrequencySeries* TDI1Signal;   /* Signal in the TDI channel 1, in the form of a Re/Im frequency series where the modes have been summed */
  struct tagReImFrequencySeries* TDI2Signal;   /* Signal in the TDI channel 2, in the form of a Re/Im frequency series where the modes have been summed */
  struct tagReImFrequencySeries* TDI3Signal;   /* Signal in the TDI channel 3, in the form of a Re/Im frequency series where the modes have been summed */
} LISASignalReIm;

typedef struct tagLISAInjectionReIm /* Storing the vectors of frequencies and noise values - We don't store the SNRs here, as we will use -1/2(h-s|h-s) for the likelihood */
{
  struct tagReImFrequencySeries* TDI1Signal;   /* Signal in the TDI channel 1, in the form of a Re/Im frequency series where the modes have been summed */
  struct tagReImFrequencySeries* TDI2Signal;   /* Signal in the TDI channel 2, in the form of a Re/Im frequency series where the modes have been summed */
  struct tagReImFrequencySeries* TDI3Signal;   /* Signal in the TDI channel 3, in the form of a Re/Im frequency series where the modes have been summed */
  gsl_vector* freq;                            /* Vector of frequencies of the injection (assumed to be the same for A,E,T) */
  gsl_vector* noisevalues1;                    /* Vector of noise values on freq for TDI channel 1 */
  gsl_vector* noisevalues2;                    /* Vector of noise values on freq for TDI channel 2 */
  gsl_vector* noisevalues3;                    /* Vector of noise values on freq for TDI channel 3 */
} LISAInjectionReIm;

typedef struct tagLISAPrior {
  double deltaT;             /* width of time prior centered on injected value (s) (default 1e5) */
  double comp_min;           /* minimum component mass (solar masses) (default 1e4) */
  double comp_max;           /* maximum component mass (solar masses) (default 1e8) */
  double mtot_min;           /* minimum total mass (solar masses) (default 5*1e4) */
  double mtot_max;           /* maximum total mass (solar masses) (default 1e8) */
  double qmax;               /* maximum asymmetric mass ratio (>=1) (default 11.98) */
  double dist_min;           /* minimum distance of source (Mpc) (default 100) */
  double dist_max;           /* maximum distance of source (Mpc) (default 40*1e3) */
  double lambda_min;         /* minimum lambda (rad, default 0) - for testing */
  double lambda_max;         /* maximum lambda (rad, default 2pi) - for testing */
  double beta_min;           /* minimum beta (rad, default 0) - for testing */
  double beta_max;           /* maximum beta (rad, default pi) - for testing */
  double phase_min;          /* minimum phase (rad, default 0) - for testing */
  double phase_max;          /* maximum phase (rad, default 2pi) - for testing */
  double pol_min;            /* minimum polarization (rad, default 0) - for testing */
  double pol_max;            /* maximum polarization (rad, default 2pi) - for testing */
  double inc_min;            /* minimum inclination (rad, default 0) - for testing */
  double inc_max;            /* maximum inclination (rad, default pi) - for testing */
  double fix_m1;
  double fix_m2;
  double fix_time;
  double fix_lambda;
  double fix_beta;
  double fix_phase;
  double fix_pol;
  double fix_dist;
  double fix_inc;
  int pin_m1;
  int pin_m2;
  int pin_time;
  int pin_lambda;
  int pin_beta;
  int pin_phase;
  int pin_pol;
  int pin_dist;
  int pin_inc;
  double snr_target;
  int rescale_distprior;
  int flat_distprior;
  int logflat_massprior;
} LISAPrior;

typedef struct tagLISARunParams {
  double eff;                /* target efficiency (default 0.1) */
  double tol;                /* logZ tolerance (default 0.5) */
  int    nlive;              /* number of live points (default 1000) */
  char   outroot[200];       /* output root (default "chains/LISAinference_") */
  int    bambi;              /* run BAMBI? (default 0) */
  int    resume;             /* resume form previous run? (default 0) */
  int    maxiter;            /* max number of iterations (default 0 - ignore) */
  char   netfile[200];       /* NN settings file (default "LISAinference.inp") */
  int    mmodal;             /* use multimodal decomposition ? */
  int    maxcls;             /* max number of modes in multimodal decomposition */
  int    nclspar;            /* number of parameters to use for multimodal decomposition - in the order of the cube */
  double ztol;               /* in multimodal decomposition, modes with lnZ lower than Ztol are ignored */
  int    seed;               /* seed the inference by setting one of the live points to the injection ? */
} LISARunParams;

/************ Functions for LISA parameters, injection, likelihood, prior ************/

/* Parse command line to initialize LISAParams, LISAGlobalParams, LISAPrior, and LISARunParams objects */
/* Masses are input in solar masses and distances in Mpc - converted in SI for the internals */
void parse_args_LISA(ssize_t argc, char **argv,
  LISAParams* params,
  LISAGlobalParams* globalparams,
  LISAPrior* prior,
  LISARunParams* run);

/* Functions to print the parameters of the run in files for reference */
int print_parameters_to_file_LISA(
  LISAParams* params,
  LISAGlobalParams* globalparams,
  LISAPrior* prior,
  LISARunParams* run);
int print_rescaleddist_to_file_LISA(
  LISAParams* params,
  LISAGlobalParams* globalparams,
  LISAPrior* prior,
  LISARunParams* run);

/* Initialization and clean-up for LISASignal structures */
void LISASignalCAmpPhase_Cleanup(LISASignalCAmpPhase* signal);
void LISASignalCAmpPhase_Init(LISASignalCAmpPhase** signal);
void LISAInjectionCAmpPhase_Cleanup(LISAInjectionCAmpPhase* signal);
void LISAInjectionCAmpPhase_Init(LISAInjectionCAmpPhase** signal);
void LISASignalReIm_Cleanup(LISASignalReIm* signal);
void LISASignalReIm_Init(LISASignalReIm** signal);
void LISAInjectionReIm_Cleanup(LISAInjectionReIm* signal);
void LISAInjectionReIm_Init(LISAInjectionReIm** signal);

/* Function generating a LISA signal as a list of modes in CAmp/Phase form, from LISA parameters */
int LISAGenerateSignalCAmpPhase(
  struct tagLISAParams* params,                 /* Input: set of LISA parameters of the signal */
  struct tagLISASignalCAmpPhase* signal);  /* Output: structure for the generated signal */
/* Function generating a LISA injection as a list of modes, given as preinterpolated splines, from LISA parameters */
int LISAGenerateInjectionCAmpPhase(
  struct tagLISAParams* injectedparams,    /* Input: set of LISA parameters of the signal */
  struct tagLISAInjectionCAmpPhase* signal);  /* Output: structure for the generated signal */
/* Function generating a LISA signal as a frequency series in Re/Im form where the modes have been summed, from LISA parameters - takes as argument the frequencies on which to evaluate */
int LISAGenerateSignalReIm(
  struct tagLISAParams* params,       /* Input: set of LISA parameters of the template */
  gsl_vector* freq,                   /* Input: frequencies on which evaluating the waveform (from the injection) */
  struct tagLISASignalReIm* signal);  /* Output: structure for the generated signal */
/* Function generating a LISA injection as a frequency series in Re/Im form where the modes have been summed, from LISA parameters - frequencies on which to evaluate are to be determined internally */
int LISAGenerateInjectionReIm(
  struct tagLISAParams* injectedparams,       /* Input: set of LISA parameters of the injection */
  double fLow,                                /* Input: starting frequency */
  int nbpts,                                  /* Input: number of frequency samples */
  int tagsampling,                            /* Input: tag for using linear (0) or logarithmic (1) sampling */
  struct tagLISAInjectionReIm* signal);       /* Output: structure for the generated signal */

/* checks prior boundaires */
int PriorBoundaryCheck(LISAPrior *prior, double *Cube);

/* Prior functions from Cube to physical parameters
   x1 is min, x2 is max when specified
   r is Cube value */
double CubeToFlatPrior(double r, double x1, double x2);
double CubeToLogFlatPrior(double r, double x1, double x2);
double CubeToPowerPrior(double p, double r, double x1, double x2);
double CubeToGaussianPrior(double r, double mean, double sigma);
double CubeToSinPrior(double r, double x1, double x2);
double CubeToCosPrior(double r, double x1, double x2);
/* Prior functions from physical parameters to Cube
   x1 is min, x2 is max when specified
   y is physical value */
double FlatPriorToCube(double y, double x1, double x2);
double LogFlatPriorToCube(double y, double x1, double x2);
double PowerPriorToCube(double p, double y, double x1, double x2);
double SinPriorToCube(double y, double x1, double x2);
double CosPriorToCube(double y, double x1, double x2);

/* log-Likelihood functions */
double CalculateLogLCAmpPhase(LISAParams *params, LISAInjectionCAmpPhase* injection);
double CalculateLogLReIm(LISAParams *params, LISAInjectionReIm* injection);

/************ Global Parameters ************/

extern LISAParams* injectedparams;
extern LISAGlobalParams* globalparams;
extern LISAPrior* priorParams;
double logZdata;

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif
