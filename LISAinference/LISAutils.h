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
#include "wip.h"
#include "likelihood.h"
#include "LISAFDresponse.h"
#include "LISAnoise.h"

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
  double fmin;               /* Minimal frequency (Hz) - when set to 0 (default), use the first frequency covered by the noise data of the detector */
  int nbmodeinj;             /* number of modes to include in the injection (starting with 22) - defaults to 5 (all modes) */
  int nbmodetemp;            /* number of modes to include in the templates (starting with 22) - defaults to 5 (all modes) */
  int tagint;                /* Tag choosing the integrator: 0 for wip (default), 1 for linear integration */
  int nbptsoverlap;          /* Number of points to use in loglinear overlaps (default 32768) */
} LISAGlobalParams;

typedef struct tagLISASignalCAmpPhase
{
  struct tagListmodesCAmpPhaseFrequencySeries* TDIASignal;   /* Signal in the 2nd generation TDI A, in the form of a list of the contribution of each mode */
  struct tagListmodesCAmpPhaseFrequencySeries* TDIESignal;   /* Signal in the 2nd generation TDI E, in the form of a list of the contribution of each mode */
  struct tagListmodesCAmpPhaseFrequencySeries* TDITSignal;   /* Signal in the 2nd generation TDI T, in the form of a list of the contribution of each mode */
  double TDIAhh;                                             /* Inner product (h|h) for TDI A */
  double TDIEhh;                                             /* Inner product (h|h) for TDI E */
  double TDIThh;                                             /* Inner product (h|h) for TDI T */
} LISASignalCAmpPhase;

typedef struct tagLISASignalReIm /* We don't store the SNRs here, as we will use -1/2(h-s|h-s) for the likelihood */
{
  struct tagReImFrequencySeries* TDIASignal;   /* Signal in the 2nd generation TDI A, in the form of a Re/Im frequency series where the modes have been summed */
  struct tagReImFrequencySeries* TDIESignal;   /* Signal in the 2nd generation TDI E, in the form of a Re/Im frequency series where the modes have been summed */
  struct tagReImFrequencySeries* TDITSignal;   /* Signal in the 2nd generation TDI T, in the form of a Re/Im frequency series where the modes have been summed */
} LISASignalReIm;

typedef struct tagLISAInjectionReIm /* Storing the vectors of frequencies and noise values - We don't store the SNRs here, as we will use -1/2(h-s|h-s) for the likelihood */
{
  struct tagReImFrequencySeries* TDIASignal;   /* Signal in the 2nd generation TDI A, in the form of a Re/Im frequency series where the modes have been summed */
  struct tagReImFrequencySeries* TDIESignal;   /* Signal in the 2nd generation TDI E, in the form of a Re/Im frequency series where the modes have been summed */
  struct tagReImFrequencySeries* TDITSignal;   /* Signal in the 2nd generation TDI T, in the form of a Re/Im frequency series where the modes have been summed */
  gsl_vector* freq;                            /* Vector of frequencies of the injection (assumed to be the same for A,E,T) */
  gsl_vector* noisevaluesA;                    /* Vector of noise values on freq for A */
  gsl_vector* noisevaluesE;                    /* Vector of noise values on freq for E */
  gsl_vector* noisevaluesT;                    /* Vector of noise values on freq for T */
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
  double fix_m1;
  double fix_m2;
  double fix_time;
  double fix_lambda;
  double fix_beta;
  double fix_phase;
  double fix_pol;
  double fix_dist;
  double fix_inc;
  double snr_target;
} LISAPrior;

typedef struct tagLISARunParams {
	double eff;                /* target efficiency (default 0.1) */
	double tol;                /* logZ tolerance (default 0.5) */
	int    nlive;              /* number of live points (default 1000) */
	char   outroot[200];       /* output root (default "chains/LISAinference_") */
	int    bambi;              /* run BAMBI? (default 0) */
	int    resume;             /* resume form previous run? (default 0) */
	char   netfile[200];       /* NN settings file (default "LISAinference.inp") */
} LISARunParams;

/************ Functions for LISA parameters, injection, likelihood, prior ************/

/* Parse command line to initialize LISAParams, LISAGlobalParams, LISAPrior, and LISARunParams objects */
/* Masses are input in solar masses and distances in Mpc - converted in SI for the internals */
void parse_args_LISA(ssize_t argc, char **argv, 
  LISAParams* params,
  LISAGlobalParams* globalparams, 
  LISAPrior* prior, 
  LISARunParams* run);

/* Initialization and clean-up for LISASignal structures */
void LISASignalCAmpPhase_Cleanup(LISASignalCAmpPhase* signal);
void LISASignalCAmpPhase_Init(LISASignalCAmpPhase** signal);
void LISASignalReIm_Cleanup(LISASignalReIm* signal);
void LISASignalReIm_Init(LISASignalReIm** signal);
void LISAInjectionReIm_Cleanup(LISAInjectionReIm* signal);
void LISAInjectionReIm_Init(LISAInjectionReIm** signal);

/* Function generating a LISA signal as a list of modes in CAmp/Phase form, from LISA parameters */
int LISAGenerateSignalCAmpPhase(
  struct tagLISAParams* params,            /* Input: set of LISA parameters of the signal */
  struct tagLISASignalCAmpPhase* signal);  /* Output: structure for the generated signal */
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
  struct tagLISAInjectionReIm* signal);          /* Output: structure for the generated signal */

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

/* log-Likelihood functions */
double CalculateLogLCAmpPhase(LISAParams *params, LISASignalCAmpPhase* injection);
double CalculateLogLReIm(LISAParams *params, LISAInjectionReIm* injection);


/************ Global Parameters ************/

extern LISAParams* injectedparams;
extern LISAGlobalParams* globalparams;
extern LISAPrior* priorParams;
double logZdata;
extern gsl_vector* noisevaluesA;
extern gsl_vector* noisevaluesE;
extern gsl_vector* noisevaluesT;

#endif
