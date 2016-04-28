#ifndef __LLVUTILS_H__
#define __LLVUTILS_H__ 1

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
#include "splinecoeffs.h"
#include "LLVFDresponse.h"
#include "LLVnoise.h"

/***************** Structures for parameters *****************/

/* Parameters for the generation of a LLV waveform (in the form of a list of modes) */
typedef struct tagLLVParams {
  double tRef;               /* reference time (s) - GPS time at the frequency representing coalescence */
  double phiRef;             /* reference phase (rad) - phase at the frequency representing coalescence (or at fRef if specified) */
  double m1;                 /* mass of companion 1 (solar masses) */
  double m2;                 /* mass of companion 2 (solar masses) */
  double distance;           /* distance of source (Mpc) */
  double ra;                 /* right ascension of the source (rad) */
  double dec;                /* declination of the source (rad) */
  double inclination;        /* inclination of L relative to line of sight (rad) */
  double polarization;       /* polarization angle (rad) */
  int nbmode;                /* number of modes to generate (starting with 22) - defaults to 5 (all modes) */
} LLVParams;

/* Global parameters for the waveform generation and overlap computation */
typedef struct tagLLVGlobalParams {
  double fRef;               /* reference frequency (Hz, default 0 which is interpreted as Mf=0.14) */
  double minf;               /* Minimal frequency (Hz) - when set to 0 (default), use the first frequency covered by the noise data of the detector */
  double maxf;               /* Maximal frequency (Hz) - when set to 0 (default), use the last frequency covered by the noise data of the detector */
  int nbmodeinj;             /* number of modes to include in the injection (starting with 22) - defaults to 5 (all modes) */
  int nbmodetemp;            /* number of modes to include in the templates (starting with 22) - defaults to 5 (all modes) */
  int tagint;                /* Tag choosing the integrator: 0 for wip (default), 1 for linear integration */
  int tagnetwork;            /* Tag choosing the network of detectors to use */
  int nbptsoverlap;          /* Number of points to use in loglinear overlaps (default 32768) */
  int constL;                /* set all logLikelihood to 0 - allows to sample from the prior for testing */
} LLVGlobalParams;

typedef struct tagLLVPrior {
	double deltaT;             /* width of time prior centered on injected value (s) (default 0.1) */
	double comp_min;           /* minimum component mass (solar masses) (default 4) */
	double comp_max;           /* maximum component mass (solar masses) (default 50) */
	double mtot_min;           /* minimum total mass (solar masses) (default 8) */
	double mtot_max;           /* maximum total mass (solar masses) (default 100) */
	double qmax;               /* maximum asymmetric mass ratio (>=1) (default 12) */
	double dist_min;           /* minimum distance of source (Mpc) (default 1) */
	double dist_max;           /* maximum distance of source (Mpc) (default 1e4) */
  double ra_min;            /* minimum ra (rad, default 0) - for testing */
  double ra_max;            /* maximum ra (rad, default 2pi) - for testing */
  double dec_min;           /* minimum dec (rad, default 0) - for testing */
  double dec_max;           /* maximum dec (rad, default pi) - for testing */
  double phase_min;          /* minimum phase (rad, default 0) - for testing */
  double phase_max;          /* maximum phase (rad, default 2pi) - for testing */
  double pol_min;            /* minimum polarization (rad, default 0) - for testing */
  double pol_max;            /* maximum polarization (rad, default 2pi) - for testing */
  double inc_min;            /* minimum inclination (rad, default 0) - for testing */
  double inc_max;            /* maximum inclination (rad, default pi) - for testing */
  double fix_m1;
  double fix_m2;
  double fix_time;
  double fix_ra;
  double fix_dec;
  double fix_phase;
  double fix_pol;
  double fix_dist;
  double fix_inc;
  int pin_m1;
  int pin_m2;
  int pin_time;
  int pin_ra;
  int pin_dec;
  int pin_phase;
  int pin_pol;
  int pin_dist;
  int pin_inc;
  double snr_target;
  int rescale_distprior;
  int flat_distprior;
} LLVPrior;

typedef struct tagLLVRunParams {
	double eff;                /* target efficiency (default 0.1) */
	double tol;                /* logZ tolerance (default 0.5) */
	int    nlive;              /* number of live points (default 1000) */
	char   outroot[200];       /* output root (default "chains/LLVinference_") */
	int    bambi;              /* run BAMBI? (default 0) */
	int    resume;             /* resume form previous run? (default 0) */
	int    maxiter;            /* max number of iterations (default 0 - ignore) */
	char   netfile[200];       /* NN settings file (default "LLVinference.inp") */
  int    mmodal;             /* use multimodal decomposition ? */
  int    maxcls;             /* max number of modes in multimodal decomposition */
  int    nclspar;            /* number of parameters to use for multimodal decomposition - in the order of the cube */
  double ztol;               /* in multimodal decomposition, modes with lnZ lower than Ztol are ignored */
  int    seed;               /* seed the inference by setting one of the live points to the injection ? */
} LLVRunParams;

/************ Structures for signals and injections ************/

// typedef struct tagLLVSignal
// {
//   struct tagListmodesCAmpPhaseFrequencySeries* LHOSignal;   /* Signal in LHO, in the form of a list of the contribution of each mode */
//   struct tagListmodesCAmpPhaseFrequencySeries* LLOSignal;   /* Signal in LLO, in the form of a list of the contribution of each mode */
//   struct tagListmodesCAmpPhaseFrequencySeries* VIRGOSignal; /* Signal in VIRGO, in the form of a list of the contribution of each mode */
//   double LHOhh;                                             /* Inner product (h|h) for LHO */
//   double LLOhh;                                             /* Inner product (h|h) for LLO */
//   double VIRGOhh;                                           /* Inner product (h|h) for VIRGO */
// } LLVSignal;

typedef struct tagLLVSignalCAmpPhase
{
  struct tagListmodesCAmpPhaseFrequencySeries* LHOSignal;   /* Signal in LHO, in the form of a list of the contribution of each mode */
  struct tagListmodesCAmpPhaseFrequencySeries* LLOSignal;   /* Signal in LLO, in the form of a list of the contribution of each mode */
  struct tagListmodesCAmpPhaseFrequencySeries* VIRGOSignal; /* Signal in VIRGO, in the form of a list of the contribution of each mode */
  double LLVhh;                                             /* Combined Inner product (h|h) for dectectors LHV */
} LLVSignalCAmpPhase;

typedef struct tagLLVInjectionCAmpPhase
{
  struct tagListmodesCAmpPhaseSpline* LHOSplines;   /* Signal in LHO, in the form of a list of splines for the contribution of each mode */
  struct tagListmodesCAmpPhaseSpline* LLOSplines;   /* Signal in LLO, in the form of a list of splines for the contribution of each mode */
  struct tagListmodesCAmpPhaseSpline* VIRGOSplines;   /* Signal in VIRGO, in the form of a list of splines for the contribution of each mode */
  double LLVss;                                   /* Combined Inner product (s|s) for dectectors LHV */
} LLVInjectionCAmpPhase;

typedef struct tagLLVSignalReIm /* We don't store the SNRs here, as we will use -1/2(h-s|h-s) for the likelihood */
{
  struct tagReImFrequencySeries* LHOSignal;   /* Signal in LHO, in the form of a Re/Im frequency series where the modes have been summed */
  struct tagReImFrequencySeries* LLOSignal;   /* Signal in LLO, in the form of a Re/Im frequency series where the modes have been summed */
  struct tagReImFrequencySeries* VIRGOSignal;   /* Signal in VIRGO, in the form of a Re/Im frequency series where the modes have been summed */
} LLVSignalReIm;

typedef struct tagLLVInjectionReIm /* Storing the vectors of frequencies and noise values - We don't store the SNRs here, as we will use -1/2(h-s|h-s) for the likelihood */
{
  struct tagReImFrequencySeries* LHOSignal;   /* Signal in LHO, in the form of a Re/Im frequency series where the modes have been summed */
  struct tagReImFrequencySeries* LLOSignal;   /* Signal in LLO, in the form of a Re/Im frequency series where the modes have been summed */
  struct tagReImFrequencySeries* VIRGOSignal;   /* Signal in VIRGO, in the form of a Re/Im frequency series where the modes have been summed */
  gsl_vector* freq;                            /* Vector of frequencies of the injection (assumed to be the same for LHO, LLO, VIRGO) */
  gsl_vector* noisevaluesLHO;                    /* Vector of noise values on freq LHO */
  gsl_vector* noisevaluesLLO;                    /* Vector of noise values on freq LLO */
  gsl_vector* noisevaluesVIRGO;                  /* Vector of noise values on freq VIRGO */
} LLVInjectionReIm;

/************ Functions for LLV parameters, injection, likelihood, prior ************/

/* Parsing parameters for the generation of a LLV waveform, from the command line */
/* Masses are input in solar masses and distances in Mpc - converted in SI for the internals */
void parse_args_LLV(ssize_t argc, char **argv,
    LLVParams* params,
    LLVGlobalParams* globalparams,
    LLVPrior* prior,
    LLVRunParams* run,
    LLVParams* addparams);

/* Functions to print the parameters of the run in files for reference */
int print_parameters_to_file_LLV(
  LLVParams* params,
  LLVGlobalParams* globalparams,
  LLVPrior* prior,
  LLVRunParams* run);
int print_rescaleddist_to_file_LLV(
  LLVParams* params,
  LLVGlobalParams* globalparams,
  LLVPrior* prior,
  LLVRunParams* run);

/* Initialization and clean-up for LLVSignal structures */
void LLVSignalCAmpPhase_Cleanup(LLVSignalCAmpPhase* signal);
void LLVSignalCAmpPhase_Init(LLVSignalCAmpPhase** signal);
void LLVInjectionCAmpPhase_Cleanup(LLVInjectionCAmpPhase* signal);
void LLVInjectionCAmpPhase_Init(LLVInjectionCAmpPhase** signal);
void LLVSignalReIm_Cleanup(LLVSignalReIm* signal);
void LLVSignalReIm_Init(LLVSignalReIm** signal);
void LLVInjectionReIm_Cleanup(LLVInjectionReIm* signal);
void LLVInjectionReIm_Init(LLVInjectionReIm** signal);

/* Function generating a LLV signal as a list of modes in CAmp/Phase form, from LLV parameters */
int LLVGenerateSignalCAmpPhase(
  struct tagLLVParams* params,                 /* Input: set of LLV parameters of the signal */
  struct tagLLVSignalCAmpPhase* signal);  /* Output: structure for the generated signal */
/* Function generating a LLV injection as a list of modes, given as preinterpolated splines, from LLV parameters */
int LLVGenerateInjectionCAmpPhase(
  struct tagLLVParams* injectedparams,    /* Input: set of LLV parameters of the signal */
  struct tagLLVInjectionCAmpPhase* signal);  /* Output: structure for the generated signal */
/* Function generating a LLV signal as a frequency series in Re/Im form where the modes have been summed, from LLV parameters - takes as argument the frequencies on which to evaluate */
int LLVGenerateSignalReIm(
  struct tagLLVParams* params,       /* Input: set of LLV parameters of the template */
  gsl_vector* freq,                   /* Input: frequencies on which evaluating the waveform (from the injection) */
  struct tagLLVSignalReIm* signal);  /* Output: structure for the generated signal */
/* Function generating a LLV injection as a frequency series in Re/Im form where the modes have been summed, from LLV parameters - frequencies on which to evaluate are to be determined internally */
int LLVGenerateInjectionReIm(
  struct tagLLVParams* injectedparams,       /* Input: set of LLV parameters of the injection */
  double fLow,                                /* Input: starting frequency (from argument minf) */
  double fHigh,                               /* Input: upper frequency (from argument maxf) */
  int nbpts,                                  /* Input: number of frequency samples */
  int tagsampling,                            /* Input: tag for using linear (0) or logarithmic (1) sampling */
  struct tagLLVInjectionReIm* signal);       /* Output: structure for the generated signal */

// checks prior boundaires
int PriorBoundaryCheck(LLVPrior *prior, double *Cube);

// Prior functions from Cube to physical parameters
// x1 is min, x2 is max when specified
// r is Cube value
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
double CalculateLogLCAmpPhase(LLVParams *params, LLVInjectionCAmpPhase* injection);
double CalculateLogLReIm(LLVParams *params, LLVInjectionReIm* injection);

/************ Global Parameters ************/

extern LLVParams* injectedparams;
extern LLVGlobalParams* globalparams;
extern LLVPrior* priorParams;
double logZdata;

#endif
