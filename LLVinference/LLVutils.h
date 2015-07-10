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
#include "LLVFDresponse.h"
#include "LLVnoise.h"

/***************** Structure definitions *****************/

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
  double fRef;               /* reference frequency (Hz) */
  int nbmode;                /* number of modes to generate (starting with 22) - defaults to 5 (all modes) */
} LLVParams;

typedef struct tagLLVSignal
{
  struct tagListmodesCAmpPhaseFrequencySeries* LHOSignal;   /* Signal in LHO, in the form of a list of the contribution of each mode */
  struct tagListmodesCAmpPhaseFrequencySeries* LLOSignal;   /* Signal in LLO, in the form of a list of the contribution of each mode */
  struct tagListmodesCAmpPhaseFrequencySeries* VIRGOSignal; /* Signal in VIRGO, in the form of a list of the contribution of each mode */
  double LHOhh;                                             /* Inner product (h|h) for LHO */
  double LLOhh;                                             /* Inner product (h|h) for LLO */
  double VIRGOhh;                                           /* Inner product (h|h) for VIRGO */
} LLVSignal;

typedef struct tagLLVPrior {
	double deltaT;             /* width of time prior centered on injected value (s) (default 0.1) */
	double comp_min;           /* minimum component mass (solar masses) (default 4) */
	double comp_max;           /* maximum component mass (solar masses) (default 50) */
	double mtot_min;           /* minimum total mass (solar masses) (default 8) */
	double mtot_max;           /* maximum total mass (solar masses) (default 100) */
	double qmax;               /* maximum asymmetric mass ratio (>=1) (default 12) */
	double dist_min;           /* minimum distance of source (pc) (default 1e6) */
	double dist_max;           /* maximum distance of source (pc) (default 10*1e9) */
} LLVPrior;

typedef struct tagLLVRunParams {
	double eff;                /* target efficiency (default 0.1) */
	double tol;                /* logZ tolerance (default 0.5) */
	int    nlive;              /* number of live points (default 1000) */
	char   outroot[200];       /* output root (default "chains/LLVinference_") */
	int    bambi;              /* run BAMBI? (default 0) */
	int    resume;             /* resume form previous run? (default 0) */
	char   netfile[200];       /* NN settings file (default "LLVinference.inp") */
} LLVRunParams;

/************ Functions for LLV parameters, injection, likelihood ************/

/* Parsing parameters for the generation of a LLV waveform, from the command line */
/* Masses are input in solar masses and distances in Mpc - converted in SI for the internals */
LLVParams* parse_args_LLV(ssize_t argc, char **argv);

void LLVSignal_Cleanup(LLVSignal* signal);
void LLVSignal_Init(LLVSignal** signal);

/* Function generating a LLV signal from LLV parameters */
int LLVGenerateSignal(
  struct tagLLVParams* params,   /* Input: set of LLV parameters of the signal */
  struct tagLLVSignal* signal);  /* Output: structure for the generated signal */

// initializes the prior boundaries
LLVPrior* LLVInitializePrior(ssize_t argc, char **argv);

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

// initialize BAMBI run parameters
LLVRunParams* LLVInitializeRunParams(ssize_t argc, char **argv);

/************ Global Parameters ************/

extern LLVParams* injectedparams;
extern LLVParams* templateparams;
extern LLVPrior* priorParams;

#endif
