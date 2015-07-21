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
  double fRef;               /* reference frequency (Hz, default 0 which is interpreted as Mf=0.14) */
  double deltatobs;          /* max duration of observation (years, default 2) - the start of the signals might be cut in time instead of cut in frequency */
  int nbmode;                /* number of modes to generate (starting with 22) - defaults to 5 (all modes) */
} LISAParams;

typedef struct tagLISASignal
{
  struct tagListmodesCAmpPhaseFrequencySeries* TDIASignal;   /* Signal in the 2nd generation TDI A, in the form of a list of the contribution of each mode */
  struct tagListmodesCAmpPhaseFrequencySeries* TDIESignal;   /* Signal in the 2nd generation TDI E, in the form of a list of the contribution of each mode */
  struct tagListmodesCAmpPhaseFrequencySeries* TDITSignal;   /* Signal in the 2nd generation TDI T, in the form of a list of the contribution of each mode */
  double TDIAhh;                                             /* Inner product (h|h) for TDI A */
  double TDIEhh;                                             /* Inner product (h|h) for TDI E */
  double TDIThh;                                             /* Inner product (h|h) for TDI T */
} LISASignal;

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

/* Parsing parameters for the generation of a LISA waveform, from the command line */
/* Masses are input in solar masses and distances in Mpc - converted in SI for the internals */
void parse_args_LISA(ssize_t argc, char **argv, 
    LISAParams* params, 
    LISAPrior *prior, 
    LISARunParams *run);

void LISASignal_Cleanup(LISASignal* signal);
void LISASignal_Init(LISASignal** signal);

/* Function generating a LISA signal from LISA parameters */
int LISAGenerateSignal(
  struct tagLISAParams* params,   /* Input: set of LISA parameters of the signal */
  struct tagLISASignal* signal);  /* Output: structure for the generated signal */

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

/* log-Likelihood function */
double CalculateLogL(LISAParams *params, LISASignal* injection);

/************ Global Parameters ************/

extern LISAParams* injectedparams;
extern LISAPrior* priorParams;
double logZdata;

#endif
