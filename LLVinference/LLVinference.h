#ifndef __LLVINFERENCE_H__
#define __LLVINFERENCE_H__ 1

#ifdef PARALLEL
#include "mpi.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <stdbool.h>

#include "constants.h"
#include "struct.h"
#include "EOBNRv2HMROMstruct.h"
#include "EOBNRv2HMROM.h"
#include "wip.h"
#include "likelihood.h"
#include "LLVFDresponse.h"
#include "LLVnoise.h"
#include "LLVInit.h"
#include "bambi.h"

#ifdef __INTEL_COMPILER 			// if the MultiNest library was compiled with ifort
       #define NESTRUN nested_mp_nestrun_
#elif defined __GNUC__ 				// if the MultiNest library was compiled with gfortran
       #define NESTRUN __nested_MOD_nestrun
#else
       #error Do not know how to link to Fortran libraries, check symbol table for your platform (nm libnest3.a | grep nestrun) & edit example_eggbox_C++/eggbox.cc
#endif

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

/************ Functions for LLV parameters, injection, likelihood ************/

/* Parsing parameters for the generation of a LLV waveform, from the command line */
/* Masses are input in solar masses and distances in Mpc - converted in SI for the internals */
LLVParams* parse_args_LLV(ssize_t argc, char **argv);

/* Function generating a LLV signal from LLV parameters */
int LLVGenerateSignal(
  struct tagLLVParams* params,   /* Input: set of LLV parameters of the signal */
  struct tagLLVSignal* signal);  /* Output: structure for the generated signal */

/***************************************** BAMBI declarations **************************************************/

extern float *omicron,tol,thL[3],logLRange;
extern double *maxsigma,logZero;
extern int nNN,nNNerr,totpar,loglcalls,ncheck,myid,nproc;
extern char root[100],networkinputs[100];
extern bool likenetinit,converged,lastconverged,netres,firstrun,discardpts;
extern int ignoredbambicalls,counter;
extern size_t nlayers,nnodes[10];
extern int doBAMBI,useNN,whitenin,whitenout,resume;

/***************************************** C Interface to BAMBI **************************************************/

/*extern void NESTRUN(int *, int *, int *, double *, double *, int *, int *, int *, int *, int *, double *,
char *, int *, int *, int *, int *, int *, int *, double *, int *, void (*Loglike)(double *, int *, int *,
double *, void *), void (*dumper)(int *, int *, int *, double **, double **, double **, double *,
double *, double *, void *), void *context);*/


void BAMBIrun(int mmodal, int ceff, int nlive, double tol, double efr, int ndims, int nPar, int nClsPar,
	int maxModes, int updInt, double Ztol, char root[], int seed, int *pWrap, int fb, int resume, int outfile,
	int initMPI, double logZero, int maxiter, void (*LogLike)(double *, int *, int *, double *, void *),
	void (*dumper)(int *, int *, int *, double **, double **, double **, double *, double *, double *, void *),
	void (*bambi)(int *, int *, double **, double *), void *context)
{
	int i;
	char rootformn[100];
	strcpy(rootformn, root);
	for (i = strlen(rootformn); i < 100; i++) rootformn[i] = ' ';

        NESTRUN(&mmodal, &ceff, &nlive, &tol, &efr, &ndims, &nPar, &nClsPar, &maxModes, &updInt, &Ztol,
        rootformn, &seed, pWrap, &fb, &resume, &outfile, &initMPI, &logZero, &maxiter, LogLike, dumper, bambi, context);
}

void BAMBIfctn(int *ndata, int *ndim, double **BAMBIData, double *lowlike)
{
	// Do "nm bambi.o | grep bambi" to find the name of the function to put here.
	// "c++filt <fctn name>" should return "bambi(...)"
	// Remove one leading underscore for the name here.

	_Z5bambiPiS_PPdS0_(ndata, ndim, BAMBIData, lowlike);
}

void LogLikeFctn(double *Cube, int *ndim, int *npars, double *lnew, void *context)
{
	// Do "nm bambi.o | grep LogLike" to find the name of the function to put here.
	// "c++filt <fctn name>" should return "LogLike(...)"
	// Remove one leading underscore for the name here.

	_Z7LogLikePdPiS0_S_Pv(Cube, ndim, npars, lnew, context);
}


/***********************************************************************************************************************/

#endif
