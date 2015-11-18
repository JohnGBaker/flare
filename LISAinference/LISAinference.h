#ifndef __LISAINFERENCE_H__
#define __LISAINFERENCE_H__ 1

#ifdef PARALLEL
#include "mpi.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <stdbool.h>

#include "LISAutils.h"
#include "bambi.h"

#ifdef __INTEL_COMPILER 			/* if the MultiNest library was compiled with ifort */
       #define NESTRUN nested_mp_nestrun_
#elif defined __GNUC__ 				/* if the MultiNest library was compiled with gfortran */
       #define NESTRUN __nested_MOD_nestrun
#else
       #error Do not know how to link to Fortran libraries, check symbol table for your platform (nm libnest3.a | grep nestrun) & edit example_eggbox_C++/eggbox.cc
#endif

/***************************************** BAMBI declarations **************************************************/

extern float *omicron,tol,thL[3],logLRange;
extern double *maxsigma,logZero;
extern int nNN,nNNerr,totpar,loglcalls,ncheck,myid,nproc;
extern char root[1000],networkinputs[1000];
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
	char rootformn[1000];
	strcpy(rootformn, root);
	for (i = strlen(rootformn); i < 1000; i++) rootformn[i] = ' ';

        NESTRUN(&mmodal, &ceff, &nlive, &tol, &efr, &ndims, &nPar, &nClsPar, &maxModes, &updInt, &Ztol,
        rootformn, &seed, pWrap, &fb, &resume, &outfile, &initMPI, &logZero, &maxiter, LogLike, dumper, bambi, context);
}

void BAMBIfctn(int *ndata, int *ndim, double **BAMBIData, double *lowlike)
{
	 /* Do "nm bambi.o | grep bambi" to find the name of the function to put here. */
	 /* "c++filt <fctn name>" should return "bambi(...)" */
	 /* Remove one leading underscore for the name here. */

	_Z5bambiPiS_PPdS0_(ndata, ndim, BAMBIData, lowlike);
}

void LogLikeFctn(double *Cube, int *ndim, int *npars, double *lnew, void *context)
{
	 /* Do "nm bambi.o | grep LogLike" to find the name of the function to put here. */
	 /* "c++filt <fctn name>" should return "LogLike(...)" */
	 /* Remove one leading underscore for the name here. */

	_Z7LogLikePdPiS0_S_Pv(Cube, ndim, npars, lnew, context);
}


/***********************************************************************************************************************/

#endif
