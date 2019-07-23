#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <stdbool.h>

#include "LISAinference_common.h"


/************************************************** Main program *******************************************************/
/* This program works compatibly with LISAinference, but takes a list of 9 param values (in order as below) and only
   computes the likelihood at that parameter point.
*/
int noMPI=1;


int main(int argc, char *argv[])
{
  int myid = 0;
  noMPI = 1; //We have to set this to avoid an MPI_Comm_rank statement in the addendum

	LISARunParams runParams = {};
	int ndim=0, nPar=0;
	int *freeparamsmap = NULL;
	void *context = NULL;
	double logZtrue;
	addendum(argc, argv, &runParams, &ndim, &nPar, &freeparamsmap, &context, &logZtrue);

  LISAInjectionCAmpPhase* injectedsignalCAmpPhase = NULL;
  LISAInjectionReIm* injectedsignalReIm = NULL;
  double logL = 0;

  LISAParams* params = NULL;
  params = (LISAParams*) malloc(sizeof(LISAParams));
  memset(params, 0, sizeof(LISAParams));

  if(!(addparams->loadparamsfile)) {
    params->m1 = addparams->m1;
    params->m2 = addparams->m2;
    params->tRef = addparams->tRef;
    params->distance = addparams->distance;
    params->phiRef = addparams->phiRef;
    params->inclination = addparams->inclination;
    params->lambda = addparams->lambda;
    params->beta = addparams->beta;
    params->polarization = addparams->polarization;
    params->nbmode = globalparams->nbmodetemp; /* Note : read from global parameters */

    /* Report params */
    printf("Template params :\n");
    report_LISAParams(params);

    /* Calculate logL of template */
    if(globalparams->tagint==0) {
      injectedsignalCAmpPhase = (LISAInjectionCAmpPhase*) context;
      logL = CalculateLogLCAmpPhase(params, injectedsignalCAmpPhase);
    }
    else if(globalparams->tagint==1) {
      injectedsignalReIm = (LISAInjectionReIm*) context;
      logL = CalculateLogLReIm(params, injectedsignalReIm);
    }
    printf("logL template = %.16e\n", logL);
  }
  else {
    int nlines = addparams->nlinesparams;

    /* Load parameters file */
    /* Format (same as in the internals): m1, m2, tRef, dist, phase, inc, lambda, beta, pol, loglike, posteriormode */
    /* Assumes not using mmodal */
    gsl_matrix* inmatrix =  gsl_matrix_alloc(nlines, 10);
    Read_Text_Matrix(addparams->indir, addparams->infile, inmatrix);

    /* Initialize output matrix */
    /* Format (same as in the internals): m1, m2, tRef, dist, phase, inc, lambda, beta, pol, loglike */
    /* Assumes not using mmodal */
    gsl_matrix* outmatrix =  gsl_matrix_alloc(nlines, 10);

    if(globalparams->tagint==0) {
      injectedsignalCAmpPhase = ((LISAInjectionCAmpPhase*) context);
    }
    else if(globalparams->tagint==1) {
      injectedsignalReIm = ((LISAInjectionReIm*) context);
    }
    double logL = 0;
    for(int i=0; i<nlines; i++) {
      //
      if(i%100 == 0) printf("Nb computed: %d/%d\n", i, nlines);

      params->m1 = gsl_matrix_get(inmatrix, i, 0);
      params->m2 = gsl_matrix_get(inmatrix, i, 1);
      params->tRef = gsl_matrix_get(inmatrix, i, 2);
      params->distance = gsl_matrix_get(inmatrix, i, 3);
      params->phiRef = gsl_matrix_get(inmatrix, i, 4);
      params->inclination = gsl_matrix_get(inmatrix, i, 5);
      params->lambda = gsl_matrix_get(inmatrix, i, 6);
      params->beta = gsl_matrix_get(inmatrix, i, 7);
      params->polarization = gsl_matrix_get(inmatrix, i, 8);
      params->nbmode = globalparams->nbmodetemp; /* Note : read from global parameters */

      if(globalparams->tagint==0) {
        ///TEST
        printf("Before CalculateLogLCAmpPhase in loop\n");
        logL = CalculateLogLCAmpPhase(params, injectedsignalCAmpPhase);
      }
      else if(globalparams->tagint==1) {
        logL = CalculateLogLReIm(params, injectedsignalReIm);
      }

      /* Set values in output matrix */
      gsl_matrix_set(outmatrix, i, 0, params->m1);
      gsl_matrix_set(outmatrix, i, 1, params->m2);
      gsl_matrix_set(outmatrix, i, 2, params->tRef);
      gsl_matrix_set(outmatrix, i, 3, params->distance);
      gsl_matrix_set(outmatrix, i, 4, params->phiRef);
      gsl_matrix_set(outmatrix, i, 5, params->inclination);
      gsl_matrix_set(outmatrix, i, 6, params->lambda);
      gsl_matrix_set(outmatrix, i, 7, params->beta);
      gsl_matrix_set(outmatrix, i, 8, params->polarization);
      gsl_matrix_set(outmatrix, i, 9, logL);
    }
    /* Output matrix */
    Write_Text_Matrix(addparams->outdir, addparams->outfile, outmatrix);

    /* Cleanup */
    free(params);
  }

  /* Cleanup */
  free(injectedparams);
  free(globalparams);
  free(addparams);
  free(priorParams);
}
