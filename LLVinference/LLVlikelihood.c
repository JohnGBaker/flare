#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <stdbool.h>

#include "LLVutils.h"


/************************************************** Main program *******************************************************/
/* This program works compatibly with LISAinference, but takes a list of 9 param values (in order as below) and only
   computes the likelihood at that parameter point.
*/


int main(int argc, char *argv[])
{
  /* Initialize structs for holding various options */
  LLVRunParams runParams;
  injectedparams = (LLVParams*) malloc(sizeof(LLVParams));
  memset(injectedparams, 0, sizeof(LLVParams));
  globalparams = (LLVGlobalParams*) malloc(sizeof(LLVGlobalParams));
  memset(globalparams, 0, sizeof(LLVGlobalParams));
  priorParams = (LLVPrior*) malloc(sizeof(LLVPrior));
  memset(priorParams, 0, sizeof(LLVPrior));
	LLVParams* addparams = (LLVParams*) malloc(sizeof(LLVParams));
  memset(addparams, 0, sizeof(LLVParams));

  /* Parse commandline to read parameters of injection - copy the number of modes demanded for the injection  */
  parse_args_LLV(argc, argv, injectedparams, globalparams, priorParams, &runParams, addparams);
  injectedparams->nbmode = globalparams->nbmodeinj;
  addparams->nbmode = globalparams->nbmodetemp;

  /* Load and initialize the detector noise */
  LLVSimFD_Noise_Init_ParsePath();

  /* Initialize the data structure for the injection */
  LLVInjectionCAmpPhase* injectedsignalCAmpPhase = NULL;
  LLVInjectionReIm* injectedsignalReIm = NULL;
  if(globalparams->tagint==0) {
    LLVInjectionCAmpPhase_Init(&injectedsignalCAmpPhase);
  }
  else if(globalparams->tagint==1) {
    LLVInjectionReIm_Init(&injectedsignalReIm);
  }

  /* Generate the injection */
  if(globalparams->tagint==0) {
    LLVGenerateInjectionCAmpPhase(injectedparams, injectedsignalCAmpPhase);
  }
  else if(globalparams->tagint==1) {
    LLVGenerateInjectionReIm(injectedparams, globalparams->minf, globalparams->maxf, globalparams->nbptsoverlap, 1, injectedsignalReIm); /* Use here logarithmic sampling as a default */
  }

  /* Print parameters */
  printf("--------------------------------------------------------\n");
  printf("Params |       Injection        |       Template        \n");
  printf("--------------------------------------------------------\n");
  printf("m1     | %.16e | %.16e\n", injectedparams->m1, addparams->m1);
  printf("m2     | %.16e | %.16e\n", injectedparams->m2, addparams->m2);
  printf("tRef   | %.16e | %.16e\n", injectedparams->tRef, addparams->tRef);
  printf("phiRef | %.16e | %.16e\n", injectedparams->phiRef, addparams->phiRef);
  printf("dist   | %.16e | %.16e\n", injectedparams->distance, addparams->distance);
  printf("ra     | %.16e | %.16e\n", injectedparams->ra, addparams->ra);
  printf("dec    | %.16e | %.16e\n", injectedparams->dec, addparams->dec);
  printf("inc    | %.16e | %.16e\n", injectedparams->inclination, addparams->inclination);
  printf("pol    | %.16e | %.16e\n", injectedparams->polarization, addparams->polarization);
  printf("nbmode |                      %d |                      %d\n", injectedparams->nbmode, addparams->nbmode);
  printf("--------------------------------------------------------\n");

  /* Calculate logL of injection */
  double logZinj = 0;
  if(globalparams->tagint==0) {
    logZinj = CalculateLogLCAmpPhase(injectedparams, injectedsignalCAmpPhase);
  }
  else if(globalparams->tagint==1) {
    logZinj = CalculateLogLReIm(injectedparams, injectedsignalReIm);
  }
  printf("logZinj  = %.16e\n", logZinj);

  /* Calculate logL of template */
  double logZtemp = 0;
  if(globalparams->tagint==0) {
    logZtemp = CalculateLogLCAmpPhase(addparams, injectedsignalCAmpPhase);
  }
  else if(globalparams->tagint==1) {
    logZtemp = CalculateLogLReIm(addparams, injectedsignalReIm);
  }
  printf("logZtemp = %.16e\n", logZtemp);

}
