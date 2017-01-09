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
  noMPI = 1; //We have to set this to avoid an MPI_Comm_rank statement in the addendum
  int i = 0;

  LISARunParams runParams={};
  int ndim=9,nPar=9;
  int *freeparamsmap = NULL;
  void *context = NULL;
  double logZtrue;
  LISAParams params;
  LISAPrior prior;
  double paramvals[nPar];
  double result;

  for (i = argc-nPar; i < argc; ++i) {
    if (i<1||strncmp(argv[i], "--xxx",2) == 0) {
      if(i<1) printf("Expected %i params at end of argument list.\n", nPar);
      else printf("Expected %i params at end of argument list but found '%s' in place of param %i.\n", nPar,argv[i],i+nPar-argc+1 );
      char* fakeargv[2]={argv[0],"--help"};
      parse_args_LISA(2, fakeargv, &params,globalparams,&prior,&runParams);
      exit(0);
    }
    double val = atof(argv[i]);
    printf("read val: %20.15g\n",val);
    paramvals[i+nPar-argc] =  atof(argv[i]);val;
  }
  addendum(argc-nPar,argv,&runParams,&ndim,&nPar,&freeparamsmap,&context,&logZtrue);

  printf("nbmode_inj   =%i\n",globalparams->nbmodeinj);
  printf("nbmode_templ =%i\n",globalparams->nbmodetemp);
  params.m1 = paramvals[0];
  params.m2 = paramvals[1];
  params.tRef = paramvals[2];
  //
  //printf("params.tRef: %g\n", params.tRef);
  params.distance = paramvals[3];
  if(!priorParams->flat_distprior)//If not using flat prior on distance then we need to transform from s(D) to D.
    params.distance= pow(paramvals[3] * pow(priorParams->dist_max, 3) + (1.0 - paramvals[3]) * pow(priorParams->dist_min, 3), 1.0 / 3.0);
  params.phiRef = paramvals[4];
  params.inclination = paramvals[5];
  params.lambda = paramvals[6];
  params.beta = paramvals[7];
  params.polarization = paramvals[8];
  params.nbmode = globalparams->nbmodetemp; /* Using the global parameter for the number of modes in templates */
  printf("Considering parameters: \n   m1 = %g\n   m2 = %g\n   t0 = %g\n",params.m1,params.m2,params.tRef);
  if(!priorParams->flat_distprior)printf(" s[d] = %g -->",paramvals[3]);
  printf(" dist = %g\n phi0 = %g\n incl = %g\n  lam = %g\n beta = %g\n  pol = %g\n",
	 params.distance,params.phiRef,params.inclination,params.lambda,params.beta,params.polarization);
  if(globalparams->tagint==0) {
    LISAInjectionCAmpPhase* injection = ((LISAInjectionCAmpPhase*) context);
    printf("Comparing with injection parameters: \n   m1 = %g\n   m2 = %g\n   t0 = %g\n dist = %g\n phi0 = %g\n incl = %g\n  lam = %g\n beta = %g\n  pol = %g\n",
	   injectedparams->m1,injectedparams->m2,injectedparams->tRef,injectedparams->distance,injectedparams->phiRef,injectedparams->inclination,injectedparams->lambda,injectedparams->beta,injectedparams->polarization);
    result = CalculateLogLCAmpPhase(&params, injection) - logZdata;
    printf("likelihood:CalculateLogL=%g.\n",result);
  }
  else if(globalparams->tagint==1) {
    LISAInjectionReIm* injection = ((LISAInjectionReIm*) context);
    result = CalculateLogLReIm(&params, injection) - logZdata;
  }

  for(int i=0;i<nPar-1;i++){
    //cout<<"i="<<i<<endl;
    printf("%20.15g,",paramvals[i]);
  }
  printf("%20.15g\n",paramvals[nPar-1]);
  printf("The result: likelihood = %23.13g\n",result);
  free(injectedparams);
  //free(priorParams);

}
