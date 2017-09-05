#include "LISAinference_common.h"

/***********************************************************************************************************************/

/*********** Addendum *************/
void addendum(int argc, char *argv[],LISARunParams *runParams, int *ndim, int *nPar,int **freeparamsmapp,void **contextp, double *logZtrue){
  /* Initialize structs for holding various options */
  //LISARunParams runParams;
   int myid = 0;
#ifdef PARALLEL
   if(noMPI==0){
     printf("noMPI=%i\n",noMPI);
     //MPI_Init(&argc,&argv);
     MPI_Comm_rank(MPI_COMM_WORLD,&myid);
   }
#endif


  injectedparams = (LISAParams*) malloc(sizeof(LISAParams));
  memset(injectedparams, 0, sizeof(LISAParams));
  globalparams = (LISAGlobalParams*) malloc(sizeof(LISAGlobalParams));
  memset(globalparams, 0, sizeof(LISAGlobalParams));
  priorParams = (LISAPrior*) malloc(sizeof(LISAPrior));
  memset(priorParams, 0, sizeof(LISAPrior));
  addparams = (LISAAddParams*) malloc(sizeof(LISAAddParams));
  memset(addparams, 0, sizeof(LISAAddParams));
  /* This structure is used only for storing precomputed values for the simple likelihood */
  simplelikelihoodinjvals = (SimpleLikelihoodPrecomputedValues*) malloc(sizeof(SimpleLikelihoodPrecomputedValues));
  memset(simplelikelihoodinjvals, 0, sizeof(SimpleLikelihoodPrecomputedValues));

  /* Parse commandline to read parameters of injection - copy the number of modes demanded for the injection */
  parse_args_LISA(argc, argv, injectedparams, globalparams, priorParams, runParams, addparams);
  injectedparams->nbmode = globalparams->nbmodeinj;

  //int notLISAlike=strstr(argv[0],"LISAlike")==0;
  if(myid == 0 && runParams->writeparams /*&& notLISAlike*/) print_parameters_to_file_LISA(injectedparams, globalparams, priorParams, runParams);
  /* Initialize the data structure for the injection */
  LISAInjectionCAmpPhase* injectedsignalCAmpPhase = NULL;
  LISAInjectionReIm* injectedsignalReIm = NULL;
  if(globalparams->tagint==0) {
    LISAInjectionCAmpPhase_Init(&injectedsignalCAmpPhase);
  }
  else if(globalparams->tagint==1) {
    LISAInjectionReIm_Init(&injectedsignalReIm);
  }

  /* Generate the injection */
  if(globalparams->tagint==0) {
    LISAGenerateInjectionCAmpPhase(injectedparams, injectedsignalCAmpPhase);
  }
  else if(globalparams->tagint==1) {
    LISAGenerateInjectionReIm(injectedparams, globalparams->minf, globalparams->nbptsoverlap, 1, injectedsignalReIm); /* Use here logarithmic sampling as a default */
  }
  printf("Injected params\n");
  report_LISAParams(injectedparams);

  /* Define SNR */
  double SNR123, SNR1, SNR2, SNR3;
  if(globalparams->tagint==0) {
    SNR123 = sqrt(injectedsignalCAmpPhase->TDI123ss);
  }
  else if(globalparams->tagint==1) {
    SNR1 = sqrt(FDOverlapReImvsReIm(injectedsignalReIm->TDI1Signal, injectedsignalReIm->TDI1Signal, injectedsignalReIm->noisevalues1));
    SNR2 = sqrt(FDOverlapReImvsReIm(injectedsignalReIm->TDI2Signal, injectedsignalReIm->TDI2Signal, injectedsignalReIm->noisevalues2));
    SNR3 = sqrt(FDOverlapReImvsReIm(injectedsignalReIm->TDI3Signal, injectedsignalReIm->TDI3Signal, injectedsignalReIm->noisevalues3));
    SNR123 = sqrt(SNR1*SNR1 + SNR2*SNR2 + SNR3*SNR3);
  }

  /* Rescale distance to match SNR */
  //
  printf("isnan(priorParams->snr_target) : %d\n", isnan(priorParams->snr_target));
  if (!isnan(priorParams->snr_target)) {
    printf("SNR=%g\n",SNR123);
    if (myid == 0) printf("Rescaling the distance to obtain a network SNR of %g\n", priorParams->snr_target);
    injectedparams->distance *= SNR123 / priorParams->snr_target;
    if (myid == 0) printf("New distance = %g Mpc\n", injectedparams->distance);
    if(priorParams->rescale_distprior) {
      priorParams->dist_min *= SNR123 / priorParams->snr_target;
      priorParams->dist_max *= SNR123 / priorParams->snr_target;
      if (myid == 0) printf("Distance prior (dist_min, dist_max) = (%g, %g) Mpc\n", priorParams->dist_min, priorParams->dist_max);
    }
    if (myid == 0 && runParams->writeparams /*&& notLISAlike*/) print_rescaleddist_to_file_LISA(injectedparams, globalparams, priorParams, runParams);
    if(globalparams->tagint==0) {
      LISAGenerateInjectionCAmpPhase(injectedparams, injectedsignalCAmpPhase);
      SNR123 = sqrt(injectedsignalCAmpPhase->TDI123ss);
    }
    else if(globalparams->tagint==1) {
      LISAGenerateInjectionReIm(injectedparams, globalparams->minf, globalparams->nbptsoverlap, 1, injectedsignalReIm); /* tagsampling fixed to 1, i.e. logarithmic sampling - could be made another global parameter */
      SNR1 = sqrt(FDOverlapReImvsReIm(injectedsignalReIm->TDI1Signal, injectedsignalReIm->TDI1Signal, injectedsignalReIm->noisevalues1));
      SNR2 = sqrt(FDOverlapReImvsReIm(injectedsignalReIm->TDI2Signal, injectedsignalReIm->TDI2Signal, injectedsignalReIm->noisevalues2));
      SNR3 = sqrt(FDOverlapReImvsReIm(injectedsignalReIm->TDI3Signal, injectedsignalReIm->TDI3Signal, injectedsignalReIm->noisevalues3));
      SNR123 = sqrt(SNR1*SNR1 + SNR2*SNR2 + SNR3*SNR3);
    }
    printf("Rescaled injected params\n");
    report_LISAParams(injectedparams);
  }

  /* If using simple likelihood, initialize precomputed values - note that the other initializations for the injection are done anyway, but will be ignored */
  /* Note: the optional distance adjustment to a given snr is done above using the response as given by responseapprox, not the simplified response */
  if(globalparams->tagsimplelikelihood) {
    LISAComputeSimpleLikelihoodPrecomputedValues(simplelikelihoodinjvals, injectedparams);
  }

  /* Print SNR */
  if (myid == 0) {
    printf("Total SNR: %g\n", SNR123);
  }

  /* Calculate logL of data */
  /*double dist_store = injectedparams->distance;
    injectedparams->distance = 1.0e9;
    logZdata = CalculateLogL(injectedparams, injectedsignal);
    printf("logZdata = %lf\n", logZdata);
    injectedparams->distance = dist_store;*/
  logZdata = 0.0;
  *logZtrue = 0.;
  if(globalparams->tagint==0) {
    *logZtrue = CalculateLogLCAmpPhase(injectedparams, injectedsignalCAmpPhase);
  }
  else if(globalparams->tagint==1) {
    *logZtrue = CalculateLogLReIm(injectedparams, injectedsignalReIm);
  }
  /* printf("Compared params\n");
  report_LISAParams(injectedparams); */
  if(myid == 0) printf("logZtrue = %lf\n", *logZtrue-logZdata);

  /* Set the context pointer */
  if((globalparams->tagint==0) && (!globalparams->tagsimplelikelihood)) {
    *contextp = injectedsignalCAmpPhase;
  }
  else if((globalparams->tagint==1) && (!globalparams->tagsimplelikelihood)) {
    *contextp = injectedsignalReIm;
  }
  else if(globalparams->tagsimplelikelihood) {
    *contextp = simplelikelihoodinjvals;
  }

  *nPar = 9;	  /* Total no. of parameters including free & derived parameters */
  *ndim = 9;  /* No. of free parameters - to be changed later if some parameters are fixed */

  /* Check for parameters pinned to injected values */
  /* Distinguish the case where we sample in Mchirp/eta instead of m1/m2 */
  if(priorParams->samplemassparams==m1m2) {
    if(priorParams->pin_m1) priorParams->fix_m1 = injectedparams->m1;
    if(priorParams->pin_m2) priorParams->fix_m2 = injectedparams->m2;
  }
  if(priorParams->samplemassparams==Mchirpeta) {
    if(priorParams->pin_Mchirp) priorParams->fix_Mchirp = Mchirpofm1m2(injectedparams->m1, injectedparams->m2);
    if(priorParams->pin_eta) priorParams->fix_eta = etaofm1m2(injectedparams->m1, injectedparams->m2);
  }
  if(priorParams->pin_dist) priorParams->fix_dist = injectedparams->distance;
  if(priorParams->pin_inc) priorParams->fix_inc = injectedparams->inclination;
  if(priorParams->pin_phase) priorParams->fix_phase = injectedparams->phiRef;
  if(priorParams->pin_pol) priorParams->fix_pol = injectedparams->polarization;
  if(priorParams->pin_lambda) priorParams->fix_lambda = injectedparams->lambda;
  if(priorParams->pin_beta) priorParams->fix_beta = injectedparams->beta;
  if(priorParams->pin_time) priorParams->fix_time = injectedparams->tRef;

  /* Check for fixed parameters, and build the map from the free cube parameters to the orignal 9 parameters */
  /* Order of the 9 original parameters (fixed): m1, m2, tRef, dist, phase, inc, lambda, beta, pol */
  /* Order of the 9 cube parameters (modified for clustering): lambda, beta, tRef, phase, pol, inc, dist, m1, m2 */
  /* Distinguish the case where one is sampling in Mchirp/eta instead of m1/m2 */
  int mapcubetophys[9] = {6, 7, 2, 4, 8, 5, 3, 0, 1};
  int freecubeparams[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
  if(!isnan(priorParams->fix_lambda)) { (*ndim)--; freecubeparams[0] = 0; }
  if(!isnan(priorParams->fix_beta))   { (*ndim)--; freecubeparams[1] = 0; }
  if(!isnan(priorParams->fix_time))   { (*ndim)--; freecubeparams[2] = 0; }
  if(!isnan(priorParams->fix_phase))  { (*ndim)--; freecubeparams[3] = 0; }
  if(!isnan(priorParams->fix_pol))    { (*ndim)--; freecubeparams[4] = 0; }
  if(!isnan(priorParams->fix_inc))    { (*ndim)--; freecubeparams[5] = 0; }
  if(!isnan(priorParams->fix_dist))   { (*ndim)--; freecubeparams[6] = 0; }
  if(priorParams->samplemassparams==m1m2) {
    if (!isnan(priorParams->fix_m1))     { (*ndim)--; freecubeparams[7] = 0; }
    if (!isnan(priorParams->fix_m2))     { (*ndim)--; freecubeparams[8] = 0; }
  }
  if(priorParams->samplemassparams==Mchirpeta) {
    if (!isnan(priorParams->fix_Mchirp)) { (*ndim)--; freecubeparams[7] = 0; }
    if (!isnan(priorParams->fix_eta))    { (*ndim)--; freecubeparams[8] = 0; }
  }

  int *freeparamsmap = malloc(*ndim*sizeof(int));

  int counter = 0;
  for(int i=0; i<*ndim; i++) {
    while(freecubeparams[counter]==0) counter++;
    freeparamsmap[i] = mapcubetophys[counter];
    counter++;
  }

  if (*ndim == 0) {
    LISAParams templateparams;
    if(priorParams->samplemassparams==m1m2) {
      templateparams.m1 = priorParams->fix_m1;
      templateparams.m2 = priorParams->fix_m2;
    }
    if(priorParams->samplemassparams==Mchirpeta) {
      double Mchirp = priorParams->fix_Mchirp;
      double eta = priorParams->fix_eta;
      templateparams.m1 = m1ofMchirpeta(Mchirp, eta);
      templateparams.m2 = m2ofMchirpeta(Mchirp, eta);
    }
    templateparams.tRef = priorParams->fix_time;
    templateparams.distance = priorParams->fix_dist;
    templateparams.phiRef = priorParams->fix_phase;
    templateparams.inclination = priorParams->fix_inc;
    templateparams.lambda = priorParams->fix_lambda;
    templateparams.beta = priorParams->fix_beta;
    templateparams.polarization = priorParams->fix_pol;
    templateparams.nbmode = globalparams->nbmodetemp; /* Using the global parameter for the number of modes in templates */

    double logL = 0.;
    if(globalparams->tagint==0) {
      logL = CalculateLogLCAmpPhase(&templateparams, injectedsignalCAmpPhase);
    }
    else if(globalparams->tagint==1) {
      logL = CalculateLogLReIm(&templateparams, injectedsignalReIm);
    }
    printf("logL = %lf\n", logL);

    free(injectedparams);
    free(priorParams);

#ifdef PARALLEL
    MPI_Finalize();
#endif

    exit(0);
  }


  *freeparamsmapp = freeparamsmap;

  return;
}
/********** End of addendum ****************/
