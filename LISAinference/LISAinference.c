#include "LISAinference.h"

/******************************************** getphysparams routine ****************************************************/

void getphysparams(double *Cube, int *ndim)
{
	int i = 0;
  double m1=0., m2=0., tRef=0., dist=0., phase=0., inc=0., lambda=0., beta=0., pol=0.;

  // component masses
  if (isnan(priorParams->fix_m1)) {
    if (isnan(priorParams->fix_m2)) {
      m1 = CubeToFlatPrior(Cube[i++], priorParams->comp_min, priorParams->comp_max);
      m2 = CubeToFlatPrior(Cube[i++], priorParams->comp_min, priorParams->comp_max);
      if (m2 > m1) {
        double tmp = m1;
        m1 = m2;
        m2 = tmp;
      }
    } else {
      m2 = priorParams->fix_m2;
      m1 = CubeToFlatPrior(Cube[i++], fmax(priorParams->comp_min,m2), priorParams->comp_max);
    }
  } else {
    m1 = priorParams->fix_m1;
    if (isnan(priorParams->fix_m2)) {
      m2 = CubeToFlatPrior(Cube[i++], priorParams->comp_min, fmin(m1,priorParams->comp_max));
    } else {
      m2 = priorParams->fix_m2;
    }
  }

  // time
  if (isnan(priorParams->fix_time)) {
    tRef = CubeToFlatPrior(Cube[i++], injectedparams->tRef - priorParams->deltaT,
           injectedparams->tRef + priorParams->deltaT);
  } else {
    tRef = priorParams->fix_time;
  }

  // distance
  if (isnan(priorParams->fix_dist)) {
    dist = CubeToPowerPrior(2.0, Cube[i++], priorParams->dist_min, priorParams->dist_max);
  } else {
    dist = priorParams->fix_dist;
  }

  // orbital phase
  if (isnan(priorParams->fix_phase)) {
    phase = CubeToFlatPrior(Cube[i++], 0.0, 2.0 * M_PI);
  } else {
    phase = priorParams->fix_phase;
  }

  // inclination
  if (isnan(priorParams->fix_inc)) {
    inc = CubeToSinPrior(Cube[i++], 0.0, M_PI);
  } else {
    inc = priorParams->fix_inc;
  }

  // sky location (lambda then beta)
  if (isnan(priorParams->fix_lambda)) {
    lambda = CubeToFlatPrior(Cube[i++], 0.0, 2.0 * M_PI);
  } else {
    lambda = priorParams->fix_lambda;
  }
  if (isnan(priorParams->fix_beta)) {
    beta = CubeToCosPrior(Cube[i++], -M_PI / 2.0, M_PI / 2.0);
  } else {
    beta = priorParams->fix_beta;
  }

  // polarization
  if (isnan(priorParams->fix_pol)) {
    pol = CubeToFlatPrior(Cube[i++], 0.0, M_PI);
  } else {
    pol = priorParams->fix_pol;
  }

  Cube[0] = m1;
  Cube[1] = m2;
  Cube[2] = tRef;
  Cube[3] = dist;
  Cube[4] = phase;
  Cube[5] = inc;
  Cube[6] = lambda;
  Cube[7] = beta;
  Cube[8] = pol;
}

/******************************************** getallparams routine ****************************************************/

void getallparams(double *Cube, int *ndim)
{
	getphysparams(Cube,ndim);
}

/******************************************** loglikelihood routine ****************************************************/

// Now an example, sample an egg box likelihood

// Input arguments
// ndim 						= dimensionality (total number of free parameters) of the problem
// npars 						= total number of free plus derived parameters
//
// Input/Output arguments
// Cube[npars] 						= on entry has the ndim parameters in unit-hypercube
//	 						on exit, the physical parameters plus copy any derived parameters you want to store with the free parameters
//
// Output arguments
// lnew 						= loglikelihood

/* Note: context must point to the LISASignal structure representing the injected signals */
void getLogLike(double *Cube, int *ndim, int *npars, double *lnew, void *context)
//void getLogLike(LISAParams *params, double *lnew, void *context)
{
  /* Convert Cube to physical parameters and check prior boundary */
  getallparams(Cube,ndim);
  if (PriorBoundaryCheck(priorParams, Cube)) {
    //
    printf("getLogLike outside of prior boundary.\n");
    *lnew = -DBL_MAX;
    return;
  }

  LISAParams templateparams;
  templateparams.m1 = Cube[0];
  templateparams.m2 = Cube[1];
  templateparams.tRef = Cube[2];
  templateparams.distance = Cube[3];
  templateparams.phiRef = Cube[4];
  templateparams.inclination = Cube[5];
  templateparams.lambda = Cube[6];
  templateparams.beta = Cube[7];
  templateparams.polarization = Cube[8];
  templateparams.nbmode = globalparams->nbmodetemp; /* Using the global parameter for the number of modes in templates */

  /* Note: context points to a LISAContext structure containing a LISASignal* */
  if(globalparams->tagint==0) {
    LISASignalCAmpPhase* injection = ((LISASignalCAmpPhase*) context);

    //TESTING
    //clock_t tbeg, tend;
    //tbeg = clock();
    *lnew = CalculateLogLCAmpPhase(&templateparams, injection) - logZdata;
    //tend = clock();
    //printf("time Likelihood: %g\n", (double) (tend-tbeg)/CLOCKS_PER_SEC);
    //
  }
  else if(globalparams->tagint==1) {
    LISAInjectionReIm* injection = ((LISAInjectionReIm*) context);

    //TESTING
    //clock_t tbeg, tend;
    //tbeg = clock();
    *lnew = CalculateLogLReIm(&templateparams, injection) - logZdata;
    //tend = clock();
    //printf("time Likelihood: %g\n", (double) (tend-tbeg)/CLOCKS_PER_SEC);
    //
  }

  //
  nblikelihoods++;
  printf("%d: %12e\n", nblikelihoods, *lnew);
  if(nblikelihoods==1) exit(0);
}


/************************************************* dumper routine ******************************************************/

// The dumper routine will be called every updInt*10 iterations
// MultiNest doesn not need to the user to do anything. User can use the arguments in whichever way he/she wants
//
//
// Arguments:
//
// nSamples 						= total number of samples in posterior distribution
// nlive 						= total number of live points
// nPar 						= total number of parameters (free + derived)
// physLive[1][nlive * (nPar + 1)] 			= 2D array containing the last set of live points (physical parameters plus derived parameters) along with their loglikelihood values
// posterior[1][nSamples * (nPar + 2)] 			= posterior distribution containing nSamples points. Each sample has nPar parameters (physical + derived) along with the their loglike value & posterior probability
// paramConstr[1][4*nPar]:
// paramConstr[0][0] to paramConstr[0][nPar - 1] 	= mean values of the parameters
// paramConstr[0][nPar] to paramConstr[0][2*nPar - 1] 	= standard deviation of the parameters
// paramConstr[0][nPar*2] to paramConstr[0][3*nPar - 1] = best-fit (maxlike) parameters
// paramConstr[0][nPar*4] to paramConstr[0][4*nPar - 1] = MAP (maximum-a-posteriori) parameters
// maxLogLike						= maximum loglikelihood value
// logZ							= log evidence value
// logZerr						= error on log evidence value
// context						void pointer, any additional information

void dumper(int *nSamples, int *nlive, int *nPar, double **physLive, double **posterior, double **paramConstr, double *maxLogLike, double *logZ, double *logZerr, void *context)
{
	// convert the 2D Fortran arrays to C arrays


	// the posterior distribution
	// postdist will have nPar parameters in the first nPar columns & loglike value & the posterior probability in the last two columns

	int i, j;

	double postdist[*nSamples][*nPar + 2];
	for( i = 0; i < *nPar + 2; i++ )
		for( j = 0; j < *nSamples; j++ )
			postdist[j][i] = posterior[0][i * (*nSamples) + j];



	// last set of live points
	// pLivePts will have nPar parameters in the first nPar columns & loglike value in the last column

	double pLivePts[*nlive][*nPar + 1];
	for( i = 0; i < *nPar + 1; i++ )
		for( j = 0; j < *nlive; j++ )
			pLivePts[j][i] = physLive[0][i * (*nlive) + j];
}

/***********************************************************************************************************************/




/************************************************** Main program *******************************************************/



int main(int argc, char *argv[])
{
  int myid = 0;
#ifdef PARALLEL
 	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
#endif

  /*********** Addendum *************/

	//
	int nblikelihoods = 0;
	int tagprint = 0;

  /* Initialize structs for holding various options */
  LISARunParams runParams;
  injectedparams = (LISAParams*) malloc(sizeof(LISAParams));
  memset(injectedparams, 0, sizeof(LISAParams));
  globalparams = (LISAGlobalParams*) malloc(sizeof(LISAGlobalParams));
  memset(globalparams, 0, sizeof(LISAGlobalParams));
  priorParams = (LISAPrior*) malloc(sizeof(LISAPrior));
  memset(priorParams, 0, sizeof(LISAPrior));
  
  /* Parse commandline to read parameters of injection - copy the number of modes demanded for the injection */
  parse_args_LISA(argc, argv, injectedparams, globalparams, priorParams, &runParams);
  injectedparams->nbmode = globalparams->nbmodeinj;

  /* Initialize the data structure for the injection */
  LISASignalCAmpPhase* injectedsignalCAmpPhase = NULL;
  LISAInjectionReIm* injectedsignalReIm = NULL;
  if(globalparams->tagint==0) {
    LISASignalCAmpPhase_Init(&injectedsignalCAmpPhase);
  }
  else if(globalparams->tagint==1) {
    LISAInjectionReIm_Init(&injectedsignalReIm);
  }

  /* Generate the injection */
  if(globalparams->tagint==0) {
    LISAGenerateSignalCAmpPhase(injectedparams, injectedsignalCAmpPhase);
  }
  else if(globalparams->tagint==1) {
    LISAGenerateInjectionReIm(injectedparams, globalparams->fmin, globalparams->nbptsoverlap, injectedsignalReIm);
  }

  /* Define SNRs */
  double SNRA, SNRE, SNRT;
  if(globalparams->tagint==0) {
    SNRA = injectedsignalCAmpPhase->TDIAhh;
    SNRE = injectedsignalCAmpPhase->TDIEhh;
    SNRT = injectedsignalCAmpPhase->TDIThh;
  }
  else if(globalparams->tagint==1) {
    SNRA = sqrt(FDOverlapReImvsReIm(injectedsignalReIm->TDIASignal, injectedsignalReIm->TDIASignal, injectedsignalReIm->noisevaluesA));
    SNRE = sqrt(FDOverlapReImvsReIm(injectedsignalReIm->TDIESignal, injectedsignalReIm->TDIESignal, injectedsignalReIm->noisevaluesE));
    SNRT = sqrt(FDOverlapReImvsReIm(injectedsignalReIm->TDITSignal, injectedsignalReIm->TDITSignal, injectedsignalReIm->noisevaluesT));
  }
  double SNR = sqrt(SNRA*SNRA + SNRE*SNRE + SNRT*SNRT);

  /* Rescale distance to match SNR */
  if (!isnan(priorParams->snr_target)) {
    if (myid == 0) printf("Rescaling the distance to obtain a network SNR of %g\n", priorParams->snr_target);
    injectedparams->distance *= SNR / priorParams->snr_target;
    if (myid == 0) printf("New distance = %g Mpc\n", injectedparams->distance);
    if(globalparams->tagint==0) {
      LISAGenerateSignalCAmpPhase(injectedparams, injectedsignalCAmpPhase);
      SNRA = injectedsignalCAmpPhase->TDIAhh;
      SNRE = injectedsignalCAmpPhase->TDIEhh;
      SNRT = injectedsignalCAmpPhase->TDIThh;
    }
    else if(globalparams->tagint==1) {
      LISAGenerateInjectionReIm(injectedparams, globalparams->fmin, globalparams->nbptsoverlap, injectedsignalReIm);
      SNRA = sqrt(FDOverlapReImvsReIm(injectedsignalReIm->TDIASignal, injectedsignalReIm->TDIASignal, injectedsignalReIm->noisevaluesA));
      SNRE = sqrt(FDOverlapReImvsReIm(injectedsignalReIm->TDIESignal, injectedsignalReIm->TDIESignal, injectedsignalReIm->noisevaluesE));
      SNRT = sqrt(FDOverlapReImvsReIm(injectedsignalReIm->TDITSignal, injectedsignalReIm->TDITSignal, injectedsignalReIm->noisevaluesT));
    }
    SNR = sqrt(SNRA*SNRA + SNRE*SNRE + SNRT*SNRT);
  }

  /* Print SNR */
  if (myid == 0) {
    printf("SNR A:     %g\n", SNRA);
    printf("SNR E:     %g\n", SNRE);
    printf("SNR T:     %g\n", SNRT);
    printf("SNR Total: %g\n", SNR);
  }

  /* Calculate logL of data */
  /*double dist_store = injectedparams->distance;
    injectedparams->distance = 1.0e9;
    logZdata = CalculateLogL(injectedparams, injectedsignal);
    printf("logZdata = %lf\n", logZdata);
    injectedparams->distance = dist_store;*/
  logZdata = 0.0;
  double logZtrue = 0.;
  if(globalparams->tagint==0) {
    logZtrue = CalculateLogLCAmpPhase(injectedparams, injectedsignalCAmpPhase);
  }
  else if(globalparams->tagint==1) {
    logZtrue = CalculateLogLReIm(injectedparams, injectedsignalReIm);
  }
  if (myid == 0) printf("logZtrue = %lf\n", logZtrue-logZdata);

  /* Set the context pointer */
  void *context = NULL;
  if(globalparams->tagint==0) {
    context = injectedsignalCAmpPhase;
  }
  else if(globalparams->tagint==1) {
    context = injectedsignalReIm;
  }

  //
  //printReImFrequencySeries(injectedsignalReIm->TDIASignal, 0, 50);
  //exit(0);
  tagprint = 1;
  //TESTING
  //exit(0);

  int ndims = 9;

  /* check for fixed parameters */
  if (!isnan(priorParams->fix_m1))
    ndims--;
  if (!isnan(priorParams->fix_m2))
    ndims--;
  if (!isnan(priorParams->fix_dist))
    ndims--;
  if (!isnan(priorParams->fix_inc))
    ndims--;
  if (!isnan(priorParams->fix_phase))
    ndims--;
  if (!isnan(priorParams->fix_pol))
    ndims--;
  if (!isnan(priorParams->fix_lambda))
    ndims--;
  if (!isnan(priorParams->fix_beta))
    ndims--;
  if (!isnan(priorParams->fix_time))
    ndims--;

  if (ndims == 0) {
    LISAParams templateparams;
    templateparams.m1 = priorParams->fix_m1;
    templateparams.m2 = priorParams->fix_m2;
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

	/********** End of addendum ****************/

	/********** Test ****************/
	/* double l; */
	/* l = CalculateLogLReIm(injectedparams, injectedsignalReIm); */
	/* printf("LogLikelihood: %g\n", l); */
	/* free(injectedparams); */
	/* LISAInjectionReIm_Cleanup(injectedsignalReIm); */
	/* exit(0); */
	/********** End of test ****************/

	int i;

	// set the MultiNest sampling parameters

	int mmodal = 0;					// do mode separation?

	int ceff = 0;					// run in constant efficiency mode?

	int nlive = runParams.nlive;				// number of live points

	double efr = runParams.eff;				// set the required efficiency

	double tol = runParams.tol;				// tol, defines the stopping criteria

	//int ndims = 9;					// dimensionality (no. of free parameters)

	int nPar = 9;					// total no. of parameters including free & derived parameters

	int nClsPar = (int) (fmin(2.,ndims));				// no. of parameters to do mode separation on

	int updInt = 50;				// after how many iterations feedback is required & the output files should be updated
							// note: posterior files are updated & dumper routine is called after every updInt*10 iterations

	double Ztol = -1E90;				// all the modes with logZ < Ztol are ignored

	int maxModes = 1;				// expected max no. of modes (used only for memory allocation)

	int pWrap[ndims];				// which parameters to have periodic boundary conditions?
	for(i = 0; i < ndims; i++) pWrap[i] = 0;
  pWrap[4] = pWrap[6] = pWrap[8] = 1;

	strcpy(root, runParams.outroot);			// root for output files
	strcpy(networkinputs, runParams.netfile);			// file with input parameters for network training

	int seed = -1;					// random no. generator seed, if < 0 then take the seed from system clock

	int fb = 1;					// need feedback on standard output?

	resume = runParams.resume;					// resume from a previous job?

	int outfile = 1;				// write output files?

	int initMPI = 0;				// initialize MPI routines?, relevant only if compiling with MPI
							// set it to F if you want your main program to handle MPI initialization

	logZero = -1E90;				// points with loglike < logZero will be ignored by MultiNest

	int maxiter = 0;				// max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it
							// has done max no. of iterations or convergence criterion (defined through tol) has been satisfied

	// void *context = 0;				// not required by MultiNest, any additional information user wants to pass

	doBAMBI = runParams.bambi;					// BAMBI?

	useNN = 0;

	// calling MultiNest

	BAMBIrun(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI,
	logZero, maxiter, LogLikeFctn, dumper, BAMBIfctn, context);

  free(injectedparams);
  free(priorParams);

#ifdef PARALLEL
 	MPI_Finalize();
#endif
}

/***********************************************************************************************************************/
