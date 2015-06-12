#include "LLVinference.h"

// Global parameters
LLVParams* injectedparams = NULL, templateparams = NULL;
LLVPrior* priorParams = NULL;

/************ Functions to initalize and clean up structure for the signals ************/

void LLVSignal_Cleanup(LLVSignal* signal) {
  if(signal->LHOSignal) ListmodesCAmpPhaseFrequencySeries_Destroy(signal->LHOSignal);
  if(signal->LLOSignal) ListmodesCAmpPhaseFrequencySeries_Destroy(signal->LLOSignal);
  if(signal->VIRGOSignal) ListmodesCAmpPhaseFrequencySeries_Destroy(signal->VIRGOSignal);
  free(signal);
}
void LLVSignal_Init(LLVSignal** signal) {
  if(!signal) exit(1);
  /* Create storage for structures */
  if(!*signal) *signal = malloc(sizeof(LLVSignal));
  else
  {
    LLVSignal_Cleanup(*signal);
  }
  (*signal)->LHOSignal = NULL;
  (*signal)->LLOSignal = NULL;
  (*signal)->VIRGOSignal = NULL;
}

/************ Functions for LLV parameters, injection, likelihood ************/

/* Parse command line and return a newly allocated LLVParams object */
/* Masses are input in solar masses and distances in Mpc - converted in SI for the internals */
LLVParams* parse_args_LLV(ssize_t argc, char **argv) {
    char help[] = ""

    ssize_t i;
    LLVParams* params;
    params = (LLVParams*) malloc(sizeof(LLVParams));
    memset(params, 0, sizeof(LLVParams));

    /* Set default values to the arguments */
    params->tRef = 0.;
    params->phiRef = 0.;
    params->m1 = 20.;
    params->m2 = 10.;
    params->distance = 100. * 1e6;
    params->ra = 0.;
    params->dec = 0.;
    params->inclination = PI/3.;
    params->polarization = 0.;
    params->fRef = 0.;
    params->nbmode = 1;

    /* Consume command line */
    for (i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--tRef") == 0) {
            params->tRef = atof(argv[++i]);
        } else if (strcmp(argv[i], "--phiRef") == 0) {
            params->phiRef = atof(argv[++i]);
        } else if (strcmp(argv[i], "--m1") == 0) {
            params->m1 = atof(argv[++i]) * MSUN_SI;
        } else if (strcmp(argv[i], "--m2") == 0) {
            params->m2 = atof(argv[++i]) * MSUN_SI;
        } else if (strcmp(argv[i], "--distance") == 0) {
            params->distance = atof(argv[++i]) * 1e6 * PC_SI;
        } else if (strcmp(argv[i], "--ra") == 0) {
            params->ra = atof(argv[++i]);
        } else if (strcmp(argv[i], "--dec") == 0) {
            params->dec = atof(argv[++i]);
        } else if (strcmp(argv[i], "--inclination") == 0) {
            params->inclination = atof(argv[++i]);
        } else if (strcmp(argv[i], "--polarization") == 0) {
            params->polarization = atof(argv[++i]);
        } else if (strcmp(argv[i], "--fRef") == 0) {
            params->fRef = atof(argv[++i]);
        } else if (strcmp(argv[i], "--nbmode") == 0) {
            params->nbmode = atof(argv[++i]);
        } else {
<<<<<<< HEAD
            //printf("Error: invalid option: %s\n", argv[i]);
            //goto fail;
=======
	  printf("Error: invalid option: %s\n", argv[i]);
	  goto fail;
>>>>>>> daf9ea2ca361b5398f3bd2e37fdb657ad901b781
        }
    }

    return params;

    /*fail:
    free(params);
    exit(1);*/
}

/* Function generating a LLV signal from LLV parameters */
int LLVGenerateSignal(
  struct tagLLVParams* params,   /* Input: set of LLV parameters of the signal */
  struct tagLLVSignal* signal)   /* Output: structure for the generated signal */
{
  ListmodesCAmpPhaseFrequencySeries* listROM = NULL;
  ListmodesCAmpPhaseFrequencySeries* listLHO = NULL;
  ListmodesCAmpPhaseFrequencySeries* listLLO = NULL;
  ListmodesCAmpPhaseFrequencySeries* listVIRGO = NULL;

  /* Should add error checking ? */
  /* Generate the waveform with the ROM */
  SimEOBNRv2HMROM(&listROM, params->nbmode, params->tRef, params->phiRef, params->fRef, (params->m1)*MSUN_SI, (params->m2)*MSUN_SI, (params->distance)*PC_SI);
  /* Process the waveform through the LLV response */
  /* TO BE MODIFIED FOR RA DEC */
  LLVSimFDResponse(&listROM, &listLHO, params->inclination, params->ra, params->dec, params->polarization, LHO);
  LLVSimFDResponse(&listROM, &listLLO, params->inclination, params->ra, params->dec, params->polarization, LLO);
  LLVSimFDResponse(&listROM, &listVIRGO, params->inclination, params->ra, params->dec, params->polarization, VIRGO);

  /* Precompute the inner products (h|h) and (s|s) */
  double LHOhh = FDListmodesOverlap(listLHO, listLHO, NoiseSnLHO, __LLVSimFD_LHONoise_fLow, __LLVSimFD_LHONoise_fHigh);
  double LLOhh = FDListmodesOverlap(listLLO, listLLO, NoiseSnLLO, __LLVSimFD_LLONoise_fLow, __LLVSimFD_LLONoise_fHigh);
  double VIRGOhh = FDListmodesOverlap(listVIRGO, listVIRGO, NoiseSnVIRGO, __LLVSimFD_VIRGONoise_fLow, __LLVSimFD_VIRGONoise_fHigh);

  /* Output and clean up */
  signal->LHOSignal = listLHO;
  signal->LLOSignal = listLLO;
  signal->VIRGOSignal = listVIRGO;
  signal->LHOhh = LHOhh;
  signal->LLOhh = LLOhh;
  signal->VIRGOhh = VIRGOhh;

  ListmodesCAmpPhaseFrequencySeries_Destroy(listROM);
  return SUCCESS;
}

/******************************************** getphysparams routine ****************************************************/

void getphysparams(double *Cube, int *ndim)
{
	int i = 0;

  // component masses
  double m1 = CubeToFlatPrior(Cube[i++], priorParams->comp_min, priorParams->comp_max);
  double m2 = CubeToFlatPrior(Cube[i++], priorParams->comp_min, priorParams->comp_max);
  if (m2 > m1) {
    double tmp = m1;
    m1 = m2;
    m2 = tmp;
  }
  templateparams->m1 = m1;
  templateparams->m2 = m2;

  // time
  templateparams->tRef = CubeToFlatPrior(Cube[i++], injectedparams->tRef - priorParams->deltaT, injectedparams->tRef + priorParams->deltaT);

  // distance
  templateparams->distance = CubeToPowerPrior(2.0, Cube[i++], priorParams->dist_min, priorParams->dist_max);

  // orbital phase
  templateparams->phiRef = CubeToFlatPrior(Cube[i++], 0.0, 2.0 * M_PI);

  // inclination
  templateparams->inclination = CubeToSinPrior(Cube[i++], 0.0, M_PI);

  // sky location
  templateparams->ra = CubeToFlatPrior(Cube[i++], 0.0, 2.0 * M_PI);
  templateparams->dec = CubeToCosPrior(Cube[i++], -M_PI / 2.0, M_PI / 2.0);

  // polarization
  templateparams->polarization = CubeToFlatPrior(Cube[i++], 0.0, M_PI);

  // put all values back in Cube as well
  Cube[0] = m1;
  Cube[1] = m2;
  Cube[2] = templateparams->tRef;
  Cube[3] = templateparams->distance;
  Cube[4] = templateparams->phiRef;
  Cube[5] = templateparams->inclination;
  Cube[6] = templateparams->ra;
  Cube[7] = templateparams-dec;
  Cube[8] = templateparams->polarization;
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

/* Note: context must point to the LLVSignal structure representing the injected signals */
//void getLogLike(double *Cube, int *ndim, int *npars, double *lnew, void *context)
void getLogLike(LLVParams *params, double *lnew, void *context)
{
  /*  */
  /*getallparams(Cube,ndim);
  if (PriorBoundaryCheck(Cube)) {
    *lnew = -DBL_MAX;
    return;
  } */

  /* Note: context points to a LLVContext structure containing a LLVSignal* */
  LLVSignal* injection = ((LLVSignal*) context);

  /* Generating the signal in the three detectors for the input parameters */
  LLVSignal* generatedsignal = NULL;
  LLVSignal_Init(&generatedsignal);
  LLVGenerateSignal(params, generatedsignal);

  /* Computing the likelihood for each detector */
  double loglikelihoodLHO = FDLogLikelihood(injection->LHOSignal, generatedsignal->LHOSignal, NoiseSnLHO, __LLVSimFD_LHONoise_fLow, __LLVSimFD_LHONoise_fHigh, injection->LHOhh, generatedsignal->LHOhh);
  double loglikelihoodLLO = FDLogLikelihood(injection->LLOSignal, generatedsignal->LLOSignal, NoiseSnLLO, __LLVSimFD_LLONoise_fLow, __LLVSimFD_LLONoise_fHigh, injection->LLOhh, generatedsignal->LLOhh);
  double loglikelihoodVIRGO = FDLogLikelihood(injection->VIRGOSignal, generatedsignal->VIRGOSignal, NoiseSnVIRGO, __LLVSimFD_VIRGONoise_fLow, __LLVSimFD_VIRGONoise_fHigh, injection->VIRGOhh, generatedsignal->VIRGOhh);

  /* Output: value of the loglikelihood for the combined signals, assuming noise independence */
  *lnew = loglikelihoodLHO + loglikelihoodLLO + loglikelihoodVIRGO;

  /* Clean up */
  LLVSignal_Cleanup(generatedsignal);
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
#ifdef PARALLEL
 	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
#endif

	/*********** Addendum *************/

	/* Parse commandline to read parameters of injection */
<<<<<<< HEAD
	printf("a\n");
	injectedparams = parse_args_LLV(argc, argv);
=======
	LLVParams* injectedparams = parse_args_LLV(argc, argv);
>>>>>>> daf9ea2ca361b5398f3bd2e37fdb657ad901b781

	/* Load and initialize the detector noise */
	LLVSimFD_Noise_Init_ParsePath();

	/* Initialize the data structure for the injection */
	LLVSignal* injectedsignal = NULL;
	LLVSignal_Init(&injectedsignal);

	/* Generate the injection */
	LLVGenerateSignal(injectedparams, injectedsignal);

	/* Set the context pointer */
	void *context = injectedsignal;

	/********** End of addendum ****************/

	/********** Test ****************/
	double l;
	getLogLike(injectedparams, &l, context);
	printf("LogLikelihood: %g\n", l);
	free(injectedparams);
	LLVSignal_Cleanup(injectedsignal);
	/********** End of test ****************/

  /* Initialize the prior */
  priorParams = LLVInitializePrior(argc, argv);

  /* Initialize Bambi run settings */

/*
	int i;

	// set the MultiNest sampling parameters

	int mmodal = 0;					// do mode separation?

	int ceff = 0;					// run in constant efficiency mode?

	int nlive = 4000;				// number of live points

	double efr = 0.8;				// set the required efficiency

	double tol = 0.5;				// tol, defines the stopping criteria

	int ndims = 2;					// dimensionality (no. of free parameters)

	int nPar = 2;					// total no. of parameters including free & derived parameters

	int nClsPar = 2;				// no. of parameters to do mode separation on

	int updInt = 4000;				// after how many iterations feedback is required & the output files should be updated
							// note: posterior files are updated & dumper routine is called after every updInt*10 iterations

	double Ztol = -1E90;				// all the modes with logZ < Ztol are ignored

	int maxModes = 1;				// expected max no. of modes (used only for memory allocation)

	int pWrap[ndims];				// which parameters to have periodic boundary conditions?
	for(i = 0; i < ndims; i++) pWrap[i] = 0;

	strcpy(root, "chains/eggboxC_");			// root for output files
	strcpy(networkinputs, "example_eggbox_C/eggbox_net.inp");			// file with input parameters for network training

	int seed = -1;					// random no. generator seed, if < 0 then take the seed from system clock

	int fb = 1;					// need feedback on standard output?

	resume = 0;					// resume from a previous job?

	int outfile = 1;				// write output files?

	int initMPI = 0;				// initialize MPI routines?, relevant only if compiling with MPI
							// set it to F if you want your main program to handle MPI initialization

	logZero = -1E90;				// points with loglike < logZero will be ignored by MultiNest

	int maxiter = 0;				// max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it
							// has done max no. of iterations or convergence criterion (defined through tol) has been satisfied

	// void *context = 0;				// not required by MultiNest, any additional information user wants to pass

	doBAMBI = 1;					// BAMBI?

	useNN = 0;

	// calling MultiNest

	BAMBIrun(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI,
	logZero, maxiter, LogLikeFctn, dumper, BAMBIfctn, context);

*/

#ifdef PARALLEL
 	MPI_Finalize();
#endif
}

/***********************************************************************************************************************/
