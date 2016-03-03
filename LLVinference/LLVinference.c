#include "LLVinference.h"

/******************************************** getphysparams routine ****************************************************/

/* Order of the parameters (fixed): m1, m2, tRef, dist, phase, inc, ra, dec, pol */
void getphysparams(double *Cube, int *ndim)
{
	int i = 0;
  double m1=0., m2=0., tRef=0., dist=0., phase=0., inc=0., ra=0., dec=0., pol=0.;

	/* Note: we use here the order for the cube parameters */
  /* Order of the 9 original parameters (fixed): m1, m2, tRef, dist, phase, inc, ra, dec, pol */
  /* Order of the 9 cube parameters (modified for clustering): ra, dec, tRef, phase, pol, inc, dist, m1, m2 */

	/* Sky location (ra then dec) */
	if (isnan(priorParams->fix_ra)) {
		ra = CubeToFlatPrior(Cube[i++], priorParams->ra_min, priorParams->ra_max);
	} else {
		ra = priorParams->fix_ra;
  }
  if (isnan(priorParams->fix_dec)) {
    dec = CubeToCosPrior(Cube[i++], priorParams->dec_min, priorParams->dec_max);
  } else {
    dec = priorParams->fix_dec;
  }

	/* Time */
  if (isnan(priorParams->fix_time)) {
    tRef = CubeToFlatPrior(Cube[i++], injectedparams->tRef - priorParams->deltaT,
           injectedparams->tRef + priorParams->deltaT);
  } else {
    tRef = priorParams->fix_time;
  }

  /* Orbital phase */
  if (isnan(priorParams->fix_phase)) {
    phase = CubeToFlatPrior(Cube[i++], priorParams->phase_min, priorParams->phase_max);
  } else {
    phase = priorParams->fix_phase;
  }

  /* Polarization */
  if (isnan(priorParams->fix_pol)) {
    pol = CubeToFlatPrior(Cube[i++], priorParams->pol_min, priorParams->pol_max);
  } else {
    pol = priorParams->fix_pol;
  }

  /* Inclination */
  if (isnan(priorParams->fix_inc)) {
    inc = CubeToSinPrior(Cube[i++], priorParams->inc_min, priorParams->inc_max);
  } else {
    inc = priorParams->fix_inc;
  }

  /* Distance - two priors allowed: r^2-weighted or flat */
  if (isnan(priorParams->fix_dist)) {
    if(priorParams->flat_distprior) {
      dist = CubeToFlatPrior(Cube[i++], priorParams->dist_min, priorParams->dist_max);
    }
    else {
      dist = CubeToPowerPrior(2.0, Cube[i++], priorParams->dist_min, priorParams->dist_max);
    }
  } else {
    dist = priorParams->fix_dist;
  }

  /* Component masses */
  if (isnan(priorParams->fix_m1)) {
    if (isnan(priorParams->fix_m2)) {
      m1 = CubeToFlatPrior(Cube[i++], priorParams->comp_min, priorParams->comp_max);
      m2 = CubeToFlatPrior(Cube[i++], priorParams->comp_min, priorParams->comp_max);
      /*if (m2 > m1) {
        double tmp = m1;
        m1 = m2;
        m2 = tmp;
      }*/
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

  Cube[0] = m1;
  Cube[1] = m2;
  Cube[2] = tRef;
  Cube[3] = dist;
  Cube[4] = phase;
  Cube[5] = inc;
  Cube[6] = ra;
  Cube[7] = dec;
  Cube[8] = pol;
}

void getcubeparams(double* Cube, int ndim, LLVParams* params, int* freeparamsmap)
{
  int i = 0;
  double m1 = params->m1;
  double m2 = params->m2;
  double tRef = params->tRef;
  double dist = params->distance;
  double phase = params->phiRef;
  double inc = params->inclination;
  double ra = params->ra;
  double dec = params->dec;
  double pol = params->polarization;

  /* Note: freeparamsmap has indices in the order of the (free) cube parameters, and values in the indices of physical parameters */
  for(int i=0; i<ndim; i++) {
    if(freeparamsmap[i]==0) Cube[i] = FlatPriorToCube(m1, priorParams->comp_min, priorParams->comp_max);
    if(freeparamsmap[i]==1) Cube[i] = FlatPriorToCube(m2, priorParams->comp_min, priorParams->comp_max);
    if(freeparamsmap[i]==2) Cube[i] = FlatPriorToCube(tRef, injectedparams->tRef - priorParams->deltaT, injectedparams->tRef + priorParams->deltaT);
    if(freeparamsmap[i]==3) Cube[i] = PowerPriorToCube(2., dist, priorParams->dist_min, priorParams->dist_max);
    if(freeparamsmap[i]==4) Cube[i] = FlatPriorToCube(phase, priorParams->phase_min, priorParams->phase_max);
    if(freeparamsmap[i]==5) Cube[i] = SinPriorToCube(inc, priorParams->inc_min, priorParams->inc_max);
    if(freeparamsmap[i]==6) Cube[i] = FlatPriorToCube(ra, priorParams->ra_min, priorParams->ra_max);
    if(freeparamsmap[i]==7) Cube[i] = CosPriorToCube(dec, priorParams->dec_min, priorParams->dec_max);
    if(freeparamsmap[i]==8) Cube[i] = FlatPriorToCube(pol, priorParams->pol_min, priorParams->pol_max);
  }
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
void getLogLike(double *Cube, int *ndim, int *npars, double *lnew, void *context)
//void getLogLike(LLVParams *params, double *lnew, void *context)
{
  /* Convert Cube to physical parameters and check prior boundary */
  getallparams(Cube,ndim);
  if (PriorBoundaryCheck(priorParams, Cube)) {
    *lnew = -DBL_MAX;
    return;
  }

  LLVParams templateparams;
  templateparams.m1 = Cube[0];
  templateparams.m2 = Cube[1];
  templateparams.tRef = Cube[2];
  templateparams.distance = Cube[3];
  templateparams.phiRef = Cube[4];
  templateparams.inclination = Cube[5];
  templateparams.ra = Cube[6];
  templateparams.dec = Cube[7];
  templateparams.polarization = Cube[8];
  templateparams.nbmode = globalparams->nbmodetemp;

  /* Note: context points to a LLVContext structure containing a LLVSignal* */
  LLVInjectionCAmpPhase* injection = ((LLVInjectionCAmpPhase*) context);

  *lnew = CalculateLogLCAmpPhase(&templateparams, injection) - logZdata;
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

	/* Initialize structs for holding various options */
  LLVRunParams runParams;
  injectedparams = (LLVParams*) malloc(sizeof(LLVParams));
  memset(injectedparams, 0, sizeof(LLVParams));
  globalparams = (LLVGlobalParams*) malloc(sizeof(LLVGlobalParams));
  memset(globalparams, 0, sizeof(LLVGlobalParams));
  priorParams = (LLVPrior*) malloc(sizeof(LLVPrior));
  memset(priorParams, 0, sizeof(LLVPrior));

  /* Parse commandline to read parameters of injection - copy the number of modes demanded for the injection  */
  parse_args_LLV(argc, argv, injectedparams, globalparams, priorParams, &runParams);
  injectedparams->nbmode = globalparams->nbmodeinj;

  /* Load and initialize the detector noise */
  LLVSimFD_Noise_Init_ParsePath();

  /* Initialize the data structure for the injection */
  LLVInjectionCAmpPhase* injectedsignal = NULL;
  LLVInjectionCAmpPhase_Init(&injectedsignal);

  /* Generate the injection */
	LLVGenerateInjectionCAmpPhase(injectedparams, injectedsignal);

  /* Rescale distance to match SNR */
  if (!isnan(priorParams->snr_target)) {
    if (myid == 0) printf("Rescaling the distance to obtain a network SNR of %g\n", priorParams->snr_target);
    injectedparams->distance *= sqrt(injectedsignal->LLVss) / priorParams->snr_target;
    if (myid == 0) printf("New distance = %g Mpc\n", injectedparams->distance);
    LLVGenerateInjectionCAmpPhase(injectedparams, injectedsignal);
  }

  /* print SNRs */
  if (myid == 0) {
    printf("SNR Network: %g\n", injectedsignal->LLVss);
  }

  /* Calculate logL of data */
  /*double dist_store = injectedparams->distance;
    injectedparams->distance = 1.0e9;
    logZdata = CalculateLogL(injectedparams, injectedsignal);
    printf("logZdata = %lf\n", logZdata);
    injectedparams->distance = dist_store;*/
  logZdata = 0.0;
  double logZtrue = CalculateLogLCAmpPhase(injectedparams, injectedsignal);
  if (myid == 0) printf("logZtrue = %lf\n", logZtrue-logZdata);

  /* Set the context pointer */
  void *context = injectedsignal;

  int ndim = 9;

  /* check for parameters pinned to injected values */
  if (priorParams->pin_m1)
    priorParams->fix_m1 = injectedparams->m1;
  if (priorParams->pin_m2)
    priorParams->fix_m2 = injectedparams->m2;
  if (priorParams->pin_dist)
    priorParams->fix_dist = injectedparams->distance;
  if (priorParams->pin_inc)
    priorParams->fix_inc = injectedparams->inclination;
  if (priorParams->pin_phase)
    priorParams->fix_phase = injectedparams->phiRef;
  if (priorParams->pin_pol)
    priorParams->fix_pol = injectedparams->polarization;
  if (priorParams->pin_ra)
    priorParams->fix_ra = injectedparams->ra;
  if (priorParams->pin_dec)
    priorParams->fix_dec = injectedparams->dec;
  if (priorParams->pin_time)
    priorParams->fix_time = injectedparams->tRef;

  /* Check for fixed parameters, and build the map from the free parameters to the orignal 9 parameters */
  /* Order of the 9 original parameters (fixed): m1, m2, tRef, dist, phase, inc, ra, dec, pol */
  int freeparams[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
  if (!isnan(priorParams->fix_m1))    { ndim--; freeparams[0] = 0; }
  if (!isnan(priorParams->fix_m2))    { ndim--; freeparams[1] = 0; }
  if (!isnan(priorParams->fix_time))  { ndim--; freeparams[2] = 0; }
  if (!isnan(priorParams->fix_dist))  { ndim--; freeparams[3] = 0; }
  if (!isnan(priorParams->fix_phase)) { ndim--; freeparams[4] = 0; }
  if (!isnan(priorParams->fix_inc))   { ndim--; freeparams[5] = 0; }
  if (!isnan(priorParams->fix_ra))    { ndim--; freeparams[6] = 0; }
  if (!isnan(priorParams->fix_dec))   { ndim--; freeparams[7] = 0; }
  if (!isnan(priorParams->fix_pol))   { ndim--; freeparams[8] = 0; }
  int* freeparamsmap = malloc(ndim*sizeof(int));
  int counter = 0;
  for(int i=0; i<ndim; i++) {
    while(freeparams[counter]==0) counter++;
    freeparamsmap[i] = counter;
    counter++;
  }

  if (ndim == 0) {
    LLVParams templateparams;
    templateparams.m1 = priorParams->fix_m1;
    templateparams.m2 = priorParams->fix_m2;
    templateparams.tRef = priorParams->fix_time;
    templateparams.distance = priorParams->fix_dist;
    templateparams.phiRef = priorParams->fix_phase;
    templateparams.inclination = priorParams->fix_inc;
    templateparams.ra = priorParams->fix_ra;
    templateparams.dec = priorParams->fix_dec;
    templateparams.polarization = priorParams->fix_pol;
    templateparams.nbmode = globalparams->nbmodetemp;

    double logL = CalculateLogLCAmpPhase(&templateparams, injectedsignal);
    if (myid == 0) printf("logL = %lf\n", logL);

    free(injectedparams);
    free(priorParams);

#ifdef PARALLEL
    MPI_Finalize();
#endif

    exit(0);
  }

	/********** End of addendum ****************/

	/********** Test ****************/
	/*double l;
	  getLogLike(injectedparams, &l, context);
	  printf("LogLikelihood: %g\n", l);
	  free(injectedparams);
	  LLVSignal_Cleanup(injectedsignal);*/
	/********** End of test ****************/

	int i;

	// set the MultiNest sampling parameters

	int mmodal = 0;					// do mode separation?

	int ceff = 0;					// run in constant efficiency mode?

	int nlive = runParams.nlive;				// number of live points

	double efr = runParams.eff;				// set the required efficiency

	double tol = runParams.tol;				// tol, defines the stopping criteria

	//int ndim = 9;					// dimensionality (no. of free parameters)

	int nPar = 13;					// total no. of parameters including free & derived parameters

	int nClsPar = (int) (fmin(2.,ndim));				// no. of parameters to do mode separation on

	int updInt = 50;				// after how many iterations feedback is required & the output files should be updated
							// note: posterior files are updated & dumper routine is called after every updInt*10 iterations

	double Ztol = -1E90;				// all the modes with logZ < Ztol are ignored

	int maxModes = 1;				// expected max no. of modes (used only for memory allocation)

	int pWrap[ndim];				// which parameters to have periodic boundary conditions?
	for(i = 0; i < ndim; i++) {
	  if(freeparamsmap[i]==4 || freeparamsmap[i]==6 || freeparamsmap[i]==8) pWrap[i] = 1;
	  else pWrap[i] = 0;
	}

	strcpy(root, runParams.outroot);		// root for output files
	strcpy(networkinputs, runParams.netfile);	// file with input parameters for network training

	int seed = -1;					// random no. generator seed, if < 0 then take the seed from system clock

	int fb = 1;					// need feedback on standard output?

	resume = runParams.resume;			// resume from a previous job?

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

	BAMBIrun(mmodal, ceff, nlive, tol, efr, ndim, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI,
	logZero, maxiter, LogLikeFctn, dumper, BAMBIfctn, context);

  free(injectedparams);
  free(priorParams);

#ifdef PARALLEL
 	MPI_Finalize();
#endif
}

/***********************************************************************************************************************/
