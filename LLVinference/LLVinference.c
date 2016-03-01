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

  Cube[9] = Cube[0] + Cube[1];
  Cube[10] = Cube[0] / Cube[1];
  Cube[11] = Cube[10] / pow(1.0 + Cube[10], 2.0);
  Cube[12] = Cube[9] * pow(Cube[11], 0.6);
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
  if(globalparams->tagint==0) {
    LLVInjectionCAmpPhase* injection = ((LLVInjectionCAmpPhase*) context);

    //TESTING
    //clock_t tbeg, tend;
    //tbeg = clock();
    *lnew = CalculateLogLCAmpPhase(&templateparams, injection);
    //tend = clock();
    //printf("time Likelihood: %g\n", (double) (tend-tbeg)/CLOCKS_PER_SEC);
    //
  }
  else if(globalparams->tagint==1) {
    LLVInjectionReIm* injection = ((LLVInjectionReIm*) context);

    //TESTING
    //clock_t tbeg, tend;
    //tbeg = clock();
    *lnew = CalculateLogLReIm(&templateparams, injection);
    //tend = clock();
    //printf("time Likelihood: %g\n", (double) (tend-tbeg)/CLOCKS_PER_SEC);
    //
  }
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
	if(myid == 0) print_parameters_to_file_LLV(injectedparams, globalparams, priorParams, &runParams);

  /* Load and initialize the detector noise */
  LLVSimFD_Noise_Init_ParsePath();

	//
	printf("%p, %p\n", injectedparams, priorParams);

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

	/* Define SNRs */
	double SNR123, SNR1, SNR2, SNR3;
	if(globalparams->tagint==0) {
		SNR123 = sqrt(injectedsignalCAmpPhase->LLVss);
	}
	else if(globalparams->tagint==1) {
		SNR1 = sqrt(FDOverlapReImvsReIm(injectedsignalReIm->LHOSignal, injectedsignalReIm->LHOSignal, injectedsignalReIm->noisevaluesLHO));
		SNR2 = sqrt(FDOverlapReImvsReIm(injectedsignalReIm->LLOSignal, injectedsignalReIm->LLOSignal, injectedsignalReIm->noisevaluesLLO));
		SNR3 = sqrt(FDOverlapReImvsReIm(injectedsignalReIm->VIRGOSignal, injectedsignalReIm->VIRGOSignal, injectedsignalReIm->noisevaluesVIRGO));
		SNR123 = sqrt(SNR1*SNR1 + SNR2*SNR2 + SNR3*SNR3);
	}

	/* Rescale distance to match SNR */
  if (!isnan(priorParams->snr_target)) {
    if (myid == 0) printf("Rescaling the distance to obtain a network SNR of %g\n", priorParams->snr_target);
    injectedparams->distance *= SNR123 / priorParams->snr_target;
    if (myid == 0) printf("New distance = %g Mpc\n", injectedparams->distance);
    if(priorParams->rescale_distprior) {
      priorParams->dist_min *= SNR123 / priorParams->snr_target;
      priorParams->dist_max *= SNR123 / priorParams->snr_target;
      if (myid == 0) printf("Distance prior (dist_min, dist_max) = (%g, %g) Mpc\n", priorParams->dist_min, priorParams->dist_max);
    }
    if (myid == 0) print_rescaleddist_to_file_LLV(injectedparams, globalparams, priorParams, &runParams);
    if(globalparams->tagint==0) {
      LLVGenerateInjectionCAmpPhase(injectedparams, injectedsignalCAmpPhase);
      SNR123 = sqrt(injectedsignalCAmpPhase->LLVss);
    }
    else if(globalparams->tagint==1) {
      LLVGenerateInjectionReIm(injectedparams, globalparams->minf, globalparams->maxf, globalparams->nbptsoverlap, 1, injectedsignalReIm); /* tagsampling fixed to 1, i.e. logarithmic sampling - could be made another global parameter */
      SNR1 = sqrt(FDOverlapReImvsReIm(injectedsignalReIm->LHOSignal, injectedsignalReIm->LHOSignal, injectedsignalReIm->noisevaluesLHO));
      SNR2 = sqrt(FDOverlapReImvsReIm(injectedsignalReIm->LLOSignal, injectedsignalReIm->LLOSignal, injectedsignalReIm->noisevaluesLLO));
      SNR3 = sqrt(FDOverlapReImvsReIm(injectedsignalReIm->VIRGOSignal, injectedsignalReIm->VIRGOSignal, injectedsignalReIm->noisevaluesVIRGO));
      SNR123 = sqrt(SNR1*SNR1 + SNR2*SNR2 + SNR3*SNR3);
    }
  }

  /* print SNRs */
  if (myid == 0) {
    printf("SNR Network: %g\n", SNR123);
  }

  /* Calculate logL of injection */
	double logZinj = 0;
	if(globalparams->tagint==0) {
    logZinj = CalculateLogLCAmpPhase(injectedparams, injectedsignalCAmpPhase);
  }
  else if(globalparams->tagint==1) {
    logZinj = CalculateLogLReIm(injectedparams, injectedsignalReIm);
  }
  if (myid == 0) printf("logZinj = %lf\n", logZinj);

	/* Set the context pointer */
  void *context = NULL;
  if(globalparams->tagint==0) {
    context = injectedsignalCAmpPhase;
  }
  else if(globalparams->tagint==1) {
    context = injectedsignalReIm;
  }

  int nPar = 9;	  /* Total no. of parameters including free & derived parameters */
  int ndim = 9;   /* No. of free parameters - to be changed later if some parameters are fixed */

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

	/* Check for fixed parameters, and build the map from the free cube parameters to the orignal 9 parameters */
	/* Order of the 9 original parameters (fixed): m1, m2, tRef, dist, phase, inc, ra, dec, pol */
	/* Order of the 9 cube parameters (modified for clustering): ra, dec, tRef, phase, pol, inc, dist, m1, m2 */
	int mapcubetophys[9] = {6, 7, 2, 4, 8, 5, 3, 0, 1};
	int freecubeparams[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
	if (!isnan(priorParams->fix_ra))     { ndim--; freecubeparams[0] = 0; }
	if (!isnan(priorParams->fix_dec))    { ndim--; freecubeparams[1] = 0; }
	if (!isnan(priorParams->fix_time))   { ndim--; freecubeparams[2] = 0; }
	if (!isnan(priorParams->fix_phase))  { ndim--; freecubeparams[3] = 0; }
	if (!isnan(priorParams->fix_pol))    { ndim--; freecubeparams[4] = 0; }
	if (!isnan(priorParams->fix_inc))    { ndim--; freecubeparams[5] = 0; }
	if (!isnan(priorParams->fix_dist))   { ndim--; freecubeparams[6] = 0; }
	if (!isnan(priorParams->fix_m1))     { ndim--; freecubeparams[7] = 0; }
	if (!isnan(priorParams->fix_m2))     { ndim--; freecubeparams[8] = 0; }
	int* freeparamsmap = malloc(ndim*sizeof(int));
	int counter = 0;
	for(int i=0; i<ndim; i++) {
		while(freecubeparams[counter]==0) counter++;
		freeparamsmap[i] = mapcubetophys[counter];
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

		double logL = 0.;
		if(globalparams->tagint==0) {
			logL = CalculateLogLCAmpPhase(&templateparams, injectedsignalCAmpPhase);
		}
		else if(globalparams->tagint==1) {
			logL = CalculateLogLReIm(&templateparams, injectedsignalReIm);
		}
    if (myid == 0) printf("logL = %lf\n", logL);

    free(injectedparams);
    free(priorParams);

#ifdef PARALLEL
    MPI_Finalize();
#endif

    exit(0);
  }

	/* If the seed option is activated, seed the initial population of live points with the injection */
	if(myid==0 && runParams.seed) {
		char* pathresume = malloc(strlen(runParams.outroot)+64);
		char* pathev = malloc(strlen(runParams.outroot)+64);
		char* pathlivepoints = malloc(strlen(runParams.outroot)+64);
		char* pathphyslivepoints = malloc(strlen(runParams.outroot)+64);
		sprintf(pathresume, "%s%s", runParams.outroot, "resume.dat");
		sprintf(pathev, "%s%s", runParams.outroot, "ev.dat");
		sprintf(pathlivepoints, "%s%s", runParams.outroot, "live.points");
		sprintf(pathphyslivepoints, "%s%s", runParams.outroot, "phys_live.points");
		int tagseed = 0; /* Tag deciding whether or not to create the seed files */
		FILE* fresume = NULL;
		FILE* fev = NULL;
		FILE* flivepoints = NULL;
		FILE* fphyslivepoints = NULL;
		if(runParams.resume) { /* Check whether files already exist */
			fresume = fopen(pathresume, "r");
			fev = fopen(pathev, "r");
			flivepoints = fopen(pathlivepoints, "r");
			fphyslivepoints = fopen(pathphyslivepoints, "r");
			if(!fresume && !fev && !flivepoints && !fphyslivepoints) { /* If none of the files already exist - create seed */
	tagseed = 1;
			}
			else if (fresume && fev && flivepoints && fphyslivepoints) { /* If all files already exist, ignore --seed and do nothing - allows to resume run */
	tagseed = 0;
	fclose(fresume);
	fclose(fev);
	fclose(flivepoints);
	fclose(fphyslivepoints);
			}
			else {
	printf("Error: when seeding, some files already exist but not all of them.");
	exit(1);
			}
		}
		else tagseed = 1; /* If the resume option is false, create the seed anyway (possibly overwriting) */

		if(tagseed) { /* Create the seeding files */
			printf("Seeding the inference with one point at the injection.\n");
			fresume = fopen(pathresume, "w");
			fev = fopen(pathev, "w");
			flivepoints = fopen(pathlivepoints, "w");
			fphyslivepoints = fopen(pathphyslivepoints, "w");
			/* Resume and ev files */
			fprintf(fresume, " T\n");
			fprintf(fev, "");
			/* Phys live points */
			fprintf(fphyslivepoints, "    %.18E", injectedparams->m1);
			fprintf(fphyslivepoints, "    %.18E", injectedparams->m2);
			fprintf(fphyslivepoints, "    %.18E", 0.); /* For templates, tRef is defined relatively to the injected value */
			fprintf(fphyslivepoints, "    %.18E", injectedparams->distance);
			fprintf(fphyslivepoints, "    %.18E", injectedparams->phiRef);
			fprintf(fphyslivepoints, "    %.18E", injectedparams->inclination);
			fprintf(fphyslivepoints, "    %.18E", injectedparams->ra);
			fprintf(fphyslivepoints, "    %.18E", injectedparams->dec);
			fprintf(fphyslivepoints, "    %.18E", injectedparams->polarization);
			fprintf(fphyslivepoints, "   %.18E", logZinj);
			fprintf(fphyslivepoints, "   %d\n", 1); /* We impose that the injection belongs to mode no. 1 */
			/* Live points - convert to values in the cube */
			double* cubevalues = malloc(ndim*sizeof(double));
			getcubeparams(cubevalues, ndim, injectedparams, freeparamsmap);
			for(int i=0; i<ndim; i++) fprintf(flivepoints, "    %.18E", cubevalues[i]);
			fprintf(flivepoints, "   %.18E\n", logZinj);
			free(cubevalues);

			/* Also set the resume option to true */
			runParams.resume = 1;

			fclose(fresume);
			fclose(fev);
			fclose(flivepoints);
			fclose(fphyslivepoints);
		}

		free(pathresume);
		free(pathev);
		free(pathlivepoints);
		free(pathphyslivepoints);
	}

	/********** End of addendum ****************/

	int i;

	/* set the MultiNest sampling parameters */

	int mmodal = runParams.mmodal;					/* do mode separation? */

	int ceff = 0;					/* run in constant efficiency mode? */

	int nlive = runParams.nlive;				/* number of live points */

	double efr = runParams.eff;				/* set the required efficiency */

	double tol = runParams.tol;				/* tol, defines the stopping criteria */

	/* int ndim = 9; */	        /* dimensionality (no. of free parameters) */
	/* int nPar = 9; */					/* total no. of parameters including free & derived parameters */

	int nClsPar = runParams.nclspar;				/* no. of parameters to do mode separation on */

	int updInt = 50;				/* after how many iterations feedback is required & the output files should be updated */
							/* note: posterior files are updated & dumper routine is called after every updInt*10 iterations */

	double Ztol = runParams.ztol;				/* all the modes with logZ < Ztol are ignored */

	int maxModes = runParams.maxcls;				/* expected max no. of modes (used only for memory allocation) */

	int pWrap[ndim];				/* which parameters to have periodic boundary conditions? */
	for(i = 0; i < ndim; i++) { /* If non-default limiting values have been set for ra,, phase, pol, do not treat them as periodic */
	  if(freeparamsmap[i]==4 && priorParams->phase_min == 0. && priorParams->phase_max == 2.*PI) pWrap[i] = 1;
	  else if(freeparamsmap[i]==6 && priorParams->ra_min == 0. && priorParams->ra_max == 2.*PI) pWrap[i] = 1;
	  else if(freeparamsmap[i]==8 && priorParams->pol_min == 0. && priorParams->pol_max == 2.*PI) pWrap[i] = 1;
	  else pWrap[i] = 0;
	}

	strcpy(root, runParams.outroot);		/* root for output files */
	strcpy(networkinputs, runParams.netfile);	/* file with input parameters for network training */

	int seed = -1;					/* random no. generator seed, if < 0 then take the seed from system clock */

	int fb = 1;					/* need feedback on standard output ? */

	resume = runParams.resume;			/* resume from a previous job? */

	int outfile = 1;				/* write output files? */

	int initMPI = 0;				/* initialize MPI routines?, relevant only if compiling with MPI */
							/* set it to F if you want your main program to handle MPI initialization */

	logZero = -1E90;				/* points with loglike < logZero will be ignored by MultiNest */

	int maxiter = runParams.maxiter;				/* max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it */
							/* has done max no. of iterations or convergence criterion (defined through tol) has been satisfied */

	/* void *context = 0; */			/* not required by MultiNest, any additional information user wants to pass */

	doBAMBI = runParams.bambi;					/* BAMBI? */

	useNN = 0;

	/* calling MultiNest */

	BAMBIrun(mmodal, ceff, nlive, tol, efr, ndim, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, LogLikeFctn, dumper, BAMBIfctn, context);

	//
	printf("%p, %p\n", injectedparams, priorParams);

  free(injectedparams);
  free(priorParams);

#ifdef PARALLEL
 	MPI_Finalize();
#endif
}

/***********************************************************************************************************************/
