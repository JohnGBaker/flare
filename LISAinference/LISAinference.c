#include "LISAinference_common.h"
#include "LISAinference.h"

/******************************************** getphysparams routine ****************************************************/

/* Order of the parameters (fixed): m1, m2, tRef, dist, phase, inc, lambda, beta, pol */
void getphysparams(double *Cube, int *ndim) /* Note: ndim not used here */
{
  int i = 0;
  double m1=0., m2=0., tRef=0., dist=0., phase=0., inc=0., lambda=0., beta=0., pol=0., Mchirp=0., eta=0.;

  /* Note: we use here the order for the cube parameters */
  /* Order of the 9 original parameters (fixed): m1, m2, tRef, dist, phase, inc, lambda, beta, pol */
  /* Order of the 9 cube parameters (modified for clustering): lambda, beta, tRef, phase, pol, inc, dist, m1, m2 */

  /* Sky location (lambda then beta) */
  if (isnan(priorParams->fix_lambda)) {
    lambda = CubeToFlatPrior(Cube[i++], priorParams->lambda_min, priorParams->lambda_max);
  } else {
    lambda = priorParams->fix_lambda;
  }
  if (isnan(priorParams->fix_beta)) {
    beta = CubeToCosPrior(Cube[i++], priorParams->beta_min, priorParams->beta_max);
  } else {
    beta = priorParams->fix_beta;
  }

  /* Time */
  /* If sampling in tL, convert tinj SSB to tLinj */
  if (isnan(priorParams->fix_time)) {
    if(priorParams->sampletimeparam==tSSB) {
      tRef = CubeToFlatPrior(Cube[i++], injectedparams->tRef - priorParams->deltaT, injectedparams->tRef + priorParams->deltaT);
    }
    else if(priorParams->sampletimeparam==tL) { /* Here tRef has the meaning of tL */
      double injectedtL = tLfromtSSB(globalparams->variant, injectedparams->tRef, injectedparams->lambda, injectedparams->beta);
      tRef = CubeToFlatPrior(Cube[i++], injectedtL - priorParams->deltaT, injectedtL + priorParams->deltaT);
    }
  } else { /* fix_time, if defined, has the sense of a SSB time */
    if(priorParams->sampletimeparam==tSSB) {
      tRef = priorParams->fix_time;
    }
    else if(priorParams->sampletimeparam==tL) { /* Set tRef (meaning tL) to the injected tL (with injected sky position) */
      tRef = tLfromtSSB(globalparams->variant, priorParams->fix_time, injectedparams->lambda, injectedparams->beta);
    }
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

  /* Mass - two priors allowed: flat or log-flat */
  double (*mass_prior) (double r, double x1, double x2);
  mass_prior = &CubeToFlatPrior;
  if(priorParams->logflat_massprior) mass_prior = &CubeToLogFlatPrior;

  /* Component masses */
  /* Branch between sampling in m1/m2 or Mchirp/eta */
  if(priorParams->samplemassparams==m1m2) {
    if (isnan(priorParams->fix_m1)) {
      if (isnan(priorParams->fix_m2)) {
        m1 = mass_prior(Cube[i++], priorParams->comp_min, priorParams->comp_max);
        m2 = mass_prior(Cube[i++], priorParams->comp_min, priorParams->comp_max);
    } else {
        m2 = priorParams->fix_m2;
        m1 = mass_prior(Cube[i++], priorParams->comp_min, priorParams->comp_max);
      }
    } else {
      m1 = priorParams->fix_m1;
      if (isnan(priorParams->fix_m2)) {
        m2 = mass_prior(Cube[i++], priorParams->comp_min, priorParams->comp_max);
      } else {
        m2 = priorParams->fix_m2;
      }
    }
  }
  else if(priorParams->samplemassparams==Mchirpeta) {
    if (isnan(priorParams->fix_Mchirp)) {
      if (isnan(priorParams->fix_eta)) {
        Mchirp = mass_prior(Cube[i++], priorParams->Mchirp_min, priorParams->Mchirp_max);
        eta = CubeToFlatPrior(Cube[i++], priorParams->eta_min, priorParams->eta_max);
      } else {
        eta = priorParams->fix_eta;
        Mchirp = mass_prior(Cube[i++], priorParams->Mchirp_min, priorParams->Mchirp_max);
      }
    } else {
      Mchirp = priorParams->fix_Mchirp;
      if (isnan(priorParams->fix_eta)) {
        eta = CubeToFlatPrior(Cube[i++], priorParams->eta_min, priorParams->eta_max);
      } else {
        eta = priorParams->fix_eta;
      }
    }
    /* Convert Mchirp/eta to m1/m2 */
    m1 = m1ofMchirpeta(Mchirp, eta);
    m2 = m2ofMchirpeta(Mchirp, eta);
  }

  /* Convert time - if sampling in tL, compute tSSB from tL using (approximate but 3e-6s accurate) inverse relation */
  /* The tRef we output has always the meaning of a SSB time */
  if(priorParams->sampletimeparam==tL) {
    tRef = tSSBfromtL(globalparams->variant, tRef, lambda, beta);
  }

  /* Note: here we output physical values in the cube (overwriting), and we keep the original order for physical parameters */
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

void getcubeparams(double* Cube, int ndim, LISAParams* params, int* freeparamsmap)
{
  int i = 0;
  double m1 = params->m1;
  double m2 = params->m2;
  double tRef = params->tRef;
  double dist = params->distance;
  double phase = params->phiRef;
  double inc = params->inclination;
  double lambda = params->lambda;
  double beta = params->beta;
  double pol = params->polarization;

  /* This is used only if sampling uses Mchirpeta */
  double Mchirp = Mchirpofm1m2(m1, m2);
  double eta = etaofm1m2(m1, m2);

  /* Note: freeparamsmap has indices in the order of the (free) cube parameters, and values in the indices of physical parameters */
  for(int i=0; i<ndim; i++) {
    if(freeparamsmap[i]==0) {
      if(priorParams->samplemassparams==m1m2) {
        if(!(priorParams->logflat_massprior)) Cube[i] = FlatPriorToCube(m1, priorParams->comp_min, priorParams->comp_max);
        else Cube[i] = LogFlatPriorToCube(m1, priorParams->comp_min, priorParams->comp_max);
      }
      else if(priorParams->samplemassparams==Mchirpeta) {
        if(!(priorParams->logflat_massprior)) Cube[i] = FlatPriorToCube(Mchirp, priorParams->Mchirp_min, priorParams->Mchirp_max);
        else Cube[i] = LogFlatPriorToCube(Mchirp, priorParams->Mchirp_min, priorParams->Mchirp_max);
      }
    }
    if(freeparamsmap[i]==1) {
      if(priorParams->samplemassparams==m1m2) {
        if(!(priorParams->logflat_massprior)) Cube[i] = FlatPriorToCube(m2, priorParams->comp_min, priorParams->comp_max);
        else Cube[i] = LogFlatPriorToCube(m2, priorParams->comp_min, priorParams->comp_max);
      }
      else if(priorParams->samplemassparams==Mchirpeta) {
        Cube[i] = FlatPriorToCube(eta, priorParams->eta_min, priorParams->eta_max);
      }
    }
    if(freeparamsmap[i]==2) Cube[i] = FlatPriorToCube(tRef, injectedparams->tRef - priorParams->deltaT, injectedparams->tRef + priorParams->deltaT);
    if(freeparamsmap[i]==3) {
      if(priorParams->flat_distprior) Cube[i] = FlatPriorToCube(dist, priorParams->dist_min, priorParams->dist_max);
      else Cube[i] = PowerPriorToCube(2., dist, priorParams->dist_min, priorParams->dist_max);
    }
    if(freeparamsmap[i]==4) Cube[i] = FlatPriorToCube(phase, priorParams->phase_min, priorParams->phase_max);
    if(freeparamsmap[i]==5) Cube[i] = SinPriorToCube(inc, priorParams->inc_min, priorParams->inc_max);
    if(freeparamsmap[i]==6) Cube[i] = FlatPriorToCube(lambda, priorParams->lambda_min, priorParams->lambda_max);
    if(freeparamsmap[i]==7) Cube[i] = CosPriorToCube(beta, priorParams->beta_min, priorParams->beta_max);
    if(freeparamsmap[i]==8) Cube[i] = FlatPriorToCube(pol, priorParams->pol_min, priorParams->pol_max);
  }
}

/******************************************** getallparams routine ****************************************************/

void getallparams(double *Cube, int *ndim)
{
	getphysparams(Cube,ndim);
}

/******************************************** loglikelihood routine ****************************************************/

 /* Now an example, sample an egg box likelihood */

 /* Input arguments */
 /* ndim 						= dimensionality (total number of free parameters) of the problem */
 /* npars 						= total number of free plus derived parameters */

 /* Input/Output arguments */
 /* Cube[npars] 						= on entry has the ndim parameters in unit-hypercube */
 /* 	 						on exit, the physical parameters plus copy any derived parameters you want to store with the free parameters */

 /* Output arguments */
 /* lnew 						= loglikelihood */

/* Note: context must point to the LISASignal structure representing the injected signals */
void getLogLike(double *Cube, int *ndim, int *npars, double *lnew, void *context)
{
  /* Convert Cube to physical parameters and check prior boundary */
  /* Branch between two different checks for m1/m2 or Mchirp/eta */
  getallparams(Cube, ndim);
  if (priorParams->samplemassparams==m1m2 && PriorBoundaryCheckm1m2(priorParams, Cube)) {
    *lnew = -DBL_MAX;
    return;
  }
  if (priorParams->samplemassparams==Mchirpeta && PriorBoundaryCheckMchirpeta(priorParams, Cube)) {
    *lnew = -DBL_MAX;
    return;
  }

  /* If zerolikelihood option is activated, simply set result to 0 for all parameters choices within prior */
  /* Useful for testing purposes, to just sample from the prior */
  if(globalparams->zerolikelihood) {
    *lnew = 0;
    return;
  }

  /* Cube values here have the order of the 9 physical parameters (fixed): m1, m2, tRef, dist, phase, inc, lambda, beta, pol */
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
  if((globalparams->tagint==0) && (!globalparams->tagsimplelikelihood)) {
    LISAInjectionCAmpPhase* injection = ((LISAInjectionCAmpPhase*) context);

    //TESTING
    //clock_t tbeg, tend;
    //tbeg = clock();
    *lnew = CalculateLogLCAmpPhase(&templateparams, injection);
    //tend = clock();
    //printf("time Likelihood: %g\n", (double) (tend-tbeg)/CLOCKS_PER_SEC);
    //
  }
  else if((globalparams->tagint==1) && (!globalparams->tagsimplelikelihood)) {
    LISAInjectionReIm* injection = ((LISAInjectionReIm*) context);

    //TESTING
    //clock_t tbeg, tend;
    //tbeg = clock();
    *lnew = CalculateLogLReIm(&templateparams, injection);
    //tend = clock();
    //printf("time Likelihood: %g\n", (double) (tend-tbeg)/CLOCKS_PER_SEC);
    //
  }
  else if(globalparams->tagsimplelikelihood) {
    SimpleLikelihoodPrecomputedValues* injection = ((SimpleLikelihoodPrecomputedValues*) context);
    *lnew = CalculateLogLSimpleLikelihood(injection, &templateparams);
  }
}


/************************************************* dumper routine ******************************************************/

 /* The dumper routine will be called every updInt*10 iterations */
 /* MultiNest doesn not need to the user to do anything. User can use the arguments in whichever way he/she wants */


 /* Arguments: */

 /* nSamples 						= total number of samples in posterior distribution */
 /* nlive 						= total number of live points */
 /* nPar 						= total number of parameters (free + derived) */
 /* physLive[1][nlive * (nPar + 1)] 			= 2D array containing the last set of live points (physical parameters plus derived parameters) along with their loglikelihood values */
 /* posterior[1][nSamples * (nPar + 2)] 			= posterior distribution containing nSamples points. Each sample has nPar parameters (physical + derived) along with the their loglike value & posterior probability */
 /* paramConstr[1][4*nPar]: */
 /* paramConstr[0][0] to paramConstr[0][nPar - 1] 	= mean values of the parameters */
 /* paramConstr[0][nPar] to paramConstr[0][2*nPar - 1] 	= standard deviation of the parameters */
 /* paramConstr[0][nPar*2] to paramConstr[0][3*nPar - 1] = best-fit (maxlike) parameters */
 /* paramConstr[0][nPar*4] to paramConstr[0][4*nPar - 1] = MAP (maximum-a-posteriori) parameters */
 /* maxLogLike						= maximum loglikelihood value */
 /* logZ							= log evidence value */
 /* logZerr						= error on log evidence value */
 /* context						void pointer, any additional information */

void dumper(int *nSamples, int *nlive, int *nPar, double **physLive, double **posterior, double **paramConstr, double *maxLogLike, double *logZ, double *logZerr, void *context)
{
  /* convert the 2D Fortran arrays to C arrays */


  /* the posterior distribution */
  /* postdist will have nPar parameters in the first nPar columns & loglike value & the posterior probability in the last two columns */
  /*
  int i, j;

  double postdist[*nSamples][*nPar + 2];
  for( i = 0; i < *nPar + 2; i++ )
    for( j = 0; j < *nSamples; j++ )
      postdist[j][i] = posterior[0][i * (*nSamples) + j];
  */


  /* last set of live points */
  /* pLivePts will have nPar parameters in the first nPar columns & loglike value in the last column */
  /*
  double pLivePts[*nlive][*nPar + 1];
  for( i = 0; i < *nPar + 1; i++ )
    for( j = 0; j < *nlive; j++ )
      pLivePts[j][i] = physLive[0][i * (*nlive) + j];
  */
}

/***********************************************************************************************************************/




/************************************************** Main program *******************************************************/



int main(int argc, char *argv[])
{
  int myid = 0;
  noMPI = 0;
#ifdef PARALLEL
 	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
#endif

	LISARunParams runParams={};
	int ndim=0,nPar=0;
	int *freeparamsmap = NULL;
	void *context = NULL;
	double logZtrue;
	addendum(argc, argv, &runParams, &ndim, &nPar, &freeparamsmap, &context, &logZtrue);

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
	    fprintf(fphyslivepoints, "    %.18E", injectedparams->tRef);
	    fprintf(fphyslivepoints, "    %.18E", injectedparams->distance);
	    fprintf(fphyslivepoints, "    %.18E", injectedparams->phiRef);
	    fprintf(fphyslivepoints, "    %.18E", injectedparams->inclination);
	    fprintf(fphyslivepoints, "    %.18E", injectedparams->lambda);
	    fprintf(fphyslivepoints, "    %.18E", injectedparams->beta);
	    fprintf(fphyslivepoints, "    %.18E", injectedparams->polarization);
	    fprintf(fphyslivepoints, "   %.18E", logZtrue);
	    fprintf(fphyslivepoints, "   %d\n", 1); /* We impose that the injection belongs to mode no. 1 */
	    /* Live points - convert to values in the cube */
	    double* cubevalues = malloc(ndim*sizeof(double));
	    getcubeparams(cubevalues, ndim, injectedparams, freeparamsmap);

	    for(int i=0; i<ndim; i++) fprintf(flivepoints, "    %.18E", cubevalues[i]);
	    fprintf(flivepoints, "   %.18E\n", logZtrue);
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

	/********** Test ****************/
	/* double l; */
	/* l = CalculateLogLReIm(injectedparams, injectedsignalReIm); */
	/* printf("LogLikelihood: %g\n", l); */
	/* free(injectedparams); */
	/* LISAInjectionReIm_Cleanup(injectedsignalReIm); */
	/* exit(0); */
	/********** End of test ****************/

	int i;

	/* set the MultiNest sampling parameters */

	int mmodal = runParams.mmodal;			/* do mode separation? */

	int ceff = 0;					/* run in constant efficiency mode? */

	int nlive = runParams.nlive;			/* number of live points */

	double efr = runParams.eff;			/* set the required efficiency */

	double tol = runParams.tol;			/* tol, defines the stopping criteria */

	/* int ndim = 9; */				        /*dimensionality (no. of free parameters) */
	/* int nPar = 9; */					/* total no. of parameters including free & derived parameters */

	int nClsPar = runParams.nclspar;            	/* no. of parameters to do mode separation on */

	int updInt = 50;				/* after how many iterations feedback is required & the output files should be updated */
							/* note: posterior files are updated & dumper routine is called after every updInt*10 iterations */

	double Ztol = runParams.ztol;	       		/* all the modes with logZ < Ztol are ignored */

	int maxModes = runParams.maxcls;		/* expected max no. of modes (used only for memory allocation) */

	int pWrap[ndim];				/* which parameters to have periodic boundary conditions? */
	for(i = 0; i < ndim; i++) { /* If non-default limiting values have been set for lambda, phase, pol, do not treat them as periodic */
	  if(freeparamsmap[i]==4 && priorParams->phase_min == 0. && priorParams->phase_max == 2.*PI) pWrap[i] = 1;
	  else if(freeparamsmap[i]==6 && priorParams->lambda_min == 0. && priorParams->lambda_max == 2.*PI) pWrap[i] = 1;
	  else if(freeparamsmap[i]==8 && priorParams->pol_min == 0. && priorParams->pol_max == PI) pWrap[i] = 1;
	  else pWrap[i] = 0;
	}

	strcpy(root, runParams.outroot);		/* root for output files */
	strcpy(networkinputs, runParams.netfile);	/* file with input parameters for network training */

	int seed = -1;					/* random no. generator seed, if < 0 then take the seed from system clock */

	int fb = 1;					/* need feedback on standard output? */

	resume = runParams.resume;			/* resume from a previous job? */

	int outfile = 1;				/* write output files? */

	int initMPI = 0;				/* initialize MPI routines?, relevant only if compiling with MPI */
							/* set it to F if you want your main program to handle MPI initialization */

	logZero = -1E90;				/* points with loglike < logZero will be ignored by MultiNest */

	int maxiter = runParams.maxiter;				/* max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it */
							/* has done max no. of iterations or convergence criterion (defined through tol) has been satisfied */

	/* void *context = 0; */				/* not required by MultiNest, any additional information user wants to pass */

	doBAMBI = runParams.bambi;			/* BAMBI? */

	useNN = 0;

	/* calling MultiNest */

	BAMBIrun(mmodal, ceff, nlive, tol, efr, ndim, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, LogLikeFctn, dumper, BAMBIfctn, context);

  free(injectedparams);
  free(priorParams);

#ifdef PARALLEL
 	MPI_Finalize();
#endif
}
