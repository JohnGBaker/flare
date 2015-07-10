#include "LLVutils.h"

/************ Global Parameters ************/

LLVParams* injectedparams = NULL;
LLVParams* templateparams = NULL;
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
LLVParams* parse_args_LLV(ssize_t argc, char **argv) {
    char help[] = "";

    ssize_t i;
    LLVParams* params;
    params = (LLVParams*) malloc(sizeof(LLVParams));
    memset(params, 0, sizeof(LLVParams));

    /* Set default values to the arguments */
    params->tRef = 0.;
    params->phiRef = 0.;
    params->m1 = 20.;
    params->m2 = 10.;
    params->distance = 100.;
    params->ra = 0.;
    params->dec = 0.;
    params->inclination = PI/3.;
    params->polarization = 0.;
    params->fRef = 0.;
    params->nbmode = 5;

    /* Consume command line */
    for (i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--tRef") == 0) {
            params->tRef = atof(argv[++i]);
        } else if (strcmp(argv[i], "--phiRef") == 0) {
            params->phiRef = atof(argv[++i]);
        } else if (strcmp(argv[i], "--m1") == 0) {
            params->m1 = atof(argv[++i]);
        } else if (strcmp(argv[i], "--m2") == 0) {
            params->m2 = atof(argv[++i]);
        } else if (strcmp(argv[i], "--distance") == 0) {
            params->distance = atof(argv[++i]);
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
            //printf("Error: invalid option: %s\n", argv[i]);
            //goto fail;
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
  int ret;
  ListmodesCAmpPhaseFrequencySeries* listROM = NULL;
  ListmodesCAmpPhaseFrequencySeries* listLHO = NULL;
  ListmodesCAmpPhaseFrequencySeries* listLLO = NULL;
  ListmodesCAmpPhaseFrequencySeries* listVIRGO = NULL;

  /* Checking that the global injectedparams has been set up */
  if (!injectedparams) {
    printf("Error: when calling LLVGenerateSignal, injectedparams points to NULL.\n");
    exit(1);
  }
  /* Should add more error checking ? */
  /* Generate the waveform with the ROM */
  /* Note: SimEOBNRv2HMROM accepts masses and distances in SI units, whereas LLV params is in solar masses and Mpc */
  ret = SimEOBNRv2HMROM(&listROM, params->nbmode, params->tRef - injectedparams->tRef, params->phiRef, params->fRef, (params->m1)*MSUN_SI, (params->m2)*MSUN_SI, (params->distance)*1e6*PC_SI);

  /* If the ROM waveform generation failed (e.g. parameters were out of bounds) return FAILURE */
  if(ret==FAILURE) return FAILURE;

  /* Process the waveform through the LLV response */
  LLVSimFDResponse(&listROM, &listLHO, params->tRef, params->ra, params->dec, params->inclination, params->polarization, LHO);
  LLVSimFDResponse(&listROM, &listLLO, params->tRef, params->ra, params->dec, params->inclination, params->polarization, LLO);
  LLVSimFDResponse(&listROM, &listVIRGO, params->tRef, params->ra, params->dec, params->inclination, params->polarization, VIRGO);

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

LLVPrior* LLVInitializePrior(ssize_t argc, char **argv)
{
    ssize_t i;
    LLVPrior* prior;
    prior = (LLVPrior*) malloc(sizeof(LLVPrior));

    // set defaults
    prior->deltaT = 0.1;
    prior->comp_min = 4.0;
    prior->comp_max = 50.0;
    prior->mtot_min = 8.0;
    prior->mtot_max = 100.0;
    prior->qmax = 12.0;
    prior->dist_min = 1.0e6;
    prior->dist_max = 1.0e10;

    /* Consume command line */
    for (i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--deltaT") == 0) {
            prior->deltaT = atof(argv[++i]);
        } else if (strcmp(argv[i], "--comp-min") == 0) {
            prior->comp_min = atof(argv[++i]);
        } else if (strcmp(argv[i], "--comp-max") == 0) {
            prior->comp_max = atof(argv[++i]);
        } else if (strcmp(argv[i], "--mtot-min") == 0) {
            prior->mtot_min = atof(argv[++i]);
        } else if (strcmp(argv[i], "--mtot-max") == 0) {
            prior->mtot_max = atof(argv[++i]);
        } else if (strcmp(argv[i], "--qmax") == 0) {
            prior->qmax = atof(argv[++i]);
        } else if (strcmp(argv[i], "--dist-min") == 0) {
            prior->dist_min = atof(argv[++i]);
        } else if (strcmp(argv[i], "--dist-max") == 0) {
            prior->dist_max = atof(argv[++i]);
        } else {
            //printf("Error: invalid option: %s\n", argv[i]);
            //goto fail;
        }
    }

    return prior;

    /*fail:
    free(prior);
    exit(1);*/
}

int PriorBoundaryCheck(LLVPrior *prior, double *Cube)
{
	if (Cube[0] < prior->comp_min || Cube[0] > prior->comp_max ||
	 	Cube[1] < prior->comp_min || Cube[1] > prior->comp_max)
	 	return 1;

	if (Cube[0] + Cube[1] < prior->mtot_min || Cube[0] + Cube[1] > prior->mtot_max)
		return 1;

	if (Cube[0] < Cube[1] || Cube[0] / Cube[1] > prior->qmax)
		return 1;

	return 0;
}

double CubeToFlatPrior(double r, double x1, double x2)
{
    return x1 + r * ( x2 - x1 );
}

double CubeToLogFlatPrior(double r, double x1, double x2)
{
    double lx1, lx2;
    lx1 = log( x1 );
    lx2 = log( x2 );
    return exp( lx1 + r * ( lx2 - lx1 ) );
}

double CubeToPowerPrior(double p, double r, double x1, double x2)
{
    double pp = p + 1.0;
    return pow(r * pow(x2, pp) + (1.0 - r) * pow(x1, pp), 1.0 / pp);
}

double CubeToGaussianPrior(double r, double mean, double sigma)
{
    return gsl_cdf_gaussian_Pinv(r,sigma) + mean;
}

double CubeToSinPrior(double r, double x1, double x2)
{
    return acos((1.0-r)*cos(x1)+r*cos(x2));
}

double CubeToCosPrior(double r, double x1, double x2)
{
    return asin((1.0-r)*sin(x1)+r*sin(x2));
}

LLVRunParams* LLVInitializeRunParams(ssize_t argc, char **argv)
{
    ssize_t i;
    LLVRunParams* runParams;
    runParams = (LLVRunParams*) malloc(sizeof(LLVRunParams));

    // set defaults
    runParams->eff = 0.1;
    runParams->tol = 0.5;
    runParams->nlive = 1000;
    strcpy(runParams->outroot, "chains/LLVinference_");
    runParams->bambi = 0;
    runParams->resume = 0;
    strcpy(runParams->netfile, "LLVinference.inp");

    /* Consume command line */
    for (i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--eff") == 0) {
            runParams->eff = atof(argv[++i]);
        } else if (strcmp(argv[i], "--tol") == 0) {
            runParams->tol = atof(argv[++i]);
        } else if (strcmp(argv[i], "--nlive") == 0) {
            runParams->nlive = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--bambi") == 0) {
            runParams->bambi = 1;
        } else if (strcmp(argv[i], "--resume") == 0) {
            runParams->resume = 1;
        } else if (strcmp(argv[i], "--outroot") == 0) {
            strcpy(runParams->outroot, argv[++i]);
        } else if (strcmp(argv[i], "--netfile") == 0) {
            strcpy(runParams->netfile, argv[++i]);
        } else {
            //printf("Error: invalid option: %s\n", argv[i]);
            //goto fail;
        }
    }

    return runParams;

    /*fail:
    free(prior);
    exit(1);*/
}
