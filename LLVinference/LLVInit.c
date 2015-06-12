#include "LLVInit.h"

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

double LALInferenceCubeToGaussianPrior(double r, double mean, double sigma)
{
    return gsl_cdf_gaussian_Pinv(r,sigma) + mean;
}

double LALInferenceCubeToSinPrior(double r, double x1, double x2)
{
    return acos((1.0-r)*cos(x1)+r*cos(x2));
}

double LALInferenceCubeToCosPrior(double r, double x1, double x2)
{
    return asin((1.0-r)*sin(x1)+r*sin(x2));
}

LLVRunParams* LLInitializeRunParams(ssize_t argc, char **argv)
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
