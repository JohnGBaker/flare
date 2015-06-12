#ifndef __LLVINIT_H__
#define __LLVINIT_H__ 1

#include <math.h>
#include <gsl/gsl_cdf.h>

typedef struct tagLLVPrior {
	double deltaT;             /* width of time prior centered on injected value (s) (default 0.1) */
	double comp_min;           /* minimum component mass (solar masses) (default 4) */
	double comp_max;           /* maximum component mass (solar masses) (default 50) */
	double mtot_min;           /* minimum total mass (solar masses) (default 8) */
	double mtot_max;           /* maximum total mass (solar masses) (default 100) */
	double qmax;               /* maximum asymmetric mass ratio (>=1) (default 12) */
	double dist_min;           /* minimum distance of source (pc) (default 1e6) */
	double dist_max;           /* maximum distance of source (pc) (default 10*1e9) */
	bool inprior;              /* bool to indicate if point is in prior */
} LLVPrior;

// initializes the prior boundaries
LLVPrior* LLVInitializePrior(ssize_t argc, char **argv);

// checks prior boundaires
int PriorBoundaryCheck(LLVPrior *prior, LLVParams *template);

// Prior functions from Cube to physical parameters
// x1 is min, x2 is max when specified
// r is Cube value
double CubeToFlatPrior(double r, double x1, double x2);
double CubeToLogFlatPrior(double r, double x1, double x2);
double CubeToPowerPrior(double p, double r, double x1, double x2);
double LALInferenceCubeToGaussianPrior(double r, double mean, double sigma);
double LALInferenceCubeToSinPrior(double r, double x1, double x2);
double LALInferenceCubeToCosPrior(double r, double x1, double x2);

#endif
