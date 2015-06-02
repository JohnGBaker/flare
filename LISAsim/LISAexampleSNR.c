/**
 * \author Sylvain Marsat, University of Maryland - NASA GSFC
 *
 * \brief C code for example generation of waveforms with the EOBNRv2HM reduced order model,
 * processing through the Fourier-domain LISA response and SNR calculation.
 *
 */


#define _XOPEN_SOURCE 500

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <unistd.h>
#include <getopt.h>
#include <stdbool.h>
#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_complex.h>

#include "constants.h"
#include "struct.h"
#include "EOBNRv2HMROMstruct.h"
#include "EOBNRv2HMROM.h"
#include "wip.h"
#include "LISAgeometry.h"
#include "LISAFDresponse.h"

/* Parameters for the generation of a ROM waveform (in the form of a list of modes) */
/* All parameters are to be given in SI units! */
typedef struct tagROMParams {
  int nbmode;                /* Number of modes to generate (starting with the 22) - defaults to 1 (22 mode only) */
  double tRef;               /* shift in time with respect to the 22-fit-removed waveform */
  double phiRef;             /* phase at fRef */
  double fRef;               /* reference frequency */
  double m1;                 /* mass of companion 1 */
  double m2;                 /* mass of companion 2 */
  double distance;           /* distance of source */
  char outname[256];         /* file to which output should be written */
} ROMParams;

/* Parameters for the generation of a ROM waveform (in the form of a list of modes) */
/* All parameters are in SI units in the internals */
/* Angle definitions are taken from the Krolak&al paper gr-qc/0401108 */
typedef struct tagLISAParams {
  int nbmode;                /* Number of modes to generate (starting with the 22) - defaults to 1 (22 mode only) */
  double tRef;               /* shift in time with respect to the 22-fit-removed waveform */
  double phiRef;             /* phase at fRef */
  double fRef;               /* reference frequency */
  double m1;                 /* mass of companion 1 */
  double m2;                 /* mass of companion 2 */
  double distance;           /* distance of source */
  double inclination;        /* inclination of L relative to line of sight */
  double lambda;             /* First angle for the position in the sky of the source */
  double beta;               /* Second angle for the position in the sky of the source */
  double psi;                /* Polarization angle */
  char outname[256];         /* file to which output should be written */
} LISAParams;

/* Parse command line and return a newly allocated ROMParams object
 * Masses are input in solar masses and distances in Mpc - converted in SI for the internals */
static ROMParams* parse_args_ROM(ssize_t argc, char **argv) {
    ssize_t i;
    ROMParams* params;
    params = (ROMParams*) malloc(sizeof(ROMParams));
    memset(params, 0, sizeof(ROMParams));

    /* Set default values to the arguments */
    params->nbmode = 1;
    params->tRef = 0.;
    params->phiRef = 0.;
    params->fRef = 0.;
    params->m1 = 1. * 1e6 * MSUN_SI;
    params->m2 = 1. * 1e6 * MSUN_SI;
    params->distance = 1. * 1e9 * PC_SI;

    /* consume command line */
    for (i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--nbmode") == 0) {
            params->nbmode = atof(argv[++i]);
        } else if (strcmp(argv[i], "--tRef") == 0) {
            params->tRef = atof(argv[++i]);
        } else if (strcmp(argv[i], "--phiRef") == 0) {
            params->phiRef = atof(argv[++i]);
        } else if (strcmp(argv[i], "--fRef") == 0) {
            params->fRef = atof(argv[++i]);
        } else if (strcmp(argv[i], "--m1") == 0) {
            params->m1 = atof(argv[++i]) * MSUN_SI;
        } else if (strcmp(argv[i], "--m2") == 0) {
            params->m2 = atof(argv[++i]) * MSUN_SI;
        } else if (strcmp(argv[i], "--distance") == 0) {
            params->distance = atof(argv[++i]) * 1e6 * PC_SI;
        } else if (strcmp(argv[i], "--outname") == 0) {
            strncpy(params->outname, argv[++i], 256);
        } else {
            printf("Error: invalid option: %s\n", argv[i]);
            goto fail;
        }
    }

    return params;

    fail:
    free(params);
    exit(1);
}

/* Parse command line and return a newly allocated LISAParams object
 * Masses are input in solar masses and distances in Mpc - converted in SI for the internals */
static LISAParams* parse_args_LISA(ssize_t argc, char **argv) {
    ssize_t i;
    LISAParams* params;
    params = (LISAParams*) malloc(sizeof(LISAParams));
    memset(params, 0, sizeof(LISAParams));

    /* Set default values to the arguments */
    params->nbmode = 1;
    params->tRef = 0.;
    params->phiRef = 0.;
    params->fRef = 0.;
    params->m1 = 1. * 1e6 * MSUN_SI;
    params->m2 = 1. * 1e6 * MSUN_SI;
    params->distance = 1. * 1e9 * PC_SI;
    params->inclination = 0.;
    params->lambda = 0.;
    params->beta = 0.;
    params->psi = 0.;
    sprintf(params->outname, "");

    /* consume command line */
    for (i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--nbmode") == 0) {
            params->nbmode = atof(argv[++i]);
        } else if (strcmp(argv[i], "--tRef") == 0) {
            params->tRef = atof(argv[++i]);
        } else if (strcmp(argv[i], "--phiRef") == 0) {
            params->phiRef = atof(argv[++i]);
        } else if (strcmp(argv[i], "--fRef") == 0) {
            params->fRef = atof(argv[++i]);
        } else if (strcmp(argv[i], "--m1") == 0) {
            params->m1 = atof(argv[++i]) * MSUN_SI;
        } else if (strcmp(argv[i], "--m2") == 0) {
            params->m2 = atof(argv[++i]) * MSUN_SI;
        } else if (strcmp(argv[i], "--distance") == 0) {
            params->distance = atof(argv[++i]) * 1e6 * PC_SI;
        } else if (strcmp(argv[i], "--inclination") == 0) {
            params->inclination = atof(argv[++i]);
        } else if (strcmp(argv[i], "--lambda") == 0) {
            params->lambda = atof(argv[++i]);
        } else if (strcmp(argv[i], "--beta") == 0) {
            params->beta = atof(argv[++i]);
        } else if (strcmp(argv[i], "--psi") == 0) {
            params->psi = atof(argv[++i]);
        } else if (strcmp(argv[i], "--outname") == 0) {
            strncpy(params->outname, argv[++i], 256);
        } else {
            printf("Error: invalid option: %s\n", argv[i]);
            goto fail;
        }
    }

    return params;

    fail:
    free(params);
    exit(1);
}

/* Function to output to a file an AmpPhaseFrequencySeries; each mode will be output to a separate file using this function */
static int Write_CAmpPhaseFrequencySeries(FILE* f, CAmpPhaseFrequencySeries* freqseries) {
    gsl_vector* freq = freqseries->freq;
    gsl_vector* amp_real = freqseries->amp_real;
    gsl_vector* amp_imag = freqseries->amp_imag;
    gsl_vector* phase = freqseries->phase;

    int len = (int) freq->size;
    /*Here, we could add a check on the length of the gsl_vectors*/

    fprintf(f, "# f amp_re amp_im phase\n");
    for (int i=0; i<len; i++) {
      fprintf(f, "%.16e %.16e %.16e %.16e\n", gsl_vector_get(freq, i), gsl_vector_get(amp_real, i), gsl_vector_get(amp_imag, i), gsl_vector_get(phase, i));
    }

    return 0;
}

/*
 * Swhitenoise
 * a simple noise function
 */
double Swhitenoise(double x){
  return 1;
};

/*
 * main
 */
int main (int argc , char **argv) {
    FILE* f;
    ListmodesCAmpPhaseFrequencySeries* listROM = NULL;
    ListmodesCAmpPhaseFrequencySeries* listA = NULL;
    ListmodesCAmpPhaseFrequencySeries* listE = NULL;
    ListmodesCAmpPhaseFrequencySeries* listT = NULL;
    //ROMParams* params;
    LISAParams* params;

    /* parse commandline */
    //params = parse_args_ROM(argc, argv);
    params = parse_args_LISA(argc, argv);

    /* Generate the waveform with the ROM */
    SimEOBNRv2HMROM(&listROM, params->nbmode, params->tRef, params->phiRef, params->fRef, params->m1, params->m2, params->distance);
    /* Process the waveform through the LISA response */
    LISASimFDResponseTDI(&listROM, &listA, &listE, &listT, params->inclination, params->lambda, params->beta, params->psi);

    ListmodesCAmpPhaseFrequencySeries* listelementA;
    ListmodesCAmpPhaseFrequencySeries* listelementE;
    ListmodesCAmpPhaseFrequencySeries* listelementT;
    int l,m;
    int status = 0;
    /* Loop over the modes */
    for(int i=0; i<params->nbmode; i++){
      l = listmode[i][0];
      m = listmode[i][1];
      listelementA = ListmodesCAmpPhaseFrequencySeries_GetMode(listA, l, m);
      listelementE = ListmodesCAmpPhaseFrequencySeries_GetMode(listE, l, m);
      listelementT = ListmodesCAmpPhaseFrequencySeries_GetMode(listT, l, m);
      /* Write files - suffix _A,E,T_lm.dat imposed to each file name, attached to the string given by outname */
      /* If outname is still the default empty string, we do not output and skip this stage */
      if(!(strcmp(params->outname, "") == 0)){
	char *filenameA = malloc(strlen(params->outname)+64);
	char *filenameE = malloc(strlen(params->outname)+64);
	char *filenameT = malloc(strlen(params->outname)+64);
	sprintf(filenameA, "%s%s%d%d%s", params->outname, "_A_", l, m, ".dat");
	sprintf(filenameE, "%s%s%d%d%s", params->outname, "_E_", l, m, ".dat");
	sprintf(filenameT, "%s%s%d%d%s", params->outname, "_T_", l, m, ".dat");
	f = fopen(filenameA, "w");
	status |= Write_CAmpPhaseFrequencySeries(f, listelementA->freqseries);
	fclose(f);
	f = fopen(filenameE, "w");
	status |= Write_CAmpPhaseFrequencySeries(f, listelementE->freqseries);
	fclose(f);
	f = fopen(filenameT, "w");
	status |= Write_CAmpPhaseFrequencySeries(f, listelementT->freqseries);
	fclose(f);
	if (status) goto fail;
      }

      /* Example SNR calculation, taking the A observable */
      double *f1=listelementA->freqseries->freq->data;
      int n1=listelementA->freqseries->freq->size;
      double *s1Ar=listelementA->freqseries->amp_real->data;
      double *s1Ai=listelementA->freqseries->amp_imag->data;
      double *s1p=listelementA->freqseries->phase->data;
      double *f2=listelementA->freqseries->freq->data;
      int n2=listelementA->freqseries->freq->size;
      double *s2Ar=listelementA->freqseries->amp_real->data;
      double *s2Ai=listelementA->freqseries->amp_imag->data;
      double *s2p=listelementA->freqseries->phase->data;

      printf("n1: %d\n", n1);
      printf("s1Ar[299]: %g\n", s1Ar[299]);

      double start=((double)clock())/CLOCKS_PER_SEC;
      double rho2= wip_phase (f1, n1, f2, n2, s1Ar, s1Ai, s1p, s2Ar, s2Ai, s2p, Swhitenoise, 1.0);
      double end=((double)clock())/CLOCKS_PER_SEC;
      printf( "SNR2 = %g,   SNR time = %g\n", rho2, end-start);
      printf( "SNR = %g,   SNR time = %g\n",sqrt(rho2), end-start);

    }

    /* clean up */
    free(params);
    ListmodesCAmpPhaseFrequencySeries_Destroy(listROM);
    ListmodesCAmpPhaseFrequencySeries_Destroy(listA);
    ListmodesCAmpPhaseFrequencySeries_Destroy(listE);
    ListmodesCAmpPhaseFrequencySeries_Destroy(listT);
    return 0;

    fail:
    free(params);
    ListmodesCAmpPhaseFrequencySeries_Destroy(listROM);
    ListmodesCAmpPhaseFrequencySeries_Destroy(listA);
    ListmodesCAmpPhaseFrequencySeries_Destroy(listE);
    ListmodesCAmpPhaseFrequencySeries_Destroy(listT);
    return 1;
}

