/**
 * \author Sylvain Marsat, AEI
 *
 * \brief C test code used to check that IMRPhenomD Amplitude/Phase gives identical results to the LDC implementation
 *
 */

// build commmand:
// gcc -c test_waveform.c -O0 -g -std=c99 -I../tools -I../IMRPhenomD -I/usr/local/include
// gcc -o test_waveform test_waveform.o IMRPhenomD_internals.o IMRPhenomD.o ../tools/struct.o -O0 -g -std=c99 -I/usr/local/include -I../tools -I../integration -I../IMRPhenomD -lm -lgsl -lgslcblas -fopenmp


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
#include "IMRPhenomD_internals.h"
#include "IMRPhenomD.h"


// /* Parameters for the generation of a ROM waveform (in the form of a list of modes) */
// /* All parameters are to be given in SI units! */
// typedef struct tagROMParams {
//   int nbmode;                /* Number of modes to generate (starting with the 22) - defaults to 1 (22 mode only) */
//   double deltatRef;               /* shift in time with respect to the approximately 22-peak-aligned waveform */
//   double phiRef;             /* phase at fRef */
//   double fRef;               /* reference frequency */
//   double m1;                 /* mass of companion 1 */
//   double m2;                 /* mass of companion 2 */
//   double distance;           /* distance of source */
//   char outname[256];         /* file to which output should be written */
// } ROMParams;
//
// /* Parse command line and return a newly allocated ROMParams object
//  * Masses are input in solar masses and distances in Mpc - converted in SI for the internals */
// static ROMParams* parse_args_ROM(ssize_t argc, char **argv) {
//     ssize_t i;
//     ROMParams* params;
//     params = (ROMParams*) malloc(sizeof(ROMParams));
//     memset(params, 0, sizeof(ROMParams));
//
//     /* Set default values to the arguments */
//     params->nbmode = 5;
//     params->deltatRef = 0.;
//     params->phiRef = 0.;
//     params->fRef = 0.;
//     params->m1 = 10.;
//     params->m2 = 10.;
//     params->distance = 200. * 1e6;
//
//     /* consume command line */
//     for (i = 1; i < argc; ++i) {
//         if (strcmp(argv[i], "--nbmode") == 0) {
//             params->nbmode = atof(argv[++i]);
//         } else if (strcmp(argv[i], "--deltatRef") == 0) {
//             params->deltatRef = atof(argv[++i]);
//         } else if (strcmp(argv[i], "--phiRef") == 0) {
//             params->phiRef = atof(argv[++i]);
//         } else if (strcmp(argv[i], "--fRef") == 0) {
//             params->fRef = atof(argv[++i]);
//         } else if (strcmp(argv[i], "--m1") == 0) {
//             params->m1 = atof(argv[++i]) * MSUN_SI;
//         } else if (strcmp(argv[i], "--m2") == 0) {
//             params->m2 = atof(argv[++i]) * MSUN_SI;
//         } else if (strcmp(argv[i], "--distance") == 0) {
//             params->distance = atof(argv[++i]) * PC_SI;
//         } else if (strcmp(argv[i], "--outname") == 0) {
//             strncpy(params->outname, argv[++i], 256);
//         } else {
//             printf("Error: invalid option: %s\n", argv[i]);
//             goto fail;
//         }
//     }
//
//     return params;
//
//     fail:
//     free(params);
//     exit(1);
// }
//
// /* Function to output to a file an AmpPhaseFrequencySeries; each mode will be output to a separate file using this function */
// static int Write_CAmpPhaseFrequencySeries(FILE* f, CAmpPhaseFrequencySeries* freqseries) {
//     gsl_vector* freq = freqseries->freq;
//     gsl_vector* amp_real = freqseries->amp_real;
//     gsl_vector* amp_imag = freqseries->amp_imag;
//     gsl_vector* phase = freqseries->phase;
//
//     int len = (int) freq->size;
//     /*Here, we could add a check on the length of the gsl_vectors*/
//
//     fprintf(f, "# f amp_re amp_im phase\n");
//     for (int i=0; i<len; i++) {
//       fprintf(f, "%.16e %.16e %.16e %.16e\n", gsl_vector_get(freq, i), gsl_vector_get(amp_real, i), gsl_vector_get(amp_imag, i), gsl_vector_get(phase, i));
//     }
//
//     return 0;
// }

/*
 * main
 */
int main (int argc , char **argv) {

    /* I/O directory */
    char dir[] = "/Users/marsat/data/flare/test/test_IMRPhenomD";

    /* Hardcoded params */
    double m1_SI = 60. * MSUN_SI;
    double m2_SI = 20. * MSUN_SI;
    double chi1 = 0.5;
    double chi2 = -0.5;
    double dist_SI = 100e6 * PC_SI;
    double fRef_in = 0.;
    double phiRef = PI/4;
    double tRef = 0.;

    /* Read the frequency vector on which to evaluate */
    int len = 1000;
    char freq_in[] = "test_LDC_f.txt";
    gsl_vector* freq_vector = gsl_vector_alloc(len);
    Read_Text_Vector(dir, freq_in, freq_vector);

    /* Generate the waveform */
    ListmodesCAmpPhaseFrequencySeries* listhlm = NULL;
    IMRPhenomDGenerateh22FDAmpPhase(
        &listhlm, /**< [out] FD waveform */
        freq_vector,                   /**< Input: frequencies (Hz) on which to evaluate h22 FD - will be copied in the output AmpPhaseFDWaveform. Frequencies exceeding max freq covered by PhenomD will be given 0 amplitude and phase. */
        tRef,                  /**< Time shift, 0 corresponds to peak approximatively set to time=0 (s) */
        phiRef,                /**< Orbital phase at fRef (rad) */
        fRef_in,               /**< reference frequency (Hz) */
        m1_SI,                 /**< Mass of companion 1 (kg) */
        m2_SI,                 /**< Mass of companion 2 (kg) */
        chi1,                  /**< Aligned-spin parameter of companion 1 */
        chi2,                  /**< Aligned-spin parameter of companion 2 */
        dist_SI               /**< Distance of source (m) */
    );

    /* Output the waveform */
    char freq_out[] = "test_flare_f.txt";
    char amp_out[] = "test_flare_amp.txt";
    char phase_out[] = "test_flare_phase.txt";
    CAmpPhaseFrequencySeries* h22 = listhlm->freqseries;
    Write_Text_Vector(dir, freq_out, h22->freq);
    Write_Text_Vector(dir, amp_out, h22->amp_real);
    Write_Text_Vector(dir, phase_out, h22->phase);

    /* Clean up */
    ListmodesCAmpPhaseFrequencySeries_Destroy(listhlm);
    return 0;
}
