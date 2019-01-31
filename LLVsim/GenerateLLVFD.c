#include "GenerateLLVFD.h"

/************ Parsing arguments function ************/

/* Function to convert string input LLV string to LLVtag */
static GenLLVFDtag ParseGenLLVFDtag(char* string) {
  GenLLVFDtag tag;
  if(strcmp(string, "LLVhlm")==0) tag = LLVhlm;
  else if(strcmp(string, "LLVFD")==0) tag = LLVFD;
  else {
    printf("Error in ParseGenLLVFDtag: string not recognized.\n");
    exit(1);
  }
  return tag;
}

/* Parse command line to initialize GenLLVFDparams object */
/* Masses are input in solar masses and distances in Mpc - converted in SI for the internals */
static void parse_args_GenerateLLVFD(ssize_t argc, char **argv, GenLLVFDparams* params)
{
    char help[] = "\
GenerateWaveform by Sylvain Marsat, John Baker, and Philip Graff\n\
Copyright July 2015\n\
\n\
This program either:\n\
(i) generates a EOBNRv2HMROM waveform, processes it through the Fourier domain response, and outputs it either in the form of Amp/Phase mode contributions or in the form of Re/Im frequency series\n\
(ii) takes as input time series for 3 LLV observables, FFT it and outputs the Re/Im frequency series. Separate output files for LLV channels\n\
Arguments are as follows:\n\
\n\
--------------------------------------------------\n\
----- Physical Parameters ------------------------\n\
--------------------------------------------------\n\
 --tRef                Time at reference frequency (sec, default=0)\n\
 --phiRef              Orbital phase at reference frequency (radians, default=0)\n\
 --fRef                Reference frequency (Hz, default=0, interpreted as Mf=0.14)\n\
 --m1                  Component mass 1 in Solar masses (larger, default=2e6)\n\
 --m2                  Component mass 2 in Solar masses (smaller, default=1e6)\n\
 --distance            Distance to source in Mpc (default=1e3)\n\
 --inclination         Inclination of source orbital plane to observer line of sight\n\
                       (radians, default=PI/3)\n\
 --lambda              First angle for the position in the sky (radians, default=0)\n\
 --beta                Second angle for the position in the sky (radians, default=0)\n\
 --polarization        Polarization of source (radians, default=0)\n\
\n\
--------------------------------------------------\n\
----- Generation Parameters ----------------------\n\
--------------------------------------------------\n\
 --nbmode              Number of modes of radiation to generate (1-5, default=5)\n\
 --minf                Minimal frequency (Hz, default=0) - when too low, use first frequency covered by the ROM\n\
 --setphiRefatfRef     Flag for adjusting the FD phase at phiRef at the given fRef, which depends also on tRef - if false, treat phiRef simply as an orbital phase shift (minus an observer phase shift) (default=1)\n\
 --tagnetwork          Tag choosing the set of detectors to use (default LHV)\n\
 --taggenwave          Tag choosing the wf format: LLVhlm (default: downsampled mode contributions to LLV in Amp/Phase form), LLVFD (hlm interpolated and summed accross modes)\n\
 --restorescaledfactor Option to restore the factors scaled out of LLV observables (default: false)\n	\
 --fromLLVtdfile       Option for loading time series for LLV observables and FFTing (default: false)\n\
 --nsamplesinfile      Number of lines of input file when loading LLV time series from file\n\
 --binaryin            Tag for loading the data in gsl binary form instead of text (default false)\n\
 --binaryout           Tag for outputting the data in gsl binary form instead of text (default false)\n\
 --indir               Input directory when loading LLV time series from file\n\
 --infile              Input file name when loading LLV time series from file\n\
 --outdir              Output directory\n\
 --outfileprefix       Output file name prefix, will output one file for each LLV channel with a built-in postfix\n\
\n";

    ssize_t i;

    /* Set default values for the physical params */
    params->tRef = 0.;
    params->phiRef = 0.;
    params->fRef = 0.;
    params->m1 = 2*1e6;
    params->m2 = 1*1e6;
    params->distance = 1e3;
    params->inclination = PI/3;
    params->ra = 0;
    params->dec = 0;
    params->polarization = 0;

    /* Set default values for the generation params */
    params->nbmode = 5;
    params->minf = 0.;
    params->setphiRefatfRef = 1;
    params->tagnetwork = LHV;
    params->taggenwave = LLVhlm;
    params->fromLLVtdfile = 0;
    params->nsamplesinfile = 0;    /* No default; has to be provided */
    strcpy(params->indir, "");   /* No default; has to be provided */
    strcpy(params->infile, "");  /* No default; has to be provided */
    strcpy(params->outdir, ".");
    strcpy(params->outfileprefix, "generated_LLVFD");

    /* Consume command line */
    for (i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--help") == 0) {
            fprintf(stdout,"%s", help);
            exit(0);
        } else if (strcmp(argv[i], "--tRef") == 0) {
            params->tRef = atof(argv[++i]);
        } else if (strcmp(argv[i], "--phiRef") == 0) {
            params->phiRef = atof(argv[++i]);
        } else if (strcmp(argv[i], "--fRef") == 0) {
            params->fRef = atof(argv[++i]);
        } else if (strcmp(argv[i], "--m1") == 0) {
            params->m1 = atof(argv[++i]);
        } else if (strcmp(argv[i], "--m2") == 0) {
            params->m2 = atof(argv[++i]);
        } else if (strcmp(argv[i], "--distance") == 0) {
            params->distance = atof(argv[++i]);
        } else if (strcmp(argv[i], "--inclination") == 0) {
            params->inclination = atof(argv[++i]);
        } else if (strcmp(argv[i], "--ra") == 0) {
            params->ra = atof(argv[++i]);
        } else if (strcmp(argv[i], "--dec") == 0) {
            params->dec = atof(argv[++i]);
        } else if (strcmp(argv[i], "--polarization") == 0) {
            params->polarization = atof(argv[++i]);
        } else if (strcmp(argv[i], "--nbmode") == 0) {
            params->nbmode = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--minf") == 0) {
            params->minf = atof(argv[++i]);
        } else if (strcmp(argv[i], "--setphiRefatfRef") == 0) {
          params->setphiRefatfRef = atof(argv[++i]);
        } else if (strcmp(argv[i], "--tagnetwork") == 0) {
	          params->tagnetwork = ParseNetworktag(argv[++i]);
        } else if (strcmp(argv[i], "--taggenwave") == 0) {
	          params->taggenwave = ParseGenLLVFDtag(argv[++i]);
        } else if (strcmp(argv[i], "--fromLLVtdfile") == 0) {
            params->fromLLVtdfile = 1;
        } else if (strcmp(argv[i], "--nsamplesinfile") == 0) {
            params->nsamplesinfile = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--binaryin") == 0) {
          params->binaryin = 1;
        } else if (strcmp(argv[i], "--binaryout") == 0) {
          params->binaryout = 1;
        } else if (strcmp(argv[i], "--indir") == 0) {
            strcpy(params->indir, argv[++i]);
        } else if (strcmp(argv[i], "--infile") == 0) {
            strcpy(params->infile, argv[++i]);
        } else if (strcmp(argv[i], "--outdir") == 0) {
            strcpy(params->outdir, argv[++i]);
        } else if (strcmp(argv[i], "--outfileprefix") == 0) {
            strcpy(params->outfileprefix, argv[++i]);
        } else {
	  printf("Error: invalid option: %s\n", argv[i]);
	  goto fail;
        }
    }

    return;

 fail:
    exit(1);
}

/* Built-in postfix for output files for LLV channels */
static char* LLVFilePostfix(int channel, int binary)
{
  if(!binary) {
      switch(channel) {
        case 1: return "_L.txt"; break;
        case 2: return "_H.txt"; break;
        case 3: return "_V.txt"; break;
      }
  }
  else {
      switch(channel) {
        case 1: return "_TDIX.dat"; break;
        case 2: return "_TDIY.dat"; break;
        case 3: return "_TDIZ.dat"; break;
      }
  }
}

/************ Functions to write waveforms to file ************/

/* Read waveform time series in Re/Im form for the 3 LLV detectors in a single file */
/* NOTE: assumes the same number of points is used to represent each mode */
static void Read_LLVTD3Det( RealTimeSeries** LLV1, RealTimeSeries** LLV2, RealTimeSeries** LLV3, const char dir[], const char file[], const int nblines, const int binary)
{
  /* Initalize and read input */
  gsl_matrix* inmatrix =  gsl_matrix_alloc(nblines, 4);
  if(!binary) Read_Text_Matrix(dir, file, inmatrix);
  else Read_Matrix(dir, file, inmatrix);

  /* Initialize structures */
  RealTimeSeries_Init(LLV1, nblines);
  RealTimeSeries_Init(LLV2, nblines);
  RealTimeSeries_Init(LLV3, nblines);

  /* Set values */
  gsl_vector_view timesview = gsl_matrix_column(inmatrix, 0);
  gsl_vector_view LLV1view = gsl_matrix_column(inmatrix, 1);
  gsl_vector_view LLV2view = gsl_matrix_column(inmatrix, 2);
  gsl_vector_view LLV3view = gsl_matrix_column(inmatrix, 3);
  gsl_vector_memcpy((*LLV1)->times, &timesview.vector);
  gsl_vector_memcpy((*LLV2)->times, &timesview.vector);
  gsl_vector_memcpy((*LLV3)->times, &timesview.vector);
  gsl_vector_memcpy((*LLV1)->h, &LLV1view.vector);
  gsl_vector_memcpy((*LLV2)->h, &LLV2view.vector);
  gsl_vector_memcpy((*LLV3)->h, &LLV3view.vector);

  /* Clean up */
  gsl_matrix_free(inmatrix);
}
/* Output waveform in frequency series form,Re/Im for hplus and hcross */
static void Write_FrequencySeries(const char dir[], const char file[], ReImFrequencySeries* freqseries, const int binary)
{
  /* Initialize output */
  /* Note: assumes hplus, hcross have same length as expected */
  int nbfreq = freqseries->freq->size;
  gsl_matrix* outmatrix = gsl_matrix_alloc(nbfreq, 3);

  /* Set output matrix */
  gsl_matrix_set_col(outmatrix, 0, freqseries->freq);
  gsl_matrix_set_col(outmatrix, 1, freqseries->h_real);
  gsl_matrix_set_col(outmatrix, 2, freqseries->h_imag);

  /* Output */
  if (!binary) Write_Text_Matrix(dir, file, outmatrix);
  else Write_Matrix(dir, file, outmatrix);
}
/* Output LLV mode contributions in downsampled form, FD AmpReal/AmpIm/Phase, all hlm modes in a single file */
/* NOTE: assumes the same number of points is used to represent each mode */
static void Write_LLVhlm(const char dir[], const char file[], ListmodesCAmpPhaseFrequencySeries* listhlm, int nbmodes, const int binary)
{
  /* Initialize output */
  /* get length from 22 mode - NOTE: assumes the same for all modes */
  int nbfreq = ListmodesCAmpPhaseFrequencySeries_GetMode(listhlm, 2, 2)->freqseries->freq->size;
  gsl_matrix* outmatrix = gsl_matrix_alloc(nbfreq, 4*nbmodes);

  /* Get data in the list of modes */
  CAmpPhaseFrequencySeries* mode;
  for(int i=0; i<nbmodes; i++) {
    mode = ListmodesCAmpPhaseFrequencySeries_GetMode(listhlm, listmode[i][0], listmode[i][1])->freqseries;
    gsl_matrix_set_col(outmatrix, 0+3*i, mode->freq);
    gsl_matrix_set_col(outmatrix, 1+3*i, mode->amp_real);
    gsl_matrix_set_col(outmatrix, 2+3*i, mode->amp_imag); /* amp_imag is important here */
    gsl_matrix_set_col(outmatrix, 3+3*i, mode->phase);
  }

  /* Output */
  if (!binary) Write_Text_Matrix(dir, file, outmatrix);
  else Write_Matrix(dir, file, outmatrix);
}

/***************** Main program *****************/

int main(int argc, char *argv[])
{
  /* Initialize structure for parameters */
  GenLLVFDparams* params;
  params = (GenLLVFDparams*) malloc(sizeof(GenLLVFDparams));
  memset(params, 0, sizeof(GenLLVFDparams));

  /* Parse commandline to read parameters */
  parse_args_GenerateLLVFD(argc, argv, params);

  if(params->fromLLVtdfile) {
    /* Load TD LLV from file */
    RealTimeSeries* LLV1 = NULL;
    RealTimeSeries* LLV2 = NULL;
    RealTimeSeries* LLV3 = NULL;
    Read_LLVTD3Det(&LLV1, &LLV2, &LLV3, params->indir, params->infile, params->nsamplesinfile, params->binaryin);

    /* Compute FFT */
    ReImFrequencySeries* LLV1FFT = NULL;
    ReImFrequencySeries* LLV2FFT = NULL;
    ReImFrequencySeries* LLV3FFT = NULL;
    gsl_vector* times = LLV1->times;
    double twindowbeg = 0.05 * (gsl_vector_get(times, times->size - 1) - gsl_vector_get(times, 0)); /* Here hardcoded relative window lengths */
    double twindowend = 0.01 * (gsl_vector_get(times, times->size - 1) - gsl_vector_get(times, 0)); /* Here hardcoded relative window lengths */
    FFTRealTimeSeries(&LLV1FFT, LLV1, twindowbeg, twindowend, 2); /* Here hardcoded 0-padding */
    FFTRealTimeSeries(&LLV2FFT, LLV2, twindowbeg, twindowend, 2); /* Here hardcoded 0-padding */
    FFTRealTimeSeries(&LLV3FFT, LLV3, twindowbeg, twindowend, 2); /* Here hardcoded 0-padding */

    /* Output */
    char *outfileLLV1 = malloc(256);
    char *outfileLLV2 = malloc(256);
    char *outfileLLV3 = malloc(256);
    sprintf(outfileLLV1, "%s%s", params->outfileprefix, LLVFilePostfix(1, params->binaryout));
    sprintf(outfileLLV2, "%s%s", params->outfileprefix, LLVFilePostfix(2, params->binaryout));
    sprintf(outfileLLV3, "%s%s", params->outfileprefix, LLVFilePostfix(3, params->binaryout));
    Write_FrequencySeries(params->outdir, outfileLLV1, LLV1FFT, params->binaryout);
    Write_FrequencySeries(params->outdir, outfileLLV2, LLV2FFT, params->binaryout);
    Write_FrequencySeries(params->outdir, outfileLLV3, LLV3FFT, params->binaryout);
    free(outfileLLV1);
    free(outfileLLV2);
    free(outfileLLV3);

    exit(0);
  }

  else {
    /* Generate Fourier-domain waveform as a list of hlm modes */
    ListmodesCAmpPhaseFrequencySeries* listROM = NULL;
    SimEOBNRv2HMROM(&listROM, params->nbmode, params->tRef, params->phiRef, params->fRef, params->m1*MSUN_SI, params->m2*MSUN_SI, params->distance*1e6*PC_SI, params->setphiRefatfRef);

    /* Process through the Fourier-domain response for LLV observables */
    ListmodesCAmpPhaseFrequencySeries* listLLV1= NULL;
    ListmodesCAmpPhaseFrequencySeries* listLLV2= NULL;
    ListmodesCAmpPhaseFrequencySeries* listLLV3= NULL;
    LLVSimFDResponse3Det(&listLLV1, &listLLV2, &listLLV3, &listROM, params->tRef, params->ra, params->dec, params->inclination, params->polarization, params->tagnetwork);

    if(params->taggenwave==LLVhlm) {

      /* Output */
      char *outfileLLV1hlm = malloc(256);
      char *outfileLLV2hlm = malloc(256);
      char *outfileLLV3hlm = malloc(256);
      sprintf(outfileLLV1hlm, "%s%s", params->outfileprefix, LLVFilePostfix(1, params->binaryout));
      sprintf(outfileLLV2hlm, "%s%s", params->outfileprefix, LLVFilePostfix(2, params->binaryout));
      sprintf(outfileLLV3hlm, "%s%s", params->outfileprefix, LLVFilePostfix(3, params->binaryout));
      Write_LLVhlm(params->outdir, outfileLLV1hlm, listLLV1, params->nbmode, params->binaryout);
      Write_LLVhlm(params->outdir, outfileLLV2hlm, listLLV2, params->nbmode, params->binaryout);
      Write_LLVhlm(params->outdir, outfileLLV3hlm, listLLV3, params->nbmode, params->binaryout);
      free(outfileLLV1hlm);
      free(outfileLLV2hlm);
      free(outfileLLV3hlm);

      exit(0);
    }

    else {

      double deltaf = 0.;
      if(params->deltaf==0.) {
        /* Determine deltaf so that N deltat = 1/deltaf > 2*tc where tc is the time to coalescence estimated from Psi22 */
        /* NOTE: assumes the TD waveform ends around t=0 */
        /* NOTE: estimate based on the 22 mode - fstartobs will be scaled from mode to mode to ensure the same deltatobs for all (except for the 21 mode, which will turn on after the others) */
        /* NOTE: the phase differs for the three detectors, and the signal can be shifted in time at worst by ~40 ms - should be ok as we take some margin */
        double tc = EstimateInitialTime(listROM, params->minf);
        double deltaf = 0.5 * 1./(2*(-tc)); /* Extra factor of 1/2 corresponding to 0-padding in TD by factor of 2 */
      }
      else deltaf = params->deltaf;

      /* Compute TDI FD */
      /* NOTE: we do note use deltatobs and fstartobs */

      /* Generate FD frequency series by summing mode contributions */
      ReImFrequencySeries* LLV1FD = NULL;
      ReImFrequencySeries* LLV2FD = NULL;
      ReImFrequencySeries* LLV3FD = NULL;
      GenerateFDReImFrequencySeries(&LLV1FD, listLLV1, params->minf, params->maxf, 0., deltaf, 0); /* Here determines the number of points from deltaf and max frequency in input list of modes */
      GenerateFDReImFrequencySeries(&LLV2FD, listLLV2, params->minf, params->maxf, 0., deltaf, 0); /* Here determines the number of points from deltaf and max frequency in input list of modes */
      GenerateFDReImFrequencySeries(&LLV3FD, listLLV3, params->minf, params->maxf, 0., deltaf, 0); /* Here determines the number of points from deltaf and max frequency in input list of modes */

      /* Output */
      char *outfileLLV1FD = malloc(256);
      char *outfileLLV2FD = malloc(256);
      char *outfileLLV3FD = malloc(256);
      sprintf(outfileLLV1FD, "%s%s", params->outfileprefix, LLVFilePostfix(1, params->binaryout));
      sprintf(outfileLLV2FD, "%s%s", params->outfileprefix, LLVFilePostfix(2, params->binaryout));
      sprintf(outfileLLV3FD, "%s%s", params->outfileprefix, LLVFilePostfix(3, params->binaryout));
      Write_FrequencySeries(params->outdir, outfileLLV1FD, LLV1FD, params->binaryout);
      Write_FrequencySeries(params->outdir, outfileLLV2FD, LLV2FD, params->binaryout);
      Write_FrequencySeries(params->outdir, outfileLLV3FD, LLV3FD, params->binaryout);
      free(outfileLLV1FD);
      free(outfileLLV2FD);
      free(outfileLLV3FD);

      exit(0);
    }
  }
}
