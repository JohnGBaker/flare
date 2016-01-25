#include "GenerateTDIFD.h"

/************ Parsing arguments function ************/

/* Function to convert string input TDI string to TDItag */
static GenTDIFDtag ParseGenTDIFDtag(char* string) {
  GenTDIFDtag tag;
  if(strcmp(string, "TDIhlm")==0) tag = TDIhlm;
  else if(strcmp(string, "TDIFD")==0) tag = TDIFD;
  else {
    printf("Error in ParseGenTDIFDtag: string not recognized.\n");
    exit(1);
  }
  return tag;
}

/* Parse command line to initialize GenTDITDparams object */
/* Masses are input in solar masses and distances in Mpc - converted in SI for the internals */
static void parse_args_GenerateTDIFD(ssize_t argc, char **argv, GenTDITDparams* params)
{
    char help[] = "\
GenerateWaveform by Sylvain Marsat, John Baker, and Philip Graff\n\
Copyright July 2015\n\
\n\
This program either:\n\
(i) generates a EOBNRv2HMROM waveform, processes it through the Fourier domain response, and outputs it either in the form of Amp/Phase mode contributions or in the form of Re/Im frequency series\n\
(ii) takes as input time series for 3 TDI observables, FFT it and outputs the Re/Im frequency series. Separate output files for TDI channels\n\
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
 --fLow                Minimal frequency (Hz, default=0) - when too low, use first frequency covered by the ROM\n\
 --tagtdi              Tag choosing the set of TDI variables to use (default TDIAETXYZ)\n\
 --taggenwave          Tag choosing the wf format: TDIhlm (default: downsampled mode contributions to TDI in Amp/Phase form), TDIFD (hlm interpolated and summed accross modes)\n\
 --restorescaledfactor Option to restore the factors scaled out of TDI observables (default: false)\n	\
 --fromtditdfile       Option for loading time series for TDI observables and FFTing (default: false)\n\
 --nlinesinfile        Number of lines of inputs file when loading TDI time series from file\n\
 --indir               Input directory when loading TDI time series from file\n\
 --infile              Input file name when loading TDI time series from file\n\
 --outdir              Output directory\n\
 --outfileprefix       Output file name prefix, will output one file for each TDI channel with a built-in postfix\n\
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
    params->lambda = 0;
    params->beta = 0;
    params->polarization = 0;

    /* Set default values for the generation params */
    params->nbmode = 5;
    params->fLow = 0.;
    params->tagtdi = TDIXYZ;
    params->taggenwave = TDIhlm;
    params->restorescaledfactor = 0;
    params->fromtditdfile = 0;
    params->nlinesinfile = 0;    /* No default; has to be provided */
    strcpy(params->indir, "");   /* No default; has to be provided */
    strcpy(params->infile, "");  /* No default; has to be provided */
    strcpy(params->outdir, ".");
    strcpy(params->outfileprefix, "generated_tdiFD");

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
        } else if (strcmp(argv[i], "--lambda") == 0) {
            params->lambda = atof(argv[++i]);
        } else if (strcmp(argv[i], "--beta") == 0) {
            params->beta = atof(argv[++i]);
        } else if (strcmp(argv[i], "--polarization") == 0) {
            params->polarization = atof(argv[++i]);
        } else if (strcmp(argv[i], "--nbmode") == 0) {
            params->nbmode = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--fLow") == 0) {
            params->fLow = atof(argv[++i]);
        } else if (strcmp(argv[i], "--tagtdi") == 0) {
	  params->tagtdi = ParseTDItag(argv[++i]);
        } else if (strcmp(argv[i], "--taggenwave") == 0) {
	  params->taggenwave = ParseGenTDIFDtag(argv[++i]);
        } else if (strcmp(argv[i], "--restorescaledfactor") == 0) {
            params->restorescaledfactor = 1;
        } else if (strcmp(argv[i], "--fromtditdfile") == 0) {
            params->fromtditdfile = 1;
        } else if (strcmp(argv[i], "--nlinesinfile") == 0) {
            params->nlinesinfile = atoi(argv[++i]);
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

/* Built-in postfix for output files for TDI channels */
static char* TDIFilePostfix(TDItag tditag, int channel)
{
  if(tditag==TDIXYZ) {
    switch(channel) {
    case 1: return "_TDIX.txt"; break;
    case 2: return "_TDIY.txt"; break;
    case 3: return "_TDIZ.txt"; break;
    }
  }
  else if(tditag==TDIAETXYZ) {
    switch(channel) {
    case 1: return "_TDIAXYZ.txt"; break;
    case 2: return "_TDIEXYZ.txt"; break;
    case 3: return "_TDITXYZ.txt"; break;
    }
  }
  else {
    printf("Error: in TDIFilePostfix TDI tag not recognized.\n");
    exit(1);
  }
}

/************ Functions to write waveforms to file ************/

/* Read waveform time series in Re/Im form for hpTD and hcTD a single file */
/* NOTE: assumes the same number of points is used to represent each mode */
static void Read_Text_TDITD3Chan( RealTimeSeries** TDI1, RealTimeSeries** TDI2, RealTimeSeries** TDI3, const char dir[], const char file[], const int nblines)
{
  /* Initalize and read input */
  gsl_matrix* inmatrix =  gsl_matrix_alloc(nblines, 4);
  Read_Text_Matrix(dir, file, inmatrix);

  /* Initialize structures */
  RealTimeSeries_Init(TDI1, nblines);
  RealTimeSeries_Init(TDI2, nblines);
  RealTimeSeries_Init(TDI3, nblines);

  /* Set values */
  gsl_vector_view timesview = gsl_matrix_column(inmatrix, 0);
  gsl_vector_view TDI1view = gsl_matrix_column(inmatrix, 1);
  gsl_vector_view TDI2view = gsl_matrix_column(inmatrix, 2);
  gsl_vector_view TDI3view = gsl_matrix_column(inmatrix, 3);
  gsl_vector_memcpy((*TDI1)->times, &timesview.vector);
  gsl_vector_memcpy((*TDI2)->times, &timesview.vector);
  gsl_vector_memcpy((*TDI3)->times, &timesview.vector);
  gsl_vector_memcpy((*TDI1)->h, &TDI1view.vector);
  gsl_vector_memcpy((*TDI2)->h, &TDI2view.vector);
  gsl_vector_memcpy((*TDI3)->h, &TDI3view.vector);

  /* Clean up */
  gsl_matrix_free(inmatrix);
}
/* Output waveform in frequency series form,Re/Im for hplus and hcross */
static void Write_Text_FrequencySeries(const char dir[], const char file[], ReImFrequencySeries* freqseries)
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
  Write_Text_Matrix(dir, file, outmatrix);
}
/* Output TDI mode contributions in downsampled form, FD AmpReal/AmpIm/Phase, all hlm modes in a single file */
/* NOTE: assumes the same number of points is used to represent each mode */
static void Write_Text_TDIhlm(const char dir[], const char file[], ListmodesCAmpPhaseFrequencySeries* listhlm, int nbmodes)
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
  Write_Text_Matrix(dir, file, outmatrix);
}

/***************** Main program *****************/

int main(int argc, char *argv[])
{
  /* Initialize structure for parameters */
  GenTDITDparams* params;
  params = (GenTDITDparams*) malloc(sizeof(GenTDITDparams));
  memset(params, 0, sizeof(GenTDITDparams));
  
  /* Parse commandline to read parameters */
  parse_args_GenerateTDIFD(argc, argv, params);

  if(params->fromtditdfile) {
    /* Load TD TDI from file */
    RealTimeSeries* TDI1 = NULL;
    RealTimeSeries* TDI2 = NULL;
    RealTimeSeries* TDI3 = NULL;
    Read_Text_TDITD3Chan(&TDI1, &TDI2, &TDI3, params->indir, params->infile, params->nlinesinfile);

    /* Compute FFT */
    ReImFrequencySeries* TDI1FFT = NULL;
    ReImFrequencySeries* TDI2FFT = NULL;
    ReImFrequencySeries* TDI3FFT = NULL;
    gsl_vector* times = TDI1->times;
    double twindowbeg = 0.05 * (gsl_vector_get(times, times->size - 1) - gsl_vector_get(times, 0)); /* Here hardcoded relative window lengths */
    double twindowend = 0.01 * (gsl_vector_get(times, times->size - 1) - gsl_vector_get(times, 0)); /* Here hardcoded relative window lengths */
    FFTTimeSeries(&TDI1FFT, TDI1, twindowbeg, twindowend, 2); /* Here hardcoded 0-padding */
    FFTTimeSeries(&TDI2FFT, TDI2, twindowbeg, twindowend, 2); /* Here hardcoded 0-padding */
    FFTTimeSeries(&TDI3FFT, TDI3, twindowbeg, twindowend, 2); /* Here hardcoded 0-padding */

    /* Output */
    char *outfileTDI1 = malloc(256);
    char *outfileTDI2 = malloc(256);
    char *outfileTDI3 = malloc(256);
    sprintf(outfileTDI1, "%s%s", params->outfileprefix, TDIFilePostfix(params->tagtdi, 1));
    sprintf(outfileTDI2, "%s%s", params->outfileprefix, TDIFilePostfix(params->tagtdi, 2));
    sprintf(outfileTDI3, "%s%s", params->outfileprefix, TDIFilePostfix(params->tagtdi, 3));
    Write_Text_FrequencySeries(params->outdir, outfileTDI1, TDI1FFT);
    Write_Text_FrequencySeries(params->outdir, outfileTDI2, TDI2FFT);
    Write_Text_FrequencySeries(params->outdir, outfileTDI3, TDI3FFT);
    free(outfileTDI1);
    free(outfileTDI2);
    free(outfileTDI3);

    exit(0);
  }

  else {
    /* Set geometric coefficients */
    SetCoeffsG(params->lambda, params->beta, params->polarization);

    /* Generate Fourier-domain waveform as a list of hlm modes */
    ListmodesCAmpPhaseFrequencySeries* listROM = NULL;
    SimEOBNRv2HMROM(&listROM, params->nbmode, params->tRef, params->phiRef, params->fRef, params->m1*MSUN_SI, params->m2*MSUN_SI, params->distance*1e6*PC_SI);

    /* Process through the Fourier-domain response for TDI observables */
    ListmodesCAmpPhaseFrequencySeries* listTDI1= NULL;
    ListmodesCAmpPhaseFrequencySeries* listTDI2= NULL;
    ListmodesCAmpPhaseFrequencySeries* listTDI3= NULL;
    LISASimFDResponseTDI3Chan(&listROM, &listTDI1, &listTDI2, &listTDI3, params->tRef, params->lambda, params->beta, params->inclination, params->polarization, params->tagtdi);

    /* If asked for it, rescale the complex amplitudes to include the factors that were scaled out of TDI observables */
    if(params->restorescaledfactor) {
      RestoreInPlaceScaledFactorTDI(listTDI1, params->tagtdi, 1);
      RestoreInPlaceScaledFactorTDI(listTDI2, params->tagtdi, 2);
      RestoreInPlaceScaledFactorTDI(listTDI3, params->tagtdi, 3);
    }

    if(params->taggenwave==TDIhlm) {
      
      /* Output */
      char *outfileTDI1hlm = malloc(256);
      char *outfileTDI2hlm = malloc(256);
      char *outfileTDI3hlm = malloc(256);
      sprintf(outfileTDI1hlm, "%s%s", params->outfileprefix, TDIFilePostfix(params->tagtdi, 1));
      sprintf(outfileTDI2hlm, "%s%s", params->outfileprefix, TDIFilePostfix(params->tagtdi, 2));
      sprintf(outfileTDI3hlm, "%s%s", params->outfileprefix, TDIFilePostfix(params->tagtdi, 3));
      Write_Text_TDIhlm(params->outdir, outfileTDI1hlm, listTDI1, params->nbmode);
      Write_Text_TDIhlm(params->outdir, outfileTDI2hlm, listTDI2, params->nbmode);
      Write_Text_TDIhlm(params->outdir, outfileTDI3hlm, listTDI3, params->nbmode);
      free(outfileTDI1hlm);
      free(outfileTDI2hlm);
      free(outfileTDI3hlm);

      exit(0);
    }

    else {

      /* Determine deltaf so that N deltat = 1/deltaf > 2*tc where tc is the time to coalescence estimated from Psi22 */
      /* Assumes the TD waveform ends around t=0 */
      /* Note: the phase is untouched by the response processing, so the phase in the 22 TDI contribution is still Psi22 and the same for all channels */
      double tc = EstimateInitialTime(listTDI1, params->fLow);
      double deltaf = 0.5 * 1./(2*(-tc)); /* Extra factor of 1/2 corresponding to 0-padding in TD by factor of 2 */
      
      /* Generate FD frequency series by summing mode contributions */
      ReImFrequencySeries* TDI1FD = NULL;
      ReImFrequencySeries* TDI2FD = NULL;
      ReImFrequencySeries* TDI3FD = NULL;
      GenerateFDReImFrequencySeries(&TDI1FD, listTDI1, params->fLow, deltaf, 0, params->nbmode); /* Here determines the number of points from deltaf and max frequency in input list of modes */
      GenerateFDReImFrequencySeries(&TDI2FD, listTDI2, params->fLow, deltaf, 0, params->nbmode); /* Here determines the number of points from deltaf and max frequency in input list of modes */
      GenerateFDReImFrequencySeries(&TDI3FD, listTDI3, params->fLow, deltaf, 0, params->nbmode); /* Here determines the number of points from deltaf and max frequency in input list of modes */

      /* Output */
      char *outfileTDI1FD = malloc(256);
      char *outfileTDI2FD = malloc(256);
      char *outfileTDI3FD = malloc(256);
      sprintf(outfileTDI1FD, "%s%s", params->outfileprefix, TDIFilePostfix(params->tagtdi, 1));
      sprintf(outfileTDI2FD, "%s%s", params->outfileprefix, TDIFilePostfix(params->tagtdi, 2));
      sprintf(outfileTDI3FD, "%s%s", params->outfileprefix, TDIFilePostfix(params->tagtdi, 3));
      Write_Text_FrequencySeries(params->outdir, outfileTDI1FD, TDI1FD);
      Write_Text_FrequencySeries(params->outdir, outfileTDI2FD, TDI2FD);
      Write_Text_FrequencySeries(params->outdir, outfileTDI3FD, TDI3FD);
      free(outfileTDI1FD);
      free(outfileTDI2FD);
      free(outfileTDI3FD);

      exit(0);
    }
  }
}
