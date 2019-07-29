#include "ComputeLISASNR.h"

/************ Parsing arguments function ************/

/* Parse command line to initialize GenTDITDparams object */
/* Masses are input in solar masses and distances in Mpc - converted in SI for the internals */
static void parse_args_ComputeLISASNR(ssize_t argc, char **argv, ComputeLISASNRparams* params)
{
  char help[] = "\
ComputeLISASNR by Sylvain Marsat, John Baker, and Philip Graff\n\
Copyright July 2015\n\
\n\
This program computes and prints the SNR for a LISA TDI waveform (currently only AET(XYZ) available); it will either:\n\
(i) generate a EOBNRv2HMROM waveform, process it through the Fourier domain response, and compute SNR using accelerated Fresnel overlap, with rescaling of the noise (follows LISAinference internals)\n\
(ii) generate a EOBNRv2HMROM waveform, process it through the Fourier domain response, and compute SNR using ordinary linear overlap, with rescaling of the noise (follows LISAinference internals)\n\
(iii) load as input time series for 3 TDI AET observables, FFT them, and compute SNR with linear integration, without rescaling of the noise.\n\
Parameters can be given either directly or in the form of a text file, in which case the output is in the form of a text file with the SNR appended at the end of the line.\n\
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
 --deltatobs           Observation duration (years, default=2)\n\
 --minf                Minimal frequency (Hz, default=0) - when set to 0, use the lowest frequency where the detector noise model is trusted __LISASimFD_Noise_fLow (set somewhat arbitrarily)\n\
 --maxf                Maximal frequency (Hz, default=0) - when set to 0, use the highest frequency where the detector noise model is trusted __LISASimFD_Noise_fHigh (set somewhat arbitrarily)\n\
 --tagextpn            Tag to allow PN extension of the waveform at low frequencies (default=1)\n\
 --tagtdi              Tag choosing the set of TDI variables to use (default TDIAETXYZ)\n\
 --tagint              Tag choosing the integrator: 0 for Fresnel (default), 1 for linear integration\n\
 --nbptsoverlap        Number of points to use for linear integration (default 32768)\n\
 --variant             String representing the variant of LISA to be applied (default LISAProposal)\n\
 --frozenLISA          Freeze the orbital configuration to the time of peak of the injection (default 0)\n\
 --responseapprox      Approximation in the GAB and orb response - choices are full (full response, default), lowfL (keep orbital delay frequency-dependence but simplify constellation response) and lowf (simplify constellation and orbital response) - WARNING : at the moment noises are not consistent, and TDI combinations from the GAB are unchanged\n\
 --delaycorrection     Include the first-order ddot delay correction in phaseRdelay (default 1) - NOTE: treated separately from frozenLISA, while strictly speaking ddot should be zero for a frozen LISA\n\
 --fromtditdfile       Option for loading time series for TDI observables and FFTing (default: false)\n\
 --nlinesinfile        Number of lines of inputs file when loading TDI time series from file\n\
 --indir               Input directory when loading TDI time series from file\n\
 --infile              Input file name when loading TDI time series from file\n\
 --loadparamsfile      Option to load physical parameters from file and to output result to file (default false)\n\
 --nlinesparams        Number of lines in params file\n\
 --paramsdir           Directory for input/output file\n\
 --paramsfile          Input file with the parameters\n\
 --outputfile          Output file\n\
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

  params->frozenLISA = 0;
  params->responseapprox = full;
  params->delaycorrection = 1;
  params->variant = &LISAProposal;

  /* Set default values for the generation params */
  params->nbmode = 5;
  params->deltatobs = 2.;
  params->minf = 0.;
  params->maxf = 0.;
  params->tagextpn = 1;
  params->Mfmatch = 0.;
  params->tagtdi = TDIAETXYZ;
  params->tagint = 0;
  params->nbptsoverlap = 32768;
  params->fromtditdfile = 0;
  params->nlinesinfile = 0;    /* No default; has to be provided */
  strcpy(params->indir, "");   /* No default; has to be provided */
  strcpy(params->infile, "");  /* No default; has to be provided */
  params->loadparamsfile = 0;
  params->nlinesparams = 0;
  strcpy(params->paramsdir, "");  /* No default; has to be provided */
  strcpy(params->paramsfile, "");  /* No default; has to be provided */
  strcpy(params->outputfile, "");  /* No default; has to be provided */

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
    } else if (strcmp(argv[i], "--deltatobs") == 0) {
      params->deltatobs = atof(argv[++i]);
    } else if (strcmp(argv[i], "--minf") == 0) {
      params->minf = atof(argv[++i]);
    } else if (strcmp(argv[i], "--maxf") == 0) {
      params->maxf = atof(argv[++i]);
    } else if (strcmp(argv[i], "--tagextpn") == 0) {
        params->tagextpn = atoi(argv[++i]);
    } else if (strcmp(argv[i], "--Mfmatch") == 0) {
        params->Mfmatch = atof(argv[++i]);
    } else if (strcmp(argv[i], "--tagtdi") == 0) {
      params->tagtdi = ParseTDItag(argv[++i]);
    } else if (strcmp(argv[i], "--tagint") == 0) {
      params->tagint = atoi(argv[++i]);
    } else if (strcmp(argv[i], "--nbptsoverlap") == 0) {
      params->nbptsoverlap = atoi(argv[++i]);
    } else if (strcmp(argv[i], "--frozenLISA") == 0) {
        params->frozenLISA = 1;
    } else if (strcmp(argv[i], "--responseapprox") == 0) {
        params->responseapprox = ParseResponseApproxtag(argv[++i]);
    } else if (strcmp(argv[i], "--delaycorrection") == 0) {
      params->delaycorrection = atoi(argv[++i]);
    } else if (strcmp(argv[i], "--variant") == 0) {
      i++;
      if (strcmp(argv[i], "LISAProposal") == 0) params->variant = &LISAProposal;
      else if (strcmp(argv[i], "LISA2017") == 0) params->variant = &LISA2017;
      else if (strcmp(argv[i], "LISA2010") == 0) params->variant = &LISA2010;
      else if (strcmp(argv[i], "fastOrbitLISA") == 0) params->variant = &fastOrbitLISA;
      else if (strcmp(argv[i], "slowOrbitLISA") == 0) params->variant = &slowOrbitLISA;
      else if (strcmp(argv[i], "tinyOrbitLISA") == 0) params->variant = &tinyOrbitLISA;
      else if (strcmp(argv[i], "bigOrbitLISA") == 0) params->variant = &bigOrbitLISA;
      else {
        printf("Error: --variant option '%s' not recognized\n",argv[i]);
        exit(1);
      }
    } else if (strcmp(argv[i], "--fromtditdfile") == 0) {
      params->fromtditdfile = 1;
    } else if (strcmp(argv[i], "--nlinesinfile") == 0) {
      params->nlinesinfile = atoi(argv[++i]);
    } else if (strcmp(argv[i], "--indir") == 0) {
      strcpy(params->indir, argv[++i]);
    } else if (strcmp(argv[i], "--infile") == 0) {
      strcpy(params->infile, argv[++i]);
    } else if (strcmp(argv[i], "--loadparamsfile") == 0) {
      params->loadparamsfile = 1;
    } else if (strcmp(argv[i], "--nlinesparams") == 0) {
      params->nlinesparams = atoi(argv[++i]);
    } else if (strcmp(argv[i], "--paramsdir") == 0) {
      strcpy(params->paramsdir, argv[++i]);
    }  else if (strcmp(argv[i], "--paramsfile") == 0) {
      strcpy(params->paramsfile, argv[++i]);
    } else if (strcmp(argv[i], "--outputfile") == 0) {
      strcpy(params->outputfile, argv[++i]);
    }  else {
      printf("Error: invalid option: %s\n", argv[i]);
      printf("argc-i=%i\n",argc-i);
      goto fail;
    }
  }

  /* Set frequency interval to default values */
  if(params->minf==0.) params->minf = __LISASimFD_Noise_fLow;
  if(params->maxf==0.) params->maxf = __LISASimFD_Noise_fHigh;

  return;

 fail:
  exit(1);
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

/***************** Function to control that the TDI tag is allowed *****************/
/* Only tdi allowed : A,E,T (built from X,Y,Z), either individually or all three together */
static int AllowedTDItag(TDItag tag) {
  int ret = 0;
  if(tag==TDIAETXYZ) ret = 1;
  else if(tag==TDIAXYZ) ret = 1;
  else if(tag==TDIEXYZ) ret = 1;
  else if(tag==TDITXYZ) ret = 1;
  else ret = 0;
  return ret;
}

/***************** Main program *****************/

int main(int argc, char *argv[])
{
  /* These global parameters are set by command line in other programs but fixed here. */
  //LISAconstellation *variant = &LISAProposal;
  int tagtRefatLISA = 0;
  int tagsimplelikelihood22 = 0;
  int tagsimplelikelihoodHM = 0;
  int zerolikelihood = 0;
  // int frozenLISA = 0;
  // ResponseApproxtag responseapprox = full;

  double SNR = 0;

  /* Initialize structure for parameters */
  ComputeLISASNRparams* params;
  params = (ComputeLISASNRparams*) malloc(sizeof(ComputeLISASNRparams));
  memset(params, 0, sizeof(ComputeLISASNRparams));

  /* Parse commandline to read parameters */
  parse_args_ComputeLISASNR(argc, argv, params);

  /* NOTE: supports only AET(XYZ) (orthogonal) */
  if(!(AllowedTDItag(params->tagtdi))) {
    printf("Error in ComputeLISASNR: TDI tag not recognized.\n");
    exit(1);
  }

  else {

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
      FFTRealTimeSeries(&TDI1FFT, TDI1, twindowbeg, twindowend, 2); /* Here hardcoded 0-padding */
      FFTRealTimeSeries(&TDI2FFT, TDI2, twindowbeg, twindowend, 2); /* Here hardcoded 0-padding */
      FFTRealTimeSeries(&TDI3FFT, TDI3, twindowbeg, twindowend, 2); /* Here hardcoded 0-padding */

      /* Restrict FFT on frequency interval of interest - no limitation put on the upper bound */
      ReImFrequencySeries* TDI1FFTrestr = NULL;
      ReImFrequencySeries* TDI2FFTrestr = NULL;
      ReImFrequencySeries* TDI3FFTrestr = NULL;
      RestrictFDReImFrequencySeries(&TDI1FFTrestr, TDI1FFT, params->minf, __LISASimFD_Noise_fHigh);
      RestrictFDReImFrequencySeries(&TDI2FFTrestr, TDI2FFT, params->minf, __LISASimFD_Noise_fHigh);
      RestrictFDReImFrequencySeries(&TDI3FFTrestr, TDI3FFT, params->minf, __LISASimFD_Noise_fHigh);

      /* Compute SNR with linear integration, weighting with non-rescaled noise functions */
      /* Note: assumes same lengths for all FD FFT frequency series */
      /* Note: non-rescaled noise functions */
      int sizeA = TDI1FFTrestr->freq->size;
      int sizeE = TDI2FFTrestr->freq->size;
      int sizeT = TDI3FFTrestr->freq->size;
      gsl_vector* noisevaluesA = gsl_vector_alloc(sizeA);
      gsl_vector* noisevaluesE = gsl_vector_alloc(sizeE);
      gsl_vector* noisevaluesT = gsl_vector_alloc(sizeT);
      for(int i=0; i<sizeA; i++) {
        gsl_vector_set(noisevaluesA, i, SnAXYZNoRescaling(params->variant, gsl_vector_get(TDI1FFTrestr->freq, i)));
      }
      for(int i=0; i<sizeE; i++) {
        gsl_vector_set(noisevaluesE, i, SnEXYZNoRescaling(params->variant, gsl_vector_get(TDI2FFTrestr->freq, i)));
      }
      for(int i=0; i<sizeT; i++) {
        gsl_vector_set(noisevaluesT, i, SnTXYZNoRescaling(params->variant, gsl_vector_get(TDI3FFTrestr->freq, i)));
      }
      double SNRA2 = FDOverlapReImvsReIm(TDI1FFTrestr, TDI1FFTrestr, noisevaluesA);
      double SNRE2 = FDOverlapReImvsReIm(TDI2FFTrestr, TDI2FFTrestr, noisevaluesE);
      double SNRT2 = FDOverlapReImvsReIm(TDI3FFTrestr, TDI3FFTrestr, noisevaluesT);
      SNR = sqrt(SNRA2 + SNRE2 + SNRT2);

      /* Print SNR to stdout */
      printf("%.8f\n", SNR);
    }

    else {

      /* Initialize the structure for LISAparams and GlobalParams and copy values */
      /* NOTE: injectedparams and globalparams are declared as extern in LISAutils.h, and used by generation functions */
      injectedparams = (LISAParams*) malloc(sizeof(LISAParams));
      memset(injectedparams, 0, sizeof(LISAParams));
      injectedparams->tRef = params->tRef;
      injectedparams->phiRef = params->phiRef;
      injectedparams->m1 = params->m1;
      injectedparams->m2 = params->m2;
      injectedparams->distance = params->distance;
      injectedparams->lambda = params->lambda;
      injectedparams->beta = params->beta;
      injectedparams->inclination = params->inclination;
      injectedparams->polarization = params->polarization;
      injectedparams->nbmode = params->nbmode;
      globalparams = (LISAGlobalParams*) malloc(sizeof(LISAGlobalParams));
      memset(globalparams, 0, sizeof(LISAGlobalParams));
      globalparams->fRef = params->fRef;
      globalparams->deltatobs = params->deltatobs; /* Default value */
      globalparams->minf = params->minf;
      globalparams->maxf = params->maxf;
      globalparams->tagextpn = params->tagextpn;
      globalparams->Mfmatch = params->Mfmatch;
      globalparams->nbmodeinj = params->nbmode;
      globalparams->nbmodetemp = params->nbmode;
      globalparams->tagint = params->tagint;
      globalparams->tagtdi = params->tagtdi;
      globalparams->nbptsoverlap = params->nbptsoverlap;
      globalparams->variant = params->variant;
      globalparams->frozenLISA = params->frozenLISA;
      globalparams->responseapprox = params->responseapprox;
      globalparams->delaycorrection = params->delaycorrection;
      /* Hardcoded */
      globalparams->tagtRefatLISA = tagtRefatLISA;
      globalparams->tagsimplelikelihood22 = tagsimplelikelihood22;
      globalparams->tagsimplelikelihoodHM = tagsimplelikelihoodHM;
      globalparams->zerolikelihood = 0;

      if(params->loadparamsfile==0) {

        /* Branch between the Fresnel or linear computation */
        if(params->tagint==0) {
          LISAInjectionCAmpPhase* injCAmpPhase = NULL;
          LISAInjectionCAmpPhase_Init(&injCAmpPhase);
          LISAGenerateInjectionCAmpPhase(injectedparams, injCAmpPhase);
          SNR = sqrt(injCAmpPhase->TDI123ss);
        }
        else if(params->tagint==1) {
          LISAInjectionReIm* injReIm = NULL;
          LISAInjectionReIm_Init(&injReIm);
          LISAGenerateInjectionReIm(injectedparams, params->minf, params->nbptsoverlap, 0, injReIm); /* Hardcoded linear sampling */

          double SNRA2 = FDOverlapReImvsReIm(injReIm->TDI1Signal, injReIm->TDI1Signal, injReIm->noisevalues1);
          double SNRE2 = FDOverlapReImvsReIm(injReIm->TDI2Signal, injReIm->TDI2Signal, injReIm->noisevalues2);
          double SNRT2 = FDOverlapReImvsReIm(injReIm->TDI3Signal, injReIm->TDI3Signal, injReIm->noisevalues3);
          SNR = sqrt(SNRA2 + SNRE2 + SNRT2);
        }
        else {
          printf("Error in ComputeLISASNR: integration tag not recognized.\n");
          exit(1);
        }

        /* Print SNR to stdout */
        printf("%.8f\n", SNR);
      }
      else {

        int nlines = params->nlinesparams;

        /* Load parameters file */
        /* Format (same as in the internals): m1, m2, tRef, dist, phase, inc, lambda, beta, pol */
        gsl_matrix* inmatrix =  gsl_matrix_alloc(nlines, 9);
        Read_Text_Matrix(params->paramsdir, params->paramsfile, inmatrix);

        /* Initialize output matrix */
        /* Format (same as in the internals): m1, m2, tRef, dist, phase, inc, lambda, beta, pol, SNR */
        gsl_matrix* outmatrix =  gsl_matrix_alloc(nlines, 10);

        for(int i=0; i<nlines; i++) {
          injectedparams->m1 = gsl_matrix_get(inmatrix, i, 0);
          injectedparams->m2 = gsl_matrix_get(inmatrix, i, 1);
          injectedparams->tRef = gsl_matrix_get(inmatrix, i, 2);
          injectedparams->distance = gsl_matrix_get(inmatrix, i, 3);
          injectedparams->phiRef = gsl_matrix_get(inmatrix, i, 4);
          injectedparams->inclination = gsl_matrix_get(inmatrix, i, 5);
          injectedparams->lambda = gsl_matrix_get(inmatrix, i, 6);
          injectedparams->beta = gsl_matrix_get(inmatrix, i, 7);
          injectedparams->polarization = gsl_matrix_get(inmatrix, i, 8);

          /* Branch between the Fresnel or linear computation */
          double SNR = 0;
          if(params->tagint==0) {
            LISAInjectionCAmpPhase* injCAmpPhase = NULL;
            LISAInjectionCAmpPhase_Init(&injCAmpPhase);
            LISAGenerateInjectionCAmpPhase(injectedparams, injCAmpPhase);
            SNR = sqrt(injCAmpPhase->TDI123ss);
          }
          else if(params->tagint==1) {
            LISAInjectionReIm* injReIm = NULL;
            LISAInjectionReIm_Init(&injReIm);
            LISAGenerateInjectionReIm(injectedparams, params->minf, params->nbptsoverlap, 0, injReIm); /* Hardcoded linear sampling */

            double SNRA2 = FDOverlapReImvsReIm(injReIm->TDI1Signal, injReIm->TDI1Signal, injReIm->noisevalues1);
            double SNRE2 = FDOverlapReImvsReIm(injReIm->TDI2Signal, injReIm->TDI2Signal, injReIm->noisevalues2);
            double SNRT2 = FDOverlapReImvsReIm(injReIm->TDI3Signal, injReIm->TDI3Signal, injReIm->noisevalues3);
            SNR = sqrt(SNRA2 + SNRE2 + SNRT2);
          }
          else {
            printf("Error in ComputeLISASNR: integration tag not recognized.\n");
            exit(1);
          }

          /* Set values in output matrix */
          gsl_matrix_set(outmatrix, i, 0, injectedparams->m1);
          gsl_matrix_set(outmatrix, i, 1, injectedparams->m2);
          gsl_matrix_set(outmatrix, i, 2, injectedparams->tRef);
          gsl_matrix_set(outmatrix, i, 3, injectedparams->distance);
          gsl_matrix_set(outmatrix, i, 4, injectedparams->phiRef);
          gsl_matrix_set(outmatrix, i, 5, injectedparams->inclination);
          gsl_matrix_set(outmatrix, i, 6, injectedparams->lambda);
          gsl_matrix_set(outmatrix, i, 7, injectedparams->beta);
          gsl_matrix_set(outmatrix, i, 8, injectedparams->polarization);
          gsl_matrix_set(outmatrix, i, 9, SNR);
        }

        /* Output matrix */
        Write_Text_Matrix(params->paramsdir, params->outputfile, outmatrix);
      }
    }
  }
}
