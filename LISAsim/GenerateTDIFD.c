#include "GenerateTDIFD.h"

/************ Parsing arguments function ************/

/* Function to convert string input TDI string to TDItag */
static GenTDIFDtag ParseGenTDIFDtag(char* string) {
  GenTDIFDtag tag;
  if(strcmp(string, "TDIhlm")==0) tag = TDIhlm;
  else if(strcmp(string, "TDIFD")==0) tag = TDIFD;
  else if(strcmp(string, "h22FFT")==0) tag = h22FFT;
  else if(strcmp(string, "yslrFFT")==0) tag = yslrFFT;
  else {
    printf("Error in ParseGenTDIFDtag: string not recognized.\n");
    exit(1);
  }
  return tag;
}

/* Parse command line to initialize GenTDIFDparams object */
/* Masses are input in solar masses and distances in Mpc - converted in SI for the internals */
static void parse_args_GenerateTDIFD(ssize_t argc, char **argv, GenTDIFDparams* params)
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
 --minf                Minimal frequency (Hz, default=0) - when set to 0, use the lowest frequency where the detector noise model is trusted __LISASimFD_Noise_fLow (set somewhat arbitrarily)\n\
 --maxf                Maximal frequency (Hz, default=0) - when set to 0, use the highest frequency where the detector noise model is trusted __LISASimFD_Noise_fHigh (set somewhat arbitrarily)\n\
 --deltatobs           Observation duration (years, default=2)\n\
 --tagextpn            Tag to allow PN extension of the waveform at low frequencies (default=1)\n\
 --Mfmatch             When PN extension allowed, geometric matching frequency: will use ROM above this value. If <=0, use ROM down to the lowest covered frequency (default=0.)\n\
 --deltaf              When generating frequency series from the mode contributions, deltaf for the output (0 to set automatically at 1/2*1/(2T))\n\
 --twindowbeg          When generating frequency series from file by FFT, twindowbeg (0 to set automatically at 0.05*duration)\n\
 --twindowend          When generating frequency series from file by FFT, twindowend (0 to set automatically at 0.01*duration)\n\
 --tagtdi              Tag choosing the set of TDI variables to use (default TDIAETXYZ)\n\
 --taggenwave          Tag choosing the wf format: TDIhlm (default: downsampled mode contributions to TDI in Amp/Phase form), TDIFD (hlm interpolated and summed accross modes), h22FFT and yslrFFT (used with fomrtditdfile -- loads time series and FFT)\n\
 --restorescaledfactor Option to restore the factors scaled out of TDI observables (default: false)\n	\
 --fromtdfile          Option for loading time series and FFTing (default: false)\n\
 --nsamplesinfile      Number of lines of input file when loading TDI time series from file\n\
 --binaryin            Tag for loading the data in gsl binary form instead of text (default false)\n\
 --binaryout           Tag for outputting the data in gsl binary form instead of text (default false)\n\
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
    params->minf = 0.;
    params->maxf = 0.;
    params->deltatobs = 2.;
    params->tagextpn = 1;
    params->Mfmatch = 0.;
    params->deltaf = 0.;
    params->twindowbeg = 0.;
    params->twindowend = 0.;
    params->tagtdi = TDIXYZ;
    params->taggenwave = TDIhlm;
    params->restorescaledfactor = 0;
    params->fromtdfile = 0;
    params->nsamplesinfile = 0;    /* No default; has to be provided */
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
        } else if (strcmp(argv[i], "--minf") == 0) {
            params->minf = atof(argv[++i]);
        } else if (strcmp(argv[i], "--maxf") == 0) {
            params->maxf = atof(argv[++i]);
        } else if (strcmp(argv[i], "--deltatobs") == 0) {
            params->deltatobs = atof(argv[++i]);
        }  else if (strcmp(argv[i], "--tagextpn") == 0) {
            params->tagextpn = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--Mfmatch") == 0) {
            params->Mfmatch = atof(argv[++i]);
        } else if (strcmp(argv[i], "--deltaf") == 0) {
            params->deltaf = atof(argv[++i]);
        } else if (strcmp(argv[i], "--twindowbeg") == 0) {
            params->twindowbeg = atof(argv[++i]);
        } else if (strcmp(argv[i], "--twindowend") == 0) {
            params->twindowend = atof(argv[++i]);
        } else if (strcmp(argv[i], "--tagtdi") == 0) {
	  params->tagtdi = ParseTDItag(argv[++i]);
        } else if (strcmp(argv[i], "--taggenwave") == 0) {
	  params->taggenwave = ParseGenTDIFDtag(argv[++i]);
        } else if (strcmp(argv[i], "--restorescaledfactor") == 0) {
            params->restorescaledfactor = 1;
        } else if (strcmp(argv[i], "--fromtdfile") == 0) {
            params->fromtdfile = 1;
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

    /* Set frequency interval to default values */
    if(params->minf==0.) params->minf = __LISASimFD_Noise_fLow;
    if(params->maxf==0.) params->maxf = __LISASimFD_Noise_fHigh;

    return;

 fail:
    exit(1);
}

/* Built-in postfix for output files for TDI channels */
static char* TDIFilePostfix(TDItag tditag, int channel, int binary)
{
  if(!binary) {
    if(tditag==y12) {
      switch(channel) {
        case 1: return "_y12.txt"; break;
        case 2: return "_null2.txt"; break;
        case 3: return "_null3.txt"; break;
      }
    }
    else if(tditag==TDIXYZ) {
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
  else {
    if(tditag==y12) {
      switch(channel) {
        case 1: return "_y12.dat"; break;
        case 2: return "_null2.dat"; break;
        case 3: return "_null3.dat"; break;
      }
    }
    else if(tditag==TDIXYZ) {
      switch(channel) {
        case 1: return "_TDIX.dat"; break;
        case 2: return "_TDIY.dat"; break;
        case 3: return "_TDIZ.dat"; break;
      }
    }
    else if(tditag==TDIAETXYZ) {
      switch(channel) {
        case 1: return "_TDIAXYZ.dat"; break;
        case 2: return "_TDIEXYZ.dat"; break;
        case 3: return "_TDITXYZ.dat"; break;
      }
    }
    else {
      printf("Error: in TDIFilePostfix TDI tag not recognized.\n");
      exit(1);
    }
  }
}

/************ Functions to write waveforms to file ************/

/* Read waveform time series in Re/Im form for hpTD and hcTD a single file */
/* NOTE: assumes the same number of points is used to represent each mode */
static void Read_TDITD3Chan( RealTimeSeries** TDI1, RealTimeSeries** TDI2, RealTimeSeries** TDI3, const char dir[], const char file[], const int nblines, const int binary)
{
  /* Initalize and read input */
  gsl_matrix* inmatrix =  gsl_matrix_alloc(nblines, 4);
  if(!binary) Read_Text_Matrix(dir, file, inmatrix);
  else Read_Matrix(dir, file, inmatrix);

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
/* Output TDI mode contributions in downsampled form, FD AmpReal/AmpIm/Phase, all hlm modes in a single file */
/* NOTE: assumes the same number of points is used to represent each mode */
static void Write_TDIhlm(const char dir[], const char file[], ListmodesCAmpPhaseFrequencySeries* listhlm, int nbmodes, const int binary)
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
  GenTDIFDparams* params;
  params = (GenTDIFDparams*) malloc(sizeof(GenTDIFDparams));
  memset(params, 0, sizeof(GenTDIFDparams));

  /* Parse commandline to read parameters */
  parse_args_GenerateTDIFD(argc, argv, params);

  if(params->fromtdfile) {

    if(params->taggenwave==TDIFD) {
      /* Load TD TDI from file */
      RealTimeSeries* TDI1 = NULL;
      RealTimeSeries* TDI2 = NULL;
      RealTimeSeries* TDI3 = NULL;
      Read_TDITD3Chan(&TDI1, &TDI2, &TDI3, params->indir, params->infile, params->nsamplesinfile, params->binaryin);

      /* Compute FFT */
      ReImFrequencySeries* TDI1FFT = NULL;
      ReImFrequencySeries* TDI2FFT = NULL;
      ReImFrequencySeries* TDI3FFT = NULL;
      gsl_vector* times = TDI1->times;
      double twindowbeg, twindowend;
      if(params->twindowbeg==0.) {
        twindowbeg = 0.05 * (gsl_vector_get(times, times->size - 1) - gsl_vector_get(times, 0)); /* Here hardcoded relative window lengths */
      }
      else {
        twindowbeg = params->twindowbeg;
      }
      if(params->twindowend==0.) {
        twindowend = 0.01 * (gsl_vector_get(times, times->size - 1) - gsl_vector_get(times, 0)); /* Here hardcoded relative window lengths */
      }
      else {
        twindowend = params->twindowend;
      }
      FFTRealTimeSeries(&TDI1FFT, TDI1, twindowbeg, twindowend, 2); /* Here hardcoded 0-padding */
      FFTRealTimeSeries(&TDI2FFT, TDI2, twindowbeg, twindowend, 2); /* Here hardcoded 0-padding */
      FFTRealTimeSeries(&TDI3FFT, TDI3, twindowbeg, twindowend, 2); /* Here hardcoded 0-padding */

      /* Output */
      char *outfileTDI1 = malloc(256);
      char *outfileTDI2 = malloc(256);
      char *outfileTDI3 = malloc(256);
      sprintf(outfileTDI1, "%s%s", params->outfileprefix, TDIFilePostfix(params->tagtdi, 1, params->binaryout));
      Write_ReImFrequencySeries(params->outdir, outfileTDI1, TDI1FFT, params->binaryout);
      if(!(params->tagtdi==y12)) {
        sprintf(outfileTDI2, "%s%s", params->outfileprefix, TDIFilePostfix(params->tagtdi, 2, params->binaryout));
        sprintf(outfileTDI3, "%s%s", params->outfileprefix, TDIFilePostfix(params->tagtdi, 3, params->binaryout));
        Write_ReImFrequencySeries(params->outdir, outfileTDI2, TDI2FFT, params->binaryout);
        Write_ReImFrequencySeries(params->outdir, outfileTDI3, TDI3FFT, params->binaryout);
      }
      free(outfileTDI1);
      free(outfileTDI2);
      free(outfileTDI3);

      exit(0);
    }
    else if(params->taggenwave==h22FFT) {
      /* Load Amp/Phase h22 from file */
      AmpPhaseTimeSeries* h22tdampphase = NULL;

      //
      printf("aaa\n");
      
      Read_AmpPhaseTimeSeries(&h22tdampphase, params->indir, params->infile, params->nsamplesinfile, params->binaryin);

      /* Convert to Re/Im form */
      ReImTimeSeries* h22td = NULL;
      AmpPhaseTimeSeries_ToReIm(&h22td, h22tdampphase);

//
printf("aaa\n");

      /* Compute FFT */
      ReImFrequencySeries* h22FFT = NULL;
      gsl_vector* times = h22td->times;
      double twindowbeg, twindowend;
      if(params->twindowbeg==0.) {
        twindowbeg = 0.05 * (gsl_vector_get(times, times->size - 1) - gsl_vector_get(times, 0)); /* Here hardcoded relative window lengths */
      }
      else {
        twindowbeg = params->twindowbeg;
      }
      if(params->twindowend==0.) {
        twindowend = 0.01 * (gsl_vector_get(times, times->size - 1) - gsl_vector_get(times, 0)); /* Here hardcoded relative window lengths */
      }
      else {
        twindowend = params->twindowend;
      }
      //
      printf("aaa\n");

      FFTTimeSeries(&h22FFT, h22td, twindowbeg, twindowend, 2); /* Here hardcoded 0-padding */

      //
      printf("aaa\n");

      /* Output */
      Write_ReImFrequencySeries(params->outdir, params->outfileprefix, h22FFT, params->binaryout);
    }
    else if(params->taggenwave==yslrFFT) {
      //
      printf("a\n");

      /* Load Re/Im h22 from file */
      RealTimeSeries* yslrtd = NULL;
      Read_RealTimeSeries(&yslrtd, params->indir, params->infile, params->nsamplesinfile, params->binaryin);

      //
      printf("a\n");
      printf("times[0,1] %g, %g\n", gsl_vector_get(yslrtd->times, 0), gsl_vector_get(yslrtd->times, 1));

      /* Compute FFT */
      ReImFrequencySeries* yslrFFT = NULL;
      gsl_vector* times = yslrtd->times;
      double twindowbeg, twindowend;
      if(params->twindowbeg==0.) {
        twindowbeg = 0.05 * (gsl_vector_get(times, times->size - 1) - gsl_vector_get(times, 0)); /* Here hardcoded relative window lengths */
      }
      else {
        twindowbeg = params->twindowbeg;
      }
      if(params->twindowend==0.) {
        twindowend = 0.01 * (gsl_vector_get(times, times->size - 1) - gsl_vector_get(times, 0)); /* Here hardcoded relative window lengths */
      }
      else {
        twindowend = params->twindowend;
      }

        //
        printf("aa\n");
      FFTRealTimeSeries(&yslrFFT, yslrtd, twindowbeg, twindowend, 2); /* Here hardcoded 0-padding */

      /* Output */
      Write_ReImFrequencySeries(params->outdir, params->outfileprefix, yslrFFT, params->binaryout);
    }
  }

  else {
    /* Set geometric coefficients */
    SetCoeffsG(params->lambda, params->beta, params->polarization);

    /* Starting frequency - takes into account the duration of observation deltatobs and the arg minf */
    double fstartobs = 0.;
    if(!(params->deltatobs==0)) fstartobs = Newtonianfoft(params->m1, params->m2, params->deltatobs);
    double fLow = fmax(params->minf, fstartobs);

    /* Generate Fourier-domain waveform as a list of hlm modes */
    /* Use TF2 extension, if required to, to arbitrarily low frequencies */
    /* NOTE: at this stage, if no extension is performed, deltatobs and minf are ignored - will start at MfROM*/
    ListmodesCAmpPhaseFrequencySeries* listROM = NULL;
    if(!(params->tagextpn)){
      //printf("Not Extending signal waveform.  mfmatch=%g\n",globalparams->mfmatch);
      SimEOBNRv2HMROM(&listROM, params->nbmode, params->tRef, params->phiRef, params->fRef, (params->m1)*MSUN_SI, (params->m2)*MSUN_SI, (params->distance)*1e6*PC_SI);
    } else {
      //printf("Extending signal waveform.  mfmatch=%g\n",globalparams->mfmatch);
      SimEOBNRv2HMROMExtTF2(&listROM, params->nbmode, params->Mfmatch, fLow, 0, params->tRef, params->phiRef, params->fRef, (params->m1)*MSUN_SI, (params->m2)*MSUN_SI, (params->distance)*1e6*PC_SI);
    }

    /* Process through the Fourier-domain response for TDI observables */
    ListmodesCAmpPhaseFrequencySeries* listTDI1= NULL;
    ListmodesCAmpPhaseFrequencySeries* listTDI2= NULL;
    ListmodesCAmpPhaseFrequencySeries* listTDI3= NULL;
    LISASimFDResponseTDI3Chan(&listROM, &listTDI1, &listTDI2, &listTDI3, params->tRef, params->lambda, params->beta, params->inclination, params->polarization, params->maxf, params->tagtdi);

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
      sprintf(outfileTDI1hlm, "%s%s", params->outfileprefix, TDIFilePostfix(params->tagtdi, 1, params->binaryout));
      Write_TDIhlm(params->outdir, outfileTDI1hlm, listTDI1, params->nbmode, params->binaryout);
      if(!(params->tagtdi==y12)) {
        sprintf(outfileTDI2hlm, "%s%s", params->outfileprefix, TDIFilePostfix(params->tagtdi, 2, params->binaryout));
        sprintf(outfileTDI3hlm, "%s%s", params->outfileprefix, TDIFilePostfix(params->tagtdi, 3, params->binaryout));
        Write_TDIhlm(params->outdir, outfileTDI2hlm, listTDI2, params->nbmode, params->binaryout);
        Write_TDIhlm(params->outdir, outfileTDI3hlm, listTDI3, params->nbmode, params->binaryout);
      }
      free(outfileTDI1hlm);
      free(outfileTDI2hlm);
      free(outfileTDI3hlm);

      exit(0);
    }

    else {

      double deltaf = 0.;
      if(params->deltaf==0.) {
        /* Determine deltaf so that N deltat = 1/deltaf > 2*tc where tc is the time to coalescence estimated from Psi22 */
        /* NOTE: assumes the TD waveform ends around t=0 */
        /* NOTE: estimate based on the 22 mode - fstartobs will be scaled from mode to mode to ensure the same deltatobs for all (except for the 21 mode, which will turn on after the others) */
        double tc = EstimateInitialTime(listROM, fLow);
        double deltaf = 0.5 * 1./(2*(-tc)); /* Extra factor of 1/2 corresponding to 0-padding in TD by factor of 2 */
      }
      else deltaf = params->deltaf;

      /* Compute TDI FD */
      /* NOTE: we do NOT use fLow in the role of fstartobs - even when not asking for a deltatobs (deltatobs=0, fstartobs=0), we could enforce that the different modes will start at a frequency that corresponds to the same starting time (again, except for the 21 mode which will turn on after the start of the waveform). We do not do that, as it is primarily useful if we plan to take an IFFT afterwards. */

      /* Generate FD frequency series by summing mode contributions */
      ReImFrequencySeries* TDI1FD = NULL;
      ReImFrequencySeries* TDI2FD = NULL;
      ReImFrequencySeries* TDI3FD = NULL;
      GenerateFDReImFrequencySeries(&TDI1FD, listTDI1, fLow, params->maxf, fstartobs, deltaf, 0); /* Here determines the number of points from deltaf and max frequency in input list of modes */
      GenerateFDReImFrequencySeries(&TDI2FD, listTDI2, fLow, params->maxf, fstartobs, deltaf, 0); /* Here determines the number of points from deltaf and max frequency in input list of modes */
      GenerateFDReImFrequencySeries(&TDI3FD, listTDI3, fLow, params->maxf, fstartobs, deltaf, 0); /* Here determines the number of points from deltaf and max frequency in input list of modes */

      /* Output */
      char *outfileTDI1FD = malloc(256);
      char *outfileTDI2FD = malloc(256);
      char *outfileTDI3FD = malloc(256);
      sprintf(outfileTDI1FD, "%s%s", params->outfileprefix, TDIFilePostfix(params->tagtdi, 1, params->binaryout));
      Write_ReImFrequencySeries(params->outdir, outfileTDI1FD, TDI1FD, params->binaryout);
      if(!(params->tagtdi==y12)) {
        sprintf(outfileTDI2FD, "%s%s", params->outfileprefix, TDIFilePostfix(params->tagtdi, 2, params->binaryout));
        sprintf(outfileTDI3FD, "%s%s", params->outfileprefix, TDIFilePostfix(params->tagtdi, 3, params->binaryout));
        Write_ReImFrequencySeries(params->outdir, outfileTDI2FD, TDI2FD, params->binaryout);
        Write_ReImFrequencySeries(params->outdir, outfileTDI3FD, TDI3FD, params->binaryout);
      }
      free(outfileTDI1FD);
      free(outfileTDI2FD);
      free(outfileTDI3FD);

      exit(0);
    }
  }
}
