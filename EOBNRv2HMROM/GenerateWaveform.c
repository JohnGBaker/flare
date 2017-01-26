#include "GenerateWaveform.h"

/************ Parsing arguments function ************/

/* Function to convert string input TDI string to TDItag */
static GenWavetag ParseGenWavetag(char* string) {
  GenWavetag tag;
  if(strcmp(string, "hlm")==0) tag = hlm;
  else if(strcmp(string, "h22TD")==0) tag = h22TD;
  else if(strcmp(string, "hphcFD")==0) tag = hphcFD;
  else if(strcmp(string, "hphcTD")==0) tag = hphcTD;
  else {
    printf("Error in ParseGenWavetag: string not recognized.\n");
    exit(1);
  }
  return tag;
}

/* Parse command line to initialize GenWaveParams object */
/* Masses are input in solar masses and distances in Mpc - converted in SI for the internals */
static void parse_args_GenerateWaveform(ssize_t argc, char **argv, GenWaveParams* params)
{
    char help[] = "\
GenerateWaveform by Sylvain Marsat, John Baker, and Philip Graff\n\
Copyright July 2015\n\
\n\
This program generates and outputs a waveform produced with EOBNRv2HMROM.\n\
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
\n\
--------------------------------------------------\n\
----- Generation Parameters ----------------------\n\
--------------------------------------------------\n\
 --nbmode              Number of modes of radiation to generate (1-5, default=1)\n\
 --minf                Minimal frequency, ignore if 0 (Hz, default=0) - will use first frequency covered by the ROM if higher\n\
 --maxf                Maximal frequency, ignore if 0 (Hz, default=0) - will use last frequency covered by the ROM if lower\n\
 --deltatobs           Observation duration (years) ignore if 0 (default=2)\n\
 --tagextpn            Tag to allow PN extension of the waveform at low frequencies (default=1)\n\
 --Mfmatch             When PN extension allowed, geometric matching frequency: will use ROM above this value. If <=0, use ROM down to the lowest covered frequency (default=0.)\n\
 --taggenwave          Tag choosing the wf format: hlm (default: downsampled modes, Amp/Phase form),  h22TD (IFFT of h22, Amp/Phase form - currently not supported for higher modes), hphcFD (hlm interpolated and summed, Re/Im form), hphcTD (IFFT of hphcFD, Re/Im form)\n\
 --f1windowbeg         If generating h22TD/hphcTD, start frequency for windowing at the beginning - set to 0 to ignore and use max(fstartobs, fLowROM, minf), where fLowROM is either the lowest frequency covered by the ROM or simply minf if PN extension is used (Hz, default=0)\n\
 --f2windowbeg         If generating h22TD/hphcTD, stop frequency for windowing at the beginning - set to 0 to ignore and use 1.1*f1windowbeg (Hz, default=0)\n\
 --f1windowend         If generating h22TD/hphcTD, start frequency for windowing at the end - set to 0 to ignore and use 0.995*f2windowend (Hz, default=0)\n\
 --f2windowend         If generating h22TD/hphcTD, stop frequency for windowing at the end - set to 0 to ignore and use min(maxf, fHighROM), where fHighROM is the highest frequency covered by the ROM (Hz, default=0)\n\
 --tagh22fromfile      Tag choosing wether to load h22 FD downsampled Amp/Phase from file (default false) - NOTE: for now low-f cut depends on fstartobs that still depends on m1, m2- so one has to re-pass m1 and m2 !\n\
 --nsamplesinfile      Number of lines of inputs file\n\
 --binaryin            Tag for loading the data in gsl binary form instead of text (default false)\n\
 --indir               Input directory\n\
 --infile              Input file name\n\
 --binaryout           Tag for writnig the data in gsl binary form instead of text (default false)\n\
 --outdir              Output directory\n\
 --outfile             Output file name\n\
\n";

    ssize_t i;

    /* set default values for the physical params */
    params->tRef = 0.;
    params->phiRef = 0.;
    params->fRef = 0.;
    params->m1 = 2*1e6;
    params->m2 = 1*1e6;
    params->distance = 1e3;
    params->inclination = PI/3;

    /* set default values for the generation params */
    params->nbmode = 1;
    params->minf = 0.;
    params->maxf = 0.;
    params->deltatobs = 0.;
    params->tagextpn = 1;
    params->Mfmatch = 0.;
    params->taggenwave = hlm;
    params->f1windowbeg = 0.;
    params->f2windowbeg = 0.;
    params->f1windowend = 0.;
    params->f2windowend = 0.;
    params->tagh22fromfile = 0;
    params->nsamplesinfile = 0; /* No default; has to be provided */
    params->binaryin = 0;
    strcpy(params->indir, ""); /* No default; has to be provided */
    strcpy(params->infile, ""); /* No default; has to be provided */
    params->binaryout = 0;
    strcpy(params->outdir, ".");
    strcpy(params->outfile, "generated_waveform.txt");

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
        } else if (strcmp(argv[i], "--nbmode") == 0) {
          params->nbmode = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--minf") == 0) {
          params->minf = atof(argv[++i]);
        } else if (strcmp(argv[i], "--maxf") == 0) {
          params->maxf = atof(argv[++i]);
        } else if (strcmp(argv[i], "--deltatobs") == 0) {
          params->deltatobs = atof(argv[++i]);
        } else if (strcmp(argv[i], "--tagextpn") == 0) {
          params->tagextpn = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--Mfmatch") == 0) {
          params->Mfmatch = atof(argv[++i]);
        } else if (strcmp(argv[i], "--taggenwave") == 0) {
          params->taggenwave = ParseGenWavetag(argv[++i]);
        } else if (strcmp(argv[i], "--f1windowbeg") == 0) {
          params->f1windowbeg = atof(argv[++i]);
        } else if (strcmp(argv[i], "--f2windowbeg") == 0) {
          params->f2windowbeg = atof(argv[++i]);
        }  else if (strcmp(argv[i], "--f1windowend") == 0) {
          params->f1windowend = atof(argv[++i]);
        } else if (strcmp(argv[i], "--f2windowend") == 0) {
          params->f2windowend = atof(argv[++i]);
        } else if (strcmp(argv[i], "--tagh22fromfile") == 0) {
          params->tagh22fromfile = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--nsamplesinfile") == 0) {
          params->nsamplesinfile = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--binaryin") == 0) {
          params->binaryin = 1;
        } else if (strcmp(argv[i], "--indir") == 0) {
            strcpy(params->indir, argv[++i]);
        } else if (strcmp(argv[i], "--infile") == 0) {
            strcpy(params->infile, argv[++i]);
        } else if (strcmp(argv[i], "--binaryout") == 0) {
          params->binaryout = 1;
        } else if (strcmp(argv[i], "--outdir") == 0) {
            strcpy(params->outdir, argv[++i]);
        } else if (strcmp(argv[i], "--outfile") == 0) {
            strcpy(params->outfile, argv[++i]);
        } else {
	  printf("Error: invalid option: %s\n", argv[i]);
	  goto fail;
        }
    }

    return;

 fail:
    exit(1);
}

/************ Functions to write waveforms to file - text or binary ************/

/* Input waveform in downsampled form, FD Amp/Pase - for now supports only single-mode h22 waveform */
/* NOTE: assumes Amp/Phase format in the file, i.e. no complex amplitude */
static void Read_Wave_h22(const char dir[], const char file[], ListmodesCAmpPhaseFrequencySeries** listhlm, int nsamples, int binary)
{
  /* Initalize and read input */
  gsl_matrix* inmatrix =  gsl_matrix_alloc(nsamples, 3);
  if (!binary) Read_Text_Matrix(dir, file, inmatrix);
  else Read_Matrix(dir, file, inmatrix);
  gsl_vector_view freqview = gsl_matrix_column(inmatrix, 0);
  gsl_vector_view ampview = gsl_matrix_column(inmatrix, 1);
  gsl_vector_view phaseview = gsl_matrix_column(inmatrix, 2);

  /* Copy in frequency series */
  CAmpPhaseFrequencySeries* h22 = NULL;
  CAmpPhaseFrequencySeries_Init(&h22, nsamples);
  gsl_vector_memcpy(h22->freq, &freqview.vector);
  gsl_vector_memcpy(h22->amp_real, &ampview.vector);
  gsl_vector_memcpy(h22->phase, &phaseview.vector);
  gsl_vector_set_zero(h22->amp_imag);

  /* Output */
  (*listhlm) = ListmodesCAmpPhaseFrequencySeries_AddModeNoCopy((*listhlm), h22, 2, 2);
}
/* Output waveform in downsampled form, FD Amp/Pase, all hlm modes in a single file */
/* NOTE: first version assumed the same number of points is used to represent each mode */
/* HACK: if not same number of freqs for the modes (due to PN extension), 0-pad at the end */
static void Write_Wave_hlm(const char dir[], const char file[], ListmodesCAmpPhaseFrequencySeries* listhlm, int nbmodes, int binary)
{
  /* Initialize output */
  /* get length from 22 mode - NOTE: assumes the same for all modes */
  int nbfreq = ListmodesCAmpPhaseFrequencySeries_GetMode(listhlm, 2, 2)->freqseries->freq->size;
  /* HACK: get the maximal length */
  for(int i=0; i<nbmodes; i++) {
    nbfreq = max(nbfreq, ListmodesCAmpPhaseFrequencySeries_GetMode(listhlm, listmode[i][0], listmode[i][1])->freqseries->freq->size);
  }
  gsl_matrix* outmatrix = gsl_matrix_alloc(nbfreq, 3*nbmodes);

  /* Get data in the list of modes */
  CAmpPhaseFrequencySeries* mode;
  for(int i=0; i<nbmodes; i++) {
    mode = ListmodesCAmpPhaseFrequencySeries_GetMode(listhlm, listmode[i][0], listmode[i][1])->freqseries;
    int nbfreqmode = mode->freq->size;
    for(int j=0; j<nbfreqmode; j++) {
      gsl_matrix_set(outmatrix, j, 0+3*i, gsl_vector_get(mode->freq, j));
      gsl_matrix_set(outmatrix, j, 1+3*i, gsl_vector_get(mode->amp_real, j)); /* amp_imag is 0 at this stage, we ignore it */
      gsl_matrix_set(outmatrix, j, 2+3*i, gsl_vector_get(mode->phase, j));
    }
    for(int j=nbfreqmode; j<nbfreq; j++) {
      gsl_matrix_set(outmatrix, j, 0+3*i, 0);
      gsl_matrix_set(outmatrix, j, 1+3*i, 0);
      gsl_matrix_set(outmatrix, j, 2+3*i, 0);
    }
  }

  /* Output */
  if (!binary) Write_Text_Matrix(dir, file, outmatrix);
  else Write_Matrix(dir, file, outmatrix);
}
/* Output waveform in frequency series form,Re/Im for hplus and hcross */
static void Write_Wave_hphcFD(const char dir[], const char file[], ReImFrequencySeries* hptilde, ReImFrequencySeries* hctilde, int binary)
{
  /* Initialize output */
  /* Note: assumes hplus, hcross have same length as expected */
  int nbfreq = hptilde->freq->size;
  gsl_matrix* outmatrix = gsl_matrix_alloc(nbfreq, 5);

  /* Set output matrix */
  gsl_matrix_set_col(outmatrix, 0, hptilde->freq);
  gsl_matrix_set_col(outmatrix, 1, hptilde->h_real);
  gsl_matrix_set_col(outmatrix, 2, hptilde->h_imag);
  gsl_matrix_set_col(outmatrix, 3, hctilde->h_real);
  gsl_matrix_set_col(outmatrix, 4, hctilde->h_imag);

  /* Output */
  if (!binary) Write_Text_Matrix(dir, file, outmatrix);
  else Write_Matrix(dir, file, outmatrix);
}
/* Output waveform in frequency series form,Re/Im for hplus and hcross */
static void Write_Wave_hphcTD(const char dir[], const char file[], RealTimeSeries* hp, RealTimeSeries* hc, int binary)
{
  /* Initialize output */
  /* Note: assumes hplus, hcross have same length as expected */
  int nbtimes = hp->times->size;
  gsl_matrix* outmatrix = gsl_matrix_alloc(nbtimes, 3);

  /* Set output matrix */
  gsl_matrix_set_col(outmatrix, 0, hp->times);
  gsl_matrix_set_col(outmatrix, 1, hp->h);
  gsl_matrix_set_col(outmatrix, 2, hc->h);

  /* Output */
  if (!binary) Write_Text_Matrix(dir, file, outmatrix);
  else Write_Matrix(dir, file, outmatrix);
}
/* Output waveform in frequency series form, Re/Im for hplus and hcross */
    // Output Re/Im because difference between numpy and C unwrapping
static void Write_Wave_h22TD(const char dir[], const char file[], AmpPhaseTimeSeries* h22, int binary)
{
  /* Initialize output */
  int nbtimes = h22->times->size;
  gsl_matrix* outmatrix = gsl_matrix_alloc(nbtimes, 3);

  /* Set output matrix */
  gsl_matrix_set_col(outmatrix, 0, h22->times);
  gsl_matrix_set_col(outmatrix, 1, h22->h_amp);
  gsl_matrix_set_col(outmatrix, 2, h22->h_phase);

  /* Output */
  if (!binary) Write_Text_Matrix(dir, file, outmatrix);
  else Write_Matrix(dir, file, outmatrix);
}

/***************** Main program *****************/

int main(int argc, char *argv[])
{
  int ret;
  /* Initialize structure for parameters */
  GenWaveParams* params;
  params = (GenWaveParams*) malloc(sizeof(GenWaveParams));
  memset(params, 0, sizeof(GenWaveParams));

  /* Parse commandline to read parameters */
  parse_args_GenerateWaveform(argc, argv, params);

  /* Generate waveform hlm FD, downsampled in amp/phase form */
  double fLowhlm = 0.;
  double fHighhlm = 0.;
  double fstartobs = 0.;
  ListmodesCAmpPhaseFrequencySeries* listhlm = NULL;
  /* If not loading h22 DS amp/phase from file, generate hlm FD */
  if(!(params->tagh22fromfile)) {
    /* Starting frequency corresponding to duration of observation deltatobs */
    if(!(params->deltatobs==0.)) fstartobs = Newtonianfoft(params->m1, params->m2, params->deltatobs);
      
    /* Generate Fourier-domain waveform as a list of hlm modes */
    /* Use TF2 extension, if required to, to arbitrarily low frequencies */
    /* NOTE: at this stage, if no extension is performed, deltatobs and minf are ignored - will start at MfROM */
    /* If extending, taking into account both fstartobs and minf */
    if(!(params->tagextpn)){
      //printf("Not Extending signal waveform.  mfmatch=%g\n",globalparams->mfmatch);
      ret = SimEOBNRv2HMROM(&listhlm, params->nbmode, params->tRef, params->phiRef, params->fRef, (params->m1)*MSUN_SI, (params->m2)*MSUN_SI, (params->distance)*1e6*PC_SI);
    } else {
      //printf("Extending signal waveform.  mfmatch=%g\n",globalparams->mfmatch);
      ret = SimEOBNRv2HMROMExtTF2(&listhlm, params->nbmode, params->Mfmatch, fmax(params->minf, fstartobs), 0, params->tRef, params->phiRef, params->fRef, (params->m1)*MSUN_SI, (params->m2)*MSUN_SI, (params->distance)*1e6*PC_SI);
    }
  }
  /* Read h22 from file - here fstartobs is ignored, use the starting frequency in the file */
  else {
    if(!(params->deltatobs)==0) {
      printf("Warning in GenerateWaveform: loading h22 FD from file, ignoring deltatobs.\n");
    }
    Read_Wave_h22(params->indir, params->infile, &listhlm, params->nsamplesinfile, params->binaryin);
  }

  /* Determine highest and lowest frequency to cover */
  /* Takes into account limited duration of obsevation, frequencies covered by the ROM and input-given minf, maxf */
  fLowhlm = ListmodesCAmpPhaseFrequencySeries_minf(listhlm);
  fHighhlm = ListmodesCAmpPhaseFrequencySeries_maxf(listhlm);
  /* If we read h22 from file, deltatobs is ignored and fstartobs is simply 0 */
  double fLow = fmax(fLowhlm, fmax(params->minf, fstartobs));
  double fHigh = fmin(fHighhlm, params->maxf);

  /* Window frequencies - used in case of an IFFT */
  double f1windowbeg = 0.;
  double f2windowbeg = 0.;
  double f1windowend = 0.;
  double f2windowend = 0.;

  if(params->taggenwave==hlm) {
    Write_Wave_hlm(params->outdir, params->outfile, listhlm, params->nbmode, params->binaryout);
    exit(0);
  }

  else if((params->taggenwave==hphcFD) || (params->taggenwave==hphcTD)) {

    /* Determine deltaf so that N deltat = 1/deltaf > 2*tc where tc is the time to coalescence estimated from Psi22 */
    /* NOTE: assumes the TD waveform ends around t=0 */
    /* NOTE: estimate based on the 22 mode - fstartobs will be scaled from mode to mode to ensure the same deltatobs for all (except for the 21 mode, which will turn on after the others) */
    double tc = EstimateInitialTime(listhlm, fLow);

    double deltaf = 0.5 * 1./(2*(-tc)); /* Extra factor of 1/2 corresponding to 0-padding in TD by factor of 2 */

    /* Compute hplus, hcross FD */
    /* NOTE: we use fLow in the role of fstartobs - even when not asking for a deltatobs (deltatobs=0, fstartobs=0), the different modes will start at a frequency that correspond to the same starting time (again, except for the 21 mode which will turn on after the start of the waveform) */
    ReImFrequencySeries* hptilde = NULL;
    ReImFrequencySeries* hctilde = NULL;
    GeneratehphcFDReImFrequencySeries(&hptilde, &hctilde, listhlm, params->minf, params->maxf, fLow, deltaf, 0, params->nbmode, params->inclination, 0., 1);

    /* Output */
    if(params->taggenwave==hphcFD) {

      Write_Wave_hphcFD(params->outdir, params->outfile, hptilde, hctilde, params->binaryout);
      exit(0);
    }

    else {

      /* Determine frequency windows */
      if(params->f1windowbeg==0.) f1windowbeg = fLow;
      else f1windowbeg = params->f1windowbeg;
      if(params->f2windowbeg==0.) f2windowbeg = 1.1*f1windowbeg;
      else f2windowbeg = params->f2windowbeg;
      if(params->f2windowend==0.) f2windowend = fHigh;
      else f2windowend = params->f2windowend;
      if(params->f1windowend==0.) f1windowend = 0.995*f2windowend;
      else f1windowend = params->f1windowend;

      /* Compute hplus, hcross TD */
      RealTimeSeries* hp = NULL;
      RealTimeSeries* hc = NULL;
      IFFTFrequencySeriesReal(&hp, hptilde, f1windowbeg, f2windowbeg, f1windowend, f2windowend, 3); /* Here hardcoded nzeropadding */
      IFFTFrequencySeriesReal(&hc, hctilde, f1windowbeg, f2windowbeg, f1windowend, f2windowend, 3); /* Here hardcoded nzeropadding */

      /* Output */
      if(params->taggenwave==hphcTD) {
        Write_Wave_hphcTD(params->outdir, params->outfile, hp, hc, params->binaryout);
        exit(0);
      }
    }
  }
  else if(params->taggenwave==h22TD) {

    /* Check nbmode */
    if(!(params->nbmode==1)) {
      printf("Error in GenerateWaveform: for taggenwave=h22TD, nbmode must be 1 (22-only waveform).\n");
      exit(1);
    }

    /* Determine deltaf so that N deltat = 1/deltaf > 2*tc where tc is the time to coalescence estimated from Psi22 */
    /* NOTE: assumes the TD waveform ends around t=0 */
    /* NOTE: estimate based on the 22 mode - see comments above for hphcFD when higher modes are included */
    double tc = EstimateInitialTime(listhlm, fLow);

    //double deltaf = 0.5 * 1./(2*(-tc)); /* Extra factor of 1/2 corresponding to 0-padding in TD by factor of 2 */
    double deltaf = 0.8 * 1./(2*(-tc)); /* Extra factor of 0.8 security margin for possible inaccuray of tf */

    /* Compute h22 FD on frequency samples to prepare FFT */
    /* NOTE: for fstartobs, see comments above for hphcFD when higher modes are included */
    ReImFrequencySeries* h22tilde = NULL;
    GenerateFDReImFrequencySeriesSingleMode(&h22tilde, listhlm, params->minf, params->maxf, fstartobs, deltaf, 0, 2, 2);

    /* Determine frequency windows */
    if(params->f1windowbeg==0.) f1windowbeg = fLow;
    else f1windowbeg = params->f1windowbeg;
    if(params->f2windowbeg==0.) f2windowbeg = 1.1*f1windowbeg;
    else f2windowbeg = params->f2windowbeg;
    if(params->f2windowend==0.) f2windowend = fHigh;
    else f2windowend = params->f2windowend;
    if(params->f1windowend==0.) f1windowend = 0.995*f2windowend;
    else f1windowend = params->f1windowend;

    /* Compute h22 TD by IFFT */
    ReImTimeSeries* h22TD = NULL;
    IFFTFrequencySeries(&h22TD, h22tilde, f1windowbeg, f2windowbeg, f1windowend, f2windowend, 1); /* Here hardcoded nzeropadding */

    /* Convert to Amp/Phase representation */
    AmpPhaseTimeSeries* h22TDAmpPhase = NULL;
    ReImTimeSeries_ToAmpPhase(&h22TDAmpPhase, h22TD);

    /* Output */
    Write_Wave_h22TD(params->outdir, params->outfile, h22TDAmpPhase, params->binaryout);
    exit(0);
  }
}
