#include "GenerateWaveform.h"

/************ Parsing arguments function ************/

/* Function to convert string input TDI string to TDItag */
static GenWavetag ParseGenWavetag(char* string) {
  GenWavetag tag;
  if(strcmp(string, "hlm")==0) tag = hlm;
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
 --nbmode              Number of modes of radiation to generate (1-5, default=5)\n\
 --fLow                Minimal frequency (Hz, default=0) - when too low, use first frequency covered by the ROM\n\
 --taggenwave          Tag choosing the wf format: hlm (default: downsampled modes in Amp/Phase form), hphcFD (hlm interpolated and summed), hphcTD (IFFT of hphcFD),\n\
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
    params->nbmode = 5;
    params->fLow = 0.;
    params->taggenwave = hlm;
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
        } else if (strcmp(argv[i], "--fLow") == 0) {
            params->fLow = atof(argv[++i]);
        } else if (strcmp(argv[i], "--taggenwave") == 0) {
	  params->taggenwave = ParseGenWavetag(argv[++i]);
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

/************ Functions to write waveforms to file ************/

/* Output waveform in downsampled form, FD Amp/Pase, all hlm modes in a single file */
/* NOTE: assumes the same number of points is used to represent each mode */
static void Write_Text_Wave_hlm(const char dir[], const char file[], ListmodesCAmpPhaseFrequencySeries* listhlm, int nbmodes)
{
  /* Initialize output */
  /* get length from 22 mode - NOTE: assumes the same for all modes */
  int nbfreq = ListmodesCAmpPhaseFrequencySeries_GetMode(listhlm, 2, 2)->freqseries->freq->size;
  gsl_matrix* outmatrix = gsl_matrix_alloc(nbfreq, 3*nbmodes); 

  /* Get data in the list of modes */
  CAmpPhaseFrequencySeries* mode;
  for(int i=0; i<nbmodes; i++) {
    mode = ListmodesCAmpPhaseFrequencySeries_GetMode(listhlm, listmode[i][0], listmode[i][1])->freqseries;
    gsl_matrix_set_col(outmatrix, 0+3*i, mode->freq);
    gsl_matrix_set_col(outmatrix, 1+3*i, mode->amp_real); /* amp_imag is 0 at this stage, we ignore it */
    gsl_matrix_set_col(outmatrix, 2+3*i, mode->phase);
  }

  /* Output */
  Write_Text_Matrix(dir, file, outmatrix);
}
/* Output waveform in frequency series form,Re/Im for hplus and hcross */
static void Write_Text_Wave_hphcFD(const char dir[], const char file[], ReImFrequencySeries* hptilde, ReImFrequencySeries* hctilde)
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
  Write_Text_Matrix(dir, file, outmatrix);
}
/* Output waveform in freqeucny sreies form,Re/Im for hplus and hcross */
static void Write_Text_Wave_hphcTD(const char dir[], const char file[], RealTimeSeries* hp, RealTimeSeries* hc)
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
  Write_Text_Matrix(dir, file, outmatrix);
}

/***************** Main program *****************/

int main(int argc, char *argv[])
{
  /* Initialize structure for parameters */
  GenWaveParams* params;
  params = (GenWaveParams*) malloc(sizeof(GenWaveParams));
  memset(params, 0, sizeof(GenWaveParams));
  
  /* Parse commandline to read parameters */
  parse_args_GenerateWaveform(argc, argv, params);

  /* Generate Fourier-domain waveform as a list of hlm modes */
  ListmodesCAmpPhaseFrequencySeries* listROM = NULL;
  SimEOBNRv2HMROM(&listROM, params->nbmode, params->tRef, params->phiRef, params->fRef, params->m1*MSUN_SI, params->m2*MSUN_SI, params->distance*1e6*PC_SI);

  if(params->taggenwave==hlm) {
    Write_Text_Wave_hlm(params->outdir, params->outfile, listROM, params->nbmode);
    exit(0);
  }

  else {

    /* Determine deltaf so that N deltat = 1/deltaf > 2*tc where tc is the time to coalescence estimated from Psi22 */
    /* Assumes the TD waveform ends around t=0 */
    double tc = EstimateInitialTime(listROM, params->fLow);
    double deltaf = 0.5 * 1./(2*(-tc)); /* Extra factor of 1/2 corresponding to 0-padding in TD by factor of 2 */

    /* Compute hplus, hcross FD */
    ReImFrequencySeries* hptilde = NULL;
    ReImFrequencySeries* hctilde = NULL;
    GeneratehphcFDReImFrequencySeries(&hptilde, &hctilde, listROM, params->fLow, deltaf, 0, params->nbmode, params->inclination, 0., 1);

    /* Output */
    if(params->taggenwave==hphcFD) {

      Write_Text_Wave_hphcFD(params->outdir, params->outfile, hptilde, hctilde);
      exit(0);
    }

    else {

      /* Determine frequency windows - min,max frequencies determined directly from the Amp/Phase hlm's */
      double fLowROM = ListmodesCAmpPhaseFrequencySeries_minf(listROM);
      double fHighROM = ListmodesCAmpPhaseFrequencySeries_maxf(listROM);
      double f1windowbeg = fmax(params->fLow, fLowROM);
      double f2windowbeg = 1.1 * f1windowbeg; /* Here hardcoded relative width of the window */
      double f2windowend = fHighROM;
      double f1windowend = 0.995 * f2windowend; /* Here hardcoded relative width of the window */
    
      /* Compute hplus, hcross TD */
      RealTimeSeries* hp = NULL;
      RealTimeSeries* hc = NULL;
      IFFTFrequencySeries(&hp, hptilde, f1windowbeg, f2windowbeg, f1windowend, f2windowend, 3); /* Here hardcoded nzeropadding */
      IFFTFrequencySeries(&hc, hctilde, f1windowbeg, f2windowbeg, f1windowend, f2windowend, 3); /* Here hardcoded nzeropadding */

      /* Output */
      if(params->taggenwave==hphcTD) {
	Write_Text_Wave_hphcTD(params->outdir, params->outfile, hp, hc);
	exit(0);
      }
    }   
  }
}
