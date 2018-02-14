#include "GenerateTDITD.h"

/************ Parsing arguments function ************/

/* Parse command line to initialize GenTDITDparams object */
/* Masses are input in solar masses and distances in Mpc - converted in SI for the internals */
static void parse_args_GenerateTDITD(ssize_t argc, char **argv, GenTDITDparams* params)
{
    char help[] = "\
GenerateWaveform by Sylvain Marsat, John Baker, and Philip Graff\n\
Copyright July 2015\n\
\n\
This program takes as input time series for hplus, hcross, generates and outputs a TDI signal in time domain.\n\
Arguments are as follows:\n\
\n\
--------------------------------------------------\n\
----- Physical Parameters ------------------------\n\
--------------------------------------------------\n\
 --phiRef              Phase (radians, default=0) -- used when input is h22\n\
 --inclination         Inclination (radians, default=PI/3) -- used when input is h22\n\
 --lambda              First angle for the position in the sky (radians, default=0)\n\
 --beta                Second angle for the position in the sky (radians, default=0)\n\
 --polarization        Polarization of source (radians, default=0)\n\
\n\
--------------------------------------------------\n\
----- Generation Parameters ----------------------\n\
--------------------------------------------------\n\
 --tagtdi              Tag choosing the set of TDI variables to use (default TDIAETXYZ)\n\
 --nsamplesinfile      Number of lines of inputs file\n\
 --binaryin            Tag for loading the data in gsl binary form instead of text (default false)\n\
 --binaryout           Tag for outputting the data in gsl binary form instead of text (default false)\n\
 --indir               Input directory\n\
 --infile              Input file name\n\
 --outdir              Output directory\n\
 --outfile             Output file name\n\
\n";

    ssize_t i;

    /* Set default values for the physical params */
    params->phiRef = 0;
    params->inclination = PI/3;
    params->lambda = 0;
    params->beta = 0;
    params->polarization = 0;

    /* Set default values for the generation params */
    params->tagtdi = TDIXYZ;
    params->nsamplesinfile = 0;    /* No default; has to be provided */
    strcpy(params->indir, "");   /* No default; has to be provided */
    strcpy(params->infile, "");  /* No default; has to be provided */
    strcpy(params->outdir, ".");
    strcpy(params->outfile, "generated_tdiTD.txt");

    /* Consume command line */
    for (i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--help") == 0) {
            fprintf(stdout,"%s", help);
            exit(0);
        } else if (strcmp(argv[i], "--phiRef") == 0) {
            params->phiRef = atof(argv[++i]);
        } else if (strcmp(argv[i], "--inclination") == 0) {
            params->inclination = atof(argv[++i]);
        } else if (strcmp(argv[i], "--lambda") == 0) {
            params->lambda = atof(argv[++i]);
        } else if (strcmp(argv[i], "--beta") == 0) {
            params->beta = atof(argv[++i]);
        } else if (strcmp(argv[i], "--polarization") == 0) {
            params->polarization = atof(argv[++i]);
        } else if (strcmp(argv[i], "--tagtdi") == 0) {
	          params->tagtdi = ParseTDItag(argv[++i]);
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

/* Read waveform time series in Re/Im form for hpTD and hcTD a single file */
/* NOTE: assumes the same number of points is used to represent each mode */
static void Read_Wave_hphcTD( RealTimeSeries** hptd, RealTimeSeries** hctd, const char dir[], const char file[], const int nsamples, int binary)
{
  /* Initalize and read input */
  gsl_matrix* inmatrix =  gsl_matrix_alloc(nsamples, 3);
  if (!binary) Read_Text_Matrix(dir, file, inmatrix);
  else Read_Matrix(dir, file, inmatrix);

  /* Initialize structures */
  RealTimeSeries_Init(hptd, nsamples);
  RealTimeSeries_Init(hctd, nsamples);

  /* Set values */
  gsl_vector_view timesview = gsl_matrix_column(inmatrix, 0);
  gsl_vector_view hptdview = gsl_matrix_column(inmatrix, 1);
  gsl_vector_view hctdview = gsl_matrix_column(inmatrix, 2);
  gsl_vector_memcpy((*hptd)->times, &timesview.vector);
  gsl_vector_memcpy((*hctd)->times, &timesview.vector);
  gsl_vector_memcpy((*hptd)->h, &hptdview.vector);
  gsl_vector_memcpy((*hctd)->h, &hctdview.vector);

  /* Clean up */
  gsl_matrix_free(inmatrix);
}
/* Output TDI 3 channels TD - real time series */
/* NOTE: assumes 3 channels with same times */
static void Write_TDITD(const char dir[], const char file[], RealTimeSeries* TDI1, RealTimeSeries* TDI2, RealTimeSeries* TDI3, int binary)
{
  /* Initialize output */
  /* NOTE: assumes identical times for all 3 TDI observables */
  int nbtimes = TDI1->times->size;
  gsl_matrix* outmatrix = gsl_matrix_alloc(nbtimes, 4);

  /* Set data */
  gsl_matrix_set_col(outmatrix, 0, TDI1->times);
  gsl_matrix_set_col(outmatrix, 1, TDI1->h);
  gsl_matrix_set_col(outmatrix, 2, TDI2->h);
  gsl_matrix_set_col(outmatrix, 3, TDI3->h);

  /* Output */
  if(!binary) Write_Text_Matrix(dir, file, outmatrix);
  else Write_Matrix(dir, file, outmatrix);
}


/***************** Main program *****************/

int main(int argc, char *argv[])
{
  int ret;

  LISAconstellation *variant = &LISAProposal;
  
  /* Initialize structure for parameters */
  GenTDITDparams* params;
  params = (GenTDITDparams*) malloc(sizeof(GenTDITDparams));
  memset(params, 0, sizeof(GenTDITDparams));

  /* Parse commandline to read parameters */
  parse_args_GenerateTDITD(argc, argv, params);

  /* Set geometric coefficients */
  SetCoeffsG(params->lambda, params->beta, params->polarization);

  /* Set aside the cases of the orbital delay and constellation response - written with TD amp/phase */
  /* If not using those tags, use the old processing that uses TD hplus, hcross interpolated */
  if(!(params->tagtdi==delayO) && !(params->tagtdi==y12L)) {

    /* Load TD hp, hc from file */
    RealTimeSeries* hptd = NULL;
    RealTimeSeries* hctd = NULL;
    Read_Wave_hphcTD(&hptd, &hctd, params->indir, params->infile, params->nsamplesinfile, params->binaryin);

    /* Interpolate hp, hc with gsl spline */
    /* NOTE: assumes indentical times vector */
    gsl_vector* times = hptd->times;
    int nbtimes = times->size;
    gsl_interp_accel* accel_hp = gsl_interp_accel_alloc();
    gsl_interp_accel* accel_hc = gsl_interp_accel_alloc();
    gsl_spline* spline_hp = gsl_spline_alloc(gsl_interp_cspline, nbtimes);
    gsl_spline* spline_hc = gsl_spline_alloc(gsl_interp_cspline, nbtimes);
    gsl_spline_init(spline_hp, gsl_vector_const_ptr(hptd->times,0), gsl_vector_const_ptr(hptd->h,0), nbtimes);
    gsl_spline_init(spline_hc, gsl_vector_const_ptr(hctd->times,0), gsl_vector_const_ptr(hctd->h,0), nbtimes);

    /* Here, hardcoded time margin that we will set to 0 on both sides, to avoid problems with delays extending past the input values */
    /* We take R/c as the maximum delay that can occur in the response, and adjust the margin accordingly */
    double maxdelay = variant->OrbitR/C_SI;
    double deltat = gsl_vector_get(times, 1) - gsl_vector_get(times, 0);
    int nptmargin = 2 * (int)(maxdelay/deltat);

    /* Evaluate TDI */
    RealTimeSeries* TDI1 = NULL;
    RealTimeSeries* TDI2 = NULL;
    RealTimeSeries* TDI3 = NULL;
    GenerateTDITD3Chanhphc(variant, &TDI1, &TDI2, &TDI3, spline_hp, spline_hc, accel_hp, accel_hc, times, nptmargin, params->tagtdi);

    /* Output */
    Write_TDITD(params->outdir, params->outfile, TDI1, TDI2, TDI3, params->binaryout);
  }

  /* Set aside the cases of the orbital delay and constellation response - written with TD amp/phase */
  /* NOTE: for now supports only single-mode (22) waveforms in amp/phase form */
  else if((params->tagtdi==delayO) || (params->tagtdi==y12L)) {
    /* Load 22-mode TD - can be the orbital-delayed 22 mode in the case of y12L */
    AmpPhaseTimeSeries* h22td = NULL;
    Read_AmpPhaseTimeSeries(&h22td, params->indir, params->infile, params->nsamplesinfile, params->binaryin);

    /* Interpolate amp, phase with gsl spline */
    gsl_vector* times = h22td->times;
    int nbtimes = times->size;
    gsl_interp_accel* accel_amp = gsl_interp_accel_alloc();
    gsl_interp_accel* accel_phase = gsl_interp_accel_alloc();
    gsl_spline* spline_amp = gsl_spline_alloc(gsl_interp_cspline, nbtimes);
    gsl_spline* spline_phase = gsl_spline_alloc(gsl_interp_cspline, nbtimes);
    gsl_spline_init(spline_amp, gsl_vector_const_ptr(h22td->times,0), gsl_vector_const_ptr(h22td->h_amp,0), nbtimes);
    gsl_spline_init(spline_phase, gsl_vector_const_ptr(h22td->times,0), gsl_vector_const_ptr(h22td->h_phase,0), nbtimes);

    /* Here, hardcoded time margin that we will set to 0 on both sides, to avoid problems with delays extending past the input values */
    /* We take R/c as the maximum delay that can occur in the response, and adjust the margin accordingly */
    double maxdelay = variant->OrbitR/C_SI;
    double deltat = gsl_vector_get(times, 1) - gsl_vector_get(times, 0);
    int nptmargin = 2 * (int)(maxdelay/deltat);

    if(params->tagtdi==delayO) {
      /* Evaluate orbital-delayed signal */
      AmpPhaseTimeSeries* h22tdO = NULL;
      Generateh22TDO(variant, &h22tdO, spline_amp, spline_phase, accel_amp, accel_phase, times, nptmargin);

      /* Output */
      Write_AmpPhaseTimeSeries(params->outdir, params->outfile, h22tdO, params->binaryout);
    }

    else if(params->tagtdi==y12L) {
      /* Evaluate y12 signal - constellation response only, input assumed to be h22tdO already */
      RealTimeSeries* y12td = NULL;
      Generatey12LTD(variant, &y12td, spline_amp, spline_phase, accel_amp, accel_phase, times, params->inclination, params->phiRef, nptmargin);

      /* Output */
      Write_RealTimeSeries(params->outdir, params->outfile, y12td, params->binaryout);
    }
  }
}
