#include "GenerateLLVTD.h"

/************ Parsing arguments function ************/

/* Parse command line to initialize GenLLVTDparams object */
/* Masses are input in solar masses and distances in Mpc - converted in SI for the internals */
static void parse_args_GenerateLLVTD(ssize_t argc, char **argv, GenLLVTDparams* params)
{
    char help[] = "\
GenerateWaveform by Sylvain Marsat, John Baker, and Philip Graff\n\
Copyright July 2015\n\
\n\
This program takes as input time series for hplus, hcross, generates and outputs a LLV signal in time domain.\n\
Arguments are as follows:\n\
\n\
--------------------------------------------------\n\
----- Physical Parameters ------------------------\n\
--------------------------------------------------\n\
 --lambda              First angle for the position in the sky (radians, default=0)\n\
 --beta                Second angle for the position in the sky (radians, default=0)\n\
 --polarization        Polarization of source (radians, default=0)\n\
\n\
--------------------------------------------------\n\
----- Generation Parameters ----------------------\n\
--------------------------------------------------\n\
 --tagnetwork          Tag choosing the detector network to use (default LLV)\n\
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
    params->lambda = 0;
    params->beta = 0;
    params->polarization = 0;

    /* Set default values for the generation params */
    params->tagnetwork = LHV;
    params->nsamplesinfile = 0;    /* No default; has to be provided */
    strcpy(params->indir, "");   /* No default; has to be provided */
    strcpy(params->infile, "");  /* No default; has to be provided */
    strcpy(params->outdir, ".");
    strcpy(params->outfile, "generated_LLVTD.txt");

    /* Consume command line */
    for (i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--help") == 0) {
            fprintf(stdout,"%s", help);
            exit(0);
        } else if (strcmp(argv[i], "--lambda") == 0) {
            params->lambda = atof(argv[++i]);
        } else if (strcmp(argv[i], "--beta") == 0) {
            params->beta = atof(argv[++i]);
        } else if (strcmp(argv[i], "--polarization") == 0) {
            params->polarization = atof(argv[++i]);
        } else if (strcmp(argv[i], "--tagnetwork") == 0) {
	  params->tagtdi = ParseNetworktag(argv[++i]); //
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
/* Output waveform in downsampled form, FD Amp/Pase, all hlm modes in a single file */
/* NOTE: assumes 3 channels with same times */
static void Write_LLVTD(const char dir[], const char file[], RealTimeSeries* LLV1, RealTimeSeries* LLV2, RealTimeSeries* LLV3, int binary)
{
  /* Initialize output */
  /* NOTE: assumes identical times for all 3 LLV observables */
  int nbtimes = LLV1->times->size;
  gsl_matrix* outmatrix = gsl_matrix_alloc(nbtimes, 4);

  /* Set data */
  gsl_matrix_set_col(outmatrix, 0, LLV1->times);
  gsl_matrix_set_col(outmatrix, 1, LLV1->h);
  gsl_matrix_set_col(outmatrix, 2, LLV2->h);
  gsl_matrix_set_col(outmatrix, 3, LLV3->h);

  /* Output */
  if(!binary) Write_Text_Matrix(dir, file, outmatrix);
  else Write_Matrix(dir, file, outmatrix);
}

/***************** Main program *****************/

int main(int argc, char *argv[])
{
  int ret;

  /* Initialize structure for parameters */
  GenLLVTDparams* params;
  params = (GenLLVTDparams*) malloc(sizeof(GenLLVTDparams));
  memset(params, 0, sizeof(GenLLVTDparams));

  /* Parse commandline to read parameters */
  parse_args_GenerateLLVTD(argc, argv, params);

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

  /* Evaluate LLV */
  RealTimeSeries* LLV1 = NULL;
  RealTimeSeries* LLV2 = NULL;
  RealTimeSeries* LLV3 = NULL;
  int nptmargin = 200; /* Here, hardcoded margin that we will set to 0 on both sides, to avoid problems with delays extending past the input values */
  GenerateTDITD3Chan(&TDI1, &TDI2, &TDI3, spline_hp, spline_hc, accel_hp, accel_hc, times, nptmargin, params->tagtdi);//



  /* Output */
  Write_TDITD(params->outdir, params->outfile, TDI1, TDI2, TDI3, params->binaryout);//
}
