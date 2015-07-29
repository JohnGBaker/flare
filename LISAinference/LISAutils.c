#include "LISAutils.h"

/************ Global Parameters ************/

LISAParams* injectedparams = NULL;
LISAGlobalParams* globalparams = NULL;
LISAParams* templateparams = NULL;
LISAPrior* priorParams = NULL;

/************ Functions to initalize and clean up structure for the signals ************/

void LISASignalCAmpPhase_Cleanup(LISASignalCAmpPhase* signal) {
  if(signal->TDIASignal) ListmodesCAmpPhaseFrequencySeries_Destroy(signal->TDIASignal);
  if(signal->TDIESignal) ListmodesCAmpPhaseFrequencySeries_Destroy(signal->TDIESignal);
  if(signal->TDITSignal) ListmodesCAmpPhaseFrequencySeries_Destroy(signal->TDITSignal);
  free(signal);
}

void LISASignalCAmpPhase_Init(LISASignalCAmpPhase** signal) {
  if(!signal) exit(1);
  /* Create storage for structures */
  if(!*signal) *signal = malloc(sizeof(LISASignalCAmpPhase));
  else
  {
    LISASignalCAmpPhase_Cleanup(*signal);
  }
  (*signal)->TDIASignal = NULL;
  (*signal)->TDIESignal = NULL;
  (*signal)->TDITSignal = NULL;
}

void LISASignalReIm_Cleanup(LISASignalReIm* signal) {
  if(signal->TDIASignal) ReImFrequencySeries_Cleanup(signal->TDIASignal);
  if(signal->TDIESignal) ReImFrequencySeries_Cleanup(signal->TDIESignal);
  if(signal->TDITSignal) ReImFrequencySeries_Cleanup(signal->TDITSignal);
  free(signal);
}

void LISASignalReIm_Init(LISASignalReIm** signal) {
  if(!signal) exit(1);
  /* Create storage for structures */
  if(!*signal) *signal = malloc(sizeof(LISASignalReIm));
  else
  {
    LISASignalReIm_Cleanup(*signal);
  }
  (*signal)->TDIASignal = NULL;
  (*signal)->TDIESignal = NULL;
  (*signal)->TDITSignal = NULL;
}

void LISAInjectionReIm_Cleanup(LISAInjectionReIm* signal) {
  if(signal->TDIASignal) ReImFrequencySeries_Cleanup(signal->TDIASignal);
  if(signal->TDIESignal) ReImFrequencySeries_Cleanup(signal->TDIESignal);
  if(signal->TDITSignal) ReImFrequencySeries_Cleanup(signal->TDITSignal);
  if(signal->freq) gsl_vector_free(signal->freq);
  if(signal->noisevaluesA) gsl_vector_free(signal->noisevaluesA);
  if(signal->noisevaluesE) gsl_vector_free(signal->noisevaluesE);
  if(signal->noisevaluesT) gsl_vector_free(signal->noisevaluesT);
  free(signal);
}

void LISAInjectionReIm_Init(LISAInjectionReIm** signal) {
  if(!signal) exit(1);
  /* Create storage for structures */
  if(!*signal) *signal = malloc(sizeof(LISAInjectionReIm));
  else
  {
    LISAInjectionReIm_Cleanup(*signal);
  }
  (*signal)->TDIASignal = NULL;
  (*signal)->TDIESignal = NULL;
  (*signal)->TDITSignal = NULL;
  (*signal)->freq = NULL;
  (*signal)->noisevaluesA = NULL;
  (*signal)->noisevaluesE = NULL;
  (*signal)->noisevaluesT = NULL;
}


/************ Parsing arguments function ************/

/* Parse command line to initialize LISAParams, LISAPrior, and LISARunParams objects */
void parse_args_LISA(ssize_t argc, char **argv, 
  LISAParams* params,
  LISAGlobalParams* globalparams, 
  LISAPrior* prior, 
  LISARunParams* run) 
{
    char help[] = "\
LISAInference by Sylvain Marsat, John Baker, and Philip Graff\n\
Copyright July 2015\n\
\n\
This program performs rapid parameter estimation for LIGO and LISA CBC sources in the no-noise case.\n\
Arguments are as follows:\n\
\n\
--------------------------------------------------\n\
----- Injected Signal Parameters -----------------\n\
--------------------------------------------------\n\
 --tRef                Time at reference frequency (sec, default=0)\n\
 --phiRef              Orbital phase at reference frequency (radians, default=0)\n\
 --m1                  Component mass 1 in Solar masses (larger, default=2e6)\n\
 --m2                  Component mass 2 in Solar masses (smaller, default=1e6)\n\
 --distance            Distance to source in Mpc (default=1e9)\n\
 --lambda              First angle for the position in the sky (radians, default=0)\n\
 --beta                Second angle for the position in the sky (radians, default=0)\n\
 --inclination         Inclination of source orbital plane to observer line of sight\n\
                       (radians, default=PI/3)\n\
 --polarization        Polarization of source (radians, default=0)\n\
 --nbmode              Number of modes of radiation to generate (1-5, default=5)\n\
 --snr                 Use a target network SNR for the injection by rescaling distance\n\
\n\
-----------------------------------------------------------------\n\
----- Global Waveform/Inner products Parameters -----------------\n\
-----------------------------------------------------------------\n\
 --fRef                Reference frequency (Hz, default=0, interpreted as Mf=0.14)\n\
 --deltatobs           Observation time (years, default=2)\n\
 --fmin                Minimal frequency (Hz, default=0) - when set to 0, use the first frequency covered by the noise data of the detector\n\
 --nbmodeinj           Number of modes of radiation to use for the injection (1-5, default=5)\n\
 --nbmodetemp          Number of modes of radiation to use for the templates (1-5, default=5)\n\
 --tagint              Tag choosing the integrator: 0 for wip (default), 1 for linear integration\n\
 --nbptsoverlap        Number of points to use for linear integration (default 32768)\n\
\n\
--------------------------------------------------\n\
----- Prior Boundary Settings --------------------\n\
--------------------------------------------------\n\
 --deltaT              Half-width of time prior (sec, default=1e5)\n\
 --comp-min            Minimum component mass in Solar masses (default=1e4)\n\
 --comp-max            Maximum component mass in Solar masses (default=1e8)\n\
 --mtot-min            Minimum total mass in Solar masses (default=5e4)\n\
 --mtot-max            Maximum total mass in Solar masses (default=1e8)\n\
 --q-max               Maximum mass ratio, m1/m2 (default=11.98, minimum is 1)\n\
 --dist-min            Minimum distance to source (Mpc, default=100)\n\
 --dist-max            Maximum distance to source (Mpc, default=40*1e3)\n\
\n\
--------------------------------------------------\n\
----- Fix Parameters In Sampling -----------------\n\
--------------------------------------------------\n\
 --fix-m1              Fix mass 1\n\
 --fix-m2              Fix mass 2\n\
 --fix-dist            Fix distance\n\
 --fix-time            Fix reference time\n\
 --fix-phase           Fix reference phase\n\
 --fix-lambda          Fix lambda\n\
 --fix-beta            Fix beta\n\
 --fix-inc             Fix inclination\n\
 --fix-pol             Fix polarization\n\
\n\
--------------------------------------------------\n\
----- BAMBI Sampler Settings ---------------------\n\
--------------------------------------------------\n\
 --eff                 Target efficiency of sampling (default=0.1)\n\
 --tol                 Tolerance for evidence calculation convergence (default=0.5)\n\
 --nlive               Number of live points for sampling (default=1000)\n\
 --bambi               Use BAMBI's neural network logL learning (no option, default off)\n\
 --resume              Resume from a previous run (no option, default off)\n\
 --outroot             Root for output files (default='chains/LISAinference_')\n\
 --netfile             Neural network settings file if using --bambi (default='LISAinference.inp')\n\
\n";

    ssize_t i;

    /* set default values for the injection params */
    params->tRef = 0.;
    params->phiRef = 0.;
    params->m1 = 2*1e6;
    params->m2 = 1*1e6;
    params->distance = 1e3;
    params->lambda = 0.;
    params->beta = 0.;
    params->inclination = PI/3.;
    params->polarization = 0.;
    params->nbmode = 5;

    /* set default values for the global params */
    globalparams->fRef = 0.;
    globalparams->deltatobs = 2.;
    globalparams->fmin = 0.;
    globalparams->nbmodeinj = 5;
    globalparams->nbmodetemp = 5;
    globalparams->tagint = 0;
    globalparams->nbptsoverlap = 32768;

    /* set default values for the prior limits */
    prior->deltaT = 1.e5;
    prior->comp_min = 1e4;
    prior->comp_max = 1e8;
    prior->mtot_min = 5e4;
    prior->mtot_max = 1e8;
    prior->qmax = 11.98;
    prior->dist_min = 100.0;
    prior->dist_max = 40*1e3;
    prior->fix_m1 = NAN;
    prior->fix_m2 = NAN;
    prior->fix_dist = NAN;
    prior->fix_time = NAN;
    prior->fix_phase = NAN;
    prior->fix_pol = NAN;
    prior->fix_lambda = NAN;
    prior->fix_beta = NAN;
    prior->fix_inc = NAN;
    prior->snr_target = NAN;

    /* set default values for the run settings */
    run->eff = 0.1;
    run->tol = 0.5;
    run->nlive = 1000;
    strcpy(run->outroot, "chains/LISAinference_");
    run->bambi = 0;
    run->resume = 0;
    strcpy(run->netfile, "LISAinference.inp");

    /* Consume command line */
    for (i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--help") == 0) {
            fprintf(stdout,"%s", help);
            exit(0);
        } else if (strcmp(argv[i], "--tRef") == 0) {
            params->tRef = atof(argv[++i]);
        } else if (strcmp(argv[i], "--phiRef") == 0) {
            params->phiRef = atof(argv[++i]);
        } else if (strcmp(argv[i], "--m1") == 0) {
            params->m1 = atof(argv[++i]);
        } else if (strcmp(argv[i], "--m2") == 0) {
            params->m2 = atof(argv[++i]);
        } else if (strcmp(argv[i], "--distance") == 0) {
            params->distance = atof(argv[++i]);
        } else if (strcmp(argv[i], "--lambda") == 0) {
            params->lambda = atof(argv[++i]);
        } else if (strcmp(argv[i], "--beta") == 0) {
            params->beta = atof(argv[++i]);
        } else if (strcmp(argv[i], "--inclination") == 0) {
            params->inclination = atof(argv[++i]);
        } else if (strcmp(argv[i], "--polarization") == 0) {
            params->polarization = atof(argv[++i]);
        } else if (strcmp(argv[i], "--fRef") == 0) {
            globalparams->fRef = atof(argv[++i]);
        } else if (strcmp(argv[i], "--deltatobs") == 0) {
            globalparams->deltatobs = atof(argv[++i]);
        } else if (strcmp(argv[i], "--fmin") == 0) {
            globalparams->fmin = atof(argv[++i]);
        } else if (strcmp(argv[i], "--nbmodeinj") == 0) {
            globalparams->nbmodeinj = atof(argv[++i]);
        } else if (strcmp(argv[i], "--nbmodetemp") == 0) {
            globalparams->nbmodetemp = atof(argv[++i]);
        } else if (strcmp(argv[i], "--tagint") == 0) {
            globalparams->tagint = atof(argv[++i]);
        } else if (strcmp(argv[i], "--nbptsoverlap") == 0) {
            globalparams->nbptsoverlap = atof(argv[++i]);
        } else if (strcmp(argv[i], "--deltaT") == 0) {
            prior->deltaT = atof(argv[++i]);
        } else if (strcmp(argv[i], "--comp-min") == 0) {
            prior->comp_min = atof(argv[++i]);
        } else if (strcmp(argv[i], "--comp-max") == 0) {
            prior->comp_max = atof(argv[++i]);
        } else if (strcmp(argv[i], "--mtot-min") == 0) {
            prior->mtot_min = atof(argv[++i]);
        } else if (strcmp(argv[i], "--mtot-max") == 0) {
            prior->mtot_max = atof(argv[++i]);
        } else if (strcmp(argv[i], "--q-max") == 0) {
            prior->qmax = atof(argv[++i]);
        } else if (strcmp(argv[i], "--dist-min") == 0) {
            prior->dist_min = atof(argv[++i]);
        } else if (strcmp(argv[i], "--dist-max") == 0) {
            prior->dist_max = atof(argv[++i]);
        } else if (strcmp(argv[i], "--fix-m1") == 0) {
            prior->fix_m1 = atof(argv[++i]);
        } else if (strcmp(argv[i], "--fix-m2") == 0) {
            prior->fix_m2 = atof(argv[++i]);
        } else if (strcmp(argv[i], "--fix-dist") == 0) {
            prior->fix_dist = atof(argv[++i]);
        } else if (strcmp(argv[i], "--fix-lambda") == 0) {
            prior->fix_lambda = atof(argv[++i]);
        } else if (strcmp(argv[i], "--fix-beta") == 0) {
            prior->fix_beta = atof(argv[++i]);
        } else if (strcmp(argv[i], "--fix-time") == 0) {
            prior->fix_time = atof(argv[++i]);
        } else if (strcmp(argv[i], "--fix-phase") == 0) {
            prior->fix_phase = atof(argv[++i]);
        } else if (strcmp(argv[i], "--fix-inc") == 0) {
            prior->fix_inc = atof(argv[++i]);
        } else if (strcmp(argv[i], "--fix-pol") == 0) {
            prior->fix_pol = atof(argv[++i]);
        } else if (strcmp(argv[i], "--snr") == 0) {
            prior->snr_target = atof(argv[++i]);
        } else if (strcmp(argv[i], "--eff") == 0) {
            run->eff = atof(argv[++i]);
        } else if (strcmp(argv[i], "--tol") == 0) {
            run->tol = atof(argv[++i]);
        } else if (strcmp(argv[i], "--nlive") == 0) {
            run->nlive = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--bambi") == 0) {
            run->bambi = 1;
        } else if (strcmp(argv[i], "--resume") == 0) {
            run->resume = 1;
        } else if (strcmp(argv[i], "--outroot") == 0) {
            strcpy(run->outroot, argv[++i]);
        } else if (strcmp(argv[i], "--netfile") == 0) {
            strcpy(run->netfile, argv[++i]);
        } else {
            printf("Error: invalid option: %s\n", argv[i]);
            goto fail;
        }
    }

    return;

    fail:
    exit(1);
}

/***************************** Functions to generate signals and compute likelihoods ******************************/

/* Function generating a LISA signal as a list of modes in CAmp/Phase form, from LISA parameters */
int LISAGenerateSignalCAmpPhase(
  struct tagLISAParams* params,            /* Input: set of LISA parameters of the signal */
  struct tagLISASignalCAmpPhase* signal)   /* Output: structure for the generated signal */
{
  int ret;
  ListmodesCAmpPhaseFrequencySeries* listROM = NULL;
  ListmodesCAmpPhaseFrequencySeries* listTDIA = NULL;
  ListmodesCAmpPhaseFrequencySeries* listTDIE = NULL;
  ListmodesCAmpPhaseFrequencySeries* listTDIT = NULL;

  /* Checking that the global injectedparams has been set up */
  if (!injectedparams) {
    printf("Error: when calling LISAGenerateSignal, injectedparams points to NULL.\n");
    exit(1);
  }
  /* Should add more error checking ? */
  /* Generate the waveform with the ROM */
  /* Note: SimEOBNRv2HMROM accepts masses and distances in SI units, whereas LISA params is in solar masses and Mpc */
  ret = SimEOBNRv2HMROM(&listROM, params->nbmode, params->tRef - injectedparams->tRef, params->phiRef, globalparams->fRef, (params->m1)*MSUN_SI, (params->m2)*MSUN_SI, (params->distance)*1e6*PC_SI);

  /* If the ROM waveform generation failed (e.g. parameters were out of bounds) return FAILURE */
  if(ret==FAILURE) return FAILURE;

  /* Process the waveform through the LISA response */
  //WARNING: tRef is ignored for now, i.e. set to 0
  //TESTING
  //clock_t tbeg, tend;
  //tbeg = clock();
  LISASimFDResponseTDIAET(&listROM, &listTDIA, &listTDIE, &listTDIT, params->tRef, params->lambda, params->beta, params->inclination, params->polarization);
  //tend = clock();
  //printf("time LISASimFDResponse: %g\n", (double) (tend-tbeg)/CLOCKS_PER_SEC);
  //

  /* Precompute the inner products (h|h) - takes into account the length of the observation with deltatobs */
  double Mfstartobs = NewtonianfoftGeom(params->m1 / params->m2, (globalparams->deltatobs * YRSID_SI) / ((params->m1 + params->m2) * MTSUN_SI));
  double fstartobs = Mfstartobs / ((params->m1 + params->m2) * MTSUN_SI);
  double fLow = fmax(__LISASimFD_Noise_fLow, globalparams->fmin);
  double fHigh = __LISASimFD_Noise_fHigh;
  //TESTING
  //tbeg = clock();
  double TDIAhh = FDListmodesOverlap(listTDIA, listTDIA, NoiseSnA, fLow, fHigh, fstartobs, fstartobs, globalparams->tagint);
  double TDIEhh = FDListmodesOverlap(listTDIE, listTDIE, NoiseSnE, fLow, fHigh, fstartobs, fstartobs, globalparams->tagint);
  double TDIThh = FDListmodesOverlap(listTDIT, listTDIT, NoiseSnT, fLow, fHigh, fstartobs, fstartobs, globalparams->tagint);
  //tend = clock();
  //printf("time SNRs: %g\n", (double) (tend-tbeg)/CLOCKS_PER_SEC);
  //

  /* Output and clean up */
  signal->TDIASignal = listTDIA;
  signal->TDIESignal = listTDIE;
  signal->TDITSignal = listTDIT;
  signal->TDIAhh = TDIAhh;
  signal->TDIEhh = TDIEhh;
  signal->TDIThh = TDIThh;

  ListmodesCAmpPhaseFrequencySeries_Destroy(listROM);
  return SUCCESS;
}

/* Function generating a LISA signal as a frequency series in Re/Im form where the modes have been summed, from LISA parameters - takes as argument the frequencies on which to evaluate */
int LISAGenerateSignalReIm(
  struct tagLISAParams* params,       /* Input: set of LISA parameters of the template */
  gsl_vector* freq,                   /* Input: frequencies on which evaluating the waveform (from the injection) */
  struct tagLISASignalReIm* signal)   /* Output: structure for the generated signal */
{
  int ret;
  ListmodesCAmpPhaseFrequencySeries* listROM = NULL;
  ListmodesCAmpPhaseFrequencySeries* listTDIA = NULL;
  ListmodesCAmpPhaseFrequencySeries* listTDIE = NULL;
  ListmodesCAmpPhaseFrequencySeries* listTDIT = NULL;

  /* Checking that the global injectedparams has been set up */
  if (!injectedparams) {
    printf("Error: when calling LISAGenerateSignal, injectedparams points to NULL.\n");
    exit(1);
  }
  /* Should add more error checking ? */
  /* Generate the waveform with the ROM */
  /* Note: SimEOBNRv2HMROM accepts masses and distances in SI units, whereas LISA params is in solar masses and Mpc */
  ret = SimEOBNRv2HMROM(&listROM, params->nbmode, params->tRef - injectedparams->tRef, params->phiRef, globalparams->fRef, (params->m1)*MSUN_SI, (params->m2)*MSUN_SI, (params->distance)*1e6*PC_SI);

  /* If the ROM waveform generation failed (e.g. parameters were out of bounds) return FAILURE */
  if(ret==FAILURE) return FAILURE;

  /* Process the waveform through the LISA response */
  //WARNING: tRef is ignored for now, i.e. set to 0
  //TESTING
  //clock_t tbeg, tend;
  //tbeg = clock();
  LISASimFDResponseTDIAET(&listROM, &listTDIA, &listTDIE, &listTDIT, params->tRef, params->lambda, params->beta, params->inclination, params->polarization);
  //tend = clock();
  //printf("time LISASimFDResponse: %g\n", (double) (tend-tbeg)/CLOCKS_PER_SEC);
  //

  /* Initialize structures for the ReIm frequency series */
  int nbpts = (int) freq->size;
  ReImFrequencySeries* TDIA = NULL;
  ReImFrequencySeries* TDIE = NULL;
  ReImFrequencySeries* TDIT = NULL;
  ReImFrequencySeries_Init(&TDIA, nbpts);
  ReImFrequencySeries_Init(&TDIE, nbpts);
  ReImFrequencySeries_Init(&TDIT, nbpts);

  /* Compute the Re/Im frequency series - takes into account the length of the observation with deltatobs */
  double Mfstartobs = NewtonianfoftGeom(params->m1 / params->m2, (globalparams->deltatobs * YRSID_SI) / ((params->m1 + params->m2) * MTSUN_SI));
  double fstartobs = Mfstartobs / ((params->m1 + params->m2) * MTSUN_SI);
  //TESTING
  //tbeg = clock();
  ReImFrequencySeries_SumListmodesCAmpPhaseFrequencySeries(TDIA, listTDIA, freq, fstartobs);
  ReImFrequencySeries_SumListmodesCAmpPhaseFrequencySeries(TDIE, listTDIE, freq, fstartobs);
  ReImFrequencySeries_SumListmodesCAmpPhaseFrequencySeries(TDIT, listTDIT, freq, fstartobs);
  //tend = clock();
  //printf("time ReIm: %g\n", (double) (tend-tbeg)/CLOCKS_PER_SEC);
  //

  /* Output and clean up */
  signal->TDIASignal = TDIA;
  signal->TDIESignal = TDIE;
  signal->TDITSignal = TDIT;

  ListmodesCAmpPhaseFrequencySeries_Destroy(listROM);
  ListmodesCAmpPhaseFrequencySeries_Destroy(listTDIA);
  ListmodesCAmpPhaseFrequencySeries_Destroy(listTDIE);
  ListmodesCAmpPhaseFrequencySeries_Destroy(listTDIT);
  return SUCCESS;
}

/* Function generating a LISA injection signal as a frequency series in Re/Im form where the modes have been summed, from LISA parameters - determines the frequencies */
int LISAGenerateInjectionReIm(
  struct tagLISAParams* injectedparams,      /* Input: set of LISA parameters of the template */
  double fLow,                               /* Input: additional lower frequency limit (argument fmin) */
  int nbpts,                                 /* Input: number of frequency samples */
  struct tagLISAInjectionReIm* injection)    /* Output: structure for the generated signal */
{
  int ret;
  ListmodesCAmpPhaseFrequencySeries* listROM = NULL;
  ListmodesCAmpPhaseFrequencySeries* listTDIA = NULL;
  ListmodesCAmpPhaseFrequencySeries* listTDIE = NULL;
  ListmodesCAmpPhaseFrequencySeries* listTDIT = NULL;

  /* Should add more error checking ? */
  /* Generate the waveform with the ROM */
  /* Note: SimEOBNRv2HMROM accepts masses and distances in SI units, whereas LISA params is in solar masses and Mpc */
  ret = SimEOBNRv2HMROM(&listROM, injectedparams->nbmode, 0., injectedparams->phiRef, globalparams->fRef, (injectedparams->m1)*MSUN_SI, (injectedparams->m2)*MSUN_SI, (injectedparams->distance)*1e6*PC_SI);

  /* If the ROM waveform generation failed (e.g. parameters were out of bounds) return FAILURE */
  if(ret==FAILURE) return FAILURE;

  /* Process the waveform through the LISA response */
  //WARNING: tRef is ignored for now, i.e. set to 0
  //TESTING
  //clock_t tbeg, tend;
  //tbeg = clock();
  LISASimFDResponseTDIAET(&listROM, &listTDIA, &listTDIE, &listTDIT, injectedparams->tRef, injectedparams->lambda, injectedparams->beta, injectedparams->inclination, injectedparams->polarization);
  //tend = clock();
  //printf("time LISASimFDResponse: %g\n", (double) (tend-tbeg)/CLOCKS_PER_SEC);
  //

  /* Determine the frequency vector - uses the fact that the detector limiting frequencies are the same in all channels - takes into account the length of the observation with deltatobs */
  gsl_vector* freq = gsl_vector_alloc(nbpts);
  double Mfstartobs = NewtonianfoftGeom(injectedparams->m1 / injectedparams->m2, (globalparams->deltatobs * YRSID_SI) / ((injectedparams->m1 + injectedparams->m2) * MTSUN_SI));
  double fstartobs = Mfstartobs / ((injectedparams->m1 + injectedparams->m2) * MTSUN_SI);
  double fLowCut = fmax(fmax(__LISASimFD_Noise_fLow, fLow), fstartobs);
  double fHigh = __LISASimFD_Noise_fHigh;
  ListmodesSetLogFrequencies(listROM, fLowCut, fHigh, nbpts, freq);

  /* Initialize structures for the ReIm frequency series */
  ReImFrequencySeries* TDIA = NULL;
  ReImFrequencySeries* TDIE = NULL;
  ReImFrequencySeries* TDIT = NULL;
  ReImFrequencySeries_Init(&TDIA, nbpts);
  ReImFrequencySeries_Init(&TDIE, nbpts);
  ReImFrequencySeries_Init(&TDIT, nbpts);

  /* Compute the Re/Im frequency series */
  //TESTING
  //tbeg = clock();
  ReImFrequencySeries_SumListmodesCAmpPhaseFrequencySeries(TDIA, listTDIA, freq, fstartobs);
  ReImFrequencySeries_SumListmodesCAmpPhaseFrequencySeries(TDIE, listTDIE, freq, fstartobs);
  ReImFrequencySeries_SumListmodesCAmpPhaseFrequencySeries(TDIT, listTDIT, freq, fstartobs);
  //tend = clock();
  //printf("time ReIm: %g\n", (double) (tend-tbeg)/CLOCKS_PER_SEC);
  //

  /* Compute the noise values */
  gsl_vector* noisevaluesA = gsl_vector_alloc(nbpts);
  gsl_vector* noisevaluesE = gsl_vector_alloc(nbpts);
  gsl_vector* noisevaluesT = gsl_vector_alloc(nbpts);
  EvaluateNoise(noisevaluesA, freq, NoiseSnA, __LISASimFD_Noise_fLow, __LISASimFD_Noise_fHigh);
  EvaluateNoise(noisevaluesE, freq, NoiseSnE, __LISASimFD_Noise_fLow, __LISASimFD_Noise_fHigh);
  EvaluateNoise(noisevaluesT, freq, NoiseSnT, __LISASimFD_Noise_fLow, __LISASimFD_Noise_fHigh);

  /* Output and clean up */
  injection->TDIASignal = TDIA;
  injection->TDIESignal = TDIE;
  injection->TDITSignal = TDIT;
  injection->freq = freq;
  injection->noisevaluesA = noisevaluesA;
  injection->noisevaluesE = noisevaluesE;
  injection->noisevaluesT = noisevaluesT;

  ListmodesCAmpPhaseFrequencySeries_Destroy(listROM);
  ListmodesCAmpPhaseFrequencySeries_Destroy(listTDIA);
  ListmodesCAmpPhaseFrequencySeries_Destroy(listTDIE);
  ListmodesCAmpPhaseFrequencySeries_Destroy(listTDIT);
  return SUCCESS;
}

/* Log-Likelihood function */

double CalculateLogLCAmpPhase(LISAParams *params, LISASignalCAmpPhase* injection)
{
  double logL = -DBL_MAX;
  int ret;

  /* Generating the signal in the three detectors for the input parameters */
  LISASignalCAmpPhase* generatedsignal = NULL;
  LISASignalCAmpPhase_Init(&generatedsignal);
  //TESTING
  //clock_t tbeg, tend;
  //tbeg = clock();
  ret = LISAGenerateSignalCAmpPhase(params, generatedsignal);
  //tend = clock();
  //printf("time GenerateSignal: %g\n", (double) (tend-tbeg)/CLOCKS_PER_SEC);
  //

  /* If LISAGenerateSignal failed (e.g. parameters out of bound), silently return -Infinity logL */
  if(ret==FAILURE) {
    logL = -DBL_MAX;
  }
  else if(ret==SUCCESS) {
    /* Computing the likelihood for each TDI channel - fstartobs is the max between the fstartobs of the injected and generated signals */
    double Mfstartobsinjected = NewtonianfoftGeom(injectedparams->m1 / injectedparams->m2, (globalparams->deltatobs * YRSID_SI) / ((injectedparams->m1 + injectedparams->m2) * MTSUN_SI));
    double Mfstartobsgenerated = NewtonianfoftGeom(params->m1 / params->m2, (globalparams->deltatobs * YRSID_SI) / ((params->m1 + params->m2) * MTSUN_SI));
    double fstartobsinjected = Mfstartobsinjected / ((injectedparams->m1 + injectedparams->m2) * MTSUN_SI);
    double fstartobsgenerated = Mfstartobsgenerated / ((params->m1 + params->m2) * MTSUN_SI);
    double fLow = fmax(__LISASimFD_Noise_fLow, globalparams->fmin);
    double fHigh = __LISASimFD_Noise_fHigh;
    //TESTING
    //tbeg = clock();
    double loglikelihoodTDIA = FDLogLikelihood(injection->TDIASignal, generatedsignal->TDIASignal, NoiseSnA, fLow, fHigh, injection->TDIAhh, generatedsignal->TDIAhh, fstartobsinjected, fstartobsgenerated, globalparams->tagint);
    double loglikelihoodTDIE = FDLogLikelihood(injection->TDIESignal, generatedsignal->TDIESignal, NoiseSnE, fLow, fHigh, injection->TDIEhh, generatedsignal->TDIEhh, fstartobsinjected, fstartobsgenerated, globalparams->tagint);
    double loglikelihoodTDIT = FDLogLikelihood(injection->TDITSignal, generatedsignal->TDITSignal, NoiseSnT, fLow, fHigh, injection->TDIThh, generatedsignal->TDIThh, fstartobsinjected, fstartobsgenerated, globalparams->tagint);
    //tend = clock();
    //printf("time Overlaps: %g\n", (double) (tend-tbeg)/CLOCKS_PER_SEC);
    //

    /* Output: value of the loglikelihood for the combined signals, assuming noise independence */
    logL = loglikelihoodTDIA + loglikelihoodTDIE + loglikelihoodTDIT;
  }

  /* Clean up */
  LISASignalCAmpPhase_Cleanup(generatedsignal);

  return logL;
}

double CalculateLogLReIm(LISAParams *params, LISAInjectionReIm* injection)
{
  double logL = -DBL_MAX;
  int ret;

  /* Frequency vector - assumes common to A,E,T, i.e. identical fLow, fHigh in all channels */
  gsl_vector* freq = injection->freq;

  /* Generating the signal in the three detectors for the input parameters */
  LISASignalReIm* generatedsignal = NULL;
  LISASignalReIm_Init(&generatedsignal);
  //TESTING
  //clock_t tbeg, tend;
  //tbeg = clock();
  ret = LISAGenerateSignalReIm(params, freq, generatedsignal);
  //tend = clock();
  //printf("time GenerateSignal: %g\n", (double) (tend-tbeg)/CLOCKS_PER_SEC);
  //
  
  //
  printf("LISAGenerateSignalReIm return code: %d\n", ret);

  /* If LISAGenerateSignal failed (e.g. parameters out of bound), silently return -Infinity logL */
  if(ret==FAILURE) {
    //
    printf("Wave generation failed.\n");
    logL = -DBL_MAX;
  }
  else if(ret==SUCCESS) {
    /* Computing the likelihood for each TDI channel - fstartobs has already been taken into account */
    //TESTING
    //tbeg = clock();
    double loglikelihoodTDIA = FDLogLikelihoodReIm(injection->TDIASignal, generatedsignal->TDIASignal, injection->noisevaluesA);
    double loglikelihoodTDIE = FDLogLikelihoodReIm(injection->TDIESignal, generatedsignal->TDIESignal, injection->noisevaluesE);
    double loglikelihoodTDIT = FDLogLikelihoodReIm(injection->TDITSignal, generatedsignal->TDITSignal, injection->noisevaluesT);
    //tend = clock();
    //printf("time Overlaps: %g\n", (double) (tend-tbeg)/CLOCKS_PER_SEC);
    //

    //
    printReImFrequencySeries(generatedsignal->TDIASignal, 1000, 1020);

    /* Output: value of the loglikelihood for the combined signals, assuming noise independence */
    logL = loglikelihoodTDIA + loglikelihoodTDIE + loglikelihoodTDIT;

    //
    printf("logL in CalculateLogLReIm: %f\n", logL);
  }

  /* Clean up */
  LISASignalReIm_Cleanup(generatedsignal);

  return logL;
}

/***************************** Functions handling the prior ******************************/

/* Function to check that returned parameter values fit in prior boundaries */
int PriorBoundaryCheck(LISAPrior *prior, double *Cube)
{
	if (Cube[0] < prior->comp_min || Cube[0] > prior->comp_max ||
	 	Cube[1] < prior->comp_min || Cube[1] > prior->comp_max)
	 	return 1;

	if (Cube[0] + Cube[1] < prior->mtot_min || Cube[0] + Cube[1] > prior->mtot_max)
		return 1;

	if (Cube[0] < Cube[1] || Cube[0] / Cube[1] > prior->qmax)
		return 1;

	return 0;
}

/* Utility prior functions to convert from Cube to common distributions */

double CubeToFlatPrior(double r, double x1, double x2)
{
    return x1 + r * ( x2 - x1 );
}

double CubeToLogFlatPrior(double r, double x1, double x2)
{
    double lx1, lx2;
    lx1 = log( x1 );
    lx2 = log( x2 );
    return exp( lx1 + r * ( lx2 - lx1 ) );
}

double CubeToPowerPrior(double p, double r, double x1, double x2)
{
    double pp = p + 1.0;
    return pow(r * pow(x2, pp) + (1.0 - r) * pow(x1, pp), 1.0 / pp);
}

double CubeToGaussianPrior(double r, double mean, double sigma)
{
    return gsl_cdf_gaussian_Pinv(r,sigma) + mean;
}

double CubeToSinPrior(double r, double x1, double x2)
{
    return acos((1.0-r)*cos(x1)+r*cos(x2));
}

double CubeToCosPrior(double r, double x1, double x2)
{
    return asin((1.0-r)*sin(x1)+r*sin(x2));
}
