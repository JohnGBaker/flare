#include "LLVutils.h"

/************ Global Parameters ************/

LLVParams* injectedparams = NULL;
LLVGlobalParams* globalparams = NULL;
LLVPrior* priorParams = NULL;

/************ Functions to initalize and clean up structure for the signals ************/

void LLVSignal_Cleanup(LLVSignal* signal) {
  if(signal->LHOSignal) ListmodesCAmpPhaseFrequencySeries_Destroy(signal->LHOSignal);
  if(signal->LLOSignal) ListmodesCAmpPhaseFrequencySeries_Destroy(signal->LLOSignal);
  if(signal->VIRGOSignal) ListmodesCAmpPhaseFrequencySeries_Destroy(signal->VIRGOSignal);
  free(signal);
}

void LLVSignal_Init(LLVSignal** signal) {
  if(!signal) exit(1);
  /* Create storage for structures */
  if(!*signal) *signal = malloc(sizeof(LLVSignal));
  else
  {
    LLVSignal_Cleanup(*signal);
  }
  (*signal)->LHOSignal = NULL;
  (*signal)->LLOSignal = NULL;
  (*signal)->VIRGOSignal = NULL;
}

void LLVSignalCAmpPhase_Cleanup(LLVSignalCAmpPhase* signal) {
  if(signal->LHOSignal) ListmodesCAmpPhaseFrequencySeries_Destroy(signal->LHOSignal);
  if(signal->LLOSignal) ListmodesCAmpPhaseFrequencySeries_Destroy(signal->LLOSignal);
  if(signal->VIRGOSignal) ListmodesCAmpPhaseFrequencySeries_Destroy(signal->VIRGOSignal);
  free(signal);
}

void LLVSignalCAmpPhase_Init(LLVSignalCAmpPhase** signal) {
  if(!signal) exit(1);
  /* Create storage for structures */
  if(!*signal) *signal = malloc(sizeof(LLVSignalCAmpPhase));
  else
  {
    LLVSignalCAmpPhase_Cleanup(*signal);
  }
  (*signal)->LHOSignal = NULL;
  (*signal)->LLOSignal = NULL;
  (*signal)->VIRGOSignal = NULL;
}

void LLVInjectionCAmpPhase_Cleanup(LLVInjectionCAmpPhase* signal) {
  if(signal->LHOSplines) ListmodesCAmpPhaseSpline_Destroy(signal->LHOSplines);
  if(signal->LLOSplines) ListmodesCAmpPhaseSpline_Destroy(signal->LLOSplines);
  if(signal->VIRGOSplines) ListmodesCAmpPhaseSpline_Destroy(signal->VIRGOSplines);
  free(signal);
}

void LLVInjectionCAmpPhase_Init(LLVInjectionCAmpPhase** signal) {
  if(!signal) exit(1);
  /* Create storage for structures */
  if(!*signal) *signal = malloc(sizeof(LLVInjectionCAmpPhase));
  else
  {
    LLVInjectionCAmpPhase_Cleanup(*signal);
  }
  (*signal)->LHOSplines = NULL;
  (*signal)->LLOSplines = NULL;
  (*signal)->VIRGOSplines = NULL;
}

void LLVSignalReIm_Cleanup(LLVSignalReIm* signal) {
  if(signal->LHOSignal) ReImFrequencySeries_Cleanup(signal->LHOSignal);
  if(signal->LLOSignal) ReImFrequencySeries_Cleanup(signal->LLOSignal);
  if(signal->VIRGOSignal) ReImFrequencySeries_Cleanup(signal->VIRGOSignal);
  free(signal);
}

void LLVSignalReIm_Init(LLVSignalReIm** signal) {
  if(!signal) exit(1);
  /* Create storage for structures */
  if(!*signal) *signal = malloc(sizeof(LLVSignalReIm));
  else
  {
    LLVSignalReIm_Cleanup(*signal);
  }
  (*signal)->LHOSignal = NULL;
  (*signal)->LLOSignal = NULL;
  (*signal)->VIRGOSignal = NULL;
}

void LLVInjectionReIm_Cleanup(LLVInjectionReIm* signal) {
  if(signal->LHOSignal) ReImFrequencySeries_Cleanup(signal->LHOSignal);
  if(signal->LLOSignal) ReImFrequencySeries_Cleanup(signal->LLOSignal);
  if(signal->VIRGOSignal) ReImFrequencySeries_Cleanup(signal->VIRGOSignal);
  if(signal->freq) gsl_vector_free(signal->freq);
  if(signal->noisevaluesLLO) gsl_vector_free(signal->noisevaluesLLO);
  if(signal->noisevaluesLHO) gsl_vector_free(signal->noisevaluesLHO);
  if(signal->noisevaluesVIRGO) gsl_vector_free(signal->noisevaluesVIRGO);
  free(signal);
}

void LLVInjectionReIm_Init(LLVInjectionReIm** signal) {
  if(!signal) exit(1);
  /* Create storage for structures */
  if(!*signal) *signal = malloc(sizeof(LLVInjectionReIm));
  else
  {
    LLVInjectionReIm_Cleanup(*signal);
  }
  (*signal)->LHOSignal = NULL;
  (*signal)->LLOSignal = NULL;
  (*signal)->VIRGOSignal = NULL;
  (*signal)->freq = NULL;
  (*signal)->noisevaluesLHO = NULL;
  (*signal)->noisevaluesLLO = NULL;
  (*signal)->noisevaluesVIRGO = NULL;
}

/************ Functions for LLV parameters, injection, likelihood, prior ************/

/* Parse command line to initialize LLVParams, LLVPrior, and LLVRunParams objects */
void parse_args_LLV(ssize_t argc, char **argv,
    LLVParams* params,
    LLVGlobalParams* globalparams,
    LLVPrior *prior,
    LLVRunParams *run)
{
    char help[] = "\
**********************************************************************\n\
LLVInference by Sylvain Marsat, John Baker, and Philip Graff\n\
Copyright July 2015\n\
**********************************************************************\n\
\n\
This program performs rapid parameter estimation for LIGO and LLV CBC sources in the no-noise case.\n\
Arguments are as follows:\n\
\n\
--------------------------------------------------\n\
----- Injected Signal Parameters -----------------\n\
--------------------------------------------------\n\
 --tRef                Time at reference frequency (sec, default J2000.0 GPS epoch)\n\
 --phiRef              Orbital phase at reference frequency (radians, default=0)\n\
 --m1                  Component mass 1 in Solar masses (larger, default=20)\n\
 --m2                  Component mass 2 in Solar masses (smaller, default=10)\n\
 --distance            Distance to source in Mpc (default=100)\n\
 --ra                  Right ascension of source sky location (radians, default=0)\n\
 --dec                 Declination of source sky location (radians, default=0)\n\
 --inclination         Inclination of source orbital plane to observer line of sight\n\
                       (radians, default=PI/3)\n\
 --polarization        Polarization of source (radians, default=0)\n\
 --fRef                Reference frequency (Hz, default=0 which is interpreted as Mf=0.14)\n\
 --nbmode              Number of modes of radiation to use (1-5, default=5)\n\
 --snr                 Use a target network SNR for the injection by rescaling distance\n\
                       (default=None)\n\
\n\
-----------------------------------------------------------------\n\
----- Global Waveform/Inner products Parameters -----------------\n\
-----------------------------------------------------------------\n\
 --fRef                Reference frequency (Hz, default=0, interpreted as Mf=0.14)\n\
 --minf                Minimal frequency (Hz, default=0) - when set to 0, use the first frequency covered by the noise data of the detector\n\
 --maxf                Maximal frequency (Hz, default=0) - when set to 0, use the last frequency covered by the noise data of the detector\n\
 --nbmodeinj           Number of modes of radiation to use for the injection (1-5, default=5)\n\
 --nbmodetemp          Number of modes of radiation to use for the templates (1-5, default=5)\n\
 --tagint              Tag choosing the integrator: 0 for Fresnel (default), 1 for linear integration\n\
 --tagnetwork          Tag choosing the network of detectors to use (default LHV)\n\
 --nbptsoverlap        Number of points to use for linear integration (default 32768)\n\
\n\
--------------------------------------------------\n\
----- Prior Boundary Settings --------------------\n\
--------------------------------------------------\n\
 --deltaT              Half-width of time prior (sec, default=0.2)\n\
 --comp-min            Minimum component mass in Solar masses (default=4)\n\
 --comp-max            Maximum component mass in Solar masses (default=50)\n\
 --mtot-min            Minimum total mass in Solar masses (default=8)\n\
 --mtot-max            Maximum total mass in Solar masses (default=100)\n\
 --q-max               Maximum mass ratio, m1/m2 (default=11.98, minimum is 1)\n\
 --dist-min            Minimum distance to source (Mpc, default=1)\n\
 --dist-max            Maximum distance to source (Mpc, default=10000)\n\
 --rescale-distprior   In case a target SNR is given with --snr, rescale dist-min and dist-max accordingly\n\
Parameters ra, dec, phase, pol, inc can also ge given min and max values (for testing)\n\
Syntax: --PARAM-min\n\
\n\
--------------------------------------------------\n\
----- Fix Parameters In Sampling -----------------\n\
--------------------------------------------------\n\
 --pin-PARAM           Pin indicated parameter to injected value\n\
 --fix-PARAM           Fix indicated parameter to specified value\n\
 Available parameter names are:\n\
   m1          Mass 1 (MSol)\n\
   m2          Mass 2 (MSol)\n\
   dist        Distance (luminosity, Mpc)\n\
   time        Reference time (GPS sec)\n\
   phase       Reference orbital phase (rad)\n\
   ra          Right ascension (rad)\n\
   dec         Declination (rad)\n\
   inc         Inclination of orbital plane to observer (rad)\n\
   pol         Polarization (rad)\n\
 Note: --pin-PARAM overrides --fix-PARAM\n\
\n\
--------------------------------------------------\n\
----- BAMBI Sampler Settings ---------------------\n\
--------------------------------------------------\n\
 --eff                 Target efficiency of sampling (default=0.1)\n\
 --tol                 Tolerance for evidence calculation convergence (default=0.5)\n\
 --nlive               Number of live points for sampling (default=1000)\n\
 --bambi               Use BAMBI's neural network logL learning (no option, default off)\n\
 --resume              Resume from a previous run (no option, default off)\n\
 --maxiter             Maximum number of iterations - if 0, use convergence criterion to stop (default 0)\n\
 --outroot             Root for output files (default='chains/LLVinference_')\n\
 --netfile             Neural network settings file if using --bambi (default='LLVinference.inp')\n\
 --mmodal              Use multimodal decomposition (no option, default off)\n\
 --maxcls              Max number of modes in multimodal decomposition (default 1)\n\
 --nclspar             Number of parameters to use for multimodal decomposition - in the order of the cube (default 1)\n\
 --ztol                In multimodal decomposition, modes with lnZ lower than ztol are ignored (default -1e90)\n\
 --seed                Seed the inference by setting one of the live points to the injection (no option, default off)\n\
\n";

    ssize_t i;

    /* set default values for the injection params */
    params->tRef = EPOCH_J2000_0_GPS;
    params->phiRef = 0.;
    params->m1 = 20.;
    params->m2 = 10.;
    params->distance = 100.;
    params->ra = 0.;
    params->dec = 0.;
    params->inclination = PI/3.;
    params->polarization = 0.;
    params->nbmode = 5;

    /* set default values for the global params */
    globalparams->fRef = 0.;
    globalparams->minf = 10.;
    globalparams->maxf = 4096.;
    globalparams->nbmodeinj = 5;
    globalparams->nbmodetemp = 5;
    globalparams->tagint = 0;
    globalparams->tagnetwork = LHV;
    globalparams->nbptsoverlap = 32768;

    /* set default values for the prior limits */
    prior->deltaT = 0.1;
    prior->comp_min = 4.0;
    prior->comp_max = 100.0;
    prior->mtot_min = 8.0;
    prior->mtot_max = 200.0;
    prior->qmax = 11.98;
    prior->dist_min = 1.0;
    prior->dist_max = 4000.0;
    prior->ra_min = 0.;
    prior->ra_max = 2*PI;
    prior->dec_min = -PI/2.;
    prior->dec_max = PI/2.;
    prior->phase_min = 0.;
    prior->phase_max = 2*PI;
    prior->pol_min = 0.;
    prior->pol_max = 2*PI;
    prior->inc_min = 0.;
    prior->inc_max = PI;
    prior->fix_m1 = NAN;
    prior->fix_m2 = NAN;
    prior->fix_dist = NAN;
    prior->fix_time = NAN;
    prior->fix_phase = NAN;
    prior->fix_pol = NAN;
    prior->fix_ra = NAN;
    prior->fix_dec = NAN;
    prior->fix_inc = NAN;
    prior->pin_m1 = 0;
    prior->pin_m2 = 0;
    prior->pin_dist = 0;
    prior->pin_time = 0;
    prior->pin_phase = 0;
    prior->pin_pol = 0;
    prior->pin_ra = 0;
    prior->pin_dec = 0;
    prior->pin_inc = 0;
    prior->snr_target = NAN;
    prior->rescale_distprior = 0;
    prior->flat_distprior = 0;

    /* set default values for the run settings */
    run->eff = 0.1;
    run->tol = 0.5;
    run->nlive = 1000;
    strcpy(run->outroot, "chains/LLVinference_");
    run->bambi = 0;
    run->resume = 0;
    run->maxiter = 0;
    strcpy(run->netfile, "LLVinference.inp");
    run->mmodal = 0;
    run->maxcls = 1;
    run->nclspar = 1;
    run->ztol = -1e90;
    run->seed = 0;

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
        } else if (strcmp(argv[i], "--ra") == 0) {
            params->ra = atof(argv[++i]);
        } else if (strcmp(argv[i], "--dec") == 0) {
            params->dec = atof(argv[++i]);
        } else if (strcmp(argv[i], "--inclination") == 0) {
            params->inclination = atof(argv[++i]);
        } else if (strcmp(argv[i], "--polarization") == 0) {
            params->polarization = atof(argv[++i]);
        } else if (strcmp(argv[i], "--fRef") == 0) {
            globalparams->fRef = atof(argv[++i]);
        } else if (strcmp(argv[i], "--minf") == 0) {
            globalparams->minf = atof(argv[++i]);
        } else if (strcmp(argv[i], "--maxf") == 0) {
            globalparams->maxf = atof(argv[++i]);
        } else if (strcmp(argv[i], "--nbmodeinj") == 0) {
            globalparams->nbmodeinj = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--nbmodetemp") == 0) {
            globalparams->nbmodetemp = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--tagint") == 0) {
            globalparams->tagint = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--tagnetwork") == 0) {
            globalparams->tagnetwork = ParseNetworktag(argv[++i]);
        } else if (strcmp(argv[i], "--nbptsoverlap") == 0) {
            globalparams->nbptsoverlap = atoi(argv[++i]);
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
        } else if (strcmp(argv[i], "--ra-min") == 0) {
            prior->ra_min = atof(argv[++i]);
        } else if (strcmp(argv[i], "--ra-max") == 0) {
            prior->ra_max = atof(argv[++i]);
        } else if (strcmp(argv[i], "--dec-min") == 0) {
            prior->dec_min = atof(argv[++i]);
        } else if (strcmp(argv[i], "--dec-max") == 0) {
            prior->dec_max = atof(argv[++i]);
        } else if (strcmp(argv[i], "--phase-min") == 0) {
            prior->phase_min = atof(argv[++i]);
        } else if (strcmp(argv[i], "--phase-max") == 0) {
            prior->phase_max = atof(argv[++i]);
        } else if (strcmp(argv[i], "--pol-min") == 0) {
            prior->pol_min = atof(argv[++i]);
        } else if (strcmp(argv[i], "--pol-max") == 0) {
            prior->pol_max = atof(argv[++i]);
        } else if (strcmp(argv[i], "--inc-min") == 0) {
            prior->inc_min = atof(argv[++i]);
        } else if (strcmp(argv[i], "--inc-max") == 0) {
            prior->inc_max = atof(argv[++i]);
        } else if (strcmp(argv[i], "--fix-m1") == 0) {
            prior->fix_m1 = atof(argv[++i]);
        } else if (strcmp(argv[i], "--fix-m2") == 0) {
            prior->fix_m2 = atof(argv[++i]);
        } else if (strcmp(argv[i], "--fix-dist") == 0) {
            prior->fix_dist = atof(argv[++i]);
        } else if (strcmp(argv[i], "--fix-ra") == 0) {
            prior->fix_ra = atof(argv[++i]);
        } else if (strcmp(argv[i], "--fix-dec") == 0) {
            prior->fix_dec = atof(argv[++i]);
        } else if (strcmp(argv[i], "--fix-time") == 0) {
            prior->fix_time = atof(argv[++i]);
        } else if (strcmp(argv[i], "--fix-phase") == 0) {
            prior->fix_phase = atof(argv[++i]);
        } else if (strcmp(argv[i], "--fix-inc") == 0) {
            prior->fix_inc = atof(argv[++i]);
        } else if (strcmp(argv[i], "--fix-pol") == 0) {
            prior->fix_pol = atof(argv[++i]);
        } else if (strcmp(argv[i], "--pin-m1") == 0) {
            prior->pin_m1 = 1;
        } else if (strcmp(argv[i], "--pin-m2") == 0) {
            prior->pin_m2 = 1;
        } else if (strcmp(argv[i], "--pin-dist") == 0) {
            prior->pin_dist = 1;
        } else if (strcmp(argv[i], "--pin-ra") == 0) {
            prior->pin_ra = 1;
        } else if (strcmp(argv[i], "--pin-dec") == 0) {
            prior->pin_dec = 1;
        } else if (strcmp(argv[i], "--pin-time") == 0) {
            prior->pin_time = 1;
        } else if (strcmp(argv[i], "--pin-phase") == 0) {
            prior->pin_phase = 1;
        } else if (strcmp(argv[i], "--pin-inc") == 0) {
            prior->pin_inc = 1;
        } else if (strcmp(argv[i], "--pin-pol") == 0) {
            prior->pin_pol = 1;
        } else if (strcmp(argv[i], "--snr") == 0) {
          prior->snr_target = atof(argv[++i]);
        } else if (strcmp(argv[i], "--rescale-distprior") == 0) {
          prior->rescale_distprior = 1;
        } else if (strcmp(argv[i], "--flat-distprior") == 0) {
            prior->flat_distprior = 1;
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
        } else if (strcmp(argv[i], "--maxiter") == 0) {
            run->maxiter = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--outroot") == 0) {
            strcpy(run->outroot, argv[++i]);
        } else if (strcmp(argv[i], "--netfile") == 0) {
            strcpy(run->netfile, argv[++i]);
        } else if (strcmp(argv[i], "--mmodal") == 0) {
            run->mmodal = 1;
        } else if (strcmp(argv[i], "--maxcls") == 0) {
            run->maxcls = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--nclspar") == 0) {
            run->nclspar = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--ztol") == 0) {
            run->ztol = atof(argv[++i]);
        } else if (strcmp(argv[i], "--seed") == 0) {
            run->seed = 1;
        } else {
            printf("Error: invalid option: %s\n", argv[i]);
            goto fail;
        }
    }

    return;

    fail:
    exit(1);
}

/************************* Functions to generate signals and compute likelihoods **************************/

/* Function generating a LLV signal as a list of modes in CAmp/Phase form, from LLV parameters */
int LLVGenerateSignalCAmpPhase(
  struct tagLLVParams* params,            /* Input: set of LLV parameters of the signal */
  struct tagLLVSignalCAmpPhase* signal)   /* Output: structure for the generated signal */
{
//
//printf("In LLVGenerateSignalCAmpPhase\n");

  int ret;
  ListmodesCAmpPhaseFrequencySeries* listROM = NULL;
  ListmodesCAmpPhaseFrequencySeries* listDet1 = NULL;
  ListmodesCAmpPhaseFrequencySeries* listDet2 = NULL;
  ListmodesCAmpPhaseFrequencySeries* listDet3 = NULL;

  /* Checking that the global injectedparams has been set up */
  if (!injectedparams) {
    printf("Error: when calling LLVGenerateSignal, injectedparams points to NULL.\n");
    exit(1);
  }
  /* Should add more error checking ? */
  /* Generate the waveform with the ROM */
  /* Note: SimEOBNRv2HMROM accepts masses and distances in SI units, whereas LLV params is in solar masses and Mpc */

//
//printf("params: (%d, %16e, %16e, %16e, %16e, %16e, %16e)\n", params->nbmode, params->tRef - injectedparams->tRef, params->phiRef, globalparams->fRef, (params->m1), (params->m2), (params->distance));

  ret = SimEOBNRv2HMROM(&listROM, params->nbmode, params->tRef - injectedparams->tRef, params->phiRef, globalparams->fRef, (params->m1)*MSUN_SI, (params->m2)*MSUN_SI, (params->distance)*1e6*PC_SI);

  //
  //printf("after SimEOBNRv2HMROM\n");

  /* If the ROM waveform generation failed (e.g. parameters were out of bounds) return FAILURE */
  if(ret==FAILURE) return FAILURE;

  /* Process the waveform through the LLV response */
  //WARNING: tRef is ignored for now, i.e. set to 0
  //TESTING
  //clock_t tbeg, tend;
  //tbeg = clock();
  LLVSimFDResponse3Det(&listROM, &listDet1, &listDet2, &listDet3, params->tRef, params->ra, params->dec, params->inclination, params->polarization, globalparams->tagnetwork);
  //tend = clock();
  //printf("time LLVSimFDResponse: %g\n", (double) (tend-tbeg)/CLOCKS_PER_SEC);
  //

  //
  //printf("after LLVSimFDResponse3Det\n");

  /* Pre-interpolate the injection, building the spline matrices */
  ListmodesCAmpPhaseSpline* listsplinesgen1 = NULL;
  ListmodesCAmpPhaseSpline* listsplinesgen2 = NULL;
  ListmodesCAmpPhaseSpline* listsplinesgen3 = NULL;
  BuildListmodesCAmpPhaseSpline(&listsplinesgen1, listDet1);
  BuildListmodesCAmpPhaseSpline(&listsplinesgen2, listDet2);
  BuildListmodesCAmpPhaseSpline(&listsplinesgen3, listDet3);

  //
  //printf("after BuildListmodesCAmpPhaseSpline\n");

  /* Precompute the inner product (h|h) */
  //TESTING
  //tbeg = clock();
  /* Note: we ignore fstartobs and assume (for the noises) that the detectors are LHO, LLO and VIRGO */
  double Det123hh = FDListmodesFresnelOverlap3Chan(listDet1, listDet2, listDet3, listsplinesgen1, listsplinesgen2, listsplinesgen3, NoiseSnLHO, NoiseSnLLO, NoiseSnVIRGO, globalparams->minf, globalparams->maxf, 0., 0.);
  //tend = clock();
  //printf("time SNRs: %g\n", (double) (tend-tbeg)/CLOCKS_PER_SEC);

  //
  //printf("after FDListmodesFresnelOverlap3Chan\n");

  /* Output and clean up */
  signal->LHOSignal = listDet1;
  signal->LLOSignal = listDet2;
  signal->VIRGOSignal = listDet3;
  signal->LLVhh = Det123hh;

  ListmodesCAmpPhaseFrequencySeries_Destroy(listROM);
  ListmodesCAmpPhaseSpline_Destroy(listsplinesgen1);
  ListmodesCAmpPhaseSpline_Destroy(listsplinesgen2);
  ListmodesCAmpPhaseSpline_Destroy(listsplinesgen3);

  //
  //printf("after cleanup\n");

  return SUCCESS;
}

/* Function generating a LLV signal as a list of modes in CAmp/Phase form, from LLV parameters */
int LLVGenerateInjectionCAmpPhase(
  struct tagLLVParams* params,       /* Input: set of LLV parameters of the signal */
  struct tagLLVInjectionCAmpPhase* signal)   /* Output: structure for the injected signal */
{
  int ret;
  ListmodesCAmpPhaseFrequencySeries* listROM = NULL;
  ListmodesCAmpPhaseFrequencySeries* listDet1 = NULL;
  ListmodesCAmpPhaseFrequencySeries* listDet2 = NULL;
  ListmodesCAmpPhaseFrequencySeries* listDet3 = NULL;

  /* Should add more error checking ? */
  /* Generate the waveform with the ROM */
  /* Note: SimEOBNRv2HMROM accepts masses and distances in SI units, whereas LLV params is in solar masses and Mpc */
  ret = SimEOBNRv2HMROM(&listROM, params->nbmode, params->tRef - injectedparams->tRef, params->phiRef, globalparams->fRef, (params->m1)*MSUN_SI, (params->m2)*MSUN_SI, (params->distance)*1e6*PC_SI);

  /* If the ROM waveform generation failed (e.g. parameters were out of bounds) return FAILURE */
  if(ret==FAILURE) return FAILURE;

  /* Process the waveform through the LLV response */
  //WARNING: tRef is ignored for now, i.e. set to 0
  //TESTING
  //clock_t tbeg, tend;
  //tbeg = clock();
  LLVSimFDResponse3Det(&listROM, &listDet1, &listDet2, &listDet3, params->tRef, params->ra, params->dec, params->inclination, params->polarization, globalparams->tagnetwork);
  //tend = clock();
  //printf("time LLVSimFDResponse: %g\n", (double) (tend-tbeg)/CLOCKS_PER_SEC);
  //

  /* Pre-interpolate the injection, building the spline matrices */
  ListmodesCAmpPhaseSpline* listsplinesinj1 = NULL;
  ListmodesCAmpPhaseSpline* listsplinesinj2 = NULL;
  ListmodesCAmpPhaseSpline* listsplinesinj3 = NULL;
  BuildListmodesCAmpPhaseSpline(&listsplinesinj1, listDet1);
  BuildListmodesCAmpPhaseSpline(&listsplinesinj2, listDet2);
  BuildListmodesCAmpPhaseSpline(&listsplinesinj3, listDet3);

  /* Precompute the inner product (h|h) - we ignore deltatobs */
  /* Note: for the noise functions we assume the detectors are L,H,V */
  //TESTING
  //tbeg = clock();
  double Det123ss = FDListmodesFresnelOverlap3Chan(listDet1, listDet2, listDet3, listsplinesinj1, listsplinesinj2, listsplinesinj3, NoiseSnLHO, NoiseSnLLO, NoiseSnVIRGO, globalparams->minf, globalparams->maxf, 0., 0.);
  //tend = clock();
  //printf("time SNRs: %g\n", (double) (tend-tbeg)/CLOCKS_PER_SEC);

  /* Output and clean up */
  signal->LHOSplines = listsplinesinj1;
  signal->LLOSplines = listsplinesinj2;
  signal->VIRGOSplines = listsplinesinj3;
  signal->LLVss = Det123ss;

  ListmodesCAmpPhaseFrequencySeries_Destroy(listROM);
  ListmodesCAmpPhaseFrequencySeries_Destroy(listDet1);
  ListmodesCAmpPhaseFrequencySeries_Destroy(listDet2);
  ListmodesCAmpPhaseFrequencySeries_Destroy(listDet3);

  return SUCCESS;
}

int LLVGenerateSignal(
  struct tagLLVParams* params,   /* Input: set of LLV parameters of the signal */
  struct tagLLVSignal* signal)   /* Output: structure for the generated signal */
{
  int ret;
  ListmodesCAmpPhaseFrequencySeries* listROM = NULL;
  ListmodesCAmpPhaseFrequencySeries* listLHO = NULL;
  ListmodesCAmpPhaseFrequencySeries* listLLO = NULL;
  ListmodesCAmpPhaseFrequencySeries* listVIRGO = NULL;

  /* Checking that the global injectedparams has been set up */
  if (!injectedparams) {
    printf("Error: when calling LLVGenerateSignal, injectedparams points to NULL.\n");
    exit(1);
  }
  /* Should add more error checking ? */
  /* Generate the waveform with the ROM */
  /* Note: SimEOBNRv2HMROM accepts masses and distances in SI units, whereas LLV params is in solar masses and Mpc */
  ret = SimEOBNRv2HMROM(&listROM, params->nbmode, params->tRef - injectedparams->tRef, params->phiRef, globalparams->fRef, (params->m1)*MSUN_SI, (params->m2)*MSUN_SI, (params->distance)*1e6*PC_SI);

  /* If the ROM waveform generation failed (e.g. parameters were out of bounds) return FAILURE */
  if(ret==FAILURE) return FAILURE;

  /* Process the waveform through the LLV response */
  LLVSimFDResponse3Det(&listROM, &listLHO, &listLLO, &listVIRGO, params->tRef, params->ra, params->dec, params->inclination, params->polarization, globalparams->tagnetwork);
  // LLVSimFDResponse(&listROM, &listLHO, params->tRef, params->ra, params->dec, params->inclination, params->polarization, LHO);
  // LLVSimFDResponse(&listROM, &listLLO, params->tRef, params->ra, params->dec, params->inclination, params->polarization, LLO);
  // LLVSimFDResponse(&listROM, &listVIRGO, params->tRef, params->ra, params->dec, params->inclination, params->polarization, VIRGO);

  /* Precompute the inner products (h|h) - 3 last args: ignore fstartobs, and use wip for overlap */
  double fLowLHO, fHighLHO, fLowLLO, fHighLLO, fLowVIRGO, fHighVIRGO;
  fLowLHO = fmax(globalparams->minf, __LLVSimFD_LHONoise_fLow);
  fLowLLO = fmax(globalparams->minf, __LLVSimFD_LLONoise_fLow);
  fLowVIRGO = fmax(globalparams->minf, __LLVSimFD_VIRGONoise_fLow);
  if(globalparams->maxf==0) {
    fHighLHO = __LLVSimFD_LHONoise_fHigh;
    fHighLLO = __LLVSimFD_LLONoise_fHigh;
    fHighVIRGO = __LLVSimFD_VIRGONoise_fHigh;
  }
  else {
    fHighLHO = fmin(globalparams->maxf, __LLVSimFD_LHONoise_fHigh);
    fHighLLO = fmin(globalparams->maxf, __LLVSimFD_LLONoise_fHigh);
    fHighVIRGO = fmin(globalparams->maxf, __LLVSimFD_VIRGONoise_fHigh);
  }
  double LHOhh = FDListmodesOverlap(listLHO, listLHO, NoiseSnLHO, fLowLHO, fHighLHO, 0., 0., 0);
  double LLOhh = FDListmodesOverlap(listLLO, listLLO, NoiseSnLLO, fLowLLO, fHighLLO, 0., 0., 0);
  double VIRGOhh = FDListmodesOverlap(listVIRGO, listVIRGO, NoiseSnVIRGO, fLowVIRGO, fHighVIRGO, 0., 0., 0);

  //
  /* Pre-interpolate the injection, building the spline matrices */
  ListmodesCAmpPhaseSpline* listsplinesgen1 = NULL;
  ListmodesCAmpPhaseSpline* listsplinesgen2 = NULL;
  ListmodesCAmpPhaseSpline* listsplinesgen3 = NULL;
  BuildListmodesCAmpPhaseSpline(&listsplinesgen1, listLHO);
  BuildListmodesCAmpPhaseSpline(&listsplinesgen2, listLLO);
  BuildListmodesCAmpPhaseSpline(&listsplinesgen3, listVIRGO);
  double Det123hh = FDListmodesFresnelOverlap3Chan(listLHO, listLLO, listVIRGO, listsplinesgen1, listsplinesgen2, listsplinesgen3, NoiseSnLHO, NoiseSnLLO, NoiseSnVIRGO, globalparams->minf, globalparams->maxf, 0., 0.);

  /* Output and clean up */
  signal->LHOSignal = listLHO;
  signal->LLOSignal = listLLO;
  signal->VIRGOSignal = listVIRGO;
  signal->LHOhh = LHOhh;
  signal->LLOhh = LLOhh;
  signal->VIRGOhh = VIRGOhh;

  ListmodesCAmpPhaseFrequencySeries_Destroy(listROM);
  ListmodesCAmpPhaseSpline_Destroy(listsplinesgen1);
  ListmodesCAmpPhaseSpline_Destroy(listsplinesgen2);
  ListmodesCAmpPhaseSpline_Destroy(listsplinesgen3);
  return SUCCESS;
}

/* Function to check that returned parameter values fit in prior boundaries */
int PriorBoundaryCheck(LLVPrior *prior, double *Cube)
{
	if (Cube[0] < Cube[1])
	 	return 1;

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

/* Log-Likelihood function */

double CalculateLogLCAmpPhase(LLVParams *params, LLVInjectionCAmpPhase* injection)
{
  double logL = -DBL_MAX;
  int ret;

  /* Generating the signal in the three detectors for the input parameters */
  LLVSignalCAmpPhase* generatedsignal = NULL;
  LLVSignalCAmpPhase_Init(&generatedsignal);
  //TESTING
  //clock_t tbeg, tend;
  //tbeg = clock();
  ret = LLVGenerateSignalCAmpPhase(params, generatedsignal);
  //tend = clock();
  //printf("time GenerateSignal: %g\n", (double) (tend-tbeg)/CLOCKS_PER_SEC);
  //

  /* If LLVGenerateSignal failed (e.g. parameters out of bound), silently return -Infinity logL */
  if(ret==FAILURE) {
    logL = -DBL_MAX;
  }
  else if(ret==SUCCESS) {
    /* Computing the likelihood for each detector - fstartobs is ignored, and we assume for the noises that the detectors are LHO, LLO, VIRGO */
    //TESTING
    //tbeg = clock();
    double overlapDet123 = FDListmodesFresnelOverlap3Chan(generatedsignal->LHOSignal, generatedsignal->LLOSignal, generatedsignal->VIRGOSignal, injection->LHOSplines, injection->LLOSplines, injection->VIRGOSplines, NoiseSnLHO, NoiseSnLLO, NoiseSnVIRGO, globalparams->minf, globalparams->maxf, 0., 0.);
    //tend = clock();
    //printf("time Overlaps: %g\n", (double) (tend-tbeg)/CLOCKS_PER_SEC);
    //

    /* Output: value of the loglikelihood for the combined signals, assuming noise independence */
    logL = overlapDet123 - 1./2*(injection->LLVss) - 1./2*(generatedsignal->LLVhh);
  }

  /* Clean up */
  LLVSignalCAmpPhase_Cleanup(generatedsignal);

  return logL;
}

double CalculateLogL(LLVParams *params, LLVInjectionCAmpPhase* injection)
{
  double logL = -DBL_MAX;
  int ret;

  /* Generating the signal in the three detectors for the input parameters */
  LLVSignalCAmpPhase* generatedsignal = NULL;
  LLVSignalCAmpPhase_Init(&generatedsignal);
  ret = LLVGenerateSignalCAmpPhase(params, generatedsignal);

  /* If LLVGenerateSignal failed (e.g. parameters out of bound), silently return -Infinity logL */
  if(ret==FAILURE) {
    logL = -DBL_MAX;
  }
  else if(ret==SUCCESS) {
    /* Computing the likelihood for each detector - last 3 args: ignoring fstartobs, and using wip for the overlap */
    // double fLowLHO, fHighLHO, fLowLLO, fHighLLO, fLowVIRGO, fHighVIRGO;
    // fLowLHO = fmax(globalparams->minf, __LLVSimFD_LHONoise_fLow);
    // fLowLLO = fmax(globalparams->minf, __LLVSimFD_LLONoise_fLow);
    // fLowVIRGO = fmax(globalparams->minf, __LLVSimFD_VIRGONoise_fLow);
    // if(globalparams->maxf==0) {
    //   fHighLHO = __LLVSimFD_LHONoise_fHigh;
    //   fHighLLO = __LLVSimFD_LLONoise_fHigh;
    //   fHighVIRGO = __LLVSimFD_VIRGONoise_fHigh;
    // }
    // else {
    //   fHighLHO = fmin(globalparams->maxf, __LLVSimFD_LHONoise_fHigh);
    //   fHighLLO = fmin(globalparams->maxf, __LLVSimFD_LLONoise_fHigh);
    //   fHighVIRGO = fmin(globalparams->maxf, __LLVSimFD_VIRGONoise_fHigh);
    // }
    // double loglikelihoodLHO = FDLogLikelihood(injection->LHOSignal, generatedsignal->LHOSignal, NoiseSnLHO, fLowLHO, fHighLHO, injection->LHOhh, generatedsignal->LHOhh, 0., 0., 0);
    // double loglikelihoodLLO = FDLogLikelihood(injection->LLOSignal, generatedsignal->LLOSignal, NoiseSnLLO, fLowLLO, fHighLLO, injection->LLOhh, generatedsignal->LLOhh, 0., 0., 0);
    // double loglikelihoodVIRGO = FDLogLikelihood(injection->VIRGOSignal, generatedsignal->VIRGOSignal, NoiseSnVIRGO, fLowVIRGO, fHighVIRGO, injection->VIRGOhh, generatedsignal->VIRGOhh, 0., 0., 0);
    double overlapDet123 = FDListmodesFresnelOverlap3Chan(generatedsignal->LHOSignal, generatedsignal->LLOSignal, generatedsignal->VIRGOSignal, injection->LHOSplines, injection->LLOSplines, injection->VIRGOSplines, NoiseSnLHO, NoiseSnLLO, NoiseSnVIRGO, globalparams->minf, globalparams->maxf, 0., 0.);
    double LLVss = injection->LLVss;
    double LLVhh = generatedsignal->LLVhh;

    /* Output: value of the loglikelihood for the combined signals, assuming noise independence */
    logL = overlapDet123 - 0.5*LLVss - 0.5*LLVhh;
  }

  /* Clean up */
  LLVSignalCAmpPhase_Cleanup(generatedsignal);

  return logL;
}
