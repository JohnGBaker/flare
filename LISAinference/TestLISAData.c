#include "ComputeLISASNR.h"

/************ Parsing arguments function ************/
extern int half_edges;

double tShift=0;

/* Masses are input in solar masses and distances in Mpc - converted in SI for the internals */
static void parse_args_ComputeLISASNR(ssize_t argc, char **argv, ComputeLISASNRparams* params)
{
  char help[] = "\
Developed from ComputeLISASNR by Sylvain Marsat, John Baker, and Philip Graff\n\
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
 --nbptsoverlap        Number of points to use for linear integration (default 32768)\n\
 --indir               Input directory when loading TDI freq series from file\n\
 --infile              Input file name when loading TDI freq series from file\n\
 --loadparamsfile      Option to load physical parameters from file and to output result to file (default false)\n\
 --nlinesparams        Number of lines in params file\n\
 --paramsdir           Directory for input/output file\n\
 --paramsfile          Input file with the parameters\n\
 --outputfile          Output file\n\
------------------------------------\n\
Test Params\n\
------------------------------------\n\
 --tShift              time to offset the data time regisration from the nominal signal registration\n\
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
      params->loadparamsfile = 1;
    } else if (strcmp(argv[i], "--outputfile") == 0) {
      strcpy(params->outputfile, argv[++i]);
    } else if (strcmp(argv[i], "--tShift") == 0) {
      tShift=atof(argv[++i]);
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
  LISAconstellation *variant = &LISAProposal;
  int tagtRefatLISA = 0;
  int tagsimplelikelihood = 0;
  int zerolikelihood = 0;
  int frozenLISA = 0;
  ResponseApproxtag responseapprox = full;

  double SNR = 0;

  /* Initialize structure for parameters */
  ComputeLISASNRparams* params;
  LISAParams* testparams=malloc(sizeof(LISAParams));
  params = (ComputeLISASNRparams*) malloc(sizeof(ComputeLISASNRparams));
  memset(params, 0, sizeof(ComputeLISASNRparams));

  /* Parse commandline to read parameters */
  parse_args_ComputeLISASNR(argc, argv, params);

    /* Initialize the structure for LISAParams and GlobalParams and copy values */
  /* NOTE: injectedparams and globalparams are declared as extern in LISAutils.h, and used by generation functions */
  injectedparams = (LISAParams*) malloc(sizeof(LISAParams));
  memset(injectedparams, 0, sizeof(LISAParams));
  injectedparams->nbmode = params->nbmode;
  if(params->loadparamsfile==0) {
    injectedparams->tRef = params->tRef;
    injectedparams->phiRef = params->phiRef;
    injectedparams->m1 = params->m1;
    injectedparams->m2 = params->m2;
    injectedparams->distance = params->distance;
    injectedparams->lambda = params->lambda;
    injectedparams->beta = params->beta;
    injectedparams->inclination = params->inclination;
    injectedparams->polarization = params->polarization;
  } else {
    printf(" params from file %s\n",params->paramsfile);
    gsl_matrix* inmatrix =  gsl_matrix_alloc(1, 9);
    Read_Text_Matrix(params->paramsdir, params->paramsfile, inmatrix);
    /* Format (same as in the internals): m1, m2, tRef, dist, phase, inc, lambda, beta, pol, SNR */
    injectedparams->m1 = gsl_matrix_get(inmatrix, 0, 0);
    injectedparams->m2 = gsl_matrix_get(inmatrix, 0, 1);
    injectedparams->tRef = gsl_matrix_get(inmatrix, 0, 2);
    injectedparams->distance = gsl_matrix_get(inmatrix, 0, 3);
    injectedparams->phiRef = gsl_matrix_get(inmatrix, 0, 4);
    injectedparams->inclination = gsl_matrix_get(inmatrix, 0, 5);
    injectedparams->lambda = gsl_matrix_get(inmatrix, 0, 6);
    injectedparams->beta = gsl_matrix_get(inmatrix, 0, 7);
    injectedparams->polarization = gsl_matrix_get(inmatrix, 0, 8);
    gsl_matrix_free(inmatrix);
  }
  memcpy(testparams,injectedparams,sizeof(LISAParams));
  
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
  /* Hardcoded */
  globalparams->variant = variant;
  globalparams->tagtRefatLISA = tagtRefatLISA;
  globalparams->frozenLISA = frozenLISA;
  globalparams->tagsimplelikelihood = tagsimplelikelihood;
  globalparams->zerolikelihood = zerolikelihood;
  globalparams->responseapprox = responseapprox;
  
  double DataSNRA2,DataSNRE2,DataSNRT2;
  double ModelSNRA2,ModelSNRE2,ModelSNRT2;
  LISADataFD * data = NULL;
  LISASignalCAmpPhase* sigCAmpPhase = NULL;
	

  /* NOTE: supports only AET(XYZ) (orthogonal) */

  // First step is that we set up the "data" either reading from file or generating anew
  
  if(!(AllowedTDItag(params->tagtdi))) {
    printf("Error in ComputeLISASNR: TDI tag not recognized.\n");
    exit(1);

  } else {

    if(strlen(params->infile)>0) {
      
      Read_LISADataFD(params->indir, params->infile, &data);
      SNR = sqrt(data->TDI123dd);
      /* Print SNR to stdout */
      printf("File Data SNR %.8f\n", SNR);
    }

    else {       /* Generate TDI Data from params */

      printf("Generating data from params\n");
      LISADataFD_Init(&data);
      double M_sec=( injectedparams->m1 + injectedparams->m2 ) * MTSUN_SI;
      double duration=globalparams->deltatobs * YRSID_SI + 1000*M_sec;
      //As a test, to mimic the effect of shifting the FT reference (by tShift) we shift the GW signal reference time and reference phase together.
      //If every thing is correct, then the only difference left between the signal and model should be the timing of the time-dependence of the response
      //Typically that isn't measurable below changes of several arc seconds, so there should be minimal effect below a few thousand seconds.
      //dataphase=phasesig(f)-phasesig(fref)+2pi*(tinj-tinj)*(f-fref)+2*phi0inj-2pi*tshift*f
      //modelphase=phasesig(f)-phasesig(fref)+2pi*(tmodel-tinj)*(f-fref)+2*phi0model
      //dphase=modelphase-datapahse=+2pi*(tmodel-tinj)*(f-fref)+2*phi0model-2*phi0inj+2pi*tshift*f
      //dphase=0 --> (tmodel-tinj)*2*pi*fref +2*(phi0model-phi0imj) = 0
      //         and (tmodel-tinj)+tshift = 0;
      //       ...-->-tshift*pi*fref + phi0model-phi0inj = 0
      double mimicShift=tShift;
      double fRef=params->fRef;
      double fRefm=0.14/M_sec;//This is currently the hard coded phi0-defining freq.
      if(fRef==0 || fRef>fRefm)fRef=fRefm;
      double phiShift=mimicShift*PI*fRef;//phiRef (in EOB...core) seems scaled as orbital freq, not GW freq, thus Pi instead of 2PI
      injectedparams->tRef+=+mimicShift;
      injectedparams->phiRef+=phiShift;
      printf("tshift,phiShift,Msec= on params = %g %g %g\n",tShift,phiShift,M_sec);
      //Now generate the reference signal.  Note that the time-shift implemented is based on the difference between testparams and injectedparams
      //LISAGenerateDataFD(testparams,0,duration+tShift,data);
      LISAGenerateDataFD(injectedparams,0,duration+tShift,data);
      //LISAGenerateDataFD(testparams,0,duration+0*tShift,data);
      SNR = sqrt(data->TDI123dd);
      
      /* Print SNR to stdout */
      printf("Generated data with params:\n");report_LISAParams(injectedparams);
      printf("Generated Data SNR %.8f\n", SNR);

      /* Output Data */
      if(params->outputfile!="")
	Write_LISADataFD(params->paramsdir, params->outputfile, data);
    }

    //Now generate the model signal for comparison. This one will be based on whichever params were last used

    /* Branch between the Fresnel or linear computation */
    if(params->tagint==0) {
      printf("&sigCAmpPhase=%p\n", &sigCAmpPhase);
      printf("*&sigCAmpPhase=%p\n", *&sigCAmpPhase);
      LISASignalCAmpPhase_Init(&sigCAmpPhase);
      clock_t tbeg, tend;
      tbeg = clock();
      LISAGenerateSignalCAmpPhase(injectedparams, sigCAmpPhase);
      tend = clock();
      double t=(double) (tend-tbeg)/CLOCKS_PER_SEC;
      printf("Model generation time=%g\n",t);
      SNR = sqrt(sigCAmpPhase->TDI123hh);
    } else {
      printf("Error in ComputeLISASNR: integration tag not recognized.\n");
      exit(1);
    }

    /* Print SNR to stdout */
    printf("Signal params:\n");report_LISAParams(injectedparams);
    printf("Model Data SNR %.8f\n", SNR);

    
  }

  double logL;
  logL = CalculateLogLDataCAmpPhase(testparams, &data,1);
  printf("Amp/Phase: logL=%g\n",logL);
  half_edges=1;
  logL = CalculateLogLDataCAmpPhase(testparams, &data,1);
  half_edges=0;
  printf("Amp/Phase: logL half=%g\n",logL);
  IntProdStyle=1;
  if(0){
  clock_t tbeg, tend;
  int ncuts=20,ntry=15;
  double cutpower=1.5;
  double times[ncuts];
  double times2[ncuts];
  for(int n=0;n<ncuts;n++)times[n]=times2[n]=0;
  for(int i=0;i<ntry;i++){
    for(int n=0;n<ncuts;n++){
      NratioCut=(int)pow(n,cutpower);
      tbeg = clock();
      logL = CalculateLogLDataCAmpPhase(testparams, &data,1);
      tend = clock();
      double t=(double) (tend-tbeg)/CLOCKS_PER_SEC;
      printf(" t=%g\n",t);
      times[n]+=t;
      times2[n]+=t*t;
    }
  }
  for(int n=0;n<ncuts;n++)
    printf("n=%i nCut = %2i: mean = %10.6g  std = %10.6g\n", n,(int)pow(n,cutpower),times[n]/ntry,sqrt((times2[n]-times[n]*times[n]/ntry)/(ntry-1)));
  }
  
  logL = CalculateLogLDataCAmpPhase(testparams, &data,1);
  printf("Amp/Phase E: logL=%g\n",logL);
  half_edges=1;
  logL = CalculateLogLDataCAmpPhase(testparams, &data,1);
  printf("Amp/Phase B half: logL=%g\n",logL);
  half_edges=0;

  LISADataFD * dataD2 = LISADataFD_Decimate2(data); 
  logL = CalculateLogLDataCAmpPhase(testparams, &dataD2,1);
  printf("Amp/Phase E (decimated by 2): logL=%g\n",logL);
  half_edges=1;
  logL = CalculateLogLDataCAmpPhase(testparams, &dataD2,1);
  half_edges=0;
  printf("Amp/Phase B half (decimated by 2): logL=%g\n",logL);

  LISADataFD * dataD4 = LISADataFD_Decimate2(dataD2); 
  //logL = CalculateLogLDataCAmpPhase(testparams, &dataD4,1);
  //printf("Amp/Phase E (decimated by 4): logL=%g\n",logL);
  //half_edges=1;
  //logL = CalculateLogLDataCAmpPhase(testparams, &dataD4,1);
  //half_edges=0;
  //printf("Amp/Phase B half (decimated by 4): logL=%g\n",logL);

  LISADataFD * dataD8 = LISADataFD_Decimate2(dataD4); 
  //logL = CalculateLogLDataCAmpPhase(testparams, &dataD8,1);
  //printf("Amp/Phase E (decimated by 8): logL=%g\n",logL);
  //half_edges=1;
  //logL = CalculateLogLDataCAmpPhase(testparams, &dataD8,1);
  //half_edges=0;
  //printf("Amp/Phase B half (decimated by 8): logL=%g\n",logL);

  LISADataFD * dataD16 = LISADataFD_Decimate2(dataD8); 
  //logL = CalculateLogLDataCAmpPhase(testparams, &dataD16,1);
  //printf("Amp/Phase E (decimated by 16): logL=%g\n",logL);
  //half_edges=1;
  //logL = CalculateLogLDataCAmpPhase(testparams, &dataD16,1);
  //half_edges=0;
  //printf("Amp/Phase B half (decimated by 16): logL=%g\n",logL);

  LISADataFD * dataD32 = LISADataFD_Decimate2(dataD16); 
  //logL = CalculateLogLDataCAmpPhase(testparams, &dataD32,1);
  //printf("Amp/Phase E (decimated by 32): logL=%g\n",logL);

  LISADataFD * dataD64 = LISADataFD_Decimate2(dataD32); 
  //logL = CalculateLogLDataCAmpPhase(testparams, &dataD64,1);
  //printf("Amp/Phase E (decimated by 64): logL=%g\n",logL);

  LISADataFD *datastack[7]={data,dataD2,dataD4,dataD8,dataD16,dataD32,dataD64};
  logL = CalculateLogLDataCAmpPhase(testparams, datastack,7);
  printf("Amp/Phase E(stacked data): logL=%g\n",logL);
  

  free(params);
  LISADataFD_Cleanup(data);
  LISASignalCAmpPhase_Cleanup(sigCAmpPhase);
}
