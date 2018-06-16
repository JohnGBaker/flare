// This is a general wrapper for black hole merger waveforms.
// The goal is to provide a relatively abstract interface which
// is independent of which instrument (LISA/LIGO/...) is assumed
// and of the specific waveform details, EOB, Taylor, Phenom,...
// If and when the code is converted to python, this may translate
// naturally to an interface in python.
//
// Authors: John Baker(2018)

#include "BBHWaveformFD.h"

//Function for parsing BBH waveform args:
void explain_BBHWaveformFDargs(char *string){
  string="";
};

//Function for parsing BBH waveform args:
int parse_BBHWaveformFDargs(int argc, char **argv, int *iarg, BBHWaveformFDargs *args){
  int nbmode;
  double fRef;
  int TF2extend;
  double TF2ext_Mfmatch;
};

//Waveform generation wrapper functions

//Generate a non-spinning EOB-ROM BBH waveform:
//
int generateEOBNonSpinBBHWaveformFD(
			ListmodesCAmpPhaseFrequencySeries** listWFmodes,   //output waveform in CAmpPhase freq series
			BBHWaveformFDparams* params,                      //input params
			SignalFraming* frame,                             //input info about signal frame
			BBHWaveformFDargs* args                           //input waveform model options
			){
  int ret;

  /* Starting frequency corresponding to duration of observation deltatobs */
  double fstartobs = 0.;
  if(!(frame->tmin==0.))fstartobs = Newtonianfoft(params->m1, params->m2, -frame->tmin/YRSID_SI);
  printf("fstartobs=%g,-frame->tmin=%g\n",fstartobs,-frame->tmin);
  /* Generate the waveform with the ROM */
  /* NOTE: SimEOBNRv2HMROM accepts masses and distances in SI units, whereas BBHWaveformFD params are in solar masses and Mpc */
  if(!(args->TF2extend)) {
    //printf("Not Extending signal waveform.  Mfmatch=%g\n",globalparams->Mfmatch);
    ret = SimEOBNRv2HMROM(listWFmodes, args->nbmode, params->tRef, params->phiRef, args->fRef, (params->m1)*MSUN_SI, (params->m2)*MSUN_SI, (params->distance)*1e6*PC_SI);
  } else {
    //printf("Extending signal waveform.  Mfmatch=%g\n",globalparams->Mfmatch);
    /* If extending, take into account signal framing */
    ret = SimEOBNRv2HMROMExtTF2(listWFmodes, args->nbmode, args->TF2ext_Mfmatch, fmax(fstartobs, frame->fmin), 0, params->tRef, params->phiRef, args->fRef, (params->m1)*MSUN_SI, (params->m2)*MSUN_SI, (params->distance)*1e6*PC_SI);
    //ret = SimEOBNRv2HMROMExtTF2(&listROM, params->nbmode, globalparams->Mfmatch, fmax(fstartobs, globalparams->minf), 0, params->tRef - injectedparams->tRef, params->phiRef, globalparams->fRef, (params->m1)*MSUN_SI, (params->m2)*MSUN_SI, (params->distance)*1e6*PC_SI);
  }
  if(ret==FAILURE){
    //printf("LISAGenerateSignalCAmpPhase: Generation of ROM for injection failed!\n");
    return FAILURE;
  }
  return SUCCESS;
};

//Generate a BBH waveform:
//
int generateBBHWaveformFD(
			int wfType,                                       //which type of waveform to generate?
			ListmodesCAmpPhaseFrequencySeries** listWFmodes,   //output waveform in CAmpPhase freq series
			BBHWaveformFDparams* params,                      //input params
			SignalFraming* frame,                             //input info about signal frame
			BBHWaveformFDargs* args                           //input waveform model options
			){
  if(wfType==0)//Non-spining EOB ROM type waveform 
    return generateEOBNonSpinBBHWaveformFD(listWFmodes, params, frame, args);
  else{
    printf("BBHWaveformFD:generateBBHWaveform: wfType not recognized.\n");
    exit(1);
  }
  return FAILURE;
};


