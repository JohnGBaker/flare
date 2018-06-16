// This is a general wrapper for black hole merger waveforms.
// The goal is to provide a relatively abstract interface which
// is independent of which instrument (LISA/LIGO/...) is assumed
// and of the specific waveform details, EOB, Taylor, Phenom,...
// If and when the code is converted to python, this may translate
// naturally to an interface in python.
//
// Authors: John Baker(2018)

#ifndef __BBHWAVEFORMFD_H__
#define __BBHWAVEFORMFD_H__ 1

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_cdf.h>

#include "constants.h"
#include "struct.h"
#include "EOBNRv2HMROMstruct.h"
#include "EOBNRv2HMROM.h"
#include "waveform.h"

#if defined(__cplusplus)
extern "C" {
#define complex _Complex
#elif 0
} /* so that editors will match preceding brace */
#endif


//Structure BBH waveform args:
typedef struct{
  int nbmode;
  double fRef;
  int TF2extend;
  double TF2ext_Mfmatch;
} BBHWaveformFDargs;

// Parameters for the generation of a BBH waveform 
typedef struct {
  double tRef;               /* reference time (s) - GPS time at the frequency representing coalescence */
  double phiRef;             /* reference phase (rad) - phase at the frequency representing coalescence (or at fRef if specified) */
  double m1;                 /* mass of companion 1 (solar masses, default 2e6) */
  double m2;                 /* mass of companion 2 (solar masses, default 1e6) */
  double distance;           /* distance of source (Mpc, default 1e3) */
  double lambda;             /* first angle for the position in the sky (rad, default 0) */
  double beta;               /* second angle for the position in the sky (rad, default 0) */
  double inclination;        /* inclination of L relative to line of sight (rad, default PI/3) */
  double polarization;       /* polarization angle (rad, default 0) */
  double chi1x;              /* x-component of scaled spin of BH1 (assume BH are on x-axis at reference moment) */
  double chi1y;              /* y-component of scaled spin of BH1 (assume BH are on x-axis at reference moment) */
  double chi1z;              /* z-component of scaled spin of BH1 (z-axis is orb ang mom dir at reference moment) */
  double chi2x;              /* x-component of scaled spin of BH1 (assume BH are on x-axis at reference moment) */
  double chi2y;              /* y-component of scaled spin of BH1 (assume BH are on x-axis at reference moment) */
  double chi2z;              /* z-component of scaled spin of BH1 (z-axis is orb ang mom dir at reference moment) */
} BBHWaveformFDparams;

//Generate a BBH waveform:
int generateBBHWaveformFD(
			int wfType,                                       //which type of waveform to generate?
			ListmodesCAmpPhaseFrequencySeries** listWFmodes,   //output waveform in CAmpPhase freq series
			BBHWaveformFDparams* params,                      //input params
			SignalFraming* frame,                             //input info about signal frame
			BBHWaveformFDargs* args                           //input waveform model options
			);

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif

