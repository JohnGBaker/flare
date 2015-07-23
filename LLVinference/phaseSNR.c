#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <stdbool.h>

#include "LLVutils.h"

int main(int argc, char *argv[])
{
  /* Initialize structs for holding various options */
  LLVRunParams runParams;
  injectedparams = (LLVParams*) malloc(sizeof(LLVParams));
  memset(injectedparams, 0, sizeof(LLVParams));
  priorParams = (LLVPrior*) malloc(sizeof(LLVPrior));
  memset(priorParams, 0, sizeof(LLVPrior));
  
  /* Parse commandline to read parameters of injection */
  parse_args_LLV(argc, argv, injectedparams, priorParams, &runParams);

  /* Load and initialize the detector noise */
  LLVSimFD_Noise_Init_ParsePath();

  /* Initialize the data structure for the injection */
  LLVSignal* injectedsignal = NULL;
  LLVSignal_Init(&injectedsignal);

  /* Generate the injections and print SNRs */
  printf("# phiRef LHO LLO Virgo Network\n");
  for (double phase = 0.0; phase <= 2.0*M_PI; phase += 2.0*M_PI/100.0)
  {
  	injectedparams->phiRef = phase;
  	LLVGenerateSignal(injectedparams, injectedsignal);
  	printf("%g %g %g %g %g\n",phase,sqrt(injectedsignal->LHOhh),sqrt(injectedsignal->LLOhh),
  		sqrt(injectedsignal->VIRGOhh),sqrt(injectedsignal->LHOhh + injectedsignal->LLOhh + injectedsignal->VIRGOhh));
  }
  
  free(injectedparams);
  free(priorParams);

  return 0;
}
