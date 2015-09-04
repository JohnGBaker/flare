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
  globalparams = (LLVGlobalParams*) malloc(sizeof(LLVGlobalParams));
  memset(globalparams, 0, sizeof(LLVGlobalParams));
  
  /* Parse commandline to read parameters of injection */
  parse_args_LLV(argc, argv, injectedparams, globalparams, priorParams, &runParams);
  injectedparams->nbmode = globalparams->nbmodeinj;

  /* Load and initialize the detector noise */
  LLVSimFD_Noise_Init_ParsePath();

  /* Initialize the data structure for the injection */
  LLVSignal* injectedsignal = NULL;
  LLVSignal_Init(&injectedsignal);

  double Mlist[7] = {10., 20., 30., 50., 100., 200., 400.};
  double Qlist[4] = {1., 2., 4., 8.};
  double Mlist2[2] = {20., 200.};
  double Qlist2[2] = {1., 4.};

  double m1, m2;
  int i, j, k=-1;

  injectedparams->distance = 100.;

  // print command line
  printf("# ");
  for (i=0; i<argc; i++) printf("%s ", argv[i]);
  printf("\n");

  for (i=0; i<7; i++)
  {
    for (j=0; j<4; j++)
    {
      k++;
      m2 = Mlist[i] / (1. + Qlist[j]);
      m1 = Mlist[i] - m2;
      injectedparams->m1 = m1;
      injectedparams->m2 = m2;
      LLVGenerateSignal(injectedparams, injectedsignal);
      printf("%d %f %f %f %f %f\n", k, Mlist[i], Qlist[j], m1, m2, 
        100.*sqrt(injectedsignal->LHOhh + injectedsignal->LLOhh + injectedsignal->VIRGOhh)/12.);
    }
  }

  for (i=0; i<2; i++)
  {
    for (j=0; j<2; j++)
    {
      k++;
      m2 = Mlist2[i] / (1. + Qlist2[j]);
      m1 = Mlist2[i] - m2;
      injectedparams->m1 = m1;
      injectedparams->m2 = m2;
      LLVGenerateSignal(injectedparams, injectedsignal);
      printf("%d %f %f %f %f %f\n", k, Mlist2[i], Qlist2[j], m1, m2, 
        100.*sqrt(injectedsignal->LHOhh + injectedsignal->LLOhh + injectedsignal->VIRGOhh)/15.);
    }
  }

  for (i=0; i<2; i++)
  {
    for (j=0; j<2; j++)
    {
      k++;
      m2 = Mlist2[i] / (1. + Qlist2[j]);
      m1 = Mlist2[i] - m2;
      injectedparams->m1 = m1;
      injectedparams->m2 = m2;
      LLVGenerateSignal(injectedparams, injectedsignal);
      printf("%d %f %f %f %f %f\n", k, Mlist2[i], Qlist2[j], m1, m2, 
        100.*sqrt(injectedsignal->LHOhh + injectedsignal->LLOhh + injectedsignal->VIRGOhh)/20.);
    }
  }
  
  free(injectedparams);
  free(priorParams);

  return 0;
}
