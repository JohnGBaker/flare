// WIP.H
//Written by John G Baker at gsfc.nasa.gov
//Weighted Inner-Product
#ifndef WIP_H
#define WIP_H

double complex wip_phase (double *f1, int n1, double *f2, int n2, double *s1Ar, double *s1Ai, double  *s1p, double *s2Ar, double*s2Ai, double *s2p, double (*Snoise)(double), double scalefactor, double min_f, double max_f);
double complex wip_adaptive_phase (double *f1, int n1, double *f2, int n2, double *s1Ar, double *s1Ai, double  *s1p, double *s2Ar, double*s2Ai, double *s2p, double (*Snoise)(double), double scalefactor, int downsample,  double errtol, double min_f, double max_f);

#endif
