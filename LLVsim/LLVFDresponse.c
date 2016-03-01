/**
 * \author Sylvain Marsat, University of Maryland - NASA GSFC
 *
 * \brief C code for the implementation of the Fourier domain response for LIGO-VIRGO detectors.
 *
 */


#define _XOPEN_SOURCE 500

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <unistd.h>
#include <getopt.h>
#include <stdbool.h>
#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_complex.h>

#include "constants.h"
#include "LLVgeometry.h"
#include "struct.h"
#include "LLVFDresponse.h"
#include "timeconversion.h"


/************************************************************/
/********* Functions setting detector geometry **************/

/* Function setting the response matrix of a given detector, in cartesian coordinates */
void SetMatrixD(
  gsl_matrix* D,                       /* Output: matrix of the detector response Dij */
  const Detectortag tag)               /* Tag identifying the detector */
{
  /* Allocating and defining the cartesian unit vectors along the two arms */
  gsl_vector* nx = gsl_vector_alloc(3);
  gsl_vector* ny = gsl_vector_alloc(3);
  if(tag==LHO) {
    gsl_vector_set(nx, 0, LAL_LHO_4K_ARM_X_DIRECTION_X);
    gsl_vector_set(nx, 1, LAL_LHO_4K_ARM_X_DIRECTION_Y);
    gsl_vector_set(nx, 2, LAL_LHO_4K_ARM_X_DIRECTION_Z);
    gsl_vector_set(ny, 0, LAL_LHO_4K_ARM_Y_DIRECTION_X);
    gsl_vector_set(ny, 1, LAL_LHO_4K_ARM_Y_DIRECTION_Y);
    gsl_vector_set(ny, 2, LAL_LHO_4K_ARM_Y_DIRECTION_Z);
  }
  else if(tag==LLO) {
    gsl_vector_set(nx, 0, LAL_LLO_4K_ARM_X_DIRECTION_X);
    gsl_vector_set(nx, 1, LAL_LLO_4K_ARM_X_DIRECTION_Y);
    gsl_vector_set(nx, 2, LAL_LLO_4K_ARM_X_DIRECTION_Z);
    gsl_vector_set(ny, 0, LAL_LLO_4K_ARM_Y_DIRECTION_X);
    gsl_vector_set(ny, 1, LAL_LLO_4K_ARM_Y_DIRECTION_Y);
    gsl_vector_set(ny, 2, LAL_LLO_4K_ARM_Y_DIRECTION_Z);
  }
  else if(tag==VIRGO) {
    gsl_vector_set(nx, 0, LAL_VIRGO_ARM_X_DIRECTION_X);
    gsl_vector_set(nx, 1, LAL_VIRGO_ARM_X_DIRECTION_Y);
    gsl_vector_set(nx, 2, LAL_VIRGO_ARM_X_DIRECTION_Z);
    gsl_vector_set(ny, 0, LAL_VIRGO_ARM_Y_DIRECTION_X);
    gsl_vector_set(ny, 1, LAL_VIRGO_ARM_Y_DIRECTION_Y);
    gsl_vector_set(ny, 2, LAL_VIRGO_ARM_Y_DIRECTION_Z);
  }
  else {
    printf("Error: invalid detector tag\n");
    exit(1);
  }

  /* Setting the components of the matrix D = 1/2 (nx nx - ny ny) */
  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++) {
      gsl_matrix_set(D, i, j, 1./2*(gsl_vector_get(nx, i)*gsl_vector_get(nx, j) - gsl_vector_get(ny, i)*gsl_vector_get(ny, j)));
    }
  }
  gsl_vector_free(nx);
  gsl_vector_free(ny);
}

/* Function setting the position of a detector, in cartesian coordinates */
void SetVectorXd(
  gsl_vector* Xd,                      /* Output: position vector of the detector */
  const Detectortag tag)               /* Tag identifying the detector */
{
  if(tag==LHO) {
    gsl_vector_set(Xd, 0, LAL_LHO_4K_VERTEX_LOCATION_X_SI);
    gsl_vector_set(Xd, 1, LAL_LHO_4K_VERTEX_LOCATION_Y_SI);
    gsl_vector_set(Xd, 2, LAL_LHO_4K_VERTEX_LOCATION_Z_SI);
  }
  else if(tag==LLO) {
    gsl_vector_set(Xd, 0, LAL_LLO_4K_VERTEX_LOCATION_X_SI);
    gsl_vector_set(Xd, 1, LAL_LLO_4K_VERTEX_LOCATION_Y_SI);
    gsl_vector_set(Xd, 2, LAL_LLO_4K_VERTEX_LOCATION_Z_SI);
  }
  else if(tag==VIRGO) {
    gsl_vector_set(Xd, 0, LAL_VIRGO_VERTEX_LOCATION_X_SI);
    gsl_vector_set(Xd, 1, LAL_VIRGO_VERTEX_LOCATION_Y_SI);
    gsl_vector_set(Xd, 2, LAL_VIRGO_VERTEX_LOCATION_Z_SI);
  }
  else {
    printf("Error: invalid detector tag\n");
    exit(1);
  }
}

/* Function setting the cartesian coordinates of the wave frame vectors (X,Y,Z), given the position in the sky and polarization */
void SetVectorsXYZ(
  gsl_vector* X,                       /* Output: cartesian vector of the wave frame unit vector X */
  gsl_vector* Y,                       /* Output: cartesian vector of the wave frame unit vector Y */
  gsl_vector* Z,                       /* Output: cartesian vector of the wave frame unit vector Z */
  const double theta,                  /* First angle for the position in the sky (Earth-based spherical angle) */
  const double phi,                    /* Second angle for the position in the sky (Earth-based spherical angle) */
  const double psi)                    /* Polarization angle */
{
  double cospsi = cos(psi);
  double sinpsi = sin(psi);
  double costheta = cos(theta);
  double sintheta = sin(theta);
  double cosphi = cos(phi);
  double sinphi = sin(phi);
  /* Unit vector X */
  gsl_vector_set(X, 0, -cospsi*sinphi+sinpsi*costheta*cosphi);
  gsl_vector_set(X, 1, cospsi*cosphi+sinpsi*costheta*sinphi);
  gsl_vector_set(X, 2, -sinpsi*sintheta);
  /* Unit vector Y */
  gsl_vector_set(Y, 0, sinpsi*sinphi+cospsi*costheta*cosphi);
  gsl_vector_set(Y, 1, -sinpsi*cosphi+cospsi*costheta*sinphi);
  gsl_vector_set(Y, 2, -cospsi*sintheta);
  /* Unit vector Z */
  gsl_vector_set(Z, 0, -sintheta*cosphi);
  gsl_vector_set(Z, 1, -sintheta*sinphi);
  gsl_vector_set(Z, 2, -costheta);
}

/***************************************/
/********* Core functions **************/

/* Function to convert string input network string to Networktag */
Networktag ParseNetworktag(char* string) {
  Networktag tag;
  if(strcmp(string, "L")==0) tag = L;
  else if(strcmp(string, "H")==0) tag = H;
  else if(strcmp(string, "V")==0) tag = V;
  else if(strcmp(string, "LH")==0) tag = LH;
  else if(strcmp(string, "LV")==0) tag = LV;
  else if(strcmp(string, "HV")==0) tag = HV;
  else if(strcmp(string, "LHV")==0) tag = LHV;
  else {
    printf("Error in ParseNetworktag: string not recognized.\n");
    exit(1);
  }
  return tag;
}

/* Function evaluating the Fourier-domain factors for LHV detectors */
/* Note: in case only one or two detectors are considered, amplitudes for the missing ones are set to 0 */
/* Contrarily to the LISA case, here returns trivially 0 or 1 */
/* (allows minimal changes from the old structure that assumed LHV - but not optimal) */
static int EvaluateDetectorFactor3Det(
  double* factor1,                       /* Output for factor for detector 1 */
  double* factor2,                       /* Output for factor for detector 2 */
  double* factor3,                       /* Output for factor for detector 3 */
  const Networktag networktag)                  /* Selector for the detector network */
{
  switch(networktag) {
    case L:
    *factor1 = 1.;
    *factor2 = 0;
    *factor3 = 0;
    break;
    case H:
    *factor1 = 0;
    *factor2 = 1.;
    *factor3 = 0;
    break;
    case V:
    *factor1 = 0;
    *factor2 = 0;
    *factor3 = 1.;
    break;
    case LH:
    *factor1 = 1.;
    *factor2 = 1.;
    *factor3 = 0;
    break;
    case LV:
    *factor1 = 1.;
    *factor2 = 0;
    *factor3 = 1.;
    break;
    case HV:
    *factor1 = 0;
    *factor2 = 1.;
    *factor3 = 1.;
    break;
    case LHV:
    *factor1 = 1.;
    *factor2 = 1.;
    *factor3 = 1.;
    break;
    default:
    printf("Error in EvaluateDetectorFactor3Det: networktag not recognized.\n");
    exit(1);
  }
  return SUCCESS;
}

/* Core function processing a signal (in the form of a list of modes) through the Fourier-domain LLV response (for a given detector), for given values of the inclination, position in the sky and polarization angle */
int LLVSimFDResponse(
  struct tagListmodesCAmpPhaseFrequencySeries **listhlm,  /* Input: list of modes in Frequency-domain amplitude and phase form as produced by the ROM */
  struct tagListmodesCAmpPhaseFrequencySeries **lists,    /* Output: list of contribution of each mode in the detector signal, in Frequency-domain amplitude and phase form, for the given detector and sky position */
  const double gpstime,                                   /* GPS time (s) when the signal at coalescence reaches geocenter */
  const double ra,                                        /* Position in the sky: J2000.0 right ascension (rad) */
  const double dec,                                       /* Position in the sky: J2000.0 declination (rad) */
  const double inclination,                               /* Inclination of the source (rad) */
  const double psi,                                       /* Polarization angle (rad) */
  const Detectortag tag)                                  /* Tag identifying the detector */
{
  /* Define matrix D and position vector Xd of the detector */
  gsl_matrix* D = gsl_matrix_alloc(3,3);
  gsl_vector* Xd = gsl_vector_alloc(3);
  SetVectorXd(Xd, tag);
  SetMatrixD(D, tag);

  /* Conversion from (ra, dec) to the Earth-based spherical angles (theta, phi) - neglecting nutation and precession, and identifying UT1 and UTC, so accurate roughly to a second of time */
  double gmst_angle = gmst_angle_from_gpstime(gpstime);
  double theta = PI/2 - dec;
  double phi = ra - gmst_angle;

  /* Define waveframe unit vectors (X,Y,Z) */
  gsl_vector* X = gsl_vector_alloc(3);
  gsl_vector* Y = gsl_vector_alloc(3);
  gsl_vector* Z = gsl_vector_alloc(3);
  SetVectorsXYZ(X, Y, Z, theta, phi, psi);

  /* Compute the delay from geocenter to the detector */
  double delaylength;
  gsl_blas_ddot(Xd, Z, &delaylength);
  double twopidelay = 2*PI*delaylength/C_SI;
  //printf("Delay: %g\n", twopidelay/(2*PI));

  /* Compute the value of pattern functions Fplus, Fcross */
  gsl_vector* DX = gsl_vector_calloc(3); /* Temporary vector D.X, initialized to 0 */
  gsl_vector* DY = gsl_vector_calloc(3); /* Temporary vector D.Y, initialized to 0 */
  gsl_blas_dgemv( CblasNoTrans, 1., D, X, 0, DX);
  gsl_blas_dgemv( CblasNoTrans, 1., D, Y, 0, DY);
  double XDX;
  double XDY;
  double YDY;
  gsl_blas_ddot(X, DX, &XDX);
  gsl_blas_ddot(X, DY, &XDY);
  gsl_blas_ddot(Y, DY, &YDY);
  double Fplus = XDX - YDY;
  double Fcross = 2*XDY;
  //printf("Fplus: %g\n", Fplus);
  //printf("Fcross: %g\n", Fcross);

  /* Main loop over the modes - goes through all the modes present, stopping when encountering NULL */
  ListmodesCAmpPhaseFrequencySeries* listelement = *listhlm;
  while(listelement) {

    /* Definitions: l,m, frequency series and length */
    int l = listelement->l;
    int m = listelement->m;
    CAmpPhaseFrequencySeries* freqseries = listelement->freqseries;
    gsl_vector* freq = freqseries->freq;
    gsl_vector* amp_real = freqseries->amp_real;
    gsl_vector* amp_imag = freqseries->amp_imag;
    gsl_vector* phase = freqseries->phase;
    int len = (int) freq->size;
    double f;
    double complex camp;
    double complex camps;

    /* Computing the Ylm combined factors for plus and cross for this mode */
    /* Capital Phi is set to 0 by convention */
    double complex Yfactorplus;
    double complex Yfactorcross;
    if (!(l%2)) {
      Yfactorplus = 1./2 * (SpinWeightedSphericalHarmonic(inclination, 0., -2, l, m) + conj(SpinWeightedSphericalHarmonic(inclination, 0., -2, l, -m)));
      Yfactorcross = I/2 * (SpinWeightedSphericalHarmonic(inclination, 0., -2, l, m) - conj(SpinWeightedSphericalHarmonic(inclination, 0., -2, l, -m)));
    }
    else {
      Yfactorplus = 1./2 * (SpinWeightedSphericalHarmonic(inclination, 0., -2, l, m) - conj(SpinWeightedSphericalHarmonic(inclination, 0., -2, l, -m)));
      Yfactorcross = I/2 * (SpinWeightedSphericalHarmonic(inclination, 0., -2, l, m) + conj(SpinWeightedSphericalHarmonic(inclination, 0., -2, l, -m)));
    }
    //if(l==2 && m==2){
    //  if (!(l%2)) {
    //printf("even\n");
    //}
    //else {printf("odd\n");}
    //printf("SpinWeightedSphericalHarmonic(inclination, 0., -2, l, m): %g+I*%g\n", creal(SpinWeightedSphericalHarmonic(inclination, 0., -2, l, m)), cimag(SpinWeightedSphericalHarmonic(inclination, 0., -2, l, m)));
    //printf("SpinWeightedSphericalHarmonic(inclination, 0., -2, l, -m): %g+I*%g\n", creal(SpinWeightedSphericalHarmonic(inclination, 0., -2, l, -m)), cimag(SpinWeightedSphericalHarmonic(inclination, 0., -2, l, -m)));
    //printf("Yfactorplus: %g + I*%g\n", creal(Yfactorplus), cimag(Yfactorplus));
    //printf("Yfactorcross: %g + I*%g\n", creal(Yfactorcross), cimag(Yfactorcross));
    //}

    /* Initializing frequency series structure for this mode, for the signal s = F+ h+ + Fx hx */
    CAmpPhaseFrequencySeries *modefreqseriess = NULL;
    CAmpPhaseFrequencySeries_Init(&modefreqseriess, len);
    gsl_vector* freqs = modefreqseriess->freq;
    gsl_vector* amp_reals = modefreqseriess->amp_real;
    gsl_vector* amp_imags = modefreqseriess->amp_imag;
    gsl_vector* phases = modefreqseriess->phase;

    /* Loop over the frequencies - multiplying by F+, Fx, and add phase due to the delay from geocenter to the detector */
    double complex factorcamp = Fplus*Yfactorplus + Fcross*Yfactorcross;
    for(int j=0; j<len; j++) {
      f = gsl_vector_get(freq, j);
      camp = gsl_vector_get(amp_real, j) + I*gsl_vector_get(amp_imag, j);
      camps = camp * factorcamp;
      gsl_vector_set(amp_reals, j, creal(camps));
      gsl_vector_set(amp_imags, j, cimag(camps));
      gsl_vector_set(phases, j, gsl_vector_get(phase, j) + twopidelay*f);
    }
    /* Copying the vectors of frequencies */
    gsl_vector_memcpy(freqs, freq);

    /* Append the modes to the ouput list-of-modes structures */
    *lists = ListmodesCAmpPhaseFrequencySeries_AddModeNoCopy(*lists, modefreqseriess, l, m);

    /* Going to the next mode in the list */
    listelement = listelement->next;
  }

  /* Cleaning */
  gsl_matrix_free(D);
  gsl_vector_free(Xd);
  gsl_vector_free(X);
  gsl_vector_free(Y);
  gsl_vector_free(Z);
  gsl_vector_free(DX);
  gsl_vector_free(DY);

  return SUCCESS;
}

/* Core function processing a signal (in the form of a list of modes) through the Fourier-domain LLV response for a given detector network, for given values of the inclination, position in the sky and polarization angle */
/* Note: as for now, asssumes the three detectors are L,H,V - amplitudes simply set to 0 in those that are not selected by networktag */
int LLVSimFDResponse3Det(
  struct tagListmodesCAmpPhaseFrequencySeries **listhlm,  /* Input: list of modes in Frequency-domain amplitude and phase form as produced by the ROM */
  struct tagListmodesCAmpPhaseFrequencySeries **list1,    /* Output: list of contributions of each mode in the signal of detector 1, in Frequency-domain amplitude and phase form */
  struct tagListmodesCAmpPhaseFrequencySeries **list2,    /* Output: list of contributions of each mode in the signal of detector 1, in Frequency-domain amplitude and phase form */
  struct tagListmodesCAmpPhaseFrequencySeries **list3,    /* Output: list of contributions of each mode in the signal of detector 1, in Frequency-domain amplitude and phase form */
  const double gpstime,                                   /* GPS time (s) when the signal at coalescence reaches geocenter */
  const double ra,                                        /* Position in the sky: J2000.0 right ascension (rad) */
  const double dec,                                       /* Position in the sky: J2000.0 declination (rad) */
  const double inclination,                               /* Inclination of the source (rad) */
  const double psi,                                       /* Polarization angle (rad) */
  const Networktag tag)                                  /* Tag identifying the detector */
{
  /* Read which detectors are to be included in the network */
  double factor1 = 0;
  double factor2 = 0;
  double factor3 = 0;
  EvaluateDetectorFactor3Det(&factor1, &factor2, &factor3, tag);

  /* Define matrix D and position vector Xd for the 3 detectors */
  gsl_matrix* D1 = gsl_matrix_alloc(3,3);
  gsl_matrix* D2 = gsl_matrix_alloc(3,3);
  gsl_matrix* D3 = gsl_matrix_alloc(3,3);
  gsl_vector* Xd1 = gsl_vector_alloc(3);
  gsl_vector* Xd2 = gsl_vector_alloc(3);
  gsl_vector* Xd3 = gsl_vector_alloc(3);
  SetVectorXd(Xd1, LHO);
  SetVectorXd(Xd2, LLO);
  SetVectorXd(Xd3, VIRGO);
  SetMatrixD(D1, LHO);
  SetMatrixD(D2, LLO);
  SetMatrixD(D3, VIRGO);

  /* Conversion from (ra, dec) to the Earth-based spherical angles (theta, phi) - neglecting nutation and precession, and identifying UT1 and UTC, so accurate roughly to a second of time */
  double gmst_angle = gmst_angle_from_gpstime(gpstime);
  double theta = PI/2 - dec;
  double phi = ra - gmst_angle;

  /* Define waveframe unit vectors (X,Y,Z) */
  gsl_vector* X = gsl_vector_alloc(3);
  gsl_vector* Y = gsl_vector_alloc(3);
  gsl_vector* Z = gsl_vector_alloc(3);
  SetVectorsXYZ(X, Y, Z, theta, phi, psi);

  /* Compute the delays from geocenter to each detector */
  double delaylength1 = 0;
  double delaylength2 = 0;
  double delaylength3 = 0;
  gsl_blas_ddot(Xd1, Z, &delaylength1);
  gsl_blas_ddot(Xd2, Z, &delaylength2);
  gsl_blas_ddot(Xd3, Z, &delaylength3);
  double twopidelay1 = 2*PI*delaylength1/C_SI;
  double twopidelay2 = 2*PI*delaylength2/C_SI;
  double twopidelay3 = 2*PI*delaylength3/C_SI;

  /* Compute the value of pattern functions Fplus, Fcross for each detector */
  gsl_vector* DX1 = gsl_vector_calloc(3); /* Temporary vector D.X, initialized to 0 */
  gsl_vector* DY1 = gsl_vector_calloc(3); /* Temporary vector D.Y, initialized to 0 */
  gsl_vector* DX2 = gsl_vector_calloc(3); /* Temporary vector D.X, initialized to 0 */
  gsl_vector* DY2 = gsl_vector_calloc(3); /* Temporary vector D.Y, initialized to 0 */
  gsl_vector* DX3 = gsl_vector_calloc(3); /* Temporary vector D.X, initialized to 0 */
  gsl_vector* DY3 = gsl_vector_calloc(3); /* Temporary vector D.Y, initialized to 0 */
  gsl_blas_dgemv( CblasNoTrans, 1., D1, X, 0, DX1);
  gsl_blas_dgemv( CblasNoTrans, 1., D1, Y, 0, DY1);
  gsl_blas_dgemv( CblasNoTrans, 1., D2, X, 0, DX2);
  gsl_blas_dgemv( CblasNoTrans, 1., D2, Y, 0, DY2);
  gsl_blas_dgemv( CblasNoTrans, 1., D3, X, 0, DX3);
  gsl_blas_dgemv( CblasNoTrans, 1., D3, Y, 0, DY3);
  double XDX1 = 0;
  double XDY1 = 0;
  double YDY1 = 0;
  double XDX2 = 0;
  double XDY2 = 0;
  double YDY2 = 0;
  double XDX3 = 0;
  double XDY3 = 0;
  double YDY3 = 0;
  gsl_blas_ddot(X, DX1, &XDX1);
  gsl_blas_ddot(X, DY1, &XDY1);
  gsl_blas_ddot(X, DX2, &XDX2);
  gsl_blas_ddot(X, DY2, &XDY2);
  gsl_blas_ddot(X, DX3, &XDX3);
  gsl_blas_ddot(X, DY3, &XDY3);
  double Fplus1 = XDX1 - YDY1;
  double Fcross1 = 2*XDY1;
  double Fplus2 = XDX2 - YDY2;
  double Fcross2 = 2*XDY2;
  double Fplus3 = XDX3 - YDY3;
  double Fcross3 = 2*XDY3;

  /* Main loop over the modes - goes through all the modes present, stopping when encountering NULL */
  ListmodesCAmpPhaseFrequencySeries* listelement = *listhlm;
  while(listelement) {

    /* Definitions: l,m, frequency series and length */
    int l = listelement->l;
    int m = listelement->m;
    CAmpPhaseFrequencySeries* freqseries = listelement->freqseries;
    gsl_vector* freq = freqseries->freq;
    gsl_vector* amp_real = freqseries->amp_real;
    gsl_vector* amp_imag = freqseries->amp_imag;
    gsl_vector* phase = freqseries->phase;
    int len = (int) freq->size;
    double f;
    double complex camp;
    double complex camp1;
    double complex camp2;
    double complex camp3;

    /* Computing the Ylm combined factors for plus and cross for this mode */
    /* Capital Phi is set to 0 by convention */
    double complex Yfactorplus;
    double complex Yfactorcross;
    if (!(l%2)) {
      Yfactorplus = 1./2 * (SpinWeightedSphericalHarmonic(inclination, 0., -2, l, m) + conj(SpinWeightedSphericalHarmonic(inclination, 0., -2, l, -m)));
      Yfactorcross = I/2 * (SpinWeightedSphericalHarmonic(inclination, 0., -2, l, m) - conj(SpinWeightedSphericalHarmonic(inclination, 0., -2, l, -m)));
    }
    else {
      Yfactorplus = 1./2 * (SpinWeightedSphericalHarmonic(inclination, 0., -2, l, m) - conj(SpinWeightedSphericalHarmonic(inclination, 0., -2, l, -m)));
      Yfactorcross = I/2 * (SpinWeightedSphericalHarmonic(inclination, 0., -2, l, m) + conj(SpinWeightedSphericalHarmonic(inclination, 0., -2, l, -m)));
    }

    /* Initializing frequency series structure for this mode, for the signal s = F+ h+ + Fx hx in each detector */
    CAmpPhaseFrequencySeries *modefreqseries1 = NULL;
    CAmpPhaseFrequencySeries *modefreqseries2 = NULL;
    CAmpPhaseFrequencySeries *modefreqseries3 = NULL;
    CAmpPhaseFrequencySeries_Init(&modefreqseries1, len);
    CAmpPhaseFrequencySeries_Init(&modefreqseries2, len);
    CAmpPhaseFrequencySeries_Init(&modefreqseries3, len);
    gsl_vector* freq1 = modefreqseries1->freq;
    gsl_vector* freq2 = modefreqseries2->freq;
    gsl_vector* freq3 = modefreqseries3->freq;
    gsl_vector* amp_real1 = modefreqseries1->amp_real;
    gsl_vector* amp_real2 = modefreqseries2->amp_real;
    gsl_vector* amp_real3 = modefreqseries3->amp_real;
    gsl_vector* amp_imag1 = modefreqseries1->amp_imag;
    gsl_vector* amp_imag2 = modefreqseries2->amp_imag;
    gsl_vector* amp_imag3 = modefreqseries3->amp_imag;
    gsl_vector* phase1 = modefreqseries1->phase;
    gsl_vector* phase2 = modefreqseries2->phase;
    gsl_vector* phase3 = modefreqseries3->phase;

    /* Loop over the frequencies - multiplying by F+, Fx, and add phase due to the delay from geocenter to the detector */
    double complex factorcamp1 = Fplus1*Yfactorplus + Fcross1*Yfactorcross;
    double complex factorcamp2 = Fplus2*Yfactorplus + Fcross2*Yfactorcross;
    double complex factorcamp3 = Fplus3*Yfactorplus + Fcross3*Yfactorcross;
    for(int j=0; j<len; j++) {
      f = gsl_vector_get(freq, j);
      camp = gsl_vector_get(amp_real, j) + I*gsl_vector_get(amp_imag, j);
      /* Note: include factor123 as detector selectors, which are simpy 0 or 1 whether the detector is in the array or not */
      camp1 = factor1 * camp * factorcamp1;
      camp2 = factor2 * camp * factorcamp2;
      camp3 = factor3 * camp * factorcamp3;
      gsl_vector_set(amp_real1, j, creal(camp1));
      gsl_vector_set(amp_real2, j, creal(camp2));
      gsl_vector_set(amp_real3, j, creal(camp3));
      gsl_vector_set(amp_imag1, j, cimag(camp1));
      gsl_vector_set(amp_imag2, j, cimag(camp2));
      gsl_vector_set(amp_imag3, j, cimag(camp3));
      gsl_vector_set(phase1, j, gsl_vector_get(phase, j) + twopidelay1*f);
      gsl_vector_set(phase2, j, gsl_vector_get(phase, j) + twopidelay2*f);
      gsl_vector_set(phase3, j, gsl_vector_get(phase, j) + twopidelay3*f);
    }
    /* Copying the vectors of frequencies */
    gsl_vector_memcpy(freq1, freq);
    gsl_vector_memcpy(freq2, freq);
    gsl_vector_memcpy(freq3, freq);

    /* Append the modes to the ouput list-of-modes structures */
    *list1 = ListmodesCAmpPhaseFrequencySeries_AddModeNoCopy(*list1, modefreqseries1, l, m);
    *list2 = ListmodesCAmpPhaseFrequencySeries_AddModeNoCopy(*list2, modefreqseries2, l, m);
    *list3 = ListmodesCAmpPhaseFrequencySeries_AddModeNoCopy(*list3, modefreqseries3, l, m);

    /* Going to the next mode in the list */
    listelement = listelement->next;
  }

  /* Cleaning */
  gsl_matrix_free(D1);
  gsl_matrix_free(D2);
  gsl_matrix_free(D3);
  gsl_vector_free(Xd1);
  gsl_vector_free(Xd2);
  gsl_vector_free(Xd3);
  gsl_vector_free(X);
  gsl_vector_free(Y);
  gsl_vector_free(Z);
  gsl_vector_free(DX1);
  gsl_vector_free(DX2);
  gsl_vector_free(DX3);
  gsl_vector_free(DY1);
  gsl_vector_free(DY2);
  gsl_vector_free(DY3);

  return SUCCESS;
}
