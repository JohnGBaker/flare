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
  const detectortag tag)               /* Tag identifying the detector */
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
  const detectortag tag)               /* Tag identifying the detector */
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

/* Core function processing a signal (in the form of a list of modes) through the Fourier-domain LLV response (for a given detector), for given values of the inclination, position in the sky and polarization angle */
int LLVSimFDResponse(
  struct tagListmodesCAmpPhaseFrequencySeries **listhlm,  /* Input: list of modes in Frequency-domain amplitude and phase form as produced by the ROM */
  struct tagListmodesCAmpPhaseFrequencySeries **lists,    /* Output: list of contribution of each mode in the detector signal, in Frequency-domain amplitude and phase form, for the given detector and sky position */
  const double gpstime,                                   /* GPS time (s) when the signal at coalescence reaches geocenter */
  const double ra,                                        /* Position in the sky: J2000.0 right ascension (rad) */
  const double dec,                                       /* Position in the sky: J2000.0 declination (rad) */
  const double inclination,                               /* Inclination of the source (rad) */
  const double psi,                                       /* Polarization angle (rad) */
  const detectortag tag)                                  /* Tag identifying the detector */
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
