/**
 * \author Sylvain Marsat, University of Maryland - NASA GSFC
 *
 * \brief C code for the conversion between time scales.
 *
 * Formulas are taken from the UNITED STATES NAVAL OBSERVATORY CIRCULAR NO. 179 by George H. Kaplan
 * See http://aa.usno.navy.mil/publications/docs/Circular_179.php
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
#include "timeconversion.h"


/***********************************************************/
/********* Core functions to compute overlaps **************/

/* Table of leap seconds - from LAL XLALLeapSeconds.h */
static const struct leaps_table { double jd; int gpssec; int taiutc; } leaps[] =
{
  {2444239.5,    -43200, 19},  /* 1980-Jan-01 */
  {2444786.5,  46828800, 20},  /* 1981-Jul-01 */
  {2445151.5,  78364801, 21},  /* 1982-Jul-01 */
  {2445516.5, 109900802, 22},  /* 1983-Jul-01 */
  {2446247.5, 173059203, 23},  /* 1985-Jul-01 */
  {2447161.5, 252028804, 24},  /* 1988-Jan-01 */
  {2447892.5, 315187205, 25},  /* 1990-Jan-01 */
  {2448257.5, 346723206, 26},  /* 1991-Jan-01 */
  {2448804.5, 393984007, 27},  /* 1992-Jul-01 */
  {2449169.5, 425520008, 28},  /* 1993-Jul-01 */
  {2449534.5, 457056009, 29},  /* 1994-Jul-01 */
  {2450083.5, 504489610, 30},  /* 1996-Jan-01 */
  {2450630.5, 551750411, 31},  /* 1997-Jul-01 */
  {2451179.5, 599184012, 32},  /* 1999-Jan-01 */
  {2453736.5, 820108813, 33},  /* 2006-Jan-01 */
  {2454832.5, 914803214, 34},  /* 2009-Jan-01 */
  {2456109.5, 1025136015, 35}, /* 2012-Jul-01 */
  {2457204.5, 1119744016, 36}  /* 2015-Jul-01 */
};
static const int numleaps = sizeof( leaps ) / sizeof( *leaps );

/* Returns the leap seconds TAI-UTC at a given GPS second */
static int leapseconds(const double gpstime)
{
  int leap;

  if ( gpstime < leaps[0].gpssec )
  {
    printf( "Error - Don't know leap seconds before GPS time %d\n", leaps[0].gpssec );
    exit(1);
  }

  /* scan leap second table and locate the appropriate interval */
  for ( leap = 1; leap < numleaps; ++leap )
    if ( gpstime < leaps[leap].gpssec )
      break;

  return leaps[leap-1].taiutc;
}

/***********************************************************/
/********* Core functions to compute overlaps **************/

/* Function computing the correction from the earth rotation angle to the gmst angle */
/* Formula (2.12) of the USNO circular */
static double gmst_correction(const double gpstime) /* gpstime in seconds */
{
  /* T is the number of julian centuries since J2000.0 */
  double T = (gpstime - EPOCH_J2000_0_GPS)/(86400.0*36525.0);
  double correction = 2*PI*(0.014506 + 4612.156534*T + 1.3915817*T*T - 0.00000044*T*T*T - 0.000029956*T*T*T*T - 0.0000000368*T*T*T*T*T)/15.0/86400.0;
  return correction;
}

/* Earth rotation angle from the UT1 seconds elapsed since J2000.0 */
/* Formula (2.11) of the USNO circular */
static double era_angle_from_ut1(const double ut1_j2000) /* ut1_j2000 in seconds */
{
  double DU = ut1_j2000/86400.0;
  double fracDU = DU - trunc(DU);
  double era_angle = 2*PI*(0.7790572732640 + 0.00273781191135448*DU + fracDU);
  return era_angle;
}

/* Function computing the gmst angle from the gps time */
/* Formula (2.12) of the USNO circular */
double gmst_angle_from_gpstime(const double gpstime) /* gpstime in seconds */
{
  /* UTC time - UTC at J2000, computed from the GPS time */
  double utc_j2000 = (gpstime - EPOCH_J2000_0_GPS) - (leapseconds(gpstime) - EPOCH_J2000_0_TAI_UTC);

  /* Earth rotation angle */
  /* BEWARE: here we assimilate UTC and UT1 times - by construction, |UTC-UT1| < 1sec */
  double era_angle = era_angle_from_ut1(utc_j2000);

  /* GMST angle */
  double gmst_angle = era_angle + gmst_correction(gpstime);
  return gmst_angle;
}




