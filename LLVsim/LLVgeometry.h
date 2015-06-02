/**
 * \author Sylvain Marsat, University of Maryland - NASA GSFC
 *
 * \brief C header for the geometric coefficients entering the response for LIGO-VIRGO detectors.
 *
 *
 */

#ifndef _LLVGEOMETRY_H
#define _LLVGEOMETRY_H

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

/********************************************************/
/**************** Type definitions **********************/

/* Enumerator for the detector selector */
typedef enum {LHO, LLO, VIRGO} detectortag;

/*************************************************************/
/**************** Geometrical constants **********************/

/**
 * \name LIGO Hanford Observatory 4km Interferometric Detector constants
 * The following constants describe the location and geometry of the
 * LIGO Hanford Observatory 4km Interferometric Detector.
 */
#define LAL_LHO_4K_DETECTOR_NAME               	"LHO_4k"	/**< LHO_4k detector name string */
#define LAL_LHO_4K_DETECTOR_PREFIX             	"H1"	/**< LHO_4k detector prefix string */
#define LAL_LHO_4K_DETECTOR_LONGITUDE_RAD      	-2.08405676917	/**< LHO_4k vertex longitude (rad) */
#define LAL_LHO_4K_DETECTOR_LATITUDE_RAD       	0.81079526383	/**< LHO_4k vertex latitude (rad) */
#define LAL_LHO_4K_DETECTOR_ELEVATION_SI       	142.554	/**< LHO_4k vertex elevation (m) */
#define LAL_LHO_4K_DETECTOR_ARM_X_AZIMUTH_RAD  	5.65487724844	/**< LHO_4k x arm azimuth (rad) */
#define LAL_LHO_4K_DETECTOR_ARM_Y_AZIMUTH_RAD  	4.08408092164	/**< LHO_4k y arm azimuth (rad) */
#define LAL_LHO_4K_DETECTOR_ARM_X_ALTITUDE_RAD 	-0.00061950000	/**< LHO_4k x arm altitude (rad) */
#define LAL_LHO_4K_DETECTOR_ARM_Y_ALTITUDE_RAD 	0.00001250000	/**< LHO_4k y arm altitude (rad) */
#define LAL_LHO_4K_DETECTOR_ARM_X_MIDPOINT_SI  	1997.54200000000	/**< LHO_4k x arm midpoint (m) */
#define LAL_LHO_4K_DETECTOR_ARM_Y_MIDPOINT_SI  	1997.52200000000	/**< LHO_4k y arm midpoint (m) */
#define LAL_LHO_4K_VERTEX_LOCATION_X_SI        	-2.16141492636e+06	/**< LHO_4k x-component of vertex location in Earth-centered frame (m) */
#define LAL_LHO_4K_VERTEX_LOCATION_Y_SI        	-3.83469517889e+06	/**< LHO_4k y-component of vertex location in Earth-centered frame (m) */
#define LAL_LHO_4K_VERTEX_LOCATION_Z_SI        	4.60035022664e+06	/**< LHO_4k z-component of vertex location in Earth-centered frame (m) */
#define LAL_LHO_4K_ARM_X_DIRECTION_X           	-0.22389266154	/**< LHO_4k x-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_LHO_4K_ARM_X_DIRECTION_Y           	0.79983062746	/**< LHO_4k y-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_LHO_4K_ARM_X_DIRECTION_Z           	0.55690487831	/**< LHO_4k z-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_LHO_4K_ARM_Y_DIRECTION_X           	-0.91397818574	/**< LHO_4k x-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_LHO_4K_ARM_Y_DIRECTION_Y           	0.02609403989	/**< LHO_4k y-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_LHO_4K_ARM_Y_DIRECTION_Z           	-0.40492342125	/**< LHO_4k z-component of unit vector pointing along y arm in Earth-centered frame */


/**
 * \name LIGO Livingston Observatory 4km Interferometric Detector constants
 * The following constants describe the location and geometry of the
 * LIGO Livingston Observatory 4km Interferometric Detector.
 */
#define LAL_LLO_4K_DETECTOR_NAME               	"LLO_4k"	/**< LLO_4k detector name string */
#define LAL_LLO_4K_DETECTOR_PREFIX             	"L1"	/**< LLO_4k detector prefix string */
#define LAL_LLO_4K_DETECTOR_LONGITUDE_RAD      	-1.58430937078	/**< LLO_4k vertex longitude (rad) */
#define LAL_LLO_4K_DETECTOR_LATITUDE_RAD       	0.53342313506	/**< LLO_4k vertex latitude (rad) */
#define LAL_LLO_4K_DETECTOR_ELEVATION_SI       	-6.574	/**< LLO_4k vertex elevation (m) */
#define LAL_LLO_4K_DETECTOR_ARM_X_AZIMUTH_RAD  	4.40317772346	/**< LLO_4k x arm azimuth (rad) */
#define LAL_LLO_4K_DETECTOR_ARM_Y_AZIMUTH_RAD  	2.83238139666	/**< LLO_4k y arm azimuth (rad) */
#define LAL_LLO_4K_DETECTOR_ARM_X_ALTITUDE_RAD 	-0.00031210000	/**< LLO_4k x arm altitude (rad) */
#define LAL_LLO_4K_DETECTOR_ARM_Y_ALTITUDE_RAD 	-0.00061070000	/**< LLO_4k y arm altitude (rad) */
#define LAL_LLO_4K_DETECTOR_ARM_X_MIDPOINT_SI  	1997.57500000000	/**< LLO_4k x arm midpoint (m) */
#define LAL_LLO_4K_DETECTOR_ARM_Y_MIDPOINT_SI  	1997.57500000000	/**< LLO_4k y arm midpoint (m) */
#define LAL_LLO_4K_VERTEX_LOCATION_X_SI        	-7.42760447238e+04	/**< LLO_4k x-component of vertex location in Earth-centered frame (m) */
#define LAL_LLO_4K_VERTEX_LOCATION_Y_SI        	-5.49628371971e+06	/**< LLO_4k y-component of vertex location in Earth-centered frame (m) */
#define LAL_LLO_4K_VERTEX_LOCATION_Z_SI        	3.22425701744e+06	/**< LLO_4k z-component of vertex location in Earth-centered frame (m) */
#define LAL_LLO_4K_ARM_X_DIRECTION_X           	-0.95457412153	/**< LLO_4k x-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_LLO_4K_ARM_X_DIRECTION_Y           	-0.14158077340	/**< LLO_4k y-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_LLO_4K_ARM_X_DIRECTION_Z           	-0.26218911324	/**< LLO_4k z-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_LLO_4K_ARM_Y_DIRECTION_X           	0.29774156894	/**< LLO_4k x-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_LLO_4K_ARM_Y_DIRECTION_Y           	-0.48791033647	/**< LLO_4k y-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_LLO_4K_ARM_Y_DIRECTION_Z           	-0.82054461286	/**< LLO_4k z-component of unit vector pointing along y arm in Earth-centered frame */

/**
 * \name VIRGO 3km Interferometric Detector constants
 * The following constants describe the location and geometry of the
 * VIRGO 3km Interferometric Detector.
 */
#define LAL_VIRGO_DETECTOR_NAME               	"VIRGO"	/**< VIRGO detector name string */
#define LAL_VIRGO_DETECTOR_PREFIX             	"V1"	/**< VIRGO detector prefix string */
#define LAL_VIRGO_DETECTOR_LONGITUDE_RAD      	0.18333805213	/**< VIRGO vertex longitude (rad) */
#define LAL_VIRGO_DETECTOR_LATITUDE_RAD       	0.76151183984	/**< VIRGO vertex latitude (rad) */
#define LAL_VIRGO_DETECTOR_ELEVATION_SI       	51.884	/**< VIRGO vertex elevation (m) */
#define LAL_VIRGO_DETECTOR_ARM_X_AZIMUTH_RAD  	0.33916285222	/**< VIRGO x arm azimuth (rad) */
#define LAL_VIRGO_DETECTOR_ARM_Y_AZIMUTH_RAD  	5.05155183261	/**< VIRGO y arm azimuth (rad) */
#define LAL_VIRGO_DETECTOR_ARM_X_ALTITUDE_RAD 	0.00000000000	/**< VIRGO x arm altitude (rad) */
#define LAL_VIRGO_DETECTOR_ARM_Y_ALTITUDE_RAD 	0.00000000000	/**< VIRGO y arm altitude (rad) */
#define LAL_VIRGO_DETECTOR_ARM_X_MIDPOINT_SI  	1500.00000000000	/**< VIRGO x arm midpoint (m) */
#define LAL_VIRGO_DETECTOR_ARM_Y_MIDPOINT_SI  	1500.00000000000	/**< VIRGO y arm midpoint (m) */
#define LAL_VIRGO_VERTEX_LOCATION_X_SI        	4.54637409900e+06	/**< VIRGO x-component of vertex location in Earth-centered frame (m) */
#define LAL_VIRGO_VERTEX_LOCATION_Y_SI        	8.42989697626e+05	/**< VIRGO y-component of vertex location in Earth-centered frame (m) */
#define LAL_VIRGO_VERTEX_LOCATION_Z_SI        	4.37857696241e+06	/**< VIRGO z-component of vertex location in Earth-centered frame (m) */
#define LAL_VIRGO_ARM_X_DIRECTION_X           	-0.70045821479	/**< VIRGO x-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_VIRGO_ARM_X_DIRECTION_Y           	0.20848948619	/**< VIRGO y-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_VIRGO_ARM_X_DIRECTION_Z           	0.68256166277	/**< VIRGO z-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_VIRGO_ARM_Y_DIRECTION_X           	-0.05379255368	/**< VIRGO x-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_VIRGO_ARM_Y_DIRECTION_Y           	-0.96908180549	/**< VIRGO y-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_VIRGO_ARM_Y_DIRECTION_Z           	0.24080451708	/**< VIRGO z-component of unit vector pointing along y arm in Earth-centered frame */

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LLVGEOMETRY_H */
