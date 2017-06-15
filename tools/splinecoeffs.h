/**
 * \author Sylvain Marsat, University of Maryland - NASA GSFC
 *
 * \brief C header for functions computing coeffcients of Not-A-Knot and Quadratic splines in matrix form.
 *
 */

#ifndef _SPLINECOEFFS_H
#define _SPLINECOEFFS_H

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
#include "struct.h"


#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

void BuildNotAKnotSpline(
  gsl_matrix* splinecoeffs,   /* Output: matrix containing all the spline coeffs (already allocated) */
  gsl_vector* vectx,          /* Input: vector x*/
  gsl_vector* vecty,          /* Input: vector y */
  int n);                      /* Size of x, y, and of output matrix */

void BuildQuadSpline(
  gsl_matrix* splinecoeffs,   /* Output: matrix containing all the spline coeffs (already allocated) */
  gsl_vector* vectx,          /* Input: vector x*/
  gsl_vector* vecty,          /* Input: vector y */
  int n);                      /* Size of x, y, and of output matrix */

void BuildSplineCoeffs(
  CAmpPhaseSpline** splines,                  /*  */
  CAmpPhaseFrequencySeries* freqseries);      /*  */

void BuildListmodesCAmpPhaseSpline(
  ListmodesCAmpPhaseSpline** listspline,              /* Output: list of modes of splines in matrix form */
  ListmodesCAmpPhaseFrequencySeries* listh);          /* Input: list of modes in amplitude/phase form */

/* Functions for spline evaluation */

/* Note: for the spines in matrix form, the first column contains the x values, so the coeffs start at 1 */
double EvalCubic(
  gsl_vector* coeffs,  /**/
  double eps,          /**/
  double eps2,         /**/
  double eps3);         /**/

double EvalQuad(
  gsl_vector* coeffs,  /**/
  double eps,          /**/
  double eps2);         /**/

void EvalCAmpPhaseSpline(
  CAmpPhaseSpline* splines,                     //input
  CAmpPhaseFrequencySeries* freqseries);  //in/out defines CAmpPhase from defined freqs  

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _SPLINECOEFFS_H */
