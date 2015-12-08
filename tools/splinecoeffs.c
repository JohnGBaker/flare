// build commmand:
// gcc -c test_fresnel.c -O2 -std=c99 -I../LISAsim -I/opt/local/include
// gcc -o test_fresnel test_fresnel.o struct.o ../LISAsim/LISAnoise.o -lgsl -lgslcblas

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_complex.h>

#include "constants.h"
#include "struct.h"
#include "splinecoeffs.h"


/* Implementation of the Thomas algorithm to solve a tridiagonal system */
/* Note: assumes numerical stability (e.g. diagonal-dominated matrix) */
static void SolveTridiagThomas(
  gsl_vector* vectx,        /* Output: solution vector, length n - already allocated */
  gsl_vector* vecta,        /* Diagonal of the matrix, length n */
  gsl_vector* vectb,        /* Lower diagonal, length n-1 */
  gsl_vector* vectc,        /* Upper diagonal, length n-1 */
  gsl_vector* vecty,        /* Right-hand-side vector, length n */
  int n)                    /* Length of vectx, vecty */
{
  /* Check lengths */
  if(!(vectx->size==n && vecty->size==n && vecta->size==n && vectb->size==n-1 && vectc->size==n-1)) {
    printf("Error: incompatible lengths in SolveTridiagThomas.\n");
    exit(1);
  }
  
  /* Note: we modify the vectors c, y in place */
  double* chat = vectc->data;
  double* yhat = vecty->data;
  double* x = vectx->data;
  double* a = vecta->data;
  double* b = vectb->data;

  /* Sweep forward, computing the chat and yhat values */
  chat[0] = chat[0] / a[0];
  yhat[0] = yhat[0] / a[0];
  double factor;
  for(int i=1; i<=n-2; i++) {
    factor = 1./(a[i] - b[i-1]*chat[i-1]);
    chat[i] = chat[i] * factor;
    yhat[i] = (yhat[i] - b[i-1]*yhat[i-1]) * factor;
  }
  factor = 1./(a[n-1] - b[n-2]*chat[n-2]);
  yhat[n-1] = (yhat[n-1] - b[n-2]*yhat[n-2]) * factor;

  /* Solve for x going backward */
  x[n-1] = yhat[n-1];
  for(int i=n-2; i>=0; i--) {
    x[i] = yhat[i] - chat[i] * x[i+1];
  }
}

void BuildNotAKnotSpline(
  gsl_matrix* splinecoeffs,   /* Output: matrix containing all the spline coeffs (already allocated) */
  gsl_vector* vectx,          /* Input: vector x*/
  gsl_vector* vecty,          /* Input: vector y */
  int n)                      /* Size of x, y, and of output matrix */
{
  /* Check lengths */
  if(!(vectx->size==n && vecty->size==n && splinecoeffs->size1==n && splinecoeffs->size2==5)) {
    printf("Error: incompatible lengths in NotAKnotSpline.\n");
    exit(1);
  }
  double* x = vectx->data;
  double* y = vecty->data;

  /* Computing vecth and vectDeltay */
  gsl_vector* vecth = gsl_vector_alloc(n-1);
  gsl_vector* vectDeltay = gsl_vector_alloc(n-1);
  gsl_vector* vectDeltayoverh = gsl_vector_alloc(n-1);
  double* h = vecth->data;
  double* Deltay = vectDeltay->data;
  double* Deltayoverh = vectDeltayoverh->data;
  for(int i=0; i<n-1; i++) {
    h[i] = x[i+1] - x[i];
    Deltay[i] = y[i+1] - y[i];
    Deltayoverh[i] = Deltay[i] / h[i];
  }

  /* Structures for the tridiagonal system */
  gsl_vector* vectY = gsl_vector_alloc(n-2);
  gsl_vector* vecta = gsl_vector_alloc(n-2);
  gsl_vector* vectb = gsl_vector_alloc(n-3);
  gsl_vector* vectc = gsl_vector_alloc(n-3);
  double* Y = vectY->data;
  double* a = vecta->data;
  double* b = vectb->data;
  double* c = vectc->data;
  for(int i=0; i<=n-3; i++) {
    Y[i] = 3.*(Deltayoverh[i+1] - Deltayoverh[i]);
    a[i] = 2.*(h[i+1] + h[i]);
  }
  for(int i=0; i<=n-4; i++) {
    b[i] = h[i+1];
    c[i] = h[i+1];
  }
  /* Adjusting for the not-a-knot condition */
  a[0] += h[0] + h[0]*h[0]/h[1];
  c[0] += -h[0]*h[0]/h[1];
  a[n-3] += h[n-2] + h[n-2]*h[n-2]/h[n-3];
  b[n-4] += -h[n-2]*h[n-2]/h[n-3];
  
  /* Solving the tridiagonal system */
  gsl_vector* vectp2 = gsl_vector_alloc(n);
  gsl_vector_view viewp2trunc = gsl_vector_subvector(vectp2, 1, n-2);
  SolveTridiagThomas(&viewp2trunc.vector, vecta, vectb, vectc, vectY, n-2);
  double* p2 = vectp2->data;
  p2[0] = p2[1] - h[0]/h[1] * (p2[2] - p2[1]);
  p2[n-1] = p2[n-2] + h[n-2]/h[n-3] * (p2[n-2] - p2[n-3]);

  /* Deducing the p1's and the p3's */
  gsl_vector* vectp1 = gsl_vector_alloc(n);
  gsl_vector* vectp3 = gsl_vector_alloc(n);
  double* p1 = vectp1->data;
  double* p3 = vectp3->data;
  for(int i=0; i<=n-2; i++) {
    p1[i] = Deltayoverh[i] - h[i]/3. * (p2[i+1] + 2.*p2[i]);
    p3[i] = (p2[i+1] - p2[i]) / (3*h[i]);
  }
  /* Note: p1[n-1], p2[n-1], p3[n-1] are set to values coherent with the derivatives of the spline at the last point, but they are not stricly speaking coefficients of the spline. */
  p1[n-1] = p1[n-2] + 2.*p2[n-2]*h[n-2] + 3.*p3[n-2]*h[n-2]*h[n-2];
  p3[n-1] = p3[n-2];

  /* Copying the results in the output matrix */
  gsl_vector_view viewx = gsl_matrix_column(splinecoeffs, 0);
  gsl_vector_view viewp0 = gsl_matrix_column(splinecoeffs, 1);
  gsl_vector_view viewp1 = gsl_matrix_column(splinecoeffs, 2);
  gsl_vector_view viewp2 = gsl_matrix_column(splinecoeffs, 3);
  gsl_vector_view viewp3 = gsl_matrix_column(splinecoeffs, 4);
  gsl_vector_memcpy(&viewx.vector, vectx);
  gsl_vector_memcpy(&viewp0.vector, vecty);
  gsl_vector_memcpy(&viewp1.vector, vectp1);
  gsl_vector_memcpy(&viewp2.vector, vectp2);
  gsl_vector_memcpy(&viewp3.vector, vectp3);

  /* Cleanup*/
  gsl_vector_free(vecth);
  gsl_vector_free(vectDeltay);
  gsl_vector_free(vectDeltayoverh);
  gsl_vector_free(vectY);
  gsl_vector_free(vecta);
  gsl_vector_free(vectb);
  gsl_vector_free(vectc);
  gsl_vector_free(vectp1);
  gsl_vector_free(vectp2);
  gsl_vector_free(vectp3);
}

void BuildQuadSpline(
  gsl_matrix* splinecoeffs,   /* Output: matrix containing all the spline coeffs (already allocated) */
  gsl_vector* vectx,          /* Input: vector x*/
  gsl_vector* vecty,          /* Input: vector y */
  int n)                      /* Size of x, y, and of output matrix */
{
  /* Check lengths */
  if(!(vectx->size==n && vecty->size==n && splinecoeffs->size1==n && splinecoeffs->size2==4)) {
    printf("Error: incompatible lengths in NotAKnotSpline.\n");
    exit(1);
  }
  double* x = vectx->data;
  double* y = vecty->data;

  /* Computing vecth and vectDeltay */
  gsl_vector* vecth = gsl_vector_alloc(n-1);
  gsl_vector* vectDeltay = gsl_vector_alloc(n-1);
  gsl_vector* vectDeltayoverh = gsl_vector_alloc(n-1);
  double* h = vecth->data;
  double* Deltay = vectDeltay->data;
  double* Deltayoverh = vectDeltayoverh->data;
  for(int i=0; i<n-1; i++) {
    h[i] = x[i+1] - x[i];
    Deltay[i] = y[i+1] - y[i];
    Deltayoverh[i] = Deltay[i] / h[i];
  }
  
  /* Solving for p1 */
  gsl_vector* vectp1 = gsl_vector_alloc(n);
  double* p1 = vectp1->data;
  double ratio = h[n-2] / h[n-3];
  p1[n-3] = ((2. + ratio)*Deltayoverh[n-3] - Deltayoverh[n-2]) / (1. + ratio);
  p1[n-2] = -p1[n-3] + 2.*Deltayoverh[n-3];
  for(int i=n-4; i>=0; i--) {
    p1[i] = -p1[i+1] + 2.*Deltayoverh[i];
  }
  p1[n-1] = (1. + ratio)*p1[n-2] - ratio*p1[n-3];

  /* Deducing the p2's */
  gsl_vector* vectp2 = gsl_vector_alloc(n);
  double* p2 = vectp2->data;
  for(int i=0; i<=n-2; i++) {
    p2[i] = (p1[i+1] - p1[i]) / (2.*h[i]);
  }
  /* Note: p2[n-1] is set to values coherent with the derivatives of the spline at the last point, but not stricly speaking a coefficient of the spline. */
  p2[n-1] = p2[n-2];

  /* Copying the results in the output matrix */
  gsl_vector_view viewx = gsl_matrix_column(splinecoeffs, 0);
  gsl_vector_view viewp0 = gsl_matrix_column(splinecoeffs, 1);
  gsl_vector_view viewp1 = gsl_matrix_column(splinecoeffs, 2);
  gsl_vector_view viewp2 = gsl_matrix_column(splinecoeffs, 3);
  gsl_vector_memcpy(&viewx.vector, vectx);
  gsl_vector_memcpy(&viewp0.vector, vecty);
  gsl_vector_memcpy(&viewp1.vector, vectp1);
  gsl_vector_memcpy(&viewp2.vector, vectp2);

  /* Cleanup*/
  gsl_vector_free(vecth);
  gsl_vector_free(vectDeltay);
  gsl_vector_free(vectDeltayoverh);
  gsl_vector_free(vectp1);
  gsl_vector_free(vectp2);
}

void BuildSplineCoeffs(
  CAmpPhaseSpline** splines,                  /*  */
  CAmpPhaseFrequencySeries* freqseries)       /*  */
{
  /* Initialize output structure */
  int n = (int) freqseries->freq->size;
  CAmpPhaseSpline_Init(splines, n);

  /* Build the splines */
  BuildNotAKnotSpline((*splines)->spline_amp_real, freqseries->freq, freqseries->amp_real, n);
  BuildNotAKnotSpline((*splines)->spline_amp_imag, freqseries->freq, freqseries->amp_imag, n);
  BuildQuadSpline((*splines)->quadspline_phase, freqseries->freq, freqseries->phase, n);
}

void BuildListmodesCAmpPhaseSpline(
  ListmodesCAmpPhaseSpline** listspline,              /* Output: list of modes of splines in matrix form */
  ListmodesCAmpPhaseFrequencySeries* listh)           /* Input: list of modes in amplitude/phase form */
{
    if(*listspline){ /* We don't allow for the case where listspline already points to something */
      printf("Error: Tried to add a mode to an already existing ListmodesCAmpPhaseSpline ");
      exit(1);
    }
    else {
      ListmodesCAmpPhaseFrequencySeries* listelementh = listh;
      while(listelementh) {
	CAmpPhaseSpline* splines = NULL;
	BuildSplineCoeffs(&splines, listelementh->freqseries);
	int l = listelementh->l;
	int m = listelementh->m;
	*listspline = ListmodesCAmpPhaseSpline_AddModeNoCopy(*listspline, splines, l, m);
	listelementh = listelementh->next;	
      }
    }
}
