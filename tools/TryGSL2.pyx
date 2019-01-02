from cython_gsl cimport *
import numpy as np
cimport numpy as np

# cdef extern from "gsl/gsl_complex.h":
#
#   ctypedef double * gsl_complex_packed_array
#   ctypedef double * gsl_complex_packed_ptr
#
#   ctypedef struct gsl_complex:
#     double dat[2]
#
#
#
#   double GSL_REAL(gsl_complex z) nogil
#
#   double GSL_IMAG(gsl_complex z) nogil
#
#   int GSL_COMPLEX_EQ(gsl_complex z1,gsl_complex z2) nogil
#
#   void GSL_SET_COMPLEX(gsl_complex * zp, double  x, double  y);
#
#   void GSL_SET_REAL(gsl_complex * zp, double x);
#
#   void GSL_SET_IMAG(gsl_complex * zp, double y);


def g_complex(a, b):
    cdef double re
    cdef double im
    re = a
    im = b
    cdef gsl_complex c
    GSL_SET_COMPLEX(&c, re, im)

    print ("We have:", GSL_REAL(c), GSL_IMAG(c))

    return(gsl_complex_abs(c), gsl_complex_arg(c))

def g_bess(double x):
    return gsl_sf_bessel_J0 (x)


# def Vec_Arr(gsl_vector* v):
#     cdef int N
#     cdef int i
#     N = v.size
#     p = np.zeros(N)
#
#     for i from 0<= i <N:
#         p[i] = gsl_vector_get (v, i)
#
#     return(p)


def Arr_Vec(np.ndarray a):

    cdef int N
    cdef int i
    cdef gsl_vector* v
    N = len(a)

    v = gsl_vector_alloc(N)

    for i in range(N):
        gsl_vector_set(v, i, a[i])

    b = np.zeros(N)

    for i in range(N):
        b[i] = gsl_vector_get (v, i)*2.0

    return(b)

def Arr_Vec2(np.ndarray[np.float_t, ndim=1] a):

    cdef int N
    cdef int i
    cdef gsl_vector_view v
    N = len(a)
    dat = <double *> a.data

    v = gsl_vector_view_array(dat,  N)

    ### Let's do some operation with vector
    for i in range(N):
        gsl_vector_set(&(v.vector), i, gsl_vector_get (&(v.vector), i)*3.0)


    b = np.zeros(N)
    cdef double [:] b_view = b
    cdef double [:] v_view = <(double)[:N]>(v.vector).data
    # cdef double [::1] v_view = <(double)[:N]>(v.vector).data

    b = np.asarray(v_view)

    # b.data = (v.vector).data

    # for i in range(N):
    #     print( gsl_vector_get (&(v.vector), i) )


    return(b)
