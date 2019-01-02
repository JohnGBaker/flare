#### Try to use cython with gsl functions/vectors matrices

cdef extern from "gsl/gsl_complex.h":

  ctypedef double * gsl_complex_packed_array
  ctypedef double * gsl_complex_packed_ptr

  ctypedef struct gsl_complex:
    double dat[2]



  double GSL_REAL(gsl_complex z) nogil

  double GSL_IMAG(gsl_complex z) nogil

  int GSL_COMPLEX_EQ(gsl_complex z1,gsl_complex z2) nogil

  # double GSL_SET_COMPLEX(gsl_complex * zp, double  x, double  y) nogil

  # double GSL_SET_REAL(gsl_complex * zp, double x) nogil

  # double GSL_SET_IMAG(gsl_complex * zp, double y) nogil
  GSL_SET_COMPLEX(gsl_complex * zp, double  x, double  y) nogil

  GSL_SET_REAL(gsl_complex * zp, double x) nogil

  GSL_SET_IMAG(gsl_complex * zp, double y) nogil


cdef extern from "gsl/gsl_sf_bessel.h":
    double gsl_sf_bessel_J0(double x)

def gsl_sf_bessel_jo(double x):
     return gsl_sf_bessel_J0 (x)

cdef extern from "gsl/gsl_block_double.h":

  ctypedef struct gsl_block:
    size_t size
    double * data

  gsl_block *  gsl_block_alloc(size_t n) nogil

  gsl_block *  gsl_block_calloc(size_t n) nogil

  void  gsl_block_free(gsl_block * b) nogil

  size_t gsl_block_size (gsl_block * b) nogil
  double * gsl_block_data (gsl_block * b) nogil


cdef extern from "gsl/gsl_block_complex_double.h":

  ctypedef struct gsl_block_complex:
    size_t size
    double * data

  gsl_block_complex *  gsl_block_complex_alloc(size_t n) nogil

  gsl_block_complex *  gsl_block_complex_calloc(size_t n) nogil

  size_t gsl_block_complex_size (gsl_block_complex * b) nogil
  double * gsl_block_complex_data (gsl_block_complex * b) nogil


cdef extern from "gsl/gsl_vector.h":

  ctypedef struct gsl_vector:
    size_t size
    size_t stride
    double *data
    gsl_block *block
    int owner

  ctypedef struct gsl_vector_view:
    gsl_vector vector

  ctypedef struct gsl_vector_const_view:
    gsl_vector vector


  # Allocation
  gsl_vector *  gsl_vector_alloc(size_t n) nogil

  gsl_vector *  gsl_vector_calloc(size_t n) nogil

  gsl_vector_alloc_from_block(gsl_block * b, size_t offset,
                              size_t n, size_t stride) nogil

  gsl_vector *gsl_vector_alloc_from_vector(gsl_vector * v,
                         size_t offset, size_t n, size_t stride) nogil

  void  gsl_vector_free(gsl_vector * v) nogil

  # Views
  gsl_vector_view  gsl_vector_view_array(double *base, size_t n) nogil

  gsl_vector_view  gsl_vector_subvector(gsl_vector *v, size_t offset, size_t n) nogil

  gsl_vector_view  gsl_vector_view_array_with_stride(double * base, size_t stride, size_t n) nogil

  gsl_vector_const_view  gsl_vector_const_view_array(double *base, size_t n) nogil

  gsl_vector_const_view  gsl_vector_const_view_array_with_stride(double * base, size_t stride, size_t n) nogil

  gsl_vector_const_view  gsl_vector_const_subvector(gsl_vector * v, size_t offset, size_t n) nogil

  gsl_vector_view  gsl_vector_subvector_with_stride(gsl_vector *v, size_t offset, size_t stride, size_t n) nogil

  gsl_vector_const_view  gsl_vector_const_subvector_with_stride(gsl_vector * v, size_t offset, size_t stride, size_t n) nogil


  # Operations
  double  gsl_vector_get(gsl_vector * v, size_t i) nogil

  void  gsl_vector_set(gsl_vector * v, size_t i, double x) nogil

  double *  gsl_vector_ptr(gsl_vector * v, size_t i) nogil

  double *  gsl_vector_const_ptr(gsl_vector * v, size_t i) nogil

  void  gsl_vector_set_zero(gsl_vector * v) nogil

  void  gsl_vector_set_all(gsl_vector * v, double x) nogil

  int  gsl_vector_set_basis(gsl_vector * v, size_t i) nogil

  # Reading and writing vectors
  # int  gsl_vector_fread(FILE * stream, gsl_vector * v) nogil
  #
  # int  gsl_vector_fwrite(FILE * stream, gsl_vector * v) nogil
  #
  # int  gsl_vector_fscanf(FILE * stream, gsl_vector * v) nogil
  #
  # int  gsl_vector_fprintf(FILE * stream, gsl_vector * v, char * format) nogil

  # Copying or exchanging elements
  int  gsl_vector_memcpy(gsl_vector * dest, gsl_vector * src) nogil

  int  gsl_vector_reverse(gsl_vector * v) nogil

  int  gsl_vector_swap(gsl_vector * v, gsl_vector * w) nogil

  int  gsl_vector_swap_elements(gsl_vector * v, size_t i, size_t j) nogil

  # Finding maximum and minimum elements of vectors

  double  gsl_vector_max(gsl_vector * v) nogil

  double  gsl_vector_min(gsl_vector * v) nogil

  void  gsl_vector_minmax(gsl_vector * v, double * min_out, double * max_out) nogil

  size_t  gsl_vector_max_index(gsl_vector * v) nogil

  size_t  gsl_vector_min_index(gsl_vector * v) nogil

  void  gsl_vector_minmax_index(gsl_vector * v, size_t * imin, size_t * imax) nogil

  # Vector operations
  int  gsl_vector_add(gsl_vector * a, gsl_vector * b) nogil

  int  gsl_vector_sub(gsl_vector * a, gsl_vector * b) nogil

  int  gsl_vector_mul(gsl_vector * a, gsl_vector * b) nogil

  int  gsl_vector_div(gsl_vector * a, gsl_vector * b) nogil

  int  gsl_vector_scale(gsl_vector * a, double x) nogil

  int  gsl_vector_add_constant(gsl_vector * a, double x) nogil

  int  gsl_vector_isnull(gsl_vector * v) nogil



cdef extern from "gsl/gsl_matrix_double.h":
    ctypedef struct gsl_matrix:
        size_t size1
        size_t size2
        size_t tda
        double * data
        gsl_block * block
        int owner

    ctypedef struct gsl_matrix_view:
        gsl_matrix matrix

    ctypedef struct gsl_matrix_const_view

    # Allocation
    gsl_matrix *  gsl_matrix_alloc(size_t n1, size_t n2) nogil

    gsl_matrix *  gsl_matrix_calloc(size_t n1, size_t n2) nogil

    gsl_matrix *  gsl_matrix_alloc_from_block(gsl_block * b,
    size_t offset, size_t n1, size_t n2, size_t d2) nogil

    gsl_matrix * gsl_matrix_alloc_from_matrix (gsl_matrix * m,  size_t k1,  size_t k2,  size_t n1,  size_t n2) nogil

    gsl_vector * gsl_vector_alloc_row_from_matrix (gsl_matrix * m,  size_t i) nogil

    gsl_vector * gsl_vector_alloc_col_from_matrix (gsl_matrix * m,  size_t j) nogil

    void  gsl_matrix_free(gsl_matrix * m) nogil

    # Views
    gsl_matrix_view  gsl_matrix_submatrix(gsl_matrix * m, size_t k1, size_t k2, size_t n1, size_t n2) nogil

    gsl_vector_view  gsl_matrix_row(gsl_matrix * m, size_t i) nogil

    gsl_vector_view  gsl_matrix_column(gsl_matrix * m, size_t j) nogil

    gsl_vector_view  gsl_matrix_diagonal(gsl_matrix * m) nogil

    gsl_vector_view  gsl_matrix_subdiagonal(gsl_matrix * m, size_t k) nogil

    gsl_vector_view  gsl_matrix_superdiagonal(gsl_matrix * m, size_t k) nogil

    gsl_matrix_view  gsl_matrix_view_array(double * base, size_t n1, size_t n2) nogil

    gsl_matrix_view  gsl_matrix_view_array_with_tda(double * base, size_t n1, size_t n2, size_t tda) nogil

    gsl_matrix_view  gsl_matrix_view_vector(gsl_vector * v, size_t n1, size_t n2) nogil

    gsl_matrix_view  gsl_matrix_view_vector_with_tda(gsl_vector * v, size_t n1, size_t n2, size_t tda) nogil

    gsl_matrix_const_view  gsl_matrix_const_submatrix(gsl_matrix * m, size_t k1, size_t k2, size_t n1, size_t n2) nogil

    gsl_vector_const_view  gsl_matrix_const_row(gsl_matrix * m, size_t i) nogil

    gsl_vector_const_view  gsl_matrix_const_column(gsl_matrix * m, size_t j) nogil

    gsl_vector_const_view  gsl_matrix_const_diagonal(gsl_matrix * m) nogil

    gsl_vector_const_view  gsl_matrix_const_subdiagonal(gsl_matrix * m, size_t k) nogil

    gsl_vector_const_view  gsl_matrix_const_superdiagonal(gsl_matrix * m, size_t k) nogil

    gsl_matrix_const_view  gsl_matrix_const_view_array(double * base, size_t n1, size_t n2) nogil

    gsl_matrix_const_view  gsl_matrix_const_view_array_with_tda(double * base, size_t n1, size_t n2, size_t tda) nogil

    gsl_matrix_const_view  gsl_matrix_const_view_vector(gsl_vector * v, size_t n1, size_t n2) nogil

    gsl_matrix_const_view  gsl_matrix_const_view_vector_with_tda(gsl_vector * v, size_t n1, size_t n2, size_t tda) nogil


    # Operations
    double  gsl_matrix_get(gsl_matrix * m, size_t i, size_t j) nogil

    void  gsl_matrix_set(gsl_matrix * m, size_t i, size_t j, double x) nogil

    double *  gsl_matrix_ptr(gsl_matrix * m, size_t i, size_t j) nogil

    double *  gsl_matrix_const_ptr(gsl_matrix * m, size_t i, size_t j) nogil

    void  gsl_matrix_set_zero(gsl_matrix * m) nogil

    void  gsl_matrix_set_identity(gsl_matrix * m) nogil

    void  gsl_matrix_set_all(gsl_matrix * m, double x) nogil

    # Reading and writing matrices
    # int  gsl_matrix_fread(FILE * stream, gsl_matrix * m) nogil
    #
    # int  gsl_matrix_fwrite(FILE * stream, gsl_matrix * m) nogil
    #
    # int  gsl_matrix_fscanf(FILE * stream, gsl_matrix * m) nogil
    #
    # int  gsl_matrix_fprintf(FILE * stream, gsl_matrix * m, char * format) nogil

    # Copying or exchanging elements
    int  gsl_matrix_memcpy(gsl_matrix * dest, gsl_matrix * src) nogil

    int  gsl_matrix_swap(gsl_matrix * m1, gsl_matrix * m2) nogil

    int  gsl_matrix_swap_rows(gsl_matrix * m, size_t i, size_t j) nogil

    int  gsl_matrix_swap_columns(gsl_matrix * m, size_t i, size_t j) nogil

    int  gsl_matrix_swap_rowcol(gsl_matrix * m, size_t i, size_t j) nogil

    int  gsl_matrix_transpose(gsl_matrix * m) nogil

    int  gsl_matrix_transpose_memcpy(gsl_matrix * dest, gsl_matrix * src) nogil

    # Finding maximum and minimum elements of matrices
    double  gsl_matrix_max(gsl_matrix * m) nogil

    double  gsl_matrix_min(gsl_matrix * m) nogil

    void  gsl_matrix_minmax(gsl_matrix * m, double * min_out, double * max_out) nogil

    void  gsl_matrix_max_index(gsl_matrix * m, size_t * imax, size_t * jmax) nogil

    void  gsl_matrix_min_index(gsl_matrix * m, size_t * imax, size_t * jmax) nogil

    void  gsl_matrix_minmax_index(gsl_matrix * m, size_t * imin, size_t * jmin, size_t * imax, size_t * jmax) nogil

    int  gsl_matrix_isnull(gsl_matrix * m) nogil

    # Matrix operations
    int  gsl_matrix_add(gsl_matrix * a, gsl_matrix * b) nogil

    int  gsl_matrix_sub(gsl_matrix * a, gsl_matrix * b) nogil

    int  gsl_matrix_mul_elements(gsl_matrix * a, gsl_matrix * b) nogil

    int  gsl_matrix_div_elements(gsl_matrix * a, gsl_matrix * b) nogil

    int  gsl_matrix_scale(gsl_matrix * a, double x) nogil

    int  gsl_matrix_add_constant(gsl_matrix * a, double x) nogil

    int gsl_matrix_add_diagonal (gsl_matrix * a,  double x) nogil

    # The functions below are obsolete
    int  gsl_matrix_get_row(gsl_vector * v, gsl_matrix * m, size_t i) nogil

    int  gsl_matrix_get_col(gsl_vector * v, gsl_matrix * m, size_t j) nogil

    int  gsl_matrix_set_row(gsl_matrix * m, size_t i, gsl_vector * v) nogil

    int  gsl_matrix_set_col(gsl_matrix * m, size_t j, gsl_vector * v) nogil
