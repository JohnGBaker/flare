//Cubic spline with derivatives
//Written by John G Baker at gsfc.nasa.gov
//SWIG code for python wrapping

%module splinecspagh
%{
#include "spline.h"
#include <stdio.h>
  //int splinecspagh_nsize;  //we use this to globally save the array size when the input is prepared

%}
//This provides the output pointer handling
%include "typemaps.i"
 
extern int spline_natural;
extern int spline_iold;

//Create a rule for setting up input arrays 
//Probably there is a better way, but if so, it wasn't easy to find.
//Might even be easier just to do the wrapping by hand... 

//%typemap(arginit) int n { 
//First we try to save the value of the int n argument for use in the conversion of the other variables
// we are hijacking "varinit" for this since the resulting wrapper code should come before 
// the conversion code in "in" so hopefully this works..
// -I only know how to do this using a global variable 
//splinecspagh_nsize=arg$argnum;
//must initialize pointers to avoid annoying warning messages
//}
//%typemap(arginit) double* {
  //varinit
  //$1=NULL;
//}
%typemap(in) double *xs, double *ys, double *zs{
  int i;
  int nsize=0;
  SWIG_AsVal_int(obj4, &nsize); //here obj4 should reference the 5th argument.  We don't do checking, but that will be done by the regular "in" method
  //printf("in: %s var %i, nsize=%i\n","$symname",$argnum,nsize);
  if (!PySequence_Check($input)) {
    PyErr_SetString(PyExc_ValueError,"Expected a sequence");
    return NULL;
  }
  $1 = (double *) malloc(nsize*sizeof(double));
  for (i = 0; i < nsize; i++) {
    PyObject *o = PySequence_GetItem($input,i);
    $1[i] = (double) PyFloat_AsDouble(o);
  }
}
//must clean up memory allocation
%typemap(freearg) double *xs, double *ys, double *zs {
       if ($1) free($1);
    }
extern void spline_intd3(double xp, const double *xs,const double *ys,const double *zs,int n,double *OUTPUT, double *OUTPUT, double *OUTPUT, double *OUTPUT);




//same as above with difference position for int n
%typemap(in) double *xs, double *ys, double *zs{
  int i;
  int nsize=0;
  SWIG_AsVal_int(obj5, &nsize); //here obj5 should reference the 6th argument.  We don't do checking, but that will be done by the regular "in" method
  //printf("in: %s var %i, nsize=%i\n","$symname",$argnum,nsize);
  if (!PySequence_Check($input)) {
    PyErr_SetString(PyExc_ValueError,"Expected a sequence");
    return NULL;
  }
  $1 = (double *) malloc(nsize*sizeof(double));
  for (i = 0; i < nsize; i++) {
    PyObject *o = PySequence_GetItem($input,i);
    $1[i] = (double) PyFloat_AsDouble(o);
  }
}
//same as above with difference position for int n and no need to copy anything in
//we can't use numinputs=0, which would leave them sensible out of input arg list, because then this call would come before
//obj1 is instantiated and we wouldn't be able to allocate the output array
%typemap(in) double *q, double *dq1, double *dq2, double *dq3 {
  int nsize=0;
  SWIG_AsVal_int(obj1, &nsize); //here obj1 should reference the 2nd argument.  We don't do checking, but that will be done by the regular "in" method
  //printf("in: %s var %i, nsize=%i\n","$symname",$argnum,nsize);
  $1 = (double *) malloc(nsize*sizeof(double));
}
%typemap(in) double *xp{
  int i;
  int nsize=0;
  SWIG_AsVal_int(obj1, &nsize); //here obj1 should reference the 2nd argument.  We don't do checking, but that will be done by the regular "in" method
  //printf("in: %s var %i, nsize=%i\n","$symname",$argnum,nsize);
  if (!PySequence_Check($input)) {
    PyErr_SetString(PyExc_ValueError,"Expected a sequence");
    return NULL;
  }
  $1 = (double *) malloc(nsize*sizeof(double));
  for (i = 0; i < nsize; i++) {
    PyObject *o = PySequence_GetItem($input,i);
    $1[i] = (double) PyFloat_AsDouble(o);
  }
}
//Create a rule for setting up non-const arrays as output argument arrays 
%typemap(argout) double *q, double *dq1, double *dq2, double *dq3 {
  int nsize=0;
  SWIG_AsVal_int(obj1, &nsize); //here obj1 should reference the 2nd argument.  We don't do checking, but that will be done by the regular "in" method
  PyObject *o = PyList_New(nsize);
  int i;
  for(i=0; i<nsize; i++)
    {
      //printf("SWIG setting z[%i]=%g\n",i,$1[i]);
      PyList_SetItem(o, i, PyFloat_FromDouble($1[i]));
    }
  //printf("o=:\n%s\n",PyString_AsString(PyObject_Str(o)));
  //printf("tuplecheck=%i\n",PyTuple_Check($result));
  if (!PyList_Check($result)) {//have to make the result a List if it is not yet
    $result = PyList_New(1);
    PyList_SetItem($result,0,o);
  } else {
    %append_output(o);
  }
  //printf("o=:\n%s\n",PyString_AsString(PyObject_Str(o)));
  //printf("result->:\n%s\n",PyString_AsString(PyObject_Str($result)));
  
}
//must clean up memory allocation
%typemap(freearg) double *xp, double *q, double *dq1, double *dq2, double *dq3 {
       if ($1) free($1);
}
void spline_intd3_vec(const double *xp,int nxp, const double *xs,const double *ys,const double *zs,int n, double *q, double *dq1, double *dq2, double *dq3);



//same as above, but we get nsize from arg 4
%typemap(in) double *xs, double *ys, double *zs{
  int i;
  int nsize=0;
  SWIG_AsVal_int(obj3, &nsize); //here obj3 should reference the 4th argument.  We don't do checking, but that will be done by the regular "in" method

  //printf("in: %s var %i, nsize=%i\n","$symname",$argnum,nsize);
  if (!PySequence_Check($input)) {
    PyErr_SetString(PyExc_ValueError,"Expected a sequence");
    return NULL;
  }
  $1 = (double *) malloc(nsize*sizeof(double));
  for (i = 0; i < nsize; i++) {
    PyObject *o = PySequence_GetItem($input,i);
    $1[i] = (double) PyFloat_AsDouble(o);
  }
}
//Create a rule for setting up non-const arrays as output argument arrays 
%typemap(argout) double *zs {
  int nsize=0;
  SWIG_AsVal_int(obj3, &nsize); //here obj3 should reference the 4th argument.  We don't do checking, but that will be done by the regular "in" method
  PyObject *o = PyList_New(nsize);
  int i;
  for(i=0; i<nsize; i++)
    {
      //printf("SWIG setting z[%i]=%g\n",i,$1[i]);
      PyList_SetItem(o, i, PyFloat_FromDouble($1[i]));
    }
  //printf("o=:\n%s\n",PyString_AsString(PyObject_Str(o)));
  $result = o;
}

extern void spline_construct(const double *xs,const double *ys,double *zs,int n);

        
