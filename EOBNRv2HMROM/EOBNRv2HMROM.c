/**
 * \author Sylvain Marsat, University of Maryland - NASA GSFC
 *
 * \brief C code for EOBNRv2HM reduced order model (non-spinning version).
 * See CQG 31 195010, 2014, arXiv:1402.4146 for details on the reduced order method.
 * See arXiv:1106.1021 for the EOBNRv2HM model.
 *  
 * Borrows from the SEOBNR ROM LAL code written by Michael Puerrer and John Veitch.
 *
 * Put the untared data in the directory designated by the environment variable ROM_DATA_PATH.
 *
 * Parameter ranges:
 *   q = 1-12 (almost)
 *   No spin
 *   Mtot >= 10Msun for fstart=8Hz
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
#include "struct.h"
#include "EOBNRv2HMROMstruct.h"
#include "EOBNRv2HMROM.h"

/*************************************************/
/********* Some general definitions **************/

/* Number and list of modes to generate - to be modified ultimately to allow for a selection of the desired mode(s) */
/* By convention the first mode of the list is used to set phiRef */
/* nbmodemax = 5 has been defined in the header */
const int listmode[nbmodemax][2] = { {2,2}, {2,1}, {3,3}, {4,4}, {5,5} };

/* Maximal mass ratio covered by the model - q=12 (almost) for now */ 
static const double q_max = 11.9894197212;
/* Minimal geometric frequency covered by the model - f=8Hz for M=10Msol for now */ 
static const double Mf_ROM_min = 0.0003940393857519091;
/* Maximal geometric frequency covered by the model - Mf=0.285 for the 55 mode - used as default */ 
static const double Mf_ROM_max = 0.285;
/* Total mass (in units of solar mass) used to generate the ROM - used to restore the correct amplitude (global factor M) when coming back to physical units */
static const double M_ROM = 10.;

/* Define the number of points in frequency used by the SVD, identical for all modes */
static const int nbfreq = 300;
/* Define the number of training waveforms used by the SVD, identical for all modes */
static const int nbwf = 301;

/******************************************************************/
/********* Definitions for the persistent structures **************/

/* SEOBNR-ROM structures are generalized to lists */
ListmodesEOBNRHMROMdata* __EOBNRv2HMROM_data_init = NULL; /* for initialization only */
ListmodesEOBNRHMROMdata** const __EOBNRv2HMROM_data = &__EOBNRv2HMROM_data_init;
ListmodesEOBNRHMROMdata_interp* __EOBNRv2HMROM_interp_init = NULL; /* for initialization only */
ListmodesEOBNRHMROMdata_interp** const __EOBNRv2HMROM_interp = &__EOBNRv2HMROM_interp_init;
/* Tag indicating whether the data has been loaded and interpolated */
int __EOBNRv2HMROM_setup = FAILURE; /* To be set to SUCCESS after initialization*/

/********************* Miscellaneous ********************/

/* Return the closest higher power of 2 */
static size_t NextPow2(const size_t n) {
  return 1 << (size_t) ceil(log2(n));
}

/* Arbitrary tuned q-dependent functions by which the frequencies for the 44 and 55 modes have been multiplied (to put the ringdowns at the same Mf). The frequencies stored in the data for the 44 and 55 modes are the rescaled ones, not the original ones. */
static double Scaling44(const double q) {
  return 1.-1./4.*exp(-(q-1.)/5.);
}
static double Scaling55(const double q) {
  return 1.-1./3.*exp(-(q-1.)/5.);
}

/* Function evaluating eta as a function of q */
static double EtaOfq(const double q) {
  return q/(1.+q)/(1.+q);
}
/* Function evaluating delta m/m = (m1-m2)/(m1+m2) as a function of q*/
static double DeltaOfq(const double q) {
  return( (q-1.)/(q+1.) );
}

/* Amplitude factors scaled out for each mode; this does not include the global amplitude factor scaled out from all modes. */
static double ModeAmpFactor(const int l, const int m, const double q) {
  double eta = EtaOfq(q);
  if( l==2 && m==2 ) return(sqrt(eta));
  else if( l==2 && m==1 ) return( sqrt(eta)*DeltaOfq(q) );
  else if( l==3 && m==3 ) return( sqrt(eta)*DeltaOfq(q) );
  else if( l==4 && m==4 ) return( sqrt(Scaling44(q))*sqrt(eta)*(1.-3.*eta) );
  else if( l==5 && m==5 ) return( pow(Scaling55(q), 1./6)*sqrt(eta)*DeltaOfq(q)*(1.-2.*eta) );
  else {
    fprintf(stderr, "Unknown mode (%d,%d) for the amplitude factor.\n", l, m);
    return(FAILURE);
  }
}

/********************* Functions to initialize and cleanup contents of data structures ********************/
void EOBNRHMROMdata_Init(EOBNRHMROMdata **data) {
  if(!data) exit(1);
  /* Create storage for structures */
  if(!*data) *data=malloc(sizeof(EOBNRHMROMdata));
  else
  {
    EOBNRHMROMdata_Cleanup(*data);
  }
  gsl_set_error_handler(&Err_Handler);
  (*data)->q = gsl_vector_alloc(nbwf);
  (*data)->freq = gsl_vector_alloc(nbfreq);
  (*data)->Camp = gsl_matrix_alloc(nk_amp,nbwf);
  (*data)->Cphi = gsl_matrix_alloc(nk_phi,nbwf);
  (*data)->Bamp = gsl_matrix_alloc(nbfreq,nk_amp);
  (*data)->Bphi = gsl_matrix_alloc(nbfreq,nk_phi);
  (*data)->shifttime = gsl_vector_alloc(nbwf);
  (*data)->shiftphase = gsl_vector_alloc(nbwf);
}
void EOBNRHMROMdata_interp_Init(EOBNRHMROMdata_interp **data_interp) {
  if(!data_interp) exit(1);
  /* Create storage for structures */
  if(!*data_interp) *data_interp=malloc(sizeof(EOBNRHMROMdata_interp));
  else
  {
    EOBNRHMROMdata_interp_Cleanup(*data_interp);
  }
  (*data_interp)->Camp_interp = NULL;
  (*data_interp)->Cphi_interp = NULL;
  (*data_interp)->shifttime_interp = NULL;
  (*data_interp)->shiftphase_interp = NULL;
}
void EOBNRHMROMdata_coeff_Init(EOBNRHMROMdata_coeff **data_coeff) {
  if(!data_coeff) exit(1);
  /* Create storage for structures */
  if(!*data_coeff) *data_coeff=malloc(sizeof(EOBNRHMROMdata_coeff));
  else
  {
    EOBNRHMROMdata_coeff_Cleanup(*data_coeff);
  }
  gsl_set_error_handler(&Err_Handler);
  (*data_coeff)->Camp_coeff = gsl_vector_alloc(nk_amp);
  (*data_coeff)->Cphi_coeff = gsl_vector_alloc(nk_phi);
  (*data_coeff)->shifttime_coeff = malloc(sizeof(double));
  (*data_coeff)->shiftphase_coeff = malloc(sizeof(double));
}
void EOBNRHMROMdata_Cleanup(EOBNRHMROMdata *data /* data to destroy */) {
  if(data->q) gsl_vector_free(data->q);
  if(data->freq) gsl_vector_free(data->freq);
  if(data->Camp) gsl_matrix_free(data->Camp);
  if(data->Cphi) gsl_matrix_free(data->Cphi);
  if(data->Bamp) gsl_matrix_free(data->Bamp);
  if(data->Bphi) gsl_matrix_free(data->Bphi);
  if(data->shifttime) gsl_vector_free(data->shifttime);
  if(data->shiftphase) gsl_vector_free(data->shiftphase);
  free(data);
}
void EOBNRHMROMdata_coeff_Cleanup(EOBNRHMROMdata_coeff *data_coeff) {
  if(data_coeff->Camp_coeff) gsl_vector_free(data_coeff->Camp_coeff);
  if(data_coeff->Cphi_coeff) gsl_vector_free(data_coeff->Cphi_coeff);
  if(data_coeff->shifttime_coeff) free(data_coeff->shifttime_coeff);
  if(data_coeff->shiftphase_coeff) free(data_coeff->shiftphase_coeff);
  free(data_coeff);
}
void EOBNRHMROMdata_interp_Cleanup(EOBNRHMROMdata_interp *data_interp) {
  if(data_interp->Camp_interp) SplineList_Destroy(data_interp->Camp_interp);
  if(data_interp->Cphi_interp) SplineList_Destroy(data_interp->Cphi_interp);
  if(data_interp->shifttime_interp) SplineList_Destroy(data_interp->shifttime_interp);
  if(data_interp->shiftphase_interp) SplineList_Destroy(data_interp->shiftphase_interp);
  free(data_interp);
}

/* Read binary ROM data for frequency vectors, coefficients matrices, basis functions matrices, and shiftvectors in time and phase */ 
int Read_Data_Mode(const char dir[], const int mode[2], EOBNRHMROMdata *data) {
  /* Load binary data for amplitude and phase spline coefficients as computed in Mathematica */
  int ret = SUCCESS;
  char* file_q = malloc(strlen(dir)+64);
  char* file_freq = malloc(strlen(dir)+64);
  char* file_Camp = malloc(strlen(dir)+64);
  char* file_Cphi = malloc(strlen(dir)+64);
  char* file_Bamp = malloc(strlen(dir)+64);
  char* file_Bphi = malloc(strlen(dir)+64);
  char* file_shifttime = malloc(strlen(dir)+64);
  char* file_shiftphase = malloc(strlen(dir)+64);
  sprintf(file_q, "%s", "EOBNRv2HMROM_q.dat"); /* The q vector is the same for all modes */
  sprintf(file_freq, "%s%d%d%s", "EOBNRv2HMROM_freq_", mode[0], mode[1], ".dat"); 
  sprintf(file_Camp, "%s%d%d%s", "EOBNRv2HMROM_Camp_", mode[0], mode[1], ".dat");
  sprintf(file_Cphi, "%s%d%d%s", "EOBNRv2HMROM_Cphi_", mode[0], mode[1], ".dat");
  sprintf(file_Bamp, "%s%d%d%s", "EOBNRv2HMROM_Bamp_", mode[0], mode[1], ".dat");
  sprintf(file_Bphi, "%s%d%d%s", "EOBNRv2HMROM_Bphi_", mode[0], mode[1], ".dat");
  sprintf(file_shifttime, "%s%d%d%s", "EOBNRv2HMROM_shifttime_", mode[0], mode[1], ".dat");
  sprintf(file_shiftphase, "%s%d%d%s", "EOBNRv2HMROM_shiftphase_", mode[0], mode[1], ".dat");
  ret |= Read_Vector(dir, file_q, data->q);
  ret |= Read_Vector(dir, file_freq, data->freq);
  ret |= Read_Matrix(dir, file_Camp, data->Camp);
  ret |= Read_Matrix(dir, file_Cphi, data->Cphi);
  ret |= Read_Matrix(dir, file_Bamp, data->Bamp);
  ret |= Read_Matrix(dir, file_Bphi, data->Bphi);
  ret |= Read_Vector(dir, file_shifttime, data->shifttime);
  ret |= Read_Vector(dir, file_shiftphase, data->shiftphase);
  free(file_q);
  free(file_freq);
  free(file_Camp);
  free(file_Cphi);
  free(file_Bamp);
  free(file_Bphi);
  free(file_shifttime);
  free(file_shiftphase);
  return(ret);
}

/* Function interpolating the data in matrix/vector form produces an interpolated data in the form of SplineLists */
int Interpolate_Spline_Data(const EOBNRHMROMdata *data, EOBNRHMROMdata_interp *data_interp) {

  gsl_set_error_handler(&Err_Handler);
  SplineList* splinelist;
  gsl_spline* spline;
  gsl_interp_accel* accel;
  gsl_vector* matrixline;
  gsl_vector* vector;
  int j;

  /* Interpolating Camp */
  splinelist = data_interp->Camp_interp;
  for (j = 0; j<nk_amp; j++) {
    matrixline = gsl_vector_alloc(nbwf);
    gsl_matrix_get_row(matrixline, data->Camp, j);

    accel = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(gsl_interp_cspline, nbwf);
    gsl_spline_init(spline, gsl_vector_const_ptr(data->q, 0), gsl_vector_const_ptr(matrixline, 0), nbwf);

    splinelist = SplineList_AddElementNoCopy(splinelist, spline,  accel, j);
    gsl_vector_free(matrixline);
  }
  data_interp->Camp_interp = splinelist;

  /* Interpolating Cphi */
  splinelist = data_interp->Cphi_interp;
  for (j = 0; j<nk_phi; j++) {
    matrixline = gsl_vector_alloc(nbwf);
    gsl_matrix_get_row(matrixline, data->Cphi, j);

    accel = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(gsl_interp_cspline, nbwf);
    gsl_spline_init(spline, gsl_vector_const_ptr(data->q, 0), gsl_vector_const_ptr(matrixline, 0), nbwf);

    splinelist = SplineList_AddElementNoCopy(splinelist, spline,  accel, j);
    gsl_vector_free(matrixline);
  }
  data_interp->Cphi_interp = splinelist;

  /* Interpolating shifttime */
  splinelist = data_interp->shifttime_interp;
  vector = data->shifttime;

  accel = gsl_interp_accel_alloc();
  spline = gsl_spline_alloc(gsl_interp_cspline, nbwf);
  gsl_spline_init(spline, gsl_vector_const_ptr(data->q, 0), gsl_vector_const_ptr(vector, 0), nbwf);

  splinelist = SplineList_AddElementNoCopy(NULL, spline,  accel, 0); /* This SplineList has only 1 element */
  data_interp->shifttime_interp = splinelist;

  /* Interpolating shiftphase */
  splinelist = data_interp->shiftphase_interp;
  vector = data->shiftphase;

  accel = gsl_interp_accel_alloc();
  spline = gsl_spline_alloc(gsl_interp_cspline, nbwf);
  gsl_spline_init(spline, gsl_vector_const_ptr(data->q, 0), gsl_vector_const_ptr(vector, 0), nbwf);

  splinelist = SplineList_AddElementNoCopy(NULL, spline,  accel, 0); /* This SplineList has only 1 element */
  data_interp->shiftphase_interp = splinelist;

  return SUCCESS;
}

/* Function taking as input interpolated data in the form of SplineLists
 * evaluates for a given q the projection coefficients and shifts in time and phase
*/
int Evaluate_Spline_Data(const double q, const EOBNRHMROMdata_interp* data_interp, EOBNRHMROMdata_coeff* data_coeff){

  SplineList* splinelist;
  /* Evaluating the vector of projection coefficients for the amplitude */
  for (int j=0; j<nk_amp; j++) {
    splinelist = SplineList_GetElement(data_interp->Camp_interp, j);
    gsl_vector_set(data_coeff->Camp_coeff, j, gsl_spline_eval(splinelist->spline, q, splinelist->accel));
  }
  /* Evaluating the vector of projection coefficients for the phase */
  for (int j=0; j<nk_phi; j++) {
    splinelist = SplineList_GetElement(data_interp->Cphi_interp, j);
    gsl_vector_set(data_coeff->Cphi_coeff, j, gsl_spline_eval(splinelist->spline, q, splinelist->accel));
  }
  /* Evaluating the shift in time */
  splinelist = SplineList_GetElement(data_interp->shifttime_interp, 0); /* This SplineList has only one element */
  *(data_coeff->shifttime_coeff) = gsl_spline_eval(splinelist->spline, q, splinelist->accel);
  /* Evaluating the shift in phase */
  splinelist = SplineList_GetElement(data_interp->shiftphase_interp, 0); /* This SplineList has only one element */
  *(data_coeff->shiftphase_coeff) = gsl_spline_eval(splinelist->spline, q, splinelist->accel);

  return SUCCESS;
}

/*
 * Core function for computing the ROM waveform.
 * Evaluates projection coefficients and shifts in time and phase at desired q.
 * Construct 1D splines for amplitude and phase.
 * Compute strain waveform from amplitude and phase.
*/
int EOBNRv2HMROMCore(
  struct tagListmodesCAmpPhaseFrequencySeries **listhlm,
  int nbmode,
  double tRef,
  double phiRef,
  double fRef,
  double Mtot_sec,
  double q,
  double distance)
{
  int ret = SUCCESS;
  int j;
  /* Check output arrays */
  if(!listhlm) exit(1);
  if(*listhlm)
  {
    printf("Error: (*listhlm) is supposed to be NULL, but got %p\n",(*listhlm));
    exit(1);
  }
  /* Check number of modes */
  if(nbmode<1 || nbmode>nbmodemax) {
    printf("Error: incorrect number of modes");
    exit(1);
  }

  /* Check if the data has been set up */
  if(__EOBNRv2HMROM_setup) {
    printf("Error: the ROM data has not been set up\n");
    exit(1);
  }
  /* Set the global pointers to data */
  ListmodesEOBNRHMROMdata* listdata = *__EOBNRv2HMROM_data;
  ListmodesEOBNRHMROMdata_interp* listdata_interp = *__EOBNRv2HMROM_interp;

  /* Global amplitude prefactor - includes total mass scaling, Fourier scaling, distance scaling, and undoing an additional arbitrary scaling */
  double Mtot_msol = Mtot_sec / MTSUN_SI; /* Mtot_msol and M_ROM in units of solar mass */
  double amp0 = (Mtot_msol/M_ROM) * Mtot_sec * 1.E-16 * 1.E6 * PC_SI / distance;

  /* Highest allowed geometric frequency for the first mode of listmode in the ROM - used for fRef
   * by convention, we use the first mode of listmode (presumably the 22) for phiref */
  ListmodesEOBNRHMROMdata* listdata_ref = ListmodesEOBNRHMROMdata_GetMode(listdata, listmode[0][0], listmode[0][1]);
  EOBNRHMROMdata* data_ref = listdata_ref->data;
  double Mf_ROM_max_ref = gsl_vector_get(data_ref->freq, nbfreq-1);
  /* Convert to geometric units the reference time and frequency */
  double tRef_geom = tRef / Mtot_sec;
  double fRef_geom = fRef * Mtot_sec;

  /* Enforce allowed geometric frequency range for fRef */
  /* In case the user asks for a reference frequency higher than covered by the ROM, we keep it that way as we will just 0-pad the waveform (and do it anyway for some modes) */
  if (fRef_geom > Mf_ROM_max_ref || fRef_geom == 0)
    fRef_geom = Mf_ROM_max_ref; /* If fRef > fhigh or 0 we reset fRef to default value of cutoff frequency for the first mode of the list (presumably the 22 mode) */
  if (0 < fRef_geom && fRef_geom < Mf_ROM_min) {
    printf("Reference frequency Mf_ref=%g is smaller than lowest frequency in ROM Mf=%g. Setting it to the lowest frequency in ROM.\n", fRef_geom, Mf_ROM_min);
    fRef_geom = Mf_ROM_min;
  }

  /* Internal storage for the projection coefficients and shifts in time and phase */
  /* Initialized only once, and reused for the different modes */
  EOBNRHMROMdata_coeff *data_coeff = NULL;
  EOBNRHMROMdata_coeff_Init(&data_coeff);

  /* The phase change imposed by phiref, from the phase of the first mode in the list - to be set in the first step of the loop on the modes */
  double phase_change_ref = 0;
  
  /* Main loop over the modes */
  for(int i=0; i<nbmode; i++ ){
    int l = listmode[i][0];
    int m = listmode[i][1];

    /* Getting the relevant modes in the lists of data */
    ListmodesEOBNRHMROMdata* listdata_mode = ListmodesEOBNRHMROMdata_GetMode(listdata, l, m);
    ListmodesEOBNRHMROMdata_interp* listdata_interp_mode = ListmodesEOBNRHMROMdata_interp_GetMode(listdata_interp, l, m);

    /* Evaluating the projection coefficients and shift in time and phase */
    ret |= Evaluate_Spline_Data(q, listdata_interp_mode->data_interp, data_coeff);

    /* Evaluating the unnormalized amplitude and unshifted phase vectors for the mode */
    /* Notice a change in convention: B matrices are transposed with respect to the B matrices in SEOBNRROM */
    /* amp_pts = Bamp . Camp_coeff */
    /* phi_pts = Bphi . Cphi_coeff */
    gsl_vector* amp_f = gsl_vector_alloc(nbfreq);
    gsl_vector* phi_f = gsl_vector_alloc(nbfreq);
    //clock_t begblas = clock();
    gsl_blas_dgemv(CblasNoTrans, 1.0, listdata_mode->data->Bamp, data_coeff->Camp_coeff, 0.0, amp_f);
    gsl_blas_dgemv(CblasNoTrans, 1.0, listdata_mode->data->Bphi, data_coeff->Cphi_coeff, 0.0, phi_f);
    //clock_t endblas = clock();
    //printf("Mode (%d,%d) Blas time: %g s\n", l, m, (double)(endblas - begblas) / CLOCKS_PER_SEC);   

     /* The downsampled frequencies for the mode - we undo the rescaling of the frequency for the 44 and 55 modes */
    gsl_vector* freq_ds = gsl_vector_alloc(nbfreq);
    gsl_vector_memcpy(freq_ds, listdata_mode->data->freq);
    if ( l==4 && m==4) gsl_vector_scale( freq_ds, 1./Scaling44(q));
    if ( l==5 && m==5) gsl_vector_scale( freq_ds, 1./Scaling55(q));

    /* Evaluating the shifts in time and phase - conditional scaling for the 44 and 55 modes */
    SplineList* shifttime_splinelist = listdata_interp_mode->data_interp->shifttime_interp;
    SplineList* shiftphase_splinelist = listdata_interp_mode->data_interp->shiftphase_interp;
    double shifttime;
    if( l==4 && m==4) {
      shifttime = gsl_spline_eval(shifttime_splinelist->spline, q, shifttime_splinelist->accel) * Scaling44(q);
    }
    else if( l==5 && m==5) {
      shifttime = gsl_spline_eval(shifttime_splinelist->spline, q, shifttime_splinelist->accel) * Scaling55(q);
    }
    else {
      shifttime = gsl_spline_eval(shifttime_splinelist->spline, q, shifttime_splinelist->accel);
    }
    double shiftphase = gsl_spline_eval(shiftphase_splinelist->spline, q, shiftphase_splinelist->accel);

    /* Adding to the internal shifttime the time shift (forward for a + sign) specified by the user with tRef - put in geometric units first */
    shifttime = shifttime + tRef_geom;

    /* If first mode in the list, set the value of phase_change_ref */
    if( i==0 ) {
      /* Setup 1d cubic spline in phase and determine phase_change_ref */
      gsl_interp_accel* accel_phi = gsl_interp_accel_alloc();
      gsl_spline* spline_phi = gsl_spline_alloc(gsl_interp_cspline, nbfreq);
      gsl_spline_init(spline_phi, gsl_vector_const_ptr(freq_ds,0), gsl_vector_const_ptr(phi_f,0), nbfreq);
      phase_change_ref = phiRef + (gsl_spline_eval(spline_phi, fRef_geom, accel_phi) - shifttime * fRef_geom - shiftphase);
      gsl_spline_free(spline_phi);
      gsl_interp_accel_free(accel_phi);
    }
    /* Propagation to this mode of the ref phase change with the proper factor of m
     * and sum with shiftphase */
    double constphaseshift = (double) m/listmode[0][1] * phase_change_ref + shiftphase;

    /* Initialize the complex series for the mode */
    CAmpPhaseFrequencySeries *modefreqseries = NULL;
    int len = (int) freq_ds->size;
    CAmpPhaseFrequencySeries_Init(&modefreqseries, len);

    /* Mode-dependent complete amplitude prefactor */
    double amp_pre = amp0 * ModeAmpFactor( l, m, q);

    /* Final result for the mode */
    /* Scale and set the amplitudes (amplitudes are real at this stage)*/
    gsl_vector_scale(amp_f, amp_pre);
    gsl_vector_memcpy(modefreqseries->amp_real, amp_f);
    gsl_vector_set_zero(modefreqseries->amp_imag); /* Amplitudes are real at this stage */
    /* Add the linear term and the constant (including the shift to phiRef), and set the phases */
    gsl_vector_scale(phi_f, -1.); /* Change the sign of the phases: ROM convention Psi=-phase */
    gsl_blas_daxpy(shifttime, freq_ds, phi_f); /*Beware: here freq_ds must still be in geometric units*/
    gsl_vector_add_constant(phi_f, constphaseshift);
    gsl_vector_memcpy(modefreqseries->phase, phi_f);
    /* Scale (to physical units) and set the frequencies */
    gsl_vector_scale(freq_ds, 1./Mtot_sec);
    gsl_vector_memcpy(modefreqseries->freq, freq_ds);

    /* Append the computed mode to the ListmodesCAmpPhaseFrequencySeries structure */
    *listhlm = ListmodesCAmpPhaseFrequencySeries_AddModeNoCopy(*listhlm, modefreqseries, l, m);

    /* Cleanup for the mode */
    gsl_vector_free(freq_ds);
    gsl_vector_free(amp_f);
    gsl_vector_free(phi_f);
  }

  /* Cleanup of the coefficients data structure */
  EOBNRHMROMdata_coeff_Cleanup(data_coeff);

  return(SUCCESS);
}

/* Compute waveform in downsampled frequency-amplitude-phase format */
int SimEOBNRv2HMROM(
  struct tagListmodesCAmpPhaseFrequencySeries **listhlm,  /* Output: list of modes in Frequency-domain amplitude and phase form */
  int nbmode,                                    /* Number of modes to generate (starting with the 22) */
  double tRef,                                   /* Time shift with respect to the 22-fit removed waveform (s) */
  double phiRef,                                 /* Phase at reference frequency */
  double fRef,                                   /* Reference frequency (Hz); 0 defaults to fLow */
  double m1SI,                                   /* Mass of companion 1 (kg) */
  double m2SI,                                   /* Mass of companion 2 (kg) */
  double distance)                               /* Distance of source (m) */
{
  /* Get masses in terms of solar mass */
  double mass1 = m1SI / MSUN_SI;
  double mass2 = m2SI / MSUN_SI;
  double Mtot = mass1 + mass2;
  double q = fmax(mass1/mass2, mass2/mass1);    /* Mass-ratio >1 by convention*/
  double Mtot_sec = Mtot * MTSUN_SI; /* Total mass in seconds */

  if ( q > q_max ) {
    printf( "Error - %s: q out of range!\nEOBNRv2HMROM is only available for a mass ratio in the range q <= %g.\n", __func__, q_max);
    exit(1);
  }

  /* Set up (load and build interpolation) ROM data if not setup already */
  //clock_t beg = clock();
  EOBNRv2HMROM_Init_DATA();
  //clock_t end = clock();
  //printf("Initialization time: %g s\n", (double)(end - beg) / CLOCKS_PER_SEC);
  
  //beg = clock();
  int retcode = EOBNRv2HMROMCore(listhlm, nbmode, tRef, phiRef, fRef, Mtot_sec, q, distance);
  //end = clock();
  //printf("ROM evaluation time: %g s\n", (double)(end - beg) / CLOCKS_PER_SEC);

  return(retcode);
}

/* Setup EOBNRv2HMROM model using data files installed in $ROM_DATA_PATH */
int EOBNRv2HMROM_Init_DATA(void) {
  if (!__EOBNRv2HMROM_setup) return SUCCESS;

  int ret=FAILURE;
  char *envpath=NULL;
  char path[32768];
  char *brkt,*word;
  envpath=getenv("ROM_DATA_PATH");
  if(!envpath) {
    printf("Error: the environment variable ROM_DATA_PATH, giving the path to the ROM data, seems undefined\n");
    return(FAILURE);
  }
  strncpy(path,envpath,sizeof(path));

  for(word=strtok_r(path,":",&brkt); word; word=strtok_r(NULL,":",&brkt))
  {
    ret = EOBNRv2HMROM_Init(word);
    if(ret == SUCCESS) break;
  }
  if(ret!=SUCCESS) {
    printf("Error: unable to find EOBNRv2HMROM data files in $ROM_DATA_PATH\n");
    exit(FAILURE);
  }
  __EOBNRv2HMROM_setup = ret;
  return(ret);
}

/* Setup EOBNRv2HMROM model using data files installed in dir */
int EOBNRv2HMROM_Init(const char dir[]) {
  if(!__EOBNRv2HMROM_setup) {
    printf("Error: EOBNRHMROMdata was already set up!");
    exit(1);
  }

  int ret = SUCCESS;
  ListmodesEOBNRHMROMdata* listdata = *__EOBNRv2HMROM_data;
  ListmodesEOBNRHMROMdata_interp* listdata_interp = *__EOBNRv2HMROM_interp;

  for (int j=0; j<nbmodemax; j++) { /* At setup, we initialize all available modes anyway */

    EOBNRHMROMdata* data = NULL;
    EOBNRHMROMdata_Init(&data);
    ret |= Read_Data_Mode( dir, listmode[j], data);
    if (!ret) {
      listdata = ListmodesEOBNRHMROMdata_AddModeNoCopy( listdata, data, listmode[j][0], listmode[j][1]);

      EOBNRHMROMdata_interp* data_interp = NULL;
      EOBNRHMROMdata_interp_Init(&data_interp);
      ret |= Interpolate_Spline_Data(data, data_interp);
      if (!ret) listdata_interp = ListmodesEOBNRHMROMdata_interp_AddModeNoCopy( listdata_interp, data_interp, listmode[j][0], listmode[j][1]);
    }
  }

  __EOBNRv2HMROM_setup = ret;
  if (!ret) {
    *__EOBNRv2HMROM_data = listdata;
    *__EOBNRv2HMROM_interp = listdata_interp;
  }
  return(ret);
}
