/**
 * \author Sylvain Marsat, University of Maryland - NASA GSFC
 *
 * \brief C code for for functions manipulating waveforms.
 *
 */

#define _XOPEN_SOURCE 500

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#include "waveform.h"

/* NOTE: uses the list of modes of EOBNRv2HMROM (listmode), to be extended when more waveform models are added */

/***************** Function estimating time to coalescence ****************/

/* Functions reading from a list of modes the minimal and maximal frequencies */
double ListmodesCAmpPhaseFrequencySeries_maxf(ListmodesCAmpPhaseFrequencySeries* listhlm)
{
  double maxf = 0;
  /* Main loop over the modes - goes through all the modes present */
  ListmodesCAmpPhaseFrequencySeries* listelementhlm = listhlm;
  CAmpPhaseFrequencySeries* hlm;
  while(listelementhlm) {
    hlm = listelementhlm->freqseries;
    maxf = fmax(maxf, gsl_vector_get(hlm->freq, hlm->freq->size - 1));
    listelementhlm = listelementhlm->next;
  }
  return maxf;
}
double ListmodesCAmpPhaseFrequencySeries_minf(ListmodesCAmpPhaseFrequencySeries* listhlm)
{
  double minf = 0;
  /* Main loop over the modes - goes through all the modes present */
  ListmodesCAmpPhaseFrequencySeries* listelementhlm = listhlm;
  CAmpPhaseFrequencySeries* hlm;
  while(listelementhlm) {
    hlm = listelementhlm->freqseries;
    minf = fmin(minf, gsl_vector_get(hlm->freq, hlm->freq->size - 1));
    listelementhlm = listelementhlm->next;
  }
  return minf;
}

/* Function estimating initial time from Psi22, according to tf_SPA = -1/(2pi)dPsi/df */
double EstimateInitialTime(ListmodesCAmpPhaseFrequencySeries* listhlm, double fLow)
{
  CAmpPhaseFrequencySeries* h22 = ListmodesCAmpPhaseFrequencySeries_GetMode(listhlm, 2, 2)->freqseries;
  double* freq = h22->freq->data;
  double* phi22 = h22->phase->data; /* Note that the waveform is given with phi = -Psi */
  int i = 0;
  while(freq[i]<fLow) i++;
  return 1./(2*PI) * (phi22[i+1] - phi22[i])/(freq[i+1] - freq[i]);
}

/***************** Functions to manipulate ReImFrequencySeries structure ****************/

/* Function evaluating a CAmpPhaseFrequencySeries on a given set of frequencies, and adding it to a ReImFrequencySeries */
void ReImFrequencySeries_AddCAmpPhaseFrequencySeries(
  struct tagReImFrequencySeries* freqseriesReIm,              /* Output Re/Im frequency series */
  struct tagCAmpPhaseFrequencySeries* freqseriesCAmpPhase,    /* Input CAmp/Phase frequency series, to be interpolated and added to the output */
  double fstartobsmode)                                       /* Starting frequency in case of limited duration of observations- assumed to have been scaled with the proper factor m/2 for this mode - set to 0 to ignore */
{
  /* Frequency vectors, and sizes */
  gsl_vector* freqin = freqseriesCAmpPhase->freq;
  gsl_vector* freqout = freqseriesReIm->freq;
  int sizein = (int) freqin->size;
  int sizeout = (int) freqout->size;
  /* Input, Output vectors */
  gsl_vector* vecampreal = freqseriesCAmpPhase->amp_real;
  gsl_vector* vecampimag = freqseriesCAmpPhase->amp_imag;
  gsl_vector* vecphase = freqseriesCAmpPhase->phase;
  gsl_vector* vechreal = freqseriesReIm->h_real;
  gsl_vector* vechimag = freqseriesReIm->h_imag;

  /* Initializing the splines */
  /* Note: since this must also apply to mode contribution after processing, real and imaginary parts of the amplitude are present - but since they should differ for LLV detectors by a constant factor, they can be interpolated */
  gsl_interp_accel* accel_ampreal = gsl_interp_accel_alloc();
  gsl_interp_accel* accel_ampimag = gsl_interp_accel_alloc();
  gsl_interp_accel* accel_phase = gsl_interp_accel_alloc();
  gsl_spline* ampreal = gsl_spline_alloc(gsl_interp_cspline, sizein);
  gsl_spline* ampimag = gsl_spline_alloc(gsl_interp_cspline, sizein);
  gsl_spline* phase = gsl_spline_alloc(gsl_interp_cspline, sizein);
  gsl_spline_init(ampreal, gsl_vector_const_ptr(freqin,0), gsl_vector_const_ptr(vecampreal,0), sizein);
  gsl_spline_init(ampimag, gsl_vector_const_ptr(freqin,0), gsl_vector_const_ptr(vecampimag,0), sizein);
  gsl_spline_init(phase, gsl_vector_const_ptr(freqin,0), gsl_vector_const_ptr(vecphase,0), sizein);

  /* First and last index of the output frequency vector that are covered by the CAmpPhase data */
  int jStart = 0;
  int jStop = sizeout - 1;
  double minfmode = fmax(gsl_vector_get(freqin, 0), fstartobsmode); /* Takes into account fstartobsmode here */
  double maxfmode = gsl_vector_get(freqin, sizein - 1);
  while(gsl_vector_get(freqout, jStart) < minfmode && jStart < sizeout-1) jStart++;
  while(gsl_vector_get(freqout, jStop) > maxfmode && jStop > 0) jStop--;
  if(jStop <= jStart) {
    printf("Error: empty range of frequencies in ReImFrequencySeries_AddCAmpPhaseFrequencySeries.\n");
    exit(1);
  }

  /* Main loop - evaluating the interpolating splines and adding to the output data */
  double f, phi;
  double complex A;
  double complex h;
  double* freqoutdata = freqout->data;
  double* hrealdata = vechreal->data;
  double* himagdata = vechimag->data;
  for(int j=jStart; j<=jStop; j++) { /* Note: loop to jStop included */
    f = freqoutdata[j];
    A = gsl_spline_eval(ampreal, f, accel_ampreal) + I*gsl_spline_eval(ampimag, f, accel_ampimag);
    phi = gsl_spline_eval(phase, f, accel_phase);
    h = A * cexp(I*phi);
    hrealdata[j] += creal(h);
    himagdata[j] += cimag(h);
  }

  /* Clean up */
  gsl_interp_accel_free(accel_ampreal);
  gsl_interp_accel_free(accel_ampimag);
  gsl_interp_accel_free(accel_phase);
  gsl_spline_free(ampreal);
  gsl_spline_free(ampimag);
  gsl_spline_free(phase);
}

/* Function evaluating a ReImFrequencySeries by interpolating wach mode of a ListmodesCAmpPhaseFrequencySeries and summing them, given a set of frequencies */
void ReImFrequencySeries_SumListmodesCAmpPhaseFrequencySeries(
  struct tagReImFrequencySeries* freqseriesReIm,                    /* Output Re/Im frequency series - already initialized */
  struct tagListmodesCAmpPhaseFrequencySeries* listmodesCAmpPhase,  /* Input CAmp/Phase frequency series, to be interpolated */
  gsl_vector* freq,                                                 /* Input set of frequencies on which evaluating */
  double fstartobs)                                                 /* Starting frequency in case of limited duration of observations - set to 0 to ignore */
{
  /* Check the sizes */
  if(freq->size != freqseriesReIm->freq->size) {
    printf("Error: incompatible sizes in ReImFrequencySeries_SumListmodesCAmpPhaseFrequencySeries.\n");
    exit(1);
  }

  /* Initialize and copy frequencies */
  gsl_vector_set_zero(freqseriesReIm->h_real);
  gsl_vector_set_zero(freqseriesReIm->h_imag);
  gsl_vector_memcpy(freqseriesReIm->freq, freq);

  /* Main loop: go through the list of modes, interpolate and add them to the output */
  ListmodesCAmpPhaseFrequencySeries* listelement = listmodesCAmpPhase;
  double fstartobsmode;
  while(listelement) {
    double fstartobsmode = fmax(fstartobs, ((double) listelement->m)/2. * fstartobs);
    ReImFrequencySeries_AddCAmpPhaseFrequencySeries(freqseriesReIm, listelement->freqseries, fstartobsmode);
    listelement = listelement->next;
  }
}

/* Helper function to add a mode to hplus, hcross in Fourier domain
 * - copies the function XLALSimAddMode, which was done only for TD structures */
int FDAddMode(
  ReImFrequencySeries* hptilde,       /* Output: frequency series for hplus */
  ReImFrequencySeries* hctilde,       /* Output: frequency series for hcross */
  ReImFrequencySeries* hlmtilde,      /* Input: frequency series for the mode hlm */
  double Theta,                       /* First angle for position in the sky of observer */
  double Phi,                         /* Second angle for position in the sky of observer */
  int l,                              /* First mode number l */
  int m,                              /* Second mode number m */
  int sym)                            /* If 1, assume planar symmetry and add also mode l,-m. Do not if set to 0. */
{
  /* Deleted the definition of the string 'func': usage ? */
  double complex Y;
  int j;
  double complex hlmtildevalue;

  double* datap_real = hptilde->h_real->data;
  double* datap_imag = hptilde->h_imag->data;
  double* datac_real = hctilde->h_real->data;
  double* datac_imag = hctilde->h_imag->data;

  int minus1l; /* (-1)^l */
  if ( l%2 ) minus1l = -1;
  else minus1l = 1;
  if ( sym ) { /* equatorial symmetry: add in -m mode */
    Y = SpinWeightedSphericalHarmonic(Theta, Phi, -2, l, m);
    double complex Ymstar = conj(SpinWeightedSphericalHarmonic(Theta, Phi, -2, l, -m));
    double complex factorp = 1./2*(Y + minus1l*Ymstar);
    double complex factorc = I/2*(Y - minus1l*Ymstar);
    for ( j = 0; j < hlmtilde->freq->size; ++j ) {
      hlmtildevalue = hlmtilde->h_real->data[j] + I*hlmtilde->h_imag->data[j];
      datap_real[j] += creal(factorp*hlmtildevalue);
      datap_imag[j] += cimag(factorp*hlmtildevalue);
      datac_real[j] += creal(factorc*hlmtildevalue);
      datac_imag[j] += cimag(factorc*hlmtildevalue);
    }
  }
  else { /* not adding in the -m mode */
    Y = SpinWeightedSphericalHarmonic(Theta, Phi, -2, l, m);
    double complex factorp = 1./2*Y;
    double complex factorc = I/2*Y;
    for ( j = 0; j < hlmtilde->freq->size; ++j ) {
      hlmtildevalue = hlmtilde->h_real->data[j] + I*hlmtilde->h_imag->data[j];
      datap_real[j] += creal(factorp*hlmtildevalue);
      datap_imag[j] += cimag(factorp*hlmtildevalue);
      datac_real[j] += creal(factorc*hlmtildevalue);
      datac_imag[j] += cimag(factorc*hlmtildevalue);
    }
  }

  return 0;
}

/* Function evaluating the FD frequency series for hplus, hcross from the modes hlm */
int GeneratehphcFDReImFrequencySeries(
  ReImFrequencySeries** hptilde,                 /* Output: frequency series for hplus */
  ReImFrequencySeries** hctilde,                 /* Output: frequency series for hcross */
  ListmodesCAmpPhaseFrequencySeries* listhlm,    /* Input: frequency series for the mode hlm  */
  double fLow,                                   /* Minimal frequency */
  double deltaf,                                 /* Frequency step */
  int nbpt,                                      /* Number of points of output - if 0, determined from deltaf and maximal frequency in input */
  int nbmode,                                    /* Number of modes to add */
  double Theta,                                  /* First angle for position in the sky of observer */
  double Phi,                                    /* Second angle for position in the sky of observer */
  int sym)                                       /* If 1, assume planar symmetry and add also mode l,-m. Do not if set to 0. */
{
  ListmodesCAmpPhaseFrequencySeries* hlm;

  /* Determine total size of frequency series if nbpt is not set */
  int nbpts;
  if(nbpt==0) {
    double fHigh = 0;
    for(int i=0; i<nbmode; i++) {
      int l = listmode[i][0];
      int m = listmode[i][1];
      hlm = ListmodesCAmpPhaseFrequencySeries_GetMode(listhlm, l, m);
      fHigh = fmax(fHigh, gsl_vector_get(hlm->freqseries->freq, hlm->freqseries->freq->size - 1));
    }
    printf("fHigh: %g\n", fHigh);
    printf("deltaf: %g\n", deltaf);
    nbpts = (int) ceil(fHigh/deltaf);
  }
  else nbpts = nbpt;

  /* Initialize hplus, hcross */
  ReImFrequencySeries_Init(hptilde, nbpts);
  ReImFrequencySeries_Init(hctilde, nbpts);
  gsl_vector_set_zero((*hptilde)->h_real);
  gsl_vector_set_zero((*hptilde)->h_imag);
  gsl_vector_set_zero((*hctilde)->h_real);
  gsl_vector_set_zero((*hctilde)->h_imag);


  /* Set values for frequency vectors */
  double f;
  for(int i=0; i<nbpts; i++) {
    f = i*deltaf;
    gsl_vector_set((*hptilde)->freq, i, f);
    gsl_vector_set((*hctilde)->freq, i, f);
  }

  /* Main loop over the modes */
  for(int i=0; i<nbmode; i++) {
    int l = listmode[i][0];
    int m = listmode[i][1];
    
    /* Initialize the complex series for the mode */
    ReImFrequencySeries* mode = NULL;
    ReImFrequencySeries_Init(&mode, nbpts);
    gsl_vector_set_zero(mode->h_real);
    gsl_vector_set_zero(mode->h_imag);

    /* Get Amp/Phase frequency series for the mode */
    ListmodesCAmpPhaseFrequencySeries* hlm = ListmodesCAmpPhaseFrequencySeries_GetMode(listhlm, l, m);
    gsl_vector* freq_ds = hlm->freqseries->freq;
    gsl_vector* amp_real = hlm->freqseries->amp_real;
    gsl_vector* amp_imag = hlm->freqseries->amp_imag;
    gsl_vector* phase = hlm->freqseries->phase;
    int nbfreq = (int) freq_ds->size;

    /* Setup 1d cubic spline for the phase and amplitude of the mode */
    gsl_interp_accel* accel_phase = gsl_interp_accel_alloc();
    gsl_interp_accel* accel_amp_real = gsl_interp_accel_alloc();
    gsl_interp_accel* accel_amp_imag = gsl_interp_accel_alloc();
    gsl_spline* spline_phase = gsl_spline_alloc(gsl_interp_cspline, nbfreq);
    gsl_spline* spline_amp_real = gsl_spline_alloc(gsl_interp_cspline, nbfreq);
    gsl_spline* spline_amp_imag = gsl_spline_alloc(gsl_interp_cspline, nbfreq);
    gsl_spline_init(spline_phase, gsl_vector_const_ptr(freq_ds,0), gsl_vector_const_ptr(phase,0), nbfreq);
    gsl_spline_init(spline_amp_real, gsl_vector_const_ptr(freq_ds,0), gsl_vector_const_ptr(amp_real,0), nbfreq);
    gsl_spline_init(spline_amp_imag, gsl_vector_const_ptr(freq_ds,0), gsl_vector_const_ptr(amp_imag,0), nbfreq);

    /* Interval in frequency covered by the ROM */
    double fLow_mode = fmax(fLow, gsl_vector_get(freq_ds, 0));
    double fHigh_mode = fmin(deltaf*(nbpts-1), gsl_vector_get(freq_ds, nbfreq-1));
    /* Initialize the loop - values outside this range in j are 0 by default */
    int jStart = ceil(fLow_mode / deltaf);
    int jStop = ceil(fHigh_mode / deltaf);
    double* modedata_real = mode->h_real->data;
    double* modedata_imag = mode->h_imag->data;
    /* Loop on the frequency samples chosen to evaluate the waveform */
    /* We set apart the first and last step to avoid falling outside of the range of the splines by numerical errors */
    double f, Ar, Ai, phi;
    double complex h;

    f = fmax(fLow_mode, jStart*deltaf);
    Ar = gsl_spline_eval(spline_amp_real, f, accel_amp_real);
    Ai = gsl_spline_eval(spline_amp_imag, f, accel_amp_imag);
    phi = gsl_spline_eval(spline_phase, f, accel_phase);
    h = (Ar + I*Ai) * cexp(I*phi);
    modedata_real[jStart] = creal(h);
    modedata_imag[jStart] = cimag(h);

    for(int j=jStart+1; j<jStop-1; j++) {
      f = j*deltaf;
      Ar = gsl_spline_eval(spline_amp_real, f, accel_amp_real);
      Ai = gsl_spline_eval(spline_amp_imag, f, accel_amp_imag);
      phi = gsl_spline_eval(spline_phase, f, accel_phase);
      h = (Ar + I*Ai) * cexp(I*phi);
      modedata_real[j] = creal(h);
      modedata_imag[j] = cimag(h);
    }

    f = fmin(fHigh_mode, (jStop-1)*deltaf);
    Ar = gsl_spline_eval(spline_amp_real, f, accel_amp_real);
    Ai = gsl_spline_eval(spline_amp_imag, f, accel_amp_imag);
    phi = gsl_spline_eval(spline_phase, f, accel_phase);
    h = (Ar + I*Ai) * cexp(I*phi);
    modedata_real[jStop-1] = creal(h);
    modedata_imag[jStop-1] = cimag(h);

    /* Add the computed mode to the SphHarmFrequencySeries structure */
    int sym;
    if ( m==0 ) sym = 0; /* We test for hypothetical m=0 modes */
    else sym = 1;

    FDAddMode(*hptilde, *hctilde, mode, Theta, Phi, l, m, sym);

    /* Cleanup for the mode */
    gsl_spline_free(spline_amp_real);
    gsl_spline_free(spline_amp_imag);
    gsl_spline_free(spline_phase);
    gsl_interp_accel_free(accel_amp_real);
    gsl_interp_accel_free(accel_amp_imag);
    gsl_interp_accel_free(accel_phase);
    ReImFrequencySeries_Cleanup(mode);
  }
}

/* Function evaluating the FD frequency series by summing mode contributions from each hlm */
/* NOTE: nbmode as arg here is a bit deceiptive, only used for max frequency */
int GenerateFDReImFrequencySeries(
  ReImFrequencySeries** freqseries,              /* Output: frequency series */
  ListmodesCAmpPhaseFrequencySeries* listhlm,    /* Input: FD modes hlm in the form AmpReal/AmpIm/Phase  */
  double fLow,                                   /* Minimal frequency */
  double deltaf,                                 /* Frequency step */
  int nbpt,                                      /* Number of points of output - if 0, determined from deltaf and maximal frequency in input */
  int nbmode)                                    /* Number of modes - used only to determine the highest frequency */
{
  ListmodesCAmpPhaseFrequencySeries* hlm;

  /* Determine total size of frequency series if nbpt is not set */
  int nbpts;
  if(nbpt==0) {
    double fHigh = 0;
    for(int i=0; i<nbmode; i++) {
      int l = listmode[i][0];
      int m = listmode[i][1];
      hlm = ListmodesCAmpPhaseFrequencySeries_GetMode(listhlm, l, m);
      fHigh = fmax(fHigh, gsl_vector_get(hlm->freqseries->freq, hlm->freqseries->freq->size - 1));
    }
    nbpts = (int) ceil(fHigh/deltaf);
  }
  else nbpts = nbpt;

  /* Initialize frequency series */
  ReImFrequencySeries_Init(freqseries, nbpts);
  gsl_vector_set_zero((*freqseries)->h_real);
  gsl_vector_set_zero((*freqseries)->h_imag);

  /* Set values for frequency vector */
  gsl_vector* freq = gsl_vector_alloc(nbpts);
  double f;
  for(int i=0; i<nbpts; i++) {
    f = i*deltaf;
    gsl_vector_set(freq, i, f);
  }

  /* Sum mode contributions */
  /* NOTE: fLow is used as fstartobs here, different meaning in spirit in the rest of the code */
  ReImFrequencySeries_SumListmodesCAmpPhaseFrequencySeries((*freqseries), listhlm, freq, fLow);

  /* Cleanup */
  gsl_vector_free(freq);
}

/* Function to restrict a frequency series (typically output of a FFT) to a given frequency range */
int RestrictFDReImFrequencySeries(
  ReImFrequencySeries** freqseriesout,           /* Output: truncated frequency series */
  ReImFrequencySeries* freqseriesin,             /* Input: frequency series */
  double fLow,                                   /* Minimal frequency */
  double fHigh)                                   /* Maximal frequency */
{
  /* Determine first and last indices */
  int imin = 0;
  int imax = freqseriesin->freq->size - 1;
  double* freqin = freqseriesin->freq->data;
  double* hrealin = freqseriesin->h_real->data;
  double* himagin = freqseriesin->h_imag->data;
  while(freqin[imin]<fLow) imin++;
  while(freqin[imax]>fHigh) imax--;

  /* Initialize output */
  int size = imax - imin + 1;
  ReImFrequencySeries_Init(freqseriesout, size);

  /* Write output */
  double* freqout = (*freqseriesout)->freq->data;
  double* hrealout = (*freqseriesout)->h_real->data;
  double* himagout = (*freqseriesout)->h_imag->data;
  for(int i=0; i<size; i++) {
    freqout[i] = freqin[imin+i];
    hrealout[i] = hrealin[imin+i];
    himagout[i] = himagin[imin+i];
  }

  return SUCCESS;
}

/***************** Other structure functions ****************/
/* Function reproducing XLALSpinWeightedSphericalHarmonic
 * - Currently only supports s=-2, l=2,3,4,5 modes */
double complex SpinWeightedSphericalHarmonic(double theta, double phi, int s, int l, int m)
{
  static const char *func = "SpinWeightedSphericalHarmonic";
  double fac;
  double complex ans;
 
  /* sanity checks ... */
  if ( l < abs(s) ) {
    printf("Error - %s: Invalid mode s=%d, l=%d, m=%d - require |s| <= l\n", func, s, l, m );
    exit(1);
  }
  if ( l < abs(m) ) {
    printf("Error - %s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", func, s, l, m );
    exit(1);
  }
 
  if ( s == -2 ) {
    if ( l == 2 ) {
      switch ( m ) {
      case -2:
	fac = sqrt( 5.0 / ( 64.0 * PI ) ) * ( 1.0 - cos( theta ))*( 1.0 - cos( theta ));
	break;
      case -1:
	fac = sqrt( 5.0 / ( 16.0 * PI ) ) * sin( theta )*( 1.0 - cos( theta ));
	break;
         
      case 0:
	fac = sqrt( 15.0 / ( 32.0 * PI ) ) * sin( theta )*sin( theta );
	break;
         
      case 1:
	fac = sqrt( 5.0 / ( 16.0 * PI ) ) * sin( theta )*( 1.0 + cos( theta ));
	break;
         
      case 2:
	fac = sqrt( 5.0 / ( 64.0 * PI ) ) * ( 1.0 + cos( theta ))*( 1.0 + cos( theta ));
	break;     
      default:
	printf("Error - %s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", func, s, l, m );
	exit(1);
	break;
      } /*  switch (m) */
    }  /* l==2*/
    else if ( l == 3 ) {
      switch ( m ) {
      case -3:
	fac = sqrt(21.0/(2.0*PI))*cos(theta/2.0)*pow(sin(theta/2.0),5.0);
	break;
      case -2:
	fac = sqrt(7.0/4.0*PI)*(2.0 + 3.0*cos(theta))*pow(sin(theta/2.0),4.0);
	break;
      case -1:
	fac = sqrt(35.0/(2.0*PI))*(sin(theta) + 4.0*sin(2.0*theta) - 3.0*sin(3.0*theta))/32.0;
	break;
      case 0:
	fac = (sqrt(105.0/(2.0*PI))*cos(theta)*pow(sin(theta),2.0))/4.0;
	break;
      case 1:
	fac = -sqrt(35.0/(2.0*PI))*(sin(theta) - 4.0*sin(2.0*theta) - 3.0*sin(3.0*theta))/32.0;
	break;
           
      case 2:
	fac = sqrt(7.0/PI)*pow(cos(theta/2.0),4.0)*(-2.0 + 3.0*cos(theta))/2.0;
	break;     
           
      case 3:
	fac = -sqrt(21.0/(2.0*PI))*pow(cos(theta/2.0),5.0)*sin(theta/2.0);
	break;     
           
      default:
	printf("Error - %s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", func, s, l, m );
	exit(1);
	break;
      } 
    }   /* l==3 */ 
    else if ( l == 4 ) {
      switch ( m ) {
      case -4:
	fac = 3.0*sqrt(7.0/PI)*pow(cos(theta/2.0),2.0)*pow(sin(theta/2.0),6.0);
	break;
      case -3:
	fac = 3.0*sqrt(7.0/(2.0*PI))*cos(theta/2.0)*(1.0 + 2.0*cos(theta))*pow(sin(theta/2.0),5.0);
	break;
         
      case -2:
	fac = (3.0*(9.0 + 14.0*cos(theta) + 7.0*cos(2.0*theta))*pow(sin(theta/2.0),4.0))/(4.0*sqrt(PI));
	break;
      case -1:
	fac = (3.0*(3.0*sin(theta) + 2.0*sin(2.0*theta) + 7.0*sin(3.0*theta) - 7.0*sin(4.0*theta)))/(32.0*sqrt(2.0*PI));
	break;
      case 0:
	fac = (3.0*sqrt(5.0/(2.0*PI))*(5.0 + 7.0*cos(2.0*theta))*pow(sin(theta),2.0))/16.0;
	break;
      case 1:
	fac = (3.0*(3.0*sin(theta) - 2.0*sin(2.0*theta) + 7.0*sin(3.0*theta) + 7.0*sin(4.0*theta)))/(32.0*sqrt(2.0*PI));
	break;
      case 2:
	fac = (3.0*pow(cos(theta/2.0),4.0)*(9.0 - 14.0*cos(theta) + 7.0*cos(2.0*theta)))/(4.0*sqrt(PI));
	break;     
      case 3:
	fac = -3.0*sqrt(7.0/(2.0*PI))*pow(cos(theta/2.0),5.0)*(-1.0 + 2.0*cos(theta))*sin(theta/2.0);
	break;     
      case 4:
	fac = 3.0*sqrt(7.0/PI)*pow(cos(theta/2.0),6.0)*pow(sin(theta/2.0),2.0);
	break;     
      default:
	printf("Error - %s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", func, s, l, m );
	exit(1);
	break;
      }
    }    /* l==4 */
    else if ( l == 5 ) {
      switch ( m ) {
      case -5:
	fac = sqrt(330.0/PI)*pow(cos(theta/2.0),3.0)*pow(sin(theta/2.0),7.0);
	break;
      case -4:
	fac = sqrt(33.0/PI)*pow(cos(theta/2.0),2.0)*(2.0 + 5.0*cos(theta))*pow(sin(theta/2.0),6.0);
	break;
      case -3:
	fac = (sqrt(33.0/(2.0*PI))*cos(theta/2.0)*(17.0 + 24.0*cos(theta) + 15.0*cos(2.0*theta))*pow(sin(theta/2.0),5.0))/4.0;
	break;
      case -2:
	fac = (sqrt(11.0/PI)*(32.0 + 57.0*cos(theta) + 36.0*cos(2.0*theta) + 15.0*cos(3.0*theta))*pow(sin(theta/2.0),4.0))/8.0;
	break;
      case -1:
	fac = (sqrt(77.0/PI)*(2.0*sin(theta) + 8.0*sin(2.0*theta) + 3.0*sin(3.0*theta) + 12.0*sin(4.0*theta) - 15.0*sin(5.0*theta)))/256.0;
	break;
      case 0:
	fac = (sqrt(1155.0/(2.0*PI))*(5.0*cos(theta) + 3.0*cos(3.0*theta))*pow(sin(theta),2.0))/32.0;
	break;
      case 1:
	fac = sqrt(77.0/PI)*(-2.0*sin(theta) + 8.0*sin(2.0*theta) - 3.0*sin(3.0*theta) + 12.0*sin(4.0*theta) + 15.0*sin(5.0*theta))/256.0;
	break;
      case 2:
	fac = sqrt(11.0/PI)*pow(cos(theta/2.0),4.0)*(-32.0 + 57.0*cos(theta) - 36.0*cos(2.0*theta) + 15.0*cos(3.0*theta))/8.0;
	break;     
      case 3:
	fac = -sqrt(33.0/(2.0*PI))*pow(cos(theta/2.0),5.0)*(17.0 - 24.0*cos(theta) + 15.0*cos(2.0*theta))*sin(theta/2.0)/4.0;
	break;     
      case 4:
	fac = sqrt(33.0/PI)*pow(cos(theta/2.0),6.0)*(-2.0 + 5.0*cos(theta))*pow(sin(theta/2.0),2.0);
	break;     
      case 5:
	fac = -sqrt(330.0/PI)*pow(cos(theta/2.0),7.0)*pow(sin(theta/2.0),3.0);
	break;     
      default:
	printf("Error - %s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", func, s, l, m );
	exit(1);
	break;
      }
    }  /* l==5 */
    else {
      printf("Error - %s: Unsupported mode l=%d (only l in [2,5] implemented)\n", func, s);
      exit(1);
    }
  }
  else {
    printf("Error - %s: Unsupported mode s=%d (only s=-2 implemented)\n", func, s);
    exit(1);
  }
  if (m)
    ans = fac*cexp(I*m*phi);
  else
    ans = fac;
  return ans;
}
