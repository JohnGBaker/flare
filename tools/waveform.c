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

/***************** Function estimating frequency corresponding to a given time to coalescence ****************/

/* Functions computing relations between chirp mass, eta, m1, m2*/
double Mchirpofm1m2(
  const double m1,   /* Mass 1 */
  const double m2)   /* Mass 2 */
{
  return pow(m1*m2, 3./5) / pow(m1+m2, 1./5);
}
double etaofm1m2(
  const double m1,   /* Mass 1 */
  const double m2)   /* Mass 2 */
{
  return m1*m2 / pow(m1+m2, 2.);
}
double m1ofMchirpeta(
  const double Mchirp,   /* Chirp mass */
  const double eta)      /* Symmetric mass ratio */
{
  double delta = sqrt(1. - 4.*eta);
  return Mchirp * pow(eta, -3./5) * (1.+delta)/2.;
}
double m2ofMchirpeta(
  const double Mchirp,   /* Chirp mass */
  const double eta)      /* Symmetric mass ratio */
{
  double delta = sqrt(1. - 4.*eta);
  return Mchirp * pow(eta, -3./5) * (1.-delta)/2.;
}

/* NOTE: these are base on Newtonian estimates, and are not accurate */
/* Newtonian estimate of the relation Mf(deltat/M) (for the 22 mode) - gives the starting geometric frequency for a given mass ratio and a given geometric duration of the observations */
double NewtonianfoftGeom(
  const double q,                      /* Mass ratio m1/m2 */
  const double t)                      /* Duration of the observations in geometric units (t/M) */
{
  if(t<=0.) return 0.;
  double nu = q/(1.+q)/(1.+q);
  return 1./PI * pow(256*nu/5. * t, -3./8);
}
/* Newtonian estimate of the relation f(deltat) (for the 22 mode) - gives the starting geometric frequency for a given mass ratio and a given geometric duration of the observations - output in Hz */
double Newtonianfoft(
  const double m1,                     /* Mass 1 (solar masses) */
  const double m2,                     /* Mass 2 (solar masses) */
  const double t)                      /* Duration of the observations in years */
{
  if(t<=0.) return 0.;
  double mtot = m1 + m2;
  double q = m1/m2;
  return NewtonianfoftGeom(q, t*YRSID_SI / (mtot*MTSUN_SI)) / (mtot*MTSUN_SI);
}
/* Newtonian estimate of the relation f(deltat) (for the 22 mode freq) - gives the starting geometric frequency for a given time to merger and chirp mass - output in Hz */
double Newtonianfoftchirp(
  const double mchirp,                 /* Chirp mass (solar masses) */
  const double t)                      /* Time in years */
{
  if(t<=0.) return 0.;
  return 1./PI * pow(mchirp*MTSUN_SI, -5./8) * pow(256.*t*YRSID_SI/5, -3./8);
}
/* Newtonian estimate of the relation deltat(f) (for the 22 mode freq) - gives the time to merger from a starting frequency for a given chirp mass - output in years */
double Newtoniantoffchirp(
  const double mchirp,                 /* Chirp mass (solar masses) */
  const double f)                      /* Freq in Hz */
{
  return 5./256 * pow(mchirp*MTSUN_SI, -5./3) * pow(PI*f, -8./3) / YRSID_SI;
}

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
    if(minf==0) minf = gsl_vector_get(hlm->freq, 0);
    else minf = fmin(minf, gsl_vector_get(hlm->freq, 0));
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
  double fLow,                                                /* Minimal frequency - set to 0 to ignore */
  double fHigh,                                               /* Maximal frequency - set to 0 to ignore */
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
  /* Takes into account fLow, fHigh and fstartobs */
  int jStart = 0;
  int jStop = sizeout - 1;
  double minfmode = fmax(fLow, fmax(gsl_vector_get(freqin, 0), fstartobsmode));
  double maxfmode =  0.;
  if(fHigh==0.) maxfmode = gsl_vector_get(freqin, sizein - 1);
  else maxfmode = fmin(fHigh, gsl_vector_get(freqin, sizein - 1));
  while(jStart<sizeout && gsl_vector_get(freqout, jStart) < minfmode) jStart++;
  while( jStop > -1 && gsl_vector_get(freqout, jStop) > maxfmode ) jStop--; /* allow continuing all the way to an impossible value to make the later loop empty if need be */
  if(jStop <= jStart) {
    printf("Warning: empty range of frequencies in ReImFrequencySeries_AddCAmpPhaseFrequencySeries.\n");
    printf(" ...maybe this is OK, continuing with jStart=%i, jStop=%i\n", jStart, jStop);
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
  double fLow,                                                      /* Minimal frequency - set to 0 to ignore */
  double fHigh,                                                     /* Maximal frequency - set to 0 to ignore */
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
    ReImFrequencySeries_AddCAmpPhaseFrequencySeries(freqseriesReIm, listelement->freqseries, fLow, fHigh, fstartobsmode);
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
/* If nbpt is specified, use (nbpt-1)*deltaf as the maximal frequency of output */
/* If not, use the last n*deltaf below min(fHigh, fHighROM) as the maximal frequency to be generated */
/* Values will be set to 0 outside [fLow, fHigh] */
int GeneratehphcFDReImFrequencySeries(
  ReImFrequencySeries** hptilde,                 /* Output: frequency series for hplus */
  ReImFrequencySeries** hctilde,                 /* Output: frequency series for hcross */
  ListmodesCAmpPhaseFrequencySeries* listhlm,    /* Input: frequency series for the mode hlm  */
  double fLow,                                   /* Minimal frequency - set to 0 to ignore */
  double fHigh,                                  /* Maximal frequency - set to 0 to ignore */
  double fstartobs,                              /* For limited duration of observation, starting frequency for the 22 mode - set to 0 to ignore */
  double deltaf,                                 /* Frequency step */
  int nbpt,                                      /* Number of points of output - if 0, determined from deltaf and maximal frequency in input */
  int nbmode,                                    /* Number of modes to add */
  double Theta,                                  /* First angle for position in the sky of observer */
  double Phi,                                    /* Second angle for position in the sky of observer */
  int sym)                                       /* If 1, assume planar symmetry and add also mode l,-m. Do not if set to 0. */
{
  ListmodesCAmpPhaseFrequencySeries* hlm;

  /* Determine total size of frequency series if nbpt is not set */
  int nbpts = 0;
  double maxf = 0.;
  if(nbpt==0) {
    double fHighROM = ListmodesCAmpPhaseFrequencySeries_maxf(listhlm);
    // NOTE : there were problems with using int for the size of very long frequency series for very low mass signals - switch to size_t ? how to get reasonable size ?
    if(fHigh==0.) maxf = fHighROM;
    else maxf = fmin(fHigh, fHighROM);
    nbpts = ceil(maxf/deltaf);
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
    double fLow_mode = 0.;
    double fHigh_mode = 0.;
    double fstartobsmode = fmax(fstartobs, m/2. * fstartobs); /* note: does not fetch frequencies < fstartobs22 for m=1*/
    if(fLow==0.) fLow_mode = fmax(fstartobsmode, gsl_vector_get(freq_ds, 0));
    else fLow_mode = fmax(fstartobsmode, fmax(fLow, gsl_vector_get(freq_ds, 0)));
    if(fHigh==0.) fHigh_mode = gsl_vector_get(freq_ds, nbfreq-1);
    else fHigh_mode = fmin(fHigh, gsl_vector_get(freq_ds, nbfreq-1));
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
/* If nbpt is specified, use (nbpt-1)*deltaf as the maximal frequency of output */
/* If not, use the last n*deltaf below min(fHigh, fHighROM) as the maximal frequency to be generated */
/* Values will be set to 0 outside [fLow, fHigh] */
int GenerateFDReImFrequencySeries(
  ReImFrequencySeries** freqseries,              /* Output: frequency series */
  ListmodesCAmpPhaseFrequencySeries* listhlm,    /* Input: FD modes hlm in the form AmpReal/AmpIm/Phase  */
  double fLow,                                   /* Minimal frequency */
  double fHigh,                                  /* Maximal frequency */
  double fstartobs,                              /* For limited duration of observation, starting frequency for the 22 mode - set to 0 to ignore */
  double deltaf,                                 /* Frequency step */
  int nbpt)                                      /* Number of points of output - if 0, determined from deltaf and maximal frequency in input */
{
  /* Determine total size of frequency series if nbpt is not set */
  int nbpts = 0;
  double maxf = 0.;
  if(nbpt==0) {
    double fHighROM = ListmodesCAmpPhaseFrequencySeries_maxf(listhlm);
    if(fHigh==0.) maxf = fHighROM;
    else maxf = fmin(fHigh, fHighROM);
    nbpts = (int) ceil(maxf/deltaf);
  }
  else nbpts = nbpt;

  /* Initialize frequency series */
  ReImFrequencySeries_Init(freqseries, nbpts);
  /* Note: initialization of freqseries hreal and himag to 0 is done in ReImFrequencySeries_SumListmodesCAmpPhaseFrequencySeries  */

  /* Set values for frequency vector */
  gsl_vector* freq = gsl_vector_alloc(nbpts);
  double f;
  for(int i=0; i<nbpts; i++) {
    f = i*deltaf;
    gsl_vector_set(freq, i, f);
  }

  /* Sum mode contributions */
  /* NOTE: fLow is used as fstartobs here, different meaning in spirit in the rest of the code */
  ReImFrequencySeries_SumListmodesCAmpPhaseFrequencySeries((*freqseries), listhlm, freq, fLow, fHigh, fstartobs);

  /* Cleanup */
  /* freq has been memcopied inside ReImFrequencySeries_SumListmodesCAmpPhaseFrequencySeries */
  gsl_vector_free(freq);
}

/* Function evaluating the FD frequency series for a single mode contribution (l,m) */
int GenerateFDReImFrequencySeriesSingleMode(
  ReImFrequencySeries** freqseries,              /* Output: frequency series */
  ListmodesCAmpPhaseFrequencySeries* listhlm,    /* Input: FD modes hlm in the form AmpReal/AmpIm/Phase  */
  double fLow,                                   /* Minimal frequency - set to 0 to ignore */
  double fHigh,                                  /* Maximal frequency - set to 0 to ignore */
  double fstartobs,                              /* For limited duration of observation, starting frequency for the 22 mode - set to 0 to ignore */
  double deltaf,                                 /* Frequency step */
  int nbpt,                                      /* Number of points of output - if 0, determined from deltaf and maximal frequency in input */
  int l,                                         /* Mode index l */
  int m)                                         /* Mode index m */
{
  /* Get relevant mode */
  CAmpPhaseFrequencySeries* hlm = ListmodesCAmpPhaseFrequencySeries_GetMode(listhlm, l, m)->freqseries;
  /* Determine total size of frequency series if nbpt is not set */
  int nbpts = 0;
  double maxf = 0.;
  if(nbpt==0) {
    double fHighROM = gsl_vector_get(hlm->freq, hlm->freq->size - 1);
    if(fHigh==0.) maxf = fHighROM;
    else maxf = fmin(fHigh, fHighROM);
    nbpts = (int) ceil(maxf/deltaf);
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

  /* Copy frequencies */
  gsl_vector_memcpy((*freqseries)->freq, freq);

  /* Take the chosen mode, interpolate it and add this single mode to the output */
  ReImFrequencySeries_AddCAmpPhaseFrequencySeries((*freqseries), hlm, fLow, fHigh, fstartobs);

  /* Clean up */
  /* freq has been memcopied */
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

/* Function for getting a phase mod 2pi (between -pi and pi) */
double mod2pi(double phase) {
  return phase - floor((phase + PI) / (2*PI)) * (2*PI);
}
/* Function for getting a phase mod pi (between 0 and pi, e.g. polarization) */
double modpi(double phase) {
  return phase - floor(phase / PI) * PI;
}

/* Function to unwrap the phase mod 2pi  - acts directly on a gsl_vector representing the phase */
/* WARNING : found significant accumulating difference with the numpy unwrap function - to be understood */
int UnwrapPhase(
  gsl_vector* phaseout,    /* Output: unwrapped phase vector, already allocated */
  gsl_vector* phasein)     /* Input: phase vector */
{
  size_t N = phasein->size;
  double* p = phasein->data;
  double* pmod = (double*) malloc(sizeof(double) * N);
  int* jumps = (int*) malloc(sizeof(int) * N);
  int* cumul = (int*) malloc(sizeof(int) * N);

  /* Compute phase mod 2pi (shifted to be between -pi and pi) */
  for(int i=0; i<N; i++) {
    pmod[i] = p[i] - floor((p[i] + PI) / (2*PI))*(2*PI);
  }

  /* Identify jumps */
  jumps[0] = 0;
  double d = 0.;
  for(int i=1; i<N; i++) {
    d = pmod[i] - pmod[i-1];
    if(d<-PI) jumps[i] = 1;
    else if(d>PI) jumps[i] = -1;
    else jumps[i] = 0;
  }

  /* Cumulative of the jump sequence */
  int c = 0;
  cumul[0] = 0;
  for(int i=1; i<N; i++) {
    c += jumps[i];
    cumul[i] = c;
  }

  /* Correct original phase series by the number of 2pi factor given by the cumulative of the jumps */
  double* pout = phaseout->data;
  for(int i=0; i<N; i++) {
    pout[i] = 2*PI*cumul[i] + p[i];
  }

  /* Cleanup */
  free(pmod);
  free(jumps);
  free(cumul);
}

/* Function to convert a time series from Re/Im form to Amp/Phase form - unwrapping the phase */
int ReImTimeSeries_ToAmpPhase(
  AmpPhaseTimeSeries** timeseriesout,             /* Output: Amp/Phase time series */
  ReImTimeSeries* timeseriesin)                   /* Input: Re/Im time series */
{
  /* Initialize output */
  int npt = (int) timeseriesin->times->size;
  AmpPhaseTimeSeries_Init(timeseriesout, npt);
  gsl_vector_memcpy((*timeseriesout)->times, timeseriesin->times);

  /* Compute amplitude and phase */
  gsl_vector* amp = gsl_vector_alloc(npt);
  gsl_vector* phase = gsl_vector_alloc(npt);
  double* ampval = amp->data;
  double* phaseval = phase->data;
  double* hreal = timeseriesin->h_real->data;
  double* himag = timeseriesin->h_imag->data;
  for(int i=0; i<npt; i++) {
    ampval[i] = cabs(hreal[i] + I*himag[i]);
    phaseval[i] = carg(hreal[i] + I*himag[i]);
  }

  /* Unwrap phase */
  UnwrapPhase((*timeseriesout)->h_phase, phase);

  /* Copy amplitude */
  gsl_vector_memcpy((*timeseriesout)->h_amp, amp);

  /* Cleanup */
  gsl_vector_free(amp);
  gsl_vector_free(phase);
}

/* Function to convert a time series from Amp/Phase form to Re/Im form */
int AmpPhaseTimeSeries_ToReIm(
  ReImTimeSeries** timeseriesout,                 /* Output: Re/Im time series */
  AmpPhaseTimeSeries* timeseriesin)               /* Input: Amp/Phase time series */
{
  /* Initialize output */
  int npt = (int) timeseriesin->times->size;
  ReImTimeSeries_Init(timeseriesout, npt);
  gsl_vector_memcpy((*timeseriesout)->times, timeseriesin->times);
  gsl_vector_set_zero((*timeseriesout)->h_real);
  gsl_vector_set_zero((*timeseriesout)->h_imag);

  /* Compute amplitude and phase */
  gsl_vector* hreal = gsl_vector_alloc(npt);
  gsl_vector* himag = gsl_vector_alloc(npt);
  double* hrealval = (*timeseriesout)->h_real->data;
  double* himagval = (*timeseriesout)->h_imag->data;
  double* hamp = timeseriesin->h_amp->data;
  double* hphase = timeseriesin->h_phase->data;
  for(int i=0; i<npt; i++) {
    hrealval[i] = hamp[i]*cos(hphase[i]);
    himagval[i] = hamp[i]*sin(hphase[i]);
  }

  return SUCCESS;
}

/* Function to compute a linear resampling at high frequencies to enforce a maximal deltaf */
/* Useful to resolve the high-f structure in the LISA response */
/* NOTE: Assumes input frequencies are logarithmic (except maybe first interval) to evaluate when to resample */
int SetMaxdeltafResampledFrequencies(
  gsl_vector** freqr,              /* Output: resampled frequencies */
  gsl_vector* freq,                /* Input: original frequencies */
  const double maxf,               /* Input: maximal frequency - set to 0. to ignore */
  const double deltaf)             /* Input: maximal deltaf aimed for - 0.002Hz appropriate for LISA */
{
  int npt = (int) freq->size;
  double* f =  freq->data;
  double ratio = f[2]/f[1]; /* The first interval might be shorter than the other ones - ignore it */
  double minf1 = f[1];
  double fHigh = 0;
  if(maxf<=0.) fHigh = f[npt-1];
  else if(maxf<f[0]) {
    printf("Error in MaxdeltafResampledFrequencies: maxf lower than first frequency present.\n");
    exit(1);
  }
  else fHigh = fmin(maxf, f[npt-1]);
  int indexrs = floor(log(deltaf/(ratio-1.)/minf1)/log(ratio) + 1);
  if((indexrs>npt-1) || (f[indexrs]>=fHigh)) { /* here no need to resample - include all samples up to and including fHigh */
    int index = 0;
    while((index+1<npt) && (f[index+1]<=fHigh)) index++;
    if(f[index]<fHigh) { /* append fHigh - note that we have fHigh<=f[npt-1] */
      *freqr = gsl_vector_alloc(index+2);
      double* fr = (*freqr)->data;
      for(int j=0; j<=index; j++) fr[j] = f[j];
      fr[index+1] = fHigh;
    }
    else { /* f[index] is already fHigh */
      *freqr = gsl_vector_alloc(index+1);
      double* fr = (*freqr)->data;
      for(int j=0; j<=index; j++) fr[j] = f[j];
    }
  }
  else { /* here f[indexrs]<fHigh, we keep original frequencies up to indexrs-1, and extend with a linspace(f[indexrs], fHigh, nrs) with nrs chosen so that the linear deltaf meets the requirement */
    int nrs = ceil((fHigh - f[indexrs])/deltaf) + 1;
    int ntot = indexrs + nrs;
    *freqr = gsl_vector_alloc(ntot);
    double* fr = (*freqr)->data;
    for(int j=0; j<indexrs; j++) fr[j] = f[j];
    double r = 1./(ntot-1-indexrs) * (fHigh - f[indexrs]);
    for(int j=indexrs; j<ntot; j++) fr[j] = fmin(f[indexrs] + (j-indexrs) * r, fHigh); /* ensure that we do not go beyond fHigh due to numerical errors */
  }

  return SUCCESS;
}

// Mathematica:
// FrequenciesLinearResampling[freqs_, deltaf_] :=
//   Module[{ratio, fmin2, fmax, indexrs, nrs},
//    ratio = freqs[[3]]/freqs[[2]];
//    fmin2 = freqs[[2]];
//    fmax = freqs[[-1]];
//    indexrs = Floor[Log[deltaf/(ratio - 1)/fmin2]/Log[ratio] + 2];
//    If[indexrs >= Length@freqs, Return@freqs];
//    nrs = Ceiling[(freqs[[-1]] - freqs[[indexrs]])/deltaf];
//    Return[
//     Join[freqs[[;; indexrs - 1]],
//      Linspace[freqs[[indexrs]], fmax, nrs]]];
//    ];

/* Function to compute a linear-in-time resampling at low frequencies to enforce a maximal deltat */
/* Useful to resolve the low-f structure in the LISA response when accumulating signal over several years */
/* Assumes input frequencies correspond to decreasing deltat */
/* Uses Newtonian estimates for t<->f, only approximate  */
/* Requires to be passed the chirpmass, as well as mode number m */
int SetMaxdeltatResampledFrequencies(
  gsl_vector** freqr,              /* Output: resampled frequencies */
  gsl_vector* freq,                /* Input: original frequencies */
  const double deltat,             /* Input: maximal deltat aimed for - fraction of a year, 1/24 (half month) appropriate for 1e-4 interpolation errors */
  const double mchirp,             /* Input: chirp mass, used for approximate t-f correspondence */
  const int m)                     /* Input: chirp mass, used for approximate t-f correspondence */
{
  /* Computed the t(f) according to Newtonian estimate */
  /* NOTE : computing at highf is wasted computation, but we try to be generic */
  int npt = (int) freq->size;
  double* f = freq->data;
  gsl_vector* times = gsl_vector_alloc(npt);
  double* t = times->data;
  for(int j=0; j<npt; j++) t[j] = Newtoniantoffchirp(mchirp, f[j]*2./m); /* include rescaling for modes other than 22 */

  /* Determine where to resample */
  int indexrs = 0;
  while((indexrs+1<npt) && ((t[indexrs] - t[indexrs+1])>deltat)) indexrs++;

  if(indexrs==0) { /* No resampling */
    *freqr = gsl_vector_alloc(npt);
    gsl_vector_memcpy(*freqr, freq);
  }
  else {
    int nrs = floor((t[0] - t[indexrs]) / deltat);
    int ntot = nrs + npt - indexrs;
    *freqr = gsl_vector_alloc(ntot);
    double* fr = (*freqr)->data;
    double r = (t[indexrs] - t[0])/nrs;
    fr[0] = f[0]; /* avoid changing the boundary by numerical errors */
    for(int j=1; j<nrs; j++) fr[j] = m/2. * Newtonianfoftchirp(mchirp, t[0] + j*r); /* include rescaling for modes other than 22 */
    for(int j=nrs; j<ntot; j++) fr[j] = f[j - nrs + indexrs];
  }
  /* Cleanup */
  gsl_vector_free(times);

  return SUCCESS;
}

//Mathematica
// funcTimeResampling[freqs_, deltat_, funct_, funcf_] :=
//   Module[{times, indexrs, nrs, timesprepend, freqsprepend},
//    times = funct /@ freqs;
//    indexrs = 1;
//    While[indexrs <
//       Length@freqs && (times[[indexrs]] - times[[indexrs + 1]] >
//        deltat), indexrs += 1];
//    If[indexrs == 1, Return@freqs];(* No need to resample *)
//
//    nrs = Floor[(times[[1]] - times[[indexrs]])/deltat];
//    timesprepend = Linspace[times[[1]], times[[indexrs]], nrs + 1];
//    freqsprepend = Drop[funcf /@ timesprepend, -1];
//    freqsprepend[[1]] = freqs[[1]];(*
//    avoid problems with numerical error *)
//
//    Return[Join[freqsprepend, freqs[[indexrs ;;]]]];
//    ];

/* Function to resample a CAmp/Phase frequency series on the specified frequencies */
int CAmpPhaseFrequencySeries_Resample(
  CAmpPhaseFrequencySeries** freqseriesout,         /* Output: CAmp/Phase freq series */
  CAmpPhaseFrequencySeries* freqseriesin,           /* Input: CAmp/Phase freq series */
  gsl_vector* freqr)                                  /* Input: freq vector to resample on */
{
  /* Definitions */
  gsl_vector* freq = freqseriesin->freq;
  gsl_vector* amp_real = freqseriesin->amp_real;
  gsl_vector* amp_imag = freqseriesin->amp_imag;
  gsl_vector* phase = freqseriesin->phase;
  int npt = (int) freq->size;
  int nptr = (int) freqr->size;

  /* Check boundaries */
  if((gsl_vector_get(freqr, 0)<gsl_vector_get(freq, 0)) || (gsl_vector_get(freqr, nptr-1)>gsl_vector_get(freq, npt-1))) {
    printf("Error in CAmpPhaseFrequencySeries_Resample: resampling exceeding bounds.\n");
    exit(1);
  }

  /* Initialize gsl splines */
  gsl_spline* spline_amp_real = gsl_spline_alloc(gsl_interp_cspline, npt);
  gsl_spline* spline_amp_imag = gsl_spline_alloc(gsl_interp_cspline, npt);
  gsl_spline* spline_phase = gsl_spline_alloc(gsl_interp_cspline, npt);
  gsl_interp_accel* accel_amp_real = gsl_interp_accel_alloc();
  gsl_interp_accel* accel_amp_imag = gsl_interp_accel_alloc();
  gsl_interp_accel* accel_phase = gsl_interp_accel_alloc();
  gsl_spline_init(spline_amp_real, gsl_vector_const_ptr(freq, 0), gsl_vector_const_ptr(amp_real, 0), npt);
  gsl_spline_init(spline_amp_imag, gsl_vector_const_ptr(freq, 0), gsl_vector_const_ptr(amp_imag, 0), npt);
  gsl_spline_init(spline_phase, gsl_vector_const_ptr(freq, 0), gsl_vector_const_ptr(phase, 0), npt);

  /* Initialize output frequency series */
  CAmpPhaseFrequencySeries_Init(freqseriesout, nptr);
  gsl_vector_memcpy((*freqseriesout)->freq, freqr);
  double* areal = (*freqseriesout)->amp_real->data;
  double* aimag = (*freqseriesout)->amp_imag->data;
  double* phi = (*freqseriesout)->phase->data;

  /* Evaluate on the new frequencies */
  /* NOTE: with resampling, some of the frequency samples are repeats of the original frequency series - so we could save some spline evaluations here */
  double f;
  for(int j=0; j<nptr; j++) {
    f = gsl_vector_get(freqr, j);
    areal[j] = gsl_spline_eval(spline_amp_real, f, accel_amp_real);
    aimag[j] = gsl_spline_eval(spline_amp_imag, f, accel_amp_imag);
    phi[j] = gsl_spline_eval(spline_phase, f, accel_phase);
  }

  /* Cleanup */
  gsl_spline_free(spline_amp_real);
  gsl_spline_free(spline_amp_imag);
  gsl_spline_free(spline_phase);
  gsl_interp_accel_free(accel_amp_real);
  gsl_interp_accel_free(accel_amp_imag);
  gsl_interp_accel_free(accel_phase);

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
