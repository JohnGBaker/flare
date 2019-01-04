from cython_gsl cimport *
import numpy as np
cimport numpy as np
import ctypes

cdef extern from "struct.h":

    ctypedef struct CAmpPhaseSpline:
      gsl_matrix* spline_amp_real;
      gsl_matrix* spline_amp_imag;
      gsl_matrix* quadspline_phase;

    ctypedef struct CAmpPhaseFrequencySeries:
      gsl_vector* freq;
      gsl_vector* amp_real;
      gsl_vector* amp_imag;
      gsl_vector* phase;


    ctypedef struct CAmpPhaseGSLSpline:
      gsl_vector* freq;
      gsl_spline* spline_amp_real;
      gsl_spline* spline_amp_imag;
      gsl_spline* spline_phase;
      gsl_interp_accel* accel_amp_real;
      gsl_interp_accel* accel_amp_imag;
      gsl_interp_accel* accel_phase;


    ctypedef struct ReImFrequencySeries:
      gsl_vector* freq;
      gsl_vector* h_real;
      gsl_vector* h_imag;

    ctypedef struct ReImTimeSeries:
      gsl_vector* times;
      gsl_vector* h_real;
      gsl_vector* h_imag;

    ctypedef struct AmpPhaseTimeSeries:
      gsl_vector* times;
      gsl_vector* h_amp;
      gsl_vector* h_phase;

    ctypedef double (*RealFunctionPtr)(double);
    ctypedef double (*RealObjectFunctionPtr)(const void *, double);

    ctypedef struct ObjectFunction:
      const void * object;
      RealObjectFunctionPtr function;

    double ObjectFunctionCall(const ObjectFunction*, double);

    void CAmpPhaseSpline_Cleanup(CAmpPhaseSpline *splines);


cdef extern from "splinecoeffs.h":
    void BuildNotAKnotSpline(
      gsl_matrix* splinecoeffs,
      gsl_vector* vectx,
      gsl_vector* vecty,
      int n);

    void BuildQuadSpline(
      gsl_matrix* splinecoeffs,
      gsl_vector* vectx,
      gsl_vector* vecty,
      int n);

    void BuildSplineCoeffs(
      CAmpPhaseSpline** splines,
      CAmpPhaseFrequencySeries* freqseries);

    # void BuildListmodesCAmpPhaseSpline(
    #   ListmodesCAmpPhaseSpline** listspline,
    #   ListmodesCAmpPhaseFrequencySeries* listh);

    double EvalCubic(
      gsl_vector* coeffs,
      double eps,
      double eps2,
      double eps3);

    double EvalQuad(
      gsl_vector* coeffs,
      double eps,
      double eps2);

    void EvalCAmpPhaseSpline(
      CAmpPhaseSpline* splines,
      CAmpPhaseFrequencySeries* freqseries)


cdef extern from "likelihoodKurz.h":

    double FDSinglemodeLogLinearOverlap(
      CAmpPhaseFrequencySeries *freqseries1,
      CAmpPhaseFrequencySeries *freqseries2,
      ObjectFunction * Snoise,
      double fLow,
      double fHigh)

    double FDOverlapReImvsReIm(
      ReImFrequencySeries *h1,
      ReImFrequencySeries *h2,
      gsl_vector* noisevalues)


    double FDSinglemodeFresnelOverlap(
      CAmpPhaseFrequencySeries *freqseries1,
      CAmpPhaseSpline *splines2,
      ObjectFunction * Snoise,
      double fLow,
      double fHigh)

    double FDSinglemodeFresnelOverlap3Chan(
      CAmpPhaseFrequencySeries *freqseries1chan1,
      CAmpPhaseFrequencySeries *freqseries1chan2,
      CAmpPhaseFrequencySeries *freqseries1chan3,
      CAmpPhaseSpline *splines2chan1,
      CAmpPhaseSpline *splines2chan2,
      CAmpPhaseSpline *splines2chan3,
      ObjectFunction * Snoisechan1,
      ObjectFunction * Snoisechan2,
      ObjectFunction * Snoisechan3,
      double fLow,
      double fHigh)


    void ComputeIntegrandValues(
      CAmpPhaseFrequencySeries** integrand,
      CAmpPhaseFrequencySeries* freqseries1,
      CAmpPhaseSpline* splines2,
      ObjectFunction * Snoise,
      double fLow,
      double fHigh)


    double FDSinglemodeWIPOverlap(
      CAmpPhaseFrequencySeries *freqseries1,
      CAmpPhaseFrequencySeries *freqseries2,
      ObjectFunction * Snoise,
      double fLow,
      double fHigh)

    double SnAProposal(double f)

    double SnEProposal(double f)

    double SnTProposal(double f)

    double SnXProposal(double f)











##############################################################################################


def BuildNotAKspline(np.ndarray[np.float_t, ndim=1] v_x, np.ndarray[np.float_t, ndim=1] v_y):

    cdef gsl_vector_view vec_x
    cdef gsl_vector_view vec_y
    cdef int N
    N = len(v_x)
    if (len(v_y) != N):
        print ("size of arrays must be the same", N, len(v_y))
        raise ValueError

    vec_x = gsl_vector_view_array(<double*> v_x.data,  N)
    vec_y = gsl_vector_view_array(<double*> v_y.data,  N)


    cdef gsl_matrix* spline_coeffs
    spline_coeffs =  gsl_matrix_alloc(N, 5)

    BuildNotAKnotSpline(spline_coeffs, &(vec_x.vector), &(vec_y.vector), N)

    ### now we need to cast it to  numoy array
    cdef double[:] spl_view = <(double)[:N*5]>spline_coeffs.data

    out = np.asarray(spl_view)
    out = out.reshape((N,5))

    return (out)

def BuildQuadspline(np.ndarray[np.float_t, ndim=1] v_x, np.ndarray[np.float_t, ndim=1] v_y):
    ### Building vectors from doubles

    # print (v_x)
    # print (v_y)
    cdef gsl_vector_view vec_x
    cdef gsl_vector_view vec_y
    cdef int N
    N = len(v_x)
    if (len(v_y) != N):
        print ("size of arrays must be the same", N, len(v_y))
        raise ValueError

    vec_x = gsl_vector_view_array(<double*> v_x.data,  N)
    vec_y = gsl_vector_view_array(<double*> v_y.data,  N)


    cdef gsl_matrix* spline_coeffs
    spline_coeffs =  gsl_matrix_alloc(N, 4)

    BuildQuadSpline(spline_coeffs, &(vec_x.vector), &(vec_y.vector), N)

    ### now we need to cast it to  numoy array

    # cdef double [:,:] v_view = <(double)[:N]>(v.vector).data
    #cdef double [:,:] spl_view = <double**>spline_coeffs.data


    cdef double[:] spl_view = <(double)[:N*4]>spline_coeffs.data

    out = np.asarray(spl_view)
    out = out.reshape((N,4))

    return (out)


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


def Arr_Matr(np.ndarray[np.float_t, ndim=2] a):

    cdef int N, M
    cdef int i
    cdef gsl_matrix_view mv
    N,M = np.shape(a)
    dat = <double *> a.data

    # print ("here", N, M, a)

    mv = gsl_matrix_view_array(dat,  N, M)

    ### Let's do some operation with vector
    for i in range(N):
        for j in range(M):
            gsl_matrix_set(&(mv.matrix), i, j, gsl_matrix_get (&(mv.matrix), i, j)*3.0)


    # b = np.zeros((N,M))
    # cdef double [:, :] b_view = b
    cdef double [:] mv_view = <(double)[:N*M]>(mv.matrix).data
    # cdef double [::1] v_view = <(double)[:N]>(v.vector).data

    b = np.asarray(mv_view)
    b = b.reshape((N,M))

    # b.data = (v.vector).data

    # for i in range(N):
    #     print( gsl_vector_get (&(v.vector), i) )


    return(b)

def GetNoiseXAET(double f):
    cdef double SA = SnAProposal(f)
    cdef double SE = SnEProposal(f)
    cdef double ST = SnTProposal(f)
    cdef double SX = SnXProposal(f)

    return (SX, SA, SE, ST)

def GetNoiseX(double f):
    cdef double SX = SnXProposal(f)
    return (SX)

def GetNoiseAE(double f):
    cdef double SA = SnAProposal(f)
    return (SA)


### try to use ObjectFunction
def GetNoise(double f, noise="X"):
    cdef ObjectFunction PSD


    PSD.object = NULL
    if (noise == "X"):
        PSD.function = <RealObjectFunctionPtr> SnXProposal
    elif (noise == "AE"):
        PSD.function = <RealObjectFunctionPtr> SnAProposal

    cdef double Sn = ObjectFunctionCall(&PSD, f)

    return (Sn)



# double FDSinglemodeWIPOverlap(
#   CAmpPhaseFrequencySeries *freqseries1,
#   CAmpPhaseFrequencySeries *freqseries2,
#   ObjectFunction * Snoise,
#   double fLow,
#   double fHigh)

def OverlapFrensel1(FD_wvf1, FD_wvf2, fmin=1.e-5, fmax=0.1, noise="X"):

    cdef int N1 = len(FD_wvf1["freq"])
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] fr1 = FD_wvf1["freq"]
    amp1 = FD_wvf1["ampl"]
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] ph1 = np.copy(FD_wvf1["phase"], order = 'C')
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] Ramp1 =  np.copy(np.real(amp1), order='C')
    # print (Ramp1.flags)
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] Iamp1 = np.copy(np.imag(amp1), order='C')

    cdef int N2 = len(FD_wvf2["freq"])
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] fr2 = FD_wvf2["freq"]
    amp2 = FD_wvf2["ampl"]
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] ph2 = np.copy(FD_wvf2["phase"], order = 'C')
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] Ramp2 =  np.copy(np.real(amp2), order='C')
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] Iamp2 = np.copy(np.imag(amp2), order='C')

    cdef double fLow = fmin
    cdef double fHigh = fmax
    cdef ObjectFunction PSD


    PSD.object = NULL
    if (noise == "X"):
        PSD.function = <RealObjectFunctionPtr> SnXProposal
    elif (noise == "AE"):
        # PSD.function = <RealObjectFunctionPtr> SnTProposal
        PSD.function = <RealObjectFunctionPtr> SnAProposal


    cdef gsl_vector_view freq1
    cdef gsl_vector_view ampl1Re
    cdef gsl_vector_view ampl1Im
    cdef gsl_vector_view phase1

    cdef double* fr1_dat = <double *>fr1.data#.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    # cdef double* am1r_dat = <double *>Ramp1.data
    cdef double* am1r_dat = <double *>Ramp1.data#.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    cdef double* am1i_dat = <double *> Iamp1.data#.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    # am1r_dat = Ramp1
    # am1i_dat = <double *>Iamp1.data
    freq1 = gsl_vector_view_array(fr1_dat,  N1)
    ampl1Re = gsl_vector_view_array(am1r_dat, N1)
    ampl1Im = gsl_vector_view_array(am1i_dat, N1)
    phase1 = gsl_vector_view_array(<double *>ph1.data, N1)
    # ampl1Re = gsl_vector_view_array(<double *>Ramp1.data, N1)
    # ampl1Im = gsl_vector_view_array(<double *>Iamp1.data, N1)
    # phase1 = gsl_vector_view_array(<double *>ph1.data, N1)

    # fout = open("TestOlap2.dat", 'w')
    # for i in range(len(fr1)):
    #     # rec = str(fr1[i]) + "  " + str(Ramp1[i]**2 + Iamp1[i]**2) + " \n"
    #     # are = gsl_vector_get(apfs1.amp_real, i)
    #     # aim = gsl_vector_get(apfs1.amp_imag, i)
    #     are  = am1r_dat[i]
    #     # are  = Ramp1[i]
    #     aim = am1i_dat[i]
    #     rec = str(fr1[i]) + "  " + str(are**2 + aim**2) + "   " + str(are) + "    " +  str(aim)  + " \n"
    #     # print (rec)
    #     fout.write(rec)
    # fout.close()

    cdef gsl_vector_view freq2
    cdef gsl_vector_view ampl2Re
    cdef gsl_vector_view ampl2Im
    cdef gsl_vector_view phase2
    freq2 = gsl_vector_view_array(<double *>fr2.data,  N2)
    ampl2Re = gsl_vector_view_array(<double *>Ramp2.data, N2)
    ampl2Im = gsl_vector_view_array(<double *>Iamp2.data, N2)
    phase2 = gsl_vector_view_array(<double *>ph2.data, N2)

    cdef CAmpPhaseFrequencySeries apfs1
    cdef CAmpPhaseFrequencySeries apfs2

    apfs1.freq = &(freq1.vector)
    apfs1.amp_real = &(ampl1Re.vector)
    apfs1.amp_imag = &(ampl1Im.vector)
    apfs1.phase = &(phase1.vector)

    apfs2.freq = &(freq2.vector)
    apfs2.amp_real = &(ampl2Re.vector)
    apfs2.amp_imag = &(ampl2Im.vector)
    apfs2.phase = &(phase2.vector)

    ### Computing spline for wvf1:
    # void BuildSplineCoeffs(
    #   CAmpPhaseSpline** splines,                  /*  */
    #   CAmpPhaseFrequencySeries* freqseries)       /*  */
    # {


    cdef CAmpPhaseSpline* splines = NULL
    BuildSplineCoeffs(&splines, &apfs1)

    cdef double Olap = 0.0
    cdef double Olap2 = 0.0

    Olap = FDSinglemodeFresnelOverlap(&apfs2, splines, &PSD, fLow, fHigh)

    # Olap2 = FDSinglemodeLogLinearOverlap(&apfs1, &apfs2, &PSD, fLow, fHigh)
    #
    # print ("Overlaps:", Olap, Olap2, np.sqrt(Olap2))

    # Olap = FDSinglemodeWIPOverlap(&apfs1, &apfs2, &PSD, fLow, fHigh)
    # ouble FDSinglemodeFresnelOverlap(
    #   struct tagCAmpPhaseFrequencySeries *freqseries1, /* First mode h1, in amplitude/phase form */
    #   struct tagCAmpPhaseSpline *splines2,             /* Second mode h2, already interpolated in matrix form */
    #   ObjectFunction * Snoise,                  /* Noise function */
    #   double fLow,                                     /* Lower bound of the frequency window for the detector */
    #   double fHigh)

    CAmpPhaseSpline_Cleanup(splines)

    return (Olap)


#### This function computes inner product assuming three noise independent channels A, E, T
def OverlapFrensel3(FD_wvf1, FD_wvf2, fmin=1.e-5, fmax=0.1, noise="AE"):

    cdef int N1 = len(FD_wvf1["freq"])
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] fr1 = FD_wvf1["freq"]
    amp1 = FD_wvf1["ampl"]
    AA1 = amp1[0, :]
    EE1 = amp1[1, :]
    TT1 = amp1[2, :]
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] ph1 = np.copy(FD_wvf1["phase"], order = 'C')
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] RA1 =  np.copy(np.real(AA1), order='C')
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] IA1 = np.copy(np.imag(AA1), order='C')
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] RE1 =  np.copy(np.real(EE1), order='C')
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] IE1 = np.copy(np.imag(EE1), order='C')
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] RT1 =  np.copy(np.real(TT1), order='C')
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] IT1 = np.copy(np.imag(TT1), order='C')
    # print (Ramp1.flags)

    cdef int N2 = len(FD_wvf2["freq"])
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] fr2 = FD_wvf2["freq"]
    amp2 = FD_wvf2["ampl"]
    AA2 = amp2[0, :]
    EE2 = amp2[1, :]
    TT2 = amp2[2, :]
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] ph2 = np.copy(FD_wvf2["phase"], order = 'C')
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] RA2 =  np.copy(np.real(AA2), order='C')
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] IA2 = np.copy(np.imag(AA2), order='C')
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] RE2 =  np.copy(np.real(EE2), order='C')
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] IE2 = np.copy(np.imag(EE2), order='C')
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] RT2 =  np.copy(np.real(TT2), order='C')
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] IT2 = np.copy(np.imag(TT2), order='C')


    cdef double fLow = fmin
    cdef double fHigh = fmax
    cdef ObjectFunction PSDA
    cdef ObjectFunction PSDT


    PSDA.object = NULL
    PSDT.object = NULL
    if (noise == "X"):
        PSDA.function = <RealObjectFunctionPtr> SnXProposal
    elif (noise == "AE"):
        PSDA.function = <RealObjectFunctionPtr> SnAProposal
        PSDT.function = <RealObjectFunctionPtr> SnTProposal

    ### First data set

    cdef gsl_vector_view freq1
    cdef gsl_vector_view A1Re
    cdef gsl_vector_view A1Im
    cdef gsl_vector_view E1Re
    cdef gsl_vector_view E1Im
    cdef gsl_vector_view T1Re
    cdef gsl_vector_view T1Im
    cdef gsl_vector_view phase1

    cdef double* fr1_dat = <double *> fr1.data
    cdef double* A1r_dat = <double *> RA1.data
    cdef double* A1i_dat = <double *> IA1.data
    cdef double* E1r_dat = <double *> RE1.data
    cdef double* E1i_dat = <double *> IE1.data
    cdef double* T1r_dat = <double *> RT1.data
    cdef double* T1i_dat = <double *> IT1.data

    freq1 = gsl_vector_view_array(fr1_dat,  N1)
    phase1 = gsl_vector_view_array(<double *>ph1.data, N1)
    A1Re = gsl_vector_view_array(A1r_dat, N1)
    A1Im = gsl_vector_view_array(A1i_dat, N1)
    E1Re = gsl_vector_view_array(E1r_dat, N1)
    E1Im = gsl_vector_view_array(E1i_dat, N1)
    T1Re = gsl_vector_view_array(T1r_dat, N1)
    T1Im = gsl_vector_view_array(T1i_dat, N1)

    ### Second data set

    cdef gsl_vector_view freq2
    cdef gsl_vector_view A2Re
    cdef gsl_vector_view A2Im
    cdef gsl_vector_view E2Re
    cdef gsl_vector_view E2Im
    cdef gsl_vector_view T2Re
    cdef gsl_vector_view T2Im
    cdef gsl_vector_view phase2
    freq2 = gsl_vector_view_array(<double *>fr2.data,  N2)
    phase2 = gsl_vector_view_array(<double *>ph2.data, N2)
    A2Re = gsl_vector_view_array(<double *>RA2.data, N2)
    A2Im = gsl_vector_view_array(<double *>IA2.data, N2)
    E2Re = gsl_vector_view_array(<double *>RE2.data, N2)
    E2Im = gsl_vector_view_array(<double *>IE2.data, N2)
    T2Re = gsl_vector_view_array(<double *>RT2.data, N2)
    T2Im = gsl_vector_view_array(<double *>IT2.data, N2)

    cdef CAmpPhaseFrequencySeries Afs1
    cdef CAmpPhaseFrequencySeries Efs1
    cdef CAmpPhaseFrequencySeries Tfs1

    cdef CAmpPhaseFrequencySeries Afs2
    cdef CAmpPhaseFrequencySeries Efs2
    cdef CAmpPhaseFrequencySeries Tfs2

    Afs1.freq = &(freq1.vector)
    Afs1.amp_real = &(A1Re.vector)
    Afs1.amp_imag = &(A1Im.vector)
    Afs1.phase = &(phase1.vector)
    Efs1.freq = &(freq1.vector)
    Efs1.amp_real = &(E1Re.vector)
    Efs1.amp_imag = &(E1Im.vector)
    Efs1.phase = &(phase1.vector)
    Tfs1.freq = &(freq1.vector)
    Tfs1.amp_real = &(T1Re.vector)
    Tfs1.amp_imag = &(T1Im.vector)
    Tfs1.phase = &(phase1.vector)

    Afs2.freq = &(freq2.vector)
    Afs2.amp_real = &(A2Re.vector)
    Afs2.amp_imag = &(A2Im.vector)
    Afs2.phase = &(phase2.vector)
    Efs2.freq = &(freq2.vector)
    Efs2.amp_real = &(E2Re.vector)
    Efs2.amp_imag = &(E2Im.vector)
    Efs2.phase = &(phase2.vector)
    Tfs2.freq = &(freq2.vector)
    Tfs2.amp_real = &(T2Re.vector)
    Tfs2.amp_imag = &(T2Im.vector)
    Tfs2.phase = &(phase2.vector)


    cdef CAmpPhaseSpline* Asplines = NULL
    cdef CAmpPhaseSpline* Esplines = NULL
    cdef CAmpPhaseSpline* Tsplines = NULL
    BuildSplineCoeffs(&Asplines, &Afs1)
    BuildSplineCoeffs(&Esplines, &Efs1)
    BuildSplineCoeffs(&Tsplines, &Tfs1)

    cdef double Olap = 0.0
    # cdef double Olap2 = 0.0

    # Olap = FDSinglemodeFresnelOverlap(&apfs2, splines, &PSD, fLow, fHigh)

    Olap =  FDSinglemodeFresnelOverlap3Chan(&Afs2, &Efs2, &Tfs2, Asplines, Esplines, Tsplines, &PSDA, &PSDA, &PSDT, fLow, fHigh)


    CAmpPhaseSpline_Cleanup(Asplines)
    CAmpPhaseSpline_Cleanup(Esplines)
    CAmpPhaseSpline_Cleanup(Tsplines)

    return (Olap)
