import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as spline
from scipy import interpolate
import Cosmology

import GenerateFD_SignalTDIs as GenTDIFD

import LISAConstants as LC
import pyIMRPhenomD
# import phase_change_dispersion




MfCUT_PhenomD = 0.2 - 1e-7

# Compute chirp mass from m1 and m2 - same units as input m1,m2
def funcMchirpofm1m2(m1, m2):
    return pow(m1*m2, 3./5) / pow(m1+m2, 1./5)

# Newtonian estimate of the relation f(deltat) (for the 22 mode freq) - gives the starting geometric frequency for a given time to merger and chirp mass - output in Hz
# Input chirp mass (solar masses), time in years
def funcNewtonianfoft(Mchirp, t):
    if(t<=0.):
        return 0.
    return 1./np.pi * pow(Mchirp*LC.MTsun, -5./8) * pow(256.*t*LC.YRSID_SI/5, -3./8)

# Newtonian estimate of the relation deltat(f) (for the 22 mode freq) - gives the time to merger from a starting frequency for a given chirp mass - output in years
# Input chirp mass (solar masses), frequency in Hz
def funcNewtoniantoffchirp(Mchirp, f):
    return 5./256 * pow(Mchirp*MTSUN_SI, -5./3) * pow(pi*f, -8./3) / YRSID_SI

### Sylvain's functions to produce a fast MBHB signal

def funcDeltaMfPowerLaw(eta, Mf, acc):
    return 3.8 * np.power(eta * acc, 1./4.) * Mf * np.power(np.pi*Mf, 5./12.)

def funcDeltaMfLog(Mfmin, Mfmax, npt, Mf):
    return Mf * (np.power(Mfmax/Mfmin, 1/(npt - 1)) - 1)

def WaveformFrequencyGridGeom(eta, Mfmin, Mfmax, acc=1e-6):
    # Get parameters for log-sampling
    #nptlogmin = 50
    nptlogmin = 250
    deltalnMfmax = 0.025/2.0
    #deltalnMfmax = 0.05
    deltalnMf = np.fmin(deltalnMfmax, np.log(Mfmax/Mfmin)/(nptlogmin - 1))
    nptlog = 1 + np.floor(np.log(Mfmax/Mfmin) / deltalnMf)

    # Will be used to ensure the last point is not too close to fmax, which poses problems for splines (notably, derivatives become meaningless)
    minsteptoMfmax = 1.e-5*Mfmax

    # Build mixed frequency grid iteratively
    Mfreq = [Mfmin]
    while Mfreq[-1] < Mfmax:
        Mf = Mfreq[-1]
        DeltaMf = np.fmin(funcDeltaMfPowerLaw(eta, Mf, acc), funcDeltaMfLog(Mfmin, Mfmax, nptlog, Mf))
        if Mf+DeltaMf < Mfmax and Mf+DeltaMf > Mfmax-minsteptoMfmax:
            Mf_append = Mfmax
        else:
            Mf_append = np.fmin(Mfmax, Mf+DeltaMf)
        Mfreq.append(Mf_append)

    # Output
    return np.array(Mfreq)

def ResponseFrequencyGrid(wf, nptlogmin=250, Deltat_max=0.02083335, Deltaf_max=0.001):
    # Min and max physical frequencies (Hz)
    freq_PhD, amp_PhD, phase_PhD = wf
    fmin, fmax = freq_PhD[0], freq_PhD[-1]

    # Sampling target for log
    Deltalnf_max = np.log(fmax/fmin)/(nptlogmin - 1)

    # Build tf function and f(t) by interpolation of the monotonous part
    tfspline = spline(freq_PhD, 1./(2.*np.pi) * phase_PhD).derivative()
    tf = tfspline(freq_PhD)
    index_cuttf = 0
    tfdiff = np.diff(tf)
    while index_cuttf<len(tfdiff)-1 and tfdiff[index_cuttf]>0:
        index_cuttf += 1
    tfr = tf[:index_cuttf+1]
    freq_PhDr = freq_PhD[:index_cuttf+1]
    ftspline = spline(tfr, freq_PhDr)

    # Will be used to ensure the last point is not too close to fmax, which poses problems for splines (notably, derivatives become meaningless)
    minsteptofmax = 1e-5*fmax

    # Build mixed frequency grid iteratively
    freq = [fmin]
    while freq[-1] < fmax:
        f = freq[-1]
        t_new = tfspline(f) + Deltat_max * LC.YRSID_SI
        if tfr[0]<=t_new and t_new<=tfr[-1]:
            Deltaf_from_time = ftspline(t_new) - f
        else:
            Deltaf_from_time = fmax - f
        Deltaf = np.min([Deltaf_max, f * Deltalnf_max, Deltaf_from_time])
        if f+Deltaf < fmax and f+Deltaf > fmax-minsteptofmax:
            f_append = fmax
        else:
            f_append = np.fmin(fmax, f+Deltaf)
        freq.append(f_append)

    # Output
    return np.array(freq)


### Auxiliary functions


AU = LC.ua
LL = 2.5e9
Reff = 2./np.pi*AU
Omega0 = 2.*np.pi/LC.YRSID_SI
armt = LL/LC.clight

def funcestimateforfomorb2nd1(f):
    return (2.*np.pi*f*Reff/LC.clight*Omega0**2)
def funcestimateforfomorb2nd2(f):
    return (2.*np.pi*f*Reff/LC.clight*Omega0)**2
def funcestimateforfomorb2nd(f):
    return max(funcestimateforfomorb2nd1(f), funcestimateforfomorb2nd2(f))
def funcestimateforfomconst2nd1(f):
    return (4.*np.pi*f*armt*Omega0**2)
def funcestimateforfomconst2nd2(f):
    return (4.*np.pi*f**2*(armt)**2*Omega0**2)
def funcestimateforfomconst2nd3(f):
    return Omega0**2
def funcestimateforfomconst2nd(f):
    return max(funcestimateforfomconst2nd1(f), funcestimateforfomconst2nd2(f), funcestimateforfomconst2nd3(f))


def funcFOMApprox(key, m1, m2, fstart):
    M = m1 + m2
    q = max(m1/m2, m2/m1)
    nu = q/(1+q)**2
    Msec = M*LC.MTsun
    def Tffunc(f):
        return 1./8*np.sqrt(5./(3*nu))*Msec*(np.pi*Msec*f)**(-11./6)
    if key=='Psiorb':
        fom = 1./2*Tffunc(fstart)**2*funcestimateforfomorb2nd(fstart)
    elif key=='Psiconst':
        fom = 1./2*Tffunc(fstart)**2*funcestimateforfomconst2nd(fstart)
    return fom


def ComputeLISASOBBH_FDTDI(p, fmax, Tobs, dt, buf=None, order=None):
    """
    This function determines the order of frensel stencil and computes SOBBH GW signal for a single source
    @param p parameter dictionary as reurned by hdf5 reading function
    @param Tobs observation time (sec)
    @param buf FrequencyArray where the computed signal will be placed
    @param order user defined order of stencil for FD response (forcing order)

    output
    frequency array
    buf FrequencyArray (3xNf) containing X, Y, Z in FD
    wfTDI debugging info dictionary containg Response functions, GW phase, amplitude, dopper delay
    """
    m1 = p['m1']
    m2 = p['m2']
    z = p['z']


    DL = p['DL']
    #print ("DL = ", DL, pars[3])

    bet = p['bet']
    lam = p['lam']
    inc = p['incl']
    psi = p['psi']

    chi1 = p['a1']
    chi2 = p['a2']

    phi0 = p['phi0']
    fstart = p['f0']
    dip_ampl=p['dip_ampl']
    wl=p['wl']
    disp_order=p['disp_order']

    if (m2>m1):
        m_tmp = m1
        m1 = m2
        m2 = m_tmp
    q = m1/m2
    eta = m1*m2/(m1+m2)**2
    Mc = funcMchirpofm1m2(m1, m2)
    fRef = fstart
    tRef = 0.0
    fny = 0.5/dt

    fend = min(fny, fmax)

    ## determine the order of stencil
    fom = funcFOMApprox('Psiorb', m1, m2, fstart)
    if fom<0.1:
        order = 3
    elif fom<0.2:
        order = 5
    elif fom<0.5:
        order=10
    else:
        order=20




    # in CI
    m1_SI = m1*LC.MsunKG
    m2_SI = m2*LC.MsunKG
    dist_SI = DL*1e6*LC.pc

    Ms=(m1+m2)*LC.MTsun#*solar mass in sec

    df = 1.0/Tobs
    acc_sampling = 1.e-5  ## hardcoded tolerance for the interpolation
    freq_PhD = 1/Ms * WaveformFrequencyGridGeom(eta, Ms*fstart, Ms*fend, acc=acc_sampling)

    ### Generating the h22 (amplitude and phase) in SSB
    wf_PhD_class = pyIMRPhenomD.IMRPhenomDh22AmpPhase(freq_PhD, phi0, fRef, m1_SI, m2_SI, chi1, chi2, dist_SI)
    freq_PhD, amp_PhD, phase_PhD = wf_PhD_class.GetWaveform()


    # if disp_order==1:
    #     phase_PhD = phase_PhD-(3./224.)*pow(eta,2/5)*dip_ampl*pow(np.pi*Mc*LC.MTsun*freq_PhD,-7./3.)\
    #     -phase_change_dispersion.reduced_disp_shift_coeff(wl,disp_order,DL,Mc)*np.log(np.pi*Mc*LC.MTsun*freq_PhD)\
    #
    # else:
    #     phase_PhD = phase_PhD-(3./224.)*pow(eta,2/5)*dip_ampl*pow(np.pi*Mc*LC.MTsun*freq_PhD,-7./3.)\
    #     -phase_change_dispersion.reduced_disp_shift_coeff(wl,disp_order,DL,Mc)*pow(np.pi*Mc*LC.MTsun*freq_PhD,disp_order-1)\


    tfspline = spline(freq_PhD, 1./(2.*np.pi) * phase_PhD).derivative()
    tf = tfspline(freq_PhD)
    tf = tf - tf[0]  ### t=0 should correspond to fstart


    # index_cuttf = 0
    # tfdiff = np.diff(tf)
    # while index_cuttf<len(tfdiff)-1 and tfdiff[index_cuttf]>0:
    #     index_cuttf += 1
    # tfr = tf[:index_cuttf+1]
    # # print ("cutoff:", index_cuttf, len(tf), tf[index_cuttf], tf[-1])
    # freq_PhDr = freq_PhD[:index_cuttf+1]
    # ftspline = spline(tfr, freq_PhDr)

    fr_end = fend
    # if the duration of the signal >Tobs, we chop it off
    ind = len(tf)-1
    if (tf[-1] > Tobs):
        ind = np.argwhere(tf > Tobs)[0][0]
        frspl = spline(tf[:ind], freq_PhD[:ind])
        fr_end = frspl(Tobs)
        ### TODO: DO I need to rederive the freq. grid again? Maybe not...

    freq_rs, amp_rs, phase_rs = np.copy(freq_PhD[:ind]), np.copy(amp_PhD[:ind]), np.copy(phase_PhD[:ind])


    settRefAtfRef = True
    if settRefAtfRef:
        fRef_inrange = min(max(freq_rs[0], fRef), freq_rs[-1])
        phase_rs += 2.*np.pi *(tRef - tfspline(fRef_inrange))*(freq_rs - fRef_inrange)  #shape problem??


    #### Computing set of frequencies  for the response
    freq_response = ResponseFrequencyGrid([freq_rs, amp_rs, phase_rs])

    #print ("SB: Return me this number", len(freq_response))

    tm_response = tfspline(freq_response)-tf[0]


    ### needed for high order
    dffphasespline = spline(freq_PhD, phase_PhD).derivative(2)  #phD or rs?
    epsTf2vec = 1./(4*np.pi*np.pi)*dffphasespline(freq_response)
    Tfvec = np.sqrt(np.abs(epsTf2vec))
    epsTfvec = np.sign(epsTf2vec)
    ### FIXME debugging!
    order = 0
    wfTDI =  GenTDIFD.JustLISAFDresponseTDI(freq_response, tm_response, Tfvec, epsTfvec, inc, lam, bet, psi, t0=0.0, \
                    order_fresnel_stencil=order)

    return(freq_rs, amp_rs, phase_rs, freq_response, wfTDI)


    #
    # df = 1.0/Tobs
    # Nf = int(fny/df)+1
    # # fr_F =  np.arange(Nf)*df ### Fourier Frequency array
    # # fr_F =  np.linspace(0, Nf-1, Nf)*df ### Fourier Frequency array
    # fr_F =  np.zeros(Nf)*df ### Fourier Frequency array
    #
    #
    # i_st = int(np.ceil(freq_rs[0]/df))  ### start index
    # i_en = int(np.floor(freq_rs[-1]/df))  ### end index
    #
    # fr_wvf = fr_F[i_st:i_en]
    #
    #
    # ### Interpolating:
    # # amp_spl = spline(freq_rs, amp_rs)
    # # ph_spl = spline(freq_rs, phase_rs)
    # #
    # # amp_wvf = amp_spl(fr_wvf)
    # # ph_wvf = ph_spl(fr_wvf)
    # #
    # # phaseRdelayspline = spline(freq_response, wfTDI['phaseRdelay'])
    # # phaseRdelay = phaseRdelayspline(fr_wvf)
    # #
    # # phasetimeshift = 2.*np.pi*tRef*fr_wvf
    # # fastPart = amp_wvf * np.exp(1j*(ph_wvf + phaseRdelay + phasetimeshift))
    # #
    # # keytrs = ['transferL1', 'transferL2', 'transferL3']
    # #
    # # for I, ky in enumerate(keytrs):
    # #     transferLRespline = spline(freq_response, np.real(wfTDI[ky]))
    # #     transferLImspline = spline(freq_response, np.imag(wfTDI[ky]))
    # #     transferLRe = transferLRespline(fr_wvf)
    # #     transferLIm = transferLImspline(fr_wvf)
    # #     if (ky == 'transferL1'):
    # #         X = np.conjugate((transferLRe+1j*transferLIm) * fastPart)
    # #         # X = np.conjugate(fastPart)
    # #     if (ky == 'transferL2'):
    # #         Y = np.conjugate((transferLRe+1j*transferLIm) * fastPart)
    # #     if (ky == 'transferL3'):
    # #         Z = np.conjugate((transferLRe+1j*transferLIm) * fastPart)
    #
    # ### TODO: Should I use FrequencyArray???
    # Xfull = np.zeros(Nf, dtype=np.complex128)
    # Yfull = np.zeros(Nf, dtype=np.complex128)
    # Zfull = np.zeros(Nf, dtype=np.complex128)
    #
    # # Xfull[i_st:i_en] = X
    # # Yfull[i_st:i_en] = Y
    # # Zfull[i_st:i_en] = Z
    # ### Output: Full Fourier frequency array and corresponding waveform
    #
    # #print(fr_F[i_en+20000])
    # return (fr_F, Xfull, Yfull, Zfull,i_st,i_en)


# pr_gr={'z': 0.055607, 'm1': 40.018, 'm2': 30.426, 'f0': 0.0124765, 'lam': 3.5064, 'bet': 0.17769632679489655, 'incl': 2.7721, 'psi': 1.0922, \
# 'phi0': 5.3862, 'DL': 258.7520806702613,'a1': 0.05420918130401683, 'a2': 0.02390622153705309, 'dip_ampl': 0., 'disp_order': 0., 'wl':0.}
#
# pr_dip={'z': 0.055607, 'm1': 40.018, 'm2': 30.426, 'f0': 0.0124765, 'lam': 3.5064, 'bet': 0.17769632679489655, 'incl': 2.7721, 'psi': 1.0922,\
#   'phi0': 5.3862, 'DL': 258.7520806702613, 'a1': 0.05420918130401683, 'a2': 0.02390622153705309,'dip_ampl':-1E-5,'disp_order': 0., 'wl':0.}
#
# pr_disp={'z': 0.055607, 'm1': 40.018, 'm2': 30.426, 'f0': 0.0124765, 'lam': 3.5064, 'bet': 0.17769632679489655, 'incl': 2.7721, 'psi': 1.0922,\
#   'phi0': 5.3862, 'DL': 258.7520806702613, 'a1': 0.05420918130401683, 'a2': 0.02390622153705309,'dip_ampl':0.,'disp_order': 0, 'wl':1E15}
#
#
# # pr['z'] = 0.055607
# # pr['m1'] = 40.018
# # pr['m2'] = 30.426
# # pr['f0'] = 1.e-2
# # pr['lam'] = 1.0
# # pr['bet'] = -0.3
# # pr['incl'] = 0.45
# # pr['psi'] = 1.7
# # chi1 = 0.7
# # chi2 = 0.4
# # pr['phi0'] = 2.1
# # pr['DL'] = 240.0
# # pr['a1'] = chi1
# # pr['a2'] = chi2
#
# Tobs = 167772160.0
# dt = 5.
# fmax = 0.018396210670471192
#
#
# fr_gr, Xf_gr, Yf_gr, Zf_gr, ist_gr, ien_gr = ComputeLISASOBBH_FDTDI(pr_gr, fmax, Tobs, dt)
#
# fr_dip, Xf_dip, Yf_dip, Zf_dip, i_st_dip, ien_dip = ComputeLISASOBBH_FDTDI(pr_dip, fmax, Tobs, dt)
#
# fr_disp, Xf_disp, Yf_disp, Zf_disp, ist_disp, iend_disp = ComputeLISASOBBH_FDTDI(pr_disp, fmax, Tobs, dt)
#
# # f_old=open("old_X.dat","r")
# #
# # old_data=np.genfromtxt(f_old,comments='#')
# #
# # fr_old=old_data[:,0]
# # Re_X_old=old_data[:,1]
#
# plt.figure()
# plt.plot(fr_gr,np.real(Xf_gr))
# plt.plot(fr_dip,np.real(Xf_dip))
#
# plt.show()
