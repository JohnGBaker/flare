import numpy as np
import sys
import flar

from scipy.interpolate import InterpolatedUnivariateSpline as spline
from scipy import interpolate


import matplotlib
import matplotlib.pyplot as plt


import tdi
import SOBBHtemplate as SOBBH
import time
import cProfile

from scipy import integrate


def GetSignal(p, fmax, obj="XYZ"):
    frW, ampW, phaseW, freq_resp, wfTDI = SOBBH.ComputeLISASOBBH_FDTDI(p, fmax, Tobs, dt)# buf=None, order=None)

    # print ("waveform freqs:", frW)
    # phR = wfTDI['phaseRdelay']
    # RX = wfTDI['transferL1']
    # RY = wfTDI['transferL2']
    # RZ = wfTDI['transferL3']
    # print ("resp freqs", freq_resp)

    # phasetimeshift = 2.*np.pi*tRef*fr_wvf
    # fastPart = amp_wvf * np.exp(1j*(ph_wvf + phaseRdelay + phasetimeshift))

    # print ("sizes Wv, Resp", len(frW), len(freq_resp))

    # plt.plot(frW, frW, 'o')
    # plt.plot(freq_resp, freq_resp, 'o')
    # plt.grid(True)
    # plt.show()

    frqs = np.unique(np.concatenate((frW, freq_resp)))
    N = len(frqs)
    # print (len(frqs))
    # print (np.diff(frqs))

    amp_spl = spline(frW, ampW)
    ph_spl = spline(frW, phaseW)
    phRspl = spline(freq_resp, wfTDI['phaseRdelay'])
    amp = amp_spl(frqs)
    ph = ph_spl(frqs) + phRspl(frqs)

    keytrs = ['transferL1', 'transferL2', 'transferL3']

    for I, ky in enumerate(keytrs):
        transferLRespline = spline(freq_resp, np.real(wfTDI[ky]))
        transferLImspline = spline(freq_resp, np.imag(wfTDI[ky]))
        transferLRe = transferLRespline(frqs)
        transferLIm = transferLImspline(frqs)
        if (ky == 'transferL1'):
                X = amp*(transferLRe+1j*transferLIm)
        # #         # X = np.conjugate(fastPart)
        if (ky == 'transferL2'):
                Y = amp*(transferLRe+1j*transferLIm)
        if (ky == 'transferL3'):
                Z = amp*(transferLRe+1j*transferLIm)

    if (obj == "XYZ"):
        return(frqs, ph, [X, Y, Z])
    elif (obj == "AET"):
        A = (Z - X)/np.sqrt(2.0)
        E = (X - 2.0*Y + Z)/np.sqrt(6.0)
        T = (X+Y+Z)/np.sqrt(3.0)
        return (frqs, ph, [A, E, T])
    else:
        raise ValueError





# phaseRdelayspline = spline(freq_response, wfTDI['phaseRdelay'])
# # phaseRdelay = phaseRdelayspline(fr_wvf)
# #
# # phasetimeshift = 2.*np.pi*tRef*fr_wvf
# # fastPart = amp_wvf * np.exp(1j*(ph_wvf + phaseRdelay + phasetimeshift))
# #


system={'z': 0.055607, 'm1': 40.018, 'm2': 30.426, 'f0': 0.0124765, 'lam': 3.5064, 'bet': 0.17769632679489655, \
        'incl': 2.7721, 'psi': 1.0922, 'phi0': 5.3862, 'DL': 258.7520806702613, 'a1': 0.05420918130401683, \
        'a2': 0.02390622153705309,'dip_ampl': 0., 'disp_order': 0, 'wl':0.}


Tobs = 167772160.0
df=1/Tobs
dt = 5.
### FIXME where did it come from?/
fmax_g = 0.018396210670471192

Ntrial = 1
st = time.time()
for i in range(Ntrial):
    frqs, phase, ampls = GetSignal(system, fmax_g, obj="AET")
en = time.time()
print ("ellapsed:", en-st, (en-st)/Ntrial)

#Compute  overlap

wvf1 = {}
wvf1["freq"] = frqs
wvf1["ampl"] = ampls[0]
wvf1["phase"] = phase

wvf2 = {}
wvf2["freq"] = frqs
wvf2["ampl"] = ampls[0]
wvf2["phase"] = phase

# fg, ax = plt.subplots(nrows =2, sharex=True)
# ax[0].plot(frqs, np.real(ampls[0]))
# ax[1].plot(frqs, np.imag(ampls[0]))
# plt.show()
#
# plt.plot(frqs, np.absolute(ampls[0]))
# plt.grid(True)
# plt.show()

#### Compute Overlap
SA = tdi.noisepsd_AE(frqs)
# AmpSpl = spline(frqs, np.absolute(ampls[0])**2/SA)
AmpSpl = spline(frqs, np.absolute(ampls[0])**2)

def Intgrnd(x):
    return(AmpSpl(x))


out = integrate.quad(Intgrnd, frqs[0], frqs[-1])

print ("Olap = ", 4.0*out[0])

#Npt = int( (frqs[-1] - frqs[0]) )

# res = 0.0
# f = frqs[0]
# while (f < frqs[-1]):
#     res += AmpSpl(f)
#     f = f+ df
#
# print ("again", 4.0*res*df)


Ntrial = 100
st = time.time()
for i in range(Ntrial):
    frqs, phase, ampls = GetSignal(system, fmax_g, obj="AET")
    wvf1 = {}
    wvf1["freq"] = frqs
    wvf1["ampl"] = ampls[0]
    wvf1["phase"] = phase

    wvf2 = {}
    wvf2["freq"] = frqs
    wvf2["ampl"] = ampls[0]
    wvf2["phase"] = phase

    olap = flar.OverlapFrensel1(wvf1, wvf2, fmin=frqs[0], fmax=frqs[-1], noise="AE")
en = time.time()
print ("ellapsed:", en-st, (en-st)/Ntrial)



print ("Olap = ", olap)

st = time.time()
for i in range(Ntrial):
    frqs, phase, ampls = GetSignal(system, fmax_g, obj="AET")
    wvf1 = {}
    wvf1["freq"] = frqs
    wvf1["ampl"] = np.array(ampls)
    wvf1["phase"] = phase

    wvf2 = {}
    wvf2["freq"] = frqs
    wvf2["ampl"] = np.array(ampls)
    wvf2["phase"] = phase

    # print (np.shape(wvf1["ampl"]), np.shape(wvf2["ampl"]))

    olap = flar.OverlapFrensel3(wvf1, wvf2, fmin=frqs[0], fmax=frqs[-1], noise="AE")
    # print (olap)
    # sys.exit(0)
en = time.time()
print ("ellapsed:", en-st, (en-st)/Ntrial)


# dat = np.genfromtxt("TestOlap1.dat")
# print (np.shape(dat))
# dat2 = np.genfromtxt("TestOlap2.dat")
# print (np.shape(dat2))
#
# fr = np.linspace(frqs[0], frqs[-1], num=100000)
# Intgr = Intgrnd(fr)
#
# # print ("Here 1 \n", ampls[0])
#
#
# # print (dat[:, 2])
# # plt.plot(frqs, np.real(ampls[0])**2 + np.imag(ampls[0])**2, '--')
# plt.plot(frqs, np.real(ampls[0]))
# # plt.plot(fr, Intgr)
# plt.plot(dat[:, 0], dat[:, 1], '--')
# # plt.plot(dat2[:, 0], dat2[:, 1], '--')
# plt.plot(dat2[:, 0], dat2[:, 2], '--')
# # plt.plot(frqs, np.real(ampls[0])**2 + np.imag(ampls[0])**2, '--')
# plt.grid(True)
# plt.show()
