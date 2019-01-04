import numpy as np
import sys, os, re
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


######################## main #####################
#### we randomly choose points around signal in 2-D parameter space (fixing other params)
################


system={'z': 0.055607, 'm1': 40.018, 'm2': 30.426, 'f0': 0.0124765, 'lam': 3.5064, 'bet': 0.17769632679489655, \
        'incl': 2.7721, 'psi': 1.0922, 'phi0': 5.3862, 'DL': 258.7520806702613, 'a1': 0.05420918130401683, \
        'a2': 0.02390622153705309,'dip_ampl': 0., 'disp_order': 0, 'wl':0.}


Tobs = 167772160.0
df=1/Tobs
dt = 5.
### FIXME where did it come from?/
fmax_g = 0.018396210670471192

frqs, phase, ampls = GetSignal(system, fmax_g, obj="AET")
wvf1 = {}
wvf1["freq"] = frqs
# wvf1["ampl"] = ampls[0]
wvf1["ampl"] = np.array(ampls)
wvf1["phase"] = phase


templ = system.copy()
m1 = templ['m1']
m2 = templ['m2']
Mc_tr = SOBBH.funcMchirpofm1m2(m1, m2)
q_tr = m1/m2
Mc_min = Mc_tr*(1.0 - 0.00001)
Mc_max = Mc_tr*(1.0 + 0.00001)
Mc_min = 30.32
Mc_max = 30.3203

q_min = q_tr*(1.0 - 0.002)
q_max = q_tr*(1.0 + 0.002)
if (q_min < 1.0):
    q_min = 1.0
q_min = 1.3
q_max = 1.33

print ("chirp mass", Mc_tr, Mc_min, Mc_max)
print ("mass ratio", q_tr, q_min, q_max)

Ntrial = 10000
np.random.seed(3)
Mcs = Mc_min  + np.random.uniform(0.0, 1.0, size=Ntrial)*(Mc_max - Mc_min)
qs = q_min + np.random.uniform(0.0, 1.0, size=Ntrial)*(q_max - q_min)

spr = "   "
fout = open("Overlaps_McqChan3.dat", 'a')
# olap = flar.OverlapFrensel1(wvf1, wvf1, fmin=frqs[0], fmax=frqs[-1], noise="AE")
olap = flar.OverlapFrensel3(wvf1, wvf1, fmin=frqs[0], fmax=frqs[-1], noise="AE")
norm = olap
fout.write(str(Mc_tr) + spr + str(q_tr) + spr + str(olap - 0.5*norm) + "  \n")
st = time.time()
for i in range(Ntrial):
    Mc = Mcs[i]
    q = qs[i]
    m1=0.5*Mc*((q/((q+1.)**2))**(-3.0/5.0) + ((q/((q+1.)**2))**(-6.0/5.0) - 4.0*(q/((q+1.)**2))**(-1.0/5.0))**(0.5))
    m2=0.5*Mc*((q/((q+1.)**2))**(-3.0/5.0) - ((q/((q+1.)**2))**(-6.0/5.0) - 4.0*(q/((q+1.)**2))**(-1.0/5.0))**(0.5))
    # print ("check:", Mc, SOBBH.funcMchirpofm1m2(m1, m2), q, m1/m2)
    templ["m1"] = m1
    templ["m2"] = m2
    frqs2, phase2, ampls2 = GetSignal(templ, fmax_g, obj="AET")
    wvf2 = {}
    wvf2["freq"] = frqs2
    # wvf2["ampl"] = ampls2[0]
    wvf2["ampl"] = np.array(ampls2)
    wvf2["phase"] = phase2
    # print (system)
    # print (templ)
    fLow = max(frqs[0], frqs2[0])
    fHigh = min(frqs[-1], frqs2[-1])

    # olap = flar.OverlapFrensel1(wvf1, wvf2, fmin=fLow, fmax=fHigh, noise="AE")
    # norm = flar.OverlapFrensel1(wvf2, wvf2, fmin=frqs2[0], fmax=frqs2[-1], noise="AE")
    olap = flar.OverlapFrensel3(wvf1, wvf2, fmin=fLow, fmax=fHigh, noise="AE")
    norm = flar.OverlapFrensel3(wvf2, wvf2, fmin=frqs2[0], fmax=frqs2[-1], noise="AE")
    fout.write(str(Mc) + spr + str(q) + spr + str(olap - 0.5*norm) + "  \n")

    # sys.exit(0)
en  = time.time()
print ("ellapsed:", en-st, (en-st)/Ntrial)

fout.close()
