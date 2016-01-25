

import sys
import time
import numpy as np
import math
import cmath
from scipy.interpolate import InterpolatedUnivariateSpline
from synthlisa import *

from argparse import ArgumentParser
parser = ArgumentParser()

parser.add_argument('--filehphc', action="store", required=True, help="file name for input time series for hplus, hcross")
parser.add_argument('--filetdi', action="store", default='./gentditd_synthlisa.txt', required=False, help="file name for output time series for TDI observable(s)")
parser.add_argument('--lambdapar', action="store", default='0.', required=False, help="ecliptic longitude lambda")
parser.add_argument('--betapar', action="store", default='0.', required=False, help="ecliptic latitide beta")
parser.add_argument('--psipar', action="store", default='0.', required=False, help="polarization psi")
parser.add_argument('--tagtdi', action="store", default='TDIX', required=False, help="tag selecting the TDI observable(s) to compute")
parser.add_argument('--tmin', action="store", default='nan', required=False, help="time where to start evaluating TDI (default: start of hp,hc time series)")
parser.add_argument('--tmax', action="store", default='nan', required=False, help="time where to end evaluating TDI (default: end of hp,hc time series)")

args = parser.parse_args()


#Parameters
lambdagen = float(args.lambdapar)
betagen = float(args.betapar)
psigen = float(args.psipar)
tmingen = float(args.tmin)
tmaxgen = float(args.tmax)

#Creating a LISA class - linear-in-e geometry used in KTV04, MLDC conventions
#Conventions xi0 = 3pi/2 and sw<0 should reproduce the MLDC conventions
lisa = CircularRotating(0., 3*np.pi/2, -1, 0.)

#Creating a WAVE class - loading and interpolating hp, hc from files
hphctd = np.loadtxt(args.filehphc);
hpint = InterpolatedUnivariateSpline(hphctd[:,0], hphctd[:,1], k=3)
hcint = InterpolatedUnivariateSpline(hphctd[:,0], hphctd[:,2], k=3)

wave = PyWave(hpint, hcint, betagen, lambdagen, psigen)

#Checks
#t = 0.
#wave.putwave(t)

#Creating a TDI class
tdi = TDIsignal(lisa, wave)

#Generating TDI XYZ

#Selecting time samples
times = hphctd[:,0]
imin = 0
imax = len(times)-1
if not tmingen==float('nan'):
    while times[imin]<tmingen:
        imin += 1
if not tmaxgen==float('nan'):
    while times[imax]>tmaxgen:
        imax -= 1
timesgen = times[imin:imax+1]

#tdi.Xm|Ym|Zm very slow, even if only computing a few thousands samples...
if args.tagtdi=='TDIX':
    XtdSynth = np.zeros((len(timesgen), 2))
    print 'Generating TDIX'
    for i in xrange(len(timesgen)):
        if i%1000 == 0:
            print 'Time samples computed: ' + str(i) + '/' + str(len(timesgen))
        XtdSynth[i,0] = timesgen[i]
        XtdSynth[i,1] = tdi.Xm(timesgen[i])
    np.savetxt(args.filetdi, XtdSynth)
elif args.tagtdi=='TDIY':
    YtdSynth = np.zeros((len(timesgen), 2))
    print 'Generating TDIY'
    for i in xrange(len(timesgen)):
        if i%1000 == 0:
            print 'Time samples computed: ' + str(i) + '/' + str(len(timesgen))
        YtdSynth[i,0] = timesgen[i]
        YtdSynth[i,1] = tdi.Ym(timesgen[i])
    np.savetxt(args.filetdi, YtdSynth)
elif args.tagtdi=='TDIZ':
    ZtdSynth = np.zeros((len(timesgen), 2))
    print 'Generating TDIZ'
    for i in xrange(len(timesgen)):
        if i%1000 == 0:
            print 'Time samples computed: ' + str(i) + '/' + str(len(timesgen))
        ZtdSynth[i,0] = timesgen[i]
        ZtdSynth[i,1] = tdi.Zm(timesgen[i])
    np.savetxt(args.filetdi, ZtdSynth)
elif args.tagtdi=='TDIXYZ':
    XYZtdSynth = np.zeros((len(timesgen), 4))
    print 'Generating TDIXYZ'
    for i in xrange(len(timesgen)):
        if i%1000 == 0:
            print 'Time samples computed: ' + str(i) + '/' + str(len(timesgen))
        XYZtdSynth[i,0] = timesgen[i]
        XYZtdSynth[i,1] = tdi.Xm(timesgen[i])
        XYZtdSynth[i,2] = tdi.Ym(timesgen[i])
        XYZtdSynth[i,3] = tdi.Zm(timesgen[i])
    np.savetxt(args.filetdi, XYZtdSynth)
else:
    print 'Error: TDI tag not recognized'
    exit()
