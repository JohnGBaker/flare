import numpy as np
import os
from scipy.interpolate import interp1d


dir=os.path.dirname(os.path.realpath(__file__))
data=np.loadtxt(dir+"/MfEOBLSO.txt")
#Mf=interp1d(data[:,0],data[:,1])
Mfeta=interp1d([x/(1.0+x)**2 for x in data[:,0]],data[:,1])

def Mf(m1,m2):
    return Mfeta(m1*m2/(m1+m2)**2)/(m1+m2)/4.9254923218988636432342917247829673e-6

"""
import matplotlib.pyplot as plt
xs=np.arange(1,12,0.0031)
eta=np.array([ x/(1.0+x)**2 for x in xs])
plt.plot(1+(1-4*eta)**.5,Mf(1,1.0/xs))
#plt.plot(eta,Mf(xs)-Mfeta(eta))
plt.show()
"""

