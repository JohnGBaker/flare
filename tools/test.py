import numpy as np
import flar

import matplotlib
import matplotlib.pyplot as plt

import tdi


# print ("Bessel:", Cy_gsl.g_bess(4.0))
#
a = np.ones(5)
print (a)

b = flar.Arr_Vec2(a)
print ("B = ", b)
#
# b = Cy_gsl.Arr_Vec2(a)
# print ("B2 = ", b)

# print ("complex:", Cy_gsl.g_complex(2.3, 4.5))

N=2
M=3
a = np.arange(N*M).reshape((N,M))
a = a*0.1
print ("a = ", a)

print ("out = ", flar.Arr_Matr(a))


x = np.arange(200)*0.01
y = np.sin(10.*x)

# plt.plot(x, y)
# plt.grid(True)
# plt.show()

cf = flar.BuildNotAKspline(x, y)

print (cf[:10, :])

cfQ = flar.BuildQuadspline(x, y)

print (cfQ[:10, :])

fr = np.logspace(-5., -1., num=1000)
SX_ldc = tdi.noisepsd_X(fr)
SA_ldc = tdi.noisepsd_AE(fr)
ST_ldc = tdi.noisepsd_T(fr)
Nf = len(fr)
SX = np.zeros(Nf)
S_ = np.zeros(Nf)
SA = np.zeros(Nf)
SE = np.zeros(Nf)
ST = np.zeros(Nf)

for i in range(Nf):
    SX[i], SA[i], SE[i], ST[i] = flar.GetNoiseXAET(fr[i])
    S_[i] = flar.GetNoise(fr[i], noise="X")


plt.loglog(fr, S_)
plt.loglog(fr, SX_ldc, "--")
plt.loglog(fr, SA)
plt.loglog(fr, SA_ldc, "--")
plt.loglog(fr, ST)
plt.loglog(fr, ST_ldc, "--")
plt.grid(True)
plt.show()
