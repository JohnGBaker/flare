import numpy as np
import pylab as pl
import flare


inj_params = dict(mass1=2e6, mass2=1e6, distance=10e3,
                  beta=-0.073604, inclination=0.6459, polarization=1.7436,
                  n_modes=5, f_low=1e-3, n_samples=32768, log_samples=True)
inj_params['lambda'] = 3.4435 # doh
injection = flare.generate_injection_reim(**inj_params)

# calculate the log-likelihood for a bunch
# of random masses around the injected masses
m1s_ = np.random.normal(inj_params['mass1'], 50, 100)
m2s_ = np.random.normal(inj_params['mass2'], 50, 100)
m1s = np.where(m1s_ > m2s_, m1s_, m2s_)
m2s = np.where(m1s_ > m2s_, m2s_, m1s_)
lls = []
for m1, m2 in zip(m1s, m2s):
    template_params = inj_params.copy()
    template_params['mass1'] = m1
    template_params['mass2'] = m2
    template_params['injection'] = injection
    lls.append(flare.calculate_logl_reim(**template_params))

# plot the likelihood
lls = np.array(lls)
sorter = np.argsort(lls)
pl.tricontourf(m1s[sorter], m2s[sorter], np.exp(lls[sorter]), 100, cmap='gnuplot2_r')
pl.colorbar()
pl.axvline(inj_params['mass1'], color='#00ff00')
pl.axhline(inj_params['mass2'], color='#00ff00')
pl.xlabel('Mass 1')
pl.ylabel('Mass 2')
pl.show()
