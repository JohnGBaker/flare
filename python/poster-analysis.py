import numpy as np
import triangle
import matplotlib.pyplot as plt
#from matplotlib.ticker import MultipleLocator
#from scipy.interpolate import interp1d
#from scipy.ndimage import gaussian_filter1d
import sys
import os.path
import re

def read_injection(file):
    pars=[]
    with open(file,'r') as f:
        line="#"
        while(not "m1:" in line): line=f.readline() #Skip comment
        while(not "---" in line): 
            val=float(line.split()[1])
            print "par val=",val
            pars.append(val)
            line=f.readline() #Skip comment
    return pars

chainfile=sys.argv[1]
injfile=re.sub("post_equal_weights",chainfile,"params")
dirname=os.path.dirname(os.path.abspath(chainfile))
basename=os.path.basename(dirname)
outpath=dirname+"/"+basename+"_corner.png"
print "reading posterior samples from file:",chainfile
print "reading injection from file:",injfile
print "output to ",outpath

pars=read_injection(injfile)
print "pars=",pars
Npar=len(pars)
data=np.loadtxt(file,usecols=range(9))
#injection=[1.5e6, 0.5e6,0, 1.04822e6, 0,1.31,1.7, 1.0471976,1.2 ]
#injection=[1.5e6, 0.5e6,0, 2.4389e6, 0,0,1.7, 1.0471976,1.2 ]
names=[r"$m_1$",r"$m_2$",r"$t_0$",r"$D$",r"$\phi_0$",r"$\iota$",r"$\lambda$",r"$\beta$",r"$\psi$"]

# Plot it.
levels = 1.0 - np.exp(-0.5 * np.linspace(1.0, 3.0, num=3) ** 2)
figure = triangle.corner(data, bins=40,labels=names,
                         truths=injection,
                         quantiles=[0.159, 0.5, 0.841],
                         show_titles=True, title_args={"fontsize": 12})
figure.gca().annotate(basename, xy=(0.5, 1.0), xycoords="figure fraction",
                      xytext=(0, -5), textcoords="offset points",
                      ha="center", va="top")
figure.savefig(outpath)

maxlike = {'NN': [.416,1.875,-.483,3.418,3421], 'RF': [.48,1.7,-5.934,6.857,4455],
    'AB': [.489,1.681,-5.95,6.682,4392]}
#maxpost = {'NN': [.487,1.677,-.516,8.855,5756], 'RF': [.316,2.132,-.677,5.049,5552],
#    'AB': [.756,1.224,-1.155,8.196,3974]}
for model in ['NN','RF','AB']:
    datafile = 'chains/realdata_varyz1_'+model+'_post_equal_weights.dat'
    data = np.loadtxt(datafile, usecols=(0,1,2,3,5))
    figure = triangle.corner(data, labels=[r'$n_0$',r'$n_1$',r'$n_2$', r'$z_1$', r'$N_{\rm tot}$'],
                             bins=50, quantiles=[0.05, 0.5, 0.95], show_titles=True, levels=levels,
                             title_args={"fontsize": 14}, verbose=False, smooth1d=1, smooth=1,
                             range=[(0,1.6),(0.7,3.2),(-6,0),(1.,10),(1500,10000)],
                             maxlike=maxlike[model], label_kwargs={"fontsize": 20})
    figure.savefig('../../paper/figures/realdata_varyz1_posterior_'+model+'.png')
    plt.close(figure)
