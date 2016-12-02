import numpy as np
import corner_with_covar as corner
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
        #The ordering in the params.txt file is out of whack...
        phi0=pars[3]
        lam=pars[5]
        beta=pars[6]
        inc=pars[7]
        pars[4]=phi0
        pars[5]=inc
        pars[6]=lam
        pars[7]=beta
        while(not "snr_target:" in line): line=f.readline() #Skip stuff
        snr=float(line.split()[1])
        while(not "dist_resc:" in line): line=f.readline() #Skip stuff
        pars[3]=float(line.split()[1])#replace original distance with the SNR-rescaled distance
    return pars,snr

#help(corner)

for chainfile in sys.argv[1:]:
    print "Processing results in ",chainfile
    run=re.search("Run\d*",chainfile).group(0)
    if(None==re.search("lm22",chainfile)):modes=""
    else: modes=" (2-2 only)"
    if(None==re.search("nl8k",chainfile)):res=""
    else: res=" (higher-res sampling)"
    if(None==re.search("nmm",chainfile)):mmodal=""
    else: mmodal=" (no modal decomp sampling)"
    print "run=",run
    injfile=re.sub("post_equal_weights.dat","params.txt",chainfile)
    dirname=os.path.dirname(os.path.abspath(chainfile))
    basename=os.path.basename(dirname)
    outpath=dirname+"/"+basename+"_corner.png"
    print "reading posterior samples from file:",chainfile
    print "reading injection from file:",injfile
    print "output to ",outpath
    pars,snr=read_injection(injfile)
    print "SNR=",snr
    print "pars=",pars
    Npar=len(pars)
    data=np.loadtxt(chainfile,usecols=range(Npar))
    names=[r"$m_1$",r"$m_2$",r"$t_0$",r"$D$",r"$\phi_0$",r"$\iota$",r"$\lambda$",r"$\beta$",r"$\psi$"]

    # Plot it.
    levels = 1.0 - np.exp(-0.5 * np.linspace(1.0, 3.0, num=3) ** 2)
    print "levels=",levels
    figure = corner.corner(data, bins=50,labels=names,levels=levels,
                             truths=pars,quantiles=[0.159, 0.5, 0.841],show_titles=True,range=[0.999 for x in pars],use_math_text=True,
                             title_args={"fontsize": 35},title_fmt='.2e',smooth1d=1,smooth=1,label_kwargs={"fontsize":30})
    figure.gca().annotate(run+"  SNR="+str(snr)+modes+res+mmodal, xy=(0.5, 1.0), xycoords="figure fraction",
                          xytext=(0, -5), textcoords="offset points",
                          ha="center", va="top",fontsize=30)
    for i in range(Npar):
        figure.gca().annotate(names[i]+"= %.3e"%+pars[i], xy=(0.75, 0.9-i/20.0), xycoords="figure fraction",fontsize=30)
    figure.savefig(outpath)

'''
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
'''
