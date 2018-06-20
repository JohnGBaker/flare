#Usage:
#python summer17analysis.py chainfile_1 chainfile_2 ... chainfile_N fisherfile_1 fisherfile_2 ... fisherfile_N

import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import corner_with_covar as corner
import matplotlib.pyplot as plt
#from matplotlib.ticker import MultipleLocator
#from scipy.interpolate import interp1d
#from scipy.ndimage import gaussian_filter1d
import sys
import os.path
import re
import astropy.units as units
#import sky_area.sky_area_clustering as skyarea
#import sets 
import math
import random
import argparse

def readCovar(file):
    pars=[]
    done=False
    trycount=0
    with open(file,'r') as f:
        line="#"
        while("#" in line): line=f.readline() #Skip comment
        for val in line.split():
            #print val
            pars.append(float(val))
        Npar=len(pars)
        while(not "#Covariance" in line):line=f.readline() #Skip until the good stuff
        covar=np.zeros((Npar,Npar))
        i=0
        for par in pars:
            line=f.readline()
            print(i,":",line)
            covar[i]=np.array(line.split())
            i+=1
        print("done")
    return covar

def read_snr(file):
    snr=0
    with open(file,'r') as f:
        line="\n"
        while(not line=="" and not ("SNR" in line)): line=f.readline() #Skip stuff
        if(not line==""):
            words=line.split()
            for i in range(len(words)-1):
                if("SNR" in words[i]):snr=float(words[i+1])
    return snr

def read_injection(file):
    pars=[]
    with open(file,'r') as f:
        line="#"
        while(not "m1:" in line): line=f.readline() #Skip comment
        while(not "---" in line): 
            val=float(line.split()[1])
            print ("par val=",val)
            pars.append(val)
            line=f.readline() #Skip comment
        #The ordering in the params.txt file is out of whack...
        phi0=pars[3]
        lam=pars[5]
        beta=pars[6]
        inc=pars[7]
        distance=pars[4]
        pars[4]=phi0
        pars[5]=inc
        pars[6]=lam
        pars[7]=beta
        while(not "snr_target:" in line): line=f.readline() #Skip stuff
        snr=0
        if(not re.search("nan",line)):
            snr=float(line.split()[1])
            while(not "dist_resc:" in line): line=f.readline() #Skip stuff
            distance=float(line.split()[1])#replace original distance with the SNR-rescaled distance
        pars[3]=distance
    return pars,snr

def read_samples(pars,sample_file,code="bambi",burnfrac=0,keeplen=-1,iskip=1,maxsamps=30000):
    if(code=="bambi"):
        print ("Reading BAMBI samples")
        Npar=len(pars)
        data=np.loadtxt(sample_file,usecols=range(Npar))
    else:
        print ("Reading PTMCMC samples")
        Npar=len(pars)
        #data=np.loadtxt(chainfile,usecols=range(5,5+Npar))
        data=np.loadtxt(sample_file,usecols=[0]+list(range(5,5+Npar)))
        print("data[-1]=",data[-1],"data[-2]=",data[-2][0])
        every=data[-1][0]-data[-2][0]
        data=data[:,1:]
        if(keeplen/iskip>maxsamps):iskip=int(keeplen/maxsamps)
        iev=int(iskip/every)
        if(iev>1):
            data=data[::iev]
            every*=iev
        print ("every=",every," iev=",iev," iskip=",iskip)            
        #print "shape=",data.shape
        keeplen=int(keeplen/every)

    print ("code=",code,"  keeplen=",keeplen,"  burnfrac=",burnfrac,"  len=",len(data) )

    if(keeplen>0):
        data=data[len(data)-keeplen:]
    else:
        data=data[int(len(data)*burnfrac):]

    if(code!="bambi"):
        outfile=re.sub("_t0.dat","_post_samples.dat",sample_file)
        print("Writing PTMCMC samples to '",outfile,"'")
        np.savetxt(outfile,data)
        
    return data

def fake_samples(pars,cov,n=1000):
    return np.random.multivariate_normal(pars,cov,n)


def cornerFisher(pars,samples,in_cov,cred_levs,iparmin=0,iparend=9):
    data=[]
    for datum in samples:
    #Sylvain:Here crop down to a limited parameter set for a smaller plot
    #the overlaid text is added with the "annotate" commands below
        #print ("data shape=",datum.shape)
        istart=iparmin;
        iend=iparend
        data.append(datum[:,istart:iend])
    pars=pars[istart:iend]
    names=parnames[istart:iend]
    if(isinstance(in_cov,list)):
        cov=[]
        for ic in in_cov:
            cov.append(ic[istart:iend,istart:iend])
    elif in_cov is not None:
        cov=in_cov[istart:iend,istart:iend]
    else:
        cov=None
    Npar=len(pars)
    fontscale=0.1+Npar/9.0
    print ("pars=",pars)
    #print "cov=",cov
    # Plot it.
    #levels = 1.0 - np.exp(-0.5 * np.linspace(1.0, 3.0, num=3) ** 2)
    levels = cred_levs
    print ("levels=",levels)
    #set plot range limits
    if(False):
        rangelimits=[corner_range_cut for x in pars] #precentile cuts
    else:
        rangelimits=[[0,0]]*len(pars) #explicit cuts.
        q = [0.5 - 0.5*corner_range_cut, 0.5 + 0.5*corner_range_cut]
        datum=data[0]
        for i in range(len(pars)):
            rangelimits[i] = list(corner.quantile(datum[:,i], q))
        for datum in data[1:]:
            for i in range(len(pars)):
                xmin, xmax =  corner.quantile(data[0][:,i], q)
                if( xmin < rangelimits[i][0] ):rangelimits[i][0]=xmin
                if( xmax > rangelimits[i][0] ):rangelimits[i][1]=xmax
        #Now enlarge just a bit
        #rangelimits=[[xmin-0.1*(xmax-xmin),xmax+0.1*(xmax-xmin)] for xmin,xmax in rangelimits ]

    figure = corner.corner(data[0], bins=50,labels=names,levels=levels,cov=cov,
                             truths=pars,quantiles=[0.159, 0.5, 0.841],show_titles=True,range=rangelimits,use_math_text=True,
                             title_args={"fontsize": 35},title_fmt='.2e',smooth1d=None,smooth=1,label_kwargs={"fontsize":30},hist_kwargs={"normed":True})
    count=0
    print ("nsamp[0] = ",len(data[0]))
    for datum in data[1:]:
        print ("nsamp[i] ",len(datum))
        figure = corner.corner(datum, bins=50,labels=names,fig=figure,color=["r","b","g","m"][count],levels=levels,cov=cov, plot_density=False,
                               truths=pars,quantiles=[0.159, 0.5, 0.841],show_titles=True,range=rangelimits,use_math_text=True,
                               title_args={"fontsize": 35},title_fmt='.2e',smooth1d=None,smooth=1,label_kwargs={"fontsize":30},hist_kwargs={"normed":True})
        count+=1
    figure.gca().annotate(run+"  SNR="+str(snr)+modes+res+mmodal, xy=(0.5, 1.0), xycoords="figure fraction",
                          xytext=(0, -5), textcoords="offset points",
                          ha="center", va="top",fontsize=30*fontscale)
    for i in range(Npar):
        figure.gca().annotate(names[i]+"= %.3e"%+pars[i], xy=(0.75, 0.9-i/20.0), xycoords="figure fraction",fontsize=30*fontscale)
    figure.savefig(outpath)
    return

#Not yet working with python3
def sky_area_compare(samples,cov,cred_levs,ilon=6,ilat=7,cred_file=None):
    nmax=2500 #samples will be cut down to (near) this number
    ntol=0.1 #this is the fractional tolerance on the number of points
    #several experiments ( with clean gaussians) indicate that there is 
    #little benefit of increasing this number beyond 1500 but considerable cost.
    #results are typically consistent within +/-5% with the calculation taking
    #about 35s,at 3000, the accuracy seems similar but it takes about 2 min.
    
    #convert to 1-D sigma scale factors: lev="percent of dist, within n sigma in each dir"
    #convert by lev  =  1 - exp(-fac^2/2)
    #thus       fac2 = sqrt(-2*log( 1 - lev )) 
    sigma_fac2s = [  -2 * math.log( 1 - x ) for x in cred_levs ]
    print ("sqrt_sigmafac2s:",[math.sqrt(x) for x in sigma_fac2s])
    if(isinstance(cov,list)):
        covar=cov[0]
    else:
        covar=cov
    #Compute sky area from Fisher covariance matrix
    val=covar[ilon][ilon]*covar[ilat][ilat]-covar[ilon][ilat]**2
    if val<0:
        dsky=float('nan')
        if -val<1e-13*covar[6][6]*covar[7][7]:dsky=0
    else: dsky=math.pi*math.sqrt(val)*math.cos(pars[ilat])
    print ("sky",val,dsky,covar[ilon][ilon],covar[ilat][ilat])
    fisher_areas=[dsky * x for x in sigma_fac2s]

    #Prepare to use Will Farr's sky area tool:
    #uses "RA and DEC in radians" 
    #we give it lon and lat
    print (samples.shape)
    pts=zip(samples[:,6],samples[:,7])
    if(len(pts)>nmax*(1+ntol)):
        print ("Trimming points for sky average comparison")
        newpts=[]
        while(len(newpts)<nmax*(1.0-ntol)):
            newpts+=[random.choice(pts) for i in range(int(nmax*(1.0+ntol)-len(newpts)))]
            newpts=list(sets.Set(newpts))
            print (len(newpts))
        pts=newpts
    print ("selected ",len(pts)," pts")
    pts=np.array(pts)
    print ("sky pts:",)
    posterior=skyarea.ClusteredSkyKDEPosterior(pts,ntrials=8,acc=0.001)
    sample_areas=posterior.sky_area(cred_levs,fast=True)
    print ("sky credibility\tfisher_est\tchain_est\tratio")
    for i in range(len(cred_levs)):
        print (cred_levs[i],fisher_areas[i],sample_areas[i],fisher_areas[i]/sample_areas[i])
    return fisher_areas,sample_areas
        

parser = argparse.ArgumentParser(description="Run standard analysis for runs started in Aug 2017");
parser.add_argument('chainfile',help="The name of the chain output file", nargs="*")
parser.add_argument('--fishfile',help="The name of the Fisher output file")
parser.add_argument('-l',help="Set the length of the effective portion of the chain at end default=-1,use last 1/4).",type=int,default=-1)
parser.add_argument('--lens',help="Set individual lengths of the effective portion of ends of the chains.",default="")
parser.add_argument('--acls',help="Individual autocorrelation lengths. May be used in defining plots. ).",default="")
parser.add_argument('--tag',help="Set a tag for the output.",default="")
parser.add_argument('--allfish',help="Use all Fisher runs together.",action="store_true")
    
args=parser.parse_args()

ncases=1

chainfiles=args.chainfile
fishfile=None
if args.fishfile is not None: fishfile=args.fishfile

parnames=[r"$m_1$",r"$m_2$",r"$t_0$",r"$D$",r"$\phi_0$",r"$\iota$",r"$\lambda$",r"$\beta$",r"$\psi$"]

#define the relevant credibility levels:
#cred_levs=[0.5,             #inter-quartile
#           0.68268949213709,#1-sigma
#           0.95449973610364,#2-sigma
#           0.99730020393674,#3-sigma
#           ]
cred_levs=[0.68268949213709,#1-sigma
           0.95449973610364,#2-sigma
           0.99730020393674,#3-sigma
           ]
corner_range_cut=0.9997
print ("chainfiles:",chainfiles)
print ("fishfile:",fishfile)

if(len(args.lens)>0):
    lens=[int(a) for a in args.lens.split()]
    print("lens=",lens)

if(len(args.acls)>0):
    acls=[float(a) for a in args.acls.split()]
    print("acls=",acls)

testing=False

if(True):
    chainfile=chainfiles[0]
    print ("Processing posterior in ",chainfile)
    print (" with Fisher results in ",fishfile)
    #run=re.search("Run\d*",chainfile).group(0)
    run=os.path.basename(chainfile)
    #run=os.path.basename(fishfile)
    if(None==re.search("lm22",chainfile)):modes=""
    else: modes=" (2-2 only)"
    if(None==re.search("nl8k",chainfile)):res=""
    else: res=" (higher-res sampling)"
    if(None==re.search("nmm",chainfile)):mmodal=""
    else: mmodal=" (no modal decomp sampling)"
    print ("run=",run)
    if("post_equal_weights.dat" in chainfile):
        code="bambi"
        burnfrac=0
        injfile=re.sub("post_equal_weights.dat","params.txt",chainfile)
        if(not os.path.exists(injfile)):
            injfile=re.sub("post_equal_weights.dat","_params.txt",chainfile)
    elif("_t0.dat" in chainfile):
        code="ptmcmc"
        burnfrac=0.75
        injfile=re.sub("_t0.dat","params.txt",chainfile)
        if(not os.path.exists(injfile)):
            injfile=re.sub("post_equal_weights.dat","_params.txt",chainfile)
    elif("fishcov.dat" in fishfile):
        code="ptmcmc"
        burnfrac=0
        injfile=re.sub("_all_fishcov.dat","params.txt",fishfile)
        injfile=re.sub("_fishcov.dat","params.txt",injfile)
        if(not os.path.exists(injfile)):
            sys.exit("inj file '"+injfile+"' not found")
    else:   sys.exit("don't know how to find inj file")
    print ("injfile=",injfile)
    dirname=os.path.dirname(os.path.abspath(chainfile))
    basename=os.path.basename(dirname)
    outbasename=basename[:]
    if(len(args.tag)>0):
        ibase=basename.find("bambi")
        if(ibase<0):ibase=basename.find("ptmcmc")
        outbasename=basename[:ibase]+args.tag
    if(testing):
        outpath=dirname+"/"+basename+"_TEST_corner.png"
        credfile=dirname+"/"+basename+"_TEST_credibility_levels.txt"
    else:
        outpath=dirname+"/"+outbasename+"_corner.png"
        outsamples=dirname+"/"+outbasename+"_samples.dat"
        credfile=dirname+"/"+basename+"_credibility_levels.txt"
    print ("reading posterior samples from file:",chainfile)
    print ("reading injection from file:",injfile)
    print ("corner output to ",outpath)
    if(args.allfish):
        cov=[]
        for tag in ['a','b','c','d','k','l']:
            newfish=fishfile.replace("_c_fishcov","_"+tag+"_fishcov")
            if(os.path.isfile(newfish)):
                print ("reading covariance from ",newfish)
                cov.append(readCovar(newfish))
    elif(fishfile is not None):
        print ("reading covariance from ",fishfile)
        cov=readCovar(fishfile)
    else:
        cov=None
    pars,snr=read_injection(injfile)
    if(snr==0):
        snrfile=re.sub("_t0.dat",".out",chainfile)
        print ("snrfile:'"+snrfile+"'")
        snr=read_snr(snrfile)
    print ("SNR=",snr)
    print ("Injected pars=",dict(zip(parnames,pars)))
    samples=[]
    if(testing):
        #cov=np.identity(len(pars))
        #pars=np.zeros(len(pars))
        samples.append(fake_samples(pars,cov,n=60000))
    i=0
    for chainfile in chainfiles:
        if("post_equal_weights.dat" in chainfile):
            code="bambi"
            burnfrac=0
            iskip=1
        elif("_t0.dat" in chainfile):
            code="ptmcmc"
            burnfrac=0.75
            iskip=50
            print("ACLS=",args.acls)
            if(len(args.acls)>0):iskip=int(50+acls[i]/100.0)
        elif("fishcov.dat" in fishfile):
            code="ptmcmc"
            burnfrac=0
            iskip
        lkeep=args.l
        if(len(args.lens)>0):lkeep=lens[i]
        samples.append(read_samples(pars,chainfile,code=code,burnfrac=burnfrac,keeplen=lkeep,iskip=iskip))
        i+=1
    print ("credfile='"+credfile+"'")
    print ("\nsamples:\n",samples)
    if(False):
        fisher_areas,sample_areas=sky_area_compare(samples,cov,cred_levs,ilon=6,ilat=7)
        with open(credfile,"w") as fd:
            fd.write("sky credibility\tfisher_est\tchain_est\tratio\n")
            for i in range(len(cred_levs)):
                fd.write(str(cred_levs[i])+"\t"+str(fisher_areas[i])+"\t"+str(sample_areas[i])+"\t"+str(fisher_areas[i]/sample_areas[i])+"\n")

    cornerFisher(pars,samples,cov,cred_levs)
    rs=[ math.sqrt(s[1]**2+s[2]**2) for s in samples[0]]
    rs.sort()
    nn=len(rs)
    rcuts=[rs[int(x*nn)] for x in cred_levs]
    print (dict(zip(cred_levs,rcuts)))
    
