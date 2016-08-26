import math
import numpy as np
import subprocess
import re
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import astropy.units as units
from astropy.cosmology import Planck15 as cosmo,z_at_value
from matplotlib.backends.backend_pdf import PdfPages

def set_flare_flags(snr,params):
    flags=""
    #Waveform model
    #flags += "--tagextpn 0" #Don't extend waveforms at low freq to allow lower masses
    #flags+="--tagint 0 " #1 for Fresnel integration(default), vs gridded quadrature"
    #flags+="--tagint 1 --nbptsoverlap 8192" #gridded quadrature
    flags+="--deltatobs 5.0 " #duration in years of LISA observation
    #flags+="--minf 1e-4 " #minimun frequency included in analysis
    flags+="--minf 3e-6 " #minimun frequency included in analysis
    #flags+="--nbmodeinj 1 --nbmodetemp 1 " #for no higher modes in injection and template
    if(snr>0):
        flags+="--snr "+str(snr)+" --rescale-distprior " #fixing SNR (rescales distance)
    flags+="--comp-min 1e5 --comp-max 1e8 " #min/max for component mass prior ranges
    flags+="--logflat-massprior " #assume prior uniform in log of masses, rather than uniform for mass."
    #flags+="--mtot-min 8e5 --mtot-max 2e8 --q-max 11.98 "  #additional prior limits on Mtot and q
    flags+="--mtot-min 1e4 --mtot-max 1e9 --q-max 11.98 "  #additional prior limits on Mtot and q
    #flags+="-dist-min 5000. --dist-max 200e4 --distance 1e5"  #prior range for distances should verify range based on distances (in Mpc).
    flags+="--dist-min 1000. --dist-max 4e5 "  #prior range for distances approx 0.2<z<33
    #set parameter flags
    m1           = params[0]
    m2           = params[1]
    tRef         = params[2]
    phiRef       = params[3]
    if(snr>0):
        dist     = 2e4
    else:
        dist     = params[4]
    lam          = params[5]
    beta         = params[6]
    inc          = params[7]
    pol          = params[8]
    flags += "--phiRef "+str(phiRef)+" "
    flags += "--m1 "+str(m1)+" "
    flags += "--m2 "+str(m2)+" "
    flags += "--tRef "+str(tRef)+" "
    flags += "--distance "+str(dist)+" "
    flags += "--lambda "+str(lam)+" "
    flags += "--beta "+str(beta)+" "
    flags += "--inclination "+str(inc)+" "
    flags += "--polarization "+str(pol)+" "
    return flags

def set_mcmc_flags(outroot,ptN):
    flags  = ""
    #MCMC basics
    flags += "--rng_seed="+str(np.random.rand())+" " 
    flags += " --outroot "+str(outroot)+" "
    flags += "--nskip=40 --info_every=10000 " #frequency of sampling/reporting
    flags += "--prop=7 --de_ni=500 --gauss_1d_frac=0.5 --de_reduce_gamma=4 " #differential evolution proposal distribution with Gaussian draws 1/2 of the time
    #Parallel Tempering setup
    flags += "--pt --pt_Tmax=1e9 "    #parallel tempering basics
    if(ptN>0):
        flags += "pt_n="+str(ptN)+" " #else default is 20
    flags += "--pt_swap_rate=0.10 "   #rate of temp swaps (or default 0.01)
    flags += "--pt_evolve_rate=0.01 " #rate at which temps are allowed to evolve 

    flags += "--pt_reboot_rate=0.0001 --pt_reboot_every=10000 --pt_reboot_grace=50000 " #Somewhat hacky trick to avoid chains getting stuck.  Not sure whether we need this.
    #stopping criteria
    flags += "--nsteps=1e7" #10 million steps may be about the most we can do
    flags += "--pt_stop_evid_err=0.05" #may terminate earlier based on evidence criterion
    return flags

def set_bambi_flags(outroot):
    flags  = "--nlive 1000 --tol 1.0 --mmodal --nclspar 2 --maxcls 10 --ztol -60 --seed "
    flags += "--outroot "+outroot+" "


def draw_params(Mtot,q):
    #we suppose fixed Mtot,q,SNR and draw the other params
    m1     = Mtot*q/(1.0+q)
    m2     = Mtot/(1.0+q)
    tRef   = np.random.randn()*1e5
    phiRef = 2*math.pi*np.random.rand()
    dist = 100*10**(np.random.rand()*math.log10(400))
    lam    = np.random.rand()*2.0*math.pi
    beta   = math.acos(np.random.rand()*2.0-1)-math.pi/2.0
    inc    = math.acos(np.random.rand()*2.0-1)
    pol    = np.random.rand()*math.pi
    params = [m1,m2,tRef,phiRef,dist,lam,beta,inc,pol]
    return params

def perform_run(name,Mtot,q,snr):
    if(BAMBI):
        cmd   = flare_dir+"/LISAinference/LISAinference"
        flags = get_bambi_flags(name)
    else:
        cmd   = flare_dir+"/LISAinference/LISAinference_ptmcmc"
        flags = get_mcmc_flags(name,60)
    params=draw_params(Mtot,q)
    flags+=set_flare_flags(snr,params)
    subprocess.call(cmd+" "+flags)

    
def SNRrun(Mtot,q,snr):
    cmd   = flare_dir+"/LISAinference/LISAinference_ptmcmc"
    flags = "--nsteps=0 --noFisher "
    params=draw_params(Mtot,q)
    flags+=set_flare_flags(snr,params)
    name="dummy"
    flags += "--rng_seed="+str(np.random.rand())+" " 
    flags += " --outroot "+str(name)+" "
    cmd += " "+flags+">"+name+".out"
    setenv = "export ROM_DATA_PATH=/Users/jgbaker/Projects/GWDA/LISA-type-response/flare/ROMdata/q1-12_Mfmin_0.0003940393857519091"
    
    print "Executing '"+cmd+"'"
    code=subprocess.call(setenv+";"+cmd,shell=True)
    print "Run completed with code(",code,")"
    with open(name+"params.txt",'r') as file:
        lines=file.read()
        #print lines
        dist=re.search("dist_resc:(.*)", lines).group(1)
        print "distance =",dist
    return float(dist)

def SNRstudy(MtotList,qList,SNRList,Navg):
    pp = PdfPages('SNRstudy.pdf')
    for q in qList:
        tags=[]
        labels=[]
        count=0
        for snr in SNRList:
            count+=1
            y1=[]
            y2=[]
            x=[]
            for Mtot in MtotList:
                print "Running SNRrun(",Mtot,",",q,",",snr,")"
                dists=np.zeros(Navg);
                for i in range(Navg):
                    dist=SNRrun(Mtot,q,snr)
                    z=z_at_value(cosmo.luminosity_distance,dist*units.Mpc,zmax=10000,ztol=1e-6)
                    print "D=",dist," z=",z
                    dists[i]=math.log10(z)
                    #dists[i]=math.log10(dist)
                mean=np.mean(dists);
                std=np.std(dists);
                print "M=",Mtot," q=",q,"dist=",mean,"+/-",std
                x.append(math.log10(Mtot/(1+10**mean)))
                #x.append(math.log10(Mtot))
                y1.append(mean-std)
                y2.append(mean+std)
            print "x=",x
            print "y1=",y1
            print "y2=",y2
            color=(0.2,0.8/math.sqrt(q),1.0/math.sqrt(count))
            plot=plt.fill_between(x, y1, y2, facecolor=color,alpha=0.3, interpolate=True)
            tags.append( Rectangle((0, 0), 1, 1, fc=color,alpha=0.3) )
            labels.append("SNR="+str(snr))
        plt.legend(tags,labels)
        plt.ylim([-1,3])
        plt.xlim([2,9])
        plt.title("SNR contours for LISA q="+str(q)+" SMBH merger")
        plt.ylabel("log(z)")
        plt.xlabel("log(M/Msun)")
        #plt.show()
        pp.savefig()
        plt.clf()
    pp.close()
    
flare_dir="../flare"
Ms=1.5e4*10**(np.arange(16)/3.0)
#Ms=2.0e5*10**(np.arange(13)/3.0)
print "Ms=",Ms
SNRstudy(Ms,[1,2,4,10],[10,100,1000],300)
#logz = np.arange(10)/2.5
#print "logz=",logz
#print [10**x for x in logz]
#logD = [cosmo.luminosity_distance(1+10**lz)/units.Mpc for lz in logz]
#print logD
#plt.clf()
#plot=plt.plot(logz,logD)
#plt.show()
