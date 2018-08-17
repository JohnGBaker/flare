import os
import subprocess
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

flare_dir=os.path.dirname(os.path.realpath(__file__))[:-7]

def callComputeSNR(f):
    times=[0]*5
    snr=0
    cmd = flare_dir+"/LISAinference/ComputeLISASNR --loadparamsfile --paramsfile q3z4_22_inc_03_id_9.par --outputfile snr.out --paramsdir . --nlinesparams 1 --nbmode 5 --verbosetmax --maxf22 "+str(f)
    p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    output = p.stdout.read()
    #print(output.decode("utf-8"))
    lines=output.decode("utf-8").split('\n')
    for line in lines:
        #print("line=",line)
        words=line.split()
        #print(words)
        #for word in words:print(word)
        if(len(words)>0):
            if(words[0]=='m'):times[int(words[3])-1]=float(words[6])
            if(words[0]=='SNR'):snr=float(words[1])
    return snr,times

fs=[]
snrs=[]
timess=[]
for logf in np.arange(-4.1,-1,0.02):
    f=10**logf
    fs.append(f)
    snr,times=callComputeSNR(f)
    snrs.append(snr)
    timess.append(times)
times=np.transpose(timess)+300


snrcuts=np.array([1/64.,1/16.,1/4.])*snrs[-1]
fcuts=np.interp(snrcuts,snrs,fs)
tcuts=np.interp(fcuts,fs,times[1])
    
print(np.transpose([snrcuts,fcuts,tcuts]))

plt.semilogx(fs,snrs,label="SNR")
plt.semilogx(fcuts,snrcuts,'*k')
plt.xlabel("f (Hz)");plt.ylabel("SNR")
plt.legend()
plt.show()
plt.loglog(fs,times[0],label="l=2,m=1")
plt.loglog(fs,times[1],label="l=2,m=2")
plt.loglog(fs,times[2],label="l=3,m=3")
plt.loglog(fs,times[3],label="l=4,m=4")
plt.loglog(fs,times[4],label="l=5,m=5")
plt.loglog(fcuts,tcuts,"*k")
plt.xlabel("f (Hz)");plt.ylabel("300-t(f) (s)")
plt.legend()
plt.show()
plt.semilogy(snrs,times[0],label="l=2,m=1")
plt.semilogy(snrs,times[1],label="l=2,m=2")
plt.semilogy(snrs,times[2],label="l=3,m=3")
plt.semilogy(snrs,times[3],label="l=4,m=4")
plt.semilogy(snrs,times[4],label="l=5,m=5")
plt.semilogy(snrcuts,tcuts,"*k")
plt.xlabel("snr");plt.ylabel("300-t(f) (s)")
plt.legend()
plt.show()

