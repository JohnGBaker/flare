import flare
import numpy as np
import sys
import os

#Usage: python FisherRun.py label snr delta param1 param2 ...

flare.flare_dir= os.path.dirname(os.path.realpath(__file__))[:-7]

label=sys.argv[1]
snr=sys.argv[2]
delta=sys.argv[3]
params=sys.argv[4:]
print "label=",label
print "snr=",snr
print "delta=",delta
print "params=",params
datum=flare.FisherRunByParams(snr,params,float(delta),label)
with open(label+"-run.dat",'w') as f:
    for val in datum:
        f.write(str(val)+"\t");
    f.write("\n")
    f.close()
