import flare
import numpy as np
import argparse
eps=0.3

parser = argparse.ArgumentParser(description="generate a discover queue script for a flare run");
parser.add_argument('label',help="The basename/directory for the run")
args=parser.parse_args()

flare.flare_dir="../flare"
#file="FisherStudy-most-HH.dat"
#file="FisherStudy.dat"
label=args.label+"/"
#label="test11.15"
#label="L3LISA-v1-sens-but-5Gm-test-wide"
#label="L3LISARef"
#label="LISA2017camp_10yr/"
#label="2arm-LISA/"
#label="tRef-redef-LISA2017/"
#label="LISA2017-Nov-flaretest/"
#label="slow-orbit-LISA/"
#label="fast-orbit-LISA/"
#label="big-orbit-LISA/"
#label="slower-orbit-LISA/"
#label="tiny-orbit-LISA/"
file=label+"FisherStudy.dat"
#flare.FisherPlot(label,9,[1.1,2.0,4.0,10.0],[10,100,1000],[0.1,0.3],file)
#flare.FisherPlot(label,0,[1.1,2,4,10],[10,100,1000],[0.1,.3],file)
#flare.FisherPlot(label,0,[2],[10],[.3],file)
flare.HorizonPlot(label,0,[2],10,eps,file,[0.001,0.003,0.01,0.03,0.10,0.30],scaled=True)
flare.HorizonPlot(label,0,[2],10,eps,file,[0.001,0.01,0.10],scaled=True,show_range=True)
flare.HorizonPlot(label,1,[2],10,eps,file,[0.001,0.003,0.01,0.03,0.10,0.30],scaled=True)
flare.HorizonPlot(label,1,[2],10,eps,file,[0.001,0.01,0.10],scaled=True,show_range=True)
#flare.HorizonPlot(label,2,[2],10,[.3],file,[0.10,0.30,1.0,3.0,10.0,30.0,100.0],scaled=True)
flare.HorizonPlot(label,3,[2],10,eps,file,[0.003,0.01,0.03,0.10,0.30],scaled=True)
flare.HorizonPlot(label,3,[2],10,eps,file,[0.01,0.10,0.5],scaled=True,show_range=True)
#flare.HorizonPlot(label,9,[2],10,[.3],file,[8.4e-7,8.4e-6,8.4e-5,3.0e-4,3.0e-3,3.0e-2],scaled=True)
#flare.HorizonPlot(label,9,[2],10,[.3],file,[100,900,3600,28800,360000],scaled=True)
#flare.HorizonPlot(label,9,[2],10,[.3],file,[100,3600,90000],scaled=True,show_range=True)
flare.HorizonPlot(label,9,[2],10,eps,file,[0.0278,0.25,1.0,9.0,100.0],scaled=True)
flare.HorizonPlot(label,9,[2],10,eps,file,[0.0278,1.0,25.0],scaled=True,show_range=True)
#flare.FisherPlot(label,0,[2],10,[.3],file,scaled=False)
