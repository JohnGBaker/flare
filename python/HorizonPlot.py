import flare
import numpy as np

flare.flare_dir="../flare"
#file="FisherStudy-most-HH.dat"
#label="test11.15"
label="L3LISA-v1-sens-but-5Gm"
file=label+"-FisherStudy.dat"
#flare.FisherPlot(label,9,[1.1,2.0,4.0,10.0],[10,100,1000],[0.1,0.3],file)
#flare.FisherPlot(label,0,[1.1,2,4,10],[10,100,1000],[0.1,.3],file)
#flare.FisherPlot(label,0,[2],[10],[.3],file)
flare.HorizonPlot(label,0,[2],[10],[.3],file,[0.001,0.003,0.01,0.03,0.10,0.30],scaled=True)
flare.HorizonPlot(label,1,[2],[10],[.3],file,[0.001,0.003,0.01,0.03,0.10,0.30],scaled=True)
flare.HorizonPlot(label,3,[2],[10],[.3],file,[0.001,0.003,0.01,0.03,0.10,0.30],scaled=True)
flare.HorizonPlot(label,9,[2],[10],[.3],file,[8.4e-7,8.4e-6,8.4e-5,3.0e-4,3.0e-3,3.0e-2],scaled=True)
#flare.FisherPlot(label,0,[2],[10],[.3],file,scaled=False)
