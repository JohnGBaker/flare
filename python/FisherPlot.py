import flare
import numpy as np

flare.flare_dir="../flare"
#file="FisherStudy-most-HH.dat"
#file="FisherStudy.dat"
#label="test11.15"
label="L3LISA-v1-sens-but-5Gm"
file=label+"-FisherStudy.dat"
#flare.FisherPlot(label,9,[1.1,2.0,4.0,10.0],[10,100,1000],[0.1,0.3],file)
#flare.FisherPlot(label,0,[1.1,2,4,10],[10,100,1000],[0.1,.3],file)
flare.FisherPlot(label,0,[2],[10],[.3],file,scaled=True)
flare.FisherPlot(label,1,[2],[10],[.3],file,scaled=True)
flare.FisherPlot(label,3,[2],[10],[.3],file,scaled=True)
flare.FisherPlot(label,9,[2],[10],[.3],file,scaled=True,targetSNR=10)
#flare.FisherPlot(label,0,[2],[10],[.3],file,scaled=False)
