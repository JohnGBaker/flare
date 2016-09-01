import flare
import numpy as np

flare.flare_dir="../flare"
#file="FisherStudy-most-HH.dat"
file="FisherStudy.dat"
label="test"
#flare.FisherPlot(label,9,[1.1,2.0,4.0,10.0],[10,100,1000],[0.1,0.3],file)
flare.FisherPlot(label,9,[1.1,2,4],[10,100,1000],[0.1,.3],file)
