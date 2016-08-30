import flare
import numpy as np

flare.flare_dir="../flare"
#Ms=1.5e4*10**(np.arange(16)/3.0)
Ms=(1.5e6*10**(np.arange(10)/3.0)).tolist()
#Ms=[1.5e4,1.5e5,1.5e5*10**0.5]+Ms
print "Ms=",Ms
#SNRstudy(Ms,[1,2,4,10],[10,100,1000],300)
#flare.FisherStudy([1.5e4],[1],[10,100,1000],[0.1,0.3],30,6)
#flare.FisherStudy([1.5e6,1.5e8],[1],[10],[0.1,0.3],12,6)
#flare.FisherStudy(Ms,[1],[10],[0.1,0.3],6,6)
flare.FisherStudy(Ms,[1.1,2.0,4.0,10.0],[10,100,1000],[0.1,0.3],60,6)
