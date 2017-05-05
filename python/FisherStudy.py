import flare
import numpy as np

#flare.flare_dir="../../LISAstudy/flaretest/flare/"
#label="LISA2017-Nov-flaretest/"
flare.flare_dir="../flare"
label="LISA2017b/"
#label="2arm-LISA/";flare.extra_flags=" --tagtdi TDIAXYZ"
#label="tRef-redef-LISA2017/"
#label="fast-orbit-LISA/"
#label="slow-orbit-LISA/"
#label="slower-orbit-LISA/"
#label="tiny-orbit-LISA/"
#Ms=1.5e4*10**(np.arange(16)/3.0)
Ms=(1.5e6*10**(np.arange(10)/3.0)).tolist()
#Ms=[1.5e3*3,1.5e4,1.5e5,1.5e5*10**0.5]+Ms
Ms=[150,470,1.5e3,4.7e3,1.5e4,1.5e5,1.5e5*10**0.5]+Ms
#Ms=[4.7e3,1.5e4,1.5e5,1.5e5*10**0.5]+Ms
#Ms=[1.5e3,4.7e3,1.5e4,1.5e5,1.5e5*10**0.5]+Ms
print "Ms=",Ms
flare.SNRstudy(label,Ms,[2.0],[10,100,1000],300)
#/////#/////flare.SNRstudy(label,Ms,[2.0],[10,100,1000],30,8)
#flare.SNRstudy(label,Ms,[2.0],[100],30,8)
##flare.SNRstudy(Ms,[1,2,4,10],[10,100,1000],300)
#flare.FisherStudy(Ms,[1],[10],[0.1,0.3],6,6)
#flare.noRun=True
#####flare.all_params_file = open(label+"FisherStudy-allparams.dat","w")
#flare.FisherStudy(Ms,[1.1,2.0,4.0,10.0],[10,100,1000],[0.1,0.3],120,6)
#Ms=(1.5e4*10**(np.arange(80)/15.0)).tolist()
#Ms=(3.162e4*10**(np.arange(20)/4.0)).tolist()
#flare.FisherStudy(label,Ms,[2.0],[10],[0.3],100,8)
#flare.FisherStudy(label,Ms,[2.0],[10],[0.1],10,8,4.0)
flare.FisherStudy(label,Ms,[2.0],[10],[0.3],30,8)
#flare.FisherStudy(label,Ms,[2.0],[100],[0.3],5,8)
#flare.all_params_file.close()
