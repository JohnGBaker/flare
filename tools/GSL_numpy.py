import Cy_gsl
import numpy as np

### let try it


def ArrayToVec(p):
    #print (p.data)
    N = len(p)
    v = Cy_gsl.gsl_vector_alloc(N)
    v.data = p.data

    return(v)
#
#
#
#
p = np.arange(10)*1.5
#
#
v = ArrayToVec(p)
