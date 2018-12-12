"""Python interface to FLARE."""

import numpy as np
from ctypes import cdll, Structure, c_double, c_int, c_void_p, byref


class LISAParams(Structure):
    _fields_ = [('tRef', c_double),
                ('phiRef', c_double),
                ('m1', c_double),
                ('m2', c_double),
                ('distance', c_double),
                ('lambda_', c_double),
                ('beta', c_double),
                ('inclination', c_double),
                ('polarization', c_double),
                ('nbmode', c_int)]

    def __init__(self, p=None):
        if type(p) is dict:
            self.m1 = p['mass1']
            self.m2 = p['mass2']
            self.distance = p['distance']
            self.lambda_ = p['lambda']
            self.beta = p['beta']
            self.inclination = p['inclination']
            self.polarization = p['polarization']
            self.nbmode = p['n_modes']


lib = cdll.LoadLibrary("libflare.so.1.0.0")

# tell ctypes the return types for the functions we'll use
lib.CalculateLogLReIm.restype = c_double
lib.CalculateOverlapReIm.restype = c_double

# must initialize flare's global internal variables first
lib.InitGlobalParams()

def generate_injection_reim(**kwa):
    """Generate a simulated signal in LISA.
    """
    # initialize a LISAInjectionReIm struct
    injection = c_void_p(None)
    lib.LISAInjectionReIm_Init(byref(injection))

    params = LISAParams(kwa)

    # realize an injected signal into the LISAInjectionReIm struct
    lib.LISAGenerateInjectionReIm(byref(params), c_double(kwa['f_low']),
                                  c_int(kwa['n_samples']),
                                  c_int(kwa['log_samples']), injection)

    return injection

def calculate_logl_reim(**kwa):
    """Calculate the log-likelihood given an injection and the parameters of a
    template.
    """
    params = LISAParams(kwa)

    return lib.CalculateLogLReIm(byref(params), kwa['injection'])

def calculate_overlap_reim(params1, params2, injection):
    """Calculate the overlap between signals of different parameter vectors.
    """
    p1 = LISAParams(params1)
    p2 = LISAParams(params2)

    return lib.CalculateOverlapReIm(p1, p2, injection)