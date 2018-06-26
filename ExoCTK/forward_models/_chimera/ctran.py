from . import _tran_module
import numpy as np
# import pdb

def CalculateTau(Xsects, RXsects, Z, P, T, Fractions, r0):
    '''
    Returns array of optical lengths at each level and wavenumber.
    '''
    Tau = np.zeros((Xsects.shape[1], Xsects.shape[0] + 1))
    #pdb.set_trace()
    _tran_module._tau_wrap(Xsects, RXsects, Z, P, T, Fractions, Tau, r0)
    return Tau

def InitXsects(xsecarr, Tgrid, Pgrid, Tavg, Pavg, wno, wno_offset, a, power):
    count = len(wno)
    numlevels = len(Tavg)+1

    numgases = xsecarr.shape[0]
    xsects_inter = np.zeros((numlevels-1, count, numgases))
    xsects_interRayleigh = np.zeros((numlevels-1, count))
    
    _tran_module._init_xsects_wrap(xsects_inter, xsects_interRayleigh, xsecarr,
                                 Tgrid, Pgrid, Tavg, Pavg, wno, wno_offset, a, power)
    return xsects_inter, xsects_interRayleigh
