import numpy as np


def MagnitudePowerLaw(mag):
    """ Hoekstra et al., Figure 1"""
    return 10**(0.36*mag - 0.36*22)
    
def HoekstraNz(path='./', survey_area=[180.,180.]):
    Mags = np.load(path+'Mags.npy')
    nzPerMag = np.load(path+'NzPerMag.npy')
    nzPerMag = (nzPerMag.T / np.sum(nzPerMag,axis=1)).T
    zs = np.load(path+'zs.npy')
    Ngals = np.array([MagnitudePowerLaw(mag) for mag in Mags])
    nobj_per_mag = Ngals * np.prod(survey_area)
    
    nzPerMag *= nobj_per_mag.reshape(-1,1)
    NumberCounts = np.sum(nzPerMag,axis=0)
    
    return NumberCounts, zs
    

