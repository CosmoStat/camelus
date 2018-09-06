
from __future__ import print_function
import scipy.integrate as integrate
from astropy.io import ascii
import sys
import math
import numpy as np
from scipy import argwhere
from scipy import interpolate
import pylab as plt
import matplotlib.cm as cm
import scipy.special as special
from scipy.interpolate import interp1d
plt.rc('text',usetex=True)
plt.rc('font', family='serif', size=15, serif='cz00')
#plt.close('all')
from astropy import cosmology

from scipy.optimize import curve_fit

fich=' ../output/test_55_z1_pk_nl.dat'

def read_pk(fich):
	data=np.loadtxt(fich)
	k=data[:,0]
	pk=data[:,1]
	return k,pk



def plot_pk(ell, p_kappa):
    plt.figure(1)
    plt.loglog(ell,p_kappa)    
    plt.xlabel('$\ell$')
    plt.ylabel('$ P_\kappa(\ell)$')
    plt.show()


