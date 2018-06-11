""" Dumb script to compute n(z) parameters per magnitude and make a bunch of plots to 
    justify(ish) our extrapolations."""

import matplotlib.pyplot as plt
import numpy as np

def PerMagNz(z, a, b, c, A):
    """ n(z) as in equation (8) of Ilbert et al., 2008
    """
    num = z**a + z**(a*b)
    den = z**b + c
    return A*num/den

# Values from Table 2, Ilbert et al., 2008
Mags = np.array([22, 22.5, 23, 23.5, 24, 24.5])+.25
ExtraMags = np.arange(20, 25, .5) + .25 
aPerMag = [0.497, 0.448, 0.372, 0.273, 0.201, 0.126]
bPerMag = [12.643, 9.251, 6.736, 5.281, 4.494, 4.146]
cPerMag = [0.381, 0.742, 1.392, 2.614, 3.932, 5.925]
APerMag = [4068.19, 9151.98, 18232.24, 35508.58, 60306.30, 103340.04]

# fit polynomials - a and c
apol = np.poly1d(np.polyfit(Mags, aPerMag, 1))
extraa = [apol(mag) for mag in ExtraMags]

cpol = np.poly1d(np.polyfit(Mags,np.log(cPerMag),1))
extrac = np.exp([cpol(mag) for mag in ExtraMags])
#cpol = np.poly1d(np.polyfit(Mags, cPerMag, 2))
#extrac = [cpol(mag) for mag in ExtraMags]

# plot em
plt.plot(Mags, aPerMag, label='a', c='midnightblue')
plt.plot(ExtraMags, extraa, '--', c='midnightblue')
plt.plot(Mags, cPerMag, label='c', c='firebrick')
plt.plot(ExtraMags, extrac, '--', c='firebrick')
plt.legend()
plt.savefig('Plots/a_c.png')
plt.close()

# fit b
bpol = np.poly1d(np.polyfit(Mags, bPerMag, 2))
extrab = [bpol(mag) for mag in ExtraMags]

# plot it
plt.plot(Mags, bPerMag, label='b', c='darkorchid')
plt.plot(ExtraMags, extrab, '--', c='darkorchid')
plt.legend()
plt.savefig('Plots/b.png')
plt.close()

# and fit A
Apol = np.poly1d(np.polyfit(Mags,np.log(APerMag),1))
extraA = np.exp([Apol(mag) for mag in ExtraMags])
#Apol = np.poly1d(np.polyfit(Mags, APerMag, 2))
#extraA = [Apol(mag) for mag in ExtraMags]

# plot it
plt.plot(Mags, APerMag, label='A', c='cornflowerblue')
plt.plot(ExtraMags, extraA, '--', c='cornflowerblue')
plt.legend()
plt.savefig('Plots/A.png')
plt.close()

# zs
zmin, zmax = 0, 2.
dz = 0.06667
bin_edges = np.arange(zmin,zmax+dz,dz)

zmid = [rf-(rf-lf)/2 for lf, rf in zip(bin_edges, bin_edges[1:])]

ExtranzPerMag = np.array([[PerMagNz(z, a, b, c, A) for z in zmid] for (a,b,c,A) in zip(extraa,extrab,extrac,extraA)])
for nz,m in zip(ExtranzPerMag, ExtraMags):
    plt.plot(zmid, nz, label='m={}'.format(m))
plt.legend()
plt.xlabel(r'$z$')
plt.ylabel(r'$n(z)$')
plt.title(r'$n(z)$ per Magnitude bin (Extrapolated parameters)')
plt.savefig('Plots/nz_paramsextrapolated.png')
plt.close()


nzPerMag = np.array([[PerMagNz(z, a, b, c, A) for z in zmid] for (a,b,c,A) in zip(aPerMag,bPerMag,cPerMag,APerMag)])
for nz,m in zip(nzPerMag, Mags):
    plt.plot(zmid, nz, label='m={}'.format(m))
plt.legend(loc=3)
plt.xlabel(r'$z$')
plt.ylabel(r'$n(z)$')
plt.title(r'$n(z)$ per Magnitude bin (Ilbert et al., 2008)')
plt.savefig('Plots/nz_COSMOS.png')
plt.close()

# ok this is not working. Let's try something else. Flip it around and plot n(z_k) as function of m
plt.rcParams['figure.figsize'] = plt.rcParamsDefault['figure.figsize']
ExtraMags = np.arange(22, 26, .5) + .25
z_splitidx = 9 # index to separate "low" and "mid" redshifts mostly for plot readability
z_splitidxhi = 30 # index to separate "mid" and "high" redshifts mostly for plot readability
extranz = []

for jz,z in enumerate(zmid[:z_splitidx]): 
    #zpol = np.poly1d(np.polyfit(Mags, nzPerMag[:,jz], 2))
    zpol = np.poly1d(np.polyfit(Mags[-2:], nzPerMag[-2:,jz], 2))
    extranz += [[max(0.,zpol(mag)) for mag in ExtraMags]]
    plt.plot(Mags, nzPerMag[:,jz], label=r'$z={}$'.format(round(z,3)), c='midnightblue', alpha=(5.+jz)/20.,
             linewidth=2-(jz/17.))
    plt.plot(ExtraMags, extranz[-1], '--', c='midnightblue', alpha=(5.+jz)/20.,
             linewidth=2-(jz/17.))
plt.legend(loc=5, bbox_to_anchor=(1.1,.5))
plt.xlabel(r'Magnitude $m$')
plt.ylabel(r'$n(z)$')
plt.savefig('Plots/lowz_highm.png')
plt.close()

for jz,z in enumerate(zmid[z_splitidx:z_splitidxhi]): 
    zpol = np.poly1d(np.polyfit(Mags[-2:], nzPerMag[-2:,jz+z_splitidx], 1))
    extranz += [[max(0.,zpol(mag)) for mag in ExtraMags]]
    plt.plot(Mags, nzPerMag[:,jz+z_splitidx], label=r'$z={}$'.format(round(z,3)), c='midnightblue', alpha=(5.+jz)/20.,
             linewidth=2-(jz/17.))
    plt.plot(ExtraMags, extranz[-1], '--', c='midnightblue', alpha=(5.+jz)/20.,
             linewidth=2-(jz/17.))
plt.legend(loc=5, bbox_to_anchor=(1.1,.5))
plt.xlabel(r'Magnitude $m$')
plt.ylabel(r'$n(z)$')
plt.savefig('Plots/midz_highm.png')
plt.close()


#ExtraMags = np.arange(23.5, 26, .5) + .25

'''for jz,z in enumerate(zmid[z_splitidxhi:]): 
    zpol = np.poly1d(np.polyfit(Mags[-2:], nzPerMag[-2:,jz+z_splitidxhi], 1))
    extranz += [[max(0.,zpol(mag)) for mag in ExtraMags]]
    plt.plot(Mags, nzPerMag[:,jz+z_splitidxhi], label=r'$z={}$'.format(round(z,3)), c='midnightblue', alpha=(5.+jz)/20.,
             linewidth=2-(jz/17.))
    plt.plot(ExtraMags, extranz[-1], '--', c='midnightblue', alpha=(5.+jz)/20.,
             linewidth=2-(jz/17.))
plt.legend(loc=5, bbox_to_anchor=(1.1,.5))
plt.xlabel(r'Magnitude $m$')
plt.ylabel(r'$n(z)$')
plt.savefig('Plots/highz_highm.png')
plt.close()'''

extranz = np.array(extranz)


for nz,m in zip(extranz.T, ExtraMags):
    plt.plot(zmid, nz, label='m={}'.format(m), c='darkorchid')
    
for nz,m in zip(nzPerMag, Mags):
    plt.plot(zmid, nz, label='m={}'.format(m), c='red')
#plt.legend()
plt.xlabel(r'$z$')
plt.ylabel(r'$n(z)$')
plt.title(r'$n(z)$ per Magnitude bin (Extrapolated n(z))s')
plt.savefig('Plots/nz_nzsextrapolated.png')
plt.close()


# merge the two
AllMags = np.arange(20, 26, .5) + .25
nzFinal = np.empty((len(AllMags),len(zmid)))

# low m's: use parameter extrapolation
#nzFinal[:4,:] = ExtranzPerMag[:4,:]
nzFinal[:10,:] = ExtranzPerMag

# mid m's: average over two interpolations
#nzFinal[4:10,:] = (ExtranzPerMag[4:,:] + extranz.T[:6])/2

# high m's: use n(z) extrapolation
nzFinal[10:] = extranz.T[6:]

# PLOOOOT
plt.rcParams['figure.figsize'] = 9,6


alphas = np.linspace(1.,.1,nzFinal.shape[0])
ls = np.linspace(1,4,nzFinal.shape[0])
for j,(nz,m) in enumerate(zip(nzFinal, AllMags)):
    alpha = alphas[j]#.1+float(j)/16
    l = ls[j] #2.8-(float(j)/10)
    plt.plot(zmid, nz, label='m={}'.format(m), c='crimson', alpha=alpha, linewidth=l)
plt.legend()
plt.xlabel(r'$z$')
plt.ylabel(r'$n(z)$')
plt.legend(loc=5, bbox_to_anchor=(1.25,.5))
plt.title(r'$n(z)$ per Magnitude bin (Extrapolated)')
plt.savefig('Plots/nz_allextrapolated.png', bbox_inches='tight')
plt.close()

# save for importation within cutoff_nz.py
np.save('./Mags.npy', AllMags)
np.save('./nzPerMag.npy', nzFinal)
