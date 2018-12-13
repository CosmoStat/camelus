### git clone 
 
from __future__ import print_function
import scipy.integrate as integrate
from astropy.io import ascii
import sys
import math
import numpy as np
from scipy import interpolate
import pylab as plt
import matplotlib.cm as cm
from matplotlib.colors import LogNorm

import scipy.integrate as integrate
import scipy.special as special
from scipy.interpolate import interp1d
import matplotlib.image as mpimg
import matplotlib.colors as colors

from matplotlib.mlab import griddata

plt.rc('text',usetex=True)
plt.rc('font', family='serif', size=12, serif='cz00')

def nz_multi2(catGal,nhalo):
	plt.close('all')
	plt.figure(1)
	nnz2=49


	for iii in range(nhalo):
		print('test {0}'.format(iii))
		catHalo2="{0}{1:03d}".format(catGal,iii+1)
		dat = np.loadtxt(catHalo2)
		z =dat[:,2]
		if iii == 0:
			zmin=min(z)
			zmax=max(z)
			zz2=np.linspace(zmin,zmax,nnz2)
			nz2=np.linspace(zmin,zmax,nnz2)*0
			dz=zz2[1]-zz2[0]

		print("z : {0} {1} ".format(min(z),max(z)))
		
		print("tot gal : {0} ".format(len(z)))
		for i in range(len(z)):
			ii=int(((z[i]-zmin)/dz))
			nz2[ii]=nz2[ii]+1

	nz2=nz2/(sum(nz2))
	plt.plot(zz2,nz2,'+')
	plt.plot(zz2,nz2,color='b')

	#plt.hist(z,bins=10,normed=True,color='b',alpha=0.5)
	zz=np.linspace(zmin,zmax,nnz2)
	nn2=nnz(zz)
	nn2=nn2/sum(nz2)/sum(nn2)
	print("tot gal : {0} ".format(sum(nn2)))
	plt.plot(zz,nn2,color='r')
	plt.show()
	return

def nz_multi(catHalo,nhalo):
	plt.close('all')
	plt.figure(1)
	nnz2=10
	zmin=0.0001 #min(z)
	zmax=2.
	zz2=np.linspace(zmin,zmax,nnz2)
	nz2=np.linspace(zmin,zmax,nnz2)*0
	dz=zz2[1]-zz2[0]

	for iii in range(nhalo):
		print('test {0}'.format(iii))
		catHalo2="{0}{1:03d}".format(catHalo,iii+1)
		dat = ascii.read(catHalo2)
		z =dat['col4']
		print("z : {0} {1} ".format(min(z),max(z)))
		Ngal_c =dat['col6']
		Ngal_s =dat['col7']
		Ntot=Ngal_c*(1.0+Ngal_s) #=ngc*(1+ngs)
		print("tot gal : {0} ".format(sum(Ntot)))
		for j in range(len(z)):
			ii=int(((z[j]-zmin)/dz)-0.5)
			nz2[ii]=nz2[ii]+Ntot[j]		
		nz2=nz2/(sum(nz2))
		plt.plot(zz2,nz2,'+')
		plt.plot(zz2,nz2,color='b')

	#plt.hist(z,bins=10,normed=True,color='b',alpha=0.5)
	zz=np.linspace(0.01,2+0.1,100)
	nn2=nnz(zz)
	nn2=nn2/(sum(nn2))
	print("tot gal : {0} ".format(sum(nn2)))
	plt.plot(zz,nn2,color='r')
	plt.show()
	return

def nz(catHalo):
	plt.close('all')
	dat = ascii.read(catHalo)
	z =dat['col4']
	print("z : {0} {1} ".format(min(z),max(z)))
	Ngal_c =dat['col6']
	Ngal_s =dat['col7']
	Ntot=Ngal_c*(1.0+Ngal_s) #=ngc*(1+ngs)
	print("tot gal : {0} ".format(sum(Ntot)))
	nnz2=30
	print("make HOD")
	zmin=min(z)
	zz2=np.linspace(zmin,max(z),nnz2)
	nz2=np.linspace(zmin,max(z),nnz2)*0
	dz=zz2[1]-zz2[0]
	for j in range(len(z)):
		ii=int(((z[j]-zmin)/dz)-0.5)
		nz2[ii]=nz2[ii]+Ntot[j]		

	plt.figure(1)
	nz2=nz2 #/(sum(nz2))
	plt.hist(z,bins=10,normed=True,color='b',alpha=0.5)
	#zz=np.linspace(0.01,max(z),len(zz2))
	nz1=nnz(zz2)
	nz1=nz1*12.*180.*180./(sum(nz1))
	#for ii in range(1,8):
	#	print("make n(z) {0}".format(ii))
	#	hh= '../build/g_random_00'+str(ii)
	#	dat = ascii.read('../build/g1')
	#	z1 =dat['col3']
	#	nz1=np.linspace(zmin,max(z),nnz2)*0
	#	for j in range(len(z1)):
	#		ii=int(((z1[j]-zmin)/dz)-0.5)
	#		nz1[ii]=nz1[ii]+1	
	#	plt.semilogy(zz2,nz1,color='r',alpha=0.4)
	plt.semilogy(zz2,nz1,color='r',alpha=0.4)
	plt.semilogy(zz2,nz2,color='b')
	plt.show()
	return	
	
def nzg(catHalo):
	plt.close('all')
	dat = ascii.read(catHalo)
	z =dat['col3']
	plt.hist(z,bins=10,normed=True,color='b',alpha=0.5)
	zz=np.linspace(0.0001,max(z),10)
	nn2=nnz(zz)
	#nn2=nn2/(sum(nn2))
	plt.plot(zz,nn2,color='r')*10.*6.
	#plt.plot(zz2,nz2,'+')
	#plt.plot(zz2,nz2,color='b')
	plt.show()
	return	

def nnz(z):
	z0=0.5
	alpha=2.
	beta=1.
	x=z/z0
	return  x**alpha * np.exp(-x**beta)

##############
def mass_function_halo(z):
		fich='../build/massFct_z{0:.3f}'.format(z)
		dat = ascii.read(fich)
		m =dat['col1']
		dn=dat['col2']
		plt.close('all')
		plt.figure(1)
		plt.xlabel('log(M) $M_{odot}/h$')
		plt.ylabel(' dn/dlogM   $Mpc/h ^{-3} $')	
		plt.semilogy(np.log10(m),dn)
		plt.title(' Halo mass function ')	
		return

def z_function_halo(fich,zbin):
		dat = ascii.read(fich)
		x= dat['col1']
		y= dat['col2']
		z= dat['col3']
		zz= dat['col4']
		mass= dat['col5']

		plt.close('all')

		zmin=min(zz)
		zmax=max(zz)
		dz=(zmax-zmin)/float(N)
		for i in range(len(zz)):
			iz = int((z-zmin)/dz)
			

		plt.figure(1)
		plt.xlabel('log(M) $M_{odot}/h$')
		plt.ylabel(' dn/dlogM   $Mpc/h ^{-3} $')	
		plt.semilogy(np.log10(m),dn)
		plt.title(' Halo mass function ')	
		return



def histogram_bias2(fich,fichbias,N):
	plt.close('all')
	ii=1

	fich2=fich+'{:03d}'.format(ii)
	dat = ascii.read(fich2)
	xmin =dat['col1']
	xmax=dat['col2']
	snr=dat['col3']
	plt.step(xmin,snr,where='post',color='crimson',alpha=0.1,label='no bias')


	fich2=fichbias+'{:03d}'.format(ii)
	dat = ascii.read(fich2)
	xmin =dat['col1']
	xmax=dat['col2']
	snr=dat['col3']
	plt.step(xmin,snr,where='post',color='crimson',alpha=0.1,label='bias')

	for ii in range(2,N):
			fich2=fich+'{:03d}'.format(ii)
			dat = ascii.read(fich2)
			xmin =dat['col1']
			xmax=dat['col2']
			snr=dat['col3']
			plt.step(xmin,snr,where='post',color='crimson',alpha=0.1)

	for ii in range(2,N):
			fich2=fichbias+'{:03d}'.format(ii)
			dat = ascii.read(fich2)
			xmin =dat['col1']
			xmax=dat['col2']
			snr=dat['col3']
			plt.step(xmin,snr,where='post',color='deepskyblue',alpha=0.1)

	plt.title('Peak abundance histogram (averaged over {0} realizations)'.format(N))
	plt.xlabel('SNR')
	plt.ylabel('Peak number')
	plt.legend()
	plt.show()
	return

def histogram_bias(fich,fichbias,N):
	plt.close('all')

	
	if (N==1):
		dat = ascii.read(fich)
		xmin =dat['col1']
		xmax=dat['col2']
		nn=dat['col3']
		nbins = xmax

		plt.figure(1)
		plt.xlabel('SNR')
		plt.ylabel('Peak number')
		plt.step(xmin,nn,where='post',color='b',alpha=1,label='no bias')
		plt.title('Peak abundance histogram ')

		dat = ascii.read(fichbias)
		xmin =dat['col1']
		xmax=dat['col2']
		nn=dat['col3']
		nbins = xmax
		plt.step(xmin,nn,where='post',color='r',alpha=1,label='bias')
		plt.legend()
		plt.show()
		return

	if(N!=1):
		fich2=fich+'001'

		dat = np.loadtxt(fich2)
		xmin =dat[:,0]
		xmax=dat[:,1]
		nn=dat[:,2]

		nbins = xmax
		mean_snr = np.copy(nn)
		mean_snr_error = np.copy(nn)*0
		
		for ii in range(2,N):
			fich2=fich+'{:03d}'.format(ii)
			dat = np.loadtxt(fich2)
			xmin =dat[:,0]
			xmax=dat[:,1]
			nn=dat[:,2]
			nbins = xmax
			mean_snr[:]=mean_snr[:]+nn[:]
		
		mean_snr[:]=mean_snr[:]/float(N)

		for ii in range(1,N):
			fich2=fich+'{:03d}'.format(ii)
			dat = np.loadtxt(fich2)
			xmin =dat[:,0]
			xmax=dat[:,1]
			nn=dat[:,2]
			mean_snr_error = mean_snr_error +(nn-mean_snr)**2.

		mean_snr_error =(1./(N-1.)*mean_snr_error)**(0.5)
		nbin_snr=9
		data=np.linspace(1,5,nbin_snr)

		plt.errorbar(data[:(nbin_snr-1)]+0.25, mean_snr[:(nbin_snr-1)], yerr=mean_snr_error[:(nbin_snr-1)],fmt='+',color='crimson',alpha=0.7)
		#plt.step(data,mean_snr,where='post',color='crimson',alpha=1,label="NoBias")
		plt.step(data,mean_snr,where='post',color='crimson',alpha=1,label="Clustered galaxies")
		#plt.step(data,mean_snr,where='post',color='crimson',alpha=1,label="Detected sources")
		
        
		fich2=fichbias+'001'
		dat = np.loadtxt(fich2)
		xmin =dat[:,0]
		xmax=dat[:,1]
		nn=dat[:,2]
		nbins = xmax
		mean_snr = np.copy(nn)
		mean_snr_error = np.copy(nn)*0
		
		for ii in range(2,N):
			fich2=fichbias+'{:03d}'.format(ii)
			dat = np.loadtxt(fich2)
			xmin =dat[:,0]
			xmax=dat[:,1]
			nn=dat[:,2]
			nbins = xmax
			mean_snr=mean_snr+nn
		mean_snr=mean_snr/float(N)

		for ii in range(1,N):
			fich2=fichbias+'{:03d}'.format(ii)
			dat = np.loadtxt(fich2)
			xmin =dat[:,0]
			xmax=dat[:,1]
			nn=dat[:,2]
			mean_snr_error = mean_snr_error +(nn-mean_snr)**2.

		mean_snr_error =(1./(N-1.)*mean_snr_error)**(0.5)
		nbin_snr=9
		data=np.linspace(1,5,nbin_snr)

		plt.errorbar(data[:(nbin_snr-1)]+0.25, mean_snr[:(nbin_snr-1)], yerr=mean_snr_error[:(nbin_snr-1)],fmt='+',color='deepskyblue',alpha=0.7)
		#plt.step(data,mean_snr,where='post',color='deepskyblue',alpha=1,label="Bias")
		plt.step(data,mean_snr,where='post',color='navy',alpha=1,label="Random galaxies")
		#plt.step(data,mean_snr,where='post',color='deepskyblue',alpha=1,label="Full population")
		plt.title('Peak abundance histogram (averaged over {0} realizations)'.format(N))
		plt.xlabel('SNR')
		plt.ylabel('Peak number')
		plt.legend()
		#plt.savefig('Plots/PeakHists.pdf')
		return


##################################

###################################
def histogram_snr(N):
	plt.close('all')

	if (N==1):
		fich='../build/peakHist'
		dat = ascii.read(fich)
		xmin =dat['col1']
		xmax=dat['col2']
		nn=dat['col3']
		nbins = xmax

		plt.figure(1)
		plt.xlabel('SNR')
		plt.ylabel('Peak number')
		plt.step(xmin,nn,where='post',color='b',alpha=1)
		plt.title('Peak abundance histogram (only one realization)')
		plt.show()
	if(N!=1):
		fich='../build/dataMat_nbFilters2_N'+str(N)
		dat = np.loadtxt(fich)
		nbin_snr=9
		mean_dat=np.ones(nbin_snr)
		err_dat=np.ones(nbin_snr)*0.
		snr=np.linspace(1,5,nbin_snr)
		for i in range(nbin_snr):
			mean_dat[i]=np.mean(dat[:,i])
			for j in range(N):
				err_dat[i]= err_dat[i]+(dat[j,i]-mean_dat[i])**2.
				if(i==1): plt.step(snr,dat[j,0:nbin_snr],where='post',color='crimson',alpha=0.05)
			err_dat[i]=(1./((np.size(dat[:,i])-1.))*err_dat[i])**(1./2.)

		plt.figure(1)
		plt.errorbar(snr[:(nbin_snr-1)]+0.25, mean_dat[:(nbin_snr-1)],  yerr=err_dat[:(nbin_snr-1)],fmt='+',color='crimson',alpha=0.7)
		plt.step(snr,mean_dat,where='post',color='crimson',alpha=1)
		plt.title('Peak abundance histogram (averaged over {0} realizations)'.format(N))
		plt.xlabel('SNR')
		plt.ylabel('Peak number')
		plt.show()
	return


##################################

def map_lensing(fich,opt):
	plt.close('all')
	dat = np.loadtxt(fich)
	if(opt == 1):
		x = dat[:,0]
		y = dat[:,1]
		z = dat[:,2]
		N = int(len(z)**.5)
		x1=np.linspace(0.5,179.5,N)
		y1=np.linspace(0.5,179.5,N)
		X1,Y1 = np.meshgrid(x1,y1)  
  		z2=np.reshape(z,(N,N))

		plt.figure(1)
		plt.pcolor(X1,Y1,(z2),cmap=cm.terrain,vmin=z2.min(),vmax=z2.max())    
		cbar =plt.colorbar()
		plt.title( "Convergence map" )
		plt.xlabel(r' $\theta_x$ (arcmin)')
		plt.ylabel(r' $\theta_y$ (arcmin)')

	if(opt == 2):
		x = dat[:,0]
		y = dat[:,1]
		g1 = dat[:,2]
		g2 = dat[:,3]
		N = int(len(x)**.5)
		x1=np.linspace(0.5,179.5,N)
		y1=np.linspace(0.5,179.5,N)
		X1,Y1 = np.meshgrid(x1,y1)  
  		z1=np.reshape(g1,(N,N))
  		z2=np.reshape(g2,(N,N))

		plt.figure(1)
		plt.subplot(1,2,1)
		plt.pcolor(X1,Y1,(z1),cmap=cm.OrRd,vmin=z1.min(),vmax=z1.max())    
		cbar =plt.colorbar()
		plt.title( '$\\rm \gamma_1$ Map')
		plt.xlabel(r' $\theta_x$ (arcmin)')
		plt.ylabel(r' $\theta_y$ (arcmin)')

		plt.subplot(1,2,2)
		plt.pcolor(X1,Y1,(z2),cmap=cm.OrRd,vmin=z2.min(),vmax=z2.max())    
		cbar =plt.colorbar()
		plt.title( '$\\rm \gamma_2$ Map')
		plt.xlabel(r' $\theta_x$ (arcmin)')
		plt.ylabel(r' $\theta_y$ (arcmin)')

	plt.show()
	return;
#####################################

def map_peak(fich,fich2,lim):
	
	dat = np.loadtxt(fich2)
	snr = dat[:,0]
	x = dat[:,1]
	y = dat[:,2]
	xx=x[snr>lim]
	yy=y[snr>lim]

	dat = np.loadtxt(fich)
	x2 = dat[:,0]
	y2 = dat[:,1]
	z = dat[:,2]

	plt.close('all')

	N = int(len(z)**.5)
	z2=np.reshape(z,(N,N))
	x1=np.linspace(min(x2),max(x2),N)
	y1=np.linspace(min(y2),max(y2),N)
	X1,Y1 = np.meshgrid(x1,y1)     
    
	plt.figure(9)
	plt.pcolor(X1,Y1,(z2),cmap=cm.terrain,vmin=z2.min(),vmax=z2.max())    
	cbar =plt.colorbar()
	cbar.set_label(r'$\kappa$', rotation=-270)
 	plt.scatter(xx,yy, marker='o',s=200, edgecolor='red', alpha=1,zorder=10,facecolor='None')
	plt.title( 'Convergence Map and peak detection for SNR $>$ {0:.1f}'.format(lim))
	plt.xlabel(r' $\theta_x$ (arcmin)')
	plt.ylabel(r' $\theta_y$ (arcmin)')
	plt.axis([X1.min(), X1.max(), Y1.min(),Y1.max()]) 
	plt.show()

	return;

##########################################
def read_h(fich1,z,dz):
	dat = np.loadtxt(fich1)
	theta_x = dat[:,0]
	theta_y = dat[:,1]
	theta_z = dat[:,3]
	M = dat[:,4]
	Nc = dat[:,5]
	Ns = dat[:,6]
	theta_x[(theta_z<(z+dz/2.))&(theta_z>(z-dz/2.))] 
	theta_y[(theta_z<(z+dz/2.))&(theta_z>(z-dz/2.))] 
	return theta_x[(theta_z<(z+dz/2.))&(theta_z>(z-dz/2.))] ,theta_y[(theta_z<(z+dz/2.))&(theta_z>(z-dz/2.))],theta_z[(theta_z<(z+dz/2.))&(theta_z>(z-dz/2.))] ,M[(theta_z<(z+dz/2.))&(theta_z>(z-dz/2.))] ,Nc[(theta_z<(z+dz/2.))&(theta_z>(z-dz/2.))],Ns[(theta_z<(z+dz/2.))&(theta_z>(z-dz/2.))]

def read_catalogue_haloANDgalaxies(fich1,fich2,z,dz):
	dat = np.loadtxt(fich1)
	theta_x = dat[:,0]
	theta_y = dat[:,1]
	theta_z = dat[:,3]
	M = dat[:,4]
	theta_x[(theta_z>(z+dz/2.))|(theta_z<(z-dz/2.))] = -100
	theta_y[(theta_z>(z+dz/2.))|(theta_z<(z-dz/2.))] = -100
	M[(theta_z>(z+dz/2.))|(theta_z<(z-dz/2.))] = 0

	dat = np.loadtxt(fich2)
	thetag_x = dat[:,0]
	thetag_y = dat[:,1]
	thetag_z = dat[:,2]
	thetag_x[(thetag_z>(z+dz/2.))|(thetag_z<(z-dz/2.))] = -100
	thetag_y[(thetag_z>(z+dz/2.))|(thetag_z<(z-dz/2.))] = -100

	plt.close('all')

	plt.figure(1)
	plt.title( ' dark matter haloes')
	#f, axarr = plt.subplots(2, sharey=True)

	plt.scatter(theta_x ,theta_y,s=M/1e12*4, color='blue', alpha=0.5)
	plt.xlim([0,180])
	plt.ylim([0,180])
	plt.xlabel(r' $\theta_x$ (arcmin)')
	plt.ylabel(r' $\theta_y$ (arcmin)')
	print("nb halo :", np.size(theta_x[(theta_z<(z+dz/2.))&(theta_z>(z-dz/2.))]))
#	print("theta x  :", max(theta_x),min(theta_x))
#	print("theta y  :", max(theta_y),min(theta_y))
	plt.show()

	#plt.figure(2)	
	#plt.title( 'galaxy source positions')
	plt.scatter(thetag_x ,thetag_y,s=5, color='red', alpha=0.25)
	plt.xlim([0,180])
	plt.ylim([0,180])
	plt.xlabel(r' $\theta_x$ (arcmin)')
	plt.ylabel(r' $\theta_y$ (arcmin)')
	plt.show()
	print("nb gal  :", np.size(thetag_x[(thetag_z<(z+dz/2.))&(thetag_z>(z-dz/2.))]))
	#print("theta x  :", max(thetag_x),min(thetag_x))
	#print("theta y  :", max(thetag_y),min(thetag_y))
	return ;

##########################################
def read_zz(fich):

	dat = np.loadtxt(fich)
	red = dat[:,2]

	plt.close('all')
	plt.figure(1)

	zz2=np.linspace(0,2,31)
	ncat=np.copy(zz2)
	for i in range(np.size(zz2)):
		if( i < np.size(zz2)-1):
			ncat[i]= sum(red[red<zz2[i]])
	

	#plt.yscale('log')
	plt.figure(1)
	plt.title( 'Redshift histogramm of galaxy sources')
	plt.ylabel(r'\LARGE{$\rm N_{sources}$}')
	plt.xlabel(r'Redshift')
	plt.plot(zz2,ncat/(180.*180.),alpha=0.80,color='green')
	plt.show()
	return

def Hoekstra(a,b,c,aa):
	zz=np.linspace(0,2,30)
	nz=aa*(zz**a+zz**(a*b))/(zz**b+c)
	return nz,zz

##########################################
def read_catalogue_galaxies1(fich):

	#galcat_dir='./HoekstraNzQuantities/'  
	#import HoekstraPop as hp  
	#Nz3, zz3 = hp.HoekstraNz(galcat_dir) 
	#Nz3, zz3 = Hoekstra(0.497,12.643,0.381,4068.19) 

	dat = np.loadtxt(fich)
	red = dat[:,2]

	plt.close('all')
	plt.figure(1)

	zz2=np.linspace(0,max(red),301)
	nn=np.copy(zz2)
	ncat=np.copy(zz2)
	for i in range(np.size(zz2)):
		if( i < np.size(zz2)-1):
			ncat[i]= np.size(red[((red<zz2[i+1])&(red>zz2[i]))])

	zz3=np.linspace(0,5,301)
	for i in range(np.size(zz3)):
		nn[i]=nnz(zz3[i])

	#plt.yscale('log')
	ng=sum(ncat)
	plt.figure(1)
	plt.title( 'Redshift histogramm of galaxy sources')
	plt.ylabel(r'\LARGE{$\rm N_{sources}$ per arcmin2} ')
	plt.xlabel(r'Redshift')
	plt.plot(zz2,ncat/180./180.,alpha=0.80,color='green')
	plt.plot(zz3,nn,alpha=0.80,color='crimson')
#	plt.plot(zz3,Nz3/((2.*60.)**2.),alpha=0.80,color='blue')


	#Ng2=sum(Nz3)
	print("\n nb galaxies : {0}, density nb/arcmin2 n={1} \n".format(ng,ng/((2.*60.)**2)))
	print("\n nb galaxies normal : {0}, density nb/arcmin2 n={1} \n".format(sum(nn),sum(nn)/sum(nn)*47.))
	#print("\n nb Hoekstra : {0}, density nb/arcmin2 n={1} \n".format(Ng2,Ng2/((2.*60.)**2)))
	plt.show()


def ngal(zz,dz):
	zz1=np.linspace(0,zz,dz)
	nn1=np.copy(zz1)
	for i in range(np.size(zz1)):
		nn1[i]=nnz(zz1[i])

	plt.figure(3)
	plt.plot(zz1,nn1,alpha=0.80,color='crimson')
	print("\n nb galaxies normal : {0}\n".format(sum(nn1)/dz))
	plt.show()
	return

def read_nz(fich,fich2):


	dat = np.loadtxt(fich)
	red = dat[:,2]

	plt.close('all')
	plt.figure(1)

	plt.hist(red,bins=30,color = 'limegreen', edgecolor = 'green',alpha=0.60,normed=True)

	dat = np.loadtxt(fich2)
	red = dat[:,2]
	plt.hist(red,bins=30,color = 'red', edgecolor = 'red',alpha=0.60,normed=True)

	plt.title( 'Redshift histogramm of galaxy sources')
	plt.ylabel(r'\LARGE{$\rm N_{sources}$}')
	plt.xlabel(r'Redshift')

	return

def read_catalogue_galaxies2(fich,fich2):
	dat = np.loadtxt(fich)
	red = dat[:,2]

	dat = np.loadtxt(fich2)
	red2 = dat[:,2]+0.03

	plt.close('all')
	plt.figure(1)

	plt.hist(red,bins=100,color = 'limegreen', edgecolor = 'green',alpha=0.60)
	plt.hist(red2,bins=100,color = 'deepskyblue', edgecolor = 'deepskyblue',alpha=0.90)

	zz=np.linspace(min(red),max(red))
	nn=np.copy(zz)
	for i in range(np.size(zz)):
		nn[i]=nnz(zz[i])*np.size(red2)
	plt.plot(zz,nn,alpha=0.30,color='crimson')
	plt.yscale('log')
	plt.title( 'Redshift histogramm of galaxy sources')
	plt.ylabel(r'\LARGE{$\rm N_{sources}$}')
	plt.xlabel(r'Redshift')

	return 

########################################

def read_catalog_halo_NM(fich):
	dat = np.loadtxt(fich)
	plt.close('all')
	theta_x = dat[:,0]
	theta_y = dat[:,1]
	ww = dat[:,2]
	zz = dat[:,3]
	mm = dat[:,4]
	nc = dat[:,5]
	ns = dat[:,6]
	rv = dat[:,7]

	Ng =nc*(1+ns)
	mm[zz>0.6]=0
	mm[zz<0.4]=0
	Ng[zz<0.4]=0
	Ng[zz>0.6]=0

	plt.figure(3)
	plt.loglog(mm,Ng,'.')
	plt.show()
	return

def read_catalog_halo(fich,opt):
	dat = np.loadtxt(fich)
	plt.close('all')
	theta_x = dat[:,0]
	theta_y = dat[:,1]
	ww = dat[:,2]
	zz = dat[:,3]
	mm = dat[:,4]
	plt.figure(1)

	if(opt==1):
		MIN=min(mm)
		MAX=max(mm)
		nbins =50
		nbins = np.linspace(np.log10(MIN), np.log10(MAX), nbins)
		plt.title( 'Halo mass histogramm')
		res =plt.hist(np.log10(mm), bins = nbins,color = 'limegreen', edgecolor = 'green',alpha=0.60, stacked=True)
		plt.yscale('log')
		ylab=r'{ $ \rm N_{halo} $}'
		xlab=r'{ $ \rm log( M_{halo} ( \rm M_{\odot}/h ) ) $}'
		axes = plt.gca()
		plt.xlabel(xlab)
		plt.ylabel(ylab)
	if(opt==2):
		res =plt.hist(zz, bins = 30,color = 'limegreen', edgecolor = 'green',alpha=0.60, stacked=True)
		plt.yscale('log')
		plt.title( 'Halo redshift histogramm')	
		ylab=r' $ \rm N_{halo} $'
		xlab=r'redshift'
		plt.xlabel(xlab)
		plt.ylabel(ylab)
	return ;
######################

def rewrite_Galaxy_catalogue(fich2,fich3):
	dat = np.loadtxt(fich2)
	theta_x = dat[:,0]
	theta_y = dat[:,1]
	red = dat[:,2]
	k  = dat[:,3]
	g1 = dat[:,4]
	g2 = dat[:,5]
	nn=np.size(k)
	m=0 #0.1 ##bias used routine   m= -0.00856*d  ## champ densite   k mass/surface 
	g1n =np.copy(g1)
	g2n =np.copy(g2)
	file3 = open(fich3,"w") 
	file2 = open(fich2,"r") 
	aa=0
	for i in range(nn):
		x = theta_x[i]
		for j in range(nn):
			y = theta_y[j]


	for line in file2: 
		file3.write(line)
		aa=aa+1
		if(aa==12): break

	file2.close()

	for i in range(nn):
		file3.write('\t {0:2.3f} \t {1:3.3f} \t {2:4.5f} \t {3:4.5f} \t {4:4.5f} \t {5:4.5f} \n'.format(dat[i,0],dat[i,1],dat[i,2],dat[i,3],g1[i],g2[i]))
  	file3.close()
	return ;

##########################
######################

def redefine_Galaxy_catalogue(fich2,fich3,zz,dz):
	datg = np.loadtxt(fich2)
	xg = datg[:,0]
	yg = datg[:,1]
	zg = datg[:,2]
	k  = datg[:,3]
	g1 = datg[:,4]
	g2 = datg[:,5]
	ng=np.size(k)

	dath = np.loadtxt(fich3)
	xh = dath[:,0]
	yh = dath[:,1]
	wh = dath[:,2]
	zh = dath[:,3]
	mh = dath[:,4]
	nh=np.size(zh)

	xh2=xh[((zh>zz-dz)&(zh<zz+dz))]
	yh2=yh[((zh>zz-dz)&(zh<zz+dz))]

	xg2=xg[((zg>zz-dz)&(zg<zz+dz))]
	yg2=yg[((zg>zz-dz)&(zg<zz+dz))]

	print("\nnb galaxies : {0}, density nb/arcmin2 n={1} \n".format(ng,ng/(180.**2)))
	print("nb halo : {0}, density nb/arcmin2 n={1} \n".format(nh,nh/(180.**2)))

	plt.close('all')

	plt.figure(1)
	plt.xlim([-10,190])
	plt.ylim([-10,190])
	plt.plot(xh2,yh2,'.',color='crimson',alpha=0.2)
	plt.show()

	plt.figure(2)
	plt.xlim([-10,190])
	plt.ylim([-10,190])
	plt.plot(xg2,yg2,'.',color='deepskyblue',alpha=0.2)
	plt.show()
	return

