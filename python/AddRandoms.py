import numpy as np
import sys
import linecache
import os
import HoekstraPop as hp
    
def AddRandoms(galcat_dir, filename, bin_edges, savestub, dz):
    # read galaxy catalog
    galcat = np.loadtxt(galcat_dir+filename)
    fullzs = galcat[:,2]
    
    # read extreme positions 
    xmin, xmax = 0., 180. #np.min(galcat[:,0]), np.max(galcat[:,0])
    ymin, ymax = 0., 180. #np.min(galcat[:,1]), np.max(galcat[:,1])
    
    # apply cutoff
    print ' > Adding randoms to galaxy catalog {}.'.format(filename)
    
    NumberCounts, popz = hp.HoekstraNz(galcat_dir)
    
    zmid = [rf-(rf-lf)/2 for lf, rf in zip(bin_edges, bin_edges[1:])]
    
    if not np.all(np.array(popz) == np.array(zmid)):
        print '/!\ REDSHIFT MISMATCH! /!\ '
    
    for nb_rand, (lf, rf) in zip(NumberCounts,zip(bin_edges, bin_edges[1:])):
        print '   > Working on z bin [{},{}]'.format(lf,rf)
        nb_popgal = len(np.where((fullzs>lf) & (fullzs<=rf))[0])
        needed_rands = int(nb_rand) - nb_popgal
        if needed_rands < 0:
            print '      > Redshift bin {}:\t no randoms added ({} populated gals vs {} in Hoekstra pop)'.format(
                    (lf,rf), nb_popgal, nb_rand)
        xs, ys = np.random.rand(2,needed_rands) * np.array([xmax, ymax]).reshape(-1,1) + np.array([
                        xmin, ymin]).reshape(-1,1)
        zs = np.ones(needed_rands) * rf-(rf-lf)/2
        dummy = np.zeros(needed_rands)
        newrows = np.vstack((xs, ys, zs, dummy, dummy, dummy)).T
        galcat = np.vstack((galcat,newrows))
    print ' > Randoms added: catalog size now {}'.format(galcat.shape[0])
    
    idnb = '_'
    for digit in [char for char in filename[-3:] if char.isdigit()]:
        idnb += digit
    np.savetxt(galcat_dir+savestub+idnb, galcat)


def main():
    """Apply bias to Camelus-generated galaxy catalogs. Syntax:
    
    > python AddRandoms.py path/to/catalogfolders/ galcat dz output_fullpop zmin zmax
    
    """
    galcat_dir, galcat_name, dz = sys.argv[1], sys.argv[2], float(sys.argv[3])
    savestub = sys.argv[4]
    zmin = float(sys.argv[5])
    zmax = float(sys.argv[6])
    bin_edges = np.arange(zmin,zmax+dz,dz)
    stub_length = len(galcat_name)
    galcat_files = [catfile for catfile in os.listdir(galcat_dir) 
                    if galcat_name == catfile[:stub_length]]
    _ = [AddRandoms(galcat_dir, filename, bin_edges, savestub, dz) for 
         filename in galcat_files]
    
if __name__ == "__main__":
    main()
