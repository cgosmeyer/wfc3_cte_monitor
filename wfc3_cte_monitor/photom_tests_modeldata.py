#! /usr/bin/env python

""" One-off script for comparing iraf.phot to photutils.aperture_photometry
on a simulated 2-d gaussian star.

Authors:

    C.M. Gosmeyer, Fall 2016

"""


from __future__ import print_function

from pyraf import iraf
from iraf import noao, digiphot, daophot

from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table
import os
import shutil
import pylab
import matplotlib.pyplot as plt
import numpy as np
from photutils_plus.apphot import apphot

from photutils.datasets import make_gaussian_sources


path_to_out = '/Users/cgosmeyer/Desktop/wfc3_repos/reports/WFC3_CTEupdate_2016/tests/photom_test_simulatedsource'
image = os.path.join(path_to_out, 'images', 'teststar.fits')


#-------------------------------------------------------------------------------# 

def do_python_phot():
    """ Runs Python photometry (from my `apphot` funtion in repo 
    detector_tools.photometry.apphot)
    """
    coofile = os.path.join(path_to_out, 'coo_files', 'teststar_python.coo')
    outname = os.path.join(path_to_out, 'mag_files', 'teststar_python_exact')
    #apertures = [1,2,3,4,5,10,12,15,20]
    apertures = np.arange(0.5,15.5,0.5) #[0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4.0, 4.5, 5.0]
    apphot(imagename=image, ext=0, coofile=coofile, colnames=['extr_xpix', 'extr_ypix'], 
           outname=outname,
           ap_radii=apertures, sep_out=False, backmethod='mean',
           backglobal=False, method='exact', subpixels=None, 
           pixelwise_error=None, effective_gain=None, mask=None,
           annulus=4., dannulus=2., verbose=False, centroid=True)
    
#-------------------------------------------------------------------------------# 

def do_iraf_phot():
    """ Runs IRAF photometry with `iraf.phot`.
    """ 
    coofile = os.path.join(path_to_out, 'coo_files', 'teststar_iraf.coo')
    outname = os.path.join(path_to_out, 'mag_files', 'teststar_iraf_centroid.mag')
    #apertures='1,2,3,4,5,10,12,15,20'
    aperture_list = np.arange(0.5,15.5,0.5) 
    # Make the aperture list IRAF-compatible.
    apertures = ''
    for ap in aperture_list:
        apertures += str(ap) + ','
    apertures[:-1]    

    print(" ")
    print(image)
    iraf.phot.unlearn()         
    iraf.phot(image=image+'[0]', 
              interactive='no', 
                verify='no', 
                coords=coofile, 
                output=outname,
                #fwhmpsf=2.5,
                #itime=exptime, 
                calgorithm='centroid', 
                #cbox=5.0,
                #annulus=12.0,
                #dannulus=10.0,
                apertures=apertures)    

#-------------------------------------------------------------------------------# 

def make_star(shape=(50, 50)):
    """ Creates a 2-d gaussian star and writes to a FITS file.
    """
    table = Table()
    table['amplitude'] = [100]
    table['x_mean'] = [26]
    table['y_mean'] = [26]
    table['x_stddev'] = [3.0]
    table['y_stddev'] = [3.0]
    table['theta'] = np.array([0.]) * np.pi / 180.
    
    star = make_gaussian_sources(shape, table)
    
    hdu = fits.PrimaryHDU(star)
    hdu.writeto(os.path.join(path_to_out, 'images', 'teststar.fits'))
    
    return star
    
#-------------------------------------------------------------------------------# 

def make_oneimage(shape=(50, 50)):
    """ Creates an image of all 1's, written to a FITS file.
    """
    im = np.ones(shape)
    hdu = fits.PrimaryHDU(im)
    hdu.writeto(os.path.join(path_to_out, 'images', 'ones.fits'))
    
    #print(im)

#-------------------------------------------------------------------------------# 

def plot_iraf_over_python(mag_iraf, mag_python, outname=''): 
    """ Plots th IRAF photometry over the Python photometry.

    """
    
    given_aps = np.arange(0.5,15.5,0.5) 
    
    # Read in iraf mag file.
     

    data_iraf = ascii.read(os.path.join(path_to_out, 'mag_files', mag_iraf))
    
    print(data_iraf.colnames)
    
    fluxes_iraf = np.array(data_iraf['col4'])
    areas_iraf = np.array(data_iraf['col3'])   


    # do this by hand. drrr astropy drrr
    """
    fluxes_iraf = np.array([100., 319.8979, 688.0577, 1156.402, 1663.54, 2211.148, 2799.446, 3312.89, 3796.934, 4238.898,
     	4606.379, 4873.84, 5106.276, 5276.05, 5400.859, 5486.187, 5486.187, 5550.457, 5589.791, 5615.728, 5632.056,
        5642.155, 5647.584, 5650.978, 5652.829, 5653.84, 5654.346, 5654.625, 5654.751, 5654.814, 5654.843])
    areas_iraf = np.array([1., 3.343146, 7.343146, 13.11146, 19.79775, 28.38807, 38.85737, 50.42309, 63.2675, 79.12653,
    	95.8515, 113.2414, 132.8508, 154.5519, 176.8093, 200.6849, 227.3683, 254.9529, 283.2889, 314.7403,
    	347.3922, 380.3317, 415.136, 452.8061, 490.9447, 530.6741, 573.2782, 616.6547, 660.3874, 706.8707])
    """
    print("IRAF COMPLETE")

    # Read in python mag file.
    data_python = ascii.read(os.path.join(path_to_out, 'mag_files', mag_python))
    fluxes_python = np.array(data_python['aperture_sum'])
    areas_python = np.pi*(np.array(given_aps))**2
    print(fluxes_python)
    print(areas_python)
    

    
    # Plot iraf and python photometry together vs given aperture.
    fig=pylab.figure(figsize=(12, 10))
    pylab.scatter(given_aps, fluxes_iraf,  color='blue', s=75, alpha=0.75, label='iraf')
    pylab.scatter(given_aps, fluxes_python,  color='orange', s=75, alpha=0.75, label='python')
    pylab.xlabel('Given Aperture Radius [pxls]', fontsize=22, weight='bold')
    pylab.ylabel('Flux', fontsize=22, weight='bold') 
    pylab.tick_params(axis='both', which='major', labelsize=20)
    pylab.title('IRAF.PHOT and PYTHON.PHOTUTILS vs Given Ap Size', fontsize=22)
    pylab.axis('tight')
    pylab.minorticks_on()
    pylab.ylim([0, 5000])
    pylab.xlim([0,5])
    pylab.legend(loc='lower right')
    pylab.savefig(os.path.join(path_to_out, 'plots', '{}_fluxVSapsize.png'.format(outname)), bbox_inches='tight') 

    # Plot ratio vs given aperture
    fig=pylab.figure(figsize=(12, 10))
    pylab.axhline(y=1.0, linewidth=2, linestyle='--', color='grey')
    pylab.axvline(x=3.0, linewidth=2, linestyle='--', color='orange')
    pylab.scatter(given_aps, fluxes_iraf/fluxes_python,  color='blue', s=75)
    #pylab.scatter(given_aps, fluxes_python,  color='orange', label='python')
    pylab.xlabel('Given Aperture Radius [pxls]', fontsize=22, weight='bold')
    pylab.ylabel('IRAF/Python Flux', fontsize=22, weight='bold') 
    pylab.tick_params(axis='both', which='major', labelsize=20)
    #pylab.title('IRAF.PHOT / PYTHON.PHOTUTILS Ratio vs Given Ap Size', fontsize=22)
    pylab.axis('tight')
    pylab.xlim([0,10])
    pylab.minorticks_on()
    pylab.savefig(os.path.join(path_to_out, 'plots', '{}_ratio_vs_given.png'.format(outname)), bbox_inches='tight') 
    
    # Plot ratio vs real aperture
    fig=pylab.figure(figsize=(12, 10))
    pylab.axhline(y=1.0, linewidth=2, linestyle='--', color='grey')
    pylab.axvline(x=2204, linewidth=2, linestyle='--', color='orange')
    pylab.scatter(fluxes_iraf, areas_iraf/areas_python, color='blue', s=75)
    pylab.annotate('Given Ap Radius = 3 pix', xy=(2280, 1.15), color='orange', size=20, weight='bold')
    pylab.xlabel('IRAF Flux [e-]', fontsize=22, weight='bold')
    pylab.ylabel('IRAF/Python Aperture Area', fontsize=22, weight='bold') 
    pylab.tick_params(axis='both', which='major', labelsize=20)
    #pylab.title('IRAF.PHOT / PYTHON.PHOTUTILS Area vs IRAF Flux', fontsize=22)
    pylab.axis('tight')
    #pylab.xlim([0,15])
    pylab.minorticks_on()
    pylab.savefig(os.path.join(path_to_out, 'plots', '{}_arearatio_vs_irafflux.png'.format(outname)), bbox_inches='tight') 
    
    
    # Plot flux IRAF vs actual radius IRAF on top of flux python vs radius python
    fig=pylab.figure(figsize=(12, 10))
    pylab.scatter(np.sqrt(areas_iraf/np.pi), fluxes_iraf,  color='blue', label='iraf fluxes and radii', s=75)
    pylab.scatter(given_aps, fluxes_python,  color='orange', label='python fluxes and radii', s=75)
    pylab.xlabel('Aperture Radius [pix]', fontsize=22, weight='bold')
    pylab.ylabel('Flux', fontsize=22, weight='bold') 
    pylab.tick_params(axis='both', which='major', labelsize=20)
    pylab.title('Flux vs Actual Radius', fontsize=22)
    pylab.axis('tight')
    pylab.xlim([0,5])
    pylab.ylim([0,5000])
    pylab.minorticks_on()
    pylab.legend(loc='lower right')
    pylab.savefig(os.path.join(path_to_out, 'plots', '{}_flux_vs_radius.png'.format(outname)), bbox_inches='tight')     
    
#-------------------------------------------------------------------------------# 
#-------------------------------------------------------------------------------# 

if __name__=='__main__':
    #star = make_star()
    #pylab.imshow(star)

    #make_oneimage(shape=(50, 50))


    #do_iraf_phot()
    #do_python_phot()

    plot_iraf_over_python(mag_iraf='teststar_iraf_centroid.mag', mag_python='teststar_python_exact.mag', outname='')


