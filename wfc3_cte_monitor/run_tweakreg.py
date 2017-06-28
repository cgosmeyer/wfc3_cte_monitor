#! /usr/bin/env python

"""
Module for running tweakreg.TweakReg and updating WCS in headers.
Runs tweakreg on the FLCs and applied the same WCS solution on the
corresponding FLTs. This is especially important in later data because
CTI gets to be so bad that tweakreg may not find enough sources in 
FLTs with short exposure time in the narrow-band filter.

This is meant to be the first prep step after the data is downloaded to 
'data/newdata', to be run first before all else. 

Author:

    C.M. Gosmeyer, April 2016
    With R. Avila's help!

Use:

    Be in the 'data/newdata' directory.

    >>>python run_tweakreg.py

    This will run tweakreg on all the FLCs in the directory and also
    correct the corresponding FLTs with the same solution.

Outputs:
    
    In the 'data/newdata' directory:
    * A WCS-corrected FLC.
    * A WCS-corrected FLT.

    Moved to the 'outputs/tweakreg' directory:
    * A matched tweakreg catalog ipppssoot_xy_catalog.match
    * Another matched tweakreg catalog ipppssoot_flc_catalog_fit.match
      (we actually use the first one, but retain this one just in case)
    * A WCS solution file ipppssoot_hlet.fits, which is used to tweak the
      FLT.
    * The tweakreg log ipppssoot_flc_tweakreg.log.
    * The diagnoistic plots hist2d_ipppssoot_flc.png
                            residuals_ipppssoot_flc.png
                            vector_ipppssoot_flc.png

    All other tweakreg outputs are removed.

"""

from __future__ import print_function

import glob
import multiprocessing
import os
import shutil
from astropy.io import fits
from drizzlepac import tweakreg
import stwcs.wcsutil.headerlet as hlet

import config


#-------------------------------------------------------------------------------# 

def run_tweakreg(imagename, updatehdr=True, sort_files=False):
    """ Runs tweakreg.TweakReg and updates WCS in headers.
    Actually first runs on the FLCs and applied the same WCS solution on
    the corresponding FLTs. This is especially important in later data 
    because CTI gets to be so bad that tweakreg may not find enough sources 
    in FLTs with short exposure time in the narrow-band filter.

    Parameters:
        imagename : string
            Name of the image.
        updatehdr : {True, False}
            Set on to update WCS in header. Usually should have this on
            unless testing stuff!!
        sort_files : {True, False}
            Turn on if want to sort the output files into proposal 
            directories.

    Returns:
        nothing.

    Outputs:
        In the 'data/newdata' directory:
        * A WCS-corrected FLC.
        * A WCS-corrected FLT.

        Moved to the 'outputs/tweakreg' directory:
        * A matched tweakreg catalog ipppssoot_xy_catalog.match
        * Another matched tweakreg catalog ipppssoot_flc_catalog_fit.match
          (we actually use the first one, but retain this one just in case)
        * A WCS solution file ipppssoot_hlet.fits, which is used to tweak 
          the FLT.
        * The tweakreg log ipppssoot_flc_tweakreg.log.
        * The diagnoistic plots hist2d_ipppssoot_flc.png
                                residuals_ipppssoot_flc.png
                                vector_ipppssoot_flc.png

        All other tweakreg outputs are removed.
    
    Notes from R. Avila: 

    By setting headerlet=True, attach=False, creates a file called ipppssoot_hlet.fits
    that contains the new solution. you then run a task called apply_headerlet 
    which can apply the headerlet to the FLT
    """
    if updatehdr:
        headerlet = True
    else:
        headerlet = False
    
    # Obtain the filter, target, expotime, and proposal number.
    priheader = fits.getheader(imagename, 0) 

    targname = (priheader['TARGNAME']).lower()
    proposid = str(priheader['PROPOSID'])
    filt = (priheader['FILTER']).lower()
    exptime = int(priheader['EXPTIME'])
    chinject = str(priheader['CHINJECT'])

     # Set the path to the proposal directory. 
    path_to_data = config.path_to_data
    proposaldir = os.path.join(path_to_data, proposid)

    # Get rootname.
    rootname = imagename.split('_fl')[0]
    imagenames = [rootname+'_flt.fits', rootname+'_flc.fits']

    # Skip the charge injections of CONT from 12348 (tweakreg can't find any sources)
    if chinject == 'CONT':
        for imagename in imagenames:
            shutil.move(imagename, proposaldir)
            print("Charge injection CONT!")
            print("Moved {} to {}".format(imagename, proposaldir))
        return ''

    # Get path to the master catalogs.
    path_to_data = config.path_to_data
    path_to_mastercat = os.path.join(path_to_data, 'master_cat')

    # Change parameters based on target, filter, and exptime.
    if '6791' in targname:
        print('ngc6791')
        master_image = os.path.join(path_to_mastercat, 'out_n6791_master_drc_sci.fits')
        master_cat = os.path.join(os.path.join(path_to_mastercat, 'final'), 'ngc6791_master.cat')
        if filt == 'f606w':
            print('f606w')
            if exptime <= 60: #short exp
                print(exptime)
                conv_width = 3.5
                threshold = 25
                peakmax = 50000
            elif exptime > 60: #long exp
                conv_width = 3.5
                threshold = 200
                peakmax = 50000
        elif filt == 'f502n':
            print('f502n')
            if exptime <= 60: #short exp
                print(exptime)
                conv_width = 3.5
                threshold = 5
                peakmax = 50000
            elif exptime > 60: #long exp
                print(exptime)
                conv_width = 3.5
                threshold = 100
                peakmax = 50000

    elif '104' in targname:
        print('ngc104')
        master_image = os.path.join(path_to_mastercat, 'out_n104_master_drc_sci.fits')
        master_cat = os.path.join(os.path.join(path_to_mastercat, 'final'), 'ngc104_master.cat')
        if filt == 'f606w':
            print('f606w')
            if exptime <= 60: #short exp
                print(exptime)
                conv_width = 3.5
                threshold = 25
                peakmax = 50000
            elif exptime > 60: #long exp
                print(exptime)
                conv_width = 3.5
                threshold = 200
                peakmax = 50000
        elif filt == 'f502n':
            print('f502n')
            if exptime <= 60: #short exp
                print(exptime)
                conv_width = 3.5
                threshold = 5
                peakmax = 50000
            elif exptime > 60: #long exp
                print(exptime)
                conv_width = 3.5
                threshold = 100
                peakmax = 50000
    elif '6583' in targname:
        print('ngc6583')
        master_image = os.path.join(path_to_mastercat, 'out_n6583_master_drc_sci.fits')
        master_cat = os.path.join(os.path.join(path_to_mastercat, 'final'), 'ngc6583_master.cat')
        conv_width = 3.5
        threshold = 25
        peakmax = 50000


    tweakreg.TweakReg(imagename,
                      updatehdr=updatehdr,
                      headerlet=headerlet,
                      attach=False,
                      wcsname='TWEAK',
                      reusename=True,
                      interactive=False,
                      conv_width=conv_width,
                      threshold=threshold, 
                      minobj=5,
                      peakmax=peakmax, 
                      refcat=master_cat,
                      refimage=master_image,
                      refxcol=2,
                      refycol=3,
                      refxyunits='pixels',
                      runfile=imagename+'_tweakreg.log')

    # Apply the WCS change to the FLT using the *hlet file.
    if updatehdr:
        hlet.apply_headerlet_as_primary(imagenames[0],
                                   '{}_hlet.fits'.format(rootname),
                                    attach=False,
                                    archive=False)
    
# how get rid of this warning?
#WARNING: AstropyDeprecationWarning: Deletion of non-existent keyword 'DISTNAME': In a future Astropy version Header.__delitem__ may be changed so that this raises a KeyError just like a dict would. Please update your code so that KeyErrors are caught and handled when deleting non-existent keywords. [astropy.io.fits.header]


    # Move the files to the proposal directory.
    if sort_files:
        for imagename in imagenames:
            shutil.move(imagename, proposaldir)
            print("Moved {} to {}".format(imagename, proposaldir))

    path_to_outputs = config.path_to_outputs
    path_to_tweakreg = os.path.join(path_to_outputs, 'tweakreg')
    # Move the *.match to outputs/tweakreg.
    match_files = glob.glob('*'+rootname+'*.match')
    for match_file in match_files:
        if os.path.isfile(os.path.join(path_to_tweakreg, match_file)):
            os.remove(os.path.join(path_to_tweakreg, match_file))   
        shutil.move(match_file, path_to_tweakreg)
        print("Moved {} to {}".format(match_file, path_to_tweakreg))

    # Move the hlet file to outputs/tweakreg
    if updatehdr:
        if os.path.isfile(os.path.join(path_to_tweakreg, 
        '{}_hlet.fits'.format(rootname))):
            os.remove(os.path.join(path_to_tweakreg, 
            '{}_hlet.fits'.format(rootname)))    
        shutil.move('{}_hlet.fits'.format(rootname), path_to_tweakreg)
        print("Moved {} to {}".format('{}_hlet.fits'.format(rootname), path_to_tweakreg))

    # Move the tweakreg png files to outputs/tweakreg
    png_files = glob.glob('*'+rootname+'*.png')
    for png_file in png_files:
        if os.path.isfile(os.path.join(path_to_tweakreg, png_file)):
            os.remove(os.path.join(path_to_tweakreg, png_file))   
        shutil.move(png_file, path_to_tweakreg)
        print("Moved {} to {}".format(png_file, path_to_tweakreg))

    # Move the log files to outputs/tweakreg
    if os.path.isfile(os.path.join(path_to_tweakreg, 
        '{}_flc.fits_tweakreg.log'.format(rootname))):
        os.remove(os.path.join(path_to_tweakreg, 
            '{}_flc.fits_tweakreg.log'.format(rootname)))
    shutil.move('{}_flc.fits_tweakreg.log'.format(rootname), path_to_tweakreg)
    print("Moved {} to {}".format('{}_flc.fits_tweakreg.log'.format(rootname), 
          path_to_tweakreg))
    


#-------------------------------------------------------------------------------# 

def run_tweakreg_parallel():
    """ Performs ``run_tweakreg`` in parallel on all FLC files in the
    current working directory.

    This does not appear to work. boo. Get error:

    The process has forked and you cannot use this CoreFoundation functionality safely. You MUST exec().
Break on __THE_PROCESS_HAS_FORKED_AND_YOU_CANNOT_USE_THIS_COREFOUNDATION_FUNCTIONALITY___YOU_MUST_EXEC__() to debug.
    """
    flc_list = glob.glob('*flc.fits')
    print(flc_list)

    pool = multiprocessing.Pool(processes=2)
    pool.map(run_tweakreg, flc_list)
    pool.close()
    pool.join()

    print("Parallel processing of run_tweakreg complete.")

    # Clean up the directory. 
    unwanted_files = glob.glob('*coo') + \
                     glob.glob('*list')
    for unwanted_file in unwanted_files:
        os.remove(unwanted_file)
        print("Removed {}".format(unwanted_file))


#-------------------------------------------------------------------------------# 

def run_tweakreg_main(updatehdr=True):
    """ Performs ``run_tweakreg`` on all FLC files in the current working 
    directory.

    Parameters:
        updatehdr : {True, False}
            Set on to update WCS in header. Usually should have this on unless
            testing stuff!!
    """
    flc_list = glob.glob('*flc.fits')
    print(flc_list)

    for flc in flc_list:
        run_tweakreg(flc, updatehdr=updatehdr)

    print(" ")
    print("Loop over run_tweakreg complete.")

    # Clean up the directory. 
    unwanted_files = glob.glob('*coo') + \
                     glob.glob('*list')
    for unwanted_file in unwanted_files:
        os.remove(unwanted_file)
        print("Removed {}".format(unwanted_file))

    print(" ")
    print("Main of run_tweakreg complete.")


#-------------------------------------------------------------------------------# 
# The Main.
#-------------------------------------------------------------------------------# 

if __name__=='__main__':
    run_tweakreg_main()
    