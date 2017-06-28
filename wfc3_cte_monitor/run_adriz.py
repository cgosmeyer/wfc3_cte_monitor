#! /usr/bin/env python

"""
Should be second prep step.

Runs astrodrizzle to 

[1] Insert the 'mdrizsky' header keyword. Potentially will use this 
    background estimate in analysis.

[2] Create a cosmic ray mask. 


Use:

    Run after the 'run_tweakreg.py' step. Be in the 'data/newdata' 
    directory. 

    >>> python run_adriz.py

Notes:

    [1] This will NOT work on a file if there are no other files in 
    the directory containing same patch of sky. So always best to 
    run this over the full visit.
    It will also run out of memory of you feed it too many files at once
    so the script automatically sorts files by proposal and visit and
    feeds AstroDrizzle the files in visit-sized chunks.

    [2] This assumes each FLC has a corresponding FLT. The same cosmic-ray
    mask generated for the FLC will be applied to the FLT, so only
    run over FLCs. (have verified that each file type generates
    same mask)

"""
from __future__ import print_function
from drizzlepac import astrodrizzle
from astropy.io import fits
from collections import OrderedDict
import config
import glob
import os
import shutil


#-------------------------------------------------------------------------------# 

def run_adriz(flc_files):
    """ Runs AstroDrizzle in order to create cosmic ray masks and to obtain
    an estimate of the global sky background in new keyword 'mdrizsky'.
    """

    # Assumes first 6 letters of a visit's images are all the same.
    common = flc_files[0][:6]
    search = '{}*flc.fits'.format(common)
    print('search term for AstroDrizzle: {}'.format(search))
    # Check that really did capture all files with the search term. 
    # You'll need do some hacking if this error should occur.
    for flc in flc_files:
        if common not in flc:
            print("Error!! File {} does not match rest of visit with {}.".format(flc, common))
            return

    astrodrizzle.AstroDrizzle(search,
                              runfile='',
                              output='',
                              preserve=False,
                              updatewcs=False,
                              skysub=True,
                              driz_cr=True,
                              driz_cr_corr=True,
                              driz_combine=False)


    # Remove unneeded files.
    unneeded_files = glob.glob('*med.fits') + glob.glob('*crclean.fits') \
        + glob.glob('*blt.fits') + glob.glob('*single_mask.fits') \
        + glob.glob('*wht.fits') + glob.glob('*sci.fits') \
        + glob.glob('*staticMask.fits') + glob.glob('*skymatch*')
    for unneeded_file in unneeded_files:
        os.remove(unneeded_file)


#-------------------------------------------------------------------------------# 

def sort_files(fits_files):
    """ Sorts FLT and FLCs into proposal directories in 'data/' and the 
    cosmic ray masks in 'outputs/crmasks'
    """

    crmask_files = glob.glob('*crmask.fits')

    path_to_crmasks = os.path.join(config.path_to_outputs, 'crmasks')

    # Sort the files 
    for fits_file in fits_files:
        proposid = str(fits.getval(fits_file,'proposid',ext=0))
        proposaldir = os.path.join(config.path_to_data, proposid)
        if not os.path.isdir(proposaldir):
            os.mkdir(proposaldir)
        shutil.move(fits_file, proposaldir)
        print("Moved {} to {}".format(fits_file, proposaldir))

    for crmask_file in crmask_files:
        shutil.move(crmask_file, path_to_crmasks)
        print("Moved {} to {}".format(crmask_file, path_to_crmasks))


#-------------------------------------------------------------------------------# 

def run_adriz_main():
    """ So AstroDrizzle doesn't explode, only feed it one visit at a time.
    """
    fits_files = glob.glob('*fl?.fits')
    flc_files = glob.glob('*flc.fits')

    file_dict = {}

    for flc_file in flc_files:
        proposid = fits.getval(flc_file,'proposid',ext=0)
        visit = fits.getval(flc_file,'rootname',ext=0)[4:6]
        if proposid not in file_dict:
            file_dict[proposid] = {}
            file_dict[proposid][visit] = []
        elif proposid in file_dict:
            if visit not in file_dict[proposid]:
                file_dict[proposid][visit] = []
    
        file_dict[proposid][visit].append(flc_file)


    for proposid in file_dict.keys():
        for visit in file_dict[proposid].keys():
            print("Running AstroDrizzle on following files: ")
            print(file_dict[proposid][visit])
            run_adriz(file_dict[proposid][visit])    


    sort_files(fits_files)


#-------------------------------------------------------------------------------# 
# The Main.
#-------------------------------------------------------------------------------# 

if __name__=='__main__':
    run_adriz_main()






