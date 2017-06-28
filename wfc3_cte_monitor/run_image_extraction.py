#! /usr/bin/env python

""" 
Module containing functions to perform the setup of the images in 
order to do analysis for the external CTE monitor.  This includes 
pulling out all an storing all relevant header information from the
images and performing source-finding, catalog-matching, and 
photometry on them. From here the FileInfo and the
Phot tables get filled. The end result is to have data from which
can measure CTE in the next step, `run_outputs.py`. 

Ideally this should only need to be run when new data is ingested.
Once you have the photometry you should be able to do analysis to 
your heart's content.

Generally called from main wrapper `uvis_external_cte.py`.
But can also be run from command line to do non-standard processing.

Authors:

    C.M. Gosmeyer, Feb. 2016

Use:
 
    Nominally execute from the main wrapper ``run_uvis_external_cte.py``.
    But can be run from the command line if user specifies desired 
    parameters in the Main.

    >>> python run_image_extraction.py

Notes:

    Translated from the IDL procedures `run_phot.pro` and
    `image_extr_org.pro`.

References:

    http://photutils.readthedocs.org/en/latest/photutils/aperture.html

    For spherematch function, see
    https://gist.github.com/eteq/4599814

To Do:
    1. Add arg parsing?

"""

from __future__ import print_function

import glob
import numpy as np
import os
import shutil
import time

from astropy.io import ascii
from astropy.io import fits
from drizzlepac import pixtosky
from drizzlepac import skytopix
from matplotlib.path import Path

# Why is the rum gone?
#from pyraf import iraf
#from iraf import noao, digiphot, daophot

import config
from set_paths_to_outputs import set_paths_to_outputs

from photutils_plus.apphot import apphot
from photutils_plus.plot_coords import plot_coords
from photutils_plus.phot_tools import read_in_photfile
from photutils_plus.phot_tools import meanclip_bkgrd

from analysis_tools.statistics.meanclip import meanclip

from database_interface import NGC104_FileInfo
from database_interface import NGC104_Phot
from database_interface import NGC6583_FileInfo
from database_interface import NGC6583_Phot
from database_interface import NGC6791_FileInfo
from database_interface import NGC6791_Phot

from database_update import update_fileinfo_table
from database_update import update_phot_table
from database_update import return_phot_dict_list
from database_update import return_fileinfo_dict

# Uncomment this line to get non-truncated printing of numpy arrays.
#np.set_printoptions(threshold=np.nan)

#-------------------------------------------------------------------------------# 

def apply_pam(imagename, param_dict):
    """
    Applies the Pixel Area Map correction to image.

    Parameters:
        imagename : string
            Name of the image.
        param_dict : dict
            Dictionary of parameters of the image.

    Returns:
        imdata_pamcorr : hdulist.data
            Data array of image with PAM correction applied.
            
    Outputs:
        nothing

    References:
        http://stackoverflow.com/questions/5830397/python-mmap-error-too-many-open-files-whats-wrong
    """

    chip = param_dict['chip']
    ext = return_extension(chip)

    # Get the chip number and read in the pixel area map.
    path_to_pams = os.path.join(config.path_to_data, 'PA_MAPS')
    if ext == 4:
        # chip 1, extension 4, dithered second exposure
        pam_file = os.path.join(path_to_pams, 'UVIS1wfc3_map.fits')
    elif ext == 1:
        #chip 2, extension 1, non-dithered first exposure
        pam_file = os.path.join(path_to_pams, 'UVIS2wfc3_map.fits')

    # Read in the extension.
    hdulist_orig = fits.open(imagename)
    imdata_orig = (hdulist_orig[ext].data) #.copy()
    hdulist_orig.close()

    # Read in the PAM correction.
    hdulist_pam = fits.open(pam_file)
    imdata_pam = (hdulist_pam[1].data) #.copy()
    hdulist_pam.close()

    # Apply the PAM correction
    imdata_pamcorr = imdata_orig * imdata_pam

    # Remove from memory.
    del imdata_orig
    del imdata_pam

    return imdata_pamcorr


#-------------------------------------------------------------------------------# 

def create_param_dict(imagename, flashlvl_desired, chargeinject='NONE', 
                       ngc104cal2=False, subdithers=False, xdithers=False):
    """
    Creates dictionary containing unique file information for images
    that will have their CTE quantified either by 
    [1] Being 'dithered' by a chip-length (~82 pixels) from chip2 to 
        chip1 in the y-direction. This yields a measure of CTE that 
        combines CTE degredation on both chips.
    [2] Being rotated 180-degrees. This yields a same-chip measure 
        of CTE.   

    Also serves to weed out undesired image types, necessary because 
    every program has special cases (see below).
    
    The dictionary will be referenced through the remainder of this 
    script and will be used to insert the image information into the 
    database.

    Parameters:
        imagename : string
            Name of the FLT or FLC file.
        flashlvl_desired : int or string
            Which post-flash level to process. Either
            0 for no post-flash. (Valid for ALL proposals and targets).
            6 for flashlvl 6 e-. (Valid only >= proposal 13083 and 
                target ngc104).
            12 for flashlvl 12 e-. (Valid only >= proposal 13083 and
                targets ngc104 and ngc6791).
            18 for flashlvl 18 e-. (Valid only >= proposal 14378 and
                target ngc104).
            24 for flashlvl 24 e-. (Same as for 18).
            33 for flashlvl 33 e- or flashexp 1.24. (Same as for 6).
            55 for flashlvl 55 e- or flashexp 1.59. (Same as for 6).
            91 for flashlvl 91 e- or flashexp 13. (Same as for 6).
            116 for flashlvl 116 e- or flashexp 21.5. (Same as for 6).
        chargeinject : string
            The level of charge injection. Changing this only relevant 
            for proposals 12348. Nominal images have chinject = 'None'.
        ngc104cal2 : {True, False}
            Switch on to process the non-nominal ngc104 field from 
            proposal 12692, targname 'NGC-104-CAL2'. False by default.
        subdithers : {True, False}
            Switch on to process the subdither exposures in proposals
            12348 and 12692. False by default.
        xdithers : {True, False}
            Switch on to process x-dithered images. False by default.

        The last four parameters are in place and set to False to speed 
        up processing. At this time there are no analysis tools set up
        for these special-case datasets. At present we are only concerned
        with nominal y-dither and 180-degree datasets of the nominal
        star fields. 

    Returns:
        param_dict : dict
            Dictionary containing image's header values and other
            useful information that will be used to fill in database
            and label plots and data files.

    Outputs:
        nothing

    Description of Programs (Everyone is special!): 
    11924 -
        3 visits of ngc6791, in F606W (30 and 360s) and F502N (60 and 420s), x AND y dithers.
    12348 -
        (charge injection test)
        4 visits of ngc104, in F502N (30 and 360s), y dithered with small 2.5 pixel dithers.
          Different charge injections used (none, continuous, 10-line, 17-line, and 25-line).
    12379 - 
        1 normal visit (01) of ngc6791, in F606W (30 and 360s) and F502N (60 and 420S)
        2 normal visits (02, 08) of ngc104, F606W (30 and 350s) and F502N (30 and 360s)
        1 ngc104 (07), F606W and F502N (348 s both), no dither
        1 ngc6791 (09) in F606W (30 and 348 s) and F502N (60 and 400s)
        2 ngc104 (11, 14) in F606W (30 and 348s) and F502N (348s), x AND y dithers.
    12692 -
        2 180-degree of ngc6583 (10, 11) in F606W (348s) and F502N (60 and 348s)
        3 ngc6791 (01,04,07) in F606W (30 and 360s) and F502N (60 and 420s), 
          where only long exptimes were y-dithered.
        3 ngc104 (02,05,08) in F606W (30 and 350s) and F502N (30 and 360s)
        3 ngc104-cal2 (03,06,09) in F606W (30 and 350s) and F502N (30 and 360s)
          This is a non-standard field of ngc104. 
        2 ngc104-cal2 (12,13) in F502N (348s) These observations are stepped in 
          smaller y-dithers than standard field.
    13083 -
        3 ngc6791 (01,03,05) in F502N (60 and 420s) and flashlvl 0 and 12.
        3 ngc104 (02,04,06) in F502N (348 s) and flashlvl 0,6,12 and flashexp 1.24,1.59,13,21.5
    13566 -
        2 ngc6791 (01,03) in F502N (60 and 420s) and flashlvl 0 and 12.
        2 ngc104 (02,04) in F502N (348 s) and flashlvl 0,6,12 and flashexp 1.24,1.59,13,21.5
    14012 -
        2 ngc6791 (01,03) in F502N (60 and 420s) and flashlvl 0 and 12.
        1 ngc6791 (05) A.Riess doing A.Riess things.
        2 ngc104 (02,04) in F502N (348 s) and flashlvl 0,6,12 and flashexp 1.24,1.59,13,21.5
    14378 -
        2 ngc6791 (01,03) in F502N (60 and 420s) and flashlvl 0 and 12.
        2 ngc104 (02,04) in F502N (348 s) and flashlvl 0,6,12,18,24 and flashexp 1.24,1.59,13,21.5 
          Added 30s in flashlvl 0 and 12. 

    Notes:
        I changed the filename scheme. In the IDL version, for example, 

        n6791_11924_55106_F606W_360_NOCI_F000_c0_1

        In this new Python version the same file is named as, 

        ngc6791_11924_55106_F606W_360_ciNO_pf000_cte0_2

        Notice that in the new naming scheme, the last extension
        refers to the chip. In the old scheme, it refered to the dither
        position.

        OLD: 1-> not dithered (postarg2=0), 2-> dithered (postar2=~82)
        NEW: 2-> chip2 observed (postarg2=0), 1-> chip1 observed (postarg2=~82)

        See the APT files of the proposals for visualization. Chip2 
        should always observe the field first and be the undithered 
        exposure.
    """

    # Fetch headers.
    priheader = fits.getheader(imagename, 0) 
    sciheader = fits.getheader(imagename, 1)

    # First check whether image was charge injected.
    # For now I will not process that dataset.
    chinject = priheader['CHINJECT']
    if chinject.upper() != chargeinject.upper():
        return 'false_chinject'

    # Second check that it is of the desired post-flash level.
    # First check FLASHSTA keyword because for some reason some
    # FLCs don't have the FLASHLVL keyword...
    flashsta = priheader['FLASHSTA']
    print(flashsta)
    if flashsta.lower() == 'successful':
        flashlvl = priheader['FLASHLVL']
        if int(flashlvl_desired) != int(flashlvl):
            return 'false_pf'
    elif int(flashlvl_desired) != 0:
        return 'false_pf'
    else: 
        flashlvl = 0
    print('flashlvl: ', flashlvl)

    targname = (priheader['TARGNAME']).split('NGC')[1]
    print('targname: ', targname)
    if targname[0] == '-':
        targname = targname[1:]

    # Third check whether using standard fields.
    # Note that ngc6583 was used in the 180-degree dataset.
    if targname == '104-CAL2' and ngc104cal2:
        # non-standard NGC104 field
        targname = '104CAL2'
    elif targname == '104-CAL2' and not ngc104cal2:
        return 'false_field'
        
    proposid = priheader['PROPOSID']

    expstart = priheader['EXPSTART']
    dateobs = str(expstart)[0:5]
    print("dateobs: ", dateobs)

    filt = priheader['FILTER']
    exptime = priheader['EXPTIME']
    naxis1 = sciheader['naxis1']
    naxis2 = sciheader['naxis2']
    try:
        mdrizsky = sciheader['mdrizsky']
    except:
        mdrizsky = 0
        print("No keyword 'MDRIZSKY' found in sci header.")

    # The postargs will tell whether the image was dithered.
    # For standard observations, want first image non-dithered and
    # second image to be y-dithered by a chip length.
    postarg1 = priheader['POSTARG1']
    postarg2 = priheader['POSTARG2']
    print("postarg1: ", postarg1)
    print("postarg2: ", postarg2)
    # In original script, postarg1 != 0 is always skipped.
    # Here I at least allow postarg1 = 0.1. All others will be skipped for 
    # now. I will retain the else cases in case we find a use for them.
    if float(postarg1) <= 0.1:
        # these 'xdither' variables don't actually do anything at this time
        xdither = 0 
    elif (float(postarg1) > 0.1 and not xdithers and targname != '6583'):
        # Add exception for target ngc6583 since that is the 180-degree
        # dataset and has unconventional dithering.
        xdither = 1
        return 'false_xdither'

    # Obtain keywords additional keywords that may be handy to have
    # in the database. 
    # (These will have nothing to do with matching images later.)
    flashdur = priheader['FLASHDUR']
    flashcur = priheader['FLASHCUR']
    shutrpos = priheader['SHUTRPOS']

    # Decide whether image is first in the pair or the second.
    if '6583' in targname:
        # The 180-degree dataset. Need hardcode which chips contain
        # the useable data...
        if any(code in imagename for code in \
            ['etq', 'euq', 'ewq', 'f2q',  'xsq', 'xtq', 'xvq', 'y1q']):
            chip = 1
        else:
            chip = 2

    else:
        # For all other datasets, the ydither tells whether the image was 
        # 'dithered' by a chip-length. Use this to find which chip is
        # the usable chip.
        if np.abs(postarg2) <= 0.1: # Not y-dithered. 
            # Nominally postarg2 = 0 or 0.1.
            chip = 2 #'y1'
        elif np.abs(postarg2) >= 80.0: # Was y-dithered. 
            # Nominally postarg2 = ~82.
            chip = 1 #'y2'
        else:
            if not subdithers:
                return 'false_subdither'
            else:
                # Small y-dithers from 12348, 12379, 12692.
                # Not sure what to do with them.
                # Retain this else in case find a use for these small dithers.
                print("Y-dither postarg2 {} cannot at this time be handled.".format(postarg2))
        
    
    if 'flc.fits' in imagename: # CTE-corrected
        cte = 1
    else: # Not CTE-corrected
        cte = 0

    # This is for sorting into directories later.
    if exptime <= 60.:
        length = 's'
    elif exptime > 60.:
        length = 'l'

    ext = return_extension(chip)

    # Add this in case want to match image's boundaries to master image.
    # I don't actually do this anymore, but doesn't hurt to keep storing
    # them...
    ra_lowerleft, dec_lowerleft = pixtosky.xy2rd('{}[{}]'.\
                                            format(imagename, ext),
                                            x=0, y=0,
                                            hms=False)
    ra_upperright, dec_upperright = pixtosky.xy2rd('{}[{}]'.\
                                            format(imagename, ext),
                                            x=4095, y=2051,
                                            hms=False)
    ra_lowerright, dec_lowerright = pixtosky.xy2rd('{}[{}]'.\
                                            format(imagename, ext),
                                            x=4095, y=0,
                                            hms=False)
    ra_upperleft, dec_upperleft = pixtosky.xy2rd('{}[{}]'.\
                                            format(imagename, ext),
                                            x=0, y=2051,
                                            hms=False)
    
    # Get the rootname (including the FLT or FLC extension)
    rootname = (imagename.split('/')[-1]).split('.fits')[0]

    param_dict = {'rootname' : rootname,
                   'targname' : 'ngc{}'.format(targname), 
                   'proposid' : proposid,
                   'dateobs' : int(dateobs), 
                   'filter' : filt.upper(), 
                   'exptime' : int(exptime),
                   'chinject' : chinject[0:2].upper(),
                   'flashlvl' : '{0:03}'.format(int(flashlvl)),
                   'ctecorr' : cte,
                   'chip' : chip,
                   'flashdur' : flashdur,
                   'flashcur' : flashcur,
                   'shutrpos' : shutrpos,
                   'postarg1' : postarg1,
                   'postarg2' : postarg2,
                   'naxis1' : naxis1,
                   'naxis2' : naxis2,
                   'length' : length,
                   'ra_lowerleft' : ra_lowerleft,
                   'dec_lowerleft' : dec_lowerleft,
                   'ra_lowerright' : ra_lowerright,
                   'dec_lowerright' : dec_lowerright,                   
                   'ra_upperright' : ra_upperright,
                   'dec_upperright' : dec_upperright,
                   'ra_upperleft' : ra_upperleft,
                   'dec_upperleft' : dec_upperleft,
                   'mnclip_bkgrd' : mdrizsky}


    return param_dict


#-------------------------------------------------------------------------------# 

def do_photom(imagename, param_dict, coofile, 
              aprads=[2,3,5,7,10,12,15,18,20,24,28,32,36,40], 
              overwrite=False, outname='default', use_iraf=False, 
              mask_bad_pxls=True, mask_crs=True):
    """
    Applies the pixel area map and runs photometry.

    Parameters:
        imagename : string
            Name of the FLT or FLC file.
        param_dict : dict
            Dictionary containing image's header values and other
            useful information that will be used to fill in database
            and label plots and data files.
        coofile : string
            Name of the coordinate file.
        aprads : list of ints/floats
            The aperture radii (pixels) on which to perform photometry.
        overwrite : {True, False}
            Switch on to overwrite pre-existing *mag file. False by 
            default.
        outname : string
            Change the default output name of the *mag file from
            "<rootname>.mag".
        use_iraf : {True, False}
            Run photometry with ``iraf.phot`` instead of ``phututils.
            aperture_photometry`` and insert into 
            'uvis_external_cte_iraf.db'.
        mask_bad_pxls : {True, False}
            Masks all non-zero pixels in the DQ array *if* doing photometry
            with Python (ie, 'use_iraf' must be False for this option to 
            take affect).

    Returns:
        magfile : Astropy Table
            Table containing photometry information.    
 
    Outputs:
        "<rootname>.mag" or "<'outname'>.mag". Contains same information 
        as the returned Table.
    """
    rootname = param_dict['rootname']

    if overwrite:
        if os.path.isfile('{}.mag'.format(rootname)):
            os.remove('{}.mag'.format(rootname))

    chip = param_dict['chip']
    ext = return_extension(chip)
    print("CHIP: {}, EXT: {}".format(chip, ext))

    # Apply the PAM correction.
    imdata_pamcorr = apply_pam(imagename, param_dict)

    if use_iraf:
        # Make the aperture list IRAF-compatible.
        apertures = ''
        for ap in aprads:
            apertures += str(ap) + ','
        apertures[:-1]

        # Make an IRAF-compatible coo file locally.
        coofile_iraf = make_coo_iraf_friendly(coofile)
        outname = '{}.mag'.format(rootname)

        # Write a PAM-corrected image locally.
        image_pamcorr = '{}_iraf.fits'.format(rootname)
        hdu = fits.PrimaryHDU(imdata_pamcorr)
        hdu.writeto(image_pamcorr, clobber=True)

        iraf.phot.unlearn()         
        iraf.phot(image='{}_iraf.fits[{}]'.format(rootname, 0), 
                  interactive='no', 
                    verify='no', 
                    coords=coofile_iraf, 
                    output=outname,
                    calgorithm='centroid', 
                    annulus=12.0,
                    dannulus=10.0,
                    apertures=apertures)
        # Remove iraf coofile.
        os.remove(coofile_iraf)
        # Remove PAM-corrected image.
        os.remove(image_pamcorr)

        magfile = '{}.mag'.format(rootname)

    elif not use_iraf:
        # Check whether to mask bad pixels, cosmic rays, both, or none.
        if mask_bad_pxls and mask_crs:
            print("Masking both bad pixels and cosmic rays.")
            mask_cr = mask_cosmicrays(rootname, chip)
            mask_badpxl = mask_bad_pixels(imagename, ext)
            mask = mask_cr + mask_badpxl
        elif mask_crs:
            print("Masking cosmic rays only.")
            mask = mask_cosmicrays(rootname, chip)
        elif mask_bad_pxls:
            print("Masking bad pixels only.")
            mask = mask_bad_pixels(imagename, ext)
        else:
            print("Masking nothing")
            mask = None

        # Finally run photometry.
        # Centroiding still does not work in photutils for multiple sources? 
        phot_tab = apphot(imdata_pamcorr, 
               ext=ext, 
               coofile=coofile, 
               colnames=['extr_xpix','extr_ypix'],
               outname=outname,
               ap_radii=aprads,
               sep_out=False, 
               backmethod='mean',
               backglobal=False, 
               method='exact',
               pixelwise_error=None, 
               effective_gain=None, 
               mask=mask,
               annulus=10., 
               dannulus=10.)

        # Clear variable
        del imdata_pamcorr

        # Since default, this should be the magfile name.
        if len(phot_tab) == 0:
            magfile = phot_tab
        else:
            magfile = '{}.mag'.format(rootname)

    return magfile


#-------------------------------------------------------------------------------# 

def extract_sources_from_master(imagename, param_dict, master_cat, out_dir=''):
    """
    Finds which sources in the master catalog are present on the desired chip
    of the image. Returns a coordinate file of the sources.

    Parameters:
        imagename : string
            Name of the FLT or FLC file.
        param_dict : dict
            Dictionary containing image's header values and other
            useful information that will be used to fill in database
            and label plots and data files.
        master_cat : string
            Name and location of the master catalog.
        out_dir : string
            Path to desired output location.

    Returns:
        coofile : string
            Name of the coordinate file.
        diagplot : string
            Name of the diagnostic JPG image, which plots all the 
            identified master sources onto the image.

    Outputs:
        "<rootname>.coo". The coordinate file.
        "<rootname>.jpg". Image with identified sources overplotted.
    """

    # read in the master catalog
    master_data = ascii.read(master_cat)
    master_id = master_data['col1'] #id
    master_ra = master_data['col4'] #ra
    master_dec = master_data['col5'] #dec

    chip = int(param_dict['chip'])
    ext = return_extension(chip)

    ra_lowerleft_pix = param_dict['ra_lowerleft']
    dec_lowerleft_pix = param_dict['dec_lowerleft']
    ra_lowerright_pix = param_dict['ra_lowerright']
    dec_lowerright_pix = param_dict['dec_lowerright']

    ra_upperright_pix = param_dict['ra_upperright']
    dec_upperright_pix = param_dict['dec_upperright']
    ra_upperleft_pix = param_dict['ra_upperleft']
    dec_upperleft_pix = param_dict['dec_upperleft']

    # Define corners of polygon. (chip in sky coords)
    verts = [(ra_lowerleft_pix, dec_lowerleft_pix),
             (ra_lowerright_pix, dec_lowerright_pix),
             (ra_upperright_pix, dec_upperright_pix),
             (ra_upperleft_pix, dec_upperleft_pix)]   

    # Create polygon
    chip_reg = Path(verts)  

    indices_in_chip = chip_reg.contains_points(list(zip(master_ra, master_dec)))
    id_retain_list = master_id[indices_in_chip]
    ra_retain_list = master_ra[indices_in_chip]   
    dec_retain_list = master_dec[indices_in_chip]

    # Transform the RA and DECs into the chip's pixel coords.
    x_retain_list = []
    y_retain_list = []
    for ra, dec in zip(ra_retain_list, dec_retain_list):
        # This part is very, very slow. Parallelize? 
        xpix, ypix = skytopix.rd2xy('{}[{}]'.\
                                    format(imagename, ext),
                                    ra, dec)
        x_retain_list.append(xpix)
        y_retain_list.append(ypix)

    # Create ID list of all found sources. 
    id_found_list = np.arange(1, len(id_retain_list)+1)

    # Write to a *coo file.
    rootname = param_dict['rootname']
    if out_dir == '':
        coofile = '{}.coo'.format(rootname)
    else:
        coofile = os.path.join(out_dir, '{}.coo'.format(rootname))

    coofile_open = open(coofile, 'w')
    coofile_open.write('#master_id\textr_id\textr_xpix\textr_ypix\textr_ra\textr_dec\n')
    for j in range(len(id_retain_list)):
        row_str = ''
        row_str = '{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                                    id_retain_list[j], 
                                    id_found_list[j],
                                    x_retain_list[j],
                                    y_retain_list[j],
                                    ra_retain_list[j],
                                    dec_retain_list[j])
        coofile_open.write(row_str)
    coofile_open.close()
    
    # Plot the coordinates over the desired chip. 
    diagplot = plot_coords(imagename=imagename, 
                           coofile=coofile, 
                           ext=ext,
                           xcol='extr_xpix', 
                           ycol='extr_ypix',
                           outname=rootname)

    return coofile, diagplot


#-------------------------------------------------------------------------------# 

def insert_rows_to_db(table_dict, table_type, targname):
    """
    Decides which cluster's table and which table type (FileInfo or Phot) 
    is to be updated.

    Parameters:
        table_dict : dict
            Dictionary of values, where keys are the columns, to be 
            inserted into table.
        table_type : string
            'phot' or 'fileinfo'
        targname : string
            Name of the target cluster.

    Returns: 
        nothing

    Outputs:
        updated Phot or FileInfo table
    """
    table_type = table_type.lower()
    targname = targname.lower()

    if table_type == 'phot' and targname == 'ngc104':
        Table = NGC104_Phot
    elif table_type == 'phot' and targname == 'ngc6791':
        Table = NGC6791_Phot
    elif table_type == 'phot' and targname == 'ngc6583':
        Table = NGC6583_Phot
    elif table_type == 'fileinfo' and targname == 'ngc104':
        Table = NGC104_FileInfo
    elif table_type == 'fileinfo' and targname == 'ngc6791':
        Table = NGC6791_FileInfo
    elif table_type == 'fileinfo' and targname == 'ngc6583':
        Table = NGC6583_FileInfo

    if table_type == 'fileinfo':
        update_fileinfo_table(Table, table_dict)

    if table_type == 'phot':
        update_phot_table(Table, table_dict)


#-------------------------------------------------------------------------------# 

def make_coo_iraf_friendly(coofile):
    """Make coo files that IRAF will accept.
    photutils is 0-based
    IRAF is 1-based
    hooray!

    Parameters:
        coofile : string
            Name of the Python coordinate file.

    Returns:
        outname : string
            Name of the IRAF coordinate file.

    Outputs:
        nothing
    """
    # Read in x,y coords
    coo_tab = ascii.read(coofile)
    # Need offset by 1 because IRAF 1-based.
    xcoords = np.array(coo_tab['extr_xpix'])-1.0
    ycoords = np.array(coo_tab['extr_ypix'])-1.0
    
    # Make outname.
    outname = (coofile.split('/')[-1]) + '_iraf'

    # Write out to new coo file.
    ascii.write([xcoords, ycoords], outname)

    return outname


#-------------------------------------------------------------------------------# 

def mask_bad_pixels(imagename, ext):
    """ Masks out all pixels not equal to 0 (nominal) or 32 (CTE trail) in 
    the DQ array.
    
    Parameters:
        imagename : string
            Name of the FLT or FLC file.
        ext : int
            The extension of the science data.

    Returns:
        mask : array of bool
            The masked image. 
    """
    # Read in the DQ array.
    hdulist = fits.open(imagename)
    dq = (hdulist[ext+2].data) 
    hdulist.close()

    # Mask out (change to True) all values not equal to 0 or 32.
    mask = np.zeros_like(dq, dtype=bool)
    for i in range(dq.shape[0]):
        for j in range(dq.shape[1]):
            if dq[i][j] != 0 and dq[i][j] != 32:
                mask[i][j] = True

    return mask


#-------------------------------------------------------------------------------# 

def mask_cosmicrays(rootname, chip):
    """ Masks out all cosmic ray-flagged pixels (equal to 0). 
    The *crmask files were generated with AstroDrizzle. 

    The FLT and FLC versions of the image share the same mask.
    
    Parameters:
        rootname : string
            Rootname of the image.
        chip : int
            The number of the chip.

    Returns:
        mask : array of bool
            The masked image. 
    """
    if int(chip) == 1:
        ext_num = 2
    elif int(chip) == 2:
        ext_num = 1

    crmask = os.path.join(config.path_to_outputs, 'crmasks', 
        '{}_sci{}_crmask.fits'.format(rootname.split('_')[0], ext_num))

    # Read in the mask.
    hdulist = fits.open(crmask)
    data = (hdulist[0].data) 
    hdulist.close()

    # Mask out (change to True) all values equal to 0.
    mask = np.zeros_like(data, dtype=bool)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            if data[i][j] == 0:
                mask[i][j] = True

    return mask


#-------------------------------------------------------------------------------# 

def return_extension(chip):
    """ 
    Returns UVIS extension corresponding to chip number.

    Parameters:
        chip : int/string
            (1,2) Chip number.
 
    Returns:
        ext : int
            Extension corrsponding to chip number.
    """

    if int(chip) == 1:
        ext = 4
    elif int(chip) == 2:
        ext = 1

    return ext

#-------------------------------------------------------------------------------# 
#-------------------------------------------------------------------------------# 

def run_image_extraction(flashlvl, programset, ctecorr, visits=[], 
    overwrite=False, insert_into_db=True, use_iraf=False, do_find=True):
    """Loops over given permutations of program number, flashlvl, and CTE 
    correction.  Generates out-filename, finds sources, and performs
    photometry.

    This step should only need be done when data first ingested, or if
    already-existing data was re-calibrated. 

    Everything generated here is needed in the analysis that follows.

    Parameters:
        flashlvl : int or string
            Which post-flash level to process. Either
            0 for no post-flash. (Valid for ALL proposals and targets).
            6 for flashlvl 6 e-. (Valid only >= proposal 13083 and 
                target ngc104).
            12 for flashlvl 12 e-. (Valid only >= proposal 13083 and
                targets ngc104 and ngc6791).
            18 for flashlvl 18 e-. (Valid only >= proposal 14378 and
                target ngc104).
            24 for flashlvl 24 e-. (Same as for 18).
            33 for flashlvl 33 e- or flashexp 1.24. (Same as for 6).
            55 for flashlvl 55 e- or flashexp 1.59. (Same as for 6).
            91 for flashlvl 91 e- or flashexp 13. (Same as for 6).
            116 for flashlvl 116 e- or flashexp 21.5. (Same as for 6).
        programset : int or string
            Which group of proposals to process. Either
            'all' - All that are available for the post-flash group.
            'last' - Newest of post-flash group available).
            <proposal number> - Specific proposal ID. Must be valid 
                for post-flash group, else will return an error.
            '180' - The 180-degree dataset from proposal 12692, visits
                10 and 11.
        ctecorr : {True, False}
            True for FLCs (pixel-based CTE-corrected)
            False for FLTs (not pixel-based CTE-corrected)
        visits : list of ints/strings
            Visit numbers to process within selected proposal. 
        overwrite : {True, False}
            Turn on to overwrite existing extracted images and *coo
            and *mag files.  *need to implement*
        insert_into_db : {True, False}
            Insert image, mag, and coo file into into FileInfo and 
            Phot tables. True by default.  
        use_iraf : {True, False}
            Run photometry with ``iraf.phot`` instead of ``phututils.
            aperture_photometry`` and insert into 
            'uvis_external_cte_iraf.db'.
        do_find : {True, False}
            Creates a new *coo file with all coordinates transformed to 
            match the Master Catalog's.  Time consuming.  If just re-running
            photometry and the coo files already exist, can shut this off.
            By default, True.

    Returns:
        nothing

    Outputs: 
        See documentation for functions called within this. Nominally,
        * "<rootname>.coo". The coordinate file.
        * "<rootname>.png". Image with identified sources overplotted.
        * "<rootname>.mag". The photometry file.
        * Updated FileInfo and Phot tables for all target clusters
            observed within range of input parameters.

    Notes:
        Analog of IDL `run_phot.pro`. 

        The pipeline allows many non-nominal data subsets to be processed, 
        but even so, in this function several datasets are thrown out
        for various reasons, since I do not foresee them ever being useful
        for the monitor. These are
        * Visit 07 of 12379 - Undithered.
        * Visit 13 of 12379 - Data quality issues, retaken in Visit 14.
        * Visit 05 of 14012 - A.Riess's added-on binning test.

        In addition, images with charge injections of 'CONT' from 12348 
        are not run through TWEAKREG (`run_twealreg.py`) since it had 
        trouble identifying sources. Proceed with caution.

    """
    # Check types.
    flashlvl = str(flashlvl)
    programset = str(programset)

    # Import globals.
    program_list = config.program_list
    path_to_data = config.path_to_data
    path_to_mastercat = os.path.join(path_to_data, 'master_cat')

    # First cut down the list to >= 13083 data if postflashed-only data selected.
    if flashlvl != '0':
        program_list = program_list[4:]

    # Second create a list of all the data directories required.
    if programset == 'all':
        datadirs = [os.path.join(path_to_data, program) for \
                    program in program_list]
    elif programset == 'last':
        datadirs = [os.path.join(path_to_data, program_list[-1])]
    elif programset == '180':
        datadirs = [os.path.join(path_to_data, '12692')]

    else:
        if programset not in program_list:
            print("Invalid proposal ID of {} for a post-flash set. (Must be >= 13083).".format(programset))
            print("Valid propals are {}".format(program_list))
            print("If you do indeed want a program set < 13083, set 'flashlvl=0'.")
            return "Invalid post-flash level for proposid."
        else:
            datadirs = [os.path.join(path_to_data, programset)]

    print(datadirs)

    # Third loop over the data directories, globbing for files in each.
    for datadir in datadirs:
        print(datadir)
        # Choose either FLTs or FLCs.
        if ctecorr:
            print("Processing FLCs.")
            image_list = glob.glob(os.path.join(datadir, '*flc.fits'))
        else:
            print("Processing FLTs.")
            image_list = glob.glob(os.path.join(datadir, '*flt.fits'))

        # Check that for proposal 14012 don't process Visit 5 (A.Riess's 
        # thing) or for proposal 12379 don't process Visit 13 (guide star 
        # failure) and Visit 07 (undithered).
        if '14012' in datadir:
            clean_image_list = []
            for image in image_list:  
                if '05' not in image:
                    clean_image_list.append(image)
            image_list = clean_image_list
        if '12379' in datadir:
            clean_image_list = []
            for image in image_list:  
                if '13' and '07' not in image:
                    clean_image_list.append(image)
            image_list = clean_image_list      
        # If programset is 180, then you are only interested in processing
        # visits 10 and 11 from program 12692.      
        if programset == '180':
            clean_image_list = []
            for image in image_list:
                if ('10' in image) or ('11' in image):
                    clean_image_list.append(image)
            image_list = clean_image_list
        # If specified only certain visits be processed, sift them out.
        if visits != []:
            clean_image_list = []
            for image in image_list:  
                for visit in visits:
                    visit = str(visit)
                    if visit in image:
                        clean_image_list.append(image)
            image_list = clean_image_list 

        # Loop over the FITS files.
        i=1
        print(image_list)
        for imagename in image_list:

            # Create the name of the processed out file.
            param_dict = create_param_dict(imagename, 
                                      flashlvl_desired=flashlvl,
                                      chargeinject='None', 
                                      ngc104cal2=False,
                                      subdithers=False, 
                                      xdithers=False)
            print(param_dict)
            if 'false' in param_dict:
                parsed_flag = param_dict
                print('{} is being skipped because of flag "{}"'.format(
                                                    imagename, 
                                                    parsed_flag))

            else:
                print(i)
                i +=1
                # Get coordinates of sources in master catalog that match to chip image.
                # Insert into a *coo file to be named ipppssoot_fl?.coo and plot the sources over 
                # the desired chip of image.
                path_to_workingdir = os.getcwd()
                out_dir = os.getcwd()
                targname = (param_dict['targname']).lower()
                master_cat = os.path.join(os.path.join(path_to_mastercat, 'final'), 
                    '{}_master.cat'.format(targname))

                # Create paths to the photometry files and plots.            
                print("path_to_workingdir: {}".format(path_to_workingdir))
                path_to_outputs = config.path_to_outputs
                print("path_to_outputs: {}".format(path_to_outputs))

                if use_iraf:
                    basedir = 'photom_iraf'
                else:
                    basedir = 'photom'

                path_to_photom = set_paths_to_outputs(path_to_outputs=path_to_outputs, 
                                                      basedir=basedir, 
                                                      flashlvl=flashlvl, 
                                                      ctecorr=ctecorr, 
                                                      mostrecent=False)
                print("path_to_photom: {}".format(path_to_photom))
                photfile_dir = '{}_{}_{}'.format(param_dict['targname'], 
                                                 param_dict['filter'], 
                                                 param_dict['length']).lower()
                print("photfile_dir: {}".format(photfile_dir))
                path_to_photdir = os.path.join(path_to_photom, photfile_dir)
                print("path_to_photdir: {}".format(path_to_photdir))
                path_to_diag = os.path.join(path_to_photom, 'diagnostics')
                print("path_to_diag: {}".format(path_to_diag))


                # Do the source-finding, ie, get the coordinates from the Master image.
                if do_find:
                    print("Extracting coordinates from Master...")
                    coofile, diagplot = extract_sources_from_master(imagename, param_dict, 
                        master_cat, out_dir=out_dir)
                elif not do_find:
                    print("NOT extracting coordinates, presumably a COO file already exists...") 
                    coofile = os.path.join(path_to_photdir, param_dict['rootname'] + '.coo')    

                print("coofile: {}".format(coofile)) 


                # Apply PAM and perform photometry.
                # Put outputs in a *mag file.  
                print("Doing photometry...")
                magfile = do_photom(imagename, param_dict, coofile, 
                    overwrite=overwrite, outname=param_dict['rootname'],
                    use_iraf=use_iraf)

                print("magfile: {}".format(magfile))


                # Return dictionaries to be put into database tables.
                if len(magfile) != 0:

                    # Check if want to insert into database.
                    if insert_into_db: 
                        fileinfo_dict = return_fileinfo_dict(
                            imagename=imagename, 
                            coofile=coofile.split('/')[-1], 
                            magfile=magfile, path_to_photom=path_to_photom, 
                            param_dict=param_dict)
                        phot_dict_list = return_phot_dict_list(
                            imagename=imagename, 
                            coofile=coofile, magfile=magfile, 
                            param_dict=param_dict,
                            use_iraf=use_iraf)

                        targname = param_dict['targname']

                        # Insert the image's entry into the FileInfo table.
                        insert_rows_to_db(fileinfo_dict, 'fileinfo', targname)
                        # Insert each source's entry into Phot table.
                        for phot_dict in phot_dict_list:
                            insert_rows_to_db(phot_dict, 'phot', targname)

                        # Clear variables.
                        del fileinfo_dict
                        del phot_dict_list

                    # Now that the photometry and file information is recorded,
                    # ready to move coo, mag, and diagnostic plots.
                    # Create directories that don't exist.

                    if not os.path.isdir(path_to_photdir):
                        os.mkdir(path_to_photdir)
                    if not os.path.isdir(path_to_diag):
                        os.mkdir(path_to_diag)

                    # Move and overwrite files to the photfile directory.
                    shutil.move(os.path.join(path_to_workingdir, magfile), 
                                os.path.join(path_to_photdir, magfile))
                    if do_find:
                        coofile = coofile.split('/')[-1]
                        shutil.move(os.path.join(coofile), 
                                    os.path.join(path_to_photdir, coofile))

                        shutil.move(os.path.join(path_to_workingdir, diagplot), 
                                    os.path.join(path_to_diag, diagplot))


                else:
                    print("The coordinate file {} and diag plot {} were not \
                          moved because no sources were found in the mag file. \
                          ".format(coofile, diagplot))
                    open_file = open('{}_error.txt'.format(imagename.split('/')[-1]), 'w')
                    open_file.write("The coordinate file {} and diag plot {} were not \
                          moved because no sources were found in the mag file. \
                          ".format(coofile, diagplot))
                    open_file.close()


#-------------------------------------------------------------------------------# 
# The Main. Use for testing.
#-------------------------------------------------------------------------------# 


if __name__=='__main__':
    """
    flashlvl=0
    programset='180' #'last'
    ctecorr=True
    run_image_extraction(flashlvl, programset, ctecorr, overwrite=True, 
        insert_into_db=True)
    """
    imagename = 'ibwb02eiq_flc.fits'
    flashlvl_desired = 0
    coofile = 'ibwb02eiq_flc.coo'


    param_dict = create_param_dict(imagename, flashlvl_desired, chargeinject='NONE', 
                       ngc104cal2=False, subdithers=False, xdithers=False)
    magfile, param_dict = do_photom(imagename, param_dict, coofile, 
                    overwrite=False, outname=coofile.split('.coo')[0],
                    use_iraf=False)
    print(param_dict)
    #phot_dict_list = return_phot_dict_list(
    #                        imagename=imagename, 
    #                        coofile=coofile, magfile=magfile, 
    #                        param_dict=param_dict,
    #                        use_iraf=True)
    #print(phot_dict_list)
    #targname = param_dict['targname']
    #for phot_dict in phot_dict_list:
    #    insert_rows_to_db(phot_dict, 'phot', targname)




    
