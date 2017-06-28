"""Module for re-populating the external CTE monitor database's
FileInfo, Master, and Phot tables, without having to redo the 
photometry calculations.  To fill out the FileInfo and Phot tables, it 
just reads the existing coo and mag files in the 'outputs' directory
and gets additional information from the headers of the FITS files in
the 'data' directory. To fill out the Master table, it reads in the 
master catalogs from the 'data/master_cat' directory.

To re-populate the Results tables with the CTE measurements you still 
need to run pipeline with parameters:

>>> python run_uvis_external_cte.py --pr 'all' --dofind 'n' --dophot 'n' --pl 'ratio_ypos'

Authors:

    C.M. Gosmeyer, May 2016 

Use:

    Can enter re-population functions into the Main for command-line
    execution. Then

    >>> python database_reset.py

To do:
    Create a function to re-populate the Results table?
    The Results table can be reset by looping over the outputs/finalresults/most_recent
    directory?

"""

from __future__ import print_function

import glob
import os
import shutil
from astropy.io import ascii
from config import path_to_data
from config import path_to_outputs
from config import target_list

from photutils_plus.phot_tools import meanclip_bkgrd

from database_query import return_tables
from database_update import update_fileinfo_table
from database_update import update_phot_table
from database_update import update_master_table
from database_update import return_fileinfo_dict
from database_update import return_phot_dict_list

from run_image_extraction import apply_pam
from run_image_extraction import create_param_dict
from run_image_extraction import do_photom

from database_interface import NGC104_FileInfo
from database_interface import NGC104_Master
from database_interface import NGC104_Phot
from database_interface import NGC104_Results
from database_interface import NGC6583_FileInfo
from database_interface import NGC6583_Master
from database_interface import NGC6583_Phot
from database_interface import NGC6583_Results
from database_interface import NGC6791_FileInfo
from database_interface import NGC6791_Master
from database_interface import NGC6791_Phot
from database_interface import NGC6791_Results


#-------------------------------------------------------------------------------# 
# Helper functions.
#-------------------------------------------------------------------------------# 

def get_proposid(imagename):
    """ 
    Fetches proposal ID of given imagename based on its ipppssoot code.
    Unfortunately, can't think of a better way to do this...

    Parameters:
        imagename : string
            Name of the image.

    Returns:
        proposid : string
            Proposal number.

    Outputs:
        nothing
    """
    if 'ibc6' in imagename:
        proposid = '11924'
    elif 'iblk' in imagename:
        proposid = '12348'
    elif 'ibnh' in imagename:
        proposid = '12379'
    elif 'ibwb' in imagename:
        proposid = '12692'
    elif 'ic5w' in imagename:
        proposid = '13083'
    elif 'ichh' in imagename:
        proposid = '13566'
    elif 'icqu' in imagename:
        proposid = '14012'
    elif 'id1v' in imagename:
        proposid = '14378'

    return proposid


#-------------------------------------------------------------------------------# 
# Functions to populate the FileInfo and Phot tables.
#-------------------------------------------------------------------------------# 

def populate_fileinfo_table(targname):
    """
    Populates FileInfo table of 'targname'. Output to be directed to
    :func:`populate_phot_table`.
    Assumes all data are in directories and just want to repopulate
    the database from scratch.

    Parameters:
        targname : string
            Name of the target cluster.

    Returns:
        param_dicts : list of dictionaries
            The 'param_dict' for each image.
        paths_to_photom : list of strings
            The path to photometry files for each image.
        magfiles : list of strings
            The name of the *mag file for each image.
        imagenames : list of strings
            Names of the images.

    Outputs:
        Repopulated FileInfo table.
    """
    print("targname: {}".format(targname))

    FileInfo, Master, Phot, Results = return_tables(targname)
    basepath_to_photom = os.path.join(path_to_outputs, 'photom')
    print("basepath_to_photom: {}".format(basepath_to_photom))
    pf_dirs = glob.glob(os.path.join(basepath_to_photom, 'pf*'))
    print("pf_dirs: {}".format(pf_dirs))

    param_dicts = []
    paths_to_photom = []
    magfiles = []
    imagenames = []

    for pf_dir in pf_dirs:
        postflash = int((pf_dir.split('pf')[1]).split('_')[0])
        print("Postflash: {}".format(postflash))

        phot_dirs = []
        phot_dirs = glob.glob(os.path.join(pf_dir, '{}*'.format(targname.lower())))
        print("phot_dirs: {}".format(phot_dirs))
        for phot_dir in phot_dirs:
            existing_magfiles = glob.glob(os.path.join(phot_dir, '*.mag'))
            for existing_magfile in existing_magfiles:
                path_to_photom = '/'.join((existing_magfile.split('/'))[0:-1])
                print("path_to_photom: {}".format(path_to_photom))
                magfile = (existing_magfile.split('/')[-1])
                print("magfile: {}".format(magfile))
                imagename = '{}.fits'.format(magfile.split('.mag')[0])
                print("imagename: {}".format(imagename))
                coofile = '{}.coo'.format(magfile.split('.mag')[0])
                print("coofile: {}".format(coofile))

                proposid = get_proposid(imagename)
                print("proposid: {}".format(proposid))
                imagepath = os.path.join(path_to_data, proposid)
                print("imagepath: {}".format(imagepath))
                imagename = os.path.join(imagepath, imagename)
                print("imagename: {}".format(imagename))

                param_dict = create_param_dict(imagename, postflash, 
                    chargeinject='NONE', ngc104cal2=False, subdithers=False, 
                    xdithers=False)

                # Recalculate the background.
                # First apply the PAM correction.
                imdata_pamcorr = apply_pam(imagename, param_dict)
                # Read in x,y coords of sources in order to make a mask.
                coo_tab = ascii.read(os.path.join(path_to_photom, coofile))
                xcoords = list(coo_tab['extr_xpix'])
                ycoords = list(coo_tab['extr_ypix'])
                # Add the average background to the param_dict

                bkgrd, bkgrd_sigma = meanclip_bkgrd(imdata_pamcorr, 
                    backmethod='mean', 
                    xcoords=xcoords, ycoords=ycoords,
                    naxis1=param_dict['naxis1'], naxis2=param_dict['naxis2'], 
                    detector='UVIS', maskrad='', multiple_sources=True) 
                print("Background: {} e-".format(bkgrd))
                param_dict['mnclip_bkgrd'] = bkgrd

                # Create the dictionary that will be inserted into db.
                fileinfo_dict = return_fileinfo_dict(imagename, coofile, 
                    magfile, path_to_photom, param_dict)
                # Update the FileInfo table in db.
                update_fileinfo_table(FileInfo, fileinfo_dict)

                param_dicts.append(param_dict)
                paths_to_photom.append(path_to_photom)
                magfiles.append(magfile)
                imagenames.append(imagename)
                print(" ")

    return param_dicts, paths_to_photom, magfiles, imagenames


#-------------------------------------------------------------------------------# 

def populate_phot_table(targname, param_dicts=[], paths_to_photom=[], 
                        magfiles=[], imagenames=[], redo_phot=False):
    """ 
    Populates Phot table of 'targname'. Uses outputs from 
    :func:`populate_fileinfo_table`.
    Assumes all data are in directories and just want to repopulate
    the database from scratch.

    Parameters:
        targname : string
            Name of the target cluster.
        param_dicts : list of dictionaries
            The 'param_dict' for each image.
        paths_to_photom : list of strings
            The path to photometry files for each image.
        magfiles : list of strings
            The name of the *mag file for each image.
        imagenames : list of strings
            Names of the images.
        redo_phot : {True, False}
            False by default. Set to True to recalculate all photometry.

    Returns:
        nothing

    Outputs:
        Repopulated Phot table.
    """
    FileInfo, Master, Phot, Results = return_tables(targname)

    for param_dict, path_to_photom, magfile, imagename in \
        zip(param_dicts, paths_to_photom, magfiles, imagenames):
         
        if redo_phot:
            print("Redoing photometry on {}!".format(imagename))
            coofile = os.path.join(path_to_photom, '{}.coo'.format(magfile.split('.mag')[0]))
            magfile, param_dict = do_photom(imagename=imagename, 
                param_dict=param_dict, coofile=coofile, 
                overwrite=True,
                outname=(coofile.split('/')[-1]).split('.coo')[0], 
                calc_global_bkgrd=False)

            phot_dict_list = return_phot_dict_list(imagename, coofile, magfile, param_dict)
            for phot_dict in phot_dict_list:
                update_phot_table(Phot, phot_dict)

            # Move magfile to replace the old.
            shutil.move(magfile, os.path.join(path_to_photom, magfile))

        else: 
            coofile = os.path.join(path_to_photom, '{}.coo'.format(magfile.split('.mag')[0]))
            magfile = os.path.join(path_to_photom, magfile)
            
            phot_dict_list = return_phot_dict_list(imagename, coofile, magfile, param_dict)
            for phot_dict in phot_dict_list:
                update_phot_table(Phot, phot_dict)


#-------------------------------------------------------------------------------# 

def populate_all_fileinfo_phot_tables():
    """ 
    Populates FileInfo and Phot tables of all target clusters.
    Assumes the two table types will always be populated together.

    Parameters:
        nothing

    Returns:
        nothing

    Outputs:
        Repopulated FileInfo and Phot tables for every target cluster.
    """

    targnames = target_list

    for targname in target_list:
        print("***Populating the FileInfo table for {}.".format(targname))
        param_dicts, paths_to_photom, magfiles, imagenames = \
            populate_fileinfo_table(targname)
        print("***Populating the Phot table for {}.".format(targname))
        populate_phot_table(targname, param_dicts, paths_to_photom, 
            magfiles, imagenames, redo_phot=True)


#-------------------------------------------------------------------------------# 
# Functions to populate the Master tables.
#-------------------------------------------------------------------------------# 

def populate_master_table(master_table, master_catalog):
    """ 
    Populates the given Master table.
    In theory, this only needs be done once, except in the case of a database
    reset. Master tables will not be updated after initial population.

    Parameters:
        master_table : Table object
            The table in which to insert information in 'master_catalog'.
        master_catalog : string
            Name of the text file containing master catalog.

    Returns:
        nothing

    Outputs:
        Populated Master table.
    """
    # Read in the master catalog
    data = ascii.read(master_catalog)
    master_id = list(data['col1'])
    xpix = list(data['col2'])
    ypix = list(data['col3'])
    ra = list(data['col4'])
    dec = list(data['col5'])

    master_dict_list = []
    # Put each row into a dictionary
    for i in range(len(master_id)):
        master_dict = {}
        master_dict['master_id'] = master_id[i]
        master_dict['xpix'] = xpix[i]
        master_dict['ypix'] = ypix[i]
        master_dict['ra'] = ra[i]
        master_dict['dec'] = dec[i]

        master_dict_list.append(master_dict)
    
    for master_dict in master_dict_list:
        update_master_table(master_table, master_dict)
        print("Updated {} with row {}".format(master_table, master_dict))


#-------------------------------------------------------------------------------# 

def populate_all_master_tables():
    """
    Populates Master tables for all target clusters.

    Parameters:
        nothing

    Returns:
        nothing

    Outputs:
        Repopulated Master tables for every target cluster.
    """
    path_to_mastercats = os.path.join(path_to_data, 'master_cat', 'final')

    for targname in target_list:
        print("***Populating the Master table for {}.".format(targname))
        FileInfo, Master, Phot, Results = return_tables(targname)
        print(master_table)
        master_catalog = os.path.join(path_to_mastercats, 
            '{}_master.cat'.format(targname))
        print(master_catalog)

        populate_master_table(Master, master_catalog)


#-------------------------------------------------------------------------------# 
# Uber-wrapper.
#-------------------------------------------------------------------------------# 

def populate_all_tables():
    """ Uber-wrapper to repopulate all the Master, FileInfo, and Phot tables.
    """
    populate_all_master_tables()
    populate_all_fileinfo_phot_tables()


#-------------------------------------------------------------------------------# 
# The main.
#-------------------------------------------------------------------------------# 

if __name__=="__main__":
    #populate_all_tables()
    #populate_all_fileinfo_phot_tables()
