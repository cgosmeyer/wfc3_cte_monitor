"""
Module containing uvis_external_cte database queries.
Add functions as needed. The idea is to hide all the messiness in this 
module, and also to minimize the scripts in which we need import the
database table objects.

Authors:

    C.M. Gosmeyer (2016)

Use:

    The functions here are meant to be imported from the analysis 
    and plotting scripts.
"""

import numpy as np

from database_interface import get_session
from database_interface import load_connection
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


session = get_session()

#-------------------------------------------------------------------------------# 
# Helper functions.
#-------------------------------------------------------------------------------# 

def return_tables(targname):
    """ Returns table objects for the given target cluster in a conventient
    shorthand.

    Parameters:
        targname : string
            Name of the target cluster.

    Returns:
        FileInfo : table object
            The FileInfo table of the target cluster.
        Master : table object 
            The Master table of the target cluster.
        Phot : table object
            The Phot table of the target cluster.
        Results : table object
            The Results table of the target cluster.

    Outputs:
        nothing
    """
    if '104' in targname:
        FileInfo = NGC104_FileInfo
        Master = NGC104_Master
        Phot = NGC104_Phot
        Results = NGC104_Results
    elif '6583' in targname:
        FileInfo = NGC6583_FileInfo
        Master = NGC6583_Master
        Phot = NGC6583_Phot
        Results = NGC6583_Results
    elif '6791' in targname:
        FileInfo = NGC6791_FileInfo
        Master = NGC6791_Master
        Phot = NGC6791_Phot
        Results = NGC6791_Results

    return FileInfo, Master, Phot, Results


#-------------------------------------------------------------------------------# 
# Queries.
#-------------------------------------------------------------------------------# 

def query_for_all_dateobss(targname):
    """Quries FileInfo table for all the dateobs for a given targname.

    Parameters:
        targname : string
            Name of the target cluster.

    Returns:
        dateobs : list of floats
            MJD dates of the observations.

    Outputs:
        nothing
    """
    FileInfo, Master, Phot, Results = return_tables(targname)  

    query = session.query(FileInfo.dateobs).all()

    dateobss = [result[0] for result in query]

    return dateobss


#-------------------------------------------------------------------------------# 

def query_for_dateobss(targname, proposid, filt, exptime):
    """Queries FileInfo table for dates of observation.

    Parameters:
        targname : string
            Name of the target cluster.
        proposid : string
            Proposal number. If equals '', then query by just filter.
        filt : string
            Name of the filter.
        exptime : float
            Expsosure time.

    Returns:
        dateobss : list of floats
            MJD dates of observations.

    Outputs:
        nothing
    """
    FileInfo, Master, Phot, Results = return_tables(targname) 
    query = session.query(FileInfo.dateobs)\
            .filter(FileInfo.proposid == proposid)\
            .filter(FileInfo.filter == filt)\
            .filter(FileInfo.exptime == exptime).all()
    dateobss = [result[0] for result in query]

    return dateobss


#-------------------------------------------------------------------------------# 

def query_for_exptimes(targname, proposid, filt, dateobs=''):
    """Queries FileInfo for exposure times. Assumes querying either by 
    proposid or dateobs. If Dateobs, will query in a range of 30 days.

    Parameters:
        targname : string
            Name of the target cluster.
        proposid : string
            Proposal number. If equals '', then query by just filter.
        filt : string
            Name of the filter.
        dateobs : float
            MJD date of observation. If equals '', then query just by 
            proposal number.

    Returns:
        exptimes : list of floats
            Exposure times.

    Outputs:
        nothing
    """
    FileInfo, Master, Phot, Results = return_tables(targname) 

    if proposid == '':
        query = session.query(FileInfo.exptime)\
           .filter(FileInfo.filter == filt).all()
    elif proposid != '' and dateobs == '':
        query = session.query(FileInfo.exptime)\
                .filter(FileInfo.proposid == proposid)\
                .filter(FileInfo.filter == filt).all()

    elif dateobs != '':
        query = session.query(FileInfo.exptime)\
                .filter(FileInfo.filter == filt)\
                .filter(FileInfo.dateobs <= dateobs+30.)\
                .filter(FileInfo.dateobs >= dateobs-30.).all()

    exptimes = [result[0] for result in query]

    return exptimes


#-------------------------------------------------------------------------------# 

def query_for_flux_by_imagename(targname, imagename, aperture, 
    bkgrd_returned='tot'):
    """Queries Phot table for fluxes, backgrounds, y-positions, and
    master ids.

    Parameters:
        targname : string
            Name of the target cluster.
        imagename : string
            Name of the image. (presumably found from `query_for_pair`)
        aperture : int
            The aperture of photometry. (3 is nominal)
        bkgrd_returned : string
            Either 'tot' for the total background or 'local' for the
            local background.

    Returns:
        fluxes : list of floats
            Fluxes of the sources for the given aperture.
        bkgrds : list of floats
            Backgrounds of the sources.
        ypix : list of ints
            Y-positions in pixels of the sources.
        master_id : list of ints
            Master IDs of the sources.

    Outputs:
        nothing
    """
    FileInfo, Master, Phot, Results = return_tables(targname) 

    flux_column_dict = {2:Phot.flux_2, 3:Phot.flux_3, 5:Phot.flux_5, 
        7:Phot.flux_7, 10:Phot.flux_10, 12:Phot.flux_12, 15:Phot.flux_15,
        18:Phot.flux_18, 20:Phot.flux_20, 24:Phot.flux_24, 28:Phot.flux_28,
        32:Phot.flux_32, 36:Phot.flux_36, 40:Phot.flux_40}
    totbkgrd_column_dict = {2:Phot.totbkgrd_2, 3:Phot.totbkgrd_3, 
        5:Phot.totbkgrd_5, 7:Phot.totbkgrd_7, 10:Phot.totbkgrd_10, 
        12:Phot.totbkgrd_12, 15:Phot.totbkgrd_15, 18:Phot.totbkgrd_18, 
        20:Phot.totbkgrd_20, 24:Phot.totbkgrd_24, 28:Phot.totbkgrd_28,
        32:Phot.totbkgrd_32, 36:Phot.totbkgrd_36, 40:Phot.totbkgrd_40}
    mnbkgrd_column_dict = {2:Phot.mnbkgrd_2, 3:Phot.mnbkgrd_3, 
        5:Phot.mnbkgrd_5, 7:Phot.mnbkgrd_7, 10:Phot.mnbkgrd_10, 
        12:Phot.mnbkgrd_12, 15:Phot.mnbkgrd_15, 18:Phot.mnbkgrd_18, 
        20:Phot.mnbkgrd_20, 24:Phot.mnbkgrd_24, 28:Phot.mnbkgrd_28,
        32:Phot.mnbkgrd_32, 36:Phot.mnbkgrd_36, 40:Phot.mnbkgrd_40}

    if bkgrd_returned == 'tot':
        bkgrd_query = totbkgrd_column_dict[aperture]
    else:
        bkgrd_query = mnbkgrd_column_dict[aperture]

    query = session.query(flux_column_dict[aperture], bkgrd_query,
            Phot.ypix, Phot.master_id)\
            .filter(Phot.imagename == imagename).all()

    fluxes = [result[0] for result in query]
    bkgrds = [result[1] for result in query]
    ypix = [result[2] for result in query]
    master_id = [result[3] for result in query]

    return fluxes, bkgrds, ypix, master_id


#-------------------------------------------------------------------------------# 

def query_for_flux_by_masterid(targname, master_id, ctecorr, filt, chinject, 
                               exptime, flashlvl, chip):
    """ Takes in an imagename (presumably found from `query_for_pair`).
    Returns list of the photometry results and background.

    Filter by imagenames (which have already been filtered by post-flash,
        ctecorr, chip1, etc.)

    Parameters:
        targname : string
            Name of the target cluster.
        master_id : int
            ID of source from Master table.
        ctecorr : int
            (0,1) CTE correction performed.
        filt : string
            Name of the filter.
        chinject : string
            Charge injection.
        exptime : float
            Expsosure time.
        flashlvl : float
            Flash level.
        chip : int
            (1,2) Number of the chip.

    Returns:
        imagenames : list of strings
            List of imagenames that match given parameters.
        fluxes_all : list of lists
            For each imagename, a list of the fluxes.

    Outputs:    
        nothing
    """
    FileInfo, Master, Phot, Results = return_tables(targname) 

    query = session.query(FileInfo.imagename)\
        .filter(FileInfo.filter == filt)\
        .filter(FileInfo.exptime == exptime)\
        .filter(FileInfo.chinject == chinject)\
        .filter(FileInfo.flashlvl == flashlvl)\
        .filter(FileInfo.ctecorr == ctecorr)\
        .filter(FileInfo.chip == 1).all()

    imagenames = [result[0] for result in query]
    print(imagenames)

    fluxes_all = []
    for imagename in imagenames:
        query = session.query(Phot.flux_10)\
            .filter(Phot.master_id == master_id)\
            .filter(Phot.imagename == imagename).all()

        fluxes = [result[0] for result in query]
        print(fluxes)
        fluxes_all.append(fluxes)

    return imagenames, fluxes_all


#-------------------------------------------------------------------------------# 

def query_for_flux_range(targname, dateobs, exptime, flashlvl, ctecorr, filt,
                         chinject, postarg1, chip, flux_lo, flux_hi, aperture,
                         subtract_background=False):
    """Queries Phot table for all fluxes in given range 'flux_lo' and 
    'flux_hi'.

    Note that dateobs searches for all observations from 30 days before
    and after given date!

    Parameters:
        targname : string
            Name of the target cluster.
        dateobs : float
            MJD date of the observations.
        exptime : float
            Expsosure time.
        flashlvl : float
            Flash level.
        ctecorr : int
            (0,1) CTE correction performed.            
        filt : string
            Name of the filter.
        chinject : string
            Charge injection.
        postarg1 : float
            Postarg1.
        chip : int
            (1,2) Number of chip.
        flux_lo : float
            Lower end of fluxbin to query.
        flux_hi : float
            Upper end of fluxbin to query.
        aperture : int
            The aperture of photometry. (3 is nominal)
        subtract_background : {True, False}
            Set on to subtract backgrounds from fluxes before
            sorting by range.
            Doesn't seem to change logflux plots. Not recommended
            since makes code 2-3x slower.
  
    Returns:
        fluxes_per_image : list of lists
            All fluxes of the image that fall in given range.
        imagenames : list of strings
            Names of images.

    Outputs: 
        nothing
    """
    FileInfo, Master, Phot, Results = return_tables(targname)  

    # First find imagenames.
    query = session.query(FileInfo.imagename)\
        .filter(FileInfo.dateobs <= dateobs+30.)\
        .filter(FileInfo.dateobs >= dateobs-30.)\
        .filter(FileInfo.exptime == exptime)\
        .filter(FileInfo.flashlvl == flashlvl)\
        .filter(FileInfo.ctecorr == ctecorr)\
        .filter(FileInfo.filter == filt)\
        .filter(FileInfo.chinject == chinject)\
        .filter(FileInfo.postarg1 == postarg1)\
        .filter(FileInfo.chip == chip).all()

    imagenames = [result[0] for result in query]


    # Second find all fluxes in that image in the desired range.
    
    flux_column_dict = {2:Phot.flux_2, 3:Phot.flux_3, 5:Phot.flux_5, 
        7:Phot.flux_7, 10:Phot.flux_10, 12:Phot.flux_12, 15:Phot.flux_15,
        18:Phot.flux_18, 20:Phot.flux_20, 24:Phot.flux_24, 28:Phot.flux_28,
        32:Phot.flux_32, 36:Phot.flux_36, 40:Phot.flux_40}

    mnbkgrd_column_dict = {2:Phot.mnbkgrd_2, 3:Phot.mnbkgrd_3, 
        5:Phot.mnbkgrd_5, 7:Phot.mnbkgrd_7, 10:Phot.mnbkgrd_10, 
        12:Phot.mnbkgrd_12, 15:Phot.mnbkgrd_15, 18:Phot.mnbkgrd_18, 
        20:Phot.mnbkgrd_20, 24:Phot.mnbkgrd_24, 28:Phot.mnbkgrd_28,
        32:Phot.mnbkgrd_32, 36:Phot.mnbkgrd_36, 40:Phot.mnbkgrd_40}

    fluxes_per_image = []

    if subtract_background:
        # This does not appear to make any difference at all in the final
        # logflux plot. Only serves to make code take 2-3x longer.
        for imagename in imagenames:
            query = session.query(flux_column_dict[aperture], mnbkgrd_column_dict[aperture])
            # Now subtract backgrounds and then filter by flux low and high.
            fluxes = [result[0] for result in query]
            bkgrds = [result[1] for result in query]
            flux_sub_bkgrds = np.array(fluxes) - np.array(bkgrds)
            good_fluxes = []
            for flux_sub_bkgrd in flux_sub_bkgrds:
                if ((flux_sub_bkgrd <= flux_hi) and (flux_sub_bkgrd >= flux_lo)):
                    good_fluxes.append(flux_sub_bkgrd)
            fluxes_per_image.append(good_fluxes)

    elif not subtract_background:
        for imagename in imagenames:
            query = session.query(flux_column_dict[aperture])\
                .filter(flux_column_dict[aperture] <= flux_hi)\
                .filter(flux_column_dict[aperture] >= flux_lo).all()

            fluxes = [result[0] for result in query]
            fluxes_per_image.append(fluxes)

    return fluxes_per_image, imagenames 


#-------------------------------------------------------------------------------# 

def query_for_fluxes_bkgrds_by_ypos(targname, proposid, filt, ctecorr, 
    exptime, chip, chinject='NO', postarg1=0, ypos=1750.):
    """Queries fluxes and backgrounds from Phot tables. 
    Return results for ALL flashlvls.

    Parameters:
        targname : string
            Name of the target cluster.
        proposid : string
            Proposal number.
        filt : string
            Name of the filter.
        ctecorr : int
            (0,1) CTE correction performed.
        exptime : float
            Expsosure time.
        chip : int
            (1,2) Number of chip.
        chinject : string
            Charge injection.
        flashlvl : float
            Flash level.
        postarg1 : float
            Postarg1.
        ypos : float/int
            Y-position in pixels above which to query.

    Returns:
        image_dict : dictionary
            Keys of imagenames. Values of lists of fluxes, local backgrounds,
            global backgrounds, and y-positions of all sources in image.

    Outputs:
        nothing
    """
    FileInfo, Master, Phot, Results = return_tables(targname) 

    # Get list of imagenames that match criteria.
    query_imagenames = session.query(FileInfo.imagename)\
        .filter(FileInfo.proposid == proposid)\
        .filter(FileInfo.filter == filt)\
        .filter(FileInfo.exptime == exptime)\
        .filter(FileInfo.chinject == chinject)\
        .filter(FileInfo.ctecorr == ctecorr)\
        .filter(FileInfo.postarg1 == postarg1)\
        .filter(FileInfo.chip == chip).all()

    imagenames = [result[0] for result in query_imagenames]

    image_dict = {}
    # Get lists of all sources that are y > 1750 pixels.
    for imagename in imagenames:
        query = session.query(Phot.flux_3, Phot.mnbkgrd_3, 
            Phot.totbkgrd_3, Phot.ypix)\
            .filter(Phot.imagename == imagename)\
            .filter(Phot.ypix >= ypos).all()

        fluxes = [result[0] for result in query]
        mnbkgrds = [result[1] for result in query]
        totbkgrds = [result[2] for result in query]
        ypixs = [result[3] for result in query]
        print(len(fluxes))
        image_dict[imagename] = [fluxes, mnbkgrds, totbkgrds, ypixs]

    return image_dict


#-------------------------------------------------------------------------------# 

def query_for_globalbkgrd(targname, imagename):
    """Queries FileInfo table for global background of given image.

    Parameters:
        targname : string
            Name of the target cluster.
        imagename : string
            Name of the image.

    Returns:
        global_bkgrds : list of floats
            Global background of image. (kept as list for consistency
            with other querying functions)

    Outputs:
        nothing
    """
    FileInfo, Master, Phot, Results = return_tables(targname)    

    # Use the FileInfo table to find the master id and imagename.
    # use that to query the Phot table for flux, etc
    query = session.query(FileInfo.mnclip_bkgrd)\
        .filter(FileInfo.imagename == imagename).all()

    global_bkgrds = [result[0] for result in query]

    return global_bkgrds


#-------------------------------------------------------------------------------# 

def query_for_matching_imagename(targname, imagename_1):
    """Queries for imagename that pairs with 'imagename_1'.

    Parameters:
        targname : string
            Name of the target cluster.
        imagename_1 : string
            Name of the image.

    Returns:
        imagename_2 : string
            Name of image that pairs with 'imagename_1'

    Outputs:
        nothing
    """
    FileInfo, Master, Phot, Results = return_tables(targname)  

    # Need to get parameters of first image in order to query for the second.
    query = session.query(FileInfo.proposid, FileInfo.dateobs, FileInfo.filter, 
        FileInfo.exptime, FileInfo.chinject, FileInfo.flashlvl, 
        FileInfo.ctecorr, FileInfo.postarg1, FileInfo.chip)\
        .filter(FileInfo.imagename == imagename_1).all()

    proposid = [result[0] for result in query][0]
    dateobs = [result[1] for result in query][0]
    filt = [result[2] for result in query][0]
    exptime = [result[3] for result in query][0] 
    chinject = [result[4] for result in query][0] 
    flashlvl = [result[5] for result in query][0] 
    ctecorr = [result[6] for result in query][0] 
    postarg1 = [results[7] for results in query][0]
    chip_1 = [result[8] for result in query][0] 

    if chip_1 == 1:
        chip_2 = 2
    elif chip_1 == 2:
        chip_2 == 1

    # Now query for the second imagename.

    query = session.query(FileInfo.imagename)\
        .filter(FileInfo.proposid == proposid)\
        .filter(FileInfo.dateobs == dateobs)\
        .filter(FileInfo.filter == filt)\
        .filter(FileInfo.exptime == exptime)\
        .filter(FileInfo.chinject == chinject)\
        .filter(FileInfo.flashlvl == flashlvl)\
        .filter(FileInfo.ctecorr == ctecorr)\
        .filter(FileInfo.postarg1 == postarg1)\
        .filter(FileInfo.chip == chip_2).all()

    if query != []:
        imagename_2 = [result[0] for result in query][0]
    else:
        imagename_2 = []

    return imagename_2


#-------------------------------------------------------------------------------# 

def query_for_pair(targname, proposid, dateobs, filt, exptime, 
                   chinject, flashlvl, ctecorr, postarg1):
    """ Queries for a matching chip1-chip2 image pair from the nominal 
    dataset. Most of the parameters are directly from the FITS headers.
    
    Parameters:
        targname : string
            Name of the target cluster.
        proposid : string
            Proposal number.
        dateobs : float
            MJD date of the observations.
        filt : string
            Name of the filter.
        exptime : float
            Expsosure time.
        chinject : string
            Charge injection.
        flashlvl : float
            Flash level.
        ctecorr : int
            (0,1) CTE correction performed.
        postarg1 : float
            Postarg1.

    Returns:
        imagenames_chip1 : list of strings
            List of imagenames that match given parameters for chip 1.
        imagenames_chip1 : list of strings
            List of imagenames that match given parameters for chip 2.

    Outputs:    
        nothing
    """
    FileInfo, Master, Phot, Results = return_tables(targname)    

    # Use the FileInfo table to find the master id and imagename.
    # use that to query the Phot table for flux, etc
    if proposid != '':
        query_chip1 = session.query(FileInfo.imagename)\
            .filter(FileInfo.proposid == proposid)\
            .filter(FileInfo.dateobs == dateobs)\
            .filter(FileInfo.filter == filt)\
            .filter(FileInfo.exptime == exptime)\
            .filter(FileInfo.chinject == chinject)\
            .filter(FileInfo.flashlvl == flashlvl)\
            .filter(FileInfo.ctecorr == ctecorr)\
            .filter(FileInfo.postarg1 == postarg1)\
            .filter(FileInfo.chip == 1).all()

        imagenames_chip1 = [result[0] for result in query_chip1]

        query_chip2 = session.query(FileInfo.imagename)\
            .filter(FileInfo.proposid == proposid)\
            .filter(FileInfo.dateobs == dateobs)\
            .filter(FileInfo.filter == filt)\
            .filter(FileInfo.exptime == exptime)\
            .filter(FileInfo.chinject == chinject)\
            .filter(FileInfo.flashlvl == flashlvl)\
            .filter(FileInfo.ctecorr == ctecorr)\
            .filter(FileInfo.postarg1 == postarg1)\
            .filter(FileInfo.chip == 2).all()
        imagenames_chip2 = [result[0] for result in query_chip2]

    else:
        query_chip1 = session.query(FileInfo.imagename)\
            .filter(FileInfo.dateobs == dateobs)\
            .filter(FileInfo.filter == filt)\
            .filter(FileInfo.exptime == exptime)\
            .filter(FileInfo.chinject == chinject)\
            .filter(FileInfo.flashlvl == flashlvl)\
            .filter(FileInfo.ctecorr == ctecorr)\
            .filter(FileInfo.postarg1 == postarg1)\
            .filter(FileInfo.chip == 1).all()

        imagenames_chip1 = [result[0] for result in query_chip1]

        query_chip2 = session.query(FileInfo.imagename)\
            .filter(FileInfo.dateobs == dateobs)\
            .filter(FileInfo.filter == filt)\
            .filter(FileInfo.exptime == exptime)\
            .filter(FileInfo.chinject == chinject)\
            .filter(FileInfo.flashlvl == flashlvl)\
            .filter(FileInfo.ctecorr == ctecorr)\
            .filter(FileInfo.postarg1 == postarg1)\
            .filter(FileInfo.chip == 2).all()
        imagenames_chip2 = [result[0] for result in query_chip2]

    # THis should return two files, one for each member of the pair
    return imagenames_chip1, imagenames_chip2


#-------------------------------------------------------------------------------# 

def query_for_180pair(targname, filt, exptime, ctecorr, chip):
    """
    Visits 10 and 11 of program 12692.  Each Visit has 7 exposures.

    The second visit, 11, is rolled 180 degrees, but is otherwise identical.
    None of them are postflashed. 

    So the pairs should match in everything but visit number and postarg1!

    Exposures 4,5,7 are dithered up so that (Chip 2) overlaps with the second visit.
    Exposures 1,2,3,6 are not dithered, so that (Chip 1) overlaps with second visit.
    (therefore can use postarg1 to distinguish pairs)
    Exposure 6 is in F606W. The others are in F502N. 


    Parameters:
        targname : string
            Name of the target cluster.
        filt : string
            Name of the filter.
        exptime : float
            Expsosure time.
        ctecorr : int
            (0,1) CTE correction performed.
        chip : int
            (1,2) Number of the chip.

    Returns:
        imagenames_first : list of strings
            List of imagenames that match given parameters for first 
            visit.
        imagenames_second : list of strings
            List of imagenames that match given parameters for second 
            visit rolled 180 degrees.

    Outputs:    
        nothing
    """
    FileInfo, Master, Phot, Results = return_tables(targname) 
    
    # Do an initial query for the imagenames and ids, then sort into 1st and 2nd
    # chips based on the id and visit (found from the imagenames)
    query = session.query(FileInfo.imagename)\
            .filter(FileInfo.filter == filt)\
            .filter(FileInfo.exptime == exptime)\
            .filter(FileInfo.ctecorr == ctecorr)\
            .filter(FileInfo.chip == chip).all()

    imagenames = [result[0] for result in query]  

    # To make sure that the pairs are correct, need to look at the
    # order of imagenames
    imagenames_sorted = sorted(imagenames)
    imagenames_first = [im for im in imagenames_sorted if '10' in im]
    imagenames_second = [im for im in imagenames_sorted if '11' in im]

    # The imagenames in both lists are of the same chip.
    return imagenames_first, imagenames_second


#-------------------------------------------------------------------------------# 

def query_results_for_slopes(targname, imagename_chip1, imagename_chip2, 
    fluxbin, aperture):
    """Queries Results table for CTE slopes.

    Parameters:
        targname : string
            Name of the target cluster.
        imagename_chip1 : string
            Name of the chip1 image.
        imagename_chip2 : string
            Name of the chip2 image.
        fluxbin : string
            Name of the fluxbin. ie, '250-500'
        aperture : int
            The aperture of photometry. (3 is nominal)

    Returns:
        slopes : list of floats
            CTE slopes for given fluxbin.
        slopesdtdevs : list of floats
            Standard deviations for CTE slopes in given fluxbin.
        numpoints : list of ints
            Number of sources in the given fluxbin.

    Outputs:
        nothing
    """
    FileInfo, Master, Phot, Results = return_tables(targname)   
    
    slope_column_dict = {'250-500':Results.slope_250_500,
                         '500-1000':Results.slope_500_1000,
                         '500-2000':Results.slope_500_2000,
                         '1000-2000':Results.slope_1000_2000,
                         '2000-4000':Results.slope_2000_4000,
                         '2000-8000':Results.slope_2000_8000,
                         '4000-8000':Results.slope_4000_8000,
                         '8000-32000':Results.slope_8000_32000}

    slopestdev_column_dict = {'250-500':Results.slopestdev_250_500,
                         '500-1000':Results.slopestdev_500_1000,
                         '500-2000':Results.slopestdev_500_2000,
                         '1000-2000':Results.slopestdev_1000_2000,
                         '2000-4000':Results.slopestdev_2000_4000,
                         '2000-8000':Results.slopestdev_2000_8000,
                         '4000-8000':Results.slopestdev_4000_8000,
                         '8000-32000':Results.slopestdev_8000_32000}

    numpoints_column_dict = {'250-500':Results.numpoints_250_500,
                         '500-1000':Results.numpoints_500_1000,
                         '500-2000':Results.numpoints_500_2000,
                         '1000-2000':Results.numpoints_1000_2000,
                         '2000-4000':Results.numpoints_2000_4000,
                         '2000-8000':Results.numpoints_2000_8000,
                         '4000-8000':Results.numpoints_4000_8000,
                         '8000-32000':Results.numpoints_8000_32000}

    if fluxbin not in slope_column_dict.keys():
        print("Invalid fluxbin of {}. Empty list returned.".format(fluxbin))
        return []
    
    print(imagename_chip1)
    print(imagename_chip2)
    print(aperture)
    query = session.query(slope_column_dict[fluxbin], 
            slopestdev_column_dict[fluxbin], 
            numpoints_column_dict[fluxbin])\
        .filter(Results.imagename_1 == imagename_chip1)\
        .filter(Results.imagename_2 == imagename_chip2)\
        .filter(Results.aperture == aperture).all()   

    slopes = [result[0] for result in query]
    slopesdtdevs = [result[1] for result in query]
    numpoints = [result[2] for result in query]

    return slopes, slopesdtdevs, numpoints


#-------------------------------------------------------------------------------# 
# Main for testing query functions.
#-------------------------------------------------------------------------------# 

if __name__=='__main__':
    

    targname='ngc104'
    proposid='12348'
    dateobs=55466
    filt='F606W'
    exptime=30
    chinject='NO'
    flashlvl=0
    ctecorr=0
    chip = 1
    postarg1=0

    #imagenames_chip1,imagenames_chip2 = query_for_pair(targname, proposid, 
    #                   dateobs, filt, exptime, chinject, flashlvl, ctecorr)

    #query_for_flux_by_imagename(targname, imagenames_chip1[0])

    #query_for_flux_by_masterid(targname, 2146, ctecorr, filt, 
    #    chinject, exptime, flashlvl, chip)

    #query_for_fluxes_bkgrds_by_ypos(targname, filt, ctecorr, exptime, chip,
    #    chinject='NO', postarg1=0)

    
    first, second = query_for_180pair(targname = 'ngc6583', filt='F606W', exptime=348.0, ctecorr=True, chip=1)
    print(first)
    print(second)
