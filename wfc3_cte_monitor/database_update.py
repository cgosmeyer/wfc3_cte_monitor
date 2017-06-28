"""
This module serves as an interface for updating the various tables of
the uvis_external_cte database, either by inserting new records, or updating
existing ones.

Authors:

    C.M. Gosmeyer, April 2016

Use:

    This module is intended to be imported from the various uvis_external_cte
    scripts, as such:

    ::

    from database_update import update_master_table
    from database_update import update_fileinfo_table
    from database_update import update_phot_table
    from database_update import update_results_table

Notes:

    Based on M.Bourque's lightcurve_pipeline `update_database.py`
    module.
"""

import datetime
import numpy as np
import os

from astropy.io import ascii

from config import path_to_data
from photutils_plus.phot_tools import read_in_photfile

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


#-------------------------------------------------------------------------------# 
# Functions to update database tables.
#-------------------------------------------------------------------------------# 

def update_fileinfo_table(FileInfo, fileinfo_dict):
    """Insert or update a record in a file information table.
    Each row needs be populated before the linked rows in the Phot
    table. They are linked by 'imagename'.

    Parameters:
        FileInfo : Table 
            The table in which to insert 'fileinfo_dict'.
        fileinfo_dict : dict
            A dictionary containing the file information. Each key of
            'fileinfo_dict' corresponds to a column in a FileInfo
            table of the database.

    Returns:
        nothing

    Outputs:
        nothing

    Notes:
        Based on M.Bourque's `update_database.py` functions.
    """

    # Get the id of the record, if it exists
    session = get_session()
    query = session.query(FileInfo.id)\
        .filter(FileInfo.imagename == fileinfo_dict['imagename']).all()
    if query == []:
        id_num = ''
    else:
        id_num = query[0][0]
    session.close()

    # If id doesn't exist then insert. If id exsits, then update
    insert_or_update(FileInfo, fileinfo_dict, id_num)


#-------------------------------------------------------------------------------# 

def update_master_table(Master, master_dict):
    """Insert or update a record in a master table.

    Unless rebuilding databse, only need to populate the master tables 
    once. 

    Parameters:
        Master : Table object
            The table in which to insert ``master_dict``.
        master_dict : dict
            A dictionary containing the master information. Each key of 
            ``master_dict`` corresponds to a column in a master
            table of the database.

    Returns:
        nothing

    Outputs:
        nothing

    Notes:
        Based on M.Bourque's `update_database.py` functions.

    """

    # Get the id of the record, if it exists
    session = get_session()
    query = session.query(Master.id)\
        .filter(Master.id == master_dict['master_id']).all()
    if query == []:
        id_num = ''
    else:
        id_num = query[0][0]
    session.close()

    # If id doesn't exist then insert.  If id exists, then update
    insert_or_update(Master, master_dict, id_num)


#-------------------------------------------------------------------------------# 

def update_phot_table(Phot, phot_dict):
    """Insert or update a record in a Phot table.

    Parameters:
        Phot : Table 
            The table in which to insert ``phot_dict``.
        phot_dict : dict
            A dictionary containing the photometry. Each key of 
            ``phot_dict`` corresponds to a column in a Phot
            table of the database.

    Returns:
        nothing

    Outputs:
        nothing

    Notes:
        Based on M.Bourque's `update_database.py` functions.
    """

    # Get the id of the record, if it exists
    session = get_session()
    query = session.query(Phot.id)\
        .filter(Phot.find_id == phot_dict['find_id'])\
        .filter(Phot.imagename == phot_dict['imagename']).all()
    if query == []:
        id_num = ''
    else:
        id_num = query[0][0]
    session.close()

    # If id doesn't exist then insert. If id exsits, then update
    insert_or_update(Phot, phot_dict, id_num)


#-------------------------------------------------------------------------------# 

def update_results_table(Results, results_dict):
    """Insert or update a record in a results table.

    Parameters:
        Results : Table 
            The table in which to insert ``results_dict``.
        results_dict : dict
            A dictionary containing the file information. Each key of
            ``results_dict`` corresponds to a column in a results
            table of the database.

    Returns:
        nothing

    Outputs:
        nothing

    Notes:
        Based on M.Bourque's `update_database.py` functions.
    """

    # Get the id of the record, if it exists
    session = get_session()
    query = session.query(Results.id)\
        .filter(Results.imagename_1 == results_dict['imagename_1'])\
        .filter(Results.imagename_2 == results_dict['imagename_2'])\
        .filter(Results.aperture == results_dict['aperture']).all()
    if query == []:
        id_num = ''
    else:
        id_num = query[0][0]
    session.close()

    # If id doesn't exist then insert. If id exsits, then update
    insert_or_update(Results, results_dict, id_num)


#-------------------------------------------------------------------------------# 

def insert_or_update(table, data, id_num):
    """
    Insert or update the given database table with the given data.
    This function performs the logic of inserting or updating an
    entry into the hstlc database; if an entry with the given
    'id_num' already exists, then the entry is updated, otherwise a
    new entry is inserted.

    Parameters:
        table : sqlalchemy.ext.declarative.api.DeclarativeMeta
            The table of the database to update
        data : dict
            A dictionary of the information to update.  Each key of the
            dictionary must be a column in the given table
        id_num : string
            The row ID to update.  If ``id_num`` is blank, then a new row
            is inserted instead.

    Returns:
        nothing

    Outputs:
        nothing

    Notes:
        Based on M.Bourque's `update_database.py` functions.
    """

    from database_interface import engine
    from database_interface import get_session

    session = get_session()
    if id_num == '':
        print("Creating entry for table {}.".format(table))
        engine.execute(table.__table__.insert(), data)
    else:
        print("Updating entry for table {}.".format(table))
        session.query(table)\
            .filter(table.id == id_num)\
            .update(data)
        session.commit()
    session.close()


#-------------------------------------------------------------------------------# 
# Functions to create dictionaries that can be passed into database tables.
#-------------------------------------------------------------------------------# 

def return_fileinfo_dict(imagename, coofile, magfile, path_to_photom, 
                         param_dict):
    """ 
    Creates dictionary of columns and values for FileInfo table.
    Each coo and mag file pair only needs one row. 

    Parameters:
        imagename : string
            Name of the image.
        coofile : string
            Name of the *coo file.
        magfile : string
            Name of the *mag file.
        path_to_photom : string
            Path to the photometry (coo and mag) files.
        param_dict : dict
            Dictionary of parameters of the image.
            
    Returns: 
        fileinfo_dict : dict
            Dictionary of column names and values for given image. This
            will be inserted into database.

    Outputs:
        nothing
    """
    parsed_name = return_parsedname(param_dict)

    fileinfo_dict = {}

    fileinfo_dict['imagename'] = imagename.split('/')[-1]
    fileinfo_dict['imagepath'] = '/'.join((imagename.split('/'))[0:-1])
    fileinfo_dict['coofile'] = coofile.split('/')[-1]
    fileinfo_dict['magfile'] = magfile.split('/')[-1]
    fileinfo_dict['photpath'] = path_to_photom
    fileinfo_dict['parsed_name'] = parsed_name
    fileinfo_dict['ingest_date'] = datetime.date.today() 
    fileinfo_dict['ra_lowerleft'] = param_dict['ra_lowerleft']
    fileinfo_dict['dec_lowerleft'] = param_dict['dec_lowerleft']
    fileinfo_dict['ra_lowerright'] = param_dict['ra_lowerright']
    fileinfo_dict['dec_lowerright'] = param_dict['dec_lowerright']
    fileinfo_dict['ra_upperleft'] = param_dict['ra_upperleft']
    fileinfo_dict['dec_upperleft'] = param_dict['dec_upperleft']
    fileinfo_dict['ra_upperright'] = param_dict['ra_upperright']
    fileinfo_dict['dec_upperright'] = param_dict['dec_upperright']
    fileinfo_dict['mnclip_bkgrd'] = param_dict['mnclip_bkgrd']
    fileinfo_dict['proposid'] = param_dict['proposid']
    fileinfo_dict['dateobs'] = param_dict['dateobs']
    fileinfo_dict['filter'] = param_dict['filter']
    fileinfo_dict['exptime'] = param_dict['exptime']
    fileinfo_dict['chinject'] = param_dict['chinject']
    fileinfo_dict['flashlvl'] = param_dict['flashlvl']
    fileinfo_dict['ctecorr'] = param_dict['ctecorr']
    fileinfo_dict['chip'] = param_dict['chip']
    fileinfo_dict['flashdur'] = param_dict['flashdur']
    fileinfo_dict['flashcur'] = param_dict['flashcur']
    fileinfo_dict['shutrpos'] = param_dict['shutrpos']
    fileinfo_dict['postarg1'] = param_dict['postarg1']
    fileinfo_dict['postarg2'] = param_dict['postarg2']

    return fileinfo_dict


#-------------------------------------------------------------------------------# 

def return_phot_dict_list(imagename, coofile, magfile, param_dict, 
    use_iraf=False):
    """ 
    Creates dictionary of columns and values for Phot table.
    Need a row for each file, for each master ID.

    Parameters:
        imagename : string
            Name of the image.
        coofile : string
            Name of the *coo file.
        magfile : string
            Name of the *mag file.
        param_dict : dict
            Dictionary of parameters of the image.
        use_iraf : {True, False}
            Run photometry with ``iraf.phot`` instead of ``phututils.
            aperture_photometry`` and insert into 
            'uvis_external_cte_iraf.db'.
            
    Returns: 
        phot_dict_list : list of dictionaries
            List of dictionaries for each source of photometry results.
            This will be inserted into database.

    Outputs:
        nothing
    """
    # read in coofile
    coo_table, coo_header = read_in_photfile(coofile, isheader=False)

    # read in magfile
    if use_iraf:
        mag_table = ascii.read(magfile, format='daophot')
    else:
        mag_table, mag_header = read_in_photfile(magfile, isheader=False)

    if not use_iraf:
        # This will be a little tricky because each source 
        # has a seperate row for each aperture.
        master_ids_present_raw = list(coo_table['master_id'])
        ra_raw = list(coo_table['extr_ra'])
        dec_raw = list(coo_table['extr_dec'])
        find_id_raw = list(mag_table['ID'])
        radius_raw = list(mag_table['radius'])
        aperture_sum_raw = list(mag_table['aperture_sum'])
        xcenter_raw = list(mag_table['xcenter'])
        ycenter_raw = list(mag_table['ycenter'])
        mean_local_bkgrd_raw = list(mag_table['mean_local_bkgrd'])
        tot_local_bkgrd_raw = list(mag_table['tot_local_bkgrd'])

        # Probs better way to do this, but for now...
        master_ids_present = []    
        find_ids = []
        radius = []
        xpixs = []
        ypixs = []
        ras = []
        decs = []
        flux_2 = []
        flux_3 = []
        flux_5 = []
        flux_7 = []
        flux_10 = []
        flux_12 = []
        flux_15 = []
        flux_18 = []
        flux_20 = []
        flux_24 = []
        flux_28 = []
        flux_32 = []
        flux_36 = []
        flux_40 = []
        mnbkgrd_2 = []
        mnbkgrd_3 = []
        mnbkgrd_5 = []
        mnbkgrd_7 = []
        mnbkgrd_10 = []
        mnbkgrd_12 = []
        mnbkgrd_15 = []
        mnbkgrd_18 = []
        mnbkgrd_20 = []
        mnbkgrd_24 = []
        mnbkgrd_28 = []
        mnbkgrd_32 = []
        mnbkgrd_36 = []
        mnbkgrd_40 = []
        totbkgrd_2 = []
        totbkgrd_3 = []
        totbkgrd_5 = []
        totbkgrd_7 = []
        totbkgrd_10 = []
        totbkgrd_12 = []
        totbkgrd_15 = []
        totbkgrd_18 = []
        totbkgrd_20 = []
        totbkgrd_24 = []
        totbkgrd_28 = []
        totbkgrd_32 = []
        totbkgrd_36 = []
        totbkgrd_40 = []

        # Compare using the 'radius' vs list of known radii.
        # when hit known_radii[-1] know to set the master_count +=1 
        flux_lists = [flux_2, flux_3, flux_5, flux_7, flux_10, flux_12, 
                      flux_15, flux_18, flux_20, flux_24, flux_28, 
                      flux_32, flux_36, flux_40]
        mnbkgrd_lists = [mnbkgrd_2, mnbkgrd_3, mnbkgrd_5, mnbkgrd_7, mnbkgrd_10, mnbkgrd_12,
                      mnbkgrd_15, mnbkgrd_18, mnbkgrd_20, mnbkgrd_24, mnbkgrd_28, 
                      mnbkgrd_32, mnbkgrd_36, mnbkgrd_40]
        totbkgrd_lists = [totbkgrd_2, totbkgrd_3, totbkgrd_5, totbkgrd_7, totbkgrd_10, totbkgrd_12,
                      totbkgrd_15, totbkgrd_18, totbkgrd_20, totbkgrd_24, totbkgrd_28, 
                      totbkgrd_32, totbkgrd_36, totbkgrd_40]
        known_radii = [2,3,5,7,10,12,15,18,20,24,28,32,36,40]


        radii_count = 0
        master_count = 0
        for i in range(len(find_id_raw)):
            print(i)

            print("radius_raw[i]: {}".format(radius_raw[i]))
            print("known_radii[radii_count]: {}".format(known_radii[radii_count]))
            if (int(radius_raw[i]) == int(known_radii[radii_count])) and \
                (radii_count != len(known_radii)-1):
                print("appending to flux list {}".format(radius_raw[i]))
                flux_lists[radii_count].append(aperture_sum_raw[i])
                mnbkgrd_lists[radii_count].append(mean_local_bkgrd_raw[i])
                totbkgrd_lists[radii_count].append(tot_local_bkgrd_raw[i])
                radii_count += 1
                print("increasing radius count to {}".format(radii_count))

            elif (int(radius_raw[i]) == int(known_radii[radii_count])) and \
                (radii_count == len(known_radii)-1):
                print("appending to flux list {}".format(radius_raw[i]))

                flux_lists[radii_count].append(aperture_sum_raw[i])
                mnbkgrd_lists[radii_count].append(mean_local_bkgrd_raw[i])
                totbkgrd_lists[radii_count].append(tot_local_bkgrd_raw[i])
                # Only need to append to the following once since the 
                # values are the same for all ap radii of given source.
                master_ids_present.append(master_ids_present_raw[master_count])
                ras.append(ra_raw[master_count])
                decs.append(dec_raw[master_count])
                find_ids.append(find_id_raw[i])
                radius.append(radius_raw[i])
                xpixs.append(xcenter_raw[i])
                ypixs.append(ycenter_raw[i])

                radii_count = 0
                master_count += 1   
                print("increasing master count to {}".format(master_count))        


    if use_iraf:
        master_ids_present = list(coo_table['master_id'])
        ras = list(coo_table['extr_ra'])
        decs = list(coo_table['extr_dec'])

        find_ids = list(mag_table['ID'])
        xpixs = list(mag_table['XCENTER'])
        ypixs = list(mag_table['YCENTER'])
        mnbkgrds = list(mag_table['MSKY'])

        flux_lists = []
        totbkgrd_lists = []
        area_lists = [] # will need multiple the mbkgrds by this to get totbkgrds
        radius_lists = []
        mnbkgrd_lists = []

        colnames = mag_table.colnames
        for col in colnames:
            if 'FLUX' in col:
                flux_lists.append(list(mag_table[col]))
            elif 'AREA' in col:
                area_lists.append(list(mag_table[col]))
            elif 'RAPERT' in col:
                radius_lists.append(list(mag_table[col]))
        
        # Calculate totbkgrds and make list of lists of mnbkrds.
        for mnbkgrd, area_list in zip(mnbkgrds, area_lists):
            totbkgrd_lists.append(np.array(area_list)*mnbkgrd)
            mnbkgrd_lists.append(list(np.ones(len(area_list))*mnbkgrd))
            

    # This portion should be the same for Python and IRAF...

    # Get lists of master IDs for master catalogs.
    master_ids = get_list_masterids(param_dict['targname'])

    # Get a list of master IDs not in this image.
    # Will loop over these seperately and set entries to NaN.
    master_ids_notpresent = list(set(master_ids) - set(master_ids_present))

    imagename = imagename.split('/')[-1]
    print(imagename)
    phot_dict_list = []

    flux_strings = ['flux_2', 'flux_3', 'flux_5', 'flux_7', 'flux_10', 
                    'flux_12', 'flux_15', 'flux_18', 'flux_20', 'flux_24', 
                    'flux_28', 'flux_32', 'flux_36', 'flux_40']
    mnbkgrd_strings = ['mnbkgrd_2', 'mnbkgrd_3', 'mnbkgrd_5', 'mnbkgrd_7', 
                     'mnbkgrd_10', 'mnbkgrd_12', 'mnbkgrd_15', 'mnbkgrd_18', 
                     'mnbkgrd_20', 'mnbkgrd_24', 'mnbkgrd_28', 'mnbkgrd_32', 
                     'mnbkgrd_36', 'mnbkgrd_40']
    totbkgrd_strings = ['totbkgrd_2', 'totbkgrd_3', 'totbkgrd_5', 'totbkgrd_7', 
                     'totbkgrd_10', 'totbkgrd_12', 'totbkgrd_15', 'totbkgrd_18', 
                     'totbkgrd_20', 'totbkgrd_24', 'totbkgrd_28', 'totbkgrd_32', 
                     'totbkgrd_36', 'totbkgrd_40']

    find_ids = npint2int(find_ids)
    print("type find_ids: ", type(find_ids[0]))

    for j in range(len(master_ids_present)):
        phot_dict = {}
        phot_dict['master_id'] = master_ids_present[j]
        phot_dict['find_id'] = find_ids[j]
        phot_dict['imagename'] = imagename
        phot_dict['ingest_date'] = datetime.date.today() 
        phot_dict['xpix'] = xpixs[j]
        phot_dict['ypix'] = ypixs[j]
        phot_dict['ra'] = ras[j]
        phot_dict['dec'] = decs[j]
        for flux_string, flux_list in zip(flux_strings, flux_lists):
            phot_dict[flux_string] = flux_list[j]
        for mnbkgrd_string, mnbkgrd_list in zip(mnbkgrd_strings, mnbkgrd_lists):
            phot_dict[mnbkgrd_string] = mnbkgrd_list[j]      
        for totbkgrd_string, totbkgrd_list in zip(totbkgrd_strings, totbkgrd_lists):
            phot_dict[totbkgrd_string] = totbkgrd_list[j]     

        phot_dict_list.append(phot_dict)


    return phot_dict_list



#-------------------------------------------------------------------------------# 

def return_results_dict(imagename_1, imagename_2, aperture, slope_file, 
    outloc=''):
    """
    Creates dictionary of columns and values for Results table.
    Need row for each image pair, for each photometry aperture that 
    outputted a CTE slope.

    Parameters:
        imagename_1 : string
            Name of the first image in pair.
        imagename_1 : string
            Name of the second image in pair.
        aperture : int
            The aperture of photometry. (3 is nominal)
        slope_file : string
            Text file containing calculated CTE slopes.
        outloc : string
            Location of output.
            
    Returns: 
        results_dict : dictionary
            Dictionary of column names and values of given image pair.
            This will be inserted into database.

    Outputs:
        nothing
    """

    results_dict = {}

    # read in slopefile
    slope_data = ascii.read(slope_file)
    slopes = list(slope_data['col1'])
    slope_stddevs = list(slope_data['col2'])
    num_points= list(slope_data['col3'])

    num_points = npint2int(num_points)


    results_dict['imagename_1'] = imagename_1.split('/')[-1]
    results_dict['imagename_2'] = imagename_2.split('/')[-1]
    results_dict['slopefile'] = slope_file.split('/')[-1]
    results_dict['slopefile_path'] = outloc
    results_dict['aperture'] = aperture
    results_dict['slope_250_500'] = slopes[0]
    results_dict['slopestdev_250_500'] = slope_stddevs[0]
    results_dict['numpoints_250_500'] = num_points[0]
    results_dict['slope_500_1000'] = slopes[1]
    results_dict['slopestdev_500_1000'] = slope_stddevs[1]
    results_dict['numpoints_500_1000'] = num_points[1]
    results_dict['slope_500_2000'] = slopes[2]
    results_dict['slopestdev_500_2000'] = slope_stddevs[2]
    results_dict['numpoints_500_2000'] = num_points[2]
    results_dict['slope_1000_2000'] = slopes[3]
    results_dict['slopestdev_1000_2000'] = slope_stddevs[3]
    results_dict['numpoints_1000_2000'] = num_points[3]
    results_dict['slope_2000_4000'] = slopes[4]
    results_dict['slopestdev_2000_4000'] = slope_stddevs[4]
    results_dict['numpoints_2000_4000'] = num_points[4]
    results_dict['slope_2000_8000'] = slopes[5]
    results_dict['slopestdev_2000_8000'] = slope_stddevs[5]
    results_dict['numpoints_2000_8000'] = num_points[5]
    results_dict['slope_4000_8000'] = slopes[6]
    results_dict['slopestdev_4000_8000'] = slope_stddevs[6]
    results_dict['numpoints_4000_8000'] = num_points[6]
    results_dict['slope_8000_32000'] = slopes[7]
    results_dict['slopestdev_8000_32000'] = slope_stddevs[7]
    results_dict['numpoints_8000_32000'] = num_points[7]

    return results_dict


#-------------------------------------------------------------------------------# 
# Helper functions.
#-------------------------------------------------------------------------------# 

def get_list_masterids(targname):
    """ Returns list of all the master IDs for target.

    Parameters:
        targname : string
            Name of the target cluster.

    Returns:
        master_id : list of ints 

    Outputs:
        nothing
    """
    path_to_mastercats = os.path.join(path_to_data, 'master_cat', 'final')
    master_catalog = os.path.join(path_to_mastercats, 
            '{}_master.cat'.format(targname.lower()))

    data = ascii.read(master_catalog)
    master_id = list(data['col1'])

    return master_id
    

#-------------------------------------------------------------------------------# 

def npint2int(npints):
    """ Changes a list of np.int64 type to plain int type.
    Because it seems sqlite is changing np.intt64 to byte type??
    """

    ints = []
    for npint in npints:
        ints.append(int(npint))

    return ints


#-------------------------------------------------------------------------------# 

def return_parsedname(param_dict):
    """ Depricated!!
    Retained only to easily compare the Python files with the original
    IDL files. The name will only be inserted into the database for 
    verification purposes.
    """
    parsed_imagename = '{}_{}_{}_{}_{}_ci{}_pf{}_cte{}_{}.fits'.format(
                          param_dict['targname'], 
                          param_dict['proposid'], 
                          param_dict['dateobs'], 
                          param_dict['filter'], 
                          param_dict['exptime'], 
                          param_dict['chinject'], 
                          param_dict['flashlvl'], 
                          param_dict['ctecorr'],
                          param_dict['chip'])

    return parsed_imagename

