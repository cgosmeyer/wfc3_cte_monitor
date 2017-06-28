
""" Module to handle all logic to name paths to outputs in which we want
flashlvl-ctecorrected-timestamp ordered directories. 

Authors:

    C.M. Gosmeyer, Feb. 2016

Outputs:


Notes:
    Translated from the IDL procedure `set_paths_to_outputs.pro`.
"""

import os
from analysis_tools.dir.make_timestamp_dir import make_timestamp_dir

def set_paths_to_outputs(path_to_outputs, basedir, flashlvl, ctecorr, 
    mostrecent, cte_vs_flashlvl=False):
    """ Handles all logic to name paths to outputs in which we want
    flashlvl-ctecorrected-timestamp ordered directories. 

    Parameters: 
	    path_to_outputs : string
	        Path to the 'outputs' directory.
	    basedir : string
	        Main subdirectory to 'outputs' whose subdirectories will be 
	        named by 'flashlvl-ctecorr'.
	    flashlvl : int/string 
	        Whether you want to use the post-flashed data.
	             0 for no post-flash. (Valid for ALL proposals).
	             6 for post-flash 6 e-. (Valid only >= proposal 13083).
	             12 for post-flash 12 e-. ""
	             33 for post-flash 33 e-. ""
	             55 for post-flash 55 e-. ""
	             91 for post-flash 91 e-. ""
	             116 for post-flash 116 e-. ""
	    ctecorr : {True, False}
	        Whether you want to use the CTE-corrected images.
	             False for FLTs.
	             True for FLCs. 
	    mostrecent : {True, False}
	        Set to True if want the directory under 'finalresults'
	        to be 'most_recent/'

	Returns:
        path_to_outputs : string
            Path to the "outputs" directory.

	Outputs:
        Directory named for 'path_to_outputs'.
    """
    
    # Check whether needs be in _ctecorr directory
    if ctecorr:
        if not cte_vs_flashlvl:
    	    pfdir = 'pf{}_ctecorr'.format(flashlvl)
        else: 
            pfdir = 'cte_vs_flashlvl_ctecorr'
    else:
        if not cte_vs_flashlvl:
    	    pfdir = 'pf{}'.format(flashlvl)
        else:
            pfdir = 'cte_vs_flashlvl'

    # Create the path name to outputs
    path_to_outputs = os.path.join(path_to_outputs, basedir, pfdir)

    # Check that the path_to_outputs exists. If not, create.
    if not os.path.isdir(path_to_outputs):
    	os.mkdir(path_to_outputs)

    if (not mostrecent) and ('finalresults' in basedir):
        path_to_outputs = make_timestamp_dir(path_to_outputs)

    elif mostrecent and ('finalresults' in basedir):
    	path_to_outputs = os.path.join(path_to_outputs, 'most_recent')

    return path_to_outputs



