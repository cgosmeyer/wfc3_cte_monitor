#! /usr/bin/env python

""" The main script for a general run of the WFC3/UVIS External CTE 
monitor's pipeline over four main scenarios:
        postfash=0, pixel-based cte-correction=False
        postfash=0, pixel-based cte-correction=True
        postfash=12, pixel-based cte-correction=False
        postfash=12, pixel-based cte-correction=True
Wraps ``run_image_extraction.py`` (performs photometry and fills in 
FileInfo and Phot tables) and ``run_outputs.py`` (matches pairs to obtain 
CTE slopes, fills in Results tables, and creates analysis plots).

Authors:

    C.M. Gosmeyer, Feb. 2015

Use:
    
    >>> python run_uvis_external_cte.py --pr <program set> --i <'y'/'n'> 
        --ap <list of aperture radii> --pl <list of plots> --ao <'y'/'n'>
        --dophot <'y'/'n'> --irafphot <'y'/'n'> --dofind <'y'/'n'>

    --pr : Which group of proposals to process. Either
        * 'all' (All that are available for the post-flash group).
        * 'last' (Newest of post-flash group available).
        * proposal number (Specific proposal ID. Must be valid 
            for post flash group. Else will return an error).
        * '180' for the 180-degree dataset, visits 10 and 11 in proposal
            12692.

    --v : (optional) Visit numbers to process within selected proposal. 

    --i : (optional) Insert image, mag, and coo file into into FileInfo 
        and Phot tables, and CTE slopes into Results tables?
        True by default. Might turn off if testing stuff.

    --ap : (optional) List of pixel-unit aperture radii (seperated by a 
        space) for which to do photometry and create outputs. By default 
        3 5.
        Note that to make the 'CTE_...' plots for a specific aperture, you 
        need to have created the slope plot ('ratio_ypos') for that 
        aperture concurrently or previously. 

    --pl : (optional) List of plots (seperated by a space) you wish to 
        create. Nomininally all should be done if new data is ingested.
        By default 'ratio_ypos' 'cte_time' 'cte_logflux'. Options also 
        include '180cte_expt' and 'cte_flashlvl'.

    --ao : (optional) Copy outputs to automated_outputs for website? 
        By default 'y'. Might turn this off if doing 180-degree dataset
        or testing stuff.

    --dofind : (optional) Do source finding and matching to Master Catalog? 
        'y'/'n'. By default, 'y'.  Can turn off if already have coo files.

    --dophot : (optional) Do photometry? 'y'/'n'. By default 'y'.

    --irafphot : (optional) Use IRAF.PHOT for photometry? 'y'/'n'.
        By default 'n'.




Examples:

    If ingesting new data of latest epoch do following. Default values
    assume this scenario. 

    >>> python run_uvis_external_cte.py --pr 'last' 

    If want to do photometry and nominal plots only on visits 03 and 04 
    from the latest proposal,

    >>> python run_uvis_external_cte.py --pr 'last' --v 03 04    

    If re-running photometry and databse ingestion over just a specific 
    program (ie, 12379) for aperture 3 and only want to recalculate the
    slope plots,

    >>> python run_uvis_external_cte.py --pr '12379' --ap 3 --pl 'ratio_ypos'

    If want to do all photometry and plots for 180-degree dataset,

    >>> python run_uvis_external_cte.py --pr '180' --pl 'ratio_ypos' '180cte_expt' --ao 'n'

    And so on.    


Outputs:

    All the things...

Notes:

    Translated from the IDL procedure `run_ext_cte.pro`.

"""

from __future__ import print_function

import argparse
import glob
import os
import shutil

from config import path_to_outputs
from config import flashlvl_list
from run_image_extraction import run_image_extraction
from run_outputs import run_outputs
from run_outputs import run_cte_flashlvl_plots
from analysis_tools.fits.sort_fits_by_keyword import sort_fits_by_keyword

#-------------------------------------------------------------------------------# 

def run_uvis_exernal_cte(programset, insert_into_db, visits=[], aprads=[3,5], 
    do_plots=['ratio_ypos', 'cte_time', 'cte_logflux'], 
    do_find=True, do_photom=True, use_iraf=False,
    copy_to_automated_outputs=True, xlims=[55000,57800], ylims=[-0.1,0.6]):
    """
    Main wrapper function of WFC3/UVIS External CTE monitor.

    Parameters:
        programset : int or string
            Which group of proposals to process. Either
            'all' - All that are available for the post-flash group.
            'last' - Newest of post-flash group available).
            <proposal number> - Specific proposal ID. Must be valid 
                for post-flash group, else will return an error.
            '180' - The 180-degree dataset from proposal 12692, visits
                10 and 11.
        insert_into_db : {True, False}
            Insert image, mag, and coo file into into FileInfo and 
            Phot tables, and CTE slopes inot Results tables. True by 
            default. Might turn off if testing stuff.
        visits : list of ints/strings
            Visit numbers to process within selected proposal. 
        aprads : list of ints/floats
            The aperture radii (pixels) on which to perform photometry.
            By default [3, 5].
            Note that to make the 'CTE_...' plots for a specific 
            aperture, you need to have created the slope plot 
            ('ratio_ypos') for that aperture concurrently or previously. 
        do_plots : list of strings
            List of plots you wish to create. Nomininally all defaults 
            should be done if new data is ingested. By default 
            'ratio_ypos' - plots flux ratio vs y-position to obtain 
                a slope, which is the CTE measurement.
            'cte_time' - plots CTE slopes vs MJD time.
            'cte_logflux' - plots CTE slopes vs log of flux. Performs
                polynomial fitting to each epoch to obtain coefficients.
            Less general options include 
            '180cte_expt' - plots CTE slopes of 180-degree data vs 
                exposure time.
            'cte_flashlvl' - plots CTE slopes vs flashlvls.
            If you don't want to run any plots, do ['none']. 
        do_find : {True, False}
            Creates a new *coo file with all coordinates transformed to 
            match the Master Catalog's.  Time consuming.  If just re-running
            photometry and the coo files already exist, can shut this off.
            By default, True.
        do_photom : {True, False}
            Set to True to perform photometry step. True by default.
        use_iraf : {True, False}
            Run photometry with ``iraf.phot`` instead of ``phututils.
            aperture_photometry`` and insert into 
            'uvis_external_cte_iraf.db'.
        copy_to_automated_outputs : {True, False}
            Set to True to copy plots to 'automated_outputs' folder for
            display on the Quicklook website.
            Might turn this off if doing 180-degree dataset or testing stuff.
        xlims : list of floats/ints
            The lower and upper bounds of hte x-axis (MJD time). Update when 
            injesting new epochs!
        ylims : list of floats
            The lower and upper bounds of the y-axis (CTE slope).


    Returns:
        nothing

    Outputs: 
        all the things

    """

    # Loop over the primary values for cte correction and posflash.
    # Perhaps want to make these values user-specified arguments with defaults
    # of below? 
    primary_ctecorr = [False, True]  # Yes/no for cte correction
    primary_flashlvl = flashlvl_list  # In e-

    # If 180-degree dataset, there is no postflash
    if programset == '180':
        primary_flashlvl = [0]

    if do_photom:  
        print("Starting photometry...")
        print(" ")
        for flashlvl in primary_flashlvl:
        # Run on flashlvl=12 first so can check initially whether
        # the proposalset selected is valid for flashlvl datatsets.
            print(">>>Starting flashlvl: {}".format(flashlvl))
            for ctecorr in primary_ctecorr:
                # It is assumed that an FLC exists for each FLT.
                print(">>>Starting ctecorr: {}".format(ctecorr))
                run_image_extraction(flashlvl=flashlvl,
                                     programset=programset,
                                     ctecorr=ctecorr,
                                     overwrite=True,
                                     insert_into_db=insert_into_db,
                                     visits=visits,
                                     use_iraf=use_iraf,
                                     do_find=do_find)
        print(" ")
        print(">>>>>>>>>><<<<<<<<<<<<")
        print("Photometry complete.")
        print(">>>>>>>>>><<<<<<<<<<<<")
        print(" ")
        print(">>>You will find photometry outputs in ")
        print("    outputs/photom/pf0/")
        print("    outputs/photom/pf0_ctecorr/")
        print("    outputs/photom/pf12/")
        print("    outputs/photom/pf12_ctecorr/")
        print(" ")

        # Take all the remaining *coo and *png files and throw them into outputs/noningested.
        path_to_noningested = os.path.join(path_to_outputs, 'noningested')
        noningested_files = glob.glob('*coo') + glob.glob('*png')
        for noningested_file in noningested_files:
            shutil.move(noningested_file, path_to_noningested)

        # Sort the noningested files by targname.
        cwd = os.getcwd()
        os.chdir(path_to_noningested)
        print("Changing directory to {} in order to sort the noningested files.".format(path_to_noningested))
        sort_fits_by_keyword(keyword='targname')
        os.chdir(cwd)
        print("Returning to cwd {}".format(cwd))


    if (bool(set(do_plots) & set(['ratio_ypos', 'cte_time', 'cte_logflux', '180cte_expt']))):
        print("Starting primary outputs generation...")
        for flashlvl in primary_flashlvl:
            print(">>>Starting flashlvl: {}".format(flashlvl))
            for ctecorr in primary_ctecorr:
                print(">>>Starting ctecorr: {}".format(ctecorr))
                run_outputs(flashlvl=flashlvl,
                            programset=programset,
                            ctecorr=ctecorr,
                            aprads=aprads,
                            insert_into_db=insert_into_db,
                            do_plots=do_plots,
                            copy_to_automated_outputs=copy_to_automated_outputs,
                            ylims=ylims,
                            use_iraf=use_iraf)

    if 'cte_flashlvl' in do_plots:
        print("Starting cte vs flashlvl plots...")
        for ctecorr in primary_ctecorr:
            run_cte_flashlvl_plots(ctecorr=ctecorr, 
                                   aprads=aprads, 
                                   copy_to_automated_outputs=False, 
                                   ylims=ylims,
                                   use_iraf=use_iraf)

    print(" ")
    print(">>>>>>>>>>>>><<<<<<<<<<<<<<<")
    print("Output generation complete.")
    print(">>>>>>>>>>>>><<<<<<<<<<<<<<<")
    print(" ")
    print(">>>You will find *xyms files with the photometry outputs listed above.")
    print(">>>You will find output plots and coefficients in ")
    print("    outputs/finalresults/pf0/most_recent/")
    print("    outputs/finalresults/pf0_ctecorr/most_recent/")
    print("    outputs/finalresults/pf12/most_recent/")
    print("    outputs/finalresults/pf12_ctecorr/most_recent/")


#-------------------------------------------------------------------------------# 

def tobool(ans):
    """Takes "Yes"/"No" or "True"/"False" and returns boolean True or False.

    Parameters:
        ans : string
            "Yes"/"No" or "True"/"False".

    Returns:
        boo : {True, False}

    Outputs:
        nothing
    """
    ans = ans.lower()

    if 'y' in ans or 't' in ans:
        boo = True
    elif 'n' in ans or 'f' in ans:
        boo = False
    else:
        boo = None

    return boo


#-------------------------------------------------------------------------------# 
    
def parse_args():
    """Parses command line arguments.
    
    Parameters:
        nothing
        
    Returns:
        args : object
            Containing the image and destination arguments.
            
    Outputs:
        nothing
    """

    programset_help = "Group of programs to process. Either 'all', 'last', "
    programset_help += "or enter program number as string."
    visits_help = "Visit numbers (seperated by space) to process within selected proposal. "
    insert_help = "Insert outputs into database? 'y'/'n'. By default 'y'."  
    insert_help += "Includes photometry and slope outputs."
    aprads_help = "List of pixel-unit aperture radii (seperated by a space) for which "
    aprads_help += "to do photometry and create outputs. By default 3 5."
    do_plots_help = "List of plots (seperated by a space) you wish to create. "
    do_plots_help += "Nomininally all should be done if new data is ingested. "
    do_plots_help += "By default 'ratio_ypos' 'cte_time' 'cte_logflux'. "
    do_plots_help += "Options also include '180cte_expt' and 'cte_flashlvl'. "  
    do_plots_help += "To skip plotting, do 'none'."
    automated_outputs_help = "Copy outputs to automated_outputs for website? 'y'/'n'."
    do_phot_help = "Do photometry? 'y'/'n'. By default 'y'."
    iraf_phot_help = "Do IRAF.PHOT photometry? 'y'/'n'. By default 'n'."
    do_find_help = "Do source finding and matching to Master Catalog? 'y'/'n'."
    do_find_help = "By default, 'y'.  Can turn off if already have coo files."

    parser = argparse.ArgumentParser()
    parser.add_argument('--pr', dest = 'programset',
                        action = 'store', type = str, required = True,
                        help = programset_help)
    parser.add_argument('--v', dest = 'visits',
                        action = 'store', type = str, required = False,
                        help = visits_help,  nargs='+', default=[])
    parser.add_argument('--i', dest = 'insert',
                        action = 'store', type = str, required = False,
                        help = insert_help, default='y')
    parser.add_argument('--ap', dest = 'aprads',
                        action = 'store', type = int, required = False,
                        help = aprads_help, nargs='+', default=[3,5])
    parser.add_argument('--pl', dest = 'do_plots',
                        action = 'store', type = str, required = False,
                        help = do_plots_help, nargs='+', 
                        default=['ratio_ypos', 'cte_time', 'cte_logflux'])
    parser.add_argument('--ao', dest = 'automated_outputs',
                        action = 'store', type = str, required = False,
                        help = automated_outputs_help, default='y')
    parser.add_argument('--dofind', dest = 'do_find',
                        action = 'store', type = str, required = False,
                        help = do_find_help, default='y')
    parser.add_argument('--dophot', dest = 'do_phot',
                        action = 'store', type = str, required = False,
                        help = do_phot_help, default='y')
    parser.add_argument('--irafphot', dest = 'iraf_phot',
                        action = 'store', type = str, required = False,
                        help = iraf_phot_help, default='n')

    args = parser.parse_args()
     
    return args


#-------------------------------------------------------------------------------# 
# The Main. Wraps argsparse. 
#-------------------------------------------------------------------------------# 

if __name__=='__main__':
    args = parse_args()
    programset = args.programset
    visits = args.visits
    insert_into_db = tobool(args.insert)
    aprads = args.aprads
    copy_to_automated_outputs = tobool(args.automated_outputs)
    do_plots = args.do_plots
    do_find = tobool(args.do_find)
    do_photom = tobool(args.do_phot)
    use_iraf = tobool(args.iraf_phot)
    # add flashlvls - list of flashlvls to loop over. config.flashlvl_list by default.
    # add use_ctecorr (or something such) - whether to loop over ctecorr=True,
    #     ctecorr=False, or both.  Both by default.


    # multi-process here?
    run_uvis_exernal_cte(programset=programset, insert_into_db=insert_into_db, 
        visits=visits, aprads=aprads, do_plots=do_plots, do_find=do_find,
        do_photom=do_photom, use_iraf=use_iraf,
        copy_to_automated_outputs=copy_to_automated_outputs)

