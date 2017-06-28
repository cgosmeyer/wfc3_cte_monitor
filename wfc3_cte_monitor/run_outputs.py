""" 
Module containing functions to generate outputs for the WFC3/UVIS External
CTE monitor, using data prepped in `run_image_extraction.py` and stored in
the FileInfo and Phot tables.

In the functions here the CTE gets measured (in call to `plot_fluxratio_vs_ypos_setup`) 
and filled into the Results table, and using database queries, various plots 
are created for analysis. 

Authors:

    C.M. Gosmeyer, Feb. 2015

Use:

    Nominally execute from the main wrapper ``run_uvis_external_cte.py``.
    But can be run from the command line if user specifies desired 
    parameters in the Main.

    >>> python run_outputs.py

Outputs:

    For each plot specificied in 'do_plots',

    * 'outputs/finalresults/<flashlvl>_<ctecorr>/<time-stamped directory>/<plot_dir>/<plot and related text files>'
    * 'outputs/finalresults/<flashlvl>_<ctecorr>/most_recent/<plot_dir>/<plot and related text files>'
        
    And if 'copy_to_automated_outputs' True,
    * '<automated outputs directory>/<plot>'
    
Notes:

    Translated from the IDL procedure `run_results.pro`.
"""

import glob
import os
import shutil
import config
from uvis_external_cte_plots import plot_fluxratio_vs_ypos_setup
from uvis_external_cte_plots import plot_cteslope_vs_time_setup
from uvis_external_cte_plots import plot_cteslope_vs_logflux_setup
from uvis_external_cte_plots import plot_cteslope_vs_flashlvl_setup
from uvis_external_cte_plots import plot_180_slope_vs_expt_setup
from database_queries import query_for_exptimes
from database_queries import query_for_dateobss
from set_paths_to_outputs import set_paths_to_outputs

# Set globals.
from config import path_to_automated_outputs
from config import path_to_outputs
from config import filter_list
from config import program_list
from config import target_list

from database_interface import get_session
from database_interface import load_connection
from database_interface import NGC104_FileInfo
from database_interface import NGC104_Results
from database_interface import NGC6583_FileInfo
from database_interface import NGC6583_Results
from database_interface import NGC6791_FileInfo
from database_interface import NGC6791_Results



#-------------------------------------------------------------------------------# 

def run_outputs(flashlvl, programset, ctecorr, aprads, insert_into_db=True,
    do_plots=['ratio_ypos', 'cte_time', 'cte_logflux'], 
    copy_to_automated_outputs=True, xlims=[55000,57800], ylims=[-0.1,0.6],
    use_iraf=False):
    """ Wraps creation of all outputs and sorts them into subdirectories
    located in "outputs" directory, and also copies them to the 
    "automated_outputs" directory for display on the Quicklook website if 
    'copy_to_automated_outputs' set to True. Which plots and related text 
    files get created are set in the 'do_plots' parameter. 

    Part of the logic for keeping flashlvl and ctecorr as parameters
    here is that the observations of all flashlvl levels for a given
    target occur in the same orbit. And of course ctecorr is just a
    different file type. 

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
        aprads : list of ints/floats
            The aperture radii (pixels) on which to perform photometry.
        insert_into_db : {True, False}
            Insert CTE slopes inot Results tables. True by default. Might 
            turn off if testing stuff.
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
            If don't want to run any plots, do ['none']. 
        copy_to_automated_outputs : {True, False}
            Set to True to copy plots to 'automated_outputs' folder for
            display on the Quicklook website.
            Might turn this off if doing 180-degree dataset or testing stuff.
        xlims : list of floats
            The lower and upper bounds of the x-axis (MJD time). 
        ylims : list of floats
            The lower and upper bounds of the y-axis (CTE slope). 
        use_iraf : {True, False}
            Run photometry with ``iraf.phot`` instead of ``phututils.
            aperture_photometry`` and insert into 
            'uvis_external_cte_iraf.db' and copy files to *_iraf
            directories in "outputs".

    Returns:
        nothing

    Outputs: 
        For each plot specificied in 'do_plots',

        * 'outputs/finalresults/<flashlvl>_<ctecorr>/<time-stamped directory>/<plot_dir>/<plot and related text files>'
        * 'outputs/finalresults/<flashlvl>_<ctecorr>/most_recent/<plot_dir>/<plot and related text files>'
            
        And if 'copy_to_automated_outputs' True,
        * '<automated outputs directory>/<plot>'
    """
    # Not sure why, but I have to set this here even though I import it above...
    # Perhaps it is because I am modifying the lists? For lists I just 
    # loop over, the above import seems to work fine.
    program_list = config.program_list 
    target_list = config.target_list

    # Check that inputs correct data types.
    flashlvl = str(flashlvl)
    programset = str(programset)

    # Determine proposids from programset, given postflash level.
    # First cut down the list to >= 13083 data if postflashed-only data selected.
    if flashlvl != '0':
        program_list = program_list[4:]

    # Second cut down by selected set.
    if programset == 'all':
        program_list = program_list
    elif programset == 'last':
        program_list = [program_list[-1]]
    elif programset == '180':
        program_list = ['12692']
        target_list = [target_list[2]]
    else:
        if programset not in program_list:
            print("Invalid proposal ID of {} for a post-flash set. (Must be >= 13083).".format(programset))
            print("Valid propals are {}".format(program_list))
            print("If you do indeed want a program set < 13083, set 'flashlvl=0'.")
            return "Invalid post-flash level for proposid."
        else:
            program_list = [programset]
    print("program_list: {}".format(program_list))

    # Create paths to the outputs for this postflash-ctecorr combination. 
    path_to_workingdir = os.getcwd()       
    print("path_to_workingdir: {}".format(path_to_workingdir))
    path_to_outputs = config.path_to_outputs

    print("path_to_outputs: {}".format(path_to_outputs))
    
    if use_iraf:
        basedir = 'finalresults_iraf'
        copy_to_automated_outputs = False
    else:
        basedir = 'finalresults'

    path_to_finalresults = set_paths_to_outputs(path_to_outputs=path_to_outputs, 
                                                basedir=basedir, 
                                                flashlvl=flashlvl, 
                                                ctecorr=ctecorr, 
                                                mostrecent=False)
    print("path_to_finalresults: {}".format(path_to_finalresults))

    path_to_mostrecent = set_paths_to_outputs(path_to_outputs=path_to_outputs, 
                                                basedir=basedir, 
                                                flashlvl=flashlvl, 
                                                ctecorr=ctecorr, 
                                                mostrecent=True)
    print("path_to_mostrecent: {}".format(path_to_mostrecent))
    if not os.path.isdir(path_to_mostrecent):
        os.mkdir(path_to_mostrecent)
        print("Created {}".format(path_to_mostrecent))

    # ***********
    # 1. Create the flux-ratio vs y-position plots and obtain slopes.
    #    The slopes will be recorded in Results tables in the database
    #    and used in all subsequent plots.
    # ***********
    if 'ratio_ypos' in do_plots:
        # Loop over all valid proposals.
        for proposid in program_list:
            print("proposid: {}".format(proposid))
            # Loop over all targets. If query turns up nothing for target, should quietly pass.
            for targname in target_list:
                print("targname: {}".format(targname))
                # Loop over all filters. If query turns up nothing for filter, should quietly pass.
                for filt in filter_list:
                    print("filter: {}".format(filt))
                    # determine exptimes for targname, proposid, and filter
                    exptimes = query_for_exptimes(targname, proposid, filt)
                    exptimes = list(set(exptimes))
                    print("exptimes found: {}".format(exptimes))
                    for exptime in exptimes:
                        print("exptime: {}".format(exptime))
                        # determine all dateobs for given targname, proposid, filter, and exptime
                        dateobss = query_for_dateobss(targname, proposid, filt, exptime)
                        dateobss = list(set(dateobss))
                        print("dateobss found: {}".format(dateobss))
                        for dateobs in dateobss:
                            print("dateobs: {}".format(dateobs))
                            # Finally! loop over given apertures
                            for aperture in aprads:
                                print("aperture: {}".format(aperture))

                                # Come up with final output location here so can record in db.
                                # This is for sorting into directories.
                                if exptime <= 60.:
                                    length = 's'
                                elif exptime > 60.:
                                    length = 'l'

                                photfile_dir = '{}_{}_{}'.format(targname, filt, length)
                                print("photfile_dir: {}".format(photfile_dir))
                                path_to_photdir = os.path.join(path_to_finalresults, photfile_dir)
                                print("path_to_photdir: {}".format(path_to_photdir))

                                # This will create slopes for the given dataset that will be used in the
                                # next plots.
                                # It also inserts slopes into Results table.
                                if programset == '180':
                                    for chip in [1,2]:
                                        plot_fluxratio_vs_ypos_setup(
                                                               targname=targname, 
                                                               proposid=proposid, 
                                                               dateobs=dateobs, 
                                                               filt=filt, 
                                                               exptime=exptime, 
                                                               chinject='NO', 
                                                               flashlvl=flashlvl,
                                                               ctecorr=ctecorr, 
                                                               chip=chip,
                                                               postarg1=0, 
                                                               aperture=aperture,
                                                               outloc=path_to_photdir,
                                                               write_ratios_to_file=False, 
                                                               insert_into_db=insert_into_db)
                                else:
                                    plot_fluxratio_vs_ypos_setup(
                                       targname=targname, 
                                       proposid=proposid, 
                                       dateobs=dateobs, 
                                       filt=filt, 
                                       exptime=exptime, 
                                       chinject='NO', 
                                       flashlvl=flashlvl, 
                                       ctecorr=ctecorr, 
                                       postarg1=0, 
                                       aperture=aperture,
                                       outloc=path_to_photdir,
                                       write_ratios_to_file=False, 
                                       insert_into_db=insert_into_db)

                        if dateobss != []:
                            # Move files in the exptime loop.
                            # Create directories that don't exist.
                            if not os.path.isdir(path_to_photdir):
                                os.mkdir(path_to_photdir)
                             
                            # Create a 'most_recent/photfile_dir' directory if it doesn't exist.
                            path_to_mostrecent_photdir = os.path.join(path_to_mostrecent, photfile_dir)
                            if not os.path.isdir(path_to_mostrecent_photdir):
                                os.mkdir(path_to_mostrecent_photdir)

                            # Move and overwrite all aperture-specific files to the photfile directory.
                            # But first copy everything to the 'most_recent/photfile_dir' dir. 
                            # glob the cwd for the *slopes.txt, *slopes.png, move to same dir
                            slope_files = glob.glob('*slopes.txt') + glob.glob('*slopes.png')
                            for slope_file in slope_files:
                                shutil.copy(slope_file, path_to_mostrecent_photdir)
                                shutil.move(slope_file, path_to_photdir)


    # ***********
    # 2. Create the slope vs MJD plots.
    # ***********
    if 'cte_time' in do_plots and programset != '180':    
        for exp_length in ['l', 's']:
            for filt in filter_list:
                plot_cteslope_vs_time_setup(
                    exp_length=exp_length, 
                    filt=filt, 
                    flashlvl=flashlvl, 
                    ctecorr=ctecorr, 
                    aperture=3,
                    chinject='NO', 
                    postarg1=0, 
                    xlims=xlims,
                    ylims=ylims,
                    write_to_file=True)

        sort_outputs(plot_dir='cte_vs_time', 
            plot_name='cteVStime', 
            extension='png', 
            path_to_finalresults=path_to_finalresults,
            path_to_mostrecent=path_to_mostrecent,
            copy_to_automated_outputs=copy_to_automated_outputs)

        sort_outputs(plot_dir='cte_vs_time', 
            plot_name='cteVStime', 
            extension='txt', 
            path_to_finalresults=path_to_finalresults,
            path_to_mostrecent=path_to_mostrecent,
            copy_to_automated_outputs=False)


    # ***********
    # 3. Create the slope vs log Flux plots, and fit the empirical model.
    # ***********
    if 'cte_logflux' in do_plots and programset != '180': 
        for exp_length in ['l', 's']:
            for filt in filter_list: 
                plot_cteslope_vs_logflux_setup(
                    targnames=target_list,
                    exp_length=exp_length, 
                    filt=filt, 
                    flashlvl=flashlvl, 
                    ctecorr=ctecorr, 
                    aperture=3, 
                    chinject='NO', 
                    postarg1=0, 
                    ylims=ylims)

        sort_outputs(plot_dir='cte_vs_logflux', plot_name='cteVSlogflux', 
            extension='png',
            path_to_finalresults=path_to_finalresults, 
            path_to_mostrecent=path_to_mostrecent,
            copy_to_automated_outputs=copy_to_automated_outputs)

        sort_outputs(plot_dir='cte_vs_logflux', plot_name='cteVSlogflux', 
            extension='txt',
            path_to_finalresults=path_to_finalresults, 
            path_to_mostrecent=path_to_mostrecent,
            copy_to_automated_outputs=False)

    elif 'cte_logflux' in do_plots and programset == '180': 
        for chip in [1,2]:
            for exp_length in ['l', 's']:
                for filt in filter_list: 
                    plot_cteslope_vs_logflux_setup(
                        targnames=target_list,
                        exp_length=exp_length, 
                        filt=filt, 
                        flashlvl=flashlvl, 
                        ctecorr=ctecorr, 
                        aperture=3, 
                        chinject='NO', 
                        postarg1=0, 
                        chip=chip,
                        ylims=ylims)

        sort_outputs(plot_dir='cte_vs_logflux', plot_name='cteVSlogflux', 
            extension='png',
            path_to_finalresults=path_to_finalresults, 
            path_to_mostrecent=path_to_mostrecent,
            copy_to_automated_outputs=False)

        sort_outputs(plot_dir='cte_vs_logflux', plot_name='cteVSlogflux', 
            extension='txt',
            path_to_finalresults=path_to_finalresults, 
            path_to_mostrecent=path_to_mostrecent,
            copy_to_automated_outputs=False)


    # ***********
    # 4. If 180 degree data...
    # ***********
    if '180cte_expt' in do_plots and programset == '180': 
        plot_180_slope_vs_expt_setup(ctecorr=ctecorr, aprads=aprads,
            chip1_invert=True)

        sort_outputs(plot_dir='180cte_expt', plot_name='180cteVSexpt', 
            extension='png',
            path_to_finalresults=path_to_finalresults, 
            path_to_mostrecent=path_to_mostrecent,
            copy_to_automated_outputs=False)


#-------------------------------------------------------------------------------# 

def run_cte_flashlvl_plots(ctecorr, aprads=[3], 
    copy_to_automated_outputs=True, xlims=[-1,120], ylims=[-0.1,0.6], 
    use_iraf=False):
    """ Wraps creater of CTE slope vs flashlvl plots.

    Parameters:
        ctecorr : {True, False}
            True for FLCs (pixel-based CTE-corrected)
            False for FLTs (not pixel-based CTE-corrected)
        aprads : list of ints/floats
            The aperture radii (pixels) on which to perform photometry.
        copy_to_automated_outputs : {True, False}
            Set to True to copy plots to 'automated_outputs' folder for
            display on the Quicklook website.
            Might turn this off if doing 180-degree dataset or testing stuff.
        xlims : list of floats/ints
            The lower and upper bounds of the x-axis (flashlvl).
        ylims : list of floats
            The lower and upper bounds of the y-axis (CTE slope).

    Returns:
        nothing

    Outputs:
        For each plot specificied in 'do_plots',

        * 'outputs/finalresults/<flashlvl>_<ctecorr>/<time-stamped directory>/<plot_dir>/<plot and related text files>'
        * 'outputs/finalresults/<flashlvl>_<ctecorr>/most_recent/<plot_dir>/<plot and related text files>'
            
        And if 'copy_to_automated_outputs' True,
        * '<automated outputs directory>/<plot>'

    Notes:
        I can't remember why this is here and not in `uvis_external_cte_plots.py`.
        Might consider moving it. 

    """  
    if use_iraf:
        basedir = 'finalresults_iraf'
    else:
        basedir = 'finalresults'
    for exp_length in ['l', 's']:
        for filt in filter_list:
            for aperture in aprads:
                plot_cteslope_vs_flashlvl_setup(
                    exp_length=exp_length, 
                    filt=filt, 
                    ctecorr=ctecorr, 
                    aperture=aperture, 
                    chinject='NO', 
                    postarg1=0, 
                    outloc='', 
                    xlims=xlims, 
                    ylims=[-0.1,0.4])

    # And really only should copy to automate_outputs the long expsore F502N ngc104
    path_to_finalresults = set_paths_to_outputs(path_to_outputs=path_to_outputs, 
                                                basedir=basedir, 
                                                flashlvl='', 
                                                ctecorr=ctecorr, 
                                                mostrecent=False,
                                                cte_vs_flashlvl=True)
    print("path_to_finalresults: {}".format(path_to_finalresults))

    path_to_mostrecent = set_paths_to_outputs(path_to_outputs=path_to_outputs, 
                                                basedir=basedir, 
                                                flashlvl='', 
                                                ctecorr=ctecorr, 
                                                mostrecent=True,
                                                cte_vs_flashlvl=True)
    print("path_to_mostrecent: {}".format(path_to_mostrecent))
    if not os.path.isdir(path_to_mostrecent):
        os.mkdir(path_to_mostrecent)
        print("Created {}".format(path_to_mostrecent))

    sort_outputs(plot_dir='', plot_name='cteVSflashlvl', 
        extension='png',
        path_to_finalresults=path_to_finalresults,
        path_to_mostrecent=path_to_mostrecent,
        copy_to_automated_outputs=copy_to_automated_outputs)


#-------------------------------------------------------------------------------# 
# Helper functions.
#-------------------------------------------------------------------------------# 

def sort_outputs(plot_dir, plot_name, extension, path_to_finalresults, 
    path_to_mostrecent, copy_to_automated_outputs):
    """ Creates time-stamped directories for outputs and sorts the 
    outputs into the directories. Also copies outputs to 'automated_outputs'
    for display on the Quicklook website.

    Parameters:
        plot_dir : string
            Name of the plot's primary directory, ie, "cte_vs_time"
        plot_name : string
            Basename of the plot, ie, "cteVStime".
        extension : string
            Extension of the file, ie, "png" or "txt".
        path_to_finalresults : string
            Path to the "finalresults" directory.
        path_to_mostrecent :
            Path to the "mostrecent" directory.
        copy_to_automated_outputs : {True, False}
            Set to True to copy plots to 'automated_outputs' folder for
            display on the Quicklook website.
            Might turn this off if doing 180-degree dataset or testing stuff.

    Returns:
        nothing

    Outputs:
        Nominally,
        * '<time-stamped directory>/<plot_dir>/<plot_name>.<extension>'
        * 'most_recent/<plot_dir>/<plot_name>.<extension>'
        If 'copy_to_automated_outputs' True,
        * '<automated outputs directory>/<plot_name>.<extension>'
    """
    # Create the time-stamped directory/plot_dir/
    path_to_plotdir = os.path.join(path_to_finalresults, plot_dir)
    print("path_to_plotdir for {}: {}".format(plot_dir, path_to_plotdir))
    if not os.path.isdir(path_to_plotdir):
        os.mkdir(path_to_plotdir)   

    # Create a 'most_recent/plot_dir' directory if it doesn't exist.
    path_to_mostrecent_plotdir = os.path.join(path_to_mostrecent, plot_dir)
    if not os.path.isdir(path_to_mostrecent_plotdir):
        os.mkdir(path_to_mostrecent_plotdir)

    slope_plots = glob.glob('{}*{}'.format(plot_name, extension))
    print("slope_plots for {}: ".format(plot_name, slope_plots))
    # Copy plots to 'most_recent/plot_dir' and automated_outputs.
    # Move them to their final locoation in 'final_results/timestampdir/plot_dir'.
    for slope_plot in slope_plots:
        shutil.copy(slope_plot, path_to_mostrecent_plotdir)
        if copy_to_automated_outputs:
            if plot_name == 'cteVSflashlvl':
                # if cteVSflaslvl, only copy the ngc104 and ngc6791, F502N, 
                # long exposure plots
                if 'ngc104' in slope_plot and '_l_' in slope_plot and \
                'F502N' in slope_plot:
                    shutil.copy(slope_plot, path_to_automated_outputs)
                    print("Copying {} to automated_outputs/daily_outputs."\
                        .format(slope_plot))
            else: 
                shutil.copy(slope_plot, path_to_automated_outputs)
                print("Copying {} to automated_outputs/daily_outputs."\
                        .format(slope_plot))
        shutil.move(slope_plot, path_to_plotdir)


#-------------------------------------------------------------------------------# 
# The Main. Use for testing.
#-------------------------------------------------------------------------------# 

if __name__=="__main__":
    flashlvl='0'
    programset='14378'
    ctecorr=True
    aprads=[3]
    #run_outputs(flashlvl, programset, ctecorr, aprads)

    exp_length = 's'
    filt = 'F502N'
    flashlvl = 0
    ctecorr=True
    plot_cteslope_vs_time_setup(exp_length, filt, flashlvl, ctecorr, aperture=3,
    chinject='NO', postarg1=0)
