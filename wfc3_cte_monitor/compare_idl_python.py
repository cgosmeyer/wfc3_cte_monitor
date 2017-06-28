#! /usr/bin/env python

"""
One-off script whose purpose is to compare the IDL pipeline's slopes to 
the Python pipeline's. 

Creates two plots:

1. Plot IDL-derived slopes on Python-derived slopes vs log10 flux bin for 
   proposal 14012 (the last proposal to be run through IDL pipeline)

2. Fraction of sources recovered vs flux bin.

Author:

    C.M. Gosmeyer, Nov. 2016

"""

from astropy.io import ascii
import config
import itertools
import numpy as np
import os
import pylab

from set_paths_to_outputs import set_paths_to_outputs
from database_queries import query_for_pair, query_results_for_slopes, query_for_dateobss


# -----------------------------------------------------------------------------

def read_idl_slope_file(targname, filt, exp_length, flashlvl, ctecorr, aperture):
    """

    outputs_old/finalresults/pf<>/most_recent/dmag_y/<target>_<filter>_<l/s>_r3_<binlo>_<binhi>.res 
    Each *res file contains all epochs for a single flux bin.
    Whereas the python files contain all flux bins for all epochs. 

    Parameters:
        targname : string
            Name of the target cluster.
        filt : string
            Name of the filter.
        exp_length : string
            Length of the exposure. Either 'l' or 's'.
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
        ctecorr : {True, False}
            True for FLCs (pixel-based CTE-corrected)
            False for FLTs (not pixel-based CTE-corrected)
        aperture : string or int
            The radius of the aperture in pixels.


    Returns:
        slopes_master : numpy array
            Slopes of given parameters.
        stderrs_master : numpy array
            Slopes' standard errors.
        epochs_master : numpy array
            Epochs of the slopes.
        fluxbins_master : numpy_array
            Fluxbins of available data.
    """
    path_to_outputs = config.path_to_outputs[:-1] + '_old'
    path_to_finalresults = os.path.join(path_to_outputs, 'finalresults')

    if ctecorr:
        ctecorr_str = '_ctecorr'
    else:
        ctecorr_str = ''
    targ_str = (targname.lower()).split('ngc')[1]

    # Get path.
    path_to_mostrecent = os.path.join(path_to_finalresults, 
        'pf{}{}'.format(flashlvl, ctecorr_str), 
        'most_recent',
        'dmag_y')

    print("IDL path_to_mostrecent: {}".format(path_to_mostrecent))

    # Initialize master lists
    mjds_master = []
    slopes_master = []
    stderrs_master = []
    fluxbins_master = []

    # Loop through the flux bins.
    for binlo, binhi in zip(config.fluxbins_lo, config.fluxbins_hi):
        # Create filename
        file_name = 'n{}_{}_{}_r{}_{}_{}.res'.format(targ_str, filt, 
            exp_length, aperture, binlo, binhi)

        print('file_name: {}'.format(file_name))

        path_to_file = os.path.join(path_to_mostrecent, file_name)
        print("path_to_file: {}".format(path_to_file))
        
        if os.path.isfile(path_to_file):
            data = ascii.read(path_to_file)
            print("data: ", data)

            mjds = np.array(data['col1'])
            slopes = np.array(data['col5'])
            stderrs = np.array(data['col6'])
            
            for mjd, slope, stderr in zip(mjds, slopes, stderrs):
                mjds_master.append(mjd)
                slopes_master.append(slope)
                stderrs_master.append(stderr)
                fluxbins_master.append('{}-{}'.format(binlo, binhi))

    print("idl slopes: {}".format(slopes_master))
    print("idl stderrs: {}".format(stderrs_master))
    print("idl mjds: {}".format(mjds_master))
    print("idl fluxbins: {}".format(fluxbins_master))

    return np.array(slopes_master), np.array(stderrs_master), np.array(mjds_master), np.array(fluxbins_master)

# -----------------------------------------------------------------------------

def query_db_python(targname, filt, exp_length, flashlvl, ctecorr, aperture, 
    epoch, proposid='14012'):
    """ 
    Alternative to :func:'read_python_slope_file'

    Parameters:
        targname : string
            Name of the target cluster.
        filt : string
            Name of the filter.
        exp_length : string
            Length of the exposure. Either 'l' or 's'.
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
        ctecorr : {True, False}
            True for FLCs (pixel-based CTE-corrected)
            False for FLTs (not pixel-based CTE-corrected)
        aperture : string or int
            The radius of the aperture in pixels.

    """
    # Initialize master lists
    mjds_master = []
    slopes_master = []
    stderrs_master = []
    fluxbins_master = []
    numpoints_master = []

    if exp_length == 'l':
        if '104' in targname:
            exptime = 348
        elif '6791' in targname:
            exptime = 420
    elif exp_length == 's':
        if '104' in targname:
            exptime = 30
        elif '6791' in targname:
            exptime = 60

    # Loop through the flux bins.
    for binlo, binhi in zip(config.fluxbins_lo, config.fluxbins_hi):
        fluxbin = '{}-{}'.format(binlo, binhi)
        # First need to query for all imagenames 
        imagenames_chip1, imagenames_chip2 = query_for_pair(
            targname=targname, proposid=proposid, dateobs=epoch, filt=filt, exptime=exptime, 
            chinject='NO', flashlvl=flashlvl, ctecorr=ctecorr, postarg1=0)

        # Now query for slopes.
        for imagename_chip1, imagename_chip2 in zip(imagenames_chip1, imagenames_chip2):
            slopes, slopesdtdevs, numpoints = query_results_for_slopes(
                targname=targname, imagename_chip1=imagename_chip1, 
                imagename_chip2=imagename_chip2, 
                fluxbin=fluxbin, aperture=aperture)
            # Don't forget to correct the slopes! 
            slopes_master.append((slopes[0]/2.)*2048.)
            stderrs_master.append(((slopesdtdevs[0]/2.)*2048) / np.sqrt(numpoints[0]))
            mjds_master.append(epoch)
            fluxbins_master.append(fluxbin)
            numpoints_master.append(numpoints[0])
    
    print("python slopes: {}".format(slopes_master))
    print("python stderrs: {}".format(stderrs_master))
    print("python mjds: {}".format(mjds_master))
    print("python fluxbins: {}".format(fluxbins_master))


    return np.array(slopes_master), np.array(stderrs_master), \
    np.array(mjds_master), np.array(fluxbins_master), np.array(numpoints_master)


# -----------------------------------------------------------------------------

def read_python_slope_file(targname, filt, exp_length, flashlvl, ctecorr, aperture):
    """ Reads in slope file of Python-derived slopes over all epochs.

    Going to match to the IDL using the epochs. 

    Parameters:
        targname : string
            Name of the target cluster.
        filt : string
            Name of the filter.
        exp_length : string
            Length of the exposure. Either 'l' or 's'.
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
        ctecorr : {True, False}
            True for FLCs (pixel-based CTE-corrected)
            False for FLTs (not pixel-based CTE-corrected)
        aperture : string or int
            The radius of the aperture in pixels.

    Returns:
        slopes_cut : numpy array
            Slopes of given parameters, including only those with given 
            targname and fluxbin.
        stderrs_cut : numpy array
            Slopes' standard errors, including only those with given 
            targname and fluxbin.
        epochs_cut : numpy array
            Epochs of the slopes, including only those with given 
            targname and fluxbin.

    Outputs: 
        nothing

    """
    # Get path.
    path_to_mostrecent = set_paths_to_outputs(
        config.path_to_outputs, 
        basedir='finalresults_iraf', 
        flashlvl=flashlvl, 
        ctecorr=ctecorr, 
        mostrecent=True)

    path_to_mostrecent = os.path.join(path_to_mostrecent, 'cte_vs_time')

    print('path_to_mostrecent: {}'.format(path_to_mostrecent))

    # Create name of file.
    file_name = 'cteVStime_{}_{}_pf{}_ctecorr{}_r{}.txt'\
    .format(filt, exp_length, flashlvl, ctecorr, aperture)

    print('file_name: {}'.format(file_name))

    path_to_file = os.path.join(path_to_mostrecent, file_name)

    # Open file if it exists and read slopes
    if os.path.isfile(path_to_file):
        data = ascii.read(path_to_file)
        print(data)
        slopes = np.array(data['slope'])
        stderrs = np.array(data['slope_stderr'])
        mjds = np.array(data['mjd'])
        targnames = np.array(data['targname'])
        fluxbins = np.array(data['fluxbin'])

        # Save only files with wanted targname and flux-bin.
        slopes_cut = np.array(slopes[np.where(targnames == targname)])
        stderrs_cut = np.array(stderrs[np.where(targnames == targname)])
        mjds_cut = np.array(mjds[np.where(targnames == targname)])   
        fluxbins_cut = np.array(fluxbins[np.where(targnames == targname)])     

        print("python slopes: {}".format(slopes_cut))
        print("python stderrs: {}".format(stderrs_cut))
        print("python mjds: {}".format(mjds_cut))
        print("python fluxbins: {}".format(fluxbins_cut))

        return slopes_cut, stderrs_cut, mjds_cut, fluxbins_cut

    else:
        print("File does not exist: {}".format(path_to_file)) 
        return [], [], [], []


# -----------------------------------------------------------------------------

def plot_numpoints(targname, filt, exp_length, flashlvl, aperture,
                      outloc, ylims):
    """ Plots fraction of non-ctecorr to ctecorr sources recovered 
    vs fluxbin, for a subset of epochs.
    """

    proposid_ls = ['12379', '12692', '13083', '13566', '14012']
    color_ls = ['red', 'orange', 'lime', 'green', 'blue'] 
    marker_ls = ['o', 'd', 's', '*', '^']
    #color_ls = itertools.cycle(colors)

    epoch_ls = []
    #proposid_ls = []
    for proposid in proposid_ls:
        if exp_length == 'l':
            if '104' in targname:
                exptime = 348
            elif '6791' in targname:
                exptime = 420
        elif exp_length == 's':
            if '104' in targname:
                exptime = 30
            elif '6791' in targname:
                exptime = 60
        epochs = query_for_dateobss(targname, proposid, filt, exptime)
        epochs = list(set(epochs))    
        epoch_ls.append(epochs[-1])
        #for epoch in epochs:
        #    proposid_ls.append(proposid)
    print("epoch_ls: {}".format(epoch_ls))
    
    # Set the figure.
    pylab.figure(figsize=(12.5,8.5))

    for epoch, proposid, color, marker in zip(epoch_ls, proposid_ls, color_ls, marker_ls):
        # ctecorr True
        slopes_ctecorr, stderrs_ctecorr, mjds_ctecorr, fluxbins_ctecorr, numpoints_ctecorr = query_db_python(
            targname, filt, exp_length, flashlvl, ctecorr=True, aperture=aperture, 
            epoch=epoch, proposid=proposid)

        mjds_ctecorr_cut = np.array(mjds_ctecorr[np.where(mjds_ctecorr == epoch)])
        fluxbins_ctecorr_cut = np.array(fluxbins_ctecorr[np.where(mjds_ctecorr == epoch)])
        numpoints_ctecorr_cut = np.array(numpoints_ctecorr[np.where(mjds_ctecorr == epoch)])

        # ctecorr False
        slopes, stderrs, mjds, fluxbins, numpoints, = query_db_python(
            targname, filt, exp_length, flashlvl, ctecorr=False, aperture=aperture, 
            epoch=epoch, proposid=proposid)

        mjds_cut = np.array(mjds[np.where(mjds == epoch)])
        fluxbins_cut = np.array(fluxbins[np.where(mjds == epoch)])
        numpoints_cut = np.array(numpoints[np.where(mjds == epoch)])

        print("fluxbins_cut: {}".format(fluxbins_cut))   
        print("numpoints_cut: {}".format(numpoints_cut))    
        print("numpoints_ctecorr_cut: {}".format(numpoints_ctecorr_cut))   

        if len(numpoints_cut) == len(numpoints_ctecorr_cut) and len(numpoints_cut) != 0:

            fluxbins_log = []
            for fluxbin in fluxbins_cut:
                fluxlo = float(fluxbin.split('-')[0])
                fluxhi = float(fluxbin.split('-')[1])
                flux_av = (fluxlo + fluxhi)/2.0
                fluxbins_log.append(np.log10(flux_av))

            # Calculate fraction recovered

            frac_recovered =(1.0 - (numpoints_ctecorr_cut.astype(float) - numpoints_cut.astype(float))/ numpoints_ctecorr_cut.astype(float))*100.

            print("frac_recovered: {}".format(frac_recovered))     

            # Set the next color in sequence.
            #color = next(color_ls)
            pylab.scatter(fluxbins_log, frac_recovered, 
                    s=120, marker=marker, color='grey', alpha=0.6, label='MJD={}'.format(epoch))

    pylab.xlabel('LOG10 Flux [e-])', fontsize=22, weight='bold') 
    pylab.ylabel('% Sources Recovered w/ CTEcorr', fontsize=22, weight='bold')
    pylab.axhline(y=100.0, linewidth=2, linestyle='--', color='grey')
    pylab.tick_params(axis='both', which='major', labelsize=20)
    title = "{} {} explen'{}' pf{} ap{}".format(targname, filt, 
        exp_length, flashlvl, aperture)
    pylab.title(title, fontsize=16)
    pylab.xlim([2.5, 4.5])
    pylab.ylim([40.,120])
    pylab.legend(scatterpoints=1, loc='lower right')
    pylab.savefig(os.path.join(outloc, '{}_{}_{}_pf{}_r{}_fracrecoverd.png'.format(
        targname, filt, exp_length, flashlvl, aperture)), 
        bbox_inches='tight')



# -----------------------------------------------------------------------------

def plot_slopes(targname, filt, exp_length, flashlvl, aperture,
                      outloc, ylims):
    """
    """
    if targname == 'ngc104':
        epoch = 57226.
    elif targname == 'ngc6791':
        epoch = 57221.
    
    # Set the figure.
    pylab.figure(figsize=(12.5,8.5))

    for ctecorr in [True, False]: 
        slopes_idl, stderrs_idl, mjds_idl, fluxbins_idl = read_idl_slope_file(
            targname, filt, exp_length, flashlvl, ctecorr, aperture)

        slopes_python, stderrs_python, mjds_python, fluxbins_python, numpoints_python = query_db_python(
            targname, filt, exp_length, flashlvl, ctecorr, aperture, 
            epoch, proposid='14012')

        #slopes_python, stderrs_python, mjds_python, fluxbins_python = read_python_slope_file(
        #    targname, filt, exp_length, flashlvl, ctecorr, aperture)

        
        mjds_idl_cut = np.array(mjds_idl[np.where(mjds_idl == epoch)])
        mjds_python_cut = np.array(mjds_python[np.where(mjds_python == epoch)])
        fluxbins_idl_cut = np.array(fluxbins_idl[np.where(mjds_idl == epoch)])
        fluxbins_python_cut = np.array(fluxbins_python[np.where(mjds_python == epoch)])
        slopes_idl_cut = np.array(slopes_idl[np.where(mjds_idl == epoch)])
        slopes_python_cut = np.array(slopes_python[np.where(mjds_python == epoch)])
        stderrs_idl_cut = np.array(stderrs_idl[np.where(mjds_idl == epoch)])
        stderrs_python_cut = np.array(stderrs_python[np.where(mjds_python == epoch)])

        print("mjds_idl_cut: {}".format(mjds_idl_cut))
        print("mjds_python_cut: {}".format(mjds_python_cut))
        print("fluxbins_idl_cut: {}".format(fluxbins_idl_cut))
        print("fluxbins_python_cut: {}".format(fluxbins_python_cut))
        print("slopes_idl_cut: {}".format(slopes_idl_cut))
        print("slopes_python_cut: {}".format(slopes_python_cut))    

        # Transform flux bins to floats, take average, and take log10, so can plot.
        fluxbins_idl_log = []
        for fluxbin in fluxbins_idl_cut:
            fluxlo = float(fluxbin.split('-')[0])
            fluxhi = float(fluxbin.split('-')[1])
            flux_av = (fluxlo + fluxhi)/2.0
            fluxbins_idl_log.append(np.log10(flux_av))

        fluxbins_python_log = []
        for fluxbin in fluxbins_python_cut:
            fluxlo = float(fluxbin.split('-')[0])
            fluxhi = float(fluxbin.split('-')[1])
            flux_av = (fluxlo + fluxhi)/2.0
            fluxbins_python_log.append(np.log10(flux_av))

        # Choose marker color.
        if ctecorr:
            idl_color = 'blue'
            python_color = 'red'
        else:
            idl_color = 'skyblue'
            python_color = 'orange'

        pylab.scatter(fluxbins_idl_log, slopes_idl_cut, 
                s=90, marker='d', color=idl_color, alpha=0.6, label='idl, ctecorr={}'.format(ctecorr))
        pylab.errorbar(fluxbins_idl_log, slopes_idl_cut,  
                    yerr=stderrs_idl_cut, 
                    color='gray', alpha=0.75, marker=None, ls='None')

        pylab.scatter(fluxbins_python_log, slopes_python_cut, 
                s=90, marker='o', color=python_color, alpha=0.6, label='python, ctecorr={}'.format(ctecorr))
        pylab.errorbar(fluxbins_python_log, slopes_python_cut,  
                    yerr=stderrs_python_cut, 
                    color='gray', alpha=0.75, marker=None, ls='None')

    # Complete the plot and save. 

    pylab.xlabel('LOG10 Flux [e-]', fontsize=22, weight='bold') 
    pylab.ylabel('CTE loss [flux / 2048 pxl]', fontsize=22, weight='bold')
    pylab.tick_params(axis='both', which='major', labelsize=20)
    for ymaj in np.arange(11)*0.1:
        pylab.axhline(y=ymaj,linestyle='--',color='grey')
    title = "{} {} explen'{}' pf{} ap{}".format(targname, filt, 
        exp_length, flashlvl, aperture)
    pylab.title(title, fontsize=16)
    pylab.annotate('MJD={}'.format(epoch), xy=(3.35, 0.55), fontsize=20, weight='bold') #1.1
    pylab.xlim([2.5, 4.5])
    pylab.ylim(ylims)
    pylab.legend(scatterpoints=1, loc='upper right')
    pylab.savefig(os.path.join(outloc, '{}_{}_{}_pf{}_r{}_idlvspy.png'.format(
        targname, filt, exp_length, flashlvl, aperture)), 
        bbox_inches='tight')


# -----------------------------------------------------------------------------
# The Main.
# -----------------------------------------------------------------------------

def compare_idl_python():
    """ Loops over combinations of targname, filt, exp_length, flashlvl, 
    and ctecorr.

    """
    outloc_compare = os.path.join(config.path_to_outputs, 'finalresults_iraf', 'compare_idl_python_pipelines')
    outloc_frac = os.path.join(config.path_to_outputs, 'finalresults_iraf', 'fraction_recovered')

    aperture = 3
    targnames = config.target_list[:2] 
    filts = ['F502N'] #config.filter_list
    exp_lengths = ['l']
    flashlvls = [0, 12] #, 6, 18, 24, 33, 55, 91, 116]
    
    for targname in targnames:
        for filt in filts:
            for exp_length in exp_lengths:
                for flashlvl in flashlvls:
                    """
                    plot_slopes(targname=targname, 
                                  filt=filt, 
                                  exp_length=exp_length, 
                                  flashlvl=flashlvl, 
                                  #ctecorr=ctecorr, 
                                  aperture=aperture,
                                  outloc=outloc_compare,
                                  ylims=[-0.1, 0.6])
                    """
                    plot_numpoints(targname='ngc6791', 
                                   filt=filt, 
                                   exp_length=exp_length, 
                                   flashlvl=flashlvl, 
                                   aperture=aperture,
                                   outloc=outloc_frac, 
                                   ylims=[0, 1.0])
                    

if __name__=='__main__':
    compare_idl_python()

